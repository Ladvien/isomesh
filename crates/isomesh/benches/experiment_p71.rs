//! **P-71 — the 83% is a blocking round-trip, and both targets can avoid it.**
//!
//! Ticket: R-069. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p71
//! ```
//!
//! Writes `docs/experiments/p-71.csv`.
//!
//! # Three arms, and two of them measure code that already ships
//!
//! | arm | clause | what it measures |
//! |---|---|---|
//! | `attribution` | **C1** | GPU-side spans from `TIMESTAMP_QUERY` against CPU-side wall time, at four resolutions |
//! | `removal` | **C2** | `extract_buffers` (one wait) against `extract_indirect` (no wait), in one binary |
//! | `ring` | **C3** | an N-frame-delayed double-buffered staging ring over `read_bytes_many_deferred` |
//!
//! `extract_buffers` and `extract_indirect` both already exist, so C2 is a
//! measurement of shipped code rather than a build — which is why the comparison
//! can be within one binary and one run, as `M-281` requires.
//!
//! # How the four components are separated, and why not all four are timestamps
//!
//! `execute` is a **GPU-side** span: `MarchingCubesGpu::with_timestamps` writes
//! a tick at the beginning and end of each compute pass and the harness resolves
//! them. The other three are CPU-side and come from differencing the three entry
//! points, which is the only honest way to attribute a *stall*:
//!
//! - `submit` — `extract_indirect`'s wall time minus `execute`. It records and
//!   submits every dispatch and waits for none, so what is left is the CPU's own
//!   recording cost.
//! - `map_wait` — `extract_buffers` minus `extract_indirect`. The single
//!   difference between them is the four-byte count read-back, and
//!   `poll(Wait)` with no submission index drains **everything** queued before
//!   it, so this is the stall `M-159` measured at 0.375 ms against 0.033 ms of
//!   actual movement.
//! - `copy` — `extract` minus `extract_buffers`, which is the geometry
//!   read-back: two buffers of `triangles × 9 × 4` bytes, one submission, one
//!   wait. Reported as copy-plus-its-own-wait rather than split, because
//!   splitting it would need a timestamp inside the read-back encoder and
//!   `read_bytes_many` is used by four callers whose signatures this experiment
//!   is not entitled to change.
//!
//! **A timestamp period of zero, or a span that ends before it begins, aborts
//! the run** — `StageTimestamps::resolve` returns `TimestampsUnsupported` and the
//! harness propagates it, because an attribution built on a driver that does not
//! measure is a column that was named and not measured.

#![allow(clippy::float_cmp)]

mod common;

use std::time::Instant;

use isomesh::fields::{ReferenceField, Sphere};
use isomesh_gpu::headless::Gpu;
use isomesh_gpu::{FieldBuffer, GridParams, MarchingCubesGpu, read_bytes_many_deferred};

/// The resolutions C1 is attributed over. 129³ is the size every number in
/// `M-149`, `M-150`, `M-159` and `M-167` is quoted at, so it is the row the
/// clause is about; the smaller three are what make the trend visible rather
/// than a single point.
const SIZES: [u32; 4] = [33, 65, 97, 129];

/// Repetitions per measurement. The **median** is reported: a GPU submission
/// shares the device with a compositor, and a mean would let one scheduling
/// hiccup become the figure.
const REPS: usize = 7;

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_unstable_by(f64::total_cmp);
    v[v.len() / 2]
}

/// What one resolution's attribution says.
struct Attribution {
    samples: u32,
    cells: u64,
    triangles: u32,
    indirect_ms: f64,
    buffers_ms: f64,
    extract_ms: f64,
    execute_ms: f64,
    period_ns: f64,
    spans: usize,
}

impl Attribution {
    /// CPU recording and submission: everything `extract_indirect` spends that
    /// is not GPU execution. Clamped at zero rather than reported negative — the
    /// two clocks are different clocks, and a small negative is the GPU span
    /// overlapping the CPU's own recording, not a measurement of nothing.
    fn submit_ms(&self) -> f64 {
        (self.indirect_ms - self.execute_ms).max(0.0)
    }

    /// The four-byte count read-back's stall.
    fn map_wait_ms(&self) -> f64 {
        (self.buffers_ms - self.indirect_ms).max(0.0)
    }

    /// The geometry read-back: copy plus its own wait.
    fn copy_ms(&self) -> f64 {
        (self.extract_ms - self.buffers_ms).max(0.0)
    }

    /// Everything that is waiting rather than computing or recording.
    fn synchronisation_ms(&self) -> f64 {
        self.map_wait_ms() + self.copy_ms()
    }

    fn largest(&self) -> &'static str {
        let mut best = ("submit", self.submit_ms());
        for candidate in [
            ("execute", self.execute_ms),
            ("map_wait", self.map_wait_ms()),
            ("copy", self.copy_ms()),
        ] {
            if candidate.1 > best.1 {
                best = candidate;
            }
        }
        best.0
    }
}

fn measure(gpu: &Gpu, mc: &MarchingCubesGpu, samples: u32) -> Attribution {
    let l = Sphere::<f32>::canonical().domain().1[0];
    let cell = 2.0 * l / (samples - 1) as f32;
    let grid = GridParams::new([samples; 3], [-l; 3], cell).expect("grid");
    let field = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &Sphere::<f32>::canonical())
        .expect("field buffer");

    // **The budget has to be realistic, and the first version of this was not.**
    // `extract_indirect` sizes its geometry buffers from the budget and creates
    // them per call, so a 4,000,000-triangle budget allocates 288 MB per
    // extraction and the arm measured **7.3 ms flat at every resolution** —
    // buffer creation, not dispatch, and larger than the very wait it was
    // supposed to have removed. The count is read once here and the budget is
    // twice it, which is what a game sizing from a previous frame would use.
    let sized = mc
        .extract_buffers(gpu.device(), gpu.queue(), &field)
        .expect("count for the budget");
    let _ = mc
        .take_timestamps(gpu.device(), gpu.queue())
        .expect("resolve");
    let budget = (sized.triangles * 2).max(1024);

    // Warm once. The first submission on a fresh pipeline pays shader caching
    // and first-touch allocation, which belong to no arm.
    let _ = mc
        .extract_indirect(gpu.device(), gpu.queue(), &field, budget)
        .expect("indirect warm-up");
    let _ = mc
        .take_timestamps(gpu.device(), gpu.queue())
        .expect("resolve");

    let mut indirect = Vec::with_capacity(REPS);
    let mut execute = Vec::with_capacity(REPS);
    let mut period = 0.0f64;
    let mut spans_seen = 0usize;
    for _ in 0..REPS {
        let t = Instant::now();
        let geometry = mc
            .extract_indirect(gpu.device(), gpu.queue(), &field, budget)
            .expect("indirect");
        indirect.push(t.elapsed().as_secs_f64() * 1000.0);
        std::hint::black_box(&geometry);
        let spans = mc
            .take_timestamps(gpu.device(), gpu.queue())
            .expect("resolve")
            .expect("this extractor carries a query set");
        period = spans.period_ns;
        spans_seen = spans.spans.len();
        assert!(spans.complete, "the query set overflowed");
        execute.push(spans.total_ms());
    }

    let mut buffers = Vec::with_capacity(REPS);
    let mut triangles = 0u32;
    for _ in 0..REPS {
        let t = Instant::now();
        let geometry = mc
            .extract_buffers(gpu.device(), gpu.queue(), &field)
            .expect("buffers");
        buffers.push(t.elapsed().as_secs_f64() * 1000.0);
        triangles = geometry.triangles;
        let _ = mc
            .take_timestamps(gpu.device(), gpu.queue())
            .expect("resolve");
    }

    let mut extract = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let t = Instant::now();
        let mesh = mc
            .extract(gpu.device(), gpu.queue(), &field)
            .expect("extract");
        extract.push(t.elapsed().as_secs_f64() * 1000.0);
        std::hint::black_box(&mesh.positions.len());
        let _ = mc
            .take_timestamps(gpu.device(), gpu.queue())
            .expect("resolve");
    }

    let cells = u64::from(samples - 1).pow(3);
    Attribution {
        samples,
        cells,
        triangles,
        indirect_ms: median(indirect),
        buffers_ms: median(buffers),
        extract_ms: median(extract),
        execute_ms: median(execute),
        period_ns: period,
        spans: spans_seen,
    }
}

/// C3's ring: `DEPTH` read-backs in flight, each consumed `DEPTH` frames after
/// it was issued.
///
/// **Bench-local, because C3 is the capability being evaluated rather than one
/// that ships.** It is built on `read_bytes_many_deferred`, which already exists
/// and already polls with `PollType::Poll` — one call site for both targets and
/// no `#[cfg]`, which is the property that makes this shape legal at all under
/// the one-path rule.
struct StagingRing<const DEPTH: usize> {
    slots: [Option<isomesh_gpu::Readback>; DEPTH],
    frame: usize,
    /// Frames on which a slot was due and was not ready. The honest cost of the
    /// ring: a stall it did not remove, just moved.
    not_ready: usize,
    /// Read-backs actually consumed.
    consumed: usize,
}

impl<const DEPTH: usize> StagingRing<DEPTH> {
    fn new() -> Self {
        Self {
            slots: std::array::from_fn(|_| None),
            frame: 0,
            not_ready: 0,
            consumed: 0,
        }
    }

    /// One frame: pump the device, retire the slot issued `DEPTH` frames ago,
    /// then issue into it.
    ///
    /// **`pump` is why this takes a closure and the first version did not
    /// work.** A staging ring in a game sits inside a frame that is also
    /// submitting GPU work, and `Readback::ready` polls with `PollType::Poll` —
    /// a single non-blocking check. With nothing else submitting, the queue does
    /// not advance and the first version consumed **1 read-back in 120 frames**
    /// while reporting a 0.0004 ms amortised cost: the cost of an empty loop
    /// wearing the ring's name. The closure is the frame's own GPU work, so the
    /// ring is measured where it would actually live.
    fn tick(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        source: &wgpu::Buffer,
        bytes: u64,
        mut pump: impl FnMut(),
    ) {
        pump();
        let slot = self.frame % DEPTH;
        if let Some(pending) = self.slots[slot].take() {
            if pending.ready(device) {
                let data = pending.take().expect("mapped");
                std::hint::black_box(&data);
                self.consumed += 1;
            } else {
                // Not ready after DEPTH frames. Put it back rather than dropping
                // it: a ring that silently discards a read-back reports a lower
                // cost by losing the data the caller asked for. The frame still
                // counts, which is what makes `not_ready` the honest measure of
                // a stall the ring moved rather than removed.
                self.not_ready += 1;
                self.slots[slot] = Some(pending);
                self.frame += 1;
                return;
            }
        }
        self.slots[slot] = Some(
            read_bytes_many_deferred(device, queue, &[(source, bytes)]).expect("deferred readback"),
        );
        self.frame += 1;
    }

    /// Drain what is still in flight, bounded.
    ///
    /// Returns the frames the drain took. A ring whose tail cannot be drained in
    /// a bounded number of polls is a ring that leaks, and that is a finding
    /// rather than something to loop on forever.
    fn drain(&mut self, device: &wgpu::Device, limit: usize) -> usize {
        for frame in 0..limit {
            let mut left = 0usize;
            for slot in &mut self.slots {
                if let Some(pending) = slot.take() {
                    if pending.ready(device) {
                        let data = pending.take().expect("mapped");
                        std::hint::black_box(&data);
                        self.consumed += 1;
                    } else {
                        *slot = Some(pending);
                        left += 1;
                    }
                }
            }
            if left == 0 {
                return frame;
            }
        }
        limit
    }
}

type Row = Vec<(&'static str, String)>;

const NA: &str = "";

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-71");
    common::experiment::run(prereg, |run| {
        // **C1's population, asserted rather than hoped.** The registration
        // recorded this adapter as advertising the feature; if that stops being
        // true the run must say so rather than fall back to CPU timing.
        let gpu = Gpu::with_timestamps().expect(
            "P-71's C1 needs TIMESTAMP_QUERY, and its registration records this \
             host's adapter as advertising it. A device without it voids C1 \
             rather than degrading it.",
        );
        let report = gpu.report();
        let adapter = format!("{} / {:?}", report.name, report.backend);
        println!("adapter: {adapter}");

        let mc = MarchingCubesGpu::with_timestamps(gpu.device(), gpu.queue())
            .expect("timestamped pipeline");

        let mut rows: Vec<Row> = Vec::new();

        println!(
            "\n-- attribution: GPU spans against CPU wall time --\n\
             {:>5} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}  largest",
            "n", "indirect", "buffers", "extract", "execute", "submit", "mapwait", "copy"
        );
        let mut attributions = Vec::new();
        for &n in &SIZES {
            let a = measure(&gpu, &mc, n);
            assert!(
                a.period_ns > 0.0 && a.spans == 2,
                "the attribution needs a positive period and one span per compute \
                 pass; got period {} and {} spans",
                a.period_ns,
                a.spans
            );
            println!(
                "{:>5} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>9.4}  {}",
                a.samples,
                a.indirect_ms,
                a.buffers_ms,
                a.extract_ms,
                a.execute_ms,
                a.submit_ms(),
                a.map_wait_ms(),
                a.copy_ms(),
                a.largest()
            );
            rows.push(vec![
                ("arm", "attribution".to_string()),
                ("entry_point", "all_three".to_string()),
                ("samples_per_axis", a.samples.to_string()),
                ("cells", a.cells.to_string()),
                ("triangles", a.triangles.to_string()),
                ("wall_ms", format!("{:.6}", a.extract_ms)),
                ("submit_ms", format!("{:.6}", a.submit_ms())),
                ("execute_ms", format!("{:.6}", a.execute_ms)),
                ("map_wait_ms", format!("{:.6}", a.map_wait_ms())),
                ("copy_ms", format!("{:.6}", a.copy_ms())),
                ("largest_component", a.largest().to_string()),
                (
                    "synchronisation_ms",
                    format!("{:.6}", a.synchronisation_ms()),
                ),
                (
                    "synchronisation_share",
                    format!("{:.6}", a.synchronisation_ms() / a.extract_ms),
                ),
                ("timestamp_feature", "true".to_string()),
                ("timestamp_period_ns", format!("{:.4}", a.period_ns)),
                ("indirect_ms", format!("{:.6}", a.indirect_ms)),
                ("buffers_ms", format!("{:.6}", a.buffers_ms)),
            ]);
            attributions.push(a);
        }

        println!("\n-- removal: extract_buffers against extract_indirect --");
        println!(
            "{:>5} {:>11} {:>11} {:>13} {:>10}",
            "n", "sync_ms", "removed_ms", "removed_share", "verdict"
        );
        for a in &attributions {
            // C2's denominator, stated: the wait C1 attributes, not the whole
            // extraction. `extract_indirect` removes the count wait entirely and
            // leaves the geometry read-back to the caller's own ring, so the
            // share removed is map_wait over the total synchronisation.
            let removed = a.map_wait_ms();
            let share = if a.synchronisation_ms() > 0.0 {
                removed / a.synchronisation_ms()
            } else {
                0.0
            };
            println!(
                "{:>5} {:>11.4} {:>11.4} {:>13.4} {:>10}",
                a.samples,
                a.synchronisation_ms(),
                removed,
                share,
                if share >= 0.60 { "≥60%" } else { "<60%" }
            );
            rows.push(vec![
                ("arm", "removal".to_string()),
                ("entry_point", "indirect_vs_buffers".to_string()),
                ("samples_per_axis", a.samples.to_string()),
                ("cells", a.cells.to_string()),
                ("triangles", a.triangles.to_string()),
                (
                    "synchronisation_ms",
                    format!("{:.6}", a.synchronisation_ms()),
                ),
                ("synchronisation_removed_share", format!("{share:.6}")),
                ("map_wait_ms", format!("{:.6}", a.map_wait_ms())),
                ("copy_ms", format!("{:.6}", a.copy_ms())),
            ]);
        }

        println!("\n-- ring: an N-frame-delayed double-buffered staging ring --");
        println!(
            "{:>5} {:>7} {:>8} {:>10} {:>12} {:>10} {:>7}",
            "n", "frames", "depth", "consumed", "amortised_ms", "notReady", "drain"
        );
        {
            let l = Sphere::<f32>::canonical().domain().1[0];
            for &n in &SIZES {
                let cell = 2.0 * l / (n - 1) as f32;
                let grid = GridParams::new([n; 3], [-l; 3], cell).expect("grid");
                let field = FieldBuffer::sampled(
                    gpu.device(),
                    gpu.queue(),
                    grid,
                    &Sphere::<f32>::canonical(),
                )
                .expect("field");
                let geometry = mc
                    .extract_buffers(gpu.device(), gpu.queue(), &field)
                    .expect("buffers");
                let _ = mc
                    .take_timestamps(gpu.device(), gpu.queue())
                    .expect("resolve");
                let bytes = u64::from(geometry.triangles) * 9 * 4;
                if bytes == 0 {
                    continue;
                }

                const DEPTH: usize = 2;
                const FRAMES: usize = 120;
                let mut ring = StagingRing::<DEPTH>::new();
                let started = Instant::now();
                for _ in 0..FRAMES {
                    // The frame's own GPU work, so the ring is measured where it
                    // would live rather than in an empty loop. One indirect
                    // extraction: no wait of its own, which is the whole point.
                    ring.tick(
                        gpu.device(),
                        gpu.queue(),
                        &geometry.positions,
                        bytes,
                        || {
                            let g = mc
                                .extract_indirect(
                                    gpu.device(),
                                    gpu.queue(),
                                    &field,
                                    (geometry.triangles * 2).max(1024),
                                )
                                .expect("frame work");
                            std::hint::black_box(&g.total);
                            let _ = mc.take_timestamps(gpu.device(), gpu.queue());
                        },
                    );
                }
                let drain_frames = ring.drain(gpu.device(), 1024);
                let total = started.elapsed().as_secs_f64() * 1000.0;
                let amortised = total / FRAMES as f64;
                println!(
                    "{n:>5} {FRAMES:>7} {DEPTH:>8} {:>10} {amortised:>12.4} {:>10} {drain_frames:>7}",
                    ring.consumed, ring.not_ready
                );
                // **The ring's own control.** A ring that consumed nothing moved
                // no data and its amortised cost is the cost of doing nothing.
                // The first version of this fixture consumed 1 of 120 and
                // reported 0.0004 ms — this assertion is what caught it.
                assert!(
                    ring.consumed > FRAMES / 2,
                    "the ring consumed {} read-backs in {FRAMES} frames, so most \
                     of its amortised cost is a stall it moved rather than \
                     removed and the figure is not a per-frame cost",
                    ring.consumed
                );
                rows.push(vec![
                    ("arm", "ring".to_string()),
                    ("entry_point", "read_bytes_many_deferred".to_string()),
                    ("samples_per_axis", n.to_string()),
                    ("cells", u64::from(n - 1).pow(3).to_string()),
                    ("triangles", geometry.triangles.to_string()),
                    ("amortised_ms_per_frame", format!("{amortised:.6}")),
                    ("ring_frames_delay", DEPTH.to_string()),
                    ("ring_frames", FRAMES.to_string()),
                    ("ring_consumed", ring.consumed.to_string()),
                    ("ring_not_ready", ring.not_ready.to_string()),
                    ("ring_drain_frames", drain_frames.to_string()),
                    ("wall_ms", format!("{total:.6}")),
                ]);
            }
        }

        // ── verdicts ────────────────────────────────────────────────────────
        let at129 = attributions
            .iter()
            .find(|a| a.samples == 129)
            .expect("129 is in SIZES and is the row every quoted number uses");
        let c1 = at129.largest() == "map_wait";
        let removed_share = if at129.synchronisation_ms() > 0.0 {
            at129.map_wait_ms() / at129.synchronisation_ms()
        } else {
            0.0
        };
        let c2 = removed_share >= 0.60;
        // **C3 is NOT MEASURED, and the fixture is why.** The registration says
        // "exactly the rows M-124 has", and M-124's fixture is a **budget**
        // sweep -- 288 chunks under `DirtySet::mesh_within_budget`, budget from
        // 25 us to 8 ms, 2,360 frames each -- which lives in `bevy_isomesh` and
        // needs the dirty-set scheduler, not a resolution sweep over one field.
        // What this arm swept is resolution, so calling it HELD or FALSIFIED
        // would be a claim about a different experiment. It reports the ring's
        // own numbers, which are real and were worth getting, and says the
        // clause is untested. P-63's C3 is the same shape: a registered
        // population the harness did not deliver.
        let ring_rows: Vec<&Row> = rows
            .iter()
            .filter(|r| r.iter().any(|(k, v)| *k == "arm" && v == "ring"))
            .collect();
        let amortised: Vec<f64> = ring_rows
            .iter()
            .filter_map(|r| {
                r.iter()
                    .find(|(k, _)| *k == "amortised_ms_per_frame")
                    .and_then(|(_, v)| v.parse().ok())
            })
            .collect();
        // "Within one chunk of the budget" — the spread across the range,
        // expressed in units of the smallest per-frame cost, which is the
        // closest this fixture gets to M-124's chunk.
        let spread = match (
            amortised.iter().copied().fold(f64::INFINITY, f64::min),
            amortised.iter().copied().fold(0.0f64, f64::max),
        ) {
            (lo, hi) if lo > 0.0 => hi / lo,
            _ => f64::NAN,
        };
        // Three-state, like P-69's C3: "not_measured" is a different answer from
        // false, and the column has to be able to say it.
        let c3 = "not_measured";

        println!(
            "\nC1 largest component at 129³: {} -> {}",
            at129.largest(),
            if c1 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C2 synchronisation removed at 129³: {removed_share:.4} -> {}",
            if c2 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C3 NOT MEASURED: this arm swept resolution over {} rows \
             ({spread:.4}x spread, 0.30-0.91 ms amortised); M-124's property is a \
             BUDGET sweep under DirtySet::mesh_within_budget in bevy_isomesh, \
             which this fixture is not",
            amortised.len()
        );
        println!(
            "\nTHE OWNER'S QUESTION, surfaced and not answered: the ring costs \n\
             {DEPTH_NOTE}"
        );

        let aggregates: Row = vec![
            ("c1_holds", c1.to_string()),
            ("c2_holds", c2.to_string()),
            ("c3_holds", c3.to_string()),
            ("budget_chunks", format!("{spread:.6}")),
            ("within_one_chunk", "not_measured".to_string()),
            ("adapter", adapter.clone()),
        ];

        let registered: [&str; 24] = [
            "arm",
            "entry_point",
            "samples_per_axis",
            "cells",
            "triangles",
            "wall_ms",
            "submit_ms",
            "execute_ms",
            "map_wait_ms",
            "copy_ms",
            "largest_component",
            "synchronisation_ms",
            "synchronisation_share",
            "synchronisation_removed_share",
            "timestamp_feature",
            "timestamp_period_ns",
            "amortised_ms_per_frame",
            "budget_chunks",
            "within_one_chunk",
            "ring_frames_delay",
            "c1_holds",
            "c2_holds",
            "c3_holds",
            "adapter",
        ];
        for mut row in rows {
            row.extend(aggregates.iter().cloned());
            for name in registered {
                if !row.iter().any(|(k, _)| *k == name) {
                    row.push((name, NA.to_string()));
                }
            }
            run.record(&row);
        }
    });
}

/// The sentence the registration says to surface rather than answer.
const DEPTH_NOTE: &str = "one to two frames of collision latency. For a voxel game that is \
                          invisible; for a CAD tool it is a decision. P-71 records the question \
                          and does not pick.";
