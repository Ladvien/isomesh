//! **P-93 — pricing `M-377`'s unpriced bill: the vertex upload against the 51× edit win.**
//!
//! Ticket: R-093. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p93
//! ```
//!
//! Writes `docs/experiments/p-93.csv`.
//!
//! # SHARE, recomputed from the committed files before a line of this was written
//!
//! `✗54`'s row is read out of `docs/measurements/gpu_vs_cpu.csv` rather than
//! quoted (`P-70`'s precedent): at 129³ `upload_ms` is **7.3236** of an
//! **8.3694** ms `gpu_total_ms`, **87.50%**, over a payload of
//! `129³ · 4 = 8,586,756` bytes — an **effective 1.1725 GB/s**.
//!
//! **That 87.50% is not the share C1 needs, and saying so is the first result
//! here.** `mesh_render.rs` states what the column is: *"field evaluation on the
//! CPU rather than anything a renderer can fix: `FieldBuffer::sampled` evaluates
//! the SDF host-side and copies the samples over"*. It is the **sample-grid**
//! upload, host evaluation included, and GPU-011a already deleted it by
//! evaluating the field in the shader (8.37 ms → 0.54 ms). `M-377`'s bill is a
//! **vertex** upload: 81,548 raw vertices at 4³ against 53,788 at 64³ for the
//! same 53,110 distinct surface points. Two different sets of bytes.
//!
//! So the reachability arithmetic, from `docs/experiments/p-72.csv`'s own
//! `ms_per_edit` and `vertices` columns:
//!
//! | | `gyroid` | `fbm_terrain` |
//! |---|---:|---:|
//! | edit cost at 4³ | 0.7826 ms | 2.6992 ms |
//! | edit cost at 64³ | 39.6726 ms | 127.4229 ms |
//! | excess vertices at 4³ | 27,760 | 18,011 |
//! | excess bytes (24 B/vertex) | 666,240 | 432,264 |
//! | excess bytes at `✗54`'s 1.1725 GB/s | 0.5682 ms | 0.3687 ms |
//! | **full-world excess re-uploads per second needed to put the crossover at 0.1/s** | **6.84** | **33.83** |
//! | the same at a realistic 5 GB/s pure-bus rate | 29.19 | 144.27 |
//!
//! **C1 is arithmetically unreachable and this is said before the run.** To put
//! the crossover at the *bottom* of the registered 0.1–100 window on `gyroid`,
//! the vertex-upload difference between 4³ and 64³ must supply
//! `0.1 · (39.6726 − 0.7826) = 3.889` ms every second. The entire excess payload
//! costs 0.57 ms to move at the committed effective rate, so the whole world's
//! *surplus* geometry would have to cross the bus **6.84 times per second** — and
//! at the registration's own predicted 12.5/s, **855 times per second**. No
//! renderer re-uploads static geometry, so the run is expected to falsify C1. It
//! is run anyway, because the number is the deliverable and because the
//! *magnitude* of the miss is what tells a consumer the granularity decision is
//! settled rather than a trade.
//!
//! # The cost model, stated because the verdict is a function of it
//!
//! Every rate-dependent column is **milliseconds per second of gameplay**:
//!
//! ```text
//! total_ms(c, r) = r · edit_ms(c) + stream_ms_per_second(c)
//!   edit_ms(c)   = (mark + remesh + upload of the re-meshed chunks) / EDITS
//!   stream_ms(c) = real uploads for chunks that ENTERED VIEW, per second
//! ```
//!
//! Both terms are measured on this machine in this run (`M-281`). The edit term
//! is `M-377`'s trace; the stream term is the only place a coarse granularity can
//! win, because a chunk's geometry is uploaded when it *changes* or when it
//! *becomes visible* and never merely because a frame happened. A model that
//! re-uploaded the resident world every frame would hand C1 its crossover by
//! fiat, and would describe no renderer.
//!
//! # Three upload paths, because the first run answered C1 with the wrong mechanism
//!
//! The first version measured streaming with **one `write_buffer` per visible
//! chunk** — which is the shape `bevy_isomesh` ships, a Bevy `Mesh` per chunk
//! entity being a GPU buffer per chunk. It reported a crossover of 0.22 edits/s
//! on `gyroid`, inside the registered window, and **could not say why**. The
//! streaming term spread **12.9×** across the six granularities while the vertex
//! data spread only **2.18×**, because that driver charges about **2 µs per
//! `write_buffer` call** and 2³ makes 38,994 of them where 64³ makes 24. A
//! crossover produced by call count is not the crossover C1 asks about: the
//! clause is denominated in *vertex data*.
//!
//! So all three are reported and the difference between them is the mechanism:
//!
//! | path | what its cost is | column |
//! |---|---|---|
//! | per-chunk `write_buffer` | bytes **+ ~2 µs per chunk** | `crossover_edits_per_second`, `c1_holds` |
//! | arena, one write per frame | bytes **+ one `submit`/`poll` per frame** | `crossover_arena_edits_per_second`, `c1_holds_arena` |
//! | **exact byte counts over the calibrated rate** | **bytes, nothing else** | `crossover_bytes_edits_per_second`, `c1_holds_bytes` |
//!
//! **The third is the one that answers the registration, and it is the only one
//! that is not a clock.** The arena path's residue is a per-frame `submit` +
//! `poll` that the 64³ arm pays across seven frames and the 4³ arm across
//! thirty-two, and seven samples of a sub-millisecond quantity on a machine
//! running eleven other builds came back with the **wrong sign** on `gyroid` —
//! `d_stream` negative, no crossover at all. The byte path is exact integers
//! (`stream_upload_bytes_per_second` is `1.5 ×` the full-world payload, to the
//! byte, because this camera's sinusoidal yaw gives every unit exactly three
//! enters and two exits) over one calibrated rate. `M-337`'s re-audit forced
//! exactly this substitution on `P-40`: gate on counts, not on a ratio that
//! swings further than the effect.
//!
//! # Fixture, inherited from `P-72` verbatim
//!
//! 128³ world cells at every granularity, `EXTENT` 4.0 centred on the fields'
//! own domain, an eleven-edit dig of 6-cell spherical brushes along a straight
//! path whose height is probed per edit at that edit's own `x`, `gyroid` and
//! `fbm_terrain`, median of three traces. Granularities 2³–64³, the range
//! `M-377` had to extend to find its interior optimum.
//!
//! # The three controls, also `P-72`'s
//!
//! - **Every arm must mesh something** — `dirty_chunks > 0`, `raw_vertices > 0`.
//! - **Every arm must produce the same surface** — the sorted multiset of vertex
//!   positions quantised to `cell_size · 1e-6`, which is what `P-61` established
//!   as the load-bearing comparison for a partition change. **The mixed arm is
//!   inside this control**, and that is where it earns its keep: a two-level
//!   scheme that dropped a seam between a coarse chunk and a subdivided
//!   neighbour is exactly the defect a speed number would hide.
//! - **The grid duplication must be the arithmetic one** — field evaluations are
//!   counted and the grid term asserted against `((c+1)/c)³` to four digits,
//!   with the normal stencil derived from the remainder rather than fudged.
//!
//! # The vacuity control the registration names
//!
//! *"the upload arm must actually upload — asserted against `✗54`'s own
//! `gpu_vs_cpu.csv` row rather than quoted"*. Four columns and four `assert!`s,
//! and the failure mode being guarded is specific: `Queue::write_buffer` stages
//! host-side and a harness that never submitted would time a `memcpy` and call it
//! an upload.
//!
//! - `upload_roundtrip_verified` — the calibration payload is uploaded through
//!   the same `Uploader` every row uses and then **read back** with
//!   `isomesh_gpu::read_bytes` and compared byte for byte. Device memory holds
//!   the bytes or the harness stops.
//! - `calibration_upload_bytes` — exactly the committed row's
//!   `samples³ · 4 = 8,586,756`, read from the file.
//! - `calibration_upload_gb_per_s` **below PCIe 4.0 ×16's 31.5 GB/s**. An upload
//!   faster than the bus did not happen; this is the assertion that a zero
//!   cannot pass.
//! - `calibration_upload_gb_per_s` **above the committed 1.1725 GB/s**, because
//!   the committed number carries host field evaluation and this path does not.
//!   Asserted against the file, so the comparison cannot drift from the artefact.
//!
//! And C3's own registered control: `visible_chunk_enters` and
//! `visible_chunk_exits` are both asserted non-zero on every arm including the
//! mixed one. A camera whose visible set only ever grew would be a frustum
//! opening, not a camera moving.
//!
//! # Three fixture defects this harness found in itself, and one it kept
//!
//! 1. **The mixed arm subdivided chunks that are not in the world.** `edit_box`
//!    pads the brush by a cell, so `mark_edit` marks fine chunks at coordinate
//!    −1 and 32; dividing those down produced coarse coordinates −1 and 2 on a
//!    world that has only 0 and 1. The arm reported **11 activations of 8 coarse
//!    chunks** and paid a full 4,096-chunk subdivision for each phantom. Caught
//!    by printing `activations/coarse_chunks` rather than `activations` alone —
//!    a count with no denominator cannot be absurd.
//! 2. **`weld_saves_upload_ms` came back negative** (−0.0093 ms at 32³ on
//!    `gyroid`), which is impossible: welding strictly removes bytes. Two upload
//!    medians of three could not resolve the delta. The saving is now derived
//!    from exact byte counts at the calibrated rate, with the measured
//!    difference kept beside it as `weld_saves_upload_ms_measured` — **and it is
//!    still negative on two `fbm_terrain` arms**, which is the honest statement
//!    of this instrument's floor: about 0.04 ms, wider than the whole saving at
//!    duplication below 1.12×.
//! 3. **C1's first crossover was produced by call count, not by vertex data** —
//!    see the three-paths section above. The decomposition into
//!    `stream_bytes_ms_per_second` and `stream_overhead_ms_per_second` is what
//!    makes that visible in the file rather than only in a reviewer's suspicion.
//!
//! **Kept deliberately:** the fixed arms' dirty sets are *not* filtered to
//! in-world chunks, because that is `M-377`'s behaviour and every integer here
//! reproduces `p-72.csv` exactly — all twelve `raw_vertices` and both
//! `distinct_surface_points`. Filtering them would be a better fixture and a
//! worse comparison.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::print_literal,
    clippy::too_many_lines
)]

mod common;

use std::cell::Cell;
use std::collections::BTreeSet;
use std::time::Instant;

use isomesh::chunk::dirty::{DirtySet, mark_edit};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::extractor::Extractor;
use isomesh::fields::{FbmTerrain, Gyroid};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::weld::{Welder, epsilon_for};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Cells per axis in the whole world, at every granularity. `P-72`'s.
const WORLD_CELLS: u32 = 128;

/// The knob, `M-377`'s range including the two granularities it had to add to
/// find an interior optimum.
const GRANULARITIES: [u32; 6] = [2, 4, 8, 16, 32, 64];

/// World extent per axis.
const EXTENT: f64 = 4.0;

/// World origin, centred on the reference fields' own domain centre.
const ORIGIN: f64 = -EXTENT * 0.5;

/// Brush radius in cells.
const BRUSH_CELLS: f64 = 6.0;

/// Edits in the trace.
const EDITS: usize = 11;

/// Traces per arm, median taken.
const REPS: usize = 3;

/// The rate sweep C1 is denominated in. `12.5` is `game_dig`'s own
/// `EDIT_PERIOD = 0.08` inverted, and the registration names it.
const RATES: [f64; 11] = [
    0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 12.5, 20.0, 50.0, 100.0,
];

/// The rate C3 is scored at: `game_dig`'s throttled stroke.
const GAME_RATE: f64 = 12.5;

/// The registered window C1's crossover has to land inside.
const RATE_MIN: f64 = 0.1;
/// Upper end of that window.
const RATE_MAX: f64 = 100.0;

/// Bytes a renderer uploads per vertex: `position` and `normal` as `f32×3`.
///
/// The interleaved layout a vertex-pulling or vertex-buffer pipeline wants, and
/// what `bevy_isomesh` builds. Indices are counted separately, at 4 bytes.
const BYTES_PER_VERTEX: usize = 24;

/// PCIe 4.0 ×16, the bus this rig's RTX 3090 sits on: 16 lanes × 1.969 GB/s.
///
/// The ceiling the calibration upload must stay under. A transfer that beats the
/// bus is a transfer that did not happen, which is the `M-44` failure this
/// experiment is most exposed to.
const PCIE4_X16_GB_PER_S: f64 = 31.5;

/// Frames the camera arm simulates per second.
const CAMERA_HZ: usize = 60;

/// Seconds of camera motion. One full yaw oscillation.
const CAMERA_SECONDS: f64 = 2.0;

/// Yaw amplitude in radians. `game_dig`'s mouse sensitivity is 0.0022 rad per
/// pixel, so this sweep's peak rate of 4.4 rad/s is a 2,000 px/s flick — fast,
/// and chosen so the world leaves the frustum completely at the extremes rather
/// than merely shifting inside it.
const CAMERA_YAW: f64 = 1.4;

/// Cone half-angle cosine, `cos(35°)`. Bevy's default perspective is a 45°
/// vertical field of view; 35° about the axis is the horizontal half-angle of a
/// 16:9 frustum with a little margin.
const CAMERA_COS_HALF: f64 = 0.819_152_044_288_991_8;

/// Upload staging capacity. The largest single payload is the calibration's
/// 8.6 MB; the largest world payload is 2³'s ~4.1 MB.
const UPLOAD_CAPACITY: u64 = 32 << 20;

/// A field wrapper that counts `sample` calls, for `P-72`'s control 3.
struct Counted<'a, F> {
    field: &'a F,
    calls: &'a Cell<u64>,
}

impl<F: Sdf<Scalar = f64>> Sdf for Counted<'_, F> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        self.calls.set(self.calls.get() + 1);
        self.field.sample(p)
    }
}

/// A sphere subtracted from a field: `max(field, -(|p - c| - r))`. `P-72`'s.
struct Dug<'a, F> {
    field: &'a F,
    centres: &'a [[f64; 3]],
    radius: f64,
}

impl<F: Sdf<Scalar = f64>> Sdf for Dug<'_, F> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        let mut v = self.field.sample(p);
        for c in self.centres {
            let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
            let sphere = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - self.radius;
            v = v.max(-sphere);
        }
        v
    }
}

/// The one path anything in this bench uploads through.
///
/// One persistent `VERTEX | INDEX | COPY_DST | COPY_SRC` buffer, written with
/// `Queue::write_buffer`, flushed with an empty command buffer and waited on with
/// `poll(Wait)`. `COPY_SRC` is there for the round-trip control and for no other
/// reason.
struct Uploader<'a> {
    device: &'a wgpu::Device,
    queue: &'a wgpu::Queue,
    buffer: wgpu::Buffer,
}

impl<'a> Uploader<'a> {
    fn new(device: &'a wgpu::Device, queue: &'a wgpu::Queue) -> Self {
        let buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("p93 vertex upload"),
            size: UPLOAD_CAPACITY,
            usage: wgpu::BufferUsages::VERTEX
                | wgpu::BufferUsages::INDEX
                | wgpu::BufferUsages::COPY_DST
                | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        Self {
            device,
            queue,
            buffer,
        }
    }

    /// Upload every payload with **one `write_buffer` per payload** and wait
    /// until the device has them. `(ms, bytes)`.
    ///
    /// One `submit` and one `poll` for the whole batch, because that is what a
    /// renderer does with a frame's worth of dirty chunks — a submission per
    /// chunk would measure the driver's synchronisation cost per chunk and call
    /// it bandwidth.
    ///
    /// This is the shape `bevy_isomesh` ships: a Bevy `Mesh` per chunk entity is
    /// a GPU buffer per chunk, so the write count is the chunk count.
    fn upload(&self, payloads: &[&[u8]]) -> (f64, u64) {
        let total: u64 = payloads.iter().map(|p| p.len() as u64).sum();
        assert!(
            total <= UPLOAD_CAPACITY,
            "payload of {total} bytes exceeds the {UPLOAD_CAPACITY}-byte staging buffer"
        );
        let started = Instant::now();
        let mut offset = 0u64;
        for payload in payloads {
            if payload.is_empty() {
                continue;
            }
            self.queue.write_buffer(&self.buffer, offset, payload);
            offset += payload.len() as u64;
        }
        self.flush();
        (started.elapsed().as_nanos() as f64 / 1e6, total)
    }

    /// The same bytes through **one** `write_buffer`, via a host arena.
    ///
    /// **This exists because the first run's C1 verdict came from the wrong
    /// mechanism and the harness could not say so.** `upload` costs about 1.3 µs
    /// per call on this driver, and a fine granularity makes tens of thousands
    /// of calls where 64³ makes eight — so the streaming term spread 12.9×
    /// across a vertex-data range of only 2.18×, and the crossover C1 asks about
    /// was being produced by **call count** rather than by the vertex volume
    /// `M-377` named. An arena renderer concatenates its visible chunks into one
    /// buffer and issues one write, at which point the cost is bytes and nothing
    /// else. Both designs are real, so both are measured and both crossovers are
    /// reported; the difference between them *is* the mechanism.
    ///
    /// The `memcpy` into the arena is inside the timer: an arena renderer pays it.
    fn upload_arena(&self, payloads: &[&[u8]], scratch: &mut Vec<u8>) -> (f64, u64) {
        let total: u64 = payloads.iter().map(|p| p.len() as u64).sum();
        assert!(
            total <= UPLOAD_CAPACITY,
            "payload of {total} bytes exceeds the {UPLOAD_CAPACITY}-byte staging buffer"
        );
        let started = Instant::now();
        scratch.clear();
        for payload in payloads {
            scratch.extend_from_slice(payload);
        }
        if !scratch.is_empty() {
            self.queue.write_buffer(&self.buffer, 0, scratch);
        }
        self.flush();
        (started.elapsed().as_nanos() as f64 / 1e6, total)
    }

    /// Submit an empty command buffer and wait. The staged writes go with it.
    fn flush(&self) {
        let encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("p93 upload flush"),
            });
        self.queue.submit(Some(encoder.finish()));
        let _ = self.device.poll(wgpu::PollType::Wait {
            submission_index: None,
            timeout: None,
        });
    }
}

/// A renderable unit: one chunk's box and one chunk's vertex payload.
///
/// The camera arm is generic over these, which is what lets the mixed arm mix
/// coarse and fine units in one list instead of needing its own visibility code.
struct Unit {
    min: [f64; 3],
    size: f64,
    bytes: Vec<u8>,
}

/// One granularity's meshed world, at the final field.
struct World {
    raw_vertices: usize,
    triangles: usize,
    surface: BTreeSet<[i64; 3]>,
    units: Vec<Unit>,
    combined: MeshBuffer<f64>,
    payload: Vec<u8>,
}

/// A group of chunks meshed at one granularity.
struct Group<'a> {
    layout: &'a ChunkLayout<f64>,
    shape: RuntimeShape3,
    ids: Vec<ChunkId>,
}

/// What one camera sweep did.
///
/// `upload_ms` and `arena_upload_ms` move **exactly the same bytes** in exactly
/// the same frames; they differ only in how many `write_buffer` calls it took.
struct Stream {
    enters: usize,
    exits: usize,
    upload_ms: f64,
    arena_upload_ms: f64,
    upload_bytes: u64,
    writes: u64,
    frames_uploading: usize,
    units: usize,
}

/// One fixed-granularity arm.
struct Arm {
    chunk_cells: u32,
    chunks: u64,
    samples_per_chunk: u64,
    field_calls: u64,
    build_vertices: u64,
    dirty_chunks: u64,
    mark_ms: f64,
    remesh_ms: f64,
    upload_dirty_ms: f64,
    upload_dirty_arena_ms: f64,
    upload_dirty_bytes: u64,
    world: World,
    weld_ms: f64,
    welded_vertices: usize,
    welded_removed: usize,
    welded_payload_bytes: u64,
    full_upload_ms: f64,
    welded_upload_ms: f64,
    stream: Stream,
}

impl Arm {
    /// Milliseconds one edit costs: mark, remesh and the upload of what was
    /// re-meshed.
    fn edit_ms(&self) -> f64 {
        (self.mark_ms + self.remesh_ms + self.upload_dirty_ms) / EDITS as f64
    }

    /// The same, on the arena upload path.
    fn edit_ms_arena(&self) -> f64 {
        (self.mark_ms + self.remesh_ms + self.upload_dirty_arena_ms) / EDITS as f64
    }

    /// Milliseconds per second of gameplay the camera's streaming costs.
    fn stream_ms_per_second(&self) -> f64 {
        self.stream.upload_ms / CAMERA_SECONDS
    }

    /// The same, on the arena upload path — bytes and nothing else.
    fn stream_arena_ms_per_second(&self) -> f64 {
        self.stream.arena_upload_ms / CAMERA_SECONDS
    }

    fn total_ms(&self, rate: f64) -> f64 {
        rate * self.edit_ms() + self.stream_ms_per_second()
    }

    fn total_arena_ms(&self, rate: f64) -> f64 {
        rate * self.edit_ms_arena() + self.stream_arena_ms_per_second()
    }

    /// What welding the world's duplicates saves on the world's upload,
    /// **derived from the byte counts at the calibrated transfer rate**.
    ///
    /// **The measured difference of two ~0.4 ms medians could not resolve it and
    /// the first run proved that**: `weld_saves_upload_ms` came back
    /// **negative** at 32³ on `gyroid`, which says the upload timer's noise floor
    /// is wider than the byte delta the weld removes. The byte counts are exact
    /// integers and the rate is one measurement shared by every row, which is the
    /// same substitution `M-337`'s re-audit forced on `P-40`: gate on counts, not
    /// on a clock that swings further than the effect.
    fn weld_saving_ms(&self, gb_per_s: f64) -> f64 {
        let saved = self.world.payload.len() as f64 - self.welded_payload_bytes as f64;
        saved / (gb_per_s * 1e6)
    }

    /// The same saving as two upload medians actually measured, for comparison.
    fn weld_saving_ms_measured(&self) -> f64 {
        self.full_upload_ms - self.welded_upload_ms
    }
}

/// The mixed arm: coarse for static chunks, fine for actively-dug ones.
struct Mixed {
    coarse_cells: u32,
    fine_cells: u32,
    activations: usize,
    coarse_chunks: usize,
    mark_ms: f64,
    remesh_ms: f64,
    subdivide_ms: f64,
    subdivide_upload_ms: f64,
    dirty_upload_ms: f64,
    dirty_upload_arena_ms: f64,
    world: World,
    stream: Stream,
}

impl Mixed {
    /// Per-edit cost including the one-off transition: subdividing a coarse
    /// chunk the first time an edit lands in it, and uploading what that
    /// produced, amortised over the eleven-edit trace the clause names.
    fn edit_ms(&self) -> f64 {
        (self.mark_ms
            + self.remesh_ms
            + self.subdivide_ms
            + self.subdivide_upload_ms
            + self.dirty_upload_ms)
            / EDITS as f64
    }

    /// Per-edit cost **after** every coarse chunk the trace touches has already
    /// been subdivided.
    ///
    /// The scheduler's best case, and the fairer number for a long session: the
    /// transition is paid once per coarse chunk and never again. Reported beside
    /// `edit_ms` rather than instead of it, because an eleven-edit trace is
    /// mostly transition and the registration denominates C3 in that trace.
    fn steady_edit_ms(&self) -> f64 {
        (self.mark_ms + self.remesh_ms + self.dirty_upload_ms) / EDITS as f64
    }

    fn stream_ms_per_second(&self) -> f64 {
        self.stream.upload_ms / CAMERA_SECONDS
    }

    fn stream_arena_ms_per_second(&self) -> f64 {
        self.stream.arena_upload_ms / CAMERA_SECONDS
    }

    fn total_ms(&self, rate: f64) -> f64 {
        rate * self.edit_ms() + self.stream_ms_per_second()
    }

    fn steady_total_ms(&self, rate: f64) -> f64 {
        rate * self.steady_edit_ms() + self.stream_ms_per_second()
    }

    fn total_arena_ms(&self, rate: f64) -> f64 {
        let edit = (self.mark_ms
            + self.remesh_ms
            + self.subdivide_ms
            + self.subdivide_upload_ms
            + self.dirty_upload_arena_ms)
            / EDITS as f64;
        rate * edit + self.stream_arena_ms_per_second()
    }
}

/// A chunk mesh in the bytes a renderer uploads.
fn payload_bytes(mesh: &MeshBuffer<f64>) -> Vec<u8> {
    assert_eq!(
        mesh.positions.len(),
        mesh.normals.len(),
        "a mesh with {} positions and {} normals is not a vertex buffer, and the byte count \
         would be a guess",
        mesh.positions.len(),
        mesh.normals.len()
    );
    let mut out =
        Vec::with_capacity(mesh.positions.len() * BYTES_PER_VERTEX + mesh.indices.len() * 4);
    for (p, n) in mesh.positions.iter().zip(&mesh.normals) {
        for v in p.iter().chain(n.iter()) {
            out.extend_from_slice(&(*v as f32).to_le_bytes());
        }
    }
    for i in &mesh.indices {
        out.extend_from_slice(&i.to_le_bytes());
    }
    out
}

fn median(mut xs: Vec<f64>) -> f64 {
    xs.sort_unstable_by(f64::total_cmp);
    xs[xs.len() / 2]
}

/// The dig path: eleven brushes along `x`, each at the surface height probed at
/// that edit's own `x`. `P-72`'s, unchanged.
fn dig_path<F: Sdf<Scalar = f64>>(field: &F) -> Vec<[f64; 3]> {
    let mid = ORIGIN + EXTENT * 0.5;
    let surface_y = |x: f64| -> f64 {
        let steps = 1024;
        let mut prev = field.sample([x, ORIGIN, mid]);
        for i in 1..=steps {
            let y = ORIGIN + EXTENT * (f64::from(i) / f64::from(steps));
            let v = field.sample([x, y, mid]);
            if (prev < 0.0) != (v < 0.0) {
                return y;
            }
            prev = v;
        }
        panic!("no surface crossing along y at x = {x}: the trace would dig in empty space");
    };
    (0..EDITS)
        .map(|i| {
            let t = (i as f64 + 0.5) / EDITS as f64;
            let x = ORIGIN + EXTENT * t;
            [x, surface_y(x), mid]
        })
        .collect()
}

/// The cell box one brush can touch, through the layout rather than by hand.
fn edit_box(layout: &ChunkLayout<f64>, centre: [f64; 3], radius: f64) -> ([i64; 3], [i64; 3]) {
    let lo_world = [0, 1, 2].map(|a| centre[a] - radius);
    let hi_world = [0, 1, 2].map(|a| centre[a] + radius);
    (
        layout.cell_of(lo_world).map(|v| v - 1),
        layout.cell_of(hi_world).map(|v| v + 1),
    )
}

/// Mesh a set of chunk groups over one field and collect everything a row needs.
fn world_from<F: Sdf<Scalar = f64>>(
    field: &F,
    mc: &mut MarchingCubes<f64>,
    groups: &[Group<'_>],
    quantum: f64,
) -> World {
    let mut surface: BTreeSet<[i64; 3]> = BTreeSet::new();
    let mut units: Vec<Unit> = Vec::new();
    let mut combined = MeshBuffer::<f64>::new();
    let mut raw_vertices = 0usize;
    let mut triangles = 0usize;
    for group in groups {
        let cell_size = group.layout.cell_size();
        let size = f64::from(group.layout.cells()) * cell_size;
        for id in &group.ids {
            let origin = group.layout.sample_origin(*id);
            let mut out = MeshBuffer::<f64>::new();
            let _ = mc.extract_into(field, &group.shape, origin, cell_size, &mut out);
            if out.positions.is_empty() {
                continue;
            }
            raw_vertices += out.positions.len();
            triangles += out.indices.len() / 3;
            for p in &out.positions {
                surface.insert([0, 1, 2].map(|a| (p[a] / quantum).round() as i64));
            }
            let base = combined.positions.len() as u32;
            combined.positions.extend_from_slice(&out.positions);
            combined.normals.extend_from_slice(&out.normals);
            combined.indices.extend(out.indices.iter().map(|i| i + base));
            units.push(Unit {
                min: origin,
                size,
                bytes: payload_bytes(&out),
            });
        }
    }
    let payload = payload_bytes(&combined);
    World {
        raw_vertices,
        triangles,
        surface,
        units,
        combined,
        payload,
    }
}

/// Is any of a unit's nine probe points inside the view cone?
fn unit_visible(unit: &Unit, eye: [f64; 3], forward: [f64; 3]) -> bool {
    let mut probes = [[0.0f64; 3]; 9];
    probes[0] = [0, 1, 2].map(|a| unit.min[a] + unit.size * 0.5);
    for corner in 0..8u32 {
        probes[1 + corner as usize] = [0, 1, 2].map(|a| {
            unit.min[a] + if corner >> a & 1 == 1 { unit.size } else { 0.0 }
        });
    }
    for q in &probes {
        let d = [0, 1, 2].map(|a| q[a] - eye[a]);
        let len2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if len2 <= 0.0 {
            return true;
        }
        let dot = d[0] * forward[0] + d[1] * forward[1] + d[2] * forward[2];
        if dot > 0.0 && dot * dot >= CAMERA_COS_HALF * CAMERA_COS_HALF * len2 {
            return true;
        }
    }
    false
}

/// Sweep the camera and upload, for real, every unit that entered view.
///
/// The eye sits above and behind the world and yaws sinusoidally, so the world
/// leaves the frustum entirely at the extremes and comes back. Units that leave
/// view cost nothing; units that arrive are uploaded in one batch per frame.
///
/// **Both upload designs are timed on the same frame's bytes**, arena first so
/// the per-chunk path gets the warm staging belt — the conservative order for
/// the claim that per-chunk call overhead is what drives the fine granularities'
/// streaming cost.
fn camera_sweep(uploader: &Uploader<'_>, units: &[Unit]) -> Stream {
    let mid = ORIGIN + EXTENT * 0.5;
    let eye = [mid, mid + 0.6 * EXTENT, mid + 1.5 * EXTENT];
    let base = {
        let d = [0.0, -0.6 * EXTENT, -1.5 * EXTENT];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        [d[0] / len, d[1] / len, d[2] / len]
    };
    let frames = (CAMERA_HZ as f64 * CAMERA_SECONDS) as usize;
    let mut prev = vec![false; units.len()];
    let mut cur = vec![false; units.len()];
    let mut scratch: Vec<u8> = Vec::new();
    let mut enters = 0usize;
    let mut exits = 0usize;
    let mut upload_ms = 0.0f64;
    let mut arena_upload_ms = 0.0f64;
    let mut upload_bytes = 0u64;
    let mut writes = 0u64;
    let mut frames_uploading = 0usize;
    for frame in 0..frames {
        let t = frame as f64 / CAMERA_HZ as f64;
        let yaw = CAMERA_YAW * (core::f64::consts::TAU * t / CAMERA_SECONDS).sin();
        let (s, c) = yaw.sin_cos();
        let forward = [
            base[0] * c + base[2] * s,
            base[1],
            -base[0] * s + base[2] * c,
        ];
        let mut arriving: Vec<&[u8]> = Vec::new();
        for (i, unit) in units.iter().enumerate() {
            let visible = unit_visible(unit, eye, forward);
            cur[i] = visible;
            if visible && !prev[i] {
                enters += 1;
                arriving.push(&unit.bytes);
            } else if !visible && prev[i] {
                exits += 1;
            }
        }
        if !arriving.is_empty() {
            arena_upload_ms += uploader.upload_arena(&arriving, &mut scratch).0;
            let (ms, bytes) = uploader.upload(&arriving);
            upload_ms += ms;
            upload_bytes += bytes;
            writes += arriving.len() as u64;
            frames_uploading += 1;
        }
        prev.copy_from_slice(&cur);
    }
    Stream {
        enters,
        exits,
        upload_ms,
        arena_upload_ms,
        upload_bytes,
        writes,
        frames_uploading,
        units: units.len(),
    }
}

/// Run one fixed-granularity arm end to end.
fn run_arm<F: Sdf<Scalar = f64>>(
    field: &F,
    chunk_cells: u32,
    uploader: &Uploader<'_>,
    centres: &[[f64; 3]],
) -> Arm {
    let cell_size = EXTENT / f64::from(WORLD_CELLS);
    let layout = ChunkLayout::<f64>::new(chunk_cells, cell_size, [ORIGIN; 3]).expect("layout");
    let shape = layout.sample_shape().expect("shape");
    let per_axis = WORLD_CELLS / chunk_cells;
    let radius = BRUSH_CELLS * cell_size;
    let mut mc = MarchingCubes::<f64>::new();

    let all_ids: Vec<ChunkId> = (0..per_axis)
        .flat_map(|cz| {
            (0..per_axis).flat_map(move |cy| {
                (0..per_axis).map(move |cx| ChunkId {
                    coords: [cx as i32, cy as i32, cz as i32],
                })
            })
        })
        .collect();

    // ── the initial build, for control 3's call accounting ───────────────────
    let calls = Cell::new(0u64);
    let mut build_vertices = 0usize;
    for id in &all_ids {
        let counted = Counted {
            field,
            calls: &calls,
        };
        let mut out = MeshBuffer::<f64>::new();
        let _ = mc.extract_into(
            &counted,
            &shape,
            layout.sample_origin(*id),
            cell_size,
            &mut out,
        );
        build_vertices += out.positions.len();
    }
    let build_calls = calls.get();

    // ── the timed trace ─────────────────────────────────────────────────────
    let mut dirty = DirtySet::new();
    let mut dirty_chunks = 0u64;
    let mut reps: Vec<(f64, f64, f64, f64, u64)> = Vec::with_capacity(REPS);
    let mut scratch: Vec<u8> = Vec::new();
    for rep in 0..REPS {
        let mut mark_ns = 0u128;
        let mut remesh_ns = 0u128;
        let mut upload_ms = 0.0f64;
        let mut arena_ms = 0.0f64;
        let mut upload_bytes = 0u64;
        for step in 0..EDITS {
            let before = Dug {
                field,
                centres: &centres[..step],
                radius,
            };
            let after = Dug {
                field,
                centres: &centres[..=step],
                radius,
            };
            let (lo, hi) = edit_box(&layout, centres[step], radius);

            let t = Instant::now();
            let _ = mark_edit(&layout, &before, &after, lo, hi, &mut dirty).expect("mark");
            mark_ns += t.elapsed().as_nanos();

            let mut fresh: Vec<MeshBuffer<f64>> = Vec::new();
            let t = Instant::now();
            let done = dirty.mesh_dirty(&layout, |_id, origin| {
                let mut out = MeshBuffer::<f64>::new();
                let _ = mc.extract_into(&after, &shape, origin, cell_size, &mut out);
                fresh.push(out);
            });
            remesh_ns += t.elapsed().as_nanos();
            if rep == 0 {
                dirty_chunks += done as u64;
            }

            // Packing is inside the upload timer because it is inside the
            // consumer's upload too: `gpu_vs_cpu.csv`'s own `upload_ms` counts
            // the `Vec<u8>` build (`FieldBuffer::from_bytes`), and a column that
            // excluded it would not be comparable to the row this experiment is
            // denominated in.
            let started = Instant::now();
            let payloads: Vec<Vec<u8>> = fresh.iter().map(payload_bytes).collect();
            let refs: Vec<&[u8]> = payloads.iter().map(Vec::as_slice).collect();
            let pack_ms = started.elapsed().as_nanos() as f64 / 1e6;
            arena_ms += pack_ms + uploader.upload_arena(&refs, &mut scratch).0;
            let (ms, bytes) = uploader.upload(&refs);
            upload_ms += pack_ms + ms;
            upload_bytes += bytes;
        }
        reps.push((
            mark_ns as f64 / 1e6,
            remesh_ns as f64 / 1e6,
            upload_ms,
            arena_ms,
            upload_bytes,
        ));
    }
    reps.sort_unstable_by(|a, b| (a.0 + a.1 + a.2).total_cmp(&(b.0 + b.1 + b.2)));
    let (mark_ms, remesh_ms, upload_dirty_ms, upload_dirty_arena_ms, upload_dirty_bytes) =
        reps[REPS / 2];

    // ── the final surface, outside every timer ──────────────────────────────
    let full = Dug {
        field,
        centres,
        radius,
    };
    let world = world_from(
        &full,
        &mut mc,
        &[Group {
            layout: &layout,
            shape,
            ids: all_ids,
        }],
        cell_size * 1e-6,
    );

    // ── C2: the weld, and the upload it saves ───────────────────────────────
    let epsilon = epsilon_for(cell_size);
    let mut welder = Welder::<f64>::new();
    let mut weld_samples: Vec<f64> = Vec::with_capacity(3);
    let mut welded = MeshBuffer::<f64>::new();
    let mut removed = 0usize;
    for _ in 0..3 {
        let mut candidate = world.combined.clone();
        let started = Instant::now();
        let report = welder.weld(&mut candidate, epsilon).expect("weld");
        weld_samples.push(started.elapsed().as_nanos() as f64 / 1e6);
        removed = report.vertices_removed();
        welded = candidate;
    }
    let weld_ms = median(weld_samples);
    let welded_payload = payload_bytes(&welded);

    // Fifteen, not three. The first run's `weld_saves_upload_ms` came back at
    // **-0.0093 ms** on `gyroid` at 32³ — a saving that cannot be negative,
    // because welding strictly removes bytes. Two medians of three, each about
    // 0.4 ms, cannot resolve a 30 kB difference on a governed machine. These
    // medians are now the *cross-check*; `weld_saving_ms` is derived from the
    // exact byte counts at the calibrated rate.
    let full_upload_ms = median(
        (0..15)
            .map(|_| uploader.upload(&[&world.payload]).0)
            .collect(),
    );
    let welded_upload_ms =
        median((0..15).map(|_| uploader.upload(&[&welded_payload]).0).collect());

    let stream = camera_sweep(uploader, &world.units);

    Arm {
        chunk_cells,
        chunks: u64::from(per_axis).pow(3),
        samples_per_chunk: u64::from(chunk_cells + 1).pow(3),
        field_calls: build_calls,
        build_vertices: build_vertices as u64,
        dirty_chunks,
        mark_ms,
        remesh_ms,
        upload_dirty_ms,
        upload_dirty_arena_ms,
        upload_dirty_bytes,
        world,
        weld_ms,
        welded_vertices: welded.positions.len(),
        welded_removed: removed,
        welded_payload_bytes: welded_payload.len() as u64,
        full_upload_ms,
        welded_upload_ms,
        stream,
    }
}

/// Run the mixed arm: `coarse_cells` everywhere until an edit lands, then
/// `fine_cells` inside every coarse chunk the trace has touched.
fn run_mixed<F: Sdf<Scalar = f64>>(
    field: &F,
    coarse_cells: u32,
    fine_cells: u32,
    uploader: &Uploader<'_>,
    centres: &[[f64; 3]],
) -> Mixed {
    let cell_size = EXTENT / f64::from(WORLD_CELLS);
    let coarse = ChunkLayout::<f64>::new(coarse_cells, cell_size, [ORIGIN; 3]).expect("coarse");
    let fine = ChunkLayout::<f64>::new(fine_cells, cell_size, [ORIGIN; 3]).expect("fine");
    let coarse_shape = coarse.sample_shape().expect("coarse shape");
    let fine_shape = fine.sample_shape().expect("fine shape");
    let ratio = (coarse_cells / fine_cells) as i32;
    let coarse_per_axis = (WORLD_CELLS / coarse_cells) as i32;
    let radius = BRUSH_CELLS * cell_size;
    let mut mc = MarchingCubes::<f64>::new();

    let fine_ids_of = |c: [i32; 3]| -> Vec<ChunkId> {
        let mut ids = Vec::with_capacity((ratio * ratio * ratio) as usize);
        for z in 0..ratio {
            for y in 0..ratio {
                for x in 0..ratio {
                    ids.push(ChunkId {
                        coords: [c[0] * ratio + x, c[1] * ratio + y, c[2] * ratio + z],
                    });
                }
            }
        }
        ids
    };

    let mut reps: Vec<(f64, f64, f64, f64, f64, f64)> = Vec::with_capacity(REPS);
    let mut last_active: BTreeSet<[i32; 3]> = BTreeSet::new();
    let mut activations = 0usize;
    let mut scratch: Vec<u8> = Vec::new();
    for _ in 0..REPS {
        // Reset per rep, or reps two and three pay no subdivision at all and the
        // scheduler's transition cost vanishes from the number that is supposed
        // to price it.
        let mut active: BTreeSet<[i32; 3]> = BTreeSet::new();
        let mut dirty = DirtySet::new();
        let mut mark_ns = 0u128;
        let mut remesh_ns = 0u128;
        let mut sub_ns = 0u128;
        let mut sub_upload_ms = 0.0f64;
        let mut dirty_upload_ms = 0.0f64;
        let mut dirty_arena_ms = 0.0f64;
        for step in 0..EDITS {
            let before = Dug {
                field,
                centres: &centres[..step],
                radius,
            };
            let after = Dug {
                field,
                centres: &centres[..=step],
                radius,
            };
            let (lo, hi) = edit_box(&fine, centres[step], radius);

            let t = Instant::now();
            let _ = mark_edit(&fine, &before, &after, lo, hi, &mut dirty).expect("mark");
            mark_ns += t.elapsed().as_nanos();

            // The coarse chunks this edit reached, derived by integer division
            // rather than by a second `mark_edit` over the same region: a
            // scheduler that marked twice would be paying twice for one fact.
            //
            // **Filtered to chunks that are inside the world, and the first run
            // is why.** `edit_box` pads the brush by one cell, so `mark_edit`
            // legitimately marks fine chunks at coordinate -1 and at 32 — and
            // dividing those down gave coarse coordinates -1 and 2 on a world
            // that has only 0 and 1. The 64³/4³ arm reported **11 activations of
            // 8 coarse chunks** and paid a full 4,096-chunk subdivision for each
            // phantom, meshing empty space outside the world and charging it to
            // the scheduler. That is a fixture defect, not a scheduler cost. The
            // fixed arms are left unfiltered on purpose, because their dirty-set
            // behaviour is `M-377`'s and this experiment inherits it.
            let touched: BTreeSet<[i32; 3]> = dirty
                .iter()
                .map(|id| [0, 1, 2].map(|a| id.coords[a].div_euclid(ratio)))
                .filter(|c| c.iter().all(|v| *v >= 0 && *v < coarse_per_axis))
                .collect();

            let mut subdivided: Vec<MeshBuffer<f64>> = Vec::new();
            let t = Instant::now();
            for c in &touched {
                if active.insert(*c) {
                    for id in fine_ids_of(*c) {
                        let mut out = MeshBuffer::<f64>::new();
                        let _ = mc.extract_into(
                            &after,
                            &fine_shape,
                            fine.sample_origin(id),
                            cell_size,
                            &mut out,
                        );
                        if !out.positions.is_empty() {
                            subdivided.push(out);
                        }
                    }
                }
            }
            sub_ns += t.elapsed().as_nanos();

            let mut fresh: Vec<MeshBuffer<f64>> = Vec::new();
            let t = Instant::now();
            let _ = dirty.mesh_dirty(&fine, |_id, origin| {
                let mut out = MeshBuffer::<f64>::new();
                let _ = mc.extract_into(&after, &fine_shape, origin, cell_size, &mut out);
                fresh.push(out);
            });
            remesh_ns += t.elapsed().as_nanos();

            // Two uploads, kept apart, because one is a one-off transition and
            // the other is what every edit costs for ever. Folding them together
            // is what made the first run's mixed arm look 25 ms/edit expensive
            // with no way to see which half it was.
            let started = Instant::now();
            let sub_payloads: Vec<Vec<u8>> = subdivided.iter().map(payload_bytes).collect();
            let sub_refs: Vec<&[u8]> = sub_payloads.iter().map(Vec::as_slice).collect();
            let sub_pack = started.elapsed().as_nanos() as f64 / 1e6;
            if !sub_refs.is_empty() {
                sub_upload_ms += sub_pack + uploader.upload(&sub_refs).0;
            }

            let started = Instant::now();
            let payloads: Vec<Vec<u8>> = fresh.iter().map(payload_bytes).collect();
            let refs: Vec<&[u8]> = payloads.iter().map(Vec::as_slice).collect();
            let pack_ms = started.elapsed().as_nanos() as f64 / 1e6;
            dirty_arena_ms += pack_ms + uploader.upload_arena(&refs, &mut scratch).0;
            dirty_upload_ms += pack_ms + uploader.upload(&refs).0;
        }
        activations = active.len();
        last_active = active;
        reps.push((
            mark_ns as f64 / 1e6,
            remesh_ns as f64 / 1e6,
            sub_ns as f64 / 1e6,
            sub_upload_ms,
            dirty_upload_ms,
            dirty_arena_ms,
        ));
    }
    reps.sort_unstable_by(|a, b| {
        (a.0 + a.1 + a.2 + a.3 + a.4).total_cmp(&(b.0 + b.1 + b.2 + b.3 + b.4))
    });
    let (mark_ms, remesh_ms, subdivide_ms, subdivide_upload_ms, dirty_upload_ms, dirty_upload_arena_ms) =
        reps[REPS / 2];

    // ── the mixed world: coarse where static, fine where dug ────────────────
    let mut coarse_ids: Vec<ChunkId> = Vec::new();
    let mut fine_ids: Vec<ChunkId> = Vec::new();
    let mut coarse_chunks = 0usize;
    for cz in 0..coarse_per_axis {
        for cy in 0..coarse_per_axis {
            for cx in 0..coarse_per_axis {
                coarse_chunks += 1;
                let c = [cx, cy, cz];
                if last_active.contains(&c) {
                    fine_ids.extend(fine_ids_of(c));
                } else {
                    coarse_ids.push(ChunkId { coords: c });
                }
            }
        }
    }
    let full = Dug {
        field,
        centres,
        radius,
    };
    let world = world_from(
        &full,
        &mut mc,
        &[
            Group {
                layout: &coarse,
                shape: coarse_shape,
                ids: coarse_ids,
            },
            Group {
                layout: &fine,
                shape: fine_shape,
                ids: fine_ids,
            },
        ],
        cell_size * 1e-6,
    );
    let stream = camera_sweep(uploader, &world.units);

    Mixed {
        coarse_cells,
        fine_cells,
        activations,
        coarse_chunks,
        mark_ms,
        remesh_ms,
        subdivide_ms,
        subdivide_upload_ms,
        dirty_upload_ms,
        dirty_upload_arena_ms,
        world,
        stream,
    }
}

type Row = Vec<(&'static str, String)>;

/// One field's three C1 crossovers and its three clause verdicts.
struct Verdict {
    field: &'static str,
    per_chunk: f64,
    arena: f64,
    bytes: f64,
    c1: bool,
    c1_arena: bool,
    c1_bytes: bool,
    c2: bool,
    c3: bool,
}

/// The committed row C1 is denominated in, read rather than quoted.
struct Committed {
    upload_ms: f64,
    total_ms: f64,
    bytes: f64,
}

fn committed_row(root: &std::path::Path) -> Committed {
    let text = std::fs::read_to_string(root.join("docs/measurements/gpu_vs_cpu.csv"))
        .expect("docs/measurements/gpu_vs_cpu.csv");
    let mut header: Vec<&str> = Vec::new();
    let mut found = None;
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let cells: Vec<&str> = line.split(',').collect();
        if header.is_empty() {
            header = cells;
            continue;
        }
        let get = |name: &str| -> Option<f64> {
            header
                .iter()
                .position(|h| *h == name)
                .and_then(|i| cells.get(i))
                .and_then(|v| v.parse().ok())
        };
        if get("samples") == Some(129.0) {
            found = Some(Committed {
                upload_ms: get("upload_ms").unwrap_or(f64::NAN),
                total_ms: get("gpu_total_ms").unwrap_or(f64::NAN),
                // The bytes that upload actually moved: an `f32` per grid
                // sample. `FieldBuffer::uploaded` writes exactly this.
                bytes: 129.0 * 129.0 * 129.0 * 4.0,
            });
        }
    }
    let row = found.expect("no 129 row in gpu_vs_cpu.csv, so the vacuity control has no anchor");
    assert!(
        row.upload_ms.is_finite() && row.upload_ms > 0.0 && row.total_ms.is_finite(),
        "gpu_vs_cpu.csv's 129 row has no usable upload_ms, so the upload arm has nothing to be \
         asserted against"
    );
    row
}

/// The CPU clock this run's milliseconds were taken on (`M-280`).
fn cpu_khz() -> String {
    std::fs::read_to_string("/sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq")
        .ok()
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| String::from("unknown"))
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-93");
    let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let committed = committed_row(&root);
    let committed_gb_per_s = committed.bytes / committed.upload_ms / 1e6;

    println!("the denominator, read from docs/measurements/gpu_vs_cpu.csv at 129³:");
    println!("  upload_ms      {:.4}", committed.upload_ms);
    println!("  gpu_total_ms   {:.4}", committed.total_ms);
    println!(
        "  share          {:.2}%  -- of the SAMPLE-GRID upload, host field evaluation included",
        100.0 * committed.upload_ms / committed.total_ms
    );
    println!("  payload        {} bytes", committed.bytes as u64);
    println!("  effective rate {committed_gb_per_s:.4} GB/s\n");

    // ── the device, and the vacuity control ─────────────────────────────────
    let gpu = isomesh_gpu::headless::Gpu::new().expect("a GPU adapter; P-93 is VOID without one");
    let report = gpu.report();
    println!(
        "adapter: {} ({:?}, {:?}) driver {}\n",
        report.name, report.backend, report.device_type, report.driver
    );
    let uploader = Uploader::new(gpu.device(), gpu.queue());

    // Exactly the committed row's payload, through exactly the path every row
    // uses, then read back off the device and compared byte for byte.
    let calib_bytes = committed.bytes as usize;
    let calib: Vec<u8> = (0..calib_bytes).map(|i| (i % 251) as u8).collect();
    let calib_ms = median((0..3).map(|_| uploader.upload(&[&calib]).0).collect());
    let returned = isomesh_gpu::read_bytes(
        gpu.device(),
        gpu.queue(),
        &uploader.buffer,
        calib_bytes as u64,
    )
    .expect("read back the calibration payload");
    let roundtrip = returned == calib;
    assert!(
        roundtrip,
        "the calibration payload did not come back off the device, so `upload_ms` is timing a \
         host memcpy and not an upload"
    );
    let calib_gb_per_s = committed.bytes / calib_ms / 1e6;
    assert!(
        calib_gb_per_s < PCIE4_X16_GB_PER_S,
        "{calib_gb_per_s:.4} GB/s beats PCIe 4.0 x16's {PCIE4_X16_GB_PER_S} GB/s, so the \
         transfer did not happen and every upload column is a zero that could not have been \
         non-zero"
    );
    assert!(
        calib_gb_per_s > committed_gb_per_s,
        "{calib_gb_per_s:.4} GB/s is slower than gpu_vs_cpu.csv's own effective \
         {committed_gb_per_s:.4} GB/s, which carries host field evaluation this path does not -- \
         so this harness is measuring something other than the bus"
    );
    println!(
        "vacuity control: {calib_bytes} bytes uploaded in {calib_ms:.4} ms = \
         {calib_gb_per_s:.4} GB/s, read back and byte-identical ({roundtrip}); \
         under the {PCIE4_X16_GB_PER_S} GB/s bus ceiling and above the committed \
         {committed_gb_per_s:.4} GB/s\n"
    );

    let khz = cpu_khz();
    let mut rows: Vec<Row> = Vec::new();
    let mut verdicts: Vec<Verdict> = Vec::new();

    for field_name in ["gyroid", "fbm_terrain"] {
        let gyroid = Gyroid::<f64>::canonical();
        let fbm = FbmTerrain::<f64>::canonical();
        let centres = if field_name == "gyroid" {
            dig_path(&gyroid)
        } else {
            dig_path(&fbm)
        };

        let mut arms: Vec<Arm> = Vec::new();
        for c in GRANULARITIES {
            let arm = if field_name == "gyroid" {
                run_arm(&gyroid, c, &uploader, &centres)
            } else {
                run_arm(&fbm, c, &uploader, &centres)
            };
            arms.push(arm);
        }
        let mixed = if field_name == "gyroid" {
            run_mixed(&gyroid, 64, 4, &uploader, &centres)
        } else {
            run_mixed(&fbm, 64, 4, &uploader, &centres)
        };
        // The registered pair is 64³/4³ and that is what C3 is scored on. This
        // second one exists because the registered pair may leave no static bulk
        // at all -- a 128³ world holds only EIGHT 64³ chunks -- and a
        // falsification that turned out to be about the fixture's chunk count
        // rather than about mixed granularity would be worthless.
        let mixed16 = if field_name == "gyroid" {
            run_mixed(&gyroid, 16, 4, &uploader, &centres)
        } else {
            run_mixed(&fbm, 16, 4, &uploader, &centres)
        };

        // ── control 1: every arm meshed something ────────────────────────────
        for a in &arms {
            assert!(
                a.dirty_chunks > 0,
                "VOID: chunk {} marked no dirty chunk in {EDITS} edits on {field_name}",
                a.chunk_cells
            );
            assert!(
                a.world.raw_vertices > 0,
                "VOID: chunk {} produced no geometry on {field_name}",
                a.chunk_cells
            );
            assert!(
                a.world.payload.len() as u64 > 0 && a.full_upload_ms > 0.0,
                "VOID: chunk {} uploaded nothing, so its upload column is not a measurement",
                a.chunk_cells
            );
            assert!(
                a.stream.upload_bytes > 0 && a.stream.upload_ms > 0.0,
                "VOID: chunk {} streamed no bytes on {field_name}, so the camera arm's cost is \
                 not an upload",
                a.chunk_cells
            );
            // C3's registered control, on every arm rather than only the mixed
            // one: a camera whose set only ever grew would be a frustum opening.
            assert!(
                a.stream.enters > 0 && a.stream.exits > 0,
                "VOID: chunk {} saw {} enters and {} exits on {field_name}; the camera did not \
                 move enough to change the visible chunk set",
                a.chunk_cells,
                a.stream.enters,
                a.stream.exits
            );
            // C2's own vacuity guard: a weld that removed nothing prices nothing.
            assert!(
                a.welded_removed > 0,
                "VOID: chunk {} welded {} vertices away, so C2 is comparing a weld cost against \
                 no saving at all",
                a.chunk_cells,
                a.welded_removed
            );
        }
        for m in [&mixed, &mixed16] {
            assert!(
                m.activations > 0,
                "VOID: the {}³/{}³ mixed arm subdivided no coarse chunk, so it IS the coarse arm \
                 and C3 is about nothing",
                m.coarse_cells,
                m.fine_cells
            );
            assert!(
                m.stream.enters > 0 && m.stream.exits > 0,
                "VOID: the {}³/{}³ mixed arm's camera saw {} enters and {} exits",
                m.coarse_cells,
                m.fine_cells,
                m.stream.enters,
                m.stream.exits
            );
        }

        // ── control 2: every arm produced the same surface ───────────────────
        let reference = &arms[0].world.surface;
        for a in &arms[1..] {
            let differing = a.world.surface.symmetric_difference(reference).count();
            assert_eq!(
                differing, 0,
                "chunk {} disagrees with chunk {} on {differing} quantised surface points of {} \
                 on {field_name}: a partition change moved the surface",
                a.chunk_cells,
                arms[0].chunk_cells,
                reference.len()
            );
        }
        for m in [&mixed, &mixed16] {
            let differing = m.world.surface.symmetric_difference(reference).count();
            assert_eq!(
                differing, 0,
                "the {}³/{}³ mixed arm disagrees with chunk {} on {differing} quantised surface \
                 points of {} on {field_name}: a two-level partition dropped a seam, which is a \
                 defect and not a speed result",
                m.coarse_cells,
                m.fine_cells,
                arms[0].chunk_cells,
                reference.len()
            );
        }

        // ── control 3: every field evaluation is accounted for ───────────────
        let baseline = u64::from(WORLD_CELLS).pow(3) as f64;
        let mut stencils: Vec<u64> = Vec::new();
        for a in &arms {
            let grid = a.chunks * a.samples_per_chunk;
            assert!(
                a.field_calls >= grid,
                "chunk {}: {} field calls is fewer than the {grid} sample-grid points",
                a.chunk_cells,
                a.field_calls
            );
            let normals = a.field_calls - grid;
            assert_eq!(
                normals % a.build_vertices,
                0,
                "chunk {}: {normals} non-grid field calls over {} vertices is not a whole stencil",
                a.chunk_cells,
                a.build_vertices
            );
            stencils.push(normals / a.build_vertices);
            let c = f64::from(a.chunk_cells);
            let predicted = ((c + 1.0) / c).powi(3);
            let measured = grid as f64 / baseline;
            assert!(
                (measured - predicted).abs() < 1e-4,
                "chunk {}: measured GRID duplication {measured:.6} against the registered \
                 ((c+1)/c)^3 = {predicted:.6}",
                a.chunk_cells
            );
        }

        // ── C1: the crossover, solved rather than searched ───────────────────
        //
        // Both totals are affine in the rate, so the crossing is one division
        // and a bisection over the sampled rates would only re-derive it less
        // exactly.
        let a4 = arms.iter().find(|a| a.chunk_cells == 4).expect("a 4 arm");
        let a64 = arms.iter().find(|a| a.chunk_cells == 64).expect("a 64 arm");
        let solve = |d_stream: f64, d_edit: f64| -> f64 {
            if d_edit > 0.0 && d_stream > 0.0 {
                d_stream / d_edit
            } else {
                f64::NAN
            }
        };
        let d_edit = a64.edit_ms() - a4.edit_ms();
        let d_stream = a4.stream_ms_per_second() - a64.stream_ms_per_second();
        let crossover = solve(d_stream, d_edit);
        let c1 = crossover.is_finite() && (RATE_MIN..=RATE_MAX).contains(&crossover);

        // The same crossing on the arena upload path, where the cost is one
        // `write_buffer` for the frame instead of one per visible chunk.
        let d_edit_arena = a64.edit_ms_arena() - a4.edit_ms_arena();
        let d_stream_arena = a4.stream_arena_ms_per_second() - a64.stream_arena_ms_per_second();
        let crossover_arena = solve(d_stream_arena, d_edit_arena);
        let c1_arena =
            crossover_arena.is_finite() && (RATE_MIN..=RATE_MAX).contains(&crossover_arena);

        // How much of each arm's streaming cost is BYTES, from the exact byte
        // counts at the calibrated rate.
        let bytes_ms =
            |a: &Arm| (a.stream.upload_bytes as f64 / (calib_gb_per_s * 1e6)) / CAMERA_SECONDS;

        // **The crossing that answers the clause as registered, and the only one
        // of the three that is not a clock.** `M-377`'s bill is 51.6% more vertex
        // DATA, so the quantity 64³ can win on is bytes. The per-chunk path's
        // crossover is contaminated by ~2 µs of driver cost per `write_buffer`,
        // which is a call count; the arena path's is contaminated the other way,
        // by a per-frame `submit` + `poll` that the 64³ arm pays across only
        // seven frames and the 4³ arm across thirty-two — a seven-sample
        // sub-millisecond measurement that came back with the WRONG SIGN on
        // `gyroid`. This one is exact integers over one measured rate, which is
        // the substitution `M-337`'s re-audit forced: gate on counts.
        let d_stream_bytes = bytes_ms(a4) - bytes_ms(a64);
        let crossover_bytes = solve(d_stream_bytes, d_edit);
        let c1_bytes =
            crossover_bytes.is_finite() && (RATE_MIN..=RATE_MAX).contains(&crossover_bytes);

        // What the streaming term would have to be for C1 to be reachable, in
        // units a reader can picture: full re-uploads of the world's EXCESS
        // vertex data per second, against what this camera actually supplies.
        let excess_bytes = a4.world.payload.len() as f64 - a64.world.payload.len() as f64;
        let excess_ms = excess_bytes / (calib_gb_per_s * 1e6);
        let needed = RATE_MIN * d_edit / excess_ms;
        let supplied = |a: &Arm| {
            (a.stream.upload_bytes as f64 / a.world.payload.len() as f64) / CAMERA_SECONDS
        };
        let shortfall = needed / supplied(a4);

        // ── C2: the weld against the upload it saves, at 4³ ──────────────────
        let c2 = a4.weld_ms > a4.weld_saving_ms(calib_gb_per_s);

        // ── C3: mixed against both fixed choices, at game_dig's own rate ─────
        let best_fixed = |r: f64| a4.total_ms(r).min(a64.total_ms(r));
        let mixed_speedup_at = |r: f64| best_fixed(r) / mixed.total_ms(r);
        let mixed_steady_speedup_at = |r: f64| best_fixed(r) / mixed.steady_total_ms(r);
        let c3 = mixed_speedup_at(GAME_RATE) >= 1.5;

        println!(
            "{field_name}: {} distinct surface points, normal stencils {stencils:?}",
            reference.len()
        );
        println!(
            "{:>6} {:>8} {:>9} {:>9} {:>8} {:>9} {:>9} {:>8} {:>7} {:>9} {:>8} {:>8}",
            "chunk",
            "mark/ed",
            "remesh/ed",
            "upload/ed",
            "edit/ed",
            "stream/s",
            "arena/s",
            "bytes/s",
            "writes",
            "verts",
            "weld",
            "saves"
        );
        for a in &arms {
            println!(
                "{:>6} {:>8.4} {:>9.4} {:>9.4} {:>8.4} {:>9.4} {:>9.4} {:>8.4} {:>7} {:>9} \
                 {:>8.2} {:>8.4}",
                a.chunk_cells,
                a.mark_ms / EDITS as f64,
                a.remesh_ms / EDITS as f64,
                a.upload_dirty_ms / EDITS as f64,
                a.edit_ms(),
                a.stream_ms_per_second(),
                a.stream_arena_ms_per_second(),
                bytes_ms(a),
                a.stream.writes,
                a.world.raw_vertices,
                a.weld_ms,
                a.weld_saving_ms(calib_gb_per_s)
            );
        }
        for m in [&mixed, &mixed16] {
            println!(
                "{:>6} {:>8.4} {:>9.4} {:>9.4} {:>8.4} {:>9.4} {:>9.4} {:>8} {:>7} {:>9} \
                 {:>8} {:>8}",
                format!("{}/{}", m.coarse_cells, m.fine_cells),
                m.mark_ms / EDITS as f64,
                m.remesh_ms / EDITS as f64,
                (m.subdivide_ms + m.subdivide_upload_ms + m.dirty_upload_ms) / EDITS as f64,
                m.edit_ms(),
                m.stream_ms_per_second(),
                m.stream_arena_ms_per_second(),
                format!("{}/{}", m.activations, m.coarse_chunks),
                m.stream.writes,
                m.world.raw_vertices,
                format!("{:.4}", m.steady_edit_ms()),
                "steady"
            );
        }
        println!(
            "  C1 per-chunk upload: d_stream {d_stream:.6} ms/s over d_edit {d_edit:.4} \
             ms/edit = {crossover:.6} edits/s -> {}",
            if c1 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "  C1 arena upload (bytes only): d_stream {d_stream_arena:.6} ms/s over d_edit \
             {d_edit_arena:.4} ms/edit = {crossover_arena:.6} edits/s -> {}",
            if c1_arena { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "  C1 exact bytes over one calibrated rate: d_stream {d_stream_bytes:.6} ms/s over \
             d_edit {d_edit:.4} ms/edit = {crossover_bytes:.6} edits/s -> {}",
            if c1_bytes { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "     the world's {excess_bytes:.0} EXCESS vertex bytes cost {excess_ms:.4} ms at \
             the calibrated {calib_gb_per_s:.4} GB/s, so they would have to be re-uploaded \
             {needed:.2} times per second to put the crossover at {RATE_MIN}/s; this camera \
             supplies {:.2} full-world uploads per second, short by {shortfall:.1}x",
            supplied(a4)
        );
        println!(
            "  C2 at 4³: weld {:.4} ms against {:.4} ms of upload saved (measured difference \
             {:.4} ms) -> weld costs {}",
            a4.weld_ms,
            a4.weld_saving_ms(calib_gb_per_s),
            a4.weld_saving_ms_measured(),
            if c2 { "MORE" } else { "LESS" }
        );
        println!(
            "  C3 at {GAME_RATE}/s: mixed {:.4} ms/s against best fixed {:.4} ms/s = {:.4}x \
             (steady state, transition amortised away: {:.4}x)",
            mixed.total_ms(GAME_RATE),
            best_fixed(GAME_RATE),
            mixed_speedup_at(GAME_RATE),
            mixed_steady_speedup_at(GAME_RATE)
        );
        println!(
            "     the fairer coarse level, {}³/{}³ with {}/{} activated: {:.4}x at \
             {GAME_RATE}/s, steady {:.4}x",
            mixed16.coarse_cells,
            mixed16.fine_cells,
            mixed16.activations,
            mixed16.coarse_chunks,
            best_fixed(GAME_RATE) / mixed16.total_ms(GAME_RATE),
            best_fixed(GAME_RATE) / mixed16.steady_total_ms(GAME_RATE)
        );
        println!();

        verdicts.push(Verdict {
            field: field_name,
            per_chunk: crossover,
            arena: crossover_arena,
            bytes: crossover_bytes,
            c1,
            c1_arena,
            c1_bytes,
            c2,
            c3,
        });

        for a in &arms {
            for rate in RATES {
                let remesh_ms = rate * (a.mark_ms + a.remesh_ms) / EDITS as f64;
                let upload_ms =
                    rate * a.upload_dirty_ms / EDITS as f64 + a.stream_ms_per_second();
                let total_ms = remesh_ms + upload_ms;
                rows.push(vec![
                    ("field", field_name.to_string()),
                    ("chunk_cells", a.chunk_cells.to_string()),
                    ("edits_per_second", format!("{rate}")),
                    ("remesh_ms", format!("{remesh_ms:.6}")),
                    ("upload_ms", format!("{upload_ms:.6}")),
                    ("weld_ms", format!("{:.6}", a.weld_ms)),
                    ("total_ms", format!("{total_ms:.6}")),
                    ("raw_vertices", a.world.raw_vertices.to_string()),
                    ("distinct_surface_points", reference.len().to_string()),
                    (
                        "vertex_duplication",
                        format!(
                            "{:.6}",
                            a.world.raw_vertices as f64 / reference.len() as f64
                        ),
                    ),
                    (
                        "crossover_edits_per_second",
                        if crossover.is_finite() {
                            format!("{crossover:.6}")
                        } else {
                            "none".to_string()
                        },
                    ),
                    (
                        "mixed_granularity_ms",
                        format!("{:.6}", mixed.total_ms(rate)),
                    ),
                    ("mixed_speedup", format!("{:.6}", mixed_speedup_at(rate))),
                    (
                        "visible_chunk_changes",
                        (a.stream.enters + a.stream.exits).to_string(),
                    ),
                    ("c1_holds", c1.to_string()),
                    ("c2_holds", c2.to_string()),
                    ("c3_holds", c3.to_string()),
                    // ── extras ──────────────────────────────────────────────
                    ("world_cells", WORLD_CELLS.to_string()),
                    ("chunks", a.chunks.to_string()),
                    ("edits", EDITS.to_string()),
                    ("triangles", a.world.triangles.to_string()),
                    ("mark_ms_per_edit", format!("{:.6}", a.mark_ms / EDITS as f64)),
                    (
                        "remesh_ms_per_edit",
                        format!("{:.6}", a.remesh_ms / EDITS as f64),
                    ),
                    (
                        "upload_dirty_ms_per_edit",
                        format!("{:.6}", a.upload_dirty_ms / EDITS as f64),
                    ),
                    (
                        "upload_dirty_bytes_per_edit",
                        (a.upload_dirty_bytes / EDITS as u64).to_string(),
                    ),
                    ("edit_ms_per_edit", format!("{:.6}", a.edit_ms())),
                    (
                        "stream_upload_ms_per_second",
                        format!("{:.6}", a.stream_ms_per_second()),
                    ),
                    (
                        "stream_upload_bytes_per_second",
                        ((a.stream.upload_bytes as f64 / CAMERA_SECONDS) as u64).to_string(),
                    ),
                    ("stream_units", a.stream.units.to_string()),
                    (
                        "stream_frames_uploading",
                        a.stream.frames_uploading.to_string(),
                    ),
                    ("visible_chunk_enters", a.stream.enters.to_string()),
                    ("visible_chunk_exits", a.stream.exits.to_string()),
                    (
                        "full_world_upload_bytes",
                        a.world.payload.len().to_string(),
                    ),
                    ("full_world_upload_ms", format!("{:.6}", a.full_upload_ms)),
                    ("welded_vertices", a.welded_vertices.to_string()),
                    ("welded_vertices_removed", a.welded_removed.to_string()),
                    (
                        "welded_upload_bytes",
                        a.welded_payload_bytes.to_string(),
                    ),
                    (
                        "welded_upload_ms",
                        format!("{:.6}", a.welded_upload_ms),
                    ),
                    (
                        "weld_saves_upload_ms",
                        format!("{:.6}", a.weld_saving_ms(calib_gb_per_s)),
                    ),
                    (
                        "weld_saves_upload_ms_measured",
                        format!("{:.6}", a.weld_saving_ms_measured()),
                    ),
                    (
                        "weld_pays",
                        (a.weld_ms < a.weld_saving_ms(calib_gb_per_s)).to_string(),
                    ),
                    (
                        "upload_share_of_total",
                        format!("{:.6}", 100.0 * upload_ms / total_ms),
                    ),
                    (
                        "excess_vertex_bytes_vs_64",
                        format!("{excess_bytes:.0}"),
                    ),
                    ("excess_vertex_bytes_ms", format!("{excess_ms:.6}")),
                    (
                        "excess_reuploads_per_second_for_crossover",
                        format!("{needed:.4}"),
                    ),
                    // The decomposition the first run could not produce: how much
                    // of the streaming cost is bytes and how much is one
                    // `write_buffer` per visible chunk.
                    (
                        "stream_writes",
                        a.stream.writes.to_string(),
                    ),
                    (
                        "stream_bytes_ms_per_second",
                        format!("{:.6}", bytes_ms(a)),
                    ),
                    (
                        "stream_overhead_ms_per_second",
                        format!("{:.6}", a.stream_ms_per_second() - bytes_ms(a)),
                    ),
                    (
                        "stream_arena_ms_per_second",
                        format!("{:.6}", a.stream_arena_ms_per_second()),
                    ),
                    (
                        "upload_arena_ms",
                        format!(
                            "{:.6}",
                            rate * a.upload_dirty_arena_ms / EDITS as f64
                                + a.stream_arena_ms_per_second()
                        ),
                    ),
                    (
                        "total_arena_ms",
                        format!("{:.6}", a.total_arena_ms(rate)),
                    ),
                    (
                        "crossover_arena_edits_per_second",
                        if crossover_arena.is_finite() {
                            format!("{crossover_arena:.6}")
                        } else {
                            "none".to_string()
                        },
                    ),
                    ("c1_holds_arena", c1_arena.to_string()),
                    (
                        "crossover_bytes_edits_per_second",
                        if crossover_bytes.is_finite() {
                            format!("{crossover_bytes:.6}")
                        } else {
                            "none".to_string()
                        },
                    ),
                    ("c1_holds_bytes", c1_bytes.to_string()),
                    (
                        "stream_full_world_uploads_per_second",
                        format!("{:.6}", supplied(a)),
                    ),
                    (
                        "crossover_streaming_shortfall",
                        format!("{shortfall:.4}"),
                    ),
                    ("grid_duplication", {
                        let c = f64::from(a.chunk_cells);
                        format!("{:.6}", ((c + 1.0) / c).powi(3))
                    }),
                    ("field_calls", a.field_calls.to_string()),
                    (
                        "normal_stencil",
                        ((a.field_calls - a.chunks * a.samples_per_chunk) / a.build_vertices)
                            .to_string(),
                    ),
                    ("mixed_coarse_cells", mixed.coarse_cells.to_string()),
                    ("mixed_fine_cells", mixed.fine_cells.to_string()),
                    ("mixed_coarse_activations", mixed.activations.to_string()),
                    ("mixed_coarse_chunks", mixed.coarse_chunks.to_string()),
                    ("mixed_edit_ms_per_edit", format!("{:.6}", mixed.edit_ms())),
                    (
                        "mixed_steady_edit_ms_per_edit",
                        format!("{:.6}", mixed.steady_edit_ms()),
                    ),
                    (
                        "mixed_steady_granularity_ms",
                        format!("{:.6}", mixed.steady_total_ms(rate)),
                    ),
                    (
                        "mixed_steady_speedup",
                        format!("{:.6}", mixed_steady_speedup_at(rate)),
                    ),
                    (
                        "mixed_subdivide_ms_total",
                        format!("{:.6}", mixed.subdivide_ms),
                    ),
                    (
                        "mixed_subdivide_upload_ms_total",
                        format!("{:.6}", mixed.subdivide_upload_ms),
                    ),
                    (
                        "mixed_arena_granularity_ms",
                        format!("{:.6}", mixed.total_arena_ms(rate)),
                    ),
                    (
                        "mixed_stream_ms_per_second",
                        format!("{:.6}", mixed.stream_ms_per_second()),
                    ),
                    ("mixed_raw_vertices", mixed.world.raw_vertices.to_string()),
                    (
                        "mixed_visible_chunk_changes",
                        (mixed.stream.enters + mixed.stream.exits).to_string(),
                    ),
                    (
                        "mixed_speedup_vs_4",
                        format!("{:.6}", a4.total_ms(rate) / mixed.total_ms(rate)),
                    ),
                    (
                        "mixed_speedup_vs_64",
                        format!("{:.6}", a64.total_ms(rate) / mixed.total_ms(rate)),
                    ),
                    (
                        "mixed16_coarse_activations",
                        mixed16.activations.to_string(),
                    ),
                    (
                        "mixed16_coarse_cells",
                        mixed16.coarse_cells.to_string(),
                    ),
                    (
                        "mixed16_coarse_chunks",
                        mixed16.coarse_chunks.to_string(),
                    ),
                    (
                        "mixed16_granularity_ms",
                        format!("{:.6}", mixed16.total_ms(rate)),
                    ),
                    (
                        "mixed16_speedup",
                        format!("{:.6}", best_fixed(rate) / mixed16.total_ms(rate)),
                    ),
                    (
                        "mixed16_steady_speedup",
                        format!("{:.6}", best_fixed(rate) / mixed16.steady_total_ms(rate)),
                    ),
                    ("committed_upload_ms", format!("{:.6}", committed.upload_ms)),
                    ("committed_upload_bytes", (committed.bytes as u64).to_string()),
                    (
                        "committed_upload_gb_per_s",
                        format!("{committed_gb_per_s:.6}"),
                    ),
                    ("calibration_upload_bytes", calib_bytes.to_string()),
                    ("calibration_upload_ms", format!("{calib_ms:.6}")),
                    (
                        "calibration_upload_gb_per_s",
                        format!("{calib_gb_per_s:.6}"),
                    ),
                    ("upload_roundtrip_verified", roundtrip.to_string()),
                    ("bus_ceiling_gb_per_s", format!("{PCIE4_X16_GB_PER_S}")),
                    ("adapter", report.name.replace([',', ' '], "_")),
                    ("backend", format!("{:?}", report.backend)),
                    ("cpu_khz", khz.clone()),
                ]);
            }
        }
    }

    println!("verdicts:");
    for v in &verdicts {
        println!(
            "  {}: C1 crossover edits/s -- per-chunk {:.6} {}, arena {:.6} {}, exact bytes \
             {:.6} {}; C2 {}; C3 {}",
            v.field,
            v.per_chunk,
            if v.c1 { "HELD" } else { "FALSIFIED" },
            v.arena,
            if v.c1_arena { "HELD" } else { "FALSIFIED" },
            v.bytes,
            if v.c1_bytes { "HELD" } else { "FALSIFIED" },
            if v.c2 { "HELD" } else { "FALSIFIED" },
            if v.c3 { "HELD" } else { "FALSIFIED" }
        );
    }

    common::experiment::run(prereg, |run| {
        for row in &rows {
            run.record(row);
        }
    });
}
