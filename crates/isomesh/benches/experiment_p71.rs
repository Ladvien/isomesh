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
//! | `budget` | **C3** | `DeferredGeometry` under `DirtySet::mesh_within_budget`, on `M-124`'s fixture |
//!
//! `extract_buffers` and `extract_indirect` both already exist, so C2 is a
//! measurement of shipped code rather than a build — which is why the comparison
//! can be within one binary and one run, as `M-281` requires.
//!
//! # C3's fixture is `M-124`'s, and that is the whole reason this arm was rewritten
//!
//! The registration says the clause fires on *"exactly the rows `M-124` has"*, and
//! `M-124` is a **budget** sweep: 288 chunks under
//! [`DirtySet::mesh_within_budget`](isomesh::chunk::dirty::DirtySet::mesh_within_budget),
//! 25 µs to 8 ms, 2,360 frames each. The first version of this arm swept
//! **resolution** over a bench-local staging ring, which measured a real
//! mechanism and could not score the clause — `c3_holds` read `not_measured` on
//! all twelve rows of `docs/experiments/p-71.csv` (`M-376`).
//!
//! Two things had to change for the clause to be answerable, and only one of them
//! is in this file. `DeferredGeometry` is now a **public type** in `isomesh-gpu`
//! (`R-071`), because the clause is about the queue *under a scheduler* and a
//! queue that exists only as bench scaffolding cannot be put under one. The
//! bench-local `StagingRing` is gone: two queues in one tree is the one-path
//! defect, and the shipped one is the one the clause is about. Its own rows are
//! historical and live in this CSV's git history.
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

use std::cell::RefCell;
use std::collections::HashMap;
use std::time::{Duration, Instant};

use isomesh::Sdf;
use isomesh::chunk::dirty::DirtySet;
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::{ReferenceField, Sphere};
use isomesh_gpu::headless::Gpu;
use isomesh_gpu::{
    DeferredGeometry, FieldBuffer, GridParams, MarchingCubesGpu, read_bytes_many_deferred,
};

/// The resolutions C1 is attributed over. 129³ is the size every number in
/// `M-149`, `M-150`, `M-159` and `M-167` is quoted at, so it is the row the
/// clause is about; the smaller three are what make the trend visible rather
/// than a single point.
const SIZES: [u32; 4] = [33, 65, 97, 129];

/// Repetitions per measurement. The **median** is reported: a GPU submission
/// shares the device with a compositor, and a mean would let one scheduling
/// hiccup become the figure.
const REPS: usize = 7;

/// ── C3's fixture, which is `M-124`'s ────────────────────────────────────────
///
/// Every constant in this block is read from
/// `bevy_isomesh/examples/game_budget.rs:58-74` rather than chosen here, because
/// the registration scores the clause on `M-124`'s own rows and a fixture that
/// merely resembles it would score a different experiment. 16 cells a chunk.
const CHUNK_CELLS: u32 = 16;
/// 0.25 world units a cell, so a chunk is 4 units on a side.
const CELL_SIZE: f32 = 0.25;
/// Chunks along x and z.
const SPAN: i32 = 12;
/// Layers in y: `game_budget`'s `LAYERS = 0..=1`.
const LAYERS: i32 = 2;
/// 12 × 2 × 12 = 288, `M-124`'s chunk count.
const CHUNKS: usize = (SPAN * LAYERS * SPAN) as usize;
/// `game_budget`'s `BUDGETS_US`. A **320× range**, which is the range the clause
/// names, and its low end is deliberately below one chunk's cost.
const BUDGETS_US: [u64; 8] = [25, 50, 200, 500, 1_000, 2_000, 4_000, 8_000];
/// Frames per `(budget, delay)` cell, `M-124`'s count.
const BUDGET_FRAMES: usize = 2_360;
/// Frames of read-back latency the queue is sized to absorb, swept rather than
/// fixed.
///
/// **A frames column, not a slots column**, which is what `ring_frames_delay`
/// has always meant here: the registration asks about *"an N-frame-delayed
/// double-buffered staging ring"* costing *"one to two frames of collision
/// latency"*, and the old ring arm's `DEPTH = 2` was two frames.
///
/// The distinction is not pedantry, it was measured. A flat slot count is a
/// **second budget**, not a latency knob: at a fixed capacity of 2 the queue held
/// one un-ready slot at every drain — the read-back submitted last in a pass is
/// microseconds old when the next frame looks at it — so the pass meshed
/// **exactly one chunk at every budget from 25 µs to 8 ms**, and `c3_holds` would
/// have been a measurement of `budget / chunk_ms` against the number 2 rather
/// than of the queue under a scheduler. [`capacity_for`] derives slots from the
/// delay and the budget, so what is swept is the latency the clause names.
const DELAYS: [usize; 3] = [1, 2, 4];

/// Slots for `delay` frames of a budget's worth of chunks, plus the frame
/// currently submitting.
///
/// `budget / chunk_ms` is what a pass buys, and `chunk_ms` is measured in this
/// same run ([`one_chunk_ms`]) rather than assumed — `M-281`. The `+ 1` frame of
/// headroom is the un-ready tail described on [`DELAYS`]: without it the queue
/// clips the last chunk of every pass.
fn capacity_for(delay: usize, budget_ms: f64, chunk_ms: f64) -> usize {
    let per_pass = (budget_ms / chunk_ms).ceil().max(1.0) as usize;
    (delay + 1) * per_pass
}

/// `game_budget`'s field, reproduced for the same reason as the constants: the
/// per-chunk triangle count is the size of the read-back the queue carries, and
/// that is a property of this surface.
///
/// A ground sheet unioned with a sphere at the world corner. Half the chunks —
/// the whole upper layer — hold no surface at all, which is `M-124`'s population
/// and not a defect: a scheduler's cost per chunk is a *mean* over a real world,
/// and `one_chunk_ms` below is measured over the same 288.
#[derive(Clone, Copy)]
struct Blobs;

impl Sdf for Blobs {
    type Scalar = f32;
    fn sample(&self, p: [f32; 3]) -> f32 {
        let ground = p[1] - 1.4 * (p[0] * 0.35).sin() * (p[2] * 0.31).cos() - 2.0;
        let r = (p[0] * p[0] + (p[1] - 2.0) * (p[1] - 2.0) + p[2] * p[2]).sqrt();
        ground.min(r - 6.0)
    }
}

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

/// What one `(budget, delay)` cell of C3's sweep measured.
struct BudgetCell {
    budget_us: u64,
    /// Frames of latency the queue was sized for. The swept column.
    delay: usize,
    /// Slots that came out of [`capacity_for`].
    capacity: usize,
    /// Frames that ran a budgeted pass. Fewer than [`BUDGET_FRAMES`] when the
    /// queue had no room, which is the latency `no_room_frames` reports.
    passes: usize,
    chunks: usize,
    /// Mean wall time of one budgeted pass, in milliseconds.
    ///
    /// **The pass, not the frame, and that is `M-124`'s quantity.** `spend`
    /// bounds `mesh_within_budget` and nothing else, so `game_budget` measures
    /// the overshoot of the budgeted pass rather than the whole Bevy frame. A
    /// mean over whole frames here would be a mean over a harness with no
    /// rendering in it, which is not a number about anything.
    mean_ms: f64,
    mean_chunks: f64,
    collected: usize,
    max_in_flight: usize,
    /// Frames where the queue was still full after the drain, so no pass ran.
    ///
    /// The honest cost of a shallow queue: a stall the queue moved rather than
    /// removed. `M-376`'s `not_ready` under a different name.
    no_room_frames: usize,
    /// Frames between a read-back's submission and its collection, summed over
    /// every collection made **during** the run, how many such collections there
    /// were, and the worst single wait.
    ///
    /// **The registration's owner-question, measured instead of assumed.** It
    /// says the queue costs *"one to two frames of collision latency"*; these are
    /// whether that is true at each delay. The tail drain's collections are
    /// excluded: a poll after the last frame is not a frame.
    latency_frames_total: usize,
    latency_samples: usize,
    max_latency_frames: usize,
    /// Frames where the dirty set was empty. Must be zero — a pass over an empty
    /// set measures an idle scheduler.
    empty_set_frames: usize,
    /// Read-backs still in flight when the cell ended, drained bounded.
    drain_frames: usize,
}

impl BudgetCell {
    fn budget_ms(&self) -> f64 {
        self.budget_us as f64 / 1000.0
    }

    /// Mean frames a read-back waited between submission and collection.
    fn mean_latency_frames(&self) -> f64 {
        if self.latency_samples == 0 {
            f64::NAN
        } else {
            self.latency_frames_total as f64 / self.latency_samples as f64
        }
    }

    /// What the pass actually spent per chunk it meshed.
    ///
    /// Reported because it is not flat, and the shape is `M-159`'s mechanism at
    /// pass granularity: `extract_buffers` waits with `poll(Wait)` and **no
    /// submission index**, which drains every dispatch queued before it — so a
    /// pass's *first* chunk pays for the previous pass's outstanding deferred
    /// copies and every later chunk does not. A one-chunk pass is all first
    /// chunk; a thirty-chunk pass amortises it thirty ways.
    fn ms_per_chunk(&self) -> f64 {
        self.mean_ms / self.mean_chunks
    }

    /// `|mean − budget|` in units of one chunk. **The number to read, rather than
    /// the boolean below.**
    ///
    /// `mesh_within_budget`'s own doc prices the never-livelock guarantee at
    /// *"overshooting by at most one chunk"*, so this ratio sits at 1.0 by design
    /// and a boolean at that boundary is a coin flip between runs — measured: the
    /// same fixture put 4 cells outside on one run and 9 on the next while every
    /// ratio stayed in the same narrow band. The ratio is the stable quantity and
    /// the entry quotes it.
    fn overshoot_chunks(&self, one_chunk_ms: f64) -> f64 {
        (self.mean_ms - self.budget_ms()).abs() / one_chunk_ms
    }

    /// C3's acceptance, which is `M-124`'s property restated: the pass tracks the
    /// budget to within one chunk.
    ///
    /// One chunk is the tolerance rather than a percentage because
    /// `mesh_within_budget` consults `spend` **after** each chunk (`dirty.rs`
    /// 105–112), so a budget below one chunk's cost overshoots by exactly one
    /// chunk *by design* and the guarantee's own price has to be inside the bar.
    fn within_one_chunk(&self, one_chunk_ms: f64) -> bool {
        self.overshoot_chunks(one_chunk_ms) <= 1.0
    }
}

/// `M-124`'s world: the chunk lattice, one resident field buffer per chunk, and
/// the camera the scheduler orders against.
///
/// The fields are built once and held, out of the frame loop deliberately: what
/// C3 is about is the **extract plus submit**, and re-uploading a chunk's samples
/// every frame would put a `write_buffer` of 19,652 bytes inside the budget the
/// clause is scored against. A game holds its chunk fields resident for the same
/// reason.
struct World {
    layout: ChunkLayout<f32>,
    fields: HashMap<ChunkId, FieldBuffer>,
    /// The world centre. The ordering is a pure function of the set's contents
    /// and the camera, so any fixed point gives a reproducible order; the centre
    /// is the one that makes nearest-first symmetric over the 12 × 2 × 12 block
    /// instead of sweeping it from a corner.
    camera: [f32; 3],
}

impl World {
    fn build(gpu: &Gpu) -> Self {
        let layout = ChunkLayout::new(CHUNK_CELLS, CELL_SIZE, [0.0f32; 3]).expect("M-124's layout");
        let mut fields = HashMap::with_capacity(CHUNKS);
        for x in 0..SPAN {
            for y in 0..LAYERS {
                for z in 0..SPAN {
                    let id = ChunkId::new([x, y, z]);
                    let params =
                        GridParams::new([CHUNK_CELLS + 1; 3], layout.sample_origin(id), CELL_SIZE)
                            .expect("a chunk grid");
                    let field = FieldBuffer::sampled(gpu.device(), gpu.queue(), params, &Blobs)
                        .expect("sample a chunk");
                    fields.insert(id, field);
                }
            }
        }
        assert_eq!(fields.len(), CHUNKS, "M-124's fixture is 288 chunks");

        let extent = SPAN as f32 * CHUNK_CELLS as f32 * CELL_SIZE;
        let height = LAYERS as f32 * CHUNK_CELLS as f32 * CELL_SIZE;
        Self {
            layout,
            fields,
            camera: [extent * 0.5, height * 0.5, extent * 0.5],
        }
    }

    fn field(&self, id: ChunkId) -> &FieldBuffer {
        self.fields.get(&id).expect("every chunk has a field")
    }

    /// A [`DirtySet`] holding every chunk.
    fn all_dirty(&self) -> DirtySet {
        let mut dirty = DirtySet::new();
        for id in self.fields.keys() {
            dirty.insert(*id);
        }
        dirty
    }

    /// The order `mesh_within_budget` meshes in, taken from the scheduler itself.
    ///
    /// Nearest-first by squared distance with ties on [`ChunkId`], so it is a pure
    /// function of the set's contents and the camera. Read from the crate rather
    /// than recomputed here: a second copy of the ordering rule would be a second
    /// thing to drift, and this also asserts all 288 are reachable in one
    /// unbounded pass.
    fn scheduler_order(&self) -> Vec<ChunkId> {
        let mut probe = self.all_dirty();
        let mut order = Vec::with_capacity(CHUNKS);
        let report =
            probe.mesh_within_budget(&self.layout, self.camera, |id, _| order.push(id), || true);
        assert_eq!(report.meshed, CHUNKS);
        assert!(report.is_drained());
        order
    }
}

/// Extract one chunk on the GPU and hand its geometry read-back to `queue`.
///
/// The unit C3 budgets. `extract_buffers` still blocks once, for the four bytes
/// of the triangle count — the queue removes the **geometry** copy, which
/// `M-376` measured at 0.7075 ms of a 1.1860 ms extraction at 129³, not the
/// 0.3177 ms count wait. A fully non-blocking path needs `extract_indirect`'s
/// device-side totals as well, which is outside `P-71`.
fn extract_and_submit(
    gpu: &Gpu,
    mc: &MarchingCubesGpu,
    field: &FieldBuffer,
    key: (ChunkId, usize),
    queue: &RefCell<DeferredGeometry<(ChunkId, usize)>>,
) {
    let geometry = mc
        .extract_buffers(gpu.device(), gpu.queue(), field)
        .expect("extract a chunk");
    let _ = mc.take_timestamps(gpu.device(), gpu.queue());
    let bytes = u64::from(geometry.triangles) * 9 * 4;
    // A chunk with no surface has no bytes to read back, and a zero-size
    // `wgpu::Buffer` is invalid -- so the request list is empty rather than a
    // request for zero bytes. The slot is still consumed and still collected:
    // "this chunk produced nothing" is an answer a caller has to receive, and
    // the upper layer of M-124's world is 144 such chunks.
    let readback = if bytes == 0 {
        read_bytes_many_deferred(gpu.device(), gpu.queue(), &[]).expect("an empty read-back")
    } else {
        read_bytes_many_deferred(gpu.device(), gpu.queue(), &[(&geometry.positions, bytes)])
            .expect("a deferred read-back")
    };
    queue
        .borrow_mut()
        .submit(key, readback)
        .expect("`spend` checked has_room() before this chunk was reached");
}

/// One `(budget, delay)` cell: `BUDGET_FRAMES` frames of drain-then-schedule.
///
/// # Why the drain comes first, and why `spend` asks the queue
///
/// `mesh_within_budget` meshes **at least one chunk** whatever its predicate says
/// — the never-livelock guarantee, `dirty.rs` 105–112 — so a room check that
/// lived only in `spend` would let the first chunk of a frame submit into a full
/// queue and lose its geometry. The frame therefore installs what completed
/// before it schedules more, which is what a consumer does anyway, and `spend`
/// stops the pass on *either* limit: the clock or the queue.
///
/// The queue is sized by [`capacity_for`] so that the clock is the limit that
/// binds and `has_room` is the guard rather than the governor — `no_room_frames`
/// is the column that says whether that held. The alternative to a guard — spin
/// inside the frame until a slot frees — is a blocking path wearing the queue's
/// name, which is the defect `M-376` recorded in this arm's first version.
///
/// # No separate pump, unlike the ring arm
///
/// `M-376`'s ring consumed 1 read-back in 120 frames because nothing else
/// submitted and `PollType::Poll` is a single non-blocking check; the fix was to
/// take the frame's own GPU work as a closure. Here the scheduler *is* the
/// frame's GPU work and every chunk's `extract_buffers` contains a `poll(Wait)`,
/// which is a strictly stronger pump than that closure was. `collected` is the
/// control that says so.
fn sweep_budget(
    gpu: &Gpu,
    mc: &MarchingCubesGpu,
    world: &World,
    budget_us: u64,
    delay: usize,
    chunk_ms: f64,
) -> BudgetCell {
    let capacity = capacity_for(delay, budget_us as f64 / 1000.0, chunk_ms);
    let budget = Duration::from_micros(budget_us);
    let queue = RefCell::new(DeferredGeometry::new(capacity).expect("capacity is at least 1"));
    let mut dirty = world.all_dirty();

    let mut span_total_ms = 0.0f64;
    let mut passes = 0usize;
    let mut chunks = 0usize;
    let mut collected = 0usize;
    let mut max_in_flight = 0usize;
    let mut no_room_frames = 0usize;
    let mut empty_set_frames = 0usize;
    let mut latency_frames_total = 0usize;
    let mut latency_samples = 0usize;
    let mut max_latency_frames = 0usize;

    for frame in 0..BUDGET_FRAMES {
        for ((id, submitted), data) in queue.borrow_mut().drain_ready(gpu.device()).expect("drain")
        {
            std::hint::black_box(&data);
            collected += 1;
            let waited = frame - submitted;
            latency_frames_total += waited;
            latency_samples += 1;
            max_latency_frames = max_latency_frames.max(waited);
            // Re-dirty what came back, so the frames measure steady state rather
            // than one drain of a world that then sits idle. It is also
            // `game_budget`'s own behaviour: the queue re-fills the moment it
            // empties, because that is the state a player carving is in.
            dirty.insert(id);
        }
        if !queue.borrow().has_room() {
            no_room_frames += 1;
            continue;
        }
        if dirty.is_empty() {
            empty_set_frames += 1;
            continue;
        }

        let started = Instant::now();
        let report = dirty.mesh_within_budget(
            &world.layout,
            world.camera,
            |id, _| {
                extract_and_submit(gpu, mc, world.field(id), (id, frame), &queue);
            },
            || started.elapsed() < budget && queue.borrow().has_room(),
        );
        span_total_ms += started.elapsed().as_secs_f64() * 1000.0;
        passes += 1;
        chunks += report.meshed;
        max_in_flight = max_in_flight.max(queue.borrow().in_flight());
    }

    // The tail, bounded. A queue whose remainder cannot be collected in a bounded
    // number of polls is leaking, and that is a finding rather than something to
    // loop on forever.
    //
    // These collections are counted but their waits are **not** added to the
    // latency stats: a poll after the run is over is not a frame, and folding it
    // in would report a collision latency the running scheduler never paid.
    let mut drain_frames = 0usize;
    let mut tail = 0usize;
    while queue.borrow().in_flight() > 0 && drain_frames < 4_096 {
        let harvest = queue
            .borrow_mut()
            .drain_ready(gpu.device())
            .expect("drain the tail")
            .len();
        collected += harvest;
        tail += harvest;
        drain_frames += 1;
    }
    assert!(
        tail <= capacity,
        "the tail held {tail} read-backs at a capacity of {capacity}"
    );

    BudgetCell {
        budget_us,
        delay,
        capacity,
        passes,
        chunks,
        mean_ms: span_total_ms / passes as f64,
        mean_chunks: chunks as f64 / passes as f64,
        collected,
        max_in_flight,
        no_room_frames,
        latency_frames_total,
        latency_samples,
        max_latency_frames,
        empty_set_frames,
        drain_frames,
    }
}

/// The scheduler's own per-pass cost, with no GPU work in the closure at all.
///
/// **A control, and it is what makes C3's verdict attributable.** `spend` is
/// consulted after each chunk, so a pass costs *one chunk plus whatever
/// `mesh_within_budget` spends on itself* — and what it spends on itself is not
/// small at 288 chunks: it computes a squared distance per chunk and sorts the
/// whole set (`dirty.rs` 137–152), every call, before meshing anything.
///
/// Without this number a cell that misses C3's bar by 0.03 ms at a 25 µs budget
/// looks like the queue traded a stall for a queue. With it, the miss is
/// attributable to the scheduler's ordering pass, which is there with or without
/// a queue and is not what C3 is about.
///
/// Measured with `spend` returning `false` immediately, so exactly one chunk is
/// visited, and with a closure that does nothing but name it — and the id is
/// re-inserted so the set stays 288 and the sort stays the sort the sweep pays.
fn scheduler_overhead_ms(world: &World) -> f64 {
    let mut dirty = world.all_dirty();
    let mut total_ms = 0.0f64;
    for _ in 0..BUDGET_FRAMES {
        let mut visited = None;
        let started = Instant::now();
        let report = dirty.mesh_within_budget(
            &world.layout,
            world.camera,
            |id, _| visited = Some(id),
            || false,
        );
        total_ms += started.elapsed().as_secs_f64() * 1000.0;
        assert_eq!(
            report.meshed, 1,
            "`spend` is consulted after the first chunk"
        );
        dirty.insert(visited.expect("one chunk was visited"));
    }
    total_ms / BUDGET_FRAMES as f64
}

/// The mean cost of one chunk's extract-plus-submit, measured in this run.
///
/// C3's tolerance, so it cannot be a constant or a figure from another run
/// (`M-281`). Three things make it the *same* quantity the sweep spends:
///
/// - **The same population, in the same order.** All 288 chunks in the
///   scheduler's own nearest-first order, so the 144 empty chunks of `M-124`'s
///   upper layer are weighted here exactly as they are there.
/// - **No drain between chunks**, which is why the queue holds all 288. A drain
///   between them polls the device and lets each read-back retire before the next
///   submission, which is not what a pass does — and it measured 0.2754 ms
///   against a sweep that spent 0.32 ms on the same chunk, a 16% underestimate on
///   the number that is C3's own bar.
/// - **Only the extract-plus-submit is timed.** The collection is the next
///   frame's work and is outside every budget.
fn one_chunk_ms(gpu: &Gpu, mc: &MarchingCubesGpu, world: &World, order: &[ChunkId]) -> f64 {
    let queue = RefCell::new(DeferredGeometry::new(order.len()).expect("one slot per chunk"));
    let mut total_ms = 0.0f64;
    for &id in order {
        let started = Instant::now();
        extract_and_submit(gpu, mc, world.field(id), (id, 0), &queue);
        total_ms += started.elapsed().as_secs_f64() * 1000.0;
    }
    let mut spins = 0;
    while queue.borrow().in_flight() > 0 {
        let _ = queue
            .borrow_mut()
            .drain_ready(gpu.device())
            .expect("drain the calibration queue");
        spins += 1;
        assert!(spins < 100_000, "a calibration read-back never completed");
    }
    total_ms / order.len() as f64
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

        // ── C3: DeferredGeometry under mesh_within_budget, on M-124's fixture ──
        println!("\n-- budget: DeferredGeometry under DirtySet::mesh_within_budget --");
        let world = World::build(&gpu);
        let order = world.scheduler_order();
        let chunk_ms = one_chunk_ms(&gpu, &mc, &world, &order);
        let overhead_ms = scheduler_overhead_ms(&world);
        println!(
            "one chunk's extract-plus-submit, over all {CHUNKS} in scheduler order: \
             {chunk_ms:.4} ms  (C3's tolerance)\n\
             mesh_within_budget's own per-pass cost over {CHUNKS} chunks, no GPU work: \
             {overhead_ms:.4} ms  (control)\n\
             so a one-chunk pass is expected at {:.4} ms",
            chunk_ms + overhead_ms
        );
        println!(
            "{:>8} {:>3} {:>4} {:>7} {:>10} {:>9} {:>8} {:>9} {:>7} {:>7} {:>8} {:>6}",
            "budgetUs",
            "dly",
            "cap",
            "passes",
            "meanPassMs",
            "chunks/p",
            "msPerChk",
            "overshoot",
            "maxFlt",
            "noRoom",
            "latency",
            "C3"
        );
        let mut cells: Vec<BudgetCell> = Vec::with_capacity(BUDGETS_US.len() * DELAYS.len());
        for &delay in &DELAYS {
            for &budget_us in &BUDGETS_US {
                let cell = sweep_budget(&gpu, &mc, &world, budget_us, delay, chunk_ms);

                // ── the controls, asserted rather than printed ────────────────
                //
                // 1. The queue must actually fill. A queue that never holds
                //    anything is measuring a blocking path wearing the queue's
                //    name -- M-376's first ring arm, which consumed 1 read-back
                //    in 120 frames and reported 0.0004 ms.
                assert!(
                    cell.max_in_flight > 0,
                    "budget {budget_us} us at delay {delay}: the queue never held a \
                     read-back, so nothing here is a deferred cost"
                );
                // 2. Geometry submitted must be geometry collected, to within one
                //    queue-full. Anything less means read-backs are being paid for
                //    and dropped.
                assert!(
                    cell.collected + cell.capacity >= cell.chunks,
                    "budget {budget_us} us at delay {delay}: {} chunks submitted, \
                     {} collected -- {} read-backs were never harvested",
                    cell.chunks,
                    cell.collected,
                    cell.chunks.saturating_sub(cell.collected + cell.capacity)
                );
                // 3. Every pass must have had a set to work on.
                assert_eq!(
                    cell.empty_set_frames, 0,
                    "budget {budget_us} us at delay {delay}: the dirty set emptied, \
                     so some frames measured an idle scheduler"
                );
                assert!(
                    cell.passes > 0,
                    "budget {budget_us} us at delay {delay}: no frame ran a budgeted \
                     pass"
                );
                assert!(
                    cell.drain_frames < 4_096,
                    "budget {budget_us} us at delay {delay}: the tail did not drain \
                     in a bounded number of polls, so the queue leaks"
                );
                // 4. The clock must be the limit that binds, not the queue. If the
                //    queue is refusing frames then `capacity_for` under-sized it
                //    and the cell measures `budget / chunk_ms` against a slot
                //    count rather than the scheduler -- the exact defect that
                //    made a flat capacity useless here.
                assert_eq!(
                    cell.no_room_frames, 0,
                    "budget {budget_us} us at delay {delay}: the queue refused {} of \
                     {BUDGET_FRAMES} frames at capacity {}, so this cell measures \
                     the queue's depth and not the budget",
                    cell.no_room_frames, cell.capacity
                );

                println!(
                    "{:>8} {:>3} {:>4} {:>7} {:>10.4} {:>9.3} {:>8.4} {:>9.3} {:>7} {:>7} {:>8} {:>6}",
                    cell.budget_us,
                    cell.delay,
                    cell.capacity,
                    cell.passes,
                    cell.mean_ms,
                    cell.mean_chunks,
                    cell.ms_per_chunk(),
                    cell.overshoot_chunks(chunk_ms),
                    cell.max_in_flight,
                    cell.no_room_frames,
                    format!(
                        "{:.2}/{}",
                        cell.mean_latency_frames(),
                        cell.max_latency_frames
                    ),
                    if cell.within_one_chunk(chunk_ms) {
                        "ok"
                    } else {
                        "OUT"
                    }
                );
                rows.push(vec![
                    ("arm", "budget".to_string()),
                    ("entry_point", "DeferredGeometry".to_string()),
                    (
                        "cells",
                        (u64::from(CHUNK_CELLS).pow(3) * CHUNKS as u64).to_string(),
                    ),
                    ("amortised_ms_per_frame", format!("{:.6}", cell.mean_ms)),
                    ("budget_chunks", format!("{:.6}", cell.mean_chunks)),
                    (
                        "within_one_chunk",
                        cell.within_one_chunk(chunk_ms).to_string(),
                    ),
                    ("ring_frames_delay", cell.delay.to_string()),
                    ("queue_capacity", cell.capacity.to_string()),
                    ("budget_us", cell.budget_us.to_string()),
                    ("budget_ms", format!("{:.6}", cell.budget_ms())),
                    ("one_chunk_ms", format!("{chunk_ms:.6}")),
                    ("scheduler_overhead_ms", format!("{overhead_ms:.6}")),
                    ("ms_per_chunk", format!("{:.6}", cell.ms_per_chunk())),
                    (
                        "overshoot_chunks",
                        format!("{:.6}", cell.overshoot_chunks(chunk_ms)),
                    ),
                    ("budget_frames", BUDGET_FRAMES.to_string()),
                    ("budget_passes", cell.passes.to_string()),
                    ("chunks_meshed", cell.chunks.to_string()),
                    ("readbacks_collected", cell.collected.to_string()),
                    ("max_in_flight", cell.max_in_flight.to_string()),
                    ("no_room_frames", cell.no_room_frames.to_string()),
                    (
                        "mean_latency_frames",
                        format!("{:.6}", cell.mean_latency_frames()),
                    ),
                    ("max_latency_frames", cell.max_latency_frames.to_string()),
                    ("drain_frames", cell.drain_frames.to_string()),
                    (
                        "wall_ms",
                        format!("{:.6}", cell.mean_ms * cell.passes as f64),
                    ),
                ]);
                cells.push(cell);
            }
        }

        // **The fifth control: chunks per pass must RISE with the budget.** A
        // scheduler ignoring its budget would be flat, and flat also satisfies the
        // tolerance at the low end, so the tolerance alone cannot see it. Asserted
        // at every delay, because `capacity_for` keeps the queue out of the way at
        // all three -- if a delay could not show the rise, control 4 above would
        // have fired first.
        for &delay in &DELAYS {
            let arm: Vec<&BudgetCell> = cells.iter().filter(|c| c.delay == delay).collect();
            let bottom = arm.first().expect("eight budgets per delay");
            let top = arm.last().expect("eight budgets per delay");
            assert!(
                top.mean_chunks > bottom.mean_chunks + 0.5,
                "delay {delay}: {:.3} chunks a pass at {} us and {:.3} at {} us -- a \
                 320x budget range moved the count by less than half a chunk, so the \
                 pass is not tracking its budget",
                bottom.mean_chunks,
                bottom.budget_us,
                top.mean_chunks,
                top.budget_us
            );
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
        // **C3 is measured now, and it is a conjunction over every row.** The
        // clause is "holds the amortised per-frame cost within one chunk of the
        // budget across a 320x range", so one cell outside the tolerance falsifies
        // it -- there is no averaging over the sweep and no scoring of a chosen
        // subset. `ring_frames_delay` is a swept column, so the conjunction runs
        // over all eight budgets at each of the three delays.
        let c3_rows = cells.len();
        let out: Vec<&BudgetCell> = cells
            .iter()
            .filter(|c| !c.within_one_chunk(chunk_ms))
            .collect();
        let c3 = if out.is_empty() { "true" } else { "false" };

        println!(
            "\nC1 largest component at 129³: {} -> {}",
            at129.largest(),
            if c1 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C2 synchronisation removed at 129³: {removed_share:.4} -> {}",
            if c2 { "HELD" } else { "FALSIFIED" }
        );
        let lo = cells
            .iter()
            .map(|c| c.overshoot_chunks(chunk_ms))
            .fold(f64::INFINITY, f64::min);
        let hi = cells
            .iter()
            .map(|c| c.overshoot_chunks(chunk_ms))
            .fold(0.0f64, f64::max);
        println!(
            "C3 {} : {} of {c3_rows} (budget, delay) cells hold |meanPassMs - budget| \
             <= one_chunk_ms = {chunk_ms:.4} ms, over a 320x budget range at delays \
             {DELAYS:?} frames.\n\
             The overshoot spans {lo:.3} to {hi:.3} chunks -- the bar is 1.000, which \
             is where `mesh_within_budget`'s own doc puts the never-livelock price, \
             so the cells straddle it and the COUNT above moves between runs while \
             this band does not.",
            if out.is_empty() { "HELD" } else { "FALSIFIED" },
            c3_rows - out.len()
        );
        for cell in &out {
            println!(
                "     out: {:>5} us at delay {} -- {:.4} ms a pass, {:.3} chunks at \
                 {:.4} ms each, overshoot {:.3} chunks",
                cell.budget_us,
                cell.delay,
                cell.mean_ms,
                cell.mean_chunks,
                cell.ms_per_chunk(),
                cell.overshoot_chunks(chunk_ms)
            );
        }
        // Where the overshoot comes from, attributed rather than guessed. The
        // scheduler's own ordering pass was the obvious suspect and the control
        // refutes it: it is measured here, not assumed.
        let single: Vec<&BudgetCell> = cells.iter().filter(|c| c.mean_chunks < 1.05).collect();
        let single_ms = single.iter().map(|c| c.ms_per_chunk()).sum::<f64>() / single.len() as f64;
        let many: Vec<&BudgetCell> = cells.iter().filter(|c| c.mean_chunks > 10.0).collect();
        let many_ms = many.iter().map(|c| c.ms_per_chunk()).sum::<f64>() / many.len() as f64;
        println!(
            "attribution: the scheduler's ordering pass is {overhead_ms:.4} ms, \
             {:.1}% of the {chunk_ms:.4} ms bar, so it is not the overshoot. A pass \
             that meshes ONE chunk spends {single_ms:.4} ms on it; a pass that meshes \
             more than ten spends {many_ms:.4} ms each, a {:.2}x premium on a pass's \
             first chunk -- which is M-159's mechanism at pass granularity, because \
             `extract_buffers` waits with no submission index and drains the previous \
             pass's outstanding deferred copies.",
            100.0 * overhead_ms / chunk_ms,
            single_ms / many_ms
        );
        // The registration's owner-question, answered as a number rather than
        // surfaced as a sentence: how many frames a chunk's geometry actually
        // waited, per delay, worst case included.
        for &delay in &DELAYS {
            let arm: Vec<&BudgetCell> = cells.iter().filter(|c| c.delay == delay).collect();
            let worst = arm.iter().map(|c| c.max_latency_frames).max().unwrap_or(0);
            let mean = arm.iter().map(|c| c.mean_latency_frames()).sum::<f64>() / arm.len() as f64;
            println!(
                "collision latency at delay {delay}: mean {mean:.2} frames, worst \
                 {worst} over the eight budgets"
            );
        }
        println!(
            "\nTHE OWNER'S QUESTION, surfaced and not answered: the queue costs \n\
             {DEPTH_NOTE}"
        );

        // `budget_chunks`, `within_one_chunk`, `amortised_ms_per_frame` and
        // `ring_frames_delay` are **per-row** now, not aggregates: `record`
        // collects into a BTreeMap and the extend below happens after each row is
        // built, so an aggregate of the same name would silently overwrite the
        // measurement it is meant to summarise.
        let aggregates: Row = vec![
            ("c1_holds", c1.to_string()),
            ("c2_holds", c2.to_string()),
            ("c3_holds", c3.to_string()),
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
