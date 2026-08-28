//! **P-95 — undo, and the checkpoint cadence nobody has measured.**
//!
//! Ticket: R-095. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p95
//! ```
//!
//! Writes `docs/experiments/p-95.csv`.
//!
//! # The two things undo can be, and why that is the whole experiment
//!
//! Edits compose rather than mutate: the field is a [`BrushStack`] over a base
//! field and carving pushes a brush, so **every field sample walks every
//! brush**. That gives exactly two ways to reach the state one edit back.
//!
//! - **The log path.** Pop the brush and re-mesh the dirty chunks against
//!   `base ⊕ log[..n-1]`. Costs one base sample plus `n-1` brush samples per
//!   grid sample. No memory beyond the log.
//! - **The checkpoint path.** Re-mesh against `checkpoint(c) ⊕ log[c..n-1]`,
//!   where the checkpoint is a dense snapshot of the fold's own value at every
//!   sample of the world lattice. Costs one trilinear lookup plus `n-1-c` brush
//!   samples. Costs `(cells+1)³` scalars of memory per checkpoint.
//!
//! The checkpoint is **exact at a lattice point**: it stores the `f64` the fold
//! produced there, and `brush::apply` is a pure function of that value and the
//! brush's own sample, so folding the remaining brushes onto it reproduces the
//! full fold bit for bit. `CELL_SIZE` is a power of two and `ORIGIN` a multiple
//! of it, so `(p - origin)/h` is exact on the lattice and the interpolation
//! weights are exactly 0 or 1 there — M-32's spacing rule, used here for the
//! same reason `game_dig` uses it.
//!
//! # SHARE, recomputed before the harness was written — and C1 is unreachable
//!
//! **C1's registered direction cannot hold, for any instrument, by
//! monotonicity.** Re-folding from a checkpoint costs one lookup plus `d` brush
//! samples where `d` is the number of edits separating the checkpoint from the
//! state being reached, so `refold(d)` is monotone non-decreasing in `d`.
//! Undoing the last edit through the log does not depend on `d` at all. So
//! `{ d : undo < refold(d) }` is an **up-set**: if the inequality ever holds it
//! holds for every larger `d`. A threshold of the registered form — "undo costs
//! strictly less *whenever fewer than* `N` edits separate them" — describes a
//! down-set and therefore cannot exist. This is `✗51`'s rule applied to a
//! direction rather than a magnitude, and it is stated here before the run.
//!
//! Worse for C1, the up-set is empty inside the log. At head `n`, undo walks
//! `base + (n-1)` brushes and `refold(d)` walks `lookup + d` brushes, so
//! `undo < refold(d)` needs `d > (n-1) + (base - lookup)/brush`, which exceeds
//! `n-1` for every field whose base sample costs more than an array index. The
//! harness sweeps `d = 0 … n-1` anyway and reports the smallest `d` that
//! crosses, or `0` for none — `P-70`'s precedent: say it before the run, then
//! produce the number.
//!
//! **And then `crossover_n` proved not to be reproducible, which is why there is
//! a second instrument beside it.** The offset `(base - lookup)/brush` is a
//! constant independent of `n`, so a crossover — if one exists — always sits
//! within a fixed number of edits of the *bottom* of the log, i.e. at the
//! boundary of any search that can be run. Two runs of this harness put
//! `sphere_r1_2`'s crossover at `d = 125` of 127 and at none at all. That is
//! `P-72`'s defect exactly: a boundary hit reported as a location. So the sweep's
//! **128 measurements are fitted** as `refold(d) = A + B·d`, and the answer is
//! reported as `crossover_n_fitted = (undo − A)/B` together with
//! `crossover_margin_edits_past_log = crossover_n_fitted − (n-1)` — how many
//! edits *past the bottom of the log* a checkpoint would have to sit before
//! undoing became cheaper. A positive margin says no crossover can exist at any
//! log length, since `d ≤ n-1` always: a checkpoint cannot be below the base
//! field. `B` is also the marginal cost of one edit, measured over 128 points
//! rather than four buckets, and the margin is a ratio of times — which is what
//! `M-280` asks for on a governed clock.
//!
//! **C2's ratio clause is reachable, and the arithmetic is here.** Least
//! squares through M-50's four bucket medians at their bucket midpoints
//! (`8, 23, 38, 53` edits → `0.158, 0.354, 0.525, 0.589` ms per re-meshed
//! chunk) gives `cost(L) = 0.10882 + 0.0097600·L` ms. With `m` dirty chunks the
//! worst single undo under cadence `k` costs `m·cost(k)`, so
//! `k ≤ (16/m − 0.10882)/0.0097600`: **1628 at `m = 1`, 399 at `m = 4`, 194 at
//! `m = 8`, 50 at `m = 27`**. The fixed per-chunk term alone only eats the frame
//! at `m ≥ 147` chunks, which is 5.4× more chunks than this world contains. So
//! a `k` satisfying the budget exists with room, and the clause is an upper
//! bound on `k` rather than an unreachable floor.
//!
//! **C2's verdict criterion, declared before the run.** `predicted_k_from_m50`
//! and `measured_k` are computed by the *same* formula from the *same* 16 ms
//! budget and the *same* measured worst-case dirty-chunk count — one from
//! M-50's four bucket medians, one from this harness's four bucket medians over
//! the same log-length buckets. C2 holds iff they agree within 2× **and** the
//! trace at the largest swept cadence not exceeding `measured_k` really keeps
//! `worst_undo_ms` under 16. 2× rather than tighter because M-50 is an `f32`
//! layout under a mouse on a different base field and this is `f64` — a bar set
//! after seeing the numbers would be worth nothing, so it is set here. The four
//! bucket ratios are reported so a reader can see whether any disagreement is a
//! scale factor (same function, different constant) or a shape difference (a
//! different function, which is C2's registered falsifier).
//!
//! # C3, and the one place a hash can actually move
//!
//! Through the log, undo-then-redo returning the same hash is an **identity**,
//! not a measurement: the field at length `L` is a pure function of `log[..L]`,
//! so no trace can refute it. That is `M-44`'s vacuous zero, and reporting it
//! would be a HELD with no instrument.
//!
//! So the trace's state is always represented the way an editor with
//! checkpoints represents it — `checkpoint(c) ⊕ log[c..L]` for the deepest
//! valid checkpoint `c` — and the hash recorded when an edit is pushed must
//! come back when it is undone and redone **across a ladder that was torn down
//! and rebuilt in between**. An undo below `c` drops that checkpoint; the redo
//! rebuilds it from the one beneath. If the rebuild were not bit-exact, or if
//! the ladder came back at a different `c`, the hash moves. Both are asserted:
//! `checkpoint_rebuild_mismatches` compares the rebuilt snapshot against the
//! dropped one scalar for scalar.
//!
//! And the instrument is shown to have teeth by the column beside it:
//! `cp_vs_log_hash_mismatches` runs the *same* hash over the *same* chunks with
//! the checkpoint path against the full fold from the base field. Positions
//! agree bit-exactly — they are computed from lattice corner values, which the
//! checkpoint stores exactly — but normals are
//! [`Sdf::gradient`](isomesh::Sdf::gradient), a central difference at
//! `p ± DIFF_STEP·max(|p|,1)`, and that lands **off** the lattice, where a
//! trilinear snapshot is not the field it was baked from. `max_normal_dev_deg`
//! is how far.
//!
//! # VACUITY CONTROL
//!
//! Registered: *the interleaved-undo trace must include undos that cross a
//! checkpoint boundary, or C3 tests only the cheap path.* Reported as
//! `undos_crossing_checkpoint` and asserted non-zero.
//!
//! **It could have been zero, and the harness reports the count that nearly
//! was** rather than only asserting the total. The regular burst schedule is
//! cadence-independent — a burst every `BURST_EVERY` edits, lengths cycling
//! through `BURST_LENS` — so whether a regular burst ever reaches a checkpoint
//! is an accident of arithmetic: a crossing needs a multiple of `k` inside
//! `[l-u+1, l]` for some burst at length `l` of depth `u`. At `k = 128`,
//! `BURST_EVERY = 25` and `u ≤ 9` exactly **one** such coincidence exists in a
//! thousand edits — the burst at `l = 900`, depth 7, over the checkpoint at 896.
//! `crossings_regular_only` is that count, and it falls 47 → 22 → 12 → 6 → 2 → 1
//! across the cadence sweep: `P-62`'s "a hair from `M-44`'s vacuous zero". The
//! first draft of this comment claimed it *was* zero at `k = 128`, and the run
//! refuted that. So the total is topped up by four *boundary bursts* placed
//! deliberately at `c + 1` for four checkpoints spread through the trace, each
//! three deep, so the second undo of each one pops the ladder at every cadence
//! by construction rather than by luck.
//!
//! # Instruments, stated plainly
//!
//! Time is `Instant` around one [`Extractor::extract_into`] per chunk, over the
//! set of chunks `mark_edit` marked for that one edit — the work
//! `DirtySet::mesh_dirty` would hand out, timed directly so the closure and the
//! set-clearing are not in the clock. That is M-50's quantity, at M-50's chunk
//! size (`16` cells) and cell size (`0.125`), so the per-chunk milliseconds are
//! comparable to its four buckets.
//! `game_dig` is an `f32` layout; this is `f64`, because
//! [`mesh_hash`](isomesh::validate::mesh_hash) is `f64`-only and C3 needs it.
//! That is a deviation, it is why C2's bar is 2× rather than tighter, and every
//! comparison inside a row is within one build and one run (`M-281`). The CPU
//! clock is on every row (`M-280`).
//!
//! # Controls, every one of them a column and an assertion
//!
//! - Every timed edit must dirty at least one chunk. `empty_dirty_edits` is
//!   reported and asserted zero: a trace that marked nothing would report a
//!   fast time for doing nothing, which is `M-44`.
//! - The three ops must all appear in the log. A log of nothing but
//!   `Subtract` reorders bit-exactly (`max` is associative and exact), so C3
//!   over one would be testing the easiest possible history. The trace cycles
//!   `Subtract, Subtract, Add, SmoothAdd`, which is semantically
//!   order-dependent (`M-37`) and numerically order-dependent (`M-38`).
//! - A checkpoint rebuilt after being dropped must be bit-identical to the one
//!   that was dropped.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::print_literal,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::brush::{Brush, BrushOp, BrushStack, apply};
use isomesh::chunk::dirty::{DirtySet, mark_edit};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::extractor::Extractor;
use isomesh::fields::{FbmTerrain, Gyroid, Sphere};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Cells per axis in the world. `48` so it divides `CHUNK_CELLS` exactly.
const WORLD_CELLS: u32 = 48;

/// Cells per chunk edge. **`game_dig`'s**, so per-chunk milliseconds here mean
/// the same thing as M-50's four buckets.
const CHUNK_CELLS: u32 = 16;

/// Cell size. `game_dig`'s `0.125`: a power of two, where M-32 measured that two
/// chunks agree on their shared sample plane bit for bit, and where the
/// checkpoint's lattice lookup is exact.
const CELL_SIZE: f64 = 0.125;

/// World origin, centred so all three fields have a surface inside.
const ORIGIN: f64 = -0.5 * CELL_SIZE * WORLD_CELLS as f64;

/// Samples per axis in a checkpoint: the world's whole cell lattice, corners
/// included, so every chunk's own sample grid is a subset of it and the seam
/// samples are shared rather than duplicated.
const LATTICE: usize = WORLD_CELLS as usize + 1;

/// Head length for the crossover sweep, and the length the cost curve is
/// measured over. `128` so that `d = k - 1` is inside the sweep for every
/// cadence, and `>= 60` so M-50's four buckets are all covered.
const HEAD: usize = 128;

/// Edits in the C3 trace. The registration's `10³`.
const TRACE_EDITS: usize = 1000;

/// Checkpoint cadences swept.
const CADENCES: [usize; 6] = [4, 8, 16, 32, 64, 128];

/// A regular undo burst every this many pushes. Cadence-independent, so the
/// rows are comparable.
const BURST_EVERY: usize = 25;

/// Regular burst depths, cycled. Coprime-ish with the cadences on purpose: at a
/// coarse cadence none of them reaches a checkpoint, which is what makes
/// `crossings_regular_only` a control rather than a decoration.
const BURST_LENS: [usize; 9] = [1, 2, 3, 5, 8, 9, 4, 6, 7];

/// Depth of each boundary burst — the four bursts placed at `c + 1` so that
/// their second undo pops the ladder at every cadence.
const BOUNDARY_BURST: usize = 3;

/// Every this many undone edits, also compare the checkpoint path against the
/// full fold from the base field. The expensive check, so it is sampled.
const CP_VS_LOG_EVERY: usize = 20;

/// Every this many undone edits, hash **every** chunk in the world rather than
/// the edit's own dirty set — a subtracting brush moves values outside its own
/// bounding box, so a dirty-set hash could in principle miss a change.
const WORLD_HASH_EVERY: usize = 60;

/// Repetitions of each timed measurement; median taken. `M-337`'s lesson: a
/// single reading on a governed CPU is the mistake that turned a registered
/// `1.25×` into `1.022×` three runs later.
const REPS: usize = 3;

/// M-50's four bucket medians, at their bucket midpoints.
///
/// `1–15 / 16–30 / 31–45 / 46–60` edits in the log against median milliseconds
/// per re-meshed chunk. Read from the registration, which quotes them, and from
/// `bevy_isomesh/examples/game_dig.rs`'s own table, which agrees exactly.
/// **There is no committed CSV for M-50** — see the report.
const M50: [(f64, f64); 4] = [
    (8.0, 0.158),
    (23.0, 0.354),
    (38.0, 0.525),
    (53.0, 0.589),
];

/// The frame C2 is denominated in.
const FRAME_MS: f64 = 16.0;

/// The edit log: one shape and one operation per entry, applied first to last.
type Log = Vec<Brush<Sphere<f64>>>;

/// A dense snapshot of the fold's own value at every sample of the world
/// lattice — what "re-fold from a checkpoint" starts from.
///
/// Exact at a lattice point by construction, and a trilinear interpolant
/// everywhere else. Outside the lattice it is extended by its boundary value;
/// that is the definition of the field a dense snapshot *is*, not a fallback
/// for one, and the trace keeps every brush well inside the lattice so the
/// extension is never load-bearing.
struct Checkpoint {
    values: Vec<f64>,
}

impl Checkpoint {
    fn bytes(&self) -> u64 {
        (self.values.len() * core::mem::size_of::<f64>()) as u64
    }
}

impl Sdf for Checkpoint {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        // `(p - ORIGIN)/CELL_SIZE` is exact on the lattice: the spacing is a
        // power of two and the origin a multiple of it. So `floor` lands on the
        // integer and the weight is exactly zero, or — on the far face — the
        // cell below with a weight of exactly one.
        let axis = |a: usize| -> (usize, f64) {
            let u = (p[a] - ORIGIN) / CELL_SIZE;
            let top = (LATTICE - 2) as f64;
            let base = u.floor().clamp(0.0, top);
            ((base as usize), (u - base).clamp(0.0, 1.0))
        };
        let (i, tx) = axis(0);
        let (j, ty) = axis(1);
        let (k, tz) = axis(2);

        let at = |di: usize, dj: usize, dk: usize| -> f64 {
            self.values[((k + dk) * LATTICE + (j + dj)) * LATTICE + (i + di)]
        };
        // `(1-t)·a + t·b`, not `a + t·(b-a)`: the first is exact at both ends
        // (`1·a + 0·b` and `0·a + 1·b`), the second only at `t = 0`.
        let mix = |a: f64, b: f64, t: f64| (1.0 - t) * a + t * b;

        let x00 = mix(at(0, 0, 0), at(1, 0, 0), tx);
        let x10 = mix(at(0, 1, 0), at(1, 1, 0), tx);
        let x01 = mix(at(0, 0, 1), at(1, 0, 1), tx);
        let x11 = mix(at(0, 1, 1), at(1, 1, 1), tx);
        mix(mix(x00, x10, ty), mix(x01, x11, ty), tz)
    }
}

/// The world position of a lattice index.
fn lattice_point(i: usize, j: usize, k: usize) -> [f64; 3] {
    [
        ORIGIN + CELL_SIZE * i as f64,
        ORIGIN + CELL_SIZE * j as f64,
        ORIGIN + CELL_SIZE * k as f64,
    ]
}

/// Bake the base field with no brushes applied — the checkpoint at `c = 0`.
fn bake<F: Sdf<Scalar = f64>>(base: &F) -> Checkpoint {
    let mut values = Vec::with_capacity(LATTICE * LATTICE * LATTICE);
    for k in 0..LATTICE {
        for j in 0..LATTICE {
            for i in 0..LATTICE {
                values.push(base.sample(lattice_point(i, j, k)));
            }
        }
    }
    Checkpoint { values }
}

/// Fold `brushes` onto an existing checkpoint, giving the checkpoint that many
/// edits later.
///
/// Bit-exact against a full fold from the base field, by induction: the stored
/// value at a lattice point *is* the fold's value there, and
/// [`apply`] is a pure function of it.
fn advance(from: &Checkpoint, brushes: &[Brush<Sphere<f64>>]) -> Checkpoint {
    let mut values = from.values.clone();
    for k in 0..LATTICE {
        for j in 0..LATTICE {
            for i in 0..LATTICE {
                let p = lattice_point(i, j, k);
                let idx = (k * LATTICE + j) * LATTICE + i;
                let mut v = values[idx];
                for b in brushes {
                    v = apply(b.op, v, b.shape.sample(p));
                }
                values[idx] = v;
            }
        }
    }
    Checkpoint { values }
}

/// The chunk grid, its sample shape, and every chunk in it.
struct World {
    layout: ChunkLayout<f64>,
    shape: RuntimeShape3,
    all: Vec<ChunkId>,
}

impl World {
    fn new() -> Self {
        let layout =
            ChunkLayout::<f64>::new(CHUNK_CELLS, CELL_SIZE, [ORIGIN; 3]).expect("chunk layout");
        let shape = layout.sample_shape().expect("sample shape");
        let per_axis = (WORLD_CELLS / CHUNK_CELLS) as i32;
        let mut all = Vec::new();
        for cz in 0..per_axis {
            for cy in 0..per_axis {
                for cx in 0..per_axis {
                    all.push(ChunkId { coords: [cx, cy, cz] });
                }
            }
        }
        Self { layout, shape, all }
    }
}

/// Extract one chunk.
fn extract<F: Sdf<Scalar = f64>>(
    world: &World,
    mc: &mut MarchingCubes<f64>,
    field: &F,
    id: ChunkId,
) -> MeshBuffer<f64> {
    let mut out = MeshBuffer::<f64>::new();
    mc.extract_into(
        field,
        &world.shape,
        world.layout.sample_origin(id),
        CELL_SIZE,
        &mut out,
    )
    .expect("extract");
    out
}

/// FNV-1a over a list of per-chunk mesh hashes. Order-dependent, so the chunk
/// list is sorted by its caller.
fn combine(hashes: &[u64]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for x in hashes {
        for byte in x.to_le_bytes() {
            h ^= u64::from(byte);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    h
}

/// The golden hash of a set of chunks under one field, and how many vertices
/// went into it.
///
/// The vertex count is the `M-44` control on C3: a hash taken over a set of
/// empty meshes is a constant, so it could not have moved, so a zero mismatch
/// count over it would mean nothing.
fn hash_chunks<F: Sdf<Scalar = f64>>(
    world: &World,
    mc: &mut MarchingCubes<f64>,
    field: &F,
    chunks: &[ChunkId],
) -> (u64, usize) {
    let mut per = Vec::with_capacity(chunks.len());
    let mut vertices = 0usize;
    for id in chunks {
        let mesh = extract(world, mc, field, *id);
        vertices += mesh.positions.len();
        per.push(mesh_hash(&mesh));
    }
    (combine(&per), vertices)
}

/// Median wall time, over [`REPS`], to re-mesh a set of chunks under one field.
///
/// This is M-50's quantity before it is divided by the chunk count: the time
/// `mesh_dirty` spends on one edit's dirty set. Every caller here uses this one
/// function, so no two timings in a row come from differently shaped loops.
fn time_extract<F: Sdf<Scalar = f64>>(
    world: &World,
    mc: &mut MarchingCubes<f64>,
    field: &F,
    chunks: &[ChunkId],
) -> f64 {
    let mut times = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let start = Instant::now();
        for id in chunks {
            std::hint::black_box(extract(world, mc, field, *id));
        }
        times.push(start.elapsed().as_nanos() as f64 / 1e6);
    }
    median(times)
}

/// The cell box one brush can touch, one cell of margin either side.
fn brush_box(world: &World, brush: &Brush<Sphere<f64>>) -> ([i64; 3], [i64; 3]) {
    let c = brush.shape.center;
    let r = brush.shape.radius;
    let lo = world.layout.cell_of([c[0] - r, c[1] - r, c[2] - r]).map(|v| v - 1);
    let hi = world.layout.cell_of([c[0] + r, c[1] + r, c[2] + r]).map(|v| v + 1);
    (lo, hi)
}

/// The chunks whose *triangles* change when brush `j` is applied on top of the
/// prefix already folded into `base`.
///
/// Through [`mark_edit`], which is the crate's own instrument for this and the
/// one M-50's E1 is defined over. Representation-independent: `mark_edit`
/// samples the global sample lattice, where a checkpoint is exact.
fn dirty_chunks<B: Sdf<Scalar = f64>>(
    world: &World,
    base: &B,
    brushes_before: &[Brush<Sphere<f64>>],
    brushes_after: &[Brush<Sphere<f64>>],
    region: ([i64; 3], [i64; 3]),
) -> Vec<ChunkId> {
    let before = BrushStack { base, brushes: brushes_before };
    let after = BrushStack { base, brushes: brushes_after };
    let mut dirty = DirtySet::new();
    mark_edit(&world.layout, &before, &after, region.0, region.1, &mut dirty).expect("mark_edit");
    let mut out: Vec<ChunkId> = dirty.iter().collect();
    out.sort_unstable_by_key(|id| id.coords);
    out
}

/// A deterministic 64-bit LCG, so the trace is the same on every run and every
/// machine. Numerical Recipes' constants.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }
}

/// Build the edit log: a cloud of overlapping spheres, every one of which
/// actually moves the surface.
///
/// Centres inside `±1.5` so every brush and every gradient probe stays well
/// within the lattice. Radii `0.30 … 0.70` world units — `2.4` to `5.6` cells,
/// so a brush straddles a chunk boundary often and dirties several chunks.
///
/// The op cycle is the point: all-`Subtract` reorders bit-exactly, so a log of
/// one would make C3's exactness claim trivial. `SmoothAdd` is not even
/// associative (`M-38`), which is the hardest history to be faithful to.
///
/// # The first run's own control rejected the first fixture, and this is the fix
///
/// With centres drawn blind, **8 of the first 128 edits dirtied no chunk at
/// all** and the `M-44` assertion refused to report a time for them. The
/// mechanism is `min`: an [`Add`](BrushOp::Add) brush sitting in bulk solid
/// computes `min(field, shape)` where `field` is already the more negative of
/// the two, so it changes nothing — not a slow edit, a *no* edit. And an undo of
/// a no-op edit is worse than useless for C3: its hashed region would be a set
/// of empty meshes, whose hash is a constant, so the check could not have
/// failed. That is `M-44` in C3's clothing.
///
/// So a candidate brush is accepted only if it **flips the sign of at least one
/// lattice sample inside its own bounding box**, tested against a running exact
/// snapshot of the fold. A sign flip somewhere guarantees a cell with mixed
/// signs on one side of the edit, hence an output-changed cell, hence a dirty
/// chunk. Rejected candidates are resampled; running out of tries is a panic
/// rather than a shrug, because a fixture that cannot place an edit is not a
/// fixture.
///
/// The log is therefore **per field**, which `P-72`'s per-field surface probe is
/// the precedent for and for the same reason. The radius and operation streams
/// come from separate generators, so those are identical on every field and only
/// the centres move.
fn build_log<F: Sdf<Scalar = f64>>(base: &F) -> Log {
    /// Tries before giving up on a position for one edit.
    const TRIES: usize = 256;

    let mut centres = Lcg(0x5EED_0095_0000_0001);
    let mut radii = Lcg(0x5EED_0095_C0FF_EE01);
    let mut grid = bake(base).values;
    let mut log = Log::with_capacity(TRACE_EDITS);

    for i in 0..TRACE_EDITS {
        let radius = radii.range(0.30, 0.70);
        let mut placed = None;
        for _ in 0..TRIES {
            let center = [
                centres.range(-1.5, 1.5),
                centres.range(-1.5, 1.5),
                centres.range(-1.5, 1.5),
            ];
            let shape = Sphere { center, radius };
            let brush = match i % 4 {
                0 | 1 => Brush::subtract(shape),
                2 => Brush::add(shape),
                _ => Brush::smooth_add(shape, 0.1),
            };
            let span = radius + CELL_SIZE;
            let bound = |v: f64| -> usize {
                (((v - ORIGIN) / CELL_SIZE).floor()).clamp(0.0, (LATTICE - 1) as f64) as usize
            };
            let lo = [0, 1, 2].map(|a| bound(center[a] - span));
            let hi = [0, 1, 2].map(|a| bound(center[a] + span) + 1);
            let mut flips = false;
            'scan: for k in lo[2]..hi[2].min(LATTICE) {
                for j in lo[1]..hi[1].min(LATTICE) {
                    for ii in lo[0]..hi[0].min(LATTICE) {
                        let p = lattice_point(ii, j, k);
                        let v = grid[(k * LATTICE + j) * LATTICE + ii];
                        let w = apply(brush.op, v, shape.sample(p));
                        if (v < 0.0) != (w < 0.0) {
                            flips = true;
                            break 'scan;
                        }
                    }
                }
            }
            if flips {
                placed = Some(brush);
                break;
            }
        }
        let brush = placed.unwrap_or_else(|| {
            panic!(
                "edit {i}: {TRIES} candidate positions in a row moved no sign, so the fixture \
                 has run out of surface to edit"
            )
        });
        // Fold it into the running snapshot over the whole lattice, not just the
        // bounding box: a `Subtract` changes the value wherever `-shape` is above
        // the field, which reaches beyond the sphere even though the *sign*
        // cannot change out there.
        for k in 0..LATTICE {
            for j in 0..LATTICE {
                for ii in 0..LATTICE {
                    let p = lattice_point(ii, j, k);
                    let idx = (k * LATTICE + j) * LATTICE + ii;
                    grid[idx] = apply(brush.op, grid[idx], brush.shape.sample(p));
                }
            }
        }
        log.push(brush);
    }
    log
}

/// Least squares through `(x, y)` pairs: returns `(intercept, slope)`.
fn fit(points: &[(f64, f64)]) -> (f64, f64) {
    let n = points.len() as f64;
    let mx = points.iter().map(|p| p.0).sum::<f64>() / n;
    let my = points.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx = points.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum::<f64>();
    let sxy = points
        .iter()
        .map(|p| (p.0 - mx) * (p.1 - my))
        .sum::<f64>();
    let slope = sxy / sxx;
    (my - slope * mx, slope)
}

fn median(mut v: Vec<f64>) -> f64 {
    assert!(!v.is_empty(), "median of nothing");
    v.sort_unstable_by(|a, b| a.partial_cmp(b).expect("finite"));
    v[v.len() / 2]
}

/// The `q`-quantile of an unsorted sample, or `0.0` for an empty one.
///
/// Empty means the comparison never ran, which the `cp_vs_log_checks` count
/// reports separately — so a zero here is readable rather than a silent default.
fn quantile(v: &[f64], q: f64) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let mut s = v.to_vec();
    s.sort_unstable_by(|a, b| a.partial_cmp(b).expect("finite"));
    let idx = ((s.len() - 1) as f64 * q).round() as usize;
    s[idx]
}

/// The largest cadence whose worst single undo fits the frame, from an affine
/// per-chunk cost curve and a worst-case dirty-chunk count.
fn cadence_from_curve(intercept: f64, slope: f64, chunks: u64) -> f64 {
    (FRAME_MS / chunks as f64 - intercept) / slope
}

/// One burst in the trace: the log length it starts from, and how deep it goes.
#[derive(Clone, Copy)]
struct Burst {
    at: usize,
    depth: usize,
    boundary: bool,
}

/// The trace's undo bursts at one cadence.
///
/// The regular bursts are cadence-independent so rows are comparable; the four
/// boundary bursts exist because at a coarse cadence no regular burst reaches a
/// checkpoint, and then the registered vacuity control would read zero.
fn bursts(k: usize) -> Vec<Burst> {
    let mut out: Vec<Burst> = Vec::new();
    let mut b = 0usize;
    let mut at = BURST_EVERY;
    while at <= TRACE_EDITS {
        out.push(Burst {
            at,
            depth: BURST_LENS[b % BURST_LENS.len()],
            boundary: false,
        });
        b += 1;
        at += BURST_EVERY;
    }
    for fifth in 1..=4usize {
        // A checkpoint sits at every multiple of `k`. Undoing from `c + 1` goes
        // `c+1 -> c` (the ladder's newest is `c`, still valid) and then
        // `c -> c-1`, which pops it. So depth 3 crosses exactly once, at any
        // cadence.
        let c = (TRACE_EDITS * fifth / 5) / k * k;
        if (BOUNDARY_BURST..TRACE_EDITS).contains(&c) {
            out.push(Burst {
                at: c + 1,
                depth: BOUNDARY_BURST,
                boundary: true,
            });
        }
    }
    // A boundary burst wins a collision with a regular one at the same length:
    // the regular schedule is comparability, the boundary one is the registered
    // vacuity control, and losing it silently is the failure the control is for.
    out.sort_unstable_by_key(|x| (x.at, !x.boundary, usize::MAX - x.depth));
    out.dedup_by_key(|x| x.at);
    out
}

/// What one (field, cadence) trace produced.
#[derive(Default)]
struct Trace {
    worst_undo_ms: f64,
    worst_undo_chunks: u64,
    worst_depth: u64,
    undos: u64,
    crossings: u64,
    crossings_regular: u64,
    hash_mismatches: u64,
    rebuilds: u64,
    rebuild_mismatches: u64,
    world_hash_checks: u64,
    cp_vs_log_checks: u64,
    cp_vs_log_hash_mismatches: u64,
    cp_vs_log_position_mismatches: u64,
    cp_vs_log_normal_mismatches: u64,
    /// Every angle between a checkpoint-path normal and its log-path twin.
    ///
    /// Kept in full rather than reduced to a maximum: a worst case over tens of
    /// thousands of vertices is one degenerate vertex, and "the worst normal
    /// flipped" and "every normal flipped" are different findings.
    normal_dev_deg: Vec<f64>,
    empty_dirty_edits: u64,
    peak_checkpoint_bytes: u64,
    checkpoints_built: u64,
    hash_sets_without_geometry: u64,
}

/// The interleaved-undo trace: `TRACE_EDITS` pushes, checkpoints every `k`, undo
/// bursts, and the golden hash checked across every one of them.
fn run_trace<F: Sdf<Scalar = f64>>(
    world: &World,
    base: &F,
    base_bake: &Checkpoint,
    log: &Log,
    k: usize,
) -> Trace {
    let mut mc = MarchingCubes::<f64>::new();
    let mut t = Trace::default();

    let schedule = bursts(k);
    // Which brush indices get undone, so the pre-edit hash is only paid for
    // those. `at` is a log *length*; undoing `depth` removes brushes
    // `at-depth … at-1`.
    let mut undone = vec![false; TRACE_EDITS];
    let mut is_boundary = vec![false; TRACE_EDITS];
    for b in &schedule {
        for j in (b.at - b.depth)..b.at {
            undone[j] = true;
            is_boundary[j] = b.boundary;
        }
    }

    // `(index, snapshot)`, ascending. Entry zero is the bake of the raw base
    // field and is never dropped.
    let mut ladder: Vec<(usize, Checkpoint)> = vec![(0, Checkpoint { values: base_bake.values.clone() })];
    let mut ladder_bytes = base_bake.bytes();
    t.peak_checkpoint_bytes = ladder_bytes;

    // Per undone brush: two regions and two hashes.
    //
    // `timed_of` is always the edit's **own dirty chunks** — M-50's quantity,
    // and the only thing `worst_undo_ms` may be measured over. `region_of` is
    // what the hash is taken over, which is periodically widened to the whole
    // world because a subtracting brush moves values outside its own bounding
    // box and a dirty-set hash could in principle miss that. Widening the hash
    // must not widen the clock: the first run reported a 27-chunk whole-world
    // re-mesh as a `worst_undo_ms`, which is a different quantity wearing
    // undo's name.
    let mut timed_of: Vec<Option<Vec<ChunkId>>> = vec![None; TRACE_EDITS];
    let mut region_of: Vec<Option<Vec<ChunkId>>> = vec![None; TRACE_EDITS];
    let mut hash_before: Vec<u64> = vec![0; TRACE_EDITS];
    let mut hash_after: Vec<u64> = vec![0; TRACE_EDITS];
    // Vertices that went into an undone edit's two hashes. Zero would mean the
    // hash was taken over empty meshes and could not have moved — M-44.
    let mut vertices_hashed: Vec<usize> = vec![0; TRACE_EDITS];

    let mut undone_seen = 0usize;
    let mut next_burst = 0usize;

    for length in 1..=TRACE_EDITS {
        let j = length - 1;
        let c = ladder.last().expect("ladder").0;
        let cp = &ladder.last().expect("ladder").1;
        let chunks = dirty_chunks(world, cp, &log[c..j], &log[c..length], brush_box(world, &log[j]));
        if chunks.is_empty() {
            t.empty_dirty_edits += 1;
        }

        if undone[j] {
            undone_seen += 1;
            let set = if undone_seen.is_multiple_of(WORLD_HASH_EVERY) {
                t.world_hash_checks += 1;
                world.all.clone()
            } else {
                chunks.clone()
            };
            let before = BrushStack { base: cp, brushes: &log[c..j] };
            let (h, verts_before) = hash_chunks(world, &mut mc, &before, &set);
            hash_before[j] = h;
            vertices_hashed[j] = verts_before;
            region_of[j] = Some(set);
            timed_of[j] = Some(chunks.clone());
        }

        // The edit is now applied: the state is `log[..length]`.
        if length.is_multiple_of(k) {
            let next = advance(cp, &log[c..length]);
            ladder_bytes += next.bytes();
            ladder.push((length, next));
            t.checkpoints_built += 1;
            t.peak_checkpoint_bytes = t.peak_checkpoint_bytes.max(ladder_bytes);
        }

        if undone[j] {
            let (c2, cp2) = ladder.last().expect("ladder");
            let after = BrushStack { base: cp2, brushes: &log[*c2..length] };
            let set = region_of[j].as_ref().expect("region").clone();
            let (h, verts_after) = hash_chunks(world, &mut mc, &after, &set);
            hash_after[j] = h;
            vertices_hashed[j] += verts_after;
            if vertices_hashed[j] == 0 {
                t.hash_sets_without_geometry += 1;
            }
        }

        if next_burst < schedule.len() && schedule[next_burst].at == length {
            let burst = schedule[next_burst];
            next_burst += 1;
            let mut dropped: Vec<(usize, Checkpoint)> = Vec::new();

            // Down: undo `depth` edits, each one timed and hashed.
            for step in 0..burst.depth {
                let target = length - step - 1;
                let j = target;
                while ladder.len() > 1 && ladder.last().expect("ladder").0 > target {
                    let gone = ladder.pop().expect("ladder");
                    ladder_bytes -= gone.1.bytes();
                    dropped.push(gone);
                    t.crossings += 1;
                    if !burst.boundary {
                        t.crossings_regular += 1;
                    }
                }
                let (c, cp) = ladder.last().expect("ladder");
                let field = BrushStack { base: cp, brushes: &log[*c..target] };
                let set = region_of[j].as_ref().expect("region");
                let timed = timed_of[j].as_ref().expect("timed");

                let ms = time_extract(world, &mut mc, &field, timed);
                t.undos += 1;
                if ms > t.worst_undo_ms {
                    t.worst_undo_ms = ms;
                    t.worst_undo_chunks = timed.len() as u64;
                }
                t.worst_depth = t.worst_depth.max((target - *c) as u64);

                if hash_chunks(world, &mut mc, &field, set).0 != hash_before[j] {
                    t.hash_mismatches += 1;
                }

                // One step per burst: a full fold from the base field at a log
                // length of several hundred is the most expensive thing here,
                // and the question it answers does not change within a burst.
                // Every boundary burst is checked at the step that crosses.
                if step == 0 && undone_seen.is_multiple_of(CP_VS_LOG_EVERY)
                    || is_boundary[j] && step == 1
                {
                    t.cp_vs_log_checks += 1;
                    let straight = BrushStack { base, brushes: &log[..target] };
                    if hash_chunks(world, &mut mc, &straight, set).0 != hash_before[j] {
                        t.cp_vs_log_hash_mismatches += 1;
                    }
                    for id in set {
                        let a = extract(world, &mut mc, &field, *id);
                        let b = extract(world, &mut mc, &straight, *id);
                        assert_eq!(
                            a.positions.len(),
                            b.positions.len(),
                            "checkpoint path changed the vertex count, which would make the \
                             normal comparison meaningless"
                        );
                        // Bit patterns, not values — "bit-exact" is what is
                        // being measured, and it is `mesh_hash`'s own rule.
                        for (pa, pb) in a.positions.iter().zip(&b.positions) {
                            if pa.map(f64::to_bits) != pb.map(f64::to_bits) {
                                t.cp_vs_log_position_mismatches += 1;
                            }
                        }
                        for (na, nb) in a.normals.iter().zip(&b.normals) {
                            if na.map(f64::to_bits) != nb.map(f64::to_bits) {
                                t.cp_vs_log_normal_mismatches += 1;
                            }
                            let dot = (na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2])
                                .clamp(-1.0, 1.0);
                            let deg = dot.acos().to_degrees();
                            t.normal_dev_deg.push(deg);
                        }
                    }
                }
            }

            // Up: redo them, rebuilding every dropped checkpoint as it is
            // reached, and requiring the hash to come back.
            for step in (0..burst.depth).rev() {
                let target = length - step;
                let j = target - 1;
                while let Some((idx, _)) = dropped.last() {
                    if *idx > target {
                        break;
                    }
                    let (idx, old) = dropped.pop().expect("dropped");
                    let (c, cp) = ladder.last().expect("ladder");
                    let rebuilt = advance(cp, &log[*c..idx]);
                    t.rebuilds += 1;
                    // Bit patterns, so a rebuilt `NaN` would still count as a
                    // match only if it were the same `NaN`.
                    let identical = rebuilt.values.len() == old.values.len()
                        && rebuilt
                            .values
                            .iter()
                            .zip(&old.values)
                            .all(|(a, b)| a.to_bits() == b.to_bits());
                    if !identical {
                        t.rebuild_mismatches += 1;
                    }
                    ladder_bytes += rebuilt.bytes();
                    ladder.push((idx, rebuilt));
                }
                let (c, cp) = ladder.last().expect("ladder");
                let field = BrushStack { base: cp, brushes: &log[*c..target] };
                let set = region_of[j].as_ref().expect("region");
                if hash_chunks(world, &mut mc, &field, set).0 != hash_after[j] {
                    t.hash_mismatches += 1;
                }
            }
            assert!(
                dropped.is_empty(),
                "a burst finished with checkpoints still torn down: the ladder did not come back"
            );
        }
    }
    t
}

/// The cost curve, the crossover sweep, and every cadence's trace, for one
/// field.
struct FieldResult {
    name: &'static str,
    buckets: [f64; 4],
    intercept: f64,
    slope: f64,
    worst_chunks_curve: u64,
    undo_ms: f64,
    refold_ms: Vec<f64>,
    sweep_chunks: u64,
    crossover_n: u64,
    refold_a: f64,
    refold_b: f64,
    crossover_fitted: f64,
    crossover_margin: f64,
    predicted_k: f64,
    measured_k: f64,
    traces: Vec<Trace>,
    ops_seen: [u64; 3],
    curve_empty_edits: u64,
    bucket_counts: Vec<usize>,
}

fn run_field<F: Sdf<Scalar = f64>>(name: &'static str, base: &F) -> FieldResult {
    let world = World::new();
    let mut mc = MarchingCubes::<f64>::new();
    let base_bake = bake(base);
    let log = &build_log(base);
    println!("  {name}: log of {} edits placed", log.len());

    // ── the cost curve, M-50's instrument on this machine ────────────────────
    //
    // Full fold from the base field: cost per re-meshed chunk against the number
    // of edits in the log, which is exactly what E-202 reported per edit.
    let mut per_chunk: Vec<(f64, f64)> = Vec::with_capacity(HEAD);
    let mut worst_chunks_curve = 0u64;
    let mut empty = 0u64;
    for length in 1..=HEAD {
        let j = length - 1;
        let chunks = dirty_chunks(
            &world,
            base,
            &log[..j],
            &log[..length],
            brush_box(&world, &log[j]),
        );
        if chunks.is_empty() {
            empty += 1;
            continue;
        }
        worst_chunks_curve = worst_chunks_curve.max(chunks.len() as u64);
        let field = BrushStack { base, brushes: &log[..length] };
        per_chunk.push((
            length as f64,
            time_extract(&world, &mut mc, &field, &chunks) / chunks.len() as f64,
        ));
    }
    assert_eq!(
        empty, 0,
        "{name}: an edit in the cost curve dirtied no chunk, so its time is the time to do \
         nothing — M-44"
    );
    // Each bucket's median must rest on enough edits to be a median. Fifteen
    // log lengths fall in each of M-50's four buckets, and every one of them
    // dirtied a chunk, so the count is fifteen — asserted rather than assumed,
    // because a bucket down to two or three samples is a number with a
    // median's name and a single reading's variance.
    let bucket_counts: Vec<usize> = [(1, 15), (16, 30), (31, 45), (46, 60)]
        .iter()
        .map(|(lo, hi)| {
            per_chunk
                .iter()
                .filter(|(l, _)| *l >= f64::from(*lo) && *l <= f64::from(*hi))
                .count()
        })
        .collect();
    assert!(
        bucket_counts.iter().all(|n| *n >= 10),
        "{name}: a bucket median rests on {bucket_counts:?} samples, and under ten is a single \
         reading wearing a median's name"
    );

    let bucket = |lo: f64, hi: f64| -> f64 {
        median(
            per_chunk
                .iter()
                .filter(|(l, _)| *l >= lo && *l <= hi)
                .map(|(_, ms)| *ms)
                .collect(),
        )
    };
    let buckets = [
        bucket(1.0, 15.0),
        bucket(16.0, 30.0),
        bucket(31.0, 45.0),
        bucket(46.0, 60.0),
    ];
    let mine: Vec<(f64, f64)> = M50
        .iter()
        .zip(&buckets)
        .map(|((mid, _), ms)| (*mid, *ms))
        .collect();
    let (intercept, slope) = fit(&mine);
    let (m50_a, m50_b) = fit(&M50);

    // ── the crossover sweep at head HEAD ────────────────────────────────────
    //
    // Undoing the last edit reaches length `HEAD - 1`. The log path walks the
    // base field and every remaining brush; the checkpoint path walks a
    // trilinear lookup and the `d` brushes above the checkpoint.
    let target = HEAD - 1;
    let chunks = dirty_chunks(
        &world,
        base,
        &log[..target],
        &log[..HEAD],
        brush_box(&world, &log[target]),
    );
    assert!(
        !chunks.is_empty(),
        "{name}: the crossover sweep's edit dirtied no chunk — M-44"
    );

    // Monomorphised, not `&dyn Sdf`: a virtual call per sample would be paid by
    // both arms and so would not bias their ratio, but it would put the sweep's
    // absolute milliseconds on a different footing from the cost curve above,
    // and C2 compares the two.
    let undo_ms = time_extract(
        &world,
        &mut mc,
        &BrushStack { base, brushes: &log[..target] },
        &chunks,
    );

    // Checkpoints at every `c` in `0 ..= target`, built incrementally — which is
    // itself the bit-exact ladder the trace uses.
    let mut refold_ms = vec![0.0f64; target + 1];
    let mut cp = Checkpoint { values: base_bake.values.clone() };
    for c in 0..=target {
        if c > 0 {
            cp = advance(&cp, &log[c - 1..c]);
        }
        let d = target - c;
        refold_ms[d] = time_extract(
            &world,
            &mut mc,
            &BrushStack { base: &cp, brushes: &log[c..target] },
            &chunks,
        );
    }
    drop(cp);

    // The registered N: the smallest separation at which undoing beats
    // re-folding. Zero for none inside the log.
    let crossover_n = (0..=target)
        .find(|d| undo_ms < refold_ms[*d])
        .map_or(0u64, |d| d as u64);

    // **The registered `N` sits at the boundary of any search that can be run,
    // and two runs of this harness disagreed about it (125 and none).** That is
    // `P-72`'s defect — a boundary minimum reported as an optimum — so the
    // question is answered from the whole sweep instead of its last two points.
    //
    // Fit `refold(d) = A + B·d` over all `target + 1` measurements. `B` is the
    // marginal cost of one more edit above the checkpoint, measured over 128
    // points rather than four buckets. The crossover in that model is
    // `d* = (undo − A)/B`, and `d* − target` is how many edits **past the
    // bottom of the log** a checkpoint would have to sit before undoing became
    // cheaper. A positive margin means no crossover can exist at any log
    // length, because `d ≤ target` always — the checkpoint cannot be below the
    // base field. It is a ratio of two measured times divided by a third, so it
    // survives a governed clock in a way milliseconds do not (`M-280`).
    let sweep: Vec<(f64, f64)> = (0..=target).map(|d| (d as f64, refold_ms[d])).collect();
    let (refold_a, refold_b) = fit(&sweep);
    let crossover_fitted = (undo_ms - refold_a) / refold_b;
    let crossover_margin = crossover_fitted - target as f64;

    let predicted_k = cadence_from_curve(m50_a, m50_b, worst_chunks_curve);
    let measured_k = cadence_from_curve(intercept, slope, worst_chunks_curve);

    let mut ops_seen = [0u64; 3];
    for b in log.iter() {
        match b.op {
            BrushOp::Add => ops_seen[0] += 1,
            BrushOp::Subtract => ops_seen[1] += 1,
            BrushOp::SmoothAdd { .. } => ops_seen[2] += 1,
        }
    }

    let traces: Vec<Trace> = CADENCES
        .iter()
        .map(|k| {
            let t = run_trace(&world, base, &base_bake, log, *k);
            println!(
                "  {name} k={k:<4} worst_undo {:.4} ms over {} chunks, depth {}, crossings {} \
                 (regular {}), hash mismatches {}, cp-vs-log {}/{}",
                t.worst_undo_ms,
                t.worst_undo_chunks,
                t.worst_depth,
                t.crossings,
                t.crossings_regular,
                t.hash_mismatches,
                t.cp_vs_log_hash_mismatches,
                t.cp_vs_log_checks
            );
            t
        })
        .collect();

    FieldResult {
        name,
        buckets,
        intercept,
        slope,
        worst_chunks_curve,
        undo_ms,
        refold_ms,
        sweep_chunks: chunks.len() as u64,
        crossover_n,
        refold_a,
        refold_b,
        crossover_fitted,
        crossover_margin,
        predicted_k,
        measured_k,
        traces,
        ops_seen,
        curve_empty_edits: empty,
        bucket_counts,
    }
}

/// The CPU's reported clock. `M-280`: on a governed CPU a nanosecond is not a
/// unit unless the clock is on the row.
fn cpu_mhz() -> String {
    std::fs::read_to_string("/proc/cpuinfo")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("cpu MHz"))
                .and_then(|l| l.split(':').nth(1))
                .map(|v| v.trim().to_string())
        })
        .unwrap_or_else(|| String::from("unknown"))
}

type Row = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|a| a == "--bench") {
        return;
    }

    let (m50_a, m50_b) = fit(&M50);
    println!(
        "M-50 fit over its four bucket midpoints: cost(L) = {m50_a:.5} + {m50_b:.7}·L ms per chunk"
    );
    println!(
        "SHARE, before the run: C1's registered direction is unreachable by monotonicity — \
         refold(d) is non-decreasing in d and undo does not depend on d, so \
         {{d : undo < refold(d)}} is an up-set and cannot be a 'fewer than N' threshold."
    );
    println!(
        "SHARE, before the run: C2 is reachable — the fixed per-chunk term alone only fills \
         {FRAME_MS} ms at {:.0} dirty chunks, and this world holds {}.",
        FRAME_MS / m50_a,
        (WORLD_CELLS / CHUNK_CELLS).pow(3)
    );
    println!();

    // Three base fields whose own sample costs differ by roughly an order of
    // magnitude — one `sqrt`, six transcendentals, four octaves of lattice
    // noise. C1's crossover condition is `base + c·brush < lookup`, so if a
    // crossover depends on the base field's cost at all, it moves across these.
    // Each field's surface passes through the `±1.5` box the brushes are drawn
    // in, which is what lets every edit move the surface.
    let results = vec![
        run_field("sphere_r1_2", &Sphere::<f64> { center: [0.0; 3], radius: 1.2 }),
        run_field("gyroid", &Gyroid::<f64>::canonical()),
        run_field("fbm_terrain", &FbmTerrain::<f64>::canonical()),
    ];

    // C1's stability clause is a statement about the set of crossovers, so it
    // can only be decided once every field has answered.
    let crossovers: Vec<u64> = results.iter().map(|r| r.crossover_n).collect();
    let all_same = crossovers.windows(2).all(|w| w[0] == w[1]);
    let stable = if all_same && crossovers[0] == 0 {
        "none_on_all_fields"
    } else if all_same {
        "true"
    } else {
        "false"
    };
    // "Undo strictly less than refold whenever FEWER than N edits separate
    // them" needs a non-empty down-set. `crossover_n == 0` is the empty up-set:
    // undo never wins, at any separation.
    let c1_holds = crossovers.iter().all(|n| *n > 0) && all_same;

    let mhz = cpu_mhz();
    let mut rows: Vec<Row> = Vec::new();

    for r in &results {
        assert!(
            r.ops_seen.iter().all(|n| *n > 0),
            "{}: the log is missing an operation, so C3 would be testing a history that \
             reorders for free",
            r.name
        );
        let ratios: Vec<f64> = M50
            .iter()
            .zip(&r.buckets)
            .map(|((_, m50_ms), ms)| ms / m50_ms)
            .collect();
        let ratio_spread = ratios.iter().copied().fold(f64::MIN, f64::max)
            / ratios.iter().copied().fold(f64::MAX, f64::min);
        let k_ratio = if r.predicted_k > 0.0 && r.measured_k > 0.0 {
            (r.predicted_k / r.measured_k).max(r.measured_k / r.predicted_k)
        } else {
            f64::INFINITY
        };

        // C2 is one claim about a derived `k`, not a claim per cadence. The
        // budget check runs at the coarsest swept cadence the derivation
        // permits: that is the trace the prediction actually endorses.
        let budget_idx = CADENCES
            .iter()
            .rposition(|k| *k as f64 <= r.measured_k);
        let budget_cadence = budget_idx.map_or(0, |i| CADENCES[i]);
        let budget_worst_ms = budget_idx.map_or(f64::INFINITY, |i| r.traces[i].worst_undo_ms);
        let c2_holds = k_ratio <= 2.0 && budget_worst_ms < FRAME_MS;

        // And the answer with no fitting in it at all: the coarsest swept
        // cadence whose trace actually kept the worst undo inside the frame.
        // `measured_k` above is the parallel of `predicted_k_from_m50` — four
        // bucket medians through one formula — and this is the observation it is
        // trying to predict.
        let trace_idx = CADENCES
            .iter()
            .enumerate()
            .rposition(|(i, _)| r.traces[i].worst_undo_ms < FRAME_MS);
        let empirical_k = trace_idx.map_or(0, |i| CADENCES[i]);
        let k_ratio_trace = if empirical_k > 0 {
            r.predicted_k / empirical_k as f64
        } else {
            f64::INFINITY
        };

        for (idx, k) in CADENCES.iter().enumerate() {
            let t = &r.traces[idx];
            // The deepest single undo a cadence of `k` can force is `k - 1`
            // edits above the checkpoint.
            let d = (k - 1).min(r.refold_ms.len() - 1);

            // `undos_crossing_checkpoint` is the registered vacuity control.
            // A zero here would mean C3 never left the cheap path.
            assert!(
                t.crossings > 0,
                "{} k={k}: no undo crossed a checkpoint boundary, so C3 tested only the cheap \
                 path — the registered vacuity control",
                r.name
            );
            assert_eq!(
                t.rebuild_mismatches, 0,
                "{} k={k}: a checkpoint rebuilt after being dropped differs from the one that \
                 was dropped, so the ladder is not a function of the log",
                r.name
            );
            assert_eq!(
                t.empty_dirty_edits, 0,
                "{} k={k}: an edit dirtied no chunk — M-44",
                r.name
            );
            assert_eq!(
                t.hash_sets_without_geometry, 0,
                "{} k={k}: an undone edit's hash was taken over a set of empty meshes, whose \
                 hash is a constant, so the C3 check over it could not have failed — M-44",
                r.name
            );

            // C3: the golden hash must come back across an undo and a redo that
            // tore the checkpoint ladder down and rebuilt it.
            let c3_holds = t.hash_mismatches == 0;

            rows.push(vec![
                ("field", r.name.to_string()),
                ("edits", TRACE_EDITS.to_string()),
                ("checkpoint_every_k", k.to_string()),
                ("undo_ms", format!("{:.6}", r.undo_ms)),
                ("refold_ms", format!("{:.6}", r.refold_ms[d])),
                ("crossover_n", r.crossover_n.to_string()),
                ("crossover_stable_across_fields", stable.to_string()),
                ("predicted_k_from_m50", format!("{:.1}", r.predicted_k)),
                ("measured_k", format!("{:.1}", r.measured_k)),
                ("worst_undo_ms", format!("{:.6}", t.worst_undo_ms)),
                (
                    "undos_crossing_checkpoint",
                    t.crossings.to_string(),
                ),
                ("hash_mismatches", t.hash_mismatches.to_string()),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_holds.to_string()),
                ("c3_holds", c3_holds.to_string()),
                // ── instrument and provenance ───────────────────────────────
                ("cpu_mhz", mhz.clone()),
                ("scalar", "f64".to_string()),
                ("world_cells", WORLD_CELLS.to_string()),
                ("chunk_cells", CHUNK_CELLS.to_string()),
                ("cell_size", format!("{CELL_SIZE}")),
                ("curve_edits", HEAD.to_string()),
                ("reps", REPS.to_string()),
                // ── C1's other side ─────────────────────────────────────────
                ("refold_ms_at_d0", format!("{:.6}", r.refold_ms[0])),
                (
                    "refold_ms_at_dmax",
                    format!("{:.6}", r.refold_ms[r.refold_ms.len() - 1]),
                ),
                (
                    "undo_over_refold_at_k",
                    format!("{:.4}", r.undo_ms / r.refold_ms[d]),
                ),
                ("crossover_search_dmax", (r.refold_ms.len() - 1).to_string()),
                ("crossover_sweep_chunks", r.sweep_chunks.to_string()),
                ("refold_fit_intercept_ms", format!("{:.6}", r.refold_a)),
                ("refold_fit_slope_ms_per_edit", format!("{:.6}", r.refold_b)),
                ("crossover_n_fitted", format!("{:.2}", r.crossover_fitted)),
                (
                    "crossover_margin_edits_past_log",
                    format!("{:.2}", r.crossover_margin),
                ),
                ("refold_ms_at_d1", format!("{:.6}", r.refold_ms[1])),
                ("refold_ms_at_d2", format!("{:.6}", r.refold_ms[2])),
                ("refold_ms_at_d4", format!("{:.6}", r.refold_ms[4])),
                ("refold_ms_at_d8", format!("{:.6}", r.refold_ms[8])),
                ("refold_ms_at_d16", format!("{:.6}", r.refold_ms[16])),
                ("refold_ms_at_d32", format!("{:.6}", r.refold_ms[32])),
                ("refold_ms_at_d64", format!("{:.6}", r.refold_ms[64])),
                (
                    "crossover_direction",
                    if r.crossover_n > 0 {
                        "undo_wins_above".to_string()
                    } else {
                        "none_undo_never_wins".to_string()
                    },
                ),
                // ── C2's other side ────────────────────────────────────────
                ("bucket1_ms_per_chunk", format!("{:.6}", r.buckets[0])),
                ("bucket2_ms_per_chunk", format!("{:.6}", r.buckets[1])),
                ("bucket3_ms_per_chunk", format!("{:.6}", r.buckets[2])),
                ("bucket4_ms_per_chunk", format!("{:.6}", r.buckets[3])),
                ("bucket1_ratio_to_m50", format!("{:.4}", ratios[0])),
                ("bucket2_ratio_to_m50", format!("{:.4}", ratios[1])),
                ("bucket3_ratio_to_m50", format!("{:.4}", ratios[2])),
                ("bucket4_ratio_to_m50", format!("{:.4}", ratios[3])),
                ("bucket_ratio_spread", format!("{ratio_spread:.4}")),
                ("fit_intercept_ms", format!("{:.6}", r.intercept)),
                ("fit_slope_ms_per_edit", format!("{:.8}", r.slope)),
                ("m50_fit_intercept_ms", format!("{m50_a:.6}")),
                ("m50_fit_slope_ms_per_edit", format!("{m50_b:.8}")),
                ("k_ratio", format!("{k_ratio:.4}")),
                ("c2_budget_cadence", budget_cadence.to_string()),
                ("c2_budget_worst_undo_ms", format!("{budget_worst_ms:.6}")),
                ("measured_k_from_trace", empirical_k.to_string()),
                ("k_ratio_m50_over_trace", format!("{k_ratio_trace:.4}")),
                (
                    "worst_dirty_chunks_curve",
                    r.worst_chunks_curve.to_string(),
                ),
                (
                    "formula_worst_undo_ms",
                    format!(
                        "{:.6}",
                        r.worst_chunks_curve as f64 * (r.intercept + r.slope * (k - 1) as f64)
                    ),
                ),
                ("worst_undo_chunks", t.worst_undo_chunks.to_string()),
                ("worst_refold_depth", t.worst_depth.to_string()),
                // ── C3's other side, and the vacuity controls ───────────────
                ("undos_total", t.undos.to_string()),
                (
                    "crossings_regular_only",
                    t.crossings_regular.to_string(),
                ),
                ("checkpoints_built", t.checkpoints_built.to_string()),
                ("checkpoint_rebuilds", t.rebuilds.to_string()),
                (
                    "checkpoint_rebuild_mismatches",
                    t.rebuild_mismatches.to_string(),
                ),
                ("world_hash_checks", t.world_hash_checks.to_string()),
                ("empty_dirty_edits", t.empty_dirty_edits.to_string()),
                ("curve_empty_dirty_edits", r.curve_empty_edits.to_string()),
                (
                    "bucket_sample_counts",
                    r.bucket_counts
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join("+"),
                ),
                (
                    "hash_sets_without_geometry",
                    t.hash_sets_without_geometry.to_string(),
                ),
                ("cp_vs_log_checks", t.cp_vs_log_checks.to_string()),
                (
                    "cp_vs_log_hash_mismatches",
                    t.cp_vs_log_hash_mismatches.to_string(),
                ),
                (
                    "cp_vs_log_position_mismatches",
                    t.cp_vs_log_position_mismatches.to_string(),
                ),
                (
                    "cp_vs_log_normal_mismatches",
                    t.cp_vs_log_normal_mismatches.to_string(),
                ),
                (
                    "cp_vs_log_normals_compared",
                    t.normal_dev_deg.len().to_string(),
                ),
                (
                    "normal_dev_p50_deg",
                    format!("{:.6}", quantile(&t.normal_dev_deg, 0.50)),
                ),
                (
                    "normal_dev_p99_deg",
                    format!("{:.6}", quantile(&t.normal_dev_deg, 0.99)),
                ),
                (
                    "max_normal_dev_deg",
                    format!("{:.6}", quantile(&t.normal_dev_deg, 1.0)),
                ),
                (
                    "peak_checkpoint_bytes",
                    t.peak_checkpoint_bytes.to_string(),
                ),
                ("log_add", r.ops_seen[0].to_string()),
                ("log_subtract", r.ops_seen[1].to_string()),
                ("log_smooth_add", r.ops_seen[2].to_string()),
            ]);
        }
    }

    common::experiment::run(isomesh::experiment!("P-95"), |run| {
        for row in rows {
            run.record(&row);
        }
    });
}
