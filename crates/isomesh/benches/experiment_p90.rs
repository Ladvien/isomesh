//! **P-90 — per-brick edit-list culling, at a brick of sixty-four cells.**
//!
//! Ticket: R-090. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p90
//! ```
//!
//! Writes `docs/experiments/p-90.csv`.
//!
//! # Hypothesis, as registered
//!
//! Dreams' hierarchical evaluator culls over 99% of naive edit evaluations by
//! maintaining a per-block edit list (Evans, *Learning from Failure*, SIGGRAPH
//! 2015 Advances in Real-Time Rendering — 1 to 100,000 edits per sculpture,
//! 10–100 M voxels evaluated per second on a PS4, culling efficiency over 99%
//! against brute force; industry slides, no DOI, cited as testimony and not as
//! measurement). This crate's equivalent is `P-39`'s Lipschitz brush pruning at
//! 3.36× median, measured at **chunk** granularity — 32 cells. `M-377` has since
//! moved the optimum chunk to **4³, sixty-four cells**. The question the two
//! results create together is whether pruning still pays when the brick is that
//! small, or whether the per-brick bookkeeping eats it.
//!
//! - **C1.** At 4³ granularity, Lipschitz pruning still gives at least **2×** on
//!   `M-50`'s 46–60 brush bucket. *Falsified by* under 2×, which would put
//!   `P-72`'s optimum and `P-39`'s win in tension.
//! - **C2.** Carrying a per-brick surviving-brush list costs under **8 bytes per
//!   brick** amortised across the dig trace. *Falsified by* above 8 — a 128³
//!   world is 32,768 bricks, so this is a real memory line item.
//! - **C3.** The two culls **compose** rather than overlap: chunk-level then
//!   brick-level pruning removes strictly more than **either alone**. *Falsified
//!   by* no additional removal.
//!
//! **VACUITY CONTROL, as registered:** the trace must contain bricks whose
//! surviving-brush set differs from their parent chunk's, reported as
//! `bricks_differing_from_parent`, or C3 cannot fire. It is asserted non-zero.
//!
//! # Why `additional_removed` is measured against the *maximum* of the two
//! single-level culls, and not against the chunk level
//!
//! A brick's survivor list is produced by pruning its parent chunk's list, so it
//! is a subsequence of it: if `bricks_differing_from_parent` is non-zero then
//! "chunk-then-brick removes more than chunk alone" is **implied by the vacuity
//! assertion itself**. Scoring C3 on that comparison would be a HELD with no
//! instrument behind it — `P-70`'s C3, the weakest row in Phase 23. The clause
//! says *"strictly more than either alone"*, so the registered column is
//!
//! ```text
//! additional_removed = removed_both − max(removed_chunk_only, removed_brick_only)
//! ```
//!
//! and `removed_brick_only` is a real fourth arm: the **full** tape pruned
//! directly over the brick box, with no chunk pass in front of it. Both
//! one-sided differences are also columns (`additional_over_chunk_only`,
//! `additional_over_brick_only`) so a reader can see which comparison each
//! number belongs to.
//!
//! # The SHARE line, recomputed before the harness was written
//!
//! The registration says C1 moves *"the field-evaluation share of remesh, which
//! `M-377` measured as essentially all of `remesh_ms`"*. Pruning removes brush
//! evaluations and **nothing else** — the base field and the marching-cubes cell
//! work are untouchable — so the share C1 can actually move is
//!
//! ```text
//! share = (t_remesh_full − t_remesh_floor) / t_remesh_full
//! ```
//!
//! where the floor is the same dirty bricks re-meshed from the base field with
//! an empty tape. That is a measured arm (`remesh_floor_ms`), not an assumption,
//! and `ceiling_speedup = 1/(1 − share)` is in the file beside the speedups.
//! Pre-run arithmetic, at ≈4 ns for a sphere or capsule distance against ≈25 ns
//! for `gyroid`'s six transcendentals: the 46–60 bucket gives `share ≈ 0.906`
//! and a ceiling of **≈10.6×** on the harder of the two bases. **C1's 2× is
//! arithmetically reachable with a factor of five to spare**, which is the
//! `✗51` check this phase requires and it passes.
//!
//! # Why the base fields are `gyroid` and `ground_slab`, and not `M-377`'s pair
//!
//! `FbmTerrain::value_bound()` is `FieldBound::Unbounded` — *"a heightfield's
//! value is a vertical distance, and fbm has no slope bound worth asserting"*.
//! With no Lipschitz constant for the base there is no enclosure to test a brush
//! against, so **the mechanism does not exist on `fbm_terrain` at any
//! granularity**: every brush survives by construction, and a timed arm would
//! report the cost of the bookkeeping against a benefit that is structurally
//! zero. That is an availability finding, asserted here rather than measured,
//! and the two bases that are measured are chosen to bracket the quantity that
//! decides the answer — the **base Lipschitz constant**:
//!
//! | base | `value_bound()` | `l` | surface |
//! |---|---|---:|---|
//! | `gyroid` | `Lipschitz { l: 2√3·scale }` | 3.4641 | everywhere |
//! | `ground_slab` (`BoxExact`) | `Exact` | 1.0000 | one sheet at `y = 0` |
//!
//! `l` is the whole mechanism: an `Add` is prunable over a box of circumradius
//! `r` only when the brush clears the running field by `(1 + l)·r`, so the
//! threshold is **0.4833 world units at 4³ and 3.8660 at 32³** on `gyroid`, and
//! `gyroid`'s own value range is ±3. The prediction that follows before any
//! timing: chunk-level pruning on `gyroid` prunes close to nothing, and brick
//! level prunes almost everything.
//!
//! # The fixture
//!
//! A 128³-cell world over `[-2, 2]³`, partitioned twice: **32,768 bricks of 4³**
//! and **64 chunks of 32³**, eight bricks to a chunk edge. The edit log is
//! `M-50`'s four buckets — 15, 30, 45 and 60 brushes — as nested prefixes of one
//! generated log, so the four rows of a field are a sweep over log length on the
//! same log and not four unrelated fixtures.
//!
//! Each log brush is a sphere or a capsule, `Add` or `Subtract`, placed **on the
//! base field's own surface** at a random column, jittered vertically by up to
//! two cells. Uniform-over-the-volume placement (P-39's choice) would leave most
//! brushes in open air or deep rock, where distance prunes them trivially; a
//! surface-following log is the harder fixture, because every brush is in the
//! band where the bricks that need re-meshing live and only *lateral* distance
//! can reject it — which is exactly the quantity a smaller box resolves better.
//!
//! The last eleven entries of every log are `P-72`'s dig trace inherited: eleven
//! spherical `Subtract`s of radius six cells along a straight path across the
//! world, each one a separate `mark_edit` cycle, the height probed per edit at
//! that edit's own `x`. So a bucket of 60 runs the trace with the log growing
//! from 50 to 60 brushes, and every step of every arm stays inside its bucket.
//!
//! # Controls
//!
//! - **The registered vacuity control**, `bricks_differing_from_parent`,
//!   asserted non-zero.
//! - **`M-44`: every edit must dirty at least one brick.** A step that marked
//!   nothing would contribute a fast time for doing nothing.
//! - **Soundness, and it is the control the speedups depend on.** Every brick is
//!   re-meshed from the full tape, from its chunk's survivors and from its
//!   brick's survivors, and all three must be **bit-identical** — positions and
//!   normals compared as bits, because `-0.0 == 0.0`, plus `mesh_hash`. A pruner
//!   that removed one brush too many would otherwise report a speedup for a
//!   wrong mesh.
//! - **The sample counts of the pruned and unpruned arms must agree.** They are
//!   counted through a wrapper, not derived: identical meshes must have made
//!   identical `Sdf::sample` sequences, and a disagreement would mean the
//!   comparison is not like-for-like.
//! - **Integers beside every clock.** `M-280`: on a governed CPU a nanosecond is
//!   not a unit. `shape_evals_*` and `eval_ratio_*` are exact counts of brush
//!   evaluations, reproduce on any machine, and are what the mechanism is
//!   actually about; the millisecond columns are the same experiment with a
//!   clock on it, all four arms interleaved within each repetition (`M-281`).

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

mod common;

use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet};
use std::hint::black_box;
use std::time::Instant;

use isomesh::brush::{Brush, BrushOp, BrushStack, Capsule};
use isomesh::chunk::dirty::{DirtySet, mark_edit};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::{BoundedSdf, BoxExact, FbmTerrain, FieldBound, Gyroid, Sphere};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf, Shape3};

/// Cells per axis in the world. `M-377`'s fixed cell count, so the granularity
/// numbers here are comparable with `p-72.csv` row for row.
const WORLD_CELLS: u32 = 128;

/// The brick: `M-377`'s optimum, sixty-four cells.
const BRICK_CELLS: u32 = 4;

/// The chunk: `P-39`'s granularity, the one the 3.36× was measured at.
const CHUNK_CELLS: u32 = 32;

/// World extent per axis, so the world is `[-2, 2]³`.
const EXTENT: f64 = 4.0;

/// World origin, centred on the reference fields' own domain centre.
const ORIGIN: f64 = -EXTENT * 0.5;

/// World units per cell.
const CELL_SIZE: f64 = EXTENT / 128.0;

/// Bricks along each axis, and chunks along each axis.
const BRICKS_PER_AXIS: u32 = WORLD_CELLS / BRICK_CELLS;
const CHUNKS_PER_AXIS: u32 = WORLD_CELLS / CHUNK_CELLS;

/// Bricks in the world. The denominator of every C2 figure, and the number the
/// registration's falsifier names.
const BRICKS_TOTAL: u64 = (BRICKS_PER_AXIS as u64).pow(3);

/// Bricks along one chunk edge.
const BRICKS_PER_CHUNK_AXIS: i32 = (CHUNK_CELLS / BRICK_CELLS) as i32;

/// `M-50`'s four log buckets, named by their top: 1–15, 16–30, 31–45, 46–60.
/// C1 is registered on the last one.
const BUCKETS: [usize; 4] = [15, 30, 45, 60];

/// Edits in the dig trace, inherited from `P-72`.
const EDITS: usize = 11;

/// Log brushes generated once per field. The largest bucket needs
/// `60 − 11 = 49`, and the smaller buckets are prefixes of the same log.
const LOG_MAX: usize = 49;

/// Dig brush radius in cells, inherited from `P-72`.
const DIG_CELLS: f64 = 6.0;

/// Repetitions of the whole trace per arm.
///
/// **Twenty-one, and the first run at nine is why.** The per-repetition ratio
/// on the registered bucket ranged 3.94 to 9.83 around a median of 6.35 — a
/// 2.5× spread, on a machine that is not quiet. `M-337`'s re-audit is the
/// precedent: a registered 1.25× floor that re-measured at 1.022 three runs
/// later. The wall clock cannot be the load-bearing number here, which is why
/// `shape_evals_*` and `eval_ratio_*` are exact integers and reproduce
/// bit-for-bit between runs; `speedup_both_min` and `speedup_both_max` put the
/// remaining spread on the row rather than in a footnote.
const REPS: usize = 21;

/// ULP of slack added to every enclosure bound, as `P-39`.
const PAD_ULPS: f64 = 16.0;

// ── the fields ──────────────────────────────────────────────────────────────

/// The base field a log is applied to.
///
/// One enum rather than a generic parameter so both bases go through the same
/// monomorphisation of every timed loop: `M-281` compares within one build, and
/// two instantiations of the same loop are not one build in the sense that
/// matters.
#[derive(Clone, Copy, Debug)]
enum Base {
    /// Surface everywhere. `l = 2√3`, the hard case for the bound.
    Gyroid(Gyroid<f64>),
    /// Solid below `y = 0` and open above it: the exact-distance stand-in for
    /// `fbm_terrain`, which declares no Lipschitz constant at all.
    Ground(BoxExact<f64>),
}

impl Sdf for Base {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Gyroid(g) => g.sample(p),
            Self::Ground(b) => b.sample(p),
        }
    }
}

impl BoundedSdf for Base {
    fn value_bound(&self) -> FieldBound {
        match self {
            Self::Gyroid(g) => g.value_bound(),
            Self::Ground(b) => b.value_bound(),
        }
    }
}

/// A brush shape. Both variants are exact distance fields and declare it
/// through the crate's own [`BoundedSdf`] rather than by a constant written
/// here.
#[derive(Clone, Copy, Debug)]
enum Shape {
    Sphere(Sphere<f64>),
    Capsule(Capsule<f64>),
}

impl Sdf for Shape {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Sphere(s) => s.sample(p),
            Self::Capsule(c) => c.sample(p),
        }
    }
}

impl BoundedSdf for Shape {
    fn value_bound(&self) -> FieldBound {
        match self {
            Self::Sphere(s) => s.value_bound(),
            Self::Capsule(c) => c.value_bound(),
        }
    }
}

/// A field wrapper that counts `sample` calls, as `P-72`'s.
///
/// `BrushStack` does not override `Sdf::gradient`, so the six central-difference
/// samples per normal are counted here exactly as the extractor makes them.
struct Counted<'a, F> {
    field: &'a F,
    calls: &'a Cell<u64>,
}

impl<F: Sdf<Scalar = f64>> Sdf for Counted<'_, F> {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        self.calls.set(self.calls.get() + 1);
        self.field.sample(p)
    }
}

/// A 64-bit LCG, so the log is the same log on every machine and every run.
struct Lcg(u64);

impl Lcg {
    const fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) as u32
    }

    /// A float in `[0, 1)`, from 24 bits so it is exactly representable.
    fn unit(&mut self) -> f64 {
        f64::from(self.next_u32() >> 8) / 16_777_216.0
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

// ── the bound ───────────────────────────────────────────────────────────────

/// The enclosure of a scalar field over a box.
#[derive(Clone, Copy, Debug)]
struct Interval {
    lo: f64,
    hi: f64,
}

/// The box a brick's or a chunk's field evaluations can touch.
#[derive(Clone, Copy, Debug)]
struct Cube {
    centre: [f64; 3],
    /// Circumradius: half the diagonal of the sampled extent, plus the margin
    /// `Sdf::gradient`'s central differences reach outside it.
    radius: f64,
}

impl Cube {
    fn new(origin: [f64; 3], span: f64) -> Self {
        let centre = [
            origin[0] + span * 0.5,
            origin[1] + span * 0.5,
            origin[2] + span * 0.5,
        ];
        let mut far = 1.0f64;
        for lo in origin {
            far = far.max(lo.abs()).max((lo + span).abs());
        }
        let margin = <f64 as Real>::DIFF_STEP * far;
        let half = span * 0.5 + margin;
        Self {
            centre,
            radius: half * 3.0f64.sqrt(),
        }
    }
}

/// Slack for one bound, in absolute units.
fn pad(value: f64, reach: f64) -> f64 {
    PAD_ULPS * f64::EPSILON * (value.abs() + reach)
}

/// `f(centre) ± l·r`, widened so it is an enclosure and not an estimate.
fn enclose<S: BoundedSdf<Scalar = f64>>(field: &S, cube: Cube) -> Interval {
    let l = field
        .value_bound()
        .lipschitz()
        .expect("every field in this fixture declares a Lipschitz constant");
    let value = field.sample(cube.centre);
    let reach = l * cube.radius;
    let slack = reach + pad(value, reach);
    Interval {
        lo: value - slack,
        hi: value + slack,
    }
}

/// Select the brushes that can still change the fold anywhere in `cube`.
///
/// `P-39`'s `Policy::Sound`, unchanged: `min` and `max` **select** an operand
/// rather than computing a new value and negation is exact, so deleting a
/// provably-losing `Add` or `Subtract` moves the result by exactly zero ULP. The
/// inequality is strict because `f64::min` may return either operand when they
/// compare equal, which is observable for `+0.0` against `-0.0`. A `SmoothAdd`
/// is never pruned — `b + (a − b)` is not bit-identical to `a` — and this
/// fixture builds none, which `assert_hard_ops` checks rather than assumes.
///
/// Order is preserved: `Add` and `Subtract` do not commute with each other.
fn prune_into(
    tape: &[Brush<Shape>],
    base: &Base,
    cube: Cube,
    out: &mut Vec<Brush<Shape>>,
) -> usize {
    out.clear();
    let mut v = enclose(base, cube);
    for brush in tape {
        let s = enclose(&brush.shape, cube);
        match brush.op {
            BrushOp::Add => {
                if s.lo > v.hi {
                    continue;
                }
                v = Interval {
                    lo: v.lo.min(s.lo),
                    hi: v.hi.min(s.hi),
                };
            }
            BrushOp::Subtract => {
                let n = Interval {
                    lo: -s.hi,
                    hi: -s.lo,
                };
                if n.hi < v.lo {
                    continue;
                }
                v = Interval {
                    lo: v.lo.max(n.lo),
                    hi: v.hi.max(n.hi),
                };
            }
            BrushOp::SmoothAdd { k } => {
                let lo = v.lo.min(s.lo) - 0.25 * k;
                v = Interval {
                    lo: lo - pad(lo, 0.25 * k),
                    hi: v.hi.min(s.hi),
                };
            }
        }
        out.push(*brush);
    }
    out.len()
}

/// The fixture contains no `SmoothAdd`, so the unprunable branch above is a
/// soundness rule and not a silent arm of the measurement.
fn assert_hard_ops(tape: &[Brush<Shape>]) {
    for b in tape {
        assert!(
            matches!(b.op, BrushOp::Add | BrushOp::Subtract),
            "the fixture must be Add/Subtract only: a SmoothAdd never prunes and \
             would enter the survivor counts as an unprunable brush without saying so"
        );
    }
}

// ── the world ───────────────────────────────────────────────────────────────

/// The two partitions of one 128³-cell world.
///
/// `SAMPLES_PER_BRICK` is derived from the layout rather than written down, and
/// checked against `(cells + 1)³`: a brick that sampled a different number of
/// points would make every `shape_evals_*` count wrong by the same factor and
/// leave the ratios looking right.
struct World {
    brick: ChunkLayout<f64>,
    chunk: ChunkLayout<f64>,
    grid: RuntimeShape3,
}

/// Grid samples in one brick, `5³`.
const SAMPLES_PER_BRICK: u64 = ((BRICK_CELLS + 1) as u64).pow(3);

impl World {
    fn new() -> Self {
        let brick = ChunkLayout::<f64>::new(BRICK_CELLS, CELL_SIZE, [ORIGIN; 3]).expect("bricks");
        let chunk = ChunkLayout::<f64>::new(CHUNK_CELLS, CELL_SIZE, [ORIGIN; 3]).expect("chunks");
        let grid = brick.sample_shape().expect("brick sample grid fits u32");
        let size = grid.size();
        assert_eq!(
            u64::from(size[0]) * u64::from(size[1]) * u64::from(size[2]),
            SAMPLES_PER_BRICK,
            "the brick sample grid is not (cells + 1)^3"
        );
        Self { brick, chunk, grid }
    }

    fn brick_cube(&self, id: ChunkId) -> Cube {
        Cube::new(
            self.brick.sample_origin(id),
            f64::from(BRICK_CELLS) * CELL_SIZE,
        )
    }

    fn chunk_cube(&self, id: ChunkId) -> Cube {
        Cube::new(
            self.chunk.sample_origin(id),
            f64::from(CHUNK_CELLS) * CELL_SIZE,
        )
    }
}

/// The chunk a brick belongs to.
fn parent_of(brick: ChunkId) -> ChunkId {
    ChunkId::new(brick.coords.map(|c| c.div_euclid(BRICKS_PER_CHUNK_AXIS)))
}

// ── the fixture ─────────────────────────────────────────────────────────────

/// The first `y` in the world where the base field changes sign at `(x, z)`.
///
/// Only the dig path uses this, and only down the world's own mid-`z` line,
/// which is `P-72`'s inherited construction and is asserted rather than assumed:
/// the panic is `P-72`'s own, and it is what stopped two of that harness's runs
/// being void.
fn surface_y(base: &Base, x: f64, z: f64) -> f64 {
    let steps = 2048i32;
    let mut prev = base.sample([x, ORIGIN, z]);
    for i in 1..=steps {
        let y = ORIGIN + EXTENT * (f64::from(i) / f64::from(steps));
        let v = base.sample([x, y, z]);
        if (prev < 0.0) != (v < 0.0) {
            return y;
        }
        prev = v;
    }
    panic!("no surface crossing along y at ({x}, {z}): the trace would dig in empty space");
}

/// Centres of the cells of a 64³ scan of the world whose corners disagree in
/// sign: the base field's surface, as a finite set.
///
/// **A per-column `y` probe cannot place a log brush and the first run said so.**
/// `gyroid`'s period is 2π against a world 4 units across, so whole columns of
/// `[-2, 2]³` never cross zero — the probe panicked at `(-0.940, -1.633)` on the
/// first attempt. Choosing the surface by *enumeration* has no failure case to
/// paper over, works for any base, and is what makes "on the surface" a
/// definition rather than a search.
fn surface_cells(base: &Base) -> Vec<[f64; 3]> {
    const N: usize = 64;
    let h = EXTENT / N as f64;
    let at = |i: usize| ORIGIN + i as f64 * h;
    let mut grid = vec![0.0f64; (N + 1) * (N + 1) * (N + 1)];
    for k in 0..=N {
        for j in 0..=N {
            for i in 0..=N {
                grid[(k * (N + 1) + j) * (N + 1) + i] = base.sample([at(i), at(j), at(k)]);
            }
        }
    }
    let idx = |i: usize, j: usize, k: usize| (k * (N + 1) + j) * (N + 1) + i;
    let mut out = Vec::new();
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                let first = grid[idx(i, j, k)] < 0.0;
                let mixed = [
                    idx(i + 1, j, k),
                    idx(i, j + 1, k),
                    idx(i + 1, j + 1, k),
                    idx(i, j, k + 1),
                    idx(i + 1, j, k + 1),
                    idx(i, j + 1, k + 1),
                    idx(i + 1, j + 1, k + 1),
                ]
                .into_iter()
                .any(|n| (grid[n] < 0.0) != first);
                if mixed {
                    out.push([
                        at(i) + h * 0.5,
                        at(j) + h * 0.5,
                        at(k) + h * 0.5,
                    ]);
                }
            }
        }
    }
    assert!(
        !out.is_empty(),
        "the base field has no surface inside the world, so no log brush could be placed \
         on it and every survivor count would be about empty space"
    );
    out
}

/// The eleven dig centres: `P-72`'s path, probed per edit at that edit's own
/// `x`, which is both what a player does and what stopped two of `P-72`'s runs
/// being void.
fn dig_centres(base: &Base) -> Vec<[f64; 3]> {
    let mid = ORIGIN + EXTENT * 0.5;
    (0..EDITS)
        .map(|i| {
            let t = (i as f64 + 0.5) / EDITS as f64;
            let x = ORIGIN + EXTENT * t;
            [x, surface_y(base, x, mid), mid]
        })
        .collect()
}

/// `LOG_MAX` surface-following brushes: the edit log a player leaves behind.
///
/// **Placement is on the surface, not uniform over the volume, and that is the
/// harder fixture.** `P-39` scattered its tape through the whole world, where
/// most brushes sit in open air or deep rock and distance rejects them for
/// free — `✗41` found 1,218 of 1,434 unnecessary survivors were more than a
/// cell clear of zero. A surface-following log puts every brush in the band
/// where the bricks that need re-meshing live, so only *lateral* distance can
/// reject it, which is exactly the quantity a smaller box resolves better and
/// therefore the quantity this experiment is about.
fn log_brushes(cells: &[[f64; 3]], seed: u64) -> Vec<Brush<Shape>> {
    let mut rng = Lcg::new(seed);
    let mut out = Vec::with_capacity(LOG_MAX);
    for _ in 0..LOG_MAX {
        let cell = cells[rng.next_u32() as usize % cells.len()];
        // Jitter within the scan cell, so the log is not on a lattice.
        let jitter = EXTENT / 64.0 * 0.5;
        let centre = [0, 1, 2].map(|a| cell[a] + rng.range(-jitter, jitter));
        let shape = if rng.next_u32() & 1 == 0 {
            Shape::Sphere(Sphere {
                center: centre,
                radius: rng.range(3.0, 8.0) * CELL_SIZE,
            })
        } else {
            let dir = [
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            ];
            let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
            let unit = if len > 1e-9 {
                [dir[0] / len, dir[1] / len, dir[2] / len]
            } else {
                [1.0, 0.0, 0.0]
            };
            let half = rng.range(2.0, 6.0) * CELL_SIZE;
            Shape::Capsule(Capsule {
                a: [
                    centre[0] - unit[0] * half,
                    centre[1] - unit[1] * half,
                    centre[2] - unit[2] * half,
                ],
                b: [
                    centre[0] + unit[0] * half,
                    centre[1] + unit[1] * half,
                    centre[2] + unit[2] * half,
                ],
                radius: rng.range(3.0, 7.0) * CELL_SIZE,
            })
        };
        let op = if rng.next_u32() & 1 == 0 {
            BrushOp::Add
        } else {
            BrushOp::Subtract
        };
        out.push(Brush { shape, op });
    }
    out
}

/// One bucket's tape: the first `total − EDITS` log brushes, then the dig.
fn tape_for(log: &[Brush<Shape>], dig: &[[f64; 3]], total: usize) -> Vec<Brush<Shape>> {
    let pre = total - EDITS;
    let mut tape: Vec<Brush<Shape>> = log[..pre].to_vec();
    for c in dig {
        tape.push(Brush {
            shape: Shape::Sphere(Sphere {
                center: *c,
                radius: DIG_CELLS * CELL_SIZE,
            }),
            op: BrushOp::Subtract,
        });
    }
    tape
}

// ── the trace ───────────────────────────────────────────────────────────────

/// One edit's re-mesh set, grouped by the chunk that owns each brick.
///
/// Grouped once, outside every timer: the chunk-level pass is *per chunk* and a
/// timed arm that re-derived the grouping would be timing a `BTreeMap`.
struct Step {
    tape_len: usize,
    chunks: Vec<(ChunkId, Vec<ChunkId>)>,
    bricks: usize,
}

fn plan_trace(world: &World, base: &Base, tape: &[Brush<Shape>], pre: usize) -> Vec<Step> {
    let mut dirty = DirtySet::new();
    let mut steps = Vec::with_capacity(EDITS);
    for step in 0..EDITS {
        let before = BrushStack {
            base: *base,
            brushes: &tape[..pre + step],
        };
        let after = BrushStack {
            base: *base,
            brushes: &tape[..=pre + step],
        };
        let c = tape[pre + step].shape;
        let radius = match c {
            Shape::Sphere(s) => s.radius,
            Shape::Capsule(cap) => cap.radius,
        };
        let centre = match c {
            Shape::Sphere(s) => s.center,
            Shape::Capsule(cap) => cap.a,
        };
        let lo_world = [0, 1, 2].map(|a| centre[a] - radius);
        let hi_world = [0, 1, 2].map(|a| centre[a] + radius);
        let lo = world.brick.cell_of(lo_world).map(|v| v - 1);
        let hi = world.brick.cell_of(hi_world).map(|v| v + 1);

        mark_edit(&world.brick, &before, &after, lo, hi, &mut dirty).expect("mark_edit");
        let mut grouped: BTreeMap<[i32; 3], Vec<ChunkId>> = BTreeMap::new();
        let mut bricks = 0usize;
        for id in dirty.iter() {
            grouped.entry(parent_of(id).coords).or_default().push(id);
            bricks += 1;
        }
        dirty.clear();
        // M-44: a step that marked nothing would contribute a fast time for
        // doing nothing, and every arm would agree on it.
        assert!(
            bricks > 0,
            "VOID: edit {step} marked no dirty brick, so the trace measures nothing"
        );
        steps.push(Step {
            tape_len: pre + step + 1,
            chunks: grouped
                .into_iter()
                .map(|(k, v)| (ChunkId::new(k), v))
                .collect(),
            bricks,
        });
    }
    steps
}

// ── the four timed arms ─────────────────────────────────────────────────────

struct Rig {
    mc: MarchingCubes<f64>,
    out: MeshBuffer<f64>,
    chunk_survivors: Vec<Brush<Shape>>,
    brick_survivors: Vec<Brush<Shape>>,
}

impl Rig {
    fn mesh(&mut self, world: &World, base: &Base, brushes: &[Brush<Shape>], id: ChunkId) {
        let field = BrushStack {
            base: *base,
            brushes,
        };
        self.out.reset();
        self.mc
            .extract(
                &field,
                &world.grid,
                world.brick.sample_origin(id),
                CELL_SIZE,
                &mut self.out,
            )
            .expect("brick extraction");
        black_box(&self.out);
    }
}

/// Every dirty brick re-meshed from the whole tape. The denominator.
fn arm_full(rig: &mut Rig, world: &World, base: &Base, tape: &[Brush<Shape>], steps: &[Step]) -> u128 {
    let t = Instant::now();
    for step in steps {
        let tk = &tape[..step.tape_len];
        for (_, bricks) in &step.chunks {
            for id in bricks {
                rig.mesh(world, base, tk, *id);
            }
        }
    }
    t.elapsed().as_nanos()
}

/// The floor: the same bricks re-meshed from the base field with an empty tape.
/// Not a competing design — the part of remesh that pruning can never remove,
/// and therefore the denominator of the SHARE.
fn arm_floor(rig: &mut Rig, world: &World, base: &Base, steps: &[Step]) -> u128 {
    let t = Instant::now();
    for step in steps {
        for (_, bricks) in &step.chunks {
            for id in bricks {
                rig.mesh(world, base, &[], *id);
            }
        }
    }
    t.elapsed().as_nanos()
}

/// `P-39` transplanted: one prune per chunk, every brick in it meshed from the
/// chunk's survivors. The pruning pass is **inside** the timer, because whether
/// the bookkeeping is affordable is the question.
fn arm_chunk(
    rig: &mut Rig,
    world: &World,
    base: &Base,
    tape: &[Brush<Shape>],
    steps: &[Step],
) -> u128 {
    let t = Instant::now();
    for step in steps {
        let tk = &tape[..step.tape_len];
        for (cid, bricks) in &step.chunks {
            prune_into(tk, base, world.chunk_cube(*cid), &mut rig.chunk_survivors);
            for id in bricks {
                let survivors = std::mem::take(&mut rig.chunk_survivors);
                rig.mesh(world, base, &survivors, *id);
                rig.chunk_survivors = survivors;
            }
        }
    }
    t.elapsed().as_nanos()
}

/// Both levels: prune per chunk, then prune the chunk's survivors per brick.
fn arm_both(
    rig: &mut Rig,
    world: &World,
    base: &Base,
    tape: &[Brush<Shape>],
    steps: &[Step],
) -> u128 {
    let t = Instant::now();
    for step in steps {
        let tk = &tape[..step.tape_len];
        for (cid, bricks) in &step.chunks {
            prune_into(tk, base, world.chunk_cube(*cid), &mut rig.chunk_survivors);
            for id in bricks {
                let chunk_list = std::mem::take(&mut rig.chunk_survivors);
                prune_into(
                    &chunk_list,
                    base,
                    world.brick_cube(*id),
                    &mut rig.brick_survivors,
                );
                let brick_list = std::mem::take(&mut rig.brick_survivors);
                rig.mesh(world, base, &brick_list, *id);
                rig.brick_survivors = brick_list;
                rig.chunk_survivors = chunk_list;
            }
        }
    }
    t.elapsed().as_nanos()
}

// ── soundness, counts and the survivor census ───────────────────────────────

/// Bit-for-bit equality of two meshes. `f64 == f64` calls `-0.0` equal to
/// `0.0`, and `-0.0` is exactly the value the selection lemma's one soft spot
/// would produce, so this compares bits.
fn bitwise_identical(a: &MeshBuffer<f64>, b: &MeshBuffer<f64>) -> bool {
    fn same(x: &[f64; 3], y: &[f64; 3]) -> bool {
        x[0].to_bits() == y[0].to_bits()
            && x[1].to_bits() == y[1].to_bits()
            && x[2].to_bits() == y[2].to_bits()
    }
    a.indices == b.indices
        && a.positions.len() == b.positions.len()
        && a.normals.len() == b.normals.len()
        && a.positions
            .iter()
            .zip(&b.positions)
            .all(|(x, y)| same(x, y))
        && a.normals.iter().zip(&b.normals).all(|(x, y)| same(x, y))
}

/// Everything the untimed census pass produces.
#[derive(Default)]
struct Census {
    remeshes: u64,
    distinct_bricks: u64,
    survivors_chunk: u64,
    survivors_brick: u64,
    survivors_brick_direct: u64,
    removed_chunk: u64,
    removed_brick_direct: u64,
    removed_both: u64,
    differing: u64,
    distinct_differing: u64,
    chunk_passes: u64,
    chunk_pass_brushes: u64,
    samples: u64,
    shape_evals_full: u64,
    shape_evals_chunk: u64,
    shape_evals_both: u64,
    meshes_identical: bool,
    hashes_identical: bool,
}

/// Re-mesh every dirty brick three ways, assert they agree bit for bit, and
/// count what each way would have evaluated.
fn census(world: &World, base: &Base, tape: &[Brush<Shape>], steps: &[Step]) -> Census {
    let mut c = Census {
        meshes_identical: true,
        hashes_identical: true,
        ..Census::default()
    };
    let mut mc = MarchingCubes::<f64>::new();
    let mut full = MeshBuffer::<f64>::new();
    let mut viac = MeshBuffer::<f64>::new();
    let mut viab = MeshBuffer::<f64>::new();
    let mut chunk_list: Vec<Brush<Shape>> = Vec::with_capacity(tape.len());
    let mut brick_list: Vec<Brush<Shape>> = Vec::with_capacity(tape.len());
    let mut direct_list: Vec<Brush<Shape>> = Vec::with_capacity(tape.len());
    let mut seen: BTreeSet<[i32; 3]> = BTreeSet::new();
    let mut seen_differing: BTreeSet<[i32; 3]> = BTreeSet::new();

    for step in steps {
        let tk = &tape[..step.tape_len];
        let l = step.tape_len as u64;
        for (cid, bricks) in &step.chunks {
            let s_chunk = prune_into(tk, base, world.chunk_cube(*cid), &mut chunk_list) as u64;
            c.chunk_passes += 1;
            c.chunk_pass_brushes += l;
            for id in bricks {
                let s_both =
                    prune_into(&chunk_list, base, world.brick_cube(*id), &mut brick_list) as u64;
                let s_direct =
                    prune_into(tk, base, world.brick_cube(*id), &mut direct_list) as u64;

                let calls = Cell::new(0u64);
                let stack_full = BrushStack {
                    base: *base,
                    brushes: tk,
                };
                let counted = Counted {
                    field: &stack_full,
                    calls: &calls,
                };
                full.reset();
                mc.extract(
                    &counted,
                    &world.grid,
                    world.brick.sample_origin(*id),
                    CELL_SIZE,
                    &mut full,
                )
                .expect("brick extraction");
                let samples_full = calls.get();

                calls.set(0);
                let stack_both = BrushStack {
                    base: *base,
                    brushes: &brick_list,
                };
                let counted = Counted {
                    field: &stack_both,
                    calls: &calls,
                };
                viab.reset();
                mc.extract(
                    &counted,
                    &world.grid,
                    world.brick.sample_origin(*id),
                    CELL_SIZE,
                    &mut viab,
                )
                .expect("brick extraction");
                assert_eq!(
                    samples_full,
                    calls.get(),
                    "the pruned and unpruned arms made different numbers of Sdf::sample \
                     calls on brick {:?}, so they are not meshing the same cells",
                    id.coords
                );

                let stack_chunk = BrushStack {
                    base: *base,
                    brushes: &chunk_list,
                };
                viac.reset();
                mc.extract(
                    &stack_chunk,
                    &world.grid,
                    world.brick.sample_origin(*id),
                    CELL_SIZE,
                    &mut viac,
                )
                .expect("brick extraction");

                let ok = bitwise_identical(&full, &viab) && bitwise_identical(&full, &viac);
                let hashed =
                    mesh_hash(&full) == mesh_hash(&viab) && mesh_hash(&full) == mesh_hash(&viac);
                c.meshes_identical &= ok;
                c.hashes_identical &= hashed;

                c.remeshes += 1;
                if seen.insert(id.coords) {
                    c.distinct_bricks += 1;
                }
                c.survivors_chunk += s_chunk;
                c.survivors_brick += s_both;
                c.survivors_brick_direct += s_direct;
                c.removed_chunk += l - s_chunk;
                c.removed_brick_direct += l - s_direct;
                c.removed_both += l - s_both;
                // The brick list is pruned *from* the chunk list, so it is a
                // subsequence of it and a length difference is exactly a set
                // difference. No element comparison is needed and none would be
                // sound: `Brush<Shape>` carries floats.
                if s_both != s_chunk {
                    c.differing += 1;
                    if seen_differing.insert(id.coords) {
                        c.distinct_differing += 1;
                    }
                }
                c.samples += samples_full;
                c.shape_evals_full += samples_full * l;
                c.shape_evals_chunk += samples_full * s_chunk;
                c.shape_evals_both += samples_full * s_both + s_chunk;
            }
        }
    }
    // Every chunk pass evaluates one shape per brush in the tape it is given.
    c.shape_evals_chunk += c.chunk_pass_brushes;
    c.shape_evals_both += c.chunk_pass_brushes;
    c
}

/// The whole world's per-brick lists at one tape length: what "carrying" costs.
struct WorldLists {
    entries: u64,
    occupied: u64,
    max_chunk_survivors: u64,
}

fn world_lists(world: &World, base: &Base, tape: &[Brush<Shape>]) -> WorldLists {
    let mut chunk_list: Vec<Brush<Shape>> = Vec::with_capacity(tape.len());
    let mut brick_list: Vec<Brush<Shape>> = Vec::with_capacity(tape.len());
    let mut out = WorldLists {
        entries: 0,
        occupied: 0,
        max_chunk_survivors: 0,
    };
    for cz in 0..CHUNKS_PER_AXIS as i32 {
        for cy in 0..CHUNKS_PER_AXIS as i32 {
            for cx in 0..CHUNKS_PER_AXIS as i32 {
                let cid = ChunkId::new([cx, cy, cz]);
                let s = prune_into(tape, base, world.chunk_cube(cid), &mut chunk_list) as u64;
                out.max_chunk_survivors = out.max_chunk_survivors.max(s);
                for bz in 0..BRICKS_PER_CHUNK_AXIS {
                    for by in 0..BRICKS_PER_CHUNK_AXIS {
                        for bx in 0..BRICKS_PER_CHUNK_AXIS {
                            let id = ChunkId::new([
                                cx * BRICKS_PER_CHUNK_AXIS + bx,
                                cy * BRICKS_PER_CHUNK_AXIS + by,
                                cz * BRICKS_PER_CHUNK_AXIS + bz,
                            ]);
                            let n =
                                prune_into(&chunk_list, base, world.brick_cube(id), &mut brick_list)
                                    as u64;
                            out.entries += n;
                            if n > 0 {
                                out.occupied += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    out
}

/// **The representation C2 is costed against, named rather than assumed.**
///
/// A CSR arena over the whole world's bricks:
///
/// - `offsets: [u32; bricks + 1]` — **4 bytes per brick**, and it is dense
///   because the structure is indexed by brick id in O(1) during re-mesh. This
///   term is paid whether or not a brick has any survivor.
/// - `indices: [u16; entries]` — **2 bytes per surviving brush**, an index into
///   the world edit log. `u16` caps the log at 65,535 edits; Dreams goes to
///   100,000, so `bytes_per_brick_u32idx` costs the version that does not.
///
/// The alternatives are costed beside it because they are what a reader would
/// otherwise assume: a per-brick `u64` bitmask over the parent chunk's survivor
/// list is the natural hierarchical encoding and is **exactly 8 bytes**, and a
/// two-level sparse form (a 512-bit occupancy mask per chunk, offsets only for
/// occupied bricks) is what wins when occupancy is low.
fn bytes_csr_u16(entries: u64) -> f64 {
    (4.0 * (BRICKS_TOTAL as f64 + 1.0) + 2.0 * entries as f64) / BRICKS_TOTAL as f64
}

fn bytes_csr_u32(entries: u64) -> f64 {
    (4.0 * (BRICKS_TOTAL as f64 + 1.0) + 4.0 * entries as f64) / BRICKS_TOTAL as f64
}

fn bytes_sparse(entries: u64, occupied: u64) -> f64 {
    let chunks = f64::from(CHUNKS_PER_AXIS).powi(3);
    // 512 bricks per chunk = 64 bytes of occupancy mask, then a u32 offset for
    // each occupied brick and a u16 index for each entry.
    (64.0 * chunks + 4.0 * occupied as f64 + 2.0 * entries as f64) / BRICKS_TOTAL as f64
}

/// The tightest encoding that still indexes a brick in O(1): a `u32` arena base
/// per chunk plus a **chunk-relative** `u16` offset per brick.
///
/// A chunk holds 512 bricks and the tape is at most 60 brushes, so a chunk's
/// arena cannot exceed 30,720 entries and the offset fits in 16 bits — half the
/// index cost of the global CSR. This is the two-level form Dreams' block
/// hierarchy would actually produce, and it is here because **C2's verdict is
/// representation-sensitive and the sensitivity is the finding**, not because
/// the registered column changed.
fn bytes_chunk_arena(entries: u64) -> f64 {
    let chunks = f64::from(CHUNKS_PER_AXIS).powi(3);
    // 513 u16 offsets per chunk (one past the end) and one u32 arena base.
    (chunks * (2.0 * 513.0 + 4.0) + 2.0 * entries as f64) / BRICKS_TOTAL as f64
}

/// The payload alone, index charged at zero: the floor under every encoding
/// that stores a `u16` brush index per surviving brush. If this exceeds 8 the
/// clause is unreachable for any representation at all.
fn bytes_payload(entries: u64) -> f64 {
    2.0 * entries as f64 / BRICKS_TOTAL as f64
}

/// `✗41`'s necessity test, run over **every brick in the world** rather than
/// over the dirty ones.
///
/// C2's registration invokes *"`✗41`'s finding that 1,507 survivors cut to 73
/// necessary"*, and that reduction is the difference between what a Lipschitz
/// bound can produce and what the mesh actually needs. `✗41`'s own Part-5 rule
/// forbids chaining the two: its headline once multiplied a per-chunk median by
/// a world-wide total and the entry had to be rewritten. So the necessary count
/// is taken over **the same 32,768 bricks** the survivor count is taken over,
/// and the two byte figures are directly comparable with nothing multiplied.
///
/// A brush is necessary when re-meshing without it changes `mesh_hash`
/// (leave-one-out). `joint_failures` then re-meshes from the necessary list
/// **alone**, every individually-unnecessary survivor removed at once, and
/// counts the bricks where that changed the mesh — the upgrade `✗41` had to
/// measure rather than infer, and which held on 64 of 64 chunks there.
///
/// **It is a count and not an assertion, because it does not have to hold.**
/// Two brushes can each be individually redundant while the other is present
/// and jointly decisive; `✗41` said so in as many words and then measured that
/// its own fixture never did it. Whether this fixture does is the question, and
/// asserting the answer would have hidden it.
struct Necessity {
    entries: u64,
    bricks_with_triangles: u64,
    bricks_reduced: u64,
    joint_failures: u64,
}

fn world_necessity(world: &World, base: &Base, tape: &[Brush<Shape>]) -> Necessity {
    let mut mc = MarchingCubes::<f64>::new();
    let mut reference = MeshBuffer::<f64>::new();
    let mut trial = MeshBuffer::<f64>::new();
    let mut chunk_list: Vec<Brush<Shape>> = Vec::with_capacity(tape.len());
    let mut brick_list: Vec<Brush<Shape>> = Vec::with_capacity(tape.len());
    let mut minus: Vec<Brush<Shape>> = Vec::with_capacity(tape.len());
    let mut needed: Vec<Brush<Shape>> = Vec::with_capacity(tape.len());
    let mut out = Necessity {
        entries: 0,
        bricks_with_triangles: 0,
        bricks_reduced: 0,
        joint_failures: 0,
    };
    for cz in 0..CHUNKS_PER_AXIS as i32 {
        for cy in 0..CHUNKS_PER_AXIS as i32 {
            for cx in 0..CHUNKS_PER_AXIS as i32 {
                let cid = ChunkId::new([cx, cy, cz]);
                prune_into(tape, base, world.chunk_cube(cid), &mut chunk_list);
                for bz in 0..BRICKS_PER_CHUNK_AXIS {
                    for by in 0..BRICKS_PER_CHUNK_AXIS {
                        for bx in 0..BRICKS_PER_CHUNK_AXIS {
                            let id = ChunkId::new([
                                cx * BRICKS_PER_CHUNK_AXIS + bx,
                                cy * BRICKS_PER_CHUNK_AXIS + by,
                                cz * BRICKS_PER_CHUNK_AXIS + bz,
                            ]);
                            prune_into(
                                &chunk_list,
                                base,
                                world.brick_cube(id),
                                &mut brick_list,
                            );
                            let origin = world.brick.sample_origin(id);
                            mesh_into(&mut mc, &mut reference, base, &brick_list, world, origin);
                            if reference.triangle_count() > 0 {
                                out.bricks_with_triangles += 1;
                            }
                            let want = mesh_hash(&reference);
                            needed.clear();
                            for drop in 0..brick_list.len() {
                                minus.clear();
                                minus.extend(
                                    brick_list
                                        .iter()
                                        .enumerate()
                                        .filter(|(i, _)| *i != drop)
                                        .map(|(_, b)| *b),
                                );
                                mesh_into(&mut mc, &mut trial, base, &minus, world, origin);
                                if mesh_hash(&trial) != want {
                                    needed.push(brick_list[drop]);
                                }
                            }
                            out.entries += needed.len() as u64;
                            if needed.len() != brick_list.len() {
                                out.bricks_reduced += 1;
                                mesh_into(&mut mc, &mut trial, base, &needed, world, origin);
                                if !bitwise_identical(&reference, &trial) {
                                    out.joint_failures += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    out
}

fn mesh_into(
    mc: &mut MarchingCubes<f64>,
    out: &mut MeshBuffer<f64>,
    base: &Base,
    brushes: &[Brush<Shape>],
    world: &World,
    origin: [f64; 3],
) {
    let field = BrushStack {
        base: *base,
        brushes,
    };
    out.reset();
    mc.extract(&field, &world.grid, origin, CELL_SIZE, out)
        .expect("brick extraction");
}

// ── one row ─────────────────────────────────────────────────────────────────

fn median_f64(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

struct Case {
    field: &'static str,
    brushes: usize,
    bucket: &'static str,
    lipschitz: f64,
    steps: Vec<Step>,
    census: Census,
    lists_final: WorldLists,
    necessity: Necessity,
    bytes_mean_trace: f64,
    ms_full: f64,
    ms_chunk: f64,
    ms_both: f64,
    ms_floor: f64,
    speedup_both_min: f64,
    speedup_both_max: f64,
}

fn run_case(
    world: &World,
    base: &Base,
    field: &'static str,
    log: &[Brush<Shape>],
    dig: &[[f64; 3]],
    total: usize,
    bucket: &'static str,
) -> Case {
    let tape = tape_for(log, dig, total);
    assert_hard_ops(&tape);
    let pre = total - EDITS;
    let steps = plan_trace(world, base, &tape, pre);
    let census = census(world, base, &tape, &steps);

    // "Carrying" costed over the whole world at every step of the trace, then
    // both readings of "amortised" reported: the peak (final tape) sizes the
    // allocation, the mean is the literal reading of the clause.
    let mut bytes_steps = Vec::with_capacity(EDITS);
    let mut lists_final = WorldLists {
        entries: 0,
        occupied: 0,
        max_chunk_survivors: 0,
    };
    for step in &steps {
        let lists = world_lists(world, base, &tape[..step.tape_len]);
        bytes_steps.push(bytes_csr_u16(lists.entries));
        lists_final = lists;
    }
    let bytes_mean_trace = bytes_steps.iter().sum::<f64>() / bytes_steps.len() as f64;

    // ✗41's necessity test over the same 32,768 bricks, at the same final tape,
    // so `bytes_per_brick` and `bytes_per_brick_necessary` are two statistics
    // over one population and neither has to be multiplied by the other.
    let necessity = world_necessity(world, base, &tape);

    let mut rig = Rig {
        mc: MarchingCubes::new(),
        out: MeshBuffer::new(),
        chunk_survivors: Vec::with_capacity(tape.len()),
        brick_survivors: Vec::with_capacity(tape.len()),
    };
    // One untimed warm-up of each arm, so no allocation and no cold branch
    // predictor is inside a timed region.
    arm_full(&mut rig, world, base, &tape, &steps);
    arm_chunk(&mut rig, world, base, &tape, &steps);
    arm_both(&mut rig, world, base, &tape, &steps);
    arm_floor(&mut rig, world, base, &steps);

    let mut full = Vec::with_capacity(REPS);
    let mut chunk = Vec::with_capacity(REPS);
    let mut both = Vec::with_capacity(REPS);
    let mut floor = Vec::with_capacity(REPS);
    let mut ratio_both = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        // M-281: the four arms are interleaved inside one repetition, so a
        // ratio is always taken across measurements a few milliseconds apart
        // rather than across two halves of a run.
        let f = arm_full(&mut rig, world, base, &tape, &steps) as f64;
        let c = arm_chunk(&mut rig, world, base, &tape, &steps) as f64;
        let b = arm_both(&mut rig, world, base, &tape, &steps) as f64;
        let z = arm_floor(&mut rig, world, base, &steps) as f64;
        ratio_both.push(f / b);
        full.push(f / 1e6);
        chunk.push(c / 1e6);
        both.push(b / 1e6);
        floor.push(z / 1e6);
    }
    ratio_both.sort_by(f64::total_cmp);

    Case {
        field,
        brushes: total,
        bucket,
        lipschitz: base
            .value_bound()
            .lipschitz()
            .expect("both bases declare a Lipschitz constant"),
        steps,
        census,
        lists_final,
        necessity,
        bytes_mean_trace,
        ms_full: median_f64(&mut full),
        ms_chunk: median_f64(&mut chunk),
        ms_both: median_f64(&mut both),
        ms_floor: median_f64(&mut floor),
        speedup_both_min: ratio_both[0],
        speedup_both_max: ratio_both[REPS - 1],
    }
}

type Row = Vec<(&'static str, String)>;

fn row_of(c: &Case) -> Row {
    let cs = &c.census;
    let remeshes = cs.remeshes;
    let speedup_chunk = c.ms_full / c.ms_chunk;
    let speedup_both = c.ms_full / c.ms_both;
    let share = (c.ms_full - c.ms_floor) / c.ms_full;
    let ceiling = c.ms_full / c.ms_floor;
    let over_chunk = cs.removed_both as i64 - cs.removed_chunk as i64;
    let over_direct = cs.removed_both as i64 - cs.removed_brick_direct as i64;
    let additional = over_chunk.min(over_direct);
    let bytes = bytes_csr_u16(c.lists_final.entries);
    let dirty_total: usize = c.steps.iter().map(|s| s.bricks).sum();

    vec![
        ("field", c.field.to_string()),
        ("log_bucket", c.bucket.to_string()),
        ("chunk_cells", CHUNK_CELLS.to_string()),
        ("brick_cells", BRICK_CELLS.to_string()),
        ("brushes", c.brushes.to_string()),
        ("bricks", BRICKS_TOTAL.to_string()),
        ("edits", EDITS.to_string()),
        ("remeshes", remeshes.to_string()),
        ("dirty_bricks_total", dirty_total.to_string()),
        ("distinct_dirty_bricks", cs.distinct_bricks.to_string()),
        ("chunk_passes", cs.chunk_passes.to_string()),
        ("survivors_chunk_level", cs.survivors_chunk.to_string()),
        ("survivors_brick_level", cs.survivors_brick.to_string()),
        (
            "survivors_brick_direct",
            cs.survivors_brick_direct.to_string(),
        ),
        ("removed_chunk_only", cs.removed_chunk.to_string()),
        ("removed_brick_only", cs.removed_brick_direct.to_string()),
        ("removed_both", cs.removed_both.to_string()),
        ("additional_removed", additional.to_string()),
        ("additional_over_chunk_only", over_chunk.to_string()),
        ("additional_over_brick_only", over_direct.to_string()),
        (
            "bricks_differing_from_parent",
            cs.differing.to_string(),
        ),
        (
            "distinct_bricks_differing",
            cs.distinct_differing.to_string(),
        ),
        (
            "mean_survivors_chunk",
            format!("{:.4}", cs.survivors_chunk as f64 / remeshes as f64),
        ),
        (
            "mean_survivors_brick",
            format!("{:.4}", cs.survivors_brick as f64 / remeshes as f64),
        ),
        ("world_brick_entries", c.lists_final.entries.to_string()),
        ("world_bricks_occupied", c.lists_final.occupied.to_string()),
        (
            "max_chunk_survivors",
            c.lists_final.max_chunk_survivors.to_string(),
        ),
        (
            "mean_survivors_world_brick",
            format!("{:.6}", c.lists_final.entries as f64 / BRICKS_TOTAL as f64),
        ),
        ("bytes_per_brick", format!("{bytes:.6}")),
        (
            "bytes_per_brick_mean_trace",
            format!("{:.6}", c.bytes_mean_trace),
        ),
        (
            "bytes_per_brick_u32idx",
            format!("{:.6}", bytes_csr_u32(c.lists_final.entries)),
        ),
        (
            "bytes_per_brick_sparse",
            format!(
                "{:.6}",
                bytes_sparse(c.lists_final.entries, c.lists_final.occupied)
            ),
        ),
        ("bytes_per_brick_bitmask64", format!("{:.6}", 8.0)),
        (
            "bytes_per_brick_chunkarena",
            format!("{:.6}", bytes_chunk_arena(c.lists_final.entries)),
        ),
        (
            "bytes_per_brick_payload",
            format!("{:.6}", bytes_payload(c.lists_final.entries)),
        ),
        (
            "world_necessary_entries",
            c.necessity.entries.to_string(),
        ),
        (
            "bytes_per_brick_necessary",
            format!("{:.6}", bytes_csr_u16(c.necessity.entries)),
        ),
        (
            "survivors_to_necessary_world",
            format!(
                "{:.4}",
                c.lists_final.entries as f64 / c.necessity.entries.max(1) as f64
            ),
        ),
        (
            "world_bricks_with_triangles",
            c.necessity.bricks_with_triangles.to_string(),
        ),
        (
            "world_bricks_leave_one_out_reduced",
            c.necessity.bricks_reduced.to_string(),
        ),
        (
            "necessary_only_bricks_differing",
            c.necessity.joint_failures.to_string(),
        ),
        (
            "necessary_only_identical",
            (c.necessity.joint_failures == 0).to_string(),
        ),
        ("remesh_full_ms", format!("{:.4}", c.ms_full)),
        ("remesh_chunk_ms", format!("{:.4}", c.ms_chunk)),
        ("remesh_both_ms", format!("{:.4}", c.ms_both)),
        ("remesh_floor_ms", format!("{:.4}", c.ms_floor)),
        ("share_prunable", format!("{share:.6}")),
        ("ceiling_speedup", format!("{ceiling:.4}")),
        ("speedup_chunk_only", format!("{speedup_chunk:.4}")),
        ("speedup_both", format!("{speedup_both:.4}")),
        ("speedup_both_min", format!("{:.4}", c.speedup_both_min)),
        ("speedup_both_max", format!("{:.4}", c.speedup_both_max)),
        ("shape_evals_full", cs.shape_evals_full.to_string()),
        ("shape_evals_chunk", cs.shape_evals_chunk.to_string()),
        ("shape_evals_both", cs.shape_evals_both.to_string()),
        (
            "eval_ratio_chunk",
            format!(
                "{:.4}",
                cs.shape_evals_full as f64 / cs.shape_evals_chunk as f64
            ),
        ),
        (
            "eval_ratio_both",
            format!(
                "{:.4}",
                cs.shape_evals_full as f64 / cs.shape_evals_both as f64
            ),
        ),
        ("samples_total", cs.samples.to_string()),
        ("samples_per_brick", SAMPLES_PER_BRICK.to_string()),
        ("base_lipschitz", format!("{:.6}", c.lipschitz)),
        ("meshes_bit_identical", cs.meshes_identical.to_string()),
        ("hashes_identical", cs.hashes_identical.to_string()),
        ("reps", REPS.to_string()),
        ("clock", "std_time_Instant_monotonic".to_string()),
        ("c1_holds", (speedup_both >= 2.0).to_string()),
        ("c2_holds", (bytes < 8.0).to_string()),
        ("c3_holds", (additional > 0).to_string()),
    ]
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-90");

    // The availability finding, asserted rather than timed: with no Lipschitz
    // constant for the base there is no enclosure to test a brush against, so
    // the mechanism does not exist on `fbm_terrain` at any granularity.
    assert!(
        FbmTerrain::<f64>::canonical()
            .value_bound()
            .lipschitz()
            .is_none(),
        "fbm_terrain now declares a Lipschitz constant, so it can carry a pruning arm \
         and this row's field set is out of date"
    );

    let world = World::new();
    let bases: [(&'static str, Base, u64); 2] = [
        ("gyroid", Base::Gyroid(Gyroid::<f64>::canonical()), 0x90_5EED_6C11_0001),
        (
            "ground_slab",
            Base::Ground(BoxExact {
                center: [0.0, -3.0, 0.0],
                half_extents: [4.0, 3.0, 4.0],
            }),
            0x90_5EED_6C11_0002,
        ),
    ];
    let bucket_names = ["1-15", "16-30", "31-45", "46-60"];

    println!(
        "{:>12} {:>7} {:>8} {:>9} {:>10} {:>10} {:>9} {:>9} {:>9} {:>8} {:>7}",
        "field",
        "brushes",
        "remeshes",
        "surv/chunk",
        "surv/brick",
        "differing",
        "x chunk",
        "x both",
        "ceiling",
        "B/brick",
        "extra"
    );

    let mut cases: Vec<Case> = Vec::new();
    for (field, base, seed) in bases {
        let dig = dig_centres(&base);
        let cells = surface_cells(&base);
        let log = log_brushes(&cells, seed);
        for (i, total) in BUCKETS.into_iter().enumerate() {
            let case = run_case(&world, &base, field, &log, &dig, total, bucket_names[i]);
            let row = row_of(&case);
            let get = |k: &str| -> String {
                row.iter()
                    .find(|(n, _)| *n == k)
                    .map(|(_, v)| v.clone())
                    .expect("column present")
            };
            println!(
                "{:>12} {:>7} {:>8} {:>9} {:>10} {:>10} {:>9} {:>9} {:>9} {:>8} {:>7}",
                case.field,
                case.brushes,
                case.census.remeshes,
                get("mean_survivors_chunk"),
                get("mean_survivors_brick"),
                case.census.differing,
                get("speedup_chunk_only"),
                get("speedup_both"),
                get("ceiling_speedup"),
                get("bytes_per_brick"),
                get("additional_removed"),
            );
            cases.push(case);
        }
    }

    let rows: Vec<Row> = cases.iter().map(row_of).collect();
    let num = |row: &Row, k: &str| -> f64 {
        row.iter()
            .find(|(n, _)| *n == k)
            .map(|(_, v)| v.parse::<f64>().expect("numeric column"))
            .expect("column present")
    };

    // ── the registered vacuity control ──────────────────────────────────────
    let differing: u64 = cases.iter().map(|c| c.census.differing).sum();
    assert!(
        differing > 0,
        "VACUITY: no brick in any trace had a surviving-brush set different from its \
         parent chunk's, so C3 could not have fired and the zero means nothing"
    );
    assert!(
        cases.iter().all(|c| c.census.meshes_identical),
        "SOUNDNESS: a pruned arm produced a different mesh, so its speedup is a speedup \
         for the wrong answer"
    );
    assert!(
        cases.iter().all(|c| c.census.hashes_identical),
        "SOUNDNESS: mesh_hash disagreed between arms"
    );

    let c1 = rows
        .iter()
        .filter(|r| {
            r.iter()
                .any(|(n, v)| *n == "log_bucket" && v == "46-60")
        })
        .map(|r| (r.iter().find(|(n, _)| *n == "field").expect("field").1.clone(), num(r, "speedup_both"), num(r, "speedup_chunk_only"), num(r, "ceiling_speedup")))
        .collect::<Vec<_>>();
    println!("\nC1, the registered 46-60 bucket:");
    for (field, both, chunk, ceiling) in &c1 {
        println!(
            "  {field:>12}  chunk-only x{chunk:.4}   both x{both:.4}   ceiling x{ceiling:.4}  -> {}",
            if *both >= 2.0 { "HELD" } else { "FALSIFIED" }
        );
    }
    println!(
        "\nC2, bytes per brick over {BRICKS_TOTAL} bricks. REGISTERED COLUMN is the CSR: \
         u32 offsets over every brick + u16 index per survivor."
    );
    for r in &rows {
        println!(
            "  {:>12} {:>2} brushes  peak {:.4} B  mean-over-trace {:.4} B | payload-only \
             {:.4}  chunk-arena {:.4}  sparse {:.4}  bitmask64 8.0000 -> {}",
            r.iter().find(|(n, _)| *n == "field").expect("field").1,
            r.iter().find(|(n, _)| *n == "brushes").expect("brushes").1,
            num(r, "bytes_per_brick"),
            num(r, "bytes_per_brick_mean_trace"),
            num(r, "bytes_per_brick_payload"),
            num(r, "bytes_per_brick_chunkarena"),
            num(r, "bytes_per_brick_sparse"),
            if num(r, "bytes_per_brick") < 8.0 {
                "HELD"
            } else {
                "FALSIFIED"
            }
        );
    }
    println!(
        "\nC2's cited premise, x41's survivors-to-necessary, over the SAME 32,768 bricks \
         (never chained onto anything):"
    );
    for r in &rows {
        println!(
            "  {:>12} {:>2} brushes  survivors {:>7}  necessary {:>7}  x{:.4}  \
             bytes-if-necessary-only {:.4}  tri bricks {:>6}  reduced bricks {:>6}  \
             JOINT FAILURES {:>6}",
            r.iter().find(|(n, _)| *n == "field").expect("field").1,
            r.iter().find(|(n, _)| *n == "brushes").expect("brushes").1,
            r.iter()
                .find(|(n, _)| *n == "world_brick_entries")
                .expect("col")
                .1,
            r.iter()
                .find(|(n, _)| *n == "world_necessary_entries")
                .expect("col")
                .1,
            num(r, "survivors_to_necessary_world"),
            num(r, "bytes_per_brick_necessary"),
            r.iter()
                .find(|(n, _)| *n == "world_bricks_with_triangles")
                .expect("col")
                .1,
            r.iter()
                .find(|(n, _)| *n == "world_bricks_leave_one_out_reduced")
                .expect("col")
                .1,
            r.iter()
                .find(|(n, _)| *n == "necessary_only_bricks_differing")
                .expect("col")
                .1,
        );
    }
    println!("\nC3, composition against the better of the two single-level culls:");
    for r in &rows {
        println!(
            "  {:>12} {:>2} brushes  removed chunk-only {:>7}  brick-only {:>7}  both {:>7}  \
             additional {:>4}  differing bricks {:>5} -> {}",
            r.iter().find(|(n, _)| *n == "field").expect("field").1,
            r.iter().find(|(n, _)| *n == "brushes").expect("brushes").1,
            r.iter()
                .find(|(n, _)| *n == "removed_chunk_only")
                .expect("col")
                .1,
            r.iter()
                .find(|(n, _)| *n == "removed_brick_only")
                .expect("col")
                .1,
            r.iter().find(|(n, _)| *n == "removed_both").expect("col").1,
            r.iter()
                .find(|(n, _)| *n == "additional_removed")
                .expect("col")
                .1,
            r.iter()
                .find(|(n, _)| *n == "bricks_differing_from_parent")
                .expect("col")
                .1,
            if num(r, "additional_removed") > 0.0 {
                "HELD"
            } else {
                "FALSIFIED"
            }
        );
    }
    println!(
        "\nexact counts (machine-independent), 46-60 bucket: \
         eval_ratio_chunk / eval_ratio_both"
    );
    for r in rows.iter().filter(|r| {
        r.iter()
            .any(|(n, v)| *n == "log_bucket" && v == "46-60")
    }) {
        println!(
            "  {:>12}  x{:.4} / x{:.4}   shape evals {} -> {} -> {}",
            r.iter().find(|(n, _)| *n == "field").expect("field").1,
            num(r, "eval_ratio_chunk"),
            num(r, "eval_ratio_both"),
            r.iter()
                .find(|(n, _)| *n == "shape_evals_full")
                .expect("col")
                .1,
            r.iter()
                .find(|(n, _)| *n == "shape_evals_chunk")
                .expect("col")
                .1,
            r.iter()
                .find(|(n, _)| *n == "shape_evals_both")
                .expect("col")
                .1,
        );
    }

    common::experiment::run(prereg, |run| {
        for r in &rows {
            run.record(r);
        }
    });
}
