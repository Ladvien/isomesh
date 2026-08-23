//! **P-39 — Lipschitz tape pruning of the brush stack.**
//!
//! Ticket: R-038. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p39
//! ```
//!
//! Writes `docs/experiments/p-39.csv`.
//!
//! # What is being measured
//!
//! [`BrushStack::sample`](isomesh::brush::BrushStack) is a linear fold over every
//! brush, and `MarchingCubes::extract` prefills the whole sample grid before any
//! cell work, so a 33³ chunk evaluates the entire edit history 35,937 times. The
//! claim is that most of that history cannot possibly matter inside any one
//! chunk, and that the ones that cannot matter can be *deleted* rather than
//! merely skipped — bit-exactly.
//!
//! No crate change is needed for this, and none is made. `BrushStack` already
//! takes a `&[Brush<S>]`, so a pruned tape is a shorter slice: the whole
//! mechanism is a selection over the slice, built here.
//!
//! # The bound, one sample per brush per chunk
//!
//! A shape with declared Lipschitz constant `l` varies by at most `l·r` over a
//! box of circumradius `r`, so `f(centre) ± l·r` encloses it. `l` is not
//! hard-coded: it comes from
//! [`BoundedSdf::value_bound`](isomesh::fields::BoundedSdf) via
//! [`FieldBound::lipschitz`](isomesh::fields::FieldBound::lipschitz), which
//! answers `1` for every exact distance field — `Sphere` and `Capsule` both.
//!
//! Two details make the enclosure an actual enclosure rather than an
//! approximate one:
//!
//! - **The box is bigger than the sample grid.** Marching cubes' normals come
//!   from `Sdf::gradient`, and `BrushStack` does not override it, so the default
//!   central differences sample `DIFF_STEP · max(|coord|, 1)` *outside* the
//!   sampled extent. A brush pruned on a bound that stopped at the grid could
//!   move a normal. The box is inflated by exactly that margin.
//! - **The bound arithmetic is rounded.** The Lipschitz inequality is about the
//!   exact function; `sample(centre)` and `f(c) ± l·r` are both floating-point,
//!   and a bound rounded the wrong way is not a bound. Every enclosure is
//!   widened by [`PAD_ULPS`] ULP of the magnitudes involved, which also covers
//!   the few-ULP evaluation error of these closed forms at the sample points
//!   the comparison is really about.
//!
//! # The lemma, and why the inequality is strict
//!
//! `apply(Add, f, s) = min(f, s)` and `apply(Subtract, f, s) = max(f, −s)`. IEEE
//! `min`/`max` **select** an operand rather than computing a new value, and
//! negation is exact, so deleting a provably-losing `Add` or `Subtract` moves
//! the result by exactly zero ULP.
//!
//! "Provably losing" is tested with a **strict** inequality, and that is
//! load-bearing rather than fastidious. `f64::min` is documented to return
//! *either* operand when they compare equal, which is only observable for
//! `+0.0` against `-0.0` — and `-0.0` is reachable here. `s.lo > v.hi` gives
//! `s(p) > v(p)` strictly at every point, so there is no tie to resolve and the
//! selection is forced. A non-strict test would leave a signed-zero hole in a
//! bit-exactness claim, for a pruning gain of measure zero.
//!
//! # The asymmetry is registered, not discovered
//!
//! `smooth_min` is **not** bit-exactly prunable in the losing direction. At
//! `h == 1` it returns `b + (a − b)`, which is not bit-identical to `a`. So the
//! registered arms contain `Add` and `Subtract` only, and
//! [`Policy::Sound`] never prunes a `SmoothAdd` at all.
//!
//! [`Policy::PruneSmoothLosers`] exists to *demonstrate* that, as a clearly
//! separated extra: the same fixture with every `Add` turned into a
//! `SmoothAdd`, pruned by the same losing test widened by `k`, and the meshes
//! compared. Its columns are prefixed `smooth_` and none of them decides a
//! registered clause.

mod common;

use std::time::Instant;

use isomesh::brush::{Brush, BrushOp, BrushStack, Capsule};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::{BoundedSdf, FieldBound, Sphere};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf, Shape3};

/// Chunks along each axis. 4³ = 64 chunks, as registered.
const CHUNKS_PER_AXIS: i32 = 4;

/// Cells per chunk axis, so 33 samples per axis.
const CELLS_PER_CHUNK: u32 = 32;

/// World units per cell. The chunk is 4 units across and the world 16.
const CELL_SIZE: f64 = 0.125;

/// The world's minimum corner, so the world is `[-8, 8]³`.
const WORLD_ORIGIN: f64 = -8.0;

/// Brushes in the tape, as registered.
const BRUSHES: usize = 64;

/// Radius of the solid the brushes carve. Its surface crosses every chunk but
/// the eight corners, so the fixture has interior, surface and empty chunks
/// without being arranged to.
const BASE_RADIUS: f64 = 6.0;

/// Join width for the `smooth_` extra, in world units — two cells.
const SMOOTH_K: f64 = 0.25;

/// Timed repetitions per arm per chunk. The median is reported; a mean would be
/// dragged by whichever run collided with a scheduler tick.
const REPS: usize = 5;

/// ULP of slack added to every enclosure bound.
///
/// Covers rounding in the bound arithmetic itself and in the evaluation of the
/// shape, at both the centre and the sample points the comparison is about.
/// Sixteen is far above the few-ULP error of a sphere or capsule distance and
/// costs a pruning decision only exactly at the boundary.
const PAD_ULPS: f64 = 16.0;

// ── the fixture ─────────────────────────────────────────────────────────────

/// A brush shape in this experiment.
///
/// One enum so the whole stack is a single `&[Brush<Shape>]` slice, which is
/// what makes a pruned tape a shorter slice of the same type. Both variants are
/// exact distance fields, so both declare `l == 1` — and they declare it through
/// the crate's own [`BoundedSdf`] rather than by a constant written here.
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

/// A 64-bit LCG, so the 64 brushes are the same 64 brushes on every machine and
/// every run.
///
/// This is fixture construction, not output: nothing measured depends on the
/// generator being good, only on it being reproducible.
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

/// The 64-brush tape: `Add` and `Subtract` over spheres and capsules, scattered
/// across the whole world.
///
/// Uniform over the world rather than concentrated on the solid's surface. That
/// is the *harder* fixture for the mechanism — surface-concentrated edits would
/// leave whole chunks touched by nothing, which prunes trivially — and it is the
/// reading of "scattered over a 4×4×4 chunk world" that cannot be accused of
/// arranging the answer.
fn tape() -> Vec<Brush<Shape>> {
    // P-39's seed, and nothing about the result depends on the value.
    let mut rng = Lcg::new(0x39_5EED_C0DE_1234);
    let mut out = Vec::with_capacity(BRUSHES);
    for _ in 0..BRUSHES {
        let centre = [
            rng.range(-6.5, 6.5),
            rng.range(-6.5, 6.5),
            rng.range(-6.5, 6.5),
        ];
        let shape = if rng.next_u32() & 1 == 0 {
            Shape::Sphere(Sphere {
                center: centre,
                radius: rng.range(0.35, 1.1),
            })
        } else {
            let dir = [
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            ];
            let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
            // A zero direction would be a sphere, which is a fine capsule but
            // not the one this row is meant to contribute.
            let unit = if len > 1e-9 {
                [dir[0] / len, dir[1] / len, dir[2] / len]
            } else {
                [1.0, 0.0, 0.0]
            };
            let half = rng.range(0.25, 1.0);
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
                radius: rng.range(0.3, 0.8),
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

/// The same tape with every `Add` replaced by a `SmoothAdd` of width
/// [`SMOOTH_K`]. Feeds the `smooth_` extra only.
fn smooth_tape(hard: &[Brush<Shape>]) -> Vec<Brush<Shape>> {
    hard.iter()
        .map(|b| match b.op {
            BrushOp::Add => Brush::smooth_add(b.shape, SMOOTH_K),
            _ => *b,
        })
        .collect()
}

// ── the bound ───────────────────────────────────────────────────────────────

/// The enclosure of a scalar field over a chunk.
#[derive(Clone, Copy, Debug)]
struct Interval {
    lo: f64,
    hi: f64,
}

/// The box a chunk's field evaluations can touch.
#[derive(Clone, Copy, Debug)]
struct ChunkBox {
    centre: [f64; 3],
    /// Circumradius: half the diagonal of the sampled extent, plus the margin
    /// `Sdf::gradient`'s central differences reach outside it.
    radius: f64,
}

impl ChunkBox {
    /// The box for a chunk whose sample grid starts at `origin` and spans `span`
    /// on every axis.
    fn new(origin: [f64; 3], span: f64) -> Self {
        let centre = [
            origin[0] + span * 0.5,
            origin[1] + span * 0.5,
            origin[2] + span * 0.5,
        ];
        // `h = DIFF_STEP * max(|p|, 1)` at the furthest corner bounds the
        // differencing reach anywhere in the box.
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
fn enclose<S: BoundedSdf<Scalar = f64>>(field: &S, chunk: ChunkBox) -> Interval {
    let l = field
        .value_bound()
        .lipschitz()
        .expect("every field in this fixture declares a Lipschitz constant");
    let value = field.sample(chunk.centre);
    let reach = l * chunk.radius;
    let slack = reach + pad(value, reach);
    Interval {
        lo: value - slack,
        hi: value + slack,
    }
}

/// Which pruning rule to apply to a `SmoothAdd`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Policy {
    /// The registered rule. `Add` and `Subtract` prune in the losing direction;
    /// a `SmoothAdd` never prunes, because `b + (a − b)` is not `a`.
    Sound,
    /// The rule the registration says must **not** be used, run to show what it
    /// costs. Prunes a `SmoothAdd` whose enclosure is `k` clear of the running
    /// one, which is exactly the `h == 1` region.
    PruneSmoothLosers,
}

/// What one pruning pass found.
#[derive(Clone, Copy, Debug, Default)]
struct PruneStats {
    survivors: usize,
    /// `Add`s that provably *win* everywhere in the chunk, so the whole tape
    /// prefix and the base field are dead. Counted, not exploited:
    /// `BrushStack` has no way to say "start from this brush" without replacing
    /// the base, which is a crate change and out of scope here.
    dominant_adds: usize,
}

/// Select the brushes that can still change the fold anywhere in `chunk`.
///
/// Order is preserved, because `Add` and `Subtract` do not commute with each
/// other — `BrushOp::commutes_with` is the crate's own statement of that.
fn prune_into(
    tape: &[Brush<Shape>],
    base: &Sphere<f64>,
    chunk: ChunkBox,
    policy: Policy,
    out: &mut Vec<Brush<Shape>>,
) -> PruneStats {
    out.clear();
    let mut stats = PruneStats::default();
    let mut v = enclose(base, chunk);
    for brush in tape {
        let s = enclose(&brush.shape, chunk);
        match brush.op {
            BrushOp::Add => {
                // Strictly above the running value everywhere, so `min` is
                // forced to select the running value at every point.
                if s.lo > v.hi {
                    continue;
                }
                if s.hi < v.lo {
                    stats.dominant_adds += 1;
                }
                v = Interval {
                    lo: v.lo.min(s.lo),
                    hi: v.hi.min(s.hi),
                };
            }
            BrushOp::Subtract => {
                // `max(v, -s)`. Negation is exact, so negating the enclosure is
                // exact too.
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
                if policy == Policy::PruneSmoothLosers && s.lo > v.hi + k {
                    continue;
                }
                // `smin(a, b) = min(a, b) − k(1 − |a − b|/k)²/4`, so the floor
                // sags by at most `k/4` below the plain minimum.
                let lo = v.lo.min(s.lo) - 0.25 * k;
                v = Interval {
                    lo: lo - pad(lo, 0.25 * k),
                    hi: v.hi.min(s.hi),
                };
            }
        }
        out.push(*brush);
        stats.survivors += 1;
    }
    stats
}

// ── meshing ─────────────────────────────────────────────────────────────────

/// One chunk's sample grid.
struct Grid {
    shape: RuntimeShape3,
    origin: [f64; 3],
    cell: f64,
    samples: f64,
}

fn extract_once<F: Sdf<Scalar = f64>>(
    mc: &mut MarchingCubes<f64>,
    out: &mut MeshBuffer<f64>,
    field: &F,
    grid: &Grid,
) {
    out.reset();
    mc.extract(field, &grid.shape, grid.origin, grid.cell, out)
        .expect("chunk extraction");
}

/// Median nanoseconds per sample over [`REPS`] timed runs, after one untimed
/// warm-up. The buffer is reused, so no allocation is inside a timed region.
fn median_ns_per_sample<F: Sdf<Scalar = f64>>(
    mc: &mut MarchingCubes<f64>,
    out: &mut MeshBuffer<f64>,
    field: &F,
    grid: &Grid,
) -> f64 {
    extract_once(mc, out, field, grid);
    let mut runs = [0.0f64; REPS];
    for slot in &mut runs {
        let t = Instant::now();
        extract_once(mc, out, field, grid);
        *slot = t.elapsed().as_secs_f64() * 1e9 / grid.samples;
    }
    runs.sort_by(f64::total_cmp);
    runs[REPS / 2]
}

/// Median nanoseconds for one whole pruning pass over the tape.
fn median_bound_ns(
    tape: &[Brush<Shape>],
    base: &Sphere<f64>,
    chunk: ChunkBox,
    out: &mut Vec<Brush<Shape>>,
) -> f64 {
    prune_into(tape, base, chunk, Policy::Sound, out);
    let mut runs = [0.0f64; REPS];
    for slot in &mut runs {
        let t = Instant::now();
        let stats = prune_into(tape, base, chunk, Policy::Sound, out);
        *slot = t.elapsed().as_secs_f64() * 1e9;
        // The pass must not be optimised away, and its result is the thing the
        // next arm consumes anyway.
        assert_eq!(stats.survivors, out.len(), "pruning lost a survivor");
    }
    runs.sort_by(f64::total_cmp);
    runs[REPS / 2]
}

/// Bit-for-bit equality of two meshes.
///
/// `f64 == f64` calls `-0.0` equal to `0.0`, and `-0.0` is exactly the value the
/// selection lemma's one soft spot would produce, so this compares bits. The
/// hash is checked separately: a hash match is necessary, this is decisive.
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

// ── the sweep ───────────────────────────────────────────────────────────────

/// Everything one chunk contributed.
#[derive(Clone, Copy, Debug)]
struct ChunkResult {
    id: [i32; 3],
    survivors: usize,
    dominant_adds: usize,
    ns_full: f64,
    ns_pruned: f64,
    bound_ns: f64,
    /// Samples in this chunk's grid, so the bound's share of a mesh is derived
    /// rather than transcribed.
    samples: f64,
    hash_full: u64,
    hash_pruned: u64,
    identical: bool,
    hash_matched: bool,
    vertices: usize,
    triangles: usize,
    smooth_survivors: usize,
    smooth_identical: bool,
}

impl ChunkResult {
    fn survivor_fraction(&self) -> f64 {
        self.survivors as f64 / BRUSHES as f64
    }

    fn speedup(&self) -> f64 {
        self.ns_full / self.ns_pruned
    }

    /// The pruning pass as a fraction of one pruned meshing of this chunk.
    fn bound_cost_fraction(&self) -> f64 {
        self.bound_ns / (self.ns_pruned * self.samples)
    }

    fn label(&self) -> String {
        format!("{}-{}-{}", self.id[0], self.id[1], self.id[2])
    }
}

/// Buffers that outlive the sweep so nothing allocates inside a timed region.
struct Rig {
    mc: MarchingCubes<f64>,
    full: MeshBuffer<f64>,
    pruned: MeshBuffer<f64>,
    survivors: Vec<Brush<Shape>>,
}

/// The whole tape set, hard and smooth, plus the solid they carve.
struct Fixture {
    base: Sphere<f64>,
    hard: Vec<Brush<Shape>>,
    smooth: Vec<Brush<Shape>>,
}

fn measure_chunk(
    rig: &mut Rig,
    fixture: &Fixture,
    layout: &ChunkLayout<f64>,
    id: ChunkId,
) -> ChunkResult {
    let origin = layout.sample_origin(id);
    let span = f64::from(layout.cells()) * layout.cell_size();
    let chunk = ChunkBox::new(origin, span);
    let shape = layout.sample_shape().expect("chunk sample grid fits u32");
    let size = shape.size();
    let grid = Grid {
        samples: f64::from(size[0]) * f64::from(size[1]) * f64::from(size[2]),
        shape,
        origin,
        cell: layout.cell_size(),
    };

    let bound_ns = median_bound_ns(&fixture.hard, &fixture.base, chunk, &mut rig.survivors);
    let stats = prune_into(
        &fixture.hard,
        &fixture.base,
        chunk,
        Policy::Sound,
        &mut rig.survivors,
    );

    let full_field = BrushStack {
        base: fixture.base,
        brushes: &fixture.hard,
    };
    let ns_full = median_ns_per_sample(&mut rig.mc, &mut rig.full, &full_field, &grid);
    let pruned_field = BrushStack {
        base: fixture.base,
        brushes: &rig.survivors,
    };
    let ns_pruned = median_ns_per_sample(&mut rig.mc, &mut rig.pruned, &pruned_field, &grid);

    let hash_full = mesh_hash(&rig.full);
    let hash_pruned = mesh_hash(&rig.pruned);
    let identical = bitwise_identical(&rig.full, &rig.pruned);
    let vertices = rig.full.vertex_count();
    let triangles = rig.full.triangle_count();

    // The extra, and the reason it is an extra: the same losing test on a
    // smooth stack is not bit-exact, and the registration says so up front.
    let smooth_stats = prune_into(
        &fixture.smooth,
        &fixture.base,
        chunk,
        Policy::PruneSmoothLosers,
        &mut rig.survivors,
    );
    let smooth_full = BrushStack {
        base: fixture.base,
        brushes: &fixture.smooth,
    };
    extract_once(&mut rig.mc, &mut rig.full, &smooth_full, &grid);
    let smooth_pruned = BrushStack {
        base: fixture.base,
        brushes: &rig.survivors,
    };
    extract_once(&mut rig.mc, &mut rig.pruned, &smooth_pruned, &grid);
    let smooth_identical = bitwise_identical(&rig.full, &rig.pruned);

    ChunkResult {
        id: id.coords,
        survivors: stats.survivors,
        dominant_adds: stats.dominant_adds,
        ns_full,
        ns_pruned,
        bound_ns,
        samples: grid.samples,
        hash_full,
        hash_pruned,
        identical,
        hash_matched: hash_full == hash_pruned,
        vertices,
        triangles,
        smooth_survivors: smooth_stats.survivors,
        smooth_identical,
    }
}

/// Upper median of a sorted-in-place copy. `sorted[len / 2]` is the same
/// convention `experiment_p38` uses for its timing medians, and for an even
/// count it is the larger of the two middles — the conservative direction for a
/// "under 0.5" clause.
fn median(values: &mut [f64]) -> f64 {
    values.sort_by(f64::total_cmp);
    values[values.len() / 2]
}

/// The aggregates every row carries.
struct Summary {
    survivor_median: f64,
    survivor_max: f64,
    survivor_min: f64,
    survivor_mean: f64,
    speedup_median: f64,
    /// The worst chunk. The registered clause is about "a chunk", so the chunk
    /// where the mechanism buys least belongs on the artefact next to the
    /// median that decides the clause.
    speedup_min: f64,
    chunks_at_or_above_target: usize,
    /// Chunks where nothing pruned at all — the shape the falsification
    /// criterion names, counted rather than left to the max column.
    chunks_all_survive: usize,
    world_ns_full: f64,
    world_ns_pruned: f64,
    mesh_identical: bool,
    hashes_matched: bool,
    empty_chunks: usize,
    smooth_differing: usize,
    smooth_survivor_median: f64,
}

impl Summary {
    fn of(rows: &[ChunkResult]) -> Self {
        let mut fractions: Vec<f64> = rows.iter().map(ChunkResult::survivor_fraction).collect();
        let mut speedups: Vec<f64> = rows.iter().map(ChunkResult::speedup).collect();
        let mut smooth: Vec<f64> = rows
            .iter()
            .map(|r| r.smooth_survivors as f64 / BRUSHES as f64)
            .collect();
        let total_full: f64 = rows.iter().map(|r| r.ns_full).sum();
        let total_pruned: f64 = rows.iter().map(|r| r.ns_pruned).sum();
        let n = rows.len() as f64;
        Self {
            survivor_median: median(&mut fractions),
            survivor_max: fractions.iter().copied().fold(f64::MIN, f64::max),
            survivor_min: fractions.iter().copied().fold(f64::MAX, f64::min),
            survivor_mean: fractions.iter().sum::<f64>() / n,
            speedup_median: median(&mut speedups),
            speedup_min: speedups.iter().copied().fold(f64::MAX, f64::min),
            chunks_at_or_above_target: speedups.iter().filter(|s| **s >= 1.25).count(),
            chunks_all_survive: rows.iter().filter(|r| r.survivors == BRUSHES).count(),
            world_ns_full: total_full / n,
            world_ns_pruned: total_pruned / n,
            mesh_identical: rows.iter().all(|r| r.identical),
            hashes_matched: rows.iter().all(|r| r.hash_matched),
            empty_chunks: rows.iter().filter(|r| r.triangles == 0).count(),
            smooth_differing: rows.iter().filter(|r| !r.smooth_identical).count(),
            smooth_survivor_median: median(&mut smooth),
        }
    }

    fn world_speedup(&self) -> f64 {
        self.world_ns_full / self.world_ns_pruned
    }
}

fn row_of(r: &ChunkResult, s: &Summary) -> Vec<(&'static str, String)> {
    let held = [
        s.survivor_median < 0.5,
        s.speedup_median >= 1.25,
        s.mesh_identical,
    ];
    vec![
        // Registered.
        ("brushes", BRUSHES.to_string()),
        ("chunks", (CHUNKS_PER_AXIS.pow(3)).to_string()),
        (
            "survivor_fraction_median",
            format!("{:.4}", s.survivor_median),
        ),
        ("survivor_fraction_max", format!("{:.4}", s.survivor_max)),
        ("ns_per_sample_full", format!("{:.4}", r.ns_full)),
        ("ns_per_sample_pruned", format!("{:.4}", r.ns_pruned)),
        ("speedup", format!("{:.4}", r.speedup())),
        ("mesh_identical", s.mesh_identical.to_string()),
        // Extra: which chunk this row is.
        ("chunk", r.label()),
        ("tape_length", BRUSHES.to_string()),
        ("survivors", r.survivors.to_string()),
        ("survivor_fraction", format!("{:.4}", r.survivor_fraction())),
        ("dominant_adds", r.dominant_adds.to_string()),
        ("vertices", r.vertices.to_string()),
        ("triangles", r.triangles.to_string()),
        ("chunk_mesh_identical", r.identical.to_string()),
        ("chunk_hash_matched", r.hash_matched.to_string()),
        ("mesh_hash_full", r.hash_full.to_string()),
        ("mesh_hash_pruned", r.hash_pruned.to_string()),
        ("bound_ns", format!("{:.1}", r.bound_ns)),
        (
            "bound_cost_fraction",
            format!("{:.8}", r.bound_cost_fraction()),
        ),
        // Extra: aggregates, repeated so every row can be read alone.
        ("survivor_fraction_min", format!("{:.4}", s.survivor_min)),
        ("survivor_fraction_mean", format!("{:.4}", s.survivor_mean)),
        (
            "mean_survivors",
            format!("{:.2}", s.survivor_mean * BRUSHES as f64),
        ),
        ("speedup_median", format!("{:.4}", s.speedup_median)),
        ("speedup_min", format!("{:.4}", s.speedup_min)),
        (
            "chunks_speedup_ge_1p25",
            s.chunks_at_or_above_target.to_string(),
        ),
        ("chunks_all_survive", s.chunks_all_survive.to_string()),
        (
            "world_ns_per_sample_full",
            format!("{:.4}", s.world_ns_full),
        ),
        (
            "world_ns_per_sample_pruned",
            format!("{:.4}", s.world_ns_pruned),
        ),
        ("world_speedup", format!("{:.4}", s.world_speedup())),
        ("empty_chunks", s.empty_chunks.to_string()),
        ("all_hashes_matched", s.hashes_matched.to_string()),
        ("c1_survivor_median_under_half", held[0].to_string()),
        ("c2_speedup_at_least_1p25", held[1].to_string()),
        ("c3_byte_identical", held[2].to_string()),
        // Extra: the registered asymmetry, demonstrated.
        (
            "smooth_losing_prune_identical",
            (s.smooth_differing == 0).to_string(),
        ),
        ("smooth_chunks_differing", s.smooth_differing.to_string()),
        (
            "smooth_survivor_fraction_median",
            format!("{:.4}", s.smooth_survivor_median),
        ),
        ("smooth_k", format!("{SMOOTH_K}")),
    ]
}

fn main() {
    let prereg = isomesh::experiment!("P-39");

    let fixture = {
        let hard = tape();
        let smooth = smooth_tape(&hard);
        Fixture {
            base: Sphere {
                center: [0.0; 3],
                radius: BASE_RADIUS,
            },
            hard,
            smooth,
        }
    };
    let adds = fixture.hard.iter().filter(|b| b.op == BrushOp::Add).count();
    let layout = ChunkLayout::new(CELLS_PER_CHUNK, CELL_SIZE, [WORLD_ORIGIN; 3])
        .expect("chunk layout is well formed");

    let mut rig = Rig {
        mc: MarchingCubes::new(),
        full: MeshBuffer::new(),
        pruned: MeshBuffer::new(),
        survivors: Vec::with_capacity(BRUSHES),
    };

    let mut rows = Vec::with_capacity((CHUNKS_PER_AXIS.pow(3)) as usize);
    for z in 0..CHUNKS_PER_AXIS {
        for y in 0..CHUNKS_PER_AXIS {
            for x in 0..CHUNKS_PER_AXIS {
                let r = measure_chunk(&mut rig, &fixture, &layout, ChunkId::new([x, y, z]));
                println!(
                    "chunk {:>6}  survivors {:>2}/{BRUSHES} ({:.3})  \
                     {:7.2} → {:7.2} ns/sample  ×{:.3}  tris {:>5}  \
                     identical {}  smooth-identical {}",
                    r.label(),
                    r.survivors,
                    r.survivor_fraction(),
                    r.ns_full,
                    r.ns_pruned,
                    r.speedup(),
                    r.triangles,
                    r.identical,
                    r.smooth_identical,
                );
                rows.push(r);
            }
        }
    }

    let summary = Summary::of(&rows);
    println!(
        "\ntape: {BRUSHES} brushes ({adds} Add, {} Subtract) over a sphere of radius {BASE_RADIUS}",
        BRUSHES - adds
    );
    println!(
        "C1 survivor fraction: median {:.4}, max {:.4}, min {:.4}, mean {:.4} → {}",
        summary.survivor_median,
        summary.survivor_max,
        summary.survivor_min,
        summary.survivor_mean,
        if summary.survivor_median < 0.5 {
            "HELD"
        } else {
            "FALSIFIED"
        }
    );
    println!(
        "   {} of {} chunks prune nothing at all (survivor fraction 1.0)",
        summary.chunks_all_survive,
        rows.len()
    );
    println!(
        "C2 speedup: median per chunk ×{:.4}, world ×{:.4}, worst chunk ×{:.4}, \
         {}/{} chunks at or above ×1.25 → {}",
        summary.speedup_median,
        summary.world_speedup(),
        summary.speedup_min,
        summary.chunks_at_or_above_target,
        rows.len(),
        if summary.speedup_median >= 1.25 {
            "HELD"
        } else {
            "FALSIFIED"
        }
    );
    println!(
        "C3 bit-exactness: {}/{} chunks byte-identical, {}/{} hashes equal → {}",
        rows.iter().filter(|r| r.identical).count(),
        rows.len(),
        rows.iter().filter(|r| r.hash_matched).count(),
        rows.len(),
        if summary.mesh_identical {
            "HELD"
        } else {
            "FALSIFIED"
        }
    );
    if !summary.mesh_identical {
        for r in rows.iter().filter(|r| !r.identical) {
            println!(
                "  DIFFERS: chunk {} survivors {} hash {} vs {}",
                r.label(),
                r.survivors,
                r.hash_full,
                r.hash_pruned
            );
        }
    }
    println!(
        "extra (not registered): smooth losing-direction prune differs on {}/{} chunks \
         at k = {SMOOTH_K}; {} empty chunks; bound cost {:.3e} of one pruned mesh",
        summary.smooth_differing,
        rows.len(),
        summary.empty_chunks,
        rows[0].bound_cost_fraction()
    );

    common::experiment::run(prereg, |run| {
        for r in &rows {
            run.record(&row_of(r, &summary));
        }
    });
}
