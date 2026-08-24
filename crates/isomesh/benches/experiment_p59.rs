//! **P-59 — how many surviving brushes are actually necessary.**
//!
//! Ticket: R-057. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p59
//! ```
//!
//! Writes `docs/experiments/p-59.csv`.
//!
//! # What is transcribed, and from where
//!
//! The fixture is P-39's, copied brush for brush out of
//! `benches/experiment_p39.rs` rather than imported. Benches in this crate do
//! not `use` one another, and a shared module would let one experiment's
//! maintenance silently move another experiment's published numbers. Copied
//! verbatim: `BRUSHES = 64` and the LCG-seeded [`tape`], the [`Shape`] enum with
//! its `Sdf` and `BoundedSdf` impls, [`Interval`], [`ChunkBox`] with the
//! central-difference margin that makes its box bigger than the sample grid,
//! [`pad`], [`enclose`], [`Policy`], [`prune_into`], and the 4³-chunk layout at
//! 32 cells and 0.125 world units per cell over a base sphere of radius 6.
//!
//! Comparability is checked rather than assumed. M-341 reports a **median of 19
//! survivors of 64**, `0.2969`. [`assert_comparable`] recomputes that median
//! from this harness's own rows and panics if it moved, because a fixture that
//! drifted by one brush makes every number below incomparable with the
//! measurement this experiment exists to extend — and reporting an
//! incomparable number under M-341's name is worse than reporting nothing.
//!
//! P-39's `Policy::PruneSmoothLosers` is deliberately **not** transcribed. It is
//! that experiment's negative control for the smooth-min asymmetry, it prunes in
//! a direction P-39's own registration calls out as not bit-exact, and it has
//! nothing to say about how tight the *sound* bound is. [`Policy`] therefore
//! carries the single variant this experiment runs, and the `SmoothAdd` arm of
//! [`prune_into`] is the sound one: never prune, widen the running enclosure by
//! the `k/4` the smooth minimum can sag.
//!
//! # What this harness owns rather than the source
//!
//! The ablation. Nothing in `src/` and nothing in P-39 measures *necessity* —
//! M-341 measured how many brushes the bound **keeps** and that the kept tape
//! meshes byte-identically, which are both statements about the pruner and
//! neither a statement about the tape it produced. Leave-one-out settles it
//! directly: for each chunk, for each brush the bound kept, drop that one brush
//! from the survivor tape *keeping the order of the rest*, re-mesh, and compare
//! `mesh_hash` against the survivors-only reference. Order is preserved because
//! `Add` and `Subtract` do not commute, so a reordered tape would be a different
//! field and the comparison would measure the reorder instead of the removal.
//!
//! No crate change is needed for this and none is made: `BrushStack` takes a
//! `&[Brush<S>]`, so an ablated tape is just a shorter slice built here.
//!
//! # Which columns decide which clause
//!
//! - **C1, the soundness control, computed and printed before any other number
//!   on the row.** `control_hash_unchanged` is `mesh_hash(full 64-brush tape) ==
//!   mesh_hash(survivors only)` — removing all `non_survivors_removed`
//!   non-survivors at once. If it is `false` on a chunk, the interval bound is
//!   unsound there and every other number on that row is void; the harness
//!   prints a `VOID` banner naming the chunk and **still writes the row**,
//!   because a suppressed row hides the exact failure the control exists to
//!   surface. C1 is the conjunction of `control_hash_unchanged` over all rows.
//! - **C2** is the median over chunks of `necessary_fraction = necessary /
//!   survivors`, against the registered `0.75`. `necessary` counts survivors
//!   whose individual removal changes `mesh_hash`. The extra column
//!   `unnecessary` carries `survivors - necessary` so C3's denominator is
//!   visible in the file rather than inferred from a subtraction.
//! - **C3** is the aggregate `unnecessary_far_from_surface / unnecessary` over
//!   all chunks, against the registered 90%. A survivor is *far from the
//!   surface* when `enclose(&brush.shape, chunk)` — the same enclosure the
//!   pruner itself consumed, not a second bound computed differently — has
//!   `lo > cell_size` or `hi < -cell_size`. `far_by_lo` and `far_by_hi` split
//!   that disjunction so the direction is readable instead of asserted.
//!
//! `unnecessary_far_fraction` is per-row and is `0` on a row with no unnecessary
//! survivors; the aggregate that decides C3 is `c3_far_fraction`, which divides
//! the summed numerator by the summed denominator and so is not a mean of
//! per-row ratios weighted by nothing.
//!
//! A chunk with zero survivors would make `necessary / survivors` undefined.
//! Rather than write a placeholder, [`measure_chunk`] panics naming the chunk.
//! P-39's committed distribution has a minimum of one survivor, so reaching that
//! panic means the fixture changed — which [`assert_comparable`] would also
//! catch, and which is the only honest thing to do about a ratio that has no
//! value.
//!
//! # Two readings the registration invites and this file refuses to assume
//!
//! **Leave-one-out does not compose, and there is a column for that.** C2's
//! wording is precise — "can be dropped **individually** with the mesh
//! bit-identical" — and that is exactly what `necessary` measures. The
//! *framing* sentence, "the bound over-keeps", invites the stronger reading that
//! the unnecessary survivors could all go at once, and leave-one-out cannot
//! support it: two brushes can each be redundant *while the other is present*
//! and jointly decisive, which is the ordinary behaviour of a `min`/`max` chain.
//! So the joint claim is measured rather than inferred.
//! `necessary_only_hash_unchanged` re-meshes each chunk from the `necessary`
//! brushes alone — every individually-unnecessary survivor removed at once,
//! order preserved — and compares against the reference. It decides no
//! registered clause; it exists so the entry that cites C2 cannot quietly
//! upgrade "individually" to "jointly".
//!
//! **An empty chunk makes every survivor unnecessary for free.** A chunk whose
//! reference mesh has no triangles hashes the same after any removal that leaves
//! it empty, so its `necessary` is zero without the bound having over-kept
//! anything. `triangles` is on every row so that confound can be divided out
//! instead of argued about.
//!
//! # Counted, not timed
//!
//! `remeshes` is an integer: the leave-one-outs plus the one control, so
//! `survivors + 1`. `ns_per_remesh` sits beside it as wall time per re-mesh and
//! **gates nothing** — M-348 is the incident where a discovery was demoted for
//! resting on a wall clock, and every clause here is decided by a hash
//! comparison or an integer ratio.

mod common;

use std::time::Instant;

use common::experiment::Run;
use isomesh::brush::{Brush, BrushOp, BrushStack, Capsule};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::{BoundedSdf, FieldBound, Sphere};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

/// Chunks along each axis. 4³ = 64 chunks, P-39's layout.
const CHUNKS_PER_AXIS: i32 = 4;

/// Cells per chunk axis, so 33 samples per axis.
const CELLS_PER_CHUNK: u32 = 32;

/// World units per cell. The chunk is 4 units across and the world 16.
const CELL_SIZE: f64 = 0.125;

/// The world's minimum corner, so the world is `[-8, 8]³`.
const WORLD_ORIGIN: f64 = -8.0;

/// Brushes in the tape. **Not** reducible: M-341's numbers are for 64, and a
/// shorter tape would break the comparability this whole file is built on.
const BRUSHES: usize = 64;

/// Radius of the solid the brushes carve.
const BASE_RADIUS: f64 = 6.0;

/// ULP of slack added to every enclosure bound. P-39's constant.
const PAD_ULPS: f64 = 16.0;

/// M-341's median survivor count over the 64 chunks, the comparability gate.
const M341_MEDIAN_SURVIVORS: usize = 19;

/// M-341's median survivor *fraction*, to four decimals as published.
const M341_MEDIAN_FRACTION: f64 = 0.2969;

/// The registered C2 bound: the median `necessary / survivors` must be at most
/// this for the "the bound over-keeps" clause to hold.
const C2_BOUND: f64 = 0.75;

/// The registered C3 bound: at least this share of unnecessary survivors must be
/// far from the surface for the over-keep to have a named cause.
const C3_BOUND: f64 = 0.90;

// ── the fixture ─────────────────────────────────────────────────────────────

/// A brush shape in this experiment.
///
/// One enum so the whole stack is a single `&[Brush<Shape>]` slice, which is
/// what makes an ablated tape a shorter slice of the same type. Both variants
/// are exact distance fields, so both declare `l == 1` — and they declare it
/// through the crate's own `BoundedSdf` rather than by a constant written here.
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
/// every run — and, because the seed below is P-39's, the same 64 brushes P-39
/// measured.
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
/// Transcribed from `experiment_p39.rs::tape`, seed included. Every constant
/// here is load-bearing for comparability, not for the mechanism.
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

// ── the bound ───────────────────────────────────────────────────────────────

/// The enclosure of a scalar field over a chunk.
#[derive(Clone, Copy, Debug)]
struct Interval {
    lo: f64,
    hi: f64,
}

impl Interval {
    /// Whether the enclosure keeps the field further than `d` from zero
    /// everywhere in the chunk, and on which side. C3's predicate, split so the
    /// direction lands in the CSV.
    fn far_from_zero(self, d: f64) -> (bool, bool) {
        (self.lo > d, self.hi < -d)
    }
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

/// Which pruning rule to apply.
///
/// P-39 has a second variant, `PruneSmoothLosers`, which is its negative control
/// for the smooth-min asymmetry and prunes in a direction that is *not*
/// bit-exact. It is not this experiment and is not transcribed, so the one
/// variant here is the registered rule.
#[derive(Clone, Copy, Debug)]
enum Policy {
    /// The registered rule. `Add` and `Subtract` prune in the losing direction;
    /// a `SmoothAdd` never prunes, because `b + (a − b)` is not `a`.
    Sound,
}

/// What one pruning pass found.
#[derive(Clone, Copy, Debug, Default)]
struct PruneStats {
    survivors: usize,
    /// `Add`s that provably *win* everywhere in the chunk, so the whole tape
    /// prefix and the base field are dead. Counted, not exploited, exactly as in
    /// P-39: `BrushStack` has no way to say "start from this brush".
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
                // `Policy::Sound` never prunes a `SmoothAdd`, because `smin` at
                // `h == 1` returns `b + (a − b)`, which is not bit-identical to
                // `a`. This binding is the compile-time statement of that: a
                // second `Policy` variant would stop building right here rather
                // than silently inherit the sound rule.
                let Policy::Sound = policy;
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
}

/// Mesh `field` over `grid` into `out` and return `mesh_hash` of the result.
fn hash_of<F: Sdf<Scalar = f64>>(
    mc: &mut MarchingCubes<f64>,
    out: &mut MeshBuffer<f64>,
    field: &F,
    grid: &Grid,
) -> u64 {
    out.reset();
    mc.extract(field, &grid.shape, grid.origin, grid.cell, out)
        .expect("chunk extraction");
    mesh_hash(out)
}

// ── the sweep ───────────────────────────────────────────────────────────────

/// Everything one chunk contributed.
#[derive(Clone, Copy, Debug)]
struct ChunkResult {
    id: [i32; 3],
    survivors: usize,
    dominant_adds: usize,
    /// Survivors whose individual removal changes `mesh_hash`.
    necessary: usize,
    /// Survivors whose enclosure over the chunk stays clear of zero by more than
    /// one cell, whether or not they proved necessary.
    far_survivors: usize,
    /// Unnecessary survivors that are far, and the split of the disjunction.
    unnecessary_far: usize,
    unnecessary_far_by_lo: usize,
    unnecessary_far_by_hi: usize,
    /// The survivors-only reference hash, the full-tape control hash, and the
    /// hash of the `necessary` brushes alone — every individually-unnecessary
    /// survivor removed at once.
    hash: u64,
    hash_full: u64,
    hash_necessary_only: u64,
    vertices: usize,
    triangles: usize,
    remeshes: usize,
    ns_per_remesh: f64,
}

impl ChunkResult {
    fn label(&self) -> String {
        format!("{}-{}-{}", self.id[0], self.id[1], self.id[2])
    }

    fn unnecessary(&self) -> usize {
        self.survivors - self.necessary
    }

    fn necessary_fraction(&self) -> f64 {
        self.necessary as f64 / self.survivors as f64
    }

    /// `unnecessary_far / unnecessary`, and `0` when nothing was unnecessary.
    ///
    /// The registered aggregate is `c3_far_fraction`, which sums numerator and
    /// denominator across chunks; this per-row ratio is the same quantity for
    /// one chunk and exists so a row can be read alone.
    fn unnecessary_far_fraction(&self) -> f64 {
        if self.unnecessary() == 0 {
            0.0
        } else {
            self.unnecessary_far as f64 / self.unnecessary() as f64
        }
    }

    fn control_hash_unchanged(&self) -> bool {
        self.hash == self.hash_full
    }

    /// Whether dropping *every* individually-unnecessary survivor at once still
    /// meshes bit-identically. Decides no registered clause; it is the guard
    /// against reading C2's "individually" as "jointly".
    fn necessary_only_hash_unchanged(&self) -> bool {
        self.hash == self.hash_necessary_only
    }

    fn non_survivors_removed(&self) -> usize {
        BRUSHES - self.survivors
    }
}

/// Buffers that outlive the sweep, so the ablation does not spend its time in
/// the allocator.
struct Rig {
    mc: MarchingCubes<f64>,
    reference: MeshBuffer<f64>,
    scratch: MeshBuffer<f64>,
    survivors: Vec<Brush<Shape>>,
    ablated: Vec<Brush<Shape>>,
    necessary_tape: Vec<Brush<Shape>>,
}

/// The tape and the solid it carves.
struct Fixture {
    base: Sphere<f64>,
    hard: Vec<Brush<Shape>>,
}

fn measure_chunk(
    rig: &mut Rig,
    fixture: &Fixture,
    layout: &ChunkLayout<f64>,
    id: ChunkId,
) -> ChunkResult {
    let Rig {
        mc,
        reference,
        scratch,
        survivors,
        ablated,
        necessary_tape,
    } = rig;

    let origin = layout.sample_origin(id);
    let span = f64::from(layout.cells()) * layout.cell_size();
    let chunk = ChunkBox::new(origin, span);
    let grid = Grid {
        shape: layout.sample_shape().expect("chunk sample grid fits u32"),
        origin,
        cell: layout.cell_size(),
    };
    let cell_size = layout.cell_size();

    let stats = prune_into(
        &fixture.hard,
        &fixture.base,
        chunk,
        Policy::Sound,
        survivors,
    );
    assert!(
        stats.survivors > 0,
        "chunk {:?} pruned to zero survivors, so `necessary / survivors` has no \
         value; the fixture is not P-39's",
        id.coords
    );

    // The survivors-only reference. Not counted as a re-mesh: it is the thing
    // every re-mesh is compared against.
    let hash = hash_of(
        mc,
        reference,
        &BrushStack {
            base: fixture.base,
            brushes: survivors.as_slice(),
        },
        &grid,
    );
    let vertices = reference.vertex_count();
    let triangles = reference.triangle_count();

    let started = Instant::now();

    // C1 first: the whole non-survivor set removed at once.
    let hash_full = hash_of(
        mc,
        scratch,
        &BrushStack {
            base: fixture.base,
            brushes: fixture.hard.as_slice(),
        },
        &grid,
    );

    // C2: leave one survivor out at a time, order of the rest preserved.
    necessary_tape.clear();
    let mut unnecessary_far = 0usize;
    let mut unnecessary_far_by_lo = 0usize;
    let mut unnecessary_far_by_hi = 0usize;
    let mut far_survivors = 0usize;
    for i in 0..survivors.len() {
        ablated.clear();
        ablated.extend_from_slice(&survivors[..i]);
        ablated.extend_from_slice(&survivors[i + 1..]);
        let ablated_hash = hash_of(
            mc,
            scratch,
            &BrushStack {
                base: fixture.base,
                brushes: ablated.as_slice(),
            },
            &grid,
        );

        // C3: the same enclosure the pruner consumed for this brush.
        let (far_lo, far_hi) = enclose(&survivors[i].shape, chunk).far_from_zero(cell_size);
        if far_lo || far_hi {
            far_survivors += 1;
        }

        if ablated_hash == hash {
            if far_lo || far_hi {
                unnecessary_far += 1;
            }
            unnecessary_far_by_lo += usize::from(far_lo);
            unnecessary_far_by_hi += usize::from(far_hi);
        } else {
            necessary_tape.push(survivors[i]);
        }
    }

    let necessary = necessary_tape.len();
    let remeshes = survivors.len() + 1;
    let ns_per_remesh = started.elapsed().as_secs_f64() * 1e9 / remeshes as f64;

    // Outside the timed region and outside `remeshes`, because it answers a
    // question the registration does not ask: are the individually-unnecessary
    // survivors *jointly* droppable? One extra mesh from the necessary brushes
    // alone, their order preserved.
    let hash_necessary_only = hash_of(
        mc,
        scratch,
        &BrushStack {
            base: fixture.base,
            brushes: necessary_tape.as_slice(),
        },
        &grid,
    );

    ChunkResult {
        id: id.coords,
        survivors: stats.survivors,
        dominant_adds: stats.dominant_adds,
        necessary,
        far_survivors,
        unnecessary_far,
        unnecessary_far_by_lo,
        unnecessary_far_by_hi,
        hash,
        hash_full,
        hash_necessary_only,
        vertices,
        triangles,
        remeshes,
        ns_per_remesh,
    }
}

/// Upper median of a sorted-in-place copy, P-39's and P-38's convention:
/// `sorted[len / 2]`, the larger of the two middles for an even count — the
/// conservative direction for an "at most 0.75" clause.
fn median(values: &mut [f64]) -> f64 {
    values.sort_by(f64::total_cmp);
    values[values.len() / 2]
}

/// The aggregates that decide the clauses, repeated on every row.
struct Summary {
    chunks: usize,
    control_all_unchanged: bool,
    control_failures: usize,
    survivors_median: f64,
    necessary_fraction_median: f64,
    necessary_fraction_min: f64,
    necessary_fraction_max: f64,
    total_survivors: usize,
    total_necessary: usize,
    total_unnecessary: usize,
    total_unnecessary_far: usize,
    total_remeshes: usize,
    far_fraction: f64,
    /// Chunks where every individually-unnecessary survivor could go at once.
    joint_drop_chunks_unchanged: usize,
    /// Chunks whose reference mesh is empty, so their `necessary == 0` is free.
    empty_chunks: usize,
    /// C2 and C3 recomputed over non-empty chunks only, so the free zeros can be
    /// divided out rather than argued about.
    necessary_fraction_median_nonempty: f64,
    far_fraction_nonempty: f64,
}

impl Summary {
    fn of(rows: &[ChunkResult]) -> Self {
        let mut fractions: Vec<f64> = rows.iter().map(ChunkResult::necessary_fraction).collect();
        let mut survivors: Vec<f64> = rows.iter().map(|r| r.survivors as f64).collect();
        let total_unnecessary: usize = rows.iter().map(ChunkResult::unnecessary).sum();
        let total_unnecessary_far: usize = rows.iter().map(|r| r.unnecessary_far).sum();
        Self {
            chunks: rows.len(),
            control_all_unchanged: rows.iter().all(ChunkResult::control_hash_unchanged),
            control_failures: rows.iter().filter(|r| !r.control_hash_unchanged()).count(),
            survivors_median: median(&mut survivors),
            necessary_fraction_median: median(&mut fractions),
            necessary_fraction_min: fractions.iter().copied().fold(f64::MAX, f64::min),
            necessary_fraction_max: fractions.iter().copied().fold(f64::MIN, f64::max),
            total_survivors: rows.iter().map(|r| r.survivors).sum(),
            total_necessary: rows.iter().map(|r| r.necessary).sum(),
            total_unnecessary,
            total_unnecessary_far,
            total_remeshes: rows.iter().map(|r| r.remeshes).sum(),
            // Summed numerator over summed denominator: the share of *brushes*,
            // not the mean of per-chunk ratios. Zero unnecessary survivors
            // anywhere would leave C3 with nothing to decide, which
            // `assert_c3_has_a_denominator` refuses to let pass silently.
            far_fraction: total_unnecessary_far as f64 / total_unnecessary as f64,
            joint_drop_chunks_unchanged: rows
                .iter()
                .filter(|r| r.necessary_only_hash_unchanged())
                .count(),
            empty_chunks: rows.iter().filter(|r| r.triangles == 0).count(),
            necessary_fraction_median_nonempty: {
                let mut nonempty: Vec<f64> = rows
                    .iter()
                    .filter(|r| r.triangles > 0)
                    .map(ChunkResult::necessary_fraction)
                    .collect();
                assert!(
                    !nonempty.is_empty(),
                    "every chunk meshed empty, so nothing here is about a surface"
                );
                median(&mut nonempty)
            },
            far_fraction_nonempty: {
                let num: usize = rows
                    .iter()
                    .filter(|r| r.triangles > 0)
                    .map(|r| r.unnecessary_far)
                    .sum();
                let den: usize = rows
                    .iter()
                    .filter(|r| r.triangles > 0)
                    .map(ChunkResult::unnecessary)
                    .sum();
                assert!(den > 0, "no non-empty chunk had an unnecessary survivor");
                num as f64 / den as f64
            },
        }
    }

    fn c1_held(&self) -> bool {
        self.control_all_unchanged
    }

    fn c2_held(&self) -> bool {
        self.necessary_fraction_median <= C2_BOUND
    }

    fn c3_held(&self) -> bool {
        self.far_fraction >= C3_BOUND
    }
}

/// The comparability gate: this fixture must be P-39's fixture.
///
/// M-341's published median is 19 survivors of 64, `0.2969`. If either number
/// moved, the ablation below is measuring a different tape and must not be
/// reported under M-341's name.
fn assert_comparable(summary: &Summary) {
    assert_eq!(
        summary.chunks,
        (CHUNKS_PER_AXIS.pow(3)) as usize,
        "the comparability gate needs all 64 chunks, and this run measured {}",
        summary.chunks
    );
    let median_survivors = summary.survivors_median as usize;
    assert_eq!(
        median_survivors, M341_MEDIAN_SURVIVORS,
        "fixture drift: median survivors is {median_survivors}, but M-341 \
         published {M341_MEDIAN_SURVIVORS} of {BRUSHES}. Every number in this \
         experiment is stated as an extension of M-341 and is meaningless \
         against a different tape"
    );
    let fraction = summary.survivors_median / BRUSHES as f64;
    assert!(
        (fraction - M341_MEDIAN_FRACTION).abs() < 5e-5,
        "fixture drift: median survivor fraction is {fraction:.4}, but M-341 \
         published {M341_MEDIAN_FRACTION}"
    );
}

/// C3 needs unnecessary survivors to be a fraction *of* something.
///
/// A run where every survivor proved necessary would make C2 read `1.0` and
/// leave C3's denominator empty. That is a coherent outcome — the bound would be
/// tight — but it is not a C3 verdict, and dividing by zero into the CSV would
/// look like one.
fn assert_c3_has_a_denominator(summary: &Summary) {
    assert!(
        summary.total_unnecessary > 0,
        "no survivor on any chunk proved unnecessary, so C3's denominator is \
         zero and C3 has no verdict; C2's median is {:.6} against {C2_BOUND}",
        summary.necessary_fraction_median
    );
}

fn row_of(r: &ChunkResult, s: &Summary) -> Vec<(&'static str, String)> {
    vec![
        // ── registered ──
        ("chunk", r.label()),
        ("brushes", BRUSHES.to_string()),
        ("survivors", r.survivors.to_string()),
        ("necessary", r.necessary.to_string()),
        (
            "necessary_fraction",
            format!("{:.9}", r.necessary_fraction()),
        ),
        (
            "non_survivors_removed",
            r.non_survivors_removed().to_string(),
        ),
        (
            "control_hash_unchanged",
            r.control_hash_unchanged().to_string(),
        ),
        (
            "unnecessary_far_from_surface",
            r.unnecessary_far.to_string(),
        ),
        (
            "unnecessary_far_fraction",
            format!("{:.9}", r.unnecessary_far_fraction()),
        ),
        ("mesh_hash", r.hash.to_string()),
        ("remeshes", r.remeshes.to_string()),
        ("ns_per_remesh", format!("{:.1}", r.ns_per_remesh)),
        // ── extra: this chunk, so a row can be read alone ──
        ("unnecessary", r.unnecessary().to_string()),
        ("far_survivors", r.far_survivors.to_string()),
        ("unnecessary_far_by_lo", r.unnecessary_far_by_lo.to_string()),
        ("unnecessary_far_by_hi", r.unnecessary_far_by_hi.to_string()),
        ("mesh_hash_full_tape", r.hash_full.to_string()),
        (
            "mesh_hash_necessary_only",
            r.hash_necessary_only.to_string(),
        ),
        (
            "necessary_only_hash_unchanged",
            r.necessary_only_hash_unchanged().to_string(),
        ),
        ("dominant_adds", r.dominant_adds.to_string()),
        ("vertices", r.vertices.to_string()),
        ("triangles", r.triangles.to_string()),
        ("cell_size", format!("{CELL_SIZE}")),
        ("samples_per_axis", (CELLS_PER_CHUNK + 1).to_string()),
        // ── extra: the aggregates that decide the clauses ──
        ("chunks_measured", s.chunks.to_string()),
        ("c1_control_all_unchanged", s.c1_held().to_string()),
        ("c1_control_failures", s.control_failures.to_string()),
        (
            "c2_necessary_fraction_median",
            format!("{:.9}", s.necessary_fraction_median),
        ),
        (
            "c2_necessary_fraction_min",
            format!("{:.9}", s.necessary_fraction_min),
        ),
        (
            "c2_necessary_fraction_max",
            format!("{:.9}", s.necessary_fraction_max),
        ),
        ("c2_bound", format!("{C2_BOUND}")),
        ("c2_median_at_most_bound", s.c2_held().to_string()),
        ("c3_far_fraction", format!("{:.9}", s.far_fraction)),
        ("c3_bound", format!("{C3_BOUND}")),
        ("c3_far_at_least_bound", s.c3_held().to_string()),
        ("total_survivors", s.total_survivors.to_string()),
        ("total_necessary", s.total_necessary.to_string()),
        ("total_unnecessary", s.total_unnecessary.to_string()),
        ("total_unnecessary_far", s.total_unnecessary_far.to_string()),
        ("total_remeshes", s.total_remeshes.to_string()),
        ("survivors_median", format!("{:.4}", s.survivors_median)),
        ("m341_median_survivors", M341_MEDIAN_SURVIVORS.to_string()),
        (
            "fixture_matches_m341",
            (s.survivors_median as usize == M341_MEDIAN_SURVIVORS).to_string(),
        ),
        (
            "joint_drop_chunks_unchanged",
            s.joint_drop_chunks_unchanged.to_string(),
        ),
        ("empty_chunks", s.empty_chunks.to_string()),
        (
            "c2_necessary_fraction_median_nonempty",
            format!("{:.9}", s.necessary_fraction_median_nonempty),
        ),
        (
            "c3_far_fraction_nonempty",
            format!("{:.9}", s.far_fraction_nonempty),
        ),
    ]
}

fn emit(run: &mut Run, rows: &[ChunkResult], summary: &Summary) {
    for r in rows {
        run.record(&row_of(r, summary));
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-59");

    let fixture = Fixture {
        base: Sphere {
            center: [0.0; 3],
            radius: BASE_RADIUS,
        },
        hard: tape(),
    };
    let adds = fixture.hard.iter().filter(|b| b.op == BrushOp::Add).count();
    let layout = ChunkLayout::new(CELLS_PER_CHUNK, CELL_SIZE, [WORLD_ORIGIN; 3])
        .expect("chunk layout is well formed");

    let mut rig = Rig {
        mc: MarchingCubes::new(),
        reference: MeshBuffer::new(),
        scratch: MeshBuffer::new(),
        survivors: Vec::with_capacity(BRUSHES),
        ablated: Vec::with_capacity(BRUSHES),
        necessary_tape: Vec::with_capacity(BRUSHES),
    };

    common::experiment::run(prereg, |run| {
        println!(
            "tape: {BRUSHES} brushes ({adds} Add, {} Subtract) over a sphere of \
             radius {BASE_RADIUS}; {}³ chunks of {CELLS_PER_CHUNK} cells at \
             {CELL_SIZE} world units",
            BRUSHES - adds,
            CHUNKS_PER_AXIS
        );
        println!(
            "{:>7}  {:>9}  {:>9}  {:>6}  {:>11}  {:>13}  {:>9}",
            "chunk", "survivors", "necessary", "frac", "unnecessary", "unnec_far/unn", "control"
        );

        let mut rows = Vec::with_capacity((CHUNKS_PER_AXIS.pow(3)) as usize);
        for z in 0..CHUNKS_PER_AXIS {
            for y in 0..CHUNKS_PER_AXIS {
                for x in 0..CHUNKS_PER_AXIS {
                    let r = measure_chunk(&mut rig, &fixture, &layout, ChunkId::new([x, y, z]));
                    if !r.control_hash_unchanged() {
                        println!(
                            "  *** VOID: chunk {} — the soundness control FAILED. \
                             survivors-only hash {} != full-tape hash {}. The \
                             interval bound is unsound on this chunk and every \
                             other number on this row is void.",
                            r.label(),
                            r.hash,
                            r.hash_full
                        );
                    }
                    println!(
                        "{:>7}  {:>4}/{BRUSHES}  {:>9}  {:>6.3}  {:>11}  \
                         {:>5}/{:<7}  {:>9}",
                        r.label(),
                        r.survivors,
                        r.necessary,
                        r.necessary_fraction(),
                        r.unnecessary(),
                        r.unnecessary_far,
                        r.unnecessary(),
                        r.control_hash_unchanged(),
                    );
                    rows.push(r);
                }
            }
        }

        let summary = Summary::of(&rows);
        assert_comparable(&summary);

        println!(
            "\ncomparability: median survivors {:.0}/{BRUSHES} ({:.4}) — M-341 \
             published {M341_MEDIAN_SURVIVORS} ({M341_MEDIAN_FRACTION})",
            summary.survivors_median,
            summary.survivors_median / BRUSHES as f64
        );
        println!(
            "C1 soundness control: {}/{} chunks byte-identical after removing all \
             non-survivors at once → {}",
            summary.chunks - summary.control_failures,
            summary.chunks,
            if summary.c1_held() {
                "HELD"
            } else {
                "FALSIFIED — every other number is void"
            }
        );
        println!(
            "C2 necessary/survivors: median {:.6}, min {:.6}, max {:.6} against \
             {C2_BOUND} → {}",
            summary.necessary_fraction_median,
            summary.necessary_fraction_min,
            summary.necessary_fraction_max,
            if summary.c2_held() {
                "HELD"
            } else {
                "FALSIFIED"
            }
        );
        println!(
            "   {} survivors over {} chunks, {} necessary, {} unnecessary",
            summary.total_survivors,
            summary.chunks,
            summary.total_necessary,
            summary.total_unnecessary
        );
        println!(
            "   {} of {} chunks meshed empty, where `necessary == 0` is free",
            summary.empty_chunks, summary.chunks
        );
        println!(
            "   non-empty chunks only: C2 median {:.6}, C3 {:.6}",
            summary.necessary_fraction_median_nonempty, summary.far_fraction_nonempty
        );
        assert_c3_has_a_denominator(&summary);
        println!(
            "C3 unnecessary far from surface: {}/{} = {:.6} against {C3_BOUND} → {}",
            summary.total_unnecessary_far,
            summary.total_unnecessary,
            summary.far_fraction,
            if summary.c3_held() {
                "HELD"
            } else {
                "FALSIFIED"
            }
        );
        println!(
            "   direction split: {} unnecessary survivors far by `lo > cell`, {} \
             by `hi < -cell`",
            rows.iter().map(|r| r.unnecessary_far_by_lo).sum::<usize>(),
            rows.iter().map(|r| r.unnecessary_far_by_hi).sum::<usize>()
        );
        println!(
            "extra, decides nothing registered: dropping every \
             individually-unnecessary survivor AT ONCE leaves the mesh \
             bit-identical on {}/{} chunks — C2 is an `individually` claim and \
             this is the `jointly` one",
            summary.joint_drop_chunks_unchanged, summary.chunks
        );
        println!(
            "extra (gates nothing): {} re-meshes total, {:.2} ms each on the \
             median chunk",
            summary.total_remeshes,
            rows[rows.len() / 2].ns_per_remesh / 1e6
        );

        emit(run, &rows, &summary);
    });
}
