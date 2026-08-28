//! **P-84 — convexity preserved rather than recovered, against `M-116`'s 241 ms.**
//!
//! Ticket: R-084. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p84
//! ```
//!
//! Writes `docs/experiments/p-84.csv`.
//!
//! # The claim, and whose it is
//!
//! Müller, Chentanez & Kim, *Real Time Dynamic Fracture with Volumetric
//! Approximate Convex Decompositions*, TOG 32(4), 2013
//! (`10.1145/2461912.2461934`): split the geometry into non-overlapping convex
//! regions **offline**, then at runtime align a convex fracture pattern to the
//! impact and intersect it against the precomputed compound. The invariant is
//! that **clipping a convex shape against a convex cell yields convex pieces**,
//! so runtime decomposition never happens. Their measurement: *"the time to
//! fracture small to average sized objects is typically negligible, i.e. below
//! 10 ms"*, staying below 50 ms throughout, reaching 20k compounds / 32k
//! convexes — on a **Core i7 @ 3.07 GHz with a GTX 680**. `A-027` already
//! proposes this construction and does not cite that paper.
//!
//! The 2026 state of the art has not improved on it and this harness is the last
//! time that literature gets tracked: `M-297` found no published convex
//! decomposition running at interactive rates, and VisACD (`arXiv:2604.04244`),
//! a **GPU** method, averages **16.97 s per model** against CoACD's 36.31 s —
//! four orders of magnitude off a frame budget.
//!
//! # The SHARE line, recomputed against `p-85.csv` before this ran
//!
//! The registration says *"SHARE: C1 moves the whole destruction-collider
//! stage."* `R-085` has since attributed the collider's 45%, and its answer is
//! against this row on cost. From `docs/experiments/p-85.csv`, over its six
//! rows: `handoff_share` **0.795854–0.815711** (the crate's own weld plus
//! `collider::readiness`), `bvh_share` **0.179151–0.199501** (parry3d's
//! `QbvhBuilder`), `copy_share` **0.000354–0.000565** — the triangle copy is
//! 0.035–0.056% of the chunk collider build. **There is no convex-decomposition
//! term anywhere in that 45%**, because the crate does not do convex
//! decomposition: the README's line is *"convex decomposition — not yet"*.
//! `R-085`'s conclusion, *"do not build `R-084` for cost"*, is therefore correct
//! about the chunk collider and says nothing about the stage this row's C1 is
//! denominated in.
//!
//! So, clause by clause, before the run:
//!
//! | clause | share it can move | reachable? |
//! |---|---|---|
//! | **C1** | `M-116`'s destruction-collider stage: `avian3d` 0.7's `convex_decomposition_from_mesh` at **240.7 / 271.8 / 249.0 ms mean per fragment**, worst 323.7 / 369.0 / 362.6 ms, over 23–24 fragments per target. 100% of that stage — there is no second path. **Not one microsecond of `R-085`'s 45%.** | yes, both ways: 240.7/10 = 24.1x and 240.7/25 = 9.6x, and the falsifier at 25 ms is inside a plausible range for a plane-clip arrangement |
//! | **C2** | **none — it is a determinism clause with no cost content.** Unaffected by `R-085`. This is the clause that survives. | yes, both ways: see below — the natural construction's vertex arithmetic is order-dependent *a priori*, and the hash's sensitivity is asserted |
//! | **C3** | **none — it is a memory bound**, a guard on C1's mechanism rather than a share of anything | yes, both ways, but only on a fixture where `convex_cells_before` is not 1; see the fixture note |
//!
//! Stated plainly: **the cost case for this row is not the 45%, it never was,
//! and after `R-085` no clause of P-84 moves any part of it.** What survives is
//! (a) C1 against `M-116`'s own stage, which is unbuilt and whose only measured
//! cost is 241–272 ms per fragment, and (b) C2, which is about replication and
//! not about time.
//!
//! # What the instrument is
//!
//! There is no convex decomposition in the crate (`A-026`/`A-027` are open, and
//! `A-026` is blocked on a scope decision), so this bench builds the substrate it
//! measures. It is *not* a convex decomposition of a mesh — it is the
//! construction `A-027` proposes, which is why it is cheap:
//!
//! 1. **The chunk supplies the outer convex cells for free.** A chunk is a grid
//!    of `chunk_cells³` cubes and a cube is convex. Müller has to compute this
//!    offline with ACD; a voxel field already has it.
//! 2. **The unedited surface contributes one plane per boundary cell.** For a
//!    cell whose eight corners straddle zero, the plane through the local zero
//!    crossing with normal `∇f/|∇f|`. `cube ∩ half-space` is one convex piece,
//!    so an unedited chunk's compound is `full cells + boundary cells` — this is
//!    what a marching-cubes-derived convex decomposition looks like, and it is
//!    why `convex_cells_before` is a real number rather than 1. **Had the
//!    partition been generated from brush planes alone, an unedited chunk would
//!    be one cell and C3's ratio would be unbounded by fixture rather than by
//!    measurement**; that is a fixture decision and it is recorded here.
//! 3. **Each brush contributes its `k`-DOP supporting planes**, from the shape's
//!    exact support function `h(d) = max_{x ∈ S} d·x` — sphere `c·d + r|d|`, box
//!    `c·d + Σ b_i|d_i|`, capsule `max(a·d, b·d) + r|d|` (Klosowski, Held,
//!    Mitchell, Sowizral & Zikan, *Efficient collision detection using bounding
//!    volume hierarchies of k-DOPs*, IEEE TVCG 4(1), 1998). The `k`-DOP is the
//!    tightest convex polytope containing the shape in that direction set, so
//!    these are tangent planes and the enclosure is exact for a box brush.
//! 4. **The partition of one cell is the arrangement of its cube by the planes
//!    that cross it**, built breadth-first: split every current piece by the next
//!    plane, drop the empty ones. Every leaf is an intersection of half-spaces,
//!    hence convex — the paper's invariant, and `convexity_violations` checks it
//!    rather than assuming it.
//! 5. **Each leaf is classified by the folded field at its centroid.** `Add` is
//!    `min`, `Subtract` is `max(-s)`; the compound is the union of solid leaves.
//!
//! A convex piece is carried as its **vertex set, each vertex tagged with the
//! sorted ids of the three planes that define it**. A convex polytope's 1-skeleton
//! is recoverable from that tagging — two vertices are adjacent exactly when they
//! share two planes — so clipping needs no face-cap reconstruction: drop the
//! outside vertices, add one vertex per straddling edge. `degenerate_vertices`
//! asserts every vertex really has three distinct planes, which is what makes the
//! adjacency rule sound.
//!
//! **C1's fracture** is Müller's runtime step exactly: the impact brush's `k`-DOP
//! is the fracture pattern, and for every cell it overlaps, every solid piece
//! already in the compound is split by the pattern's planes. A leaf inside the
//! pattern is a **fragment** — already convex, no decomposition, which is the
//! whole claim. No field is sampled during the timed region: a leaf inherits its
//! parent's solid classification, which is what the invariant buys.
//!
//! **C2 has two design forks, not one, and the second is the finding.** The
//! registration anticipated only the first.
//!
//! *Fork one: how a crossing vertex is computed.*
//!
//! - `vertex_mode = lerp` — Sutherland–Hodgman's own arithmetic, `v = a +
//!   t(b−a)` with `t = d_a/(d_a − d_b)`. The endpoints `a` and `b` were
//!   themselves produced by earlier clips, so the value depends on **which
//!   plane was applied first**. This is the natural implementation and the one
//!   `A-027` describes.
//! - `vertex_mode = solve` — the vertex is the exact solution of its three plane
//!   equations by Cramer's rule, with the triple in ascending id order. The
//!   plane ids are the **stable edit ids** (each edit's own identity, not its
//!   position in the log — an edit log entry is identifiable, which is the
//!   premise `M-36` is about), so the position is a function of the triple alone.
//!
//! *Fork two: whether coincident planes are collapsed.* Two brushes can
//! contribute the **same plane**, and `M-36`'s eight-brush fixture does it three
//! times per chunk (`coincident_plane_pairs = 3`, in three cells). A half-space
//! appearing twice adds nothing to an intersection of half-spaces, so it cannot
//! move the geometry — but it leaves a vertex on that plane with **two ids it
//! could equally be tagged with**, and which one it gets is decided by the
//! traversal.
//!
//! - `planes_mode = as_logged` — every brush's DOP planes as the log gives them.
//! - `planes_mode = canonical` — a cut plane is dropped when a *smaller stable
//!   id* in the same cell bounds the same plane in either orientation. Order-free
//!   because the ids are, and it deliberately does **not** reorder the survivors:
//!   `arrange` still applies them in the log's order, so the sweep still asks
//!   whether the arrangement depends on that order instead of being handed the
//!   answer. `distinct_partitions_raw = 18720` in every arm is the proof that it
//!   was still asked — the traversal really did take 18,720 distinct shapes.
//!
//! All four arms run all 40,320 orderings. Whether any is 1 is the measurement.
//!
//! # Determinism instrumentation, and why the zero could have been non-zero
//!
//! `partition_hash` is FNV-1a over, per hashed cell: the cell index, the piece
//! count, and the sorted multiset of per-piece hashes, each piece hashed over its
//! vertices sorted by plane triple — three `u16` ids and three `f64` bit
//! patterns each — plus the piece's solid flag. Sorting isolates *arithmetic*
//! order-dependence from emission order, which is the question C2 asks; the
//! emission-order hash is reported beside it as `partition_hash_raw` /
//! `distinct_partitions_raw` so the byte-level replication question is not
//! reconstructed from the other one.
//!
//! **A count of distinct partitions above one names no mechanism**, so the hash
//! is also reported as three nested projections, each adding exactly one kind of
//! information to the one before it:
//!
//! - `distinct_topology` — the leaf plane-triple sets alone, no coordinates and
//!   no solid flags. Is the *arrangement* order-free?
//! - `distinct_positions` — those plus the vertex bit patterns. Is the
//!   *arithmetic* order-free?
//! - `distinct_solid_flags` — the triples plus the solid classification. Is the
//!   *union* order-free?
//!
//! Beside them, `topo_class_lo_pieces` / `topo_class_hi_pieces` and their solid
//! counterparts rebuild a representative ordering of each extreme topology class
//! and count its leaves, because "the same leaves labelled differently" and "the
//! arrangement lost a leaf" are different failures and the worse one has to be
//! nameable.
//!
//! Two `M-44` controls, because `distinct` and `distinct_topology` are zeros over
//! two different instruments and one control cannot license both.
//! `perturbed_hash_differs` displaces brush 0's centre by `1e-9`: far below a
//! cell, so it moves no plane across a cell boundary and can only show up in the
//! coordinates. `perturbed_topology_differs` displaces it by a quarter of a cell,
//! which does change which cells the planes cross. Without the second,
//! `distinct_topology == 1` would be `P-70`'s C3 — a held clause with no
//! instrument.
//!
//! # Controls
//!
//! - **`orderings == 40320`**, asserted, and the 40,320 permutations asserted
//!   pairwise distinct (`permutations_distinct`). The registered vacuity control:
//!   the arm must actually reach every ordering, which is `M-36`'s own fixture.
//! - **`partition_volume_error_rel`** — for every cell, the leaf volumes must sum
//!   to the cube's volume. That is the paper's *non-overlapping* half of the
//!   invariant: overlapping pieces would sum high, a gap would sum low. Asserted
//!   below `1e-9` relative.
//! - **`convexity_violations == 0`** — every vertex of every piece is on the
//!   piece's own side of every plane incident to it, to `1e-9`. That is the
//!   *convex* half: this is the assertion that "clipping a convex shape against a
//!   convex cell yields convex pieces" held on this data rather than in the
//!   abstract.
//! - **`degenerate_vertices == 0`** and **`singular_triples == 0`** — no vertex
//!   carries a repeated plane, and no plane triple was near-singular, so the
//!   adjacency rule and Cramer's rule were both applicable everywhere.
//! - **`fragments > 0`**, **`convex_cells_before > 0`**, **`cells_cut_by_brushes
//!   > 0`** — C1 and C3 both had a non-empty population to fail on.
//! - **`volume_error_rel`** — the compound's volume against the folded field's
//!   volume sampled on a `VOLUME_SAMPLES³` lattice. **This is the price of the
//!   comparison with `M-116` and it is not free**: `avian3d` decomposes an
//!   arbitrary fragment mesh, while this substrate approximates the solid to the
//!   fidelity of a `k`-DOP-cut grid. The column says how much.
//! - **`cycles` / `ghz`** — `M-280`: on a governed CPU a millisecond is not a
//!   unit. The counted fracture rep's cycles are reported beside the clock, and
//!   `cycles_per_fragment` is clock-independent.
//! - **`coincident_plane_pairs`, `nonsimple_points`, `nonsimple_cells`,
//!   `max_planes_at_a_point`** — defect 4's mechanism, measured instead of
//!   assumed. Scanned per cell over the cube's six faces and the cell's cut set:
//!   pairs bounding the same plane, and points in the cube where four or more
//!   planes meet. Order-independent by construction — it reads the cut *set*,
//!   never a traversal — so it runs once per arm, outside the sweep.
//!
//! # Deviations, and the six defects these controls caught
//!
//! Every one of these was found by a control in this file, not by inspection,
//! and each is recorded because the first three make the *first* version of
//! this harness's numbers worthless rather than merely noisy.
//!
//! 1. **`cube_verts` aliased its cube-face plane ids.** Corners were tagged
//!    `i*2, 2+j*2, 4+k*2`, which puts the x-max face on id 2 — the same id
//!    `set_cube` writes the y-min face to — the y-max face on the z-min's id,
//!    and the z-max face on `SURFACE_ID`. Every corner on a far face then
//!    carried a repeated plane id, `shared_two` reported false edges between
//!    unrelated vertices, and the false vertices bred: one piece reached
//!    **8.4 million vertices in a single 268 MB allocation**. Caught by
//!    `degenerate_vertices`.
//! 2. **A split half was pruned with the wrong orientation.** The first repair
//!    for (1) tested each candidate against *every* plane the piece's tags
//!    mention as `n·x + d <= 0`. The outside half of a split is bounded by the
//!    cut plane the other way round, so that test deleted it whole. Caught by
//!    the per-cell volume closure, at a relative error of exactly 1.0.
//! 3. **The compound did not carry each cell's surface plane.** Ids `0..=6` are
//!    rewritten per cell; `fracture` restored the cube but not `SURFACE_ID`, so
//!    every piece bounded by the surface was measured against another cell's
//!    surface plane and pruned to nothing. Caught by `fragments > 0`.
//! 4. **The 14-DOP arm is not in this file, and the volume control is why.**
//!    With the eight cube-diagonal directions added, cells lost leaves: one
//!    16³ `fbm_terrain` cell closed at `1.3495e-2` against a cube of
//!    `1.5625e-2` — **13.6% of the cube missing** — and another emitted
//!    duplicated leaves. The vertex-tag clipper's "two vertices are adjacent
//!    exactly when they share two planes" rule is sound for a **simple**
//!    polytope, and an oblique 14-plane set makes non-simple ones. Rather than
//!    report a piece count from an arrangement that does not tile its own cube,
//!    the arm was cut. **Every row in `p-84.csv` is `dop_dirs = 6`** — the
//!    axis-aligned DOP, which is the *exact* brush for a box and a conservative
//!    enclosure for a sphere or a capsule. A 14-DOP arm needs a clipper that
//!    derives vertices from plane triples rather than from edges, and that is
//!    the change, not a tolerance.
//! 5. **The canonical sort ran after the sums it was supposed to canonicalise.**
//!    The first version sorted each leaf's vertices by plane triple *for the
//!    hash only*, leaving `centroid` and `volume` summing in the breadth-first
//!    split's emission order. Both are naive float sums, so both were
//!    ULP-order-dependent, and both feed a decision: the centroid picks the
//!    leaf's solid flag through `fold`, the volume picks whether the leaf is a
//!    sliver at all. That is two order dependencies **manufactured by the
//!    instrument**, standing in front of the exact question C2 asks. The sort now
//!    runs first. Caught by adding the three projections and finding
//!    `distinct_positions` and `distinct_solid_flags` disagreeing for reasons the
//!    single number could not distinguish.
//! 6. **The 6-DOP arm is not immune to defect 4 either, and the fixture proves
//!    it.** `simplicity_scan` finds three coincident plane pairs in `M-36`'s
//!    eight-brush fixture, in three of the 182 cut cells, and 36 points where
//!    four planes meet — `max_planes_at_a_point = 4`. Those 36 points are exactly
//!    the 36 zero-width leaves `slivers_dropped` reports, and collapsing the
//!    coincident planes takes all four counters to 0. Axis-aligned does not mean
//!    generic: two brushes agreeing on one tangent plane is not a coincidence
//!    worth calling unlikely when brush centres and radii are authored by hand.
//!
//! One fixture caveat that is not a defect: `BRUSH_RADIUS_CELLS` is in **cells**,
//! so the 32³ arms dig brushes of half the world radius the 16³ arms dig. That
//! is the natural scaling for a voxel editor — a brush is so many cells wide —
//! but it means C3's two resolutions are not the same hole, and the ratio must
//! be read with `cell_size_world` beside it.
//!
//! A second fixture caveat, and it is C1's: **`fracture_ms_worst` is a
//! single-piece clock on a machine with three sibling builds running**, and it
//! swung 0.055 → 1.66 ms across two runs of the same binary while the median per
//! fragment moved 0.0028 → 0.0050. The medians and the integer counts reproduce;
//! `fracture_ms_worst` is the one column here that should be re-read from a quiet
//! tree (`M-281`).
//!
//! # References
//!
//! Müller, Chentanez & Kim 2013, `10.1145/2461912.2461934`. Klosowski et al.
//! 1998 for the `k`-DOP and its support function. Sutherland & Hodgman,
//! *Reentrant polygon clipping*, CACM 17(1), 1974, for the `lerp` arm's
//! arithmetic. `M-116` (E-204) for the 241–272 ms this is measured against,
//! `M-36` (G-003) for the 40,320-ordering fixture, `M-50` (E-202) for the
//! 46–60-brush bucket, `M-297` and `p-85.csv` for what this row cannot claim.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::needless_range_loop,
    clippy::print_literal,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::Sdf;
use isomesh::brush::{Brush, BrushOp, Capsule, apply};
use isomesh::fields::{BoxExact, FbmTerrain, Gyroid, Sphere};

use crate::common::counters::{MIN_TIME_RATIO, Probe};

/// Brushes in `M-36`'s fixture.
const BRUSHES: usize = 8;

/// `8!`. The registered vacuity control is a count against this.
const ORDERINGS: usize = 40_320;

/// `M-50`'s largest bucket: 46–60 brushes on one chunk. C3's fixture.
const SIXTY: usize = 60;

/// C1's registered bar, in milliseconds per fragment. Müller's own number.
const C1_BAR_MS: f64 = 10.0;

/// C1's registered falsifier, in milliseconds per fragment.
const C1_FALSIFIER_MS: f64 = 25.0;

/// C3's registered bar: cells after / cells before.
const C3_BAR: f64 = 4.0;

/// `M-116`'s mean cost per fragment, in ms — the low end of its three targets
/// (wall 240.7, hollow shell 271.8, spiral 249.0).
const M116_MEAN_LOW_MS: f64 = 240.7;

/// `M-116`'s mean cost per fragment, in ms — the high end.
const M116_MEAN_HIGH_MS: f64 = 271.8;

/// `M-116`'s worst single fragment, in ms — the worst of its three targets.
const M116_WORST_MS: f64 = 369.0;

/// Repetitions of the timed fracture. Odd, so the median is a sample.
const REPS: usize = 9;

/// Lattice per axis for the sampled true solid volume, for `volume_error_rel`.
const VOLUME_SAMPLES: u32 = 96;

/// A leaf whose volume is below this fraction of its cube's is a sliver from a
/// plane grazing a corner, not a piece of the compound.
const SLIVER_FRACTION: f64 = 1e-12;

/// Tolerance on the two invariant checks, relative to the cell size.
const GEOM_EPS: f64 = 1e-9;

/// Relative slack on the per-cell volume closure. A cell's leaves must tile its
/// cube to this, or the arrangement lost or duplicated a piece.
const VOLUME_CLOSURE_REL: f64 = 1e-9;

/// Leaves one cell's arrangement may hold at once before the build is declared
/// broken rather than slow. The arrangement of `p` planes in a cube has at most
/// `C(p,3)+C(p,2)+p+1` cells, so this is reached only if the split is
/// manufacturing zero-volume parts and doubling them.
const LEAF_CAP: usize = 50_000;

/// Threads. Capped so the harness's numbers do not depend on the host's core
/// count — only its runtime does.
const MAX_THREADS: usize = 12;

/// Brush radius for the 60-brush arm, in cells.
const BRUSH_RADIUS_CELLS: f64 = 1.5;

/// Impact-brush radius for C1's fracture, in cells.
const IMPACT_RADIUS_CELLS: f64 = 3.0;

// ---------------------------------------------------------------------------
// Half-spaces and convex pieces
// ---------------------------------------------------------------------------

/// A half-space. **Inside is `n·x + d <= 0`.**
#[derive(Clone, Copy)]
struct Plane {
    n: [f64; 3],
    d: f64,
}

impl Plane {
    /// `n·p + d`. Negative inside.
    #[inline]
    fn at(&self, p: [f64; 3]) -> f64 {
        self.n[0] * p[0] + self.n[1] * p[1] + self.n[2] * p[2] + self.d
    }

    /// `|n|`, for turning `at` into a distance.
    #[inline]
    fn norm(&self) -> f64 {
        (self.n[0] * self.n[0] + self.n[1] * self.n[1] + self.n[2] * self.n[2]).sqrt()
    }
}

/// A vertex of a convex piece, tagged with the three planes that define it.
///
/// The tag is what makes the 1-skeleton recoverable without storing faces: two
/// vertices of a convex polytope are adjacent exactly when they share two
/// planes.
#[derive(Clone, Copy)]
struct Vert {
    p: [f64; 3],
    /// Plane ids, ascending. Ids are **stable** — a cube face, the surface
    /// plane, or `(edit id, DOP direction)` — never a position in the log.
    t: [u16; 3],
}

/// How a crossing vertex's position is computed. C2's fork.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    /// Sutherland–Hodgman: lerp along the straddling edge. Order-dependent by
    /// construction, which is the thing being measured.
    Lerp,
    /// Cramer's rule on the vertex's own three planes, triple in ascending
    /// stable-id order.
    Solve,
}

impl Mode {
    const fn name(self) -> &'static str {
        match self {
            Self::Lerp => "lerp",
            Self::Solve => "solve",
        }
    }
}

/// Whether a cell's plane list is deduplicated by geometric plane.
///
/// The second design fork in `A-027`, and the one the measurement found. Two
/// brushes can contribute the *same* plane — the eight-brush fixture does it
/// three times per chunk — and a half-space appearing twice adds nothing to an
/// intersection of half-spaces. What it does add is an **ambiguity in the name
/// of a vertex**: a vertex on that plane can be tagged with either id, and which
/// one it gets is decided by the traversal.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Planes {
    /// Every brush's DOP planes as the log gives them.
    AsLogged,
    /// Coincident planes collapsed to the smallest **stable** id of the group,
    /// which is order-free because the ids are.
    Canonical,
}

impl Planes {
    const fn name(self) -> &'static str {
        match self {
            Self::AsLogged => "as_logged",
            Self::Canonical => "canonical",
        }
    }
}

/// The two planes shared by two vertices, if there are exactly two.
///
/// Exactly two means the pair spans an edge. Three means the two vertices are
/// the same vertex; fewer means they are not adjacent.
#[inline]
fn shared_two(a: [u16; 3], b: [u16; 3]) -> Option<[u16; 2]> {
    let mut out = [0u16; 2];
    let mut n = 0usize;
    for x in a {
        if b.contains(&x) {
            if n < 2 {
                out[n] = x;
            }
            n += 1;
        }
    }
    if n == 2 { Some(out) } else { None }
}

/// `det` of a matrix given by rows.
#[inline]
fn det3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

/// The point where three planes meet, by Cramer's rule.
///
/// `None` when the triple is near-singular. Counted rather than worked around:
/// `singular_triples` is asserted zero, because a singular triple means the
/// fixture put three planes through one line and the arrangement is not the
/// arrangement anyone meant.
#[inline]
fn solve3(a: &Plane, b: &Plane, c: &Plane) -> Option<[f64; 3]> {
    let det = det3(a.n, b.n, c.n);
    if det.abs() < 1e-12 {
        return None;
    }
    let r = [-a.d, -b.d, -c.d];
    let x = det3(
        [r[0], a.n[1], a.n[2]],
        [r[1], b.n[1], b.n[2]],
        [r[2], c.n[1], c.n[2]],
    );
    let y = det3(
        [a.n[0], r[0], a.n[2]],
        [b.n[0], r[1], b.n[2]],
        [c.n[0], r[2], c.n[2]],
    );
    let z = det3(
        [a.n[0], a.n[1], r[0]],
        [b.n[0], b.n[1], r[1]],
        [c.n[0], c.n[1], r[2]],
    );
    Some([x / det, y / det, z / det])
}

/// Whether two half-spaces bound the **same plane**, in either orientation.
///
/// Same orientation means the second half-space is redundant outright; opposite
/// means the pair contributes one zero-width leaf and nothing else. Both leave a
/// vertex on that plane with two ids it could equally be tagged with, so both
/// count. `eps` is a distance, so the normals are compared after normalisation.
#[inline]
fn coincident(a: &Plane, b: &Plane, eps: f64) -> bool {
    let (la, lb) = (a.norm(), b.norm());
    if la <= 0.0 || lb <= 0.0 {
        return false;
    }
    let same = (0..3).all(|c| (a.n[c] / la - b.n[c] / lb).abs() <= 1e-12)
        && (a.d / la - b.d / lb).abs() <= eps;
    let flip = (0..3).all(|c| (a.n[c] / la + b.n[c] / lb).abs() <= 1e-12)
        && (a.d / la + b.d / lb).abs() <= eps;
    same || flip
}

/// Counters a partition build accumulates and then has asserted.
#[derive(Clone, Copy, Default)]
struct Flags {
    degenerate_verts: u64,
    singular_triples: u64,
    slivers: u64,
    convexity_violations: u64,
    max_verts: usize,
    worst_volume_error: f64,
    /// Most leaves any one cell's arrangement held at once.
    max_leaves: usize,
    /// Most planes any one cell was cut by.
    max_cuts: usize,
}

impl Flags {
    fn merge(&mut self, o: Self) {
        self.degenerate_verts += o.degenerate_verts;
        self.singular_triples += o.singular_triples;
        self.slivers += o.slivers;
        self.convexity_violations += o.convexity_violations;
        self.max_verts = self.max_verts.max(o.max_verts);
        self.worst_volume_error = self.worst_volume_error.max(o.worst_volume_error);
        self.max_leaves = self.max_leaves.max(o.max_leaves);
        self.max_cuts = self.max_cuts.max(o.max_cuts);
    }
}

/// Reduce a split half to the vertex set of a convex polytope.
///
/// Two repairs, and the arrangement does not terminate without the second.
///
/// **Every vertex must lie on all three of its own planes**, to `eps` in
/// distance. A tag is a claim; this is the check. Nothing here tests an
/// *inequality*, and that is deliberate: the two halves of a split bound
/// themselves with opposite signs of the cut plane, so a single "inside every
/// plane the tags mention" test silently deletes the whole outside half. It
/// did — one 16³ `fbm_terrain` cell reported its single surviving leaf at
/// `8.6e-5` against a cube of `1.5625e-2`, a relative volume error of 1.0.
///
/// **Dedup by plane triple.** A vertex of a convex polytope *is* its triple of
/// planes, so a triple appearing twice is one vertex written twice. The
/// "sharing two planes means adjacent" shortcut is sound for a **simple**
/// polytope and not for a degenerate one: when four planes meet at a point,
/// three tagged vertices sit on one plane pair, two of them at the same
/// position, and both chords to the third produce the *same* crossing vertex.
/// Left in, each duplicate spawns a further pair at the next plane and the
/// piece doubles per cut — measured as a single 268 MB, 8.4-million-vertex
/// allocation on the first `fbm_terrain` arm before this existed. With the
/// dedup a piece cut by `p` planes cannot exceed `C(p,3)` vertices.
fn tidy(table: &[Plane], v: &mut Vec<Vert>, eps: f64) {
    v.retain(|x| {
        x.t.iter().all(|&q| {
            let pl = &table[q as usize];
            pl.at(x.p).abs() <= eps * pl.norm()
        })
    });
    v.sort_unstable_by(|a, b| {
        a.t.cmp(&b.t)
            .then_with(|| a.p[0].total_cmp(&b.p[0]))
            .then_with(|| a.p[1].total_cmp(&b.p[1]))
            .then_with(|| a.p[2].total_cmp(&b.p[2]))
    });
    v.dedup_by(|a, b| a.t == b.t);
}

/// Split a convex piece by one half-space, into the two convex pieces.
///
/// The `<=` / `>=` split keeps a vertex that lands exactly on the plane in both
/// halves with its original tag, so no duplicate vertex is manufactured there.
/// Crossing vertices are created once per straddling edge and appear in both
/// halves, which is correct: they are on the shared face. Both halves go
/// through [`tidy`] before they are returned.
fn split(
    table: &[Plane],
    src: &[Vert],
    id: u16,
    mode: Mode,
    eps: f64,
    sd: &mut Vec<f64>,
    inside: &mut Vec<Vert>,
    outside: &mut Vec<Vert>,
    flags: &mut Flags,
) {
    let pl = &table[id as usize];
    // Snap to the plane inside `eps` of it. A vertex that already lies on this
    // plane has a signed distance of a few ULP either way, and reading that as
    // a strict sign manufactures a "crossing" whose triple repeats a plane —
    // 5,697 of them on the first `fbm_terrain` arm before this snap existed.
    let tol = eps * pl.norm();
    sd.clear();
    for v in src {
        let s = pl.at(v.p);
        sd.push(if s.abs() <= tol { 0.0 } else { s });
    }
    inside.clear();
    outside.clear();
    for (i, v) in src.iter().enumerate() {
        if sd[i] <= 0.0 {
            inside.push(*v);
        }
        if sd[i] >= 0.0 {
            outside.push(*v);
        }
    }
    for i in 0..src.len() {
        if sd[i] >= 0.0 {
            continue;
        }
        for j in 0..src.len() {
            if sd[j] <= 0.0 {
                continue;
            }
            let Some(sh) = shared_two(src[i].t, src[j].t) else {
                continue;
            };
            let mut t = [sh[0], sh[1], id];
            t.sort_unstable();
            let p = match mode {
                Mode::Lerp => {
                    let f = sd[i] / (sd[i] - sd[j]);
                    [
                        src[i].p[0] + f * (src[j].p[0] - src[i].p[0]),
                        src[i].p[1] + f * (src[j].p[1] - src[i].p[1]),
                        src[i].p[2] + f * (src[j].p[2] - src[i].p[2]),
                    ]
                }
                Mode::Solve => {
                    match solve3(
                        &table[t[0] as usize],
                        &table[t[1] as usize],
                        &table[t[2] as usize],
                    ) {
                        Some(p) => p,
                        None => {
                            flags.singular_triples += 1;
                            let f = sd[i] / (sd[i] - sd[j]);
                            [
                                src[i].p[0] + f * (src[j].p[0] - src[i].p[0]),
                                src[i].p[1] + f * (src[j].p[1] - src[i].p[1]),
                                src[i].p[2] + f * (src[j].p[2] - src[i].p[2]),
                            ]
                        }
                    }
                }
            };
            if t[0] == t[1] || t[1] == t[2] {
                flags.degenerate_verts += 1;
            }
            inside.push(Vert { p, t });
            outside.push(Vert { p, t });
        }
    }
    tidy(table, inside, eps);
    tidy(table, outside, eps);
    flags.max_verts = flags.max_verts.max(inside.len()).max(outside.len());
}

/// Centroid of a vertex set. Interior for a convex piece.
fn centroid(v: &[Vert]) -> [f64; 3] {
    let n = v.len() as f64;
    let mut c = [0.0f64; 3];
    for x in v {
        c[0] += x.p[0];
        c[1] += x.p[1];
        c[2] += x.p[2];
    }
    [c[0] / n, c[1] / n, c[2] / n]
}

/// Volume of a convex piece, by the divergence theorem over its faces.
///
/// A face is the vertex subset sharing one plane id; the piece is convex, so
/// sorting that subset by angle about its own centroid in the plane's basis is
/// its polygon order, and the shoelace cross-product sum is twice its area. The
/// centroid of the piece is interior, so every face's distance to it is positive
/// and no outward orientation has to be tracked: `V = Σ A_f h_f / 3`.
fn volume(table: &[Plane], v: &[Vert], scratch: &mut Vec<(f64, [f64; 3])>) -> f64 {
    if v.len() < 4 {
        return 0.0;
    }
    let g = centroid(v);
    let mut ids: Vec<u16> = v.iter().flat_map(|x| x.t).collect();
    ids.sort_unstable();
    ids.dedup();
    let mut vol = 0.0f64;
    for id in ids {
        let pl = &table[id as usize];
        let nl = pl.norm();
        if nl <= 0.0 {
            continue;
        }
        let nh = [pl.n[0] / nl, pl.n[1] / nl, pl.n[2] / nl];
        // A basis of the plane, from whichever axis is least aligned with it.
        let axis = if nh[0].abs() <= nh[1].abs() && nh[0].abs() <= nh[2].abs() {
            [1.0, 0.0, 0.0]
        } else if nh[1].abs() <= nh[2].abs() {
            [0.0, 1.0, 0.0]
        } else {
            [0.0, 0.0, 1.0]
        };
        let dot = nh[0] * axis[0] + nh[1] * axis[1] + nh[2] * axis[2];
        let mut u = [
            axis[0] - dot * nh[0],
            axis[1] - dot * nh[1],
            axis[2] - dot * nh[2],
        ];
        let ul = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        u = [u[0] / ul, u[1] / ul, u[2] / ul];
        let w = [
            nh[1] * u[2] - nh[2] * u[1],
            nh[2] * u[0] - nh[0] * u[2],
            nh[0] * u[1] - nh[1] * u[0],
        ];
        scratch.clear();
        for x in v {
            if !x.t.contains(&id) {
                continue;
            }
            scratch.push((0.0, x.p));
        }
        if scratch.len() < 3 {
            continue;
        }
        let mut fc = [0.0f64; 3];
        for (_, p) in scratch.iter() {
            fc[0] += p[0];
            fc[1] += p[1];
            fc[2] += p[2];
        }
        let m = scratch.len() as f64;
        fc = [fc[0] / m, fc[1] / m, fc[2] / m];
        for e in scratch.iter_mut() {
            let r = [e.1[0] - fc[0], e.1[1] - fc[1], e.1[2] - fc[2]];
            let a = r[0] * u[0] + r[1] * u[1] + r[2] * u[2];
            let b = r[0] * w[0] + r[1] * w[1] + r[2] * w[2];
            e.0 = b.atan2(a);
        }
        scratch.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        let mut area2 = 0.0f64;
        for i in 0..scratch.len() {
            let p0 = scratch[i].1;
            let p1 = scratch[(i + 1) % scratch.len()].1;
            let a = [p0[0] - fc[0], p0[1] - fc[1], p0[2] - fc[2]];
            let b = [p1[0] - fc[0], p1[1] - fc[1], p1[2] - fc[2]];
            let cr = [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ];
            area2 += (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
        }
        let h = pl.at(g).abs() / nl;
        vol += area2 * 0.5 * h / 3.0;
    }
    vol
}

/// Every vertex of a piece on the piece's own side of every plane incident to
/// it. The paper's convexity invariant, checked rather than assumed.
fn convexity_violations(table: &[Plane], v: &[Vert], eps: f64) -> u64 {
    let g = centroid(v);
    let mut ids: Vec<u16> = v.iter().flat_map(|x| x.t).collect();
    ids.sort_unstable();
    ids.dedup();
    let mut bad = 0u64;
    for id in ids {
        let pl = &table[id as usize];
        let nl = pl.norm();
        if nl <= 0.0 {
            continue;
        }
        let side = pl.at(g) / nl;
        for x in v {
            let s = pl.at(x.p) / nl;
            if side < 0.0 {
                if s > eps {
                    bad += 1;
                }
            } else if s < -eps {
                bad += 1;
            }
        }
    }
    bad
}

// ---------------------------------------------------------------------------
// Shapes, support functions, k-DOPs
// ---------------------------------------------------------------------------

/// A brush shape, as one type so a log can hold a mixture.
///
/// The three the crate has. Identical in kind to `brush::tests::Shape`, which is
/// private to that module.
#[derive(Clone, Copy)]
enum Shape {
    Sphere(Sphere<f64>),
    Cube(BoxExact<f64>),
    Capsule(Capsule<f64>),
}

impl Sdf for Shape {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Sphere(s) => s.sample(p),
            Self::Cube(b) => b.sample(p),
            Self::Capsule(c) => c.sample(p),
        }
    }
}

impl Shape {
    /// `h(d) = max_{x ∈ S} d·x`, exact for all three shapes.
    ///
    /// Klosowski et al. 1998: the `k`-DOP is `⋂_i { x : d_i·x ≤ h(d_i) }`, the
    /// tightest convex polytope containing `S` in that direction set.
    fn support(&self, d: [f64; 3]) -> f64 {
        let dl = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        match self {
            Self::Sphere(s) => {
                s.center[0] * d[0] + s.center[1] * d[1] + s.center[2] * d[2] + s.radius * dl
            }
            Self::Cube(b) => {
                b.center[0] * d[0]
                    + b.center[1] * d[1]
                    + b.center[2] * d[2]
                    + b.half_extents[0] * d[0].abs()
                    + b.half_extents[1] * d[1].abs()
                    + b.half_extents[2] * d[2].abs()
            }
            Self::Capsule(c) => {
                let a = c.a[0] * d[0] + c.a[1] * d[1] + c.a[2] * d[2];
                let b = c.b[0] * d[0] + c.b[1] * d[1] + c.b[2] * d[2];
                a.max(b) + c.radius * dl
            }
        }
    }

    /// Axis-aligned bounds of the shape, for cell overlap rejection.
    fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let mut lo = [0.0f64; 3];
        let mut hi = [0.0f64; 3];
        for a in 0..3 {
            let mut d = [0.0f64; 3];
            d[a] = 1.0;
            hi[a] = self.support(d);
            d[a] = -1.0;
            lo[a] = -self.support(d);
        }
        (lo, hi)
    }
}

/// The 6-DOP directions: the coordinate axes. Gives the exact AABB, and for a
/// box brush the exact brush.
const DIRS6: [[f64; 3]; 6] = [
    [1.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, -1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, -1.0],
];

/// The eight cube-diagonal directions, which with [`DIRS6`] make a 14-DOP.
const DIRS8: [[f64; 3]; 8] = [
    [1.0, 1.0, 1.0],
    [1.0, 1.0, -1.0],
    [1.0, -1.0, 1.0],
    [1.0, -1.0, -1.0],
    [-1.0, 1.0, 1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [-1.0, -1.0, -1.0],
];

/// The direction set for a `k`-DOP with `k` planes, normalised.
fn dirs(k: usize) -> Vec<[f64; 3]> {
    let mut out: Vec<[f64; 3]> = DIRS6.to_vec();
    if k == 14 {
        let s = 1.0 / 3.0f64.sqrt();
        out.extend(DIRS8.iter().map(|d| [d[0] * s, d[1] * s, d[2] * s]));
    } else {
        assert_eq!(k, 6, "only the 6-DOP and the 14-DOP are fixtures here");
    }
    out
}

/// The `k`-DOP's supporting half-spaces for one shape.
fn dop_planes(shape: &Shape, ds: &[[f64; 3]]) -> Vec<Plane> {
    ds.iter()
        .map(|d| Plane {
            n: *d,
            d: -shape.support(*d),
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Base fields
// ---------------------------------------------------------------------------

/// The base field a chunk is cut out of.
#[derive(Clone, Copy)]
enum Base {
    /// `M-38`'s `BoxExact::canonical()`, for `M-36`'s fixture.
    Box(BoxExact<f64>),
    Fbm(FbmTerrain<f64>),
    Gyroid(Gyroid<f64>),
}

impl Sdf for Base {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Box(b) => b.sample(p),
            Self::Fbm(f) => f.sample(p),
            Self::Gyroid(g) => g.sample(p),
        }
    }

    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        match self {
            Self::Box(b) => b.gradient(p),
            Self::Fbm(f) => f.gradient(p),
            Self::Gyroid(g) => g.gradient(p),
        }
    }
}

impl Base {
    const fn name(&self) -> &'static str {
        match self {
            Self::Box(_) => "box_exact",
            Self::Fbm(_) => "fbm_terrain",
            Self::Gyroid(_) => "gyroid",
        }
    }
}

/// One edit: a shape, an op, and a **stable id**.
///
/// The id is the edit's own identity, assigned by whoever made the edit. It is
/// not its position in the log, and that distinction is the whole of C2: `M-36`
/// is about eight identifiable edits arriving in any order.
#[derive(Clone, Copy)]
struct Edit {
    brush: Brush<Shape>,
    id: u16,
}

/// Fold the base and a log, in the log's order, through the crate's own `apply`.
#[inline]
fn fold(base: &Base, log: &[Edit], order: &[usize], p: [f64; 3]) -> f64 {
    let mut v = base.sample(p);
    for &i in order {
        let e = &log[i];
        v = apply(e.brush.op, v, e.brush.shape.sample(p));
    }
    v
}

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

/// `M-36`'s eight brushes, all three shapes represented, deliberately
/// overlapping so the order could matter — as **`Subtract`**, which is one of
/// the two same-kind cases `M-36` reports one distinct result for.
///
/// The centres, the radii, the half-extents and the `i % 3` shape rotation are
/// copied from `brush::tests::eight` and *are* the fixture; changing any of them
/// would stop this from being `M-36`'s measurement.
fn eight() -> Vec<Edit> {
    let centres = [
        [0.30, 0.10, -0.20],
        [-0.25, 0.35, 0.15],
        [0.05, -0.30, 0.25],
        [-0.15, -0.10, -0.35],
        [0.40, 0.25, 0.05],
        [-0.35, 0.05, 0.30],
        [0.20, -0.40, -0.10],
        [-0.05, 0.20, -0.30],
    ];
    centres
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let shape = match i % 3 {
                0 => Shape::Sphere(Sphere {
                    center: *c,
                    radius: 0.30 + 0.02 * i as f64,
                }),
                1 => Shape::Cube(BoxExact {
                    center: *c,
                    half_extents: [0.22 + 0.01 * i as f64; 3],
                }),
                _ => Shape::Capsule(Capsule {
                    a: *c,
                    b: [c[0] + 0.25, c[1] - 0.15, c[2] + 0.1],
                    radius: 0.16,
                }),
            };
            Edit {
                brush: Brush {
                    shape,
                    op: BrushOp::Subtract,
                },
                id: i as u16,
            }
        })
        .collect()
}

/// A chunk: a cubical region, a cell count, and the cell size that follows.
#[derive(Clone, Copy)]
struct Chunk {
    origin: [f64; 3],
    extent: f64,
    cells: u32,
}

impl Chunk {
    fn cell_size(&self) -> f64 {
        self.extent / f64::from(self.cells)
    }

    fn cell_lo(&self, i: u32, j: u32, k: u32) -> [f64; 3] {
        let h = self.cell_size();
        [
            self.origin[0] + f64::from(i) * h,
            self.origin[1] + f64::from(j) * h,
            self.origin[2] + f64::from(k) * h,
        ]
    }
}

/// `n` digging brushes, hugging the base surface, deterministic.
///
/// `(x, z)` walk a golden-angle lattice over the middle 70% of the chunk — the
/// low-discrepancy sequence, so the brushes spread rather than clump — and `y`
/// is the first zero crossing found from the top of the chunk by a 64-step scan
/// and 40 bisections. A candidate column with no crossing is skipped, which is
/// a filter rather than a fallback: the count is asserted at the end.
fn dig_log(base: &Base, chunk: Chunk, n: usize, radius: f64) -> Vec<Edit> {
    const PHI: f64 = 1.618_033_988_749_895;
    let mut out = Vec::with_capacity(n);
    let mut cand = 0usize;
    while out.len() < n {
        cand += 1;
        assert!(cand < 100_000, "no surface found in the chunk to dig into");
        let u = (cand as f64 * PHI).fract();
        let v = (cand as f64 * PHI * PHI).fract();
        let x = chunk.origin[0] + chunk.extent * (0.15 + 0.7 * u);
        let z = chunk.origin[2] + chunk.extent * (0.15 + 0.7 * v);
        let y0 = chunk.origin[1] + chunk.extent * 0.9;
        let y1 = chunk.origin[1] + chunk.extent * 0.1;
        let mut lo = y0;
        let mut f_lo = base.sample([x, lo, z]);
        let mut hit = None;
        for s in 1..=64 {
            let y = y0 + (y1 - y0) * f64::from(s) / 64.0;
            let f = base.sample([x, y, z]);
            if (f_lo > 0.0) != (f > 0.0) {
                hit = Some((lo, y));
                break;
            }
            lo = y;
            f_lo = f;
        }
        let Some((mut a, mut b)) = hit else { continue };
        for _ in 0..40 {
            let m = 0.5 * (a + b);
            if (base.sample([x, m, z]) > 0.0) == (base.sample([x, a, z]) > 0.0) {
                a = m;
            } else {
                b = m;
            }
        }
        let y = 0.5 * (a + b);
        let i = out.len();
        let c = [x, y, z];
        let shape = match i % 3 {
            0 => Shape::Sphere(Sphere { center: c, radius }),
            1 => Shape::Cube(BoxExact {
                center: c,
                half_extents: [radius * 0.85; 3],
            }),
            _ => Shape::Capsule(Capsule {
                a: c,
                b: [c[0] + radius, c[1] - radius * 0.5, c[2] + radius * 0.3],
                radius: radius * 0.6,
            }),
        };
        out.push(Edit {
            brush: Brush {
                shape,
                op: BrushOp::Subtract,
            },
            id: i as u16,
        });
    }
    out
}

// ---------------------------------------------------------------------------
// The plane table, and one cell's arrangement
// ---------------------------------------------------------------------------

/// Reserved ids: six cube faces then the base-surface plane.
const RESERVED: usize = 7;

/// The surface plane's id.
const SURFACE_ID: u16 = 6;

/// The whole chunk's plane table.
///
/// Ids `0..6` are the current cell's cube faces and `6` its surface plane, both
/// rewritten per cell; from `RESERVED` on, one block of `k` planes per edit,
/// indexed by the edit's **stable id**. That is why reordering the log does not
/// renumber a plane, and it is what the `solve` arm's determinism rests on.
struct Table {
    planes: Vec<Plane>,
    k: usize,
}

impl Table {
    fn build(log: &[Edit], ds: &[[f64; 3]]) -> Self {
        let k = ds.len();
        let n = log.iter().map(|e| usize::from(e.id)).max().unwrap_or(0) + 1;
        let mut planes = vec![
            Plane {
                n: [0.0, 0.0, 0.0],
                d: 0.0,
            };
            RESERVED + n * k
        ];
        for e in log {
            let base = RESERVED + usize::from(e.id) * k;
            for (j, p) in dop_planes(&e.brush.shape, ds).into_iter().enumerate() {
                planes[base + j] = p;
            }
        }
        Self { planes, k }
    }

    fn brush_plane_ids(&self, id: u16) -> std::ops::Range<u16> {
        let b = (RESERVED + usize::from(id) * self.k) as u16;
        b..b + self.k as u16
    }

    /// Rewrite ids `0..6` for one cell's cube.
    fn set_cube(&mut self, lo: [f64; 3], h: f64) {
        for a in 0..3 {
            let mut nm = [0.0f64; 3];
            nm[a] = -1.0;
            self.planes[a * 2] = Plane { n: nm, d: lo[a] };
            let mut np = [0.0f64; 3];
            np[a] = 1.0;
            self.planes[a * 2 + 1] = Plane {
                n: np,
                d: -(lo[a] + h),
            };
        }
    }
}

/// The cube's eight corners, tagged with their three face planes.
fn cube_verts(lo: [f64; 3], h: f64) -> Vec<Vert> {
    let mut out = Vec::with_capacity(8);
    for k in 0..2u32 {
        for j in 0..2u32 {
            for i in 0..2u32 {
                // Axis `a`'s two faces are ids `2a` and `2a+1`, which is what
                // `set_cube` writes. Tagging them `i*2`, `2+j*2`, `4+k*2`
                // aliased the x-max face onto the y-min face, the y-max onto
                // the z-min, and the z-max onto `SURFACE_ID` — every corner on
                // a far face then carried a repeated plane id, and a repeated
                // id makes `shared_two` report a false edge.
                let t = [i as u16, (2 + j) as u16, (4 + k) as u16];
                out.push(Vert {
                    p: [
                        lo[0] + f64::from(i) * h,
                        lo[1] + f64::from(j) * h,
                        lo[2] + f64::from(k) * h,
                    ],
                    t,
                });
            }
        }
    }
    out
}

/// Buffers one worker reuses across every cell it handles.
struct Scratch {
    sd: Vec<f64>,
    cur: Vec<Vec<Vert>>,
    next: Vec<Vec<Vert>>,
    pool: Vec<Vec<Vert>>,
    face: Vec<(f64, [f64; 3])>,
    cuts: Vec<u16>,
    hashes: Vec<u64>,
    htopo: Vec<u64>,
    hpos: Vec<u64>,
    hsolid: Vec<u64>,
}

impl Scratch {
    fn new() -> Self {
        Self {
            sd: Vec::with_capacity(64),
            cur: Vec::with_capacity(64),
            next: Vec::with_capacity(64),
            pool: Vec::with_capacity(64),
            face: Vec::with_capacity(32),
            cuts: Vec::with_capacity(64),
            hashes: Vec::with_capacity(64),
            htopo: Vec::with_capacity(64),
            hpos: Vec::with_capacity(64),
            hsolid: Vec::with_capacity(64),
        }
    }

    fn take(&mut self) -> Vec<Vert> {
        self.pool.pop().map_or_else(
            || Vec::with_capacity(32),
            |mut v| {
                v.clear();
                v
            },
        )
    }
}

/// The arrangement of a seed convex piece by a list of planes, breadth-first.
///
/// Every leaf is `seed ∩ ⋂_i H_i^{s_i}` for some sign vector, hence convex. The
/// leaf **set** is a property of the plane set; the order the planes are applied
/// in affects only the arithmetic path each vertex took, which is exactly what
/// C2's two modes separate.
fn arrange(
    table: &Table,
    seed: Vec<Vert>,
    mode: Mode,
    eps: f64,
    sc: &mut Scratch,
    flags: &mut Flags,
) {
    for v in sc.cur.drain(..) {
        sc.pool.push(v);
    }
    sc.cur.push(seed);
    flags.max_cuts = flags.max_cuts.max(sc.cuts.len());
    for ci in 0..sc.cuts.len() {
        let id = sc.cuts[ci];
        sc.next.clear();
        while let Some(piece) = sc.cur.pop() {
            let mut inside = sc.take();
            let mut outside = sc.take();
            split(
                &table.planes,
                &piece,
                id,
                mode,
                eps,
                &mut sc.sd,
                &mut inside,
                &mut outside,
                flags,
            );
            sc.pool.push(piece);
            for part in [inside, outside] {
                if part.len() >= 4 {
                    sc.next.push(part);
                } else {
                    sc.pool.push(part);
                }
            }
            assert!(
                sc.next.len() <= LEAF_CAP,
                "one cell's arrangement passed {LEAF_CAP} leaves at plane {} of {}, which is not \
                 an arrangement any more: the split is manufacturing degenerate parts",
                ci + 1,
                sc.cuts.len(),
            );
        }
        flags.max_leaves = flags.max_leaves.max(sc.next.len());
        std::mem::swap(&mut sc.cur, &mut sc.next);
    }
}

/// What one chunk partition measured.
#[derive(Clone, Default)]
struct Partition {
    solid_pieces: usize,
    total_pieces: usize,
    cells_full: usize,
    cells_boundary: usize,
    cells_cut_by_brushes: usize,
    max_pieces_in_a_cell: usize,
    compound_volume: f64,
    flags: Flags,
    /// Xor-folded over hashed cells: the canonical hash, the emission-order
    /// hash, and the three projections that say *which part* of a partition two
    /// orderings disagreed about.
    ///
    /// `hash` is everything. `hash_topo` is the leaf plane-triple sets alone —
    /// no coordinates, no solid flags — so it answers "is the *arrangement*
    /// order-free?". `hash_pos` adds the vertex bit patterns, `hash_solid`
    /// replaces them with the solid classification. A `distinct` above one whose
    /// `distinct_topology` is one is a floating-point finding; one whose
    /// topology also splits is a combinatorial one, and those are different
    /// bugs. Without the split, `distinct_partitions = 6` names no mechanism.
    hash: u64,
    hash_raw: u64,
    hash_topo: u64,
    hash_pos: u64,
    hash_solid: u64,
    cells_hashed: usize,
}

impl Partition {
    fn merge(&mut self, o: &Self) {
        self.solid_pieces += o.solid_pieces;
        self.total_pieces += o.total_pieces;
        self.cells_full += o.cells_full;
        self.cells_boundary += o.cells_boundary;
        self.cells_cut_by_brushes += o.cells_cut_by_brushes;
        self.max_pieces_in_a_cell = self.max_pieces_in_a_cell.max(o.max_pieces_in_a_cell);
        self.compound_volume += o.compound_volume;
        self.flags.merge(o.flags);
        self.cells_hashed += o.cells_hashed;
        // Commutative fold: workers own disjoint cell ranges and the per-cell
        // hashes already carry the cell index, so xor-and-add is order-free and
        // still position-sensitive.
        self.hash ^= o.hash;
        self.hash_raw ^= o.hash_raw;
        self.hash_topo ^= o.hash_topo;
        self.hash_pos ^= o.hash_pos;
        self.hash_solid ^= o.hash_solid;
    }
}

/// FNV-1a, one `u64` at a time.
#[inline]
fn fnv(h: &mut u64, x: u64) {
    for b in x.to_le_bytes() {
        *h ^= u64::from(b);
        *h = h.wrapping_mul(0x0100_0000_01b3);
    }
}

/// Fold one cell's per-piece hashes into a cell hash, order-free.
///
/// Sorted, so the value is the piece **multiset**: which order the breadth-first
/// split emitted the leaves in is the question `hash_raw` asks, and it must not
/// leak into the questions the other four ask.
fn cell_fold(lin: u32, hs: &mut [u64]) -> u64 {
    hs.sort_unstable();
    let mut ch = 0u64;
    fnv(&mut ch, u64::from(lin));
    fnv(&mut ch, hs.len() as u64);
    for &x in hs.iter() {
        fnv(&mut ch, x);
    }
    ch
}

/// Build a chunk's convex-cell partition, over the cells `range`.
///
/// `hash_cells` selects the cells that go into the hash: `None` hashes every
/// cell, `Some(set)` only those. C2's arm hashes only the cells at least one
/// brush plane crosses — the rest receive no brush plane at all, so their piece
/// set is literally the same objects for every ordering, and hashing them would
/// pad the measurement with a constant.
fn partition(
    base: &Base,
    log: &[Edit],
    order: &[usize],
    chunk: Chunk,
    table: &mut Table,
    mode: Mode,
    planes: Planes,
    range: std::ops::Range<u32>,
    only_cut_cells: bool,
    sc: &mut Scratch,
) -> Partition {
    let h = chunk.cell_size();
    let cells = chunk.cells;
    let cube_vol = h * h * h;
    let eps = GEOM_EPS * h;
    let aabbs: Vec<([f64; 3], [f64; 3])> = log.iter().map(|e| e.brush.shape.aabb()).collect();
    let mut out = Partition::default();
    for lin in range {
        let i = lin % cells;
        let j = (lin / cells) % cells;
        let k = lin / (cells * cells);
        let lo = chunk.cell_lo(i, j, k);
        let hi = [lo[0] + h, lo[1] + h, lo[2] + h];
        table.set_cube(lo, h);
        let seed = cube_verts(lo, h);

        let (centre, cut_by_brushes) = cell_cuts(
            base,
            log,
            order,
            &aabbs,
            lo,
            hi,
            h,
            &seed,
            table,
            planes,
            &mut sc.cuts,
        );
        if cut_by_brushes {
            out.cells_cut_by_brushes += 1;
        }
        if only_cut_cells && !cut_by_brushes {
            continue;
        }

        if sc.cuts.is_empty() {
            // No plane crosses: the cube is wholly solid or wholly empty.
            let solid = fold(base, log, order, centre) <= 0.0;
            if solid {
                out.cells_full += 1;
                out.solid_pieces += 1;
                out.total_pieces += 1;
                out.compound_volume += cube_vol;
                out.max_pieces_in_a_cell = out.max_pieces_in_a_cell.max(1);
                let mut topo = 0u64;
                for v in &seed {
                    fnv(&mut topo, u64::from(v.t[0]));
                    fnv(&mut topo, u64::from(v.t[1]));
                    fnv(&mut topo, u64::from(v.t[2]));
                }
                let mut pos = topo;
                for v in &seed {
                    for a in 0..3 {
                        fnv(&mut pos, v.p[a].to_bits());
                    }
                }
                let mut sol = topo;
                fnv(&mut sol, 1);
                let mut ph = pos;
                fnv(&mut ph, 1);
                out.hash ^= cell_fold(lin, &mut [ph]);
                out.hash_raw ^= cell_fold(lin, &mut [ph]);
                out.hash_topo ^= cell_fold(lin, &mut [topo]);
                out.hash_pos ^= cell_fold(lin, &mut [pos]);
                out.hash_solid ^= cell_fold(lin, &mut [sol]);
                out.cells_hashed += 1;
            }
            continue;
        }
        out.cells_boundary += 1;

        arrange(table, seed, mode, eps, sc, &mut out.flags);

        // Volume closure: the leaves must tile the cube.
        let mut leaf_vol = 0.0f64;
        let mut kept = 0usize;
        sc.hashes.clear();
        sc.htopo.clear();
        sc.hpos.clear();
        sc.hsolid.clear();
        let mut raw = 0u64;
        fnv(&mut raw, u64::from(lin));
        for piece in &mut sc.cur {
            // Canonical vertex order **before** anything is summed over it. A
            // vertex of a convex polytope *is* its plane triple — `tidy` dedups
            // by triple, so the triple is a unique key within a piece — and
            // sorting by it makes every quantity below a function of the leaf's
            // vertex *set* rather than of the order the breadth-first split
            // happened to emit it in.
            //
            // The first version of this file sorted only for the hash and left
            // `centroid` and `volume` summing in emission order. Both are naive
            // float sums, so both were ULP-order-dependent, and both feed a
            // decision: the centroid picks the leaf's solid flag through `fold`,
            // the volume picks whether the leaf is a sliver at all. That is two
            // order dependencies manufactured by the instrument, sitting in
            // front of the exact question C2 asks.
            piece.sort_unstable_by(|a, b| {
                a.t.cmp(&b.t)
                    .then_with(|| a.p[0].total_cmp(&b.p[0]))
                    .then_with(|| a.p[1].total_cmp(&b.p[1]))
                    .then_with(|| a.p[2].total_cmp(&b.p[2]))
            });
            let vol = volume(&table.planes, piece, &mut sc.face);
            leaf_vol += vol;
            if vol < SLIVER_FRACTION * cube_vol {
                out.flags.slivers += 1;
                piece.clear();
                continue;
            }
            kept += 1;
            let g = centroid(piece);
            let solid = fold(base, log, order, g) <= 0.0;
            out.flags.convexity_violations += convexity_violations(&table.planes, piece, eps);
            if solid {
                out.solid_pieces += 1;
                out.compound_volume += vol;
            }
            // The three projections, nested so each adds exactly one kind of
            // information to the one before it.
            let mut topo = 0u64;
            for v in piece.iter() {
                fnv(&mut topo, u64::from(v.t[0]));
                fnv(&mut topo, u64::from(v.t[1]));
                fnv(&mut topo, u64::from(v.t[2]));
            }
            let mut pos = topo;
            for v in piece.iter() {
                for a in 0..3 {
                    fnv(&mut pos, v.p[a].to_bits());
                }
            }
            let mut sol = topo;
            fnv(&mut sol, u64::from(u8::from(solid)));
            let mut ph = pos;
            fnv(&mut ph, u64::from(u8::from(solid)));
            sc.hashes.push(ph);
            sc.htopo.push(topo);
            sc.hpos.push(pos);
            sc.hsolid.push(sol);
            fnv(&mut raw, ph);
        }
        out.total_pieces += kept;
        out.max_pieces_in_a_cell = out.max_pieces_in_a_cell.max(kept);
        let err = (leaf_vol - cube_vol).abs() / cube_vol;
        // The paper's *non-overlapping* half, checked per cell rather than in
        // aggregate so the failure names the cell and its plane set. Leaves
        // that overlap sum high; a leaf that was never generated sums low.
        if err >= VOLUME_CLOSURE_REL {
            let mut s = String::new();
            for piece in &sc.cur {
                use std::fmt::Write as _;
                let vv = volume(&table.planes, piece, &mut sc.face);
                write!(s, "({},{vv:.3e})", piece.len()).unwrap();
            }
            panic!(
                "cell {lin}'s leaves did not tile its cube: cuts={:?} leaf_vol={leaf_vol:e} \
                 cube={cube_vol:e} pieces={s}",
                sc.cuts
            );
        }
        out.flags.worst_volume_error = out.flags.worst_volume_error.max(err);
        out.hash ^= cell_fold(lin, &mut sc.hashes);
        out.hash_topo ^= cell_fold(lin, &mut sc.htopo);
        out.hash_pos ^= cell_fold(lin, &mut sc.hpos);
        out.hash_solid ^= cell_fold(lin, &mut sc.hsolid);
        out.hash_raw ^= raw;
        out.cells_hashed += 1;
    }
    out
}

/// Fill `cuts` with the plane ids that cross one cell, and write that cell's
/// base-surface plane into the table.
///
/// Returns the cell centre and whether any brush plane was added. Shared by the
/// partition builder and the compound collector, so the two cannot disagree
/// about what a cell's cut set is.
fn cell_cuts(
    base: &Base,
    log: &[Edit],
    order: &[usize],
    aabbs: &[([f64; 3], [f64; 3])],
    lo: [f64; 3],
    hi: [f64; 3],
    h: f64,
    seed: &[Vert],
    table: &mut Table,
    planes: Planes,
    cuts: &mut Vec<u16>,
) -> ([f64; 3], bool) {
    cuts.clear();
    // The base surface's plane, if the cube straddles it.
    let mut any_neg = false;
    let mut any_pos = false;
    for v in seed {
        if base.sample(v.p) <= 0.0 {
            any_neg = true;
        } else {
            any_pos = true;
        }
    }
    let centre = [lo[0] + h * 0.5, lo[1] + h * 0.5, lo[2] + h * 0.5];
    if any_neg && any_pos {
        let g = base.gradient(centre);
        let gl = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        if gl > 1e-12 {
            let nh = [g[0] / gl, g[1] / gl, g[2] / gl];
            let v = base.sample(centre);
            // Solid is f <= 0 and ∇f points away from it, so the half-space is
            // n̂·x + d <= 0 with the plane through centre − (v/|∇f|)·n̂.
            let foot = [
                centre[0] - v / gl * nh[0],
                centre[1] - v / gl * nh[1],
                centre[2] - v / gl * nh[2],
            ];
            table.planes[SURFACE_ID as usize] = Plane {
                n: nh,
                d: -(nh[0] * foot[0] + nh[1] * foot[1] + nh[2] * foot[2]),
            };
            cuts.push(SURFACE_ID);
        }
    }
    let surface_cuts = cuts.len();

    // Brush planes, in the log's order, that actually cross this cube.
    for &oi in order {
        let (blo, bhi) = aabbs[oi];
        if bhi[0] <= lo[0]
            || blo[0] >= hi[0]
            || bhi[1] <= lo[1]
            || blo[1] >= hi[1]
            || bhi[2] <= lo[2]
            || blo[2] >= hi[2]
        {
            continue;
        }
        for id in table.brush_plane_ids(log[oi].id) {
            let pl = &table.planes[id as usize];
            let mut lo_v = f64::INFINITY;
            let mut hi_v = f64::NEG_INFINITY;
            for v in seed {
                let s = pl.at(v.p);
                lo_v = lo_v.min(s);
                hi_v = hi_v.max(s);
            }
            if lo_v < 0.0 && hi_v > 0.0 {
                cuts.push(id);
            }
        }
    }
    if planes == Planes::Canonical {
        // A half-space appearing twice adds nothing to an intersection of
        // half-spaces, and a half-space appearing beside its own complement adds
        // only a zero-width leaf — but *which* of the two ids a vertex on that
        // plane gets tagged with is decided by the traversal, and the
        // measurement says that is the whole of the `solve` arm's residual.
        //
        // Keeping the smallest **stable** id of each coincident group is
        // order-free because the ids are. It deliberately does **not** reorder
        // the survivors: `arrange` still applies them in the log's order, so the
        // sweep is still asking whether the arrangement depends on that order
        // rather than being handed the answer.
        let eps = GEOM_EPS * h;
        let n = cuts.len();
        let mut w = 0usize;
        for i in 0..n {
            let x = cuts[i];
            let px = table.planes[x as usize];
            let mut redundant = false;
            for j in 0..n {
                let y = cuts[j];
                if y < x && coincident(&px, &table.planes[y as usize], eps) {
                    redundant = true;
                    break;
                }
            }
            if !redundant {
                cuts[w] = x;
                w += 1;
            }
        }
        cuts.truncate(w);
    }
    (centre, cuts.len() > surface_cuts)
}

/// Partition a whole chunk, over `threads` workers on disjoint cell ranges.
fn partition_all(
    base: &Base,
    log: &[Edit],
    order: &[usize],
    chunk: Chunk,
    ds: &[[f64; 3]],
    mode: Mode,
    planes: Planes,
    only_cut_cells: bool,
    threads: usize,
) -> Partition {
    let total = chunk.cells * chunk.cells * chunk.cells;
    let per = total.div_ceil(threads as u32);
    let parts: Vec<Partition> = std::thread::scope(|s| {
        let mut hs = Vec::with_capacity(threads);
        for t in 0..threads as u32 {
            let lo = t * per;
            let hi = ((t + 1) * per).min(total);
            hs.push(s.spawn(move || {
                let mut table = Table::build(log, ds);
                let mut sc = Scratch::new();
                partition(
                    base,
                    log,
                    order,
                    chunk,
                    &mut table,
                    mode,
                    planes,
                    lo..hi,
                    only_cut_cells,
                    &mut sc,
                )
            }));
        }
        hs.into_iter().map(|h| h.join().expect("worker")).collect()
    });
    let mut out = Partition::default();
    for p in &parts {
        out.merge(p);
    }
    out
}

// ---------------------------------------------------------------------------
// C1: the fracture
// ---------------------------------------------------------------------------

/// The solid pieces of a compound, grouped by the cell they came from.
struct Compound {
    cells: Vec<Vec<Vec<Vert>>>,
    lo: Vec<[f64; 3]>,
    /// Each cell's base-surface plane. Ids `0..=6` in the table are **rewritten
    /// per cell**, so a piece tagged with `SURFACE_ID` is meaningless against
    /// another cell's surface plane: the fracture restores this before it
    /// touches the cell's pieces. Without it every such piece is pruned to
    /// nothing and the fracture reports zero fragments.
    surf: Vec<Plane>,
}

/// Collect the compound's solid pieces for the cells the impact overlaps.
///
/// This is the **offline** half of Müller's method, and it is deliberately
/// outside C1's timed region: the compound is precomputed, which is the whole
/// premise. It goes through the same [`cell_cuts`] the partition builder uses,
/// so the pieces C1 fractures are the pieces C3 counted.
fn compound_near(
    base: &Base,
    log: &[Edit],
    chunk: Chunk,
    ds: &[[f64; 3]],
    impact: &Shape,
    mode: Mode,
    planes: Planes,
) -> (Compound, Table) {
    let order: Vec<usize> = (0..log.len()).collect();
    let aabbs: Vec<([f64; 3], [f64; 3])> = log.iter().map(|e| e.brush.shape.aabb()).collect();
    let mut table = Table::build(log, ds);
    let mut sc = Scratch::new();
    let mut flags = Flags::default();
    let h = chunk.cell_size();
    let cube_vol = h * h * h;
    let (ilo, ihi) = impact.aabb();
    let mut cells = Vec::new();
    let mut los = Vec::new();
    let mut surfs = Vec::new();
    let n = chunk.cells;
    for lin in 0..n * n * n {
        let i = lin % n;
        let j = (lin / n) % n;
        let k = lin / (n * n);
        let lo = chunk.cell_lo(i, j, k);
        let hi = [lo[0] + h, lo[1] + h, lo[2] + h];
        if ihi[0] <= lo[0]
            || ilo[0] >= hi[0]
            || ihi[1] <= lo[1]
            || ilo[1] >= hi[1]
            || ihi[2] <= lo[2]
            || ilo[2] >= hi[2]
        {
            continue;
        }
        table.set_cube(lo, h);
        let seed = cube_verts(lo, h);
        let (centre, _) = cell_cuts(
            base,
            log,
            &order,
            &aabbs,
            lo,
            hi,
            h,
            &seed,
            &mut table,
            planes,
            &mut sc.cuts,
        );
        let mut keep: Vec<Vec<Vert>> = Vec::new();
        if sc.cuts.is_empty() {
            if fold(base, log, &order, centre) <= 0.0 {
                keep.push(seed);
            }
        } else {
            arrange(&table, seed, mode, GEOM_EPS * h, &mut sc, &mut flags);
            for piece in &sc.cur {
                if volume(&table.planes, piece, &mut sc.face) < SLIVER_FRACTION * cube_vol {
                    continue;
                }
                if fold(base, log, &order, centroid(piece)) <= 0.0 {
                    keep.push(piece.clone());
                }
            }
        }
        if !keep.is_empty() {
            cells.push(keep);
            los.push(lo);
            surfs.push(table.planes[SURFACE_ID as usize]);
        }
    }
    (
        Compound {
            cells,
            lo: los,
            surf: surfs,
        },
        table,
    )
}

/// What one timed fracture produced.
struct Fracture {
    fragments: usize,
    remainder: usize,
    total_ms: f64,
    worst_per_fragment_ms: f64,
    cells: usize,
}

/// Müller's runtime step: intersect the impact's convex pattern against the
/// precomputed compound.
///
/// The pattern is the impact brush's `k`-DOP, a convex polytope. Every solid
/// piece of every overlapped cell is split by the pattern's planes; a leaf
/// inside every pattern plane is a **fragment**, and it is already convex. No
/// field is sampled and nothing is decomposed — that is the claim.
fn fracture(
    comp: &Compound,
    table: &mut Table,
    impact_planes: &[Plane],
    impact_base_id: u16,
    chunk: Chunk,
    mode: Mode,
    sc: &mut Scratch,
) -> Fracture {
    let h = chunk.cell_size();
    let cube_vol = h * h * h;
    for (j, p) in impact_planes.iter().enumerate() {
        table.planes[usize::from(impact_base_id) + j] = *p;
    }
    let mut fragments = 0usize;
    let mut remainder = 0usize;
    let mut total = 0.0f64;
    let mut worst = 0.0f64;
    let mut flags = Flags::default();
    for (ci, pieces) in comp.cells.iter().enumerate() {
        table.set_cube(comp.lo[ci], h);
        table.planes[SURFACE_ID as usize] = comp.surf[ci];
        let t = Instant::now();
        let mut cell_frags = 0usize;
        let mut cell_rem = 0usize;
        for piece in pieces {
            sc.cuts.clear();
            for j in 0..impact_planes.len() {
                sc.cuts.push(impact_base_id + j as u16);
            }
            arrange(table, piece.clone(), mode, GEOM_EPS * h, sc, &mut flags);
            for leaf in &sc.cur {
                let vol = volume(&table.planes, leaf, &mut sc.face);
                if vol < SLIVER_FRACTION * cube_vol {
                    continue;
                }
                let g = centroid(leaf);
                let inside = impact_planes.iter().all(|p| p.at(g) <= 0.0);
                if inside {
                    cell_frags += 1;
                } else {
                    cell_rem += 1;
                }
            }
        }
        let ms = t.elapsed().as_secs_f64() * 1e3;
        total += ms;
        if cell_frags > 0 {
            worst = worst.max(ms / cell_frags as f64);
        }
        fragments += cell_frags;
        remainder += cell_rem;
    }
    Fracture {
        fragments,
        remainder,
        total_ms: total,
        worst_per_fragment_ms: worst,
        cells: comp.cells.len(),
    }
}

// ---------------------------------------------------------------------------
// C2: all 40,320 orderings
// ---------------------------------------------------------------------------

/// Every permutation of `0..BRUSHES`, by Heap's algorithm.
///
/// Generated rather than sampled, exactly as `M-36` does it. The registered
/// vacuity control is a count against [`ORDERINGS`], asserted below.
fn permutations() -> Vec<[usize; BRUSHES]> {
    let mut out = Vec::with_capacity(ORDERINGS);
    let mut items = [0usize; BRUSHES];
    for (i, slot) in items.iter_mut().enumerate() {
        *slot = i;
    }
    let mut counters = [0usize; BRUSHES];
    out.push(items);
    let mut i = 0;
    while i < BRUSHES {
        if counters[i] < i {
            if i % 2 == 0 {
                items.swap(0, i);
            } else {
                items.swap(counters[i], i);
            }
            out.push(items);
            counters[i] += 1;
            i = 0;
        } else {
            counters[i] = 0;
            i += 1;
        }
    }
    out
}

/// Where the vertex-tag clipper's adjacency rule is unsound, measured.
///
/// The rule — two vertices of a convex piece are adjacent exactly when they
/// share two planes — is sound for a **simple** polytope, one whose every vertex
/// lies on exactly three facets. Four planes through one point breaks it: three
/// tagged vertices then share a plane pair, two of them sit at the same
/// position, and which leaf survives can depend on the order the splits were
/// applied in. That is defect 4's mechanism, and the 6-DOP arm is not immune to
/// it just because it is axis-aligned — a brush plane landing exactly on a cell
/// corner puts four planes through that corner.
///
/// So the fixture is asked rather than assumed. Order-independent by
/// construction: it reads the cell's cut **set** and the cube's six faces, never
/// a traversal, so it runs once per arm outside the 40,320-ordering sweep.
#[derive(Default)]
struct Simplicity {
    /// Pairs among a cell's planes with the same point set, either orientation.
    /// A duplicated plane makes zero-width leaves and every vertex on it
    /// non-simple.
    coincident_pairs: usize,
    /// Distinct points in or on a cell's cube where four or more of its planes
    /// meet.
    nonsimple_points: usize,
    /// The most planes found through any one point.
    max_planes_at_a_point: usize,
    /// Cells holding at least one such point.
    nonsimple_cells: usize,
}

fn simplicity_scan(
    base: &Base,
    log: &[Edit],
    chunk: Chunk,
    ds: &[[f64; 3]],
    planes: Planes,
) -> Simplicity {
    let h = chunk.cell_size();
    let eps = GEOM_EPS * h;
    let cells = chunk.cells;
    let order: Vec<usize> = (0..log.len()).collect();
    let aabbs: Vec<([f64; 3], [f64; 3])> = log.iter().map(|e| e.brush.shape.aabb()).collect();
    let mut table = Table::build(log, ds);
    let mut cuts: Vec<u16> = Vec::new();
    let mut all: Vec<u16> = Vec::new();
    let mut pts: Vec<[f64; 3]> = Vec::new();
    let mut out = Simplicity::default();
    for lin in 0..cells * cells * cells {
        let i = lin % cells;
        let j = (lin / cells) % cells;
        let k = lin / (cells * cells);
        let lo = chunk.cell_lo(i, j, k);
        let hi = [lo[0] + h, lo[1] + h, lo[2] + h];
        table.set_cube(lo, h);
        let seed = cube_verts(lo, h);
        cell_cuts(
            base, log, &order, &aabbs, lo, hi, h, &seed, &mut table, planes, &mut cuts,
        );
        if cuts.is_empty() {
            continue;
        }
        // The cube's own six faces bound the seed, so they are part of the
        // arrangement even though `cell_cuts` never lists them.
        all.clear();
        all.extend(0..RESERVED as u16 - 1);
        all.extend(cuts.iter().copied());

        for a in 0..all.len() {
            let pa = table.planes[all[a] as usize];
            for b in a + 1..all.len() {
                if coincident(&pa, &table.planes[all[b] as usize], eps) {
                    out.coincident_pairs += 1;
                }
            }
        }

        pts.clear();
        for a in 0..all.len() {
            for b in a + 1..all.len() {
                for c in b + 1..all.len() {
                    let Some(p) = solve3(
                        &table.planes[all[a] as usize],
                        &table.planes[all[b] as usize],
                        &table.planes[all[c] as usize],
                    ) else {
                        continue;
                    };
                    if (0..3).any(|x| p[x] < lo[x] - eps || p[x] > hi[x] + eps) {
                        continue;
                    }
                    let mut through = 0usize;
                    for &id in &all {
                        let pl = &table.planes[id as usize];
                        let l = pl.norm();
                        if l > 0.0 && (pl.at(p) / l).abs() <= eps {
                            through += 1;
                        }
                    }
                    if through <= 3 {
                        continue;
                    }
                    out.max_planes_at_a_point = out.max_planes_at_a_point.max(through);
                    if !pts
                        .iter()
                        .any(|q| (0..3).all(|x| (q[x] - p[x]).abs() <= eps))
                    {
                        pts.push(p);
                    }
                }
            }
        }
        if !pts.is_empty() {
            out.nonsimple_cells += 1;
            out.nonsimple_points += pts.len();
        }
    }
    out
}

/// The log with brush 0's centre displaced along `x` by `dx`.
///
/// The `M-44` control's only moving part. Shape-preserving: a sphere stays a
/// sphere, so the plane count and the `k`-DOP direction set are untouched and
/// the only thing that can move the hash is the displacement.
fn nudge(log: &[Edit], dx: f64) -> Vec<Edit> {
    let mut out: Vec<Edit> = log.to_vec();
    out[0].brush.shape = match out[0].brush.shape {
        Shape::Sphere(s) => Shape::Sphere(Sphere {
            center: [s.center[0] + dx, s.center[1], s.center[2]],
            radius: s.radius,
        }),
        Shape::Cube(b) => Shape::Cube(BoxExact {
            center: [b.center[0] + dx, b.center[1], b.center[2]],
            half_extents: b.half_extents,
        }),
        Shape::Capsule(c) => Shape::Capsule(Capsule {
            a: [c.a[0] + dx, c.a[1], c.a[2]],
            b: c.b,
            radius: c.radius,
        }),
    };
    out
}

/// What the 40,320-ordering sweep found for one `(fixture, mode)` arm.
struct Sweep {
    orderings: usize,
    distinct: usize,
    distinct_raw: usize,
    /// The leaf plane-triple sets alone: is the *arrangement* order-free?
    distinct_topo: usize,
    /// Those plus the vertex bit patterns: is the *arithmetic* order-free?
    distinct_pos: usize,
    /// Those plus the solid classification instead of the coordinates: is the
    /// *union* order-free?
    distinct_solid: usize,
    hash: u64,
    hash_raw: u64,
    cells_hashed: usize,
    solid_pieces: usize,
    total_pieces: usize,
    flags: Flags,
    /// Leaves and solid leaves of the two extreme `hash_topo` classes. Equal
    /// counts with different hashes means the same leaves cut differently;
    /// unequal counts means the arrangement lost or gained one, which is a
    /// different failure and a worse one.
    topo_lo_pieces: usize,
    topo_hi_pieces: usize,
    topo_lo_solid: usize,
    topo_hi_solid: usize,
    perturbed_differs: bool,
    perturbed_topo_differs: bool,
    seconds: f64,
}

/// Run every ordering of `log` and count distinct partitions.
fn sweep(
    base: &Base,
    log: &[Edit],
    chunk: Chunk,
    ds: &[[f64; 3]],
    mode: Mode,
    planes: Planes,
    perms: &[[usize; BRUSHES]],
    threads: usize,
) -> Sweep {
    let t0 = Instant::now();
    let total = chunk.cells * chunk.cells * chunk.cells;
    let chunks: Vec<Vec<[u64; 5]>> = std::thread::scope(|s| {
        let per = perms.len().div_ceil(threads);
        let mut hs = Vec::with_capacity(threads);
        for t in 0..threads {
            let slice = &perms[(t * per).min(perms.len())..((t + 1) * per).min(perms.len())];
            hs.push(s.spawn(move || {
                let mut table = Table::build(log, ds);
                let mut sc = Scratch::new();
                let mut out = Vec::with_capacity(slice.len());
                for perm in slice {
                    let p = partition(
                        base,
                        log,
                        perm,
                        chunk,
                        &mut table,
                        mode,
                        planes,
                        0..total,
                        true,
                        &mut sc,
                    );
                    out.push([p.hash, p.hash_raw, p.hash_topo, p.hash_pos, p.hash_solid]);
                }
                out
            }));
        }
        hs.into_iter().map(|h| h.join().expect("worker")).collect()
    });
    // One column per projection, each counted independently: a `distinct` above
    // one is only a finding once it says *which* projection split.
    let mut cols: [Vec<u64>; 5] = std::array::from_fn(|_| Vec::with_capacity(perms.len()));
    for c in &chunks {
        for row in c {
            for (col, &h) in cols.iter_mut().zip(row.iter()) {
                col.push(h);
            }
        }
    }
    let orderings = cols[0].len();
    assert_eq!(
        orderings,
        perms.len(),
        "every ordering must contribute a hash"
    );
    let mut distinct = [0usize; 5];
    let mut lowest = [0u64; 5];
    for (i, col) in cols.iter_mut().enumerate() {
        col.sort_unstable();
        lowest[i] = *col.first().expect("at least one ordering");
        col.dedup();
        distinct[i] = col.len();
    }

    // Name the two extreme topology classes by rebuilding a representative
    // ordering of each. When `distinct_topo == 1` both representatives are the
    // same partition and the two pairs agree, which is the reading that says the
    // comparison was live and found nothing.
    let flat: Vec<[u64; 5]> = chunks.iter().flatten().copied().collect();
    let rep = |target: u64| {
        let idx = flat
            .iter()
            .position(|r| r[2] == target)
            .expect("a hash in the column came from some ordering");
        let ord: Vec<usize> = perms[idx].to_vec();
        partition_all(base, log, &ord, chunk, ds, mode, planes, true, threads)
    };
    let topo_lo = rep(lowest[2]);
    let topo_hi = rep(*cols[2].last().expect("at least one ordering"));

    // The identity ordering again, for the shape counters and the M-44 control.
    let identity: Vec<usize> = (0..log.len()).collect();
    let id_part = partition_all(base, log, &identity, chunk, ds, mode, planes, true, threads);

    // `M-44`, twice over, because `distinct` and `distinct_topo` are zeros over
    // two different instruments and one control cannot license both.
    //
    // The **fine** nudge displaces brush 0's centre by `1e-9`. That is far below
    // a cell, so it moves no plane across a cell boundary and the arrangement is
    // combinatorially identical — it can only show up in the coordinates. It is
    // therefore the control for `distinct` and `distinct_pos` and it is *not* a
    // control for `distinct_topo`.
    //
    // The **coarse** nudge displaces it by a quarter of a cell, which does change
    // which cells the brush's planes cross, so the leaf plane-triple sets change.
    // Without it, `distinct_topo == 1` would be a zero over an instrument nobody
    // ever showed could return anything else — `P-70`'s C3, a held clause with no
    // instrument, which is the worst outcome available here.
    let fine = nudge(log, 1e-9);
    let fine_part = partition_all(
        base, &fine, &identity, chunk, ds, mode, planes, true, threads,
    );
    let coarse = nudge(log, 0.25 * chunk.cell_size());
    let coarse_part = partition_all(
        base, &coarse, &identity, chunk, ds, mode, planes, true, threads,
    );

    Sweep {
        orderings,
        distinct: distinct[0],
        distinct_raw: distinct[1],
        distinct_topo: distinct[2],
        distinct_pos: distinct[3],
        distinct_solid: distinct[4],
        hash: lowest[0],
        hash_raw: lowest[1],
        cells_hashed: id_part.cells_hashed,
        solid_pieces: id_part.solid_pieces,
        total_pieces: id_part.total_pieces,
        flags: id_part.flags,
        topo_lo_pieces: topo_lo.total_pieces,
        topo_hi_pieces: topo_hi.total_pieces,
        topo_lo_solid: topo_lo.solid_pieces,
        topo_hi_solid: topo_hi.solid_pieces,
        perturbed_differs: fine_part.hash != id_part.hash,
        perturbed_topo_differs: coarse_part.hash_topo != id_part.hash_topo,
        seconds: t0.elapsed().as_secs_f64(),
    }
}

// ---------------------------------------------------------------------------
// Sampled volume, for the fidelity column
// ---------------------------------------------------------------------------

/// The folded field's solid volume inside the chunk, by lattice sampling.
fn sampled_volume(base: &Base, log: &[Edit], chunk: Chunk, n: u32) -> f64 {
    let order: Vec<usize> = (0..log.len()).collect();
    let h = chunk.extent / f64::from(n);
    let inside: u64 = std::thread::scope(|s| {
        let per = n.div_ceil(MAX_THREADS as u32);
        let mut hs = Vec::new();
        for t in 0..MAX_THREADS as u32 {
            let klo = t * per;
            let khi = ((t + 1) * per).min(n);
            let order = order.clone();
            hs.push(s.spawn(move || {
                let mut c = 0u64;
                for k in klo..khi {
                    for j in 0..n {
                        for i in 0..n {
                            let p = [
                                chunk.origin[0] + (f64::from(i) + 0.5) * h,
                                chunk.origin[1] + (f64::from(j) + 0.5) * h,
                                chunk.origin[2] + (f64::from(k) + 0.5) * h,
                            ];
                            if fold(base, log, &order, p) <= 0.0 {
                                c += 1;
                            }
                        }
                    }
                }
                c
            }));
        }
        hs.into_iter().map(|h| h.join().expect("worker")).sum()
    });
    inside as f64 * h * h * h
}

// ---------------------------------------------------------------------------
// Rows
// ---------------------------------------------------------------------------

type Row = Vec<(&'static str, String)>;

/// One `(field, chunk_cells, dop)` arm of C1 and C3.
struct C13 {
    base: Base,
    chunk: Chunk,
    dop: usize,
}

fn run_c13(a: &C13, planes: Planes, threads: usize, probe: &mut Probe) -> Row {
    let ds = dirs(a.dop);
    let h = a.chunk.cell_size();
    let identity_empty: Vec<usize> = Vec::new();
    let empty_log: Vec<Edit> = Vec::new();

    // C3's "before": the unedited chunk's convex-cell compound.
    let before = partition_all(
        &a.base,
        &empty_log,
        &identity_empty,
        a.chunk,
        &ds,
        Mode::Lerp,
        planes,
        false,
        threads,
    );

    // C3's "after": `M-50`'s largest bucket.
    let log = dig_log(&a.base, a.chunk, SIXTY, BRUSH_RADIUS_CELLS * h);
    let order: Vec<usize> = (0..log.len()).collect();
    let after = partition_all(
        &a.base,
        &log,
        &order,
        a.chunk,
        &ds,
        Mode::Lerp,
        planes,
        false,
        threads,
    );

    // C1: one more brush arrives on the finished compound. It lands on the
    // surface at the first column `dig_log`'s golden-angle walk visits, so the
    // impact site is in the fixture rather than chosen after seeing a number.
    let impact_log = dig_log(&a.base, a.chunk, 1, IMPACT_RADIUS_CELLS * h);
    let impact = Shape::Sphere(Sphere {
        center: match impact_log[0].brush.shape {
            Shape::Sphere(s) => s.center,
            Shape::Cube(b) => b.center,
            Shape::Capsule(c) => c.a,
        },
        radius: IMPACT_RADIUS_CELLS * h,
    });
    let impact_planes = dop_planes(&impact, &ds);

    let (comp, mut table) = compound_near(&a.base, &log, a.chunk, &ds, &impact, Mode::Lerp, planes);
    // The impact's own plane block: one past the last edit's.
    let impact_base_id = (RESERVED + log.len() * ds.len()) as u16;
    table.planes.resize(
        usize::from(impact_base_id) + ds.len(),
        Plane {
            n: [0.0; 3],
            d: 0.0,
        },
    );

    let mut sc = Scratch::new();
    let mut totals: Vec<f64> = Vec::with_capacity(REPS);
    let mut worst = 0.0f64;
    let mut frags = 0usize;
    let mut rem = 0usize;
    let mut cells = 0usize;
    for _ in 0..REPS {
        let f = fracture(
            &comp,
            &mut table,
            &impact_planes,
            impact_base_id,
            a.chunk,
            Mode::Lerp,
            &mut sc,
        );
        totals.push(f.total_ms);
        worst = worst.max(f.worst_per_fragment_ms);
        frags = f.fragments;
        rem = f.remainder;
        cells = f.cells;
    }
    totals.sort_unstable_by(f64::total_cmp);
    let median_ms = totals[REPS / 2];

    // The counted rep, for `M-280`'s cycles.
    probe.reset_and_enable();
    let t = Instant::now();
    let counted = fracture(
        &comp,
        &mut table,
        &impact_planes,
        impact_base_id,
        a.chunk,
        Mode::Lerp,
        &mut sc,
    );
    let counted_ns = t.elapsed().as_secs_f64() * 1e9;
    probe.disable();
    let counts = probe.read();
    assert!(
        counts.worst_ratio() >= MIN_TIME_RATIO,
        "a counter was multiplexed at ratio {:.4}, so `ghz` would be an extrapolation",
        counts.worst_ratio()
    );

    assert!(
        frags > 0,
        "the fracture produced no fragments, so C1 has nothing to fail on"
    );
    assert_eq!(
        counted.fragments, frags,
        "the counted rep produced a different fragment count from the timed reps"
    );
    assert!(
        before.solid_pieces > 0,
        "the unedited chunk has no convex cells, so C3's ratio is vacuous"
    );
    assert!(
        after.cells_cut_by_brushes > 0,
        "no cell was cut by a brush plane, so the 60-brush arm measured nothing"
    );
    assert_eq!(
        after.flags.degenerate_verts, 0,
        "a vertex carried a repeated plane"
    );
    assert_eq!(
        after.flags.singular_triples, 0,
        "a plane triple was near-singular"
    );
    assert_eq!(
        after.flags.convexity_violations, 0,
        "a piece was not convex, which is the invariant this whole construction rests on"
    );
    assert!(
        after.flags.worst_volume_error < 1e-9,
        "a cell's leaves did not tile its cube: relative volume error {:.3e}",
        after.flags.worst_volume_error
    );

    let per_frag = median_ms / frags as f64;
    let ratio = after.solid_pieces as f64 / before.solid_pieces as f64;
    let c1 = per_frag < C1_BAR_MS;
    let c3 = ratio < C3_BAR;
    let field_vol = sampled_volume(&a.base, &log, a.chunk, VOLUME_SAMPLES);

    vec![
        ("arm", "c1c3".to_string()),
        ("field", a.base.name().to_string()),
        ("chunk_cells", a.chunk.cells.to_string()),
        ("cell_size_world", format!("{h:.9}")),
        ("dop_dirs", a.dop.to_string()),
        ("vertex_mode", Mode::Lerp.name().to_string()),
        ("planes_mode", planes.name().to_string()),
        ("brushes", SIXTY.to_string()),
        ("brush_radius_cells", format!("{BRUSH_RADIUS_CELLS:.2}")),
        ("impact_radius_cells", format!("{IMPACT_RADIUS_CELLS:.2}")),
        ("fragments", frags.to_string()),
        ("fracture_remainder_pieces", rem.to_string()),
        ("fracture_cells", cells.to_string()),
        ("fracture_ms_total", format!("{median_ms:.6}")),
        ("fracture_ms_per_fragment", format!("{per_frag:.6}")),
        ("fracture_ms_worst", format!("{worst:.6}")),
        ("reps", REPS.to_string()),
        ("cycles", counts.cycles.count.to_string()),
        (
            "ghz",
            format!("{:.4}", counts.cycles.count as f64 / counted_ns),
        ),
        (
            "cycles_per_fragment",
            format!("{:.1}", counts.cycles.count as f64 / frags as f64),
        ),
        ("c1_bar_ms", format!("{C1_BAR_MS:.1}")),
        ("c1_falsifier_ms", format!("{C1_FALSIFIER_MS:.1}")),
        (
            "m116_speedup_mean_low",
            format!("{:.1}", M116_MEAN_LOW_MS / per_frag),
        ),
        (
            "m116_speedup_mean_high",
            format!("{:.1}", M116_MEAN_HIGH_MS / per_frag),
        ),
        (
            "m116_speedup_worst",
            format!("{:.1}", M116_WORST_MS / worst.max(f64::MIN_POSITIVE)),
        ),
        ("convex_cells_before", before.solid_pieces.to_string()),
        ("convex_cells_after", after.solid_pieces.to_string()),
        ("cell_growth_ratio", format!("{ratio:.6}")),
        ("cells_full_before", before.cells_full.to_string()),
        ("cells_boundary_before", before.cells_boundary.to_string()),
        (
            "cells_cut_by_brushes",
            after.cells_cut_by_brushes.to_string(),
        ),
        ("total_pieces_before", before.total_pieces.to_string()),
        ("total_pieces_after", after.total_pieces.to_string()),
        (
            "max_pieces_in_a_cell",
            after.max_pieces_in_a_cell.to_string(),
        ),
        ("max_cuts_in_a_cell", after.flags.max_cuts.to_string()),
        ("max_leaves_in_a_cell", after.flags.max_leaves.to_string()),
        ("max_verts_seen", after.flags.max_verts.to_string()),
        (
            "degenerate_vertices",
            after.flags.degenerate_verts.to_string(),
        ),
        ("singular_triples", after.flags.singular_triples.to_string()),
        ("slivers_dropped", after.flags.slivers.to_string()),
        (
            "convexity_violations",
            after.flags.convexity_violations.to_string(),
        ),
        (
            "partition_volume_error_rel",
            format!("{:.3e}", after.flags.worst_volume_error),
        ),
        (
            "compound_volume_after",
            format!("{:.6}", after.compound_volume),
        ),
        ("field_volume_sampled", format!("{field_vol:.6}")),
        (
            "volume_error_rel",
            format!(
                "{:.6}",
                (after.compound_volume - field_vol).abs() / field_vol
            ),
        ),
        ("volume_samples", VOLUME_SAMPLES.to_string()),
        ("threads", threads.to_string()),
        ("orderings", "NA".to_string()),
        ("distinct_partitions", "NA".to_string()),
        ("distinct_partitions_raw", "NA".to_string()),
        ("distinct_topology", "NA".to_string()),
        ("distinct_positions", "NA".to_string()),
        ("distinct_solid_flags", "NA".to_string()),
        ("topo_class_lo_pieces", "NA".to_string()),
        ("topo_class_hi_pieces", "NA".to_string()),
        ("topo_class_lo_solid", "NA".to_string()),
        ("topo_class_hi_solid", "NA".to_string()),
        ("coincident_plane_pairs", "NA".to_string()),
        ("nonsimple_points", "NA".to_string()),
        ("nonsimple_cells", "NA".to_string()),
        ("max_planes_at_a_point", "NA".to_string()),
        ("partition_hash", "NA".to_string()),
        ("partition_hash_raw", "NA".to_string()),
        ("perturbed_hash_differs", "NA".to_string()),
        ("perturbed_topology_differs", "NA".to_string()),
        ("permutations_distinct", "NA".to_string()),
        ("cells_hashed", "NA".to_string()),
        ("sweep_seconds", "NA".to_string()),
        ("c1_holds", c1.to_string()),
        ("c2_holds", "NA".to_string()),
        ("c3_holds", c3.to_string()),
    ]
}

fn run_c2(
    base: &Base,
    log: &[Edit],
    chunk: Chunk,
    dop: usize,
    mode: Mode,
    planes: Planes,
    perms: &[[usize; BRUSHES]],
    threads: usize,
) -> Row {
    let ds = dirs(dop);
    let s = sweep(base, log, chunk, &ds, mode, planes, perms, threads);
    let simp = simplicity_scan(base, log, chunk, &ds, planes);
    assert_eq!(
        s.orderings, ORDERINGS,
        "the registered vacuity control: the arm must reach all 8! orderings"
    );
    assert!(
        s.cells_hashed > 0,
        "no cell was cut by a brush plane, so the hash covers nothing"
    );
    assert!(
        s.perturbed_differs,
        "displacing a brush centre by 1e-9 did not change the partition hash, so \
         distinct_partitions is a zero over an instrument that cannot say `different`"
    );
    assert!(
        s.perturbed_topo_differs,
        "displacing a brush centre by a quarter of a cell did not change the leaf plane-triple \
         sets, so distinct_topology is a zero over an instrument that cannot say `different`"
    );
    assert_eq!(
        s.flags.degenerate_verts, 0,
        "a vertex carried a repeated plane"
    );
    assert_eq!(
        s.flags.singular_triples, 0,
        "a plane triple was near-singular"
    );
    assert_eq!(
        s.flags.convexity_violations, 0,
        "a piece was not convex, which is the invariant this construction rests on"
    );
    assert!(
        s.flags.worst_volume_error < 1e-9,
        "a cell's leaves did not tile its cube: relative volume error {:.3e}",
        s.flags.worst_volume_error
    );
    let c2 = s.distinct == 1;
    vec![
        ("arm", "c2".to_string()),
        ("field", base.name().to_string()),
        ("chunk_cells", chunk.cells.to_string()),
        ("cell_size_world", format!("{:.9}", chunk.cell_size())),
        ("dop_dirs", dop.to_string()),
        ("vertex_mode", mode.name().to_string()),
        ("planes_mode", planes.name().to_string()),
        ("brushes", log.len().to_string()),
        ("brush_radius_cells", "NA".to_string()),
        ("impact_radius_cells", "NA".to_string()),
        ("fragments", "NA".to_string()),
        ("fracture_remainder_pieces", "NA".to_string()),
        ("fracture_cells", "NA".to_string()),
        ("fracture_ms_total", "NA".to_string()),
        ("fracture_ms_per_fragment", "NA".to_string()),
        ("fracture_ms_worst", "NA".to_string()),
        ("reps", "NA".to_string()),
        ("cycles", "NA".to_string()),
        ("ghz", "NA".to_string()),
        ("cycles_per_fragment", "NA".to_string()),
        ("c1_bar_ms", format!("{C1_BAR_MS:.1}")),
        ("c1_falsifier_ms", format!("{C1_FALSIFIER_MS:.1}")),
        ("m116_speedup_mean_low", "NA".to_string()),
        ("m116_speedup_mean_high", "NA".to_string()),
        ("m116_speedup_worst", "NA".to_string()),
        ("convex_cells_before", "NA".to_string()),
        ("convex_cells_after", s.solid_pieces.to_string()),
        ("cell_growth_ratio", "NA".to_string()),
        ("cells_full_before", "NA".to_string()),
        ("cells_boundary_before", "NA".to_string()),
        ("cells_cut_by_brushes", s.cells_hashed.to_string()),
        ("total_pieces_before", "NA".to_string()),
        ("total_pieces_after", s.total_pieces.to_string()),
        ("max_pieces_in_a_cell", "NA".to_string()),
        ("max_cuts_in_a_cell", s.flags.max_cuts.to_string()),
        ("max_leaves_in_a_cell", s.flags.max_leaves.to_string()),
        ("max_verts_seen", s.flags.max_verts.to_string()),
        ("degenerate_vertices", s.flags.degenerate_verts.to_string()),
        ("singular_triples", s.flags.singular_triples.to_string()),
        ("slivers_dropped", s.flags.slivers.to_string()),
        (
            "convexity_violations",
            s.flags.convexity_violations.to_string(),
        ),
        (
            "partition_volume_error_rel",
            format!("{:.3e}", s.flags.worst_volume_error),
        ),
        ("compound_volume_after", "NA".to_string()),
        ("field_volume_sampled", "NA".to_string()),
        ("volume_error_rel", "NA".to_string()),
        ("volume_samples", "NA".to_string()),
        ("threads", threads.to_string()),
        ("orderings", s.orderings.to_string()),
        ("distinct_partitions", s.distinct.to_string()),
        ("distinct_partitions_raw", s.distinct_raw.to_string()),
        ("distinct_topology", s.distinct_topo.to_string()),
        ("distinct_positions", s.distinct_pos.to_string()),
        ("distinct_solid_flags", s.distinct_solid.to_string()),
        ("topo_class_lo_pieces", s.topo_lo_pieces.to_string()),
        ("topo_class_hi_pieces", s.topo_hi_pieces.to_string()),
        ("topo_class_lo_solid", s.topo_lo_solid.to_string()),
        ("topo_class_hi_solid", s.topo_hi_solid.to_string()),
        ("coincident_plane_pairs", simp.coincident_pairs.to_string()),
        ("nonsimple_points", simp.nonsimple_points.to_string()),
        ("nonsimple_cells", simp.nonsimple_cells.to_string()),
        (
            "max_planes_at_a_point",
            simp.max_planes_at_a_point.to_string(),
        ),
        ("partition_hash", format!("{:016x}", s.hash)),
        ("partition_hash_raw", format!("{:016x}", s.hash_raw)),
        ("perturbed_hash_differs", s.perturbed_differs.to_string()),
        (
            "perturbed_topology_differs",
            s.perturbed_topo_differs.to_string(),
        ),
        ("permutations_distinct", perms.len().to_string()),
        ("cells_hashed", s.cells_hashed.to_string()),
        ("sweep_seconds", format!("{:.1}", s.seconds)),
        ("c1_holds", "NA".to_string()),
        ("c2_holds", c2.to_string()),
        ("c3_holds", "NA".to_string()),
    ]
}

fn main() {
    if !std::env::args().any(|a| a == "--bench") {
        return;
    }
    common::experiment::run(isomesh::experiment!("P-84"), |run| {
        let threads = std::thread::available_parallelism()
            .map_or(1, std::num::NonZeroUsize::get)
            .min(MAX_THREADS);

        // The registered vacuity control, before anything else runs.
        let perms = permutations();
        assert_eq!(
            perms.len(),
            ORDERINGS,
            "the 40,320-ordering arm must actually reach every ordering"
        );
        let mut seen = perms.clone();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(
            seen.len(),
            ORDERINGS,
            "the permutations must be pairwise distinct, or the count is not a count of orderings"
        );

        let mut rows: Vec<Row> = Vec::new();
        let mut probe = Probe::open();

        for arm in [
            C13 {
                base: Base::Fbm(FbmTerrain::canonical()),
                chunk: Chunk {
                    origin: [-2.0, -2.0, -2.0],
                    extent: 4.0,
                    cells: 16,
                },
                dop: 6,
            },
            C13 {
                base: Base::Fbm(FbmTerrain::canonical()),
                chunk: Chunk {
                    origin: [-2.0, -2.0, -2.0],
                    extent: 4.0,
                    cells: 32,
                },
                dop: 6,
            },
            C13 {
                base: Base::Gyroid(Gyroid::canonical()),
                chunk: Chunk {
                    origin: [-2.0, -2.0, -2.0],
                    extent: 4.0,
                    cells: 16,
                },
                dop: 6,
            },
            C13 {
                base: Base::Gyroid(Gyroid::canonical()),
                chunk: Chunk {
                    origin: [-2.0, -2.0, -2.0],
                    extent: 4.0,
                    cells: 32,
                },
                dop: 6,
            },
        ] {
            for planes in [Planes::AsLogged, Planes::Canonical] {
                let row = run_c13(&arm, planes, threads, &mut probe);
                println!(
                    "P-84 c1c3 {} at {}³, {}-DOP, {} planes: {} fragments, {} ms/fragment, \
                     worst {} ms; cells {} -> {} = {}x",
                    arm.base.name(),
                    arm.chunk.cells,
                    arm.dop,
                    planes.name(),
                    find(&row, "fragments"),
                    find(&row, "fracture_ms_per_fragment"),
                    find(&row, "fracture_ms_worst"),
                    find(&row, "convex_cells_before"),
                    find(&row, "convex_cells_after"),
                    find(&row, "cell_growth_ratio"),
                );
                rows.push(row);
            }
        }

        // C2: `M-36`'s own fixture, on `M-38`'s base cube, across both forks —
        // how a crossing vertex is computed, and whether coincident planes are
        // collapsed. Four arms, because the measurement found the second fork
        // matters and the registration only anticipated the first.
        let m36 = eight();
        let c2_chunk = Chunk {
            origin: [-1.0, -1.0, -1.0],
            extent: 2.0,
            cells: 8,
        };
        let base = Base::Box(BoxExact::canonical());
        for mode in [Mode::Lerp, Mode::Solve] {
            for planes in [Planes::AsLogged, Planes::Canonical] {
                let row = run_c2(&base, &m36, c2_chunk, 6, mode, planes, &perms, threads);
                println!(
                    "P-84 c2 {} vertices, {} planes: {} orderings, {} distinct partitions \
                     ({} raw); topology {}, positions {}, solid flags {}; coincident pairs {}, \
                     non-simple points {}; hash {}, perturbation seen fine/coarse: {}/{}, {} s",
                    mode.name(),
                    planes.name(),
                    find(&row, "orderings"),
                    find(&row, "distinct_partitions"),
                    find(&row, "distinct_partitions_raw"),
                    find(&row, "distinct_topology"),
                    find(&row, "distinct_positions"),
                    find(&row, "distinct_solid_flags"),
                    find(&row, "coincident_plane_pairs"),
                    find(&row, "nonsimple_points"),
                    find(&row, "partition_hash"),
                    find(&row, "perturbed_hash_differs"),
                    find(&row, "perturbed_topology_differs"),
                    find(&row, "sweep_seconds"),
                );
                rows.push(row);
            }
        }

        for row in rows {
            run.record(&row);
        }
    });
}

/// One column of a built row, for the progress lines.
fn find(row: &Row, key: &str) -> String {
    row.iter()
        .find(|(k, _)| *k == key)
        .map_or_else(|| "?".to_string(), |(_, v)| v.clone())
}
