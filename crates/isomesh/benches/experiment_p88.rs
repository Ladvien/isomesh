//! **P-88 — a clearance lower bound read off an octree of free cells, with no
//! medial axis, and whether it can ride `P-87`'s local repair.**
//!
//! Ticket: R-088. Pre-registered before this harness existed. Rides `P-87`.
//!
//! ```bash
//! cargo bench --bench experiment_p88
//! ```
//!
//! Writes `docs/experiments/p-88.csv`.
//!
//! # What is being asked
//!
//! CALIBRE gives every creature a half-width `λ` and defines reachable space as
//! the connected component of `{r ≥ λ}`, where `r` is clearance. The docs say
//! the hard sub-problem is maintaining the medial axis `ρ` under material
//! removal, and the defining paper (Chazal & Lieutier,
//! `10.1016/j.gmod.2005.01.002`) is unobtainable. The registration's claim is
//! that an octree of free cells sidesteps `ρ` entirely: *"a cell of side `s`
//! that is entirely empty admits a sphere of radius `s/2`, and `P-87`'s merged
//! convex regions give a larger bound directly"*.
//!
//! # SHARE, recomputed before a line of this was written
//!
//! The registration's SHARE line is *"this is an accuracy and a rate, not a
//! ratio of a total"*, and that is **correct as written**. `✗51`'s ceiling
//! `1/(1 − s + s/factor)` binds a *speedup* denominated in a share `s` of some
//! larger total; none of these three clauses is a speedup, so no share can
//! throttle them:
//!
//! - **C1 is an accuracy rate.** `within_one_voxel / samples`, denominator the
//!   sample set itself — share 1.0 by construction. Reachable **upward**: a
//!   bound that is exact scores 1.0. Reachable **downward**: a bound that
//!   reports the inradius of a one-cell-thick slab inside an eight-cell channel
//!   is 7.5 cells short, so a fixture with a fat channel drives the fraction
//!   toward 0. Both directions are live before the run.
//! - **C2 is a ratio of two directly counted sets** — λ-membership flips over
//!   sign-flipped samples — plus a **subset test**. Reachable upward: dissolving
//!   one merged region rewrites the bound of every cell in it, and a region can
//!   hold thousands of cells against a few hundred flipped samples. Reachable
//!   downward: an edit whose repair set is a handful of leaves flips a handful
//!   of cells. The subset half is reachable in both directions and this harness
//!   *demonstrates* the non-zero direction rather than asserting it is possible
//!   — see the controls below.
//! - **C3 is a one-sided count** over ≥ 10⁶ samples, and a count has no
//!   denominator to be diluted by. Reachable downward (a sound bound scores 0);
//!   reachable upward is **demonstrated on the same samples** by the half-voxel
//!   probe.
//!
//! What is *not* reachable, and is recorded here rather than discovered
//! afterwards: no clause here is a comparison against another machine's number,
//! so `M-281` does not bind. `clock_mhz` is on every row anyway (`M-280`)
//! because `bound_ns_per_query` is on every row.
//!
//! # The ambiguity in "the octree clearance lower bound", and how it is resolved
//!
//! The registration's own sentence contains **two** estimators and they are not
//! the same function:
//!
//! | id | bound at `p` | cost |
//! |---|---|---|
//! | `region` | `d(p, ∂B)` for `B` the merged convex region containing `p` | one leaf walk, `O(depth)` |
//! | `ball` | `d(p, ⋃ non-free cells)`, i.e. the largest ball centred at `p` inside free space | branch-and-bound descent of the same pyramid |
//!
//! `region` is the literal reading — *"merged convex regions give a larger bound
//! directly"* — and it is what a navigation layer gets for free once `P-87` has
//! run. `ball` is the reading CALIBRE's own definition invites, because
//! λ-membership *is* a thresholded ball query. `ball ≥ region` always, since the
//! region box is a subset of the free-cell union.
//!
//! **Both are run, over the same samples, the same trace and the same truth, and
//! both get a row.** Scoring only the literal one would answer the question with
//! a straw estimator; scoring only the tight one would be `P-70`'s C3 — quietly
//! answering an easier question. The `estimator` column says which is which and
//! the verdicts are per row.
//!
//! # The truth, and why there are two of them
//!
//! - **`analytic`** — `M-346`'s fixture, inherited literally. A `BoxExact` big
//!   enough to swallow the domain with a `Capsule` of radius `r` subtracted
//!   along `x`: subtraction is `max(field, −shape)`, so inside the channel the
//!   box term is ≤ −2 and the value *is* `r − ρ`, exactly, `Capsule` being an
//!   exact distance field. `M-346` measured **zero** error against that truth at
//!   `r ∈ {2, 4, 8}` cells, which is why it is the fixture C1 names.
//! - **`refined`** — on `fbm_terrain` and `gyroid` no closed-form distance
//!   exists, so the truth is the same query run on a **4× refined** sign field
//!   of the same edited field ([`REFINE`]). This is one-sided evidence in
//!   exactly `M-347`'s sense: refinement can *catch* a coarse cell that is
//!   all-air while the surface passes through it (`M-347` counted 4,063 such
//!   cells across eight fields, 97 on `gyroid` and 94 on `fbm_terrain` at 33³)
//!   and cannot prove there are none.
//!
//!   **The refinement factor is the strictness, and it is arithmetic rather
//!   than taste.** A coarse non-free cell at distance `D` is guaranteed to
//!   contain a non-free *fine* cell only in the corner holding its solid corner
//!   sample, so the refined truth can exceed the coarse bound by at most
//!   `√3·(1 − 1/f)`: **0.866** cells at `f = 2` and **1.299** at `f = 4`. At
//!   `f = 2` no reference row could ever miss the one-voxel bar — the first run
//!   of this harness scored `within_fraction = 1.000000` and
//!   `max_understatement_cells = 0.866025` on every one of them, which is that
//!   cap and not a measurement of the surface. `f = 4` puts the cap past one
//!   voxel so the bar can be missed, and the analytic rows remain the ones C1
//!   is registered on.
//!
//! Both truths are clamped by the distance to the world boundary, because both
//! estimators are, and for the same reason: nothing outside the sampled world is
//! known, so a navigation layer must not claim clearance through its own face.
//!
//! # Controls, each an assertion rather than a printed number
//!
//! - **VACUITY (registered), first half.** `known_clearance_samples` — points
//!   whose true clearance is known analytically. Asserted `> 0` on the slab rows
//!   and in the total.
//! - **VACUITY (registered), second half, and it is the one that matters.**
//!   `narrow_passage_samples` — points in **dug** passages narrower than two
//!   cells, i.e. sign-flipped from solid to air by an edit *and* of true
//!   clearance under one cell. A bound tested only in open space is a bound
//!   tested nowhere interesting. The fixture is built so it can fire: every
//!   trace ends with capsule cuts of radius **0.9 cells** whose axis runs down a
//!   line of **cell centres**, so the four surrounding sample lines (at
//!   `√2/2 = 0.707` cells) go to air and the next ring (at `1.58` cells) does
//!   not — a passage exactly one cell wide. Asserted `> 0` on every row.
//! - **C3's zero is proved capable of being non-zero (`M-44`).**
//!   `probe_overstatements` re-runs the identical comparison against
//!   `bound + 0.5` cells — the classic half-voxel error of measuring to an
//!   obstacle cell's *centre* instead of its box — and is asserted `> 0` on
//!   every row. A zero that could not have been non-zero is not a measurement.
//! - **`flips_outside_repair_set`'s zero, likewise.** Asserted: there were flips
//!   at all (`lambda_flips > 0`) and the repair set is a strict subset of the
//!   world (`repair_set_cells < world_cells`), so a flip *could* have landed
//!   outside; and, globally, that the identical predicate returns **non-zero**
//!   somewhere in this run — which it does, on the `ball` rows, and that is the
//!   result rather than the control.
//! - **The regions tile the free leaves, and each is an exactly tiled box.**
//!   `P-87`'s `Nav::audit`, kept verbatim in spirit: Σ leaf volumes equals the
//!   world, Σ region volumes equals Σ free-leaf volumes, and every region's box
//!   volume equals its own members' volumes. That control caught a leaf-id
//!   recycling bug on `P-87`'s first run and it is re-run here after the build
//!   and after every edit.
//! - **The dirty set is the crate's, not the harness's.** The harness's count of
//!   cells incident to a flipped sign bit is asserted equal to
//!   `mark_edit`'s `sign_changed_cells` over the same box, per edit — two
//!   independent walks, which is the defect P-72 found in itself.
//! - **The bound never exceeds the tight bound.** `region ≤ ball` asserted per
//!   sample; a merged region that is not inside the free-cell union breaks it.
//!
//! # What is deliberately not here
//!
//! `P-87`'s region **adjacency graph**. No P-88 clause reads an edge: the
//! clearance bound is a property of a region's box and the λ-filter is a
//! threshold on it. Building the graph would be 150 lines that no column
//! depends on.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::needless_range_loop,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::many_single_char_names
)]

mod common;

use std::collections::BTreeMap;
use std::time::Instant;

use isomesh::Sdf;
use isomesh::brush::{Brush, BrushStack, Capsule};
use isomesh::chunk::ChunkLayout;
use isomesh::chunk::dirty::{DirtySet, mark_edit};
use isomesh::fields::{BoxExact, FbmTerrain, Gyroid, ReferenceField, Sphere};

// ── the world ───────────────────────────────────────────────────────────────

/// World extent per axis. `BoxExact::canonical().domain()` is `[−2, 2]³`, which
/// is `M-346`'s grid and `P-87`'s world at once.
const EXTENT: f64 = 4.0;
/// World origin, the same centring `P-72` and `P-87` record.
const ORIGIN: f64 = -EXTENT * 0.5;
/// Chunk granularity in cells per axis, `P-87`'s.
const CHUNK_CELLS: u32 = 4;

/// `M-346`'s resolution: 65 samples over `[−2, 2]`, so 64 cells.
const SLAB_CELLS: u32 = 64;
/// `M-346`'s channel radii, in cells.
const SLAB_RADII: [u32; 3] = [2, 4, 8];
/// World sizes for the reference fields.
const REF_WORLDS: [u32; 2] = [64, 128];
/// Refinement factor for the reference fields' truth.
///
/// Four rather than two: the refined truth can exceed the coarse bound by at
/// most `√3·(1 − 1/f)`, which is 0.866 cells at `f = 2` — under the one-voxel
/// bar, so a two-times grid cannot fail C1 no matter what the field does. At
/// `f = 4` the cap is 1.299 and the bar is reachable in both directions.
const REFINE: u32 = 4;

/// Wide sphere brushes in a reference trace. `P-72`'s dig, at `P-87`'s radius.
const WIDE_EDITS: usize = 8;
/// Narrow capsule cuts closing every trace: the passage the vacuity control
/// needs.
const NARROW_EDITS: usize = 3;
/// Narrow cuts in a slab trace, where there is no wide dig.
const SLAB_NARROW_EDITS: usize = 11;
/// Wide brush radius in cells — `M-311`'s protocol and `P-87`'s.
const WIDE_CELLS: f64 = 6.0;
/// Narrow cut radius in cells.
///
/// Under 1 so that only the four sample lines at `√2/2 = 0.7071` cells from a
/// cell-centre axis flip, and over `√2/2` so that they all do. The result is a
/// passage exactly one cell across — narrower than the two cells the registered
/// vacuity control names.
const NARROW_CELLS: f64 = 0.9;

/// Samples per fixture. Seven fixtures clears 10⁶ for each estimator.
const SAMPLES_PER_FIXTURE: usize = 160_000;
/// Share of a reference fixture's samples forced into the narrow cut.
///
/// The tube is ~0.05% of a 64³ world by volume, so uniform sampling would put a
/// two-digit count in the column the registration says is the one that matters.
/// `P-62`'s precedent: add the stratum, report it, and report the uniform-only
/// accuracy beside it so the stratification cannot hide inside a fraction.
const NARROW_SHARE: f64 = 0.05;

/// Creature half-widths swept, in cells. The registered columns carry the worst.
const LAMBDAS: [f64; 5] = [0.5, 1.0, 1.5, 2.0, 3.0];

/// Column names for the per-λ sweep. Parallel to [`LAMBDAS`], and asserted so
/// at run time, because a silently mismatched pair would relabel every number.
const FLIP_FACTOR_COLS: [&str; LAMBDAS.len()] = [
    "flip_factor_l05",
    "flip_factor_l10",
    "flip_factor_l15",
    "flip_factor_l20",
    "flip_factor_l30",
];
/// Per-λ counterpart of [`FLIP_FACTOR_COLS`] for the subset half of C2.
const OUTSIDE_COLS: [&str; LAMBDAS.len()] = [
    "flips_outside_repair_set_l05",
    "flips_outside_repair_set_l10",
    "flips_outside_repair_set_l15",
    "flips_outside_repair_set_l20",
    "flips_outside_repair_set_l30",
];

/// No leaf, no region.
const NONE: u32 = u32::MAX;

/// Every corner sample of the cell is air.
const AIR: u8 = 0;
/// Every corner sample of the cell is solid.
const SOLID: u8 = 1;
/// The surface passes through: at the leaf there is nothing left to split, so it
/// is an obstacle.
const MIXED: u8 = 2;

// ── the sign field ──────────────────────────────────────────────────────────

/// One bit per sample, packed 64-to-a-word along `x`. `P-87`'s, and `dual.rs`'s
/// `inside` prepass layout (R-039): `is_inside(value)` rather than the IEEE sign
/// bit, because `-0.0 < 0.0` is false and exactly zero is outside.
struct Signs {
    words: Vec<u64>,
    bit_row: usize,
    dims: [u32; 3],
}

impl Signs {
    fn new(dims: [u32; 3]) -> Self {
        let bit_row = (dims[0] as usize).div_ceil(64);
        let words = vec![0u64; bit_row * dims[1] as usize * dims[2] as usize];
        Self {
            words,
            bit_row,
            dims,
        }
    }

    #[inline]
    fn slot(&self, s: [u32; 3]) -> (usize, u64) {
        let row = (s[2] as usize * self.dims[1] as usize + s[1] as usize) * self.bit_row;
        (row + (s[0] as usize >> 6), 1u64 << (s[0] & 63))
    }

    #[inline]
    fn solid(&self, s: [u32; 3]) -> bool {
        let (w, m) = self.slot(s);
        self.words[w] & m != 0
    }

    #[inline]
    fn set(&mut self, s: [u32; 3], solid: bool) {
        let (w, m) = self.slot(s);
        if solid {
            self.words[w] |= m;
        } else {
            self.words[w] &= !m;
        }
    }

    #[inline]
    fn cell_state(&self, c: [u32; 3]) -> u8 {
        let mut solid = 0u32;
        for dz in 0..2 {
            for dy in 0..2 {
                for dx in 0..2 {
                    if self.solid([c[0] + dx, c[1] + dy, c[2] + dz]) {
                        solid += 1;
                    }
                }
            }
        }
        match solid {
            0 => AIR,
            8 => SOLID,
            _ => MIXED,
        }
    }

    /// Sample the field over the whole grid.
    fn fill<F: Sdf<Scalar = f64>>(&mut self, field: &F, layout: &ChunkLayout<f64>) {
        let dims = self.dims;
        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    let p = layout.world_of_sample([i64::from(x), i64::from(y), i64::from(z)]);
                    self.set([x, y, z], field.sample(p) < 0.0);
                }
            }
        }
    }
}

// ── the octree, as a state pyramid over cells ───────────────────────────────

/// `levels[0]` is one state per cell; `levels[k]` is one state per `2^k` node.
/// `P-87`'s pyramid: a node's state is a pure function of its eight children, so
/// an edit's update is a bottom-up walk over its own ancestors and nothing else.
struct Tree {
    levels: Vec<Vec<u8>>,
    n: u32,
    depth: usize,
}

impl Tree {
    fn new(n: u32) -> Self {
        assert!(n.is_power_of_two(), "world must be a power of two in cells");
        let depth = n.trailing_zeros() as usize;
        let mut levels = Vec::with_capacity(depth + 1);
        for k in 0..=depth {
            let m = (n >> k) as usize;
            levels.push(vec![MIXED; m * m * m]);
        }
        Self { levels, n, depth }
    }

    #[inline]
    fn side(&self, level: usize) -> u32 {
        self.n >> level
    }

    #[inline]
    fn idx(&self, level: usize, c: [u32; 3]) -> usize {
        let m = (self.n >> level) as usize;
        (c[2] as usize * m + c[1] as usize) * m + c[0] as usize
    }

    #[inline]
    fn state(&self, level: usize, c: [u32; 3]) -> u8 {
        self.levels[level][self.idx(level, c)]
    }

    fn rebuild(&mut self, signs: &Signs) {
        let n = self.n as usize;
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    self.levels[0][(z * n + y) * n + x] =
                        signs.cell_state([x as u32, y as u32, z as u32]);
                }
            }
        }
        for l in 1..=self.depth {
            let m = (self.n >> l) as usize;
            let mc = m * 2;
            let (lower, upper) = self.levels.split_at_mut(l);
            let child = &lower[l - 1];
            let parent = &mut upper[0];
            for z in 0..m {
                for y in 0..m {
                    for x in 0..m {
                        let mut kids = [MIXED; 8];
                        for dz in 0..2 {
                            for dy in 0..2 {
                                for dx in 0..2 {
                                    kids[(dz * 2 + dy) * 2 + dx] = child
                                        [((2 * z + dz) * mc + (2 * y + dy)) * mc + (2 * x + dx)];
                                }
                            }
                        }
                        parent[(z * m + y) * m + x] = combine(kids);
                    }
                }
            }
        }
    }
}

/// A node is uniform only when all eight children agree and are not `MIXED`.
#[inline]
fn combine(kids: [u8; 8]) -> u8 {
    let f = kids[0];
    if f != MIXED && kids.iter().all(|&s| s == f) {
        f
    } else {
        MIXED
    }
}

// ── the tight bound: nearest non-free cell, by branch and bound ─────────────

/// Squared distance from `q` (cell coordinates) to an axis-aligned box.
#[inline]
fn box_dist2(q: [f64; 3], lo: [f64; 3], hi: [f64; 3]) -> f64 {
    let mut s = 0.0;
    for a in 0..3 {
        let d = if q[a] < lo[a] {
            lo[a] - q[a]
        } else if q[a] > hi[a] {
            q[a] - hi[a]
        } else {
            0.0
        };
        s += d * d;
    }
    s
}

/// `d(q, ⋃ non-free cells)`, in cells, clamped by `best` (an upper bound, in
/// cells — always the distance to the world boundary here).
///
/// Exact rather than approximate: a node whose state is `AIR` contains no
/// non-free cell and is pruned outright; any other node is descended in order of
/// its own distance to `q`, and a level-0 node reached this way *is* a non-free
/// unit cell whose box distance is the answer. This is the octree used the way
/// an octree is used for a nearest-obstacle query — no medial axis, no distance
/// transform, and no dependence on the world size beyond `depth`.
fn nearest_non_free(
    tree: &Tree,
    q: [f64; 3],
    best: f64,
    stack: &mut Vec<(f64, u8, [u32; 3])>,
) -> f64 {
    let mut best2 = best * best;
    if tree.state(tree.depth, [0; 3]) == AIR {
        return best;
    }
    stack.clear();
    stack.push((0.0, tree.depth as u8, [0; 3]));
    let mut kids: [(f64, [u32; 3]); 8] = [(0.0, [0; 3]); 8];
    while let Some((d2, l, c)) = stack.pop() {
        if d2 >= best2 {
            continue;
        }
        let l = l as usize;
        if l == 0 {
            best2 = d2;
            continue;
        }
        let child = l - 1;
        let side = f64::from(1u32 << child);
        let mut m = 0usize;
        for dz in 0..2 {
            for dy in 0..2 {
                for dx in 0..2 {
                    let cc = [2 * c[0] + dx, 2 * c[1] + dy, 2 * c[2] + dz];
                    if tree.state(child, cc) == AIR {
                        continue;
                    }
                    let lo = [
                        f64::from(cc[0]) * side,
                        f64::from(cc[1]) * side,
                        f64::from(cc[2]) * side,
                    ];
                    let hi = [lo[0] + side, lo[1] + side, lo[2] + side];
                    let dd = box_dist2(q, lo, hi);
                    if dd < best2 {
                        kids[m] = (dd, cc);
                        m += 1;
                    }
                }
            }
        }
        // Descending order, so the nearest child is popped first and prunes its
        // siblings before they are ever expanded.
        kids[..m].sort_unstable_by(|a, b| b.0.total_cmp(&a.0));
        for &(dd, cc) in &kids[..m] {
            stack.push((dd, child as u8, cc));
        }
    }
    best2.sqrt()
}

// ── boxes, leaves, regions ──────────────────────────────────────────────────

/// A half-open box in cell coordinates. `hi` is exclusive.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Box3 {
    lo: [u32; 3],
    hi: [u32; 3],
}

impl Box3 {
    #[inline]
    fn volume(&self) -> u64 {
        let mut v = 1u64;
        for a in 0..3 {
            v *= u64::from(self.hi[a] - self.lo[a]);
        }
        v
    }

    /// Distance from `q` to the boundary, i.e. the largest ball centred at `q`
    /// inside this box. Zero if `q` is outside.
    #[inline]
    fn inradius_at(&self, q: [f64; 3]) -> f64 {
        let mut r = f64::INFINITY;
        for a in 0..3 {
            r = r.min(q[a] - f64::from(self.lo[a]));
            r = r.min(f64::from(self.hi[a]) - q[a]);
        }
        r.max(0.0)
    }
}

/// One octree leaf: a uniform node, or a `MIXED` unit cell.
struct Leaf {
    level: u8,
    state: u8,
    coords: [u32; 3],
    region: u32,
    alive: bool,
}

/// One merged convex navigation cell.
struct Region {
    bx: Box3,
    leaves: Vec<u32>,
    alive: bool,
}

// ── the navigation structure ────────────────────────────────────────────────

/// Octree leaves and the merged convex regions over the free ones.
///
/// `P-87`'s, minus the adjacency graph, which no P-88 clause reads.
struct Nav {
    tree: Tree,
    leaf_at: Vec<Vec<u32>>,
    leaves: Vec<Leaf>,
    leaf_free: Vec<u32>,
    /// Leaf ids killed by the repair in progress, held back from `leaf_free`
    /// until it finishes.
    ///
    /// `P-87`'s region-tiling audit caught this on its first run: recycling an
    /// id inside one repair makes a dissolved region's member list point at a
    /// different node, and a `SOLID` leaf then reaches the free-space merge.
    leaf_retire: Vec<u32>,
    regions: Vec<Region>,
    region_free: Vec<u32>,
    live_leaves: u64,
    live_regions: u64,
}

impl Nav {
    fn new(n: u32) -> Self {
        let tree = Tree::new(n);
        let mut leaf_at = Vec::with_capacity(tree.depth + 1);
        for k in 0..=tree.depth {
            let m = (n >> k) as usize;
            leaf_at.push(vec![NONE; m * m * m]);
        }
        Self {
            tree,
            leaf_at,
            leaves: Vec::new(),
            leaf_free: Vec::new(),
            leaf_retire: Vec::new(),
            regions: Vec::new(),
            region_free: Vec::new(),
            live_leaves: 0,
            live_regions: 0,
        }
    }

    fn clear(&mut self) {
        for l in &mut self.leaf_at {
            l.fill(NONE);
        }
        self.leaves.clear();
        self.leaf_free.clear();
        self.leaf_retire.clear();
        self.regions.clear();
        self.region_free.clear();
        self.live_leaves = 0;
        self.live_regions = 0;
    }

    #[inline]
    fn leaf_box(&self, id: u32) -> Box3 {
        let leaf = &self.leaves[id as usize];
        let s = 1u32 << leaf.level;
        Box3 {
            lo: [leaf.coords[0] * s, leaf.coords[1] * s, leaf.coords[2] * s],
            hi: [
                (leaf.coords[0] + 1) * s,
                (leaf.coords[1] + 1) * s,
                (leaf.coords[2] + 1) * s,
            ],
        }
    }

    fn new_leaf(&mut self, level: usize, coords: [u32; 3], state: u8) -> u32 {
        let id = match self.leaf_free.pop() {
            Some(id) => {
                self.leaves[id as usize] = Leaf {
                    level: level as u8,
                    state,
                    coords,
                    region: NONE,
                    alive: true,
                };
                id
            }
            None => {
                self.leaves.push(Leaf {
                    level: level as u8,
                    state,
                    coords,
                    region: NONE,
                    alive: true,
                });
                (self.leaves.len() - 1) as u32
            }
        };
        let i = self.tree.idx(level, coords);
        self.leaf_at[level][i] = id;
        self.live_leaves += 1;
        id
    }

    fn kill_leaf(&mut self, id: u32) {
        let (level, coords) = {
            let leaf = &self.leaves[id as usize];
            (leaf.level as usize, leaf.coords)
        };
        let i = self.tree.idx(level, coords);
        self.leaf_at[level][i] = NONE;
        self.leaves[id as usize].alive = false;
        self.leaves[id as usize].region = NONE;
        self.leaf_retire.push(id);
        self.live_leaves -= 1;
    }

    /// Create leaves under `(level, coords)`, stopping at the first uniform node.
    fn derive_under(&mut self, level: usize, coords: [u32; 3], added: &mut Vec<u32>) {
        let mut stack = vec![(level, coords)];
        while let Some((l, c)) = stack.pop() {
            let st = self.tree.state(l, c);
            if st != MIXED || l == 0 {
                let id = self.new_leaf(l, c, st);
                added.push(id);
                continue;
            }
            for dz in 0..2 {
                for dy in 0..2 {
                    for dx in 0..2 {
                        stack.push((l - 1, [2 * c[0] + dx, 2 * c[1] + dy, 2 * c[2] + dz]));
                    }
                }
            }
        }
    }

    /// The existing leaves covering the subtree at `(level, coords)`.
    fn collect_under(&self, level: usize, coords: [u32; 3], out: &mut Vec<u32>) {
        let mut stack = vec![(level, coords)];
        while let Some((l, c)) = stack.pop() {
            let id = self.leaf_at[l][self.tree.idx(l, c)];
            if id != NONE {
                out.push(id);
                continue;
            }
            assert!(
                l > 0,
                "a level-0 node with no leaf and no uniform ancestor: the leaf structure is not \
                 a partition, which means the repair is wrong"
            );
            for dz in 0..2 {
                for dy in 0..2 {
                    for dx in 0..2 {
                        stack.push((l - 1, [2 * c[0] + dx, 2 * c[1] + dy, 2 * c[2] + dz]));
                    }
                }
            }
        }
    }

    fn new_region(&mut self, bx: Box3, leaves: Vec<u32>) -> u32 {
        let id = match self.region_free.pop() {
            Some(id) => {
                self.regions[id as usize] = Region {
                    bx,
                    leaves,
                    alive: true,
                };
                id
            }
            None => {
                self.regions.push(Region {
                    bx,
                    leaves,
                    alive: true,
                });
                (self.regions.len() - 1) as u32
            }
        };
        let members = core::mem::take(&mut self.regions[id as usize].leaves);
        for &l in &members {
            self.leaves[l as usize].region = id;
        }
        self.regions[id as usize].leaves = members;
        self.live_regions += 1;
        id
    }

    fn kill_region(&mut self, id: u32) {
        let members = core::mem::take(&mut self.regions[id as usize].leaves);
        for l in members {
            if self.leaves[l as usize].alive && self.leaves[l as usize].region == id {
                self.leaves[l as usize].region = NONE;
            }
        }
        self.regions[id as usize].alive = false;
        self.region_free.push(id);
        self.live_regions -= 1;
    }

    /// Merge a pool of free leaves into convex regions. Returns the new ids.
    fn merge_pool(&mut self, pool: &[u32], made: &mut Vec<u32>) {
        let mut items: Vec<(Box3, Vec<u32>)> =
            pool.iter().map(|&l| (self.leaf_box(l), vec![l])).collect();
        greedy_merge(&mut items);
        made.clear();
        for (bx, leaves) in items {
            let id = self.new_region(bx, leaves);
            made.push(id);
        }
    }

    /// Full build from the pyramid: leaves, then a global merge.
    fn build(&mut self) {
        self.clear();
        let mut added = Vec::new();
        let (depth, root) = (self.tree.depth, [0u32; 3]);
        self.derive_under(depth, root, &mut added);
        let pool: Vec<u32> = added
            .into_iter()
            .filter(|&l| self.leaves[l as usize].state == AIR)
            .collect();
        let mut made = Vec::new();
        self.merge_pool(&pool, &mut made);
    }

    /// The leaf covering cell `c`, or [`NONE`].
    #[inline]
    fn leaf_of_cell(&self, c: [u32; 3]) -> u32 {
        let mut l = 0usize;
        let mut cc = c;
        loop {
            let id = self.leaf_at[l][self.tree.idx(l, cc)];
            if id != NONE {
                return id;
            }
            if l == self.tree.depth {
                return NONE;
            }
            l += 1;
            cc = [cc[0] >> 1, cc[1] >> 1, cc[2] >> 1];
        }
    }

    /// `P-87`'s audit, kept: every invariant the structure has, checked at once.
    fn audit(&self, where_: &str) {
        let mut leaf_volume = 0u64;
        let mut free_volume = 0u64;
        for (i, leaf) in self.leaves.iter().enumerate() {
            if !leaf.alive {
                continue;
            }
            let v = self.leaf_box(i as u32).volume();
            leaf_volume += v;
            if leaf.state == AIR {
                free_volume += v;
                assert!(
                    leaf.region != NONE && self.regions[leaf.region as usize].alive,
                    "{where_}: free leaf {i} has no live region"
                );
            }
        }
        let world = u64::from(self.tree.n).pow(3);
        assert_eq!(
            leaf_volume, world,
            "{where_}: leaves cover {leaf_volume} of {world} cells, so they are not a partition"
        );
        let mut region_volume = 0u64;
        for (i, region) in self.regions.iter().enumerate() {
            if !region.alive {
                continue;
            }
            let members: u64 = region
                .leaves
                .iter()
                .map(|&l| self.leaf_box(l).volume())
                .sum();
            assert_eq!(
                region.bx.volume(),
                members,
                "{where_}: region {i} has box volume {} against member volume {members}: the \
                 greedy merge produced a region that is not an exactly tiled box, so it is not \
                 convexity-preserving",
                region.bx.volume()
            );
            region_volume += members;
        }
        assert_eq!(
            region_volume, free_volume,
            "{where_}: regions cover {region_volume} of {free_volume} free cells"
        );
    }
}

/// Hertel–Mehlhorn-inspired greedy merge, restricted to axis-aligned boxes.
///
/// `P-87`'s, verbatim in effect: the union of two boxes is convex **iff** they
/// are contiguous on one axis and their extents agree on the other two, so the
/// pass groups by cross-section and merges maximal contiguous runs, repeating
/// over the three axes until a whole sweep merges nothing. Deterministic,
/// because the groups come out of a `BTreeMap` in key order.
fn greedy_merge(items: &mut Vec<(Box3, Vec<u32>)>) {
    loop {
        let mut merged_any = false;
        for axis in 0..3 {
            let b = (axis + 1) % 3;
            let c = (axis + 2) % 3;
            let mut groups: BTreeMap<[u32; 4], Vec<usize>> = BTreeMap::new();
            for (i, it) in items.iter().enumerate() {
                groups
                    .entry([it.0.lo[b], it.0.hi[b], it.0.lo[c], it.0.hi[c]])
                    .or_default()
                    .push(i);
            }
            let mut dead = vec![false; items.len()];
            let mut any_here = false;
            for (_, mut group) in groups {
                if group.len() < 2 {
                    continue;
                }
                group.sort_unstable_by_key(|&i| items[i].0.lo[axis]);
                let mut head = group[0];
                for &i in &group[1..] {
                    if items[head].0.hi[axis] == items[i].0.lo[axis] {
                        items[head].0.hi[axis] = items[i].0.hi[axis];
                        let taken = core::mem::take(&mut items[i].1);
                        items[head].1.extend(taken);
                        dead[i] = true;
                        any_here = true;
                    } else {
                        head = i;
                    }
                }
            }
            if any_here {
                merged_any = true;
                let old = core::mem::take(items);
                let mut out = Vec::with_capacity(old.len());
                for (i, it) in old.into_iter().enumerate() {
                    if !dead[i] {
                        out.push(it);
                    }
                }
                *items = out;
            }
        }
        if !merged_any {
            return;
        }
    }
}

// ── the local repair, and the cell set it touched ───────────────────────────

/// What one edit's repair did, and where.
struct Touched {
    /// Cells covered by leaves removed, leaves added, regions dissolved and
    /// regions created — **`P-87`'s repair set**, expanded from nav cells to the
    /// world cells they cover, because a λ-membership flip is a per-cell event.
    stamp_all: Vec<u32>,
    /// The same, restricted to leaves. The shrunken set exists so the
    /// containment predicate is exercised in the direction that can fail.
    stamp_leaf: Vec<u32>,
    tag: u32,
    cells_all: u64,
    cells_leaf: u64,
    flipped: u64,
    pattern_cells: u64,
    repair_nav_cells: u64,
}

impl Touched {
    fn new(n: u32) -> Self {
        let count = (n as usize).pow(3);
        Self {
            stamp_all: vec![0; count],
            stamp_leaf: vec![0; count],
            tag: 0,
            cells_all: 0,
            cells_leaf: 0,
            flipped: 0,
            pattern_cells: 0,
            repair_nav_cells: 0,
        }
    }

    fn begin(&mut self) {
        self.tag += 1;
        self.cells_all = 0;
        self.cells_leaf = 0;
        self.flipped = 0;
        self.pattern_cells = 0;
        self.repair_nav_cells = 0;
    }
}

/// Mark every cell of `bx` in `stamp`, counting the ones that were not already
/// marked.
fn stamp_box(stamp: &mut [u32], tag: u32, n: u32, bx: Box3, count: &mut u64) {
    let n = n as usize;
    for z in bx.lo[2]..bx.hi[2] {
        for y in bx.lo[1]..bx.hi[1] {
            let row = (z as usize * n + y as usize) * n;
            for x in bx.lo[0]..bx.hi[0] {
                let i = row + x as usize;
                if stamp[i] != tag {
                    stamp[i] = tag;
                    *count += 1;
                }
            }
        }
    }
}

/// Apply one edit's new sign bits and repair the octree locally.
///
/// `P-87`'s repair, with the repair set recorded as a cell stamp so C2 can ask
/// whether a λ-flip landed inside it.
fn repair(
    nav: &mut Nav,
    signs: &mut Signs,
    lo: [u32; 3],
    ext: [u32; 3],
    new_solid: &[bool],
    flipped: &mut Vec<[u32; 3]>,
    touched: &mut Touched,
) {
    let n = nav.tree.n;
    let chunks_per_axis = n / CHUNK_CELLS;
    touched.begin();

    flipped.clear();
    for z in 0..ext[2] {
        for y in 0..ext[1] {
            for x in 0..ext[0] {
                let s = [lo[0] + x, lo[1] + y, lo[2] + z];
                let want = new_solid[((z * ext[1] + y) * ext[0] + x) as usize];
                if signs.solid(s) != want {
                    signs.set(s, want);
                    flipped.push(s);
                }
            }
        }
    }
    touched.flipped = flipped.len() as u64;

    let mut chunk_ids: Vec<u32> = Vec::new();
    let mut pattern: Vec<u32> = Vec::new();
    for &s in flipped.iter() {
        for dz in 0..2u32 {
            for dy in 0..2u32 {
                for dx in 0..2u32 {
                    let mut cell = [0u32; 3];
                    let mut ok = true;
                    for (a, d) in [dx, dy, dz].into_iter().enumerate() {
                        let v = s[a] + d;
                        if v == 0 || v > n {
                            ok = false;
                            break;
                        }
                        cell[a] = v - 1;
                    }
                    if !ok {
                        continue;
                    }
                    pattern.push(nav.tree.idx(0, cell) as u32);
                    let ch = [
                        cell[0] / CHUNK_CELLS,
                        cell[1] / CHUNK_CELLS,
                        cell[2] / CHUNK_CELLS,
                    ];
                    chunk_ids.push((ch[2] * chunks_per_axis + ch[1]) * chunks_per_axis + ch[0]);
                }
            }
        }
    }
    pattern.sort_unstable();
    pattern.dedup();
    touched.pattern_cells = pattern.len() as u64;
    chunk_ids.sort_unstable();
    chunk_ids.dedup();

    // Reclassify whole chunks: that is the granularity a chunked world hands its
    // consumers, and it is P-87's.
    let mut changed: Vec<u32> = Vec::new();
    for &ch in &chunk_ids {
        let cz = ch / (chunks_per_axis * chunks_per_axis);
        let cy = (ch / chunks_per_axis) % chunks_per_axis;
        let cx = ch % chunks_per_axis;
        let base = [cx * CHUNK_CELLS, cy * CHUNK_CELLS, cz * CHUNK_CELLS];
        for dz in 0..CHUNK_CELLS {
            for dy in 0..CHUNK_CELLS {
                for dx in 0..CHUNK_CELLS {
                    let cell = [base[0] + dx, base[1] + dy, base[2] + dz];
                    let st = signs.cell_state(cell);
                    let i = nav.tree.idx(0, cell);
                    if nav.tree.levels[0][i] != st {
                        nav.tree.levels[0][i] = st;
                        changed.push(i as u32);
                    }
                }
            }
        }
    }

    let depth = nav.tree.depth;
    let mut changed_levels: Vec<Vec<u32>> = vec![Vec::new(); depth + 1];
    changed_levels[0] = changed;
    for l in 1..=depth {
        let m = nav.tree.side(l) as usize;
        let mc = m * 2;
        let mut cand: Vec<u32> = Vec::with_capacity(changed_levels[l - 1].len());
        for &i in &changed_levels[l - 1] {
            let i = i as usize;
            let x = i % mc;
            let y = (i / mc) % mc;
            let z = i / (mc * mc);
            cand.push((((z / 2) * m + y / 2) * m + x / 2) as u32);
        }
        cand.sort_unstable();
        cand.dedup();
        let mut moved = Vec::with_capacity(cand.len());
        for &i in &cand {
            let i = i as usize;
            let x = i % m;
            let y = (i / m) % m;
            let z = i / (m * m);
            let mut kids = [MIXED; 8];
            for dz in 0..2 {
                for dy in 0..2 {
                    for dx in 0..2 {
                        kids[(dz * 2 + dy) * 2 + dx] = nav.tree.levels[l - 1]
                            [((2 * z + dz) * mc + (2 * y + dy)) * mc + (2 * x + dx)];
                    }
                }
            }
            let st = combine(kids);
            if nav.tree.levels[l][i] != st {
                nav.tree.levels[l][i] = st;
                moved.push(i as u32);
            }
        }
        changed_levels[l] = moved;
    }

    // The top-most nodes whose state moved: everything below them follows.
    let mut tops: Vec<(usize, [u32; 3])> = Vec::new();
    for l in 0..=depth {
        let m = nav.tree.side(l) as usize;
        for &i in &changed_levels[l] {
            let i = i as usize;
            let x = i % m;
            let y = (i / m) % m;
            let z = i / (m * m);
            if l < depth {
                let mp = m / 2;
                let pi = (((z / 2) * mp + y / 2) * mp + x / 2) as u32;
                if changed_levels[l + 1].binary_search(&pi).is_ok() {
                    continue;
                }
            }
            tops.push((l, [x as u32, y as u32, z as u32]));
        }
    }

    let mut old_leaves: Vec<u32> = Vec::new();
    for &(l, c) in &tops {
        nav.collect_under(l, c, &mut old_leaves);
    }
    old_leaves.sort_unstable();
    old_leaves.dedup();
    touched.repair_nav_cells += old_leaves.len() as u64;

    let mut doomed: Vec<u32> = Vec::new();
    for &l in &old_leaves {
        let r = nav.leaves[l as usize].region;
        if r != NONE {
            doomed.push(r);
        }
    }
    doomed.sort_unstable();
    doomed.dedup();
    touched.repair_nav_cells += doomed.len() as u64;

    let tag = touched.tag;
    for &l in &old_leaves {
        let bx = nav.leaf_box(l);
        stamp_box(&mut touched.stamp_all, tag, n, bx, &mut touched.cells_all);
        stamp_box(&mut touched.stamp_leaf, tag, n, bx, &mut touched.cells_leaf);
    }
    for &r in &doomed {
        let bx = nav.regions[r as usize].bx;
        stamp_box(&mut touched.stamp_all, tag, n, bx, &mut touched.cells_all);
    }

    for &l in &old_leaves {
        nav.kill_leaf(l);
    }
    let mut new_leaves: Vec<u32> = Vec::new();
    for &(l, c) in &tops {
        nav.derive_under(l, c, &mut new_leaves);
    }
    touched.repair_nav_cells += new_leaves.len() as u64;
    for &l in &new_leaves {
        let bx = nav.leaf_box(l);
        stamp_box(&mut touched.stamp_all, tag, n, bx, &mut touched.cells_all);
        stamp_box(&mut touched.stamp_leaf, tag, n, bx, &mut touched.cells_leaf);
    }

    let mut pool: Vec<u32> = new_leaves
        .iter()
        .copied()
        .filter(|&l| nav.leaves[l as usize].state == AIR)
        .collect();
    for &r in &doomed {
        if !nav.regions[r as usize].alive {
            continue;
        }
        for l in nav.regions[r as usize].leaves.clone() {
            if nav.leaves[l as usize].alive {
                pool.push(l);
            }
        }
        nav.kill_region(r);
    }
    pool.sort_unstable();
    pool.dedup();
    let mut made = Vec::new();
    nav.merge_pool(&pool, &mut made);
    touched.repair_nav_cells += made.len() as u64;
    for &r in &made {
        let bx = nav.regions[r as usize].bx;
        stamp_box(&mut touched.stamp_all, tag, n, bx, &mut touched.cells_all);
    }

    // Only now may killed ids be handed out again: within one repair a dissolved
    // region's member list still names them.
    let mut retired = core::mem::take(&mut nav.leaf_retire);
    nav.leaf_free.append(&mut retired);
    nav.leaf_retire = retired;
}

// ── the edits ───────────────────────────────────────────────────────────────

/// A shape removed from a field. Two of them, because a trace needs a wide dig
/// and a passage narrower than two cells.
#[derive(Clone, Copy)]
enum Cut {
    Ball(Sphere<f64>),
    Tube(Capsule<f64>),
}

impl Sdf for Cut {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Ball(s) => s.sample(p),
            Self::Tube(c) => c.sample(p),
        }
    }
}

impl Cut {
    /// The world box the shape can touch.
    fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        match self {
            Self::Ball(s) => (
                [0, 1, 2].map(|a| s.center[a] - s.radius),
                [0, 1, 2].map(|a| s.center[a] + s.radius),
            ),
            Self::Tube(c) => (
                [0, 1, 2].map(|a| c.a[a].min(c.b[a]) - c.radius),
                [0, 1, 2].map(|a| c.a[a].max(c.b[a]) + c.radius),
            ),
        }
    }
}

/// An axis-aligned cylinder along `x`, in world units. Both the fixture's
/// channels and the sampling strata are one of these.
#[derive(Clone, Copy)]
struct Channel {
    x0: f64,
    x1: f64,
    y: f64,
    z: f64,
    radius: f64,
}

impl Channel {
    fn volume(&self) -> f64 {
        std::f64::consts::PI * self.radius * self.radius * (self.x1 - self.x0)
    }
}

// ── sampling ────────────────────────────────────────────────────────────────

/// SplitMix64. A fixed seed per fixture, so the sample set is reproducible
/// (`M-36`): the accuracy fraction has to be a property of the bound rather than
/// of a lucky draw.
struct Rng(u64);

impl Rng {
    #[inline]
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    #[inline]
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }

    #[inline]
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

/// One evaluation point.
struct Sample {
    /// Cell coordinates in the coarse world.
    q: [f64; 3],
    /// True clearance in coarse cells.
    truth: f64,
    /// Air here, solid before the trace: the point is in dug material.
    dug: bool,
    /// The truth came from a closed form rather than a refined grid.
    known: bool,
    /// From the uniform stratum rather than a forced one.
    uniform: bool,
}

/// How a fixture's sample set is drawn.
enum SampleMode {
    /// Two cylinders, sampled in proportion to their volumes. On the slab this
    /// is the whole of the air, so the set is unbiased over it.
    Channels { wide: Channel, narrow: Channel },
    /// Uniform over the world, rejecting solid, plus a forced narrow stratum.
    WorldAir { narrow: Channel },
}

/// Where a fixture's truth comes from.
enum TruthMode {
    /// The field value is the exact distance to the solid, so the clearance is
    /// read straight off it. `M-346`'s fixture and its argument.
    Analytic,
    /// The same nearest-non-free query on a refined sign field.
    Refined { tree: Tree, factor: u32 },
}

// ── measurement ─────────────────────────────────────────────────────────────

/// One estimator's accuracy over one fixture's samples.
#[derive(Default)]
struct Accuracy {
    within: u64,
    within_half: u64,
    within_two: u64,
    over: u64,
    max_over: f64,
    max_under: f64,
    sum_err: f64,
    probe_over: u64,
    zero_bound: u64,
    uniform: u64,
    uniform_within: u64,
    narrow_within: u64,
    ns: u128,
    errors: Vec<f32>,
}

/// One estimator's λ accounting over one fixture's trace.
#[derive(Clone, Copy)]
struct Locality {
    /// Per λ: flips, flips outside the repair set, flips outside the leaf-only
    /// subset, and the worst single edit's flip count.
    flips: [u64; LAMBDAS.len()],
    outside_all: [u64; LAMBDAS.len()],
    outside_leaf: [u64; LAMBDAS.len()],
    worst_edit_flips: [u64; LAMBDAS.len()],
    worst_edit_changed: [u64; LAMBDAS.len()],
    worst_edit_index: [usize; LAMBDAS.len()],
}

impl Locality {
    fn new() -> Self {
        Self {
            flips: [0; LAMBDAS.len()],
            outside_all: [0; LAMBDAS.len()],
            outside_leaf: [0; LAMBDAS.len()],
            worst_edit_flips: [0; LAMBDAS.len()],
            worst_edit_changed: [1; LAMBDAS.len()],
            worst_edit_index: [0; LAMBDAS.len()],
        }
    }
}

/// Everything one fixture measured, for one estimator.
struct Out {
    fixture: String,
    estimator: &'static str,
    field: &'static str,
    world: u32,
    depth: usize,
    truth_source: &'static str,
    refine: u32,
    samples: u64,
    known: u64,
    narrow: u64,
    narrow_within: u64,
    acc: Accuracy,
    changed_samples: u64,
    lambda: f64,
    flips: u64,
    outside_all: u64,
    outside_leaf: u64,
    worst_edit_factor: f64,
    worst_edit_index: usize,
    repair_set_cells: u64,
    repair_nav_cells: u64,
    world_cells: u64,
    free_leaves: u64,
    regions: u64,
    all_leaves: u64,
    free_cells: u64,
    edits: usize,
    wide_edits: usize,
    narrow_edits: usize,
    loc: Locality,
}

/// Fill `out` with the region-box bound at every cell centre.
fn bounds_region(nav: &Nav, out: &mut [f32]) {
    out.fill(0.0);
    let n = nav.tree.n as usize;
    for region in &nav.regions {
        if !region.alive {
            continue;
        }
        let bx = region.bx;
        for z in bx.lo[2]..bx.hi[2] {
            let dz = f64::from(z - bx.lo[2] + 1) - 0.5;
            let uz = f64::from(bx.hi[2] - z) - 0.5;
            let mz = dz.min(uz);
            for y in bx.lo[1]..bx.hi[1] {
                let dy = f64::from(y - bx.lo[1] + 1) - 0.5;
                let uy = f64::from(bx.hi[1] - y) - 0.5;
                let mzy = mz.min(dy).min(uy);
                let row = (z as usize * n + y as usize) * n;
                for x in bx.lo[0]..bx.hi[0] {
                    let dx = f64::from(x - bx.lo[0] + 1) - 0.5;
                    let ux = f64::from(bx.hi[0] - x) - 0.5;
                    out[row + x as usize] = mzy.min(dx).min(ux) as f32;
                }
            }
        }
    }
}

/// Fill `out` with the nearest-non-free bound at every cell centre.
fn bounds_ball(tree: &Tree, out: &mut [f32], stack: &mut Vec<(f64, u8, [u32; 3])>) {
    let n = tree.n;
    let ns = n as usize;
    let fn_ = f64::from(n);
    for z in 0..n {
        for y in 0..n {
            let row = (z as usize * ns + y as usize) * ns;
            for x in 0..n {
                let i = row + x as usize;
                if tree.levels[0][i] != AIR {
                    out[i] = 0.0;
                    continue;
                }
                let q = [f64::from(x) + 0.5, f64::from(y) + 0.5, f64::from(z) + 0.5];
                let wall = q
                    .iter()
                    .map(|&v| v.min(fn_ - v))
                    .fold(f64::INFINITY, f64::min);
                out[i] = nearest_non_free(tree, q, wall, stack) as f32;
            }
        }
    }
}

/// Count λ-membership flips between two bound arrays and locate them against the
/// repair set.
fn accumulate_flips(
    prev: &[f32],
    cur: &[f32],
    touched: &Touched,
    edit: usize,
    changed: u64,
    loc: &mut Locality,
) {
    for (li, &lambda) in LAMBDAS.iter().enumerate() {
        let l = lambda as f32;
        let mut flips = 0u64;
        let mut out_all = 0u64;
        let mut out_leaf = 0u64;
        for i in 0..prev.len() {
            if (prev[i] >= l) == (cur[i] >= l) {
                continue;
            }
            flips += 1;
            if touched.stamp_all[i] != touched.tag {
                out_all += 1;
            }
            if touched.stamp_leaf[i] != touched.tag {
                out_leaf += 1;
            }
        }
        loc.flips[li] += flips;
        loc.outside_all[li] += out_all;
        loc.outside_leaf[li] += out_leaf;
        if flips > loc.worst_edit_flips[li] {
            loc.worst_edit_flips[li] = flips;
            loc.worst_edit_changed[li] = changed.max(1);
            loc.worst_edit_index[li] = edit;
        }
    }
}

/// Score one estimator's bound against the truth over a fixture's samples.
fn score(
    samples: &[Sample],
    bound_of: &mut dyn FnMut(&Sample) -> f64,
    other: Option<&[f64]>,
    acc: &mut Accuracy,
) -> Vec<f64> {
    let mut bounds = Vec::with_capacity(samples.len());
    let t = Instant::now();
    for s in samples {
        bounds.push(bound_of(s));
    }
    acc.ns = t.elapsed().as_nanos();
    acc.errors.reserve(samples.len());
    for (i, s) in samples.iter().enumerate() {
        let b = bounds[i];
        if let Some(tight) = other {
            assert!(
                b <= tight[i] + 1e-9,
                "the merged-region bound {b} exceeds the free-cell-union bound {} at cell \
                 ({:.3}, {:.3}, {:.3}): a region that is not inside the free set",
                tight[i],
                s.q[0],
                s.q[1],
                s.q[2]
            );
        }
        let err = s.truth - b;
        acc.errors.push(err as f32);
        acc.sum_err += err;
        if err.abs() <= 1.0 {
            acc.within += 1;
            if s.uniform {
                acc.uniform_within += 1;
            }
            if s.dug && s.truth < 1.0 {
                acc.narrow_within += 1;
            }
        }
        if err.abs() <= 0.5 {
            acc.within_half += 1;
        }
        if err.abs() <= 2.0 {
            acc.within_two += 1;
        }
        if s.uniform {
            acc.uniform += 1;
        }
        if b > s.truth {
            acc.over += 1;
            acc.max_over = acc.max_over.max(b - s.truth);
        } else {
            acc.max_under = acc.max_under.max(s.truth - b);
        }
        // M-44: the identical comparison against a bound half a voxel too
        // optimistic — the error of measuring to an obstacle cell's centre
        // instead of to its box.
        if b + 0.5 > s.truth {
            acc.probe_over += 1;
        }
        if b <= 0.0 {
            acc.zero_bound += 1;
        }
    }
    bounds
}

/// A percentile of an unsorted error list.
fn percentile(errors: &mut [f32], p: f64) -> f64 {
    if errors.is_empty() {
        return 0.0;
    }
    errors.sort_unstable_by(f32::total_cmp);
    let k = ((errors.len() - 1) as f64 * p).round() as usize;
    f64::from(errors[k])
}

// ── the fixture ─────────────────────────────────────────────────────────────

/// Build, dig, and measure one fixture; returns one row per estimator.
fn run_fixture<F: Sdf<Scalar = f64>>(
    base: &F,
    fixture: String,
    field_name: &'static str,
    n: u32,
    base_cuts: Vec<Brush<Cut>>,
    trace: Vec<Brush<Cut>>,
    truth_mode: TruthMode,
    sample_mode: SampleMode,
    seed: u64,
    wide_edits: usize,
    narrow_edits: usize,
) -> [Out; 2] {
    let h = EXTENT / f64::from(n);
    let layout = ChunkLayout::<f64>::new(CHUNK_CELLS, h, [ORIGIN; 3]).expect("layout");
    let sdim = n + 1;
    let cells = (n as usize).pow(3);

    let mut all_cuts = base_cuts.clone();
    all_cuts.extend_from_slice(&trace);
    let pristine = BrushStack {
        base,
        brushes: &base_cuts[..],
    };
    let final_field = BrushStack {
        base,
        brushes: &all_cuts[..],
    };

    // ── the pristine world, and the structure over it ───────────────────────
    let mut signs = Signs::new([sdim; 3]);
    signs.fill(&pristine, &layout);
    let mut nav = Nav::new(n);
    nav.tree.rebuild(&signs);
    nav.build();
    nav.audit("after build");

    let free_leaves = nav
        .leaves
        .iter()
        .filter(|l| l.alive && l.state == AIR)
        .count() as u64;
    let free_cells: u64 = nav
        .leaves
        .iter()
        .enumerate()
        .filter(|(_, l)| l.alive && l.state == AIR)
        .map(|(i, _)| nav.leaf_box(i as u32).volume())
        .sum();
    let all_leaves = nav.live_leaves;
    let static_regions = nav.live_regions;
    let depth = nav.tree.depth;

    // ── the trace, with both estimators' λ accounting ───────────────────────
    let mut touched = Touched::new(n);
    let mut flipped: Vec<[u32; 3]> = Vec::new();
    let mut stack: Vec<(f64, u8, [u32; 3])> = Vec::new();
    let mut dirty = DirtySet::new();

    let mut prev_region = vec![0f32; cells];
    let mut cur_region = vec![0f32; cells];
    let mut prev_ball = vec![0f32; cells];
    let mut cur_ball = vec![0f32; cells];
    bounds_region(&nav, &mut prev_region);
    bounds_ball(&nav.tree, &mut prev_ball, &mut stack);

    let mut loc_region = Locality::new();
    let mut loc_ball = Locality::new();
    let mut changed_samples = 0u64;
    let mut repair_set_cells = 0u64;
    let mut repair_nav_cells = 0u64;

    for step in 0..trace.len() {
        let k = base_cuts.len() + step;
        let before = BrushStack {
            base,
            brushes: &all_cuts[..k],
        };
        let after = BrushStack {
            base,
            brushes: &all_cuts[..=k],
        };
        let (lo_world, hi_world) = all_cuts[k].shape.aabb();
        let lo_i = layout
            .cell_of(lo_world)
            .map(|v| (v - 1).clamp(0, i64::from(n) - 1));
        let hi_i = layout
            .cell_of(hi_world)
            .map(|v| (v + 1).clamp(0, i64::from(n) - 1));

        // Untimed control: the crate's own instrument over the same box.
        let report = mark_edit(&layout, &before, &after, lo_i, hi_i, &mut dirty).expect("mark");
        dirty.clear();

        let lo = [0, 1, 2].map(|a| lo_i[a] as u32);
        let ext = [0, 1, 2].map(|a| (hi_i[a] - lo_i[a] + 2) as u32);
        let mut new_solid = vec![false; (ext[0] * ext[1] * ext[2]) as usize];
        for z in 0..ext[2] {
            for y in 0..ext[1] {
                for x in 0..ext[0] {
                    let s = [
                        i64::from(lo[0] + x),
                        i64::from(lo[1] + y),
                        i64::from(lo[2] + z),
                    ];
                    let v = after.sample(layout.world_of_sample(s));
                    new_solid[((z * ext[1] + y) * ext[0] + x) as usize] = v < 0.0;
                }
            }
        }

        repair(
            &mut nav,
            &mut signs,
            lo,
            ext,
            &new_solid,
            &mut flipped,
            &mut touched,
        );

        assert_eq!(
            touched.pattern_cells, report.sign_changed_cells,
            "{fixture}, edit {step}: the harness found {} cells incident to a flipped sign bit \
             and `mark_edit` found {} over the same box, so one of the two is looking at the \
             wrong region",
            touched.pattern_cells, report.sign_changed_cells
        );
        nav.audit("after edit");

        changed_samples += touched.flipped;
        repair_set_cells += touched.cells_all;
        repair_nav_cells += touched.repair_nav_cells;

        bounds_region(&nav, &mut cur_region);
        bounds_ball(&nav.tree, &mut cur_ball, &mut stack);
        accumulate_flips(
            &prev_region,
            &cur_region,
            &touched,
            step,
            touched.flipped,
            &mut loc_region,
        );
        accumulate_flips(
            &prev_ball,
            &cur_ball,
            &touched,
            step,
            touched.flipped,
            &mut loc_ball,
        );
        core::mem::swap(&mut prev_region, &mut cur_region);
        core::mem::swap(&mut prev_ball, &mut cur_ball);
    }

    // ── the sample set, on the final world ──────────────────────────────────
    let mut rng = Rng(seed);
    let mut samples: Vec<Sample> = Vec::with_capacity(SAMPLES_PER_FIXTURE);
    let to_cell = |p: [f64; 3]| [0, 1, 2].map(|a| (p[a] - ORIGIN) / h);
    let fnn = f64::from(n);
    let wall_of = |q: [f64; 3]| {
        q.iter()
            .map(|&v| v.min(fnn - v))
            .fold(f64::INFINITY, f64::min)
    };
    let mut truth_of = |q: [f64; 3], p: [f64; 3]| -> f64 {
        let wall = wall_of(q);
        match &truth_mode {
            TruthMode::Analytic => (final_field.sample(p) / h).min(wall),
            TruthMode::Refined { tree, factor } => {
                let f = f64::from(*factor);
                let qf = [q[0] * f, q[1] * f, q[2] * f];
                nearest_non_free(tree, qf, wall * f, &mut stack) / f
            }
        }
    };

    match &sample_mode {
        SampleMode::Channels { wide, narrow } => {
            let vw = wide.volume();
            let vn = narrow.volume();
            let p_narrow = vn / (vw + vn);
            for _ in 0..SAMPLES_PER_FIXTURE {
                let in_narrow = rng.unit() < p_narrow;
                let tube = if in_narrow { narrow } else { wide };
                let p = point_in_tube(&mut rng, tube);
                let q = to_cell(p);
                let dug = pristine.sample(p) <= 0.0;
                samples.push(Sample {
                    q,
                    truth: truth_of(q, p),
                    dug,
                    known: true,
                    uniform: !in_narrow,
                });
            }
        }
        SampleMode::WorldAir { narrow } => {
            let forced = (SAMPLES_PER_FIXTURE as f64 * NARROW_SHARE) as usize;
            for _ in 0..forced {
                let p = point_in_tube(&mut rng, narrow);
                let q = to_cell(p);
                if final_field.sample(p) <= 0.0 {
                    continue;
                }
                let dug = pristine.sample(p) <= 0.0;
                samples.push(Sample {
                    q,
                    truth: truth_of(q, p),
                    dug,
                    known: false,
                    uniform: false,
                });
            }
            let mut attempts = 0u64;
            while samples.len() < SAMPLES_PER_FIXTURE {
                attempts += 1;
                assert!(
                    attempts < 200 * SAMPLES_PER_FIXTURE as u64,
                    "{fixture}: rejection sampling cannot find air; the world is solid"
                );
                let p = [
                    rng.range(ORIGIN, ORIGIN + EXTENT),
                    rng.range(ORIGIN, ORIGIN + EXTENT),
                    rng.range(ORIGIN, ORIGIN + EXTENT),
                ];
                if final_field.sample(p) <= 0.0 {
                    continue;
                }
                let q = to_cell(p);
                let dug = pristine.sample(p) <= 0.0;
                samples.push(Sample {
                    q,
                    truth: truth_of(q, p),
                    dug,
                    known: false,
                    uniform: true,
                });
            }
        }
    }

    let known = samples.iter().filter(|s| s.known).count() as u64;
    let narrow = samples.iter().filter(|s| s.dug && s.truth < 1.0).count() as u64;

    // ── score both estimators over the same samples ─────────────────────────
    let mut acc_ball = Accuracy::default();
    let ball_bounds = {
        let nav_tree = &nav.tree;
        let stack = &mut stack;
        let mut f = |s: &Sample| {
            let wall = wall_of(s.q);
            nearest_non_free(nav_tree, s.q, wall, stack)
        };
        score(&samples, &mut f, None, &mut acc_ball)
    };
    let mut acc_region = Accuracy::default();
    {
        let navr = &nav;
        let mut f = |s: &Sample| {
            let c = [0, 1, 2].map(|a| (s.q[a].floor() as i64).clamp(0, i64::from(n) - 1) as u32);
            let leaf = navr.leaf_of_cell(c);
            if leaf == NONE || navr.leaves[leaf as usize].state != AIR {
                return 0.0;
            }
            let r = navr.leaves[leaf as usize].region;
            if r == NONE {
                return 0.0;
            }
            navr.regions[r as usize].bx.inradius_at(s.q)
        };
        score(&samples, &mut f, Some(&ball_bounds), &mut acc_region);
    }

    let world_cells = u64::from(n).pow(3);
    let truth_source = match &truth_mode {
        TruthMode::Analytic => "analytic",
        TruthMode::Refined { .. } => "refined",
    };
    let edits = trace.len();

    let finish = |estimator: &'static str, acc: Accuracy, loc: Locality| -> Out {
        // The registered columns carry the worst creature in the sweep. "Worst"
        // is ordered by `flips_outside_repair_set` first and flip count second,
        // because C2 is falsified *either* by the factor *or* by a flip set that
        // is not a subset of the repair set, and a λ that breaks the subset half
        // has to reach the row rather than be averaged away. Every λ's numbers
        // are on the row anyway, as `flip_factor_l*` and `outside_l*`.
        let (li, _) = loc
            .flips
            .iter()
            .enumerate()
            .max_by_key(|&(i, f)| (loc.outside_all[i], *f))
            .expect("non-empty sweep");
        Out {
            fixture: fixture.clone(),
            estimator,
            field: field_name,
            world: n,
            depth,
            truth_source,
            refine: match &truth_mode {
                TruthMode::Analytic => 1,
                TruthMode::Refined { factor, .. } => *factor,
            },
            samples: samples.len() as u64,
            known,
            narrow,
            narrow_within: acc.narrow_within,
            acc,
            changed_samples,
            lambda: LAMBDAS[li],
            flips: loc.flips[li],
            outside_all: loc.outside_all[li],
            outside_leaf: loc.outside_leaf[li],
            worst_edit_factor: loc.worst_edit_flips[li] as f64 / loc.worst_edit_changed[li] as f64,
            worst_edit_index: loc.worst_edit_index[li],
            repair_set_cells,
            repair_nav_cells,
            world_cells,
            free_leaves,
            regions: static_regions,
            all_leaves,
            free_cells,
            edits,
            wide_edits,
            narrow_edits,
            loc,
        }
    };

    [
        finish("region", acc_region, loc_region),
        finish("ball", acc_ball, loc_ball),
    ]
}

/// Uniform point inside an `x`-aligned cylinder.
fn point_in_tube(rng: &mut Rng, tube: &Channel) -> [f64; 3] {
    let x = rng.range(tube.x0, tube.x1);
    loop {
        let dy = rng.range(-tube.radius, tube.radius);
        let dz = rng.range(-tube.radius, tube.radius);
        if dy * dy + dz * dz <= tube.radius * tube.radius {
            return [x, tube.y + dy, tube.z + dz];
        }
    }
}

// ── the fixtures ────────────────────────────────────────────────────────────

/// Split an `x`-aligned tube into `count` capsule cuts that meet end to end.
fn tube_cuts(tube: &Channel, count: usize) -> Vec<Brush<Cut>> {
    (0..count)
        .map(|i| {
            let a = tube.x0 + (tube.x1 - tube.x0) * (i as f64) / count as f64;
            let b = tube.x0 + (tube.x1 - tube.x0) * ((i + 1) as f64) / count as f64;
            Brush::subtract(Cut::Tube(Capsule {
                a: [a, tube.y, tube.z],
                b: [b, tube.y, tube.z],
                radius: tube.radius,
            }))
        })
        .collect()
}

/// `M-346`'s drilled slab, plus a narrow cut dug through solid beside it.
///
/// The wide channel is `M-346`'s exact fixture and is part of the pristine
/// world; the narrow channel is the **edit trace**, because the registered
/// vacuity control asks for points in *dug* passages.
fn slab_fixture(radius_cells: u32) -> [Out; 2] {
    let n = SLAB_CELLS;
    let h = EXTENT / f64::from(n);
    let base = BoxExact::<f64> {
        center: [0.0; 3],
        half_extents: [4.0; 3],
    };
    let (lo, hi) = BoxExact::<f64>::canonical().domain();
    assert!(
        (lo[0] - ORIGIN).abs() < 1e-12 && (hi[0] - ORIGIN - EXTENT).abs() < 1e-12,
        "the slab fixture assumes M-346's [-2, 2] domain"
    );

    let wide = Channel {
        x0: ORIGIN,
        x1: ORIGIN + EXTENT,
        y: 0.0,
        z: 0.0,
        radius: h * f64::from(radius_cells),
    };
    // M-346's capsule runs from -4 to 4, so inside the domain it is a straight
    // cylinder rather than a capped one: the truth is `r - rho` with no end cap.
    let base_cuts = vec![Brush::subtract(Cut::Tube(Capsule {
        a: [-4.0, wide.y, wide.z],
        b: [4.0, wide.y, wide.z],
        radius: wide.radius,
    }))];

    // The narrow channel: axis down a line of cell centres, far enough away that
    // the escape ray from either channel never enters the other, so each
    // channel's clearance is exactly its own `r - rho`.
    let narrow = Channel {
        x0: ORIGIN + 4.0 * h,
        x1: ORIGIN + EXTENT - 4.0 * h,
        y: 1.0 + 0.5 * h,
        z: 0.5 * h,
        radius: h * NARROW_CELLS,
    };
    let separation = ((narrow.y - wide.y).powi(2) + (narrow.z - wide.z).powi(2)).sqrt();
    assert!(
        separation > wide.radius + narrow.radius + 4.0 * h,
        "the two channels are {separation} apart and the analytic truth needs them clear of \
         each other by more than their radii"
    );

    let trace = tube_cuts(&narrow, SLAB_NARROW_EDITS);
    run_fixture(
        &base,
        format!("slab_r{radius_cells}"),
        "slab",
        n,
        base_cuts,
        trace,
        TruthMode::Analytic,
        SampleMode::Channels { wide, narrow },
        0x5157_0000 + u64::from(radius_cells),
        0,
        SLAB_NARROW_EDITS,
    )
}

/// `P-72`'s dig path: straight across `x` through the middle of the world, with
/// the height probed **per edit at that edit's own `x`**. `P-87`'s, inherited.
fn wide_path<F: Sdf<Scalar = f64>>(field: &F, count: usize) -> Vec<[f64; 3]> {
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
    (0..count)
        .map(|i| {
            let t = (i as f64 + 0.5) / count as f64;
            let x = ORIGIN + EXTENT * t;
            [x, surface_y(x), mid]
        })
        .collect()
}

/// The cell-centre line through the most solid material, searched rather than
/// assumed.
///
/// `P-87`'s rule, applied to a different fixture: a narrow cut down a line that
/// is already air removes nothing, and the registered vacuity column would then
/// be a zero that could not have been non-zero.
fn narrow_line<F: Sdf<Scalar = f64>>(
    field: &F,
    layout: &ChunkLayout<f64>,
    n: u32,
    x0: f64,
    x1: f64,
) -> (f64, f64) {
    let h = EXTENT / f64::from(n);
    let probes = 64u32;
    let mut best: Option<(u32, u32, u32)> = None;
    let step = (n / 32).max(1);
    let mut j = 2;
    while j + 2 < n {
        let mut k = 2;
        while k + 2 < n {
            let y = layout.world_of_sample([0, i64::from(j), 0])[1] + 0.5 * h;
            let z = layout.world_of_sample([0, 0, i64::from(k)])[2] + 0.5 * h;
            let mut solid = 0u32;
            for i in 0..probes {
                let x = x0 + (x1 - x0) * (f64::from(i) + 0.5) / f64::from(probes);
                // The four sample lines a 0.9-cell cut actually reaches.
                let all = [(0.5, 0.5), (0.5, -0.5), (-0.5, 0.5), (-0.5, -0.5)]
                    .into_iter()
                    .all(|(dy, dz)| field.sample([x, y + dy * h, z + dz * h]) < 0.0);
                if all {
                    solid += 1;
                }
            }
            if best.is_none_or(|(s, _, _)| solid > s) {
                best = Some((solid, j, k));
            }
            k += step;
        }
        j += step;
    }
    let (solid, j, k) = best.expect("a candidate line");
    assert!(
        solid * 2 > probes,
        "the best cell-centre line is solid at only {solid} of {probes} probes, so the narrow \
         cut would dig mostly in air and `narrow_passage_samples` would be a vacuous count"
    );
    (
        layout.world_of_sample([0, i64::from(j), 0])[1] + 0.5 * h,
        layout.world_of_sample([0, 0, i64::from(k)])[2] + 0.5 * h,
    )
}

/// A reference field, dug with `P-87`'s wide brush and then cracked open with a
/// passage one cell across.
fn reference_fixture<F: Sdf<Scalar = f64> + Copy>(
    field: F,
    name: &'static str,
    n: u32,
    seed: u64,
) -> [Out; 2] {
    let h = EXTENT / f64::from(n);
    let layout = ChunkLayout::<f64>::new(CHUNK_CELLS, h, [ORIGIN; 3]).expect("layout");

    let wide: Vec<Brush<Cut>> = wide_path(&field, WIDE_EDITS)
        .into_iter()
        .map(|c| {
            Brush::subtract(Cut::Ball(Sphere {
                center: c,
                radius: h * WIDE_CELLS,
            }))
        })
        .collect();

    let x0 = ORIGIN + 4.0 * h;
    let x1 = ORIGIN + EXTENT - 4.0 * h;
    // Searched on the field the crack will actually be cut into: after the wide
    // dig, not before it.
    let dug_so_far = BrushStack {
        base: &field,
        brushes: &wide[..],
    };
    let (y, z) = narrow_line(&dug_so_far, &layout, n, x0, x1);
    let narrow = Channel {
        x0,
        x1,
        y,
        z,
        radius: h * NARROW_CELLS,
    };
    let mut trace = wide;
    trace.extend(tube_cuts(&narrow, NARROW_EDITS));

    // The refined truth: the same edited field on a 2x grid.
    let fine_n = n * REFINE;
    let fine_layout = ChunkLayout::<f64>::new(CHUNK_CELLS, EXTENT / f64::from(fine_n), [ORIGIN; 3])
        .expect("fine layout");
    let mut fine_signs = Signs::new([fine_n + 1; 3]);
    let full = BrushStack {
        base: &field,
        brushes: &trace[..],
    };
    fine_signs.fill(&full, &fine_layout);
    let mut fine_tree = Tree::new(fine_n);
    fine_tree.rebuild(&fine_signs);
    drop(fine_signs);

    run_fixture(
        &field,
        format!("{name}_{n}"),
        name,
        n,
        Vec::new(),
        trace,
        TruthMode::Refined {
            tree: fine_tree,
            factor: REFINE,
        },
        SampleMode::WorldAir { narrow },
        seed,
        WIDE_EDITS,
        NARROW_EDITS,
    )
}

/// The clock, on the row. `M-280`: on a governed CPU a nanosecond is not a unit.
fn clock_mhz() -> f64 {
    let text = std::fs::read_to_string("/proc/cpuinfo").unwrap_or_default();
    let mut best = 0f64;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("cpu MHz") else {
            continue;
        };
        let Some(value) = rest.split(':').nth(1) else {
            continue;
        };
        let Ok(mhz) = value.trim().parse::<f64>() else {
            continue;
        };
        if mhz > best {
            best = mhz;
        }
    }
    best
}

type Row = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|a| a == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-88");
    let clock = clock_mhz();

    println!(
        "SHARE: the registration says 'this is an accuracy and a rate, not a ratio of a total', \
         and that is correct. C1 is an accuracy rate whose denominator IS the sample set (share \
         1.0); C2 is a ratio of two directly counted sets plus a subset test; C3 is a one-sided \
         count over 10^6 samples. None is a speedup, so the 1/(1 - s + s/factor) ceiling of \
         \u{2717}51 cannot bind. Reachability, both directions, before the run: C1 upward by an \
         exact bound, downward by a one-cell-thick region inside an eight-cell channel (7.5 \
         cells short); C2 upward by dissolving one large region, downward by an edit that \
         touches a handful of leaves; C3 downward by a sound bound, and upward is DEMONSTRATED \
         on the same samples by the half-voxel probe rather than argued.\n"
    );

    let mut outs: Vec<Out> = Vec::new();
    for r in SLAB_RADII {
        outs.extend(slab_fixture(r));
    }
    for n in REF_WORLDS {
        outs.extend(reference_fixture(
            FbmTerrain::<f64>::canonical(),
            "fbm_terrain",
            n,
            0xFB00_0000 + u64::from(n),
        ));
        outs.extend(reference_fixture(
            Gyroid::<f64>::canonical(),
            "gyroid",
            n,
            0x6900_0000 + u64::from(n),
        ));
    }

    println!(
        "{:>14} {:>8} {:>9} {:>8} {:>8} {:>7} {:>8} {:>7} {:>7} {:>8} {:>8}",
        "fixture",
        "est",
        "samples",
        "narrow",
        "within",
        "frac",
        "over",
        "lambda",
        "flips",
        "factor",
        "outside"
    );
    for o in &outs {
        println!(
            "{:>14} {:>8} {:>9} {:>8} {:>8} {:>7.4} {:>8} {:>7.1} {:>7} {:>8.3} {:>8}",
            o.fixture,
            o.estimator,
            o.samples,
            o.narrow,
            o.acc.within,
            o.acc.within as f64 / o.samples as f64,
            o.acc.over,
            o.lambda,
            o.flips,
            o.flips as f64 / o.changed_samples.max(1) as f64,
            o.outside_all
        );
    }

    // C1 is registered "across M-346's fixtures", so its verdict is the pooled
    // rate over the analytic rows and nothing else. Pooled here rather than in
    // a reader's head.
    for est in ["region", "ball"] {
        let (w, s): (u64, u64) = outs
            .iter()
            .filter(|o| o.estimator == est && o.truth_source == "analytic")
            .fold((0, 0), |(w, s), o| (w + o.acc.within, s + o.samples));
        println!(
            "C1 pooled over M-346's fixtures, {est}: {w} of {s} within one voxel = {:.6} \
             against a bar of 0.90",
            w as f64 / s as f64
        );
    }

    // ── the registered vacuity control, both halves, asserted ───────────────
    let total_known: u64 = outs
        .iter()
        .filter(|o| o.estimator == "region")
        .map(|o| o.known)
        .sum();
    let total_narrow: u64 = outs
        .iter()
        .filter(|o| o.estimator == "region")
        .map(|o| o.narrow)
        .sum();
    assert!(
        total_known > 0,
        "VACUOUS: not one sample has an analytically known clearance, so C1 has no fixture at \
         all — M-346's slab did not reach the sample set"
    );
    assert!(
        total_narrow > 0,
        "VACUOUS: not one sample lies in a dug passage narrower than two cells, so the bound is \
         tested only in open space, which is nowhere interesting (M-44)"
    );
    for o in &outs {
        assert!(
            o.known > 0 || o.truth_source != "analytic",
            "VACUOUS: {} claims an analytic truth and reports no known-clearance sample",
            o.fixture
        );
        assert!(
            o.narrow > 0,
            "VACUOUS: {} has no sample in a dug passage narrower than two cells. The narrow cut \
             is {} edits of radius {NARROW_CELLS} cells down a cell-centre line and it removed \
             nothing that the sample set found; the fixture cannot show the effect and must not \
             be scored (M-44).",
            o.fixture,
            o.narrow_edits
        );
        // M-44 for C3's one-sided zero: the identical comparison, half a voxel
        // more optimistic, over the identical samples.
        assert!(
            o.acc.probe_over > 0,
            "VACUOUS: on {}/{} a bound half a voxel too optimistic still never overstates the \
             truth, so `overstatements = {}` is a zero that could not have been non-zero",
            o.fixture,
            o.estimator,
            o.acc.over
        );
        // M-44 for C2's zero: there were flips, and most of the world is outside
        // the repair set, so a flip could have landed outside it.
        assert!(
            o.flips > 0,
            "VACUOUS: {}/{} recorded no lambda-membership flip at any of {} creature sizes over \
             {} edits, so `flips_outside_repair_set` is a zero over an empty set",
            o.fixture,
            o.estimator,
            LAMBDAS.len(),
            o.edits
        );
        assert!(
            o.repair_set_cells < o.world_cells,
            "VACUOUS: {}/{} stamped {} cells over {} edits against a world of {}, so the repair \
             set is not a strict subset and the containment test cannot fail",
            o.fixture,
            o.estimator,
            o.repair_set_cells,
            o.edits,
            o.world_cells
        );
    }
    // And the containment predicate is demonstrated firing, on this run's own
    // data, rather than argued to be capable of it.
    let ball_outside: u64 = outs
        .iter()
        .filter(|o| o.estimator == "ball")
        .map(|o| o.outside_all)
        .sum();
    assert!(
        ball_outside > 0,
        "VACUOUS: the `flips_outside_repair_set` predicate never returned non-zero anywhere in \
         this run, on either estimator, so a zero on the region rows proves nothing about the \
         predicate (M-44)"
    );
    for est in ["region", "ball"] {
        let total: u64 = outs
            .iter()
            .filter(|o| o.estimator == est)
            .map(|o| o.samples)
            .sum();
        assert!(
            total >= 1_000_000,
            "C3 is registered over 10^6 samples and the `{est}` estimator saw only {total}"
        );
    }

    // ── rows ────────────────────────────────────────────────────────────────
    let mut rows: Vec<Row> = Vec::new();
    for o in &mut outs {
        let samples = o.samples as f64;
        let within_fraction = o.acc.within as f64 / samples;
        let flip_factor = o.flips as f64 / o.changed_samples.max(1) as f64;
        let c1 = within_fraction >= 0.90;
        let c2 = flip_factor <= 4.0 && o.outside_all == 0;
        let c3 = o.acc.over == 0;
        let mean_err = o.acc.sum_err / samples;
        let median = percentile(&mut o.acc.errors, 0.50);
        let p90 = percentile(&mut o.acc.errors, 0.90);
        let mut row: Row = vec![
            ("fixture", o.fixture.clone()),
            ("samples", o.samples.to_string()),
            ("known_clearance_samples", o.known.to_string()),
            ("narrow_passage_samples", o.narrow.to_string()),
            ("within_one_voxel", o.acc.within.to_string()),
            ("within_fraction", format!("{within_fraction:.6}")),
            ("overstatements", o.acc.over.to_string()),
            ("max_overstatement_cells", format!("{:.9}", o.acc.max_over)),
            ("lambda_flips", o.flips.to_string()),
            ("changed_samples", o.changed_samples.to_string()),
            ("flip_factor", format!("{flip_factor:.4}")),
            ("flips_outside_repair_set", o.outside_all.to_string()),
            ("c1_holds", c1.to_string()),
            ("c2_holds", c2.to_string()),
            ("c3_holds", c3.to_string()),
            // ── extras ──────────────────────────────────────────────────────
            ("estimator", o.estimator.to_string()),
            ("field", o.field.to_string()),
            ("world_cells", o.world.to_string()),
            ("octree_depth", o.depth.to_string()),
            ("truth_source", o.truth_source.to_string()),
            ("refine_factor", o.refine.to_string()),
            (
                "cell_size_world",
                format!("{:.8}", EXTENT / f64::from(o.world)),
            ),
            ("edits", o.edits.to_string()),
            ("wide_edits", o.wide_edits.to_string()),
            ("narrow_edits", o.narrow_edits.to_string()),
            ("wide_brush_cells", format!("{WIDE_CELLS:.1}")),
            ("narrow_cut_cells", format!("{NARROW_CELLS:.2}")),
            ("lambda_cells", format!("{:.1}", o.lambda)),
            ("within_half_voxel", o.acc.within_half.to_string()),
            ("within_two_voxels", o.acc.within_two.to_string()),
            (
                "within_fraction_uniform",
                format!(
                    "{:.6}",
                    o.acc.uniform_within as f64 / o.acc.uniform.max(1) as f64
                ),
            ),
            ("uniform_samples", o.acc.uniform.to_string()),
            ("narrow_within_one_voxel", o.narrow_within.to_string()),
            ("mean_error_cells", format!("{mean_err:.6}")),
            ("median_error_cells", format!("{median:.6}")),
            ("p90_error_cells", format!("{p90:.6}")),
            (
                "max_understatement_cells",
                format!("{:.6}", o.acc.max_under),
            ),
            ("zero_bound_samples", o.acc.zero_bound.to_string()),
            ("probe_overstatements", o.acc.probe_over.to_string()),
            ("flips_outside_leaf_subset", o.outside_leaf.to_string()),
            (
                "worst_edit_flip_factor",
                format!("{:.4}", o.worst_edit_factor),
            ),
            ("worst_edit_index", o.worst_edit_index.to_string()),
            ("repair_set_cells", o.repair_set_cells.to_string()),
            ("repair_nav_cells", o.repair_nav_cells.to_string()),
            (
                "repair_set_fraction",
                format!("{:.6}", o.repair_set_cells as f64 / o.world_cells as f64),
            ),
            ("world_total_cells", o.world_cells.to_string()),
            ("free_leaves", o.free_leaves.to_string()),
            ("static_regions", o.regions.to_string()),
            ("all_leaves", o.all_leaves.to_string()),
            ("free_cells", o.free_cells.to_string()),
            (
                "bound_ns_per_query",
                format!("{:.1}", o.acc.ns as f64 / samples),
            ),
            ("clock_mhz", format!("{clock:.0}")),
        ];
        // The whole λ sweep, so the row's own choice of "worst creature" can be
        // checked rather than trusted.
        for i in 0..LAMBDAS.len() {
            row.push((
                FLIP_FACTOR_COLS[i],
                format!(
                    "{:.4}",
                    o.loc.flips[i] as f64 / o.changed_samples.max(1) as f64
                ),
            ));
            row.push((OUTSIDE_COLS[i], o.loc.outside_all[i].to_string()));
        }
        rows.push(row);
    }

    common::experiment::run(prereg, |run| {
        for row in &rows {
            run.record(row);
        }
    });
}
