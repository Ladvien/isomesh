//! **P-96 — how far apart smooth union's 40,317 answers actually are.**
//!
//! Ticket: R-096. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p96
//! ```
//!
//! Writes `docs/experiments/p-96.csv`.
//!
//! # The question, and why counting could not answer it
//!
//! `M-38` folded eight `SmoothAdd { k: 0.15 }` brushes in all `8! = 40,320`
//! orderings, hashed the resulting field at 64 probe points, and counted
//! **40,317 distinct results**. That is a count of distinct *bit patterns*. It
//! says nothing about how far apart the resulting *surfaces* are, and the
//! protocol decision — must a networked editor impose a total order on smooth
//! edits? — turns entirely on the distance, not the count.
//!
//! # Reachability, computed before the run (the SHARE line)
//!
//! The registration says "SHARE: this measures a distance, not a ratio", so
//! there is no share of a total to name. What still has to be checked is that
//! the 0.1-cell bar is reachable from both sides.
//!
//! `smooth_min(a, b, k)` satisfies `min(a,b) − k/4 ≤ smin ≤ min(a,b)` — the
//! penalty term is `k·h·(1−h)`, maximised at `h = 1/2`. A stack of `n` operands
//! folds `n − 1` times, so **every ordering's field lies within `(n−1)k/4` of
//! the ordering-invariant hard min**, and two orderings differ by at most
//! `(n−1)k/2`. Here `n = 9` — the base box plus eight brushes — so the fold is
//! within `2k` of hard `min` and two orderings are within `4k` of each other. On
//! a near-unit-gradient SDF a field perturbation `δ` moves the isosurface by
//! about `δ`, so the worst-case surface spread is `4k / cell_size` cells:
//!
//! | k | worst case, `m38` cells (0.078125) | vs the 0.1-cell bar |
//! |---:|---:|---:|
//! | 0.3 | 15.36 | 154x |
//! | 0.15 | 7.68 | 77x |
//! | 0.075 | 3.84 | 38x |
//! | 0.0375 | 1.92 | 19x |
//! | 0.01 | 0.512 | 5.1x |
//! | 0.001 | 0.0512 | **0.51x — bar guaranteed** |
//!
//! So C1 is **not** decided by arithmetic at `M-38`'s own `k`: the bound sits 77x
//! above the bar and the clause could fail. It *is* decided for `k ≲ 0.002`,
//! which is why the sweep runs down to `1e-6` — those arms confirm the bound
//! rather than test the clause.
//!
//! # The pairwise reduction, stated rather than assumed
//!
//! C1 names `max_{i,j} d_H(M_i, M_j)` over 40,320 meshes — 812,838,080 pairs.
//! Symmetric Hausdorff is a metric on compact sets, so for any reference `R` the
//! triangle inequality gives
//!
//! ```text
//! max_{i,j} d_H(M_i, M_j) ≤ max_i d_H(M_i, R) + max_j d_H(M_j, R)
//! ```
//!
//! and since `i ≠ j` may be assumed, the **two largest** distances-to-reference
//! bound the diameter: `top1 + top2`. That is `max_hausdorff_bound_cells`, it is
//! computed from **all** 40,320 orderings against the identity ordering, and it
//! is the column C1 is scored on — a bound under the bar proves the registered
//! quantity is under the bar. `max_hausdorff_cells` carries the **exact**
//! all-pairs maximum over a sampled subset of at most 128 orderings, built from
//! the 64 largest distances-to-reference plus 64 evenly strided, so it is a
//! lower bound biased towards the extremes. Both are reported and neither is
//! substituted for the other.
//!
//! # Two fixtures, and the second one is the finding
//!
//! `m38` is `M-38`'s fixture verbatim: `BoxExact::canonical()` — the `[-1,1]³`
//! cube — plus the eight overlapping brushes of `brush::tests::eight`, and the
//! same 4³ probe lattice. It reproduces `M-38` exactly, which is the registered
//! vacuity control.
//!
//! It also cannot answer the question, and the harness's own controls say so.
//! Every brush lies **strictly inside** the base cube, so the base is the
//! smallest of the nine operands everywhere and the runner-up is far away:
//! `min_margin_world` is the smallest gap between the two smallest operands
//! anywhere on the grid, and it exceeds `k = 0.15`. `smooth_min` clamps `h` to
//! `0` or `1` whenever `|b − a| ≥ k`, so at that `k` **not one fold step on the
//! whole grid enters the smooth branch** — `blended_samples` is zero. The 40,317
//! distinct results are the cancellation residue of evaluating
//! `(b + (a − b)·1) − k·1·0` instead of returning `a`: algebraically the hard
//! min, numerically a few ULP. `m38`'s C1 answer is therefore a zero that could
//! not have been anything else, and it is scored **vacuous** rather than held.
//!
//! `m38_exposed` keeps all eight brushes, all eight `SmoothAdd` ops, the probe
//! lattice and the permutation machinery, and shrinks the base cube's
//! half-extent from `1.0` to `0.35` so the brush union protrudes and brush-brush
//! seams reach the surface. The base is not one of the permuted objects — the
//! fold seeds with it and permutes the eight brushes — so this is the smallest
//! edit that makes the blend reachable while leaving the permuted fixture
//! untouched. It is a deviation from the registration and is labelled as one in
//! the `fixture` column.
//!
//! # Controls
//!
//! - **`orderings == 40320`**, and the 40,320 permutations are asserted pairwise
//!   distinct, inheriting `there_really_are_40320_orderings`.
//! - **`distinct_results == 40317` and `coincident_orderings == 3` on the `m38`
//!   fixture at `k = 0.15`.** The registered vacuity control: a harness that
//!   does not reproduce `M-38`'s count *and* its three coincidences is not
//!   measuring `M-38`'s fixture.
//! - **`stack_agreement == 64`.** The precomputed nine-operand fold is compared
//!   bit-for-bit against `BrushStack::sample` at every probe for the identity
//!   ordering. Sampling the eight shapes once and re-folding is the reduction
//!   that makes 40,320 orderings affordable, and it is only legitimate if it is
//!   the same arithmetic; this asserts that it is.
//! - **`grid_index_roundtrip == grid_samples`.** The cached-grid `Sdf` recovers
//!   its index from the position Marching Cubes asks for, verified exhaustively
//!   per arm, because one mis-rounded index would silently mesh a wrong field.
//! - **`boundary_min_world > 0`.** The minimum, over every grid-boundary sample
//!   and every one of the 40,320 orderings, of the folded field. Positive means
//!   the surface is closed inside the grid, so no Hausdorff distance is
//!   measuring a clipped mesh edge.
//! - **`ref_vertices > 0` and `ref_triangles > 0`.** A non-empty population on
//!   both sides of every Hausdorff query.
//! - **Sweep-level `M-44` controls**, asserted after every arm has run: at least
//!   one arm reaches the smooth branch (`blended_samples > 0`), at least one arm
//!   reports a Hausdorff distance above 0.01 cells — the instrument can return
//!   non-zero — and at least one arm has both a non-empty outside-shell
//!   population and spread inside the shell, so C3 has somewhere to fail.
//! - **The reduction is asserted sound**: the exact sampled-subset maximum must
//!   not exceed `top1 + top2`.
//!
//! # How a clause is scored `vacuous`
//!
//! - **C1** on `active_samples == 0`. C1 asks how far apart two *orderings* are,
//!   so the zero it could have avoided is a field that depends on the ordering
//!   somewhere above `spread_floor_world`. A blend every ordering evaluates
//!   identically leaves it nothing to measure — which is exactly `m38` at
//!   `k = 0.3`: 777 `blended_samples`, zero `active_samples`.
//! - **C2** on fewer than three arms above `SLOPE_FLOOR_CELLS`, since a slope
//!   through two points is not a measurement of scaling.
//! - **C3** on `outside_shell_population == 0` — the `10k` shell has swallowed
//!   the grid, so `points_outside_shell` cannot be non-zero — or on
//!   `inside_shell_with_spread == 0`, where there is no spread to confine.
//!
//! # C3's columns, and which reading of the registration they take
//!
//! `max_deviation_from_seam_cells` is read as *the largest distance from a seam
//! at which ordering spread was observed*, measured in the memo's own
//! coordinate — the margin `f₍₂₎ − f₍₁₎` between the two smallest operands,
//! which `docs/research/2026-08-23-unmined-mathematics-for-meshing.md` §3 names
//! as "the coordinate that measures it". `shell_10k_cells` is `10k` in the same
//! unit and `points_outside_shell` counts grid samples beyond it that still move
//! with the ordering. That triple instruments C3.
//!
//! The memo's *other* quantity — `max |smin_k − min|`, measured there at 0.135
//! against a bound of 0.208 (`k·ln n` at `k = 0.1`, `n = 8`) — is emitted as
//! `max_deviation_world` beside `deviation_bound_kln9_world` (`k·ln 9`, since
//! this fixture folds nine operands) and `deviation_bound_poly_world`
//! (`(n−1)k/4 = 2k`, the exact bound for *this* smooth-min, and tighter than
//! `k·ln 9` for every `n ≤ 9`). Both readings are in the file so neither has to
//! be reconstructed.
//!
//! A sample "moves with the ordering" when its spread over all 40,320 orderings
//! exceeds `spread_floor_world = 1e-12`. That is three orders above the f64
//! cancellation residue this fixture produces (1.4e-15) and eleven orders below
//! the geometric scale, so the threshold separates arithmetic from geometry
//! rather than tuning the answer.
//!
//! # References
//!
//! Quilez's polynomial smooth minimum, as implemented in `brush::smooth_min`.
//! The log-semiring / Maslov-dequantisation reading of `k` as a length, and the
//! `|smin_k − min| ≤ k·ln n` bound, are
//! `docs/research/2026-08-23-unmined-mathematics-for-meshing.md` §3.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::print_literal,
    clippy::too_many_lines
)]

mod common;

use std::num::NonZeroUsize;

use isomesh::brush::{Brush, BrushOp, BrushStack, Capsule, apply};
use isomesh::fields::{BoxExact, Sphere};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Brushes in the fixture. `M-38`'s eight.
const BRUSHES: usize = 8;

/// `8!`. The registered vacuity control is a count against this.
const ORDERINGS: usize = 40_320;

/// Probes in `M-38`'s signature lattice: `4³`.
const PROBES: usize = 64;

/// Operands in the fold: the base field plus the eight brushes.
const OPERANDS: usize = BRUSHES + 1;

/// `M-38`'s own blend radius. The arm the registered assertions apply to.
const M38_K: f64 = 0.15;

/// The blend-radius sweep, ascending.
///
/// `0.15` is `M-38`'s. `0.3`, `0.075` and `0.0375` are its octave neighbours;
/// `0.01`, `1e-3` and `1e-6` walk down towards the regime where the `4k` bound
/// alone guarantees the 0.1-cell bar. `0.0` is the tropical limit and takes
/// `smooth_min`'s `k <= 0` branch, which is an ordinary `min` — the only arm
/// where `M-36`'s single distinct result can be recovered.
const KS: [f64; 8] = [0.0, 1e-6, 1e-3, 0.01, 0.0375, 0.075, M38_K, 0.3];

/// Cells per axis in the extraction grid.
const GRID_CELLS: u32 = 32;

/// C1's bar, in cells.
const C1_BAR: f64 = 0.1;

/// Below this, a per-sample spread is f64 cancellation rather than geometry.
const SPREAD_FLOOR: f64 = 1e-12;

/// A Hausdorff distance below this, in cells, is excluded from C2's log-log fit.
const SLOPE_FLOOR_CELLS: f64 = 1e-6;

/// Orderings taken into the exact all-pairs subset, from each of two sources:
/// the largest distances-to-reference, and an even stride.
const SUBSET_TOP: usize = 64;

/// Threads. Capped so the harness's numbers do not depend on the host's core
/// count — only its runtime does.
const MAX_THREADS: usize = 12;

/// A brush shape, as one type so a stack can hold a mixture.
///
/// Identical to `brush::tests::Shape`, which is private to that module.
#[derive(Clone, Copy, Debug)]
enum Shape {
    Sphere(Sphere<f64>),
    Cube(BoxExact<f64>),
    Capsule(Capsule<f64>),
}

impl Sdf for Shape {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Sphere(s) => s.sample(p),
            Self::Cube(b) => b.sample(p),
            Self::Capsule(c) => c.sample(p),
        }
    }
}

/// `M-38`'s eight brushes, all three shapes represented, deliberately
/// overlapping so the order could matter.
///
/// Copied from `brush::tests::eight` — the centres, the radii, the half-extents
/// and the `i % 3` shape rotation *are* the fixture, and changing any of them
/// would stop this from being `M-38`'s measurement.
fn eight(k: f64) -> Vec<Brush<Shape>> {
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
            Brush {
                shape,
                op: BrushOp::SmoothAdd { k },
            }
        })
        .collect()
}

/// `M-38`'s signature lattice: `4³` points spanning the region the brushes
/// occupy.
fn probes() -> Vec<[f64; 3]> {
    let mut out = Vec::with_capacity(PROBES);
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..4 {
                out.push([
                    -0.6 + 0.4 * f64::from(x),
                    -0.6 + 0.4 * f64::from(y),
                    -0.6 + 0.4 * f64::from(z),
                ]);
            }
        }
    }
    out
}

/// Every permutation of `0..BRUSHES`, by Heap's algorithm.
///
/// Generated rather than sampled, exactly as `M-38` does it: "we tried a
/// thousand random orderings" is a different and weaker claim.
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

/// One of the two fixtures.
#[derive(Clone, Copy)]
struct Fixture {
    /// The `fixture` column.
    name: &'static str,
    /// Half-extent of the base cube. `1.0` is `M-38`'s `BoxExact::canonical()`.
    base_half: f64,
    /// Lowest grid corner, on every axis.
    origin: f64,
    /// Grid extent, on every axis.
    extent: f64,
}

/// `m38` verbatim, then the smallest edit that exposes a brush-brush seam.
const FIXTURES: [Fixture; 2] = [
    Fixture {
        name: "m38",
        base_half: 1.0,
        origin: -1.25,
        extent: 2.5,
    },
    Fixture {
        name: "m38_exposed",
        base_half: 0.35,
        origin: -1.0,
        extent: 2.0,
    },
];

/// Cell size of a fixture's extraction grid.
fn cell_of(fx: Fixture) -> f64 {
    fx.extent / f64::from(GRID_CELLS)
}

/// The nine operand values at every grid sample and every probe, sampled once.
///
/// The eight shapes and the base do not depend on the ordering — only the fold
/// does — so this is the whole reason 40,320 orderings fit inside a bench.
struct Operands {
    size: [u32; 3],
    origin: [f64; 3],
    cell: f64,
    /// `OPERANDS` values per grid sample, base first.
    grid: Vec<[f64; OPERANDS]>,
    /// `OPERANDS` values per probe, base first.
    probe: Vec<[f64; OPERANDS]>,
    /// Ordering-invariant hard `min` over the nine operands, per grid sample.
    hard: Vec<f64>,
    /// `f₍₂₎ − f₍₁₎`, the memo's seam coordinate, per grid sample.
    margin: Vec<f64>,
}

/// Sample the nine operands over the grid and the probe lattice.
fn operands(fx: Fixture, brushes: &[Brush<Shape>]) -> Operands {
    let n = GRID_CELLS + 1;
    let cell = cell_of(fx);
    let origin = [fx.origin; 3];
    let base = BoxExact::<f64> {
        center: [0.0; 3],
        half_extents: [fx.base_half; 3],
    };

    let count = (n as usize).pow(3);
    let mut grid = Vec::with_capacity(count);
    let mut hard = Vec::with_capacity(count);
    let mut margin = Vec::with_capacity(count);
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let p = [
                    origin[0] + cell * f64::from(x),
                    origin[1] + cell * f64::from(y),
                    origin[2] + cell * f64::from(z),
                ];
                let v = nine(&base, brushes, p);
                let (lo, second) = two_smallest(&v);
                grid.push(v);
                hard.push(lo);
                margin.push(second - lo);
            }
        }
    }

    let probe = probes().iter().map(|p| nine(&base, brushes, *p)).collect();

    Operands {
        size: [n; 3],
        origin,
        cell,
        grid,
        probe,
        hard,
        margin,
    }
}

/// The nine operand values at one point, base first.
fn nine(base: &BoxExact<f64>, brushes: &[Brush<Shape>], p: [f64; 3]) -> [f64; OPERANDS] {
    let mut v = [0.0f64; OPERANDS];
    v[0] = base.sample(p);
    for (slot, brush) in v[1..].iter_mut().zip(brushes) {
        *slot = brush.shape.sample(p);
    }
    v
}

/// The two smallest of the nine operands, in order.
fn two_smallest(v: &[f64; OPERANDS]) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut second = f64::INFINITY;
    for &x in v {
        if x < lo {
            second = lo;
            lo = x;
        } else if x < second {
            second = x;
        }
    }
    (lo, second)
}

/// Fold the nine operands in one ordering, through the crate's own `apply`.
#[inline]
fn fold(v: &[f64; OPERANDS], order: &[usize; BRUSHES], op: BrushOp) -> f64 {
    let mut x = v[0];
    for &i in order {
        x = apply(op, x, v[1 + i]);
    }
    x
}

/// A grid of already-folded field values, presented as an `Sdf`.
///
/// Marching Cubes samples only at grid corners — `sdf::sample_grid`, and
/// `crossing_refinement` is `0` in the default configuration — so recovering the
/// index by rounding is exact. `grid_index_roundtrip` verifies that
/// exhaustively rather than trusting it.
struct Cached<'a> {
    values: &'a [f64],
    size: [u32; 3],
    origin: [f64; 3],
    inv_cell: f64,
}

impl Sdf for Cached<'_> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        let ix = ((p[0] - self.origin[0]) * self.inv_cell).round() as usize;
        let iy = ((p[1] - self.origin[1]) * self.inv_cell).round() as usize;
        let iz = ((p[2] - self.origin[2]) * self.inv_cell).round() as usize;
        let nx = self.size[0] as usize;
        let ny = self.size[1] as usize;
        self.values[ix + nx * (iy + ny * iz)]
    }
}

/// A bucketed point set, for **exact** nearest-neighbour queries.
///
/// Buckets are the extraction grid's own cells, so a Marching Cubes vertex and
/// its counterpart from another ordering share a bucket unless the field moved
/// by a whole cell.
struct PointIndex {
    dims: [usize; 3],
    origin: [f64; 3],
    cell: f64,
    starts: Vec<u32>,
    cursor: Vec<u32>,
    items: Vec<u32>,
    pts: Vec<[f64; 3]>,
}

impl PointIndex {
    fn new(dims: [usize; 3], origin: [f64; 3], cell: f64) -> Self {
        let buckets = dims[0] * dims[1] * dims[2];
        Self {
            dims,
            origin,
            cell,
            starts: vec![0; buckets + 1],
            cursor: vec![0; buckets + 1],
            items: Vec::new(),
            pts: Vec::new(),
        }
    }

    #[inline]
    fn bucket(&self, p: [f64; 3]) -> [usize; 3] {
        let mut b = [0usize; 3];
        for (axis, slot) in b.iter_mut().enumerate() {
            let t = (p[axis] - self.origin[axis]) / self.cell;
            let clamped = if t < 0.0 { 0.0 } else { t };
            *slot = (clamped as usize).min(self.dims[axis] - 1);
        }
        b
    }

    /// Rebuild over `pts`, reusing every allocation.
    fn build(&mut self, pts: &[[f64; 3]]) {
        self.pts.clear();
        self.pts.extend_from_slice(pts);
        self.starts.fill(0);
        let nx = self.dims[0];
        let ny = self.dims[1];
        for p in pts {
            let b = self.bucket(*p);
            self.starts[b[0] + nx * (b[1] + ny * b[2])] += 1;
        }
        let mut acc = 0u32;
        for s in &mut self.starts {
            let c = *s;
            *s = acc;
            acc += c;
        }
        self.cursor.copy_from_slice(&self.starts);
        self.items.clear();
        self.items.resize(pts.len(), 0);
        for (i, p) in pts.iter().enumerate() {
            let b = self.bucket(*p);
            let cell = b[0] + nx * (b[1] + ny * b[2]);
            self.items[self.cursor[cell] as usize] = i as u32;
            self.cursor[cell] += 1;
        }
    }

    /// Exact squared distance from `p` to the nearest indexed point.
    ///
    /// After scanning every bucket within Chebyshev radius `r`, an unscanned
    /// point lies in a bucket at radius `≥ r + 1`, and `p` lies inside the centre
    /// bucket, so that point is at least `r · cell` away. `best ≤ (r · cell)²` is
    /// therefore a proof of exactness rather than a heuristic — and the loop
    /// panics instead of returning an unproven answer.
    fn nearest_sq(&self, p: [f64; 3]) -> f64 {
        let b = self.bucket(p);
        let nx = self.dims[0];
        let ny = self.dims[1];
        let limit = self.dims[0].max(self.dims[1]).max(self.dims[2]);
        let mut best = f64::INFINITY;
        for r in 1..=limit {
            let x0 = b[0].saturating_sub(r);
            let x1 = (b[0] + r).min(self.dims[0] - 1);
            let y0 = b[1].saturating_sub(r);
            let y1 = (b[1] + r).min(self.dims[1] - 1);
            let z0 = b[2].saturating_sub(r);
            let z1 = (b[2] + r).min(self.dims[2] - 1);
            for z in z0..=z1 {
                for y in y0..=y1 {
                    let row = nx * (y + ny * z);
                    for x in x0..=x1 {
                        let cell = row + x;
                        let lo = self.starts[cell] as usize;
                        let hi = self.starts[cell + 1] as usize;
                        for &item in &self.items[lo..hi] {
                            let q = self.pts[item as usize];
                            let dx = p[0] - q[0];
                            let dy = p[1] - q[1];
                            let dz = p[2] - q[2];
                            let d = dx * dx + dy * dy + dz * dz;
                            if d < best {
                                best = d;
                            }
                        }
                    }
                }
            }
            let proven = r as f64 * self.cell;
            if best <= proven * proven {
                return best;
            }
        }
        panic!("nearest-neighbour search exhausted the grid: the indexed point set is empty");
    }
}

/// Symmetric Hausdorff distance between two indexed vertex sets, in world units.
fn symmetric_hausdorff(a: &PointIndex, b: &PointIndex) -> f64 {
    let mut worst = 0.0f64;
    for p in &a.pts {
        let d = b.nearest_sq(*p);
        if d > worst {
            worst = d;
        }
    }
    for q in &b.pts {
        let d = a.nearest_sq(*q);
        if d > worst {
            worst = d;
        }
    }
    worst.sqrt()
}

/// Buffers one worker reuses across every ordering it handles.
struct Scratch {
    mc: MarchingCubes<f64>,
    field: Vec<f64>,
    mesh: MeshBuffer<f64>,
    index: PointIndex,
}

impl Scratch {
    fn new(ops: &Operands, dims: [usize; 3]) -> Self {
        Self {
            mc: MarchingCubes::<f64>::new(),
            field: Vec::with_capacity(ops.hard.len()),
            mesh: MeshBuffer::<f64>::new(),
            index: PointIndex::new(dims, ops.origin, ops.cell),
        }
    }

    /// Fold one ordering onto the grid, mesh it, and index its vertices.
    fn mesh_ordering(
        &mut self,
        ops: &Operands,
        shape: &RuntimeShape3,
        order: &[usize; BRUSHES],
        op: BrushOp,
    ) {
        self.field.clear();
        self.field
            .extend(ops.grid.iter().map(|v| fold(v, order, op)));
        let cached = Cached {
            values: &self.field,
            size: ops.size,
            origin: ops.origin,
            inv_cell: 1.0 / ops.cell,
        };
        self.mesh.reset();
        self.mc
            .extract(&cached, shape, ops.origin, ops.cell, &mut self.mesh)
            .expect("marching cubes over a finite cached grid");
        self.index.build(&self.mesh.positions);
    }
}

/// What one ordering contributed.
#[derive(Clone)]
struct Slot {
    /// Symmetric Hausdorff to the reference mesh, in world units.
    hausdorff: f64,
    vertices: u32,
    triangles: u32,
    /// `M-38`'s signature: the field at the 64 probes, as raw bits.
    sig: [u64; PROBES],
}

/// Per-worker accumulators, merged after the parallel loop.
struct Acc {
    lo: Vec<f64>,
    hi: Vec<f64>,
}

/// Everything one `(fixture, k)` arm measured, before the sweep-level verdicts.
struct Arm {
    k: f64,
    distinct: usize,
    coincident: usize,
    groups: usize,
    largest_group: usize,
    stack_agreement: usize,
    roundtrip: usize,
    boundary_min: f64,
    min_margin: f64,
    ref_vertices: usize,
    ref_triangles: usize,
    max_vertices: usize,
    max_to_ref_cells: f64,
    bound_cells: f64,
    exact_pair_max_cells: f64,
    exact_pairs: usize,
    subset: usize,
    mean_cells: f64,
    p99_cells: f64,
    blended: usize,
    active: usize,
    sign_flips: usize,
    topology_changes: usize,
    max_spread_world: f64,
    max_dev_world: f64,
    max_dev_from_seam_cells: f64,
    shell_cells: f64,
    outside_pop: usize,
    outside_with_spread: usize,
    inside_with_spread: usize,
    max_spread_outside_world: f64,
}

/// Run one `(fixture, k)` arm over all 40,320 orderings.
fn run_arm(fx: Fixture, k: f64, perms: &[[usize; BRUSHES]], threads: usize) -> Arm {
    let brushes = eight(k);
    let ops = operands(fx, &brushes);
    let op = BrushOp::SmoothAdd { k };
    let samples = ops.hard.len();
    let shape = RuntimeShape3::new(ops.size).expect("a 33-cubed grid fits u32");
    let dims = [GRID_CELLS as usize; 3];

    // ── the fold reduction is the crate's own arithmetic ────────────────────
    let stack = BrushStack {
        base: BoxExact::<f64> {
            center: [0.0; 3],
            half_extents: [fx.base_half; 3],
        },
        brushes: &brushes,
    };
    let stack_agreement = probes()
        .iter()
        .zip(&ops.probe)
        .filter(|(p, v)| fold(v, &perms[0], op).to_bits() == stack.sample(**p).to_bits())
        .count();

    // ── reference mesh: the identity ordering ───────────────────────────────
    let mut scratch = Scratch::new(&ops, dims);
    scratch.mesh_ordering(&ops, &shape, &perms[0], op);
    let ref_index = {
        let mut ix = PointIndex::new(dims, ops.origin, ops.cell);
        ix.build(&scratch.mesh.positions);
        ix
    };
    let ref_vertices = scratch.mesh.positions.len();
    let ref_triangles = scratch.mesh.indices.len() / 3;

    // ── the cached grid returns what it was asked for, at every sample ──────
    let roundtrip = {
        let cached = Cached {
            values: &scratch.field,
            size: ops.size,
            origin: ops.origin,
            inv_cell: 1.0 / ops.cell,
        };
        let n = ops.size[0] as usize;
        let mut ok = 0usize;
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let p = [
                        ops.origin[0] + ops.cell * x as f64,
                        ops.origin[1] + ops.cell * y as f64,
                        ops.origin[2] + ops.cell * z as f64,
                    ];
                    if cached.sample(p) == scratch.field[x + n * (y + n * z)] {
                        ok += 1;
                    }
                }
            }
        }
        ok
    };

    // ── all 40,320 orderings ────────────────────────────────────────────────
    let mut slots = vec![
        Slot {
            hausdorff: 0.0,
            vertices: 0,
            triangles: 0,
            sig: [0; PROBES],
        };
        ORDERINGS
    ];
    let chunk = ORDERINGS.div_ceil(threads);
    let accs: Vec<Acc> = std::thread::scope(|s| {
        let mut handles = Vec::new();
        for (slot_chunk, perm_chunk) in slots.chunks_mut(chunk).zip(perms.chunks(chunk)) {
            let ops = &ops;
            let shape = &shape;
            let ref_index = &ref_index;
            handles.push(s.spawn(move || {
                let mut acc = Acc {
                    lo: vec![f64::INFINITY; samples],
                    hi: vec![f64::NEG_INFINITY; samples],
                };
                let mut scratch = Scratch::new(ops, dims);
                for (slot, order) in slot_chunk.iter_mut().zip(perm_chunk) {
                    for (bits, v) in slot.sig.iter_mut().zip(&ops.probe) {
                        *bits = fold(v, order, op).to_bits();
                    }
                    scratch.mesh_ordering(ops, shape, order, op);
                    for ((x, lo), hi) in scratch
                        .field
                        .iter()
                        .zip(acc.lo.iter_mut())
                        .zip(acc.hi.iter_mut())
                    {
                        if *x < *lo {
                            *lo = *x;
                        }
                        if *x > *hi {
                            *hi = *x;
                        }
                    }
                    slot.vertices = scratch.mesh.positions.len() as u32;
                    slot.triangles = (scratch.mesh.indices.len() / 3) as u32;
                    slot.hausdorff = symmetric_hausdorff(&scratch.index, ref_index);
                }
                acc
            }));
        }
        handles
            .into_iter()
            .map(|h| h.join().expect("ordering worker finished"))
            .collect()
    });

    let mut lo = vec![f64::INFINITY; samples];
    let mut hi = vec![f64::NEG_INFINITY; samples];
    for acc in &accs {
        for ((l, h), (al, ah)) in lo
            .iter_mut()
            .zip(hi.iter_mut())
            .zip(acc.lo.iter().zip(acc.hi.iter()))
        {
            if *al < *l {
                *l = *al;
            }
            if *ah > *h {
                *h = *ah;
            }
        }
    }

    // ── M-38's own measurement, reproduced ──────────────────────────────────
    let mut sigs: Vec<[u64; PROBES]> = slots.iter().map(|s| s.sig).collect();
    sigs.sort_unstable();
    let mut distinct = 0usize;
    let mut groups = 0usize;
    let mut largest_group = 1usize;
    let mut i = 0usize;
    while i < sigs.len() {
        let mut j = i + 1;
        while j < sigs.len() && sigs[j] == sigs[i] {
            j += 1;
        }
        let count = j - i;
        distinct += 1;
        if count > 1 {
            groups += 1;
        }
        if count > largest_group {
            largest_group = count;
        }
        i = j;
    }
    let coincident = ORDERINGS - distinct;

    // ── the distance-to-reference distribution, and the diameter bound ──────
    let dh: Vec<f64> = slots.iter().map(|s| s.hausdorff / ops.cell).collect();
    let mean_cells = dh.iter().sum::<f64>() / ORDERINGS as f64;
    let mut ranked = dh.clone();
    ranked.sort_by(f64::total_cmp);
    let p99_cells = ranked[(0.99 * ORDERINGS as f64) as usize];
    let max_to_ref_cells = ranked[ORDERINGS - 1];
    let bound_cells = ranked[ORDERINGS - 1] + ranked[ORDERINGS - 2];

    // ── the exact all-pairs maximum over a subset biased to the extremes ────
    let mut order_by_distance: Vec<usize> = (0..ORDERINGS).collect();
    order_by_distance.sort_by(|&a, &b| dh[b].total_cmp(&dh[a]));
    let mut subset: Vec<usize> = order_by_distance[..SUBSET_TOP].to_vec();
    let stride = ORDERINGS / SUBSET_TOP;
    for m in 0..SUBSET_TOP {
        subset.push(m * stride);
    }
    subset.sort_unstable();
    subset.dedup();
    let subset_indices: Vec<PointIndex> = subset
        .iter()
        .map(|&o| {
            scratch.mesh_ordering(&ops, &shape, &perms[o], op);
            let mut ix = PointIndex::new(dims, ops.origin, ops.cell);
            ix.build(&scratch.mesh.positions);
            ix
        })
        .collect();
    let exact_pair_max_cells = std::thread::scope(|s| {
        let mut handles = Vec::new();
        for t in 0..threads {
            let members = &subset_indices;
            handles.push(s.spawn(move || {
                let mut best = 0.0f64;
                let mut a = t;
                while a < members.len() {
                    for b in (a + 1)..members.len() {
                        let d = symmetric_hausdorff(&members[a], &members[b]);
                        if d > best {
                            best = d;
                        }
                    }
                    a += threads;
                }
                best
            }));
        }
        handles
            .into_iter()
            .map(|h| h.join().expect("pair worker finished"))
            .fold(0.0f64, f64::max)
    }) / ops.cell;
    let exact_pairs = subset_indices.len() * (subset_indices.len() - 1) / 2;

    // ── C3: where the spread lives, in the memo's margin coordinate ─────────
    let shell = 10.0 * k;
    let mut blended = 0usize;
    let mut active = 0usize;
    let mut sign_flips = 0usize;
    let mut max_spread_world = 0.0f64;
    let mut max_dev_world = 0.0f64;
    let mut max_dev_from_seam_cells = 0.0f64;
    let mut outside_pop = 0usize;
    let mut outside_with_spread = 0usize;
    let mut inside_with_spread = 0usize;
    let mut max_spread_outside_world = 0.0f64;
    for idx in 0..samples {
        let spread = hi[idx] - lo[idx];
        let dev = (hi[idx] - ops.hard[idx])
            .abs()
            .max((lo[idx] - ops.hard[idx]).abs());
        if dev > SPREAD_FLOOR {
            blended += 1;
        }
        if spread > max_spread_world {
            max_spread_world = spread;
        }
        if dev > max_dev_world {
            max_dev_world = dev;
        }
        if (lo[idx] < 0.0) != (hi[idx] < 0.0) {
            sign_flips += 1;
        }
        let moves = spread > SPREAD_FLOOR;
        if moves {
            active += 1;
            let from_seam = ops.margin[idx] / ops.cell;
            if from_seam > max_dev_from_seam_cells {
                max_dev_from_seam_cells = from_seam;
            }
        }
        if ops.margin[idx] > shell {
            outside_pop += 1;
            if moves {
                outside_with_spread += 1;
            }
            if spread > max_spread_outside_world {
                max_spread_outside_world = spread;
            }
        } else if moves {
            inside_with_spread += 1;
        }
    }

    // ── the mesh is closed inside the grid, for every ordering ──────────────
    let n = ops.size[0] as usize;
    let mut boundary_min = f64::INFINITY;
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                if x == 0 || y == 0 || z == 0 || x == n - 1 || y == n - 1 || z == n - 1 {
                    let v = lo[x + n * (y + n * z)];
                    if v < boundary_min {
                        boundary_min = v;
                    }
                }
            }
        }
    }

    let min_margin = ops.margin.iter().copied().fold(f64::INFINITY, f64::min);
    let max_vertices = slots.iter().map(|s| s.vertices as usize).max().unwrap_or(0);
    let topology_changes = slots
        .iter()
        .filter(|s| s.triangles as usize != ref_triangles)
        .count();

    println!(
        "  {:<12} k={:<9} distinct={distinct:<6} coincident={coincident} \
         blended={blended:<6} maxHd_to_ref={max_to_ref_cells:.6}c bound={bound_cells:.6}c \
         exact_pairs_max={exact_pair_max_cells:.6}c seam={max_dev_from_seam_cells:.3}c/\
         {:.3}c outside={outside_with_spread}",
        fx.name,
        k,
        10.0 * k / ops.cell
    );

    Arm {
        k,
        distinct,
        coincident,
        groups,
        largest_group,
        stack_agreement,
        roundtrip,
        boundary_min,
        min_margin,
        ref_vertices,
        ref_triangles,
        max_vertices,
        max_to_ref_cells,
        bound_cells,
        exact_pair_max_cells,
        exact_pairs,
        subset: subset_indices.len(),
        mean_cells,
        p99_cells,
        blended,
        active,
        sign_flips,
        topology_changes,
        max_spread_world,
        max_dev_world,
        max_dev_from_seam_cells,
        shell_cells: 10.0 * k / ops.cell,
        outside_pop,
        outside_with_spread,
        inside_with_spread,
        max_spread_outside_world,
    }
}

/// Ordinary least squares slope of `y` on `x`.
fn least_squares(pts: &[(f64, f64)]) -> f64 {
    let n = pts.len() as f64;
    let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
    let num: f64 = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let den: f64 = pts.iter().map(|p| (p.0 - mx).powi(2)).sum();
    num / den
}

/// One CSV row.
type Row = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|a| a == "--bench") {
        return;
    }

    common::experiment::run(isomesh::experiment!("P-96"), |run| {
        let perms = permutations();
        assert_eq!(perms.len(), ORDERINGS, "8! orderings");
        let mut sorted = perms.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(
            sorted.len(),
            ORDERINGS,
            "the 40,320 permutations must be distinct"
        );

        let threads = std::thread::available_parallelism()
            .map_or(1, NonZeroUsize::get)
            .min(MAX_THREADS);
        println!("  threads: {threads}\n");

        let grid_samples = ((GRID_CELLS + 1) as usize).pow(3);
        let mut rows: Vec<Row> = Vec::new();
        let mut any_blend = false;
        let mut any_distance = false;
        let mut any_c3_reachable = false;

        for fx in FIXTURES {
            let arms: Vec<Arm> = KS
                .iter()
                .map(|&k| run_arm(fx, k, &perms, threads))
                .collect();

            // C2 is a property of the sweep, not of one arm, so its slope and
            // its verdict are computed here and repeated on every row of this
            // fixture. Repeated deliberately: `M-377`'s defect was one arm's
            // value on another arm's row, and this is the opposite — a
            // sweep-level scalar, named as one.
            let fit: Vec<(f64, f64)> = arms
                .iter()
                .filter(|a| a.k > 0.0 && a.max_to_ref_cells > SLOPE_FLOOR_CELLS)
                .map(|a| (a.k.ln(), a.max_to_ref_cells.ln()))
                .collect();
            let slope = if fit.len() >= 3 {
                least_squares(&fit)
            } else {
                0.0
            };
            // The same fit on the *field* spread rather than the mesh's. The
            // mesh number is what C1 is denominated in, but it is not a power
            // law -- a sample that flips sign adds a whole vertex, so the mesh
            // curve has a knee where the field curve does not. Both slopes are
            // reported because C2's claim is about the field's dependence on k
            // and C1's bar is about the mesh's.
            let field_fit: Vec<(f64, f64)> = arms
                .iter()
                .filter(|a| a.k > 0.0 && a.max_spread_world / cell_of(fx) > SLOPE_FLOOR_CELLS)
                .map(|a| (a.k.ln(), (a.max_spread_world / cell_of(fx)).ln()))
                .collect();
            let field_slope = if field_fit.len() >= 3 {
                least_squares(&field_fit)
            } else {
                0.0
            };
            let monotone = arms
                .windows(2)
                .all(|w| w[1].max_to_ref_cells >= w[0].max_to_ref_cells - 1e-15);
            let limit_recovers_m36 = arms.iter().any(|a| a.k == 0.0 && a.distinct == 1);
            let c2 = if fit.len() < 3 {
                "vacuous"
            } else if slope >= 0.5 && monotone && limit_recovers_m36 {
                "true"
            } else {
                "false"
            };

            for arm in &arms {
                any_blend |= arm.blended > 0;
                any_distance |= arm.max_to_ref_cells > 0.01;
                any_c3_reachable |= arm.outside_pop > 0 && arm.inside_with_spread > 0;

                assert_eq!(
                    arm.stack_agreement, PROBES,
                    "{}: the precomputed fold must equal BrushStack::sample bit-for-bit at every \
                     probe, or the reduction that makes 40,320 orderings affordable is not the \
                     crate's arithmetic",
                    fx.name
                );
                assert_eq!(
                    arm.roundtrip, grid_samples,
                    "{}: the cached grid must return the value Marching Cubes asked for at every \
                     sample",
                    fx.name
                );
                assert!(
                    arm.boundary_min > 0.0,
                    "{}: the folded field reaches {} on the grid boundary, so the surface is \
                     clipped and the Hausdorff distances measure a cut edge",
                    fx.name,
                    arm.boundary_min
                );
                assert!(
                    arm.ref_vertices > 0 && arm.ref_triangles > 0,
                    "{}: the reference mesh is empty, so no Hausdorff distance means anything",
                    fx.name
                );
                assert!(
                    arm.exact_pair_max_cells <= arm.bound_cells + 1e-12,
                    "{}: the exact sampled maximum {} exceeds the top1+top2 bound {}, so the \
                     triangle-inequality reduction is wrong",
                    fx.name,
                    arm.exact_pair_max_cells,
                    arm.bound_cells
                );
                if fx.name == "m38" && arm.k == M38_K {
                    assert_eq!(
                        arm.distinct, 40_317,
                        "M-38's own count must reproduce, or this is not M-38's fixture"
                    );
                    assert_eq!(
                        arm.coincident, 3,
                        "M-38's three coincident orderings must still coincide"
                    );
                }

                // Gated on `active_samples`, not on `blended_samples`. C1 asks
                // how far apart two *orderings* are, so the zero it could have
                // avoided is a field that depends on the ordering somewhere. A
                // blend that every ordering evaluates identically -- which is
                // what `m38` at k = 0.3 is, 777 blended samples and zero
                // active ones -- leaves C1 with nothing to measure.
                let c1 = if arm.active == 0 {
                    "vacuous"
                } else if arm.bound_cells < C1_BAR {
                    "true"
                } else {
                    "false"
                };
                let c3 = if arm.outside_pop == 0 || arm.inside_with_spread == 0 {
                    "vacuous"
                } else if arm.outside_with_spread == 0 {
                    "true"
                } else {
                    "false"
                };
                let cell = cell_of(fx);

                rows.push(vec![
                    ("fixture", fx.name.to_string()),
                    ("base_half_extent", format!("{:.4}", fx.base_half)),
                    ("blend_radius_k", format!("{:.9}", arm.k)),
                    ("orderings", ORDERINGS.to_string()),
                    ("distinct_results", arm.distinct.to_string()),
                    ("coincident_orderings", arm.coincident.to_string()),
                    ("coincident_groups", arm.groups.to_string()),
                    ("largest_coincident_group", arm.largest_group.to_string()),
                    (
                        "max_hausdorff_cells",
                        format!("{:.9}", arm.exact_pair_max_cells),
                    ),
                    (
                        "max_hausdorff_bound_cells",
                        format!("{:.9}", arm.bound_cells),
                    ),
                    (
                        "max_hausdorff_to_ref_cells",
                        format!("{:.9}", arm.max_to_ref_cells),
                    ),
                    (
                        // The fixed-point columns above round a few-ULP
                        // difference to a wall of zeros, and "exactly zero" and
                        // "1e-14 cells" are different findings.
                        "max_hausdorff_to_ref_world",
                        format!("{:.6e}", arm.max_to_ref_cells * cell),
                    ),
                    ("mean_hausdorff_cells", format!("{:.9}", arm.mean_cells)),
                    ("p99_hausdorff_cells", format!("{:.9}", arm.p99_cells)),
                    ("exact_pairs_sampled", arm.exact_pairs.to_string()),
                    ("subset_orderings", arm.subset.to_string()),
                    ("spread_vs_k_slope", format!("{slope:.6}")),
                    ("spread_vs_k_slope_field", format!("{field_slope:.6}")),
                    (
                        "max_deviation_from_seam_cells",
                        format!("{:.6}", arm.max_dev_from_seam_cells),
                    ),
                    ("shell_10k_cells", format!("{:.6}", arm.shell_cells)),
                    ("points_outside_shell", arm.outside_with_spread.to_string()),
                    ("outside_shell_population", arm.outside_pop.to_string()),
                    (
                        "inside_shell_with_spread",
                        arm.inside_with_spread.to_string(),
                    ),
                    (
                        "max_spread_outside_shell_world",
                        format!("{:.3e}", arm.max_spread_outside_world),
                    ),
                    ("blended_samples", arm.blended.to_string()),
                    ("active_samples", arm.active.to_string()),
                    ("sign_flips", arm.sign_flips.to_string()),
                    ("topology_changes", arm.topology_changes.to_string()),
                    ("max_spread_world", format!("{:.6e}", arm.max_spread_world)),
                    (
                        "max_spread_cells",
                        format!("{:.9}", arm.max_spread_world / cell),
                    ),
                    ("max_deviation_world", format!("{:.6e}", arm.max_dev_world)),
                    (
                        "max_deviation_cells",
                        format!("{:.6}", arm.max_dev_world / cell),
                    ),
                    (
                        "deviation_bound_kln9_world",
                        format!("{:.6}", arm.k * (OPERANDS as f64).ln()),
                    ),
                    (
                        "deviation_bound_poly_world",
                        format!("{:.6}", arm.k * (OPERANDS - 1) as f64 / 4.0),
                    ),
                    ("min_margin_world", format!("{:.6}", arm.min_margin)),
                    ("min_margin_cells", format!("{:.6}", arm.min_margin / cell)),
                    ("boundary_min_world", format!("{:.6}", arm.boundary_min)),
                    ("grid_cells", GRID_CELLS.to_string()),
                    ("grid_samples", grid_samples.to_string()),
                    ("grid_index_roundtrip", arm.roundtrip.to_string()),
                    ("cell_size_world", format!("{cell:.9}")),
                    ("ref_vertices", arm.ref_vertices.to_string()),
                    ("ref_triangles", arm.ref_triangles.to_string()),
                    ("max_vertices", arm.max_vertices.to_string()),
                    ("stack_agreement", arm.stack_agreement.to_string()),
                    ("spread_floor_world", format!("{SPREAD_FLOOR:.0e}")),
                    ("c1_holds", c1.to_string()),
                    ("c2_holds", c2.to_string()),
                    ("c3_holds", c3.to_string()),
                ]);
            }
        }

        // The M-44 controls this harness needed on top of the registered one.
        assert!(
            any_blend,
            "no arm reached smooth_min's smooth branch, so every zero in this run is vacuous"
        );
        assert!(
            any_distance,
            "no arm reported a Hausdorff distance above 0.01 cells, so the instrument was never \
             shown able to return non-zero"
        );
        assert!(
            any_c3_reachable,
            "no arm had both a non-empty outside-shell population and spread inside the shell, so \
             C3 could not have failed anywhere"
        );

        for row in rows {
            run.record(&row);
        }
    });
}
