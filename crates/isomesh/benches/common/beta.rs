//! Jones' `beta`-numbers: multiscale flatness, per cell and summed over scales.
//!
//! Ticket: R-152, which owns this file. Consumed unchanged by R-153 (`beta`
//! against curvature against camera distance as a refinement criterion) and
//! R-154 (a triangle budget from the Traveling Salesman sum).
//!
//! # The definition
//!
//! For a set `E` and a cube `Q`, Jones' `beta_inf(E, Q)` is
//!
//! ```text
//!     beta_inf(E, Q) = inf over affine planes L of
//!                          sup { dist(x, L) : x in E n Q } / diam(Q)
//! ```
//!
//! — the half-width of the thinnest slab containing the piece of `E` inside
//! `Q`, made scale-free by dividing through by `diam Q`. It is zero exactly
//! when the patch is planar and it grows with the patch's departure from any
//! plane. Jones introduced it for curves in the plane; the surface case, and
//! the theorem that
//!
//! ```text
//!     sum over dyadic Q of beta_inf(Q)^2 * diam(Q)^d
//! ```
//!
//! is finite precisely for `d`-rectifiable sets, is Azzam & Schul's
//! higher-dimensional Traveling Salesman Theorem, `arXiv:1609.02892`. Here
//! `d = 2`, because the thing being measured is a surface. That sum is what
//! [`beta_sum`] accumulates and what gives R-154 a cost estimate computed from
//! the field alone, before any mesh exists.
//!
//! Why this is worth a module: Dual Contouring minimises squared distance to
//! the tangent planes of a cell's crossing points, which is an *assumption*
//! that the patch is nearly planar. `beta_inf(Q)` is the name of exactly that
//! assumption, and it is a quantity, not a hope.
//!
//! # `beta` here is an UPPER BOUND, and the budget is named
//!
//! The exact thinnest slab through `n` points is a minimum-width enclosing-slab
//! problem: the optimum is supported by a small subset of the points and
//! finding it exactly means enumerating candidate supports, which is cubic in
//! `n` at best. This module does not do that. It computes a **feasible** slab
//! and therefore an **upper bound** on `beta_inf`:
//!
//! 1. the total-least-squares plane — the eigenvector of the point set's
//!    covariance matrix belonging to its smallest eigenvalue — gives a first
//!    normal. TLS minimises the *sum of squares* of the residuals, not their
//!    maximum, so it is a good guess and not the answer;
//! 2. that normal is then refined by a deterministic local search:
//!    [`SLAB_ROUNDS`] rounds, each evaluating the max-minus-min objective on a
//!    Fibonacci-spiral spherical cap of [`SLAB_DIRECTIONS`] directions around
//!    the incumbent, the cap's half-angle starting at [`SLAB_CAP_DEGREES`] and
//!    shrinking by [`SLAB_CAP_SHRINK`] each round. The incumbent moves only
//!    between rounds, and ties go to the earlier candidate, so the result is a
//!    pure function of the point set.
//!
//! Total budget: `1 + SLAB_ROUNDS * SLAB_DIRECTIONS` = 73 objective
//! evaluations per patch, each a dot product per point. Every reported `beta`
//! is `>=` the true `beta_inf`, never below it — which is the honest direction
//! for a flatness statistic, since it can only ever over-report roughness and
//! so cannot manufacture a claim of flatness the geometry does not support.
//! Any bench quoting a `beta` must quote the budget with it.
//!
//! # The patch
//!
//! A patch is a point set, and the point set is the one Marching Cubes would
//! place: the field is sampled on a regular grid and every grid edge whose ends
//! disagree in sign contributes its linearly interpolated zero. No extractor is
//! involved, so `beta` is measurable independently of the thing it is meant to
//! explain.
//!
//! - [`beta_per_cell`] uses the extraction grid's own cells, so a cell's patch
//!   is its up-to-twelve edge crossings — byte-for-byte the point set the QEF
//!   is fitted to. That is what makes R-152's correlation against the QEF
//!   residual a comparison of two statistics of one point set.
//! - [`beta_sum`] cuts the box into `n^3` sub-boxes and samples each on
//!   [`SUB_SAMPLES`]`^3` points, up to `3 * (SUB_SAMPLES - 1) * SUB_SAMPLES^2`
//!   = 300 crossings per sub-box. The sub-grids of adjacent sub-boxes share
//!   their common face, so the field is sampled once on a global grid of
//!   `(n * (SUB_SAMPLES - 1) + 1)^3` points — `(4n + 1)^3` `f64`s of memory,
//!   and identical crossing points on shared faces, which a per-sub-box grid
//!   would not give.
//!
//! # `beta` is EXACTLY zero on most cells, for a combinatorial reason
//!
//! Three points are coplanar. A Marching Cubes cell whose sign pattern isolates
//! a single corner has exactly three crossings, so its patch lies in a plane
//! identically and its `beta` is zero — not small, zero to the rounding floor.
//! Measured on this crate's own fields at `33^3`: of the surface cells, 944 of
//! 1160 on `sphere`, 888 of 1128 on `torus`, **all 512** on `thin_plate` and
//! 4101 of 6176 on `noise_cavity` carry a patch that is planar by construction.
//!
//! This is a property of the mechanism and not a defect, but a bench must know
//! it before correlating `beta` against anything: the `beta` column is mostly
//! exact zeros, ranks over it are heavily tied, and a Spearman coefficient
//! computed over the full cell set is dominated by that tie block. Reporting
//! the count of non-planar patches beside the coefficient is the honest form.
//! [`rank_correlation`] averages tied ranks precisely so this stays a defined
//! number rather than an accident of sort order.
//!
//! # Measured cost, so a cost clause is not measuring the wrong thing
//!
//! On a 5900X, release, `f64`, best of five, `sphere` on its own domain:
//!
//! | stage | `65^3` | `129^3` |
//! |---|---|---|
//! | sample the grid + [`SampleGrid::straddles`] over every cell | 2.2 ms | 18.9 ms |
//! | + slab fits on the surface cells | 10.1 ms | 60.6 ms |
//! | Marching Cubes over the same grid, for scale | 4.5 ms | 35.4 ms |
//!
//! So the slab fit is about 2 us per surface cell at the shipped budget, and
//! the whole of [`beta_per_cell`] lands at 1.3x-1.7x an extraction rather than
//! a fraction of it. That number is the mechanism's, not an accident: before
//! the sign-change prepass existed the same walk cost 57 ns on *every* cell and
//! the figure was 4.2x. A cost clause on `beta` must therefore say which cost
//! it means — the standalone pass measured here, or the marginal cost inside an
//! extractor that has already computed the crossings, which is the slab fit
//! alone.
//!
//! What the refinement buys, same fields at `33^3`, over patches with five or
//! more points and a width above the rounding floor: **0.0%** on `sphere`
//! (a symmetric patch's total-least-squares plane already is its min-max
//! plane), mean **5.7%** and up to 21.7% on `torus`, mean **3.5%** and up to
//! 14.5% on `csg_difference`, mean **9.3%** and up to 29.5% on `noise_cavity`.
//! The seventy-two candidates are therefore not decoration on an asymmetric
//! patch, and they are pure overhead on a symmetric one.
//!
//! # Duplication, deliberately
//!
//! The 3x3 symmetric eigensolver below (cyclic Jacobi, forty lines) is also
//! wanted by `common::metric`, which is R-146's file. It is duplicated rather
//! than shared: the two modules are authored concurrently by different tickets,
//! and a dependency between them would be a coupling neither owner could
//! verify. The duplication is scoped to `jacobi_eigen` and to nothing else.

use std::cmp::Ordering;
use std::sync::LazyLock;

use isomesh::{Sdf, Shape3};

/// Samples per axis inside one sub-box of [`beta_sum`]'s dyadic cut.
///
/// Five gives `5^3 = 125` field samples and `3 * 4 * 25 = 300` candidate edges
/// per sub-box, which is enough that a curved patch's extremes are found rather
/// than missed, and small enough that the whole sweep is a few hundred thousand
/// samples at the resolutions the crate meshes at.
pub(crate) const SUB_SAMPLES: u32 = 5;

/// The surface dimension `d` in `sum beta(Q)^2 * diam(Q)^d`.
///
/// Two. The exponent is what makes the sum an area, and it is the reason the
/// sum has the units of a triangle budget rather than of a length.
pub(crate) const SURFACE_DIMENSION: i32 = 2;

/// Candidate directions per refinement round of [`thinnest_slab`].
pub(crate) const SLAB_DIRECTIONS: usize = 24;

/// Refinement rounds of [`thinnest_slab`]. Round zero searches the widest cap;
/// each later round re-centres on the incumbent and shrinks.
pub(crate) const SLAB_ROUNDS: usize = 3;

/// Half-angle, in degrees, of the first refinement cap around the
/// total-least-squares normal.
pub(crate) const SLAB_CAP_DEGREES: f64 = 20.0;

/// Factor by which the cap's half-angle shrinks per round: 20 deg, 5 deg,
/// 1.25 deg.
pub(crate) const SLAB_CAP_SHRINK: f64 = 4.0;

/// Total objective evaluations per patch: the TLS seed plus every candidate.
///
/// Named so a bench can report the search budget beside the `beta` it bought,
/// which is required — an upper bound without its budget is not a measurement.
pub(crate) const SLAB_EVALUATIONS: usize = 1 + SLAB_ROUNDS * SLAB_DIRECTIONS;

/// Jones' `beta_infinity` for a point set inside a cube: the half-width of the
/// thinnest slab containing every point, divided by the cube's diameter.
///
/// Reported as an **upper bound** — see the module header for the search
/// budget. Returns `0.0` for fewer than three points, which fit in a plane
/// trivially, and for a non-positive `diam`, which names no cube.
pub(crate) fn beta_infinity(points: &[[f64; 3]], diam: f64) -> f64 {
    if points.len() < 3 || diam <= 0.0 {
        return 0.0;
    }
    thinnest_slab(points).1 / diam
}

/// The unnormalised slab half-width: `beta_infinity * diam`, in world units.
///
/// This is the quantity R-152 correlates against the QEF residual, because the
/// residual is also in world units and dividing one of them by `diam Q` would
/// make the pair scale-inconsistent. Zero for fewer than three points.
pub(crate) fn beta_times_diam(points: &[[f64; 3]]) -> f64 {
    if points.len() < 3 {
        return 0.0;
    }
    thinnest_slab(points).1
}

/// The thinnest slab found: `(unit normal, half-width, slab centre)`.
///
/// The returned triple is a slab in the literal sense — every point of
/// `points` satisfies `|(p - centre) . normal| <= half_width` — and the normal
/// is the searched direction minimising the maximum absolute signed distance
/// from that direction's best plane, which is the mid-plane between the two
/// extreme points along it. `centre` is therefore the point centroid displaced
/// along `normal` onto the mid-plane, not the bare centroid: the bare centroid
/// would give a self-consistent slab only for a symmetric patch, and would
/// report a needlessly loose bound for every other one.
///
/// The result is an upper bound on the true minimum width (module header). The
/// normal's sign is canonicalised so its largest-magnitude component is
/// non-negative, which leaves the slab unchanged and makes the output
/// reproducible.
///
/// Degenerate input — fewer than three points — yields `([0, 0, 1], 0.0,
/// centroid)`, with `[0, 0, 0]` for the centroid of an empty set.
pub(crate) fn thinnest_slab(points: &[[f64; 3]]) -> ([f64; 3], f64, [f64; 3]) {
    let centroid = centroid(points);
    if points.len() < 3 {
        return ([0.0, 0.0, 1.0], 0.0, centroid);
    }

    // The patch is NOT copied into a centred buffer. Folding the three
    // subtractions into the objective's inner loop costs three flops per point
    // per direction; hoisting them costs one heap allocation per patch, and at
    // twenty thousand surface cells the allocation measured the more expensive
    // of the two. Subtracting the centroid inside the loop also keeps the dot
    // products the size of the patch rather than the size of the patch's
    // distance from the origin, which matters on the domains that reach +-8:
    // a truly planar patch there still reports a width near the rounding
    // floor rather than near `eps * 8`.

    // Step 1: the total-least-squares plane. Its normal is the covariance
    // matrix's eigenvector for the smallest eigenvalue -- the direction of
    // least variance, which is the least-squares answer to a question we are
    // asking in the maximum norm, hence a seed and not a solution.
    let (values, vectors) = jacobi_eigen(covariance(points, centroid));
    let smallest = (0..3)
        .min_by(|&i, &j| values[i].total_cmp(&values[j]))
        .unwrap_or(0);
    let mut best_normal = normalise([
        vectors[0][smallest],
        vectors[1][smallest],
        vectors[2][smallest],
    ]);
    let (mut best_half, mut best_offset) = slab_about(points, centroid, best_normal);

    // Step 2: the documented local search. Each round tilts the incumbent by
    // every direction of a fixed Fibonacci cap; the cap shrinks; ties go to the
    // earlier candidate so the walk is deterministic.
    for cap in SLAB_CAPS.iter() {
        let (e1, e2) = frame(best_normal);
        let mut round_normal = best_normal;
        let mut round_half = best_half;
        let mut round_offset = best_offset;
        for dir in cap {
            let candidate = [
                e1[0] * dir[0] + e2[0] * dir[1] + best_normal[0] * dir[2],
                e1[1] * dir[0] + e2[1] * dir[1] + best_normal[1] * dir[2],
                e1[2] * dir[0] + e2[2] * dir[1] + best_normal[2] * dir[2],
            ];
            let (half, offset) = slab_about(points, centroid, candidate);
            if half < round_half {
                round_normal = candidate;
                round_half = half;
                round_offset = offset;
            }
        }
        best_normal = round_normal;
        best_half = round_half;
        best_offset = round_offset;
    }

    // Canonical sign: the slab is the same object either way, so pick one.
    let widest = (0..3)
        .max_by(|&i, &j| best_normal[i].abs().total_cmp(&best_normal[j].abs()))
        .unwrap_or(0);
    if best_normal[widest] < 0.0 {
        best_normal = [-best_normal[0], -best_normal[1], -best_normal[2]];
        best_offset = -best_offset;
    }

    let centre = [
        centroid[0] + best_normal[0] * best_offset,
        centroid[1] + best_normal[1] * best_offset,
        centroid[2] + best_normal[2] * best_offset,
    ];
    (best_normal, best_half, centre)
}

/// The Traveling Salesman sum, `sum over Q of beta(Q)^2 * diam(Q)^d`,
/// accumulated over dyadic scales.
///
/// `d` is [`SURFACE_DIMENSION`]. The theorem's content is that this converges
/// for a rectifiable surface and diverges otherwise, so the convergence column
/// is a statement about the field, not about the harness.
#[derive(Clone, Debug, Default)]
pub(crate) struct BetaSum {
    /// The running sum over every scale accumulated so far.
    pub(crate) total: f64,
    /// One entry per scale, in the order the scales were given:
    /// `(n, this scale's contribution, sub-boxes that carried a patch)`, where
    /// `n` is the scale as passed to [`beta_sum`] — the number of sub-boxes per
    /// axis, so the level holds `n^3` of them.
    pub(crate) per_scale: Vec<(u32, f64, u64)>,
}

impl BetaSum {
    /// Does the sum converge as scales are added? True when the last scale's
    /// increment is under `rel_tol` of the running total.
    ///
    /// This is the observable form of the theorem's hypothesis: a rectifiable
    /// surface's dyadic contributions must fall away, so a finite tail is the
    /// signature of convergence and a contribution that stays a fixed fraction
    /// of the total at every level is the signature of divergence. False with
    /// fewer than two scales, because one increment is the whole sum and
    /// nothing has been observed, and false for a non-positive or non-finite
    /// total, where the ratio means nothing.
    pub(crate) fn converges(&self, rel_tol: f64) -> bool {
        if self.per_scale.len() < 2 || !self.total.is_finite() || self.total <= 0.0 {
            return false;
        }
        match self.per_scale.last() {
            Some(&(_, last, _)) => last <= rel_tol * self.total,
            None => false,
        }
    }

    /// How many scales were accumulated. R-154's vacuity control requires at
    /// least four, or convergence cannot be observed.
    pub(crate) fn scales_used(&self) -> usize {
        self.per_scale.len()
    }
}

/// Accumulate the beta-sum over `scales` dyadic resolutions of `sdf` on the box
/// `(lo, hi)`.
///
/// At each scale `n` the box is cut into `n^3 ` sub-boxes; a sub-box's point set
/// is the surface samples inside it, found by marching the field on a shared
/// sub-grid of [`SUB_SAMPLES`] samples per axis per sub-box (module header).
/// `diam(Q)` is the sub-box's true space diagonal, so a non-cubic box is
/// handled correctly and a cubic one gives the familiar `side * sqrt(3)`.
///
/// A sub-box with fewer than three crossings contributes nothing and is not
/// counted: it holds no patch a plane could be fitted to. Scales are visited in
/// the order given, which is the order `per_scale` reports and the order
/// [`BetaSum::converges`] reads.
pub(crate) fn beta_sum<S: Sdf<Scalar = f64>>(
    sdf: &S,
    lo: [f64; 3],
    hi: [f64; 3],
    scales: &[u32],
) -> BetaSum {
    let mut out = BetaSum::default();
    let span = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];

    for &n in scales {
        if n == 0 {
            continue;
        }
        let stride = usize::try_from(SUB_SAMPLES - 1).unwrap_or(1);
        let per_axis = usize::try_from(n).unwrap_or(1) * stride + 1;
        let divisions = f64::from(n);
        let step = [
            span[0] / (divisions * stride as f64),
            span[1] / (divisions * stride as f64),
            span[2] / (divisions * stride as f64),
        ];
        let grid = SampleGrid::sampled(sdf, lo, step, [per_axis; 3]);

        // The sub-box diameter at this scale: the space diagonal of one cut.
        let side = [
            span[0] / divisions,
            span[1] / divisions,
            span[2] / divisions,
        ];
        let diam = (side[0] * side[0] + side[1] * side[1] + side[2] * side[2]).sqrt();
        let weight = diam.powi(SURFACE_DIMENSION);

        let block = [stride + 1; 3];
        let mut points = Vec::with_capacity(300);
        let mut partial = 0.0;
        let mut counted = 0u64;
        let cuts = usize::try_from(n).unwrap_or(1);
        for bz in 0..cuts {
            for by in 0..cuts {
                for bx in 0..cuts {
                    let base = [bx * stride, by * stride, bz * stride];
                    if !grid.straddles(base, block) {
                        continue;
                    }
                    points.clear();
                    grid.crossings(base, block, &mut points);
                    if points.len() < 3 {
                        continue;
                    }
                    counted += 1;
                    let beta = thinnest_slab(&points).1 / diam;
                    partial += beta * beta * weight;
                }
            }
        }
        out.total += partial;
        out.per_scale.push((n, partial, counted));
    }
    out
}

/// Per-cell beta over one grid: one entry per cell in x-fastest order, `None`
/// where the cell holds no surface.
///
/// `shape` counts **samples**, so an `n`-sample axis spans `n - 1` cells and the
/// returned vector has `(sx - 1) * (sy - 1) * (sz - 1)` entries indexed
/// `cx + cy * (sx - 1) + cz * (sx - 1) * (sy - 1)` — the same x-fastest
/// convention [`Shape3`] uses for samples. An axis with fewer than two samples
/// spans no cells and yields an empty vector.
///
/// A cell's patch is its own edge crossings, at most twelve, which is exactly
/// the point set Dual Contouring fits its plane to. `None` marks a cell with
/// fewer than three crossings; because a closed sign pattern on a cube produces
/// at least three crossings whenever it produces any, that is the same set of
/// cells as "no surface here".
///
/// Expect most `Some` values to be **exactly zero**: a three-crossing cell's
/// patch is three points and three points are coplanar. See the module header
/// for the measured proportion per field — it is the majority on every one of
/// them — and for why a correlation against this column must report the
/// non-planar count beside the coefficient.
pub(crate) fn beta_per_cell<S: Sdf<Scalar = f64>>(
    sdf: &S,
    shape: &impl Shape3,
    origin: [f64; 3],
    cell_size: f64,
) -> Vec<Option<f64>> {
    let size = shape.size();
    if size[0] < 2 || size[1] < 2 || size[2] < 2 {
        return Vec::new();
    }
    let dim = [
        usize::try_from(size[0]).unwrap_or(0),
        usize::try_from(size[1]).unwrap_or(0),
        usize::try_from(size[2]).unwrap_or(0),
    ];
    let grid = SampleGrid::sampled(sdf, origin, [cell_size; 3], dim);

    // A cube of side h has diameter h*sqrt(3); that is the normaliser Jones'
    // definition divides by, and it is what makes beta dimensionless.
    let diam = cell_size * 3.0_f64.sqrt();
    let cells = [dim[0] - 1, dim[1] - 1, dim[2] - 1];
    let mut out = Vec::with_capacity(cells[0] * cells[1] * cells[2]);
    let mut points = Vec::with_capacity(12);
    for cz in 0..cells[2] {
        for cy in 0..cells[1] {
            for cx in 0..cells[0] {
                if diam <= 0.0 || !grid.straddles([cx, cy, cz], [2; 3]) {
                    out.push(None);
                    continue;
                }
                points.clear();
                grid.crossings([cx, cy, cz], [2; 3], &mut points);
                if points.len() < 3 {
                    out.push(None);
                } else {
                    out.push(Some(thinnest_slab(&points).1 / diam));
                }
            }
        }
    }
    out
}

/// Spearman rank correlation of two equal-length samples, ties averaged.
///
/// Spearman rather than Pearson because R-152 asks whether `beta` and the QEF
/// residual *order* cells the same way, not whether they are linearly related:
/// the two have different units and no reason to be affine in one another.
/// Ranks are assigned by [`f64::total_cmp`], so the result does not depend on
/// the input's order or on how the platform compares NaN. Returns `0.0` for
/// fewer than three pairs and for a sample with no rank variance at all, where
/// the correlation would be against a constant.
///
/// # Panics
///
/// If the two samples differ in length; a bench that mismatched them has
/// paired the wrong columns and no correlation of the two is meaningful.
pub(crate) fn rank_correlation(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(
        a.len(),
        b.len(),
        "rank_correlation needs paired samples of equal length"
    );
    if a.len() < 3 {
        return 0.0;
    }
    let (ra, rb) = (ranks(a), ranks(b));
    pearson(&ra, &rb)
}

// ── the patch: field samples, and the crossings Marching Cubes would place ──

/// A scalar field sampled on an axis-aligned regular grid, x-fastest.
///
/// One grid serves a whole scale: sub-boxes read overlapping blocks of it, so
/// their shared faces carry identical crossing points, which is what makes the
/// dyadic sum a sum over a tiling rather than over 8/27 disagreeing copies.
#[derive(Clone, Debug)]
struct SampleGrid {
    values: Vec<f64>,
    dim: [usize; 3],
    origin: [f64; 3],
    step: [f64; 3],
}

impl SampleGrid {
    /// Sample `sdf` on `dim` points per axis starting at `origin`, spacing
    /// `step`. One pass, x innermost, which is the crate's own convention.
    fn sampled<S: Sdf<Scalar = f64>>(
        sdf: &S,
        origin: [f64; 3],
        step: [f64; 3],
        dim: [usize; 3],
    ) -> Self {
        let mut values = Vec::with_capacity(dim[0] * dim[1] * dim[2]);
        for z in 0..dim[2] {
            let pz = origin[2] + step[2] * z as f64;
            for y in 0..dim[1] {
                let py = origin[1] + step[1] * y as f64;
                for x in 0..dim[0] {
                    let px = origin[0] + step[0] * x as f64;
                    values.push(sdf.sample([px, py, pz]));
                }
            }
        }
        Self {
            values,
            dim,
            origin,
            step,
        }
    }

    /// Field value at a sample index.
    fn value(&self, p: [usize; 3]) -> f64 {
        self.values[p[0] + self.dim[0] * (p[1] + self.dim[1] * p[2])]
    }

    /// World position of a sample index.
    fn position(&self, p: [usize; 3]) -> [f64; 3] {
        [
            self.origin[0] + self.step[0] * p[0] as f64,
            self.origin[1] + self.step[1] * p[1] as f64,
            self.origin[2] + self.step[2] * p[2] as f64,
        ]
    }

    /// Does the block of `span` samples per axis at `base` straddle the
    /// isosurface?
    ///
    /// Exact, not conservative: a grid block is edge-connected, so if its
    /// samples do not all share a sign then some path from a negative sample to
    /// a positive one crosses an edge of the block, and that edge produces a
    /// crossing. Contrapositively, one sign everywhere means no crossing
    /// anywhere.
    ///
    /// This is the prepass that makes `beta` affordable. Without it every cell
    /// of the grid walks its twelve edges — two loads and a comparison each —
    /// where the overwhelming majority of cells hold no surface at all. At
    /// `129^3` the walk cost 57 ns per cell over 2.1 M cells, which is three
    /// times the whole of Marching Cubes; the eight-corner test answers the
    /// same question in eight loads and quits on the first disagreement.
    fn straddles(&self, base: [usize; 3], span: [usize; 3]) -> bool {
        let first = self.value(base) < 0.0;
        for k in 0..span[2] {
            let z = base[2] + k;
            if z >= self.dim[2] {
                break;
            }
            for j in 0..span[1] {
                let y = base[1] + j;
                if y >= self.dim[1] {
                    break;
                }
                let row = self.dim[0] * (y + self.dim[1] * z);
                let hi = (base[0] + span[0]).min(self.dim[0]);
                for value in &self.values[row + base[0]..row + hi] {
                    if (*value < 0.0) != first {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// Append the surface crossings of the block of `span` samples per axis
    /// whose low corner is the sample `base`.
    ///
    /// Every axis-aligned edge inside the block whose ends disagree in sign
    /// contributes its linearly interpolated zero — the same point Marching
    /// Cubes places on that edge, with the same `t = a / (a - b)`. A sample of
    /// exactly zero counts as outside, so it produces no crossing on its own
    /// edges; that is the `=`-corner convention, and it matters only on a
    /// measure-zero set of fields.
    fn crossings(&self, base: [usize; 3], span: [usize; 3], out: &mut Vec<[f64; 3]>) {
        for axis in 0..3 {
            let mut delta = [0usize; 3];
            delta[axis] = 1;
            let mut extent = span;
            extent[axis] -= 1;
            for k in 0..extent[2] {
                for j in 0..extent[1] {
                    for i in 0..extent[0] {
                        let a = [base[0] + i, base[1] + j, base[2] + k];
                        let b = [a[0] + delta[0], a[1] + delta[1], a[2] + delta[2]];
                        if b[axis] >= self.dim[axis] {
                            continue;
                        }
                        let (va, vb) = (self.value(a), self.value(b));
                        if (va < 0.0) == (vb < 0.0) {
                            continue;
                        }
                        let denominator = va - vb;
                        let t = if denominator.abs() > 0.0 {
                            (va / denominator).clamp(0.0, 1.0)
                        } else {
                            0.5
                        };
                        let (pa, pb) = (self.position(a), self.position(b));
                        out.push([
                            pa[0] + (pb[0] - pa[0]) * t,
                            pa[1] + (pb[1] - pa[1]) * t,
                            pa[2] + (pb[2] - pa[2]) * t,
                        ]);
                    }
                }
            }
        }
    }
}

// ── the slab objective and its deterministic direction search ──

/// Centroid of a point set; the origin for an empty one.
fn centroid(points: &[[f64; 3]]) -> [f64; 3] {
    if points.is_empty() {
        return [0.0; 3];
    }
    let mut sum = [0.0; 3];
    for p in points {
        sum[0] += p[0];
        sum[1] += p[1];
        sum[2] += p[2];
    }
    let n = points.len() as f64;
    [sum[0] / n, sum[1] / n, sum[2] / n]
}

/// Covariance matrix of the points about `centre`, unnormalised — the scale
/// factor cannot move an eigenvector, so dividing by `n` would only cost a
/// division.
fn covariance(points: &[[f64; 3]], centre: [f64; 3]) -> [[f64; 3]; 3] {
    let mut m = [[0.0; 3]; 3];
    for p in points {
        let d = [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]];
        for r in 0..3 {
            for c in 0..3 {
                m[r][c] += d[r] * d[c];
            }
        }
    }
    m
}

/// The objective: `(half-width, offset)` of the thinnest slab with normal `n`.
///
/// `offset` is where that slab's mid-plane sits along `n` relative to `centre`.
/// Taking `(max - min) / 2` about the mid-plane rather than `max |s|` about
/// `centre` is what makes this the best slab *for this direction*, so the
/// search compares directions and not centrings.
///
/// `min`/`max` rather than a branch per comparison: this runs
/// [`SLAB_EVALUATIONS`] times per patch, the branch is unpredictable by
/// construction, and the pair compiles to two instructions with no branch at
/// all — which on a patch's worth of points is a third of the cost.
fn slab_about(points: &[[f64; 3]], centre: [f64; 3], n: [f64; 3]) -> (f64, f64) {
    if points.is_empty() {
        return (0.0, 0.0);
    }
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for p in points {
        let s = (p[0] - centre[0]) * n[0] + (p[1] - centre[1]) * n[1] + (p[2] - centre[2]) * n[2];
        lo = lo.min(s);
        hi = hi.max(s);
    }
    ((hi - lo) / 2.0, (hi + lo) / 2.0)
}

/// The refinement caps, in a canonical frame whose `z` is the incumbent normal.
///
/// Computed once and shared: they are a pure function of the three budget
/// constants, so recomputing them per patch would buy nothing but trig. Round
/// `r` samples the cap of half-angle `SLAB_CAP_DEGREES / SLAB_CAP_SHRINK^r`
/// with a Fibonacci spiral — `cos(theta)` uniform over the cap, azimuth
/// advancing by the golden angle — which spreads the directions without a
/// random number in sight.
static SLAB_CAPS: LazyLock<[[[f64; 3]; SLAB_DIRECTIONS]; SLAB_ROUNDS]> = LazyLock::new(|| {
    // 1 / phi, the golden-ratio conjugate: the azimuth increment that keeps
    // consecutive samples maximally apart.
    const GOLDEN: f64 = 0.618_033_988_749_894_9;
    let mut caps = [[[0.0; 3]; SLAB_DIRECTIONS]; SLAB_ROUNDS];
    let mut half_angle = SLAB_CAP_DEGREES.to_radians();
    for cap in &mut caps {
        let cos_min = half_angle.cos();
        for (k, dir) in cap.iter_mut().enumerate() {
            let fraction = (k as f64 + 0.5) / SLAB_DIRECTIONS as f64;
            let cos_theta = 1.0 - fraction * (1.0 - cos_min);
            let radius = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
            let azimuth = std::f64::consts::TAU * (k as f64 * GOLDEN).fract();
            *dir = [radius * azimuth.cos(), radius * azimuth.sin(), cos_theta];
        }
        half_angle /= SLAB_CAP_SHRINK;
    }
    caps
});

/// Two unit vectors completing `n` to a right-handed orthonormal frame.
///
/// The seed axis is the coordinate axis least aligned with `n`, whose component
/// is at most `1 / sqrt(3)` in magnitude, so the cross product's length is at
/// least `sqrt(2/3)` and the normalisation is never near-singular.
fn frame(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let axis = (0..3)
        .min_by(|&i, &j| n[i].abs().total_cmp(&n[j].abs()))
        .unwrap_or(0);
    let mut seed = [0.0; 3];
    seed[axis] = 1.0;
    let e1 = normalise(cross(n, seed));
    let e2 = cross(n, e1);
    (e1, e2)
}

/// Cross product.
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Unit vector along `v`; `+z` for a vector with no length, so that a
/// degenerate point set still yields a well-defined slab rather than a NaN.
fn normalise(v: [f64; 3]) -> [f64; 3] {
    let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if length > 0.0 && length.is_finite() {
        [v[0] / length, v[1] / length, v[2] / length]
    } else {
        [0.0, 0.0, 1.0]
    }
}

/// Eigenvalues and eigenvectors of a symmetric 3x3 matrix, by cyclic Jacobi.
///
/// Returns `(values, vectors)` with `vectors[row][k]` the `row`-th component of
/// the eigenvector belonging to `values[k]`. Jacobi rather than the closed-form
/// cubic because the cubic's discriminant loses most of its digits on a nearly
/// planar patch, which is precisely the case `beta` is computed for: a patch
/// whose smallest eigenvalue is `1e-12` of its largest is the flat cell, and
/// the flat cell is the one whose normal must still come out right.
///
/// Deliberately duplicated from `common::metric` — see the module header.
fn jacobi_eigen(mut a: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut v = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    // Six sweeps of three rotations is far beyond convergence for 3x3: each
    // sweep squares the off-diagonal norm, so this terminates on the tolerance
    // long before the count, and the count only bounds a pathological input.
    for _ in 0..6 {
        let off = a[0][1].abs() + a[0][2].abs() + a[1][2].abs();
        let scale = a[0][0].abs() + a[1][1].abs() + a[2][2].abs();
        if off <= f64::EPSILON * scale {
            break;
        }
        for &(p, q) in &[(0usize, 1usize), (0, 2), (1, 2)] {
            if a[p][q].abs() <= f64::EPSILON * scale {
                continue;
            }
            // The rotation that zeroes a[p][q]: theta parameterises the pivot
            // ratio, and taking the smaller root of t^2 + 2*theta*t - 1 keeps
            // |t| <= 1 and the rotation numerically benign.
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = if theta >= 0.0 {
                1.0 / (theta + (theta * theta + 1.0).sqrt())
            } else {
                -1.0 / (-theta + (theta * theta + 1.0).sqrt())
            };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            let (app, aqq, apq) = (a[p][p], a[q][q], a[p][q]);
            a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
            a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
            a[p][q] = 0.0;
            a[q][p] = 0.0;
            let r = 3 - p - q; // the row this rotation leaves alone
            let (arp, arq) = (a[r][p], a[r][q]);
            a[r][p] = c * arp - s * arq;
            a[p][r] = a[r][p];
            a[r][q] = s * arp + c * arq;
            a[q][r] = a[r][q];
            for row in &mut v {
                let (vp, vq) = (row[p], row[q]);
                row[p] = c * vp - s * vq;
                row[q] = s * vp + c * vq;
            }
        }
    }
    ([a[0][0], a[1][1], a[2][2]], v)
}

// ── Spearman's two halves ──

/// Ranks of `values`, one-based, tied values sharing the mean of their block.
///
/// Ordered by [`f64::total_cmp`] and broken by index, so the ranking is a pure
/// function of the slice's contents and position — no NaN-dependent branch and
/// no reliance on the sort's stability.
fn ranks(values: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&i, &j| values[i].total_cmp(&values[j]).then(i.cmp(&j)));
    let mut out = vec![0.0; values.len()];
    let mut start = 0;
    while start < order.len() {
        let mut end = start + 1;
        while end < order.len()
            && values[order[end]].total_cmp(&values[order[start]]) == Ordering::Equal
        {
            end += 1;
        }
        // Mean of the one-based ranks start+1 ..= end.
        let mean = (start + 1 + end) as f64 / 2.0;
        for &i in &order[start..end] {
            out[i] = mean;
        }
        start = end;
    }
    out
}

/// Pearson correlation; `0.0` when either sample has no variance, which for
/// ranks means every value tied and no ordering to correlate.
fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let mean_a = a.iter().sum::<f64>() / n;
    let mean_b = b.iter().sum::<f64>() / n;
    let mut sab = 0.0;
    let mut saa = 0.0;
    let mut sbb = 0.0;
    for (x, y) in a.iter().zip(b.iter()) {
        let dx = x - mean_a;
        let dy = y - mean_b;
        sab += dx * dy;
        saa += dx * dx;
        sbb += dy * dy;
    }
    let denominator = (saa * sbb).sqrt();
    if denominator > 0.0 {
        sab / denominator
    } else {
        0.0
    }
}
