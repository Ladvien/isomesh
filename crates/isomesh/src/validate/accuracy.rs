//! How far a mesh sits from the surface it claims to represent.
//!
//! [`validate`](super::validate) answers "is this mesh well formed"; this module
//! answers "is it in the right place". They are independent: a mesh can be
//! perfectly manifold, correctly oriented, closed, and describe the wrong
//! surface entirely.
//!
//! # What is measured
//!
//! Both one-sided distances, and the symmetric Hausdorff distance as the larger
//! of the two:
//!
//! - **mesh → field**, sampled at every referenced vertex and every triangle
//!   centroid, which is the set the accuracy ticket names.
//! - **field → mesh**, sampled by projecting the extractor's own lattice points
//!   onto the surface. This is the direction that notices *missing* geometry;
//!   the forward direction only notices *misplaced* geometry.
//!
//! Both are sampled, so both are **lower bounds** on the true Hausdorff
//! distance. The same limitation applies to Metro (Cignoni, Rocchini &
//! Scopigno, *Metro: measuring error on simplified surfaces*, Computer Graphics
//! Forum 17(2), 1998, `10.1111/1467-8659.00236`), which is the standard tool for
//! this measurement and samples for the same reason: the exact quantity requires
//! a continuous supremum.
//!
//! # Why distance is measured geometrically and not from `|f(p)|`
//!
//! One unconditional path, for all seven reference fields and any other `Sdf`.
//!
//! `|f(p)|` is a distance only when `|∇f| == 1`, which is false for two of the
//! seven reference fields. Branching on
//! [`is_exact_distance`](crate::fields::ReferenceField::is_exact_distance) would
//! put two execution paths in the harness, and the crate's rule is one.
//!
//! The obvious branch-free repair, `|f(p)|/|∇f(p)|`, is rejected for a reason
//! the field module already states: dividing by `|∇f|` "introduces a second
//! failure mode wherever `|∇g| → 0`". It is also only first order.
//!
//! So the distance reported is the length of the **gradient-flow chord**: project
//! the sample onto the surface by Newton iteration along `∇f`, and measure how
//! far it moved.
//!
//! ```text
//! p ← p − f(p)·∇f(p)/|∇f(p)|²        d(p) = |p − p_projected|
//! ```
//!
//! This is exact for the five fields that are true distance fields — the first
//! step lands on the surface and its length is `|f(p)|` — and it iterates to the
//! same accuracy on the other two. Note that the first step's length *is*
//! `|f|/|∇f|`, so the rejected estimate is this routine's first iteration; the
//! difference is that this one keeps going. The field module's own instruction
//! is followed rather than worked around: "an accuracy harness should measure
//! against geometry rather than against `|sample|`".
//!
//! Whenever the projection converges the projected point is on the surface, so
//! `d(p) ≥ dist(p, S)`: the reported error is an **over**-estimate, which makes
//! "max error below X" a conservative pass rather than an optimistic one. Near a
//! concave seam — `csg_difference` is the case — the flow can land further away
//! than the true nearest point, and that is the direction the bias runs.
//!
//! # Cost
//!
//! Dominated by the seed band filter, which evaluates the field and its gradient
//! once per lattice point. Like [`validate_indexed`](super::validate_indexed),
//! this is a measurement and not a hot path: it allocates freely.

use alloc::vec::Vec;
use core::fmt;

use super::ValidateConfig;
use super::tri_grid::TriangleGrid;
use crate::{Real, Sdf, Shape3, vec3};

/// Thresholds and iteration limits for an accuracy measurement.
///
/// Private fields and one checked constructor, for the reason
/// [`ValidateConfig`] gives: every threshold here is relative to the grid
/// spacing, so a config that exists is one whose thresholds mean something.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AccuracyConfig {
    cell_size: f64,
    band_radius: f64,
    residual_tolerance: f64,
}

impl AccuracyConfig {
    /// Seed band half-width, relative to `cell_size`.
    ///
    /// A surface point's furthest cell corner is `h·√3/2 ≈ 0.87h` away, so `1.0`
    /// guarantees that every surface-crossing cell contributes at least one
    /// seed, with margin for fields whose `|∇f|` is not 1.
    pub const BAND_RADIUS_REL: f64 = 1.0;

    /// Newton is converged when its step is shorter than this, times `cell_size`.
    ///
    /// Relative to the spacing rather than absolute, so it scales with the
    /// quantity being measured. At `1e-4·h` it sits two orders below the errors
    /// an extractor actually produces (`~1e-3·h`) and still well above the `f32`
    /// rounding floor at ordinary coordinate magnitudes, so it is reachable in
    /// single precision rather than stagnating.
    pub const RESIDUAL_TOLERANCE_REL: f64 = 1e-4;

    /// Newton iteration cap.
    ///
    /// One step is exact for a true distance field, and a second confirms it, so
    /// five of the seven reference fields converge in two. The cap exists for
    /// the other two and, more importantly, to bound the cost of points that
    /// never converge — a field critical point, or a concave seam where the
    /// iteration oscillates between two operands.
    pub const MAX_NEWTON_ITERATIONS: u32 = 8;

    /// The only constructor.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidCellSize`](crate::Error::InvalidCellSize) if `cell_size`
    /// is not finite and positive.
    pub fn from_cell_size(cell_size: f64) -> crate::Result<Self> {
        if !cell_size.is_finite() || cell_size <= 0.0 {
            return Err(crate::Error::InvalidCellSize { value: cell_size });
        }
        Ok(Self {
            cell_size,
            band_radius: cell_size * Self::BAND_RADIUS_REL,
            residual_tolerance: cell_size * Self::RESIDUAL_TOLERANCE_REL,
        })
    }

    /// Grid spacing `h`, in world units.
    #[must_use]
    pub fn cell_size(&self) -> f64 {
        self.cell_size
    }

    /// Lattice points further than this from the surface are not seeded.
    #[must_use]
    pub fn band_radius(&self) -> f64 {
        self.band_radius
    }

    /// Newton stops when its step falls to this length.
    #[must_use]
    pub fn residual_tolerance(&self) -> f64 {
        self.residual_tolerance
    }

    /// `cell_size · √3`, the diagonal of one cell.
    ///
    /// Provided so that the usual "within one cell diagonal" acceptance test is
    /// written once here rather than re-derived at each call site.
    #[must_use]
    pub fn cell_diagonal(&self) -> f64 {
        self.cell_size * libm::sqrt(3.0)
    }
}

/// Max and mean over one set of distance samples.
///
/// Generic over `R` because [`Real`] offers no widening to `f64`, so a distance
/// cannot be stored more precisely than the mesh it came from.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DistanceStats<R: Real> {
    /// Samples that produced a distance.
    pub samples: u64,
    /// The largest distance found. Zero when `samples == 0`.
    pub max: R,
    /// The mean distance. Every distance here is non-negative, so this is the
    /// mean *absolute* error over this set. Zero when `samples == 0`.
    pub mean: R,
}

impl<R: Real> DistanceStats<R> {
    const EMPTY: Self = Self {
        samples: 0,
        max: R::ZERO,
        mean: R::ZERO,
    };
}

/// What an accuracy measurement found.
///
/// Every field is a count or a distance; nothing here asserts. The loud path is
/// the opt-in [`panic_if_worse_than`](Self::panic_if_worse_than).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AccuracyReport<R: Real> {
    /// Mesh to field: every referenced vertex and every usable centroid.
    pub mesh_to_field: DistanceStats<R>,
    /// Field to mesh: lattice seeds projected onto the surface.
    pub field_to_mesh: DistanceStats<R>,

    /// Referenced vertices sampled.
    pub vertex_samples: u64,
    /// Triangle centroids sampled. Equal to `triangles`.
    pub centroid_samples: u64,
    /// Mesh samples whose projection did not converge, excluded from
    /// `mesh_to_field`. A property of the field at that point, not of the mesh.
    pub unconverged_mesh_samples: u64,

    /// Lattice points considered. The five buckets below plus
    /// `field_to_mesh.samples` sum to exactly this.
    pub seeds: u64,
    /// Seeds whose first Newton step exceeded the band radius: too far from the
    /// surface to be worth projecting.
    pub seeds_out_of_band: u64,
    /// Seeds at a field critical point — `∇f` zero or not finite.
    pub seeds_non_finite: u64,
    /// Seeds that hit the iteration cap.
    pub seeds_unconverged: u64,
    /// Seeds that projected to a point outside the measurable interior.
    pub seeds_outside_domain: u64,

    /// Triangles that passed the usability filter.
    pub triangles: u64,
    /// Faces rejected for an out-of-range or repeated index.
    pub faces_skipped: u64,
    /// Faces rejected as degenerate. Reported, not an error: Marching Cubes
    /// emits slivers whenever a corner value sits near zero.
    pub degenerate_triangles: u64,
    /// Faces rejected for a non-finite coordinate.
    pub non_finite_positions: u64,

    /// The config this was measured with.
    pub config: AccuracyConfig,
}

impl<R: Real> AccuracyReport<R> {
    /// `max(mesh_to_field.max, field_to_mesh.max)`.
    ///
    /// Sampled on both sides, so a lower bound on the true symmetric Hausdorff
    /// distance.
    #[must_use]
    pub fn symmetric_hausdorff(&self) -> R {
        self.mesh_to_field.max.max(self.field_to_mesh.max)
    }

    /// `mesh_to_field.mean`.
    ///
    /// Deliberately the forward set and not a blend of the two. Its sampling
    /// density is a property of the mesh alone, whereas the reverse set's size
    /// depends on [`AccuracyConfig::BAND_RADIUS_REL`] and the seed lattice — a
    /// mean over it would move when a tuning constant moved, which is useless as
    /// a regression number. "Mean of the symmetric distance" is not a defined
    /// quantity, so it is not invented here; both means are reported separately.
    #[must_use]
    pub fn mean_absolute_error(&self) -> R {
        self.mesh_to_field.mean
    }

    /// Both directions produced samples.
    ///
    /// False for an empty mesh, and false when no seed survived — which is
    /// itself the signal that the extractor missed the surface entirely.
    #[must_use]
    pub fn has_coverage(&self) -> bool {
        self.mesh_to_field.samples > 0 && self.field_to_mesh.samples > 0
    }

    /// Panic unless the mesh is within `limit` of the surface.
    ///
    /// Opt in, mirroring
    /// [`MeshReport::panic_if_invalid`](super::MeshReport::panic_if_invalid).
    ///
    /// # Panics
    ///
    /// If [`has_coverage`](Self::has_coverage) is false, or if
    /// [`symmetric_hausdorff`](Self::symmetric_hausdorff) exceeds `limit`.
    pub fn panic_if_worse_than(&self, limit: R) {
        assert!(self.has_coverage(), "no accuracy coverage\n{self}");
        assert!(
            self.symmetric_hausdorff() <= limit,
            "mesh is further than {limit:?} from the surface\n{self}"
        );
    }
}

impl<R: Real> fmt::Display for AccuracyReport<R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "accuracy")?;
        writeln!(f, "  h              {:?}", self.config.cell_size())?;
        writeln!(f, "  ---------------------------")?;
        writeln!(
            f,
            "  mesh->field    max {:?}  mean {:?}  n {}",
            self.mesh_to_field.max, self.mesh_to_field.mean, self.mesh_to_field.samples
        )?;
        writeln!(
            f,
            "  field->mesh    max {:?}  mean {:?}  n {}",
            self.field_to_mesh.max, self.field_to_mesh.mean, self.field_to_mesh.samples
        )?;
        writeln!(f, "  symmetric      {:?}", self.symmetric_hausdorff())?;
        writeln!(f, "  ---------------------------")?;
        writeln!(
            f,
            "  mesh samples   {} vertices + {} centroids, {} unconverged",
            self.vertex_samples, self.centroid_samples, self.unconverged_mesh_samples
        )?;
        writeln!(
            f,
            "  seeds          {} total: {} out of band, {} non-finite, {} unconverged, {} outside",
            self.seeds,
            self.seeds_out_of_band,
            self.seeds_non_finite,
            self.seeds_unconverged,
            self.seeds_outside_domain
        )?;
        writeln!(
            f,
            "  triangles      {} used, {} skipped, {} degenerate, {} non-finite",
            self.triangles,
            self.faces_skipped,
            self.degenerate_triangles,
            self.non_finite_positions
        )?;
        if self.has_coverage() {
            write!(f, "  measured both directions")
        } else {
            write!(
                f,
                "  !! NO COVERAGE — one or both directions sampled nothing"
            )
        }
    }
}

/// Measure how far an indexed mesh sits from an analytic field.
///
/// `shape`, `origin` and `cfg.cell_size()` describe the **seed lattice** for the
/// reverse direction. Normally they are the same grid the mesh was extracted on,
/// which is why the parameter list mirrors
/// [`MarchingCubes::extract`](crate::mc::MarchingCubes::extract) — the caller
/// already has those values in hand. Nothing requires it: passing a coarser
/// shape samples the reverse direction more cheaply.
///
/// There is no `cell_size` parameter because the spacing lives in `cfg` and
/// nowhere else. One source of truth makes a spacing/threshold mismatch
/// unrepresentable rather than merely unlikely.
///
/// # Errors
///
/// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if any axis of `shape`
/// has fewer than two samples.
///
/// [`Error::CellSizeMismatch`](crate::Error::CellSizeMismatch) if the spacing
/// does not describe this mesh: one triangle spanning more than 512 grid cells,
/// or a grid exceeding 2²² cells in total. Both mean `cfg.cell_size()` and this
/// mesh do not belong to each other.
pub fn accuracy<R, S>(
    positions: &[[R; 3]],
    indices: &[u32],
    sdf: &S,
    shape: &impl Shape3,
    origin: [R; 3],
    cfg: &AccuracyConfig,
) -> crate::Result<AccuracyReport<R>>
where
    R: Real,
    S: Sdf<Scalar = R>,
{
    let size = shape.size();
    if size.iter().any(|&n| n < 2) {
        return Err(crate::Error::GridTooSmall { size });
    }

    let h = R::from_f64(cfg.cell_size());
    let band = R::from_f64(cfg.band_radius());
    let tol = R::from_f64(cfg.residual_tolerance());

    let mut report = AccuracyReport {
        mesh_to_field: DistanceStats::EMPTY,
        field_to_mesh: DistanceStats::EMPTY,
        vertex_samples: 0,
        centroid_samples: 0,
        unconverged_mesh_samples: 0,
        seeds: 0,
        seeds_out_of_band: 0,
        seeds_non_finite: 0,
        seeds_unconverged: 0,
        seeds_outside_domain: 0,
        triangles: 0,
        faces_skipped: 0,
        degenerate_triangles: 0,
        non_finite_positions: 0,
        config: *cfg,
    };

    // ── the usable triangle set ─────────────────────────────────────────────
    //
    // One filter, applied once. Its output feeds both the centroid sample set
    // and the spatial index, so the two cannot disagree about what a triangle
    // is.
    let two_area_limit =
        R::from_f64(2.0 * ValidateConfig::AREA_EPSILON_REL * cfg.cell_size() * cfg.cell_size());
    let limit_sq = two_area_limit * two_area_limit;

    let whole = indices.len() - indices.len() % 3;
    // A trailing partial group is one face that cannot be read.
    report.faces_skipped = u64::from(indices.len() % 3 != 0);
    let mut tris: Vec<[u32; 3]> = Vec::with_capacity(whole / 3);
    let mut referenced = alloc::vec![false; positions.len()];
    for face in indices[..whole].chunks_exact(3) {
        let in_range = face.iter().all(|&i| (i as usize) < positions.len());
        let distinct = face[0] != face[1] && face[1] != face[2] && face[0] != face[2];
        if !in_range || !distinct {
            report.faces_skipped += 1;
            continue;
        }
        let a = positions[face[0] as usize];
        let b = positions[face[1] as usize];
        let c = positions[face[2] as usize];
        if !finite(a) || !finite(b) || !finite(c) {
            report.non_finite_positions += 1;
            continue;
        }
        let cross = vec3::cross(vec3::sub(b, a), vec3::sub(c, a));
        if vec3::length_squared(cross) <= limit_sq {
            report.degenerate_triangles += 1;
            continue;
        }
        for &i in face {
            referenced[i as usize] = true;
        }
        tris.push([face[0], face[1], face[2]]);
    }
    report.triangles = tris.len() as u64;

    // ── mesh → field ────────────────────────────────────────────────────────
    let mut forward = Accumulator::<R>::new();
    for (i, &p) in positions.iter().enumerate() {
        if !referenced[i] {
            continue;
        }
        report.vertex_samples += 1;
        match project(
            sdf,
            p,
            R::INFINITY,
            tol,
            AccuracyConfig::MAX_NEWTON_ITERATIONS,
        ) {
            Projected::Converged(q) => forward.push(vec3::length(vec3::sub(p, q))),
            _ => report.unconverged_mesh_samples += 1,
        }
    }
    let third = R::ONE / R::from_f64(3.0);
    for t in &tris {
        let a = positions[t[0] as usize];
        let b = positions[t[1] as usize];
        let c = positions[t[2] as usize];
        let centroid = vec3::scale(
            [a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2]],
            third,
        );
        report.centroid_samples += 1;
        match project(
            sdf,
            centroid,
            R::INFINITY,
            tol,
            AccuracyConfig::MAX_NEWTON_ITERATIONS,
        ) {
            Projected::Converged(q) => forward.push(vec3::length(vec3::sub(centroid, q))),
            _ => report.unconverged_mesh_samples += 1,
        }
    }
    report.mesh_to_field = forward.finish();

    // ── field → mesh ────────────────────────────────────────────────────────
    let grid = TriangleGrid::build(positions, &tris, cfg.cell_size())?;

    // Accept only projections landing one cell inside the lattice box. For an
    // open field the mesh is clipped at the wall while the surface is not, so
    // wall-adjacent seeds would measure the domain rather than the extractor.
    // For a closed field the surface is strictly interior and this costs
    // nothing. A uniform rule, so no branch on the field.
    let lo_accept = [origin[0] + h, origin[1] + h, origin[2] + h];
    let hi_accept = [
        origin[0] + h * R::from_f64(f64::from(size[0] - 2)),
        origin[1] + h * R::from_f64(f64::from(size[1] - 2)),
        origin[2] + h * R::from_f64(f64::from(size[2] - 2)),
    ];

    let mut reverse = Accumulator::<R>::new();
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                report.seeds += 1;
                let seed = [
                    origin[0] + h * R::from_f64(f64::from(x)),
                    origin[1] + h * R::from_f64(f64::from(y)),
                    origin[2] + h * R::from_f64(f64::from(z)),
                ];
                let q = match project(sdf, seed, band, tol, AccuracyConfig::MAX_NEWTON_ITERATIONS) {
                    Projected::Converged(q) => q,
                    Projected::OutOfBand => {
                        report.seeds_out_of_band += 1;
                        continue;
                    }
                    Projected::NonFinite => {
                        report.seeds_non_finite += 1;
                        continue;
                    }
                    Projected::Unconverged => {
                        report.seeds_unconverged += 1;
                        continue;
                    }
                };
                let inside = (0..3).all(|a| q[a] >= lo_accept[a] && q[a] <= hi_accept[a]);
                if !inside {
                    report.seeds_outside_domain += 1;
                    continue;
                }
                let d2 = grid.nearest_distance_squared(q, positions, &tris);
                if d2.is_finite() {
                    reverse.push(d2.sqrt());
                } else {
                    // No triangle anywhere: an empty mesh. Counted as outside
                    // rather than silently dropped, so the buckets still sum.
                    report.seeds_outside_domain += 1;
                }
            }
        }
    }
    report.field_to_mesh = reverse.finish();

    debug_assert_eq!(
        report.seeds,
        report.seeds_out_of_band
            + report.seeds_non_finite
            + report.seeds_unconverged
            + report.seeds_outside_domain
            + report.field_to_mesh.samples,
        "a seed was dropped without landing in a bucket"
    );

    Ok(report)
}

/// Every coordinate finite.
#[inline]
fn finite<R: Real>(p: [R; 3]) -> bool {
    p[0].is_finite() && p[1].is_finite() && p[2].is_finite()
}

/// The outcome of projecting one point onto the zero set.
enum Projected<R: Real> {
    Converged([R; 3]),
    /// The first step was longer than the band radius.
    OutOfBand,
    /// `∇f` vanished or went non-finite.
    NonFinite,
    /// Hit the iteration cap.
    Unconverged,
}

/// Newton iteration onto the zero set, along the gradient.
///
/// `max_first_step` bounds the *first* step only, which is what makes the band
/// test free: for a true distance field that step length is the distance to the
/// surface, so the test costs no extra field evaluation.
fn project<R, S>(
    sdf: &S,
    start: [R; 3],
    max_first_step: R,
    tol: R,
    max_iterations: u32,
) -> Projected<R>
where
    R: Real,
    S: Sdf<Scalar = R>,
{
    let mut p = start;
    for i in 0..max_iterations {
        let f = sdf.sample(p);
        let g = sdf.gradient(p);
        let gg = vec3::dot(g, g);
        // `is_finite` carries the NaN case, so this needs no negated comparison:
        // a NaN `gg` fails the finiteness test rather than the ordering one.
        if gg <= R::ZERO || !gg.is_finite() || !f.is_finite() {
            return Projected::NonFinite;
        }
        let step = vec3::scale(g, f / gg);
        let len = vec3::length(step);
        if !len.is_finite() {
            return Projected::NonFinite;
        }
        if i == 0 && len > max_first_step {
            return Projected::OutOfBand;
        }
        p = vec3::sub(p, step);
        if !finite(p) {
            return Projected::NonFinite;
        }
        if len <= tol {
            return Projected::Converged(p);
        }
    }
    Projected::Unconverged
}

/// Max and compensated mean over a stream of non-negative distances.
///
/// [`Real`] has no widening to `f64`, so the sum happens at the mesh's own
/// precision. Naive `f32` accumulation over `~3·10⁴` positive terms carries a
/// worst-case relative error around `n·ε ≈ 4e-3`, which is the same order as the
/// quantity being measured. Neumaier compensation costs three flops per sample
/// and removes the objection; it is a fixed arithmetic sequence, so it stays
/// deterministic.
struct Accumulator<R: Real> {
    count: u64,
    max: R,
    sum: R,
    compensation: R,
}

impl<R: Real> Accumulator<R> {
    fn new() -> Self {
        Self {
            count: 0,
            max: R::ZERO,
            sum: R::ZERO,
            compensation: R::ZERO,
        }
    }

    fn push(&mut self, d: R) {
        self.count += 1;
        if d > self.max {
            self.max = d;
        }
        let t = self.sum + d;
        if self.sum.abs() >= d.abs() {
            self.compensation += (self.sum - t) + d;
        } else {
            self.compensation += (d - t) + self.sum;
        }
        self.sum = t;
    }

    fn finish(self) -> DistanceStats<R> {
        if self.count == 0 {
            return DistanceStats::EMPTY;
        }
        let total = self.sum + self.compensation;
        DistanceStats {
            samples: self.count,
            max: self.max,
            mean: total / R::from_f64(self.count as f64),
        }
    }
}

#[cfg(test)]
mod tests;
