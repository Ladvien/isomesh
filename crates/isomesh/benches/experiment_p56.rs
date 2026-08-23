//! **P-56 — the seam cone: a central difference across a `min`/`max` is wrong by
//! at most half the angle the two branches span.**
//!
//! Ticket: R-051. Pre-registered before this file existed.
//!
//! ```bash
//! cargo bench --bench experiment_p56
//! ```
//!
//! Writes `docs/experiments/p-56.csv`, one row per (dihedral, resolution).
//!
//! # What P-47 left behind
//!
//! P-47's accuracy clause died by three orders of magnitude, and its own
//! artefact said why: `bulk_mean_angular_error_deg` is `1.9e-8` while **one**
//! vertex in 57,470 carries `4.365` degrees, and `worst_stencil_straddles_seam`
//! is true from 32 brushes upward. The surviving claim is mechanical rather than
//! statistical, so it needs a fixture where the seam angle is a **dial** and not
//! whatever a random brush pile happened to produce.
//!
//! # The fixture, derived
//!
//! Two spheres, `A` of radius `r₁` centred at the origin and `B` of radius `r₂`
//! centred at `(dist, 0, 0)`, composed as [`Difference`] — `f = max(f_A, −f_B)`,
//! the crate's own `A − B`.
//!
//! Let `p` be a point of the intersection circle, so `|p| = r₁` and
//! `|p − c| = r₂`. The two surface sheets meeting there are the sphere-`A` sheet,
//! whose field gradient is the outward normal `n₁ = p / r₁`, and the sphere-`B`
//! cavity wall, whose field gradient is `−n₂` for `n₂ = (p − c) / r₂`. The
//! **dihedral** `theta` of the solid at that edge is the angle between the two
//! tangent planes measured through the material, and for two planes with outward
//! normals `u` and `v` that is `180° − angle(u, v)`. So
//!
//! ```text
//! cos theta = −(n₁ · (−n₂)) = n₁ · n₂
//! ```
//!
//! and `n₁ · n₂` is the cosine of the angle at `p` in the triangle
//! `(0, c, p)`, whose sides are `r₁`, `r₂` and `dist`:
//!
//! ```text
//! n₁ · n₂ = (r₁² + r₂² − dist²) / (2 · r₁ · r₂)
//! ```
//!
//! Setting the two equal and solving gives the law of cosines itself:
//!
//! ```text
//! dist = sqrt(r₁² + r₂² − 2 · r₁ · r₂ · cos theta)
//! ```
//!
//! — `theta = 180°` puts the spheres externally tangent (a flat, degenerate
//! seam), `theta = 90°` makes them orthogonal, `theta → 0` makes them internally
//! tangent. The seam circle then sits at `x = (r₁² + dist² − r₂²) / (2 · dist)`
//! with radius `r₁ · r₂ · sin theta / dist`.
//!
//! **That derivation is checked before anything is measured.** `seam_residual`
//! is `|angle(n₁, −n₂) − (180 − theta)|` at eight points around each circle; the
//! run aborts above `1e-9` degrees, because a wrong fixture voids every number
//! downstream.
//!
//! # Why the bound is `(180 − theta) / 2`, and why it is a theorem
//!
//! Write `s(q) = f_A(q) + f_B(q)`. Then `max(f_A, −f_B) = f_A − min(s, 0)`, so
//! branch `A` is active exactly where `s ≥ 0` — and `s = 0` is the **prolate
//! spheroid** `|q| + |q − c| = r₁ + r₂` with foci at the two centres, which
//! contains the seam circle. That surface, not the seam curve, is what a stencil
//! straddles.
//!
//! Per axis the central difference of `f` is
//!
//! ```text
//! g_i = ∂_i f_A − λ_i · ∂_i s,   λ_i = (fraction of [p − h·e_i, p + h·e_i] with s < 0)
//! ```
//!
//! so with `u = n₁`, `v = −n₂` and `w = u − v`, the returned direction is
//! `g = u − Λw` for a **diagonal** `Λ` with entries in `[0, 1]`. Over the
//! stencil `s` is affine to `O(h²)`, so `s(p) ≥ 0` forces every `λ_i ≤ ½` and
//! `s(p) < 0` forces every `λ_i ≥ ½`. In the first case
//! `|Λw|² = Σ λ_i² w_i² ≤ |w|² / 4`, and `|w| = 2·sin(angle(u, v) / 2)`, so the
//! deviation from `u` has norm at most `sin(alpha/2)` for
//! `alpha = angle(u, v) = 180 − theta`. The largest angle between a unit vector
//! and `u − z` over `|z| ≤ R < 1` is `asin(R)`, hence
//!
//! ```text
//! angle(g, active branch gradient) ≤ asin(sin(alpha / 2)) = (180 − theta) / 2
//! ```
//!
//! and symmetrically when branch `B` is active. Equality is reached at `s(p) = 0`
//! with all three `λ_i = ½`, where `g` is exactly `(u + v) / 2`. So this is not a
//! fitted envelope: it is tight, and `worst_over_bound_ratio` measures how close
//! the sampled vertices got to the corner case rather than how generous the
//! envelope is.
//!
//! # The window, and why the canonical domain cannot see this at `f64`
//!
//! `Sdf::gradient`'s step is `DIFF_STEP · max(|pₓ|, |p_y|, |p_z|, 1)`, which at
//! `f64` is `6.06e-6` — the registration says `DIFF_STEP · |p|`, which is the
//! max-norm and not the Euclidean one, and on this fixture every seam point lies
//! on the unit sphere so the factor is exactly `1` either way. A vertex straddles
//! only within about that distance of the spheroid, while dual contouring places
//! a seam vertex on the intersection of two *tangent planes* fitted from Hermite
//! samples up to half a cell away, which misses the true crease by roughly
//! `(cell/2)² / (2r)`. Over the fields' canonical `[-2, 2]³` domain that is
//! `1e-3` at 129³ against a `6e-6` stencil, so the straddle is a coincidence and
//! P-47 duly saw exactly one of them.
//!
//! This experiment therefore meshes a **cube window of side [`WINDOW`] centred on
//! a seam point** instead of the whole solid, so the cell is `1.25e-3` at the
//! coarsest resolution here rather than `0.06`. The window is a magnifying glass
//! on the seam, not a change to the mechanism: the field, the composition and the
//! differencing step are the crate's own.
//!
//! [`WINDOW`] was chosen from `(cell/2)² / (2r)`, which predicts a tangent-plane
//! miss of `2.8e-7` at 33³ — twelve stencils inside the step, i.e. every seam
//! cell already saturated at every resolution. **That prediction is wrong, and
//! the run says so**: `seam_offset_stencils` reads `75` at 30°/33³ and `48` at
//! 175°/33³. The model left out the conditioning of the QEF itself — two planes
//! meeting at `theta` locate their intersection with an error amplified by
//! `1 / sin theta`, which is `2` at 30° and `11.5` at 175°, and it is worst at
//! exactly the two ends of the sweep. Near 120° the model is right and the
//! diagnostic reads `0.74` down to `0.08`. The estimate is left in place rather
//! than retro-fitted: it is why the window is `0.04` and not something else, and
//! the column that refutes it is in the CSV.
//!
//! # `seam_offset_stencils`, and what it does to C2's fit
//!
//! C2 is registered as a **fitted exponent** of `straddling_vertices` against
//! `n`, required under `1.5`. That statistic turns out not to be a property of
//! the mechanism: it depends on whether the grid has resolved the crease to
//! within one differencing step, which is a joint function of cell size, `theta`
//! and `DIFF_STEP` and has nothing to do with the seam's dimensionality. Below
//! saturation the count is still climbing toward its ceiling and the fit reads
//! the climb; at saturation it reads `1`.
//!
//! So the ceiling is recorded directly beside it, and it is the registered
//! sentence — *"the count of such vertices scales like the seam's length in
//! cells"* — in a form no crossover can move: `straddling_per_seam_cell`, the
//! straddling count over the number of grid cells the seam circle passes
//! through, and `straddling_share_of_vertices`. The first cannot exceed about
//! `1` because a straddling vertex is within one stencil of the spheroid and the
//! spheroid meets the surface only at the seam; the second must fall like `1/n`
//! if the population is a curve on an `n²` surface. Both are ratios of counts.
//! This is post-hoc, in M-288's sense: the registered fit stays and is reported
//! as `scaling_exponent`, and these are measurements offered beside it.
//!
//! # Dual contouring, not surface nets
//!
//! Surface nets places a cell's vertex at the centroid of its edge crossings,
//! which in a seam cell is a point in the cell's interior with no relation to the
//! crease — it would report zero straddling vertices at any resolution and
//! measure nothing. Dual contouring solves the QEF, and on a crease that solution
//! *is* the crease. The straddling population exists because dual contouring is
//! good at sharp features, which is the honest statement of who this finding is
//! about.
//!
//! # Angles are measured with `atan2`, not `acos`
//!
//! `acos` near `1` has resolution `sqrt(2·eps) = 2.1e-8` radians, i.e. `1.2e-6`
//! degrees — larger than C3's whole threshold, so an `acos` control column would
//! be reporting its own conditioning. [`angle_deg`] uses
//! `2·atan2(|û − v̂|, |û + v̂|)`, which is accurate at both ends.

mod common;

use core::f64::consts::TAU;

use isomesh::dual_contouring::DualContouring;
use isomesh::fields::{Difference, Sphere};
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

/// The scalar the whole experiment runs in.
///
/// `f64`, because C3 asks for a control under `1e-6` degrees and an `f32`
/// central difference is wrong by `1.4e-3` degrees on a *smooth* patch. The
/// scalar is not a free parameter here; C3 picks it.
type Scalar = f64;

/// `A − B` for two spheres: the crate's own `max(f_A, −f_B)`.
type Fixture = Difference<Sphere<Scalar>, Sphere<Scalar>>;

/// Radius of the solid being cut.
const R1: Scalar = 1.0;

/// Radius of the tool. Deliberately different from [`R1`]: equal radii make the
/// seam symmetric about the mid-plane and would let a sign error cancel.
const R2: Scalar = 0.7;

/// Side of the cube window centred on a seam point.
///
/// Sized from the tangent-plane estimate in the module header, which turns out
/// to hold only near 120°: the crease is resolved to well inside a differencing
/// step there at every resolution, and to `3.3` steps at 30°/257³. Also smaller
/// than the smallest seam radius in the sweep (`0.036` at 175°), so the window
/// never swallows the far side of the circle.
const WINDOW: Scalar = 0.04;

/// The swept dihedrals, in degrees, as registered: 30 through 175.
const DIHEDRALS: [Scalar; 7] = [30.0, 60.0, 90.0, 120.0, 150.0, 165.0, 175.0];

/// Samples per axis. Four points for C2's fit, doubling each time.
const RESOLUTIONS: [u32; 4] = [33, 65, 129, 257];

/// Points used to walk the seam circle when counting the cells it enters.
///
/// The finest cell in the sweep is `1.56e-4` and the longest circle is `4.2`
/// across, so this is a step of `4.2e-6` — at least thirty samples per cell
/// everywhere, which is what makes the count a count and not a sampling artefact.
const CIRCLE_SAMPLES: usize = 1_000_000;

/// How far the measured seam angle may sit from the constructed one, in degrees.
const SEAM_TOLERANCE_DEG: Scalar = 1e-9;

/// C2's ceiling on the fitted exponent.
const EXPONENT_CEILING: Scalar = 1.5;

/// Below this median tightness the registration calls C1 a vacuous pass.
const VACUITY_THRESHOLD: Scalar = 0.1;

fn dot(u: [Scalar; 3], v: [Scalar; 3]) -> Scalar {
    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
}

fn norm(u: [Scalar; 3]) -> Scalar {
    dot(u, u).sqrt()
}

/// `u`, scaled to unit length. `None` when it has no length to scale.
fn unit(u: [Scalar; 3]) -> Option<[Scalar; 3]> {
    let n = norm(u);
    if n > 0.0 && n.is_finite() {
        let inv = n.recip();
        Some([u[0] * inv, u[1] * inv, u[2] * inv])
    } else {
        None
    }
}

/// The angle between two directions, in degrees, accurate at both ends.
///
/// `2·atan2(|û − v̂|, |û + v̂|)`. See the module header for why not `acos`.
fn angle_deg(u: [Scalar; 3], v: [Scalar; 3]) -> Option<Scalar> {
    let (a, b) = (unit(u)?, unit(v)?);
    let diff = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let sum = [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
    Some((2.0 * norm(diff).atan2(norm(sum))).to_degrees())
}

/// One fixture: the two spheres, and where their seam is.
struct Seam {
    dihedral_deg: Scalar,
    /// Distance between the two centres, from the law of cosines.
    separation: Scalar,
    /// `x` of the plane the seam circle lies in.
    plane_x: Scalar,
    /// Radius of the seam circle.
    radius: Scalar,
    field: Fixture,
}

impl Seam {
    /// The fixture whose seam dihedral is `dihedral_deg`.
    fn new(dihedral_deg: Scalar) -> Self {
        let cos_theta = dihedral_deg.to_radians().cos();
        let separation = (R1 * R1 + R2 * R2 - 2.0 * R1 * R2 * cos_theta).sqrt();
        let plane_x = (R1 * R1 + separation * separation - R2 * R2) / (2.0 * separation);
        // From `r₁² = plane_x² + radius²`. Written this way rather than as
        // `r₁·r₂·sin(theta)/separation` so the circle is exactly on sphere A to
        // the last bit, which is what the residual check is looking at.
        let radius = (R1 * R1 - plane_x * plane_x).max(0.0).sqrt();
        Self {
            dihedral_deg,
            separation,
            plane_x,
            radius,
            field: Difference {
                a: Sphere {
                    center: [0.0; 3],
                    radius: R1,
                },
                b: Sphere {
                    center: [separation, 0.0, 0.0],
                    radius: R2,
                },
            },
        }
    }

    /// The seam point at parameter `t`.
    fn point(&self, t: Scalar) -> [Scalar; 3] {
        [self.plane_x, self.radius * t.cos(), self.radius * t.sin()]
    }

    /// Worst `|angle(n₁, −n₂) − (180 − theta)|` over eight points of the circle.
    ///
    /// The construction check. A fixture that fails this has the wrong
    /// separation, and every number measured on it is about a different angle
    /// than the one in its own row.
    fn residual_deg(&self) -> Scalar {
        let target = 180.0 - self.dihedral_deg;
        let mut worst: Scalar = 0.0;
        for k in 0..8 {
            let p = self.point(TAU * f64::from(k) / 8.0);
            let n1 = self.field.a.gradient(p);
            let g2 = self.field.b.gradient(p);
            let Some(measured) = angle_deg(n1, [-g2[0], -g2[1], -g2[2]]) else {
                return Scalar::INFINITY;
            };
            worst = worst.max((measured - target).abs());
        }
        worst
    }

    /// Origin of the cube window, centred on the seam point at `t = 0`.
    fn window_origin(&self) -> [Scalar; 3] {
        let c = self.point(0.0);
        [
            c[0] - WINDOW * 0.5,
            c[1] - WINDOW * 0.5,
            c[2] - WINDOW * 0.5,
        ]
    }

    /// Cells of the window grid that the seam circle enters.
    fn seam_cells(&self, origin: [Scalar; 3], cell: Scalar, cells: u32) -> usize {
        let extent = Scalar::from(cells);
        let mut hit: Vec<[i32; 3]> = Vec::new();
        for k in 0..CIRCLE_SAMPLES {
            let p = self.point(TAU * k as Scalar / CIRCLE_SAMPLES as Scalar);
            let mut index = [0i32; 3];
            let mut inside = true;
            for axis in 0..3 {
                let f = (p[axis] - origin[axis]) / cell;
                if f < 0.0 || f >= extent {
                    inside = false;
                    break;
                }
                index[axis] = f as i32;
            }
            if inside {
                hit.push(index);
            }
        }
        hit.sort_unstable();
        hit.dedup();
        hit.len()
    }
}

/// The differencing step `Sdf::gradient` uses at `p`.
fn step(p: [Scalar; 3]) -> Scalar {
    <Scalar as Real>::DIFF_STEP * p[0].abs().max(p[1].abs()).max(p[2].abs()).max(1.0)
}

/// The composed value at `q`, and whether branch `A` was the active one.
///
/// The value is derived from the two branch samples rather than read from
/// `Difference::sample`, and is bit-identical to it: that function *is*
/// `if fa >= -fb { fa } else { -fb }`. Computing both here is what makes the
/// straddle test exact instead of a Lipschitz bound on the margin.
fn value(field: &Fixture, q: [Scalar; 3]) -> (Scalar, bool) {
    let fa = field.a.sample(q);
    let fb = field.b.sample(q);
    if fa >= -fb { (fa, true) } else { (-fb, false) }
}

/// What the six-sample stencil at `p` returned, and whether it crossed branches.
struct Probe {
    gradient: [Scalar; 3],
    straddles: bool,
}

/// The crate's default central difference, re-implemented here.
///
/// It has to be re-implemented: [`Difference`] **overrides** `Sdf::gradient` with
/// the analytic active-branch gradient, so calling `field.gradient` would return
/// the thing this experiment is measuring *against*. The formula is `sdf.rs`'s,
/// step for step.
fn probe(field: &Fixture, p: [Scalar; 3]) -> Probe {
    let h = step(p);
    let inv = (2.0 * h).recip();
    let mut gradient = [0.0; 3];
    let mut branch = [false; 6];
    for axis in 0..3 {
        let mut lo = p;
        let mut hi = p;
        lo[axis] -= h;
        hi[axis] += h;
        let (v_lo, a_lo) = value(field, lo);
        let (v_hi, a_hi) = value(field, hi);
        gradient[axis] = (v_hi - v_lo) * inv;
        branch[2 * axis] = a_lo;
        branch[2 * axis + 1] = a_hi;
    }
    Probe {
        gradient,
        straddles: branch.iter().any(|b| *b != branch[0]),
    }
}

/// Least-squares slope of `ln y` against `ln x`.
///
/// `None` when fewer than two points survive, which is the honest answer for a
/// dihedral whose straddling population never left zero.
fn fit_exponent(points: &[(Scalar, Scalar)]) -> Option<Scalar> {
    if points.len() < 2 {
        return None;
    }
    let n = points.len() as Scalar;
    let mean_x = points.iter().map(|p| p.0.ln()).sum::<Scalar>() / n;
    let mean_y = points.iter().map(|p| p.1.ln()).sum::<Scalar>() / n;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    for &(x, y) in points {
        let dx = x.ln() - mean_x;
        sxy += dx * (y.ln() - mean_y);
        sxx += dx * dx;
    }
    if sxx > 0.0 { Some(sxy / sxx) } else { None }
}

/// Median of a sample, by `total_cmp` so a NaN surfaces rather than vanishing.
fn median(mut v: Vec<Scalar>) -> Option<Scalar> {
    if v.is_empty() {
        return None;
    }
    v.sort_unstable_by(Scalar::total_cmp);
    let mid = v.len() / 2;
    if v.len().is_multiple_of(2) {
        Some((v[mid - 1] + v[mid]) * 0.5)
    } else {
        Some(v[mid])
    }
}

/// One (dihedral, resolution) measurement.
struct Outcome {
    dihedral_deg: Scalar,
    samples: u32,
    cell_size: Scalar,
    seam_cells: usize,
    vertices: usize,
    straddling: usize,
    straddling_max_deg: Option<Scalar>,
    straddling_mean_deg: Option<Scalar>,
    non_straddling: usize,
    non_straddling_mean_deg: Option<Scalar>,
    degenerate: usize,
    /// Median distance from the branch surface, in stencil widths, over the
    /// `seam_cells` vertices closest to it. Under `1` the stencil reaches the
    /// branch surface and the vertex straddles; well over `1` the grid has not
    /// yet resolved the crease to within a differencing step. This is the
    /// crossover diagnostic C2's fit turns out to need.
    seam_offset_stencils: Option<Scalar>,
}

impl Outcome {
    fn bound_deg(&self) -> Scalar {
        (180.0 - self.dihedral_deg) * 0.5
    }

    fn ratio(&self) -> Option<Scalar> {
        self.straddling_max_deg.map(|e| e / self.bound_deg())
    }

    /// The registered claim in its direct form: is the straddling population
    /// the seam's length in cells, rather than the surface's area in cells?
    fn per_seam_cell(&self) -> Option<Scalar> {
        (self.seam_cells > 0).then(|| self.straddling as Scalar / self.seam_cells as Scalar)
    }

    fn share_of_vertices(&self) -> Option<Scalar> {
        (self.vertices > 0).then(|| self.straddling as Scalar / self.vertices as Scalar)
    }
}

/// How far `p` sits from the branch surface `s = 0`, measured in stencil widths.
///
/// `s(p) / (h · max_i |∂_i s|)`, absolute. The straddle test flips an axis when
/// this drops below `1`, since `s` is affine to `O(h²)` over a stencil this
/// small, so it is the same quantity the classification uses — reported as a
/// continuum instead of a bit.
fn seam_offset_stencils(field: &Fixture, p: [Scalar; 3]) -> Option<Scalar> {
    let n1 = field.a.gradient(p);
    let n2 = field.b.gradient(p);
    let grad_s = [n1[0] + n2[0], n1[1] + n2[1], n1[2] + n2[2]];
    let widest = grad_s[0].abs().max(grad_s[1].abs()).max(grad_s[2].abs());
    if widest > 0.0 {
        let s = field.a.sample(p) + field.b.sample(p);
        Some(s.abs() / (step(p) * widest))
    } else {
        None
    }
}

/// Mesh the window and classify every vertex.
fn measure(seam: &Seam, samples: u32) -> Outcome {
    let origin = seam.window_origin();
    let cell_size = WINDOW / Scalar::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("window grid fits u32");
    let mut out = MeshBuffer::<Scalar>::new();
    DualContouring::<Scalar>::new()
        .extract(&seam.field, &shape, origin, cell_size, &mut out)
        .expect("dual contouring on a two-sphere difference");
    let seam_cells = seam.seam_cells(origin, cell_size, samples - 1);

    let mut straddling = Vec::new();
    let mut smooth = Vec::new();
    let mut offsets = Vec::new();
    let mut degenerate = 0usize;
    for &p in &out.positions {
        let probed = probe(&seam.field, p);
        let analytic = seam.field.gradient(p);
        let Some(error) = angle_deg(probed.gradient, analytic) else {
            degenerate += 1;
            continue;
        };
        if let Some(offset) = seam_offset_stencils(&seam.field, p) {
            offsets.push(offset);
        }
        if probed.straddles {
            straddling.push(error);
        } else {
            smooth.push(error);
        }
    }

    let mean = |v: &[Scalar]| -> Option<Scalar> {
        if v.is_empty() {
            None
        } else {
            Some(v.iter().sum::<Scalar>() / v.len() as Scalar)
        }
    };
    let worst = straddling.iter().copied().max_by(Scalar::total_cmp);

    // The `seam_cells` vertices nearest the branch surface are the seam
    // population; anything past that is the smooth sheet and would drown the
    // statistic in numbers about the wrong vertices.
    offsets.sort_unstable_by(Scalar::total_cmp);
    offsets.truncate(seam_cells);

    Outcome {
        dihedral_deg: seam.dihedral_deg,
        samples,
        cell_size,
        seam_cells,
        vertices: out.positions.len(),
        straddling: straddling.len(),
        straddling_max_deg: worst,
        straddling_mean_deg: mean(&straddling),
        non_straddling: smooth.len(),
        non_straddling_mean_deg: mean(&smooth),
        degenerate,
        seam_offset_stencils: median(offsets),
    }
}

fn show(v: Option<Scalar>) -> String {
    v.map_or_else(|| String::from("n/a"), |x| format!("{x:.6e}"))
}

fn main() {
    let prereg = isomesh::experiment!("P-56");

    common::experiment::run(prereg, |run| {
        // ── the construction check, first, because it can void everything ──
        println!("fixture check — angle(n1, -n2) at the seam against 180 - theta");
        let seams: Vec<Seam> = DIHEDRALS.iter().copied().map(Seam::new).collect();
        let mut worst_residual: Scalar = 0.0;
        for seam in &seams {
            let residual = seam.residual_deg();
            worst_residual = worst_residual.max(residual);
            println!(
                "  theta {:>6.1}°  dist {:.9}  seam x {:.9} r {:.9}  residual {residual:.3e}°",
                seam.dihedral_deg, seam.separation, seam.plane_x, seam.radius
            );
        }
        assert!(
            worst_residual < SEAM_TOLERANCE_DEG,
            "fixture is wrong: worst seam-angle residual {worst_residual:.3e}° \
             exceeds {SEAM_TOLERANCE_DEG:e}°, so no downstream number is about \
             the dihedral its row claims"
        );
        println!("  worst residual {worst_residual:.3e}° — construction verified\n");

        // ── the sweep ──
        let mut outcomes: Vec<Outcome> = Vec::new();
        let mut exponents: Vec<(Scalar, Option<Scalar>)> = Vec::new();
        for seam in &seams {
            let mut here: Vec<Outcome> = Vec::new();
            for samples in RESOLUTIONS {
                let o = measure(seam, samples);
                println!(
                    "theta {:>6.1}°  {samples:>4}³  cell {:.3e}  seam cells {:>5}  \
                     vertices {:>7}  straddling {:>5}  per seam cell {:>10}  \
                     seam offset {:>10} stencils  worst {:>12}  bound {:>7.4}°  \
                     ratio {:>10}  control mean {:>12}",
                    o.dihedral_deg,
                    o.cell_size,
                    o.seam_cells,
                    o.vertices,
                    o.straddling,
                    show(o.per_seam_cell()),
                    show(o.seam_offset_stencils),
                    show(o.straddling_max_deg),
                    o.bound_deg(),
                    show(o.ratio()),
                    show(o.non_straddling_mean_deg),
                );
                here.push(o);
            }
            let fit: Vec<(Scalar, Scalar)> = here
                .iter()
                .filter(|o| o.straddling > 0)
                .map(|o| (Scalar::from(o.samples), o.straddling as Scalar))
                .collect();
            let exponent = fit_exponent(&fit);
            println!(
                "  → straddling scaling exponent against n: {}\n",
                exponent.map_or_else(|| String::from("n/a"), |e| format!("{e:.4}"))
            );
            exponents.push((seam.dihedral_deg, exponent));
            outcomes.extend(here);
        }

        // ── verdicts ──
        let ratios: Vec<Scalar> = outcomes.iter().filter_map(Outcome::ratio).collect();
        let sweep_median = median(ratios.clone());
        let min_ratio = ratios.iter().copied().min_by(Scalar::total_cmp);
        let max_ratio = ratios.iter().copied().max_by(Scalar::total_cmp);
        let measured_rows = outcomes.iter().filter(|o| o.straddling > 0).count();
        let breaches = outcomes
            .iter()
            .filter(|o| o.straddling_max_deg.is_some_and(|e| e > o.bound_deg()))
            .count();
        let worst_exponent = exponents
            .iter()
            .filter_map(|&(_, e)| e)
            .max_by(Scalar::total_cmp);
        let worst_control = outcomes
            .iter()
            .filter_map(|o| o.non_straddling_mean_deg)
            .max_by(Scalar::total_cmp);

        println!(
            "C1  rows with a straddling population: {measured_rows}/{}",
            outcomes.len()
        );
        println!("C1  bound breaches: {breaches}");
        println!(
            "C1  tightness  min {}  median {}  max {}",
            show(min_ratio),
            show(sweep_median),
            show(max_ratio)
        );
        println!(
            "C1  {}",
            match sweep_median {
                Some(m) if m < VACUITY_THRESHOLD =>
                    format!("VACUOUS PASS — median tightness {m:.4} < {VACUITY_THRESHOLD}"),
                Some(m) => format!("tightness median {m:.4} ≥ {VACUITY_THRESHOLD}, not vacuous"),
                None => String::from("UNMEASURED — no straddling vertex anywhere"),
            }
        );
        let worst_per_seam_cell = outcomes
            .iter()
            .filter_map(Outcome::per_seam_cell)
            .max_by(Scalar::total_cmp);
        let shares: Vec<Option<Scalar>> = RESOLUTIONS
            .iter()
            .map(|&n| {
                outcomes
                    .iter()
                    .filter(|o| o.samples == n)
                    .filter_map(Outcome::share_of_vertices)
                    .max_by(Scalar::total_cmp)
            })
            .collect();
        println!(
            "C2  worst fitted exponent {} (ceiling {EXPONENT_CEILING})",
            show(worst_exponent)
        );
        println!(
            "C2  per-dihedral exponents: {}",
            exponents
                .iter()
                .map(|&(t, e)| format!(
                    "{t:.0}°:{}",
                    e.map_or_else(|| String::from("n/a"), |x| format!("{x:.3}"))
                ))
                .collect::<Vec<_>>()
                .join("  ")
        );
        println!(
            "C2  direct form — worst straddling / seam_cells over all rows: {}",
            show(worst_per_seam_cell)
        );
        println!(
            "C2  worst straddling / vertices by resolution: {}",
            RESOLUTIONS
                .iter()
                .zip(&shares)
                .map(|(n, s)| format!("{n}³:{}", show(*s)))
                .collect::<Vec<_>>()
                .join("  ")
        );
        println!(
            "C3  worst non-straddling mean error {}°",
            show(worst_control)
        );

        let sweep_median_text = show(sweep_median);
        let exponent_of = |theta: Scalar| -> String {
            exponents
                .iter()
                .find(|&&(t, _)| (t - theta).abs() < 1e-12)
                .and_then(|&(_, e)| e)
                .map_or_else(|| String::from("n/a"), |e| format!("{e:.6}"))
        };
        for o in &outcomes {
            run.record(&[
                ("dihedral_deg", format!("{:.1}", o.dihedral_deg)),
                ("samples_per_axis", o.samples.to_string()),
                ("seam_cells", o.seam_cells.to_string()),
                ("vertices", o.vertices.to_string()),
                ("straddling_vertices", o.straddling.to_string()),
                ("straddling_max_error_deg", show(o.straddling_max_deg)),
                ("predicted_bound_deg", format!("{:.6}", o.bound_deg())),
                ("worst_over_bound_ratio", show(o.ratio())),
                (
                    "within_bound",
                    o.straddling_max_deg
                        .map_or_else(|| String::from("n/a"), |e| (e <= o.bound_deg()).to_string()),
                ),
                (
                    "non_straddling_mean_error_deg",
                    show(o.non_straddling_mean_deg),
                ),
                ("scaling_exponent", exponent_of(o.dihedral_deg)),
                // Extras: the fixture, the window, and the population the
                // control was taken over.
                (
                    "separation",
                    format!("{:.9}", Seam::new(o.dihedral_deg).separation),
                ),
                (
                    "seam_radius",
                    format!("{:.9}", Seam::new(o.dihedral_deg).radius),
                ),
                ("window_side", format!("{WINDOW}")),
                ("cell_size", format!("{:.6e}", o.cell_size)),
                ("straddling_mean_error_deg", show(o.straddling_mean_deg)),
                ("non_straddling_vertices", o.non_straddling.to_string()),
                ("degenerate_normals", o.degenerate.to_string()),
                ("sweep_median_ratio", sweep_median_text.clone()),
                ("seam_angle_residual_deg", format!("{worst_residual:.3e}")),
                ("straddling_per_seam_cell", show(o.per_seam_cell())),
                ("straddling_share_of_vertices", show(o.share_of_vertices())),
                ("seam_offset_stencils", show(o.seam_offset_stencils)),
            ]);
        }
    });
}
