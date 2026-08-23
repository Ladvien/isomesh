//! E-311 — the seam cone: half the crease angle, and it is *attained*.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example seam_normal_bound --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 20-50x slower, and this one
//! extracts a 129³ grid and then re-differences the field at every vertex it
//! produced.
//!
//! `1` `2` `3` freeze on the three ledger rows the startup self-check
//! reproduces — 90°, 120°, 175°. `X` restarts the sweep. The rest are the
//! shared keys — `W` wireframe, `G` the window box, `Space` freeze,
//! `F12` screenshot.
//!
//! Under `ISOMESH_CAPTURE` it needs no keyboard: the dihedral advances one
//! schedule entry per captured frame, and the schedule is exactly 80 entries
//! long — the harness's own default `ISOMESH_CAPTURE_FRAMES`, so a default
//! capture is exactly one sweep from 30° to 175°. `ISOMESH_FIELD=0|1|2` pins
//! 90°, 120° or 175° for a still.
//!
//! ```bash
//! # FPS=12 rather than the script's 20: 80 frames at 20 fps is four seconds
//! # for a sweep with seven ledger stops in it, and the numbers are the demo.
//! # The camera never moves and only the wedge and the bar do, which is what a
//! # GIF's inter-frame compression wants.
//! FPS=12 ./scripts/record_gif.sh seam_normal_bound docs/gifs/e311.gif
//! ```
//!
//! Demonstrates **M-350 / P-56** (`docs/experiments/p-56.csv`).
//!
//! # The claim, and why it is a theorem rather than an envelope
//!
//! At a `min`/`max` CSG seam the field is C0 and not C1. A six-sample central
//! difference whose stencil **straddles** the seam averages two different
//! branch gradients, so the direction it returns lies inside the cone the two
//! branches span — and is therefore wrong by at most **half** the angle between
//! them. With `u` and `v` the two outward branch normals and
//! `alpha = angle(u, v) = 180 - theta` for a seam dihedral `theta`, the
//! returned direction is `u - Λ(u - v)` for a diagonal `Λ` with entries in
//! `[0, 1]`, every entry on one side of `½`. That bounds the deviation by
//! `sin(alpha/2)` in norm, hence
//!
//! ```text
//! angle(central difference, active branch gradient) <= (180 - theta) / 2
//! ```
//!
//! Equality is reached when all three entries are exactly `½`, where the
//! returned direction is exactly `(u + v) / 2`. So the interesting number is not
//! *does it hold* — it must — but **how close the sampled vertices get to that
//! corner**. That is `worst_over_bound_ratio`, and the ledger's median over its
//! 24 measured rows is `0.9748`. The bar along the bottom of this window is that
//! ratio, live.
//!
//! # What is on screen
//!
//! - **The surface** — two spheres composed by [`Difference`], meshed by Dual
//!   Contouring. Inside the window it is two sheets meeting along one crease:
//!   sphere `A`'s outer surface and the cavity wall sphere `B` cut into it. The
//!   angle of that crease *is* `theta`, and it is what sweeps.
//! - **White polyline** — the seam circle itself, the exact crease, from the
//!   fixture's closed form rather than from the mesh.
//! - **Cyan rays** — the analytic active-branch gradient at a straddling vertex,
//!   which is the right answer.
//! - **Pink rays** — what the central difference returned there, which is not.
//! - **Amber arcs** — the angle between the two, shaded. That angle is the
//!   error, and the widest one in the frame is the number on the bar.
//! - **White dot with the long pair of rays** — the *worst* straddling vertex,
//!   the one `straddling_max_error_deg` is measured at.
//! - **Grey box** — the sampled window, on `G`.
//!
//! Non-straddling vertices are drawn as nothing at all. Their mean angular error
//! at 129³ is `4.2e-10` to `6.5e-10` degrees across the sweep, which is `f64`
//! round-off and not a finding; a marker per vertex would bury the ones that
//! matter under 28,000 that do not.
//!
//! # Why the view is so tight, and why it has to be
//!
//! `Sdf::gradient`'s step is `DIFF_STEP · max(|pₓ|, |p_y|, |p_z|, 1)`, which at
//! `f64` is `6.06e-6`. A vertex straddles only within about that distance of the
//! branch surface `f_A + f_B = 0`. Dual contouring, meanwhile, places a seam
//! vertex on the intersection of two tangent planes fitted from Hermite samples
//! up to half a cell away, which misses the true crease by roughly
//! `(cell/2)² / (2r)` — `1e-3` at 129³ over the fields' canonical `[-2, 2]³`
//! domain, against a `6e-6` stencil. On the canonical domain a straddle is
//! therefore a **coincidence**: P-47 saw exactly one, in 57,470 vertices.
//!
//! So this meshes a cube window of side [`WINDOW`] centred on one seam point
//! instead of the whole solid, which takes the cell to `3.125e-4` at 129³. The
//! window is a magnifying glass on the seam, not a change to the mechanism: the
//! field, the composition and the differencing step are the crate's own, and
//! `0.04` is smaller than the smallest seam radius in the sweep (`0.036` at
//! 175°) so the window never swallows the far side of the circle.
//!
//! # The fixture is derived, not fitted
//!
//! Two spheres, `A` of radius [`R1`] at the origin and `B` of radius [`R2`] at
//! `(dist, 0, 0)`, composed as `A − B` — the crate's own `max(f_A, −f_B)`. At a
//! point `p` of the intersection circle the two sheets have outward normals
//! `n₁ = p / r₁` and `−n₂` for `n₂ = (p − c) / r₂`, and the dihedral through the
//! material is `180° − angle(n₁, −n₂)`, so `cos theta = n₁ · n₂`. That dot
//! product is the cosine of the angle at `p` in the triangle `(0, c, p)`, giving
//! the law of cosines itself:
//!
//! ```text
//! dist = sqrt(r₁² + r₂² − 2 · r₁ · r₂ · cos theta)
//! ```
//!
//! `theta = 180°` puts the spheres externally tangent (flat), `theta = 90°`
//! makes them orthogonal, `theta → 0` makes them internally tangent — a knife
//! edge. The startup check measures `|angle(n₁, −n₂) − (180 − theta)|` at eight
//! points of each of the seven ledger circles and holds it against the CSV's own
//! `seam_angle_residual_deg` of `1.243e-13`, because a wrong fixture voids every
//! number downstream.
//!
//! # Angles are measured with `atan2`, not `acos`
//!
//! `acos` near `1` resolves only to `sqrt(2·eps) = 2.1e-8` radians, i.e.
//! `1.2e-6` degrees — larger than the whole non-straddling control, so an `acos`
//! column would be reporting its own conditioning rather than the field's.
//! [`angle_deg`] uses `2·atan2(|û − v̂|, |û + v̂|)`, accurate at both ends.
//!
//! # The bar dips at 30°, and that is a measurement
//!
//! At `theta = 30` the ledger's ratio at 129³ is `0.2648`, not `0.99`. The
//! reason is in the CSV beside it: `seam_offset_stencils` reads `8.14` there,
//! meaning the nearest vertices sit eight differencing steps from the branch
//! surface, so the grid has not resolved the crease to within a stencil and only
//! ten vertices straddle at all — none of them near the corner case. Two planes
//! meeting at `theta` locate their intersection with an error amplified by
//! `1 / sin theta`, which is worst at both ends of the sweep. The `offset` line
//! on the HUD is that explanation, live, and the bar is left to dip rather than
//! rescaled to look full.
//!
//! # The view frame is the seam's own
//!
//! The picture is drawn in the frame `(e₁, e₂, ẑ)` where `e₁` bisects the two
//! outward normals and `ẑ` is the crease tangent, scaled so the window is a unit
//! cube. That is a rigid rotation about `z` plus a translation and a uniform
//! scale, so **every angle on screen is the angle in the field**; it exists so
//! the crease holds still while `theta` sweeps, and so the geometry is order-one
//! in `f32` rather than a `0.04`-wide detail at `x ≈ 0.9`.
//!
//! # `f64`, and Dual Contouring only
//!
//! M-350 asks for a control under `1e-6` degrees and an `f32` central difference
//! is wrong by `1.4e-3` degrees on a *smooth* patch, so the scalar is not a free
//! parameter here. The surface is cast to `f32` on its way into the [`Mesh`]
//! asset and nothing but the picture depends on that.
//!
//! Dual Contouring rather than Surface Nets because Surface Nets places a cell's
//! vertex at the centroid of its edge crossings, which in a seam cell is a point
//! in the cell's interior with no relation to the crease — it would report zero
//! straddling vertices at any resolution and measure nothing. Dual contouring
//! solves the QEF, and on a crease that solution *is* the crease. The straddling
//! population exists because dual contouring is good at sharp features, which is
//! the honest statement of who this finding is about.

mod common;

use core::f64::consts::TAU;
use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{
    Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags, samples_override,
};
use isomesh::dual_contouring::DualContouring;
use isomesh::fields::{Difference, Sphere as SdfSphere};
use isomesh::{MeshBuffer, MeshSink, Real, RuntimeShape3, Sdf};

// ─── the registered fixture, copied from `benches/experiment_p56.rs` ────────

/// `A − B` for two spheres: the crate's own `max(f_A, −f_B)`.
type Fixture = Difference<SdfSphere<f64>, SdfSphere<f64>>;

/// Radius of the solid being cut.
const R1: f64 = 1.0;

/// Radius of the tool. Deliberately different from [`R1`]: equal radii make the
/// seam symmetric about the mid-plane and would let a sign error cancel.
const R2: f64 = 0.7;

/// Side of the cube window centred on a seam point.
///
/// Smaller than the smallest seam radius in the sweep (`0.036` at 175°), so the
/// window never swallows the far side of the circle. See the module docs for why
/// the canonical domain cannot see this effect at `f64` at all.
const WINDOW: f64 = 0.04;

/// P-56's resolution for the row this example reproduces. `129` samples span
/// `128` cells across the window, so the cell is `3.125e-4`.
const LEDGER_SAMPLES: u32 = 129;

/// Below this the window has no interior cell for the crease to cross.
const MIN_SAMPLES: u32 = 9;

/// Above this a rebuild costs more than the sweep is worth.
const MAX_SAMPLES: u32 = 257;

/// Points used to walk the seam circle when counting the cells it enters.
///
/// The bench's number, unchanged, because the count is only a count if the walk
/// out-samples the grid: the finest cell here is `1.56e-4` and the longest
/// circle is `4.2` across, so this is a step of `4.2e-6`.
const CIRCLE_SAMPLES: usize = 1_000_000;

/// How far the measured seam angle may sit from the constructed one, in degrees.
const SEAM_TOLERANCE_DEG: f64 = 1e-9;

/// The seven dihedrals P-56 swept, in the CSV's order. The fixture check walks
/// all of them, so a wrong construction is loud even at a dihedral the sweep
/// happens not to be sitting on.
const DIHEDRALS: [f64; 7] = [30.0, 60.0, 90.0, 120.0, 150.0, 165.0, 175.0];

/// The three rows the startup self-check reproduces, in E-311's own order:
/// the flattest crease first, so the log reads from the tightest bound down to
/// the widest.
const SELF_CHECK: [f64; 3] = [175.0, 120.0, 90.0];

/// The same three rows in ascending order, which is what `1`, `2`, `3` and
/// `ISOMESH_FIELD=0|1|2` select.
///
/// Ascending rather than sharing [`SELF_CHECK`]'s order because a key row and a
/// log row are read by different people: `1` next to `2` next to `3` has to walk
/// the crease one way, and having `1` mean the *last* thing the log printed is
/// the kind of off-by-one nothing on screen would report.
const LEDGER_PINS: [f64; 3] = [90.0, 120.0, 175.0];

fn dot(u: [f64; 3], v: [f64; 3]) -> f64 {
    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
}

fn norm(u: [f64; 3]) -> f64 {
    dot(u, u).sqrt()
}

/// `u`, scaled to unit length. `None` when it has no length to scale.
fn unit(u: [f64; 3]) -> Option<[f64; 3]> {
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
/// `2·atan2(|û − v̂|, |û + v̂|)`. See the module docs for why not `acos`.
fn angle_deg(u: [f64; 3], v: [f64; 3]) -> Option<f64> {
    let (a, b) = (unit(u)?, unit(v)?);
    let diff = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let sum = [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
    Some((2.0 * norm(diff).atan2(norm(sum))).to_degrees())
}

/// One fixture: the two spheres, and where their seam is.
struct Seam {
    dihedral_deg: f64,
    /// Distance between the two centres, from the law of cosines.
    separation: f64,
    /// `x` of the plane the seam circle lies in.
    plane_x: f64,
    /// Radius of the seam circle.
    radius: f64,
    field: Fixture,
}

impl Seam {
    /// The fixture whose seam dihedral is `dihedral_deg`.
    fn new(dihedral_deg: f64) -> Self {
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
                a: SdfSphere {
                    center: [0.0; 3],
                    radius: R1,
                },
                b: SdfSphere {
                    center: [separation, 0.0, 0.0],
                    radius: R2,
                },
            },
        }
    }

    /// The seam point at parameter `t`.
    fn point(&self, t: f64) -> [f64; 3] {
        [self.plane_x, self.radius * t.cos(), self.radius * t.sin()]
    }

    /// Worst `|angle(n₁, −n₂) − (180 − theta)|` over eight points of the circle.
    fn residual_deg(&self) -> f64 {
        let target = 180.0 - self.dihedral_deg;
        let mut worst: f64 = 0.0;
        for k in 0..8 {
            let p = self.point(TAU * f64::from(k) / 8.0);
            let n1 = self.field.a.gradient(p);
            let g2 = self.field.b.gradient(p);
            let Some(measured) = angle_deg(n1, [-g2[0], -g2[1], -g2[2]]) else {
                return f64::INFINITY;
            };
            worst = worst.max((measured - target).abs());
        }
        worst
    }

    /// Origin of the cube window, centred on the seam point at `t = 0`.
    fn window_origin(&self) -> [f64; 3] {
        let c = self.point(0.0);
        [
            c[0] - WINDOW * 0.5,
            c[1] - WINDOW * 0.5,
            c[2] - WINDOW * 0.5,
        ]
    }

    /// Half-width, in the circle's own parameter, of the arc that can possibly
    /// land inside the window.
    ///
    /// Two **necessary** conditions on the same `t`, both from the window being
    /// a cube of side [`WINDOW`] centred on `point(0.0)`: the `y` span forces
    /// `radius·cos t >= radius − WINDOW/2` and the `z` span forces
    /// `|radius·sin t| <= WINDOW/2`. Necessary, so a sample outside this arc
    /// would have failed the `inside` test anyway — which is what makes
    /// [`Self::seam_cells`] identical to the bench's full-circle walk rather
    /// than an approximation of it. Checked, not argued: the startup self-check
    /// holds `seam_cells` against the committed CSV.
    fn window_arc(&self) -> f64 {
        let half = WINDOW * 0.5;
        if self.radius <= half {
            // The circle is no wider than the window's own half-width, so no
            // `t` is excluded and the walk is the whole circle.
            return TAU;
        }
        let by_y = (1.0 - half / self.radius).acos();
        let by_z = (half / self.radius).asin();
        by_y.min(by_z)
    }

    /// Cells of the window grid that the seam circle enters.
    fn seam_cells(&self, origin: [f64; 3], cell: f64, cells: u32) -> usize {
        let extent = f64::from(cells);
        let arc = self.window_arc();
        let step = TAU / CIRCLE_SAMPLES as f64;
        // `+ 1` of slack so the truncation to whole samples can only ever widen
        // the arc, never clip a sample the full walk would have taken.
        let reach = (arc / step).ceil() as usize + 1;

        let mut hit: Vec<[i32; 3]> = Vec::new();
        let mut visit = |k: usize| {
            let p = self.point(TAU * k as f64 / CIRCLE_SAMPLES as f64);
            let mut index = [0i32; 3];
            for axis in 0..3 {
                let f = (p[axis] - origin[axis]) / cell;
                if f < 0.0 || f >= extent {
                    return;
                }
                index[axis] = f as i32;
            }
            hit.push(index);
        };
        if reach * 2 >= CIRCLE_SAMPLES {
            for k in 0..CIRCLE_SAMPLES {
                visit(k);
            }
        } else {
            for k in 0..=reach {
                visit(k);
            }
            for k in (CIRCLE_SAMPLES - reach)..CIRCLE_SAMPLES {
                visit(k);
            }
        }
        hit.sort_unstable();
        hit.dedup();
        hit.len()
    }
}

/// The differencing step `Sdf::gradient` uses at `p`.
fn step(p: [f64; 3]) -> f64 {
    <f64 as Real>::DIFF_STEP * p[0].abs().max(p[1].abs()).max(p[2].abs()).max(1.0)
}

/// The composed value at `q`, and whether branch `A` was the active one.
///
/// Derived from the two branch samples rather than read from
/// `Difference::sample`, and bit-identical to it: that function *is*
/// `if fa >= -fb { fa } else { -fb }`. Computing both here is what makes the
/// straddle test exact instead of a Lipschitz bound on the margin.
fn value(field: &Fixture, q: [f64; 3]) -> (f64, bool) {
    let fa = field.a.sample(q);
    let fb = field.b.sample(q);
    if fa >= -fb { (fa, true) } else { (-fb, false) }
}

/// What the six-sample stencil at `p` returned, and whether it crossed branches.
struct Probe {
    gradient: [f64; 3],
    straddles: bool,
}

/// The crate's default central difference, re-implemented here.
///
/// It has to be re-implemented: [`Difference`] **overrides** `Sdf::gradient`
/// with the analytic active-branch gradient, so calling `field.gradient` would
/// return the thing this example is measuring *against*. The formula is
/// `sdf.rs`'s, step for step.
fn probe(field: &Fixture, p: [f64; 3]) -> Probe {
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

/// How far `p` sits from the branch surface `s = f_A + f_B = 0`, in stencil
/// widths.
///
/// The straddle test flips an axis when this drops below `1`, since `s` is
/// affine to `O(h²)` over a stencil this small — so it is the same quantity the
/// classification uses, reported as a continuum instead of a bit. It is what
/// explains a low ratio: well above `1` the grid has not resolved the crease to
/// within a differencing step, and no vertex is near the corner case.
fn seam_offset_stencils(field: &Fixture, p: [f64; 3]) -> Option<f64> {
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

/// Median of a sample, by `total_cmp` so a NaN surfaces rather than vanishing.
fn median(mut v: Vec<f64>) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_unstable_by(f64::total_cmp);
    let mid = v.len() / 2;
    if v.len().is_multiple_of(2) {
        Some((v[mid - 1] + v[mid]) * 0.5)
    } else {
        Some(v[mid])
    }
}

/// The bench's `show`: the format every floating-point column of the CSV is
/// written in, so a comparison can be made as text.
fn show(v: Option<f64>) -> String {
    v.map_or_else(|| String::from("n/a"), |x| format!("{x:.6e}"))
}

// ─── one measurement ────────────────────────────────────────────────────────

/// One (dihedral, resolution) row, measured live.
#[derive(Default, Clone)]
struct Row {
    dihedral_deg: f64,
    samples: u32,
    cell_size: f64,
    seam_cells: usize,
    vertices: usize,
    triangles: usize,
    straddling: usize,
    straddling_max_deg: Option<f64>,
    straddling_mean_deg: Option<f64>,
    non_straddling: usize,
    non_straddling_mean_deg: Option<f64>,
    degenerate: usize,
    seam_offset_stencils: Option<f64>,
    separation: f64,
    seam_radius: f64,
    extract_ms: f64,
}

impl Row {
    fn bound_deg(&self) -> f64 {
        (180.0 - self.dihedral_deg) * 0.5
    }

    fn ratio(&self) -> Option<f64> {
        self.straddling_max_deg.map(|e| e / self.bound_deg())
    }

    fn per_seam_cell(&self) -> Option<f64> {
        (self.seam_cells > 0).then(|| self.straddling as f64 / self.seam_cells as f64)
    }

    fn share_of_vertices(&self) -> Option<f64> {
        (self.vertices > 0).then(|| self.straddling as f64 / self.vertices as f64)
    }

    /// Every column of `docs/experiments/p-56.csv` this row can reproduce,
    /// formatted exactly the way the bench wrote it.
    ///
    /// Four columns are deliberately absent. `scaling_exponent` is a fit across
    /// four resolutions, `sweep_median_ratio` is a median over all 28 rows and
    /// `seam_angle_residual_deg` is the worst over all seven fixtures — none of
    /// which one row can produce, and the first two would need a 257³ sweep at
    /// startup. The residual is checked separately, over the same seven
    /// dihedrals the bench used.
    fn ledger_fields(&self) -> Vec<(&'static str, String)> {
        vec![
            ("seam_cells", self.seam_cells.to_string()),
            ("vertices", self.vertices.to_string()),
            ("straddling_vertices", self.straddling.to_string()),
            ("straddling_max_error_deg", show(self.straddling_max_deg)),
            ("predicted_bound_deg", format!("{:.6}", self.bound_deg())),
            ("worst_over_bound_ratio", show(self.ratio())),
            (
                "within_bound",
                self.straddling_max_deg.map_or_else(
                    || String::from("n/a"),
                    |e| (e <= self.bound_deg()).to_string(),
                ),
            ),
            (
                "non_straddling_mean_error_deg",
                show(self.non_straddling_mean_deg),
            ),
            ("non_straddling_vertices", self.non_straddling.to_string()),
            ("straddling_mean_error_deg", show(self.straddling_mean_deg)),
            ("degenerate_normals", self.degenerate.to_string()),
            ("cell_size", format!("{:.6e}", self.cell_size)),
            ("separation", format!("{:.9}", self.separation)),
            ("seam_radius", format!("{:.9}", self.seam_radius)),
            ("seam_offset_stencils", show(self.seam_offset_stencils)),
            ("straddling_per_seam_cell", show(self.per_seam_cell())),
            (
                "straddling_share_of_vertices",
                show(self.share_of_vertices()),
            ),
            ("window_side", format!("{WINDOW}")),
        ]
    }
}

/// One straddling vertex, already in the view frame so drawing costs nothing.
struct Straddle {
    at: Vec3,
    /// The analytic active-branch gradient — the right answer.
    analytic: Vec3,
    /// What the central difference returned — the wrong one.
    difference: Vec3,
    error_deg: f64,
}

/// The rigid frame the picture is drawn in: `e₁` bisects the two outward
/// normals, `e₃` is the crease tangent `ẑ`, and the window becomes a unit cube.
///
/// Both branch normals lie in the `z = 0` plane at the seam point `t = 0`, so
/// this is a rotation about `z`, a translation and a uniform scale — nothing
/// that could change an angle.
#[derive(Clone, Copy)]
struct View {
    centre: [f64; 3],
    e1: [f64; 3],
    e2: [f64; 3],
}

/// How much the window is magnified: `1 / WINDOW`, so the window is a unit cube.
const VIEW_SCALE: f64 = 1.0 / WINDOW;

impl View {
    /// The frame the seam point at `t = 0` defines.
    fn of(seam: &Seam) -> Self {
        let centre = seam.point(0.0);
        let n1 = seam.field.a.gradient(centre);
        let g2 = seam.field.b.gradient(centre);
        let outward = [-g2[0], -g2[1], -g2[2]];
        let bisector = unit([n1[0] + outward[0], n1[1] + outward[1], n1[2] + outward[2]])
            .unwrap_or([1.0, 0.0, 0.0]);
        let e1 = unit([bisector[0], bisector[1], 0.0]).unwrap_or([1.0, 0.0, 0.0]);
        // `ẑ × e₁`, so `(e₁, e₂, ẑ)` is right-handed.
        let e2 = [-e1[1], e1[0], 0.0];
        Self { centre, e1, e2 }
    }

    fn point(&self, p: [f64; 3]) -> Vec3 {
        let d = [
            p[0] - self.centre[0],
            p[1] - self.centre[1],
            p[2] - self.centre[2],
        ];
        Vec3::new(
            (dot(d, self.e1) * VIEW_SCALE) as f32,
            (dot(d, self.e2) * VIEW_SCALE) as f32,
            (d[2] * VIEW_SCALE) as f32,
        )
    }

    fn direction(&self, v: [f64; 3]) -> Vec3 {
        Vec3::new(dot(v, self.e1) as f32, dot(v, self.e2) as f32, v[2] as f32).normalize_or_zero()
    }
}

impl Default for View {
    /// The identity frame. Only ever seen by [`Live`]'s default, which no
    /// system reads: `rebuild` runs ahead of every consumer in the chain.
    fn default() -> Self {
        Self {
            centre: [0.0; 3],
            e1: [1.0, 0.0, 0.0],
            e2: [0.0, 1.0, 0.0],
        }
    }
}

/// Everything one rebuild produced.
#[derive(Default)]
struct Built {
    row: Row,
    straddles: Vec<Straddle>,
    /// Index into `straddles` of the vertex `straddling_max_error_deg` is at.
    worst: usize,
    /// The exact crease, clipped to the window, in the view frame.
    crease: Vec<Vec3>,
    /// The eight corners of the sampled window, in the view frame.
    corners: [Vec3; 8],
    /// The frame the mesh, the rays, the crease and the box are all drawn in.
    view: View,
}

/// Mesh the window, classify every vertex, and lay out the picture.
fn build(
    seam: &Seam,
    samples: u32,
    dc: &mut DualContouring<f64>,
    buffer: &mut MeshBuffer<f64>,
) -> Built {
    let origin = seam.window_origin();
    let cell_size = WINDOW / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("window grid fits u32");
    let view = View::of(seam);

    buffer.reset();
    let started = Instant::now();
    if let Err(e) = dc.extract(&seam.field, &shape, origin, cell_size, buffer) {
        error!("dual contouring failed on the seam window: {e}");
    }
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;
    let seam_cells = seam.seam_cells(origin, cell_size, samples - 1);

    let mut straddles: Vec<Straddle> = Vec::new();
    let mut straddling_sum = 0.0;
    let mut smooth_sum = 0.0;
    let mut smooth_count = 0usize;
    let mut degenerate = 0usize;
    let mut offsets: Vec<f64> = Vec::with_capacity(buffer.positions.len());
    for &p in &buffer.positions {
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
            straddling_sum += error;
            straddles.push(Straddle {
                at: view.point(p),
                analytic: view.direction(analytic),
                difference: view.direction(probed.gradient),
                error_deg: error,
            });
        } else {
            smooth_sum += error;
            smooth_count += 1;
        }
    }

    // The `seam_cells` vertices nearest the branch surface are the seam
    // population; anything past that is the smooth sheet and would drown the
    // statistic in numbers about the wrong vertices.
    offsets.sort_unstable_by(f64::total_cmp);
    offsets.truncate(seam_cells);

    let straddling = straddles.len();
    let worst = straddles
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.error_deg.total_cmp(&b.1.error_deg))
        .map_or(0, |(i, _)| i);

    let row = Row {
        dihedral_deg: seam.dihedral_deg,
        samples,
        cell_size,
        seam_cells,
        vertices: buffer.positions.len(),
        triangles: buffer.indices.len() / 3,
        straddling,
        straddling_max_deg: straddles.get(worst).map(|s| s.error_deg),
        straddling_mean_deg: (straddling > 0).then(|| straddling_sum / straddling as f64),
        non_straddling: smooth_count,
        non_straddling_mean_deg: (smooth_count > 0).then(|| smooth_sum / smooth_count as f64),
        degenerate,
        seam_offset_stencils: median(offsets),
        separation: seam.separation,
        seam_radius: seam.radius,
        extract_ms,
    };

    Built {
        row,
        straddles,
        worst,
        crease: crease_polyline(seam, &view, origin, cell_size, samples - 1),
        corners: window_corners(&view, origin),
        view,
    }
}

/// The exact seam circle, sampled over the arc that reaches the window and
/// clipped to it.
///
/// From the fixture's closed form rather than from the mesh: the white line has
/// to be the crease the bound is about, not the extractor's guess at it.
fn crease_polyline(seam: &Seam, view: &View, origin: [f64; 3], cell: f64, cells: u32) -> Vec<Vec3> {
    const STEPS: usize = 192;
    let extent = f64::from(cells);
    let arc = seam.window_arc().min(TAU * 0.5);
    let mut out = Vec::with_capacity(STEPS + 1);
    for k in 0..=STEPS {
        let t = -arc + 2.0 * arc * k as f64 / STEPS as f64;
        let p = seam.point(t);
        let inside = (0..3).all(|axis| {
            let f = (p[axis] - origin[axis]) / cell;
            f >= 0.0 && f < extent
        });
        if inside {
            out.push(view.point(p));
        }
    }
    out
}

/// The eight corners of the sampled window, in the view frame.
///
/// Corner `i` has bit `a` of `i` selecting the far side on axis `a`, matching
/// the extractors' own corner indexing.
fn window_corners(view: &View, origin: [f64; 3]) -> [Vec3; 8] {
    let mut out = [Vec3::ZERO; 8];
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = view.point([
            origin[0] + if i & 1 == 0 { 0.0 } else { WINDOW },
            origin[1] + if i & 2 == 0 { 0.0 } else { WINDOW },
            origin[2] + if i & 4 == 0 { 0.0 } else { WINDOW },
        ]);
    }
    out
}

// ─── the ledger, compiled in ────────────────────────────────────────────────

/// P-56's committed artefact, embedded at compile time.
///
/// `include_str!` rather than a runtime read, and rather than transcribed
/// constants: the path is resolved against this source file so the check cannot
/// be broken by a working directory, and a number that lives only here could
/// drift from the CSV without anything saying so. Costs 6 KB in the binary.
const LEDGER_CSV: &str = include_str!("../../docs/experiments/p-56.csv");

/// The parsed CSV. Rows are `&'static str` slices of [`LEDGER_CSV`], so this
/// owns no text.
#[derive(Resource)]
struct Ledger {
    header: Vec<&'static str>,
    rows: Vec<Vec<&'static str>>,
    /// `sweep_median_ratio`, the ledger's median tightness over its measured
    /// rows — the reference tick on the bar.
    median_ratio: f64,
    median_text: String,
    measured_rows: usize,
}

impl Ledger {
    fn parse() -> Self {
        let mut lines = LEDGER_CSV
            .lines()
            .map(str::trim_end)
            .filter(|l| !l.is_empty() && !l.starts_with('#'));
        let header: Vec<&'static str> = lines
            .next()
            .map(|l| l.split(',').collect())
            .unwrap_or_default();
        let rows: Vec<Vec<&'static str>> = lines.map(|l| l.split(',').collect()).collect();

        let mut out = Self {
            header,
            rows,
            median_ratio: 0.0,
            median_text: String::from("n/a"),
            measured_rows: 0,
        };
        // Taken from the artefact rather than hardcoded, for the same reason the
        // CSV is embedded rather than transcribed. `sweep_median_ratio` is a
        // sweep-wide statistic repeated on every row, so the first row carries
        // it.
        let median = out
            .rows
            .first()
            .and_then(|r| out.field(r, "sweep_median_ratio"));
        if let Some(text) = median {
            out.median_text = text.to_string();
            out.median_ratio = text.parse().unwrap_or(0.0);
        }
        let measured = out
            .rows
            .iter()
            .filter(|r| {
                out.field(r, "worst_over_bound_ratio")
                    .is_some_and(|v| v != "n/a")
            })
            .count();
        out.measured_rows = measured;
        out
    }

    fn field(&self, row: &[&'static str], name: &str) -> Option<&'static str> {
        let index = self.header.iter().position(|h| *h == name)?;
        row.get(index).copied()
    }

    /// The row for one (dihedral, resolution), if the ledger measured it.
    fn row(&self, dihedral_deg: f64, samples: u32) -> Option<&Vec<&'static str>> {
        let want_theta = format!("{dihedral_deg:.1}");
        let want_samples = samples.to_string();
        self.rows.iter().find(|r| {
            self.field(r, "dihedral_deg") == Some(want_theta.as_str())
                && self.field(r, "samples_per_axis") == Some(want_samples.as_str())
        })
    }

    /// Every column of `row` that `measured` also produces, compared as text.
    ///
    /// Text rather than a tolerance, because "reproduces the ledger" either
    /// means the digits or it means nothing.
    fn diff(&self, row: &[&'static str], measured: &Row) -> Vec<(&'static str, String, String)> {
        measured
            .ledger_fields()
            .into_iter()
            .filter_map(|(name, mine)| {
                let theirs = self.field(row, name)?;
                (theirs != mine).then(|| (name, mine, theirs.to_string()))
            })
            .collect()
    }
}

// ─── the startup self-check ─────────────────────────────────────────────────

/// The three 129³ rows of `docs/experiments/p-56.csv`, re-measured before the
/// window opens and compared against the CSV column by column.
///
/// A wrong answer must be **loud and must not take the window down with it**:
/// it logs at `error!` and the HUD says so, because a demo a stranger runs is
/// not the place for an assertion.
#[derive(Resource)]
struct Replay {
    /// Every compared column of every row matched to the digit.
    exact: bool,
    /// One HUD line naming the three headline numbers and the verdict.
    line: String,
    /// The fixture check: worst `|angle(n₁, −n₂) − (180 − theta)|` over all
    /// seven ledger dihedrals, and whether it matches the CSV's own column.
    residual_deg: f64,
    residual_matches: bool,
}

impl Replay {
    fn run(ledger: &Ledger) -> Self {
        // ── the construction check, first, because it can void everything ──
        let mut residual_deg: f64 = 0.0;
        for theta in DIHEDRALS {
            let seam = Seam::new(theta);
            let r = seam.residual_deg();
            info!(
                "E-311 fixture  theta {theta:>6.1} deg  dist {:.9}  seam x {:.9} r {:.9}  \
                 residual {r:.3e} deg",
                seam.separation, seam.plane_x, seam.radius
            );
            residual_deg = residual_deg.max(r);
        }
        let residual_text = format!("{residual_deg:.3e}");
        let ledger_residual = ledger
            .rows
            .first()
            .and_then(|r| ledger.field(r, "seam_angle_residual_deg"));
        let residual_matches = ledger_residual == Some(residual_text.as_str());
        if residual_deg >= SEAM_TOLERANCE_DEG {
            error!(
                "E-311 fixture is wrong: worst seam-angle residual {residual_text} deg exceeds \
                 {SEAM_TOLERANCE_DEG:e} deg, so no number below is about the dihedral it claims"
            );
        }
        if !residual_matches {
            error!(
                "E-311 seam_angle_residual_deg: measured {residual_text}, \
                 p-56.csv says {}",
                ledger_residual.unwrap_or("<absent>")
            );
        }
        info!(
            "E-311 fixture  worst residual {residual_text} deg over {} dihedrals, \
             p-56.csv seam_angle_residual_deg {} - {}",
            DIHEDRALS.len(),
            ledger_residual.unwrap_or("<absent>"),
            if residual_matches {
                "MATCH"
            } else {
                "DISAGREES"
            },
        );

        // ── the three rows ──
        let mut dc = DualContouring::<f64>::new();
        let mut buffer = MeshBuffer::<f64>::new();
        let mut exact = residual_matches && residual_deg < SEAM_TOLERANCE_DEG;
        let mut headline: Vec<String> = Vec::new();
        let mut compared = 0usize;
        for theta in SELF_CHECK {
            let seam = Seam::new(theta);
            let built = build(&seam, LEDGER_SAMPLES, &mut dc, &mut buffer);
            let row = &built.row;
            let Some(ledger_row) = ledger.row(theta, LEDGER_SAMPLES) else {
                error!("E-311 p-56.csv has no {theta:.1} deg row at {LEDGER_SAMPLES}^3");
                exact = false;
                continue;
            };
            let mismatches = ledger.diff(ledger_row, row);
            let columns = row.ledger_fields().len();
            compared += columns;
            info!(
                "E-311 replay   theta {theta:>6.1} deg  {LEDGER_SAMPLES}^3  \
                 worst {:.6} deg  bound {:.6} deg  ratio {:.6}  straddling {}/{}  \
                 {} of {columns} columns match p-56.csv",
                row.straddling_max_deg.unwrap_or(f64::NAN),
                row.bound_deg(),
                row.ratio().unwrap_or(f64::NAN),
                row.straddling,
                row.vertices,
                columns - mismatches.len(),
            );
            for (name, mine, theirs) in &mismatches {
                error!("E-311 replay   theta {theta:.1} {name}: measured {mine}, csv {theirs}");
                exact = false;
            }
            headline.push(format!(
                "{theta:.0}:{:.6}/{:.1}",
                row.straddling_max_deg.unwrap_or(f64::NAN),
                row.bound_deg(),
            ));
        }

        // Two lines, not one: the HUD gets the verdict and the count, the log
        // gets the three headline numbers. The one-line version ran to 108
        // characters and reached across the frame into the ray fan, and a HUD
        // line lying over the subject is the same defect as a wrapped one.
        let line = format!(
            "check     p-56.csv {LEDGER_SAMPLES}^3 replay {}, {compared} columns",
            if exact {
                "EXACT to the digit"
            } else {
                "FAILED - see the log"
            },
        );
        info!("E-311 {line}   {}", headline.join("  "));
        Self {
            exact,
            line,
            residual_deg,
            residual_matches,
        }
    }
}

// ─── the sweep ──────────────────────────────────────────────────────────────

/// Degrees between consecutive dihedrals on the ramp.
///
/// `2.5` and not a rounder number: `30 + 2.5k` lands **exactly** on all seven
/// dihedrals P-56 measured — 30, 60, 90, 120, 150, 165 and 175 — so the sweep
/// passes through every row of the ledger rather than near them. 59 entries.
const RAMP_STEP: f64 = 2.5;

/// Extra schedule entries spent on each of the seven ledger dihedrals.
///
/// `59 + 7 · 3 = 80`, which is the harness's own default
/// `ISOMESH_CAPTURE_FRAMES` — so a default capture is exactly one sweep and
/// stops for four frames at each measured row.
const LEDGER_HOLD: usize = 3;

/// Seconds for one pass through the schedule, when nobody is capturing.
const SWEEP_SECONDS: f32 = 20.0;

/// The dihedrals the sweep visits, in order.
#[derive(Resource)]
struct Schedule(Vec<f64>);

impl Schedule {
    fn build() -> Self {
        let mut out = Vec::new();
        let steps = ((175.0 - 30.0) / RAMP_STEP).round() as usize;
        for k in 0..=steps {
            let theta = 30.0 + RAMP_STEP * k as f64;
            out.push(theta);
            if DIHEDRALS.iter().any(|d| (d - theta).abs() < 1e-12) {
                for _ in 0..LEDGER_HOLD {
                    out.push(theta);
                }
            }
        }
        Self(out)
    }

    fn at(&self, index: usize) -> f64 {
        self.0
            .get(index % self.0.len().max(1))
            .copied()
            .unwrap_or(90.0)
    }

    /// First schedule entry sitting on `dihedral_deg`.
    fn find(&self, dihedral_deg: f64) -> usize {
        self.0
            .iter()
            .position(|d| (d - dihedral_deg).abs() < 1e-12)
            .unwrap_or(0)
    }
}

/// Samples per axis, fixed at startup from `ISOMESH_SAMPLES`.
#[derive(Resource)]
struct Resolution(u32);

/// A ledger row pinned by `ISOMESH_FIELD`, which overrides the sweep.
///
/// The harness's contract is that anything a capture depends on is reachable
/// from the environment; without this a still of one row could only be produced
/// by holding a key down.
#[derive(Resource)]
struct Pinned(Option<f64>);

/// Where the sweep is.
#[derive(Resource, Default)]
struct Cursor {
    index: usize,
    seconds: f32,
}

/// The current measurement and its picture.
#[derive(Resource, Default)]
struct Live(Built);

/// The rig.
#[derive(Resource)]
struct Demo {
    dc: DualContouring<f64>,
    buffer: MeshBuffer<f64>,
    mesh: Option<Handle<Mesh>>,
    surface: Entity,
    /// `(dihedral bits, samples)` of the last rebuild, so a held frame is free.
    last: Option<(u64, u32)>,
}

/// The comb of ray pairs, the crease and the window box. Its own group so the
/// depth bias can pull it in front of the surface without dragging the shared
/// wireframe along.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct CombGizmos;

/// The worst straddling vertex, biased harder still and drawn wider, so the
/// vertex the bar is about is never lost in the comb.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct HeroGizmos;

// ─── framing ────────────────────────────────────────────────────────────────

/// Radians the camera swings off the wedge's bisector, toward `e₂`.
///
/// Straight down the bisector the two sheets are equally lit and equally
/// foreshortened, which makes the crease hard to read as a fold; swung too far
/// and the far sheet goes edge-on and there is only one surface on screen.
/// Measured on a 1280x720 capture at all three ledger rows: `0.32` keeps both
/// sheets front-facing at 90° and 120° — cosines `0.51` and `0.27` against the
/// two normals — and is the widest swing that still shows both at 175°.
const CAMERA_SWING: f32 = 0.32;

/// Radians the camera is tilted off the `e₁e₂` plane toward the crease.
///
/// **Mostly along the crease, and that is the whole composition.** The pair of
/// rays at a straddling vertex lies in the plane spanned by the two branch
/// normals, which is exactly the plane perpendicular to the crease — so looking
/// along the crease is the only view in which the angle between them projects to
/// its true size. It costs a foreshortened crease, `cos` of this off the crease
/// axis, and [`CAMERA_RADIUS`] pays for that by zooming in.
const CAMERA_TILT: f32 = 0.95;

/// Orbit radius, in window widths.
///
/// The window is a unit cube in the view frame and Bevy's default vertical FOV
/// is 45°, so the half-height in view units is `radius · tan(22.5°)`. The
/// subject is the fan of ray pairs along the crease, whose projected extent is
/// `0.58` of the window plus a ray length either side — about `0.8` — and this
/// puts that at some 60% of the frame height. The window box then overflows the
/// frame top and bottom, deliberately: the box is context and the fan is the
/// finding.
const CAMERA_RADIUS: f32 = 1.60;

/// Where the subject sits in frame, as a fraction of the orbit radius: right,
/// and up out of the bar's way.
///
/// **The HUD is fifteen lines in the upper left and the bar is along the
/// bottom.** Centring the wedge photographs the argument with its evidence
/// hidden, which is E-112's lesson.
const SUBJECT_OFFSET: Vec2 = Vec2::new(0.20, -0.04);

/// Length of a comb ray, in window widths.
const COMB_LEN: f32 = 0.11;

/// Length of the worst vertex's rays.
const HERO_LEN: f32 = 0.26;

/// How many ray pairs the comb draws.
///
/// At 129³ there are up to 132 straddling vertices along a crease one window
/// wide, i.e. one every `1/128` of the frame — a fan that dense is a filled
/// triangle rather than a set of pairs, and the angle it is supposed to show
/// disappears. Every straddling vertex is counted on the HUD; every
/// `stride`-th one is drawn, and the HUD says which.
const COMB_RAYS: usize = 20;

const ANALYTIC: Color = Color::srgb(0.25, 0.88, 0.95);
const DIFFERENCE: Color = Color::srgb(0.99, 0.42, 0.86);
const WEDGE: Color = Color::srgb(0.99, 0.78, 0.30);
const CREASE: Color = Color::srgb(0.96, 0.97, 1.0);
const BOX_EDGE: Color = Color::srgb(0.34, 0.37, 0.44);

// ─── the bar ────────────────────────────────────────────────────────────────

/// Width and height of the bar's track, in logical pixels.
///
/// 460 of 1280 is 36% of the frame, which survives `record_gif.sh`'s downscale
/// to 900 wide as 323 px — wide enough that the 2.5% gap between the ledger's
/// median tightness and a full bar is still a visible gap.
const TRACK_W: f32 = 460.0;
const TRACK_H: f32 = 26.0;

const FILL_GOOD: Color = Color::srgb(0.30, 0.85, 0.42);
const FILL_SLACK: Color = Color::srgb(0.99, 0.74, 0.22);
const FILL_BREACH: Color = Color::srgb(0.94, 0.26, 0.22);
const TRACK_BG: Color = Color::srgb(0.11, 0.12, 0.16);

/// Above this the bar is green: the sampled vertices reached the corner case.
const TIGHT: f64 = 0.9;

/// Everything the bar is made of, so `nohud` hides all of it.
#[derive(Component)]
struct HudPanel;

/// The filled part of the track.
#[derive(Component)]
struct MeterFill;

/// The big number under the track.
#[derive(Component)]
struct MeterValue;

// ─── app ────────────────────────────────────────────────────────────────────

fn main() {
    let samples = samples_override()
        .unwrap_or(LEDGER_SAMPLES)
        .clamp(MIN_SAMPLES, MAX_SAMPLES);
    let ledger = Ledger::parse();
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-311 seam normal bound".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<CombGizmos>()
        .init_gizmo_group::<HeroGizmos>()
        .insert_resource(Resolution(samples))
        // Replayed here rather than in `setup`, so every system can take it as a
        // plain `Res`, and after `add_plugins` so the log subscriber is already
        // installed and a failed replay is heard.
        .insert_resource(Replay::run(&ledger))
        .insert_resource(ledger)
        .insert_resource(Schedule::build())
        .insert_resource(Pinned(pinned_row()))
        .init_resource::<Cursor>()
        .init_resource::<Live>()
        .add_systems(Startup, setup)
        // **`PreUpdate`, and that is load-bearing rather than a preference.**
        // The harness's `update_hud` renders `DemoStats` and its
        // `capture_sequence` both takes the screenshot and advances
        // `Capture::taken`, and `Update` gives no ordering against either — so
        // in `Update` the HUD renders a frame-old row while the mesh is current
        // and the sweep reads `taken` on either side of the increment. For a
        // demo whose whole claim is "the number on screen is this picture" that
        // is the one defect that matters. After `InputSystems` so a keypress is
        // seen in the frame it happened.
        .add_systems(
            PreUpdate,
            (advance, rebuild, draw_scene, paint_meter, report)
                .chain()
                .after(bevy::input::InputSystems),
        )
        .run();
}

/// The ledger row `ISOMESH_FIELD` asks for, if it asks for one.
fn pinned_row() -> Option<f64> {
    let raw = std::env::var("ISOMESH_FIELD").ok()?;
    match raw.trim().parse::<usize>() {
        Ok(index) if index < LEDGER_PINS.len() => Some(LEDGER_PINS[index]),
        _ => {
            error!(
                "ISOMESH_FIELD={raw} is not one of 0..{} - 90, 120, 175 degrees",
                LEDGER_PINS.len()
            );
            None
        }
    }
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut config: ResMut<GizmoConfigStore>,
    mut camera: Query<&mut OrbitCamera>,
    ledger: Res<Ledger>,
) {
    let (comb, _) = config.config_mut::<CombGizmos>();
    comb.line.width = 1.5;
    comb.depth_bias = -0.25;
    let (hero, _) = config.config_mut::<HeroGizmos>();
    hero.line.width = 3.4;
    hero.depth_bias = -0.7;

    // Dark, and darker than this repo's usual surface grey. The two sheets fill
    // most of the frame at every theta in the sweep, and the rays and the HUD
    // text are the evidence -- so the rock loses. `double_sided` because the
    // camera crosses the plane of one sheet as the dihedral closes.
    let surface_material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.24, 0.27, 0.34),
        perceptual_roughness: 0.68,
        metallic: 0.05,
        double_sided: true,
        cull_mode: None,
        ..default()
    });

    // `Mesh3d::default()` names no asset, so nothing is drawn and nothing is
    // uploaded until the first rebuild. An empty mesh would be worse than
    // nothing: `MeshAllocator` skips a zero-byte vertex buffer and then copies
    // into it anyway, once per frame, in red.
    let surface = commands
        .spawn((
            Mesh3d::default(),
            MeshMaterial3d(surface_material),
            DemoMesh,
        ))
        .id();

    // The subject is always the view frame's origin, so the camera is set once
    // and never moves again -- which is also what makes the GIF compress.
    let dir = Vec3::new(
        CAMERA_TILT.cos() * CAMERA_SWING.cos(),
        CAMERA_TILT.cos() * CAMERA_SWING.sin(),
        CAMERA_TILT.sin(),
    );
    for mut orbit in &mut camera {
        orbit.yaw = dir.z.atan2(dir.x);
        orbit.pitch = dir.y.asin();
        orbit.radius = CAMERA_RADIUS;
        // The camera's own basis, from the same yaw/pitch the harness's
        // `orbit_camera` builds its transform from. It places the eye at
        // `focus + dir * radius`, so the view direction is `-dir` and a focus
        // moved along `-right` puts the subject right of centre.
        let forward = -dir;
        let right = forward.cross(Vec3::Y).normalize_or_zero();
        let up = right.cross(forward).normalize_or_zero();
        orbit.focus =
            -right * (SUBJECT_OFFSET.x * CAMERA_RADIUS) + up * (SUBJECT_OFFSET.y * CAMERA_RADIUS);
    }

    spawn_meter(&mut commands, ledger.median_ratio);

    commands.insert_resource(Demo {
        dc: DualContouring::<f64>::new(),
        buffer: MeshBuffer::<f64>::new(),
        mesh: None,
        surface,
        last: None,
    });
}

/// Decide which dihedral this frame is about.
///
/// Under capture the schedule advances one entry per captured frame, so a clip
/// of any length is a sweep rather than a still. Interactively it advances on
/// wall-clock time, and the digits freeze it on a ledger row.
fn advance(
    time: Res<Time>,
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    schedule: Res<Schedule>,
    pinned: Res<Pinned>,
    mut flags: ResMut<ViewFlags>,
    mut cursor: ResMut<Cursor>,
) {
    if let Some(theta) = pinned.0 {
        cursor.index = schedule.find(theta);
        return;
    }

    for (key, theta) in [KeyCode::Digit1, KeyCode::Digit2, KeyCode::Digit3]
        .into_iter()
        .zip(LEDGER_PINS)
    {
        if keys.just_pressed(key) {
            cursor.index = schedule.find(theta);
            cursor.seconds = 0.0;
            flags.paused = true;
        }
    }
    if keys.just_pressed(KeyCode::KeyX) {
        cursor.index = 0;
        cursor.seconds = 0.0;
        flags.paused = false;
    }
    if flags.paused {
        return;
    }

    if capture.is_active() {
        cursor.index = capture.taken as usize;
        return;
    }
    let per_entry = SWEEP_SECONDS / schedule.0.len().max(1) as f32;
    cursor.seconds += time.delta_secs();
    while cursor.seconds >= per_entry {
        cursor.seconds -= per_entry;
        cursor.index += 1;
    }
}

/// Mesh the window and classify its vertices — only when the answer would
/// change.
#[allow(clippy::too_many_arguments)]
fn rebuild(
    schedule: Res<Schedule>,
    resolution: Res<Resolution>,
    cursor: Res<Cursor>,
    flags: Res<ViewFlags>,
    mut demo: ResMut<Demo>,
    mut live: ResMut<Live>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut commands: Commands,
) {
    let dihedral = schedule.at(cursor.index);
    let key = (dihedral.to_bits(), resolution.0);
    if demo.last == Some(key) && !flags.remesh_requested {
        return;
    }
    demo.last = Some(key);

    let seam = Seam::new(dihedral);
    let Demo {
        dc, buffer, mesh, ..
    } = &mut *demo;
    let built = build(&seam, resolution.0, dc, buffer);

    let handle = meshes.add(to_mesh(buffer, &built.view));
    if let Some(old) = mesh.replace(handle.clone()) {
        meshes.remove(&old);
    }
    commands.entity(demo.surface).insert(Mesh3d(handle));
    live.0 = built;
}

/// The `f64` extraction as a Bevy mesh, in the view frame.
///
/// Cast rather than re-extracted in `f32`: the numbers on the HUD are M-350's
/// and they are `f64` numbers, so the mesh the picture is drawn from has to be
/// the one they were computed on. The view frame is applied here rather than as
/// a `Transform` so the vertices handed to `f32` are order-one instead of a
/// `0.04`-wide detail at `x ≈ 0.9`, which is the difference between a crease and
/// a row of z-fighting.
fn to_mesh(buffer: &MeshBuffer<f64>, view: &View) -> Mesh {
    let mut builder = MeshBuilder::new();
    for (p, n) in buffer.positions.iter().zip(&buffer.normals) {
        builder.vertex(view.point(*p).to_array(), view.direction(*n).to_array());
    }
    for t in buffer.indices.chunks_exact(3) {
        let (Some(a), Some(b), Some(c)) = (t.first(), t.get(1), t.get(2)) else {
            continue;
        };
        builder.triangle(*a, *b, *c);
    }
    builder.into_mesh()
}

/// The comb of ray pairs, the crease, and the window box.
fn draw_scene(
    live: Res<Live>,
    flags: Res<ViewFlags>,
    mut comb: Gizmos<CombGizmos>,
    mut hero: Gizmos<HeroGizmos>,
) {
    if flags.grid {
        for i in 0..8usize {
            for bit in [1usize, 2, 4] {
                if i & bit == 0 {
                    comb.line(live.0.corners[i], live.0.corners[i | bit], BOX_EDGE);
                }
            }
        }
    }

    for pair in live.0.crease.windows(2) {
        comb.line(pair[0], pair[1], CREASE);
    }

    let total = live.0.straddles.len();
    if total == 0 {
        return;
    }
    let stride = total.div_ceil(COMB_RAYS).max(1);
    for (i, s) in live.0.straddles.iter().enumerate() {
        if i == live.0.worst || !i.is_multiple_of(stride) {
            continue;
        }
        comb.line(s.at, s.at + s.analytic * COMB_LEN, ANALYTIC);
        comb.line(s.at, s.at + s.difference * COMB_LEN, DIFFERENCE);
        // Two arcs rather than one, so the angle reads as a shaded wedge
        // instead of as a stray curve.
        for scale in [0.62f32, 0.86] {
            comb.short_arc_3d_between(
                s.at,
                s.at + s.analytic * (COMB_LEN * scale),
                s.at + s.difference * (COMB_LEN * scale),
                WEDGE.with_alpha(0.55),
            )
            .resolution(12);
        }
    }

    let Some(s) = live.0.straddles.get(live.0.worst) else {
        return;
    };
    hero.line(s.at, s.at + s.analytic * HERO_LEN, ANALYTIC);
    hero.line(s.at, s.at + s.difference * HERO_LEN, DIFFERENCE);
    for scale in [0.42f32, 0.58, 0.74, 0.90] {
        hero.short_arc_3d_between(
            s.at,
            s.at + s.analytic * (HERO_LEN * scale),
            s.at + s.difference * (HERO_LEN * scale),
            WEDGE,
        )
        .resolution(32);
    }
    hero.sphere(Isometry3d::from_translation(s.at), 0.014, CREASE);
}

// ─── the bar ────────────────────────────────────────────────────────────────

/// The bound meter, centred along the bottom edge where a GIF viewer looks.
///
/// `median_ratio` places the reference tick, and it is read from the CSV rather
/// than written here: a tick drawn at a literal is exactly the kind of number
/// that drifts away from the artefact it claims to quote.
fn spawn_meter(commands: &mut Commands, median_ratio: f64) {
    commands
        .spawn((
            Node {
                position_type: PositionType::Absolute,
                left: Val::Px(0.0),
                right: Val::Px(0.0),
                bottom: Val::Px(14.0),
                justify_content: JustifyContent::Center,
                ..default()
            },
            GlobalZIndex(4),
            HudPanel,
        ))
        .with_children(|parent| {
            parent
                .spawn((
                    Node {
                        flex_direction: FlexDirection::Column,
                        align_items: AlignItems::Center,
                        row_gap: Val::Px(5.0),
                        padding: UiRect::axes(Val::Px(12.0), Val::Px(9.0)),
                        ..default()
                    },
                    BackgroundColor(Color::srgba(0.03, 0.04, 0.07, 0.88)),
                ))
                .with_children(|column| {
                    column.spawn((
                        Text::new(
                            "worst straddling error against the predicted (180 - theta)/2\n\
                             cyan tick: the ledger's median tightness over its measured rows",
                        ),
                        TextFont {
                            font_size: FontSize::Px(12.0),
                            ..default()
                        },
                        // `NoWrap`: in a centring flex column the text measure
                        // is handed the container's whole width while the node's
                        // own height resolves before the wrap, so a soft wrap
                        // pushes the number off the bottom of the frame.
                        TextLayout {
                            linebreak: bevy::text::LineBreak::NoWrap,
                            ..default()
                        },
                        TextColor(Color::srgb(0.80, 0.84, 0.90)),
                    ));
                    column
                        .spawn((
                            Node {
                                width: Val::Px(TRACK_W),
                                height: Val::Px(TRACK_H),
                                ..default()
                            },
                            BackgroundColor(TRACK_BG),
                        ))
                        .with_children(|track| {
                            track.spawn((
                                Node {
                                    position_type: PositionType::Absolute,
                                    left: Val::Px(0.0),
                                    top: Val::Px(0.0),
                                    width: Val::Px(0.0),
                                    height: Val::Px(TRACK_H),
                                    ..default()
                                },
                                BackgroundColor(FILL_SLACK),
                                MeterFill,
                            ));
                            // The ledger's median tightness. Drawn after the
                            // fill so it stays visible once the bar passes it,
                            // which is the whole point of having it -- and
                            // overhanging the track top and bottom, because at
                            // `record_gif.sh`'s downscale to 900 wide a 3 px
                            // notch inside a green bar is one pixel of cyan on
                            // green and cannot be seen at all.
                            track.spawn((
                                Node {
                                    position_type: PositionType::Absolute,
                                    left: Val::Px(median_ratio as f32 * TRACK_W),
                                    top: Val::Px(-6.0),
                                    width: Val::Px(3.0),
                                    height: Val::Px(TRACK_H + 12.0),
                                    ..default()
                                },
                                BackgroundColor(ANALYTIC),
                            ));
                        });
                    column.spawn((
                        Text::new(""),
                        TextFont {
                            font_size: FontSize::Px(18.0),
                            ..default()
                        },
                        TextLayout {
                            linebreak: bevy::text::LineBreak::NoWrap,
                            ..default()
                        },
                        TextColor(FILL_SLACK),
                        MeterValue,
                    ));
                });
        });
}

/// Size the bar, colour it, and write the number under it.
fn paint_meter(
    live: Res<Live>,
    flags: Res<ViewFlags>,
    mut fills: Query<(&mut Node, &mut BackgroundColor), With<MeterFill>>,
    mut values: Query<(&mut Text, &mut TextColor), With<MeterValue>>,
    mut panels: Query<&mut Visibility, With<HudPanel>>,
) {
    for mut visibility in &mut panels {
        *visibility = if flags.hud {
            Visibility::Inherited
        } else {
            Visibility::Hidden
        };
    }

    let row = &live.0.row;
    let bound = row.bound_deg();
    let (colour, fraction, text) = match (row.straddling_max_deg, row.ratio()) {
        (Some(worst), Some(ratio)) => (
            if ratio > 1.0 {
                FILL_BREACH
            } else if ratio >= TIGHT {
                FILL_GOOD
            } else {
                FILL_SLACK
            },
            ratio.clamp(0.0, 1.0) as f32,
            format!(
                "{worst:.6} deg measured   {bound:.6} deg bound   {ratio:.6}x{}",
                if ratio > 1.0 { "  BOUND BREACHED" } else { "" },
            ),
        ),
        _ => (
            TRACK_BG,
            0.0,
            format!(
                "no straddling vertex   {bound:.6} deg bound   crease unresolved \
                 at this cell"
            ),
        ),
    };

    for (mut node, mut background) in &mut fills {
        node.width = Val::Px(fraction * TRACK_W);
        background.0 = colour;
    }
    for (mut target, mut tint) in &mut values {
        if target.0 != text {
            target.0.clone_from(&text);
        }
        tint.0 = colour;
    }
}

// ─── the HUD ────────────────────────────────────────────────────────────────

/// The numbers are the demo.
///
/// Every line is kept inside 76 characters. At the harness's 13 px font that is
/// about 600 logical pixels, so nothing wraps at the 640x360 a smoke capture
/// uses -- and a wrapped line in a GIF reads as a bug.
fn report(live: Res<Live>, ledger: Res<Ledger>, replay: Res<Replay>, mut stats: ResMut<DemoStats>) {
    let row = &live.0.row;
    if row.samples == 0 {
        return;
    }
    let theta = row.dihedral_deg;
    let samples = row.samples;
    let drawn = {
        let total = live.0.straddles.len();
        if total == 0 {
            String::from("nothing to draw")
        } else {
            let stride = total.div_ceil(COMB_RAYS).max(1);
            format!("drawing every {stride} of them, worst in white")
        }
    };

    // The verdict is in the title as well as in the `check` line: a viewer
    // reading a GIF looks at the top-left and the bar, and a failed replay must
    // not be something you have to scroll fifteen lines to find.
    stats.title = format!(
        "E-311  seam normal bound - dihedral {theta:.1} deg  {samples}^3  \
         [1-3] ledger rows{}",
        if replay.exact {
            ""
        } else {
            "  - SELF-CHECK FAILED"
        },
    );
    stats.vertices = row.vertices;
    stats.triangles = row.triangles;
    stats.extract_ms = row.extract_ms;

    stats.extra = vec![
        format!(
            "fixture   A(r {R1:.3}) - B(r {R2:.3})  separation {:.6}  seam r {:.6}",
            row.separation, row.seam_radius,
        ),
        format!(
            "window    {WINDOW:.3} cube on one seam point, cell {:.6e} - a vertex",
            row.cell_size,
        ),
        String::from("          straddles only within the 6.06e-6 stencil of the branch"),
        String::from("          surface, and the canonical domain misses the crease by 1e-3"),
        String::new(),
        format!(
            "seam      {} cells crossed   {} vertices   {} degenerate normals",
            row.seam_cells, row.vertices, row.degenerate,
        ),
        format!("straddle  {} vertices - {drawn}", row.straddling),
        format!(
            "worst     {} deg   mean {} deg",
            row.straddling_max_deg
                .map_or_else(|| String::from("n/a"), |e| format!("{e:.6}")),
            row.straddling_mean_deg
                .map_or_else(|| String::from("n/a"), |e| format!("{e:.6}")),
        ),
        format!(
            "bound     (180 - {theta:.1})/2 = {:.6} deg   worst/bound {}",
            row.bound_deg(),
            row.ratio()
                .map_or_else(|| String::from("n/a"), |r| format!("{r:.6}x")),
        ),
        format!(
            "quiet     {} non-straddling   mean error {} deg",
            row.non_straddling,
            show(row.non_straddling_mean_deg),
        ),
        format!(
            "offset    median {} stencils from the branch surface",
            show(row.seam_offset_stencils),
        ),
        String::new(),
        ledger_line(&ledger, row),
        format!(
            "median    ledger tightness {} over {} measured rows",
            ledger.median_text, ledger.measured_rows,
        ),
        replay.line.clone(),
        format!(
            "          seam angle residual {:.3e} deg, p-56.csv column {}",
            replay.residual_deg,
            if replay.residual_matches {
                "MATCHES"
            } else {
                "DISAGREES"
            },
        ),
    ];
}

/// The live row against the committed one, when the sweep is standing on a row
/// the ledger measured.
///
/// The ramp lands exactly on all seven of P-56's dihedrals, so this fires seven
/// times per sweep at any resolution the CSV carries — and a disagreement is
/// reported rather than smoothed, because a live run that does not reproduce the
/// artefact is a finding about the artefact.
fn ledger_line(ledger: &Ledger, row: &Row) -> String {
    let Some(csv) = ledger.row(row.dihedral_deg, row.samples) else {
        return String::from("ledger    between measured dihedrals - no p-56.csv row here");
    };
    let mismatches = ledger.diff(csv, row);
    let columns = row.ledger_fields().len();
    format!(
        "ledger    p-56.csv worst {}  ratio {}  {}",
        ledger
            .field(csv, "straddling_max_error_deg")
            .unwrap_or("<absent>"),
        ledger
            .field(csv, "worst_over_bound_ratio")
            .unwrap_or("<absent>"),
        if mismatches.is_empty() {
            format!("all {columns} columns MATCH")
        } else {
            format!(
                "{} of {columns} DISAGREE ({})",
                mismatches.len(),
                mismatches[0].0
            )
        },
    )
}
