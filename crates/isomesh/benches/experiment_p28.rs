//! **P-28 — the medial identity against derived truths, not a second mollifier.**
//!
//! Ticket: R-030. Pre-registered in the commit before this one; the instrument
//! it replaces is recorded at V-46.
//!
//! ```bash
//! cargo bench --bench experiment_p28
//! ```
//!
//! Writes `docs/experiments/p-28.csv`.
//!
//! # The design in one line
//!
//! Measure `r_f = ρ·√(1 − ‖∇ρ‖²)` through the crate's real gradient path and
//! compare it against **derived truths**: exact closed forms where the field is
//! piecewise linear (slab, wedge, prism — the two- and three-point
//! closest-point cases), and an analytic zero where the true inscribed radius
//! vanishes (a capsule's off-axis interior), so the residual isolates the
//! discrete gradient and a wrong closed form has somewhere to go red.
//!
//! # The three registered clauses
//!
//! - **C1 (form).** On the piecewise-linear fixtures, per-axis voxel-step
//!   central differences of `min` of linear functions are computable exactly —
//!   `d_j(x ± h·eᵢ) = d_j(x) ± h·n̂_jᵢ`, then the min and the difference — so
//!   `r̃` is derived without touching [`central_difference`] or the field's
//!   sample path. Measured `r_f` must match within `1e-9` of the gap on 100%
//!   of medial-band samples.
//! - **C2 (curvature floor).** On the capsule's mid-axis band the field is an
//!   exact cylinder, every band sample has a *unique* closest point, and the
//!   true inscribed radius is **zero** — so `r_f` there is the formula's own
//!   noise: `‖CD‖² = 1 − O(h²/s²)` at radius `s`, √-amplified to
//!   `r_f ≈ ρ·h/s · c_or`. First order in `h`. The registered convergence is
//!   the 33³→129³ world-unit median ratio ≤ 0.35; ≥ 0.7 (an h-independent
//!   floor) is the falsifier that kills the three dependent mechanics.
//! - **C3 (clearance envelope).** Axis-aligned slabs of gap `W ∈ {3h, 6h,
//!   10h}` at 8 controlled sub-voxel phases: the band-max `r_f` must sit in
//!   `[√3/2·(W − h/2), W]`, from the profile `(W − φh)·√(1 − φ²)` minimised
//!   over `φ ∈ [0, ½]`.
//!
//! # The inversion, run before the verdict is read
//!
//! The wrong form `ρ·(1 − ‖∇ρ‖)` agrees with the right one at `‖∇ρ‖ ∈ {0, 1}`
//! and differs by up to `0.41·ρ` at `0.7` — the wedge's mid-range band hits
//! that by construction, and the wrong form must fail C1 on ≥ 30% of band
//! samples by more than `0.1` of the gap. If it does not go red, V-46 applies
//! to this instrument too.
//!
//! # Counted and measured against closed forms, not timed
//!
//! No timing A/B exists here, so M-197's interleaving rule does not apply, and
//! this note exists so its absence is not read as an oversight. Every fixture
//! constant is registered below; there is no randomness anywhere.

mod common;

use isomesh::Sdf;
use isomesh::normals::central_difference;

/// Domain half-width; grids are `n³` samples over `[-2, 2]`.
const DOMAIN: f64 = 2.0;
/// The three query resolutions (samples per axis).
const RESOLUTIONS: [u32; 3] = [33, 65, 129];
/// C1 tolerance scale: fixtures are built at gap scale `G = 0.4`, and the
/// registered float-dust tolerance is `1e-9·G`.
const GAP_SCALE: f64 = 0.4;
/// Samples closer to the boundary than `2h` are excluded everywhere — the
/// identity is about off-surface points, and `ρ ≥ 2h` is the registered
/// population cut.
const MIN_RHO_CELLS: f64 = 2.0;

/// A field that is the minimum of affine distance functions — air inside.
///
/// `d_j(x) = n_j · x + c_j`, `f = min_j d_j`. This one shape expresses the
/// tilted slab (two opposed planes), the wedge (two planes through an edge)
/// and the triangular prism (three planes around an axis), which are exactly
/// the fixtures whose voxel-step mollified truth is derivable in closed form.
struct MinLinear {
    planes: Vec<([f64; 3], f64)>,
}

impl MinLinear {
    fn distances(&self, p: [f64; 3]) -> Vec<f64> {
        self.planes
            .iter()
            .map(|(n, c)| n[0] * p[0] + n[1] * p[1] + n[2] * p[2] + c)
            .collect()
    }

    /// The derived mollified truth: per-axis central differences of a min of
    /// affine functions, evaluated from the plane data alone — no call into
    /// the field's sample path and none into [`central_difference`].
    fn rtilde(&self, p: [f64; 3], h: f64) -> f64 {
        let d = self.distances(p);
        let rho = d.iter().copied().fold(f64::INFINITY, f64::min);
        let mut norm2 = 0.0;
        for axis in 0..3 {
            let plus = self
                .planes
                .iter()
                .zip(&d)
                .map(|((n, _), dj)| dj + h * n[axis])
                .fold(f64::INFINITY, f64::min);
            let minus = self
                .planes
                .iter()
                .zip(&d)
                .map(|((n, _), dj)| dj - h * n[axis])
                .fold(f64::INFINITY, f64::min);
            let g = (plus - minus) / (2.0 * h);
            norm2 += g * g;
        }
        rho * (1.0 - norm2).max(0.0).sqrt()
    }
}

impl Sdf for MinLinear {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        self.distances(p).into_iter().fold(f64::INFINITY, f64::min)
    }
}

/// A capsule-shaped air pocket: `f = R − |p − closest point on segment|`.
/// Positive inside (air), the crate's sign convention.
struct AirCapsule {
    a: [f64; 3],
    b: [f64; 3],
    radius: f64,
}

impl AirCapsule {
    /// Axis parameter of the closest point on the segment, in `[0, 1]`.
    fn axis_t(&self, p: [f64; 3]) -> f64 {
        let ab = sub(self.b, self.a);
        let ap = sub(p, self.a);
        (dot(ap, ab) / dot(ab, ab)).clamp(0.0, 1.0)
    }

    /// Radial distance from the axis line (unclamped — valid mid-segment).
    fn radial(&self, p: [f64; 3]) -> f64 {
        let t = self.axis_t(p);
        let ab = sub(self.b, self.a);
        let foot = [
            self.a[0] + t * ab[0],
            self.a[1] + t * ab[1],
            self.a[2] + t * ab[2],
        ];
        len(sub(p, foot))
    }
}

impl Sdf for AirCapsule {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        self.radius - self.radial(p)
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn normalize(a: [f64; 3]) -> [f64; 3] {
    let l = len(a);
    assert!(l > 0.0, "registered fixture vector must be non-zero");
    [a[0] / l, a[1] / l, a[2] / l]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Gram–Schmidt: `v` made orthonormal to unit `e`.
fn ortho_to(v: [f64; 3], e: [f64; 3]) -> [f64; 3] {
    let d = dot(v, e);
    normalize([v[0] - d * e[0], v[1] - d * e[1], v[2] - d * e[2]])
}

/// The measured quantity: `r_f = ρ·√(clamp(1 − ‖∇ρ‖², 0, 1))` through the
/// crate's real gradient path. Returns `(r_f, clamped)`.
fn r_formula(field: &impl Sdf<Scalar = f64>, p: [f64; 3], h: f64) -> (f64, bool) {
    let rho = field.sample(p);
    let g = central_difference(field, p, h);
    let radicand = 1.0 - (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
    (rho * radicand.max(0.0).sqrt(), radicand < 0.0)
}

/// The deliberately wrong form for the inversion: `ρ·(1 − ‖∇ρ‖)`, clamped.
fn r_wrong(field: &impl Sdf<Scalar = f64>, p: [f64; 3], h: f64) -> f64 {
    let rho = field.sample(p);
    let g = central_difference(field, p, h);
    let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    rho * (1.0 - mag).max(0.0)
}

fn median(values: &mut [f64]) -> f64 {
    assert!(!values.is_empty(), "median of an empty band — reachability");
    values.sort_unstable_by(f64::total_cmp);
    values[values.len() / 2]
}

/// One registered row.
struct Row {
    fixture: String,
    n: u32,
    band: u64,
    within_pct: f64,
    median_world: f64,
    clearance_true_voxels: f64,
    clearance_est_voxels: f64,
    clamped: u64,
}

/// Scan a `MinLinear` fixture at `n³`: band = samples with derived `r̃ > 0`
/// and `ρ ≥ 2h`. Residual is `|r_f − r̃|`; `within` is the C1 float-dust
/// tolerance `1e-9·G`.
fn scan_min_linear(name: &str, field: &MinLinear, n: u32) -> Row {
    let h = 2.0 * DOMAIN / f64::from(n - 1);
    let tol = 1e-9 * GAP_SCALE;
    let mut residuals = Vec::new();
    let mut within = 0u64;
    let mut clamped = 0u64;
    let mut rt_max = 0.0f64;
    let mut rf_max = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let p = [
                    -DOMAIN + f64::from(i) * h,
                    -DOMAIN + f64::from(j) * h,
                    -DOMAIN + f64::from(k) * h,
                ];
                let rho = field.sample(p);
                if rho < MIN_RHO_CELLS * h {
                    continue;
                }
                let rt = field.rtilde(p, h);
                if rt <= 0.0 {
                    continue;
                }
                let (rf, was_clamped) = r_formula(field, p, h);
                if was_clamped {
                    clamped += 1;
                }
                let residual = (rf - rt).abs();
                if residual <= tol {
                    within += 1;
                }
                residuals.push(residual);
                rt_max = rt_max.max(rt);
                rf_max = rf_max.max(rf);
            }
        }
    }
    let band = residuals.len() as u64;
    assert!(band > 0, "{name} at {n}: empty medial band — reachability");
    Row {
        fixture: name.to_string(),
        n,
        band,
        within_pct: 100.0 * within as f64 / band as f64,
        median_world: median(&mut residuals),
        clearance_true_voxels: rt_max / h,
        clearance_est_voxels: rf_max / h,
        clamped,
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-28");
    common::experiment::run(prereg, |run| {
        // ---- registered fixture constants -------------------------------
        // Tilted slab: gap half-width 0.4, generic normal and offset.
        let slab_n = normalize([0.31, 0.92, 0.24]);
        let slab_x0 = [0.0137, -0.0731, 0.0413];
        let slab = MinLinear {
            planes: vec![
                (
                    [-slab_n[0], -slab_n[1], -slab_n[2]],
                    GAP_SCALE + dot(slab_n, slab_x0),
                ),
                (slab_n, GAP_SCALE - dot(slab_n, slab_x0)),
            ],
        };

        // Wedge: interior half-angle γ = 30°, so the on-bisector true
        // ‖∇ρ‖ = sin γ = 0.5 — the mid-range the inversion needs. Walls
        // through the edge line at x₀ along ê; normals sin γ·b̂ ± cos γ·t̂.
        let gamma = 30.0f64.to_radians();
        let e_hat = normalize([0.27, 0.4, 0.88]);
        let b_hat = ortho_to([0.9, -0.2, -0.3], e_hat);
        let t_hat = cross(e_hat, b_hat);
        let wedge_x0 = [0.0213, 0.0117, -0.0331];
        let mut wedge_planes = Vec::new();
        for sign in [1.0, -1.0] {
            let n = [
                gamma.sin() * b_hat[0] + sign * gamma.cos() * t_hat[0],
                gamma.sin() * b_hat[1] + sign * gamma.cos() * t_hat[1],
                gamma.sin() * b_hat[2] + sign * gamma.cos() * t_hat[2],
            ];
            wedge_planes.push((n, -dot(n, wedge_x0)));
        }
        let wedge = MinLinear {
            planes: wedge_planes,
        };

        // Prism: equilateral triangular air tube, inradius 0.45, tilted axis.
        // d_k = 0.45 − (x − x₀)·m̂_k with m̂_k the three outward normals.
        let a_hat = normalize([0.2, 0.86, 0.47]);
        let p_hat = ortho_to([1.0, 0.1, -0.4], a_hat);
        let q_hat = cross(a_hat, p_hat);
        let prism_x0 = [-0.0173, 0.0293, 0.0119];
        let inradius = 0.45;
        let mut prism_planes = Vec::new();
        for k in 0..3 {
            let ang = 2.0 * std::f64::consts::PI * f64::from(k) / 3.0;
            let m = [
                ang.cos() * p_hat[0] + ang.sin() * q_hat[0],
                ang.cos() * p_hat[1] + ang.sin() * q_hat[1],
                ang.cos() * p_hat[2] + ang.sin() * q_hat[2],
            ];
            prism_planes.push(([-m[0], -m[1], -m[2]], inradius + dot(m, prism_x0)));
        }
        let prism = MinLinear {
            planes: prism_planes,
        };

        // Capsule: air pocket, generic position. The C2 band is the exact
        // cylinder region: axis parameter in the middle [0.3, 0.7] and radius
        // s ∈ [0.3R, 0.7R], where every point has a unique closest boundary
        // point and the true inscribed radius is zero.
        let capsule = AirCapsule {
            a: [-0.83, -0.12, 0.07],
            b: [0.79, 0.15, -0.11],
            radius: 0.55,
        };

        let mut rows: Vec<Row> = Vec::new();

        // ---- C1: form, on the three piecewise-linear fixtures -----------
        for n in RESOLUTIONS {
            rows.push(scan_min_linear("slab_tilted", &slab, n));
            rows.push(scan_min_linear("wedge30", &wedge, n));
            rows.push(scan_min_linear("prism", &prism, n));
        }
        let c1_rows = rows.len();
        let c1_held = rows.iter().all(|r| r.within_pct >= 100.0);
        let prism_band: u64 = rows
            .iter()
            .filter(|r| r.fixture == "prism")
            .map(|r| r.band)
            .sum();
        assert!(
            prism_band > 0,
            "the |Π| = 3 case was not exercised — reachability"
        );

        // ---- C2: the curvature floor on the capsule ----------------------
        let mut c2_medians = Vec::new();
        for n in RESOLUTIONS {
            let h = 2.0 * DOMAIN / f64::from(n - 1);
            let mut floor = Vec::new();
            let mut clamped = 0u64;
            let mut over_cap = 0u64;
            for i in 0..n {
                for j in 0..n {
                    for k in 0..n {
                        let p = [
                            -DOMAIN + f64::from(i) * h,
                            -DOMAIN + f64::from(j) * h,
                            -DOMAIN + f64::from(k) * h,
                        ];
                        let t = capsule.axis_t(p);
                        let s = capsule.radial(p);
                        if !(0.3..=0.7).contains(&t)
                            || s < 0.3 * capsule.radius
                            || s > 0.7 * capsule.radius
                        {
                            continue;
                        }
                        let rho = capsule.sample(p);
                        if rho < MIN_RHO_CELLS * h {
                            continue;
                        }
                        let (rf, was_clamped) = r_formula(&capsule, p, h);
                        if was_clamped {
                            clamped += 1;
                        }
                        // Recorded (not gating): the derived pointwise cap.
                        if rf > 1.5 * rho * h / s {
                            over_cap += 1;
                        }
                        floor.push(rf);
                    }
                }
            }
            let band = floor.len() as u64;
            assert!(band > 0, "capsule band empty at {n} — reachability");
            let mut sorted = floor.clone();
            let med = median(&mut sorted);
            let rf_max = floor.iter().copied().fold(0.0f64, f64::max);
            c2_medians.push((n, med));
            rows.push(Row {
                fixture: "capsule_floor".to_string(),
                n,
                band,
                within_pct: 100.0 * (band - over_cap) as f64 / band as f64,
                median_world: med,
                clearance_true_voxels: 0.0,
                clearance_est_voxels: rf_max / h,
                clamped,
            });
        }
        let end_to_end = c2_medians[2].1 / c2_medians[0].1;

        // ---- C3: clearance envelope, axis-aligned slabs at 65³ ------------
        let n65 = 65u32;
        let h65 = 2.0 * DOMAIN / f64::from(n65 - 1);
        let mut c3_ok = 0u32;
        let mut c3_total = 0u32;
        for gap_cells in [3.0f64, 6.0, 10.0] {
            let w = gap_cells * h65;
            for phase16 in 0..8u32 {
                let phi = f64::from(phase16) / 16.0;
                let y0 = phi * h65;
                let n_axis = [0.0, 1.0, 0.0];
                let fixture = MinLinear {
                    planes: vec![([0.0, -1.0, 0.0], w + y0), (n_axis, w - y0)],
                };
                let row = scan_min_linear(&format!("slab_w{gap_cells}_p{phase16}"), &fixture, n65);
                let est = row.clearance_est_voxels * h65;
                let lo = 3.0f64.sqrt() / 2.0 * (w - h65 / 2.0);
                let inside = est >= lo && est <= w * (1.0 + 1e-12);
                c3_total += 1;
                if inside {
                    c3_ok += 1;
                }
                rows.push(Row {
                    fixture: format!("c3_w{gap_cells}_phi{phase16}of16"),
                    clearance_true_voxels: gap_cells,
                    ..row
                });
            }
        }

        // ---- Inversion: the wrong form must go red on C1 -----------------
        let h = 2.0 * DOMAIN / 64.0;
        let mut band = 0u64;
        let mut failed = 0u64;
        for i in 0..65u32 {
            for j in 0..65u32 {
                for k in 0..65u32 {
                    let p = [
                        -DOMAIN + f64::from(i) * h,
                        -DOMAIN + f64::from(j) * h,
                        -DOMAIN + f64::from(k) * h,
                    ];
                    let rho = wedge.sample(p);
                    if rho < MIN_RHO_CELLS * h {
                        continue;
                    }
                    let rt = wedge.rtilde(p, h);
                    if rt <= 0.0 {
                        continue;
                    }
                    band += 1;
                    if (r_wrong(&wedge, p, h) - rt).abs() > 0.1 * GAP_SCALE {
                        failed += 1;
                    }
                }
            }
        }
        let wrong_fail = failed as f64 / band as f64;
        assert!(
            wrong_fail >= 0.30,
            "inversion: the wrong form failed C1 on only {:.1}% of the wedge \
             band — the instrument has not been shown able to go red, and \
             V-46 applies to this design too",
            wrong_fail * 100.0
        );

        // ---- emit ---------------------------------------------------------
        println!(
            "{:>18} {:>4} {:>7} {:>12} {:>13} {:>10} {:>10} {:>8}",
            "fixture", "n", "band", "within%", "median", "true(vx)", "est(vx)", "clamped"
        );
        for r in &rows {
            println!(
                "{:>18} {:>4} {:>7} {:>12.4} {:>13.3e} {:>10.3} {:>10.3} {:>8}",
                r.fixture,
                r.n,
                r.band,
                r.within_pct,
                r.median_world,
                r.clearance_true_voxels,
                r.clearance_est_voxels,
                r.clamped
            );
            run.record(&[
                ("fixture", r.fixture.clone()),
                ("samples_per_axis", r.n.to_string()),
                ("band_samples", r.band.to_string()),
                ("within_tol_pct", format!("{:.4}", r.within_pct)),
                (
                    "band_median_residual_world",
                    format!("{:.6e}", r.median_world),
                ),
                (
                    "clearance_true_voxels",
                    format!("{:.4}", r.clearance_true_voxels),
                ),
                (
                    "clearance_est_voxels",
                    format!("{:.4}", r.clearance_est_voxels),
                ),
                ("clamped", r.clamped.to_string()),
            ]);
        }

        println!();
        println!(
            "C1 (form): {}/{c1_rows} piecewise-linear rows at 100% within 1e-9·G -- {}",
            rows.iter()
                .take(c1_rows)
                .filter(|r| r.within_pct >= 100.0)
                .count(),
            if c1_held { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C2 (curvature floor): medians {:.4e} / {:.4e} / {:.4e} world; 33->129 ratio {:.3} -- {} \
             (H says <= 0.35; falsified at >= 0.7)",
            c2_medians[0].1,
            c2_medians[1].1,
            c2_medians[2].1,
            end_to_end,
            if end_to_end <= 0.35 {
                "HELD"
            } else if end_to_end >= 0.7 {
                "FALSIFIED -- h-independent floor, the three mechanics die"
            } else {
                "UNDECIDED, loudly -- report as falsification-shaped"
            }
        );
        println!(
            "C3 (clearance envelope): {c3_ok}/{c3_total} rows inside [sqrt(3)/2*(W-h/2), W] -- {}",
            if c3_ok == c3_total {
                "HELD"
            } else {
                "FALSIFIED"
            }
        );
        println!(
            "inversion: wrong form fails C1 on {:.1}% of the wedge band (>= 30% required) -- RED \
             demonstrated",
            wrong_fail * 100.0
        );
    });
}
