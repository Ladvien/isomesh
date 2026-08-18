//! **P-33 — the arch golden values: a feasibility solver held to 0.1075 and 15.84°.**
//!
//! Ticket: R-034a. Pre-registered in the commit before this one; the solver
//! was corrected from FISTA to alternating projection *before any run*, on
//! convergence arithmetic recorded in the registration.
//!
//! ```bash
//! cargo bench --bench experiment_p33
//! ```
//!
//! Writes `docs/experiments/p-33.csv`.
//!
//! # The program (Whiting, Ochsendorf & Durand 2009, READ)
//!
//! Rigid blocks; a 3-D force at each interface vertex, expressed in the
//! interface frame as one axial and one in-plane component (2-D here, per
//! unit depth); equilibrium `A_eq·f = −w`, three rows per block; friction
//! `|f_t| ≤ μ·f_n` with μ = 0.7 (the paper's typical value); compression-only
//! `f_n ≥ 0`. A structure stands iff a force solution exists. The paper
//! validated on exactly the two numbers this harness bisects for: minimum
//! thickness/centerline-radius **0.1075** (Milankovitch 1907; their solver:
//! 0.10746 at 100 blocks) and critical ground tilt **15.84°** at t/r = 0.20
//! (Ochsendorf 2002).
//!
//! # The two decisions with committed reasons
//!
//! - **Weights at exact annular-sector centroids.** Area `Δθ·R·t`, centroid
//!   radius `(2 sin β / 3β)·(R_o³−R_i³)/(R_o²−R_i²)` — a centerline-lumped
//!   weight reproduces Heyman's 0.106, not Milankovitch's 0.1075, and the
//!   third decimal is the point. Torque rows are taken about each block's
//!   centroid with the weight applied there, so the centroid enters through
//!   every interface arm.
//! - **Alternating projection with the affine side exact.** One prefactored
//!   dense Cholesky of `A·Aᵀ` makes the equilibrium projection exact per
//!   iteration; the friction-cone projection is closed-form per vertex.
//!   Feasible probes converge into the intersection; infeasible probes stall
//!   at the minimal-distance pair. The decision reads the cone-side
//!   iterate's equilibrium residual per unit weight — **feasible below 1e-5,
//!   infeasible above 1e-4**, band asserted never hit — and both bisections
//!   stop with a classification margin (4e-4 in t/r, 0.02°) inside the
//!   registered tolerances.
//!
//! # Instrument checks, before any verdict
//!
//! Bracket ends must classify (0.08 infeasible, 0.16 feasible; 10° stands,
//! 20° falls) — both branches demonstrated reachable — and doubling gravity
//! must not move any classification: the program is scale-invariant or it is
//! wrong. No timing A/B exists; M-197 does not apply.

mod common;

/// Friction coefficient (the paper's typical value).
const MU: f64 = 0.7;
/// Alternating-projection iterations, fixed.
const ITERS: usize = 20_000;
/// Decision thresholds on ‖A·z + w‖ / W.
const FEASIBLE_BELOW: f64 = 1e-5;
const INFEASIBLE_ABOVE: f64 = 1e-4;
/// Centerline radius.
const R: f64 = 1.0;

struct Problem {
    /// Dense equilibrium matrix, rows-major: 3·blocks × 4·(blocks+1).
    a: Vec<f64>,
    rows: usize,
    cols: usize,
    /// Right-hand side: −(weight terms), per row.
    b: Vec<f64>,
    /// Total weight, the residual scale.
    total_w: f64,
}

/// Build the semicircular arch: `n` radial blocks, thickness ratio `t_over_r`,
/// gravity tilted by `tilt` radians (the ground-tilt experiment rotates
/// gravity in the arch frame), gravity magnitude `g`.
#[allow(
    clippy::needless_range_loop,
    reason = "dense matrix kernels index two arrays per loop; iterator forms obscure the algebra"
)]
fn arch(n: usize, t_over_r: f64, tilt: f64, g: f64) -> Problem {
    let t = t_over_r * R;
    let ri = R - t / 2.0;
    let ro = R + t / 2.0;
    let dth = std::f64::consts::PI / n as f64;
    let rows = 3 * n;
    let cols = 4 * (n + 1);
    let mut a = vec![0.0f64; rows * cols];
    let mut b = vec![0.0f64; rows];
    let gv = [g * tilt.sin(), -g * tilt.cos()];

    // Interface k sits at θ_k = k·Δθ; its two vertices at radii ri and ro;
    // normal along +θ, in-plane tangent radial.
    let iface = |k: usize| -> ([f64; 2], [f64; 2], [[f64; 2]; 2]) {
        let th = dth * k as f64;
        let nrm = [-th.sin(), th.cos()];
        let tan = [th.cos(), th.sin()];
        let vlo = [ri * th.cos(), ri * th.sin()];
        let vhi = [ro * th.cos(), ro * th.sin()];
        (nrm, tan, [vlo, vhi])
    };

    let mut total_w = 0.0;
    for j in 0..n {
        // Annular-sector centroid, exact.
        let beta = dth / 2.0;
        let rbar = (2.0 * beta.sin() / (3.0 * beta)) * (ro.powi(3) - ri.powi(3))
            / (ro.powi(2) - ri.powi(2));
        let phi = dth * (j as f64 + 0.5);
        let c = [rbar * phi.cos(), rbar * phi.sin()];
        let w = dth * R * t * g;
        total_w += w;

        let fx = 3 * j;
        let fy = 3 * j + 1;
        let tq = 3 * j + 2;
        b[fx] = -w * gv[0] / g;
        b[fy] = -w * gv[1] / g;
        // (weight acts at the centroid; its torque about the centroid is 0.)

        // Interface j acts on block j with +F, interface j+1 with −F.
        for (k, sign) in [(j, 1.0f64), (j + 1, -1.0f64)] {
            let (nrm, tan, verts) = iface(k);
            for (v, p) in verts.iter().enumerate() {
                let col_n = 4 * k + 2 * v;
                let col_t = col_n + 1;
                a[fx * cols + col_n] += sign * nrm[0];
                a[fx * cols + col_t] += sign * tan[0];
                a[fy * cols + col_n] += sign * nrm[1];
                a[fy * cols + col_t] += sign * tan[1];
                let arm = [p[0] - c[0], p[1] - c[1]];
                a[tq * cols + col_n] += sign * (arm[0] * nrm[1] - arm[1] * nrm[0]);
                a[tq * cols + col_t] += sign * (arm[0] * tan[1] - arm[1] * tan[0]);
            }
        }
    }
    Problem {
        a,
        rows,
        cols,
        b,
        total_w,
    }
}

/// Dense Cholesky of A·Aᵀ, lower factor in place. Aborts loudly if rank
/// drops — that is a formulation failure, not a case to smooth over.
struct AffineProjector {
    rows: usize,
    cols: usize,
    a: Vec<f64>,
    l: Vec<f64>,
}

#[allow(
    clippy::needless_range_loop,
    reason = "dense matrix kernels index two arrays per loop; iterator forms obscure the algebra"
)]
impl AffineProjector {
    fn new(p: &Problem) -> Self {
        let (rows, cols) = (p.rows, p.cols);
        let mut m = vec![0.0f64; rows * rows];
        for i in 0..rows {
            for j in 0..=i {
                let mut s = 0.0;
                for k in 0..cols {
                    s += p.a[i * cols + k] * p.a[j * cols + k];
                }
                m[i * rows + j] = s;
                m[j * rows + i] = s;
            }
        }
        for j in 0..rows {
            for k in 0..j {
                let ljk = m[j * rows + k];
                if ljk == 0.0 {
                    continue;
                }
                for i in j..rows {
                    m[i * rows + j] -= m[i * rows + k] * ljk;
                }
            }
            let d = m[j * rows + j];
            assert!(d > 1e-12, "A·Aᵀ lost rank at row {j} — formulation failure");
            let root = d.sqrt();
            for i in j..rows {
                m[i * rows + j] /= root;
            }
        }
        AffineProjector {
            rows,
            cols,
            a: p.a.clone(),
            l: m,
        }
    }

    /// f ← f − Aᵀ·(A·Aᵀ)⁻¹·(A·f − b): exact projection onto {A·f = b}.
    fn project(&self, f: &mut [f64], b: &[f64], scratch: &mut [f64]) {
        let (rows, cols) = (self.rows, self.cols);
        for i in 0..rows {
            let mut s = -b[i];
            for k in 0..cols {
                s += self.a[i * cols + k] * f[k];
            }
            scratch[i] = s;
        }
        for i in 0..rows {
            let mut s = scratch[i];
            for k in 0..i {
                s -= self.l[i * rows + k] * scratch[k];
            }
            scratch[i] = s / self.l[i * rows + i];
        }
        for i in (0..rows).rev() {
            let mut s = scratch[i];
            for k in i + 1..rows {
                s -= self.l[k * rows + i] * scratch[k];
            }
            scratch[i] = s / self.l[i * rows + i];
        }
        for k in 0..cols {
            let mut s = 0.0;
            for i in 0..rows {
                s += self.a[i * cols + k] * scratch[i];
            }
            f[k] -= s;
        }
    }
}

/// Closed-form projection of each vertex's (f_n, f_t) onto the friction cone.
fn project_cone(f: &mut [f64]) {
    for pair in f.chunks_exact_mut(2) {
        let n = pair[0];
        let t = pair[1];
        if t.abs() <= MU * n {
            continue;
        }
        if MU * t.abs() <= -n {
            pair[0] = 0.0;
            pair[1] = 0.0;
        } else {
            let nn = (n + MU * t.abs()) / (1.0 + MU * MU);
            pair[0] = nn;
            pair[1] = t.signum() * MU * nn;
        }
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Verdict {
    Feasible,
    Infeasible,
}

/// Run the alternating projection and classify. Returns the verdict and the
/// cone-side residual per unit weight.
#[allow(
    clippy::needless_range_loop,
    reason = "dense matrix kernels index two arrays per loop; iterator forms obscure the algebra"
)]
fn classify(p: &Problem) -> (Verdict, f64) {
    let proj = AffineProjector::new(p);
    let mut f = vec![0.0f64; p.cols];
    let mut scratch = vec![0.0f64; p.rows];
    for _ in 0..ITERS {
        proj.project(&mut f, &p.b, &mut scratch);
        project_cone(&mut f);
    }
    let mut e = 0.0f64;
    for i in 0..p.rows {
        let mut s = -p.b[i];
        for k in 0..p.cols {
            s += p.a[i * p.cols + k] * f[k];
        }
        e += s * s;
    }
    let rel = e.sqrt() / p.total_w;
    let verdict = if rel < FEASIBLE_BELOW {
        Verdict::Feasible
    } else if rel > INFEASIBLE_ABOVE {
        Verdict::Infeasible
    } else {
        panic!(
            "residual {rel:.3e} landed in the undecided band ({FEASIBLE_BELOW:e}, \
             {INFEASIBLE_ABOVE:e}] — the instrument refuses to guess"
        );
    };
    (verdict, rel)
}

/// Bisect thickness/radius for the minimum standing thickness at `n` blocks.
/// Returns (threshold, last feasible residual, last infeasible residual).
fn min_thickness(n: usize) -> (f64, f64, f64) {
    let (mut lo, mut hi) = (0.08f64, 0.16f64);
    let (vlo, rlo) = classify(&arch(n, lo, 0.0, 1.0));
    let (vhi, rhi) = classify(&arch(n, hi, 0.0, 1.0));
    assert!(
        vlo == Verdict::Infeasible && vhi == Verdict::Feasible,
        "bracket ends failed to classify as registered (0.08 infeasible, 0.16 feasible)"
    );
    let (mut r_feas, mut r_inf) = (rhi, rlo);
    while hi - lo > 4e-4 {
        let mid = 0.5 * (lo + hi);
        let (v, r) = classify(&arch(n, mid, 0.0, 1.0));
        match v {
            Verdict::Feasible => {
                hi = mid;
                r_feas = r;
            }
            Verdict::Infeasible => {
                lo = mid;
                r_inf = r;
            }
        }
    }
    (0.5 * (lo + hi), r_feas, r_inf)
}

/// Bisect the ground tilt (gravity rotation) at t/r = 0.20, 100 blocks.
fn critical_tilt(n: usize) -> (f64, f64, f64) {
    let (mut lo, mut hi) = (10.0f64.to_radians(), 20.0f64.to_radians());
    let (vlo, rlo) = classify(&arch(n, 0.20, lo, 1.0));
    let (vhi, rhi) = classify(&arch(n, 0.20, hi, 1.0));
    assert!(
        vlo == Verdict::Feasible && vhi == Verdict::Infeasible,
        "tilt bracket ends failed to classify as registered (10° stands, 20° falls)"
    );
    let (mut r_feas, mut r_inf) = (rlo, rhi);
    while hi - lo > 0.02f64.to_radians() {
        let mid = 0.5 * (lo + hi);
        let (v, r) = classify(&arch(n, 0.20, mid, 1.0));
        match v {
            Verdict::Feasible => {
                lo = mid;
                r_feas = r;
            }
            Verdict::Infeasible => {
                hi = mid;
                r_inf = r;
            }
        }
    }
    (0.5 * (lo + hi).to_degrees(), r_feas, r_inf)
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-33");
    common::experiment::run(prereg, |run| {
        // Scale invariance: doubling gravity must not move any verdict.
        for tr in [0.10f64, 0.12] {
            let (v1, _) = classify(&arch(100, tr, 0.0, 1.0));
            let (v2, _) = classify(&arch(100, tr, 0.0, 2.0));
            assert!(
                v1 == v2,
                "gravity doubling moved the verdict at t/r = {tr}: {v1:?} vs {v2:?}"
            );
        }
        println!("scale invariance: verdicts unmoved under doubled gravity at t/r = 0.10, 0.12");

        type Row = (String, usize, f64, f64, f64, bool, f64, f64);
        let mut rows: Vec<Row> = Vec::new();

        let (thr100, rf, ri) = min_thickness(100);
        let err = (thr100 - 0.1075).abs();
        rows.push((
            "milankovitch_100".into(),
            100,
            thr100,
            0.1075,
            err,
            err <= 0.0010,
            rf,
            ri,
        ));

        let (tilt, tf, ti) = critical_tilt(100);
        let terr = (tilt - 15.84).abs();
        rows.push((
            "tilt_100".into(),
            100,
            tilt,
            15.84,
            terr,
            terr <= 0.05,
            tf,
            ti,
        ));

        for n in [25usize, 50, 200] {
            let (thr, rf, ri) = classify_threshold_for_coarsening(n);
            // Recorded, not registered: the paper's direction — coarser
            // over-estimates stability, so thr(coarse) should sit below
            // thr(100).
            rows.push((
                format!("coarsen_{n}"),
                n,
                thr,
                thr100,
                thr - thr100,
                thr < thr100 || n > 100,
                rf,
                ri,
            ));
        }

        println!(
            "\n{:>18} {:>7} {:>10} {:>9} {:>10} {:>7} {:>10} {:>10}",
            "test", "blocks", "value", "target", "err", "within", "res_feas", "res_inf"
        );
        for (test, blocks, value, target, err, within, rf, ri) in &rows {
            println!(
                "{:>18} {:>7} {:>10.5} {:>9.4} {:>10.5} {:>7} {:>10.2e} {:>10.2e}",
                test, blocks, value, target, err, within, rf, ri
            );
            run.record(&[
                ("test", test.clone()),
                ("blocks", blocks.to_string()),
                ("value", format!("{value:.5}")),
                ("target", format!("{target:.5}")),
                ("abs_error", format!("{err:.5}")),
                ("within_tolerance", within.to_string()),
                ("residual_feasible", format!("{rf:.3e}")),
                ("residual_infeasible", format!("{ri:.3e}")),
            ]);
        }

        println!();
        let gold_thickness = rows[0].5;
        let gold_tilt = rows[1].5;
        println!(
            "Milankovitch: {:.5} vs 0.1075 (Whiting: 0.10746) -- {}",
            rows[0].2,
            if gold_thickness { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "Ochsendorf tilt: {:.3}° vs 15.84° -- {}",
            rows[1].2,
            if gold_tilt { "HELD" } else { "FALSIFIED" }
        );
        let coarser_below = rows[2].2 < thr100 && rows[3].2 < thr100;
        println!(
            "coarsening (recorded): thr(25) = {:.5}, thr(50) = {:.5}, thr(200) = {:.5} against \
             thr(100) = {:.5} -- coarser-sits-below {}",
            rows[2].2,
            rows[3].2,
            rows[4].2,
            thr100,
            if coarser_below {
                "as the paper warns"
            } else {
                "NOT observed — the paper's direction did not reproduce here"
            }
        );
    });
}

/// The same thickness bisection, named for the coarsening sweep.
fn classify_threshold_for_coarsening(n: usize) -> (f64, f64, f64) {
    min_thickness(n)
}
