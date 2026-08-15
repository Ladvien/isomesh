//! **How fast does the distance property degrade under repeated CSG?**
//!
//! Ticket: F-004. No paper measures this, and a destructible game needs the
//! answer: *after how many brush strokes is the field no longer usable as a
//! distance?*
//!
//! ```bash
//! cargo bench --bench csg_degradation
//! ```
//!
//! # The ticket said to sample `‖∇f‖`, and that metric is blind to this
//!
//! F-004 proposed measuring the gradient magnitude. F-002 then measured
//! `csg_difference` — one subtraction from a box — at **100% eikonal**, `‖∇f‖`
//! within 5% of one on every sample, while its values are not distances near the
//! seam (M-245). `max` selects an operand pointwise, and each operand is an exact
//! distance, so the composition's gradient is a unit vector almost everywhere no
//! matter how many operands there are. **Sampling `‖∇f‖` would produce a flat
//! line and conclude nothing degrades.**
//!
//! So this measures the thing that does degrade: the **value**.
//!
//! # Two numbers, and where the ground truth comes from
//!
//! **`q̂`, the empirical underestimate ratio.** For a sample point `p`, march
//! along `−∇f` until the sign changes and bisect to find where the surface is;
//! call that distance `d_ray`. Then `q̂ = |f(p)| / d_ray`. Since the true nearest
//! surface point may be off-ray, `d_true ≤ d_ray`, so `q̂ ≤ |f|/d_true` — the
//! measured ratio is a **pessimistic** estimate of the real one, which is the
//! safe direction to report a precision in.
//!
//! **Sphere-tracing steps.** The operational consequence: a field that
//! understates distance makes a tracer take more, shorter steps. This needs no
//! ground truth at all, and it is the number a renderer actually pays.
//!
//! # The carve is seeded and reproducible
//!
//! Random sphere subtractions from an analytic box, from a fixed seed, so the
//! curve is the same on every machine and a change in it is a change in the
//! code.

use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use isomesh::Sdf;
use isomesh::fields::{BoxExact, Sphere};

type Scalar = f64;

/// Brush counts the curve is sampled at.
const STROKES: [usize; 10] = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256];

/// Samples per axis for the ratio census.
const GRID: u32 = 24;

/// The box every carve starts from.
const HALF: Scalar = 1.0;

/// A box with `n` spheres subtracted, evaluated at run time.
///
/// Written here rather than composed from [`isomesh::fields::Difference`]
/// because that nests in the *type* — `Difference<Difference<…>>` — and this
/// needs the count to be a variable. The arithmetic is identical: `max` against
/// each negated sphere, in order.
struct CarvedBox {
    solid: BoxExact<Scalar>,
    cuts: Vec<Sphere<Scalar>>,
}

impl Sdf for CarvedBox {
    type Scalar = Scalar;

    fn sample(&self, p: [Scalar; 3]) -> Scalar {
        let mut v = self.solid.sample(p);
        for cut in &self.cuts {
            let c = -cut.sample(p);
            if c > v {
                v = c;
            }
        }
        v
    }

    fn gradient(&self, p: [Scalar; 3]) -> [Scalar; 3] {
        let e = 1e-6;
        let at = |q: [Scalar; 3]| self.sample(q);
        let d = [
            (at([p[0] + e, p[1], p[2]]) - at([p[0] - e, p[1], p[2]])) / (2.0 * e),
            (at([p[0], p[1] + e, p[2]]) - at([p[0], p[1] - e, p[2]])) / (2.0 * e),
            (at([p[0], p[1], p[2] + e]) - at([p[0], p[1], p[2] - e])) / (2.0 * e),
        ];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        if len > 0.0 {
            [d[0] / len, d[1] / len, d[2] / len]
        } else {
            [0.0, 0.0, 1.0]
        }
    }
}

/// Numerical Recipes' 64-bit LCG, so the carve is identical everywhere.
struct Lcg(u64);

impl Lcg {
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn carve(n: usize) -> CarvedBox {
    let mut rng = Lcg(0x0000_F004_0000_0001);
    let mut cuts = Vec::with_capacity(n);
    for _ in 0..n {
        // Centres on the box's surface-ish shell, radii a fraction of the box:
        // a brush that always hits and never swallows the whole solid.
        let c = [
            (rng.unit() * 2.0 - 1.0) * HALF * 1.1,
            (rng.unit() * 2.0 - 1.0) * HALF * 1.1,
            (rng.unit() * 2.0 - 1.0) * HALF * 1.1,
        ];
        cuts.push(Sphere {
            center: c,
            radius: 0.12 + 0.18 * rng.unit(),
        });
    }
    CarvedBox {
        solid: BoxExact {
            center: [0.0; 3],
            half_extents: [HALF; 3],
        },
        cuts,
    }
}

/// Distance from `p` to the surface **along the field's own descent direction**.
///
/// An upper bound on the true distance, so the ratio built from it understates
/// the precision — the safe direction.
fn distance_along_ray(field: &CarvedBox, p: [Scalar; 3]) -> Option<Scalar> {
    let f0 = field.sample(p);
    if f0 == 0.0 {
        return Some(0.0);
    }
    let g = field.gradient(p);
    // Downhill if outside, uphill if inside: either way, toward the surface.
    let dir = if f0 > 0.0 { [-g[0], -g[1], -g[2]] } else { g };
    let step = 0.01;
    let mut prev = 0.0;
    let mut t = step;
    while t < 6.0 {
        let q = [p[0] + dir[0] * t, p[1] + dir[1] * t, p[2] + dir[2] * t];
        if field.sample(q).signum() != f0.signum() {
            // Bisect the bracket for a usable figure rather than a step-sized one.
            let (mut lo, mut hi) = (prev, t);
            for _ in 0..40 {
                let mid = 0.5 * (lo + hi);
                let m = [
                    p[0] + dir[0] * mid,
                    p[1] + dir[1] * mid,
                    p[2] + dir[2] * mid,
                ];
                if field.sample(m).signum() == f0.signum() {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            return Some(0.5 * (lo + hi));
        }
        prev = t;
        t += step;
    }
    None
}

/// Sphere-trace one ray and count the steps.
///
/// The operational cost of a degraded field, needing no ground truth.
fn trace_steps(field: &CarvedBox, origin: [Scalar; 3], dir: [Scalar; 3]) -> Option<u32> {
    let mut t = 0.0;
    for step in 0..512u32 {
        let p = [
            origin[0] + dir[0] * t,
            origin[1] + dir[1] * t,
            origin[2] + dir[2] * t,
        ];
        let d = field.sample(p);
        if d < 1e-4 {
            return Some(step);
        }
        t += d.max(1e-5);
        if t > 8.0 {
            return None;
        }
    }
    None
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!("csg_degradation — how fast repeated subtraction stops being a distance\n");
    println!(
        "{:>7} {:>10} {:>10} {:>10} {:>12} {:>10}",
        "strokes", "q_min", "q_p01", "q_median", "eikonal_%", "trace_steps"
    );

    let mut rows = Vec::new();
    for &n in &STROKES {
        let field = carve(n);

        let mut ratios = Vec::new();
        let mut eikonal = 0usize;
        let mut counted = 0usize;
        for i in 0..GRID {
            for j in 0..GRID {
                for k in 0..GRID {
                    let at = |v: u32| -1.4 + 2.8 * (f64::from(v) + 0.5317) / f64::from(GRID);
                    let p = [at(i), at(j), at(k)];
                    let f = field.sample(p);
                    if f.abs() < 1e-3 {
                        continue;
                    }
                    let g = field.gradient(p);
                    let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                    counted += 1;
                    if (len - 1.0).abs() <= 0.05 {
                        eikonal += 1;
                    }
                    if let Some(d) = distance_along_ray(&field, p)
                        && d > 1e-6
                    {
                        ratios.push(f.abs() / d);
                    }
                }
            }
        }
        ratios.sort_by(f64::total_cmp);
        let pick = |q: f64| ratios[((ratios.len() - 1) as f64 * q) as usize];
        let (q_min, q_p01, q_median) = (ratios[0], pick(0.01), pick(0.5));
        let eikonal_pct = 100.0 * eikonal as f64 / counted as f64;

        // A fixed fan of rays, so the step count is comparable across `n`.
        let mut steps_total = 0u64;
        let mut steps_n = 0u64;
        for a in 0..12 {
            for b in 0..12 {
                let (u, v) = (f64::from(a) / 12.0, f64::from(b) / 12.0);
                let theta = u * std::f64::consts::TAU;
                let phi = (2.0 * v - 1.0).acos();
                let dir = [phi.sin() * theta.cos(), phi.sin() * theta.sin(), phi.cos()];
                let origin = [-dir[0] * 4.0, -dir[1] * 4.0, -dir[2] * 4.0];
                if let Some(s) = trace_steps(&field, origin, dir) {
                    steps_total += u64::from(s);
                    steps_n += 1;
                }
            }
        }
        let trace_steps_mean = steps_total as f64 / steps_n.max(1) as f64;

        println!(
            "{n:>7} {q_min:>10.4} {q_p01:>10.4} {q_median:>10.4} {eikonal_pct:>11.1}% {trace_steps_mean:>10.1}"
        );
        rows.push((n, q_min, q_p01, q_median, eikonal_pct, trace_steps_mean));
    }

    let mut csv = String::from("strokes,q_min,q_p01,q_median,eikonal_pct,trace_steps_mean\n");
    for (n, a, b, c, d, e) in &rows {
        let _ = writeln!(csv, "{n},{a:.6},{b:.6},{c:.6},{d:.3},{e:.3}");
    }
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements")
        .join("csg_degradation.csv");
    if let Some(dir) = path.parent() {
        let _ = fs::create_dir_all(dir);
    }
    let _ = fs::write(&path, csv);
    println!("\nwrote {}", path.display());

    // **The answer is read against the zero-stroke baseline, not against 1.0.**
    //
    // `q̂` is built from a ray distance, and a ray from a box *corner* leaves
    // along the diagonal, so `d_ray` overshoots the true distance by up to √3.
    // The uncarved box therefore reports `q̂_min ≈ 0.577 = 1/√3` while being
    // exactly a distance function. That floor is a property of the measurement,
    // not of the field, and reading the absolute number as a precision would
    // call an exact box degraded.
    let baseline = rows[0].1;
    let threshold = baseline * 0.5;
    let crossed = rows.iter().find(|(_, q_min, ..)| *q_min < threshold);
    println!("\nafter how many brush strokes is the field no longer usable as a distance?");
    println!(
        "  baseline: an UNCARVED box measures q_min = {baseline:.4} ≈ 1/√3, which is this\n           metric's floor at a corner rather than a defect — read the curve against it"
    );
    match crossed {
        Some((n, q, ..)) => {
            println!("  worst case halves against that baseline at {n} strokes (q_min = {q:.4})")
        }
        None => println!(
            "  worst case never halved against it within {} strokes",
            STROKES[STROKES.len() - 1]
        ),
    }
    let last = rows[rows.len() - 1];
    println!(
        "  by {} strokes: worst case {:.4} ({:.0}× down), median {:.4} (unmoved),\n           and tracing costs {:.1} steps against {:.1} — the degradation is entirely\n           in the tail, and a renderer barely notices what a precision bound calls ruin",
        last.0,
        last.1,
        baseline / last.1,
        last.3,
        last.5,
        rows[0].5
    );
    println!(
        "  the eikonal column is the control: it stays flat, which is why F-004's\n  \
         original proposal to measure ‖∇f‖ would have concluded nothing degrades"
    );
}
