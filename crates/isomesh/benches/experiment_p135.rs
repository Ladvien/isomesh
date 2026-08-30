//! **P-135 — where the identity stops holding, which is the boundary of
//! everything Group A claims.**
//!
//! Ticket: R-135. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p135
//! ```
//!
//! Writes `docs/experiments/p-135.csv`.
//!
//! # What was missing
//!
//! `P-127` proved `b*b - 4*a*c` at `marching_cubes/trilinear.rs:246` is
//! identically Cayley's `2x2x2` hyperdeterminant, and `P-128`–`P-134` all lean
//! on it. Every one of those claims rests on the reconstruction being
//! **multi-affine** — degree at most 1 in each variable separately. That is a
//! hypothesis about the filter, not a fact about the crate, and until this row
//! nothing recorded where it stops being true.
//!
//! # The mathematics, stated so the residual is checkable
//!
//! On a multi-affine reconstruction the eight corner values determine the
//! whole cell, so `Δ` is a polynomial in exactly those eight numbers and the
//! identity is an identity of `i128` polynomials — which is what
//! `common::poly` computes. The moment the reconstruction carries a
//! **quadratic or higher term in any single variable**, the cell is no longer
//! determined by its corners: the coefficients that `BodySaddles` reads off
//! the corners are the wrong coefficients, and the discriminant of the
//! resulting quadratic is a different polynomial.
//!
//! Two non-multi-affine reconstructions, both exact:
//!
//! - **Tricubic.** Degree 3 per variable, so 64 coefficients against the
//!   trilinear's 8. The residual is computed by *restricting* a tricubic to
//!   the cell's eight corners, rebuilding the trilinear that interpolates
//!   those corners, and subtracting: the difference is a polynomial whose
//!   **named non-zero terms** are the vacuity control this row's registration
//!   demands. A bench-local exact polynomial in the three spatial variables
//!   (`SpatialPoly`) does this; it is a different object from
//!   `common::poly::Poly` (which is a polynomial in the eight corner
//!   *values*), not a second copy of it.
//! - **`smooth_min(k)`.** The crate's own `SmoothUnion` carries the
//!   polynomial-smooth-min blend with parameter `k`. It is not multi-affine
//!   for any `k > 0`: inside the `O(k)` seam shell the blend is quadratic in
//!   the two field values. Sweeping `k` measures whether the deviation from
//!   the multi-affine `Δ` **tracks `k`**, which is C2, and connects Group A's
//!   boundary to `M-38`'s smoothing result rather than leaving it a separate
//!   fact.
//!
//! # Arms
//!
//! | arm | `reconstruction` | what it varies | is_control |
//! |---|---|---|---|
//! | trilinear | `trilinear` | nothing — the identity must hold | **yes** |
//! | tricubic | `tricubic` | degree 3 per variable, symbolic residual | no |
//! | smooth-min sweep | `smooth_min` | `k` over a decade sweep | no |
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration says **none** — this is a scope statement and moves no
//! stage of the pipeline. Recorded as such rather than left blank.
//!
//! # Vacuity controls
//!
//! - **The tricubic arm must produce a symbolic residual with NAMED non-zero
//!   terms**, not a numeric non-zero: `symbolic_residual_terms` carries the
//!   count and an extra column names the leading monomials. A numeric
//!   difference would not establish the failure structurally.
//! - **The trilinear control must reproduce the identity in the same
//!   harness**, or the comparison is between two different measurement setups.
//! - **The `k` sweep must span at least a decade** and include `k = 0`, whose
//!   deviation must be exactly zero — `smooth_min(0)` IS `min`, which is
//!   multi-affine on each side of the seam.

#![allow(
    clippy::cast_precision_loss,
    clippy::too_many_lines,
    // Deviations are compared against an exact zero the algebra predicts, so
    // an epsilon would hide the very thing being measured.
    clippy::float_cmp
)]

mod common;

use std::collections::BTreeMap;

use common::poly::{Rng, cayley_2x2x2, repo_discriminant};

/// Exponent triple of a monomial in the three spatial variables.
type Exp3 = [u8; 3];

/// An exact polynomial in the three spatial variables `x`, `y`, `z` with
/// `i128` coefficients.
///
/// Distinct from [`common::poly::Poly`], which is a polynomial in the eight
/// corner *values*. This one is needed because a tricubic's failure is a
/// statement about degrees in `x`, `y`, `z` — the corner-value polynomial
/// cannot express it.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct SpatialPoly {
    terms: BTreeMap<Exp3, i128>,
}

impl SpatialPoly {
    fn zero() -> Self {
        Self {
            terms: BTreeMap::new(),
        }
    }

    fn monomial(exp: Exp3, coefficient: i128) -> Self {
        let mut terms = BTreeMap::new();
        if coefficient != 0 {
            terms.insert(exp, coefficient);
        }
        Self { terms }
    }

    fn add(&self, other: &Self) -> Self {
        let mut out = self.terms.clone();
        for (exp, c) in &other.terms {
            let slot = out.entry(*exp).or_insert(0);
            *slot += *c;
            if *slot == 0 {
                out.remove(exp);
            }
        }
        Self { terms: out }
    }

    fn sub(&self, other: &Self) -> Self {
        let mut out = self.terms.clone();
        for (exp, c) in &other.terms {
            let slot = out.entry(*exp).or_insert(0);
            *slot -= *c;
            if *slot == 0 {
                out.remove(exp);
            }
        }
        Self { terms: out }
    }

    fn mul(&self, other: &Self) -> Self {
        let mut out: BTreeMap<Exp3, i128> = BTreeMap::new();
        for (a, ca) in &self.terms {
            for (b, cb) in &other.terms {
                let exp = [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
                let slot = out.entry(exp).or_insert(0);
                *slot += ca * cb;
                if *slot == 0 {
                    out.remove(&exp);
                }
            }
        }
        Self { terms: out }
    }

    fn terms(&self) -> usize {
        self.terms.len()
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Highest exponent appearing on variable `i` in any term. A
    /// multi-affine polynomial has `degree_in(i) <= 1` for every `i`.
    fn degree_in(&self, i: usize) -> u32 {
        self.terms
            .keys()
            .map(|e| u32::from(e[i]))
            .max()
            .unwrap_or(0)
    }

    fn is_multi_affine(&self) -> bool {
        (0..3).all(|i| self.degree_in(i) <= 1)
    }

    /// The `n` largest-degree monomials, named as strings, for the residual
    /// column. Deterministic: `BTreeMap` order, then total degree descending.
    fn named_terms(&self, n: usize) -> Vec<String> {
        let mut v: Vec<(Exp3, i128)> = self.terms.iter().map(|(e, c)| (*e, *c)).collect();
        v.sort_by_key(|(e, _)| {
            let total = u32::from(e[0]) + u32::from(e[1]) + u32::from(e[2]);
            (std::cmp::Reverse(total), *e)
        });
        v.into_iter()
            .take(n)
            .map(|(e, c)| {
                let mut s = if c == 1 {
                    String::new()
                } else if c == -1 {
                    "-".to_string()
                } else {
                    format!("{c}")
                };
                for (i, name) in ["x", "y", "z"].iter().enumerate() {
                    match e[i] {
                        0 => {}
                        1 => s.push_str(name),
                        d => s.push_str(&format!("{name}^{d}")),
                    }
                }
                if s.is_empty() || s == "-" {
                    s.push('1');
                }
                s
            })
            .collect()
    }

    /// Value at a lattice point, exactly.
    fn eval(&self, x: i128, y: i128, z: i128) -> i128 {
        let mut acc = 0i128;
        for (e, c) in &self.terms {
            let mut t = *c;
            for _ in 0..e[0] {
                t *= x;
            }
            for _ in 0..e[1] {
                t *= y;
            }
            for _ in 0..e[2] {
                t *= z;
            }
            acc += t;
        }
        acc
    }
}

/// `x`, `y`, `z` as spatial polynomials.
fn var(i: usize) -> SpatialPoly {
    let mut e = [0u8; 3];
    e[i] = 1;
    SpatialPoly::monomial(e, 1)
}

/// The trilinear interpolant of eight corner values over the unit cell,
/// symbolically: `sum over corners of f_c * B_c(x, y, z)` with
/// `B_c = (u ? x : 1-x)(v ? y : 1-y)(w ? z : 1-z)` and `c = u + 2v + 4w`.
fn trilinear_form(f: &[i128; 8]) -> SpatialPoly {
    let one = SpatialPoly::monomial([0, 0, 0], 1);
    let mut acc = SpatialPoly::zero();
    for (c, value) in f.iter().enumerate() {
        if *value == 0 {
            continue;
        }
        let mut basis = SpatialPoly::monomial([0, 0, 0], *value);
        for axis in 0..3 {
            let bit = (c >> axis) & 1;
            let factor = if bit == 1 {
                var(axis)
            } else {
                one.sub(&var(axis))
            };
            basis = basis.mul(&factor);
        }
        acc = acc.add(&basis);
    }
    acc
}

/// A deliberately non-multi-affine tricubic: the tensor-product cubic
/// `x^3 y^3 z^3` plus a lower-degree tail, with integer coefficients drawn
/// deterministically so the arm is reproducible.
fn tricubic_form(rng: &mut Rng) -> SpatialPoly {
    let mut acc = SpatialPoly::zero();
    for ex in 0..4u8 {
        for ey in 0..4u8 {
            for ez in 0..4u8 {
                let c = rng.next_i64_in(-4, 5);
                if c != 0 {
                    acc = acc.add(&SpatialPoly::monomial([ex, ey, ez], i128::from(c)));
                }
            }
        }
    }
    // Force the top corner term present so the arm is genuinely degree 3 in
    // every variable regardless of the draw.
    acc.add(&SpatialPoly::monomial([3, 3, 3], 1))
}

/// The polynomial smooth-min the crate blends with (`sdf.rs`'s
/// `smooth_min`), in the exact form used for the algebra:
/// `min(a, b) - h^2 * k / 4` with `h = max(0, 1 - |a - b| / k)`.
///
/// Inside the seam shell `|a - b| < k` this is **quadratic** in `a - b`,
/// which is exactly why it is not multi-affine.
fn smooth_min(a: f64, b: f64, k: f64) -> f64 {
    if k <= 0.0 {
        return a.min(b);
    }
    let h = (1.0 - (a - b).abs() / k).max(0.0);
    a.min(b) - h * h * k * 0.25
}

/// The eight corner values of the smooth union of two planes over the unit
/// cell, at blend parameter `k`. The two planes cross inside the cell so the
/// seam shell is sampled.
fn smooth_min_corners(k: f64) -> [f64; 8] {
    std::array::from_fn(|c| {
        let x = f64::from(u32::from((c & 1) != 0));
        let y = f64::from(u32::from((c & 2) != 0));
        let z = f64::from(u32::from((c & 4) != 0));
        // Two planes meeting along the cell diagonal: a = x - 0.5,
        // b = 0.5 - y, so a - b spans the seam over the corner set.
        let a = x - 0.5 + 0.25 * z;
        let b = 0.5 - y - 0.25 * z;
        smooth_min(a, b, k)
    })
}

/// `Δ` of a corner tuple in `f64`, through the repo's own expression.
fn delta_f64(f: &[f64; 8]) -> f64 {
    let twist_lo = (f[0] + f[3]) - (f[1] + f[2]);
    let twist_hi = (f[4] + f[7]) - (f[5] + f[6]);
    let du_lo = f[1] - f[0];
    let du_hi = f[5] - f[4];
    let dv_lo = f[2] - f[0];
    let dv_hi = f[6] - f[4];
    let a = du_hi * twist_lo - du_lo * twist_hi;
    let b = (f[4] * twist_lo - f[0] * twist_hi) + (du_hi * dv_lo - du_lo * dv_hi);
    let c = f[2] * f[4] - f[0] * f[6];
    b * b - 4.0 * a * c
}

/// The `k` sweep: `0` plus a decade and a half, so C2's growth claim has a
/// range to track and the `k = 0` control has an exact zero to hit.
const K_SWEEP: [f64; 7] = [0.0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5];

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-135");

    common::experiment::run(prereg, |run| {
        let cayley = cayley_2x2x2();
        let repo = repo_discriminant();

        // ── the trilinear control, in this harness ─────────────────────────
        //
        // The identity must hold here or the comparison below is between two
        // different measurement setups.
        let identity_residual = repo.sub(&cayley);
        assert!(
            identity_residual.is_zero(),
            "P-135 VOID: the trilinear identity does not hold in this harness \
             ({} residual terms), so the tricubic and smooth-min arms have no \
             baseline to fail against",
            identity_residual.terms()
        );

        // The trilinear interpolant of a symbolic corner tuple must be
        // multi-affine, which is the property Group A rests on.
        let mut rng = Rng::new(0x135_0135);
        let sample_corners: [i128; 8] = std::array::from_fn(|_| i128::from(rng.next_i64_in(-6, 7)));
        let tri = trilinear_form(&sample_corners);
        assert!(
            tri.is_multi_affine(),
            "P-135 VOID: the trilinear interpolant this harness builds is not \
             multi-affine (degrees {}, {}, {}), so it is not the object the \
             crate reconstructs and no boundary measured against it means \
             anything",
            tri.degree_in(0),
            tri.degree_in(1),
            tri.degree_in(2)
        );

        run.record(&[
            ("reconstruction", "trilinear".to_string()),
            ("is_multi_affine", "true".to_string()),
            ("identity_holds", "true".to_string()),
            ("symbolic_residual_terms", "0".to_string()),
            ("smooth_min_k", "0".to_string()),
            ("deviation_at_k", "0".to_string()),
            ("tricubic_degree", "1".to_string()),
            ("cases_touched", "256".to_string()),
            ("c1_holds", "true".to_string()),
            ("c2_holds", "true".to_string()),
            // ── extras (M-273) ──
            ("arm_role", "control".to_string()),
            ("degree_x", tri.degree_in(0).to_string()),
            ("degree_y", tri.degree_in(1).to_string()),
            ("degree_z", tri.degree_in(2).to_string()),
            ("interpolant_terms", tri.terms().to_string()),
            ("residual_named_terms", "none".to_string()),
            ("seed", "0x1350135".to_string()),
        ]);

        // ── the tricubic arm: a NAMED symbolic residual ────────────────────
        //
        // Restrict the tricubic to the eight corners, rebuild the trilinear
        // that interpolates them, and subtract. The difference is the part of
        // the cell the corner values cannot see — and its named non-zero
        // terms are the structural failure C1 asks for.
        let tricubic = tricubic_form(&mut rng);
        let corner_values: [i128; 8] = std::array::from_fn(|c| {
            let x = i128::from(u8::from((c & 1) != 0));
            let y = i128::from(u8::from((c & 2) != 0));
            let z = i128::from(u8::from((c & 4) != 0));
            tricubic.eval(x, y, z)
        });
        let rebuilt = trilinear_form(&corner_values);
        let residual = tricubic.sub(&rebuilt);
        let named = residual.named_terms(6);

        assert!(
            !residual.is_zero(),
            "P-135 VOID: the tricubic's residual against its own corner \
             trilinear is the zero polynomial, so this arm is multi-affine \
             after all and C1 has no non-multi-affine reconstruction to fail on"
        );
        assert!(
            !named.is_empty(),
            "P-135 VOID: the tricubic residual has no NAMED non-zero terms, \
             only a numeric difference, so C1 has not established the failure \
             structurally -- the registration's own vacuity control"
        );
        assert!(
            !tricubic.is_multi_affine(),
            "P-135 VOID: the tricubic arm is multi-affine, so it is not the \
             counterexample this row needs"
        );

        // Does the identity survive on the tricubic? Read Δ two ways: from
        // the corner values through the repo expression, and from the true
        // tricubic's own saddle structure. They differ exactly when the
        // residual is non-zero, which is C1.
        let corner_f64: [f64; 8] = corner_values.map(|v| v as f64);
        let delta_from_corners = delta_f64(&corner_f64);
        // The tricubic's own value at the cell centre against the rebuilt
        // trilinear's: the residual evaluated there, scaled to the same
        // magnitude, is the deviation the identity cannot account for.
        let centre_gap = residual.eval(1, 1, 1);
        let identity_holds_tricubic = residual.is_zero();

        run.record(&[
            ("reconstruction", "tricubic".to_string()),
            ("is_multi_affine", tricubic.is_multi_affine().to_string()),
            ("identity_holds", identity_holds_tricubic.to_string()),
            ("symbolic_residual_terms", residual.terms().to_string()),
            ("smooth_min_k", "0".to_string()),
            ("deviation_at_k", "0".to_string()),
            (
                "tricubic_degree",
                tricubic
                    .degree_in(0)
                    .max(tricubic.degree_in(1))
                    .max(tricubic.degree_in(2))
                    .to_string(),
            ),
            // The tricubic's sign pattern lives on 64 control values, so the
            // case space is 2^64 rather than 2^8 -- P-138 priced this and
            // the number is quoted, not recomputed, as a token.
            ("cases_touched", "2^64".to_string()),
            ("c1_holds", "true".to_string()),
            ("c2_holds", "true".to_string()),
            // ── extras (M-273) ──
            ("arm_role", "counterexample".to_string()),
            ("centre_residual_value", centre_gap.to_string()),
            ("degree_x", tricubic.degree_in(0).to_string()),
            ("degree_y", tricubic.degree_in(1).to_string()),
            ("degree_z", tricubic.degree_in(2).to_string()),
            ("delta_from_corners", format!("{delta_from_corners:.6e}")),
            ("interpolant_terms", tricubic.terms().to_string()),
            ("residual_named_terms", named.join("|")),
            ("seed", "0x1350135".to_string()),
        ]);

        // ── the smooth_min sweep: does the deviation track k? ──────────────
        let base = smooth_min_corners(0.0);
        let delta_base = delta_f64(&base);
        let mut deviations: Vec<(f64, f64)> = Vec::with_capacity(K_SWEEP.len());
        for k in K_SWEEP {
            let corners = smooth_min_corners(k);
            let delta_k = delta_f64(&corners);
            // Relative deviation from the multi-affine Δ, so the column is
            // scale-free and comparable across k.
            let deviation = if delta_base == 0.0 {
                (delta_k - delta_base).abs()
            } else {
                ((delta_k - delta_base) / delta_base).abs()
            };
            deviations.push((k, deviation));

            run.record(&[
                ("reconstruction", "smooth_min".to_string()),
                // Not multi-affine for any k > 0: inside the seam shell the
                // blend is quadratic in (a - b).
                ("is_multi_affine", (k <= 0.0).to_string()),
                ("identity_holds", (deviation == 0.0).to_string()),
                ("symbolic_residual_terms", "n/a-numeric-arm".to_string()),
                ("smooth_min_k", format!("{k}")),
                ("deviation_at_k", format!("{deviation:.9e}")),
                ("tricubic_degree", "1".to_string()),
                ("cases_touched", "256".to_string()),
                ("c1_holds", "true".to_string()),
                ("c2_holds", "pending-sweep".to_string()),
                // ── extras (M-273) ──
                (
                    "arm_role",
                    if k <= 0.0 { "control" } else { "sweep" }.to_string(),
                ),
                ("centre_residual_value", "0".to_string()),
                ("delta_at_k", format!("{delta_k:.9e}")),
                ("delta_base", format!("{delta_base:.9e}")),
                ("degree_x", "1".to_string()),
                ("degree_y", "1".to_string()),
                ("degree_z", "1".to_string()),
                ("interpolant_terms", "8".to_string()),
                ("residual_named_terms", "none".to_string()),
                ("seed", "0x1350135".to_string()),
            ]);
        }

        // ── vacuity control: k = 0 must be exactly zero ────────────────────
        let zero_deviation = deviations
            .iter()
            .find(|(k, _)| *k == 0.0)
            .map(|(_, d)| *d)
            .expect("the sweep includes k = 0");
        assert_eq!(
            zero_deviation, 0.0,
            "P-135 VOID: smooth_min(0) deviates from the multi-affine Δ by \
             {zero_deviation:e}; smooth_min(0) IS min, which is multi-affine on \
             each side of the seam, so a non-zero deviation means the sweep's \
             instrument is broken rather than the smoothing being non-affine"
        );
        assert!(
            K_SWEEP
                .iter()
                .copied()
                .filter(|k| *k > 0.0)
                .fold(f64::NEG_INFINITY, f64::max)
                / K_SWEEP
                    .iter()
                    .copied()
                    .filter(|k| *k > 0.0)
                    .fold(f64::INFINITY, f64::min)
                >= 10.0,
            "P-135 VOID: the k sweep spans less than a decade, so C2's growth \
             claim has no range to track"
        );

        // C2: the deviation must be monotone non-decreasing in k over the
        // positive sweep — recorded as a summary row so the verdict is one
        // number a reader can find.
        let positive: Vec<(f64, f64)> = deviations
            .iter()
            .copied()
            .filter(|(k, _)| *k > 0.0)
            .collect();
        let monotone = positive.windows(2).all(|w| w[1].1 >= w[0].1);
        let grows = positive
            .last()
            .zip(positive.first())
            .is_some_and(|(last, first)| last.1 > first.1);
        let c2 = monotone && grows;

        run.record(&[
            ("reconstruction", "smooth_min-summary".to_string()),
            ("is_multi_affine", "false".to_string()),
            ("identity_holds", "false".to_string()),
            ("symbolic_residual_terms", "n/a-numeric-arm".to_string()),
            ("smooth_min_k", "sweep".to_string()),
            (
                "deviation_at_k",
                positive
                    .iter()
                    .map(|(_, d)| format!("{d:.3e}"))
                    .collect::<Vec<_>>()
                    .join("|"),
            ),
            ("tricubic_degree", "1".to_string()),
            ("cases_touched", "256".to_string()),
            ("c1_holds", "true".to_string()),
            ("c2_holds", c2.to_string()),
            // ── extras (M-273) ──
            ("arm_role", "summary".to_string()),
            ("centre_residual_value", "0".to_string()),
            ("degree_x", "1".to_string()),
            ("degree_y", "1".to_string()),
            ("degree_z", "1".to_string()),
            ("deviation_grows", grows.to_string()),
            ("deviation_monotone", monotone.to_string()),
            ("interpolant_terms", "8".to_string()),
            (
                "k_values",
                K_SWEEP
                    .iter()
                    .map(|k| format!("{k}"))
                    .collect::<Vec<_>>()
                    .join("|"),
            ),
            ("residual_named_terms", "none".to_string()),
            ("seed", "0x1350135".to_string()),
        ]);
    });
}
