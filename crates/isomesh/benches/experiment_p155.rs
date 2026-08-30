//! **P-155 — predicting `M-12`'s fitted constant instead of fitting it.**
//!
//! Ticket: R-155. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p155
//! ```
//!
//! Writes `docs/experiments/p-155.csv`.
//!
//! # What was missing
//!
//! `M-12` (`FINDINGS.md:1128`) is *"Marching Cubes' error falls like `h²`,
//! measured"* — mean error `2.7168e-3` at 32³ against `6.5015e-4` at 64³, a
//! ratio of **4.179** against the ideal `((4/31)/(4/63))² = 4.13`. It is a **fit**
//! on **one field**: its test is `the_error_falls_like_h_squared`
//! (`crates/isomesh/src/validate/accuracy/tests.rs:588-605`), which extracts a
//! `Sphere::<f64>::canonical()` at two resolutions and asserts only that the
//! ratio is at least 3.0. The exponent is asserted; the **constant** is not
//! recorded anywhere, and no second field has ever been asked for it.
//!
//! `M-113` then found that a fitted accuracy constant *does not survive*: at
//! `ulp/h = 4` the same quantity read **4.0000** in one alignment and **5.8564**
//! in another. A constant that moves when the grid phase moves is not a property
//! of the algorithm.
//!
//! Strang–Fix says the **order** is a property of the reconstruction filter and
//! nothing else, and it gives the leading error term in closed form. So the
//! trilinear's order-2 should be a derivation rather than a measurement, and the
//! constant should be predictable from the field's own second derivative rather
//! than fitted from the answer. That is this row.
//!
//! Two neighbours depend on the same algebra and neither has it yet: `✗42`
//! measured a **knot shift** at fixed order and found the gain maps to root
//! position only as a lottery over the crossing position, and `P-157` proposes
//! raising the order itself. Neither can be priced without the order-2 baseline
//! this row derives, which is why `P-157`'s own vacuity control is *"the
//! trilinear arm must reproduce exponent 2 in the same harness"*.
//!
//! # The derivation, written out before any number
//!
//! **The order.** The Strang–Fix condition: a filter `phi` reproduces every
//! polynomial of total degree `<= m - 1` iff its Fourier transform has a zero of
//! order `m` at every non-zero point of the dual lattice. The trilinear
//! interpolant is the tensor-product linear B-spline — the hat,
//! `beta1(x) = max(0, 1 - |x|)` — placed on the integer lattice, whose transform
//! is
//!
//! ```text
//! beta1_hat(w) = sinc^2(w/2) = 2 (1 - cos w) / w^2
//! ```
//!
//! (the two forms are identical: `4 sin^2(w/2) = 2(1 - cos w)`). `sin(w/2)`
//! vanishes at every `w = 2 pi k`, and it is **squared**, so every non-zero dual
//! point carries a zero of order exactly **2**. The tensor product
//! `beta_hat(w) = prod_i beta1_hat(w_i)` has, at a dual point `2 pi k` with `r`
//! non-zero components, a zero of order `2r`; the minimum over `k != 0` is at
//! `r = 1`, so `m = 2`. The hat reproduces constants and linears and fails on
//! quadratics.
//!
//! This harness does not quote that. It **measures both halves**:
//!
//! * the Fourier half numerically — `beta1_hat` is cross-checked against a
//!   262,144-point trapezoid quadrature of its own defining integral, its value
//!   on all 26 dual points `k in {-1,0,1}^3 \ {0}` is recorded, and the order of
//!   each zero is read off a log-log slope at `eps = 1e-2` and `1e-3`;
//! * the reproduction half **exactly and symbolically**, in `common::poly`'s
//!   `i128` arithmetic. Variables `0,1,2` are the local cell coordinates
//!   `(u,v,w)` and variables `3,4,5` are the integer cell origin `(i,j,k)`, so
//!   the residual `T_h p - p` is a polynomial in six variables and
//!   `Poly::is_zero()` on it is a statement about **every cell of the lattice at
//!   once**, not about one sampled cell. Degrees 0 and 1 must give the zero
//!   polynomial; degree 2 must not.
//!
//! `measured_order` is `1 +` the largest degree whose every monomial reproduces;
//! `strang_fix_order` is the minimum dual-lattice zero order. C1 asks the two to
//! agree on 2, which is a two-instrument agreement rather than a restatement.
//!
//! **A caveat that is a real property and not a hedge:** the trilinear reproduces
//! the whole *multi-affine* space exactly, including `xy`, `xz`, `yz` and `xyz`
//! — total degree 3. Strang–Fix order is about the **total-degree** space, where
//! `x^2` is the first failure, so the order is 2 and not 4.
//! `c1_multiaffine_reproduced` records the other fact so no reader has to
//! rediscover it.
//!
//! **The constant.** On a grid edge of length `h` in direction `e_i`, the hat is
//! linear interpolation, and the classical remainder is exact at the midpoint:
//!
//! ```text
//! L(1/2) - f(mid) = (h^2 / 8) * d2f/dx_i^2 (xi)      for some xi on the edge
//! ```
//!
//! The extracted vertex sits where `L` vanishes, so its **distance to the true
//! surface** is that residual divided by `|grad f|` — and for a field whose
//! `bound()` is `Exact` the denominator is 1 by definition. Taking the supremum
//! over the surface, and noting that as `h -> 0` the number of crossing edges
//! grows like `h^-2` so a crossing does eventually land at the midpoint of an
//! edge at the worst point of the surface, the `L^inf` asymptotic constant is
//!
//! ```text
//! C = (1/8) * sup_surface max_i |d2f / dx_i^2|
//! ```
//!
//! — the `1/8` the task names, with the field's own second derivative in the
//! sampling directions as the only field-dependent factor.
//!
//! `predicted_constant` is that supremum, measured **without any differencing
//! step**: for every sign-changing axis edge the harness evaluates the three
//! collinear samples `f(p)`, `f(p + h e)` and `f(p + h e / 2)` and takes
//! `|midpoint(f0, f1) - f_mid| / (|grad f| * h^2)`, whose limit is exactly
//! `(1/8)|d2f/dx_i^2|`. It is a prediction in the only sense that matters here:
//! it is computed from the **field and the grid alone**, before any mesh exists
//! and without ever looking at an extractor's output.
//!
//! `analytic_constant` is the closed form where one exists, as an independent
//! check on that supremum. For a true distance field the Hessian at a surface
//! point is `diag(k1, k2, 0)` in the principal frame, so the sup is
//! `max(|k1|, |k2|)` whenever a principal direction is axis-aligned — which it is
//! for every canonical shape here: `sphere` `k = 1` gives `1/8 = 0.125`, and
//! `torus` `k = 1/0.3` at the outer equator, whose meridian direction is the `z`
//! axis, gives `1/(8*0.3) = 0.4166667`. For `box_exact` and `thin_plate` the
//! closed form on the *smooth part* is **0** — flat faces, and the trilinear
//! reproduces a linear field exactly — with **no** second derivative at the
//! creases; that zero is recorded rather than hidden, because the gap between it
//! and the measured supremum is the whole story on those two fields. The four
//! fields whose SDF is not a distance function have no closed-form curvature and
//! record `nan`.
//!
//! **`box_exact` is the row where that goes furthest, and the arithmetic is
//! worth stating before it is measured.** Along any axis edge that *crosses* the
//! box's surface, the box SDF is exactly linear — near the `+x` face, and with
//! the other two `q` components negative, `f` is identically `x - 1` — so the
//! hat's residual on that edge is not small, it is **zero**. Every crossing edge
//! is such an edge, because the regions where two `q` components are positive
//! are strictly outside the solid and produce no sign change. So the measured
//! supremum is zero to `f64` round-off at every resolution, the closed form
//! agrees at zero, and `ratio` and `predicted_stability_ratio` are consequently
//! `inf` by division rather than by any defect. That infinity is the sharpest
//! sentence this experiment produces: the filter theory predicts **no error at
//! all** on a polyhedron and the extractor still has some, so whatever is left
//! is not reconstruction error. `predicted_series` carries the round-off floor
//! so the zero can be read as a zero rather than as a missing measurement, and
//! `thin_plate` — the same code path, the same family of field — returns a
//! non-zero supremum from its rim, which is what licenses `box_exact`'s zero as
//! a measured zero rather than an instrument that cannot report one (`M-44`).
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | polynomial reproduction, degrees 0–3, all 20 monomials, exact over `(u,v,w,i,j,k)` | the total degree | the degree-2 row **is** the control: it must fail, or the order test cannot fire |
//! | Fourier zero order at the 26 dual points `k in {-1,0,1}^3 \ {0}` | the dual point | `beta_hat` against its own quadrature, `beta_hat(0) = 1` included |
//! | one row per reference field, five resolutions each | the field | the four fields whose `bound()` is not `Exact` are recorded as skips |
//! | `M-12`'s own 32³ and 64³ on `sphere`, mean metric | nothing | **yes** — the calibration that says this harness measures what `M-12` measured |
//!
//! Every row carries the same `strang_fix_order`, `measured_order`, `c1_holds`
//! and `c2_holds`: **both clauses are global**, C1 because the filter is one
//! object and C2 because it counts fields. The per-field half of C2 is `c2_hit`.
//!
//! # The resolution ladder, and the arithmetic that chose it
//!
//! `19 | 27 | 35 | 47 | 63` samples per axis over each field's own domain — five
//! resolutions spanning `3.44x` in `h` and `11.9x` in `h^2`, which is what a
//! two-parameter fit needs and more than the registration's *"at least four"*.
//! Every one has `n - 1 = 2m` with **`m` odd**, and that is load-bearing in both
//! directions:
//!
//! * **`n - 1` even** puts `y = 0` on a sample plane, which is the only thing
//!   that makes `thin_plate` extractable at all at these resolutions. Its
//!   half-thickness is `0.4 * 0.5 * 4/64 = 0.0125`, far under every `h` here, so
//!   with `y = 0` *between* two sample planes every corner is outside the plate
//!   and Marching Cubes correctly emits **nothing** (`fields/mod.rs:579-591`).
//!   `M-266` is the note that `y = 0` being a sample plane puts a whole plane of
//!   corners inside it at every level.
//! * **`m` odd** keeps every grid sample off the surface. A sample sits at
//!   `2(k - m)/m`, so `sphere`'s `x^2+y^2+z^2 = 1` needs `4 sum (k_i - m)^2 = m^2`
//!   — even equals odd, impossible; `box_exact`'s face `x = 1` needs `2(k-m) = m`,
//!   impossible; `torus`' outer equator `x = 13/10` needs `20(k-m) = 13m`,
//!   impossible; the plate's rim and face likewise. An exactly-zero sample is
//!   `M-48`'s degenerate crossing and would put an `O(h)` term into a fit about
//!   an `O(h^2)` one.
//!
//! `zero_valued_samples` counts them over the whole ladder, and the argument
//! above covers exactly the **four measurable fields** — it is an argument about
//! four algebraic surfaces, not a property of the ladder. It says nothing about
//! `gyroid`, whose nodal set `sin x cos y + ... = 0` contains lattice points by
//! construction, or about the two noise fields; those are skipped rows and the
//! column is what reports the difference rather than a claim that is quietly
//! false on half the table.
//!
//! `M-12`'s own 32 and 64 have `n - 1` odd and so cannot resolve the plate; they
//! are therefore run as a **separate control** on `M-12`'s own field, rather than
//! folded into a ladder they would break.
//!
//! # The fit
//!
//! `fitted_constant` is the constant of **the same law**, `error = C h^2`, fitted
//! with the exponent held at the derived order: `log C_i = log e_i - 2 log h_i`,
//! and `C = exp(mean log C_i)`. `fitted_constant_ci` is the Student `t` 95%
//! interval on that mean — `t = 2.776` at four degrees of freedom — so a field
//! whose errors do not actually follow `h^2` pays for it with a wide interval
//! rather than with a silently wrong point estimate. `ci_relative_halfwidth` sits
//! beside it precisely so that a `within_ci = true` bought by a useless interval
//! can be told from one that means something. `fitted_exponent` is the free
//! two-parameter slope, recorded so the shape of the law is visible rather than
//! assumed.
//!
//! The error metric is **symmetric Hausdorff** from `validate::accuracy`, which
//! measures the gradient-flow chord and not `|f|` (`validate/accuracy.rs:26-63`).
//! It is comparable against a distance-function prediction only where `bound()`
//! is `Exact`, which is **four of the eight** fields; the other four are recorded
//! as skips with `skip_reason = bound_not_exact`. Their `predicted_constant`,
//! gradient range and crossing-edge count are still measured and recorded —
//! those are properties of the field, not of the extractor — so a skipped row
//! still says something falsifiable.
//!
//! `M-12`'s own metric is the **mean** rather than the max, so it is fitted in
//! parallel and carried as `fitted_constant_mean` / `ratio_mean` /
//! `within_ci_mean`. On `sphere` the mean fit is directly comparable with
//! `M-12`'s `2.7168e-3` and `6.5015e-4`, which the control reproduces in this
//! harness as `m12_mean_32` and `m12_mean_64`.
//!
//! # What C2 can arithmetically reach, computed before the run
//!
//! C2 needs the prediction to land inside the fitted interval on **four of
//! eight** fields. The instrument is valid on four, so C2 needs **four of four**
//! — and two of those four, `box_exact` and `thin_plate`, are polyhedra whose
//! surfaces have **no second derivative** at the creases that dominate their
//! error and whose closed-form constant on the smooth part is zero. The
//! remaining two are `sphere` and `torus`. So C2 is predicted **FALSIFIED**, and
//! the informative part is the `ratio` column rather than the boolean: `M-10`
//! already puts the unit sphere's symmetric Hausdorff at `1.380e-3` at 64³, i.e.
//! `1.380e-3 / (4/63)^2 = 0.342`, against a predicted `0.125` — a ratio near
//! `2.7`, because a symmetric Hausdorff over triangle centroids and lattice seeds
//! also carries the **chordal** deviation of flat triangles across a curved
//! surface, which is `O(h^2)` as well and which no reconstruction-filter argument
//! bounds. That is exactly the registered negative — *"`M-12`'s law is empirical
//! rather than asymptotic and should not be extrapolated"* — and the sharper
//! version this row can state is **which** `O(h^2)` term the filter theory does
//! and does not own.
//!
//! # SHARE, recomputed before the numbers
//!
//! **None, and the registration says none:** *"this converts a fit into a
//! prediction"*. Nothing here changes an extraction path, nothing is timed, and
//! no clause is a ratio of wall clocks — so `M-280`'s `1.45x` governor swing has
//! nothing to bite on. Every number is either an exact integer over an enumerated
//! population (C1) or a distance in world units (C2).
//!
//! # Vacuity controls
//!
//! * **The registration's own:** the fitted constants must differ across the
//!   measurable fields by at least `2x`, or the prediction is matching a
//!   universal constant by accident. Recorded as `fitted_spread`.
//! * **C1 could fail:** `c1_residual_deg2` must be non-zero, or the hat
//!   reproduces everything put to it and `measured_order` is a number that could
//!   not have come out any other way (`M-44`). Predicted `0.25`, from `u - u^2`
//!   at `u = 1/2`.
//! * **C1 was asked anything:** `c1_monomials_tested` must be 20.
//! * **The derived order is derived:** the closed-form transform must agree with
//!   the quadrature of its own defining integral, or the object being
//!   differentiated is not the hat's transform.
//! * **The instrument is calibrated:** `M-12`'s own ratio, reproduced here on
//!   `M-12`'s field, metric and resolutions, must clear `M-12`'s own bar of 3.0.
//!   If this harness cannot reproduce the finding whose constant it claims to
//!   predict, nothing downstream of it means anything.
//! * **Every fit has a population:** every measurable field must report coverage
//!   in both directions at all five resolutions and a strictly positive, finite
//!   error at each, or `log e_i` is undefined and the fit is fitted to nothing.
//! * **Every prediction has a population:** `crossing_edges_finest` must be
//!   non-zero on every field, or its supremum was taken over an empty set.
//! * **At least one field is measurable at all**, or C2 has no denominator.

#![allow(clippy::cast_precision_loss, clippy::float_cmp, clippy::too_many_lines)]

mod common;

use isomesh::fields::{FieldBound, ReferenceField, Sphere};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_cubes::table::is_inside;
use isomesh::validate::{AccuracyConfig, accuracy};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf, Shape3};

use common::poly::{Poly, VARS};

// ─── the registered constants ───────────────────────────────────────────────

/// The order Strang–Fix gives for the tensor-product hat, and the exponent the
/// error law is fitted with. Derived in the header, measured two ways below.
const STRANG_FIX_ORDER: u32 = 2;

/// Samples per axis. `n - 1 = 2m` with `m` odd on every one — see the header for
/// why both halves of that are load-bearing.
const LADDER: [u32; 5] = [19, 27, 35, 47, 63];

/// `M-12`'s coarse resolution, on `M-12`'s own field.
const M12_COARSE: u32 = 32;

/// `M-12`'s fine resolution.
const M12_FINE: u32 = 64;

/// The mean error `M-12` published at 32³, for the calibration column to be read
/// against.
const M12_PUBLISHED_COARSE: f64 = 2.7168e-3;

/// The mean error `M-12` published at 64³.
const M12_PUBLISHED_FINE: f64 = 6.5015e-4;

/// `M-12`'s own acceptance bar on the 32³/64³ ratio
/// (`validate/accuracy/tests.rs:601-604`).
const M12_RATIO_BAR: f64 = 3.0;

/// A gradient shorter than this is a field critical point rather than a surface
/// normal, and dividing a residual by it would manufacture a supremum. The same
/// value `common::metric`'s `GRAD_FLOOR` uses for the same question.
const GRADIENT_FLOOR: f64 = 1e-12;

/// The registration's vacuity bar on the spread of fitted constants.
const SPREAD_BAR: f64 = 2.0;

/// Two-sided 95% Student `t`, indexed by degrees of freedom. Index 0 is unused
/// and is `NAN`, so a fit with one point cannot quietly produce an interval.
const T95: [f64; 9] = [
    f64::NAN,
    12.706,
    4.303,
    3.182,
    2.776,
    2.571,
    2.447,
    2.365,
    2.306,
];

// ─── C1, half one: exact polynomial reproduction ────────────────────────────

/// `common::poly` variable slot for the local cell coordinate `u`.
///
/// `u, v, w` are the local cell coordinates and `i, j, k` the integer cell
/// origin, so a residual that is the zero polynomial in all six is a statement
/// about every cell of the lattice at once.
const U: usize = 0;
/// Local `v`.
const V: usize = 1;
/// Local `w`.
const W: usize = 2;
/// Cell origin `i`.
const I: usize = 3;
/// Cell origin `j`.
const J: usize = 4;
/// Cell origin `k`.
const K: usize = 5;

/// The residual supremum is read on this sub-lattice of the unit cell. `1/8` is
/// fine enough to land exactly on `u = 1/2`, where the degree-2 residual attains
/// its true maximum of `1/4`.
const RESIDUAL_LATTICE: i64 = 8;

/// The world coordinate `local + origin`, as an exact polynomial.
fn world_axis(local: usize, origin: usize) -> Poly {
    Poly::var(local).add(&Poly::var(origin))
}

/// The world coordinate of one corner of the cell: `bit + origin`.
fn corner_axis(bit: u8, origin: usize) -> Poly {
    Poly::constant(i128::from(bit)).add(&Poly::var(origin))
}

/// One factor of the trilinear weight: `t` for the far corner, `1 - t` for the
/// near one. These eight products are the hat `beta(x - k)` restricted to the
/// cell, which is why trilinear interpolation *is* hat quasi-interpolation.
fn weight_axis(bit: u8, local: usize) -> Poly {
    if bit == 1 {
        Poly::var(local)
    } else {
        Poly::constant(1).sub(&Poly::var(local))
    }
}

/// The exact residual `T_h p - p` of the tensor-product hat for the monomial
/// `x^a y^b z^c`, as a polynomial in `(u, v, w, i, j, k)`.
fn hat_residual(exponents: [u32; 3]) -> Poly {
    let [a, b, c] = exponents;
    let local = world_axis(U, I)
        .pow(a)
        .mul(&world_axis(V, J).pow(b))
        .mul(&world_axis(W, K).pow(c));

    let mut interpolant = Poly::zero();
    for corner in 0..8u8 {
        let bits = [corner & 1, (corner >> 1) & 1, (corner >> 2) & 1];
        let value = corner_axis(bits[0], I)
            .pow(a)
            .mul(&corner_axis(bits[1], J).pow(b))
            .mul(&corner_axis(bits[2], K).pow(c));
        let weight = weight_axis(bits[0], U)
            .mul(&weight_axis(bits[1], V))
            .mul(&weight_axis(bits[2], W));
        interpolant = interpolant.add(&value.mul(&weight));
    }
    interpolant.sub(&local)
}

/// Every monomial of exactly this total degree, in `(x, y, z)`.
fn monomials_of_degree(degree: u32) -> Vec<[u32; 3]> {
    let mut out = Vec::new();
    for a in 0..=degree {
        for b in 0..=(degree - a) {
            out.push([a, b, degree - a - b]);
        }
    }
    out
}

/// `sup |residual|` over the unit cell `[0,1]^3` at cell origin `(0,0,0)`, read
/// exactly on the `1/8` sub-lattice and only then converted to `f64`.
///
/// The cell origin is pinned because the degree-3 residual grows with it — `x^3`
/// leaves `3i(u - u^2) + (u - u^3)` — so "the maximum over the lattice" would be
/// unbounded and would say nothing. The unit cell is the normalisation.
fn residual_supremum(residual: &Poly) -> f64 {
    if residual.is_zero() {
        return 0.0;
    }
    let mut den = [1i64; VARS];
    den[U] = RESIDUAL_LATTICE;
    den[V] = RESIDUAL_LATTICE;
    den[W] = RESIDUAL_LATTICE;
    let mut worst = 0.0f64;
    for nu in 0..=RESIDUAL_LATTICE {
        for nv in 0..=RESIDUAL_LATTICE {
            for nw in 0..=RESIDUAL_LATTICE {
                let mut num = [0i64; VARS];
                num[U] = nu;
                num[V] = nv;
                num[W] = nw;
                let (n, d) = residual.eval_ratio(&num, &den);
                worst = worst.max((n as f64 / d as f64).abs());
            }
        }
    }
    worst
}

// ─── C1, half two: the Fourier side ─────────────────────────────────────────

/// `sin t / t`, with its removable singularity filled in.
fn sinc(t: f64) -> f64 {
    if t == 0.0 { 1.0 } else { t.sin() / t }
}

/// `beta1_hat(w) = sinc^2(w/2)`, the transform of the hat `max(0, 1 - |x|)`.
fn hat_transform_1d(w: f64) -> f64 {
    let s = sinc(0.5 * w);
    s * s
}

/// The tensor-product transform, whose zeros on the dual lattice are what
/// Strang–Fix counts.
fn hat_transform(w: [f64; 3]) -> f64 {
    hat_transform_1d(w[0]) * hat_transform_1d(w[1]) * hat_transform_1d(w[2])
}

/// The same transform by quadrature of its own defining integral,
/// `int_{-1}^{1} (1 - |x|) cos(w x) dx` — the sine half vanishes because the hat
/// is even. The step is a negative power of two, so every node, including the
/// kink at `x = 0`, is exact in binary.
fn hat_transform_quadrature(w: f64) -> f64 {
    const STEPS: u32 = 262_144;
    let step = 2.0 / f64::from(STEPS);
    let mut sum = 0.0f64;
    for k in 0..=STEPS {
        let x = -1.0 + f64::from(k) * step;
        let value = (1.0 - x.abs()) * (w * x).cos();
        sum += if k == 0 || k == STEPS {
            0.5 * value
        } else {
            value
        };
    }
    sum * step
}

/// What the Fourier arm found.
#[derive(Debug)]
struct Fourier {
    /// The smallest zero order over the 26 non-zero dual points — the
    /// Strang–Fix order.
    min_zero_order: u32,
    /// The largest `|beta_hat|` actually seen on those points: a zero, up to the
    /// `f64` representation of `2 pi k`.
    max_dual_value: f64,
    /// Worst disagreement between the closed form and the quadrature.
    max_quadrature_deviation: f64,
}

/// Measure the transform, its zeros and their orders.
fn fourier() -> Fourier {
    // The closed form against its own integral, at frequencies on and off the
    // dual lattice. `beta_hat(0) = 1` is the partition-of-unity end of it.
    let mut max_quadrature_deviation = 0.0f64;
    for step in 0..13u32 {
        let w = f64::from(step) * 0.5;
        let deviation = (hat_transform_1d(w) - hat_transform_quadrature(w)).abs();
        max_quadrature_deviation = max_quadrature_deviation.max(deviation);
    }

    let two_pi = core::f64::consts::PI * 2.0;
    let direction = 1.0 / 3.0f64.sqrt();
    let mut min_zero_order = u32::MAX;
    let mut max_dual_value = 0.0f64;
    for kx in -1i32..=1 {
        for ky in -1i32..=1 {
            for kz in -1i32..=1 {
                if kx == 0 && ky == 0 && kz == 0 {
                    continue;
                }
                let base = [
                    two_pi * f64::from(kx),
                    two_pi * f64::from(ky),
                    two_pi * f64::from(kz),
                ];
                max_dual_value = max_dual_value.max(hat_transform(base).abs());

                // The order of the zero, off a log-log slope along
                // (1,1,1)/sqrt(3), which is a direction no dual point is special
                // in.
                let at = |eps: f64| -> f64 {
                    let off = eps * direction;
                    hat_transform([base[0] + off, base[1] + off, base[2] + off]).abs()
                };
                let order = (at(1e-2) / at(1e-3)).ln() / 10.0f64.ln();
                let rounded = order.round();
                assert!(
                    (order - rounded).abs() < 1e-2,
                    "VOID: the zero order at the dual point ({kx}, {ky}, {kz}) measured \
                     {order}, which is not an integer, so the log-log slope is not reading \
                     the order of a zero and the Fourier half of C1 is unmeasured"
                );
                min_zero_order = min_zero_order.min(rounded as u32);
            }
        }
    }
    Fourier {
        min_zero_order,
        max_dual_value,
        max_quadrature_deviation,
    }
}

/// What the reproduction arm found.
#[derive(Debug)]
struct Reproduction {
    /// `sup |T_h p - p|` over the unit cell, worst monomial of each degree
    /// `0..=3`.
    residual_by_degree: [f64; 4],
    /// `1 +` the largest degree reproduced exactly — the measured approximation
    /// order.
    measured_order: u32,
    /// Monomials put through the test.
    monomials: usize,
    /// Whether `xy`, `xz`, `yz` and `xyz` are all reproduced exactly, which is
    /// true and is not the same question as the order.
    multiaffine: bool,
}

/// Exact polynomial reproduction, degrees 0 through 3.
fn reproduction() -> Reproduction {
    let mut residual_by_degree = [0.0f64; 4];
    let mut exact_by_degree = [true; 4];
    let mut monomials = 0usize;
    for (degree, slot) in residual_by_degree.iter_mut().enumerate() {
        for exponents in monomials_of_degree(degree as u32) {
            let residual = hat_residual(exponents);
            if !residual.is_zero() {
                exact_by_degree[degree] = false;
            }
            *slot = slot.max(residual_supremum(&residual));
            monomials += 1;
        }
    }

    let mut reproduced_through = 0u32;
    for (degree, exact) in exact_by_degree.iter().enumerate() {
        if !*exact {
            break;
        }
        reproduced_through = degree as u32;
    }

    let multiaffine = [[1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
        .iter()
        .all(|exponents| hat_residual(*exponents).is_zero());

    Reproduction {
        residual_by_degree,
        measured_order: reproduced_through + 1,
        monomials,
        multiaffine,
    }
}

// ─── C2: the fields ─────────────────────────────────────────────────────────

/// What one grid says about the hat's residual on one field, before any mesh
/// exists.
#[derive(Debug)]
struct Scan {
    /// `sup |residual| / (|grad f| h^2)` over the sign-changing axis edges.
    constant: f64,
    /// Edges the supremum was taken over.
    edges: u64,
    /// Edges whose midpoint gradient was not usable, and which therefore did not
    /// enter the supremum.
    dropped: u64,
    /// Grid samples that are exactly zero — `M-48`'s degenerate crossing.
    zero_samples: u64,
    /// Shortest usable gradient on the sampled surface.
    grad_min: f64,
    /// Longest one. For an `Exact` field both are 1 to within the central
    /// difference.
    grad_max: f64,
}

/// Scan one field on one grid.
fn scan<F>(field: &F, samples: u32) -> Scan
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, h) = common::grid::<f64, _>(field, samples);
    let size = shape.size();
    let at = |x: u32, y: u32, z: u32| -> [f64; 3] {
        [
            origin[0] + f64::from(x) * h,
            origin[1] + f64::from(y) * h,
            origin[2] + f64::from(z) * h,
        ]
    };

    let mut values = vec![0.0f64; shape.element_count()];
    let mut zero_samples = 0u64;
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let value = field.sample(at(x, y, z));
                if value == 0.0 {
                    zero_samples += 1;
                }
                values[shape.linearize([x, y, z]) as usize] = value;
            }
        }
    }

    let mut constant = 0.0f64;
    let mut edges = 0u64;
    let mut dropped = 0u64;
    let mut grad_min = f64::INFINITY;
    let mut grad_max = 0.0f64;
    for axis in 0..3usize {
        let mut step = [0u32; 3];
        step[axis] = 1;
        for z in 0..size[2] - step[2] {
            for y in 0..size[1] - step[1] {
                for x in 0..size[0] - step[0] {
                    let near = values[shape.linearize([x, y, z]) as usize];
                    let far =
                        values[shape.linearize([x + step[0], y + step[1], z + step[2]]) as usize];
                    if is_inside(near) == is_inside(far) {
                        continue;
                    }
                    let mut mid = at(x, y, z);
                    mid[axis] += 0.5 * h;
                    let residual = (f64::midpoint(near, far) - field.sample(mid)).abs();
                    let g = field.gradient(mid);
                    let norm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                    if !norm.is_finite() || norm < GRADIENT_FLOOR {
                        dropped += 1;
                        continue;
                    }
                    edges += 1;
                    grad_min = grad_min.min(norm);
                    grad_max = grad_max.max(norm);
                    constant = constant.max(residual / (norm * h * h));
                }
            }
        }
    }

    Scan {
        constant,
        edges,
        dropped,
        zero_samples,
        grad_min,
        grad_max,
    }
}

/// Marching Cubes at its defaults — the configuration `M-12` measured, from
/// `validate/accuracy/tests.rs:68-72`.
fn extract<F>(field: &F, shape: &RuntimeShape3, origin: [f64; 3], h: f64) -> MeshBuffer<f64>
where
    F: Sdf<Scalar = f64>,
{
    let mut mc = MarchingCubes::<f64>::new();
    let mut mesh = MeshBuffer::<f64>::new();
    mc.extract(field, shape, origin, h, &mut mesh)
        .expect("every ladder grid has at least two samples per axis");
    mesh
}

/// Both error metrics at one resolution: the symmetric Hausdorff this row is
/// scored on, and the mean `M-12` used.
fn errors_at<F>(field: &F, samples: u32) -> (f64, f64, bool)
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, h) = common::grid::<f64, _>(field, samples);
    let mesh = extract(field, &shape, origin, h);
    let cfg = AccuracyConfig::from_cell_size(h).expect("every ladder cell size is positive");
    let report = accuracy(&mesh.positions, &mesh.indices, field, &shape, origin, &cfg)
        .expect("the mesh and the grid it came from belong to each other");
    (
        report.symmetric_hausdorff(),
        report.mean_absolute_error(),
        report.has_coverage(),
    )
}

/// A one-parameter fit of `error = C h^order`, plus the free two-parameter
/// slope beside it.
#[derive(Debug)]
struct Fit {
    /// `exp(mean(log e_i - order * log h_i))`.
    constant: f64,
    /// Lower end of the Student `t` 95% interval on that mean, exponentiated.
    ci_lo: f64,
    /// The upper end.
    ci_hi: f64,
    /// `(hi - lo) / (2 C)`, so a wide interval cannot pass for a tight one.
    half_width_rel: f64,
    /// The unconstrained slope of `log e` against `log h`.
    exponent: f64,
}

impl Fit {
    /// Whether a prediction lands inside the interval.
    fn contains(&self, predicted: f64) -> bool {
        predicted >= self.ci_lo && predicted <= self.ci_hi
    }
}

/// Fit the law at the derived order, and separately measure the order the data
/// actually shows.
fn fit_law(cell_sizes: &[f64], errors: &[f64]) -> Fit {
    let n = cell_sizes.len();
    assert_eq!(n, errors.len(), "one error per cell size");
    assert!(
        (4..=8).contains(&n),
        "the fit needs at least four resolutions and its t-table stops at eight"
    );
    let count = n as f64;
    let logs: Vec<(f64, f64)> = cell_sizes
        .iter()
        .zip(errors)
        .map(|(h, e)| (h.ln(), e.ln()))
        .collect();

    let mean_x = logs.iter().map(|(x, _)| *x).sum::<f64>() / count;
    let mean_y = logs.iter().map(|(_, y)| *y).sum::<f64>() / count;
    let sxy: f64 = logs.iter().map(|(x, y)| (x - mean_x) * (y - mean_y)).sum();
    let sxx: f64 = logs.iter().map(|(x, _)| (x - mean_x) * (x - mean_x)).sum();

    let order = f64::from(STRANG_FIX_ORDER);
    let held: Vec<f64> = logs.iter().map(|(x, y)| y - order * x).collect();
    let mean_c = held.iter().sum::<f64>() / count;
    let variance = held
        .iter()
        .map(|c| (c - mean_c) * (c - mean_c))
        .sum::<f64>()
        / (count - 1.0);
    let half = T95[n - 1] * (variance / count).sqrt();

    let constant = mean_c.exp();
    let ci_lo = (mean_c - half).exp();
    let ci_hi = (mean_c + half).exp();
    Fit {
        constant,
        ci_lo,
        ci_hi,
        half_width_rel: (ci_hi - ci_lo) / (2.0 * constant),
        exponent: sxy / sxx,
    }
}

/// The ladder, where the field's `bound()` lets it mean anything.
#[derive(Debug)]
struct Ladder {
    /// Symmetric Hausdorff at each resolution.
    hausdorff: Vec<f64>,
    /// Whether every resolution reported coverage in both directions.
    covered: bool,
    /// The fit this row is scored on.
    haus: Fit,
    /// The same fit on `M-12`'s metric.
    mean: Fit,
}

/// Everything one reference field contributes.
#[derive(Debug)]
struct Measured {
    /// The field's own name.
    name: &'static str,
    /// Its declared bound, which decides whether the Hausdorff arm runs.
    bound: FieldBound,
    /// `predicted_constant` at each ladder resolution.
    predicted: Vec<f64>,
    /// Crossing edges the finest supremum was taken over.
    edges: u64,
    /// Edges dropped at the finest resolution for an unusable gradient.
    dropped: u64,
    /// Exactly-zero samples over the whole ladder.
    zero_samples: u64,
    /// Shortest surface gradient at the finest resolution.
    grad_min: f64,
    /// Longest one.
    grad_max: f64,
    /// The closed form where one exists, `NAN` otherwise.
    analytic: f64,
    /// `None` where `bound()` is not `Exact`.
    ladder: Option<Ladder>,
}

impl Measured {
    /// The headline prediction: the supremum at the finest resolution, which is
    /// the closest this ladder gets to the asymptotic one.
    fn predicted_constant(&self) -> f64 {
        *self
            .predicted
            .last()
            .expect("the ladder has at least one resolution")
    }

    /// How much the prediction moved across the ladder. Near `1` for a `C^2`
    /// surface, whose second derivative converges; larger says it does not,
    /// which is what a crease looks like. `inf` where the supremum is an exact
    /// zero at some resolution, which is `box_exact` and is derived in the
    /// header rather than being a division that got away.
    fn prediction_stability(&self) -> f64 {
        let lo = self.predicted.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = self.predicted.iter().copied().fold(0.0f64, f64::max);
        hi / lo
    }
}

/// The closed-form asymptotic constant, where the surface has one.
///
/// `sup max_i |d2f/dx_i^2| / 8`, with the Hessian of a distance function at a
/// surface point being `diag(k1, k2, 0)` in the principal frame and a principal
/// direction axis-aligned on every canonical shape here. `box_exact` and
/// `thin_plate` are polyhedra: zero on the smooth part, undefined at the creases.
fn analytic_constant(name: &str) -> f64 {
    match name {
        "sphere" => 1.0 / 8.0,
        "torus" => 1.0 / (0.3 * 8.0),
        "box_exact" | "thin_plate" => 0.0,
        _ => f64::NAN,
    }
}

/// The name of a bound, for the CSV.
fn bound_name(bound: FieldBound) -> &'static str {
    match bound {
        FieldBound::Exact => "exact",
        FieldBound::Lipschitz { .. } => "lipschitz",
        FieldBound::Underestimate { .. } => "underestimate",
        FieldBound::Unbounded => "unbounded",
    }
}

/// Measure one reference field: the prediction on every ladder grid, and the fit
/// where the error metric is meaningful.
fn measure<F>(name: &'static str, field: &F) -> Measured
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let scans: Vec<Scan> = LADDER.iter().map(|n| scan(field, *n)).collect();
    let finest = scans.last().expect("the ladder is not empty");
    let bound = field.bound();

    let ladder = if bound.is_exact() {
        let (lo, hi) = field.domain();
        let cell_sizes: Vec<f64> = LADDER
            .iter()
            .map(|n| (hi[0] - lo[0]) / f64::from(n - 1))
            .collect();
        let mut hausdorff = Vec::with_capacity(LADDER.len());
        let mut mean = Vec::with_capacity(LADDER.len());
        let mut covered = true;
        for samples in LADDER {
            let (worst, average, has_coverage) = errors_at(field, samples);
            hausdorff.push(worst);
            mean.push(average);
            covered &= has_coverage;
        }
        let haus = fit_law(&cell_sizes, &hausdorff);
        let mean_fit = fit_law(&cell_sizes, &mean);
        Some(Ladder {
            hausdorff,
            covered,
            haus,
            mean: mean_fit,
        })
    } else {
        None
    };

    Measured {
        name,
        bound,
        predicted: scans.iter().map(|s| s.constant).collect(),
        edges: finest.edges,
        dropped: finest.dropped,
        zero_samples: scans.iter().map(|s| s.zero_samples).sum(),
        grad_min: finest.grad_min,
        grad_max: finest.grad_max,
        analytic: analytic_constant(name),
        ladder,
    }
}

/// `M-12`'s own measurement, in this harness: its field, its metric, its two
/// resolutions.
fn m12_control() -> (f64, f64) {
    let field = Sphere::<f64>::canonical();
    (
        errors_at(&field, M12_COARSE).1,
        errors_at(&field, M12_FINE).1,
    )
}

// ─── formatting ─────────────────────────────────────────────────────────────

/// A value that can span orders of magnitude. `nan` is written as itself rather
/// than as a zero, because a skipped measurement is not a measured zero.
fn num(value: f64) -> String {
    if value.is_nan() {
        String::from("nan")
    } else {
        format!("{value:.6e}")
    }
}

/// A ratio or an order, where fixed point reads better.
fn plain(value: f64) -> String {
    if value.is_nan() {
        String::from("nan")
    } else {
        format!("{value:.6}")
    }
}

/// A list of numbers as one CSV-safe token.
fn series(values: &[f64]) -> String {
    values.iter().map(|v| num(*v)).collect::<Vec<_>>().join("|")
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-155");

    common::experiment::run(prereg, |run| {
        // ── C1, both instruments ────────────────────────────────────────────
        let spectrum = fourier();
        let repro = reproduction();
        let strang_fix_order = spectrum.min_zero_order;
        let c1 = strang_fix_order == STRANG_FIX_ORDER && repro.measured_order == STRANG_FIX_ORDER;

        println!(
            "C1  dual-lattice min zero order {}  (max |beta_hat| on the 26 points {:e}, \
             worst quadrature deviation {:e})",
            spectrum.min_zero_order, spectrum.max_dual_value, spectrum.max_quadrature_deviation
        );
        println!(
            "C1  exact reproduction, sup residual on the unit cell by degree: \
             0 {:.6}  1 {:.6}  2 {:.6}  3 {:.6}  over {} monomials -> measured order {}",
            repro.residual_by_degree[0],
            repro.residual_by_degree[1],
            repro.residual_by_degree[2],
            repro.residual_by_degree[3],
            repro.monomials,
            repro.measured_order
        );

        // ── C2, the eight fields ────────────────────────────────────────────
        let mut fields: Vec<Measured> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            fields.push(measure(name, &field));
        });

        // ── the calibration control, on M-12's own terms ────────────────────
        let (m12_coarse, m12_fine) = m12_control();
        let m12_ratio = m12_coarse / m12_fine;
        println!(
            "M-12 control  mean error 32^3 {m12_coarse:.9} (published {M12_PUBLISHED_COARSE:.9}), \
             64^3 {m12_fine:.9} (published {M12_PUBLISHED_FINE:.9}), ratio {m12_ratio:.3}"
        );

        // ── the vacuity controls, before any row is written ─────────────────
        assert!(
            repro.residual_by_degree[2] > 0.0,
            "VOID: the hat reproduced every quadratic put to it, so `measured_order` is a \
             number that could not have come out any other way (M-44) and C1 is unmeasured"
        );
        assert_eq!(
            repro.monomials, 20,
            "VOID: {} monomials were tested rather than the 20 of degrees 0 through 3, so \
             the reproduction arm did not ask what it says it asked",
            repro.monomials
        );
        assert!(
            spectrum.max_quadrature_deviation < 1e-6,
            "VOID: the closed-form transform and the quadrature of its own defining integral \
             disagree by {:e}, so `sinc^2(w/2)` is not the object whose zeros are being \
             counted and the derived order is a quotation rather than a measurement",
            spectrum.max_quadrature_deviation
        );
        assert!(
            m12_ratio >= M12_RATIO_BAR,
            "VOID: M-12's own 32^3/64^3 ratio reads {m12_ratio} here against its own bar of \
             {M12_RATIO_BAR}, so this harness does not reproduce the finding whose constant \
             it claims to predict and every ratio below is against an unmoored baseline"
        );
        for f in &fields {
            assert!(
                f.edges > 0,
                "VOID: {}: no sign-changing axis edge at the finest resolution, so its \
                 predicted constant is a supremum over an empty set",
                f.name
            );
            if let Some(ladder) = &f.ladder {
                assert!(
                    ladder.covered,
                    "VOID: {}: accuracy reported no coverage at some ladder resolution, so \
                     one of the five points of its fit is not a measurement of anything",
                    f.name
                );
                assert!(
                    ladder.hausdorff.iter().all(|e| *e > 0.0 && e.is_finite()),
                    "VOID: {}: a zero or non-finite symmetric Hausdorff on the ladder, so \
                     `log e` is undefined and the fit is fitted to nothing: {:?}",
                    f.name,
                    ladder.hausdorff
                );
            }
        }
        let measurable: Vec<&Measured> = fields.iter().filter(|f| f.ladder.is_some()).collect();
        assert!(
            !measurable.is_empty(),
            "VOID: not one field has an Exact bound, so C2 has no denominator and the \
             symmetric Hausdorff means nothing anywhere in this run"
        );

        let constants: Vec<f64> = measurable
            .iter()
            .filter_map(|f| f.ladder.as_ref().map(|l| l.haus.constant))
            .collect();
        let spread = constants.iter().copied().fold(0.0f64, f64::max)
            / constants.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            spread >= SPREAD_BAR,
            "VOID: the fitted constants span only {spread}x across the {} measurable fields, \
             under the registered bar of {SPREAD_BAR}x, so any agreement with the prediction \
             is agreement with a universal constant: {constants:?}",
            measurable.len()
        );

        // ── the verdicts ────────────────────────────────────────────────────
        let hits = fields
            .iter()
            .filter(|f| {
                f.ladder
                    .as_ref()
                    .is_some_and(|l| l.haus.contains(f.predicted_constant()))
            })
            .count();
        let c2 = hits >= 4;

        // ── the rows ────────────────────────────────────────────────────────
        let resolutions = LADDER
            .iter()
            .map(u32::to_string)
            .collect::<Vec<_>>()
            .join("|");
        for f in &fields {
            let predicted = f.predicted_constant();
            let (fitted, ci_lo, ci_hi, ratio, within, exponent, half_width) = match &f.ladder {
                Some(l) => (
                    l.haus.constant,
                    l.haus.ci_lo,
                    l.haus.ci_hi,
                    l.haus.constant / predicted,
                    l.haus.contains(predicted),
                    l.haus.exponent,
                    l.haus.half_width_rel,
                ),
                None => (
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    false,
                    f64::NAN,
                    f64::NAN,
                ),
            };
            let (mean_constant, mean_ci_lo, mean_ci_hi, mean_ratio, mean_within) = match &f.ladder {
                Some(l) => (
                    l.mean.constant,
                    l.mean.ci_lo,
                    l.mean.ci_hi,
                    l.mean.constant / predicted,
                    l.mean.contains(predicted),
                ),
                None => (f64::NAN, f64::NAN, f64::NAN, f64::NAN, false),
            };
            let hausdorff = f
                .ladder
                .as_ref()
                .map_or_else(|| String::from("nan"), |l| series(&l.hausdorff));

            println!(
                "{:>15} {:>14}  predicted {:.6e} (analytic {})  fitted {:.6e} \
                 [{:.6e}, {:.6e}]  ratio {}  exponent {}  within_ci {}",
                f.name,
                bound_name(f.bound),
                predicted,
                plain(f.analytic),
                fitted,
                ci_lo,
                ci_hi,
                plain(ratio),
                plain(exponent),
                within
            );

            run.record(&[
                ("field", f.name.to_string()),
                ("resolution_series", resolutions.clone()),
                ("fitted_constant", num(fitted)),
                (
                    "fitted_constant_ci",
                    format!("{}|{}", num(ci_lo), num(ci_hi)),
                ),
                ("predicted_constant", num(predicted)),
                ("ratio", plain(ratio)),
                ("within_ci", within.to_string()),
                ("strang_fix_order", strang_fix_order.to_string()),
                ("measured_order", repro.measured_order.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                // ── extras (M-273) ──────────────────────────────────────────
                ("analytic_constant", plain(f.analytic)),
                ("bound", bound_name(f.bound).to_string()),
                ("c1_monomials_tested", repro.monomials.to_string()),
                ("c1_multiaffine_reproduced", repro.multiaffine.to_string()),
                ("c1_residual_deg0", plain(repro.residual_by_degree[0])),
                ("c1_residual_deg1", plain(repro.residual_by_degree[1])),
                ("c1_residual_deg2", plain(repro.residual_by_degree[2])),
                ("c1_residual_deg3", plain(repro.residual_by_degree[3])),
                ("c2_hit", within.to_string()),
                ("c2_hits", hits.to_string()),
                ("c2_measurable_fields", measurable.len().to_string()),
                ("ci_relative_halfwidth", plain(half_width)),
                ("crossing_edges_finest", f.edges.to_string()),
                ("edges_dropped_finest", f.dropped.to_string()),
                ("fitted_constant_mean", num(mean_constant)),
                (
                    "fitted_constant_mean_ci",
                    format!("{}|{}", num(mean_ci_lo), num(mean_ci_hi)),
                ),
                ("fitted_exponent", plain(exponent)),
                ("fitted_spread", plain(spread)),
                (
                    "fourier_dual_lattice_max_value",
                    num(spectrum.max_dual_value),
                ),
                (
                    "fourier_quadrature_max_deviation",
                    num(spectrum.max_quadrature_deviation),
                ),
                (
                    "fourier_zero_order_min",
                    spectrum.min_zero_order.to_string(),
                ),
                (
                    "grad_norm_range",
                    format!("{}|{}", plain(f.grad_min), plain(f.grad_max)),
                ),
                ("hausdorff_series", hausdorff),
                ("hausdorff_valid", f.bound.is_exact().to_string()),
                ("m12_mean_32", num(m12_coarse)),
                ("m12_mean_64", num(m12_fine)),
                ("m12_ratio", plain(m12_ratio)),
                ("predicted_series", series(&f.predicted)),
                ("predicted_stability_ratio", plain(f.prediction_stability())),
                ("ratio_mean", plain(mean_ratio)),
                (
                    "skip_reason",
                    String::from(if f.bound.is_exact() {
                        "none"
                    } else {
                        "bound_not_exact"
                    }),
                ),
                ("within_ci_mean", mean_within.to_string()),
                ("zero_valued_samples", f.zero_samples.to_string()),
            ]);
        }

        println!();
        println!(
            "C1 {}: derived order {strang_fix_order} from the dual-lattice zeros, measured \
             order {} from exact reproduction, against the registered {STRANG_FIX_ORDER}",
            if c1 { "HELD" } else { "FALSIFIED" },
            repro.measured_order
        );
        println!(
            "C2 {}: the prediction lands inside the fitted interval on {hits} of 8 fields \
             against a bar of 4 — and the instrument is valid on {} of them, so {} was the \
             most C2 could reach",
            if c2 { "HELD" } else { "FALSIFIED" },
            measurable.len(),
            measurable.len()
        );
        println!(
            "vacuity: the fitted constants span {spread:.3}x across the measurable fields, \
             against the registered bar of {SPREAD_BAR}x"
        );
    });
}
