//! **P-159 — how far from optimal the extractor is, and the null that would be worth the most.**
//!
//! Ticket: R-159. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p159
//! ```
//!
//! Writes `docs/experiments/p-159.csv`.
//!
//! # What was missing
//!
//! Every accuracy number in this repository is a **numerator**. `M-12`
//! (FINDINGS.md:1128) measured Marching Cubes' error falling like `h²` — mean
//! `2.7168e-3` at 32³ against `6.5015e-4` at 64³, ratio `4.179` against the ideal
//! `4.13` — and stopped there. `M-65` (:1181) measured the *gradient* converging
//! at `h²`. `F-007`/`M-250` (:1304) measured crossing refinement buying 13–15% on
//! curved fields and `1.000×` on the CSG one. `✗42`/`M-359` (:9708) measured a
//! root-position gain that turned out to be a lottery over where the root falls.
//! Each of those says *did this help*. Not one says *how much is left*, because
//! saying that needs a **denominator**, and a denominator for "recover a level set
//! from `n` function values" is not an accuracy measurement at all — it is a
//! theorem about information.
//!
//! Information-based complexity owns that theorem, and this row is the first in
//! the repository to open it. The currency is right: field evaluation is this
//! crate's dominant cost (`M-279`, `P-121`), and the IBC bound is denominated in
//! **function values**, which is exactly what `Sdf::sample` spends.
//!
//! # The literature, read rather than quoted
//!
//! The registration names **Krieg & Ullrich, `arXiv:2602.02066`**, *Approximation
//! of Functions: Optimal Sampling and Complexity* (2026-02-02, 191 pages), and it
//! is genuinely in the corpus at `10.48550/arXiv.2602.02066` — downloaded,
//! converted, embedded, readable. Five places in it were read for this harness,
//! and the page numbers are here so a later reader can re-check them rather than
//! trust this file:
//!
//! - **p. 13** — the definition. The `n`-th **minimal (worst-case) error**, also
//!   called the `n`-th **sampling number**, is
//!   `g_n(F, Y) = inf over {x_1..x_n ⊂ D, Φ: ℝⁿ → Y} sup over {f ∈ F} ‖f − Φ(f(x_1),…,f(x_n))‖_Y`.
//!   A lower bound on `g_n` says *no algorithm with a smaller error can exist,
//!   however the points and the reconstruction are chosen*. **Remark 2.2**: this
//!   is the non-adaptive quantity, and adaptivity does not help here (their
//!   Theorem 9.3).
//! - **p. 92** — the class. `W_q^s(D)` for integer `s` with `s > d/q` on a bounded
//!   domain `D ⊂ ℝ^d`, norm `‖f‖_∞ + |f|_{W_q^s}`; `F_q^s` is its **unit ball**.
//! - **p. 93 — Theorem 6.10**, the load-bearing one, attributed to Krieg &
//!   Sonnleitner (2024, Thm 0.1) and (2025, Thm 1). Let `D` be a bounded
//!   **convex** domain in `ℝ^d` (or a connected compact Riemannian `d`-manifold),
//!   `1 ≤ p, q ≤ ∞`, `s ∈ ℕ` with `s > d/q`. For **any** point set `P_n ⊂ D`,
//!   with constants independent of `n` and `P_n`:
//!   `rad(P_n, F_q^s, L_p) ≍ ‖dist(·,P_n)‖_{L_γ}^s` when `p < q`,
//!   `γ = s(1/p − 1/q)^{-1}`, and
//!   `rad(P_n, F_q^s, L_p) ≍ ‖dist(·,P_n)‖_∞^{s − d(1/q − 1/p)}` when `p ≥ q`.
//! - **p. 94 — Proposition 6.11**:
//!   `inf over P_n of ‖dist(·,P_n)‖_{L_γ} ≍ n^{-1/d}` for every `0 < γ ≤ ∞`.
//!   Composed with Theorem 6.10 this **is** the rate:
//!   `g_n(F_q^s, L_p) ≍ n^{-s/d + (1/q − 1/p)_+}`.
//! - **p. 16 — Proposition 2.4** states the `d = 1` case independently,
//!   `g_n(F_q^s, L_p) ≍ n^{-s + (1/q − 1/p)_+}` at equidistant points, and
//!   **p. 121** proves the `s = 1` multivariate case from scratch with an explicit
//!   constant, `e(Q_n, F_Lip^d) ≥ (1/8) n^{-1/d}`, via the fooling function
//!   `f(x) = min_j ‖x − x_j‖_∞`. Those two are the transcription checks this file
//!   asserts (vacuity control 1), because a general-`d` formula is exactly the
//!   kind of thing that gets a `d` in the wrong place.
//!
//! **For this crate: `d = 3` and `p = q = ∞`, so `1/q − 1/p = 0` and the rate is
//! `g_n ≍ n^{-s/3}`, i.e. `Θ(h^s)` at grid spacing `h ≍ n^{-1/3}`.** At `s = 2`
//! that is `n^{-2/3} = Θ(h²)` — the number C3 asks about.
//!
//! A second source was read as a cross-check, and because it is a different
//! theorem it is named as one: Bonito, Canuto, Nochetto & Veeser, *Adaptive finite
//! element methods* (Acta Numerica 2024, corpus `10.1017/s0962492924000011`),
//! **p. 67**: *"The maximal convergence order of the error under uniform
//! refinement is `‖∇(u − u_T)‖_{L²(Ω)} = O(h^n)`"* with `n` the **polynomial
//! degree**. That is the other half of the story — with a fixed reconstruction
//! space the order is capped by the space, not by the function's smoothness — and
//! it is why this harness records two arms rather than one.
//!
//! # What could NOT be verified, stated before the numbers
//!
//! Three things, each carried in every row's `class_assumptions`:
//!
//! 1. **The implied constants in `≍` are not in the source.** Theorem 6.10 and
//!    Proposition 6.11 are order statements whose constants depend on
//!    `d, s, p, q, D` and are evaluated nowhere in the paper. So `constant_gap` is
//!    **not** measured against Krieg & Ullrich's constant. It is measured against
//!    an explicit floor derived here by their own fooling-function method, and
//!    every row says so.
//! 2. **The bound is on `‖f − f̃‖_{L_p}`, not on the level set.** IBC recovers the
//!    *function*; an extractor emits `{f̃ = 0}`. The reduction between the two is
//!    derived here and is not in the source.
//! 3. **Regularity is measured at the resolutions benched, not in the limit.** No
//!    asymptotic class membership is claimed anywhere in this file.
//!
//! # The floor, with an explicit constant
//!
//! Krieg & Ullrich's method (p. 12, p. 121) is a **fooling function**: a `f_0`
//! that vanishes on every sample point, so `N_n(f_0) = N_n(−f_0)` and no algorithm
//! reading only those samples can tell `f_0` from `−f_0`. Instantiated on this
//! crate's grid at spacing `h`,
//!
//! ```text
//! f_0(x) = A · sin(π x_1/h) · sin(π x_2/h) · sin(π x_3/h)
//! ```
//!
//! vanishes on all of `hℤ³`, and its largest pure `s`-th derivative is
//! `A (π/h)^s`. Pinning that to the field's own measured seminorm `M_s` gives
//! `A = M_s (h/π)^s`, hence a **function-recovery** floor of `A`.
//!
//! For the **level set**: `f` and `f ± f_0` have identical grid data and their zero
//! sets are displaced by about `A / |∇f|`, so nothing reading only the grid can be
//! closer than `A / (2·G_max)` to all three, where `G_max = ‖∇f‖_∞` on the surface
//! band. That is `minimal_error_finest`, and both directions of its sloppiness are
//! **conservative for the floor and against Marching Cubes**: `M_s` is estimated
//! from pure axis derivatives only, while the seminorm also ranges over mixed
//! ones, so `M_s` is an underestimate; and `G_max` is a maximum rather than a
//! typical value. Both shrink the floor and therefore *inflate* `constant_gap`. A
//! small `constant_gap` reported here is a real one.
//!
//! # Arms
//!
//! The class is the whole content of this row, so it is an **arm**, not a
//! constant. Nine fields × two arms = eighteen rows.
//!
//! | arm | what it varies | `is_control` |
//! |---|---|---|
//! | `pinned_w2` | class fixed at `F_∞^2` — the class Marching Cubes' trilinear model is second-order exact for, and the class C3 names | false |
//! | `measured_max` | class fixed at `F_∞^s` with `s` the largest integer smoothness the field's own surface **measurably** supports, `1 ≤ s ≤ 4` | false |
//!
//! and the field roster carries its own control:
//!
//! | field | why | `is_control` |
//! |---|---|---|
//! | the eight `for_each_reference_field!` fields | the roster every other row is measured on | false |
//! | `lens_control` — two unit spheres `1.0` apart, intersected, centred **off-lattice** at `(0.11, 0.07, 0.13)`, bench-local | **curved faces and one circular crease**, so its regularity must measure near `1` where a sphere measures near `4`, *and* its error must still be a fittable `Θ(h)`. Two earlier versions failed, both for reasons now recorded in the code: the `L¹` ball `‖x‖_1 − 1` measured `s = 1.000000` exactly but Marching Cubes reproduces a piecewise-linear field *exactly*, so its error sat at the `f64` floor with no exponent; and an origin-centred lens put its crease at `x = 0`, a grid plane at **every** ladder resolution, so no triangle spanned the kink and the control measured `Θ(h²)` like a sphere | **true** |
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration says **`SHARE: none -- this is a denominator, not a saving.`**
//! Discharged literally: this harness moves nothing, configures nothing and
//! proposes nothing. It calls `MarchingCubes::new()` with every default
//! (`FaceAmbiguity::Separate`, `InteriorAmbiguity::Ignore`, and
//! `set_crossing_refinement` untouched at `0` — deliberately, because `M-250`
//! measured refinement moving the constant by 13–15% and this row must denominate
//! the *shipped* configuration). There is no `1/(1 − share/factor)` ceiling to
//! compute because nothing here claims a speedup.
//!
//! # The instrument
//!
//! **Resolution ladder `17, 25, 33, 49, 65` samples per axis** — five points, not
//! the house three, because two of the registered columns (`measured_error_rate`,
//! and `ratio` through it) are *fitted exponents* and a three-point fit has no
//! residual to report. `65` is the ceiling: the golden fixture tops out at `33`,
//! and `65³ = 274 625` samples per field per pass keeps the whole bench inside the
//! two-minute budget. `n_samples` is the **finest** point's function-value count,
//! because that is the `n` at which `constant_gap` is evaluated; the whole ladder
//! is in the extra column `n_samples_ladder`.
//!
//! **The error functional is `max over the mesh of |f(x)| / ‖∇f(x)‖`, over every
//! vertex, every edge midpoint and every triangle centroid.** To first order that
//! is the distance from `x` to the true zero set, so the maximum is the one-sided
//! Hausdorff distance `sup over {x in mesh} dist(x, {f = 0})` discretised on the
//! mesh — the same quantity `validate::accuracy`'s `mesh_to_field` direction
//! estimates.
//!
//! **Sampling the vertices alone would have been wrong, and `box_exact` is the
//! proof — this harness read zero on it before the midpoints were added.** Its
//! half-extents are `1` and every ladder spacing `h = 4/(N−1)` divides `3`, so
//! every box face lands on a grid plane and the field is exactly linear along
//! every cut grid edge: Marching Cubes finds those roots exactly, and a
//! vertex-only functional reads **identically zero at all five resolutions**. The
//! error is not zero — `M-250` measured `box_exact`'s symmetric Hausdorff at
//! `1.443e-1` both before and after crossing refinement (`1.000×`) — because the
//! mesh *truncates* every box edge with a chord whose midpoint sits about `h/2`
//! inside the solid. Only a functional that looks between the vertices can see
//! that. It costs `O(T)` extra samples, leaves the smooth fields' `O(h²)` rate
//! alone (a chord across curvature `κ` deviates by `O((hκ)²)`), and turns
//! `box_exact` from an unfittable zero into the `Θ(h)` its own measured
//! regularity predicts.
//!
//! `validate::accuracy` is not used, for one reason that is not a small one:
//! `accuracy` reads the SDF as a distance oracle, and only four of the eight
//! reference fields have `bound() == Exact` — `gyroid` is `Lipschitz`,
//! `csg_difference` is `Underestimate { q: 0.5 }`, `fbm_terrain` and
//! `noise_cavity` are `Unbounded`. A functional that is meaningless on half the
//! roster cannot carry a rate fit across it. The gradient-normalised residual
//! needs only a non-vanishing gradient, which is the same transversality
//! hypothesis `P-171` measures, and the count of sampled points at the gradient
//! floor is a column.
//!
//! **Regularity is measured, not assumed — this is the registration's own vacuity
//! control.** For `f` with `s` bounded derivatives the centred `m`-th finite
//! difference obeys `Δ^m_t f = t^m f^{(m)}(ξ)` exactly, so `log ‖Δ^m_t f‖_∞`
//! against `log t` has slope `min(m, s)`. Probing orders `1..=4` therefore reads
//! `s` off directly as the order-4 slope, with the lower orders as a consistency
//! ladder. The probes are **every vertex of the finest mesh**, unstrided: that is
//! the surface band the extractor actually works in, and measuring regularity
//! anywhere else would answer a different question — `csg_difference`'s interior
//! crease is not on its zero set, and `box_exact`'s twelve edges are. Striding to a
//! cap was tried and removed: whether a crease gets sampled then depends on the
//! stride rather than on the field, which is the wrong instrument for a `‖·‖_∞`.
//! The step ladder is the **cell sizes of the resolution ladder**, so "the class"
//! means the class at the scales benched. Both a max reduction (theory-matched:
//! `q = ∞` is a sup norm) and a median reduction are fitted, and the median slope
//! is a column, because a max over probes is one probe's opinion. The two are
//! floored **independently** — on a flat surface the median difference is
//! identically zero while the max is a healthy `O(t)`.
//!
//! `f64` throughout: an order-4 difference at `t = 0.0625` is about `1e-4` for a
//! smooth field against a cancellation floor near `4e-14`, four orders of
//! headroom. Order 6 would not have had it, which is why
//! `MAX_DIFFERENCE_ORDER` stops at four and is itself a column.
//!
//! # What the two arms mean, and why C3 needs both
//!
//! `g_n` is a property of a **class**, never of a function: a single smooth field
//! belongs to `F_∞^s` for every `s`, and an algorithm that already knows it has
//! zero error. So "how far from optimal is Marching Cubes" has no answer until
//! somebody commits to a regularity class, and **the commitment is the answer**.
//! `pinned_w2` records C3's registered branch — against the class its own
//! trilinear model targets, `Θ(h²)` *is* the floor and only the constant is left.
//! `measured_max` records the other branch — against the smoother class the fields
//! measurably inhabit, the floor is `Θ(h^s)` with `s > 2` and Marching Cubes is
//! not order-optimal at all, which *licenses* `P-157` rather than capping it. The
//! registration says C3 *"is not falsifiable -- it is a branch, and both branches
//! are recorded"*, so both are rows, and vacuity control 7 refuses to let the file
//! be written with only one leg of `order_optimal` present.
//!
//! `c1_holds` and `c3_holds` are **global** verdicts and carry the same value on
//! every row: C1 is a statement about the literature, and C3 is answered by the
//! pair of arms existing. `c2_holds` is **per row** — `false` exactly when the
//! row's ratio is not trustworthy, which is the registration's own C2 falsifier
//! (*"a ratio that is not computable because the class assumptions are
//! unverifiable on procedural fields"*): either the log-log fit is worse than
//! [`FIT_R2_FLOOR`], or the measured regularity does not reach the arm's pinned
//! `s`.
//!
//! # Vacuity controls
//!
//! All seven run before the first `run.record`, and every panic starts `VOID: `.
//! Controls 5a and 6a fire inside [`measure`] because a bad ladder would otherwise
//! panic inside the fit with an uninformative message; the rest fire in `main`.
//!
//! 1. **The literature is transcribed faithfully.** The general-`d` rate must
//!    reproduce Proposition 2.4's independently stated `d = 1` form (p. 16), the
//!    `s = 1` form `n^{-1/d}` proved from scratch on p. 121, and Theorem 6.10's
//!    `p ≥ q` covering-radius exponent as literally printed (p. 93). Proves the
//!    formula this whole row is denominated in. Names `minimal_error_rate`.
//! 2. **Regularity exists and is finite for every field** — the registration's own
//!    control, verbatim. Names `regularity`.
//! 3. **The regularity estimator discriminates.** `lens_control` must
//!    measure at least [`REGULARITY_SEPARATION`] rougher than the **smoothest**
//!    field on the roster, or the estimator returns about the same order for a
//!    creased polyhedron as for a sphere and every class argument in the file is
//!    vacuous. Against the smoothest and not the roughest, deliberately:
//!    `box_exact` is genuinely creased along twelve edges and measures about as
//!    rough as the control, so the roughest-field form of this clause fires on a
//!    *correct* measurement. Names `regularity`.
//! 4. **The probes are on the surface, and measurably closer to it than a random
//!    point.** Two clauses: the worst probe's first-order distance
//!    `|f|/‖∇f‖` must be inside [`PROBE_BAND_CELLS`] cells, **and** must be below
//!    the median distance of a deterministic [`SCATTER_COUNT`]-point scatter over
//!    the field's own domain. The comparative clause is the one with teeth — it
//!    needs no constant chosen after seeing the data. A raw `|f|` bar was tried
//!    first and was wrong: `capped_gyroid` is `Lipschitz { l: 3.464 }` and
//!    `noise_cavity` is `Unbounded`, so a field value is not a length on this
//!    roster. Names `probe_max_distance` and `scatter_median_distance`.
//! 5. **The gradient does not vanish on the band** (`G_max > 0`), or the level-set
//!    reduction divides by zero and `constant_gap` measures nothing. Names
//!    `grad_max`. **5a**: the probe set is non-empty at all.
//! 6. **The instrument sees a falling error on at least
//!    [`MIN_FALLING_FIELDS`] fields.** Per field this is *data*, not a
//!    precondition, and `thin_plate` is why: it is
//!    `ThinPlate::for_cell_size(4/64)` with `THICKNESS_IN_CELLS = 0.4`, so its
//!    half-thickness is exactly `0.0125` and the plate is **sub-voxel at every
//!    resolution on this ladder by construction** (M-266, M-72). The mesh always
//!    bridges its two sheets inside one cell, an edge midpoint lands on the
//!    mid-plane where `|f|` is the half-thickness, and the error is pinned at
//!    `1.25e-2` independently of `h`. It cannot fall before `h < 0.025`, i.e.
//!    `N ≥ 161` — `161³ = 4 173 281` samples per pass against the `65³ = 274 625`
//!    this bench can afford. That exponent is **arithmetically unreachable** and is
//!    recorded as one, `error_falls = false` and `c2_holds = false`, which is the
//!    registration's own C2 falsifier verbatim. Names `error_falls`, `error_span`
//!    and `measured_error_rate`. **6a**: at least [`MIN_LADDER_POINTS`]
//!    resolutions meshed and every recorded error is strictly positive.
//! 7. **Both C3 branches are present.** `order_optimal` must take both values
//!    across the eighteen rows. Names `order_optimal`.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::fields::{FieldBound, ReferenceField};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, Sdf};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// The resolution ladder, samples per axis.
///
/// Five points rather than the house three: two registered columns are fitted
/// exponents, and a three-point fit reports no residual. `65` is the ceiling for
/// the two-minute budget.
const RESOLUTIONS: [u32; 5] = [17, 25, 33, 49, 65];

/// Which ladder point donates the probe set. Must be one of [`RESOLUTIONS`], and
/// is the last of them so the probes come from the finest surface available.
const REFERENCE_RESOLUTION: u32 = 65;

/// Ladder points a field must reach for its fit to be recorded.
const MIN_LADDER_POINTS: usize = 4;

/// How many of the nine fields must show a falling error for the instrument to be
/// believed. A majority: one deliberately sub-voxel fixture (`thin_plate`) has an
/// arithmetically unreachable exponent and says so in its own row, but a roster
/// that is flat everywhere means the functional is blind.
const MIN_FALLING_FIELDS: usize = 5;

// No cap on surface probes: the sup norm is taken over **every** vertex of the
// reference mesh, which is why no `PROBE_COUNT` constant appears here.
//
// An earlier version strided to 8192, and striding is the wrong instrument for a
// `‖·‖_∞`: whether a crease gets sampled then depends on the stride rather than on
// the field, so the measured `W_∞^s` seminorm moved with a constant nobody chose
// for a reason. The whole bench runs in about a second, so there is no case for
// approximating the maximum.

/// Points in vacuity control 4's negative-control scatter over the field's domain.
const SCATTER_COUNT: usize = 4096;

/// The scatter's seed. Stated because a control that moves between runs is not a
/// control. `common::poly::Rng` is the shared SplitMix64.
const SCATTER_SEED: u64 = 0x5EED_1234_A159;

/// How many cells wide the surface band may be before the probes are not on the
/// surface. A Marching Cubes vertex lies on a cut grid edge, so its true distance
/// to the zero set is at most one cell diagonal `h√3 ≈ 1.73 h`; the first-order
/// estimate overshoots that where the field curves inside a cell, which
/// `noise_cavity` is deliberately built to do (its features are `1/3.45 ≈ 0.29`
/// across, M-209). Four cells admits that and still sits far below a typical
/// domain point, which is why the comparative clause beside this one carries the
/// weight.
const PROBE_BAND_CELLS: f64 = 4.0;

/// Highest finite-difference order probed, and therefore the highest integer
/// smoothness this instrument can certify. That it is a **choice** rather than a
/// measurement is the finding, so it is a column.
const MAX_DIFFERENCE_ORDER: usize = 4;

/// The domain dimension `d` of Theorem 6.10.
const DOMAIN_DIMENSION: u32 = 3;

/// The `p` of `L_p`. A mesh error is a worst-case displacement, so `p = ∞`.
const LP_EXPONENT: f64 = f64::INFINITY;

/// The `q` of `W_q^s`. Bounded derivatives, so `q = ∞`, and Theorem 6.10's
/// embedding hypothesis `s > d/q` reduces to `s > 0`.
const LQ_EXPONENT: f64 = f64::INFINITY;

/// The `s` of the class Marching Cubes' trilinear model is second-order exact
/// for. The `pinned_w2` arm's whole content.
const TRILINEAR_MODEL_ORDER: u32 = 2;

/// Gradient magnitudes below this are treated as the floor rather than divided
/// by, and counted in `grad_floor_vertices`.
const GRAD_FLOOR: f64 = 1e-12;

/// A finite difference below this is cancellation noise, not a measurement: an
/// order-4 difference of an `O(10)` field cancels at about `4e-14`.
const DIFFERENCE_FLOOR: f64 = 1e-12;

/// Below this the log-log fit is not a rate, and `c2_holds` is `false`.
const FIT_R2_FLOOR: f64 = 0.90;

/// `|ratio − 1|` within this counts as order-optimal.
const ORDER_OPTIMAL_TOLERANCE: f64 = 0.15;

/// Slack on `regularity` when asking whether a field is in an arm's class.
const REGULARITY_TOLERANCE: f64 = 0.25;

/// How much rougher the control must measure than every reference field.
const REGULARITY_SEPARATION: f64 = 0.5;

/// Exact-arithmetic slack for vacuity control 1's transcription checks.
const TRANSCRIPTION_TOLERANCE: f64 = 1e-12;

/// Pascal's triangle, rows `0..=MAX_DIFFERENCE_ORDER`.
const BINOMIALS: [[f64; MAX_DIFFERENCE_ORDER + 1]; MAX_DIFFERENCE_ORDER + 1] = [
    [1.0, 0.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 2.0, 1.0, 0.0, 0.0],
    [1.0, 3.0, 3.0, 1.0, 0.0],
    [1.0, 4.0, 6.0, 4.0, 1.0],
];

// ---------------------------------------------------------------------------
// The literature, transcribed
// ---------------------------------------------------------------------------

/// The `n`-th minimal error exponent: `g_n(F_q^s, L_p) ≍ n^{-α}` on a bounded
/// convex `D ⊂ ℝ^d`.
///
/// Theorem 6.10 (Krieg & Ullrich p. 93) composed with Proposition 6.11 (p. 94):
/// the optimal point set has `‖dist(·,P_n)‖_{L_γ} ≍ n^{-1/d}` for every
/// `0 < γ ≤ ∞`, and Theorem 6.10 raises that to the power `s` when `p < q`, and to
/// `s − d(1/q − 1/p)` when `p ≥ q`. Dividing the second exponent by `d` cancels
/// the `d` in front of the bracket, which is the whole reason vacuity control 1
/// exists.
///
/// # Panics
///
/// If the Sobolev embedding hypothesis `s > d/q` fails, because outside it
/// `W_q^s(D)` is not a space of continuous functions and `f(x_i)` is not defined.
fn minimal_rate_exponent(s: u32, d: u32, p: f64, q: f64) -> f64 {
    assert!(s >= 1, "Theorem 6.10 needs s in N, got {s}");
    assert!(d >= 1, "Theorem 6.10 needs d in N, got {d}");
    assert!(
        f64::from(s) > f64::from(d) / q,
        "Theorem 6.10 needs s > d/q for the embedding into C(D): s={s} d={d} q={q}"
    );
    f64::from(s) / f64::from(d) - (1.0 / q - 1.0 / p).max(0.0)
}

/// Proposition 2.4 as printed on p. 16, independently: at `d = 1` and equidistant
/// points, `g_n(F_q^s, L_p) ≍ n^{-s + (1/q − 1/p)_+}`.
fn prop_2_4_exponent(s: u32, p: f64, q: f64) -> f64 {
    f64::from(s) - (1.0 / q - 1.0 / p).max(0.0)
}

/// Theorem 6.10's `p ≥ q` branch as printed on p. 93: the exponent **of the
/// covering radius**, `s − d(1/q − 1/p)`, before Proposition 6.11 turns the
/// covering radius into `n^{-1/d}`.
fn theorem_6_10_covering_exponent(s: u32, d: u32, p: f64, q: f64) -> f64 {
    f64::from(s) - f64::from(d) * (1.0 / q - 1.0 / p)
}

/// The three transcription residuals of vacuity control 1, worst of each kind.
///
/// 1. Against Proposition 2.4 (p. 16) at `d = 1`, over a grid of `(s, p, q)`.
/// 2. Against the `s = 1`, `p = q = ∞` case proved from scratch on p. 121, whose
///    rate is `n^{-1/d}`.
/// 3. Against Theorem 6.10's own `p ≥ q` covering-radius exponent (p. 93),
///    divided by `d`.
///
/// The grid skips cells where the paper's own hypotheses do not hold, which is
/// why the `s > d/q` test appears here as well as inside
/// [`minimal_rate_exponent`].
fn transcription_residuals() -> (f64, f64, f64) {
    let exponents = [1.0, 2.0, 4.0, f64::INFINITY];
    let mut prop_2_4 = 0.0_f64;
    let mut covering = 0.0_f64;

    for s in 1..=4_u32 {
        for &p in &exponents {
            for &q in &exponents {
                if f64::from(s) <= 1.0 / q {
                    continue;
                }
                let mine = minimal_rate_exponent(s, 1, p, q);
                prop_2_4 = prop_2_4.max((mine - prop_2_4_exponent(s, p, q)).abs());

                for d in 1..=5_u32 {
                    if p < q || f64::from(s) <= f64::from(d) / q {
                        continue;
                    }
                    let general = minimal_rate_exponent(s, d, p, q);
                    let printed = theorem_6_10_covering_exponent(s, d, p, q) / f64::from(d);
                    covering = covering.max((general - printed).abs());
                }
            }
        }
    }

    let mut lipschitz = 0.0_f64;
    for d in 1..=5_u32 {
        let mine = minimal_rate_exponent(1, d, f64::INFINITY, f64::INFINITY);
        lipschitz = lipschitz.max((mine - 1.0 / f64::from(d)).abs());
    }

    (prop_2_4, lipschitz, covering)
}

/// The explicit fooling-function floor on the **level-set** error at spacing `h`,
/// or `NaN` when the field's `s`-th seminorm is not measurable.
///
/// `f_0(x) = A ∏ sin(π x_a/h)` vanishes on all of `hℤ³` and has largest pure
/// `s`-th derivative `A (π/h)^s`; pinning that to `seminorm` gives
/// `A = seminorm · (h/π)^s`. `f` and `f ± f_0` share their grid data and their zero
/// sets are displaced by about `A/|∇f|`, so nothing reading only the grid can be
/// closer than `A/(2·grad_max)` to all three.
///
/// **A seminorm at the difference floor returns `NaN`, not zero, and
/// `box_exact` is why.** Its surface is flat, so an order-2 difference over its
/// own mesh vertices is *identically zero* at some steps: the measured
/// `|f|_{W_∞^2}` is `0`, the floor collapses to `0`, and `measured/floor` is `inf`
/// — a number that reads like "infinitely far from optimal" when what happened is
/// that `F_∞^2` is the wrong class for a piecewise-linear field. `NaN` says *not
/// computed*, and `c2_holds` is `false` on that row, which is the registration's
/// own C2 falsifier.
fn minimal_level_set_error(seminorm: f64, s: u32, h: f64, grad_max: f64) -> f64 {
    if seminorm.is_nan() || seminorm <= DIFFERENCE_FLOOR {
        return f64::NAN;
    }
    let amplitude = seminorm * (h / std::f64::consts::PI).powi(s as i32);
    amplitude / (2.0 * grad_max)
}

// ---------------------------------------------------------------------------
// The roughness control field
// ---------------------------------------------------------------------------

/// The intersection of two unit spheres centred at `(∓0.5, 0, 0)`: a lens with
/// **curved faces and one circular crease** at `x = 0`, `y² + z² = 0.75`.
///
/// The roughness control. It has to fail in two directions at once, and finding a
/// field that does took two attempts:
///
/// - Its **regularity** must measure near `1` where a sphere measures near `4`, or
///   the estimator is not responding to smoothness and every `field_class` in this
///   file is a class we merely hope applies (vacuity control 3). The `max(a, b)` of
///   two smooth branches is only Lipschitz across the crease, so every difference
///   order collapses to slope `1` there.
/// - Its **error** must still be a measurable number, and this is where the first
///   attempt died. `‖x‖_1 − 1` — the `L¹` unit ball — is a perfect regularity
///   control (it measured `s = 1.000000` exactly) and a useless error control,
///   because it is **linear inside every octant** and Marching Cubes reproduces a
///   piecewise-linear field *exactly*: its measured error was `~1e-14`, at the
///   `f64` floor, with no fittable exponent. A control whose own rows are noise
///   cannot demonstrate the thing being controlled for. The same mechanism is why
///   `box_exact`'s order-2 seminorm vanishes above.
///
/// The lens keeps the crease and adds curvature, so the error is a genuine
/// `Θ(h)` — `O(h²)` on the curved faces, `O(h)` at the chord across the crease —
/// and `order_optimal` can come out either way on its rows rather than always
/// `false`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Lens;

impl Lens {
    /// Half the distance between the two sphere centres. Below [`Self::RADIUS`],
    /// so the spheres overlap and the lens is non-empty.
    const OFFSET: f64 = 0.5;
    /// Both spheres' radius.
    const RADIUS: f64 = 1.0;
    /// The lens's centre, and it is **off-lattice on purpose**.
    ///
    /// The crease is the perpendicular bisector plane of the two centres, so it
    /// lies at `x = CENTRE[0]`. Centred at the origin the crease would sit at
    /// `x = 0`, which is a **grid plane at every resolution on this ladder** —
    /// `2/h` is `(N−1)/2`, an integer for every odd `N` — and this harness measured
    /// the consequence: with the crease exactly on a cell boundary each cell sees
    /// only one smooth branch, no triangle spans a chord across the kink, and the
    /// control's error came out `Θ(h²)` like a sphere's instead of the `Θ(h)` a
    /// sharp edge must give. It is `box_exact`'s alignment trap in a second
    /// costume, and `ThinPlate`'s own docs record the third (M-266).
    ///
    /// `0.11` puts the crease at `(2 + 0.11)/h` cells from the domain edge:
    /// `8.44`, `12.66`, `16.88`, `25.32`, `33.76` at `N = 17, 25, 33, 49, 65`. The
    /// closest approach to a grid plane is `0.12` of a cell, so the crease is
    /// strictly interior to a cell at every level. `y` and `z` are offset too, so
    /// the curved faces are not aligned either.
    const CENTRE: [f64; 3] = [0.11, 0.07, 0.13];
}

impl Sdf for Lens {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let ball = |offset: f64| {
            let d = [
                p[0] - (Self::CENTRE[0] + offset),
                p[1] - Self::CENTRE[1],
                p[2] - Self::CENTRE[2],
            ];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - Self::RADIUS
        };
        ball(-Self::OFFSET).max(ball(Self::OFFSET))
    }
}

impl ReferenceField for Lens {
    const NAME: &'static str = "lens_control";

    fn domain(&self) -> ([f64; 3], [f64; 3]) {
        ([-2.0; 3], [2.0; 3])
    }

    fn closed_in_domain(&self) -> bool {
        true
    }

    fn expected_euler(&self) -> Option<i64> {
        Some(2)
    }

    fn bound(&self) -> FieldBound {
        // Both branches are exact unit-gradient sphere distances and `max`
        // preserves a Lipschitz constant, so `l = 1`. It is *not* `Exact`: outside
        // the lens near the crease the true distance is to the crease circle, which
        // is shorter than either branch.
        FieldBound::Lipschitz { l: 1.0 }
    }
}

// ---------------------------------------------------------------------------
// Fitting
// ---------------------------------------------------------------------------

/// A least-squares fit of `ln y = intercept + slope · ln x`.
#[derive(Clone, Copy, Debug)]
struct Fit {
    slope: f64,
    r2: f64,
}

/// Fit `(x, y)` pairs in log-log. `r2` is the squared correlation, which for a
/// single-predictor least-squares line is the coefficient of determination.
///
/// # Panics
///
/// If fewer than two points are supplied, or any coordinate is not positive.
/// Every caller has already asserted both, so reaching either is a bug here
/// rather than a property of a field.
fn fit_log_log(points: &[(f64, f64)]) -> Fit {
    assert!(points.len() >= 2, "a log-log fit needs two points");
    let count = points.len() as f64;
    let logs: Vec<(f64, f64)> = points
        .iter()
        .map(|&(x, y)| {
            assert!(
                x > 0.0 && y > 0.0,
                "a log-log fit needs positive data, got ({x:e}, {y:e})"
            );
            (x.ln(), y.ln())
        })
        .collect();
    let xbar = logs.iter().map(|&(x, _)| x).sum::<f64>() / count;
    let ybar = logs.iter().map(|&(_, y)| y).sum::<f64>() / count;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    for &(x, y) in &logs {
        sxy += (x - xbar) * (y - ybar);
        sxx += (x - xbar) * (x - xbar);
        syy += (y - ybar) * (y - ybar);
    }
    let slope = sxy / sxx;
    let r2 = if syy > 0.0 {
        sxy * sxy / (sxx * syy)
    } else {
        0.0
    };
    Fit { slope, r2 }
}

// ---------------------------------------------------------------------------
// Measuring one field
// ---------------------------------------------------------------------------

/// The centred `order`-th finite difference of `field` along `axis` at step `t`.
///
/// `Δ^m_t f = Σ_j (−1)^{m−j} C(m,j) f(p + (j − m/2) t e_axis)`, which is exactly
/// `t^m f^{(m)}(ξ)` for some `ξ` in the stencil when `f ∈ C^m`. Left unscaled on
/// purpose: dividing by `t^m` is the seminorm estimate, and the caller does it.
///
/// The stencil is free to leave the field's declared `domain()`. Every field in
/// the roster is a closed-form expression defined on all of `ℝ³`, and clamping
/// would silently turn an `m`-th difference into something else.
fn nth_difference<S>(field: &S, p: [f64; 3], axis: usize, t: f64, order: usize) -> f64
where
    S: Sdf<Scalar = f64>,
{
    let centre = order as f64 / 2.0;
    let mut acc = 0.0;
    for (j, &binomial) in BINOMIALS[order][..=order].iter().enumerate() {
        let mut q = p;
        q[axis] += (j as f64 - centre) * t;
        let sign = if (order - j).is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        acc += sign * binomial * field.sample(q);
    }
    acc
}

/// Per-order regularity readings over one field's surface band.
#[derive(Clone, Copy, Debug)]
struct Roughness {
    /// Slope of `log max|Δ^m_t f|` against `log t`, per order: `min(m, s)`.
    slope_max: [f64; MAX_DIFFERENCE_ORDER + 1],
    /// The same slope from the median rather than the max reduction.
    slope_median: [f64; MAX_DIFFERENCE_ORDER + 1],
    /// `max|Δ^m_t f| / t^m` at the smallest step: an estimate of `‖f^{(m)}‖_∞`
    /// over pure axis derivatives only, hence an underestimate of `|f|_{W_∞^m}`.
    seminorm: [f64; MAX_DIFFERENCE_ORDER + 1],
    /// Set when some step drove the **max** reduction of an order below
    /// [`DIFFERENCE_FLOOR`], so that order's `slope_max` is cancellation noise and
    /// is reported as `NaN`.
    below_floor: [bool; MAX_DIFFERENCE_ORDER + 1],
    /// The same, independently, for the median reduction. Separate because the two
    /// genuinely disagree on a piecewise-linear surface — see [`roughness`].
    below_floor_median: [bool; MAX_DIFFERENCE_ORDER + 1],
}

/// Fit the difference ladder over `probes` at every step in `steps`.
fn roughness<S>(field: &S, probes: &[[f64; 3]], steps: &[f64]) -> Roughness
where
    S: Sdf<Scalar = f64>,
{
    let mut out = Roughness {
        slope_max: [f64::NAN; MAX_DIFFERENCE_ORDER + 1],
        slope_median: [f64::NAN; MAX_DIFFERENCE_ORDER + 1],
        seminorm: [f64::NAN; MAX_DIFFERENCE_ORDER + 1],
        below_floor: [false; MAX_DIFFERENCE_ORDER + 1],
        below_floor_median: [false; MAX_DIFFERENCE_ORDER + 1],
    };
    let mut buffer: Vec<f64> = Vec::with_capacity(probes.len() * 3);

    for order in 1..=MAX_DIFFERENCE_ORDER {
        let mut maxima: Vec<(f64, f64)> = Vec::with_capacity(steps.len());
        let mut medians: Vec<(f64, f64)> = Vec::with_capacity(steps.len());
        // The two reductions are floored **independently**, and coupling them was
        // a real bug this harness hit on `box_exact`: its probes sit on box faces
        // where the field is exactly linear, so the *median* second difference is
        // identically `0.0` while the *max* — driven by the stencils that straddle
        // a box edge — is a healthy `O(t)`. One shared flag threw away the
        // theory-matched reading because the robustness readout beside it was
        // degenerate.
        let mut floored_max = false;
        let mut floored_median = false;

        for &t in steps {
            buffer.clear();
            for &p in probes {
                for axis in 0..3 {
                    buffer.push(nth_difference(field, p, axis, t, order).abs());
                }
            }
            let worst = buffer.iter().copied().fold(0.0_f64, f64::max);
            buffer.sort_unstable_by(f64::total_cmp);
            let middle = buffer[buffer.len() / 2];

            floored_max |= worst < DIFFERENCE_FLOOR;
            floored_median |= middle < DIFFERENCE_FLOOR;
            maxima.push((t, worst));
            medians.push((t, middle));
        }

        // The seminorm is read at the smallest step, found by value rather than by
        // position so it does not depend on how `steps` happens to be ordered.
        let finest = maxima
            .iter()
            .copied()
            .reduce(|a, b| if b.0 < a.0 { b } else { a })
            .expect("the step ladder is non-empty");
        out.seminorm[order] = finest.1 / finest.0.powi(order as i32);

        out.below_floor[order] = floored_max;
        out.below_floor_median[order] = floored_median;
        if !floored_max {
            out.slope_max[order] = fit_log_log(&maxima).slope;
        }
        if !floored_median {
            out.slope_median[order] = fit_log_log(&medians).slope;
        }
    }

    out
}

/// Everything measured about one field, before any class is pinned.
#[derive(Clone, Debug)]
struct Measured {
    name: &'static str,
    is_control: bool,
    /// `(samples, cell_size, level-set error)` for each ladder point that meshed.
    ladder: Vec<(u32, f64, f64)>,
    /// The error against the function-value count `n = samples³` — IBC's currency.
    fit_n: Fit,
    /// The same error against the cell size `h`, for readability only.
    fit_h: Fit,
    rough: Roughness,
    grad_max: f64,
    probe_count: usize,
    /// The largest raw `|f|` over the probes. Informative only: on a non-unit
    /// gradient field it is not a length.
    probe_max_abs_value: f64,
    /// The largest first-order **distance** `|f|/‖∇f‖` over the probes. This is
    /// what vacuity control 4 reads.
    probe_max_distance: f64,
    /// The same statistic's median over a deterministic scatter of the field's own
    /// domain: control 4's negative control.
    scatter_median_distance: f64,
    /// The cell size the probes were taken at.
    probe_cell_size: f64,
    grad_floor_vertices: u64,
    vertices_finest: usize,
    triangles_finest: usize,
    wall_seconds: f64,
}

impl Measured {
    /// The headline regularity: the order-[`MAX_DIFFERENCE_ORDER`] slope, which is
    /// `min(MAX_DIFFERENCE_ORDER, s)`.
    fn regularity(&self) -> f64 {
        self.rough.slope_max[MAX_DIFFERENCE_ORDER]
    }

    /// The finest ladder point, where `constant_gap` is evaluated.
    fn finest(&self) -> (u32, f64, f64) {
        *self.ladder.last().expect("the ladder is non-empty")
    }

    /// The coarsest ladder point, so a growing `constant_gap` is visible.
    fn coarsest(&self) -> (u32, f64, f64) {
        *self.ladder.first().expect("the ladder is non-empty")
    }

    /// Whether the error fell at all across the ladder.
    ///
    /// `false` is a **result**, not a failure: it means this field's exponent is
    /// not measurable at the resolutions benched, which is the registration's own
    /// C2 falsifier. `thin_plate` is the field it fires on.
    fn error_falls(&self) -> bool {
        self.finest().2 < self.coarsest().2
    }

    /// `error_coarsest / error_finest`. `1.0` is flat; `h²` over this ladder is
    /// about `(64/16)² = 16`. One number that says whether an exponent exists.
    fn error_span(&self) -> f64 {
        self.coarsest().2 / self.finest().2
    }
}

/// The first-order distance from `x` to `{f = 0}`, and whether the gradient there
/// sat at the floor.
///
/// `|f(x)| / ‖∇f(x)‖`. Dividing by the gradient is what makes the quantity a
/// **length** rather than a field value, and that matters on this roster: only
/// four of the eight reference fields are unit-gradient. `capped_gyroid` is
/// `Lipschitz { l: 3.464 }` and `noise_cavity` is `Unbounded`, so comparing a raw
/// `|f|` against a cell size asks a scale-dependent question — a bug this harness
/// hit in its own vacuity control 4 before the normalisation was added.
fn surface_distance<S>(field: &S, x: [f64; 3]) -> (f64, bool)
where
    S: Sdf<Scalar = f64>,
{
    let g = field.gradient(x);
    let norm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    (
        field.sample(x).abs() / norm.max(GRAD_FLOOR),
        norm < GRAD_FLOOR,
    )
}

/// The extractor's level-set error over one mesh, and how many of its sampled
/// points sat at the gradient floor.
///
/// `|f(x)| / ‖∇f(x)‖` is, to first order, the distance from `x` to the true zero
/// set, and the functional is its **maximum over the mesh** — every vertex, every
/// edge midpoint and every triangle centroid. That is the one-sided Hausdorff
/// distance `sup over {x in mesh} dist(x, {f = 0})`, discretised on the mesh, and
/// it is what `validate::accuracy`'s `mesh_to_field` direction estimates.
///
/// **Sampling the vertices alone is wrong, and `box_exact` is the proof.** Its
/// half-extents are `1` and every ladder spacing `h = 4/(N−1)` divides `3`, so
/// every box face lands on a grid plane, every box edge on a grid line, and the
/// field is *exactly linear* along every cut grid edge — Marching Cubes finds
/// those roots exactly and a vertex-only functional reads **identically zero at
/// all five resolutions**. It is not zero: `M-250` measured `box_exact`'s
/// symmetric Hausdorff at `1.443e-1` before and after crossing refinement,
/// because the mesh *truncates* every box edge with a chord across it. The chord's
/// midpoint is about `h/2` inside the solid, and only a functional that looks
/// between the vertices can see it. Adding midpoints and centroids costs `O(T)`
/// samples, leaves the `O(h²)` rate of the smooth fields alone (a chord across
/// curvature `κ` deviates by `O((hκ)²)`), and turns `box_exact` from an
/// unfittable zero into the `Θ(h)` its own measured regularity predicts.
///
/// It is used in preference to `validate::accuracy` for one reason, and it is not
/// a small one: `accuracy` reads the SDF as a distance oracle, and only four of
/// the eight reference fields have `bound() == Exact` — `gyroid` is `Lipschitz`,
/// `csg_difference` is `Underestimate { q: 0.5 }`, `fbm_terrain` and
/// `noise_cavity` are `Unbounded`. A functional that is meaningless on half the
/// roster cannot carry a rate fit across it. The gradient-normalised residual
/// needs only a non-vanishing gradient, which is the same transversality
/// hypothesis `P-171` measures.
///
/// Shared edges are sampled twice and that is deliberate: a maximum is
/// idempotent, so deduplicating them would buy nothing but a hash set.
fn level_set_error<S>(field: &S, mesh: &MeshBuffer<f64>) -> (f64, u64)
where
    S: Sdf<Scalar = f64>,
{
    let mut worst = 0.0_f64;
    let mut at_floor = 0_u64;
    let probe = |x: [f64; 3], worst: &mut f64, at_floor: &mut u64| {
        let (distance, floored) = surface_distance(field, x);
        if floored {
            *at_floor += 1;
        }
        *worst = worst.max(distance);
    };

    for &v in &mesh.positions {
        probe(v, &mut worst, &mut at_floor);
    }
    for tri in mesh.indices.as_chunks::<3>().0 {
        let a = mesh.positions[tri[0] as usize];
        let b = mesh.positions[tri[1] as usize];
        let c = mesh.positions[tri[2] as usize];
        for (p, q) in [(a, b), (b, c), (c, a)] {
            probe(
                [
                    (p[0] + q[0]) / 2.0,
                    (p[1] + q[1]) / 2.0,
                    (p[2] + q[2]) / 2.0,
                ],
                &mut worst,
                &mut at_floor,
            );
        }
        probe(
            [
                (a[0] + b[0] + c[0]) / 3.0,
                (a[1] + b[1] + c[1]) / 3.0,
                (a[2] + b[2] + c[2]) / 3.0,
            ],
            &mut worst,
            &mut at_floor,
        );
    }
    (worst, at_floor)
}

/// Walk the resolution ladder, then the difference ladder on the finest mesh's own
/// vertices.
///
/// # Panics
///
/// Vacuity controls 5a and 6a: `VOID: ` if the reference resolution produced no
/// probes, if fewer than [`MIN_LADDER_POINTS`] resolutions meshed, or if any
/// recorded error is not strictly positive. All three would otherwise surface as
/// an uninformative panic inside [`fit_log_log`].
fn measure<F>(name: &'static str, field: &F, is_control: bool) -> Measured
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let started = Instant::now();
    let mut extractor = MarchingCubes::<f64>::new();
    let mut mesh = MeshBuffer::<f64>::new();

    let mut ladder: Vec<(u32, f64, f64)> = Vec::with_capacity(RESOLUTIONS.len());
    let mut grad_floor_vertices = 0_u64;
    let mut probes: Vec<[f64; 3]> = Vec::new();
    let mut probe_cell_size = 0.0_f64;
    let mut vertices_finest = 0;
    let mut triangles_finest = 0;

    for samples in RESOLUTIONS {
        let (shape, origin, cell_size) = common::grid::<f64, _>(field, samples);
        mesh.reset();
        extractor
            .extract(field, &shape, origin, cell_size, &mut mesh)
            .expect("every ladder grid has at least two samples per axis");

        if mesh.indices.is_empty() {
            continue;
        }
        let (error, at_floor) = level_set_error(field, &mesh);
        grad_floor_vertices += at_floor;
        ladder.push((samples, cell_size, error));

        if samples == REFERENCE_RESOLUTION {
            probe_cell_size = cell_size;
            vertices_finest = mesh.vertex_count();
            triangles_finest = mesh.triangle_count();
            probes = mesh.positions.clone();
        }
    }

    assert!(
        !probes.is_empty() && probe_cell_size > 0.0,
        "VOID: {name} produced no mesh at the reference resolution \
         {REFERENCE_RESOLUTION}, so there is no surface band to measure regularity on and its \
         class membership could only be assumed"
    );
    assert!(
        ladder.len() >= MIN_LADDER_POINTS,
        "VOID: {name} meshed at only {} of {} ladder points, fewer than the \
         {MIN_LADDER_POINTS} a reportable exponent fit needs",
        ladder.len(),
        RESOLUTIONS.len()
    );
    assert!(
        ladder.iter().all(|&(_, _, e)| e > 0.0),
        "VOID: {name} records a level-set error of exactly zero at some resolution, so the \
         extractor reproduces this field inside its own model space, there is no exponent to \
         fit and there is nothing left to denominate"
    );

    // The step ladder is the resolution ladder's own cell sizes: "the class" here
    // means the class at the scales the extractor was benched at, never a limit.
    let steps: Vec<f64> = ladder.iter().map(|&(_, h, _)| h).collect();
    let rough = roughness(field, &probes, &steps);

    let mut grad_max = 0.0_f64;
    let mut probe_max_abs_value = 0.0_f64;
    let mut probe_max_distance = 0.0_f64;
    for &p in &probes {
        let g = field.gradient(p);
        grad_max = grad_max.max((g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt());
        probe_max_abs_value = probe_max_abs_value.max(field.sample(p).abs());
        probe_max_distance = probe_max_distance.max(surface_distance(field, p).0);
    }

    // The negative control for vacuity control 4: the same statistic over a
    // deterministic scatter of the field's own domain. A probe set that is *not*
    // on the surface cannot be much closer to it than a typical domain point is,
    // so `probe_max_distance < scatter_median_distance` is a comparison the
    // instrument cannot satisfy by accident — and unlike an absolute bar in cells,
    // it needs no constant chosen after seeing the data. `common::poly::Rng` is
    // the shared SplitMix64; writing a second generator here would be a second
    // path to one answer.
    let (lo, hi) = field.domain();
    let mut rng = common::poly::Rng::new(SCATTER_SEED);
    let mut scatter: Vec<f64> = Vec::with_capacity(SCATTER_COUNT);
    for _ in 0..SCATTER_COUNT {
        let x = [
            lo[0] + (hi[0] - lo[0]) * rng.next_f64_unit(),
            lo[1] + (hi[1] - lo[1]) * rng.next_f64_unit(),
            lo[2] + (hi[2] - lo[2]) * rng.next_f64_unit(),
        ];
        scatter.push(surface_distance(field, x).0);
    }
    scatter.sort_unstable_by(f64::total_cmp);
    let scatter_median_distance = scatter[scatter.len() / 2];

    let against_n: Vec<(f64, f64)> = ladder
        .iter()
        .map(|&(samples, _, e)| (f64::from(samples).powi(3), e))
        .collect();
    let against_h: Vec<(f64, f64)> = ladder.iter().map(|&(_, h, e)| (h, e)).collect();

    Measured {
        name,
        is_control,
        fit_n: fit_log_log(&against_n),
        fit_h: fit_log_log(&against_h),
        ladder,
        rough,
        grad_max,
        probe_count: probes.len(),
        probe_max_abs_value,
        probe_max_distance,
        scatter_median_distance,
        probe_cell_size,
        grad_floor_vertices,
        vertices_finest,
        triangles_finest,
        wall_seconds: started.elapsed().as_secs_f64(),
    }
}

// ---------------------------------------------------------------------------
// Rows
// ---------------------------------------------------------------------------

/// One recorded row: one field, one pinned class.
#[derive(Clone, Copy, Debug)]
struct Row {
    field: &'static str,
    arm: &'static str,
    is_control: bool,
    s_pinned: u32,
    in_class: bool,
    minimal_rate: f64,
    measured_rate: f64,
    ratio: f64,
    order_optimal: bool,
    minimal_error_finest: f64,
    constant_gap: f64,
    constant_gap_coarsest: f64,
    c2_holds: bool,
}

/// Build both arms for one field.
///
/// `pinned_w2` fixes the class at [`TRILINEAR_MODEL_ORDER`] — the class C3 names.
/// `measured_max` fixes it at the largest integer smoothness the field's own
/// surface measurably supports, which is what the registration's vacuity control
/// demands the class be argued from.
fn rows_for(m: &Measured) -> [Row; 2] {
    let s_hat = m.regularity();
    // Theorem 6.10 needs `s in N`, so the measured order is rounded to an integer
    // and the rounding is visible in `regularity` beside it.
    let s_measured = (s_hat.round().max(1.0) as u32).min(MAX_DIFFERENCE_ORDER as u32);
    [
        row_for(m, "pinned_w2", TRILINEAR_MODEL_ORDER, s_hat),
        row_for(m, "measured_max", s_measured, s_hat),
    ]
}

/// One arm.
fn row_for(m: &Measured, arm: &'static str, s_pinned: u32, s_hat: f64) -> Row {
    let minimal_rate = minimal_rate_exponent(s_pinned, DOMAIN_DIMENSION, LP_EXPONENT, LQ_EXPONENT);
    let measured_rate = -m.fit_n.slope;
    let ratio = measured_rate / minimal_rate;

    let (_, h_fine, error_fine) = m.finest();
    let (_, h_coarse, error_coarse) = m.coarsest();
    let seminorm = m.rough.seminorm[s_pinned as usize];
    let floor_fine = minimal_level_set_error(seminorm, s_pinned, h_fine, m.grad_max);
    let floor_coarse = minimal_level_set_error(seminorm, s_pinned, h_coarse, m.grad_max);

    let in_class = s_hat >= f64::from(s_pinned) - REGULARITY_TOLERANCE;

    Row {
        field: m.name,
        arm,
        is_control: m.is_control,
        s_pinned,
        in_class,
        minimal_rate,
        measured_rate,
        ratio,
        order_optimal: (ratio - 1.0).abs() <= ORDER_OPTIMAL_TOLERANCE,
        minimal_error_finest: floor_fine,
        constant_gap: error_fine / floor_fine,
        constant_gap_coarsest: error_coarse / floor_coarse,
        // C2's registered falsifier, per row: the ratio is not trustworthy when the
        // error does not fall at all (`thin_plate`, sub-voxel by construction), when
        // the fit is not a rate, when the field is not in the pinned class, or when
        // the field's own seminorm for that class is unmeasurable (`box_exact` at
        // `s = 2`) so the floor -- and the gap against it -- does not exist.
        c2_holds: m.error_falls()
            && m.fit_n.r2 >= FIT_R2_FLOOR
            && in_class
            && ratio.is_finite()
            && (error_fine / floor_fine).is_finite(),
    }
}

/// The class token, one CSV-safe word.
fn class_token(s_pinned: u32) -> String {
    format!("W_inf^{s_pinned}_unit_ball_on_bounded_convex_D_in_R3")
}

/// Every hypothesis the comparison rests on, and every one it could not check.
///
/// This is C1's *"written down with all its hypotheses"* discharged into the data
/// rather than only into this file's prose, and it is where the three unverifiable
/// items from the header live. Pipe-separated because the CSV writer refuses
/// commas.
fn class_assumptions(row: &Row, s_hat: f64) -> String {
    let s = row.s_pinned;
    let in_class = row.in_class;
    format!(
        "source=KriegUllrich_arXiv2602.02066|thm6.10_p93|prop6.11_p94|prop2.4_p16|gn_def_p13\
         |crosscheck=BonitoCanutoNochettoVeeser_10.1017-s0962492924000011_p67\
         |class=W_inf^{s}_unit_ball|s_integer=true|s_gt_d_over_q=true|d=3|p=inf|q=inf\
         |domain=bounded_convex|information=n_function_values|adaptive=no_by_thm9.3\
         |UNVERIFIED_implied_constants_absent_from_source\
         |UNVERIFIED_error_is_Lp_of_f_not_of_level_set\
         |level_set_reduction=derived_here_not_in_source\
         |explicit_floor=own_fooling_function_sin_product\
         |seminorm=pure_axis_derivatives_only_so_underestimated\
         |regularity_probed_to_order={MAX_DIFFERENCE_ORDER}|regularity_at_benched_h_only\
         |s_measured={s_hat:.4}|field_in_class={in_class}"
    )
}

/// A wide-range number the CSV writer will accept.
fn sci(v: f64) -> String {
    format!("{v:.6e}")
}

/// A bounded number — a rate, a ratio, an exponent.
fn dec(v: f64) -> String {
    format!("{v:.6}")
}

/// Pipe-join, because the CSV writer refuses commas.
fn joined(parts: &[String]) -> String {
    parts.join("|")
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-159");

    common::experiment::run(prereg, |run| {
        // ---- measure everything before asserting anything -------------------
        let mut measured: Vec<Measured> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            // The macro inlines this body once per field, so no `return`, no
            // `break` and no `continue` may appear in it (M-199, M-253).
            measured.push(measure(name, &field, false));
        });
        measured.push(measure(Lens::NAME, &Lens, true));

        // ---- vacuity control 1: the literature is transcribed faithfully ----
        let (prop_2_4, lipschitz, covering) = transcription_residuals();
        assert!(
            prop_2_4 <= TRANSCRIPTION_TOLERANCE,
            "VOID: the general-d rate disagrees with Krieg & Ullrich Proposition 2.4 (p. 16) at \
             d=1 by {prop_2_4:e}, so `minimal_error_rate` is not the quantity the source states \
             and every ratio in this file is denominated in a formula nobody wrote"
        );
        assert!(
            lipschitz <= TRANSCRIPTION_TOLERANCE,
            "VOID: the general-d rate disagrees with the s=1 p=q=inf case proved from scratch on \
             p. 121 -- e(Q_n, F_Lip^d) >= (1/8) n^(-1/d) -- by {lipschitz:e}, so the d is in the \
             wrong place"
        );
        assert!(
            covering <= TRANSCRIPTION_TOLERANCE,
            "VOID: the general-d rate disagrees with Theorem 6.10's own p>=q covering-radius \
             exponent s - d(1/q - 1/p) as printed on p. 93 by {covering:e}"
        );

        // ---- per-field vacuity controls 2, 4, 5, 6 --------------------------
        for m in &measured {
            let s_hat = m.regularity();
            assert!(
                s_hat.is_finite() && s_hat > 0.0 && s_hat <= MAX_DIFFERENCE_ORDER as f64 + 0.5,
                "VOID: {}'s regularity estimate is {s_hat}, not a finite positive order, so its \
                 class membership would be assumed rather than argued from a measured property -- \
                 which is this registration's own vacuity control verbatim. below_floor={:?}",
                m.name,
                m.rough.below_floor
            );
            assert!(
                m.probe_max_distance <= PROBE_BAND_CELLS * m.probe_cell_size,
                "VOID: {}'s {} surface probes reach {:e} away from the zero set against the cell \
                 size {:e} they were taken at -- over {PROBE_BAND_CELLS} cells, so regularity was \
                 measured somewhere the extractor never looks and the class argument is about a \
                 different set",
                m.name,
                m.probe_count,
                m.probe_max_distance,
                m.probe_cell_size
            );
            assert!(
                m.probe_max_distance < m.scatter_median_distance,
                "VOID: {}'s worst surface probe is {:e} from the zero set while a *typical point \
                 of its domain* is only {:e} away, so the probe set is no closer to the surface \
                 than a random scatter and nothing licenses calling it a surface band",
                m.name,
                m.probe_max_distance,
                m.scatter_median_distance
            );
            assert!(
                m.grad_max > GRAD_FLOOR,
                "VOID: {}'s gradient maximum over the surface band is {:e}, so the level-set \
                 reduction from the L_inf function bound divides by zero and `constant_gap` \
                 measures nothing",
                m.name,
                m.grad_max
            );
            // Whether *this* field's error falls is data, not a precondition: see
            // the roster-level clause below.
        }

        // ---- vacuity control 6: the instrument sees a falling error at all ----
        // Per field this is *data*, and this harness learned that the hard way.
        // `thin_plate` is `ThinPlate::for_cell_size(4/64)` with
        // `THICKNESS_IN_CELLS = 0.4`, so its half-thickness is exactly
        // `0.0625 * 0.4 * 0.5 = 0.0125` and the plate is **sub-voxel at every
        // resolution on this ladder by construction** (M-266, M-72). The mesh
        // therefore always bridges its two sheets inside one cell, an edge midpoint
        // lands on the plate's mid-plane where `|f|` is the half-thickness, and the
        // error is pinned at `1.25e-2` independently of `h`. It cannot begin to
        // fall before `h < 0.025`, i.e. `N >= 161` -- `161^3 = 4_173_281` samples
        // per pass against the `65^3 = 274_625` this bench can afford, fifteen
        // times over. That is an **arithmetically unreachable** exponent and it is
        // recorded as one, with `error_falls = false` and `c2_holds = false`, which
        // is this registration's own C2 falsifier verbatim: *a ratio that is not
        // computable*. Aborting the whole file over it would have dropped seventeen
        // other rows to report one field's arithmetic.
        //
        // What must still hold is that the instrument is not uniformly blind.
        let falling = measured.iter().filter(|m| m.error_falls()).count();
        assert!(
            falling >= MIN_FALLING_FIELDS,
            "VOID: the error falls on only {falling} of {} fields, fewer than the \
             {MIN_FALLING_FIELDS} a working instrument must show -- so `measured_error_rate` is a \
             slope through noise across the roster rather than on one deliberately sub-voxel \
             fixture, and no ratio in this file means anything",
            measured.len()
        );

        // ---- vacuity control 3: the estimator discriminates -------------------
        // Against the **smoothest** field, not the roughest. Comparing the control
        // against the roughest reference field was this harness's second real bug:
        // `box_exact` is genuinely creased along twelve edges and measures about as
        // rough as the lens control, so that form fired on a *correct* measurement.
        // What has to be shown is that the estimator RESPONDS to smoothness — that
        // something on the roster reads well above the creased control.
        let control = measured
            .iter()
            .find(|m| m.is_control)
            .expect("the roughness control is in the roster");
        let smoothest = measured
            .iter()
            .map(Measured::regularity)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            control.regularity() + REGULARITY_SEPARATION <= smoothest,
            "VOID: the roughness control `{}` measures regularity {:.4} against the smoothest \
             field on the roster at {smoothest:.4}, a spread below {REGULARITY_SEPARATION} -- so \
             the estimator returns about the same order for a creased polyhedron as for a smooth \
             field, and every `field_class` in this file is a class we merely hope applies",
            control.name,
            control.regularity()
        );

        // ---- build the rows, then control 7 before recording any of them -----
        let rows: Vec<(Row, f64)> = measured
            .iter()
            .flat_map(|m| rows_for(m).map(|r| (r, m.regularity())))
            .collect();

        assert!(
            rows.iter().any(|(r, _)| r.order_optimal) && rows.iter().any(|(r, _)| !r.order_optimal),
            "VOID: `order_optimal` is {} on all {} rows, so only one leg of C3 was recorded -- and \
             C3's registration is explicitly a branch whose *both* branches are recorded. A file \
             with one leg cannot tell a reader whether Marching Cubes is order-optimal or whether \
             the class was chosen to make it so",
            rows[0].0.order_optimal,
            rows.len()
        );

        // C1: the n-th minimal error rate for a stated class is established from a
        // source that was read, its hypotheses are written into every row, and the
        // transcription is checked against two independently stated special cases
        // in the same paper. Vacuity control 1 is what earns this.
        let c1_holds = true;
        // C3: not falsifiable, a branch. Answered by both legs of `order_optimal`
        // being present, which vacuity control 7 is what earns.
        let c3_holds = true;

        println!(
            "{:<20} {:<13} {:>3} {:>9} {:>10} {:>7} {:>6} {:>11} {:>6}",
            "field", "arm", "s", "min_rate", "meas_rate", "ratio", "opt", "const_gap", "c2"
        );

        for (row, s_hat) in &rows {
            println!(
                "{:<20} {:<13} {:>3} {:>9.4} {:>10.4} {:>7.3} {:>6} {:>11.3e} {:>6}",
                row.field,
                row.arm,
                row.s_pinned,
                row.minimal_rate,
                row.measured_rate,
                row.ratio,
                row.order_optimal,
                row.constant_gap,
                row.c2_holds
            );

            let m = measured
                .iter()
                .find(|m| m.name == row.field)
                .expect("every row came from a measured field");
            let (_, h_fine, error_fine) = m.finest();
            let resolutions: Vec<String> =
                m.ladder.iter().map(|&(s, _, _)| s.to_string()).collect();
            let ladder_n: Vec<String> = m
                .ladder
                .iter()
                .map(|&(s, _, _)| u64::from(s).pow(3).to_string())
                .collect();
            let cell_sizes: Vec<String> = m.ladder.iter().map(|&(_, h, _)| sci(h)).collect();
            let errors: Vec<String> = m.ladder.iter().map(|&(_, _, e)| sci(e)).collect();

            run.record(&[
                ("field_class", class_token(row.s_pinned)),
                ("regularity", dec(*s_hat)),
                (
                    "n_samples",
                    u64::from(REFERENCE_RESOLUTION).pow(3).to_string(),
                ),
                ("minimal_error_rate", dec(row.minimal_rate)),
                ("measured_error_rate", dec(row.measured_rate)),
                ("ratio", dec(row.ratio)),
                ("order_optimal", row.order_optimal.to_string()),
                ("constant_gap", sci(row.constant_gap)),
                ("class_assumptions", class_assumptions(row, *s_hat)),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", row.c2_holds.to_string()),
                ("c3_holds", c3_holds.to_string()),
                // ── extras (M-273) ──
                ("arm", row.arm.to_string()),
                ("cell_sizes", joined(&cell_sizes)),
                ("constant_gap_coarsest", sci(row.constant_gap_coarsest)),
                (
                    "difference_ladder_below_floor",
                    format!("{:?}", m.rough.below_floor).replace(", ", "|"),
                ),
                (
                    "difference_ladder_below_floor_median",
                    format!("{:?}", m.rough.below_floor_median).replace(", ", "|"),
                ),
                ("error_coarsest", sci(m.coarsest().2)),
                ("error_falls", m.error_falls().to_string()),
                ("error_finest", sci(error_fine)),
                ("error_span", sci(m.error_span())),
                (
                    "error_functional",
                    String::from("max_abs_f_over_grad_norm_on_vertices_edge_midpoints_centroids"),
                ),
                ("errors", joined(&errors)),
                (
                    "extractor",
                    String::from("marching_cubes_defaults_refinement_0"),
                ),
                ("field", row.field.to_string()),
                ("field_in_pinned_class", row.in_class.to_string()),
                ("fit_r2", dec(m.fit_n.r2)),
                ("fit_r2_h", dec(m.fit_h.r2)),
                ("grad_floor_vertices", m.grad_floor_vertices.to_string()),
                ("grad_max", sci(m.grad_max)),
                ("h_finest", sci(h_fine)),
                ("is_control", row.is_control.to_string()),
                ("ladder_points", m.ladder.len().to_string()),
                ("max_difference_order", MAX_DIFFERENCE_ORDER.to_string()),
                // `err ∝ h^{+p}` while `err ∝ n^{-α}`, so this slope is NOT negated
                // the way `measured_error_rate` is. Directly comparable to
                // `minimal_h_exponent`, which is `s`.
                ("measured_h_exponent", dec(m.fit_h.slope)),
                ("minimal_error_finest", sci(row.minimal_error_finest)),
                ("minimal_h_exponent", dec(f64::from(row.s_pinned))),
                ("n_samples_ladder", joined(&ladder_n)),
                ("probe_cell_size", sci(m.probe_cell_size)),
                ("probe_count", m.probe_count.to_string()),
                ("probe_max_abs_value", sci(m.probe_max_abs_value)),
                ("probe_max_distance", sci(m.probe_max_distance)),
                ("scatter_median_distance", sci(m.scatter_median_distance)),
                (
                    "regularity_median",
                    dec(m.rough.slope_median[MAX_DIFFERENCE_ORDER]),
                ),
                ("resolutions", joined(&resolutions)),
                ("s_pinned", row.s_pinned.to_string()),
                (
                    "seminorm_pinned",
                    sci(m.rough.seminorm[row.s_pinned as usize]),
                ),
                ("slope_order1", dec(m.rough.slope_max[1])),
                ("slope_order2", dec(m.rough.slope_max[2])),
                ("slope_order3", dec(m.rough.slope_max[3])),
                ("slope_order4", dec(m.rough.slope_max[4])),
                ("triangles_finest", m.triangles_finest.to_string()),
                ("vertices_finest", m.vertices_finest.to_string()),
                ("wall_seconds", dec(m.wall_seconds)),
            ]);
        }
    });
}
