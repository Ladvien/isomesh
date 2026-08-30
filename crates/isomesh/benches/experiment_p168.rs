//! **P-168 — noise stability as a robustness number, tested against the wrong
//! noise model on purpose.**
//!
//! Ticket: R-168. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p168
//! ```
//!
//! Writes `docs/experiments/p-168.csv`.
//!
//! # What was missing
//!
//! `P-167` took the first Walsh–Hadamard transform this repository has ever
//! taken, of the 256-case table read as a Boolean function of eight corner signs
//! (`crates/isomesh/src/marching_cubes/table.rs:180`, `corner_inside` at
//! `cube.rs`). Its result is committed at `docs/experiments/p-167.csv` and this
//! row is built on it rather than on a re-derivation:
//!
//! - **The spectrum is dense.** `triangle_counts.bit1` and `.bit2` each have
//!   `spectral_terms = 256` — every one of the 256 Fourier coefficients is
//!   non-negligible — with `spectral_concentration` (weight at or below degree 2)
//!   of `0.243408203` and `0.428710938`. `bit0` is exactly
//!   [`Bool8::parity`](common::boolean::Bool8::parity): one coefficient, all of
//!   it at degree 8. `bit3` is the constant zero. `edge_masks.bit0` is a single
//!   coefficient at degree 2.
//! - **`total_influence` is `8.0` / `4.0` / `3.5` / `0.0` for the four triangle
//!   count bits and `2.0` for an edge bit**, and all eight per-corner influences
//!   are equal on every bit of `triangle_counts`.
//! - `P-167`'s C3 was falsified on every row: the branchless spectral evaluation
//!   loses to the lookup by `886×` on the dense bits.
//!
//! What `P-167` did **not** do is ask whether the spectrum describes any
//! perturbation this crate suffers. That is the whole of this row, and the
//! registration is explicit that the expected answer is no.
//!
//! # C1: the flip rate, derived rather than quoted
//!
//! Write the bit in the Fourier convention, `chi = 1 - 2f: {0,1}^8 -> {-1,+1}`.
//! Let `y` be the `rho`-correlated copy of a uniform `x`: independently for each
//! coordinate, `y_i = x_i` with probability `rho` and uniform otherwise, so
//!
//! ```text
//! Pr[y_i = x_i] = rho + (1 - rho)/2 = (1 + rho)/2,     eps := Pr[y_i != x_i] = (1 - rho)/2
//! ```
//!
//! and therefore `rho = 1 - 2*eps`. Noise stability is `Stab_rho = E[chi(x)chi(y)]`
//! and the standard spectral formula is `Stab_rho = sum_S rho^|S| fhat(S)^2`,
//! which is what [`Bool8::noise_stability`](common::boolean::Bool8::noise_stability)
//! computes. Because `chi` is `+/-1`-valued the product `chi(x)chi(y)` is itself
//! `+/-1`, so
//!
//! ```text
//! Stab_rho = Pr[same] - Pr[differ] = 1 - 2*Pr[differ]
//!   =>   predicted_flip_rate = Pr[chi(y) != chi(x)] = (1 - Stab_rho) / 2
//! ```
//!
//! That is the whole of C1 and there is no free parameter in it.
//!
//! **The dense spectrum makes the curve fall off fast.** `Stab_rho` weights
//! degree `k` by `rho^k`, so a function whose mass sits high in degree loses
//! stability quickly as `rho` drops. Computed here from the committed
//! `fourier_weight_by_degree` of `p-167.csv`, `predicted_flip_rate` for
//! `triangle_counts.bit1` is `0.019617` at `rho = 0.99`, `0.090882` at `0.95`,
//! `0.165627` at `0.90`, `0.277348` at `0.80` and `0.436706` at `0.50` — one
//! corner sign in two hundred going wrong already costs a two-percent output
//! flip rate. For the parity bit, where the closed form is exactly `rho^8`, the
//! same five points are `0.038628`, `0.168290`, `0.284766`, `0.416114`,
//! `0.498047`.
//!
//! **The small-`eps` linearisation names the constant.** Expanding
//! `1 - (1 - 2eps)^k = 2k*eps + O(eps^2)`,
//!
//! ```text
//! predicted_flip_rate = (1/2) * sum_k w_k * (1 - rho^k)  ~  eps * sum_k k*w_k  =  eps * I(f)
//! ```
//!
//! where `I(f) = sum_S |S| fhat(S)^2` is the total influence. So the independent
//! model's prediction, to first order, is **`I(f)` output flips per unit
//! per-corner flip rate** — `4.0` for `bit1`, `3.5` for `bit2`, `8.0` for parity,
//! `0` for the constant. `predicted_flip_rate_linear` records `eps * I(f)` beside
//! the exact figure, and `measured_response_per_corner_flip` records the same
//! ratio taken from the measurement, so the disagreement can be read as a single
//! dimensionless number: `response_ratio`.
//!
//! # C2: the wrong noise model, and why "round to `f32`" is not it
//!
//! `f32` rounding is not eight independent coin flips, and the registration says
//! so before any number exists. Two mechanisms are in play and only the second
//! can move a sign at all:
//!
//! 1. **Rounding the *value* cannot flip its sign.** `v as f32` has relative
//!    error at most `2^-24`, so `v` and `v as f32` have the same sign for every
//!    `v` whose magnitude is above the `f32` subnormal floor. This harness
//!    measures that rather than asserting it: `f32_roundtrip_sign_flips` counts
//!    the grid points where `is_inside(v)` differs from
//!    `is_inside(v as f32 as f64)`, and it is expected to be **zero on every
//!    field**. A sweep over ULP-multiples *of the result* would therefore be a
//!    sweep over eight orders of magnitude of exactly zero, which is precisely
//!    the vacuous comparison the registration's vacuity control forbids.
//! 2. **Rounding the *arithmetic* can.** The error in evaluating a signed
//!    distance is one ULP of the **operands**, not of the result. `Sphere`'s
//!    `sample` is `sqrt(x^2+y^2+z^2) - 1`: near the surface the subtraction
//!    cancels to `~0` while the absolute error stays at `~ULP(1) = 6e-8`. So the
//!    sign of a corner whose true magnitude is below that is not determined by
//!    the field, it is determined by the rounding — and that error field is
//!    **deterministic** (the same point always gives the same answer) and
//!    **spatially correlated** (it can only matter in a thin shell about the zero
//!    set, which is where every cell that produces a triangle lives).
//!
//! So the perturbation swept here is mechanism 2, amplified:
//!
//! ```text
//! e(p)          = f32_field.sample(p as f32) - f64_field.sample(p)     // the crate's own error
//! v_pert(p, k)  = f64_field.sample(p) + k * e(p)                       // ulp_perturbation = k
//! ```
//!
//! `k = 1` is the genuine artefact: the value a caller running the crate at `f32`
//! actually gets. `k > 1` scales the amplitude while keeping the *shape* of the
//! real error field, which is the one-parameter family the registration asks for
//! — "sweep the perturbation magnitude (ULP multiples) so the comparison is a
//! curve rather than a point". `k` is a multiple of the crate's own single
//! precision error, and since that error is `~ULP` of the operands it is also a
//! multiple of an ULP, which is why the column keeps its registered name.
//!
//! Three honest notes about the `f32` instance, none of which is a defect:
//!
//! - The grid points are the same real points at both precisions. Every domain
//!   half-extent is `2`, `7` or `8` and every cell size is a power-of-two
//!   fraction of it, so `p as f32` is exact; `grid_point_cast_exact` measures
//!   that over all eight fields rather than assuming it. The measured `e` is
//!   therefore evaluation precision and not a different sample location.
//! - `noise_cavity`, `fbm_terrain` and `torus` carry decimal constants
//!   (`frequency = 3.45`, `iso = 0.25`, `minor = 0.3`) that round differently in
//!   the two precisions, so their `f32` instance is a very slightly different
//!   field as well as a less precise evaluator. That is exactly what a caller who
//!   instantiates `FbmTerrain::<f32>` receives, so it belongs in the measurement.
//!   The lattice hash is pure integer arithmetic on `floor(p)` (`fields/noise.rs`),
//!   so both precisions select the same lattice cell and the same gradients.
//! - **A field can be immune to the whole sweep, and that is a measurement.**
//!   The perturbation is `k * e`, so a point where `e` is *exactly* zero cannot
//!   change sign at any `k` whatsoever. An exact polyhedral distance sampled on a
//!   dyadic grid is precisely that case near its own faces: the value there is a
//!   difference of exactly representable numbers and both precisions return it
//!   bit-identically, so the error lives only where the field is far from zero
//!   and no amount of amplification reaches the sign. `err_nonzero_points` and
//!   `err_min_abs_value_where_nonzero` are how a reader tells that structural
//!   zero from a sweep that merely stopped too early: the second is the smallest
//!   `|v|` at which the two precisions disagree at all, and a sign can only move
//!   where `k*|e|` exceeds `|v|`.
//!
//! # The comparison is calibrated, which is what makes it a test
//!
//! Comparing a measured flip rate against the model at some arbitrary `rho` would
//! measure the choice of `rho`. Instead every row **measures its own `eps`**: the
//! per-corner sign-flip rate over the cells' corner slots, `noise_rate`. Then
//! `rho = 1 - 2*noise_rate` and the model is evaluated *there*. The two sides now
//! differ in exactly one respect — the model's eight flips are independent and
//! uniform over the 256 cases, the measurement's are neither — so the gap
//! measures the modelling error and nothing else.
//!
//! Because `noise_rate` is measured, the `rho` sweep and the ULP sweep are the
//! same sweep: `eps` runs from `0` at `k = 1` up towards the saturating end at
//! `k = 10^7`, giving one `rho` per (field, `k`) pair — forty-eight of them,
//! every one a value the crate actually produced rather than a value chosen for
//! the plot. The a-priori model curve is recorded as well, at the fixed
//! [`RHO_GRID`], in `rho_grid`, `stability_at_rho_grid` and
//! `predicted_at_rho_grid` (`|`-joined, because the CSV writer refuses a comma),
//! so C1's "the predicted case-flip rate at a given `rho` is stated" is stated on
//! every row.
//!
//! **The disagreement is decomposed, not just reported.** Two of its causes are
//! separable and both are recorded:
//!
//! - *Distribution.* The model draws `x` uniformly from 256 cases; a real grid is
//!   overwhelmingly case `0` and case `255` (`empty_case_share`).
//!   `predicted_flip_rate_case_weighted` re-runs the **same independent-flip
//!   model** with `x` drawn from the field's measured case histogram instead, so
//!   `case_weighted_agreement` is the gap that survives after the distributional
//!   mismatch is removed.
//! - *Amplitude.* `eps_needed_for_measured` inverts the model curve: the `eps` at
//!   which independent noise would have predicted the rate actually measured.
//!   `eps_needed_ratio` is that over `noise_rate`, so a value far from `1` says
//!   how many times more independent noise the model needs to reproduce what
//!   correlated rounding did.
//!
//! # Arms
//!
//! One row per (output bit, reference field, `k`) — 5 × 8 × 6 = 240 rows.
//!
//! | arm | what it varies | `is_control` |
//! |---|---|---|
//! | `triangle_counts.bit1`, `.bit2` | the two informative bits of the primary reading — **the rows to quote** | no |
//! | `triangle_counts.bit0` | parity: `Stab_rho = rho^8` in closed form, asserted | yes |
//! | `triangle_counts.bit3` | the constant zero: `Stab_rho = 1`, so `predicted_flip_rate = 0` at every `rho` | yes |
//! | `edge_masks.bit0` | a single degree-2 coefficient: `Stab_rho = rho^2`, a closed form at a degree neither of the other two occupies | yes |
//! | `ulp_perturbation = 1` | the crate's actual `f32` field against its `f64` field | no |
//! | `ulp_perturbation = 10^2 … 10^7` | the same error field, amplified | no |
//!
//! The three control bits are not a separate fixture: they are bits of the same
//! two shipped tables, and each one's stability has a closed form that this file
//! computes independently and asserts. Rows whose model prediction cannot move
//! (`role = constant`) or whose measured `eps` is zero carry `is_degenerate` and
//! are excluded from the run-level C2 verdict, since a zero compared against a
//! zero is the failure M-44 names.
//!
//! # SHARE, recomputed before the numbers
//!
//! Registered: *"SHARE: none — this is a modelling check."* Discharged as
//! registered, and the arithmetic says why it must be none. The two candidate
//! consumers of a robustness number would be a refinement heuristic and a
//! precision choice, and both need the number to be *about* the perturbation they
//! face. `predicted_flip_rate ~ eps * I(f)` with `I = 4`, and the `eps` the `f32`
//! path actually delivers is bounded by the fraction of corners within `~6e-8` of
//! the zero set — on a 65³ grid over a domain of span 4 that is on the order of
//! `1e-7` of the corner slots. A stage-moving recommendation cannot be built on a
//! prediction whose input is seven orders of magnitude away from the measurement,
//! and saying so quantitatively is the deliverable. No stage moves.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record` and every panic starts `VOID: `.
//!
//! - **The measured flip rate must be non-zero at some perturbation magnitude.**
//!   The registration's own control, verbatim. `vacuity_max_measured_flip_rate`
//!   is the maximum of `measured_flip_rate_f32` over the whole sweep.
//! - **The measured per-corner flip rate must be non-zero somewhere**, or `rho`
//!   is identically `1`, the prediction is identically zero, and C2 compares two
//!   zeros for a second reason.
//! - **At least one non-degenerate row must exist** — a row with
//!   `noise_rate > 0` **and** `predicted_flip_rate > 0`. Without it the C2
//!   verdict is taken over an empty set.
//! - **The `f32` error field must be non-zero on every field.** If `e` were
//!   identically zero the whole one-parameter family would be the identity map
//!   and `k` would sweep nothing. `perturbation_max_abs` at `k = 1`.
//! - **The cells must not be uniformly empty.** At least one field must have
//!   `empty_case_share < 1`, or the output bit is constant over the population
//!   for reasons that have nothing to do with the clause.
//! - **The grid points must be the same points at both precisions.**
//!   `grid_point_cast_exact`, so `e` is evaluation error and not a relocated
//!   sample.
//! - **The instrument, against closed forms computed in this file.** For every
//!   bit whose spectrum is a single coefficient, `Stab_rho` must equal `rho^deg`
//!   to [`NEGLIGIBLE`] across [`RHO_GRID`] — `rho^8` for parity, `rho^0 = 1` for
//!   the constant, `rho^2` for an edge bit. `closed_form`, `closed_form_max_gap`.
//! - **Stability must be computable two ways.** C1's falsifier is *"stability not
//!   being computable from the spectrum"*, so the spectral sum is checked against
//!   an exact combinatorial recomputation that contains no transform at all:
//!   `Stab = sum_x w_x sum_y eps^d (1-eps)^(8-d) chi(x) chi(y)` with
//!   `d = popcount(x xor y)`, all 65,536 pairs. Asserted across [`RHO_GRID`] for
//!   every bit and recorded per row as `stability_gap_vs_combinatorial`.
//! - **`common::boolean::self_check()`** against `1`, `8*(70/256)^2` and
//!   `8*35/128` computed here rather than transcribed, matching `P-167`'s control.
//! - **The informative bits must be non-constant**, or the model predicts zero
//!   flips whatever the noise is.
//!
//! # What this row does not claim
//!
//! It does not claim that the amplified error field is a *model* of anything. It
//! is the crate's own error scaled, which is the only one-parameter family whose
//! `k = 1` member is the artefact under discussion. It does not claim the eight
//! fields are a sample of anything wider than themselves. And it does not touch
//! `crates/isomesh/src/**`: `common::boolean` is `R-167`'s module, consumed
//! unchanged.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use isomesh::marching_cubes::table::is_inside;
use isomesh::{Sdf, for_each_reference_field};

use crate::common::boolean::{
    Bool8, NEGLIGIBLE, self_check, shipped_edge_masks, shipped_triangle_counts,
};

// ─── clause constants ───────────────────────────────────────────────────────

/// Samples per axis, so `SAMPLES - 1` cells per axis.
///
/// 65, giving 274,625 samples and 262,144 cells per field. The resolution is a
/// *detector* here rather than a fidelity choice: the measured `eps` is the
/// fraction of corner slots the perturbation flips, and a flip can only happen
/// within `k * |e|` of the zero set, so the smallest `eps` this harness can
/// resolve is one slot in `8 * 262,144 = 2,097,152`. 65³ is also the largest of
/// the three golden resolutions, so the case histogram is the one the rest of the
/// repository's numbers are taken over.
const SAMPLES: u32 = 65;

/// The amplification of the crate's own single-precision error, swept.
///
/// `k = 1` is the artefact: the `f32` field against the `f64` field, no scaling.
/// The rest are decades. Two decades of gap at the bottom because nothing is
/// expected to move there — `|e| ~ 6e-8` and a cell is `0.0625` across, so a
/// perturbation of `6e-6` reaches a shell `1e-4` cells deep — and single decades
/// at the top, where the curve does its work. `10^7` puts `k*|e|` at order `1` on
/// a unit-scale field, which saturates and is recorded as `model_saturated`.
const K_SWEEP: [u64; 6] = [1, 100, 10_000, 100_000, 1_000_000, 10_000_000];

/// The `rho` values at which the a-priori model curve is stated on every row.
///
/// Descending, and stopping at `0` (`eps = 0.5`, every corner independently
/// re-randomised) because that is where the model saturates: `Stab_0 =
/// fhat(empty)^2` and the flip rate cannot exceed `(1 - fhat(empty)^2)/2`.
const RHO_GRID: [f64; 7] = [0.99, 0.95, 0.90, 0.80, 0.50, 0.20, 0.00];

/// The relative gap above which the two flip rates count as disagreeing.
///
/// A quarter. The comparison is between two probabilities that may both be very
/// small, so an absolute bar would call `1e-9` against `1e-6` an agreement; the
/// relative form asks the question the clause asks, which is whether the model
/// gets the *size* of the effect right. A quarter is generous to the model on
/// purpose: C2 predicts disagreement, so the bar is set where the model would
/// have to be badly wrong to clear it, not where it would have to be perfect.
const AGREEMENT_RELATIVE: f64 = 0.25;

/// Floor on the denominator of a relative gap, so two exact zeros give `0`
/// rather than `NaN`. Such a row is `is_degenerate` and excluded from the
/// verdict; the floor only keeps the column numeric.
const AGREEMENT_FLOOR: f64 = 1e-12;

/// Bisection steps used to invert the model curve for `eps_needed_for_measured`.
///
/// Sixty halvings of `[0, 0.5]` reach `4e-19`, well below the resolution of the
/// rate being inverted. The curve is monotone in `eps` on this interval: `rho =
/// 1 - 2*eps` decreases, every `rho^k` with `k >= 1` decreases with it, so
/// `Stab` decreases and the flip rate increases.
const BISECTION_STEPS: usize = 60;

/// The uniform distribution on the 256 case indices — the measure the spectral
/// formula is taken over, used to check the combinatorial recomputation against
/// it.
const UNIFORM: [f64; 256] = [1.0 / 256.0; 256];

// ─── exact recomputation of the same quantity, with no transform in it ──────

/// Noise stability under independent per-coordinate flips at rate `eps`, summed
/// directly over all 65,536 input pairs.
///
/// `Stab = sum_x w_x sum_y Pr[y|x] chi(x) chi(y)` with
/// `Pr[y|x] = eps^d (1-eps)^(8-d)` and `d = popcount(x xor y)`. Nothing here
/// knows what a Fourier coefficient is, which is the point: C1's falsifier is
/// *"stability not being computable from the spectrum"*, and the only way to
/// answer it is to compute the same number a second way.
///
/// `weight` is a probability distribution over the 256 case indices. With
/// [`UNIFORM`] this is exactly `sum_S rho^|S| fhat(S)^2`; with a field's measured
/// case histogram it is the same independent-flip model asked about the
/// population that actually occurs, which is how the distributional part of C2's
/// disagreement is separated from the correlational part.
fn stability_weighted(chi: &[f64; 256], eps: f64, weight: &[f64; 256]) -> f64 {
    let keep = 1.0 - eps;
    let mut by_distance = [0.0f64; 9];
    for (d, slot) in by_distance.iter_mut().enumerate() {
        let mut w = 1.0f64;
        for _ in 0..d {
            w *= eps;
        }
        for _ in d..8 {
            w *= keep;
        }
        *slot = w;
    }

    let mut acc = 0.0f64;
    for (x, (&wx, &cx)) in weight.iter().zip(chi.iter()).enumerate() {
        if wx <= 0.0 {
            continue;
        }
        let mut inner = 0.0f64;
        for (y, &cy) in chi.iter().enumerate() {
            inner += by_distance[(x ^ y).count_ones() as usize] * cy;
        }
        acc += wx * cx * inner;
    }
    acc
}

/// The flip rate the independent model predicts at `eps`, from the spectrum.
fn flip_rate_at(bit: &Bool8, eps: f64) -> f64 {
    0.5 * (1.0 - bit.noise_stability(1.0 - 2.0 * eps))
}

/// The `eps` at which independent noise would predict `target`, and whether the
/// model saturated before reaching it.
///
/// Saturation is a real outcome and not an error: the model's flip rate is
/// bounded above by `(1 - fhat(empty)^2)/2`, so a measurement above that ceiling
/// is one the independent-flip model cannot produce at any noise rate. Reported
/// as `eps = 0.5` with `model_saturated = true` rather than extrapolated.
fn eps_for_flip_rate(bit: &Bool8, target: f64) -> (f64, bool) {
    if target <= 0.0 {
        return (0.0, false);
    }
    if target >= flip_rate_at(bit, 0.5) {
        return (0.5, true);
    }
    let mut lo = 0.0f64;
    let mut hi = 0.5f64;
    for _ in 0..BISECTION_STEPS {
        let mid = 0.5 * (lo + hi);
        if flip_rate_at(bit, mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (0.5 * (lo + hi), false)
}

/// `|a - b|` over the larger of the two, floored so two zeros give zero.
fn relative_gap(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs()).max(AGREEMENT_FLOOR)
}

/// Nine or seven values as `a|b|c`, because `Run::record` refuses a comma.
fn joined(values: &[f64]) -> String {
    let parts: Vec<String> = values.iter().map(|v| format!("{v:.9}")).collect();
    parts.join("|")
}

// ─── one output bit, everything the spectrum says about it ──────────────────

/// One bit of one shipped table, transformed once and then reused across all 48
/// (field, `k`) measurements.
struct BitAnalysis {
    /// Which shipped table the bit came out of.
    reading: &'static str,
    /// `reading.bitN`, the row's identity.
    label: String,
    /// The `0/1` truth table, for counting output flips against a case index.
    table: [u8; 256],
    /// The `+/-1` truth table, for the combinatorial recomputation.
    chi: [f64; 256],
    /// The bit itself, for [`Bool8::noise_stability`].
    bit: Bool8,
    /// Bit index within the reading.
    index: u32,
    /// Non-negligible Fourier coefficients, out of 256. `P-167` measured 256 for
    /// the two informative bits and 1 for the three controls.
    spectral_terms: usize,
    /// Fourier degree.
    max_degree: u32,
    /// `sum_S |S| fhat(S)^2`, which is the small-`eps` slope of the model curve.
    total_influence: f64,
    /// Derived from the measurement, never transcribed.
    role: &'static str,
    /// `rho^deg` when the spectrum is a single coefficient, else `none`.
    closed_form: String,
    /// Worst deviation from that closed form across [`RHO_GRID`]; `0` when there
    /// is none, and the column is only meaningful when `closed_form != none`.
    closed_form_max_gap: f64,
    /// Worst spectral-versus-combinatorial gap across [`RHO_GRID`].
    grid_gap_vs_combinatorial: f64,
    /// The a-priori model curve, stated on every row.
    stability_at_grid: Vec<f64>,
    /// `(1 - Stab)/2` at the same points.
    predicted_at_grid: Vec<f64>,
}

/// Transform one bit and derive everything that does not depend on a field.
///
/// The closed form is *derived*, not looked up: by Parseval a spectrum with a
/// single non-negligible coefficient has that coefficient squared equal to `1`,
/// so `Stab_rho = rho^deg` exactly. Parity (`deg = 8`), the constant (`deg = 0`)
/// and a two-corner edge parity (`deg = 2`) all fall out of the same line, which
/// is why three different calibrations cost one branch.
fn analyse(values: &[u32; 256], reading: &'static str, index: u32) -> BitAnalysis {
    let bit = Bool8::from_values(values, index);
    let spectrum = bit.fourier();
    let table = bit.0;
    let mut chi = [0.0f64; 256];
    for (slot, &t) in chi.iter_mut().zip(table.iter()) {
        *slot = 1.0 - 2.0 * f64::from(t);
    }

    let spectral_terms = spectrum.iter().filter(|c| c.abs() > NEGLIGIBLE).count();
    let max_degree = bit.max_degree();
    let single = spectral_terms == 1;
    let role = if max_degree == 0 {
        "constant"
    } else if single && max_degree == 8 {
        "parity"
    } else if single {
        "single_coefficient"
    } else {
        "informative"
    };

    let mut stability_at_grid = Vec::with_capacity(RHO_GRID.len());
    let mut predicted_at_grid = Vec::with_capacity(RHO_GRID.len());
    let mut closed_form_max_gap = 0.0f64;
    let mut grid_gap_vs_combinatorial = 0.0f64;
    for &rho in &RHO_GRID {
        let stability = bit.noise_stability(rho);
        let combinatorial = stability_weighted(&chi, 0.5 * (1.0 - rho), &UNIFORM);
        grid_gap_vs_combinatorial =
            grid_gap_vs_combinatorial.max((stability - combinatorial).abs());
        if single {
            let mut closed = 1.0f64;
            for _ in 0..max_degree {
                closed *= rho;
            }
            closed_form_max_gap = closed_form_max_gap.max((stability - closed).abs());
        }
        stability_at_grid.push(stability);
        predicted_at_grid.push(0.5 * (1.0 - stability));
    }

    BitAnalysis {
        reading,
        label: format!("{reading}.bit{index}"),
        table,
        chi,
        total_influence: bit.total_influence(),
        bit,
        index,
        spectral_terms,
        max_degree,
        role,
        closed_form: if single {
            format!("rho^{max_degree}")
        } else {
            String::from("none")
        },
        closed_form_max_gap,
        grid_gap_vs_combinatorial,
        stability_at_grid,
        predicted_at_grid,
    }
}

// ─── the fields, sampled once at both precisions ────────────────────────────

/// One reference field, sampled at `f64` and at `f32`, with the case indices the
/// shipped classifier assigns to the unperturbed grid.
struct FieldSamples {
    /// The `ReferenceField` name.
    name: &'static str,
    /// The `f64` value at every grid point, `x` fastest.
    v64: Vec<f64>,
    /// `f32_field.sample(p as f32) - f64_field.sample(p)`: the crate's own
    /// single-precision error, absolute rather than relative, at every point.
    err: Vec<f64>,
    /// Grid points whose sign changes under a plain `f64 -> f32` round trip of
    /// the value. Measured rather than assumed; expected zero everywhere.
    roundtrip_sign_flips: u64,
    /// Grid points where the two precisions disagree at all.
    ///
    /// The perturbation is `k * e`, so a point with `e == 0` can never change
    /// sign at any `k`. This column and the next are how a reader tells a field
    /// whose error is structurally absent near its own zero set from a sweep that
    /// is simply too small.
    err_nonzero_points: u64,
    /// The smallest `|v64|` among the points where `e != 0`.
    ///
    /// Threshold-free, and it decides the whole question for a field: a sign can
    /// only move where `k*|e| > |v64|`, so if the field is evaluated *exactly*
    /// everywhere near its own zero set this number is large and no `k` will ever
    /// flip a corner. An exact polyhedral distance sampled on a dyadic grid is
    /// the case that does this — its value near a face is a difference of exactly
    /// representable numbers and both precisions get it bit-identically.
    err_min_abs_value_where_nonzero: f64,
    /// Unperturbed case index of every cell.
    case64: Vec<u8>,
    /// The unperturbed case histogram, normalised to a probability distribution.
    weight: [f64; 256],
    /// Share of cells that are case `0` or case `255`.
    empty_case_share: f64,
}

/// Sample a field at `f64` over the grid `common::grid` defines.
fn sample_at<F: Sdf<Scalar = f64>>(field: &F, origin: [f64; 3], h: f64) -> Vec<f64> {
    let n = SAMPLES as usize;
    let mut out = Vec::with_capacity(n * n * n);
    for z in 0..n {
        let pz = origin[2] + h * z as f64;
        for y in 0..n {
            let py = origin[1] + h * y as f64;
            for x in 0..n {
                out.push(field.sample([origin[0] + h * x as f64, py, pz]));
            }
        }
    }
    out
}

/// Sample the `f32` instance of a field at the `f32` cast of the same points.
///
/// Returns the values widened back to `f64` — widening is exact, so no
/// information is added or lost — and whether every coordinate cast was exact,
/// which is what makes the difference against [`sample_at`] evaluation precision
/// rather than a relocated sample. Compared through `to_bits` because an exact
/// float equality is the one comparison a float equality is right for.
fn sample_at_f32<F: Sdf<Scalar = f32>>(field: &F, origin: [f64; 3], h: f64) -> (Vec<f64>, bool) {
    let n = SAMPLES as usize;
    let mut out = Vec::with_capacity(n * n * n);
    let mut exact = true;
    for z in 0..n {
        let pz = origin[2] + h * z as f64;
        for y in 0..n {
            let py = origin[1] + h * y as f64;
            for x in 0..n {
                let px = origin[0] + h * x as f64;
                for c in [px, py, pz] {
                    exact &= f64::from(c as f32).to_bits() == c.to_bits();
                }
                out.push(f64::from(field.sample([px as f32, py as f32, pz as f32])));
            }
        }
    }
    (out, exact)
}

/// The case index of every cell, from the per-sample inside flags.
///
/// Bit `i` of the case is corner `i`, at local coordinate
/// `(i & 1, (i >> 1) & 1, (i >> 2) & 1)` — the numbering `corner_inside` uses and
/// the numbering `common::boolean` builds its `Bool8` variables from, so a case
/// index computed here indexes a truth table built there.
fn case_indices(inside: &[bool]) -> Vec<u8> {
    let n = SAMPLES as usize;
    let cells = n - 1;
    let mut out = Vec::with_capacity(cells * cells * cells);
    for cz in 0..cells {
        for cy in 0..cells {
            for cx in 0..cells {
                let mut case = 0u8;
                for corner in 0..8usize {
                    let x = cx + (corner & 1);
                    let y = cy + ((corner >> 1) & 1);
                    let z = cz + ((corner >> 2) & 1);
                    if inside[x + y * n + z * n * n] {
                        case |= 1u8 << corner;
                    }
                }
                out.push(case);
            }
        }
    }
    out
}

/// Sample all eight reference fields at both precisions.
///
/// Two passes, because the macro inlines its body once per field and the eight
/// fields are eight distinct types at each precision — there is no way to hold an
/// `f32` and an `f64` instance of the same field in one block. The passes are
/// joined by position and the name is asserted at every position, so a change to
/// the macro's order cannot silently pair `gyroid`'s `f64` values with
/// `fbm_terrain`'s `f32` ones.
fn sample_fields() -> (Vec<FieldSamples>, bool) {
    let mut geometry: Vec<(&'static str, [f64; 3], f64)> = Vec::new();
    let mut v64s: Vec<Vec<f64>> = Vec::new();
    for_each_reference_field!(f64, |name, field| {
        let (_shape, origin, h) = common::grid::<f64, _>(&field, SAMPLES);
        geometry.push((name, origin, h));
        v64s.push(sample_at(&field, origin, h));
    });

    let mut v32s: Vec<Vec<f64>> = Vec::new();
    let mut cast_exact = true;
    for_each_reference_field!(f32, |name, field| {
        let position = v32s.len();
        let (expected, origin, h) = geometry[position];
        assert_eq!(
            expected, name,
            "the two passes over the eight reference fields disagree at position {position}"
        );
        let (values, exact) = sample_at_f32(&field, origin, h);
        cast_exact &= exact;
        v32s.push(values);
    });

    let cells = {
        let c = SAMPLES as usize - 1;
        c * c * c
    };
    let mut out = Vec::with_capacity(v64s.len());
    for ((name, _origin, _h), (v64, v32)) in geometry.into_iter().zip(v64s.into_iter().zip(v32s)) {
        let err: Vec<f64> = v64.iter().zip(v32.iter()).map(|(&a, &b)| b - a).collect();
        let roundtrip_sign_flips = v64
            .iter()
            .filter(|&&v| is_inside(v) != is_inside(f64::from(v as f32)))
            .count() as u64;

        let mut err_nonzero_points = 0u64;
        let mut err_min_abs_value_where_nonzero = f64::INFINITY;
        for (&v, &e) in v64.iter().zip(err.iter()) {
            if e.abs() > 0.0 {
                err_nonzero_points += 1;
                err_min_abs_value_where_nonzero = err_min_abs_value_where_nonzero.min(v.abs());
            }
        }

        let inside: Vec<bool> = v64.iter().map(|&v| is_inside(v)).collect();
        let case64 = case_indices(&inside);
        let mut histogram = [0u64; 256];
        for &case in &case64 {
            histogram[usize::from(case)] += 1;
        }
        let mut weight = [0.0f64; 256];
        for (slot, &count) in weight.iter_mut().zip(histogram.iter()) {
            *slot = count as f64 / cells as f64;
        }
        let empty_case_share = (histogram[0] + histogram[255]) as f64 / cells as f64;

        out.push(FieldSamples {
            name,
            v64,
            err,
            roundtrip_sign_flips,
            err_nonzero_points,
            err_min_abs_value_where_nonzero,
            case64,
            weight,
            empty_case_share,
        });
    }
    (out, cast_exact)
}

// ─── one perturbation magnitude ─────────────────────────────────────────────

/// The field re-classified under `v + k*e`, and what moved.
struct Perturbed {
    /// The amplification.
    k: u64,
    /// Per-cell case index under the perturbation.
    case: Vec<u8>,
    /// Corner slots whose sign changed, out of `8 * cells`. This is the `eps` the
    /// model is calibrated to, and it is `sum_cells popcount(case64 xor case)`
    /// exactly, because bit `i` of a case index *is* the sign of corner `i`.
    corner_sign_flips: u64,
    /// Cells whose case index changed at all — the triangulation for that cell is
    /// a different entry of the table.
    topology_changes: u64,
    /// `max |k*e|` in field units, so a reader can see when the "perturbation"
    /// has stopped being small.
    max_abs: f64,
    /// `rms |k*e|` over all grid points.
    rms: f64,
}

/// Re-classify one field at one amplification.
fn perturb(samples: &FieldSamples, k: u64) -> Perturbed {
    let scale = k as f64;
    let mut max_abs = 0.0f64;
    let mut square_sum = 0.0f64;
    let inside: Vec<bool> = samples
        .v64
        .iter()
        .zip(samples.err.iter())
        .map(|(&v, &e)| {
            let delta = scale * e;
            max_abs = max_abs.max(delta.abs());
            square_sum += delta * delta;
            is_inside(v + delta)
        })
        .collect();

    let case = case_indices(&inside);
    let mut corner_sign_flips = 0u64;
    let mut topology_changes = 0u64;
    for (&before, &after) in samples.case64.iter().zip(case.iter()) {
        let moved = before ^ after;
        corner_sign_flips += u64::from(moved.count_ones());
        if moved != 0 {
            topology_changes += 1;
        }
    }

    Perturbed {
        k,
        case,
        corner_sign_flips,
        topology_changes,
        max_abs,
        rms: (square_sum / samples.v64.len() as f64).sqrt(),
    }
}

// ─── one measured row ───────────────────────────────────────────────────────

/// One (bit, field, `k`) measurement: the model at the measured `eps`, the
/// measurement, and the gap.
struct Row {
    /// Index into the analysed bits.
    bit: usize,
    /// Reference field name.
    field: &'static str,
    /// The amplification.
    k: u64,
    /// Measured per-corner sign-flip rate — the model's `eps`.
    noise_rate: f64,
    /// `1 - 2*noise_rate`.
    rho: f64,
    /// `Stab_rho` from the spectrum.
    stability: f64,
    /// `Stab_rho` recomputed over all 65,536 pairs, no transform involved.
    stability_combinatorial: f64,
    /// `(1 - Stab_rho)/2`.
    predicted: f64,
    /// The same independent-flip model, over the field's own case histogram.
    predicted_case_weighted: f64,
    /// `noise_rate * total_influence`, the small-`eps` linearisation.
    predicted_linear: f64,
    /// Cells whose output bit changed, over cells.
    measured: f64,
    /// The same as a count.
    output_flips: u64,
    /// Cells whose case index changed at all.
    topology_changes: u64,
    /// Corner slots whose sign changed.
    corner_sign_flips: u64,
    /// `measured / noise_rate` — the empirical total influence.
    measured_response: f64,
    /// That over `total_influence`.
    response_ratio: f64,
    /// The `eps` at which the model would have predicted `measured`.
    eps_needed: f64,
    /// Whether the model saturated before reaching `measured`.
    model_saturated: bool,
    /// `eps_needed / noise_rate`.
    eps_needed_ratio: f64,
    /// `max |k*e|`.
    max_abs: f64,
    /// `rms |k*e|`.
    rms: f64,
    /// Share of cells that are case 0 or 255.
    empty_case_share: f64,
    /// The naive reading's flip count, per field.
    roundtrip_sign_flips: u64,
    /// Grid points where the two precisions disagree at all, per field.
    err_nonzero_points: u64,
    /// Smallest `|v64|` among those points, per field.
    err_min_abs_value_where_nonzero: f64,
    /// Relative gap between `measured` and `predicted`.
    agreement: f64,
    /// Absolute gap between the same two.
    agreement_absolute: f64,
    /// Relative gap after the distributional mismatch is removed.
    case_weighted_agreement: f64,
    /// Spectral against combinatorial, at this row's own `rho`.
    gap_vs_combinatorial: f64,
    /// The model cannot move on this row, so C2 would compare two zeros.
    degenerate: bool,
    /// Stability was computed from the spectrum and verified two other ways.
    c1: bool,
    /// The two flip rates disagree, as registered.
    c2: bool,
}

/// Cells per field, `(SAMPLES - 1)^3`.
fn cell_count() -> usize {
    let c = SAMPLES as usize - 1;
    c * c * c
}

/// Measure one bit against one perturbed field.
fn measure(analysis: &BitAnalysis, bit: usize, samples: &FieldSamples, pert: &Perturbed) -> Row {
    let cells = cell_count();
    let slots = (cells * 8) as f64;
    let noise_rate = pert.corner_sign_flips as f64 / slots;
    let rho = 1.0 - 2.0 * noise_rate;

    let stability = analysis.bit.noise_stability(rho);
    let stability_combinatorial = stability_weighted(&analysis.chi, noise_rate, &UNIFORM);
    let predicted = 0.5 * (1.0 - stability);
    let predicted_case_weighted =
        0.5 * (1.0 - stability_weighted(&analysis.chi, noise_rate, &samples.weight));
    let predicted_linear = noise_rate * analysis.total_influence;

    let output_flips = samples
        .case64
        .iter()
        .zip(pert.case.iter())
        .filter(|&(&before, &after)| {
            analysis.table[usize::from(before)] != analysis.table[usize::from(after)]
        })
        .count() as u64;
    let measured = output_flips as f64 / cells as f64;

    let (eps_needed, model_saturated) = eps_for_flip_rate(&analysis.bit, measured);
    let gap_vs_combinatorial = (stability - stability_combinatorial).abs();
    let agreement = relative_gap(measured, predicted);
    let degenerate = noise_rate <= 0.0 || predicted <= 0.0;

    let c1 = gap_vs_combinatorial <= NEGLIGIBLE
        && analysis.grid_gap_vs_combinatorial <= NEGLIGIBLE
        && analysis.closed_form_max_gap <= NEGLIGIBLE;

    Row {
        bit,
        field: samples.name,
        k: pert.k,
        noise_rate,
        rho,
        stability,
        stability_combinatorial,
        predicted,
        predicted_case_weighted,
        predicted_linear,
        measured,
        output_flips,
        topology_changes: pert.topology_changes,
        corner_sign_flips: pert.corner_sign_flips,
        measured_response: if noise_rate > 0.0 {
            measured / noise_rate
        } else {
            0.0
        },
        response_ratio: if noise_rate > 0.0 && analysis.total_influence > 0.0 {
            measured / noise_rate / analysis.total_influence
        } else {
            0.0
        },
        eps_needed,
        model_saturated,
        eps_needed_ratio: if noise_rate > 0.0 {
            eps_needed / noise_rate
        } else {
            0.0
        },
        max_abs: pert.max_abs,
        rms: pert.rms,
        empty_case_share: samples.empty_case_share,
        roundtrip_sign_flips: samples.roundtrip_sign_flips,
        err_nonzero_points: samples.err_nonzero_points,
        err_min_abs_value_where_nonzero: samples.err_min_abs_value_where_nonzero,
        agreement,
        agreement_absolute: (measured - predicted).abs(),
        case_weighted_agreement: relative_gap(measured, predicted_case_weighted),
        gap_vs_combinatorial,
        degenerate,
        c1,
        c2: !degenerate && agreement > AGREEMENT_RELATIVE,
    }
}

// ─── the run ────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-168");

    common::experiment::run(prereg, |run| {
        // ── the instrument, before anything is measured with it ─────────────
        let (parity_weight8, majority_weight1, majority_influence) = self_check();
        let majority_weight1_closed = 8.0 * (70.0 / 256.0) * (70.0 / 256.0);
        let majority_influence_closed = 8.0 * 35.0 / 128.0;
        assert!(
            (parity_weight8 - 1.0).abs() <= NEGLIGIBLE
                && (majority_weight1 - majority_weight1_closed).abs() <= NEGLIGIBLE
                && (majority_influence - majority_influence_closed).abs() <= NEGLIGIBLE,
            "VOID: common::boolean's calibration does not reproduce the closed forms computed \
             here — parity weight at degree 8 {parity_weight8} against 1, majority weight at \
             degree 1 {majority_weight1} against {majority_weight1_closed}, majority total \
             influence {majority_influence} against {majority_influence_closed}. The transform \
             is unvalidated and every stability below is a number from a broken instrument"
        );

        let triangles = shipped_triangle_counts();
        let edges = shipped_edge_masks();
        let bits = vec![
            analyse(&triangles, "triangle_counts", 0),
            analyse(&triangles, "triangle_counts", 1),
            analyse(&triangles, "triangle_counts", 2),
            analyse(&triangles, "triangle_counts", 3),
            analyse(&edges, "edge_masks", 0),
        ];

        for analysis in &bits {
            assert!(
                analysis.grid_gap_vs_combinatorial <= NEGLIGIBLE,
                "VOID: {}'s spectral noise stability disagrees with the exact 65,536-pair \
                 recomputation by {} across the rho grid, so stability is not computable from \
                 the spectrum and C1 is answered before C2 is asked",
                analysis.label,
                analysis.grid_gap_vs_combinatorial
            );
            assert!(
                analysis.closed_form_max_gap <= NEGLIGIBLE,
                "VOID: {}'s stability deviates from its closed form {} by {} across the rho \
                 grid; a single-coefficient spectrum has Stab_rho = rho^deg exactly by Parseval, \
                 so the transform or the stability sum is wrong",
                analysis.label,
                analysis.closed_form,
                analysis.closed_form_max_gap
            );
        }

        let informative = bits.iter().filter(|b| b.max_degree > 0).count();
        assert!(
            informative >= 2,
            "VOID: only {informative} of the analysed bits is non-constant, so the model \
             predicts zero flips whatever the noise is and C2 compares two zeros by \
             construction"
        );

        // ── the fields, and the perturbation family built on their own error ─
        let (fields, cast_exact) = sample_fields();
        assert!(
            cast_exact,
            "VOID: a grid coordinate does not cast exactly to f32, so the measured difference \
             between the two precisions mixes evaluation error with a relocated sample point \
             and does not isolate rounding"
        );
        for samples in &fields {
            let worst = samples.err.iter().fold(0.0f64, |acc, e| acc.max(e.abs()));
            assert!(
                worst > 0.0,
                "VOID: the f32 error field is identically zero on {}, so v + k*e is the \
                 identity for every k and the whole perturbation sweep measures nothing",
                samples.name
            );
        }
        assert!(
            fields.iter().any(|f| f.empty_case_share < 1.0),
            "VOID: every cell of every field is case 0 or case 255, so no output bit can move \
             for reasons that have nothing to do with the noise model"
        );

        // ── every row, before any of them is recorded ────────────────────────
        let mut rows: Vec<Row> = Vec::with_capacity(bits.len() * fields.len() * K_SWEEP.len());
        for samples in &fields {
            for &k in &K_SWEEP {
                let pert = perturb(samples, k);
                for (index, analysis) in bits.iter().enumerate() {
                    rows.push(measure(analysis, index, samples, &pert));
                }
            }
        }

        let max_measured = rows.iter().fold(0.0f64, |acc, r| acc.max(r.measured));
        assert!(
            max_measured > 0.0,
            "VOID: the measured flip rate is zero at every perturbation magnitude on every \
             field and every bit, so C2 compares two zeros — this is the registration's own \
             vacuity control verbatim"
        );
        let max_noise_rate = rows.iter().fold(0.0f64, |acc, r| acc.max(r.noise_rate));
        assert!(
            max_noise_rate > 0.0,
            "VOID: no corner sign flipped at any perturbation magnitude, so rho is identically \
             1, the predicted flip rate is identically 0, and the comparison is vacuous from \
             the model's side as well"
        );
        let counted = rows.iter().filter(|r| !r.degenerate).count();
        assert!(
            counted > 0,
            "VOID: every row is degenerate (measured eps of zero, or a model prediction that \
             cannot move), so the C2 verdict would be taken over an empty set"
        );
        let disagreeing = rows.iter().filter(|r| r.c2).count();

        // ── the rows ────────────────────────────────────────────────────────
        let rho_grid = joined(&RHO_GRID);
        let cells = cell_count();
        for row in &rows {
            let analysis = &bits[row.bit];
            run.record(&[
                ("noise_rate", format!("{:.12}", row.noise_rate)),
                ("noise_stability", format!("{:.12}", row.stability)),
                ("predicted_flip_rate", format!("{:.12}", row.predicted)),
                ("measured_flip_rate_f32", format!("{:.12}", row.measured)),
                ("agreement", format!("{:.9}", row.agreement)),
                ("ulp_perturbation", row.k.to_string()),
                ("topology_changes", row.topology_changes.to_string()),
                ("field", row.field.to_string()),
                ("c1_holds", row.c1.to_string()),
                ("c2_holds", row.c2.to_string()),
                // ── extras (M-273) ──
                ("output_bit", analysis.label.clone()),
                ("reading", analysis.reading.to_string()),
                ("bit_index", analysis.index.to_string()),
                ("role", analysis.role.to_string()),
                ("closed_form", analysis.closed_form.clone()),
                (
                    "closed_form_max_gap",
                    format!("{:.15}", analysis.closed_form_max_gap),
                ),
                ("spectral_terms", analysis.spectral_terms.to_string()),
                ("max_degree", analysis.max_degree.to_string()),
                (
                    "total_influence",
                    format!("{:.9}", analysis.total_influence),
                ),
                ("rho", format!("{:.12}", row.rho)),
                ("rho_grid", rho_grid.clone()),
                ("stability_at_rho_grid", joined(&analysis.stability_at_grid)),
                ("predicted_at_rho_grid", joined(&analysis.predicted_at_grid)),
                (
                    "stability_combinatorial",
                    format!("{:.12}", row.stability_combinatorial),
                ),
                (
                    "stability_gap_vs_combinatorial",
                    format!("{:.15}", row.gap_vs_combinatorial),
                ),
                (
                    "predicted_flip_rate_case_weighted",
                    format!("{:.12}", row.predicted_case_weighted),
                ),
                (
                    "predicted_flip_rate_linear",
                    format!("{:.12}", row.predicted_linear),
                ),
                (
                    "case_weighted_agreement",
                    format!("{:.9}", row.case_weighted_agreement),
                ),
                (
                    "agreement_absolute",
                    format!("{:.12}", row.agreement_absolute),
                ),
                ("agreement_relative_bar", format!("{AGREEMENT_RELATIVE:.2}")),
                (
                    "measured_response_per_corner_flip",
                    format!("{:.9}", row.measured_response),
                ),
                ("response_ratio", format!("{:.9}", row.response_ratio)),
                ("eps_needed_for_measured", format!("{:.12}", row.eps_needed)),
                ("eps_needed_ratio", format!("{:.6}", row.eps_needed_ratio)),
                ("model_saturated", row.model_saturated.to_string()),
                ("output_flips", row.output_flips.to_string()),
                ("corner_sign_flips", row.corner_sign_flips.to_string()),
                ("corner_slots", (cells * 8).to_string()),
                ("cells", cells.to_string()),
                (
                    "case_index_change_rate",
                    format!("{:.12}", row.topology_changes as f64 / cells as f64),
                ),
                ("perturbation_max_abs", format!("{:.12e}", row.max_abs)),
                ("perturbation_rms", format!("{:.12e}", row.rms)),
                (
                    "f32_roundtrip_sign_flips",
                    row.roundtrip_sign_flips.to_string(),
                ),
                ("err_nonzero_points", row.err_nonzero_points.to_string()),
                (
                    "err_min_abs_value_where_nonzero",
                    format!("{:.12e}", row.err_min_abs_value_where_nonzero),
                ),
                ("empty_case_share", format!("{:.9}", row.empty_case_share)),
                ("grid_point_cast_exact", cast_exact.to_string()),
                ("resolution", SAMPLES.to_string()),
                ("is_degenerate", row.degenerate.to_string()),
                ("c2_rows_counted", counted.to_string()),
                ("c2_rows_disagreeing", disagreeing.to_string()),
                (
                    "vacuity_max_measured_flip_rate",
                    format!("{max_measured:.12}"),
                ),
                (
                    "calibration_parity_weight_degree8",
                    format!("{parity_weight8:.12}"),
                ),
                (
                    "calibration_majority_weight_degree1",
                    format!("{majority_weight1:.12}"),
                ),
                (
                    "calibration_majority_total_influence",
                    format!("{majority_influence:.12}"),
                ),
            ]);
        }

        // ── the run-level verdict ───────────────────────────────────────────
        //
        // Over the non-degenerate rows only. A constant bit predicts zero flips
        // at every rho and a row whose measured eps is zero calibrates the model
        // to rho = 1, and in both cases the two sides are zero for structural
        // reasons rather than for the reason C2 is about.
        let c1_held = rows.iter().filter(|r| r.c1).count();
        println!();
        println!(
            "C1 held on {c1_held} of {} rows: the spectral stability agrees with the exact \
             65,536-pair recomputation, and every single-coefficient bit reproduces its closed \
             form rho^deg",
            rows.len()
        );
        println!(
            "C2: {disagreeing} of {counted} non-degenerate rows disagree by more than \
             {:.0}% relative; {} rows are degenerate (zero measured eps, or a model that \
             cannot move)",
            AGREEMENT_RELATIVE * 100.0,
            rows.len() - counted
        );

        let primary: Vec<&Row> = rows
            .iter()
            .filter(|r| bits[r.bit].role == "informative" && !r.degenerate)
            .collect();
        let worst = primary
            .iter()
            .fold(0.0f64, |acc, r| acc.max(r.eps_needed_ratio));
        let best = primary
            .iter()
            .fold(f64::INFINITY, |acc, r| acc.min(r.eps_needed_ratio));
        println!(
            "over the {} non-degenerate rows of the two informative bits, the independent model \
             needs between {best:.3}x and {worst:.3}x the measured per-corner flip rate to \
             reproduce the measured output flip rate",
            primary.len()
        );
        let response = primary
            .iter()
            .fold(0.0f64, |acc, r| acc.max(r.response_ratio));
        println!(
            "the measured output response per corner flip reaches {response:.3} of the total \
             influence the independent model predicts as its small-eps slope"
        );
        println!(
            "the naive reading — round the value itself to f32 — flipped {} signs in total \
             across all eight fields, which is why the sweep amplifies the crate's own \
             evaluation error instead of nudging result ULPs",
            fields.iter().map(|f| f.roundtrip_sign_flips).sum::<u64>()
        );
    });
}
