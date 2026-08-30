//! **P-127 — the body-saddle discriminant is Cayley's hyperdeterminant, and the
//! proof is a test rather than a comment.**
//!
//! Ticket: R-127. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p127
//! ```
//!
//! Writes `docs/experiments/p-127.csv`.
//!
//! # What was missing
//!
//! `crates/isomesh/src/marching_cubes/trilinear.rs:246` computes
//! `b*b - R::TWO * R::TWO * a * c` from the three coefficients built at
//! `:199-214`. Every id the crate has ever attached to that expression treats it
//! as **a transcription of Grosso's quadratic**:
//!
//! - `M-207` defends the *root extraction* around it — `(-b +- sqrt(d))/2a`
//!   loses the linear root when `a` is exactly zero and drops an exact tangency,
//!   so this crate solves the smaller polynomial as a smaller polynomial. That
//!   is a finding about the **solver**, and it says nothing about what the
//!   discriminant *is*.
//! - `M-206` is the odd one out and is the reason this row exists:
//!   *"two independently derived constructions locate the same body saddles, to
//!   `1.1e-12`"*. `interior::SweptFaces` solves a quadratic in the sweep
//!   **height**; `trilinear::BodySaddles` solves a quadratic in the face
//!   **coordinate**. `M-206` records that they *share no arithmetic, no
//!   coefficient and no parametrisation* and then agree to twelve digits. An
//!   unexplained correct number is exactly the artefact this row converts into a
//!   theorem: they are two of the three slicings of **one pencil**, and the
//!   invariant they both compute is the same degree-4 form.
//! - `M-214`, `M-215`, `M-216` and `M-217` all build on the saddle count without
//!   ever asking what the sign of that discriminant is invariant under.
//!
//! Nothing in the repository states the identity, and nothing gates it. The one
//! artefact that proves it —
//! `docs/research/2026-08-29-phase-27-hyperdeterminant-identity.py:45-67` —
//! needs `sympy`, which `scripts/preflight.sh:98-101` does not install into
//! `~/.venvs/isomesh`, and it is wired into no gate. So the identity has been
//! *known* for exactly as long as it has been *unenforced*, which is zero
//! coverage in the only sense that matters.
//!
//! `crates/isomesh/benches/common/poly.rs` (this ticket owns it) is the gate: a
//! degree-4 form in eight variables has at most 330 monomials, so exact
//! expansion over `i128` is a hundred lines and needs no dependency. This
//! harness is what reads the gate's answer into a file.
//!
//! # What the object is, quoted rather than paraphrased
//!
//! de Silva & Lim, *Tensor rank and the ill-posedness of the best low-rank
//! approximation problem*, `arXiv:math/0607647`, corpus `doc_id`
//! `10.48550_arXiv.math_0607647`, §6, verbatim: *"the rank of a tensor is 2 on
//! the set `{A | Det_2,2,2(A) > 0}` and 3 on the set `{A | Det_2,2,2(A) < 0}`."*
//! Their sign convention is the discriminant of `det(l1*A1 + l2*A2)`, which is
//! the normalisation `common::poly` uses — `c1^2 - 4*c0*c2`, no leading minus —
//! and therefore the crate's. That agreement is not assumed here: it is what
//! `symbolic_difference_is_zero` measures, and a sign flip would show up as a
//! twelve-term residual rather than as a zero.
//!
//! # Arms
//!
//! Three rows, one per axis pairing, because `pencil_axis_pairing` and
//! `pencil_matches` are the only registered columns that vary within the run.
//! Everything else is a **global** quantity and carries the same value on every
//! row; the header says so here so that a reader of the CSV does not mistake
//! three identical `f32_sign_disagreements` for three measurements.
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | pencil `0123\|4567` | the cube is split along `w` | no |
//! | pencil `0145\|2367` | the cube is split along `v` | no |
//! | pencil `0246\|1357` | the cube is split along `u` | no |
//! | corrupted-Cayley | the `f0^2*f7^2` coefficient is bumped `1 -> 2` | **yes**, an assert rather than a row |
//! | general rational stratum | 3,000 dyadic 8-tuples, `f_i` in `[-1024, 1024]` at spacing down to `1/8` | no |
//! | near-zero stratum | `f7` bisected onto a root of the discriminant | **yes**, the stratum that licenses the `f32` count |
//!
//! The corrupted arm is a control and not a row on purpose: its whole job is to
//! show that `Poly::sub(..).is_zero()` — the single predicate carrying both C1
//! and C2 — is *capable of reading `false`*. A control that produced a CSV row
//! would invite that row being quoted as a result.
//!
//! # The three clauses, and how each is decided
//!
//! **C1 is symbolic and exact.** `repo_discriminant()` is transcribed line for
//! line from `BodySaddles::coefficients`; `cayley_2x2x2()` is the standard
//! explicit twelve-term normalisation, written out rather than derived from the
//! pencil precisely so that C2's three pairings are three genuine checks instead
//! of one construction and two. C1 holds when all four of its own sentences do:
//! twelve terms on each side, total degree 4, and an identically zero
//! difference.
//!
//! **C2 is the same predicate three times.** `pencil_discriminant(p)` builds
//! `det(A0 + lambda*A1)` for pairing `p`, reads off `c1^2 - 4*c0*c2`, and is
//! compared to Cayley. `pencil_matches` is the **running** count over the rows,
//! so the last row reads 3 when every pairing agrees; `pencil_matches_here` is
//! that row's own boolean, which is the column to read if one pairing fails.
//!
//! **C3 is numeric, and its whole difficulty is making the `f32` count mean
//! something.** Two design decisions, both load-bearing:
//!
//! 1. **Every corner value is exactly representable in `f32`.** Each `f_i` is a
//!    dyadic rational `n / 2^k` with `|n| < 2^24`, so `num as f32 / den as f32`
//!    is exact and `f64::from(that) == num as f64 / den as f64`. That equality is
//!    *checked on every trial* (`inexact_f32_inputs`, asserted zero) rather than
//!    argued for. Without it a "the `f32` sign is wrong" count would be
//!    contaminated by the input having been rounded before any arithmetic
//!    happened, and `P-134` would inherit a number about representation instead
//!    of about cancellation.
//! 2. **The exact sign is the sign of an `i128` numerator.** `Poly::eval_ratio`
//!    clears denominators over the polynomial's own `degree_in` — which is 2 in
//!    every corner, this being a quadratic invariant and not a multi-affine one —
//!    and reduces. So "ratio exactly 1" is decided by two reduced pairs being
//!    *identical*, never by forming a quotient, and `max_abs_ratio_deviation` is
//!    reported as an exact integer ratio `p/q` rather than as a rounded decimal.
//!
//! `f32_sign_disagreements` counts trials where the `f32` evaluation of *either*
//! polynomial disagrees in sign with the exact value. Because C1 makes the two
//! polynomials the same `BTreeMap`, `Poly::eval_f32` on them is bit-identical and
//! `f32_cross_disagreements` is structurally 0 — that column is a witness that
//! C1 held all the way down to the float path, not an independent measurement.
//! A `f32` result of `+-0.0` against a non-zero exact value **counts as a
//! disagreement**, and it is the most interesting kind: `trilinear.rs:250`
//! branches on `discriminant == R::ZERO` and takes the double-root path, so a
//! flush to zero is a cell silently re-classified as a tangency.
//! `f32_flush_to_zero_trials` is that subset.
//!
//! # The near-zero stratum, and why it is bisection rather than luck
//!
//! The registration's vacuity control demands at least 50 trials within `1e-6`
//! of zero. Drawing them is hopeless — a random 8-tuple with corners of size
//! `10^3` has `|Delta|` of order `10^12` — so they are **constructed**, and the
//! construction is exact:
//!
//! Fix `f0..f6`. The discriminant is then a *quadratic in `f7` alone*:
//! `Delta(f7) = f0^2*f7^2 + (4*f1*f2*f4 - 2*f0*f3*f4 - 2*f0*f2*f5 -
//! 2*f0*f1*f6)*f7 + C`. So a sign change between two consecutive `f7` samples
//! brackets a root, and bisection converges on it. The harness scans `f7` over
//! `[-4, 4]` at spacing `1/16`, takes the first bracket, and halves — evaluating
//! the **exact rational** discriminant at every step, so the bracket is never
//! lost to rounding.
//!
//! **The depth cap is derived, not chosen.** Bisection stops when one more
//! halving would push `|numerator|` past `2^24` and take the trial off the grid
//! `f32` represents exactly — decision 1 above is worth more than another bit of
//! `f7`. `f0..f6` are drawn in `[-1, 1]` at spacing `1/16` for the same reason it
//! matters: with every corner bounded by 1 the slope `|dDelta/df7|` is at most
//! `4 + 2 + 2 + 2 + 2*|f0^2*f7|`, of order ten, so a final bracket of width
//! `2^-21` puts `|Delta|` in the `10^-7`–`10^-6` band that the control asks for
//! and a `f32` evaluation of twelve `O(1)` terms cannot resolve.
//!
//! `near_zero_trials` is the honest census — `|Delta| < 1e-6` over the **whole**
//! population, general stratum included — and `near_zero_constructed` is how
//! many the bisection produced. `near_zero_brackets_found` separates "the
//! bisection failed to converge" from "no sign change existed in `[-4, 4]` for
//! that draw", which is a real and uninteresting outcome for a quadratic whose
//! roots are complex.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration says: *"C1 and C2 move no runtime cost at all -- this row is
//! a renaming plus two theorems, and its whole value is what Group A's other rows
//! can then assume."* Discharged: **zero**, and zero by construction rather than
//! by measurement. This harness reads no field, extracts no mesh, and proposes
//! no change to `crates/isomesh/src/**`. `trilinear.rs:246` computes the same
//! twelve-term form after this row as before it.
//!
//! What stands in a share's place is what each clause licenses downstream, with
//! an exact denominator:
//!
//! | clause | quantity | denominator | exact because |
//! |---|---|---|---|
//! | C1 | `symbolic_residual_terms` | 12 | `Poly` prunes cancelled terms, so the count is of genuinely non-zero monomials |
//! | C2 | `pencil_matches` | 3 | `PAIRINGS` is the three ways to split a cube into opposite faces, and there are three |
//! | C3 | `random_rational_trials - ratio_pairs_equal_trials` | `random_rational_trials` | every trial is a reduced-pair comparison, decided or not at all |
//!
//! Concretely, what C1 and C2 buy is that `R-128` may assert the `GL(2)^3`
//! weight is a *square* and therefore that the body-saddle count cannot depend
//! on cell aspect ratio; that `R-129` may replace an octahedral sweep with
//! algebra; that `R-130` may read tensor rank off the sign via de Silva & Lim;
//! and that `R-133`'s exact sign has something to be exact *about*. All four are
//! wave-2 rows and none of them is measured here.
//!
//! # Vacuity controls
//!
//! `M-44`: a zero that could not have been non-zero is not a measurement. Every
//! control runs before the first `run.record` and every panic message starts
//! `VOID: `.
//!
//! - **`positive_trials > 0` and `negative_trials > 0`** — the registration's own
//!   words, *"the trial set must contain 8-tuples of both discriminant signs"*.
//!   The `f32` disagreement count is a claim about signs; a one-signed population
//!   could report it without ever having had a sign to get wrong.
//! - **`near_zero_trials >= 50`** — the registration's bar verbatim. Without it
//!   C3 samples only the stratum where `|Delta|` dwarfs `f32` epsilon by twelve
//!   orders, `f32_sign_disagreements` reads 0 for a reason that has nothing to do
//!   with the arithmetic, and the falsifier fires on an artefact of the fixture.
//! - **`random_rational_trials >= 3000`** — the registered population size.
//! - **`inexact_f32_inputs == 0`** — every corner value round-trips `f32 -> f64`
//!   exactly, so a sign disagreement is the *evaluation* and not the input. This
//!   is the control that makes `f32_sign_disagreements` the number `P-134` can
//!   act on rather than a restatement of `f32`'s mantissa width.
//! - **`corrupted_control_residual_terms == 1`** — Cayley with its `f0^2*f7^2`
//!   coefficient bumped from 1 to 2 must differ from the repo discriminant by
//!   exactly one monomial. This is the control on the predicate itself:
//!   `sub(..).is_zero()` carries C1 *and* all of C2, and a predicate that
//!   answered `true` unconditionally would pass every clause in the row.
//! - **each pencil is non-zero with twelve terms** — recorded per row as
//!   `pencil_terms`, so `pencil_matches` counting to 3 cannot be three
//!   comparisons of the zero polynomial against itself.
//!
//! # Determinism
//!
//! One thread. No wall clock gates anything: `wall_ns` is recorded because it is
//! interesting and is read by nothing, which is the only safe status for a
//! nanosecond on a host whose governor swings the same binary 1.45x (`M-280`).
//! The PRNG is `common::poly::Rng`, a SplitMix64 seeded by [`SEED`] and recorded
//! in the `seed` column, so the disagreement counts are the same on every host
//! and every re-run — a count that moves between runs is not a measurement.
//! Every polynomial is a `BTreeMap`, so `expression` is byte-identical across
//! runs; every exact comparison is on `i128` integers; every float sign goes
//! through [`float_sign`], which maps both `+0.0` and `-0.0` to `0` because
//! `trilinear.rs:250` does.

#![allow(
    clippy::float_cmp,
    reason = "C3's whole point is exact float equality: the f32-exactness control asserts \
              f64::from(x as f32) == x, and the sign test compares against exactly zero \
              because trilinear.rs:250 does"
)]

mod common;

use std::cmp::Ordering;
use std::time::Instant;

use crate::common::poly::{self, Poly};

// ─── the registered constants ───────────────────────────────────────────────

/// The PRNG seed, recorded in the `seed` column.
///
/// Any value would do; this one is fixed so the counts are reproducible, which
/// is the only property a seed needs to have.
const SEED: u64 = 0x0127_5AD1_E0C4_A1E9;

/// Trials in the general stratum. The registration's floor is 3,000 over the
/// whole population, and the general stratum alone meets it.
const GENERAL_TRIALS: usize = 3_000;

/// Attempts at a constructed near-zero trial. Not every draw brackets a root in
/// `[-4, 4]` — a quadratic with complex roots has none to bracket — so this is
/// an attempt count and `near_zero_brackets_found` is the yield.
const NEAR_ZERO_ATTEMPTS: usize = 800;

/// The registration's near-zero band: `|Delta| < 1e-6`.
const NEAR_ZERO_THRESHOLD: f64 = 1e-6;

/// The registration's near-zero floor: at least this many.
const NEAR_ZERO_MIN: usize = 50;

/// The registration's population floor.
const MIN_TRIALS: usize = 3_000;

/// Both sides of C1 must have exactly this many non-zero monomials.
const EXPECTED_TERMS: usize = 12;

/// Both sides of C1 must have exactly this total degree.
const EXPECTED_DEGREE: u32 = 4;

/// The general stratum's numerator span: `f_i` in `[-1024, 1024]`.
///
/// Bounded so `eval_ratio`'s cleared numerator stays inside `i128` with room to
/// spare: the worst monomial is `4 * num^4 * den^12 <= 2^2 * 2^40 * 2^36`, and
/// twelve of them sum below `2^82`.
const GENERAL_NUM_SPAN: i64 = 1_024;

/// The general stratum's denominators are `2^0 .. 2^3`, i.e. `1, 2, 4, 8`.
///
/// Powers of two rather than arbitrary integers because decision 1 in the header
/// requires every corner value to be exactly an `f32`, and a dyadic rational
/// with a small numerator is.
const GENERAL_DEN_LOG2_MAX: i64 = 3;

/// The near-zero stratum draws `f0..f6` at spacing `2^-4`, in `[-1, 1]`.
const NEAR_ZERO_DEN_LOG2: u32 = 4;

/// The near-zero stratum scans `f7` over `[-4, 4]` at the same spacing, so
/// `[-BRACKET_SPAN, BRACKET_SPAN]` in units of `2^-NEAR_ZERO_DEN_LOG2`.
const BRACKET_SPAN: i64 = 64;

/// `f32` carries 24 significant bits, so `n / 2^k` is exactly an `f32` while
/// `|n| < 2^24`. Bisection stops at the last halving that keeps it there — the
/// cap is derived from this, never chosen.
const F32_MANTISSA_LIMIT: i64 = 1 << 24;

// ─── exact and float signs ──────────────────────────────────────────────────

/// The sign of an exact numerator, as `-1`, `0` or `1`.
///
/// `Poly::eval_ratio` returns a reduced pair with a positive denominator, so the
/// numerator's sign *is* the rational's sign.
fn exact_sign(n: i128) -> i32 {
    match n.cmp(&0) {
        Ordering::Greater => 1,
        Ordering::Less => -1,
        Ordering::Equal => 0,
    }
}

/// The sign of a float, with **both** zeros mapping to `0`.
///
/// Not `total_cmp` against `0.0`, which orders `-0.0` below `+0.0` and would
/// report a negative sign for a value the extractor treats as zero:
/// `trilinear.rs:250` tests `discriminant == R::ZERO`, and `-0.0 == 0.0`. An
/// `f32` that flushed to zero has *lost* the sign, which is why that case is
/// counted as a disagreement rather than resolved into one.
fn float_sign(x: f64) -> i32 {
    if x == 0.0 {
        return 0;
    }
    if x.is_sign_negative() { -1 } else { 1 }
}

// ─── one trial ──────────────────────────────────────────────────────────────

/// One 8-tuple of corner values, as exact dyadic rationals `num[i] / den[i]`.
///
/// Kept as a numerator/denominator pair rather than as floats because the exact
/// arm needs `eval_ratio`'s `i64` inputs and the float arms need a value that is
/// *representable*, and one dyadic pair is both.
#[derive(Clone, Copy, Debug)]
struct Trial {
    /// Numerators, one per corner.
    num: [i64; poly::VARS],
    /// Denominators, one per corner. Always a power of two.
    den: [i64; poly::VARS],
    /// Whether the bisection built this trial, as opposed to the general draw.
    constructed: bool,
    /// `log2` of `den[7]` — the bisection depth reached, `0` for a general draw.
    depth: u32,
}

impl Trial {
    /// The corner values in `f64`. Exact: a small numerator over a power of two.
    fn as_f64(&self) -> [f64; poly::VARS] {
        std::array::from_fn(|i| self.num[i] as f64 / self.den[i] as f64)
    }

    /// The corner values in `f32`. Exact for the same reason, and
    /// [`Trial::inputs_are_f32_exact`] is the check that it really is.
    fn as_f32(&self) -> [f32; poly::VARS] {
        std::array::from_fn(|i| self.num[i] as f32 / self.den[i] as f32)
    }

    /// Whether every corner value survives `f64 -> f32 -> f64` unchanged.
    ///
    /// Measured, not argued: this is the control that makes a `f32` sign
    /// disagreement a fact about the *arithmetic*.
    fn inputs_are_f32_exact(&self) -> bool {
        let wide = self.as_f64();
        let narrow = self.as_f32();
        (0..poly::VARS).all(|i| wide[i] == f64::from(narrow[i]))
    }
}

/// A general-stratum draw: eight independent dyadic rationals.
fn draw_general(rng: &mut poly::Rng) -> Trial {
    let mut num = [0i64; poly::VARS];
    let mut den = [1i64; poly::VARS];
    for (n, d) in num.iter_mut().zip(den.iter_mut()) {
        *n = rng.next_i64_in(-GENERAL_NUM_SPAN, GENERAL_NUM_SPAN + 1);
        *d = 1i64 << rng.next_i64_in(0, GENERAL_DEN_LOG2_MAX + 1);
    }
    Trial {
        num,
        den,
        constructed: false,
        depth: 0,
    }
}

/// Assemble a near-zero trial: `f0..f6 = base[i] / base_den`, `f7 = n7 / 2^k7`.
fn trial_from(base: &[i64; 7], base_den: i64, n7: i64, k7: u32) -> Trial {
    let mut num = [0i64; poly::VARS];
    let mut den = [base_den; poly::VARS];
    num[..7].copy_from_slice(base);
    num[7] = n7;
    den[7] = 1i64 << k7;
    Trial {
        num,
        den,
        constructed: true,
        depth: k7,
    }
}

/// The exact discriminant at one `f7` candidate: the trial, its sign, its value.
fn eval_at(p: &Poly, base: &[i64; 7], base_den: i64, n7: i64, k7: u32) -> (Trial, i32, f64) {
    let trial = trial_from(base, base_den, n7, k7);
    let (n, d) = p.eval_ratio(&trial.num, &trial.den);
    (trial, exact_sign(n), n as f64 / d as f64)
}

/// What one bisection attempt produced.
#[derive(Clone, Copy, Debug)]
enum Attempt {
    /// No sign change over the scanned range: nothing to bisect. A quadratic in
    /// `f7` with complex roots has none, which is an ordinary outcome and not a
    /// failure of the method.
    NoBracket,
    /// A bracket was found and bisected. Carries the best point seen; its
    /// discriminant is re-evaluated by `Tally::score`, which is the one place
    /// every trial's exact value is read, so the value is not carried out.
    Bisected(Trial),
}

/// Drive the exact discriminant toward zero by bisecting one corner.
///
/// `f0..f6` are drawn in `[-1, 1]` at spacing `2^-NEAR_ZERO_DEN_LOG2`; the
/// discriminant is then a quadratic in `f7` alone, so a sign change between two
/// consecutive scan points brackets a root. Every evaluation is the **exact**
/// rational value, so the bracket cannot be lost to rounding, and every
/// evaluation is counted into `evals`.
///
/// Stops at the first of: `|Delta| < NEAR_ZERO_THRESHOLD`, an exact root, or the
/// last halving that keeps `|numerator| < 2^24` and the trial exactly
/// representable in `f32`.
fn draw_near_zero(disc: &Poly, rng: &mut poly::Rng, evals: &mut u64) -> Attempt {
    let base_den = 1i64 << NEAR_ZERO_DEN_LOG2;
    let mut base = [0i64; 7];
    for b in &mut base {
        *b = rng.next_i64_in(-base_den, base_den + 1);
    }

    // Scan for a sign change. An exact root on the scan grid is already inside
    // the band and is taken as-is.
    let mut previous: Option<(i64, i32)> = None;
    let mut bracket: Option<(i64, i64)> = None;
    for j in -BRACKET_SPAN..=BRACKET_SPAN {
        let (trial, sign, _) = eval_at(disc, &base, base_den, j, NEAR_ZERO_DEN_LOG2);
        *evals += 1;
        if sign == 0 {
            return Attempt::Bisected(trial);
        }
        if let Some((previous_j, previous_sign)) = previous
            && previous_sign != sign
        {
            bracket = Some((previous_j, j));
            break;
        }
        previous = Some((j, sign));
    }
    let Some((mut lo, mut hi)) = bracket else {
        return Attempt::NoBracket;
    };

    let mut depth = NEAR_ZERO_DEN_LOG2;
    let (lo_trial, lo_sign, lo_value) = eval_at(disc, &base, base_den, lo, depth);
    let (hi_trial, _, hi_value) = eval_at(disc, &base, base_den, hi, depth);
    *evals += 2;
    let mut best = if lo_value.abs() <= hi_value.abs() {
        (lo_trial, lo_value)
    } else {
        (hi_trial, hi_value)
    };

    while best.1.abs() >= NEAR_ZERO_THRESHOLD
        && lo.abs().max(hi.abs()).saturating_mul(2).saturating_add(1) < F32_MANTISSA_LIMIT
    {
        lo *= 2;
        hi *= 2;
        depth += 1;
        let mid = lo + 1;
        let (trial, sign, value) = eval_at(disc, &base, base_den, mid, depth);
        *evals += 1;
        if value.abs() < best.1.abs() {
            best = (trial, value);
        }
        if sign == 0 {
            break;
        }
        if sign == lo_sign {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    Attempt::Bisected(best.0)
}

// ─── the exact ratio deviation ──────────────────────────────────────────────

/// `ratio - 1` for one trial, as an exact integer ratio.
///
/// The registration asks for `max_abs_ratio_deviation`, which must be exactly
/// zero. Reporting it as a decimal would make an exact zero indistinguishable
/// from a rounded one, so it is reported as `p/q` — `0/1` when every trial's two
/// reduced pairs are identical.
///
/// A denominator of zero is the honest encoding of an *undefined* ratio, which
/// is what a zero Cayley evaluation against a non-zero repo evaluation is. Its
/// magnitude is then infinite and it wins the maximum, which is the right
/// verdict for a trial whose ratio does not exist.
#[derive(Clone, Copy, Debug)]
struct Deviation {
    /// The numerator of `ratio - 1`.
    num: i128,
    /// The denominator. Zero encodes an undefined ratio.
    den: i128,
}

impl Deviation {
    /// Exactly zero: the value the registration predicts on every trial.
    const ZERO: Self = Self { num: 0, den: 1 };

    /// `|num / den|`, used **only** to order deviations against each other. The
    /// value reported is always the exact pair.
    fn magnitude(self) -> f64 {
        (self.num as f64 / self.den as f64).abs()
    }

    /// The CSV token. No comma, no quote, no newline.
    fn text(self) -> String {
        format!("{}/{}", self.num, self.den)
    }
}

// ─── the tally over the whole population ────────────────────────────────────

/// Everything C3 counts, over both strata.
#[derive(Clone, Copy, Debug, Default)]
struct Tally {
    /// Trials scored. The registration's floor is 3,000.
    trials: usize,
    /// Trials from the general draw.
    general: usize,
    /// Trials the bisection constructed.
    constructed: usize,
    /// `|Delta| < NEAR_ZERO_THRESHOLD`, over the whole population.
    near_zero: usize,
    /// Exact sign `+1`.
    positive: usize,
    /// Exact sign `-1`.
    negative: usize,
    /// Exact value identically zero.
    zero: usize,
    /// Trials whose two reduced pairs were identical, i.e. ratio exactly 1.
    pairs_equal: usize,
    /// Trials where the Cayley evaluation was zero and the repo's was not.
    undefined_ratio: usize,
    /// The largest deviation seen, exactly.
    worst: Option<Deviation>,
    /// `f64` sign of the repo discriminant differs from the exact sign.
    f64_disagree_disc: usize,
    /// `f64` sign of Cayley differs from the exact sign.
    f64_disagree_cayley: usize,
    /// Either `f64` evaluation disagrees. This is `f64_sign_disagreements`.
    f64_disagree: usize,
    /// `f32` sign of the repo discriminant differs from the exact sign.
    f32_disagree_disc: usize,
    /// `f32` sign of Cayley differs from the exact sign.
    f32_disagree_cayley: usize,
    /// Either `f32` evaluation disagrees. This is `f32_sign_disagreements`.
    f32_disagree: usize,
    /// The two `f32` evaluations disagree with **each other**. Structurally zero
    /// once C1 holds, because they are then the same polynomial.
    f32_cross: usize,
    /// The `f32` evaluation was `+-0.0` against a non-zero exact value — the
    /// cell `trilinear.rs:250` would re-classify as a tangency.
    f32_flush: usize,
    /// `f32` disagreements inside the constructed stratum.
    f32_disagree_constructed: usize,
    /// `f32` disagreements inside the general stratum.
    f32_disagree_general: usize,
    /// Trials whose corner values did **not** round-trip `f32` exactly.
    inexact_inputs: usize,
    /// Non-finite float evaluations. Expected zero: the inputs are bounded.
    nonfinite: usize,
    /// Deepest bisection reached, as `log2(den[7])`.
    max_depth: u32,
    /// `eval_ratio` calls, bisection included.
    exact_evals: u64,
}

impl Tally {
    /// Score one trial against both polynomials, exactly and in both precisions.
    fn score(&mut self, trial: &Trial, disc: &Poly, cayley: &Poly) {
        self.trials += 1;
        if trial.constructed {
            self.constructed += 1;
            self.max_depth = self.max_depth.max(trial.depth);
        } else {
            self.general += 1;
        }
        if !trial.inputs_are_f32_exact() {
            self.inexact_inputs += 1;
        }

        // ── exact ──
        let (disc_num, disc_den) = disc.eval_ratio(&trial.num, &trial.den);
        let (cayley_num, cayley_den) = cayley.eval_ratio(&trial.num, &trial.den);
        self.exact_evals += 2;
        let sign = exact_sign(disc_num);
        let value = disc_num as f64 / disc_den as f64;

        match sign {
            1 => self.positive += 1,
            -1 => self.negative += 1,
            _ => self.zero += 1,
        }
        if value.abs() < NEAR_ZERO_THRESHOLD {
            self.near_zero += 1;
        }

        let deviation = if (disc_num, disc_den) == (cayley_num, cayley_den) {
            self.pairs_equal += 1;
            Deviation::ZERO
        } else if cayley_num == 0 {
            self.undefined_ratio += 1;
            Deviation {
                num: i128::from(sign),
                den: 0,
            }
        } else {
            let cross = disc_num
                .checked_mul(cayley_den)
                .and_then(|left| {
                    cayley_num
                        .checked_mul(disc_den)
                        .and_then(|right| left.checked_sub(right))
                })
                .expect("the exact ratio deviation's numerator fits i128");
            let den = cayley_num
                .checked_mul(disc_den)
                .expect("the exact ratio deviation's denominator fits i128");
            Deviation { num: cross, den }
        };
        let beats = self
            .worst
            .is_none_or(|current| deviation.magnitude() > current.magnitude());
        if beats {
            self.worst = Some(deviation);
        }

        // ── f64 and f32, on exactly representable inputs ──
        let wide = trial.as_f64();
        let narrow = trial.as_f32();
        let disc_f64 = disc.eval_f64(&wide);
        let cayley_f64 = cayley.eval_f64(&wide);
        let disc_f32 = disc.eval_f32(&narrow);
        let cayley_f32 = cayley.eval_f32(&narrow);
        if !disc_f64.is_finite()
            || !cayley_f64.is_finite()
            || !disc_f32.is_finite()
            || !cayley_f32.is_finite()
        {
            self.nonfinite += 1;
        }

        let disc_f64_sign = float_sign(disc_f64);
        let cayley_f64_sign = float_sign(cayley_f64);
        let disc_f32_sign = float_sign(f64::from(disc_f32));
        let cayley_f32_sign = float_sign(f64::from(cayley_f32));

        if disc_f64_sign != sign {
            self.f64_disagree_disc += 1;
        }
        if cayley_f64_sign != sign {
            self.f64_disagree_cayley += 1;
        }
        if disc_f64_sign != sign || cayley_f64_sign != sign {
            self.f64_disagree += 1;
        }
        if disc_f32_sign != sign {
            self.f32_disagree_disc += 1;
        }
        if cayley_f32_sign != sign {
            self.f32_disagree_cayley += 1;
        }
        if disc_f32_sign != cayley_f32_sign {
            self.f32_cross += 1;
        }
        if disc_f32_sign == 0 && sign != 0 {
            self.f32_flush += 1;
        }
        if disc_f32_sign != sign || cayley_f32_sign != sign {
            self.f32_disagree += 1;
            if trial.constructed {
                self.f32_disagree_constructed += 1;
            } else {
                self.f32_disagree_general += 1;
            }
        }
    }

    /// The largest deviation, or exact zero if nothing was scored.
    fn worst_deviation(&self) -> Deviation {
        self.worst.unwrap_or(Deviation::ZERO)
    }
}

// ─── the pairing name, split into its two readable halves ───────────────────

/// The corner-set token of a pairing name, e.g. `0123|4567`.
fn pairing_token(name: &str) -> &str {
    name.rsplit(' ')
        .next()
        .expect("a pairing name always has a last token")
}

/// The axis half of a pairing name, e.g. `w-slices`.
fn pairing_axis(name: &str) -> &str {
    name.split(' ')
        .next()
        .expect("a pairing name always has a first token")
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-127");

    common::experiment::run(prereg, |run| {
        // ── C1: the symbolic identity, exactly ──────────────────────────────
        let disc = poly::repo_discriminant();
        let cayley = poly::cayley_2x2x2();
        let residual = disc.sub(&cayley);

        let terms_disc = disc.terms();
        let terms_cayley = cayley.terms();
        let total_degree = disc.total_degree();
        let difference_is_zero = residual.is_zero();
        let c1 = difference_is_zero
            && terms_disc == EXPECTED_TERMS
            && terms_cayley == EXPECTED_TERMS
            && total_degree == EXPECTED_DEGREE
            && cayley.total_degree() == EXPECTED_DEGREE;

        // The expanded form as one CSV-safe token: `Poly`'s own `Display` with
        // the spaces removed, so the string is derived from the algebra rather
        // than transcribed beside it.
        let expression = disc.to_string().replace(' ', "");
        let expression_cayley = cayley.to_string().replace(' ', "");

        // ── C2: the three pencils ───────────────────────────────────────────
        let pencils: [Poly; 3] = std::array::from_fn(poly::pencil_discriminant);
        let pencil_residual: [Poly; 3] = std::array::from_fn(|p| pencils[p].sub(&cayley));
        let pencil_agrees: [bool; 3] = std::array::from_fn(|p| pencil_residual[p].is_zero());
        let pencil_total = pencil_agrees.iter().filter(|ok| **ok).count();
        let c2 = pencil_total == pencils.len();

        // ── the controls on the predicate itself, before any measurement ────
        //
        // `sub(..).is_zero()` carries C1 and all three of C2's comparisons. A
        // predicate that could only answer `true` would pass both clauses
        // without reading anything, so it is shown answering `false` first.
        assert!(
            !disc.is_zero(),
            "VOID: the repo discriminant expanded to the zero polynomial, so C1's \
             zero difference would be two zeros agreeing and would say nothing \
             about trilinear.rs:246"
        );
        assert!(
            !cayley.is_zero(),
            "VOID: Cayley's form expanded to the zero polynomial, so C1 and all of \
             C2 would compare against nothing"
        );
        let mut bump = [0u8; poly::VARS];
        bump[0] = 2;
        bump[7] = 2;
        assert_eq!(
            disc.coefficient(bump),
            1,
            "P-127: the f0^2*f7^2 coefficient pins the normalisation at c1^2 - 4*c0*c2 \
             with no leading sign; a value other than 1 means the transcription or the \
             corner indexing moved"
        );
        let corrupted = cayley.add(&Poly::monomial(bump, 1));
        let control_residual = disc.sub(&corrupted);
        assert!(
            !control_residual.is_zero(),
            "VOID: bumping Cayley's f0^2*f7^2 coefficient from 1 to 2 left the \
             difference zero, so Poly::sub(..).is_zero() cannot read false and \
             neither C1 nor C2 is falsifiable by it"
        );
        assert_eq!(
            control_residual.terms(),
            1,
            "VOID: the corrupted control must differ by exactly the one monomial that \
             was corrupted, got {} terms",
            control_residual.terms()
        );
        for (p, pencil) in pencils.iter().enumerate() {
            assert!(
                !pencil.is_zero() && pencil.terms() == EXPECTED_TERMS,
                "VOID: pencil {p} ({}) expanded to {} terms, so its agreement with \
                 Cayley would not be a comparison of two twelve-term forms",
                poly::PAIRING_NAMES[p],
                pencil.terms()
            );
        }

        // ── C3: the exact-rational population ───────────────────────────────
        let started = Instant::now();
        let mut rng = poly::Rng::new(SEED);
        let mut tally = Tally::default();
        let mut bracket_evals: u64 = 0;
        let mut brackets_found = 0usize;

        for _ in 0..GENERAL_TRIALS {
            let trial = draw_general(&mut rng);
            tally.score(&trial, &disc, &cayley);
        }
        for _ in 0..NEAR_ZERO_ATTEMPTS {
            match draw_near_zero(&disc, &mut rng, &mut bracket_evals) {
                Attempt::NoBracket => {}
                Attempt::Bisected(trial) => {
                    brackets_found += 1;
                    tally.score(&trial, &disc, &cayley);
                }
            }
        }
        tally.exact_evals += bracket_evals;
        let wall_ns = started.elapsed().as_nanos();

        let worst = tally.worst_deviation();
        let ratio_exactly_one = tally.pairs_equal == tally.trials;
        let c3 = tally.trials >= MIN_TRIALS
            && ratio_exactly_one
            && tally.f64_disagree == 0
            && tally.f32_disagree > 0;

        // ── the registration's vacuity controls ─────────────────────────────
        assert!(
            tally.trials >= MIN_TRIALS,
            "VOID: {} trials against the registered floor of {MIN_TRIALS}, so C3's \
             agreement is claimed over a smaller population than was registered",
            tally.trials
        );
        assert!(
            tally.positive > 0,
            "VOID: no trial had a positive discriminant, so the sign-disagreement \
             counts are claims about signs over a one-signed population (the \
             registration: 'the trial set must contain 8-tuples of both \
             discriminant signs')"
        );
        assert!(
            tally.negative > 0,
            "VOID: no trial had a negative discriminant, so the sign-disagreement \
             counts are claims about signs over a one-signed population (the \
             registration: 'the trial set must contain 8-tuples of both \
             discriminant signs')"
        );
        assert!(
            tally.near_zero >= NEAR_ZERO_MIN,
            "VOID: only {} trials inside |Delta| < {NEAR_ZERO_THRESHOLD:e} against \
             the registered floor of {NEAR_ZERO_MIN}, so C3 is sampling only the easy \
             stratum where |Delta| dwarfs f32 epsilon and a zero f32 disagreement \
             count would be an artefact of the fixture ({} brackets found in \
             {NEAR_ZERO_ATTEMPTS} attempts)",
            tally.near_zero,
            brackets_found
        );
        assert_eq!(
            tally.inexact_inputs, 0,
            "VOID: {} trials carried a corner value that does not round-trip f32 \
             exactly, so a sign disagreement could be the input having been rounded \
             before any arithmetic happened rather than the arithmetic getting the \
             sign wrong -- which is not the number P-134 exists to act on",
            tally.inexact_inputs
        );
        assert_eq!(
            tally.nonfinite, 0,
            "VOID: {} float evaluations were not finite, so a sign comparison against \
             them is a comparison against an infinity and the disagreement count is \
             measuring overflow rather than cancellation",
            tally.nonfinite
        );

        // ── three rows, one per axis pairing ────────────────────────────────
        let mut running = 0usize;
        for (p, pencil) in pencils.iter().enumerate() {
            let agrees = pencil_agrees[p];
            running += usize::from(agrees);
            let name = poly::PAIRING_NAMES[p];
            run.record(&[
                // ── the registered columns, in registration order ──────────
                ("expression", expression.clone()),
                ("terms_disc", terms_disc.to_string()),
                ("terms_cayley", terms_cayley.to_string()),
                ("total_degree", total_degree.to_string()),
                (
                    "symbolic_difference_is_zero",
                    difference_is_zero.to_string(),
                ),
                ("pencil_axis_pairing", pairing_token(name).to_string()),
                ("pencil_matches", running.to_string()),
                ("random_rational_trials", tally.trials.to_string()),
                ("max_abs_ratio_deviation", worst.text()),
                ("f32_sign_disagreements", tally.f32_disagree.to_string()),
                ("f64_sign_disagreements", tally.f64_disagree.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // ── extras (M-273) ─────────────────────────────────────────
                //
                // C1's algebra, in enough detail that a non-zero residual is
                // readable from the file rather than only from a re-run.
                ("expression_cayley", expression_cayley.clone()),
                (
                    "expression_matches",
                    (expression == expression_cayley).to_string(),
                ),
                ("symbolic_residual_terms", residual.terms().to_string()),
                ("total_degree_cayley", cayley.total_degree().to_string()),
                (
                    "degree_per_corner",
                    (0..poly::VARS)
                        .map(|i| disc.degree_in(i).to_string())
                        .collect::<Vec<_>>()
                        .join("|"),
                ),
                ("disc_is_multi_affine", disc.is_multi_affine().to_string()),
                ("coefficient_f0sq_f7sq", disc.coefficient(bump).to_string()),
                (
                    "corrupted_control_residual_terms",
                    control_residual.terms().to_string(),
                ),
                // C2, per row and in total.
                ("pencil_split_axis", pairing_axis(name).to_string()),
                ("pencil_matches_here", agrees.to_string()),
                ("pencil_terms", pencil.terms().to_string()),
                ("pencil_total_degree", pencil.total_degree().to_string()),
                (
                    "pencil_residual_terms",
                    pencil_residual[p].terms().to_string(),
                ),
                ("pencil_pairings_checked", pencils.len().to_string()),
                ("pencil_matches_total", pencil_total.to_string()),
                // C3's population, and the strata it is made of.
                ("seed", format!("{SEED:#018x}")),
                ("general_trials", tally.general.to_string()),
                ("near_zero_trials", tally.near_zero.to_string()),
                ("near_zero_constructed", tally.constructed.to_string()),
                ("near_zero_brackets_found", brackets_found.to_string()),
                ("near_zero_attempts", NEAR_ZERO_ATTEMPTS.to_string()),
                ("near_zero_floor", NEAR_ZERO_MIN.to_string()),
                ("near_zero_threshold", format!("{NEAR_ZERO_THRESHOLD:e}")),
                ("positive_trials", tally.positive.to_string()),
                ("negative_trials", tally.negative.to_string()),
                ("zero_trials", tally.zero.to_string()),
                ("bisection_max_depth", tally.max_depth.to_string()),
                ("f32_mantissa_limit", F32_MANTISSA_LIMIT.to_string()),
                // C3's exact half.
                ("ratio_pairs_equal_trials", tally.pairs_equal.to_string()),
                ("ratio_exactly_one", ratio_exactly_one.to_string()),
                ("undefined_ratio_trials", tally.undefined_ratio.to_string()),
                (
                    "max_abs_ratio_deviation_f64",
                    format!("{:.6e}", worst.magnitude()),
                ),
                ("exact_evals", tally.exact_evals.to_string()),
                // C3's float half, split so P-134 inherits a decomposition
                // rather than one number.
                (
                    "f32_sign_disagreements_disc",
                    tally.f32_disagree_disc.to_string(),
                ),
                (
                    "f32_sign_disagreements_cayley",
                    tally.f32_disagree_cayley.to_string(),
                ),
                ("f32_cross_disagreements", tally.f32_cross.to_string()),
                ("f32_flush_to_zero_trials", tally.f32_flush.to_string()),
                (
                    "f32_disagreements_near_zero",
                    tally.f32_disagree_constructed.to_string(),
                ),
                (
                    "f32_disagreements_general",
                    tally.f32_disagree_general.to_string(),
                ),
                (
                    "f64_sign_disagreements_disc",
                    tally.f64_disagree_disc.to_string(),
                ),
                (
                    "f64_sign_disagreements_cayley",
                    tally.f64_disagree_cayley.to_string(),
                ),
                ("inexact_f32_inputs", tally.inexact_inputs.to_string()),
                ("nonfinite_evaluations", tally.nonfinite.to_string()),
                // Time, recorded beside the verdicts and read by nothing.
                ("wall_ns", wall_ns.to_string()),
            ]);
        }
    });
}
