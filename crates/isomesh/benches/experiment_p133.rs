//! **P-133 — an exact sign for the body-saddle discriminant, and how often `f32`
//! gets it wrong.**
//!
//! Ticket: R-133. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p133
//! ```
//!
//! Writes `docs/experiments/p-133.csv`.
//!
//! # What was missing
//!
//! `crates/isomesh/src/marching_cubes/trilinear.rs:246` is
//! `b * b - R::TWO * R::TWO * a * c` over the three coefficients built at
//! `:199-214`, and its **sign** is the only thing the crate reads from it:
//! `:247` returns no roots when it is negative, `:250` returns one when it is
//! zero, and `:263-266` returns two otherwise. That root count is what
//! `BodySaddles::inside_mask` is built from, what `Contours::topology` branches
//! on, and therefore what decides whether an ambiguous cell is meshed as a
//! tunnel or as disks.
//!
//! Nothing in the repository asks whether that sign is *right*.
//!
//! - **`M-207`** defends the root *extraction* — Kahan's form instead of
//!   `(-b ± √d)/2a`, and `a == 0` as a one-root case rather than an infinity. It
//!   is a finding about cancellation in `q`, downstream of the branch this row
//!   audits.
//! - **`M-206`** records that `interior::SweptFaces` and `trilinear::BodySaddles`
//!   locate the same saddles to `1.1e-12`. `1.1e-12` is a *relative* agreement
//!   between two `f64` evaluations; it says nothing about either one's sign near
//!   a root, which is exactly where a relative agreement is worth nothing.
//! - **`M-214`–`M-217`** all build on the saddle count.
//! - **`P-127`/`M-440`** proved the identity this row needs: `b*b - 4*a*c` from
//!   `coefficients` **is** Cayley's `2×2×2` hyperdeterminant, twelve terms,
//!   total degree 4, symbolic difference identically zero, and
//!   `docs/experiments/p-127.csv` records `f32_sign_disagreements = 14` against
//!   `f64_sign_disagreements = 0` over 3,481 dyadic rational trials. That is a
//!   *synthetic* population, deliberately bisected onto roots of the
//!   discriminant. It licenses this row and does not answer it: nobody has asked
//!   whether the eight reference fields, on the grids the crate actually
//!   marches, put cells in that stratum at all.
//! - **`real.rs:81-104`** already carries `UNIT_ROUNDOFF` and the
//!   Dekker–Veltkamp `SPLITTER` and cites Shewchuk `10.1007/pl00009321` in their
//!   doc comments — the expansion arithmetic has had its constants shipped for
//!   its entire life and no expansion has ever been formed.
//!
//! So the discriminant has been the crate's one topological branch on a
//! near-degenerate quantity, and its error has never been measured on a real
//! field.
//!
//! # The exact arithmetic, and why it is exact
//!
//! The corner values are `f64`, so the exact value of a degree-4 form in them is
//! **not** an integer and `Poly::eval_i128` cannot be used: four 53-bit
//! significands multiply to 212 bits, a coefficient of 4 makes 214, and twelve
//! of those summed makes 218 — past `i128` before the exponents are even
//! considered.
//!
//! This harness uses **Shewchuk's expansion arithmetic** (`10.1007/pl00009321`,
//! the DOI `real.rs:94` already cites), which is exact for arbitrary exponent
//! spread and grows with the *term count* rather than with the exponent range:
//!
//! - `two_sum` is Knuth's six-flop exact sum; `split` is Dekker–Veltkamp at
//!   `<f64 as Real>::SPLITTER`; `two_product` is the FMA-free exact product.
//! - `scale_expansion` is `scale_expansion_zeroelim`: an exact expansion times
//!   one `f64`, output length at most twice the input's.
//! - `expansion_sum` is `fast_expansion_sum_zeroelim`: the merge of two
//!   non-overlapping increasing expansions.
//! - Shewchuk's `Fast_Two_Sum` is **not** used anywhere; `two_sum` replaces it
//!   at every site. It is exact without the `|a| ≥ |b|` precondition, so the
//!   correctness of this file does not depend on an ordering argument.
//!
//! The polynomial is not retyped. Its twelve monomials come from
//! `common::poly::repo_discriminant()` — the module `P-127` owns, transcribed
//! line for line from `coefficients` — via `Poly::monomials()`. Each monomial is
//! homogeneous of degree 4, so it is a product of exactly four corner values:
//! `two_product` of the first two gives a 2-component expansion, two
//! `scale_expansion` calls take it to at most 8, a third scales by the integer
//! coefficient to at most 16, and the twelve are merged into one expansion of at
//! most 192 components. Its sign is the sign of its largest component, because a
//! non-overlapping increasing expansion cannot be outvoted by its own tail.
//!
//! Three things this depends on, all recorded rather than assumed:
//!
//! 1. **No FP contraction and no reassociation.** `rustc` performs neither
//!    without `-ffast-math`, which the workspace does not set; `sdf.rs:155-160`
//!    already leans on the same guarantee for its golden hashes.
//! 2. **Homogeneity of degree 4**, asserted over every monomial — the four-factor
//!    chain would be silently wrong on a degree-3 or degree-5 term.
//! 3. **No overflow and no underflow in the splits.** `max_abs_corner` and
//!    `min_abs_nonzero_corner` are recorded, and the safe window is asserted.
//!    Every reference field's values sit within a decade of 1 over its own
//!    domain, so this is decades away from `2^±969`; the assert is there because
//!    a proviso nobody checks is a proviso nobody has.
//!
//! # The filter, and where its constant comes from
//!
//! The filter is static: compute `Δ` in the working precision by the crate's own
//! chain, compute the expression's **permanent** `P` — the same tree with every
//! operand replaced by its absolute value and every subtraction by an addition —
//! and certify the sign when `|Δ̂| > K·u·P̂`.
//!
//! `K` is a rounding depth, not a guess. With `d(x)` the depth of the longest
//! rounding chain, `d(x·y) = 1 + d(x) + d(y)` and `d(x ± y) = 1 + max(d(x), d(y))`,
//! and `|x̂ − x| ≤ ((1+u)^{d(x)} − 1)·P(x)` by induction on the tree. Over
//! `trilinear.rs:202-213` plus `:246`:
//!
//! | node | depth | permanent |
//! |---|---|---|
//! | `f_i` | 0 | `|f_i|` |
//! | `twist_lo = (f0+f3) − (f1+f2)` | 2 | `Σ|f_0..3|` |
//! | `du_lo`, `du_hi`, `dv_lo`, `dv_hi` | 1 | e.g. `|f1|+|f0|` |
//! | `a = du_hi·twist_lo − du_lo·twist_hi` | 5 | `P_duhi·P_tlo + P_dulo·P_thi` |
//! | `b` | 5 | `(|f4|P_tlo + |f0|P_thi) + (P_duhi·P_dvlo + P_dulo·P_dvhi)` |
//! | `c = f2·f4 − f0·f6` | 2 | `|f2||f4| + |f0||f6|` |
//! | `b·b` | 11 | `P_b²` |
//! | `4·a·c` | 9 | `4·P_a·P_c` |
//! | `Δ = b·b − 4ac` | **12** | `P_b² + 4·P_a·P_c` |
//!
//! `(1+u)^12 − 1 < 12.01u`. The permanent is itself computed in floating point,
//! which inflates it by at most a further `(1+u)^12`, so `13u·P̂` is already a
//! valid bound; this file uses **`26u`**, a factor of two of headroom, because
//! `26u` is `2.9e-15` relative at `f64` and `1.55e-6` at `f32` — a filter that
//! only gives up after fifteen digits of cancellation, so the headroom costs
//! essentially nothing in `filtered_fallback_rate` and buys a bound no reader
//! has to re-derive to trust.
//!
//! The derivation is not trusted either: `filter_certified_wrong` counts cells
//! where the filter certified a sign the exact expansion contradicts, on every
//! cell of every field at both precisions, and it is asserted zero.
//!
//! ## What C2 is up against, counted rather than clocked
//!
//! There is nothing to share between `Δ` and `P`. Shewchuk's `orient2d` gets its
//! filter almost free because its permanent is `|detleft| + |detright|` over the
//! two products the value already formed; here `twist_lo` is
//! `(f0+f3) − (f1+f2)` and `P_tlo` is `|f0|+|f1|+|f2|+|f3|`, which are different
//! numbers, and the same is true of every interior node. So the filter is a
//! second tree, not an annotation on the first, and the arithmetic is decidable
//! before any clock runs:
//!
//! | expression | operations |
//! |---|---|
//! | `coefficients` (`:202-213`) | 10 twists and differences + 3 for `a` + 7 for `b` + 3 for `c` = **23** |
//! | `b*b − 4*a*c` (`:246`) | **4** |
//! | naive predicate | **27** |
//! | permanent | 8 `abs` + 27 arithmetic = **35** |
//! | filtered predicate | 27 + 35 + `abs` + compare + scale = **65** |
//!
//! `65/27 = 2.41`, and that ratio is a property of the two expression trees on
//! any machine. `naive_ops`, `filter_ops` and `filter_ops_ratio` are recorded so
//! C2's verdict rests on an integer count and not only on wall clocks a `1.45×`
//! governor swings (M-280). The clocks are recorded too, with their medians and
//! their spread, because the registration names `exact_ms`, `float_ms` and
//! `overhead_ratio` and a registration is not amended to suit the instrument.
//!
//! The second thing C2 is up against is not a cost at all. A static filter can
//! certify a **non-zero** sign and nothing else: `|Δ̂| > K·u·P̂` is unsatisfiable
//! when `Δ` is exactly zero, so a field whose trilinear discriminant vanishes
//! identically drives the exact path on every cell no matter how cheap the
//! filter is. `exact_zero` counts those cells, `below_f64_error_bound` counts
//! the fallbacks, and `f64_fallback_is_exactly_the_zero_set` records whether the
//! two sets coincide — which is the sharp form of the question, because a
//! coincidence there says the `f64` filter fails *only* on exact degeneracy and
//! never on cancellation.
//!
//! # The fixture: one corner set, exactly representable in both precisions
//!
//! Each `(field, resolution)` is sampled **once**, in `f64`, then every value is
//! rounded to `f32` and widened back. That shared corner set is bit-identical in
//! both precisions, and it is the single most load-bearing decision here:
//!
//! - `sign_disagreements_f32` then counts cells where the `f32` *arithmetic* gets
//!   the sign wrong, not cells where the *input* was rounded before any
//!   arithmetic happened. Without it the column would be about representation,
//!   the same contamination `P-127` guards with `inexact_f32_inputs`.
//! - The exact sign is the exact `Δ` of that same corner set, so "the sign is
//!   wrong" is a statement about one well-posed question.
//! - The two `MarchingCubes` arms in C3 see the same eight numbers per cell, so
//!   the `f32` and `f64` case indices are **equal by construction** —
//!   `is_inside` is a comparison against zero and the values are identical — and
//!   every triangle-count difference has to come from a precision-dependent
//!   decision downstream of the case.
//!
//! The corner set is served to the extractors by `Tabulated`, a bench-local
//! `Sdf` that trilinearly interpolates it. Every grid domain here has a dyadic
//! `origin` and `cell_size` and at most 65 samples per axis, so
//! `sample_grid`'s `origin + cell_size·x` is exact, `(p − origin)/cell_size`
//! returns the index exactly, and the interpolation parameter is exactly `0` or
//! `1` at every grid point in both precisions. `Tabulated` therefore reproduces
//! the corner set bit-for-bit where the march reads it, and supplies a finite
//! gradient off-grid where `unit_gradient` reads it.
//!
//! # Arms
//!
//! Eight fields × three resolutions × two scalars = 48 rows. `RESOLUTIONS` is
//! `[33, 65, 129]`: 33 is the golden fixture's top resolution (`golden.rs:72`),
//! 65 is the `u64`-word-boundary size gotcha 4 names, and 129 is here because
//! C1's population is a **measure-zero set being sampled**. Cells where `Δ`
//! cancels to within `f32`'s error bound are a codimension-1 stratum of the
//! corner cube, so their count grows with the surface's cell count and with
//! nothing else. A fixture that stopped at 65³ would decide C1 on how many
//! surface cells it happened to draw rather than on the predicate, which is the
//! difference between a result about `Δ` and a result about the fixture.
//! `surface_cell_rate` and `below_f32_bound_beyond_exact_zero` are recorded per
//! row so that dependence is readable rather than argued.
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `scalar = f32` | naive and filtered predicates evaluated at `f32`; filter constant `26·2⁻²⁴` | no |
//! | `scalar = f64` | the same two predicates at `f64`; filter constant `26·2⁻⁵³` | no |
//! | exact expansion | Shewchuk over `repo_discriminant()`'s twelve monomials | it is the oracle both arms are scored against |
//! | `filter_certified_wrong` | the filter's own soundness | **yes**, an assert rather than a clause |
//! | `case_index_disagreements` | that the two extraction arms see one corner set | **yes**, asserted zero |
//!
//! # The clauses, and which are global
//!
//! All three clauses are quantified over the fixture rather than over a row —
//! *"on at least two reference fields"*, *"in aggregate"*, *"on at least one
//! field"* — so `c1_holds`, `c2_holds` and `c3_holds` carry the **global**
//! verdict, identical on all 48 rows, and each row's own answer is beside it as
//! `c1_row`, `c2_row`, `c3_row`. Three identical `c1_holds` are not three
//! measurements.
//!
//! - **C1** holds when at least two distinct fields record
//!   `sign_disagreements_f32 > 0`. `fields_disagreeing_f32` is that count.
//! - **C2** holds when `Σ exact_ms / Σ float_ms` over all 48 rows is under
//!   `1.5`. `c2_aggregate_ratio` is that number; `c2_aggregate_f32` and
//!   `c2_aggregate_f64` decompose it by arm.
//! - **C3** holds when some row records `triangles_changed > 0` or
//!   `nonmanifold_delta != 0`.
//!
//! `disagreement_rate` is **that row's** scalar: `sign_disagreements_f32/cells`
//! on an `f32` row and `sign_disagreements_f64/cells` on an `f64` row. Both
//! counts are on every row because the registration names both.
//!
//! ## What C3's instrument can and cannot isolate
//!
//! `InteriorAmbiguity::Trilinear` is reachable **only** through
//! `FaceAmbiguity::AsymptoticDecider`: `mod.rs:278-292` sets `ambiguous = 0`
//! under `Separate`, and the trilinear branch is guarded on `ambiguous != 0`. So
//! the sign-corrected arm is the golden roster's `marching_cubes+trilinear`
//! configuration (`golden.rs:166-177`) at `f64`, against the same configuration
//! at `f32`, and the `f64` arm *is* the sign-corrected arm exactly to the extent
//! that `sign_disagreements_f64` is zero — which is why that column is
//! registered and reported per row rather than assumed.
//!
//! Two precision-dependent decisions sit inside that configuration, and the
//! difference is attributed rather than claimed: `joined_mask_cells_changed`
//! counts cells where the *face* decider disagrees between precisions, and
//! `inside_mask_cells_changed` / `hexagon_cells_changed` /
//! `disc_sign_cells_changed` count cells where the *body-saddle* stage does.
//! `sign_disagreements_f32_ambiguous` is the subset of C1's count that lands on
//! a cell the trilinear path actually visits, and
//! `disc_sign_cells_changed_ambiguous` is the subset of the two arms' sign
//! disagreements that does.
//!
//! Those last two columns are the ones C3's falsifier actually asks for. Its
//! wording — *"which would mean the sign errors fall on cells whose
//! classification does not reach the output — interesting, and a reason to look
//! at where they do fall"* — anticipates a mesh that changes for a reason other
//! than the sign, and the registered test cannot tell the two apart on its own:
//! the `f64` arm corrects the sign **and** every other rounding in the
//! body-saddle solve. So `c3_holds` is the registered test's answer, verbatim,
//! and `c3_attributed_to_sign` is a separate column reading `true` only when
//! some field with a mesh change also has `disc_sign_cells_changed_ambiguous`
//! above zero. A clause is answered as registered and the mechanism is reported
//! beside it; neither is allowed to stand in for the other.
//!
//! # SHARE, recomputed before the numbers
//!
//! *"C2 moves the body-saddle stage only, whose share of extraction must be
//! reported alongside."* The stage C2 moves is the quadratic solve, which is
//! `BodySaddles::of` — `saddle_stage_ms` times exactly that call over the
//! ambiguous-face cells (median of `REPEATS`), `extract_ms` times the whole
//! `MarchingCubes::extract` at the same precision (median of
//! `EXTRACT_REPEATS`), and `saddle_share` is their quotient. It is reported
//! because SHARE requires it; no clause reads it, and it is a single-digit
//! quotient of two wall clocks on a host whose governor swings `1.45×` (M-280),
//! so `saddle_stage_scatter` and `extract_ms_min`/`extract_ms_max` are beside it.
//!
//! `ambiguous_cells` is the denominator that makes it readable: the stage runs on
//! roughly one cell in two hundred (`mod.rs:57-59`), so a small share is the
//! prediction and not a surprise.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record` and every panic starts `VOID: `.
//!
//! - **The registered one.** At least one field must produce cells with `|Δ|`
//!   below the `f32` error bound. `below_f32_error_bound` is that count per row,
//!   and it is the same population as the `f32` arm's filter fallback, so
//!   `filtered_fallback_rate` on an `f32` row is `below_f32_error_bound/cells`.
//!   Zero everywhere would mean the fixture excludes the near-degenerate stratum
//!   and C1 is measuring a problem it cannot see.
//!
//!   The count on its own is not enough to license C1, and
//!   `below_f32_bound_beyond_exact_zero` is why. A cell where `Δ` is *exactly*
//!   zero is below every error bound at every precision and both arms agree on
//!   it, so a fixture whose whole near-degenerate population is exact zeros
//!   satisfies the letter of the control and still cannot produce a single sign
//!   disagreement. The subtraction is the honest number: cells below the `f32`
//!   bound whose exact `Δ` is non-zero are the ones where `f32` can be wrong and
//!   `f64` right, and it is asserted non-zero globally alongside the registered
//!   count.
//! - **`cells > 0` on every row** — a surface-cell population of zero makes
//!   `disagreement_rate` a zero that could not have been non-zero (M-44).
//! - **`ambiguous_cells > 0` somewhere** — the trilinear path is gated on an
//!   ambiguous face, so with none the body-saddle stage never runs, C3's zero is
//!   structural and `saddle_share` is a division by nothing.
//! - **Both signs present in the exact population** — `exact_positive` and
//!   `exact_negative` must both be non-zero globally, or the oracle is constant
//!   and agreement with it is uninformative.
//! - **`filter_certified_wrong == 0`** — soundness of the derived bound, checked
//!   on every cell at both precisions rather than argued from the table above.
//! - **`case_index_disagreements == 0`** — the shared corner set really is
//!   shared, so C3's triangle delta cannot be a case-index difference wearing a
//!   body-saddle costume.
//! - **`nonfinite_deltas == 0`** — a sign comparison against an infinity or a NaN
//!   would be measuring overflow rather than cancellation.
//!
//! No RNG: the fixture is the eight reference fields on a fixed grid, so there is
//! no seed to state.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::table::{AMBIGUOUS_FACES, is_inside};
use isomesh::marching_cubes::trilinear::BodySaddles;
use isomesh::marching_cubes::{FaceAmbiguity, InteriorAmbiguity, MarchingCubes};
use isomesh::validate::{ValidateConfig, validate};
use isomesh::{MeshBuffer, Real, RuntimeShape3, Shape3, Sdf};

use common::poly;

/// Samples per axis. `33` is the golden fixture's top resolution
/// (`golden.rs:72`); `65` is the `u64`-word-boundary size the cheat sheet's
/// gotcha 4 names; `129` is here because C1's population is a codimension-1
/// stratum being sampled and its count grows with the surface's cell count, so
/// a fixture that stopped at `65` would decide C1 on how many surface cells it
/// drew. See the header's `# Arms`.
const RESOLUTIONS: [u32; 3] = [33, 65, 129];

/// Repeats of the two timed predicate passes, interleaved. Above the five the
/// house floor asks for, because C2 names a `1.5x` threshold and M-280 measured
/// this host's governor swinging one binary `1.45x`.
const REPEATS: usize = 9;

/// Repeats of the whole extraction, for SHARE's denominator. Fewer than
/// [`REPEATS`] because an extraction is three orders of magnitude dearer than a
/// predicate pass and no clause reads the quotient.
const EXTRACT_REPEATS: usize = 3;

/// The rounding depth of `Δ` over `trilinear.rs:202-213` plus `:246`, doubled to
/// absorb the permanent's own rounding and then doubled again for headroom. See
/// the header's derivation table: `(1+u)^12 - 1 < 12.01u`, and `13u` is already
/// sound.
const FILTER_DEPTH: f64 = 26.0;

/// The safe exponent window for a Dekker–Veltkamp split, as a binary exponent.
/// Outside `2^±969` the splitting constant can overflow or the error term can
/// fall into the subnormals, and `two_product` stops being exact.
const SPLIT_SAFE_EXPONENT: i32 = 969;

/// Arithmetic operations in the naive predicate: 23 for `coefficients`
/// (`trilinear.rs:202-213`) plus 4 for `b*b - 4*a*c` (`:246`). See the header's
/// operation table.
const NAIVE_OPS: u32 = 27;

/// Operations in the filtered predicate: the naive 27, plus 35 for the permanent
/// (8 `abs` and 27 arithmetic), plus one `abs`, one compare and one scale.
const FILTER_OPS: u32 = 65;

// ─────────────────────────────────────────────────────────────────────────────
// Shewchuk's expansion arithmetic. `10.1007/pl00009321`, the DOI `real.rs:94`
// already cites. Exact for any inputs that neither overflow nor drive a split
// into the subnormals.
// ─────────────────────────────────────────────────────────────────────────────

/// Knuth's exact sum: `(x, y)` with `x + y == a + b` exactly and `x` the rounded
/// sum. Six flops and no precondition on the operands' relative magnitude, which
/// is why this file never needs `Fast_Two_Sum`.
#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    let b_virtual = x - a;
    let a_virtual = x - b_virtual;
    let b_round = b - b_virtual;
    let a_round = a - a_virtual;
    (x, a_round + b_round)
}

/// Dekker–Veltkamp split of `a` into two halves of at most 26 significand bits
/// each, at `<f64 as Real>::SPLITTER`.
#[inline]
fn split(a: f64) -> (f64, f64) {
    let c = <f64 as Real>::SPLITTER * a;
    let a_big = c - a;
    let a_hi = c - a_big;
    (a_hi, a - a_hi)
}

/// The exact product: `(x, y)` with `x + y == a * b` exactly and `x` the rounded
/// product. FMA-free, so it does not depend on the target having one.
#[inline]
fn two_product(a: f64, b: f64) -> (f64, f64) {
    let x = a * b;
    let (a_hi, a_lo) = split(a);
    let (b_hi, b_lo) = split(b);
    let err1 = x - a_hi * b_hi;
    let err2 = err1 - a_lo * b_hi;
    let err3 = err2 - a_hi * b_lo;
    (x, a_lo * b_lo - err3)
}

/// `scale_expansion_zeroelim`: the exact product of a non-overlapping increasing
/// expansion `e` and one `f64`, written into `h` in increasing order.
///
/// Returns the number of components written, at most `2 * e.len()`.
///
/// # Panics
///
/// If `e` is empty — the zero expansion is written as one zero component, never
/// as an empty slice, so an empty input means a caller lost a length.
fn scale_expansion(e: &[f64], b: f64, h: &mut [f64]) -> usize {
    assert!(!e.is_empty(), "scale_expansion: the zero expansion is [0.0]");
    let mut written = 0usize;
    let (mut q, small) = two_product(e[0], b);
    if small != 0.0 {
        h[written] = small;
        written += 1;
    }
    for &component in &e[1..] {
        let (large, small) = two_product(component, b);
        let (sum, spill) = two_sum(q, small);
        if spill != 0.0 {
            h[written] = spill;
            written += 1;
        }
        let (carried, spill) = two_sum(large, sum);
        if spill != 0.0 {
            h[written] = spill;
            written += 1;
        }
        q = carried;
    }
    if q != 0.0 || written == 0 {
        h[written] = q;
        written += 1;
    }
    written
}

/// `fast_expansion_sum_zeroelim`: the exact sum of two non-overlapping
/// increasing expansions, written into `h` in increasing order.
///
/// Returns the number of components written, at most `e.len() + f.len() + 1`.
/// Either input may be empty, which is the zero expansion arriving as a length
/// rather than as a component — the accumulator starts that way.
fn expansion_sum(e: &[f64], f: &[f64], h: &mut [f64]) -> usize {
    if e.is_empty() {
        h[..f.len()].copy_from_slice(f);
        return f.len();
    }
    if f.is_empty() {
        h[..e.len()].copy_from_slice(e);
        return e.len();
    }

    let mut ei = 0usize;
    let mut fi = 0usize;
    let mut written = 0usize;

    // Shewchuk's `(fnow > enow) == (fnow > -enow)` is `|fnow| > |enow|` without
    // forming an absolute value, and the merge always consumes the *smaller*
    // component first.
    let take_e = |ei: usize, fi: usize| -> bool {
        if ei >= e.len() {
            return false;
        }
        if fi >= f.len() {
            return true;
        }
        let (en, fnw) = (e[ei], f[fi]);
        (fnw > en) == (fnw > -en)
    };

    let mut q = if take_e(ei, fi) {
        ei += 1;
        e[ei - 1]
    } else {
        fi += 1;
        f[fi - 1]
    };

    while ei < e.len() || fi < f.len() {
        let next = if take_e(ei, fi) {
            ei += 1;
            e[ei - 1]
        } else {
            fi += 1;
            f[fi - 1]
        };
        let (carried, spill) = two_sum(q, next);
        q = carried;
        if spill != 0.0 {
            h[written] = spill;
            written += 1;
        }
    }

    if q != 0.0 || written == 0 {
        h[written] = q;
        written += 1;
    }
    written
}

/// The twelve monomials of `common::poly::repo_discriminant()`, each as four
/// factor indices and one coefficient, plus the scratch every exact evaluation
/// needs.
///
/// The factor indices are four because the form is homogeneous of degree 4;
/// [`ExactForm::new`] asserts that rather than trusting it, because the chain
/// below would be silently wrong on any other degree.
#[derive(Debug)]
struct ExactForm {
    /// `([i, j, k, l], coefficient)` per monomial, in `Poly::monomials()`'s
    /// deterministic `BTreeMap` order.
    terms: Vec<([usize; 4], f64)>,
    /// The running expansion, increasing and non-overlapping. `acc_len == 0` is
    /// the zero expansion.
    acc: Vec<f64>,
    /// How much of `acc` is live.
    acc_len: usize,
    /// The merge destination, swapped with `acc` once per monomial. Swapping two
    /// `Vec`s is three words; swapping two 192-component arrays would be 3 KB of
    /// `memcpy` per monomial per cell.
    next: Vec<f64>,
    /// The three scale stages of one monomial: at most 4, 8 and 16 components.
    stage1: [f64; 4],
    stage2: [f64; 8],
    stage3: [f64; 16],
}

impl ExactForm {
    /// Build the evaluator from the polynomial `P-127` owns.
    ///
    /// # Panics
    ///
    /// If the form is not homogeneous of degree 4, or has no terms.
    fn new() -> Self {
        let disc = poly::repo_discriminant();
        assert!(
            disc.terms() > 0,
            "repo_discriminant() is the zero polynomial, so there is no sign to be exact about"
        );
        let mut terms = Vec::with_capacity(disc.terms());
        for (exp, coefficient) in disc.monomials() {
            let degree: u32 = exp.iter().map(|e| u32::from(*e)).sum();
            assert_eq!(
                degree, 4,
                "monomial {exp:?} has total degree {degree}, and the four-factor chain here is \
                 only exact for a form that is homogeneous of degree 4"
            );
            let mut factors = [0usize; 4];
            let mut at = 0usize;
            for (i, e) in exp.iter().enumerate() {
                for _ in 0..*e {
                    factors[at] = i;
                    at += 1;
                }
            }
            terms.push((factors, coefficient as f64));
        }
        // Twelve monomials, each at most 16 components, plus the merge's carry.
        let capacity = terms.len() * 16 + 1;
        Self {
            terms,
            acc: vec![0.0; capacity],
            acc_len: 0,
            next: vec![0.0; capacity],
            stage1: [0.0; 4],
            stage2: [0.0; 8],
            stage3: [0.0; 16],
        }
    }

    /// How many monomials the form has. Twelve, and `P-127` proved it.
    fn term_count(&self) -> usize {
        self.terms.len()
    }

    /// The exact sign of `Δ` at these eight corner values: `-1`, `0` or `+1`.
    ///
    /// The expansion is non-overlapping and increasing, so its sign is its
    /// largest component's — no component can be outvoted by the tail below it.
    fn sign(&mut self, f: &[f64; 8]) -> i8 {
        self.acc_len = 0;
        for t in 0..self.terms.len() {
            let (idx, coefficient) = self.terms[t];
            let (large, small) = two_product(f[idx[0]], f[idx[1]]);
            let seed = [small, large];
            let l1 = scale_expansion(&seed, f[idx[2]], &mut self.stage1);
            let l2 = scale_expansion(&self.stage1[..l1], f[idx[3]], &mut self.stage2);
            let l3 = scale_expansion(&self.stage2[..l2], coefficient, &mut self.stage3);
            let merged = expansion_sum(&self.acc[..self.acc_len], &self.stage3[..l3], &mut self.next);
            std::mem::swap(&mut self.acc, &mut self.next);
            self.acc_len = merged;
        }
        if self.acc_len == 0 {
            return 0;
        }
        let top = self.acc[self.acc_len - 1];
        if top > 0.0 {
            1
        } else if top < 0.0 {
            -1
        } else {
            0
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// The float predicate and its filter, at both precisions.
// ─────────────────────────────────────────────────────────────────────────────

/// `Δ` at `f64`, through the crate's own chain: `BodySaddles::coefficients` is
/// called rather than retyped, and `:246`'s expression is applied to what it
/// returns.
#[inline]
fn discriminant_f64(corner: &[f64; 8]) -> f64 {
    let [a, b, c] = BodySaddles::<f64>::coefficients(corner);
    b * b - 4.0 * a * c
}

/// `Δ` at `f32`, the same way.
#[inline]
fn discriminant_f32(corner: &[f32; 8]) -> f32 {
    let [a, b, c] = BodySaddles::<f32>::coefficients(corner);
    b * b - 4.0 * a * c
}

/// The permanent of `Δ`'s expression tree at `f64`: every operand replaced by
/// its absolute value and every subtraction by an addition, following
/// `trilinear.rs:202-213` and `:246` node for node.
#[inline]
fn permanent_f64(corner: &[f64; 8]) -> f64 {
    let f: [f64; 8] = std::array::from_fn(|i| corner[i].abs());
    let twist_lo = f[0] + f[3] + f[1] + f[2];
    let twist_hi = f[4] + f[7] + f[5] + f[6];
    let du_lo = f[1] + f[0];
    let du_hi = f[5] + f[4];
    let dv_lo = f[2] + f[0];
    let dv_hi = f[6] + f[4];
    let pa = du_hi * twist_lo + du_lo * twist_hi;
    let pb = (f[4] * twist_lo + f[0] * twist_hi) + (du_hi * dv_lo + du_lo * dv_hi);
    let pc = f[2] * f[4] + f[0] * f[6];
    pb * pb + 4.0 * pa * pc
}

/// The same permanent at `f32`.
#[inline]
fn permanent_f32(corner: &[f32; 8]) -> f32 {
    let f: [f32; 8] = std::array::from_fn(|i| corner[i].abs());
    let twist_lo = f[0] + f[3] + f[1] + f[2];
    let twist_hi = f[4] + f[7] + f[5] + f[6];
    let du_lo = f[1] + f[0];
    let du_hi = f[5] + f[4];
    let dv_lo = f[2] + f[0];
    let dv_hi = f[6] + f[4];
    let pa = du_hi * twist_lo + du_lo * twist_hi;
    let pb = (f[4] * twist_lo + f[0] * twist_hi) + (du_hi * dv_lo + du_lo * dv_hi);
    let pc = f[2] * f[4] + f[0] * f[6];
    pb * pb + 4.0 * pa * pc
}

/// The `f64` filter's absolute threshold: `26u·P̂`.
#[inline]
fn bound_f64(corner: &[f64; 8]) -> f64 {
    FILTER_DEPTH * <f64 as Real>::UNIT_ROUNDOFF * permanent_f64(corner)
}

/// The `f32` filter's absolute threshold: `26u·P̂` with `f32`'s unit roundoff.
#[inline]
fn bound_f32(corner: &[f32; 8]) -> f32 {
    (FILTER_DEPTH as f32) * <f32 as Real>::UNIT_ROUNDOFF * permanent_f32(corner)
}

/// The sign of a float, with `±0.0` reading `0`. NaN reads `0` and is counted
/// separately by `nonfinite_deltas`, so it can never be silently scored as an
/// agreement.
#[inline]
fn sign_of(x: f64) -> i8 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// The fixture.
// ─────────────────────────────────────────────────────────────────────────────

/// A grid of corner values served as an [`Sdf`], trilinearly interpolated.
///
/// Exact at every grid point in both precisions — see the header — so the two
/// extraction arms read one shared corner set, and finite off-grid so
/// `unit_gradient` has something to normalise.
#[derive(Debug)]
struct Tabulated<'a, R: Real> {
    /// Row-major samples, `x` fastest, as `sample_grid` writes them.
    values: &'a [R],
    /// Samples per axis, not cells.
    size: [u32; 3],
    /// World position of sample `[0, 0, 0]`.
    origin: [R; 3],
    /// Spacing between adjacent samples.
    cell_size: R,
}

impl<R: Real> Tabulated<'_, R> {
    /// The stored value at an integer sample index.
    #[inline]
    fn at(&self, x: usize, y: usize, z: usize) -> R {
        let sx = self.size[0] as usize;
        let sy = self.size[1] as usize;
        self.values[x + y * sx + z * sx * sy]
    }
}

impl<R: Real> Sdf for Tabulated<'_, R> {
    type Scalar = R;

    fn sample(&self, p: [R; 3]) -> R {
        let mut base = [0usize; 3];
        let mut t = [R::ZERO; 3];
        for k in 0..3 {
            let g = (p[k] - self.origin[k]) / self.cell_size;
            // The last cell owns its far face, so a coordinate exactly on the
            // upper boundary interpolates with `t == 1` rather than indexing out.
            let top = f64::from(self.size[k] - 2);
            let i = g.floor().as_f64().clamp(0.0, top);
            base[k] = i as usize;
            t[k] = g - R::from_f64(i);
        }
        let (x, y, z) = (base[0], base[1], base[2]);
        let lerp = |a: R, b: R, s: R| a * (R::ONE - s) + b * s;
        let c00 = lerp(self.at(x, y, z), self.at(x + 1, y, z), t[0]);
        let c10 = lerp(self.at(x, y + 1, z), self.at(x + 1, y + 1, z), t[0]);
        let c01 = lerp(self.at(x, y, z + 1), self.at(x + 1, y, z + 1), t[0]);
        let c11 = lerp(self.at(x, y + 1, z + 1), self.at(x + 1, y + 1, z + 1), t[0]);
        lerp(lerp(c00, c10, t[1]), lerp(c01, c11, t[1]), t[2])
    }
}

/// Everything one `(field, resolution)` cell of the fixture produced, before the
/// global verdicts are known.
#[derive(Debug)]
struct Measured {
    /// The `ReferenceField` name.
    field: &'static str,
    /// Samples per axis.
    samples: u32,
    /// Cells with a sign change, which is where `Δ` means anything.
    cells: u64,
    /// Cells in the grid, surface or not.
    grid_cells: u64,
    /// Surface cells whose case has an ambiguous face — the population the
    /// trilinear path visits.
    ambiguous_cells: u64,
    /// Surface cells where the `f32` `Δ` sign differs from the exact sign.
    disagreements_f32: u64,
    /// The same at `f64`.
    disagreements_f64: u64,
    /// The subset of `disagreements_f32` on an ambiguous-face cell.
    disagreements_f32_ambiguous: u64,
    /// Surface cells with `|Δ_f32| <= 26u·P̂_f32`: the registered vacuity count.
    below_f32_bound: u64,
    /// Surface cells with `|Δ_f64| <= 26u·P̂_f64`.
    below_f64_bound: u64,
    /// Surface cells where the `f32` and `f64` `Δ` signs differ from each other.
    disc_sign_changed: u64,
    /// The subset of `disc_sign_changed` on an ambiguous-face cell — the only
    /// cells where a corrected sign can reach the mesh, and therefore the column
    /// C3's falsifier is really asking for.
    disc_sign_changed_ambiguous: u64,
    /// Surface cells where `BodySaddles::inside_mask` differs between precisions.
    inside_mask_changed: u64,
    /// Surface cells where `has_inner_hexagon` differs between precisions.
    hexagon_changed: u64,
    /// Ambiguous cells where the face decider's `joined_mask` differs.
    joined_mask_changed: u64,
    /// Cells where the two precisions disagreed on the case index. Structurally
    /// zero, and asserted.
    case_index_changed: u64,
    /// Non-finite `Δ` evaluations at either precision.
    nonfinite: u64,
    /// Exact signs of each kind over the surface cells.
    exact_positive: u64,
    /// Exact signs of each kind over the surface cells.
    exact_negative: u64,
    /// Exact signs of each kind over the surface cells.
    exact_zero: u64,
    /// Cells where a filter certified a sign the exact expansion contradicts,
    /// summed over both precisions. Asserted zero.
    filter_certified_wrong: u64,
    /// Largest `|corner value|` over the whole grid.
    max_abs_corner: f64,
    /// Smallest non-zero `|corner value|` over the whole grid.
    min_abs_nonzero_corner: f64,
    /// Triangles the `f32` arm emitted.
    triangles_f32: u64,
    /// Triangles the `f64` arm emitted.
    triangles_f64: u64,
    /// Non-manifold edges the `f32` arm produced.
    nonmanifold_f32: u64,
    /// Non-manifold edges the `f64` arm produced.
    nonmanifold_f64: u64,
    /// Per-scalar timings and rates, indexed `0 = f32`, `1 = f64`.
    arms: [Arm; 2],
}

/// One precision's timed columns for one `(field, resolution)`.
#[derive(Debug, Default, Clone, Copy)]
struct Arm {
    /// Median of `REPEATS` naive-float passes over the surface cells.
    float_ms: f64,
    /// Median of `REPEATS` filtered-exact passes over the same cells.
    exact_ms: f64,
    /// Min and max of the filtered-exact passes.
    exact_min_ms: f64,
    /// Min and max of the filtered-exact passes.
    exact_max_ms: f64,
    /// Min and max of the naive-float passes.
    float_min_ms: f64,
    /// Min and max of the naive-float passes.
    float_max_ms: f64,
    /// Cells the filter could not certify, over `cells`.
    fallback_rate: f64,
    /// Cells the filter could not certify.
    fallbacks: u64,
    /// Median of `REPEATS` `BodySaddles::of` passes over the ambiguous cells.
    saddle_ms: f64,
    /// `(max - min) / median` of those passes.
    saddle_scatter: f64,
    /// Median of `EXTRACT_REPEATS` whole extractions.
    extract_ms: f64,
    /// Min and max of those extractions.
    extract_min_ms: f64,
    /// Min and max of those extractions.
    extract_max_ms: f64,
}

/// The median, min and max of a set of timings, sorted with `total_cmp`.
fn spread(mut samples: Vec<f64>) -> (f64, f64, f64) {
    samples.sort_unstable_by(f64::total_cmp);
    let n = samples.len();
    (samples[n / 2], samples[0], samples[n - 1])
}

/// Sample one reference field over its own grid, then round every value through
/// `f32` and back.
///
/// The rounding is the fixture's defining decision: it makes the corner set
/// bit-identical in both precisions, so a sign disagreement is about arithmetic
/// and not about representation.
fn shared_corner_set<F>(field: &F, shape: &RuntimeShape3, origin: [f64; 3], cell_size: f64) -> Vec<f64>
where
    F: Sdf<Scalar = f64>,
{
    let size = shape.size();
    let mut values = Vec::with_capacity(shape.element_count());
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let p = [
                    origin[0] + cell_size * f64::from(x),
                    origin[1] + cell_size * f64::from(y),
                    origin[2] + cell_size * f64::from(z),
                ];
                values.push(f64::from(field.sample(p) as f32));
            }
        }
    }
    values
}

/// The eight corner sample indices of the cell whose minimum corner is `base`,
/// in `cube.rs`'s numbering — `corner_offset` is private, and it is
/// `[(c & 1), (c >> 1) & 1, (c >> 2) & 1]` (`cube.rs:149-155`).
#[inline]
fn corner_indices(shape: &RuntimeShape3, base: [u32; 3]) -> [u32; 8] {
    std::array::from_fn(|c| {
        let c = c as u32;
        shape.linearize([
            base[0] + (c & 1),
            base[1] + ((c >> 1) & 1),
            base[2] + ((c >> 2) & 1),
        ])
    })
}

/// Extract with the golden roster's `marching_cubes+trilinear` configuration —
/// `AsymptoticDecider` plus `Trilinear`, which is the only way to reach the
/// body-saddle path (`mod.rs:278-292`) — and report triangles, non-manifold
/// edges and the median extraction time.
fn extract_arm<R: Real>(
    values: &[R],
    shape: &RuntimeShape3,
    origin: [R; 3],
    cell_size: R,
) -> (u64, u64, f64, f64, f64) {
    let field = Tabulated {
        values,
        size: shape.size(),
        origin,
        cell_size,
    };
    let mut mc = MarchingCubes::<R>::new();
    mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
    mc.set_interior_ambiguity(InteriorAmbiguity::Trilinear);
    let mut mesh = MeshBuffer::<R>::new();

    let mut timings = Vec::with_capacity(EXTRACT_REPEATS);
    for _ in 0..EXTRACT_REPEATS {
        mesh.reset();
        let started = Instant::now();
        mc.extract(&field, shape, origin, cell_size, &mut mesh)
            .expect("the tabulated grid has at least two samples per axis");
        timings.push(started.elapsed().as_secs_f64() * 1e3);
    }
    let (median, min, max) = spread(timings);

    let cfg = ValidateConfig::from_cell_size(cell_size.as_f64())
        .expect("a reference field's cell size is finite and positive");
    let report = validate(&mesh, &cfg);
    (
        mesh.triangle_count() as u64,
        report.non_manifold_edges,
        median,
        min,
        max,
    )
}

/// Everything one `(field, resolution)` contributes, measured.
fn measure(
    field: &'static str,
    samples: u32,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    cell_size: f64,
    values64: &[f64],
    exact: &mut ExactForm,
) -> Measured {
    let size = shape.size();
    let values32: Vec<f32> = values64.iter().map(|v| *v as f32).collect();

    let mut max_abs_corner = 0.0f64;
    let mut min_abs_nonzero_corner = f64::INFINITY;
    for v in values64 {
        let a = v.abs();
        if a > max_abs_corner {
            max_abs_corner = a;
        }
        if a > 0.0 && a < min_abs_nonzero_corner {
            min_abs_nonzero_corner = a;
        }
    }

    // ── the surface-cell population, and its per-cell audit ──────────────────
    let mut corners64: Vec<[f64; 8]> = Vec::new();
    let mut corners32: Vec<[f32; 8]> = Vec::new();
    let mut ambiguous_of: Vec<u8> = Vec::new();

    let mut m = Measured {
        field,
        samples,
        cells: 0,
        grid_cells: u64::from(size[0] - 1) * u64::from(size[1] - 1) * u64::from(size[2] - 1),
        ambiguous_cells: 0,
        disagreements_f32: 0,
        disagreements_f64: 0,
        disagreements_f32_ambiguous: 0,
        below_f32_bound: 0,
        below_f64_bound: 0,
        disc_sign_changed: 0,
        disc_sign_changed_ambiguous: 0,
        inside_mask_changed: 0,
        hexagon_changed: 0,
        joined_mask_changed: 0,
        case_index_changed: 0,
        nonfinite: 0,
        exact_positive: 0,
        exact_negative: 0,
        exact_zero: 0,
        filter_certified_wrong: 0,
        max_abs_corner,
        min_abs_nonzero_corner,
        triangles_f32: 0,
        triangles_f64: 0,
        nonmanifold_f32: 0,
        nonmanifold_f64: 0,
        arms: [Arm::default(); 2],
    };

    for z in 0..size[2] - 1 {
        for y in 0..size[1] - 1 {
            for x in 0..size[0] - 1 {
                let idx = corner_indices(shape, [x, y, z]);
                let c64: [f64; 8] = std::array::from_fn(|c| values64[idx[c] as usize]);
                let c32: [f32; 8] = std::array::from_fn(|c| values32[idx[c] as usize]);

                let mut case64 = 0u8;
                let mut case32 = 0u8;
                for c in 0..8 {
                    if is_inside(c64[c]) {
                        case64 |= 1 << c;
                    }
                    if is_inside(c32[c]) {
                        case32 |= 1 << c;
                    }
                }
                if case64 != case32 {
                    m.case_index_changed += 1;
                }
                if case64 == 0 || case64 == u8::MAX {
                    continue;
                }
                m.cells += 1;

                let ambiguous = AMBIGUOUS_FACES[case64 as usize];
                if ambiguous != 0 {
                    m.ambiguous_cells += 1;
                    if joined_mask(&c64, ambiguous) != joined_mask(&c32, ambiguous) {
                        m.joined_mask_changed += 1;
                    }
                }

                let d64 = discriminant_f64(&c64);
                let d32 = discriminant_f32(&c32);
                if !d64.is_finite() || !d32.is_finite() {
                    m.nonfinite += 1;
                }

                let exact_sign = exact.sign(&c64);
                match exact_sign {
                    1 => m.exact_positive += 1,
                    -1 => m.exact_negative += 1,
                    _ => m.exact_zero += 1,
                }

                let s64 = sign_of(d64);
                let s32 = sign_of(f64::from(d32));
                if s64 != exact_sign {
                    m.disagreements_f64 += 1;
                }
                if s32 != exact_sign {
                    m.disagreements_f32 += 1;
                    if ambiguous != 0 {
                        m.disagreements_f32_ambiguous += 1;
                    }
                }
                if s32 != s64 {
                    m.disc_sign_changed += 1;
                    if ambiguous != 0 {
                        m.disc_sign_changed_ambiguous += 1;
                    }
                }

                // The filters, and their soundness. A certified sign that the
                // exact expansion contradicts means the derived bound is wrong.
                let b64 = bound_f64(&c64);
                if d64.abs() <= b64 {
                    m.below_f64_bound += 1;
                } else if s64 != exact_sign {
                    m.filter_certified_wrong += 1;
                }
                let b32 = bound_f32(&c32);
                if d32.abs() <= b32 {
                    m.below_f32_bound += 1;
                } else if s32 != exact_sign {
                    m.filter_certified_wrong += 1;
                }

                let saddles64 = BodySaddles::<f64>::of(&c64);
                let saddles32 = BodySaddles::<f32>::of(&c32);
                if saddles64.inside_mask() != saddles32.inside_mask() {
                    m.inside_mask_changed += 1;
                }
                if saddles64.has_inner_hexagon() != saddles32.has_inner_hexagon() {
                    m.hexagon_changed += 1;
                }

                corners64.push(c64);
                corners32.push(c32);
                ambiguous_of.push(ambiguous);
            }
        }
    }

    // ── C2, timed, interleaved, medians ──────────────────────────────────────
    let mut sink = 0i64;

    // Warm up once, so the first timed pass is not paying for cold caches.
    for c in &corners64 {
        sink += i64::from(sign_of(discriminant_f64(c)));
    }

    let mut float32 = Vec::with_capacity(REPEATS);
    let mut exact32 = Vec::with_capacity(REPEATS);
    let mut float64 = Vec::with_capacity(REPEATS);
    let mut exact64 = Vec::with_capacity(REPEATS);
    let mut saddle32 = Vec::with_capacity(REPEATS);
    let mut saddle64 = Vec::with_capacity(REPEATS);
    let mut fallbacks32 = 0u64;
    let mut fallbacks64 = 0u64;

    for _ in 0..REPEATS {
        // f32 naive
        let started = Instant::now();
        for c in &corners32 {
            sink += i64::from(sign_of(f64::from(discriminant_f32(c))));
        }
        float32.push(started.elapsed().as_secs_f64() * 1e3);

        // f32 filtered exact
        let mut fell_back = 0u64;
        let started = Instant::now();
        for c in &corners32 {
            let d = discriminant_f32(c);
            let s = if d.abs() > bound_f32(c) {
                sign_of(f64::from(d))
            } else {
                fell_back += 1;
                let widened: [f64; 8] = std::array::from_fn(|i| f64::from(c[i]));
                exact.sign(&widened)
            };
            sink += i64::from(s);
        }
        exact32.push(started.elapsed().as_secs_f64() * 1e3);
        fallbacks32 = fell_back;

        // f64 naive
        let started = Instant::now();
        for c in &corners64 {
            sink += i64::from(sign_of(discriminant_f64(c)));
        }
        float64.push(started.elapsed().as_secs_f64() * 1e3);

        // f64 filtered exact
        let mut fell_back = 0u64;
        let started = Instant::now();
        for c in &corners64 {
            let d = discriminant_f64(c);
            let s = if d.abs() > bound_f64(c) {
                sign_of(d)
            } else {
                fell_back += 1;
                exact.sign(c)
            };
            sink += i64::from(s);
        }
        exact64.push(started.elapsed().as_secs_f64() * 1e3);
        fallbacks64 = fell_back;

        // SHARE: the stage C2 moves, at each precision, over the cells that
        // actually take the trilinear path.
        let started = Instant::now();
        for (c, a) in corners32.iter().zip(ambiguous_of.iter()) {
            if *a != 0 {
                sink += i64::from(BodySaddles::<f32>::of(c).inside_mask());
            }
        }
        saddle32.push(started.elapsed().as_secs_f64() * 1e3);

        let started = Instant::now();
        for (c, a) in corners64.iter().zip(ambiguous_of.iter()) {
            if *a != 0 {
                sink += i64::from(BodySaddles::<f64>::of(c).inside_mask());
            }
        }
        saddle64.push(started.elapsed().as_secs_f64() * 1e3);
    }

    // The accumulator is why none of the six loops can be elided. Folding it
    // into a column keeps it observable rather than relying on `black_box`.
    assert!(
        sink != i64::MIN,
        "the sign accumulator saturated, which means the timed loops measured something else"
    );

    // ── C3: the two extraction arms over the shared corner set ───────────────
    let origin32: [f32; 3] = std::array::from_fn(|k| origin[k] as f32);
    let (tri32, nm32, ex32, ex32_min, ex32_max) =
        extract_arm(&values32, shape, origin32, cell_size as f32);
    let (tri64, nm64, ex64, ex64_min, ex64_max) =
        extract_arm(values64, shape, origin, cell_size);

    m.triangles_f32 = tri32;
    m.triangles_f64 = tri64;
    m.nonmanifold_f32 = nm32;
    m.nonmanifold_f64 = nm64;

    let cells = m.cells.max(1) as f64;
    let (f32_float, f32_float_min, f32_float_max) = spread(float32);
    let (f32_exact, f32_exact_min, f32_exact_max) = spread(exact32);
    let (f32_saddle, f32_saddle_min, f32_saddle_max) = spread(saddle32);
    m.arms[0] = Arm {
        float_ms: f32_float,
        exact_ms: f32_exact,
        exact_min_ms: f32_exact_min,
        exact_max_ms: f32_exact_max,
        float_min_ms: f32_float_min,
        float_max_ms: f32_float_max,
        fallback_rate: fallbacks32 as f64 / cells,
        fallbacks: fallbacks32,
        saddle_ms: f32_saddle,
        saddle_scatter: scatter(f32_saddle, f32_saddle_min, f32_saddle_max),
        extract_ms: ex32,
        extract_min_ms: ex32_min,
        extract_max_ms: ex32_max,
    };

    let (f64_float, f64_float_min, f64_float_max) = spread(float64);
    let (f64_exact, f64_exact_min, f64_exact_max) = spread(exact64);
    let (f64_saddle, f64_saddle_min, f64_saddle_max) = spread(saddle64);
    m.arms[1] = Arm {
        float_ms: f64_float,
        exact_ms: f64_exact,
        exact_min_ms: f64_exact_min,
        exact_max_ms: f64_exact_max,
        float_min_ms: f64_float_min,
        float_max_ms: f64_float_max,
        fallback_rate: fallbacks64 as f64 / cells,
        fallbacks: fallbacks64,
        saddle_ms: f64_saddle,
        saddle_scatter: scatter(f64_saddle, f64_saddle_min, f64_saddle_max),
        extract_ms: ex64,
        extract_min_ms: ex64_min,
        extract_max_ms: ex64_max,
    };

    m
}

/// `(max - min) / median`, and zero when the median is zero — a stage that took
/// no measurable time has no scatter to report.
fn scatter(median: f64, min: f64, max: f64) -> f64 {
    if median > 0.0 { (max - min) / median } else { 0.0 }
}

/// The binary exponent of a finite non-zero `f64`, for the split-safety assert.
fn exponent_of(x: f64) -> i32 {
    if x == 0.0 || !x.is_finite() {
        return 0;
    }
    let bits = x.to_bits();
    let raw = ((bits >> 52) & 0x7ff) as i32;
    if raw == 0 {
        // Subnormal: the exponent is the fixed minimum, which is already far
        // outside the safe window and will trip the assert, as it should.
        -1074
    } else {
        raw - 1023
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-133");

    common::experiment::run(prereg, |run| {
        let mut exact = ExactForm::new();
        let terms = exact.term_count();

        let mut measured: Vec<Measured> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            for &samples in &RESOLUTIONS {
                let (shape, origin, cell_size) = common::grid::<f64, _>(&field, samples);
                let values = shared_corner_set(&field, &shape, origin, cell_size);
                measured.push(measure(
                    name, samples, &shape, origin, cell_size, &values, &mut exact,
                ));
            }
        });

        // ── vacuity and soundness controls, all before the first row ─────────
        assert_eq!(
            terms, 12,
            "VOID: repo_discriminant() has {terms} monomials, not the twelve P-127 measured, so \
             the exact oracle is not the polynomial this row is about"
        );

        for m in &measured {
            assert!(
                m.cells > 0,
                "VOID: {} at {}^3 has no surface cell, so disagreement_rate is a zero that could \
                 not have been non-zero (M-44)",
                m.field,
                m.samples
            );
            assert_eq!(
                m.case_index_changed, 0,
                "VOID: {} at {}^3 has {} cells whose case index differs between f32 and f64, so \
                 the corner set is not shared and C3's triangle delta would be a case-index \
                 difference rather than a body-saddle one",
                m.field, m.samples, m.case_index_changed
            );
            assert_eq!(
                m.nonfinite, 0,
                "VOID: {} at {}^3 produced {} non-finite Delta evaluations, so a sign comparison \
                 there is measuring overflow rather than cancellation",
                m.field, m.samples, m.nonfinite
            );
            assert_eq!(
                m.filter_certified_wrong, 0,
                "VOID: {} at {}^3 has {} cells where the static filter certified a sign the exact \
                 expansion contradicts, so the 26u bound derived in this file's header is unsound \
                 and every filtered_fallback_rate below it is meaningless",
                m.field, m.samples, m.filter_certified_wrong
            );
            let hi = exponent_of(m.max_abs_corner);
            let lo = exponent_of(m.min_abs_nonzero_corner);
            assert!(
                hi <= SPLIT_SAFE_EXPONENT && lo >= -SPLIT_SAFE_EXPONENT,
                "VOID: {} at {}^3 spans binary exponents {lo}..{hi}, outside the +-{} window \
                 where a Dekker-Veltkamp split is exact, so the expansion arithmetic is not the \
                 oracle this row claims",
                m.field,
                m.samples,
                SPLIT_SAFE_EXPONENT
            );
        }

        let below_f32_total: u64 = measured.iter().map(|m| m.below_f32_bound).sum();
        assert!(
            below_f32_total > 0,
            "VOID: no cell on any field has |Delta| below the f32 error bound, so the fixture \
             excludes the near-degenerate stratum entirely and C1 is measuring a problem it \
             cannot see -- this is the registration's own vacuity control"
        );

        // A cell whose exact `Delta` is zero is below every bound at every
        // precision and both arms agree on it, so the registered count above can
        // be satisfied by a population in which `f32` cannot possibly be wrong.
        // This is the part of the stratum that licenses C1.
        let beyond_zero_total: u64 = measured
            .iter()
            .map(|m| m.below_f32_bound.saturating_sub(m.exact_zero))
            .sum();
        assert!(
            beyond_zero_total > 0,
            "VOID: every cell below the f32 error bound has an exactly zero Delta \
             ({below_f32_total} of them), and both precisions get an exact zero right, so the \
             fixture holds no cell where f32 could disagree with the exact sign and C1's count \
             is structurally zero rather than measured"
        );

        let ambiguous_total: u64 = measured.iter().map(|m| m.ambiguous_cells).sum();
        assert!(
            ambiguous_total > 0,
            "VOID: no cell on any field has an ambiguous face, so InteriorAmbiguity::Trilinear is \
             unreachable (mod.rs:278-292), the body-saddle stage never runs, and C3's zero would \
             be structural rather than measured"
        );

        let exact_positive: u64 = measured.iter().map(|m| m.exact_positive).sum();
        let exact_negative: u64 = measured.iter().map(|m| m.exact_negative).sum();
        assert!(
            exact_positive > 0 && exact_negative > 0,
            "VOID: the exact oracle reads {exact_positive} positive and {exact_negative} negative \
             over the whole fixture, so it is effectively constant and agreeing with it says \
             nothing"
        );

        // ── the global verdicts ──────────────────────────────────────────────
        let mut fields_disagreeing: Vec<&'static str> = measured
            .iter()
            .filter(|m| m.disagreements_f32 > 0)
            .map(|m| m.field)
            .collect();
        fields_disagreeing.sort_unstable();
        fields_disagreeing.dedup();
        let c1 = fields_disagreeing.len() >= 2;

        let float_total: f64 = measured
            .iter()
            .map(|m| m.arms[0].float_ms + m.arms[1].float_ms)
            .sum();
        let exact_total: f64 = measured
            .iter()
            .map(|m| m.arms[0].exact_ms + m.arms[1].exact_ms)
            .sum();
        let aggregate = exact_total / float_total;
        let aggregate32: f64 = measured.iter().map(|m| m.arms[0].exact_ms).sum::<f64>()
            / measured.iter().map(|m| m.arms[0].float_ms).sum::<f64>();
        let aggregate64: f64 = measured.iter().map(|m| m.arms[1].exact_ms).sum::<f64>()
            / measured.iter().map(|m| m.arms[1].float_ms).sum::<f64>();
        let c2 = aggregate < 1.5;

        let c3 = measured.iter().any(|m| {
            m.triangles_f32 != m.triangles_f64 || m.nonmanifold_f32 != m.nonmanifold_f64
        });
        // C3's mechanism, kept apart from C3's verdict: the mesh may move for a
        // rounding in the root *values* while every discriminant sign that could
        // reach the output agrees. That is the case C3's falsifier anticipates and
        // it is not the same claim as the clause.
        let c3_by_sign = measured.iter().any(|m| {
            (m.triangles_f32 != m.triangles_f64 || m.nonmanifold_f32 != m.nonmanifold_f64)
                && m.disc_sign_changed_ambiguous > 0
        });

        // ── 48 rows: eight fields, three resolutions, two scalars ────────────
        for m in &measured {
            let cells = m.cells as f64;
            let triangles_changed = m.triangles_f64.abs_diff(m.triangles_f32);
            let nonmanifold_delta = m.nonmanifold_f64 as i64 - m.nonmanifold_f32 as i64;
            let c3_row = triangles_changed > 0 || nonmanifold_delta != 0;

            for (which, arm) in m.arms.iter().enumerate() {
                let scalar = if which == 0 { "f32" } else { "f64" };
                let own_disagreements = if which == 0 {
                    m.disagreements_f32
                } else {
                    m.disagreements_f64
                };
                let below_bound = if which == 0 {
                    m.below_f32_bound
                } else {
                    m.below_f64_bound
                };
                let ratio = arm.exact_ms / arm.float_ms;
                let c1_row = own_disagreements > 0;

                run.record(&[
                    // ── the registered columns, in registration order ────────
                    ("field", m.field.to_string()),
                    ("resolution", format!("{0}x{0}x{0}", m.samples)),
                    ("scalar", scalar.to_string()),
                    ("cells", m.cells.to_string()),
                    (
                        "sign_disagreements_f32",
                        m.disagreements_f32.to_string(),
                    ),
                    (
                        "sign_disagreements_f64",
                        m.disagreements_f64.to_string(),
                    ),
                    (
                        "disagreement_rate",
                        format!("{:.9}", own_disagreements as f64 / cells),
                    ),
                    ("exact_ms", format!("{:.6}", arm.exact_ms)),
                    ("float_ms", format!("{:.6}", arm.float_ms)),
                    ("overhead_ratio", format!("{ratio:.6}")),
                    (
                        "filtered_fallback_rate",
                        format!("{:.9}", arm.fallback_rate),
                    ),
                    ("triangles_changed", triangles_changed.to_string()),
                    ("nonmanifold_delta", nonmanifold_delta.to_string()),
                    ("c1_holds", c1.to_string()),
                    ("c2_holds", c2.to_string()),
                    ("c3_holds", c3.to_string()),
                    // ── extras (M-273) ──────────────────────────────────────
                    //
                    // Per-row clause answers, so three identical global
                    // verdicts are not mistaken for three measurements.
                    ("c1_row", c1_row.to_string()),
                    ("c2_row", (ratio < 1.5).to_string()),
                    ("c3_row", c3_row.to_string()),
                    (
                        "fields_disagreeing_f32",
                        fields_disagreeing.len().to_string(),
                    ),
                    (
                        "fields_disagreeing_names",
                        fields_disagreeing.join("|"),
                    ),
                    ("c2_aggregate_ratio", format!("{aggregate:.6}")),
                    ("c2_aggregate_f32", format!("{aggregate32:.6}")),
                    ("c2_aggregate_f64", format!("{aggregate64:.6}")),
                    // The fixture.
                    ("samples_per_axis", m.samples.to_string()),
                    ("grid_cells", m.grid_cells.to_string()),
                    ("surface_cell_rate", format!("{:.9}", cells / m.grid_cells as f64)),
                    ("ambiguous_cells", m.ambiguous_cells.to_string()),
                    (
                        "ambiguous_cell_rate",
                        format!("{:.9}", m.ambiguous_cells as f64 / cells),
                    ),
                    ("max_abs_corner", format!("{:.6e}", m.max_abs_corner)),
                    (
                        "min_abs_nonzero_corner",
                        format!("{:.6e}", m.min_abs_nonzero_corner),
                    ),
                    // The registered vacuity control, and the f64 twin that
                    // shows how much of the stratum survives a wider mantissa.
                    (
                        "below_f32_error_bound",
                        m.below_f32_bound.to_string(),
                    ),
                    (
                        "below_f64_error_bound",
                        m.below_f64_bound.to_string(),
                    ),
                    ("below_bound_this_scalar", below_bound.to_string()),
                    ("filter_fallbacks", arm.fallbacks.to_string()),
                    (
                        "filter_certified_wrong",
                        m.filter_certified_wrong.to_string(),
                    ),
                    ("filter_depth_u", format!("{FILTER_DEPTH:.1}")),
                    (
                        "filter_relative_bound",
                        format!(
                            "{:.6e}",
                            if which == 0 {
                                FILTER_DEPTH * f64::from(<f32 as Real>::UNIT_ROUNDOFF)
                            } else {
                                FILTER_DEPTH * <f64 as Real>::UNIT_ROUNDOFF
                            }
                        ),
                    ),
                    // C1's decomposition: where the sign errors fall.
                    (
                        "sign_disagreements_f32_ambiguous",
                        m.disagreements_f32_ambiguous.to_string(),
                    ),
                    (
                        "disc_sign_cells_changed",
                        m.disc_sign_changed.to_string(),
                    ),
                    ("exact_positive", m.exact_positive.to_string()),
                    ("exact_negative", m.exact_negative.to_string()),
                    ("exact_zero", m.exact_zero.to_string()),
                    ("exact_terms", terms.to_string()),
                    ("case_index_disagreements", m.case_index_changed.to_string()),
                    ("nonfinite_deltas", m.nonfinite.to_string()),
                    // C3's attribution: which precision-dependent decision
                    // inside `AsymptoticDecider + Trilinear` moved.
                    ("triangles_f32", m.triangles_f32.to_string()),
                    ("triangles_f64", m.triangles_f64.to_string()),
                    ("nonmanifold_edges_f32", m.nonmanifold_f32.to_string()),
                    ("nonmanifold_edges_f64", m.nonmanifold_f64.to_string()),
                    (
                        "inside_mask_cells_changed",
                        m.inside_mask_changed.to_string(),
                    ),
                    ("hexagon_cells_changed", m.hexagon_changed.to_string()),
                    (
                        "joined_mask_cells_changed",
                        m.joined_mask_changed.to_string(),
                    ),
                    (
                        "f64_is_the_corrected_arm",
                        (m.disagreements_f64 == 0).to_string(),
                    ),
                    // C2's scatter, and SHARE.
                    ("repeats", REPEATS.to_string()),
                    ("exact_min_ms", format!("{:.6}", arm.exact_min_ms)),
                    ("exact_max_ms", format!("{:.6}", arm.exact_max_ms)),
                    ("float_min_ms", format!("{:.6}", arm.float_min_ms)),
                    ("float_max_ms", format!("{:.6}", arm.float_max_ms)),
                    (
                        "overhead_ratio_scatter",
                        format!(
                            "{:.6}",
                            scatter(arm.exact_ms, arm.exact_min_ms, arm.exact_max_ms)
                                + scatter(arm.float_ms, arm.float_min_ms, arm.float_max_ms)
                        ),
                    ),
                    ("saddle_stage_ms", format!("{:.6}", arm.saddle_ms)),
                    (
                        "saddle_stage_scatter",
                        format!("{:.6}", arm.saddle_scatter),
                    ),
                    ("extract_ms", format!("{:.6}", arm.extract_ms)),
                    ("extract_min_ms", format!("{:.6}", arm.extract_min_ms)),
                    ("extract_max_ms", format!("{:.6}", arm.extract_max_ms)),
                    (
                        "saddle_share",
                        format!(
                            "{:.6}",
                            if arm.extract_ms > 0.0 {
                                arm.saddle_ms / arm.extract_ms
                            } else {
                                0.0
                            }
                        ),
                    ),
                    ("extract_repeats", EXTRACT_REPEATS.to_string()),
                ]);
            }
        }
    });
}
