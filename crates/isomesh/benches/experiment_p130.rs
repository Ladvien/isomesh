//! **P-130 — `Delta > 0` means real tensor rank 2, and the census says whether
//! rank predicts anything.**
//!
//! Ticket: R-130. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p130
//! ```
//!
//! Writes `docs/experiments/p-130.csv`.
//!
//! # What was missing
//!
//! `P-127` established the identity this row stands on: the quantity
//! `crates/isomesh/src/marching_cubes/trilinear.rs:246` computes as
//! `b*b - R::TWO*R::TWO*a*c`, from the coefficients built at `:199-214`, **is**
//! Cayley's `2x2x2` hyperdeterminant of the cell's eight corner values.
//! `common::poly`'s `repo_discriminant().sub(&cayley_2x2x2()).is_zero()` is that
//! proof, and `docs/experiments/p-127.csv` is where it was measured — twelve
//! terms each side, total degree 4, all three axis pairings agreeing.
//!
//! de Silva & Lim, *Tensor rank and the ill-posedness of the best low-rank
//! approximation problem*, `arXiv:math/0607647`, §6, verbatim: *"the rank of a
//! tensor is 2 on the set `{A | Det_2,2,2(A) > 0}` and 3 on the set
//! `{A | Det_2,2,2(A) < 0}`."* Their normalisation is the discriminant of
//! `det(l1*A1 + l2*A2)` — `c1^2 - 4*c0*c2`, no leading minus — which is exactly
//! `common::poly`'s and therefore exactly the crate's. `P-127` measured that
//! agreement rather than assuming it.
//!
//! So the crate has computed a real-tensor-rank certificate on every ambiguous
//! cell since `M-206`, and nothing has ever asked how that certificate is
//! distributed. `M-214`, `M-215`, `M-216` and `M-217` all consume the saddle
//! *count*; none of them censuses the sign. This row is the census, and it is a
//! census of counts rather than of timings — nothing here is a stopwatch.
//!
//! # Which `Delta` this row measures, and why it is the crate's and not the
//! polynomial's
//!
//! Two `f64` routes to the same polynomial:
//!
//! 1. **The crate's.** `BodySaddles::coefficients(&corner)` then
//!    `b*b - 4.0*a*c`, associating exactly as `trilinear.rs:246` does
//!    (`R::TWO*R::TWO*a*c` groups as `((4*a)*c)`, and so does `4.0 * a * c`).
//! 2. **The polynomial's.** `common::poly::cayley_2x2x2().eval_f64(&corner)`,
//!    twelve terms summed over the *expanded* form
//!    (`common/poly.rs:296-308`).
//!
//! `P-127` proved these are the same polynomial exactly. They are **not** the
//! same `f64` computation, and the difference is not academic: on a cell where
//! the field is affine, route 1 returns a **hard zero** (`twist_lo`, `twist_hi`,
//! `dv_lo` and `dv_hi` are all exactly `0.0`, so `a`, `b` and `c` are, so
//! `Delta` is) while route 2 sums twelve terms of equal magnitude and opposite
//! sign whose multiplication orders differ — `f0^2*f7^2` evaluates as
//! `1*A*A*B*B` and `f0*f1*f6*f7` as `1*A*B*A*B` — and need not cancel to `0.0`.
//!
//! The registered clause says *"`Delta = 0` is rare — under 0.1% of surface
//! cells **at `f64`**"*, so the precision is part of the question and the route
//! has to be named. **Route 1 is the measurement**: C1's subject is "a sign the
//! crate already computes", and a stratification taken off a polynomial the
//! crate does not evaluate would be a finding about `common::poly`. Route 2 is
//! recorded as a cross-check in `cayley_sign_disagreements`,
//! `cayley_zero_only_cells` and `crate_zero_only_cells`, which — because the
//! two are the same polynomial, asserted here before any row is written — can
//! only be rounding. An exact sign over exact inputs is `P-133`'s row, not this
//! one.
//!
//! # Rank, computed rather than renamed
//!
//! `Delta != 0` fixes the rank by the theorem. `Delta == 0` does **not**: the
//! degenerate stratum carries rank 0, rank 1, rank 2 *and* rank 3 tensors, and
//! the rank-3 orbit inside it (the tangential one, `W`-state up to `GL(2)^3` —
//! `P-131`'s fixture) is the whole reason the ill-posedness paper exists. So
//! this harness decides the stratum instead of assuming it, in one pass, with
//! every decision an exact `f64` comparison so that "zero" means zero:
//!
//! - all eight corners `0.0` -> **rank 0**;
//! - all three flattenings of rank `<= 1` (every `2x2` minor of each `2x4`
//!   slice pair vanishing) -> **rank 1**;
//! - the pencil quadratic `det(x*M + y*N) = alpha*x^2 + beta*x*y + gamma*y^2`
//!   identically zero -> **rank 2**. Every matrix in the pencil is singular, so
//!   the pencil is a compression space, `T = u (x) W` or `T = W (x) u`, and
//!   `rank T = rank W <= 2`;
//! - otherwise the quadratic has a double root `(x0, y0)`, and
//!   `M0 = x0*M + y0*N` is the singular slice there. Normalising the other slice
//!   to the identity turns `T` into `(B^-1 M0, I)`; that matrix is singular with
//!   a double eigenvalue, hence nilpotent, hence non-diagonalisable **unless it
//!   is zero**. So `M0 == 0` -> **rank 2** (`T = N (x) v`, `det N != 0`) and
//!   `M0 != 0` -> **rank 3**, the tangential orbit.
//!
//! `rank_two_cells >= delta_positive` and `rank_three_cells >= delta_negative`
//! are therefore inequalities and not identities, and both are asserted on
//! every row. `rank_zero_cells` and `rank_one_cells` are extras; the four rank
//! counts sum to `cells` exactly, asserted.
//!
//! The double-root test is taken in pairing `0` (`w`-slices `0123|4567`) and
//! **cross-checked in pairings 1 and 2** on every `Delta == 0` cell. Rank is a
//! tensor invariant, so a disagreement is `f64` incoherence and nothing else;
//! `degenerate_pairing_disagreements` is that count.
//!
//! # C2 is arithmetically unreachable, and the arithmetic is recorded
//!
//! C2 asks that `I(rank; ambiguous)` **exceed** `I(rank; case index)`. In this
//! crate ambiguity is a *lookup on the case index* —
//! `AMBIGUOUS_FACES[case] != 0`, `marching_cubes/table.rs:202-238`, derived
//! from the case and from nothing else. Write `A = g(K)`. The data-processing
//! inequality then gives
//!
//! ```text
//! I(R; A) = I(R; g(K)) <= I(R; K)     for every field, resolution and margin,
//! ```
//!
//! with equality exactly when `g` is a sufficient statistic of `K` for `R`.
//! C2 as registered cannot hold. Per the phase's rule, it is **recorded as
//! unreachable with the arithmetic** and the other two clauses still run:
//! `c2_holds` is `false` on every row, `c2_unreachable` is `true`,
//! `mi_gap_bits = I(R;K) - I(R;A)` is recorded and asserted `>= 0` (a numeric
//! check on the inequality rather than a claim about it), and
//! `c2_margin_required_bits` states the margin the clause would have needed.
//!
//! What *is* reachable is C2's own falsifier — *"rank being a function of the
//! case index, which would make it a renaming with no new signal"* — and that is
//! what `rank_vs_case_index_agreement` decides. It is the accuracy of the best
//! predictor of rank from the case index, `sum_k max_r N(k, r) / cells`, so it
//! reads exactly `1.0` if and only if rank is a function of the case index.
//! `rank_majority_share` (the constant predictor) and `case_index_lift` (the
//! difference) are beside it, because on a population where one rank dominates
//! the agreement is high for a reason that has nothing to do with the case
//! index.
//!
//! **Mutual information is in bits, `log2`**, with the standard convention that
//! a zero joint cell contributes zero. `mi_units` carries the word `bits_log2`
//! in every row so the CSV states its own convention. `mutual_information` is
//! the registered column and is `I(R; A)`; `mi_rank_case_index_bits` is
//! `I(R; K)`. Both are recorded because C2 is a comparison and one number
//! cannot express it. The three marginal entropies are recorded too:
//! `I(R;A) = 0` on a row where `entropy_ambiguous_bits` is `0` is not a finding
//! about rank.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | eight reference fields x 17/33/65 samples | field and resolution — 24 rows | no |
//! | `sphere` at all three resolutions | the field with no ambiguous cell | **yes**, C3's control |
//! | `common::poly` identity re-affirmed before any row | nothing; it licenses the cross-check | **yes**, an assert |
//! | Cayley `eval_f64` beside the crate's `b*b-4ac` | the evaluation route, not the geometry | **yes**, three columns |
//! | pairings 1 and 2 on the `Delta == 0` stratum | which pencil decides the orbit | **yes**, one column |
//!
//! Resolutions are `17`, `33`, `65` **samples** per axis, so `16^3`, `32^3` and
//! `64^3` cells. The registration says "three resolutions" without naming them;
//! 17 and 33 are two of `golden.rs:72`'s `RESOLUTIONS`, so two thirds of the
//! sweep sit on the grids the rest of the repo already reports, and 65 straddles
//! the `u64` word boundary the dual path's bitmap is built on. `129` is left to
//! rows that need it: this one visits 2.4 million cells as it is and nothing
//! here is timed.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE reads *"C1 and C2 are counts, not timings, and move
//! nothing; they decide whether the rank stratification is worth a stage."* That
//! is discharged as written: there is no `Instant` in this file, every clause is
//! decided by integer counts and by mutual informations computed from integer
//! counts, and no shipped path changes. What the row hands forward is a decision
//! for `P-131` and `P-132` — the size of the `Delta = 0` stratum, and whether
//! the case index already knows what rank knows.
//!
//! # The verdict columns, and which of them are global
//!
//! `c1_holds` and `c3_holds` are **global** verdicts, stamped identically on all
//! 24 rows, because neither clause is a property of one row: C1's stability half
//! quantifies over the three resolutions of a field and C3 quantifies over
//! `sphere`. The per-row and per-field facts behind them are separate columns —
//! `c1_zero_bar_holds` and `delta_zero_share` per row, `c1_stability_holds` and
//! `partition_share_range` per field, `c3_scope` and `rank_two_share` for C3.
//! `c2_holds` is `false` on every row for the structural reason above.
//!
//! Both bars are stated in the CSV rather than only here:
//!
//! - `c1_zero_bar = 0.001` — the registration's own "under 0.1% of surface
//!   cells", per row.
//! - `c1_stability_bar = 0.05` — "the partition is stable" read as: for each of
//!   the three classes, the share of surface cells varies by no more than five
//!   percentage points across a field's three resolutions.
//!   `partition_share_range` is the worst of the three ranges, so the bar is
//!   checked against a single number per field.
//! - `c3_rank_two_bar = 0.99` — "rank is 2 on essentially every surface cell",
//!   one cell in a hundred, on every `sphere` row, together with
//!   `ambiguous_cells == 0` there, which is C3's own premise.
//!
//! # The registered vacuity control is half unsatisfiable, and that is measured
//!
//! The registration's control reads *"`gyroid` and `csg_difference` must
//! contribute a non-zero ambiguous-cell count, or C2's mutual information is
//! computed against a constant."* **`csg_difference` cannot contribute one.**
//! Censused before this harness was finished, at seven resolutions, counting
//! surface cells with `AMBIGUOUS_FACES[case] != 0`:
//!
//! | samples | 17 | 25 | 33 | 49 | 65 | 97 | 129 |
//! |---|---|---|---|---|---|---|---|
//! | `sphere` | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
//! | `torus` | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
//! | `box_exact` | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
//! | `csg_difference` | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
//! | `thin_plate` | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
//! | `gyroid` | 52 | 61 | 27 | 93 | 132 | 132 | 153 |
//! | `fbm_terrain` | 16 | 18 | 30 | 38 | 58 | 48 | 52 |
//! | `noise_cavity` | 193 | 409 | 502 | 556 | 567 | 606 | 569 |
//!
//! It is not a resolution artefact and it is not going to become non-zero. An
//! ambiguous face needs **two diagonal** inside corners on one cell face, which
//! means the surface cuts that face into two components.
//! `csg_difference = BoxExact - Sphere{c: 0.6, r: 0.75}` is planar on the box
//! faces, convex on the spherical bite, and joins the two along a **single
//! monotone concave seam curve**; a curve crossing a square cuts it into two
//! *adjacent* pieces, never a diagonal pair. Both its constituents read zero
//! too. The three fields that carry ambiguity are exactly the three whose
//! surfaces are not built from planes and spheres.
//!
//! So the control **fails as registered**, and this harness records that failure
//! with its arithmetic rather than aborting. Aborting would discard 24 rows of
//! C1 and C3 — neither of which the control is about — over a control whose only
//! subject, C2, is separately and structurally unreachable. The failure is in
//! the artefact and not only in this header:
//! `vacuity_control_as_registered_holds` is `false` on every row,
//! `vacuity_csg_difference_ambiguous_cells` carries the zero,
//! `vacuity_gyroid_ambiguous_cells` carries the count that does exist, and
//! `mi_against_a_constant` marks, per row, every row whose ambiguity variable is
//! constant and whose `mutual_information` is therefore structurally zero.
//!
//! # Vacuity controls, as this harness enforces them
//!
//! - **`gyroid` must contribute a non-zero ambiguous-cell count** — the half of
//!   the registered control that is satisfiable. Column:
//!   `vacuity_gyroid_ambiguous_cells`.
//! - **Some row must have a non-constant ambiguity variable**, i.e.
//!   `0 < ambiguous_cells < cells`. Columns: `vacuity_rows_with_ambiguity`,
//!   `mi_against_a_constant`, `entropy_ambiguous_bits`. This is the registered
//!   control's *stated purpose* — "or C2's mutual information is computed
//!   against a constant" — asserted directly and per row, which is strictly
//!   sharper than asserting it of two named fields, because a zero
//!   `mutual_information` on a row with `entropy_ambiguous_bits = 0` is not a
//!   fact about rank (`M-44`).
//! - **At least two distinct ranks must occur over the run.** Column:
//!   `rank_classes_present`. `I(R; anything)` is zero when `R` is constant, and
//!   that zero would say nothing about ambiguity.
//! - **Every row must have at least one surface cell.** Column: `cells`. It is
//!   the denominator of every share and the population of both mutual
//!   informations.
//! - **`common::poly`'s identity must still hold**, i.e.
//!   `repo_discriminant() - cayley_2x2x2()` is the zero polynomial. Columns:
//!   `cayley_sign_disagreements`, `cayley_zero_only_cells`,
//!   `crate_zero_only_cells`. Those three are only interpretable as rounding if
//!   the two routes are algebraically identical; if they were not, the columns
//!   would be measuring an algebra difference while being read as a precision
//!   one.
//! - Not a vacuity control but a precondition, asserted plainly:
//!   `non_finite_samples` is `0` on every row. A `NaN` corner would make every
//!   sign comparison below meaningless in a way no count would reveal.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]

mod common;

use common::poly;
use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::table::{AMBIGUOUS_FACES, is_inside};
use isomesh::marching_cubes::trilinear::BodySaddles;
use isomesh::{Sdf, Shape3, for_each_reference_field};

/// Samples per axis. `n` samples span `n - 1` cells; see the header for why
/// these three.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// Rank classes tracked, `0..=3`. A real `2x2x2` tensor has rank at most 3.
const RANKS: usize = 4;

/// The registration's "under 0.1% of surface cells", per row.
const C1_ZERO_BAR: f64 = 0.001;

/// The worst per-class share range a field may show across its three
/// resolutions and still count as a stable partition.
const C1_STABILITY_BAR: f64 = 0.05;

/// The margin C2 would have had to clear, in bits. Recorded so the CSV states
/// the bar C2 is unreachable against, rather than only this header.
const C2_MARGIN_BITS: f64 = 0.01;

/// "Rank is 2 on essentially every surface cell" on `sphere`: one cell in a
/// hundred.
const C3_RANK_TWO_BAR: f64 = 0.99;

/// The three axis pairings of the eight corners into a `2x2` matrix pencil, and
/// simultaneously the three flattenings of the cell tensor.
///
/// Transcribed from `common/poly.rs:523-527`, which builds
/// `pencil_discriminant` from the same three index sets. That copy is private
/// and exact-`i128` symbolic; this one is the `f64` numeric twin, needed because
/// the orbit test on the degenerate stratum is a numeric decision about one
/// cell's corner values and not a polynomial identity. `0` splits along `w`, `1`
/// along `v`, `2` along `u`, under the crate's corner indexing `f[u + 2v + 4w]`.
const PAIRINGS: [([usize; 4], [usize; 4]); 3] = [
    ([0, 1, 2, 3], [4, 5, 6, 7]),
    ([0, 1, 4, 5], [2, 3, 6, 7]),
    ([0, 2, 4, 6], [1, 3, 5, 7]),
];

/// The three partition classes as accessors, so the stability scan is one loop
/// and the CSV's three shares cannot drift from the three ranges.
const CLASSES: [fn(&Census) -> u64; 3] =
    [|c| c.delta_positive, |c| c.delta_negative, |c| c.delta_zero];

/// One `(field, resolution)` census. Every field is a count or a quantity
/// derived from counts; nothing here is timed.
#[derive(Clone, Debug)]
struct Census {
    /// `(samples - 1)^3`, the whole grid.
    grid_cells: u64,
    /// Surface cells: `case` neither `0` nor `255`. The population of every
    /// count, share and mutual information in this row.
    cells: u64,
    /// Surface cells with `Delta > 0`.
    delta_positive: u64,
    /// Surface cells with `Delta < 0`.
    delta_negative: u64,
    /// Surface cells with `Delta` exactly `0.0`.
    delta_zero: u64,
    /// Surface cells by computed real tensor rank, indexed `0..RANKS`.
    rank_cells: [u64; RANKS],
    /// Surface cells with `AMBIGUOUS_FACES[case] != 0`.
    ambiguous_cells: u64,
    /// How many of the four rank classes are non-empty here.
    rank_classes_present: u32,
    /// `sum_k max_r N(k, r) / cells`: the accuracy of the best predictor of
    /// rank from the case index. Exactly `1.0` iff rank is a function of it.
    agreement: f64,
    /// The accuracy of the best *constant* predictor of rank, the baseline
    /// `agreement` has to beat to mean anything.
    majority_share: f64,
    /// `I(rank; ambiguous)` in bits.
    mi_ambiguous_bits: f64,
    /// `I(rank; case index)` in bits.
    mi_case_bits: f64,
    /// `H(rank)` in bits.
    entropy_rank_bits: f64,
    /// `H(ambiguous)` in bits.
    entropy_ambiguous_bits: f64,
    /// `H(case index)` in bits.
    entropy_case_bits: f64,
    /// Surface cells where the Cayley `eval_f64` and the crate's `b*b-4ac`
    /// disagree in sign class. Pure rounding: the two are the same polynomial.
    cayley_sign_disagreements: u64,
    /// Surface cells where Cayley reads `0.0` and the crate does not.
    cayley_zero_only: u64,
    /// Surface cells where the crate reads `0.0` and Cayley does not.
    crate_zero_only: u64,
    /// `Delta == 0` cells where the three axis pairings do not agree on the
    /// rank. Rank is an invariant, so this is `f64` incoherence.
    degenerate_pairing_disagreements: u64,
    /// Grid samples that were not finite. A precondition, asserted zero.
    non_finite_samples: u64,
    /// The largest `|Delta|` seen, as scale for the zeros.
    max_abs_delta: f64,
    /// The smallest non-zero `|Delta|` seen: how close `f64` came to the
    /// stratum without landing on it.
    min_abs_nonzero_delta: f64,
}

/// `-1`, `0` or `+1`. The sign *class*, which is what the stratification is
/// about; `-0.0` and `+0.0` are both the zero class.
fn sign_class(v: f64) -> i8 {
    if v > 0.0 {
        1
    } else if v < 0.0 {
        -1
    } else {
        0
    }
}

/// `det(x*M + y*N) = alpha*x^2 + beta*x*y + gamma*y^2` for one axis pairing,
/// each index set read row-major into a `2x2` matrix.
///
/// Written in the same association `common/poly.rs:554-560` uses, so the two
/// differ only by `f64` rounding and not by grouping.
fn pencil_quadratic(corner: &[f64; 8], pairing: usize) -> [f64; 3] {
    let (s0, s1) = PAIRINGS[pairing];
    let m: [f64; 4] = std::array::from_fn(|k| corner[s0[k]]);
    let n: [f64; 4] = std::array::from_fn(|k| corner[s1[k]]);

    let alpha = m[0] * m[3] - m[1] * m[2];
    let beta = m[0] * n[3] + n[0] * m[3] - m[1] * n[2] - n[1] * m[2];
    let gamma = n[0] * n[3] - n[1] * n[2];
    [alpha, beta, gamma]
}

/// Whether the `2x4` flattening named by `pairing` has matrix rank at most 1,
/// decided by its six `2x2` minors.
fn flattening_is_rank_one_or_less(corner: &[f64; 8], pairing: usize) -> bool {
    let (s0, s1) = PAIRINGS[pairing];
    for i in 0..4 {
        for j in (i + 1)..4 {
            if corner[s0[i]] * corner[s1[j]] - corner[s0[j]] * corner[s1[i]] != 0.0 {
                return false;
            }
        }
    }
    true
}

/// The real rank of the cell's `2x2x2` corner tensor, in `0..=3`.
///
/// `delta` is the crate's discriminant for this cell. Off the `Delta = 0`
/// stratum the theorem decides it and `pairing` is not consulted; on the
/// stratum the orbit is decided in `pairing`. See the header for the derivation
/// of each branch.
fn real_rank(corner: &[f64; 8], delta: f64, pairing: usize) -> usize {
    if delta > 0.0 {
        return 2;
    }
    if delta < 0.0 {
        return 3;
    }

    if corner.iter().all(|v| *v == 0.0) {
        return 0;
    }
    if (0..3).all(|p| flattening_is_rank_one_or_less(corner, p)) {
        return 1;
    }

    let [alpha, beta, gamma] = pencil_quadratic(corner, pairing);
    if alpha == 0.0 && beta == 0.0 && gamma == 0.0 {
        // Every matrix in the pencil is singular: a compression space, so the
        // tensor factors through a 2x2 matrix and its rank is that matrix's.
        return 2;
    }

    // The double root of `alpha*x^2 + beta*x*y + gamma*y^2`. With `alpha == 0`
    // the quadratic's own discriminant is `beta^2`, so the root is at
    // `(1, 0)` and the singular slice is `M` itself.
    let (x0, y0) = if alpha == 0.0 {
        (1.0, 0.0)
    } else {
        (-beta, 2.0 * alpha)
    };
    let (s0, s1) = PAIRINGS[pairing];
    let singular_slice_is_zero = (0..4).all(|k| x0 * corner[s0[k]] + y0 * corner[s1[k]] == 0.0);
    if singular_slice_is_zero { 2 } else { 3 }
}

/// Shannon entropy of a count vector, in bits.
fn entropy_bits(counts: &[u64]) -> f64 {
    let total: u64 = counts.iter().sum();
    assert!(total > 0, "entropy over an empty population");
    let total_f = total as f64;
    let mut h = 0.0;
    for n in counts {
        if *n == 0 {
            continue;
        }
        let p = *n as f64 / total_f;
        h -= p * p.log2();
    }
    h
}

/// Discrete mutual information in **bits** (`log2`) from a joint count table
/// stored row-major with `cols` columns.
///
/// `I(X;Y) = sum p(x,y) * log2(p(x,y) / (p(x)*p(y)))`, a zero joint cell
/// contributing zero.
fn mutual_information_bits(joint: &[u64], cols: usize) -> f64 {
    assert!(cols > 0, "a joint table needs at least one column");
    assert!(
        joint.len().is_multiple_of(cols),
        "joint table of {} entries is not a multiple of {cols} columns",
        joint.len()
    );
    let total: u64 = joint.iter().sum();
    assert!(total > 0, "mutual information over an empty population");
    let total_f = total as f64;

    let mut col_sums = vec![0u64; cols];
    for row in joint.chunks_exact(cols) {
        for (sum, n) in col_sums.iter_mut().zip(row.iter()) {
            *sum += *n;
        }
    }

    let mut mi = 0.0;
    for row in joint.chunks_exact(cols) {
        let row_sum: u64 = row.iter().sum();
        if row_sum == 0 {
            continue;
        }
        let px = row_sum as f64 / total_f;
        for (n, col_sum) in row.iter().zip(col_sums.iter()) {
            if *n == 0 {
                continue;
            }
            let p = *n as f64 / total_f;
            let py = *col_sum as f64 / total_f;
            mi += p * (p / (px * py)).log2();
        }
    }
    mi
}

/// Census one reference field at one resolution.
fn census<F>(field: &F, samples: u32, cayley: &poly::Poly) -> Census
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, cell_size) = common::grid::<f64, _>(field, samples);
    let size = shape.size();

    // The grid, sampled exactly as `sdf::sample_grid` does at sdf.rs:180-193 —
    // `x` innermost, position `origin + cell_size * index`.
    let mut values = Vec::with_capacity(shape.element_count());
    let mut non_finite_samples = 0u64;
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let v = field.sample([
                    origin[0] + cell_size * f64::from(x),
                    origin[1] + cell_size * f64::from(y),
                    origin[2] + cell_size * f64::from(z),
                ]);
                if !v.is_finite() {
                    non_finite_samples += 1;
                }
                values.push(v);
            }
        }
    }

    let mut joint_case = vec![0u64; 256 * RANKS];
    let mut joint_ambiguous = vec![0u64; 2 * RANKS];
    let mut cells = 0u64;
    let mut delta_positive = 0u64;
    let mut delta_negative = 0u64;
    let mut delta_zero = 0u64;
    let mut ambiguous_cells = 0u64;
    let mut cayley_sign_disagreements = 0u64;
    let mut cayley_zero_only = 0u64;
    let mut crate_zero_only = 0u64;
    let mut degenerate_pairing_disagreements = 0u64;
    let mut max_abs_delta = 0.0f64;
    let mut min_abs_nonzero_delta = f64::INFINITY;

    for z in 0..size[2] - 1 {
        for y in 0..size[1] - 1 {
            for x in 0..size[0] - 1 {
                // The crate's own gather and case construction, corner `k` at
                // offset `(k & 1, (k >> 1) & 1, (k >> 2) & 1)` — the indexing
                // `f[u + 2v + 4w]` that `common::poly` is written in.
                let mut corner = [0.0f64; 8];
                let mut case = 0u8;
                for (k, slot) in corner.iter_mut().enumerate() {
                    let p = [
                        x + (k as u32 & 1),
                        y + ((k as u32 >> 1) & 1),
                        z + ((k as u32 >> 2) & 1),
                    ];
                    let v = values[shape.linearize(p) as usize];
                    *slot = v;
                    if is_inside(v) {
                        case |= 1u8 << k;
                    }
                }
                if case == 0 || case == 255 {
                    continue;
                }
                cells += 1;

                let [a, b, c] = BodySaddles::coefficients(&corner);
                let delta = b * b - 4.0 * a * c;
                let delta_class = sign_class(delta);
                match delta_class {
                    1 => delta_positive += 1,
                    -1 => delta_negative += 1,
                    _ => delta_zero += 1,
                }
                let magnitude = delta.abs();
                if magnitude > max_abs_delta {
                    max_abs_delta = magnitude;
                }
                if delta_class != 0 && magnitude < min_abs_nonzero_delta {
                    min_abs_nonzero_delta = magnitude;
                }

                let cayley_class = sign_class(cayley.eval_f64(&corner));
                if cayley_class != delta_class {
                    cayley_sign_disagreements += 1;
                    if cayley_class == 0 {
                        cayley_zero_only += 1;
                    }
                    if delta_class == 0 {
                        crate_zero_only += 1;
                    }
                }

                let rank = real_rank(&corner, delta, 0);
                if delta_class == 0
                    && (real_rank(&corner, delta, 1) != rank
                        || real_rank(&corner, delta, 2) != rank)
                {
                    degenerate_pairing_disagreements += 1;
                }

                let ambiguous = usize::from(AMBIGUOUS_FACES[case as usize] != 0);
                ambiguous_cells += ambiguous as u64;
                joint_case[case as usize * RANKS + rank] += 1;
                joint_ambiguous[ambiguous * RANKS + rank] += 1;
            }
        }
    }

    assert!(
        cells > 0,
        "no surface cell at {samples}^3, so every count below is over an empty \
         population"
    );

    // `as_chunks` over the const `RANKS` rather than `chunks_exact`: the row
    // type is then `&[u64; RANKS]` and the column count cannot drift from the
    // rank count the joint tables were built with.
    let case_rows = joint_case.as_chunks::<RANKS>().0;
    let ambiguous_rows = joint_ambiguous.as_chunks::<RANKS>().0;

    let mut rank_cells = [0u64; RANKS];
    for (r, slot) in rank_cells.iter_mut().enumerate() {
        *slot = case_rows.iter().map(|row| row[r]).sum::<u64>();
    }
    let case_counts: Vec<u64> = case_rows.iter().map(|row| row.iter().sum()).collect();
    let ambiguous_counts: Vec<u64> = ambiguous_rows.iter().map(|row| row.iter().sum()).collect();

    let best_per_case: u64 = case_rows
        .iter()
        .map(|row| row.iter().copied().max().unwrap_or(0))
        .sum();
    let cells_f = cells as f64;

    Census {
        grid_cells: u64::from(samples - 1).pow(3),
        cells,
        delta_positive,
        delta_negative,
        delta_zero,
        rank_cells,
        ambiguous_cells,
        rank_classes_present: rank_cells.iter().filter(|n| **n > 0).count() as u32,
        agreement: best_per_case as f64 / cells_f,
        majority_share: rank_cells.iter().copied().max().unwrap_or(0) as f64 / cells_f,
        mi_ambiguous_bits: mutual_information_bits(&joint_ambiguous, RANKS),
        mi_case_bits: mutual_information_bits(&joint_case, RANKS),
        entropy_rank_bits: entropy_bits(&rank_cells),
        entropy_ambiguous_bits: entropy_bits(&ambiguous_counts),
        entropy_case_bits: entropy_bits(&case_counts),
        cayley_sign_disagreements,
        cayley_zero_only,
        crate_zero_only,
        degenerate_pairing_disagreements,
        non_finite_samples,
        max_abs_delta,
        min_abs_nonzero_delta,
    }
}

/// Total ambiguous cells one field contributed across its three resolutions.
fn ambiguous_total(rows: &[(&'static str, u32, Census)], field: &str) -> u64 {
    rows.iter()
        .filter(|(name, _, _)| *name == field)
        .map(|(_, _, c)| c.ambiguous_cells)
        .sum()
}

/// The worst share range, over the three partition classes, across one field's
/// three resolutions. This is what `C1_STABILITY_BAR` is checked against.
fn partition_share_range(rows: &[(&'static str, u32, Census)], field: &str) -> f64 {
    let mut worst = 0.0f64;
    for class in CLASSES {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for (_, _, c) in rows.iter().filter(|(name, _, _)| *name == field) {
            let share = class(c) as f64 / c.cells as f64;
            lo = lo.min(share);
            hi = hi.max(share);
        }
        worst = worst.max(hi - lo);
    }
    worst
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-130");

    common::experiment::run(prereg, |run| {
        let cayley = poly::cayley_2x2x2();

        // ── vacuity control: the identity the cross-check columns rest on ────
        assert!(
            poly::repo_discriminant().sub(&cayley).is_zero(),
            "VOID: `repo_discriminant - cayley_2x2x2` is not the zero polynomial, so \
             `cayley_sign_disagreements` would be counting an algebra difference while \
             being read as a rounding one, and P-127's identity — which is the whole \
             licence for calling the crate's sign a tensor-rank certificate — does not \
             hold on this tree"
        );

        let mut rows: Vec<(&'static str, u32, Census)> = Vec::new();
        for_each_reference_field!(f64, |name, field| {
            for samples in RESOLUTIONS {
                rows.push((name, samples, census(&field, samples, &cayley)));
            }
        });

        // ── vacuity controls, all before the first record ────────────────────
        for (name, samples, c) in &rows {
            assert_eq!(
                c.non_finite_samples, 0,
                "{name} at {samples}^3 produced {} non-finite samples, so every sign \
                 comparison in this row is meaningless",
                c.non_finite_samples
            );
            assert!(
                c.cells > 0,
                "VOID: {name} at {samples}^3 has no surface cell, so its shares have no \
                 denominator and both mutual informations are over an empty population"
            );
        }
        // The registered control names two fields. `gyroid` satisfies it;
        // `csg_difference` is measured to be incapable of satisfying it at any
        // resolution, so its half is **recorded as failed with its arithmetic**
        // rather than aborting 24 rows of C1 and C3. See the header.
        let gyroid_ambiguous = ambiguous_total(&rows, "gyroid");
        let csg_ambiguous = ambiguous_total(&rows, "csg_difference");
        assert!(
            gyroid_ambiguous > 0,
            "VOID: gyroid contributed no ambiguous cell across {RESOLUTIONS:?}, so no row \
             in this run has a non-constant ambiguity variable and every \
             `mutual_information` here is a zero that could not have been non-zero (M-44)"
        );
        let rows_with_ambiguity = rows
            .iter()
            .filter(|(_, _, c)| c.ambiguous_cells > 0 && c.ambiguous_cells < c.cells)
            .count();
        assert!(
            rows_with_ambiguity > 0,
            "VOID: no row has a non-constant ambiguity variable, so `I(rank; ambiguous)` is \
             identically zero for a reason that has nothing to do with rank (M-44)"
        );
        let ranks_seen: usize = (0..RANKS)
            .filter(|r| rows.iter().any(|(_, _, c)| c.rank_cells[*r] > 0))
            .count();
        assert!(
            ranks_seen >= 2,
            "VOID: only {ranks_seen} rank class occurs anywhere in the run, so `I(rank; ·)` \
             is zero because rank is constant and says nothing about ambiguity"
        );

        // ── clause verdicts ─────────────────────────────────────────────────
        let stability: Vec<(&'static str, f64)> = rows
            .iter()
            .map(|(name, _, _)| *name)
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .map(|name| (name, partition_share_range(&rows, name)))
            .collect();
        let stable = |field: &str| -> f64 {
            stability
                .iter()
                .find(|(name, _)| *name == field)
                .map(|(_, range)| *range)
                .expect("every field was censused")
        };

        let c1 = rows.iter().all(|(name, _, c)| {
            (c.delta_zero as f64 / c.cells as f64) < C1_ZERO_BAR && stable(name) <= C1_STABILITY_BAR
        });
        // C2 is unreachable: `ambiguous = g(case)`, so `I(R;A) <= I(R;K)` by the
        // data-processing inequality and no margin can be cleared.
        let c2 = false;
        let c3 = rows
            .iter()
            .filter(|(name, _, _)| *name == "sphere")
            .all(|(_, _, c)| {
                c.ambiguous_cells == 0
                    && (c.rank_cells[2] as f64 / c.cells as f64) >= C3_RANK_TWO_BAR
            });

        for (name, samples, c) in &rows {
            let cells_f = c.cells as f64;
            let mi_gap = c.mi_case_bits - c.mi_ambiguous_bits;
            assert!(
                mi_gap >= -1e-9,
                "{name} at {samples}^3: I(rank; ambiguous) = {} exceeds \
                 I(rank; case index) = {}, which the data-processing inequality forbids \
                 because ambiguity is a lookup on the case index — the mutual-information \
                 code is wrong, not the theorem",
                c.mi_ambiguous_bits,
                c.mi_case_bits
            );
            assert_eq!(
                c.delta_positive + c.delta_negative + c.delta_zero,
                c.cells,
                "{name} at {samples}^3: the sign classes do not partition the surface cells"
            );
            assert_eq!(
                c.rank_cells.iter().sum::<u64>(),
                c.cells,
                "{name} at {samples}^3: the rank classes do not partition the surface cells"
            );
            assert!(
                c.rank_cells[2] >= c.delta_positive,
                "{name} at {samples}^3: {} rank-2 cells against {} with Delta > 0, but the \
                 theorem makes every Delta > 0 cell rank 2",
                c.rank_cells[2],
                c.delta_positive
            );
            assert!(
                c.rank_cells[3] >= c.delta_negative,
                "{name} at {samples}^3: {} rank-3 cells against {} with Delta < 0, but the \
                 theorem makes every Delta < 0 cell rank 3",
                c.rank_cells[3],
                c.delta_negative
            );

            let share_of_case = if c.mi_case_bits > 0.0 {
                c.mi_ambiguous_bits / c.mi_case_bits
            } else {
                // `I(R;K) == 0` happens exactly when rank is constant on the row,
                // which `rank_classes_present` reads; the share is then 0 by
                // definition rather than undefined.
                0.0
            };

            run.record(&[
                ("field", (*name).to_string()),
                ("resolution", samples.to_string()),
                ("cells", c.cells.to_string()),
                ("delta_positive", c.delta_positive.to_string()),
                ("delta_negative", c.delta_negative.to_string()),
                ("delta_zero", c.delta_zero.to_string()),
                ("rank_two_cells", c.rank_cells[2].to_string()),
                ("rank_three_cells", c.rank_cells[3].to_string()),
                ("ambiguous_cells", c.ambiguous_cells.to_string()),
                (
                    "rank_vs_case_index_agreement",
                    format!("{:.6}", c.agreement),
                ),
                ("mutual_information", format!("{:.6}", c.mi_ambiguous_bits)),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // ── extras (M-273) ──
                (
                    "ambiguous_share",
                    format!("{:.6}", c.ambiguous_cells as f64 / cells_f),
                ),
                ("c1_stability_bar", format!("{C1_STABILITY_BAR:.6}")),
                (
                    "c1_stability_holds",
                    (stable(name) <= C1_STABILITY_BAR).to_string(),
                ),
                ("c1_zero_bar", format!("{C1_ZERO_BAR:.6}")),
                (
                    "c1_zero_bar_holds",
                    ((c.delta_zero as f64 / cells_f) < C1_ZERO_BAR).to_string(),
                ),
                ("c2_margin_required_bits", format!("{C2_MARGIN_BITS:.6}")),
                ("c2_unreachable", true.to_string()),
                ("c3_rank_two_bar", format!("{C3_RANK_TWO_BAR:.6}")),
                (
                    "c3_scope",
                    if *name == "sphere" { "sphere" } else { "other" }.to_string(),
                ),
                (
                    "case_index_lift",
                    format!("{:.6}", c.agreement - c.majority_share),
                ),
                (
                    "cayley_sign_disagreements",
                    c.cayley_sign_disagreements.to_string(),
                ),
                ("cayley_zero_only_cells", c.cayley_zero_only.to_string()),
                ("crate_zero_only_cells", c.crate_zero_only.to_string()),
                (
                    "degenerate_pairing_disagreements",
                    c.degenerate_pairing_disagreements.to_string(),
                ),
                (
                    "delta_negative_share",
                    format!("{:.6}", c.delta_negative as f64 / cells_f),
                ),
                (
                    "delta_positive_share",
                    format!("{:.6}", c.delta_positive as f64 / cells_f),
                ),
                (
                    "delta_zero_share",
                    format!("{:.6}", c.delta_zero as f64 / cells_f),
                ),
                (
                    "entropy_ambiguous_bits",
                    format!("{:.6}", c.entropy_ambiguous_bits),
                ),
                (
                    "entropy_case_index_bits",
                    format!("{:.6}", c.entropy_case_bits),
                ),
                ("entropy_rank_bits", format!("{:.6}", c.entropy_rank_bits)),
                ("grid_cells", c.grid_cells.to_string()),
                ("max_abs_delta", format!("{:.6e}", c.max_abs_delta)),
                ("mi_ambiguous_share_of_case", format!("{share_of_case:.6}")),
                ("mi_gap_bits", format!("{mi_gap:.6}")),
                ("mi_rank_case_index_bits", format!("{:.6}", c.mi_case_bits)),
                ("mi_units", String::from("bits_log2")),
                (
                    "min_abs_nonzero_delta",
                    format!("{:.6e}", c.min_abs_nonzero_delta),
                ),
                ("non_finite_samples", c.non_finite_samples.to_string()),
                ("partition_share_range", format!("{:.6}", stable(name))),
                ("rank_classes_present", c.rank_classes_present.to_string()),
                (
                    "rank_is_a_function_of_case_index",
                    (c.agreement == 1.0).to_string(),
                ),
                ("rank_majority_share", format!("{:.6}", c.majority_share)),
                ("rank_one_cells", c.rank_cells[1].to_string()),
                (
                    "rank_three_share",
                    format!("{:.6}", c.rank_cells[3] as f64 / cells_f),
                ),
                (
                    "rank_two_share",
                    format!("{:.6}", c.rank_cells[2] as f64 / cells_f),
                ),
                ("rank_zero_cells", c.rank_cells[0].to_string()),
                (
                    "surface_share",
                    format!("{:.6}", cells_f / c.grid_cells as f64),
                ),
                (
                    "mi_against_a_constant",
                    (c.ambiguous_cells == 0 || c.ambiguous_cells == c.cells).to_string(),
                ),
                (
                    "vacuity_control_as_registered_holds",
                    (gyroid_ambiguous > 0 && csg_ambiguous > 0).to_string(),
                ),
                (
                    "vacuity_csg_difference_ambiguous_cells",
                    csg_ambiguous.to_string(),
                ),
                (
                    "vacuity_gyroid_ambiguous_cells",
                    gyroid_ambiguous.to_string(),
                ),
                (
                    "vacuity_rows_with_ambiguity",
                    rows_with_ambiguity.to_string(),
                ),
            ]);
        }
    });
}
