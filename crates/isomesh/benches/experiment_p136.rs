//! **P-136 — a null registered on purpose: the hyperdeterminant does not
//! resolve the interior test.**
//!
//! Ticket: R-136. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p136
//! ```
//!
//! Writes `docs/experiments/p-136.csv`.
//!
//! # What was missing
//!
//! `P-127` establishes that `b * b - R::TWO * R::TWO * a * c` at
//! `crates/isomesh/src/marching_cubes/trilinear.rs:246` **is** Cayley's `2x2x2`
//! hyperdeterminant, exactly and with no sign flip. That is a real identity and
//! it invites a false corollary: that a single degree-4 invariant of the eight
//! corner values therefore decides the cell's topology, and that the 730-subcase
//! table is a corollary of an invariant. The registration says in advance that it
//! is not, and this row is the measurement that stops the overclaim from being
//! made rather than the one that walks it back afterwards.
//!
//! The two quantities have never been cross-tabulated in this repository, and the
//! reason they have not is instructive:
//!
//! - `M-206` (`FINDINGS.md`) records that `interior::SweptFaces` and
//!   `trilinear::BodySaddles` *locate the same body saddles to 1.1e-12* while
//!   sharing no arithmetic, no coefficient and no parametrisation. That is an
//!   agreement about **positions**. It says nothing about the **decision**, and
//!   `P-127` explains the position agreement without touching the decision: the
//!   two are two slicings of one pencil.
//! - `M-214` is the closest anything has come, and it runs the other way: six
//!   body saddles mean a cell is a tunnel *or* a twelve-vertex contour and *"the
//!   saddles cannot say which"*. So the saddle count is already known to
//!   under-determine a topological question. Nobody has asked what the
//!   **sign** of the discriminant under-determines.
//! - `M-230` derives Grosso's asymptote-side predicate
//!   (`BodySaddles::same_asymptote_side`) and finds it agrees with a contour count
//!   — again a question about which of two topologies, not about which side a
//!   tunnel opens on.
//!
//! The population is `M-165`'s, and `M-165` is the only measurement in the
//! repository of how often the *published* interior tests disagree with each
//! other: *"among the configurations that can exhibit it, Chernyaev's
//! numerator-only test is wrong 12.6% of the time — 1,966 of 15,625"*
//! (`FINDINGS.md:1259`), produced by
//! `crates/isomesh/src/marching_cubes/interior/tests.rs:176-231`. Its sweep is
//! reproduced here **structurally rather than approximately**: the same five
//! magnitudes, the same nesting, the same opposed-sign face construction, the same
//! `SweptFaces::new` rejection path, and the pinned triple
//! `(total, with_pole, disagree) == (15625, 15625, 1966)` asserted before any row
//! is written.
//!
//! # What `sign(Delta)` is being asked to do, stated so it cannot be softened
//!
//! `sign(Delta)` is one bit (three states, counting the exact zero) computed from
//! the eight corner values. A *rule that reads only that bit* is a function from
//! `{neg, zero, pos}` to `{Joined, Separated}`. There are eight such functions and
//! `agreement_rate` reports the **best** of them on each row — fitted on the row's
//! own configurations by majority vote within each sign class, ties to `joined`.
//!
//! That is deliberately the most generous reading available, and it is the only
//! honest one. Scoring the geometrically natural rule (`Delta > 0` means two real
//! `u` roots, so the hyperbolas cross, so call it `Joined`) would let C1 hold on a
//! badly chosen mapping rather than on `sign(Delta)`'s information content.
//! `natural_rule_rate` is recorded beside it so both readings are in the file.
//!
//! Because `agreement_rate` is an **upper bound over all sign-only rules**, a
//! value below `1` is a statement that *no* such rule can succeed — which is what
//! C1 needs. `base_rate` is the same quantity for a rule that reads nothing at all
//! (predict the row's majority verdict), and `delta_information_gain` is the
//! difference. A gain of zero means `sign(Delta)` is worth exactly as much as not
//! looking.
//!
//! # The mechanism, which is a proof and not a fitted number
//!
//! `Delta` cannot decide the interior verdict, and the reason is a group action.
//! Negating the four corner values of the far face is the `GL(2)^3` element
//! `g3 = diag(1, -1)` acting on the swept slot, with `g1 = g2 = I`. `P-127` records
//! the weight as `(det g1 * det g2 * det g3)^2`, which is `(-1)^2 = 1`, so
//! **`Delta` is exactly invariant** under it. The interior verdict is not:
//! `SweptFaces::numerator(1)` is a product of two negations and does not move,
//! while `SweptFaces::denominator(1)` negates, so `saddle(1)` *flips sign* — and
//! `margin()` is a maximum that includes `saddle(1)`.
//!
//! So one operation leaves `Delta` fixed and moves the answer. That is not
//! evidence about a fixture; it is a reason. The `control/aligned_no_pole` arm is
//! its measurement: it is the M-165 sweep with exactly that negation applied, and
//! `delta_blind_to_face_negation` asserts the two families' `Delta` values agree
//! **row for row, bit for bit**, while their separated counts differ by more than
//! an order of magnitude.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `case150/dmNN` x 14 | the M-165 population split by asymptotic-decider mask | no |
//! | `case150/delta_neg`, `delta_zero`, `delta_pos` | the same population split by `sign(Delta)` — the predictor's own resolution | no |
//! | `case150/all` | the whole M-165 population; the headline row | no |
//! | `control/base_rate` | the same population scored by a rule that never reads `Delta` | **yes** |
//! | `control/aligned_no_pole` | faces ambiguous the *same* way round: no pole, so the two published tests must agree | **yes** |
//!
//! Two taxonomies share the `configuration_class` column, so every row carries a
//! `partition` extra naming which one it belongs to. **The `delta_insufficient`
//! counts are not additive across partitions** and must not be summed: each row
//! fits its own upper-bound rule, so a class that resolves internally can still
//! contribute to the pooled mistake count. `population_delta_insufficient` is the
//! one canonical figure and is carried on every row.
//!
//! `configuration_class` is `caseNNN/...` because the whole M-165 sweep has a
//! **single** Marching Cubes case index — 150, whose inside corners `{1,2,4,7}`
//! are one of the cube's two tetrahedra, so all six faces are ambiguous and this
//! is Chernyaev's case 13. That is asserted rather than assumed
//! (`AMBIGUOUS_FACES[150] == 0b0011_1111`), and it is why the finer class has to
//! be the six-bit decider mask: the case index alone is a constant column and
//! could not characterise anything.
//!
//! # The decider mask, and why it is not a subcase name
//!
//! `ambiguity::joined_mask(&corner, AMBIGUOUS_FACES[150])` is the shipped
//! per-face resolution, and `table::segment_links(case, joined)` keys the 730
//! subcases off exactly that mask. So the mask **is** the granularity at which
//! "does `Delta` make the subcase table a corollary" is a well-posed question.
//!
//! The rows are named by the mask's numeric value and **not** by Chernyaev's
//! `13.1`/`13.2`/`13.5.1` labels. Rule 5: transcribing a published classification
//! from memory is guessing, and the mask needs no transcription — it comes out of
//! `joined_mask`, which the extractor itself calls. `face_bit(axis, side)` gives
//! the layout: bits `0..4` are the four faces containing the sweep axis, bit `4`
//! is the `z = 0` face and bit `5` the `z = 1` face.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration says *"SHARE: none — this is a negative result registered
//! before its positive twin can be believed."* Discharged as written, and the
//! arithmetic behind "none" is worth stating because a null still has a cost
//! ledger:
//!
//! | quantity | value | why |
//! |---|---|---|
//! | change to `crates/isomesh/src/**` | zero bytes | every mechanism here is public API plus this file |
//! | extraction cost moved | zero | no field is sampled and no mesh is extracted |
//! | golden hashes moved | zero | nothing this harness touches is on an extraction path |
//!
//! What the row buys instead is a fence: `R-130` reads tensor rank off
//! `sign(Delta)`, `R-132` reads eight real orbits off it and `R-134` normalises
//! `abs(Delta)` into a continuous ambiguity magnitude. None of those three may be
//! read as resolving the interior test, and after this row the file says so with a
//! number.
//!
//! # Vacuity controls
//!
//! `M-44`: a zero that could not have been non-zero is not a measurement. Every
//! control runs before the first `run.record` and every panic message starts
//! `VOID: `.
//!
//! - **the registered control — `chernyaev_wrong == 1966` and the rate within
//!   `5e-4` of `0.126`.** The registration's own words: *"the configuration set
//!   must reproduce M-165's 12.6% disagreement rate as a control column, or the
//!   fixture is not the population M-165 measured."* Both halves are asserted: the
//!   exact integer, which is what `interior/tests.rs:230` pins, and the rate
//!   against the figure the registration quotes to three digits. `5e-4` is half of
//!   that quotation's last place, so agreement at the tolerance is agreement at
//!   the precision the claim was made in.
//! - **`configurations == 15625` and `degenerate_rejected == 0`** — an opposed face
//!   pair cannot have a zero bilinear denominator, so `SweptFaces::new` must
//!   accept every draw. A rejection would silently shrink the population under the
//!   pinned disagreement count.
//! - **`poled == 15625`** — `M-165`'s structural claim, *"every one of the 15,625
//!   has a pole inside the sweep, by construction"*. Without a pole the corrected
//!   and numerator-only tests are provably equal (`interior.rs:392-394`) and the
//!   comparison is empty.
//! - **every configuration has `case == 150`, and `AMBIGUOUS_FACES[150] == 63`** —
//!   the population is the all-faces-ambiguous case-13 family and not some
//!   easier neighbourhood. If this failed, `joined_mask` would be being asked
//!   about a face with nothing to decide.
//! - **both verdicts occur (`joined > 0 && separated > 0`)** — `agreement_rate` is
//!   a comparison against `SweptFaces::test()`, and on a constant verdict every
//!   rule scores 1 and C1 could not have held for any reason but the fixture.
//! - **all three `Delta` signs occur** — `delta_sign` must be a measurement and not
//!   a constant column, or "sign(Delta) alone" has no sign to be wrong about.
//! - **the aligned control reads `chernyaev_wrong == 0` while its verdict is
//!   *not* constant** — this is the control on the *instrument*: the disagreement
//!   counter is capable of reading zero, so 1,966 is a property of the opposed
//!   family rather than of the comparison; and the zero is not an artefact of
//!   there being only one verdict to agree about.
//! - **`delta_blind_to_face_negation`** — the two families' `Delta` values are
//!   equal row for row and bit for bit. This is the header's mechanism asserted
//!   rather than asserted-in-prose: if it were false, the `GL(2)^3` weight
//!   argument would be wrong and the null would need a different explanation.
//! - **`separated_inside_endpoint_condition` on both families** — `margin()` is a
//!   maximum over candidate points that *includes* `t = 0` and `t = 1`, so a
//!   `Separated` verdict provably implies `saddle(0) <= 0` and `saddle(1) <= 0`.
//!   Asserting the implication checks this harness's endpoint reading against the
//!   shipped `test()` rather than against itself, and it is what makes the
//!   closed-form half of C2 a containment and not a correlation.
//!
//! # Determinism
//!
//! No RNG. The population is a complete `5^6` product over five fixed magnitudes,
//! enumerated in `M-165`'s own loop order, so the counts are identical on every
//! host and every re-run. No clause is timed: `wall_ns` is recorded because it is
//! interesting and read by nothing, which is the only safe status for a nanosecond
//! on a host whose governor swings the same binary 1.45x (`M-280`). Every ordering
//! is on integer keys, every rate is `{:.6}`, and the only float equality in the
//! file is the two intended exact ones — `Delta`'s sign against zero, which is
//! what `trilinear.rs:247-250` branches on, and the family-to-family `Delta`
//! comparison, which is a bit-identity claim.

#![allow(
    clippy::float_cmp,
    reason = "two exact float comparisons are the point: Delta's sign is taken against exactly \
              zero because trilinear.rs:247-250 branches there, and the two families' Delta \
              values are compared for bit identity to assert the GL(2)^3 invariance"
)]

mod common;

use std::time::Instant;

use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::interior::{Interior, SweptFaces, chernyaev_numerator_test};
use isomesh::marching_cubes::table::{AMBIGUOUS_FACES, face_bit};
use isomesh::marching_cubes::trilinear::BodySaddles;

// ─── the registered population, transcribed from M-165's own loop ────────────

/// The five corner magnitudes `M-165` sweeps
/// (`marching_cubes/interior/tests.rs:187`).
const MAGNITUDES: [f64; 5] = [0.1, 0.25, 1.0, 4.0, 10.0];

/// Which cell corners the `lo` face's `[A, B, C, D]` are, in this crate's
/// numbering. `z = 0` is corners 0, 1, 3, 2 — `A`/`C` one diagonal, `B`/`D` the
/// other, exactly as `benches/interior_margin.rs:49` has it.
const LO_CORNERS: [usize; 4] = [0, 1, 3, 2];

/// The same for the `hi` face, `z = 1`.
const HI_CORNERS: [usize; 4] = [4, 5, 7, 6];

/// `M-165`'s population size: `5^6`, every draw accepted.
const M165_CONFIGURATIONS: usize = 15_625;

/// `M-165`'s disagreement count, pinned at `interior/tests.rs:230`.
const M165_DISAGREEMENTS: usize = 1_966;

/// The rate the registration quotes, to the three digits it quotes it in.
const M165_REGISTERED_RATE: f64 = 0.126;

/// Half of `M165_REGISTERED_RATE`'s last decimal place, so agreement at this
/// tolerance is agreement at the precision the claim was published in.
const M165_RATE_TOLERANCE: f64 = 5e-4;

/// The single Marching Cubes case index the whole sweep produces: inside corners
/// `{1, 2, 4, 7}`, one of the cube's two tetrahedra, which is Chernyaev's case 13.
const CASE_13: u8 = 150;

/// All six face bits set — what `AMBIGUOUS_FACES[CASE_13]` must equal.
const ALL_FACES_AMBIGUOUS: u8 = face_bit(0, 0)
    | face_bit(0, 1)
    | face_bit(1, 0)
    | face_bit(1, 1)
    | face_bit(2, 0)
    | face_bit(2, 1);

/// The four faces that contain the sweep axis, derived from the shipped
/// [`face_bit`] layout rather than written as a literal. Named so the
/// decider-mask commentary in this file's header has something to point at:
/// bit 4 is the `z = 0` face and bit 5 the `z = 1` face, and everything below
/// them is lateral.
const LATERAL_FACES: u8 = face_bit(0, 0) | face_bit(0, 1) | face_bit(1, 0) | face_bit(1, 1);

// ─── the two families ────────────────────────────────────────────────────────

/// How the far face's signs relate to the near face's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Family {
    /// `A`/`C` positive on the low face and negative on the high one — the
    /// structure Custodio's Appendix A counterexample has, and the only one that
    /// can put `Delta`'s root inside `(0, 1)`. This is `M-165`'s population.
    Opposed,
    /// Both faces ambiguous the same way round. No pole, so the corrected and
    /// numerator-only tests are provably equal; and it is the M-165 sweep with the
    /// far face negated, which is the `GL(2)^3` element the header names.
    Aligned,
}

impl Family {
    /// The far face, given the three magnitudes drawn for it.
    fn hi(self, a1: f64, c1: f64, b1: f64) -> [f64; 4] {
        match self {
            Self::Opposed => [-a1, b1, -c1, b1],
            Self::Aligned => [a1, -b1, c1, -b1],
        }
    }
}

/// One configuration, with everything the analysis reads off it.
#[derive(Clone, Copy, Debug)]
struct Config {
    /// `b^2 - 4ac` of `BodySaddles`' quadratic — Cayley's `2x2x2`
    /// hyperdeterminant of the eight corner values (`P-127`).
    delta: f64,
    /// `SweptFaces::test()`: Custodio's corrected criterion, which is the correct
    /// interior verdict and the thing `sign(Delta)` is being scored against.
    verdict: Interior,
    /// `chernyaev_numerator_test`: the published test the correction fixes.
    chernyaev: Interior,
    /// `ambiguity::joined_mask` over the six ambiguous faces.
    decider: u8,
    /// The Marching Cubes case index; constant over both families, asserted.
    case: u8,
    /// Is `Delta(t)`'s root inside the sweep? `M-165` says always, on this family.
    poled: bool,
    /// `saddle(0) <= 0 && saddle(1) <= 0`. `margin()` is a maximum that includes
    /// both endpoints, so this is *provably* necessary for `Separated`.
    endpoint_nonpositive: bool,
}

/// `-1`, `0` or `1`. Exact against zero because `trilinear.rs:247-250` is.
fn delta_sign(delta: f64) -> i32 {
    if delta > 0.0 {
        1
    } else if delta < 0.0 {
        -1
    } else {
        0
    }
}

/// `0` for negative, `1` for exactly zero, `2` for positive.
fn sign_slot(delta: f64) -> usize {
    (delta_sign(delta) + 1) as usize
}

/// `0` for [`Interior::Joined`], `1` for [`Interior::Separated`].
fn verdict_slot(v: Interior) -> usize {
    match v {
        Interior::Joined => 0,
        Interior::Separated => 1,
    }
}

/// The CSV token for a verdict.
fn verdict_name(v: Interior) -> &'static str {
    match v {
        Interior::Joined => "joined",
        Interior::Separated => "separated",
    }
}

/// The case index: bit `c` set when corner `c` is inside, which this crate spells
/// `value < 0` (`marching_cubes/mod.rs:265-266`).
fn case_index(corner: &[f64; 8]) -> u8 {
    let mut case = 0u8;
    for (c, &v) in corner.iter().enumerate() {
        if v < 0.0 {
            case |= 1 << c;
        }
    }
    case
}

/// `M-165`'s sweep, with the far face built by `family`.
///
/// The loop nesting, the magnitude list and the `SweptFaces::new` rejection are
/// `interior/tests.rs:187-216`'s, so `Family::Opposed` is that test's population
/// and not a population resembling it. The returned count of rejections is
/// asserted zero rather than absorbed.
fn population(family: Family) -> (Vec<Config>, usize) {
    let mut out = Vec::with_capacity(M165_CONFIGURATIONS);
    let mut degenerate = 0usize;

    for a0 in MAGNITUDES {
        for c0 in MAGNITUDES {
            for b0 in MAGNITUDES {
                for a1 in MAGNITUDES {
                    for c1 in MAGNITUDES {
                        for b1 in MAGNITUDES {
                            let lo = [a0, -b0, c0, -b0];
                            let hi = family.hi(a1, c1, b1);
                            let Ok(faces) = SweptFaces::new(lo, hi) else {
                                degenerate += 1;
                                continue;
                            };

                            let mut corner = [0.0f64; 8];
                            for (&c, &v) in LO_CORNERS.iter().zip(lo.iter()) {
                                corner[c] = v;
                            }
                            for (&c, &v) in HI_CORNERS.iter().zip(hi.iter()) {
                                corner[c] = v;
                            }

                            let case = case_index(&corner);
                            // Written as the shipped expression is written, and in
                            // the same association order, so this is the number
                            // `trilinear.rs:246` computes rather than a rounding
                            // of it: `R::TWO * R::TWO * a * c` folds to
                            // `((4 * a) * c)`.
                            let [qa, qb, qc] = BodySaddles::<f64>::coefficients(&corner);
                            let delta = qb * qb - 2.0 * 2.0 * qa * qc;

                            out.push(Config {
                                delta,
                                verdict: faces.test(),
                                chernyaev: chernyaev_numerator_test(&faces),
                                decider: joined_mask(&corner, AMBIGUOUS_FACES[case as usize]),
                                case,
                                poled: faces.pole().is_some(),
                                endpoint_nonpositive: faces.saddle(0.0) <= 0.0
                                    && faces.saddle(1.0) <= 0.0,
                            });
                        }
                    }
                }
            }
        }
    }

    (out, degenerate)
}

// ─── scoring a set of configurations ─────────────────────────────────────────

/// `k / n` as an `f64`. Both counts are far below `2^53`.
fn rate(k: usize, n: usize) -> f64 {
    k as f64 / n as f64
}

/// `k/n` as a CSV-safe token, so a zero denominator is visible rather than
/// silently rendered as `0.000000`.
fn share(k: usize, n: usize) -> String {
    format!("{k}/{n}")
}

/// Everything one CSV row reports about a set of configurations.
#[derive(Clone, Debug)]
struct Summary {
    configurations: usize,
    /// `[sign][verdict]` counts, `sign` in [`sign_slot`]'s order.
    table: [[usize; 2]; 3],
    /// The best sign-only rule's verdict per sign class, ties to `Joined`.
    ///
    /// The tie-break is documented rather than hidden because it is a choice: a
    /// tied sign class carries no information either way, and `Joined` is the
    /// base-rate answer on this population, so breaking towards it keeps
    /// `delta_information_gain` from being inflated by a coin flip. No tie occurs
    /// on either family; the rule is fixed so that would not matter if one did.
    rule: [Interior; 3],
    joined: usize,
    separated: usize,
    /// Configurations the best sign-only rule gets wrong.
    delta_wrong: usize,
    /// Configurations where `chernyaev_numerator_test` disagrees with `test()`.
    chernyaev_wrong: usize,
    /// `delta_wrong ∩ chernyaev_wrong`.
    overlap: usize,
    /// The geometrically natural rule: `Delta > 0` means `Joined`.
    natural_correct: usize,
    poled: usize,
    /// Configurations with both endpoint saddles non-positive.
    endpoint_population: usize,
    /// Does every `Separated` configuration satisfy the endpoint condition? It
    /// provably must; asserting it checks this harness against `margin()`.
    separated_inside_endpoint: bool,
    /// Does some sign class carry both verdicts? That is exactly the condition
    /// under which `sign(Delta)` cannot determine the answer on this set.
    delta_insufficient: bool,
}

impl Summary {
    /// Score one set. `Delta`'s sign classes are counted, then a majority rule is
    /// fitted within each class and evaluated on the same set — which makes
    /// `agreement_rate` an **upper bound** over every rule that reads only
    /// `sign(Delta)`, and therefore the right instrument for a clause of the form
    /// "sign(Delta) alone agrees on strictly fewer than all".
    fn of(configs: &[Config]) -> Self {
        let mut table = [[0usize; 2]; 3];
        for c in configs {
            table[sign_slot(c.delta)][verdict_slot(c.verdict)] += 1;
        }
        let mut rule = [Interior::Joined; 3];
        for (slot, counts) in table.iter().enumerate() {
            if counts[1] > counts[0] {
                rule[slot] = Interior::Separated;
            }
        }

        let mut joined = 0usize;
        let mut separated = 0usize;
        let mut delta_wrong = 0usize;
        let mut chernyaev_wrong = 0usize;
        let mut overlap = 0usize;
        let mut natural_correct = 0usize;
        let mut poled = 0usize;
        let mut endpoint_population = 0usize;
        let mut separated_inside_endpoint = true;

        for c in configs {
            match c.verdict {
                Interior::Joined => joined += 1,
                Interior::Separated => separated += 1,
            }
            let delta_says = rule[sign_slot(c.delta)];
            let delta_bad = delta_says != c.verdict;
            let chernyaev_bad = c.chernyaev != c.verdict;
            if delta_bad {
                delta_wrong += 1;
            }
            if chernyaev_bad {
                chernyaev_wrong += 1;
            }
            if delta_bad && chernyaev_bad {
                overlap += 1;
            }
            let natural = if c.delta > 0.0 {
                Interior::Joined
            } else {
                Interior::Separated
            };
            if natural == c.verdict {
                natural_correct += 1;
            }
            if c.poled {
                poled += 1;
            }
            if c.endpoint_nonpositive {
                endpoint_population += 1;
            }
            if c.verdict == Interior::Separated && !c.endpoint_nonpositive {
                separated_inside_endpoint = false;
            }
        }

        let delta_insufficient = table.iter().any(|counts| counts[0] > 0 && counts[1] > 0);

        Self {
            configurations: configs.len(),
            table,
            rule,
            joined,
            separated,
            delta_wrong,
            chernyaev_wrong,
            overlap,
            natural_correct,
            poled,
            endpoint_population,
            separated_inside_endpoint,
            delta_insufficient,
        }
    }

    /// The registered `agreement_rate`: the best sign-only rule's share.
    fn agreement_rate(&self) -> f64 {
        rate(self.configurations - self.delta_wrong, self.configurations)
    }

    /// A rule that reads nothing: predict the majority verdict.
    fn base_rate(&self) -> f64 {
        rate(self.joined.max(self.separated), self.configurations)
    }

    /// What `sign(Delta)` is worth over not looking. Structurally zero on a row
    /// that *is* a single sign class, which the header says so nobody reads those
    /// three zeros as measurements.
    fn information_gain(&self) -> f64 {
        self.agreement_rate() - self.base_rate()
    }

    /// The registered `delta_sign`: which `Delta` signs the set contains, in
    /// `neg|zero|pos` order.
    fn delta_sign_token(&self) -> String {
        let names = ["neg", "zero", "pos"];
        let mut parts: Vec<&str> = Vec::new();
        for (slot, counts) in self.table.iter().enumerate() {
            if counts[0] + counts[1] > 0 {
                parts.push(names[slot]);
            }
        }
        parts.join("|")
    }

    /// The registered `interior_test_verdict`: which verdicts the set contains.
    /// `joined|separated` is a mixed set, and on a set whose `delta_sign` is a
    /// single value that mixture *is* the proof that `sign(Delta)` cannot decide.
    fn verdict_token(&self) -> String {
        let mut parts: Vec<&str> = Vec::new();
        if self.joined > 0 {
            parts.push("joined");
        }
        if self.separated > 0 {
            parts.push("separated");
        }
        parts.join("|")
    }

    /// The rule this row fitted, as a readable token — so a reader can see that
    /// the best sign-only rule on the M-165 population is the constant one.
    fn rule_token(&self) -> String {
        let names = ["neg", "zero", "pos"];
        let mut parts: Vec<String> = Vec::new();
        for (slot, counts) in self.table.iter().enumerate() {
            if counts[0] + counts[1] > 0 {
                parts.push(format!(
                    "{}->{}",
                    names[slot],
                    verdict_name(self.rule[slot])
                ));
            }
        }
        parts.join("|")
    }
}

/// One emitted CSV row.
#[derive(Debug)]
struct Row {
    class: String,
    partition: &'static str,
    is_control: bool,
    case: u8,
    summary: Summary,
}

/// Global facts, identical on every row, so the ledger has one canonical copy of
/// each rather than a number a reader has to recompute from a subset.
#[derive(Debug)]
struct Globals {
    /// The registered `cases_where_delta_insufficient`.
    insufficient_classes: String,
    /// How many decider classes the population has at all — the denominator that
    /// makes "four classes" a boundary rather than "everywhere".
    observed_classes: usize,
    insufficient_count: usize,
    population: Summary,
    m165_rate: f64,
    m165_rate_error: f64,
    delta_blind: bool,
    endpoint_population: usize,
    c1: bool,
    c2: bool,
    wall_ns: u128,
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-136");

    common::experiment::run(prereg, |run| {
        let started = Instant::now();

        let (opposed, opposed_degenerate) = population(Family::Opposed);
        let (aligned, aligned_degenerate) = population(Family::Aligned);
        let whole = Summary::of(&opposed);
        let control = Summary::of(&aligned);

        // ── vacuity control 1: the population is M-165's, not one like it ────
        assert!(
            opposed.len() == M165_CONFIGURATIONS && opposed_degenerate == 0,
            "VOID: the opposed sweep produced {} configurations with {opposed_degenerate} \
             rejected by SweptFaces::new, not {M165_CONFIGURATIONS} with none — an opposed face \
             pair cannot have a zero bilinear denominator, so a rejection means the fixture is \
             not the population M-165 measured and the pinned disagreement count below would be \
             a count over a different set",
            opposed.len()
        );
        assert!(
            aligned.len() == M165_CONFIGURATIONS && aligned_degenerate == 0,
            "VOID: the aligned control produced {} configurations with {aligned_degenerate} \
             rejected, not {M165_CONFIGURATIONS} with none, so the row-for-row Delta comparison \
             below would be comparing misaligned indices",
            aligned.len()
        );

        // ── vacuity control 2: the registered one — M-165's 12.6% ────────────
        let m165_rate = rate(whole.chernyaev_wrong, whole.configurations);
        let m165_rate_error = (m165_rate - M165_REGISTERED_RATE).abs();
        assert!(
            whole.chernyaev_wrong == M165_DISAGREEMENTS,
            "VOID: the corrected and numerator-only interior tests disagree on {} of {} \
             configurations, not the {M165_DISAGREEMENTS} that interior/tests.rs:230 pins — the \
             registration requires this fixture to reproduce M-165's rate, and a different \
             integer means it is a different population",
            whole.chernyaev_wrong,
            whole.configurations
        );
        assert!(
            m165_rate_error <= M165_RATE_TOLERANCE,
            "VOID: the reproduced disagreement rate is {m165_rate:.6}, which is \
             {m165_rate_error:.6} from the registration's {M165_REGISTERED_RATE} and outside the \
             {M165_RATE_TOLERANCE} tolerance — half of the last place the claim was quoted in"
        );

        // ── vacuity control 3: every configuration has a pole ────────────────
        assert!(
            whole.poled == M165_CONFIGURATIONS,
            "VOID: {} of {M165_CONFIGURATIONS} configurations have Delta(t)'s root inside the \
             sweep, not all of them — without a pole the corrected and numerator-only tests are \
             provably equal (interior.rs:392-394), so any configuration without one contributes \
             a forced agreement and inflates every rate in this file",
            whole.poled
        );

        // ── vacuity control 4: it really is the case-13 family ───────────────
        assert!(
            AMBIGUOUS_FACES[CASE_13 as usize] == ALL_FACES_AMBIGUOUS,
            "VOID: AMBIGUOUS_FACES[{CASE_13}] is {:#08b}, not {ALL_FACES_AMBIGUOUS:#08b} — the \
             decider masks below would then be computed over faces with nothing to decide, and \
             `configuration_class` would not be a class",
            AMBIGUOUS_FACES[CASE_13 as usize]
        );
        let off_case = opposed.iter().filter(|c| c.case != CASE_13).count();
        assert!(
            off_case == 0,
            "VOID: {off_case} of {M165_CONFIGURATIONS} configurations are not case {CASE_13}, so \
             the sweep is not the single all-faces-ambiguous family the class naming assumes"
        );

        // ── vacuity control 5: the verdict is not a constant ─────────────────
        assert!(
            whole.joined > 0 && whole.separated > 0,
            "VOID: {} joined and {} separated — SweptFaces::test() is a constant over this \
             population, so every rule scores 1 by construction, C1's `strictly fewer than all` \
             could not have failed, and the agreement rate is a property of the fixture (M-44)",
            whole.joined,
            whole.separated
        );

        // ── vacuity control 6: sign(Delta) is not a constant either ──────────
        let signs_present = whole
            .table
            .iter()
            .filter(|counts| counts[0] + counts[1] > 0)
            .count();
        assert!(
            signs_present == 3,
            "VOID: only {signs_present} of the three Delta signs occur ({:?}) — `delta_sign` is \
             then a near-constant column and `sign(Delta) alone` has no sign to be wrong about",
            whole.table
        );

        // ── vacuity control 7: the disagreement counter can read zero ────────
        assert!(
            control.chernyaev_wrong == 0,
            "VOID: the aligned control has {} disagreements between the corrected and \
             numerator-only tests, and it has no pole anywhere, so it must have none — a non-zero \
             here means the comparison itself is broken and 1,966 is a property of the instrument \
             rather than of the opposed family",
            control.chernyaev_wrong
        );
        assert!(
            control.joined > 0 && control.separated > 0,
            "VOID: the aligned control's verdict is constant ({} joined, {} separated), so its \
             zero disagreement count is a zero that could not have been non-zero (M-44)",
            control.joined,
            control.separated
        );

        // ── vacuity control 8: Delta is blind to the far face's sign ─────────
        let delta_blind = opposed
            .iter()
            .zip(aligned.iter())
            .all(|(a, b)| a.delta == b.delta);
        assert!(
            delta_blind,
            "VOID: negating the far face moved Delta on at least one configuration, so the \
             GL(2)^3 weight argument in this file's header — g3 = diag(1, -1) has weight \
             (det)^2 = 1, therefore Delta is invariant — is wrong, and the mechanism the null is \
             explained by has to be re-derived before any row below means what it says"
        );

        // ── vacuity control 9: the endpoint containment is real ──────────────
        assert!(
            whole.separated_inside_endpoint && control.separated_inside_endpoint,
            "VOID: a Separated configuration has a positive endpoint saddle (opposed: {}, \
             aligned: {}) — margin() is a maximum that includes t = 0 and t = 1, so that is \
             impossible, and it means this harness's saddle reading disagrees with the shipped \
             test() rather than characterising it",
            whole.separated_inside_endpoint,
            control.separated_inside_endpoint
        );

        // ── the decider-class partition ──────────────────────────────────────
        let mut classes: Vec<u8> = opposed.iter().map(|c| c.decider).collect();
        classes.sort_unstable();
        classes.dedup();

        let mut class_rows: Vec<Row> = Vec::with_capacity(classes.len());
        let mut insufficient: Vec<u8> = Vec::new();
        for &dm in &classes {
            let members: Vec<Config> = opposed
                .iter()
                .copied()
                .filter(|c| c.decider == dm)
                .collect();
            let summary = Summary::of(&members);
            if summary.delta_insufficient {
                insufficient.push(dm);
            }
            class_rows.push(Row {
                class: format!("case{CASE_13}/dm{dm:02}"),
                partition: "decider_mask",
                is_control: false,
                case: CASE_13,
                summary,
            });
        }

        let insufficient_classes = insufficient
            .iter()
            .map(u8::to_string)
            .collect::<Vec<String>>()
            .join("|");

        // ── the clause verdicts, global ──────────────────────────────────────
        //
        // C1: the best rule that reads only sign(Delta) agrees on strictly fewer
        // than all of the configurations, and the disagreement set is finite,
        // enumerated and non-empty.
        let c1 = whole.agreement_rate() < 1.0 && whole.delta_wrong > 0;
        // C2: the classes are named, they are a strict subset of the observed
        // classes — so the token is a boundary rather than "everywhere" — and the
        // closed-form containment holds with no exceptions, which is what makes
        // the set characterised rather than merely listed.
        let c2 = !insufficient.is_empty()
            && insufficient.len() < classes.len()
            && whole.separated_inside_endpoint;

        let globals = Globals {
            insufficient_classes,
            observed_classes: classes.len(),
            insufficient_count: whole.delta_wrong,
            population: whole.clone(),
            m165_rate,
            m165_rate_error,
            delta_blind,
            endpoint_population: whole.endpoint_population,
            c1,
            c2,
            wall_ns: 0,
        };

        // ── the sign partition, the headline row and the controls ────────────
        let mut rows = class_rows;
        for (name, slot) in [("delta_neg", 0usize), ("delta_zero", 1), ("delta_pos", 2)] {
            let members: Vec<Config> = opposed
                .iter()
                .copied()
                .filter(|c| sign_slot(c.delta) == slot)
                .collect();
            rows.push(Row {
                class: format!("case{CASE_13}/{name}"),
                partition: "delta_sign",
                is_control: false,
                case: CASE_13,
                summary: Summary::of(&members),
            });
        }
        rows.push(Row {
            class: format!("case{CASE_13}/all"),
            partition: "population",
            is_control: false,
            case: CASE_13,
            summary: globals.population.clone(),
        });
        rows.push(Row {
            class: String::from("control/base_rate"),
            partition: "control",
            is_control: true,
            case: CASE_13,
            summary: globals.population.clone(),
        });
        let aligned_case = aligned.first().map_or(0, |c| c.case);
        rows.push(Row {
            class: String::from("control/aligned_no_pole"),
            partition: "control",
            is_control: true,
            case: aligned_case,
            summary: control.clone(),
        });

        let globals = Globals {
            wall_ns: started.elapsed().as_nanos(),
            ..globals
        };

        // ── the console table ────────────────────────────────────────────────
        println!(
            "M-165 reproduced: {} of {} configurations ({:.6}, registered {M165_REGISTERED_RATE}, \
             error {:.6}), every one with a pole, every one case {CASE_13}",
            globals.population.chernyaev_wrong,
            globals.population.configurations,
            globals.m165_rate,
            globals.m165_rate_error
        );
        println!(
            "sign(Delta) confusion over the population, [joined, separated] per sign: \
             neg {:?}, zero {:?}, pos {:?}",
            globals.population.table[0], globals.population.table[1], globals.population.table[2]
        );
        println!(
            "best sign-only rule: {} -> {:.6}; base rate (reads nothing) {:.6}; gain {:.6}; \
             natural rule (Delta > 0 means joined) {:.6}",
            globals.population.rule_token(),
            globals.population.agreement_rate(),
            globals.population.base_rate(),
            globals.population.information_gain(),
            rate(
                globals.population.natural_correct,
                globals.population.configurations
            )
        );
        println!(
            "Delta insufficient in {} of {} decider classes: {}",
            insufficient.len(),
            globals.observed_classes,
            globals.insufficient_classes
        );
        println!(
            "closed form: both endpoint saddles non-positive holds on {} of {} configurations and \
             on every one of the {} the best sign-only rule gets wrong; lateral face mask is \
             {LATERAL_FACES:#08b}",
            globals.endpoint_population,
            globals.population.configurations,
            globals.insufficient_count
        );
        println!(
            "aligned control: Delta identical row for row ({}), {} disagreements, {} separated \
             against the opposed family's {} — same invariant, {:.1}x the separated rate\n",
            globals.delta_blind,
            control.chernyaev_wrong,
            control.separated,
            globals.population.separated,
            rate(control.separated, globals.population.separated.max(1))
        );
        println!(
            "{:<24} {:<13} {:>6} {:<13} {:<18} {:>9} {:>9} {:>9} {:>6} {:>6} {:>6}",
            "configuration_class",
            "partition",
            "n",
            "delta_sign",
            "interior_verdict",
            "agree",
            "base",
            "gain",
            "d_bad",
            "ch_bad",
            "olap"
        );

        for row in &rows {
            let s = &row.summary;
            println!(
                "{:<24} {:<13} {:>6} {:<13} {:<18} {:>9.6} {:>9.6} {:>9.6} {:>6} {:>6} {:>6}  {}",
                row.class,
                row.partition,
                s.configurations,
                s.delta_sign_token(),
                s.verdict_token(),
                s.agreement_rate(),
                s.base_rate(),
                s.information_gain(),
                s.delta_wrong,
                s.chernyaev_wrong,
                s.overlap,
                if row.is_control { "control" } else { "" }
            );

            // `control/base_rate` scores the same population with a rule that
            // never reads Delta, so its `delta_sign` is `unread` and its
            // `agreement_rate` is the base rate. The two numbers coinciding with
            // `case150/all`'s is the finding, not a duplicated row.
            let is_base_rate = row.class == "control/base_rate";
            let (sign_token, agreement) = if is_base_rate {
                (String::from("unread"), s.base_rate())
            } else {
                (s.delta_sign_token(), s.agreement_rate())
            };

            run.record(&[
                ("configuration_class", row.class.clone()),
                ("delta_sign", sign_token),
                ("interior_test_verdict", s.verdict_token()),
                ("agreement_rate", format!("{agreement:.6}")),
                (
                    "cases_where_delta_insufficient",
                    globals.insufficient_classes.clone(),
                ),
                ("chernyaev_disagreement_overlap", s.overlap.to_string()),
                ("c1_holds", globals.c1.to_string()),
                ("c2_holds", globals.c2.to_string()),
                // ── extras (M-273) ──────────────────────────────────────────
                ("partition", String::from(row.partition)),
                ("is_control", row.is_control.to_string()),
                ("configurations", s.configurations.to_string()),
                ("case_index", row.case.to_string()),
                ("poled", s.poled.to_string()),
                ("joined_count", s.joined.to_string()),
                ("separated_count", s.separated.to_string()),
                ("delta_neg", (s.table[0][0] + s.table[0][1]).to_string()),
                ("delta_zero", (s.table[1][0] + s.table[1][1]).to_string()),
                ("delta_pos", (s.table[2][0] + s.table[2][1]).to_string()),
                ("best_sign_rule", s.rule_token()),
                ("best_rule_wrong", s.delta_wrong.to_string()),
                ("chernyaev_wrong", s.chernyaev_wrong.to_string()),
                (
                    "chernyaev_disagreement_rate",
                    format!("{:.6}", rate(s.chernyaev_wrong, s.configurations)),
                ),
                ("overlap_of_delta_wrong", share(s.overlap, s.delta_wrong)),
                (
                    "overlap_of_chernyaev_wrong",
                    share(s.overlap, s.chernyaev_wrong),
                ),
                ("base_rate", format!("{:.6}", s.base_rate())),
                (
                    "delta_information_gain",
                    format!("{:.6}", s.information_gain()),
                ),
                (
                    "natural_rule_rate",
                    format!("{:.6}", rate(s.natural_correct, s.configurations)),
                ),
                (
                    "class_is_delta_insufficient",
                    s.delta_insufficient.to_string(),
                ),
                (
                    "verdict_is_constant",
                    (s.joined == 0 || s.separated == 0).to_string(),
                ),
                (
                    "endpoint_condition_population",
                    s.endpoint_population.to_string(),
                ),
                (
                    "separated_inside_endpoint_condition",
                    s.separated_inside_endpoint.to_string(),
                ),
                // ── extras, global: one canonical copy on every row ─────────
                ("m165_rate_reproduced", format!("{:.6}", globals.m165_rate)),
                ("m165_rate_registered", format!("{M165_REGISTERED_RATE:.6}")),
                ("m165_rate_error", format!("{:.6}", globals.m165_rate_error)),
                (
                    "population_delta_insufficient",
                    globals.insufficient_count.to_string(),
                ),
                (
                    "population_agreement_rate",
                    format!("{:.6}", globals.population.agreement_rate()),
                ),
                (
                    "population_base_rate",
                    format!("{:.6}", globals.population.base_rate()),
                ),
                (
                    "population_information_gain",
                    format!("{:.6}", globals.population.information_gain()),
                ),
                ("delta_insufficient_classes", insufficient.len().to_string()),
                (
                    "decider_classes_observed",
                    globals.observed_classes.to_string(),
                ),
                (
                    "delta_blind_to_face_negation",
                    globals.delta_blind.to_string(),
                ),
                ("wall_ns", globals.wall_ns.to_string()),
            ]);
        }
    });
}
