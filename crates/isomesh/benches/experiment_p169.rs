//! **P-169 — which corner matters, as an influence rather than an intuition.**
//!
//! Ticket: R-169. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p169
//! ```
//!
//! Writes `docs/experiments/p-169.csv`.
//!
//! # What was missing
//!
//! `validate_table()` (`crates/isomesh/src/marching_cubes/mod.rs:836`) audits the
//! 256-case table **combinatorially** — cut-edge agreement, face consistency,
//! manifoldness. It never asks a *quantitative* question about the function, and
//! in particular it never asks whether the eight corners are treated alike. The
//! influence of a variable, `Inf_i(f) = Pr[f(x) != f(x + e_i)]`, is exactly that
//! question, and it is a theorem rather than a heuristic: a table that is
//! octahedrally correct **must** give all eight corners the same influence, so an
//! unequal influence localises a defect.
//!
//! `R-167` built the instrument (`benches/common/boolean.rs`, which that ticket
//! owns and this one consumes unchanged) and `docs/experiments/p-167.csv` already
//! carries the influences as extra columns. This row is the clause that reads
//! them as a check, and the second half of the registration — is the number any
//! use to a refinement heuristic — has never been asked at all.
//!
//! # Two measured facts from `p-167.csv` that decide what the clauses mean
//!
//! Both are read from the committed CSV, not from a summary, and both are
//! **re-measured here** so this artefact stands on its own.
//!
//! **1. `corner_symmetry_classes` is `1`.** Every row of `p-167.csv` carries it.
//! `common::boolean::corner_symmetry_classes()` generates the orbits of the
//! 48-element cube group acting on the eight corner indices by flood fill, and
//! the answer is a **single orbit** — the three axis sign-flips alone already
//! carry corner 0 to all eight. So C1's *"equal within each octahedral symmetry
//! class"* is not a four-against-four constraint with a gap a defect could hide
//! in. It is the **strongest reading the clause could have: all eight influences
//! must be equal to one another**, and any inequality at all is a violation. That
//! is worth saying out loud, because had the corners split into two classes the
//! clause would have been half as sharp and the same words would have looked
//! identical.
//!
//! **2. The influences are all-equal on `triangle_counts` and unequal on
//! `edge_masks`.** `p-167.csv` reports `influences_all_equal = true` on all four
//! bits of `triangle_counts` (`influence_by_corner` is `1.0` × 8, `0.5` × 8,
//! `0.4375` × 8, `0.0` × 8, and `total_influence` is `8.0`, `4.0`, `3.5`, `0.0`)
//! and `false` on `edge_masks.bit0`–`.bit11` (each `1|1|0|0|0|0|0|0` up to a
//! relabelling, `total_influence = 2.0`).
//!
//! **The edge-mask inequality is not a defect, and reading C1 against it would be
//! a false positive.** Here is the arithmetic. A cube symmetry `pi` relabels the
//! corner signs, `s'[i] = s[pi[i]]`, and it maps edges to edges: the pair
//! `{pi[a], pi[b]}` of edge `e = {a, b}` is again an edge `sigma(e)`. So
//! `mask(x') = sigma(mask(x))` — the mask is **equivariant**, and the group acts
//! on the twelve edge *labels*. A single edge bit is therefore not an invariant
//! function: `bit_e` is literally `x_a xor x_b`, attached to two named corners,
//! and its influence vector is `1` on those two and `0` on the other six under
//! *any* measure. Only the **multiset over bits** is invariant.
//!
//! Both halves of that are measured rather than argued, using
//! `common::poly::octahedral_relabellings()` — an independently generated copy of
//! the same 48-element group, asserted distinct and closed under composition by
//! its own module:
//!
//! - `octahedral_violations` counts `(pi, x)` with
//!   `triangle_counts[x'] != triangle_counts[x]`. **Predicted 0**: the primary
//!   reading is invariant as a labelled function, which is *why* C1's equality is
//!   a consequence of anything.
//! - `edge_mask_labelled_violations` counts `(pi, x)` with
//!   `edge_masks[x'] != edge_masks[x]`. **Predicted large and non-zero**: the
//!   labelled reading is not invariant.
//! - `edge_mask_equivariant_violations` counts `(pi, x)` with
//!   `edge_masks[x'] != sigma(edge_masks[x])`. **Predicted 0**: the failure above
//!   is exactly the edge relabelling and nothing else.
//!
//! So **C1 is read against `shipped_triangle_counts`**, and the reason is a
//! measured property of that reading rather than a preference. The multiset-level
//! statement the edge masks *do* satisfy is measured too, as its own arm — see
//! *The aggregate arms* below.
//!
//! # Three readings, and which is which
//!
//! `common::boolean` exposes three readings of `CASES`; all three are run, and
//! each arm's `reading`, `role` and `is_degenerate` columns say what it turned out
//! to be. Bit counts follow P-167's rule, `bit_length(max) + 1`, so the two CSVs
//! analyse the same functions and join on `output_bit`. The extra top bit is the
//! **constant-zero witness**: its measured degeneracy proves the reading's output
//! really stops where the arithmetic says it does.
//!
//! | reading | bits | what its influences are | C1's reading? |
//! |---|---|---|---|
//! | `shipped_triangle_counts` | 4 | all eight equal on every bit | **yes — the primary** |
//! | `shipped_edge_masks` | 13 | `1` on two named corners, `0` on six | no — equivariant, not invariant |
//! | `shipped_centroid_counts` | 1 | all zero; the reading is the constant | no — vacuously equal |
//!
//! ## The aggregate arms
//!
//! Each reading also gets an `aggregate` arm: the per-corner sum of the influence
//! over all of the reading's bits. That is not a convenience, it is **the
//! vector-valued influence** — `E[number of output bits that flip when corner i
//! flips]`, the expected Hamming displacement of the integer output — and it is
//! the level at which the group actually acts, so it is the honest form of the
//! multiset statement.
//!
//! Its value for the edge masks is known in closed form before it is computed:
//! corner `i` gets `1` from each of the twelve edge bits incident to it, and every
//! corner of a cube has [`CUBE_VERTEX_DEGREE`] `= 3` incident edges, so the
//! aggregate is **exactly `3.0` on all eight corners** and the total is `24 = 2 ×
//! 12`. That makes it a *calibration* — a third function with a known answer
//! beside R-167's parity and majority — and it is asserted, not merely reported.
//!
//! # C1, and what its verdict is taken over
//!
//! C1 asks that all eight influences be computed and be equal within class. Both
//! halves are measured per row: `influence` is each corner's number and
//! `influence_equal_within_class` is the equality **for that arm**. Because a
//! clause verdict must not read differently on rows the clause is not about,
//! `c1_holds` is a **global** verdict — the same value on every row, as the
//! authoring contract provides for — taken over the primary reading's
//! non-degenerate bit arms plus its aggregate. Degenerate arms are excluded
//! because an equality of zeros is not an equality (M-44), and the exclusion is
//! guarded: a vacuity control requires a non-degenerate primary bit to exist.
//!
//! `influence` is cross-checked against `influence_combinatorial` on every corner
//! of every arm and the agreement is recorded as
//! `influence_agrees_combinatorial`. The module asserts it internally; recording
//! it puts the check in the artefact instead of in a panic that never fired.
//!
//! # C2, and the two independent ways the number is uninformative
//!
//! C2 asks that total influence be *reported and compared against the average
//! sensitivity a refinement heuristic would need*. Total influence **is** average
//! sensitivity — `sum_S |S| fhat(S)^2`, the expected number of pivotal corners at
//! a random input — so the comparison has to be against a heuristic's *actual*
//! quantity, and there are two of those. Both are measured.
//!
//! **(a) Does the per-corner influence rank the corners the way a refinement
//! heuristic would?** `refinement_priority_correlation` is
//! `common::beta::rank_correlation` — Spearman, ties averaged, `total_cmp`
//! ordering — between the arm's eight influences and a per-corner **refinement
//! priority** accumulated over the eight reference fields at
//! [`RESOLUTIONS`]. The priority is the registration's own suggestion, the
//! extracted vertex's displacement: on each cut edge `{a, b}` the Marching Cubes
//! vertex sits at `t = v_a / (v_a - v_b)`, so
//!
//! ```text
//! priority[a] += min(1, d * |v_b| / (v_a - v_b)^2)
//! priority[b] += min(1, d * |v_a| / (v_a - v_b)^2),   d = PERTURBATION * cell_size
//! ```
//!
//! which is the linearised movement of that vertex when corner `a`'s sample is
//! perturbed by the amount one step of grid refinement would remove, clamped to
//! the edge it cannot leave. `v_a - v_b` is strictly non-zero on a cut edge — one
//! endpoint is strictly negative and the other is not — so there is no singular
//! case and no second path; the clamp is part of the definition, and it is also
//! what keeps a near-tangential crossing from owning the sum.
//!
//! **The clause's verdict is decided by the arithmetic of C1 and it is worth
//! seeing that before the run.** `rank_correlation` returns `0.0` when either
//! sample has no rank variance (`beta.rs:854-874`: a zero denominator). The
//! primary reading's influences are all *equal*, so that sample has no rank
//! variance and the correlation is **exactly zero on every primary arm** — not
//! small, zero, and for a reason. **C1 holding is C2 failing.** The table's
//! octahedral correctness is precisely what makes per-corner influence useless as
//! a refinement signal, and that is the registration's own expectation: *"C2 by
//! total influence being uninformative, which is the expected outcome for a
//! symmetric table."* `c2_holds` is `|correlation| >=`
//! [`INFORMATIVE_CORRELATION`] on the primary reading, and is predicted **false**.
//!
//! A zero that could not have been non-zero would be worth nothing, so the
//! priority side is controlled three ways: it must have at least
//! [`MIN_DISTINCT_PRIORITIES`] distinct values, and correlating it with itself and
//! with its own negation must return `+1` and `-1` through the same function.
//! With those, a zero is attributable to the influence side and to nothing else.
//!
//! The priority is also *reported*, with its spread, because there is a second
//! result hiding in it: a per-corner statistic accumulated over a translation-
//! invariant grid is **forced** to be nearly corner-symmetric. The same sample
//! plays all eight corner roles in the eight cells that meet it, so `priority[i]`
//! is a sum of one octant's worth of response per sample and the eight sums can
//! differ only by the fields' own asymmetry.
//!
//! **That residual does not wash out with refinement, and the first reading of
//! this was wrong.** It looks like a boundary-layer effect, which would shrink
//! like `1/N`; it is not. Measured, `refinement_priority_spread` as a share of
//! the vector's own mean is **1.52% at `33³`** and **2.68% at `65³`** — it *grows*
//! — while the pooled vector over both resolutions sits at **0.58%**, the two
//! grids' asymmetries partly cancelling. So the near-degeneracy is a property of
//! the eight fields rather than of the grid's edge, and a finer grid does not buy
//! a ranking. Both per-resolution spreads are printed and the pooled one is a
//! column, so a reader can see exactly how near-constant the sample the
//! correlation is taken against actually is. Nothing here is dressed up as signal.
//!
//! **(b) Is the uniform-measure number the number a heuristic would face?**
//! Influence is defined under the uniform measure on `{0,1}^8`. A refinement
//! heuristic runs on real cells, whose case histogram is nothing like uniform —
//! it is dominated by the empty and full cases. So `empirical_influence` is the
//! same flip count re-weighted by the measured histogram of all eight reference
//! fields at both resolutions, and `empirical_total_influence` is the average
//! sensitivity a heuristic would actually see. Predicted, and the prediction is
//! the interesting part: any function that is a `xor` of corner signs has a
//! **measure-independent** influence, so parity (`triangle_counts.bit0`, `8.0`)
//! and every edge bit (`2.0`) must agree exactly, while the two *informative*
//! bits of the primary reading — the only ones C2 could have used — must not. The
//! quantity is doubly uninformative: it does not rank the corners, and it is not
//! even the right magnitude for the cells a caller has.
//!
//! # Arms
//!
//! | arm | what it varies | `is_control` |
//! |---|---|---|
//! | `triangle_counts.bit0` … `.bit3` | the primary reading, one arm per output bit | no |
//! | `triangle_counts.aggregate` | the primary reading's vector-valued influence | no |
//! | `edge_masks.bit0` … `.bit12` | a reading that is equivariant but not invariant | no |
//! | `edge_masks.aggregate` | the multiset statement, closed form `3.0` per corner | no |
//! | `centroid_counts.bit0`, `.aggregate` | the constant-zero reading | no |
//! | `corrupt.bit1@case37`, `corrupt.bit2@case7` | one flipped case, on a bit the check can see | **yes** |
//! | `corrupt.bit0@case37`, `corrupt.bit3@case37` | one flipped case, on a bit it structurally cannot | **yes** |
//!
//! Eight rows per arm, one per corner, so `corner_index` is the row key and every
//! per-corner column is a genuine per-row value rather than a `|`-joined blob.
//! The `|`-joined vector is *also* carried on every row as `influence_by_corner`,
//! which is how `p-167.csv` reports it, so the two artefacts can be diffed
//! directly.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record` and every panic starts `VOID: `.
//!
//! - **A corrupted table must produce unequal within-class influences**, the
//!   registration's own control. `common::boolean::corrupt` flips one case, and
//!   the module's measurements name `(37, delta 2)` on bit 1 and `(7, delta 4)` on
//!   bit 2 as corruptions that split the eight influences. Both must report
//!   `influence_equal_within_class = false`, or C1 cannot detect the defect it
//!   exists for. Column `control_detected`.
//! - **The corruption must have landed on all four control arms, including the
//!   two the check is blind to.** A single-case flip moves `Inf_i` by exactly
//!   `+/- 2/256` = [`SINGLE_FLIP_SHIFT`] for *every* `i` — only the pair
//!   `(x0, x0 xor e_i)` is touched — so the blind spot is a property of the
//!   **sign pattern**, not of the corruption failing to apply. Every corner of
//!   every control arm must show `influence_shift_abs` equal to `2/256`. Without
//!   this, "still equal" could mean "`corrupt` did nothing", and bits 0 and 3
//!   would look like a working check rather than a measured blind spot. Bit 0 is
//!   parity, so every neighbour always differed and all eight influences fall
//!   together; bit 3 is the constant zero, so none ever differed and all eight
//!   rise together.
//! - **The primary reading must be non-degenerate.** At least one of its bits must
//!   have a strictly positive influence, or C1's equality is an equality of zeros.
//! - **The primary reading must be octahedrally invariant**, `octahedral_violations
//!   = 0`. Without it, an unequal influence would localise nothing and C1 would not
//!   be a check on `validate_table()` but a coincidence.
//! - **The edge-mask aggregate must be the cube's vertex degree**, `3.0` on all
//!   eight corners — a closed-form answer through the same instrument, so a
//!   transform or an accumulation error is visible against arithmetic rather than
//!   against a transcribed decimal.
//! - **`influence` and `influence_combinatorial` must agree** on every corner of
//!   every arm. Two independent computations of one rational number.
//! - **`influence_equal_within_class` must be able to come out both ways over the
//!   census**, or it is a constant column and not a measurement (M-44). The
//!   primary reading supplies `true`, the labelled edge bits supply `false`.
//! - **The refinement priority must have rank variance**, at least
//!   [`MIN_DISTINCT_PRIORITIES`] distinct values, and `rank_correlation` must
//!   return `+1` against it and `-1` against its negation. Only then is a zero
//!   correlation attributable to the influence side.
//! - **Every scanned grid must contain a cut cell**, or that field contributes a
//!   zero priority vector and the pooled ranking is over fewer fields than the
//!   header claims.
//! - **The primary reading's per-bit total influence must reproduce
//!   `p-167.csv`'s committed values** `8.0`, `4.0`, `3.5`, `0.0`
//!   ([`P167_PRIMARY_TOTAL_INFLUENCE`]). Two wave-27 CSVs measure the same
//!   quantity through the same module; if they disagree, one of them is wrong and
//!   silence is the worst outcome.
//!
//! # SHARE, recomputed before the numbers
//!
//! Registered: **`SHARE: none`**, and that is correct rather than an omission.
//! C1 is a *check* on a table already gated by `validate_table()`, so holding
//! moves nothing — it adds a cheap independent instrument, and the value of a
//! passing check is that it could have failed. C2 is registered expecting a null,
//! and a null cannot move a stage. The one branch that would have carried a SHARE
//! is C1 *failing*, which would have been a table defect and the whole ticket by
//! itself; a registration cannot promise a share for an outcome it hopes not to
//! see. So: nothing to land, and the finding says so.
//!
//! # Determinism and cost
//!
//! **No clock is read.** The registration names no ratio and no threshold on time,
//! so every column is an exact rational or an integer count — which M-280's
//! `amd-pstate-epp` scatter is the standing reason to prefer. No RNG either: the
//! influence census is exhaustive over all 256 cases and the field scans are the
//! crate's own deterministic grids at [`RESOLUTIONS`] `= 33³` and `65³`, the two
//! the authoring contract names. Two resolutions rather than one so the pooled
//! ranking cannot be one grid's artefact, and so the refinement priority's
//! residual asymmetry can be seen not to shrink with refinement — a single grid
//! could not distinguish a boundary effect from the fields' own.
//!
//! Whole run: **4.2 s** wall including the 2,359,296 cells of the sixteen scans,
//! against a two-minute budget. 200 rows, 45 columns.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::table::{EDGE_CORNERS, is_inside};
use isomesh::{Sdf, Shape3, for_each_reference_field};

use crate::common::beta::rank_correlation;
use crate::common::boolean::{
    Bool8, NEGLIGIBLE, corner_symmetry_classes, corrupt, shipped_centroid_counts,
    shipped_edge_masks, shipped_triangle_counts,
};
use crate::common::poly::{octahedral_relabellings, relabel};

// ─── clause constants ───────────────────────────────────────────────────────

/// A rank correlation at least this large, in absolute value, counts as
/// informative for C2.
///
/// One half. Below it the per-corner ordering agrees with the refinement
/// priority on fewer than three quarters of concordant pairs, which is not a
/// signal a caller could spend a sample on. The number does not decide this
/// clause's verdict in practice — the primary reading's correlation is exactly
/// zero, for the structural reason the header sets out — but a threshold read off
/// after the fact is not a threshold.
const INFORMATIVE_CORRELATION: f64 = 0.5;

/// The corner-sample perturbation the refinement priority is the response to,
/// as a fraction of the cell size.
///
/// One part in a thousand. The quantity is a *linearised* displacement, so it
/// wants a perturbation small enough that the linearisation is the answer, and
/// the clamp at one edge length is what makes the near-tangential crossing
/// finite rather than what makes it small. A fraction of `cell_size` rather than
/// of the cell's own value scale so the clamp engages at the same field-value
/// threshold everywhere on a grid instead of drifting with local magnitude.
const PERTURBATION: f64 = 1e-3;

/// The grids the refinement priority and the case histogram are accumulated
/// over.
///
/// `33³` and `65³` — the two the authoring contract names, and `n` samples span
/// `n - 1` cells, so `32³` and `64³` cells per field. Two rather than one so the
/// pooled priority is not one grid's artefact, and so its residual corner
/// asymmetry can be seen not to shrink with refinement: measured `1.52%` of the
/// vector's own mean at `33³` and `2.68%` at `65³`, which rules out the
/// boundary-layer reading a single grid would have left open.
const RESOLUTIONS: [u32; 2] = [33, 65];

/// How many distinct values the pooled refinement priority must take before a
/// correlation against it means anything.
///
/// Three: Spearman on eight points with fewer than three distinct values is a
/// comparison against something very close to a constant, and
/// `rank_correlation`'s own floor is three pairs.
const MIN_DISTINCT_PRIORITIES: usize = 3;

/// The exact amount a single-case flip moves every one of the eight influences.
///
/// `2/256`. Flipping the output at one input `x0` touches only the pair
/// `(x0, x0 xor e_i)` in the flip count for variable `i`, and that pair
/// contributes `2` when unequal and `0` when equal — so the change is `+/- 2/256`
/// for every `i`, always. The *sign pattern* is what decides whether the eight
/// influences split, and that is the whole content of the blind spot at bits 0
/// and 3.
const SINGLE_FLIP_SHIFT: f64 = 2.0 / 256.0;

/// Every corner of a cube has three incident edges, so the edge-mask aggregate
/// influence is `3.0` on every corner.
///
/// The cube graph is 3-regular. An edge bit is `x_a xor x_b`, whose influence is
/// `1` on `a` and `b` and `0` elsewhere; summing over the twelve bits gives each
/// corner its degree. A closed form through the same instrument, which is what
/// makes it a calibration.
const CUBE_VERTEX_DEGREE: f64 = 3.0;

/// The primary reading's per-bit total influence, as committed in
/// `docs/experiments/p-167.csv`.
///
/// `8.0` for bit 0 (parity: every flip always flips the output), `4.0` for bit 1,
/// `3.5` for bit 2, `0.0` for bit 3 (the constant-zero witness). Read from that
/// CSV's `total_influence` column, not from a summary. Both rows measure the same
/// quantity through the same module, so a disagreement means one of the two
/// artefacts is wrong.
const P167_PRIMARY_TOTAL_INFLUENCE: [f64; 4] = [8.0, 4.0, 3.5, 0.0];

/// The absolute gap allowed between an arm's uniform-measure total influence and
/// its field-measure counterpart before the two are called different numbers.
///
/// `0.8`, which is a tenth of the maximum total influence an 8-variable function
/// can have. Absolute rather than relative so a constant arm — both numbers zero
/// — is compared without a division.
const SENSITIVITY_TOLERANCE: f64 = 0.8;

/// One deliberately corrupted table: which output bit is analysed, which case is
/// altered, and by what `xor` delta.
#[derive(Clone, Copy, Debug)]
struct Corruption {
    bit: u32,
    case: usize,
    delta: u32,
}

/// Corruptions of the primary reading whose defect the influence check must see.
///
/// `common::boolean::corrupt`'s own measurements over all 256 single-case flips
/// name these two: bit 1 splits on 236 of 256 cases and bit 2 on 204, and `37`
/// and `7` are cases with a mixed neighbourhood. The delta is the bit's own mask
/// because the alteration is an `xor` and must move exactly one bit.
const CORRUPT_DETECTING: [Corruption; 2] = [
    Corruption {
        bit: 1,
        case: 37,
        delta: 2,
    },
    Corruption {
        bit: 2,
        case: 7,
        delta: 4,
    },
];

/// Corruptions the influence check is structurally blind to, run on purpose.
///
/// Bit 0 of the triangle count is exactly `Bool8::parity` and bit 3 is the
/// constant zero; for parity every neighbour of `x0` already differed and for a
/// constant none did, so in both cases all eight influences move the same way and
/// equality survives. The module measures 0 of 256 flips detected on each. A
/// bench that corrupted bit 0 and concluded "C1 cannot detect a defect" would
/// have measured the wrong bit, so the blind spot is recorded rather than
/// avoided — with `influence_shift_abs` proving the corruption landed.
const CORRUPT_BLIND: [Corruption; 2] = [
    Corruption {
        bit: 0,
        case: 37,
        delta: 1,
    },
    Corruption {
        bit: 3,
        case: 37,
        delta: 8,
    },
];

// ─── the cube's own layout ──────────────────────────────────────────────────

/// Corner `i`'s offset in cells: `[(i & 1), ((i >> 1) & 1), ((i >> 2) & 1)]`.
///
/// The crate's own layout (`crates/isomesh/src/cube.rs:149-153`, pinned by the
/// test at `cube.rs:315-320`), restated here because `cube::corner_offset` is
/// `pub(crate)` and out of a bench's reach. It is the same numbering
/// `common::boolean` documents for its variables and the same one
/// `octahedral_relabellings` rebuilds indices with, which is what lets an
/// influence, a case index and a relabelling all mean the corner they say.
const fn corner_offset(corner: usize) -> [u32; 3] {
    [
        (corner & 1) as u32,
        ((corner >> 1) & 1) as u32,
        ((corner >> 2) & 1) as u32,
    ]
}

/// Relabel a case index by a corner permutation: the new sign of corner `i` is
/// the old sign of corner `perm[i]`.
///
/// `relabel` is `common::poly`'s, so the convention is the one its 48
/// relabellings were generated with and the two cannot drift apart.
fn relabel_case(perm: &[u8; 8], case: usize) -> usize {
    let signs: [u8; 8] = std::array::from_fn(|i| ((case >> i) & 1) as u8);
    let moved = relabel(perm, &signs);
    moved
        .iter()
        .enumerate()
        .fold(0usize, |acc, (i, &b)| acc | (usize::from(b) << i))
}

/// The permutation of the twelve edge labels induced by a corner permutation.
///
/// Bit `e` of `mask(x')` is set when `s[perm[a]] != s[perm[b]]` for
/// `{a, b} = EDGE_CORNERS[e]`, and `{perm[a], perm[b]}` is again an edge of the
/// cube because a cube symmetry maps edges to edges. So `mask(x')` reads bit
/// `sigma[e]` of `mask(x)`, and this returns `sigma`.
///
/// # Panics
///
/// If some corner pair's image is not an edge, or if the twelve images are not a
/// permutation — either would mean the relabellings are not the cube group and
/// every equivariance count below would be about the wrong set.
fn induced_edge_permutation(perm: &[u8; 8]) -> [u8; 12] {
    let mut out = [u8::MAX; 12];
    for (e, corners) in EDGE_CORNERS.iter().enumerate() {
        let pa = perm[usize::from(corners[0])];
        let pb = perm[usize::from(corners[1])];
        let (lo, hi) = if pa < pb { (pa, pb) } else { (pb, pa) };
        let target = EDGE_CORNERS
            .iter()
            .position(|c| c[0] == lo && c[1] == hi)
            .expect("a cube symmetry maps an edge of the cube to an edge of the cube");
        out[e] = target as u8;
    }
    let mut seen = [false; 12];
    for &e in &out {
        assert!(
            !seen[usize::from(e)],
            "the induced edge map is not a permutation of the twelve edges"
        );
        seen[usize::from(e)] = true;
    }
    out
}

/// Apply an edge permutation to a twelve-bit mask: bit `e` of the result is bit
/// `sigma[e]` of the input.
fn permute_mask(mask: u32, sigma: &[u8; 12]) -> u32 {
    sigma
        .iter()
        .enumerate()
        .fold(0u32, |acc, (e, &s)| acc | (((mask >> s) & 1) << e))
}

/// How the three readings behave under the 48 corner relabellings.
#[derive(Clone, Copy, Debug)]
struct Invariance {
    /// `(pi, x)` pairs where the primary reading's output moved. Predicted 0.
    triangle_labelled: usize,
    /// `(pi, x)` pairs where the edge mask moved as a labelled function.
    /// Predicted large — this is what makes a per-bit influence non-invariant.
    edge_labelled: usize,
    /// `(pi, x)` pairs where the edge mask moved by anything other than the
    /// induced edge permutation. Predicted 0 — the failure above is exactly the
    /// relabelling.
    edge_equivariant: usize,
    /// `(pi, x)` pairs where the centroid reading's output moved. Predicted 0; a
    /// constant is invariant under everything.
    centroid_labelled: usize,
}

/// Push all three readings through all 48 corner relabellings.
fn invariance(triangles: &[u32; 256], masks: &[u32; 256], centroids: &[u32; 256]) -> Invariance {
    let perms = octahedral_relabellings();
    let mut out = Invariance {
        triangle_labelled: 0,
        edge_labelled: 0,
        edge_equivariant: 0,
        centroid_labelled: 0,
    };
    for perm in &perms {
        let sigma = induced_edge_permutation(perm);
        for case in 0..256usize {
            let moved = relabel_case(perm, case);
            if triangles[moved] != triangles[case] {
                out.triangle_labelled += 1;
            }
            if centroids[moved] != centroids[case] {
                out.centroid_labelled += 1;
            }
            if masks[moved] != masks[case] {
                out.edge_labelled += 1;
            }
            if masks[moved] != permute_mask(masks[case], &sigma) {
                out.edge_equivariant += 1;
            }
        }
    }
    out
}

// ─── the readings ───────────────────────────────────────────────────────────

/// One reading of the case table's integer output, and how many bits it is
/// analysed in.
#[derive(Clone, Debug)]
struct Reading {
    name: &'static str,
    values: [u32; 256],
    bits: u32,
}

impl Reading {
    /// P-167's bit rule, reused verbatim so the two CSVs analyse the same
    /// functions: `bit_length(max)` bits carry the output, and one more is the
    /// constant-zero witness whose measured degeneracy proves the reading stops
    /// where the arithmetic says it does.
    fn new(name: &'static str, values: [u32; 256]) -> Self {
        let max = values.iter().fold(0u32, |a, &b| a.max(b));
        let bits = 32 - max.leading_zeros() + 1;
        Self { name, values, bits }
    }
}

// ─── the field side: case histogram and refinement priority ─────────────────

/// One scanned grid: which cases its cells landed in, and how much each corner's
/// sample moves the extracted geometry.
#[derive(Clone, Debug)]
struct Scan {
    label: String,
    /// Per-corner refinement priority, normalised to sum to one so a field with
    /// more cut cells does not simply outvote the others.
    priority: [f64; 8],
    cells: u64,
    cut_cells: u64,
}

/// Sample one reference field on one grid, accumulate the case histogram, and
/// return the per-corner refinement priority.
///
/// The priority is the linearised displacement of the Marching Cubes vertex under
/// a [`PERTURBATION`]-of-a-cell perturbation of that corner's sample, clamped to
/// the edge the vertex cannot leave. On a cut edge `v_a - v_b` is strictly
/// non-zero — `is_inside` is `v < 0`, so one endpoint is strictly negative and the
/// other is not — so the expression is total; the clamp is part of the definition
/// and keeps one near-tangential crossing from owning the sum.
///
/// # Panics
///
/// If the grid contains no cut cell, which would make this field's contribution a
/// zero vector and the pooled ranking a ranking over fewer fields than claimed.
fn scan<F>(name: &str, field: &F, samples: u32, histogram: &mut [u64; 256]) -> Scan
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, h) = common::grid::<f64, _>(field, samples);
    let n = samples as usize;
    let mut values = Vec::with_capacity(n * n * n);
    for z in 0..n {
        let pz = origin[2] + h * z as f64;
        for y in 0..n {
            let py = origin[1] + h * y as f64;
            for x in 0..n {
                values.push(field.sample([origin[0] + h * x as f64, py, pz]));
            }
        }
    }

    let d = PERTURBATION * h;
    let mut priority = [0.0f64; 8];
    let mut cut_cells = 0u64;
    let mut cells = 0u64;
    for cz in 0..samples - 1 {
        for cy in 0..samples - 1 {
            for cx in 0..samples - 1 {
                let corner: [f64; 8] = std::array::from_fn(|i| {
                    let o = corner_offset(i);
                    values[shape.linearize([cx + o[0], cy + o[1], cz + o[2]]) as usize]
                });

                let mut case = 0u8;
                for (i, &v) in corner.iter().enumerate() {
                    if is_inside(v) {
                        case |= 1 << i;
                    }
                }
                histogram[usize::from(case)] += 1;
                cells += 1;

                let mut cut = false;
                for [a, b] in EDGE_CORNERS {
                    let va = corner[usize::from(a)];
                    let vb = corner[usize::from(b)];
                    if is_inside(va) == is_inside(vb) {
                        continue;
                    }
                    cut = true;
                    let gap = va - vb;
                    let response = d / (gap * gap);
                    priority[usize::from(a)] += (response * vb.abs()).min(1.0);
                    priority[usize::from(b)] += (response * va.abs()).min(1.0);
                }
                if cut {
                    cut_cells += 1;
                }
            }
        }
    }

    let total: f64 = priority.iter().sum();
    assert!(
        cut_cells > 0 && total > 0.0,
        "VOID: {name} at {samples}^3 has {cut_cells} cut cells and a priority sum of {total}, so \
         it contributes nothing to the refinement ranking and the pooled priority is over fewer \
         fields than this experiment claims"
    );
    for slot in &mut priority {
        *slot /= total;
    }

    Scan {
        label: format!("{name}@{samples}"),
        priority,
        cells,
        cut_cells,
    }
}

// ─── one measured arm ───────────────────────────────────────────────────────

/// Everything this experiment measures about one Boolean function of the eight
/// corner signs.
#[derive(Clone, Debug)]
struct Arm {
    /// `triangle_counts.bit1`, and the same key `p-167.csv` uses for
    /// `output_bit`, so the two artefacts join.
    label: String,
    reading: &'static str,
    bit_label: String,
    control: Option<Corruption>,
    /// Spectral influence, `sum_{S ni i} fhat(S)^2`.
    influences: [f64; 8],
    /// The same number counted directly over the 256 inputs.
    combinatorial: [f64; 8],
    /// `sum_S |S| fhat(S)^2`, from the module's own spectral form.
    total_influence: f64,
    /// The same total, re-added from the eight per-corner numbers.
    total_from_corners: f64,
    /// Influence under the measured case histogram of the reference fields
    /// instead of under the uniform measure.
    empirical: [f64; 8],
    empirical_total: f64,
    /// Distance from the uncorrupted same-bit influence, per corner. `None` off
    /// the control arms, where an uncorrupted table has no distance from itself
    /// and a zero would read as a measurement; exactly [`SINGLE_FLIP_SHIFT`] on
    /// every corner of every control arm, which is what proves the corruption
    /// landed even where the equality check is blind to it.
    shift_from_shipped: Option<[f64; 8]>,
    role: &'static str,
    correlation: f64,
    correlation_per_scan: Vec<f64>,
}

impl Arm {
    /// Analyse one output bit of one reading.
    ///
    /// `control` is `Some((corruption, the uncorrupted same-bit influence))` on a
    /// control arm and `None` otherwise, so the corruption and the baseline it is
    /// measured against cannot be supplied one without the other — a shipped arm
    /// has no distance from itself, and passing a zero baseline for one would
    /// silently report the influence as if it were a shift.
    fn bit(
        reading: &'static str,
        values: &[u32; 256],
        bit: u32,
        control: Option<(Corruption, [f64; 8])>,
        histogram: &[u64; 256],
        cells: u64,
    ) -> Self {
        let f = Bool8::from_values(values, bit);
        let influences: [f64; 8] = std::array::from_fn(|i| f.influence(i));
        let combinatorial: [f64; 8] = std::array::from_fn(|i| f.influence_combinatorial(i));
        let empirical = empirical_bit(&f, histogram, cells);
        let label = match control {
            Some((c, _)) => format!("corrupt.bit{}@case{}", c.bit, c.case),
            None => format!("{reading}.bit{bit}"),
        };
        Self {
            label,
            reading,
            bit_label: format!("bit{bit}"),
            control: control.map(|(c, _)| c),
            influences,
            combinatorial,
            total_influence: f.total_influence(),
            total_from_corners: influences.iter().sum(),
            empirical,
            empirical_total: empirical.iter().sum(),
            shift_from_shipped: control
                .map(|(_, base)| std::array::from_fn(|i| (influences[i] - base[i]).abs())),
            role: role(&influences),
            correlation: 0.0,
            correlation_per_scan: Vec::new(),
        }
    }

    /// The reading's **vector-valued** influence: per corner, the expected number
    /// of output bits that flip when that corner's sign flips. This is the level
    /// at which the cube group acts on an equivariant reading, so it is the
    /// honest form of "only the multiset over bits is invariant".
    fn aggregate(reading: &'static str, r: &Reading, histogram: &[u64; 256], cells: u64) -> Self {
        let mut influences = [0.0f64; 8];
        let mut combinatorial = [0.0f64; 8];
        let mut empirical = [0.0f64; 8];
        let mut total_influence = 0.0f64;
        for bit in 0..r.bits {
            let f = Bool8::from_values(&r.values, bit);
            for (i, slot) in influences.iter_mut().enumerate() {
                *slot += f.influence(i);
            }
            for (i, slot) in combinatorial.iter_mut().enumerate() {
                *slot += f.influence_combinatorial(i);
            }
            let e = empirical_bit(&f, histogram, cells);
            for (slot, v) in empirical.iter_mut().zip(e.iter()) {
                *slot += v;
            }
            total_influence += f.total_influence();
        }
        Self {
            label: format!("{reading}.aggregate"),
            reading,
            bit_label: String::from("aggregate"),
            control: None,
            influences,
            combinatorial,
            total_influence,
            total_from_corners: influences.iter().sum(),
            empirical,
            empirical_total: empirical.iter().sum(),
            shift_from_shipped: None,
            role: role(&influences),
            correlation: 0.0,
            correlation_per_scan: Vec::new(),
        }
    }

    /// Is this arm's influence equal, within [`NEGLIGIBLE`], across every pair of
    /// corners sharing a symmetry class? With one class — which is what
    /// `corner_symmetry_classes()` generates — this is "all eight equal".
    fn equal_within_class(&self, classes: &[u8; 8]) -> bool {
        (0..8).all(|i| {
            (0..8).all(|j| {
                classes[i] != classes[j]
                    || (self.influences[i] - self.influences[j]).abs() <= NEGLIGIBLE
            })
        })
    }

    /// A constant arm's equality is an equality of zeros, which is not an
    /// equality (M-44), so it does not carry a clause verdict.
    fn is_degenerate(&self) -> bool {
        self.role == "constant"
    }

    fn spread(&self) -> f64 {
        let lo = self.influences.iter().copied().fold(f64::MAX, f64::min);
        let hi = self.influences.iter().copied().fold(f64::MIN, f64::max);
        hi - lo
    }

    /// Do the spectral and combinatorial influences agree on every corner? The
    /// module asserts it internally; this puts it in the artefact.
    fn agrees_combinatorially(&self) -> bool {
        self.influences
            .iter()
            .zip(self.combinatorial.iter())
            .all(|(a, b)| (a - b).abs() <= NEGLIGIBLE)
    }

    /// Does the uniform-measure average sensitivity match the one a refinement
    /// heuristic would face on the reference fields' own cells?
    fn sensitivity_agrees(&self) -> bool {
        (self.total_influence - self.empirical_total).abs() <= SENSITIVITY_TOLERANCE
    }
}

/// Influence under the measured case histogram rather than the uniform measure.
///
/// The same flip count, re-weighted: `Pr[f(x) != f(x + e_i)]` when `x` is drawn
/// from the cells of the eight reference fields instead of uniformly from
/// `{0,1}^8`. A `xor` of corner signs is pivotal at every input, so its influence
/// is measure-independent and the two numbers must agree exactly; anything else
/// is free to disagree, and the informative bits of the primary reading are
/// predicted to.
fn empirical_bit(f: &Bool8, histogram: &[u64; 256], cells: u64) -> [f64; 8] {
    std::array::from_fn(|i| {
        let bit = 1usize << i;
        let flips: u64 = (0..256usize)
            .filter(|&x| f.0[x] != f.0[x ^ bit])
            .map(|x| histogram[x])
            .sum();
        flips as f64 / cells as f64
    })
}

/// What the influence vector's shape is, derived from the measurement rather than
/// transcribed from the module's documentation.
///
/// `parity` is a theorem and not a guess: an influence of `1` on every variable
/// means every flip always flips the output, which on `{0,1}^n` is parity or its
/// complement. `edge_pair` is the `(1, 1, 0, 0, 0, 0, 0, 0)` shape of a
/// two-corner `xor` up to relabelling — equivariant but not invariant, and the
/// shape C1 must not be read against. `orbit_equal` is all eight equal without
/// being either. `split` is the shape a corrupted table is supposed to produce.
fn role(influences: &[f64; 8]) -> &'static str {
    let equal = influences
        .iter()
        .all(|v| (v - influences[0]).abs() <= NEGLIGIBLE);
    let zeros = influences.iter().filter(|v| v.abs() <= NEGLIGIBLE).count();
    let ones = influences
        .iter()
        .filter(|v| (*v - 1.0).abs() <= NEGLIGIBLE)
        .count();
    if zeros == 8 {
        "constant"
    } else if equal && ones == 8 {
        "parity"
    } else if equal {
        "orbit_equal"
    } else if ones == 2 && zeros == 6 {
        "edge_pair"
    } else {
        "split"
    }
}

// ─── formatting ─────────────────────────────────────────────────────────────

/// A vector as `a|b|c`, because `Run::record` refuses a comma: the writer does
/// not quote, so a comma would shift every later column.
fn joined(values: &[f64]) -> String {
    values
        .iter()
        .map(|v| format!("{v:.9}"))
        .collect::<Vec<String>>()
        .join("|")
}

/// Labels as `a|b|c`, for the per-scan correlation column's key.
fn joined_labels(labels: &[String]) -> String {
    labels.join("|")
}

// ─── the run ────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-169");

    common::experiment::run(prereg, |run| {
        // ── the symmetry classes, generated ────────────────────────────────
        let (classes, class_count) = corner_symmetry_classes();

        // ── the three readings ─────────────────────────────────────────────
        let triangles = Reading::new("triangle_counts", shipped_triangle_counts());
        let masks = Reading::new("edge_masks", shipped_edge_masks());
        let centroids = Reading::new("centroid_counts", shipped_centroid_counts());

        // ── vacuity control 1: the primary reading is octahedrally invariant,
        //    and the edge masks are equivariant without being invariant ──────
        let inv = invariance(&triangles.values, &masks.values, &centroids.values);
        assert!(
            inv.triangle_labelled == 0,
            "VOID: the primary reading moved on {} of the 48*256 corner relabellings, so it is \
             not octahedrally invariant, an unequal influence would localise nothing, and C1 is \
             not a check on validate_table() but a coincidence",
            inv.triangle_labelled
        );
        assert!(
            inv.edge_equivariant == 0,
            "VOID: the edge mask failed to follow the induced edge permutation on {} of the \
             48*256 relabellings, so its per-bit influence being unequal is not explained by the \
             relabelling and this experiment's reason for excluding it from C1 is wrong",
            inv.edge_equivariant
        );
        assert!(
            inv.edge_labelled > 0,
            "VOID: the edge mask is invariant as a labelled function, so the labelled/multiset \
             distinction this experiment is built on does not exist and C1 could have been read \
             against it after all"
        );

        // ── the field side ─────────────────────────────────────────────────
        let mut histogram = [0u64; 256];
        let mut scans: Vec<Scan> = Vec::new();
        for_each_reference_field!(f64, |name, field| {
            for samples in RESOLUTIONS {
                scans.push(scan(name, &field, samples, &mut histogram));
            }
        });
        let scanned_cells: u64 = scans.iter().map(|s| s.cells).sum();
        let cut_cells: u64 = scans.iter().map(|s| s.cut_cells).sum();
        assert!(
            scanned_cells == histogram.iter().sum::<u64>(),
            "VOID: {scanned_cells} cells were walked but the case histogram holds {}, so the \
             empirical influence is normalised against the wrong denominator",
            histogram.iter().sum::<u64>()
        );

        // Each scan's priority already sums to one, so the pooled vector sums to
        // the number of scans; normalise it back so the column is readable as a
        // share and the spread is a share of a share.
        let mut priority = [0.0f64; 8];
        for s in &scans {
            for (slot, v) in priority.iter_mut().zip(s.priority.iter()) {
                *slot += v;
            }
        }
        for slot in &mut priority {
            *slot /= scans.len() as f64;
        }
        let priority_lo = priority.iter().copied().fold(f64::MAX, f64::min);
        let priority_hi = priority.iter().copied().fold(f64::MIN, f64::max);
        // The mean is 1/8 by construction, so the spread is directly comparable
        // across resolutions.
        let priority_spread = (priority_hi - priority_lo) * 8.0;

        // ── vacuity control 2: the priority side could have carried a ranking ─
        let mut distinct: Vec<f64> = priority.to_vec();
        distinct.sort_by(f64::total_cmp);
        distinct.dedup_by(|a, b| (*a - *b).abs() <= NEGLIGIBLE);
        assert!(
            distinct.len() >= MIN_DISTINCT_PRIORITIES,
            "VOID: the pooled refinement priority takes only {} distinct values, so a zero \
             correlation would be a property of the priority rather than of the influences and \
             C2's verdict would be unattributable",
            distinct.len()
        );
        let control_self = rank_correlation(&priority, &priority);
        let negated: [f64; 8] = std::array::from_fn(|i| -priority[i]);
        let control_negated = rank_correlation(&priority, &negated);
        assert!(
            (control_self - 1.0).abs() <= NEGLIGIBLE && (control_negated + 1.0).abs() <= NEGLIGIBLE,
            "VOID: rank_correlation returns {control_self} against the priority itself and \
             {control_negated} against its negation, not +1 and -1, so the instrument cannot \
             report a correlation over this very sample and a zero means nothing"
        );

        // ── the arms ───────────────────────────────────────────────────────
        let mut arms: Vec<Arm> = Vec::new();
        for r in [&triangles, &masks, &centroids] {
            for bit in 0..r.bits {
                arms.push(Arm::bit(
                    r.name,
                    &r.values,
                    bit,
                    None,
                    &histogram,
                    scanned_cells,
                ));
            }
            arms.push(Arm::aggregate(r.name, r, &histogram, scanned_cells));
        }

        // ── vacuity control 3: the primary reading is not all constants ─────
        let primary_positive = arms
            .iter()
            .filter(|a| a.reading == triangles.name && !a.is_degenerate())
            .count();
        assert!(
            primary_positive > 0,
            "VOID: every arm of the primary reading has all eight influences zero, so C1's \
             equality is an equality of zeros and its verdict is a property of the fixture (M-44)"
        );

        // ── vacuity control 4: the primary reading reproduces p-167.csv ─────
        for (bit, &expected) in P167_PRIMARY_TOTAL_INFLUENCE.iter().enumerate() {
            let measured = arms
                .iter()
                .find(|a| a.label == format!("{}.bit{bit}", triangles.name))
                .map_or(f64::NAN, |a| a.total_influence);
            assert!(
                (measured - expected).abs() <= NEGLIGIBLE,
                "VOID: {}.bit{bit} has total influence {measured} but p-167.csv committed \
                 {expected} for the same quantity through the same module, so one of the two \
                 artefacts is wrong",
                triangles.name
            );
        }

        // ── vacuity control 5: the edge-mask aggregate is the cube's degree ──
        let edge_aggregate = arms
            .iter()
            .find(|a| a.label == format!("{}.aggregate", masks.name))
            .map_or([f64::NAN; 8], |a| a.influences);
        assert!(
            edge_aggregate
                .iter()
                .all(|v| (v - CUBE_VERTEX_DEGREE).abs() <= NEGLIGIBLE),
            "VOID: the edge-mask aggregate influence is {} and not {CUBE_VERTEX_DEGREE} on every \
             corner, which is the cube graph's vertex degree in closed form — so either the \
             transform or this accumulation is wrong",
            joined(&edge_aggregate)
        );

        // ── the corrupt control arms ───────────────────────────────────────
        let shipped_bit = |bit: u32| {
            arms.iter()
                .find(|a| a.label == format!("{}.bit{bit}", triangles.name))
                .map_or([f64::NAN; 8], |a| a.influences)
        };
        let mut controls: Vec<Arm> = Vec::new();
        for c in CORRUPT_DETECTING.into_iter().chain(CORRUPT_BLIND) {
            let table = corrupt(&triangles.values, c.case, c.delta);
            let baseline = shipped_bit(c.bit);
            controls.push(Arm::bit(
                triangles.name,
                &table,
                c.bit,
                Some((c, baseline)),
                &histogram,
                scanned_cells,
            ));
        }

        // ── vacuity control 6: the registration's own control ──────────────
        for c in CORRUPT_DETECTING {
            let arm = controls
                .iter()
                .find(|a| {
                    a.control
                        .is_some_and(|k| k.bit == c.bit && k.case == c.case)
                })
                .expect("every corruption produced an arm");
            assert!(
                !arm.equal_within_class(&classes),
                "VOID: flipping case {} of the primary reading by {} left all eight influences of \
                 bit {} equal, so C1 cannot detect the defect it is designed for and a held C1 \
                 says nothing about the table",
                c.case,
                c.delta,
                c.bit
            );
        }

        // ── vacuity control 7: every corruption landed, blind ones included ──
        for arm in &controls {
            assert!(
                arm.shift_from_shipped.as_ref().is_some_and(|shift| {
                    shift
                        .iter()
                        .all(|v| (v - SINGLE_FLIP_SHIFT).abs() <= NEGLIGIBLE)
                }),
                "VOID: {} moved the eight influences by {} rather than by exactly \
                 {SINGLE_FLIP_SHIFT} on every corner — a single-case flip touches exactly one \
                 pair per variable, so either the corruption did not land or the influence is \
                 not counting what it claims, and an unsplit control would then look like a \
                 blind spot instead of a no-op",
                arm.label,
                joined(arm.shift_from_shipped.as_ref().unwrap_or(&arm.influences))
            );
        }
        arms.extend(controls);

        // ── vacuity control 8: the two independent influences agree ────────
        for arm in &arms {
            assert!(
                arm.agrees_combinatorially(),
                "VOID: {} has spectral influences {} and flip-counted influences {}, so the \
                 Walsh-Hadamard transform and the direct count disagree and neither number is a \
                 measurement",
                arm.label,
                joined(&arm.influences),
                joined(&arm.combinatorial)
            );
            assert!(
                (arm.total_influence - arm.total_from_corners).abs() <= NEGLIGIBLE,
                "VOID: {} reports total influence {} spectrally but {} when the eight per-corner \
                 numbers are re-added",
                arm.label,
                arm.total_influence,
                arm.total_from_corners
            );
        }

        // ── the correlations ───────────────────────────────────────────────
        for arm in &mut arms {
            arm.correlation = rank_correlation(&arm.influences, &priority);
            arm.correlation_per_scan = scans
                .iter()
                .map(|s| rank_correlation(&arm.influences, &s.priority))
                .collect();
        }
        let scan_labels = joined_labels(
            &scans
                .iter()
                .map(|s| s.label.clone())
                .collect::<Vec<String>>(),
        );

        // ── vacuity control 9: the equality column is a measurement ────────
        let equal_arms = arms
            .iter()
            .filter(|a| a.equal_within_class(&classes))
            .count();
        assert!(
            equal_arms > 0 && equal_arms < arms.len(),
            "VOID: influence_equal_within_class is {} on all {} arms, so it is a constant column \
             over this census and could not have come out the other way (M-44)",
            equal_arms > 0,
            arms.len()
        );

        // ── the two global verdicts ────────────────────────────────────────
        //
        // Both clauses are claims about the shipped table read on the reading
        // they are about, so both are global and both go on every row; the
        // per-arm measurements live in the registered per-row columns. The
        // primary reading's degenerate bits are excluded from either verdict --
        // an equality of zeros is not an equality, and a constant's correlation
        // is not a correlation.
        let primary: Vec<&Arm> = arms
            .iter()
            .filter(|a| a.reading == triangles.name && a.control.is_none() && !a.is_degenerate())
            .collect();
        let c1 = primary.iter().all(|a| a.equal_within_class(&classes));
        let c2_correlation = primary
            .iter()
            .map(|a| a.correlation.abs())
            .fold(0.0f64, f64::max);
        let c2 = c2_correlation >= INFORMATIVE_CORRELATION;

        // ── what came out ──────────────────────────────────────────────────
        println!(
            "the cube group's orbits on the eight corners: {class_count} class(es), assignment \
             {classes:?}\nso C1's 'equal within each octahedral symmetry class' is its strongest \
             reading: all eight influences must be equal"
        );
        println!(
            "over 48 relabellings x 256 cases -- triangle_counts moved {} times (invariant), \
             edge_masks moved {} times as a labelled function but {} times against its induced \
             edge permutation (equivariant, not invariant), centroid_counts moved {} times",
            inv.triangle_labelled, inv.edge_labelled, inv.edge_equivariant, inv.centroid_labelled
        );
        println!(
            "{} grids over {} reference fields at {RESOLUTIONS:?}: {scanned_cells} cells, \
             {cut_cells} of them cut",
            scans.len(),
            scans.len() / RESOLUTIONS.len()
        );
        println!(
            "pooled refinement priority {} -- spread {priority_spread:.6} of its own mean, \
             {} distinct values; rank_correlation against itself {control_self:.6} and against \
             its negation {control_negated:.6}",
            joined(&priority),
            distinct.len()
        );
        for samples in RESOLUTIONS {
            let suffix = format!("@{samples}");
            let mut per: [f64; 8] = [0.0; 8];
            let mut count = 0usize;
            for s in scans.iter().filter(|s| s.label.ends_with(&suffix)) {
                for (slot, v) in per.iter_mut().zip(s.priority.iter()) {
                    *slot += v;
                }
                count += 1;
            }
            for slot in &mut per {
                *slot /= count as f64;
            }
            let lo = per.iter().copied().fold(f64::MAX, f64::min);
            let hi = per.iter().copied().fold(f64::MIN, f64::max);
            println!("  priority spread at {samples}^3: {:.6}", (hi - lo) * 8.0);
        }
        println!();
        println!(
            "{:<26} {:<12} {:>9} {:>5} {:>8} {:>9}  influence_by_corner",
            "arm", "role", "total", "eq", "corr", "empirical"
        );
        for arm in &arms {
            println!(
                "{:<26} {:<12} {:>9.6} {:>5} {:>8.4} {:>9.6}  [{}]{}",
                arm.label,
                arm.role,
                arm.total_influence,
                arm.equal_within_class(&classes),
                arm.correlation,
                arm.empirical_total,
                arm.influences
                    .iter()
                    .map(|v| format!("{v:.4}"))
                    .collect::<Vec<String>>()
                    .join(" "),
                if arm.control.is_some() {
                    " control"
                } else {
                    ""
                }
            );
        }

        // ── the rows ───────────────────────────────────────────────────────
        for arm in &arms {
            let equal = arm.equal_within_class(&classes);
            let is_control = arm.control.is_some();
            for corner in 0..8usize {
                run.record(&[
                    ("corner_index", corner.to_string()),
                    ("influence", format!("{:.9}", arm.influences[corner])),
                    ("total_influence", format!("{:.9}", arm.total_influence)),
                    (
                        "symmetry_class",
                        if is_control {
                            String::from("control")
                        } else {
                            format!("class{}", classes[corner])
                        },
                    ),
                    ("influence_equal_within_class", equal.to_string()),
                    (
                        "refinement_priority_correlation",
                        format!("{:.6}", arm.correlation),
                    ),
                    ("c1_holds", c1.to_string()),
                    ("c2_holds", c2.to_string()),
                    // ── extras (M-273) ──────────────────────────────────────
                    ("output_bit", arm.label.clone()),
                    ("reading", String::from(arm.reading)),
                    ("bit_index", arm.bit_label.clone()),
                    ("role", String::from(arm.role)),
                    ("is_control", is_control.to_string()),
                    ("is_degenerate", arm.is_degenerate().to_string()),
                    (
                        "influence_combinatorial",
                        format!("{:.9}", arm.combinatorial[corner]),
                    ),
                    (
                        "influence_agrees_combinatorial",
                        arm.agrees_combinatorially().to_string(),
                    ),
                    ("influence_by_corner", joined(&arm.influences)),
                    ("influence_spread", format!("{:.9}", arm.spread())),
                    (
                        "total_influence_from_corners",
                        format!("{:.9}", arm.total_from_corners),
                    ),
                    ("symmetry_class_index", classes[corner].to_string()),
                    ("symmetry_class_count", class_count.to_string()),
                    ("refinement_priority", format!("{:.9}", priority[corner])),
                    ("refinement_priority_by_corner", joined(&priority)),
                    (
                        "refinement_priority_spread",
                        format!("{priority_spread:.9}"),
                    ),
                    ("refinement_priority_distinct", distinct.len().to_string()),
                    (
                        "refinement_priority_correlation_per_scan",
                        joined(&arm.correlation_per_scan),
                    ),
                    ("refinement_priority_scans", scan_labels.clone()),
                    ("control_correlation_self", format!("{control_self:.9}")),
                    (
                        "control_correlation_negated",
                        format!("{control_negated:.9}"),
                    ),
                    (
                        "informative_correlation",
                        format!("{INFORMATIVE_CORRELATION:.2}"),
                    ),
                    (
                        "empirical_influence",
                        format!("{:.9}", arm.empirical[corner]),
                    ),
                    (
                        "empirical_total_influence",
                        format!("{:.9}", arm.empirical_total),
                    ),
                    (
                        "uniform_minus_empirical",
                        format!("{:.9}", arm.total_influence - arm.empirical_total),
                    ),
                    (
                        "empirical_agrees_with_uniform",
                        arm.sensitivity_agrees().to_string(),
                    ),
                    (
                        "influence_shift_abs",
                        // Undefined off a control arm by construction: a
                        // shipped arm has no distance from itself, and a zero
                        // here would read as "no shift" rather than as "not a
                        // shift-shaped question".
                        arm.shift_from_shipped.as_ref().map_or_else(
                            || String::from("undefined"),
                            |s| format!("{:.9}", s[corner]),
                        ),
                    ),
                    (
                        "control_case",
                        arm.control
                            .map_or_else(|| String::from("none"), |c| c.case.to_string()),
                    ),
                    (
                        "control_delta",
                        arm.control
                            .map_or_else(|| String::from("none"), |c| c.delta.to_string()),
                    ),
                    (
                        "control_detected",
                        arm.control
                            .map_or_else(|| String::from("none"), |_| (!equal).to_string()),
                    ),
                    ("octahedral_violations", inv.triangle_labelled.to_string()),
                    (
                        "edge_mask_labelled_violations",
                        inv.edge_labelled.to_string(),
                    ),
                    (
                        "edge_mask_equivariant_violations",
                        inv.edge_equivariant.to_string(),
                    ),
                    ("scanned_cells", scanned_cells.to_string()),
                    ("scanned_cut_cells", cut_cells.to_string()),
                    (
                        "scan_resolutions",
                        RESOLUTIONS
                            .iter()
                            .map(u32::to_string)
                            .collect::<Vec<String>>()
                            .join("|"),
                    ),
                    ("perturbation", format!("{PERTURBATION:e}")),
                ]);
            }
        }

        // ── the verdicts, spelled out ──────────────────────────────────────
        println!();
        println!(
            "C1 over the {} non-degenerate arms of the primary reading: {}",
            primary.len(),
            if c1 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C2: the largest |rank correlation| between a primary arm's eight influences and the \
             refinement priority is {c2_correlation:.6}, against a threshold of \
             {INFORMATIVE_CORRELATION} -- {}",
            if c2 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "and the two verdicts are the same fact: the eight influences are equal, so the \
             sample has no rank variance and rank_correlation is exactly zero. C1 holding is \
             what makes C2 fail."
        );
        let disagreeing: Vec<&str> = arms
            .iter()
            .filter(|a| !a.sensitivity_agrees())
            .map(|a| a.label.as_str())
            .collect();
        println!(
            "average sensitivity under the fields' own case histogram differs from the \
             uniform-measure total influence by more than {SENSITIVITY_TOLERANCE} on {} of {} \
             arms: {:?}",
            disagreeing.len(),
            arms.len(),
            disagreeing
        );
        println!(
            "the corrupt controls: {} of {} split the eight influences; the other {} moved every \
             influence by exactly {SINGLE_FLIP_SHIFT} and stayed equal, which is the check's \
             measured blind spot at parity and at the constant bit",
            arms.iter()
                .filter(|a| a.control.is_some() && !a.equal_within_class(&classes))
                .count(),
            arms.iter().filter(|a| a.control.is_some()).count(),
            arms.iter()
                .filter(|a| a.control.is_some() && a.equal_within_class(&classes))
                .count()
        );
    });
}
