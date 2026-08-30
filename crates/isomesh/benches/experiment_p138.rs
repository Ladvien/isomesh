//! **P-138 — what the case explosion would cost at tricubic, priced before anyone proposes it.**
//!
//! Ticket: R-138. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p138
//! ```
//!
//! Writes `docs/experiments/p-138.csv`.
//!
//! # What was missing
//!
//! **The 2026-08-23 memo killed Bernstein–Bézier form for the trilinear by four
//! separate arguments and then wrote down the one condition under which all four
//! die** (`docs/research/2026-08-23-unmined-mathematics-for-meshing.md:197-225`,
//! and the standing prohibition at `:265` — *"Bernstein–Bézier form — **Dead. Do
//! not attempt.**"*). The four, verbatim in substance:
//!
//! 1. *"For a degree-(1,1,1) tensor-product function the Bernstein basis is
//!    `B₀(t) = 1−t`, `B₁(t) = t`, so the **eight Bernstein coefficients are
//!    identically the eight corner values**"* (`:202-203`). The control net is
//!    the data; the change of basis is the identity map and buys nothing.
//! 2. *"De Casteljau subdivision of a multi-affine function produces sub-cell
//!    coefficients that are exactly `f` at the sub-cell corners — subdivision
//!    **is** resampling"* (`:203-205`), so *"the convex-hull enclosure … is
//!    **exact at every level** — Garloff's vertex condition holds unconditionally
//!    for multi-affine functions, so overestimation is zero"* (`:207-208`).
//! 3. *"The **sign pattern of the coefficients is the Marching Cubes case
//!    index**"* (`:211`) — the object already known to be insufficient to resolve
//!    the ambiguity, so variation-diminishing restates the ambiguity rather than
//!    resolving it.
//! 4. *"The sharpness theorems are vacuous here. The canonical pair is linear in
//!    1/degree under degree elevation (Rivlin 1970) and quadratic in box width
//!    under subdivision (Stahl 1995) — both describe overestimation that is
//!    already zero"* (`:214-216`).
//!
//! And then: ***"Revisit the entire hypothesis if the crate ever adopts a
//! tricubic reconstruction filter. At degree ≥ 2 per variable every part of it
//! becomes true, useful, and as far as I can tell unpublished for isosurface
//! topology"*** (`:223-225`). That is a prediction, in a document, with nothing
//! behind it. **`P-157` proposes exactly the adoption** — *"Raising the order,
//! which is the lever `✗42` did not pull"* (`FINDINGS.md:25370-25374`), a filter
//! of approximation order 4, and it records a `cases_invalidated` column without
//! anywhere saying how many cases there would be to invalidate. `✗42 / M-359`
//! (`FINDINGS.md:9708`) moved the *knot* at fixed order; nothing in this
//! repository has ever moved the *degree*, and nothing has priced it.
//!
//! This row is the price tag, and it is only the price tag. It proposes no
//! filter, changes no shipped path and makes no recommendation: it computes three
//! numbers per candidate degree — the BKK critical-point bound, the case-space
//! size, and the Bernstein verdict with its reason — so that `P-157` argues
//! against arithmetic instead of against a memory of a memo.
//!
//! The trilinear's own bound is `2`, and it is not taken on trust either:
//! `docs/research/2026-08-29-phase-27-axes-and-vocabulary-v2.md:195` records
//! *"`MV = 8 − 6 + 0 = 2`: a trilinear has at most **2** critical points in the
//! complex torus"*, `P-137` computes it by convex hull, and **this harness
//! recomputes it from the same inclusion–exclusion as the calibration of the
//! whole pipeline**. A tricubic number produced by a pipeline that cannot
//! reproduce `2` is not a measurement.
//!
//! # `TriPoly`, and why it is not a second copy of `common::poly::Poly`
//!
//! `common::poly::Poly` is a polynomial in the **eight corner values** `f0..f7`
//! — the object `P-127`'s hyperdeterminant identity lives in, where the *cell's
//! data* are the variables. Everything in this row is about the **three spatial
//! variables** `x, y, z`, where a tricubic reconstruction has `4³ = 64`
//! coefficients and the corner values are not variables at all. Those are
//! different objects, not two spellings of one, so this bench carries its own
//! `TriPoly`: a flat `Vec<i128>` on a `(d+1)³` exponent cube with `d ≤ 3`. It has
//! no `add`, no `mul` and no `sub`, because this row never needs them.
//!
//! What *is* consumed from `common::poly`, unchanged: `Rng` (the seeded
//! SplitMix64 every Phase 27 fixture draws from, so the instances are the same on
//! every host) and `octahedral_relabellings` + `relabel` (the generated
//! order-48 group, used for the trilinear case-class count that calibrates the
//! symmetry-quotient arithmetic).
//!
//! # Arms
//!
//! | arm | `reconstruction` | degree per variable | control values | is_control |
//! |---|---|---|---|---|
//! | 1 | `trilinear` | 1 | 8 | **yes** — the registered calibration |
//! | 2 | `triquadratic` | 2 | 27 | no — **extra arm**, not registered |
//! | 3 | `tricubic` | 3 | 64 | no — the row |
//!
//! **The triquadratic arm is an addition and is marked as one.** The registration
//! names two degrees; the memo's sentence is *"at degree ≥ 2 per variable"*, and
//! degree 2 is the cheapest place that sentence can be tested. It costs one more
//! row and it turns out to carry the only tractable case space of the three,
//! which is the sharpest thing this row has to give `P-157`. It is recorded with
//! `is_control=false` and is never used to decide `c1_holds` or `c2_holds`.
//!
//! # C1 — the mixed volume, by two exact routes that must agree
//!
//! For a tensor-product `f` of degree `d` in each variable, `∂f/∂x` has degree
//! `d−1` in `x` and `d` in `y` and `z`, so the three Newton polytopes of
//! `∇f = 0` are the axis-aligned boxes
//!
//! ```text
//! P₁ = [0,d−1] × [0,d]   × [0,d]
//! P₂ = [0,d]   × [0,d−1] × [0,d]
//! P₃ = [0,d]   × [0,d]   × [0,d−1]
//! ```
//!
//! **Route A, polarisation / inclusion–exclusion** — the same identity `P-137`
//! applies, and a Minkowski sum of boxes is a box whose sides add:
//!
//! ```text
//! MV(P₁,P₂,P₃) = Σ_{∅≠S⊆{1,2,3}} (−1)^(3−|S|) · Vol(Σ_{i∈S} Pᵢ)
//! ```
//!
//! At `d = 1` that is `0 − 6 + 8 = 2`, exactly the axes document's `8 − 6 + 0`
//! read in the other order. At `d = 3`: `54 − 450 + 512 = 116`.
//!
//! **Route B, the permanent.** The mixed volume of `n` boxes in `ℝⁿ` is the
//! permanent of their side-length matrix. For `d = 3` that matrix is
//! `[[2,3,3],[3,2,3],[3,3,2]]` and its permanent is `8 + 3·18 + 2·27 = 116`.
//! Route B shares no line of code with route A — one is a signed sum over `2³−1`
//! subsets of boxes, the other an unsigned sum over `3!` permutations of columns
//! — so `mv_routes_agree` is a real cross-check and not a restatement.
//!
//! By Bernstein's theorem the mixed volume **is** the bound on isolated
//! solutions in `(ℂ*)³`, so `critical_points_bound == mixed_volume` on every row
//! by definition. The two columns are both recorded anyway, because the
//! registration names both and because a future sharper *real* bound would move
//! one without moving the other.
//!
//! # C1's falsifier, and where the hopelessness actually lives
//!
//! The falsifier reads: *"C1 by a mixed volume that makes the tricubic case space
//! computationally hopeless"*. Those are two different quantities and the row
//! records both, with the bar fixed here, before the numbers:
//!
//! - **The bound.** `116` isolated critical points per cell against the
//!   trilinear's `2` — a `58×` rise in the root-isolation budget of a
//!   Plantinga–Vegter-style subdivision solver, which is expensive and is not
//!   hopeless. `mv_vs_trilinear` is that ratio.
//! - **The case space.** The sign pattern is over the `(d+1)³` control values, so
//!   the case count is `2^((d+1)³)`: `2^8 = 256` for the trilinear, `2^27` for
//!   the triquadratic, `2^64` for the tricubic. **`2^64` is recorded as the token
//!   `2^64` and never as a float** — the decimal is in `case_count_decimal` as an
//!   exact `u128`. The bar for "tabulable" is this machine's RAM, `32 GiB = 2^35`
//!   bytes at one byte per case: `2^8` is 256 B, `2^27` is 128 MiB, `2^64` is
//!   18.4 EB. `case_space_tractable` is `net_sign_bits <= ram_bytes_log2`, and
//!   both sides are columns so the verdict is auditable from the file alone.
//!
//! So the sharper statement, which the registration did not anticipate, is that
//! **the hopelessness is in the table and not in the bound** — and that the
//! middle degree the memo's *"≥ 2"* also covers has a case space of `2^27`, or
//! about 1.4 M classes after the order-96 quotient, which is tabulable on a
//! laptop. `c1_holds` therefore carries the clause the registration actually
//! wrote — *"is computed and reported"* — gated on the calibration and the
//! two-route agreement, and the judgement half is in the file as numbers rather
//! than as an adjective.
//!
//! `case_classes_lower_bound` is `⌈2^((d+1)³) / 96⌉`, the orbit-counting floor
//! under the order-48 octahedral group on the control net times complementation.
//! It is a floor and is labelled one; `trilinear_case_classes` is the **exact**
//! enumerated orbit count of the 256 corner patterns under that same group,
//! recorded identically on every row as the calibration that the group and the
//! quotient are the ones being claimed.
//!
//! # C2 — the four arguments, re-checked rather than re-read
//!
//! Each argument is decided by a measured integer over `INSTANCES` seeded
//! instances per arm, not by an appeal to the memo. An instance is a **Bernstein
//! control net** with non-zero integer entries in `[−9, 9]` — the natural object
//! for a reconstruction filter, and the one direction of the basis change that
//! needs no division, so every number below is exact `i128` arithmetic with no
//! rational type anywhere.
//!
//! | argument | measured as | trilinear | tricubic |
//! |---|---|---|---|
//! | 1 — the net *is* the samples | `control_points_matching_samples` of `control_points_total`: `b_ijk` against `f(i/d, j/d, k/d)`, exactly | all | corners only |
//! | 2 — subdivision *is* resampling | `subdivision_sites_matching` of total: the 8 de Casteljau half-nets at `t = ½` against `f` at the sub-cell abscissae, exactly | all | corners only |
//! | 2′ — the enclosure is exact at every level | `enclosure_narrowed_instances`: does the union of the 8 sub-enclosures sit **strictly inside** the parent's `[min b, max b]` | none | all |
//! | 3 — the signs *are* the case index | `net_sign_bits` and `hidden_edge_sign_changes` | 8, 0 | 64, >0 |
//! | 4 — the sharpness theorems are vacuous | `sharpness_theorems_vacuous`, i.e. no overestimation for Rivlin's and Stahl's rates to describe | true | false |
//!
//! Argument 2′ is what makes "strictly" a measurement rather than a hope.
//! Bernstein enclosures are nested under subdivision, so the union of the
//! children can only be narrower or equal; **narrower is a certificate that the
//! parent enclosure was a strict over-estimate**, with no need to bound the true
//! range of `f` from above. For a multi-affine `f` the range over the cell is
//! exactly `[min corner, max corner]` and every sub-net coefficient is a value of
//! `f` inside it, so the union equals the parent and `narrowed` is `false` — the
//! trilinear arm is the proof that the instrument can report `false`. Argument 4
//! is derived from 2′ and not measured separately: both named theorems bound
//! *overestimation*, so a measured overestimation of zero makes both vacuous and
//! a measured non-zero one makes both bite.
//!
//! `hidden_edge_sign_changes` is argument 3's substance. It counts, over the 12
//! cell edges, lines of the control net whose **two endpoints agree in sign** but
//! which carry a sign change between adjacent control values somewhere along
//! their interior — information the 8-corner case index structurally cannot hold.
//! At `d = 1` a cell-edge line has two entries and both are endpoints, so the
//! count is `0` by construction and the case index *is* the sign pattern. Zeros
//! never enter the comparison: instance coefficients are drawn non-zero, and
//! "opposite signs" is `a.signum() * b.signum() < 0`.
//!
//! `bernstein_route_alive` is the token `alive` or `dead`, and it is **derived
//! from the four measured booleans** so it cannot disagree with the columns
//! beside it: `alive` iff the net is not the samples **and** subdivision is not
//! resampling **and** the signs are not the case index **and** the enclosure
//! narrows. `bernstein_reason` spells the four out as a `|`-joined token built
//! from the same booleans. `subdivision_exactness` is
//! `exact_resampling` / `inexact_strict_overestimate` / `inexact_but_not_narrowed`
//! — the third is a registered outcome that this harness can report and is not
//! expected to.
//!
//! `c2_holds` is the tricubic verdict being `alive`, which is the memo's own
//! prediction; a `dead` there falsifies C2 and contradicts the memo, which is
//! precisely why the prediction was worth registering.
//!
//! # Both verdict columns are global, and carry the same value on every row
//!
//! C1 is *"the mixed volume … is computed and reported"* and C2 is *"the
//! Bernstein-form arguments … are re-checked at degree 3"*. Neither is a
//! per-arm property: C1 needs the trilinear arm to calibrate the tricubic one,
//! and C2 is a statement about degree 3 specifically. So `c1_holds` and
//! `c2_holds` are computed once across all arms and written identically to every
//! row, and the per-arm content lives in the per-arm columns. There is no
//! `c3_holds`; the registration names two clauses.
//!
//! # SHARE, recomputed before the numbers
//!
//! **`SHARE: none — this prices a decision, it does not make one.`** Discharged
//! as written, and recomputed rather than copied: this harness calls no
//! extractor, samples no field, touches no `crates/isomesh/src/**` path and
//! proposes no landing, so there is no stage whose share of an extraction could
//! be taken. What stands in a share's place is one exact integer per clause over
//! a denominator that is exact by construction:
//!
//! | clause | quantity | denominator | exact because |
//! |---|---|---|---|
//! | C1 | `mixed_volume` | — | a permanent of a `3×3` integer matrix |
//! | C1 | `case_count_decimal` | `2^((d+1)³)` | the control count |
//! | C2 | `control_points_matching_samples` | `control_points_total` | `INSTANCES · (d+1)³` |
//! | C2 | `enclosure_narrowed_instances` | `INSTANCES` | the fixture |
//!
//! No wall clock is read anywhere in this file, so `M-280`'s 1.45× governor
//! swing cannot reach any column. Every number is an integer, a ratio of
//! integers, or a token.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record` and every panic starts `VOID: `.
//!
//! - **`mixed_volume` on the trilinear arm must be exactly `2`** — the
//!   registration's own control, verbatim: *"the trilinear must be run through
//!   the same pipeline and reproduce `MV = 2`, or the tricubic number has no
//!   calibration"*. Column: `mixed_volume`.
//! - **`estimated_case_count` on the trilinear arm must be `2^8`, and
//!   `case_count_decimal` exactly `256`** — the second half of the same
//!   calibration, and the reason the tricubic count is a token rather than a
//!   float. Columns: `estimated_case_count`, `case_count_decimal`.
//! - **The two mixed-volume routes must agree on every arm.** A tricubic `116`
//!   from one route is a number; from two disjoint routes it is a measurement.
//!   Column: `mv_routes_agree`.
//! - **`bernstein_route_alive` on the trilinear arm must be `dead`.** This is the
//!   control that licenses the tricubic `alive`: the same four checks, the same
//!   code, the same instances generator, reaching the *opposite* verdict one
//!   degree down. Without it an `alive` at degree 3 is a property of the pipeline
//!   rather than of the degree. Column: `bernstein_route_alive`.
//! - **On the trilinear arm, all four arguments must reproduce the memo:** every
//!   control point equal to its sample, every sub-net coefficient equal to its
//!   sub-cell sample, `hidden_edge_sign_changes == 0`, and
//!   `enclosure_narrowed_instances == 0`. Four columns, and a pipeline that fails
//!   any of them is not measuring the memo's arguments.
//! - **The corner sites must match on *every* arm.** Bernstein endpoint
//!   interpolation makes `b` at an index vertex exactly `f` at the cell corner
//!   for any degree, so `corner_sites_matching == 8 · INSTANCES` is an
//!   independent exactness check on the Bernstein→monomial conversion. It is what
//!   licenses reading a *mismatch* at the interior sites as a fact about the
//!   degree rather than a bug in the change of basis. Column:
//!   `corner_sites_matching`.
//! - **No arm may be measured on a degenerate fixture.** `degenerate_nets` counts
//!   instances whose parent enclosure has zero width, for which the narrowing
//!   comparison says nothing; it is asserted strictly below `INSTANCES`. Column:
//!   `degenerate_nets`.
//! - **`case_classes_lower_bound` must not exceed the enumerated trilinear
//!   count** on the trilinear arm — `⌈256/96⌉ = 3 ≤ 14`. A floor above the exact
//!   value would mean the group order used for the quotient is wrong, which would
//!   silently corrupt the tricubic floor where no exact value can be enumerated.
//!
//! # Determinism
//!
//! One thread, no clock, no float in any gated quantity, no map iteration. The
//! instance stream is `common::poly::Rng` (SplitMix64) seeded `SEED ^ degree`,
//! stated in the `seed` column. Every arithmetic path is `i128`; the largest
//! absolute intermediate is tracked per arm and asserted below `2^100`, so the
//! word "exactly" above is checked rather than asserted by prose.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use std::collections::BTreeSet;

use common::poly::{Rng, octahedral_relabellings, relabel};

/// Instances per arm. The fixture is exhaustive in structure and random only in
/// coefficients, so this is a witness count rather than a sample size: every
/// argument below is decided by "does any instance narrow" or "do all sites
/// match", and 2048 non-zero nets is far past the point where either flips.
const INSTANCES: usize = 2048;

/// The instance stream's base seed. `"P138"` twice, so the number in the `seed`
/// column is recognisable in the file.
const SEED: u64 = 0x5031_3338_5031_3338;

/// The three degrees, ascending, with the calibration first.
const DEGREES: [usize; 3] = [1, 2, 3];

/// Arm names, parallel to [`DEGREES`].
const NAMES: [&str; 3] = ["trilinear", "triquadratic", "tricubic"];

/// The order-48 octahedral group on the control net, times complementation.
///
/// The quotient the case-class floor divides by. `48` is
/// [`octahedral_relabellings`]' own count and `2` is `f → −f`.
const GROUP_ORDER: u128 = 96;

/// `log2` of this machine's RAM in bytes: 32 GiB.
///
/// The bar `case_space_tractable` is measured against, at one byte per case.
/// Fixed here, before any number, so the word "hopeless" in the registration's
/// falsifier is decided by arithmetic against a stated quantity.
const RAM_BYTES_LOG2: usize = 35;

/// Every exact intermediate must stay this far inside `i128`.
///
/// Asserted rather than argued: the header claims exactness throughout, and a
/// silent release-profile wrap would make that claim false while every column
/// still looked plausible.
const EXACTNESS_CEILING: i128 = 1 << 100;

/// `n choose k`, exact.
///
/// The incremental form `acc * (n − i) / (i + 1)` is exact at every step because
/// after `i + 1` factors the accumulator is `C(n, i+1)`, an integer.
fn binom(n: usize, k: usize) -> i128 {
    if k > n {
        return 0;
    }
    let mut acc: i128 = 1;
    for i in 0..k {
        acc = acc * (n - i) as i128 / (i + 1) as i128;
    }
    acc
}

/// The Bernstein-to-monomial matrix for degree `d`, row-major `[out·side + inp]`.
///
/// `Σ_j b_j C(d,j) t^j (1−t)^(d−j)`, so the coefficient of `t^m` is
/// `Σ_{j≤m} b_j C(d,j) C(d−j, m−j) (−1)^(m−j)`. **Pure integers, no division** —
/// which is why instances are generated as control nets and converted this way
/// rather than the other, and why nothing in this bench needs a rational type.
fn bernstein_to_monomial(d: usize) -> Vec<i128> {
    let side = d + 1;
    let mut m = vec![0i128; side * side];
    for out in 0..side {
        for inp in 0..=out {
            let sign = if (out - inp).is_multiple_of(2) { 1 } else { -1 };
            m[out * side + inp] = sign * binom(d, inp) * binom(d - inp, out - inp);
        }
    }
    m
}

/// The two de Casteljau half-nets at `t = ½`, each over the common denominator
/// `2^d`, row-major `[out·side + inp]`.
///
/// `b^(r)_i = 2^(−r) Σ_m C(r,m) b_(i+m)`, so the left net's entry `r` is
/// `2^(d−r) Σ_m C(r,m) b_m` and the right net's entry `r` is
/// `2^r Σ_m C(d−r,m) b_(r+m)`. Both are integral because `r ≤ d`, so subdivision
/// stays exact with a single shared power-of-two denominator instead of a
/// per-entry one.
fn de_casteljau_halves(d: usize) -> (Vec<i128>, Vec<i128>) {
    let side = d + 1;
    let mut left = vec![0i128; side * side];
    let mut right = vec![0i128; side * side];
    for r in 0..side {
        let lscale = 1i128 << (d - r);
        for m in 0..=r {
            left[r * side + m] = lscale * binom(r, m);
        }
        let rscale = 1i128 << r;
        for m in 0..=(d - r) {
            right[r * side + (r + m)] = rscale * binom(d - r, m);
        }
    }
    (left, right)
}

/// Flat index into a `side³` cube with `t` along `axis` and `(a, b)` the other
/// two coordinates in increasing axis order.
///
/// The cube is `(i·side + j)·side + k` throughout this file, `i` along `x`.
fn flat(side: usize, axis: usize, a: usize, b: usize, t: usize) -> usize {
    let (i, j, k) = match axis {
        0 => (t, a, b),
        1 => (a, t, b),
        _ => (a, b, t),
    };
    (i * side + j) * side + k
}

/// Apply one `side × side` integer matrix along each of the three axes.
///
/// A tensor-product basis change is a per-axis one, so this single function does
/// Bernstein→monomial (the same matrix three times) and de Casteljau subdivision
/// into one of the eight sub-boxes (a per-axis choice of left or right half).
fn apply_per_axis(side: usize, cube: &[i128], mats: [&[i128]; 3]) -> Vec<i128> {
    let mut cur = cube.to_vec();
    for (axis, m) in mats.into_iter().enumerate() {
        let mut next = vec![0i128; cur.len()];
        for a in 0..side {
            for b in 0..side {
                for out in 0..side {
                    let mut acc = 0i128;
                    for inp in 0..side {
                        acc += m[out * side + inp] * cur[flat(side, axis, a, b, inp)];
                    }
                    next[flat(side, axis, a, b, out)] = acc;
                }
            }
        }
        cur = next;
    }
    cur
}

/// `den^(3d) · f(num₀/den, num₁/den, num₂/den)` for a monomial cube, exactly.
///
/// Every term carries `den^(3d − p − q − r)` with a non-negative exponent
/// because `p, q, r ≤ d`, so the whole evaluation is one integer sum and no
/// rational arithmetic appears anywhere. `den_powers[e]` must be `den^e` for
/// `e` in `0..=3d`.
fn eval_scaled(side: usize, mono: &[i128], num: [i128; 3], den_powers: &[i128]) -> i128 {
    let d = side - 1;
    let mut np = [[1i128; 4]; 3];
    for (axis, row) in np.iter_mut().enumerate() {
        for i in 1..side {
            row[i] = row[i - 1] * num[axis];
        }
    }
    let mut total = 0i128;
    for p in 0..side {
        for q in 0..side {
            for r in 0..side {
                let c = mono[(p * side + q) * side + r];
                if c != 0 {
                    total += c * np[0][p] * np[1][q] * np[2][r] * den_powers[3 * d - p - q - r];
                }
            }
        }
    }
    total
}

/// `den^e` for `e` in `0..=3d`.
fn power_table(den: i128, d: usize) -> Vec<i128> {
    let mut out = vec![1i128; 3 * d + 1];
    for e in 1..out.len() {
        out[e] = out[e - 1] * den;
    }
    out
}

/// The three Newton polytopes of `∇f = 0` for degree `d` per variable, as boxes
/// given by their side lengths.
///
/// `∂f/∂x_i` drops one degree in `x_i` and keeps `d` in the other two, so box
/// `i` is `d` everywhere except `d − 1` on the diagonal.
fn gradient_newton_boxes(d: usize) -> [[u64; 3]; 3] {
    let mut out = [[d as u64; 3]; 3];
    for (i, row) in out.iter_mut().enumerate() {
        row[i] = d as u64 - 1;
    }
    out
}

/// The three Newton polytopes as one CSV-safe token: side lengths joined by
/// `x` within a box and by `|` between boxes.
fn newton_boxes_token(d: usize) -> String {
    let boxes = gradient_newton_boxes(d);
    let one = |r: [u64; 3]| format!("{}x{}x{}", r[0], r[1], r[2]);
    format!("{}|{}|{}", one(boxes[0]), one(boxes[1]), one(boxes[2]))
}

/// Mixed volume by polarisation: the signed inclusion–exclusion over subsets.
///
/// `MV = Σ_{∅≠S} (−1)^(3−|S|) Vol(Σ_{i∈S} Pᵢ)`, and a Minkowski sum of boxes is
/// the box whose sides add — so `Vol` is one product of three integers and the
/// whole computation is exact. Route A of two.
fn mixed_volume_by_polarisation(boxes: &[[u64; 3]; 3]) -> i128 {
    let mut total: i128 = 0;
    for mask in 1u32..8 {
        let mut sides = [0u64; 3];
        let mut chosen = 0u32;
        for (i, row) in boxes.iter().enumerate() {
            if mask & (1 << i) != 0 {
                chosen += 1;
                for (a, s) in sides.iter_mut().enumerate() {
                    *s += row[a];
                }
            }
        }
        let vol = i128::from(sides[0]) * i128::from(sides[1]) * i128::from(sides[2]);
        if (3 - chosen).is_multiple_of(2) {
            total += vol;
        } else {
            total -= vol;
        }
    }
    total
}

/// Mixed volume as the permanent of the side-length matrix.
///
/// A standard identity for boxes, and route B of two: an unsigned sum over the
/// `3!` column permutations, sharing no line with route A's signed sum over the
/// `2³ − 1` subsets of boxes. Disagreement between the two is a `VOID`.
fn mixed_volume_by_permanent(boxes: &[[u64; 3]; 3]) -> i128 {
    const PERMS: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    PERMS
        .iter()
        .map(|p| (0..3).map(|i| i128::from(boxes[i][p[i]])).product::<i128>())
        .sum()
}

/// The exact orbit count of the 256 corner sign patterns under the order-48
/// octahedral group times complementation.
///
/// Enumerated, not transcribed: the folklore numbers for the Marching Cubes case
/// classes (14, 15, 22, 23) each name a *different* group, so the count is
/// computed here and the group it belongs to is [`GROUP_ORDER`], recorded beside
/// it. Its only job is to calibrate the symmetry quotient that produces
/// `case_classes_lower_bound` at degrees where no orbit count can be enumerated.
fn corner_sign_pattern_classes() -> usize {
    let perms = octahedral_relabellings();
    let mut canonical = BTreeSet::new();
    for mask in 0u16..256 {
        let bits: [u8; 8] = std::array::from_fn(|i| ((mask >> i) & 1) as u8);
        let mut best = u16::MAX;
        for perm in &perms {
            let moved = relabel(perm, &bits);
            let mut m: u16 = 0;
            for (i, b) in moved.iter().enumerate() {
                m |= u16::from(*b) << i;
            }
            best = best.min(m).min(m ^ 0xFF);
        }
        canonical.insert(best);
    }
    canonical.len()
}

/// Everything one arm measures, before any verdict exists.
struct Arm {
    /// `trilinear` / `triquadratic` / `tricubic`.
    name: &'static str,
    /// Degree per variable: 1, 2, 3.
    degree: usize,
    /// Route A.
    mv_polarisation: i128,
    /// Route B.
    mv_permanent: i128,
    /// `(d+1)³` — the control count, and the width of the sign pattern.
    net_sign_bits: usize,
    /// `2^net_sign_bits`, exact.
    case_count_decimal: u128,
    /// `⌈case_count / GROUP_ORDER⌉`, an orbit-counting floor.
    case_classes_lower_bound: u128,
    /// Argument 1: `b_ijk` against `f(i/d, j/d, k/d)`.
    sites_matching: u64,
    /// `INSTANCES · (d+1)³`.
    sites_total: u64,
    /// The subset of [`Arm::sites_matching`] at the 8 index vertices.
    corner_sites_matching: u64,
    /// Argument 2: sub-net coefficients against `f` at the sub-cell abscissae.
    sub_sites_matching: u64,
    /// `INSTANCES · 8 · (d+1)³`.
    sub_sites_total: u64,
    /// Argument 2′: instances whose child enclosures sit strictly inside the
    /// parent's.
    narrowed: usize,
    /// Instances whose parent enclosure has zero width, for which 2′ says
    /// nothing.
    degenerate: usize,
    /// Narrowing as parts per million of the parent width, over the
    /// non-degenerate instances.
    shrink_ppm_min: i128,
    /// The other end of the same range.
    shrink_ppm_max: i128,
    /// Argument 3: cell-edge control lines whose endpoints agree in sign but
    /// whose interior carries a sign change.
    hidden_edge_sign_changes: u64,
    /// The largest absolute exact intermediate this arm produced.
    max_magnitude: i128,
    /// This arm's instance-stream seed.
    seed: u64,
}

impl Arm {
    /// Argument 1's verdict: is the control net literally the samples?
    fn net_is_the_samples(&self) -> bool {
        self.sites_matching == self.sites_total
    }

    /// Argument 2's verdict: is de Casteljau subdivision resampling?
    fn subdivision_is_resampling(&self) -> bool {
        self.sub_sites_matching == self.sub_sites_total
    }

    /// Argument 2′'s verdict: does the enclosure strictly over-estimate?
    fn enclosure_overestimates(&self) -> bool {
        self.narrowed > 0
    }

    /// Argument 4's verdict, derived from 2′: Rivlin's and Stahl's rates both
    /// bound overestimation, so a measured zero makes both vacuous.
    fn sharpness_vacuous(&self) -> bool {
        !self.enclosure_overestimates()
    }

    /// Argument 3's verdict: is the sign pattern of the coefficients the
    /// Marching Cubes case index? Eight bits **and** nowhere for extra sign
    /// information to hide.
    fn signs_are_the_case_index(&self) -> bool {
        self.net_sign_bits == 8 && self.hidden_edge_sign_changes == 0
    }

    /// All four arguments dead, which is the route being alive.
    fn alive(&self) -> bool {
        !self.net_is_the_samples()
            && !self.subdivision_is_resampling()
            && !self.signs_are_the_case_index()
            && self.enclosure_overestimates()
    }

    /// The four verdicts as a `|`-joined CSV-safe token, built from the same
    /// booleans as [`Arm::alive`] so the prose cannot drift from the columns.
    fn reason(&self) -> String {
        let parts = [
            if self.net_is_the_samples() {
                "coefficients_are_the_samples"
            } else {
                "control_net_is_not_the_samples"
            },
            if self.subdivision_is_resampling() {
                "subdivision_is_resampling"
            } else {
                "subdivision_strictly_refines"
            },
            if self.signs_are_the_case_index() {
                "signs_are_the_case_index"
            } else {
                "signs_finer_than_the_case_index"
            },
            if self.sharpness_vacuous() {
                "sharpness_theorems_vacuous"
            } else {
                "sharpness_theorems_bite"
            },
        ];
        format!(
            "{}|{}",
            if self.alive() { "alive" } else { "dead" },
            parts.join("|")
        )
    }

    /// `exact_resampling` when subdivision reproduces `f`; otherwise whether the
    /// enclosure was shown to be a strict over-estimate. The third token is a
    /// registered outcome this harness can report and does not expect.
    fn subdivision_exactness(&self) -> &'static str {
        if self.subdivision_is_resampling() {
            "exact_resampling"
        } else if self.enclosure_overestimates() {
            "inexact_strict_overestimate"
        } else {
            "inexact_but_not_narrowed"
        }
    }

    /// `2^n` as an exponent-bearing token. A case count of `2^64` cannot be a
    /// float and must not be one.
    fn case_count_token(&self) -> String {
        format!("2^{}", self.net_sign_bits)
    }

    /// One byte per case against [`RAM_BYTES_LOG2`].
    fn case_space_tractable(&self) -> bool {
        self.net_sign_bits <= RAM_BYTES_LOG2
    }
}

/// Measure one arm: the mixed volume, the case space, and the four Bernstein
/// arguments over [`INSTANCES`] seeded control nets.
fn measure(degree: usize, name: &'static str) -> Arm {
    let side = degree + 1;
    assert!(
        (1..=3).contains(&degree),
        "the fixed-size exponent arrays in `eval_scaled` cover degree 1..=3, not {degree}"
    );
    let cells = side * side * side;
    assert!(
        cells < 128,
        "`case_count_decimal` is a `u128`, so the control count must stay below 128"
    );

    let boxes = gradient_newton_boxes(degree);
    let b2m = bernstein_to_monomial(degree);
    let (left, right) = de_casteljau_halves(degree);

    // Abscissae for argument 1 are `i/d`; for argument 2 they are `i/(2d)` in the
    // parent's coordinates, and the sub-nets share the denominator `2^(3d)`.
    let abscissa_powers = power_table(degree as i128, degree);
    let sub_powers = power_table(2 * degree as i128, degree);
    let sub_denominator = 1i128 << (3 * degree);

    let seed = SEED ^ degree as u64;
    let mut rng = Rng::new(seed);

    let mut arm = Arm {
        name,
        degree,
        mv_polarisation: mixed_volume_by_polarisation(&boxes),
        mv_permanent: mixed_volume_by_permanent(&boxes),
        net_sign_bits: cells,
        case_count_decimal: 1u128 << cells,
        case_classes_lower_bound: (1u128 << cells).div_ceil(GROUP_ORDER),
        sites_matching: 0,
        sites_total: 0,
        corner_sites_matching: 0,
        sub_sites_matching: 0,
        sub_sites_total: 0,
        narrowed: 0,
        degenerate: 0,
        shrink_ppm_min: i128::MAX,
        shrink_ppm_max: i128::MIN,
        hidden_edge_sign_changes: 0,
        max_magnitude: 0,
        seed,
    };

    let mut net = vec![0i128; cells];
    for _ in 0..INSTANCES {
        // Non-zero coefficients in [−9, 9], symmetric about zero: a zero control
        // value would make "opposite signs" undefined for argument 3 and would
        // put the fixture on the degenerate stratum M-48 already owns.
        for c in &mut net {
            let v = rng.next_i64_in(0, 18);
            *c = i128::from(if v < 9 { v - 9 } else { v - 8 });
        }

        let mono = apply_per_axis(side, &net, [&b2m, &b2m, &b2m]);
        for c in &mono {
            arm.max_magnitude = arm.max_magnitude.max(c.abs());
        }

        // ── argument 1: is the net the samples? ──
        for i in 0..side {
            for j in 0..side {
                for k in 0..side {
                    let num = [i as i128, j as i128, k as i128];
                    let exact = eval_scaled(side, &mono, num, &abscissa_powers);
                    let claimed = net[(i * side + j) * side + k] * abscissa_powers[3 * degree];
                    arm.max_magnitude = arm.max_magnitude.max(exact.abs()).max(claimed.abs());
                    arm.sites_total += 1;
                    let is_vertex = [i, j, k].iter().all(|t| *t == 0 || *t == degree);
                    if exact == claimed {
                        arm.sites_matching += 1;
                        if is_vertex {
                            arm.corner_sites_matching += 1;
                        }
                    }
                }
            }
        }

        // ── argument 3: sign information the 8-corner case index cannot hold ──
        for axis in 0..3 {
            let (oa, ob) = match axis {
                0 => (1, 2),
                1 => (0, 2),
                _ => (0, 1),
            };
            for ea in [0usize, degree] {
                for eb in [0usize, degree] {
                    let mut coord = [0usize; 3];
                    coord[oa] = ea;
                    coord[ob] = eb;
                    let mut line = [0i128; 4];
                    for (t, cell) in line.iter_mut().enumerate().take(side) {
                        coord[axis] = t;
                        *cell = net[(coord[0] * side + coord[1]) * side + coord[2]];
                    }
                    let ends_agree = line[0].signum() * line[degree].signum() > 0;
                    let interior_flips =
                        (0..degree).any(|t| line[t].signum() * line[t + 1].signum() < 0);
                    if ends_agree && interior_flips {
                        arm.hidden_edge_sign_changes += 1;
                    }
                }
            }
        }

        // ── arguments 2 and 2′: subdivision, resampling and the enclosure ──
        let mut parent_lo = i128::MAX;
        let mut parent_hi = i128::MIN;
        for c in &net {
            let scaled = c * sub_denominator;
            parent_lo = parent_lo.min(scaled);
            parent_hi = parent_hi.max(scaled);
        }
        let mut union_lo = i128::MAX;
        let mut union_hi = i128::MIN;
        for mask in 0usize..8 {
            let mats: [&[i128]; 3] = std::array::from_fn(|a| {
                if ((mask >> a) & 1) == 1 {
                    right.as_slice()
                } else {
                    left.as_slice()
                }
            });
            let sub = apply_per_axis(side, &net, mats);
            for c in &sub {
                union_lo = union_lo.min(*c);
                union_hi = union_hi.max(*c);
                arm.max_magnitude = arm.max_magnitude.max(c.abs());
            }
            for i in 0..side {
                for j in 0..side {
                    for k in 0..side {
                        let num = [
                            (i + degree * (mask & 1)) as i128,
                            (j + degree * ((mask >> 1) & 1)) as i128,
                            (k + degree * ((mask >> 2) & 1)) as i128,
                        ];
                        let exact = eval_scaled(side, &mono, num, &sub_powers);
                        let idx = (i * side + j) * side + k;
                        // `sub[idx] / 2^(3d)` against `exact / (2d)^(3d)`, cross
                        // multiplied so the comparison stays in integers.
                        let lhs = sub[idx] * sub_powers[3 * degree];
                        let rhs = exact * sub_denominator;
                        arm.max_magnitude = arm.max_magnitude.max(lhs.abs()).max(rhs.abs());
                        arm.sub_sites_total += 1;
                        if lhs == rhs {
                            arm.sub_sites_matching += 1;
                        }
                    }
                }
            }
        }

        // Bernstein enclosures are nested under subdivision, so the children can
        // only be narrower or equal. If that fails, the subdivision matrices are
        // wrong and every argument-2 number is meaningless.
        assert!(
            union_lo >= parent_lo && union_hi <= parent_hi,
            "VOID: the {name} child enclosure [{union_lo}, {union_hi}] escapes its parent \
             [{parent_lo}, {parent_hi}], so the de Casteljau matrices are not a subdivision \
             and no argument-2 column means anything"
        );

        let parent_width = parent_hi - parent_lo;
        if parent_width == 0 {
            arm.degenerate += 1;
        } else {
            let union_width = union_hi - union_lo;
            if union_width < parent_width {
                arm.narrowed += 1;
            }
            let ppm = (parent_width - union_width) * 1_000_000 / parent_width;
            arm.shrink_ppm_min = arm.shrink_ppm_min.min(ppm);
            arm.shrink_ppm_max = arm.shrink_ppm_max.max(ppm);
        }
    }

    assert!(
        arm.max_magnitude < EXACTNESS_CEILING,
        "VOID: the {name} arm reached |{}|, which is within a factor of 2^28 of `i128`'s range, \
         so the exact arithmetic this row's every column rests on may have wrapped",
        arm.max_magnitude
    );
    arm
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-138");

    common::experiment::run(prereg, |run| {
        let arms: Vec<Arm> = DEGREES
            .iter()
            .zip(NAMES)
            .map(|(d, name)| measure(*d, name))
            .collect();
        let trilinear_classes = corner_sign_pattern_classes();

        let calibration = arms
            .iter()
            .find(|a| a.degree == 1)
            .expect("the trilinear arm is the first of DEGREES");
        let subject = arms
            .iter()
            .find(|a| a.degree == 3)
            .expect("the tricubic arm is the last of DEGREES");

        // ── vacuity controls ──

        assert!(
            calibration.mv_polarisation == 2,
            "VOID: the trilinear arm gives MV = {}, not 2, so the tricubic mixed volume has no \
             calibration and the registration's own control has failed \
             (docs/research/2026-08-29-phase-27-axes-and-vocabulary-v2.md:195)",
            calibration.mv_polarisation
        );
        assert!(
            calibration.case_count_token() == "2^8" && calibration.case_count_decimal == 256,
            "VOID: the trilinear case space came out {} = {}, not 2^8 = 256, so the pipeline that \
             produced the tricubic 2^64 cannot reproduce the count everyone already knows",
            calibration.case_count_token(),
            calibration.case_count_decimal
        );
        for a in &arms {
            assert!(
                a.mv_polarisation == a.mv_permanent,
                "VOID: the {} arm's two mixed-volume routes disagree, {} by polarisation against \
                 {} by permanent, so neither number is a measurement",
                a.name,
                a.mv_polarisation,
                a.mv_permanent
            );
            assert!(
                a.corner_sites_matching == 8 * INSTANCES as u64,
                "VOID: the {} arm matched {} of {} index-vertex sites, but Bernstein endpoint \
                 interpolation makes the control net equal f at a cell corner at every degree - \
                 so the Bernstein-to-monomial conversion is wrong and an interior mismatch says \
                 nothing about the degree",
                a.name,
                a.corner_sites_matching,
                8 * INSTANCES as u64
            );
            assert!(
                a.degenerate < INSTANCES,
                "VOID: all {INSTANCES} of the {} arm's nets have a zero-width enclosure, so the \
                 narrowing comparison that decides arguments 2' and 4 is being taken over a \
                 constant fixture",
                a.name
            );
            assert!(
                u128::from(a.net_sign_bits as u64) <= 128,
                "VOID: the {} arm's control count {} overflows the exact case-count arithmetic",
                a.name,
                a.net_sign_bits
            );
        }
        assert!(
            !calibration.alive(),
            "VOID: the trilinear arm came out `{}`, but the 2026-08-23 memo killed all four \
             Bernstein arguments at degree 1 \
             (docs/research/2026-08-23-unmined-mathematics-for-meshing.md:197-225). An `alive` \
             here would mean the tricubic verdict is a property of this pipeline rather than of \
             the degree",
            calibration.reason()
        );
        assert!(
            calibration.net_is_the_samples(),
            "VOID: the trilinear arm matched {} of {} sites, but the memo's argument 1 is that the \
             eight Bernstein coefficients ARE the eight corner values - a pipeline that cannot \
             reproduce that is not measuring the memo's arguments",
            calibration.sites_matching,
            calibration.sites_total
        );
        assert!(
            calibration.subdivision_is_resampling(),
            "VOID: the trilinear arm matched {} of {} subdivision sites, but the memo's argument 2 \
             is that de Casteljau subdivision of a multi-affine function IS resampling",
            calibration.sub_sites_matching,
            calibration.sub_sites_total
        );
        assert!(
            calibration.hidden_edge_sign_changes == 0,
            "VOID: the trilinear arm found {} hidden cell-edge sign changes, but a degree-1 \
             control line has two entries and both are endpoints, so the count is 0 by \
             construction - a non-zero means the line enumeration is not walking cell edges",
            calibration.hidden_edge_sign_changes
        );
        assert!(
            !calibration.enclosure_overestimates(),
            "VOID: the trilinear arm narrowed on {} of {INSTANCES} instances, but Garloff's vertex \
             condition holds unconditionally for a multi-affine function so the enclosure is exact \
             at every level and cannot narrow \
             (docs/research/2026-08-23-unmined-mathematics-for-meshing.md:207-208)",
            calibration.narrowed
        );
        assert!(
            calibration.case_classes_lower_bound <= trilinear_classes as u128,
            "VOID: the orbit-counting floor {} exceeds the enumerated trilinear class count {}, so \
             GROUP_ORDER = {GROUP_ORDER} is not the group being quotiented by and the tricubic \
             floor - where no orbit count can be enumerated - is corrupt",
            calibration.case_classes_lower_bound,
            trilinear_classes
        );
        assert!(
            subject.sites_total > 0 && subject.sub_sites_total > 0,
            "VOID: the tricubic arm measured no site at all, so C2 was never checked at degree 3"
        );

        // ── the two global verdicts ──

        // C1 is "is computed and reported": the calibration reproduces 2, both
        // routes agree on every arm, and the tricubic number is exact. What that
        // number *means* is reported as `mv_vs_trilinear`, `case_space_tractable`
        // and `case_count_decimal` rather than folded into a boolean.
        let c1 = calibration.mv_polarisation == 2
            && arms.iter().all(|a| a.mv_polarisation == a.mv_permanent)
            && subject.mv_polarisation > 0;
        // C2 is the memo's own prediction: alive at degree 3.
        let c2 = subject.alive();

        let trilinear_mv = calibration.mv_polarisation as f64;
        for a in &arms {
            // `degenerate < INSTANCES` is asserted above, so at least one
            // instance set the range and the sentinels are always overwritten.
            run.record(&[
                ("reconstruction", a.name.to_string()),
                ("degree_per_variable", a.degree.to_string()),
                ("mixed_volume", a.mv_polarisation.to_string()),
                ("critical_points_bound", a.mv_polarisation.to_string()),
                (
                    "bernstein_route_alive",
                    if a.alive() { "alive" } else { "dead" }.to_string(),
                ),
                (
                    "subdivision_exactness",
                    a.subdivision_exactness().to_string(),
                ),
                ("estimated_case_count", a.case_count_token()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                // ── extras (M-273) ──
                ("bernstein_reason", a.reason()),
                (
                    "case_classes_lower_bound",
                    a.case_classes_lower_bound.to_string(),
                ),
                ("case_count_decimal", a.case_count_decimal.to_string()),
                ("case_space_tractable", a.case_space_tractable().to_string()),
                ("ram_bytes_log2", RAM_BYTES_LOG2.to_string()),
                (
                    "control_points_matching_samples",
                    a.sites_matching.to_string(),
                ),
                ("control_points_total", a.sites_total.to_string()),
                ("corner_sites_matching", a.corner_sites_matching.to_string()),
                ("degenerate_nets", a.degenerate.to_string()),
                ("enclosure_narrowed_instances", a.narrowed.to_string()),
                ("enclosure_shrink_ppm_max", a.shrink_ppm_max.to_string()),
                ("enclosure_shrink_ppm_min", a.shrink_ppm_min.to_string()),
                (
                    "hidden_edge_sign_changes",
                    a.hidden_edge_sign_changes.to_string(),
                ),
                ("instances", INSTANCES.to_string()),
                ("is_control", (a.degree == 1).to_string()),
                ("max_exact_magnitude", a.max_magnitude.to_string()),
                ("mv_permanent", a.mv_permanent.to_string()),
                (
                    "mv_routes_agree",
                    (a.mv_polarisation == a.mv_permanent).to_string(),
                ),
                (
                    "mv_vs_trilinear",
                    format!("{:.6}", a.mv_polarisation as f64 / trilinear_mv),
                ),
                ("net_is_the_samples", a.net_is_the_samples().to_string()),
                ("net_sign_bits", a.net_sign_bits.to_string()),
                ("newton_boxes", newton_boxes_token(a.degree)),
                ("seed", a.seed.to_string()),
                (
                    "sharpness_theorems_vacuous",
                    a.sharpness_vacuous().to_string(),
                ),
                (
                    "signs_are_the_case_index",
                    a.signs_are_the_case_index().to_string(),
                ),
                (
                    "subdivision_sites_matching",
                    a.sub_sites_matching.to_string(),
                ),
                ("subdivision_sites_total", a.sub_sites_total.to_string()),
                ("trilinear_case_classes", trilinear_classes.to_string()),
                ("trilinear_case_group_order", GROUP_ORDER.to_string()),
            ]);
        }
    });
}
