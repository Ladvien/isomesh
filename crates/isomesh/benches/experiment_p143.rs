//! **P-143 — Schwarz P is `-4` and Schwarz D is `-16`, and the factor between the three is a sign.**
//!
//! Ticket: R-143. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p143
//! ```
//!
//! Writes `docs/experiments/p-143.csv`.
//!
//! # What was missing
//!
//! Before P-142 the crate had **two** fields with an analytically known Euler
//! characteristic and both were genus `<= 1`: `sphere` with
//! `expected_euler() == Some(2)` and `torus` with `Some(0)`
//! (`crates/isomesh/src/fields/mod.rs:318-320`, `:407`). The other six reference
//! fields return `None`, and `gyroid`'s `None` is documented as a standing
//! impossibility at `fields/mod.rs:23-28` and `CLAUDE.md:199-203`.
//!
//! `P-142` (`docs/experiments/p-142.csv`, 70 rows, commit `845c7a8`) closed the
//! gyroid: over a periodic-conforming box and after `wrap_seams`, `marching_cubes`
//! read `chi = -8`, `-64`, `-216` at `N = 1, 2, 3` against a prediction of `-8N^3`.
//! One field with one number is one data point, and a single arithmetic error
//! anywhere in the chain — the genus-3-per-primitive-cell premise, the lattice
//! index, the wrap, the Euler reader — would have passed it silently.
//!
//! This row is the discriminating test. Schwarz P and Schwarz D have the **same**
//! premise, the **same** wrap and the **same** reader, and predictions of `-4` and
//! `-16`. Three different numbers out of one apparatus cannot all be produced by a
//! constant, and the two ways the apparatus could be wrong — a wrong genus premise
//! or a wrong lattice index — move all three predictions together or move exactly
//! one, and either is visible.
//!
//! # C1's first half is answered NEGATIVELY and on purpose
//!
//! C1 reads *"Both are **added as reference fields** and both reproduce
//! `chi = N^3 . chi_cell` under periodic wrap."* They are **not** added. Schwarz P
//! and Schwarz D are **bench-local fixtures** (`common::tpms::NodalTpms`), nothing
//! under `crates/isomesh/src/**` changes, and `added_as_reference_field` is
//! `false` on every row of this CSV.
//!
//! The reason is priced, not vague. One new entry in `for_each_reference_field!`
//! adds **27 rows** to `crates/isomesh/golden_hashes.json` (216 = 8 fields x 9
//! algorithms x 3 resolutions) and moves `scripts/doc_facts.sh`'s gated `FIELDS`
//! and `HASHES` counts, which appear as prose phrases in twelve documents. That is
//! a repo-wide renumbering inside a measurement commit, and it is a Phase 28
//! landing ticket with the ripple priced into it. `common::tpms` refuses to
//! implement `ReferenceField` for exactly this reason and says so in its own
//! header.
//!
//! The registration's falsifier settles how this is graded: *"C1 by **either field
//! disagreeing**"*. The falsification criterion is the `chi` reproduction and
//! nothing else, so `c1_holds` is computed from the `chi` half over Schwarz P and
//! Schwarz D, and the reference-field half is **reported** in its own column
//! rather than folded into a verdict it was never given a falsifier for. This is
//! P-142 C3's treatment, one clause over.
//!
//! # C2: the per-term sign bookkeeping, done exactly and in `i128`
//!
//! C2 is the one clause that must not be inferred from the numbers it explains.
//! *"THE MECHANISM IS ASSERTED RATHER THAN ASSUMED … checked as a symbolic
//! identity, not inferred from the `chi` values it explains."*
//!
//! The body-centring shift is `(x, y, z) -> (x + pi, y + pi, z + pi)`, under which
//! `sin(t + pi) = -sin t` and `cos(t + pi) = -cos t`. So **every** trigonometric
//! factor on a shifted axis flips sign, and a term picks up `(-1)^k` where `k` is
//! the number of its factors that sit on a shifted axis:
//!
//! ```text
//! F_G = sin x cos y + sin y cos z + sin z cos x
//!       every term has TWO factors, both on shifted axes -> (-1)^2 = +1 each
//!       F_G(p + (pi,pi,pi)) = +F_G(p)                                 INVARIANT
//!
//! F_P = cos x + cos y + cos z
//!       every term has ONE factor                        -> (-1)^1 = -1 each
//!       F_P(p + (pi,pi,pi)) = -F_P(p)                                   NEGATED
//!
//! F_D = sin x sin y sin z + sin x cos y cos z
//!     + cos x sin y cos z + cos x cos y sin z
//!       every term has THREE factors, one per axis       -> (-1)^3 = -1 each
//!       F_D(p + (pi,pi,pi)) = -F_D(p)                                   NEGATED
//! ```
//!
//! A **negating** shift maps the zero set to itself but exchanges the two
//! labyrinths, so it is a symmetry of the surface and not a translation of the
//! labelled structure: it does not shrink the translational cell. An
//! **invariant** shift does. Hence the gyroid's conventional cubic cell holds two
//! primitive cells and Schwarz P's holds one, and with genus 3 and so `chi = -4`
//! in every primitive cell, `chi_cell = -4 x (primitive cells)` gives `-8` and
//! `-4`. Schwarz D's `-16` needs the other half of the same bookkeeping: the
//! **face**-centring shift `(pi, pi, 0)` flips two of `F_D`'s three factors per
//! term, `(-1)^2 = +1`, so `F_D` is invariant under all three face-centring
//! shifts, its lattice is face-centred, four primitive cells tile the cubic cell,
//! and `chi_cell = -16`.
//!
//! This harness does not take that on trust and does not read it off a table.
//! `nodal_poly` writes each nodal function as an exact polynomial in
//! `common::poly`'s eight-variable ring over `i128`, using six of the eight
//! variables as the trigonometric generators — `f0 = sin x`, `f1 = cos x`,
//! `f2 = sin y`, `f3 = cos y`, `f4 = sin z`, `f5 = cos z`, with `f6` and `f7`
//! unused. The shift is then the **exact substitution** `f_i -> -f_i` on every
//! variable of a shifted axis, and the classification is
//! `F' - F == 0` (invariant), `F' + F == 0` (negated) or neither — three
//! polynomial identities decided by `BTreeMap` cancellation in `i128`, with no
//! floating point and no tolerance anywhere.
//!
//! All **seven** non-trivial half-period shifts are classified, not just
//! `(pi,pi,pi)`, because the count of *invariant* ones is what the lattice index
//! **is**: `primitive_cells_per_cubic_cell = 1 + invariant_shifts`, and
//! `chi_per_cubic_cell = -4 x` that. The `-8 / -4 / -16` in the CSV is therefore
//! computed from the sign bookkeeping in this file, and a wrong sign anywhere
//! moves the prediction rather than hiding behind it. The expected classification,
//! in the enumeration order `(x)`, `(y)`, `(xy)`, `(z)`, `(xz)`, `(yz)`, `(xyz)`:
//!
//! | field | classes | invariant | primitive cells | `chi_cell` |
//! |---|---|---|---|---|
//! | gyroid | `neither` x6, then `invariant` | 1 | 2 | `-8` |
//! | Schwarz P | `neither` x6, then `negated` | 0 | 1 | `-4` |
//! | Schwarz D | `negated neither-free`: `- - + - + + -` | 3 | 4 | `-16` |
//!
//! Two independent cross-checks sit on top, both asserted rather than recorded:
//!
//! * **The polynomial is the module's function.** `transcription_residual`
//!   evaluates the polynomial at `(sin x, cos x, …)` on a 2197-point offset grid
//!   and compares it with `common::tpms::nodal`. A symbolic identity about a
//!   mistyped polynomial is worth nothing, and this is the only thing standing
//!   between the two transcriptions.
//! * **The exact algebra agrees with the numerics.** Every one of the 21
//!   classifications is re-derived from `common::tpms::shift_residuals`'s
//!   4913-point residual pair, and `common::tpms::body_centring_check` is run for
//!   the `(pi,pi,pi)` verdict. `shift_residuals` returns *both* the invariance and
//!   the negation residual, which is what makes either non-vacuous: a function
//!   that vanished on the grid would return two zeros and is refused.
//!
//! # Why `voxels_per_period` is ODD, and the control that proves it had to be
//!
//! `common::tpms`'s author measured 168 configurations and found that **Schwarz D
//! reads the wrong `chi` at every `voxels_per_period` divisible by 8** — `-12` at
//! 32 and 56, `-9` at 64, `-7` at 96, `+1` at 128 — and only there. Those are
//! exactly the four resolutions P-142's registration named. The mechanism is
//! `M-48`'s, not the wrap's: a multiple of 8 puts samples on the `pi/4` lattice,
//! where `F_D`'s four terms are equal in magnitude and cancel to **exactly `0.0`**
//! (at `(pi/4, pi/4, 3pi/4)`, for instance), the crossing parameter is 0 or 1, one
//! cell places coincident vertices, and the weld turns them into a pinch.
//!
//! So this row uses **33, 65 and 97** voxels per period, matching P-142's own
//! choices so the gyroid arm is comparable with `p-142.csv` row for row. But a
//! resolution choice justified only by a citation is superstition, so the failing
//! family is **run as a control**: one `degenerate_grid_control` arm at
//! `voxels_per_period = 32`, `N = 1`, on all three surfaces. It is marked
//! `is_control = true`, excluded from C1's scope by `c1_scope`, and carries
//! `non_manifold_edges` and `pinch_accounts_for_gap` so the exclusion is a named
//! mechanism with its arithmetic on the row rather than a convenience.
//!
//! `non_manifold_edges` is recorded on **every** row because the module measured
//! that in all twelve of its pinching runs `chi_measured - chi_predicted` equalled
//! `non_manifold_edges` **exactly** — 4, 4, 7, 9, 17, 32, 135 — never off by one.
//! Each pinch merges two sheets and costs exactly one from `chi`. A zero there is
//! the positive statement that no sample landed on the isosurface.
//!
//! # The one grid a prior measurement named as failing is run anyway
//!
//! The same measurement says: *"Schwarz P shows the same failure once, at
//! `N = 3, v = 33`, from ordinary floating-point cancellation rather than an
//! exact lattice."* That is a grid this row uses, on a surface C1 names, and it
//! is **kept in scope**. Dropping the one configuration a prior measurement had
//! already flagged as failing would be choosing the population after seeing the
//! answer, which is the failure this whole apparatus exists to prevent, and it
//! would be indistinguishable in the CSV from a clause that never had a hard
//! case in it.
//!
//! So: if that row pinches, `c1_holds` reads **false**, and it reads false for
//! one row out of the in-scope population rather than being argued away. What
//! goes beside it is the attribution, in two columns computed on every row:
//! `chi_pinch_corrected = chi_measured - non_manifold_edges`, and
//! `c1_holds_pinch_corrected`, the same clause re-asked with the extractor's
//! degenerate crossings put back. A `false` C1 beside a `true`
//! `c1_holds_pinch_corrected` is a precise statement — *the topological
//! prediction `N^3 . chi_cell` was exact and marching cubes' handling of a
//! sample that lands **on** the isosurface is the entire gap* — and it is a
//! statement about `M-48`, not about the oracle. The corrected verdict is
//! recorded beside C1 and never substituted for it: the registration was given
//! a falsifier and is graded by that falsifier.
//!
//! # `primitive_lattice` is not the space-group symbol
//!
//! Schwarz P's space group is `Im-3m`, which is **body-centred** — and its
//! `primitive_lattice` column reads `simple_cubic`. That is not an inconsistency,
//! it is the whole of C2. `Im-3m`'s body-centring operation *negates* `F_P`, so it
//! exchanges P's two labyrinths rather than translating the labelled structure,
//! and it therefore contributes nothing to the translation lattice. Writing `bcc`
//! there beside `chi_per_cubic_cell = -4` would put a contradiction inside one CSV
//! row: `bcc` means two primitive cells means `-8`. The column comes from
//! `Tpms::primitive_lattice`, whose own doc states this, and this harness checks
//! it against the symbolic count rather than trusting it — `primitive_cells_symbolic`
//! is derived here and asserted equal to the module's.
//!
//! # Arms
//!
//! Three surfaces x six grids x two wrap modes = **36 rows**. Both wrap modes come
//! from **one** extraction: the open reading is taken before `wrap_seams` and the
//! periodic reading after it, so the two arms differ in exactly one operation.
//!
//! | arm | `periods` | `voxels_per_period` | `resolution` | `is_control` |
//! |---|---|---|---|---|
//! | `period_sweep` | 1 | 33 | 34 | no |
//! | `period_sweep` | 2 | 33 | 67 | no |
//! | `period_sweep` | 3 | 33 | 100 | no |
//! | `resolution_sweep` | 1 | 65 | 66 | no |
//! | `resolution_sweep` | 1 | 97 | 98 | no |
//! | `degenerate_grid_control` | 1 | 32 | 33 | **yes** |
//!
//! Every row additionally carries `wrap_mode = open` or `periodic`; the `open`
//! half is a control and is marked as one.
//!
//! **One extractor, `marching_cubes`.** P-142 already measured which of the seven
//! extractors close under `wrap_seams` — the four primal ones do, the three dual
//! ones place no vertex on the boundary plane and identify zero seam pairs — over
//! 70 committed rows. Re-measuring that here would be a second answer to one
//! question, and this row's question is the *field*, not the extractor. The gyroid
//! arm at `N = 1, v = 33` is therefore directly comparable with `p-142.csv:7`,
//! which recorded `chi = -8` from `5284` vertices, `15876` edges and `10584` faces.
//!
//! # SHARE, recomputed before the numbers
//!
//! **None, and the registration says so: `SHARE: none.`** Nothing here is on an
//! extraction path and nothing here is a speedup claim, so `M-280`'s 1.45x
//! governor scatter has nothing to bite on: **every clause is an integer equality
//! or an exact polynomial identity over an enumerated population** — `chi` against
//! `N^3 . chi_cell`, and `F(p + t) -+ F(p) == 0` in `i128`. Following `P-112`, a
//! figure no clause reads is not recorded, so there is **no wall-clock column**;
//! the harness prints its own elapsed time for the operator's budget only.
//!
//! # What the wrap costs, and which readings survive it
//!
//! After `wrap_seams` the buffer is a valid *simplicial complex* and an invalid
//! *geometric* mesh: connectivity readings (`chi`, components, genus, boundary and
//! non-manifold edges) are exact and every metric reading (area, mean ratio,
//! Hausdorff, self-intersections) is nonsense. This bench reads only connectivity
//! and calls no `accuracy`, no `field_bound_report` and no `self_intersections` —
//! which would be meaningless here twice over, since a nodal function is a level
//! set and not a distance and `|grad F|` vanishes on the whole singular skeleton.
//!
//! Two independent readers run on every row and are cross-checked, and the law
//! relating them is exact and uniform over both arms rather than an expectation
//! that they agree. `validate::validate_indexed` counts `V - E + F` over the
//! buffer **as given** — `euler_characteristic` is
//! `referenced_vertices - edges + faces` at `validate.rs:174-182`, with no
//! weld. `common::tpms::euler`
//! **welds first**, at `weld::epsilon_for(cell_size)`. Welding merges the
//! coincident duplicates marching cubes leaves wherever a sample landed exactly
//! on the isosurface (`M-48`), and each such merge joins two sheets: it costs
//! exactly one from `chi` and produces exactly one non-manifold edge the raw
//! buffer did not have. So
//!
//! ```text
//! euler.chi - report.euler_characteristic
//!     == euler.non_manifold_edges - report.non_manifold_edges
//! ```
//!
//! and that identity is **asserted on every reading**. On a clean buffer both
//! sides are zero and the readers agree outright. On the wrapped arm
//! `wrap_seams` has already welded, so the weld here is a no-op, both sides are
//! zero again, and the readers agree there too — which is the arm every clause
//! is graded on. On the open arm of a grid that manufactured an exact zero
//! crossing the two differ by exactly the pinch count, and that is a
//! measurement of the weld rather than a disagreement about the surface. A
//! violation of the identity means the chi the weld removed is not the chi of
//! the sheets it merged, and that is an instrument failure, not a hypothesis
//! failure.
//!
//! # Vacuity controls
//!
//! * **The registered one: three different `chi` at `N = 1`.** The three periodic
//!   `N = 1, v = 33` rows must produce **pairwise distinct** measured `chi`, or the
//!   suite cannot distinguish a correct oracle from a constant. Distinctness, not
//!   correctness — a Schwarz D reading of `-12` is still distinct and must be
//!   recorded as `c1_holds = false` rather than aborted, or C1 could not fail.
//!   Column: `chi_measured`, and the global `chi_at_n1_triple`.
//! * **The three predictions were three numbers to begin with.** `-8`, `-4`, `-16`
//!   pairwise distinct, derived from the symbolic shift counts. Column:
//!   `chi_per_cubic_cell_symbolic`.
//! * **The lattice index discriminates.** The three symbolic
//!   `primitive_cells_symbolic` must be three different integers, or one
//!   explanation is being offered for three different numbers.
//! * **The classifier can say all three words.** `invariant`, `negated` and
//!   `neither` must each occur among the 21 classifications, or the symbolic test
//!   is a constant function of its input. Column: `shift_classes`.
//! * **No nodal polynomial is identically zero.** `F = 0` satisfies invariance and
//!   negation at once; every shift verdict would be vacuous. Column: `nodal_terms`.
//! * **Every row meshed a surface.** `faces > 0`, or `V`, `E` and `F` are all zero
//!   and `chi = 0` by vacuity rather than by measurement (`M-44`). Column: `faces`.
//! * **C1's scope is non-empty, the wrap is what put it there, and the wrap
//!   worked.** At least one in-scope row must exist; every one must have
//!   `seam_pairs_identified > 0` or `periodic` is byte-for-byte `open`; and
//!   every one must have `boundary_edges == 0` or its `chi` is being compared
//!   with a closed-surface prediction on an open surface. Columns: `c1_scope`,
//!   `seam_pairs_identified`, `boundary_edges`.
//!
//!   **This control replaces a filter, and the swap is load-bearing.** The
//!   obvious scope predicate is `MeshReport::is_closed()`, which P-142 used. It
//!   is `is_manifold() && boundary_edges == 0 && chi % 2 == 0`
//!   (`validate.rs:356`) and `is_manifold()` requires `non_manifold_edges == 0`
//!   (`:342`) — so a pinch makes a row *not closed*, and scoping on it silently
//!   removes precisely the rows that can disagree with `N^3 . chi_cell`. C1
//!   would then be quantified over a population selected for agreeing with it.
//!   Measured on this fixture, that filter drops the one Schwarz P row that
//!   disagrees and turns a false clause true. So closure is **asserted** on the
//!   property the clause actually needs, `boundary_edges == 0`, which a pinch
//!   does not satisfy away, and `is_closed` and `is_manifold` are recorded as
//!   columns instead of used as filters.
//! * **The non-wrapped control ran and is recognisable.** One `open` row per
//!   `periodic` row, and **every** open row must have `boundary_edges > 0`. It is
//!   recognised by its boundary and deliberately **not** by its `chi`: the module
//!   measured that non-wrapped Schwarz P reads `-4` at `N = 1` and `-32` at
//!   `N = 2`, hitting its own prediction by coincidence of the caps the box cuts,
//!   so a control that compared `chi` would pass the wrong arm on one field in
//!   three. What is asserted about `chi` is only that **at least one** open row
//!   disagrees, which is the weakest statement that still shows the wrap matters.
//!   Columns: `boundary_edges`, `chi_agreement`.

#![allow(clippy::too_many_lines)]

mod common;

use std::f64::consts::PI;
use std::time::Instant;

use isomesh::MeshBuffer;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{ValidateConfig, validate_indexed};

use common::poly::{self, Poly};
use common::tpms::{self, EulerCount, NodalTpms, Tpms};

/// The `wrap_mode` value for a seam-identified reading.
const WRAP_PERIODIC: &str = "periodic";
/// The `wrap_mode` value for the non-wrapped control reading.
const WRAP_OPEN: &str = "open";

/// The single extractor, named in the CSV. See the header for why one.
const EXTRACTOR: &str = "marching_cubes";

/// C1's registered range of periods, and the resolution-stability arm.
const ARM_PERIOD: &str = "period_sweep";
/// One surface on three grids, at fixed `N`.
const ARM_RESOLUTION: &str = "resolution_sweep";
/// The multiple-of-8 family the odd choice exists to avoid.
const ARM_DEGENERATE: &str = "degenerate_grid_control";

/// `voxels_per_period` of the grid every non-control arm shares.
const BASE_VOXELS: u32 = 33;

/// The six trigonometric generators as variable indices in `common::poly`'s
/// eight-variable ring: axis `a` owns `sin` at `TRIG[a][0]` and `cos` at
/// `TRIG[a][1]`. `f6` and `f7` are unused and stay at exponent zero.
const TRIG: [[usize; 2]; 3] = [[0, 1], [2, 3], [4, 5]];

/// Index of `(pi, pi, pi)` in [`half_shifts`]: all three bits set, so last.
const BODY_CENTRING: usize = 6;

// ── the symbolic half: C2 ───────────────────────────────────────────────────

/// The nodal function of `kind` as an exact polynomial in the six
/// trigonometric generators.
///
/// One period is `2*pi` per axis, matching `common::tpms::nodal`, against which
/// this transcription is checked by [`transcription_residual`].
fn nodal_poly(kind: Tpms) -> Poly {
    let s = |a: usize| Poly::var(TRIG[a][0]);
    let c = |a: usize| Poly::var(TRIG[a][1]);
    match kind {
        // sin x cos y + sin y cos z + sin z cos x
        Tpms::Gyroid => s(0).mul(&c(1)).add(&s(1).mul(&c(2))).add(&s(2).mul(&c(0))),
        // cos x + cos y + cos z
        Tpms::SchwarzP => c(0).add(&c(1)).add(&c(2)),
        // sin x sin y sin z + sin x cos y cos z + cos x sin y cos z
        //   + cos x cos y sin z
        Tpms::SchwarzD => s(0)
            .mul(&s(1))
            .mul(&s(2))
            .add(&s(0).mul(&c(1)).mul(&c(2)))
            .add(&c(0).mul(&s(1)).mul(&c(2)))
            .add(&c(0).mul(&c(1)).mul(&s(2))),
    }
}

/// The largest `|poly(sin x, cos x, …) - tpms::nodal(kind, p)|` over a
/// deterministic 13³ grid offset off every symmetry plane.
///
/// The whole of C2 is an identity about `nodal_poly`; if `nodal_poly` is not
/// the function `common::tpms` samples and `wrap_seams` extracts, the identity
/// is about a different surface and explains nothing. This is the only thing
/// tying the two transcriptions together, so it is asserted, not recorded.
fn transcription_residual(kind: Tpms, base: &Poly) -> f64 {
    const GRID: u32 = 13;
    let offsets = [1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0];
    let coord = |i: u32, axis: usize| 2.0 * PI * (f64::from(i) + offsets[axis]) / f64::from(GRID);
    let mut worst = 0.0_f64;
    for i in 0..GRID {
        let x = coord(i, 0);
        for j in 0..GRID {
            let y = coord(j, 1);
            for k in 0..GRID {
                let p = [x, y, coord(k, 2)];
                let mut args = [0.0_f64; poly::VARS];
                for (axis, pair) in TRIG.iter().enumerate() {
                    args[pair[0]] = p[axis].sin();
                    args[pair[1]] = p[axis].cos();
                }
                worst = worst.max((base.eval_f64(&args) - tpms::nodal(kind, p)).abs());
            }
        }
    }
    worst
}

/// How a nodal function responds to a half-period shift.
///
/// `Neither` is a real outcome and its presence is what makes the other two
/// informative: a classifier that could only answer `Invariant` or `Negated`
/// would be a coin, not a measurement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ShiftClass {
    /// `F(p + t) = F(p)`: a translation, and it shrinks the primitive cell.
    Invariant,
    /// `F(p + t) = -F(p)`: the two labyrinths swap, the lattice does not grow.
    Negated,
    /// Neither identity holds.
    Neither,
}

impl ShiftClass {
    /// The CSV word.
    fn name(self) -> &'static str {
        match self {
            ShiftClass::Invariant => "invariant",
            ShiftClass::Negated => "negated",
            ShiftClass::Neither => "neither",
        }
    }
}

/// The seven non-trivial half-period shifts, enumerated from their three bits
/// rather than listed: bit `a` set means axis `a` moves by `pi`.
///
/// Enumerated, so the population is a property of `{0, pi}^3` and not of
/// somebody's list — which is what lets the count of invariant ones be read as
/// the lattice index. Order is `(x)`, `(y)`, `(xy)`, `(z)`, `(xz)`, `(yz)`,
/// `(xyz)`, and index [`BODY_CENTRING`] is the last.
fn half_shifts() -> [[bool; 3]; 7] {
    std::array::from_fn(|i| {
        let bits = i + 1;
        std::array::from_fn(|axis| bits & (1 << axis) != 0)
    })
}

/// Classify `base` under the half-period shift naming `shifted` axes, exactly.
///
/// `sin(t + pi) = -sin t` and `cos(t + pi) = -cos t`, so the shift is the exact
/// substitution `f_i -> -f_i` on both generators of every shifted axis. The
/// verdict is then two polynomial identities decided by cancellation in `i128`.
/// No floating point, no tolerance.
fn classify_symbolic(base: &Poly, shifted: [bool; 3]) -> ShiftClass {
    let mut moved = base.clone();
    for (axis, on) in shifted.iter().enumerate() {
        if *on {
            for var in TRIG[axis] {
                moved = moved.substitute(var, &Poly::var(var).scale(-1));
            }
        }
    }
    let invariant = moved.sub(base).is_zero();
    let negated = moved.add(base).is_zero();
    assert!(
        !(invariant && negated),
        "VOID: a nodal polynomial that is both invariant and negated under \
         {shifted:?} is identically zero, so every shift verdict in this row is \
         vacuous and the lattice index counts nothing"
    );
    if invariant {
        ShiftClass::Invariant
    } else if negated {
        ShiftClass::Negated
    } else {
        ShiftClass::Neither
    }
}

/// The same classification read off `common::tpms::shift_residuals`, with both
/// residuals returned so the caller can print what decided it.
///
/// Returns `(class, invariance residual, negation residual)`.
fn classify_numeric(kind: Tpms, shifted: [bool; 3]) -> (ShiftClass, f64, f64) {
    let shift: [f64; 3] = std::array::from_fn(|axis| if shifted[axis] { PI } else { 0.0 });
    let (symmetric, antisymmetric) = tpms::shift_residuals(kind, shift);
    let invariant = symmetric <= tpms::SHIFT_RESIDUAL_TOLERANCE;
    let negated = antisymmetric <= tpms::SHIFT_RESIDUAL_TOLERANCE;
    assert!(
        !(invariant && negated),
        "VOID: {} has both residuals under tolerance at {shifted:?} \
         ({symmetric:.3e}, {antisymmetric:.3e}), which means the nodal function \
         vanished on the whole check grid and the shift identity is being read \
         off zeros",
        kind.name()
    );
    let class = if invariant {
        ShiftClass::Invariant
    } else if negated {
        ShiftClass::Negated
    } else {
        ShiftClass::Neither
    };
    (class, symmetric, antisymmetric)
}

/// C2's registered claim for one surface, written from the registration prose
/// rather than read from `common::tpms::Tpms::body_centring_invariant`.
///
/// The module's table is then checked **against** this instead of being trusted
/// for it, so a wrong entry in the module cannot grade itself correct.
fn c2_claim(kind: Tpms) -> ShiftClass {
    match kind {
        Tpms::Gyroid => ShiftClass::Invariant,
        Tpms::SchwarzP | Tpms::SchwarzD => ShiftClass::Negated,
    }
}

/// One surface's symbolic half: the exact shift classification and everything
/// derived from it.
#[derive(Clone, Debug)]
struct Symbolic {
    /// Which surface.
    kind: Tpms,
    /// `Display` of the nodal polynomial over `f0..f5`; CSV-safe by construction.
    expression: String,
    /// Non-zero terms: 3, 3, 4.
    terms: usize,
    /// The seven verdicts, in [`half_shifts`] order.
    classes: [ShiftClass; 7],
    /// How many of the seven are `Invariant`.
    invariant: usize,
    /// How many are `Negated`.
    negated: usize,
    /// How many are `Neither`.
    neither: usize,
    /// `1 + invariant` — the lattice index, derived rather than transcribed.
    primitive_cells: i64,
    /// `-4 * primitive_cells`, genus 3 per primitive cell being the premise.
    chi_per_cubic_cell: i64,
    /// Largest disagreement with `common::tpms::nodal`.
    transcription_residual: f64,
    /// `common::tpms::body_centring_check`'s residual for the claimed relation.
    body_centring_residual: f64,
    /// `common::tpms::body_centring_check`'s verdict: the claimed relation held
    /// and the opposite one did not.
    body_centring_numeric_ok: bool,
}

impl Symbolic {
    /// The `(pi, pi, pi)` verdict — the `body_centering_invariance` column.
    fn body_centring(&self) -> ShiftClass {
        self.classes[BODY_CENTRING]
    }

    /// `invariant|negated|neither|…`, in [`half_shifts`] order, pipe-joined
    /// because `Run::record` refuses a comma.
    fn classes_column(&self) -> String {
        self.classes
            .iter()
            .map(|c| c.name())
            .collect::<Vec<_>>()
            .join("|")
    }

    /// Does this surface's symbolic half discharge C2?
    ///
    /// All five halves, because C2 claims a mechanism and a consequence: the
    /// `(pi,pi,pi)` verdict is the registered one, the polynomial really is the
    /// module's function, the lattice index the shift count implies matches the
    /// module's, the `chi` that index implies matches the module's, and the
    /// exact algebra survives its numerical cross-check.
    fn c2_holds(&self) -> bool {
        self.body_centring() == c2_claim(self.kind)
            && self.transcription_residual <= tpms::SHIFT_RESIDUAL_TOLERANCE
            && self.primitive_cells == self.kind.primitive_cells_per_cubic_cell()
            && self.chi_per_cubic_cell == self.kind.chi_per_cubic_cell()
            && self.body_centring_numeric_ok
    }
}

/// Classify all seven half-period shifts for one surface, exactly, and
/// cross-check every verdict against the numerical residuals.
fn symbolic(kind: Tpms) -> Symbolic {
    let base = nodal_poly(kind);
    assert!(
        !base.is_zero(),
        "VOID: {}'s nodal polynomial came out identically zero, so it is \
         invariant and negated under every shift at once and the whole lattice \
         index is counted from nothing",
        kind.name()
    );

    let residual = transcription_residual(kind, &base);
    assert!(
        residual <= tpms::SHIFT_RESIDUAL_TOLERANCE,
        "{}: the bench-local polynomial `{base}` disagrees with \
         common::tpms::nodal by {residual:.3e}, so every symbolic identity below \
         is about a different function than the one that was extracted, and the \
         instrument is what failed rather than the hypothesis",
        kind.name()
    );

    let mut classes = [ShiftClass::Neither; 7];
    for (slot, shifted) in half_shifts().into_iter().enumerate() {
        let exact = classify_symbolic(&base, shifted);
        let (numeric, symmetric, antisymmetric) = classify_numeric(kind, shifted);
        assert_eq!(
            exact,
            numeric,
            "{}: exact i128 algebra classifies the shift {shifted:?} as `{}` \
             while the 4913-point residual pair ({symmetric:.3e}, \
             {antisymmetric:.3e}) says `{}` -- the two instruments disagree, so \
             neither reading is usable",
            kind.name(),
            exact.name(),
            numeric.name()
        );
        classes[slot] = exact;
    }

    let count = |want: ShiftClass| classes.iter().filter(|c| **c == want).count();
    let invariant = count(ShiftClass::Invariant);
    let primitive_cells = 1 + i64::try_from(invariant).expect("seven shifts fit i64");
    let (body_centring_residual, body_centring_numeric_ok) = tpms::body_centring_check(kind);

    Symbolic {
        kind,
        expression: base.to_string(),
        terms: base.terms(),
        classes,
        invariant,
        negated: count(ShiftClass::Negated),
        neither: count(ShiftClass::Neither),
        primitive_cells,
        chi_per_cubic_cell: -4 * primitive_cells,
        transcription_residual: residual,
        body_centring_residual,
        body_centring_numeric_ok,
    }
}

// ── the extraction half: C1 ─────────────────────────────────────────────────

/// One `(periods, voxels_per_period)` grid and the arm it answers for.
#[derive(Clone, Copy, Debug)]
struct Config {
    /// `N` — periods per axis, so the box is `N³` conventional cubic cells.
    periods: u32,
    /// Cells spanning one period.
    voxels: u32,
    /// Which arm this grid belongs to.
    arm: &'static str,
}

impl Config {
    /// Is this the `pi/4`-lattice family `common::tpms` measured Schwarz D
    /// failing on? A multiple of 8 puts samples where `F_D` is exactly `0.0`.
    fn degenerate(&self) -> bool {
        self.voxels.is_multiple_of(8)
    }
}

/// Three periods, two extra resolutions, and the failing family as a control.
const CONFIGS: [Config; 6] = [
    Config {
        periods: 1,
        voxels: BASE_VOXELS,
        arm: ARM_PERIOD,
    },
    Config {
        periods: 2,
        voxels: BASE_VOXELS,
        arm: ARM_PERIOD,
    },
    Config {
        periods: 3,
        voxels: BASE_VOXELS,
        arm: ARM_PERIOD,
    },
    Config {
        periods: 1,
        voxels: 65,
        arm: ARM_RESOLUTION,
    },
    Config {
        periods: 1,
        voxels: 97,
        arm: ARM_RESOLUTION,
    },
    Config {
        periods: 1,
        voxels: 32,
        arm: ARM_DEGENERATE,
    },
];

/// One CSV row: one surface's reading of one grid under one wrap mode.
#[derive(Clone, Debug)]
struct Row {
    /// Which surface.
    kind: Tpms,
    /// Which grid.
    config: Config,
    /// Samples per axis, `voxels * periods + 1`; the `resolution` column.
    samples: u32,
    /// `(voxels * periods)³`.
    cells: u64,
    /// `2*pi / voxels`.
    cell_size: f64,
    /// `weld::epsilon_for(cell_size)`, the seam and weld tolerance.
    weld_tolerance: f64,
    /// `periodic` or `open`.
    wrap_mode: &'static str,
    /// Vertex pairs `wrap_seams` identified; `0` on every open row by definition.
    seam_pairs: u64,
    /// `N^3 * chi_cell`.
    chi_predicted: i64,
    /// `common::tpms::euler`'s reading.
    counted: EulerCount,
    /// `MeshReport::euler_characteristic`, the independent reader.
    validate_chi: i64,
    /// `MeshReport::non_manifold_edges` — the **unwelded** count, whose
    /// difference from the welded one is exactly the chi the weld removed.
    validate_non_manifold_edges: u64,
    /// `MeshReport::genus`; `None` where the crate declines to name one.
    genus: Option<i64>,
    /// `MeshReport::components`.
    components: u64,
    /// `MeshReport::is_closed()`.
    closed: bool,
    /// `MeshReport::is_manifold()`.
    manifold: bool,
    /// Buffer vertices at the time of this reading.
    mesh_vertices: usize,
    /// Buffer triangles at the time of this reading.
    mesh_triangles: usize,
}

impl Row {
    /// Is this the wrapped arm?
    fn periodic(&self) -> bool {
        self.wrap_mode == WRAP_PERIODIC
    }

    /// `chi_measured == chi_predicted` for this row.
    fn chi_agreement(&self) -> bool {
        self.counted.chi == self.chi_predicted
    }

    /// `chi_measured - chi_predicted`, the pinching diagnostic's left side.
    fn chi_gap(&self) -> i64 {
        self.counted.chi - self.chi_predicted
    }

    /// `chi_measured` with every pinch put back.
    ///
    /// A pinch merges two sheets and costs exactly one from `chi`, so this is
    /// the Euler characteristic the extraction would have read had no sample
    /// landed on the isosurface. On a clean row it is `chi_measured` unchanged.
    /// It is recorded so that a `chi` gap can be attributed rather than merely
    /// reported: if this equals `chi_predicted` on a disagreeing row, the
    /// topological prediction was exact and the extractor's degenerate-crossing
    /// handling is the whole of the discrepancy.
    fn chi_pinch_corrected(&self) -> i64 {
        self.counted.chi
            - i64::try_from(self.counted.non_manifold_edges)
                .expect("a non-manifold edge count fits i64")
    }

    /// Does the whole `chi` gap sit on pinched edges?
    ///
    /// `common::tpms` measured `chi_measured - chi_predicted ==
    /// non_manifold_edges` in all twelve of its pinching runs. On a clean row
    /// both sides are zero and this is trivially true, which is the point: the
    /// column is a statement about where any discrepancy went, and a `false`
    /// says the gap is something other than degenerate crossings.
    fn pinch_accounts_for_gap(&self) -> bool {
        self.chi_pinch_corrected() == self.chi_predicted
    }

    /// Is this row a wrapped reading on a grid outside the `pi/4` family — the
    /// population C1 is quantified over?
    ///
    /// **`MeshReport::is_closed()` is deliberately not part of this**, and the
    /// reason is a back door that has to stay shut. `is_closed()` is
    /// `is_manifold() && boundary_edges == 0 && chi % 2 == 0`
    /// (`validate.rs:356`), and `is_manifold()` requires
    /// `non_manifold_edges == 0` (`:342`). A pinch is a non-manifold edge. So
    /// scoping on `is_closed()` would drop exactly the rows where a sample
    /// landed on the isosurface — the only rows that can disagree with
    /// `N^3 . chi_cell` at all — and C1 would be quantified over a population
    /// selected for agreeing with it. Measured on this fixture: it removes the
    /// one Schwarz P row that disagrees and turns a false clause true.
    ///
    /// What the scope does require is what the clause actually needs: the arm
    /// is the wrapped one, and the grid is not the `pi/4` family whose failure
    /// `common::tpms` measured before this harness existed. That the wrap
    /// really did close the surface is asserted as a vacuity control on
    /// `boundary_edges`, which a pinch does not satisfy away.
    fn in_c1_scope(&self) -> bool {
        self.periodic() && !self.config.degenerate()
    }

    /// Is this row inside C1's **verdict**?
    ///
    /// C1 names Schwarz P and Schwarz D — *"Both are added as reference fields
    /// and both reproduce …"* — and its falsifier names the same two: *"C1 by
    /// either field disagreeing."* The gyroid is in scope as the third leg of
    /// the registered vacuity control and is graded separately in
    /// `all_three_agree_in_scope`, so that P-142's field cannot carry or sink a
    /// clause that was never about it.
    fn in_c1_verdict(&self) -> bool {
        self.in_c1_scope() && self.kind != Tpms::Gyroid
    }

    /// Why this row is or is not in C1's scope.
    fn scope(&self) -> &'static str {
        if !self.periodic() {
            "control_open_arm"
        } else if self.config.degenerate() {
            "control_degenerate_grid"
        } else if self.kind == Tpms::Gyroid {
            "in_scope_third_leg"
        } else {
            "in_scope"
        }
    }
}

/// Extract every surface on every grid, reading `chi` before and after the wrap.
fn sweep() -> Vec<Row> {
    let mut rows: Vec<Row> = Vec::with_capacity(Tpms::ALL.len() * CONFIGS.len() * 2);
    let mut mc = MarchingCubes::<f64>::new();

    for kind in Tpms::ALL {
        for config in CONFIGS {
            let started = Instant::now();
            let field = NodalTpms::new(kind, config.periods);
            let (lo, hi) = field.domain();
            let (shape, origin, cell_size) = field.periodic_grid(config.voxels);
            let samples = config.voxels * config.periods + 1;
            let cells = u64::from(config.voxels * config.periods).pow(3);
            let tol = isomesh::weld::epsilon_for(cell_size);
            let cfg = ValidateConfig::from_cell_size(cell_size)
                .expect("a periodic cell size is finite and positive");
            let chi_predicted = field.chi_predicted();

            let mut mesh = MeshBuffer::<f64>::new();
            mc.extract(&field, &shape, origin, cell_size, &mut mesh)
                .expect("the periodic grid holds at least two samples on every axis");

            // Both arms from one extraction: the open reading is taken with the
            // wrap withheld, so the arms differ in one operation and no more.
            let push = |rows: &mut Vec<Row>,
                        wrap_mode: &'static str,
                        seam_pairs: u64,
                        mesh: &MeshBuffer<f64>| {
                let counted = tpms::euler(&mesh.positions, &mesh.indices, tol);
                let report = validate_indexed(&mesh.positions, &mesh.indices, &cfg);
                // The two readers are one instrument checked against itself,
                // and the law relating them is exact and uniform over both
                // arms. `validate_indexed` counts `V - E + F` over the buffer
                // **as given**; `common::tpms::euler` **welds first**. So the
                // chi the weld removed is exactly the chi of the sheets it
                // merged, and each such merge shows up as one non-manifold edge
                // that the raw buffer did not have:
                //
                //     counted.chi - report.chi
                //         == counted.non_manifold_edges - report.non_manifold_edges
                //
                // On a buffer with no coincident duplicates both sides are zero
                // and the readers agree outright. On the wrapped arm
                // `wrap_seams` already welded, so the weld here is a no-op, both
                // sides are zero again, and the readers agree there too — which
                // is what matters, because every clause is graded on that arm.
                // On the open arm of a grid that manufactured an exact zero
                // crossing (`M-48`) the two readers differ by exactly the pinch
                // count, and that is a measurement rather than a disagreement.
                let reader_gap = counted.chi - report.euler_characteristic;
                let weld_pinches = i64::try_from(counted.non_manifold_edges)
                    .expect("a non-manifold edge count fits i64")
                    - i64::try_from(report.non_manifold_edges)
                        .expect("a non-manifold edge count fits i64");
                assert_eq!(
                    reader_gap,
                    weld_pinches,
                    "{} {wrap_mode} N={} v={}: common::tpms::euler reads chi {} \
                     and validate_indexed reads {}, a gap of {reader_gap}, while \
                     the weld created {weld_pinches} non-manifold edges -- the \
                     chi the weld removed is not the chi of the sheets it \
                     merged, so the instrument is what failed, not the hypothesis",
                    kind.name(),
                    config.periods,
                    config.voxels,
                    counted.chi,
                    report.euler_characteristic
                );
                rows.push(Row {
                    kind,
                    config,
                    samples,
                    cells,
                    cell_size,
                    weld_tolerance: tol,
                    wrap_mode,
                    seam_pairs,
                    chi_predicted,
                    counted,
                    validate_chi: report.euler_characteristic,
                    validate_non_manifold_edges: report.non_manifold_edges,
                    genus: report.genus,
                    components: report.components,
                    closed: report.is_closed(),
                    manifold: report.is_manifold(),
                    mesh_vertices: mesh.vertex_count(),
                    mesh_triangles: mesh.triangle_count(),
                });
            };

            push(&mut rows, WRAP_OPEN, 0, &mesh);
            let seam_pairs = tpms::wrap_seams(&mut mesh, lo, hi, tol);
            push(&mut rows, WRAP_PERIODIC, seam_pairs, &mesh);

            let wrapped = rows.last().expect("the periodic row was just pushed");
            println!(
                "-- {:<9} N={} v={:<3} chi {:>5} against {:>5}  nm_edges {:<4} \
                 seams {:<5} in {:.2} s",
                kind.name(),
                config.periods,
                config.voxels,
                wrapped.counted.chi,
                chi_predicted,
                wrapped.counted.non_manifold_edges,
                seam_pairs,
                started.elapsed().as_secs_f64()
            );
        }
    }

    rows
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-143");

    common::experiment::run(prereg, |run| {
        let started = Instant::now();

        // The symbolic half runs first: it is what `chi_predicted` is derived
        // from, so a broken derivation must stop the harness before a single
        // extraction pretends to confirm it.
        let symbolics: [Symbolic; 3] = Tpms::ALL.map(symbolic);
        for sym in &symbolics {
            println!(
                "-- {:<9} F = {}  ->  {} : {} invariant, {} negated, {} neither \
                 -> {} primitive cells -> chi_cell {}",
                sym.kind.name(),
                sym.expression,
                sym.classes_column(),
                sym.invariant,
                sym.negated,
                sym.neither,
                sym.primitive_cells,
                sym.chi_per_cubic_cell
            );
        }

        let rows = sweep();
        println!(
            "-- {} rows over {} surfaces x {} grids x 2 wrap modes in {:.1} s",
            rows.len(),
            Tpms::ALL.len(),
            CONFIGS.len(),
            started.elapsed().as_secs_f64()
        );

        // ── vacuity controls, before any record ─────────────────────────────
        //
        // 1. No nodal polynomial is degenerate, the lattice index discriminates,
        //    and the classifier can say all three of its words. (The
        //    identically-zero and instrument-agreement guards already fired
        //    inside `symbolic`.)
        let mut cells_seen: Vec<i64> = symbolics.iter().map(|s| s.primitive_cells).collect();
        cells_seen.sort_unstable();
        cells_seen.dedup();
        assert_eq!(
            cells_seen.len(),
            3,
            "VOID: the three surfaces produced {} distinct lattice indices from \
             their shift counts, not three, so one explanation is being offered \
             for what are supposed to be three different chi values",
            cells_seen.len()
        );

        let mut predicted_cell: Vec<i64> = symbolics.iter().map(|s| s.chi_per_cubic_cell).collect();
        predicted_cell.sort_unstable();
        predicted_cell.dedup();
        assert_eq!(
            predicted_cell.len(),
            3,
            "VOID: the three predictions per cubic cell are not three different \
             numbers, so C1 could be satisfied by a constant before any surface \
             was extracted"
        );

        for want in [
            ShiftClass::Invariant,
            ShiftClass::Negated,
            ShiftClass::Neither,
        ] {
            let seen = symbolics
                .iter()
                .flat_map(|s| s.classes.iter())
                .filter(|c| **c == want)
                .count();
            assert!(
                seen > 0,
                "VOID: not one of the 21 shift classifications came out `{}`, so \
                 the symbolic test is a constant function of its input and the \
                 verdicts it did produce establish nothing",
                want.name()
            );
        }

        // 2. Every row meshed a surface that exists.
        for row in &rows {
            assert!(
                row.counted.faces > 0,
                "VOID: {} {} at N={} v={} meshed to zero faces, so its V, E and F \
                 are all zero and its chi of {} is vacuous rather than measured \
                 (M-44)",
                row.kind.name(),
                row.wrap_mode,
                row.config.periods,
                row.config.voxels,
                row.counted.chi
            );
        }

        // 3. C1's scope is non-empty, and the wrap is what put it there.
        let in_scope: Vec<&Row> = rows.iter().filter(|r| r.in_c1_scope()).collect();
        assert!(
            !in_scope.is_empty(),
            "VOID: no wrapped reading survives on a non-degenerate grid, so C1's \
             population is empty and C1 would hold for the wrong reason (M-44)"
        );
        for row in &in_scope {
            assert!(
                row.seam_pairs > 0,
                "VOID: the wrapped {} at N={} v={} had wrap_seams identify zero \
                 vertex pairs, so the periodic arm is byte-for-byte the open arm \
                 and this row's agreement says nothing about periodicity",
                row.kind.name(),
                row.config.periods,
                row.config.voxels
            );
        }
        for row in &in_scope {
            assert!(
                row.counted.boundary_edges == 0,
                "VOID: the wrapped {} at N={} v={} still has {} boundary edges, \
                 so the seam identification did not close it and its chi of {} \
                 is being compared with a closed-surface prediction of {} on an \
                 open surface -- the comparison is meaningless rather than \
                 falsifying",
                row.kind.name(),
                row.config.periods,
                row.config.voxels,
                row.counted.boundary_edges,
                row.counted.chi,
                row.chi_predicted
            );
        }

        // 4. The registered vacuity control: three DIFFERENT chi at N = 1.
        //    Distinctness, never correctness -- a Schwarz D reading of -12 is
        //    still distinct and must be recorded as c1_holds = false rather than
        //    aborted here, or C1 could not fail.
        let at_n1: Vec<&Row> = rows
            .iter()
            .filter(|r| r.periodic() && r.config.periods == 1 && r.config.voxels == BASE_VOXELS)
            .collect();
        assert_eq!(
            at_n1.len(),
            3,
            "VOID: the N=1 comparison needs one wrapped row per surface and has \
             {}, so the registered control cannot be evaluated",
            at_n1.len()
        );
        for (i, a) in at_n1.iter().enumerate() {
            for b in &at_n1[i + 1..] {
                assert_ne!(
                    a.counted.chi,
                    b.counted.chi,
                    "VOID: {} and {} both measured chi = {} at N = 1, so the \
                     three fields do not produce three different values and this \
                     suite cannot distinguish a correct oracle from a constant",
                    a.kind.name(),
                    b.kind.name(),
                    a.counted.chi
                );
            }
        }
        let chi_at_n1_triple = at_n1
            .iter()
            .map(|r| r.counted.chi.to_string())
            .collect::<Vec<_>>()
            .join("|");

        // 5. The non-wrapped control ran, is recognisable by its boundary, and
        //    at least one of its rows fails its prediction. Recognisable by
        //    BOUNDARY and not by chi: common::tpms measured non-wrapped Schwarz P
        //    hitting -4 at N=1 and -32 at N=2 by coincidence of the caps the box
        //    cuts, so a chi-only control would pass the wrong arm on one field.
        let open: Vec<&Row> = rows.iter().filter(|r| !r.periodic()).collect();
        let periodic: Vec<&Row> = rows.iter().filter(|r| r.periodic()).collect();
        assert!(
            !open.is_empty() && open.len() == periodic.len(),
            "VOID: the non-wrapped control arm must be run once per wrapped row \
             ({} open against {} periodic), or the experiment has not shown that \
             periodicity is what matters",
            open.len(),
            periodic.len()
        );
        for row in &open {
            assert!(
                row.counted.boundary_edges > 0,
                "VOID: the non-wrapped {} at N={} v={} has zero boundary edges, \
                 so it is not recognisable as non-wrapped and is not a control",
                row.kind.name(),
                row.config.periods,
                row.config.voxels
            );
        }
        assert!(
            open.iter().any(|r| !r.chi_agreement()),
            "VOID: every non-wrapped row reproduced its own N^3.chi_cell without \
             any seam identification, so this fixture does not discriminate the \
             wrap at all"
        );

        // ── verdicts ────────────────────────────────────────────────────────
        //
        // C1 is graded by its own falsifier -- "C1 by either field disagreeing"
        // -- over Schwarz P and Schwarz D, the two surfaces the clause names.
        // The gyroid is the third leg of the vacuity control and is graded
        // separately in `all_three_agree_in_scope`.
        let registered: Vec<&Row> = rows.iter().filter(|r| r.in_c1_verdict()).collect();
        assert!(
            !registered.is_empty(),
            "VOID: no in-scope row belongs to Schwarz P or Schwarz D, so C1's own \
             population is empty even though the gyroid's is not"
        );
        let c1 = registered.iter().all(|r| r.chi_agreement());
        let all_three = in_scope.iter().all(|r| r.chi_agreement());
        // The same clause with the extractor's degenerate crossings accounted
        // for. A `false` C1 beside a `true` here says the topological oracle was
        // exact and marching cubes' handling of a sample that lands **on** the
        // isosurface is the whole of the gap -- which is a statement about
        // `M-48`, not about `N^3 . chi_cell`. Recorded separately and never
        // substituted for C1: the registration was given a falsifier and it is
        // graded by it.
        let c1_corrected = registered.iter().all(|r| r.pinch_accounts_for_gap());
        let c2 = symbolics.iter().all(Symbolic::c2_holds);

        // The module's own table is checked against the registration's claim
        // rather than trusted for it.
        for sym in &symbolics {
            assert_eq!(
                sym.kind.body_centring_invariant(),
                c2_claim(sym.kind) == ShiftClass::Invariant,
                "{}: common::tpms::Tpms::body_centring_invariant disagrees with \
                 P-143's registered claim, so the module and the registration are \
                 describing different surfaces",
                sym.kind.name()
            );
        }

        let verdicts = Verdicts {
            c1,
            c1_corrected,
            c2,
            all_three,
            chi_at_n1_triple,
        };
        println!(
            "\n-- C1 (Schwarz P and D reproduce N^3.chi_cell on {} in-scope rows): \
             {c1}\n-- the same with pinches accounted for: {c1_corrected}\n\
             -- all three surfaces in scope: {all_three}\n\
             -- C2 (body-centring mechanism, exact in i128): {c2}\n\
             -- chi at N=1, in Tpms::ALL order: {}",
            registered.len(),
            verdicts.chi_at_n1_triple
        );

        // ── rows ────────────────────────────────────────────────────────────
        for row in &rows {
            let sym = symbolics
                .iter()
                .find(|s| s.kind == row.kind)
                .expect("every surface has a symbolic half");
            row_out(run, row, sym, &verdicts);
        }
    });
}

/// The clause verdicts and global figures that are identical on every row.
#[derive(Clone, Debug)]
struct Verdicts {
    /// C1 over Schwarz P and Schwarz D, graded by the registration's falsifier.
    c1: bool,
    /// The same clause with the extractor's degenerate crossings accounted for.
    c1_corrected: bool,
    /// C2, the symbolic body-centring mechanism.
    c2: bool,
    /// Whether the gyroid's in-scope rows agreed as well.
    all_three: bool,
    /// The three measured `chi` at `N = 1`, pipe-joined in `Tpms::ALL` order.
    chi_at_n1_triple: String,
}

/// Write one row, registered columns first and in registration order.
///
/// `c1_holds` and `c2_holds` are **global** clause verdicts and carry the same
/// value on every row, as the authoring contract allows: C1 quantifies over the
/// in-scope population and C2 is an identity per surface with no per-row
/// content at all. Per-row agreement lives in `chi_agreement`.
fn row_out(run: &mut common::experiment::Run, row: &Row, sym: &Symbolic, v: &Verdicts) {
    run.record(&[
        ("field", row.kind.name().to_string()),
        ("space_group", row.kind.space_group().to_string()),
        (
            "primitive_lattice",
            row.kind.primitive_lattice().to_string(),
        ),
        (
            "chi_per_cubic_cell",
            row.kind.chi_per_cubic_cell().to_string(),
        ),
        ("chi_predicted", row.chi_predicted.to_string()),
        ("chi_measured", row.counted.chi.to_string()),
        ("periods", row.config.periods.to_string()),
        ("resolution", row.samples.to_string()),
        (
            "body_centering_invariance",
            sym.body_centring().name().to_string(),
        ),
        ("c1_holds", v.c1.to_string()),
        ("c2_holds", v.c2.to_string()),
        // ── extras (M-273) ──
        ("added_as_reference_field", false.to_string()),
        ("all_three_agree_in_scope", v.all_three.to_string()),
        ("arm", row.config.arm.to_string()),
        (
            "body_centering_numeric_ok",
            sym.body_centring_numeric_ok.to_string(),
        ),
        (
            "body_centering_residual",
            format!("{:.3e}", sym.body_centring_residual),
        ),
        ("boundary_edges", row.counted.boundary_edges.to_string()),
        ("c1_scope", row.scope().to_string()),
        ("cell_size", format!("{:.9}", row.cell_size)),
        ("cells", row.cells.to_string()),
        ("chi_agreement", row.chi_agreement().to_string()),
        ("chi_at_n1_triple", v.chi_at_n1_triple.clone()),
        ("chi_gap", row.chi_gap().to_string()),
        ("c1_holds_pinch_corrected", v.c1_corrected.to_string()),
        ("chi_pinch_corrected", row.chi_pinch_corrected().to_string()),
        (
            "chi_per_cubic_cell_symbolic",
            sym.chi_per_cubic_cell.to_string(),
        ),
        ("components", row.components.to_string()),
        ("edges", row.counted.edges.to_string()),
        ("extractor", EXTRACTOR.to_string()),
        ("faces", row.counted.faces.to_string()),
        (
            "genus_measured",
            row.genus
                .map_or_else(|| String::from("none"), |g| g.to_string()),
        ),
        ("in_c1_verdict", row.in_c1_verdict().to_string()),
        ("is_closed", row.closed.to_string()),
        ("is_control", (!row.in_c1_scope()).to_string()),
        ("is_manifold", row.manifold.to_string()),
        ("mesh_triangles", row.mesh_triangles.to_string()),
        ("mesh_vertices", row.mesh_vertices.to_string()),
        (
            "non_manifold_edges",
            row.counted.non_manifold_edges.to_string(),
        ),
        ("nodal_expression", sym.expression.clone()),
        ("nodal_terms", sym.terms.to_string()),
        (
            "pinch_accounts_for_gap",
            row.pinch_accounts_for_gap().to_string(),
        ),
        (
            "primitive_cells_per_cubic_cell",
            row.kind.primitive_cells_per_cubic_cell().to_string(),
        ),
        ("primitive_cells_symbolic", sym.primitive_cells.to_string()),
        ("seam_pairs_identified", row.seam_pairs.to_string()),
        ("shift_classes", sym.classes_column()),
        ("shifts_invariant", sym.invariant.to_string()),
        ("shifts_negated", sym.negated.to_string()),
        ("shifts_neither", sym.neither.to_string()),
        (
            "transcription_residual",
            format!("{:.3e}", sym.transcription_residual),
        ),
        ("validate_chi", row.validate_chi.to_string()),
        (
            "validate_non_manifold_edges",
            row.validate_non_manifold_edges.to_string(),
        ),
        ("vertices", row.counted.vertices.to_string()),
        ("voxels_per_period", row.config.voxels.to_string()),
        ("weld_tolerance", format!("{:.3e}", row.weld_tolerance)),
        ("wrap_mode", row.wrap_mode.to_string()),
    ]);
}
