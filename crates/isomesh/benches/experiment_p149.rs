//! **P-149 — the features where a flat direction exists, isolated and measured
//! on their own.**
//!
//! Ticket: R-149. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p149
//! ```
//!
//! Writes `docs/experiments/p-149.csv`.
//!
//! # What was missing
//!
//! `P-146` measured the **aggregate** — one triangle ratio per field, read off a
//! matched-error fit — and its C1 was **FALSIFIED**: `docs/experiments/p-146.csv`
//! carries `c1_holds=false` on all forty rows, with `c1_population=4` (only
//! `FieldBound::Exact` admits a Hausdorff measurement, `fields/mod.rs:83-84`) and
//! `c1_winners=0`. Nothing in that file says **why**, because an aggregate cannot.
//! This row measures the **mechanism**: it bins the surface cells of each field by
//! their own principal-curvature ratio and asks whether the triangle saving lives
//! where the anisotropy is.
//!
//! The four numbers from `p-146.csv` this row is built against, quoted rather
//! than remembered:
//!
//! | field | `axis_ratio` (17³ → 65³) | `ratio` at matched Hausdorff | what the anisotropic arm did |
//! |---|---|---|---|
//! | `sphere`, `box_exact`, `csg_difference`, `gyroid`, `noise_cavity` | `1.000000` at every rung | `1.000000` where measurable | **nothing** — `triangles_anisotropic == triangles_isotropic` on every rung |
//! | `torus` | `1.117647 → 1.095238` | `1.153211` | spent **15.3% more** triangles for the same error |
//! | `thin_plate` | `7.666667 → 30.238095` | `2.298935` | 65³: `4088 → 536` triangles, but Hausdorff `0.045927933 → 0.141556797`, **3.08× worse** |
//! | `fbm_terrain` | `6.200000 → 47.400000` | unmeasurable (`bound=Unbounded`) | 65³: `16884 → 115174` triangles, **6.8× more** |
//!
//! So six of the eight reference fields make the two arms **the same grid**, and
//! that fact governs everything below: on those six every bin's saving is exactly
//! zero and could not have been anything else, which is `M-44`'s definition of a
//! non-measurement. They are recorded, with `arms_identical=true`, and their C1
//! and C2 verdicts read `unmeasurable:arms_identical` rather than `false` — a
//! falsification fitted through a constant is not a falsification.
//!
//! # The construction, restated from P-146 on purpose
//!
//! The task of this row is to **decompose** `P-146`'s saving, so it must measure
//! the same two meshes `P-146` measured. Every part of the anisotropic arm below
//! is restated from `crates/isomesh/benches/experiment_p146.rs` — [`band_points`]
//! (:413), [`axis_weights`] (its `census_of` weights, :548-550), [`round_odd`]
//! (:587), [`anisotropic_grid`] (:620) and [`Stretched`] (:676) — deliberately
//! and verbatim in behaviour, so the two rows are commensurable. **It is restated
//! and not shared**: `benches/common/` is owned by the module authors and R-146
//! kept these five inside its own file, so a consumer either restates them or
//! edits a file it does not own. This paragraph is the statement that they were
//! restated.
//!
//! The arm, in one sentence: the metric prescribes a point density
//! `√(e_aᵀ M e_a)` along each world axis, averaged over the surface band; the
//! anisotropic arm spends the **same total sample budget** `N³` on per-axis counts
//! `(n_x, n_y, n_z)` proportional to those weights, rounded to odd (`M-266`),
//! reached through the shipped extractor by the [`Stretched`] coordinate warp
//! because every `extract` in the crate takes a scalar `cell_size`
//! (`marching_cubes/mod.rs:193`). It is **per-axis global anisotropy, not per-cell
//! anisotropy** — `P-146`'s header says so at length and it is why six fields
//! come out isotropic.
//!
//! # The mechanism, in closed form, and where it is in the data
//!
//! The whole of Group D rests on the AM–GM gap: the anisotropic error constant
//! beats the isotropic one by the factor by which the geometric mean of the
//! curvature magnitudes falls below their arithmetic mean, and **that gap
//! collapses toward zero exactly where one principal curvature vanishes**. In two
//! surface dimensions the gap is a function of the ratio and of nothing else. With
//! `r = |κ₁| / |κ₀| ≥ 1`,
//!
//! ```text
//!     GM / AM  =  √(|κ₀ κ₁|) / ((|κ₀| + |κ₁|)/2)  =  2 √r / (1 + r),
//! ```
//!
//! which takes the value `1.000000` at an umbilic (`r = 1`), `0.993808` at the
//! `near_umbilic` edge `r = 1.25`, `0.942809` at `r = 2`, `0.745356` at `r = 5`,
//! `0.277297` at the `flat_direction` edge `r = 50` and `0.002000` at the
//! reported cap `r = 10⁶`. It is a **strictly decreasing** function of
//! `r` on `[1, ∞)`, so it is the same statement as the registration's monotonicity
//! read backwards. Both halves are recorded: `am_gm_gap_measured` is the median of
//! the per-cell `√(|κ₀κ₁|) / ((|κ₀|+|κ₁|)/2)` inside the bin and
//! `am_gm_gap_predicted` is `2√R/(1+R)` at the bin's median ratio `R`. They agree
//! to the width of the bin or the mechanism is not the one claimed.
//!
//! **`common::metric::am_gm_gap` is deliberately *not* used per bin**, and the
//! module's own doc says why (`benches/common/metric.rs:596-606`): it is
//! `‖√|det H|‖ / ‖tr|H|/d‖` over a population and carries **one power of a
//! curvature**, so comparing it across bins that differ in curvature *magnitude*
//! as well as in ratio would confound the two. The dimensionless GM/AM above
//! isolates the ratio exactly. That module function is R-147's instrument, used
//! across fields at one scale, which is the comparison it is homogeneous for.
//!
//! # The orientation of `principal_curvature_ratio`, decided in the open
//!
//! `common::metric::principal_curvatures` returns the two curvatures **ascending
//! by magnitude**, so `|κ₀| ≤ |κ₁|`, and its doc comment at `metric.rs:666` says
//! *"P-149's `principal_curvature_ratio` is `|κ₀| / |κ₁|`"* — a number in `(0, 1]`.
//! **This harness records the reciprocal, `|κ₁| / |κ₀| ≥ 1`, and the registration
//! is what forces that.** Under the module's orientation a spherical cap is still
//! `1`, but the saving would *decrease* with the ratio, the Spearman correlation
//! would sit near `−1`, and C1's own test — *"correlation above 0.7"* — would
//! record a **falsification on a mechanism that held perfectly**. The registration
//! pairs "monotone" with "above 0.7" and names "a non-monotone relationship" as
//! the falsifier, so the orientation it intends is the one in which a confirmed
//! mechanism reads positive. Neither convention is hidden: `curvature_ratio_
//! reciprocal` carries the module's orientation on every row and
//! `saving_vs_reciprocal_correlation` carries its correlation, which must be the
//! sign-flipped twin of the registered one.
//!
//! # The bins
//!
//! Six classes plus a residue, and every cell of every grid lands in exactly one.
//! A cell is classified when it **straddles** the isotropic grid's zero level and
//! `principal_curvatures` at its centre returns `Some`.
//!
//! | `feature_class` | window on `r` | what it is |
//! |---|---|---|
//! | `planar` | `‖κ‖∞ < 1e-6` | flat in **both** directions — a box face. Umbilic at zero, and not a spherical cap; excluded from the correlation because it has no curvature to be a ratio of |
//! | `umbilic` | `1 ≤ r < 1.25` | **C2's control**: spherical caps, and minimal-surface saddles — see below |
//! | `near_umbilic` | `1.25 ≤ r < 2` | |
//! | `mild_anisotropy` | `2 ≤ r < 5` | |
//! | `strong_anisotropy` | `5 ≤ r < 50` | |
//! | `flat_direction` | `r ≥ 50` | one curvature at or below the floor: the cylinder, the ridge, the box edge |
//! | `unclassified` | — | straddles but has no tangent plane (`‖∇f‖ ≤ GRAD_FLOOR`), **or** received a triangle from an arm in a cell the isotropic grid does not see the surface in. Recorded so no triangle is silently dropped; never in the correlation |
//!
//! `|κ₀|` is floored at `K_FLOOR = 1e-6` before the division and the ratio is
//! capped at `1e6`. The module states the reference fields' genuine curvatures run
//! `1e-2` to `1e2` (`metric.rs:62-64`), and a numerically-flat direction differs
//! at second order in `h` from `f64` round-off — of order `1e-14` — so `1e-6` sits
//! four orders below the smallest real curvature and eight above the noise.
//!
//! **The ratio is blind to sign, and one field exploits that.** A minimal surface
//! has `κ₁ = −κ₀` exactly, so a gyroid saddle has `r = 1` and lands in `umbilic`
//! beside a sphere. That is correct for the *mechanism* — a perfect saddle has no
//! flat direction and nothing to exploit — and wrong for the *word*. Every row
//! carries `elliptic_fraction`, the fraction of the bin's cells with `κ₀κ₁ > 0`,
//! which is `1` for a cap and `0` for a saddle, so a reader can tell which
//! `umbilic` bin they are looking at.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | **isotropic** | `(N, N, N)` samples, `h` equal on every axis | **yes** — uniform refinement, the baseline every `saving` is against |
//! | **anisotropic** | `(n_x, n_y, n_z)` from the metric, `∏ n_a ≈ N³` | no — this is P-146's arm, restated |
//! | **the `umbilic` bin** | cells with `r` near 1 inside the *same* two meshes | **yes** — C2's control, and the one that says the effect is anisotropy and not refinement |
//! | **the `planar` bin** | cells with no curvature at all | **yes** — a second, independent no-anisotropy control, reported and not registered |
//!
//! Both arms are extracted by one `MarchingCubes::new()` at its shipped defaults,
//! and both are binned on **one** classification — the isotropic grid's cells, by
//! the world-space centroid of each triangle. A comparison in which each arm is
//! binned by its own grid is not a comparison.
//!
//! # `hausdorff_delta`, which is local and says so
//!
//! `validate::accuracy` returns one scalar for a whole mesh and cannot be
//! restricted to a bin, so this row measures the two-sided distance itself, over
//! the bin's own cells, and only where it is a distance at all:
//!
//! - **mesh → surface.** For a `FieldBound::Exact` field `|f(p)|` **is** the
//!   distance to the zero set (`fields/mod.rs:83-84`), so the forward term is
//!   `max |f(v)|` over the vertices of the triangles binned here.
//! - **surface → mesh.** The bin's cell centres, Newton-projected onto the zero
//!   set along `∇f` in two steps, each measured to the **nearest vertex of the
//!   whole arm mesh** — not of the bin. That is what makes a hole visible: an
//!   anisotropic arm that emits nothing here is charged the distance to wherever
//!   it did emit, instead of quietly scoring a 100% saving. `probes` and
//!   `probe_residual_max` report the projection so the reader can see it worked.
//!   Nearest-*vertex* over-estimates nearest-point-on-mesh; it is applied
//!   identically to both arms, so the **delta** is fair and the level is an upper
//!   bound. Probes are strided to at most 128 per bin.
//!
//! `hausdorff_delta = max(forward, backward)_anisotropic − max(...)_isotropic`.
//! Positive means the saving was bought with error. The four non-`Exact` fields
//! read `unmeasurable:bound=…`, exactly as `p-146.csv` does.
//!
//! # Resolutions
//!
//! Two rungs, `33³` and `65³` — the authoring contract's default pair. The
//! registration names no ladder, and a bin population is a census rather than a
//! convergence study; two rungs are enough to show that a bin's monotonicity is
//! not an artefact of one grid phase, and cheap enough to keep the row inside two
//! minutes with the nearest-vertex search in it.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"C1 decomposes `P-146`'s C1 saving by feature
//! class; the bin populations are reported so each saving has a denominator."*
//! Discharged twice over. `P-146`'s C1 saving is, in the committed file,
//! **zero on all four measurable fields** (`c1_winners=0`), so what is decomposed
//! is a zero, and the decomposition's job is to say whether it is zero
//! *everywhere* or a cancellation. Every row carries `cells` — the bin's own
//! population — beside its `saving`, and `cells_all_fields` — that class's
//! population across the entire sweep, which is the vacuity control's own count.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record`, and every panic starts `"VOID: "`.
//! `M-44`: a zero that could not have been non-zero is not a measurement.
//!
//! - **The registered control.** Every one of the five ordered curvature-ratio
//!   bins must be non-empty on at least one field, or the monotonicity is fitted
//!   through gaps. Columns: `cells`, and `cells_all_fields` on every row.
//! - **The anisotropic arm must be anisotropic somewhere.** `axis_ratio > 1.5` on
//!   at least one group, or every `saving` in the file is identically zero by
//!   construction and the correlation is a correlation through a constant. Column:
//!   `axis_ratio`. Expected witnesses from `p-146.csv`: `thin_plate` and
//!   `fbm_terrain`.
//! - **The two arms must be on one budget.** `budget_ratio = ∏n_a / N³ ∈ [0.25, 4]`
//!   on every group, or a "saving" is a smaller sample count wearing a metric's
//!   clothes. Column: `budget_ratio`.
//! - **The band must be non-empty** on every group, or the per-axis weights the
//!   anisotropic grid is derived from are statistics of nothing.
//! - **At least one group must reach three populated ordered bins**, or C1 has no
//!   instrument anywhere rather than a narrow one. Column: `correlation_bins`.
//!
//! # Predicted verdicts, written before the harness ever ran
//!
//! C1 is predicted **FALSIFIED where it is measurable at all**. Six of eight
//! fields have `axis_ratio = 1.000000` in `p-146.csv`, so twelve of the sixteen
//! groups are `unmeasurable:arms_identical`; the four that remain are
//! `thin_plate` and `fbm_terrain` at two rungs each, and on both of those the
//! anisotropic arm's global behaviour (`thin_plate` `4088 → 536` with 3.08× the
//! error, `fbm_terrain` `16884 → 115174`) is dominated by a single pinned axis
//! rather than by a per-cell curvature ratio — a per-axis grid cannot spend a
//! flat direction that rotates, and cannot decline to spend one where the cell is
//! umbilic. The predicted finding is that **the mechanism is not measurable
//! through P-146's construction**, and the closed-form `2√r/(1+r)` agreement
//! between `am_gm_gap_measured` and `am_gm_gap_predicted` is the part of the
//! mechanism that *is* measured here, on all eight fields and both rungs.
//!
//! C2 is predicted **HELD wherever it is measurable**: `thin_plate`'s and
//! `fbm_terrain`'s `umbilic` bins should show `|saving| ≤ 0.05` if the effect is
//! anisotropy, and on the six identical-arm fields it reads
//! `unmeasurable:arms_identical` because a zero there was not free to be non-zero.
//!
//! # Determinism
//!
//! One thread, no PRNG, `f64` throughout. Sorting is [`f64::total_cmp`]. Grids are
//! swept `z`, `y`, `x` with `x` innermost, the crate's order. Probe striding is by
//! integer index. The only seeded object anywhere near this row is `FbmTerrain`'s
//! committed `0x5EED_1234`, which is the field's and not the harness's. No column
//! here is a wall clock: `R-151` owns the cost question and this row is a census.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

mod common;

use isomesh::fields::{FieldBound, ReferenceField};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

use crate::common::beta::rank_correlation;
use crate::common::metric::{GRAD_FLOOR, Sym3, hessian, metric_lp, principal_curvatures};

/// Samples per axis for the isotropic arm. Odd throughout (`M-266`).
const RUNGS: [u32; 2] = [33, 65];

/// The norm `M_Lp` is optimised for. One value, and `L²`, matching `P-146`;
/// `R-150` owns the norm sweep.
const P_NORM: f64 = 2.0;

/// A grid sample joins the surface band when `|f| <= BAND_CELLS · h`. One cell,
/// as in `P-146`.
const BAND_CELLS: f64 = 1.0;

/// Fewest samples any axis of the anisotropic grid may carry. Odd, and five
/// rather than two, exactly as `P-146` sets it.
const MIN_SAMPLES: u32 = 5;

/// `∏ n_a / N³` must land inside this, or the arms are not budget-matched.
const BUDGET_BAND: [f64; 2] = [0.25, 4.0];

/// The anisotropic arm must be anisotropic somewhere, by at least this factor.
const AXIS_RATIO_FLOOR: f64 = 1.5;

/// The absolute floor on `|κ|` below which a direction counts as flat.
///
/// The reference fields' genuine curvatures run `1e-2` to `1e2`
/// (`benches/common/metric.rs:62-64`); a numerically flat direction reads at
/// `f64` round-off divided by `h²`, of order `1e-14`. This sits four orders below
/// the smallest real curvature and eight above the noise.
const K_FLOOR: f64 = 1e-6;

/// The reported ratio is capped here, so `flat_direction`'s statistics are a
/// number rather than an infinity. `|κ₁|` at `1` over `K_FLOOR` is exactly this.
const RATIO_CAP: f64 = 1e6;

/// C1's bar on the Spearman correlation.
const CORRELATION_BAR: f64 = 0.7;

/// Fewest populated ordered bins a group needs before its correlation is a
/// correlation rather than two points and a line through them.
const C1_MIN_BINS: usize = 3;

/// C2's bar: the `umbilic` bin's saving must be inside this of zero.
const C2_TOLERANCE: f64 = 0.05;

/// Most Newton probes taken per bin for the surface → mesh direction.
const PROBE_CAP: usize = 128;

/// Newton steps used to put a cell centre on the zero set.
const NEWTON_STEPS: usize = 2;

/// `for_each_reference_field!` yields eight (`fields/mod.rs:211-255`).
const FIELDS: usize = 8;

/// Bin index of the flat-in-every-direction class.
const PLANAR: usize = 0;

/// Bin index of C2's control, the cells with `r` near one.
const UMBILIC: usize = 1;

/// First and last ordered bin, inclusive. `UMBILIC ..= FLAT_DIRECTION` are the
/// five bins the monotonicity is read over.
const FLAT_DIRECTION: usize = 5;

/// Bin index of the residue: no tangent plane, or a triangle in a cell the
/// isotropic grid does not see the surface in.
const UNCLASSIFIED: usize = 6;

/// How many bins a group carries.
const BINS: usize = 7;

/// The class names, in row-emission order. CSV-safe tokens.
const CLASS_NAMES: [&str; BINS] = [
    "planar",
    "umbilic",
    "near_umbilic",
    "mild_anisotropy",
    "strong_anisotropy",
    "flat_direction",
    "unclassified",
];

/// Inclusive lower edge of each ordered bin, in `UMBILIC ..= FLAT_DIRECTION`
/// order. Ascending, which is what makes [`class_of_ratio`]'s last-match search
/// correct.
const ORDERED_LO: [f64; 5] = [1.0, 1.25, 2.0, 5.0, 50.0];

/// Exclusive upper edge of each ordered bin, parallel to [`ORDERED_LO`].
const ORDERED_HI: [f64; 5] = [1.25, 2.0, 5.0, 50.0, f64::INFINITY];

// ─── curvature classification ────────────────────────────────────────────────

/// What one classified cell contributes.
#[derive(Clone, Copy)]
struct CellKind {
    /// Index into [`CLASS_NAMES`].
    class: usize,
    /// `|κ₁| / max(|κ₀|, K_FLOOR)`, capped at [`RATIO_CAP`]; `1` for a plane.
    ratio: f64,
    /// `|κ₀| / |κ₁|`, the orientation `metric.rs:666` names; `1` for a plane.
    reciprocal: f64,
    /// `√(|κ₀κ₁|) / ((|κ₀| + |κ₁|)/2)`, the dimensionless AM–GM gap; `1` for a
    /// plane, since both curvatures are equal there.
    am_gm: f64,
    /// `κ₀ κ₁ > 0` — a cap rather than a saddle.
    elliptic: bool,
}

/// The ordered bin a ratio falls in.
///
/// The last edge at or below `ratio` wins, which is well defined because
/// [`ORDERED_LO`] ascends and `ratio ≥ 1 = ORDERED_LO[0]` by construction.
fn class_of_ratio(ratio: f64) -> usize {
    let mut class = UMBILIC;
    for (offset, &edge) in ORDERED_LO.iter().enumerate() {
        if ratio >= edge {
            class = UMBILIC + offset;
        }
    }
    class
}

/// Classify one cell from its two principal curvatures.
fn classify(kappa: [f64; 2]) -> CellKind {
    let small = kappa[0].abs();
    let large = kappa[1].abs();
    let elliptic = kappa[0] * kappa[1] > 0.0;

    if large < K_FLOOR {
        // Flat in both directions. Umbilic at zero and not a spherical cap: the
        // ratio is `0/0`, defined here as `1` because the two curvatures are
        // equal, and the class keeps it out of the correlation.
        return CellKind {
            class: PLANAR,
            ratio: 1.0,
            reciprocal: 1.0,
            am_gm: 1.0,
            elliptic,
        };
    }

    let ratio = (large / small.max(K_FLOOR)).min(RATIO_CAP);
    CellKind {
        class: class_of_ratio(ratio),
        ratio,
        reciprocal: small / large,
        am_gm: (small * large).sqrt() / (0.5 * (small + large)),
        elliptic,
    }
}

/// The closed-form AM–GM gap at a ratio: `2√r / (1 + r)`.
///
/// Strictly decreasing on `[1, ∞)`, `1` at an umbilic and `→ 0` as one curvature
/// vanishes. This is the mechanism the registration's first sentence names, and
/// `am_gm_gap_measured` is compared against it on every row.
fn am_gm_of_ratio(ratio: f64) -> f64 {
    2.0 * ratio.sqrt() / (1.0 + ratio)
}

// ─── the anisotropic arm, restated from P-146 ────────────────────────────────

/// Nearest odd integer, at least one. Ties go up, deterministically.
///
/// Restated from `experiment_p146.rs:587`.
fn round_odd(x: f64) -> u32 {
    let half = ((x - 1.0) * 0.5).round();
    (2.0f64.mul_add(half, 1.0)).max(1.0) as u32
}

/// Per-axis sample counts from the metric's per-axis point densities, at the
/// isotropic arm's total budget.
///
/// Restated from `experiment_p146.rs:620`, behaviour for behaviour: `n_a ∝
/// weights[a]` with `∏ n_a = N³`; one **lower** clamp at [`MIN_SAMPLES`], after
/// which the remaining budget is re-solved over the axes still free, at most
/// three rounds because every round that changes anything pins one more axis.
/// There is no upper clamp, and R-146's header records the measured reason: a
/// ceiling can bind on two axes in the same round and double-count the budget.
///
/// Returns the counts and how many axes were pinned at the floor.
fn anisotropic_grid(weights: [f64; 3], samples: u32) -> ([u32; 3], usize) {
    let budget = f64::from(samples).powi(3);
    let mut pinned = [false; 3];
    let mut n = [samples; 3];

    for _round in 0..3 {
        let free: Vec<usize> = (0..3).filter(|&axis| !pinned[axis]).collect();
        if free.is_empty() {
            break;
        }
        let mut held = 1.0f64;
        for (axis, &fixed) in pinned.iter().enumerate() {
            if fixed {
                held *= f64::from(n[axis]);
            }
        }
        let count = free.len() as f64;
        let target = budget / held;
        let mut logsum = 0.0f64;
        for &axis in &free {
            logsum += weights[axis].ln();
        }
        let geometric_mean = (logsum / count).exp();
        let scale = target.powf(1.0 / count) / geometric_mean;

        let mut newly_pinned = false;
        for &axis in &free {
            let raw = round_odd(scale * weights[axis]);
            if raw < MIN_SAMPLES {
                n[axis] = MIN_SAMPLES;
                pinned[axis] = true;
                newly_pinned = true;
            } else {
                n[axis] = raw;
            }
        }
        if !newly_pinned {
            break;
        }
    }

    (n, pinned.iter().filter(|fixed| **fixed).count())
}

/// The field seen through a per-axis coordinate stretch, `sample(q) = f(lo + q ⊙ s)`.
///
/// Restated from `experiment_p146.rs:676`. Extracting this on a **cubic** grid of
/// `cell_size = h` in `q` is exactly extracting `f` on a rectilinear grid whose
/// physical spacings are `h · s`, which is the only way to reach an anisotropic
/// grid through an `extract` that takes a scalar `cell_size`
/// (`marching_cubes/mod.rs:193`). Positions are mapped back to world space by the
/// caller before anything measures them.
struct Stretched<'a, F> {
    field: &'a F,
    lo: [f64; 3],
    s: [f64; 3],
}

impl<F> Sdf for Stretched<'_, F>
where
    F: Sdf<Scalar = f64>,
{
    type Scalar = f64;

    fn sample(&self, q: [f64; 3]) -> f64 {
        self.field.sample([
            q[0].mul_add(self.s[0], self.lo[0]),
            q[1].mul_add(self.s[1], self.lo[1]),
            q[2].mul_add(self.s[2], self.lo[2]),
        ])
    }
}

/// The metric's mean point density along each world axis, over the surface band.
///
/// Restated from `experiment_p146.rs:548-560`: `mean √(e_aᵀ M e_a)`, and it is the
/// only thing the anisotropic split is derived from.
fn axis_weights<F>(field: &F, points: &[[f64; 3]], h: f64) -> [f64; 3]
where
    F: Sdf<Scalar = f64>,
{
    let mut weights = [0.0f64; 3];
    for &p in points {
        let metric: Sym3 = metric_lp(&hessian(field, p, h), P_NORM);
        for (axis, weight) in weights.iter_mut().enumerate() {
            *weight += metric.get(axis, axis).sqrt();
        }
    }
    let n = points.len() as f64;
    for weight in &mut weights {
        *weight /= n;
    }
    weights
}

// ─── the sampled grid ────────────────────────────────────────────────────────

/// Central-difference gradient at the same step the Hessian uses.
///
/// Not `Sdf::gradient`, whose default step is `Real::DIFF_STEP · max(|pᵢ|, 1)` and
/// therefore a different discrete object from the one the curvature was measured
/// with (`metric.rs:668-673`).
fn central_gradient<F>(field: &F, p: [f64; 3], h: f64) -> [f64; 3]
where
    F: Sdf<Scalar = f64>,
{
    let mut gradient = [0.0f64; 3];
    for (axis, slot) in gradient.iter_mut().enumerate() {
        let mut plus = p;
        plus[axis] += h;
        let mut minus = p;
        minus[axis] -= h;
        *slot = (field.sample(plus) - field.sample(minus)) / (2.0 * h);
    }
    gradient
}

/// The isotropic grid's values, swept `z`, `y`, `x` with `x` innermost.
fn sample_grid<F>(field: &F, origin: [f64; 3], h: f64, samples: u32) -> Vec<f64>
where
    F: Sdf<Scalar = f64>,
{
    let n = samples as usize;
    let mut values = Vec::with_capacity(n * n * n);
    for k in 0..samples {
        for j in 0..samples {
            for i in 0..samples {
                values.push(field.sample([
                    origin[0] + f64::from(i) * h,
                    origin[1] + f64::from(j) * h,
                    origin[2] + f64::from(k) * h,
                ]));
            }
        }
    }
    values
}

/// The grid samples within [`BAND_CELLS`] cells of the surface, read off the
/// value grid rather than re-sampled. Restated from `experiment_p146.rs:413`.
fn band_points(values: &[f64], origin: [f64; 3], h: f64, samples: u32) -> Vec<[f64; 3]> {
    let band = BAND_CELLS * h;
    let mut out = Vec::new();
    let mut index = 0usize;
    for k in 0..samples {
        for j in 0..samples {
            for i in 0..samples {
                if values[index].abs() <= band {
                    out.push([
                        origin[0] + f64::from(i) * h,
                        origin[1] + f64::from(j) * h,
                        origin[2] + f64::from(k) * h,
                    ]);
                }
                index += 1;
            }
        }
    }
    out
}

/// The cell of the isotropic grid containing a world point, clamped to the grid.
fn cell_of(p: [f64; 3], lo: [f64; 3], h: f64, cells: u32) -> usize {
    let last = f64::from(cells - 1);
    let mut index = [0u32; 3];
    for (axis, slot) in index.iter_mut().enumerate() {
        *slot = ((p[axis] - lo[axis]) / h).floor().clamp(0.0, last) as u32;
    }
    (index[0] + index[1] * cells + index[2] * cells * cells) as usize
}

/// Classify every cell of the isotropic grid that straddles the zero level.
///
/// Returns one slot per cell in `x`-fastest order: `Some(kind)` for a straddling
/// cell with a tangent plane, `None` otherwise. A cell whose eight corners agree
/// in sign is not a surface cell and is not asked.
fn classify_cells<F>(
    field: &F,
    values: &[f64],
    origin: [f64; 3],
    h: f64,
    samples: u32,
) -> (Vec<Option<CellKind>>, Vec<bool>)
where
    F: Sdf<Scalar = f64>,
{
    let cells = samples - 1;
    let n = samples as usize;
    let count = (cells as usize).pow(3);
    let mut kinds: Vec<Option<CellKind>> = Vec::with_capacity(count);
    let mut residue = vec![false; count];

    for cz in 0..cells {
        for cy in 0..cells {
            for cx in 0..cells {
                let base = cx as usize + cy as usize * n + cz as usize * n * n;
                let mut negatives = 0u32;
                for dz in 0..2usize {
                    for dy in 0..2usize {
                        for dx in 0..2usize {
                            if values[base + dx + dy * n + dz * n * n] < 0.0 {
                                negatives += 1;
                            }
                        }
                    }
                }
                if negatives == 0 || negatives == 8 {
                    kinds.push(None);
                    continue;
                }
                let centre = [
                    origin[0] + (f64::from(cx) + 0.5) * h,
                    origin[1] + (f64::from(cy) + 0.5) * h,
                    origin[2] + (f64::from(cz) + 0.5) * h,
                ];
                match principal_curvatures(field, centre, h) {
                    Some(kappa) => kinds.push(Some(classify(kappa))),
                    None => {
                        // Straddles, but the level set has no tangent plane
                        // there. Counted, never silently dropped.
                        residue[kinds.len()] = true;
                        kinds.push(None);
                    }
                }
            }
        }
    }

    (kinds, residue)
}

// ─── the bins ────────────────────────────────────────────────────────────────

/// One class's accumulation inside one `(field, resolution)` group.
struct Bin {
    cells: u64,
    ratios: Vec<f64>,
    reciprocals: Vec<f64>,
    am_gms: Vec<f64>,
    centres: Vec<[f64; 3]>,
    elliptic: u64,
    tri_iso: u64,
    tri_aniso: u64,
    /// `max |f(v)|` over the vertices of the triangles binned here, per arm.
    forward: [f64; 2],
    /// `max` over the bin's Newton-projected cell centres of the distance to the
    /// nearest vertex of the whole arm mesh.
    backward: [f64; 2],
    probes: usize,
    probe_residual: f64,
}

impl Bin {
    fn new() -> Self {
        Self {
            cells: 0,
            ratios: Vec::new(),
            reciprocals: Vec::new(),
            am_gms: Vec::new(),
            centres: Vec::new(),
            elliptic: 0,
            tri_iso: 0,
            tri_aniso: 0,
            forward: [0.0; 2],
            backward: [0.0; 2],
            probes: 0,
            probe_residual: 0.0,
        }
    }

    /// `1 − aniso / iso`, or `None` where the isotropic arm put no triangle here
    /// and the ratio has no denominator.
    fn saving(&self) -> Option<f64> {
        if self.tri_iso == 0 {
            return None;
        }
        Some(1.0 - self.tri_aniso as f64 / self.tri_iso as f64)
    }
}

/// Median of a sample, or `None` for an empty one. Sorted by
/// [`f64::total_cmp`], upper middle on an even count — the same rule
/// `experiment_p146.rs:400-407` uses, so the two files report one statistic.
fn median(xs: &[f64]) -> Option<f64> {
    if xs.is_empty() {
        return None;
    }
    let mut sorted = xs.to_vec();
    sorted.sort_by(f64::total_cmp);
    Some(sorted[sorted.len() / 2])
}

/// Arithmetic mean of a sample, summed in slice order, or `None` for an empty one.
fn mean(xs: &[f64]) -> Option<f64> {
    if xs.is_empty() {
        return None;
    }
    Some(xs.iter().sum::<f64>() / xs.len() as f64)
}

// ─── one group ───────────────────────────────────────────────────────────────

/// Everything one `(field, resolution)` produced.
struct Group {
    field: &'static str,
    bound: &'static str,
    /// The token a distance column carries when this field's `bound()` does not
    /// make `|f|` a distance. Same vocabulary `p-146.csv` uses.
    reason: &'static str,
    exact: bool,
    samples: u32,
    grid: [u32; 3],
    pinned: usize,
    axis_ratio: f64,
    budget_ratio: f64,
    band: usize,
    /// The anisotropic grid is the isotropic grid, so every saving in this group
    /// is zero and could not have been anything else (`M-44`).
    arms_identical: bool,
    verts_iso: usize,
    verts_aniso: usize,
    bins: Vec<Bin>,
    correlation: Option<f64>,
    reciprocal_correlation: Option<f64>,
    correlation_bins: usize,
    c1: String,
    c2: String,
}

impl Group {
    /// This group's per-bin `(median ratio, saving)` pairs, over the ordered bins
    /// that have both. The pairs C1's correlation is taken over.
    fn correlation_pairs(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let mut ratios = Vec::new();
        let mut reciprocals = Vec::new();
        let mut savings = Vec::new();
        for index in UMBILIC..=FLAT_DIRECTION {
            let bin = &self.bins[index];
            let (Some(ratio), Some(reciprocal), Some(saving)) = (
                median(&bin.ratios),
                median(&bin.reciprocals),
                bin.saving(),
            ) else {
                continue;
            };
            ratios.push(ratio);
            reciprocals.push(reciprocal);
            savings.push(saving);
        }
        (ratios, reciprocals, savings)
    }
}

/// Why a field's distances cannot be measured, from its declared bound.
///
/// Same tokens `p-146.csv` uses, so the two files can be joined.
fn unmeasurable_reason(bound: FieldBound) -> &'static str {
    match bound {
        FieldBound::Exact => "measurable",
        FieldBound::Lipschitz { .. } => "unmeasurable:bound=Lipschitz",
        FieldBound::Underestimate { .. } => "unmeasurable:bound=Underestimate",
        FieldBound::Unbounded => "unmeasurable:bound=Unbounded",
    }
}

/// The name of a `FieldBound` variant, without its parameters — the CSV writer
/// refuses a `,` inside a value and `Lipschitz { l: 3.46 }` has one.
fn bound_name(bound: FieldBound) -> &'static str {
    match bound {
        FieldBound::Exact => "Exact",
        FieldBound::Lipschitz { .. } => "Lipschitz",
        FieldBound::Underestimate { .. } => "Underestimate",
        FieldBound::Unbounded => "Unbounded",
    }
}

/// Two Newton steps along `∇f` from a cell centre onto the zero set.
fn project<F>(field: &F, start: [f64; 3], h: f64) -> [f64; 3]
where
    F: Sdf<Scalar = f64>,
{
    let mut q = start;
    for _step in 0..NEWTON_STEPS {
        let g = central_gradient(field, q, h);
        let norm2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
        if !norm2.is_finite() || norm2.sqrt() <= GRAD_FLOOR {
            // No direction to move along. The residual column reports it.
            break;
        }
        let value = field.sample(q);
        for (axis, slot) in q.iter_mut().enumerate() {
            *slot -= value * g[axis] / norm2;
        }
    }
    q
}

/// Distance from a point to the nearest vertex of a mesh, or `None` for an empty
/// mesh.
fn nearest_vertex(point: [f64; 3], positions: &[[f64; 3]]) -> Option<f64> {
    let mut best = f64::INFINITY;
    for v in positions {
        let dx = v[0] - point[0];
        let dy = v[1] - point[1];
        let dz = v[2] - point[2];
        let d2 = dx.mul_add(dx, dy.mul_add(dy, dz * dz));
        if d2 < best {
            best = d2;
        }
    }
    if best.is_finite() { Some(best.sqrt()) } else { None }
}

/// Bin one arm's triangles by the isotropic grid cell their centroid falls in,
/// and accumulate that arm's forward `max |f(v)|` per bin.
fn bin_triangles<F>(
    field: &F,
    mesh: &MeshBuffer<f64>,
    kinds: &[Option<CellKind>],
    residue: &mut [bool],
    bins: &mut [Bin],
    arm: usize,
    lo: [f64; 3],
    h: f64,
    cells: u32,
    exact: bool,
) where
    F: Sdf<Scalar = f64>,
{
    for triangle in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[triangle[0] as usize];
        let b = mesh.positions[triangle[1] as usize];
        let c = mesh.positions[triangle[2] as usize];
        let centroid = [
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        ];
        let cell = cell_of(centroid, lo, h, cells);
        let class = match kinds[cell] {
            Some(kind) => kind.class,
            None => {
                // Either a straddling cell without a tangent plane or, for the
                // anisotropic arm, a cell the isotropic grid does not see the
                // surface in at all. Both are the residue and both are counted.
                residue[cell] = true;
                UNCLASSIFIED
            }
        };
        if arm == 0 {
            bins[class].tri_iso += 1;
        } else {
            bins[class].tri_aniso += 1;
        }
        if exact {
            for vertex in [a, b, c] {
                let value = field.sample(vertex).abs();
                if value > bins[class].forward[arm] {
                    bins[class].forward[arm] = value;
                }
            }
        }
    }
}

/// Measure one reference field at one resolution.
fn measure<F>(field: &F, name: &'static str, samples: u32) -> Group
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let bound = field.bound();
    let exact = bound.is_exact();
    let (shape, origin, h) = common::grid::<f64, _>(field, samples);
    let (lo, hi) = field.domain();
    let extent = hi[0] - lo[0];
    let cells = samples - 1;

    // ── the classification, on the isotropic grid ────────────────────────────
    let values = sample_grid(field, origin, h, samples);
    let (kinds, mut residue) = classify_cells(field, &values, origin, h, samples);

    // ── the anisotropic grid, from P-146's per-axis weights ──────────────────
    let band = band_points(&values, origin, h, samples);
    assert!(
        !band.is_empty(),
        "VOID: {name} at {samples}^3 put no grid sample within {BAND_CELLS} cell of its surface, \
         so the per-axis weights the anisotropic arm is built from are statistics of an empty \
         population and every saving below is a saving over nothing (M-44)"
    );
    let weights = axis_weights(field, &band, h);
    let (grid, pinned) = anisotropic_grid(weights, samples);

    // ── both arms, one extractor at its shipped defaults ─────────────────────
    let mut mc = MarchingCubes::<f64>::new();
    let mut iso = MeshBuffer::<f64>::new();
    mc.extract(field, &shape, origin, h, &mut iso)
        .expect("isotropic extraction over the reference grid");

    let stretch = [
        extent / f64::from(grid[0] - 1) / h,
        extent / f64::from(grid[1] - 1) / h,
        extent / f64::from(grid[2] - 1) / h,
    ];
    let stretched = Stretched {
        field,
        lo,
        s: stretch,
    };
    let aniso_shape = RuntimeShape3::new(grid).expect("anisotropic grid fits u32");
    let mut aniso = MeshBuffer::<f64>::new();
    mc.extract(&stretched, &aniso_shape, [0.0; 3], h, &mut aniso)
        .expect("anisotropic extraction over the stretched grid");
    for p in &mut aniso.positions {
        p[0] = p[0].mul_add(stretch[0], lo[0]);
        p[1] = p[1].mul_add(stretch[1], lo[1]);
        p[2] = p[2].mul_add(stretch[2], lo[2]);
    }

    // ── the bins ─────────────────────────────────────────────────────────────
    let mut bins: Vec<Bin> = (0..BINS).map(|_| Bin::new()).collect();
    bin_triangles(
        field,
        &iso,
        &kinds,
        &mut residue,
        &mut bins,
        0,
        origin,
        h,
        cells,
        exact,
    );
    bin_triangles(
        field,
        &aniso,
        &kinds,
        &mut residue,
        &mut bins,
        1,
        origin,
        h,
        cells,
        exact,
    );

    for (cell, kind) in kinds.iter().enumerate() {
        let Some(kind) = kind else {
            if residue[cell] {
                bins[UNCLASSIFIED].cells += 1;
            }
            continue;
        };
        let bin = &mut bins[kind.class];
        bin.cells += 1;
        bin.ratios.push(kind.ratio);
        bin.reciprocals.push(kind.reciprocal);
        bin.am_gms.push(kind.am_gm);
        bin.elliptic += u64::from(kind.elliptic);
        let z = cell / (cells as usize * cells as usize);
        let y = (cell / cells as usize) % cells as usize;
        let x = cell % cells as usize;
        bin.centres.push([
            origin[0] + (x as f64 + 0.5) * h,
            origin[1] + (y as f64 + 0.5) * h,
            origin[2] + (z as f64 + 0.5) * h,
        ]);
    }

    // ── the surface → mesh direction, only where `|f|` is a distance ─────────
    if exact && !iso.positions.is_empty() && !aniso.positions.is_empty() {
        for bin in &mut bins {
            if bin.centres.is_empty() {
                continue;
            }
            let stride = bin.centres.len().div_ceil(PROBE_CAP);
            for centre in bin.centres.iter().step_by(stride) {
                let q = project(field, *centre, h);
                let residual = field.sample(q).abs();
                if residual > bin.probe_residual {
                    bin.probe_residual = residual;
                }
                bin.probes += 1;
                for (arm, positions) in [&iso.positions, &aniso.positions].into_iter().enumerate() {
                    let Some(distance) = nearest_vertex(q, positions) else {
                        continue;
                    };
                    if distance > bin.backward[arm] {
                        bin.backward[arm] = distance;
                    }
                }
            }
        }
    }

    // ── the verdicts this group decides ──────────────────────────────────────
    let budget_ratio = (f64::from(grid[0]) * f64::from(grid[1]) * f64::from(grid[2]))
        / f64::from(samples).powi(3);
    let axis_hi = f64::from(grid.iter().copied().max().unwrap_or(samples));
    let axis_lo = f64::from(grid.iter().copied().min().unwrap_or(samples));
    let arms_identical = grid == [samples; 3];

    let mut group = Group {
        field: name,
        bound: bound_name(bound),
        reason: unmeasurable_reason(bound),
        exact,
        samples,
        grid,
        pinned,
        axis_ratio: axis_hi / axis_lo,
        budget_ratio,
        band: band.len(),
        arms_identical,
        verts_iso: iso.positions.len(),
        verts_aniso: aniso.positions.len(),
        bins,
        correlation: None,
        reciprocal_correlation: None,
        correlation_bins: 0,
        c1: String::new(),
        c2: String::new(),
    };

    let (ratios, reciprocals, savings) = group.correlation_pairs();
    group.correlation_bins = ratios.len();
    if ratios.len() >= C1_MIN_BINS {
        group.correlation = Some(rank_correlation(&ratios, &savings));
        group.reciprocal_correlation = Some(rank_correlation(&reciprocals, &savings));
    }

    group.c1 = if arms_identical {
        String::from("unmeasurable:arms_identical")
    } else {
        match group.correlation {
            Some(rho) => (rho >= CORRELATION_BAR).to_string(),
            None => format!("unmeasurable:bins={}<{C1_MIN_BINS}", group.correlation_bins),
        }
    };
    group.c2 = if arms_identical {
        String::from("unmeasurable:arms_identical")
    } else {
        match group.bins[UMBILIC].saving() {
            Some(saving) => (saving.abs() <= C2_TOLERANCE).to_string(),
            None => String::from("unmeasurable:no_umbilic_bin"),
        }
    };

    group
}

// ─── reporting ───────────────────────────────────────────────────────────────

/// One group's block on the console.
fn report(group: &Group) {
    println!(
        "{:<15} {:>3}^3  grid {}x{}x{} (axis_ratio {:.3}, budget {:.3}, pinned {}) bound {} \
         verts {} -> {}",
        group.field,
        group.samples,
        group.grid[0],
        group.grid[1],
        group.grid[2],
        group.axis_ratio,
        group.budget_ratio,
        group.pinned,
        group.bound,
        group.verts_iso,
        group.verts_aniso
    );
    for (index, bin) in group.bins.iter().enumerate() {
        println!(
            "    {:<18} cells {:>7}  r_med {:>12}  gm/am {:>8} (pred {:>8})  tri {:>7} -> {:>7}  \
             saving {:>10}  elliptic {:>8}",
            CLASS_NAMES[index],
            bin.cells,
            median(&bin.ratios).map_or_else(|| String::from("-"), |r| format!("{r:.4e}")),
            median(&bin.am_gms).map_or_else(|| String::from("-"), |g| format!("{g:.5}")),
            median(&bin.ratios)
                .map_or_else(|| String::from("-"), |r| format!("{:.5}", am_gm_of_ratio(r))),
            bin.tri_iso,
            bin.tri_aniso,
            bin.saving()
                .map_or_else(|| String::from("-"), |s| format!("{s:.5}")),
            if bin.cells == 0 {
                String::from("-")
            } else {
                format!("{:.4}", bin.elliptic as f64 / bin.cells as f64)
            },
        );
    }
    println!(
        "    correlation over {} ordered bins: {} (reciprocal {})  C1 {}  C2 {}\n",
        group.correlation_bins,
        group
            .correlation
            .map_or_else(|| String::from("-"), |r| format!("{r:+.4}")),
        group
            .reciprocal_correlation
            .map_or_else(|| String::from("-"), |r| format!("{r:+.4}")),
        group.c1,
        group.c2
    );
}

/// A value, or the token that says why there is none.
fn or_undefined(value: Option<f64>, digits: usize) -> String {
    value.map_or_else(|| String::from("undefined"), |v| format!("{v:.digits$}"))
}

// ─── the run ─────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-149");

    common::experiment::run(prereg, |run| {
        println!(
            "construction: P-146's per-axis global anisotropic GRID, restated verbatim in \
             behaviour, binned per cell.\n  ratio = |k1|/max(|k0|, {K_FLOOR:e}) >= 1, capped at \
             {RATIO_CAP:e}; bins {ORDERED_LO:?} .. {ORDERED_HI:?}\n  mechanism: gm/am = \
             2*sqrt(r)/(1+r), collapsing to 0 as one curvature vanishes\n  ladder {RUNGS:?} \
             samples/axis, metric M_Lp p = {P_NORM}, band |f| <= {BAND_CELLS}h\n"
        );

        let mut groups: Vec<Group> = Vec::with_capacity(FIELDS * RUNGS.len());
        isomesh::for_each_reference_field!(f64, |name, field| {
            for samples in RUNGS {
                groups.push(measure(&field, name, samples));
            }
        });
        assert_eq!(
            groups.len(),
            FIELDS * RUNGS.len(),
            "P-149: for_each_reference_field! must yield {FIELDS} fields at {} rungs",
            RUNGS.len()
        );

        for group in &groups {
            report(group);
        }

        // ── vacuity controls, all before the first record ────────────────────
        //
        // M-44: a zero that could not have been non-zero is not a measurement.
        let mut class_cells = [0u64; BINS];
        for group in &groups {
            for (index, bin) in group.bins.iter().enumerate() {
                class_cells[index] += bin.cells;
            }
        }
        for index in UMBILIC..=FLAT_DIRECTION {
            assert!(
                class_cells[index] > 0,
                "VOID: the curvature-ratio bin `{}` ({} <= r < {}) is empty on every one of the \
                 {FIELDS} reference fields at every rung, so C1's monotonicity would be fitted \
                 through a gap. That is the registration's own vacuity control. Populations: {}",
                CLASS_NAMES[index],
                ORDERED_LO[index - UMBILIC],
                ORDERED_HI[index - UMBILIC],
                (0..BINS)
                    .map(|b| format!("{}={}", CLASS_NAMES[b], class_cells[b]))
                    .collect::<Vec<_>>()
                    .join(" ")
            );
        }

        let axis_ratio = groups
            .iter()
            .map(|group| group.axis_ratio)
            .fold(0.0f64, f64::max);
        assert!(
            axis_ratio > AXIS_RATIO_FLOOR,
            "VOID: the most anisotropic grid the metric asked for anywhere is {axis_ratio:.4}:1, \
             below {AXIS_RATIO_FLOOR}:1. Both arms are then the same grid in every group, every \
             `saving` in the file is identically zero by construction, and the correlation C1 \
             reads is a correlation through a constant"
        );

        for group in &groups {
            assert!(
                group.budget_ratio >= BUDGET_BAND[0] && group.budget_ratio <= BUDGET_BAND[1],
                "VOID: {} at {}^3 spent {:.4}x the isotropic arm's sample budget on {}x{}x{}, \
                 outside [{}, {}]. A saving over an arm that was handed fewer samples is not a \
                 saving from anisotropy",
                group.field,
                group.samples,
                group.budget_ratio,
                group.grid[0],
                group.grid[1],
                group.grid[2],
                BUDGET_BAND[0],
                BUDGET_BAND[1]
            );
        }

        let measurable_groups = groups
            .iter()
            .filter(|group| group.correlation_bins >= C1_MIN_BINS)
            .count();
        assert!(
            measurable_groups > 0,
            "VOID: no group anywhere in the sweep reached {C1_MIN_BINS} populated ordered bins, \
             so C1 has no instrument at all rather than a narrow one. Best was {} bins",
            groups
                .iter()
                .map(|group| group.correlation_bins)
                .max()
                .unwrap_or(0)
        );

        // ── the pooled correlation, per resolution ───────────────────────────
        //
        // Every (field, ordered bin) pair at one rung, which is the widest read
        // of the same question and the one that does not depend on any single
        // field having five populated bins.
        let mut pooled: Vec<(u32, Option<f64>, usize)> = Vec::with_capacity(RUNGS.len());
        for samples in RUNGS {
            let mut ratios = Vec::new();
            let mut savings = Vec::new();
            for group in groups.iter().filter(|group| group.samples == samples) {
                let (group_ratios, _, group_savings) = group.correlation_pairs();
                ratios.extend(group_ratios);
                savings.extend(group_savings);
            }
            let count = ratios.len();
            let rho = if count >= C1_MIN_BINS {
                Some(rank_correlation(&ratios, &savings))
            } else {
                None
            };
            pooled.push((samples, rho, count));
        }

        // ── the global verdicts ──────────────────────────────────────────────
        let c1_measurable = groups
            .iter()
            .filter(|group| group.c1 == "true" || group.c1 == "false")
            .count();
        let c1_holding = groups.iter().filter(|group| group.c1 == "true").count();
        let c1_global = c1_measurable > 0 && c1_holding == c1_measurable;

        let c2_measurable = groups
            .iter()
            .filter(|group| group.c2 == "true" || group.c2 == "false")
            .count();
        let c2_holding = groups.iter().filter(|group| group.c2 == "true").count();
        let c2_global = c2_measurable > 0 && c2_holding == c2_measurable;

        println!(
            "bin populations across the sweep: {}",
            (0..BINS)
                .map(|b| format!("{}={}", CLASS_NAMES[b], class_cells[b]))
                .collect::<Vec<_>>()
                .join(" ")
        );
        for (samples, rho, count) in &pooled {
            println!(
                "pooled correlation at {samples}^3 over {count} (field, bin) pairs: {}",
                rho.map_or_else(|| String::from("undefined"), |r| format!("{r:+.6}"))
            );
        }
        println!(
            "C1: {c1_holding} of {c1_measurable} measurable groups reach rho >= \
             {CORRELATION_BAR} -> {c1_global}   ({} of {} groups had identical arms)",
            groups.iter().filter(|group| group.arms_identical).count(),
            groups.len()
        );
        println!(
            "C2: {c2_holding} of {c2_measurable} measurable groups keep the umbilic bin inside \
             +/-{C2_TOLERANCE} -> {c2_global}\n"
        );

        // ── the rows ─────────────────────────────────────────────────────────
        for group in &groups {
            let reason = group.reason;
            let meshes_present = group.verts_iso > 0 && group.verts_aniso > 0;
            let pooled_rho = pooled
                .iter()
                .find(|(samples, _, _)| *samples == group.samples)
                .and_then(|(_, rho, _)| *rho);

            for (index, bin) in group.bins.iter().enumerate() {
                let ratio_median = median(&bin.ratios);
                let in_correlation = (UMBILIC..=FLAT_DIRECTION).contains(&index)
                    && ratio_median.is_some()
                    && bin.saving().is_some();

                let hausdorff = |arm: usize| bin.forward[arm].max(bin.backward[arm]);
                let hausdorff_delta = if !group.exact {
                    reason.to_string()
                } else if !meshes_present {
                    String::from("unmeasurable:empty_mesh")
                } else if bin.cells == 0 && bin.tri_iso == 0 && bin.tri_aniso == 0 {
                    String::from("undefined:empty_bin")
                } else {
                    format!("{:.9}", hausdorff(1) - hausdorff(0))
                };
                let hausdorff_arm = |arm: usize| {
                    if group.exact && meshes_present {
                        format!("{:.9}", hausdorff(arm))
                    } else if group.exact {
                        String::from("unmeasurable:empty_mesh")
                    } else {
                        reason.to_string()
                    }
                };

                run.record(&[
                    ("feature_class", CLASS_NAMES[index].to_string()),
                    (
                        "principal_curvature_ratio",
                        or_undefined(ratio_median, 6),
                    ),
                    ("cells", bin.cells.to_string()),
                    ("triangles_isotropic", bin.tri_iso.to_string()),
                    ("triangles_anisotropic", bin.tri_aniso.to_string()),
                    (
                        "saving",
                        bin.saving().map_or_else(
                            || String::from("undefined:no_isotropic_triangles"),
                            |s| format!("{s:.6}"),
                        ),
                    ),
                    ("hausdorff_delta", hausdorff_delta),
                    (
                        "saving_vs_curvature_ratio_correlation",
                        group.correlation.map_or_else(
                            || format!("unmeasurable:bins={}<{C1_MIN_BINS}", group.correlation_bins),
                            |r| format!("{r:.6}"),
                        ),
                    ),
                    ("c1_holds", group.c1.clone()),
                    ("c2_holds", group.c2.clone()),
                    // ── extras (M-273) ──
                    (
                        "am_gm_gap_measured",
                        or_undefined(median(&bin.am_gms), 6),
                    ),
                    (
                        "am_gm_gap_predicted",
                        or_undefined(ratio_median.map(am_gm_of_ratio), 6),
                    ),
                    ("arms_identical", group.arms_identical.to_string()),
                    ("axis_ratio", format!("{:.6}", group.axis_ratio)),
                    ("axes_pinned", group.pinned.to_string()),
                    ("band_points", group.band.to_string()),
                    ("budget_ratio", format!("{:.6}", group.budget_ratio)),
                    ("c1_global_holds", c1_global.to_string()),
                    ("c1_groups_holding", c1_holding.to_string()),
                    ("c1_groups_measurable", c1_measurable.to_string()),
                    ("c2_global_holds", c2_global.to_string()),
                    (
                        "c2_umbilic_saving",
                        or_undefined(group.bins[UMBILIC].saving(), 6),
                    ),
                    ("cells_all_fields", class_cells[index].to_string()),
                    ("correlation_bins", group.correlation_bins.to_string()),
                    (
                        "curvature_ratio_mean",
                        or_undefined(mean(&bin.ratios), 6),
                    ),
                    (
                        "curvature_ratio_reciprocal",
                        or_undefined(median(&bin.reciprocals), 9),
                    ),
                    (
                        "elliptic_fraction",
                        or_undefined(
                            (bin.cells > 0).then(|| bin.elliptic as f64 / bin.cells as f64),
                            6,
                        ),
                    ),
                    ("field", group.field.to_string()),
                    ("field_bound", group.bound.to_string()),
                    (
                        "grid_anisotropic",
                        format!("{}x{}x{}", group.grid[0], group.grid[1], group.grid[2]),
                    ),
                    ("grid_isotropic", format!("{0}x{0}x{0}", group.samples)),
                    ("hausdorff_anisotropic_local", hausdorff_arm(1)),
                    ("hausdorff_isotropic_local", hausdorff_arm(0)),
                    ("in_correlation", in_correlation.to_string()),
                    (
                        "pooled_correlation",
                        or_undefined(pooled_rho, 6),
                    ),
                    (
                        "probe_residual_max",
                        if group.exact {
                            format!("{:.6e}", bin.probe_residual)
                        } else {
                            reason.to_string()
                        },
                    ),
                    ("probes", bin.probes.to_string()),
                    (
                        "ratio_bin_hi",
                        if (UMBILIC..=FLAT_DIRECTION).contains(&index) {
                            format!("{:.6}", ORDERED_HI[index - UMBILIC])
                        } else {
                            String::from("undefined")
                        },
                    ),
                    (
                        "ratio_bin_lo",
                        if (UMBILIC..=FLAT_DIRECTION).contains(&index) {
                            format!("{:.6}", ORDERED_LO[index - UMBILIC])
                        } else {
                            String::from("undefined")
                        },
                    ),
                    ("resolution", group.samples.to_string()),
                    (
                        "saving_vs_reciprocal_correlation",
                        group.reciprocal_correlation.map_or_else(
                            || format!("unmeasurable:bins={}<{C1_MIN_BINS}", group.correlation_bins),
                            |r| format!("{r:.6}"),
                        ),
                    ),
                    ("vertices_anisotropic", group.verts_aniso.to_string()),
                    ("vertices_isotropic", group.verts_iso.to_string()),
                ]);
            }
        }
    });
}
