//! **P-150 — a null registered on purpose: the `L^1` anomaly in the source's
//! own table, and whether our own fields reproduce it.**
//!
//! Ticket: R-150. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p150
//! ```
//!
//! Writes `docs/experiments/p-150.csv`.
//!
//! # What was missing
//!
//! Group D's headline number is a comparison **across two norms**, and the
//! registration exists to say so before anyone quotes it. Cao's Table 2
//! (`10.1090/S0025-5718-07-01981-3`, in corpus) shows `5.43e-7` against
//! `2.79e-6` at `N_e ~ 16,000` — a real `5.14x` spread — but the two metrics in
//! that comparison are optimised for **different** norms, in `H^1` the ranking
//! reverses, and Cao himself writes that his theory holds *"in all the cases
//! except the `L^1`-error of quadratic interpolation"*. The `L^1` column is the
//! one place his own optimal metric does **not** win. Building on it is the
//! weakest available ground.
//!
//! So this row is a methodology guard, not a win: it asks whether **the metric
//! optimised for a given norm wins in that norm** on our eight reference fields,
//! reports every comparison **inside** one norm, and characterises anything that
//! violates it.
//!
//! It consumes `common::metric` unchanged — `hessian`, `metric_lp`,
//! `complexity`, `aspect_ratio` and `Sym3::eigen` — owned by R-146
//! (`benches/common/metric.rs`). The metric is
//!
//! ```text
//!     M_Lp = D_Lp · det(|H_u|)^(−1/(2p + 3)) · |H_u|,
//! ```
//!
//! and its primary is **NASA NTRS 20200003084**, which restates Loseille &
//! Alauzet verbatim; the two SIAM originals are paywalled and the corpus holds
//! only their landing pages. A finding must cite the restatement.
//!
//! ## What P-146 measured about this machinery, quoted from its CSV
//!
//! Read from `docs/experiments/p-146.csv` (commit `77e71cb`), not from a
//! summary:
//!
//! - **P-146 C1 and C3 were both FALSIFIED.** The metric-driven arm did not
//!   reach 25% fewer triangles at matched Hausdorff, and `metric_share` ran to
//!   `2.023210` on `gyroid` at `65³` — 202% of extraction against a 15% bar. C2
//!   reads `unmeasurable` on twenty of forty rows because `validate::accuracy`
//!   is meaningless where `field.bound()` is not `Exact`. **So this row does not
//!   inherit a win from R-146 and does not claim one.** What it inherits is the
//!   metric field itself, which P-146 established is computable and anisotropic.
//! - **`aspect_ratio_max` is mostly `H_FLOOR` talking, and P-146's own numbers
//!   prove it.** At `65³`: `box_exact` reports `aspect_ratio_max = 2.4e10` with
//!   `aspect_ratio_max_off_floor = 3.824847e0` and `at_floor_cells = 18042` of
//!   `band_points = 18082` (`at_floor_fraction = 0.997788`); `fbm_terrain`
//!   reports `aspect_ratio_max = 8.544427e9` with
//!   `aspect_ratio_max_off_floor = 0.000000e0` and `at_floor_fraction =
//!   1.000000`, because a heightfield SDF is exactly linear in `y`; `gyroid`
//!   reports `at_floor_cells = 7` of `21655` with a genuine
//!   `aspect_ratio_max_off_floor = 3.511176e5`.
//!
//!   `local://shared-modules.md` states the rule that follows and this file
//!   obeys it: **every row that carries an aspect ratio also carries
//!   `at_floor_cells` and `at_floor_fraction`**, or the maximum is a restatement
//!   of `H_FLOOR = 1e-9` rather than a measurement
//!   (`benches/common/metric.rs:67-74`).
//!
//! # The construction
//!
//! The framework is the continuous mesh: a metric field `M` is a mesh, its
//! **complexity** `C(M) = ∫ √det M` stands in for the element count, and the
//! local interpolation error `e(M, x)` of the unit element at `x` is integrated
//! to give the global error
//!
//! ```text
//!     E_p(M) = ( ∫_Ω e(M, x)^p dx )^{1/p}.
//! ```
//!
//! Four arms are built, one per exponent `q ∈ {1, 2, 3, ∞}`, each **scaled to
//! the same complexity** so the comparison is at a fixed element count. Scaling
//! is the one thing `common::metric` leaves to the caller: it folds `D_Lp` to
//! `1` because every quantity P-146 reports divides it out
//! (`benches/common/metric.rs:44-53`). Here it cannot be folded — matching two
//! metrics on budget is exactly what `D_Lp` is for — so it is set per arm by
//! `C(s·M) = s^{3/2} C(M)`, i.e.
//!
//! ```text
//!     s_q = (N / C(M_q))^{2/3},        N = the band's sample count,
//! ```
//!
//! and `complexity_rel_error` on every row proves the match landed. It lands to
//! `1e-8`, six orders tighter than the `1.4%` Cao's own table manages — which is
//! the caveat the registration flags in his numbers, closed in ours.
//!
//! ## `e_K`: the element, and the twenty-one points its error is taken over
//!
//! At band sample `x` with scaled metric `M`, eigendecomposed by
//! `Sym3::eigen`, the **unit-metric tetrahedron** has edge vectors
//! `a_i = v_i / √λ_i` (so `a_iᵀ M a_i = 1`) and is centred on `x`. Its linear
//! interpolant `Π f` is the affine function agreeing with `f` at the four
//! vertices, evaluated in barycentric coordinates as `Σ b_i f(P_i)` — no linear
//! system is solved. Then
//!
//! ```text
//!     e_K = max over the 21-point barycentric stencil of |f(Σ b_i P_i) − Σ b_i f(P_i)|.
//! ```
//!
//! The stencil is the six edge midpoints, the four face centroids, the centroid,
//! and the ten strictly interior points of the barycentric-`6` lattice
//! (`(3,1,1,1)/6` and `(2,2,1,1)/6` up to permutation). It is a **lower bound**
//! on the element's supremum error and the header says so rather than claiming
//! exactness; it is the *same* stencil on every arm, so the comparison the row
//! is about is unaffected. Edge midpoints and the centroid are in the set
//! because for a quadratic `f` the error is
//! `−½ Σ_{i<j} b_i b_j (P_i − P_j)ᵀ H (P_i − P_j)`, whose extrema over a simplex
//! sit at an edge midpoint when one pair dominates and at the centroid when the
//! six pairs are comparable.
//!
//! That identity is also recorded, as `quadratic_model_ratio`: the same maximum
//! taken over the **cell's own quadratic model** rather than over `f`, needing no
//! further samples, divided by the measured one and reduced to a median over the
//! band. It is the direct test of the hypothesis the whole metric rests on — that
//! the cell's Hessian describes `f` over the element the metric prescribed.
//!
//! ## How each norm is discretised
//!
//! Every band sample owns one cell of volume `h³`, so the `h³` cancels out of
//! `(Σ e^p h³ / Σ h³)^{1/p}` and each norm is exactly a root-mean-power over the
//! `n` band samples:
//!
//! | `norm` | discretisation | note |
//! |---|---|---|
//! | `L^1` | `(1/n) Σ_K e_K` | the arithmetic mean of the element errors |
//! | `L^2` | `((1/n) Σ_K e_K²)^{1/2}` | |
//! | `L^3` | `((1/n) Σ_K e_K³)^{1/3}` | a fourth exponent, so the trend in `p` is visible and not just its endpoints |
//! | `L^inf` | `max_K e_K` | **one** cell decides it, so `top_cell_share` is `1` by construction |
//!
//! `Ω` is the **surface band** (`|f| ≤ BAND_CELLS · h`, the same one-cell shell
//! P-146 censuses at `experiment_p146.rs:286-290`) and not the whole box. Away
//! from the surface an SDF is near-eikonal, `|H| → 0`, every eigenvalue floors,
//! and `det|H|_fl` goes constant — which, by the derivative below, makes every
//! arm tie. Integrating over the box would therefore replace the measurement
//! with `H_FLOOR`. The band is where the extractor's mesh is, and it is the
//! domain the mesh occupies that the framework integrates over.
//!
//! # Why the analytic version of C1 would be a tautology, and what the
//! # discrete one adds
//!
//! This has to be stated before the numbers, because otherwise a `true` here
//! looks like evidence and is not.
//!
//! In the continuous model `M_q = c_q(x)|H|_fl` with `c_q = s_q det|H|_fl^{−α_q}`
//! and `α_q = 1/(2q + 3)`, so `e(M_q, x) = trace(M_q^{-1/2} |H|_fl M_q^{-1/2})
//! = 3 det|H|_fl^{α_q} / s_q` exactly. Writing `D = det|H|_fl`, `u = α_q` and
//! using `q·α_q = (1 − 3u)/2`,
//!
//! ```text
//!     Ê_p(M_q) = 3 h² · A(u)^{2/3} · B(u)^{1/p},
//!     A(u) = mean D^{(1−3u)/2},    B(u) = mean D^{pu},
//!     d(ln Ê_p)/du = E_B[ln D] − E_A[ln D],
//! ```
//!
//! where `E_A`, `E_B` are the tilted means of `ln D` at exponents `(1−3u)/2` and
//! `pu`. At `u = α_p` those two exponents are both `p/(2p + 3)`, the derivative
//! is **identically zero**, and since `E_B` rises with `pu` (its derivative is
//! `Var(ln D) ≥ 0`) while `E_A` falls with `u`, the derivative is increasing and
//! the stationary point is the unique **minimum**. So `q = p` wins the continuous
//! comparison by arithmetic, strictly — unless `Var(ln D) = 0`, in which case
//! every arm ties exactly.
//!
//! `Ê_p` is computed and recorded anyway, as `continuous_error`, with
//! `continuous_wins` and `continuous_agrees`, because it is the cleanest possible
//! separation of *theory on our data* from *measurement on our data*. And
//! `det_floored_log_variance` is recorded because the derivative above says that
//! is the number which decides whether the continuous ranking can discriminate
//! at all — a field whose band is entirely at `H_FLOOR` has `D` constant, `Var =
//! 0`, and a four-way tie.
//!
//! Four things the **discrete** measurement can break that the closed form
//! cannot, and they are the vocabulary `violation_mechanism` speaks:
//!
//! 1. **`H_FLOOR`.** `D` is a floored determinant, so where a direction is
//!    genuinely flat the exponent `−α_q` acts on `1e-9` and not on a curvature.
//! 2. **The element leaves the region its Hessian describes.** `a_0 = v_0/√λ_0`
//!    with `λ_0` at the floor is a very long edge; `element_extent_relative`
//!    reports it against the domain's own extent.
//! 3. **The Riemann sum over a one-cell band** is not `∫_Ω`.
//! 4. **A twenty-one-point maximum** is not a supremum.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `M_L1` | `q = 1` — the metric optimised for `L^1` | no |
//! | `M_L2` | `q = 2` | no |
//! | `M_L3` | `q = 3` | no |
//! | `M_Linf` | `q = ∞`, where `det^(−0) = 1` and the metric is `\|H\|_fl` itself: no grading at all | **yes** — the ungraded metric is what "optimal" is optimal against |
//!
//! There is no isotropic arm and there should not be: C1 is a statement about
//! **which exponent** wins in **which norm**, and every arm shares one
//! eigenvector frame and one set of eigenvalue *ratios* because `M_q` differs
//! from `M_2` only by the scalar `c_q(x)`. So the arms differ purely in
//! **grading** — how the fixed budget is spread between high- and low-curvature
//! band samples — which is exactly the free parameter the theory picks. It also
//! means `aspect_ratio_max` is the same number on every arm of a group up to the
//! eigensolver's round-off, and it is reported per row from that row's own arm.
//!
//! # The grouping, which is C2's structure and not its prose
//!
//! C2 says the comparison is *"always reported within a norm, never across
//! two"*. That is made structural rather than promised: `group_key` is
//! `field|resolution|norm`, `group_arms` is the number of metrics compared inside
//! it, the winner is an `argmin` over that group and nothing else, and the run
//! asserts every group holds exactly one row per metric. A reader can partition
//! `p-150.csv` on `group_key` and recover every comparison this file made.
//!
//! `optimal_metric_wins`, `anomaly_reproduced`, `cao_table_agreement`, `c1_holds`
//! and `c2_holds` are therefore **group** verdicts, identical on the four rows of
//! a group; `error_measured`, `element_count`, `quadratic_model_ratio`,
//! `element_extent_*`, `complexity_*` and `metric_scale` are per-arm. The global
//! verdicts are `c1_global_holds`, `c2_global_holds` and
//! `anomaly_reproduced_global`, identical on every row.
//!
//! `cao_table_agreement` is one of three words, and the three are the only
//! outcomes his paper admits:
//!
//! | value | meaning |
//! |---|---|
//! | `theory_holds` | the metric optimised for this norm won in it |
//! | `l1_anomaly_reproduced` | it lost, **and the norm is `L^1`** — his own stated exception, reproduced on different data |
//! | `violation_outside_l1` | it lost in a norm that is not `L^1` — a violation his theory does not admit |
//!
//! `violation_mechanism` is a fixed ordered classification, evaluated only for a
//! group that violates, with every threshold a named constant: `budget_mismatch`
//! → `unresolved_tie` → `hessian_floor` → `element_exceeds_domain` →
//! `nonquadratic_element` → `tail_dominated` → `uncharacterised`. **C2 holds for
//! a group exactly when its mechanism is not `uncharacterised`**, so the one way
//! to falsify C2 is to find a violation none of those five explains.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE line is *"none — this is a methodology guard"*, and
//! that is discharged rather than skipped: **nothing here proposes a change to
//! `crates/isomesh/src/`, and that is the finding.** What the row produces is a
//! constraint on how Group D may be quoted — an error comparison between two
//! metrics is meaningless unless both are read in the same norm at the same
//! complexity — plus the two columns that let a reader check it,
//! `det_floored_log_variance` and `continuous_agrees`. A guard that shipped a
//! landing would be a guard that had drifted into being a feature.
//!
//! # Vacuity controls
//!
//! Every one runs **before** the first `run.record` and every panic starts
//! `"VOID: "`. `M-44`: a zero that could not have been non-zero is not a
//! measurement.
//!
//! - **The registered control — at least three norms must be *run*.** Counted
//!   from the groups that actually exist (`seen`), not from `LEVELS.len()`, and
//!   `L^1`, `L^2` and `L^inf` must each be among them by name. Without three,
//!   *"within a norm"* is not a constraint. Column: `group_key`, whose third
//!   field is the norm.
//! - **At least two metrics per group**, asserted as
//!   `group.errors.len() == LEVEL_COUNT ≥ 2` on every group. With one arm
//!   *"the optimal metric wins"* is a comparison against nothing and is true for
//!   free. Column: `group_arms`.
//! - **Every arm's error population must contain a positive error.** An arm whose
//!   `e_K` is zero on every band sample has an affine field on every element and
//!   admits no ranking at all. Asserted where the errors are produced. Column:
//!   `error_measured`.
//! - **The budget must actually be matched**, `complexity_rel_error ≤ 1e-8` on
//!   every arm of every rung. Unmatched, the arms are meshes of different sizes
//!   and the whole comparison is Cao's `1.4%` caveat made worse. Column:
//!   `complexity_rel_error`.
//! - **At least one group's ranking must be resolved** — best and second-best
//!   apart by more than `MIN_SEPARATION = 1e-3`. If every group is a four-way
//!   tie then both a "win" and a "violation" are arithmetic noise and C1 cannot
//!   be falsified. Columns: `runner_up_margin_relative`, `ranking_resolved`,
//!   `groups_resolved`.
//! - **The band must be non-empty** on every field and rung, or every statistic
//!   above is a statistic of nothing.
//!
//! `norm_dependent_winner_rungs` is deliberately **not** a control. A zero there
//! would mean one metric minimises every norm, which is C1 falsified in the
//! strongest possible way and a registered outcome; aborting on it would suppress
//! the finding rather than report it.
//!
//! # Determinism
//!
//! One thread, no PRNG, no map iteration, `f64` throughout. Sorting is
//! [`f64::total_cmp`]. The band is swept `z`, `y`, `x` with `x` innermost, the
//! crate's order, and every sum is taken in slice order. The only seeded objects
//! anywhere near this row are `FbmTerrain`'s and `NoiseVolume`'s committed
//! `0x5EED_1234`, which are the fields' and not the harness's. **No wall clock
//! is read and no column is machine-dependent** — the registration names no cost
//! threshold, so M-280's `1.45x` governor swing has nothing to swing here.
//!
//! Three rungs, `17 → 33 → 65` samples per axis: a factor of `4` in `h` and so
//! roughly `16x` in a second-order error, all odd (`M-266`), and the authoring
//! contract's preferred `33³` and `65³` with one coarser rung so a
//! resolution-dependent verdict is visible rather than assumed away.

#![allow(clippy::too_many_lines)]

mod common;

use isomesh::Sdf;
use isomesh::fields::ReferenceField;

use crate::common::metric::{H_FLOOR, Sym3, aspect_ratio, complexity, hessian, metric_lp};

/// Samples per axis. Odd throughout (`M-266`).
const RUNGS: [u32; 3] = [17, 33, 65];

/// `for_each_reference_field!` yields eight (`fields/mod.rs:211-256`).
const FIELDS: usize = 8;

/// How many exponents [`LEVELS`] carries.
const LEVEL_COUNT: usize = 4;

/// The exponents this row runs, as `(norm label, metric label, exponent)`.
///
/// **One list and not two**, because "the metric optimised for a given norm" is
/// then index identity and there is no lookup that can pair a norm with the
/// wrong metric — which is the single mistake that would turn C1 into noise.
const LEVELS: [(&str, &str, f64); LEVEL_COUNT] = [
    ("L^1", "M_L1", 1.0),
    ("L^2", "M_L2", 2.0),
    ("L^3", "M_L3", 3.0),
    ("L^inf", "M_Linf", f64::INFINITY),
];

/// The registration's vacuity control: this many norms must be **run**.
const MIN_NORMS: usize = 3;

/// ...and these three must be among them, by name.
const REQUIRED_NORMS: [&str; 3] = ["L^1", "L^2", "L^inf"];

/// A grid sample joins the surface band when `|f| <= BAND_CELLS · h`.
///
/// One cell, the same shell P-146 censuses (`experiment_p146.rs:286-290`), so
/// the two rows' at-floor fractions are comparable numbers rather than two
/// conventions wearing one name.
const BAND_CELLS: f64 = 1.0;

/// Two arms tie when their errors are within this, relative.
///
/// This is a **numerical** tolerance, not an evidential one: the arms share one
/// cell set and one stencil, so their difference is `f64` arithmetic over a sum
/// of at most `2.2e4` terms, whose round-off is of order `1e-13`. Three orders
/// above that is a tie that cannot be an artefact of summation order.
const RANK_RESOLUTION: f64 = 1e-9;

/// A group's ranking counts as resolved when best and second-best are apart by
/// more than this, relative.
///
/// This is the **evidential** threshold and is deliberately six orders looser
/// than [`RANK_RESOLUTION`]: a `0.1%` gap is a difference a reader can defend
/// against the stencil's own crudeness, and a smaller one is not.
const MIN_SEPARATION: f64 = 1e-3;

/// `|C(s·M)/N − 1|` may not exceed this on any arm.
///
/// `s^{3/2}` is applied and then `C` is recomputed by the same cofactor
/// determinant, so the amplified cancellation at a `1e9` aspect ratio cancels in
/// the ratio and only the two `powf` roundings survive.
const BUDGET_TOLERANCE: f64 = 1e-8;

/// Above this at-floor fraction, `det|H|_fl` is the product of `H_FLOOR`s over
/// most of the band and the continuous ranking is degenerate by construction.
const FLOOR_SHARE: f64 = 0.5;

/// `element_extent_relative` above this means the prescribed element is longer
/// than the domain, so the cell's Hessian is being asked about `f` far outside
/// the cell it was differenced in.
const EXTENT_LIMIT: f64 = 1.0;

/// `quadratic_model_ratio` outside this band means the quadratic bound `M_Lp` is
/// derived from does not predict the measured error to within a factor of two.
const MODEL_BAND: [f64; 2] = [0.5, 2.0];

/// A single band sample carrying more than this share of `Σ e^p` decides the
/// norm on its own.
const TAIL_SHARE: f64 = 0.5;

/// Cao's Table 2, `L^1` column, his optimal metric — the smaller of the pair.
const CAO_L1_OPTIMAL: f64 = 5.43e-7;

/// Cao's Table 2, `L^1` column, the metric optimised for another norm.
const CAO_L1_OTHER: f64 = 2.79e-6;

/// How closely Cao's two element counts match at `N_e ~ 16,000`: `1.4%`, which
/// the registration names as *"not exactly"*.
const CAO_ELEMENT_MATCH: f64 = 0.014;

/// `1/3` — the barycentric weight of a face centroid.
const THIRD: f64 = 1.0 / 3.0;

/// `1/6` — the unit of the interior barycentric-`6` lattice.
const SIXTH: f64 = 1.0 / 6.0;

/// How many points [`STENCIL`] carries.
const STENCIL_POINTS: usize = 21;

/// Barycentric weights of the points `e_K` is maximised over.
///
/// Six edge midpoints, four face centroids, the centroid, and the ten strictly
/// interior points of the barycentric-`6` lattice. See the header for why the
/// midpoints and the centroid have to be in the set.
const STENCIL: [[f64; 4]; STENCIL_POINTS] = [
    // six edge midpoints
    [0.5, 0.5, 0.0, 0.0],
    [0.5, 0.0, 0.5, 0.0],
    [0.5, 0.0, 0.0, 0.5],
    [0.0, 0.5, 0.5, 0.0],
    [0.0, 0.5, 0.0, 0.5],
    [0.0, 0.0, 0.5, 0.5],
    // four face centroids
    [THIRD, THIRD, THIRD, 0.0],
    [THIRD, THIRD, 0.0, THIRD],
    [THIRD, 0.0, THIRD, THIRD],
    [0.0, THIRD, THIRD, THIRD],
    // the centroid
    [0.25, 0.25, 0.25, 0.25],
    // (3,1,1,1)/6
    [0.5, SIXTH, SIXTH, SIXTH],
    [SIXTH, 0.5, SIXTH, SIXTH],
    [SIXTH, SIXTH, 0.5, SIXTH],
    [SIXTH, SIXTH, SIXTH, 0.5],
    // (2,2,1,1)/6
    [THIRD, THIRD, SIXTH, SIXTH],
    [THIRD, SIXTH, THIRD, SIXTH],
    [THIRD, SIXTH, SIXTH, THIRD],
    [SIXTH, THIRD, THIRD, SIXTH],
    [SIXTH, THIRD, SIXTH, THIRD],
    [SIXTH, SIXTH, THIRD, THIRD],
];

/// The six vertex pairs of a tetrahedron, in the order the difference vectors
/// are built in [`element_error`].
const PAIRS: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// `a − b`.
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// `‖v‖`.
fn length(v: [f64; 3]) -> f64 {
    v[0].mul_add(v[0], v[1].mul_add(v[1], v[2] * v[2])).sqrt()
}

/// `vᵀ M v`.
///
/// Three lines of arithmetic over `Sym3::get`, not a second copy of a mechanism:
/// `common::metric`'s own `quadratic_form` is private and is not part of the
/// module's API.
fn form(m: &Sym3, v: [f64; 3]) -> f64 {
    let mut sum = 0.0f64;
    for (i, &vi) in v.iter().enumerate() {
        for (j, &vj) in v.iter().enumerate() {
            sum += m.get(i, j) * vi * vj;
        }
    }
    sum
}

/// `Π max(|λᵢ|, H_FLOOR)` and whether any `|λᵢ|` sat at the floor.
///
/// This is the determinant `metric_lp`'s exponent acts on
/// (`benches/common/metric.rs:492-497`), read back out so the row can report the
/// floor rather than be quietly shaped by it.
fn spectrum(hess: &Sym3) -> (f64, bool) {
    let (values, _) = hess.eigen();
    let mut determinant = 1.0f64;
    let mut at_floor = false;
    for value in values {
        let magnitude = value.abs();
        if magnitude <= H_FLOOR {
            at_floor = true;
        }
        determinant *= magnitude.max(H_FLOOR);
    }
    (determinant, at_floor)
}

/// The three edge vectors of the unit-metric tetrahedron: `a_i = v_i / √λ_i`.
fn element_edges(metric: &Sym3) -> [[f64; 3]; 3] {
    let (values, vectors) = metric.eigen();
    let mut out = [[0.0f64; 3]; 3];
    for (column, edge) in out.iter_mut().enumerate() {
        let lambda = values[column];
        assert!(
            lambda > 0.0,
            "P-150: metric eigenvalue {lambda} is not positive, so the element it prescribes has \
             no edge length. metric_lp floors every eigenvalue at D_LP * H_FLOOR = {H_FLOOR:e} \
             and cannot produce this (benches/common/metric.rs:485-509)"
        );
        let scale = 1.0 / lambda.sqrt();
        for (row, component) in edge.iter_mut().enumerate() {
            *component = vectors[row][column] * scale;
        }
    }
    out
}

/// `(measured e_K, the same maximum over the cell's own quadratic model)`.
///
/// The model half is free: for a quadratic the error is
/// `−½ Σ_{i<j} b_i b_j (P_i − P_j)ᵀ H (P_i − P_j)`, which needs only `H` and the
/// six pair differences and no further samples of `f`.
fn element_error<F>(field: &F, centre: [f64; 3], edges: &[[f64; 3]; 3], hess: &Sym3) -> (f64, f64)
where
    F: Sdf<Scalar = f64>,
{
    // The tetrahedron's centroid is the band sample, so the corner sits at
    // `x − (a_0 + a_1 + a_2)/4`.
    let mut corner = centre;
    for edge in edges {
        for (slot, &component) in corner.iter_mut().zip(edge) {
            *slot -= 0.25 * component;
        }
    }
    let mut verts = [corner; 4];
    for (index, edge) in edges.iter().enumerate() {
        for (slot, &component) in verts[index + 1].iter_mut().zip(edge) {
            *slot += component;
        }
    }
    let values = verts.map(|vertex| field.sample(vertex));

    // `P_i − P_j` for [`PAIRS`], in that order: the three edges, then their
    // three differences.
    let differences = [
        edges[0],
        edges[1],
        edges[2],
        sub(edges[1], edges[0]),
        sub(edges[2], edges[0]),
        sub(edges[2], edges[1]),
    ];
    let curvatures = differences.map(|difference| form(hess, difference));

    let mut measured = 0.0f64;
    let mut modelled = 0.0f64;
    for weights in STENCIL {
        let mut point = [0.0f64; 3];
        let mut interpolated = 0.0f64;
        for (index, &weight) in weights.iter().enumerate() {
            for (slot, component) in point.iter_mut().zip(verts[index]) {
                *slot += weight * component;
            }
            interpolated += weight * values[index];
        }
        measured = measured.max((field.sample(point) - interpolated).abs());

        let mut quadratic = 0.0f64;
        for (&(i, j), &curvature) in PAIRS.iter().zip(curvatures.iter()) {
            quadratic += weights[i] * weights[j] * curvature;
        }
        modelled = modelled.max((0.5 * quadratic).abs());
    }

    (measured, modelled)
}

/// `‖e‖_p` as a root-mean-power over the band, and the largest single sample's
/// share of `Σ e^p`.
///
/// Every band sample owns the same `h³`, so the cell volume cancels out of the
/// Riemann sum and this is exact rather than an approximation of the header's
/// formula. `common::metric`'s `l_tau_norm` is private, so this is written here.
fn norm_of(errors: &[f64], p: f64) -> (f64, f64) {
    assert!(
        !errors.is_empty(),
        "P-150: a norm over an empty band is unasked, not zero"
    );
    if p.is_infinite() {
        let mut worst = 0.0f64;
        for &error in errors {
            worst = worst.max(error);
        }
        return (worst, 1.0);
    }
    let mut sum = 0.0f64;
    let mut worst = 0.0f64;
    for &error in errors {
        let term = error.powf(p);
        sum += term;
        worst = worst.max(term);
    }
    let mean = sum / errors.len() as f64;
    (mean.powf(1.0 / p), worst / sum)
}

/// Median of a non-empty population, by [`f64::total_cmp`].
fn median(mut xs: Vec<f64>) -> f64 {
    assert!(!xs.is_empty(), "P-150: median of an empty population");
    xs.sort_by(f64::total_cmp);
    xs[xs.len() / 2]
}

/// `Var(ln D)` over the band — the quantity the header's derivative names as the
/// discriminating power of the whole continuous comparison.
fn log_variance(values: &[f64]) -> f64 {
    let count = values.len() as f64;
    let mut mean = 0.0f64;
    for &value in values {
        mean += value.ln();
    }
    mean /= count;
    let mut variance = 0.0f64;
    for &value in values {
        let deviation = value.ln() - mean;
        variance += deviation * deviation;
    }
    variance / count
}

/// Index of the smallest entry; the first on a tie, so the answer is one
/// ordering rather than whichever the comparison happened to leave.
fn argmin(values: &[f64]) -> usize {
    let mut best = 0usize;
    for (index, value) in values.iter().enumerate() {
        if *value < values[best] {
            best = index;
        }
    }
    best
}

/// Second-smallest entry of a slice of at least two.
fn second_smallest(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    sorted[1]
}

/// The grid samples of the surface band, swept `z`, `y`, `x` with `x` innermost.
fn band_samples<F>(field: &F, origin: [f64; 3], h: f64, samples: u32) -> Vec<[f64; 3]>
where
    F: Sdf<Scalar = f64>,
{
    let band = BAND_CELLS * h;
    let mut out = Vec::new();
    for k in 0..samples {
        for j in 0..samples {
            for i in 0..samples {
                let point = [
                    origin[0] + f64::from(i) * h,
                    origin[1] + f64::from(j) * h,
                    origin[2] + f64::from(k) * h,
                ];
                if field.sample(point).abs() <= band {
                    out.push(point);
                }
            }
        }
    }
    out
}

/// What one `(field, resolution)` hands to every arm, unchanged.
struct RungInput<'a> {
    /// The band's grid samples.
    points: &'a [[f64; 3]],
    /// One Hessian per band sample, differenced at the cell size.
    hessians: &'a [Sym3],
    /// Whether that sample has a floored `|H|` eigenvalue.
    at_floor: &'a [bool],
    /// `h³`, the volume each band sample owns.
    cell_volume: f64,
    /// The complexity every arm is scaled to: the band's sample count.
    target: f64,
    /// `hi[0] − lo[0]`, the domain's own extent.
    extent: f64,
}

/// One exponent's arm on one rung.
struct Arm {
    /// The `metric` column: `M_L1`, `M_L2`, `M_L3` or `M_Linf`.
    label: &'static str,
    /// `q`.
    exponent: f64,
    /// The `D_Lp` this arm needed to hit the budget: `(N / C(M_q))^{2/3}`.
    scale: f64,
    /// `C(s·M_q)`, recomputed after scaling.
    complexity_measured: f64,
    /// `|C(s·M_q)/N − 1|`.
    complexity_rel_error: f64,
    /// Longest prescribed element edge, in world units.
    extent_max: f64,
    /// The same, over the domain's extent.
    extent_relative: f64,
    /// `max` and off-floor `max` of `metric::aspect_ratio` over the band.
    aspect_max: f64,
    /// `aspect_max` restricted to band samples with no floored eigenvalue.
    aspect_max_off_floor: f64,
    /// Median of `modelled / measured` over the band samples with a positive
    /// measured error.
    model_ratio: f64,
    /// `e_K` per band sample, in band order.
    errors: Vec<f64>,
}

/// One `(field, resolution, norm)` — the group C2 requires every comparison to
/// happen inside.
struct Group {
    /// The `norm` column.
    norm: &'static str,
    /// `p`.
    exponent: f64,
    /// `E_p(M_q)` per arm, in [`LEVELS`] order.
    errors: Vec<f64>,
    /// The largest band sample's share of `Σ e^p`, per arm.
    top_share: Vec<f64>,
    /// `Ê_p(M_q)`, the continuous-model error, per arm.
    continuous: Vec<f64>,
    /// `argmin` over `errors`.
    winner: usize,
    /// `argmin` over `continuous`.
    continuous_winner: usize,
    /// Did the metric optimised for **this** norm win in it?
    wins: bool,
    /// The same question asked of the closed form.
    continuous_wins: bool,
    /// `(own − best) / best`; zero when the own arm wins.
    loss_margin: f64,
    /// `(second-best − best) / best`.
    runner_up_margin: f64,
    /// `runner_up_margin > MIN_SEPARATION`.
    resolved: bool,
    /// `max(errors) / min(errors)` — the number directly comparable to Cao's
    /// `5.14x`.
    spread: f64,
    /// The ordered classification of a violation; `none` when there is none.
    mechanism: &'static str,
}

/// Everything one `(field, resolution)` produced.
struct Rung {
    /// Samples per axis.
    samples: u32,
    /// Band sample count.
    band: usize,
    /// `h`.
    cell_size: f64,
    /// `hi[0] − lo[0]`.
    extent: f64,
    /// Band samples with a floored `|H|` eigenvalue.
    at_floor: usize,
    /// `Var(ln det|H|_fl)` over the band.
    det_log_variance: f64,
    /// One per [`LEVELS`] entry.
    arms: Vec<Arm>,
    /// One per [`LEVELS`] entry.
    groups: Vec<Group>,
}

/// Build and measure one arm.
fn measure_arm<F>(field: &F, input: &RungInput<'_>, level: (&str, &'static str, f64)) -> Arm
where
    F: Sdf<Scalar = f64>,
{
    let (_, label, exponent) = level;

    let unscaled: Vec<Sym3> = input
        .hessians
        .iter()
        .map(|hess| metric_lp(hess, exponent))
        .collect();
    let raw = complexity(&unscaled, input.cell_volume);
    let scale = (input.target / raw).powf(2.0 / 3.0);
    let metrics: Vec<Sym3> = unscaled.iter().map(|metric| metric.scale(scale)).collect();
    let complexity_measured = complexity(&metrics, input.cell_volume);

    let mut errors = Vec::with_capacity(input.points.len());
    let mut ratios = Vec::with_capacity(input.points.len());
    let mut extent_max = 0.0f64;
    let mut aspect_max = 0.0f64;
    let mut aspect_max_off_floor = 0.0f64;

    for (index, metric) in metrics.iter().enumerate() {
        let aspect = aspect_ratio(metric);
        aspect_max = aspect_max.max(aspect);
        if !input.at_floor[index] {
            aspect_max_off_floor = aspect_max_off_floor.max(aspect);
        }

        let edges = element_edges(metric);
        for &edge in &edges {
            extent_max = extent_max.max(length(edge));
        }

        let (measured, modelled) =
            element_error(field, input.points[index], &edges, &input.hessians[index]);
        if measured > 0.0 {
            ratios.push(modelled / measured);
        }
        errors.push(measured);
    }

    assert!(
        !ratios.is_empty(),
        "VOID: arm {label} produced e_K = 0 on every one of {} band samples, so the field is \
         affine on every element the metric prescribed, there is no interpolation error to rank \
         and every zero in this arm's column is a zero that could not have been non-zero (M-44)",
        input.points.len()
    );

    Arm {
        label,
        exponent,
        scale,
        complexity_measured,
        complexity_rel_error: (complexity_measured / input.target - 1.0).abs(),
        extent_max,
        extent_relative: extent_max / input.extent,
        aspect_max,
        aspect_max_off_floor,
        model_ratio: median(ratios),
        errors,
    }
}

/// The ordered classification of one group's violation. `none` when it does not
/// violate; `uncharacterised` is the one value that falsifies C2.
fn mechanism_of(group_wins: bool, arms: &[Arm], group: &GroupDraft<'_>) -> &'static str {
    if group_wins {
        return "none";
    }
    let budget = arms
        .iter()
        .map(|arm| arm.complexity_rel_error)
        .fold(0.0f64, f64::max);
    let own = &arms[group.own];
    let winner = &arms[group.winner];

    if budget > BUDGET_TOLERANCE {
        "budget_mismatch"
    } else if !group.resolved {
        "unresolved_tie"
    } else if group.at_floor_fraction > FLOOR_SHARE {
        "hessian_floor"
    } else if own.extent_relative.max(winner.extent_relative) > EXTENT_LIMIT {
        "element_exceeds_domain"
    } else if own.model_ratio < MODEL_BAND[0] || own.model_ratio > MODEL_BAND[1] {
        "nonquadratic_element"
    } else if group.top_share_own.max(group.top_share_winner) > TAIL_SHARE {
        "tail_dominated"
    } else {
        "uncharacterised"
    }
}

/// The part of a group [`mechanism_of`] needs, assembled before the verdict.
struct GroupDraft<'a> {
    /// Index of the arm optimised for this group's norm.
    own: usize,
    /// Index of the arm that actually minimised the norm.
    winner: usize,
    /// `runner_up_margin > MIN_SEPARATION`.
    resolved: bool,
    /// Band samples at the floor, over the band.
    at_floor_fraction: f64,
    /// The own arm's largest-sample share of `Σ e^p`.
    top_share_own: f64,
    /// The winning arm's.
    top_share_winner: f64,
    /// Borrow marker so the struct carries the rung's lifetime.
    norm: &'a str,
}

/// One group per [`LEVELS`] entry, each comparing every arm **inside** one norm.
fn build_groups(arms: &[Arm], det_floored: &[f64], at_floor_fraction: f64) -> Vec<Group> {
    // `det|H|_fl^{α_q}` per arm, hoisted out of the norm loop: the tilt depends
    // on the arm and not on the norm.
    let tilted: Vec<Vec<f64>> = arms
        .iter()
        .map(|arm| {
            let alpha = 1.0 / (2.0 * arm.exponent + 3.0);
            det_floored
                .iter()
                .map(|determinant| determinant.powf(alpha))
                .collect()
        })
        .collect();

    let mut out = Vec::with_capacity(LEVEL_COUNT);
    for (own, &(norm, _, p)) in LEVELS.iter().enumerate() {
        let mut errors = Vec::with_capacity(LEVEL_COUNT);
        let mut top_share = Vec::with_capacity(LEVEL_COUNT);
        let mut continuous = Vec::with_capacity(LEVEL_COUNT);
        for (arm, tilt) in arms.iter().zip(tilted.iter()) {
            let (value, share) = norm_of(&arm.errors, p);
            errors.push(value);
            top_share.push(share);
            // Ê_p(M_q) = (3 / s_q) · ‖det|H|_fl^{α_q}‖_p, from the header.
            continuous.push(3.0 / arm.scale * norm_of(tilt, p).0);
        }

        let winner = argmin(&errors);
        let best = errors[winner];
        let own_error = errors[own];
        let wins = own_error <= best * (1.0 + RANK_RESOLUTION);
        let loss_margin = (own_error - best) / best;
        let runner_up_margin = (second_smallest(&errors) - best) / best;
        let resolved = runner_up_margin > MIN_SEPARATION;
        let worst = errors.iter().copied().fold(0.0f64, f64::max);
        let continuous_winner = argmin(&continuous);
        let continuous_wins =
            continuous[own] <= continuous[continuous_winner] * (1.0 + RANK_RESOLUTION);

        let draft = GroupDraft {
            own,
            winner,
            resolved,
            at_floor_fraction,
            top_share_own: top_share[own],
            top_share_winner: top_share[winner],
            norm,
        };
        let mechanism = mechanism_of(wins, arms, &draft);

        out.push(Group {
            norm: LEVELS[own].0,
            exponent: p,
            errors,
            top_share,
            continuous,
            winner,
            continuous_winner,
            wins,
            continuous_wins,
            loss_margin,
            runner_up_margin,
            resolved,
            spread: worst / best,
            mechanism,
        });
        // `draft.norm` is the same string the group carries; asserting it here
        // is the cheapest possible proof that the group was assembled for the
        // norm it reports, which is C2's structural half.
        assert_eq!(
            draft.norm, out[own].norm,
            "P-150: group {own} was assembled for a different norm than it reports"
        );
    }
    out
}

/// Measure one `(field, resolution)`.
fn measure_rung<F>(field: &F, name: &str, samples: u32) -> Rung
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (_, origin, h) = common::grid::<f64, _>(field, samples);
    let (domain_lo, domain_hi) = field.domain();
    let extent = domain_hi[0] - domain_lo[0];

    let points = band_samples(field, origin, h, samples);
    assert!(
        !points.is_empty(),
        "VOID: {name} at {samples}^3 put no grid sample within {BAND_CELLS} cell of its surface, \
         so every error norm, aspect ratio and at-floor count below is a statistic of an empty \
         population and every zero among them a zero that could not have been non-zero (M-44)"
    );

    // One Hessian and one floored spectrum per band sample, read once and shared
    // by every arm: four arms differencing the same field at the same step would
    // be four answers to one question.
    let mut hessians = Vec::with_capacity(points.len());
    let mut det_floored = Vec::with_capacity(points.len());
    let mut at_floor_flags = Vec::with_capacity(points.len());
    let mut at_floor = 0usize;
    for &point in &points {
        let hess = hessian(field, point, h);
        let (determinant, floored) = spectrum(&hess);
        if floored {
            at_floor += 1;
        }
        det_floored.push(determinant);
        at_floor_flags.push(floored);
        hessians.push(hess);
    }

    let input = RungInput {
        points: &points,
        hessians: &hessians,
        at_floor: &at_floor_flags,
        cell_volume: h * h * h,
        target: points.len() as f64,
        extent,
    };

    let arms: Vec<Arm> = LEVELS
        .iter()
        .map(|&level| measure_arm(field, &input, level))
        .collect();
    let at_floor_fraction = at_floor as f64 / points.len() as f64;
    let groups = build_groups(&arms, &det_floored, at_floor_fraction);

    Rung {
        samples,
        band: points.len(),
        cell_size: h,
        extent,
        at_floor,
        det_log_variance: log_variance(&det_floored),
        arms,
        groups,
    }
}

/// Measure one reference field across the whole ladder.
fn measure_field<F>(field: &F, name: &str) -> Vec<Rung>
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let mut out = Vec::with_capacity(RUNGS.len());
    for samples in RUNGS {
        let rung = measure_rung(field, name, samples);
        let violations = rung.groups.iter().filter(|group| !group.wins).count();
        println!(
            "  {name:<15} {:>3}^3  band {:>6}  at_floor {:.4}  Var(ln det) {:>10.4}  \
             violations {violations}/{}",
            rung.samples,
            rung.band,
            rung.at_floor as f64 / rung.band as f64,
            rung.det_log_variance,
            rung.groups.len()
        );
        out.push(rung);
    }
    out
}

/// Which of Cao's three outcomes this group landed on.
fn agreement(norm: &str, wins: bool) -> &'static str {
    if wins {
        "theory_holds"
    } else if norm == REQUIRED_NORMS[0] {
        "l1_anomaly_reproduced"
    } else {
        "violation_outside_l1"
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-150");

    common::experiment::run(prereg, |run| {
        println!(
            "construction: within-a-norm ranking of {LEVEL_COUNT} L^p metrics at MATCHED \
             complexity.\n  e_K = max over {STENCIL_POINTS} barycentric points of |f - Pi f| on \
             the unit-metric tetrahedron\n  norms {:?}, band |f| <= {BAND_CELLS}h, H_FLOOR = \
             {H_FLOOR:e}, ladder {RUNGS:?} samples/axis\n  Cao Table 2 L^1: {CAO_L1_OPTIMAL:e} vs \
             {CAO_L1_OTHER:e} = {:.4}x at element counts matching to {CAO_ELEMENT_MATCH}\n",
            LEVELS.map(|level| level.0),
            CAO_L1_OTHER / CAO_L1_OPTIMAL
        );

        let mut fields: Vec<(&'static str, Vec<Rung>)> = Vec::with_capacity(FIELDS);
        isomesh::for_each_reference_field!(f64, |name, field| {
            fields.push((name, measure_field(&field, name)));
        });
        assert_eq!(
            fields.len(),
            FIELDS,
            "P-150: for_each_reference_field! must yield {FIELDS} fields"
        );

        // ── vacuity controls ─────────────────────────────────────────────────
        // The registered one: at least three norms must have been RUN, counted
        // from the groups that exist rather than from LEVELS.len().
        let mut seen: Vec<&str> = Vec::new();
        for (_, rungs) in &fields {
            for rung in rungs {
                for group in &rung.groups {
                    if !seen.contains(&group.norm) {
                        seen.push(group.norm);
                    }
                }
            }
        }
        assert!(
            seen.len() >= MIN_NORMS,
            "VOID: only {} norm(s) were run ({seen:?}), below the registration's own control of \
             {MIN_NORMS}. With fewer than three, 'the comparison is always reported within a \
             norm' is not a constraint on anything",
            seen.len()
        );
        for required in REQUIRED_NORMS {
            assert!(
                seen.contains(&required),
                "VOID: {required} was not among the norms run ({seen:?}), and the registration \
                 names L^1, L^2 and L^inf as the minimum. L^1 in particular is the whole subject \
                 of this row: it is the column Cao's own theory fails in"
            );
        }

        // Every group must compare at least two metrics, and the arms must be
        // budget-matched, or a "win" is free. Counted from the groups that were
        // built, for the same reason the norm count is: a control read off a
        // constant is a control that cannot fail.
        let thinnest_group = fields
            .iter()
            .flat_map(|(_, rungs)| rungs.iter())
            .flat_map(|rung| rung.groups.iter())
            .map(|group| group.errors.len())
            .min()
            .expect("at least one group was built, or the field loop above did nothing");
        assert!(
            thinnest_group >= 2,
            "VOID: the thinnest group in the sweep compares {thinnest_group} metric(s), so \
             'the metric optimised for a given norm wins in that norm' is a comparison against \
             nothing and is true for free"
        );
        for (name, rungs) in &fields {
            for rung in rungs {
                assert_eq!(
                    rung.groups.len(),
                    LEVEL_COUNT,
                    "P-150: {name} at {}^3 built {} groups, not {LEVEL_COUNT}",
                    rung.samples,
                    rung.groups.len()
                );
                for group in &rung.groups {
                    assert_eq!(
                        group.errors.len(),
                        LEVEL_COUNT,
                        "VOID: {name} at {}^3, norm {}, compared {} arms rather than \
                         {LEVEL_COUNT}. 'The optimal metric wins in that norm' is a comparison \
                         with nothing when a group holds one arm",
                        rung.samples,
                        group.norm,
                        group.errors.len()
                    );
                }
                for arm in &rung.arms {
                    assert!(
                        arm.complexity_rel_error <= BUDGET_TOLERANCE,
                        "VOID: {name} at {}^3, arm {}, landed C(s.M) at {:.6e} against a target \
                         of {} -- a relative miss of {:.3e}, above {BUDGET_TOLERANCE:e}. The arms \
                         are then meshes of different sizes and comparing their errors is exactly \
                         the mismatch the registration flags in Cao's own table",
                        rung.samples,
                        arm.label,
                        arm.complexity_measured,
                        rung.band,
                        arm.complexity_rel_error
                    );
                }
            }
        }

        // At least one group's ranking has to be resolved, or every verdict in
        // the file is a tie wearing a verdict's clothes.
        let groups_total: usize = fields
            .iter()
            .map(|(_, rungs)| rungs.iter().map(|rung| rung.groups.len()).sum::<usize>())
            .sum();
        let groups_resolved: usize = fields
            .iter()
            .map(|(_, rungs)| {
                rungs
                    .iter()
                    .map(|rung| rung.groups.iter().filter(|group| group.resolved).count())
                    .sum::<usize>()
            })
            .sum();
        let widest = fields
            .iter()
            .flat_map(|(_, rungs)| rungs.iter())
            .flat_map(|rung| rung.groups.iter())
            .map(|group| group.runner_up_margin)
            .fold(0.0f64, f64::max);
        assert!(
            widest > MIN_SEPARATION,
            "VOID: the widest gap between the best and second-best metric in any of \
             {groups_total} groups is {widest:.3e}, below {MIN_SEPARATION:e}. Every ranking in \
             the file is then a tie, a 'win' costs nothing and a 'violation' is arithmetic noise, \
             so C1 could not have been falsified (M-44)"
        );

        // ── the verdicts ─────────────────────────────────────────────────────
        let violations: Vec<(&str, u32, &str, f64, &str)> = fields
            .iter()
            .flat_map(|(name, rungs)| {
                rungs.iter().flat_map(move |rung| {
                    rung.groups
                        .iter()
                        .filter(|group| !group.wins)
                        .map(move |group| {
                            (
                                *name,
                                rung.samples,
                                group.norm,
                                group.loss_margin,
                                group.mechanism,
                            )
                        })
                })
            })
            .collect();
        let violations_in_l1 = violations
            .iter()
            .filter(|entry| entry.2 == REQUIRED_NORMS[0])
            .count();
        let violations_outside_l1 = violations.len() - violations_in_l1;
        let c1_global = violations.is_empty();
        let anomaly_global = violations_in_l1 > 0;
        let c2_global = violations.iter().all(|entry| entry.4 != "uncharacterised");

        // How many rungs had the winner change with the norm. A characterisation
        // column, not a control: see the header.
        let norm_dependent: usize = fields
            .iter()
            .flat_map(|(_, rungs)| rungs.iter())
            .filter(|rung| {
                let first = rung.groups[0].winner;
                rung.groups.iter().any(|group| group.winner != first)
            })
            .count();

        println!(
            "\nC1: the metric optimised for a norm won in it on {} of {groups_total} groups -> \
             {c1_global}",
            groups_total - violations.len()
        );
        println!(
            "C2: {} violation(s), {violations_in_l1} in L^1 and {violations_outside_l1} outside \
             it; every one characterised -> {c2_global}",
            violations.len()
        );
        println!(
            "anomaly_reproduced (Cao's L^1 exception, on our fields) -> {anomaly_global}\n\
             groups_resolved {groups_resolved}/{groups_total}, \
             norm_dependent_winner_rungs {norm_dependent}/{}",
            fields.len() * RUNGS.len()
        );
        for (name, samples, norm, margin, mechanism) in &violations {
            println!("  violation: {name} {samples}^3 {norm} loses by {margin:.4e} -- {mechanism}");
        }
        println!();

        // ── the rows ─────────────────────────────────────────────────────────
        for (name, rungs) in &fields {
            for rung in rungs {
                let at_floor_fraction = rung.at_floor as f64 / rung.band as f64;
                for group in &rung.groups {
                    let agreement = agreement(group.norm, group.wins);
                    for (index, arm) in rung.arms.iter().enumerate() {
                        run.record(&[
                            ("norm", group.norm.to_string()),
                            ("metric", arm.label.to_string()),
                            ("element_count", format!("{:.0}", arm.complexity_measured)),
                            ("error_measured", format!("{:.6e}", group.errors[index])),
                            ("optimal_metric_wins", group.wins.to_string()),
                            (
                                "anomaly_reproduced",
                                (!group.wins && group.norm == REQUIRED_NORMS[0]).to_string(),
                            ),
                            ("cao_table_agreement", agreement.to_string()),
                            ("c1_holds", group.wins.to_string()),
                            (
                                "c2_holds",
                                (group.mechanism != "uncharacterised").to_string(),
                            ),
                            // ── extras (M-273) ──
                            ("anomaly_reproduced_global", anomaly_global.to_string()),
                            ("aspect_ratio_max", format!("{:.6e}", arm.aspect_max)),
                            (
                                "aspect_ratio_max_off_floor",
                                format!("{:.6e}", arm.aspect_max_off_floor),
                            ),
                            ("at_floor_cells", rung.at_floor.to_string()),
                            ("at_floor_fraction", format!("{at_floor_fraction:.6}")),
                            ("band_samples", rung.band.to_string()),
                            (
                                "cao_element_match_reference",
                                format!("{CAO_ELEMENT_MATCH:.6}"),
                            ),
                            (
                                "cao_spread_reference",
                                format!("{:.6}", CAO_L1_OTHER / CAO_L1_OPTIMAL),
                            ),
                            ("c1_global_holds", c1_global.to_string()),
                            ("c2_global_holds", c2_global.to_string()),
                            ("cell_size", format!("{:.9}", rung.cell_size)),
                            (
                                "complexity_measured",
                                format!("{:.6e}", arm.complexity_measured),
                            ),
                            (
                                "complexity_rel_error",
                                format!("{:.3e}", arm.complexity_rel_error),
                            ),
                            ("complexity_target", rung.band.to_string()),
                            (
                                "continuous_agrees",
                                (group.continuous_wins == group.wins).to_string(),
                            ),
                            (
                                "continuous_error",
                                format!("{:.6e}", group.continuous[index]),
                            ),
                            ("continuous_wins", group.continuous_wins.to_string()),
                            (
                                "continuous_winner_metric",
                                LEVELS[group.continuous_winner].1.to_string(),
                            ),
                            (
                                "det_floored_log_variance",
                                format!("{:.6e}", rung.det_log_variance),
                            ),
                            ("domain_extent", format!("{:.6}", rung.extent)),
                            ("element_extent_max", format!("{:.6e}", arm.extent_max)),
                            (
                                "element_extent_relative",
                                format!("{:.6e}", arm.extent_relative),
                            ),
                            ("field", (*name).to_string()),
                            ("group_arms", rung.arms.len().to_string()),
                            (
                                "group_key",
                                format!("{name}|{}|{}", rung.samples, group.norm),
                            ),
                            ("groups_resolved", groups_resolved.to_string()),
                            ("groups_total", groups_total.to_string()),
                            ("loss_margin_relative", format!("{:.6e}", group.loss_margin)),
                            ("metric_exponent", format!("{}", arm.exponent)),
                            ("metric_scale", format!("{:.6e}", arm.scale)),
                            ("norm_dependent_winner_rungs", norm_dependent.to_string()),
                            ("p_exponent", format!("{}", group.exponent)),
                            ("quadratic_model_ratio", format!("{:.6e}", arm.model_ratio)),
                            ("ranking_resolved", group.resolved.to_string()),
                            ("resolution", rung.samples.to_string()),
                            (
                                "runner_up_margin_relative",
                                format!("{:.6e}", group.runner_up_margin),
                            ),
                            ("spread_in_norm", format!("{:.6}", group.spread)),
                            ("stencil_points", STENCIL_POINTS.to_string()),
                            ("top_cell_share", format!("{:.6}", group.top_share[index])),
                            ("violation_mechanism", group.mechanism.to_string()),
                            ("violations_in_l1", violations_in_l1.to_string()),
                            ("violations_outside_l1", violations_outside_l1.to_string()),
                            ("violations_total", violations.len().to_string()),
                            ("winner_metric", LEVELS[group.winner].1.to_string()),
                        ]);
                    }
                }
            }
        }
    });
}
