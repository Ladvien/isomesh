//! **P-152 — the `beta`-number is the planarity assumption the QEF never checks.**
//!
//! Ticket: R-152. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p152
//! ```
//!
//! Writes `docs/experiments/p-152.csv`.
//!
//! # What was missing
//!
//! Dual Contouring places a cell's vertex where the tangent planes of its edge
//! crossings best agree. The shipped rule is four lines
//! (`crates/isomesh/src/dual_contouring/solve.rs:7-11`):
//!
//! ```text
//! c = centroid of the crossing positions        d_i = n_i.(p_i - c)
//! M = sum n_i n_i^T        g = sum d_i n_i        lambda = 0.01
//! x = c + adj(M + lambda*I).g / det(M + lambda*I)
//! ```
//!
//! **`solve` returns the point and throws the objective away** (`solve.rs:248`,
//! `:272`). Nothing in `crates/isomesh/src/**` computes the value of the
//! quadratic it just minimised — a `grep` for `residual` over the whole of
//! `src/` finds only `subgrid`'s residual *points*, an unrelated integer
//! construction at `subgrid/curves.rs:47-51`. So the crate has shipped a
//! least-squares plane fit for its entire life and has never once recorded
//! **how badly the planes disagreed**.
//!
//! Nor has it ever asked, per cell, whether the fit was determined by the data
//! at all. `lambda = 0.01` (`solve.rs:85`) is described exactly as *"the
//! Tikhonov regularizer that stops an under-determined cell — a flat region,
//! where `M` is rank 1 — from flying off"* (`solve.rs:44-47`), and `✗12` records
//! the deliberate removal of the rank branch that used to decide it
//! (`solve.rs:20-34`: the branch disagreed after a rotation in 454 of 20,000
//! trials and moved the vertex a median of 2.13 cells when it flipped). The
//! constant is therefore a *fixed answer* to a question — is this cell
//! under-determined? — that the crate has never evaluated. `singular_cells`
//! below is that evaluation, and it uses the crate's own number as the bar.
//!
//! `V-18` is the third part of the same gap: `M = A^T A` squares the condition
//! number, DC's own paper measures `b^T b` reaching `~1e6` on a `256^3` grid,
//! and that is why `Real` spans `f64` at all (`solve.rs:59-71`). A
//! condition-number argument about `M` with no per-cell measurement of `M`'s
//! smallest eigenvalue is an argument about a matrix nobody has looked at.
//!
//! What Jones' `beta_inf(Q)` adds is a *name with theorems* for the assumption
//! the fit makes. `benches/common/beta.rs` — R-152's own module, consumed
//! unchanged here — computes it as the half-width of the thinnest slab
//! containing the cell's patch over `diam Q`, and
//! [`beta_per_cell`](common::beta::beta_per_cell) uses the cell's own eight
//! corners and twelve edges: **byte-for-byte the point set the QEF is fitted
//! to** (`beta.rs:72-75`, `:421-425`). That is what makes this row a comparison
//! of two statistics of *one* point set rather than of two point sets. The
//! surface case of the Traveling Salesman Theorem is Azzam & Schul,
//! `arXiv:1609.02892`; the QEF is Ju, Losasso, Schaefer & Warren,
//! `10.1145/566570.566586` §2.3, cited at `hermite.rs:12-16`.
//!
//! `beta` is reported as an **upper bound** on Jones' infimum, at a budget of
//! [`SLAB_EVALUATIONS`](common::beta::SLAB_EVALUATIONS)` = 73` objective
//! evaluations per patch (`beta.rs:37-62`). Every `beta` in this CSV is
//! `>=` the true `beta_inf`, never below it, so it can only over-report
//! roughness — the honest direction for a flatness statistic, and the reason the
//! budget is recorded beside it in `slab_evaluations_per_patch`.
//!
//! # The QEF residual, and why it is the shipped one rather than a bench copy
//!
//! The brief this harness was cut from assumed the crate's solver is private and
//! that the residual would need a bench-local `3x3` normal-equations solve.
//! **It is not private.** `isomesh::hermite::HermiteCell` and
//! `isomesh::dual_contouring::solve::{LAMBDA, solve}` are both public and both
//! already used from benches (`experiment_p52.rs:217-220`), so this file
//! measures the residual of the vertex the *shipped* path actually places rather
//! than the residual of a second solver that would have to be trusted:
//!
//! ```text
//! corners  -> HermiteCell::from_corners(field, &corners, cell_origin, h)
//! x        -> solve(&cell)                       // solve.rs:248, lambda = 0.01
//! residual -> sum over crossings of (n_i . (x - p_i))^2
//! ```
//!
//! Three properties of that make it the right quantity and not merely an
//! available one:
//!
//! - **It is the data misfit and excludes the prior.** `solve` minimises
//!   `sum (n_i.(x - p_i))^2 + lambda*|x - c|^2`; the Tikhonov term is a
//!   statement about where the vertex ought to be, not about how well the planes
//!   agree, so it is not in the residual. C2 pairs `beta` against the
//!   *plane-fit* residual and that is what this is.
//! - **It is the unclamped solve.** `solve` is documented as *not* clamped
//!   (`solve.rs:242-246`); `Clamp::ToCell` is applied afterwards by the rule
//!   (`dual_contouring.rs:168`, `:212`). A residual measured at a clamped vertex
//!   would be the residual of a point the QEF did not choose.
//! - **The population is literally `beta`'s population.** `HermiteCell` cuts an
//!   edge on `is_inside(a) != is_inside(b)` with `is_inside(v) = v < 0`
//!   (`cube.rs:171-173`), which is character-for-character `beta.rs:622`'s test,
//!   over the same twelve cube edges. `population_mismatch` is asserted **0** on
//!   every row, and `HermiteCell::len()` is asserted against a crossing count
//!   derived here from the public `EDGE_CORNERS` — the check
//!   `bevy_isomesh/examples/hermite_debug.rs:146-152` exists for, because
//!   `from_corners`' corner order lives in a private module and a permuted guess
//!   would silently land the crossings on different edges.
//!
//! The two crossing *positions* differ in their last bits and nowhere else:
//! `beta.rs:626-630` uses `t = a/(a - b)` from the low corner and the crate uses
//! the centred `d = ((a + b)/2)/(a - b)` that `P-61` installed for exact
//! reflection antisymmetry (`cube.rs:175-225`). Algebraically
//! `0.5 + d = a/(a - b)`, so the point sets are the same set and only the
//! rounding of the last ulp differs. No clause here is decided at that scale —
//! but the *planarity floor* below is, so it is measured rather than assumed.
//!
//! # `beta` on a planar patch is rounding noise, not zero — and the gap is eight decades
//!
//! `beta.rs:84-99` is right that a three-crossing cell's patch is planar *by
//! construction*: three points are coplanar, so the true `beta_inf` is exactly
//! `0`. It calls the computed value *"zero to the rounding floor"*, which is the
//! accurate phrase — the value is **not** exactly `0.0`. `thinnest_slab` takes
//! the eigenvector of a covariance matrix by cyclic Jacobi and then a
//! max-minus-min of dot products, and both steps round.
//!
//! Measured out-of-tree on this exact code path over all sixteen
//! `(field, resolution)` pairs, the `beta_inf` column is **bimodal with a seven
//! decade dead band**: every occupied decade is at or below `1e-16`, or at or
//! above `1e-8`. Decades `1e-15` through `1e-9` are empty on every field at both
//! resolutions. So
//!
//! ```text
//! PLANAR_FLOOR = 1e-12
//! ```
//!
//! is not a tolerance invented to make a clause pass; it is the middle of a band
//! no cell occupies, three decades above the largest noise and four below the
//! smallest signal. Two columns re-derive that on every run — `beta_max_below_floor`
//! and `beta_min_above_floor` — and a control asserts the separation is at least
//! three decades wide, so a future field that lands *in* the band fails the run
//! instead of quietly having its ranks decided by Jacobi's last bits.
//!
//! **Why this matters, and why it is not bookkeeping: the floor decides C2.**
//! Without it, 96% of `sphere`'s surface cells carry a "non-zero" `beta` of order
//! `1e-17` and `rank_correlation` spends most of its ranks ordering rounding
//! noise. Measured out-of-tree on this code path at `65^3`, the unfloored column
//! gives `thin_plate` a coefficient of **0.996** — over a `beta` column whose
//! every entry lies below `1e-15`, on a field `beta.rs:90` itself calls planar in
//! all 512 of its cells — and pushes `csg_difference` from below C2's bar to
//! above it. Six of eight fields clear `0.7` unfloored and four do floored, so
//! the two readings return **opposite verdicts on C2**, and one of them is a
//! verdict about Jacobi's last bits.
//!
//! The registered `rank_correlation` is therefore computed on the **floored**
//! column, where every sub-floor cell is collapsed to one honest tie block.
//! `rank_correlation_unfloored` is recorded beside it as the sensitivity strip,
//! in the shape `experiment_p124.rs:307-313` uses, so the entry must name which
//! reading it quotes rather than inherit one silently.
//!
//! # A correction to `common::beta`'s own census, reported rather than edited
//!
//! `beta.rs:89-91` states: *"of the surface cells, 944 of 1160 on `sphere`, 888
//! of 1128 on `torus`, all 512 on `thin_plate` and 4101 of 6176 on
//! `noise_cavity` carry a patch that is planar by construction."* Those four
//! numbers are, to the cell, **the count of cells with at most four crossings** —
//! and four crossings are coplanar only on a field that is affine across the
//! cell, or by an exact symmetry. The census's own justification (*"Three points
//! are coplanar"*) licenses the three-crossing cells and no more.
//!
//! This harness records the correction instead of asserting it, in four columns
//! per row: `three_crossing_cells`, `three_crossing_above_floor` (asserted
//! **0** — the coplanarity theorem, on every field, at every resolution),
//! `four_crossing_cells` and `four_crossing_below_floor`. The last of those is
//! the whole disagreement: on a piecewise-planar field it is the entire
//! four-crossing population, and on a curved one it is the handful of cells a
//! symmetry plane makes exactly flat. `crates/isomesh/src/**` and
//! `benches/common/**` are read-only for this row, so the module keeps its
//! header and the finding names it.
//!
//! # Arms
//!
//! One build, one grid, one shared population per row (`M-281`). The "arms" are
//! the two timed passes over that grid and the two statistics over that one
//! point set.
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `beta` pass | `common::beta::beta_per_cell` over every cell of the grid | no |
//! | extraction | `DualContouring::<f64>::new()` over the same grid, default `lambda` and `Clamp::ToCell` | no — it is C1's denominator |
//! | crossing-count cross-check | the cut count derived from public `EDGE_CORNERS`, against `HermiteCell::len()` and against `beta`'s `Some` set | **yes** — `population_mismatch` |
//! | three-crossing check | every cut-3 cell must land below the floor | **yes** — `three_crossing_above_floor` |
//! | floor sensitivity | `rank_correlation` floored, and unfloored beside it | **yes** — `rank_correlation_unfloored` |
//! | singular threshold strip | `t = LAMBDA` as registered, `t = 1e-9` beside it | **yes** — a count identical at both is decided by neither |
//!
//! `resolution` counts **samples**, so `n` samples span `n - 1` cells
//! (`benches/common/mod.rs:40-43`), and `cells = (n - 1)^3` is that arithmetic
//! and nothing else. `surface_cells` is the paired population.
//!
//! **Two resolutions, 33 and 65, and the registration names none.** The
//! authoring contract's default ladder is `33^3` and `65^3`; `129^3` is left out
//! because `beta.rs:101-118` already publishes the committed `129^3`
//! decomposition of the `beta` pass on this same machine — 18.9 ms of prepass
//! inside 60.6 ms of pass, against 35.4 ms for Marching Cubes — and C1's verdict
//! turns on an order of magnitude rather than on a factor of two, so a third
//! resolution would buy scatter and no verdict.
//!
//! # Which reading carries each clause
//!
//! - **C1 is per-row.** `beta_share = beta_ms / extract_ms` on that row, against
//!   the registered bar `0.10`.
//! - **C2 is global**, because the clause counts fields: *at least six of
//!   eight*. The same boolean is written on every row, decided at
//!   `resolution = 65` — the finer of the two grids. `c2_fields_above_bar`
//!   records the count that decided it, and `rank_correlation` is per-row
//!   throughout, so the file settles the other resolution too.
//! - **C3 is per-row**: `confidently_wrong_cells > 0` on that row.
//!
//! # C3's thresholds, stated because they are percentiles
//!
//! C3 asks for *"cells with low residual and high `beta`"* — the cells where the
//! plane fit is confidently wrong. Both sides are quantiles of that row's own
//! surface-cell population, so the thresholds are per-field and per-resolution
//! and both are written into the CSV (`low_residual_threshold`,
//! `high_beta_threshold`):
//!
//! ```text
//! low residual  :  residual <= p25(residual)
//! high beta     :  beta*diam >  p75(beta*diam)      // floored column
//! ```
//!
//! The asymmetry between `<=` and `>` is forced by the mechanism rather than
//! chosen. On the `beta` side the floored column has a tie block at exactly `0`
//! that runs from a third of the population to all of it, so on a field where it
//! exceeds three quarters `p75` **is** `0.0` and `>` then reads "departs from
//! planarity at all", which is the only available sense of "high". On the
//! residual side the same choice inverted: `<=` keeps the exactly-zero
//! residuals, which are the cells where the fit is *perfect*, and those are
//! precisely the ones C3 is about. Using `<` on both sides would let a row whose
//! `p25(residual)` is `0.0` report an empty C3 set for a reason about the
//! threshold convention rather than about the geometry — a zero that could not
//! have been non-zero, which is `M-44`'s failure.
//!
//! Quantiles are nearest-rank on a `total_cmp` sort:
//! `index = floor(q * (len - 1))`. No interpolation, so a threshold is always a
//! value some cell actually took.
//!
//! # `singular_cells`, from Sylvester rather than from an eigensolver
//!
//! `M = sum n_i n_i^T` is assembled here from the same crossings.
//! "Rank deficient" is read at the crate's own scale — **the direction the
//! regulariser decides rather than the data**:
//!
//! ```text
//! singular(cell)  <=>  lambda_min(M) <= LAMBDA = 0.01
//! ```
//!
//! which is exactly when `lambda*I` dominates `M` along some direction of
//! `M + lambda*I`. The two constants the corpus circulates agree here to the
//! digit: DC's own SVD truncation is `sigma = 0.1` on `A` (`solve.rs:47-49`,
//! `:56-57`), and `sigma^2 = 0.01 = LAMBDA` is the corresponding eigenvalue of
//! `M = A^T A`. The bar is not invented by this harness.
//!
//! The test is **Sylvester's criterion on `M - t*I`**, not an eigenvalue:
//! `lambda_min(M) > t` iff `M - t*I` is positive definite iff its three leading
//! principal minors are all strictly positive. That is exact for a symmetric
//! matrix, costs about fifteen flops, and sidesteps the failure
//! `beta.rs:766-773` names — the closed-form cubic loses most of its digits on
//! precisely the nearly-planar cell whose smallest eigenvalue is being measured.
//! `singular_cells_tight` runs the same test at `t = 1e-9`, so the CSV carries a
//! strip rather than one threshold.
//!
//! No eigensolver is written here and none is borrowed: `common::beta`'s
//! `jacobi_eigen` is private to that module, and `common::metric`'s `Sym3::eigen`
//! belongs to R-146, whose own header explains why the two concurrently authored
//! tickets deliberately do not depend on one another (`beta.rs:128-134`).
//!
//! # `spearman_ceiling`: what C2's bar can reach given the tie structure
//!
//! `beta`'s tie block is not a nuisance, it is a cap on Spearman. Take `a` to be
//! the fraction of the paired population at `beta = 0` after flooring. Against an
//! **untied** companion the largest attainable coefficient is the concordant
//! ordering's, and it is closed-form: writing `r(u) = a/2` for `u < a` and
//! `r(u) = u` above, both `Var(r)` and `Cov(r, u)` come to `(1 - a^3)/12` while
//! `Var(u) = 1/12`, so the ratio is `sqrt(1 - a^3)`. Inverting at C2's bar,
//!
//! ```text
//! sqrt(1 - a^3) = 0.7   <=>   a = (1 - 0.49)^(1/3) = 0.798936...
//! ```
//!
//! — so against an untied residual, a planar fraction above 79.89% puts C2's
//! `0.7` out of reach before the run starts.
//!
//! **That closed form is not what is recorded, because the residual is not
//! guaranteed untied.** A residual carrying its own tie block in the same cells
//! can exceed `sqrt(1 - a^3)`, and a piecewise-planar field is exactly the case
//! where it does: every flat cell fits perfectly, so the residual has a tie
//! block too and the two block structures agree. So the recorded
//! `spearman_ceiling` is the **exact** bound over both observed tie structures,
//! obtained without a new mechanism: sort the two columns independently and ask
//! `common::beta::rank_correlation` for the coefficient of that pairing. By the
//! rearrangement inequality the concordant pairing maximises the covariance of
//! two fixed multisets while leaving both variances alone, so that coefficient is
//! the largest Spearman any pairing of these two columns could produce.
//! `spearman_ceiling_admits_bar` is the boolean, and a field whose measured
//! coefficient falls short of `0.7` while its ceiling also falls short was
//! refused by arithmetic rather than by measurement — which is the distinction
//! the entry has to draw.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"C1 moves the vertex-placement stage, whose cost
//! `M-25` puts at 3% over Surface Nets."* Discharged, and it discharges against
//! C1 rather than with it: the stage the pass would have to be folded into is
//! worth **3%** of an extraction, and `beta.rs:101-118`'s committed measurement
//! puts the whole `beta` pass at **1.3x-1.7x** an extraction. A pass costing more
//! than an extraction cannot be folded into 3% of one. There is no share to move,
//! and C1's bar of `0.10` is more than an order of magnitude away from the
//! mechanism's own published cost.
//!
//! Which is why the count columns matter more here than the milliseconds.
//! `slab_fits` is the number of thinnest-slab searches performed and equals
//! `surface_cells` exactly; `slab_evaluations = slab_fits * 73` is the objective
//! count. Both are machine-independent integers, and unlike a ratio they cannot
//! be moved by this host's `amd-pstate-epp` governor, which `M-280` measured
//! swinging the same binary 1.45x between runs. `beta_ms` and `extract_ms` are
//! medians of five timed repeats after one warm-up, with min and max beside them,
//! and `beta_share_min`/`beta_share_max` bracket the share over the extreme
//! repeat pairing — so a row whose repeats disagree by more than C1's own bar is
//! visible as scatter rather than averaged into a verdict.
//!
//! # Vacuity controls
//!
//! `M-44`: a zero that could not have been non-zero is not a measurement. Each of
//! these runs before the first `run.record` and panics with a message starting
//! `VOID: `.
//!
//! - **The registered one.** `qef_residual_variance > 0` on **every** field and
//!   resolution. Recorded as a column, because the registration's own falsifier
//!   says the correlation is against a constant otherwise.
//! - **`population_mismatch == 0` on every row** — the count of cells where
//!   `beta`'s `Some` set and this file's `cut >= 3` set disagree. Non-zero there
//!   means the two statistics are of two different point sets and the comparison
//!   is void, not merely wrong.
//! - **`three_crossing_above_floor == 0` on every row** — three points are
//!   coplanar, so a cut-3 cell above `PLANAR_FLOOR` means the floor is below the
//!   noise it was chosen to sit above, and every planar/non-planar split in the
//!   row is then an artefact.
//! - **The dead band is at least three decades wide on every row** that has cells
//!   on both sides of the floor: `beta_min_above_floor > 1e3 *
//!   beta_max_below_floor`. This is what licenses the floor as a measurement
//!   rather than a tolerance.
//! - **`extract_ms > 0` on every row** — C1's share is a division by it, and a
//!   zero denominator would report `inf` as a pass.
//! - **`triangles > 0` on every row** — an extraction that produced no mesh is
//!   not the extraction cost of anything.
//! - **`nonplanar_cells > 0` somewhere in the run** (global, not per row): a
//!   piecewise-planar field's `beta` column is identically zero, which is a
//!   measured fact and not a fixture defect, so a per-row assertion would abort
//!   on the answer. The global form still licenses every `rank_correlation` of
//!   `0.0`: it proves `beta` has variance somewhere, so a zero coefficient on a
//!   given row is that row's geometry rather than a dead instrument.
//! - **The singular strip must be monotone** in its threshold, or
//!   `rank_deficient` is not testing what it claims.
//!
//! # Determinism
//!
//! One thread, no PRNG, no map iteration, `f64` throughout. `beta_per_cell` is a
//! pure function of the point set at a fixed 73-evaluation budget
//! (`beta.rs:49-56`), the sweep order is `z`, `y`, `x` fixed, and every sort is
//! [`f64::total_cmp`] — a total order, so a NaN would sort into view rather than
//! be dropped by a partial comparison.

#![allow(clippy::cast_precision_loss, clippy::too_many_lines)]

mod common;

use std::time::Instant;

use isomesh::dual_contouring::DualContouring;
use isomesh::dual_contouring::solve::{LAMBDA, solve};
use isomesh::fields::ReferenceField;
use isomesh::hermite::HermiteCell;
use isomesh::marching_cubes::table::{EDGE_CORNERS, is_inside};
use isomesh::{MeshBuffer, Sdf};

/// Samples per axis. Two, and the header says why `129` is not a third.
const RESOLUTIONS: [u32; 2] = [33, 65];

/// Timed repeats per arm, after one untimed warm-up. Five is the authoring
/// contract's floor for a clause with a ratio threshold.
const REPEATS: usize = 5;

/// Below this, a `beta_inf` is the rounding noise of a patch that is planar by
/// construction rather than a departure from planarity.
///
/// Read off the distribution, not chosen: measured over all sixteen
/// `(field, resolution)` pairs on this code path, every occupied decade of the
/// `beta_inf` column is at or below `1e-16` or at or above `1e-8`, and the seven
/// decades between are empty on every field at both resolutions. This sits in
/// the middle of that band. `beta_max_below_floor` and `beta_min_above_floor`
/// re-derive the separation on every run and a control asserts it.
const PLANAR_FLOOR: f64 = 1e-12;

/// Decades of separation the dead band must show for [`PLANAR_FLOOR`] to be a
/// measurement rather than a tolerance.
const FLOOR_SEPARATION: f64 = 1e3;

/// C1's bar: `beta_ms / extract_ms` must be below this.
const C1_BAR: f64 = 0.10;

/// C2's bar on the Spearman coefficient.
const C2_BAR: f64 = 0.7;

/// C2's count: at least this many of the eight fields must clear [`C2_BAR`].
const C2_MIN_FIELDS: usize = 6;

/// The resolution C2's global verdict is read at — the finer of the two.
const C2_RESOLUTION: u32 = 65;

/// C3's "low residual" quantile.
const LOW_RESIDUAL_QUANTILE: f64 = 0.25;

/// C3's "high `beta`" quantile.
const HIGH_BETA_QUANTILE: f64 = 0.75;

/// The tight arm of the `singular_cells` threshold strip: numerically rank
/// deficient, as against dominated by the shipped regulariser.
const SINGULAR_TIGHT: f64 = 1e-9;

/// `for_each_reference_field!` yields eight (`fields/mod.rs:195`).
const FIELDS: usize = 8;

/// `true` when `lambda_min(m) <= t`, by Sylvester's criterion on `m - t*I`.
///
/// `m` is the upper triangle of a symmetric `3x3` in the order
/// `[xx, xy, xz, yy, yz, zz]`. A symmetric matrix is positive definite exactly
/// when its three leading principal minors are strictly positive, so
/// `lambda_min(m) > t` iff `m - t*I` passes that test — and the negation is the
/// rank deficiency this row counts. Exact, fifteen flops, and no cubic
/// discriminant to lose digits on the flat cell that is the whole point
/// (`benches/common/beta.rs:766-773`).
fn rank_deficient(m: &[f64; 6], t: f64) -> bool {
    let (xx, xy, xz) = (m[0] - t, m[1], m[2]);
    let (yy, yz, zz) = (m[3] - t, m[4], m[5] - t);
    let minor1 = xx;
    let minor2 = xx * yy - xy * xy;
    let minor3 = xx * (yy * zz - yz * yz) - xy * (xy * zz - yz * xz) + xz * (xy * yz - yy * xz);
    !(minor1 > 0.0 && minor2 > 0.0 && minor3 > 0.0)
}

/// Nearest-rank quantile of an already-`total_cmp`-sorted slice.
///
/// `index = floor(q * (len - 1))`, no interpolation, so the value returned is one
/// some cell actually took. `0.0` for an empty slice, which cannot reach a
/// recorded row: `surface_cells >= 3` is asserted first.
fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let last = sorted.len() - 1;
    let index = (q * last as f64).floor().max(0.0) as usize;
    sorted[index.min(last)]
}

/// Arithmetic mean; `0.0` for an empty slice.
fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().sum::<f64>() / values.len() as f64
}

/// Population variance. This is the registered vacuity control's quantity.
fn variance(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let m = mean(values);
    values.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / values.len() as f64
}

/// Median, minimum and maximum of a set of timings, in that order.
fn spread(samples: &mut [f64]) -> (f64, f64, f64) {
    samples.sort_unstable_by(f64::total_cmp);
    (
        quantile(samples, 0.5),
        samples[0],
        samples[samples.len() - 1],
    )
}

/// One `(field, resolution)` measurement, before the cross-row verdicts exist.
struct Row {
    field: &'static str,
    samples: u32,
    cells: u64,
    surface_cells: usize,
    nonplanar_cells: usize,
    population_mismatch: usize,
    three_crossing_cells: usize,
    three_crossing_above_floor: usize,
    four_crossing_cells: usize,
    four_crossing_below_floor: usize,
    beta_max_below_floor: f64,
    /// `f64::INFINITY` when no cell in the row clears the floor.
    beta_min_above_floor: f64,
    /// Mean, median and max of the floored `beta_inf * diam(Q)`, in world units.
    beta_times_diam_mean: f64,
    beta_times_diam_median: f64,
    beta_times_diam_max: f64,
    /// `diam(Q) = h * sqrt(3)`, so the dimensionless statistics are these
    /// divided by it.
    diam: f64,
    residual_mean: f64,
    residual_median: f64,
    residual_max: f64,
    residual_variance: f64,
    residual_zeros: usize,
    rank_correlation: f64,
    rank_correlation_unfloored: f64,
    spearman_ceiling: f64,
    beta_ms: f64,
    beta_ms_min: f64,
    beta_ms_max: f64,
    extract_ms: f64,
    extract_ms_min: f64,
    extract_ms_max: f64,
    singular_cells: usize,
    singular_cells_tight: usize,
    low_residual_threshold: f64,
    high_beta_threshold: f64,
    confidently_wrong_cells: usize,
    vertices: usize,
    triangles: usize,
}

impl Row {
    /// C1's registered quantity: the standalone `beta` pass over an extraction.
    fn beta_share(&self) -> f64 {
        self.beta_ms / self.extract_ms
    }

    /// The share at its most favourable repeat pairing, and at its least.
    fn share_bounds(&self) -> (f64, f64) {
        (
            self.beta_ms_min / self.extract_ms_max,
            self.beta_ms_max / self.extract_ms_min,
        )
    }

    fn planar_fraction(&self) -> f64 {
        if self.surface_cells == 0 {
            return 1.0;
        }
        (self.surface_cells - self.nonplanar_cells) as f64 / self.surface_cells as f64
    }

    /// Decades between the largest sub-floor `beta` and the smallest above it.
    ///
    /// `0.0` when one side of the floor is empty, or when every sub-floor value
    /// is exactly `0.0` — in both cases the separation is not a finite number of
    /// decades and the control below tests it directly rather than through this.
    fn floor_gap_decades(&self) -> f64 {
        if self.beta_min_above_floor.is_finite() && self.beta_max_below_floor > 0.0 {
            (self.beta_min_above_floor / self.beta_max_below_floor).log10()
        } else {
            0.0
        }
    }
}

/// Measure one `(field, resolution)`.
///
/// Order matters: the `beta` pass and the extraction are each timed against a
/// warm cache, and the untimed QEF pass runs last so its own grid allocation
/// cannot land inside either timer.
fn measure<F>(name: &'static str, field: &F, samples: u32) -> Row
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, h) = common::grid::<f64, F>(field, samples);
    let n = samples as usize;
    let span = n - 1;
    let cells = (span as u64).pow(3);
    // A cube of side h has diameter h*sqrt(3), which is what Jones' definition
    // divides by and what `beta_per_cell` has already divided by
    // (`beta.rs:449-451`). Multiplying back is the same product
    // `beta_times_diam` would recompute, without a second seventy-three
    // evaluation slab search.
    let diam = h * 3.0_f64.sqrt();

    // ── arm 1: the beta pass, warmed then timed ──
    let betas = common::beta::beta_per_cell(field, &shape, origin, h);
    assert_eq!(
        betas.len(),
        cells as usize,
        "{name} at {samples}: beta_per_cell returned {} entries for {cells} cells",
        betas.len()
    );
    let mut beta_samples = Vec::with_capacity(REPEATS);
    for _ in 0..REPEATS {
        let start = Instant::now();
        let again = common::beta::beta_per_cell(field, &shape, origin, h);
        beta_samples.push(start.elapsed().as_secs_f64() * 1e3);
        std::hint::black_box(&again);
    }
    let (beta_ms, beta_ms_min, beta_ms_max) = spread(&mut beta_samples);

    // ── arm 2: the extraction, C1's denominator ──
    let mut mesher = DualContouring::<f64>::new();
    let mut mesh = MeshBuffer::<f64>::new();
    mesher
        .extract(field, &shape, origin, h, &mut mesh)
        .expect("dual contouring extracts on a benchmark grid");
    let (vertices, triangles) = (mesh.vertex_count(), mesh.triangle_count());
    let mut extract_samples = Vec::with_capacity(REPEATS);
    for _ in 0..REPEATS {
        // Outside the timer, and `reset` keeps the capacity (`mesh.rs:125`), so
        // the timed region is extraction and not allocation.
        mesh.reset();
        let start = Instant::now();
        mesher
            .extract(field, &shape, origin, h, &mut mesh)
            .expect("dual contouring extracts on a benchmark grid");
        extract_samples.push(start.elapsed().as_secs_f64() * 1e3);
    }
    let (extract_ms, extract_ms_min, extract_ms_max) = spread(&mut extract_samples);

    // ── the untimed QEF pass, over beta's own population ──
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

    let mut beta_times_diam: Vec<f64> = Vec::new();
    let mut beta_times_diam_raw: Vec<f64> = Vec::new();
    let mut residual: Vec<f64> = Vec::new();
    let mut nonplanar_cells = 0usize;
    let mut population_mismatch = 0usize;
    let mut residual_zeros = 0usize;
    let mut singular_cells = 0usize;
    let mut singular_cells_tight = 0usize;
    let mut three_crossing_cells = 0usize;
    let mut three_crossing_above_floor = 0usize;
    let mut four_crossing_cells = 0usize;
    let mut four_crossing_below_floor = 0usize;
    let mut beta_max_below_floor = 0.0f64;
    let mut beta_min_above_floor = f64::INFINITY;

    for cz in 0..span {
        for cy in 0..span {
            for cx in 0..span {
                let cell = cx + span * (cy + span * cz);
                // Corner `c` sits at `[c & 1, (c >> 1) & 1, (c >> 2) & 1]`.
                // `cube::corner_offset` is private and this is the three-line
                // replacement the API notes prescribe; `HermiteCell::len()` is
                // asserted against `cut` below, which is the check that catches
                // a permuted guess.
                let mut corners = [0.0f64; 8];
                for (c, slot) in corners.iter_mut().enumerate() {
                    let x = cx + (c & 1);
                    let y = cy + ((c >> 1) & 1);
                    let z = cz + ((c >> 2) & 1);
                    *slot = values[x + n * (y + n * z)];
                }

                let mut cut = 0usize;
                for [lo, hi] in EDGE_CORNERS {
                    if is_inside(corners[lo as usize]) != is_inside(corners[hi as usize]) {
                        cut += 1;
                    }
                }
                let entry = betas[cell];
                if (cut >= 3) != entry.is_some() {
                    population_mismatch += 1;
                }
                let Some(beta) = entry else {
                    continue;
                };

                let planar = beta <= PLANAR_FLOOR;
                if planar {
                    beta_max_below_floor = beta_max_below_floor.max(beta);
                } else {
                    beta_min_above_floor = beta_min_above_floor.min(beta);
                    nonplanar_cells += 1;
                }
                if cut == 3 {
                    three_crossing_cells += 1;
                    if !planar {
                        three_crossing_above_floor += 1;
                    }
                }
                if cut == 4 {
                    four_crossing_cells += 1;
                    if planar {
                        four_crossing_below_floor += 1;
                    }
                }

                let cell_origin = [
                    origin[0] + h * cx as f64,
                    origin[1] + h * cy as f64,
                    origin[2] + h * cz as f64,
                ];
                let hermite = HermiteCell::from_corners(field, &corners, cell_origin, h);
                assert_eq!(
                    hermite.len(),
                    cut,
                    "{name} at {samples}: HermiteCell cut {} edges where the corner signs cut \
                     {cut} — the corner order handed to `from_corners` is not the crate's",
                    hermite.len()
                );
                // `solve` returns `None` only for a cell with no crossings
                // (`solve.rs:237-238`), and this one has at least three.
                let x = solve(&hermite).expect("a three-crossing cell is not empty");

                let mut sum = 0.0;
                let mut m = [0.0f64; 6];
                for crossing in hermite.iter() {
                    let (nrm, p) = (crossing.normal, crossing.position);
                    let d =
                        nrm[0] * (x[0] - p[0]) + nrm[1] * (x[1] - p[1]) + nrm[2] * (x[2] - p[2]);
                    sum += d * d;
                    m[0] += nrm[0] * nrm[0];
                    m[1] += nrm[0] * nrm[1];
                    m[2] += nrm[0] * nrm[2];
                    m[3] += nrm[1] * nrm[1];
                    m[4] += nrm[1] * nrm[2];
                    m[5] += nrm[2] * nrm[2];
                }

                if rank_deficient(&m, LAMBDA) {
                    singular_cells += 1;
                }
                if rank_deficient(&m, SINGULAR_TIGHT) {
                    singular_cells_tight += 1;
                }
                // A sum of squares is non-negative, so `<= 0.0` is `== 0.0`
                // without asking clippy's `float_cmp` for an exemption.
                if sum <= 0.0 {
                    residual_zeros += 1;
                }
                beta_times_diam.push(if planar { 0.0 } else { beta * diam });
                beta_times_diam_raw.push(beta * diam);
                residual.push(sum);
            }
        }
    }

    let surface_cells = beta_times_diam.len();
    let rank_correlation = common::beta::rank_correlation(&beta_times_diam, &residual);
    let rank_correlation_unfloored =
        common::beta::rank_correlation(&beta_times_diam_raw, &residual);

    let mut beta_sorted = beta_times_diam.clone();
    beta_sorted.sort_unstable_by(f64::total_cmp);
    let mut residual_sorted = residual.clone();
    residual_sorted.sort_unstable_by(f64::total_cmp);

    // The exact ceiling over both observed tie structures. Both columns sorted
    // ascending is the concordant pairing, which maximises the covariance of two
    // fixed multisets (rearrangement inequality) and leaves both variances
    // untouched — so this is the largest Spearman any pairing could produce.
    let spearman_ceiling = common::beta::rank_correlation(&beta_sorted, &residual_sorted);

    let low_residual_threshold = quantile(&residual_sorted, LOW_RESIDUAL_QUANTILE);
    let high_beta_threshold = quantile(&beta_sorted, HIGH_BETA_QUANTILE);
    let confidently_wrong_cells = residual
        .iter()
        .zip(&beta_times_diam)
        .filter(|(r, b)| **r <= low_residual_threshold && **b > high_beta_threshold)
        .count();

    Row {
        field: name,
        samples,
        cells,
        surface_cells,
        nonplanar_cells,
        population_mismatch,
        three_crossing_cells,
        three_crossing_above_floor,
        four_crossing_cells,
        four_crossing_below_floor,
        beta_max_below_floor,
        beta_min_above_floor,
        beta_times_diam_mean: mean(&beta_times_diam),
        beta_times_diam_median: quantile(&beta_sorted, 0.5),
        beta_times_diam_max: quantile(&beta_sorted, 1.0),
        diam,
        residual_mean: mean(&residual),
        residual_median: quantile(&residual_sorted, 0.5),
        residual_max: quantile(&residual_sorted, 1.0),
        residual_variance: variance(&residual),
        residual_zeros,
        rank_correlation,
        rank_correlation_unfloored,
        spearman_ceiling,
        beta_ms,
        beta_ms_min,
        beta_ms_max,
        extract_ms,
        extract_ms_min,
        extract_ms_max,
        singular_cells,
        singular_cells_tight,
        low_residual_threshold,
        high_beta_threshold,
        confidently_wrong_cells,
        vertices,
        triangles,
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-152");

    common::experiment::run(prereg, |run| {
        let mut rows: Vec<Row> = Vec::with_capacity(FIELDS * RESOLUTIONS.len());
        for samples in RESOLUTIONS {
            // Inline block per field, not a closure, so no `return` in here
            // (M-253).
            isomesh::for_each_reference_field!(f64, |name, field| {
                let row = measure(name, &field, samples);
                println!(
                    "  {:<14} {:>4}^3  surface {:>6}  planar {:>7.4}  cut3 {:>6}  \
                     cut4-flat {:>6}/{:<6}  rho {:>7.4}  ceiling {:>7.4}  \
                     beta {:>9.3} ms  dc {:>9.3} ms  share {:>7.3}  singular {:>6}  c3 {:>6}",
                    row.field,
                    row.samples,
                    row.surface_cells,
                    row.planar_fraction(),
                    row.three_crossing_cells,
                    row.four_crossing_below_floor,
                    row.four_crossing_cells,
                    row.rank_correlation,
                    row.spearman_ceiling,
                    row.beta_ms,
                    row.extract_ms,
                    row.beta_share(),
                    row.singular_cells,
                    row.confidently_wrong_cells,
                );
                rows.push(row);
            });
        }

        // ── vacuity controls, all of them, before the first record ──

        assert_eq!(
            rows.len(),
            FIELDS * RESOLUTIONS.len(),
            "VOID: expected {} rows ({FIELDS} fields x {} resolutions) and built {}, so some \
             field or resolution was never measured and its clauses were never at risk",
            FIELDS * RESOLUTIONS.len(),
            RESOLUTIONS.len(),
            rows.len()
        );

        for row in &rows {
            assert!(
                row.surface_cells >= 3,
                "VOID: {} at {}^3 has {} paired cells, and `common::beta` defines a rank \
                 correlation over fewer than three pairs to be 0.0 rather than measuring one",
                row.field,
                row.samples,
                row.surface_cells
            );
            // The registration's own vacuity control, verbatim: "the residual
            // column must have non-zero variance on every field, or the
            // correlation is against a constant".
            assert!(
                row.residual_variance > 0.0,
                "VOID: {} at {}^3 has a QEF residual of zero variance over {} cells, so C2's \
                 rank correlation is taken against a constant and its value is an accident of \
                 tie-breaking rather than a measurement",
                row.field,
                row.samples,
                row.surface_cells
            );
            assert_eq!(
                row.population_mismatch, 0,
                "VOID: {} at {}^3 disagrees with `common::beta` about {} cells' surface \
                 membership, so `beta` and the QEF residual are statistics of two different \
                 point sets and the whole comparison is void",
                row.field, row.samples, row.population_mismatch
            );
            assert_eq!(
                row.three_crossing_above_floor, 0,
                "VOID: {} at {}^3 puts {} of its {} three-crossing cells above PLANAR_FLOOR = \
                 {PLANAR_FLOOR:e}, and three points are coplanar — so the floor is inside the \
                 noise it was chosen to sit above and every planar split in this row is an \
                 artefact",
                row.field, row.samples, row.three_crossing_above_floor, row.three_crossing_cells
            );
            assert!(
                !row.beta_min_above_floor.is_finite()
                    || row.beta_min_above_floor > FLOOR_SEPARATION * row.beta_max_below_floor,
                "VOID: {} at {}^3 has its largest sub-floor beta at {:e} and its smallest \
                 above-floor beta at {:e}, less than {FLOOR_SEPARATION:e} apart — the floor is \
                 cutting through the distribution rather than through a dead band, so \
                 planar_fraction and rank_correlation are threshold artefacts",
                row.field,
                row.samples,
                row.beta_max_below_floor,
                row.beta_min_above_floor
            );
            assert!(
                row.extract_ms > 0.0,
                "VOID: {} at {}^3 measured an extraction time of {} ms, and C1's share is a \
                 division by it — a zero denominator would report `inf` as a verdict",
                row.field,
                row.samples,
                row.extract_ms
            );
            assert!(
                row.triangles > 0,
                "VOID: {} at {}^3 extracted {} triangles, so C1's denominator is the cost of \
                 producing no mesh",
                row.field,
                row.samples,
                row.triangles
            );
        }

        // Global rather than per row, and the header says why: a piecewise-planar
        // field's beta column is identically zero, which is the answer and not a
        // fixture defect. What has to be shown is that the instrument can read a
        // non-planar patch at all.
        let nonplanar_total: usize = rows.iter().map(|r| r.nonplanar_cells).sum();
        assert!(
            nonplanar_total > 0,
            "VOID: not one cell in the whole run clears PLANAR_FLOOR, so every rank correlation \
             is against a constant column and every zero is a zero that could not have been \
             non-zero (M-44)"
        );

        // The threshold strip has to move, or `singular_cells` is reporting a
        // floor rather than a rank.
        let singular_at_lambda: usize = rows.iter().map(|r| r.singular_cells).sum();
        let singular_at_tight: usize = rows.iter().map(|r| r.singular_cells_tight).sum();
        assert!(
            singular_at_lambda >= singular_at_tight,
            "VOID: the singular count at LAMBDA ({singular_at_lambda}) is below the count at \
             {SINGULAR_TIGHT:e} ({singular_at_tight}), which is impossible for a monotone \
             threshold and means `rank_deficient` is not testing what it claims"
        );

        // ── C2's global verdict ──
        let c2_fields_above_bar = rows
            .iter()
            .filter(|r| r.samples == C2_RESOLUTION && r.rank_correlation > C2_BAR)
            .count();
        let c2_holds = c2_fields_above_bar >= C2_MIN_FIELDS;
        println!(
            "\n  C2 at {C2_RESOLUTION}^3: {c2_fields_above_bar} of {FIELDS} fields above \
             {C2_BAR} — needs {C2_MIN_FIELDS}, so C2 {}",
            if c2_holds { "HOLDS" } else { "is FALSIFIED" }
        );
        for row in rows.iter().filter(|r| r.samples == C2_RESOLUTION) {
            if row.spearman_ceiling <= C2_BAR {
                println!(
                    "    {:<14} planar {:.4} puts the exact concordance ceiling at {:.4}, so \
                     {C2_BAR} was unreachable before the run",
                    row.field,
                    row.planar_fraction(),
                    row.spearman_ceiling
                );
            }
        }

        for row in &rows {
            let share = row.beta_share();
            let (share_min, share_max) = row.share_bounds();
            let c1_holds = share < C1_BAR;
            let c3_holds = row.confidently_wrong_cells > 0;
            let slab_fits = row.surface_cells;
            let slab_evaluations = slab_fits * common::beta::SLAB_EVALUATIONS;
            let min_above = if row.beta_min_above_floor.is_finite() {
                row.beta_min_above_floor
            } else {
                0.0
            };

            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.samples.to_string()),
                (
                    "beta_infinity",
                    format!("{:.9}", row.beta_times_diam_mean / row.diam),
                ),
                (
                    "beta_times_diam",
                    format!("{:.9}", row.beta_times_diam_mean),
                ),
                ("qef_residual", format!("{:.9e}", row.residual_mean)),
                ("rank_correlation", format!("{:.6}", row.rank_correlation)),
                ("cells", row.cells.to_string()),
                ("beta_ms", format!("{:.4}", row.beta_ms)),
                ("beta_share", format!("{share:.6}")),
                ("singular_cells", row.singular_cells.to_string()),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_holds.to_string()),
                ("c3_holds", c3_holds.to_string()),
                // ── extras (M-273) ──
                (
                    "beta_infinity_max",
                    format!("{:.9}", row.beta_times_diam_max / row.diam),
                ),
                (
                    "beta_infinity_median",
                    format!("{:.9}", row.beta_times_diam_median / row.diam),
                ),
                ("beta_floor", format!("{PLANAR_FLOOR:e}")),
                (
                    "beta_floor_gap_decades",
                    format!("{:.3}", row.floor_gap_decades()),
                ),
                (
                    "beta_max_below_floor",
                    format!("{:.6e}", row.beta_max_below_floor),
                ),
                ("beta_min_above_floor", format!("{min_above:.6e}")),
                ("beta_ms_max", format!("{:.4}", row.beta_ms_max)),
                ("beta_ms_min", format!("{:.4}", row.beta_ms_min)),
                (
                    "beta_ns_per_surface_cell",
                    format!("{:.1}", row.beta_ms * 1e6 / slab_fits.max(1) as f64),
                ),
                ("beta_share_max", format!("{share_max:.6}")),
                ("beta_share_min", format!("{share_min:.6}")),
                (
                    "beta_times_diam_max",
                    format!("{:.9}", row.beta_times_diam_max),
                ),
                (
                    "beta_times_diam_median",
                    format!("{:.9}", row.beta_times_diam_median),
                ),
                ("c1_bar", format!("{C1_BAR:.2}")),
                ("c2_bar", format!("{C2_BAR:.2}")),
                ("c2_fields_above_bar", c2_fields_above_bar.to_string()),
                ("c2_resolution", C2_RESOLUTION.to_string()),
                (
                    "confidently_wrong_cells",
                    row.confidently_wrong_cells.to_string(),
                ),
                ("diam", format!("{:.9}", row.diam)),
                ("extract_ms", format!("{:.4}", row.extract_ms)),
                ("extract_ms_max", format!("{:.4}", row.extract_ms_max)),
                ("extract_ms_min", format!("{:.4}", row.extract_ms_min)),
                (
                    "four_crossing_below_floor",
                    row.four_crossing_below_floor.to_string(),
                ),
                ("four_crossing_cells", row.four_crossing_cells.to_string()),
                (
                    "high_beta_threshold",
                    format!("{:.9}", row.high_beta_threshold),
                ),
                ("lambda", format!("{LAMBDA:.4}")),
                (
                    "low_residual_threshold",
                    format!("{:.9e}", row.low_residual_threshold),
                ),
                ("nonplanar_cells", row.nonplanar_cells.to_string()),
                ("planar_fraction", format!("{:.6}", row.planar_fraction())),
                ("population_mismatch", row.population_mismatch.to_string()),
                ("qef_residual_max", format!("{:.9e}", row.residual_max)),
                (
                    "qef_residual_median",
                    format!("{:.9e}", row.residual_median),
                ),
                (
                    "qef_residual_variance",
                    format!("{:.9e}", row.residual_variance),
                ),
                ("qef_residual_zeros", row.residual_zeros.to_string()),
                (
                    "rank_correlation_unfloored",
                    format!("{:.6}", row.rank_correlation_unfloored),
                ),
                ("repeats", REPEATS.to_string()),
                ("singular_bar", format!("{LAMBDA:.4}")),
                ("singular_cells_tight", row.singular_cells_tight.to_string()),
                ("slab_evaluations", slab_evaluations.to_string()),
                (
                    "slab_evaluations_per_patch",
                    common::beta::SLAB_EVALUATIONS.to_string(),
                ),
                ("slab_fits", slab_fits.to_string()),
                ("spearman_ceiling", format!("{:.6}", row.spearman_ceiling)),
                (
                    "spearman_ceiling_admits_bar",
                    (row.spearman_ceiling > C2_BAR).to_string(),
                ),
                ("surface_cells", row.surface_cells.to_string()),
                (
                    "three_crossing_above_floor",
                    row.three_crossing_above_floor.to_string(),
                ),
                ("three_crossing_cells", row.three_crossing_cells.to_string()),
                ("triangles", row.triangles.to_string()),
                ("vertices", row.vertices.to_string()),
            ]);
        }
    });
}
