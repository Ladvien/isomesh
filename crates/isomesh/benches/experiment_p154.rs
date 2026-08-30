//! **P-154 — a triangle budget computed before meshing, from a convergent sum.**
//!
//! Ticket: R-154. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p154
//! ```
//!
//! Writes `docs/experiments/p-154.csv`.
//!
//! # What was missing
//!
//! `MeshSink` has a capacity hook — `fn reserve(&mut self, vertices: usize,
//! triangles: usize)` (`crates/isomesh/src/mesh.rs:225-232`) — and **no
//! extractor in the crate ever calls it.** A `grep` for `reserve(` over
//! `src/marching_cubes/`, `src/dual.rs` and `src/greedy_quads.rs` finds
//! `GreedyQuads`' `solid.reserve(cell_count)` (`greedy_quads.rs:181`), which is
//! a *cell* count and not a triangle count, and `DualMesher`'s
//! `slot_vertex.reserve(slot_position.len())` (`dual.rs:638`), which reserves
//! against work already done. The mesh sink is never sized in advance, because
//! nothing in the crate can say how big the mesh will be.
//!
//! The one budget the crate does own is denominated in **time**:
//! `DirtySet::mesh_within_budget` (`chunk/dirty.rs:121`) takes a caller-supplied
//! `spend: FnMut() -> bool` predicate and meshes nearest-first until it returns
//! false (`dirty.rs:98-103`: *"a `no_std` crate cannot read a clock"*). It
//! answers *"stop when the frame is gone"*, never *"this chunk will cost N
//! triangles"*.
//!
//! What exists instead is `M-13`, and it is explicitly **after** the fact
//! (`FINDINGS.md:1129`):
//!
//! > **Surface cells ≈ `1.5·A/h²`, not `A/h²`.** Measured `1.450` (25³), `1.442`
//! > (33³), `1.517` (64³) on the unit sphere. The constant is derivable: a plane
//! > of unit normal `n` crosses `(|nx|+|ny|+|nz|)/h²` cells per unit area, and
//! > `E[|nx|] = ½` over the sphere, so an isotropic surface gives
//! > `E[Σ|nᵢ|] = 3/2`. Predicted 6,430 triangles at 64³ from `A/h²` and measured
//! > **9,452** — a 1.47× miss, which is this factor.
//!
//! `A` is the mesh's own area, so `M-13` prices a mesh that already exists. Its
//! own distilled lesson says as much (`FINDINGS.md:1815`): *"Estimate a count
//! from the geometry, then measure it before writing it down."*
//!
//! Azzam & Schul's higher-dimensional Traveling Salesman Theorem
//! (`arXiv:1609.02892`, acquired for `R-152`) offers the other direction. Its
//! content is that
//!
//! ```text
//!     S = sum over dyadic Q of beta_inf(Q)^2 * diam(Q)^d
//! ```
//!
//! is finite precisely for `d`-rectifiable sets, with `d = 2` here
//! (`benches/common/beta.rs:22-30`, `:149-153`). `S` has the units of an area
//! and is computed **from the field alone**, so `S / h^2` is a count and needs
//! no mesh. `common::beta::beta_sum` (R-152's module, consumed unchanged,
//! `beta.rs:293-410`) accumulates it and reports `per_scale` increments,
//! `converges(rel_tol)` and `scales_used()` (`beta.rs:310-336`).
//!
//! # What P-152 measured, quoted from its committed CSV
//!
//! `docs/experiments/p-152.csv`, 16 rows × 61 columns, commit `98c5309`. All
//! three of its clauses were falsified and two of those results are load-bearing
//! here:
//!
//! - **`beta` costs more than the extraction it describes.** `beta_share` runs
//!   `0.824`–`1.694` against C1's bar of `0.10` — `gyroid` at `65^3` is
//!   `0.964901`, `thin_plate` at `33^3` is `1.693981`. The `beta` pass is not a
//!   fraction of a mesh; it is a mesh's worth of work.
//! - **`beta` is exactly zero on the piecewise-planar fields.** `beta_infinity`
//!   is `0.000000000` for `box_exact` and `thin_plate` at both resolutions, with
//!   `planar_fraction = 1.000000` and `nonplanar_cells = 0`. Per *cell* — over
//!   the cell's own twelve edge crossings — those two fields carry no departure
//!   from planarity at all.
//!
//! That second number is the reason this row cuts the box at **dyadic sub-boxes
//! of the whole domain** rather than at extraction cells. A cell-sized patch of
//! `box_exact` is flat; a domain-quarter-sized sub-box of `box_exact` contains an
//! edge or a corner, and no plane fits that. `beta_sum` is therefore non-zero on
//! all eight fields, which is what lets the prediction be *wrong* rather than
//! *undefined* — the difference between a falsification and a vacuity.
//!
//! P-152's per-cell columns are not re-copied into this CSV. They are committed,
//! and one answer belongs in one file.
//!
//! One column here *is* comparable to P-152's, and it is used as a cross-check
//! rather than a copy: `surface_cells`. P-152 counts it through
//! `HermiteCell`/`beta_per_cell`; this file counts it from its own grid pass with
//! the crate's `is_inside`. Measured out of tree, the two agree **exactly on all
//! sixteen shared `(field, resolution)` pairs** — `1160 / 1128 / 1352 / 1388 /
//! 512 / 5240 / 1958 / 6176` at `33^3` and `4760 / 4208 / 5768 / 6014 / 2048 /
//! 21432 / 8413 / 28375` at `65^3`. Two independent implementations of "which
//! cells hold surface" landing on the same integer sixteen times is what licenses
//! `M-13`'s cells-level column below.
//!
//! # The prediction, and the one parameter it spends
//!
//! `S` is an area and a triangle count is an area over `h^2`, so the only
//! dimensionally admissible one-parameter form is
//!
//! ```text
//!     predicted_triangles = k * beta_sum / h^2
//! ```
//!
//! **`k` is fitted, and this is therefore a one-parameter prediction and not a
//! parameter-free one.** It is fitted once, globally, over all
//! `8 x 3 = 24` rows — never per field and never per resolution — as the
//! geometric mean of `actual / (beta_sum / h^2)`. That is the exact minimiser of
//! the sum of squared *log* errors, which is the right loss for a clause phrased
//! as a relative bound: "within 25%" is a statement about a ratio, so the fit
//! that C1 should be judged against is the one that centres the ratio.
//! `dof_spent = 1`, `fit_rows = 24`, `dof_residual = 23`, and all three are
//! columns rather than prose.
//!
//! `k` is **not** chosen to maximise the number of rows inside C1's bar. Tuning
//! the estimator to the verdict is the failure this apparatus exists to prevent,
//! so the honest fit is recorded — and then the question "would a luckier
//! constant have passed?" is answered exactly rather than dodged:
//!
//! ```text
//!     row i is inside the bar  <=>  k in [(1 - 0.25) * r_i, (1 + 0.25) * r_i]
//!     where r_i = actual_i / (beta_sum_i / h_i^2)
//! ```
//!
//! A field passes C1 when one `k` puts **all three** of its resolutions inside
//! the bar, i.e. when the intersection of its three intervals is non-empty
//! (`field_constant_lo`, `field_constant_hi`, `field_feasible_any_constant`).
//! `max_fields_any_constant` is the largest number of fields any single constant
//! can satisfy at once, found by evaluating the count at every interval endpoint
//! — the optimum of a closed-interval stabbing problem is attained there. If
//! that number is below C1's five, C1 is refused by arithmetic and not by the
//! fit, and the entry has to say so.
//!
//! # The incumbent, priced three ways
//!
//! C1's falsifier says a miss *"would leave `M-13`'s area law the better
//! predictor"*, so `M-13` is measured on the same 24 rows rather than cited:
//!
//! | column | form | fitted parameters |
//! |---|---|---|
//! | `m13_predicted_cells` | `1.5 * A / h^2` — `M-13` as written, against `surface_cells` | 0 |
//! | `m13_predicted_triangles` | `2.0 * 1.5 * A / h^2` — `M-13`'s own triangle arithmetic | 0 |
//! | `m13_fitted_predicted_triangles` | `k13 * A / h^2`, `k13` the geometric mean over the same 24 rows | 1 |
//! | `m13_anisotropic_predicted_triangles` | `2.0 * mean_l1_normal * A / h^2` | 0 |
//!
//! The `2.0` is `M-13`'s, not this file's: its published figure of *"6,430
//! triangles at 64³"* for the unit sphere is `2 * 4*pi / (4/64)^2 = 6434`, so
//! `M-13` converts cells to triangles at two per cell. `triangles_per_surface_cell`
//! re-measures that on every row, and the *cells* column tests `M-13`'s literal
//! claim without the conversion at all.
//!
//! The last row of that table is `M-13`'s own derivation with its own constant
//! **measured instead of assumed**. `M-13` derives `Σ|nᵢ|` and then evaluates it
//! at `3/2` for an isotropic normal distribution. An axis-aligned plane has
//! `Σ|nᵢ| = 1`, so on `box_exact` and `thin_plate` the assumed constant is 50%
//! high by construction. `mean_l1_normal` is the area-weighted mean of
//! `|nx|+|ny|+|nz|` over the mesh's triangles, computed as
//! `sum |cross|_1 / sum |cross|_2` — which needs no normalisation, because the
//! cross product's `L1` norm over its `L2` norm *is* the unit normal's `L1` norm.
//! It is still a-posteriori; it is `M-13` sharpened, not replaced.
//!
//! # `beta_sum` per unit area is the diagnostic that decides C1
//!
//! `beta_sum / mesh_area` is the conversion factor a single `k` has to be
//! constant against. It is recorded per row, and its extreme ratio over the run
//! is recorded as `beta_sum_spread` on every row. `S` is a *bending* functional —
//! `beta_inf(Q) ~ curvature * diam(Q)` on a smooth patch, so
//! `beta^2 * diam^2 ~ curvature^2 * diam^4` — while a uniform extractor's
//! triangle count is an *area* functional. Two functionals that disagree on the
//! sphere-versus-plane axis cannot be related by one constant, and the spread is
//! the size of that disagreement.
//!
//! # Arms
//!
//! One field, one grid, one `beta_sum` per row. The "arms" are the four
//! predictors above laid against one measured triangle count, plus the
//! convergence reading.
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `beta`-sum budget | `k * beta_sum / h^2`, one global `k` | no — it is C1 |
//! | `M-13` as published | `3 * A / h^2`, no fitted parameter | **yes** — the incumbent C1's falsifier names |
//! | `M-13` refitted | `k13 * A / h^2`, one global `k13` | **yes** — the same degrees of freedom as C1's |
//! | `M-13` anisotropic | `2 * mean_l1_normal * A / h^2` | **yes** — `M-13`'s constant measured, not assumed |
//! | best-case constant | the exact interval analysis over every `k` | **yes** — `max_fields_any_constant` |
//! | convergence | the dyadic increments at that row's ladder | no — it is C2 |
//!
//! # Three resolutions, and why 17 / 33 / 65
//!
//! The registration says *"at three resolutions"* and names none. The ladder is
//! forced from both ends:
//!
//! - The **vacuity control demands at least four dyadic scales**, and this row's
//!   ladder is the dyadic scales the mesh actually resolves — `2, 4, ..., n-1`
//!   samples-minus-one. `17^3` spans `16 = 2^4` cells and gives exactly four; a
//!   coarser resolution gives three and cannot observe convergence.
//! - `129^3` would extend the ladder to `n = 128`, whose sub-grid alone is
//!   `(4 * 128 + 1)^3 = 513^3 = 135` million field samples per field. The sum's
//!   own convergence is what makes that level irrelevant, and `increment_ratio_last`
//!   is the column that proves it rather than assuming it.
//!
//! So `17`, `33`, `65`, giving `scales_used` of `4`, `5`, `6`. `resolution`
//! counts **samples**, so `n` samples span `n - 1` cells
//! (`benches/common/mod.rs:40-43`).
//!
//! # Which reading carries each clause, and the two convergence columns
//!
//! **C1 is global**, because the clause counts fields: *at least five of eight*,
//! and a field passes only when all three of its resolutions do. The same
//! boolean is therefore written on every row; `c1_fields_within_bar` records the
//! count that decided it, `row_within_bar` is that row's own reading, and
//! `max_fields_any_constant` says whether any constant could have done better.
//!
//! **C2 is per row**, and the reading it takes is *divergence*, not a tolerance:
//!
//! ```text
//!     c2_holds  <=>  increment_ratio_last < 1
//! ```
//!
//! — the last dyadic increment is smaller than the one before it, so the tail is
//! decaying and the sum has somewhere to converge to. That is the theorem's own
//! dichotomy: rectifiable means the increments fall away, and the registration's
//! falsifier is *"C2 by divergence"* rather than by a slow approach.
//!
//! `sum_converges` is the **module's** instrument — `BetaSum::converges(0.10)`,
//! last increment at or under a tenth of the running total (`beta.rs:310-329`) —
//! and it is recorded beside `c2_holds` rather than instead of it, because the
//! two disagree exactly where the *truncation depth* decides rather than the
//! field. Measured out of tree, `box_exact` at `17^3` has a cleanly decaying tail
//! (`increment_ratio_last = 0.515`) whose last increment is still **12.2%** of
//! the total, so the module says `false` and the divergence test says `true`;
//! at `33^3` and `65^3` the same field reads `true` on both. A row where they
//! disagree has not diverged, it has been cut short — and `scales` and
//! `increment_last` are on the row so a reader can see which.
//!
//! **One field is expected to diverge, and the mechanism is nameable in advance.**
//! `thin_plate` is a slab of half-thickness `0.4 * (4/64) / 2 = 0.0125`
//! (`fields/mod.rs:617-627`), so its two faces are `0.025` apart while the finest
//! sub-box in this ladder is `4/64 = 0.0625` across. **Every sub-box in the
//! ladder therefore contains both faces**, and the thinnest slab containing both
//! has a half-width set by the *plate*, not by the sub-box: `w ~ 0.0125`
//! independent of `diam(Q)`, so
//!
//! ```text
//!     beta(Q)^2 * diam(Q)^2 = (w / diam)^2 * diam^2 = w^2 = constant per box
//! ```
//!
//! while the boxes carrying that patch grow with `n`. Measured out of tree, its
//! increments are `0.00031, 0.00045, 0.00075, 0.00153, 0.00283, 0.00557` — a
//! ratio approaching **2** at every level, i.e. contribution `∝ n`, i.e. the
//! truncated sum grows without bound as scales are added. Contrast `box_exact`,
//! whose `beta` also comes from a non-smooth feature but whose dihedral edge
//! gives `w ∝ diam`, so its increments *halve* (`0.515, 0.507, 0.504`).
//!
//! The plate is a smooth compact surface and is of course rectifiable. What
//! diverges is the sum **truncated above the field's own smallest length scale**,
//! and that is the honest content of the registration's phrase *"not rectifiable
//! at the resolutions we mesh at"*: convergence of the TST sum is a statement
//! about scales finer than the thinnest feature, and this crate meshes
//! `thin_plate` at a cell size 2.5× its thickness on purpose (`M-266`).
//!
//! # Cost, in integers rather than in milliseconds
//!
//! No clause here is a cost clause, so nothing is timed — this host's
//! `amd-pstate-epp` governor swings the same binary 1.45× between runs (`M-280`)
//! and a wall clock in this CSV would only invite misquotation. The cost of an
//! a-priori budget is nonetheless the reason to want one, so it is recorded as
//! machine-independent counts:
//!
//! - `field_samples` — field evaluations the whole ladder performs, summed over
//!   scales as `(4n + 1)^3` (`beta.rs:77-82`);
//! - `extract_samples` — `resolution^3`, what the extraction it predicts costs;
//! - `samples_ratio` — the first over the second;
//! - `slab_fits` — thinnest-slab searches, the sum of `per_scale`'s carried
//!   sub-box counts, at `slab_evaluations_per_patch = 73` objective evaluations
//!   each (`beta.rs:170-174`).
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"none — this predicts a cost, it does not
//! change one"*. Discharged as registered: nothing in `crates/isomesh/src/**` is
//! touched, no stage moves, and the row's whole product is a number that would
//! be *read* by a caller sizing `MeshSink::reserve`. The `samples_ratio` column
//! is what a reader needs in order to decide whether reading it is worth it, and
//! it is reported whichever way C1 goes.
//!
//! # Vacuity controls
//!
//! `M-44`: a zero that could not have been non-zero is not a measurement. Every
//! one of these runs **before** the first `run.record` and panics with a message
//! starting `VOID: `.
//!
//! - **The registered one.** `scales_used >= 4` on every row. Recorded as a
//!   column, because with fewer than four dyadic scales `converges` is reading
//!   one or two increments and convergence has not been observed.
//! - **Every scale must have carried a patch**: each `per_scale` entry's sub-box
//!   count and contribution are both strictly positive. A level that found
//!   nothing contributes a zero increment, and `increment_ratio_last` over a zero
//!   is not a convergence reading.
//! - **`beta_sum > 0` on every row.** A zero sum predicts zero triangles for
//!   every `k`, so `prediction_error_ratio` would be `1.0` for a reason about the
//!   instrument rather than about the field — and C1 would be falsified by an
//!   undefined prediction instead of by a wrong one.
//! - **`actual_triangles > 0` on every row.** A budget for an empty mesh is not a
//!   budget, and it is C1's denominator.
//! - **`mesh_area > 0` on every row** — `M-13`'s predictions divide by nothing
//!   but they are meaningless over an empty surface, and `mean_l1_normal` is a
//!   ratio of areas.
//! - **`beta_sum_spread > 2`** across the run: the conversion factor a single `k`
//!   must be constant against has to actually vary, or C1's failure would be a
//!   statement about a constant fixture rather than about the fields. This is the
//!   `M-44` control for a *negative* result — a falsification that could not have
//!   been a pass is not a falsification.
//! - **`fitted_constant` finite and strictly positive**, or the prediction column
//!   is not a prediction.
//! - **All 24 rows present**, so no field or resolution escaped its clause.
//!
//! # Determinism
//!
//! One thread, no PRNG, no clock in any column. `beta_sum` is a pure function of
//! the field at a fixed 73-evaluation-per-patch budget (`beta.rs:49-62`), the
//! ladder is a fixed list of powers of two, the field order is
//! `for_each_reference_field!`'s, and the only sort is the interval-endpoint
//! sweep, taken with [`f64::total_cmp`] on values that are all finite by the
//! controls above.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_cubes::table::is_inside;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf, Shape3};

/// Samples per axis. Three, and the header says why these three.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// `for_each_reference_field!` yields eight (`fields/mod.rs:195`).
const FIELDS: usize = 8;

/// C1's bar: `|predicted - actual| / actual` must be at or below this.
const C1_BAR: f64 = 0.25;

/// C1's count: this many of the eight fields must clear [`C1_BAR`] at **all
/// three** resolutions.
const C1_MIN_FIELDS: usize = 5;

/// The registered vacuity floor on dyadic scales.
const MIN_SCALES: usize = 4;

/// The relative tolerance [`BetaSum::converges`] is asked at.
///
/// The increments decay geometrically at a ratio the CSV records, so a last
/// increment under a tenth of the total bounds the unsummed tail at roughly the
/// same tenth — comfortably inside C1's own 25% bar, which is what makes this a
/// statement about truncation rather than a free parameter.
///
/// [`BetaSum::converges`]: common::beta::BetaSum::converges
const CONVERGENCE_TOLERANCE: f64 = 0.10;

/// `M-13`'s surface-cell constant: `surface cells ~ 1.5 * A / h^2`
/// (`FINDINGS.md:1129`).
const M13_CELLS_PER_AREA: f64 = 1.5;

/// `M-13`'s own cells-to-triangles conversion.
///
/// Not fitted here and not chosen here: `M-13` quotes *"6,430 triangles at 64³"*
/// for the unit sphere, and `2 * 4*pi / (4/64)^2 = 6434`. Two per cell is the
/// arithmetic behind its published figure. `triangles_per_surface_cell`
/// re-measures it on every row.
const M13_TRIANGLES_PER_CELL: f64 = 2.0;

/// The minimum spread of `beta_sum / mesh_area` the run must show, or C1's
/// verdict is about a constant fixture rather than about the fields.
const MIN_SPREAD: f64 = 2.0;

/// Objective evaluations per thinnest-slab search (`beta.rs:170-174`).
const SLAB_EVALUATIONS_PER_PATCH: u64 = 73;

/// The dyadic scales a mesh of `cells` cells per axis resolves: `2, 4, ...,
/// cells`.
///
/// The ladder ends at the extraction grid's own cell count, so the finest
/// sub-box in the sum is the coarsest thing the mesh cannot see inside of. The
/// `while` bound is `cells` and not `cells / 2`, so `16` cells gives exactly the
/// four levels the vacuity control requires.
fn dyadic_ladder(cells: u32) -> Vec<u32> {
    let mut out = Vec::new();
    let mut n = 2u32;
    while n <= cells {
        out.push(n);
        n *= 2;
    }
    out
}

/// Mesh area, and the area-weighted mean of `|nx| + |ny| + |nz|` over its
/// triangles.
///
/// The mean is `sum |cross|_1 / sum |cross|_2`: each triangle's cross product has
/// `L2` norm twice its area and direction its normal, so the ratio is exactly the
/// area-weighted mean of the unit normal's `L1` norm and no normalisation is
/// performed. `(0.0, 0.0)` for a mesh with no triangles, which the controls
/// refuse.
fn area_and_mean_l1(mesh: &MeshBuffer<f64>) -> (f64, f64) {
    let mut l2 = 0.0;
    let mut l1 = 0.0;
    // `as_chunks` over `chunks_exact` is the crate's own convention for walking
    // an index buffer (`common/tpms.rs:618-621`); it drops a ragged tail the
    // same way.
    for t in mesh.indices.as_chunks::<3>().0 {
        let p = mesh.positions[t[0] as usize];
        let q = mesh.positions[t[1] as usize];
        let r = mesh.positions[t[2] as usize];
        let u = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
        let v = [r[0] - p[0], r[1] - p[1], r[2] - p[2]];
        let c = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        l2 += (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        l1 += c[0].abs() + c[1].abs() + c[2].abs();
    }
    if l2 > 0.0 {
        (0.5 * l2, l1 / l2)
    } else {
        (0.0, 0.0)
    }
}

/// Cells of the extraction grid whose eight corners do not all agree in sign.
///
/// `M-13`'s claim is about this quantity and not about triangles, so it is
/// counted rather than inferred. The sign rule is the crate's own
/// [`is_inside`] — `value < 0`, exactly zero is outside (`cube.rs:157-173`) —
/// so this population is the extractor's, not a second convention.
fn surface_cells<S: Sdf<Scalar = f64>>(
    sdf: &S,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    h: f64,
) -> u64 {
    let size = shape.size();
    let mut values = Vec::with_capacity(shape.element_count());
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                values.push(sdf.sample([
                    origin[0] + f64::from(x) * h,
                    origin[1] + f64::from(y) * h,
                    origin[2] + f64::from(z) * h,
                ]));
            }
        }
    }

    let sx = size[0] as usize;
    let sy = size[1] as usize;
    let mut cells = 0u64;
    for z in 0..size[2] as usize - 1 {
        for y in 0..sy - 1 {
            for x in 0..sx - 1 {
                let base = x + y * sx + z * sx * sy;
                let mut inside = 0u8;
                for corner in 0u8..8 {
                    let i = base
                        + usize::from(corner & 1)
                        + usize::from((corner >> 1) & 1) * sx
                        + usize::from((corner >> 2) & 1) * sx * sy;
                    if is_inside(values[i]) {
                        inside += 1;
                    }
                }
                if inside != 0 && inside != 8 {
                    cells += 1;
                }
            }
        }
    }
    cells
}

/// Geometric mean of a set of strictly positive ratios.
///
/// The minimiser of the sum of squared log errors, which is the fit a relative
/// bound should be judged against. `0.0` for an empty slice, which the controls
/// refuse.
fn geometric_mean(ratios: &[f64]) -> f64 {
    if ratios.is_empty() {
        return 0.0;
    }
    let mut sum = 0.0;
    for &r in ratios {
        sum += r.ln();
    }
    (sum / ratios.len() as f64).exp()
}

/// One `(field, resolution)` measurement, before any cross-row fit exists.
struct Row {
    field: &'static str,
    samples: u32,
    cell_size: f64,
    beta_sum: f64,
    /// `(n, this scale's contribution, sub-boxes that carried a patch)`.
    per_scale: Vec<(u32, f64, u64)>,
    sum_converges: bool,
    actual_triangles: u64,
    surface_cells: u64,
    mesh_area: f64,
    mean_l1_normal: f64,
    field_samples: u64,
}

impl Row {
    /// The area-like sum turned into a count: `beta_sum / h^2`.
    fn beta_over_h2(&self) -> f64 {
        self.beta_sum / (self.cell_size * self.cell_size)
    }

    /// `M-13`'s quantity: `A / h^2`.
    fn area_over_h2(&self) -> f64 {
        self.mesh_area / (self.cell_size * self.cell_size)
    }

    /// The constant this row alone would need: `actual / (beta_sum / h^2)`.
    fn required_constant(&self) -> f64 {
        self.actual_triangles as f64 / self.beta_over_h2()
    }

    /// The constant this row alone would need for `M-13`'s form.
    fn required_m13_constant(&self) -> f64 {
        self.actual_triangles as f64 / self.area_over_h2()
    }

    /// The closed interval of constants that put this row inside C1's bar.
    fn constant_interval(&self) -> (f64, f64) {
        let r = self.required_constant();
        ((1.0 - C1_BAR) * r, (1.0 + C1_BAR) * r)
    }

    /// `beta_sum` per unit of measured surface area — the factor a single `k`
    /// has to be constant against.
    fn beta_per_area(&self) -> f64 {
        self.beta_sum / self.mesh_area
    }

    /// The last dyadic increment, and its ratio to the one before it.
    ///
    /// A ratio below one is a decaying tail, which is convergence; at or above
    /// one the truncated sum is still growing as scales are added, which is the
    /// divergence C2 is about.
    fn increment(&self) -> (f64, f64) {
        let n = self.per_scale.len();
        let last = self.per_scale[n - 1].1;
        let previous = self.per_scale[n - 2].1;
        (last, last / previous)
    }

    /// Thinnest-slab searches performed over the whole ladder.
    fn slab_fits(&self) -> u64 {
        self.per_scale.iter().map(|&(_, _, c)| c).sum()
    }

    /// Relative error of a prediction against this row's measured count.
    fn error_ratio(&self, predicted: i64) -> f64 {
        (predicted - self.actual_triangles as i64).abs() as f64 / self.actual_triangles as f64
    }
}

/// Measure one `(field, resolution)`.
fn measure<F>(name: &'static str, field: &F, samples: u32) -> Row
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();
    let (shape, origin, h) = common::grid::<f64, _>(field, samples);

    let scales = dyadic_ladder(samples - 1);
    let sum = common::beta::beta_sum(field, lo, hi, &scales);
    let sub_samples = u64::from(common::beta::SUB_SAMPLES - 1);
    let field_samples = scales
        .iter()
        .map(|&n| {
            let per_axis = u64::from(n) * sub_samples + 1;
            per_axis * per_axis * per_axis
        })
        .sum();

    let mut mesh = MeshBuffer::<f64>::new();
    let mut mc = MarchingCubes::<f64>::new();
    mc.extract(field, &shape, origin, h, &mut mesh)
        .expect("marching cubes over a reference field's own grid");
    let (mesh_area, mean_l1_normal) = area_and_mean_l1(&mesh);

    Row {
        field: name,
        samples,
        cell_size: h,
        beta_sum: sum.total,
        sum_converges: sum.converges(CONVERGENCE_TOLERANCE),
        per_scale: sum.per_scale.clone(),
        actual_triangles: mesh.triangle_count() as u64,
        surface_cells: surface_cells(field, &shape, origin, h),
        mesh_area,
        mean_l1_normal,
        field_samples,
    }
}

/// The largest number of fields any single constant can put inside C1's bar at
/// all three resolutions, and the interval that achieves it.
///
/// Each field contributes the intersection of its three rows' admissible
/// intervals; the answer is the maximum stabbing count over those closed
/// intervals, which is attained at some interval's left endpoint.
fn best_constant(per_field: &[(f64, f64)]) -> (usize, f64, f64) {
    let mut best = (0usize, f64::NAN, f64::NAN);
    for &(candidate, _) in per_field {
        if !candidate.is_finite() {
            continue;
        }
        let mut count = 0usize;
        let mut lo = f64::NEG_INFINITY;
        let mut hi = f64::INFINITY;
        for &(l, h) in per_field {
            if l <= candidate && candidate <= h {
                count += 1;
                lo = lo.max(l);
                hi = hi.min(h);
            }
        }
        if count > best.0 {
            best = (count, lo, hi);
        }
    }
    best
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-154");

    common::experiment::run(prereg, |run| {
        let mut rows: Vec<Row> = Vec::with_capacity(FIELDS * RESOLUTIONS.len());
        for samples in RESOLUTIONS {
            // Inline block per field, not a closure, so no `return` in here
            // (M-253).
            isomesh::for_each_reference_field!(f64, |name, field| {
                let row = measure(name, &field, samples);
                let (last, ratio) = row.increment();
                println!(
                    "  {:<15} {:>4}^3  scales {}  beta_sum {:>12.6}  last {:>11.6}  \
                     ratio {:>6.3}  conv {:<5}  tris {:>7}  cells {:>7}  area {:>10.4}  \
                     l1 {:>6.4}  k_row {:>10.3}",
                    row.field,
                    row.samples,
                    row.per_scale.len(),
                    row.beta_sum,
                    last,
                    ratio,
                    row.sum_converges,
                    row.actual_triangles,
                    row.surface_cells,
                    row.mesh_area,
                    row.mean_l1_normal,
                    row.required_constant(),
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
            // The registration's own vacuity control, verbatim: "at least four
            // dyadic scales must be summed, or convergence cannot be observed".
            assert!(
                row.per_scale.len() >= MIN_SCALES,
                "VOID: {} at {}^3 summed {} dyadic scales and the registration requires at \
                 least {MIN_SCALES}; with fewer, `converges` is reading one or two increments \
                 and convergence has not been observed",
                row.field,
                row.samples,
                row.per_scale.len()
            );
            for &(n, partial, counted) in &row.per_scale {
                assert!(
                    counted > 0,
                    "VOID: {} at {}^3 found no sub-box carrying a patch at scale n = {n}, so \
                     that level contributes a structural zero and the increment ratio through \
                     it is not a convergence reading",
                    row.field,
                    row.samples
                );
                assert!(
                    partial > 0.0,
                    "VOID: {} at {}^3 contributed {partial} at scale n = {n} over {counted} \
                     sub-boxes, and a zero increment makes the convergence ratio undefined \
                     rather than small",
                    row.field,
                    row.samples
                );
            }
            assert!(
                row.beta_sum > 0.0,
                "VOID: {} at {}^3 has beta_sum = {}, which predicts zero triangles for every \
                 constant — C1 would then be falsified by an undefined prediction rather than \
                 by a wrong one",
                row.field,
                row.samples,
                row.beta_sum
            );
            assert!(
                row.actual_triangles > 0,
                "VOID: {} at {}^3 extracted no triangles, so C1's denominator is zero and a \
                 budget for an empty mesh is not a budget",
                row.field,
                row.samples
            );
            assert!(
                row.mesh_area > 0.0,
                "VOID: {} at {}^3 has a mesh of zero area, so M-13's A/h^2 is zero and \
                 mean_l1_normal is a ratio of zeros",
                row.field,
                row.samples
            );
        }

        let mut spread_lo = f64::INFINITY;
        let mut spread_hi = 0.0f64;
        for row in &rows {
            spread_lo = spread_lo.min(row.beta_per_area());
            spread_hi = spread_hi.max(row.beta_per_area());
        }
        let beta_sum_spread = spread_hi / spread_lo;
        assert!(
            beta_sum_spread > MIN_SPREAD,
            "VOID: beta_sum per unit area spans only {beta_sum_spread:.6} across the run \
             ({spread_lo:e} to {spread_hi:e}), so a single constant converts one into the other \
             trivially and C1's verdict is a statement about a constant fixture rather than \
             about the eight fields"
        );

        // ── the two global fits, each one parameter over all 24 rows ──

        let beta_ratios: Vec<f64> = rows.iter().map(Row::required_constant).collect();
        let m13_ratios: Vec<f64> = rows.iter().map(Row::required_m13_constant).collect();
        let fitted_constant = geometric_mean(&beta_ratios);
        let m13_fitted_constant = geometric_mean(&m13_ratios);
        assert!(
            fitted_constant.is_finite() && fitted_constant > 0.0,
            "VOID: the fitted constant is {fitted_constant}, so the predicted_triangles column \
             is not a prediction"
        );
        assert!(
            m13_fitted_constant.is_finite() && m13_fitted_constant > 0.0,
            "VOID: M-13's refitted constant is {m13_fitted_constant}, so the matched-degrees-of\
             -freedom comparison has no incumbent"
        );

        // ── the four verdicts, each counted over fields rather than rows ──

        let names: Vec<&'static str> = {
            let mut seen: Vec<&'static str> = Vec::with_capacity(FIELDS);
            for row in &rows {
                if !seen.contains(&row.field) {
                    seen.push(row.field);
                }
            }
            seen
        };
        assert_eq!(
            names.len(),
            FIELDS,
            "VOID: the run saw {} distinct fields and the reference roster has {FIELDS}, so a \
             field-counting clause is being decided over the wrong population",
            names.len()
        );

        let predicted = |row: &Row| (fitted_constant * row.beta_over_h2()).round() as i64;
        let m13_published = |row: &Row| {
            (M13_TRIANGLES_PER_CELL * M13_CELLS_PER_AREA * row.area_over_h2()).round() as i64
        };
        let m13_refitted = |row: &Row| (m13_fitted_constant * row.area_over_h2()).round() as i64;
        let m13_anisotropic = |row: &Row| {
            (M13_TRIANGLES_PER_CELL * row.mean_l1_normal * row.area_over_h2()).round() as i64
        };

        let count_fields = |f: &dyn Fn(&Row) -> i64| {
            names
                .iter()
                .filter(|name| {
                    rows.iter()
                        .filter(|row| row.field == **name)
                        .all(|row| row.error_ratio(f(row)) <= C1_BAR)
                })
                .count()
        };
        let c1_fields = count_fields(&predicted);
        let m13_fields = count_fields(&m13_published);
        let m13_fitted_fields = count_fields(&m13_refitted);
        let m13_anisotropic_fields = count_fields(&m13_anisotropic);
        let c1_holds = c1_fields >= C1_MIN_FIELDS;

        // Per field: the intersection of its three rows' admissible intervals.
        let field_interval = |name: &'static str| {
            let mut lo = f64::NEG_INFINITY;
            let mut hi = f64::INFINITY;
            for row in rows.iter().filter(|row| row.field == name) {
                let (l, h) = row.constant_interval();
                lo = lo.max(l);
                hi = hi.min(h);
            }
            (lo, hi)
        };
        let per_field: Vec<(f64, f64)> = names.iter().map(|&name| field_interval(name)).collect();
        let (max_fields, best_lo, best_hi) = best_constant(&per_field);

        println!(
            "\n  fit: k = {fitted_constant:.6} over {} rows, 1 parameter, {} residual degrees \
             of freedom",
            rows.len(),
            rows.len() - 1
        );
        println!(
            "  beta-sum budget      within {C1_BAR:.2} on {c1_fields}/{FIELDS} fields \
             (needs {C1_MIN_FIELDS})  -> C1 {}",
            if c1_holds { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "  best possible k      any single constant reaches {max_fields}/{FIELDS} fields, \
             k in [{best_lo:.4}, {best_hi:.4}]"
        );
        println!("  M-13 as published    3*A/h^2, 0 parameters, {m13_fields}/{FIELDS} fields");
        println!(
            "  M-13 refitted        k13 = {m13_fitted_constant:.6}, 1 parameter, \
             {m13_fitted_fields}/{FIELDS} fields"
        );
        println!(
            "  M-13 anisotropic     2*mean_l1*A/h^2, 0 parameters, \
             {m13_anisotropic_fields}/{FIELDS} fields"
        );
        println!(
            "  beta_sum/area spread {beta_sum_spread:.3}x  ({spread_lo:.6} to {spread_hi:.6})"
        );

        for row in &rows {
            let (lo, hi) = field_interval(row.field);
            let (last, ratio) = row.increment();
            let pred = predicted(row);
            let m13 = m13_published(row);
            let m13_fit = m13_refitted(row);
            let m13_aniso = m13_anisotropic(row);
            let actual = row.actual_triangles as i64;
            // C2 is per row: the tail decays, or the truncated sum is still
            // growing as scales are added, which is divergence.
            let c2 = ratio < 1.0;
            let scales: Vec<String> = row
                .per_scale
                .iter()
                .map(|&(n, _, _)| n.to_string())
                .collect();
            let slab_fits = row.slab_fits();
            let extract_samples = u64::from(row.samples).pow(3);

            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.samples.to_string()),
                ("beta_sum", format!("{:.9}", row.beta_sum)),
                ("predicted_triangles", pred.to_string()),
                ("actual_triangles", actual.to_string()),
                ("prediction_error", (pred - actual).to_string()),
                (
                    "prediction_error_ratio",
                    format!("{:.6}", row.error_ratio(pred)),
                ),
                ("sum_converges", row.sum_converges.to_string()),
                ("scales_used", row.per_scale.len().to_string()),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2.to_string()),
                // ── extras (M-273) ──
                ("area_over_h2", format!("{:.6}", row.area_over_h2())),
                ("beta_over_h2", format!("{:.6}", row.beta_over_h2())),
                ("beta_sum_per_area", format!("{:.9}", row.beta_per_area())),
                ("beta_sum_spread", format!("{beta_sum_spread:.6}")),
                ("c1_bar", format!("{C1_BAR:.2}")),
                ("c1_fields_within_bar", c1_fields.to_string()),
                ("c1_min_fields", C1_MIN_FIELDS.to_string()),
                ("cell_size", format!("{:.9}", row.cell_size)),
                ("cells", (u64::from(row.samples - 1).pow(3)).to_string()),
                (
                    "convergence_tolerance",
                    format!("{CONVERGENCE_TOLERANCE:.2}"),
                ),
                ("dof_residual", (rows.len() - 1).to_string()),
                ("dof_spent", 1.to_string()),
                ("extract_samples", extract_samples.to_string()),
                ("field_constant_hi", format!("{hi:.6}")),
                ("field_constant_lo", format!("{lo:.6}")),
                ("field_feasible_any_constant", (lo <= hi).to_string()),
                ("field_samples", row.field_samples.to_string()),
                ("fit_rows", rows.len().to_string()),
                ("fitted_constant", format!("{fitted_constant:.6}")),
                ("increment_last", format!("{last:.9}")),
                ("increment_ratio_last", format!("{ratio:.6}")),
                ("max_fields_any_constant", max_fields.to_string()),
                (
                    "max_fields_constant_hi",
                    if best_hi.is_finite() {
                        format!("{best_hi:.6}")
                    } else {
                        String::from("none")
                    },
                ),
                (
                    "max_fields_constant_lo",
                    if best_lo.is_finite() {
                        format!("{best_lo:.6}")
                    } else {
                        String::from("none")
                    },
                ),
                ("mean_l1_normal", format!("{:.6}", row.mean_l1_normal)),
                ("mesh_area", format!("{:.6}", row.mesh_area)),
                (
                    "m13_anisotropic_error_ratio",
                    format!("{:.6}", row.error_ratio(m13_aniso)),
                ),
                (
                    "m13_anisotropic_fields_within_bar",
                    m13_anisotropic_fields.to_string(),
                ),
                ("m13_anisotropic_predicted_triangles", m13_aniso.to_string()),
                (
                    "m13_cells_constant",
                    format!("{:.6}", row.surface_cells as f64 / row.area_over_h2()),
                ),
                (
                    "m13_cells_error_ratio",
                    format!(
                        "{:.6}",
                        ((M13_CELLS_PER_AREA * row.area_over_h2()).round()
                            - row.surface_cells as f64)
                            .abs()
                            / row.surface_cells as f64
                    ),
                ),
                ("m13_error_ratio", format!("{:.6}", row.error_ratio(m13))),
                ("m13_fields_within_bar", m13_fields.to_string()),
                ("m13_fitted_constant", format!("{m13_fitted_constant:.6}")),
                (
                    "m13_fitted_error_ratio",
                    format!("{:.6}", row.error_ratio(m13_fit)),
                ),
                (
                    "m13_fitted_fields_within_bar",
                    m13_fitted_fields.to_string(),
                ),
                ("m13_fitted_predicted_triangles", m13_fit.to_string()),
                (
                    "m13_predicted_cells",
                    (M13_CELLS_PER_AREA * row.area_over_h2())
                        .round()
                        .to_string(),
                ),
                ("m13_predicted_triangles", m13.to_string()),
                (
                    "required_constant",
                    format!("{:.6}", row.required_constant()),
                ),
                (
                    "row_within_bar",
                    (row.error_ratio(pred) <= C1_BAR).to_string(),
                ),
                (
                    "samples_ratio",
                    format!("{:.3}", row.field_samples as f64 / extract_samples as f64),
                ),
                ("scales", scales.join("|")),
                (
                    "slab_evaluations",
                    (slab_fits * SLAB_EVALUATIONS_PER_PATCH).to_string(),
                ),
                (
                    "slab_evaluations_per_patch",
                    SLAB_EVALUATIONS_PER_PATCH.to_string(),
                ),
                ("slab_fits", slab_fits.to_string()),
                ("surface_cells", row.surface_cells.to_string()),
                (
                    "triangles_per_surface_cell",
                    format!(
                        "{:.6}",
                        row.actual_triangles as f64 / row.surface_cells as f64
                    ),
                ),
            ]);
        }
    });
}
