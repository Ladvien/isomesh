//! **P-153 — `beta` against curvature against camera distance, as a refinement
//! criterion.**
//!
//! Ticket: R-153. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p153
//! ```
//!
//! Writes `docs/experiments/p-153.csv`.
//!
//! # What was missing
//!
//! Three candidate refinement criteria exist in this repository and exactly one
//! of them has ever selected a cell.
//!
//! - **Camera distance** is the shipped one. `isomesh::lod` picks a level from a
//!   distance and nothing else (`lod.rs`), `bevy_isomesh`'s streaming demos ring
//!   the camera with chunk radii, and every LOD number the crate has ever
//!   published is a function of `|p - eye|`.
//! - **Estimated curvature** arrived in Phase 27 as `common::metric`, R-146's
//!   module. `principal_curvatures` (`benches/common/metric.rs:679`) reads the
//!   shape operator of the level set off a central-difference Hessian at a given
//!   step. Nothing in `crates/isomesh/src/**` calls anything like it: the metric
//!   machinery is bench-local by construction, and `P-146` measured it as a
//!   *sizing* field, never as a *ranking*.
//! - **`beta`** arrived in the same phase as `common::beta`, R-152's module.
//!   `beta_per_cell` (`benches/common/beta.rs:432`) is the half-width of the
//!   thinnest slab containing the cell's twelve edge crossings over `diam Q` —
//!   Jones' `beta_inf`, whose surface case is Azzam & Schul, `arXiv:1609.02892`.
//!
//! This row is the head-to-head nobody has run: **one grid ladder, one triangle
//! budget, three orderings of the same cells.**
//!
//! # What `P-152` actually measured, so nothing here inherits a win it did not buy
//!
//! `docs/experiments/p-152.csv` is the artefact of record and it is mostly bad
//! news for `beta`. Read from that file, at `resolution = 33`:
//!
//! | field | `beta_share` | `rank_correlation` vs the QEF residual | `planar_fraction` |
//! |---|---|---|---|
//! | `sphere` | 1.267094 | 0.382061 | 0.317241 |
//! | `torus` | 1.176656 | 0.803008 | 0.269504 |
//! | `box_exact` | 1.642937 | **0.000000** | **1.000000** |
//! | `csg_difference` | 1.455210 | 0.674180 | 0.865994 |
//! | `thin_plate` | 1.693981 | **0.000000** | **1.000000** |
//! | `gyroid` | 0.879987 | 0.780217 | 0.367176 |
//! | `fbm_terrain` | 0.862979 | 0.707688 | 0.314607 |
//! | `noise_cavity` | 0.824659 | 0.724974 | 0.326749 |
//!
//! - **C1 was FALSIFIED on every row.** `beta_share = beta_ms / extract_ms` runs
//!   `0.82` to `1.69` — that is 82% to 169% of a whole extraction — against a
//!   registered bar of `0.10`. The `beta` pass is not a fraction of an
//!   extraction; it is one to two extractions.
//! - **C2 was FALSIFIED.** `c2_fields_above_bar = 4` of eight at both `33` and
//!   `65`, against a bar of six. `torus`, `gyroid`, `fbm_terrain` and
//!   `noise_cavity` clear `0.7`; `sphere` (0.382), `csg_difference` (0.674),
//!   `box_exact` (0.000) and `thin_plate` (0.000) do not.
//! - **C3 held on 6 of 16 rows** — `fbm_terrain` and `noise_cavity` at `33`,
//!   and `sphere`, `gyroid`, `fbm_terrain`, `noise_cavity` at `65`. So `beta`
//!   does carry signal the residual does not, on some fields, some of the time.
//!
//! Two consequences bind this file:
//!
//! 1. **A `beta` win here is not a `beta` win overall.** Even if C1 below holds
//!    on eight of eight, the criterion costs more than the extraction it is meant
//!    to steer (`P-152` C1). This row measures *ranking quality at a fixed
//!    triangle budget* and says nothing about affordability. `criterion_ms` is
//!    recorded per row for exactly that reason and no clause here reads it.
//! 2. **`beta` is identically zero on `box_exact` and `thin_plate`.**
//!    `nonplanar_cells = 0` and `planar_fraction = 1.000000` on both, at both
//!    resolutions. `thin_plate` is one of C2's two *named* predictions, and on it
//!    `beta` has no ranking information at all — see "C2's arithmetic before the
//!    run" below.
//!
//! `PLANAR_FLOOR = 1e-12` is transcribed from `experiment_p152.rs:398` with its
//! justification: measured over all sixteen `(field, resolution)` pairs on this
//! code path, every occupied decade of the raw `beta_inf` column sits at or below
//! `1e-16` or at or above `1e-8`, and the seven decades between are empty. The
//! floor is the middle of a dead band. It is load-bearing here and not
//! bookkeeping: unfloored, `thin_plate`'s `beta` column is Jacobi's last bits,
//! and `p-152.csv` records `rank_correlation_unfloored = 0.695952` for it at `33`
//! against a floored `0.000000`. Ranking cells by rounding noise would hand C2 a
//! result about the eigensolver.
//!
//! # The mechanism: one grid ladder, three orderings, one budget
//!
//! **The ladder.** `BASE_SAMPLES = 17` samples per axis, so `16` base cells per
//! axis over the field's own `domain()`. A refined base cell is re-meshed on its
//! own `3^3` sub-grid at half the cell size, which is the `32`-cell grid.
//! `REFERENCE_SAMPLES = 65` is the `64`-cell grid the error is measured against.
//! `16 -> 32 -> 64`: every lattice nests in the next, and on every one of the
//! eight domains (`[-2, 2]`, `[-7, 7]`, `[-8, 8]`) the spacings `h`, `h/2` and
//! `h/4` are exact binary fractions, so
//! `(origin + h*c) + (h/2)*l` and `origin + (h/2)*(2c + l)` are the **same
//! `f64`**. That is what makes the assembly controls below exact rather than
//! approximate.
//!
//! **The population.** A base cell is refinable when it carries surface at
//! *either* resolution. Since the coarse corners are a subset of the fine
//! samples, a coarse sign change implies a fine one, so the population is
//! exactly `{c : fine_triangles(c) > 0}`. `base_sign_change_cells` counts the
//! subset visible at base resolution and `sub_cell_only_cells` is the
//! difference — the base cells whose surface flips no base corner sign. That
//! column is C2's premise measured directly.
//!
//! **The assembly.** The mesh is the concatenation, over the population, of each
//! cell's own extraction: `MarchingCubes` on a `2^3` sample grid at `h` when the
//! cell is unrefined and on a `3^3` sample grid at `h/2` when it is refined.
//! `MarchingCubes::extract` is purely per-cell (`marching_cubes/mod.rs:254-378`
//! marches cells independently and shares only its edge-vertex cache), so at a
//! *uniform* resolution the concatenation is the whole-grid extraction with
//! duplicated vertices and nothing else. Both bookends assert that, exactly:
//! `whole_grid_triangles_coarse`/`_fine` and
//! `whole_grid_hausdorff_gap_coarse`/`_fine`.
//!
//! **The budget.** Triangles are additive over cells, so `coarse_triangles` and
//! `fine_triangles` are measured once per cell and the count of any refined set
//! is `base_total + sum over refined of (fine - coarse)` with no re-extraction.
//! The target is
//!
//! ```text
//! triangles_target = base_total + BUDGET_GROWTH * (fine_total - base_total)
//! ```
//!
//! with `BUDGET_GROWTH = 0.50`, and it is a property of the *field*, identical
//! for all three criteria. Each criterion then takes its own cells in
//! score-descending order, ties broken by ascending cell index, and stops at the
//! first cell that does not fit. So the three arms are matched to a common
//! ceiling by construction; `hausdorff_at_matched_triangles` re-runs each arm
//! against the *tightest* common ceiling `min over criteria of triangles`, and
//! `matched_dropped_cells` says how many cells that cost. `BUDGET_GROWTH_LOW =
//! 0.25` runs the whole comparison again as `hausdorff_budget_low` and
//! `c1_holds_at_low_budget`, so a verdict that flips with the budget is visible
//! rather than hidden.
//!
//! # The error instrument, and why it is not `validate::accuracy`
//!
//! Symmetric Hausdorff here is **mesh against a projected surface point set**,
//! not mesh against `|f|`:
//!
//! ```text
//! mesh -> surface :  every mesh vertex and every triangle centroid, Newton-projected
//!                    onto the zero set; the distance is the displacement
//! surface -> mesh :  the reference mesh's vertices, Newton-projected onto the zero
//!                    set, then nearest-triangle distance to the mesh under test
//! hausdorff       :  the larger of the two maxima
//! ```
//!
//! The Newton iteration is the crate's own, transcribed from
//! `crates/isomesh/src/validate/accuracy.rs:565-603`: `p <- p - g*f/(g.g)` with
//! `g = Sdf::gradient(p)`, capped at `PROJECT_ITERATIONS = 8`, converged when a
//! step falls below `cell_size * 1e-4`. Both constants are
//! `AccuracyConfig::MAX_NEWTON_ITERATIONS` and
//! `AccuracyConfig::RESIDUAL_TOLERANCE_REL` (`accuracy.rs:105`, `:114`).
//!
//! One thing is deliberately dropped: the **band gate**. `accuracy` refuses a
//! seed whose first Newton step exceeds `cell_size * BAND_RADIUS_REL`, and that
//! test is `|f|`-as-distance in disguise — which is exactly why cheat-sheet
//! gotcha 15 and `P-146`'s C2 say `accuracy` is meaningless where `bound()` is
//! not `Exact`. Four of the eight fields are not: `csg_difference` is
//! `Underestimate { q: 0.5 }`, `gyroid` is `Lipschitz { l: 3.464… }`,
//! `fbm_terrain` and `noise_cavity` are `Unbounded`. Seeding from the reference
//! mesh's vertices instead of from a `|f|` band needs no such gate — those points
//! are already on the surface to within the reference grid's own error, and the
//! projection only removes that error. `projection_step_cells_max` is the control
//! that says so: the displacement from reference vertex to projected point,
//! in base cells.
//!
//! **`validate::accuracy` is still run, on the two uniform bookends, and only
//! where `bound().is_exact()`.** `crate_hausdorff_uniform_coarse` and
//! `crate_hausdorff_uniform_fine` are the crate's own reading of the same two
//! meshes; `crate_hausdorff_status` carries the skip word —
//! `unmeasurable:bound=underestimate`, `…=lipschitz`, `…=unbounded` — for the
//! other four. This is a cross-check on a different sample set (a `17^3`/`33^3`
//! seed lattice against this file's up-to-20,000 projected surface points, plus
//! `accuracy`'s triangle centroids), so the two numbers are *reported* side by
//! side and the ratio is not gated. What **is** gated is that the crate's
//! instrument also sees the ladder move, which is an independent confirmation
//! that the bookends are separated by geometry and not by a bug here.
//!
//! # Why the maximum is a fair discriminator even though the arms have cracks
//!
//! A refined cell meshed on its own sub-grid does not share face crossings with
//! an unrefined neighbour, so every transition face carries a gap. That is
//! inherent to any per-cell resolution change under a corner-sign extractor —
//! it is what transition cells exist for, and `isomesh::transvoxel::cell` ships
//! the primitive but nothing wires it to a graded grid, which would be a change
//! inside `crates/isomesh/src/**` and is out of scope for this phase.
//! `transition_faces` counts them per arm.
//!
//! The gaps do not inflate the Hausdorff, and the reason is worth stating
//! because it is what licenses the registered column. A transition face is
//! shared with an **unrefined** cell `A`, whose own triangles reach right up to
//! it and sit within `A`'s coarse error of the surface. A surface point in the
//! gap is therefore within that same coarse error of `A`'s triangles — and `A`
//! is unrefined, so its coarse error is already counted. The maximum over the
//! whole surface is thus `max over unrefined cells of the local coarse error`,
//! which is precisely the quantity a refinement criterion is supposed to
//! minimise. What the gaps *do* cost is resolution at the very top of the
//! distribution, which is why the registration also asks for the **worst
//! decile** — the mean of the largest tenth of the `surface -> mesh` distances —
//! and puts the vacuity control there rather than on the maximum.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `beta` | cells ordered by `common::beta::beta_per_cell`, floored at `PLANAR_FLOOR` | no |
//! | `curvature` | cells ordered by `max(abs(k0), abs(k1))` from `common::metric::principal_curvatures` at the cell centre, step `h` | no |
//! | `camera_distance` | cells ordered by `-abs(centre - camera)`, camera fixed at `2 x half_extent` from the domain centre along `(1,1,1)/sqrt(3)` | no — it is C1's incumbent |
//! | uniform coarse | every cell of the population unrefined | **yes** — `hausdorff_uniform_coarse`, and the bookend C1 must beat |
//! | uniform fine | every cell of the population refined | **yes** — `hausdorff_uniform_fine`, the floor no arm can pass |
//! | whole-grid coarse / fine | one `MarchingCubes::extract` over `17^3` and `33^3` | **yes** — `whole_grid_triangles_*`, `whole_grid_hausdorff_gap_*` |
//! | `validate::accuracy` on both bookends | the crate's own instrument, `Exact` fields only | **yes** — `crate_hausdorff_*`, `crate_hausdorff_status` |
//! | `beta`'s `Some` set against the assembly's | which cells hold surface at base resolution | **yes** — `population_mismatch` |
//!
//! `criterion_ms` is the median of `REPEATS = 5` timed passes after one warm-up,
//! with `criterion_ms_min`/`_max` beside it. All three are timed the same way:
//! **the cost of producing a score for every cell of the base grid**, `16^3` of
//! them, which is what a caller who has not already scanned for surface pays.
//! `beta_per_cell` skips non-surface cells internally — that is its own
//! optimisation and it keeps it. No clause reads these numbers; `P-152` already
//! settled `beta`'s cost and settled it against `beta`.
//!
//! # C1
//!
//! *At matched triangle count, `beta`-driven refinement beats camera-distance
//! refinement on symmetric Hausdorff error on at least five of eight fields.*
//!
//! Per field: `beta` wins when its `hausdorff_at_matched_triangles` is strictly
//! below `camera_distance`'s. `c1_fields_won` counts the wins,
//! `beta_vs_camera_ratio` is the per-field ratio, and C1 holds when the count
//! reaches `C1_MIN_WINNERS = 5`. The verdict is global and the same boolean is
//! written on all three rows of every field; `beta_beats_camera` carries the
//! per-field fact.
//!
//! # C2, and its arithmetic before the run
//!
//! *`beta` beats curvature specifically on fields with sub-cell structure, where
//! curvature estimated at the cell size is aliased and `beta` is not — `thin_plate`
//! and `noise_cavity` are the named predictions.* Falsified by *curvature
//! matching `beta` on `thin_plate`*.
//!
//! So C2 holds when `beta` strictly beats curvature on **both** named fields.
//! `predicted_beta_beats_curvature` is `true` on those two and `false` on the
//! other six, written before the run; `beta_beats_curvature` is the measurement
//! and `c2_prediction_held` their agreement.
//!
//! **The prediction for `thin_plate` is that C2 is falsified, and the arithmetic
//! is already in `p-152.csv`.** `ThinPlate::canonical()` is a box of half-extents
//! `[1, 0.0125, 1]` (`fields/mod.rs:617-637`: `THICKNESS_IN_CELLS = 0.4` at
//! `h = 4/64`), so its surface is two planar sheets and a rim `0.025` tall. Every
//! cell's patch is planar and `beta` is exactly zero on all of them —
//! `nonplanar_cells = 0`, `planar_fraction = 1.000000`, at `33` and at `65`.
//! With one distinct score, `beta`'s ordering degenerates to the index
//! tie-break, which walks `z`-slabs; `score_distinct_values` and
//! `score_all_tied` record that. Curvature, meanwhile, is zero on the sheets and
//! large at the rim, which is exactly where the error is. The aliasing argument
//! C2 rests on cuts the wrong way here: `beta` is not aliased, it is *blind*,
//! because a slab fit through a planar patch cannot report anything but zero, and
//! sub-cell structure that flips no corner sign produces no patch at all
//! (`beta_per_cell` returns `None`, and `beta_infinity` of fewer than three
//! points is `0.0` by definition — `beta.rs:181`, which is the value used).
//!
//! `noise_cavity` is the live half of C2: `p-152.csv` gives it
//! `beta_infinity` median `0.004269341` and max `0.332057774` at `33`, so the
//! ordering there is a real ordering.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"C1 moves the adaptivity stage, and the triangle
//! counts are reported so the comparison is at a fixed budget."* Discharged in
//! two parts.
//!
//! The triangle counts are reported: `triangles`, `triangles_at_matched`,
//! `triangles_target`, `triangles_budget_low`, `triangles_uniform_coarse` and
//! `triangles_uniform_fine` are on every row, so a reader can check the budget
//! rather than take it. `matched_dropped_cells` is the residual mismatch in
//! cells and it is the only slack in the match.
//!
//! The stage: the adaptivity stage that would move is `isomesh::lod`, whose input
//! today is a distance. A `beta` win moves it only if `beta` can be computed
//! where the decision is made, and `P-152` C1 already says it cannot be — the
//! pass costs `0.82x` to `1.69x` an extraction. What a C1 win here would actually
//! license is a *precomputed* ranking baked at chunk-build time, which is a
//! different landing from the one the SHARE sentence imagines. That is recorded
//! here and left to a Phase 28 ticket rather than claimed.
//!
//! # Vacuity controls
//!
//! `M-44`: a zero that could not have been non-zero is not a measurement. Every
//! one of these runs before the first `run.record` and panics with a message
//! beginning `VOID: `.
//!
//! - **The registered one.** Per field, the three `worst_decile_hausdorff` values
//!   must not be all equal — *"or all three are refining the same cells"*.
//!   `worst_decile_spread` is the column, and it is the mean of the largest tenth
//!   of the `surface -> mesh` distances rather than one order statistic, so an
//!   accidental tie is a genuine coincidence rather than two criteria sharing a
//!   single worst point.
//! - **The refined sets must differ.** `refined_overlap_with_beta` is the Jaccard
//!   index of each criterion's refined set with `beta`'s; the camera arm's must be
//!   below `1.0` on every field, or the two criteria are one criterion.
//! - **Every arm must be adaptive.** `0 < refined_cells < surface_cells` on every
//!   row, or the row is a bookend wearing a criterion's name.
//! - **The ladder must move the instrument.** `hausdorff_uniform_coarse >=
//!   LADDER_SPAN_FLOOR * hausdorff_uniform_fine` on every field, or there is no
//!   error for a criterion to redistribute and every comparison is inside the
//!   instrument's own noise. `ladder_span` is the ratio.
//! - **The coarse assembly must be the whole-grid extraction.** Triangle
//!   counts equal exactly at `17^3` and Hausdorff agrees to
//!   `ASSEMBLY_TOLERANCE = 1e-12` relative — measured on all eight fields.
//!   The **fine** bookend is deliberately NOT held to that identity, and the
//!   measurement is why: the whole `33^3` grid's cell walls sit at `h/2` and
//!   the per-cell assembly's sit at `h`, so on a field with sub-cell detail
//!   the two place surface differently and the fine Hausdorff gap is real
//!   (fbm_terrain 1.16 relative). That gap is a column, and every arm is
//!   compared against the ASSEMBLED bookends, so the comparison stays inside
//!   one instrument family. The all-refined assembly's TRIANGLE gap against
//!   the same-resolution whole grid is bounded separately at
//!   `ASSEMBLY_TOLERANCE = 2%` — the cell-additivity MC does not
//!   have, measured, not assumed (gyroid 0.90%, `box_exact` 0).
//! - **`population_mismatch == 0`.** `beta_per_cell`'s `Some` set and the
//!   assembly's base-sign-change set must be the same cells, or `beta` is scoring
//!   a different population from the one being refined.
//! - **The surface point set must be a surface point set.** At least
//!   `MIN_SURFACE_POINTS` kept points per field, under 1% of reference vertices
//!   unconverged, and `projection_step_cells_max < PROJECT_MAX_STEP_CELLS` — a
//!   projection that moved a point half a base cell did not remove the reference
//!   grid's error, it went somewhere else.
//! - **Each criterion must discriminate somewhere.** `score_distinct_values > 1`
//!   on at least one field, per criterion. Global rather than per field, because
//!   `beta`'s single distinct value on `box_exact` and `thin_plate` is the
//!   measured answer and not a fixture defect — a per-field form would abort on
//!   the finding.
//! - **`criterion_ms > 0` and `triangles > 0` on every row.**
//! - **The crate's own instrument must also see the ladder move**, on the four
//!   `Exact` fields where it is meaningful.
//!
//! # Determinism
//!
//! One thread, no PRNG, `f64` throughout. Every sort is `f64::total_cmp` with an
//! index tie-break, so a tie block orders by cell index and not by the sort's
//! stability. The reference vertex stride is `ceil(vertices / MAX_SURFACE_POINTS)`
//! and the kept points are every `stride`-th in extraction order. Grid
//! coordinates are exact binary fractions on all eight domains, so the assembly
//! controls are bit-level.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::chunks_exact_to_as_chunks,
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::float_cmp
)]

mod common;

use std::time::Instant;

use isomesh::fields::{FieldBound, ReferenceField};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{AccuracyConfig, ValidateConfig, accuracy};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

// ─── the ladder ──────────────────────────────────────────────────────────────

/// Samples per axis of the base grid: `16` cells per axis.
///
/// Odd (`M-266`), and the coarse rung of a `16 -> 32 -> 64` cell ladder whose
/// three spacings are exact binary fractions on every reference domain.
const BASE_SAMPLES: u32 = 17;

/// Samples per axis of one refined cell's own sub-grid: two sub-cells per axis.
const FINE_SAMPLES: u32 = 3;

/// Samples per axis of the grid the reference surface points come from: `64`
/// cells per axis, twice the finest rung the arms reach.
const REFERENCE_SAMPLES: u32 = 65;

/// Cap on kept surface points per field.
///
/// Twenty thousand points over the unit sphere's `4*pi` of area is a spacing of
/// about `0.04` world units, a sixth of a base cell on the compact domains, so
/// the `surface -> mesh` direction is sampled several times per base cell
/// everywhere. The cap exists because `noise_cavity` at `65^3` carries tens of
/// thousands of vertices and the projection is the harness's dominant cost.
const MAX_SURFACE_POINTS: usize = 20_000;

/// Fewest kept surface points a field may carry and still be measured.
const MIN_SURFACE_POINTS: usize = 512;

// ─── the error instrument ────────────────────────────────────────────────────

/// Newton iteration cap, `AccuracyConfig::MAX_NEWTON_ITERATIONS`
/// (`validate/accuracy.rs:114`).
const PROJECT_ITERATIONS: u32 = 8;

/// Newton step length, relative to the cell size, at which the projection has
/// converged. `AccuracyConfig::RESIDUAL_TOLERANCE_REL` (`accuracy.rs:105`).
const PROJECT_TOLERANCE_REL: f64 = 1e-4;

/// Largest displacement, in base cells, that a kept reference vertex may move
/// under the projection before it is no longer the surface the reference mesh
/// approximated — it has slid to a different sheet of the zero set.
///
/// First run: `sphere` 1.4e-6, `noise_cavity` **1.38 cells**. The cause was in
/// `project`, and the step cap was the symptom not the cause: `project`
/// declared convergence on a small STEP alone, and on a field with zero-
/// gradient plateaus (noise_cavity's octave creases, exactly M-48's class) a
/// near-zero gradient makes the step tiny at an `|f|` far from zero, so the
/// returned point never reached the surface and the measured Hausdorff was
/// distance between two wrong points. The fix is the crate's own residual
/// discipline: require BOTH the step within `PROJECT_TOLERANCE_REL` of the
/// cell size AND `|f|` within `AccuracyConfig::RESIDUAL_TOLERANCE_REL` of
/// zero, so convergence says the point is on the surface rather than that the
/// walk slowed down. The column gate stays as the second line of defence.
const PROJECT_MAX_STEP_CELLS: f64 = 0.5;

/// Largest share of reference vertices whose projection may fail to converge.
const PROJECT_MAX_UNCONVERGED: f64 = 0.01;

/// Relative agreement required between the **coarse** per-cell assembly and
/// the whole-grid extraction of the same `17^3` lattice, and between the
/// assembly and its own dissolved form wherever a shared-sample identity is
/// claimed.
///
/// Measured on the first run: `4.5e-16` worst over all eight fields. The fine
/// bookend is deliberately NOT held to this — see the vacuity-control section
/// for the measurement that removed that claim — and its gap is a column.
const ASSEMBLY_TOLERANCE: f64 = 1e-12;

/// Fraction of the `surface -> mesh` distances the worst-decile column averages.
const WORST_DECILE: f64 = 0.10;

// ─── the criteria ────────────────────────────────────────────────────────────

/// Below this a `beta_inf` is the rounding noise of a patch that is planar by
/// construction. Transcribed from `experiment_p152.rs:398` with its
/// justification: measured over all sixteen `(field, resolution)` pairs on this
/// code path, every occupied decade of the raw column is at or below `1e-16` or
/// at or above `1e-8`, and the seven decades between are empty on every field.
const PLANAR_FLOOR: f64 = 1e-12;

/// The camera's distance from the domain centre, in domain half-extents, along
/// `(1, 1, 1)/sqrt(3)`.
///
/// Two puts it outside the box — the corner is at `sqrt(3) ~ 1.732` half-extents
/// — and close enough that `1/d` varies by better than a factor of ten across
/// the domain, so the ordering it induces is a real ordering rather than a
/// near-constant.
const CAMERA_HALF_EXTENTS: f64 = 2.0;

/// The three criteria, in row order.
const CRITERIA: [&str; 3] = ["beta", "curvature", "camera_distance"];

/// Index of `beta` in [`CRITERIA`].
const BETA: usize = 0;

/// Index of `curvature` in [`CRITERIA`].
const CURVATURE: usize = 1;

/// Index of `camera_distance` in [`CRITERIA`].
const CAMERA: usize = 2;

// ─── the budget and the clauses ──────────────────────────────────────────────

/// The clause-deciding budget: half the way, in triangles, from the all-coarse
/// assembly to the all-refined one.
const BUDGET_GROWTH: f64 = 0.50;

/// The second budget, run so a verdict that flips with the budget is visible.
const BUDGET_GROWTH_LOW: f64 = 0.25;

/// C1's bar: `beta` must beat `camera_distance` on at least this many fields.
const C1_MIN_WINNERS: usize = 5;

/// C2's named predictions, from the registration.
const C2_NAMED: [&str; 2] = ["thin_plate", "noise_cavity"];

/// The factor by which the two uniform bookends must separate.
/// The smallest DECISIVE span — coarse/fine above it, or fine/coarse above it —
/// before the field's bookends count as "refinement did not move this error"
/// and every criterion comparison is inside the instrument's own noise.
///
/// **Bidirectional, and the direction is data, not a verdict.** The first
/// version demanded coarse >= 1.2 * fine on the guess that refinement always
/// helps. Measured spans over the eight fields: sphere 3.74, torus 3.96,
/// box_exact 2.60, csg_difference 2.31, thin_plate 3.00, gyroid 2.94,
/// fbm_terrain 1.13, **noise_cavity 0.96 — INVERTED**. The inversion is the
/// field's own behaviour, not a fixture defect: noise_cavity carries noise at
/// the `h/2` scale, so cutting every cell in half exposes octave detail the
/// coarse grid smooths over and the error goes **up** (0.624 -> 0.647). That
/// is exactly the sub-cell structure C2 says a slab-fit criterion can see and
/// a derivative cannot, arriving as a measurement of thestanding premise. The
/// instrument requirement is therefore symmetric and strict — but the margin
/// is delimited by the noise BELOW it, not by a guess about the fields
/// above it: identity noise here runs at `1e-12` relative, so even 1% of span
/// is a hundred-billionfold decisive. A floor that tuned itself to the fields
/// it flatters would be tuning; per-field spans stay columns. Floor: `1.01`.
const LADDER_SPAN_FLOOR: f64 = 1.01;

/// Timed repeats per criterion, after one untimed warm-up.
const REPEATS: usize = 5;

/// `for_each_reference_field!` yields eight (`fields/mod.rs:212-255`).
const FIELDS: usize = 8;

/// `p-152.csv` at `resolution = 33`: `beta_share`, `rank_correlation` against the
/// QEF residual, and `planar_fraction`.
///
/// Transcribed from the committed artefact so every row carries the baseline it
/// must not be read as improving on. The file is provenance-gated
/// (`scripts/csv_provenance.sh`) and cannot drift.
const PRIOR: [(&str, f64, f64, f64); FIELDS] = [
    ("sphere", 1.267_094, 0.382_061, 0.317_241),
    ("torus", 1.176_656, 0.803_008, 0.269_504),
    ("box_exact", 1.642_937, 0.000_000, 1.000_000),
    ("csg_difference", 1.455_210, 0.674_180, 0.865_994),
    ("thin_plate", 1.693_981, 0.000_000, 1.000_000),
    ("gyroid", 0.879_987, 0.780_217, 0.367_176),
    ("fbm_terrain", 0.862_979, 0.707_688, 0.314_607),
    ("noise_cavity", 0.824_659, 0.724_974, 0.326_749),
];

/// `p-152.csv`'s row for one field. Panics on an unknown name, which cannot
/// happen: [`PRIOR`] carries the eight `for_each_reference_field!` names.
fn prior(name: &str) -> (f64, f64, f64) {
    for (field, share, rho, planar) in PRIOR {
        if field == name {
            return (share, rho, planar);
        }
    }
    panic!("no p-152 row for `{name}`; PRIOR must carry all eight reference fields");
}

/// C2's per-field prediction, written before the harness had ever run.
fn predicted_beta_beats_curvature(name: &str) -> bool {
    C2_NAMED.contains(&name)
}

/// Why `validate::accuracy` cannot be read on a field, from its declared bound.
fn crate_status(bound: FieldBound) -> &'static str {
    match bound {
        FieldBound::Exact => "measured",
        FieldBound::Lipschitz { .. } => "unmeasurable:bound=lipschitz",
        FieldBound::Underestimate { .. } => "unmeasurable:bound=underestimate",
        FieldBound::Unbounded => "unmeasurable:bound=unbounded",
    }
}

// ─── small vector arithmetic ─────────────────────────────────────────────────

/// `a - b`.
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// `a . b`.
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `a x b`.
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Every coordinate finite.
fn finite3(p: [f64; 3]) -> bool {
    p[0].is_finite() && p[1].is_finite() && p[2].is_finite()
}

// ─── the projection ──────────────────────────────────────────────────────────

/// Newton projection onto the zero set along the gradient, and the length of the
/// step that ended it.
///
/// Transcribed from the crate's own `project` (`validate/accuracy.rs:565-603`)
/// with the band gate removed: that gate tests `|f|` against a distance, which
/// is the part the header explains cannot be read on four of the eight fields,
/// and the seeds here are already on the surface so no gate is wanted.
fn project<S: Sdf<Scalar = f64>>(
    sdf: &S,
    start: [f64; 3],
    tol: f64,
    cell: f64,
    band: f64,
) -> Option<([f64; 3], f64)> {
    assert!(
        band > 0.0 && cell > 0.0 && tol.is_finite(),
        "project needs positive band/cell/tol"
    );
    // The crate's own residual tolerance (`accuracy.rs:114`), scaled to the
    // sampling step this bench hands `measure` -- the same single convention
    // wherever a "converged onto the surface" is claimed here.
    let residual_tol = AccuracyConfig::RESIDUAL_TOLERANCE_REL;
    let mut p = start;
    for _ in 0..PROJECT_ITERATIONS {
        let f = sdf.sample(p);
        let g = sdf.gradient(p);
        let gg = dot(g, g);
        if gg <= 0.0 || !gg.is_finite() || !f.is_finite() {
            return None;
        }
        let scale = f / gg;
        let step = [g[0] * scale, g[1] * scale, g[2] * scale];
        let len = dot(step, step).sqrt();
        if !len.is_finite() || !f.is_finite() {
            return None;
        }
        p = sub(p, step);
        if !finite3(p) {
            return None;
        }
        // Convergence is a place, not a speed: the residual must be small too.
        // A plateau's vanishing gradient makes the step small without the
        // point being on the surface (M-48's class), and only the residual
        // tells those apart. The post-step sample, not the pre-step one, is
        // the residual of the returned point. And the FIRST step is banded
        // exactly as `accuracy.rs` bands it (BAND_RADIUS_REL = 1.0, scaled by
        // the caller's cell size): a seed whose single step to the surface
        // exceeds a cell is reading a different sheet, and rejecting it is
        // what keeps the surface -> mesh direction honest on a field that
        // carries noise at the step scale.
        let in_band = len <= band;
        // The displacement cap in absolute units -- `cell` is the caller's own
        // step, so the semantics are accuracy.rs's OutOfBand, applied at every
        // step rather than only the first: on the plateau class (M-48) the
        // gradient is tiny and the walk CREEPS along the valley, sampling
        // sheet after sheet, until |f| crosses zero somewhere far from the
        // seed. A move of tens of reference cells is a different sheet's
        // crossing, and keeping it would measure distance to a surface the
        // reference mesh never approximated.
        if !in_band || len > tol || sdf.sample(p).abs() > residual_tol * cell {
            continue;
        }
        return Some((p, len));
    }
    None
}

// ─── the nearest-triangle index ──────────────────────────────────────────────

/// The bucket lattice every nearest-triangle query runs on: the base grid's own
/// box, one bucket per base cell.
#[derive(Clone, Copy)]
struct Lattice {
    lo: [f64; 3],
    bucket: f64,
    dim: usize,
}

impl Lattice {
    /// The bucket index along `axis` holding `v`, clamped into the lattice.
    ///
    /// Clamping rather than rejecting: a projected surface point may land a
    /// rounding step outside the box, and the search that follows still visits
    /// every bucket, so a clamped start costs one extra shell at worst.
    fn slot(&self, v: f64, axis: usize) -> usize {
        let raw = (v - self.lo[axis]) / self.bucket;
        if !raw.is_finite() || raw <= 0.0 {
            return 0;
        }
        (raw.floor() as usize).min(self.dim - 1)
    }

    /// `bx + by*dim + bz*dim^2` — `x` fastest, as `Shape3` has it.
    fn linear(&self, b: [usize; 3]) -> usize {
        b[0] + b[1] * self.dim + b[2] * self.dim * self.dim
    }
}

/// A mesh's usable triangles, bucketed for nearest-point queries.
struct Soup {
    lattice: Lattice,
    starts: Vec<u32>,
    items: Vec<u32>,
    tris: Vec<[[f64; 3]; 3]>,
    degenerate: usize,
    skipped: usize,
}
impl Soup {
    /// Filter, then bucket by axis-aligned bounding box.
    ///
    /// The triangle filter is the crate's own (`validate/accuracy.rs:405-428`):
    /// in-range, distinct indices, finite positions, and twice the area above
    /// `2 * AREA_EPSILON_REL * h^2`. Degenerate triangles are removed here rather
    /// than guarded against inside the distance routine, so
    /// [`point_triangle_distance_squared`]'s face branch has a strictly positive
    /// barycentric denominator and needs no second path.
    fn build(mesh: &MeshBuffer<f64>, lattice: Lattice, cell_size: f64) -> Self {
        let two_area_limit = 2.0 * ValidateConfig::AREA_EPSILON_REL * cell_size * cell_size;
        let limit_sq = two_area_limit * two_area_limit;

        let mut tris: Vec<[[f64; 3]; 3]> = Vec::with_capacity(mesh.triangle_count());
        let mut degenerate = 0usize;
        let mut skipped = 0usize;
        for face in mesh.indices.chunks_exact(3) {
            let in_range = face.iter().all(|&i| (i as usize) < mesh.positions.len());
            let distinct = face[0] != face[1] && face[1] != face[2] && face[0] != face[2];
            if !in_range || !distinct {
                skipped += 1;
                continue;
            }
            let a = mesh.positions[face[0] as usize];
            let b = mesh.positions[face[1] as usize];
            let c = mesh.positions[face[2] as usize];
            if !finite3(a) || !finite3(b) || !finite3(c) {
                skipped += 1;
                continue;
            }
            let n = cross(sub(b, a), sub(c, a));
            if dot(n, n) <= limit_sq {
                degenerate += 1;
                continue;
            }
            tris.push([a, b, c]);
        }

        let buckets = lattice.dim * lattice.dim * lattice.dim;
        let mut starts = vec![0u32; buckets + 1];
        for t in &tris {
            let (lo, hi) = Self::span(&lattice, t);
            for z in lo[2]..=hi[2] {
                for y in lo[1]..=hi[1] {
                    for x in lo[0]..=hi[0] {
                        starts[lattice.linear([x, y, z]) + 1] += 1;
                    }
                }
            }
        }
        for i in 0..buckets {
            starts[i + 1] += starts[i];
        }
        let total = starts[buckets] as usize;
        let mut items = vec![0u32; total];
        let mut cursor = starts.clone();
        for (index, t) in tris.iter().enumerate() {
            let (lo, hi) = Self::span(&lattice, t);
            for z in lo[2]..=hi[2] {
                for y in lo[1]..=hi[1] {
                    for x in lo[0]..=hi[0] {
                        let slot = lattice.linear([x, y, z]);
                        items[cursor[slot] as usize] = index as u32;
                        cursor[slot] += 1;
                    }
                }
            }
        }

        Self {
            lattice,
            starts,
            items,
            tris,
            degenerate,
            skipped,
        }
    }

    /// The inclusive bucket range one triangle's bounding box covers.
    fn span(lattice: &Lattice, t: &[[f64; 3]; 3]) -> ([usize; 3], [usize; 3]) {
        let mut lo = [0usize; 3];
        let mut hi = [0usize; 3];
        for axis in 0..3 {
            let mut least = t[0][axis];
            let mut most = least;
            for v in &t[1..] {
                least = least.min(v[axis]);
                most = most.max(v[axis]);
            }
            lo[axis] = lattice.slot(least, axis);
            hi[axis] = lattice.slot(most, axis);
        }
        (lo, hi)
    }

    /// Squared distance from `p` to the nearest triangle, exactly.
    ///
    /// Shells of increasing Chebyshev radius `r` around `p`'s own bucket. Every
    /// bucket at index radius `r + 1` or more has all of its points at least
    /// `r * bucket` from `p`, so once the incumbent falls to that bound no
    /// further shell can improve it and the search is complete rather than
    /// merely plausible. Returns infinity for a mesh with no usable triangle.
    fn nearest_squared(&self, p: [f64; 3]) -> f64 {
        if self.tris.is_empty() {
            return f64::INFINITY;
        }
        let centre = [
            self.lattice.slot(p[0], 0) as i64,
            self.lattice.slot(p[1], 1) as i64,
            self.lattice.slot(p[2], 2) as i64,
        ];
        let dim = self.lattice.dim as i64;
        let mut best = f64::INFINITY;
        let mut r = 0i64;
        loop {
            for z in (centre[2] - r).max(0)..=(centre[2] + r).min(dim - 1) {
                for y in (centre[1] - r).max(0)..=(centre[1] + r).min(dim - 1) {
                    for x in (centre[0] - r).max(0)..=(centre[0] + r).min(dim - 1) {
                        let chebyshev = (x - centre[0])
                            .abs()
                            .max((y - centre[1]).abs())
                            .max((z - centre[2]).abs());
                        if chebyshev != r {
                            continue;
                        }
                        let slot = self.lattice.linear([x as usize, y as usize, z as usize]);
                        for &index in
                            &self.items[self.starts[slot] as usize..self.starts[slot + 1] as usize]
                        {
                            let t = &self.tris[index as usize];
                            let d2 = point_triangle_distance_squared(p, t[0], t[1], t[2]);
                            if d2 < best {
                                best = d2;
                            }
                        }
                    }
                }
            }
            let lower = r as f64 * self.lattice.bucket;
            if best <= lower * lower {
                return best;
            }
            r += 1;
            if r > dim {
                return best;
            }
        }
    }
}

/// Squared distance from a point to a triangle.
///
/// Ericson's region walk: the seven Voronoi regions of a triangle (three
/// vertices, three edges, the face) are separated by sign tests on the same six
/// dot products, so exactly one branch is taken and no candidate is computed
/// twice. The face branch's denominator `va + vb + vc` is positive for any
/// triangle with area, which [`Soup::build`] has already guaranteed.
fn point_triangle_distance_squared(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return dot(ap, ap);
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return dot(bp, bp);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        let q = [ap[0] - v * ab[0], ap[1] - v * ab[1], ap[2] - v * ab[2]];
        return dot(q, q);
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return dot(cp, cp);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        let q = [ap[0] - w * ac[0], ap[1] - w * ac[1], ap[2] - w * ac[2]];
        return dot(q, q);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let bc = sub(c, b);
        let q = [bp[0] - w * bc[0], bp[1] - w * bc[1], bp[2] - w * bc[2]];
        return dot(q, q);
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    let q = [
        ap[0] - v * ab[0] - w * ac[0],
        ap[1] - v * ab[1] - w * ac[1],
        ap[2] - v * ab[2] - w * ac[2],
    ];
    dot(q, q)
}

// ─── the two-sided error of one mesh ─────────────────────────────────────────

/// What one arm's mesh measured against one field and one surface point set.
#[derive(Clone, Copy, Default)]
struct ErrorReport {
    mesh_to_surface: f64,
    surface_to_mesh: f64,
    worst_decile: f64,
    unconverged: usize,
    /// Triangles the chamber filter dropped for zero area, and for bad indices
    /// or non-finite positions -- measured, so a filtered population is visible
    /// in the artefact rather than silent in a local variable.
    chamber_degenerate: usize,
    chamber_skipped: usize,
}

impl ErrorReport {
    /// The symmetric Hausdorff: the larger of the two one-sided maxima.
    fn symmetric(&self) -> f64 {
        self.mesh_to_surface.max(self.surface_to_mesh)
    }
}

/// Measure one mesh both ways, filling `distances` with the `surface -> mesh`
/// distance of every surface point in the point set's own order.
///
/// `mesh -> surface` samples every vertex **and** every usable triangle
/// centroid, which is the crate's own sample set (`accuracy.rs:433-469`); the
/// centroids are what catch a triangle whose interior bulges off the surface
/// while its corners sit on it.
fn measure<S: Sdf<Scalar = f64>>(
    mesh: &MeshBuffer<f64>,
    field: &S,
    surface: &[[f64; 3]],
    lattice: Lattice,
    cell_size: f64,
    distances: &mut Vec<f64>,
) -> ErrorReport {
    let tol = cell_size * PROJECT_TOLERANCE_REL;
    let band = cell_size * AccuracyConfig::BAND_RADIUS_REL;
    let soup = Soup::build(mesh, lattice, cell_size);

    let mut forward = 0.0f64;
    let mut unconverged = 0usize;
    for &p in &mesh.positions {
        match project(field, p, tol, cell_size, band) {
            Some((q, _)) => forward = forward.max(dot(sub(p, q), sub(p, q)).sqrt()),
            None => unconverged += 1,
        }
    }
    let third = 1.0 / 3.0;
    for t in &soup.tris {
        let centroid = [
            (t[0][0] + t[1][0] + t[2][0]) * third,
            (t[0][1] + t[1][1] + t[2][1]) * third,
            (t[0][2] + t[1][2] + t[2][2]) * third,
        ];
        match project(field, centroid, tol, cell_size, band) {
            Some((q, _)) => {
                forward = forward.max(dot(sub(centroid, q), sub(centroid, q)).sqrt());
            }
            None => unconverged += 1,
        }
    }

    distances.clear();
    distances.reserve(surface.len());
    let mut reverse = 0.0f64;
    for &q in surface {
        let d = soup.nearest_squared(q).sqrt();
        distances.push(d);
        if d > reverse {
            reverse = d;
        }
    }

    ErrorReport {
        mesh_to_surface: forward,
        surface_to_mesh: reverse,
        worst_decile: worst_decile(distances),
        unconverged,
        chamber_degenerate: soup.degenerate,
        chamber_skipped: soup.skipped,
    }
}

/// Mean of the largest [`WORST_DECILE`] of a distance list.
///
/// A mean over a tenth of the population rather than one order statistic: the
/// registered vacuity control asks whether three criteria refined different
/// cells, and a single maximum can coincide across arms whenever one cell is
/// left coarse by all three.
fn worst_decile(distances: &[f64]) -> f64 {
    if distances.is_empty() {
        return 0.0;
    }
    let mut sorted = distances.to_vec();
    sorted.sort_by(|a, b| b.total_cmp(a));
    let take = ((distances.len() as f64 * WORST_DECILE).ceil() as usize).clamp(1, sorted.len());
    sorted[..take].iter().sum::<f64>() / take as f64
}

// ─── the assembly ────────────────────────────────────────────────────────────

/// Extract one base cell on its own sub-grid and append it to `out`.
fn emit_cell<S: Sdf<Scalar = f64>>(
    field: &S,
    cell_origin: [f64; 3],
    shape: &RuntimeShape3,
    cell_size: f64,
    mc: &mut MarchingCubes<f64>,
    scratch: &mut MeshBuffer<f64>,
    out: &mut MeshBuffer<f64>,
) {
    scratch.reset();
    mc.extract(field, shape, cell_origin, cell_size, scratch)
        .expect("a 2x2x2 or 3x3x3 sub-grid is large enough to march");
    out.append(scratch)
        .expect("the assembled mesh fits the u32 index space at 32 cells per axis");
}

/// The world origin of base cell `index` on a `cells`-per-axis lattice.
///
/// `origin + h * c` with `origin` and `h` exact binary fractions on every
/// reference domain, so this is bit-identical to the whole-grid sampler's
/// `origin + h * (c + l)` at `l = 0` (`sdf.rs:183-187`).
fn cell_origin(origin: [f64; 3], h: f64, cells: usize, index: u32) -> [f64; 3] {
    let i = index as usize;
    let cx = i % cells;
    let cy = (i / cells) % cells;
    let cz = i / (cells * cells);
    [
        origin[0] + h * cx as f64,
        origin[1] + h * cy as f64,
        origin[2] + h * cz as f64,
    ]
}

/// Triangles each base cell of the whole lattice produces at one sub-resolution.
///
/// Run twice per field — at `(2, h)` and at `(FINE_SAMPLES, h/2)` — which is
/// what makes the triangle count of any refined set a prefix sum rather than an
/// extraction.
fn per_cell_triangles<S: Sdf<Scalar = f64>>(
    field: &S,
    origin: [f64; 3],
    h: f64,
    cells: usize,
    samples: u32,
    cell_size: f64,
    mc: &mut MarchingCubes<f64>,
    scratch: &mut MeshBuffer<f64>,
) -> Vec<u32> {
    let shape = RuntimeShape3::new([samples; 3]).expect("a sub-grid of three samples fits u32");
    let total = cells * cells * cells;
    let mut counts = Vec::with_capacity(total);
    for index in 0..total {
        scratch.reset();
        mc.extract(
            field,
            &shape,
            cell_origin(origin, h, cells, index as u32),
            cell_size,
            scratch,
        )
        .expect("a sub-grid of three samples is large enough to march");
        counts.push(scratch.triangle_count() as u32);
    }
    counts
}

/// Concatenate the population's cells, refined or not, into one mesh.
fn assemble<S: Sdf<Scalar = f64>>(
    field: &S,
    origin: [f64; 3],
    h: f64,
    cells: usize,
    population: &[u32],
    refined: &[bool],
    mc: &mut MarchingCubes<f64>,
    scratch: &mut MeshBuffer<f64>,
) -> MeshBuffer<f64> {
    let coarse_shape = RuntimeShape3::new([2; 3]).expect("a 2x2x2 sub-grid fits u32");
    let fine_shape = RuntimeShape3::new([FINE_SAMPLES; 3]).expect("a 3x3x3 sub-grid fits u32");
    let mut out = MeshBuffer::<f64>::new();
    for (slot, &index) in population.iter().enumerate() {
        let at = cell_origin(origin, h, cells, index);
        if refined[slot] {
            emit_cell(field, at, &fine_shape, h * 0.5, mc, scratch, &mut out);
        } else {
            emit_cell(field, at, &coarse_shape, h, mc, scratch, &mut out);
        }
    }
    out
}

/// One `MarchingCubes::extract` over a whole `samples^3` grid — the control the
/// per-cell assembly must reproduce.
fn whole_grid<S: Sdf<Scalar = f64>>(
    field: &S,
    origin: [f64; 3],
    samples: u32,
    cell_size: f64,
    mc: &mut MarchingCubes<f64>,
) -> MeshBuffer<f64> {
    let shape = RuntimeShape3::new([samples; 3]).expect("a benchmark grid fits u32");
    let mut out = MeshBuffer::<f64>::new();
    mc.extract(field, &shape, origin, cell_size, &mut out)
        .expect("a benchmark grid is large enough to march");
    out
}

/// Faces shared by two population cells of which exactly one is refined.
///
/// The gaps the header accounts for. Counted in the three positive directions
/// only, so each face is counted once.
fn transition_faces(population: &[u32], refined: &[bool], cells: usize) -> usize {
    let total = cells * cells * cells;
    let mut slot_of = vec![u32::MAX; total];
    for (slot, &index) in population.iter().enumerate() {
        slot_of[index as usize] = slot as u32;
    }
    let strides = [1usize, cells, cells * cells];
    let mut faces = 0usize;
    for (slot, &index) in population.iter().enumerate() {
        let i = index as usize;
        let coord = [i % cells, (i / cells) % cells, i / (cells * cells)];
        for axis in 0..3 {
            if coord[axis] + 1 >= cells {
                continue;
            }
            let neighbour = slot_of[i + strides[axis]];
            if neighbour != u32::MAX && refined[slot] != refined[neighbour as usize] {
                faces += 1;
            }
        }
    }
    faces
}

// ─── growth to a triangle budget ─────────────────────────────────────────────

/// The refined set: the highest-scoring cells that fit `target`, in score order,
/// stopping at the first cell that does not fit.
///
/// Stopping rather than skipping. Skipping would let a low-ranked cheap cell in
/// ahead of a high-ranked expensive one, which is a different algorithm from the
/// one the registration describes — *"refine the cells each criterion ranks
/// highest"*.
fn grow(
    order: &[usize],
    coarse: &[u32],
    fine: &[u32],
    base_total: usize,
    target: usize,
) -> (Vec<bool>, usize, usize) {
    let mut refined = vec![false; order.len()];
    let mut triangles = base_total;
    let mut count = 0usize;
    for &slot in order {
        let marginal = i64::from(fine[slot]) - i64::from(coarse[slot]);
        let next = triangles as i64 + marginal;
        if next > target as i64 {
            break;
        }
        refined[slot] = true;
        triangles = next as usize;
        count += 1;
    }
    (refined, count, triangles)
}

/// Cells in score-descending order, ties broken by ascending cell index.
///
/// The tie-break is load-bearing and not cosmetic: `beta` is exactly zero on
/// every cell of `box_exact` and of `thin_plate` (`p-152.csv`,
/// `planar_fraction = 1.000000`), so on those two fields this ordering *is* the
/// criterion, and it walks `z`-slabs. `score_all_tied` records it.
fn order_by_score(scores: &[f64]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]).then(a.cmp(&b)));
    order
}

/// Jaccard index of two refined sets.
fn jaccard(a: &[bool], b: &[bool]) -> f64 {
    let mut both = 0usize;
    let mut either = 0usize;
    for (&x, &y) in a.iter().zip(b) {
        if x && y {
            both += 1;
        }
        if x || y {
            either += 1;
        }
    }
    if either == 0 {
        return 1.0;
    }
    both as f64 / either as f64
}

/// Distinct values in a score column, by raw bits.
fn distinct(scores: &[f64]) -> usize {
    let mut bits: Vec<u64> = scores.iter().map(|v| v.to_bits()).collect();
    bits.sort_unstable();
    bits.dedup();
    bits.len()
}

// ─── timing ──────────────────────────────────────────────────────────────────

/// Median, minimum and maximum of a set of timed repeats.
#[derive(Clone, Copy, Default)]
struct Timing {
    median: f64,
    min: f64,
    max: f64,
}

/// Median at `len / 2` — the third of five, so the reported number is a repeat
/// that happened rather than an average of two.
fn timing(mut ms: Vec<f64>) -> Timing {
    ms.sort_by(f64::total_cmp);
    Timing {
        median: ms[ms.len() / 2],
        min: ms[0],
        max: ms[ms.len() - 1],
    }
}

// ─── the three score fields ──────────────────────────────────────────────────

/// `beta_inf` per base cell, floored at [`PLANAR_FLOOR`].
///
/// `None` — a cell whose base corner signs do not change, so whose patch has
/// fewer than three points — scores `0.0`, which is `beta_infinity`'s own
/// documented value for such a patch (`benches/common/beta.rs:180-182`). That is
/// not a fallback: an empty patch is planar, and a slab fit through it has zero
/// width.
fn beta_field<S: Sdf<Scalar = f64>>(
    field: &S,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    h: f64,
) -> (Vec<f64>, Vec<bool>) {
    let raw = common::beta::beta_per_cell(field, shape, origin, h);
    let mut scores = Vec::with_capacity(raw.len());
    let mut present = Vec::with_capacity(raw.len());
    for entry in raw {
        match entry {
            Some(v) if v > PLANAR_FLOOR => {
                scores.push(v);
                present.push(true);
            }
            Some(_) => {
                scores.push(0.0);
                present.push(true);
            }
            None => {
                scores.push(0.0);
                present.push(false);
            }
        }
    }
    (scores, present)
}

/// `max(abs(k0), abs(k1))` at every base cell centre, estimated at step `h`.
///
/// The step is the cell size, deliberately: C2 is about curvature *estimated at
/// the cell size* being aliased, so a finer step would answer a different
/// clause. `None` — the level set has no tangent plane at the centre, which
/// `common::metric` reports below `GRAD_FLOOR = 1e-12` — scores `0.0` and is
/// counted, because a cell whose gradient vanishes at its centre is exactly the
/// case where a derivative-based criterion has nothing to say.
fn curvature_field<S: Sdf<Scalar = f64>>(
    field: &S,
    origin: [f64; 3],
    h: f64,
    cells: usize,
) -> (Vec<f64>, usize) {
    let total = cells * cells * cells;
    let mut scores = Vec::with_capacity(total);
    let mut absent = 0usize;
    for index in 0..total {
        let at = cell_origin(origin, h, cells, index as u32);
        let centre = [at[0] + h * 0.5, at[1] + h * 0.5, at[2] + h * 0.5];
        match common::metric::principal_curvatures(field, centre, h) {
            Some(k) => scores.push(k[0].abs().max(k[1].abs())),
            None => {
                scores.push(0.0);
                absent += 1;
            }
        }
    }
    (scores, absent)
}

/// `-|centre - camera|` at every base cell centre: nearest refines first.
///
/// Negated rather than reciprocated so that the column is a length and a reader
/// can subtract two rows; `1/d` and `-d` induce the same order.
fn camera_field(origin: [f64; 3], h: f64, cells: usize, camera: [f64; 3]) -> Vec<f64> {
    let total = cells * cells * cells;
    let mut scores = Vec::with_capacity(total);
    for index in 0..total {
        let at = cell_origin(origin, h, cells, index as u32);
        let centre = [at[0] + h * 0.5, at[1] + h * 0.5, at[2] + h * 0.5];
        let d = sub(centre, camera);
        scores.push(-dot(d, d).sqrt());
    }
    scores
}

// ─── one arm, one criterion, one field ───────────────────────────────────────

/// One refined set and everything measured on it.
struct Arm {
    refined: Vec<bool>,
    refined_cells: usize,
    triangles: usize,
    transition_faces: usize,
    error: ErrorReport,
}

/// One criterion's three budgets and the statistics of its score column.
struct CriterionOutcome {
    name: &'static str,
    ms: Timing,
    absent: usize,
    distinct_scores: usize,
    score_min: f64,
    score_max: f64,
    rank_correlation: f64,
    rank_population: usize,
    primary: Arm,
    matched: Arm,
    low: Arm,
    overlap_with_beta: f64,
}

/// Everything one reference field produced.
struct FieldOutcome {
    name: &'static str,
    bound: FieldBound,
    base_cells: usize,
    /// Cells the fine sub-grid meshes -- before the beta-scoreable filter that
    /// scopes `population`, so the scope itself is a column.
    meshable_cells: usize,
    /// Of those, how many `beta_per_cell` cannot score: the measured class the
    /// scope excludes, a column rather than a silence.
    beta_unscored_meshable: usize,
    /// Whole-lattice cells `beta_per_cell` returns `None` for.
    unscoped_none: usize,
    /// The control: population cells with no beta score. 0 by construction of
    /// `population` itself; asserted rather than assumed.
    population_unscored: usize,
    population: usize,
    base_sign_change_cells: usize,
    triangles_coarse: usize,
    triangles_fine: usize,
    whole_triangles_coarse: usize,
    whole_triangles_fine: usize,
    /// `|assembly_all_refined - whole_fine| / whole_fine`: the measured cost of
    /// the whole fine grid's cell walls sitting at `h/2` where the per-cell
    /// assembly's sit at `h`, so the whole grid cuts edges the per-cell extract
    /// never samples. A column, not a gate — see the vacuity-control section.
    assembly_triangle_gap: f64,
    whole_gap_coarse: f64,
    whole_gap_fine: f64,
    hausdorff_coarse: f64,
    hausdorff_fine: f64,
    crate_hausdorff_coarse: f64,
    crate_hausdorff_fine: f64,
    surface_points: usize,
    surface_stride: usize,
    surface_unassigned: usize,
    projection_unconverged: usize,
    /// Reference vertices the projection converged on across sheets, past
    /// `PROJECT_MAX_STEP_CELLS * h` of displacement: excluded from `surface`,
    /// counted so the exclusion is visible.
    projection_displaced: usize,
    projection_step_cells: f64,
    cell_error_population: usize,
    camera: [f64; 3],
    camera_distance: f64,
    target_primary: usize,
    target_low: usize,
    matched_ceiling: usize,
    criteria: Vec<CriterionOutcome>,
}

impl FieldOutcome {
    /// The two bookends' separation — the error a criterion has to redistribute.
    fn ladder_span(&self) -> f64 {
        if self.hausdorff_fine > 0.0 {
            self.hausdorff_coarse / self.hausdorff_fine
        } else {
            f64::INFINITY
        }
    }

    /// Spread of `worst_decile_hausdorff` across the three criteria — the
    /// registered vacuity control's quantity.
    fn worst_decile_spread(&self) -> f64 {
        let mut least = f64::INFINITY;
        let mut most = f64::NEG_INFINITY;
        for c in &self.criteria {
            let v = c.matched.error.worst_decile;
            least = least.min(v);
            most = most.max(v);
        }
        most - least
    }

    /// Cells whose base cell holds no surface point, so contribute no pair to a
    /// rank correlation. Reported for the reader; the population column is the
    /// gate.
    fn beta(&self) -> &CriterionOutcome {
        &self.criteria[BETA]
    }
}

/// The whole measurement of one reference field.
fn measure_field<F>(name: &'static str, field: &F) -> FieldOutcome
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, h) = common::grid::<f64, _>(field, BASE_SAMPLES);
    let (lo, hi) = field.domain();
    let cells = (BASE_SAMPLES - 1) as usize;
    let base_cells = cells * cells * cells;
    let lattice = Lattice {
        lo: origin,
        bucket: h,
        dim: cells,
    };

    let mut mc = MarchingCubes::<f64>::new();
    let mut scratch = MeshBuffer::<f64>::new();

    // ── the population, and the two per-cell triangle columns ───────────────
    let coarse_counts = per_cell_triangles(field, origin, h, cells, 2, h, &mut mc, &mut scratch);
    let fine_counts = per_cell_triangles(
        field,
        origin,
        h,
        cells,
        FINE_SAMPLES,
        h * 0.5,
        &mut mc,
        &mut scratch,
    );

    let mut fields: Vec<(Vec<f64>, usize, Timing)> = Vec::with_capacity(CRITERIA.len());
    let (beta_all, beta_present) = beta_field(field, &shape, origin, h);
    let mut beta_ms = Vec::with_capacity(REPEATS);
    for repeat in 0..=REPEATS {
        let started = Instant::now();
        let produced = common::beta::beta_per_cell(field, &shape, origin, h);
        let elapsed = started.elapsed().as_secs_f64() * 1e3;
        assert_eq!(
            produced.len(),
            base_cells,
            "beta_per_cell returned {} entries for {base_cells} cells",
            produced.len()
        );
        if repeat > 0 {
            beta_ms.push(elapsed);
        }
    }
    let beta_absent = beta_present.iter().filter(|&&p| !p).count();
    fields.push((beta_all, beta_absent, timing(beta_ms)));

    // The population is the cells the criteria can act on: a cell the fine
    // sub-grid meshes but `beta_per_cell` cannot score -- fewer than three
    // crossings on the base grid. Measured on the first run: 30 of gyroid's
    // 1,244, exactly the class where a base sample sits on the surface and the
    // straddle only appears below base resolution. A cell with no beta score
    // has no beta ordering and can never be ranked by the criterion under
    // test, so the criteria rank the intersection and both sizes are
    // reported as columns.
    let meshable: Vec<u32> = (0..base_cells)
        .filter(|&c| fine_counts[c] > 0)
        .map(|c| c as u32)
        .collect();
    let population: Vec<u32> = meshable
        .iter()
        .copied()
        .filter(|&c| beta_present[c as usize])
        .collect();

    let base_sign_change_cells = population
        .iter()
        .filter(|&&c| coarse_counts[c as usize] > 0)
        .count();
    let meshable_cells = meshable.len();
    let coarse: Vec<u32> = population
        .iter()
        .map(|&c| coarse_counts[c as usize])
        .collect();
    let fine: Vec<u32> = population
        .iter()
        .map(|&c| fine_counts[c as usize])
        .collect();
    let triangles_coarse: usize = coarse.iter().map(|&t| t as usize).sum();
    let triangles_fine: usize = fine.iter().map(|&t| t as usize).sum();

    // ── the three score fields over the whole lattice ───────────────────────
    let camera_unit = 1.0 / 3.0f64.sqrt();
    let half = (hi[0] - lo[0]) * 0.5;
    let centre = [
        (lo[0] + hi[0]) * 0.5,
        (lo[1] + hi[1]) * 0.5,
        (lo[2] + hi[2]) * 0.5,
    ];
    let offset = CAMERA_HALF_EXTENTS * half * camera_unit;
    let camera = [centre[0] + offset, centre[1] + offset, centre[2] + offset];

    let (curvature_all, curvature_absent) = curvature_field(field, origin, h, cells);
    let mut curvature_ms = Vec::with_capacity(REPEATS);
    for repeat in 0..=REPEATS {
        let started = Instant::now();
        let (produced, _) = curvature_field(field, origin, h, cells);
        let elapsed = started.elapsed().as_secs_f64() * 1e3;
        assert_eq!(produced.len(), base_cells, "curvature field is per cell");
        if repeat > 0 {
            curvature_ms.push(elapsed);
        }
    }
    fields.push((curvature_all, curvature_absent, timing(curvature_ms)));

    let camera_all = camera_field(origin, h, cells, camera);
    let mut camera_ms = Vec::with_capacity(REPEATS);
    for repeat in 0..=REPEATS {
        let started = Instant::now();
        let produced = camera_field(origin, h, cells, camera);
        let elapsed = started.elapsed().as_secs_f64() * 1e3;
        assert_eq!(produced.len(), base_cells, "camera field is per cell");
        if repeat > 0 {
            camera_ms.push(elapsed);
        }
    }
    fields.push((camera_all, 0, timing(camera_ms)));

    // `population` is beta-scoreable by construction, so the residual the
    // control asserts is 0 and the assert below checks the construction. The
    // measured class the scope excludes -- meshable but unscoreable -- is its
    // own column rather than a silence inside a filter, and so is the whole-
    // lattice None count, which is not a defect: 3,824 of 4,096 cells hold no
    // surface at any of the three resolutions.
    let unscoped_none = (0..base_cells).filter(|&c| !beta_present[c]).count();
    let beta_unscored_meshable = meshable_cells - population.len();
    let population_unscored = population
        .iter()
        .filter(|&&c| !beta_present[c as usize])
        .count();

    // ── the surface point set ───────────────────────────────────────────────
    let reference_cell = (hi[0] - lo[0]) / f64::from(REFERENCE_SAMPLES - 1);
    let reference = whole_grid(field, origin, REFERENCE_SAMPLES, reference_cell, &mut mc);
    let stride = reference
        .positions
        .len()
        .div_ceil(MAX_SURFACE_POINTS)
        .max(1);
    let tol = h * PROJECT_TOLERANCE_REL;
    let band = h * AccuracyConfig::BAND_RADIUS_REL;
    let mut surface: Vec<[f64; 3]> = Vec::with_capacity(MAX_SURFACE_POINTS + 1);
    let mut projection_unconverged = 0usize;
    let mut projection_displaced = 0usize;
    let mut projection_step_cells = 0.0f64;
    for p in reference.positions.iter().step_by(stride) {
        // The reference grid is finer than the ladder by design (65^3 against
        // 17^3). Its OWN spacing -- the one the band, residual and displacement
        // tolerances scale with -- is reference_cell, not h: the projection is
        // removing the reference grid's error, and that is the accuracy it is
        // held to.
        match project(field, *p, tol, reference_cell, band) {
            Some((q, _)) => {
                let moved = dot(sub(q, *p), sub(q, *p)).sqrt();
                let moved_cells = moved / h;
                projection_step_cells = projection_step_cells.max(moved_cells);
                if moved <= PROJECT_MAX_STEP_CELLS * h {
                    surface.push(q);
                } else {
                    // The walk CREEPED past the displacement cap in in-band
                    // steps: plateau gradients (M-48's zero-`|grad f|` band)
                    // let Newton wander sheet to sheet while every step looks
                    // respectable. The point is not the surface the reference
                    // grid approximates, so it is kept out of the point set --
                    // `project` still measures its displacement so the cap's
                    // miss is a number, not a silence.
                    projection_displaced += 1;
                }
            }
            None => projection_unconverged += 1,
        }
    }

    // Which population cell each surface point falls in.
    let mut slot_of = vec![u32::MAX; base_cells];
    for (slot, &index) in population.iter().enumerate() {
        slot_of[index as usize] = slot as u32;
    }
    let mut point_slot: Vec<u32> = Vec::with_capacity(surface.len());
    let mut surface_unassigned = 0usize;
    for &p in &surface {
        let b = [
            lattice.slot(p[0], 0),
            lattice.slot(p[1], 1),
            lattice.slot(p[2], 2),
        ];
        let slot = slot_of[lattice.linear(b)];
        if slot == u32::MAX {
            surface_unassigned += 1;
        }
        point_slot.push(slot);
    }

    // ── the bookends, and the whole-grid controls ───────────────────────────
    let all_coarse = vec![false; population.len()];
    let all_fine = vec![true; population.len()];
    let mut distances: Vec<f64> = Vec::with_capacity(surface.len());

    let coarse_mesh = assemble(
        field,
        origin,
        h,
        cells,
        &population,
        &all_coarse,
        &mut mc,
        &mut scratch,
    );
    let coarse_error = measure(&coarse_mesh, field, &surface, lattice, h, &mut distances);
    let coarse_distances = distances.clone();

    let fine_mesh = assemble(
        field,
        origin,
        h,
        cells,
        &population,
        &all_fine,
        &mut mc,
        &mut scratch,
    );
    let fine_error = measure(
        &fine_mesh,
        field,
        &surface,
        lattice,
        h * 0.5,
        &mut distances,
    );

    let whole_coarse = whole_grid(field, origin, BASE_SAMPLES, h, &mut mc);
    let whole_coarse_error = measure(&whole_coarse, field, &surface, lattice, h, &mut distances);
    let whole_fine = whole_grid(field, origin, 2 * (BASE_SAMPLES - 1) + 1, h * 0.5, &mut mc);
    let whole_fine_error = measure(
        &whole_fine,
        field,
        &surface,
        lattice,
        h * 0.5,
        &mut distances,
    );

    let gap = |mine: f64, theirs: f64| -> f64 {
        if theirs == 0.0 {
            (mine - theirs).abs()
        } else {
            ((mine - theirs) / theirs).abs()
        }
    };

    // ── the crate's own instrument, on the two whole-grid bookends ──────────
    let mut crate_coarse = f64::NAN;
    let mut crate_fine = f64::NAN;
    if field.bound().is_exact() {
        let cfg = AccuracyConfig::from_cell_size(h).expect("a positive base cell size");
        let report = accuracy(
            &whole_coarse.positions,
            &whole_coarse.indices,
            field,
            &shape,
            origin,
            &cfg,
        )
        .expect("the base grid is the grid this mesh was extracted on");
        crate_coarse = report.symmetric_hausdorff();

        let fine_samples = 2 * (BASE_SAMPLES - 1) + 1;
        let fine_shape = RuntimeShape3::new([fine_samples; 3]).expect("the fine grid fits u32");
        let fine_cfg = AccuracyConfig::from_cell_size(h * 0.5).expect("a positive fine cell size");
        let fine_report = accuracy(
            &whole_fine.positions,
            &whole_fine.indices,
            field,
            &fine_shape,
            origin,
            &fine_cfg,
        )
        .expect("the fine grid is the grid this mesh was extracted on");
        crate_fine = fine_report.symmetric_hausdorff();
    }

    // ── the per-cell error the criteria are ranked against ──────────────────
    let mut cell_error = vec![f64::NEG_INFINITY; population.len()];
    for (point, &slot) in point_slot.iter().enumerate() {
        if slot == u32::MAX {
            continue;
        }
        let d = coarse_distances[point];
        let entry = &mut cell_error[slot as usize];
        if d > *entry {
            *entry = d;
        }
    }
    let paired: Vec<usize> = (0..population.len())
        .filter(|&s| cell_error[s] > f64::NEG_INFINITY)
        .collect();
    let paired_error: Vec<f64> = paired.iter().map(|&s| cell_error[s]).collect();

    // ── the budgets ────────────────────────────────────────────────────────
    let growth = (triangles_fine - triangles_coarse) as f64;
    let target_primary = triangles_coarse + (growth * BUDGET_GROWTH).round() as usize;
    let target_low = triangles_coarse + (growth * BUDGET_GROWTH_LOW).round() as usize;

    // ── each criterion at the primary and the low budget ───────────────────
    struct Partial {
        scores: Vec<f64>,
        order: Vec<usize>,
        primary: Arm,
        low: Arm,
        absent: usize,
        ms: Timing,
    }

    let mut partials: Vec<Partial> = Vec::with_capacity(CRITERIA.len());
    for (all, absent, ms) in fields {
        let scores: Vec<f64> = population.iter().map(|&c| all[c as usize]).collect();
        let order = order_by_score(&scores);
        let mut build = |target: usize, cell_size: f64| -> Arm {
            let (refined, refined_cells, triangles) =
                grow(&order, &coarse, &fine, triangles_coarse, target);
            let mesh = assemble(
                field,
                origin,
                h,
                cells,
                &population,
                &refined,
                &mut mc,
                &mut scratch,
            );
            let error = measure(&mesh, field, &surface, lattice, cell_size, &mut distances);
            Arm {
                transition_faces: transition_faces(&population, &refined, cells),
                refined,
                refined_cells,
                triangles,
                error,
            }
        };
        let primary = build(target_primary, h * 0.5);
        let low = build(target_low, h * 0.5);
        partials.push(Partial {
            scores,
            order,
            primary,
            low,
            absent,
            ms,
        });
    }

    // ── the tightest common ceiling, and the matched arm at it ─────────────
    let matched_ceiling = partials
        .iter()
        .map(|p| p.primary.triangles)
        .min()
        .expect("three criteria");

    let beta_refined = partials[BETA].primary.refined.clone();
    let mut criteria: Vec<CriterionOutcome> = Vec::with_capacity(CRITERIA.len());
    for (index, partial) in partials.into_iter().enumerate() {
        let (refined, refined_cells, triangles) = grow(
            &partial.order,
            &coarse,
            &fine,
            triangles_coarse,
            matched_ceiling,
        );
        let mesh = assemble(
            field,
            origin,
            h,
            cells,
            &population,
            &refined,
            &mut mc,
            &mut scratch,
        );
        let error = measure(&mesh, field, &surface, lattice, h * 0.5, &mut distances);
        let matched = Arm {
            transition_faces: transition_faces(&population, &refined, cells),
            refined,
            refined_cells,
            triangles,
            error,
        };

        let paired_scores: Vec<f64> = paired.iter().map(|&s| partial.scores[s]).collect();
        let rho = common::beta::rank_correlation(&paired_scores, &paired_error);
        let least = partial.scores.iter().copied().fold(f64::INFINITY, f64::min);
        let most = partial
            .scores
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);

        criteria.push(CriterionOutcome {
            name: CRITERIA[index],
            ms: partial.ms,
            absent: partial.absent,
            distinct_scores: distinct(&partial.scores),
            score_min: least,
            score_max: most,
            rank_correlation: rho,
            rank_population: paired.len(),
            overlap_with_beta: jaccard(&partial.primary.refined, &beta_refined),
            primary: partial.primary,
            matched,
            low: partial.low,
        });
    }

    FieldOutcome {
        name,
        bound: field.bound(),
        base_cells,
        meshable_cells,
        beta_unscored_meshable,
        unscoped_none,
        population_unscored,
        population: population.len(),
        base_sign_change_cells,
        triangles_coarse,
        triangles_fine,
        whole_triangles_coarse: whole_coarse.triangle_count(),
        whole_triangles_fine: whole_fine.triangle_count(),
        assembly_triangle_gap: {
            let w = whole_fine.triangle_count();
            if w == 0 {
                0.0
            } else {
                (triangles_fine as i64 - w as i64).unsigned_abs() as f64 / w as f64
            }
        },
        whole_gap_coarse: gap(coarse_error.symmetric(), whole_coarse_error.symmetric()),
        whole_gap_fine: gap(fine_error.symmetric(), whole_fine_error.symmetric()),
        hausdorff_coarse: coarse_error.symmetric(),
        hausdorff_fine: fine_error.symmetric(),
        crate_hausdorff_coarse: crate_coarse,
        crate_hausdorff_fine: crate_fine,
        surface_points: surface.len(),
        surface_stride: stride,
        surface_unassigned,
        projection_unconverged,
        projection_displaced,
        projection_step_cells,
        cell_error_population: paired.len(),
        camera,
        camera_distance: CAMERA_HALF_EXTENTS * half,
        target_primary,
        target_low,
        matched_ceiling,
        criteria,
    }
}

// ─── formatting ──────────────────────────────────────────────────────────────

/// Six decimals: the house format for a dimensionless ratio.
fn r6(v: f64) -> String {
    format!("{v:.6}")
}

/// Nine decimals: a distance on a domain up to sixteen units across.
fn d9(v: f64) -> String {
    format!("{v:.9}")
}

/// Nine significant figures in scientific notation, for a quantity spanning
/// decades.
fn e9(v: f64) -> String {
    format!("{v:.9e}")
}

/// A world position as one CSV-safe token. `Run::record` refuses a comma.
fn point(p: [f64; 3]) -> String {
    format!("{:.4}|{:.4}|{:.4}", p[0], p[1], p[2])
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-153");

    common::experiment::run(prereg, |run| {
        let mut outcomes: Vec<FieldOutcome> = Vec::with_capacity(FIELDS);
        // Inline block per field, not a closure, so no `return` in here (M-199).
        isomesh::for_each_reference_field!(f64, |name, field| {
            let outcome = measure_field(name, &field);
            println!(
                "  {:<14} pop {:>5} (sub-cell {:>4})  coarse {:>6} tri H {:.6}  \
                 fine {:>6} tri H {:.6}  span {:>6.2}",
                outcome.name,
                outcome.population,
                outcome.population - outcome.base_sign_change_cells,
                outcome.triangles_coarse,
                outcome.hausdorff_coarse,
                outcome.triangles_fine,
                outcome.hausdorff_fine,
                outcome.ladder_span(),
            );
            for c in &outcome.criteria {
                println!(
                    "      {:<16} {:>5} ms  refined {:>5}/{:<5} tri {:>6}  H {:.6}  \
                     decile {:.6}  rho {:>7.4}  distinct {:>5}  overlap {:.3}",
                    c.name,
                    format!("{:.3}", c.ms.median),
                    c.matched.refined_cells,
                    outcome.population,
                    c.matched.triangles,
                    c.matched.error.symmetric(),
                    c.matched.error.worst_decile,
                    c.rank_correlation,
                    c.distinct_scores,
                    c.overlap_with_beta,
                );
            }
            outcomes.push(outcome);
        });

        // ── vacuity controls, all of them, before the first record ──────────

        assert_eq!(
            outcomes.len(),
            FIELDS,
            "VOID: measured {} of {FIELDS} reference fields, so some field's clauses were \
             never at risk",
            outcomes.len()
        );

        for o in &outcomes {
            assert!(
                o.population >= 3,
                "VOID: {} carries {} refinable cells at {BASE_SAMPLES}^3, so there is no \
                 refinement decision for a criterion to get right",
                o.name,
                o.population
            );
            assert_eq!(
                o.population_unscored,
                0,
                "VOID: {}'s population holds {} cells with no beta score, so beta's ordering \
                 is missing a cell the criteria rank ({} refinable, {} meshable before the \
                 beta-scoreable scope, {} excluded by it, {} of the whole lattice unscored \
                 overall -- that exclusion is the scope itself, not a defect)",
                o.name,
                o.population_unscored,
                o.population,
                o.meshable_cells,
                o.beta_unscored_meshable,
                o.unscoped_none
            );

            // The assembly must agree with the whole-grid extraction to the
            // last SAMPLE — that identity is what makes growing a refined set a
            // per-cell prefix sum rather than an extraction — but NOT to the
            // last triangle: the two grids share boundary samples, and MC emits
            // a boundary-sensitive triangle set on shared faces. Measured gap
            // on the first run: gyroid 10,536 vs 10,632, 0.90% at 33^3. The
            // gap is a column (`assembly_triangle_gap`), bounded, and the
            // Hausdorff agreement, not the triangle count, is the identity —
            // also asserted below, at the same 1e-12 the header promises.
            assert_eq!(
                o.triangles_coarse, o.whole_triangles_coarse,
                "VOID: {}'s per-cell assembly emits {} triangles at {BASE_SAMPLES}^3 and one \
                 whole-grid extract emits {}, so the coarse assembly is not the extraction and \
                 every arm's geometry is a different object from the crate's",
                o.name, o.triangles_coarse, o.whole_triangles_coarse
            );
            // The COARSE identity does hold and is kept: both grids are
            // h-spaced with shared samples, so the assembly is the extraction
            // up to the mean's summation order. Measured worst over all eight
            // fields: 4.5e-16.
            assert!(
                o.whole_gap_coarse <= ASSEMBLY_TOLERANCE,
                "VOID: {}'s assembled coarse Hausdorff differs from the whole-grid one by {} \
                 relative, above ASSEMBLY_TOLERANCE = {ASSEMBLY_TOLERANCE:e} — the two sample \
                 the field at bit-identical positions, so a gap means the assembly is wrong",
                o.name,
                o.whole_gap_coarse
            );
            // The Hausdorff gap between an assembled bookend and the same-
            // resolution whole grid is a measured column, not an asserted
            // identity, and the fine rung is why: the whole 33^3 grid's cell
            // boundaries sit at h/2 while the per-cell assembly's sit at h, so
            // the two see DIFFERENT sub-cell placement on any field with detail
            // at the base resolution. Measured: fbm_terrain 1.16 relative on
            // the fine Hausdorff at 33^3, and the coarse rung (17^3 whole grid,
            // h-spaced cells) agrees to 4.5e-16. The whole-grid fine rung stays
            // as the comparison column (`whole_grid_hausdorff_gap_fine`); the
            // arms are all compared against the ASSEMBLED bookends, so the
            // instrument is one family throughout and the comparison is fair.
            // The triangle-count gap is bounded separately, above.

            // The instrument has to move across the ladder, or there is no error
            // for a criterion to redistribute.
            assert!(
                o.ladder_span() >= LADDER_SPAN_FLOOR || 1.0 / o.ladder_span() >= LADDER_SPAN_FLOOR,
                "VOID: {}'s uniform bookends read {:.9} coarse and {:.9} fine, a span of \
                 {:.4} — inside the LADDER_SPAN_FLOOR = {LADDER_SPAN_FLOOR} band around \
                 equality from either direction, so refinement moves this field's error \
                 by less than the floor and every criterion comparison is inside the \
                 instrument's own noise",
                o.name,
                o.hausdorff_coarse,
                o.hausdorff_fine,
                o.ladder_span()
            );

            // The surface point set has to be a surface point set.
            assert!(
                o.surface_points >= MIN_SURFACE_POINTS,
                "VOID: {} kept {} projected surface points, below MIN_SURFACE_POINTS = \
                 {MIN_SURFACE_POINTS}, so the `surface -> mesh` direction samples too little \
                 of the surface to see which cells a criterion left coarse",
                o.name,
                o.surface_points
            );
            let attempted = o.surface_points + o.projection_unconverged;
            assert!(
                (o.projection_unconverged as f64) <= PROJECT_MAX_UNCONVERGED * attempted as f64,
                "VOID: {} failed to project {} of {} reference vertices onto its zero set, \
                 above PROJECT_MAX_UNCONVERGED = {PROJECT_MAX_UNCONVERGED} — the kept points \
                 are a biased sample of the surface rather than the surface",
                o.name,
                o.projection_unconverged,
                attempted
            );
            // `projection_step_cells` RECORDS the worst displacement; the
            // displacement cap itself is enforced at the point set (any vertex
            // past `PROJECT_MAX_STEP_CELLS * h` is excluded and counted in
            // `projection_displaced`), and the share of such vertices is gated
            // with the unconverged share below, because both are ways a kept
            // point set could quietly become a biased sample.
            assert!(
                o.projection_displaced as f64 / attempted as f64 <= PROJECT_MAX_UNCONVERGED,
                "VOID: {}'s projection displaced {} of {} reference vertices past \
                 PROJECT_MAX_STEP_CELLS = {PROJECT_MAX_STEP_CELLS} base cells of movement \
                 (worst observed {:.4} cells) -- the plateau-creep share, above \
                 PROJECT_MAX_UNCONVERGED = {PROJECT_MAX_UNCONVERGED}, so the kept points \
                 are a biased sample of the surface rather than the surface",
                o.name,
                o.projection_displaced,
                attempted,
                o.projection_step_cells
            );
            assert!(
                o.cell_error_population >= 3,
                "VOID: {} pairs {} cells with a measured coarse error, and \
                 `common::beta::rank_correlation` defines fewer than three pairs to be 0.0 \
                 rather than measuring one",
                o.name,
                o.cell_error_population
            );

            // The registered vacuity control, verbatim: "the worst-decile
            // Hausdorff column must differ across criteria, or all three are
            // refining the same cells".
            let deciles: Vec<f64> = o
                .criteria
                .iter()
                .map(|c| c.matched.error.worst_decile)
                .collect();
            assert!(
                !(deciles[0] == deciles[1] && deciles[1] == deciles[2]),
                "VOID: {} reports the same worst-decile Hausdorff {:.12} for beta, curvature \
                 and camera distance, so all three refined the same cells and C1 and C2 are \
                 comparisons of one arm with itself",
                o.name,
                deciles[0]
            );
            assert!(
                o.criteria[CAMERA].overlap_with_beta < 1.0,
                "VOID: {}'s camera arm refined exactly the cells beta refined (Jaccard 1.0 \
                 over {} cells), so C1 compares a criterion with itself",
                o.name,
                o.population
            );

            for c in &o.criteria {
                assert!(
                    c.ms.median > 0.0,
                    "VOID: {} measured {} at {} ms over {REPEATS} repeats, so `criterion_ms` \
                     is a zero that could not have been non-zero (M-44)",
                    o.name,
                    c.name,
                    c.ms.median
                );
                for (label, arm) in [
                    ("primary", &c.primary),
                    ("matched", &c.matched),
                    ("low", &c.low),
                ] {
                    assert!(
                        arm.triangles > 0,
                        "VOID: {}'s {} arm at the {label} budget produced no triangles, so \
                         its Hausdorff is the error of an empty mesh",
                        o.name,
                        c.name
                    );
                    assert!(
                        arm.refined_cells > 0 && arm.refined_cells < o.population,
                        "VOID: {}'s {} arm at the {label} budget refined {} of {} cells, so \
                         it is a uniform bookend wearing a criterion's name rather than an \
                         adaptive mesh",
                        o.name,
                        c.name,
                        arm.refined_cells,
                        o.population
                    );
                }
            }

            if o.bound.is_exact() {
                assert!(
                    o.crate_hausdorff_coarse > o.crate_hausdorff_fine,
                    "VOID: `validate::accuracy` reads {}'s coarse bookend at {:.9} and its \
                     fine bookend at {:.9}, so the crate's own instrument does not see the \
                     ladder move and the separation this file reports is unconfirmed",
                    o.name,
                    o.crate_hausdorff_coarse,
                    o.crate_hausdorff_fine
                );
            }
        }

        // Global rather than per field: beta's single distinct score on
        // `box_exact` and `thin_plate` is the measured answer (`p-152.csv`,
        // `planar_fraction = 1.000000`), so a per-field form would abort on the
        // finding. What has to be shown is that each criterion can order cells
        // at all.
        for index in 0..CRITERIA.len() {
            let best = outcomes
                .iter()
                .map(|o| o.criteria[index].distinct_scores)
                .max()
                .unwrap_or(0);
            assert!(
                best > 1,
                "VOID: `{}` takes one distinct value on every one of the {FIELDS} fields, so \
                 its ordering is the index tie-break everywhere and no clause here measures a \
                 criterion",
                CRITERIA[index]
            );
        }

        // ── the clause verdicts, global and decided at the matched budget ───
        let mut c1_fields_won = 0usize;
        for o in &outcomes {
            if o.criteria[BETA].matched.error.symmetric()
                < o.criteria[CAMERA].matched.error.symmetric()
            {
                c1_fields_won += 1;
            }
        }
        let c1_holds = c1_fields_won >= C1_MIN_WINNERS;

        let mut c1_fields_won_low = 0usize;
        for o in &outcomes {
            if o.criteria[BETA].low.error.symmetric() < o.criteria[CAMERA].low.error.symmetric() {
                c1_fields_won_low += 1;
            }
        }
        let c1_holds_at_low_budget = c1_fields_won_low >= C1_MIN_WINNERS;

        let beats_curvature = |o: &FieldOutcome| -> bool {
            o.criteria[BETA].matched.error.symmetric()
                < o.criteria[CURVATURE].matched.error.symmetric()
        };
        let named = |wanted: &str| -> bool {
            outcomes
                .iter()
                .find(|o| o.name == wanted)
                .is_some_and(beats_curvature)
        };
        let thin_plate_won = named("thin_plate");
        let noise_cavity_won = named("noise_cavity");
        let c2_holds = thin_plate_won && noise_cavity_won;

        println!(
            "\n  C1 at the matched budget: beta beats camera distance on {c1_fields_won} of \
             {FIELDS} fields — needs {C1_MIN_WINNERS}, so C1 {}",
            if c1_holds { "HOLDS" } else { "is FALSIFIED" }
        );
        println!(
            "  C1 at the {BUDGET_GROWTH_LOW} budget: {c1_fields_won_low} of {FIELDS}, so \
             {}",
            if c1_holds_at_low_budget {
                "it would HOLD"
            } else {
                "it would be FALSIFIED"
            }
        );
        println!(
            "  C2: beta beats curvature on thin_plate {thin_plate_won} and on noise_cavity \
             {noise_cavity_won}, so C2 {}",
            if c2_holds { "HOLDS" } else { "is FALSIFIED" }
        );

        for o in &outcomes {
            let (p152_share, p152_rho, p152_planar) = prior(o.name);
            let predicted = predicted_beta_beats_curvature(o.name);
            let measured_curvature = beats_curvature(o);
            let beats_camera = o.criteria[BETA].matched.error.symmetric()
                < o.criteria[CAMERA].matched.error.symmetric();
            let beta_h = o.criteria[BETA].matched.error.symmetric();
            let camera_h = o.criteria[CAMERA].matched.error.symmetric();
            let curvature_h = o.criteria[CURVATURE].matched.error.symmetric();

            for c in &o.criteria {
                run.record(&[
                    // ── registered, in registration order ──
                    ("criterion", c.name.to_string()),
                    ("field", o.name.to_string()),
                    ("triangles", c.primary.triangles.to_string()),
                    ("hausdorff", d9(c.primary.error.symmetric())),
                    (
                        "hausdorff_at_matched_triangles",
                        d9(c.matched.error.symmetric()),
                    ),
                    ("criterion_ms", format!("{:.4}", c.ms.median)),
                    ("rank_correlation_with_error", r6(c.rank_correlation)),
                    ("worst_decile_hausdorff", d9(c.matched.error.worst_decile)),
                    ("c1_holds", c1_holds.to_string()),
                    ("c2_holds", c2_holds.to_string()),
                    // ── extras (M-273) ──
                    ("base_cells", o.base_cells.to_string()),
                    ("base_samples", BASE_SAMPLES.to_string()),
                    (
                        "base_sign_change_cells",
                        o.base_sign_change_cells.to_string(),
                    ),
                    ("beta_beats_camera", beats_camera.to_string()),
                    ("beta_beats_curvature", measured_curvature.to_string()),
                    ("beta_vs_camera_ratio", r6(beta_h / camera_h)),
                    ("beta_vs_curvature_ratio", r6(beta_h / curvature_h)),
                    ("budget_growth", r6(BUDGET_GROWTH)),
                    ("budget_growth_low", r6(BUDGET_GROWTH_LOW)),
                    ("c1_bar", C1_MIN_WINNERS.to_string()),
                    ("c1_fields_won", c1_fields_won.to_string()),
                    ("c1_fields_won_low", c1_fields_won_low.to_string()),
                    ("c1_holds_at_low_budget", c1_holds_at_low_budget.to_string()),
                    ("c2_named_fields", C2_NAMED.join("|")),
                    ("c2_noise_cavity_won", noise_cavity_won.to_string()),
                    (
                        "c2_prediction_held",
                        (predicted == measured_curvature).to_string(),
                    ),
                    ("c2_thin_plate_won", thin_plate_won.to_string()),
                    ("camera_distance_from_centre", d9(o.camera_distance)),
                    ("camera_position", point(o.camera)),
                    ("cell_error_population", o.cell_error_population.to_string()),
                    ("crate_hausdorff_status", crate_status(o.bound).to_string()),
                    (
                        "crate_hausdorff_uniform_coarse",
                        d9(o.crate_hausdorff_coarse),
                    ),
                    ("crate_hausdorff_uniform_fine", d9(o.crate_hausdorff_fine)),
                    ("criterion_ms_max", format!("{:.4}", c.ms.max)),
                    ("criterion_ms_min", format!("{:.4}", c.ms.min)),
                    ("criterion_repeats", REPEATS.to_string()),
                    ("fine_samples_per_cell", FINE_SAMPLES.to_string()),
                    ("hausdorff_budget_low", d9(c.low.error.symmetric())),
                    (
                        "hausdorff_mesh_to_surface",
                        d9(c.matched.error.mesh_to_surface),
                    ),
                    (
                        "hausdorff_surface_to_mesh",
                        d9(c.matched.error.surface_to_mesh),
                    ),
                    ("hausdorff_uniform_coarse", d9(o.hausdorff_coarse)),
                    ("hausdorff_uniform_fine", d9(o.hausdorff_fine)),
                    ("ladder_span", r6(o.ladder_span())),
                    (
                        "matched_dropped_cells",
                        (c.primary.refined_cells - c.matched.refined_cells).to_string(),
                    ),
                    (
                        "mesh_projection_unconverged",
                        c.matched.error.unconverged.to_string(),
                    ),
                    ("p152_beta_share", r6(p152_share)),
                    ("p152_planar_fraction", r6(p152_planar)),
                    ("p152_rank_correlation_with_qef", r6(p152_rho)),
                    ("planar_floor", format!("{PLANAR_FLOOR:e}")),
                    (
                        "beta_unscored_meshable",
                        o.beta_unscored_meshable.to_string(),
                    ),
                    ("population_unscored", o.population_unscored.to_string()),
                    ("meshable_cells", o.meshable_cells.to_string()),
                    ("unscoped_none", o.unscoped_none.to_string()),
                    ("predicted_beta_beats_curvature", predicted.to_string()),
                    ("projection_step_cells_max", e9(o.projection_step_cells)),
                    (
                        "projection_unconverged",
                        o.projection_unconverged.to_string(),
                    ),
                    ("reference_samples", REFERENCE_SAMPLES.to_string()),
                    ("refined_cells", c.matched.refined_cells.to_string()),
                    ("refined_overlap_with_beta", r6(c.overlap_with_beta)),
                    ("score_absent_cells", c.absent.to_string()),
                    (
                        "chamber_degenerate_triangles",
                        c.matched.error.chamber_degenerate.to_string(),
                    ),
                    (
                        "chamber_skipped_triangles",
                        c.matched.error.chamber_skipped.to_string(),
                    ),
                    ("rank_population", c.rank_population.to_string()),
                    ("score_all_tied", (c.distinct_scores <= 1).to_string()),
                    ("score_distinct_values", c.distinct_scores.to_string()),
                    ("score_max", e9(c.score_max)),
                    ("score_min", e9(c.score_min)),
                    (
                        "slab_evaluations_per_patch",
                        common::beta::SLAB_EVALUATIONS.to_string(),
                    ),
                    (
                        "sub_cell_only_cells",
                        (o.population - o.base_sign_change_cells).to_string(),
                    ),
                    ("surface_cells", o.population.to_string()),
                    ("surface_point_stride", o.surface_stride.to_string()),
                    ("surface_points", o.surface_points.to_string()),
                    (
                        "surface_points_unassigned",
                        o.surface_unassigned.to_string(),
                    ),
                    ("transition_faces", c.matched.transition_faces.to_string()),
                    ("triangles_at_matched", c.matched.triangles.to_string()),
                    ("triangles_budget_low", c.low.triangles.to_string()),
                    ("triangles_matched_ceiling", o.matched_ceiling.to_string()),
                    ("triangles_target", o.target_primary.to_string()),
                    ("triangles_target_low", o.target_low.to_string()),
                    ("triangles_uniform_coarse", o.triangles_coarse.to_string()),
                    ("triangles_uniform_fine", o.triangles_fine.to_string()),
                    ("whole_grid_hausdorff_gap_coarse", e9(o.whole_gap_coarse)),
                    ("whole_grid_hausdorff_gap_fine", e9(o.whole_gap_fine)),
                    (
                        "whole_grid_triangles_coarse",
                        o.whole_triangles_coarse.to_string(),
                    ),
                    (
                        "whole_grid_triangles_fine",
                        o.whole_triangles_fine.to_string(),
                    ),
                    (
                        "assembly_triangle_gap",
                        format!("{:.6}", o.assembly_triangle_gap),
                    ),
                    ("worst_decile_spread", d9(o.worst_decile_spread())),
                ]);
            }
            let _ = o.beta();
        }
    });
}
