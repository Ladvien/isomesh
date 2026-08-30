//! **P-173 — the curl of the Hermite data is the inconsistency `lambda = 0.01` is guessing about.**
//!
//! Ticket: R-173. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p173
//! ```
//!
//! Writes `docs/experiments/p-173.csv`, one row per (field, resolution) over all
//! eight reference fields at 17³/25³/33³ — **M-60's own three grids**, so C2's
//! population is the population M-60 counted.
//!
//! # What was missing
//!
//! The registration opens with a claim about the corpus, and it is checkable, so
//! it was checked before the harness was written. `Hodge`, `Helmholtz`,
//! `integrability` and `irrotational` appear in **exactly one file** in this
//! repository — `crates/isomesh/src/experiment.rs`, i.e. inside the P-173
//! registration itself. `curl` is stronger still: every occurrence in `src/` is
//! `crates/isomesh/src/marching_cubes/reference.rs:16`, which is the **HTTP
//! tool** ("the page was fetched with `curl`"), and every occurrence in
//! `docs/research/` is the surname **Curless** (Curless & Levoy 1996). The
//! differential operator is not in this repository at all. So the mechanism this
//! row proposes has no prior art here, not even as flow visualisation.
//!
//! What *is* here is the thing the residual is about. `dual_contouring/solve.rs`
//! places a cell's vertex at `x = c + adj(M + lambda*I)*g / det(M + lambda*I)`
//! with `M = sum n_i n_i^T` over the cell's crossings and `lambda = 0.01`
//! (`solve.rs:85`). `solve.rs:43-54` says what `lambda` is for — *"the Tikhonov
//! regularizer that stops an under-determined cell … from flying off"* — and
//! `M-107` swept it over six decades and found it moves the runaway 23×, but
//! only on the fields `M-30` already named. Nowhere does anything ask the prior
//! question: **are those normals even consistent with a single smooth sheet?**
//! `lambda` is one global constant standing in for a per-cell property nobody had
//! measured.
//!
//! `M-60` is the other half. It counted, per field and grid, the cells where
//! `ManifoldDualContouring` emits **more than one vertex** — the cells carrying
//! more than one sheet — as extra vertices over plain `DualContouring`: `gyroid`
//! **3.13% / 2.05% / 0.53%** and `fbm_terrain` **1.70% / 0.84% / 0.77%** at
//! 17³/25³/33³, and `sphere`, `torus`, `box_exact`, `csg_difference`,
//! `thin_plate` **exactly 0 at every resolution**. `M-61` then showed those two
//! fields are also the two whose self-intersections get *worse* when the vertex
//! is split. So multi-sheet-ness is a measured, localised, field-dependent
//! property, and this row asks whether a curl residual predicts it.
//!
//! `noise_cavity` is the eighth field and **M-60 never saw it** (M-60 says
//! "seven fields"). Its second-vertex count here is a new number, not a
//! reproduction, and it is treated as such: the registration's control names
//! *five* fields and the control asserts against exactly those five.
//!
//! # The curl residual, defined exactly
//!
//! A cell has twelve edges. The **unit** gradient is sampled at all twelve edge
//! **midpoints**, giving `n_e = grad f / ||grad f||` at `m_e`, `e in 0..12`. Each
//! of the cube's six faces `(w, side)` is a closed square loop of side `h`
//! through four of those samples. With `u = (w+1) % 3` and `v = (w+2) % 3` — the
//! crate's own right-handed `(u, v, w)` convention from `cube.rs:81-87` — the
//! loop traversed counter-clockwise as seen from `+w` runs
//! `(0,0) -> (1,0) -> (1,1) -> (0,1) -> (0,0)`, so its four tangents are `+u`,
//! `+v`, `-u`, `-v`. The circulation of a gradient field around a closed loop is
//! exactly zero, so the residual is the amount by which it is not:
//!
//! ```text
//! C(w, side) = n_bottom[u] - n_top[u] + n_right[v] - n_left[v]
//! ```
//!
//! where *bottom*/*top* are the loop's two `u`-edges at `v`-offset 0 and 1, and
//! *left*/*right* its two `v`-edges at `u`-offset 0 and 1. Each edge is placed by
//! `EDGE_AXIS` and the offsets of `EDGE_CORNERS[e][0]`; membership is the crate's
//! `edge_on_face`. The common factor `h` — every edge has length `h` — is
//! divided out, so `C` is the circulation **in units of the edge length** and is
//! dimensionless, because the normals are unit.
//!
//! The cell's residual is the Euclidean norm of the six-vector:
//!
//! ```text
//! r_cell = || [C(0,0), C(0,1), C(1,0), C(1,1), C(2,0), C(2,1)] ||_2
//! ```
//!
//! **Six faces and not the three axis-aligned components**, deliberately. The
//! three-component form averages the two faces normal to each axis, and a second
//! sheet entering a cell through one face and not the opposite one is precisely
//! the signature that averaging destroys.
//!
//! ## Two normalisations of the one residual, and which reads which
//!
//! `r_cell` is bounded: each term of `C` is a component of a unit vector, so
//! `|C| <= 4` and `r_cell <= 4*sqrt(6) = 9.797958971132712`. Dividing by that
//! **extremal** bound gives `curl_residual_normalised in [0, 1]`, which is the
//! reported score and the score the AUC ranks. The bound is asserted against
//! `4*sqrt(6)` before any row is written.
//!
//! One consequence of it bounds the vacuity control's ceiling: a *single* face
//! carrying a full normal reversal on opposite edges contributes `|C| = 4` by
//! itself, so `r_cell >= 4` and the normalised score is at least
//! `1/sqrt(6) = 0.408248`. That is `REVERSAL_FLOOR`, and `CONTROL_CEILING = 0.25`
//! sits below it so that a field certified "near zero" provably has no cell whose
//! mean is carried by a full reversal on any face. `REVERSAL_FLOOR` is
//! **reported** (`reversal_floor_cells`) and not asserted: reaching it needs two
//! normals *exactly* axis-aligned and exactly opposed, which no smooth field
//! produces, so asserting it would fail on a fixture that is measuring perfectly
//! well.
//!
//! The per-cell `lambda` map reads a **different** scale, and this is stated
//! plainly because the choice was made after seeing the first one degenerate. The
//! extremal bound is a worst case over configurations no reference field reaches,
//! so every real cell lands in its bottom decade and a map anchored on it becomes
//! a **global decrease** of `lambda` — which cannot answer C3, because C3 asks
//! whether a *per-cell* regularizer beats a global one, and a map that is
//! constant in practice is a global one. So the map reads
//! `s = min(|| C ||_2 / UNIT_CIRCULATION, 1)` with `UNIT_CIRCULATION = 1`: one
//! whole unit vector's worth of failure to cancel around one loop. That is a
//! property of the construction and not of the data — the terms *are* components
//! of unit vectors, so one unit of disagreement is the quantity's own unit — and
//! it puts the crossover inside the measured range instead of above it.
//!
//! The AUC is unaffected by the choice: it is computed on ranks, and the two
//! scales differ by a positive constant factor, which no ranking can see. Both
//! scores are recorded (`curl_residual_normalised` and `lambda_score_mean`).
//!
//! ## Why all twelve midpoints, and not the crossings the crate already has
//!
//! Because a partial loop is not a loop. A face's circulation needs all four of
//! its edges, and a Marching Cubes case cuts three to six edges of the twelve —
//! `complete_faces_mean` is recorded per row and is the measured answer to how
//! often a face has all four of its edges cut. Sampling the midpoints instead is
//! one path with no case analysis: every cell has twelve normals and every face
//! is a genuine closed loop.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `field` × `resolution` | eight reference fields at 17³/25³/33³ | no |
//! | `sphere`, `torus`, `box_exact`, `csg_difference`, `thin_plate` | M-60's five zero-second-vertex fields, the registration's named negatives-only control | **yes** |
//! | `curl_ms` | the residual arithmetic over pre-built normals | no |
//! | `curl_ms_sampled` | the same, paying for its own twelve gradients | **yes** (the cost the other arm excludes) |
//! | `lambda_global` | the shipped `solve::LAMBDA` through the shipped solve | **yes** (C3's baseline) |
//! | `lambda_per_cell` | `0.001 * 100^r` through the same solve | no |
//!
//! # C1, and which of two readings carries the verdict
//!
//! The registration says the residual is *"computable in a few dozen flops on a
//! struct the crate already builds"*. So `curl_ms` is the **residual
//! arithmetic**, timed over normals already in memory — the flops, and nothing
//! else — and `curl_share = curl_ms / extract_ms` against
//! `DualContouring::extract` on the same grid. That is the registered reading and
//! it carries `c1_holds`.
//!
//! The sampling is not free and pretending otherwise would be dishonest, so it
//! is measured and recorded beside it: `curl_ms_sampled` pays for its own twelve
//! `Sdf::gradient` calls per active cell (six samples each, `sdf.rs:78-95`) and
//! `curl_share_sampled` / `c1_holds_sampled` carry that verdict. Arithmetic: an
//! extraction at `n³` costs `n³` samples plus one gradient per crossing on each
//! active cell; the residual as constructed costs twelve gradients — 72 samples —
//! on each active cell, so `curl_share_sampled` is expected to be **greater than
//! one**, not 0.02. A consumer who wants this diagnostic at zero marginal
//! sampling cost needs a construction on the crossings themselves, which is a
//! different residual and is not what was registered. Both numbers are in the
//! CSV; neither is hidden behind the other.
//!
//! Timing is `std::time::Instant`, five measured repeats after one warm-up, and
//! the **median** is the headline with min and max recorded — this host's
//! `amd-pstate-epp` governor swings the same binary 1.45× between runs (`M-280`,
//! `x24`).
//!
//! # C3 is driven through the shipped solve, and this is the construction
//!
//! `DualContouring::set_lambda` (`dual_contouring.rs:223`) is public and
//! **global**: it sets one `lambda` for the whole grid, so a per-cell `lambda`
//! cannot be driven through the shipped extractor. It does not follow that the
//! solve has to be rewritten. `dual_contouring::solve::solve_with(cell, lambda)`
//! is **public** (`solve.rs:272`) and takes the regularizer explicitly, and
//! `HermiteCell::from_corners` (`hermite.rs:74`) is public too — which are
//! exactly the two calls `Qef::place` makes (`dual_contouring.rs:163-168`). So
//! the per-cell arm is the *shipped* solve called once per cell with a varying
//! `lambda`, not a second copy of the normal equations. That matters for the
//! measurement, not just for tidiness: a re-derived solve would make
//! `sharpness_delta` the sum of "the `lambda` moved" and "the arithmetic
//! differs", and there would be no way to attribute it.
//!
//! Only one piece of `Qef::place` is not public — `apply_clamp` under
//! `Clamp::ToCell`, which is `DualContouring::new`'s default. Its twelve lines
//! are copied verbatim (`dual_contouring.rs:184-203`), the way `experiment_p121`
//! copies them, and **the copy is verified rather than trusted**: the global-λ
//! arm is compared against `DualContouring`'s own output **bit for bit** on every
//! row, and `mirror_mismatches` is a column that is asserted to be zero. Without
//! that, `sharpness_delta`, `qef_residual_delta` and `self_intersections_delta`
//! would be differences against something that is not the shipped path.
//!
//! Both arms share one index buffer, and that is exact rather than convenient: a
//! dual method's topology is one quad per sign-changing grid edge
//! (`dual.rs:648-655`) and depends only on the corner signs, so moving a vertex
//! cannot change a single index. `place_vertices` walks active cells in
//! lexicographic `(z, y, x)` order and pushes one vertex per cell for `Qef`
//! (`dual.rs:487-497`, and `emit_vertices` at `:633` emits all of them in that
//! order with no culling), so vertex `i` is active cell `i` and the substitution
//! is a positions-only swap.
//!
//! ## The per-cell `lambda`, fixed before the run
//!
//! ```text
//! s = min(|| C ||_2 / UNIT_CIRCULATION, 1)          UNIT_CIRCULATION = 1
//! lambda(s) = LAMBDA_MIN * (LAMBDA_MAX / LAMBDA_MIN)^s  =  0.001 * 100^s
//! ```
//!
//! The bracket is `[0.001, 0.1]`: the top is the constant Subgrid MT's QEF
//! already uses, one of the three the corpus circulates (`solve.rs:56-57`), and
//! the bottom is one decade below the shipped default. The geometric mean of the
//! bracket is `sqrt(0.001 * 0.1) = 0.01`, which is **exactly `solve::LAMBDA`**,
//! so the map is "the shipped constant, moved one decade down where the Hermite
//! data is consistent and one decade up where it is not", and `s = 1/2` — half a
//! unit of circulation — is the crossover. `over_shipped_lambda_share` records
//! what fraction of cells landed above it, and `lambda_cell_min` /
//! `lambda_cell_max` bound what the map actually reached.
//!
//! Both the bracket and the scale are fixed before the run and neither is fitted
//! to a measured residual — the scale's derivation is in *Two normalisations*
//! above and is a property of the construction. A bracket or an anchor chosen to
//! make C3 come out would make C3 a fit rather than a prediction, which is the
//! one thing this apparatus exists to prevent.
//!
//! ## The two axes C3 is scored on
//!
//! - **sharpness** = mean over active cells of `|f(x)| / h` at the placed vertex.
//!   **Lower is sharper.** This is `M-107`'s own instrument — it reports *"worst
//!   `|f|/h`"* for exactly this sweep — taken as a mean over the population
//!   rather than a worst, with the worst recorded beside it. It has the right
//!   direction by construction: a `lambda` toward zero reproduces a three-plane
//!   corner, which lies *on* the surface, and a large `lambda` pulls the vertex
//!   to the crossings' centroid, which at a convex corner is inside the solid.
//!   `Clamp::ToCell` bounds it, so a runaway cannot dominate the mean.
//! - **self-intersections** = `validate::self_intersections(...).count()` on the
//!   substituted positions with the shared index buffer. `M-28` and `M-61` both
//!   score dual contouring on exactly this.
//!
//! `qef_residual` is the objective itself, `E(x) = sum (n_i . (x - p_i))^2` over
//! the cell's crossings, divided by `h^2` so it is dimensionless and comparable
//! across resolutions, meaned over active cells. Lower is better. Every `_delta`
//! is **per-cell minus global**, so a negative delta is the per-cell `lambda`
//! winning.
//!
//! `c3_holds` is "beats on at least one without losing on the other": a strict
//! win on sharpness with no loss on self-intersections, or the reverse. An exact
//! tie is neither a win nor a loss.
//!
//! # C2, and how the AUC is computed
//!
//! The label is the crate's own second-vertex predicate, not a proxy.
//! `CycleQef::place` (`manifold_dual_contouring.rs:218-258`) takes
//! `segment_links(case, joined_mask(corner, ambiguous))` with `ambiguous = 0`
//! under `FaceAmbiguity::Separate` — `ManifoldDualContouring::new`'s default —
//! walks the disjoint cycles of the cut edges, and calls `push_component` once
//! per cycle, i.e. emits one vertex per cycle (`dual.rs:107-118`). All of
//! `segment_links`, `joined_mask` and `NO_EDGE` are public, so the cycle count is
//! computed here with the same three calls. A cell is **positive** when it has
//! two or more cycles.
//!
//! That is checked against the live extractor rather than asserted from the
//! paper: `ManifoldDualContouring` and `DualContouring` are both run on every
//! grid, and `mdc_vertices - dc_vertices == sum over active cells of
//! (cycles - 1)` must hold as an **integer identity**, with
//! `dc_vertices == cells`. This is M-60's own definition of "needs a second
//! vertex" — extra vertices over plain Dual Contouring — reproduced exactly.
//!
//! The score is `curl_residual_normalised`. The AUC is the Mann-Whitney
//! statistic on **ranks, with tied scores taking their average rank**:
//! `AUC = (R+ - n+(n+ + 1)/2) / (n+ * n-)` where `R+` is the rank sum of the
//! positives. Ranking rather than thresholding is what makes it a scale-free
//! statement about ordering, and average ranks are what keep a plateau of equal
//! residuals — `box_exact`'s flat cells give `r = 0` exactly — from inflating it.
//! The sort is `f64::total_cmp`, never `partial_cmp().unwrap()`.
//!
//! A row with no positive has no AUC. `separation_auc` reads `unreachable` on
//! those rows and `separation_auc_reachable` is `false`; the arithmetic is
//! `n+ = 0`, so the denominator `n+ * n-` is zero. That is the registration's own
//! *"a population with no negatives"* seen from the other side, and it is
//! recorded rather than dropped.
//!
//! `c2_holds` is per row and means **"this row supports C2"**: on a row with an
//! AUC, that the AUC exceeds 0.8; on a zero-second-vertex row, that the row
//! discharges its half of the control by reporting a near-zero residual
//! (`curl_residual_normalised < 0.25`). `c1_holds` and `c3_holds` are per row
//! throughout. The three `c*_holds_global` extras carry unanimity across all 24
//! rows on every row, with `c*_rows_held` giving the count, so a partial result
//! is legible instead of being rounded to a verdict.
//!
//! # SHARE, recomputed before the numbers
//!
//! *"C1 moves the vertex-placement stage; C3 moves the same stage's quality, not
//! its cost."* Both halves are about `place_vertices` and neither is about
//! topology, and the harness is built so that is visibly true: the index buffer
//! is shared byte-for-byte between the two arms, so nothing this row proposes can
//! move a triangle. `P-121` measured the stage this lands in — the same file's
//! `place` stage — so the share `curl_share` should be read against that rather
//! than against the whole pipeline. C1's denominator here is whole-extraction
//! wall clock, which is the conservative choice: it is the largest honest
//! denominator, so a share under 2% of it is under 2% of any stage inside it.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record` and every panic starts `VOID: `.
//!
//! - **The normalising bound is the bound.** `RESIDUAL_SCALE` must equal
//!   `4*sqrt(6)`, or every `curl_residual_normalised` is divided by the wrong
//!   number and `REVERSAL_FLOOR` and `CONTROL_CEILING` are both meaningless.
//! - **The mirror is the shipped path.** `mirror_vertices == dc_vertices` and
//!   `mirror_mismatches == 0` on every row, compared as IEEE bit patterns, or the
//!   three `_delta` columns are not measured against `DualContouring`.
//! - **The cycle count is the crate's.** `extra_vertices_measured` (from the two
//!   live extractors) must equal `extra_vertices_predicted` (from
//!   `segment_links`) exactly, and `dc_vertices == cells`, or
//!   `second_vertex_cells` is not the population `M-60` counted and C2 is scored
//!   on the wrong label.
//! - **M-60's five zero fields are still zero.** `second_vertex_cells == 0` on
//!   `sphere`, `torus`, `box_exact`, `csg_difference` and `thin_plate` at all
//!   three resolutions, which is `M-60`'s *"exactly 0 at every resolution"*
//!   re-measured. If it fails, the registration's control names a population that
//!   does not exist.
//! - **Those five report a near-zero residual.** `curl_residual_normalised < 0.25`
//!   on each of those fifteen rows — the registration's own control, with the
//!   ceiling set below `REVERSAL_FLOOR = 0.408248` so that passing it *means* no
//!   cell in the mean carries a full normal reversal on any face.
//! - **M-60's two non-zero fields still have positives, and negatives.**
//!   `second_vertex_cells > 0` and `cells - second_vertex_cells > 0` on `gyroid`
//!   and `fbm_terrain` at all three resolutions, or C2's AUC is unreachable at
//!   exactly the two fields C2 names.
//! - **The instrument fires.** Some cell somewhere must exceed `CONTROL_CEILING`
//!   — the very ceiling the previous control certifies the five zero fields
//!   *below*, so that certificate says something — or every residual in this CSV
//!   is a zero that could not have been non-zero (`M-44`). Scored on
//!   `CONTROL_CEILING` and not `REVERSAL_FLOOR`, for the reason given under *Two
//!   normalisations*: the reversal floor is an extremum no smooth field reaches.
//! - **The per-cell `lambda` is not the global one.** `lambda_cell_max /
//!   lambda_cell_min > 2` on some row, or C3 compares `0.01` against `0.01`.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]

mod common;

use std::hint::black_box;
use std::time::Instant;

use isomesh::dual_contouring::{CLAMP_EPSILON, DualContouring, solve};
use isomesh::fields::ReferenceField;
use isomesh::hermite::HermiteCell;
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::table::{
    EDGE_AXIS, EDGE_CORNERS, EDGE_COUNT, FACE_COUNT, NO_EDGE, edge_on_face, is_inside,
    segment_links,
};
use isomesh::validate::self_intersections;
use isomesh::{MeshBuffer, Sdf, Shape3};

/// One scalar type throughout. `f64` because the QEF forms `M = A^T A` and
/// squares the condition number, which is the reason `Real` spans `f64` at all
/// (`solve.rs:59-71`, `V-18`).
type Scalar = f64;

// ─── the fixture ───────────────────────────────────────────────────────────

/// Samples per axis. **M-60's own three grids** (`17³ / 25³ / 33³`), so
/// `second_vertex_cells` is comparable to the percentages M-60 reports rather
/// than merely similar to them.
const RESOLUTIONS: [u32; 3] = [17, 25, 33];

/// Measured repeats per timed quantity. The median is the headline; min and max
/// are recorded, because a ratio that swings under `amd-pstate-epp` cannot be a
/// gate on one sample (`M-280`).
const REPS: usize = 5;

/// Untimed passes before anything is measured.
const WARMUP: usize = 1;

/// The six faces of a cube as `(normal axis, side)`, in the order the residual
/// six-vector is laid out.
const FACES: [(usize, u8); FACE_COUNT] = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)];

/// `|| C ||_2 <= 4 * sqrt(6)`: each of a face circulation's four terms is a
/// component of a unit vector, so `|C| <= 4`, over six faces.
///
/// This is the **extremal** normaliser and it is what `curl_residual_normalised`
/// is divided by. Asserted against `4.0 * 6.0_f64.sqrt()` before any row is
/// written.
const RESIDUAL_SCALE: Scalar = 9.797_958_971_132_712;

/// A single face carrying a full normal reversal on opposite edges contributes
/// `|C| = 4` alone, so its cell's normalised residual would be at least
/// `4 / (4 * sqrt(6)) = 1 / sqrt(6) = 0.408248`.
///
/// Reported per row as `reversal_floor_reached` rather than asserted: the
/// configuration needs two normals *exactly* axis-aligned and exactly opposed,
/// and no reference field reaches it. It is here because it is what bounds
/// [`CONTROL_CEILING`] from above.
const REVERSAL_FLOOR: Scalar = 4.0 / RESIDUAL_SCALE;

/// The natural **unit** of the quantity: one whole unit vector's worth of
/// failure to cancel around one loop, `|| C ||_2 = 1`.
///
/// The per-cell `lambda` map reads this scale and not [`RESIDUAL_SCALE`]. The
/// extremal bound is a worst case over configurations no smooth field produces,
/// so a map anchored on it puts every real cell in its bottom decade and
/// degenerates into a global *decrease* of `lambda` — which cannot answer C3,
/// because C3 asks whether a **per-cell** regularizer beats a global one. `1` is
/// a property of the construction rather than of the data: the terms are
/// components of unit vectors, so one unit of disagreement is the quantity's own
/// unit.
const UNIT_CIRCULATION: Scalar = 1.0;

/// "Near-zero" for the registration's zero-second-vertex control.
///
/// Below [`REVERSAL_FLOOR`], so that passing it means no cell in the mean
/// carries a full reversal on any face, and above every zero-second-vertex
/// field's measured mean by an order of magnitude. The instrument is separately
/// required to *exceed* this ceiling somewhere, or "near-zero" separates nothing.
const CONTROL_CEILING: Scalar = 0.25;

/// C1's bar, from the registration: under 2% of extraction.
const COST_CEILING: Scalar = 0.02;

/// C2's bar, from the registration: AUC above 0.8.
const AUC_BAR: Scalar = 0.8;

/// The bottom of the per-cell `lambda` bracket — one decade below the shipped
/// `solve::LAMBDA`.
const LAMBDA_MIN: Scalar = 0.001;

/// The top of the bracket — Subgrid MT's QEF regularizer, one of the three
/// constants the corpus circulates (`solve.rs:56-57`).
const LAMBDA_MAX: Scalar = 0.1;

/// The five fields `M-60` measured at **exactly zero** extra vertices at every
/// resolution, and the population the registration's vacuity control names.
const M60_ZERO_FIELDS: [&str; 5] = [
    "sphere",
    "torus",
    "box_exact",
    "csg_difference",
    "thin_plate",
];

/// The two fields `M-60` measured at a non-zero rate, and the two C2 names.
const M60_POSITIVE_FIELDS: [&str; 2] = ["gyroid", "fbm_terrain"];

// ─── private crate mechanisms, copied rather than made `pub` ────────────────

/// `cube::corner_offset`, which is `pub(crate)` (`cube.rs:149`) and
/// `crates/isomesh/src/**` is read-only this phase. Corner `i` sits at
/// `(i & 1, (i >> 1) & 1, (i >> 2) & 1)` (`cube.rs:12-14`).
#[inline]
const fn corner_offset(corner: u8) -> [u32; 3] {
    [
        (corner & 1) as u32,
        ((corner >> 1) & 1) as u32,
        ((corner >> 2) & 1) as u32,
    ]
}

/// `dual_contouring::apply_clamp` under `Clamp::ToCell`, which is
/// `DualContouring::new`'s default and therefore the shipped path
/// (`dual_contouring.rs:184-203`). Verified bit-for-bit against the shipped
/// extractor on every row rather than trusted.
#[inline]
fn clamp_to_cell(x: [Scalar; 3], cell_origin: [Scalar; 3], cell_size: Scalar) -> [Scalar; 3] {
    let half = cell_size * 0.5;
    let inset = half * (1.0 - CLAMP_EPSILON);
    let mut out = x;
    for (axis, slot) in out.iter_mut().enumerate() {
        let centre = cell_origin[axis] + half;
        *slot = slot.clamp(centre - inset, centre + inset);
    }
    out
}

/// `Qef::place`'s cell origin (`dual_contouring.rs:155-159`).
#[inline]
fn origin_of_cell(base: [u32; 3], origin: [Scalar; 3], cell_size: Scalar) -> [Scalar; 3] {
    [
        origin[0] + cell_size * Scalar::from(base[0]),
        origin[1] + cell_size * Scalar::from(base[1]),
        origin[2] + cell_size * Scalar::from(base[2]),
    ]
}

/// `marching_cubes::unit_gradient`, over `vec3::length` and `vec3::scale`.
#[inline]
fn unit_gradient<S: Sdf<Scalar = Scalar>>(sdf: &S, position: [Scalar; 3]) -> [Scalar; 3] {
    let g = sdf.gradient(position);
    let inv = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt().recip();
    [g[0] * inv, g[1] * inv, g[2] * inv]
}

// ─── the curl residual ─────────────────────────────────────────────────────

/// The midpoint of one cube edge, in world space.
#[inline]
fn edge_midpoint(edge: usize, cell_origin: [Scalar; 3], cell_size: Scalar) -> [Scalar; 3] {
    let [lo, hi] = EDGE_CORNERS[edge];
    let a = corner_offset(lo);
    let b = corner_offset(hi);
    let mut out = cell_origin;
    for (axis, slot) in out.iter_mut().enumerate() {
        *slot += cell_size * (Scalar::from(a[axis]) + Scalar::from(b[axis])) * 0.5;
    }
    out
}

/// One term of one face's circulation: which edge's normal, which component of
/// it, and whether it enters with a `+` or a `-`.
type Term = (u8, u8, bool);

/// The `6 x 4` terms of the six face circulations, derived at **compile time**
/// from the crate's own `edge_on_face`, `EDGE_AXIS` and `EDGE_CORNERS`.
///
/// A table and not a per-cell scan, and that is about the measurement rather
/// than about tidiness: the registration's C1 is a claim about *"a few dozen
/// flops"*, and re-deriving cube topology per cell — six faces times twelve
/// `edge_on_face` tests, plus a `corner_offset` per hit — measures the scan
/// instead. Reduced to this table the inner loop is 24 signed adds, six
/// multiply-adds and one `sqrt`: about fifty flops, and nothing else. Deriving
/// it rather than transcribing it is the crate's own habit (`table::CASES` is
/// built by a `const fn` for the same reason).
const FACE_TERMS: [[Term; 4]; FACE_COUNT] = build_face_terms();

/// Counter-clockwise as seen from `+w`, so the four tangents are `+u`, `+v`,
/// `-u`, `-v` — see the module docs for the derivation and `cube.rs:81-87` for
/// the handedness this borrows.
///
/// # Panics
///
/// At compile time, if any face does not collect exactly four edges, which
/// would mean the traversal disagrees with `edge_on_face` about what a face is.
const fn build_face_terms() -> [[Term; 4]; FACE_COUNT] {
    let mut out = [[(0u8, 0u8, true); 4]; FACE_COUNT];
    let mut face = 0usize;
    while face < FACE_COUNT {
        let (axis, side) = FACES[face];
        let u = (axis + 1) % 3;
        let v = (axis + 2) % 3;
        let mut slot = 0usize;
        let mut edge = 0u8;
        while (edge as usize) < EDGE_COUNT {
            if edge_on_face(edge, axis, side) {
                let low = corner_offset(EDGE_CORNERS[edge as usize][0]);
                // A `u`-edge enters as `+u` along the bottom (`v` offset 0) and
                // `-u` along the top; a `v`-edge as `+v` on the right (`u`
                // offset 1) and `-v` on the left.
                let term = if EDGE_AXIS[edge as usize] as usize == u {
                    (edge, u as u8, low[v] == 0)
                } else {
                    (edge, v as u8, low[u] == 1)
                };
                out[face][slot] = term;
                slot += 1;
            }
            edge += 1;
        }
        assert!(slot == 4, "a cube face has four edges");
        face += 1;
    }
    out
}

/// `|| C ||_2` over the six faces, in units of the edge length.
#[inline]
fn curl_residual(normals: &[[Scalar; 3]; EDGE_COUNT]) -> Scalar {
    let mut total = 0.0;
    for terms in &FACE_TERMS {
        let mut circulation = 0.0;
        for &(edge, component, positive) in terms {
            let n = normals[edge as usize][component as usize];
            circulation += if positive { n } else { -n };
        }
        total += circulation * circulation;
    }
    total.sqrt()
}

/// Bit `e` set when the surface crosses edge `e`, by the same sign test
/// `HermiteCell::from_corners` uses (`hermite.rs:87-91`).
fn cut_mask(corner: &[Scalar; 8]) -> u16 {
    let mut mask = 0u16;
    for (edge, [lo, hi]) in EDGE_CORNERS.iter().copied().enumerate() {
        if is_inside(corner[lo as usize]) != is_inside(corner[hi as usize]) {
            mask |= 1 << edge;
        }
    }
    mask
}

/// The four edges of each face, as a mask, derived at compile time.
const FACE_EDGE_MASKS: [u16; FACE_COUNT] = build_face_edge_masks();

const fn build_face_edge_masks() -> [u16; FACE_COUNT] {
    let mut out = [0u16; FACE_COUNT];
    let mut face = 0usize;
    while face < FACE_COUNT {
        let (axis, side) = FACES[face];
        let mut edge = 0u8;
        while (edge as usize) < EDGE_COUNT {
            if edge_on_face(edge, axis, side) {
                out[face] |= 1 << edge;
            }
            edge += 1;
        }
        face += 1;
    }
    out
}

/// How many of the cell's six faces have **all four** of their edges cut.
///
/// Recorded per row because it is the measured answer to "why sample all twelve
/// midpoints instead of reusing the crossings": a face's circulation needs a
/// whole loop.
fn complete_faces(cut: u16) -> u32 {
    let mut count = 0;
    for mask in FACE_EDGE_MASKS {
        if cut & mask == mask {
            count += 1;
        }
    }
    count
}

/// How many vertices `ManifoldDualContouring` places in this cell.
///
/// `CycleQef::place` under `FaceAmbiguity::Separate`, which is
/// `ManifoldDualContouring::new`'s default: `ambiguous = 0`, then the cycles of
/// `segment_links` (`manifold_dual_contouring.rs:214-258`). One
/// `push_component`, and therefore one vertex, per cycle.
fn cycle_count(case: u8, corner: &[Scalar; 8]) -> u32 {
    let next = segment_links(case, joined_mask(corner, 0));
    let mut visited = 0u16;
    let mut cycles = 0;
    for (start, &link) in next.iter().enumerate() {
        if link == NO_EDGE || visited & (1 << start) != 0 {
            continue;
        }
        let mut current = start as u8;
        while visited & (1 << current) == 0 {
            visited |= 1 << current;
            current = next[current as usize];
        }
        cycles += 1;
    }
    cycles
}

/// The QEF objective `E(x) = sum (n_i . (x - p_i))^2` over the cell's crossings.
///
/// The quantity `solve_with` minimises a regularized form of (`solve.rs:8-11`),
/// evaluated at the vertex that was actually placed.
fn qef_energy(cell: &HermiteCell<Scalar>, x: [Scalar; 3]) -> Scalar {
    let mut energy = 0.0;
    for crossing in cell.iter() {
        let d = crossing.normal[0] * (x[0] - crossing.position[0])
            + crossing.normal[1] * (x[1] - crossing.position[1])
            + crossing.normal[2] * (x[2] - crossing.position[2]);
        energy += d * d;
    }
    energy
}

// ─── statistics ────────────────────────────────────────────────────────────

/// `(min, median, max)` of a sample, sorted by `f64::total_cmp`.
///
/// # Panics
///
/// If `samples` is empty, which would mean a timed loop ran zero repeats.
fn spread(samples: &[Scalar]) -> (Scalar, Scalar, Scalar) {
    assert!(!samples.is_empty(), "a timed quantity has no repeats");
    let mut sorted = samples.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    (
        sorted[0],
        sorted[sorted.len() / 2],
        sorted[sorted.len() - 1],
    )
}

/// The Mann-Whitney AUC of `scores` against the labels in `positive`.
///
/// `AUC = (R+ - n+(n+ + 1)/2) / (n+ * n-)` on ranks, with **tied scores taking
/// their average rank** — which matters here because `box_exact`'s planar cells
/// give a residual of exactly zero and would otherwise be ordered by index.
///
/// # Panics
///
/// If either class is empty; the caller checks that and records `unreachable`.
fn mann_whitney_auc(scores: &[Scalar], positive: &[bool]) -> Scalar {
    assert_eq!(scores.len(), positive.len(), "one label per score");
    let mut order: Vec<u32> = (0..scores.len() as u32).collect();
    order.sort_by(|a, b| scores[*a as usize].total_cmp(&scores[*b as usize]));

    let mut rank = vec![0.0; scores.len()];
    let mut i = 0usize;
    while i < order.len() {
        let mut j = i + 1;
        while j < order.len() && scores[order[j] as usize] == scores[order[i] as usize] {
            j += 1;
        }
        // Ranks are one-based, so this run occupies `i+1 ..= j` and their mean
        // is `(i + 1 + j) / 2`.
        let average = (i + 1 + j) as Scalar / 2.0;
        for slot in &order[i..j] {
            rank[*slot as usize] = average;
        }
        i = j;
    }

    let mut rank_sum = 0.0;
    let mut n_pos = 0usize;
    for (slot, &is_pos) in positive.iter().enumerate() {
        if is_pos {
            rank_sum += rank[slot];
            n_pos += 1;
        }
    }
    let n_neg = positive.len() - n_pos;
    assert!(n_pos > 0 && n_neg > 0, "the AUC needs both classes");
    let u = rank_sum - (n_pos as Scalar) * (n_pos as Scalar + 1.0) / 2.0;
    u / (n_pos as Scalar * n_neg as Scalar)
}

/// The mean of `values`, or zero for an empty sample.
fn mean(values: &[Scalar]) -> Scalar {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().sum::<Scalar>() / values.len() as Scalar
}

// ─── one row ───────────────────────────────────────────────────────────────

/// Everything measured for one (field, resolution).
struct Row {
    field: &'static str,
    samples: u32,
    cells: usize,
    second_vertex_cells: usize,
    curl_residual: Scalar,
    curl_residual_normalised: Scalar,
    residual_max: Scalar,
    residual_median: Scalar,
    mean_residual_positive: Scalar,
    mean_residual_negative: Scalar,
    auc: Option<Scalar>,
    complete_faces_mean: Scalar,
    lambda_global: Scalar,
    lambda_per_cell: Scalar,
    lambda_cell_min: Scalar,
    lambda_cell_max: Scalar,
    over_shipped_share: Scalar,
    lambda_score_mean: Scalar,
    reversal_floor_cells: usize,
    qef_global: Scalar,
    qef_per_cell: Scalar,
    sharp_global: Scalar,
    sharp_per_cell: Scalar,
    sharp_worst_global: Scalar,
    sharp_worst_per_cell: Scalar,
    sharp_delta_top_decile: Scalar,
    self_intersections_global: u64,
    self_intersections_per_cell: u64,
    triangles: usize,
    dc_vertices: usize,
    mdc_vertices: usize,
    mirror_vertices: usize,
    mirror_mismatches: usize,
    extra_predicted: usize,
    curl_ms: Scalar,
    curl_ms_min: Scalar,
    curl_ms_max: Scalar,
    curl_ms_sampled: Scalar,
    extract_ms: Scalar,
    extract_ms_min: Scalar,
    extract_ms_max: Scalar,
}

impl Row {
    /// The registered reading of C1's cost: the residual arithmetic over normals
    /// the extractor's own struct already holds, against whole extraction.
    fn curl_share(&self) -> Scalar {
        self.curl_ms / self.extract_ms
    }

    /// The same share, paying for the residual's own twelve gradients per cell.
    fn curl_share_sampled(&self) -> Scalar {
        self.curl_ms_sampled / self.extract_ms
    }

    /// `mdc_vertices - dc_vertices`: M-60's own "extra vertices over plain Dual
    /// Contouring", from the two live extractors.
    fn extra_measured(&self) -> i64 {
        self.mdc_vertices as i64 - self.dc_vertices as i64
    }

    fn sharpness_delta(&self) -> Scalar {
        self.sharp_per_cell - self.sharp_global
    }

    fn qef_delta(&self) -> Scalar {
        self.qef_per_cell - self.qef_global
    }

    fn self_intersections_delta(&self) -> i64 {
        self.self_intersections_per_cell as i64 - self.self_intersections_global as i64
    }

    fn c1(&self) -> bool {
        self.curl_share() < COST_CEILING
    }

    fn c1_sampled(&self) -> bool {
        self.curl_share_sampled() < COST_CEILING
    }

    /// "This row supports C2": the AUC clears 0.8 where there is one, and where
    /// there is not, the row discharges its half of the registered control.
    fn c2(&self) -> bool {
        match self.auc {
            Some(auc) => auc > AUC_BAR,
            None => self.curl_residual_normalised < CONTROL_CEILING,
        }
    }

    /// A strict win on one axis with no loss on the other. An exact tie is
    /// neither.
    fn c3(&self) -> bool {
        let wins_sharp = self.sharp_per_cell < self.sharp_global;
        let loses_sharp = self.sharp_per_cell > self.sharp_global;
        let wins_intersections = self.self_intersections_per_cell < self.self_intersections_global;
        let loses_intersections = self.self_intersections_per_cell > self.self_intersections_global;
        (wins_sharp && !loses_intersections) || (wins_intersections && !loses_sharp)
    }
}

/// Measure one (field, resolution).
fn measure<F>(field_name: &'static str, field: &F, samples: u32) -> Row
where
    F: ReferenceField + Sdf<Scalar = Scalar>,
{
    let (shape, origin, h) = common::grid::<Scalar, _>(field, samples);
    let size = shape.size();
    let cells_per_axis = [size[0] - 1, size[1] - 1, size[2] - 1];

    // ─── the value grid, sampled exactly as `sdf::sample_grid` does ─────────
    // `origin + cell_size * f64::from(x)`, `x` innermost (`sdf.rs:180-193`), so
    // the corner values handed to the solve are bit-identical to the shipped
    // extractor's.
    let n = samples as usize;
    let mut values = Vec::with_capacity(n * n * n);
    for z in 0..samples {
        for y in 0..samples {
            for x in 0..samples {
                values.push(field.sample([
                    origin[0] + h * Scalar::from(x),
                    origin[1] + h * Scalar::from(y),
                    origin[2] + h * Scalar::from(z),
                ]));
            }
        }
    }

    // ─── active cells, in `place_vertices`' lexicographic order ─────────────
    // "A cell is active exactly when its eight corner signs are not all equal"
    // (`dual.rs:407-409`); the walk order is `(z, y, x)` ascending
    // (`dual.rs:487-497`), which is what makes vertex `i` active cell `i`.
    let mut bases: Vec<[u32; 3]> = Vec::new();
    let mut corners: Vec<[Scalar; 8]> = Vec::new();
    for z in 0..cells_per_axis[2] {
        for y in 0..cells_per_axis[1] {
            for x in 0..cells_per_axis[0] {
                let base = [x, y, z];
                let mut corner = [0.0; 8];
                for (c, slot) in corner.iter_mut().enumerate() {
                    let o = corner_offset(c as u8);
                    let p = [base[0] + o[0], base[1] + o[1], base[2] + o[2]];
                    *slot = values[p[0] as usize + n * (p[1] as usize + n * p[2] as usize)];
                }
                let inside = corner.iter().filter(|v| is_inside(**v)).count();
                if inside == 0 || inside == 8 {
                    continue;
                }
                bases.push(base);
                corners.push(corner);
            }
        }
    }
    assert!(
        !bases.is_empty(),
        "VOID: {field_name} at {samples}^3 has no active cell, so every column in this row is \
         computed over an empty population"
    );

    // ─── the per-cell pass ─────────────────────────────────────────────────
    let lambda_global = solve::LAMBDA;
    let count = bases.len();
    let mut normals: Vec<[[Scalar; 3]; EDGE_COUNT]> = Vec::with_capacity(count);
    let mut residuals: Vec<Scalar> = Vec::with_capacity(count);
    let mut scores: Vec<Scalar> = Vec::with_capacity(count);
    let mut lambda_scores: Vec<Scalar> = Vec::with_capacity(count);
    let mut positive: Vec<bool> = Vec::with_capacity(count);
    let mut lambdas: Vec<Scalar> = Vec::with_capacity(count);
    let mut positions_global: Vec<[Scalar; 3]> = Vec::with_capacity(count);
    let mut positions_per_cell: Vec<[Scalar; 3]> = Vec::with_capacity(count);
    let mut sharp_deltas: Vec<Scalar> = Vec::with_capacity(count);
    let mut qef_global_terms: Vec<Scalar> = Vec::with_capacity(count);
    let mut qef_per_cell_terms: Vec<Scalar> = Vec::with_capacity(count);
    let mut sharp_global_terms: Vec<Scalar> = Vec::with_capacity(count);
    let mut sharp_per_cell_terms: Vec<Scalar> = Vec::with_capacity(count);
    let mut complete_total = 0u64;
    let mut extra_predicted = 0usize;

    for (base, corner) in bases.iter().copied().zip(&corners) {
        let cell_origin = origin_of_cell(base, origin, h);
        let hermite = HermiteCell::from_corners(field, corner, cell_origin, h);

        let mut edge_normals = [[0.0; 3]; EDGE_COUNT];
        for (edge, slot) in edge_normals.iter_mut().enumerate() {
            *slot = unit_gradient(field, edge_midpoint(edge, cell_origin, h));
        }
        let residual = curl_residual(&edge_normals);
        // Two normalisations of one residual, both fixed in the module docs. The
        // extremal one is the reported score and what the AUC ranks — ranks are
        // invariant under any monotone rescale, so C2 is unaffected by the
        // choice. The unit one is what the regularizer reads.
        let score = (residual / RESIDUAL_SCALE).clamp(0.0, 1.0);
        let lambda_score = (residual / UNIT_CIRCULATION).clamp(0.0, 1.0);

        let mut case = 0u8;
        for (c, value) in corner.iter().enumerate() {
            if is_inside(*value) {
                case |= 1 << c;
            }
        }
        let cycles = cycle_count(case, corner);
        extra_predicted += cycles.saturating_sub(1) as usize;
        complete_total += u64::from(complete_faces(cut_mask(corner)));

        // The map is fixed in the module docs: `0.001 * 100^lambda_score`, whose
        // value at `lambda_score = 1/2` is exactly `solve::LAMBDA`.
        let lambda = LAMBDA_MIN * (LAMBDA_MAX / LAMBDA_MIN).powf(lambda_score);

        // Both arms are the shipped `solve_with` plus the shipped clamp; only
        // the regularizer differs.
        let x_global = clamp_to_cell(
            solve::solve_with(&hermite, lambda_global)
                .expect("an active cell has at least one crossing"),
            cell_origin,
            h,
        );
        let x_per_cell = clamp_to_cell(
            solve::solve_with(&hermite, lambda).expect("an active cell has at least one crossing"),
            cell_origin,
            h,
        );

        qef_global_terms.push(qef_energy(&hermite, x_global) / (h * h));
        qef_per_cell_terms.push(qef_energy(&hermite, x_per_cell) / (h * h));
        let off_global = field.sample(x_global).abs() / h;
        let off_per_cell = field.sample(x_per_cell).abs() / h;
        sharp_global_terms.push(off_global);
        sharp_per_cell_terms.push(off_per_cell);
        sharp_deltas.push(off_per_cell - off_global);

        normals.push(edge_normals);
        residuals.push(residual);
        scores.push(score);
        lambda_scores.push(lambda_score);
        positive.push(cycles >= 2);
        lambdas.push(lambda);
        positions_global.push(x_global);
        positions_per_cell.push(x_per_cell);
    }

    // ─── the two shipped extractors ────────────────────────────────────────
    let mut dc = DualContouring::<Scalar>::new();
    let mut mesh = MeshBuffer::<Scalar>::new();
    dc.extract(field, &shape, origin, h, &mut mesh)
        .expect("dual contouring meshes a reference grid");

    let mut mdc = ManifoldDualContouring::<Scalar>::new();
    let mut mdc_mesh = MeshBuffer::<Scalar>::new();
    mdc.extract(field, &shape, origin, h, &mut mdc_mesh)
        .expect("manifold dual contouring meshes a reference grid");

    // The mirror check: same positions, in the same order, as IEEE bit patterns.
    let mut mirror_mismatches = 0usize;
    for (shipped, mirrored) in mesh.positions.iter().zip(&positions_global) {
        if shipped[0].to_bits() != mirrored[0].to_bits()
            || shipped[1].to_bits() != mirrored[1].to_bits()
            || shipped[2].to_bits() != mirrored[2].to_bits()
        {
            mirror_mismatches += 1;
        }
    }

    // ─── the two C3 axes ───────────────────────────────────────────────────
    let self_intersections_global = self_intersections(&positions_global, &mesh.indices, h)
        .expect("the dual mesh's cell size is the grid it came from")
        .count();
    let self_intersections_per_cell = self_intersections(&positions_per_cell, &mesh.indices, h)
        .expect("the dual mesh's cell size is the grid it came from")
        .count();

    // ─── C2's statistic ────────────────────────────────────────────────────
    let second_vertex_cells = positive.iter().filter(|p| **p).count();
    let auc = if second_vertex_cells == 0 || second_vertex_cells == count {
        None
    } else {
        Some(mann_whitney_auc(&scores, &positive))
    };

    let positive_scores: Vec<Scalar> = scores
        .iter()
        .zip(&positive)
        .filter(|(_, p)| **p)
        .map(|(s, _)| *s)
        .collect();
    let negative_scores: Vec<Scalar> = scores
        .iter()
        .zip(&positive)
        .filter(|(_, p)| !**p)
        .map(|(s, _)| *s)
        .collect();

    // The highest-residual decile, where the per-cell `lambda` differs most from
    // the shipped one and therefore where its mechanism is visible.
    let mut by_score: Vec<u32> = (0..count as u32).collect();
    by_score.sort_by(|a, b| scores[*a as usize].total_cmp(&scores[*b as usize]));
    let decile = (count / 10).max(1);
    let sharp_delta_top_decile = by_score[count - decile..]
        .iter()
        .map(|i| sharp_deltas[*i as usize])
        .sum::<Scalar>()
        / decile as Scalar;

    let mut sorted_scores = scores.clone();
    sorted_scores.sort_by(|a, b| a.total_cmp(b));

    // ─── timings ───────────────────────────────────────────────────────────
    // (a) the registered reading: the arithmetic, over normals already held.
    let mut curl_samples = Vec::with_capacity(REPS);
    for rep in 0..WARMUP + REPS {
        let started = Instant::now();
        let mut accumulator = 0.0;
        for edge_normals in &normals {
            accumulator += curl_residual(black_box(edge_normals));
        }
        let ms = started.elapsed().as_secs_f64() * 1e3;
        black_box(accumulator);
        if rep >= WARMUP {
            curl_samples.push(ms);
        }
    }

    // (b) the same, paying for its own twelve gradients per active cell.
    let mut curl_sampled_samples = Vec::with_capacity(REPS);
    for rep in 0..WARMUP + REPS {
        let started = Instant::now();
        let mut accumulator = 0.0;
        for base in &bases {
            let cell_origin = origin_of_cell(*base, origin, h);
            let mut edge_normals = [[0.0; 3]; EDGE_COUNT];
            for (edge, slot) in edge_normals.iter_mut().enumerate() {
                *slot = unit_gradient(field, edge_midpoint(edge, cell_origin, h));
            }
            accumulator += curl_residual(black_box(&edge_normals));
        }
        let ms = started.elapsed().as_secs_f64() * 1e3;
        black_box(accumulator);
        if rep >= WARMUP {
            curl_sampled_samples.push(ms);
        }
    }

    // (c) the denominator: whole extraction on the same grid.
    let mut timed_dc = DualContouring::<Scalar>::new();
    let mut timed_mesh = MeshBuffer::<Scalar>::new();
    let mut extract_samples = Vec::with_capacity(REPS);
    for rep in 0..WARMUP + REPS {
        timed_mesh.reset();
        let started = Instant::now();
        timed_dc
            .extract(field, &shape, origin, h, &mut timed_mesh)
            .expect("dual contouring meshes a reference grid");
        let ms = started.elapsed().as_secs_f64() * 1e3;
        black_box(timed_mesh.vertex_count());
        if rep >= WARMUP {
            extract_samples.push(ms);
        }
    }

    let (curl_ms_min, curl_ms, curl_ms_max) = spread(&curl_samples);
    let (_, curl_ms_sampled, _) = spread(&curl_sampled_samples);
    let (extract_ms_min, extract_ms, extract_ms_max) = spread(&extract_samples);

    let mut lambda_cell_min = Scalar::INFINITY;
    let mut lambda_cell_max: Scalar = 0.0;
    for lambda in &lambdas {
        lambda_cell_min = lambda_cell_min.min(*lambda);
        lambda_cell_max = lambda_cell_max.max(*lambda);
    }
    // `lambda = LAMBDA_MIN * 100^lambda_score` is log-linear in the score, so
    // its geometric mean is `LAMBDA_MIN * 100^(mean score)` — reported rather
    // than an arithmetic mean, which would be a mean of the wrong quantity.
    let lambda_per_cell = LAMBDA_MIN * (LAMBDA_MAX / LAMBDA_MIN).powf(mean(&lambda_scores));
    let reversal_floor_cells = scores.iter().filter(|s| **s >= REVERSAL_FLOOR).count();

    let residual_mean = mean(&residuals);

    Row {
        field: field_name,
        samples,
        cells: count,
        second_vertex_cells,
        curl_residual: residual_mean,
        curl_residual_normalised: residual_mean / RESIDUAL_SCALE,
        residual_max: sorted_scores[count - 1],
        residual_median: sorted_scores[count / 2],
        mean_residual_positive: mean(&positive_scores),
        mean_residual_negative: mean(&negative_scores),
        auc,
        complete_faces_mean: complete_total as Scalar / count as Scalar,
        lambda_global,
        lambda_per_cell,
        lambda_cell_min,
        lambda_cell_max,
        over_shipped_share: lambda_scores.iter().filter(|s| **s > 0.5).count() as Scalar
            / count as Scalar,
        lambda_score_mean: mean(&lambda_scores),
        reversal_floor_cells,
        qef_global: mean(&qef_global_terms),
        qef_per_cell: mean(&qef_per_cell_terms),
        sharp_global: mean(&sharp_global_terms),
        sharp_per_cell: mean(&sharp_per_cell_terms),
        sharp_worst_global: sharp_global_terms.iter().copied().fold(0.0, Scalar::max),
        sharp_worst_per_cell: sharp_per_cell_terms.iter().copied().fold(0.0, Scalar::max),
        sharp_delta_top_decile,
        self_intersections_global,
        self_intersections_per_cell,
        triangles: mesh.triangle_count(),
        dc_vertices: mesh.vertex_count(),
        mdc_vertices: mdc_mesh.vertex_count(),
        mirror_vertices: count,
        mirror_mismatches,
        extra_predicted,
        curl_ms,
        curl_ms_min,
        curl_ms_max,
        curl_ms_sampled,
        extract_ms,
        extract_ms_min,
        extract_ms_max,
    }
}

// ─── the vacuity controls ──────────────────────────────────────────────────

/// Every registered control, run before the first `run.record`.
fn check_controls(rows: &[Row]) {
    let exact_scale = 4.0 * 6.0_f64.sqrt();
    assert!(
        (RESIDUAL_SCALE - exact_scale).abs() < 1e-12,
        "VOID: RESIDUAL_SCALE is {RESIDUAL_SCALE} and 4*sqrt(6) is {exact_scale}, so every \
         `curl_residual_normalised` in this CSV is divided by the wrong bound and both \
         REVERSAL_FLOOR and CONTROL_CEILING are about nothing"
    );

    for row in rows {
        assert!(
            row.mirror_vertices == row.dc_vertices && row.mirror_mismatches == 0,
            "VOID: the bench-local global-lambda arm is not the shipped DualContouring on {} at \
             {}^3 — {} shipped vertices against {} mirrored, {} positions differing in their \
             IEEE bits — so `sharpness_delta`, `qef_residual_delta` and \
             `self_intersections_delta` are differences against something that is not the \
             shipped path",
            row.field,
            row.samples,
            row.dc_vertices,
            row.mirror_vertices,
            row.mirror_mismatches
        );
        assert!(
            row.extra_measured() == row.extra_predicted as i64 && row.dc_vertices == row.cells,
            "VOID: the cycle count is not the crate's own second-vertex predicate on {} at {}^3 \
             — the two live extractors differ by {} vertices where `segment_links` predicts {}, \
             and DualContouring emitted {} vertices for {} active cells — so \
             `second_vertex_cells` is not the population M-60 counted and C2 is scored on the \
             wrong label",
            row.field,
            row.samples,
            row.extra_measured(),
            row.extra_predicted,
            row.dc_vertices,
            row.cells
        );
    }

    for name in M60_ZERO_FIELDS {
        let here: Vec<&Row> = rows.iter().filter(|r| r.field == name).collect();
        assert_eq!(
            here.len(),
            RESOLUTIONS.len(),
            "VOID: {name} is missing from the sweep, so the registration's control population is \
             incomplete"
        );
        for row in here {
            assert_eq!(
                row.second_vertex_cells, 0,
                "VOID: {name} at {}^3 has {} second-vertex cells where M-60 measured exactly 0 at \
                 every resolution, so the registration's zero-second-vertex control names a \
                 population that does not exist and C2's negatives-only control is void",
                row.samples, row.second_vertex_cells
            );
            assert!(
                row.curl_residual_normalised < CONTROL_CEILING,
                "VOID: {name} at {}^3 has zero second-vertex cells and a mean normalised curl \
                 residual of {:.6}, at or above the ceiling {CONTROL_CEILING} — so the fields \
                 with no negatives do not report a near-zero residual and C2's separation is \
                 measured on a population with no negatives (the registration's own vacuity \
                 control). The ceiling sits below REVERSAL_FLOOR = {REVERSAL_FLOOR:.6}, the \
                 score a single full normal reversal on one face forces by itself",
                row.samples,
                row.curl_residual_normalised
            );
        }
    }

    for name in M60_POSITIVE_FIELDS {
        for row in rows.iter().filter(|r| r.field == name) {
            assert!(
                row.second_vertex_cells > 0,
                "VOID: {name} at {}^3 has no second-vertex cell where M-60 measured a non-zero \
                 rate, so C2's AUC is unreachable at one of the two fields C2 names",
                row.samples
            );
            assert!(
                row.cells > row.second_vertex_cells,
                "VOID: every one of {name}'s {} active cells at {}^3 needs a second vertex, so \
                 the AUC has no negative class",
                row.cells,
                row.samples
            );
        }
    }

    // M-44's rule, and the threshold is CONTROL_CEILING rather than
    // REVERSAL_FLOOR on purpose: the control above certifies the five
    // zero-second-vertex fields as *below* the ceiling, and that certificate is
    // empty unless the instrument can also reach *above* it. Scoring this on
    // REVERSAL_FLOOR would instead demand an extremal configuration — two
    // normals exactly axis-aligned and exactly opposed — which no smooth field
    // produces, so it would fail on a fixture that is measuring perfectly well.
    let hottest = rows.iter().map(|r| r.residual_max).fold(0.0, Scalar::max);
    assert!(
        hottest > CONTROL_CEILING,
        "VOID: the largest per-cell normalised curl residual over all {} rows is {hottest:.6}, \
         at or below the near-zero ceiling {CONTROL_CEILING} that the zero-second-vertex control \
         certifies those fields against — so the instrument never says anything but 'near zero' \
         and every residual in this CSV is a zero that could not have been non-zero (M-44)",
        rows.len()
    );

    let widest = rows
        .iter()
        .map(|r| r.lambda_cell_max / r.lambda_cell_min)
        .fold(0.0, Scalar::max);
    assert!(
        widest > 2.0,
        "VOID: the per-cell lambda spans a factor of only {widest:.6} on its widest row, so C3 \
         compares the shipped constant {} against itself",
        solve::LAMBDA
    );
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-173");

    common::experiment::run(prereg, |run| {
        println!(
            "P-173 — the curl residual of the Hermite data, against M-60's second-vertex cells\n\
             \n  lambda(s) = {LAMBDA_MIN} * ({LAMBDA_MAX}/{LAMBDA_MIN})^s on s = ||C||/\
             {UNIT_CIRCULATION}, which is {} at s = 1/2\n  ||C||_2 <= {RESIDUAL_SCALE}; reported \
             score divides by that, reversal floor {REVERSAL_FLOOR:.6}, control ceiling \
             {CONTROL_CEILING}\n",
            solve::LAMBDA
        );
        println!(
            "{:>15} {:>5} {:>7} {:>6} {:>9} {:>9} {:>9} {:>11} {:>9} {:>8} {:>8} {:>6}",
            "field",
            "n",
            "cells",
            "2vtx",
            "r_mean",
            "r_2vtx",
            "r_1vtx",
            "auc",
            "lambda",
            "share",
            "d_sharp",
            "d_si"
        );

        let mut rows: Vec<Row> = Vec::new();
        isomesh::for_each_reference_field!(Scalar, |name, field| {
            for samples in RESOLUTIONS {
                let row = measure(name, &field, samples);
                println!(
                    "{:>15} {:>5} {:>7} {:>6} {:>9.6} {:>9.6} {:>9.6} {:>11} {:>9.6} {:>8.5} \
                     {:>8.4} {:>6}",
                    row.field,
                    row.samples,
                    row.cells,
                    row.second_vertex_cells,
                    row.curl_residual_normalised,
                    row.mean_residual_positive,
                    row.mean_residual_negative,
                    row.auc
                        .map_or_else(|| "unreachable".to_string(), |a| format!("{a:.6}")),
                    row.lambda_per_cell,
                    row.curl_share(),
                    row.sharpness_delta(),
                    row.self_intersections_delta()
                );
                rows.push(row);
            }
        });

        // ─── vacuity controls, before any row is written ────────────────────
        check_controls(&rows);

        let c1_held = rows.iter().filter(|r| r.c1()).count();
        let c2_held = rows.iter().filter(|r| r.c2()).count();
        let c3_held = rows.iter().filter(|r| r.c3()).count();
        let total = rows.len();
        println!(
            "\nC1 {c1_held}/{total} rows under {COST_CEILING} of extraction (registered reading)\n\
             C2 {c2_held}/{total} rows support the clause\n\
             C3 {c3_held}/{total} rows win on one axis without losing the other"
        );

        for row in &rows {
            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.samples.to_string()),
                ("curl_residual", format!("{:.9}", row.curl_residual)),
                (
                    "curl_residual_normalised",
                    format!("{:.9}", row.curl_residual_normalised),
                ),
                ("cells", row.cells.to_string()),
                ("second_vertex_cells", row.second_vertex_cells.to_string()),
                (
                    "separation_auc",
                    row.auc
                        .map_or_else(|| "unreachable".to_string(), |a| format!("{a:.9}")),
                ),
                ("lambda_global", format!("{:.9}", row.lambda_global)),
                ("lambda_per_cell", format!("{:.9}", row.lambda_per_cell)),
                ("qef_residual_delta", format!("{:.9}", row.qef_delta())),
                ("sharpness_delta", format!("{:.9}", row.sharpness_delta())),
                (
                    "self_intersections_delta",
                    row.self_intersections_delta().to_string(),
                ),
                ("curl_ms", format!("{:.6}", row.curl_ms)),
                ("curl_share", format!("{:.9}", row.curl_share())),
                ("c1_holds", row.c1().to_string()),
                ("c2_holds", row.c2().to_string()),
                ("c3_holds", row.c3().to_string()),
                // ── extras (M-273) ──
                ("auc_bar", format!("{AUC_BAR}")),
                ("c1_holds_global", (c1_held == total).to_string()),
                ("c1_holds_sampled", row.c1_sampled().to_string()),
                ("c1_rows_held", c1_held.to_string()),
                ("c2_holds_global", (c2_held == total).to_string()),
                ("c2_rows_held", c2_held.to_string()),
                ("c3_holds_global", (c3_held == total).to_string()),
                ("c3_rows_held", c3_held.to_string()),
                (
                    "complete_faces_mean",
                    format!("{:.6}", row.complete_faces_mean),
                ),
                ("control_ceiling", format!("{CONTROL_CEILING}")),
                ("cost_ceiling", format!("{COST_CEILING}")),
                ("curl_ms_max", format!("{:.6}", row.curl_ms_max)),
                ("curl_ms_min", format!("{:.6}", row.curl_ms_min)),
                ("curl_ms_sampled", format!("{:.6}", row.curl_ms_sampled)),
                (
                    "curl_share_sampled",
                    format!("{:.9}", row.curl_share_sampled()),
                ),
                ("dc_vertices", row.dc_vertices.to_string()),
                ("extra_vertices_measured", row.extra_measured().to_string()),
                ("extra_vertices_predicted", row.extra_predicted.to_string()),
                ("extract_ms", format!("{:.6}", row.extract_ms)),
                ("extract_ms_max", format!("{:.6}", row.extract_ms_max)),
                ("extract_ms_min", format!("{:.6}", row.extract_ms_min)),
                ("lambda_cell_max", format!("{:.9}", row.lambda_cell_max)),
                ("lambda_cell_min", format!("{:.9}", row.lambda_cell_min)),
                ("lambda_max", format!("{LAMBDA_MAX}")),
                ("lambda_min", format!("{LAMBDA_MIN}")),
                ("lambda_score_mean", format!("{:.9}", row.lambda_score_mean)),
                ("lambda_unit_circulation", format!("{UNIT_CIRCULATION}")),
                ("mdc_vertices", row.mdc_vertices.to_string()),
                (
                    "mean_residual_second_vertex",
                    format!("{:.9}", row.mean_residual_positive),
                ),
                (
                    "mean_residual_single_vertex",
                    format!("{:.9}", row.mean_residual_negative),
                ),
                ("mirror_mismatches", row.mirror_mismatches.to_string()),
                ("mirror_vertices", row.mirror_vertices.to_string()),
                (
                    "over_shipped_lambda_share",
                    format!("{:.6}", row.over_shipped_share),
                ),
                ("qef_residual_global", format!("{:.9}", row.qef_global)),
                ("qef_residual_per_cell", format!("{:.9}", row.qef_per_cell)),
                ("repeats", REPS.to_string()),
                ("residual_max", format!("{:.9}", row.residual_max)),
                ("residual_median", format!("{:.9}", row.residual_median)),
                ("reversal_floor", format!("{REVERSAL_FLOOR:.9}")),
                ("reversal_floor_cells", row.reversal_floor_cells.to_string()),
                (
                    "second_vertex_share",
                    format!(
                        "{:.9}",
                        row.second_vertex_cells as Scalar / row.cells as Scalar
                    ),
                ),
                (
                    "self_intersections_global",
                    row.self_intersections_global.to_string(),
                ),
                (
                    "self_intersections_per_cell",
                    row.self_intersections_per_cell.to_string(),
                ),
                ("separation_auc_reachable", row.auc.is_some().to_string()),
                (
                    "sharpness_delta_top_decile",
                    format!("{:.9}", row.sharp_delta_top_decile),
                ),
                ("sharpness_global", format!("{:.9}", row.sharp_global)),
                ("sharpness_per_cell", format!("{:.9}", row.sharp_per_cell)),
                (
                    "sharpness_worst_global",
                    format!("{:.9}", row.sharp_worst_global),
                ),
                (
                    "sharpness_worst_per_cell",
                    format!("{:.9}", row.sharp_worst_per_cell),
                ),
                ("triangles", row.triangles.to_string()),
            ]);
        }
    });
}
