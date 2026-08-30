//! E-317 — the `b*b - 4*a*c` the body-saddle solver has always computed is
//! Cayley's `2x2x2` hyperdeterminant, and this is what its magnitude looks like
//! spread over a grid.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example hyperdeterminant_cells --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 20-50x slower, and this one
//! samples the field twice: once for the census and once inside the extractor.
//!
//! Keys: `1`-`7` field, `[` `]` resolution, `X` the `f32` sign-disagreement
//! cages, `Z` the cells where the `Delta = 0` branch fires, `V` the surface.
//! The rest are the shared keys -- `W` wireframe, `N` normals, `G` domain box,
//! `H` HUD, `Space` pause, `R` re-mesh, `F12` screenshot, `Esc` quit.
//!
//! Under `ISOMESH_CAPTURE` it needs no keyboard: the sequence walks all eight
//! reference fields and, inside each, both rungs of the resolution ladder, so
//! `record_gif.sh`'s default 80 frames is exactly one pass through the whole
//! sixteen-row census and the clip loops. `ISOMESH_FIELD=0..7` pins one field
//! and takes it out of the walk; `ISOMESH_SAMPLES=33` or `=65` pins the rung and
//! takes that out; `ISOMESH_SPIN` adds yaw.
//!
//! ```bash
//! # The clip size is not measured here: this host has no display, so the
//! # orchestrator records it. Every counting/comparing clip keeps its HUD,
//! # which is the whole evidence in this one.
//! ISOMESH_SPIN=0.003 scripts/record_gif.sh hyperdeterminant_cells docs/gifs/e317.gif
//! ```
//!
//! # The identity, which is P-127 and is the phase's headline
//!
//! `crates/isomesh/src/marching_cubes/trilinear.rs:246` reads
//! `b * b - R::TWO * R::TWO * a * c`, over the coefficients built at `:199-214`.
//! Read as a quadratic's discriminant it is unremarkable. Read as a polynomial
//! in the cell's eight corner values it **is** Cayley's `2x2x2`
//! hyperdeterminant, `Det(2,2,2)`, under the indexing `f[u + 2v + 4w]` --
//! `docs/experiments/p-127.csv`, all three clauses held:
//!
//! - **C1** — twelve terms on each side, total degree 4, symbolic difference
//!   exactly the zero polynomial. A committed symbolic check, not a numeric
//!   sample.
//! - **C2** — the same polynomial is `disc(det(A0 + lambda*A1))` for the two
//!   opposite-face `2x2` corner matrices, and for **all three** axis pairings
//!   (`0123|4567`, `0145|2367`, `0246|1357`). That is the mechanism behind
//!   M-206.
//! - **C3** — over 3,481 random 8-tuples in exact rational arithmetic the two
//!   expressions agree with ratio exactly `1` and zero sign disagreements; in
//!   `f32` **14 of those 3,481 flip sign**.
//!
//! And de Silva & Lim (`arXiv:math/0607647`, §6) then fix a meaning for the
//! sign: real tensor rank is **2** on `{Det > 0}` and **3** on `{Det < 0}`. So
//! the crate has computed a real-tensor-rank certificate since M-206 and has
//! never looked at how it is distributed. That distribution is what is on
//! screen.
//!
//! # What is on screen
//!
//! - **The surface, coloured per cell** by `abs(Delta) / max(abs(f_i))^4`. The
//!   ramp is logarithmic over the five decades `1e-8 .. 1e0` -- P-134's own
//!   `threshold_sweep` decades -- from dark blue through blue, green and yellow
//!   to red.
//! - **Flat light grey** is the `Delta = 0` stratum. It is not the cold end of
//!   the ramp; a logarithm has no place to put a zero, so the stratum is drawn
//!   off the ramp entirely. `box_exact` is grey from edge to edge, which is the
//!   honest picture of a polyhedron (below).
//! - **Magenta cages** are cells where the `f32` sign of `Delta` disagrees with
//!   the `f64` reference sign. `X` toggles them; they are on by default because
//!   there are never many -- measured, 4 cells on `csg_difference` at 33 samples
//!   and none at all on six of the eight fields.
//! - **Cyan cages** are cells where `a != 0` and `Delta == 0`, so the branch at
//!   `trilinear.rs:250` actually fires -- P-131's population, which is **every**
//!   cell and not only the surface ones. `Z` toggles them, and `V` hides the
//!   surface, since a cage can sit inside solid rock. Measured: 17 cells on
//!   `csg_difference` at 33 samples and 41 at 65, all of them surface cells, and
//!   **zero on every other reference field**.
//!
//! # Where `:246` actually runs, which is not everywhere
//!
//! `MarchingCubes::new()` defaults to `FaceAmbiguity::Separate` and
//! `InteriorAmbiguity::Ignore`, and in that configuration **line 246 never
//! executes**: `BodySaddles::of` is reached only from `emit_trilinear`, which is
//! gated on `ambiguous != 0 && interior_ambiguity == Trilinear`. So this demo
//! extracts with `AsymptoticDecider` + `Trilinear`, the only configuration in
//! which its subject is live code, and the HUD reports how many cells have an
//! ambiguous face: **27 of 5,240 on `gyroid` at 33 samples**, and zero on four
//! of the eight fields.
//!
//! The heatmap is therefore a census of the *polynomial* over every surface
//! cell, not a trace of the *branch* over the cells that run it. That is exactly
//! the population P-127, P-130, P-133 and P-134 measured, and the two numbers
//! are on the HUD side by side so the difference cannot be read past.
//!
//! # What the number means
//!
//! `Delta` is homogeneous of degree 4 in the eight corner values, so dividing by
//! `max(abs(f_i))^4` is exact rather than approximate: scale all eight corners
//! by `s` and numerator and denominator both scale by `s^4`. P-134 measured that
//! invariance and **C1 held** -- the worst scale error over the eight fields is
//! `2.7e-15`, and it is exact on dyadic scalings.
//!
//! The denominator is formed as `(M*M)*(M*M)` and not as `M*M*M*M`, which is
//! `p-134.csv`'s own association and its stated reason: on a dyadic rescaling
//! that leaves every mantissa in the computation untouched, so "the
//! normalisation is exact" is an exact-zero prediction rather than a tolerance.
//!
//! The number is a *magnitude of ambiguity* where the crate previously had only
//! a sign. It is not a distance and not an error.
//!
//! # Two evaluation routes, one polynomial
//!
//! The panel puts the live `min`/`mean`/`max` of the normalised magnitude
//! directly under `p-134.csv`'s three columns, and they do not agree bit for
//! bit. The gap is a **route** and not a defect, and the panel reports it as a
//! measured relative gap rather than as a verdict word:
//!
//! - `p-134.csv` evaluates `Delta` through `benches/common::poly`'s **expanded**
//!   twelve-term Cayley form, summing terms of equal magnitude and opposite
//!   sign whose multiplication orders differ.
//! - this demo evaluates it through the shipped `BodySaddles::coefficients`,
//!   because the shipped route is its subject -- and it is also P-130's route,
//!   which is what lets the stratum counts be compared for equality.
//!
//! P-127 proved those are the same polynomial *exactly*. P-130's header records
//! that they are **not** the same `f64` computation, which is precisely why it
//! names its route. Measured here over all sixteen rows: the worst relative gap
//! is `7.25e-8`, on `gyroid` at 65 samples where the minimum sits on one cell,
//! and thirteen of the sixteen rows agree to about `1e-10` -- the CSV's own
//! ten-digit rendering. A gap of that size cannot move a cell across a decade of
//! the ramp; a gap of order `1e-1` would be a disagreement about the polynomial,
//! and that is what the bound in [`MAGNITUDE_GAP_BOUND`] separates.
//!
//! The stratum counts and the branch counts carry no such caveat.
//! **All sixteen rows reproduce `p-130.csv`'s partition and `p-131.csv`'s branch
//! counts cell for cell**, and
//! `every_row_reproduces_p130s_partition_and_p131s_branch_counts` asserts it.
//!
//! # What this shows, and what it does not
//!
//! Three separate negative results, and a demo that hid any of them would be
//! worse than no demo.
//!
//! **P-130's C1 and C2 were both FALSIFIED.**
//!
//! - C1 asked for a partition that is stable across resolutions *and* a
//!   `Delta = 0` stratum under 0.1% of surface cells. The zero stratum is
//!   **100% of them on `box_exact`**, 98.4% on `thin_plate` and 83.6% on
//!   `csg_difference`. Nowhere near 0.1%.
//! - C2 asked that real tensor rank carry information about ambiguity that the
//!   8-bit case index does not. It **cannot**, and the reason is arithmetic
//!   rather than empirical: in this crate a cell is ambiguous exactly when
//!   `AMBIGUOUS_FACES[case] != 0` (`marching_cubes/table.rs:202`), so
//!   `ambiguous = g(case)` is a lookup on the case index, and the
//!   data-processing inequality gives `I(R; g(K)) <= I(R; K)` for every field,
//!   resolution and margin. C2 as registered is **unreachable**, not merely
//!   unsupported, and `p-130.csv` carries `c2_unreachable = true` on all 24 rows
//!   with `mi_gap_bits >= 0` asserted as a numeric check on the inequality.
//!
//! What did hold is C3, the control: on `sphere` every surface cell is rank 2
//! and no cell is ambiguous.
//!
//! **P-134's C2 was FALSIFIED.** It asked for a rank correlation above `0.5`
//! between the normalised magnitude and per-cell symmetric Hausdorff error on at
//! least four of eight fields. It cleared the bar on **three** --
//! `gyroid 0.743`, `csg_difference 0.674`, `noise_cavity 0.600` at 65 samples --
//! and `box_exact` is **arithmetically excluded from the population**: its
//! normalised magnitude has *zero variance* across all **1,352** surface cells
//! at 33 samples and all 5,768 at 65, because the field is a polyhedron and
//! every cell of it is affine along the axes the coefficients difference. A
//! correlation against a constant is not a correlation, so `p-134.csv` records
//! `in_correlation_population = false` and
//! `exclusion_reason = constant-normalised-magnitude-polyhedral-field` on those
//! rows rather than counting a zero as a failure to clear the bar. Ambiguity
//! magnitude and geometric error are, on this evidence, different phenomena.
//!
//! **P-133's C2 was FALSIFIED too**, and it is the one that would have made the
//! sign fixable cheaply: a filtered exact predicate was to cost under `1.5x` the
//! naive float evaluation, and it costs **3.8x to 300x** (`p-133.csv`,
//! `overhead_ratio`). C1 and C3 held -- the `f32` sign does disagree with the
//! exact sign on `csg_difference` and `fbm_terrain`, and correcting it moves 6
//! to 12 triangles on `fbm_terrain`.
//!
//! **P-131 held, and it contradicts a comment.** `trilinear.rs:251-254` explains
//! the `discriminant == 0` branch as *"the two hyperbolas touch rather than
//! cross"*. The branch does fire -- 17 cells of 32,768 on `csg_difference` at 33
//! samples, 41 of 262,144 at 65, and **zero on every other reference field** --
//! and where it fires P-131 classified every hit as **border rank 2**, with two
//! cells at 129 samples carrying **true rank 3** (the `W`-state orbit up to
//! `GL(2)^3`). A tangential touch is not what those cells are.
//!
//! # What the reference sign is, and what it is not
//!
//! The `f32` arm is compared against the **`f64` evaluation of the same
//! polynomial on the same corner values**, not against a certified exact
//! predicate. P-133 is the row that built the exact one -- a Shewchuk-style
//! filter at `26u` times the term-magnitude permanent, with an exact fallback --
//! and it measured `sign_disagreements_f64 = 0` on 23 of its 24 rows, the
//! exception being one cell on `csg_difference` at 65 samples. So `f64` is a
//! reference the measurement licenses rather than one this file assumes; the HUD
//! names it `f64 reference` and never calls it exact.
//!
//! The live count also is not expected to equal P-133's. P-133 sampled its grid
//! **at `f32` and widened back to `f64`** so that its two arms share a
//! bit-identical corner set; this demo samples at `f64` because its stratum
//! census has to reproduce P-130's, which did. The two lattices agree on *which*
//! cells are surface cells -- P-133's `cells` column equals P-130's on all
//! sixteen rows this demo can show -- and differ by one `f32` rounding in the
//! corner values, which is precisely the quantity a sign disagreement is
//! sensitive to.
//!
//! Measured, live, over the sixteen rows: **4** disagreements on
//! `csg_difference` at 33 samples and **13** at 65, **1** on `fbm_terrain` at
//! 65, and zero everywhere else. `p-133.csv`'s `f32` rows read 4, 14 and 0 for
//! those three. So the phenomenon, its concentration on one field and its order
//! of magnitude all carry across the lattice change, and the individual counts
//! do not -- which is what a stratum of cells whose `Delta` cancels to within
//! `f32`'s error bound should do. The `X` overlay is therefore empty on six of
//! the eight fields, and that is P-133's C1 rather than a hole in it:
//! `fields_disagreeing_f32` is **2** of 8.
//!
//! # The lattice, and why the HUD can say EQUAL
//!
//! `origin` is the field's own `domain()` minimum, `cell_size` is
//! `(max - min) / (samples - 1)`, and a sample sits at
//! `origin + cell_size * index` with `x` innermost -- `sdf::sample_grid`'s
//! arithmetic, which is also `benches/common::grid`'s and therefore P-130's and
//! P-131's. A cell is a surface cell when its case byte is neither `0x00` nor
//! `0xFF`, from `cube::is_inside`'s strict `value < 0.0`. Nothing here is a
//! re-implementation: `is_inside`, `AMBIGUOUS_FACES` and
//! `BodySaddles::coefficients` are all imported.
//!
//! So the live census is not merely comparable to `p-130.csv`, it is the same
//! computation over the same numbers, and the HUD's `EQUAL` is a gate: a
//! difference of one cell means this file, that CSV or the shipped coefficients
//! have moved.
//!
//! # Colouring cells through a mesh whose vertices are shared
//!
//! Marching Cubes caches one vertex per crossed grid edge, so a vertex belongs
//! to up to four cells and cannot carry a cell's colour. A *triangle* can: every
//! vertex a cell emits lies in that cell's closed box, so a candidate host is
//! any cell whose closed box covers the triangle's whole span. The mesh is
//! therefore rebuilt as a soup, three vertices per triangle, each carrying its
//! cell's colour.
//!
//! **The candidate is not always unique, and the obvious tie-break is wrong.**
//! A triangle lying exactly in a face shared by two cells belongs to both boxes.
//! Picking the far side by `floor(midpoint)` was the first implementation here
//! and `box_exact` measured it wrong: that field is planar on its own faces, the
//! far cell is routinely one whose eight corners are all outside and which
//! emitted nothing, and its surface came out in **two** colours instead of the
//! one flat stratum it is -- 588 of its 1,352 surface cells reading "emitted no
//! triangle" while 588 triangles were painted from a cell that was never
//! censused. Only a surface cell emits a triangle, so the candidate set is
//! enumerated -- at most two per axis, at most eight in all -- and the surface
//! cell among them is the host. That is a constraint from the extractor, not a
//! tie-break, and with it `box_exact` reads one colour and zero cells without a
//! triangle.
//!
//! Three counts on the HUD keep the rest of it honest: `unattributed`, triangles
//! no candidate could host, which is `0` on all sixteen rows; `in a shared face`,
//! triangles where **two** surface cells were candidates and either is a correct
//! answer -- measured 4 on `gyroid` at 33 samples, 62 on `noise_cavity` at 65,
//! and zero on five of the eight fields; and `cells emit none`, surface cells the
//! extractor produced no triangle for, so the picture cannot show their heat --
//! measured `0` on fourteen of the sixteen rows and 8 and 13 on `noise_cavity`.

mod common;

use bevy::platform::time::Instant;
use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{
    BoxExact, FbmTerrain, ReferenceField, Sphere, ThinPlate, Torus, capped_gyroid, csg_difference,
    noise_cavity,
};
use isomesh::marching_cubes::table::{AMBIGUOUS_FACES, is_inside};
use isomesh::marching_cubes::trilinear::BodySaddles;
use isomesh::marching_cubes::{FaceAmbiguity, InteriorAmbiguity, MarchingCubes};
use isomesh::{MeshBuffer, MeshSink, Real, RuntimeShape3, Sdf};

// ─── the two axes the demo moves along ──────────────────────────────────────

/// The eight reference fields, in the order the digit keys select them.
///
/// The seven most informative sit on `1`-`7`; `torus` is last because it reads
/// the same as `sphere` here -- every surface cell rank 2 -- and index 7 is
/// reachable only from `ISOMESH_FIELD`, which is the harness's contract
/// (`common/mod.rs:634-644` maps digits to `0..=6`).
const FIELD_COUNT: usize = 8;

/// Samples per axis, and the ladder `[` and `]` walk.
///
/// **These two rungs are not a taste.** P-130, P-131, P-133 and P-134 all have a
/// committed row at 33 and at 65 samples per axis and only there: P-130 stops at
/// 65, P-133 starts at 33, P-134 measured exactly this pair. A third rung would
/// leave at least one citation on the HUD without a row to quote.
const LADDER: [u32; 2] = [33, 65];

/// The bottom of the colour ramp, as `log10` of the normalised magnitude.
///
/// `1e-8` and `1e0` are the outer two of P-134's own `threshold_sweep` decades,
/// so the ramp's anchors are the thresholds the finding was reported against
/// rather than round numbers picked to look good.
const RAMP_LO: f64 = -8.0;

/// The top of the colour ramp, as `log10` of the normalised magnitude.
const RAMP_HI: f64 = 0.0;

/// P-134's `headline_threshold`: the decade the HUD counts cells above.
const HEADLINE_THRESHOLD: f64 = 1e-4;

/// Slack, in cells, on the check that a triangle lies inside the cell it was
/// attributed to.
///
/// The quantity checked is `(position - origin) / cell_size`, whose relative
/// error is a few ULP of at most `16 / 0.0625 = 256` on the widest domain here,
/// so about `1e-13`. `1e-6` is seven orders of margin above that and still eight
/// orders below "the neighbouring cell", which is the only thing the check has
/// to be able to tell apart.
const CELL_TOLERANCE: f64 = 1e-6;

/// The most cells of one class the overlay will cage.
///
/// `draw_normals` caps its gizmo lines for the reason this does: twelve lines
/// per cage means an uncapped class would let the frame-time readout measure the
/// overlay rather than the thing under test (`common/mod.rs:808-843`). It has
/// never bound on a reference field -- the two caged classes are at most 41 and
/// 33 cells on the whole sixteen-row census -- and the truncated count is
/// reported so a field that did bind it could not do so silently.
const MAX_CAGES: usize = 3_000;

// ─── the committed citations ────────────────────────────────────────────────

/// One committed row, quoted from the CSVs beside the live numbers.
///
/// The house rule is `game_dig.rs:2946-2952`: a figure that is not recomputed
/// live is written into the HUD **as a citation naming its finding id**. These
/// are those figures, and `the_citation_table_is_the_committed_csvs` re-reads all
/// four CSVs and checks every column of every row, so the transcription cannot
/// rot without a test going red.
struct Citation {
    /// `ReferenceField::NAME`, which is also the CSVs' `field` column.
    field: &'static str,
    /// Samples per axis, which is the CSVs' `resolution` column.
    samples: u32,
    /// `p-130.csv`: `delta_positive`, `delta_negative`, `delta_zero`,
    /// `ambiguous_cells`, over surface cells.
    partition: [u64; 4],
    /// `p-131.csv`: `a_zero_hits` and `discriminant_zero_hits`, over **every**
    /// cell rather than only the surface ones -- that is the population P-131
    /// censused, because `roots()` runs on every cell the mesher hands it.
    branch: [u64; 2],
    /// `p-133.csv`, the `f32` row: `sign_disagreements_f32`.
    f32_disagreements: u64,
    /// `p-134.csv`: `rank_correlation_with_hausdorff`.
    hausdorff_rho: f64,
    /// `p-134.csv`: `delta_magnitude_min`, `delta_magnitude_mean` and
    /// `delta_magnitude_max`, the normalised magnitude over the same surface
    /// cells this demo censuses.
    ///
    /// These three are quoted **and reproduced**: the live census computes the
    /// same quotient over the same lattice, so the panel puts the two rows on
    /// top of each other with a verdict, and
    /// `every_row_reproduces_p130s_partition_and_p131s_branch_counts` asserts
    /// the agreement on all sixteen rows.
    magnitude: [f64; 3],
    /// `p-134.csv`: `in_correlation_population`. `false` on `box_exact`, whose
    /// normalised magnitude has zero variance.
    in_population: bool,
}

/// P-134's `c2_bar`: the rank correlation a field had to clear.
const C2_BAR: f64 = 0.50;

/// P-134's `c2_fields_above_bar` at its `c2_resolution` of 65.
const C2_FIELDS_ABOVE_BAR: u32 = 3;

/// P-134's `c2_min_fields`: how many it needed.
const C2_MIN_FIELDS: u32 = 4;

/// P-127's `random_rational_trials`.
const P127_TRIALS: u64 = 3_481;

/// P-127's `f32_sign_disagreements` over those trials.
const P127_F32_DISAGREEMENTS: u64 = 14;

/// P-127's `terms_disc`, which is also its `terms_cayley`.
const P127_TERMS: u32 = 12;

/// P-127's `pencil_matches_total` out of `pencil_pairings_checked`.
const P127_PENCILS: u32 = 3;

/// Every committed row this demo can show, in `field * LADDER.len() + rung`
/// order so the lookup is arithmetic rather than a search that could miss.
///
/// A `static` rather than a `const` because [`Census::cited`] hands out a
/// `&'static Citation`: a `const` is inlined at each use, so borrowing one
/// borrows a temporary. `AMBIGUOUS_FACES` is a `static` for the same reason.
///
/// Sources, all committed: `docs/experiments/p-130.csv`, `p-131.csv`,
/// `p-133.csv`, `p-134.csv`.
static CITED: [Citation; FIELD_COUNT * LADDER.len()] = [
    Citation {
        field: "gyroid",
        samples: 33,
        partition: [4743, 497, 0, 27],
        branch: [0, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.689_276,
        magnitude: [3.444552358e-7, 7.436178186e-2, 4.493030175e0],
        in_population: true,
    },
    Citation {
        field: "gyroid",
        samples: 65,
        partition: [20392, 1040, 0, 132],
        branch: [0, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.742_686,
        magnitude: [9.585891404e-9, 2.854461852e-2, 3.255925886e0],
        in_population: true,
    },
    Citation {
        field: "noise_cavity",
        samples: 33,
        partition: [5018, 1158, 0, 502],
        branch: [6, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.475_842,
        magnitude: [2.238136486e-7, 2.477387172e-1, 6.424062620e0],
        in_population: true,
    },
    Citation {
        field: "noise_cavity",
        samples: 65,
        partition: [23338, 5037, 0, 567],
        branch: [30, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.599_755,
        magnitude: [3.219107767e-9, 1.068246551e-1, 4.852256740e0],
        in_population: true,
    },
    Citation {
        field: "csg_difference",
        samples: 33,
        partition: [201, 27, 1160, 0],
        branch: [27206, 17],
        f32_disagreements: 4,
        hausdorff_rho: 0.491_836,
        magnitude: [0.000000000e0, 1.554114470e-3, 1.149606044e-1],
        in_population: true,
    },
    Citation {
        field: "csg_difference",
        samples: 65,
        partition: [933, 59, 5022, 0],
        branch: [218_120, 41],
        f32_disagreements: 14,
        hausdorff_rho: 0.674_217,
        magnitude: [0.000000000e0, 9.842778897e-4, 2.454373847e-1],
        in_population: true,
    },
    Citation {
        field: "box_exact",
        samples: 33,
        partition: [0, 0, 1352, 0],
        branch: [28672, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.0,
        magnitude: [0.000000000e0, 0.000000000e0, 0.000000000e0],
        in_population: false,
    },
    Citation {
        field: "box_exact",
        samples: 65,
        partition: [0, 0, 5768, 0],
        branch: [229_376, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.0,
        magnitude: [0.000000000e0, 0.000000000e0, 0.000000000e0],
        in_population: false,
    },
    Citation {
        field: "thin_plate",
        samples: 33,
        partition: [8, 0, 504, 0],
        branch: [24576, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.071_906,
        magnitude: [0.000000000e0, 1.929012346e-4, 1.234567901e-2],
        in_population: true,
    },
    Citation {
        field: "thin_plate",
        samples: 65,
        partition: [8, 0, 2040, 0],
        branch: [196_608, 0],
        f32_disagreements: 0,
        hausdorff_rho: -0.926_920,
        magnitude: [0.000000000e0, 2.441406250e-4, 6.250000000e-2],
        in_population: true,
    },
    Citation {
        field: "fbm_terrain",
        samples: 33,
        partition: [1958, 0, 0, 30],
        branch: [28755, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.165_041,
        magnitude: [2.348812127e-7, 7.845582343e-2, 1.480288276e0],
        in_population: true,
    },
    Citation {
        field: "fbm_terrain",
        samples: 65,
        partition: [8413, 0, 0, 58],
        branch: [245_314, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.249_949,
        magnitude: [1.543015463e-9, 5.197687315e-2, 1.631874567e0],
        in_population: true,
    },
    Citation {
        field: "sphere",
        samples: 33,
        partition: [1160, 0, 0, 0],
        branch: [0, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.198_941,
        magnitude: [2.748432540e-6, 9.976194379e-4, 5.724661581e-3],
        in_population: true,
    },
    Citation {
        field: "sphere",
        samples: 65,
        partition: [4760, 0, 0, 0],
        branch: [0, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.186_149,
        magnitude: [3.825221637e-8, 2.482375455e-4, 1.917420416e-3],
        in_population: true,
    },
    Citation {
        field: "torus",
        samples: 33,
        partition: [1128, 0, 0, 0],
        branch: [0, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.257_079,
        magnitude: [9.629397257e-5, 1.736637197e-2, 1.635731727e-1],
        in_population: true,
    },
    Citation {
        field: "torus",
        samples: 65,
        partition: [4208, 0, 0, 0],
        branch: [0, 0],
        f32_disagreements: 0,
        hausdorff_rho: 0.211_601,
        magnitude: [6.248944416e-6, 4.942329936e-3, 3.105287130e-2],
        in_population: true,
    },
];

// ─── colour ─────────────────────────────────────────────────────────────────

/// The ramp's anchor colours in sRGB, one per decade from `1e-8` to `1e0`.
const RAMP: [[f32; 3]; 5] = [
    [0.07, 0.11, 0.42],
    [0.10, 0.55, 0.88],
    [0.25, 0.80, 0.35],
    [0.96, 0.84, 0.20],
    [0.94, 0.20, 0.10],
];

/// The `Delta = 0` stratum, drawn off the ramp because a logarithm has nowhere
/// to put a zero.
const ZERO_STRATUM: [f32; 3] = [0.82, 0.83, 0.86];

/// A triangle whose cell could not be established, or a normalisation that did
/// not produce a finite number.
///
/// Loud on purpose, and paired with an `error!`: a colour nobody chose is a
/// colour a reader would rationalise, so this one is a colour nobody would.
const REFUSED: [f32; 3] = [1.0, 0.0, 0.65];

/// sRGB as a human picks it into the linear RGBA [`Mesh::ATTRIBUTE_COLOR`]
/// wants. Feeding sRGB in raw renders it washed out (E-208).
fn linear(srgb: [f32; 3]) -> [f32; 4] {
    Color::srgb(srgb[0], srgb[1], srgb[2])
        .to_linear()
        .to_f32_array()
}

/// The ramp colour for one normalised magnitude, as linear RGBA.
///
/// A magnitude of `0.0` takes the ramp's **coldest** colour rather than the
/// [`ZERO_STRATUM`] grey, and the distinction is the point: `log10(0.0)` is
/// `-inf`, which clamps to the bottom of the ramp, and that is the right answer
/// for a cell whose `Delta` is non-zero but whose quotient underflowed. Only
/// [`cell_colour`] can reach the stratum, and only from the sign class -- so
/// grey on screen means one thing.
///
/// A non-finite magnitude gets [`REFUSED`] and is counted, since it means the
/// denominator under- or overflowed rather than that the cell is uninteresting.
fn ramp(magnitude: f64) -> [f32; 4] {
    if !magnitude.is_finite() {
        return linear(REFUSED);
    }
    let t = ((magnitude.log10() - RAMP_LO) / (RAMP_HI - RAMP_LO)).clamp(0.0, 1.0);
    let last = RAMP.len() - 1;
    let scaled = t * last as f64;
    let step = (scaled as usize).min(last - 1);
    let f = (scaled - step as f64) as f32;
    let (lo, hi) = (RAMP[step], RAMP[step + 1]);
    linear([
        lo[0] + (hi[0] - lo[0]) * f,
        lo[1] + (hi[1] - lo[1]) * f,
        lo[2] + (hi[2] - lo[2]) * f,
    ])
}

// ─── the grid ───────────────────────────────────────────────────────────────

/// The lattice the census and the extraction share.
///
/// One definition, so the two cannot disagree about which numbers they looked
/// at. The arithmetic is `sdf::sample_grid`'s at `sdf.rs:180-193`, which is what
/// makes the live counts comparable to `p-130.csv` and `p-131.csv` rather than
/// merely similar to them.
struct Grid {
    /// World position of sample `[0, 0, 0]`: the field's `domain()` minimum.
    origin: [f64; 3],
    /// `(max - min) / (samples - 1)`.
    cell_size: f64,
    /// Samples per axis.
    samples: u32,
    /// Cells per axis: `samples - 1`.
    cells: u32,
}

impl Grid {
    /// `i = x + y*n + z*n*n`, the crate's order, with `x` innermost.
    fn sample_index(&self, x: usize, y: usize, z: usize) -> usize {
        let n = self.samples as usize;
        x + y * n + z * n * n
    }

    /// The same order over cells rather than samples.
    fn cell_index(&self, cell: [u32; 3]) -> usize {
        let c = self.cells as usize;
        cell[0] as usize + cell[1] as usize * c + cell[2] as usize * c * c
    }

    /// Cells in the whole grid.
    fn cell_count(&self) -> u64 {
        u64::from(self.cells).pow(3)
    }

    /// The world position of one sample.
    fn point(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        [
            self.origin[0] + self.cell_size * x as f64,
            self.origin[1] + self.cell_size * y as f64,
            self.origin[2] + self.cell_size * z as f64,
        ]
    }

    /// The minimum corner of one cell, in world space.
    fn cell_corner(&self, cell: [u32; 3]) -> Vec3 {
        let p = self.point(cell[0] as usize, cell[1] as usize, cell[2] as usize);
        Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32)
    }

    /// Which cell emitted a triangle, as a linear cell index, and whether the
    /// triangle had more than one candidate.
    ///
    /// Every vertex a cell emits lies inside that cell's closed box -- on one of
    /// its twelve edges for a cached edge vertex, strictly inside for A-015's
    /// cycle centroids and A-002h's inner hexagon -- so a candidate is any cell
    /// `i` with `i <= lo` and `hi <= i + 1` on every axis, where `lo` and `hi`
    /// bound the triangle's span in cell units. Since the span is at most one
    /// cell wide there are at most two candidates per axis and so at most eight
    /// in all, and two on an axis only when the triangle lies exactly in the
    /// plane the pair shares.
    ///
    /// **The pair is not interchangeable, and assuming it was is what this
    /// signature exists to stop.** Picking the far side of the plane by
    /// `floor(midpoint)` looked right and was measured wrong on `box_exact`,
    /// which is planar on its own faces: the far cell is routinely one whose
    /// eight corners are all outside, which emitted nothing, and its whole
    /// surface came out in two colours instead of the one flat stratum it is.
    /// Only a surface cell emits a triangle, so among the candidates the surface
    /// cell is the one that did -- that is a constraint from the extractor and
    /// not a tie-break.
    ///
    /// `None` when no candidate is a surface cell, which is a measurement rather
    /// than a case to paper over: [`Census::unattributed`] carries the count and
    /// [`rebuild`] shouts. Where *two* candidates are surface cells the triangle
    /// lies in a face they share and either is a correct answer to "which cell
    /// emitted this"; [`Census::face_coplanar`] counts those.
    fn cell_of_triangle(
        &self,
        triangle: &[[f64; 3]; 3],
        heat: &[CellHeat],
    ) -> Option<(usize, bool)> {
        let last = f64::from(self.cells.saturating_sub(1));
        let mut span = [[0u32; 2]; 3];
        for (axis, slot) in span.iter_mut().enumerate() {
            let at = |v: &[f64; 3]| (v[axis] - self.origin[axis]) / self.cell_size;
            let mut lo = at(&triangle[0]);
            let mut hi = lo;
            for vertex in triangle.iter().skip(1) {
                let t = at(vertex);
                lo = lo.min(t);
                hi = hi.max(t);
            }
            // Smallest `i` with `hi <= i + 1`, and largest with `i <= lo`.
            let first = (hi - 1.0 - CELL_TOLERANCE).ceil().clamp(0.0, last);
            let final_ = (lo + CELL_TOLERANCE).floor().clamp(0.0, last);
            if first > final_ {
                return None;
            }
            *slot = [first as u32, final_ as u32];
        }

        let mut found = None;
        let mut candidates = 0u32;
        for z in span[2][0]..=span[2][1] {
            for y in span[1][0]..=span[1][1] {
                for x in span[0][0]..=span[0][1] {
                    let slot = self.cell_index([x, y, z]);
                    if heat.get(slot).is_some_and(|cell| cell.surface) {
                        candidates += 1;
                        found = found.or(Some(slot));
                    }
                }
            }
        }
        found.map(|slot| (slot, candidates > 1))
    }
}

// ─── per-cell state ─────────────────────────────────────────────────────────

/// What the heatmap needs to know about one cell.
#[derive(Clone, Copy)]
struct CellHeat {
    /// Linear RGBA for every triangle this cell emits.
    colour: [f32; 4],
    /// Triangles the extractor's output was attributed to this cell.
    triangles: u32,
    /// Whether the cell's corner signs are mixed.
    ///
    /// This is what makes [`Grid::cell_of_triangle`] able to choose between the
    /// two cells that share a face a coplanar triangle lies in: only a surface
    /// cell emits a triangle, so only a surface cell can be its host.
    surface: bool,
}

impl Default for CellHeat {
    /// A cell that is not a surface cell emits nothing, so this colour is
    /// unreachable -- `cell_of_triangle` will not attribute a triangle to a cell
    /// with `surface: false`. It is [`REFUSED`] rather than black so that
    /// reaching it would be visible instead of plausible.
    fn default() -> Self {
        Self {
            colour: linear(REFUSED),
            triangles: 0,
            surface: false,
        }
    }
}

/// `+1`, `-1` or `0`, so comparing two signs compares three states.
///
/// A `NaN` lands on `0` and is counted separately by
/// [`Census::non_finite_samples`]; without that split a `NaN` corner would make
/// every comparison below meaningless in a way no count would reveal.
fn sign_class<R: Real>(value: R) -> i8 {
    if value > R::ZERO {
        1
    } else if value < R::ZERO {
        -1
    } else {
        0
    }
}

// ─── resources ──────────────────────────────────────────────────────────────

/// Which field is showing.
#[derive(Resource, Default)]
struct Field(usize);

/// Which rung of [`LADDER`] is showing.
#[derive(Resource, Default)]
struct Rung(usize);

/// A field pinned by `ISOMESH_FIELD`, which takes it out of the capture walk.
#[derive(Resource)]
struct PinnedField(Option<usize>);

/// A rung pinned by `ISOMESH_SAMPLES`, which takes it out of the capture walk.
#[derive(Resource)]
struct PinnedRung(Option<usize>);

/// Which overlays are drawn.
#[derive(Resource)]
struct Show {
    /// The magenta cages: cells where the `f32` sign of `Delta` disagrees with
    /// the `f64` reference. On by default -- there are never many, and they are
    /// the point of P-133's C1.
    disagreements: bool,
    /// The cyan cages: cells where `a != 0` and `Delta == 0`, so
    /// `trilinear.rs:250` fires. Off by default because on most fields the class
    /// is empty and an always-on empty overlay teaches nothing.
    zero_branch: bool,
    /// The coloured surface itself. A cage can sit inside solid rock.
    surface: bool,
}

impl Default for Show {
    fn default() -> Self {
        Self {
            disagreements: true,
            zero_branch: false,
            surface: true,
        }
    }
}

/// Cell boxes to draw, resolved to world space once per rebuild.
#[derive(Resource, Default)]
struct Overlay {
    /// Cell size in world units, which is also the cage size.
    cell_size: f32,
    /// Minimum corners of the `f32` sign-disagreement cells.
    disagreements: Vec<Vec3>,
    /// Minimum corners of the cells where `trilinear.rs:250` fires.
    zero_branch: Vec<Vec3>,
}

/// Where the camera sits for the field on screen.
#[derive(Resource, Default)]
struct Framing {
    /// Domain centre.
    focus: Vec3,
    /// Orbit radius, from the field's own extent. `0.0` until the first rebuild.
    radius: f32,
}

/// The material every field's surface shares.
#[derive(Resource)]
struct SurfaceMaterial(Handle<StandardMaterial>);

/// Everything the census measured, and everything the HUD reads.
#[derive(Resource, Default)]
struct Census {
    /// `ReferenceField::NAME`.
    field_name: &'static str,
    /// Index into the field list, so the citation lookup is arithmetic.
    field: usize,
    /// Index into [`LADDER`].
    rung: usize,
    /// Samples per axis.
    samples: u32,
    /// Cells in the whole grid, which is P-131's population.
    grid_cells: u64,
    /// Cells whose case byte is neither `0x00` nor `0xFF`, which is P-130's.
    surface_cells: u64,
    /// Surface cells with `Delta > 0`: real tensor rank 2.
    delta_positive: u64,
    /// Surface cells with `Delta < 0`: real tensor rank 3.
    delta_negative: u64,
    /// Surface cells with `Delta == 0`: rank not fixed by the sign.
    delta_zero: u64,
    /// Surface cells with `AMBIGUOUS_FACES[case] != 0` -- the only cells on
    /// which the shipped `trilinear.rs:246` runs at all.
    ambiguous_cells: u64,
    /// Smallest normalised magnitude over surface cells.
    magnitude_min: f64,
    /// Mean normalised magnitude over surface cells.
    magnitude_mean: f64,
    /// Largest normalised magnitude over surface cells.
    magnitude_max: f64,
    /// Surface cells above [`HEADLINE_THRESHOLD`].
    above_headline: u64,
    /// Surface cells whose `max(abs(f_i))^4` did not produce a finite quotient.
    unnormalised: u64,
    /// Surface cells where the `f32` sign of `Delta` differs from the `f64` one.
    f32_disagreements: u64,
    /// Surface cells whose case byte differs between `f64` and `f32`.
    ///
    /// A coherence guard: a non-zero here means the two arms censused different
    /// cell sets and the disagreement count above is not about arithmetic.
    f32_case_flips: u64,
    /// Cells with `a == 0`, over the whole grid: `roots()` solves the linear
    /// polynomial and never reaches `:246`.
    branch_a_zero: u64,
    /// Cells with `a != 0` and `Delta == 0`, over the whole grid: the branch at
    /// `trilinear.rs:250` fires.
    branch_disc_zero: u64,
    /// How many of those are surface cells.
    branch_disc_zero_surface: u64,
    /// Non-finite samples. Anything but zero invalidates every count above.
    non_finite_samples: u64,
    /// Vertices in the extracted, edge-shared mesh.
    vertices: usize,
    /// Triangles in the extracted mesh.
    triangles: usize,
    /// Triangles whose cell could not be established.
    unattributed: u64,
    /// Triangles that lie exactly in a face two surface cells share, so both are
    /// correct answers to "which cell emitted this".
    ///
    /// Non-zero only on a field with axis-aligned planar geometry, and it is the
    /// count that says how much of the picture rests on that choice.
    face_coplanar: u64,
    /// Surface cells the extractor emitted no triangle for, so the coloured
    /// surface does not show their heat.
    cells_without_triangles: u64,
    /// Cages dropped by [`MAX_CAGES`].
    cages_dropped: u64,
    /// Extraction time, in milliseconds.
    extract_ms: f64,
    /// Census time, in milliseconds: sampling plus the per-cell pass.
    census_ms: f64,
}

impl Census {
    /// The committed row for what is on screen.
    fn cited(&self) -> Option<&'static Citation> {
        CITED.get(self.field * LADDER.len() + self.rung)
    }

    /// Whether the live partition reproduces `p-130.csv` exactly.
    fn matches_p130(&self, cited: &Citation) -> bool {
        cited.partition
            == [
                self.delta_positive,
                self.delta_negative,
                self.delta_zero,
                self.ambiguous_cells,
            ]
    }

    /// The worst relative gap between the live normalised magnitude and
    /// `p-134.csv`'s three columns.
    ///
    /// **Not zero, and that is not an error.** `p-134.csv` evaluates `Delta`
    /// through `benches/common::poly`'s **expanded** twelve-term Cayley form;
    /// this demo evaluates it through the shipped `BodySaddles::coefficients`,
    /// because the shipped route is its subject. P-127 proved those are the same
    /// polynomial *exactly*; P-130's header records that they are **not** the
    /// same `f64` computation and names its route for that reason. So a gap of
    /// order `1e-8` is the two routes rounding differently on the one cell that
    /// holds the extremum, and a gap of order `1e-1` would be a disagreement
    /// about the polynomial.
    ///
    /// A committed exact zero is compared exactly: `box_exact`'s whole surface
    /// and `csg_difference`'s minimum are hard zeros in both routes, and a live
    /// near-zero against one of those would be a different finding rather than a
    /// rounding of the same one.
    fn magnitude_gap(&self, cited: &Citation) -> f64 {
        [self.magnitude_min, self.magnitude_mean, self.magnitude_max]
            .iter()
            .zip(cited.magnitude.iter())
            .fold(0.0f64, |worst, (live, quoted)| {
                let gap = if *quoted == 0.0 {
                    if *live == 0.0 { 0.0 } else { f64::INFINITY }
                } else {
                    (*live - *quoted).abs() / quoted.abs()
                };
                worst.max(gap)
            })
    }
}

/// The most the live normalised magnitude may differ from `p-134.csv`'s, as a
/// relative gap.
///
/// The two differ only because the two evaluation routes round differently --
/// see [`Census::magnitude_gap`] -- and the observed worst over all sixteen
/// rows is `7.25e-8`, on `gyroid` at 65 samples, where the minimum sits on a
/// single cell. `1e-6` is more than an order above that and six orders below one
/// decade of the ramp, so a gap inside this bound cannot move a cell to a
/// different colour on screen, and a gap outside it is not a rounding.
const MAGNITUDE_GAP_BOUND: f64 = 1e-6;

/// Cages get their own group so they can be pulled in front of the opaque
/// surface without dragging the shared wireframe's width along with them.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct CellGizmos;

// ─── app ────────────────────────────────────────────────────────────────────

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-317 hyperdeterminant cells".into(),
                // Web only, inert on native: bind to the 1280x720 canvas the
                // page supplies rather than letting Bevy append its own. The HUD
                // panel is laid out in pixels for that size, so the canvas is
                // fixed and CSS scales it -- `fit_canvas_to_parent` stays at its
                // `false` default for the same reason.
                canvas: Some("#isomesh-canvas".into()),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<CellGizmos>()
        .insert_resource(PinnedField(pinned_field()))
        .insert_resource(PinnedRung(pinned_rung()))
        .init_resource::<Field>()
        .init_resource::<Rung>()
        .init_resource::<Show>()
        .init_resource::<Overlay>()
        .init_resource::<Framing>()
        .init_resource::<Census>()
        .add_systems(Startup, setup)
        // `PreUpdate` for E-306's reason: the harness's `update_hud` renders
        // `DemoStats` and its `capture_sequence` advances `Capture::taken`, both
        // in `Update` with no ordering against an example's own systems. In
        // `Update` the HUD would render a frame-old readout beside a current
        // caption, which on this demo means a cited partition sitting next to
        // the previous field's live one -- the exact comparison the panel exists
        // to make (E-312).
        // `.after(InputSystems)` because `controls` reads `just_pressed`, and
        // `bevy_input`'s own systems live in `PreUpdate` too: unordered, a
        // keypress would be read on the frame before it arrived about half the
        // time. `affine_rejection.rs:1211-1216` orders its `PreUpdate` chain the
        // same way.
        .add_systems(
            PreUpdate,
            (controls, rebuild, frame_camera, report)
                .chain()
                .after(bevy::input::InputSystems),
        )
        .add_systems(Update, draw_overlay)
        .run();
}

/// The field `ISOMESH_FIELD` asks for, if it asks for one this demo has.
fn pinned_field() -> Option<usize> {
    let raw = std::env::var("ISOMESH_FIELD").ok()?;
    match raw.trim().parse::<usize>() {
        Ok(index) if index < FIELD_COUNT => Some(index),
        _ => {
            error!("ISOMESH_FIELD={raw} is not one of 0..{FIELD_COUNT}");
            None
        }
    }
}

/// The rung `ISOMESH_SAMPLES` asks for, if it names one of [`LADDER`]'s.
///
/// Refused rather than clamped or rounded to the nearest rung: an off-ladder
/// resolution has no row in any of the four CSVs this demo quotes, so the panel
/// would show live numbers beside four blanks and the `EQUAL` gate would have
/// nothing to compare against.
fn pinned_rung() -> Option<usize> {
    let samples = common::samples_override()?;
    match LADDER.iter().position(|rung| *rung == samples) {
        Some(rung) => Some(rung),
        None => {
            error!(
                "ISOMESH_SAMPLES={samples} is not one of {LADDER:?}. P-130, P-131, P-133 and \
                 P-134 all have a committed row at 33 and at 65 samples per axis and at no \
                 other resolution the four share, so an off-ladder grid would leave every \
                 citation on the HUD without a row to quote."
            );
            None
        }
    }
}

/// Captured frames spent on each (field, rung) pair.
///
/// Read from the environment rather than from [`Capture`], which keeps its
/// length private, and the alternative is editing the harness. Sixteen stops at
/// `frames/16` each means `record_gif.sh`'s default 80 frames is exactly one
/// pass through the whole census and a short smoke capture still visits every
/// stop.
fn capture_frames_per_stop() -> u32 {
    let frames: u32 = std::env::var("ISOMESH_CAPTURE_FRAMES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(60);
    (frames / (FIELD_COUNT * LADDER.len()) as u32).max(1)
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut gizmo_config: ResMut<GizmoConfigStore>,
    mut camera: Query<&mut OrbitCamera>,
) {
    // Off both axes, so a cage reads as a box rather than as a square. A square
    // hides what is inside it, which is E-304's lesson about the same overlay.
    for mut orbit in &mut camera {
        orbit.yaw = 0.68;
        orbit.pitch = 0.34;
    }

    let (cells, _) = gizmo_config.config_mut::<CellGizmos>();
    cells.line.width = 1.8;
    // Harder than E-304's -0.4, because the surface here is opaque rather than
    // translucent and a cage inside it would otherwise be invisible exactly
    // where it matters.
    cells.depth_bias = -0.7;

    // `Color::WHITE`, so `StandardMaterial`'s multiply leaves the vertex colour
    // alone -- E-301's rule for a mesh whose whole content is its colour.
    // Double-sided because `gyroid` and `noise_cavity` are cave-like and a
    // reader looking into one would otherwise see through the far wall.
    commands.insert_resource(SurfaceMaterial(materials.add(StandardMaterial {
        base_color: Color::WHITE,
        perceptual_roughness: 0.85,
        metallic: 0.0,
        double_sided: true,
        cull_mode: None,
        ..default()
    })));

    // Spawned rather than assumed: `draw_domain` queries for it, and without one
    // the `G` toggle silently does nothing. Filled in by the first rebuild.
    commands.spawn(DemoDomain {
        min: Vec3::splat(-1.0),
        max: Vec3::splat(1.0),
    });
}

/// Field, resolution and the three overlay toggles.
///
/// Under capture the field and the rung both advance on the captured-frame
/// counter, because an example whose subject only changes on a keypress captures
/// as a still frame. A pinned field or rung wins over the walk, which is the
/// harness's contract: anything a committed capture depends on has to be
/// reachable from the environment.
#[allow(clippy::too_many_arguments)]
fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    pinned_field: Res<PinnedField>,
    pinned_rung: Res<PinnedRung>,
    mut field: ResMut<Field>,
    mut rung: ResMut<Rung>,
    mut show: ResMut<Show>,
) {
    let walking = capture.is_active();
    let stop = if walking {
        (capture.taken / capture_frames_per_stop()) as usize % (FIELD_COUNT * LADDER.len())
    } else {
        0
    };

    field.0 = match (pinned_field.0, walking) {
        (Some(index), _) => index,
        (None, true) => stop / LADDER.len(),
        (None, false) => flags.field.min(FIELD_COUNT - 1),
    };

    match (pinned_rung.0, walking) {
        (Some(index), _) => rung.0 = index,
        (None, true) => rung.0 = stop % LADDER.len(),
        (None, false) => {
            if keys.just_pressed(KeyCode::BracketRight) && rung.0 + 1 < LADDER.len() {
                rung.0 += 1;
            }
            if keys.just_pressed(KeyCode::BracketLeft) && rung.0 > 0 {
                rung.0 -= 1;
            }
        }
    }

    if keys.just_pressed(KeyCode::KeyX) {
        show.disagreements = !show.disagreements;
    }
    if keys.just_pressed(KeyCode::KeyZ) {
        show.zero_branch = !show.zero_branch;
    }
    if keys.just_pressed(KeyCode::KeyV) {
        show.surface = !show.surface;
    }
}

/// Census the grid, extract the surface, and colour it per cell.
#[allow(clippy::too_many_arguments)]
fn rebuild(
    field: Res<Field>,
    rung: Res<Rung>,
    flags: Res<ViewFlags>,
    mut census: ResMut<Census>,
    mut overlay: ResMut<Overlay>,
    mut framing: ResMut<Framing>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut query: Query<&mut Mesh3d, With<DemoMesh>>,
    mut domain: Query<&mut DemoDomain>,
    material: Res<SurfaceMaterial>,
    mut commands: Commands,
    mut last: Local<Option<(usize, usize)>>,
) {
    let key = (field.0, rung.0);
    if *last == Some(key) && !flags.remesh_requested {
        return;
    }
    *last = Some(key);

    let Some(built) = build(field.0, rung.0) else {
        return;
    };

    for mut edge in &mut domain {
        edge.min = built.domain_min;
        edge.max = built.domain_max;
    }
    *framing = built.framing;
    *overlay = built.overlay;

    // Loud, not decorative. Each of these means a number on the panel is not
    // measuring what its label says, and a colour alone would let a reader
    // rationalise it -- E-301's rule.
    if built.census.non_finite_samples > 0 {
        error!(
            "{} at {}^3: {} non-finite samples, so every sign comparison on the panel is \
             meaningless",
            built.census.field_name, built.census.samples, built.census.non_finite_samples
        );
    }
    if built.census.unattributed > 0 {
        error!(
            "{} at {}^3: {} triangles do not lie inside any single cell, so the heatmap is \
             colouring them from no cell's magnitude",
            built.census.field_name, built.census.samples, built.census.unattributed
        );
    }
    if built.census.f32_case_flips > 0 {
        error!(
            "{} at {}^3: {} cells have different case bytes at f64 and f32, so the two arms \
             censused different cell sets and the disagreement count is not about arithmetic",
            built.census.field_name, built.census.samples, built.census.f32_case_flips
        );
    }
    if built.census.unnormalised > 0 {
        error!(
            "{} at {}^3: {} cells have no finite abs(Delta)/max(abs(f_i))^4, so the \
             normalisation P-134 measured as exact did not hold here",
            built.census.field_name, built.census.samples, built.census.unnormalised
        );
    }

    // The HUD is the evidence and a headless capture has no HUD to read. One
    // line per rebuild, so `ISOMESH_CAPTURE` leaves the census in the log where
    // a script can hold it against p-130.csv -- E-203 learned this the hard
    // way: a measurement that only exists on screen cannot be verified from a
    // terminal.
    let cited = built.census.cited();
    info!(
        "{} at {}^3: {} surface cells of {} = {} pos + {} neg + {} zero, {} ambiguous; \
         p-130.csv {}; f32 sign disagreements {} ({} case flips); :250 fires on {} of {} cells \
         ({} on the surface), a == 0 on {}; magnitude {:.3e}..{:.3e}; census {:.1} ms, \
         extract {:.1} ms",
        built.census.field_name,
        built.census.samples,
        built.census.surface_cells,
        built.census.grid_cells,
        built.census.delta_positive,
        built.census.delta_negative,
        built.census.delta_zero,
        built.census.ambiguous_cells,
        match cited {
            Some(row) if built.census.matches_p130(row) => "EQUAL",
            Some(_) => "DIFFERS",
            None => "not cited",
        },
        built.census.f32_disagreements,
        built.census.f32_case_flips,
        built.census.branch_disc_zero,
        built.census.grid_cells,
        built.census.branch_disc_zero_surface,
        built.census.branch_a_zero,
        built.census.magnitude_min,
        built.census.magnitude_max,
        built.census.census_ms,
        built.census.extract_ms,
    );

    if let Some(row) = cited
        && !built.census.matches_p130(row)
    {
        error!(
            "{} at {}^3: the live partition is {}/{}/{}/{} and p-130.csv says {}/{}/{}/{}. \
             This file, that CSV or BodySaddles::coefficients has moved.",
            built.census.field_name,
            built.census.samples,
            built.census.delta_positive,
            built.census.delta_negative,
            built.census.delta_zero,
            built.census.ambiguous_cells,
            row.partition[0],
            row.partition[1],
            row.partition[2],
            row.partition[3],
        );
    }

    *census = built.census;

    // `Mesh3d::default()` names no asset and draws nothing, which is what an
    // empty extraction actually wants: adding an empty mesh to `Assets<Mesh>`
    // produces `Use-after-free: attempted to copy element data for an
    // unallocated key` twice a frame, forever (E-307).
    let handle = match built.mesh {
        Some(mesh) => Mesh3d(meshes.add(mesh)),
        None => Mesh3d::default(),
    };
    if query.is_empty() {
        commands.spawn((handle, MeshMaterial3d(material.0.clone()), DemoMesh));
    } else {
        for mut mesh in &mut query {
            *mesh = handle.clone();
        }
    }
}

/// Everything one (field, rung) pair produced.
struct Built {
    /// The coloured triangle soup, or `None` when the extraction was empty.
    mesh: Option<Mesh>,
    /// Every number the HUD reads.
    census: Census,
    /// The cell boxes to cage.
    overlay: Overlay,
    /// Where the camera goes.
    framing: Framing,
    /// The `G` box's minimum corner.
    domain_min: Vec3,
    /// The `G` box's maximum corner.
    domain_max: Vec3,
}

/// Dispatch on the field index, then do the work once in [`measure`].
///
/// The eight reference fields are eight different types, so
/// `for_each_reference_field!` cannot serve a runtime choice -- the index is
/// matched here instead, the same shape `critical_cells` and `manifold_check`
/// use.
fn build(field: usize, rung: usize) -> Option<Built> {
    let Some(&samples) = LADDER.get(rung) else {
        error!("rung {rung} is not one of 0..{}", LADDER.len());
        return None;
    };
    let Some(cited) = CITED.get(field * LADDER.len() + rung) else {
        error!("no committed citation for field {field} at rung {rung}");
        return None;
    };
    match field {
        0 => measure(&capped_gyroid::<f64>(), samples, field, rung, cited),
        1 => measure(&noise_cavity::<f64>(), samples, field, rung, cited),
        2 => measure(&csg_difference::<f64>(), samples, field, rung, cited),
        3 => measure(&BoxExact::<f64>::canonical(), samples, field, rung, cited),
        4 => measure(&ThinPlate::<f64>::canonical(), samples, field, rung, cited),
        5 => measure(&FbmTerrain::<f64>::canonical(), samples, field, rung, cited),
        6 => measure(&Sphere::<f64>::canonical(), samples, field, rung, cited),
        _ => measure(&Torus::<f64>::canonical(), samples, field, rung, cited),
    }
}

/// Sample the lattice, evaluate `Delta` on every cell, extract, colour.
///
/// The field is sampled twice -- once here and once inside the extractor -- and
/// that is deliberate: the census has to see the corner values as `f64` numbers
/// it can put through [`BodySaddles::coefficients`], and `MarchingCubes` owns
/// its own sample buffer. `critical_cells` pays the same price for the same
/// reason. The two passes cannot disagree, because [`Grid::point`] is
/// `sdf::sample_grid`'s arithmetic.
fn measure<F>(
    field: &F,
    samples: u32,
    field_index: usize,
    rung: usize,
    cited: &'static Citation,
) -> Option<Built>
where
    F: Sdf<Scalar = f64> + ReferenceField,
{
    if cited.field != F::NAME || cited.samples != samples {
        error!(
            "the citation at field {field_index} rung {rung} is {} at {}^3 and the field being \
             measured is {} at {samples}^3, so every quoted number on the panel belongs to a \
             different row",
            cited.field,
            cited.samples,
            F::NAME,
        );
    }

    let (min, max) = field.domain();
    let cells = samples.saturating_sub(1);
    if cells == 0 {
        error!("{samples} samples per axis leaves no cell to census");
        return None;
    }
    let grid = Grid {
        origin: min,
        cell_size: (max[0] - min[0]) / f64::from(cells),
        samples,
        cells,
    };

    let mut census = Census {
        field_name: F::NAME,
        field: field_index,
        rung,
        samples,
        grid_cells: grid.cell_count(),
        magnitude_min: f64::INFINITY,
        ..Census::default()
    };

    // ── sample, in `sample_grid`'s order and with its arithmetic ────────────
    let census_started = Instant::now();
    let n = samples as usize;
    let mut values = Vec::with_capacity(n * n * n);
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let v = field.sample(grid.point(x, y, z));
                if !v.is_finite() {
                    census.non_finite_samples += 1;
                }
                values.push(v);
            }
        }
    }

    // ── every cell: the branch census; every surface cell: the strata ───────
    let c = cells as usize;
    let mut heat = vec![CellHeat::default(); c * c * c];
    let mut surface_slots: Vec<u32> = Vec::new();
    let mut overlay = Overlay {
        cell_size: grid.cell_size as f32,
        ..Overlay::default()
    };
    let mut magnitude_total = 0.0f64;

    for cz in 0..c {
        for cy in 0..c {
            for cx in 0..c {
                // The crate's own gather: corner `k` at offset
                // `(k & 1, (k >> 1) & 1, (k >> 2) & 1)`, which is the indexing
                // `f[u + 2v + 4w]` P-127's polynomials are written in.
                let mut corner = [0.0f64; 8];
                let mut case = 0u8;
                for (k, slot) in corner.iter_mut().enumerate() {
                    let v = values
                        [grid.sample_index(cx + (k & 1), cy + ((k >> 1) & 1), cz + ((k >> 2) & 1))];
                    *slot = v;
                    if is_inside(v) {
                        case |= 1 << k;
                    }
                }

                // `4.0 * a * c` associates as `((4*a)*c)`, and so does
                // `R::TWO * R::TWO * a * c` at `trilinear.rs:246`. Same
                // groupings, same roundings, same number.
                let [a, b, quad_c] = BodySaddles::<f64>::coefficients(&corner);
                let delta = b * b - 4.0 * a * quad_c;
                let class = sign_class(delta);

                // P-131's population is every cell, because `roots()` runs on
                // every cell the mesher hands it -- so this half of the census
                // is taken before the surface filter, not after.
                let fires = if a == 0.0 {
                    census.branch_a_zero += 1;
                    false
                } else if class == 0 {
                    census.branch_disc_zero += 1;
                    true
                } else {
                    false
                };
                if fires {
                    if overlay.zero_branch.len() < MAX_CAGES {
                        overlay
                            .zero_branch
                            .push(grid.cell_corner([cx as u32, cy as u32, cz as u32]));
                    } else {
                        census.cages_dropped += 1;
                    }
                }

                if case == 0 || case == 255 {
                    continue;
                }
                census.surface_cells += 1;
                if fires {
                    census.branch_disc_zero_surface += 1;
                }
                match class {
                    1 => census.delta_positive += 1,
                    -1 => census.delta_negative += 1,
                    _ => census.delta_zero += 1,
                }
                if AMBIGUOUS_FACES[case as usize] != 0 {
                    census.ambiguous_cells += 1;
                }

                // `Delta` is homogeneous of degree 4, so this quotient is
                // scale-free exactly rather than approximately -- P-134's C1.
                //
                // `(M*M)*(M*M)` and not `M*M*M*M`, which is `p-134.csv`'s own
                // association and its stated reason: a *dyadic* rescaling of the
                // eight corners then leaves every mantissa in the computation
                // untouched, so "the normalisation is exact" is an exact-zero
                // prediction rather than a tolerance.
                let peak = corner.iter().fold(0.0f64, |m, v| m.max(v.abs()));
                let square = peak * peak;
                let magnitude = delta.abs() / (square * square);
                if magnitude.is_finite() {
                    census.magnitude_min = census.magnitude_min.min(magnitude);
                    census.magnitude_max = census.magnitude_max.max(magnitude);
                    magnitude_total += magnitude;
                    if magnitude > HEADLINE_THRESHOLD {
                        census.above_headline += 1;
                    }
                } else {
                    census.unnormalised += 1;
                }

                // The `f32` arm, on the same corner values narrowed. The case
                // byte is recomputed rather than reused, because a corner that
                // narrows to `-0.0` is outside by `is_inside`'s strict `< 0`
                // and would move the cell out of the population.
                let corner32: [f32; 8] = std::array::from_fn(|i| corner[i] as f32);
                let mut case32 = 0u8;
                for (k, v) in corner32.iter().enumerate() {
                    if is_inside(*v) {
                        case32 |= 1 << k;
                    }
                }
                if case32 != case {
                    census.f32_case_flips += 1;
                }
                let [a32, b32, quad_c32] = BodySaddles::<f32>::coefficients(&corner32);
                let delta32 = b32 * b32 - 4.0 * a32 * quad_c32;
                let disagrees = sign_class(delta32) != class;
                if disagrees {
                    census.f32_disagreements += 1;
                    if overlay.disagreements.len() < MAX_CAGES {
                        overlay
                            .disagreements
                            .push(grid.cell_corner([cx as u32, cy as u32, cz as u32]));
                    } else {
                        census.cages_dropped += 1;
                    }
                }

                let slot = grid.cell_index([cx as u32, cy as u32, cz as u32]);
                if let Some(cell) = heat.get_mut(slot) {
                    cell.colour = cell_colour(magnitude, class);
                    cell.surface = true;
                }
                surface_slots.push(slot as u32);
            }
        }
    }
    census.census_ms = census_started.elapsed().as_secs_f64() * 1000.0;
    census.magnitude_mean = if census.surface_cells == 0 {
        0.0
    } else {
        magnitude_total / census.surface_cells as f64
    };
    if census.surface_cells == 0 {
        census.magnitude_min = 0.0;
    }

    // ── the surface, in the one configuration that runs `:246` ──────────────
    let shape = match RuntimeShape3::new([samples; 3]) {
        Ok(shape) => shape,
        Err(error) => {
            error!("grid {samples}^3 rejected: {error}");
            return None;
        }
    };
    let mut mesher = MarchingCubes::<f64>::new();
    // Without both of these `trilinear.rs:246` never executes: the default is
    // `Separate` + `Ignore`, and `emit_trilinear` is gated on an ambiguous face
    // *and* on `Trilinear`. A demo about that line has to run it.
    mesher.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
    mesher.set_interior_ambiguity(InteriorAmbiguity::Trilinear);
    let mut buffer = MeshBuffer::<f64>::new();
    let extract_started = Instant::now();
    if let Err(error) = mesher.extract(field, &shape, grid.origin, grid.cell_size, &mut buffer) {
        error!(
            "marching cubes failed on {} at {samples}^3: {error}",
            F::NAME
        );
        return None;
    }
    census.extract_ms = extract_started.elapsed().as_secs_f64() * 1000.0;
    census.vertices = buffer.vertex_count();
    census.triangles = buffer.triangle_count();

    // ── the soup: three vertices per triangle, carrying its cell's colour ───
    let mut builder = MeshBuilder::new();
    for triangle in buffer.indices.as_chunks::<3>().0 {
        let (Some(p0), Some(p1), Some(p2)) = (
            buffer.positions.get(triangle[0] as usize),
            buffer.positions.get(triangle[1] as usize),
            buffer.positions.get(triangle[2] as usize),
        ) else {
            continue;
        };
        let (Some(n0), Some(n1), Some(n2)) = (
            buffer.normals.get(triangle[0] as usize),
            buffer.normals.get(triangle[1] as usize),
            buffer.normals.get(triangle[2] as usize),
        ) else {
            continue;
        };
        let corners = [*p0, *p1, *p2];
        let colour = match grid.cell_of_triangle(&corners, &heat) {
            Some((slot, shared)) => {
                if shared {
                    census.face_coplanar += 1;
                }
                match heat.get_mut(slot) {
                    Some(cell) => {
                        cell.triangles += 1;
                        cell.colour
                    }
                    None => {
                        census.unattributed += 1;
                        linear(REFUSED)
                    }
                }
            }
            None => {
                census.unattributed += 1;
                linear(REFUSED)
            }
        };
        let a = builder.vertex(narrow(*p0), narrow(*n0));
        let b = builder.vertex(narrow(*p1), narrow(*n1));
        let c = builder.vertex(narrow(*p2), narrow(*n2));
        builder.triangle(a, b, c);
        let colours = builder.colors_mut();
        colours.push(colour);
        colours.push(colour);
        colours.push(colour);
    }

    census.cells_without_triangles = surface_slots
        .iter()
        .filter(|slot| heat.get(**slot as usize).is_some_and(|c| c.triangles == 0))
        .count() as u64;

    // Read before the builder is moved into `into_mesh`, and the reason it is a
    // question at all is E-307's: adding an empty mesh to `Assets<Mesh>`
    // produces `Use-after-free: attempted to copy element data for an
    // unallocated key` twice a frame, forever.
    let has_geometry = builder.vertex_count() > 0;

    // ── framing, from the field's own domain ────────────────────────────────
    //
    // Never a hardcoded radius: the eight domains differ by 4x -- five are
    // half-extent 2, the capped gyroid is 7 and `fbm_terrain` is 8 -- and a
    // fixed number puts the camera comfortably *inside* the gyroid. E-110 found
    // E-109's committed screenshot was a picture of an inner wall for exactly
    // that reason.
    let domain_min = Vec3::new(min[0] as f32, min[1] as f32, min[2] as f32);
    let domain_max = Vec3::new(max[0] as f32, max[1] as f32, max[2] as f32);
    let extent = domain_max.x - domain_min.x;

    Some(Built {
        mesh: has_geometry.then(|| builder.into_mesh()),
        census,
        overlay,
        framing: Framing {
            focus: (domain_min + domain_max) * 0.5,
            radius: extent * VIEW_EXTENTS,
        },
        domain_min,
        domain_max,
    })
}

/// `f64` triple to `f32`, for the [`Mesh`] the picture is drawn from.
///
/// Cast rather than re-extracted in `f32`: the numbers on the panel are `f64`
/// numbers, so the mesh they are painted onto has to be the one they were
/// computed on.
fn narrow(v: [f64; 3]) -> [f32; 3] {
    [v[0] as f32, v[1] as f32, v[2] as f32]
}

/// The cell's colour: the ramp for `Delta != 0`, the flat stratum otherwise.
///
/// The sign class decides, not the magnitude: a `Delta` of exactly zero and a
/// `Delta` that rounded to a magnitude of zero are different facts, and only the
/// first is the stratum P-131 is about.
fn cell_colour(magnitude: f64, class: i8) -> [f32; 4] {
    if class == 0 {
        linear(ZERO_STRATUM)
    } else {
        ramp(magnitude)
    }
}

/// Orbit radius as a multiple of the domain's extent.
///
/// The whole surface has to be in frame -- a heatmap read on one face is not a
/// heatmap -- and this is the closest that keeps every field's silhouette inside
/// a 1280x720 viewport.
const VIEW_EXTENTS: f32 = 1.35;

/// Where the subject sits in frame, as a fraction of the orbit radius, right and
/// down from the centre.
///
/// **The HUD is twenty-odd lines in the upper left and the coloured surface is
/// the subject.** Centring it photographs the argument with its evidence hidden,
/// which is E-112's lesson and E-109's committed screenshot. Applied in the
/// camera's own basis rather than as a world offset, so it holds while
/// `ISOMESH_SPIN` yaws.
const SUBJECT_OFFSET: Vec2 = Vec2::new(0.22, 0.10);

/// Point the orbit camera at the field on screen.
///
/// The radius is written **only when the field changes**, so scroll-zoom
/// survives a resolution change and a re-mesh; the focus is written every frame
/// from the current radius, so the off-centre nudge stays one screen-space
/// offset however far the spin has turned.
fn frame_camera(
    field: Res<Field>,
    framing: Res<Framing>,
    mut camera: Query<&mut OrbitCamera>,
    mut last: Local<Option<usize>>,
) {
    if framing.radius <= 0.0 {
        return;
    }
    let moved = *last != Some(field.0);
    *last = Some(field.0);

    for mut orbit in &mut camera {
        if moved {
            orbit.radius = framing.radius;
        }
        // The camera's own basis, from the same yaw/pitch the harness's
        // `orbit_camera` builds its transform from. It places the eye at
        // `focus + dir * radius`, so the view direction is `-dir` and a focus
        // moved along `-right` puts the subject right of centre.
        let dir = Vec3::new(
            orbit.yaw.cos() * orbit.pitch.cos(),
            orbit.pitch.sin(),
            orbit.yaw.sin() * orbit.pitch.cos(),
        );
        let forward = -dir;
        let right = forward.cross(Vec3::Y).normalize_or_zero();
        let up = right.cross(forward).normalize_or_zero();
        orbit.focus = framing.focus - right * (SUBJECT_OFFSET.x * orbit.radius)
            + up * (SUBJECT_OFFSET.y * orbit.radius);
    }
}

// ─── the panel ──────────────────────────────────────────────────────────────

/// Write the census, and the committed rows it is meant to be read against, to
/// [`DemoStats`].
///
/// Separate from [`rebuild`] so it can be driven as a one-shot system by a test
/// on a machine with no display, which is the only way to see this demo's screen
/// here.
fn report(census: Res<Census>, show: Res<Show>, mut stats: ResMut<DemoStats>) {
    stats.title = format!(
        "E-317  hyperdeterminant cells - {} at {}^3",
        census.field_name, census.samples
    );
    stats.vertices = census.vertices;
    stats.triangles = census.triangles;
    stats.extract_ms = census.extract_ms;
    stats.keys = Some(
        "[1-7] field   [ ] resolution   [X] f32 cages   [Z] Delta=0 cages   [V] surface\n\
         [W] wire   [N] normals   [G] grid   [Space] pause   [R] re-mesh   [H] HUD   [Esc] quit"
            .into(),
    );
    stats.hint = Some("[H] HUD   [X] f32 cages   [Z] Delta=0 cages   [V] surface".into());

    let Some(cited) = census.cited() else {
        stats.banner = Some((
            format!(
                "no committed row for field {} at rung {}",
                census.field, census.rung
            ),
            Color::srgb(1.0, 0.45, 0.25),
        ));
        stats.extra = vec![
            "Every number this panel is built to compare comes from a committed CSV row, and"
                .into(),
            "there is none for what is on screen. Nothing below would be a comparison.".into(),
        ];
        return;
    };

    let equal = census.matches_p130(cited);
    stats.banner = Some((
        format!(
            "P-127: b*b-4ac IS Det(2,2,2)   |   strata {} / {} / {}   {} p-130.csv",
            census.delta_positive,
            census.delta_negative,
            census.delta_zero,
            if equal { "EQUAL to" } else { "DIFFER from" }
        ),
        if equal {
            Color::srgb(0.45, 0.95, 0.55)
        } else {
            Color::srgb(1.0, 0.35, 0.30)
        },
    ));

    let share = |n: u64| {
        if census.surface_cells == 0 {
            0.0
        } else {
            100.0 * n as f64 / census.surface_cells as f64
        }
    };
    let gap = census.magnitude_gap(cited);

    stats.extra = vec![
        format!(
            "{:>9} grid cells   {:>8} surface cells   {:.1} ms census",
            census.grid_cells, census.surface_cells, census.census_ms
        ),
        String::new(),
        " stratum            live     p-130     real tensor rank".into(),
        format!(
            " Delta > 0      {:>8}  {:>8}     2   (de Silva & Lim s6)",
            census.delta_positive, cited.partition[0]
        ),
        format!(
            " Delta < 0      {:>8}  {:>8}     3",
            census.delta_negative, cited.partition[1]
        ),
        format!(
            " Delta = 0      {:>8}  {:>8}     not fixed by the sign alone",
            census.delta_zero, cited.partition[2]
        ),
        format!(
            " ambiguous face {:>8}  {:>8}     the only cells :246 runs on",
            census.ambiguous_cells, cited.partition[3]
        ),
        format!(
            "  {:>7.3}% of surface cells are degenerate       {}",
            share(census.delta_zero),
            if equal {
                "EQUAL <- P-130's own lattice"
            } else {
                "DIFFERS <- see the log"
            }
        ),
        String::new(),
        format!(
            " abs(Delta)/max(abs(f_i))^4  min {:>9.3e}  mean {:>9.3e}  max {:>9.3e}",
            census.magnitude_min, census.magnitude_mean, census.magnitude_max
        ),
        format!(
            "  P-134 measured             min {:>9.3e}  mean {:>9.3e}  max {:>9.3e}",
            cited.magnitude[0], cited.magnitude[1], cited.magnitude[2],
        ),
        format!(
            "  worst relative gap {:>9.2e}   {}   two evaluation routes, one polynomial",
            gap,
            if gap <= MAGNITUDE_GAP_BOUND {
                "within 1e-6"
            } else {
                "OUTSIDE 1e-6"
            }
        ),
        format!(
            "  {:>8} of {} above 1e-4    degree-4 homogeneous, so exactly scale-free",
            census.above_headline, census.surface_cells
        ),
        format!(
            "  P-134  rho(magnitude, Hausdorff) {:>7.3} vs a {:.2} bar   in population: {}",
            cited.hausdorff_rho,
            C2_BAR,
            if cited.in_population {
                "yes"
            } else {
                "no, zero variance"
            }
        ),
        format!(
            "         C2 FALSIFIED: {C2_FIELDS_ABOVE_BAR} of the {C2_MIN_FIELDS} fields it \
             needed, box_exact excluded"
        ),
        String::new(),
        format!(
            " f32 sign vs the f64 reference {:>8} of {}   {} case flips",
            census.f32_disagreements, census.surface_cells, census.f32_case_flips
        ),
        format!(
            "  P-133 measured {:>7} here; it samples at f32, so not these corner values",
            cited.f32_disagreements
        ),
        format!(
            "  P-127  {P127_TERMS} terms each side, symbolic difference 0, \
             {P127_PENCILS}/{P127_PENCILS} pencil pairings;"
        ),
        format!(
            "         {P127_F32_DISAGREEMENTS} of {P127_TRIALS} exact-rational 8-tuples flip \
             sign at f32"
        ),
        String::new(),
        format!(
            "{:<36}a == 0 {:>9}   :250 fires {:>6}",
            format!(" Delta = 0, over all {} cells:", census.grid_cells),
            census.branch_a_zero,
            census.branch_disc_zero
        ),
        format!(
            "{:<36}a == 0 {:>9}   :250 fires {:>6}",
            "  P-131 measured", cited.branch[0], cited.branch[1]
        ),
        format!(
            "         {} of those are surface cells; border rank 2 where it fires,",
            census.branch_disc_zero_surface
        ),
        "         not the tangential touch trilinear.rs:251 describes".into(),
        String::new(),
        format!(
            " soup {:>8} tris   {} unattributed   {} in a shared face   {} cells emit none",
            census.triangles,
            census.unattributed,
            census.face_coplanar,
            census.cells_without_triangles,
        ),
        format!(
            "         cages: {}   {}",
            if show.disagreements {
                "X f32 disagreements on"
            } else {
                "X f32 disagreements off"
            },
            if show.zero_branch {
                "Z the :250 branch on"
            } else {
                "Z the :250 branch off"
            },
        ),
    ];

    if census.cages_dropped > 0 {
        stats.extra.push(format!(
            "         {} cages dropped by the {MAX_CAGES}-per-class overlay cap",
            census.cages_dropped
        ));
    }
}

/// Draw the cages, and hide the surface when `V` says so.
///
/// Runs every frame; the overlay changes only on rebuild.
fn draw_overlay(
    overlay: Res<Overlay>,
    show: Res<Show>,
    mut visibility: Query<&mut Visibility, With<DemoMesh>>,
    mut gizmos: Gizmos<CellGizmos>,
) {
    /// The `f32` sign-disagreement cages.
    const MAGENTA: Color = Color::srgb(1.0, 0.25, 0.90);
    /// The cages where `trilinear.rs:250` fires.
    const CYAN: Color = Color::srgb(0.15, 0.85, 1.0);

    // Written only when it differs. A `*visible = ...` every frame marks the
    // component changed every frame, and Bevy's visibility propagation is
    // change-driven -- so an unconditional write turns a toggle nobody pressed
    // into per-frame work on every descendant.
    let wanted = if show.surface {
        Visibility::Inherited
    } else {
        Visibility::Hidden
    };
    for mut visible in &mut visibility {
        if *visible != wanted {
            *visible = wanted;
        }
    }

    if show.disagreements {
        for corner in &overlay.disagreements {
            cage(&mut gizmos, *corner, overlay.cell_size, MAGENTA);
        }
    }
    if show.zero_branch {
        for corner in &overlay.zero_branch {
            cage(&mut gizmos, *corner, overlay.cell_size, CYAN);
        }
    }
}

/// The twelve edges of one cell, at its exact bounds.
///
/// Exact rather than inflated: a cage larger than its cell would make every cell
/// look like the one next to it, which is the one thing this picture must not
/// do. Corner indexing matches the extractor's -- bit `i` of the corner index is
/// axis `i`, the same convention `common::draw_domain` uses.
fn cage(gizmos: &mut Gizmos<CellGizmos>, min: Vec3, size: f32, colour: Color) {
    let corner = |i: usize| {
        min + Vec3::new(
            if i & 1 == 0 { 0.0 } else { size },
            if i & 2 == 0 { 0.0 } else { size },
            if i & 4 == 0 { 0.0 } else { size },
        )
    };
    for i in 0..8usize {
        for axis in 0..3usize {
            let bit = 1 << axis;
            if i & bit == 0 {
                gizmos.line(corner(i), corner(i | bit), colour);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use core::time::Duration;

    use bevy::asset::AssetApp;
    use bevy::ecs::system::RunSystemOnce;
    use bevy::mesh::VertexAttributeValues;
    use bevy::time::TimeUpdateStrategy;

    use super::*;

    /// A frame, fixed, so nothing here depends on how long the test machine
    /// took.
    const FRAME: Duration = Duration::from_millis(16);

    /// `p-130.csv`, read at compile time so the transcription in [`CITED`]
    /// cannot drift from the artefact it quotes.
    const P130_CSV: &str = include_str!("../../docs/experiments/p-130.csv");

    /// `p-131.csv`, likewise.
    const P131_CSV: &str = include_str!("../../docs/experiments/p-131.csv");

    /// `p-133.csv`, likewise.
    const P133_CSV: &str = include_str!("../../docs/experiments/p-133.csv");

    /// `p-134.csv`, likewise.
    const P134_CSV: &str = include_str!("../../docs/experiments/p-134.csv");

    /// One CSV as a header row plus data rows, comment lines dropped.
    ///
    /// The experiment CSVs carry the hypothesis, the falsifier and the
    /// provenance as `#` lines above the header, which is why this cannot be a
    /// `lines().next()`.
    fn table(csv: &str) -> (Vec<&str>, Vec<Vec<&str>>) {
        let mut rows = csv.lines().filter(|l| !l.starts_with('#'));
        let header: Vec<&str> = rows
            .next()
            .expect("the CSV has a header row")
            .split(',')
            .collect();
        let data = rows
            .filter(|l| !l.trim().is_empty())
            .map(|l| l.split(',').collect())
            .collect();
        (header, data)
    }

    /// Whether a transcribed constant is the number the CSV holds.
    ///
    /// Not `==`: the artefact holds a decimal rendering of an `f64` and the
    /// constant is a decimal literal, so a bit-exact comparison would be a test
    /// of two parsers agreeing. Zero is compared exactly, because a committed
    /// `0.000000000e0` came from a hard zero -- `box_exact`'s whole surface --
    /// and a near-zero beside it would be a different figure rather than a
    /// rounding of the same one.
    fn reproduces(live: f64, cited: f64) -> bool {
        if cited == 0.0 {
            live == 0.0
        } else {
            (live - cited).abs() <= 1e-9 * cited.abs()
        }
    }

    /// One cell of one CSV, by column name, from the first row matching every
    /// `(column, value)` in `where_`.
    fn cell(csv: &str, where_: &[(&str, &str)], column: &str) -> String {
        let (header, rows) = table(csv);
        let at = |name: &str| {
            header
                .iter()
                .position(|h| *h == name)
                .unwrap_or_else(|| panic!("{name} is not a column of this CSV"))
        };
        let wanted = at(column);
        for row in &rows {
            if where_
                .iter()
                .all(|(k, v)| row.get(at(k)).is_some_and(|got| got == v))
            {
                return row
                    .get(wanted)
                    .unwrap_or_else(|| panic!("row is short of {column}"))
                    .to_string();
            }
        }
        panic!("no row matching {where_:?}");
    }

    /// The demo's headless app: its own systems, no window and no renderer.
    ///
    /// `setup` is left out because it wants an `Assets<StandardMaterial>` and a
    /// `GizmoConfigStore`; the one thing it produces that `rebuild` needs is the
    /// surface material handle, which is inserted here by hand. `report` is left
    /// out too and run as a one-shot below, which is the same system the demo
    /// runs every frame.
    ///
    /// No stall-detecting drain: this demo's rebuild is synchronous, so one
    /// stepped frame produces the whole census. There is no queue to drain and
    /// no command flush to wait on -- the census is a resource written in
    /// `PreUpdate`, not geometry published through `Commands`.
    fn harness(field: usize, rung: usize) -> App {
        let mut app = App::new();
        app.add_plugins(MinimalPlugins)
            .add_plugins(bevy::asset::AssetPlugin::default())
            .init_asset::<Mesh>()
            .insert_resource(TimeUpdateStrategy::ManualDuration(FRAME))
            .init_resource::<ButtonInput<KeyCode>>()
            .init_resource::<ViewFlags>()
            .insert_resource(Capture::default())
            // Pinned rather than left to `ViewFlags::field`, which reads
            // `ISOMESH_FIELD` from the environment and would otherwise decide
            // what this test measures. `unsafe_code = "forbid"` means a test
            // here cannot set that variable, so it has to out-rank it.
            .insert_resource(PinnedField(Some(field)))
            .insert_resource(PinnedRung(Some(rung)))
            .insert_resource(Field(field))
            .insert_resource(Rung(rung))
            .init_resource::<Show>()
            .init_resource::<Overlay>()
            .init_resource::<Framing>()
            .init_resource::<Census>()
            .init_resource::<DemoStats>()
            .insert_resource(SurfaceMaterial(Handle::default()))
            .add_systems(Update, (controls, rebuild).chain());
        app
    }

    /// One frame, with the input clearing `InputPlugin` would have done.
    fn step(app: &mut App) {
        app.update();
        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .clear();
    }

    /// Run the demo once and hand back the panel a reader would read.
    fn panel(field: usize, rung: usize) -> (String, Vec<String>, Option<(String, Color)>) {
        let mut app = harness(field, rung);
        step(&mut app);
        app.world_mut()
            .run_system_once(report)
            .expect("the HUD system");
        let stats = app.world().resource::<DemoStats>();
        (
            stats.title.clone(),
            stats.extra.clone(),
            stats.banner.clone(),
        )
    }

    /// Find the one panel line containing `needle`.
    fn line<'a>(lines: &'a [String], needle: &str) -> &'a str {
        lines
            .iter()
            .find(|l| l.contains(needle))
            .map(String::as_str)
            .unwrap_or_else(|| panic!("no panel line contains {needle:?}"))
    }

    /// The transcription in [`CITED`] is what the four CSVs say.
    ///
    /// **This is the test that makes every quoted number on the panel a citation
    /// rather than a remembered figure.** The `include_str!`s above are
    /// `#[cfg(test)]` so the shipped example carries the house-style transcribed
    /// constants (`game_dig.rs:2946-2952`) while the artefacts still get the
    /// last word: rename or re-run a CSV and this goes red before anything can
    /// quote a figure that is no longer in it.
    #[test]
    fn the_citation_table_is_the_committed_csvs() {
        for cited in &CITED {
            let res = cited.samples.to_string();
            let key: [(&str, &str); 2] = [("field", cited.field), ("resolution", res.as_str())];

            for (column, quoted) in [
                ("delta_positive", cited.partition[0]),
                ("delta_negative", cited.partition[1]),
                ("delta_zero", cited.partition[2]),
                ("ambiguous_cells", cited.partition[3]),
            ] {
                assert_eq!(
                    cell(P130_CSV, &key, column),
                    quoted.to_string(),
                    "p-130.csv {column} for {} at {}^3",
                    cited.field,
                    cited.samples
                );
            }

            for (column, quoted) in [
                ("a_zero_hits", cited.branch[0]),
                ("discriminant_zero_hits", cited.branch[1]),
            ] {
                assert_eq!(
                    cell(P131_CSV, &key, column),
                    quoted.to_string(),
                    "p-131.csv {column} for {} at {}^3",
                    cited.field,
                    cited.samples
                );
            }

            let long = format!("{res}x{res}x{res}");
            let p133_key: [(&str, &str); 3] = [
                ("field", cited.field),
                ("resolution", long.as_str()),
                ("scalar", "f32"),
            ];
            assert_eq!(
                cell(P133_CSV, &p133_key, "sign_disagreements_f32"),
                cited.f32_disagreements.to_string(),
                "p-133.csv sign_disagreements_f32 for {} at {}^3",
                cited.field,
                cited.samples
            );

            let rho: f64 = cell(P134_CSV, &key, "rank_correlation_with_hausdorff")
                .parse()
                .expect("p-134.csv rank_correlation_with_hausdorff is a number");
            assert!(
                (rho - cited.hausdorff_rho).abs() < 5e-7,
                "p-134.csv rank_correlation_with_hausdorff for {} at {}^3 is {rho}, quoted {}",
                cited.field,
                cited.samples,
                cited.hausdorff_rho
            );
            assert_eq!(
                cell(P134_CSV, &key, "in_correlation_population"),
                cited.in_population.to_string(),
                "p-134.csv in_correlation_population for {} at {}^3",
                cited.field,
                cited.samples
            );

            // The three magnitude columns, parsed rather than string-compared:
            // the panel reproduces them, so they have to be the same *numbers*
            // and not the same rendering of them.
            for (column, quoted) in [
                ("delta_magnitude_min", cited.magnitude[0]),
                ("delta_magnitude_mean", cited.magnitude[1]),
                ("delta_magnitude_max", cited.magnitude[2]),
            ] {
                let raw = cell(P134_CSV, &key, column);
                let got: f64 = raw
                    .parse()
                    .unwrap_or_else(|_| panic!("p-134.csv {column} is not a number: {raw}"));
                assert!(
                    reproduces(quoted, got),
                    "p-134.csv {column} for {} at {}^3 is {got:e}, quoted {quoted:e}",
                    cited.field,
                    cited.samples
                );
            }
        }

        // P-127's globals, which are the same on all three of its rows.
        let any: [(&str, &str); 0] = [];
        let p127 = include_str!("../../docs/experiments/p-127.csv");
        assert_eq!(cell(p127, &any, "terms_disc"), P127_TERMS.to_string());
        assert_eq!(cell(p127, &any, "terms_cayley"), P127_TERMS.to_string());
        assert_eq!(cell(p127, &any, "symbolic_difference_is_zero"), "true");
        assert_eq!(
            cell(p127, &any, "pencil_matches_total"),
            P127_PENCILS.to_string()
        );
        assert_eq!(
            cell(p127, &any, "pencil_pairings_checked"),
            P127_PENCILS.to_string()
        );
        assert_eq!(
            cell(p127, &any, "random_rational_trials"),
            P127_TRIALS.to_string()
        );
        assert_eq!(
            cell(p127, &any, "f32_sign_disagreements"),
            P127_F32_DISAGREEMENTS.to_string()
        );

        // P-134's C2 verdict, which is the negative result the header states.
        let bar: [(&str, &str); 2] = [("field", "gyroid"), ("resolution", "65")];
        assert_eq!(cell(P134_CSV, &bar, "c2_holds"), "false");
        assert_eq!(
            cell(P134_CSV, &bar, "c2_fields_above_bar"),
            C2_FIELDS_ABOVE_BAR.to_string()
        );
        assert_eq!(
            cell(P134_CSV, &bar, "c2_min_fields"),
            C2_MIN_FIELDS.to_string()
        );

        // P-130's C1 and C2, likewise.
        assert_eq!(cell(P130_CSV, &bar, "c1_holds"), "false");
        assert_eq!(cell(P130_CSV, &bar, "c2_holds"), "false");
        assert_eq!(cell(P130_CSV, &bar, "c2_unreachable"), "true");
    }

    /// [`CITED`] is indexed arithmetically, so its order has to be the demo's.
    #[test]
    fn the_citation_table_is_in_the_order_the_lookup_assumes() {
        for (rung, samples) in LADDER.iter().enumerate() {
            for field in 0..FIELD_COUNT {
                let cited = CITED
                    .get(field * LADDER.len() + rung)
                    .expect("a row per field and rung");
                assert_eq!(
                    cited.samples,
                    *samples,
                    "row {} of CITED is at {} samples, not {samples}",
                    field * LADDER.len() + rung,
                    cited.samples
                );
            }
        }
        // Two rows per field, adjacent, and eight distinct field names.
        let mut names: Vec<&str> = CITED
            .iter()
            .step_by(LADDER.len())
            .map(|c| c.field)
            .collect();
        assert_eq!(names.len(), FIELD_COUNT);
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), FIELD_COUNT, "CITED repeats a field name");
    }

    /// The panel says what the census found, beside what P-130 measured.
    ///
    /// **This is the only way to see this demo's screen on a machine with no
    /// display**, and every line it checks is a line a reader is meant to read.
    /// `gyroid` at 33 samples per axis is the row with both signs present and a
    /// non-zero ambiguous-face count, so it exercises the whole panel rather
    /// than a degenerate corner of it.
    #[test]
    fn the_panel_reports_the_strata_beside_p130s_partition() {
        let (title, lines, banner) = panel(0, 0);
        for line in &lines {
            println!("{line}");
        }

        assert!(
            title.contains("E-317") && title.contains("gyroid") && title.contains("33^3"),
            "the title lost its ticket, field or resolution: {title}"
        );

        // The three strata, live and cited, and the verdict that says the live
        // census reproduced the committed one rather than merely resembling it.
        let positive = line(&lines, "Delta > 0");
        assert!(
            positive.contains("4743"),
            "the Delta > 0 row does not carry p-130.csv's 4743: {positive}"
        );
        let negative = line(&lines, "Delta < 0");
        assert!(
            negative.contains("497"),
            "the Delta < 0 row does not carry p-130.csv's 497: {negative}"
        );
        let ambiguous = line(&lines, "ambiguous face");
        assert!(
            ambiguous.contains("27") && ambiguous.contains(":246"),
            "the ambiguous-face row lost its count or the line it names: {ambiguous}"
        );
        assert!(
            line(&lines, "EQUAL").contains("P-130"),
            "the verdict line no longer names whose lattice it matched"
        );

        let (headline, _) = banner.expect("a banner");
        assert!(
            headline.contains("Det(2,2,2)") && headline.contains("EQUAL to p-130.csv"),
            "the banner lost the identity or the verdict: {headline}"
        );

        // The normalisation, named in full so a reader can check the exponent.
        let magnitude = line(&lines, "abs(Delta)/max(abs(f_i))^4");
        assert!(
            magnitude.contains("min") && magnitude.contains("max"),
            "the magnitude row lost its range: {magnitude}"
        );
        assert!(
            line(&lines, "above 1e-4").contains("scale-free"),
            "the threshold row no longer says why the normalisation is exact"
        );
        let quoted = line(&lines, "P-134 measured");
        assert!(
            quoted.contains("3.445e-7") && quoted.contains("4.493e0"),
            "the panel stopped quoting p-134.csv's magnitude range beside the live one: \
             {quoted}"
        );
        let gap = line(&lines, "worst relative gap");
        assert!(
            gap.contains("within 1e-6") && gap.contains("two evaluation routes"),
            "the panel stopped saying how far the live magnitude sits from p-134.csv's, or \
             why it is not zero: {gap}"
        );

        // The f32 arm, and the two rows that stop it being read as exact.
        let f32_row = line(&lines, "f32 sign vs the f64 reference");
        assert!(
            f32_row.contains("case flips"),
            "the f32 row dropped its case-index coherence guard: {f32_row}"
        );
        assert!(
            line(&lines, "P-133").contains("samples at f32"),
            "the P-133 citation stopped saying its corner set is a different one"
        );
        assert!(
            line(&lines, "exact-rational 8-tuples").contains("14 of 3481"),
            "the P-127 citation stopped quoting its f32 disagreement count"
        );

        // The Delta = 0 branch, over P-131's population rather than P-130's.
        let branch = line(&lines, "over all");
        assert!(
            branch.contains("over all 32768 cells") && branch.contains(":250 fires"),
            "the branch row stopped naming the population it counts over: {branch}"
        );
        assert!(
            line(&lines, "P-131").contains("a == 0"),
            "the P-131 citation lost the linear-branch count beside the fired one"
        );
        assert!(
            line(&lines, "border rank 2").contains("surface cells"),
            "the panel stopped separating the branch's surface cells from its hits"
        );
        assert!(
            line(&lines, "tangential touch").contains("trilinear.rs:251"),
            "the panel stopped contradicting the comment P-131 falsified"
        );

        // The falsified clause, stated on screen and not only in the header.
        assert!(
            line(&lines, "rho(magnitude").contains("bar"),
            "the P-134 citation lost the bar its correlation is read against"
        );
        assert!(
            line(&lines, "C2 FALSIFIED").contains("3 of the 4 fields"),
            "the panel stopped saying how far short P-134's C2 fell"
        );

        // The heatmap's own honesty line.
        let soup = line(&lines, "unattributed");
        assert!(
            soup.contains("in a shared face") && soup.contains("cells emit none"),
            "the soup row lost the count of cells the picture cannot show, or the count of \
             triangles whose host cell was a choice between two: {soup}"
        );
    }

    /// `box_exact` is the falsified control, and the panel shows it as one.
    ///
    /// Every one of its surface cells is `Delta = 0`, so the surface is grey from
    /// edge to edge and P-134's correlation population excludes it outright. A
    /// demo that only ever showed a field with a spread of magnitudes would
    /// misrepresent how often the interesting case arrives, which is E-102's
    /// lesson about the same field picker.
    #[test]
    fn the_panel_shows_box_exacts_whole_surface_in_the_zero_stratum() {
        let (_, lines, banner) = panel(3, 0);
        for line in &lines {
            println!("{line}");
        }

        let zero = line(&lines, "not fixed by the sign");
        assert!(
            zero.contains("1352"),
            "box_exact's zero stratum is not p-130.csv's whole 1352 surface cells: {zero}"
        );
        assert!(
            line(&lines, "of surface cells are degenerate").contains("100.000%"),
            "the share row does not read 100% on a field where every cell is degenerate"
        );
        assert!(
            line(&lines, "in population:").contains("no, zero variance"),
            "the P-134 citation stopped recording box_exact's exclusion"
        );

        let (headline, _) = banner.expect("a banner");
        assert!(
            headline.contains("EQUAL to p-130.csv"),
            "the census no longer reproduces p-130.csv on box_exact: {headline}"
        );

        // `a == 0` on every cell of a polyhedron, and the `:250` branch never
        // reached. The two numbers are different facts and the panel keeps them
        // apart.
        let branch = line(&lines, "over all");
        assert!(
            branch.contains("28672"),
            "box_exact's a == 0 count is not p-131.csv's 28672: {branch}"
        );
    }

    /// Every one of the sixteen rows the demo can show reproduces `p-130.csv`'s
    /// partition and `p-131.csv`'s branch counts, cell for cell.
    ///
    /// **This is the gate the whole panel rests on.** The census here is not
    /// merely comparable to those CSVs, it is the same computation over the same
    /// numbers: the same lattice arithmetic, the same `cube::is_inside`, the same
    /// `BodySaddles::coefficients`, the same `4*a*c` association. So a
    /// difference of one cell on one row means this file, one of those CSVs, or
    /// the shipped coefficients have moved, and every quoted figure beside a live
    /// one has stopped being a comparison.
    ///
    /// It is cheap because the workspace's dev profile is optimised: sixteen
    /// rows, the largest visiting 262,144 cells, run in well under a second.
    #[test]
    fn every_row_reproduces_p130s_partition_and_p131s_branch_counts() {
        for rung in 0..LADDER.len() {
            for field in 0..FIELD_COUNT {
                let mut app = harness(field, rung);
                step(&mut app);
                let census = app.world().resource::<Census>();
                let cited = census.cited().expect("a committed row");
                assert_eq!(
                    census.field_name, cited.field,
                    "field {field} measured {} and the citation is {}",
                    census.field_name, cited.field
                );
                let gap = census.magnitude_gap(cited);
                println!(
                    "{:>14} {:>3}^3  {:>6} surface  {:>6}/{:>5}/{:>6}/{:>4}  a==0 {:>7}  \
                     :250 {:>3} ({} on the surface)  f32 dis {:>3}  no-tri {:>5}  \
                     coplanar {:>5}  mag {:.3e}..{:.3e}  gap {gap:.2e}",
                    census.field_name,
                    census.samples,
                    census.surface_cells,
                    census.delta_positive,
                    census.delta_negative,
                    census.delta_zero,
                    census.ambiguous_cells,
                    census.branch_a_zero,
                    census.branch_disc_zero,
                    census.branch_disc_zero_surface,
                    census.f32_disagreements,
                    census.cells_without_triangles,
                    census.face_coplanar,
                    census.magnitude_min,
                    census.magnitude_max,
                );
                assert!(
                    census.matches_p130(cited),
                    "{} at {}^3: the live partition is {}/{}/{}/{} and p-130.csv says \
                     {}/{}/{}/{}",
                    census.field_name,
                    census.samples,
                    census.delta_positive,
                    census.delta_negative,
                    census.delta_zero,
                    census.ambiguous_cells,
                    cited.partition[0],
                    cited.partition[1],
                    cited.partition[2],
                    cited.partition[3],
                );
                assert_eq!(
                    [census.branch_a_zero, census.branch_disc_zero],
                    cited.branch,
                    "{} at {}^3: the branch census disagrees with p-131.csv",
                    census.field_name,
                    census.samples
                );
                assert!(
                    gap <= MAGNITUDE_GAP_BOUND,
                    "{} at {}^3: the live normalised magnitude is {:e}/{:e}/{:e}, p-134.csv \
                     says {:e}/{:e}/{:e}, and the worst relative gap {gap:e} is outside \
                     {MAGNITUDE_GAP_BOUND:e} -- which is a disagreement about the polynomial \
                     rather than the two evaluation routes rounding differently",
                    census.field_name,
                    census.samples,
                    census.magnitude_min,
                    census.magnitude_mean,
                    census.magnitude_max,
                    cited.magnitude[0],
                    cited.magnitude[1],
                    cited.magnitude[2],
                );
                assert_eq!(
                    census.delta_positive + census.delta_negative + census.delta_zero,
                    census.surface_cells,
                    "{} at {}^3: the sign classes do not partition the surface cells",
                    census.field_name,
                    census.samples
                );
                assert_eq!(
                    census.non_finite_samples, 0,
                    "{} at {}^3 produced a non-finite sample",
                    census.field_name, census.samples
                );
                assert_eq!(
                    census.unattributed, 0,
                    "{} at {}^3 emitted triangles that lie in no single cell",
                    census.field_name, census.samples
                );
                assert_eq!(
                    census.f32_case_flips, 0,
                    "{} at {}^3 has cells whose case byte differs at f32",
                    census.field_name, census.samples
                );
                assert_eq!(
                    census.unnormalised, 0,
                    "{} at {}^3 has cells with no finite normalised magnitude",
                    census.field_name, census.samples
                );
                assert_eq!(
                    census.cages_dropped, 0,
                    "{} at {}^3 bound the {MAX_CAGES}-per-class cage cap, so the overlay is \
                     no longer showing the whole class",
                    census.field_name, census.samples
                );
            }
        }
    }

    /// The `f32` disagreements land on `csg_difference` and nowhere else, which
    /// is the honest shape of P-133's C1 rather than a demo-wide effect.
    ///
    /// P-133 measured `fields_disagreeing_f32 = 2` -- `csg_difference` and
    /// `fbm_terrain`, the latter only at 129 samples, which is off this demo's
    /// ladder. So a zero on `gyroid` here is the finding and not a hole in it,
    /// and the assertion is that the *population* is non-empty somewhere rather
    /// than that it is distributed any particular way.
    #[test]
    fn the_f32_sign_disagreements_are_rare_and_concentrated() {
        let mut clean = harness(0, 0);
        step(&mut clean);
        assert_eq!(
            clean.world().resource::<Census>().f32_disagreements,
            0,
            "gyroid at 33^3 gained an f32 sign disagreement; p-133.csv measures none there"
        );

        let mut sharp = harness(2, 0);
        step(&mut sharp);
        let census = sharp.world().resource::<Census>();
        assert!(
            census.f32_disagreements > 0,
            "csg_difference at 33^3 has no f32 sign disagreement, so the X overlay is empty on \
             every field and P-133's C1 has nothing on screen to show. p-133.csv measures 4 \
             there on a lattice sampled at f32."
        );
        assert!(
            census.f32_disagreements * 100 < census.surface_cells,
            "f32 sign disagreements are {} of {} surface cells on csg_difference, which is not \
             the rare stratum p-133.csv measured (0.288%)",
            census.f32_disagreements,
            census.surface_cells
        );
    }

    /// The `:250` branch fires on `csg_difference` and on no other reference
    /// field, and the live count is p-131.csv's.
    ///
    /// This is the population P-131's C1 is about, and it is over **every** cell
    /// rather than only the surface ones -- so a demo that censused the surface
    /// alone would report a different number under the same name.
    #[test]
    fn the_discriminant_zero_branch_fires_where_p131_measured_it() {
        let mut sharp = harness(2, 0);
        step(&mut sharp);
        let census = sharp.world().resource::<Census>();
        let cited = census.cited().expect("a committed row");
        assert_eq!(
            census.branch_disc_zero, cited.branch[1],
            "the :250 branch fired on {} cells of csg_difference at 33^3 and p-131.csv says {}",
            census.branch_disc_zero, cited.branch[1]
        );
        assert_eq!(
            census.branch_a_zero, cited.branch[0],
            "the a == 0 count is {} and p-131.csv says {}",
            census.branch_a_zero, cited.branch[0]
        );
        assert!(
            census.branch_disc_zero > 0,
            "p-131.csv's C1 is that this branch fires; a zero here would mean the demo cannot \
             show it"
        );

        let mut clean = harness(0, 0);
        step(&mut clean);
        assert_eq!(
            clean.world().resource::<Census>().branch_disc_zero,
            0,
            "the :250 branch fired on gyroid, where p-131.csv measures no hit"
        );
    }

    /// The `Delta = 0` stratum is off the ramp, and the ramp's anchors are
    /// P-134's decades.
    ///
    /// The stratum's colour has to be reachable **only** from the sign class:
    /// otherwise a magnitude that underflowed to zero on a cell with a perfectly
    /// good sign would be drawn as an exact zero, and the grey on screen would
    /// mean two different things. The anchors are checked because they are the
    /// thresholds `p-134.csv`'s `threshold_sweep` reports against rather than
    /// round numbers chosen to look good, so a reader can read a decade off the
    /// picture.
    #[test]
    fn the_zero_stratum_is_not_the_cold_end_of_the_ramp() {
        assert_ne!(
            ramp(0.0),
            cell_colour(0.0, 0),
            "the Delta = 0 stratum is drawn in the ramp's coldest colour, so an underflowed \
             magnitude and an exact zero are indistinguishable on screen"
        );
        assert_eq!(
            cell_colour(0.0, 1),
            ramp(0.0),
            "a signed cell whose magnitude underflowed no longer takes the ramp"
        );
        assert_eq!(
            cell_colour(1e-3, 0),
            cell_colour(0.0, 0),
            "the stratum's colour depends on the magnitude, so it is not a stratum"
        );

        // Clamped at both ends, exactly: the `clamp` lands on `0.0` and `1.0`
        // whatever the logarithm returned.
        assert_eq!(
            ramp(1e-30),
            ramp(1e-9),
            "the ramp does not clamp below 1e-8"
        );
        assert_eq!(ramp(1e6), ramp(1e1), "the ramp does not clamp above 1e0");

        // Each decade lands on its anchor. Not exact equality: `1e-6f64.log10()`
        // is not exactly `-6.0`, so the interpolation weight is near zero rather
        // than zero, and a test that demanded exactness would be testing libm.
        for (decade, anchor) in [1e-8, 1e-6, 1e-4, 1e-2, 1e0].iter().zip(RAMP.iter()) {
            let got = ramp(*decade);
            let want = linear(*anchor);
            for channel in 0..3usize {
                assert!(
                    (got[channel] - want[channel]).abs() < 1e-3,
                    "the ramp at {decade:e} is {got:?} and its anchor is {want:?}"
                );
            }
        }

        // Five distinct anchors, or the ramp cannot separate five decades.
        let mut seen: Vec<[f32; 4]> = RAMP.iter().map(|c| linear(*c)).collect();
        seen.push(linear(ZERO_STRATUM));
        seen.push(linear(REFUSED));
        for (i, a) in seen.iter().enumerate() {
            for b in seen.iter().skip(i + 1) {
                assert_ne!(a, b, "two of the palette's colours are the same");
            }
        }
    }

    /// The distinct linear colours on the surface asset, and its vertex count.
    fn painted(field: usize, rung: usize) -> (usize, usize, usize) {
        let mut app = harness(field, rung);
        step(&mut app);
        let mut query = app.world_mut().query_filtered::<&Mesh3d, With<DemoMesh>>();
        let handles: Vec<Handle<Mesh>> =
            query.iter(app.world()).map(|mesh| mesh.0.clone()).collect();
        assert_eq!(
            handles.len(),
            1,
            "field {field} at rung {rung} did not spawn exactly one surface"
        );
        let assets = app.world().resource::<Assets<Mesh>>();
        let mesh = assets
            .get(&handles[0])
            .expect("the surface asset the demo just added");
        let triangles = mesh.indices().map_or(0, |i| i.len() / 3);
        let colours = match mesh.attribute(Mesh::ATTRIBUTE_COLOR) {
            Some(VertexAttributeValues::Float32x4(values)) => values,
            Some(other) => panic!(
                "ATTRIBUTE_COLOR is not Float32x4 but {}",
                other.enum_variant_name()
            ),
            None => panic!(
                "the surface carries no ATTRIBUTE_COLOR, so the heatmap is not on screen at all \
                 -- MeshBuilder::into_mesh omits the attribute when colors_mut was never filled"
            ),
        };
        assert_eq!(
            colours.len(),
            mesh.count_vertices(),
            "one colour per vertex is `into_mesh`'s contract"
        );
        let mut distinct: Vec<[u32; 4]> = colours
            .iter()
            .map(|c| [c[0], c[1], c[2], c[3]].map(f32::to_bits))
            .collect();
        distinct.sort_unstable();
        distinct.dedup();
        (mesh.count_vertices(), triangles, distinct.len())
    }

    /// The surface really is a coloured soup, and its colours really vary with
    /// the magnitude.
    ///
    /// **This is the picture, asserted.** A `MeshBuilder` whose `colors_mut` was
    /// never filled omits `ATTRIBUTE_COLOR` entirely and renders a plain white
    /// surface -- so a heatmap that silently stopped painting would look like a
    /// demo about nothing and pass every count-based test above.
    ///
    /// The pair of fields is the whole argument. `gyroid`'s magnitudes span
    /// seven decades, so its surface must carry many colours; `box_exact`'s are
    /// all exactly zero, so its surface must carry **one** -- P-134 excluded it
    /// from the correlation population for precisely that reason, and the screen
    /// has to say so too.
    #[test]
    fn the_surface_is_a_coloured_soup_and_box_exact_is_one_flat_colour() {
        let (vertices, triangles, spread) = painted(0, 0);
        assert_eq!(
            vertices,
            triangles * 3,
            "the surface is not a soup: a shared vertex cannot carry one cell's colour"
        );
        assert!(
            spread > 100,
            "gyroid's surface carries only {spread} distinct colours over {triangles} \
             triangles, and its normalised magnitude spans seven decades"
        );

        let (flat_vertices, flat_triangles, flat_spread) = painted(3, 0);
        assert_eq!(flat_vertices, flat_triangles * 3);
        assert_eq!(
            flat_spread, 1,
            "box_exact's surface carries {flat_spread} distinct colours; every one of its \
             1352 surface cells has Delta exactly zero, so the whole surface is one flat \
             stratum colour and that is the falsified result the picture has to show"
        );
        assert_eq!(painted(3, 0).2, 1, "the flat reading is not reproducible");
    }
}
