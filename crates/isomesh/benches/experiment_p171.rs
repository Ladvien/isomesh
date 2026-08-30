//! **P-171 — how often the hypothesis every correctness theorem needs actually fails.**
//!
//! Ticket: R-171. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p171
//! ```
//!
//! Writes `docs/experiments/p-171.csv`.
//!
//! # What was missing
//!
//! Every PL-correctness theorem this crate has ever cited needs two hypotheses
//! before it says anything: the isovalue must be a **regular value** of `f`, and
//! the surface must meet each simplex **transversally**. Nothing in this
//! repository has ever counted how often either fails, so nothing has ever known
//! which cells are inside the theory and which are outside it.
//!
//! The nearest prior art is `isomesh::validate::cell_is_certified`
//! (`validate/isotopy.rs:126-140`) and the sweep over it, `isotopy_report`
//! (`:188`), which report `uncertified` active cells. That is a **different**
//! question and answering it does not answer this one. `cell_is_certified` is
//! Plantinga & Vegter's *interval* condition on the **trilinear interpolant**:
//! clause one is `0 not in box F(C)`, clause two is `<box grad F, box grad F> > 0`
//! over the cell, and its verdict is *"the trilinear surface inside this cell is
//! isotopic to a single facet"* (`isotopy.rs:117-124`). It is a sufficient
//! condition on `F`, evaluated by interval arithmetic over eight corner values,
//! and — read the first clause — **it passes every inactive cell trivially**
//! (`isotopy.rs:123-124`, `:129-132`). So it is structurally unable to see the
//! failure this row is mostly about: a cell the sign test calls *inactive* that
//! the surface nonetheless touches.
//!
//! The mechanism of that failure is a decision the crate made once, globally, and
//! wrote down at `cube.rs:157-173`: **a sample of exactly zero is outside**,
//! after Lengyel (2010 dissertation, §3.1.1), because strict `< 0` is what makes
//! `a - b` strictly non-zero on a cut edge and removes the interpolation's
//! division-by-zero guard entirely. That choice is correct and this row does not
//! argue with it. What it means, though, is that a corner sample landing exactly
//! on the isovalue is *silently absorbed*: `is_inside(0.0) == false`
//! (`cube.rs:171-173`), the corner is booked as outside, and if its neighbours
//! are also outside the cell reports case `0` and emits nothing while the surface
//! is demonstrably passing through one of its corners. That is a transversality
//! failure against the cell's **0-skeleton** and it is invisible to every
//! instrument the crate ships.
//!
//! It is also not hypothetical. `fields/mod.rs:606-617` records `M-266`, where a
//! claim that *"no grid phase can ever put a corner inside the plate"* was found
//! false by the canonical grid itself: the plate is centred at `y = 0`, every
//! grid this crate measures on has an odd sample count, so `y = 0` **is** a
//! sample plane. Lattice alignment is the mechanism, the crate already knew it in
//! one place, and nobody had swept for it.
//!
//! # What "non-transverse" means in this instrument, exactly
//!
//! A cell is **non-transverse** when either of two stated conditions holds. Both
//! are computed and both are recorded, because they are different failures with
//! different consequences and averaging them into one fraction without saying so
//! would make the fraction meaningless.
//!
//! **(a) Regular-value violation, `regular_value_violations`.** At least one of
//! the cell's eight corner samples is **exactly** the isovalue: `value == 0.0`.
//! The comparison is exact and deliberately so — this is not a tolerance test, it
//! is the test for the one input the sign rule cannot classify. `-0.0 == 0.0` in
//! IEEE 754, and `is_inside(-0.0)` is also `false`, so negative zero is counted
//! and is counted correctly. **This criterion is evaluated on every cell,
//! including inactive ones**, and it is the only criterion here that can flag a
//! cell the extractor never looked at. That is the whole point of it.
//!
//! **(b) Vanishing gradient, `gradient_norm_min`.** The field's own gradient
//! norm, `||grad f(p)||`, falls below [`GRADIENT_FLOOR`] at a point where the
//! surface passes through the cell — so `0` is not a regular value there and no
//! implicit-function-theorem constant survives. `grad f` is `Sdf::gradient`
//! (`sdf.rs:79-92`), which is the field's analytic gradient wherever one is
//! implemented and a six-sample central difference otherwise; it is the gradient
//! the crate itself uses, not a second one written here.
//!
//! The probe points are the **crossings on the cell's cut edges**, and that
//! choice is not arbitrary either: `cube::edge_offset` (`cube.rs:216-225`) is
//! *the* definition of where a crossing is — `P-61` made it so and
//! `cube.rs:214-220` says in as many words that every extractor routes through
//! it and that a second copy of the formula is what `x39` cost 216 golden hashes.
//! So the point set at which transversality is interrogated is byte-for-byte the
//! point set at which the crate places vertices.
//!
//! `non_transverse_cells` is `|(a) union (b)|` and `non_transverse_fraction` is
//! that over `(n - 1)^3`.
//!
//! ## The floor, and the promise not to ship an inert guard
//!
//! `x36 / M-351` and `P-124` are the standing warning here: `P-124`'s tolerance
//! was scaled by a quantity that was **exactly 0.0** on `box_exact` at every
//! resolution and never exceeded `7.63e-13` anywhere against a minimum deciding
//! value of `1.408e-3` — *nine orders of magnitude of inert guard*, and nobody
//! could see it from the CSV. So this row refuses to hide the same thing:
//! `gradient_floor` and `gradient_margin` (`= gradient_norm_min /
//! GRADIENT_FLOOR`) are **both columns**. If the margin is enormous on every row
//! then criterion (b) is inert on this fixture, that is a reported result rather
//! than a silence, and the reader can re-threshold from `gradient_norm_min`
//! without re-running anything.
//!
//! `gradient_norm_min` is written in scientific notation, `{:.6e}`. The quantity
//! spans orders of magnitude — an exact distance field reads `1.0` and a
//! near-critical point reads whatever it reads — and a fixed-point `0.000000`
//! would erase exactly the evidence the column exists to carry.
//!
//! ## The crossing coordinate, as an independent second view of (a)
//!
//! `edge_offset` returns `d = ((a + b)/2)/(a - b)`, a signed offset from the edge
//! **midpoint** in units of the edge, so `d` is in `[-1/2, +1/2]`. `|d| = 1/2`
//! means the crossing sits exactly on a grid **vertex**, which is the maximal
//! transversality failure against the 0-skeleton, and `d = 0` is a crossing at
//! the midpoint, maximally transverse. Algebraically `|d| = 1/2` requires an
//! endpoint to be exactly zero, and with `a == 0.0` the arithmetic is exact:
//! `((0 + b) * 0.5) / (0 - b)` is `-0.5` bit for bit for every finite non-zero
//! `b`. `max_abs_edge_offset` and `edges_at_half_offset` are therefore a second,
//! independent view of criterion (a) taken through the extractor's own
//! coordinate.
//!
//! They are **not** asserted equal to (a), and the reason is worth writing down:
//! in floating point `|d|` can reach exactly `1/2` without an exactly-zero
//! endpoint, because an endpoint negligible against its partner rounds
//! `(a + b) * 0.5` to `0.5 * (a - b)` — `a = -1e-20, b = 1.0` gives `d = -0.5`
//! exactly. That is *also* a transversality failure, and a stricter one than it
//! looks: the extractor pins the crossing to a grid vertex up to rounding. So the
//! two counts are reported side by side and neither is derived from the other.
//!
//! # The two overlap statistics, and why they are not the same statistic
//!
//! **C2, `overlap_with_ambiguous`, is the Jaccard index** `|N cap A| / |N cup A|`,
//! which is what the registration asks for by name ("set overlap reported"). The
//! ambiguous set comes from the crate and is not re-derived:
//! `marching_cubes::table::AMBIGUOUS_FACES` (`table.rs:196-231`) is a
//! `pub static [u8; 256]` built at compile time, and a cell is ambiguous when
//! `AMBIGUOUS_FACES[case] != 0`. The case index is the same `u8` the shipped
//! extractor computes, for the reason given under *Agreement* below.
//!
//! C2's verdict is `jaccard >= C2_JACCARD_FLOOR` = `0.10`, stated here because
//! "substantially" is not a number. Chance-level enrichment is also recorded, as
//! `ambiguous_lift` = `|N cap A| * cells / (|N| * |A|)`, so that a small Jaccard
//! can be read as *"both sets are small"* (lift about 1) or as *"these
//! populations avoid each other"* (lift < 1) rather than being one
//! undifferentiated zero. `c2_reachable` is `|N| > 0 && |A| > 0`: where it is
//! false, C2's zero is a zero that could not have been non-zero (`M-44`) and must
//! not be read as a falsification.
//!
//! ## Why C2 is unreachable on some fields, derived rather than discovered
//!
//! `AMBIGUOUS_FACES[case] != 0` requires a cube face whose four corner signs
//! **alternate** around the ring (`table.rs:198-201` and `build_ambiguous_faces`
//! at `:213-220`) — that is, the two inside corners are **diagonal**. For three
//! of the eight reference fields that is impossible, and the argument is exact
//! rather than empirical, so `ambiguous_cells == 0` there is a consequence of the
//! fixture's geometry and not a broken instrument.
//!
//! **A box.** `box_sample` (`fields/mod.rs:438-443`) is
//! `length(max(q, 0)) + min(max q, 0)`, which is negative exactly when every
//! `q_i < 0` — a product of open intervals. On an axis-aligned cube face two
//! coordinates range over `{x0, x1}` and `{y0, y1}`; if the diagonal pair
//! `(x0, y0)` and `(x1, y1)` are both inside then `x0`, `x1`, `y0` and `y1` each
//! pass their own interval test independently, so all four corners are inside. No
//! alternating face exists. That is `box_exact` and `thin_plate`, which share the
//! same `box_sample` (`fields/mod.rs:541`, `:651`).
//!
//! **A ball.** Inside is `|p - c| < R`. For a square with centre `O`, either
//! diagonal pair sums to the same thing from any point: `A - O = -(C - O)`, so
//! the cross terms cancel and `|c - A|^2 + |c - C|^2 = 2|c - O|^2 + 2m^2` where
//! `m` is the half-diagonal — identical for `B, D`. So `A, C` inside forces a sum
//! `< 2R^2` while `B, D` outside forces `>= 2R^2`, and the two sums are equal. No
//! alternating face exists. That is `sphere`.
//!
//! Nothing analogous holds for a difference, an intersection, a gyroid or a noise
//! volume, so those are the fields on which C2 can be asked at all. Which rows
//! those turned out to be is what `c2_reachable` records.
//!
//! **C3, `overlap_with_defects`, is a containment and not a Jaccard**:
//! `|N cap D| / |D|`, the share of defect cells that are non-transverse. `D` is
//! tiny beside `N` by construction, so a Jaccard over the union would be pinned
//! near zero by the size ratio alone and could not distinguish *"every defect is
//! non-transverse"* from *"no defect is"* — which is precisely the distinction
//! C3 asks for. The raw counts are also columns (`non_transverse_and_defect`,
//! `defect_cells`) so nothing is hidden behind the ratio, and `c3_reachable` =
//! `|D| > 0` marks the rows where a zero is unreachable rather than false.
//!
//! `D` is built from the shipped validators, again not re-derived:
//! `validate::validate_features` (`validate.rs:653`) for non-manifold **edges**
//! and **vertices** and for inconsistently oriented edges,
//! `validate::self_intersections` (`validate/self_intersection.rs:153`) for
//! intersecting triangle pairs, and **zero-area triangles**. `boundary_edges` is
//! deliberately **excluded**: `validate.rs:605-607` states it is not a defect on
//! an open field, and `fbm_terrain` is open (`closed_in_domain() == false`), so
//! counting it would make one field's whole surface a defect.
//!
//! ## Zero-area triangles, which is the defect criterion (a) actually predicts
//!
//! Criterion (a) is not an abstract complaint. A corner sample of exactly zero
//! forces `d = ±1/2` on **every** cut edge that ends at that corner, so all of
//! those crossings are placed on the same point — the grid vertex itself. Two of
//! them in one triangle and the triangle has zero area. So a degenerate triangle
//! is the mechanical consequence of a regular-value violation, and it is the one
//! defect C3 has a *reason* to expect. Leaving it out of `D` would have been
//! asking C3 about defect populations criterion (a) has no mechanism to cause.
//!
//! Neither validator returns *which* triangles are degenerate, only how many, so
//! the census locates them itself — and then proves it located the right ones.
//! The predicate is `self_intersection.rs:191-198`'s: `length(cross(b − a,
//! c − a))` is not strictly positive and finite. `vec3::{sub, cross, dot,
//! length}` are `pub(crate)` (`vec3.rs:17`, `:27`, `:32`, `:46`) so the
//! expressions are written out in the crate's own order of operations, and
//! [`census`] **asserts** the located count against
//! `SelfIntersectionReport::degenerate_triangles` on every row. That makes it a
//! mirror, not a second rule — `P-121`'s discipline again.
//!
//! The mirror also surfaced something nobody had recorded: **the crate's two
//! degeneracy predicates do not agree.** `MeshReport::degenerate_triangles`
//! (`validate_features`, an area test relative to `area_epsilon_rel = 1e-6`) and
//! `SelfIntersectionReport::degenerate_triangles` (this exact-zero-length test)
//! are both columns here — `near_degenerate_triangles` and
//! `degenerate_triangles` — because they are answers to different questions and
//! only the second one can be located. `D` uses the second.
//!
//! ## Mapping a defect back to a cell, and what that costs
//!
//! A defect is a vertex index or a triangle index; a cell is a base grid index.
//! The map is stated rather than assumed: `floor((p - origin) / h)` per axis,
//! clamped to `[0, n - 2]`. A marching-cubes vertex sits **on a grid edge**, so
//! it lies on the shared boundary of up to four cells and the floor rule
//! attributes it to exactly one of them — the lowest-indexed. That
//! under-attributes, and rather than pretend otherwise the census records the
//! **dilated** overlap too: `defect_cells_dilated` and
//! `overlap_with_defects_dilated` take the defect set out to the `2x2x2` block of
//! base indices `{i-1, i} x {j-1, j} x {k-1, k}`, which is a superset of every
//! cell whose closed box can contain the vertex and needs no floating-point
//! equality test to compute. If the strict and dilated numbers disagree, the
//! attribution rule is doing the work and the finding must say so.
//!
//! # Agreement: why this bench's grid is the shipped extractor's grid
//!
//! Two facts, and the second is measured rather than argued.
//!
//! The value grid is built with `origin[k] + cell_size * f64::from(index)`, which
//! is the expression `crate::sdf::sample_grid` uses at `sdf.rs:183-187` — the
//! function `MarchingCubes::extract` calls at `marching_cubes/mod.rs:240-247`.
//! Same operands, same order, so the arrays are bit-identical and the case index
//! at `mod.rs:259-268` is the case index computed here.
//!
//! And **`mesh_vertices` must equal `cut_grid_edges`** on every row, asserted.
//! Plain Marching Cubes places exactly one vertex per crossed grid edge —
//! `table.rs:128-131` states it and names `x1 / M-2 / M-22`'s `V_mc = C` — so a
//! disagreement means this bench and the extractor are not looking at the same
//! grid, and every number below would be about a fixture the crate never saw.
//! `P-121`'s rule: a new instrument's first job is to agree with the old one
//! where they overlap. Both counts are columns, not just an assert.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `<field>` at 17/33/65 samples | eight reference fields x three grids — the 24 census rows | no |
//! | `box_exact` at every resolution | the registration's own named control: the field whose faces lie exactly on the lattice | **yes** |
//! | `(samples - 1) % 4 == 0` | the resolution arithmetic that lets a corner land on a box face at all | **yes** |
//! | `mesh_vertices` vs `cut_grid_edges` | this bench's grid against the shipped extractor's | **yes** |
//! | `ambiguous_cells` per row | C2's second set; zero makes C2 unreachable on that row | **yes** |
//! | `defect_cells` per row | C3's second set; zero makes C3 unreachable on that row | **yes** |
//! | `gradient_margin` | how inert [`GRADIENT_FLOOR`] is, as a number | **yes** |
//!
//! The three resolutions are `17`, `33`, `65`. Each is `2^k + 1`, so
//! `(n - 1) % 4 == 0` and `box_exact`'s faces at `+/-1` land on samples at all
//! three — see the second vacuity control, which is where that arithmetic earns
//! its keep. They double, so the scaling of a fraction whose numerator counts
//! cells touching a **surface** and whose denominator counts a **volume** is
//! readable straight off three rows: it should fall like `1/n`.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE line is *"none — this partitions the cells into
//! in-theory and out-of-theory, which nothing currently does"*, and that is still
//! the right answer after writing the harness. Nothing here proposes a source
//! change, moves a nanosecond or touches a golden hash: `crates/isomesh/src/**`
//! is untouched, the mesh is the default `MarchingCubes` mesh, and the only write
//! is `docs/experiments/p-171.csv`.
//!
//! What the row does produce is a partition that did not exist, and one thing
//! that can be acted on later: `regular_value_violations` is a **grid** property,
//! not a field property. A half-cell shift of the origin removes every
//! exactly-equal corner sample on a lattice-aligned feature — which is precisely
//! the remedy `fields/mod.rs:614-616` already records for the thin plate, where
//! shifting by half a cell removes the surface entirely. Whether the extractor
//! should ever be told to do that is a Phase 28 question and is not asked here.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record`, and every message begins
//! `VOID: `.
//!
//! - **`box_exact` is in the census.** The registration names it as the control
//!   because its faces have exactly-equal samples; excluding it would exclude the
//!   phenomenon. Proved by the column `field` carrying it on three rows.
//! - **Every resolution can land a sample on a box face.** `box_exact` spans
//!   `[-1, 1]^3` inside a `[-2, 2]^3` domain, so the faces sit at sample index
//!   `(n - 1)/4` and `3(n - 1)/4`. If `(n - 1) % 4 != 0` no corner can be exactly
//!   on the surface and `regular_value_violations` would read zero **by the
//!   choice of resolution rather than by measurement**. `17`, `33` and `65` give
//!   `16`, `32`, `64`; the control asserts the arithmetic instead of trusting the
//!   constants.
//! - **`box_exact` reports a non-zero violation count at every resolution.**
//!   Proved by `regular_value_violations`. This is the assert that would fire if
//!   the sign rule, the sampler or the exact comparison were wrong.
//! - **Something in the fixture is ambiguous.** Proved by `ambiguous_cells`
//!   summed over the census. Without it C2's Jaccard is computed against an empty
//!   set everywhere and every zero is a zero that could not have been non-zero.
//! - **Every row has an active surface.** Proved by `active_cells`. A row with no
//!   surface reports zero for every clause and none of those zeros mean anything.
//! - **This bench's grid is the extractor's grid.** Proved by
//!   `mesh_vertices == cut_grid_edges`, and by `intersections.triangles ==
//!   mesh.triangle_count()`, which is what licenses reading
//!   `self_intersections`' triangle indices as indices into `mesh.indices`.
//! - **The located degenerate triangles are the crate's degenerate triangles.**
//!   Proved by `degenerate_triangles`, asserted equal to
//!   `SelfIntersectionReport::degenerate_triangles`. Without it `D` could be a
//!   population of this bench's own invention and C3 would be answering a
//!   question about nothing.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_cubes::table::{
    AMBIGUOUS_FACES, EDGE_AXIS, EDGE_CORNERS, EDGE_COUNT, edge_offset, is_inside,
};
use isomesh::validate::{ValidateConfig, self_intersections, validate_features};
use isomesh::{MeshBuffer, Sdf, Shape3, for_each_reference_field};

/// The three grids the census runs on, in **samples** per axis.
///
/// Each is `2^k + 1`, so `(n - 1) % 4 == 0` and `box_exact`'s faces land on
/// samples — the second vacuity control asserts that rather than assuming it.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// Below this, `||grad f||` at a crossing counts as vanished and the isovalue is
/// not a regular value there.
///
/// Six orders of magnitude below the eikonal value: every field in the fixture is
/// a signed distance or is built from signed distances, so `||grad f||` is `1`
/// wherever it is smooth. A norm of `1e-6` means the level set's local spacing is
/// a million times the value scale, and no theorem's constant survives that.
///
/// The floor and the distance from it are **both columns** (`gradient_floor`,
/// `gradient_margin`), so an inert guard is a number rather than a silence —
/// `P-124`'s nine orders of inert tolerance is the reason that rule exists.
const GRADIENT_FLOOR: f64 = 1e-6;

/// C2's bar: the Jaccard index at or above which the two populations count as
/// "substantially" overlapping.
///
/// The registration says "substantially" and "substantially" is not a number, so
/// one is fixed here, before the run, and `ambiguous_lift` is recorded beside it
/// so a reader can re-decide.
const C2_JACCARD_FLOOR: f64 = 0.10;

/// Flag bits over the cell array.
mod flag {
    /// The cell fails criterion (a) or criterion (b).
    pub(crate) const NON_TRANSVERSE: u8 = 1 << 0;
    /// `AMBIGUOUS_FACES[case] != 0`.
    pub(crate) const AMBIGUOUS: u8 = 1 << 1;
    /// A defect vertex floors into this cell.
    pub(crate) const DEFECT: u8 = 1 << 2;
    /// A defect vertex lies in this cell's closed box, up to the `2x2x2`
    /// dilation.
    pub(crate) const DEFECT_DILATED: u8 = 1 << 3;
}

/// One census row: one field at one resolution.
#[derive(Clone, Debug)]
struct Row {
    /// The reference field's name, as `for_each_reference_field!` gives it.
    field: &'static str,
    /// Samples per axis. Cells per axis is one less.
    resolution: u32,
    /// `(resolution - 1)^3`.
    cells: u64,
    /// Cells failing criterion (a) or (b).
    non_transverse_cells: u64,
    /// Cells with at least one corner sample exactly equal to the isovalue —
    /// criterion (a). Counted on inactive cells too.
    regular_value_violations: u64,
    /// Active cells whose minimum crossing gradient norm is below
    /// [`GRADIENT_FLOOR`] — criterion (b).
    non_transverse_gradient_cells: u64,
    /// Minimum `||grad f||` over every crossing on the grid.
    gradient_norm_min: f64,
    /// Grid samples exactly equal to the isovalue. The primitive count; each one
    /// touches up to eight cells.
    exact_zero_samples: u64,
    /// Grid edges with a sign change. Equals the mesh's vertex count.
    cut_grid_edges: u64,
    /// `max |d|` over every crossing, `d` from `edge_offset`. `0.5` is a crossing
    /// pinned to a grid vertex.
    max_abs_edge_offset: f64,
    /// Crossings with `|d| >= 0.5`.
    edges_at_half_offset: u64,
    /// Cells with a sign change.
    active_cells: u64,
    /// Cells with `AMBIGUOUS_FACES[case] != 0`.
    ambiguous_cells: u64,
    /// `|N cap A|`.
    non_transverse_and_ambiguous: u64,
    /// `|N cup A|`.
    non_transverse_or_ambiguous: u64,
    /// `|D|` under the floor rule.
    defect_cells: u64,
    /// `|D|` under the `2x2x2` dilation.
    defect_cells_dilated: u64,
    /// `|N cap D|` under the floor rule.
    non_transverse_and_defect: u64,
    /// `|N cap D|` under the dilation.
    non_transverse_and_defect_dilated: u64,
    /// Vertices in the default marching-cubes mesh.
    mesh_vertices: u64,
    /// Triangles in the default marching-cubes mesh.
    triangles: u64,
    /// `MeshReport::non_manifold_edges`.
    non_manifold_edges: u64,
    /// `MeshReport::non_manifold_vertices`.
    non_manifold_vertices: u64,
    /// `MeshReport::inconsistently_oriented_edges`.
    inconsistent_edges: u64,
    /// `SelfIntersectionReport::count()`.
    self_intersecting_pairs: u64,
    /// `SelfIntersectionReport::degenerate_triangles` — zero-area by the exact
    /// `length(cross(..)) > 0` test. These are the ones the census locates and
    /// puts in `D`.
    degenerate_triangles: u64,
    /// `MeshReport::degenerate_triangles` — the *area-epsilon* test, a different
    /// predicate that returns a different number and cannot be located.
    near_degenerate_triangles: u64,
    /// Wall clock for this row. Read by no clause; recorded so the phase can see
    /// where the census spends its time. `M-280`: this host's governor swings the
    /// same binary 1.45x, so it is not evidence about anything.
    wall_ms: f64,
}

impl Row {
    /// `non_transverse_cells / cells`.
    fn non_transverse_fraction(&self) -> f64 {
        self.non_transverse_cells as f64 / self.cells as f64
    }

    /// C2's statistic: `|N cap A| / |N cup A|`, zero when the union is empty.
    fn jaccard_with_ambiguous(&self) -> f64 {
        if self.non_transverse_or_ambiguous == 0 {
            0.0
        } else {
            self.non_transverse_and_ambiguous as f64 / self.non_transverse_or_ambiguous as f64
        }
    }

    /// Observed `|N cap A|` over what independence would give. Zero when either
    /// set is empty, in which case `c2_reachable` is false and this number says
    /// nothing.
    fn ambiguous_lift(&self) -> f64 {
        if self.non_transverse_cells == 0 || self.ambiguous_cells == 0 {
            return 0.0;
        }
        let expected =
            (self.non_transverse_cells as f64 * self.ambiguous_cells as f64) / self.cells as f64;
        self.non_transverse_and_ambiguous as f64 / expected
    }

    /// C3's statistic: the share of defect cells that are non-transverse.
    fn defect_containment(&self) -> f64 {
        if self.defect_cells == 0 {
            0.0
        } else {
            self.non_transverse_and_defect as f64 / self.defect_cells as f64
        }
    }

    /// The same containment under the `2x2x2` dilation of the defect set.
    fn defect_containment_dilated(&self) -> f64 {
        if self.defect_cells_dilated == 0 {
            0.0
        } else {
            self.non_transverse_and_defect_dilated as f64 / self.defect_cells_dilated as f64
        }
    }

    /// `|N| > 0 && |A| > 0`. False means C2's zero could not have been non-zero.
    fn c2_reachable(&self) -> bool {
        self.non_transverse_cells > 0 && self.ambiguous_cells > 0
    }

    /// `|D| > 0`. False means C3's zero could not have been non-zero.
    fn c3_reachable(&self) -> bool {
        self.defect_cells > 0
    }

    /// C2's verdict for this row.
    fn c2_holds(&self) -> bool {
        self.jaccard_with_ambiguous() >= C2_JACCARD_FLOOR
    }

    /// C3's verdict for this row: the non-transverse cells do overlap the defect
    /// cells.
    fn c3_holds(&self) -> bool {
        self.non_transverse_and_defect > 0
    }
}

/// Where the crossing goes on the segment `lo..hi`, given `d` from
/// [`edge_offset`].
///
/// Copied verbatim from the crate's `cube::place` (`cube.rs:233-235`), which is
/// `pub(crate)` and so unreachable from a bench. One line, and the centred frame
/// has to have exactly one spelling or `P-61`'s exact antisymmetry is lost in the
/// copy.
fn place(lo: f64, hi: f64, d: f64) -> f64 {
    (lo + hi) * 0.5 + (hi - lo) * d
}

/// The `(i & 1, (i >> 1) & 1, (i >> 2) & 1)` corner offset.
///
/// `cube::corner_offset` is private and is not in `table.rs:88-91`'s re-export
/// list, so the three lines are written here. The numbering they implement is
/// stated at `cube.rs:12-14`.
fn corner_offset(corner: u8) -> [u32; 3] {
    [
        u32::from(corner & 1),
        u32::from((corner >> 1) & 1),
        u32::from((corner >> 2) & 1),
    ]
}

/// Linear cell id from a base index, `x` fastest — the crate's index order.
fn cell_index(c: [u32; 3], cells_per_axis: u32) -> usize {
    let n = cells_per_axis as usize;
    c[0] as usize + c[1] as usize * n + c[2] as usize * n * n
}

/// The cell a mesh position floors into, clamped into the grid.
///
/// Stated rather than assumed: a marching-cubes vertex lies on a grid edge, so it
/// sits on the shared boundary of up to four cells and this names the
/// lowest-indexed one. The dilated set in [`census`] is what covers the rest.
fn cell_of(p: [f64; 3], origin: [f64; 3], h: f64, cells_per_axis: u32) -> [u32; 3] {
    let mut out = [0u32; 3];
    for k in 0..3 {
        let t = ((p[k] - origin[k]) / h).floor() as i64;
        out[k] = t.clamp(0, i64::from(cells_per_axis) - 1) as u32;
    }
    out
}

/// Is this triangle zero-area, by the predicate `self_intersections` uses?
///
/// `validate/self_intersection.rs:191-198` builds `cross(sub(b, a), sub(c, a))`,
/// takes its `length`, and excludes the triangle when that length is not strictly
/// positive and finite. `vec3::{sub, cross, dot, length}` are `pub(crate)`
/// (`vec3.rs:17`, `:27`, `:32`, `:46`), so the expressions are written out here
/// in the crate's own order of operations — `dot` is
/// `a0*b0 + a1*b1 + a2*b2` and `length` is its `sqrt`, associating left to right,
/// which is what makes the mirror bit-exact rather than merely close.
///
/// [`census`] asserts this predicate's count against
/// `SelfIntersectionReport::degenerate_triangles` on every row, so it is a mirror
/// of the shipped test and not a second opinion about degeneracy.
fn is_degenerate(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> bool {
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let length = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    !(length > 0.0 && length.is_finite())
}

/// The whole census for one field at one resolution.
fn census<F>(name: &'static str, field: &F, samples: u32) -> Row
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let started = Instant::now();

    let (shape, origin, h) = common::grid::<f64, _>(field, samples);
    let sample_count = (samples as usize).pow(3);
    let cells_per_axis = samples - 1;
    let cell_count = (cells_per_axis as usize).pow(3);

    // ── the value grid ──────────────────────────────────────────────────────
    //
    // `origin[k] + cell_size * f64::from(index)` is `sdf.rs:183-187` verbatim, so
    // this array is bit-identical to the one `MarchingCubes::extract` builds at
    // `mod.rs:240-247`. That is what makes the case indices below the shipped
    // extractor's case indices rather than a lookalike.
    let mut values = vec![0.0f64; sample_count];
    for z in 0..samples {
        for y in 0..samples {
            for x in 0..samples {
                let p = [
                    origin[0] + h * f64::from(x),
                    origin[1] + h * f64::from(y),
                    origin[2] + h * f64::from(z),
                ];
                values[shape.linearize([x, y, z]) as usize] = field.sample(p);
            }
        }
    }
    let exact_zero_samples = values.iter().filter(|v| **v == 0.0).count() as u64;

    // ── pass A: one gradient probe per cut grid edge ─────────────────────────
    //
    // Keyed the way `marching_cubes` keys its own edge cache (`mod.rs:250-251`
    // sizes it `sample_count * 3`): `3 * linearize(lower sample) + axis`. Probing
    // per *grid* edge rather than per *cell* edge is a 4x saving — a cut grid edge
    // is shared by four cells — and it makes `cut_grid_edges` directly comparable
    // with the mesh's vertex count.
    let mut edge_norm = vec![f64::INFINITY; sample_count * 3];
    let mut cut_grid_edges = 0u64;
    let mut gradient_norm_min = f64::INFINITY;
    let mut max_abs_edge_offset = 0.0f64;
    let mut edges_at_half_offset = 0u64;

    for axis in 0..3usize {
        let mut limit = [samples; 3];
        limit[axis] -= 1;
        for z in 0..limit[2] {
            for y in 0..limit[1] {
                for x in 0..limit[0] {
                    let lo_i = [x, y, z];
                    let mut hi_i = lo_i;
                    hi_i[axis] += 1;
                    let a = values[shape.linearize(lo_i) as usize];
                    let b = values[shape.linearize(hi_i) as usize];
                    if is_inside(a) == is_inside(b) {
                        continue;
                    }
                    cut_grid_edges += 1;

                    let d = edge_offset(a, b);
                    let mut p = [
                        origin[0] + h * f64::from(lo_i[0]),
                        origin[1] + h * f64::from(lo_i[1]),
                        origin[2] + h * f64::from(lo_i[2]),
                    ];
                    let hi_w = origin[axis] + h * f64::from(hi_i[axis]);
                    p[axis] = place(p[axis], hi_w, d);

                    let g = field.gradient(p);
                    let norm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                    edge_norm[3 * shape.linearize(lo_i) as usize + axis] = norm;
                    gradient_norm_min = gradient_norm_min.min(norm);

                    let abs_d = d.abs();
                    if abs_d > max_abs_edge_offset {
                        max_abs_edge_offset = abs_d;
                    }
                    if abs_d >= 0.5 {
                        edges_at_half_offset += 1;
                    }
                }
            }
        }
    }

    // ── pass B: the cell census ──────────────────────────────────────────────
    let mut flags = vec![0u8; cell_count];
    let mut regular_value_violations = 0u64;
    let mut non_transverse_gradient_cells = 0u64;
    let mut active_cells = 0u64;

    for k in 0..cells_per_axis {
        for j in 0..cells_per_axis {
            for i in 0..cells_per_axis {
                let mut corner = [0.0f64; 8];
                let mut case = 0u8;
                let mut exactly_zero = false;
                for c in 0..8u8 {
                    let o = corner_offset(c);
                    let v = values[shape.linearize([i + o[0], j + o[1], k + o[2]]) as usize];
                    corner[c as usize] = v;
                    if is_inside(v) {
                        case |= 1 << c;
                    }
                    if v == 0.0 {
                        exactly_zero = true;
                    }
                }

                let mut non_transverse = exactly_zero;
                if exactly_zero {
                    regular_value_violations += 1;
                }

                // Criterion (b) is a statement about where the surface passes, so
                // it is only asked of a cell the surface passes through. An active
                // cell always has at least one cut cube edge: two corners of
                // opposite sign are joined by an edge path, and some adjacent pair
                // on it must change sign.
                if case != 0 && case != 255 {
                    active_cells += 1;
                    let mut cell_min = f64::INFINITY;
                    for e in 0..EDGE_COUNT {
                        let lo_c = EDGE_CORNERS[e][0];
                        let hi_c = EDGE_CORNERS[e][1];
                        if is_inside(corner[lo_c as usize]) == is_inside(corner[hi_c as usize]) {
                            continue;
                        }
                        let o = corner_offset(lo_c);
                        let s = shape.linearize([i + o[0], j + o[1], k + o[2]]) as usize;
                        cell_min = cell_min.min(edge_norm[3 * s + EDGE_AXIS[e] as usize]);
                    }
                    if cell_min < GRADIENT_FLOOR {
                        non_transverse_gradient_cells += 1;
                        non_transverse = true;
                    }
                }

                let id = cell_index([i, j, k], cells_per_axis);
                if non_transverse {
                    flags[id] |= flag::NON_TRANSVERSE;
                }
                if AMBIGUOUS_FACES[case as usize] != 0 {
                    flags[id] |= flag::AMBIGUOUS;
                }
            }
        }
    }

    // ── the shipped mesh, and T-001's defect sets over it ────────────────────
    let mut mesh = MeshBuffer::<f64>::new();
    let mut mc = MarchingCubes::<f64>::new();
    mc.extract(field, &shape, origin, h, &mut mesh)
        .expect("the census grid has at least two samples on every axis");

    assert_eq!(
        mesh.vertex_count() as u64,
        cut_grid_edges,
        "VOID: {name} at {samples} samples: the default MarchingCubes mesh has {} \
         vertices and this bench counts {cut_grid_edges} cut grid edges. Plain \
         Marching Cubes places exactly one vertex per crossed grid edge \
         (table.rs:128-131, x1/M-2/M-22), so a disagreement means this census is \
         not looking at the extractor's grid and every fraction below is about a \
         fixture the crate never saw",
        mesh.vertex_count()
    );

    let cfg = ValidateConfig::from_cell_size(h)
        .expect("the census cell size is finite and positive by construction");
    let (report, features) = validate_features(&mesh.positions, &mesh.indices, &cfg);
    let intersections = self_intersections(&mesh.positions, &mesh.indices, h)
        .expect("a marching-cubes triangle never spans more than the broadphase guard's cells");

    assert_eq!(
        intersections.triangles,
        mesh.triangle_count() as u64,
        "VOID: {name} at {samples} samples: self_intersections filtered {} of {} \
         triangles, so its pair indices are indices into its own filtered list \
         and reading them as indices into mesh.indices would attribute defects to \
         the wrong cells",
        mesh.triangle_count() as u64 - intersections.triangles,
        mesh.triangle_count()
    );

    // Every vertex that participates in a defect. `boundary_edges` is excluded on
    // purpose (validate.rs:605-607: not a defect on an open field, and
    // `fbm_terrain` is open).
    let mut defect_vertices: Vec<u32> = Vec::new();
    defect_vertices.extend_from_slice(&features.vertices);
    for edge in &features.edges {
        defect_vertices.extend_from_slice(edge);
    }
    for edge in &features.inconsistently_oriented_edges {
        defect_vertices.extend_from_slice(edge);
    }
    for pair in &intersections.pairs {
        for &triangle in pair {
            let base = triangle as usize * 3;
            defect_vertices.extend_from_slice(&mesh.indices[base..base + 3]);
        }
    }

    // Zero-area triangles: the one defect criterion (a) has a mechanism to cause,
    // because an exactly-zero corner pins every crossing on its incident cut
    // edges to that grid vertex. Neither validator hands back which triangles
    // they are, so they are located here and the count is asserted against the
    // crate's.
    let mut degenerate_triangles = 0u64;
    for triangle in mesh.indices.as_chunks::<3>().0 {
        if is_degenerate(
            mesh.positions[triangle[0] as usize],
            mesh.positions[triangle[1] as usize],
            mesh.positions[triangle[2] as usize],
        ) {
            degenerate_triangles += 1;
            defect_vertices.extend_from_slice(triangle);
        }
    }
    assert_eq!(
        degenerate_triangles, intersections.degenerate_triangles,
        "VOID: {name} at {samples} samples: this bench locates \
         {degenerate_triangles} zero-area triangles and self_intersections counts \
         {}, so the located population is not the crate's and every cell C3 \
         attributes a degeneracy to is this bench's own invention \
         (self_intersection.rs:191-198 is the predicate being mirrored)",
        intersections.degenerate_triangles
    );

    for &v in &defect_vertices {
        let c = cell_of(mesh.positions[v as usize], origin, h, cells_per_axis);
        flags[cell_index(c, cells_per_axis)] |= flag::DEFECT;
        for dz in 0..2u32 {
            for dy in 0..2u32 {
                for dx in 0..2u32 {
                    let near = [
                        c[0].saturating_sub(dx),
                        c[1].saturating_sub(dy),
                        c[2].saturating_sub(dz),
                    ];
                    flags[cell_index(near, cells_per_axis)] |= flag::DEFECT_DILATED;
                }
            }
        }
    }

    // ── tally ────────────────────────────────────────────────────────────────
    let mut non_transverse_cells = 0u64;
    let mut ambiguous_cells = 0u64;
    let mut non_transverse_and_ambiguous = 0u64;
    let mut non_transverse_or_ambiguous = 0u64;
    let mut defect_cells = 0u64;
    let mut defect_cells_dilated = 0u64;
    let mut non_transverse_and_defect = 0u64;
    let mut non_transverse_and_defect_dilated = 0u64;
    for &f in &flags {
        let nt = f & flag::NON_TRANSVERSE != 0;
        let am = f & flag::AMBIGUOUS != 0;
        let de = f & flag::DEFECT != 0;
        let dd = f & flag::DEFECT_DILATED != 0;
        if nt {
            non_transverse_cells += 1;
        }
        if am {
            ambiguous_cells += 1;
        }
        if nt && am {
            non_transverse_and_ambiguous += 1;
        }
        if nt || am {
            non_transverse_or_ambiguous += 1;
        }
        if de {
            defect_cells += 1;
        }
        if dd {
            defect_cells_dilated += 1;
        }
        if nt && de {
            non_transverse_and_defect += 1;
        }
        if nt && dd {
            non_transverse_and_defect_dilated += 1;
        }
    }

    Row {
        field: name,
        resolution: samples,
        cells: cell_count as u64,
        non_transverse_cells,
        regular_value_violations,
        non_transverse_gradient_cells,
        gradient_norm_min,
        exact_zero_samples,
        cut_grid_edges,
        max_abs_edge_offset,
        edges_at_half_offset,
        active_cells,
        ambiguous_cells,
        non_transverse_and_ambiguous,
        non_transverse_or_ambiguous,
        defect_cells,
        defect_cells_dilated,
        non_transverse_and_defect,
        non_transverse_and_defect_dilated,
        mesh_vertices: mesh.vertex_count() as u64,
        triangles: mesh.triangle_count() as u64,
        non_manifold_edges: report.non_manifold_edges,
        non_manifold_vertices: report.non_manifold_vertices,
        inconsistent_edges: report.inconsistently_oriented_edges,
        self_intersecting_pairs: intersections.count(),
        degenerate_triangles,
        near_degenerate_triangles: report.degenerate_triangles,
        wall_ms: started.elapsed().as_secs_f64() * 1e3,
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-171");

    common::experiment::run(prereg, |run| {
        // ── measure ─────────────────────────────────────────────────────────
        //
        // Everything is measured before anything is asserted or recorded, because
        // C1 is a clause about the census as a whole — "non-zero on at least two
        // fields" cannot be decided from inside one field's block.
        // `for_each_reference_field!` inlines its body once per field with a
        // different concrete type each time, so there is no other way to hold all
        // eight; and no `return` appears in the body (M-253).
        let mut rows: Vec<Row> = Vec::new();
        for_each_reference_field!(f64, |name, field| {
            for samples in RESOLUTIONS {
                rows.push(census(name, &field, samples));
            }
        });

        // ── vacuity controls ────────────────────────────────────────────────
        assert!(
            rows.iter().any(|r| r.field == "box_exact"),
            "VOID: box_exact is not in the census, and the registration names it \
             as the control precisely because its faces have exactly-equal \
             samples. Without it the field most likely to violate transversality \
             has been excluded and every fraction reported here is a fraction \
             over a fixture that excluded the phenomenon"
        );

        for samples in RESOLUTIONS {
            let cells = samples - 1;
            assert!(
                cells % 4 == 0,
                "VOID: at {samples} samples the cell size over box_exact's \
                 [-2, 2] domain is 4/{cells}, so its faces at +/-1 sit at sample \
                 index {cells}/4 which is not an integer. No corner can land on \
                 the surface and regular_value_violations would read zero by the \
                 choice of resolution rather than by measurement"
            );
        }

        for row in rows.iter().filter(|r| r.field == "box_exact") {
            assert!(
                row.regular_value_violations > 0,
                "VOID: box_exact at {} samples reports no corner sample exactly \
                 equal to the isovalue, yet its faces at +/-1 land on sample \
                 index {} of {}. Either the sign rule, the sampler or the exact \
                 comparison is wrong, and the census cannot distinguish a clean \
                 fixture from a broken instrument",
                row.resolution,
                (row.resolution - 1) / 4,
                row.resolution - 1
            );
        }

        let ambiguous_total: u64 = rows.iter().map(|r| r.ambiguous_cells).sum();
        assert!(
            ambiguous_total > 0,
            "VOID: no cell anywhere in the census has an ambiguous face, so C2's \
             set overlap is computed against an empty set on all {} rows and \
             every zero is a zero that could not have been non-zero (M-44)",
            rows.len()
        );

        for row in &rows {
            assert!(
                row.active_cells > 0,
                "VOID: {} at {} samples has no cell with a sign change, so there \
                 is no surface, every clause reports zero and none of those zeros \
                 mean anything",
                row.field,
                row.resolution
            );
        }

        // ── C1, which is global ─────────────────────────────────────────────
        let mut fields_with_failures: Vec<&'static str> = rows
            .iter()
            .filter(|r| r.non_transverse_cells > 0)
            .map(|r| r.field)
            .collect();
        fields_with_failures.sort_unstable();
        fields_with_failures.dedup();
        let c1_fields_nonzero = fields_with_failures.len() as u64;
        let c1_holds = c1_fields_nonzero >= 2;

        // ── report ──────────────────────────────────────────────────────────
        println!(
            "gradient floor {GRADIENT_FLOOR:e}, C2 Jaccard bar {C2_JACCARD_FLOOR:.2}, \
             resolutions {RESOLUTIONS:?}"
        );
        println!(
            "C1: {c1_fields_nonzero} of 8 fields carry a non-transverse cell ({}) -> {}",
            fields_with_failures.join("|"),
            if c1_holds { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "{:>15} {:>5} {:>8} {:>7} {:>11} {:>6} {:>5} {:>10} {:>6} {:>8} {:>6} {:>6}",
            "field",
            "n",
            "cells",
            "non-tr",
            "fraction",
            "rvv",
            "grad",
            "gradmin",
            "ambig",
            "jaccard",
            "defec",
            "C2/C3"
        );

        for row in &rows {
            let jaccard = row.jaccard_with_ambiguous();
            let c2_holds = row.c2_holds();
            let c3_holds = row.c3_holds();
            println!(
                "{:>15} {:>5} {:>8} {:>7} {:>11.8} {:>6} {:>5} {:>10.2e} {:>6} {:>8.5} {:>6} {:>6}",
                row.field,
                row.resolution,
                row.cells,
                row.non_transverse_cells,
                row.non_transverse_fraction(),
                row.regular_value_violations,
                row.non_transverse_gradient_cells,
                row.gradient_norm_min,
                row.ambiguous_cells,
                jaccard,
                row.defect_cells,
                format!(
                    "{}{}/{}{}",
                    if c2_holds { "H" } else { "F" },
                    if row.c2_reachable() { "" } else { "*" },
                    if c3_holds { "H" } else { "F" },
                    if row.c3_reachable() { "" } else { "*" }
                )
            );

            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.resolution.to_string()),
                ("cells", row.cells.to_string()),
                ("non_transverse_cells", row.non_transverse_cells.to_string()),
                (
                    "non_transverse_fraction",
                    format!("{:.9}", row.non_transverse_fraction()),
                ),
                (
                    "gradient_norm_min",
                    format!("{:.6e}", row.gradient_norm_min),
                ),
                (
                    "regular_value_violations",
                    row.regular_value_violations.to_string(),
                ),
                ("overlap_with_ambiguous", format!("{jaccard:.9}")),
                (
                    "overlap_with_defects",
                    format!("{:.9}", row.defect_containment()),
                ),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_holds.to_string()),
                ("c3_holds", c3_holds.to_string()),
                // ── extras (M-273) ──
                ("active_cells", row.active_cells.to_string()),
                ("ambiguous_cells", row.ambiguous_cells.to_string()),
                ("ambiguous_lift", format!("{:.6}", row.ambiguous_lift())),
                ("c1_fields_nonzero", c1_fields_nonzero.to_string()),
                ("c2_reachable", row.c2_reachable().to_string()),
                ("c3_reachable", row.c3_reachable().to_string()),
                ("cut_grid_edges", row.cut_grid_edges.to_string()),
                ("defect_cells", row.defect_cells.to_string()),
                ("defect_cells_dilated", row.defect_cells_dilated.to_string()),
                ("degenerate_triangles", row.degenerate_triangles.to_string()),
                ("edges_at_half_offset", row.edges_at_half_offset.to_string()),
                ("exact_zero_samples", row.exact_zero_samples.to_string()),
                ("gradient_floor", format!("{GRADIENT_FLOOR:.6e}")),
                (
                    "gradient_margin",
                    format!("{:.6e}", row.gradient_norm_min / GRADIENT_FLOOR),
                ),
                ("inconsistent_edges", row.inconsistent_edges.to_string()),
                (
                    "max_abs_edge_offset",
                    format!("{:.9}", row.max_abs_edge_offset),
                ),
                ("mesh_vertices", row.mesh_vertices.to_string()),
                (
                    "near_degenerate_triangles",
                    row.near_degenerate_triangles.to_string(),
                ),
                ("non_manifold_edges", row.non_manifold_edges.to_string()),
                (
                    "non_manifold_vertices",
                    row.non_manifold_vertices.to_string(),
                ),
                (
                    "non_transverse_and_ambiguous",
                    row.non_transverse_and_ambiguous.to_string(),
                ),
                (
                    "non_transverse_and_defect",
                    row.non_transverse_and_defect.to_string(),
                ),
                (
                    "non_transverse_gradient_cells",
                    row.non_transverse_gradient_cells.to_string(),
                ),
                (
                    "overlap_with_defects_dilated",
                    format!("{:.9}", row.defect_containment_dilated()),
                ),
                (
                    "self_intersecting_pairs",
                    row.self_intersecting_pairs.to_string(),
                ),
                ("triangles", row.triangles.to_string()),
                ("wall_ms", format!("{:.3}", row.wall_ms)),
            ]);
        }

        // A `*` on a verdict above means that clause was unreachable on that row:
        // C2 needs both populations non-empty, C3 needs a defect to exist at all.
        // An unreachable clause is recorded as unreachable, with its arithmetic —
        // it is not a falsification and it is not a silence.
        let c2_reachable = rows.iter().filter(|r| r.c2_reachable()).count();
        let c3_reachable = rows.iter().filter(|r| r.c3_reachable()).count();
        println!(
            "\nC2 reachable on {c2_reachable} of {} rows, held on {}; C3 reachable on \
             {c3_reachable} of {} rows, held on {}",
            rows.len(),
            rows.iter().filter(|r| r.c2_holds()).count(),
            rows.len(),
            rows.iter().filter(|r| r.c3_holds()).count()
        );
        println!(
            "total wall {:.1} ms over {} rows",
            rows.iter().map(|r| r.wall_ms).sum::<f64>(),
            rows.len()
        );
    });
}
