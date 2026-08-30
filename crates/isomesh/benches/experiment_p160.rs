//! **P-160 — whether adaptivity beats uniform sampling, measured because the
//! famous theorem does not apply.**
//!
//! Ticket: R-160. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p160
//! ```
//!
//! Writes `docs/experiments/p-160.csv`.
//!
//! # What was missing
//!
//! Two things, and they pull in opposite directions.
//!
//! The first is the **theorem everyone reaches for**. Gal & Micchelli and Novak
//! (`10.1006/jcom.1996.0015`) prove that adaption improves worst-case error by at
//! most a factor of two. It is the single most quotable result in
//! information-based complexity and it is quoted at problems it does not cover,
//! because its four hypotheses are usually left in the paper: the class of
//! problem instances must be **convex**, it must be **symmetric**, the error must
//! be **worst-case**, and the solution operator must be a **continuous linear
//! mapping**. Nothing in this repository had ever written those four down beside
//! the crate's own setting, so nothing stopped a future reader from importing the
//! bound and concluding that the whole octree direction is capped at 2×. C3
//! exists to make that impossible, and it is answered with measurements rather
//! than with prose — see *C3, answered four times* below.
//!
//! The second is that **the crate has no adaptive sampler to measure**, and this
//! row could not add one: `crates/isomesh/src/**` is read-only for Phase 27.
//! What the crate does have is the *predicate* an octree is built out of.
//! `SubgridMarchingTetrahedra::cell_is_provably_empty` (subgrid/extract.rs:494-508)
//! is exactly
//!
//! ```text
//! sdf.sample(centre).abs() > l * (cell_size * 0.5 * 1.7320508075688772)
//! ```
//!
//! — a one-evaluation certificate that a cell holds no surface, with `l` from
//! [`FieldBound::lipschitz`](isomesh::fields::FieldBound::lipschitz) and a
//! **strict** comparison so that equality subdivides. It prunes cells of one
//! fixed 16³ subgrid; it is never recursed. `isomesh::lod` is a *downsampling*
//! module (lod.rs:13-15) and explicitly not a re-sampling one, so it refines
//! nothing. `P-48`'s own registration already records that this crate has no
//! balanced octree (experiment.rs:2386). So the octree here is bench-local, and
//! its refinement test is the crate's shipped predicate **negated and recursed
//! over levels** rather than a new criterion invented for the occasion.
//!
//! `P-159` (the group's null: Marching Cubes is already order-optimal and only
//! the constant is in play) and `P-161` (which approximation class `A^s` each
//! field is in) are the two rows either side of this one in Group G of
//! `docs/research/2026-08-29-phase-27-fifty-experiments-from-unmined-mathematics.md:390`.
//! This row is the empirical one: it does not ask what the rate is, it asks
//! whether spending a fixed number of field evaluations *adaptively* buys
//! anything, and by how much.
//!
//! # The refinement criterion, stated exactly, and the one the wording would
//! have given
//!
//! A cell at level `l` spans `edge = extent / 2^l`. It is **subdivided** iff
//!
//! ```text
//! signs_disagree(8 corners)  ||  |f(centre)| <= l_field * edge * 0.5 * sqrt(3)
//! ```
//!
//! and becomes a leaf otherwise. `l_field` is not a tuning knob: it is
//! `field.bound().lipschitz()`, which is `Some(1.0)` for every field this row
//! measures, and the crate's own instruction is to pass the constant rather than
//! guess it (subgrid/extract.rs:468-475, M-244).
//!
//! **Why this is sound, in full, because everything downstream rests on it.**
//! Every point of a cell lies within half the space diagonal of the cell centre.
//! For an `l_field`-Lipschitz field, `|f(centre)| > l_field * half_diagonal`
//! therefore gives `|f(p)| >= |f(centre)| - l_field * |p - centre| > 0` for every
//! `p` in the **closed** cell. A pruned cell contains no zero of the field, on
//! its boundary faces included.
//!
//! The second term subsumes the first: a sign disagreement among the corners
//! implies a zero inside the cell, which implies `|f(centre)| <= half_diagonal`.
//! `active_by_sign_only` counts the cells where the sign term fires and the
//! centre term does not, and the prediction is **zero**. It is a column rather
//! than a claim.
//!
//! **The registration's wording is not this test, and the difference is
//! measured rather than argued.** Read literally, "the field's range over the
//! corners exceeds the cell diagonal times a Lipschitz slack" is a *variation*
//! test. `corner_range_over_diagonal_max` is exactly that quantity — the worst
//! `(max corner − min corner) / (edge * sqrt(3))` over every cell this harness
//! ever tested — and for a 1-Lipschitz field it is **bounded above by 1 by
//! definition**, because two corners are at most one diagonal apart. So a
//! range-over-diagonal test with slack `>= 1` can never fire and one with slack
//! `< 1` fires on any cell whose diagonal happens to lie along the gradient,
//! which is most of them far from the surface: neither follows a surface. The
//! column reports the number; the criterion used is the crate's.
//!
//! `corner_range_over_diagonal_max <= l_field + 1e-9` is **asserted**, not
//! merely recorded. It is the pruning proof's only assumption, and a field that
//! violated its declared constant would hide surface in a coarse leaf and make
//! the adaptive arm's error something other than a measurement of adaptivity.
//!
//! # Why there are no cracks at the level transitions, and the column that
//! witnesses it
//!
//! Octree isosurface extraction is famous for cracks, and a crack in the
//! adaptive arm would be an error charged to adaptivity that is really an error
//! of stitching. There are none here, and the reason is the soundness argument
//! above rather than a repair pass: **a coarse leaf contains no zero of the
//! field**, so it emits no triangle, and a fine cell's face shared with a coarse
//! leaf lies inside that leaf's zero-free closure and therefore carries no
//! crossing either. Every triangle in the adaptive mesh is emitted by a cell at
//! the deepest level, whose six neighbours across a face are also at the deepest
//! level. The mesh is Marching Cubes over a locally uniform fine grid, restricted
//! to the shell where the surface is — which is the whole point.
//!
//! No dilation, no 2:1 balancing and no transition cells are needed, and none is
//! used. `adaptive_boundary_edges` is the witness — `boundary_edges` from
//! [`validate_indexed`](isomesh::validate::validate_indexed) — and it only
//! becomes one **after** the mesh is welded, which is the correction the first
//! run of this harness forced. `validate_indexed` keys its edges on vertex
//! *index*, not on position (validate.rs:736-744), and per-leaf Marching Cubes
//! emits every shared vertex once per leaf, so on the raw buffer every triangle
//! is an island and `boundary_edges` read 358,360 across these rows: a
//! measurement of unweldedness, saying nothing about cracks. Both arms now go
//! through the crate's own [`Welder`](isomesh::weld::Welder) at its own
//! [`epsilon_for`](isomesh::weld::epsilon_for) first — one instrument, applied
//! identically — and the coincident vertices merge because adjacent leaves
//! compute them from the same two corner values at bit-identical positions
//! (`extent` and every `cell_adaptive` here are exact powers of two times a
//! dyadic origin, so `origin + h*(i+1)` and `(origin + h*i) + h` are the same
//! `f64`). `adaptive_vertices_removed` and `uniform_vertices_removed` say what
//! each arm lost, `adaptive_triangles_collapsed` says whether the weld also
//! destroyed a sliver, and `uniform_boundary_edges` is the same crack number for
//! the other arm so neither is reported alone.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `uniform` | a single `n³` grid over the whole domain, `n` the smallest count with `n³ >= sample_budget` | **yes** — the non-adaptive information the theorem's `n` evaluations buy |
//! | `adaptive` | an octree refined by the crate's emptiness predicate to `max_level`, Marching Cubes over its deepest leaves | no |
//!
//! Both arms are Marching Cubes with the shipped defaults, both on `f64`, both
//! measured by the same [`accuracy`](isomesh::validate::accuracy) call against
//! the same reference lattice. The only thing that differs is **where the field
//! was evaluated**.
//!
//! # How the budget is matched, and in whose favour
//!
//! `sample_budget` is the number of **distinct field-evaluation positions** the
//! adaptive arm needed — the miss count of a cache keyed on the deepest lattice.
//! Every octree corner and every cell centre at every level lands on that
//! lattice (a level-`l` corner at a multiple of `res >> l`, a level-`l` centre at
//! an odd multiple of `res >> (l+1)`), so the cache is a flat array indexed by
//! integer coordinates: exact dedup, no hashing, no float comparison, and a miss
//! count that is the information the arm consumed. That is the right currency
//! for this question, because "restricted information" in the Novak sense is a
//! count of evaluations and nothing else.
//!
//! The uniform arm is then given the **smallest odd** `n` with `n³ >=
//! sample_budget`, so it receives *at least* the adaptive arm's budget and
//! usually a little more. `budget_ratio = n³ / sample_budget >= 1` is on every
//! row and reaches about `1.3` where the odd step is coarse. Rounding the other
//! way would have made any adaptive win partly a budget win; rounding this way
//! means a win is a win, and a `1.3` budget advantage is at most a `1.10`
//! cell-size advantage against gains this row reports in the single and double
//! digits. `evaluations_uncached` records what the harness actually spent
//! (`8 * leaves_finest` for the per-leaf Marching Cubes passes, which re-samples
//! rather than reading the cache) so the difference between the *information* and
//! the *implementation* is visible rather than implied.
//!
//! **Odd is load-bearing, and a first run proved it.** An octree of depth `L`
//! has `2^L + 1` lattice points per axis, so the centre of the domain is a
//! sample plane on every axis, always. A uniform grid of `n` samples has the
//! centre as a sample plane **iff `n` is odd** — for `n = 2m + 1` the spacing is
//! `extent / 2m` and index `m` lands exactly on the centre. With `n` free to be
//! even, the first run of this harness gave `thin_plate` at `max_level = 6` a
//! uniform arm of `n = 20`, whose vertical grid edges all straddle the plate
//! without either endpoint entering it — and classic Marching Cubes correctly
//! emitted **nothing at all**, 0 triangles against the adaptive arm's 4,088.
//! That is a *registered property of the field* and not a discovery
//! (fields/mod.rs:579-583, and M-266 at fields/mod.rs:606-617: "the plate is
//! centred at `y = 0` and every grid this crate measures on has an *odd* sample
//! count"). An infinite gain from it would have been a measurement of lattice
//! phase wearing adaptivity's clothes. The whole point of holding the budget
//! fixed is that **one** thing differs between the arms, so the phase is matched
//! and the house convention — `golden.rs:72`'s `[17, 25, 33]`, every
//! `experiment_pNN` resolution — is followed rather than accidentally broken.
//!
//! Three depths, `max_level` 5, 6 and 7 — finest cells of `extent/32`,
//! `extent/64` and `extent/128`. Not one depth, because the interesting quantity
//! is how the advantage *grows*: the octree pays for the surface area at its
//! finest spacing while a uniform grid pays for the volume, so the cell-size
//! ratio the adaptive arm buys at matched budget grows like `2^(L/3)` and one
//! point cannot show that. `cell_uniform` and `cell_adaptive` are columns.
//!
//! # How to read `gain`, derived before the numbers rather than after
//!
//! Two different laws are at work and confusing them would misread the CSV.
//!
//! On a **smooth** field the Hausdorff error of Marching Cubes is the linear
//! interpolation error on a curved surface, `O(h²)`, so at matched budget
//! `gain ≈ (cell_uniform / cell_adaptive)²`. `cell_ratio` is a column precisely
//! so the reader can square it and compare. `sphere` and `torus` are these.
//!
//! On a **polyhedral** field the surface is flat and the error is not in the
//! interpolation but in the corners, so it is `O(h)` and `gain ≈ cell_ratio`.
//! `box_exact` is this one, and it carries a **factor-of-two lattice-phase term
//! that the other three fields do not**, which is worth deriving here so that its
//! rows are not over-read. When the cube's face `x = 1` lands on a sample plane,
//! `f` is exactly `0` at those corners, `is_inside` (cube.rs:171, `v < 0`) calls
//! them outside, and the cell at the cube's vertex `(1,1,1)` has exactly one
//! corner inside. The three crossings all interpolate to `t = 1`, so the emitted
//! triangle is the plane `x + y + z = 3 − 2h` and the cube's own vertex sits
//! `(3 − (3 − 2h)) / sqrt(3) = 2h/sqrt(3) ≈ 1.1547 h` off the mesh. When the face
//! does not land on a sample plane the same corner truncation is half as deep,
//! `h/sqrt(3) ≈ 0.5774 h`.
//!
//! The adaptive arm is **always** in the first case: `x = 1` sits at lattice
//! index `3 * 2^(L−2)` of the octree's `2^L + 1` points, an integer for every
//! depth here, so `hausdorff_adaptive` on `box_exact` should read exactly
//! `2 * cell_adaptive / sqrt(3)` at all three depths. The uniform arm alternates
//! with `n`, so `box_exact`'s `gain` should read about `2 * cell_ratio` at some
//! depths and about `cell_ratio` at others — the same phase effect the odd-`n`
//! rule above removes from *whether* a surface is found, and cannot remove from
//! *where a polyhedron's corners fall*. It is a property of the fixture, it is
//! bounded by two, and `box_exact` is the only one of the four it touches.
//!
//! # The reference lattice, and the four fields that cannot have one
//!
//! Both arms of a row are measured by one `accuracy` call each against **one
//! shared lattice**, so the reverse (field→mesh) direction is sampled identically
//! for both. `reference_cells` is the largest power of two in `[16, 64]` whose
//! spacing is at least a sixth of the coarser arm's cell — the guard inside
//! `validate::accuracy`'s triangle grid rejects a triangle whose bounding box
//! spans more than 512 lattice cells (validate/tri_grid.rs:236), and a triangle
//! confined to a cell `k` lattice cells across spans at most `(k+1)³`, so `k <= 6`
//! keeps 343 under 512 with room. It is asserted, not hoped.
//!
//! `symmetric_hausdorff` is `max(mesh→field, field→mesh)` from the crate's own
//! report and is not recomputed here. In practice the forward direction
//! dominates, and that direction is **independent of the lattice**: it projects
//! every referenced vertex and every triangle centroid onto the zero set by
//! Newton and takes the largest displacement. Note the one bias this leaves, and
//! it runs against C1: the adaptive mesh has several times as many vertices as
//! the uniform mesh at matched budget, so it has several times as many chances to
//! find its own worst point. `mae_uniform` and `mae_adaptive` are beside the
//! maxima so a reader can see both.
//!
//! **`accuracy` is meaningless where `field.bound()` is not `Exact`** — the
//! value is not the distance, the band radius gates seeds on a first Newton step
//! that is not a distance, and a number would be a fiction. `csg_difference`
//! (`Underestimate { q: 0.5 }`), `gyroid` (`Lipschitz { l: 3.4641… }`),
//! `fbm_terrain` and `noise_cavity` (`Unbounded`) are therefore **skipped, and
//! the skip is a row**: `measured = false`, `bound` naming which of the four
//! kinds it is, and zeros in the distance columns. They are excluded from C1 and
//! C2 and from nothing else. Skipping them silently is the failure mode this
//! apparatus exists to prevent; four rows saying *why* cost nothing.
//!
//! The same gate is load-bearing twice over: the octree's pruning **is** the
//! Lipschitz certificate, so a field with no declared constant has no sound
//! octree either. The two exclusions have one cause.
//!
//! # C3, answered four times, by measurement
//!
//! The registration asks for "a row of booleans". Booleans that are only
//! asserted are prose in a numeric column, so each of the four is decided by a
//! fixture that runs in this bench and reports its own evidence. All four are
//! **global** — properties of the setting, not of a field — and are written
//! identically onto every row.
//!
//! - **`class_convex`.** Two unit spheres at `±(1.2, 0, 0)`. Both are members:
//!   exact distance fields with a closed non-empty surface. Their **midpoint**
//!   `(f + g)/2` is not, and not marginally: `min(|p+c| + |p−c|) = 2 * 1.2 = 2.4`
//!   along the segment between the centres, so the average is
//!   `(2.4 − 2)/2 = 0.2 > 0` there and larger everywhere else. The midpoint field
//!   has **no zero set at all**. `convex_midpoint_min_value` is that minimum
//!   over a 65³ grid whose origin is a lattice point (predicted `0.200000`), and
//!   `convex_midpoint_triangles` is the triangle count of a Marching Cubes pass
//!   over the same grid (predicted `0`), beside
//!   `convex_endpoint_triangles_min` for the two endpoints (predicted `> 0`, and
//!   asserted, or the midpoint's zero is M-44's vacuous zero).
//!   `class_convex = convex_midpoint_triangles > 0` → predicted **false**.
//!
//! - **`class_symmetric`.** Symmetric means `f ∈ F ⟹ −f ∈ F`. The canonical unit
//!   sphere and its negation are both extracted over `[-2, 2]³` at 65³. The
//!   surfaces coincide; the **winding** does not, because the crate's convention
//!   is negative-inside with normals away from the solid (lib.rs:56-67), so `−f`
//!   describes the unbounded complement. `symmetry_signed_volume` and
//!   `symmetry_signed_volume_negated` are the divergence-theorem volumes of the
//!   two meshes (predicted `+4.18879` and `−4.18879`, i.e. `±4π/3`): a member of
//!   a class of bounded solids cannot enclose a negative volume.
//!   `symmetry_hash_differs` ([`mesh_hash`](isomesh::validate::mesh_hash) of the
//!   two) proves the two extractions are genuinely different artefacts and not
//!   one computation reported twice.
//!   `class_symmetric = (v > 0) == (v_negated > 0)` → predicted **false**.
//!
//! - **`operator_linear`.** Decided by **positive homogeneity**, which is the
//!   cheapest half of linearity and fails exactly. Marching Cubes' crossing is
//!   `t = a/(a − b)`, invariant under scaling the field by a positive constant,
//!   and doubling is exact in IEEE-754 — so `S(2f)` is the mesh `S(f)`
//!   *bit for bit*, while a linear operator would have to produce `2 * S(f)`.
//!   `linear_scaled_positions_equal` compares the two position arrays by raw bit
//!   pattern (predicted **true**) and `linear_max_abs_position` is the largest
//!   absolute coordinate of `S(f)` (predicted `≈ 1`, asserted `> 0`, since
//!   `S(f) = 2 S(f)` would otherwise be satisfiable by a mesh at the origin).
//!   `operator_linear = !(scaled_positions_equal && max_abs_position > 0)` →
//!   predicted **false**.
//!
//!   The word *continuous* in the hypothesis is answered beside it and does not
//!   feed the verdict, because one mechanism per column: two unit spheres centred
//!   at `±(1, 0, 0)`, touching at the origin, offset by `ε = 0.08`. At `+ε` the
//!   level set is the boundary of two overlapping balls — one component. At `−ε`
//!   it is two disjoint balls — two components. The two fields are `2ε = 0.16`
//!   apart in the sup norm and the solution's component count is not the same
//!   number. `continuity_components_grown` (predicted `1`),
//!   `continuity_components_shrunk` (predicted `2`) and `continuity_epsilon`
//!   record it; `ε` is the smallest the fixture's 97³ grid resolves (the gap at
//!   `−ε` is `0.16`, about three cells), not a claim that the jump needs a
//!   perturbation that large.
//!
//! - **`error_worst_case`** (an extra column; the registration names three
//!   booleans and the fourth hypothesis belongs on the row too). This one **is
//!   satisfied**, and saying so is the point: symmetric Hausdorff is a sup-norm
//!   criterion and every number in the distance columns is a maximum, not a mean.
//!   `mae_uniform` and `mae_adaptive` sit beside them so the distinction is
//!   visible. Recorded **true**.
//!
//! `hypotheses_satisfied` counts the four (predicted `1`). The theorem needs all
//! four, so `c3_holds = hypotheses_satisfied < 4`: the four are checked, at least
//! one fails, and the bound cannot be imported.
//!
//! **The falsifier's strict reading is recorded separately rather than
//! reinterpreted.** *"C3 by any hypothesis being arguably satisfied"*, read
//! literally, is falsified by `error_worst_case = true`. `c3_strict_holds =
//! hypotheses_satisfied == 0` is therefore a column and is predicted **false**,
//! naming which hypothesis does it. Both readings are on the row; neither is
//! hidden, and a later reader can score C3 either way without re-running
//! anything.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"C1 moves the field-evaluation stage at fixed
//! budget"*, and that is a claim about a stage share, so it gets recomputed here
//! rather than repeated. What C1 moves is not the *cost* of the field-evaluation
//! stage — the budget is held fixed, so the stage costs the same by construction —
//! it is the **error at that cost**. There is no `1/(1 − share/factor)` ceiling
//! to compute, because the denominator is not a time. The reachable claim is
//! bounded instead by the cell-size ratio: `cell_uniform / cell_adaptive` is a
//! column, and no error ratio can exceed that ratio raised to the extractor's
//! order. `P-159` is the row that fixes the order, and this row records the ratio
//! it would be raised to rather than assuming one. A landing is a Phase 28
//! ticket against the shipped extractor, not a claim made here.
//!
//! # Which clauses are per-row and which are global
//!
//! `Run::record` writes rows, so a global verdict appears identically on all of
//! them. Stated once here so no reader has to infer it:
//!
//! - **`c1_holds` is global**: every measured row has `hausdorff_adaptive <
//!   hausdorff_uniform`. The per-row comparison is `c1_row`, and it is `false` on
//!   a skipped field because a field with no measurable distance is not a
//!   counterexample to anything.
//! - **`gain_exceeds_two` is per-row** and is the disjunct C2 is scored on;
//!   **`c2_holds` is global** and is its `OR` over the measured rows.
//! - **`c3_holds`, `c3_strict_holds`, `hypotheses_satisfied` and the three
//!   registered hypothesis booleans are global** — one setting, measured once.
//!
//! **Both clause comparisons are multiplications, never divisions**, and that is
//! deliberate. `hausdorff_adaptive` can be **exactly zero**: `box_exact`'s faces
//! land on lattice planes at every power-of-two resolution, so the octree arm can
//! reproduce the cube exactly, and P-100's harness already records a Hausdorff of
//! exactly 0 on `box_exact` at 33³. So C1 is `h_adaptive < h_uniform` and
//! `gain_exceeds_two` is `h_uniform > 2 * h_adaptive`, both total and exact at
//! zero. `gain` is reported as the ratio `h_uniform / h_adaptive` and will read
//! `inf` on such a row, with `gain_defined = false` beside it. `inf` in an
//! experiment CSV is house practice (p-114, p-53, p-54, p-73 all carry it) and no
//! clause reads the column.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `record`, and each names the column that
//! proves its fixture could have failed:
//!
//! - **The registered control**: `density_ratio > 1` on every measured row, per
//!   `(field, max_level)`. It is `8^(deepest leaf level − shallowest leaf level)`,
//!   the exact ratio of samples per unit volume between the finest and coarsest
//!   leaf. A ratio of 1 means the octree refined everything to one level and
//!   *both arms are uniform grids*, at which point every gain on the row is a
//!   measurement of the harness. `coarsest_leaf_level`, `leaves` and
//!   `leaves_finest` are beside it.
//! - **Both arms produced a surface**: `adaptive_triangles > 0` and
//!   `uniform_triangles > 0`.
//! - **Both arms were measured in both directions**: `has_coverage()` on both
//!   accuracy reports, with `reverse_samples_adaptive` and
//!   `reverse_samples_uniform` as the counts — a Hausdorff over an empty sample
//!   set is a zero that could not have been non-zero (M-44).
//! - **At least one field was measurable at all**: a run in which every field is
//!   skipped scores C1 and C2 over nothing.
//! - **The convexity fixture's endpoints have surfaces**:
//!   `convex_endpoint_triangles_min > 0`, or the midpoint's zero triangles prove
//!   nothing about convexity.
//! - **The homogeneity fixture's mesh is not at the origin**:
//!   `linear_max_abs_position > 0`, or `S(f) = 2 S(f)` is satisfiable and the
//!   contradiction evaporates.
//! - **The discontinuity fixture produced two meshes**:
//!   `continuity_components_grown > 0` and `continuity_components_shrunk > 0`,
//!   or the inequality between them is an inequality between two absences.
//!
//! Determinism: no RNG, no wall clock, no clause on a timing. The octree's cell
//! order is fixed by construction (level by level, children in corner-index
//! order from a parent list that starts as one root), the cache is keyed on
//! integers, and no comparison sorts floats.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]

mod common;

use isomesh::fields::{FieldBound, ReferenceField, Sphere};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_cubes::table::is_inside;
use isomesh::validate::{AccuracyConfig, ValidateConfig, accuracy, mesh_hash, validate_indexed};
use isomesh::weld::{Welder, epsilon_for};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

// ─── the fixture constants ───────────────────────────────────────────────────

/// The three octree depths, in levels below the domain cube.
///
/// Level `L` has a finest cell of `extent / 2^L`, so these are finest
/// resolutions of 32, 64 and 128 cells per axis. Three rather than one because
/// the quantity of interest is how the advantage **grows** with the budget: the
/// octree pays for surface area at its finest spacing, a uniform grid pays for
/// volume, so the cell-size ratio at matched budget grows like `2^(L/3)`.
const MAX_LEVELS: [u32; 3] = [5, 6, 7];

/// Half the space diagonal of a unit cube, as the crate spells it.
///
/// Copied byte for byte from `subgrid/extract.rs:504` so this harness's
/// emptiness certificate and the shipped one cannot drift apart in the last bit.
const HALF_DIAGONAL_UNIT: f64 = 1.732_050_807_568_877_2;

/// Cells per axis of the finest reference lattice both arms may be measured on.
const REFERENCE_CELLS_MAX: u32 = 64;

/// Cells per axis of the coarsest reference lattice this harness will accept.
///
/// Below this the reverse (field→mesh) direction has too few seeds to be worth
/// reporting, and a row that needed a coarser one would be a row whose arms are
/// coarser than the harness was designed for. It is asserted rather than
/// silently relaxed.
const REFERENCE_CELLS_MIN: u32 = 16;

/// Largest ratio of a mesh cell to a reference cell that
/// [`accuracy`](isomesh::validate::accuracy) will accept.
///
/// Its triangle grid refuses a triangle whose bounding box spans more than 512
/// reference cells (validate/tri_grid.rs:236). A triangle confined to a cell `k`
/// reference cells across spans at most `(k + 1)³`, so `k = 6` gives 343 and
/// leaves margin; `k = 7` would sit exactly on the limit.
const MAX_CELL_IN_REFERENCE_CELLS: f64 = 6.0;

/// How far above its declared constant a field's measured corner range may sit
/// before the octree's pruning is unsound.
///
/// Pure floating-point slack: the quantity is a difference of two samples over a
/// diagonal of at most 7 world units, so the rounding is nine orders below this.
const LIPSCHITZ_TOLERANCE: f64 = 1e-9;

/// Centre offset of the two spheres whose average refutes convexity.
///
/// `1.2 > 1`, so the two balls are disjoint and `|p + c| + |p − c| >= 2 * 1.2`
/// everywhere with equality on the segment between the centres. The average of
/// the two exact distance fields is therefore `(2.4 − 2) / 2 = 0.2` at its
/// minimum: strictly positive, so the midpoint of two members has no zero set.
const CONVEX_OFFSET: f64 = 1.2;

/// Half-extent of the convexity fixture's grid. `1.2 + 1 = 2.2` plus margin.
const CONVEX_HALF: f64 = 2.5;

/// Samples per axis for the convexity fixture.
///
/// Odd and centred over a symmetric box, so the origin — where the midpoint
/// field attains its minimum — is a lattice point and the reported minimum is
/// the true one rather than a nearby one.
const CONVEX_SAMPLES: u32 = 65;

/// Half-extent for the symmetry and homogeneity fixtures: the canonical unit
/// sphere's own domain.
const SPHERE_HALF: f64 = 2.0;

/// Samples per axis for the symmetry and homogeneity fixtures.
const SPHERE_SAMPLES: u32 = 65;

/// The factor the homogeneity fixture scales the field by.
///
/// A power of two, so `2a`, `2b` and `2a − 2b` are all exact and the crossing
/// parameter `a / (a − b)` is bit-identical to the unscaled one. Any other
/// factor would leave a rounding difference and turn an exact contradiction into
/// an approximate one.
const HOMOGENEITY_FACTOR: f64 = 2.0;

/// The perturbation that changes the component count of the touching-spheres
/// level set.
///
/// The two surfaces are `2 * epsilon = 0.16` apart in the sup norm, about three
/// cells of the fixture's grid — the smallest separation that grid resolves,
/// not a claim that the topology needs a perturbation this large.
const TOUCH_EPSILON: f64 = 0.08;

/// Half-extent of the touching-spheres grid. The grown balls reach `2.08`.
const TOUCH_HALF: f64 = 2.5;

/// Samples per axis for the touching-spheres fixture. Odd and centred, so the
/// touching point is a lattice point.
const TOUCH_SAMPLES: u32 = 97;

// ─── the fixtures that are not reference fields ──────────────────────────────

/// The pointwise average of two fields.
///
/// The midpoint of the segment between them in any linear space of fields, and
/// therefore a member of any **convex** class that contains both.
struct Averaged<A, B> {
    a: A,
    b: B,
}

impl<A, B> Sdf for Averaged<A, B>
where
    A: Sdf<Scalar = f64>,
    B: Sdf<Scalar = f64>,
{
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        f64::midpoint(self.a.sample(p), self.b.sample(p))
    }
}

/// `k · f`, for the positive homogeneity a linear operator would have to respect.
struct Scaled<F> {
    inner: F,
    k: f64,
}

impl<F: Sdf<Scalar = f64>> Sdf for Scaled<F> {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        self.k * self.inner.sample(p)
    }
}

/// `−f`, the member a **symmetric** class would have to contain.
struct Negated<F> {
    inner: F,
}

impl<F: Sdf<Scalar = f64>> Sdf for Negated<F> {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        -self.inner.sample(p)
    }
}

/// Two unit spheres centred at `±(1, 0, 0)`, so their surfaces touch at the
/// origin, with the whole field offset by `epsilon`.
///
/// The extracted level set is `{ min(d₁, d₂) = epsilon }`: two balls of radius
/// `1 + epsilon`. Positive `epsilon` overlaps them into one component, negative
/// separates them into two, and the two fields differ by `2 * epsilon` in the
/// sup norm.
struct TouchingSpheres {
    epsilon: f64,
}

impl Sdf for TouchingSpheres {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let ball = |cx: f64| {
            let q = [p[0] - cx, p[1], p[2]];
            (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt() - 1.0
        };
        ball(1.0).min(ball(-1.0)) - self.epsilon
    }
}

// ─── the sample cache: its miss count is the adaptive arm's budget ───────────

/// The octree's field-sample cache, keyed on the deepest level's lattice.
///
/// Every corner and every cell centre the refinement asks for lands on this
/// lattice, so a flat array indexed by integer coordinates dedups **exactly**:
/// no hashing, no float comparison, no tolerance. [`Cache::evaluations`] is the
/// number of distinct positions the field was evaluated at, which is the
/// information the adaptive arm consumed and therefore its sample budget.
struct Cache {
    /// Lattice points per axis, `(1 << max_level) + 1`.
    stride: usize,
    values: Vec<f64>,
    known: Vec<bool>,
    lo: [f64; 3],
    /// Spacing of the deepest lattice.
    h: f64,
    /// Distinct positions evaluated. The budget.
    evaluations: u64,
    /// Every query, cached or not. Reported so the caching claim is auditable.
    lookups: u64,
}

impl Cache {
    fn new(lo: [f64; 3], extent: f64, max_level: u32) -> Self {
        let res = 1usize << max_level;
        let stride = res + 1;
        let points = stride * stride * stride;
        Self {
            stride,
            values: vec![0.0; points],
            known: vec![false; points],
            lo,
            h: extent / res as f64,
            evaluations: 0,
            lookups: 0,
        }
    }

    fn position(&self, i: [u32; 3]) -> [f64; 3] {
        [
            self.lo[0] + self.h * f64::from(i[0]),
            self.lo[1] + self.h * f64::from(i[1]),
            self.lo[2] + self.h * f64::from(i[2]),
        ]
    }

    fn value<F: Sdf<Scalar = f64>>(&mut self, field: &F, i: [u32; 3]) -> f64 {
        let k = i[0] as usize + self.stride * (i[1] as usize + self.stride * i[2] as usize);
        self.lookups += 1;
        if !self.known[k] {
            self.values[k] = field.sample(self.position(i));
            self.known[k] = true;
            self.evaluations += 1;
        }
        self.values[k]
    }
}

// ─── the octree ──────────────────────────────────────────────────────────────

/// What the refinement produced.
struct Octree {
    /// Leaves at the deepest level, in that level's integer cell coordinates.
    ///
    /// The only cells that can carry surface, by the pruning argument in this
    /// file's header, and therefore the only ones Marching Cubes is run on.
    finest: Vec<[u32; 3]>,
    /// Every leaf, at every level.
    leaves: u64,
    /// Shallowest and deepest level any leaf sits at. `u32::MAX` and `0` when
    /// there are no leaves at all, which the triangle-count control catches.
    shallowest_leaf_level: u32,
    deepest_leaf_level: u32,
    /// Cells whose activity was tested, i.e. every cell at every level above the
    /// deepest.
    cells_tested: u64,
    /// Cells the corner-sign term made active while the centre term did not.
    ///
    /// Predicted **zero**: a sign disagreement implies a zero inside the cell,
    /// which implies the centre is within half a diagonal of the surface. A
    /// column rather than an argument.
    active_by_sign_only: u64,
    /// `(max corner − min corner) / (edge * sqrt(3))`, worst over every cell
    /// tested. Bounded by the field's declared Lipschitz constant, and that
    /// bound is the pruning proof's only assumption.
    corner_range_over_diagonal_max: f64,
}

impl Octree {
    /// `8^(deepest − shallowest)`: the exact ratio of samples per unit volume
    /// between the finest and the coarsest leaf. `0` when there are no leaves.
    fn density_ratio(&self) -> f64 {
        if self.shallowest_leaf_level > self.deepest_leaf_level {
            return 0.0;
        }
        let levels = self.deepest_leaf_level - self.shallowest_leaf_level;
        (1u64 << (3 * levels)) as f64
    }
}

/// Refine the domain cube until every cell either sits at `max_level` or is
/// certified surface-free.
///
/// The certificate is the crate's own `cell_is_provably_empty`
/// (subgrid/extract.rs:494-508) with `lipschitz` from the field's declared
/// bound, recursed over levels rather than applied once.
fn refine<F: Sdf<Scalar = f64>>(
    field: &F,
    cache: &mut Cache,
    max_level: u32,
    lipschitz: f64,
) -> Octree {
    let res = 1u32 << max_level;
    let mut current: Vec<[u32; 3]> = vec![[0, 0, 0]];
    let mut leaves = 0u64;
    let mut shallowest = u32::MAX;
    let mut deepest = 0u32;
    let mut cells_tested = 0u64;
    let mut sign_only = 0u64;
    let mut range_max = 0.0_f64;

    for level in 0..max_level {
        let step = res >> level;
        let half_step = step / 2;
        let edge = cache.h * f64::from(step);
        let diagonal = edge * HALF_DIAGONAL_UNIT * 2.0;
        let radius = edge * 0.5 * HALF_DIAGONAL_UNIT;
        let mut next: Vec<[u32; 3]> = Vec::new();

        for &c in &current {
            cells_tested += 1;
            let base = [c[0] * step, c[1] * step, c[2] * step];

            let mut lowest = f64::INFINITY;
            let mut highest = f64::NEG_INFINITY;
            let mut inside = 0u32;
            for k in 0..8u32 {
                let v = cache.value(
                    field,
                    [
                        base[0] + (k & 1) * step,
                        base[1] + ((k >> 1) & 1) * step,
                        base[2] + ((k >> 2) & 1) * step,
                    ],
                );
                lowest = lowest.min(v);
                highest = highest.max(v);
                if is_inside(v) {
                    inside += 1;
                }
            }
            range_max = range_max.max((highest - lowest) / diagonal);

            let signs_disagree = inside > 0 && inside < 8;
            let centre = cache.value(
                field,
                [
                    base[0] + half_step,
                    base[1] + half_step,
                    base[2] + half_step,
                ],
            );
            // The shipped predicate is `> l * radius` for *provably empty*, and
            // it is strict so equality subdivides. This is its negation.
            let near_surface = centre.abs() <= lipschitz * radius;
            if signs_disagree && !near_surface {
                sign_only += 1;
            }

            if signs_disagree || near_surface {
                for k in 0..8u32 {
                    next.push([
                        c[0] * 2 + (k & 1),
                        c[1] * 2 + ((k >> 1) & 1),
                        c[2] * 2 + ((k >> 2) & 1),
                    ]);
                }
            } else {
                leaves += 1;
                shallowest = shallowest.min(level);
                deepest = deepest.max(level);
            }
        }
        current = next;
    }

    if !current.is_empty() {
        leaves += current.len() as u64;
        shallowest = shallowest.min(max_level);
        deepest = deepest.max(max_level);
    }

    Octree {
        finest: current,
        leaves,
        shallowest_leaf_level: shallowest,
        deepest_leaf_level: deepest,
        cells_tested,
        active_by_sign_only: sign_only,
        corner_range_over_diagonal_max: range_max,
    }
}

// ─── the two arms ────────────────────────────────────────────────────────────

/// Marching Cubes over the octree's deepest leaves, accumulated into one buffer.
///
/// One `extract` per leaf on a `[2; 3]` grid — a single cell — into the same
/// sink, so `MeshSink::vertex`'s returned indices are already global and nothing
/// has to be re-offset. The extractor never clears its sink
/// (marching_cubes/mod.rs:237-252 resets only its own scratch), so the passes
/// accumulate. Vertices shared between adjacent leaves are emitted twice and are
/// bit-identical, which is what lets `validate_indexed`'s weld see the mesh as
/// one surface.
fn adaptive_mesh<F: Sdf<Scalar = f64>>(
    field: &F,
    lo: [f64; 3],
    h: f64,
    finest: &[[u32; 3]],
) -> MeshBuffer<f64> {
    let cell = RuntimeShape3::new([2; 3]).expect("a single cell is a legal grid");
    let mut mc = MarchingCubes::<f64>::new();
    let mut mesh = MeshBuffer::<f64>::new();
    for &c in finest {
        let origin = [
            lo[0] + h * f64::from(c[0]),
            lo[1] + h * f64::from(c[1]),
            lo[2] + h * f64::from(c[2]),
        ];
        mc.extract(field, &cell, origin, h, &mut mesh)
            .expect("marching cubes over one octree leaf");
    }
    mesh
}

/// Marching Cubes over one `samples³` grid spanning the field's whole domain.
fn uniform_mesh<F>(field: &F, samples: u32) -> MeshBuffer<f64>
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, h) = common::grid::<f64, _>(field, samples);
    let mut mc = MarchingCubes::<f64>::new();
    let mut mesh = MeshBuffer::<f64>::new();
    mc.extract(field, &shape, origin, h, &mut mesh)
        .expect("marching cubes over the uniform grid");
    mesh
}

/// The lattice both arms of one row are measured against.
struct Reference {
    shape: RuntimeShape3,
    origin: [f64; 3],
    config: AccuracyConfig,
    /// Cells per axis. On the row, so a reader can see how many seeds the
    /// reverse direction had.
    cells: u32,
}

/// One arm's geometry, measured against the shared reference lattice.
struct Arm {
    hausdorff: f64,
    mean_absolute_error: f64,
    coverage: bool,
    reverse_samples: u64,
    triangles: u64,
    vertices: u64,
    vertices_removed: u64,
    triangles_collapsed: u64,
    boundary_edges: u64,
    components: u64,
}

/// Weld the arm's mesh, then measure it.
///
/// **The weld is load-bearing and the first run of this harness proved it.**
/// `validate_indexed` keys its edges on vertex *index* and not on position
/// (validate.rs:736-744); `weld_epsilon` feeds only its `duplicate_vertices` and
/// `weld_buckets` counters. Per-leaf Marching Cubes emits every shared vertex
/// once per leaf, so on an unwelded adaptive mesh every triangle is an island and
/// `boundary_edges` counts `3 * triangles` — 358,360 of them across these rows on
/// the first run, which measures unweldedness and says nothing whatever about
/// cracks. So both arms go through the crate's own [`Welder`] at the crate's own
/// [`epsilon_for`] before anything is read off them: one instrument, applied
/// identically, and `vertices_removed` says how much each arm had to lose (very
/// nearly nothing on the uniform arm, which Marching Cubes already welds per grid
/// edge inside one `extract`). Only then is `boundary_edges` a crack witness.
fn measure<F: Sdf<Scalar = f64>>(
    field: &F,
    mesh: &mut MeshBuffer<f64>,
    cell_size: f64,
    reference: &Reference,
) -> Arm {
    let mut welder = Welder::<f64>::new();
    let weld = welder
        .weld(mesh, epsilon_for(cell_size))
        .expect("marching cubes produced a well-formed buffer");

    let acc = accuracy(
        &mesh.positions,
        &mesh.indices,
        field,
        &reference.shape,
        reference.origin,
        &reference.config,
    )
    .expect("the reference lattice describes this mesh");
    let vcfg = ValidateConfig::from_cell_size(cell_size).expect("positive cell size");
    let report = validate_indexed(&mesh.positions, &mesh.indices, &vcfg);
    Arm {
        hausdorff: acc.symmetric_hausdorff(),
        mean_absolute_error: acc.mean_absolute_error(),
        coverage: acc.has_coverage(),
        reverse_samples: acc.field_to_mesh.samples,
        triangles: acc.triangles,
        vertices: weld.vertices_after as u64,
        vertices_removed: weld.vertices_removed() as u64,
        triangles_collapsed: weld.triangles_collapsed as u64,
        boundary_edges: report.boundary_edges,
        components: report.components,
    }
}

// ─── C3: the four hypotheses, each answered by a fixture ─────────────────────

/// Marching Cubes over a centred cube grid, for the fixtures that are not
/// reference fields. Returns the mesh and the grid's cell size.
fn mesh_on_cube<F: Sdf<Scalar = f64>>(
    field: &F,
    half: f64,
    samples: u32,
) -> (MeshBuffer<f64>, f64) {
    let h = 2.0 * half / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("fixture grid fits u32");
    let mut mc = MarchingCubes::<f64>::new();
    let mut mesh = MeshBuffer::<f64>::new();
    mc.extract(field, &shape, [-half; 3], h, &mut mesh)
        .expect("marching cubes over the fixture grid");
    (mesh, h)
}

/// The smallest field value over a centred cube grid.
fn min_on_cube<F: Sdf<Scalar = f64>>(field: &F, half: f64, samples: u32) -> f64 {
    let h = 2.0 * half / f64::from(samples - 1);
    let mut lowest = f64::INFINITY;
    for z in 0..samples {
        for y in 0..samples {
            for x in 0..samples {
                let p = [
                    -half + h * f64::from(x),
                    -half + h * f64::from(y),
                    -half + h * f64::from(z),
                ];
                lowest = lowest.min(field.sample(p));
            }
        }
    }
    lowest
}

/// The volume a closed triangle mesh encloses, by the divergence theorem.
///
/// Positive for the crate's convention — negative inside, faces wound
/// counter-clockwise seen from outside (lib.rs:56-67) — and negative for a mesh
/// whose winding is inverted, which is exactly what negating the field produces.
fn signed_volume(mesh: &MeshBuffer<f64>) -> f64 {
    let mut six = 0.0_f64;
    for face in mesh.indices.as_chunks::<3>().0 {
        let a = mesh.positions[face[0] as usize];
        let b = mesh.positions[face[1] as usize];
        let c = mesh.positions[face[2] as usize];
        six += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
    six / 6.0
}

/// Edge-connected components of a mesh's face set, welded at the cell size.
fn components_of(mesh: &MeshBuffer<f64>, cell_size: f64) -> u64 {
    let cfg = ValidateConfig::from_cell_size(cell_size).expect("positive cell size");
    validate_indexed(&mesh.positions, &mesh.indices, &cfg).components
}

/// The four hypotheses of the factor-two theorem, and the evidence for each.
struct Hypotheses {
    class_convex: bool,
    class_symmetric: bool,
    operator_linear: bool,
    error_worst_case: bool,

    convex_midpoint_min_value: f64,
    convex_midpoint_triangles: u64,
    convex_endpoint_triangles_min: u64,

    symmetry_signed_volume: f64,
    symmetry_signed_volume_negated: f64,
    symmetry_hash_differs: bool,

    linear_scaled_positions_equal: bool,
    linear_max_abs_position: f64,

    continuity_epsilon: f64,
    continuity_components_grown: u64,
    continuity_components_shrunk: u64,
}

fn hypotheses() -> Hypotheses {
    // ── convexity: the midpoint of two members has no zero set ──────────────
    let left = Sphere::<f64> {
        center: [-CONVEX_OFFSET, 0.0, 0.0],
        radius: 1.0,
    };
    let right = Sphere::<f64> {
        center: [CONVEX_OFFSET, 0.0, 0.0],
        radius: 1.0,
    };
    let midpoint = Averaged { a: left, b: right };
    let convex_midpoint_min_value = min_on_cube(&midpoint, CONVEX_HALF, CONVEX_SAMPLES);
    let (mid_mesh, _) = mesh_on_cube(&midpoint, CONVEX_HALF, CONVEX_SAMPLES);
    let (left_mesh, _) = mesh_on_cube(&left, CONVEX_HALF, CONVEX_SAMPLES);
    let (right_mesh, _) = mesh_on_cube(&right, CONVEX_HALF, CONVEX_SAMPLES);
    let convex_midpoint_triangles = mid_mesh.triangle_count() as u64;
    let convex_endpoint_triangles_min =
        left_mesh.triangle_count().min(right_mesh.triangle_count()) as u64;
    let class_convex = convex_midpoint_triangles > 0;

    // ── symmetry: −f encloses a negative volume ─────────────────────────────
    let sphere = Sphere::<f64>::canonical();
    let (plus, _) = mesh_on_cube(&sphere, SPHERE_HALF, SPHERE_SAMPLES);
    let (minus, _) = mesh_on_cube(&Negated { inner: sphere }, SPHERE_HALF, SPHERE_SAMPLES);
    let symmetry_signed_volume = signed_volume(&plus);
    let symmetry_signed_volume_negated = signed_volume(&minus);
    let symmetry_hash_differs = mesh_hash(&plus) != mesh_hash(&minus);
    let class_symmetric = (symmetry_signed_volume > 0.0) == (symmetry_signed_volume_negated > 0.0);

    // ── linearity: S(2f) is S(f) bit for bit, and 2·S(f) is not ─────────────
    let (doubled, _) = mesh_on_cube(
        &Scaled {
            inner: sphere,
            k: HOMOGENEITY_FACTOR,
        },
        SPHERE_HALF,
        SPHERE_SAMPLES,
    );
    let linear_scaled_positions_equal = plus.positions.len() == doubled.positions.len()
        && plus.indices == doubled.indices
        && plus.positions.iter().zip(&doubled.positions).all(|(p, q)| {
            p.iter()
                .zip(q.iter())
                .all(|(x, y)| x.to_bits() == y.to_bits())
        });
    let linear_max_abs_position = plus
        .positions
        .iter()
        .flatten()
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    let operator_linear = !(linear_scaled_positions_equal && linear_max_abs_position > 0.0);

    // ── continuity: the component count jumps across epsilon ────────────────
    let (grown, grown_h) = mesh_on_cube(
        &TouchingSpheres {
            epsilon: TOUCH_EPSILON,
        },
        TOUCH_HALF,
        TOUCH_SAMPLES,
    );
    let (shrunk, shrunk_h) = mesh_on_cube(
        &TouchingSpheres {
            epsilon: -TOUCH_EPSILON,
        },
        TOUCH_HALF,
        TOUCH_SAMPLES,
    );
    let continuity_components_grown = components_of(&grown, grown_h);
    let continuity_components_shrunk = components_of(&shrunk, shrunk_h);

    Hypotheses {
        class_convex,
        class_symmetric,
        operator_linear,
        // The one hypothesis our setting satisfies: symmetric Hausdorff is a
        // sup-norm criterion and every distance column is a maximum.
        error_worst_case: true,
        convex_midpoint_min_value,
        convex_midpoint_triangles,
        convex_endpoint_triangles_min,
        symmetry_signed_volume,
        symmetry_signed_volume_negated,
        symmetry_hash_differs,
        linear_scaled_positions_equal,
        linear_max_abs_position,
        continuity_epsilon: TOUCH_EPSILON,
        continuity_components_grown,
        continuity_components_shrunk,
    }
}

// ─── one row ─────────────────────────────────────────────────────────────────

/// One `(field, max_level)` pair, or one skipped field.
#[derive(Default)]
struct Row {
    field: &'static str,
    bound: &'static str,
    measured: bool,
    max_level: u32,

    sample_budget: u64,
    uniform_samples: u64,
    uniform_samples_per_axis: u32,
    budget_ratio: f64,
    evaluations_uncached: u64,
    cache_lookups: u64,

    cell_uniform: f64,
    cell_adaptive: f64,
    cell_ratio: f64,
    reference_cells: u32,

    hausdorff_uniform: f64,
    hausdorff_adaptive: f64,
    mae_uniform: f64,
    mae_adaptive: f64,
    gain: f64,
    gain_defined: bool,
    gain_exceeds_two: bool,
    c1_row: bool,

    density_ratio: f64,
    leaves: u64,
    leaves_finest: u64,
    coarsest_leaf_level: u32,
    cells_tested: u64,
    active_by_sign_only: u64,
    corner_range_over_diagonal_max: f64,

    adaptive_triangles: u64,
    uniform_triangles: u64,
    adaptive_boundary_edges: u64,
    uniform_boundary_edges: u64,
    adaptive_components: u64,
    uniform_components: u64,
    adaptive_vertices: u64,
    uniform_vertices: u64,
    adaptive_vertices_removed: u64,
    uniform_vertices_removed: u64,
    adaptive_triangles_collapsed: u64,
    uniform_triangles_collapsed: u64,
    reverse_samples_adaptive: u64,
    reverse_samples_uniform: u64,
}

/// Which of the four kinds a field's declared bound is, as one CSV word.
fn bound_name(bound: FieldBound) -> &'static str {
    match bound {
        FieldBound::Exact => "exact",
        FieldBound::Lipschitz { .. } => "lipschitz",
        FieldBound::Underestimate { .. } => "underestimate",
        FieldBound::Unbounded => "unbounded",
    }
}

/// Every row for one reference field: three budgets if its bound is `Exact`, one
/// skip row otherwise.
fn collect<F>(name: &'static str, field: &F, rows: &mut Vec<Row>)
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();
    let extent = hi[0] - lo[0];
    let bound = bound_name(field.bound());

    // `accuracy` compares |f| against a true distance, and the octree's pruning
    // *is* a Lipschitz certificate. A field with no declared constant has
    // neither, and one cause excludes it from both.
    if !field.bound().is_exact() {
        println!("  {name:<15} bound {bound} — skipped, no Hausdorff is defined against it");
        rows.push(Row {
            field: name,
            bound,
            ..Row::default()
        });
        return;
    }
    let lipschitz = field
        .bound()
        .lipschitz()
        .expect("an exact bound declares a Lipschitz constant of one");

    assert!(
        (hi[1] - lo[1] - extent).abs() <= extent * f64::EPSILON * 8.0
            && (hi[2] - lo[2] - extent).abs() <= extent * f64::EPSILON * 8.0,
        "P-160: {name} has a non-cubic domain and the octree subdivides a cube"
    );

    for &max_level in &MAX_LEVELS {
        let mut cache = Cache::new(lo, extent, max_level);
        let tree = refine(field, &mut cache, max_level, lipschitz);
        let cell_adaptive = cache.h;

        assert!(
            tree.corner_range_over_diagonal_max <= lipschitz + LIPSCHITZ_TOLERANCE,
            "P-160: {name} at level {max_level} spreads {:.12} of a cell diagonal across its \
             corners, above the declared Lipschitz constant {lipschitz} the octree's pruning \
             rests on — a coarse leaf may hide surface and the adaptive arm's error would not \
             be a measurement of adaptivity",
            tree.corner_range_over_diagonal_max
        );

        let density_ratio = tree.density_ratio();
        assert!(
            density_ratio > 1.0,
            "VOID: {name} at level {max_level} put every leaf at one level, so the adaptive \
             arm is a uniform grid under another name and every gain on this row would be a \
             measurement of the harness (density_ratio {density_ratio:.6}, leaves {}, \
             finest {})",
            tree.leaves,
            tree.finest.len()
        );

        // The smallest **odd** cube of samples that is at least the adaptive
        // budget, so the uniform arm is never short-changed and its lattice has
        // the octree's phase. Odd is not a nicety: see the header. The octree's
        // deepest lattice is `2^L + 1` points, so the domain centre is always a
        // sample plane on every axis; an odd `n` puts the uniform arm's centre
        // plane in the same place, and an even one does not.
        let budget = cache.evaluations;
        let mut n = (budget as f64).cbrt().floor() as u32;
        if n < 3 {
            n = 3;
        }
        if n.is_multiple_of(2) {
            n += 1;
        }
        while u64::from(n).pow(3) < budget {
            n += 2;
        }
        let uniform_samples = u64::from(n).pow(3);
        let cell_uniform = extent / f64::from(n - 1);

        // The reference lattice: as fine as the accuracy harness's triangle-span
        // guard allows, and shared by both arms of this row.
        let coarsest = cell_uniform.max(cell_adaptive);
        let mut reference_cells = REFERENCE_CELLS_MAX;
        while reference_cells > REFERENCE_CELLS_MIN
            && f64::from(reference_cells) * coarsest > MAX_CELL_IN_REFERENCE_CELLS * extent
        {
            reference_cells /= 2;
        }
        assert!(
            f64::from(reference_cells) * coarsest <= MAX_CELL_IN_REFERENCE_CELLS * extent,
            "P-160: {name} at level {max_level} has a coarsest cell of {coarsest:.6} against a \
             domain of {extent:.6}, so even a {REFERENCE_CELLS_MIN}-cell reference lattice \
             would trip the accuracy harness's 512-cell triangle-span guard"
        );
        let reference_h = extent / f64::from(reference_cells);
        let reference = Reference {
            shape: RuntimeShape3::new([reference_cells + 1; 3])
                .expect("reference lattice fits u32"),
            origin: lo,
            config: AccuracyConfig::from_cell_size(reference_h)
                .expect("positive reference cell size"),
            cells: reference_cells,
        };

        let mut adaptive = adaptive_mesh(field, lo, cell_adaptive, &tree.finest);
        let mut uniform = uniform_mesh(field, n);

        let a = measure(field, &mut adaptive, cell_adaptive, &reference);
        let u = measure(field, &mut uniform, cell_uniform, &reference);

        assert!(
            a.triangles > 0 && u.triangles > 0,
            "VOID: {name} at level {max_level} produced {} adaptive and {} uniform usable \
             triangles, so at least one arm has no surface and its Hausdorff is a zero that \
             could not have been non-zero (M-44)",
            a.triangles,
            u.triangles
        );
        assert!(
            a.coverage && u.coverage,
            "VOID: {name} at level {max_level} was measured in one direction only (adaptive \
             reverse {}, uniform reverse {}), so the maximum is taken over a set the fixture \
             could not populate (M-44)",
            a.reverse_samples,
            u.reverse_samples
        );

        // Both clause comparisons are multiplications: `hausdorff_adaptive` can
        // be exactly zero on a lattice-aligned field, and a division there is a
        // reported `inf` rather than a verdict.
        let gain = u.hausdorff / a.hausdorff;
        let gain_defined = a.hausdorff > 0.0;
        let c1_row = a.hausdorff < u.hausdorff;
        let gain_exceeds_two = u.hausdorff > 2.0 * a.hausdorff;

        println!(
            "  {name:<15} L{max_level}  budget {budget:>8}  n {n:>3}  h_u {cell_uniform:.6} \
             h_a {cell_adaptive:.6}  H_u {:.9} H_a {:.9}  gain {gain:.4}  density {density_ratio:.0}",
            u.hausdorff, a.hausdorff
        );

        rows.push(Row {
            field: name,
            bound,
            measured: true,
            max_level,
            sample_budget: budget,
            uniform_samples,
            uniform_samples_per_axis: n,
            budget_ratio: uniform_samples as f64 / budget as f64,
            evaluations_uncached: 8 * tree.finest.len() as u64,
            cache_lookups: cache.lookups,
            cell_uniform,
            cell_adaptive,
            cell_ratio: cell_uniform / cell_adaptive,
            reference_cells: reference.cells,
            hausdorff_uniform: u.hausdorff,
            hausdorff_adaptive: a.hausdorff,
            mae_uniform: u.mean_absolute_error,
            mae_adaptive: a.mean_absolute_error,
            gain,
            gain_defined,
            gain_exceeds_two,
            c1_row,
            density_ratio,
            leaves: tree.leaves,
            leaves_finest: tree.finest.len() as u64,
            coarsest_leaf_level: tree.shallowest_leaf_level,
            cells_tested: tree.cells_tested,
            active_by_sign_only: tree.active_by_sign_only,
            corner_range_over_diagonal_max: tree.corner_range_over_diagonal_max,
            adaptive_triangles: a.triangles,
            uniform_triangles: u.triangles,
            adaptive_boundary_edges: a.boundary_edges,
            uniform_boundary_edges: u.boundary_edges,
            adaptive_components: a.components,
            uniform_components: u.components,
            adaptive_vertices: a.vertices,
            uniform_vertices: u.vertices,
            adaptive_vertices_removed: a.vertices_removed,
            uniform_vertices_removed: u.vertices_removed,
            adaptive_triangles_collapsed: a.triangles_collapsed,
            uniform_triangles_collapsed: u.triangles_collapsed,
            reverse_samples_adaptive: a.reverse_samples,
            reverse_samples_uniform: u.reverse_samples,
        });
    }
}

// ─── recording ───────────────────────────────────────────────────────────────

/// Six decimals: the house format for a dimensionless ratio.
fn r6(v: f64) -> String {
    format!("{v:.6}")
}

/// Nine decimals: distances, where the sixth is already inside the noise of a
/// 128-cell grid and a reader recomputing a ratio needs the digits.
fn d9(v: f64) -> String {
    format!("{v:.9}")
}

/// The verdicts and the C3 evidence, all global, carried onto every row because
/// `Run::record` writes rows and not a preamble.
struct Globals {
    c1: bool,
    c2: bool,
    c3: bool,
    c3_strict: bool,
    hypotheses_satisfied: u32,
    h: Hypotheses,
}

fn record(run: &mut common::experiment::Run, row: &Row, g: &Globals) {
    run.record(&[
        // ── registered, in registration order ──
        ("field", row.field.to_string()),
        ("sample_budget", row.sample_budget.to_string()),
        ("hausdorff_uniform", d9(row.hausdorff_uniform)),
        ("hausdorff_adaptive", d9(row.hausdorff_adaptive)),
        ("gain", r6(row.gain)),
        ("gain_exceeds_two", row.gain_exceeds_two.to_string()),
        ("class_convex", g.h.class_convex.to_string()),
        ("class_symmetric", g.h.class_symmetric.to_string()),
        ("operator_linear", g.h.operator_linear.to_string()),
        ("c1_holds", g.c1.to_string()),
        ("c2_holds", g.c2.to_string()),
        ("c3_holds", g.c3.to_string()),
        // ── extras (M-273) ──
        ("bound", row.bound.to_string()),
        ("measured", row.measured.to_string()),
        ("max_level", row.max_level.to_string()),
        ("uniform_samples", row.uniform_samples.to_string()),
        (
            "uniform_samples_per_axis",
            row.uniform_samples_per_axis.to_string(),
        ),
        ("budget_ratio", r6(row.budget_ratio)),
        ("evaluations_uncached", row.evaluations_uncached.to_string()),
        ("cache_lookups", row.cache_lookups.to_string()),
        ("cell_uniform", d9(row.cell_uniform)),
        ("cell_adaptive", d9(row.cell_adaptive)),
        ("cell_ratio", r6(row.cell_ratio)),
        ("reference_cells", row.reference_cells.to_string()),
        ("mae_uniform", d9(row.mae_uniform)),
        ("mae_adaptive", d9(row.mae_adaptive)),
        ("gain_defined", row.gain_defined.to_string()),
        ("c1_row", row.c1_row.to_string()),
        ("density_ratio", r6(row.density_ratio)),
        ("leaves", row.leaves.to_string()),
        ("leaves_finest", row.leaves_finest.to_string()),
        ("coarsest_leaf_level", row.coarsest_leaf_level.to_string()),
        ("cells_tested", row.cells_tested.to_string()),
        ("active_by_sign_only", row.active_by_sign_only.to_string()),
        (
            "corner_range_over_diagonal_max",
            d9(row.corner_range_over_diagonal_max),
        ),
        ("adaptive_triangles", row.adaptive_triangles.to_string()),
        ("uniform_triangles", row.uniform_triangles.to_string()),
        (
            "adaptive_boundary_edges",
            row.adaptive_boundary_edges.to_string(),
        ),
        (
            "uniform_boundary_edges",
            row.uniform_boundary_edges.to_string(),
        ),
        ("adaptive_components", row.adaptive_components.to_string()),
        ("uniform_components", row.uniform_components.to_string()),
        ("adaptive_vertices", row.adaptive_vertices.to_string()),
        ("uniform_vertices", row.uniform_vertices.to_string()),
        (
            "adaptive_vertices_removed",
            row.adaptive_vertices_removed.to_string(),
        ),
        (
            "uniform_vertices_removed",
            row.uniform_vertices_removed.to_string(),
        ),
        (
            "adaptive_triangles_collapsed",
            row.adaptive_triangles_collapsed.to_string(),
        ),
        (
            "uniform_triangles_collapsed",
            row.uniform_triangles_collapsed.to_string(),
        ),
        (
            "reverse_samples_adaptive",
            row.reverse_samples_adaptive.to_string(),
        ),
        (
            "reverse_samples_uniform",
            row.reverse_samples_uniform.to_string(),
        ),
        // ── extras (M-273): the C3 evidence, global ──
        ("error_worst_case", g.h.error_worst_case.to_string()),
        ("hypotheses_satisfied", g.hypotheses_satisfied.to_string()),
        ("c3_strict_holds", g.c3_strict.to_string()),
        (
            "convex_midpoint_min_value",
            d9(g.h.convex_midpoint_min_value),
        ),
        (
            "convex_midpoint_triangles",
            g.h.convex_midpoint_triangles.to_string(),
        ),
        (
            "convex_endpoint_triangles_min",
            g.h.convex_endpoint_triangles_min.to_string(),
        ),
        ("symmetry_signed_volume", d9(g.h.symmetry_signed_volume)),
        (
            "symmetry_signed_volume_negated",
            d9(g.h.symmetry_signed_volume_negated),
        ),
        (
            "symmetry_hash_differs",
            g.h.symmetry_hash_differs.to_string(),
        ),
        (
            "linear_scaled_positions_equal",
            g.h.linear_scaled_positions_equal.to_string(),
        ),
        ("linear_max_abs_position", d9(g.h.linear_max_abs_position)),
        ("continuity_epsilon", d9(g.h.continuity_epsilon)),
        (
            "continuity_components_grown",
            g.h.continuity_components_grown.to_string(),
        ),
        (
            "continuity_components_shrunk",
            g.h.continuity_components_shrunk.to_string(),
        ),
    ]);
}

// ─── main ────────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-160");

    common::experiment::run(prereg, |run| {
        // ── C3 first: the four hypotheses, measured once ────────────────────
        println!("C3 — the four hypotheses of the factor-two theorem, each measured:");
        let h = hypotheses();

        assert!(
            h.convex_endpoint_triangles_min > 0,
            "VOID: the convexity fixture's two endpoint spheres produced {} triangles between \
             them at their worst, so the midpoint field's {} triangles prove nothing about \
             convexity — the grid missed all three fields, not just the average",
            h.convex_endpoint_triangles_min,
            h.convex_midpoint_triangles
        );
        assert!(
            h.linear_max_abs_position > 0.0,
            "VOID: the homogeneity fixture's mesh sits at the origin, where S(f) = 2 S(f) is \
             satisfiable, so the bit-identical positions are not a contradiction"
        );
        assert!(
            h.continuity_components_grown > 0 && h.continuity_components_shrunk > 0,
            "VOID: the discontinuity fixture produced {} and {} components, so the inequality \
             between them would be an inequality between two absences",
            h.continuity_components_grown,
            h.continuity_components_shrunk
        );

        let hypotheses_satisfied = u32::from(h.class_convex)
            + u32::from(h.class_symmetric)
            + u32::from(h.operator_linear)
            + u32::from(h.error_worst_case);

        println!(
            "  class_convex      {}   midpoint min {:.9} over {} triangles \
             (endpoints at worst {})",
            h.class_convex,
            h.convex_midpoint_min_value,
            h.convex_midpoint_triangles,
            h.convex_endpoint_triangles_min
        );
        println!(
            "  class_symmetric   {}   enclosed volume {:+.6} against {:+.6} negated, \
             hash differs {}",
            h.class_symmetric,
            h.symmetry_signed_volume,
            h.symmetry_signed_volume_negated,
            h.symmetry_hash_differs
        );
        println!(
            "  operator_linear   {}   S(2f) == S(f) bit for bit: {}, max |position| {:.6}",
            h.operator_linear, h.linear_scaled_positions_equal, h.linear_max_abs_position
        );
        println!(
            "  error_worst_case  {}   and the topology jumps at eps {:.4}: {} component(s) \
             grown against {} shrunk",
            h.error_worst_case,
            h.continuity_epsilon,
            h.continuity_components_grown,
            h.continuity_components_shrunk
        );
        println!("  {hypotheses_satisfied} of 4 hypotheses satisfied\n");

        // ── the two arms, per field, per budget ─────────────────────────────
        println!("C1/C2 — uniform against octree-adaptive at matched sample budget:");
        let mut rows: Vec<Row> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            collect(name, &field, &mut rows);
        });

        let measured = rows.iter().filter(|r| r.measured).count();
        assert!(
            measured > 0,
            "VOID: every one of the eight reference fields was skipped, so C1 and C2 would be \
             scored over no measurement at all"
        );

        let c1 = rows.iter().filter(|r| r.measured).all(|r| r.c1_row);
        let c2 = rows
            .iter()
            .filter(|r| r.measured)
            .any(|r| r.gain_exceeds_two);
        let c3 = hypotheses_satisfied < 4;
        let c3_strict = hypotheses_satisfied == 0;

        let best = rows
            .iter()
            .filter(|r| r.measured && r.gain_defined)
            .map(|r| r.gain)
            .fold(0.0_f64, f64::max);
        println!(
            "\n  {measured} measured rows, {} skipped; best defined gain {best:.4}",
            rows.len() - measured
        );
        println!("  C1 {c1}   C2 {c2}   C3 {c3}   C3 (strict reading) {c3_strict}");

        // The witnesses the pruning argument stands or falls on, summarised so a
        // reader of the run log does not have to open the CSV to see them.
        let worst_range = rows
            .iter()
            .map(|r| r.corner_range_over_diagonal_max)
            .fold(0.0_f64, f64::max);
        let sign_only: u64 = rows.iter().map(|r| r.active_by_sign_only).sum();
        let adaptive_boundary: u64 = rows.iter().map(|r| r.adaptive_boundary_edges).sum();
        let uniform_boundary: u64 = rows.iter().map(|r| r.uniform_boundary_edges).sum();
        let collapsed: u64 = rows
            .iter()
            .map(|r| r.adaptive_triangles_collapsed + r.uniform_triangles_collapsed)
            .sum();
        println!(
            "  soundness: worst corner range over cell diagonal {worst_range:.9}; cells made \
             active by the sign term alone {sign_only}; boundary edges {adaptive_boundary} \
             adaptive against {uniform_boundary} uniform; {collapsed} triangles collapsed by \
             the weld"
        );

        let globals = Globals {
            c1,
            c2,
            c3,
            c3_strict,
            hypotheses_satisfied,
            h,
        };
        for row in &rows {
            record(run, row, &globals);
        }
    });
}
