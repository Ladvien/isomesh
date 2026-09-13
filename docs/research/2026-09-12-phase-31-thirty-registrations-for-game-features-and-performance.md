# Thirty registrations for novel game features and performance, from fourteen fields the ledger has not used

**Date:** 2026-09-12 · **Repo state:** `329747e` on `big` (`main` = `origin/main`), 606 findings entries, 76
open tickets
**Corpus:** home-still, 9,723 documents, 296,176 chunks, distill on CUDA
**Companion vocabulary:** `docs/research/2026-09-04-phase-29-axes-and-vocabulary-v3.md` is the latest on
`main` — v4 (Axes 19–23) has not landed — and the words for this phase's fields are in each group header.

**Numbering.** Phase 29 closed at `P-232` / `R-237`. Drafted as `P-277`–`P-306` against tickets
`R-282`–`R-311` on the `mac-air` checkout, numbered after a Phase 30 drafted there and ending at
`P-276`/`R-281`, then shifted **by −44** at the 2026-09-12 landing on `big`'s `main`, where **Phase 30 had
not landed** — landed numbers never move, the arriving text does, as `M-440`'s phase did twice and Phase 29
did once. Every artefact in this phase — this document, the `FINDINGS.md` section and the thirty backlog
rows — carries the final numbers **`P-233`–`P-262`** / **`R-238`–`R-267`**. Phase 30 will renumber above this
block when it arrives.

**Phase 30 is not on `main`, so its eleven cross-references are by title and carry no number.** The draft
cited `P-231`, `P-235`, `P-237`, `P-238`, `P-246`, `P-257`, `P-258`, `P-261`, `P-270`, `P-271` and `P-276`;
seven of those numbers now belong to *this* phase's rows, so a number here would be worse than no number.
Each reads as "Phase 30's ⟨instrument⟩ registration", and the record columns that named them read
`thickness_phase30_depth0`, `local_thickness_min_phase30` and `crease_edges_flagged_phase30_indicator`.
Back-linking them is Phase 30's landing's job, not this one's. Cross-references to Phase 29 rows were
resolved **by title, not by arithmetic**; the corrections paragraph below lists the three the draft got wrong
by content as well.

**What this phase is for.** Phases 27–30 were grades and guarantees — how to know a mesh is right. This phase
is the other half of the memory file: **features a game would notice** and **speed it would feel**, each from
a mathematical result the ledger has not touched, each with a hook into a number already measured, each
falsifiable on the crate's own roster. Tag per group: **[feature]**, **[perf]**, or both.

**Corpus status.** Probed at landing; §1 is the record. Every group header still names the field's canonical
result and the search terms, so a ticket can re-probe.

**Protocol.** Phase 15's rules. Every `P-` registered in `experiment.rs` before any harness commit; every
ratio-of-total carries SHARE; every row carries a VACUITY CONTROL; timings are minima over repetitions with
the count and CSV hash (Phase 30's timing-protocol registration: minima over repetitions with the count and
CSV hash). Falsified rows record which of the falsifier's meanings applies and end there.

**Eight registered nulls:** `P-235`, `P-238`, `P-243`, `P-246`, `P-249`, `P-254`, `P-256`, `P-261`.

**Corrections owed to the draft.** (1) The draft cited `✗42` for *crossing positions vary with evaluation
order within one machine*; `✗42` / `M-359` is `P-60`'s root-position gain and says nothing about order. The
within-machine order sensitivity is `M-32` (two expressions for one seam point, `1.57e-16` world units) and
`M-175` (a sum's result depends on tie order); Group F cites `M-32`. (2) The draft cited `P-183` for a
mesh-free Gaussian-curvature / mean-breadth instrument; on the landed numbering `P-183` is Phase 29's area
estimator, and the instrument described is `P-185` (`R-190`), integrated mean curvature from the sign grid.
`R-248` is blocked on it. (3) The draft's Group I proposes a Barnes–Hut dipole expansion after Barill et al.
2018. S-007's own backlog row already records the opposite verdict, from the 2026 Antipodal paper
(`10.1145/3811323`): the order-0/order-1 expansions are *"very imprecise… not useful for applications."*
Antipodal and Xie, Hafner & Wojtan (`10.1145/3811339`) are both in the corpus and both exact. `R-258`
measures those; Barnes–Hut is the baseline the ledger already distrusts, not the proposal.

---

## 1. What the fourteen probes found, in one table

| Area | Corpus | Repo cites | Verdict |
|---|---|---|---|
| Simplicial refinement (Freudenthal–Bey, Kuhn simplices, congruence classes) | ABSENT (0.62) → **acquired** Maubach; Bey 1995 and Bey 2000 **paywalled** | 11 (5 in `crates/`) | Group A |
| Multi-material meshing (multi-label level sets, non-manifold PL complexes) | ABSENT (0.60); Wu & Sullivan **paywalled** | 6 | Group B |
| Curvilinear grids on the sphere (cubed sphere, gnomonic projection) | ABSENT (0.62); Ronchi–Iacono–Paolucci **paywalled** | 0 | Group C |
| Singularity theory (apparent contours, folds and cusps) | ABSENT (0.59); Koenderink and Whitney **paywalled** | 0 | Group D |
| Mathematical morphology (offset surfaces, erosion and dilation, Lipschitz bounds) | ABSENT (0.58), no canonical paper | 4 (1 in `crates/`) | Group E |
| Floating-point determinism (reproducibility, FMA contraction, reduction order) | **PRESENT 0.66** Revol & Théveny | 30 (24 in `crates/`) | Group F — re-scoped: `M-32`, `M-175`, `M-31` |
| Discrete conformal geometry (conformal equivalence, cone singularities) | **PRESENT 0.64** Bobenko–Pinkall–Springborn; Springborn–Schröder–Pinkall 2008 **paywalled** | 5 | Group G |
| Heat-kernel PDE (the heat method, prefactored Laplacian) | **PRESENT 0.68** Crane–Weischedel–Wardetzky; Caissard et al. 2018 at 0.6766 | 16 (5 in `crates/`) | Group H |
| Generalized winding numbers (solid angle, inside–outside segmentation) | **PRESENT 0.72** Jacobson–Kavan–Sorkine-Hornung | 16 (7 in `crates/`) | Group I — re-scoped: S-006, S-007, `M-262`, `M-299`, `M-303` |
| Sparse voxel DAGs (hash-consing, subtree deduplication) | **PRESENT 0.77** Kämpe–Sintorn–Assarsson | 6 (1 in `crates/`) | Group J |
| Vertex quantisation (rate–distortion, uniform scalar quantiser) | ABSENT (0.62), no canonical paper | 70 (36 in `crates/`) | Group K |
| Space-filling curve analysis (Hilbert against Morton, clustering) | **PRESENT 0.64** Haverkort & van Walderveen; Moon et al. 2001 **paywalled** | 8 | Group L |
| Power diagrams (Laguerre tessellation, prescribed volumes) | ABSENT (0.58); Aurenhammer–Hoffmann–Aronov **paywalled** | 8 | Group M |
| Anisotropic diffusion (Perona–Malik, edge-stopping conductance) | ABSENT (0.53); Perona & Malik **paywalled** | 0 | Group N |

**The calibration, and the top hits that decided the absences.** Phase 29's rule was applied unchanged:
present at a top score of 0.65 or above with a relevant top hit, absent at 0.62 or below or with an unrelated
top hit, and the 0.62–0.65 band decided by the top hit's relevance with the score recorded either way. Six
areas came back present and eight absent. The absences are absences of the *canonical result*, not of the
field: Group A's top hit was Stevenson 2007's bisection-completion paper (`10.1090/s0025-5718-07-01959-x`) —
adjacent to Bey and not Bey; Group B's was *Sphere Carving*; Group C's was Springborn's ideal hyperbolic
polyhedra; Group D's was `10.1109_tvcg.2003.1207447`; Group E's was `10.1145_258734.258781`; Group K's was
*Bounded-Distortion Piecewise Mesh Parameterization*, with *High-Pass Quantization for Mesh Encoding* just
behind at 0.6144; Group M's was *Diagrammes de puissance restreint sur le GPU* — the right field and the
wrong result; and Group N's was a set of `eth-cgl-various-Pap13` slides. Two present verdicts fell inside the
band and were decided on relevance: Group G at 0.6382 (Bobenko–Pinkall–Springborn 2010,
`10.48550/arXiv.1005.2698`, whose canonical SSP 2008 companion is a different document) and Group L at 0.6416
(Haverkort & van Walderveen 2008, `10.48550/arXiv.0806.4787`, with Moon et al. 2001 absent). Repo cites were
counted over `crates/ docs/ bevy_isomesh/ isomesh_web/ FINDINGS.md BACKLOG.md BACKLOG_ARCHIVE.md`, `*.rs` and
`*.md`, case-insensitive, over files matching either of a group's two search terms — 478 files scanned.

**Two groups are re-scoped, because the probe came back present *and* the repo already measures part of the
hypothesis.** Group F: 24 files under `crates/` name `mul_add`, and `M-32`, `M-175` and `M-31`/`T-007`'s 378
golden hashes are the determinism entries already on the books; what is unmeasured is the thread count and a
cross-machine hash of the same edit. Group I: the crate ships both directions of the mesh ↔ field conversion
(`crates/isomesh/src/construct/winding.rs`, S-007; `crates/isomesh/src/construct/from_mesh.rs`, S-006) and
`M-262`, `M-299` and `M-303` priced them; what is unmeasured is S-006's own stated acceptance — mesh a
sphere, convert back to a field, re-mesh, compare — which has no Hausdorff or χ number in `FINDINGS.md`.
Groups A and K are heavily cited but their probes came back absent, so the re-scope rule does not fire; their
hooks already cite `M-452` and Phase 25 respectively.

**Acquired and in the corpus at time of writing (1).** Maubach 1995, *Local bisection refinement for
n-simplicial grids generated by reflection*, `10.1137/0916014` (downloaded and converted this session).

**Blocked acquisitions — paywalled or unreachable by the resolver this session.** Bey 1995 *Tetrahedral grid
refinement* `10.1007/bf02238487`; Bey 2000 *Simplicial grid refinement* `10.1007/s002110050475`; Wu &
Sullivan *Multiple material marching cubes algorithm* `10.1002/nme.775`; Ronchi–Iacono–Paolucci *The "Cubed
Sphere"* `10.1006/jcph.1996.0047`; Koenderink *What does the occluding contour tell us about solid shape?*
`10.1068/p130321`; Whitney *On singularities of mappings of Euclidean spaces. I* `10.2307/1970070`;
Springborn–Schröder–Pinkall *Conformal equivalence of triangle meshes* `10.1145/1360612.1360676`;
Moon–Jagadish–Faloutsos–Saltz *Analysis of the clustering properties of the Hilbert space-filling curve*
`10.1109/69.908985`; Aurenhammer–Hoffmann–Aronov *Minkowski-type theorems and least-squares clustering*
`10.1007/pl00009187`; Perona & Malik *Scale-space and edge detection using anisotropic diffusion*
`10.1109/34.56205`. Every title was matched against the canonical title in `paper_search` before the
identifier was used; no DOI or arXiv id was supplied from memory.

---

# Group A — LOD without a stitching table: Freudenthal–Bey refinement of the Kuhn triangulation (finite-element mesh refinement) [feature][perf]

**Field and words.** Bey 1995 *Tetrahedral grid refinement*; Bey 2000 *Simplicial grid refinement: on
Freudenthal's algorithm*; Maubach 1995 bisection on Kuhn simplices; Rivara longest-edge bisection. Terms:
`Freudenthal refinement`, `Kuhn simplex`, `congruence classes`, `conforming refinement`, `hanging node`,
`red-green refinement`. **The result [classical]:** a Kuhn tetrahedron refines into eight Kuhn-type
tetrahedra with a bounded number of congruence classes, so recursive refinement never degrades shape; and
bisection-based refinement of Kuhn simplices produces *conforming* meshes — no hanging nodes — with local
refinement closed by a bounded number of extra bisections. **The hook:** `M-474`'s LOD stitches levels with
a Transvoxel-style transition table; the Kuhn triangulation `M-452` already uses admits a conforming
refinement that needs no transition cases at all.

#### P-233 — registered for R-238, before the harness: Bey refinement of the crate's six-tet cell stays in a bounded congruence class set, so needle ratio does not grow with LOD depth

**Ticket:** `R-238` (S). **Records:** `refinement_depth`, `tet_count`, `congruence_classes_observed`, `congruence_classes_bound`, `thickness_min`, `thickness_phase30_depth0`, `needle_ratio_extracted`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) Refining each Kuhn tet by Freudenthal's rule to depth 4 produces at most the classical bound of congruence classes (three for `d = 3`), verified by canonicalising each tet's edge-length multiset. (C2) `thickness_min` at every depth equals the depth-0 value of Phase 30's tetrahedron-thickness registration (the `thickness_min` of the six-tet cell) to precision — quality is *exactly* preserved, not merely bounded — and the extracted needle ratio on `sphere` at mixed depth is within **10%** of uniform-depth extraction. SHARE: none.

**Falsified by.** C1 by more classes, which would mean the crate's six-tet split is the chiral variant Bey's rule does not preserve, and the row names which. C2 by degrading thickness, which would say the implementation bisected the wrong edge — Maubach's rule picks a specific one. VACUITY CONTROL: longest-edge bisection of a *non*-Kuhn regular tet must produce new congruence classes.

#### P-234 — registered for R-239, before the harness: Conforming Kuhn refinement replaces the LOD transition table, with no cracks and no special cases

**Ticket:** `R-239` (M). **Records:** `field`, `lod_layout`, `transition_faces`, `hanging_nodes`, `crack_edges` (boundary edges of the merged mesh), `transition_table_entries_used`, `chi_merged`, `chi_true_p145`, `hausdorff_vs_uniform`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** (C1) With coarse cells refined by bisection closure to meet fine neighbours, `hanging_nodes = 0` and `crack_edges = 0` on every field at every LOD layout the harness generates — conformity is a property of the triangulation, not of a case table. (C2) `transition_table_entries_used = 0`: the Transvoxel path is bypassed entirely and Marching Tetrahedra runs unchanged on the refined complex. (C3) `chi_merged = chi_true` wherever `M-474`'s stitched mesh also had it right, and `hausdorff_vs_uniform` is within **2×** of `M-474`'s. SHARE: none.

**Falsified by.** C1 by a crack, which would be a closure bug — the number of closure bisections is bounded and the row reports where the bound was exceeded. C3 by worse Hausdorff, which would say bisection puts vertices on cell edges where the crossing is poorly resolved — the price of conformity, and worth a number. VACUITY CONTROL: skipping the closure step must produce `crack_edges > 0`.

#### P-235 — registered for R-240, before the harness: A null registered on purpose — conforming refinement costs more triangles than Transvoxel at the transition, and the ratio is bounded by the closure depth

**Ticket:** `R-240` (S). **Records:** `field`, `triangles_conforming`, `triangles_transvoxel_m474`, `ratio`, `closure_bisections_per_transition_cell`, `c1_holds`.

**Hypothesis.** The registered null: (C1) `ratio ∈ [1.0, 1.5]` on every field — conformity is paid for in transition-cell triangles, and the price is bounded by the closure depth. SHARE: `ratio` is over the *whole* mesh; the transition cells' share of triangles is reported so a 1.5 on 5% of cells is read correctly.

**Falsified by.** C1 by `ratio < 1`, a strictly better result (Marching Tets on refined Kuhn cells emitting fewer triangles than a transition table is plausible on flat regions), or `> 1.5`, which prices the feature. VACUITY CONTROL: with no LOD (uniform depth) the ratio must be exactly 1.

---

# Group B — multi-material meshing (multi-label level sets, non-manifold PL complexes) [feature]

**Field and words.** Wu & Sullivan 2003 *Multiple material marching cubes*; multi-label / multi-phase level
sets; `non-manifold cell complex`, `material interface`, `triple junction`, `label configuration`. **The
result [classical]:** with `k` labels per corner, the cell has `k⁸` label configurations; the interface is a
non-manifold 2-complex whose 1-skeleton (triple-junction curves) is where the manifold case table fails, and
the number of *orbits* under `O_h × S_k` (cube symmetries times label permutations) is what a table must
enumerate. **The hook:** the roster is single-material; a game's terrain is not (rock / dirt / grass / ore),
and the `M-160` extraction has no seam path at all.

#### P-236 — registered for R-241, before the harness: The multi-label configuration count under `O_h × S_k` is small enough to table for `k ≤ 4`, and Burnside gives the number before enumeration

**Ticket:** `R-241` (S). **Records:** `k`, `configurations_total`, `orbits_burnside`, `orbits_enumerated`, `agree`, `orbits_with_triple_junction`, `table_bytes`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) Burnside's lemma over `O_h × S_k` (48·k! elements) predicts the orbit count for `k = 2, 3, 4`, and the exhaustive enumeration (the `A-002f` machinery extended to labels) agrees — the `k = 2` case must reproduce the 22 orbits of MC with complement. (C2) For `k = 4` the orbit count is under **5,000** and `table_bytes` under 1 MB — tableable, not a runtime decision. SHARE: none.

**Falsified by.** C1 by disagreement, a bug in the group action (label permutation composed with cube symmetry has a non-obvious order). C2 by a larger count, which prices multi-material at a per-cell solve rather than a lookup, and the row records the crossover `k`. VACUITY CONTROL: `k = 1` must give exactly one orbit.

#### P-237 — registered for R-242, before the harness: Material seams are closed curves on the mesh, and the mesh-free `χ` of each material's solid grades them

**Ticket:** `R-242` (M). **Records:** `fixture`, `k`, `resolution`, `material_id`, `chi_material_mesh`, `chi_material_oracle_p145`, `seam_curves`, `seam_curves_closed`, `triple_junction_vertices`, `nonmanifold_edges_expected`, `nonmanifold_edges_found`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) On a three-material fixture (three overlapping spheres, each cell's label by nearest centre, a triple junction where all three meet), each material's own solid — `label = i` as a binary field — has `chi_material_mesh = chi_material_oracle` at 33³ and above, using `P-145`'s oracle per label. (C2) Every seam curve is closed, and the non-manifold edges of the merged complex are exactly the triple-junction edges (`found = expected`) — the complex is non-manifold *only* where the geometry is. SHARE: none.

**Falsified by.** C1 by a χ mismatch on one material, which locates a wrong orbit in `P-236`'s table — the oracle finds the table's bug. C2 by non-manifold edges away from junctions, a seam-stitching bug. VACUITY CONTROL: two materials with no contact must produce `seam_curves = 0`.

#### P-238 — registered for R-243, before the harness: A null registered on purpose — multi-material extraction costs under 1.3× single-material on the same grid

**Ticket:** `R-243` (S). **Records:** `k`, `resolution`, `wall_ms_multi`, `wall_ms_single_m160`, `ratio`, `mixed_cells_share`, `c1_holds`.

**Hypothesis.** The registered null: (C1) with a lookup table, the per-cell cost is a wider index and one extra branch, so `ratio ≤ 1.3` at `k = 3` on the three-sphere fixture. SHARE: `mixed_cells_share` (cells with more than one label) is reported so the ratio is read against how often the multi-label path runs.

**Falsified by.** C1 by a larger ratio, which would say the table's cache footprint (`P-236`'s bytes) is the cost, and the row names the `k` where it starts. VACUITY CONTROL: `k = 1` must give `ratio = 1.00 ± 0.02`.

---

# Group C — planets: extraction on the cubed sphere and under diffeomorphism (curvilinear grids, differential topology) [feature]

**Field and words.** Ronchi–Iacono–Paolucci 1996 *The cubed sphere*; gnomonic projection; `curvilinear grid`,
`Jacobian of the map`, `topological invariance of sign structure`. **The result [classical]:** a
homeomorphism of the parameter domain preserves the sign of the field at every corner, so the *combinatorial*
output of the case table is invariant under the warp; only the *geometry* changes, by the map's Jacobian.
The trilinear is exact only for affine maps, so on the cubed sphere the extracted surface is the image of the
parameter-domain PL surface, with Hausdorff error bounded by `‖J‖·h²`-type terms. **The hook:** games want
spherical worlds; the crate meshes boxes. The cubed sphere has six cube-face charts, each a near-uniform grid,
and the seams between charts are the interesting part.

#### P-239 — registered for R-244, before the harness: Case indices are invariant under the cubed-sphere warp, so topology on the planet equals topology in the parameter cube

**Ticket:** `R-244` (S). **Records:** `field`, `resolution`, `chart`, `case_index_hash_param`, `case_index_hash_warped`, `equal`, `chi_param`, `chi_warped`, `c1_holds`.

**Hypothesis.** (C1) For a field defined on the sphere (`sphere` with radius varying by `fbm_terrain`'s heights, sampled through the gnomonic map), the per-cell case indices computed in the parameter cube and after the warp are identical on every chart, and `chi_warped = chi_param` — a theorem checked. SHARE: none.

**Falsified by.** C1 by any difference, which can only come from the field being sampled at *different* points (the warp applied to corners after sampling rather than before) — a fixture bug, but exactly the bug a planet implementation will have. VACUITY CONTROL: a warp with a fold (non-injective) must change case indices.

#### P-240 — registered for R-245, before the harness: Hausdorff error on the planet scales with the chart's Jacobian, and is worst at the cube's corners by the gnomonic distortion factor

**Ticket:** `R-245` (M). **Records:** `field`, `resolution`, `chart`, `cell`, `jacobian_norm`, `jacobian_norm_max_over_min`, `hausdorff_cell`, `hausdorff_over_h2`, `correlation_with_jacobian`, `corner_over_centre`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) Per-cell Hausdorff error correlates with `jacobian_norm` at `ρ ≥ 0.7` on a planet-sized `sphere`. (C2) The ratio of error at chart corners to chart centres is within **30%** of the gnomonic distortion ratio (classically about `1.3–1.4` for the equal-angle cubed sphere) — the error is the map, not the extractor. SHARE: none.

**Falsified by.** C1 by low correlation, which would locate the error in the warp of the *crossing positions* (interpolated in the parameter domain, then mapped) — fixable by mapping corners first and interpolating in world space, and the row measures both orders. C2 by a larger ratio, which would say the extractor's own anisotropy (`P-186`) compounds the map's. VACUITY CONTROL: the identity map must give `corner_over_centre = 1`.

#### P-241 — registered for R-246, before the harness: Chart seams are conforming if the six charts share their edge samples, and the seam's crack count is zero without any stitching code

**Ticket:** `R-246` (M). **Records:** `resolution`, `seam_edges`, `seam_vertices_shared`, `crack_edges`, `duplicate_vertices`, `chi_planet`, `chi_true`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) When adjacent charts sample the same world points on their shared edge (the cubed sphere's charts meet at cube edges with a common parameterisation), the two charts' meshes share vertices by construction and `crack_edges = 0`. (C2) `chi_planet = 2` for a planet with no caves and `chi_true` (from `P-145` run in world space) for one with. SHARE: none.

**Falsified by.** C1 by cracks at seams, which come from the two charts interpolating the crossing along the shared edge from *different* corner pairs — the cubed sphere's known subtlety at the eight cube corners where three charts meet, and the row reports whether the cracks are only there. VACUITY CONTROL: offsetting one chart's samples by `h/2` must produce cracks.

---

# Group D — silhouettes from the field: apparent contours (singularity theory) [feature]

**Field and words.** Whitney 1955 (folds and cusps of maps from surfaces to the plane); Koenderink 1984 *What
does the occluding contour tell us about solid shape?*; `apparent contour`, `contour generator`, `fold`,
`cusp`, `T-junction`. **The result [classical]:** the silhouette of a smooth surface under a generic view is
the image of the contour generator `{x : n(x)·v = 0}`, whose singularities are only folds and cusps
(Whitney), and the sign of the silhouette's curvature equals the sign of the surface's Gaussian curvature at
the generator (Koenderink). **The hook:** outline rendering and shadow silhouettes walk mesh edges looking
for normal-sign changes; the generator is a level set of `∇f·v` and the crate has a level-set extractor.

#### P-242 — registered for R-247, before the harness: The contour generator is the zero set of `∇f·v` on the surface, extractable by the same machinery as a curve, and its cusp count matches Whitney's genericity

**Ticket:** `R-247` (M). **Records:** `field`, `view_direction`, `resolution`, `generator_curves`, `generator_closed`, `cusps_found`, `cusps_by_koenderink_sign_change`, `mesh_silhouette_edges`, `hausdorff_generator_vs_mesh_silhouette`, `cost_share`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) Extracting `{f = 0} ∩ {∇f·v = 0}` as the intersection of two level sets (a curve on the surface — per-cell, both signs, a 1-D case table) gives closed generator curves on `sphere`, `torus`, `gyroid`, with `cusps_found` equal to the number of points where the generator's tangent aligns with `v`. (C2) The generator's projection is within `2h` Hausdorff of the mesh's silhouette edges, at a cost under **10%** of extraction (the second field is a dot product on already-computed gradients). SHARE: C2's cost is against one extraction.

**Falsified by.** C1 by open generator curves, which come from the two level sets being extracted from *different* interpolants — the per-cell curve must use the same trilinear as the surface. C2 by a larger cost, which would say the gradients were not already available (they are for DC, not for MC) and the row prices both paths. VACUITY CONTROL: `v` along a symmetry axis of `sphere` must give exactly one generator curve (the equator).

#### P-243 — registered for R-248, before the harness: A null registered on purpose — Koenderink's sign theorem holds on the mesh's silhouette only where `P-185`'s Gaussian curvature estimate has a sign

**Ticket:** `R-248` (S). **Records:** `field`, `view_direction`, `silhouette_vertex`, `silhouette_curvature_sign`, `gaussian_curvature_sign_p185`, `agree`, `agree_share`, `indeterminate_share`, `c1_holds`.

**Hypothesis.** The registered null: (C1) on `torus`, silhouette curvature sign agrees with the mesh-free Gaussian curvature sign at `≥ 90%` of silhouette vertices where `P-185`'s estimate is not within its own error of zero — a theorem about smooth surfaces surviving discretisation on the vertices where the discretisation has an opinion. SHARE: `agree_share` is of *determinate* silhouette vertices; `indeterminate_share` is reported.

**Falsified by.** C1 by low agreement, which would say the silhouette walked on the mesh is not the contour generator — the mesh's silhouette includes staircase edges (`P-186`) that have no smooth counterpart — and would make `P-242`'s field-derived generator the *only* correct outline. VACUITY CONTROL: `sphere` must give 100% positive on both.

---

# Group E — conservative offset meshes: occluders and hulls by isovalue shift (mathematical morphology, Lipschitz analysis) [feature]

**Field and words.** Minkowski erosion/dilation, offset surfaces, `1-Lipschitz`, `conservative occluder`,
`bounding hull`. **The result [classical]:** for a signed distance function, `{f ≤ −ε}` is exactly the
`ε`-erosion of the solid and `{f ≤ +ε}` exactly its `ε`-dilation; for an `L`-Lipschitz field that is not a
distance function, `{f ≤ −Lε}` is *contained in* the erosion and `{f ≤ +Lε}` *contains* the dilation. **The
hook:** software occlusion culling needs a mesh strictly inside the solid; broad-phase collision needs one
strictly outside; both are one isovalue away, and the Lipschitz bound the crate already computes for
empty-cell rejection is the constant.

#### P-244 — registered for R-249, before the harness: The `−Lε` isosurface is a conservative occluder — inside the true solid everywhere — and the `+Lε` one a conservative hull, verified by the oracle

**Ticket:** `R-249` (M). **Records:** `field`, `is_sdf`, `lipschitz_L`, `epsilon`, `resolution`, `occluder_vertices_outside_solid`, `hull_vertices_inside_solid`, `occluder_volume_over_true`, `hull_volume_over_true`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) On every field, extracting at `−Lε` with `ε = 2h` gives a mesh with `occluder_vertices_outside_solid = 0` (each vertex tested against the field, and against the mesh's own Hausdorff bound from `M-472`), and at `+Lε` a mesh with `hull_vertices_inside_solid = 0`. (C2) On the SDF fields (`is_sdf = true`) the occluder's volume ratio is within **5%** of the analytic erosion volume (`V − ε·A + ε²·M − …`, the Steiner formula the mean-breadth instrument `P-185` already uses), so the shift is not merely conservative but *tight*. SHARE: none.

**Falsified by.** C1 by a vertex on the wrong side, which would mean `ε < 2h` is below the extractor's own error — the row raises `ε` until it holds and reports the multiple of `h` that suffices, which is the conservative offset a game must use. C2 by a loose volume on an SDF, which would say the crossing interpolation biases inward (Phase 30's signed-defect registration's inward bias) more than `ε` absorbs. VACUITY CONTROL: `ε = 0` must reproduce the ordinary mesh and fail C1 on at least one field.

#### P-245 — registered for R-250, before the harness: The occluder needs fewer triangles than the render mesh at the same resolution, because erosion removes features smaller than `ε`

**Ticket:** `R-250` (S). **Records:** `field`, `epsilon_over_h`, `triangles_occluder`, `triangles_render`, `ratio`, `components_occluder`, `components_render`, `local_thickness_min_phase30`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) On `noise_cavity` and `thin_plate`, `ratio ≤ 0.7` at `ε = 2h` — thin features vanish under erosion and take their triangles with them. (C2) `components_occluder ≤ components_render`, with the difference equal to the number of components whose local thickness (Phase 30's local-thickness registration) is below `2ε`. SHARE: `ratio` is over whole meshes.

**Falsified by.** C1 by `ratio ≈ 1`, which would say the roster's features are all thicker than `2ε` — a statement about the roster, and `P-223` records it. C2 by more components after erosion, which is possible (erosion can disconnect) and the row records it as the feature's known behaviour. VACUITY CONTROL: `ε = 0` must give `ratio = 1`.

---

# Group F — bit-identical extraction across threads and architectures (floating-point determinism) [feature][perf]

**Field and words.** IEEE 754 reproducibility; `fused multiply-add contraction`, `reduction order`,
`canonical order`, `lockstep determinism`. **The result [classical]:** floating-point addition is not
associative, so any parallel reduction whose order depends on thread scheduling is non-deterministic; and
`arm64` and `x86-64` differ in whether `a*b+c` contracts to an FMA unless the code forbids it. **The hook:**
lockstep multiplayer and replay require every machine to produce the *same* mesh from the same edits; `M-32`
showed the crate's seam vertices differ by up to `1.57e-16` world units when the same point is reached by
two algebraically equal expressions *within one machine*. `M-32`, `M-175` and `M-31`'s 378 cross-platform
golden hashes are what the repo has already measured — two expressions for one seam point, a sum that
depends on tie order, and bit-identity of one machine's output across macOS/arm64 and Linux/x86-64. None of
them varies the **thread count** or hashes the **same edit on two different machines**, which is what this
group's two rows do.

#### P-246 — registered for R-251, before the harness: A null registered on purpose — with a canonical (Morton) emission order the mesh hash is invariant across thread counts on one machine

**Ticket:** `R-251` (S). **Records:** `field`, `resolution`, `thread_count`, `mesh_hash`, `hashes_equal_across_threads`, `c1_holds`.

**Hypothesis.** The registered null: (C1) the current extractor, run at 1, 4, 8 and 16 threads with output sorted into Morton order before hashing, produces identical hashes on every field — thread count does not reach the arithmetic. SHARE: none.

**Falsified by.** C1 by differing hashes, which would locate a reduction whose order depends on scheduling (a QEF accumulation, a vertex-dedup hash-map iteration) and would be a real bug for any multiplayer use. VACUITY CONTROL: hashing *without* the canonical sort must differ across thread counts, or the sort is not the thing making them equal.

#### P-247 — registered for R-252, before the harness: Across `mac-air` (arm64) and `big` (x86-64) the hashes differ, and forbidding FMA contraction plus a fixed reduction order makes them equal at under 5% cost

**Ticket:** `R-252` (M). **Records:** `field`, `resolution`, `hash_arm64`, `hash_x86`, `equal_default`, `equal_no_fma_fixed_order`, `differing_vertices_default`, `max_vertex_delta_default`, `wall_ms_default`, `wall_ms_deterministic`, `cost_share`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** (C1) With default settings the two machines' hashes differ on at least one field, and `max_vertex_delta_default` is at the `ulp` scale of `M-32` — not a bug, an architecture. (C2) With `mul_add` replaced by explicit `mul` then `add` (or all arithmetic routed through `mul_add`) *and* the crossing evaluation order fixed as `P-224`'s backward-error row prescribes, `equal_no_fma_fixed_order = true` on every field. (C3) `cost_share ≤ 5%` of extraction. SHARE: C3 is against `M-160`'s per-chunk time on each machine separately.

**Falsified by.** C1 by equal hashes, which would say the crate already avoids contraction — good, and the row documents why. C2 by remaining differences, which would locate them in a library (`libm` transcendental in the noise fields) and the row names the function; a game then ships its own. C3 by a larger cost, which prices determinism as a build flag rather than a default. VACUITY CONTROL: a deliberate `mul_add` in the crossing path must break equality.

---

# Group G — automatic UVs for voxel terrain: discrete conformal equivalence (discrete conformal geometry) [feature]

**Field and words.** Springborn–Schröder–Pinkall 2008 *Conformal equivalence of triangle meshes*; discrete
Ricci flow (Gu–Yau); `discrete conformal map`, `quasi-conformal dilatation`, `cone singularity`, `seam`,
`atlas`. **The result [classical]:** a triangle mesh admits a discrete conformal flattening obtained by
minimising a convex energy in per-vertex log scale factors, with prescribed cone angles; the flattening's
distortion is measured by the per-triangle dilatation. **The hook:** voxel terrain is textured by tri-planar
projection (three samples, blended) because it has no UVs; an atlas per chunk with bounded dilatation would
cost one sample.

#### P-248 — registered for R-253, before the harness: A per-chunk discrete conformal flattening of the extracted mesh has median dilatation under 1.2 and is computed at under 3× extraction cost

**Ticket:** `R-253` (L). **Records:** `field`, `chunk`, `resolution`, `cone_vertices`, `flattening_converged`, `dilatation_median`, `dilatation_p95`, `seam_length_over_perimeter`, `wall_ms_flatten`, `wall_ms_extract`, `cost_share`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** (C1) On every chunk of every field, the convex energy minimisation converges (Newton on the SSP energy) with cone singularities placed by Gaussian-curvature clustering, and `dilatation_median ≤ 1.2`, `p95 ≤ 2.0`. (C2) `seam_length_over_perimeter ≤ 1.5` — the atlas cuts are mostly the chunk boundary, which is a seam already. (C3) `cost_share ≤ 3×` extraction. SHARE: C3 is a multiple of one chunk's extraction.

**Falsified by.** C1 by non-convergence, which happens when the mesh has near-degenerate triangles — `M-309`'s needle ratio on `box_exact` is the suspect, and the row reports whether Phase 30's edge-flip registration's post-flip mesh converges where the raw one does not. C2 by long seams, which prices the atlas against tri-planar's zero seams. C3 by a larger cost, which moves flattening to a background pass. VACUITY CONTROL: a planar chunk must flatten with dilatation exactly 1.

#### P-249 — registered for R-254, before the harness: A null registered on purpose — one-sample atlas texturing and three-sample tri-planar differ by under 10% in fragment cost on the render path

**Ticket:** `R-254` (S). **Records:** `field`, `scene`, `fragment_ms_triplanar`, `fragment_ms_atlas`, `ratio`, `seam_visible_pixels_share`, `c1_holds`.

**Hypothesis.** The registered null: (C1) on the Bevy example scene, `ratio ≥ 0.9` — the GPU's texture cache makes the three samples nearly free, and the atlas's benefit is *quality* (no blend smearing), not speed. SHARE: `fragment_ms` is the fragment stage alone, from GPU timestamps.

**Falsified by.** C1 by `ratio < 0.9`, which would make the atlas a performance feature as well, and the row records the bandwidth-bound scene where it happens. VACUITY CONTROL: a scene with the texture cache defeated (random UV offsets) must show `ratio ≈ 1/3`.

---

# Group H — geodesic distance as an extraction by-product: the heat method (heat-kernel PDE) [feature]

**Field and words.** Crane–Weischedel–Wardetzky 2013 *Geodesics in heat*; Varadhan's formula; `heat method`,
`prefactored Laplacian`, `geodesic distance field`. **The result [classical]:** geodesic distance from a
source on a mesh is obtained by one short-time heat diffusion, normalising the gradient, and one Poisson solve
— two sparse linear solves with a matrix that can be prefactored once per mesh. **The hook:** pathfinding and
"distance to the player" over terrain currently run on a navmesh graph (`✗67`'s octree); the extracted mesh's
cotan Laplacian (already built for Phase 30's cotan-Laplacian registration) gives geodesic distance in two
solves per query.

#### P-250 — registered for R-255, before the harness: Heat-method geodesics on the extracted mesh are within 3% of true geodesics on `sphere` and `torus`, and a prefactored chunk answers a query in under 1 ms

**Ticket:** `R-255` (M). **Records:** `field`, `resolution`, `source`, `target`, `geodesic_heat`, `geodesic_true`, `relative_error`, `t_parameter_over_h2`, `prefactor_ms`, `query_ms`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) On `sphere` (great-circle distance) and `torus` (numerically integrated geodesic), `relative_error ≤ 3%` at 65³ with `t = h²`. (C2) `query_ms ≤ 1` per chunk after a `prefactor_ms` reported once — a distance field per frame is affordable. SHARE: none.

**Falsified by.** C1 by larger error, which on an MC mesh is the needle-triangle conditioning of the cotan matrix (`M-309`) and the row reports error on Phase 30's edge-flip registration's post-flip mesh alongside. C2 by slow queries, which prices the feature per chunk size. VACUITY CONTROL: source equal to target must give distance `0 ± h`.

#### P-251 — registered for R-256, before the harness: Geodesic distance from the *field* — the heat method on the digital surface — agrees with the mesh's, so navigation needs no mesh either

**Ticket:** `R-256` (M). **Records:** `field`, `resolution`, `source`, `target`, `geodesic_mesh_p250`, `geodesic_digital`, `relative_difference`, `c1_holds`.

**Hypothesis.** [abstract — Caissard et al. 2018, in the corpus at `10.1007/s10851-018-0839-4`, 0.6766 at this phase's probe] The digital Laplace–Beltrami operator of Phase 30's digital-operator registration supports the same two solves on the voxel boundary. (C1) `relative_difference ≤ 5%` on every field at 65³ — a navigation distance field that is a by-product of the sign grid, before any mesh exists, which is what a server without a renderer wants. SHARE: none.

**Falsified by.** C1 by larger differences on the ungradeable fields, which is a finding about which of the two is closer to the truth — and on `sphere` the row can say which. VACUITY CONTROL: the digital operator with the integration radius set to zero must fail to converge.

---

# Group I — importing meshes into the voxel world: generalized winding numbers (computational topology) [feature][perf]

**Field and words.** Jacobson–Kavan–Sorkine 2013 *Robust inside-outside segmentation using generalized winding
numbers*; Barill et al. 2018 *Fast winding numbers for soups and clouds*; `generalized winding number`, `solid
angle`, `Barnes–Hut dipole expansion`. **The result [classical]:** the winding number `w(p) = (1/4π) Σ solid
angles` is an integer inside a closed mesh, fractional near holes, and robust to self-intersections and open
boundaries — a signed indicator that degrades gracefully; with a dipole-tree expansion it is `O(log n)` per
query. **The hook:** the crate goes field → mesh; a game also needs mesh → field (drop a prop, carve it into
the world), and `P-230`'s volume test is the oracle for whether the round trip closes. Re-scoped at this
phase's probe: the crate already ships both directions — `crates/isomesh/src/construct/winding.rs` (S-007) and
`crates/isomesh/src/construct/from_mesh.rs` (S-006) — and `M-262`, `M-299` and `M-303` priced them (the
winding number beats the pseudonormal on holed meshes; the on-demand split costs a factor of N; the crossover
is N² with per-point cost linear in boundary edges). What is unmeasured is S-006's own acceptance criterion —
mesh a sphere, convert back to a field, re-mesh, compare — which has no Hausdorff or χ number in
`FINDINGS.md`, and which is this group's first row.

#### P-252 — registered for R-257, before the harness: The round trip `field → mesh → winding number → field → mesh` closes to within `2h` Hausdorff and preserves `χ`

**Ticket:** `R-257` (M). **Records:** `field`, `resolution`, `hausdorff_roundtrip`, `hausdorff_over_h`, `chi_original`, `chi_roundtrip`, `winding_fractional_cells_share`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) On every closed-mesh field, re-voxelising the extracted mesh by `w(p) ≥ ½` and re-extracting gives `hausdorff_over_h ≤ 2` and `chi_roundtrip = chi_original`. (C2) `winding_fractional_cells_share` (cells whose corner values are not near 0 or 1) is under **2%** — the extracted mesh is closed enough that the winding number is essentially an indicator. SHARE: none.

**Falsified by.** C1 by χ changing, which would locate a non-closed mesh (a chunk seam with a crack, `✗-class` bug) that the winding number *tolerates* — which is the point of using it, and the row reports where the fraction lives. C2 by a large fractional share, the same finding stated the other way. VACUITY CONTROL: deleting 5% of triangles before re-voxelising must raise the fractional share and *not* change χ, or the winding number is not doing its job.

#### P-253 — registered for R-258, before the harness: The Barnes–Hut winding number voxelises a 100k-triangle prop into a 65³ chunk in under 50 ms, against seconds for the exact sum

**Ticket:** `R-258` (S). **Records:** `mesh_triangles`, `resolution`, `wall_ms_exact`, `wall_ms_fast`, `speedup`, `max_abs_winding_error`, `sign_flips_vs_exact`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) `speedup ≥ 20×` at 100k triangles and 65³. (C2) `sign_flips_vs_exact = 0` — the approximation never changes a corner's side, because the error is far from `½` away from the surface and the near-field is evaluated exactly. SHARE: none.

**Falsified by.** C1 by less speedup, which prices the tree's construction against the query count and the row reports the crossover. C2 by sign flips, which are at cells within `h` of the surface and the row reports whether the exact near-field radius was set below `h`. VACUITY CONTROL: a mesh of 100 triangles must give `speedup < 2`.

---

# Group J — compressing the sign grid: sparse voxel DAGs against the entropy floor (data structures, information theory) [perf]

**Field and words.** Kämpe–Sintorn–Assarsson 2013 *High resolution sparse voxel DAGs*; `hash-consing`,
`subtree deduplication`, `symmetry-aware DAG`. **The result [classical]:** an octree of a binary volume whose
identical subtrees are merged into a DAG compresses geometry by orders of magnitude where the volume is
self-similar. **The hook:** `P-204` (`R-209`) measured the sign grid's entropy rate as the lossless floor for
bit-packing; a DAG is the *other* way to approach that floor, and it is queryable in place.

#### P-254 — registered for R-259, before the harness: A null registered on purpose — the DAG compresses the sign grid to within 4× of the entropy-rate floor on smooth fields and cannot approach it on noise

**Ticket:** `R-259` (M). **Records:** `field`, `resolution`, `bits_dag`, `bits_entropy_floor_p204`, `ratio_over_floor`, `subtrees_merged_share`, `bits_bitpacked_phase25`, `c1_holds`, `c2_holds`.

**Hypothesis.** The registered null: (C1) on `sphere`, `torus`, `box_exact`, `ratio_over_floor ≤ 4` and below Phase 25's bit-packed size. (C2) On `noise_cavity` the DAG is *larger* than the bit-packed grid — dedup finds nothing and the pointers are overhead. SHARE: none.

**Falsified by.** C1 by a larger ratio on smooth fields, which would say the octree's node overhead dominates at 65³ and the row reports the resolution where the DAG wins. C2 by the DAG winning on noise, which would be surprising and would say value noise at one octave has more repeated subtrees than its entropy rate suggests. VACUITY CONTROL: a constant field must compress to one node.

#### P-255 — registered for R-260, before the harness: In-place surface-cell enumeration on the DAG is faster than the bit-packed scan at 257³ because it skips merged empty subtrees

**Ticket:** `R-260` (S). **Records:** `field`, `resolution`, `surface_cells`, `wall_ms_dag_enumerate`, `wall_ms_bitpacked_scan_m306`, `ratio`, `nodes_visited`, `nodes_total`, `c1_holds`.

**Hypothesis.** (C1) At 257³ on `sphere`, enumerating surface cells by DAG traversal visits under **10%** of nodes and is faster than `M-306`'s rejection scan by at least **2×**. SHARE: the ratio is of the enumeration stage.

**Falsified by.** C1 by no speedup, which would say `M-306`'s coarse rejection is already skipping the same volume and the DAG's pointer chasing costs what the skip saves. VACUITY CONTROL: at 33³ the DAG must *lose*, or the traversal is not paying its overhead.

---

# Group K — a 6-byte vertex: quantised crossing parameters (rate–distortion, lattice quantisation) [perf]

**Field and words.** `rate–distortion`, `uniform scalar quantiser`, `edge-indexed vertex`, `vertex shader
decode`. **The result [classical]:** a crossing position on grid edge `e` at parameter `t ∈ [0,1]` quantised
to `b` bits has error `≤ h/2^{b+1}` along the edge and zero across it; a vertex is then `(edge_id: u32, t:
u16)` = 6 bytes, decoded in the vertex shader, against 12 bytes for `f32×3` (plus 12 for a normal). **The
hook:** Phase 25 compressed the *input*; the *output* vertex buffer is the bandwidth the GPU actually moves,
and its size has never had a floor.

#### P-256 — registered for R-261, before the harness: A null registered on purpose — 16-bit `t` is below the extractor's own error at every resolution on the roster

**Ticket:** `R-261` (S). **Records:** `field`, `resolution`, `bits_t`, `quantisation_error_bound`, `hausdorff_m472`, `bound_over_hausdorff`, `hausdorff_quantised`, `c1_holds`.

**Hypothesis.** The registered null: (C1) `quantisation_error_bound / hausdorff_m472 ≤ 0.01` at `b = 16` on every row, and `hausdorff_quantised` equals `hausdorff_m472` to three digits — quantisation is invisible. SHARE: none.

**Falsified by.** C1 by a visible change, which can only happen at resolutions above 4096³ where `h/2^17` competes with `h²` — and the row states that crossover. VACUITY CONTROL: `b = 4` must raise Hausdorff measurably.

#### P-257 — registered for R-262, before the harness: The 6-byte vertex halves vertex-buffer bandwidth and does not slow the vertex stage

**Ticket:** `R-262` (M). **Records:** `field`, `scene`, `vertex_bytes_f32`, `vertex_bytes_quantised`, `ratio`, `vertex_stage_ms_f32`, `vertex_stage_ms_quantised`, `decode_cost_ratio`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) `ratio ≤ 0.5` including a 4-byte octahedral-encoded normal in the quantised format against `f32×3 + f32×3`. (C2) `decode_cost_ratio ≤ 1.1` — the two integer-to-float conversions in the shader are hidden under memory latency. SHARE: `vertex_stage_ms` is from GPU timestamps, the vertex stage alone.

**Falsified by.** C1 by a ratio above 0.5, an arithmetic error. C2 by a slower vertex stage, which would say the scene is ALU-bound, and the row names it. VACUITY CONTROL: a scene with one triangle must show `decode_cost_ratio ≈ 1`.

---

# Group L — traversal order with a theorem: Hilbert against Morton (space-filling curve analysis) [perf]

**Field and words.** Moon–Jagadish–Faloutsos–Saltz 2001 *Analysis of the clustering properties of the Hilbert
space-filling curve*; `clustering`, `number of runs`, `locality`. **The result [classical]:** for a query
region, the expected number of contiguous runs along the Hilbert curve is provably lower than along the
Z-order (Morton) curve, with the constant derived. **The hook:** `M-306`'s scan and Phase 30's
memory-traffic-floor registration are measured on Morton; Hilbert has a theorem saying it should touch fewer
cache lines.

#### P-258 — registered for R-263, before the harness: Hilbert order reduces cache-line runs per chunk by the Moon et al. constant and wall time by less than 5%

**Ticket:** `R-263` (S). **Records:** `field`, `resolution`, `order`, `runs_per_chunk`, `runs_ratio_hilbert_over_morton`, `moon_predicted_ratio`, `l2_misses`, `wall_ms`, `wall_ratio`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) `runs_ratio` is within **20%** of the theorem's prediction for a `2×2×2` stencil window. (C2) `wall_ratio ∈ [0.95, 1.0]` — the theorem is right about runs and the runs do not matter at chunk size, because the chunk fits in L2 (Phase 30's "the chunk fits in L2" reading). SHARE: none.

**Falsified by.** C1 by a ratio far from prediction, which would mean the `runs` were counted over the wrong window. C2 by a real speedup, which would say the working set exceeds L2 at 65³ and the row reports the resolution boundary — exactly the number that registration needs. VACUITY CONTROL: lexicographic order must have more runs than either.

---

# Group M — destruction with prescribed fragment volumes: power diagrams (computational geometry, semi-discrete optimal transport) [feature]

**Field and words.** Aurenhammer–Hoffmann–Aronov 1998 *Minkowski-type theorems and least-squares clustering*;
`power diagram`, `Laguerre tessellation`, `weights`, `prescribed volumes`. **The result [classical]:** for any
set of sites and any prescribed volumes summing to the solid's volume, there exist weights making the power
diagram's cells have exactly those volumes, found by a concave maximisation. **The hook:** Voronoi fracture
gives fragments whose sizes are whatever they are; a designer wants "one big chunk and five small ones". OT
was mined in Phase 27 for a different purpose (transport of the error field); the semi-discrete case as a
*fracture* control is new to the ledger and the row says so.

#### P-259 — registered for R-264, before the harness: Fragment volumes are controllable to 2% by power-diagram weights, and each fragment's mesh is the solid clipped by the cell

**Ticket:** `R-264` (M). **Records:** `field`, `sites`, `prescribed_volumes`, `weights_converged`, `newton_iterations`, `fragment_volume_p230`, `volume_error_share`, `fragment_chi`, `fragments_closed`, `total_volume_conserved`, `c1_holds`, `c2_holds`.

**Hypothesis.** (C1) On `sphere` with 8 sites and volumes `(50%, 7×~7%)`, the concave maximisation converges in under 30 Newton steps and every fragment's `P-230` volume is within **2%** of prescribed. (C2) Each fragment (field `∧` power cell, extracted with the cell's planes as additional half-space fields) is a closed mesh with `chi = 2` and the fragment volumes sum to the solid's within extraction error. SHARE: `volume_error_share` is per fragment against its own prescribed volume.

**Falsified by.** C1 by non-convergence, which happens when a prescribed volume exceeds what its site can reach — the theorem needs volumes in the feasible set, and the row records the feasibility test. C2 by a non-closed fragment, a clipping seam bug at the power-cell face. VACUITY CONTROL: equal prescribed volumes with equal weights must reproduce the plain Voronoi fracture.

#### P-260 — registered for R-265, before the harness: The clipped-fragment extraction costs under 2× the unfragmented mesh for 8 fragments, because each cell's planes touch few grid cells

**Ticket:** `R-265` (S). **Records:** `sites`, `wall_ms_fragmented`, `wall_ms_whole`, `ratio`, `cells_touched_by_any_plane_share`, `c1_holds`.

**Hypothesis.** (C1) `ratio ≤ 2` at 8 sites and 65³; `cells_touched_by_any_plane_share ≤ 15%` — the planes are sparse in the grid and the rest of the extraction is unchanged. SHARE: the ratio is of whole extractions; the touched share explains it.

**Falsified by.** C1 by a larger ratio, which prices fragmentation per site and the row fits the linear model. VACUITY CONTROL: one site must give `ratio = 1.00 ± 0.02`.

---

# Group N — ridge-preserving pre-filter: anisotropic diffusion of the field (PDE image processing) [feature]

**Field and words.** Perona–Malik 1990 *Scale-space and edge detection using anisotropic diffusion*; `edge-
stopping function`, `anisotropic diffusion`, `scale space`. **The result [classical]:** diffusing a field with
conductance that vanishes at large gradient smooths noise inside regions while sharpening, not blurring, the
boundaries between them. **The hook:** `✗116` found the tricubic pre-filter blurs `box_exact`'s creases; Phase
30's WENO crossing registration attacks the crossing; this attacks the *field*, so terrain noise smooths and
cliff edges stay.

#### P-261 — registered for R-266, before the harness: A null registered on purpose — Perona–Malik on `fbm_terrain` removes octave-3 noise triangles and does not move the ridge lines

**Ticket:** `R-266` (M). **Records:** `field`, `iterations`, `kappa`, `triangles_before`, `triangles_after`, `triangle_ratio`, `ridge_displacement_max`, `ridge_displacement_over_h`, `hausdorff_to_unfiltered`, `c1_holds`, `c2_holds`.

**Hypothesis.** The registered null: (C1) 10 iterations at `κ` set to the octave-2 gradient magnitude reduce triangle count by at least **20%** on `fbm_terrain`. (C2) Ridge lines (the crest lines of `P-242`'s machinery with `v` vertical) move by `ridge_displacement_over_h ≤ 0.5` — the filter is anisotropic where it matters. SHARE: `triangle_ratio` is over the whole mesh.

**Falsified by.** C1 by a small reduction, which would say the field's noise is below `κ` everywhere and the filter is isotropic in practice — a `κ` calibration finding. C2 by ridges moving, which would say `κ` was set below the ridge gradient and the row reports the `κ` at which they stop. VACUITY CONTROL: `κ = ∞` must reduce to Gaussian blur and move ridges by more than `h`.

#### P-262 — registered for R-267, before the harness: On `box_exact` the diffused field's creases are sharper than the tricubic's, measured by Phase 30's smoothness indicator

**Ticket:** `R-267` (S). **Records:** `filter`, `crease_edges_flagged_phase30_indicator`, `crease_edges_true`, `recall`, `crease_sharpness_beta_mean`, `hausdorff_at_crease_cells`, `c1_holds`.

**Hypothesis.** (C1) After Perona–Malik, `recall` of true crease edges by Phase 30's smoothness-indicator registration is above **0.9** and `hausdorff_at_crease_cells` is below the tricubic filter's (`✗116`) by at least **30%** — the crease survives the pre-filter that `✗116`'s did not. SHARE: none.

**Falsified by.** C1 by no improvement, which would say the diffusion's stopping function saturates at the SDF's constant gradient magnitude (an SDF has `|∇f| = 1` everywhere, so gradient *magnitude* cannot distinguish crease from face) — and the row switches the conductance to the *Hessian* norm, which does. That is the expected outcome on an SDF, and it is registered here because it is the right diagnosis. VACUITY CONTROL: `sphere` must be unchanged by the filter to within `0.1h`.

---

## What this phase does not claim

- No row claims the feature is *shippable*; each claims a measurable property of a prototype on the roster.
- The corpus probe was run at landing (§1); a group whose probe came back present-and-cited in the repo is
  re-scoped there to what the repo did not measure.
- Group M's field overlaps Phase 27's optimal transport row in name only; the row states the difference.
- Groups G and H produce assets a *renderer* and a *navigation* system consume; they are registered here
  because each is a by-product of the mesh the crate already builds, and each is graded against a truth the
  crate already computes (Phase 30's cotan Laplacian, `sphere`'s great circles).
- Group F is the only group whose failure is a bug rather than a finding, and `P-246` is registered as a
  null for that reason: it is the cheapest test in the phase and the most expensive to be wrong about.
