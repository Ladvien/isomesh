# isomesh — findings ledger

**Started:** 2026-08-11 · **Append-only.** Entries are never deleted, only re-tiered when new evidence
arrives — with the old verdict left visible.

This is the project's epistemic state: what we believe, how strongly, and **on what evidence**. It
exists because this project has already been wrong six times in ways that would have silently
propagated into code, and because the research corpus contains several published figures that failed
verification. A belief with no recorded falsification method is not a finding, it's a preference.

## Machines

Every measured figure in this repo comes from one of these two, and machine-specific figures say
which (the README and demo pages lean on this block by reference; added at D-003, sourced from
`docs/research/2026-08-13-measured-comparison.md`).

| | |
|---|---|
| **M5** | Apple M5 · macOS · arm64 · Metal — `docs/measurements/resolution_sweep.csv` |
| **Zen 3** | AMD Ryzen 9 5900X · Linux · x86-64 · NVIDIA RTX 3090 (Vulkan) — `docs/measurements/resolution_sweep-ryzen9-5900x.csv`, and every GPU figure |

## How to use this file

- **Before acting on a "known" fact, look for it here.** If it's not here, it hasn't been checked.
- **When a measurement contradicts something written down, that's an entry** — the contradiction is
  the finding, not an inconvenience.
- **Falsified entries stay.** They're the most valuable rows, because they tell you which *sources* to
  distrust, not just which facts.
- Every entry names how it could be shown wrong. If you can't write that line, you have an opinion.
- **An experiment that was run and reverted is an entry too** — Part 4b. The verdict line is the
  point: a variant that was built, measured and put back is the most expensive thing to rediscover,
  because nothing left in the tree records that it was ever tried.

## Index

<!-- BEGIN GENERATED INDEX -- scripts/findings_index.sh -->

**358 entries** — 19 falsified, 285 measured, 34 verified, 16 open, 4 experiments. Regenerate with `scripts/findings_index.sh`; CI fails if this is stale.

| # | Claim |
|---|---|
| `✗1` | "Surface Nets produces substantially fewer triangles than Marching Cubes" |
| `✗2` | "You can have a manifold mesh or an intersection-free one, not both" |
| `✗3` | "Every interior Surface Nets vertex has four neighbours" |
| `✗4` | "Dual Contouring is absent from the home-still corpus" |
| `✗5` | "naga_oil is the shader composition path for Bevy" |
| `✗6` | "Mesh shaders aren't reachable from inside Bevy" |
| `✗7` | "CBT 5.78 → 0.40 ms is from Dupuy 2020" |
| `✗8` | Velo3D "93M vertices / 31M faces" |
| `✗9` | "Marching Cubes' cost inside the volumetric loop was never measured" / "navmesh rebuild cost was never measured" |
| `✗10` | "glam should be the internal math library from day one" |
| `✗11` | "Plain Marching Cubes has ambiguous faces and produces holes" |
| `✗12` | "The equivariant vertex rule needs a fast three-plane path with a fallback" |
| `✗13` | "Real::as_f32 is the only narrowing operation in the crate, and the crate itself never calls it" |
| `✗14` | "Surface Nets is the cheapest thing in the family and the natural default" |
| `✗15` | "Marching Cubes is unconditionally manifold" |
| `✗16` | "glam 0.32 lands with A-007's vertex solve" |
| `✗17` | "Only the interior test can make Marching Cubes 33 non-manifold, so the face decider alone cannot" |
| `✗18` | "A hairline seam difference is a crack no weld can close" |
| `✗19` | "Manifold Dual Contouring's uniform-grid surface is always a manifold" |
| `M-1` | surface cells = crossed edges + χ |
| `M-2` | V_sn = V_mc + χ, F_sn = F_mc + 2χ |
| `M-3` | Surface Nets max vertex degree 10; Marching Cubes 9 |
| `M-4` | Surface Nets is non-manifold where one cell carries two sheets: 48 non-manifold edges on capped gyroid, 15 on fbm_terrai… |
| `M-5` | On box_exact, Surface Nets' nearest vertex to the corner (1,1,1) is 1.15 cells away |
| `M-6` | libm::sqrtf lowers to hardware fsqrt (aarch64+neon) / sqrtss (x86-64+sse2) |
| `M-7` | dev-dependencies do not propagate: consumer resolves 3 packages, the crate's own lockfile has 137 |
| `M-8` | 317 packages, both 29.0.4 and 30.0.0 |
| `M-9` | Workspace feature unification leaks: -p isomesh alone gives glam libm; whole-workspace gives it std, serde, bytemuck, en… |
| `M-10` | Unit sphere at 64³ (h = 4/63), symmetric Hausdorff: Marching Cubes 1.380e-3, Surface Nets 2.288e-3. |
| `M-11` | T-003's own acceptance criterion is loose by ~80×. |
| `M-12` | Marching Cubes' error falls like h², measured. |
| `M-13` | Surface cells ≈ 1.5·A/h², not A/h². |
| `M-14` | The reverse direction finds defects the forward direction structurally cannot. |
| `M-15` | Surface Nets' non-manifoldness is a resolution effect, not a topology effect. |
| `M-16` | The even-χ parity check is not independent of manifoldness — it is a corollary of it. |
| `M-17` | A case-table entry naming an uncut edge is caught inside the crate, before any mesh exists |
| `M-18` | (refined by T-008 — the arithmetic below is about adjacent cells; the effect on a real mesh is gradual, because a mesh w… |
| `M-19` | There is no meaningful fixed cost on the CPU extraction path, and the prediction saying so was written down before the r… |
| `M-20` | Marching Cubes' marginal cost is 4.75 ns/sample — 211 M samples/s, single-threaded, f32, Apple M5. |
| `M-21` | Surface Nets is not O(n³) over this range; Marching Cubes is. |
| `M-22` | ✗1's identity holds at every resolution to 256³ |
| `M-23` | f64 costs 8–10% on extraction paths with no matrix solve in them. |
| `M-24` | Bit-exact lattice equivariance needs magnitude-ordered products, not just sums. |
| `M-25` | The sharp-feature solve is nearly free: Dual Contouring costs 3% more than Surface Nets. |
| `M-26` | Dual Contouring reaches a box corner to 0.01 cells where Surface Nets stops at 0.58 |
| `M-27` | The two dual methods differ only at features, with a 14-order-of-magnitude gap and nothing in between. |
| `M-28` | The cell clamp eliminates placement-caused self-intersections entirely, and costs nothing in sharpness. |
| `M-29` | The literature's two branches both fire, on disjoint fields — which is a sharper answer than either alone. |
| `M-30` | An unclamped solve can fling a vertex 3.18 cells out of its own cell |
| `M-31` | libm delivers the bit-identical cross-platform meshes it was chosen for — verified, not reasoned. |
| `M-32` | Chunk seams are bit-exact only when the cell size is a power of two. |
| `M-33` | E1, the unpublished number: a brush changes 15–36% of the cells in its own bounding box. |
| `M-34` | Counting value changes overstates the re-mesh set by 2.8–3.7×. |
| `M-35` | A brush stepping two cells sweeps 6–14% of its sign-changed cells entirely through. |
| `M-36` | The multiplayer story survives, with a boundary. Eight brushes in all 8! = 40,320 orderings give exactly one result — bi… |
| `M-37` | Mixed add and subtract do not commute: 11 distinct results from the same 40,320 orderings. |
| `M-38` | Smooth union destroys reordering almost completely: 40,317 distinct results from 40,320 orderings. |
| `M-39` | Fixed-point storage is unnecessary for the guarantee it was proposed for. |
| `M-40` | The ambiguous face is rarer than the literature suggests — on five of the seven reference fields it never occurs at all. |
| `M-41` | 88 of the 256 cases change their Euler characteristic when their ambiguous faces are joined. |
| `M-42` | The asymptotic decider is free to within a few percent, which is the first time this repo's "~free" claim has had a benc… |
| `M-43` | The decider needs no division and no epsilon, and the brief's "guard the denominator" is unnecessary. |
| `M-44` | The decider does not widen M-32's chunk-seam problem, and there is margin to spare. |
| `M-45` | ✗14 reproduces on a second machine and gets worse; its crossover does not reproduce at all. |
| `M-46` | A chord is only collidable when its two cut edges share a cube face, and that makes the manifold fix nearly free. |
| `M-47` | The validator's duplicate_vertices is an upper bound on what a weld removes, not the count. |
| `M-48` | The edge-vertex cache does not share everything, and welding removes a class of sliver nobody expected it to. |
| `M-49` | ChunkLayout::cell_of inverts world_of_sample inside a cell and not reliably on its corner — M-32 in a second place. |
| `M-50` | E1 and M-34's ratio both reproduce live, under a mouse. |
| `M-51` | Marching Tetrahedra costs ~3× the triangles for ~4% worse geometry — and the literature's 2–3× is too low. |
| `M-52` | The Marching Tetrahedra ratio is 4.0 when the surface normal lies in one octant and 2.0 when it changes sign, and P-1's… |
| `M-53` | The five algorithms fill three of the four corners of manifold × intersection-free, and Marching Cubes is the only one i… |
| `M-54` | Dual Contouring is 101× more accurate than Marching Cubes on a sharp field, and indistinguishable on a smooth one. |
| `M-55` | O-14 falsified: Marching Tetrahedra's accuracy penalty is 4.3%, not 86%, and it beats Surface Nets rather than losing to… |
| `M-56` | Greedy meshing's 2.76× saving over face culling is a property of one scene, not of the algorithm: measured 1.70× to 256×… |
| `M-57` | Greedy merging manufactures T-junctions, and no weld can remove them. |
| `M-58` | A-010's vertex splitting removes the one-vertex-per-cell pinch completely, and the ticket named a field that was never a… |
| `M-59` | The dual of a manifold surface is a manifold complex, and an indexed triangle mesh cannot always represent it. This is a… |
| `M-60` | Only two of seven fields ever need a second vertex in a cell, and the rate falls with resolution — so Nielson's "about 1… |
| `M-61` | Splitting the vertex makes self-intersection worse, not better — ✗2's report is confirmed and the natural reading of M-2… |
| `M-62` | The t = a + b·n³ fit had been printing NaN since the day the algorithm names were spelled out, and once it prints number… |
| `M-63` | Both papers docs/research/ lists as "genuinely absent, blocking" are in the home-still corpus, so the acquisition lists… |
| `M-64` | A Transvoxel lateral face does not always cross the resolution boundary, and the exception is what transition cells are… |
| `M-65` | Central differences at the cell size cost under half a degree of normal direction, and converge at h². |
| `M-66` | On a sharp field the geometry and the field disagree by an angle that does not fall with resolution. |
| `M-67` | A sign test cannot distinguish 95.6% of the configurations a tetrahedron can actually be in. |
| `M-68` | parry3d's constructor is not a validity check: the only mesh it refuses is one with no triangles. |
| `M-69` | A chunk seam costs 72 boundary edges, and welding removes exactly those and nothing else. |
| `M-70` | Field-derived LOD is exact, not approximate: a coarse sample position is bit-identical to the fine one it sits on. |
| `M-71` | Cells fall by 8 per LOD level and triangles by 4 — and the 4 degrades exactly where the grid stops resolving the surface… |
| `M-72` | A sub-cell feature does not vanish under coarsening — it aliases, which is worse. |
| `M-73` | A transition cell that computes its sample positions by offsetting from a face origin puts a hairline crack in the seam,… |
| `M-74` | A zero-width transition cell stitches the hole and has no normal at all — which is what "severe shading problems" means,… |
| `M-75` | The transition width is what makes the patch's winding a measurable question at all, and the answer is unanimous. |
| `M-76` | Two blocks at differing resolution leave 88 unmatched boundary edges in the seam plane, and transition cells take it to… |
| `M-77` | Lengyel's Equation 4.2 loses its level index entirely when written in the block's own cells, and the seam holds at a rea… |
| `M-78` | A no_std crate cannot take a millisecond budget, and the honest API is a predicate the caller owns. |
| `M-79` | Subgrid MT's conformity is locality plus a shared vertex ordering — it is emphatically not invariance under relabelling,… |
| `M-80` | Every point §3.1 leaves unpaired is forced by arithmetic, not chosen — which is what makes the procedure implementable w… |
| `M-81` | §3.2.1's first two cases cover the whole of classic Marching Tetrahedra, and the rest of the machinery exists entirely f… |
| `M-82` | Theorems B.4 and B.6 hold against a reconstruction that does not use them — 27 of 27 (d₁, d₂) patterns, exactly. |
| `M-83` | T-002 is structurally blind to a fold inside a Steiner fan, and §3.2 is built entirely from Steiner fans — so "zero self… |
| `M-84` | The Figure-13 subdivision stencil closes: all four sub-tets it produces are themselves normal configurations, on all 27… |
| `M-85` | Property II is a statement about Γ_normal's residual, not about the tetrahedron's — and taking it for the tet's silently… |
| `M-86` | Conforming does not mean "no boundary at a shared face" — it means "boundary exactly where an open curve was discarded",… |
| `M-87` | §3.2 as implemented meshes 3,394 of 4,096 configurations, the entire remaining gap is one construction, and Property II… |
| `M-88` | §3.2.2's two labelling rules are stated separately and are checkable against each other, which turns a transcription int… |
| `M-89` | On thin_plate — the field this whole track exists for — 93.75% of tets already mesh, and every one of the remaining 6.25… |
| `M-90` | A scoop belongs to one face, but its two endpoints lie on an edge that belongs to two — so "both endpoints are on this f… |
| `M-91` | The contractible spanning disk introduces no vertices at all, and coverage moves 3,394 → 3,808 of 4,096. |
| `M-92` | §3.2 is complete over the tested space: 4,096 of 4,096 configurations mesh, with NoPattern and Inconsistent both zero —… |
| `M-93` | Subdivision's output reported 30 self-intersections in 52 triangles, and every one was an artefact of vertex duplication… |
| `M-94` | 1D root finding resolves a slab at 1/1000 of the edge length — and the fixture that appears to demonstrate its limit pas… |
| `M-95` | A-014c's acceptance is met with a number rather than a hedge: thin_plate returns 4,328 triangles at 33³ where greedy qua… |
| `M-96` | Orienting each triangle against the gradient at its own centroid is sufficient — the welded output is closed, manifold a… |
| `M-97` | §3.2.3's immersion is not hypothetical, and it is the single cause of every violation the subgrid extractor produces: 4… |
| `M-98` | Subgrid Marching Tetrahedra costs 70× classic Marching Tetrahedra and 196× Marching Cubes, and the ratio is field evalua… |
| `M-99` | The subgrid mesh's connectivity is provably manifold and my weld is what breaks it — the same mechanism as M-59, in a se… |
| `M-100` | A demo can be broken in a way that looks entirely correct: E-108's letters, centred on z = 0, are never lost by Marching… |
| `M-101` | §3.2.3's inset cannot be reconstructed from its prose, and two attempts made the measured result worse rather than bette… |
| `M-102` | ChunkId orders lexicographically on [x, y, z], so a residency sweep must iterate x slowest — and the natural z-outer loo… |
| `M-103` | rustdoc had never run on bevy_isomesh, and two doc links were already broken — the third time an excluded workspace turn… |
| `M-104` | A radius-based residency rule loads a ball of chunks, and for a heightfield that is 4x too many with most of them meshin… |
| `M-105` | E-203's first run reported 439 holes and a verdict of "G-001's overlap is wrong" — and the bug was in the test, one oper… |
| `M-106` | The acid test passes with a margin, and the margin is the interesting part: across 495 seam crossings the worst vertical… |
| `M-107` | λ is the sharpness/stability trade in one number, and swept over six decades it moves the runaway by a factor of 23 — bu… |
| `M-108` | The prediction was registered before the sweep and came out half right, and the half that failed is the more interesting… |
| `M-109` | M-61's "splitting the vertex makes self-intersection worse" is a 33³ fact, not a property — across nine resolutions it r… |
| `M-110` | The self-intersection counter never tests 99% of the pairs on a dual contouring mesh, and the skipped count is identical… |
| `M-111` | bevy_isomesh's four async drain tests do not pass on Linux, and one of them never passes — B-003's acceptance property i… |
| `M-112` | E-112's premise is wrong by an order of magnitude, and the failure it names is really two failures with two different la… |
| `M-113` | M-112's two laws reproduce in the committed example, and its fitted constant does not survive — the accuracy figure depe… |
| `M-114` | HermiteCell::from_corners is public and its contract is not: the corner order it requires lives in a private module, so… |
| `M-115` | A moving body is stopped harder and more often by ordinary terrain than by a chunk join — the same shape as M-106's answ… |
| `M-116` | Runtime convex decomposition costs 241–272 ms per fragment, which is 14–22 whole frames — so correct destruction collide… |
| `M-117` | A metric that reports a defect can be measuring the absence of a floor. E-204's first run accused 15 of 23 working colli… |
| `M-118` | B-003's async tests were a wall-clock race, not a plugin defect — and the number of iterations they ran was irrelevant t… |
| `M-119` | The obvious repair would have gutted the assertion it was repairing, and the mutation test is what proves the replacemen… |
| `M-120` | Transvoxel's mirrored seam works, and it had never been run — E-107 only ever meshed the fine block on the low-x side. |
| `M-121` | A level change moves the surface by up to 3.14 cells, which is the pop nobody measures. |
| `M-122` | Re-extracting a whole LOD ladder on every level change costs 12–23 ms and hitches; re-extracting only the blocks that ch… |
| `M-123` | Aliasing on a composed field came from the field, not the extractor — and swapping Surface Nets for Dual Contouring is w… |
| `M-124` | The amortized cost per frame is the number a game spends, and it tracks the budget to within one chunk across a 320x ran… |
| `M-125` | The never-livelock guarantee costs one chunk, exactly, and below one chunk's cost it degrades to precisely 1.00 chunks p… |
| `M-126` | BrushOp::commutes_with is sound and tight — on eight overlapping brushes it called all seven adjacent pairs correctly, w… |
| `M-127` | A fixture built to demonstrate that order matters reported 0 of 7 pairs where order mattered. |
| `M-128` | The dual methods do not tile across chunk boundaries and Marching Cubes does — measured on two fields, after the defect… |
| `M-129` | The seam counter's own exclusion list was missing an axis, and it accused Marching Cubes of a defect it does not have. |
| `M-130` | On a concave edge Dual Contouring's advantage is 3.6x in the mean and only 1.56x in the worst case — nothing like the 10… |
| `M-131` | The cell clamp does not bind on a concave edge either, which extends M-28 rather than limiting it — and the prediction t… |
| `M-132` | Subgrid Marching Tetrahedra does tile across a chunk boundary — 0 open edges in 20 configurations — and the 1 that said… |
| `M-133` | The dual methods are not reliably seam-closing rather than reliably open, which is a weaker and more useful claim than t… |
| `M-134` | M-21's negative intercept reproduces in direction on a second machine and a third of the range, and its magnitude is 40x… |
| `M-135` | The contour is 29% of a usable mesh, not 54% — and the largest stage is the collider check at 45%. |
| `M-136` | "What fraction is the contour" has no single answer — it is 13.1% to 74.3% across seven fields, and the variable is how… |
| `M-137` | Paint that lives in the field does not move when the geometry under it is destroyed — drift is exactly 0.000000, not sma… |
| `M-138` | The price of exact paint is that every sample walks the edit log, and it is sub-linear in practice: 2.33x the cost per c… |
| `M-139` | The wgpu pin holds, and the number it holds at is 29.0.4 rather than the 29.0.3 written down. |
| `M-140` | The GPU this repository's numbers will be measured on, and the two limits that bound a dispatch. |
| `M-141` | Include-once and cycle detection interact, and the obvious ordering silently accepts a circular graph. |
| `M-142` | GPU and CPU Marching Cubes agree on every triangle and disagree on 6% of vertices by exactly one ULP — and the cause is… |
| `M-143` | GridParams::sample_position was a mul_add where isomesh is origin + h·i, and a power-of-two cell size hid it through an… |
| `M-144` | GPU/CPU bit-identity is a property of the cell size, not of the port — 93.8% at h = 0.125, 1.1% at h = 0.1, 98.0% at h =… |
| `M-145` | The GPU extraction itself is essentially free and almost flat — 0.05 ms at 17³ and 0.13 ms at 129³, over 420× the cells… |
| `M-146` | Mesh shaders are advertised by this adapter and cannot be enabled from this workspace, because ExperimentalFeatures::ena… |
| `M-147` | Bevy's device already has mesh shaders enabled, with no configuration and no unsafe anywhere in this repository — so GPU… |
| `M-148` | A WGSL mesh shader's output block is a fixed contract that naga defines and no document states — and it is derivable fro… |
| `M-149` | A mesh-shader draw removes the smallest of the GPU path's three data-movement costs — 6.7% at 129³ — and the doc comment… |
| `M-150` | Moving the prefix sum onto the GPU took the extraction path from 15.01 ms to 9.65 ms at 129³ — 1.56× — and the stage it… |
| `M-151` | A hierarchical GPU scan is exactly the kind of code that is wrong invisibly, so it is asserted element-for-element again… |
| `M-152` | "The upload is field evaluation" is 57% true: of the 8.40 ms upload at 129³, evaluation is 2.65 ms, converting Vec<f32>… |
| `M-153` | Eliminating the "redundant" copy into GPU memory costs 1.6×, and the copy that was redundant was a different one. |
| `M-154` | GPU field evaluation does not match libm bit-for-bit on any field, including the ones with no transcendentals — and it d… |
| `M-155` | Evaluating the field on the GPU takes the 129³ path from 8.37 ms to 0.54 ms — 15.5× — and 37× ahead of a single-threaded… |
| `M-156` | The mesh-shader coverage hole is accepted rather than closed, and the deciding constraint is a stated goal rather than a… |
| `M-157` | A GPU-interpreted edit log reproduces BrushStack to under 8.4e-7 across every shape, every op and a twelve-brush mixed s… |
| `M-158` | The whole pipeline runs on the GPU and never touches the bus: field, edit log, extraction and draw, at 0.41 ms for 16,43… |
| `M-159` | The last four bytes cost 0.033 ms to move and 0.375 ms to wait for — because poll(Wait) drains every dispatch queued bef… |
| `M-160` | A zero-synchronisation extraction makes CPU time independent of grid size: flat at ~0.17 ms from 33³ to 129³, where the… |
| `M-161` | §3.2.3's immersion is real, measured on five of seven fields — and it is not the defect. |
| `M-162` | A-014d's blocking question is answered, and the answer is no: the inset needs neighbour information, and the neighbours… |
| `M-163` | csg_difference's three surviving non-manifold edges are each 4 faces but only 3 distinct polygons, from exactly 2 tetrah… |
| `M-164` | gyroid's 138 inconsistently-oriented edges are not A-014d's, and the ticket's stated mechanism reaches 8% of them. |
| `M-165` | Among the configurations that can exhibit it, Chernyaev's numerator-only test is wrong 12.6% of the time — 1,966 of 15,6… |
| `M-166` | The interior test does not inherit ambiguity's crack-freeness argument, and the reason is one operator. |
| `M-167` | Across the whole GPU series the arithmetic never moved and was never the point: synchronisation was 83% of an extraction… |
| `M-168` | Giving a crossing an identity instead of a position removes 5.01× of the subgrid extractor's vertices and changes no tri… |
| `M-169` | Identity-based sharing is complete exactly when no root lands on a grid sample point, and the correlation is not approxi… |
| `M-170` | A GPU adapter is a finite resource and the test harness was treating it as free: 67 tests opened 67 devices, and once an… |
| `M-171` | A shipped example panicked on every run, and three siblings carried the same race unfired. |
| `M-172` | A signed distance field's gradient is exactly zero on its medial axis, and for a slab the medial plane is precisely wher… |
| `M-173` | Two coincident surfaces are not a rendering bug to be biased away, they are a modelling statement — and the fix is to st… |
| `M-174` | CI has been red on every push of the GPU series, and behind the one accepted failure three unrelated ones accumulated un… |
| `M-175` | Sorting by magnitude is not enough to make a sum permutation-invariant, and the comment that said it was had been wrong… |
| `M-176` | Zero-padding a reduction is transparent, negative zero included — and the reason is the accumulator's seed, not the sort… |
| `M-177` | Reordering cannot buy negation equivariance, and the obstruction is structural rather than a missing tie-break. |
| `M-178` | A-016 moved 34 of 168 golden rows, not the 42 predicted, and the 8 that held identify their own mechanism. |
| `M-179` | A-014h's stated mechanism finds almost nothing, and the reason is a deliberate choice in the root finder. |
| `M-180` | The population M-169 named as one defect is two, and only the smaller one is an identity problem — so A-014h cannot reac… |
| `M-181` | T-009's own premise was wrong: normalising four scattered weld epsilons onto one moved nothing, because the weld's answe… |
| `M-182` | P-7 confirmed, and the check that confirmed it found thin_plate failing at two of the three resolutions nobody was testi… |
| `M-184` | Naming a root by the grid point it lies on completes the identity outright: on every field the extractor now emits each… |
| `M-185` | Completing the identity turns a sliver into a repeated-index triangle, and the extractor now declines to emit those. |
| `M-186` | M-162 is falsified, and A-014d's blocker with it: after A-014h no coincident polygon has a foreign tetrahedron on its bo… |
| `M-187` | Orienting each connected component from its most confident triangle drives the inconsistently-oriented-edge count to exa… |
| `M-188` | Two of the seven reference fields are non-manifold at a resolution the suite never tested, and one of them is torus. |
| `M-189` | Running the subgrid validity census at the three resolutions its own gate names doubles the number of known defects, and… |
| `M-192` | I-007's own finding happened to I-007, in the commit that filed it. |
| `M-193` | The polygon type A-014d has to inset is the one type its figure could not be read for — measured, not assumed. |
| `M-194` | A-014d's acceptance criterion asks §3.2.3 to fix a defect §3.2.3 does not cause, and this repository's own instrument ha… |
| `M-195` | The mirrored seam's crack is exactly invisible to a seam-plane counter — 0 edges on the plane, 28 on each side of it. |
| `M-196` | Consuming fields:: instead of a hand-rolled copy is a speed-up, not a tidy-up — and a dispatch layer can throw the whole… |
| `M-197` | GPU-013's two halves both misnamed their target, and the measurement says the prize is waits rather than submissions. |
| `M-198` | The publish job ran for the first time in its life, and failed on a secret that was never set — so the release pipeline… |
| `M-199` | Subgrid Marching Tetrahedra is already intersection-free on every reference field at every resolution — so §3.2.3, whose… |
| `M-200` | Two of A-014i's three recorded defects are not what the review called them: one is real and its obvious fix is measurabl… |
| `M-201` | The subdivision path has a grid-level fixture at last, and it says the orphaned vertices are real — 9.8% of them — which… |
| `M-202` | Subdivision only fires where the grid badly under-resolves the field, so "T-001 validity on the subdivision fixture" can… |
| `M-203` | A-014i's orphaned vertices are fixed by giving children their parent's names instead of new ones, and the mesh is bit-fo… |
| `M-204` | M-166 is closed by two parentheses, and the test that was already asserting the property could not have caught it. |
| `M-205` | A-014i's "orphaned bigons" is real, and it is a law rather than a defect: unreferenced positions arrive in pairs, one pa… |
| `M-206` | Two independently derived constructions locate the same body saddles, to 1.1e-12. |
| `M-207` | The reference implementation loses a root this one keeps, in two places, and both are the textbook quadratic formula. |
| `M-208` | Interior ambiguity was unreachable by this crate's entire test suite: 0 of 68,385 reference-field surface cells have six… |
| `M-209` | The eighth reference field is searched, not designed — 97 of 610 candidates qualify, and the one taken is the gentlest r… |
| `M-210` | The natural level set of gradient noise is its own lattice, which is degenerate — so a noise reference field may not use… |
| `M-211` | A-010's zero was conditional, and nobody knew on what: Manifold Dual Contouring is not manifold on a field with interior… |
| `M-212` | Welding two coincident vertices can create a non-manifold edge, which a test comment had asserted was impossible. |
| `M-213` | Orientation can raise the flipped-edge count, and after <= before was not a law. |
| `M-214` | The tunnel and the twelve-vertex contour are told apart by counting rings, and both are reachable: 2,053 and 173 in 396,… |
| `M-215` | The u pair's two lines are crossed relative to the other two, and the inner hexagon is what proves it. |
| `M-216` | The interior vertex is a transcription, and the check on it is geometric: 149,803 of them, every one on the level set to… |
| `M-217` | The disk path costs no new budget: worst case 12 triangles and 1 interior vertex per cell, against MAX_TRIANGLES = 12 an… |
| `M-218` | The twelve-vertex contour needs a closing step the tunnel does not, and a manifoldness test found it rather than a readi… |
| `M-219` | The reference implementation has a transcription-grade typo in the detached-ring test, and it is one line from six. |
| `M-220` | The singular face is an artifact of quantised data, and this crate cannot reach it: 0 of 1,838 on eight reference fields… |
| `M-221` | 0 × NaN is NaN, and it took the extractor to find that — not 400,000 random cells. |
| `M-222` | χ falls by exactly two per tunnel and by nothing else — the interior rule's topology change is arithmetic, not approxima… |
| `M-223` | The interior rule costs 1.95% at 33³ and 0.14% at 65³, on the only field that exercises it. |
| `M-224` | Manifold Dual Contouring's non-manifoldness has nothing to do with tunnels, which is what A-017 assumed, and it survives… |
| `M-225` | A-017's mechanism, and the grid predicts the mesh defect exactly — 30, 64, 8 and 40, with zero error. |
| `M-226` | Subgrid Marching Tetrahedra's output does not need welding and is damaged by it — the doc had been telling consumers to… |
| `M-227` | Orientation now reaches zero on every reference field at every resolution, and M-187's caveat was about the walk rather… |
| `M-228` | Grosso's tunnel triangulation has an undefined case, it is reachable, and 400,000 random cells could not reach it. |
| `M-229` | The contour-count discriminator misclassifies case 13, and the misclassified cells are exactly the ones with no triangul… |
| `M-230` | Corollary 6 is the tunnel test; Proposition 1 is not, and the derivation that shows it also validates itself. |
| `M-231` | The [9,3] cell is not a topological subcase; it is a singular face that the strict interior test lets through. |
| `M-232` | A singular face is unreachable from continuous data and routine from quantised data — and the rate at u8 density matches… |
| `M-233` | A-002i's blocker is not the vertex cache, it is that a singular face needs a third routing and the resolution mask has t… |
| `M-234` | A-017 closed by decision rather than by code, and the limit is stated as an identity so it cannot rot. |
| `M-235` | The example capture rig renders headlessly but cannot control its own window size, and it fails silently — which is the… |
| `M-236` | The shootout's header claimed the wrong counts for both of the things it was counting, and nothing could have caught it. |
| `M-237` | The QEF buys 2× accuracy on smooth fields and 100× on sharp ones, and pays for all of it in self-intersections — measure… |
| `M-238` | Probabilistic quadrics are this crate's existing solve with a different regularizer, and the derivation is the finding —… |
| `M-239` | A window can resize itself without a window manager, and two obvious fixes both fail silently before the third works. |
| `M-240` | The crate could subtract and intersect but not union, and nothing noticed because no reference field unions anything. |
| `M-241` | Two demo GIFs were shipped that did not show what their captions claimed, and a single frame is not an inspection. |
| `M-242` | The shootout's 112 rows came back structurally identical after X-001 rewrote how it enumerates, which is the refactor's… |
| `M-243` | An unquoted heredoc executed the words in this script's own Python comments, and it failed silently for a value of "fail… |
| `M-244` | A declared Lipschitz constant was wrong on the first run, and the test that caught it was written in the same commit by… |
| `M-245` | The eikonal condition cannot tell a CSG underestimate from a true distance, which is why a bound needs two numbers rathe… |
| `M-246` | min is not better than max; each is exact in one region and wrong in the other, by the same amount. |
| `M-247` | Repeated CSG destroys the worst case and leaves the typical case untouched, and a renderer barely notices what a precisi… |
| `M-248` | One field evaluation replaces 576, and buys between 1.1× and 11.8× depending on how much of the volume the surface reach… |
| `M-249` | A directional Lipschitz bound buys exactly nothing on five of six fields, and 1.80× on the sixth — the null result F-006… |
| `M-250` | Refining an edge crossing on the real field helps curved fields by 13–15% and does nothing at all for the CSG one — the… |
| `M-251` | The exact distance transform agrees with brute force to the last bit, and is exactly one sample spacing from the analyti… |
| `M-252` | Fast sweeping beats the exact transform everywhere, including where it was predicted to lose — because the seeding, not… |
| `M-255` | naive reinitialisation moves the zero set, measured before it was fixed (S-004) |
| `M-256` | the narrow band's cost claim needs the march, not the sweep (S-004) |
| `M-257` | the approximate GPU method beats both exact CPU ones (S-005) |
| `M-258` | a u32 followed by vec3<u32> is 32 bytes, not 16 (S-005) |
| `M-259` | the round trip the crate did not have (S-006) |
| `M-260` | a uniform grid over the sample cells lost to a flat box reject, 3.9× (S-006) |
| `M-261` | Real gained acos, and libm is why (S-006) |
| `M-262` | the winding number beats the pseudonormal on holed meshes, by a widening margin (S-007) |
| `M-263` | the boundary must be counted with multiplicity, not as a boolean (S-007) |
| `M-264` | the uncertified set is a curve, not a resolution failure (T-015) |
| `M-265` | decimation is re-sampling on a nested grid, so the literature's rule has no bite here (T-016) |
| `M-266` | M-72's aliasing is alignment, not chance (T-016) |
| `M-267` | the sampled gradient supremum is not even monotone in sampling density (T-017) |
| `M-268` | field quality is now in the regression gate, and every column is exact (T-017) |
| `M-269` | a grid-aligned ray double-counts shared edges, and the constructor shootout is what found it (T-018) |
| `M-270` | a benchmark that hands a repair algorithm perfect input measures the benchmark (T-018) |
| `M-271` | the exact transform is the hungriest constructor, and the mesh-based pair are free (T-019) |
| `M-272` | the pre-registration gate is a const assertion, and that is the whole of it (R-000) |
| `M-273` | the first thing done with the pre-registration mechanism was amend a registration to fit the code (R-002) |
| `M-274` | the fixture never contained the configuration both experiments were about (R-002) |
| `M-275` | FALSIFIED on its first clause, and the residue is an edge defect (R-003) |
| `M-276` | the dual methods' non-manifold edges are the ambiguous face, all 314 of them (A-021) |
| `M-277` | the index generator was blind to the shape its own file had moved to (R-004) |
| `M-278` | HELD, and the crack budget is not where the ticket assumed (R-004) |
| `M-279` | the mechanism is FALSIFIED, and the registered falsifier could not have caught it (R-005) |
| `M-280` | this harness's nanoseconds are not a unit, and the committed Zen 3 sweep is 1.45× stale (R-005) |
| `M-281` | a timing here is a property of the binary, not only of the code (M-001) |
| `M-282` | the whole family, in one binary and one run (M-001) |
| `M-283` | FALSIFIED, and the fixture that confirmed it agreed to four decimal places (R-006) |
| `M-284` | HELD, and the dual's IPC wall is one function (R-007) |
| `M-285` | the dual's axis had to be a constant, and that was 82% of its cycles (A-023) |
| `M-286` | the misses M-279 measured as free were hidden behind the stall, and now they cost (A-023) |
| `M-287` | one bit of the row length was a 3.4× tax at the chunk size everybody uses (A-024) |
| `M-288` | FALSIFIED, and the registered definition did not implement the registered claim (R-008) |
| `M-289` | the reference gradient was noise, and it falsified two hypotheses that were true (R-009) |
| `M-290` | the ambiguous face in the dual path, with the source read at last (A-022) |
| `V-1` | wgpu / wgpu-types / naga 29.0.3, glam 0.32.0, encase 0.12 |
| `V-2` | Bevy 0.19 removed RenderGraph; passes are systems in ECS schedules; non-camera work targets the RenderGraph schedule |
| `V-3` | Marching Cubes peak: 5.42 G voxel/s, 330 M tri/s (RTX 2080 Ti). DMC costs 1.52–3.50×; FlexiCubes 2.77–3.92× |
| `V-4` | Contouring 68 ms vs halfedge construction 58 ms |
| `V-5` | On unstructured grids, Delaunay/MT ratio 15.3×–81.5× — contouring is 1–2% of the pipeline |
| `V-6` | 73% of FlexiCubes' 64³ Marching Cubes timing is fixed launch overhead |
| `V-7` | Cross-paper reproducibility floor is ~1.5× in opposite directions: TetWeave re-measured FlexiCubes at 128³ as 9.63/15.25… |
| `V-8` | 10.7× more bandwidth bought ~1.7× more throughput |
| `V-9` | Same Marching Cubes, compute shader → mesh shader: 114.2 → 2679.4 fps (23.4×) |
| `V-10` | CBT sum-reduction, atomics → LDS staging: 5.78 → 0.40 ms |
| `V-11` | Meshlet compression: 15.5 M tri in 0.59 ms (RX 7900 XTX) |
| `V-12` | Work graphs: 79,710 instances in 3.74 ms — but 2.8–3.4× slower on classification workloads |
| `V-13` | nvblox: meshing is the least GPU-accelerable stage, ×3–13 vs fusion's ×174–177 |
| `V-14` | Aokana renders 10¹⁰ voxels at 6 ms, 5% resident, RTX 3060 Ti — explicitly not editable |
| `V-15` | CoACD vs V-HACD: 49% → 80% downstream manipulation success |
| `V-16` | Dimforge migrated parry (0.26.0) and rapier (0.32) off nalgebra onto glam, citing rust-gpu support; performance "nothing… |
| `V-17` | No paper since 2020 benchmarks Marching Cubes vs Surface Nets vs Dual Contouring against each other. |
| `V-18` | Dual Contouring's own paper quantifies the f32 QEF failure. |
| `V-19` | Dual Contouring's topology is Surface Nets' topology. |
| `V-20` | A QEF is stored as AᵀA (symmetric 3×3), Aᵀb (3-vector) and bᵀb (scalar) — 10 floats — rather than as A and b |
| `V-21` | The corner/contractible spanning disk is supposed to come out degenerate, and mistaking that for a bug would send A-014b… |
| `V-22` | (settled at E-206b: adopted as a bevy_isomesh dev-dependency only, avian3d 0.7.0, which compiles clean against Bevy 0.19… |
| `V-23` | CLAUDE.md's Metal contradiction is resolved, and both sources were right about different layers. |
| `V-24` | Custodio §5.1's correction, verified from the paper and re-derived rather than transcribed. |
| `V-25` | "Lewiner's reference implementation omits disambiguation for cases 10 and 12 entirely" is true of the code and false of… |
| `V-26` | A follow-up exists that builds MC33's triangulation with no lookup table at all, and A-002b's blocker did not know about… |
| `V-27` | §3.2.3's prose is fully retrievable and its triangulation patterns are still only a picture — and the paper offers a sec… |
| `V-28` | A-014d's rule-5 stop is lifted, and not by the figure — by the authors' own reference implementation, shipped as arXiv a… |
| `V-29` | Two sources this repo recorded as unobtainable were obtainable, and one of them removes the largest single piece of A-00… |
| `V-30` | Grosso's quadratic re-derived, and now agreeing three ways rather than two. |
| `V-31` | The reference implementation both Grosso papers cite has been deleted from GitHub, and it survives — the paper's own lis… |
| `V-32` | Two more rows of meshing-library-target.md marked PAYWALL are in the corpus, and the table's status code has now been wr… |
| `V-33` | A paper_download that reports success can return a landing page, and this is the third producer of that signature. |
| `V-34` | Manifold Dual Contouring's uniform-grid criterion is one vertex per cycle of a decider-modified Marching Cubes table, an… |
| `O-1` | Settled at G-002 (M-33, M-34), and confirmed live under a mouse at E-202 (M-50). |
| `O-2` | Settled at A-009 (M-28, M-29): not entirely, and the residue names its own mechanism. |
| `O-3` | Marching Cubes vs Surface Nets vs Dual Contouring vs MT — actual relative speed on one machine? |
| `O-4` | Settled at G-003 (M-36, M-37, M-38): conditionally, and the condition is narrow. |
| `O-5` | Do mesh shaders work on macOS/Metal? |
| `O-6` | What is amortized meshing cost per frame under continuous editing? |
| `O-7` | What fraction of our pipeline is contouring vs everything else? |
| `O-8` | Does Dual Contouring's vertex placement need f64 in practice, or is f32 enough? |
| `O-9` | How much does T-003's gradient-flow chord over-estimate distance at a concave seam? |
| `O-10` | What is Surface Nets' non-manifold rate as a function of feature thickness over h? First curve measured at A-010 (M-60),… |
| `O-11` | (ANSWERED at R-007, M-284: the dual carries a fourth stage Marching Cubes has not — emit_quads, which walks every grid e… |
| `O-12` | Is Marching Cubes unconditionally manifold now? |
| `O-13` | Pre-registered: Marching Tetrahedra vertex count = 3.0× Marching Cubes, converging from above Confirmed at A-003/M-001,… |
| `O-14` | Pre-registered: Marching Tetrahedra symmetric Hausdorff at 64³ ≈ 2.6e-3, about 1.86× Marching Cubes, i.e. slightly worse… |
| `O-15` | Answered at A-003 (M-52): the normal's sign pattern, not its direction. |
| `O-16` | Can the parallel dual-edge collapse (M-59) be removed without giving up the cycle partition? |
| `E×1` | Surface Nets' centroid as Dual Contouring's vertex rule |
| `E×2` | A separate probabilistic-quadric solver |
| `E×3` | Crossing-count-scaled regularizer |
| `E×4` | Weld gated on the pairwise link condition, rejected pairs left split |

<!-- END GENERATED INDEX -->

## When this file splits

**Stated in advance, so the decision is not made under pressure by whoever is holding it when it
becomes unbearable.** At **500 entries or 600 KB**, whichever comes first, `FINDINGS.md` becomes the
index plus Parts 5 and 6 — the method rules and the how-to — and Parts 1 through 4 move to
`findings/` split by *axis*, not by date: `findings/correctness.md`, `findings/performance.md`,
`findings/sources.md`, `findings/method.md`. Date order is how the file grew and is useless for
lookup; nobody has ever wanted "everything we learned in August".

Two things do not change at the split. **The index stays generated** — `scripts/findings_index.sh`
gains the directory and keeps one table across all files, because an index per file is four places
to look and therefore none. And **entry numbering stays global**: `M-231` keeps its identifier
wherever it lives, since every cross-reference in the backlog, the archive and a year of commit
messages names it by number and a renumber would silently break all of them.

The index above carries the current count; it is generated, so it cannot drift from the file.

## Confidence tiers

| Tier | Meaning | Bar |
|---|---|---|
| **M** | **Measured here** | We ran it. Code and numbers in this repo, reproducible by checkout. |
| **V** | **Verified externally** | We read the primary source ourselves. DOI or file attached. |
| **R** | **Reported** | A credible source asserts it; we have not independently checked. |
| **F** | **Folklore** | Widely repeated, no verified primary source found. |
| **✗** | **Falsified** | Tested and found false. Includes what we believed and why. |

**Never cite an R or F claim as justification for a design decision without saying which tier it is.**

---

## Part 1 — Falsified

The most valuable section. Each row includes where the wrong belief came from, because provenance
predicts the next error.

### ✗1 — "Surface Nets produces substantially fewer triangles than Marching Cubes"

**Believed because:** stated in this repo's own implementation brief, and near-universal folklore.
**Falsified by:** A-004 measurement, then derivation, then independent numeric check.
**What's true instead:** on a closed manifold surface meshed on the same grid, both counts are pinned
by Euler's formula:

```
V_sn = V_mc + χ            F_sn = F_mc + 2χ
```

Exactly, at every resolution. `V_mc = C` (one vertex per crossed grid edge), `F_sn = 2C` (two
triangles per crossed edge), `F = 2V − 2χ` on any closed triangulated manifold. The middle step is a
combinatorial identity worth naming on its own: **surface cells = crossed edges + χ.**

| field | χ | crossed edges C | surface cells S | S−C |
|---|---|---|---|---|
| sphere 17³ / 49³ | 2 | 414 / 4206 | 416 / 4208 | 2 / 2 |
| torus 33³ / 49³ | 0 | 1296 / 3216 | 1296 / 3216 | 0 / 0 |
| box 49³ | 2 | 4374 | 4376 | 2 |
| two disjoint spheres 49³ | 4 | 1300 | 1304 | 4 |

**Consequence:** Surface Nets' case must rest on quad connectivity and inner-loop cost, not output size. In
M-001 the count columns are a **checksum with a predicted value**, not a result.
**Would be shown wrong by:** any closed-manifold pair on the same grid where the difference ≠ 2χ.
**Legitimately breaks at:** boundary-clipped meshes (**incoming at G-001 chunking** — expect the
assertion to fail there and do not "fix" it), A-013 welding, Marching Cubes vs Marching Cubes 33 differing in χ.
**A-013 landed and sharpened the welding clause (M-48):** a weld lowers `V_mc` only where a grid
*sample* lands on the isosurface, which is a sampling accident rather than a property of welding —
measured at 48 of 654 vertices on `sphere` at 25³, and **zero** at 17³ and 33³ and on four of the
seven fields entirely. So the identity is not broken by A-013 existing; it is broken by a field that
happens to touch a lattice point, and by chunking.

### ✗2 — "You can have a manifold mesh or an intersection-free one, not both"

**Believed because:** folklore, repeated in several secondary sources.
**Falsified by:** literature review round 1. Manson & Schaefer 2010 achieved both. ODC (2024) measured
Manifold Dual Contouring at **100% of models self-intersecting** against ODC at **0 of 1500**.
**Consequence:** guaranteed intersection-free extraction is on the table, which is the premise under
A-009 and the runtime-convex-decomposition opportunity.

### ✗3 — "Every interior Surface Nets vertex has four neighbours"

**Believed because:** written into isomesh's own module docs before measuring.
**Falsified by:** A-004. Measured max degree **10** — higher than Marching Cubes' **9**.
**Method rule earned:** a doc comment the test suite disproves is worse than no doc comment.

### ✗4 — "Dual Contouring is absent from the home-still corpus"

**Believed because:** `distill_search` returned nothing.
**Falsified by:** `catalog_read` on `10.1145_566570.566586` — present, with zero Qdrant chunks.
**Root cause:** **342 documents are readable but invisible to `distill_search`.**
**Method rule earned:** presence in the corpus is decided by `catalog_read` / `catalog_list`, **never**
by search. Any claim of absence made from a search result is unfounded.

### ✗5 — "naga_oil is the shader composition path for Bevy"

**Believed because:** it was, and it's the name everyone knows. Carried into this project's premise.
**Falsified by:** naga_oil's last release is **v0.18, 2025-06-26** (14 months stale), and Bevy 0.19's
own release notes say they will *"port our existing internal shaders to use WESL, and endorse it as
the shader language of choice for Bevy."*
**Consequence:** GPU-002 uses a ~40-line preprocessor; WESL revisited when shader count justifies it.

### ✗6 — "Mesh shaders aren't reachable from inside Bevy"

**Believed because:** an imprecise claim made in this project's own speed analysis.
**Falsified by:** wgpu **v28** shipped `Features::EXPERIMENTAL_MESH_SHADER`; Bevy 0.19 pins wgpu
29.0.3; `WgpuSettings.features` requests it at device creation and
`RenderDevice::wgpu_device() -> &wgpu::Device` gives raw access.
**What's actually true:** wgpu's *maturity* is the blocker, not Bevy. Experimental, redesign issue
open, no browser support, **Metal status contradictory** (see O-5).

### ✗7 — "CBT 5.78 → 0.40 ms is from Dupuy 2020"

**Numbers confirmed. Attribution wrong** — it's the Unity SIGGRAPH 2021 talk.
**Method rule earned:** verify the *source* separately from the *number*. A correct figure with a
wrong citation propagates as an uncheckable claim.

### ✗8 — Velo3D "93M vertices / 31M faces"

**Unsupported.** No primary source located. Removed from the catalog rather than softened.

### ✗9 — "Marching Cubes' cost inside the volumetric loop was never measured" / "navmesh rebuild cost was never measured"

Both false. Dong 2018 measured meshing at **76.5–89.6%** of the pipeline; van Toll 2012 measured
navmesh rebuild. Both were asserted as gaps in this project's v2 catalog and corrected in round 1.
**Method rule earned:** "nobody measured X" is a claim requiring the same evidence as any other.

### ✗10 — "glam should be the internal math library from day one"

**Believed because:** stated in this repo's architecture doc.
**Refined by:** the agent at I-001. The core took **libm only**; glam is deferred to A-007, where the
3×3 solve is the first thing that actually needs matrix math.
**Why the refinement is better:** the public API is arrays, so glam was never load-bearing; deferring
it means the crate carries zero exposure to glam's ~quarterly breaking releases until it buys
something. Consumers need no conversion impls either — `glam::Vec3::from([f32; 3])` already exists.

### ✗11 — "Plain Marching Cubes has ambiguous faces and produces holes"

**Believed because:** stated in this repo's own implementation brief (Stage 2, "Plain Marching Cubes has ambiguous
faces and produces holes"), carried into `BACKLOG.md`'s A-002 acceptance criterion, and near-universal
folklore about Marching Cubes.
**Tested by:** `validate_table()` (`marching_cubes/mod.rs:319`), which checks all 256 cases structurally, and the
assertion `assert_eq!(report.face_disagreements, 0)` at `marching_cubes/tests.rs:30`.
**Result:** zero face disagreements, across all 256 cases.
**What's true instead:** holes require two cells sharing a face to *disagree* about how the surface
crosses it. In this implementation a face's segments are a function of that face's own four corner
signs and nothing else — the two cells meeting on a face read the same four corners, so they cannot
disagree. The property is structural, not empirical, and it falls out of the table being **derived at
compile time by walking each face counter-clockwise** rather than transcribed from a diagram.

The folklore is not wrong about Marching Cubes in general; it is wrong about *this* Marching Cubes.
Lorensen & Cline's original table was transcribed per-case and its ambiguous cases were resolved
inconsistently between complementary configurations, which is where the holes came from.

**Consequence:** A-002's acceptance criterion was unsatisfiable and has been re-scoped. Marching Cubes 33's
remaining value is topological agreement with the trilinear interpolant — a genuinely different
surface on ambiguous faces, measurable as a **χ difference** — not crack-fixing. The research is
explicit that a game wants *consistency* over topological fidelity, so the `L` slot is now spent
knowingly rather than against a test that could never go green.
**Would be shown wrong by:** any field producing `boundary_edges > 0` from A-001 on a closed field, or
any non-zero `face_disagreements`.

### ✗12 — "The equivariant vertex rule needs a fast three-plane path with a fallback"

**Believed because:** `BACKLOG.md` split it into A-007 ("the three-plane rotation-equivariant rule…
falls through to A-008 when the triple product is near zero") and A-008 ("for >3 planes and degenerate
cells"), and the brief says "falls back when the triple product is near zero (near-parallel planes)".
The crate architecture doc says the same in a third phrasing: "fall back to the regularized
normal-equation form only for >3 planes".
**Tested by:** reading the audit doc all three of them cite —
`docs/research/2026-08-10-adjacent-math-transfer-audit.md:182-219`.
**Result:** the audit gives the Tikhonov-adjugate form as the *production* form and describes it as
"branch-free, handles all degeneracies", closing with "**no eigendecomposition, no SVD, no iteration,
no data-dependent branch**". It is a single unconditional path in the source, not a fallback arm.

Worse, the audit's diagnosis of *why Dual Contouring pops* is the branch itself: Dual Contouring's hard SVD
truncation at σ < 0.1 is a discontinuous branch, and over 20,000 trials seeded at the threshold in f32
the rank branch disagreed after a rotation in **454 cases**, with `‖f(Rx) − Rf(x)‖` median **2.13** and
max **9.10** — a several-cell vertex pop from an infinitesimal rotation. A triple-product threshold is
the same construction with a different discriminant, so the split would have reintroduced the exact
failure the rule exists to remove.

The measured equivariance residual (f32, coordinates in [0,256], 4000 random cells) also shows the
"fast path" is not the accurate one:

| rule | median | p99 | max |
|---|---:|---:|---:|
| Dual Contouring normal equations | 6.80e−05 | 2.48e−01 | 5.6e+02 |
| dual basis (Cramer) | 1.61e−05 | 7.23e−04 | 3.6e−01 |
| **Tikhonov adjugate** | **1.59e−05** | **1.81e−04** | **6.4e−04** |

Tikhonov dominates Cramer on both tail columns, so nothing is traded away by dropping the three-plane
form. The two paths also do not agree to within noise, which means the branch would have been
*observable* in the output.

**Consequence:** A-007 and A-008 merged into one ticket with one unconditional path. Two requirements
the audit states and no ticket had recorded are now in it: **magnitude-sorted 3-term dot products**
(4328/9600 equivariance failures unsorted, **0/9600** sorted — the guarantee does not hold in f32
without this), and the derivation of **λ = 0.01** as the value that reproduces Dual Contouring's σ = 0.1 truncation
smoothly. The corpus circulates three constants — 0.01, 0.1, and σ=0.1 — and an implementer reading
only the algorithm catalog would have picked 0.1.
**Would be shown wrong by:** a measured configuration where the adjugate form is less accurate or less
equivariant than the Cramer form, or where `det(M + λI)` is small enough at λ = 0.01 to matter.

### ✗13 — "`Real::as_f32` is the only narrowing operation in the crate, and the crate itself never calls it"

**Believed because:** stated in `real.rs`'s own doc comment for `as_f32`, which frames it as an output
convenience — "it exists for consumers writing into an `f32` vertex buffer".
**Tested by:** `grep -rn "as_f32()" crates/isomesh/src/`.
**Result:** **three** call sites, none of them a consumer and none of them writing a vertex buffer:

| site | purpose |
|---|---|
| `fields/noise.rs:82` `lattice_index` | Perlin lattice coordinate |
| `validate.rs:810` `quantise` | duplicate-vertex bucketing |
| `validate/tri_grid.rs:60` `cell_of` | spatial-grid binning |

All three are the same idiom, `f.as_f32() as iN` — a float→integer narrowing. Two predate this
session, so the doc comment has been false since I-002.

**Why this is more than a tidiness point.** The doc's framing hides the property that actually matters
at these call sites: `as_f32` on an `f64` is exact only up to `2²⁴`, so using it as an integer
narrowing step carries a silent lattice cliff. `noise.rs` knows this and guards with a
`LATTICE_LIMIT` debug assertion; the other two did not, and nothing said they should. See M-18 for
where that bites.
**Consequence:** the doc comment corrected in place, and T-008 opened for the one exposed call site.
**Method rule earned:** ✗3's, again — a doc comment the code disproves is worse than no doc comment.
Twice now the false comment has been an *architectural* claim ("nothing calls this", "every vertex has
four neighbours") rather than a factual slip, which is the kind that shapes later decisions.

### ✗16 — "glam 0.32 lands with A-007's vertex solve"

**Believed because:** stated in four places — `CLAUDE.md`'s crate layout ("Deps: libm today; glam 0.32
joins it at A-007. Nothing else, ever.") and its dependency-justification section, `BACKLOG.md`'s A-007
ticket, `BACKLOG_ARCHIVE.md`'s I-001 note, and ✗10 itself, which deferred glam to A-007 on the grounds
that "the 3×3 solve is the first thing that actually needs matrix math".
**Tested by:** reading glam 0.32.1's source before adding the dependency.
**Result:** **glam has no scalar abstraction.** `Mat3` lives in `src/f32/mat3.rs` and `DMat3` in
`src/f64/dmat3.rs` as separate concrete types, `lib.rs` re-exports per-scalar modules, and there is no
generic `Mat3<T>` and no trait spanning the two. The only `pub trait` in the crate is `FloatExt`, a
scalar extension trait.

**What's true instead:** the premise was right and the conclusion did not follow. The solve *does* need
matrix math — but this crate is generic over `Real`, which spans `f32` **and** `f64`, and glam's types
do not. Using it would mean a bridge trait with two impls forwarding every operation, which is more
code than the 3×3 adjugate it would wrap, adds a dependency, and puts two float backends inside one
solve — the exact thing the `libm` justification rejects.

So the 3×3 lives in `dual_contouring/solve.rs` as a six-entry symmetric matrix over `[R; 3]`, about 40 lines, and
**the crate stays at one dependency.** The "as light as possible" pitch survives A-007 intact.
**Would be shown wrong by:** glam gaining a generic scalar parameter, or this crate dropping `f64`.
**Note this is ✗10's second correction.** ✗10 moved glam from "day one" to "A-007"; the deferral target
was wrong too. The recurring error is reasoning about glam from its reputation rather than its API.

### ✗15 — "Marching Cubes is unconditionally manifold"

**Believed because:** every measurement in this repo said so. M-4 contrasts Surface Nets'
non-manifoldness against *"Marching Cubes' zero at every resolution"*, the README says *"Marching Cubes
stays manifold"*, and the mechanism looked airtight — Marching Cubes places vertices on grid **edges** rather than
one per cell, so the multi-sheet argument that sinks Surface Nets does not apply. `SurfaceGate::Closed`'s own doc
comment asserted it.
**Falsified by:** T-005b's `marching_cubes_meshes_sphere_unions`, on its first run against a fresh
proptest seed — during T-006, which is a nice demonstration that the property tests keep working after
the ticket that wrote them.
**Result:** a union of three spheres at `h = 2/3` gives **2 non-manifold edges and 3 non-manifold
vertices** on a mesh that is otherwise perfect: closed, `χ = 2`, one component, consistently oriented,
zero boundary edges.

**What's true instead:** Marching Cubes is manifold when the grid **resolves** the surface. Where the surface
*pinches* inside a single cell — two lobes of a union meeting at sub-cell scale — the shared grid edge
ends up carrying four faces. Refinement fixes it, and sharply:

| n | 7 | 9 | 13 | 17 | 25 | 33 | 49 | 65 |
|---|---|---|---|---|---|---|---|---|
| non-manifold edges | **2** | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

**Consequence:** the property suite's gate is renamed from `ClosedAllowingMultiSheet` to
**`ClosedAllowingUnresolvedTopology`** and now covers Marching Cubes on generated fields too, because the condition
was never about the algorithm — it is about whether the grid resolves the field. The strict gate is
still asserted where it is actually true, on the seven reference fields in `marching_cubes/tests.rs`. The exact
counts are pinned in both directions by
`an_under_resolved_pinch_makes_marching_cubes_non_manifold`, following M-4's precedent: the defect is
an assertion, not an exclusion, so it fails if it spreads *and* if it silently disappears.
**Would be shown wrong by:** the same field at `h = 2/3` coming back manifold, or any *resolved* field
coming back non-manifold under Marching Cubes.
**Worth noting for G-001:** a chunk boundary is a place where a surface can be under-resolved relative
to the chunk it lands in. This is a plausible source of seam defects later.

**The counterexample stands; the mechanism was wrong, and A-015 proved it.** The falsification
condition above — "the same field at `h = 2/3` coming back manifold" — **is now met**, by a change that
touched no geometry at all: the union of three spheres reports **0 non-manifold edges and 0 non-manifold
vertices** once cycles are fanned from a chord-safe apex.

That settles the attribution, because *a triangulation cannot repair a pinch*. If two sheets genuinely
met inside a cell they would still meet however each is triangulated. What was actually happening is
the fan chord of ✗17: two cells choosing the same interior diagonal, four faces on one mesh edge. The
refinement table above is real but was reading a proxy — refinement changed which sign patterns
occurred, not whether the grid "resolved" anything.

The reading of χ was wrong too, and in a way worth keeping: a collided mesh edge is counted **once** in
`E` while carrying four faces, so `E` was short by exactly the two collisions and `χ` was long by two.
The old `χ = 2` was not a topology measurement. The same fixture now reports `χ = 0`, genus 1 — three
lobes meeting in a genuine handle at that spacing.

**Consequence:** the property suite's Marching Cubes gate goes back to the strict `SurfaceGate::Closed` it was
waived from, and passes 8,000 generated cases where before A-015 it failed on the first fresh seed.
Whether Marching Cubes is *unconditionally* manifold is now **open, not settled** — see O-12. Nothing
here proves it; only that the one mechanism ever exhibited is gone.

### ✗17 — "Only the interior test can make Marching Cubes 33 non-manifold, so the face decider alone cannot"

**Believed because:** it was written down as a prediction in A-002's plan, and the reasoning looked
tight. Custodio et al. (2013) §6.2 give exactly one mechanism for Marching Cubes 33's non-manifold edges — *"two
adjacent voxels that share an ambiguous face have tunnels in the voxel interior"* — and tunnels come
from the interior test, which A-002 deliberately does not implement. The face decider only re-pairs cut
edges; it moves no vertex and creates no new one.
**Falsified by:** `every_closed_reference_field_meshes_cleanly_under_the_decider`, on its first run.
**Result:** the capped gyroid at 25³ gives **2 non-manifold edges and 3 non-manifold vertices** under
the decider against plain Marching Cubes' **0**, on a mesh with zero boundary edges and zero
inconsistently oriented edges. At 17³ and 33³ the two rules agree exactly.

**What's true instead — and it is not the decider's fault.** Both offending mesh edges carry **four**
faces. Inspected rather than reasoned about: each joins two cut cube edges that lie on a shared cube
face but are *not* connected by a segment there, and each is emitted **twice by each of the two
adjacent cells**. They are **fan chords**. `marching_cubes/table.rs`'s `triangulate` fans each cycle from its first
edge, and any triangulation of a `k`-gon that adds no vertices has `k − 3` interior chords; nothing
local stops two neighbouring cells choosing the same one.

An exhaustive two-cell search settles the attribution. Two cells stacked along z share a face and have
twelve samples between them, so all 4,096 sign patterns fit in a loop:

| rule | worst faces on one mesh edge | patterns affected |
|---|---|---|
| separated (A-001) | **4** | **12 / 4096** |
| decider (A-002) | **4** | **12 / 4096** |

Identical. **Plain Marching Cubes has this defect at exactly the same rate.** The decider only changes
*which* sign patterns are reached, and on the gyroid at 25³ it happens to reach one.

**Consequence:** the defect is ticketed against Marching Cubes as **A-015**, not against Marching Cubes 33, and the
fix is a per-cycle centroid vertex — the only chord-free triangulation of a polygon. That breaks
✗1/M-2/M-22's `V_mc = C` identity and re-baselines every golden hash, which is why it is its own
ticket rather than a fix folded into this one. Meanwhile `is_closed()` is deliberately *not* the gate
in the decider's reference-field sweep, because it folds in manifoldness and manifoldness is owned by
the fan; the census is pinned exactly instead, following M-4.
**Fixed at A-015, and for far less than it looked like it would cost.** The obvious repair — a centroid
for every cycle of four or more — works and is ruinous: **+73% vertices and +74% triangles** across the
seven reference fields, up to +99.7% on `box_exact`, because almost every cycle is long enough to
qualify. The cheap repair follows from asking *which* chords can collide rather than eliminating all of
them: only a cell containing **both** of a chord's cut edges can emit it, and two cells share a pair of
cube edges only when those edges share a cube **face**. So a fan is safe when none of its chords joins
two edges of one face, and that is a local, `const`-evaluable test. Measured over all 256 cases and
every canonical mask, **a safe apex exists for every cycle of length 3–7 and for 48 of the 60 length-8
cycles**; only lengths 9 and 12 have none, and plain Marching Cubes tops out at 7. Final cost across the
whole golden fixture: **one row's counts changed** — `marching_cubes+decider/gyroid/25`, +6 vertices and +12 triangles.
`V_mc = C` is intact.

**Would be shown wrong by:** the two-cell search returning a different count for the two rules, or the
offending mesh edges turning out to be face segments rather than chords.
**Method note:** the prediction was written into the plan *before* the code ran, which is the only
reason this reads as a falsification rather than as a bug that was quietly fixed.

### ✗14 — "Surface Nets is the cheapest thing in the family and the natural default"

**Believed because:** this repo's own algorithm catalog states it as the engine verdict —
`docs/research/2026-08-10-meshing-algorithm-catalog-v2.md:163`, *"cheapest thing in the family and the
natural default"* — reinforced by the same folklore ✗1 already corrected once.
**Tested by:** T-006's resolution sweep, `cargo bench --bench resolution_sweep`. Sphere, `f32`, single
thread, median of 5 timed runs after 2 warm-ups, identical grid and reused output buffers for both
algorithms. Raw data committed at `docs/measurements/resolution_sweep.csv`.
**Result:** Surface Nets is cheaper only below about 48³. The crossover sits between 48³ and 64³, and
past it Surface Nets loses steadily and then sharply:

| n | Marching Cubes ms | Surface Nets ms | Surface Nets/Marching Cubes | Marching Cubes ns/sample | Surface Nets ns/sample |
|---:|---:|---:|---:|---:|---:|
| 16 | 0.090 | 0.038 | **0.42** | 21.88 | 9.17 |
| 32 | 0.601 | 0.296 | **0.49** | 18.33 | 9.04 |
| 48 | 1.251 | 0.976 | **0.78** | 11.31 | 8.82 |
| 64 | 2.246 | 2.425 | 1.08 | 8.57 | 9.25 |
| 128 | 10.195 | 20.006 | 1.96 | 4.86 | 9.54 |
| 192 | 33.898 | 70.432 | 2.08 | 4.79 | 9.95 |
| 256 | 80.257 | **221.223** | **2.76** | **4.78** | **13.19** |

**What's true instead:** the two curves are not parallel — one converges and the other degrades.
Marching Cubes' per-sample cost *falls* from 21.9 to 4.78 ns as the `O(n²)` surface term amortises
away, then is flat from 128³ on. Surface Nets' is flat at ~9 ns to 128³ and then *rises*, reaching
13.19 ns/sample at 256³ — a 33% jump from 192³ alone.

**Consequence:** taken with ✗1, which showed Surface Nets emits `2χ` *more* triangles rather than
fewer, **both halves of the case for Surface Nets as the default are now falsified by measurement in
this repo.** What it actually retains is quad connectivity and one vertex per cell — a topology
argument, not a cost one. M-001's shootout must not present it as the cheap baseline, and the choice of
a default extractor for the game should be revisited.

**Would be shown wrong by:** a machine where Surface Nets' per-sample cost stays flat to 256³, which
would localise this to one cache hierarchy rather than to the algorithm.

**Amended again the same day (M-285), and this time the conclusion itself moves.** A-023 made
`emit_quads`' axis a `const` generic — three monomorphisations of one function, same order, byte-identical
mesh — and Surface Nets went **693.8 → 224.4 ms** at 256³, `SN/MC` **5.43× → 1.72×**, IPC **1.20 → 2.77**.
At 48³ Surface Nets is now **faster** than Marching Cubes, 30.7 cycles per sample against 33.1. So the
sentence *"Surface Nets never wins on Zen 3, at any resolution"* is false, and the cost half of the case
against Surface Nets as a default is much weaker than this entry says: 1.7× at the largest grid, a win at
the smallest. **✗1's half is untouched** — Surface Nets still emits `2χ` more triangles — so the entry is
narrowed rather than reversed, and the residual gap is the cache story of M-286 rather than the algorithm.

**Amended 2026-08-16 (M-282), and the conclusion is stronger while every number here is superseded.**
Measured in one binary and one run on the Ryzen: `surface_nets / marching_cubes` at 256³ is **5.43×**
against the 3.72× above, and **3.19×** at 16³. Marching Cubes got 1.74× faster between `d2ab82a` and
now while Surface Nets got 1.17×, which is where the widening comes from. The mechanism is also no
longer "per-sample cost rises": M-279 measured that the whole family divides on **IPC** — everything
table-driven runs at 3.7–4.2 and everything on `DualMesher` at 1.20–1.42 — and that the dual's rise
with `n` is a 16% IPC decline on an instruction stream that is flat per sample. And per M-281, the
absolute milliseconds in this entry are comparable only against others from the same binary.

**Confirmed on a second machine at O-11, and harder — but the crossover is not (M-45).** On an AMD
Ryzen 9 5900X, Surface Nets' per-sample cost climbs `37.4 → 49.1 ns` over the same sweep, so the
falsification condition above is not met: the verdict is a property of the algorithm rather than of one
cache hierarchy. The ratio is *worse* there — **3.72× at 256³** against the M5's 2.65×.

What does **not** transfer is the sub-claim "Surface Nets is cheaper only below about 48³". That
crossover exists on the M5 only because Marching Cubes *starts* expensive there — 24.99 ns/sample at
16³ — and converges. On Zen 3 Marching Cubes is flat at 13–15 ns from 16³ up, so it has nothing to
converge from and **Surface Nets is behind at every resolution measured**, 2.46× even at 16³. Any
statement of the form "Surface Nets wins at small grids" is an Apple-silicon statement.

**Narrowed by A-007.** Dual Contouring, which shares the dual topology and differs only in vertex
placement, shows the *same* curve: `218.9 ms` at 256³ against Surface Nets' `212.5`, and the same
negative fitted intercept. So the superlinearity is a property of the **shared dual engine** — its
sampling and its strided quad walk — and not of either vertex rule. That is a real narrowing of O-11's
search: the suspect is the gather in `emit_quads`, whose `z`-stride is `n²` cells apart.

### ✗18 — "A hairline seam difference is a crack no weld can close"

**Source: this repository — M-73, and `transvoxel::cell`'s own module docs, which say it twice.**

> *"A weld cannot rescue it: welding merges vertices it can see are the same, and these two differ in
> the last bit."* (M-73)
>
> *"The weld can close a seam it can see; it cannot invent a shared vertex where the two sides
> disagree in the last bit."* (`cell.rs`, module docs)

**Falsified at R-004 (M-278).** *"Vertices it can see are the same"* reads as bit-identity, and the
welder's rule is not that: it is **first fit within `epsilon`**, and the one policy is
`epsilon_for(h) = h · 1e-4`. At `h = 1/12` that is `8.33e-6` against a worst measured seam
disagreement of `1.440e-15` — **nine orders of magnitude inside the tolerance**, and the 27-cell
lattice probe means straddling a bucket boundary does not save it either.

Measured on the two-resolution fixture, offset arithmetic, all five spacings and both LOD pairs:
**0 seam-plane boundary edges in all twenty rows under the crate's own weld** — including the twelve
where a bit-identity merge leaves **63 to 348**.

**Two things survive, and they are why the M-73 fix was still right.**

- **A consumer that does not weld gets the crack in full.** M-69 is that consumer on the record: *"a
  seam of unshared vertices that a renderer draws correctly and a collider reads as a hole."* Sharing
  by construction is structural; sharing by tolerance is a policy someone can change.
- **The weld is not free to lean on.** R-001 measured that it can *create* a non-manifold edge
  (M-226) and R-002 that its result depends on input order. *"The weld will catch it"* is a bad
  thing to design against even when it is true.

**And the claim understates the damage in the other direction, in 2 of 12 rows.** The offset
arithmetic is not always a hairline: on `torus` at `h = 1/12` the widest bit-identity crack is
**2.076 cells** at LOD 0–1 and **1.053** at 1–2, because a sample perturbed in its last bit crossed
zero and the two sides then disagreed about whether an edge was cut at all. That is a hole, not a
rounding error, and no weld epsilon short of a cell would close it.

### ✗19 — "Manifold Dual Contouring's uniform-grid surface is always a manifold"

**Source: the paper itself, read at A-022 (V-34).** Schaefer, Ju & Warren, §3:

> *"To create polygons, the algorithm constructs one polygon connecting the vertices associated with
> that edge in the four adjacent cells. … **this surface is always a manifold** because the original
> MC algorithm always constructs a manifold and the dual preserves the topology of the surface."*

Unconditional, and at the uniform-grid level — before any octree, any simplification, any clustering.

**Falsified at A-022 (M-290).** Eight reference fields × three resolutions, no chunking, no weld:

| extractor | non-manifold edges | non-manifold vertices |
|---|---|---|
| `marching_cubes` | **0** | **0** |
| `marching_cubes` + asymptotic decider | **0** | **0** |
| `manifold_dual_contouring`, unmodified table | 143 | 286 |
| `manifold_dual_contouring`, **decider-modified table** — the paper's own construction | **114** | **222** |
| `surface_nets` / `dual_contouring` | 1,128 | 2,156 |

**The interesting part is which half of the argument fails.** The premise is *"the original MC
algorithm always constructs a manifold"* — ✗15 records that as false in general, so the obvious guess
is that the premise is what breaks. It is not: **Marching Cubes measures 0 on all 24 configurations
here, under both face rules**, A-015 having removed the fan chords that made ✗15's counterexample.
The premise holds and the conclusion still fails, so what is false is *"the dual preserves the
topology of the surface"*.

**And the mechanism was already written down in this crate, as something the paper did not claim.**
`manifold_dual_contouring`'s module docs say the residue is *"two crossed edges of one **shared** face
resolving to the same cycle pair"*, producing four quads on one dual edge, and then: *"Schaefer, Ju
and Warren separate sheets within a cell and **never claimed** to handle [this], so this is outside
what they guarantee rather than a defect against it."* That charitable reading was written while the
paper was unobtainable. **They do claim it**, in one sentence, unconditionally. So the residue is a
defect against the paper's stated guarantee and not outside it, and the crate's own docs need the
correction — which they now have.

**What is not falsified.** The paper's *contribution* is the octree vertex-clustering algorithm and
its topology-preservation proofs, none of which this touches. What fails is the one-sentence argument
for the uniform-grid base case it builds on, and it fails on `noise_cavity` only — the field A-002e
added because none of the other seven produces a cell with an interior ambiguity (M-208).

---

## Part 2 — Measured here (tier M)

| # | Finding | Evidence |
|---|---|---|
| M-1 | **surface cells = crossed edges + χ** | 4 topologies × 3 resolutions, table in ✗1 |
| M-2 | `V_sn = V_mc + χ`, `F_sn = F_mc + 2χ` | A-004 tests, all four clean fields |
| M-3 | Surface Nets max vertex degree **10**; Marching Cubes **9** | A-004 |
| M-4 | Surface Nets is non-manifold where one cell carries two sheets: **48** non-manifold edges on capped gyroid, **15** on fbm_terrain at 33³ | A-004; pinned as non-zero assertions, not excluded silently |
| M-5 | On `box_exact`, Surface Nets' nearest vertex to the corner (1,1,1) is **1.15 cells** away | A-004 — this gap is what E-104 exists to show |
| M-6 | `libm::sqrtf` lowers to hardware `fsqrt` (aarch64+neon) / `sqrtss` (x86-64+sse2) | libm 0.2.16 source: `src/math/arch/aarch64.rs` raw asm, dispatched by `select_implementation!` on `target_feature` |
| M-7 | dev-dependencies do not propagate: consumer resolves **3 packages**, the crate's own lockfile has **137** | Experiment, cloud container |
| M-8 | Cargo silently co-resolves two wgpu majors — **317 packages, both 29.0.4 and 30.0.0**, no resolution error; fails later as `expected TextureFormat, found a different TextureFormat` | Experiment |
| M-9 | Workspace feature unification leaks: `-p isomesh` alone gives glam `libm`; whole-workspace gives it `std`, `serde`, `bytemuck`, `encase`, `rand` | Experiment — the reason `bevy_isomesh` is excluded |
| M-10 | **Unit sphere at 64³ (`h = 4/63`), symmetric Hausdorff: Marching Cubes `1.380e-3`, Surface Nets `2.288e-3`.** Mean absolute error Marching Cubes `6.50e-4`, Surface Nets `1.367e-3`. Surface Nets is **1.66×** worse than Marching Cubes on both | T-003, `a_unit_sphere_at_64_cubed_is_within_one_cell_diagonal` |
| M-11 | **T-003's own acceptance criterion is loose by ~80×.** One cell diagonal is `0.10997`; Marching Cubes measures `0.00138`. A harness returning a constant `0.01` would pass it | T-003 — which is why the ticket also ships a convergence-order test and closed-form fixtures |
| M-12 | **Marching Cubes' error falls like `h²`, measured.** Mean error `2.7168e-3` at 32³ against `6.5015e-4` at 64³ — a ratio of **4.179**, against the ideal `((4/31)/(4/63))² = 4.13` | T-003, `the_error_falls_like_h_squared` |
| M-13 | **Surface cells ≈ `1.5·A/h²`, not `A/h²`.** Measured `1.450` (25³), `1.442` (33³), `1.517` (64³) on the unit sphere. The constant is derivable: a plane of unit normal `n` crosses `(\|nₓ\|+\|n_y\|+\|n_z\|)/h²` cells per unit area, and `E[\|nₓ\|] = ½` over the sphere, so an isotropic surface gives `E[Σ\|nᵢ\|] = 3/2` | T-003. Predicted 6,430 triangles at 64³ from `A/h²` and measured **9,452** — a 1.47× miss, which is this factor |
| M-14 | **The reverse direction finds defects the forward direction structurally cannot.** `box_exact` at 33³: forward `0.0833`, reverse `0.1443` — the reverse number is Marching Cubes' rounding of the sharp corner. `thin_plate` at 33³: forward `0.0083`, reverse `0.0893` — an under-resolved plate | T-003. Deleting one face of an octahedron leaves `mesh_to_field` bit-identical and moves `field_to_mesh` to `√(3/2 − 2/√3)` |
| M-15 | **Surface Nets' non-manifoldness is a resolution effect, not a topology effect.** M-4 measured it on `gyroid` (48 edges) and `fbm_terrain` (15) and read it as a high-genus / open-field property. T-005b finds it on a randomly generated **convex body** — 1–2 non-manifold edges, 3–4 non-manifold vertices, zero boundary edges. Any feature thinner than one cell forces two sheets through it | T-005b, `surface_nets_meshes_convex_bodies`. This is why the sweep has a named `SurfaceGate` rather than a per-field exception |
| M-16 | **The even-`χ` parity check is not independent of manifoldness — it is a corollary of it.** `χ = 2 − 2g` holds for a closed *orientable manifold*, so a gate that waives manifoldness and keeps parity is incoherent. Measured: Surface Nets on a generated convex body gives **`χ = 1`** with one non-manifold edge and zero boundary edges | T-005b. Cost one wrong gate before it was noticed; `SurfaceGate::ClosedAllowingMultiSheet` now documents the omission rather than leaving it to be rediscovered |
| M-17 | **A case-table entry naming an *uncut* edge is caught inside the crate, before any mesh exists** — `edge_crossing`'s `is_inside(a) != is_inside(b)` precondition fires. Real defence, but it is a `debug_assert`, so it is absent from a release build | T-005b. The mutation check therefore confines its wrong-edge corruption to *cut* edges, which is both the plausible transcription error and the one that actually reaches the validity gate |
| M-18 | *(refined by T-008 — the arithmetic below is about **adjacent cells**; the effect on a real mesh is **gradual**, because a mesh whose vertices span many lattice cells keeps resolving them for a while after neighbours have merged. Measured on a 1158-vertex sphere at `h = 0.125`: no collapse at all at `1e4`, and **1158 → 918 buckets**, a 21% loss, at `1e6`. Fixed by anchoring the lattice to the mesh's own minimum corner, so the scale depends on the mesh's extent rather than its position)* **`quantise`'s weld lattice collapses beyond ~105 world units, and it is a performance cliff rather than a correctness one.** It scales absolute coordinates by `1/(h·1e-4)` — 160,000 at `h = 0.0625` — so it passes `f32`'s exact-integer range at `2²⁴·weld_epsilon ≈ 104.86`. Measured: at `p = 104` two cells one apart stay distinct; at `p = 105` they collapse; by `p = 1000` a whole region is one bucket. **Correctness survives** — coarsening only *merges* buckets, and the 27-neighbour probe plus exact distance test still finds every duplicate — but the scan degrades toward quadratic, silently, at exactly the coordinates G-001 chunking and G-007 streaming produce. `TriangleGrid` is immune because it quantises *relative* to its own AABB origin, which is the fix pattern | T-005b follow-up, ✗13. Ticketed as T-008 |
| M-19 | **There is no meaningful fixed cost on the CPU extraction path, and the prediction saying so was written down before the run.** Marching Cubes' fitted `a` is `0.49 ms` against a largest measured run of `80.3 ms` — **0.61%**. V-6's "73% of a published 64³ figure was fixed launch overhead" is a *GPU dispatch* property and does not transfer; the "stop trusting single-grid numbers" rule belongs to Phase 6, not here. **Caveat that matters:** `a` is 543% of the *smallest* measured run, so the fit must not be extrapolated below 16³ — down there the `O(n²)` surface term dominates | T-006, `benches/resolution_sweep.rs`. The prediction is in that file's module docs, committed before the first measurement |
| M-20 | **Marching Cubes' marginal cost is `4.75 ns/sample` — `211 M samples/s`, single-threaded, `f32`, Apple M5.** Per-sample cost is flat within 2% from 128³ upward, `r² = 0.99986` | T-006. Against V-3's `5.42 G voxel/s` on an RTX 2080 Ti that is a **~26× gap**, which is the number the Phase 6 GPU decision should be argued from rather than from folklore |
| M-21 | **Surface Nets is not `O(n³)` over this range; Marching Cubes is.** Surface Nets' fitted intercept is **negative** — `−3.13 ms` full sweep, `−7.32 ms` on the tail — which is physically impossible and is the signature of a curve convex in `n³`. `r² = 0.9899` against Marching Cubes' `0.99986`. Per-sample cost rises `9.0 → 13.19 ns` while Marching Cubes' falls and flattens | T-006. Cause unmeasured — see O-11. This is why ✗14's gap widens rather than staying constant | **Amended at O-11: the intercept-sign diagnostic is machine-specific; the per-sample rise it stands in for is not.** On Zen 3 *both* fitted intercepts come back negative and both are numerically negligible — −0.04% of the largest run either way, `r² = 0.99999` — so the sign diagnoses nothing there. What reproduces is the underlying effect: Surface Nets' per-sample cost rises **31%** across the sweep (37.4 → 49.1 ns) while Marching Cubes' *falls* 13% (15.2 → 13.2). Report the per-sample curve, not the intercept. **Zen 3 figures re-measured 2026-08-16 (M-282): 29.63 → 41.35 ns for Surface Nets, a 40% rise, and 9.29 → 7.62 for Marching Cubes, an 18% fall** — same shape, both faster, and the gap wider. **Then A-023 moved the Surface Nets curve bodily (M-285): 7.92 → 13.37 ns over the same range**, so the rise is still there and the constant is 3× smaller. What remains rising is the cache term M-286 can now see. The M5 half of this row has not been re-run (M-005)
| M-22 | **✗1's identity holds at every resolution to 256³**: `V_sn − V_mc = 2` and `F_sn − F_mc = 4` exactly, nine resolutions, `χ = 2`. The original table topped out at 49³, so this is corroboration at **5× the resolution** and 16.8 M samples | T-006's sweep, which records vertex and triangle counts alongside the timings |
| M-23 | **`f64` costs 8–10% on extraction paths with no matrix solve in them.** At 65³ on a sphere: Marching Cubes `1.3928 ms` (f32) against `1.5083 ms` (f64), **+8.3%**; Surface Nets `2.3625` against `2.6036`, **+10.2%**. Not the 2× a naive "twice the bytes" guess suggests, because the work is dominated by field evaluation and branchy table lookup rather than by memory bandwidth | T-006, `benches/extract.rs`, the `precision` group. **Partially answers O-8** for the non-QEF paths; A-007's solve is where `AᵀA` squares the condition number and the answer may differ |
| M-24 | **Bit-exact lattice equivariance needs magnitude-ordered *products*, not just sums.** The audit prescribes "magnitude-sorted 3-term dot products", which is necessary and **not sufficient**: a cofactor expansion of `det(M+λI)` along a fixed row selects three of the six entries *by position*, so relabelling the axes evaluates a different expression. Measured **19 ULP** disagreement under a cyclic permutation, on all three fixtures, with the dots already sorted. Fixed by the symmetric determinant form with magnitude-ordered 3-factor products — FP multiplication is commutative but not associative, so `(a·b)·c ≠ (b·c)·a`. Now **72/72** rotation×fixture cases are bit-identical | A-007, `the_vertex_is_bit_exactly_equivariant_under_lattice_rotations`, which failed before the fix |
| M-25 | **The sharp-feature solve is nearly free: Dual Contouring costs 3% more than Surface Nets.** At 256³ on a sphere, `218.9 ms` against `212.5 ms`; marginal `78.1` against `80.4 M samples/s`. A full 3×3 regularized solve per surface cell, and it barely registers — because both methods are dominated by the *shared* dual topology (sampling and the quad walk), not by vertex placement | A-007, `benches/resolution_sweep.rs`, now sweeping three algorithms | **Confirmed on a second machine at O-11**, at roughly double the fraction: on Zen 3 Dual Contouring costs **6.5%** more than Surface Nets (`877.4` against `823.4 ms` at 256³). Small on both, so the conclusion holds and the exact percentage is a machine property
| M-26 | **Dual Contouring reaches a box corner to `0.01` cells where Surface Nets stops at `0.58`** — measured at 27³ on `box_exact`, `0.0009` against `0.0888` in world units. The resolution is deliberately **not** grid-aligned; on an aligned grid this measures the zero-classification rule instead (E-103's trap) | A-007, `the_corner_is_sharper_than_surface_nets`. This is E-104's money shot, measured before the example exists |
| M-27 | **The two dual methods differ *only* at features, with a 14-order-of-magnitude gap and nothing in between.** On `box_exact` at 27³: **864** of 1016 vertices agree with Surface Nets to within `2e-15` cells, **152** move by `0.35`–`0.57` cells, **0** land between. Exact reason: on a planar patch every crossing lies in the plane, so `pᵢ − c ⊥ n`, every `dᵢ` is exactly zero, `g` is exactly zero and the solve returns the centroid | A-007, `dual_contouring_moves_only_the_feature_vertices`. Consequence: E-104's side-by-side measures the feature and nothing else. Note the agreement is to *rounding*, not to the bit — the two centroids are computed by different expressions |
| M-31 | **`libm` delivers the bit-identical cross-platform meshes it was chosen for — verified, not reasoned.** T-007's 63 golden hashes were generated on **macOS / arm64** and pass unchanged on **Linux / x86-64** in CI. Every position, normal and index bit-for-bit equal across two architectures, two operating systems and two libm-vs-hardware float paths | I-001 chose `libm` over `std` unconditionally on the grounds that `std`'s `sin`/`cos` are the platform's own and differ between macOS and Linux, and recorded that T-007's committed hashes would be the proof. They are. The claim moves from a design argument to a measurement |
| M-32 | **Chunk seams are bit-exact only when the cell size is a power of two.** Two adjacent chunks meshed independently agree on **16 of 16** shared-plane vertices bit-for-bit at `h = 0.125`, and on **0 of 14** at `h = 4/35` — worst gap `1.57e-16` world units, `1.37e-15` cells. Cause: an extractor computes `origin + h·local`, so chunk `c`'s last plane is `(o + h·cn) + h·n` while `c+1`'s first is `o + h·(c+1)n` — equal by algebra, not by IEEE. **22% of 200,000 random `(origin, h, cells, chunk)` combinations disagree**, by one or two ulp | G-001. A rounding error rather than a crack, but it decides whether a chunked world can claim the bit-identity T-007 and M-31 establish for a single volume. **Recommend power-of-two cell sizes for chunked worlds**, and note G-003's commutativity acceptance will depend on this |
| M-33 | **E1, the unpublished number: a brush changes `15–36%` of the cells in its own bounding box.** Measured for a sphere brush carved two cells deeper, `h = 0.0625`: radius `0.25` → **36.0%**, `0.50` → **29.2%**, `0.75` → **24.4%**, `1.00` → **14.9%**. It falls roughly as `1/r`, and the reason is geometric — the changed set is a shell of fixed thickness, so it grows as `r²` inside a bounding box that grows as `r³` | G-002. The research doc asks for this and records that nobody has published it. **The incremental case is real but not overwhelming**: a bounding-box re-mesh does 3–7× the necessary work, not 30× |
| M-34 | **Counting *value* changes overstates the re-mesh set by 2.8–3.7×.** An SDF brush perturbs values throughout its support — carving a sphere deeper moves the field everywhere the sphere term dominates, including deep inside the solid — but a cell wholly inside or wholly outside emits no triangles either way, so its **output** is unchanged. Measured at radius `0.25`: **100%** of cells changed value, **36%** changed output. At `1.00`: 55% against 15% | G-002. The first version of E1 counted value changes and read `100%`, which would have been published as "incremental meshing buys nothing" — the exact opposite of what the data says |
| M-35 | **A brush stepping two cells sweeps `6–14%` of its sign-changed cells entirely through.** Every corner flips from outside to inside in one edit, so the cell is a surface cell neither before nor after and needs no re-mesh at all. Measured 167 of 1155 at radius `0.25`, 1040 of 8660 at `1.00`. Consequence: **a sign change does not imply an output change**, which looks obviously false and is why the assertion caught it | G-002. The signature of an edit moving further than one cell per step, which is where a stroke starts skipping geometry between frames |
| M-36 | **The multiplayer story survives, with a boundary. Eight brushes in all `8! = 40,320` orderings give exactly *one* result — bit-identical — when they are all `Add`, and again when all `Subtract`.** `min` and `max` are commutative *and* associative in IEEE, and they introduce no rounding at all: they select an argument rather than computing a value. Concurrent clients may reorder a run of same-kind hard edits freely and converge bit for bit | G-003, the ticket's acceptance criterion. **No fixed-point storage was needed to achieve it**, which is worth recording because the ticket proposed `i32` storage for exactly this guarantee — see M-39 |
| M-37 | **Mixed add and subtract do *not* commute: 11 distinct results from the same 40,320 orderings.** Carving a hole and then filling it is a different solid from filling and then carving. **Semantic, not numerical** — no storage format or arithmetic repairs it, and a concurrent-editing protocol must preserve order across an add/subtract boundary while remaining free to merge within a run | G-003. `BrushOp::commutes_with` returns the honest answer rather than the optimistic one |
| M-38 | **Smooth union destroys reordering almost completely: 40,317 distinct results from 40,320 orderings.** Smooth-min is **not associative** — measured `1.694e-2` apart on `smin(smin(a,b),c)` against `smin(a,smin(b,c))`, four orders above rounding — and it is not *bit*-commutative either, disagreeing by **1 ulp** when its arguments are swapped. So "all the same operation" is not sufficient for reordering; only "all the same *hard* operation" is | G-003. The two failures are independent and have different fixes: fixed-point would repair the 1-ulp one and nothing repairs the associativity one |
| M-39 | **Fixed-point storage is unnecessary for the guarantee it was proposed for.** G-003's ticket asks for an `i32` option "for bit-exact determinism", but `min` and `max` are already exact in IEEE — the ordering guarantee in M-36 holds in `f64` with nothing added. The only place fixed-point would help is smooth-min's 1-ulp commutativity gap (M-38), and that is the lesser of its two reordering failures | G-003. Not implemented, on the evidence: `Real` is a sealed trait and adding a fixed-point scalar is a large change to buy a property already held |
| M-28 | **The cell clamp eliminates placement-caused self-intersections entirely, and costs nothing in sharpness.** λ (pairs per 1,000 triangles) at 33³, clamp off → on: `torus` **2.66 → 0**, `gyroid` **71.43 → 3.12**, `fbm_terrain` **189.46 → 13.84**; `sphere`, `box_exact`, `csg_difference` and `thin_plate` were already 0. Corner gap on `box_exact` at 27³: **0.0057 cells either way** — identical, because a convex corner's solution is interior to its own cell, so the constraint never binds where the feature is | A-009, `the_clamp_measured_on_every_reference_field`. Default is now `Clamp::ToCell`, chosen by this measurement rather than by preference |
| M-29 | **The literature's two branches both fire, on disjoint fields — which is a sharper answer than either alone.** The review states the rule in advance: λ→0 means placement was the cause, λ unchanged means the defect is *connectivity* and needs A-010. Measured: λ → **exactly 0** on five of seven fields, and drops **23×** and **13.7×** on `gyroid` and `fbm_terrain` without reaching it. Those two are precisely the fields with multi-sheet cells (M-4, M-15). So the clamp removes the placement failure completely and the residue is exactly A-010's problem, with nothing left unaccounted for | A-009 |
| M-30 | **An unclamped solve can fling a vertex 3.18 cells out of its own cell** — measured max displacement on `gyroid` at 33³, with 618 of 5240 vertices outside; `fbm_terrain` 2.17 cells and 1097 of 1958. On the smooth closed fields it never leaves at all: `sphere`, `box_exact` and `thin_plate` have **zero** vertices outside | A-009. This is the failure mode the clamp exists for, quantified rather than asserted |
| M-40 | **The ambiguous face is rarer than the literature suggests — on five of the seven reference fields it never occurs at all.** At 33³: `sphere` 0 of 1160 surface cells, `torus` 0 of 1128, `box_exact` 0 of 1352, `csg_difference` 0 of 1388, `thin_plate` 0 of 512. Only `gyroid` (**27 of 5240, 0.515%**) and `fbm_terrain` (**30 of 1958, 1.532%**) reach it, and the decider joins roughly half the ambiguous faces it finds — 12 and 18 respectively. **So Marching Cubes 33 and Marching Cubes are bit-identical on five of seven fields at every resolution tested**, which the 84-row golden fixture now pins | A-002. Verifies Custodio et al. 2013's *"the vast majority of Marching Cubes cases match the non-ambiguous configurations"* (tier R) against this crate's own fields and finds it understated. Consequence for E-102: an Marching Cubes-vs-Marching Cubes 33 example must use `gyroid` or `fbm_terrain` or it will show two identical meshes |
| M-41 | **88 of the 256 cases change their Euler characteristic when their ambiguous faces are joined.** The smallest is case 6 — corners 1 and 2 inside, diagonally opposite on the `z = 0` face and on no other face together, so exactly one ambiguous face: separated it is two discs (`χ = 2`), joined it is one (`χ = 1`). 136 of the 256 cases have no ambiguous face at all, so the rule cannot reach them. The decider's worst cell is **10 triangles** against the separated table's 5 — predicted before running, from the longest cycle that can use all twelve cut edges | A-002's acceptance criterion, `the_decider_and_marching_cubes_disagree_about_chi`. Searched over all 256 cases rather than picked, and pinned in both directions |
| M-42 | **The asymptotic decider is free to within a few percent, which is the first time this repo's "~free" claim has had a benchmark behind it.** Median extraction, f32, Apple M5: `sphere` 33³ **206.25 → 205.65 µs** (−0.3%, i.e. noise), `sphere` 65³ **1.4954 → 1.5189 ms** (+1.6%), `gyroid` 33³ **786.89 → 795.44 µs** (+1.1%), `gyroid` 65³ **5.6378 → 5.8236 ms** (+3.3%). `sphere` has no ambiguous face at all (M-40), so its difference is the price of *asking* — one table lookup and a branch per surface cell; `gyroid`'s extra ~1.7 points is the price of *answering*, building the cell's triangulation at run time instead of reading it | A-002, `cargo bench --bench extract -- decider`. Confirms the v1 catalog's "~free" (tier R) for the decider, against its "730 subcases in the LUT" for the guaranteed version |
| M-43 | **The decider needs no division and no epsilon, and the brief's "guard the denominator" is unnecessary.** On an ambiguous face one diagonal is strictly negative and the other non-negative — a sample of exactly zero is outside — so `v0 + v2 − v1 − v3` is strictly non-zero *by the sign rule alone*. Only `sign(S)` is wanted and the denominator's sign is already known, so the whole test is **`joined ⟺ d_in > d_out`** on the two diagonal products. Both branches of the derivation reduce to the same comparison, and it is invariant under rotation and reflection of the corner order because IEEE multiplication is commutative and correctly rounded — which is what makes two adjacent cells agree bit for bit | A-002. Structurally the same argument as `edge_crossing`'s missing epsilon: strictness in the sign rule pays for itself twice |
| M-44 | **The decider does not widen M-32's chunk-seam problem, and there is margin to spare.** Over 217 seam planes where the two chunk expressions differ bit for bit and 499,968 faces lying in them: **0** where the ulp moved a corner across zero (which would be a crack for plain Marching Cubes too, not just for the decider), 205 ambiguous faces, **0** decision flips. The closest any ambiguous seam face came to its own decision boundary was a relative margin of **1.535e-2** — about fourteen orders of magnitude above the `~1e-16` perturbation the seam arithmetic introduces | A-002, `the_decider_at_a_chunk_seam_is_measured`. A count of zero says nothing about how nearly it happened, which is why the margin is recorded alongside it. The first sweep found **0** ambiguous seam faces and the test's own reachability gate caught it — the fixture trap for the third time |
| M-45 | **✗14 reproduces on a second machine and gets worse; its crossover does not reproduce at all.** Same sweep, same field, same commit, AMD Ryzen 9 5900X (Zen 3, x86-64, single thread) against the Apple M5. Per-sample cost, `16³ → 256³`: **Marching Cubes M5 24.99 → 4.78 ns, Marching Cubes Zen 3 15.18 → 13.19; Surface Nets M5 8.40 → 12.66, Surface Nets Zen 3 37.38 → 49.08.** So Surface Nets degrades on both — the effect is the algorithm's memory pattern, not one cache hierarchy — and `Surface Nets/Marching Cubes` at 256³ is **3.72× on Zen 3 against 2.65× on the M5**. But Surface Nets never wins on Zen 3, at any resolution: `2.46×` behind even at 16³. The M5 crossover exists only because Marching Cubes starts expensive there and converges, which Zen 3's Marching Cubes does not do because it is flat from the start | O-11, `cargo bench --bench resolution_sweep` on `big` at commit `d2ab82a`. Raw data committed as `docs/measurements/resolution_sweep-ryzen9-5900x.csv` beside the M5's. Also: the M5 is **2.76× faster than the Ryzen on Marching Cubes at 256³** (80.2 vs 221.4 ms) single-threaded, while the Ryzen is faster below about 32³ | **Amended 2026-08-16: that cross-machine figure is withdrawn and the Zen 3 half is superseded (M-280, M-282).** A worktree at `d2ab82a` re-run on the same Ryzen tonight reproduces this row to within 1.8% (222.649 ms against 221.363), so the numbers were sound — and Marching Cubes has since got **1.74× faster** (127.8 ms) while Surface Nets got 1.17×. The M5 half was not re-run, so `2.76×` compares a current Ryzen against a stale Mac and must not be quoted; the within-machine claims — Surface Nets degrading on both, no crossover on Zen 3 — are unaffected. **Except the second one, which A-023 falsified the same day (M-285): Surface Nets is now faster than Marching Cubes at 48³ on Zen 3, 30.7 cycles per sample against 33.1, so the crossover this row says does not exist here now does** |
| M-46 | **A chord is only collidable when its two cut edges share a cube face, and that makes the manifold fix nearly free.** The naive repair — centroid-fan every cycle of four or more — costs **+73.1% vertices and +73.8% triangles** over the seven reference fields at three resolutions (23,034 → 39,881 and 45,662 → 79,356), worst on `box_exact` at **+99.7%**, because nearly every cycle qualifies. Restricting it to cycles with no chord-safe apex costs **+6 vertices and +12 triangles, on one row of eighty-four**. The enabling fact is local: only a cell containing both of a chord's cut edges can emit that mesh edge, and two cells share a pair of cube edges only if the edges share a face. **A safe apex exists for every cycle of length 3–7 and 48 of the 60 length-8 cycles; plain Marching Cubes never exceeds length 7**, so it never pays anything and `V_mc = C` survives | A-015. The naive version was implemented and measured first, which is the only reason the cheap one was looked for — the ticket had been written expecting to re-baseline ✗1/M-2/M-22 and the whole golden fixture |
| M-47 | **The validator's `duplicate_vertices` is an upper bound on what a weld removes, not the count.** It asks whether *any* earlier vertex is within ε; the welder asks for the lowest-indexed *kept* one, and the two part company wherever a chain of near-misses exists — the validator counts the middle of a chain as a duplicate of its start, the welder leaves the end of the chain unwelded. Predicted before running and both halves held: **equal on a real chunk seam** (14 and 14), where duplicates are pairs an ulp apart with no chains, and **different on a constructed chain** (2 against 1). Also measured, two chunks of a unit sphere at `h = 4/35`: `273 → 259` vertices, boundary edges `85 → 59`, `duplicate_vertices 14 → 0`, `χ 2 → 1` — two discs glued along an arc | A-013. The two share one `Lattice` so they cannot disagree about *which* cells are probed; what differs is the question, and that difference is now measured rather than assumed |
| M-48 | **The edge-vertex cache does not share everything, and welding removes a class of sliver nobody expected it to.** The cache shares vertices between cells meeting on a grid *edge*; when a grid **sample** lands on the isosurface, `t` is 0 or 1 and the crossing sits *at that sample*, so every cut edge meeting there places its own vertex at the same point and nothing shares them. Whole-volume weld census over the seven fields at 17/25/33³: `sphere` 25³ **48 vertices and 96 triangles**, `gyroid` **2 and 4** at all three resolutions, `fbm_terrain` 33³ **1 and 2**, everything else **zero**. The 96 is exactly the degenerate-sliver count A-001 measured at that resolution from the 30 lattice points sitting exactly on the unit sphere — the same 96, so **welding is a fix for that class of sliver** | A-013. Falsified a claim written in `weld.rs`'s own module docs the same day — "every reference mesh reports `duplicate_vertices == 0`" — which was asserted from the edge cache's design rather than measured. The test that disproved it was written to *confirm* it |
| M-49 | **`ChunkLayout::cell_of` inverts `world_of_sample` inside a cell and not reliably on its corner — M-32 in a second place.** `world_of_sample` computes `origin + h·sample`, `cell_of` computes `floor((p − origin)/h)`; inverse by algebra, not by IEEE. Measured over three cells: **3 of 3 corners round-trip at `h = 0.125`, 1 of 3 at `h = 4/35`**, where the division lands a hair under the integer and `floor` takes it down. The interior of a cell is unambiguous at any spacing. **No epsilon was added** — a point exactly on a cell boundary belongs to either cell by convention, and at a non-power-of-two spacing "exactly on" is not a decidable question, so snapping would trade a visible ambiguity for an invisible one | A-013's `cell_of`, added for E-202. Callers needing a cell *range* pad it; E-202 pads by the cell size. The test was written to assert a round trip and immediately caught that it does not hold |
| M-50 | **E1 and M-34's ratio both reproduce live, under a mouse.** E-202 carves with a brush and re-meshes only the dirty chunks, reporting per edit: a typical carve is **265 of 1,728 cells in the brush's bounding box = E1 15.3%**, against **756 cells whose sample value moved** — a ratio of **2.85×**, inside M-34's measured 2.8–3.7×, and E1 inside M-33's 15–36%. Over a scripted 60-carve run E1 ranges **0.6% to 27.3%**. Cost per re-meshed chunk against edit-log length, median: **0.158 / 0.354 / 0.525 / 0.589 ms** for logs of 1–15 / 16–30 / 31–45 / 46–60 — **3.7× for 7× the log, and flattening.** So the `BrushStack` walk is a real cost and *not* proportional at these lengths, which is weaker than "every sample walks every brush" suggests | E-202, `ISOMESH_AUTOCARVE=60`. The first offline measurements of E1 (M-33) and of the value-versus-output ratio (M-34) were made on synthetic edits; this is the first time either has been measured on the interactive path they were written to justify |
| M-51 | **Marching Tetrahedra costs ~3× the triangles for ~4% worse geometry — and the literature's `2–3×` is too low.** Vertex and triangle ratio against Marching Cubes on identical grids, seven reference fields at 33³ and 49³: `gyroid` and `fbm_terrain` **2.87×**, `sphere` and `torus` **3.04×**, `csg_difference` **3.83×**, `thin_plate` **3.84×**, `box_exact` **3.91×**. The tier-R figure from `10.1109/2945.485620` covers only the two roughest fields. On the other side, Lewiner et al. 2003's *"weaker geometrical accuracy… the vertex position cannot be adjusted to fit the geometrical trilinear approximation"* measures **4.3%**: symmetric Hausdorff on a unit sphere at 64³, marching cubes `1.3798e-3` against marching tetrahedra `1.4386e-3`. Directionally right and far weaker than it reads | A-003. The marching cubes figure reproduces **M-10's recorded `1.380e-3` exactly**, so the harness is measuring what it measured before. **P-1's `2.992` is confirmed on the smooth closed fields** and is not the whole story — see M-52 |
| M-52 | **The Marching Tetrahedra ratio is `4.0` when the surface normal lies in one octant and `2.0` when it changes sign, and P-1's `2.992` is the average of the two.** Written out for a single plane of normal `n`, the crossings are `Σ\|nᵢ\|` on the three axis families, `Σ\|nᵢ+nⱼ\|` on the three face diagonals and `\|nₓ+n_y+n_z\|` on the body diagonal. With every component the same sign nothing cancels and those sum to **exactly `4·Σ\|nᵢ\|`** — so the ratio is `4.0` for *any* orientation inside one octant, which is why a plane at four different orientations measured `3.919 / 3.939 / 3.945 / 3.943`. Across a sign change the diagonal terms cancel to `2.0`; measured `1.980 / 2.265 / 2.267`. Integrating over the sphere gives **2.9916**, reproducing P-1 to four figures | A-003, O-15. **This explains the whole reference-field spread with no new mechanism:** `box_exact`'s faces are axis-aligned one-octant normals (3.91), a sphere samples every octant (3.04), `gyroid` sits just below the isotropic average because its normals favour the cancelling ones (2.87). Two earlier hypotheses of mine — orientation, then curvature — were tested and killed first; the second failure is what forced doing the algebra instead of guessing a third time |
| M-53 | **The five algorithms fill three of the four corners of manifold × intersection-free, and Marching Cubes is the only one in the good one.** Seven reference fields, two grids, one process, one run: `marching_cubes` and `marching_cubes+decider` **0 non-manifold edges and 0 self-intersections**; `marching_tetrahedra` **0 non-manifold but 3.405 per 1k** on `csg_difference`; `surface_nets` **128 non-manifold and 0 self-intersections**; `dual_contouring` **128 and 13.837 per 1k** on `fbm_terrain`. So each of the three non-Marching-Cubes methods fails exactly one property or both, and the method the folklore treats as the crude baseline is the only one that fails neither | M-001. Cross-checks that the numbers are the same numbers: `dual_contouring`'s 13.837 and 3.118 reproduce **M-28**'s clamped `fbm_terrain` 13.84 and `gyroid` 3.12 exactly, and Surface Nets' triangle ratio of `0.977–1.001` is **✗1**'s `F_sn = F_mc + 2χ` seen from the other side |
| M-54 | **Dual Contouring is 101× more accurate than Marching Cubes on a sharp field, and indistinguishable on a smooth one.** Symmetric Hausdorff at 65³: `box_exact` **7.217e-2 → 7.145e-4 (101×)**, `thin_plate` **4.593e-2 → 5.892e-4 (77.9×)**, `csg_difference` **7.655e-2 → 2.057e-2 (3.7×)** — against `sphere` **1.2×** and `torus` **1.6×**. Marching Tetrahedra sits within 6% of Marching Cubes on the smooth fields and *better* on the sharp ones (`box_exact` 5.103e-2 against 7.217e-2), because its extra edge families sample the corner from more directions | M-001. M-26 measured this as a corner *gap* — 0.01 cells against 0.58 — and this is the same result as a whole-surface distance, which is the form that transfers to a field whose features are not corners. It also puts a number on the sentence the crate's pitch rests on: the sharp-feature solve is worth two orders of magnitude exactly where the features are sharp, and nothing at all where they are not |
| M-55 | **O-14 falsified: Marching Tetrahedra's accuracy penalty is 4.3%, not 86%, and it beats Surface Nets rather than losing to it.** Symmetric Hausdorff on a unit sphere at 64³: Marching Cubes **1.3798e-3**, Marching Tetrahedra **1.4386e-3** (`1.043×`, against a pre-registered `2.6e-3` and `1.86×`), Surface Nets **2.251e-3** (`1.69×`). And on the sharp fields Marching Tetrahedra is *better* than Marching Cubes — `box_exact` **5.103e-2** against **7.217e-2** — because its extra edge families sample a corner from more directions | M-001b. **The prediction was registered before the measurement and is wrong in its most interesting clause**: "more vertices and worse accuracy" was flagged as the counterintuitive part, and the accuracy half does not hold. Lewiner et al. 2003's underlying claim survives in direction and not in magnitude — see M-51 |
| M-56 | **Greedy meshing's `2.76×` saving over face culling is a property of one scene, not of the algorithm: measured `1.70×` to `256×`.** Same occupancy, merge on against merge off, seven reference fields at 33³: `gyroid` **1.70×**, `sphere` **1.94×**, `torus` **2.69×**, `fbm_terrain` **4.60×**, `csg_difference` **10.64×**, `box_exact` **256×**. Merging pays for flat runs, so a grid-aligned box collapses to **six quads at every resolution** — 12 triangles at 17³, 33³ and 65³ alike — while a sphere's staircase surface barely merges at all. The published figure (tier R, the UE5 benchmark) happens to sit beside `torus` | A-005. **Predicted before running** that it would not reproduce as a constant and that `box_exact` would collapse while `sphere` would not, for exactly this reason. Against Marching Cubes the blocky path costs `0.004×` the triangles on `box_exact` and `0.32–0.58×` elsewhere, which is the budget end of the tradeoff table with numbers on it |
| M-57 | **Greedy merging manufactures T-junctions, and no weld can remove them.** A blocky mesh carries split vertices on purpose — a cube corner has three faces at three normals — so its index buffer describes an open surface: a merged quad contributes five edges of which four are unshared. A-013's weld fixes that wherever quads meet corner to corner, and cannot fix a **T-junction**, because where a long quad butts against several short ones the vertex they meet at does not exist on the long quad's edge and there is nothing to merge it with. Measured on `sphere` at 33³: `2568 → 848` vertices, boundary edges `2568 → 768`, so the weld closes **70%** and the remainder is T-junctions. On `box_exact`, where every face merges to one quad and no T-junction can arise, `24 → 0` and the result is closed with `χ = 2` | A-005. The box is the control that makes this a mechanism rather than an observation. It is also the concrete form of the catalog's *"no LOD, no seam story"*: the same missing-vertex problem that breaks a greedy mesh internally is what breaks it at a chunk seam |
| M-58 | **A-010's vertex splitting removes the one-vertex-per-cell pinch completely, and the ticket named a field that was never a counterexample.** Manifold Dual Contouring against plain Dual Contouring, identical grids: `gyroid` at 33³ **15 non-manifold edges and 40 non-manifold vertices → 0 and 0**, at 49³ **48 and 99 → 0 and 0**. All seven reference fields at 17³ and 33³ come out with `non_manifold_edges == 0`, `non_manifold_vertices == 0` and `χ` equal to Marching Cubes'. But `csg_difference` — the second field the ticket names as one *"plain Dual Contouring will not manage"* — measures **0 non-manifold edges under plain Dual Contouring already**, at 33³ and 49³ | A-010. Half the acceptance criterion was vacuous as written, which is why the test asserts the comparison is non-vacuous (`any_pinched`) rather than only asserting the new number. M-53's `128` is a total over every field *and* resolution, and it is `gyroid` and `fbm_terrain` that supply all of it |
| M-59 | **The dual of a manifold surface is a manifold *complex*, and an indexed triangle mesh cannot always represent it. This is a second non-manifold mechanism, unrelated to the one A-010 fixes.** Where two cells share a face carrying **two** surface segments and each cell puts both segments in the *same* cycle, the two dual edges have the same two endpoints; an index buffer has no way to keep them apart, so they collapse into one edge with four faces. Exhibited on the ✗15 fixture — the same three-sphere union at `h = 2/3`, reached independently by proptest shrinking — where Marching Cubes is clean (`χ = 0`, zero non-manifold edges, since A-015) and the dual reports **1 non-manifold edge, 2 non-manifold vertices, `χ = 1`**. The mechanism is identified by arithmetic rather than inspection: a collapse costs exactly one edge and nothing else, so **`χ_dual − χ_mc == non_manifold_edges`** must hold, and it does — `1 − 0 == 1`. Gone by 9³ and at every resolution above, where `χ_dual == χ_mc` exactly | A-010, `the_parallel_dual_edge_collapse_is_the_only_residue`. **This bounds Nielson's guarantee**, quoted in the module docs as *"always a manifold because the original MC algorithm always constructs a manifold and the dual preserves the topology"* — the dual preserves the topology of the *complex*, and the index buffer is where it is lost. Same coarse-grid shape as ✗15, and the reason the property gate is `ClosedAllowingUnresolvedTopology` rather than `Closed`. See O-16 |
| M-60 | **Only two of seven fields ever need a second vertex in a cell, and the rate *falls* with resolution — so Nielson's "about 1.3%" is a statement about the case table, not about a scene.** Extra vertices over plain Dual Contouring, by field and grid: `gyroid` **3.13% / 2.05% / 0.53%** at 17³/25³/33³, `fbm_terrain` **1.70% / 0.84% / 0.77%**, and `sphere`, `torus`, `box_exact`, `csg_difference`, `thin_plate` **exactly 0 at every resolution**. So the cost of the manifold guarantee is zero on five of seven fields and under one percent on the other two once the grid resolves them **And it costs ~5% of the run time**: median `1.046×` plain Dual Contouring over the shootout's 14 (field, grid) points, range `1.007×` (`fbm_terrain` 33³) to `1.178×` (`gyroid` 33³) — the timing column is the noisy one, so read the median rather than the ends. Triangle counts are **identical** to Surface Nets' and Dual Contouring's on every field (`0.977–1.001×` Marching Cubes), because splitting moves vertices without adding quads | A-010. The falling rate is M-15 seen from the other side — *"any feature thinner than one cell forces two sheets through it"*, so refining removes the multi-sheet cells rather than adding them — and it is the **first curve for O-10**, which asked for exactly this rate as a function of resolution. Nielson's *"typically comprise about 1.3% of all configurations"* counts configurations in the 256-case table; on a real field the answer is field-dependent and usually zero. Fourth figure in this repo to turn out field-dependent, after ✗14, M-51 and M-56 |
| M-61 | **Splitting the vertex makes self-intersection worse, not better — ✗2's report is confirmed and the natural reading of M-29 is falsified.** Self-intersections per 1,000 triangles at 33³, plain → manifold Dual Contouring: `gyroid` **3.118 → 5.669 (1.82×)**, `fbm_terrain` **13.837 → 15.434 (1.12×)**; `sphere`, `torus`, `box_exact`, `csg_difference`, `thin_plate` **0 → 0**. The two fields that get worse are exactly the two that split (M-60), so the extra intersections are caused by the splitting itself | A-010. **The prediction registered before running was the opposite**, reasoning from M-29's *"the residue is exactly A-010's problem, with nothing left unaccounted for"* that removing the shared vertex would remove the residue. M-29's attribution is right — only multi-sheet cells have a non-zero count — and the inference from it was wrong. The mechanism is Manson & Schaefer's: the cell clamp's partition argument assumes **one** vertex per cell, and two vertices in one cell is precisely the assumption being dropped. So ✗2's ODC figure (Manifold Dual Contouring 100% of models self-intersecting) reproduces in direction here, and the tradeoff is real rather than an artefact of their implementation |
| M-62 | **The `t = a + b·n³` fit had been printing `NaN` since the day the algorithm names were spelled out, and once it prints numbers it falsifies the ticket that asked for it.** `report()` in `resolution_sweep` filtered on `["mc", "sn", "dc"]` while `Extractor::NAME` had become `marching_cubes` / `surface_nets` / `dual_contouring`; every selection came back empty and `fit` divided by zero. Fixed by deriving the list from the rows. Fitted on the **committed** CSV: `marching_cubes` **a = 0.5118 ms, b = 4.7389 ns/sample, r² = 0.99976**, and `a` is **0.64%** of the largest run — so there is no meaningful fixed cost, against M-002's *"expect a large fixed cost at small grids."* Worse for the dual methods: `surface_nets` **a = −2.746 ms (r² = 0.9923)** and `dual_contouring` **a = −2.449 ms (r² = 0.9928)**. A negative intercept is not physically possible, so the two-term model does not describe them at all — their cost grows **faster than `n³`** over this range | The `report()` fix, found while adding A-010's row. Two lessons, one of them a method rule: a list that names a thing it does not derive from will drift from it, and **nothing in this repo asserts on a benchmark's stdout**, so it drifted silently for however long. The `a < 0` result is **O-11 stated by the bench itself** rather than inferred from a ratio table, which is a stronger form of the same evidence |
| M-63 | **Both papers `docs/research/` lists as "genuinely absent, blocking" are in the home-still corpus, so the acquisition lists are stale rather than the corpus thin.** `catalog-v2.md:629` names *"Transvoxel (Lengyel 2010 …); Manifold Dual Contouring"* as blocking absences, and `catalog.md:711` and `meshing-library-target.md:97` repeat the second. Both were retrieved and read this session: `10.1109/TVCG.2007.1012` as stem `dualsimp_tvcg`, and the Lengyel dissertation as `transvoxel_dissertation_lengyel2010`. Two out of two | A-010, A-011a. **Method rule, added to Part 5:** search home-still before believing a doc that says a paper is missing. The cost of not doing so is high and asymmetric — A-010's ticket named the wrong paper for its own algorithm precisely because nobody had opened the right one |
| M-64 | **A Transvoxel lateral face does not always cross the resolution boundary, and the exception is what transition cells are *for*.** Written as an assertion — "a lateral link always joins a full-resolution sub-edge to a half-resolution one" — and falsified by the case with only the midpoint sample inside: both fine sub-edges are cut, the coarse edge is not, and the lateral links **fine to fine**, capping the feature off entirely on the fine side. Correct, because the coarse neighbour has both endpoints outside and cannot represent that feature at all. The rule that does hold is sharper: a lateral link crosses the resolution boundary **iff its half-resolution edge is cut**. Over all 512 cases and every ambiguity mask: **2,080 links stitch the seam, 1,128 cap a sub-coarse feature** **And the dissertation prescribes exactly this, which was found only after the assertion failed** — §4.3: *"In the two configurations for which the sample states alternate, either inside-outside-inside or outside-inside-outside, a mesh edge is placed on the boundary edge between the lateral face and the full-resolution face, and it is thus **not connected to the half-resolution face**."* So the derivation reproduces a rule the paper states, without the rule having been read into it | A-011a. Also measured, and both needed by A-011b: the longest cycle is **12** edges and a cell yields at most **4** cycles — the same slot budget the cube needed at A-010, so `CellVertices` transfers unchanged |
| M-65 | **Central differences at the cell size cost under half a degree of normal direction, and converge at `h²`.** Analytic gradient against `CentralDifference { step: h }` on `sphere`, worst and mean angle between the two normals: **0.460° / 0.299°** at 17³, **0.121° / 0.079°** at 33³, **0.031° / 0.020°** at 65³. Successive mean ratios **3.76** and **3.92** — `h²`, which is what a central difference must be, asserted as a range rather than only printed | A-012. This is the number a game without an analytic field actually gets, because a sampled voxel buffer has nothing finer than its own spacing to difference over. Same convergence order M-12 measured for *position* error, now for direction, and measured independently of it |
| M-66 | **On a sharp field the geometry and the field disagree by an angle that does not fall with resolution.** Mean angle between area-weighted face normals and the analytic gradient: `sphere` **3.25° → 2.16° → 1.08°** and `torus` **11.65° → 6.07° → 2.45°** at 17/33/65³, both falling; `box_exact`'s mean falls **13.55° → 6.73° → 3.34°** but its **worst is 35.796° at all three resolutions, identical to six figures**. Refining a grid does not soften a corner — the disagreement there is geometry, not discretisation, where on a smooth field it is discretisation and shrinks | A-012. So "which way does the surface face" and "which way does the field increase" are different questions wherever the surface has a crease, and that is the whole reason the strategy is selectable rather than fixed. Asserted as the mechanism (box worst resolution-invariant, sphere worst falling) rather than as a pinned constant, since the constant is a property of the box's corner angle and would move with a different field | **Amended at R-006 (M-283), and the caution at the end was right for the wrong reason.** The constant is not a property of the corner *angle*: on a wedge with one crease and no corners it does not track the dihedral at all — 5° of crease gives 87.9° — and 149 of 168 rows sit at or above 60° whatever the dihedral. What *is* confirmed is the mechanism: with the crease removed (a 180° wedge) the disagreement is **0.0000° worst and mean at every resolution**, because Marching Cubes is exact on a linear field. And the two halves of this row are one thing: the **median is 0.000° everywhere**, so the disagreement is confined to a one-dimensional crease, which dilutes out of the mean as the grid refines and leaves the worst untouched |
| M-67 | **A sign test cannot distinguish 95.6% of the configurations a tetrahedron can actually be in.** Over every edge-coordinate vector with counts 0–3 on a tet's six edges, **181** satisfy normal surface theory's two conditions (even sum, triangle inequality) and only **8** are *classic* — every edge carrying at most one crossing. The other **173** put two or more crossings on some edge, where a sign test reads the parity alone and returns the classic configuration with the same parities. And classic Marching Tetrahedra is exactly the 0/1 corner of this encoding: taking `eᵢⱼ = 1` where corner signs differ reproduces A-003's triangle count on **all 96 (tet, configuration) pairs**, as 48 corner cuts and 36 diagonal cuts | A-014a, from Baktash, Gillespie & Crane `10.48550/arXiv.2606.00454` §2. The paper's own framing is that marching tetrahedra *"reinvented a small piece of this story"* — this puts a number on how small. It is also the quantitative form of A-005's `thin_plate` result: a feature thinner than a cell does not exist to a sign-based method, and 173 of 181 is how much else does not either |
| M-68 | **`parry3d`'s constructor is not a validity check: the only mesh it refuses is one with no triangles.** `TriMesh::new` returns `Result`, and its documented failure is *"the index buffer is empty (at least one triangle is required)."* Measured: it accepts a single zero-area triangle (three collinear points) and it accepts a two-chunk mesh with an unwelded seam, both without complaint. What *does* check is `set_flags(TriMeshFlags::HALF_EDGE_TOPOLOGY)`, which builds the half-edge adjacency and returns `TopologyError` when the mesh cannot support one, and `TriMeshFlags::ORIENTED`, which needs a closed consistently-oriented surface to compute pseudo-normals from | G-005, parry3d 0.30.2. So a caller who treats "the constructor took it" as "the mesh is fine" has checked nothing, which is the gap `collider::ColliderReadiness` exists to fill. The carved acceptance case — `csg_difference` at 41³, 4,484 triangles — passes both flags |
| M-69 | **A chunk seam costs 72 boundary edges, and welding removes exactly those and nothing else.** Two adjacent 16-cell chunks of a torus, meshed independently and concatenated: **36 duplicate vertices and 180 boundary edges**. After a weld at `1e-4` cells: **0 duplicate vertices and 108 boundary edges.** The 108 that remain are the two-chunk slab's own outer border, which is legitimately open — the surface leaves through the sides — so the weld closed the 72 that were the seam and left the rest alone | G-005. A renderer draws the unwelded version correctly and a physics engine reads those 72 as a hole, which is the concrete form of G-005's ticket note *"a chunked collider must be welded first or parry sees a seam of unshared vertices."* M-46 measured the same seam at A-013 as 80 boundary edges and 40 duplicated vertices on a different chunk pair; the mechanism is the same and the count is a property of the pair |
| M-70 | **Field-derived LOD is exact, not approximate: a coarse sample position is bit-identical to the fine one it sits on.** Level `k` doubles the spacing `k` times, so a level-`k` sample at index `s` and a level-0 sample at index `2^k·s` must land on the same world point — and they do, **bit for bit**, over cell sizes `0.125`, `4/35`, `0.1` and `1/3` and levels 0–3. Doubling is exact in IEEE and so is doubling a small integer, so `(h·2^k)·s` and `h·(2^k·s)` are the same real rounded the same way | G-004. Asserted rather than argued, because **M-32 and M-49 both caught this crate assuming an algebraic identity IEEE did not honour** — `cell_of` round-trips 3 of 3 cell corners at `h = 0.125` and 1 of 3 at `h = 4/35`. This one holds at every spacing tried, including those two. It is the precondition A-011b rests on: no coordinate drift can open a crack at an LOD boundary before transition cells get a chance to be wrong |
| M-71 | **Cells fall by 8 per LOD level and triangles by 4 — and the 4 degrades exactly where the grid stops resolving the surface.** Unit sphere over a fixed world extent, levels 0–3: **262,144 / 32,768 / 4,096 / 512 cells** (exactly `8×` each, by construction) against **9,512 / 2,312 / 536 / 104 triangles**, ratios **4.114, 4.313, 5.154**. A surface is two-dimensional so its triangle count tracks `area / h²`, but that is a *continuum* claim: by level 3 the sphere is four cells across and a 104-triangle staircase is not approximating anything smoothly | G-004. So each level buys back `8×` the sampling work and only `4×` the rendering, which is the whole economics of LOD and the reason the ticket's own acceptance figure is about cells. The tight `3.8–4.6` bound is asserted only on the two steps where the premise holds, and the upward drift is asserted as a *direction* rather than waived as tolerance |
| M-72 | **A sub-cell feature does not vanish under coarsening — it aliases, which is worse.** `thin_plate` across LOD 0–3: **4,088 → 1,016 → 248 → 56** triangles, still 56 at `h = 0.5` where the plate is a fraction of a cell thick. The test was written asserting it would be *gone* by the coarsest level and that assertion is what failed. Marching Cubes samples **corners** and cuts **edges**, so whichever edges happen to straddle a thin slab still register a sign change and what comes back is a partial, holey remnant | G-004. **The contrast is the mechanism:** A-005 measured the same field returning **zero** triangles under greedy quads, which asks one question per cell *centre* and therefore misses it cleanly. For a streamed world the aliasing is the worse behaviour — a feature that vanishes at a known distance can be faded, one that disintegrates into a resolution-dependent scatter pops. It is also the cost M-67 quantified from the other side |
| M-73 | **A transition cell that computes its sample positions by offsetting from a face origin puts a hairline crack in the seam, and no weld can close it.** The first version of `TransitionCell::sample` took the face's world origin and added local offsets. At `h = 4/14` a half-resolution crossing came out at `y = -1.11e-16` where the coarse mesh had the same vertex at exactly `0` — because `(origin + h·i) + h ≠ origin + h·(i + 1)` in IEEE at a spacing that is not a power of two. Indexing from the **grid** origin instead, `origin + h·index`, which is exactly `ChunkLayout::world_of_sample`'s expression, makes the coarse `origin + (2h)·c` and the fine `origin + h·(2c)` bit-identical by M-70. Measured after the fix: **56 half-resolution crossings over 256 transition faces on a sphere at `h = 1/8`, and 24 on a torus at `h = 4/14`, every one matching a coarse mesh vertex exactly** | A-011b. **Third time this crate has assumed an algebraic identity IEEE does not honour** — M-32 at a chunk seam, M-49 in `cell_of`, now here — and the first where the consequence is a visible hole rather than a classification wobble. ~~A weld cannot rescue it: welding merges vertices it can see are the same, and these two differ in the last bit.~~ **That clause is ✗18, falsified at R-004 (M-278):** the welder's rule is first fit within `epsilon_for(h) = h · 1e-4`, not bit-identity, and it closes every one of these — 0 seam-plane boundary edges in all twenty offset rows. The headline claim stands and the reason changes: what the offset arithmetic costs is **sharing by construction**, which is what an unwelded consumer gets (M-69's collider), plus a tail where an ulp flips a sign and opens a hole 1.05–2.08 cells wide. The test searches an actual coarse vertex buffer rather than re-deriving the interpolation, because two copies of a formula agreeing proves only that they were written on the same day |
| M-74 | **A zero-width transition cell stitches the hole and has no normal at all — which is what "severe shading problems" means, precisely.** All nine of a transition cell's samples lie in the transition *face*, so at zero width every crossing does too and every triangle it emits is **coplanar with that face**. Measured: the worst `\|cos\|` between a patch face's normal and the field's gradient is **exactly 0** over 136 faces — the patch stands perpendicular to the surface it is stitching. A winding test against the gradient therefore cannot mean anything, which is how this was found: it reported 136 of 136 faces inward, and reversing the fan reported the same 136 | A-011b. Lengyel 2010 §4.3 says a zero width *"still produce\[s\] results that seamlessly stitch multiresolution meshes together, but this width leads to **severe shading problems**"* — the stitch closes the hole and shades as a hard crease **because it is one**, a flat wall standing edge-on. So the transition width is not a polish item to defer after the gap closes: **it is what gives the patch a normal**, and it is a hard dependency of any correct winding rather than a nicety. A-011c owns it |
| M-75 | **The transition width is what makes the patch's winding a measurable question at all, and the answer is unanimous.** With Lengyel's `w = 2^(k−2)` — half the adjacent full-resolution cell — the patch stops being coplanar with the transition face and becomes a ribbon: best `\|cos\|` against the surface normal **1.000**, where M-74 measured **0.000** at zero width. The orientation is then decided by measurement rather than convention: one fan order faces away from the solid on **144 of 144** faces and the other on **none**. Also measured, and the property the seam depends on: the width displaces the half-resolution face **and nothing else** — over 256 transition faces, **88 fine crossings unmoved bit-for-bit and 56 coarse crossings displaced by exactly the width along the face normal**, with no in-plane movement | A-011b. Two winding attempts were made before this one and both were meaningless, reporting the identical count in either direction, because the quantity being tested was `dot(face_normal, gradient)` on faces exactly perpendicular to the gradient. **The lesson is the method rule, not the sign:** a test that returns the same answer when you invert the thing it is testing is not measuring that thing |
| M-76 | **Two blocks at differing resolution leave 88 unmatched boundary edges in the seam plane, and transition cells take it to 0.** A full-resolution block over `x ∈ [−2, 0]` at `h = 1/8` against a half-resolution one over `x ∈ [0, 2]` at `2h`, both meshed with Marching Cubes and welded: the fine side ends on a contour of 32×32 sub-squares, the coarse on 16×16, and **88 boundary edges lying wholly in the seam plane** are the crack between them. Adding one transition cell per coarse cell face — **28 of the 256 were cut** — takes that count to **exactly 0** | A-011b's acceptance. Asserted in both directions: the `before > 0` half is what stops the test passing vacuously if the two resolutions ever stopped disagreeing. Counted *in the seam plane only*, because both blocks are legitimately open at their outer borders — a global boundary count says nothing about a seam, which is the same distinction `collider::ColliderReadiness` draws for a chunk |
| M-77 | **Lengyel's Equation 4.2 loses its level index entirely when written in the block's own cells, and the seam holds at a real width.** The published form is `Δx = (1 − 2^−k·x)·w(k)` for `x < 2^k`, zero in the middle, `(s − 1 − 2^−k·x)·w(k)` beyond `2^k(s−1)` — in **level-0** cell units, which is where the `2^−k` comes from. Substituting `x = c·2^k` with `c` the coordinate in *this block's own* cells, `k` cancels: `Δ = (1 − c)·w` for `c < 1`, `0` for `1 ≤ c ≤ s − 1`, `(s − 1 − c)·w` beyond. A linear taper across the first and last cell and nothing between. Measured with `w = 2^(k−2)`: **0 gaps at the fine plane, 0 at the coarse plane**, and the patch's best `\|cos\|` against the surface normal is **1.000** where M-74 measured **0.000** at zero width | A-011c. The cancellation is what makes it a post-pass over a finished `MeshBuffer` rather than something extraction has to know about — a block does not know which of its neighbours are coarser when it is meshed, and that can change while it is resident. Lengyel stores two positions per boundary vertex for exactly that reason; here the primary mesh is the un-inset one and the secondary is a tapered copy. **A vertex at `c = 0` moves by exactly `w`**, which is precisely the displacement `TransitionCell` gives its half-resolution face, and that coincidence is the whole reason the seam survives the width |
| M-78 | **A `no_std` crate cannot take a millisecond budget, and the honest API is a predicate the caller owns.** G-006 asked for `mesh_within_budget(ms)`. `core` has no `Instant`, so the only ways to accept a duration are a `std` feature — two paths, which this crate does not do — or a number the crate cannot compare against anything. The shipped signature takes `spend: FnMut() -> bool`, and a caller with a clock writes `\|\| start.elapsed() < budget`. **And the predicate is consulted *after* each chunk, never before:** checking first means a budget too small for one chunk drains nothing, forever, while the queue grows — a livelock that presents as a memory leak. Overshooting by at most one chunk is the price, and it is the right way round | G-006. The `no_std` constraint is doing real design work here rather than being an inconvenience: "how long does the queue take" is the question `mesh_dirty` answers and the wrong one, and being unable to phrase the ticket's version forced the API that matches what a frame actually does. Ordering is nearest-first by *squared* distance — same order, and the square root is not free per chunk per frame — with ties broken by `ChunkId` so the schedule is a pure function of the set's contents and the camera, never of insertion order |
| M-79 | **Subgrid MT's conformity is locality plus a shared vertex ordering — it is emphatically *not* invariance under relabelling, and the difference is testable.** §3.1 promises *"identical curves on triangles shared by neighboring tetrahedra"*, and the obvious reading — that the construction is symmetric in a face's corners — is **false**. Step 3(b) skips *"the first residual point along oriented edge `ij` (assuming a canonical orientation `i < j`)"*, so swapping `i` and `j` skips the point at the other end: measured on `e = (3, 0, 0)`, one labelling pairs crossings `{1, 2}` and the swapped one pairs `{0, 1}`. What is actually true is **locality**: a face's segments are a function of that face's own three edge coordinates and nothing else, verified over 4,096 configurations against all 64 fillings of the three edges a neighbouring tet would differ on. Two tets share the face's *global* vertex indices, so they agree on `i < j` without communicating | A-014b. The first version of the test asserted the symmetry and failed, which is how this was found. **The construction depends on vertex order by design and that is safe** — but only because the ordering is global. A mesh that renumbered vertices per tet would crack along every shared face, and nothing in the paper's phrasing warns you |
| M-105 | **E-203's first run reported 439 holes and a verdict of "G-001's overlap is wrong" — and the bug was in the test, one operator wide.** A probe misses only counts as a hole once **every** chunk layer that could hold the surface has been meshed; the guard said `||` where it must say `&&`, so a column whose surface sat in the one layer still in the meshing queue reported a gap. Diagnosed by printing the field at the offending points rather than by reasoning: `f(y = 0) = +0.42`, `f(y = −4) = −3.58` put the surface at `y ≈ −0.42`, in layer −1, and `colliders[−1] = false` — that layer had not meshed yet. Corrected, the same run gives **495 seam crossings, 0 holes** | E-203. **The rule this instance argues for is already in CLAUDE.md** — *"when investigating whether an issue is fixed, actually inspect the underlying data/code first; do not assume a file is broken"* — and it applies with more force when the test is *newer than the code it accuses*. G-001 has M-32, M-46 and M-69 behind it; a two-hour-old probe sweep does not, and the prior should be weighted accordingly |
| V-22 | *(settled at E-206b: adopted as a `bevy_isomesh` **dev-dependency only**, `avian3d` 0.7.0, which compiles clean against Bevy 0.19. **It costs a duplicate `glam`** — Avian pulls `parry3d` 0.27, which brings `glam` 0.33 alongside Bevy's pinned 0.32, exactly the silent double-compile `CLAUDE.md` warns about. Confined to `[dev-dependencies]`, so no consumer inherits it.)* **Avian has not been evaluated for this project — it appears once in the research, incidentally, as evidence about *release cadence* rather than as a physics choice.** `docs/research/2026-08-11-meshing-crate-architecture.md:133` cites `avian3d` 0.7.0 shipping the day after Bevy 0.19 to argue that a Bevy wrapper rides Bevy's train, and its comparison table lists **`bevy_rapier3d` (55,916)** as the Bevy physics data point. Measured 2026-08-13: **`avian3d` 0.7.0 at 127,279** 90-day downloads against **`bevy_rapier3d` 0.36.0 at 56,538** — so the doc's table names the *less* used of the two by a factor of 2.25, and `parry3d` itself is 596,889, an order above both | Queried from crates.io, and the repository the user cited (`github.com/avianphysics/avian`) is the one crates.io lists — verified rather than recalled. **What does not change:** Avian is built on parry, so the geometry contract is the one G-005 already encodes — welded, manifold, correctly wound — and M-69 measured its cost at 72 boundary edges per seam. **What it would change:** E-204's *"debris becomes rigid bodies"*, E-207 and E-209 need simulation, which a raycast cannot provide, and E-206b needs a capsule that slides. **Where it must not go:** `crates/isomesh` (rule 2), and arguably not `bevy_isomesh`'s `[dependencies]` either — the plugin stops at `Handle<Mesh>` precisely so a consumer picks their own renderer, and picking their physics engine for them is the same mistake |
| M-129 | **The seam counter's own exclusion list was missing an axis, and it accused Marching Cubes of a defect it does not have.** The two-chunk block has six walls and a boundary edge on any of them is the block ending rather than a seam failing — but the first version excluded only `x` and `z`. On `blobs`, whose spheres extend past the chunk's `y` range, the surface clipped at the top and one of those edges landed in the seam plane: **`marching_cubes` reported 1 open edge** on a seam it closes perfectly. Adding `y` takes it to 0 on both fields | B-006. **Ninth instance this session, and the second time this exact omission has been made** — E-205's crack counter needed the same fix for the same reason four commits earlier, and I wrote that one too. The pinning test carries the note so the next person adding an axis-aligned exclusion does not make it a third time |
| M-141 | **Include-once and cycle detection interact, and the obvious ordering silently accepts a circular graph.** GPU-002's preprocessor pastes each module at most once, because WGSL has no forward declarations and a duplicated function is a hard error — so "two modules both include the shared header" has to work. Written the natural way, the include site checks *have I already included this* before *am I currently inside this*, and then `a` includes `b` includes `a` finds `a` already included, skips, and **reports success** on a genuinely circular graph. Swapping the two checks fixes it | GPU-002. Found by a test written from the module docs' own claim that *"a cycle is still an error rather than being absorbed by that rule"* — the sentence was written first, the code did not implement it, and the test that asserted the sentence is what noticed. **Second failure in the same commit from the same cause:** a bare `#ifdef` with no symbol fell through the `"#ifdef "`-with-a-space matcher and was emitted as *text*, so the error surfaced one line later at the orphaned `#endif`. Splitting `#keyword` from its argument without requiring the space makes a malformed use of a known keyword an error where it occurs, while leaving an unrecognised `#word` to pass through untouched |
| M-142 | **GPU and CPU Marching Cubes agree on every triangle and disagree on 6% of vertices by exactly one ULP — and the cause is float contraction, not the algorithm.** Same 33³ grid, same `f32` samples, and the *same case table* (the shader's is uploaded from `isomesh`'s `CASES` rather than transcribed). Triangle counts match exactly on `sphere`, `torus` and `box_exact`. Of 6,936 emitted vertices: **6,507 bit-identical to a CPU vertex, 429 within one ULP per axis, 0 further** — worst separation `5.96e-8`, which is `2⁻²⁴`, one ULP at magnitude 1 | GPU-004, RTX 3090 / Vulkan. **Bit-equality was the expectation and it is wrong for a reason worth writing down.** Both sides evaluate the same expressions — `t = a / (a − b)`, then `lo + (hi − lo)·t` over `origin + h·index` — so the arithmetic is identical on identical inputs. WGSL *permits* a multiply-add to be contracted into a fused one, and this driver takes that permission: an FMA rounds once where the CPU rounds twice. **This is the divergence GPU-005's ticket anticipates**, measured before that ticket starts rather than discovered during it, and the assertion is the bound rather than a tolerance — two ULPs would fail the test. Diagnosed by measuring the *magnitude* of the miss rather than by reasoning about it: 429 mismatches could have been a formula error, and 1 ULP could not |
| M-143 | **`GridParams::sample_position` was a `mul_add` where `isomesh` is `origin + h·i`, and a power-of-two cell size hid it through an entire GPU test suite.** The GPU crate's own sample-position helper decides where the field is evaluated *before upload*, so it must be the same expression `marching_cubes::corner_position` evaluates — a multiply and then an add, rounding twice. `mul_add` rounds once. At `h = 0.125` the two forms are bit-identical because `h·i` is exact, so GPU-004's whole suite passed at that spacing and reported 0 vertices further than one ULP. E-301 ran at `h = 0.1` and **8,951 vertices came back further than a ULP from any CPU vertex** — the GPU was reading a field sampled at different points | GPU-005. **Fifth instance in this repository of an algebraic identity IEEE does not honour**, after M-32, M-49, M-70 and M-73, and the first where the wrong form was *chosen* rather than assumed. The guard is a test at `h = 0.1` that asserts the expression **and asserts that the spacing can tell the two forms apart** — at a power of two it cannot, and a fixture that cannot fail is M-44's rule wearing a new hat |
| M-144 | **GPU/CPU bit-identity is a property of the cell size, not of the port — 93.8% at `h = 0.125`, 1.1% at `h = 0.1`, 98.0% at `h = 0.0625` — and what holds at every spacing is a *distance*, under `2.1e-6` cells.** Same sphere, same samples, same uploaded case table, RTX 3090 over Vulkan through Bevy's own device | GPU-005, `cargo run --example gpu_compute_mc`. Full table, and **triangle counts are equal at every row with zero vertices moved**: `33³/0.125` **6,507 of 6,936 bit-identical**, worst offset `4.77e-7` cells; `41³/0.1` **121 of 11,112**, worst `2.02e-6`; `65³/0.0625` **27,954 of 28,536**, worst `9.54e-7`. **So GPU-005's "bit-identical, or document the divergence" resolves to the second, and the reason is not the algorithm.** Both sides evaluate the same expressions on the same samples; WGSL permits a multiply-add to be contracted into a fused one and this driver takes it, which is invisible wherever `h·i` is exact and near-universal wherever it is not. Reporting the bit-identical *fraction* without the spacing would be reporting a property of the fixture. ~~**And the honest timing, which GPU-006 owns: the GPU is slower here.** `11.54 / 11.95 / 14.36 ms` against the CPU's `0.46 / 0.88 / 2.82 ms` — 5× to 25× *behind*.~~ **Falsified at GPU-006 (M-145): those were cold numbers and the conclusion drawn from them was wrong.** Each was the *first* extraction in its process, so one call absorbed shader compilation, the first submit and driver initialisation — measured at **10.76 ms landing on a single `counts_readback`**, against `0.04 ms` for the same call warm. Warmed and taken as a median of three, 65³ is **1.95 ms against the CPU's 2.36**, i.e. the GPU *ahead*. The instrument was a single cold shot; the lesson is the one this repo keeps relearning about single samples |
| M-145 | **The GPU extraction itself is essentially free and almost flat — `0.05 ms` at 17³ and `0.13 ms` at 129³, over 420× the cells — and everything the GPU path costs is data movement around it.** Warmed, median of three, sphere, RTX 3090 over Vulkan through Bevy's device | GPU-006, `docs/measurements/gpu_vs_cpu.csv`. Total against a single-threaded CPU: **17³ 0.22 vs 0.05 ms (4.4× behind), 33³ 0.47 vs 0.33 (1.4× behind), 49³ 1.01 vs 1.04 (level), 65³ 2.09 vs 2.33, 97³ 6.38 vs 8.28, 129³ 15.03 vs 19.02 (1.27× ahead)**. **So GPU-006's expectation was right and the crossover is at about 40³** — the gap closes at small grids and inverts. **But the interesting number is the other one:** with the mesh left on the GPU, `count + emit` is **0.13 ms against 19.02**, and the ratio runs 1.04× → **0.01×**. The extraction is ~150× faster; the whole of the loss is getting data in and out. **The ceiling is the upload, and it is architectural.** At 129³ upload is **8.61 of 15.03 ms (57%)** — and `FieldBuffer::sampled` evaluates the SDF **on the CPU**, so this design does not remove field evaluation from the CPU's budget, it adds a copy to it. M-136 measured field evaluation at 65–74% of the whole job on `fbm_terrain`: on exactly the workload a GPU would help most, this design helps least. Evaluating the field in the shader is what changes that. Read-back is the rest, falling from 72% to 43% of the GPU path as the grid grows |
| M-146 | **Mesh shaders are advertised by this adapter and cannot be enabled from this workspace, because `ExperimentalFeatures::enabled()` is `unsafe` and the workspace forbids `unsafe`.** Probe over every adapter on every backend: **NVIDIA GeForce RTX 3090, Vulkan, DiscreteGpu — `EXPERIMENTAL_MESH_SHADER` advertised, `_MULTIVIEW` advertised, `_POINTS` advertised**. But requesting any of them at device creation needs an `ExperimentalFeatures` token, and its constructor is a `const unsafe fn` in `wgpu-types` 29.0.4 `src/tokens.rs`, with `disabled()` documented as *"uses of `Features` prefixed with EXPERIMENTAL are disallowed"* | GPU-007, `cargo run -p isomesh-gpu --example mesh_shader_probe`. ~~**This makes GPU-008 a policy decision rather than a task** — `unsafe_code = "forbid"` is set workspace-wide and in `bevy_isomesh`, so nothing in this repository can open a mesh-shader device without that being changed deliberately.~~ **Falsified within the hour by M-147: the conclusion was drawn from this crate's constraint without checking the consumer that actually needs it.** `isomesh-gpu` never opens a device — its API takes `&wgpu::Device`, which is GPU-001's rule — so the token is not its problem at all. Bevy writes that `unsafe` itself and already ships it, and E-303 is a Bevy example. What genuinely needs `unsafe` is only `headless::Gpu` opening its *own* mesh-shader device, which is a test convenience and not the ticket. **The first version of the probe reported `usable: false` for this adapter and it was measuring its own configuration**: it requested the feature while passing `ExperimentalFeatures::disabled()`, then reported the refusal as a property of the hardware. The field was removed rather than fixed, because without `unsafe` there is nothing honest for it to say |
| V-23 | **`CLAUDE.md`'s Metal contradiction is resolved, and both sources were right about different layers.** It records that *"wgpu's spec table lists MSL as planned while the tracking issue says the Metal HAL backend merged"*. `wgpu-types` 29.0.4's own source settles it on `EXPERIMENTAL_MESH_SHADER`: *"Supported platforms: Vulkan (with `VK_EXT_mesh_shader`), DX12, Metal"* followed by **"Naga is only supported on vulkan. On other platforms you will have to use passthrough shaders."** | Read from `wgpu-types-29.0.4/src/features.rs`, GPU-007. So the **feature** reaches Metal and the **WGSL compiler** does not: on Metal a caller hands wgpu pre-compiled MSL rather than the WGSL this crate composes. Mesh shaders are therefore a *fork in the shader pipeline*, not a flag on it, which is the fact GPU-008 has to be designed around. **Metal remains unmeasured here** — this machine is Linux with one Vulkan adapter, and the probe says so rather than inferring. **naga 29 does implement the WGSL side on Vulkan**, verified by parsing rather than by reading: `enable wgpu_mesh_shader;` is accepted as *implemented* (not merely named), `@task` parses, and `@mesh` requires `@mesh(<global>)` naming an output variable in the **`workgroup`** address space (`front/wgsl/parse/mod.rs:1922`, `valid/interface.rs:1531`) — which is how the first test failed, having written a bare `@mesh` and a `dispatchMeshWorkgroups` call that does not exist; the task dispatch is a builtin **output**, `@builtin(mesh_task_size)` |
| M-147 | **Bevy's device already has mesh shaders enabled, with no configuration and no `unsafe` anywhere in this repository — so GPU-008 was never blocked, and the block was my inference rather than a measurement.** Measured on Bevy 0.19's `RenderDevice` in a running example: **`mesh_shader=true multiview=true points=true`** | GPU-007's correction, RTX 3090 / Vulkan. Three facts compose to it, each read from source. `WgpuSettings`' default priority is **`Functionality`** (`bevy_render-0.19.0/src/settings.rs:89`), which sets `features = adapter.features()` — *every* feature the adapter advertises. `features \|= options.features` then honours anything a caller adds explicitly (`renderer/mod.rs:318`). And Bevy passes the experimental token itself: **`experimental_features: unsafe { wgpu::ExperimentalFeatures::enabled() }`** at `renderer/mod.rs:335`, carrying its own `SAFETY: TODO` and a tracking issue. **So the crates here stay 100% safe Rust and still reach mesh shaders**, because `isomesh-gpu` takes `&wgpu::Device` and never creates one — GPU-001's API rule paying a dividend it was not designed for. **The error worth recording is the reasoning, not the fact.** M-146 measured a real constraint on `isomesh-gpu`'s *own* headless device and generalised it to "nothing in this repository can open a mesh-shader device", without checking the path E-303 actually takes. A blocker asserted from one crate's configuration is a claim about the consumer, and it needed the consumer to be looked at. **It is also worth knowing that the `unsafe` guards nothing**: `ExperimentalFeatures::enabled()`'s entire body is `Self { enabled: true }` — a struct literal with a `bool`. The `unsafe` is a hand-raise token acknowledging that experimental *features* may carry UB, not a memory operation | **Amended on independent review, four qualifications, three of which narrow the claim.** **(1) The free device is contingent on an unresolved upstream question.** Bevy's comment at that line is `// SAFETY: TODO, see https://github.com/bevyengine/bevy/issues/22082` — an *admission that a justification is owed*, with an open issue, not a justification. Describing it as "its own SAFETY note" read as more settled than it is. If #22082 lands as opt-in, the default path loses mesh shaders. **(2) "E-303 gets a device for free" holds on this machine, with default settings, today — and is one of three branches.** `WgpuSettings::default()` calls `settings_priority_from_env()` first, so **`WGPU_SETTINGS_PRIO` overrides the default** and consumers do set it (wasm); under any priority other than `Functionality`, `features` starts at `wgpu::Features::empty()` and gains only `options.features` (`renderer/mod.rs:296`); and `adapter.features()` is machine-dependent, measured here on one RTX 3090. **So GPU-008's runtime capability probe is load-bearing rather than belt-and-braces, and the ticket was right to specify it** — the correction over-generalised in the opposite direction from M-146. **(3) "The `unsafe` guards nothing" undersells it.** A `Self { enabled: true }` body does not make the contract vacuous; it makes it **compiler-unenforceable**, which is exactly what an `unsafe fn` used as a hand-raise token is. wgpu marked it unsafe because experimental features can produce driver UB and someone must accept that. *"There is no risk"* and *"we did not take the risk"* are different claims and only the second is true. **(4) The headless gap is a coverage hole, not a convenience** — reclassified and ticketed as GPU-009 |
| M-148 | **A WGSL mesh shader's output block is a fixed contract that naga defines and no document states — and it is derivable from the source rather than guessable, which is the difference between writing one and inventing one.** The mesh output must be a `var<workgroup>` **struct** whose members carry exactly four builtins: `@builtin(vertex_count) u32`, `@builtin(primitive_count) u32`, `@builtin(vertices) array<V, N>` and `@builtin(primitives) array<P, M>`, both arrays **constant**-sized (`ArraySize::Dynamic` is rejected), with the topology inferred from the primitive struct — `@builtin(triangle_indices)` → triangles, `point_index` → points, `line_indices` → lines | GPU-008a, read from `naga-29.0.4` `proc/mod.rs:723` and `valid/interface.rs:1531`, then confirmed by validating. Derived this way the first candidate **parsed and validated on the first attempt**, where three earlier guesses at the syntax had each failed. **And one real language restriction fell out of it:** a function parameter may not be a pointer into the `storage` address space — naga reports `InvalidArgumentPointerSpace` — so a shared `read_vec3(ptr, i)` helper is impossible and two near-identical readers are the language's requirement rather than a preference. **The consequence for coverage is the useful part:** the shader validates with **no GPU**, so GPU-003's naga sweep covers it in CI even though GPU-009 records that it cannot be *executed* headlessly. Validity is covered; behaviour is not; and those are now separable claims |
| M-149 | **A mesh-shader draw removes the *smallest* of the GPU path's three data-movement costs — 6.7% at 129³ — and the doc comment claiming otherwise was written before the arithmetic was done.** Breaking `gpu_vs_cpu.csv`'s 129³ row down: **upload 8.63 ms (57%), CPU prefix sum 3.27 ms (22%), counts read-back 1.97 ms (13%), geometry read-back 1.00 ms (6.7%)** — and the geometry read-back is the only one a mesh-shader render path eliminates | GPU-008a, from `docs/measurements/gpu_vs_cpu.csv`. **The share falls with resolution, from 40.5% at 17³ to 6.7% at 129³**, because it scales with triangle count while upload scales with cells. So the saving is largest exactly where the whole path is cheapest and least worth optimising. **The counts read-back is bigger than the geometry read-back** (1.97 against 1.00) — 4 bytes per *cell* beats 36 bytes per *triangle* at any resolution where the surface is sparse in the volume, which is all of them. **What would actually remove the movement** is a GPU prefix scan feeding `draw_mesh_tasks_indirect` (verified present in `wgpu-29.0.4` `api/render_pass.rs:303`): the counts never come home, the total never comes home, and the draw reads its own workgroup count from GPU memory. Ticketed as GPU-010. **The error being recorded is an overclaim in prose that the repo's own committed CSV contradicted** — `mesh_render.wgsl` said "there is no read-back" and "removes the largest remaining piece" while `docs/measurements/gpu_vs_cpu.csv` had the numbers to disprove it. Checking a claim against data already in the repository costs one query, and this one was not made until the ticket after |
| M-150 | **Moving the prefix sum onto the GPU took the extraction path from 15.01 ms to 9.65 ms at 129³ — 1.56× — and the stage it replaced went from 5.24 ms to 0.37 ms, a factor of 14.** Predicted 34.9% from M-149's breakdown, **measured 35.7%** | GPU-010a, `docs/measurements/gpu_vs_cpu.csv`, RTX 3090 / Vulkan, warmed, median of three. Against a single-threaded CPU the whole path goes from **0.79× to 0.51×** — 1.27× ahead to **1.95× ahead** — and every row from 49³ up improves. **What it removed was not computation, it was a copy:** the CPU prefix sum needed every per-cell count, which is 4 bytes × 2,097,152 cells = **8.4 MB** read home at 129³ purely to add up. The scan leaves **four bytes** crossing the bus, the grand total, needed only to size the geometry buffer. **It is a small loss at tiny grids and that is worth saying: 17³ went 0.24 → 0.28 ms (0.85×)**, because two extra dispatches cost more than summing four thousand numbers on the CPU. The crossover is around 25³, well below the 40³ where the GPU path as a whole starts winning, so it never has to be chosen against. **The consequence for what comes next: upload is now 87% of the GPU path** (was 57%), and upload is `FieldBuffer::sampled` evaluating the SDF **on the CPU** and copying samples over. Every remaining stage together is 13%. Evaluating the field in the shader is the largest lever left, and it is not ticketed | **Refined at M-152: "the only lever" was too strong.** Only **57%** of that upload is `field.sample()`; the other 43% is a `Vec<f32>` → `Vec<u8>` conversion and a `write_buffer` copy, which is a second, independent lever and the only one that can help a consumer who already has samples |
| M-151 | **A hierarchical GPU scan is exactly the kind of code that is wrong invisibly, so it is asserted element-for-element against a CPU reference at all three depths — and the depth is asserted too.** `ScanOutput::levels` is exposed for that reason: a suite that only ran inputs under 65,536 elements would exercise one and two levels and never touch the cross-block add that three levels need, while reporting full coverage | GPU-010a. Tested: 200 elements (1 level), exactly 256, 1000 (2 levels), 300,000 (**3 levels**, the depth a 129³ extraction uses), all-zeros, all-ones, single element, and a **sparse clustered distribution** shaped like real Marching Cubes counts — zero almost everywhere with short runs where a surface crosses — because a scan tested only on smooth data has not met its actual input. Every case compared offset-by-offset against `cpu_prefix_sum`, not by total. **A total-only comparison would pass on a scan that shuffled the offsets**, and shuffled offsets produce triangles written to each other's slots: a mesh that renders, has the right triangle count, and is silently corrupt. The end-to-end evidence is separate and stronger — the pre-existing CPU-comparison tests still report identical triangle counts and vertices within one ULP with the scan in the path |
| M-152 | **"The upload is field evaluation" is 57% true: of the 8.40 ms upload at 129³, evaluation is 2.65 ms, converting `Vec<f32>` to `Vec<u8>` is 1.14 ms and `write_buffer` is 0.86 ms.** Timed separately, release build, warmed, RTX 3090 host | GPU-011a/GPU-012's premise, checked before either was written. **So there are two independent levers, not one**, and they serve different consumers: evaluating the field in a shader removes the whole upload for a caller with an analytic SDF and can do nothing for a caller who *already has samples* — a voxel volume or an imported scan — whose upload is 100% copy. Writing samples straight into a mapped buffer helps that caller and is worth ~43% of the upload either way. **Two smaller corrections fell out of the same measurement.** `create_buffer` for 8.4 MB measures **0.00 ms** — wgpu allocates lazily, so "the upload allocates a fresh buffer every call" is not a cost worth avoiding. And the first run of this measurement reported `eval 98.76 ms`, which is a **debug build**: `cargo test` does not inherit `--release`, and a timing taken that way is 37× off and looks like a catastrophic finding rather than a mistake. The repo already knows this — CLAUDE.md says *"always `--release` for examples"* — and the same rule needed applying to a test used as an instrument |
| M-153 | **Eliminating the "redundant" copy into GPU memory costs 1.6×, and the copy that *was* redundant was a different one.** GPU-012 predicted 43% of the upload recoverable by writing samples straight into a `mapped_at_creation` buffer. Three variants, 129³ upload: **`Vec<u8>` + `write_buffer` 8.40 ms** (the code that was already there), **`mapped_at_creation` + per-element `write_iter` 13.47 ms**, **`mapped_at_creation` + bulk `copy_from_slice` 8.62 ms** | GPU-012. **So the mapping is not the problem and per-element writes are.** wgpu says why in a comment on the type itself: `BufferViewMut` deliberately does *not* deref to `[u8]` because *"mapped memory may be write-combining memory"* — where a stream of small writes is far worse than one memcpy. A bulk copy into a mapping merely ties `write_buffer`, so `write_buffer` was already doing the efficient thing and there was nothing to reclaim there at all. **What did work was the copy nobody had named:** `sampled` built a `Vec<f32>` and then converted it to bytes in a second full pass over 8.4 MB. Evaluating straight into the byte layout takes the upload to **7.04 ms** and the whole path from **9.65 to 8.15 ms — 1.18×**, and from 1.95× to **2.33× ahead** of a single-threaded CPU. **The prediction was 43% of the upload and the answer was 16%**, because half the target turned out to be an irreducible 8.4 MB memcpy. Registering the number first is what made that visible as a miss rather than a success |
| M-154 | **GPU field evaluation does not match `libm` bit-for-bit on *any* field, including the ones with no transcendentals — and it does not change the mesh.** Same grid, same expressions in the same order, `h = 0.125` from origin `-2.0` so the sample *positions* are bit-identical and only the field arithmetic differs. Bit-exact samples of 35,937: **`sphere` 26,009, `torus` 25,449, `box_exact` 30,637, `gyroid` 3,873**; worst absolute deviation **2.38e-7 / 4.77e-7 / 1.19e-7 / 6.56e-7** | GPU-011a. **The hypothesis was that sqrt-only fields would agree exactly, and it is wrong**: a GPU may contract `x*x + y*y + z*z` into fused multiply-adds, rounding once where `libm` rounds twice — M-142's cause, one layer up from the vertex positions where it was first found. **`box_exact` comes closest for a structural reason worth generalising:** its expression is `abs`/`max`/`min`, exact in IEEE, over a `sqrt` whose argument is identically zero *inside* the box, so only the exterior shell has anything to round — 1 ULP worst against `sphere`'s and `torus`'s hundreds. **So the drift is a property of the expression, not of the GPU**, which is the fact GPU-011b has to design around. **And the consequence that actually matters is measured rather than argued: the mesh is unchanged.** Extracting from a CPU-uploaded field and a GPU-evaluated one gives **identical triangle counts** on all four fields, including `gyroid` at a non-power-of-two spacing where the positions differ too — no sample within `7e-7` of zero crossed the sign test. **ULP counts are reported but not gated**, because a signed distance field is zero on its own surface and ULP distance is unbounded near a zero crossing; absolute deviation is the quantity that decides whether a crossing moves |
| M-155 | **Evaluating the field on the GPU takes the 129³ path from 8.37 ms to 0.54 ms — 15.5× — and 37× ahead of a single-threaded CPU. The prediction was "under 2 ms and ≥4×", registered before the work.** Warmed, median of three, RTX 3090 / Vulkan | GPU-011a, `docs/measurements/gpu_vs_cpu.csv`. Full sweep, GPU-native against the uploaded path: **17³ 0.22 vs 0.30, 33³ 0.22 vs 0.42, 49³ 0.26 vs 0.78, 65³ 0.27 vs 1.35, 97³ 0.33 vs 3.72, 129³ 0.54 vs 8.37**. **The shape is the finding as much as the number: it is nearly flat** — 0.22 → 0.54 ms across a **420×** increase in cells — the same not-remotely-saturated signature M-145 found in `count + emit`, now for the whole path. What was 86% of the run has become nothing, because the samples are produced where they are read and the bus is never touched. **The crossover moves to about 25³**: at 17³ the GPU is still 3.7× behind on ~0.22 ms of fixed overhead, and from 33³ up it wins, by 37× at 129³. **Cumulative across the three tickets that touched this path: 15.01 ms at GPU-010a's start, 0.54 ms now — 28×** — and none of it came from making the extractor faster. `count + emit` never moved: it was 0.11 ms before and 0.04-0.05 ms now. Every gain was data movement removed |
| M-156 | **The mesh-shader coverage hole is accepted rather than closed, and the deciding constraint is a stated goal rather than a technical one.** `headless::Gpu` cannot open a device with `EXPERIMENTAL_MESH_SHADER` without an `unsafe` block, so a mesh-shader kernel can be *validated* headlessly (GPU-003's `naga` sweep needs no device) and *refused* headlessly (`MeshShaderRenderer::new` errors, asserted) but never **run** | GPU-009. The alternative was to relax `unsafe_code = "forbid"` for one crate and confine the `unsafe` to one documented call — which would work, and whose entire body is `Self { enabled: true }`. It is not taken because the user's stated goal is a crate that can be described, accurately, as safe Rust, and one `#[allow(unsafe_code)]` costs exactly that description for a test-only convenience. **What is bought and what is given up, separated:** shader *validity* stays covered in CI; mesh-shader *behaviour* does not, and must be asserted inside an application whose device asked for the feature — Bevy's does by default (M-147). Compute kernels are entirely unaffected: Marching Cubes, the scan and the field sampler all run and are tested headlessly. **The doc is pinned by a test** asserting the device really has no experimental features, so if wgpu stops gating them or the lint policy changes, the claim fails loudly instead of quietly becoming false |
| M-157 | **A GPU-interpreted edit log reproduces `BrushStack` to under `8.4e-7` across every shape, every op and a twelve-brush mixed sequence — and extracts the same mesh.** Worst deviation from `isomesh`'s CPU fold, 25³ grid: **each op alone 1.19e-7**, **each shape alone 2.38e-7**, **a 12-brush mixed log 8.34e-7**; an empty log reproduces the base field **bit-identically**; and a three-brush log at 49³ gives **the same triangle count** as the CPU-sampled equivalent | GPU-011b. **The mixed-log case is the one that matters and it is not redundant.** Mixed adds and subtracts do not commute and a smooth add is not associative (M-36..M-38), so folding first-to-last is part of the *value*, not an implementation detail — a shader that batched or reordered brushes for parallelism would pass every single-op test and fail only this one. **The design was chosen over composing a consumer's WGSL for a reason measurement supports rather than taste:** M-154 established that arithmetic drift is a property of the expression rather than of the GPU, so an interpreter evaluating the same primitives drifts the same as compiled WGSL would — option (b) buys no accuracy, while option (a) survives editing without recompiling a shader per carve. What it costs is stated rather than hidden: the log is bounded by sphere, box and capsule, so a consumer with an arbitrary analytic SDF is not served, and that is the CAD half of the audience |
| M-158 | **The whole pipeline runs on the GPU and never touches the bus: field, edit log, extraction and draw, at 0.41 ms for 16,436 triangles at 65³, 60 fps.** E-303 evaluates a base sphere with three moving brushes folded over it (GPU-011a/b), extracts with the GPU-scanned Marching Cubes (GPU-010a), and draws straight out of the position and normal buffers with a mesh shader (GPU-008a). Per frame the CPU sends **a camera matrix and three brushes** | GPU-008b, RTX 3090 / Vulkan through Bevy. **Two integration facts cost a run each and are worth writing down.** *A render pipeline must match the pass's MSAA sample count*, and wgpu refuses the draw outright rather than degrading — *"the RenderPass uses textures with sample count 4 but the RenderPipeline uses attachments with format 1"*. Bevy defaults to 4×, so `MeshShaderRenderer::new` gained a `samples` parameter; reading it from `ViewDepthTexture::texture.sample_count()` is better than assuming, because the view owns it. And *Bevy 0.19 has no `Node` trait*: `Core3d` is a **schedule** with `Prepass → MainPass → EarlyPostProcess → PostProcess`, so a custom pass is a system in one of those sets — `EarlyPostProcess` here, which orders after the opaque pass by the chain without naming a Bevy system and depending on it staying public. **The refusal branch is forced and screenshotted rather than assumed**: `ISOMESH_NO_MESH_SHADER=1` makes the unsupported path reachable on hardware that supports the feature, because a branch nothing can execute is precisely the failure this project keeps finding, and *"never panics on an unsupported adapter"* was otherwise an untested claim |
| M-192 | **I-007's own finding happened to I-007, in the commit that filed it.** The new `msrv` leg for `bevy_isomesh` went red immediately, and not on a Rust version: `wayland-sys`'s build script panicked because `libwayland-dev` is absent from the runner. Bevy 0.19 enables Wayland by default on Linux, so its build scripts need the system packages **even to type-check** — the `bevy` job has carried an `apt-get` step for exactly this since it was written, complete with a comment recording that *"none of these exist on macOS, which is why a local verification could not have caught it"* | I-007. `cargo +1.95 check --all-targets` passes locally in that directory, so the code was never the problem; the *job* was, in the one way the excluded workspace always breaks jobs. **That is the fourth instance of the pattern I-007 filed** — after rustfmt, rustdoc and MSRV itself — and it arrived inside the fix for the third. A rule earned twice in one commit is worth more than the incident: **any new CI step that enters `bevy_isomesh` needs that directory's whole environment, not just its Rust toolchain**, and the only reliable way to know is to copy the `bevy` job's prologue rather than to reason about what a step "should" need |
| V-27 | **§3.2.3's prose is fully retrievable and its triangulation patterns are still only a picture — and the paper offers a second option the ticket never recorded.** Retrieved from the source (`10.48550/arXiv.2606.00454` §3.2.3, via home-still): the inset construction is stated in full — *"we insert the midpoints of all polygon edges contained in any edge of `f`, and move them a small distance in the inward normal direction. For pentagons, we also insert the midpoint of the edge opposite the distinguished corner."* — and is immediately followed by *"The triangulation patterns in **Figure 15, right** ensure that inset regions stay within the tetrahedron."* The patterns themselves appear nowhere in the text. **The sentence A-014d did not have** is the one after it: *"This procedure is needed only when taking the union of the two polygons would yield a nonmanifold edge. **If manifold connectivity is not required, one can also just discard one of the polygons.**"* | A-014d's rule-5 check. So there are two published routes, not one, and the paper is explicit that the cheap one costs manifold connectivity. **This crate cannot take it**: `the_validity_suite_over_every_reference_field` gates on manifoldness and `csg_difference`'s pinned `(3, 6, 6)` is precisely a manifoldness row, so discarding a polygon would trade the defect the ticket owns for a different one. The rule-5 stop therefore stands exactly where the ticket predicted — at the quad/pentagon/hexagon triangulations, which must be derived or the figure read, and **not invented** |, not in a hot path — the plugin never shades at all. What was real was the doc comment telling consumers to do it.** The ticket recorded *"`bevy_isomesh/src/mesh.rs`'s shading path allocates a scratch buffer per call to split the borrow, which at re-meshing scale is an allocation per chunk per edit."* There is no such path: `isomesh::paint::shade` appears in `bevy_isomesh/src` exactly twice, once in an `ignore`d doc example and once in a test that runs once. `grep -n "colors\|shade" src/plugin.rs` returns nothing — the plugin builds geometry and never touches colour | B-008. **What was genuinely wrong is one line of advice**: the doc said the borrow *"needs splitting in real code — read the positions into the call, or shade into a scratch buffer and swap"*, and a consumer following it pays exactly the allocation the ticket described, in their code rather than in ours. The fix is therefore an API rather than a pool: `positions_and_colors_mut` hands out both borrows at once from disjoint fields, so the natural line compiles and the pass writes straight into the array the `Mesh` will own. **A pool would have been the wrong answer even where the premise held** — the allocation is avoidable outright, and reusing one is what you do when it is not |, and nothing had checked it — the third hole the excluded workspace has produced.** The `msrv` job reads `cargo metadata --no-deps` from the repository root, and the root workspace *excludes* `bevy_isomesh` by design (feature unification), so the `rust-version = "1.95"` that landed with M-174 was derived from Bevy 0.19's floor and enforced by nothing. Verified from the resolved graph: the highest `rust_version` anywhere in `bevy_isomesh`'s dependencies is **1.95.0**, declared by `bevy` and `bevy_ecs`, with `bevy_math` and `bevy_input_focus` at 1.94.0 and the isomesh crates at 1.89 — so the claim was right, by luck of having been copied from the correct source | I-007. **The pattern is what matters, not this instance.** rustfmt, rustdoc and now MSRV have each had to be patched into the `bevy` job separately, always after the gap was found rather than before. **Every check that reads the root workspace needs a deliberate answer for `bevy_isomesh`, and "it is excluded" is not one** — the exclusion is a deliberate trade for a pristine root lockfile (CLAUDE.md), and its cost is that every gate has to be written twice or explicitly waived |
| V-28 | **A-014d's rule-5 stop is lifted, and not by the figure — by the authors' own reference implementation, shipped as arXiv ancillary material that nobody had looked for.** V-27 concluded the triangulation patterns *"appear nowhere in the text"* and left the ticket stopped on *"derive the patterns or look at the figure"*. Both halves have now been done, and a third route turned out to exist. **(a) The figure was read.** Figure 15 is on page 9; rendered with `pdftoppm` at 1200–2400 DPI it is fully legible, and the quad and hexagon patterns transcribe exactly. The figure's own colouring is the key the prose never states: **thick blue is a segment of `γ`, along which the two coincident copies stay glued — which is why two quads make an annulus and the paper can call the pair a "tube" — and thin black is a polygon edge contained in an edge of `f`**, the ones §3.2.3 inserts a midpoint into. Measured rather than eyeballed, because the pentagon panel turns on it: sampling the rendered PPM gives blue-excess **+109** on the two chord edges and **0, 0, +14** on the three `∂f` edges. The crate already carries this distinction — `Arc::Chord` against `Arc::Edge` in `subgrid/surface.rs`. **(b) The pentagon could not be read**, and it is the only type A-014d needs (M-193). Its structure is unambiguous — a fan from one interior point over the boundary ring — but which edge that point is the midpoint of cannot be told apart from the drawing, and the prose's *"insert the midpoints of **all** polygon edges contained in any edge of `f`"* contradicts a figure that inserts one. **(c) The arXiv listing has ancillary files** — `anc/SubgridSupplement.pdf` and `anc/subgrid-tet.html` — and the second is a standalone JS implementation with a **"Simplicial Embedding"** checkbox. `combFaceToPolygonSoup(combFace, tetPositions, edgeIsectTs, midScoopVertices, bulge)` is §3.2.3, from its authors, in executable form | A-014d. **What the implementation settles that neither prose nor figure did.** *Ordinary crossings never move* — `if (!cv.scoopSteiner && !cv.scoopInteriorSteiner) return pos;` — so conformity across tets is preserved by construction, which is the exact worry M-101 raised twice and the reason the inset is safe at all. An inserted point is the **midpoint of two consecutive crossings on the same tet edge**, displaced `lerp(midPos, otherEdgeMid, bulge)` toward the midpoint of the **opposite** tet edge; two tets sharing `f` have different fourth corners, so their copies move apart and each into its own tet. `INTERIOR_HEXAGON` is 9 vertices and 7 triangles, `[8,2,5] [2,8,0] [2,0,1] [5,2,3] [5,3,4] [8,5,6] [8,6,7]` — inserted at ring positions 2, 5, 8, central triangle on all three — **which reproduces the hexagon transcribed from the figure triangle-for-triangle**, two independent sources agreeing. `INTERIOR_CORNER_TYPE` inserts **one** Steiner point at ring position 1, at `lerp(mid(v0,v1), mid(commonI, opposite), bulge)` where `commonI` is the tet corner shared by the first two vertices' edges, and **fans the whole polygon from it**. So the figure was right and the prose was the loose one. **The method rule is the finding:** V-27 asked whether a paper stated something and stopped at the paper. An arXiv `abs` page lists ancillary files in the HTML; one `curl` would have found an executable answer at any point in this ticket's life |
| M-193 | **The polygon type A-014d has to inset is the one type its figure could not be read for — measured, not assumed.** §3.2.3 names three (*"a quad, a hexagon, and a pentagon made at a corner"*), and implementing a type no field produces would be unverifiable code, so `which_polygon_types_coincide_across_a_shared_face` counts them. A coincident polygon is a `Region` of one face seen by two tets; classifying by `(loop kind, chords, edges)` — `Arc::Chord` against `Arc::Edge`, the figure's blue against black — gives, at 17³ over all seven fields: **`(contractible, 1, 1)` 384, `(corner, 1, 1)` 126, `(corner, 1, 2)` 14, `(corner, 2, 3)` 5**. `sphere` and `torus` have none of any kind | A-014d. **Three things follow.** (1) `(corner, 2, 3)` is 5 arcs from a corner-type loop — **exactly Figure 15's "pentagon made at a corner"** — and it is precisely `csg_difference`'s **3** plus `thin_plate`'s **2**, which is the ticket's whole target restated from the other side. The quad and hexagon patterns, the two that transcribed cleanly, are **unreachable on every reference field**. (2) The 510 `(*, 1, 1)` entries are **bigons** — one chord, one edge — the zero-area scoop V-21 predicted. They are invisible to the triangle-keyed census that pins `csg_difference` at 3, because a 2-node region emits no triangles at all: `for k in 1..corner.len()-1` is an empty range. Two instruments, two different denominators, and neither was wrong. (3) `box_exact` carries 84 coincident bigons and **zero** of anything else, which is why M-161 called it the worst field and M-186 found it at 0 — it was never §3.2.3's immersion |
| M-205 | **A-014i's "orphaned bigons" is real, and it is a law rather than a defect: unreferenced positions arrive in pairs, one pair per bigon nothing else reaches.** A bigon is a two-arc region — one chord, one edge piece, the zero-area scoop V-21 predicted and M-193 counted 510 of — and it emits no triangle at all, since a fan over fewer than three nodes is an empty range. Where its two crossings are reached by no other region, the fill records two positions and nothing uses them. Measured per **patch**, because `cell_tet` hands every position to the sink and a neighbour may then reference it: `sphere` and `torus` have **0 bigons and 0 orphans**; the other five have both, and **every affected tetrahedron leaves exactly two** | A-014i, third defect, characterised and not fixed because there is nothing to fix. **The apparent exception is the best part.** `fbm_terrain` leaves 302 across 148 tetrahedra where twice-per-tetrahedron gives 296, which read as a second mechanism until the histogram was printed: **145 tetrahedra with 2 and three with 4**, i.e. two bigons in one tetrahedron. `145·2 + 3·4 = 302`. So the law is exact on all seven fields once stated as *pairs* rather than as *two per tetrahedron*, and the assertion is now evenness — which is ✗1's rule, assert the identity and the counterexample explains itself. **Why this is not a defect.** A bigon legitimately emits nothing: it has zero area, and A-014d closed as *not required* precisely because there is no geometry to inset into it. Its two positions carry real global crossing identities, so a neighbouring tetrahedron consumes them — which is exactly what M-201 measured, where dropping such positions changes vertex order and not vertex count. Unlike the subdivision orphans, which nothing consumed and which M-203 fixed |
| M-204 | **M-166 is closed by two parentheses, and the test that was already asserting the property could not have caught it.** The interior test's denominator was `v[0] + v[2] - v[1] - v[3]`, which evaluates as `((A + C) − B) − D` — a fixed subtraction order that a rotation **permutes**, and IEEE addition is not associative. So two cells reading one shared face could evaluate the same sweep to different bits and disagree about a tunnel, which is why A-002b was blocked on this. Grouping each diagonal, `(A + C) − (B + D)`, makes every symmetry of the face exact: a two-corner rotation swaps the operands *within* each parenthesis and IEEE addition is commutative, so the result is bit-identical; a one-corner rotation or a reflection exchanges the two parentheses and IEEE subtraction is exactly antisymmetric, so the result is the exact negation — which cancels against the numerator's own exact negation in `saddle`. Over 20,166 corner quadruples the old order disagrees on **2,764** and the grouped one on **zero** | A-002b's named prerequisite. **The existing test asserted this property and passed with the defect present**, which is the part worth keeping. `the_test_is_deterministic_and_independent_of_diagonal_order` rotates by two corners and compares saddles, and strengthening its tolerance to bit-identity **still passes** — `opposed()`'s corners are all of similar magnitude, so both orders round the same way. That is the fixture trap for the **fifth** time (M-32, M-38, M-44, G-003, here), and the fixture that exposes it had to be *searched* for: `(1, 1, 1, 10⁻⁸)`, giving `0.99999999` against `0.9999999900000001`. The new test asserts first that the fixture still separates the two orders, then that the sweep agrees bit-for-bit — and it fails on the old expression and passes on the new one, verified in both directions. **A tolerance was doing the work of a fixture search**: the original comment correctly identified the risk and deferred it, and the deferral outlived the two lines that fix it |
| M-203 | **A-014i's orphaned vertices are fixed by giving children their parent's names instead of new ones, and the mesh is bit-for-bit the same mesh.** M-201 measured 16,296 of 166,591 vertices referenced by nothing on the subdividing plane fixture, and named the cause: a subdivided parent records a position for every crossing its cycles name, and its four children then re-derive those same points as **child-local positions beyond `crossings.len()`, with no identity at all**. `fill_append` now takes an `inherited` map, and `fill_subdivision` builds one per child from the parent's own crossing indices — forward or reversed exactly as the parameters are. Result: the bare `(4, 2)` tetrahedron goes **112 positions → 64 with zero unreferenced**, and the grid fixture **166,591 → 104,767 vertices with 46 unreferenced**, a 37% reduction | A-014i. **The triangle count is identical in both — 52 and 219,084 — which is the evidence that this renamed vertices rather than changing geometry.** All 523 tests pass including the golden hashes, which is the second guard: subdivision is unreachable on the seven reference fields, so a change confined to it must leave every golden alone, and it does. **The detector for the `1 − t` re-lerp moved and had to be re-checked rather than assumed.** M-200's bare-tetrahedron test worked by observing that children *duplicate* parent crossings, and inheritance removes exactly that — `rediscovered` goes 24 → 0, which is now the invariant it asserts. Inverting the branch is still caught, but by the **grid fixture** instead, verified by re-running the mutation after the fix. Had that not been checked, the fix would have silently deleted the only test with power over a defect this ticket had just finished proving was untested. **The 46 remaining are pinned, not chased**: 0.04%, and crossings the parent named that no child's cycle reaches at all, which is a different mechanism from the systematic duplication measured here |
| M-202 | **Subdivision only fires where the grid badly under-resolves the field, so "T-001 validity on the subdivision fixture" cannot be a clean gate — the fixture is compromised by construction.** The `k = 96` plane fixture welds to 73,855 vertices with **249 non-manifold edges, 295 non-manifold vertices and 131,815 inconsistently oriented**, which reads as a subdivision defect and is not one. Holding the algorithm and the grid fixed and moving only the frequency: the same plane family at `k = 8`, which this 9³ grid resolves and which fires **no** subdivision at all, gives **0, 0, 0** — clean, with 948 boundary edges, which is a plane family leaving through the walls of a finite box and is the field rather than the mesher. Same extractor, same grid, same field family, opposite verdicts | A-014i, and ✗15's condition — *"the grid resolves the surface"* — measured for this extractor rather than inherited. **The structural consequence is the finding.** Case (5) needs a residual with `g > 1` and `ℓ > 8`, which needs many crossings per edge, which needs the surface to oscillate several times inside one cell. That *is* under-resolution, by definition. So every fixture that exercises subdivision is necessarily in the regime where manifoldness fails for reasons that have nothing to do with subdivision — and it is why the seven reference fields, all comfortably resolved at 17³–33³, never reach it. A-014i's "T-001 validity on the fixture" has to be re-read as a comparison against the resolved case, not as an assertion of zero. **Determinism, by contrast, is a clean gate and passes**: two extractions of the subdividing field are byte-identical, which T-004 could never have checked because it only runs the reference fields. **The obvious control was the wrong instrument**, and worth recording as such: classic Marching Tetrahedra on the same field returns **2,113** vertices against subgrid's 73,855 and looks far cleaner, because M-67 measured a sign test as blind to 95.6% of tetrahedron configurations — it is reconstructing a much simpler surface, not the same one more carefully |
| M-201 | **The subdivision path has a grid-level fixture at last, and it says the orphaned vertices are real — 9.8% of them — which falsifies the inference M-200 was heading toward.** Every measurement of case (5) until now was on a bare tetrahedron built from hand-chosen edge coordinates, the one configuration where a crossing has no neighbour to be shared with. M-200 observed that dropping unreferenced positions on the *reference* fields changes vertex order but not vertex **count**, and the natural inference — that such positions are always consumed by an adjacent tetrahedron, so the bare-tet orphans are a fixture artefact — is **wrong**. Searched rather than assumed: `sin(k(x+y+z))` at `k = 96` fires case (5) **276 times** on a 9³ grid, and extracting it gives **166,591 vertices, 219,084 triangles and 16,296 referenced by nothing** | A-014i. **Why the neighbour argument holds on the reference fields and fails here.** A parent's crossings are real grid crossings carrying global identities, so an adjacent tetrahedron that *does not* subdivide reuses them — which is exactly why they only reorder on the seven reference fields, where `subdivision == 0` is pinned. A subdividing parent's children instead re-derive those same points as **child-local positions beyond `crossings.len()`, with no identity at all**, so where both sides of a face subdivide nobody consumes the parent's registration. **That names the fix and rules out the tempting one:** children should inherit the parent's crossing identities on parent edges rather than re-deriving anonymous copies; dropping the orphans afterwards is the shape that does not work, because it removes positions the reference fields legitimately share. **Finding a usable field was itself the hard part.** `sin(kx)·sin(ky)·sin(kz)` reaches case (5) most readily of anything tried — subdivision wants two transverse families of cuts — and its gradient **vanishes where those families meet**, so the extractor refuses it with `DegenerateNormal`. The obvious way to write the configuration is also the degenerate way; the plane family works because on `sin = 0` the cosine is `±1`, so `‖∇f‖ = k√3` everywhere on the zero set |
| M-200 | **Two of A-014i's three recorded defects are not what the review called them: one is real and its obvious fix is measurably wrong, and one is not a code fault at all but a hole in the test suite.** The ticket's own instruction was that the three names are *"one-line summaries, not verified re-derivations"*, so each was reproduced before anything was touched. **The fixture is derivable rather than searched for:** case (5) is taken only when `ℓ ∉ {4, 8}` **and** `loop_count > 1`, so with `ℓ = 4(d₁+d₂)/g` and `g = gcd` it needs `g > 1` and `d₁ + d₂ > 2g`. The smallest is **(4, 2)** — `g = 2`, `ℓ = 12`, two loops, `4(d₁+d₂) = 24` intersections, which is exactly the *"≥24 edge intersections"* the ticket quotes. It is also the **only** one of the twenty patterns in `every_implemented_case_emits_an_intersection_free_patch`'s sweep that reaches case (5), so that sweep has covered this path all along and nothing said so | A-014i. **Orphaned vertices: real, 24 of 112 positions for 52 triangles**, and the 24 is exactly the parent's own crossing count — `fill_append` records a position per named crossing because every case *except* subdivision triangulates over exactly those, and subdivision hands the work to four children that re-derive their own. **The obvious fix is wrong, measured:** dropping unreferenced positions at the end of `fill` and remapping takes the fixture to 88 positions and 0 orphans, and drifts six golden rows — `gyroid` and `fbm_terrain` at all three resolutions — with **identical vertex and triangle counts and a different hash**. So on the reference fields it removes nothing and only reorders: an unreferenced position there still carries a global crossing identity and seeds the map a neighbouring tetrahedron reuses, so dropping it moves which tetrahedron emits that vertex first. A-014i's own guard is that reference-field output must not change, and that guard is what caught it. Count pinned at 24, per M-4's rule. **The `1 − t` re-lerp: the code is correct and nothing tested it.** `fill_subdivision` already reverses the list *and* maps `t ↦ 1 − t`, and **inverting that branch passes all 32 tests in the module**, including the intersection-free sweep — because the `child.coordinates() == predicted` guard counts crossings per edge and a mirror leaves every count identical. The assertion with power is that a parent crossing must be re-derived **bit-identically** by whichever child expresses it: 24 of 24 with the code as written, **12 of 24 under the mutation**, verified in both directions. That is ✗-tier "verify a property test can actually fail" applied to a defect report rather than to a test |
| M-199 | **Subgrid Marching Tetrahedra is already intersection-free on every reference field at every resolution — so §3.2.3, whose entire purpose is to remove self-intersections, has nothing to remove and A-014d has no reachable target.** §3.2.3 turns an immersed Δ-complex into an **embedded** simplicial complex; the paper is explicit that the coincident pairs it addresses *"define manifold connectivity, but degenerate geometry"*. M-194 falsified the manifoldness half of A-014d's acceptance. This measures the half that was actually left: T-002's counter over the extractor's own welded output, which had never been run there — `subgrid/surface/tests.rs` runs it per tetrahedron on synthetic patterns, a different question. Result: **0 intersecting pairs on all 21 rows**, seven fields × 17³/25³/33³, up to 32,364 triangles on `gyroid` at 33³ | A-014d, closed as **not required** rather than blocked. **The control is what makes the twenty-one zeroes mean anything.** Dual contouring, through the *identical* pipeline — same fields, same resolution, same `epsilon_for(cell)` weld, same counter — reports **`gyroid` 18 pairs (7.143 per 1k)** and **`fbm_terrain` 21 pairs (24.249 per 1k)**, which is A-009's residue (M-28, M-29) showing up exactly where it was measured before. So the counter is live on this path and subgrid's zero is a property of the surface. **The first control tried was worse than useless and the test caught it**: counting the *unwelded* soup, on M-93's note that neighbours sharing a position but not an index are counted, reads zero too — face-touching is not crossing — and the reachability assertion failed rather than passing quietly. **What would reopen this ticket** is a field or resolution where this census goes non-zero; nothing else. Cost: the census adds ~120 s, dominated by subgrid extraction rather than by the intersection test |
| M-198 | **The publish job ran for the first time in its life, and failed on a secret that was never set — so the release pipeline has still never published anything, and now the reason is known rather than inferred.** M-174 established that `0.0.3` and `0.0.4` reached crates.io **by hand**, and that the job had been `skipped` on every run since it was written because the suite was red. Closing the last of that red (I-007, and this session's four tickets) made the whole suite green for the first time, which is the only condition under which `needs: [lint, test, bevy, msrv, package]` lets the job execute. It executed, and died on `CARGO_REGISTRY_TOKEN is not set` — the secret exists neither at repository level nor in the `crates-io` environment the job declares. **It never reached `scripts/publish.sh`**, which would have found `0.0.4` already on the registry for both crates, skipped both and exited 0 | GPU-013's push, and M-174's rule collecting its third instance. **The check was in the wrong place, and the script's own header said so.** `publish.sh` opens by explaining that it is version-driven precisely because *"the alternative … would leave main permanently red and train everyone to ignore it"* — and the workflow then required a token **unconditionally**, on every push to a green main, including the ones that upload nothing. The check now lives next to the upload it guards, and both directions are verified rather than argued: with no token and nothing to publish it exits **0** (`published 0, skipped 2`), and with the workspace version bumped so an upload is imminent it exits **1** with the loud error **before** `cargo publish` runs. **What is still owed is the user's:** the secret has to exist before any real release, and nothing in CI can tell you it is missing until a version bump makes an upload imminent — which is the correct trade, but worth knowing in advance rather than at release time |
| M-197 | **GPU-013's two halves both misnamed their target, and the measurement says the prize is waits rather than submissions.** **(a) The duplicated prologue is not between the functions the ticket names.** `extract` has no prologue at all — it delegates to `extract_buffers` and reads two buffers home. The ~30 duplicated lines are between `extract_indirect` and `extract_buffers`, and the ticket's stated justification is also wrong: the indirect-draw budget-clamp fix (`99415af`) landed **15 and 74 lines past the end** of the shared region. The real justification is that the copies had **already drifted** — `counts` carried `STORAGE \| COPY_SRC` in one and `STORAGE \| COPY_DST \| COPY_SRC` in the other, and both explanatory comments existed in only one. The unified version drops `COPY_DST`: nothing in the crate copies or writes into `counts`, so the flag was dead where it existed rather than missing where it did not. **(b) Batching the two geometry read-backs is worth 8–9%**, measured on an RTX 3090 at `gpu+read`: 97³ **0.835 → 0.764 ms**, 129³ **1.341 → 1.223 ms**, with after-run spread ±0.005. The saving is roughly *constant* rather than proportional to the bytes moved, which is the signature of removing a fixed synchronisation rather than bandwidth — `poll(Wait { submission_index: None })` drains the entire queue whether it carries four bytes or five megabytes | GPU-013. **The number that decides the rest of the ticket** comes from the harness this added, because `gpu_vs_cpu` had never timed `extract_indirect` at all and a claim about its submissions cannot be checked against a path nothing measures. The zero-read-back path costs **0.155 / 0.169 / 0.188 / 0.265 / 0.264 / 0.304 ms** at 17³…129³ — nearly flat across a **435×** rise in cell count, which refines M-160's "~0.17 ms" and confirms it is pure CPU recording. All **eight** of its `queue.submit` calls live inside that 0.304 ms. Against them, `extract_buffers` costs 0.552 ms at the same size, and **the entire 0.25 ms difference is one `poll(Wait)`** — the counts read-back that sizes the geometry buffers, which cannot be batched away because its value is what sizes them. So batching submissions can recover at most a fraction of 0.3 ms and removes no wait; batching read-backs removes a whole wait and measured 8–9%. **Submissions are not where this crate's time goes, waits are** — M-159 and M-167 restated with the two paths side by side. The dispatch-submit batching is therefore answered rather than deferred, and it carries a cost besides: `ExtractTimings`' per-stage `count_ms`/`scan_ms`/`emit_ms` are CPU wall-clock around each submission, so collapsing them into one encoder would destroy the attribution M-149 and M-160 both rest on, and the alternative — GPU timestamp queries — needs a device feature this crate does not request |
| M-196 | **Consuming `fields::` instead of a hand-rolled copy is a speed-up, not a tidy-up — and a dispatch layer can throw the whole benefit away silently.** `Sdf::gradient`'s default is central differences: **six extra `sample` calls per normal**. `fields::BoxExact`, `Sphere`, `Difference` and `Intersection` all override it analytically; the examples' private re-implementations implemented `sample` alone, so every normal they produced cost seven field evaluations instead of one. Measured on `game_csg_props` at 41³, extract time per frame: dual contouring **5.81 → 5.35 ms**, dual contouring with the clamp off **5.82 → 5.25 ms**, surface nets **2.82 → 2.76 ms** | E-212. **The trap is the dispatch enum.** `game_destruction`'s `Solid` wraps three fields and implemented only `sample`, which means swapping its variants for the crate's types would have bought **nothing** — the default `gradient` would have been reimposed one layer above the fields that carry the real one. Forwarding `gradient` through the enum is what makes the swap real, and nothing about the types would have warned that it was missing. **Two things this measurement cannot say.** The demo's accuracy columns are driven by a *moving* reflex edge on a wall clock, so a fixed frame index samples a different phase run to run — the baseline read `dc worst 0.3700` where the changed build read `0.2255`, and neither is evidence about the other. What *is* checkable is that the changed build reproduces the figure the demo's own HUD prose has recorded all along (`0.2255`, and surface nets `2.34×` further out), which the baseline run did not reach by frame 61. A like-for-like accuracy comparison needs a deterministic sweep this demo does not have |
| M-195 | **The mirrored seam's crack is exactly invisible to a seam-plane counter — 0 edges on the plane, 28 on each side of it.** E-211 said `game_lod_flyover`'s HUD scans the seam plane while Eq-4.2's taper puts transition geometry at `seam ± w`, and estimated the blind spot from the 44 boundary edges the 2026-08-14 review found. `the_mirrored_seam_closes_only_with_the_signed_width` now measures it directly, headless and in CI, on the configuration E-107 never ran — coarse block **below** the seam, its high face tapered to `seam − w`. With the signed width the three planes `(seam, seam − w, seam + w)` report **`(0, 0, 0)`**. With the unsigned width — the one-token revert of the fix — they report **`(0, 28, 28)`** | E-211. **The middle zero is the finding.** 56 boundary edges are open and *not one of them lies on the seam plane*, so the blind spot is not a matter of tolerance but of geometry: the taper moves every open edge off the plane the counter watches, and `w ≥ BASE_H` against a `BASE_H · 0.25` window puts the band **four tolerances** outside it. An instrument reading 0 on a mesh with 56 holes in it is the strongest possible statement of E-208's rule, and this is the first time this project has had the number for a case where the instrument was structurally blind rather than merely starved. The test carries its own positive control — the unsigned run is asserted non-zero — because M-44 and E-208 both say a zero that cannot go non-zero is not a measurement, and here the zero being reported *was* the bug. **Confirmed end-to-end in the example too**, run on this machine's RTX 3090: with the band-scanning counter the fixed build reports `open_low 0, open_high 0`, and reverting the one token `sample_width` → `width` makes it report **`open_low 44`** — the review's own figure, recovered by the instrument that could not see it before, and only on the rows where `seams_low` is 1, which is the mirrored configuration and no other |
| M-194 | **A-014d's acceptance criterion asks §3.2.3 to fix a defect §3.2.3 does not cause, and this repository's own instrument has been saying so, in a passing test, since it was written.** With the source finally in hand (V-28) the inset was implemented three ways and measured on `csg_difference` at 17³ against a baseline that is **closed** — `(3, 6, 6)`, 0 boundary edges, 1 component, χ = 5. **(a)** A midpoint on every polygon edge contained in an edge of `f`, per the prose: `(0, 12, 6)` — but **48 boundary edges, 9 boundary loops, 4 components**. **(b)** One Steiner point on the segment across the omitted distinguished corner, per the reference implementation's corner-type branch: `(0, 12, 0)`, and the same 48/9/4. Both tear, and for one reason: in this decomposition a region is *one face's share* of a disk, so **every part of its boundary is already shared** — a chord with the coincident copy in the neighbouring tet, an edge piece with the same tet's region on the adjacent face, the corner cut with §3.2.2's cap. The paper's comb face spans the whole tet boundary, where those segments really are free. **(c)** Displacement moved to the only place this decomposition has room for it, the interior — one Steiner at the ring's centroid pushed toward the corner the shared face does not carry, coning over the ring. Closed again (χ = 5, 1 component, 0 boundary edges) and **`nm-edges` unchanged at exactly 3** | A-014d. **(c) is the result.** It separates the two copies perfectly — different apexes, no shared interior edge — and the non-manifold count does not move, because **duplication was never the cause**. `the_defects_traced_back_to_the_tetrahedra_that_made_them` already reports `csg_difference` as **`duplication-only 0, three-distinct 3`**, and its own comment states the criterion: *"if fewer than three remain, duplication is the whole story and §3.2.3's inset is the fix. If three or more distinct polygons genuinely meet, **it is not**."* All 3 are genuine three-sheet junctions — 4 faces, 3 distinct triangles, 2 tets — so one tet contributes two sheets and the other two, and the pair being identical is incidental. The paper agrees in advance: *"Pairs of identical polygons … **define manifold connectivity**, but degenerate geometry"* — §3.2.3 buys intersection-freedom, not manifoldness, and T-002 is the counter it should be graded against. **Where the inset does have a target:** `fbm_terrain` is `duplication-only 2` of 4 bad edges. That is the whole reachable prize, and it is on a different field from the one the ticket names. **This is ✗11's rule happening a second time** — *"a ticket's acceptance criterion is itself a claim about the code; check it against the code before starting"* — and the check cost one `cargo test` against a counter that was already green |
| M-189 | **Running the subgrid validity census at the three resolutions its own gate names doubles the number of known defects, and every new one is on a field the suite called clean.** T-010 widened `the_validity_suite_over_every_reference_field` from `n = 17` to `[17, 25, 33]`. `(non-manifold edges, non-manifold vertices, inconsistently oriented)`, the rows that are not all-zero: **17³** `csg_difference` (3, 6, 6), `gyroid` (0, 0, 138), `fbm_terrain` (4, 6, 19); **25³** **`torus` (4, 6, 6)**, `thin_plate` (0, 0, 8), `gyroid` (0, 0, 150), `fbm_terrain` (2, 3, 29); **33³** `csg_difference` (0, 0, 36), `thin_plate` (0, 0, 6), `gyroid` (0, 0, 330), `fbm_terrain` (8, 12, 53) | T-010. **Three of the seven fields were believed clean and are not**: `torus` is non-manifold at 25³ only, `thin_plate` is incoherently wound at 25³ and 33³ only, and `csg_difference` — whose 17³ row is A-014d's entire target — is **clean at 25³ and non-manifold at neither finer size**. The defect that ticket exists to remove occurs at one of three resolutions and was never checked at the other two. **The lesson is not that 17³ is a bad choice, it is that one sample of a stated range is not a sample of it**: every one of these fields is defect-free at some resolution and defective at another, so any single size gives an answer that is true and unrepresentative |
| M-187 | **Orienting each connected component from its most confident triangle drives the inconsistently-oriented-edge count to *exactly zero* on every mesh whose edges are all manifold — and the residue where it does not is non-manifoldness, not a limit of propagation.** Before → after, all seven fields at 17³/25³/33³: **`gyroid` 138 → 0, 150 → 0, 330 → 0**; **`thin_plate` 0 → 0, 8 → 0, 6 → 0**; `csg_difference` 6 → 6, 0 → 0, **36 → 0**; `fbm_terrain` 19 → 6, 29 → 3, 53 → 12; `torus` 0 → 0, 6 → 6, 0 → 0; `sphere` and `box_exact` 0 throughout. **The law is exact and holds on all 21 rows: every row with `non_manifold_edges == 0` comes out at 0, and every row with a residue has non-manifold edges** | A-014f, `every_reference_field_comes_out_coherently_oriented`. The mechanism is forced rather than incidental: an edge carrying four faces cannot have them pairwise oppose, so no assignment of windings makes it consistent — that residue is A-014d's defect wearing an orientation costume. Vertex counts, face counts, manifoldness and boundary counts are asserted unchanged across the pass, so it moves winding and nothing else. **The seed is the most *confident* triangle rather than the first**, which is load-bearing and has its own test: a triangle whose gradient lies in its own plane has a vote of ~0, and seeding from one flips a whole component on a coin toss |
| M-188 | **Two of the seven reference fields are non-manifold at a resolution the suite never tested, and one of them is `torus`.** Measured while establishing M-187's law, which needed the non-manifold count per row: `torus` reports **4 non-manifold edges at 25³** and 0 at 17³ and 33³; `csg_difference` reports 3 at 17³ and **0 at both 25³ and 33³**; `fbm_terrain` 4, 2, 8. So `csg_difference`'s pinned `(3, 6, 6)` — the row A-014d exists to drive to zero — is a **17³-only** defect that the field does not exhibit at either finer resolution | A-014f. This is the same hole as M-182 and it is worse than it looked: the subgrid validity census runs at `n = 17` alone (T-010), and 17³ is not merely *a* sample of the three, it is the resolution that decides which defects this project believes it has. A-014d's target exists at one of the three resolutions it was never checked at the others of |
| M-184 | **Naming a root by the grid point it lies on completes the identity outright: on every field the extractor now emits each vertex exactly once, and a positional weld removes nothing.** Raw → welded, all seven at 17³: `sphere` **812 → 812**, `torus` 912 → 912, `box_exact` **338 → 338**, `csg_difference` **482 → 482**, `thin_plate` **422 → 422**, `gyroid` **4014 → 4014**, `fbm_terrain` 1758 → 1758. Every one of those raw counts is exactly what a positional weld used to reach, so M-169's gap is closed and M-180's ceiling is beside the point — it bounded rules that leave positions alone, and this moves them | A-014h. The rule is not a cleverer key: an endpoint root is **emitted at the grid point's own position**, `origin + cell_size · index`, the expression `corners` already uses and every cell computes identically (M-32). Detection is `f(G) == 0` exactly — no tolerance — because that is the definition of the surface passing through `G` and it is a property of `G` alone, so every one of the up-to-24 tetrahedra meeting there reaches the same answer without consulting the others. **The cost is 15 of 168 golden rows, and its shape is the check that it is right**: 14,171 vertices removed and **28,342 triangles — exactly 2 per vertex, on all 15 rows independently**. That is the edge-collapse signature (`ΔV = −1, ΔE = −3, ΔF = −2`), which preserves the Euler characteristic, and the validity suite's χ = 2 rows confirm it did. **The 6 rows that held are `torus` and `fbm_terrain` at all three resolutions — precisely the two fields M-169 measured with zero vertices on a grid point.** The correlation that opened the ticket also closes it |
| M-185 | **Completing the identity turns a sliver into a repeated-index triangle, and the extractor now declines to emit those.** When the surface passes exactly through a grid point, the two tetrahedron edges meeting there each carry a root *at* it; naming both by the point makes them one vertex, and any triangle spanning both becomes a triangle with two equal indices. `validate_indexed` counts that as a **structural error** — 36 of them on `sphere` at 17³, `INVALID: 36 violations` — not as the recorded metric it treats a near-zero-area sliver as | A-014h. Declining to emit them is not dropping geometry: they carried no area before the merge either, and the alternative is shipping an index buffer this crate's own validator rejects. **What confirms it rather than merely excusing it**: with them gone, `sphere` validates 0/0/0/0 with χ = 2 and one component, and every pinned defect census in `the_validity_suite_over_every_reference_field` is **unchanged** — `csg_difference` (3, 6, 6), `gyroid` (0, 0, 138), `fbm_terrain` (4, 6, 19). A geometry change that moved 15 golden rows disturbed not one of the numbers those rows exist to protect |
| M-186 | **M-162 is falsified, and A-014d's blocker with it: after A-014h no coincident polygon has a foreign tetrahedron on its boundary edges, and none is cross-cell.** M-162 is what blocked A-014d — *"the inset needs information no per-tetrahedron and no per-cell rule carries"* — measured as `csg_difference` carrying **33 coincident polygons, 27 of them with foreign tetrahedra standing on their boundary edges, 312 such triangles and 150 of them in a different cell**. Re-measured on the same instrument with the identity complete: **`csg_difference` 3 coincident, 0 with foreign edge users, 0 cross-cell, all 3 within one cell**; `fbm_terrain` 4, likewise all within one cell; **every other field 0, `box_exact` included** — the field M-161 called *"the most coincidence of any field, 30 polygons and 348 triangles"* | A-014h, `which_polygons_coincide_across_a_shared_face`. **So most of what M-162 measured was never §3.2.3's immersion at all** — it was one grid point wearing many names, and it disappeared when the names were fixed rather than when any geometry was inset. `the_defects_traced_back_to_the_tetrahedra_that_made_them` agrees from the other side: `thin_plate`'s 2 soup non-manifold edges → **0**, `csg_difference`'s 6 → **3**, and the 3 that went are exactly the *duplication-only* ones, leaving the 3 M-163 identified as genuinely three-distinct polygons. **A-014d is now a per-cell problem on 3 edges of one field**, which is a different ticket from the one that was blocked |
| M-182 | **P-7 confirmed, and the check that confirmed it found `thin_plate` failing at two of the three resolutions nobody was testing.** P-7's limb (b) holds exactly as registered: welded `thin_plate` is **one component, 0 boundary edges, χ = 2, 0 non-manifold edges and 0 non-manifold vertices at 17³, 25³ and 33³** — a closed genus-0 surface, so its top and bottom are joined at a rim and are one component *by topology rather than by tolerance*. Both remaining falsifiers fail to fire: the component count is 1, and it stays 1 at every epsilon from `h·10⁻⁹` to `h·10⁻²`. **But the inconsistently-oriented-edge census is `17³ → 0`, `25³ → 8`, `33³ → 6`** | A-014f/P-7, `the_plate_is_one_closed_orientable_component_at_every_tolerance`. **The reason nobody had seen it is the interesting part: `the_validity_suite_over_every_reference_field` runs at `n = 17` and nothing else**, while Phase 1's own gate says *"no unexplained violations on all seven test fields **at three resolutions**"*. The subgrid extractor has been meeting one third of its stated gate, and 17³ is the resolution where the plate happens to be clean. **This strengthens the case for A-014f's remedy rather than weakening it**: the surface *is* orientable and there is a coherent orientation to propagate to — A-014e's per-triangle vote simply is not finding it, on the very field the per-triangle rule was written to protect. P-7's worry was that propagation would invert a sheet; the measurement says the plate is a single closed orientable component, which is the precondition under which propagation is safe |
| M-181 | **T-009's own premise was wrong: normalising four scattered weld epsilons onto one moved nothing, because the weld's answer is flat across *six* orders of magnitude.** The ticket said *"normalizing shifts pinned non-manifold censuses — which is exactly why this is a ticket rather than a drive-by."* Four policies were in use (`h·1e-4` via the constant, `h·1e-4` as a literal, `h·1e-6`, `h·1e-5`) across 38 call sites. Routing all of them through one changed **not a single pinned count** anywhere in the suite. Welding the subgrid output at `h·10⁻⁹ … h·10⁻³`: every field returns the *same* vertex count at all seven — `sphere` 812, `torus` 912, `box_exact` 338, `csg_difference` 482, `thin_plate` 422, `gyroid` 4014, `fbm_terrain` 1758 | T-009, `the_weld_answer_is_flat_across_the_range_the_four_policies_spanned`. **The plateau's edge is the useful half.** First factor whose answer differs from the policy: `h·10⁻²` on `gyroid` (4014 → 3932) and `fbm_terrain` (1758 → 1735) — the two multi-sheet/rough fields — and `h·10⁻¹` on `sphere`, `torus`, `csg_difference` and `thin_plate`. `box_exact` never changes, even at `h/2`, because its geometry lives on lattice planes and its genuine neighbours are a whole cell apart. **So the policy sits 100× below the nearest tolerance that changes any answer, and 1000× below on five of seven fields.** That margin is why four different numbers coexisted for months without disagreeing — and exactly why they were dangerous: they agree until a field puts two real vertices `10⁻⁵·h` apart, and then they stop agreeing silently |
| M-179 | **A-014h's stated mechanism finds almost nothing, and the reason is a deliberate choice in the root finder.** The ticket — and M-169 before it — says *"a crossing at parameter 0 or 1 should be named by the grid point it lies on."* Measured over every tetrahedron edge of every cell at 17³, across all seven fields: **no root anywhere reports `t == 0`**, and **36 report `t == 1`, all on `gyroid`**, out of 52,488 roots. `refine` returns the **upper** end of its final bracket by design, so it can keep the ascending-and-distinct contract when a root sits on a sample — which means a root at an edge's *lower* endpoint comes back as a tiny positive parameter and never as zero. The position behaves asymmetrically too, and that was not predicted either: `a + (b − a)·t` **never** rounds back onto `corners[a]` on any field, while it lands exactly on `corners[b]` 1,692 times on `box_exact`, 1,284 on `csg_difference`, 112 on `thin_plate`, 36 on `gyroid` and 18 on `sphere` | A-014h, `no_root_reports_parameter_zero_and_almost_none_reports_one`. The two endpoints were assumed to behave alike; only one does. A `t == 0.0 \|\| t == 1.0` test would have been written, passed its own unit test, and moved nothing |
| M-180 | **The population M-169 named as one defect is two, and only the smaller one is an identity problem — so A-014h cannot reach its stated acceptance.** M-169 measured that a positional weld removes vertices identity-based sharing could not, and concluded the fix was a second keying rule. That assumes the copies *are the same point*. Counting raw vertices, raw vertices **distinct by bit pattern**, and vertices after the weld: `box_exact` **1262 → 1028 → 338**, `csg_difference` **1280 → 1094 → 482**, `sphere` **830 → 830 → 812**, `thin_plate` **450 → 444 → 422**, `gyroid` **4020 → 4014 → 4014**, `torus` and `fbm_terrain` unchanged throughout. **The middle column is the ceiling on any rule that leaves positions where the extractor put them**, and on `box_exact` it is 1028 against the weld's 338: an exact rule can remove 234 of the 924 duplicates, and 690 are positions IEEE calls different numbers | A-014h, `how_much_of_the_positional_weld_an_exact_identity_could_ever_reach`. Two tetrahedra reach one grid point along different edges and `a + (b − a)·t` rounds differently on each, so merging them is not *naming* a shared vertex, it is **moving** two distinct ones together. **The split runs by field and is not uniform**: `gyroid`'s 6 merges are entirely exact, so an identity closes it completely; `sphere`'s 18 are entirely inexact, so an identity closes none of it. Same symptom, opposite causes. What A-014h now has to decide is whether the extractor may canonicalise an endpoint root onto the grid point's own position — a geometry change of a few ulps, not a keying change |
| M-175 | **Sorting by magnitude is not enough to make a sum permutation-invariant, and the comment that said it was had been wrong since ✗12.** `sort_by_magnitude`'s doc argued that *"ties need no tie-break, since IEEE addition and multiplication are both commutative — only associativity fails — so two terms of equal magnitude give the same answer either way round."* True of `a + b` against `b + a`; false of every longer sum, which is the only kind that occurs here. Witness, measured: smallest-magnitude-first, `[1e-16, +1, −1]` sums to **0** and `[1e-16, −1, +1]` to **1.11e-16**. Both are correctly sorted by magnitude — the `±1` pair *ties*, so a stable insertion sort left their order to the order they arrived in, which in this crate is the axis or edge labelling that a lattice rotation permutes | A-016, `sorting_on_magnitude_alone_would_still_depend_on_the_arrival_order`. Found by writing the padding test, not by reading the code — the first three attempts at the test all failed for this reason and were read as fixture bugs before the third one was traced. **The fix is to compare the signed values on a tie** (`−c` before `+c`), which makes the sequence a function of the multiset alone. **Why the original claim survived ✗12's 9600 trials:** `dot_equivariant`'s terms are *products of two components*, and a signed permutation flips both signs inside each product, so its terms permute **without negating** and a tie group is always repeats of one value. The claim was true in the only place it had been tested |
| M-176 | **Zero-padding a reduction is transparent, negative zero included — and the reason is the accumulator's seed, not the sort.** A-016 reduces over twelve slots keyed by edge label with absent edges at `R::ZERO`, so the padding question is load-bearing: a twelve-slot reduction over three real terms must give the same bits as a three-slot one or every cell with fewer than twelve crossings would answer differently for nothing. It does. The one place an exception was predicted — every real term `−0.0`, where the padded form was expected to return `+0.0` against an unpadded `−0.0` — **does not occur**, because `sum_equivariant` seeds at `R::ZERO` and `+0.0 + (−0.0)` is `+0.0` under round-to-nearest, so *both* return `+0.0` | A-016, `negative_zero_survives_padding_because_the_accumulator_starts_positive`. The predicted exception came from reasoning about the sorted sequence and forgetting the seed. Pinned in both directions rather than deleted, because the golden hash deliberately distinguishes signed zero and this is exactly the kind of claim that gets re-derived wrongly later |
| M-177 | **Reordering cannot buy negation equivariance, and the obstruction is structural rather than a missing tie-break.** A-016's brief was to make the QEF accumulation "a function of the set of terms". Permutation invariance is achievable and is delivered. `φ(−S) = −φ(S)` is **not achievable by any ordering rule**: a magnitude tie group holds `m` copies of `−c` and `n` of `+c`, negation swaps those counts, and no order that is a function of the multiset can map the group onto its own reverse unless `m == n`. Witness: `[1e-16, +1, −1]` sums to `1.11e-16`, its negation to `0` | A-016, `negation_equivariance_is_not_achievable_by_ordering_and_here_is_the_witness`. **This bounds what ✗12's "bit-exact equivariance" can mean for a sum of *components*.** A lattice rotation can negate a component, so a sum of positions or of normal-outer-product entries is not bit-exactly equivariant under the full octahedral group by this route — where the three-term dot product is, for M-175's reason. The property the extractors actually need and now have is the weaker, sufficient one: the vertex is a function of the crossings, not of the edge labels they arrived under. Anything stronger needs exact summation, not a better sort |
| M-178 | **A-016 moved 34 of 168 golden rows, not the 42 predicted, and the 8 that held identify their own mechanism.** All 21 dual-contouring rows and all 21 manifold rows were expected to shift, since both share the solve path. 17 of 21 did, twice over: **34 hashes changed, 0 vertex or triangle counts changed**, and no algorithm outside those two was touched. The four that held are `box_exact` and `thin_plate` at **17³ and 33³** — the two axis-aligned fields at the two resolutions whose cell size is a binary fraction of the domain (`16` and `32` cells; `25³`'s `24` is not). On axis-aligned geometry at an exactly-representable spacing the crossing coordinates are exact, their sums are exact, and reordering cannot change a result that was never rounded — which is why the same two fields *did* move at 25³ | A-016. **This is M-94's fixture trap in its useful direction**: the fields that are easiest to reason about are also the ones least able to detect a floating-point defect, and a suite that only used them would have reported this fix as a no-op |
| M-172 | **A signed distance field's gradient is exactly zero on its medial axis, and for a slab the medial plane is precisely where you would aim.** `game_destruction` fires charges at `Vec3::new(u * 1.9, v * 1.5, 0.0)` — `z = 0` is the wall slab's mid-plane, equidistant from both faces — and `BrushStack::gradient` there returns `[0.0, 0.0, 0.0]`, measured on two consecutive shots. On a third, off-centre shot it returned `[0.17, 0.98, 0.00]`: non-zero, and pointing **along** the wall rather than out of its face, because the nearest surface to that charge was the top edge | E-204's z-fighting fix. "Which way is out of the solid" reads as the obvious use for a gradient and is wrong twice over here — it is undefined exactly where the geometry is symmetric, and where it is defined it names the nearest surface, which is not necessarily the one you meant. **The general rule is that a gradient answers "where is the nearest surface", never "which way did the projectile go"**; the second question needs a second piece of information, and in this case the camera already had it |
| M-173 | **Two coincident surfaces are not a rendering bug to be biased away, they are a modelling statement — and the fix is to stop the geometry coinciding.** `game_destruction`'s fragment is `solid ∩ charge` and its crater is `solid − charge`, so the two share a surface *exactly*: the same sphere, extracted twice from two different fields. A fragment left at rest in the hole it came from renders as a chunk fused into the wall | E-204. **Velocity alone did not fix it, measured**: some fragments leave at 9 m/s while others sit at the impact with a velocity of **0.004 and no response to gravity**, on bodies whose mass is a healthy **1.19** and whose sleeping is disabled — so the cause is inside avian's solver and not somewhere the example can reach. What does fix it is spawning the fragment at the **mouth** of the crater rather than inside it: a fragment that fails to move then fails *clear of the wall*, which is the only place its faces are not coincident by construction. **Depth bias would have been the wrong instrument** — it would hide the overlap at one camera angle and not another, where moving the geometry removes the overlap at every angle |
| M-171 | **A shipped example panicked on every run, and three siblings carried the same race unfired.** `game_terrain_stream` dies with *"Entity despawned: the entity with ID 404v0 is invalid; its index now has generation 1"* — `attach_meshes` queues `commands.entity(e).insert(Mesh3d…)` for a chunk in the same frame the streamer queues `commands.entity(e).despawn()`, and whichever order the queue flushes, the insert can land on a dead entity | Found while re-recording the README's GIFs, not by a test — nothing in CI opens a window, so no example is ever *run*. `game_walk`, `game_showcase` and `game_capsule_walk` have the identical pattern and had simply never streamed aggressively enough to hit it. Fixed with `try_insert` at all five sites. **That is not a swallowed error and the distinction matters here**: attaching a mesh to a chunk that no longer exists is *vacuous*, not degraded — there is no other thing it could correctly do — where a fallback would be picking a different, worse outcome for the same input. **The gap this exposes is coverage, not the fix**: 30 examples compile in CI and none of them runs, so a panic on frame 200 is invisible until a human launches it |
| M-174 | **CI has been red on every push of the GPU series, and behind the one accepted failure three unrelated ones accumulated unseen — including the two that block the release itself.** The `test (ubuntu)` job runs `cargo test --workspace --all-targets` and GitHub's ubuntu runners expose no GPU adapter, so every isomesh-gpu device test has panicked there since GPU-001 (`headless.rs` refuses software rasterisers by design — GPU-009, M-147 — so there is nothing to fall back to; macOS runners expose Metal and stayed green). With that job red as the norm, nobody saw: (1) the lint job's rule-2 gate `grep -ril bevy crates/` matching **11 files** — every one a comment, README or manifest note *explaining* the wgpu-follows-Bevy pin, zero dependencies (verified: non-comment `.rs` matches = 0, comment-stripped manifest matches = 0, resolved-graph bevy count = 0); (2) the msrv job failing because `rust-version = "1.85"` was never true once isomesh-gpu existed — **wgpu/naga 29.0.3 and 29.0.4 both declare 1.87** (registry index), and the dev graph's `wide`/`safe_arch` (via parry3d) declare **1.89**, while `bevy_isomesh` claimed 1.85 against Bevy 0.19's own **1.95**; (3) a one-line rustdoc `redundant explicit link` in `bevy_isomesh/src/mesh.rs`. The `publish` job `needs: [lint, test, bevy, msrv, package]` and reports `publish to crates.io: skipped` on every run — four of its five gates red. **This finding's first draft concluded that 0.0.3 and 0.0.4 were therefore never uploaded, and checking it against the registry falsified that outright**: both are on crates.io (0.0.3 at `11:37:52Z`, 0.0.4 at `12:17:57Z`, 2026-08-14), uploaded **by hand from the developer machine seconds after each push, and before the publish job existed at all** — that job landed at `d01756f`, `12:22Z`, after both. So the release process was never dead; it was **manual**, and the automation written to replace it has still never run once. That is the worse of the two states and much the harder to see, because *versions appearing on crates.io is exactly the evidence you would use to conclude the pipeline worked* | Review of 2026-08-14, from `gh run list` (failures back through the whole visible GPU series) + `--log-failed` on run 31796991183, each cause reproduced locally; the upload claim re-checked against the crates.io API and `gh run view --json jobs` on 2026-08-14, which is what falsified it. Fixes in the same session: the ubuntu leg now tests the CPU crate (the one the cross-OS golden-hash claim is about) and macOS keeps the full workspace; rule 2's gate checks dependency positions and non-comment code instead of prose; `rust-version = "1.89"` verified by `cargo +1.89 check --workspace --all-targets` clean. **The method failure is the interesting part: an accepted red is indistinguishable from a new red.** GPU-009 accepted a *coverage* hole and left it expressed as a permanent *failure*, which converted CI from a signal into noise for the entire series |
| M-170 | **A GPU adapter is a finite resource and the test harness was treating it as free: 67 tests opened 67 devices, and once another process held 5.4 GB of the card, `request_adapter` stopped answering partway through the run.** Deterministic, not flaky — the same 11 tests failed on two consecutive runs, **every one of them passes in isolation**, and every failure is `Error::NoAdapter` | Found while verifying 0.0.3 for release. The boundary is visible in the log: device tests pass up to #38 and every device test after it fails, while pure-CPU tests keep passing throughout. **The crate's own no-software-fallback rule is what made this legible** — a fallback adapter would have quietly run eleven tests on a CPU rasteriser and reported success. **The fix is one shared device behind a `OnceLock`, and it is the more accurate harness rather than a workaround**: this crate's API takes `&wgpu::Device` and never opens one, so nothing under test wants a *fresh* device; only `Gpu::new` itself does, and `headless::tests` still calls it directly because that is its subject. **67 passed in 0.85 s against 6.72 s before — 8× faster and green** |
| M-168 | **Giving a crossing an identity instead of a position removes 5.01× of the subgrid extractor's vertices and changes no triangle at all.** Across the 21 golden rows: **389,146 → 77,605 vertices**, per-row **4.76× (`fbm_terrain`) to 5.08× (`torus`)**, and **0 of 168 golden rows changed a triangle count** | A-014g. The key is `(the tet edge's two grid points, the root's ordinal)`. It works because `TETS[t]` orders corners by inclusion, offsets are the grid point minus the cell origin, and subtracting the same origin preserves the componentwise comparison — so the edge runs the same direction from whichever cell it is viewed, and `FacePoint` already counts `index` from the lower corner for exactly this reason. **The point is not the vertex count.** A vertex now has an identity rather than only a location, so a later stage can move one and every triangle standing on it moves too. Welding by position cannot do that — after the move the copies no longer coincide — and that is precisely the tearing M-101 measured twice and M-162 traced to 150 triangles in neighbouring cells |
| M-169 | **Identity-based sharing is complete exactly when no root lands on a grid sample point, and the correlation is not approximate.** Raw vertices → welded, with the count sitting on a lattice point: `torus` **912 → 912, 0 on lattice**; `fbm_terrain` **1758 → 1758, 0**; `sphere` 830 → 812, 24; `gyroid` 4020 → 4014, 7; `thin_plate` 450 → 422, 58; `csg_difference` 1280 → 482, 1087; **`box_exact` 1262 → 338, and all 1262 are on a lattice point** | A-014g, `how_complete_the_shared_table_is_against_a_positional_weld`. **The two fields where the weld finds nothing to do are exactly the two with nothing on a grid point**, and the field where *every* vertex is on one loses 73%. The mechanism is that a root at a tetrahedron edge's **endpoint** sits on a grid point that up to 24 tet edges meet at, so it has one position and a different `(edge, index)` on each — one point wearing many names, which no correct sharing under this key can merge. That is a property of how the field sits on the grid, not of the algorithm: `box_exact`'s faces are axis-aligned and land on sample planes, which is **M-94's fixture trap showing up structurally rather than by accident**. **The remedy is still an identity and not a weld** — a crossing at parameter 0 or 1 should be named by the grid point it lies on rather than by the edge it was found along — and it is named rather than done, because it is a second keying rule and belongs in its own commit |
| M-167 | **Across the whole GPU series the arithmetic never moved and was never the point: synchronisation was 83% of an extraction and the payload 7%.** Breaking the 0.454 ms `extract_buffers` at 129³: **count pass 0.017 ms (3.7%), the four-byte read-back's own transfer 0.033 (7.3%), and the `poll(Wait)` around it 0.375 (82.6%)** | GPU-010a, GPU-012, GPU-011a, GPU-010b together. The cumulative figure is **15.01 ms → 0.54 ms, 28×**, and `count + emit` — the actual Marching Cubes — was **0.11 ms at the start and 0.04 ms now**. Every gain came from data movement or a stall removed: a CPU prefix sum and an 8.4 MB counts read-back (1.56×), a second pass over 8.4 MB during evaluation (1.18×), the entire upload (15.5×), and finally the wait itself. **Attribution matters here and the obvious citation is the wrong one.** V-10's atomics → LDS staging result is directionally the same thing. **V-9 is not support for it** — that is the compute-shader → mesh-shader 23.4× figure, and **M-149 measured in this repo that the mesh-shader draw removes the *smallest* of the three movement costs**, so quoting V-9 as confirmed here would contradict a correction this file already carries. The thesis stands on this series' own numbers |
| V-24 | **Custodio §5.1's correction, verified from the paper and re-derived rather than transcribed.** The swept saddle value is `f(x_c(t)) = (A_t·C_t − B_t·D_t) / (A_t + C_t − B_t − D_t) = F(t)/Δ(t)`, with **F quadratic in t and Δ linear**. Chernyaev's test tracks `F` alone: *"the polynomial F(t) … is a second order equation in t and thus can only allow for two sign changes. Therefore, the sign tracked by the MC33 algorithm will not match the expected one at some point"* | A-002c, `10.1016/j.cag.2013.04.004` §5.1, read this session. The quotient of a quadratic by a linear is a **hyperbola**, so where Δ has a root the saddle has a pole and the sign jumps across it — three sign changes on (0,1), which no quadratic can track. Their Figure 6 is a case 13.5.2 read as 13.5.1 for exactly this reason. **The three coefficients were re-derived here before being compared** — expanding `A_t = A₀ + tα` gives `a = αγ − βδ`, `b = C₀α + A₀γ − D₀β − B₀δ`, `c = A₀C₀ − B₀D₀` — and they agree with the paper's term for term, which is how rule 5 is satisfied without transcribing a formula |
| V-25 | **"Lewiner's reference implementation omits disambiguation for cases 10 and 12 entirely" is true of the *code* and false of the *literature*, and A-002b's row stated it as the latter.** The paper's §5.4 calls it *"a missing step in the implementation"* and then gives the rule: **10.1.1** positive nodes separated on both faces and at the cube diagonals; **10.1.2** separated on both faces but not at the diagonals; **10.2** separated on top and connected on the bottom face | A-002c, same paper §5.4 and §6.4. §6.4 is explicit that the fix is a code change — *"We fixed the MC33 implementation by adding the lines 16-20, which in the original implementation were replaced by the result case 10.1.1"* — so what was missing was never the published disambiguation. **This weakens A-002b's rule-5 stop from "no correct source exists" to "the source is prose and a figure rather than a table"**, which is a different and much smaller obstacle |
| V-26 | **A follow-up exists that builds MC33's triangulation with no lookup table at all, and A-002b's blocker did not know about it.** Custodio, Pesco & Silva, *An extended triangulation to the Marching Cubes 33 algorithm*, `10.1186/s13173-019-0086-6` (2019): vertices are labelled `+`/`−`/`=`, connected into **groups of vertices and groups of edges**, and each group's geometry is the boundary of a **convex hull** with the cube-face triangles removed | A-002c, read this session. Three categories — *simple leaves* (discs), *tunnel* (cylinder), *interior point* (disc plus a cell-centre vertex) — with case 13.3, 13.5.1 and 13.5.2 given as named combinations of them. The authors state the consequence directly: *"the triangulation is generated without the need of a look up table."* **So A-002b's "there is no correct published table to transcribe" is answered by a construction that needs no table**, and the `=` label additionally removes the degenerate slivers Marching Cubes emits by construction. The same paper also confirms the non-manifold remedy is a *grid* change rather than a table change: split the two cells at the shared ambiguous face's critical point |
| M-165 | **Among the configurations that can exhibit it, Chernyaev's numerator-only test is wrong 12.6% of the time — 1,966 of 15,625.** Every one of the 15,625 has a pole inside the sweep, by construction | A-002c, `how_often_the_correction_changes_the_answer`. The family swept is face pairs **opposed in sign** — `A`/`C` positive on the low face and negative on the high one — which is the structure Appendix A's counterexample has and the only one that can put Δ's root inside `(0, 1)`, since Δ is linear and a linear function positive at both ends is positive throughout. **This is not Custodio's rate and must not be quoted as one.** They report once in 10,000 random 5×5×5 fields and six times across 50 isosurfaces of the Skull dataset; that is a rate over *fields*, and it is small because the opposed family is itself rare. Within the family the correction is not a rounding detail |
| V-29 | **Two sources this repo recorded as unobtainable were obtainable, and one of them removes the largest single piece of A-002b.** `docs/research/2026-08-10-meshing-library-target.md:76-77` lists Grosso 2016 (`10.1111/cgf.12975`) and Grosso 2017 (`10.1145/3095140.3095179`) as **PAYWALL**. Both are in home-still, converted and indexed — `10.1111_cgf.12975` (16 chunks) and `10.1145_3095140.3095179` (5 pages, 11 chunks) — and have been since 2026-08-10, the same sweep that acquired the papers A-002c was built on | A-002d, read this session. **What the availability changes is the ticket, not just the citation.** A-002b was sized around Custodio's non-manifold remedy, a *grid* preprocessing pass that splits both cells at a shared ambiguous face's critical point (V-26). Grosso §6 says outright that this is not the only route: *"The algorithm proposed by Custodio et al. [CEPS13] … **do not really solve the manifold problem for case 13. Instead, they split neighbor cells** to avoid consistency problems across cells borders. Our method solves manifoldness across cells based on the asymptotic decider adding additional interior vertices."* So the ticket's constraint (b) — *"the non-manifold remedy is a grid change, not a table change"* — is true of one route and false of the other, and the route that avoids it keeps unambiguous cells on the existing table path. **The method rule is V-28's, one ticket later and still not learned:** an acquisition status recorded once was trusted for four days without being re-checked against the library that had meanwhile acquired it |
| V-30 | **Grosso's quadratic re-derived, and now agreeing three ways rather than two.** The coefficients of `a·u² + b·u + c`, whose roots are where two opposite faces' hyperbolas meet, were expanded here from `(i₀ − g₀)(g̃₁ − g̃₀) − (i₀ − g̃₀)(g₁ − g₀)` **before** being compared with anything: `a = (f₅−f₄)(f₀+f₃−f₁−f₂) − (f₁−f₀)(f₄+f₇−f₅−f₆)`, `b = (i₀−f₀)(f₄+f₇−f₅−f₆) − (f₁−f₀)(f₆−f₄) − (i₀−f₄)(f₀+f₃−f₁−f₂) + (f₅−f₄)(f₂−f₀)`, `c = (i₀−f₀)(f₆−f₄) − (i₀−f₄)(f₂−f₀)`. They agree term for term with the paper's eq. (4) **and** with the authors' own program, which writes them in a different but algebraically identical grouping | A-002d, `the_coefficients_reproduce_the_face_hyperbola_difference` — worst disagreement against a direct evaluation **1.33e-15** over 20,000 cells. V-24 had a derivation agreeing with a paper; this has a derivation, a paper and an implementation. **Also settled: Grosso's corner numbering is this crate's**, `v₀ = (0,0,0)`, `v₁ = (1,0,0)`, `v₂ = (0,1,0)`, `v₄ = (0,0,1)` … is exactly `cube.rs`'s `(i&1, (i>>1)&1, (i>>2)&1)`, so every formula transcribes with no permutation. That is a coincidence and not a design, so `grosso_corner_numbering_is_ours` pins it |
| V-31 | **The reference implementation both Grosso papers cite has been deleted from GitHub, and it survives — the paper's own listing of the inner hexagon is corrupt, so recovering it was load-bearing rather than tidy.** `github.com/rogrosso/tmc` 404s and the author's public repositories list only `infovis`. **Software Heritage has a full snapshot**: `swh:1:snp:8b33006176ec6abd716a2f259933f4f14f42b46d`, visited 2024-04-06, revision `swh:1:rev:b2048b9081bb2a097fbbcabf500784afa726490d`, from which `mc/MarchingCubes.cpp` (`swh:1:cnt:b77a45b1f2be63e877521746be225e5903fc45df`, 2,226 lines, MIT) reads out intact. The Wayback Machine has only the repository's root page — its subdirectory listings 404 — so it would have answered nothing | A-002d. **Three things the program settles that the prose does not.** **(a)** The paper's own hexagon listing, in the copy this project holds, assigns `p₁` and `p₂` the same triple, which cannot be a ring; the program's order is `(u₀,v₀,w₀), (u₀,v₀,w₁), (u₁,v₀,w₁), (u₁,v₁,w₁), (u₁,v₁,w₀), (u₀,v₁,w₀)` — consecutive vertices differing in exactly one coordinate, verified here as a closed axis-parallel ring over 1,146 hexagons. **(b)** *"Three quadratic equations"* is one quadratic: *"It is enough to compute a pair of solutions for one face. The other solutions are obtained by evaluating the equations for the common variable."* **(c)** Tunnel-versus-twelve-vertex-contour is **not** decided by the asymptote-side predicate of Proposition 1 and Corollary 1, whose precise form the prose never pins down. The program does not evaluate one; it branches on the number of closed contours. ~~**So the rule-5 stop A-002h was pre-emptively flagged for never had to be taken**~~ **— amended at A-020, and claim (c) was too strong (M-229).** Branching on the contour count is what the program does and it is *not* equivalent to Proposition 1: it calls a case-13 cell with contours of nine and three a tunnel, and that cell's inside region is two separate blobs rather than one cylinder. Those are exactly the cells whose contour edges span three inner-hexagon steps, which the construction has no rule for. **So the asymptote-side predicate is load-bearing after all**, and the reference's shortcut has a hole in it. The stop was still correctly not taken for A-002h's hexagon, which the program does settle — and the method rule is V-28's exactly: a deleted repository is not a missing source until an archive has been asked ~~**So the asymptote-side predicate is load-bearing after all**~~ **— amended again at A-020, and claim (c) was right the first time (M-230).** The predicate is derivable from the paper's own normal form and is now implemented (`BodySaddles::same_asymptote_side`); measured against the contour count over 400,000 random cells it agrees **exactly** — false precisely when there is one contour, which is Corollary 1. So the program's shortcut *is* equivalent to Proposition 1, and following it was sound. What it is not is the whole tunnel test: Corollary 6's length bound is what excludes the `[9,3]` cells, and neither Proposition 1 nor the contour count sees them. **The method note survives all three revisions**: the prose was called unpinnable without the derivation being attempted, and one derivation settled in an afternoon what two amendments had guessed at. |
| V-32 | **Two more rows of `meshing-library-target.md` marked `PAYWALL` are in the corpus, and the table's status code has now been wrong four times.** Nielson 2003 *On Marching Cubes* (`10.1109/tvcg.2003.1207437`) and Lopes & Brodlie 2003 (`10.1109/tvcg.2003.1175094`) both answer to a corpus query and both were marked as having no open-access route | A-020, checked while looking for a second source on the case-13 tunnel criterion. **This is V-29's finding repeating**: that entry corrected the two Grosso rows on 2026-08-14 after they had been indexed since 2026-08-10, and the same mistake was sitting two rows away. `PAYWALL` in that table records *what the resolver could not reach*, which is a statement about the tooling and reads as a statement about the library. **The doc now says so in a callout above the legend**, because correcting the individual rows is what was done last time and it did not stop the next pair. **The method rule is V-28's, one level up**: ask the artefact before believing the index that describes it — and if an index has misled you twice, fix the index rather than the entry |
| V-33 | **A `paper_download` that reports success can return a landing page, and this is the third producer of that signature.** `10.1142_s0218195912600060` (Attali, Lieutier & Salinas) reports downloaded and converts to French HAL UI chrome plus an abstract — no body. **Purge and re-fetch from `hal.science/hal-00785082/document`**, the direct PDF path | noted during the Phase 8–15 backlog build-out. **A second identifier to protect, and the failure mode is the opposite one:** Dey, Edelsbrunner, Guha & Nekhayev, *Topology preserving edge contraction* (1999) — the primary source for the link condition R-001 rests on — **has no DOI at all**. Semantic Scholar carries it with `doi: null`; it is free at emis.de. **Nothing may invent one.** A fabricated DOI is worse than a missing one because it resolves to something, and the citation then looks checked. **The method rule is V-28's, applied to acquisition rather than to archives**: a tool reporting success is not evidence the artefact arrived — open it and look |
| M-206 | **Two independently derived constructions locate the same body saddles, to 1.1e-12.** `interior::SweptFaces` sweeps a plane between two opposite faces and solves a quadratic in the sweep **height** `t` for the heights at which the plane's bilinear saddle sits on the surface. `trilinear::BodySaddles` intersects the two faces' hyperbolas, solving a quadratic in the face **coordinate** `u`, and interpolates the other two coordinates linearly. They share no arithmetic, no coefficient and no parametrisation. Over 200,000 random cells, 1,076 reached the comparison and the worst gap between the two answers was **1.13e-12** | A-002d, `the_body_saddle_heights_agree_with_the_swept_saddle_roots`. **This relationship was not predicted; it was measured, and the reasoning that preceded it was wrong.** The two were expected to differ, because a swept-plane saddle has `∂F/∂u = ∂F/∂v = 0` while a hexagon vertex is the crossing of two axis-parallel lines and so has a *different* pair of partials vanishing. They coincide anyway. Two supporting measurements from the same sweep: all **6,618** hexagon vertices lie on the level set to within **1.8e-13**, which is Grosso's equation (5) checked rather than assumed; and over 200,000 cells relabelled by a cyclic axis permutation there are **0** disagreements about the saddle count, which is M-204's rotation property holding for a construction that is deliberately asymmetric in the axes. **Six body saddles occur in 1,140 of 200,000 uniformly random cells (0.57%)** — the population A-002e has to find a field to reach |
| M-207 | **The reference implementation loses a root this one keeps, in two places, and both are the textbook quadratic formula.** `(−b ± √d)/2a` divides by `2a`, so where `a` is exactly zero the equation is linear with a real root `−c/b` and the formula returns infinities instead; and the reference guards with `d > 0`, so an exact tangency — a double root, one intersection point — is dropped rather than counted once | A-002d, `a_linear_equation_keeps_its_root` and `a_double_root_is_one_point_not_two`. Both matter downstream rather than here: Grosso §5.3 selects a cell's interior vertex by counting face pairs with a **single** solution, a count the textbook form cannot produce. This module solves the linear case as a linear equation and reports a tangency as **one** point, not two — the second guarding against a degenerate zero-area "hexagon" claiming six saddles. **Both fixtures are constructed, not searched, and the reason is the fixture trap in a new direction:** a bisection onto a tangency lands *near* the crossing, where the discriminant is a small non-zero number and two nearly-equal roots is the correct answer, so a searched fixture would have tested the branch next to the one it names. Constructing the coefficients backwards out of exact binary fractions (`a, b, c = 1, −1, ¼`, double root at `u = ½`) makes the discriminant exactly `+0.0`. **Also carried over from A-002c: Kahan's stable form rather than the textbook one**, for the reason `interior.rs` already records — `a` is a difference of near-equal products with nothing keeping it away from zero |
| M-214 | **The tunnel and the twelve-vertex contour are told apart by counting rings, and both are reachable: 2,053 and 173 in 396,877 random surface cells.** Six body saddles mean a cell is one or the other and the saddles cannot say which (M-206). Grosso states an asymptote-side criterion for the split (Proposition 1, Corollary 1) whose precise predicate the prose never pins down; the authors' implementation does not evaluate one, branching on the number of closed contours instead, and so does this | A-002f, `six_saddles_split_into_tunnels_and_twelve_vertex_contours`, both branches asserted reachable so neither is dead code. **0.52% tunnels and 0.044% twelve-vertex contours**, against Grosso's own 2,057 and 7 on a 512²×641 CT skull — the same order for tunnels and a ratio two orders apart for the rarer case, which is what uniformly random corner values buy over real data. **Two exhaustive results alongside it, over all 256 cases × 64 masks (16,384 combinations).** Every cut edge lies on exactly one closed ring and the rings follow the links at the wrap, so the walk is the link structure and not a re-derivation of it. And the ring lengths produced are exactly **{3, 4, 5, 6, 7, 8, 9, 12}** — the set the paper lists, with **ten and eleven absent**, which is the half worth noticing: an eleven-ring would leave one cut edge of a twelve-edge cell over, and the face walk's parity forbids it. At most **4** rings and longest **12**, both pinned so an over-sized buffer is as visible as an under-sized one |
| M-215 | **The `u` pair's two lines are crossed relative to the other two, and the inner hexagon is what proves it.** A line joining one pair of opposite faces is fixed by the other two coordinates, and the obvious reading is that solution `k` contributes the line at `(v[k], w[k])`. That is right for the `v` and `w` pairs and **wrong for `u`**: reading the hexagon's ring — consecutive vertices differ in exactly one coordinate, so a pair differing only in `u` *is* a `u`-line — gives `(v₀, w₁)` and `(v₁, w₀)`, indices `0,1` and `1,0` | A-002g, `the_line_inventory_agrees_with_the_inner_hexagon`, checked over **1,083** hexagons against the ring rather than against any transcription. The authors' implementation agrees (`fc3 = fs[1][0]*fs[2][1] + fs[1][1]*fs[2][0]`, crossed, against `fc1` and `fc2` uncrossed), **and the disagreement with the obvious reading is what made it worth deriving** — a copied crossed index looks like a typo to fix, and fixing it would have put the interior vertex off the surface in a third of cells |
| M-216 | **The interior vertex is a transcription, and the check on it is geometric: 149,803 of them, every one on the level set to 6.7e-12.** Grosso §5.3 says only *"if there are only two such lines, the additional vertex is the intersection point. If there are three such lines, the additional inner vertex is the midpoint between the two saddle points"* — which does not determine the four-line case at all, and does not say which of three lines' **three** pairwise intersections the "two saddle points" are. So the branch selection comes from the authors' program (V-31) | A-002g, `the_interior_vertex_lies_on_the_level_set`. **A transcription cannot be checked by re-reading it, only by a property it must have**: the point is a saddle of the trilinear interpolant, so it lies *on* the surface, and an index taken from the wrong axis lands it somewhere the interpolant is not zero. All three reachable branches fire — **99,777** two-line, **32,916** three-line, **17,110** four-line — and the five-and-six-line arm is separately asserted unreachable from the disk path over 397,666 cells, so its `None` is not a silent wrong answer |
| M-217 | **The disk path costs no new budget: worst case 12 triangles and 1 interior vertex per cell, against `MAX_TRIANGLES = 12` and `MAX_CENTROIDS = 3` already in place.** A ring of `k` fanned from an interior vertex costs `k` triangles against a vertex fan's `k − 2`, so the twelve-vertex ring goes 10 → 12 — which is exactly the bound A-015 already raised it to, for the same reason (a centroid fan emits one triangle per cycle edge) | A-002g, `the_worst_case_triangle_count_is_pinned`, exhaustive over all 16,384 `(case, mask)` pairs. The ticket asked for the per-cell vertex bound at `mod.rs:140-141` to be **re-derived rather than reused**, and the re-derivation is that §5.3 adds at most **one** cell-local vertex where A-015 budgeted three. ~~A-002h's tunnel needs exactly three, which is the same bound again — so the `u32` index-space check does not move for either~~ **— amended at A-002h and wrong in both halves (M-218).** Three inner vertices is Grosso *2016*'s collapse; **2017 §4.1 keeps all six**, and the reference implements 2017. Measured: a tunnel names **6** interior vertices and reaches **22** triangles, so A-002b must raise both `MAX_CENTROIDS` and `MAX_TRIANGLES`. The claim was made from the earlier paper without checking which one the construction being implemented came from |
| M-218 | **The twelve-vertex contour needs a closing step the tunnel does not, and a manifoldness test found it rather than a reading of the paper.** A tunnel is an *annulus*: its two contours pass the inner hexagon from opposite sides, so every hexagon edge is traversed twice and the patch closes itself. A twelve-vertex contour is a **disk** — one ring, circling the hexagon once — so it leaves the hexagon as a six-edge hole. Without a closure the cell emits a patch with a boundary no neighbour shares, which is a hole in the mesh | A-002h, found by `the_tunnel_patch_is_manifold_inside_the_cell` on its first run, not by reading §5.2 — which does say *"at the end, vertices of the inner hexagon are merged into three inner vertices which generates the last triangle"*, and which I had implemented as far as the ring walk and no further. **The winding cannot be fixed either**: a twelve-vertex ring circles the hexagon in either direction depending on the configuration, so the closing fan's orientation is derived from the hexagon edges the ring walk actually laid. The first two attempts were a fixed winding in each direction and the directed-edge check rejected both. **Measured: 2,297 tunnel and twelve-vertex patches, every one manifold *and* consistently wound inside the cell** — each non-contour edge traversed exactly once in each direction, each contour edge exactly once, because its second face belongs to the neighbouring cell. That is the property Chernyaev's tunnel triangulation fails and the entire reason this route needs no grid subdivision. Costs, pinned: **22** triangles and **6** interior vertices worst case, both above the budgets A-015 set |
| M-219 | **The reference implementation has a transcription-grade typo in the detached-ring test, and it is one line from six.** Finding which of three contours is *not* part of a tunnel means comparing the ring's span in `u` against the interval between the two `u` roots. The reference computes that span with `umax = (u_e2 > umax) ? u_e1 : umax` — comparing the third vertex and then storing the **second** | A-002h, `10.1111/cgf.12975`'s implementation, recovered per V-31. The other five lines of the same six-line block are written correctly, and a max that can return a value smaller than the one it just compared is not a max. Written correctly here. **The point is not the typo but what it says about the method**: this is the one file that settled three questions the papers left open (V-31), and it is also a file with a bug in it. Transcribing from a reference implementation and verifying by a geometric or topological property — as A-002g and A-002h both do — is the only combination that survives either source being wrong |
| M-220 | **The singular face is an artifact of quantised data, and this crate cannot reach it: 0 of 1,838 on eight reference fields and 0 of 299,215 over 400,000 random cells.** A face is singular when its bilinear saddle sits *exactly* on the level set — the two hyperbola branches degenerate into crossing straight lines, and the asymptotic decider's binary choice is between two answers that are both wrong because the surface passes through the saddle. Grosso 2017 exists for this case and Table 1 counts **8, 58 and 20** of them per 512²×~700 CT volume | A-002i's reachability, measured before implementing it, `how_often_a_face_is_singular`. **The difference is the data, not the code.** A singular face needs `v₀·v₂` and `v₁·v₃` to be *bit-identical* `f64`s. CT voxels are quantised integers and collide into that readily; a continuous signed-distance field essentially never does. **Two consequences.** A-002b is *not* blocked by A-002i, which is how the series was originally sequenced: a case no reference field can produce cannot change any mesh the acceptance criteria measure, so it cannot be a prerequisite for them. And the case is still worth implementing eventually, because **a consumer feeding quantised density — `u8` voxels, the normal thing in a voxel game — reaches it immediately**, which is exactly the audience this crate is for. The fixture will have to be constructed rather than sampled, as ✗22's was |
| M-221 | **`0 × NaN` is `NaN`, and it took the extractor to find that — not 400,000 random cells.** The interior-vertex construction selects among candidate coordinates by multiplying each by a `0`/`1` flag and summing, which is a selection *for finite values only*. A coordinate the mask marks as out of range **can be non-finite**: the linear solve for `v` or `w` divides by a difference that vanishes wherever the interpolant does not vary along that axis. So a zero weight does not suppress it | A-002b, found the first time `extract` ran the trilinear path, on `noise_cavity` at 17³ — mask `0b011010`, giving `[0.758, NaN, 0.347]`. **The reference dodges this and this crate deliberately did not**: it substitutes `−1` for every non-finite coordinate on creation, where `BodySaddles` documents the mask as the single authority on which numbers mean anything and stores whatever the arithmetic produced. Keeping that is right — a sentinel is a second representation of "absent" — but it makes the *consumers* responsible, and the selection has to be a branch rather than a product. **Why 400,000 random cells missed it**: `the_interior_vertex_lies_on_the_level_set` checks `trilinear(f, p) ≈ 0`, and a `NaN` coordinate makes that comparison `false`, so the assertion `worst < 1e-9` never fired on it — `NaN.max(x)` returns `x`. A test that tracks a **maximum** cannot see a `NaN`; the extractor's debug assertion could |
| M-222 | **χ falls by exactly two per tunnel and by nothing else — the interior rule's topology change is arithmetic, not approximate.** A tunnel is a handle and a handle costs a closed surface exactly two. Everything else the rule does is topology-neutral by construction: giving an ambiguous contour an interior vertex and fanning from it adds one vertex, three edges and two faces, and `1 − 3 + 2 = 0`. Measured across all eight reference fields at three resolutions, against a tunnel count taken from the classifier rather than from the mesh: **`noise_cavity` 17³ 3 tunnels χ 62 → 56; 25³ 4 tunnels 72 → 64; 33³ 2 tunnels 10 → 6**, and χ unchanged on all seven other fields | A-002b, `chi_falls_by_two_for_every_tunnel_and_by_nothing_else`. **The 33³ row is what makes this a test rather than a tautology.** That resolution has **four** six-saddle cells but only **two** tunnels: the other two are twelve-vertex contours, which have the same six saddles and are *disks*, not handles. A rule that treated all six-saddle cells alike would miss by exactly four, and the identity catches it. Confirms A-002f's ring-count discriminator end to end, from a direction it was not built from |
| M-223 | **The interior rule costs 1.95% at 33³ and 0.14% at 65³, on the only field that exercises it.** `noise_cavity`, f32, median of criterion's samples, Ryzen 9 5900X / Linux: decider **2.8492 ms → 2.9048 ms** at 33³ and **18.630 ms → 18.656 ms** at 65³ | A-002b, `cargo bench --bench extract -- noise_cavity`. **The cost falls as the grid refines, which is the same relationship M-209 records from the other end**: tunnels are an undersampling phenomenon, so a finer grid has proportionally fewer cells for the rule to do anything in — 3, 4, 4 six-saddle cells at 17/25/33 and none at all by 65³ on this field. The rule is confined to cells with an ambiguous face, so on the five reference fields that have none it is **byte-identical** to the decider, which the golden fixture pins: 15 of the 24 new rows share a hash with `marching_cubes+decider` and 9 differ, and the 9 are exactly `gyroid`, `fbm_terrain` and `noise_cavity` (M-40's two, plus the field added to reach the third case) |
| M-224 | **Manifold Dual Contouring's non-manifoldness has nothing to do with tunnels, which is what A-017 assumed, and it survives the correct face rule.** Three measurements on `noise_cavity`, each ruling out an explanation. **(a) Not tunnels.** Exactly **one** of the 30 offending edges at 17³ and one of the 64 at 33³ lies within `1.5h` of a tunnel cell, while **all** of them lie within `1.5h` of an *ambiguous* cell — the field has 193 and 502 ambiguous cells against 3 and 2 tunnels, so the correlation is with ambiguity in general at a rate of about 13–15% of ambiguous cells. **(b) Not duplication.** Every offending edge carries exactly **four faces and four *distinct* triangles**, uniformly — no threes, no fives, nothing emitted twice. Four distinct triangles on one edge is two sheets genuinely meeting, which is the opposite of ✗17's finding for Marching Cubes' fan chords, where the four faces were two cells each emitting the same triangle twice. **(c) Not the face pairing.** Under `FaceAmbiguity::AsymptoticDecider` the count falls from 30 to 8 at 17³ but does not reach zero — and the same setting **introduces** 3 offending edges on `gyroid` at 25³ where the default `Separate` gives none | A-017, `the_manifold_dual_contouring_defect_is_four_distinct_faces_on_one_edge`, pinned in both directions. **So A-017's own framing was wrong in both halves**: it asked whether this is a defect in this crate's construction *or* the published guarantee not covering a cell whose interior the interpolant joins, and the answer is neither, because interior topology is not involved. What is left is the quad walk around a crossed grid edge. **The gyroid row is the sharpest part** — it means no face rule is uniformly better, so this cannot be fixed by choosing one |
| M-225 | **A-017's mechanism, and the grid predicts the mesh defect exactly — 30, 64, 8 and 40, with zero error.** An ambiguous face has **all four** of its edges cut. Manifold Dual Contouring places one vertex per cycle per cell, so when all four of those edges lie in *one* cycle in each of the two cells sharing the face, all four dual quads — one per crossed grid edge — connect the **same pair** of cell vertices. Four quads on one dual edge, and a quad gives exactly one of its two triangles to each side, which is the *four distinct faces* M-224 measured. The identity is `non_manifold_edges == shared ambiguous faces whose four cut edges lie in one cycle on both sides`, computed from the **grid** with no mesh involved against a count the validator takes from the **mesh** with no grid involved | A-017, `the_defect_count_is_predicted_from_the_grid_alone`. Exact at 17³ and 33³ under **both** face rules: `Separate` 30 and 64, `AsymptoticDecider` 8 and 40. **This is a limit of one-vertex-per-cycle rather than a defect in this transcription.** Schaefer, Ju & Warren's argument separates sheets *within* a cell and nothing in it stops two different crossed edges of one shared face resolving to the same pair of cycles — which is only possible at all because an ambiguous face has four crossings rather than two. **It also kills the obvious fix**: the decider reduces the count without reaching zero and *raises* it on `gyroid`, so no choice of face rule makes the output manifold |
| M-226 | **Subgrid Marching Tetrahedra's output does not need welding and is damaged by it — the doc had been telling consumers to do the harmful thing.** The module said *"vertices are emitted per tetrahedron and are not shared … before welding, the output is a triangle soup: every edge looks like a boundary edge"*, and told callers to weld before validating. That was true when M-93 and M-96 were written; **A-014h ended it** by giving every crossing a global identity, and the doc did not follow. Measured on all eight reference fields at 17³: the **raw** output already has `boundary_edges == 0` on every closed field — it is a surface, not a soup — and the weld is a no-op on seven of them, identical vertex count and identical topology | A-018. **On `noise_cavity` it is worse than a no-op**: it merges exactly one pair of vertices that are coincident *by position* and distinct *by identity*, fusing two sheets and taking non-manifold edges from **288 to 290** and vertices from 392 to 395. Sharing by identity is strictly finer than sharing by position, so once identity sharing is complete the coarser rule can only contribute mistakes. **The general statement is the finding**: a positional weld is not a way to *make* a mesh a surface, it is a way to *join* two surfaces, and applying it to output that is already shared trades correctness for nothing. The subgrid validity suite now judges the extractor's own output, which is both more honest and two edges kinder; the weld-specific tests keep welding, because measuring the weld is what they are for |
| M-227 | **Orientation now reaches zero on every reference field at every resolution, and M-187's caveat was about the walk rather than the mesh.** M-213 had found `after <= before` was not a law — `noise_cavity` went 1,580 → **2,422** at 25³. The cause is propagation *crossing* a four-face edge: with no well-defined neighbour there, the flood fill takes whichever it reaches first and carries a **consistent** winding across a patch that is consistent with the wrong side. Stopping at such an edge instead — one `if neighbours.len() > 2 { continue; }` — leaves each sheet to be seeded and oriented on its own evidence | A-019, the flipped-edge census, now **all zeros**. **The result is larger than the ticket asked for.** It wanted the growth removed; what it got is the residue *eliminated*: `csg_difference` 6 → 0, `torus` 6 → 0, `thin_plate` 8 → 0, `gyroid` 138 → 0, `fbm_terrain` 19 → 0, `noise_cavity` 1,629 → 0 and 1,580 → 0. **So M-187's second half was a property of the walk, not of the mesh.** It said *"no assignment of windings can make four faces pairwise oppose across one edge"*, which is true and does not imply a residue — that edge only contributes one if the walk *visits* it, and `inconsistently_oriented_edges` counts runs of exactly two faces, which a four-face edge is not. The non-manifold edges are still there and still counted as non-manifold; what has gone is the orientation damage they were causing elsewhere |
| M-228 | **Grosso's tunnel triangulation has an undefined case, it is reachable, and 400,000 random cells could not reach it.** The rule closes each contour edge by how many steps its two endpoints are apart around the inner hexagon — one triangle for zero, two for one, three for two. **Three steps has no rule.** The paper gives none and the authors' implementation has no branch for it: its `switch` runs `case 0/1/2` with no `default`, so it silently emits nothing and leaves a hole | A-002h's gap, found at E-213 while searching for a presentable demo configuration. **The fixture trap, for the sixth time in this repository** (M-32, M-38, M-44, G-003, M-204, here): `the_tunnel_patch_is_manifold_inside_the_cell` asserted the three-step count was zero over 2,297 tunnel patches and passed, because uniform random corner values never produce the shape. **Where it lives:** Marching Cubes' **case 13** — the four alternating corners, the only case with all six faces ambiguous — at particular face resolutions, giving a tunnel whose two contours are **nine and three** vertices. That is also outside Grosso's **Corollary 6**, which says a tunnel's contours are at most six and three; the random sweep produces only `[3,3] [3,4] [3,5] [3,6] [4,4]` and `[3,3,6]`, all compliant, which is why the corollary looked verified. **The response is a refusal, not a repair.** `extract` returns `Error::UnresolvedTunnel` rather than emitting the hole — inventing the missing triangulation is exactly what rule 5 forbids, and a hole reachable only on case 13 is the kind of defect that arrives in a consumer's collider before anyone notices. A-020 owns deriving it |
| M-229 | **The contour-count discriminator misclassifies case 13, and the misclassified cells are exactly the ones with no triangulation rule.** `Contours::topology` calls a six-saddle cell a tunnel when it has two or more contours — which is what the authors' implementation branches on, and what V-31 recorded as making Proposition 1's asymptote-side predicate unnecessary. Measured by flood-filling the cell's inside region on a 96³ grid and counting how many components its inside *corners* fall into, a computation sharing nothing with the classifier: a genuine tunnel joins its same-signed corners through the interior and lands in **one** component, and both shipped tunnel fixtures do. The `[9,3]` case-13 cells land in **two** — their inside region is two separate blobs, which is not a cylinder | A-020, `a_nine_and_three_cell_is_not_one_connected_tunnel`. **Grosso's Corollary 6 was right and was being read as a description rather than a test**: it says a tunnel's contours are at most six and three, so a `[9,3]` cell is excluded from the tunnel case *by the corollary itself*, and the ring count admits it anyway. **The two findings are one finding.** M-228's undefined three-step edge and this misclassification occur on the same cells, and in that order: the cell is called a tunnel it is not, sent to a triangulation built for a cylinder, and that triangulation then meets a configuration it has no rule for. So A-020 is a **classification** problem first — which is cheap, and is where its ticket now says to start — and a triangulation problem only if some genuinely-tunnel cell also reaches three steps, which nothing has yet shown |
| M-208 | **Interior ambiguity was unreachable by this crate's entire test suite: 0 of 68,385 reference-field surface cells have six body saddles.** Swept over all seven reference fields at 17³, 33³ and 65³ with A-002d's classifier. `sphere`, `box_exact` and `thin_plate` reach **zero** saddles in *every* surface cell; `torus` and `fbm_terrain` top out at three; `gyroid` and `csg_difference` reach five, five times between them in the whole sweep. Six — the configuration that carries a tunnel or a twelve-vertex contour — never happens | A-002e, `how_often_the_reference_fields_reach_six_body_saddles`. **This is M-40 one level deeper and much stronger.** That entry found the ambiguous *face* absent from five of seven fields; this finds the ambiguous *cell* absent from all seven. The consequence is what matters: every gate the A-002 series could have leaned on — the T-001 validity sweep, the golden fixture, the property suite — was structurally incapable of exercising the code the series exists to write, and a per-cell proof was all it could ever have got. **Smooth analytic solids do not produce tunnels because they are too smooth**: a tunnel needs the field to reverse twice across one cell, and a sphere never does. Both papers behind MC33's corrections used randomly generated scalar fields for precisely this reason |
| M-230 | **Corollary 6 is the tunnel test; Proposition 1 is not, and the derivation that shows it also validates itself.** Grosso's Proposition 1 splits tunnel from twelve-vertex contour by whether both roots of the quadratic lie *at the same side of the corresponding asymptotes*, which V-31 recorded as unpinnable prose. It is not: §3 gives the face bilinear the normal form `G = α + η(s−s_c)(t−t_c)` and names the asymptotes, and matching coefficients against `a + (b−a)s + (c−a)t + (a+d−b−c)st` gives `η·s_c = a − c` directly. Because a point on the hyperbola has `η(s−s_c)(t−t_c)` fixed and non-zero, one coordinate per face decides it, and the division-free form `η·s − (a−c)` needs no guard for an untwisted face | A-020, `the_asymptote_side_predicate_is_corollary_1`, `corollary_6s_length_bound_agrees_with_the_flood_fill`, `the_tunnel_contour_shapes_are_pinned`. **Three results, and the first is what makes the other two usable.** **(a) The derivation is right, graded by something that shares no arithmetic with it.** Over 400,000 random cells the predicate is false *exactly* when the cell has one contour — zero disagreements — which is Corollary 1 stated the other way round. A derived predicate and a face-segment walk agreeing cell-for-cell is evidence in the shape V-30 established. **(b) It is therefore not what A-020 needed.** Every multi-contour shape passes it, the `[9,3]` cells included, so Proposition 1 admits exactly the cells Corollary 6 excludes and the plan that expected it to separate them was wrong. The fix is the corollary's length bound read contrapositively — a contour past six cannot belong to a tunnel — which is one comparison against the existing `Contours::longest`. **(c) Half of Corollary 6 is false.** It reads *"one of the contours can have at most 6 vertices and the other 3 vertices"*; the measured tunnel shapes are `[3,3] [3,4] [3,5] [3,6] [4,4]` and `[3,3,6]`, and `[4,4]` has no three-vertex contour at all. Only the bound is used, and only as a necessary condition. **The independent grade on the bound is M-229's flood fill, restricted to two-contour cells** — a three-contour tunnel's detached ring caps its own corner group and so *correctly* gives two components, and a twelve-vertex contour gives one while being a disk, so the oracle is decisive nowhere else. On two-contour cells the bound and the flood fill agree on all 400 sampled. **The three-step hole M-228 found is now unreachable from the tunnel path** and its error is kept as a live guard rather than deleted, because "nothing has reached it" is a sample, not a proof |
| M-231 | **The `[9,3]` cell is not a topological subcase; it is a singular face that the strict interior test lets through.** A-020 classified a six-saddle cell with a contour past Corollary 6's bound as `SeparateDisks` and made `extract` refuse it, and A-020b was written to derive the triangulation it needs. **There is nothing to derive.** Measured over 2,000,000 random cells per sweep: **continuous corner values produce zero such cells** among 11,354 six-saddle ones, and **every one that quantised values produce has a body saddle within `1e-12` of a cell face** — 261 of 261 across quanta of 0.1, 0.25, 0.05 and 1/3. A saddle *on* a face is Grosso 2017 §4.2's singular configuration, which `has_inner_hexagon`'s strict `0 < x < 1` test admits because the arithmetic misses the face by a few ulps | A-020b, `every_separate_disks_cell_has_a_saddle_on_a_face`. **The control is what makes it a finding rather than an observation about rounding.** The same degeneracy among six-saddle cells that are *not* `SeparateDisks` runs at 8.0% at quantum 0.05 and 79.4% at 1/3 — it swings by an order of magnitude with the quantum while the `SeparateDisks` rate stays pinned at 100%. So the association is with this configuration, not with quantisation in general. **Two things it costs.** A-020b is blocked on A-002i and will most likely be closed by it. And A-002i's own reachability claim needs reading carefully: it measured **0 of 1,838** ambiguous faces singular on the reference fields by a *bit-exact* product comparison, which is true and narrower than it sounds — 86–100% of these cells have a bit-exact singular face and the rest are the same configuration a rounding away, so the exact test undercounts the phenomenon it is testing for. **The method note is the fixture trap again, in its eighth appearance, and from the other direction**: the shipped `[9,3]` fixtures were found by a search over *rounded* corner values, and rounding is precisely what creates the degeneracy. A fixture found by search carries the search's bias, and here the bias *was* the phenomenon |
| M-232 | **A singular face is unreachable from continuous data and routine from quantised data — and the rate at `u8` density matches Grosso's CT counts.** A-002i measured **0 of 1,838** ambiguous faces on the eight reference fields and **0 of 299,215** over 400,000 random cells, and read that as the case being rare. It is not rare; it is *conditional*, and the condition is the one this crate's audience meets. Quantise the same 400,000 random cells and the count goes to **6,658** at a quantum of 0.1 and **20** at 1/255 — `u8` density — against **0** unquantised | A-002i, `a_singular_face_needs_quantised_data`, and the pre-existing `how_often_a_face_is_singular` census is unchanged at 0 and 0 because it now calls the crate's own `singular_face_mask` rather than a copy of the predicate, so the two cannot drift. **20 per 400,000 cells is the number to compare**: Grosso 2017 Table 1 counts 8, 58 and 20 singular faces per 512²×~700 CT volume (tier V), the same order, from data quantised the same way. **The method note is the one M-208 earned, applied to a rate rather than a property**: a count of zero is a statement about the sample, and the sample here was continuous `f64` while the stated audience feeds integers. Nothing was wrong with the original measurement; it was narrower than the conclusion drawn from it |
| M-233 | **A-002i's blocker is not the vertex cache, it is that a singular face needs a third routing and the resolution mask has two bits of room.** The BLOCKED note said the work was *"a new cache keyed on faces, not a branch in a decider"*. The cache half is right and is determined — Grosso 2017 §4.2 says *"three saddle points will lie on a singular face, but only **one** will be shared with the neighbor cell"*, so one slot per grid face suffices, and a grid face is named by min-corner sample plus normal axis exactly as an edge is named by lower sample plus direction axis. **The other half was missed.** Definition 3.2 requires that a singular face *not* divide the surface into two branches, and on a singular ambiguous face all four edges are cut with the level set being the two crossing asymptotes — so the four cut edges must meet at the hyperbola **centre**, a four-valent junction. `segment_links` takes `joined` as **one bit per face**: exactly two routings exist, and both are permutations of the cut edges, asserted over all **384** (case, ambiguous face, bit) combinations with no exceptions | A-002i, `a_face_has_only_two_routings_so_a_singular_face_has_nowhere_to_go`. **So the blocker is representational and sits in `Contours`**, which the whole trilinear path and A-002's 16,384-pair decider validation rest on — a larger and more delicate change than a second cache, and the reason the ticket's size is still wrong. **A third obstacle is a rule-5 stop rather than work**: §4.2's fewer-than-six-saddle arm is fully specified, but its six-saddle arm says *"the other two points are **slightly moved** towards the interior of the cell"* and gives no distance. The recovered reference is the 2016 code (V-31), whose singular handling is the face-pairing choice of the 2016 paper rather than §4.2's construction — verified by reading it — so **no artefact supplies the constant**. **The method note**: the ticket's blocker had been written from the *symptom* (two cells emit coincident vertices) rather than from the construction, and a blocker written from a symptom under-describes the work. Reading §4.2 and Definition 3.2 against the code took an hour and moved the ticket from `M` to something that needs splitting |
| M-234 | **A-017 closed by decision rather than by code, and the limit is stated as an identity so it cannot rot.** Manifold Dual Contouring's headline property — *"the entry that takes the non-manifold count to zero"* — was conditional on the seven original reference fields being unable to produce an interior ambiguity (M-208). The user chose option 1 of three on 2026-08-15: document the limit, keep the censuses, stop making the unqualified claim | A-017. **Options 2 and 3 are recorded as rejected, not deferred**: splitting the cycle needs a second vertex where Schaefer, Ju and Warren give one and stops being their algorithm, and taking the dual of the trilinear surface changes what A-010 claims. **The documentation is placed where the claim was made, not in a notes file** — the module header states the precondition beside the guarantee, the crate README's tradeoff row carries the exception, and A-010's archive row is amended. **The method note is what makes this cheap rather than lossy**: the limit is written as the identity `non_manifold_edges == shared ambiguous faces whose four cut edges lie in one cycle on both sides`, which is asserted grid-against-mesh with zero error under both face rules (30/64 and 8/40). A documented limit that is also an executable prediction fails loudly if it is ever wrong, which a paragraph of prose does not |
| M-235 | **The example capture rig renders headlessly but cannot control its own window size, and it fails silently — which is the guarantee it was written to provide.** `examples/common/mod.rs`'s `size_window` reads `ISOMESH_WINDOW=WIDTHxHEIGHT` and its doc comment says it exists so *"two captures of the same example"* do not *"come back different shapes"*. Measured on a bare X server with no window manager: `1280x720` and `1600x900` both produce **exactly 836×1356**, portrait, with no error logged | E-214. **The half that works is worth recording too**: `ISOMESH_SCREENSHOT` produced a correct 264 KB still of `marching_cubes_tunnel` — HUD, contours, gold hexagon — because it goes through Bevy's own screenshot path and reads back from the GPU rather than scraping the root window. **So the rig needs a compositor for geometry and not for rendering**, which is the opposite of what one would guess, and it is why an `Xvfb` display fixes it while the whole thing already works over an unmapped window. **The method note is the silent-default one**: `size_window` returns early on a missing variable and logs on a malformed one, but has no branch for *applied and ignored* — the case that actually occurred. A setting that cannot verify it took effect is a setting that reports success when it did nothing. ~~So the rig needs a compositor for geometry~~ **— wrong, and corrected at E-214 the same day (M-239).** The window was always resizable; `size_window` merely ran in `PreStartup`, before `bevy_winit` creates the OS window, so the write landed on an entity and was overwritten. No compositor and no `Xvfb` were ever needed |
| M-236 | **The shootout's header claimed the wrong counts for both of the things it was counting, and nothing could have caught it.** `benches/shootout.rs` printed *"seven reference fields, five algorithms"* while the run enumerated **eight and seven** — `noise_cavity` landed at A-002e and the `marching_cubes+decider` row was added later, and neither edit reached the string | X-001. **This is the failure mode the ticket predicted, found in the ticket's own subject.** X-001's premise was that hand-enumeration makes adding an algorithm an `O(N)` edit; the sharper cost is that the edit *fails silently* — a bench that measures seven algorithms while its header says five compiles, runs, and writes a CSV nobody doubts. Both counts now derive from the lists. **The registry's order was checked against the committed evidence rather than assumed**: `docs/measurements/shootout.csv`'s algorithm column matches `ALL_EXTRACTORS` element for element, which is what says the conversion preserved the file that ✗14, M-19 through M-22 and O-11 all quote |
| M-237 | **The QEF buys 2× accuracy on smooth fields and 100× on sharp ones, and pays for all of it in self-intersections — measured with one algorithm and one rule swapped.** X-002's seam makes `DualContouring` generic over its `VertexRule`, so `Qef` and `Centroid` run through the *same* cell classification, quad walk, sign conventions and buffers. Symmetric Hausdorff at 65³, as a ratio of QEF to centroid: **sphere 0.486, torus 0.457, csg_difference 0.255, box_exact 0.010, thin_plate 0.010** — a factor of two where the surface is smooth and a factor of a hundred where it has a feature, which is precisely what A-007 exists for. Self-intersections per 1,000 triangles at 33³: QEF **3.118** on gyroid, **13.837** on fbm_terrain, **29.745** on noise_cavity; centroid **0.000** on all three | X-002, `benches/ablation.rs`, `docs/measurements/ablation.csv`. **Both halves were pre-registered in the bench's own `report` before the run** and both hold. **The controlled part is what makes it a measurement of the rule.** Vertex counts, triangle counts and non-manifold-edge counts are *identical* between the two arms on every field and resolution — gyroid 33/33 and 15/15 and 69/69, noise_cavity 131/131 and 297/297 and 322/322 — so the rule provably does not reach the topology. `the_ablation_arms_differ_only_in_position` pins the same thing behaviourally: 680 vertices, byte-identical index buffers, and **all 680 positions different**. A comparison between `SurfaceNets` and `DualContouring` cannot say this, because those are two structs with two `extract` methods and any difference between them is a difference between two implementations |
| M-238 | **Probabilistic quadrics are this crate's existing solve with a different regularizer, and the derivation is the finding — X-004's premise that they *supersede* it is falsified.** Trettner & Kobbelt (`10.1111/cgf.13933`) §3.1 gives the probabilistic plane quadric under Gaussian normal noise as `A = Σnᵢnᵢᵀ + N·Σₙ`, `b = Σnᵢnᵢᵀqᵢ + Σₙ·Σqᵢ`. This crate already solves in **centroid-relative** coordinates, so writing `x = c + Δ` and `qᵢ = c + rᵢ` with isotropic `Σₙ = σ²I` turns that system into `(M + Nσ²I)Δ = Σnᵢdᵢ + σ²Σrᵢ` — and **`Σrᵢ ≡ 0`, because `c` is defined as the arithmetic mean of the `qᵢ`**. The term the paper adds to `b` vanishes and what remains is `solve_with` at `λ = Nσ²` | X-004, `the_probabilistic_quadric_is_the_existing_solve`. **Checked numerically, not left as algebra**: a direct assembly of the paper's equations in world coordinates, sharing no line with the crate's solve, agrees to **1.110e-16 over 296 cells**. **So no `ProbabilisticQuadric` solver was written** — it would be a second execution path computing numbers the existing path already computes, which the one-path rule forbids. What shipped is the part that genuinely differs: the regularizer scales with the **crossing count** where `Qef` applies one fixed λ to every cell. **Measured over all eight fields at 65³** (`docs/measurements/ablation.csv`): Hausdorff ratio of scaled to fixed is **1.0000 sphere, 0.9957 torus, 0.9992 csg_difference, 0.7519 box_exact, 0.7519 thin_plate** — no worse anywhere, **25% better on both sharp-feature fields**, and identical on both of them to four figures, which is itself unexplained. Self-intersections per 1k at 33³ fall too: gyroid **3.118 → 2.551**, fbm_terrain **13.837 → 13.571**, noise_cavity **29.745 → 28.749**. **The paper's headline does not apply here and that is worth saying**: its selling point is robustness without an SVD, and this crate's solve never used one — it is a 3×3 adjugate (✗16). Its real novelty is *anisotropic* `Σₙ`, which needs a noise model that analytic fields with exact gradients do not have |
| M-239 | **A window can resize itself without a window manager, and two obvious fixes both fail silently before the third works.** E-214 was blocked on the belief that controlling window size needs a compositor (M-235). It does not. The fault was *when* the request was made | E-214, `examples/common/mod.rs`. **Three attempts, and the first two are the instructive part.** (1) `PreStartup` — the OS window does not exist yet, so the write lands on an entity and is overwritten by whatever the window comes back as: `1280x720` and `1600x900` both gave **836×1356**. (2) `Update`, latching once the window reports the requested size back — **also 836×1356**, because `Window::resolution` reads back the value *this system wrote* rather than what the platform granted, so it latched on its own echo. (3) `Update`, re-applying across the first 30 frames, then stopping — **1280x720 → 1493×840 and 1600x900 → 1866×1050**, both exactly the request at the display's 1.1666 scale factor, and 1493/840 is 16:9 to four figures. **The method note is that a blocker is a hypothesis.** M-235 recorded "needs a compositor for geometry" from two measurements that were both consistent with it and with the true cause, and the ticket sat blocked on a system package nobody needed to install. The cheap disproof — move the system and re-run — took one attempt |
| M-240 | **The crate could subtract and intersect but not union, and nothing noticed because no reference field unions anything.** `Difference` exists because `csg_difference` needs it and `Intersection` exists because `Gyroid` needs capping — every combinator in the crate was there because a *fixture* asked for it, so the most basic CSG operation was the one missing. Found by trying to write an authoring demo, which is the first thing that ever wanted to put two primitives together | E-216. **The fixture trap again, in its ninth appearance and from a new direction**: the previous eight were properties that held because no fixture could falsify them (M-208's five at once, M-228's three-step edge, M-231's `[9,3]` cells). This is the same mechanism applied to *API surface* rather than to a claim — not a wrong belief but an absent capability, invisible for exactly as long as nothing exercised it. **`Union` is also the safe direction and that is worth having written down**: `min` of two 1-Lipschitz functions never overestimates distance, so a sphere tracer stepping by it can only under-step, where the `max` in `Intersection` and `Difference` overestimates near concave seams — the direction that lets a tracer step through a surface. That asymmetry is what Phase 11's `F-001` gives a type |
| M-241 | **Two demo GIFs were shipped that did not show what their captions claimed, and a single frame is not an inspection.** Both were checked by reading one still and both passed that check while failing their purpose | E-214, E-216, raised by the user looking at them. **`marching_cubes_tunnel` had no `Mesh3d` at all** — every triangle was three `gizmos.line` calls, with a per-triangle normal tick added *"so the two sides read as surfaces rather than as a cage of lines"*, which concedes the problem in a comment. Interactively it is a good diagram, because you can see through it to the hexagon and the contours; as a picture it fails, because the claim is about **surfaces** — two discs against one cylinder — and an outline shows no surface. Filled, the tunnel's hole is unmistakable. **`sdf_authoring` never turned its own knob**: 79 frames at a fixed `k = 0.12`, advertising a blend sweep it never performed, when `Capture::taken` is exposed *specifically* so a sweep runs in step with the capture. **The method note is that a GIF has a failure mode a still cannot expose**, and it is the one that matters: a recording of something that does not change is a recording of nothing. Check the first frame against the last, or check the parameter the caption names |
| M-242 | **The shootout's 112 rows came back structurally identical after X-001 rewrote how it enumerates, which is the refactor's independent check.** Converting `benches/shootout.rs` from a hand-written `Algorithm` enum to the `for_each_extractor!` registry touched every measurement the file produces. Re-run and compared row for row against the previously committed CSV: **112 rows, zero differences in vertices, triangles or non-manifold edges** | T-011. **This is what the baseline is for and it earned its keep on the first run.** A refactor of the code that *drives* the measurements is exactly the change golden hashes cannot see — T-007 pins meshes from the extractors, not from the bench — and "the diff looked mechanical" is not evidence. **The regression gate's tolerances follow the same split**: vertices, triangles and non-manifold edges are compared **exactly**, because they are deterministic and a tolerance there could only hide a real change; `hausdorff` gets 2% and `self_intersections_per_1k` 5% for cross-architecture float accumulation; `median_ms` gets 60%, which is a tripwire for a doubling rather than a benchmark. Baselines are per machine and the machine is in the filename, so a run on another host finds none and says so instead of comparing against numbers that never applied to it |
| M-243 | **An unquoted heredoc executed the words in this script's own Python comments, and it failed silently for a value of "failed" that exits zero.** `scripts/regress.sh` embeds Python in a `<<PYEOF` heredoc. Unquoted, bash expands the body — so a backtick in a docstring becomes command substitution, and a comment reading *"it found no row with a `median_ms`"* made bash run `median_ms` and print `command not found` on **every invocation** | T-014, found while wiring provenance into baselines. **The gate still passed.** Command substitution failing does not fail the script; the noise went to stderr and the exit code was zero, so a CI step would have stayed green while printing an error on every run — the kind of thing that survives for months because it looks like someone else's warning. **Fixed structurally rather than by escaping**: the heredoc is quoted and the four values it needs come through the environment, which removes the whole class rather than the instance. **The method note is that prose inside a heredoc is code**, and this repository writes unusually long comments, so it is unusually exposed to it |
| M-244 | **A declared Lipschitz constant was wrong on the first run, and the test that caught it was written in the same commit by the same person.** F-001 replaced `is_exact_distance() -> bool` with a declared `FieldBound`, and the first draft gave `noise_cavity` `Lipschitz { l: 2.598 }`, reasoned loosely from Perlin noise's gradient bound. `every_field_meets_the_bound_it_declares` measured `|∇f|` reaching **7.734** and failed it immediately | F-001. **The defect is identical to the one the ticket existed to remove** — `csg_difference` declaring `true` with `// away from the seam` beside it — reproduced by the person removing it, one screen away, in the same sitting. A declaration is not evidence about the thing declared, however careful the declarer. **The gyroid's constant was guessed too and survived only by luck**: `√3` was declared, the sampled maximum came out **1.695**, and the test passed. Deriving it properly gives `|∂g/∂x| = |cos a cos b − sin c sin a| ≤ 2`, hence `|∇g| ≤ 2√3 ≈ 3.464` — twice the guess, and the guess would have been an *under*-declared bound, which is the dangerous direction: a sphere tracer dividing by a constant smaller than the true one steps through the surface. **`noise_cavity` is now `Unbounded` rather than measured**, because a sampled maximum is a lower bound on a supremum and declaring one as a Lipschitz constant would be unsound in exactly that direction. F-002 is where a real constant gets established |
| M-245 | **The eikonal condition cannot tell a CSG underestimate from a true distance, which is why a bound needs two numbers rather than one.** `csg_difference` is `max(box, −sphere)` and measures **100.0% eikonal** — `‖∇f‖` within 5% of one on every sample — while its values are *not* distances near a concave seam. Away from the seam the active operand is an exact distance, so `‖∇f‖ = 1` there; the seam is a measure-zero set a grid essentially never lands on | F-002, `every_reference_field_meets_its_declared_bound`. **So `|∇f| ≈ 1` is necessary and not sufficient**, and the obvious single-number validator — *check the field is eikonal* — would have passed the exact field this whole phase exists because of. `FieldBound` carries `l` for how fast the field changes and `q` for how far its value is from the distance precisely because the first cannot detect a failure of the second. **The census across all eight fields**, `sup ‖∇f‖` and eikonal fraction: sphere/torus/box_exact/thin_plate/csg_difference all **1.000 and 100%**; gyroid **1.727, 67.0%**; fbm_terrain **2.850, 11.7%**; noise_cavity **7.748, 85.6%**. Note noise_cavity's pairing — a high eikonal fraction *and* a gradient reaching 7.7 — which is a second way the fraction alone misleads |
| M-246 | **`min` is not better than `max`; each is exact in one region and wrong in the other, by the same amount.** F-003 predicted an asymmetry — that composing two exact fields with `min` *"yields a declared bound F-002 confirms"* while `max` yields *"a strictly weaker one"* — and asked for a test asserting it. Measured against closed-form ground truth on two overlapping unit-ish balls, over 64,000 sample points: **union is exact *outside* (0.000e0) and wrong inside by 6.989e-1; intersection is exact *inside* (0.000e0) and wrong outside by 6.995e-1** | F-003, `min_and_max_are_both_inexact_and_the_asymmetry_is_by_region`. **A ratio of 0.999 — they are mirror images.** The prediction was directionally wrong in an interesting way: there *is* an asymmetry and it is by **region**, not by operator. `min` gets the outside right because the nearest surface of a union is the nearest surface of one of its parts; `max` gets the inside right for the dual reason. Neither operator is better than the other; they are wrong about opposite halves of space. **The test asserts the ratio stays within 0.5–2.0**, so if one operator ever does become strictly weaker the prediction is revived rather than silently forgotten |
| M-247 | **Repeated CSG destroys the worst case and leaves the typical case untouched, and a renderer barely notices what a precision bound calls ruin.** F-004 asked how fast a field stops being a distance under repeated subtraction. Carving a box with `N` seeded random spheres and measuring the empirical underestimate ratio `q̂ = |f| / d_ray` over 13,824 points: worst case **0.5774 → 0.1815 → 0.0726 → 0.0040** at 0, 4, 32 and 256 strokes — **143× down** — while the **median stays at 1.0000 at every count**, and mean sphere-tracing steps go **5.2 → 5.5 → 9.9** | F-004, `benches/csg_degradation.rs`, `docs/measurements/csg_degradation.csv`. **So "the field is no longer a distance" and "the field is no longer usable" are different questions with different answers**: the precision bound is wrecked at a handful of strokes and the practical cost has not quite doubled at 256. A destructible game should watch the tracing cost, not the bound. **The eikonal column is the control and it reads 100.0% at every stroke count** — confirming M-245 operationally, and confirming that F-004's own proposal to sample `‖∇f‖` would have produced a flat line and concluded nothing degrades. **One caveat is in the measurement rather than the field**: `q̂` uses a ray distance, and a ray from a box corner leaves along the diagonal, so an **uncarved** box reports `q̂_min = 0.5774 ≈ 1/√3`. That is this metric's floor, not a defect, and the curve is read against it |
| M-248 | **One field evaluation replaces 576, and buys between 1.1× and 11.8× depending on how much of the volume the surface reaches.** F-005 applies Hart's bound to subgrid Marching Tetrahedra: with Lipschitz constant `l`, a cell whose centre satisfies `|f| > l·(√3/2)·h` cannot contain the surface, so one sample settles what M-98 measured as 576 evaluations. At 33³ with 16 samples per edge, cost against Marching Cubes **before → after**: sphere 204× → 36.7×, torus 294× → 45.6×, box_exact 381× → 88.9×, csg_difference 400× → 91.4×, thin_plate 389× → **32.9×**, gyroid 449× → 295× | F-005, `empty_cell_rejection_is_measured_per_field`. **The spread is the finding, and it is geometric.** `thin_plate` gains **11.8×** because a plate leaves almost all of its domain empty; `gyroid` gains **1.5×** because a triply periodic surface reaches nearly every cell, *and* because its Lipschitz constant is 3.46 rather than 1, which inflates the rejection radius by the same factor and disqualifies cells a distance field would have rejected. **So the optimisation pays in proportion to emptiness and is taxed by a loose constant** — which is a second reason F-001's constants had to be derived rather than guessed, beyond soundness. **Bit-identical output is the safety property and it is checked, not argued**: a rejection test that is ever wrong produces a *hole*, and a hole is invisible to every validity gate — the mesh is simply missing a piece and remains perfectly manifold |
| M-249 | **A directional Lipschitz bound buys exactly nothing on five of six fields, and 1.80× on the sixth — the null result F-006 asked for.** Galin et al.'s segment tracing marches by a bound computed along the ray rather than over all directions. **The paper states the condition for its own failure**: *"when the implicit objects have an almost uniform distribution of primitives and a uniform Lipschitz bound over their support Ω, the benefit is limited or negative in terms of speed."* Five of this crate's fields are exactly that — `sphere`, `torus`, `box_exact`, `csg_difference` and `thin_plate` are 1-Lipschitz everywhere, so no directional bound can be smaller than 1 | F-006, `a_directional_bound_helps_only_where_the_global_one_is_loose`. **Measured, sphere-tracing steps global → directional**: 228 → 228, 360 → 360, 108 → 108, 108 → 108, 148 → 148 — **identical, asserted as equality rather than as a small difference**. The gyroid is the one field that can gain, and its directional bound is *derived* rather than sampled: along a coordinate axis the directional derivative is a single partial, `|∂g/∂x| ≤ 2`, against the global `|∇g| ≤ 2√3`. **Steps 2184 → 1213, a gain of 1.80× from a bound 1.73× tighter** — slightly more than the ratio, because fewer steps also accumulate less conservatism. **So the technique is worth exactly what the field's global bound is loose by**, which for an analytic distance primitive is nothing, and this crate's tracing costs are dominated by fields whose constants are already tight |
| M-250 | **Refining an edge crossing on the real field helps *curved* fields by 13–15% and does nothing at all for the CSG one — the reverse of F-007's prediction, and for a reason the ticket did not consider.** The ticket argued that `min`/`max` kink the field along an edge crossing a seam, so linear interpolation misses the root, and set its acceptance as *"`csg_difference`'s Hausdorff improves"*. The reasoning about the kink is correct; the conclusion about the measurement is not | F-007, `refinement_helps_curved_fields_and_not_the_kinked_one`. **Measured at 33³ with 24 bisection steps, Hausdorff before → after**: sphere **5.361e-3 → 4.561e-3 (0.851×)**, torus **1.379e-2 → 1.201e-2 (0.871×)**, box_exact **1.443e-1 → 1.443e-1 (1.000×)**, thin_plate **8.927e-2 → 8.927e-2 (1.000×)**, csg_difference **1.515e-1 → 1.515e-1 (1.001×)**. **The mechanism: `csg_difference` is `max(box, −sphere)` and a box is planar**, so its field is *exactly* linear along an axis-aligned edge and there is nothing for a root-finder to find. Its Hausdorff is dominated by the box's own edges and corners, which Marching Cubes cannot represent at any resolution — that is A-007's problem and no amount of refinement touches it. **The fields that improve are the ones whose error genuinely is interpolation error**: a sphere and a torus are curved, so `f` along an edge is a curve and a straight line through its endpoints really does miss. **The method note is that "where the field is kinked" and "where the error is" were assumed to be the same place and are not** — the ticket reasoned from the defect's mechanism to the metric without checking what the metric was measuring |
| M-251 | **The exact distance transform agrees with brute force to the last bit, and is exactly one sample spacing from the analytic sphere — which is the sampling limit, not an error.** Felzenszwalb & Huttenlocher's separable `O(n)` transform, checked two independent ways: against an `O(n²)` search sharing no line of code with it (**1,287 samples, exact equality**, not a tolerance), and against the closed-form sphere (**worst 0.10000 against a spacing of 0.1 = 1.0000 cells**) | S-001, `agrees_with_brute_force_exactly`, `matches_the_analytic_sphere_within_one_spacing`. **The two checks answer different questions and both are needed**: brute force proves the algorithm computes what it intends — the lower-envelope pass is subtle and its failure mode is plausible numbers — and the analytic comparison proves the intent is a distance field. **The one-cell gap is structural.** The transform answers with the distance to the nearest opposite-signed *sample*, and the surface lies between samples, so a point whose nearest crossing falls mid-cell is off by up to a full spacing. Landing exactly on the bound rather than inside it is the expected result, and the test compares against `1 + ε` for that reason rather than for slack. **The grid is deliberately not a cube** — 11×9×13 — because an axis transposition survives a cube and dies here |
| M-252 | **Fast sweeping beats the exact transform everywhere, including where it was predicted to lose — because the seeding, not the solver, is where the accuracy is.** S-002's docs asserted that sweeping *"accumulates: a value ten cells away is the result of ten Godunov updates, each a first-order approximation"*, and would therefore lose at distance. Measured on a sphere at 41³ against the analytic field, worst error by band: **within two cells of the surface, swept 0.0333 against the exact transform's 0.1000; beyond eight cells, swept 0.0933 against 0.1000** | S-002, `sweeping_and_the_exact_transform_trade_places_with_distance`. **Sweeping wins by 3× near the surface and still wins, narrowly, far from it.** The mechanism is the seed rather than the sweep: the exact transform answers with the distance to the nearest opposite-signed *sample* and is therefore quantised to the grid, floor **one full spacing** (M-251), while sweeping seeds from the *interpolated* crossing and starts sub-cell. A sphere's characteristics are radial straight lines, the eight-orthant sweep follows them, and that head start survives to the edge of the domain. **The prediction was mine and was written into a doc comment as though established** — the same failure M-244 and M-250 record, which is three times in one day that a stated expectation went into prose before a measurement went into a test. The assertion is now "does not lose", so a field that flips the ordering fails loudly rather than quietly vindicating the original claim |
| M-209 | **The eighth reference field is searched, not designed — 97 of 610 candidates qualify, and the one taken is the gentlest rather than the richest.** The criterion was three-fold and all of it measured: the meshed result must be **closed** at 17³, 25³ and 33³, it must reach six body saddles at **all three** rather than at a lucky one, and extraction must not trip the zero-gradient assertion. `noise_cavity` — one octave of gradient noise at frequency `3.45`, iso `0.25`, capped to a sphere of radius `1.5` over `[-2, 2]³` — reaches **3, 4, 4** six-saddle cells and thins to **zero** by 65³ | A-002e, `the_tunnel_field_actually_contains_tunnels`, pinned in both directions. **97 qualifying combinations is the important number**: this is a plateau, not a knife edge, unlike the gyroid at high iso where the same search found the property appearing and vanishing between `iso = 1.44` and `1.45`. **The gentlest was taken deliberately.** Frequency `4.9` reaches four times as many six-saddle cells (4/13/16, still 3 at 65³) and roughly **quadruples** Subgrid Marching Tetrahedra's damage — 713 non-manifold and 5,803 flipped edges at 17³ against `3.45`'s 290 and 1,629. That extra is collateral from sampling noise near its feature size, not coverage of the interior rule, and the ticket asked for a field that reaches the configuration rather than one that maximises it. **The undersampling is not a defect but the mechanism**: Grosso reports refining a volume with 16 tunnels twice leaves one, and it reproduces here from the other direction |
| M-210 | **The natural level set of gradient noise is its own lattice, which is degenerate — so a noise reference field may not use `iso = 0`.** Perlin noise is *exactly* `0.000e0` at every point of its integer lattice, measured at six lattice points against `1.5e-1`-ish just off it. The zero level set therefore **contains the whole lattice**, and the surface is pinned through a regular array of the noise's own critical points | A-002e. The first version of `noise_cavity` used `iso = 0` on the grounds that it was the field's own balanced surface and no constant would be doing hidden work — an argument that is aesthetically right and geometrically wrong. It tripped `unit_gradient`'s zero-gradient assertion during extraction. Any non-zero level avoids the lattice entirely. **The method note is that the reasoning was backwards**: "zero is the least arbitrary choice" is a statement about the *number*, and what mattered was a property of the *generator* that one `println!` of six samples settled |
| M-211 | **A-010's zero was conditional, and nobody knew on what: Manifold Dual Contouring is not manifold on a field with interior ambiguity.** On `noise_cavity` it returns **30 non-manifold edges / 60 vertices at 17³** and **64 / 128 at 33³**, where all seven original fields give zero at every resolution. Its Euler characteristic parts company with Marching Cubes' too — `−30` against `−96` at 17³ — which is the sharper half, because χ agreement was P-5, *pre-registered before A-010 ran* | A-002e, pinned as `MDC_NON_MANIFOLD_CENSUS` and `MDC_CHI_CENSUS`, owned by **A-017**. **Note the exact factor of two in every row**: each offending edge carries exactly two offending vertices, its own endpoints, so the defect is a property of edges and the vertex count is a corollary. **Why P-5 was always going to fail here rather than merely happening to**: the dual identity `V' = F, E' = E, F' = V` needs the two methods to be describing the same surface, and on a cell whose interior is ambiguous they are not — Marching Cubes separates a tunnel that the trilinear interpolant joins. The prediction was sound and its precondition was invisible while no field could reach the configuration (M-208) |
| M-212 | **Welding two coincident vertices can *create* a non-manifold edge, which a test comment had asserted was impossible.** `subgrid/extract/tests.rs` reasoned that "sharing by identity is a subset of sharing by position, so every counter must be at least as good after welding as before". True of the sharing, false of the counters: merging two vertices that are geometrically coincident and topologically **distinct** fuses two sheets, which is exactly the pinch A-010 exists to avoid. On `noise_cavity` it costs **one** merged pair — 5,567 → 5,566 — and that pair adds **2 non-manifold edges and 3 non-manifold vertices** | A-002e, owned by **A-018**. **One pair explains three separate failures**, which is how it was identified as one mechanism rather than three: the over-merge itself; P-7 limb (a)'s weld plateau, `[5567, 5567, 5566]` across the four tolerances T-009 merged; and A-014h's completeness claim, where the exact identity rule removes **0** of the weld's 1. That last one is the confirmation rather than a fourth problem — the exact rule leaves the pair alone *correctly*, because they are different crossings, so the disagreement points at the weld and not at the identity rule. **The margin is inverted, not small**: the pair merges at exactly the `1e-4` policy and stays separate at every finer tolerance, so the usual "how far to the nearest disagreement" ratio comes out **below one**. A relaxed floor would have reported a direct hit as a near miss |
| M-213 | **Orientation can *raise* the flipped-edge count, and `after <= before` was not a law.** Propagation was expected to only ever agree edges it reaches, so the residue could shrink but not grow. On `noise_cavity` it shrinks at 17³ (**1,629 → 1,015**) and grows at 25³ (**1,580 → 2,422**) and 33³ (**1,477 → 3,341**) | A-002e, pinned in the flipped-edge census, owned by **A-019**. The mechanism is the same four-face edge M-187 already names as the residue's cause, reached often enough to change what propagation does rather than merely what it cannot fix: with 318 non-manifold edges the flood fill crosses one, commits to one side's winding, and carries a *consistent* orientation across a whole patch that is consistent with the wrong neighbour. Nothing is more wrong afterwards — a different and larger set of edges is now the one that disagrees. **M-187's law survives exactly where it was stated** (zero non-manifold edges ⟹ zero residue); what fails is the monotonicity that had been assumed alongside it and never separately tested |
| M-166 | **The interior test does not inherit `ambiguity`'s crack-freeness argument, and the reason is one operator.** The face decider reduces to comparing two *products*, and IEEE multiplication is The face decider reduces to comparing two *products*, and IEEE multiplication is commutative, so two cells reading a shared face in rotated order get bit-identical answers. The interior test's denominator is `((A + C) − B) − D` — a fixed subtraction order that a rotation **permutes** — and floating-point addition is not associative | A-002c, found while writing `the_test_is_deterministic_and_independent_of_diagonal_order`, which originally asserted bit-identical agreement and could not justify it. So two cells meeting on a face could evaluate the same sweep to different bits and disagree about a tunnel. **M-32's caveat in a second place**: equal by algebra is not equal by IEEE, and only the same expression is safe. The test now asserts agreement to a tolerance and says why; establishing the exact-agreement property is A-002b's, before any of this is wired into an extractor |
| ✗22 | **Appendix A's counterexample could not be recovered from the converted PDF, and fitting the numbers until they matched would have been rule 5's exact failure.** The OCR lists eight `(_0, _1)` value pairs with their letters scrambled, and no grouping of them reproduces the paper's own `c = A₀C₀ − B₀D₀ = 0.1701` | A-002c. Two candidate groupings give 0.0417 and 0.0598; a third that lands within 0.0002 requires reading `D₀ = −0.1075` where the text says `D₀ = −0.1692`. **That last one is the trap** — it is close enough to feel like a win and is arrived at by searching for the answer, which is precisely how a wrong case table gets into a repository. The fixture used instead is **derived**: `lo = [0.1, −2, 10, −2]`, `hi = [−10, 2, −0.1, 2]`, both faces ambiguous, `F` convex and negative at both ends so negative throughout, Δ running +14.1 → −14.1. It has the same sign *structure* as Appendix A's — which is itself corroboration, since that structure was not chosen to match but fell out of requiring a pole |
| M-161 | **§3.2.3's immersion is real, measured on five of seven fields — and it is not the defect.** Coincident polygons, counted by exact position over the unwelded soup with every triangle traced to the tetrahedron that emitted it: **`box_exact` 30, `csg_difference` 33, `thin_plate` 4, `fbm_terrain` 4, `gyroid` 1; `sphere` and `torus` none** | A-014d's unblocking measurement, `which_polygons_coincide_across_a_shared_face`. **`box_exact` carries the most coincidence of any field and validates `(0, 0, 0)`.** If duplication were the defect those two columns would track each other and they invert. Two further contradictions of the section as this repo restated it: **every coincident pair is wound the *same* way, not oppositely** — A-014e imposes winding from the gradient at the triangle's centroid, and two coincident triangles share a centroid, so the pair survives the weld as a doubled face rather than as the two sides of a sheet — and **three pairs on `box_exact`, three on `csg_difference` and two on `fbm_terrain` come from tetrahedra that do not share a face at all**, which §3.2.3's framing ("two adjacent tetrahedra sharing a face `f`") does not cover |
| M-162 | **A-014d's blocking question is answered, and the answer is no: the inset needs neighbour information, and the neighbours are not all in this cell.** Of `csg_difference`'s 33 coincident polygons, **27 have other tetrahedra standing on their boundary edges — 312 triangles in total, up to 12 on a single polygon — and 150 of those 312 (48%) belong to a *different cell*.** `box_exact`: 348 foreign users, 168 cross-cell, the same 48% | A-014d, same instrument. §3.2.3 moves the midpoints of polygon edges lying in an edge of the shared face; every triangle standing on those edges must move with it or the disk detaches, which is precisely what M-101's two reverted attempts measured. **So the condition is not expressible per tetrahedron — and a per-cell rewrite would not reach it either**, because half the affected triangles are in another cell. That makes the inset an architectural change to the extractor, on a number rather than on the suspicion the ticket recorded. **Two fields escape**: `thin_plate`'s 4 pairs and `fbm_terrain`'s 4 have **zero** foreign edge users, so those could be inset locally — the obstruction is not universal, which is worse than if it were, since a rule that works on some fields and tears on others is the failure mode this repo's one-path rule exists to prevent |
| M-163 | **`csg_difference`'s three surviving non-manifold edges are each 4 faces but only 3 *distinct* polygons, from exactly 2 tetrahedra in one cell.** `([11, 11, 7], tets 4 & 5)`, `([11, 7, 11], 2 & 3)`, `([7, 11, 11], 0 & 1)` — the three symmetric images of one feature. The soup carries **6** bad edges; the 3 that are duplication-only do **not** survive the weld and the 3 with three distinct polygons do | A-014d, `the_surviving_non_manifold_edges_are_not_duplicated_polygons`, matched by exact position from the welded mesh's `validate_features` back to the soup. **This is the one row where §3.2.3 is the right remedy:** removing the duplicated pair from such an edge leaves exactly 2 faces, which is manifold. It is also the row that shows deduplication would *not* be enough — dropping one copy leaves 3 faces and the edge is still bad. The two tetrahedra being siblings in one cell is what makes it look locally fixable; M-162 is why it is not |
| M-164 | **`gyroid`'s 138 inconsistently-oriented edges are not A-014d's, and the ticket's stated mechanism reaches 8% of them.** 186 triangles stand on those edges: **0 zero-area, 15 with the gradient lying in the triangle's own plane, 171 with a decisive vote** | A-014d, `what_gyroids_flipped_edges_are_standing_on`. The ticket claims insetting *"gives those polygons area and a normal transverse to the gradient, which fixes both symptoms at once."* The first half is real — those 15 have no answer to A-014e's `dot(face_normal, gradient(centroid))` and a transverse normal would settle them. The other **171 have an unambiguous answer and disagree with their neighbour anyway**: two triangles of one sheet, each correctly oriented against the gradient at its *own* centroid, wound opposite ways. That is A-014e's per-triangle vote having no consistency guarantee between neighbours, not §3.2.3's immersion, and no inset reaches it — `gyroid` has **1** coincident polygon in 8,088 triangles, so there is almost nothing there to inset. Moved to A-014f. (Zero-area reads 0 rather than the soup's 12 because the weld removes a triangle with a repeated vertex before the validator sees it) |
| M-159 | **The last four bytes cost 0.033 ms to move and 0.375 ms to wait for — because `poll(Wait)` drains every dispatch queued before it, not just the copy.** Breaking `extract_buffers` down at 129³: **total 0.454 ms, count pass 0.017, scan-plus-four-byte-read-back 0.408** — while a **bare** four-byte read-back with nothing queued measures **0.033 ms** | GPU-010b's premise, checked before building it. **So "remove the read-back" was the wrong framing and the ticket's own wording invited it.** The bytes are free; what costs is the *synchronisation point*, and what it waits for is the GPU doing work the caller asked for. Removing it does not delete that work — it stops the CPU standing still through it. Same discipline as GPU-012: the premise was checked and it turned out to be about something else |
| M-160 | **A zero-synchronisation extraction makes CPU time independent of grid size: flat at ~0.17 ms from 33³ to 129³, where the synchronising path grows to 0.722 ms — 4.2× at the top and widening.** Median of seven, warmed, same field and grid for both | GPU-010b. Measured: **33³ 0.170 → 0.167, 65³ 0.245 → 0.206, 129³ 0.722 → 0.171**. **The flatness is the result, not the ratio.** `extract_indirect` never waits, so what remains is CPU *recording* — buffer creation, bind groups, submits — which does not depend on how many cells there are. The synchronising path's cost is the GPU's execution time surfaced at the first wait, so it tracks the workload. **The GPU still does the work, and saying otherwise would be M-149's mistake again**: a frame is no shorter, the CPU is simply free during it. That is worth little to a demo extracting one surface and a great deal to a caller meshing many chunks, where the stalls serialise what could overlap. **What it costs is stated as a contract rather than hidden:** without the total the geometry buffers are sized from a *budget*. Worst case would be `MAX_TRIANGLES` = 12 per cell — **906 MB per buffer at 129³**, 1.8 GB for the pair, against a measured 38,456 triangles, a factor of **190×**. So the budget is the caller's, a surface exceeding it is truncated, and `IndirectGeometry::total` holds the real count for a caller who chooses to check — the same shape as `collider::readiness` being a check the caller runs rather than a promise the extractor makes |
| M-139 | **The wgpu pin holds, and the number it holds at is 29.0.4 rather than the 29.0.3 written down.** `CLAUDE.md` requires `wgpu` to match Bevy 0.19 exactly, on the grounds that Cargo resolves two majors side by side with no resolution error. Checked rather than assumed, against the two independent lockfiles this repo keeps: the root workspace and `bevy_isomesh` both resolve **`wgpu 29.0.4` and `wgpu-types 29.0.4`** — the same patch in both, from `"29.0.3"` requirements that are caret ranges on both sides | GPU-001. So the pin is doing its job and the *documented* version is one patch stale; the thing that matters is that the two agree, and they do. **Five API details differ from what the architecture doc's era would suggest**, all found by the compiler rather than by reading: `PollType::Wait` is a struct variant carrying `submission_index` and `timeout`, `InstanceDescriptor` has no `Default` (it has `new_without_display_handle_from_env`), `Instance::new` takes the descriptor by value, `DeviceDescriptor` gained a required `experimental_features` field, and `Limits::max_storage_buffer_binding_size` is `u64` where an `AdapterReport` written from memory would type it `u32` |
| M-140 | **The GPU this repository's numbers will be measured on, and the two limits that bound a dispatch.** `headless::Gpu` on the Ryzen 9 5900X box reports **NVIDIA GeForce RTX 3090, Vulkan backend, DiscreteGpu**, with **`max_storage_buffer_binding_size` 2,147,483,644 bytes** and **`max_compute_workgroups_per_dimension` 65,535** | GPU-001, `cargo test -p isomesh-gpu -- --nocapture`. The storage limit is the one GPU-004 has to size against: at 4 bytes per `f32` sample it caps a single binding at **536,870,911 samples, about 812³** — so a 1024³ grid does not fit in one binding on this adapter and will need splitting, which is worth knowing before writing the dispatch rather than after. **The adapter type is asserted, not merely printed:** a test fails if `device_type` is `Cpu`, because `force_fallback_adapter: false` is a promise that a software rasteriser would silently break, and a lavapipe run reports timings orders of magnitude off while looking merely slow |
| M-137 | **Paint that lives in the field does not move when the geometry under it is destroyed — drift is exactly 0.000000, not small — and the L²-nearest attribute transfer the research prices this feature on is not cheap here, it is unnecessary.** `docs/research/2026-08-11-novel-gameplay-opportunities.md:16` builds row 4 on *"common subdivision + L²-nearest attribute transfer"* (Integer Coordinates §3.8, §4.4), on the reasoning that the shared tet grid is the common refinement so the expensive half comes free. That machinery moves attributes from an old mesh to a new one and is needed exactly when the attributes live **on the mesh**. This crate's world is a base field plus an ordered edit log, so paint goes *in the log* and is a function of world position: a carve moves the surface, and the paint was never on the surface. Measured over a scripted run of 40 sprays and 2 carves, 320 recorded world points: **drift 0.000000 on both carves**, with **20 probes losing their surface entirely** (320 → 308 → 300 still within a cell of solid) | E-208, `bevy_isomesh/examples/game_paint.rs`. **The control is the load-bearing half, and the first version of this measurement did not have one.** The script originally sprayed a single colour, so overpainting red with red was indistinguishable from paint staying put and the instrument reported 0.000000 either way — a number that could not have come out otherwise. Cycling the palette makes it sensitive: **27 of 40 sprays now register non-zero drift, up to 0.886**, and the two carves still read exactly zero. A transfer-based implementation could only ever report a tolerance; this reports an equality, and there is a core test asserting the same thing bit-for-bit in `paint::tests` |
| M-138 | **The price of exact paint is that every sample walks the edit log, and it is sub-linear in practice: 2.33x the cost per chunk for 40x the log.** Median milliseconds per re-meshed chunk over a scripted 42-step run, binned by how many edits were in the log at the time: **0.1375 (edits 1–10), 0.1748 (11–20), 0.2450 (21–30), 0.3200 (31–40)** | E-208. Same shape as `game_dig`'s measured 3.7x-for-7x on brushes alone (E-202), and the mechanism is the same — the walk is real cost and not the dominant one at these lengths. **Two design choices keep it this cheap and both are asserted by test.** `PaintStack`'s `Sdf::sample` skips sprays *without evaluating their shapes*, so meshing a painted world costs one `match` arm more than meshing an unpainted one and not a field evaluation more — pinned by a test that a spray does not move the surface, bit-for-bit. And a log containing no sprays samples **bit-identically** to the equivalent `BrushStack`, because both fold through the same `brush::apply`; that test exists so the painted and unpainted paths cannot silently diverge, which is the failure rule 1 is actually about |
| M-135 | **The contour is 29% of a usable mesh, not 54% — and the largest stage is the *collider check* at 45%.** Seven reference fields at 33³ and 65³, `f64`, one process, median of seven timed runs after two warmups. Mean share per stage: **contour 29.0%, weld 25.5%, collider 45.0%, normals 0.4%**. The published comparison M-003 cites has contouring at 68 ms against halfedge construction at 58 ms — 54% — so this repository's own ratio is **half** that, and **optimising the extractor alone can buy at most 1.41x** on the whole job | M-003, `docs/measurements/stage_breakdown.csv`. **The normals stage is effectively free (0.4%)**, which matters because it is the one every extractor here already does for nothing: all seven reference fields override `Sdf::gradient`, so the `AreaWeightedFaces` pass measured is the cost a consumer pays only when their field has no analytic gradient — and even then it is noise. **The collider figure needs one qualification and keeps its force after it:** what is timed is `collider::readiness`, which *validates*, not a physics engine building its own structure. A shipping game may validate once rather than per chunk. But G-005 makes welded-manifold-correctly-wound the contract, and the cost of checking that contract is the single biggest line in the budget |
| M-136 | **"What fraction is the contour" has no single answer — it is 13.1% to 74.3% across seven fields, and the variable is how expensive the *field* is to sample, not anything about the mesher.** Ordered: `sphere` 13.1%, `box_exact` 15.1%, `csg_difference` 16.1%, `torus` 16.4%, `gyroid` 17.5% (all at 33³), rising through `thin_plate` 33.0% and 44.0% to **`fbm_terrain` 65.2% and 74.3%**. A 5.7x spread on the same extractor, the same post-passes and the same machine | M-003. **`fbm_terrain` is the outlier because fBm noise is expensive per sample**, so its contour swamps the fixed post-processing, while a cheap analytic field leaves the post-passes dominant. **Consequence for every optimisation decision here:** the ratio is a property of the *workload*, so "is it worth optimising the contour" cannot be answered in general and must be answered per field. A game on procedural terrain is in `fbm_terrain`'s regime and should optimise the extractor; one on authored CSG is in `sphere`'s and should not. The share also **rises with resolution on every field** — the contour is `O(n³)` while the post-passes scale with triangle count, which is `O(n²)` |
| M-134 | **M-21's negative intercept reproduces in *direction* on a second machine and a third of the range, and its magnitude is 40x smaller — so the sign is the finding and the value never was.** Fitting `t = a + b·n³` live over 17³–89³ on a Ryzen 9 5900X, median of three timings per point, eight independent runs: **Surface Nets' `a` is negative in 7 of 8** (range −0.133 … +0.009), Marching Cubes positive in 7 of 8 (−0.020 … 0.155), Dual Contouring positive in **8 of 8** (0.042 … 0.171). Every value is under **0.18 ms in absolute value**, against M-21's committed **−3.13 ms** for Surface Nets over a sweep reaching 256³ | M-002. **Both halves matter.** M-62 concluded from the committed CSV that there is no meaningful fixed cost; this reaches the same conclusion from a different machine, a shorter range and a different timing method, which is as close to independent replication as this repository gets. And M-21's rule — *"a physically impossible fitted parameter is the model telling you it is wrong"* — survives being measured at 1/40th the magnitude, because what carries the meaning is the sign, not the number. **`r²` inverts, which is worth knowing:** over this short range Marching Cubes fits *worst* (0.909–0.996) and the dual methods best (0.9986–0.99997), the reverse of M-21's full-sweep ordering, because MC's timings are smaller and therefore noisier relative to themselves |
| M-132 | **Subgrid Marching Tetrahedra *does* tile across a chunk boundary — 0 open edges in 20 configurations — and the 1 that said otherwise was a wrong bound in my own instrument.** Two adjacent chunks, welded, counting boundary edges in the shared plane: **0 on both fields, and 0 across six field phases at three sampling resolutions each**. The mechanism holds up in the source: `TETS[t]` is ordered by inclusion, so a tet edge always runs from the lower **cube-corner index** to the higher, which makes traversal a property of the grid. **M-79's warning is about a different renumbering** — it says a mesh renumbering vertices *per tet* would crack every shared face, and chunking renumbers per *chunk* while leaving the within-cell corner order untouched | B-007. **The phantom edge was the block's own `z` wall.** A chunk here is 16 cells at 0.25 — **4.0 deep** — and the exclusion list gave `z` an upper bound of **8.0**, so that wall was never excluded and one of its edges landed in the seam plane. **Tenth instance this session, and the second in the same helper**: B-006 fixed this exclusion once already, for omitting `y` entirely, and left `z`'s bound wrong in the same edit. Both hand-written versions were wrong in different ways, so the bounds are now *derived from the layout* and there is nothing left to write out |
| M-133 | **The dual methods are not *reliably* seam-closing rather than reliably open, which is a weaker and more useful claim than the one first published.** With the bounds corrected, `surface_nets` measures **5 on `waves` and 0 on `blobs`**, and `dual_contouring` **4 and 1**. So a gap depends on how a particular surface meets the seam, and either method can come out clean on a given field — what is missing is the *guarantee*, not the geometry on any one frame | B-007. `ChunkSeams::Gapped` now says exactly that, because "leaves gaps at every chunk boundary" — which is what B-006 shipped — is falsified by its own second field, and a consumer reading it would have concluded the opposite of the truth on half their scenes |
| M-130 | **On a *concave* edge Dual Contouring's advantage is 3.6x in the mean and only 1.56x in the worst case — nothing like the 101x E-104 measures on a convex corner.** A block with a quarter-space subtracted leaves a reflex dihedral whose position is known in closed form; the cutter orbits so the edge sweeps the grid continuously, and the distance from the exact edge to the nearest vertex is measured every frame at 41³ over thousands of frames. Mean / worst, in cells: **dual contouring 0.0798 / 0.3481, surface nets 0.2895 / 0.5426, marching cubes 0.3412 / 0.5528**. So the sharp-feature method is decisively better *typically* and converges toward the others at the extreme | E-209. **This is why the sweep matters.** E-104 measured one static configuration and needed a rule about which resolutions to skip to defend it; a moving feature is unaligned almost everywhere and aligned occasionally, and reporting the worst over the sweep is the only way that occasional case gets counted. A single-position measurement of a sharp feature is a measurement of that position |
| M-131 | **The cell clamp does not bind on a concave edge either, which extends M-28 rather than limiting it — and the prediction that it would was registered first and falsified.** The hypothesis, written into E-209's HUD before the run: a convex corner's solution is interior to its own cell so the clamp never binds (M-28), while a reflex edge's solution "wants to sit outside" and would therefore be capped by it. Measured: clamped and unclamped Dual Contouring are **identical** on this fixture — mean 0.0798 against 0.0801, worst **0.3481 against 0.3481**. The reasoning was wrong in a simple way: a reflex edge passing *through* a cell has its QEF solution inside that same cell, so there is nothing for the clamp to catch | E-209. **What does cap the worst case is deliberately left unestablished.** It coincides with the edge lying near a sample plane — E-104's alignment trap arriving continuously rather than at chosen resolutions — but this demo does not isolate it, and the HUD says so rather than asserting a mechanism it did not test |
| M-128 | **The dual methods do not tile across chunk boundaries and Marching Cubes does — measured on two fields, after the defect shipped in a README GIF.** Two adjacent chunks meshed independently and welded, 16 cells at 0.25, counting boundary edges lying in the shared plane and excluding the two-chunk block's own **six** walls. `waves` / `blobs`: **`marching_cubes` 0 / 0, `surface_nets` 5 / 0, `dual_contouring` 4 / 1, `subgrid(4)` 0 / 0** *(corrected at B-007 — the figures first published here were inflated by a wrong `z` bound, see M-132)*. The mechanism is the one-vertex-per-cell rule: Marching Cubes places vertices on grid **edges**, which two neighbouring chunks compute from identical corner values and therefore agree on exactly, while a dual method places one vertex per cell **interior** and a boundary quad needs the *neighbour's* vertex, which the chunk does not have | E-210, found because the user looked at the hero GIF and said the gaps were broken. **It shipped.** The showcase was switched to `DualContouring` during an unrelated aliasing investigation and never switched back, so the README's lead image was a cracked world for one commit. Switching to `MarchingCubes` removed every slash. **Consequence:** `IsomeshPlugin` exposes `Extractor::DualContouring` and `Extractor::SurfaceNets` on a chunked `VoxelVolume` with no warning, and choosing either silently produces a cracked world — B-006. Every existing streaming demo uses Marching Cubes, which is why nothing had caught it |
| M-126 | **`BrushOp::commutes_with` is sound *and tight* — on eight overlapping brushes it called all seven adjacent pairs correctly, with zero unsound and zero conservative answers.** Audited by swapping every adjacent pair in an edit log, re-folding the world and comparing an FNV-1a hash over every position and index against what the predicate promised: **5 pairs said commute and the mesh was bit-identical** (M-36's case), **2 said no and the mesh moved** (M-37's case), **0 promised commutation and moved** — the only genuinely wrong answer — and **0 refused a commutation that actually held**. Undo/redo round-trips bit-identically over the same hash, and the check asserts the undone state *differed* so it cannot pass on a no-op | E-207. **The predicate answers on operations alone and cannot see shapes, which is where its conservatism lives.** The first fixture moved each brush ~3 units against radii of ~1, so every pair was disjoint — and disjoint brushes commute whatever their ops are, so the two add/subtract boundaries came back **conservative**: `commutes_with` said no, the mesh was unchanged. **This refines M-37 rather than contradicting it:** order matters across an add/subtract boundary *when the shapes overlap*, and M-37's 40,320-ordering sweep was on overlapping brushes |
| M-127 | **A fixture built to demonstrate that order matters reported 0 of 7 pairs where order mattered.** E-207's first scripted sequence moved each brush's centre by roughly 3 units between ops while the radii were about 1, so consecutive brushes never touched and every add/subtract boundary commuted anyway. The audit was working perfectly and had nothing to find. Fixed by lowering the path frequencies until consecutive centres sit ~1.5 apart, inside the sum of the radii, at which point 2 of 7 pairs stop commuting | E-207. **Eighth instance this session of a fixture that could not exhibit the property it was chosen for** (M-107, M-112, M-113, M-115, M-117, M-119, M-120's control), and the second where the *demo of a defect* was the thing that could not show the defect. Part 5's rule — search for a fixture that exhibits the property rather than picking one that looks like it should — is now the most-cited rule in the file |
| M-124 | **The amortized cost per frame is the number a game spends, and it tracks the budget to within one chunk across a 320x range of budgets.** 288 chunks re-meshed under `DirtySet::mesh_within_budget`, budget swept from 25 us to 8 ms, 2,360 frames each. Chunks per frame **1.00, 1.02, 3.49, 15.15, 28.80, 96.02**; mean cost per frame **0.085, 0.087, 0.251, 1.016, 1.905, 6.269 ms** against budgets of 0.025, 0.05, 0.2, 1, 2 and 8 ms. Above a few chunks the mean lands within 5% under the budget asked for, which is the figure a frame actually pays and the one no paper reports — papers report throughput, and a game cannot spend throughput | E-206. **The control is what makes the flat line mean anything:** the identical queue through `mesh_dirty` costs **20.62 ms in a single frame**, missing a 16.7 ms frame outright, against a budgeted **peak of 2.10 ms** — a **9.8x lower peak for exactly the same total work**. A budget does not make meshing cheaper; it decides which frame pays |
| M-125 | **The never-livelock guarantee costs one chunk, exactly, and below one chunk's cost it degrades to precisely 1.00 chunks per frame rather than to zero.** `mesh_within_budget` consults its predicate *after* each chunk, which its doc prices as "overshooting by at most one chunk". Measured: one chunk here is ~0.072 ms, and the worst overshoot past the budget is **0.157, 0.172, 0.123, 0.099, 0.117, 0.122 ms** at budgets from 0.025 to 8 ms — bounded by one chunk's cost and **flat across a 320x budget range**, exactly as claimed. At a 25 us budget, a third of one chunk, the rate is **1.00 chunks per frame**: never zero, never two | E-206. **Both halves of the doc's argument are now numbers rather than assertions.** The guarantee is real, its price is one chunk, and the alternative it rejects — a queue that cannot progress while it grows — never happens even when the budget is a third of the smallest indivisible unit of work |
| M-123 | **Aliasing on a composed field came from the *field*, not the extractor — and swapping Surface Nets for Dual Contouring is what proved it.** Intersecting a terrain half-space with a thickened gyroid, `max(p.y − height, |g| − t)`, staircases badly along the rim where the two surfaces cross. The obvious suspect is the extractor rounding a crease, so Surface Nets was replaced with Dual Contouring — the crate's sharp-feature method — and **the artefact was unchanged**, which rules the extractor out in one run. Two field-side fixes removed it: a **smooth** intersection (`smooth_max(a, b) = −smooth_min(−a, −b)`, blend `0.7`) rounds the knife edge to wider than a cell, and raising the shell from `0.55` to `0.85` takes the sheet from **2.2 cells to 3.4 cells** thick | E-210. **This is M-72's aliasing arriving from the other direction.** M-72 measured Marching Cubes disintegrating a feature near the sampling limit and read it as a property of the extractor; the same failure is reachable purely by *composing* two smooth fields into one whose interesting region is sub-cell. A hard `max` is not a smooth field, and nothing warns you: both operands are exact distance functions and the result is not. **Rule of thumb this earns:** a CSG intersection's rim is a feature, and it has to clear the sampling limit like any other |
| M-120 | **Transvoxel's mirrored seam works, and it had never been run — E-107 only ever meshed the fine block on the low-`x` side.** A camera flying out *and back* along an LOD ladder puts coarse blocks on both sides of itself, so half the seams are `inset_boundary(face_bit(0, 1))` with transition cells sampled at the fine block's **first** `x` index instead of its last. Measured across a full out-and-back over 12 blocks at levels 0–2: **0 open edges on both sides at every position**, with up to 2 seams below the camera and 2 above. The mirror is the classic place for an inside-out winding, which no manifold or Euler check can see — so this is the configuration most likely to have been silently wrong, and it is not | E-205. Proven by a negative control rather than by the zero: with `ISOMESH_TRANSITIONS=0` the same worlds report **71 open edges on the low side and 102–111 on the high side**, so the counter demonstrably reads failure |
| M-121 | **A level change moves the surface by up to 3.14 cells, which is the pop nobody measures.** Meshing a block at both its old and its new level at the instant it switches, and taking the worst vertex-to-nearest-vertex distance, over a full flight of an LOD ladder: **worst 3.136 cells**, typically 0.6–1.6. A pop is not a defect — a coarser mesh genuinely *is* a different surface — but its size is what decides whether it can be hidden by a fade, a morph, or nothing at all, and no figure for it exists in the literature review | E-205. **Method note:** it is only measurable at the instant of the switch, because after it the old mesh is gone; the demo therefore extracts the block twice on that frame rather than storing every block at every level |
| M-122 | **Re-extracting a whole LOD ladder on every level change costs 12–23 ms and hitches; re-extracting only the blocks that changed costs 4.6–12.4 ms and does not.** Twelve blocks, one or two of them changing level at a time. Caching each block's *un-inset* extraction and rebuilding only the changed ones takes the same world from **12–23 ms to 4.6–12.4 ms**, inside a 16.7 ms frame. Validation is a further **5–7 ms** and is the demo's own instrument, not a cost a game pays on a level change — timed separately for exactly that reason | E-205. **What must be cached is the extraction *before* `inset_boundary`**: the taper mutates positions in place and which faces need it depends on the neighbours' levels, which move independently of the block's own |
| M-118 | **B-003's async tests were a wall-clock race, not a plugin defect — and the number of *iterations* they ran was irrelevant to whether they passed.** `spawning_the_work_does_not_do_the_work` spun 500 `app.update()` calls with no yield, completing them in **~20 ms** against extractions needing **~49 ms**, and drained 2 chunks of 12; it failed **8 runs of 8** in isolation and the suite failed 1–3 tests per run. The discriminating probe moved wall time while *reducing* iterations — one millisecond of sleep per pass — and it drained all 12 in **39 iterations and 42 ms**. Twelve times fewer iterations, twice the wall clock, and it passed. `IsomeshPlugin` was correct throughout: the pool ran (2 chunks did complete), and `apply_finished_meshes` polls every task each frame rather than one | B-005. **Five tests carried the same defect, not four** — `the_subgrid_extractor_is_reachable_through_the_component` failed rarely enough to be missed in the first survey and surfaced only under 20 consecutive runs. All now wait on a **deadline** rather than an iteration count. **Verified 0 failures in 20 consecutive full-suite runs and 20 consecutive isolated runs**, which is the acceptance and is also the only evidence that would have distinguished a fix from a coincidence |
| M-119 | **The obvious repair would have gutted the assertion it was repairing, and the mutation test is what proves the replacement is real.** The acceptance property compared `spawn_cost * 4 < total`, where `total` was the drain loop's wall time — so *any* fix that let the loop wait would pad `total` with sleep and make a main-thread implementation pass. In the probe, **39 ms of the 42 ms total was sleep**. The denominator is now the measured cost of performing **one** extraction synchronously in the test, which removes the scheduler from the assertion entirely: the claim becomes "queuing twelve of these cost far less than doing one." **Mutation-tested by moving the work onto the main thread** — `let built = extract_chunk(..); pool.spawn(async move { built })` — which fails with *"the frame that queued twelve extractions took 149.18 ms, against 16.43 ms to perform ONE of them here"* | B-005. The reference extraction also carries its own negative control, asserting it produced triangles: a field that was free to sample would time no work and let the comparison pass for any implementation. **This is the sixth instance this session of a measurement whose denominator was the wrong quantity** (M-107, M-112, M-113, M-115, M-117) |
| M-116 | **Runtime convex decomposition costs 241–272 ms per fragment, which is 14–22 whole frames — so correct destruction colliders are viable, but never synchronously.** Measured on a Ryzen 9 5900X, `avian3d` 0.7's `convex_decomposition_from_mesh` over 23–24 fragments per target, each fragment meshed at 21³ from the intersection of the solid with the charge: **wall** 23 fragments, 211 convex parts, mean **240.7 ms**, worst **323.7 ms**; **hollow shell** 24 fragments, 305 parts, mean **271.8 ms**, worst **369.0 ms**; **spiral** 23 fragments, 224 parts, mean **249.0 ms**, worst **362.6 ms**. A 60 fps frame is 16.7 ms, so **one fragment's decomposition is 14–22× the entire frame budget**, and it falls on the frame the shot lands — the worst frame available | E-204. **The shape ranking is the predicted one, which is the useful part:** the hollow shell — whose convex *hull* is a solid ball, so a hull collider fills the cavity that makes it a shell — needs the most pieces, **12.7 per fragment** against the wall's 9.2 and the spiral's 9.7, and costs the most time. **Zero fragments failed to get a collider** on any of the three. **Consequence:** this belongs on a worker with the finished collider swapped in later, exactly as G-006's budget already treats meshing; a synchronous call is a 300 ms hitch per fragment. The README's "convex decomposition — not yet" is now a measured position rather than an unexamined one |
| M-117 | **A metric that reports a defect can be measuring the absence of a floor. E-204's first run accused 15 of 23 working colliders of tunnelling.** The demo counts fragments that fall below a plane as having gone *through* something — and the scene had no ground, so every fragment eventually fell past it and the HUD read *"15 fragment(s) left the world — they fell through something"*. Adding a static floor takes it to **1 of 23** on the wall and **0 of 24** and **0 of 23** on the shell and spiral. The one survivor is a real tunnelling case and is reported as such | E-204. **Fifth instance this session of the same shape** (M-107, M-112's controls, M-113's grader, M-115's mismatched windows): a number that looked like a finding was a property of the fixture. The rule that catches it is the one already in Part 5 — require every non-zero to prove it *could* have been zero, which here means asking what the fragment was supposed to land on |
| M-115 | **A moving body is stopped harder and more often by ordinary terrain than by a chunk join — the same shape as M-106's answer, from the test M-106 could not run.** A dynamic capsule driven at 7 m/s along E-203's path across the same streamed `fbm_terrain`, measured over 66 s and 441 m of travel: **97.0% of commanded distance actually covered**, and stalls — frames where the body advanced less than half what it was asked to — split **13 at a seam against 70 inside a chunk**, with worst single-frame shortfall **0.809 at a seam against 1.000 inside**. So the interior of a chunk both stalls a body five times as often and stops it *completely* at least once, while no seam ever did | E-206b. **This is E-203's comparison re-run with a body instead of a ray**, and it lands the same way: M-106 measured a worst vertical step of 0.412 cells at a seam against 0.539 within a chunk, and the honest reading of both is the *ratio*, not either number. A fixed stall threshold would have been measuring the landscape. **What is new is that it could fail differently.** A ray hits or misses; it is never *caught*. A seam lip that reads as a 0.2-unit step is, to a capsule, either nothing or a wall depending on the capsule, its speed and its approach angle — and nothing in E-203 could have distinguished those. **Two instrument bugs found on the way, both mine**: `commanded` accumulated from frame zero while `travelled` accumulated only after settling, reporting 12.8% progress on a body walking perfectly well; and the body spawned at y = 24 so most of any short run was a measurement of gravity. Both are the same error — comparing two quantities gathered over different windows |
| M-114 | **`HermiteCell::from_corners` is public and its contract is not: the corner order it requires lives in a private module, so no consumer can satisfy it without guessing.** The signature takes `corner_values: &[R; 8]` *"in this crate's corner order"*; `cube` is declared `mod cube;` in `lib.rs`, so `corner_offset`, `EDGE_CORNERS` and `CORNER_COUNT` are all unreachable from outside. The layout is `[c & 1, (c >> 1) & 1, (c >> 2) & 1]` and is pinned by `corner_offset_matches_the_bit_layout`, but that test is internal too. **A public function whose precondition can only be met by reading private source is an API defect, not a documentation gap** | E-114, which is the first consumer to call `from_corners` from outside the crate and therefore the first to hit it. **Worked around rather than papered over:** the example duplicates the layout and then *verifies* it at startup against the crate itself, by building a cell from a plane whose crossing position is known in closed form. **The check is mutation-tested, and the obvious version of it does not work**: swapping x and y in the duplicated `corner_offset` — the likeliest transcription slip — still produces **four** crossings, so a check that counted them would pass; it is caught only by the position, worst x error `7.5e-1`. Fixing this properly means exporting the corner order, which is an API decision and is not E-114's to make |
| M-113 | **M-112's two laws reproduce in the committed example, and its fitted constant does not survive — the accuracy figure depends on how the offset sits on the `f32` lattice, not on `ulp/h` alone.** M-112 offered `worst \|f\|/h ≈ 1.4 · ulp/h` from a standalone probe at decimal offsets. `precision_f32_vs_f64` sweeps **powers of two** instead, and at `ulp/h = 4` it measures **4.0000** where the probe's `6e6` — same `ulp/h`, same `ulp = 0.5`, different alignment — measured **5.8564**. The two agree exactly wherever the offsets agree (`ulp/h = 1 → 1.3808`, `2 → 3.1010`), so nothing is contradicted except the formula. **The law that is exact is a different one: at a fixed offset, halving `h` exactly doubles the error measured in cells** — `3.1010 → 6.2020`, `4.0000 → 8.0000`, `5.8564 → 11.7128` from 33³ to 65³. The absolute error is set by `ulp(offset)` alone, and expressing it in cells is the entire reason a finer grid looks worse | E-112. **Both of M-112's structural claims hold in the shipped code**: 33³ and 65³ are clean at `2²²` and torn at `2²³` (χ 2 → 1, 42 and 78 boundary edges), so the topology threshold is absolute; and the discriminating fixture reproduces — 65³ at `ulp/h = 8` is clean while 33³ at `ulp/h = 8` is torn. **Method:** the fitted constant came from four points that all happened to be decimal, and Part 5 already says the tidy formula is usually missing a constant. Here it was missing a *variable* |
| M-112 | **E-112's premise is wrong by an order of magnitude, and the failure it names is really two failures with two different laws.** The ticket says *"same field at ~1e6 offsets; f32 cracks, f64 doesn't"*. Measured on a unit sphere translated to `offset` on every axis, dual contouring, 33³ and 65³ over a ±2 box: **at 1e6 `f32` does not crack.** Topology is perfect — same 1,160 vertices and 2,316 triangles as `f64`, `χ = 2`, 0 non-manifold and 0 boundary edges — and only *accuracy* moves, worst radial error **0.0362 → 0.6603 cells**. What actually breaks: **(1) accuracy degrades in proportion to `ulp(offset)/h`**, measured `0.66, 1.38, 3.10, 5.86` cells at `ulp/h = 0.5, 1, 2, 4`, i.e. ≈1.4·`ulp/h` — a *relative* law that depends on cell size; **(2) topology breaks at an absolute threshold, `offset ≥ 2²³ ≈ 8.39e6`, where `ulp(f32) ≥ 1.0`** — `χ` drops 2 → 1, vertices collapse 1,160 → 271, and **54 boundary edges** appear at 33³ (102 at 65³). `f64` is untouched at every offset tried, 0.0362 cells at 1e7, identical to its reading at the origin | E-112 research, before implementing. **The two laws are independent, and that is the load-bearing part**: 65³ at `ulp/h = 8` (offset 6e6) is topologically **clean**, while 33³ at the same `ulp/h = 8` (offset 1e7) **cracks** — so the crack cannot be predicted from `ulp/h`, and the accuracy loss cannot be predicted from the absolute offset. **Two controls, both run before believing any of it.** *Is the crack in the mesh or in the validator?* `validate.rs` quantises with `as_f32() as iN` and T-008 anchors that lattice to the mesh minimum, so the same `f32` vertices were re-validated after moving them to the origin **in f64** — reports are **bit-identical**, `χ = 1` and 54 boundary edges either way, so the defect is in the mesh and T-008's anchoring is not the limit. *Is it the fixture's gradient?* `Sdf::gradient`'s default is a central difference with a step scaled by `\|p\|`, and its own doc warns that a field whose characteristic length is far from 1 should override it — so the probe was re-run with an **analytic** sphere normal. **Bit-identical again**, at every offset, which rules the gradient out entirely |
| M-111 | **`bevy_isomesh`'s four async drain tests do not pass on Linux, and one of them never passes — B-003's acceptance property is currently unverified rather than verified.** Measured 2026-08-13 at `4369e3c` on a 12-core Ryzen 5900X under CachyOS, with **no local changes to `bevy_isomesh/src`**, so this is HEAD's behaviour and not a working-tree artefact. `plugin::tests::spawning_the_work_does_not_do_the_work` fails **8 of 8** runs in isolation, meshing **2 of 12** chunks across 500 `app.update()` calls; the whole test returns in ~20 ms, so the updates are not blocking and the spawned extractions simply never complete. Three consecutive full-suite runs gave **1, 3 and 2 failures** out of the same four drain tests, so the other three are intermittent rather than sound | Found while running the E-110 definition-of-done checks, which is the only reason it was run at all. **Not caused by E-110** — that ticket adds one file under `examples/`, which the lib tests do not compile. **What this contradicts:** the previous session's handoff reported "14 plugin" tests green, and the count is right while the result is not; a suite that passes on one platform and fails on another is exactly the rule CI already earned once (`libwayland-dev`). **The obvious suspect is unchecked:** `TaskPoolPlugin::default()` sizes its pools from available parallelism and `cargo test` builds one `App` per test thread, so oversubscription would explain both the total failure in isolation and the intermittency in the suite — but it is a hypothesis, and B-005 owns settling it before anything is changed |
| M-108 | **The prediction was registered before the sweep and came out half right, and the half that failed is the more interesting one: the clamped residue falls with resolution on `gyroid` and does not move at all on `fbm_terrain`.** Clamped dual contouring, pairs per 1,000 triangles, at 17³→49³ in steps of 4. `gyroid`: **7.14, 3.02, 7.14, 2.70, 3.12, 5.26, 0.72, 0.76, 1.12** — a 6.4× fall end to end, non-monotone in the middle. `fbm_terrain`: **25.40, 20.20, 9.82, 24.01, 13.84, 19.03, 21.87, 17.46, 20.43** — 25.40 at the coarsest and 20.43 at the finest, wandering without direction in between. The registered prediction was that *both* would fall, reasoning from M-29 (the residue is multi-sheet cells) plus M-15 (multi-sheet is a resolution effect, not a topological one) | E-110, `examples/qef_clamp.rs` module header, which states the prediction and the falsification condition before the numbers. **M-15's mechanism survives; the prediction was too broad.** M-15 says *"any feature thinner than one cell forces two sheets through it"*, and that silently assumes a **fixed** stock of features to resolve. `gyroid` has one, so finer cells retire them. `fbm_terrain` is fractional Brownian motion and has structure at every scale by construction, so refining the grid uncovers new sub-cell features exactly as fast as it resolves the old ones — the residue is scale-invariant because the *field* is. **Consequence:** "raise the resolution" is not a mitigation for self-intersection on procedural terrain, which is the one field in this set a game would actually ship |
| M-109 | **M-61's "splitting the vertex makes self-intersection worse" is a 33³ fact, not a property — across nine resolutions it reverses three times and ties four.** Clamped, pairs per 1,000, dual contouring against manifold dual contouring, 17³→49³ in steps of 4. On `gyroid` MDC is worse at eight of nine grids and **better at 29³, 2.311 against 2.696**. On `fbm_terrain` the two are **bit-identical at 25³, 37³, 41³ and 45³**, MDC is better at 29³ (23.305 vs 24.011) and 49³ (19.545 vs 20.428), and worse only at 17³, 21³ and 33³. So M-61's 1.82× and 1.12× are single-grid samples of a quantity whose sign is not stable, and on the field where it matters most MDC is indistinguishable from plain DC on four grids out of nine | E-110's resolution sweep. **What does not change:** M-61's direction *at the default grid* reproduces here exactly — 3.118 → 5.669 and 13.837 → 15.434 — and its mechanism, that the clamp's partition argument assumes one vertex per cell, is untouched. What changes is the strength of claim those two numbers can carry. **Method:** the third time a figure measured at a single resolution has turned out to be resolution-dependent (V-6, M-19). A single-grid number belongs in the ledger *as* a single-grid number |
| M-110 | **The self-intersection counter never tests 99% of the pairs on a dual contouring mesh, and the skipped count is identical clamped and unclamped — which corroborates M-29 from the counter's own bookkeeping.** `gyroid` at 33³: **756 pairs found against 71,748 skipped for sharing a vertex** with the clamp off, and **33 found against the same 71,748 skipped** with it on. M-83 established the blind spot on subgrid MT's Steiner fans; on dual contouring it is larger still, because one vertex per cell means every quad shares vertices with its neighbours across every cell face. **The invariant skip count is the corroboration**: `adjacent_pairs_skipped` is a function of connectivity alone, so a clamp that changed *which triangles are adjacent* would move it. It does not move by one, which is M-29's "the clamp fixes placement, not connectivity" visible as an equality rather than an inference | E-110. Also the first count of **distinct offending triangles** rather than pairs — 483 of 10,584 with the clamp off and **54** with it on, so the red is 4.6% of the mesh and then 0.5% of it. Nobody had counted triangles; M-28 counted pairs, and a pair-count cannot say how much of a mesh is affected |
| M-107 | **λ is the sharpness/stability trade in one number, and swept over six decades it moves the runaway by a factor of 23 — but only on the fields M-30 already named.** `gyroid` at 25³, clamp off, as λ goes `1e-6 → 7e-5 → 4e-3 → 5e-1`: worst distance a vertex has to be dragged back into its own cell reads **18.03 → 8.99 → 3.41 → 0.78 cells**, and worst `|f|/h` reads **10.54 → 3.38 → 2.24 → 0.60**. On **`box_exact` the same sweep gives 0.000 for the runaway at every λ**, which is not a defect in the demo — M-30 measured the runaway at 3.18 cells on `gyroid` and 2.17 on `fbm_terrain` and said explicitly that *"sphere, box_exact and thin_plate have zero vertices outside"* | E-109, which made λ configurable (`solve_with`, `Qef::lambda`, `DualContouring::set_lambda`) after finding the ticket named a "normal-deviation threshold" this implementation does not have. **Two lessons, both about choosing the instrument before running it.** `|f|/h` cannot see the runaway *where it matters most*: a flat cell is rank 1 and its unconstrained directions lie **within the surface**, so an unheld vertex slides along the plane and stays exactly on it — `box_exact` at λ = 1e-6 reads `0.000` off-surface while being several cells from where it belongs. And the *field* has to be chosen as deliberately as the metric: opening this demo on `box_exact` and sweeping λ shows exactly half the story, which is the mistake the field list is now ordered to prevent |
| M-106 | **The acid test passes with a margin, and the margin is the interesting part: across 495 seam crossings the worst vertical discontinuity *at a seam* is 0.412 cells, against 0.539 cells *within* a single chunk.** So the seams are measurably **smoother than the terrain itself** — which is the right comparison, because a fixed height threshold would have been measuring the landscape rather than the joins. Zero probes hit nothing. Measured by casting rays against the meshed triangles through `parry3d`, not against the field: whether two independently meshed chunks *meet* is decided by G-001's overlap, and only the triangles know it | E-203. **What this does not test, and the ticket asked for:** a character *controller*. The sweep finds holes and lips; it does not find an invisible wall that only a moving capsule reveals, because nothing here slides along a surface. Named as E-206b rather than left implied |
| M-104 | **A radius-based residency rule loads a *ball* of chunks, and for a heightfield that is 4x too many with most of them meshing to nothing — 952 resident and 606 permanently waiting, against 234 and 0 once the vertical extent is bounded.** G-007's `ChunkStream` is radius-based because that is what a general residency rule is, and it is right. What is wrong is using it unfiltered for terrain: `fbm_terrain` has bounded amplitude, so every chunk more than one layer off the surface is empty air or solid rock, and each costs a **full extraction** to discover that. Filtering to two vertical layers in the *example* — which is what a real game does — takes the resident set from 952 to 234, the permanent backlog from 606 to 0, and the visible holes in the terrain to none | E-201. **The queue that never drains is the symptom to recognise:** `waiting` pinned at hundreds while `applied` trickles is not a budget being respected, it is a working set larger than the budget can ever serve, and it renders as holes that never fill. **A second effect, not root-caused:** with the oversized set, `bevy_render`'s slab allocator logged **704** use-after-free errors per run, and with the bounded set it logs **0** — consistent with meshes being freed while still queued for upload, but the race was not isolated, and the honest statement is that the churn was removed rather than the bug found |
| M-103 | **`rustdoc` had never run on `bevy_isomesh`, and two doc links were already broken — the third time an excluded workspace turned out to be excluded from a CI step nobody noticed was missing.** The root lint job runs fmt, clippy *and* rustdoc; the bevy job ran fmt, check, clippy and test. So `cargo doc -D warnings` had never been executed against that crate, and it fails immediately: a module doc linking `[`Mesh`]` with no import path, and a public function linking `[`triplanar_uv`]`, a private item. Both fixed, and the step added to the job | B-003. **E-111 is the same finding and said so at the time** — *"a workspace that is excluded from the root is excluded from the root's CI commands too. Check each one separately"* — and it was recorded after `fmt` turned out to be missing from exactly this job. The rule was right and incomplete: `fmt` was added, the *list* was not audited, and rustdoc sat missing for another twenty tickets. **The sharper rule is to enumerate the steps rather than patch the one that bit you**, because the gap is always the step nobody thinks of as linting — A-002 recorded that phrasing and this is its third instance |
| M-102 | **`ChunkId` orders lexicographically on `[x, y, z]`, so a residency sweep must iterate `x` *slowest* — and the natural `z`-outer loop produces exactly the wrong order.** G-007's diff between one frame's resident set and the next is a sorted merge, which on unsorted input is silently wrong rather than loudly: it emits spurious loads and unloads and the set drifts. Caught on the first run by a `debug_assert!` that the swept candidates come out ascending, placed there precisely because a merge cannot detect its own precondition | G-007. **The general shape is worth more than the instance.** The bug is invisible in every unit that produces a *set* — membership, counts, and "is chunk X resident" all pass — and only appears in the *diff*, which is the thing consumers act on. Asserting the precondition at the point where it is relied on cost one line and one run; discovering it downstream would have looked like chunks flickering in a demo, three abstractions away from the loop nest that caused it |
| M-101 | **§3.2.3's inset cannot be reconstructed from its prose, and two attempts made the measured result worse rather than better — the missing half is *where* it applies, not just how to triangulate it.** The point-insertion rules are fully specified (*"insert the midpoints of all polygon edges contained in any edge of `f`, and move them a small distance in the inward normal direction"*), so the first implementation applied them to every boundary-disk region: `box_exact` went from a clean `(0, 0, 0)` to **(84 non-manifold edges, 84 non-manifold vertices, 180 mis-oriented)**. Restricting it to loops confined to a single face — the per-tet expression of *"two adjacent tetrahedra sharing a face `f`"*, which needs no communication with the neighbour — still gives `box_exact` **(42, 84, 180)**. Both attempts reverted | A-014d, attempted and abandoned rather than shipped. **The mechanism the failures point at:** a polygon edge lying along an edge of `f` is shared with whatever else meets at that tet edge, and inserting a midpoint splits it — so the disk detaches from every neighbour that did *not* also insert one. Which polygons are duplicated, and therefore which edges may be split, is a property of a **pair** of tetrahedra; this implementation is strictly per-tet, and the guess that "confined to one face" captures the pair condition is what both numbers falsify. **Rule 5 is the right call here and the numbers say so**: a construction whose specification is a figure was guessed at twice, and both guesses produced meshes that would render fine and are measurably worse. The prose half is not enough |
| M-100 | **A demo can be broken in a way that looks entirely correct: E-108's letters, centred on `z = 0`, are never lost by Marching Cubes at any thickness — 576 triangles at 0.45, 0.30 *and* 0.15 voxels, the same number every time.** `z = 0` is a grid plane at every odd resolution, so a sheet centred there always contains a whole plane of nodes; every vertical edge through one registers a sign change, and the extractor the demo exists to embarrass finds the letters perfectly. The constant triangle count across a 3× thickness change is the tell, and it is the kind of tell that is only visible if you look at the number rather than the picture — the render looked exactly as intended. Fixed by offsetting the mid-plane to **0.41 of a cell** (off the lattice, and not a round fraction that a resolution change lands back on) and tilting it by `0.03` in `y`. Measured after: Marching Cubes gives **220 → 0 → 0 → 0** triangles at 0.70, 0.45, 0.25, 0.10 voxels while subgrid gives **1706 → 1380 → 1304 → 1252** | E-108. **Fifth occurrence of the fixture trap** (M-32, M-38, M-44, M-94), and the first in an *example* rather than a test — which is worse, because a demo has no assertion to fail. The tilt earns its place separately: without it the failure is binary, every node inside or none, and M-72's middle phase — the holey remnant that is what a streamed world actually suffers — never appears at all |
| M-99 | **The subgrid mesh's connectivity is provably manifold and my weld is what breaks it — the same mechanism as M-59, in a second algorithm. Unwelded: 0 non-manifold edges and 0 non-manifold vertices on every field. Welded: 3 and 6 on `csg_difference`, 4 and 6 on `fbm_terrain`.** Appendix A is explicit — Theorem A.1: *"In the even sum case, both the primal and dual algorithms yield manifold connectivity"*; Lemma A.2: *"The primal polygons constructed **before** Section 3.2.3 form an edge-manifold cell complex."* The precondition was checked rather than assumed: over all seven reference fields at 17³, **0 of 98,304 tetrahedron faces violate the even-sum condition** — so the theorem applies everywhere and the guarantee is in force. What §3.2.3 actually says is that the pair of coincident polygons *"define manifold connectivity, but degenerate geometry"*, and a weld by position cannot tell two coincident-but-distinct polygons apart | A-014d, measured after retrieving Appendix A through home-still, which no amount of §3.2's own text would have supplied. **This re-aims A-014d.** Its job is not "make the output manifold" — it already is, as a complex. It is "separate coincident geometry far enough that an indexed, welded triangle mesh can still represent it", which is exactly the inset. And **M-59 is the precedent, not a coincidence**: that entry recorded the same shape for Manifold Dual Contouring — a manifold complex an index buffer cannot express. Two algorithms, one representational limit |
| M-98 | **Subgrid Marching Tetrahedra costs 70× classic Marching Tetrahedra and 196× Marching Cubes, and the ratio is field evaluations rather than anything algorithmic.** On `sphere`, f32, Apple M5: at 33³ Marching Cubes **115 µs**, Marching Tetrahedra **323 µs**, subgrid **22.5 ms**; at 65³ **965 µs**, **2.68 ms**, **173 ms**. Other fields at 33³: `torus` 32.5 ms, `box_exact` 36.8 ms. Scaling 33³→65³ is **7.7×** against 8× the cells, so it is linear in cells like the others — the constant is the whole story. **The constant is not mysterious.** At 16 samples per edge, a cell costs `6 tets × 6 edges × 16` = 576 field evaluations before any bisection refinement, against Marching Cubes' 8 shared corner samples; ~72× is the ratio that predicts the measured 70× against Marching Tetrahedra | A-014c, `benches/extract.rs`, `cargo bench --bench extract`. **Two consequences.** This is a *correctness-first* implementation and says so: every cell re-finds the roots on edges its neighbours already found, deliberately, because identical endpoints through a deterministic root finder is what makes conformity hold without a cache — a grid-edge cache is the obvious optimisation and the redundancy is large, but it has a correctness precondition and is unmeasured. And **the trade is real rather than a defect**: 22.5 ms buys geometry no sign-based method produces at any grid resolution (M-95), so the comparison to make is not "subgrid vs Marching Cubes at 33³" but "subgrid at 33³ vs Marching Cubes at whatever resolution resolves the same feature", which for `thin_plate` is *no resolution at all* |
| M-97 | **§3.2.3's immersion is not hypothetical, and it is the single cause of every violation the subgrid extractor produces: 4 of the 7 reference fields are clean, 3 are not, and both failure modes trace to the same missing inset.** At 17³ with 16 samples per edge, welded, `validate_indexed` reports **(non-manifold edges, non-manifold vertices, inconsistently oriented)** of `sphere (0,0,0)`, `torus (0,0,0)`, `box_exact (0,0,0)`, `thin_plate (0,0,0)`, **`csg_difference (3,6,6)`**, **`gyroid (0,0,138)`**, **`fbm_terrain (4,6,19)`**. The two shapes of failure look unrelated and are not. `csg_difference` and `fbm_terrain` show the *non-manifold* form §3.2.3 names outright — two adjacent tetrahedra each emitting an oppositely-oriented copy of a polygon bound by the same non-normal loop, whose union pinches. `gyroid` shows the *orientation* form, and the mechanism was measured rather than guessed: of its 8,088 triangles, **12 have zero area and 24 have the gradient lying exactly in the triangle's plane** — worst `|cos| = 0.00e0` — so the orient-against-the-gradient vote is not merely close, it is **undecided**, on precisely the degenerate boundary-disk triangles V-21 says are supposed to be there | A-014c/A-014e, `the_validity_suite_over_every_reference_field`. **Both are downstream of one fix**: §3.2.3 pushes those polygons into the tet interior, which gives them area and a normal transverse to the gradient — so the union stops pinching *and* the orientation vote becomes decisive. Pinned as exact counts rather than asserted to zero, per the Phase 1 amendment: a known defect with a number and a ticket that owns it satisfies the gate, and **A-014d now owns it with two measured fixtures instead of a hypothesis**. **Sharpened by M-99:** the counts are not a manifoldness failure at all — connectivity is provably manifold under Theorem A.1, whose even-sum precondition holds on all seven fields, and welding is what collapses the coincident pair |
| M-96 | **Orienting each triangle against the gradient at its own centroid is sufficient — the welded output is closed, manifold and consistently oriented with `χ = 2` on both `thin_plate` and `sphere`. But the weld is a precondition, not a tidy-up, and the numbers before it are worthless.** Unwelded, `thin_plate` at 17³ reports **2,240 boundary edges** across 896 triangles, because the extractor emits each tetrahedron's vertices independently and no two tetrahedra share one — it is a triangle soup with no topology to check, and the orientation figure of 8 was measured only over edges *interior to a single tetrahedron*. Welded: **840 triangles, 0 inconsistently oriented, 0 non-manifold edges, 0 non-manifold vertices, 0 boundary edges, χ = 2**, and the 56 degenerate triangles collapse away in the weld | A-014e, `the_welded_output_is_a_closed_consistently_oriented_manifold`. **Per triangle, not per patch:** `thin_plate`'s two faces are 0.4 cells apart and routinely land in the same tetrahedron facing opposite ways, so a per-patch decision would get one of them backwards. **Second time this session that a measurement on unwelded subgrid output was meaningless** — M-93 was the self-intersection count — which is why the extractor's doc comment now says it outright rather than leaving it to be rediscovered a third time |
| M-95 | **A-014c's acceptance is met with a number rather than a hedge: `thin_plate` returns 4,328 triangles at 33³ where greedy quads returns 0 on the same grid, and every vertex is on the surface to within `1e-9`.** At 17³ it is 896 triangles from 2,248 vertices. The contrast is the whole point of the track — A-005 measured greedy quads' zero because no cell *centre* is inside a plate 0.4 cells thick, and M-72 measured Marching Cubes aliasing the same field into a resolution-dependent scatter. Sharper still, with the **grid held fixed**: a slab 1/20 of a cell thick gives 0 triangles at 2 samples per edge and a mesh at 256, which no sign-based method can do at any sample count | A-014c, `thin_plate_comes_back_where_greedy_quads_returns_nothing` and `the_extractor_resolves_a_feature_the_grid_cannot`. **And a determinism boundary worth having in writing.** Raising `samples` from 16 to 32 leaves the topology *identical* — same vertex count, same index buffer — but moves the positions, by under `1e-12` and by more than zero. Bisection converges to *an* ulp of a root and which one depends on the bracket it started from. So the guarantee is **"same arguments, same output", not "same field, same output"**, and a golden hash over this extractor must pin `samples` alongside the grid. The first version of that test asserted bit-equality across sample counts and failed |
| M-94 | **1D root finding resolves a slab at 1/1000 of the edge length — and the fixture that appears to demonstrate its limit passes for the wrong reason unless the slab is moved off the sample lattice.** `all_roots` finds exactly two crossings on a slab of half-width `w` centred mid-edge for `w = 0.25, 0.05, 0.01, 0.001`, each to within `1e-12`, while **both endpoints read outside** — which is A-005's zero-triangle result at its origin: a sign test asks the endpoints and learns nothing. The stated limit is real too: 8 samples step over a slab `2e-4` wide, and 100,000 find it, which is §1.3's *"1D marching can of course miss intersections, \[but\] we are no worse off than classic marching"* | A-014c, `a_slab_thinner_than_the_grid_still_has_two_roots` and its companion. **The trap:** the miss-it test was first written with the slab centred at `0.5`, where 8 samples put one sample *exactly* on the centre — it reads inside, the slab is found, and the assertion that it should be missed fails. Centring at `0.5137` fixes it. **Fourth occurrence of the fixture trap** (M-32, M-38, M-44) and the first where the fixture was written for this ticket rather than inherited |
| M-93 | **Subdivision's output reported 30 self-intersections in 52 triangles, and every one was an artefact of vertex duplication rather than a fold. Welding first takes it to zero.** Each sub-tet appends its own copy of the crossings on a face it shares with a sibling, so two triangles that legitimately *meet* there carry different vertex indices — and `self_intersections` only skips pairs sharing an index (M-83), so it compares them and finds them touching. Diagnosed by elimination rather than by inspection: all four children of `(4, 2)` fill to **0** self-intersections standalone (8 + 8 + 18 + 18 = 52, matching the total), so nothing was wrong *inside* a child; and all four children are **normal** curves, which ruled out the first hypothesis — §3.2.3's immersed Δ-complex, which only arises from non-normal loops on a shared face | A-014b, `every_implemented_case_emits_an_intersection_free_patch`. **Consequence for every future measurement on this extractor:** self-intersection counts are meaningless on an unwelded multi-tet patch, and the weld must be by *exact* position, which is sound only because both sub-tets compute a shared crossing from the same corners and the same parameter and get bit-identical results. Same reasoning as M-32's, and the same reason a tolerance would hide the failure rather than absorb it |
| M-92 | **§3.2 is complete over the tested space: 4,096 of 4,096 configurations mesh, with `NoPattern` and `Inconsistent` both zero — and 105 of them correctly triangulate to *nothing*.** Corner type was the last of §3.2.2's three, needing the disk plus two extras: the vertex coinciding with the inside corner omitted from its polygon, and a triangle capping that corner, built from the crossings nearest it along its three incident edges. With it in, every configuration with 0–3 crossings per edge fills. **The 105 is the interesting number.** Those carry a non-normal loop whose inside region is a **bigon** — one chord and one edge piece, two nodes — and a fan over two nodes is empty. The minimal case is `e = (2, 0, 0, 0, 0, 0)`, a single scoop pair. That is V-21 at its most extreme and it is right: the region has zero area until A-014d insets it into the tet interior | A-014b, `all_three_non_normal_types_now_have_a_spanning_disk` and `the_coverage_of_the_implemented_cases_is_pinned`. **The assertion this replaced was "every non-normal loop produces triangles", which failed** — and rewriting it as a pinned count rather than deleting it is the point: a bigon disk is a real state of the algorithm, and a test that merely tolerated it would stop distinguishing "degenerate by construction" from "emitted nothing by mistake". `SingleLoop` and `Subdivision` stay at 0 here because they remain unreachable at these counts, not because they are covered |
| M-91 | **The contractible spanning disk introduces no vertices at all, and coverage moves 3,394 → 3,808 of 4,096.** §3.2.2's disk for corner and contractible loops is built *in the tet boundary*, so every one of its vertices is a crossing the tet already has — no Steiner point, no centroid, nothing new. That makes it **cheaper than any of the four Steiner cases** despite being the one that needed a planar arrangement to write, which is the opposite of what the case list's difficulty ordering suggests. With it wired in, the sweep over 0–3 crossings per edge moves from 3,394 filled / 702 refused to **3,808 filled / 288 refused**, and every one of the remaining 288 carries a **corner-type** loop — the only kind left, and the only one needing §3.2.2's two extras (omit the vertex coinciding with the inside corner, add a triangle at that corner) | A-014b, `the_coverage_of_the_implemented_cases_is_pinned`. A `Node::Corner` reaching the emit is treated as [`Unfilled::Inconsistent`] rather than fanned over: a contractible loop marks all four corners outside so no qualifying region can contain one, and for a corner-type loop fanning over it would produce a plausible wrong disk instead of the omission rule |
| M-90 | **A scoop belongs to one face, but its two endpoints lie on an edge that belongs to two — so "both endpoints are on this face" hands the same scoop to both, and a crossing ends up bounding three regions instead of two.** Building §3.2.2's `σ \ γ` needs to know which face each of `γ`'s segments lies in. For an ordinary segment that is determined: its two endpoints are on two different edges, and two distinct edges of a tet lie in a common face only if that face is unique. **A scoop breaks it** — both endpoints on one edge, and every edge is shared by exactly two faces, so the endpoint test cannot tell which one owns it. Measured on `e = (2, 2, 0, 0, 0, 0)`, face 3: crossing `(edge 1, index 1)` appeared in **3** regions. §3.1 already has the answer, because `face_segments` builds segments per face in the first place, so the fix is to take the face assignment from there rather than re-derive it | A-014b, `a_face_split_along_a_loop_partitions_into_regions`. **Found only because the weaker checks passed.** Edge-piece conservation and `regions == chords/2 + 1` both held with the bug present — they count things without noticing *which* region got them. The check that caught it is order-sensitive: a crossing lies between exactly two regions and a corner belongs to exactly one. The same assertion also killed a second mutation the counting checks had survived (peeling a chord that is not innermost), so one order-sensitive property was worth more than two conservation laws |
| V-21 | **The corner/contractible spanning disk is *supposed* to come out degenerate, and mistaking that for a bug would send A-014b chasing a non-problem.** §3.2.2 builds this disk *"contained mostly in the tet boundary, rather than its interior"* — so its triangles are coplanar with the tet's faces by construction. Worse-looking still, a **scoop** (a segment joining two crossings on the *same* edge) realises as a straight chord lying *along* that edge, so the region it bounds with the edge pieces between its endpoints has **zero area**. Both are expected: §3.2.3 says the mesh at this point *"has manifold connectivity, but may be an immersed Δ-complex rather than an embedded simplicial complex"*, and the repair is exactly to *"push polygons into the tet interior: we insert the midpoints of all polygon edges contained in any edge of `f`, and move them a small distance in the inward normal direction"* | Read from `10.48550/arXiv.2606.00454` §3.2.2–3.2.3 at A-014b, before implementing, specifically to find out how the case should be judged. **Two consequences for the work.** A-014b's corner/contractible output must be graded on *connectivity*, not on area or on self-intersection — CLAUDE.md already treats degenerate triangles as a recorded metric rather than a gate, and here they are the specified output. And the arrangement should be built **combinatorially** rather than from straight-line geometry, since a geometric build collapses every scoop before A-014d gets the chance to inset it |
| M-89 | **On `thin_plate` — the field this whole track exists for — 93.75% of tets already mesh, and every one of the remaining 6.25% is the *contractible* case. Not corner, not diagonal: contractible, all 192 of them.** Sampled at 33³, decomposed into the 6-tet cube, with all roots along each edge found by dense 1D sampling: **3,072 tets carry crossings**, of which **2,880 fill completely** and **192 stop at a non-normal loop**. The curve census is **0 open / 3,060 normal / 192 non-normal**, and the non-normal breakdown is **192 contractible / 0 diagonal / 0 corner**. **620 edges carry more than one root**, which is the sub-voxel signal itself — the thing a sign test cannot represent, and the reason A-005 measured greedy quads returning *zero* triangles on this field | A-014b, measured to decide whether §3.2.2's planar-arrangement case was on the critical path for A-014c. **It is**: 192 holes is not a mesh. But the useful half is that it is one construction and the *simpler* one — contractible loops mark all four vertices outside, so they need no distinguished vertex, no corner triangle and no vertex omission, which are the three complications the corner case adds. **Provenance, resolved at A-014c:** first measured with a throwaway probe carrying its own ad-hoc root finder, deleted rather than committed because a second root finder is the two-path failure the crate forbids. `the_thin_plate_census_reproduces_with_the_real_root_finder` now re-runs it through `all_roots` and gets **all three numbers unchanged** — 3,072, 620, and 0/3,060/192 — so the probe and the shipped root finder agree on counts, and the census is a committed test rather than a figure that only existed in a commit message |
| M-88 | **§3.2.2's two labelling rules are stated separately and are checkable against each other, which turns a transcription into a verified one.** The vertex rule assigns sides to the tet's four corners (corner type: the distinguished vertex inside, the rest outside; contractible: all outside; diagonal: no distinction). The edge rule walks an edge from its lower corner `i`, takes `i`'s side, and flips at every crossing — **and never looks at the far corner `j`**. So the side it arrives at is a prediction, and it must equal the side `j` was independently assigned. The two agree only if crossing `γ` an odd number of times means opposite sides, which is exactly what the parity bit `b_ij` is defined to say — so the check couples the parity rule to the vertex rule, and nothing else does. Verified over all 4,096 configurations with 0–3 crossings per edge, every corner-type and contractible loop, every edge: **no disagreement**. Mutating the vertex rule so a corner loop marks nothing inside fails at `e = (1, 1, 1, 2, 0, 0)` | A-014b, `the_edge_labelling_closes_on_the_far_corner_it_was_not_told_about`. **The distinguished corner is derivable rather than given**: it is the vertex whose three incident edges all have odd parity, which is consistent with the classification's `p = b₀₁ + b₀₂ + b₀₃` giving 3 when it is corner 0 and 1 otherwise |
| M-87 | **§3.2 as implemented meshes 3,394 of 4,096 configurations, the entire remaining gap is one construction, and Property II held every single time.** Over all edge coordinates with 0–3 crossings on each of six edges: **3,394 fill completely (82.9%)**, **702 stop at a non-normal loop**, and `SingleLoop`, `Subdivision`, `NoPattern` and `Inconsistent` are all **0**. Three of those zeros mean different things and the distinction is the point. `Inconsistent` is a bug and never a case, so zero is a pass. `SingleLoop` and `Subdivision` are zero because they are *unreachable* at these counts — both need a residual pattern with `ℓ > 8`, which wants larger `d₁, d₂` than three crossings per edge can supply — so they are covered by targeted fixtures instead, and this sweep says nothing about them. **`NoPattern` at zero is the strong one**: Property II held on `Γ_normal`'s own residual in all 4,096 configurations, which is an empirical check of Theorem B.3 across the whole space rather than on the 27 constructed patterns of M-82 | A-014b, `the_coverage_of_the_implemented_cases_is_pinned`, pinned as exact counts so a regression moves a number rather than rounding away. **The useful consequence for planning:** the entire 702-configuration gap is §3.2.2's corner and contractible types — one construction, not four — so what is left of A-014b is narrower than the case list suggests |
| M-86 | **Conforming does not mean "no boundary at a shared face" — it means "boundary exactly where an open curve was discarded", and the difference is measurable.** Two tets sharing a face, filled independently with no communication and no second pass, merged and welded: on `e_shared = (1, 1, 2)` with differing apex coordinates, the shared face carries **one** edge used by a single triangle. That is not a crack. §3.1 discards curves with a degree-1 vertex, and whether a given shared-face segment lands in an open or a closed curve depends on the *whole* tet — so one tet can legitimately drop a segment its neighbour keeps, which is exactly what the paper means by *"such curves are discarded, but their segments may still appear in neighboring tets as part of the mesh boundary."* Every one-sided shared-face edge across the fixture is accounted for by a discarded open curve, none is used by more than two triangles, and the excuse fires at least once so it is not a clause that passes by never applying | A-014b, `the_shared_face_of_two_tets_carries_no_crack`, the ticket's acceptance criterion. **The first version of this test asserted zero one-sided edges outright and failed** — and the failure is the finding, because "conforming" had been written down as the stronger and wrong property. The right gate is the conditional one, and it is what a full grid needs: a tet interior to the grid has every neighbour present, so its open curves close through them |
| M-85 | **Property II is a statement about `Γ_normal`'s residual, not about the tetrahedron's — and taking it for the tet's silently disables §3.2.2 on every configuration that needs it.** Subtracting only the corner cuts from the tet's own edge coordinates leaves in the points belonging to *non-normal* loops, which are in neither set, so `Pattern::of` fails for a reason that has nothing to do with §3.2.1 — and because that failure is an early return, the non-normal handling below it never ran. Exhibited on `e = (0, 0, 2, 0, 0, 0)`: one non-normal loop, no normal residual at all, reported `NoPattern` where the loop should have been triangulated. The fix is to sum the residual coordinates from the residual *normal* loops directly, which cannot make the mistake, and to run §3.2.2 before §3.2.1's dispatch rather than after it — the two sections partition `Γ`, and one configuration can carry both | A-014b, found by `a_diagonal_non_normal_loop_fans_around_its_centre_of_mass`. **The paper's structure is the tell and it was there to read**: §3.2.1 opens *"we first consider the subset of normal curves `Γ_normal ⊂ Γ`"* and §3.2.2 *"we next consider the subset of non-normal loops `Γ_nonnormal ⊂ Γ`"*. Two subsets, two independent procedures. Implementing them as sequential stages of one pipeline, where the first can bail out, is what created the coupling |
| M-84 | **The Figure-13 subdivision stencil closes: all four sub-tets it produces are themselves normal configurations, on all 27 patterns tried, and the margin is zero.** The stencil is asymmetric — `e_ai = 2d₂` against `e_aj = d₁`, `e_ak = d₂`, `e_al = d₁ − d₂` — so the recursion is only well-founded if every tet it makes is a legal input to the same procedure, and the paper states that nowhere. Checked directly: over `1 ≤ d₁ ≤ 6`, `0 ≤ d₂ ≤ d₁`, every one of the 108 sub-tets satisfies both normality conditions. **What makes this a real check rather than a formality is that it is tight**: three of each sub-tet's four faces meet the triangle inequality with *equality* (`2d₂ ≤ d₂ + d₂`, `d₁ + d₂ ≤ 2d₂ + (d₁ − d₂)`, `d₁ ≤ (d₁ − d₂) + d₂`), so there is no slack for an off-by-one to hide in — and mutating `2d₂` to either `d₂` or `d₁ + d₂` fails the sweep | A-014b, `the_subdivision_stencil_produces_four_normal_sub_tets`. **Also: the labelling is load-bearing, not cosmetic.** Because the stencil is asymmetric, which corner is called `i` decides the entire subdivision, so `Subdivision::label` searches for the `i, j, k, l` Property II guarantees instead of trusting corner order to supply it, and returns the lexicographically smallest when `d₁ = d₂` makes several equivalent |
| M-83 | **T-002 is structurally blind to a fold *inside* a Steiner fan, and §3.2 is built entirely from Steiner fans — so "zero self-intersections" means less here than it does for Marching Cubes, and exactly how much less is now pinned.** The counter excludes any triangle pair sharing a vertex index (`self_intersection.rs:266`), which its own docs already price as *"a fold that pinches exactly at a shared vertex is not counted"*. A fan gives **every** triangle of a loop the same apex, so no two triangles within one loop are ever compared. Found by a reachability check that failed: collapsing all three Steiner points of a `(3, 3)` octagon pattern onto one apex — a definitively wrong mesh — reports **0 intersections**, because the collapse is what creates the sharing. **The boundary is sharp and the useful half survives:** distinct loops share no crossing vertex and no Steiner point, so *cross-loop* intersections are fully visible, and reversing the innermost-to-outermost assignment on the same pattern does report a non-zero count | A-014b, `reversing_the_steiner_assignment_is_visibly_wrong`. **Method consequence:** a mutation test on fan geometry must permute the apexes, never merge them, or it silently tests nothing. The first version of this test merged them and passed for the wrong reason — the fourth occurrence of the fixture trap (M-32, M-38, M-44), and the second caught by an assertion rather than a reviewer. **Taken to its conclusion, the single-loop case's zero is entirely vacuous**: `m = 1` is one fan, so on a `(3, 2)` pattern — 20 triangles, 190 pairs — the counter skips **190 of 190** and reports zero having compared nothing. `the_single_loop_zero_is_vacuous_and_this_is_what_makes_that_visible` asserts that equality rather than the zero, so the green tick cannot be mistaken for coverage, and it fails loudly if the counter ever learns to see inside a fan |
| M-81 | **§3.2.1's first two cases cover the whole of classic Marching Tetrahedra, and the rest of the machinery exists entirely for what a sign test cannot see.** All sixteen sign patterns fill completely with corner cuts (`ℓ = 3`) and quads (`ℓ = 4`) alone — nothing reaches the octagon, single-loop, subdivision or non-normal branches, and the triangle counts come out 0 / 1 / 2 exactly as A-003's table has them. This is M-67 from the other end: the 8 classic configurations of 181 are also the 8 easiest to triangulate, so **every line of §3.2 past the quad case is there to serve the 173 that a sign test reads as one of those 8** | A-014b, `every_classic_configuration_fills_completely`. Useful as a staging property, not just a curiosity: the two implemented cases already subsume the extractor this crate ships today, so the remaining cases can land incrementally without a half-working path existing at any point |
| M-82 | **Theorems B.4 and B.6 hold against a reconstruction that does not use them — 27 of 27 `(d₁, d₂)` patterns, exactly.** The appendix predicts the number of boundary components as `gcd(d₁, d₂)` (Theorem B.4) and every component's length as `4(d₁ + d₂) / gcd(d₁, d₂)` (Corollary B.6). §3.1's segment reconstruction computes neither — it pairs points face by face from the edge coordinates — yet over every `1 ≤ d₁ ≤ 6`, `0 ≤ d₂ ≤ d₁` the component count and every component's length match the formulas with no exceptions, up to a single loop of length 20 at `(3, 2)` and two octagons at `(2, 2)` | A-014b, `the_appendix_formulas_predict_the_curves_section_3_1_actually_finds`. **Both directions are load-bearing**, so the test is mutation-checked: replacing `gcd` with 1 fails at `(2, 0)`, and `4(d₁+d₂)/g` with `8(d₁+d₂)/g` fails at `(1, 0)`. The length formula is partly forced by arithmetic — `m · ℓ` must equal the `4(d₁ + d₂)` total crossings — but **the component count is not**, and that is the part that says the curve reconstruction agrees with the topology the paper derives |
| M-80 | **Every point §3.1 leaves unpaired is forced by arithmetic, not chosen — which is what makes the procedure implementable without inventing a convention.** Step 2 handles an odd sum by subtracting one from each edge coordinate, *"effectively creating three open endpoints"*, without saying *which* three. It does not need to: the corner cuts claim `cᵢ` from one end of an edge and `cⱼ` from the other, and on the reduced coordinates `cᵢ + cⱼ = eᵢⱼ − 1`, so exactly one point — the one at index `cᵢ` — is unclaimed. Verified over 2,875 odd-sum faces, each leaving exactly `e − 1` endpoints used per edge. The same holds in step 3: the residual run on the long edge is `[cᵢ, eᵢⱼ − cⱼ)`, of size exactly `r = eᵢⱼ − eⱼₖ − eₖᵢ` | A-014b. This mattered for rule 5: an under-specified step is normally where a construction gets guessed at, and here the specification is complete once the arithmetic is done. Measured alongside: 8,424 even and 14,256 odd residual cases exercised, and over configurations with counts to 3 the curves come out **5,568 open, 3,682 normal, 771 non-normal** — all three of §3.1's kinds reachable |

---

## Part 3 — Verified from primary sources (tier V)

| # | Finding | Source |
|---|---|---|
| V-1 | Bevy 0.19 pins **wgpu / wgpu-types / naga 29.0.3, glam 0.32.0, encase 0.12** | `bevy_render/Cargo.toml` @ v0.19.0 |
| V-2 | Bevy 0.19 **removed `RenderGraph`**; passes are systems in ECS schedules; non-camera work targets the `RenderGraph` schedule | 0.18→0.19 migration guide |
| V-3 | Marching Cubes peak: **5.42 G voxel/s, 330 M tri/s** (RTX 2080 Ti). DMC costs 1.52–3.50×; FlexiCubes 2.77–3.92× | Grosso & Zint, `10.1007/s00371-021-02139-w` |
| V-4 | **Contouring 68 ms vs halfedge construction 58 ms** — extraction is 54% of a usable mesh | same, Table 5 |
| V-5 | On unstructured grids, Delaunay/MT ratio **15.3×–81.5×** — contouring is 1–2% of the pipeline | TetWeave Table 3 |
| V-6 | **73% of FlexiCubes' 64³ Marching Cubes timing is fixed launch overhead** (fitted a ≈ 1.88 ms) | fit over FlexiCubes' own resolution series |
| V-7 | Cross-paper reproducibility floor is ~1.5× **in opposite directions**: TetWeave re-measured FlexiCubes at 128³ as 9.63/15.25 ms vs FlexiCubes' own 14.06/9.53 | both papers |
| V-8 | GPU Marching Cubes throughput has not tracked hardware: **10.7× more bandwidth bought ~1.7× more throughput** (GTS 450 → 2080 Ti) | speed analysis |
| V-9 | Same Marching Cubes, compute shader → mesh shader: **114.2 → 2679.4 fps (23.4×)** | Elliott MSc, Waikato 2022 |
| V-10 | CBT sum-reduction, atomics → LDS staging: **5.78 → 0.40 ms** | Unity SIGGRAPH 2021 (see ✗7) |
| V-11 | Meshlet compression: **15.5 M tri in 0.59 ms** (RX 7900 XTX) | `10.2312/vmv20241204` |
| V-12 | Work graphs: 79,710 instances in 3.74 ms — **but 2.8–3.4× slower** on classification workloads | `10.1145/3675376` + independent profile |
| V-13 | nvblox: meshing is the least GPU-accelerable stage, **×3–13 vs fusion's ×174–177** | nvblox |
| V-14 | Aokana renders **10¹⁰ voxels at 6 ms**, 5% resident, RTX 3060 Ti — **explicitly not editable** | Aokana |
| V-15 | CoACD vs V-HACD: **49% → 80%** downstream manipulation success | CoACD |
| V-16 | Dimforge migrated parry (0.26.0) and rapier (0.32) **off nalgebra onto glam**, citing rust-gpu support; performance *"nothing changed, at all"* | dimforge.com, 2026-01-09 |
| V-17 | **No paper since 2020 benchmarks Marching Cubes vs Surface Nets vs Dual Contouring against each other.** Surface Nets has no credible published timings at all | literature review round 1 |
| V-18 | **Dual Contouring's own paper quantifies the f32 QEF failure.** At 256³, `bᵀb` reaches ~10⁶; f32 carries six decimal digits, so `E[x]` evaluated on a flat region — where it should be zero — has error **on the order of 1**. The paper's own remedy is double precision | Ju, Losasso, Schaefer & Warren 2002, `10.1145/566570.566586`, §2.3, read this session |
| V-19 | **Dual Contouring's topology is Surface Nets' topology.** The paper's algorithm is literally: vertex at the QEF minimizer for each sign-changing cube, quad joining the four cubes of each sign-changing edge. Only vertex *placement* differs | same, §2.2 |
| V-20 | A QEF is stored as `AᵀA` (symmetric 3×3), `Aᵀb` (3-vector) and `bᵀb` (scalar) — 10 floats — rather than as `A` and `b` | same, §2.3 |
| V-34 | **Manifold Dual Contouring's uniform-grid criterion is one vertex per cycle of a *decider-modified* Marching Cubes table, and the paper claims that surface is unconditionally manifold.** §3, *Contouring on a Uniform Grid*: *"On a uniform grid, **DC leads to nonmanifold vertices and edges for all of the ambiguous sign configurations** in the original MC algorithm. … Nielson associates one vertex with each cycle of a **modified MC table [26]**. … each edge is associated with exactly one vertex. … **this surface is always a manifold** because the original MC algorithm always constructs a manifold and the dual preserves the topology of the surface."* **Reference [26] is Nielson & Hamann, *The Asymptotic Decider*, VIS 1991**, and [13] is Nielson's *Dual Marching Cubes*, VIS 2004 — so the face ambiguity is resolved **upstream, inside the table the cycles are read from**, and the dual walk needs no rule of its own. That is the question A-022 was blocked on and could not answer. **The claim itself is ✗19** | Schaefer, Ju & Warren, `10.1109/TVCG.2007.1012`, read at A-022 from `cs.wustl.edu/~taoju/research/dualsimp_tvcg.pdf` |

---

## Part 4 — Open questions

Each has the test that would settle it. **An open question with no proposed test is a wish.**

| # | Question | Settled by | Why it matters |
|---|---|---|---|
| O-1 | ~~What fraction of cells actually change per brush stroke?~~ **Settled at G-002 (M-33, M-34), and confirmed live under a mouse at E-202 (M-50).** **15–36%** of cells in the brush's bounding box move a *triangle*, against ~100% whose sample value moves — the distinction is the whole finding | G-002 instrumentation; hash cell slabs, log per stroke | **Unpublished**, and the ceiling on every incremental-repair idea in the opportunities doc. The ceiling is high enough to be worth having: incremental re-meshing does 15–36% of the work, not 100% |
| O-2 | ~~Does clamping the QEF vertex to (1−ε) inside its cell eliminate self-intersections?~~ **Settled at A-009 (M-28, M-29): not entirely, and the residue names its own mechanism.** λ → **exactly 0** on five of seven fields; `gyroid` and `fbm_terrain` drop **23×** and **13.7×** without reaching it, and those two are precisely the multi-sheet-cell fields — so what the clamp cannot fix is *connectivity*, which is A-010's problem, not placement | A-009: measure per 1,000 triangles, clamp on vs off, all seven fields | Decided whether guaranteed intersection-free extraction is free → whether runtime convex decomposition can stop failing. **Answer: free for placement, not sufficient overall** |
| O-3 | Marching Cubes vs Surface Nets vs Dual Contouring vs MT — actual relative speed on one machine? | M-001 | The comparison does not exist (V-17). We'd have the only apples-to-apples measurement |
| O-4 | ~~Do brush operations commute?~~ **Settled at G-003 (M-36, M-37, M-38): conditionally, and the condition is narrow.** One result from 40,320 orderings for a run of same-kind *hard* edits; **11** across an add/subtract boundary (semantic, unrepairable); **40,317** for smooth union. `BrushOp::commutes_with` returns that honest answer | G-003: 8 ops × 40,320 orderings, count distinct results. Expect 1 | The coordination-free multiplayer story survives inside a run and dies at every boundary. **That is a protocol's problem, not this crate's** — networked editing is closed out in `BACKLOG.md` on exactly this evidence |
| O-5 | Do mesh shaders work on macOS/Metal? | GPU-007 capability probe | **Sources contradict:** wgpu's spec table lists MSL as *planned*; the tracking issue says the Metal HAL backend merged. Neither is trustworthy until probed |
| O-6 | What is amortized meshing cost per frame under continuous editing? | E-206 under a deliberately overloaded queue | The only number a game cares about, and no paper reports it |
| O-7 | What fraction of *our* pipeline is contouring vs everything else? | M-003 | V-4 says 54% for someone else's code with no physics. Ours is probably worse |
| O-8 | Does Dual Contouring's vertex placement need f64 in practice, or is f32 enough? | E-112, with the QEF condition number in the HUD | `M = AᵀA` squares the condition number. **Half answered by V-18**: the original paper measures f32 error ~1 on flat regions at 256³, and recommends f64. **Partially answered by M-23**: on extraction paths with no solve, `f64` costs only 8–10% of wall time, so precision is cheap where there is no QEF. Still open for the vertex solve itself, and for *our* fields at *our* resolutions — which sidesteps `AᵀA` entirely and may not degrade the same way |
| O-9 | How much does T-003's gradient-flow chord **over**-estimate distance at a concave seam? | A comparison against nearest-point search over a dense surface point cloud, or E-104 once Dual Contouring lands | The chord follows `∇f` to the zero set, which near `csg_difference`'s seam can land further away than the true nearest point. The bias direction is known and safe for a "below X" gate; the *magnitude* is not measured, and M-001's shootout column would inherit it. `csg_difference` measured forward `0.0833` at 33³ — how much of that is seam bias is unknown |
| O-10 | ~~What is Surface Nets' non-manifold **rate** as a function of feature thickness over `h`?~~ **First curve measured at A-010 (M-60)**, as the multi-sheet-cell rate: `gyroid` 3.13% → 2.05% → 0.53% and `fbm_terrain` 1.70% → 0.84% → 0.77% at 17³/25³/33³, and exactly zero on the other five fields at every resolution. Still open only as the *slab* sweep, which would give thickness-over-`h` directly rather than resolution-at-fixed-field | A-010 drove it to zero, which was the ticket's job; a sweep over a slab of shrinking thickness would give the parametrised form | M-15 established it is a resolution effect rather than a topology one, and M-4 has counts at two resolutions on two fields. It decides whether Surface Nets is usable at game resolutions or needs A-010 first |
| O-11 | *(**ANSWERED at R-007, M-284**: the dual carries a fourth stage Marching Cubes has not — `emit_quads`, which walks every grid edge on all three axes and loads both endpoints **before** the sign test that would let it skip. It is `O(n³)` where the surface is `O(n²)`, runs at **IPC 0.72** against the rest of the mesher's 3.8–6.5, and is **82% of the cycles**. Every other candidate was excluded by measurement first — see M-279. The remedy is **A-023**, and it is not a vertex rule.)* **Why does the dual topology go superlinear in `n³` while Marching Cubes does not?** *(Half-answered at M-45: it is not one machine's cache hierarchy. Surface Nets degrades on Zen 3 too — 37.4 → 49.1 ns/sample — and the `Surface Nets/Marching Cubes` ratio is worse there than on the M5. What remains open is the mechanism, not whether the effect is real.)* | A profile or cache-miss counter at 192³ vs 256³. The cross-machine experiment is **done**; a second one would not add anything | The working-set hypothesis survives and is strengthened: Surface Nets gathers the four cells around each crossed edge with one stride `n²` cells apart, and that stride is architecture-independent, which is exactly the kind of cost that would reproduce across microarchitectures. Note both machines show a per-sample **spike at 128³** (M5 Surface Nets 9.35, Zen 3 Surface Nets 53.84 against 45.6 at 96³ and 47.3 at 192³) — a working-set effect at one specific grid size on two unrelated cache hierarchies, which is itself a clue nobody has followed. **Narrowed at R-005 (M-279), and the working-set hypothesis in this cell is wrong.** The counters are in: the gather is `O(n²)` and the cost is `O(n³)`, and a field with **no surface at all** costs the same to within 0.9%, so it is not the gather. Nor branches (they fall), nor allocation (0 page faults), nor the TLB. What is left is a 16% IPC decline on an instruction stream that is flat per sample — and at 16.7 M samples a **2.4× swing in misses moves the cycles by 0.4%**, so the miss column is not the driver either. The 128³ spike **is** resolved: 127³ and 129³ are normal and only 128 has a 64 KiB plane stride, so it is conflict aliasing, and it survives on the empty field. Residue is **R-007** |
| O-12 | **Is Marching Cubes unconditionally manifold now?** ✗15's only counterexample was the fan chord and A-015 removed it; the strict gate passes 8,000 generated cases where it used to fail on the first seed. But nothing proves a second mechanism does not exist | An exhaustive search over configurations spanning more than two cells — the two-cell sweep is exhaustive and the vertex-link case is not covered by it at all. Or a proof that a cell-local cycle triangulation plus shared face segments cannot produce a non-manifold **vertex** | The strict gate is now asserted, so if a second mechanism exists CI will find it on some future seed. That is the intended outcome: a failure there is a finding, not a regression, and the failing case would be the first example of whatever the mechanism is |
| O-13 | ~~**Pre-registered:** Marching Tetrahedra vertex count = **3.0× Marching Cubes**, converging from above~~ **Confirmed at A-003/M-001, exactly and including the convergence.** Measured on `sphere`: 33³ **3.036**, 49³ **3.026**, 65³ **3.003** — from above, onto 3.0 | *(closed)* | And M-52 supplies the mechanism the prediction did not need but turns out to have: the ratio is `4.0` in one octant and `2.0` across a sign change, so `2.992` is an average hiding a factor-of-two spread. That is why the shootout CSV carries every field |
| O-14 | ~~**Pre-registered:** Marching Tetrahedra symmetric Hausdorff at 64³ ≈ **2.6e-3**, about **1.86×** Marching Cubes, i.e. slightly worse than Surface Nets~~ **Falsified at A-003/M-001 (M-55): measured 1.4386e-3, which is 1.043×.** Not slightly worse than Surface Nets — **better by 1.6×** (Surface Nets is 2.251e-3, 1.69× Marching Cubes) | *(closed)* | The prediction's stated counterintuitive part, *"more vertices **and** worse accuracy"*, is the half that is wrong. Marching Tetrahedra buys 3× the vertices for 4% worse accuracy on smooth fields and **better** accuracy on sharp ones |
| O-15 | ~~Why does a plane cost `3.94×` and a sphere `3.00×` when both are locally flat at cell scale?~~ **Answered at A-003 (M-52): the normal's sign pattern, not its direction.** One octant gives `4.0` exactly, a sign change gives `2.0`, and the isotropic average is P-1's `2.992`. A plane has one normal and a sphere has all of them | *(closed)* | What remains is small and not worth a ticket: the mixed-sign measurement spreads `1.98–2.27` against a predicted flat `2.0`, so the continuum model gets the mechanism exactly and carries a discretisation term it does not describe |
| O-16 | **Can the parallel dual-edge collapse (M-59) be removed without giving up the cycle partition?** The dual is a manifold complex; the index buffer is where it stops being a manifold mesh. A finer split — one vertex per *face-segment adjacency* rather than per cycle — would separate the endpoints, but it is not obviously topology-preserving and would cost vertices on every field rather than the two that need them | Enumerate, over all `(case, joined)` pairs on both sides of a shared face, the configurations where both cells put both of that face's segments in one cycle. That is a finite two-cell sweep of the same shape as ✗17's, and it would say whether the collapse is rare-and-coarse-only or merely unobserved on the seven fields | Bounded in practice: zero on all seven reference fields at every tested resolution, one edge on the ✗15 fixture at `h = 2/3`, gone by `h = 1/2`. So it is the same "coarse grid does not resolve the surface" regime as ✗15 and ✗17, and the same answer probably applies — refine, or accept it and pin the count |

### Pre-registered predictions

Registering the predicted value *before* running the measurement is the point. A prediction that first
appears after the number is known is not evidence, and this project has already caught itself writing
expectations into docs that the measurement then disproved (✗1, ✗3).

**P-1 — Marching Tetrahedra / Marching Cubes vertex ratio = 3.0×.** Derived, then confirmed
numerically, before implementation. **Outcome at A-003: confirmed exactly, and the spread explained** — `sphere` 3.04, `torus` 3.04 against a predicted 2.99, and `box_exact` 3.91 with a bare plane 3.94 — because the ratio is `4.0` inside one octant and `2.0`
across a sign change, and `2.992` is the integral of those over the sphere. The derivation was not
just right, it was the whole answer; what it did not say out loud is that the average hides a factor-of-two
spread over individual orientations. See M-51 and M-52; O-15 is closed.

The 6-tet decomposition is the **Kuhn / Freudenthal triangulation** — Freudenthal 1942, Kuhn 1960,
decades older than Doi & Koide 1991 and far better documented. Verified: the six monotone 000→111
paths give six tets, each of volume exactly 1/6, summing to 1, all sharing the main diagonal.

Its tet mesh uses **7 edge families per cell** against Marching Cubes' 3 — 3 cube-axis (|e|=1),
3 face diagonals (√2), 1 body diagonal (√3). Weighting each by crossing probability `E[|n·e|]` over
uniformly random surface orientations: 4.4877 / 1.4999 = **2.992**. Counted directly on a sphere SDF:

| grid | Marching Cubes (3 axis) | Marching Tetrahedra (all 7) | ratio |
|---|---|---|---|
| 33³ | 1830 | 5582 | 3.050 |
| 65³ | 7470 | 22394 | 2.998 |
| 129³ | 30078 | 89978 | 2.991 |

**Falsified by** a converged ratio outside ~2.95–3.05 — which would indicate the implementation is not
emitting one vertex per crossed tet-mesh edge, i.e. an edge-cache or decomposition defect.

**P-2 — Marching Tetrahedra is *less* accurate despite having 3× the vertices.** Mechanism: linear
interpolation places a vertex at the zero of the linear approximation along an edge, so its error
scales with `|e|²`. Every Marching Cubes crossed edge has `|e|² = 1`. Weighting Marching Tetrahedra's
seven families by crossing probability gives mean `|e|² = 1.859`, worst case **3** (body diagonal).
More vertices, each individually further from the true surface. Against M-10's `1.380e-3`, that
predicts ≈ **2.6e-3** — above Surface Nets' `2.288e-3`.

This is consistent with M-12's measured `h²` convergence: the error term is `O(|e|²·κ)`, and Marching
Tetrahedra simply draws from a longer-edge distribution.

> **Scoping caveat — do not file this as confirming Lewiner.** Lewiner et al. 2003 says tetrahedral
> vertices *"cannot be adjusted to fit the geometrical trilinear approximation as we do with cubes"* —
> he is comparing against **his trilinear-fitted MC33**, not linear-interpolation Marching Cubes.
> A-003 therefore tests **"Marching Tetrahedra vs linear-interp Marching Cubes"**, the weaker half.
> The trilinear half needs A-002's decider *plus* trilinear-aware vertex placement to compare against,
> which A-001 does not do. Record the verdict against the comparison actually run.

**P-3 — the conforming property is assertable, not just assumable.** Kuhn tiles space face-to-face
**only if every cell picks the same main-diagonal direction**; the 5-tet decomposition needs
alternating orientation instead, which is why it is the wrong choice here. Verified: cell A's `z=1`
face is split by its local diagonal `001–111`, cell B's `z=0` face by `000–110`, and those are the
same world-space segment. Assert it directly — compute the shared face's diagonal from each side and
require equality — rather than trusting the construction. Same shape as `face_disagreements`.

**P-4 — the case table is 16 entries and needs no source paper.** A tet has 4 corners, so 2⁴ = 16 sign
configurations: 2 trivial, **8** one-vs-three (the isolated corner's 3 edges cross → 1 triangle), **6**
two-vs-two (4 edges cross → 2 triangles). Generate it and prove it exhaustively. Doi & Koide being
unobtainable costs nothing here — unlike A-001's 256 cases, this table cannot be mistyped because it
is not typed.

**P-5 — Manifold Dual Contouring carries Marching Cubes' Euler characteristic across exactly.**
Registered before running, from the construction rather than from a paper: the output is the dual of
this crate's Marching Cubes, and a dual has `V' = F`, `E' = E`, `F' = V`, so `χ' = χ` identically.
**Outcome at A-010: confirmed on every closed field at 17³, 25³ and 33³** — and the one place it fails
is what found M-59. On the ✗15 fixture at `h = 2/3` the dual reports `χ = 1` against Marching Cubes'
`0`, and the discrepancy is exactly the number of collapsed parallel edges. So the prediction held, and
its single failure localised a mechanism rather than merely being wrong.

**P-7 — `thin_plate` welds into exactly one connected component, and orienting it as one would
*not* flip a face.** Registered before A-014f is attempted, and before the measurement, because
A-014f's proposed remedy — orient each connected component from its most confident triangle and
propagate — is a per-patch decision, and M-96 is on record that `thin_plate`'s two faces sit 0.4 cells
apart and land in one tetrahedron facing opposite ways. The concern raised against it is that a weld
could merge the two sheets *through their thickness* and let propagation silently invert one of them,
on the one field that justifies the whole subgrid track (M-95: 4,328 triangles where greedy quads
returns zero).

The prediction is that this specific failure does not occur, for two reasons that should be checked
separately because they can fail separately. **(a)** The validity suite welds at `cell · 1e-6`, six
orders of magnitude below the 0.4-cell separation, so no through-thickness merge is possible at the
epsilon actually in use — but the margin, not the outcome, is what should be recorded, since it is
the quantity that would decide a coarser weld. **(b)** More importantly, `thin_plate` is **closed** in
its domain and reports zero boundary edges, so its top and bottom are joined at a rim and are
*already* one component by topology rather than by tolerance. A closed orientable surface has a
coherent orientation across that rim, so component-wise propagation is the right answer there rather
than a risk to it.

**What would falsify it:** more than one component after the weld, or a component count that changes
with the epsilon, or an orientable-surface check that fails on the welded plate. Any of those means
the remedy is unsafe on the field it most needs to be safe on, and A-014f must find another.

**Outcome — limb (a) confirmed at T-009, with the margin it asked for (M-181).** The concern was that
a weld could merge `thin_plate`'s two sheets *through their thickness*, and the prediction was that the
epsilon in use sits far below the 0.4-cell separation. Measured rather than argued: welding
`thin_plate`'s subgrid output at every factor from `h·10⁻⁹` to `h·10⁻³` returns **422 vertices at all
seven**, and the first tolerance that changes it at all is `h·10⁻¹` — **1000× the one policy** — where
it collapses to 128, which is the through-thickness merge finally happening. So the margin is three
orders of magnitude, not the "six orders below the separation" the prediction reasoned from, and the
quantity is now a number rather than an inference. The same test covers the second falsifier outright:
**the count does not change with the epsilon**, on any of the seven fields, across the whole range the
crate's four historical policies spanned.

**Outcome — limb (b) confirmed too, and both remaining falsifiers fail to fire (M-182).** Welded
`thin_plate` is **one component with 0 boundary edges and χ = 2 at 17³, 25³ and 33³** — closed, genus 0,
joined at a rim, one component by topology and not by tolerance, exactly as registered. The component
count stays 1 at every epsilon from `h·10⁻⁹` to `h·10⁻²`.

**So P-7 holds in full, and the remedy it was registered to de-risk is safe on the field that
justifies the subgrid track.** The prediction was worth making: it named its falsifiers in advance, it
was checked against them, and the answer arrived before a line of A-014f was written.

**The check also found something P-7 was not looking for.** `thin_plate`'s inconsistently-oriented
edges run `17³ → 0`, `25³ → 8`, `33³ → 6`, and nobody had seen it because
`the_validity_suite_over_every_reference_field` runs at `n = 17` alone while Phase 1's gate says three
resolutions. That is A-014f's defect on A-014e's protected field, and it *supports* propagation: the
surface is orientable and the local vote is failing to find the orientation that exists.

**P-6 — splitting the vertex reduces self-intersection on the two multi-sheet fields.** Registered
before running, reasoning from M-29 (*"the residue is exactly A-010's problem"*): if the clamp removed
every placement-caused intersection and the remainder is caused by two sheets sharing a vertex, giving
them separate vertices should remove the remainder. **Falsified at A-010 (M-61): it goes up**, `gyroid`
3.118 → 5.669 and `fbm_terrain` 13.837 → 15.434. The premise was right and the inference was wrong —
two vertices in one cell is exactly what breaks the within-cell partition the clamp's guarantee rests
on. ✗2's competing figure (ODC measuring Manifold Dual Contouring at 100% self-intersecting) was on
the record the whole time and should have been weighted against M-29 before registering this.


### P-8 … P-13 — registered before anything was written, and enforced by the compiler (R-000)

**The practice became a compile error.** P-1 through P-7 were pre-registered in prose, and the record
is good — P-6 was falsified and says so, P-7 held. The weakness was that nothing stopped an
experiment from being *written first and predicted afterwards*, which is the exact failure the
practice exists to prevent and is **invisible in the artefact**: a prediction recorded after the
numbers came in reads identically to one recorded before.

So the predictions now live in `crates/isomesh/src/experiment.rs`, and `isomesh::experiment!("P-n")`
**does not compile** for an id that is not in `PREREGISTERED`. Registering is a commit, so git carries
the ordering that prose cannot. `scripts/backlog_gate.sh` checks every registered id reaches this
file — the Rust is the source, this is the elaboration, and the gate stops the two from parting
company.

`Preregistration` has a `falsified_by` field and no way to omit it, which is the other half: a
hypothesis with no stated refutation is a description in the future tense.

| id | ticket | prediction | falsified by |
|---|---|---|---|
| **P-8** | R-001 | A weld gated on `Lk u ∩ Lk v = ∅`, leaving rejected pairs split, yields **exactly 0** non-manifold edges and vertices on all eight fields × all extractors, where the unconditional weld yields `N > 0` | The gated weld still producing non-manifold output — which would prove the surface link condition insufficient for index-buffer realisation, and is the more interesting result |
| **P-9** | R-002 | For buckets of ≥3 coincident vertices, at least one field yields **≥2 distinct outputs** across seeded permutations of within-bucket merge order | All permutations byte-identical on every field, meaning the k-way weld is confluent and no canonical order is needed |
| **P-10** | R-003 | Vertex inflation from gated-weld-plus-split is **< 1%**, and self-intersections per 1k are unchanged | Inflation > 1% (a real trade-off needing a stated policy), or self-intersections rising (M-93's duplication artefact returns) |
| **P-11** | R-004 | With one canonical `world_of_sample` rather than offset-and-add, seam cracks fall to **0 for all cell sizes** — not only powers of two — and M-73's hairline disappears with no change to Transvoxel | Cracks surviving canonical reconstruction, which localises the defect back in Transvoxel |
| **P-12** | R-005 | The dual's superlinear cost is the four-cells-around-a-crossed-edge gather at stride `n²`: cache-miss count per sample **rises with `n`** for Surface Nets and stays flat for Marching Cubes | Flat miss rates, pointing at branch misprediction or allocation instead |
| **P-13** | R-006 | M-66's non-convergent angle is **bounded below by the dihedral angle** of the feature, so it is a property of sharp edges rather than of resolution | The angle failing to track the dihedral prediction, which makes it a defect with a location |

Each experiment's numbers land in `docs/experiments/p-n.csv`, stamped with the git SHA, the machine
and the time — and flagged **`WORKING TREE DIRTY`** when the tree was not clean, because numbers that
correspond to no commit have to say so on the artefact rather than in someone's memory.


### P-8 — FALSIFIED, in both clauses, and the gate makes things worse (R-001)

**M.** Eight reference fields × seven extractors × a 2×2×2 block of independently meshed chunks at
`h = 4/35`, both welds run on the same meshes in one pass. `docs/experiments/p-8.csv`.

> **These numbers replace a first run whose fixture was broken (M-274).** The block was placed at
> `-2.0` with 8-cell chunks, which spans 1.83 of the 4-unit domain and clipped a corner off every
> field. It produced **no bucket with more than two members**, so the mechanism the first write-up
> asserted was about a configuration the run never contained. The block is now centred and 18-cell,
> and P-9 reports 6–32 buckets of `k ≥ 3` per configuration.

**Clause one — "where the unconditional weld yields `N > 0`" — is false in 47 of 56 configurations.**
`Welder` produces zero non-manifold edges and zero non-manifold vertices everywhere except nine rows,
all under **dual or subgrid** extractors:

| field | extractor | ungated e/v | gated e/v | rejected |
|---|---|---|---|---|
| `torus` | `dual_contouring` | 0 / 4 | 0 / 4 | 0 |
| `torus` | `manifold_dual_contouring` | 0 / 4 | 0 / 4 | 0 |
| `thin_plate` | `subgrid_marching_tetrahedra` | 4 / 6 | **0** / 135 | 143 |
| `fbm_terrain` | `surface_nets` | 1 / 2 | 1 / 2 | 0 |
| `fbm_terrain` | `dual_contouring` | 1 / 2 | 1 / 2 | 0 |
| `noise_cavity` | `surface_nets` | 276 / 536 | 276 / 536 | 1 |
| `noise_cavity` | `dual_contouring` | 277 / 539 | 276 / 537 | 8 |
| `noise_cavity` | `manifold_dual_contouring` | 66 / 145 | 65 / 143 | 9 |
| `noise_cavity` | `subgrid_marching_tetrahedra` | 222 / 301 | 218 / **1,092** | 928 |

**Clause two — "yields exactly 0" — is false, and the trade is catastrophic.** Across all 56
configurations the gate removes **at most 4** non-manifold edges (`thin_plate` + subgrid, 4 → 0) and
adds **up to 791** non-manifold vertices (`noise_cavity` + subgrid, 301 → 1,092). On every **primal**
extractor — where the weld was perfectly clean — it manufactures non-manifoldness from nothing:
`sphere` 0 → 96, `torus` 0 → 108, `box_exact` 0 → 102, `noise_cavity` 0 → 407.

**The mechanism.** `vertex_delta` equals `rejected_merges` in all 56 rows, so every refusal leaves
exactly one extra vertex — the inflation is arithmetic and harmless. The damage is that a coincidence
of `k` vertices is manifold **only when all `k` merge**: refusing one leaves the representative
carrying cones from some copies and not others, which is a bowtie — two cones sharing an apex, every
edge with exactly two faces, χ intact. That is precisely what `validate`'s link walk exists to catch
and what an edge count cannot see, and it is why the damage shows up in the *vertex* column while the
edge column barely moves.

So the pre-registration's stated falsifier — *"proving the surface link condition insufficient for
index-buffer realisation"* — has the verdict right and the reason half right. The **pairwise**
condition is not insufficient for a pairwise merge. It is being applied greedily inside a `k`-way
group, and Dey, Fan & Wang's decomposition into `k − 1` pairwise merges *in the intermediate complex*
is exactly what does not commute with rejecting one of them. Rejecting a merge is not a safe no-op —
it is a decision to split, and a split part-way through a group is worse than either whole.

**Two consequences downstream.** R-003 asked whether splitting the unsafe merges is free; it is not,
and the cost is not the vertex inflation it anticipated but new non-manifoldness. And the residue —
276 edges and 536 vertices on `noise_cavity` under `surface_nets`, which the gate rejects **one**
merge against and changes not at all — is an **extractor** defect rather than a weld defect.

### P-9 — HELD, and the driver is not `k` (R-002)

**M.** Eight seeded permutations of within-bucket order per configuration, re-welded, compared by
byte-identity. Same fixture as P-8. `docs/experiments/p-9.csv`.

**Nine of 56 configurations produce more than one distinct output**, so H holds: the `k`-way weld is
**not** confluent and `CLAUDE.md`'s byte-identity guarantee rests on the input order being
deterministic rather than on the weld being order-free.

**But the premise is wrong about why.** H said *"for buckets of ≥3 coincident vertices"*. Measured, the
correlation runs the other way:

| field | extractor | distinct | vertex spread | buckets `k ≥ 3` |
|---|---|---|---|---|
| `sphere` | `marching_cubes` | **1** | 0 | **6** |
| `box_exact` | `marching_tetrahedra` | **1** | 0 | **6** |
| `torus` | `dual_contouring` | **5** | 0 | **0** |
| `torus` | `manifold_dual_contouring` | **5** | 0 | **0** |
| `noise_cavity` | `dual_contouring` | **8** | **4** | 1 |
| `noise_cavity` | `manifold_dual_contouring` | **8** | **2** | 2 |
| `noise_cavity` | `marching_tetrahedra` | 6 | 0 | 21 |

Every primal extractor has 6–32 buckets of `k ≥ 3` and is **perfectly confluent**. `torus` under the
duals has **none** and yields five different outputs.

**The driver is whether the coincident vertices are bit-identical, not how many there are.** A `k ≥ 3`
bucket on a primal extractor comes from M-48 — a grid sample landing exactly on the isosurface, so
every cut edge meeting there places its vertex at the *same computed point*, within one chunk, by the
same expression. Those are bit-identical, so which survives cannot change a byte. A seam bucket spans
two chunks that computed the plane by different expressions, and M-32 measured those disagreeing by
an ulp — so which survives changes the output exactly.

**And two rows change the vertex *count*, which is the part that matters.** `noise_cavity` under
`dual_contouring` spans **4** vertices across the eight permutations, and under
`manifold_dual_contouring` **2**. Epsilon-closeness is not transitive and `Welder` uses first fit, so
permuting the order changes **how many representatives are elected**, not merely which. The weld's own
documentation predicted this in prose (*"a chain `a ~ b ~ c` with `a` and `c` further than `ε` apart
yields two representatives"*); this is the first measurement of it.

**What it means for a consumer, which is not what the ticket feared.** R-002 worried that gating would
break determinism; gating is reverted (P-8), so that threat is gone. The live one is that a chunked
consumer **appending chunks in a different order** gets a different vertex count on `noise_cavity`,
not merely different bytes. The order is deterministic today because the loops are; nothing enforces
it.

---

## Part 4b — Experiments run, and what happened to them

**The slot this file did not have (T-013).** `M-` records a measurement, `✗` a belief that died, `O-`
an open question and a `P-` prediction rides inside the `O-` row it settles. None of those is the
right home for *"we built the variant, measured it, and put it back"* — and with Phase 8's ablation
seam, most experiments will end that way. A reverted experiment whose numbers were never written
down gets re-run in six weeks by someone who has forgotten, which is the specific waste this section
exists to prevent.

**An entry is owed whenever an ablation runs**, whether it was kept or not. The verdict line is the
point; the numbers are what make it re-checkable.

```markdown
| E×n | **<the variant, in one line>** | H: <the hypothesis, as pre-registered> | <numbers, both arms> | **KEPT** / **REVERTED** — <why> | <ticket, harness> |
```

| # | Variant | Hypothesis | Both arms | Verdict | Where |
|---|---|---|---|---|---|
| E×1 | **Surface Nets' centroid as Dual Contouring's vertex rule** | The QEF is worth its cost on sharp fields and not on smooth ones, and pays for it in self-intersection | Hausdorff at 65³, QEF ÷ centroid: sphere **0.486**, torus **0.457**, csg_difference **0.255**, box_exact **0.010**, thin_plate **0.010**. Self-intersections per 1k at 33³: QEF **3.118 / 13.837 / 29.745** on gyroid, fbm_terrain, noise_cavity against centroid's **0.000** | **KEPT as an ablation, not as a default.** Both arms are real answers to different questions and neither dominates — 100× accuracy on sharp features against zero self-intersections. The seam stays so the comparison can be re-run; `Qef` stays the default because sharp-feature recovery is what A-007 exists for | X-002, `benches/ablation.rs`, M-237 |
| E×2 | **A separate probabilistic-quadric solver** (Trettner & Kobbelt, `10.1111/cgf.13933`) | It supersedes the Tikhonov regularizer and is more robust on near-singular cells | Never measured as a separate solver, because it was shown identical first: a direct assembly of the paper's equations agrees with `solve_with` at `λ = Nσ²` to **1.110e-16 over 296 cells** | **REVERTED before it was written.** In this crate's centroid-relative coordinates the paper's extra term is `σ²Σrᵢ`, and `Σrᵢ ≡ 0` because the centroid *is* the mean of the crossings. A second solver would have been a second execution path computing identical numbers. **Do not re-attempt for isotropic noise**; the open door is anisotropic `Σₙ`, which needs a noise model analytic fields do not have | X-004, `the_probabilistic_quadric_is_the_existing_solve`, M-238 |
| E×3 | **Crossing-count-scaled regularizer** (`λ = Nσ²`, the part of E×2 that *is* different) | Scaling λ with the number of planes beats one fixed λ per cell | Hausdorff at 65³, scaled ÷ fixed: sphere **1.0000**, torus **0.9957**, csg_difference **0.9992**, box_exact **0.7519**, thin_plate **0.7519**. Self-intersections at 33³ fall on all three noisy fields: **3.118→2.551, 13.837→13.571, 29.745→28.749** | **KEPT behind `experimental`, not made default.** Never worse and 25% better on both sharp fields, which is a real improvement — but the default carries T-007's committed golden hashes and 112 baseline rows, and moving those for a 25% gain on two of eight fields is a decision with evidence attached, not a tidy-up. Promoting it is its own ticket | X-004, `crates/isomesh/src/experimental.rs`, M-238 |
| E×4 | **Weld gated on the pairwise link condition, rejected pairs left split** | H (P-8): exactly 0 non-manifold edges and vertices on all eight fields × all extractors, where the unconditional weld yields N > 0 | 56 configurations on a centred 2×2×2 block. Ungated is **0/0 in 47 of them**. Across all 56 the gate removes **at most 4** non-manifold edges and adds **up to 791** non-manifold vertices — `noise_cavity` + subgrid goes 301 → **1,092**, and `sphere` + Marching Cubes goes 0 → 96 | **REVERTED, and it was never merged.** Strictly worse: it fixes almost nothing where there was something to fix and manufactures non-manifoldness where there was none. A `k`-way coincidence is manifold only if all `k` merge; refusing one leaves the representative a bowtie, which is why the damage is in the vertex column and the edge column barely moves | R-001, `benches/experiment_p8.rs`, `docs/experiments/p-8.csv`, P-8 |

---

## Part 5 — Method rules, and the failure that earned each

Rules with no incident behind them get ignored. These all have one.

| Rule | Earned from |
|---|---|
| **Derive a measurement's bounds from the thing being measured. A hand-written bound is wrong twice before anyone notices** | B-006 and B-007 — the same seam-counting helper had its exclusion list written out by hand twice. The first version omitted the `y` axis entirely and accused Marching Cubes of a seam defect it does not have; the fix for that left `z` with a bound of `8.0` on a chunk `4.0` deep, and the uncounted wall produced a phantom open edge that made subgrid look non-conforming and got shipped as `Unverified` in a public API. Both readings were plausible and both were the fixture. The bounds now come from `layout.cell_size() * layout.cells()`, which cannot disagree with the chunks it is measuring |
| **A counter that is only populated on the success path cannot report the failure. Register what is being checked *before* the thing that fixes it** | E-205 — the crack counter built its list of seam planes inside `if transitions`, so running with transitions **off** left nothing to compare against and the demo reported a confident **0 cracks** on a world with 182 open edges in it. The zero was not wrong about the geometry; it was computed over an empty set. Moving the seam-plane scan out of the conditional makes the control real: 71 low and 102–111 high with transitions off, 0 and 0 with them on. **Seventh instance in one session** of a number that was a property of the fixture rather than of the code |
| **A trap that has to be dodged by choosing a constant will be walked into again. Put the dodge in the step, not in the default** | E-104 found that `box_exact` is exactly zero on its whole boundary, so a grid aligned to the box faces lets the sign convention decide instead of the algorithm — over the ±2 domain, whenever `n − 1` is a multiple of 4. E-114 then defaulted to **13**, which is one of them, and opened on the degenerate case it exists to explain: corner 7 sampled `-0.0000`. E-104's own note says the dodge belongs *"in the code rather than left as a warning"*, and it was — **in E-104's code**. The fix that survives is arithmetic rather than vigilance: step by 4 from an odd base so `n − 1 ≡ 2 (mod 4)` throughout and no reachable resolution is aligned |
| **Frame the camera from the field's own extent. A fixed orbit radius ships a screenshot of the *inside* of the mesh, and it looks like a rendering rather than a framing bug** | E-110, looking at E-109's committed image for a style reference. `sharp_features.rs:131` sets `orbit.radius = 7.0` for every field, and the capped gyroid's domain extent is **14** — so `docs/screenshots/e109-sharp-features-gyroid.png` is a picture of an inner wall, with the HUD legible over a flat beige surface and none of the geometry the demo exists to show. `manifold_check.rs:256-261` already had the fix and the comment explaining it (*"a fixed radius puts the camera comfortably inside the gyroid"*), written one ticket earlier; E-109 did not reuse it. **The image passed review because it is not obviously wrong** — nothing is missing, nothing is inverted, the numbers are correct. That is the failure mode: a framing bug produces a plausible picture, so only comparing against what the field *should* look like catches it |
| **A falsifier has to separate the hypothesis from its rivals, not merely be capable of failing** | M-279 — P-12 registered *"flat miss rates"* as its refutation. Miss rates were not flat, so by its own falsifier the hypothesis survived, and its mechanism was false anyway: *"the crossed-edge gather misses"* and *"the `O(n³)` scan misses"* both predict rising misses on a growing grid, so the stated observation could never have told them apart. The control that settled it — the same sweep on a field with **no surface**, where the gather runs zero times and the dense state is unchanged — was not in the registration. **When registering, ask what else would produce the same reading** |
| **A control run where it cannot discriminate reports "no effect", and reads exactly like a real negative** | M-279 — the axis-order control was first run at 4.3 M samples, a 17 MB array inside this machine's 32 MB L3, where no traversal order can miss. All three orders came out within 20% and the honest-looking conclusion was *"orientation does not matter"*. Re-run at 16.7 M it is a **2.4× spread**. The fixture is now run at **both** sizes, the small one as a control on the control. This is G-003's rule at a different scale: there the fixture's *value* sat in the degenerate region, here its *size* did |
| **Check a new harness against a committed measurement of the same thing before believing its new columns** | M-279 — `experiment_p12` forgot `MeshBuffer::reset()`, so the output buffer grew by a whole mesh every run and later runs paid reallocation the extraction did not cause. Every exotic column looked plausible; the tell was the boring one, `triangles`, which was **not monotone in `n`** — 145900 at 112³, 190060 at 128³, 144708 at 144³. `resolution_sweep-ryzen9-5900x.csv` has 5180 triangles at 48³ and the fixed harness reproduces it exactly. **A new instrument's first job is to agree with the old one where they overlap** |
| **On a governed CPU a nanosecond is not a unit. Report cycles, and put the clock on the row** | M-280 — the same binary reported Marching Cubes at 48³ as 8.13 and 14.66 ns/sample with cycles/sample unchanged at ~34, because `amd-pstate-epp` on `powersave` spans 1.96–5.62 GHz. Nothing on the face of either number said which clock it was. Every row now carries `ghz`, computed as cycles ÷ nanoseconds, so the artefact states it rather than inviting the inference |
| **A reference implementation used as ground truth needs the same scrutiny as the thing it checks — and when a measurement is impossible, suspect the instrument before the world** | M-289 — R-006 and R-008 both compared a mesh against an analytic gradient, and every control they carried was about whether the *mesh* was being measured fairly: rotate the fixture, offset the apex, add a no-crease case, check the vertex is on the surface. **None asked whether the gradient was right**, and it was wrong at exactly the points being measured — normalising a cancellation residue for any point epsilon-outside the surface, which is about half of them. Two hypotheses were reported falsified and both were true. The tell was there from the start and was misread twice: an area-weighted normal cannot leave the cone its faces span, so a past-90° reading was **arithmetically impossible**, and M-283 recorded the impossibility and went looking for strange geometry instead of a broken instrument |
| **A fixture can exhibit the property too perfectly. Exactness is the tell, not the confirmation** | M-283 — a wedge whose bisector lay on a grid axis reproduced P-13's predicted angle to **four decimal places** at three dihedrals and three resolutions: 75.0000, 60.0000, 45.0000. Turning the same wedge 17° or 37° about its own crease gives 20.1–128.0° for the same dihedrals. The exact agreement is a property of the symmetry between the crease and the sampling — the worst vertex is then the symmetric one — and nothing about `75.0000` at three resolutions looks like an artefact. *(The numbers in this row were first taken from the run M-289 corrects; the aligned fixture's four decimals are unchanged by that correction, which is the part the rule is about.)* Part 5 already says to *search* for a fixture that exhibits the property; this is the other half — **when a measurement matches a prediction exactly, vary the fixture's orientation before believing it** |
| **A millisecond is a property of the binary. Compare within one build and one run, or compare ratios** | M-281 — two of this repo's benches measured Marching Cubes on the same field at the same resolutions with the same median rule and disagreed by a **uniform 1.24–1.36×**, including at 16³ where the whole run is 40 µs. Both loop shapes in **one** binary are identical (0.991–1.002), and adding **one unrelated function** to `resolution_sweep.rs` moved its own 256³ row from 152.5 to 130.8 ms. Layout bias, with a paper — Mytkowicz et al., ASPLOS 2009 (`10.1145/1508284.1508275`). `benches/layout_bias` is the standing check, and it asserts rather than prints |
| **A generator that recognises one shape stops counting when the shape changes — and its staleness check cannot see that** | M-277 — `findings_index.sh` matched Part 2's table rows and Part 1's `✗` headings. Measurements became `###` sections at M-255 and **twenty-two of them fell out of the index while `--check` stayed green**, because a staleness check compares the file against the generator and the two agreed. The count read *"249 measured"* against 271 present. **Check a generator against its source's own vocabulary, not against its previous output** — one `grep -c '^### M-'` beside `grep -c '^| M-'` is the whole test, and it is the check that was never written |
| **A file that records state drifts from the state unless something checks it. Write the check the second time, not the third** | E-113 — the demo shipped at `a0859e8`, the README referenced it, and its row sat in `BACKLOG.md` for four more commits; the header counts drifted from the row counts separately. Both were found by audit, both were caused by editing one file with several scripts in one turn where a later write clobbered an earlier one. The same defect was sitting in `FINDINGS.md` unnoticed: **O-1, O-2 and O-4 were still listed as open questions** after G-002, A-009 and G-003 had all landed and answered them. `scripts/backlog_gate.sh` is the check, and it is mutation-tested against all eight ways the two files can disagree — because a gate that has only ever passed is indistinguishable from one that cannot fail (M-44) |
| **An instrument that cannot report the failure has not reported the success. Show it producing a non-zero before trusting its zero** | E-208 — the paint-drift readout measured "did the colour at this point change", and the scripted run sprayed one colour, so repainting red over red was numerically identical to paint that never moved. It printed **0.000000 at every step**, which is the answer the ticket wanted, and it would have printed the same on an implementation that smeared. Cycling the palette turned the same instrument sensitive — **27 of 40 sprays register drift, up to 0.886** — and only then does the **0.000000 across both carves** mean anything. This is M-75's rule in a different costume (*"a test that returns the same answer when you invert the thing it is testing is not measuring that thing"*), and the reason it earns its own row is that here the instrument was not inverted, it was **starved**: the input never varied in the dimension being measured |
| A typed error at the call site is louder than an abort — make the invalid state unrepresentable where you can, report it where you can't, and never substitute a default | The no-panic rule, reconciled with "fail loudly": `ValidateConfig` has private fields and one checked constructor, so the validator needs no runtime guard at all |
| Corpus presence is decided by `catalog_read`, never by `distill_search` | ✗4 — 342 documents readable but unsearchable |
| **`for_each_reference_field!` looks like a closure and is not. A `return` in its body exits the whole test** | M-199 — the macro takes `\|name, field\|` and **inlines its body once per field as a plain block**, because each field is a different type and no single closure can take all seven. So `if name != "gyroid" { return; }` returned from the *test function* on `sphere`, and the test **passed while running neither control nor the assertion that both had run**. Use `if name == …` instead. The `continue`/`break` uses elsewhere in the tree are all inside genuine inner loops and are fine — a bare one would not compile, which is the only reason this trap is spelt `return` and nothing else. The macro's own doc now says so |
| **A precondition check belongs next to the action it guards, not at the top of the job.** Hoisted, it fires on every run that never reaches the action | M-198 — the publish job demanded `CARGO_REGISTRY_TOKEN` before running `publish.sh`, so a push that uploads nothing still needed a secret to stay green. `publish.sh`'s own header had already written the rule down — it is version-driven rather than push-driven because *"the alternative … would leave main permanently red and train everyone to ignore it"* — and the workflow reintroduced exactly that, one layer up, where the script could not see it. **The guard is not weakened by moving it**: a release push with no token still fails loudly, and now it says which crate and version it was about to upload. Both directions were checked, because a guard nobody has watched fail is decoration |
| **When two paths are timed in one process, the one that runs second pays. Re-order before believing either number** | M-197 — batching `extract`'s read-backs appeared to make `extract_buffers` **75% slower** (0.54 → 0.95 ms at 129³), consistently across four runs, on a path the change does not touch. Collapsing the two staging buffers into one allocation did not move it, so it was not allocation count. Swapping the order of the two measurements in `gpu_vs_cpu` put `extract_buffers` back at 0.553 and made the *read-back* path slow instead — whichever runs second pays, because the first leaves the device allocator warmed to a different peak. The real comparison was sound (both builds measured the read-back path first) and the apparent regression was the harness. **A one-line reordering is the cheapest way to find out, and it is worth doing before writing any number down** |
| **A property that holds on every fixture you have is a property of the fixtures until one of them can fail it** | M-208 through M-213 — five pre-registered claims (Manifold Dual Contouring's manifoldness, P-5's dual-χ identity, the weld's never-adds-a-defect, P-7's weld plateau, orientation's monotonicity) were all true on seven reference fields and all false on the eighth. Not one was wrong about the seven; every one was **narrower than it was stated**, and the missing precondition was the same in every case: *no field in this crate could produce a cell with an interior ambiguity* (M-208). **The tell was available in advance and nobody looked**: each claim is about what happens on a hard configuration, and no one had asked whether any fixture reached it. Before trusting a green suite as evidence for a property, ask what the fixtures *cannot* express |
| **A deleted repository is not a missing source. Ask an archive before taking a rule-5 stop** | V-31 — `github.com/rogrosso/tmc`, cited by both Grosso papers, 404s. **Software Heritage had a full git snapshot**, and `mc/MarchingCubes.cpp` reads out of it intact by content hash. The Wayback Machine had only the root page and would have answered nothing, so "I checked the archive" has to name *which* — Wayback archives pages, Software Heritage archives repositories. This mattered rather than being tidy: the paper's own listing of the inner hexagon is corrupt in the copy held here, and the program is the only surviving source for it. **Sibling of the rule below, and the same failure one ticket later**: the question was never "does the text say it" but "does any artefact say it" |
| **An acquisition status is a measurement with a date on it. Re-check it before it justifies a decision** | V-29 — `meshing-library-target.md` recorded Grosso 2016 and 2017 as `PAYWALL`. Both were in home-still, converted and indexed, acquired by the *same sweep* that supplied the papers A-002c was built on. The stale row was trusted for four days and it had sized a ticket: A-002b's largest single piece, a grid-subdivision preprocessing pass, exists only on the route the unread paper replaces. **A `PAYWALL` row is evidence about one attempt, not a property of the paper** |
| **A paper is not only its text. Check its ancillary files before concluding it does not state something** | V-28 — A-014d sat rule-5 stopped for two tickets on *"the triangulation patterns appear nowhere in the text"*, which was true and not the question. The arXiv `abs` page lists `anc/` files in its HTML; this paper ships a supplement **and a standalone JS implementation with a "Simplicial Embedding" checkbox** — §3.2.3, from its authors, executable. One `curl` at any point would have found it. The stop was correctly *taken* — inventing a case table is the failure rule 5 exists to prevent — and incorrectly *maintained*, because "the prose does not say" was never checked against "the artefact does not say". **Figures and code are sources; only text had been searched** |
| **Sample the pixels rather than trusting the eye. A figure read by eye is a measurement with no error bar** | V-28 — Figure 15's whole meaning hangs on which edges are drawn blue, and at a glance the pentagon's bottom edge looks navy against a blue-tinted fill. Sampling the rendered PPM along each edge gives blue-excess **+109** on the two real chords and **+14** on the one that merely looked it — a 7× separation that settles in one command what staring at the panel could not. The same read also caught a segment sampling as *fill* rather than as a line, which is how the boundary was shown to bend at an inserted vertex |
| **Search home-still before believing a doc that says a paper is missing** | M-63 — `docs/research/` lists Manifold Dual Contouring and the Transvoxel dissertation as "genuinely absent, blocking"; both are in the corpus and both were read in one session. A-010's ticket named the wrong paper for its own algorithm because nobody had opened the right one |
| **A test that gives the same answer when you invert what it tests is not measuring it** | M-75 — two winding tests on a zero-width transition patch reported the identical count in both fan orders, because the patch was exactly perpendicular to the gradient they were dotted against. Flipping the code and re-running is the cheapest check that a test has any power at all |
| **Assert the property you believe, not the one that is easy — and when the assertion fails, the counterexample is usually the finding** | M-64, and M-59 before it. "A lateral link always crosses the resolution boundary" was false, and the case that broke it is precisely what transition cells exist to do. "Manifold Dual Contouring is manifold" was too strong, and the case that broke it is a second mechanism nobody had named |
| **Never guess a DOI or arXiv ID.** Look it up or stop | A subagent guessed an ID from memory and downloaded an unrelated condensed-matter physics paper under a meshing DOI |
| Verify the *source* separately from the *number* | ✗7 — right figure, wrong attribution |
| "Nobody measured X" needs the same evidence as any other claim | ✗9 — asserted twice, false twice |
| No performance number without the benchmark that produced it, in the repo | The corpus contains several published figures that failed verification |
| A doc comment the test suite disproves is worse than no doc comment | ✗3 |
| Assert the identity, not the inequality — a weak assertion hides a strong fact | ✗1 |
| **Verify that a property test can actually fail.** Corrupt an input and confirm red | A test that cannot fail is decoration, not evidence |
| Record an assertion's break conditions *next to it* | ✗1 — G-001 chunking will break it correctly, and it will look like a regression |
| A command queued against a streamed entity must tolerate the entity being gone | M-171 — one example panicked on every run for want of `try_insert`, and three more were one busy frame away from it |
| **An accepted CI failure is indistinguishable from a new one. Encode the acceptance — skip, gate, or split the job — and keep the checkmark green, or the checkmark stops meaning anything** | M-174 — the ubuntu test job's known GPU red concealed a broken release gate, a false MSRV claim and a rule-2 gate that had been failing since GPU-001, for the whole GPU series |
| **A deliberate exclusion buys something and costs a gate. Write down which gates it exempts, or rediscover them one incident at a time** | M-190 — `bevy_isomesh` is excluded from the root workspace for a good reason, and rustfmt, rustdoc and MSRV each had to be patched in separately after the gap was found. Three incidents, one cause |
| **A suite that samples one point of a stated range is reporting on that point, not on the range. Check that each gate is run everywhere it claims to be** | M-182 — the subgrid validity census runs at 17³ while Phase 1's gate says three resolutions, and 17³ is exactly where `thin_plate` is clean. 25³ and 33³ carry 8 and 6 flipped edges, on the one field the whole subgrid track exists for |
| **Before implementing a remedy, measure the ceiling it could reach. A remedy whose maximum is known to fall short is a re-scope, not a task** | M-180 — A-014h's rule was implementable in an afternoon and would have removed 234 of `box_exact`'s 924 duplicates while its acceptance asked for all of them. One `HashSet` of bit patterns, written before the code, said so |
| **A claim that a property "needs no special case" is a prediction, and it gets a test like any other. Write the test that would fail if it did** | M-175 — "ties need no tie-break" sat in a doc comment for the whole of ✗12 and was false the entire time. It survived 9600 measured trials because those trials were the one case where it happened to be true, and it fell to the first test written *about the tie itself* |
| **A shipped artifact is not evidence that the pipeline shipped it. Check the pipeline's own record, not the outcome it was supposed to produce** | M-174 — 0.0.3 and 0.0.4 were on crates.io, which is precisely what a working `publish` job would look like from outside. It had been `skipped` on every run since it was written, and both versions went up by hand. **The finding that said so was itself wrong the first time in the opposite direction** ("never uploaded"), from reasoning about the gate instead of asking the registry |
| Pin known defects as non-zero assertions rather than excluding them | M-4 — the numbers only move when someone means them to |
| Single-grid timings measure dispatch latency; sweep resolution and report the fixed cost | V-6 |
| Treat any published cross-paper ratio below ~2× as noise | V-7 |
| A green local run on one platform is not a green build. CI is the first real test of anything platform-shaped, and it will find things a local pass structurally cannot | First push: every job passed except `bevy_isomesh` on Linux, where Bevy 0.19's default Wayland backend needs `libwayland-dev` / `libxkbcommon-dev`. No such package exists on macOS, so no amount of local verification could have caught it |
| **A ticket's acceptance criterion is itself a claim about the code. Check it against the code before starting the ticket, not after.** | ✗11 — A-002 carried an `L`-sized acceptance criterion that the existing test suite had already made unsatisfiable. Nothing flagged it, because acceptance criteria are read as instructions rather than as assertions to verify |
| A property that falls out of *how a table was constructed* outranks folklore about the algorithm the table implements | ✗11 — "Marching Cubes produces holes" is true of a transcribed table and false of a derived one; the distinction is invisible if you reason about "Marching Cubes" rather than about this code |
| **When a ticket paraphrases a research doc, re-read the doc.** A paraphrase can invert the property that made the technique worth adopting | ✗12 — "branch-free, handles all degeneracies" became "falls through when the triple product is near zero" across three documents, turning the rule's central guarantee into its opposite |
| **When an acceptance criterion passes by two orders of magnitude, it is not the test — find the one that fails.** Ship it anyway, and ship the real one beside it | M-11 — T-003's stated criterion passes with 80× margin, so a constant-returning harness satisfies it. The convergence-order test and the closed-form fixtures are what actually constrain the code |
| Estimate a count from the geometry, then **measure it before writing it down** — the tidy formula is usually missing a constant | M-13 — `A/h²` under-predicted the triangle count by 1.47×, because a surface crosses `3/2` cells per unit area, not one |
| **A test double is evidence only if it is pinned to the thing it stands in for.** Write the equivalence test first, and make it bit-exact | T-005b — the case-table mutation check runs a corrupted table through a local marcher. Without `the_double_reproduces_marching_cubes` comparing the two bit-for-bit on the *uncorrupted* table, a corrupted-table failure would be indistinguishable from the double having drifted |
| **Before waiving a property in a gate, check what else was resting on it.** Derived checks fail with their premise | M-16 — a gate that waived manifoldness kept the even-`χ` assertion, which is a *corollary* of manifoldness. It failed on the first non-manifold mesh it saw, and the assertion was the bug, not the mesh |
| When a mutation test passes for the wrong reason, the message tells you: **check where the panic came from, not just that there was one** | M-17 — the wrong-edge corruption first tripped `edge_crossing`'s precondition deep inside the crate, so the validity gate under test was never reached. `should_panic` alone would have called that a pass |
| **Write the prediction into the benchmark before the first run**, in the file, committed. Then the result cannot be rationalised afterwards | M-19 — T-006 predicted `a ≈ 0` on the CPU path from the fact that V-6's 73% figure is GPU dispatch overhead. It came out at 0.61%. Had it come out large, the prediction being on record is what would have forced the awkward question instead of a tidy story |
| **A fitted coefficient means nothing until it is compared to the data's own range.** Report it against both ends | M-19 — Marching Cubes' `a` is 0.61% of the largest run and 543% of the smallest. Either number alone tells a different and misleading story; the pair says "negligible at scale, do not extrapolate below the range" |
| **A physically impossible fitted parameter is the model telling you it is wrong.** Do not report it as a value | M-21 — Surface Nets' fitted fixed cost is *negative*. Reported as "there is no fixed cost" it would be nonsense; read correctly it says the cost grows faster than `n³` and the whole two-term model does not apply |
| **A property that has held in every measurement so far is still a hypothesis, not a mechanism.** Say which condition it depends on | ✗15 — "Marching Cubes is manifold" held on seven reference fields at every resolution ever tried, and the mechanism offered for it (vertices on edges, not one per cell) was real but insufficient. The true condition is "the grid resolves the surface", which nothing had stated |
| **A fixture chosen by intuition can sit in the degenerate region where the property being tested does not apply.** Search for one | G-003 — the smooth-min associativity fixture used values `0.5` apart with `k = 0.4`; past `|a − b| ≥ k` the blend saturates and smooth-min *is* just `min`, which is associative. The test asserting non-associativity was exercising the associative case. **This is the second time in two tickets** (M-32 was the first), which is why it is now a rule rather than an anecdote |
| **Count what changes in the *output*, not what changes in the input.** They can differ by an order of magnitude | M-34 — E1 first counted cells whose samples moved and read 100%, which says incremental meshing is pointless. Counting cells whose *triangles* move gives 15–36%, which says the opposite. An SDF edit perturbs a whole solid; it re-shapes only a shell |
| **Choose a test fixture by searching for one that exhibits the property, not by picking one that looks like it should.** | M-32 — the non-power-of-two seam test first used `h = 4/33`, which *looks* irregular and lands in the 78% of cases that happen to agree exactly. It passed while proving nothing about the case it was named after. A search over `(origin, h, cells, chunk)` found the 22% that disagree, and the fixture now comes from that search with an assertion that the two expressions still differ |
| **A workspace that is excluded from the root is excluded from the root's CI commands too.** Check each one separately | E-111 — the lint job runs `cargo fmt --all --check` from the root, which excludes `bevy_isomesh`, so that crate's formatting had never been checked in 20 tickets and an example was committed unformatted. `cargo check`/`clippy`/`test` had their own steps in the bevy job; `fmt` was the one nobody noticed was missing |
| **Count what the claim is about, not what is easy to count.** A bitwise "did it move" is not a measure of "was it wrong" | A-009 — the clamp box is the cell shrunk by `ε`, so a vertex sitting exactly *on* a cell face is nudged by `5e-5` cells and a bitwise count calls it clamped. On a grid-aligned `box_exact` that read as **1176 of 1352 vertices displaced**; the honest count, against a `1e-3`-cell threshold, is **zero**. The first number would have gone into a commit message and been wrong |
| **A remedy stated for one operation does not cover the pipeline.** If a property is claimed end-to-end, check every reduction in it | M-24 — the audit's "magnitude-sorted dot products" is real and insufficient; the determinant needed the same treatment, and nothing said so. The equivariance test caught it because it asserted bit-equality rather than a tolerance |
| **Read a dependency's API before believing what it is for.** Reputation is not a type signature | ✗16 — glam is "the" Rust math library and was written into four documents as A-007's dependency. It has no generic scalar, so it cannot serve a crate generic over `f32` and `f64` |
| Before believing a performance verdict, ask **how many machines it has run on.** One is a hypothesis | ✗14 — Surface Nets loses to Marching Cubes by 2.76× at 256³ on an Apple M5, and the mechanism is probably cache. That is a strong result and a weak generalisation until it runs somewhere else (O-11) |
| **When a new feature shows a defect, check whether the old one has it too before attributing it.** The cheapest version of that check is usually an exhaustive small search | ✗17 — the decider produced 2 non-manifold edges where plain Marching Cubes produced 0, which reads unambiguously as "the decider broke it". An exhaustive sweep over all 4,096 two-cell sign patterns found **12 affected under each rule** — the defect is A-001's fan, and the decider only changes which patterns are reached. Attributing it to A-002 would have put the fix in the wrong ticket and left Marching Cubes' own version of it undiscovered |
| **A measurement that comes back zero has to prove it could have come back non-zero.** Put the reachability check in the test | M-44 — the first chunk-seam sweep reported 0 decision flips and 0 ambiguous seam faces, which is a pass that means nothing. The assertion `ambiguous_faces > 100` failed and forced the sweep to be retuned until it actually reached the configuration. **Third occurrence of the fixture trap** (M-32, M-38), and the first where a test caught it rather than a reviewer |
| **Run every step of a CI job locally, not the ones you remember it having.** Name them in the definition of done so the list is not held in memory | A-002 — a public doc comment linked to a `pub(crate)` item, which `cargo doc` under `-D warnings` rejects and which clippy and fmt both pass. Two of the lint job's three steps were run locally and the third was not, so a green local run pushed a red CI. Same shape as E-111's missing `fmt` on the excluded workspace: the gap is always the step nobody thinks of as linting |
| **Implement the expensive fix, measure it, and only then look for the cheap one.** The measurement is what tells you a cheap one is worth hunting | A-015 — the ticket was written expecting to re-baseline ✗1, M-2, M-22 and all 84 golden hashes, and the naive centroid fix duly cost +73% vertices. That number was so much worse than the "a vertex and two triangles per long cycle" estimate that it forced the question "which chords can *actually* collide?", whose answer is local and made the fix free. Estimating the cost instead of measuring it would have shipped the expensive version or abandoned the ticket |
| **Record the margin, not just the verdict.** "It did not happen" and "it came within an ulp of happening" are the same count | M-44 — zero seam decisions flipped, but the number that makes that trustworthy is the closest observed margin, `1.535e-2`, against a perturbation of `~1e-16`. Without it, the zero could have been luck |

---

## Part 6 — Adding an entry

```markdown
### <tier><n> — <one-line claim>

**Believed because:** <where it came from — and be specific, "folklore" is a real answer>
**Tested by:** <the command, test name, or source read. Must be repeatable.>
**Result:** <the numbers>
**Consequence:** <what changed as a result — a decision, a ticket, an assertion>
**Would be shown wrong by:** <the observation that would falsify this>
```

**An experiment that was reverted still gets an entry.** Part 4b is the slot, and the verdict line is
the load-bearing part — a variant that was built, measured and put back is the single most expensive
thing to rediscover, because nothing in the tree records that it was ever tried.

Two things that make this file worth keeping rather than a chore:

**Re-tier rather than rewrite.** When an R becomes an M, leave the R text and add the measurement
below it. The gap between what was reported and what was measured is itself data.

**Record the ones we got right for the wrong reason.** ✗10 is in the falsified section even though
the outcome was fine, because the reasoning was wrong and the reasoning is what generalizes.

### M-255 — naive reinitialisation moves the zero set, measured before it was fixed (S-004)

**M.** Sussman & Fatemi's warning is not folklore and not small. A narrow-band reinitialiser that
seeds from the interpolated crossing, solves, and keeps the solved value everywhere drifts the
surface **0.152 of a cell over twenty applications** on `Sphere` at 33³, `h = 0.125`, `band = 3` —
tracked as the crossing fraction along all 386 cut x-edges, none of which appeared or vanished.

At roughly 0.0076 cells per application that is invisible per edit and visible after a hundred. In a
destructible game it is a wall changing shape while nobody edits it.

**Freezing the seeds inside one call does not fix it**, which was the first thing tried and the
reason this entry exists: the *next* call recomputes the seeds from the previous call's output, so
the freeze holds for one round and the drift resumes. The fix is to restore the **input** values at
every sample adjacent to a sign change, which makes the crossing fraction bit-identical rather than
close — measured drift is then exactly `0.0`, and the assertion is `assert_eq!(worst, 0.0)` rather
than a tolerance.

**Why a tolerance would have been the wrong assertion.** A per-application bound of 0.01 cells is
satisfied by a steady creep of 0.0076. The bound has to be on the *total after many applications*,
which is why the test runs twenty and not one.

`construct::tests::reinitialisation_does_not_move_the_zero_set`.

### M-256 — the narrow band's cost claim needs the march, not the sweep (S-004)

**M.** The ticket's premise is that reinitialisation should cost **edited surface area** rather than
chunk volume. Delegating to `signed_distance_field_swept` and discarding what falls outside the band
satisfies every accuracy assertion and none of that premise: fast sweeping visits every sample on
every one of its eight passes regardless of value, so bounding it saves the *update* and not the
*visit*.

Fast marching bounds exactly. It finalises in increasing order of distance, so the first dequeued
value above the limit proves every remaining one is too, and `break` there is exact rather than a
heuristic cutoff. Measured on `Sphere` at 33³, `band = 3`: **4,802 of 35,937 samples finalised —
13.4%**, against 100% for both unbounded constructors.

So `march()` takes a `limit` and `signed_distance_field_marched` passes `far()`. One implementation,
one path; the bound is a parameter, not a second algorithm.

The test asserts the share stays below 25% for the same reason M-255's assertion is on the total: a
solve that quietly touched the whole grid would pass the drift check and defeat the ticket.

### M-257 — the approximate GPU method beats both exact CPU ones (S-005)

**M.** Jump flooding (Rong & Tan 2006, `10.1145/1111411.1111431`) against `construct::
signed_distance_field` (S-001, "exact") and `construct::signed_distance_field_swept` (S-002), on the
four reference fields whose analytic value is a distance, at 17³/33³/65³. Worst error against the
**analytic** field, at 65³:

| field | jump flood | swept | exact transform |
|---|---|---|---|
| `sphere` | **0.02048** | 0.09271 | 0.06250 |
| `torus` | **0.02370** | 0.10527 | 0.06233 |
| `box_exact` | **0.08839** | 0.15861 | 0.10825 |
| `thin_plate` | **0.06250** | 0.19221 | 0.08927 |

The flood wins all twelve rows. `docs/measurements/jump_flood.csv`.

**Why, and why it is not a paradox.** "Exact" in S-001 means exact to the nearest *sample*, so it
quantises every distance to the grid — its error is pinned at exactly `h/4` on `sphere` and
`box_exact` at every resolution, which is the signature of quantisation rather than of algorithm.
The sweep is exact in the limit but its Godunov update is first-order. The flood is seeded from the
same sub-cell crossings the sweep uses and then measures a **true Euclidean distance** to that seed,
so on these fields the sub-cell seeding buys more than the 27-offset lattice restriction costs.

**So the ticket's stated acceptance was the wrong assertion, and it took two failures to see it.**
S-005 asks for "error against S-001". Written literally, the gate asserts the flood *agree* with a
CPU constructor — which asserts it reproduce that constructor's error. Both attempts failed on that
basis while the flood was the most accurate of the three:

- against S-001: `sphere` at 17³ disagreed by a full cell (0.250) while sitting 0.082 from truth
  against the transform's 0.250.
- against S-002: `thin_plate` at 17³ disagreed by 0.523 cells while sitting 0.250 from truth against
  the sweep's 0.407.

The gate is `flood_err <= xform_err && flood_err <= swept_err`, against the analytic field. The
agreement figures stay recorded, because *how far the GPU and CPU paths drift apart* is a real
question for a consumer meshing on one and colliding on the other — it is just not a correctness
bound.

**Rule 5 note.** The seeding is a deliberate departure from the paper, stated rather than slipped in:
Rong & Tan seed each boundary *sample* with its own position. Seeding the interpolated crossing is
what makes the comparison measure the flood rather than the seeding, and it is why these numbers do
not transfer to a textbook JFA implementation.

### M-258 — a `u32` followed by `vec3<u32>` is 32 bytes, not 16 (S-005)

**V + M.** std140 aligns `vec3<u32>` to 16, so `struct { stride: u32, pad: vec3<u32> }` puts the pad
at offset 16 and the struct is **32 bytes**. Sizing the uniform buffer at 16 produced
*"the buffer bound at binding index 1 is bound with size 16 where the shader expects 32"* — caught by
wgpu, but only at dispatch, not at pipeline creation. `struct { stride: vec4<u32> }` with `.x` used
is 16 and says the same thing. The existing `GridParams` shader header already carried this rule in
prose ("two vec4s rather than a struct of scalars so std140 and std430 agree"); it was not followed
in the new file.

### M-259 — the round trip the crate did not have (S-006)

**M.** Sphere → Marching Cubes → `signed_distance_from_mesh` → Marching Cubes, on a 33³ grid at
`h = 0.125`. Vertex and triangle counts are **identical** across the trip (1,158 vertices, 2,312
triangles), χ is 2 both times, boundary and non-manifold edges zero both times, and the worst
`|analytic field|` over the output vertices goes from `0.001701` to `0.003803` — 2.2× degradation,
against a bound of half a cell (`0.0625`).

The **sign** agrees with the analytic field at every sample more than half a cell from the surface —
zero disagreements out of 35,937. That is Bærentzen & Aanæs Theorem 1's whole claim, tested rather
than trusted, and it is the assertion that catches a pseudonormal read from the wrong feature: a face
normal where a vertex one was needed is correct for most samples and wrong near a crease.

This is the first end-to-end test in the crate. It exercises the extractor, the pseudonormal sign,
the closest-point classification and the acceleration structure *against each other* — a fault in any
one shows up as geometry that moved.

**It needed a piece that did not exist.** Every constructor in `construct` returns a `Vec<R>` and
every extractor consumes an `Sdf`, so there was no way to mesh what had just been built.
`construct::SampledField` is that adapter, and it interpolates **trilinearly** rather than by nearest
sample because Marching Cubes' case table is derived from the trilinear interpolant — a field
wrapping the same samples any other way disagrees with the mesher about where the surface is.

### M-260 — a uniform grid over the sample cells lost to a flat box reject, 3.9× (S-006)

**M.** Mesh-to-SDF needs a closest-triangle query per sample. The first implementation binned
triangles into the sample grid's own cells and searched expanding shells, stopping when the shell's
own lower bound exceeded the best distance found. It measured **3,457 ms** against an unaccelerated
scan's **892 ms** on 9,261 samples × 872 triangles — 3.9× *slower* than the thing it was accelerating.

**Why, and it is not an implementation detail.** Reaching radius `k` costs `O(k³)` bins. A sample at
the corner of a 21³ grid around a unit sphere is twelve cells from the surface, so it walks
essentially the whole grid before it finds anything — and most samples in any grid are far ones.
Fixing the shell iteration from `O(k³)` to `O(k²)` per shell does not change the conclusion, because
the *total* over all shells is `O(k³)` either way.

Replaced with a two-level axis-aligned box reject: a box per triangle, and a box per block of 64
consecutive triangles. **394 ms — 2.3× faster than brute force**, and the whole `from_mesh` test
module dropped from 31.6 s to 3.0 s.

The reject is exact — a box whose nearest corner is beyond the current best cannot contain a closer
point — so the accelerated and unaccelerated paths are asserted **bit-identical**, not merely close.

**What makes the blocks work is a property of the input, stated rather than assumed:** Marching Cubes
emits triangles in grid order, so 64 consecutive ones are spatially close and their block box is
tight. A mesh whose triangles arrive in arbitrary order degrades to the unaccelerated scan and stays
correct. This is **not** a BVH, which is what the paper uses and what would make the query `O(log m)`
rather than `O(m)`.

### M-261 — `Real` gained `acos`, and `libm` is why (S-006)

**V.** The angle-weighted pseudonormal weights each face normal by the incident angle, which needs
an arc cosine `Real` did not have. Added as `libm::acosf` / `libm::acos`, unconditionally, for the
reason in CLAUDE.md's `libm` justification: `std`'s `acos` is the platform's and differs between
macOS and Linux, and T-007's 63 golden hashes are committed. A platform-dependent angle would make
this crate's mesh-to-field path produce different geometry on the dev machine and in CI.

### M-262 — the winding number beats the pseudonormal on holed meshes, by a widening margin (S-007)

**M.** A meshed sphere at 17³ with a cap removed, scored against the analytic sphere on the 4,491
samples more than one cell from the surface. Four hole sizes:

| cut | triangles removed | boundary edges | pseudonormal wrong | winding wrong | mean `\|w−½\|` |
|---|---|---|---|---|---|
| 0.6 | 104 | 24 | 5 | **0** | 0.0000 |
| 0.3 | 212 | 28 | 131 | **27** | 0.0674 |
| 0.0 | 268 | 28 | 327 | **50** | 0.1225 |
| −0.5 | 396 | 28 | 1,435 | **88** | 0.2650 |

The margin widens from 5-versus-0 to **16×**. On the closed mesh the winding number is exactly
`1.000000000` inside and `0.000000000` outside — the correction term is identically zero when there
are no boundary edges, so that case reduces to plain ray parity and calibrates the intersection code
separately from the construction.

**The one-hole version of this test was not evidence.** At `cut = 0.6` the two methods differ by five
samples out of 4,491, which is indistinguishable from noise in a discretisation. The sweep is what
makes it a finding rather than an anecdote.

**Why the winding number is not perfect either, and why that is correct.** The `mean |w−½|` column
grows with the hole. That is not the measure becoming confidently wrong — it is the *question*
becoming wrong. Once half the sphere is deleted the mesh no longer encloses the points under the
hole, so calling them outside is the right answer for the surface that actually exists, and scoring
against the analytic sphere is scoring against geometry that was removed. Jacobson et al.'s framing
is that a GWN measures *"how confident we can be that a point is inside"*; the column is that
confidence at exactly the samples this counts as wrong. So the assertion is a **3× margin over the
pseudonormal at every hole size**, plus exactness on the smallest hole, not zero errors everywhere.

**Sources checked rather than assumed.** Xie, Hafner & Wojtan (`10.1145/3811339`) give the
construction verbatim: `w_M(q) = Σᵢ sgn(r·nᵢ) − (1/4π) Σⱼ Ωⱼ`, with the cone apex *"directly behind
the ray"* so the cone contributes no forward intersections. They are also the source for **not**
citing Barill et al. 2018 as state of the art: Barnes–Hut summation *"trades off accuracy for
computational speed… the resulting values are merely approximations."* Martens & Bessmeltsev supply
the grid optimisation used here — *"to compute voxelizations of resolution N³, we only need to shoot
N² rays"* — which is why this casts one ray per row and sums the intersections beyond each sample.

### M-263 — the boundary must be counted with multiplicity, not as a boolean (S-007)

**V.** The closing cone needs one triangle per *net* directed boundary edge. Treating "is this a
boundary edge" as a boolean is correct on a manifold mesh with boundary and wrong on the triangle
soup this exists to handle: an edge with three incident faces has a net of one and needs one closing
triangle, while a boolean either drops it (if any partner exists) or double-counts it. The net is
also what makes the sign right without a separate orientation pass — a net of `+n` for `u → v` means
`n` copies of the triangle `(v, u, apex)`, the reverse edge, which is what makes `M + C` consistently
oriented.

`Real` gained `atan2` for van Oosterom & Strackee's solid angle. The four-quadrant form is
load-bearing: the denominator goes negative for a triangle subtending more than a hemisphere, and a
plain `atan` of the quotient loses the half-turn.

### M-264 — the uncertified set is a curve, not a resolution failure (T-015)

**M.** Plantinga & Vegter's per-cell condition evaluated on all eight reference fields at 17³, 33³
and 65³. Share of **active** cells carrying the isotopy guarantee:

| field | 17³ | 33³ | 65³ |
|---|---|---|---|
| `sphere` | 100% | 100% | 100% |
| `torus` | 100% | 100% | 100% |
| `thin_plate` | 100% | 100% | 100% |
| `fbm_terrain` | 100% | 100% | 100% |
| `box_exact` | 72.97% | 86.98% | 93.62% |
| `csg_difference` | 69.26% | 85.66% | 92.75% |
| `gyroid` | 84.84% | 92.94% | 96.36% |
| `noise_cavity` | 35.89% | 48.59% | 80.02% |

`docs/measurements/isotopy.csv`.

**The interesting number is not the fraction, it is how it scales.** Halving `h` multiplies a
`d`-dimensional set's cell count by `2^d`, so two doublings give `2^(2d)` and
`d = log₂(growth) / 2`. Measured 17³ → 65³:

| field | active `d` | uncertified `d` |
|---|---|---|
| `box_exact` | 2.14 | **1.10** |
| `csg_difference` | 2.17 | **1.13** |
| `gyroid` | 2.07 | **1.04** |
| `noise_cavity` | 2.33 | **1.49** |

The active set is two-dimensional, because a surface is. **The uncertified set is one-dimensional.**
So those cells are not scattered failures of an under-resolved field — they trace a *feature curve*:
the box's edges, the CSG seam, the gyroid's high-curvature ridges. No spacing will certify them,
because the feature is genuinely not smooth, and the fraction climbing toward 100% is just the 1-D
set thinning against a 2-D one rather than the problem going away.

`noise_cavity` at 1.49 sits between the two, which is the honest reading of a field that is *both*
under-resolved at 17³ and creased: its fraction jumps 36% → 49% → 80%, far faster than the others.
That is consistent with M-244, where its gradient was measured at 7.73 against a first draft's
declared 2.598.

**The certificate is exact here, and the reason it can be is worth recording.** The general condition
needs interval arithmetic over an arbitrary `F`, which this crate cannot do — an `Sdf` returns point
values, and a sampled gradient hull *underestimates* variation, so the predicate would pass where the
truth fails. That is the one direction a certificate must never err in. But the surface Marching
Cubes approximates is the **trilinear interpolant**, and for that the bounds are exact and
closed-form: `∂F/∂x` is bilinear in `(y, z)`, hence a convex combination of the four `x`-edge
differences, so its exact range is their min and max. `F` is a convex combination of the eight
corners, so `0 ∉ □F(C)` is exactly *"the cell is inactive"* — the first clause is free and the
predicate reduces to the inner product. No interval library, no sampling, and `h²` factors out of the
sign test on isotropic cells.

**Two configurations this crate already cares about are refused, and correctly.** The paper's own
Figure 1 — alternating corner signs — is refused, which is the ticket's "engineered violator". So is
a face-ambiguous cell: a configuration whose topology depends on a tie-break is by definition not
determined by the corners, and A-002's entire series exists because of it.

### M-265 — decimation *is* re-sampling on a nested grid, so the literature's rule has no bite here (T-016)

**M.** Four downsampling operators against re-sampling, on all eight reference fields, 65³ → 9³ in
four levels. `docs/measurements/downsample.csv`.

**Decimation matches re-sampling bit-for-bit** — every triangle count, every vertex position, every
error figure, on every field at every level. Asserted as an equality, not a tolerance. The reason is
elementary once stated: a grid that was sampled from the field, decimated, keeps a *subset* of those
same points, and decimation composes, so level 3 is still a subset of the field's own samples.

That falsifies the framing this ticket inherited. *"You do not downsample, you re-sample"* (Frisken,
Koschier) is sound advice about **filtering** operators; for decimation on a `2ᵏ + 1` nested grid it
is a distinction with no difference, and the free operator is already the recommended one. It stops
being free the moment a level is *edited* rather than sampled — which is the case a destructible
game is in, and is where the rule earns its keep.

**The filters are strictly worse, both ways at once.** `Mean` (box, `[1,1,1]/3`) and `Tent`
(`[1,2,1]/4`) produce **fewer** triangles *and* **larger** geometric error at every level — worst
`|analytic|` on `sphere` at level 3 is `0.07196` for the mean against re-sampling's `0.01430`, 5×.
And at level 3 the mean **deletes the entire torus**: 0 triangles against 128. A filter that removes
a whole object one level before the grid stops resolving it is not "correctly discarding a sub-cell
feature".

**`Min` behaves exactly as designed and the price is now known.** It never loses anything — more
triangles at every level on every field — because taking a minimum can only pull values negative, so
solid grows and never shrinks. The cost is a systematic dilation of about one fine cell per level:
worst error reaches `4.02` on `fbm_terrain` at level 3 against re-sampling's `0.93`, and `0.75` on
`sphere` against `0.014`.

**Recommendation, with the number behind it:** decimate. It is free, it is identical to re-sampling
while the fine level is field-derived, and both alternatives trade accuracy for something that is not
worth it — the filters for smoothing that deletes objects, `min` for a conservatism that costs
several cells of dilation.

### M-266 — M-72's aliasing is *alignment*, not chance (T-016)

**M, and it refines M-72 rather than contradicting it.** M-72 recorded `thin_plate` at
**4,088 → 1,016 → 248 → 56** triangles across LOD 0–3 and explained the survival as *"whichever edges
happen to straddle a thin slab still register a sign change"* — which reads as chance.

It is not chance. The plate is `0.4 × h` thick and centred at `y = 0`; every level here has an **odd**
sample count (65, 33, 17, 9), so `y = 0` is a sample plane at every level and the feature is exactly
aligned with the sampling. Shifting the plate by **half a cell** produces **zero triangles at every
level, including the finest**, where it is still only 0.4 cells thick.

So the sequence is not a sub-cell feature disintegrating; it is a perfectly aligned one being found
every time, and the falling counts are just its `x`/`z` extent being sampled more coarsely.

**And the ticket's attribution was backwards.** It predicted that 4,088 → 1,016 → 248 → 56 is what
*box-filter averaging* produces, and that re-sampling would make the plate *"correctly disappear"*.
Measured: that sequence is **re-sampling's**, and the box filter is what makes the plate vanish —
completely, from level 1. Both halves inverted.

The mechanism is clear once separated: re-sampling asks a **point** question, so a corner inside the
plate registers a sign change however thin it is; a filter asks a **neighbourhood** question, so a
thin negative sliver averaged against a large positive surround comes out positive. Neither is
"correct" — they answer different questions, and which one a consumer wants depends on whether a
feature that vanishes at a known distance is better than one that survives by luck of alignment.

**This also makes `ThinPlate::THICKNESS_IN_CELLS`'s doc comment wrong.** It reads *"Below `1.0` so
that no grid phase can ever put a corner inside the plate, with margin"*, and the canonical grid puts
a whole **plane** of corners inside it. Thin enough is not the same as off-lattice.

### M-267 — the sampled gradient supremum is not even monotone in sampling density (T-017)

**M.** `sup ‖∇f‖` re-measured at 8³, 16³, 32³ and 64³ samples over each field's own domain:

| field | n=8 | n=16 | n=32 | n=64 | still climbing |
|---|---|---|---|---|---|
| `sphere`, `torus`, `box_exact`, `csg_difference`, `thin_plate` | 1.00000 | 1.00000 | 1.00000 | 1.00000 | no |
| `gyroid` | 1.54383 | 1.72726 | 1.73003 | 1.73131 | no |
| `fbm_terrain` | 2.35846 | 2.84980 | 2.90822 | 3.04828 | **yes** |
| `noise_cavity` | 8.84299 | **7.74791** | 9.37619 | 9.72480 | **yes** |

**`noise_cavity` goes down and then up.** M-244 recorded its gradient reaching 7.734 and used that to
justify declaring the field `Unbounded`; the figure was measured at `n = 16`, and at `n = 8` it was
already **8.84** — higher. So 7.734 was not a lower bound that later runs merely improved on. A
sampled maximum over one point set is not nested inside a sampled maximum over a denser one, because
the denser set is not a superset — F-002's grid is offset by a fixed nudge and rescaled with `n`.

That is a sharper argument than M-244's for the same conclusion. M-244 said a sampled maximum is a
lower bound on a supremum and so cannot be declared as a Lipschitz constant. This says the sampled
maximum is not even a *stable* lower bound: no single run establishes anything, and taking the max
over several runs is the least one would have to do to claim a number.

**The five exact fields read exactly `1.00000` at every density**, which is the eikonal property
holding to five figures across four independent samplings — a stronger confirmation than any one run.

**The gyroid's declared bound is loose by 2×, and that is now recorded rather than suspected.** It
declares `Lipschitz { l = 2√3 ≈ 3.4641 }`, derived at M-244 from `|∂g/∂x| ≤ 2`. Measured supremum
converges to `1.731` — `bound_gap = 0.4994`. The derivation is sound and the bound is therefore
*correct*; it is not *tight*, because the three terms of `∇g` cannot be extremal simultaneously.
Sampling cannot settle whether `2√3` is attainable, which is why the declaration stands. The cost is
paid in F-005's empty-cell rejection, where a loose constant inflates the rejection radius and
disqualifies cells a tighter one would skip — M-248 measured exactly that, at 1.5× against
`thin_plate`'s 11.8×.

### M-268 — field quality is now in the regression gate, and every column is exact (T-017)

**M.** `docs/measurements/field_quality.csv`, baselined per machine and diffed by
`scripts/regress.sh field_quality`. Columns: `sup_grad`, `inf_grad`, `eikonal_pct`, `bound_gap`,
`certified_pct`, keyed on `(field, samples)`.

**All five are compared exactly**, with `None` tolerance, and the reason is structural rather than
optimistic: every one is a **max, a min, or a ratio of counts** over deterministic `libm` arithmetic.
None is a sum, so there is no accumulation order to differ across architectures — unlike `hausdorff`
and `self_intersections_per_1k`, which carry stated tolerances for exactly that reason. A field is
not supposed to move at all; if one of these shifts, either a field's definition changed or a
constant it depends on did, and both should require an acknowledgement.

`Unbounded` fields write **empty cells** for `declared_bound` and `bound_gap` rather than zeros.
`regress.sh` skips an empty value, and a literal `0` would read as *"the gradient is nowhere near the
bound"* — the opposite of *"there is no bound"*.

### M-269 — a grid-aligned ray double-counts shared edges, and the constructor shootout is what found it (T-018)

**M, and it is a defect in S-007 that S-007's own tests missed.** T-018 put the pseudonormal and the
winding number side by side on the same mesh for the first time. Their **magnitudes are computed by
the same code**, so any difference is a sign difference — and on `torus` the winding number read a
worst error of `1.223` against the pseudonormal's `0.014`, and `1.168` against `0.004` at 65³.

The cause: `intersect_x_ray` tested `u ∈ [0,1]` and `u + v ∈ [0,1]` **inclusively on both sides**,
while its own doc comment claimed a half-open convention. A ray passing exactly through an edge
shared by two triangles is therefore counted **twice**, which flips the parity and so the sign.

**And a grid-aligned ray hits shared edges constantly, not rarely.** The samples lie on grid lines,
and Marching Cubes places every vertex on a grid *edge*, so a `+x` ray from a sample runs along a
line that triangle edges also lie on. This is not a floating-point near-miss; it is exact
coincidence, by construction.

A half-open barycentric rule does not fix it in 3D: two triangles meeting at an edge do not agree on
which of their own `u`, `v`, `w` that edge is, so there is no local rule that makes exactly one of
them own it. The fix is Plantinga & Vegter's device, and they state the licence — *"We can take ε
arbitrarily small using a symbolic perturbation"* — so the ray **line** is nudged off the lattice by
`~1e-5` of a cell, with different fractions on `y` and `z` so it cannot land back on a cell diagonal.

**The perturbation has to move the query point too, and getting that wrong is a second bug the tests
caught immediately.** The first attempt nudged only the line and left `q` and the cone apex on the
sample: `χ` was then counted along one line while the correction was subtracted for a cone whose apex
sat on another, which is not the construction. It turned 0 misclassifications into **12** on the
smallest hole. With the query point moved with the line, `torus` matches the pseudonormal to the last
digit at both resolutions and the hole sweep is unchanged.

**The method note is the point.** S-007 shipped with four passing tests including a closed-sphere
calibration that read exactly `1.000000000` and `0.000000000`. A sphere is convex, so a grid-aligned
ray leaves through one triangle and the double-count has nothing to double. **It took a second
constructor computing the same magnitudes by different means to expose it** — which is what a
head-to-head harness is for, and is the argument for T-018 existing at all.

### M-270 — a benchmark that hands a repair algorithm perfect input measures the benchmark (T-018)

**M.** The first version of this shootout ranked all four field-based constructors in one table and
put **`band` first on both speed and accuracy** — 5.02 ms against the exact transform's 5.11, and a
worst error of `0.0729` against `0.1206`, `0.1623` and `0.1453`.

That was an artefact and the algorithm deserved no credit for it. S-004 **keeps the input** outside
its band by design; the harness fed it the analytic field, so most of its output was ground truth it
had merely copied. The other three read only signs and sub-cell crossings and discard magnitudes, so
they gained nothing from the same input — which is exactly why the comparison was unfair without
looking unfair.

Refed with the analytic field **scaled by two** — same zero set, twice the gradient, the shape a CSG
scale or a careless composition produces — the honest numbers are:

| | worst | worst near surface |
|---|---|---|
| degraded input | 2.26132 | 0.18740 |
| after reinitialisation | 2.26132 | **0.11325** |

It repairs the near band and **leaves the far field exactly as it found it**. That is what a narrow
band *is*, and stating it as a number rather than as a design intent is the difference between a
measurement and a description.

**The taxonomy is the fix, not a tolerance.** Three families now: *signs → field* (exact, swept,
marched), *field → field* (band, on degraded input), *mesh → field* (pseudonormal, winding). A
constructor may only be ranked against ones answering the same question.

### M-271 — the *exact* transform is the hungriest constructor, and the mesh-based pair are free (T-019)

**M.** Peak resident set per constructor, measured out of process — one process per constructor via
`--only`, reading its own `VmHWM` from `/proc/self/status`, minus a `baseline` run that builds the
input and stops. `sphere`, KiB above baseline:

| constructor | 33³ | 65³ | growth |
|---|---|---|---|
| `swept` | 376 | **4,228** | 11.2× |
| `marched` | 404 | 4,616 | 11.4× |
| `exact` | 608 | 6,212 | 10.2× |
| `band` | 784 | 6,704 | 8.6× |
| `pseudonormal` | 292 | **0** above the mesh | — |
| `winding` | 208 | **0** above the mesh | — |

Eight times the samples costs 8.6–11.4× the memory, so every one is essentially linear in sample
count with different constants.

**The ordering is the finding.** `signed_distance_field` is the *exact* transform and the only `O(n)`
one, and it is the most memory-hungry of the three — **47% more than sweeping**. Felzenszwalb &
Huttenlocher's separable algorithm needs per-row envelope state and a squared-distance staging array
on top of its output, and none of that is visible from the complexity class. "Exact and linear" does
not mean lean, and a chunked consumer choosing a constructor by asymptotics would have chosen wrong.

**`band` is the most expensive of all**, which is the opposite of what a *narrow-band* method sounds
like. It runs the bounded march and additionally holds the input it is repairing and an `on_surface`
flag array — M-256 measured it finalising only 13.4% of samples, and that is a saving in *work*, not
in *footprint*.

**The two mesh-based constructors cost nothing above the mesh they were handed** at 65³. Meshing
allocates more transiently than either of them does, so their own peak sits underneath it. That is
the honest reading of a zero here, not a measurement failure — and it means a consumer already
holding a mesh pays only for the output buffer.

**On the instrument.** T-018 declared this gap because a counting `GlobalAlloc` needs `unsafe impl`
and the workspace sets `unsafe_code = "forbid"`. Reading `/proc/self/status` is a plain file read and
needed no such thing; what it needed was leaving the *process*, since `VmHWM` is a high-water mark
over a whole process life and would otherwise attribute every constructor's peak to whichever ran
first. So the rule cost a process boundary, not a measurement. Linux only — macOS publishes no
`/proc` and its `ru_maxrss` needs `libc` — and CI runs Linux, which is where the figure is checked.

### M-272 — the pre-registration gate is a `const` assertion, and that is the whole of it (R-000)

**V.** R-000 asked for *"the feedback loop is currently a discipline; make it a compile error"*, with
the acceptance *"an experiment without a pre-registered `P-` fails to build"*.

It needs no proc macro and no build script. `experiment!("P-n")` expands to

```rust
const _CHECK: () = assert!(isomesh::experiment::is_preregistered("P-n"), "…");
```

and `is_preregistered` is a `const fn` doing a byte-wise `str` comparison, because `str`'s
`PartialEq` is not const and neither is `<[u8]>::eq`. A `const` item is evaluated whether or not its
value is read, so an unregistered id cannot reach a run. **The acceptance is proved by a
`compile_fail` doctest**, which is the only kind of test that can assert a compile error and is
therefore the only honest way to check this ticket.

`Preregistration` is `#[non_exhaustive]`, so a consumer crate cannot construct one. The harness takes
`&'static Preregistration` and there is no other way to obtain one, which makes the macro the sole
entrance rather than the polite one.

**Two smaller things the design gets for free.** `falsified_by` is a required field, so a hypothesis
with no stated refutation is unrepresentable rather than discouraged. And `records` names the columns
in advance, with `Run::record` panicking on any other set — so a metric that was predicted and then
quietly not measured is a failure, not a silence. That is F-002's one-sidedness in another costume:
the instrument has to be able to report the bad news.

The one thing it cannot enforce is that the registration was written *first*. Nothing in a source
file can. What it does is make registering a **commit**, so git carries the ordering, and
`scripts/backlog_gate.sh` fails if a registered id never reaches `FINDINGS.md`.

### P-14 — registered for R-003, before it was measured

**The residual non-manifold output under the dual extractors is the one-vertex-per-cell rule meeting
a cell with more than one surface component.** Almost all of it should sit in cells where Manifold
Dual Contouring emits more than one vertex, and MDC's own count should be strictly lower on every
field where either is non-zero.

**Falsified by** a substantial share of non-manifold vertices sitting in cells where MDC emits
exactly one vertex — which would mean the defect is not the one-vertex-per-cell rule and is somewhere
else entirely.

Records: `non_manifold_vertices`, `non_manifold_edges`, `multi_vertex_cells`,
`nm_vertices_in_single_vertex_cells`, `worst_link_components`. **Answered above — falsified on its
first clause.**

### M-275 / P-14 — FALSIFIED on its first clause, and the residue is an *edge* defect (R-003)

**M.** One 49³ grid per field, one extraction per dual rule, no chunking and no weld anywhere — so a
seam cannot be a second explanation. `docs/experiments/p-14.csv`.

| field | extractor | nm vertices | nm edges | cells MDC split | nm vertices in **un**split cells | link components | worst link degree |
|---|---|---|---|---|---|---|---|
| `gyroid` | `surface_nets` | 99 | 48 | 54 | **45** | 2 | 4 |
| `gyroid` | `manifold_dual_contouring` | **0** | 0 | 54 | 0 | — | — |
| `fbm_terrain` | `surface_nets` | 38 | 19 | 19 | **19** | 1 | 4 |
| `fbm_terrain` | `manifold_dual_contouring` | **0** | 0 | 19 | 0 | — | — |
| `noise_cavity` | `surface_nets` | 597 | 314 | 290 | **307** | 2 | 4 |
| `noise_cavity` | `dual_contouring` | 597 | 314 | 290 | 307 | 2 | 4 |
| `noise_cavity` | `manifold_dual_contouring` | 106 | 53 | 290 | **106** | 1 | 4 |

`sphere`, `torus`, `box_exact`, `csg_difference` and `thin_plate` are clean under all three.

**Clause two held.** Manifold Dual Contouring's count is strictly lower on every field where either
is non-zero — 99 → 0, 38 → 0, 597 → 106. Splitting multi-component cells removes the defect entirely
on two of three fields.

**Clause one is false, by about half.** *"Almost all of them sit in cells where MDC emits more than
one vertex"* — measured, **45 of 99, 19 of 38 and 307 of 597 sit in cells MDC did not split**, and
MDC's own residual 106 are **106 of 106** in unsplit cells. So the one-vertex-per-cell rule explains
roughly half the defect on the dual rules and none of what survives its own remedy.

**What the other half is, and it reframes the ticket.** The vertex-to-edge ratio is 597/314, 99/48,
38/19 and 106/53 — **exactly two to one, on every row**. Every non-manifold vertex is an *endpoint of
a non-manifold edge*. So there are not 597 defects plus 314 defects; there are 314 defects, each
counted twice more by the vertex walk. The residue is an **edge with three or more incident faces**,
which is a different failure from a vertex whose link falls apart, and no amount of cell splitting
addresses it — MDC splits the cells and the edges remain.

**A second shape of vertex defect, which the first run of this experiment could not see.**
`worst_link_components = 1` on a vertex `validate` flagged looks like a contradiction and is not: a
link can be **connected and still not a simple cycle**, with one link vertex reached by four link
edges instead of two. That is a **pinch**, not a bowtie. The added `worst_link_vertex_degree` column
reads **4** on every non-clean row, and separates them: `fbm_terrain` and MDC's residue on
`noise_cavity` are *purely* pinches (components 1, degree 4); `gyroid` and `noise_cavity` under
Surface Nets and Dual Contouring carry both (components 2 **and** degree 4).

**The instrument was wrong before the hypothesis was.** `link_components` alone would have reported
`fbm_terrain` as having non-manifold vertices with one-component links and left that unexplained.
Adding a column rather than editing the registration is what M-273's relaxation made possible.

### M-273 — the first thing done with the pre-registration mechanism was amend a registration to fit the code (R-002)

**✗ against my own practice, caught by the mechanism itself.** R-000 shipped a harness whose
`Run::record` demanded the row's keys be **exactly** the registered `records`. The first experiment
needed a `field` column to identify its rows — and the path of least resistance was to edit `P-8`'s
registration and add it. Which is amending a prediction to fit the code, one commit after building
the thing whose entire purpose is to stop that.

It was caught on the next experiment, when P-9 panicked with *"`field` is not one of the records this
experiment registered"* and the same edit beckoned again.

**The design was wrong, not only the act.** `records` is a list of **metrics that must be reported**,
not a schema for the file. Row keys are not anybody's hypothesis. `Run::record` now panics only on a
**missing** metric and permits extra columns, which are written after the registered ones, sorted —
so identifying a row never requires touching a prediction. Both amendments were reverted to what
R-000 committed.

**The general rule, which is Part 5 material.** A gate that is stricter than the thing it protects
creates pressure to weaken the gate, and the weakening is indistinguishable from the abuse it was
built to stop. Check exactly the property you mean: *"every registered metric appears"* was the
property; *"the columns are exactly these"* was a stronger claim nobody needed and could not keep.

### M-274 — the fixture never contained the configuration both experiments were about (R-002)

**✗.** P-8's first run used a 2×2×2 block of **8-cell** chunks at `h = 4/35` anchored at `-2.0`. That
spans `16 × 4/35 = 1.83` units of a 4-unit domain, so the block covered a **corner** of every
reference field and never reached the origin. P-8 recorded numbers, they were internally consistent,
and the write-up asserted a mechanism about `k`-way merges.

P-9 measured the same fixture and reported **`buckets_of_three_or_more = 0` in all 49 rows it produced**. There
was no `k`-way merge anywhere in it. The mechanism was a description of something the run did not
contain, arrived at by reasoning from the algorithm rather than from the data — which is the failure
mode this file exists to catch, committed while writing an entry for this file.

Centring the block and taking 18-cell chunks — `36 × 4/35 = 4.11`, covering `[-2, 2]` — puts the
eight-chunk corner at the origin and the four-chunk edges through the surface. `k ≥ 3` buckets appear
immediately: **6–32 per configuration** on every primal extractor. Both experiments were re-run and
P-8's entry replaced; the conclusion survived and got much stronger (the gate's damage went from
"up to 42 non-manifold vertices" to "up to 791"), which is luck rather than vindication.

**The rule.** A fixture that produces zero of the thing being studied looks exactly like one that
produces some, because every other column is populated. **Count the phenomenon, not just the
outcome** — P-9's `buckets_of_three_or_more` column is what exposed this, and it existed only because
the pre-registration demanded it.

### M-276 — the dual methods' non-manifold edges are the ambiguous face, all 314 of them (A-021)

**M + V.** `noise_cavity` at 49³ under Surface Nets, no chunking and no weld. For each non-manifold
edge, the grid face its two cells share, and how many of that face's four boundary edges carry a sign
change:

| population | crossed boundary edges | count |
|---|---|---|
| **non-manifold** edges | **4** | **314** |
| manifold (2-face) edges | 2 | 30,891 |
| manifold (2-face) edges | 0 | 8 |

**Every one. No exceptions, and no overlap between the populations.** Zero non-manifold edges join
cells that are only *diagonally* adjacent, so the quad's triangulation diagonal is never the culprit
either.

**The mechanism, derived from the code and then confirmed by a control.** `DualMesher::emit_quads`
winds a quad as `(0,1,2)` and `(0,2,3)`, so each quad **side** lands in exactly one triangle and only
the diagonal `0–2` lands in two. The mesh edge between two face-adjacent cells is a side of the quad
of **every crossed boundary edge of their shared grid face**. So the face count on that mesh edge
*equals* the number of crossed boundary edges: two is manifold, four is not.

A square whose four boundary edges are all crossed is one whose corner signs **alternate** — and that
is the literature's *ambiguous face*, verbatim (`10.1145/195826.195828`):

> *"Each cell face has four edges, an even number of which will contain intersection points. If two
> edges of a face contain intersection points, then the only choice is to connect them with an edge…
> However, if a cell face has intersection points in all four of its edges, there are choices of how
> to connect up pairs… which justifies the name ambiguous face."*

**So this is the A-002 series arriving from the dual side.** Grosso's trilinear work, the asymptotic
decider, `FaceAmbiguity::AsymptoticDecider` — all of it exists to make *one* choice on an ambiguous
face in the Marching Cubes path. The dual path makes **both choices at once**: it emits a quad for
each of the four crossed edges, and all four pass through the same mesh edge.

**It is therefore not inherent to "one quad per crossed edge".** The ambiguous face has a choice by
definition, and a decider that picked one pairing would emit two quads there instead of four. That is
a specific remedy with a precedent inside this crate, and it is **A-022** rather than a claim made
here — the dual path has no face-ambiguity setting today, and asserting what MDC's own criterion is
without reading it would be rule 5's failure.

**It also closes P-14's loose end exactly.** P-14 found MDC's residual 106 non-manifold vertices
sitting **106 of 106** in cells MDC did *not* split. Now explained: MDC splits cells with more than
one surface component, and an ambiguous *face* between two cells that each hold a single component
does not meet that criterion. The two measurements interlock — one says the residue is not about
cells, the other says what it is about.

**Method note: the reasoning was wrong before it was right, and the probe is what caught it.** Reading
`emit_quads` first produced "a mesh edge gets a quad per crossed boundary edge, and a face with a sign
change has two, so every dual mesh is non-manifold everywhere" — which contradicts every measurement
in this file. The error was treating a quad as contributing two faces to each of its sides. Printing
the face-count histogram for a **plain half-space** (`{1: 8, 2: 8}` — no edge above two) settled it in
one run and cost less than the third attempt at the algebra.

### M-277 — the index generator was blind to the shape its own file had moved to (R-004)

**M.** `scripts/findings_index.sh` matched exactly two things: table rows `| M-n | … |` in Part 2, and
`### ✗n` headings in Part 1. Measurements stopped being table rows at **M-255** — the entries had
outgrown a cell — and the generator never noticed. Counted directly:

| | count |
|---|---|
| `^\| M-n \|` table rows | 249 |
| `^### M-n —` heading sections | 22 |
| what the index said | **"319 entries — 249 measured"** |

So every measurement from S-004 onward, **twenty-two of them including all of R-001, R-002, R-003 and
A-021**, was in `FINDINGS.md` and absent from the index a reader is told to consult first.

**`--check` was green throughout, and correctly so.** A staleness check compares the artefact against
what the generator produces, so it cannot see the generator's own blind spot: both agreed, and the
one gate guarding the index was structurally unable to report this. The script's own header says an
index that is allowed to drift *"answers confidently and wrongly"*; it was doing that about itself.

Found while adding R-004's entry, by going to copy the row format and finding that the last
twenty-two entries do not use one.

Fixed by also matching `### (M|V|O)-n — …` headings, with `M-275 / P-14`'s two-id form resolved to
the first id and a table row winning over a heading for the same id. **`P-` headings stay unindexed
deliberately** — `### P-8 … P-13` names six predictions in one heading and would collide with each of
their own outcome sections. Index after the fix: **341 entries, 271 measured.**

### M-278 / P-11 — HELD, and the crack budget is not where the ticket assumed (R-004)

**M.** `benches/experiment_p11.rs`, `docs/experiments/p-11.csv`, 40 rows: `sphere` and `torus` ×
five level-0 spacings × two LOD pairs × two arithmetics. A fine block meets a half-resolution one
across a seam, stitched with Lengyel transition cells at zero width, and every block starts twelve
cells from the grid origin so that both arithmetics can differ at all.

| arm | sample position |
|---|---|
| `canonical` | `o + h·(base + local)` — one multiply from a global integer index |
| `offset` | `(o + h·base) + h·local` — a step added to the block's own corner |

**Both come from the shipped functions.** `TransitionCell::sample` already takes `(origin, base)`, so
handing it the face's own world origin with a zero base *is* the pre-M-73 code. Marching Cubes has no
integer base, so its canonical arm roots the extraction at the grid origin and clips to the block's
cells; `clip_agrees_with_the_block` checks that clip at every power-of-two spacing, where the two
arms compute identical positions and a difference could only be the clip.

#### The result, in one table — seam-plane boundary edges

| pipeline | `canonical` | `offset` |
|---|---|---|
| **weld + transition cells** (what the crate ships) | **0**, all 20 rows | **0**, all 20 rows |
| **bit-identity + transition cells** | **0**, all 20 rows | **0** at `h = 0.125` and `0.0625`; **63–348** at `0.1`, `1/12`, `1/14` |
| **weld, transition cells removed** | 32–184 | 32–184, differing from `canonical` in 5 of 20 |

**P-11 HELD, in both clauses.** Canonical reconstruction gives **zero** seam cracks at every cell
size — not only the powers of two — under *both* merge policies, and every one of the fine block's
seam vertices is bit-identical to its partner (`shared == pairs`, 20 of 20 rows, up to 124/124). The
registered falsifier — *"cracks surviving canonical reconstruction"* — did not fire anywhere, so the
defect does not localise back into Transvoxel.

#### What the ticket assumed, and what is actually true

R-004 is titled *"quantify the crack budget: arithmetic vs algorithm"*, and the split is sharper than
expected but on a different axis:

- **The algorithm owns the whole of the *visible* budget.** Remove the transition cells and the seam
  opens to 32–184 boundary edges with a widest hole of **1.03–3.01 cells** — a hole you can see the
  sky through. Put them back and it is 0. Arithmetic changes none of that.
- **The arithmetic owns the whole of the *invisible* one.** Its worst positional disagreement is
  `1.440e-15` world units, and the weld epsilon is `h · 1e-4`. So under the shipped pipeline the
  offset arm is indistinguishable from the canonical one, and under bit-identity it is 63–348 cracks
  wide. That is ✗18.

#### Three things the sweep found that were not predicted

**1. The seam *plane* is the wrong unit of analysis, which sharpens M-32 rather than falsifying it.**
M-32 is confirmed where it speaks: at `h = 0.125` and `0.0625` the two arms are bit-identical
everywhere and every crack column is 0, so *"recommend power-of-two cell sizes for chunked worlds"*
holds. But at `h = 0.1` and `1/12` the seam plane's own coordinate agreed **bit for bit**
(`seam_plane_delta = 0`) while only **24 of 92** vertices did. The disagreement is not confined to
the coordinate along the seam normal; it is *every* coordinate that two blocks with different bases
reach by different groupings — the fine block's `(o + h·12) + h·i`, the patch's
`(o + h·(12 + 2j)) + h·s`, the coarse block's `(o + 2h·6) + 2h·j`. Only `h = 1/14` disagreed on the
plane itself, at `4.441e-16`.

**2. Bit-exact sharing degrades continuously, and one spacing loses it entirely.** Partners that are
bit-identical, offset arm: `0.125` and `0.0625` → **100%**; `0.1` → 26/76, 14/48, 6/36, 4/24;
`1/12` → 24/92, 13/56, 12/44, 9/24; `1/14` → **0 of 108, 0 of 72, 0 of 52, 0 of 40**. There is no
threshold here to sit above — a spacing either happens to keep the grouping exact or it does not.

**3. An ulp of coordinate can be a cell of geometry.** The bit-identity crack is a hairline
(`~1e-14` cells) in 8 of the 12 non-power-of-two offset rows, `8.1e-3` and `1.7e-2` cells in two
more, and **1.053 and 2.076 cells** on `torus` at `h = 1/12`. A sample perturbed in its last bit
crossed zero, the two sides disagreed about whether an edge was cut, and a patch that should have
been emitted was not. The same effect shows in the control: with transition cells removed, the offset
arm's crack count differs from the canonical one in **5 of 20 rows** (68→81, 44→48, 37→41, 152→160,
72→80). So "the arithmetic only moves things by an ulp" is true on average and false in the tail.

#### The consequence for the crate, which is a ticket rather than a claim

`ChunkLayout::world_of_sample`'s doc calls itself *"the single place a sample's world position is
defined. Everything else routes through it"* — and **no extractor does.** `Extractor::extract_into`
takes `origin: [R; 3]` and computes `origin + h·local`, so a chunk at a non-zero base gets the offset
arithmetic by construction and there is no argument a caller can pass to avoid it. That is why this
experiment had to reach the canonical arm by clipping an over-extended extraction. **X-005** owns it;
it is an API change to the crate's central trait and therefore a decision, not a fix.

### P-15 — registered for R-007, before it was measured

**The residue M-279 leaves, stated so it can fail.** Surface Nets executes 207 instructions per
sample at IPC 1.22 where Marching Cubes executes 132 at 4.04 — 1.57× the instructions and 5.24× the
cycles — and none of the obvious candidates survived: not the crossed-edge gather (a field with no
surface costs the same), not branches (they fall), not allocation (zero page faults), not the TLB,
and not, at 16.7 M samples, the misses (a 2.4× swing moves cycles by 0.4%).

> **H.** More than half of the dual mesher's cycles per sample are spent in `emit_quads`, which is
> three unconditional `O(n³)` sweeps over the sample grid rather than work proportional to the
> surface.

The reasoning, so the prediction is a claim and not a shrug: `place_vertices` reads the same eight
corners per cell that Marching Cubes' march does, and `sample` is identical, so the dual's extra ~75
instructions per sample have to come from `emit_quads` — which walks every grid edge on all three
axes and loads two samples *before* the sign test that would let it skip. Two of its three passes
have their innermost loop on the wrong axis (`v = (axis + 2) % 3`), giving strides of `nx` and
`nx·ny`. Whether that is where the *cycles* go rather than merely the instructions is exactly what is
not known.

**Falsified by** `emit_quads` accounting for half or less, which puts the cost in `sample` or
`place_vertices` — work Marching Cubes does too, at four times the IPC — and means the dual's IPC is
lost to something other than its extra traversal.

**Records** `stage`, `cycles_per_sample`, `instructions_per_sample`, `ipc`, `samples`.

**A note on the instrument, because the obvious one is unavailable.** `STALLED_CYCLES_BACKEND` would
answer this directly and AMD does not map it (`perf_event_open` → ENOENT). An intra-extractor
breakdown therefore needs either an ablation seam in `DualMesher` or per-stage counter windows inside
it, and which of those is acceptable under the one-path rule is part of R-007.

### M-279 / P-12 — the mechanism is FALSIFIED, and the registered falsifier could not have caught it (R-005)

**M.** `benches/experiment_p12.rs`, `docs/experiments/p-12.csv`, 90 rows on the Ryzen 9 5900X.
Six hardware events and one software event through `perf_event_open`; every counter ran the whole
window (`counter_time_ratio` 1.0 on all 90 rows) and the clock held at **4.148–4.188 GHz**, which the
rows carry rather than imply (M-280).

#### What P-12 predicted, and what the numbers say

`sphere`, `f32`, the resolution sweep's own grid:

| n | MC cyc/sample | MC instr | MC IPC | MC miss | SN cyc/sample | SN instr | SN IPC | SN miss |
|---|---|---|---|---|---|---|---|---|
| 48 | 33.7 | 138.8 | 4.12 | 0.104 | 144.0 | 210.5 | 1.46 | 1.635 |
| 128 | 32.4 | 133.2 | 4.11 | 0.022 | **178.4** | 207.7 | **1.16** | **5.465** |
| 256 | 32.5 | 131.5 | 4.04 | 0.020 | 170.2 | 206.8 | 1.22 | 4.258 |

Both of P-12's *observable* clauses came out its way. Marching Cubes' miss rate stays flat
(0.104 → 0.020, falling); Surface Nets' rises 2.6× (1.635 → 4.258), though not monotonically — 112³
and 144³ sit *below* 64³.

#### The mechanism clause is false, and one control settles it

P-12 says the cost **is** the four-cells-around-a-crossed-edge gather. That gather runs once per
crossed edge, which is `O(n²)`; the cost is `O(n³)`. So the sweep is run a second time on a field
with **no surface at all** — a sphere of radius 10 over a domain of half-extent 2, every sample
negative, no edge crossed, no vertex placed, no quad walked — while the dense per-cell state is
allocated, filled and scanned exactly as before.

At 256³, Surface Nets: **168.6 cycles/sample and 4.277 misses/sample with 0 triangles**, against
**170.2 and 4.258 with 153,552**. A difference of 0.9%. Removing every gather the hypothesis blames
changes nothing, so **the hypothesis is false**. Dual Contouring, which shares the engine, lands on
the same number (167.4 cycles, 4.240 misses, 0 triangles).

#### The other two candidates the registration named are excluded too

**Branch misprediction cannot produce a rising cost because it falls**: Surface Nets 0.0436 → 0.0267
misses per sample across the sweep, Marching Cubes 0.0343 → 0.0115. **Allocation is not it either**:
page faults after the two warmup runs are **0 on all 90 rows**. And a candidate the registration did
not name, added because "the working set grew" and "the page walk grew" have the same symptom:
**dTLB read misses peak at 0.104 per sample and are typically ~1e-3**, transparent huge pages being
`always` on this host.

#### What it actually is: IPC, on an instruction stream that does not grow

Surface Nets executes **1.57×** Marching Cubes' instructions and takes **5.24×** the cycles. Its
instruction count per sample **falls** across the sweep (210.5 → 206.8) while its cycles **rise**
(144.0 → 170.2), so the whole of the superlinearity is a **16% decline in IPC** (1.46 → 1.22) on a
constant instruction stream. Marching Cubes' IPC is flat at 4.04–4.28.

#### And the miss column does not explain the cycles either

Three grids of **16.7 M samples**, no surface, same code, differing only in axis order:

| shape | SN misses/sample | SN cycles/sample |
|---|---|---|
| `68×496×496` | 3.274 | 151.49 |
| `496×68×496` | **1.362** | 151.74 |
| `496×496×68` | 3.360 | 151.14 |

**A 2.4× spread in misses buys 0.4% of cycles.** At this size the misses are streaming and the
prefetcher hides them, so the LLC-miss column tracks the dual's cost without accounting for it.

**The one place misses do cost is the 128³ spike, and it is a stride effect.** 127³ / 128³ / 129³ on
`sphere`: **152.5 / 178.4 / 152.1** cycles and **2.468 / 5.465 / 2.413** misses. Their working sets
differ by 2%; only 128 has a plane stride of exactly 64 KiB. The same spike is on the empty field
(149.1 / 176.4 / 148.8 cycles), so it is not the surface. These misses are conflict misses and cost
about 10 cycles each; the streaming ones above cost nothing. That is the clue M-45 recorded on two
machines and nobody had followed.

#### The registered falsifier was insufficient, which is the method finding

`falsified_by` was *"flat miss rates, pointing at branch misprediction or allocation instead."* Miss
rates were **not** flat, so **by its own stated falsifier P-12 survives** — and its mechanism is
still false. The falsifier named an observation the hypothesis could fail, but not one that could
separate it from its rival: "the gather misses" and "the `O(n³)` scan misses" both predict rising
misses on a growing grid. The control that settled it — a field with no surface — was not in the
registration and was added while writing the harness.

#### O-11 is narrowed rather than closed

It asked for *"a profile or cache-miss counter at 192³ vs 256³"* and hypothesised the stride-`n²`
gather. The counters are in and the gather is out. What remains is **where the dual's IPC goes** —
not the gather, not branches, not allocation, not the TLB, and not, at this size, the misses.
`STALLED_CYCLES_BACKEND` is the event that would name it and **AMD does not map it**;
`perf_event_open` answers ENOENT. **R-007** owns the residue.

### M-280 — this harness's nanoseconds are not a unit, and the committed Zen 3 sweep is 1.45× stale (R-005)

**M.** Two runs of the same binary at the same commit reported Marching Cubes at 48³ as **8.13 and
14.66 ns/sample** while cycles/sample held at ~34 both times. This host runs `amd-pstate-epp` on the
`powersave` governor over **1.96–5.62 GHz**, so that is a 1.8× clock excursion and nothing about the
code. `experiment_p12` now leads with `cyc/samp` and carries a `ghz` column — cycles ÷ nanoseconds —
so a row states the clock it was taken at instead of leaving it to be inferred from a number that
looks like a measurement.

Two full 90-row runs after that change: **ns median |Δ| 0.43%, cycles 0.31%**, clock 4.15–4.20 GHz
throughout. The nanosecond is fine *when the clock is*; the column is what says which.

**Which then found something bigger.** Re-running `benches/resolution_sweep` at the current commit,
with the clock sampled every five seconds and steady at **4.20 GHz**, five times:

| | committed `resolution_sweep-ryzen9-5900x.csv` | measured now | |
|---|---|---|---|
| Marching Cubes 256³ | 221.363 ms | **152.2–153.3 ms** | **1.45× faster** |
| Surface Nets 256³ | 823.442 ms | **693–721 ms** | 1.18× faster |

**Consequences for what is written down.** ✗14's Zen 3 ratio was `SN/MC = 3.72×` at 256³; measured in
one binary and one run it is **5.43×** (M-282), so ✗14's conclusion is strengthened and its number is
wrong. M-45's *"the M5 is 2.76× faster than the Ryzen on Marching Cubes at 256³ (80.2 vs 221.4 ms)"*
**cannot be quoted at all**: its Ryzen half is superseded by a real code change and its M5 half was
not re-run, so the figure is unmeasured rather than merely stale.

**One excursion in six runs, unexplained.** One of the six gave **264.714 ms** for Marching Cubes at
256³ at the same steady 4.20 GHz and the same binary — 1.74× the other five, which cluster within
0.7%. Surface Nets did not move in that run (697.8 against 693–721). No cause identified; recorded so
that a single ms figure from this host is not read as a measurement.

**Corrected the same night, and the correction matters (M-281).** The table above attributes the gap
to the code, and that attribution was not earned when it was written: the "measured now" column came
from the *current* `resolution_sweep` binary, and adding one unrelated function to that binary moves
its own Marching Cubes 256³ row from 152.5 to **130.8 ms**. So 152.5 is a property of a build, not of
the algorithm.

**Settled by re-running the old commit rather than by arguing.** A worktree at `d2ab82a` — the commit
M-45 cites — built and run on this machine tonight reproduces its committed CSV to within 1.8%:
Marching Cubes 256³ **222.649 ms** against 221.363, Surface Nets **813.964** against 823.442, Dual
Contouring **861.370** against 877.359. **The committed Zen 3 figures are sound and reproducible**, so
the difference *is* a real code change and not the clock and not the machine — but its size is
`222.6 → 127.8`, **1.74×**, measured against M-282's family run rather than against a binary whose
layout costs it 17%. Read the direction from this entry and the magnitude from M-282.

**The CSV was restored rather than overwritten.** ✗14, M-19, M-20, M-21, M-22 and O-11 all quote
exact figures from it, and the archive says three separate times that re-measuring the family belongs
to **M-001**. M-001 is referenced **nineteen times** across `FINDINGS.md` and `BACKLOG_ARCHIVE.md`
and **had no ticket row in either file** — so the designated home for this re-measurement did not
exist and it could never have been taken off the queue. Now filed.

### M-281 — a timing here is a property of the binary, not only of the code (M-001)

**M.** `benches/family` and `benches/resolution_sweep` measure Marching Cubes on the same field, the
same resolutions, the same warmup and median rule, on this machine with the clock held at 4.18 GHz —
and disagreed by a **uniform 1.24–1.36×** at every resolution:

| n | `resolution_sweep` ms | `family` ms | ratio |
|---|---|---|---|
| 16 | 0.047 | 0.038 | 1.24 |
| 64 | 2.681 | 2.011 | 1.33 |
| 128 | 20.972 | 15.461 | 1.36 |
| 256 | 172.477 | 127.605 | 1.35 |

A uniform ratio **at 16³**, where the whole run is 40 µs and nothing is under memory pressure, rules
out cache, TLB and allocation at a stroke.

**What it is not.** `benches/layout_bias` puts both loop shapes in **one binary** — the sweep's single
drained loop with `black_box(triangle_count())` against `family`'s split loops with
`black_box(&mesh)` — in both orders, so M-197's *"whichever runs second pays"* cannot apply. They come
out **identical, ratio 0.991–1.002**, and the file asserts that rather than printing it. Nor is it the
counters (`family` without them: 127.6 ms against 127.1 with), nor extractor lifetime (rebuilt per
resolution: 127.1 against 131.6 carried).

**What it is.** Adding **one unrelated function** to `resolution_sweep.rs` — a copy of the probe,
called once before the sweep — moved that binary's own Marching Cubes 256³ row from **152.5 ms to
130.8 ms**. Same measured code, same loop, same machine, same clock; 17% from linking something else
in beside it. This is layout bias, and it has a paper: Mytkowicz, Diwan, Hauswirth & Sweeney,
*Producing Wrong Data Without Doing Anything Obviously Wrong!*, ASPLOS 2009
(`10.1145/1508284.1508275`), whose finding is that link order and environment size move measurements
by enough to invert a conclusion.

**The consequence for this repository, which is the point.** A millisecond figure is comparable
**only against other figures from the same binary and the same build**. Ratios measured side by side
in one run survive; absolute numbers quoted across benches do not. That is precisely M-001's design —
one bench, one process, one run — and it is now measured rather than asserted.

### M-282 — the whole family, in one binary and one run (M-001)

**M.** `benches/family`, `docs/measurements/family.csv`. `sphere`, `f32`, the sweep's resolution set,
median of 5 after 2 warmups, clock 4.17–4.25 GHz on every row, Ryzen 9 5900X, single thread.

| algorithm | 256³ ms | × MC | ns/sample | cycles/sample | **IPC** | triangles |
|---|---|---|---|---|---|---|
| `marching_cubes` | **127.8** | 1.00 | 7.62 | 31.8 | **4.13** | 153,548 |
| `marching_cubes+decider` | 130.8 | 1.02 | 7.80 | 32.6 | **4.16** | 153,548 |
| `marching_tetrahedra` | 261.8 | 2.05 | 15.60 | 65.3 | **4.19** | 458,568 |
| `surface_nets` | 693.8 | 5.43 | 41.35 | 172.7 | **1.20** | 153,552 |
| `dual_contouring` | 751.1 | 5.88 | 44.77 | 187.0 | **1.35** | 153,552 |
| `manifold_dual_contouring` | 771.0 | 6.03 | 45.96 | 192.0 | **1.42** | 153,552 |
| `subgrid_marching_tetrahedra` | *budget* | — | 1575.0 | 6583.2 | **3.69** | 113,568 (at 128³) |

**The strongest pattern in the table is the IPC column, and it is a clean partition.** Everything
table-driven runs at **3.7–4.2** instructions per cycle; everything built on `DualMesher` runs at
**1.20–1.42**. Three algorithms with three different vertex rules — a centroid, a QEF, a QEF with
cell splitting — land within 18% of each other and a factor of three below the rest, which says the
cost is the shared scaffolding rather than any of the rules. That is R-007's target, arrived at from
a second direction.

**Four things nobody had priced.**

- **The asymptotic decider costs 2.4%** — 130.8 against 127.8 ms, same 153,548 triangles. ✗11 and the
  whole A-002 series argue about what it *fixes*; this is what it costs.
- **Marching Tetrahedra is 2.05× Marching Cubes in time and 2.99× in triangles** (458,568 against
  153,548). P-1 predicted the *3.0× ratio* and it lands at 2.99 — but the time ratio is only 2.05,
  because Marching Tetrahedra has the **highest IPC in the family**. Three times the triangles for
  twice the time.
- **Manifold Dual Contouring's guarantee costs 2.6% over Dual Contouring** — 771.0 against 751.1 ms.
  Given A-010 and M-59, this is the price of the entry that takes the manifold zero, and it is
  approximately free.
- **Subgrid Marching Tetrahedra is 100.7× classic Marching Tetrahedra and 215× Marching Cubes**, at
  128³ on Zen 3. M-98 measured 70× and 196× on the Apple M5 at 33³/65³, and said *"the constant is
  the whole story"* — so the constant is also machine-dependent, by 1.44× between these two hosts.
  It is the one entry the 2000 ms budget stops, at 128³, which the CSV records rather than omitting.

**✗14's ratio has widened and its numbers are superseded.** `surface_nets / marching_cubes` at 256³
is **5.43×** here against the 3.72× on record, and **3.19×** at 16³ against M-45's *"2.46× behind even
at 16³"*. Both are within-binary ratios, which M-281 says is the comparison that survives.

### M-283 / P-13 — ~~FALSIFIED~~, and the fixture that confirmed it agreed to four decimal places (R-006)

> **⚠ THE VERDICT IN THIS ENTRY IS WRONG. P-13 HELD — see M-289.** The reference gradient this
> experiment compared against was normalising a cancellation residue at points epsilon-outside the
> surface, so roughly half the vertices near a face were measured against a random unit vector. The
> corrected numbers are in M-289: 6,959 past-90° vertices become **472**, the worst angle **tracks**
> `(180° − θ)/2` and is bounded by it in 136 of 168 rows, and it is resolution-invariant to four
> significant figures. The entry is kept whole because what it got right — the `θ = 180°` control,
> the median of zero, the four-decimal fixture — it got right for the right reasons, and because a
> falsification that was itself falsified is the most useful kind of thing this file can hold.


**M.** `benches/experiment_p13.rs`, `docs/experiments/p-13.csv`, 384 rows: an exact convex wedge of
controllable dihedral × 4 resolutions × 2 apex alignments × 3 rotations × Marching Cubes and Dual
Contouring. The metric is M-66's: the angle between a vertex's **area-weighted face normal** and the
field's own gradient there, over vertices more than one cell from the domain boundary.

#### The prediction, and the fixture that appeared to confirm it exactly

Across a crease of dihedral `θ` the surface normal turns by `180° − θ`, so a vertex given the average
of the two faces should sit at most `(180° − θ)/2` from either. With the wedge's bisector on a grid
axis — the obvious way to build the fixture — that is what comes out, **to four decimal places**:

| θ | predicted | 17³ | 33³ | 65³ | 129³ |
|---|---|---|---|---|---|
| 30° | 75.0 | **75.0000** | **75.0000** | **75.0000** | **75.0000** |
| 60° | 60.0 | **60.0000** | **60.0000** | **60.0000** | **60.0000** |
| 90° | 45.0 | **45.0000** | **45.0000** | **45.0000** | 90.0000 |

Four significant figures, three resolutions, three dihedrals. It is wrong.

#### Rotate the wedge by 17° and the agreement is gone

Same dihedrals, same resolutions, apex off the lattice, wedge turned 17° about its own crease
(129³, and the p99 is there because one vertex should not set a conclusion):

| θ | predicted | MC worst | MC p99 | MC median | DC worst | DC p99 | DC median |
|---|---|---|---|---|---|---|---|
| 30° | 75.0 | 90.96 | 89.04 | 0.000 | 75.00 | 75.00 | 0.045 |
| 45° | 67.5 | 85.93 | 85.93 | 0.000 | 98.46 | 41.88 | 0.036 |
| 60° | 60.0 | 92.13 | 60.00 | 0.000 | 60.05 | 60.00 | 0.035 |
| 90° | 45.0 | 71.57 | 71.57 | 0.000 | 90.00 | 71.57 | 0.028 |
| 120° | 30.0 | 93.43 | 60.00 | 0.000 | 75.12 | 74.95 | 0.048 |
| 150° | 15.0 | 90.96 | 78.43 | 0.000 | 78.43 | 78.14 | 0.019 |
| **170°** | **5.0** | **87.88** | **87.88** | 0.000 | **91.91** | **91.42** | 0.396 |
| **180°** | **0.0** | **0.0000** | 0.0000 | 0.000 | 0.0028 | 0.0000 | 0.000 |

**A five-degree crease produces an eighty-eight-degree disagreement.** Across all 168 non-control
rows, Marching Cubes' worst runs 5.00°–128.00° and **149 of them are at or above 60°** whatever the
dihedral. The registered falsifier — *"the angle failing to track the dihedral prediction, which makes
it a defect with a location rather than a property"* — fires.

#### One clause of P-13 held, and it is the one the `θ = 180°` control settles

With the crease removed entirely the disagreement is **0.0000° worst and 0.0000° mean for Marching
Cubes at every resolution**, and 0.0028/0.0003 for Dual Contouring. Marching Cubes is *exact* on a
linear field — every crossing lands on the plane — so this is not a floor, it is zero. So the error
does require a sharp feature and does not require a coarse grid, which is P-13's first clause. What
is false is the second: the magnitude is **not** a function of the dihedral and is **not** predictable
from the field.

#### And it reconciles the two halves of M-66 that looked contradictory

M-66 reported a *mean* that falls with resolution and a *worst* that does not. The median here is
**0.000° in every row** — the disagreement is confined to the crease, which is a one-dimensional set
in a two-dimensional mesh. So refining the grid adds smooth vertices that are exactly right and
dilutes the mean, while leaving the crease's own error untouched. Both halves, one mechanism.

#### Dual Contouring does not fix it, and on one axis is worse

A-007 exists to recover sharp features, so the dual is the natural remedy. It is not one. At `θ = 90°`
its worst **rises** with resolution where Marching Cubes' does not, and so does its mean:

| | 17³ | 33³ | 65³ | 129³ |
|---|---|---|---|---|
| Marching Cubes worst / mean | 90.00 / 13.28 | 90.00 / 8.51 | 90.00 / 9.13 | 71.57 / **8.45** |
| Dual Contouring worst / mean | 36.17 / 2.52 | 44.96 / 3.33 | 71.57 / 5.06 | 90.00 / **6.79** |

The dual starts three times better and converges *toward* Marching Cubes from below. Its median stays
an order of magnitude worse than Marching Cubes' exact 0.000, which is the QEF placing vertices off
the plane in the smooth region — the trade ✗12 and E×1 already describe, now visible in the normals.

#### The past-90° vertices are real, and they are on the surface

6,959 vertices under Marching Cubes and 4,868 under Dual Contouring, over the non-control rows, have
an area-weighted normal **more than 90° from the field gradient** — pointing into the solid. Worst
128.00° and 125.30°.

That should be impossible: an area-weighted normal is a convex combination of incident face normals,
so it lies in the cone they span, which for two planes is `(180° − θ)/2 ≤ 75°` wide. The obvious
escape is that the vertex is not on the surface — Marching Cubes interpolates linearly and a wedge is
not linear near its apex. **Measured rather than assumed, and the escape is closed:** over the rows
that have one, the median `|f(v)| / h` at those vertices is **0.0000** and the largest is 0.28. They
are on the isosurface.

So the incident faces are not confined to the two planes' cone: **the crease is bridged by triangles
that face somewhere else**, and for an acute wedge a whole band along it is thinner than a cell, which
is M-15's *"any feature thinner than one cell forces two sheets through it"* arriving on a sharp
field. Whether that is inherent to one vertex per crossed edge or a defect with a fix is **R-008**.

**R-008 answered it and this reading is narrowed (M-288): the crease accounts for 20% of them.**
536 of 2,657 past-90° vertices sit under a cell from the crease; the other 2,121 are 11 to 94 cells
away. The counts are exact multiples of `n − 2` and the wedge is extruded along `z`, so there are
**one to three offending locations in the whole cross-section** and neither classifier found them.
**R-009** owns that.

#### The method finding, which is the one worth keeping

**The fixture agreed with the prediction to four decimal places and was wrong.** A wedge whose
bisector lies on a grid axis puts its crease in a symmetric relationship to the sampling, and the
worst vertex is then exactly the symmetric one — the average of two faces, at exactly half the turn.
The number is a property of that symmetry, not of the dihedral, and nothing about `75.0000` at three
resolutions looks like a fixture artefact. Only turning the fixture 17° showed it. Part 5 already
carries *"choose a fixture by searching for one that exhibits the property"*; this is the harder
version — **a fixture can exhibit the property too perfectly, and exactness is the tell rather than
the confirmation.**

### M-284 / P-15 — HELD, and the dual's IPC wall is one function (R-007)

**M.** `benches/experiment_p15.rs`, `docs/experiments/p-15.csv`, 16 rows, Ryzen 9 5900X, `f32`,
Surface Nets, no surface. **P-15 predicted more than half the cycles; the answer is 82%.**

#### The instrument, which needed nothing from the library

R-007 offered two ways into a private function and asked which was allowed: counter windows inside
`DualMesher::extract`, or an ablation seam. Neither was used. The first is impossible before it is
undesirable — `isomesh` is `no_std` and cannot make a Linux system call — and the second is public
API for one experiment.

The third way is that the stages have **different iteration counts, and the counts depend on the
grid's shape rather than only its size**: `sample` runs `S = ∏ size`, the two resizes and
`place_vertices` run `C = ∏ (size − 1)`, and `emit_quads` runs
`Q = Σ_axis (size[axis] − 1)(cells[u] − 1)(cells[v] − 1)`. On a cube `Q/C = 2.97`; on a slab two
samples deep it is `1.00`; on a **rod two samples deep on both minor axes the inner loops are empty
and it is exactly 0**, while `S/C` runs 1 → 4 across the same shapes. Thirteen shapes make the design
matrix separable and least squares reads the stages off it.

| shape | cells | `S/C` | `Q/C` | cycles/cell | instructions/cell |
|---|---|---|---|---|---|
| 193×193×193 | 7,077,888 | 1.02 | 2.97 | 156.04 | 205.29 |
| 513×513×9 | 2,097,152 | 1.13 | 2.74 | 153.21 | 199.29 |
| 1025×1025×3 | 2,097,152 | 1.50 | 2.00 | 126.06 | 179.88 |
| 1449×1449×2 | 2,096,704 | 2.00 | **1.00** | **66.63** | 145.92 |
| 500001×2×2 | 500,000 | 4.00 | **0.00** | **31.32** | 145.88 |

#### The decomposition, and the number that answers R-007

| per iteration | instructions | cycles | **IPC** |
|---|---|---|---|
| `sample` | 16.31 | 2.51 | 6.49 |
| `place_vertices` + resizes | 95.82 | 24.99 | **3.83** |
| **`emit_quads`** | 31.24 | 43.26 | **0.72** |

`r² = 0.9954` on instructions and `0.9990` on cycles, and the model **predicts rather than
describes**: a shape held out of the fit entirely, `385×385×17`, comes in at **+0.12% on instructions
and +0.50% on cycles**.

**At 193³, `emit_quads` is 45% of the instructions and 82% of the cycles.** So P-15 held, and by more
than it claimed.

**The punchline is the IPC column.** `place_vertices` — the dual's cell loop, which reads the same
eight corners Marching Cubes' march does — runs at **3.83**, which is Marching Cubes' own 4.04 (M-282).
`sample` runs at 6.49. Only `emit_quads` is slow, at **0.72**, and it is slow enough to drag the whole
mesher to 1.20. **The dual does not have an IPC problem; one of its four stages does**, and that
stage is 45% of the instructions and 82% of the time.

#### The fit-free cross-check, and where it disagrees

`1025×1025×3` and `1449×1449×2` have cell counts agreeing to **0.02%** and `Q/C` differing by exactly
1, so their difference in cost per cell is one `emit_quads` iteration with only a small `S/C`
correction. It gives **60.8 cycles and 42.2 instructions** against the fit's 43.3 and 31.2.

That is a 40% disagreement and it is reported rather than smoothed: 60.8 × 2.97 = 181 cycles per cell
would **exceed** the cube's measured 156, so the two-row estimate over-attributes — the two slabs
differ in more than `Q/C`. The fit is the conservative reading and the one with a validated
out-of-sample prediction, so **82% is a floor**.

#### A limitation the sweep found in itself

`500001×2×2` and `2×2×500001` have identical `S`, `C` and `Q` and differ by **62% in instructions**
(145.88 against 236.88). The model has no term for loop overhead, which dominates when the innermost
trip count is 2 — and `sample` iterates `x` innermost, so which axis is short decides how many times
the loop preamble runs. It does not affect the conclusion (both rods anchor `Q = 0` and the two
`500001×2×2`-shaped rows agree exactly at 145.88), and it bounds how far the coefficients should be
trusted on degenerate shapes.

#### O-11 is answered

*"Why does the dual topology go superlinear in `n³` while Marching Cubes does not?"* — asked at T-006,
half-answered at M-45, narrowed at M-279, and closed here. The dual carries a fourth stage Marching
Cubes does not have, `emit_quads`, which walks every grid edge on all three axes and loads both
endpoint samples **before** the sign test that would let it skip. It is `O(n³)` where the surface is
`O(n²)`, it runs at one sixth of the rest of the mesher's IPC, and it is 82% of the cost. The
superlinearity M-21 saw is that stage's working set growing past the caches while the other three
stay flat — and the remedy is **A-023**, not a vertex rule.

**A-023 landed and O-11 is now a smaller question rather than a closed one (M-285, M-286).** The
constant is gone: the dual is 2.5–3.1× faster, byte-identically, and Surface Nets beats Marching
Cubes at 48³. The *curve* is not: per-sample cost still rises 7.92 → 13.37 ns over 16³…256³, and
that residue is now visibly the cache — misses and cycles track at 8.4 cycles per miss where before
they did not track at all. What is left of O-11 is **A-024**, and specifically the **2.6× penalty at
exactly 128³** that a 64 KiB plane stride buys.

**A-024 landed and O-11 is closed (M-287).** Forcing the row length odd removed the aliasing at
every size — 128³ against its neighbours goes 3.37× → **1.01×** — and Surface Nets' per-sample
cost across 16³…256³ is now **8.71 → 9.70 ns, +11%**, against the +40% this question was raised
about. Two causes, both removed, and neither of them a vertex rule or an algorithm: a dynamically
indexed coordinate array and a power-of-two row stride.

### M-285 — the dual's axis had to be a constant, and that was 82% of its cycles (A-023)

**M.** `benches/family`, `benches/experiment_p15`, `benches/experiment_p12`, all re-run at one commit
with the clock at 4.14–4.26 GHz.

`DualMesher::emit_quads` took `axis`, `u` and `v` as **runtime** values, so `p[axis] = a` was a
dynamically indexed store and `p` could not live in registers: every iteration wrote three coordinates
to the stack and `linearize` read them straight back — a store-to-load chain the scheduler cannot
break. R-007 had already priced the stage at **82% of the dual's cycles at IPC 0.72** (M-284), beside
a cell loop doing more work per iteration at 3.83.

It is now three monomorphisations of one function with the axis a `const` generic. **Same three
passes, same order, same bounds, same triangles in the same sequence.**

| | before | after |
|---|---|---|
| `emit_quads` per iteration | 43.26 cycles, 31.24 instructions | **3.33 cycles, 12.02 instructions** |
| `emit_quads`' share at 193³ | 82% of cycles | **27–32%** |
| the whole mesher at 193³ | 156.0 cycles/cell | **37.1 cycles/cell** |

At 256³, one binary, one run:

| algorithm | before | after | speedup | IPC | × Marching Cubes |
|---|---|---|---|---|---|
| `surface_nets` | 693.8 ms | **224.4 ms** | **3.09×** | 1.20 → **2.77** | 5.43 → **1.72** |
| `dual_contouring` | 751.1 ms | **279.5 ms** | **2.69×** | 1.35 → **2.88** | 5.88 → **2.14** |
| `manifold_dual_contouring` | 771.0 ms | **305.6 ms** | **2.52×** | 1.42 → **2.90** | 6.03 → **2.34** |

**Marching Cubes, its decider variant, Marching Tetrahedra and subgrid Marching Tetrahedra are
unchanged to within noise**, which is the control: the change is inside `DualMesher` and only
`DualMesher`'s three entries moved.

**The mesh is byte-identical.** `golden_hashes_are_unchanged` passes untouched, and every triangle
count across the family matches the run before it. That was the acceptance: an optimisation that
changes the mesh is a bug in the optimisation.

**And Surface Nets is now faster than Marching Cubes at 48³** — 30.7 cycles per sample against 33.1,
at **IPC 5.29** against 4.19. ✗14 said Surface Nets *"never wins on Zen 3, at any resolution: 2.46×
behind even at 16³"*. That is no longer true, and the reason it was true was one missing `const`.

### M-286 — the misses M-279 measured as free were hidden behind the stall, and now they cost (A-023)

**M.** M-279 measured, on the pre-A-023 dual, that **a 2.4× swing in cache misses moved cycles by
0.4%** — three grids of 16.7 M samples differing only in axis order — and concluded that the LLC-miss
column tracked the dual's cost without accounting for it. **That conclusion was true of that code and
is false of this code**, because the store-to-load chain in `emit_quads` was stalling the loop deeply
enough to cover every miss underneath it.

The same control, re-run: misses **2.12 / 1.35 / 1.95** and cycles **37.4 / 33.4 / 35.7**. They now
track, at about 5 cycles per miss.

And the cubic sweep, which before showed cycles rising while misses did nothing useful:

| n | Marching Cubes cycles | Surface Nets cycles | Surface Nets IPC | Surface Nets misses |
|---|---|---|---|---|
| 48 | 33.1 | **30.7** | **5.29** | 0.612 |
| 96 | 30.8 | 36.8 | 4.28 | 1.658 |
| 127 | 31.0 | 33.2 | 4.72 | 1.496 |
| **128** | 31.2 | **84.3** | **1.86** | **4.724** |
| 129 | 31.2 | 33.0 | 4.75 | 1.483 |
| 192 | 31.6 | 43.1 | 3.60 | 2.540 |
| 256 | 32.0 | 56.7 | 2.73 | 3.720 |

Cycles and misses now move together: 48³ → 256³ is `+26.0` cycles for `+3.11` misses, **8.4 cycles per
miss**; 129³ → 128³ is `+51.3` for `+3.24`, **15.8 cycles per miss**.

**So the residue of O-11 is now, and only now, a cache story.** What A-023 removed was a constant, not
the curve: Surface Nets' per-sample cost still rises with `n`, from 7.92 ns at 16³ to 13.37 at 256³,
and that rise is the dense per-cell state outgrowing the caches exactly as M-279's working-set
reasoning said — the reasoning was right and the *evidence* for it was masked.

**The 128³ spike is now the largest single feature of the curve and is unambiguous.** 127³ and 129³
cost 33.2 and 33.0 cycles per sample with 1.50 and 1.48 misses; 128³ costs **84.3 with 4.72** —
a **2.6× penalty at one resolution**, on working sets 2% apart, because `n²·4` is exactly 64 KiB and
that is a cache-set aliasing period on this machine. Before A-023 the same spike was 1.24× and looked
like noise on a curve. **A-024** owns it, and it matters more than its size suggests: 128³ is the
canonical chunk size in voxel engines.

**A-024 has since landed and removed it (M-287).** 128³ now costs 36.65 cycles per sample against
127³'s 36.01 and 129³'s 36.20, and the miss rate at 256³ falls from 3.72 to **1.56** per sample —
so the residual cache term this entry identifies was, at every size and not only at 128, largely
the same aliasing.

**Method note, and it is the third time tonight.** M-279's null result — *"misses do not cost"* — was
a correct measurement of a machine whose bottleneck was somewhere else. A null measured under a
dominant confound is a statement about the confound. The confound has to be removed before the null
means what it says.

### M-287 — one bit of the row length was a 3.4× tax at the chunk size everybody uses (A-024)

**M.** `benches/a024_aliasing`, `docs/measurements/a024-aliasing.csv`, plus `family`,
`experiment_p12` and `experiment_p15` re-run.

`DualMesher::values` was laid out by the caller's shape, so its **row** stride was `size[0]·4` bytes
and its **plane** stride `size[0]·size[1]·4`. At 128 those are **512 bytes and exactly 64 KiB** —
a cache-set aliasing period twice over on this machine.

#### Diagnosed before anything was changed, by letting the caller arrange the pad

The aliasing depends only on the shape, and the shape is an argument. So adding **one sample** moves
the stride while changing the work by under 1%:

| | plane bytes | cycles/sample | vs 127³/129³ |
|---|---|---|---|
| 127³ | 64,516 | 33.10 | — |
| **128³** | **65,536 = 2¹⁶** | **108.51** | **3.37×** |
| 129³ | 66,564 | 31.39 | — |
| 129×128×128 (pad `x`) | 66,048 | 31.48 | **0.98×** |
| 128×129×128 (pad `y`) | 66,048 | 36.45 | 1.13× |
| **128×128×129 (pad `z`, control)** | **65,536 = 2¹⁶** | **107.89** | **3.35×** |
| 128×131×131 (512-byte rows only) | 67,072 | 36.74 | 1.14× |
| 256³ | 262,144 = 2¹⁸ | 54.07 | **1.39×** vs 255³/257³ |

**The `z` control is what makes this a measurement rather than a story.** It adds the same 0.8% of
work and does *not* touch either stride, and it keeps the entire penalty. And the two periods separate
cleanly: 512-byte rows alone cost 1.14×, the 64 KiB plane costs the remaining ~3×, and only padding
the **fastest** axis fixes both because `size[0]` appears in both.

#### The fix is `size[0] | 1`

Unconditional, because a pad applied only when the stride looks bad is a second layout reachable from
one call — the shape the one-path rule forbids. And **idempotent, which is the part that matters**: a
*fixed* pad of one would be worse than nothing, since it maps every `size[0] = 2ᵏ − 1` onto the stride
it is trying to avoid. `| 1` has no such image.

It cannot reintroduce either period: `4·odd` is never a multiple of 512, and `4·odd·size[1]` is a
multiple of 65,536 only if `size[1]` is a multiple of 16,384. The cost is one float per row when the
row is even — **0.8% of `values` at 128³** — and nothing at run time, because the multiply is by a
different constant rather than an extra one.

#### After

| | before | after |
|---|---|---|
| 128³ against its neighbours | 3.37× | **1.01×** |
| 256³ against its neighbours | 1.39× | **0.92×** |
| Surface Nets at 128³ | 48.51 ms | **18.35 ms** (2.64×) |
| Surface Nets at 256³ | 221.4 ms | **162.7 ms** (1.36×) |
| Surface Nets IPC at 256³ | 2.80 | **4.09** |
| cache misses per sample at 256³ | 3.72 | **1.56** |

Marching Cubes and Marching Tetrahedra move 1.02× and 1.01×, which is the control: the change is
inside `DualMesher`. **Triangle counts across the family are identical and
`golden_hashes_are_unchanged` passes.**

#### Taken with A-023, and this is the headline

| | this morning | now |
|---|---|---|
| Surface Nets at 256³ | 693.8 ms | **162.7 ms** — **4.26×** |
| Surface Nets IPC | 1.20 | **4.09** |
| `SN / MC` at 256³ | 5.43× | **1.26×** |
| per-sample cost, 16³ → 256³ | 29.63 → 41.35 ns (+40%) | **8.71 → 9.70 ns (+11%)** |

Not one triangle changed. And Surface Nets is now **faster than Marching Cubes** on small grids:

| n | 16 | 24 | 32 | 48 | 64 | 96 | 128 | 192 | 256 |
|---|---|---|---|---|---|---|---|---|---|
| Marching Cubes ns/sample | 9.37 | 8.87 | 8.29 | 8.03 | 7.74 | 7.41 | 7.48 | 7.62 | 7.67 |
| Surface Nets ns/sample | **8.71** | **8.59** | **8.05** | 8.24 | 8.42 | 8.53 | 8.75 | 9.53 | 9.70 |
| ratio | **0.93** | **0.97** | **0.97** | 1.03 | 1.09 | 1.15 | 1.17 | 1.25 | 1.26 |

**O-11 is closed.** *"Why does the dual go superlinear where Marching Cubes does not?"* Two reasons,
both now removed: a dynamically indexed coordinate array that stalled `emit_quads` on store-to-load
forwarding (A-023), and a power-of-two row stride that put the sample array on one cache set
(A-024). What is left is **+11% across four octaves of grid**, against Marching Cubes' own −18% then
flat. The question was asked at T-006 on 2026-08-13 and it took the counters to answer, but the answer
was two keywords.

**And it re-opens a design question rather than settling one.** ✗14 exists to say Surface Nets is not
the cheap default. Its triangle-count half (✗1, `2χ` more triangles) is untouched; its **cost half is
now 1.26× at the largest grid and a win below 48³**. Whether that changes the default is not a
measurement, and this file does not decide it.

### P-16 — registered for R-008, before it was measured

**What M-283 left, stated so it can fail.** 6,959 vertices under Marching Cubes and 4,868 under Dual
Contouring, on an exact convex wedge, carry an area-weighted normal **more than 90° from the field
gradient** — worst 128.0°. An area-weighted normal is a convex combination of incident face normals,
so it lies inside the cone those faces span, and for two planes meeting at `θ` that cone is at most
`(180° − θ)/2 ≤ 75°` wide. The obvious escape was measured and closed: the median `|f(v)|/h` at those
vertices is **0.0000**, so they are on the isosurface.

> **H.** Every such vertex lies on a grid edge at least one of whose incident cells **straddles the
> crease** — its eight corners do not all have the same nearer plane — so the phenomenon is two faces
> meeting inside one cell rather than a winding or ordering defect.

If that holds, it is M-15's *"any feature thinner than one cell forces two sheets through it"*
arriving on a sharp field, and it is inherent to one vertex per crossed edge rather than fixable.

**Falsified by** more than 5% of past-90° vertices whose incident cells all lie on one side of the
crease — which would make it a defect with a location and a fix, and a much better outcome.

**Records** `dihedral_deg`, `samples`, `past90_vertices`, `past90_in_straddling_cell`,
`offending_faces_per_past90_vertex`.

**The threshold is 5% and not 0% on purpose.** A cell is classified by its eight corners' nearer
plane, and a corner sitting within rounding of the bisector can be classified either way; a hard zero
would make the hypothesis fail on arithmetic rather than on geometry. 5% is far above that noise and
far below "a substantial share".

### M-288 / P-16 — ~~FALSIFIED~~, and the registered definition did not implement the registered claim (R-008)

> **⚠ THE VERDICT IN THIS ENTRY IS WRONG. P-16 HELD, at 0% — see M-289.** Every number below is
> counted over vertices whose reference gradient was noise. Corrected: **442 offenders in 6 rows,
> 100% of them on the crease**, only at `θ = 30°` and `60°`. What survives is the observation about
> the registered definition — it splits by the bisector rather than the crease, and it still agrees
> with the corrected classifier on every row.


**M.** `benches/experiment_p16.rs`, `docs/experiments/p-16.csv`, 63 rows: an exact convex wedge ×
7 dihedrals × 3 rotations × 3 resolutions, Marching Cubes, area-weighted normals.

**23 rows contain a past-90° vertex, 2,657 in total, and 536 of them — 20% — are on the crease.**
The registered falsifier was *"more than 5% whose incident cells all lie on one side"*. It is 80%.

| | rows | median distance to the crease |
|---|---|---|
| offenders **on** the crease | 6 | **0.69–0.90 cells** |
| offenders **elsewhere** | 17 | **11.14–93.81 cells** |

So the mechanism P-16 named is real and is a **minority**. Two faces meeting inside one cell does
produce past-90° normals — the six rows where it happens are unambiguous, the offenders sit under a
cell from the crease, and the control (vertices under 90°) is on the crease only 0–12.8% of the time,
so the classification is not vacuous. It is simply not what most of them are.

#### The registration's operational definition was wrong, and it did not matter

P-16 reads *"straddles the crease — **its eight corners do not all have the same nearer plane**"*. The
clause after the dash does not implement the claim before it: `d0 ≥ d1` splits by the **bisector**
plane, which runs from the apex through the *interior* of the solid, not by the crease. A cell deep
inside the wedge and nowhere near the surface straddles the bisector; a cell holding the crease need
not.

The registered test was kept and reported anyway, because a registration is not edited after its
experiment runs, and a corrected test — does the cell hold the crease line, which rotation does not
move — was reported beside it and labelled post-hoc. **They agree on the verdict of every one of the
23 rows.** The error was real and changed nothing, which is worth knowing in both directions: the
result does not depend on the mistake, and a wrong definition can survive a whole experiment
undetected because it happened to correlate.

#### What the counts localise, which is more than either classifier managed

Offender counts are **31, 62 at n = 33; 63, 126 at n = 65; 127, 254, 381 at n = 129** — every one an
exact multiple of `n − 2`. The wedge is extruded along `z`, so the sampling is identical in every
layer and any feature of the two-dimensional cross-section repeats once per layer. **There are
therefore one, two or three offending locations in the entire cross-section**, and they are 11 to 94
cells from the crease.

That is a small, specific, reproducible configuration rather than a diffuse effect — the same shape
as A-021's 314 non-manifold edges before M-276 named them. Finding it wants a constructed minimal
fixture, not a wider sweep, and that is **R-009**.

#### What R-008 concludes

The question was *"is the crease bridged by triangles that face somewhere else, and is that
inherent?"* The answer is that **the crease does that and accounts for a fifth of it**; the other
four fifths are somewhere else, at one to three points per cross-section, and this experiment did not
find them. M-283's reading — that the bridging is M-15's thin-feature mechanism on a sharp field — is
**narrowed rather than falsified**: it is one of at least two causes and not the main one.

### M-289 — the reference gradient was noise, and it falsified two hypotheses that were true (R-009)

**M.** `benches/r009_locate.rs`, `docs/measurements/r009-locate.csv`, and `experiment_p13` and
`experiment_p16` re-run.

**This entry reverses M-283's and M-288's verdicts. Read it before either of them.**

#### How it was found, which is the method and not an accident

R-008 left 80% of the past-90° normals unlocated and two classifiers saying only where they were not.
What it *had* bounded was tight: the counts were exact multiples of `n − 2` and the wedge is extruded
along `z`, so there were **one to three offending locations in the entire cross-section**. A-021 is
the model for that situation — it found its answer by printing a face-count histogram for a plain
half-space, not by widening a census — so R-009 dumped one configuration.

The answer was in the first six lines. Two cells per cross-section, **six** incident faces each,
**every face lying exactly on a plane** (worst face `0.00°`), no slivers, and a stored vertex normal
**exactly equal to a plane normal**. The mesh was correct in every respect. And the gradient it was
being compared against came back `[-0.5156, 0.8568, 0]`, which is neither plane normal.

#### The bug

`Wedge::gradient`'s exterior branch computes `away = q − dir·t` and normalises it. On a point lying
**on** a ray those two vectors are equal, so `away` is not zero — it is a cancellation residue of
order `ε·|q|` — and normalising it returns a **random unit vector**. The guard was `e > 0.0`.

Every Marching Cubes vertex is on the surface to within an ulp, and about half of them land
epsilon-*outside*, so **half the vertices near a face were compared against noise**. The fix is a
threshold relative to `|q|`, falling back to the plane normal, which is what the exterior gradient
converges to approaching a ray along the surface — the right answer rather than a tolerance.

#### What it reversed

| | as reported | corrected |
|---|---|---|
| past-90° vertices, Marching Cubes | 6,959 | **472** |
| past-90° vertices, Dual Contouring | 4,868 | **232** |
| rows containing one (of 168) | 75 | **8** |
| **P-13** | FALSIFIED | **HELD** |
| **P-16** | FALSIFIED at 4× its threshold | **HELD at 0%** |

**P-13 holds.** With the corrected gradient the worst angle tracks `(180° − θ)/2` and is bounded by
it in **136 of 168 rows**, falling monotonically with the dihedral — at rotation 17°, `n = 129`:
θ = 30° → 58.0 against 75 predicted, 45° → 50.5 against 67.5, 90° → 33.3 against 45, 120° → 25.3
against 30, 150° → 10.3 against 15, 170° → 2.8 against 5, and 180° → **0.00**. And it is
**resolution-invariant to four significant figures** — 58.00, 58.00, 58.00, 58.00 across 17³…129³ —
which is M-66's non-convergence, reproduced cleanly and now predictable from the field. That is
exactly what P-13 claimed.

**P-16 holds, at 0%.** Its falsifier was *"more than 5% of past-90° vertices whose incident cells all
lie on one side of the crease"*. The corrected measurement is **442 offenders in 6 rows and 100% of
them on the crease**, median 0.69–0.73 cells from it — and they occur only at **θ = 30° and 60°**,
the acute wedges, where the wedge is thin enough for two sheets to share a cell. That is M-15 on a
sharp field, exactly as registered.

#### What survives from M-283 unchanged

The `θ = 180°` control (**0.0000°** worst and mean at every resolution, because Marching Cubes is
exact on a linear field), the median of `0.000°` everywhere (the disagreement is confined to a
one-dimensional crease, which is what reconciles M-66's falling mean with its flat worst), and the
grid-aligned fixture reporting the prediction to **four decimal places** — 75.0000, 67.5000, 60.0000,
45.0000 — while rotated ones give 20.1 to 128.0. The fixture is still too perfect; what changed is
that the rotated values now *bracket* the prediction instead of contradicting it.

#### The rule

**A reference implementation used as ground truth needs the same scrutiny as the thing it checks.**
Every control in R-006 pointed at the mesh: rotate the wedge, offset the apex, add a no-crease case,
check the vertex is on the surface. All of them were about whether *Marching Cubes* was being measured
fairly. **None of them asked whether the `gradient` on the other side of the comparison was right**,
and it was wrong at precisely the points being measured — which is the only place a reference is ever
wrong in a way that matters.

The tell was available and was misread twice: an area-weighted normal cannot leave the cone its faces
span, so a past-90° reading was **arithmetically impossible** from the start. M-283 recorded that
impossibility, chased the one escape it could think of (the vertex being off-surface), measured it
closed, and concluded the geometry must be strange — rather than concluding that one of the two
quantities in the comparison had to be wrong. **When a measurement is impossible, suspect the
instrument before the world.**

### M-290 — the ambiguous face in the dual path, with the source read at last (A-022)

**M.** `benches/a022_decider.rs`, `docs/measurements/a022-decider.csv`, 144 rows: eight reference
fields × three resolutions × six extractor configurations, `f64`, no chunking and no weld (matching
A-021, because the weld can create a non-manifold edge — M-226).

#### The block was a lookup failure, not a paywall

A-022 sat blocked on *"the source is paywalled and not in the corpus"*, with `paper_download`
reporting *"No open-access PDF found"* for `10.1109/TVCG.2007.1012`. The paper is at
**`cs.wustl.edu/~taoju/research/dualsimp_tvcg.pdf`** — Tao Ju's own publications page, and **the exact
filename A-022's own text already named**. One web search found it.

Two obstacles were in the way and neither was the publisher. `paper_download` resolves a DOI through
arXiv, Unpaywall and provider resolvers; an author's personal copy is indexed by none of them. And the
server presents a valid InCommon certificate **without its intermediate**, so every TLS client refuses
it — fetched by supplying the missing intermediate from the certificate's own AIA URL, with
verification left on.

**It is still not in the corpus**, and `paper_download` still answers *"No open-access PDF found"* —
re-checked after reading it. So the recipe is recorded rather than left to be rediscovered:

```bash
curl -sSo inter.crt http://repository.emsign.com/certs/EEEMIncommonDVG2C.crt
openssl x509 -inform DER -in inter.crt -out inter.pem
cat /etc/ssl/certs/ca-certificates.crt inter.pem > chain.pem
curl -sSL --cacert chain.pem -o dualsimp_tvcg.pdf \\
  https://www.cs.wustl.edu/~taoju/research/dualsimp_tvcg.pdf
```

`--cacert` **adds** the missing intermediate to the system bundle; nothing is disabled, and the root
(`emSign Root TLS CA - G1`) is one the system already trusts.

#### What it says, and it answers the question A-022 could not

V-34 has the quotations. The criterion is **one vertex per cycle of a table whose ambiguous faces the
asymptotic decider has already resolved**, each edge owned by exactly one vertex. So A-022's framing —
*"whether MDC's own criterion is face-based or component-based decides whether this is a new rule or a
bug in the existing one"* — has a third answer: it is **component-based over a face-resolved table**,
and the dual walk needs no rule of its own because the ambiguity is settled upstream.

#### A-022's acceptance was unreachable, and the paper says so

It asked for M-276's 314 non-manifold edges to go to **zero under Surface Nets and Dual Contouring**.
Those are one vertex per cell, and §3 says *"DC leads to nonmanifold vertices and edges for **all** of
the ambiguous sign configurations"*. Measured here: 1,128 edges across the 24 configurations, against
Marching Cubes' 0. **The 314 is the literature's own prediction and not a defect to remove** — removing
it means splitting the cell's vertex, which is Manifold Dual Contouring, which this crate already has.

#### The decider helps by a fifth and does not eliminate

`manifold_dual_contouring` **defaults to `FaceAmbiguity::Separate`**, which is not the table the paper
specifies. Switching to the decider-modified one takes it from **143 to 114** non-manifold edges over
the 24 configurations — `noise_cavity` at 33³ from 64 to 40, at 49³ from 53 to 49, at 65³ from 26 to
25. Better, sourced, and still not zero, which is ✗19.

**And the residue is one field.** Manifold Dual Contouring is manifold on **seven of the eight**
reference fields under both rules; every one of the 143 is `noise_cavity`. That is the field A-002e
added precisely because none of the other seven produces a cell with an **interior** ambiguity
(M-208), which a *face* decider cannot see by construction. That is the next question and it is
**A-025**, not this ticket.

#### The rule

**A paywalled DOI is not a missing paper. Check the author's own page before taking a rule-5 stop.**
V-31 earned the sibling of this — *"a deleted repository is not a missing source; ask an archive, and
say which archive"* — one ticket earlier, and this is the same failure with a different resolver:
`paper_download` answers a question about **DOI resolvers**, not about whether the paper is readable.
A-022 was blocked for a day on a file that a search engine returns first, under the filename the
ticket had already written down.
