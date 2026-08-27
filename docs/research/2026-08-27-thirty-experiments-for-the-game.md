# Phase 23 read back, and thirty experiments for the game

**Written:** 2026-08-27 · **Reads:** `FINDINGS.md` Phase 23 in full (lines 10194–12114), the twelve `docs/experiments/p-6*.csv` and `p-7*.csv`, `BACKLOG.md`'s 22 open rows, the six gameplay research docs, and roughly 220 external fetches across rendering, physics, navigation and multiplayer.

**Two halves.** Part 1 reads Phase 23 back — what held, what the falsifications actually mean for a game, and five defects the entries carry. Part 2 onward is what was asked for: **thirty registrable experiments aimed at making a meshing-based game better**, numbered `P-73` … `P-102`.

**The lens has changed and the filter with it.** Phase 23 was about mathematics and logic. This set is about a game, so every row below is scored on *what a player would notice*, and a row that improves a number nobody can see is not here. Everything was checked against the 2026-08-18 novelty table's fifteen rows and its rejected list, the mechanics dossier's Tier-3 losers, `BACKLOG.md`'s "deliberately not in scope", and the 2026-08-23 dossier's twenty-five foreclosures. Where a row overlaps something already proposed, it says so and states what is new about it.

---

## Part 1 — Phase 23, read back

Ten of twelve ran. **Four clean wins, one real bug caught, five informative falsifications, two blocked on acquisition.** The process result is the headline: every Phase-23 CSV resolves to an ancestor of HEAD and every one was written against a clean tree, which is the exact debt the 2026-08-26 audit opened on — and `scripts/csv_provenance.sh` is now a `preflight.sh` step rather than a recommendation.

### What landed

**`M-372` — Marching Cubes is now bit-exactly equivariant under all 48 octahedral elements**, `worst_component_ulp` **0** on all 16 rows, `gyroid` and `noise_cavity` included. The crossing moved from `a/(a−b)` anchored at the lower corner to a signed offset from the edge midpoint, `cube::edge_crossing` became `cube::edge_offset`, and the 216 golden hashes were rebaselined in the same commit. For the game this is not a rendering change — it is a **test oracle with 48 independent relations per fixture** where there were six, and it is the reason `P-100` and `P-101` below are worth registering.

Both of C1 and C2 were falsified, and both falsifications are more informative than the clauses were. C1 asked for all four primal extractors: `marching_tetrahedra` improves from {6, 12} to a flat **12 of 48 and stops**, because a six-tetrahedron decomposition of a cell is not octahedrally invariant — its diagonals cut different edges after a relabelling, and **no placement rule can reach that**. And `subgrid_marching_tetrahedra` did not move at all, because it **does not use the linear crossing**: `subgrid::roots::all_roots` bisects the real field. That was found by C2, not by reading the code.

**`M-374` — `O-12` is settled for Marching Cubes, by exhaustion, in 21.5 seconds.** The third mechanism `O-12` named — a vertex whose two face groups sit in different cells — **does not exist for Marching Cubes**. Every face incident to an edge vertex comes from one of the four cells sharing that edge; those four cells are 18 corners; 18 corners is 262,144 patterns; zero link defects on four independent magnitude seeds, with the pre-fix single-apex fan reproducing **5,302** defects in the same walk. What remains open is the dual family, which needs the 3 × 3 × 3 block at 2²⁷ — filed as `R-072`.

**`✗50` / `M-373` — the harness found a panic before it measured anything.** `MAX_PATCH_TRIANGLES = 24` was a *sampled* maximum plus two; `fan_tunnel` buffered its own output at the derived `40`; the caller's array in `mod.rs` used the smaller one, and **the smaller one was the one that indexed**. Sign pattern `0x96` — MC case 13 — emits 26 triangles. Index out of bounds, in release, on an ordinary trilinear cell. **The margin was the entire safety argument and it was two triangles wide.** In a shipped game that is a crash on a cave wall.

**`M-379` — the case table is proved indexable**, over all 256 sign patterns in 2.04 s and all 16,384 `(case, mask)` pairs in 161 s, with four properties: shape, every code nameable, every named edge cut, no degenerate triangle. The interesting part is the failure that preceded it: pointing Kani at `triangulate(segment_links(case, joined))` exhausted **32 GB** — 909,347 VCs, 885 s of symbolic execution — because `segment_links`' twelve-edge walk is data-dependent and symbolic execution cannot merge paths. The registration's falsifier said "a property that cannot be expressed against the sign abstraction means the abstraction is wrong". It was a third shape: **the abstraction was right and the function was the wrong thing to bit-blast.**

**`M-377` — 51× from one knob.** Chunk granularity swept 2³–64³ at a fixed 128³ world: 4³ costs **8.61 ms** against 64³'s **439.51 ms** on `gyroid`, and the optimum is *interior* — 4³ beats 2³ by 13.6% and 8³ by 11.5%. Had the sweep stopped at the registered 8³ floor it would have reported the boundary and been wrong about the shape. The mechanism is arithmetic: at 64³ eleven edits re-mesh 8,126,464 cells, at 4³ they re-mesh 25,344 — **321× fewer cells for the same edits**. GVDB's 2048³ regime made a large volume look like the reason for its 256× spread; it is not. The cost is re-meshing a whole brick to service a local edit, and that does not care how big the volume is.

**This is the single most game-relevant result in the phase**, and it comes with a bill the registration did not name: 4³ carries **51.6% more vertex data** than 64³ for the same 53,110 distinct surface points. `P-93` below is that trade-off, registered.

**`M-378` — the crate can now say where it is correct.** The Plantinga–Vegter gradient predicate certifies **zero unsound cells over 2,389 tunnel and twelve-vertex cells**, while refusing **95.07%** of an adversarial random population. It is not a predicate that says yes to everything. C3 failed at 21% of extraction, but two thirds of that is **re-reading corners the extractor already has** — a fused version pays 0.0658, and the entry correctly refuses to substitute that for the verdict. `P-98` registers the fused version.

### What the falsifications mean for the game

**`✗51` — autovectorisation is dead here, and Amdahl killed it, not LLVM.** `total %ymm` across all eleven monomorphisations of both loop shapes is **zero**. But the clause was unreachable before that mattered: the sample loop is **1.23 ns/sample against an extraction marginal of 10.68** — 11.5% of the quantity C1 was denominated in, so halving it gives 1.06×, not 2×. The new Part 5 rule is the durable half: *a clause stated as a ratio of a total must name the share of that total it can move.* Every speed row in Part 2 below states its share.

The other durable half was not in any clause: there were **three copies of the sample loop** (`marching_cubes/mod.rs`, `dual.rs`, `marching_tetrahedra.rs`), now one `sdf::sample_grid`.

**`✗52` — the 83% is real and the decomposition was wrong.** Synchronisation reproduces at **86.5%** and `M-160`'s flat CPU time at 0.16 ms over a 60× cell range. But the largest component is not the count wait — it is the **geometry copy at 0.7075 ms**, 60% of the whole extraction. `extract_indirect` removes the count wait entirely, which is 31% and not 60%, and it cannot remove the copy because it removes it by *not delivering the geometry to the CPU at all*. The consequence splits cleanly by consumer, and this is the sentence that matters for a game: **for a consumer that only draws, `extract_indirect` removes 100% of the synchronisation it needs to; for a consumer that needs the bytes — collider generation, `M-116`'s decomposition — the copy is unavoidable.** Which is why `P-85` below points the same attribution instrument at the collider.

**`✗54` — the subgroup scan is real, reproduces the literature, and is not worth landing.** 1.4376× at 2²¹ elements on an RTX 3090, sitting in the middle of the six-device band Smith, Levien & Owens published — an independent reproduction on a seventh device. Applied to the shipped path it moves `gpu_total_ms` by **1.33%**, because the upload is 87.50% and the scan is 4.37%. Not landed. Two naga realities strengthen it: `enable subgroups;` is *rejected* by naga 29.0.4 while the intrinsics work without it, so the shader compiles on `wgpu` today and **would not compile on a spec-conforming WGSL implementation**.

**`✗53` — the under-resolution witness is sound at `k = 17` and not at `k = 5`**, and all eight false negatives at `k = 5` are on `noise_cavity`, the only field with genuine sub-cell noise. Two registration defects are worth carrying forward. C2's second half asked for a false-positive rate that *falls* with `k`, which is **mechanically impossible** — adding a sample point can only add an opportunity to disagree, so `flagged` is monotone non-decreasing. And C3's metric was taken over *all* grid edges, so it measures how much of the **volume** has turning gradient, not how badly the **surface** is under-resolved — which is why `thin_plate`, a thin plate in a mostly-empty box, ranked eighth. Its real signature is that its false-positive rate is **the only one that rises with resolution**. `P-99` registers the reformulation the data supports.

**`✗55` — the coefficient is 3, and the derivation says why.** `fl(a + b)` carries a relative error up to `u`, so an absolute error up to `u·|a + b|`; dividing by `(a − b)` makes that `u·|a + b|/|a − b|`, and `|a + b|/|a − b|` **is** `2|d|`. The numerator contributes two units, the quotient one. The registration counted roundings instead of propagating them. Ten violations at coefficient 2, zero at 3, over 19,415 crossings — and all ten are on `gyroid`, `noise_cavity` and `fbm_terrain`, none on the five smooth fields, which is a mechanism rather than a tail.

### Five defects in the entries

Verified against the CSVs. None reverses a verdict; all are the kind of thing this ledger exists to catch.

**D-1. `✗49` states "zero on 9.2 million straddling pairs".** The premeasure block's own `pairs` column sums to **5,800,000** (2 × 2,000,000 + 6 × 300,000), and the entry prints those row counts three lines above. Off by 59%. The claim it supports is fully sound; the headline count is not.

**D-2. `✗49` states "2,285 of 28,124 cut edges moved — 7.3% of them".** Both counts reproduce exactly and the ratio is **8.12%**. (7.24% is the 33³-only subset; the mean of the per-row fraction column is 10.9%.)

**D-3. `✗52`'s C1 and C2 figures are not in the file its `M.` line names.** `p-71.csv` at `2fc75b4` reads copy **0.666295**, map-wait 0.319133, `synchronisation_removed_share` **0.323852**; the entry quotes 0.7075 and 0.3098 from the superseded `d3b79e7` run. The two-run split *is* declared — but the reconciliation sentence claims the amendment "reproduced them to within 0.5%", and its own two cited numbers are **−5.83%** and **+4.55%**. The stated 0.31–0.45 band is likewise nowhere in `docs/experiments/`; the committed file's band is 0.30–0.33.

**D-4. `P-70`'s C3 is a HELD with no instrument behind it.** "The fallback is exercised, not merely present" is recorded as the bare word HELD, and `p-70.csv` has **no column that could have made it false** — no wasm arm, no forced-`SUBGROUP`-off arm, only `hillis` and `subgroup` rows which both always run. The registered falsifier was explicitly *"`M-44`'s zero-that-could-not-have-been-non-zero in a new place"*. In a phase this careful about `M-44`, this is the weakest thing in it. C2 is a milder case of the same: registered as "bit-identical output on all eight fields — same triangles, same order", measured as prefix-scan output at four sizes, scored HELD. The entry declares the gap and scores it anyway.

**D-5. `✗43`'s 8,064-case sweep is still unreproducible and now has no ticket.** The audit offered two fixes; the repo took the second and downgraded the entry's own text to say the sweep exists in no bench, no test and no CSV. Honest, and it leaves the evidence that the per-ring apex fix *generalises* permanently absent. `P-102` registers it properly.

**Minor, for completeness:** `✗51`'s "11.6%" is 11.53%; `✗52`'s "86.5%" is 86.44%; `M-377`'s `spread` column carries gyroid's 51.045955 on the `fbm_terrain` rows, so the entry's 47.19× is derivable but not recorded; `p-61.csv`'s obviously-named `c1_rows_at_48` (36 of 112) is not the clause's column (`c1_can_fail_rows_at_48`, 28 of 98) and wants a rename. And the phase header's claim that all twelve ids were registered "each in its own commit, before their harnesses exist" is not true of `P-65` and `P-67`, which landed together in `d3b79e7` after the commit that declared them blocked. Nothing about falsifiability is compromised — neither ever ran.

### The two acquisitions that blocked

`P-65` (MCPro) and `P-67` (reduced affine arithmetic) are both blocked on a paper, not on a ticket, and the MCPro blocker produced the most useful method finding of the phase. **`M-371`: `catalog_read` is not a presence oracle either.** MCPro is in the corpus as a 279,882-byte SciTePress landing page whose markdown is **383 characters** — a title, a topic list and a page range — and `catalog_read` reported `markdown_path` set, conversion complete, `chunks_indexed: 1`. Every field a present paper has. Three mechanical discriminators are now known: markdown length (383 against Finken's 37,165), `chunks_indexed` of 1 against a real paper's 12, and `pdf_path` not ending in `.pdf` — the third being exact rather than statistical. That check reported sound for thirteen days on the acquisition gating the largest item in the phase.

---

## Part 2 — Where these thirty come from

Three sources, in descending order of how much they contribute.

**The gap by omission.** Across all six gameplay research docs there is **no proposal touching global illumination, shadowing, volumetrics, atmospherics, or a material system beyond the phase field, the pseudonormal, the off-surface shading canary and decals-from-the-log-map**. Lighting appears twice: once as a failure symptom (dual-contouring self-intersection reads as broken lighting) and once as an analogy (a static acoustic bake is like static lighting). That is a large hole in a document set that has costed reverb, structural collapse, hydrology, speleothems and modal audio to three decimal places — and it is a hole in the half of the engine a player looks at for every frame of the game. Eight of the thirty are there.

It is also the place where this crate's constraints buy the most and where the literature has measured the least. Every cost in the SDF-rendering literature is dominated by *constructing* the field: RTSDF pays 2.09 ms per frame jump-flooding a 128³ SDF before it traces anything, and nvblox pays 0.4 ms integrating a TSDF to reach 3 ns per distance query. **This crate already has the field; it is the source of truth the mesh was made from, and `isomesh-gpu` already ships a `jump_flood` module.** AMD's Brixelizer GI is structurally the same idea — cascades of sparse SDFs, voxelised every frame — and has **zero published numbers**. That is a benchmark-shaped absence.

**Two corrected citations, both of which change what a planned feature should be.** The mechanics dossier proposes the **angle-weighted pseudonormal** as a ~40-line improvement. Jin, Lewis & West (`10.1007/s00371-004-0271-1`, 140 citations) measured six connectivity weightings against the analytic normal and found that on **marching-tetrahedra output the median discrepancy is 5–20° for every one of them, even at the highest resolution**, and say plainly that none does particularly well — they attribute it to spatial aliasing in irregular implicit tessellation, which is the same cause as `M-72`. `P-73` puts that to the crate's own fields. Separately, `A-027`'s cut-and-assign is Müller, Chentanez & Kim's 2013 preserve-convexity-by-construction, which the ticket does not cite and which reports **sub-10 ms fracture on 2013 hardware** — a number worth having next to `M-116`'s 241–272 ms.

**Phase 23's own leftovers.** Five rows the entries named as follow-ons and did not file, plus the audit item that was downgraded rather than fixed.

### What was checked and deliberately not re-proposed

The banked list (Sabine reverb, modal analysis, masonry statics, heat-method/CPM routing, medial axis and weak-feature-size certificates, contour-tree maintenance, persistence-thresholded ambiguity, union-find air connectivity); the novelty table's fifteen rows and its rejected list; the mechanics dossier's Tier-3 losers (fracture modes on carved geometry, full-bandwidth modal audio, handholds from solid thickness, ligament severance, 3D repose-angle relaxation, buckling eigenproblems, topology whorls, freeze–thaw); `BACKLOG.md`'s "deliberately not in scope" (`O-17`, `O-18`, Nanite-style mesh-space cluster simplification, networked/concurrent editing, neural extraction); and the 2026-08-23 dossier's twenty-five foreclosures. Where a row below sits next to one of those, its opening line says which and what is different.

**One boundary worth restating, because three rows below sit near it.** `BACKLOG.md` closes networked and concurrent editing out — *"closed out, not deferred… a networking layer needs sockets, clocks and a session model — none of which belong in a `no_std` crate whose public API is `[f32; 3]`."* `P-94`, `P-96` and `P-97` do not propose a networking layer. They measure properties of the **edit log**, which the crate owns and which serves saves and undo whether or not anything is ever sent over a wire.

---

## Part 3 — The thirty

Numbering continues from `P-72`. Protocol as always: the id goes into `crates/isomesh/src/experiment.rs` in its own commit before any harness, `falsified_by` is required by the type, numbers land in `docs/experiments/p-n.csv` with a resolving SHA against a clean tree, and an `E×n` row is owed either way. Per `✗51`'s new Part 5 rule, **every clause stated as a ratio names the share of the total it can move.**

### Group A — rendering and shading, the gap by omission

#### P-73 — the angle-weighted pseudonormal is the wrong fix, and the field's gradient is both better and the only bit-exact option

**The hook is a proposal in this repository's own docs.** The mechanics dossier proposes the angle-weighted pseudonormal at "<2° difference on good triangles, >15° on radius-ratio <0.15, ~40 lines". Jin, Lewis & West 2005 scored six connectivity weightings — equal, angle, sine-and-edge-length-reciprocal, adjacent-triangle-area, edge-length-reciprocal, sqrt-edge-length-reciprocal — against the analytic normal, and on **marching-tetrahedra output every one of them has a median discrepancy of 5–20°**, at the highest resolution they tested. Their diagnosis is spatial aliasing in irregular implicit tessellation, which is `M-72`'s mechanism seen from the shading side. `NormalStrategy` already has `AnalyticGradient` and `CentralDifference`; this measures whether a third is worth having.

**C1.** On eight reference fields at 33³/65³, every connectivity weighting has a median angular error against the analytic gradient of at least **3×** `CentralDifference`'s at the cell size, and `thin_plate` and `noise_cavity` are the worst. *Falsified by:* any weighting within 3×, which would mean the 2005 result does not transfer to trilinear output and the dossier's proposal is right.

**C2.** `mean |f(v)|` — the off-surface shading canary the dossier proposes in two lines — **predicts** the per-vertex angular error, rank correlation above 0.7 on at least six of eight fields. *Falsified by:* below 0.7 on three or more, which would mean the canary and the error are different phenomena and the canary is not a proxy for shimmer.

**C3 (the determinism argument, and it is independent of quality).** Angle weighting needs `acos`. `libm`'s `acos` has no architecture selection, so a connectivity-weighted normal is a golden-hash liability on a crate that commits 216 of them; a central difference is a fixed sequence of subtractions and one normalise. Measured: the connectivity route moves at least one hash between the M5 and Zen 3, and the gradient route moves none. *Falsified by:* zero hash movement on the connectivity route, which would be the more interesting result and would remove the determinism objection entirely.

#### P-74 — ambient occlusion and soft shadows traced against the field the mesh came from

**Nothing in six gameplay documents proposes a lighting feature.** This is the cheapest one available, because the expensive half is already paid. RTSDF (`10.5220/0010996200003124`) measures, on an RTX 2080 Ti at 1024²: **jump flood 2.09 ms, ray trace 4.60 ms, total frame 10.22 ms / 97 fps** — and 2.09 ms of that is *building* a 128³ SDF the game does not have. nvblox (`10.1109/icra57147.2024.10611532`) measures **0.8–7.3 billion distance queries per second, about 3 ns per query**, on an RTX 3090 Ti, after paying 0.4 ms per frame to integrate the TSDF. This crate's field is analytic, resident, and already evaluated on the GPU at `M-155`'s 0.54 ms for a 129³ grid; `isomesh-gpu/src/jump_flood` already exists.

**C1.** Screen-space AO by sphere-tracing the resident field costs under **2.0 ms** at 1920×1080 on the Zen 3 / RTX 3090 rig, with 8 rays per pixel and a 16-step cone march. *Falsified by:* above 2.0 ms — the share it must move is the whole AO budget, which a mesh-based SSAO pass currently pays in full.

**C2 (the quality claim, and the one worth the experiment).** Against a mesh-based SSAO baseline on `gyroid` and `fbm_terrain`, field-traced AO has **no** darkening error at chunk seams and no haloing at silhouettes, measured as mean absolute difference against a 512-ray offline reference. SSAO's seam and halo error is non-zero on both fields. *Falsified by:* SSAO matching it, which would mean the field buys nothing a depth buffer does not already have.

**C3 (the mechanism that makes it worth a paper).** The construction cost RTSDF pays is **zero** here — measured as the difference between total frame cost with the field already resident and the same trace after a jump-flood rebuild. Predicted: the rebuild is 2 ms and the resident path is 0 ms of it. *Falsified by:* a non-zero resident cost, which would mean the field is not actually reusable in the form the tracer wants and there is a conversion nobody costed.

#### P-75 — material weights carried at the vertex, against the edit log walked at the fragment

**`M-138` measured the price of exact paint: every sample walks the edit log, at 2.33× the cost per chunk.** That is per *sample*. A shading pass pays it per *fragment*, every frame, at screen resolution. The alternative is bounded per-vertex material weights — an `[f32; 4]`, L1-normalised, the de facto standard in voxel renderers — computed once during extraction from the same walk. Nobody has published the comparison; the only public arithmetic for a real pipeline is a blog post costing 4 materials × 3 triplanar planes × 5 maps = **60 texture fetches per fragment**, with no timing behind it.

**C1.** Per-vertex weights, computed inside the existing extraction walk, cost under **5%** of extraction — the share is the whole vertex-attribute stage, currently zero. *Falsified by:* above 5%.

**C2.** At 1920×1080 on `game_dig`'s scene, per-vertex weights beat a per-fragment log walk by at least **4×** in frame time, and the gap widens with edit-log length across `M-50`'s four log buckets (1–15, 16–30, 31–45, 46–60 brushes). *Falsified by:* under 4×, or a gap that does not widen — the second would mean the fragment path is not paying `M-138`'s cost and the premise is wrong.

**C3 (what it costs).** Interpolating weights across a triangle is not exact: measure the material misclassification rate against the field at 10⁶ surface points. Predicted under **2%**, concentrated within one cell of a material boundary. *Falsified by:* above 2%, or misclassification that is not boundary-local — the latter would mean the weights are not a resampling of the field but a different quantity.

#### P-76 — one stochastic triplanar plane instead of three, and what it costs the temporal budget

**Triplanar and stochastic filtering compose multiplicatively and nobody states the product.** Three planes × three stochastic taps is nine fetches per map; `bevy_isomesh/examples/triplanar.wgsl` already pays the three. Stochastic Texture Filtering (`arXiv:2305.05810`) and Heitz–Neyret histogram-preserving blending (`10.1145/3233304`, ">20× faster" than procedural-noise state of the art, **hardware not named**) both trade fetches for temporal accumulation. The conflict is the finding worth registering: **a destructible world has already spent its temporal budget**, because geometry that changes rejects history.

**C1.** Selecting one triplanar plane per pixel stochastically drops the fetch count by exactly **3×** and the fragment cost by at least 2× at 1920×1080. *Falsified by:* under 2× — the share is the material stage, which `P-75`'s C2 will have measured.

**C2 (the clause the whole row exists for).** With TAA resolving the stochastic choice, the ghosting cost under an active dig is **worse** than the fetch saving is worth: measured as mean absolute error against a 3-plane reference in the 8 frames following an edit, above the error of the same scene with no digging. *Falsified by:* no difference between digging and static, which would mean history rejection is not the bottleneck and stochastic filtering is free here after all.

**C3.** Biplanar mapping — two planes rather than three, no stochastic term, no temporal debt — gets **at least half** the saving of C1 with none of C2's cost. *Falsified by:* under half, which makes the honest recommendation "keep three planes".

#### P-77 — how much temporal history a dig destroys, and whether 0.2 ms buys it back

**The one measured mitigation in this area is cheap and the problem it mitigates has never been measured under destruction.** k-DOP clipping (`10.1145/3681758.3697996`) replaces TAA's axis-aligned neighbourhood clamp with a tighter k-discrete-oriented polytope at **0.2 ms overhead** (GPU not named in the abstract), and its framing is that prior bounding-box methods are "only situationally effective". Separately, Epic's account of virtual shadow maps in Fortnite says animated deformation invalidates shadow pages and that they **abandoned caching for directional sun shadows entirely** — with no invalidation cost published anywhere.

**C1.** In `game_dig`, the fraction of TAA history samples rejected in the frame after a brush stroke is at least **5×** the steady-state rejection rate, and the elevated rate persists for at least 3 frames. *Falsified by:* under 5×, or a single-frame spike — either would mean destruction is not a temporal problem and this whole row and `P-76`'s C2 close together.

**C2.** k-DOP clipping recovers at least **half** the rejected samples at under 0.3 ms. *Falsified by:* under half, which would mean the rejections are genuine disocclusions rather than clamp conservatism, and no clipping scheme can help.

**C3 (the one a player sees).** The rejection is **spatially concentrated at the brush**, not global: over 80% of rejected samples fall within the brush's screen-space bounding box dilated by one dig radius. *Falsified by:* under 80%, which would mean an edit disturbs the whole frame and a localised mitigation cannot work.

#### P-78 — how many light probes one dig invalidates

**Probe-based GI is the only GI family that tolerates geometry changing every frame, and its cost under a *destructible* world is unmeasured.** The best rigorous measurement available is a bachelor's thesis (Dell'Ova, BTH 2025, **DOI unverified**) finding that **87–96% of radiance-cascade frame time is the gather pass** — i.e. cost is dominated by tracing the scene, not by the cascade structure. That matters here for one reason: this crate can trace an SDF instead of a BVH, which is what `P-74` measures. But before any of that, the question a game asks first is *how much of the cache does a dig dirty*, and the crate already owns the instrument — `M-311`'s dirty-cell set and `M-314`'s edit-proportionality decomposition.

**C1.** The probe invalidation set is **edit-proportional**, not volume-proportional: the number of probes whose visibility changes tracks the brush's dilated support with a constant factor under 4, across `M-50`'s four edit-log buckets and three probe densities. *Falsified by:* a factor above 4, or a count that grows with world size at fixed edit size — the second kills probe GI for this game outright and is worth knowing in an afternoon.

**C2.** Tracing the field beats tracing the extracted mesh for probe updates by at least **3×**, at equal ray count. *Falsified by:* under 3× — the share is the gather pass, which the thesis puts at 87–96% of the total.

**C3.** A dig that opens a new air component (`M-311`'s union-find merge, the banked breakthrough event) invalidates **strictly more** probes than a dig of the same volume that does not. *Falsified by:* no difference, which would mean the topological event has no lighting signature and the banked breakthrough event cannot drive a lighting response.

#### P-79 — shadow-map page invalidation per edit, which nobody has published

**Epic's own answer to deforming geometry was to stop caching.** The Fortnite VSM write-up gives exactly one quantitative figure — a light-loop optimisation from 1.56 ms to 1.08 ms — and describes invalidation qualitatively: sun movement causes "quite significant shadow page table changes frame to frame", animated deformation invalidates pages, and directional sun shadows were left uncached. **For a world where the geometry is being dug away, that is the relevant precedent and the number is missing from the literature entirely.**

**C1.** In `game_dig` with a cached shadow atlas, one brush stroke invalidates a page count proportional to the brush's **projected** area from the light, with a constant under 3 — not to the brush volume and not to the scene. *Falsified by:* a constant above 3, or scene-proportional invalidation.

**C2.** Invalidating only the pages the brush's light-space bounding volume touches produces a **pixel-identical** shadow to a full re-render, on all eight reference fields. *Falsified by:* any difference, which would mean the conservative bound is wrong and localised invalidation is unsound.

**C3.** The saving is worth having: cached-with-invalidation beats uncached by at least **2×** in shadow cost under a continuous dig at 12.5 strokes/second (`game_dig`'s throttle). *Falsified by:* under 2×, which vindicates Epic's decision and closes the direction with a number they did not publish.

#### P-80 — the LOD residual as a normal map, so coarse chunks keep their detail as shading

**`M-72` is the finding this is built on and it is the one nobody has exploited.** Sub-cell features do not vanish under coarsening, they **alias** — `thin_plate` goes 4,088 → 1,016 → 248 → 56 triangles across LOD 0–3 — and the dossier's response is to make the disappearance *authorable*. There is a second response: the difference between the fine surface and the coarse one is a **vector field the crate can evaluate**, because it owns the analytic field at both scales. Bake it into a tangent-space normal map per coarse chunk and the geometry fades while the shading does not. The sweep found **no measured study of normal-mapping detail onto a coarse isosurface** at all.

**C1.** For each coarse-LOD vertex, the direction from the coarse surface to the nearest fine-surface point is computable from the field alone — one gradient and one root refinement, no fine mesh — and agrees with the true nearest point to within **0.1 cells** on 95% of vertices at LOD 1 and 2. *Falsified by:* under 95%, or an error that grows with LOD level faster than the cell size.

**C2.** A normal map baked from that residual reduces the perceived detail loss: mean angular difference between LOD-2 shading-with-map and LOD-0 shading is under **10°** on `fbm_terrain`, against over 25° for LOD-2 without. *Falsified by:* over 10° with the map, or under 25° without — the second would mean coarsening does not cost shading detail on this field and the premise is wrong for terrain.

**C3 (where it must fail, stated in advance).** It cannot work on `thin_plate` or `gyroid` at LOD 3, because a normal map cannot restore a *silhouette* and both fields lose topology rather than curvature. Predicted: angular difference above 25° on both, with the map. *Falsified by:* the map working there, which would be a much stronger result than the row claims and would need explaining.

### Group B — collision and physics

The framing number for this whole group: **the collider check is 45% of the pipeline, larger than contouring**, and it is described in the crate's own docs as "the unexamined 45%". Two independent voxel engines have now landed on the same shape — Godot Voxel Tools' performance docs say *"creating a collider from a mesh is actually much more expensive than meshing itself (about 3 to 5 times)"*, with no absolute timings and no hardware. That is corroboration, not measurement.

#### P-81 — a capsule against the field, instead of a capsule against the triangles

**The docs already propose SDF collision and cite the wrong paper for the game's actual case.** The cited source is Liu et al. `10.1016/j.cagd.2024.102305`, which is SDF-**vs**-SDF for two comparable solids — a strictly harder problem needing interval arithmetic, which is a determinism liability. The game's case is a small analytic proxy against one enormous static field, and that is Macklin, Erleben, Müller, Chentanez, Jeschke & Corse, `10.1145/3384538`: per-element local optimisation between an SDF isosurface and mesh elements, comparing projected gradient descent, Frank–Wolfe and golden-section search, with **GSS winning on the 1-D edge problem**. Their decisive line, on 129k-triangle rigid shells: *"this mesh-based collision took approximately 15 ms per-step, compared to < 0.5 ms using SDF-based contact"* — 48–445 µs per timestep across scenarios, CUDA on a GTX 2080 Ti.

**Discount that honestly at registration.** Their 30× is GPU-parallel against a deep BVH over 129k triangles; this crate's baseline is a much smaller per-chunk `parry3d` `TriMesh` on CPU. Single-digit× is the prediction, not 30×.

**C1.** Capsule-vs-field by GSS costs under **20 µs** per query on the Zen 3, and beats the shipped `TriMesh` path by at least **3×** on `fbm_terrain` at 33³ chunks. The share it moves is the 45%. *Falsified by:* under 3×, which would mean the `TriMesh` path's cost is collider *construction* rather than query and this attacks the wrong half — in which case `P-85` is the row that matters.

**C2 (determinism, which is the reason GSS and not gradient descent).** A fixed iteration count with no data-dependent branching gives bit-identical contact points across the M5 and the Zen 3, on 10⁶ queries. *Falsified by:* one differing contact, which would put the method behind the same wall as every other float-sensitive path and needs the golden-hash treatment before it can ship.

**C3 (the class it deletes rather than mitigates).** Ghost contacts — collisions against internal edges between adjacent triangles — are **structurally absent** from the field query, because there are no internal edges. Measured: non-zero ghost-contact count on the `TriMesh` path over a 495-seam-crossing walk (`M-106`'s fixture), exactly zero on the field path. *Falsified by:* zero on both, which would mean `avian3d` already handles it and Jolt's v5.0.0 internal-edge-removal work and `avian#612` are solving a problem this crate does not have.

#### P-82 — tunnelling through a wall you just made thin

**A digger is the fastest thing in the game and the wall in front of it is getting thinner every frame.** That is the exact configuration discrete collision cannot survive, and there is now a CPU-only measured answer: Pelletier-Guénette, Mercier-Aubin & Andrews, `10.1145/3747862` (SCA 2025) — spatio-temporal local optimisation for first time of impact, a modified Frank–Wolfe with golden-section search over barycentric coordinates **and time**, with adaptive triangle subdivision giving multiple contacts per triangle. **Intel i9, single CPU thread, CPU-only: an 888-triangle shuriken at 0.4 ms total, 11.13 µs mean per triangle; 100K triangles at 0.96 µs/tri.** It beats Macklin 2020's discrete detection on all their tests — self-adjudicated, since no standard triangle-SDF CCD benchmark exists.

**C1.** At the speed a `game_dig` projectile travels, the discrete path tunnels through a wall of thickness `t < v·Δt` and the CCD path does not, over 10⁴ randomised shots at a wall swept from 2 cells down to `t/h = 0.05` (`subgrid`'s floor, `M-95`). *Falsified by:* the discrete path not tunnelling, which would mean the game's speeds are below the threshold and this is premature.

**C2.** Cost per moving element stays under **25 µs** and is linear in element count. *Falsified by:* superlinear, or above 25 µs — the share is a single dynamic body's collision budget, not the whole frame.

**C3.** A capsule approximated by 8 swept spheres reaches the same first-time-of-impact as a 200-triangle capsule mesh, to within one cell, at **10×** less cost. *Falsified by:* a disagreement above one cell, which would mean the proxy is too coarse and the cheap version of this does not exist.

#### P-83 — mass, centre of mass and inertia from the triangles already emitted

**Nothing in the literature computes mass properties from an SDF analytically, and nothing benchmarks it.** The nearest usable result is Hartmann & Ewougsi Tekeu, `10.1007/s00707-025-04419-1`: inserting **T** = **x** ⊗ **x** into the divergence theorem converts the volume integral for the inertia tensor into a **pure surface integral**, evaluated analytically per triangle with linear shape functions. No wall-clock timings, but the operation count is below tetrahedral volume discretisation at equal resolution. This crate emits the surface; the integral is one pass over data it just produced, with no mesh round-trip to `parry3d`.

It is also the only row in this group that is **core-crate eligible** — pure arithmetic, `no_std`, allocation-free, generic over `Real`.

**C1.** Volume, centre of mass and the inertia tensor agree with a dense voxel reference to **1e-4 relative** on all eight reference fields at 33³, and the error falls at `h²`. *Falsified by:* worse than 1e-4, or a convergence order below 1.5 — the second would mean the surface integral is not seeing the geometry the volume integral sees.

**C2 (determinism, and the trap the paper names).** Float summation is non-associative, so a fixed triangle iteration order is required for hash stability; and the paper's own discretisation yields a **non-symmetric** inertia tensor despite the continuous form being symmetric, which they fix by averaging off-diagonals. Both are measured: fixed order gives bit-identical tensors across machines, and the pre-symmetrisation asymmetry is non-zero. *Falsified by:* zero asymmetry, which would mean the crate's triangles are better-behaved than the paper's and the symmetrisation step can be dropped.

**C3.** It costs under **2%** of extraction — the share is the whole mass-properties stage, currently paid by `parry3d` after a mesh handoff. *Falsified by:* above 2%.

#### P-84 — convexity preserved rather than recovered, against `M-116`'s 241 ms

**`A-027` already proposes this and does not cite the paper that measured it.** Müller, Chentanez & Kim, `10.1145/2461912.2461934` (TOG 2013): split the geometry into non-overlapping convex regions *offline*, then at runtime align a convex fracture pattern to the impact and intersect it against the precomputed compound. The invariant is that **clipping a convex shape against a convex cell yields convex pieces**, so runtime decomposition never happens. Measured: *"the time to fracture small to average sized objects is typically negligible, i.e. below 10 ms"*, staying below 50 ms throughout, reaching 20k compounds / 32k convexes — on a **Core i7 @ 3.07 GHz with a GTX 680**.

**And the 2026 state of the art has not improved on it.** `M-297` said no published convex decomposition runs at interactive rates; that holds. VisACD (`arXiv:2604.04244`) is a **GPU** method averaging **16.97 seconds per model** against CoACD's 36.31 s. Four orders of magnitude off a frame budget. Stop tracking that literature.

**C1.** A brush intersected against a convex-cell partition of a chunk produces convex fragments in under **10 ms** per fragment, against `M-116`'s 241–272 ms. *Falsified by:* above 25 ms, which is still 10× better and would still be worth landing — the clause is set at 10 to make the 2013 number the bar rather than a comfortable one.

**C2 (the multiplayer clause, and the one that makes this crate's version better than Müller's).** The convex-cell partition is a **pure function of the edit log**: eight same-kind brushes in all 40,320 orderings give one partition, bit-for-bit, reproducing `M-36`'s result at the cell level. *Falsified by:* more than one partition, which would mean the decomposition introduces order-dependence the field does not have and breaks the coordination-free story `M-36` bought.

**C3.** Piece count stays bounded: a chunk that has taken 60 brushes (`M-50`'s largest bucket) has under **4×** the convex cells of an unedited one. *Falsified by:* above 4×, which is the failure mode that turns this into a memory problem instead of a time problem.

#### P-85 — point the attribution instrument at the 45%

**`✗52` is the template and it is the most transferable thing Phase 23 produced.** `P-71` predicted the GPU synchronisation would be dominated by the count wait, measured it, and found the geometry copy at 60% of the whole extraction instead. **The collider stage has never had that treatment.** It is 45% of the pipeline, larger than contouring, and nobody knows whether it is BVH construction, triangle copying, `parry3d`'s constructor, or the weld that `M-69` and `✗18` argue about.

This is the cheapest row in the document and it decides which of `P-81` and `P-84` is worth doing.

**C1.** Decomposed into copy / construct / BVH-build / handoff, **one stage is over 50%** of the collider cost at 33³ and 65³ chunks on `fbm_terrain` and `gyroid`. *Falsified by:* four stages each under 50%, which would mean the cost is diffuse and there is no single lever — a null worth having, and it would redirect the group toward `P-81`'s query-side attack rather than a construction-side one.

**C2.** The dominant stage is **not** the one the docs assume. Godot's docs blame BVH/octree construction; predicted here: it is the triangle copy, by analogy with `✗52`'s finding that the copy dominated on the GPU side for the same structural reason — the bytes have to move. *Falsified by:* BVH construction dominating, which vindicates the folklore and is equally useful.

**C3.** The cost is **not** proportional to triangle count alone: at fixed triangle count, `gyroid`'s collider costs at least 1.5× `sphere`'s, because seam boundary edges (`M-69`'s 72 per seam) and degenerate slivers are per-field. *Falsified by:* proportionality to triangle count, which would make the whole stage predictable from a number the extractor already reports and is the best possible outcome for the scheduler.

#### P-86 — how many of a character's stops are geometry artefacts

**`M-115` found that a moving body is stopped harder and more often by ordinary terrain than by a chunk join, and left the harder question open: how many of those stops are *real*.** The crate records degenerate near-zero-area triangles as a metric rather than a gate, because Marching Cubes genuinely emits slivers whenever a grid corner sits near zero — that is the algorithm, not a bug. But a sliver is exactly what a capsule controller catches on, and `M-185` already found that completing the crossing identity turned a sliver into a repeated-index triangle the extractor now declines to emit. Nobody has connected the recorded metric to the gameplay symptom.

This is also where the literature is thinnest: **there is essentially no peer-reviewed work on capsule-vs-isosurface locomotion at chunk seams**, and this repository's own `M-115` is better evidence than anything published.

**C1.** Over `game_capsule_walk`'s 495 seam crossings plus a 10⁴-step randomised walk on `fbm_terrain`, at least **20%** of controller stops occur on a triangle in the bottom decile of the aspect-ratio distribution. *Falsified by:* under 20%, which would mean slivers are not what stops a character and the recorded metric has no gameplay consequence — a genuine null and worth the afternoon.

**C2.** Stops are **not** concentrated at chunk seams, confirming `M-115` from the controller's side rather than the body's: seam-adjacent triangles carry no more stops per triangle than interior ones. *Falsified by:* a seam excess, which reopens `M-133`'s "not reliably seam-closing" as a gameplay defect rather than a topology one.

**C3.** The field query path from `P-81` removes **all** of C1's stops, because it never sees a triangle. *Falsified by:* any surviving stop, which would mean the terrain genuinely stops the character there and C1 was measuring geometry that is correct.

### Group C — navigation

#### P-87 — an octree that repairs itself locally, against a Morse–Smale route this project already costed at 0.36 seconds

**The navigation direction in this repository is costed and the cost is fatal.** The docs propose deriving the encounter graph from field topology via Morse–Smale segmentation, and record the price: **PLMSS at 256³ takes 4.40 s single-threaded, 0.36 s on 24 threads**, which the docs themselves call *"20–40× over a 16 ms frame"*. Recast is worse in a different way — a 2026 devlog measures Unreal's navmesh generation at **a constant ~5 ms and a ~10 FPS drop for a relatively simple sublevel**, and Recast rebuilds a tile by voxelising collision geometry, so this game would voxelise a mesh it generated from voxels.

**There is a third option nobody in these docs has considered, and it is measured.** Massonnat & Verbrugge, `10.1109/CoG60054.2024.10645669` (IEEE CoG 2024): an octree splits cells containing obstacles, adjacent cells merge by a Hertel–Mehlhorn-inspired greedy method **while preserving convexity**, and A* runs on the resulting coarse graph with visibility pruning and a 3-D funnel. Extended to dynamic environments with **local** octree and graph updates plus a cell-repairing strategy. Measured on an **Intel Core i7-12700H laptop: octree update 0.22–1.36 ms, local graph update 0.03 ms, ≈1 ms total**, with cell reduction up to an order of magnitude (28,190 → 303).

The convexity-preserving merge is the same trick as `P-84`'s fracture invariant. Two subfields arriving at it independently is a signal.

**C1.** Built on the sign field and the active-cell bitmap — never on triangles — the local repair after `P-72`'s eleven-edit dig trace costs under **2 ms** at 4³ chunk granularity. *Falsified by:* above 2 ms, which puts it in Recast's band and removes the reason to prefer it.

**C2.** The repair set is **edit-proportional**, tracking `M-311`'s 925 dirty cells with a constant factor under 3, and does not grow with world size at fixed edit size. *Falsified by:* world-proportional growth — the same falsifier as `P-78`'s C1 and for the same reason.

**C3.** Convexity-preserving merge gives at least **5×** cell reduction on `fbm_terrain` and `gyroid`, and `gyroid` — which is topologically the hardest of the eight, tunnelling in all three axes by construction — is the one where it does worst. *Falsified by:* under 5× on either, or `gyroid` not being the worst, which would mean the reduction tracks something other than topological complexity.

#### P-88 — clearance for free from the octree, without maintaining a medial axis

**The highest-value gameplay feature in the mechanics dossier is blocked on a sub-problem the dossier itself names as the hard one.** CALIBRE — every creature carries a half-width λ, reachable space is the connected component of `{r ≥ λ}` — needs `ρ` maintained under material removal, and the docs say plainly: *"the genuinely hard sub-problem is maintaining `ρ` itself under material removal"*, with the defining paper (Chazal & Lieutier `10.1016/j.gmod.2005.01.002`) unobtainable and `O-19` still open on which direction the λ-filter runs.

**An octree of free cells gives a clearance lower bound with no medial axis at all**: a cell of side `s` that is entirely empty admits a sphere of radius `s/2`, and the merged convex regions of `P-87` give a larger bound directly. `M-346` already established the crate can state a clearance rather than only a connection, with exactly zero error where the answer is known. This asks whether the octree route reaches the same number cheaply enough to ship.

**C1.** The octree clearance lower bound is within **1 voxel** of the true clearance on ≥90% of sampled points across `M-346`'s fixtures — the same bar the dossier set for CALIBRE's oracle check. *Falsified by:* under 90%, which means the bound is too loose to gate a passage and CALIBRE still needs `ρ`.

**C2.** λ-membership flips after an edit are **at most 4×** the changed samples — again the dossier's own registered bar — and the flip set is computable from `P-87`'s repair set with no extra traversal. *Falsified by:* above 4×, or a flip set that is not a subset of the repair set.

**C3.** It is conservative in the safe direction: the bound **never overstates** clearance, over 10⁶ samples. A creature is never told it fits where it does not. *Falsified by:* one overstatement — this is one-sided and a single violation kills it, because a gameplay gate that lies in the permissive direction is worse than no gate.

### Group D — streaming, LOD and memory

#### P-89 — below 2³

**`M-377`'s own open question, and it is an hour of work.** The granularity curve is still rising at 2³ on both fields, so 1³ should be worse — the sample-duplication penalty `((c+1)/c)³` is **8.0 at c = 1** against 3.375 at c = 2 — but it is untested, and `M-377` is the entry that established a boundary minimum is not an answer.

**C1.** 1³ is worse than 2³ on both `gyroid` and `fbm_terrain`, confirming the interior optimum at 4³. *Falsified by:* 1³ winning, which would mean the duplication model is wrong and the whole granularity result needs re-deriving.

**C2.** The measured cost at 1³ matches the two-term model `M-377` derived (which gets a stencil of exactly 6 on all twelve arms) to within **10%**. *Falsified by:* above 10%, which means the model does not extrapolate and the 4³ optimum is a coincidence of the sampled points.

#### P-90 — per-block edit-list culling, at a block of 64 cells

**Dreams' hierarchical evaluator culls over 99% of naive edit evaluations by maintaining a per-block edit list** (Evans, SIGGRAPH 2015 Advances — 1 to 100,000 edits per sculpture, 10–100 M voxels evaluated per second on a PS4, culling efficiency over 99% against brute force). This crate's equivalent is `P-39`'s Lipschitz brush pruning at **3.36× median**, measured at chunk granularity. `M-377` has now moved the optimum chunk to **4³ — sixty-four cells**. The question the two results create together is whether pruning still pays when the brick is that small, or whether the per-brick bookkeeping eats it.

**C1.** At 4³ granularity, Lipschitz pruning still gives at least **2×** on `M-50`'s 46–60 brush bucket. *Falsified by:* under 2×, which would mean `P-72`'s optimum and `P-39`'s win are in tension and the granularity choice has a hidden cost.

**C2.** Carrying a per-brick surviving-brush list costs under **8 bytes per brick** amortised across the dig trace, using `✗41`'s finding that 1,507 survivors cut to 73 necessary. *Falsified by:* above 8 bytes, which at 4³ granularity is a real memory line item — a 128³ world is 32,768 bricks.

**C3.** The two culls compose rather than overlap: pruning at chunk granularity then at brick granularity removes strictly more than either alone, on the eleven-edit trace. *Falsified by:* no additional removal, which would mean the chunk-level cull already found everything and the brick level is bookkeeping for nothing.

#### P-91 — geomorph against dither, in units of pixels

**Two rows of the existing docs meet here and neither names the head-to-head.** Item 3.6 proposes Lengyel's one-scalar geomorph at 6 bytes/vertex instead of 12, and predicts it **fails on `gyroid`** and succeeds on `fbm_terrain` — "geomorph on terrain, not on caves". Item 3.8 proposes switching LOD at the distance where p99 pixels-of-pop falls under one pixel, against `M-121`'s measured pop of up to **3.14 cells**. The comparison nobody has run is geomorph against the other standard answer, dithered or alpha-to-coverage cross-fade.

**The sweep found no benchmark for any dithered LOD blend, anywhere.** The one measured argument in the literature favours morphing, and it is a *bandwidth* argument rather than an aesthetic one: Haydel, Yuksel & Seiler (`10.1145/3618359`) chose morphing over stochastic LOD because stochastic causes *"significantly increased data movement"* spikes at transitions — though their hardware is a **cycle-accurate simulation of TRaX, not a real GPU**, which is why the argument that survives the port is the structural one and not the number.

**C1.** Geomorph gets p99 pixels-of-pop under **1.0** on `fbm_terrain` at the current switch distance and fails on `gyroid`, reproducing 3.6's prediction. *Falsified by:* either half — geomorph failing on terrain would kill it outright, and succeeding on `gyroid` would be a better result than predicted and needs explaining.

**C2.** Dithered cross-fade gets a *lower* p99 pop than geomorph on `gyroid` — because it does not need a coarse-mesh counterpart along `p + t·n`, which is exactly what a triply-periodic field denies it — and a **higher** cost, measured as the temporal-history rejection of `P-77`'s instrument. *Falsified by:* dither also failing on `gyroid`, which would mean neither transition works on cave topology and the honest answer is a hard switch at a distance where nobody looks.

**C3 (the clause that connects the two).** Dither's cost lands in the same budget destruction has already spent: rejection under dithered LOD plus active digging is **more than the sum** of each alone. *Falsified by:* additivity, which would mean the two effects are independent and dither is affordable after all.

#### P-92 — regenerate, do not transmit

**The bar is public and low.** meshoptimizer decodes sponza — 184k vertices, 262k triangles — in **1.92 ms on a Core i7-8650U @ 2 GHz**, against Draco's 169 ms. That is **≈7.3 ns per triangle, single-threaded, to decode**, and decoding produces geometry and nothing else. This crate's extraction produces the geometry *and* the material weights *and* the gradient normals from a field plus an edit log that is far smaller than any encoding of 440k triangles. **And `M-31`'s bit-exact cross-platform determinism turns that from a gamble into a protocol: both endpoints extract the same mesh from the same field, and the golden hashes prove it.** The sweep found nothing published on this.

**C1.** Marginal extraction cost per triangle is under **7.3 ns** on the Zen 3 for `marching_cubes` at 33³ chunks, i.e. re-extracting is cheaper than decoding an encoded mesh of the same output. *Falsified by:* above 7.3 ns — `M-20`'s 4.75 ns/sample and `✗51`'s 10.68 ns/sample extraction marginal make this genuinely uncertain, which is why it is worth measuring rather than asserting.

**C2.** The field plus edit log for a dug chunk is at least **20×** smaller than meshopt's encoding of that chunk's triangles, across `M-50`'s four log buckets. *Falsified by:* under 20×, which would mean long edit logs approach the size of the geometry they produce and the argument has a crossover the log-trimming ticket must respect.

**C3.** Extraction from the same field and log on the M5 and the Zen 3 produces byte-identical output over a 10⁴-edit trace — `M-31` at 216 fixtures extended to a long replay, which is a different regime. *Falsified by:* one differing byte, which is the single most important negative result available here and would end the direction.

#### P-93 — what the 4³ optimum costs the GPU upload

**`M-377` named this cost and did not price it.** At 4³ the trace carries **51.6% more vertex data** than 64³ for the same 53,110 distinct surface points — 81,548 against 53,788 — with the duplication ratio running 2.2122× at 2³ down to 1.0128× at 64³. The entry says plainly that this is *"the number a consumer paying for GPU upload should weigh against the 51× edit win"*, and `✗54` has since measured what upload actually costs: **7.3236 ms of an 8.3694 ms GPU total, 87.50%**.

So the weigh-off is now computable and it has never been computed.

**C1.** There is a crossover: below some edits-per-second rate, 64³ wins on total frame cost despite the 51× edit penalty, because the upload dominates. Predicted to exist between 1 and 20 edits/second — `game_dig` throttles the stroke to **12.5/s**, which is inside that window. *Falsified by:* no crossover in 0.1–100 edits/s, which would mean 4³ wins unconditionally and the granularity decision is settled rather than a trade.

**C2.** Welding closes the duplication (`M-377` says it does) at a cost that is itself granularity-dependent, and at 4³ the weld costs more than the upload it saves. *Falsified by:* the weld being cheaper, which resolves the trade in 4³'s favour and makes C1's crossover a curiosity.

**C3.** A per-chunk granularity — coarse for static chunks, fine for actively-dug ones — beats both fixed choices by at least **1.5×** on the eleven-edit trace with a moving camera. *Falsified by:* under 1.5×, in which case one global granularity is the right answer and the scheduler stays simple.

### Group E — the edit log, which is the save file, the undo stack and the wire format

**These four measure the edit log, not a network layer.** `BACKLOG.md` closes networked and concurrent editing out and that is not being reopened. But the log is the crate's own data structure, it is what a save file is, it is what undo rewinds, and its properties have never been measured. The one public datum for a destructible world anywhere is Teardown's **~1 Mbit per client**, from a blog post — and Teardown replicates *commands*, not voxels, for the reason Gustafsson states directly: *"commands are the same regardless of object size."*

#### P-94 — how big is a dig, in bytes

**No peer-reviewed paper measures bandwidth or storage for a destructible world.** Voxel Plugin and Godot Voxel Tools both document edit-replication designs with no numbers. The 08-11 doc estimates 100k operations × 48 B = **4.80 MB** against 1.7 GB per-voxel, from CALM; that is an estimate, not a measurement, and the crate now has the machinery to make it one.

**C1.** A one-hour `game_dig` session at the throttled 12.5 strokes/second produces an edit log under **2 MB** uncompressed, and the per-edit cost is **constant in world size and in chunk granularity** — the property that makes the log the right representation. *Falsified by:* either — size growth with world size would mean the log is carrying position precision it does not need.

**C2.** Coaxial and nested edits collapse: a tunnel dug as 200 overlapping capsules along one axis compresses under `max` to fewer than **40** surviving brushes, reproducing `✗41`'s survivor arithmetic (1,507 → 73) on a realistic trace rather than a fixture. *Falsified by:* fewer than 3× collapse, which would mean the log grows with play time in a way trimming cannot fix.

**C3.** Entropy-coding the log beats a general-purpose compressor by at least **2×**, because brush parameters are strongly correlated along a stroke. *Falsified by:* under 2×, which means `zstd` is the answer and no bespoke format is worth writing.

#### P-95 — undo, and the checkpoint cadence nobody has measured

**The 08-11 doc names this gap in one line: "out-of-order arrival forces a re-fold from checkpoint; needs checkpoint cadence measured."** It is also the feature a volumetric editor is unusable without, and `game_editor` exists.

**C1.** Undoing the last edit costs strictly less than re-folding from a checkpoint whenever fewer than **N** edits separate them, and `N` is measurable and stable across fields. *Falsified by:* no crossover, which would mean undo is always a re-fold and the checkpoint cadence is the only knob.

**C2.** `M-50`'s cost curve gives the cadence directly: checkpoints every **k** edits hold worst-case undo under 16 ms, and `k` falls out of the 0.158 / 0.354 / 0.525 / 0.589 ms bucket measurements rather than being tuned. *Falsified by:* a `k` that does not predict from `M-50` — i.e. undo cost is not the same function as evaluation cost, which would be a finding about the fold rather than about undo.

**C3.** Undo is **bit-exact**: undoing then redoing an edit returns the golden hash to its previous value, over a 10³-edit trace with interleaved undo. *Falsified by:* one mismatch, which would mean the log is not a faithful history and the save format is lossy.

#### P-96 — how far apart are smooth union's 40,317 answers

**`M-38` is the measurement that decides whether smooth blending can ever be order-free, and it only counted.** Smooth union gives **40,317 distinct results from 40,320 orderings** — against hard union's 1 (`M-36`) and mixed add/subtract's 11 (`M-37`). Counting distinct results says the operator is order-dependent. It does not say the results are *different in any way a player could see*, and that is a completely different question with a completely different consequence: if the spread is sub-voxel, then smooth union is effectively commutative at gameplay tolerance even though it is not bit-commutative, and the ordering authority that `M-38` seemed to demand is not needed.

This is a one-day experiment on a fixture that already exists and it changes a protocol decision.

**C1.** Across all 40,320 orderings of the eight-brush fixture, the maximum symmetric Hausdorff distance between any two resulting meshes is under **0.1 cells**. *Falsified by:* above 0.1 cells, which confirms `M-38`'s implication and means smooth-union edits carry a hard ordering requirement.

**C2.** The spread scales with the blend radius `k` and vanishes as `k → 0`, recovering `M-36`'s exact commutativity in the limit — so `k` is the knob, and there is a `k` below which smooth union is gameplay-commutative. *Falsified by:* a spread that does not fall with `k`, which would mean the order-dependence is not the blend and is something else entirely.

**C3.** The spread is **spatially confined** to within `10k` of a seam between two brushes, matching the `|smin_k − min| ≤ k·ln n` bound the 2026-08-23 memo derived and measured at max deviation 0.135 against a bound of 0.208. *Falsified by:* spread outside that shell, which would mean the bound does not localise and the memo's `k`-as-a-length interpretation is wrong.

#### P-97 — determinism at replay scale, not at fixture scale

**`M-31` is 216 golden hashes on eight reference fields. A save file is a hundred thousand edits.** These are different regimes and nothing has tested the second. Teardown's team is the only public account, and their arc is a warning: they **first rewrote destruction in fixed-point integer arithmetic**, then found floating point workable with precautions — *"floating point operations were considered unsafe for deterministic purposes. That is still true to some extent, but the picture is more nuanced."* They were not committing cross-platform hashes. This crate is, so its bar is higher and its evidence should be too.

**C1.** A 10⁵-edit trace replayed on the M5 and the Zen 3 produces byte-identical meshes on all eight fields. *Falsified by:* one differing byte — and the value of this experiment is entirely in that outcome, because a divergence at edit 60,000 that fixture testing cannot reach is the single worst bug this crate could ship.

**C2.** If C1 fails, the divergence is **localised**: the first differing edit is identifiable by bisection and its brush parameters name the operation responsible. *Falsified by:* a divergence that bisection cannot localise, which would mean the fold accumulates rather than diverges and no per-operation fix exists.

**C3.** Replay cost is linear in log length with a constant under **1.2×** the sum of per-edit costs `M-50` measured — i.e. replay does not have a hidden superlinearity. *Falsified by:* above 1.2×, which would put a ceiling on session length independent of correctness.

### Group F — Phase 23's own leftovers

Five follow-ons the entries named and did not file, plus the audit item that was downgraded rather than fixed. All five are cheap and four of them are the direct continuation of an experiment that has already run.

#### P-98 — the fused certificate, at 0.0658

**`M-378` named this and explicitly refused to substitute it for the verdict.** The Plantinga–Vegter predicate costs 21% standalone, against a registered 5% ceiling — but the **bare eight-corner gather is 0.51 ms on every single arm**, and two thirds of the standalone cost is re-reading corners the extractor has already read. A version fused into the extractor's existing gather would pay **0.0658**, which is still above 5% but by 1.3 points rather than 16. The entry states that as a number and files nothing.

**C1.** Fused into the extraction gather, the predicate costs at most **0.0658** of extraction at 65³ on all eight fields — the figure `M-378` derived. *Falsified by:* above 0.0658, which would mean the derivation double-counted something and the fused cost has a term the decomposition missed.

**C2.** The certification result is **identical** to the standalone predicate on all 25 rows of `p-62.csv` plus the 400,000-cell random arm — same certified set, same zero unsound. *Falsified by:* any difference, which would mean fusing changed the arithmetic and the soundness proof does not carry over.

**C3.** At 0.0658 it is shippable as a default rather than a debug gate: the certified fraction is reportable per chunk at no extra pass, and a consumer can read "this chunk is 87% certified" from the existing report. *Falsified by:* needing a second pass to aggregate, which puts it back above the ceiling.

#### P-99 — the under-resolution witness, on the metric the data supports

**`✗53` falsified all three clauses and said what the next registration should be about.** Two candidate reformulations survive the data: the non-monotonic rate **among single-root edges** — which ranks `thin_plate` first and is the only quantity in the run that fails to converge — and the false-negative count at fixed `k`, which needs the oracle and is therefore offline. The entry registers neither, on purpose: *"Neither is registered, so neither is claimed."*

**C1.** The non-monotonic rate among single-root edges ranks `thin_plate` **first** of eight at 33³ and 65³, where the all-edges rate ranked it eighth. *Falsified by:* any other ranking, which would mean the denominator was not the problem and the whole witness family is measuring the volume rather than the surface.

**C2.** It **fails to converge** on `thin_plate` and converges on the other seven, reproducing `✗53`'s observation that `thin_plate`'s false-positive rate is the only one that rises with resolution (0.3889 → 0.4412 → 0.4697). That non-convergence is the signal an LOD decision wants. *Falsified by:* convergence on `thin_plate`, or non-convergence on a second field.

**C3 (the clause `✗53` got backwards, stated correctly this time).** False **negatives** fall with `k` — 322 / 94 / 8 / 1 / 0 at `k` ∈ {2,3,5,9,17} — and false positives are monotone **non-decreasing** in `k` by construction, because adding a sample point can only add an opportunity to disagree. Both are asserted in the direction the predicate's form allows. *Falsified by:* a non-monotone false-positive count, which would mean the implementation is not the predicate.

#### P-100 — is there a cell decomposition that is octahedrally invariant

**`M-372` found the obstruction and named it as structural.** Marching Tetrahedra improves to a flat **12 of 48 and stops**, because a six-tetrahedron decomposition of a cell is not octahedrally invariant — its diagonals cut different edges after a relabelling — and **no placement rule reaches that**. That is a statement about the *decomposition*, not about the extractor, and decompositions that are invariant exist: the barycentric (Kuhn/Freudenthal) subdivision into 24 tetrahedra is invariant under the full octahedral group by construction.

`P-3` already verified that Kuhn tiles face-to-face only if every cell picks the same main diagonal, and `M-51` measured classic MT at ~3× the triangles for ~4% worse geometry. So the cost side is known and the question is whether invariance is worth paying for.

**C1.** A 24-tetrahedron barycentric decomposition reaches **48 of 48** on all eight fields at both resolutions, where the six-tet split reaches 12. *Falsified by:* under 48, which would mean the obstruction is not the decomposition and `M-372`'s diagnosis needs revising.

**C2.** It costs at most **2×** the six-tet split in triangle count — the share is `M-51`'s already-measured 3× penalty over Marching Cubes, so the combined penalty must stay under 6×. *Falsified by:* above 2×, which prices invariance out for a game and leaves it a CAD-only option.

**C3.** It still tiles across chunk seams with zero open edges, which is the property `M-132` measured for the subgrid extractor and the reason Marching Tetrahedra is in this crate at all. *Falsified by:* any open edge, which kills it outright regardless of C1.

#### P-101 — can the duals be made equivariant by canonicalising the accumulation

**`M-372` says the duals accumulate crossings in an order that axis relabelling permutes, and points at `M-177`.** But `M-177` is about **negation** equivariance — *"reordering cannot buy negation equivariance, and the obstruction is structural"* — and octahedral relabelling is a different group action. The evidence that the two are not the same case is in `p-61.csv`: `dual_contouring` and `manifold_dual_contouring` already reach 48 on **3 of 16 rows**, which a structural obstruction would not permit.

So the open question is narrow and cheap: does accumulating in a **relabelling-invariant order** — sorting the crossings by a key that is a function of the cell's geometry rather than of the axis indices — take the duals from 3 of 16 to 16 of 16?

**C1.** With crossings accumulated in an order keyed on `(|value|, |offset|)` rather than on edge index, `dual_contouring` reaches 48 on at least **12 of 16** rows. *Falsified by:* under 12, which would mean `M-177`'s obstruction does cover this case and the duals are structurally out of reach.

**C2.** It is free: no golden hash moves except through the accumulation order itself, and the Hausdorff error is unchanged to 1e-12 — the same shape as `✗49`'s C3. *Falsified by:* any geometric change, which would mean the reordering is not a reordering.

**C3.** `manifold_dual_contouring` does **worse** than `dual_contouring`, because its cycle partition introduces a second order-dependence the accumulation key cannot reach. *Falsified by:* the two matching, which would be a stronger result and would mean the cycle partition is already invariant.

#### P-102 — ✗43's rate, this time with an artefact

**The one audit item that was downgraded rather than fixed.** `✗43` claims "2 of 8,064 not closed before the fix, 0 of 8,064 after" and the sweep that produced it exists in no bench, no test and no CSV. The entry now says so honestly, which leaves the claim that the per-ring apex fix **generalises** with no evidence behind it, and no ticket to produce any. `P-63` has since built exactly the machinery this needs — an exhaustive sign sweep with a fixture-can-fail control and a magnitude-seed protocol.

**C1.** Over all 8,064 configurations (16 × 8 × 9 × 7 sizes from 6³ to 12³), the post-fix extractor produces **zero** unclosed meshes, on four independent magnitude seeds reported per seed and never pooled — `M-374`'s protocol. *Falsified by:* any unclosed mesh, which reopens `✗43`.

**C2 (the control the original lacked).** The pre-fix single-apex fan produces exactly **2** at 6³ and zero elsewhere, reproducing the claim the entry could not support. *Falsified by:* a different count, which would mean the original figure was wrong as well as unreproducible.

**C3.** The `±1` magnitude arm is a bad fixture here for the reason `M-374` found — `has_inner_hexagon`'s strict `0 < x < 1` rejects, so interior vertices go to zero — and the harness reports it as VOID rather than as a pass. *Falsified by:* the unit arm reporting a pass, which is `M-44` in a new place.

---

## Part 4 — Foreclosed, and why

Seven new rows. Same framing as always: each looked transferable under its lens and is not, and recording it stops the next sweep paying for it again.

**1. Cluster LOD / Nanite-style virtual geometry, and the reason is now a number rather than an argument.** `BACKLOG.md` already rejects it on "no local validity certificate", and the measurement supports that harder than the argument did. Nothing published builds cluster hierarchies for meshes that change every frame — every paper from Nanite (SIGGRAPH 2021 Advances, course notes, no DOI) through V3DG (`10.1145/3721238.3730602`) splits into an offline build and an online selection. The build cost is measured: zeux's DAG builder does **1.64 billion triangles in 2 m 35 s on a Ryzen 7950X, 16c/32t** — ≈331k triangles per second per thread. **At 440k triangles over 376 chunks, ≈1,170 triangles per chunk, that is ≈3.5 ms on one worker to re-DAG a single edited chunk**, which is 21% of a 60 fps frame, and it overstates throughput because zeux's figure comes from 30 M-triangle meshes where fixed overheads amortise. *The chunks already are the clusters and extraction already is the LOD generator.* The one paper worth reading anyway is Ladeuil, Trabucato, Vaisse & Faraj, `10.1111/cgf.70380` (CGF 2026), which attacks locked cluster boundaries with progressive boundary portions in bijection with neighbours — the seam mechanism, not the hierarchy. Read it for `P-91`, not for a DAG.

**2. Meshlet compression in a mesh shader.** All three of the strong results — `arXiv:2404.06359` (GTS-Reuse ≈16:1, 5.9 bits/triangle, 15.56 M triangles decompressed in **0.59 ms on a Radeon RX 7900 XTX**), `10.1111/cgf.15002`, and DGF `10.1145/3675383` — decode **in a mesh shader**. `V-49` establishes that WGSL mesh shaders on Metal are `unimplemented!()` in naga 29 and land in `wgpu` v30. And `P-92` is the argument that this crate should not be transmitting geometry at all.

**3. NanoVDB, GVDB, fVDB and the SVDAG family as a storage layer.** NanoVDB is disqualified by its own documentation — *"values can be modified in a NanoVDB grid, its tree topology cannot"* — and a brush that carves new empty space outside the current narrow band **is** a topology change. GVDB is CUDA-only, dead against the Metal target. fVDB (`10.1145/3658226`) is a deep-learning framework whose construction path is optimised for point-cloud ingestion. Aokana is a raymarched renderer whose DAG deduplication makes edits expensive by construction, which is the same blocker the 08-11 doc already recorded. **`P-72` took GVDB's brick-size finding without taking any of the libraries, and that was the whole transferable content.**

**4. Breaking Good and the learned-fracture family, for *carved* geometry.** The mechanics dossier already puts fracture modes in Tier 3 at 0.5–12 s per mode; the sweep confirms it from the source and adds the structural reasons. `10.1145/3549540` explicitly produces **non-convex** fragments — the paper contrasts itself favourably against Voronoi's "convex fragments with perfectly flat sides" — so every fragment lands back in `M-116`'s 241–272 ms decomposition; precompute is per *shape*, and in a live-edited world every chunk is a new shape after every brush; and even runtime projection is 10–100 ms. DeepFracture (`10.1111/cgf.70002`) is worse on every axis: **3 days of data generation and ~24 h of training per shape**, 0.27–0.31 s inference. Its useful contribution is a baseline table putting classical Voronoi pre-fracture at **31–339 ms**, which is the honest number for the technique this project would otherwise reach for. **`P-84` is the only route, and it is 2013.**

**5. Fluids, granular material and soft bodies against a changing field, as a research direction.** Everything credible is GPU-resident and global — XPBI (`10.1145/3680528.3687577`) publishes **no timings, no particle counts and no hardware**; the convex MPM/rigid coupling (`arXiv:2503.05046`) is 21.7 ms per timestep on an **RTX 4090** with a global convex solve; the SPH–SDF coupling paper (`10.3390/math14111845`) uses three analytic primitives, never isolates the coupling cost, and does not name its GPU. And the correct primitive — project particles out along `∇φ`, which the field gives free — is a handful of lines and needs no citation. **There is nothing to transfer; there is only something to write.**

**6. Recast/Detour as a navigation substrate.** It voxelises collision geometry to rebuild a tile, so this game would voxelise a mesh it generated from voxels; the agent radius is **baked**, so N agent sizes is N navmeshes; it is 2.5D; and a 2026 devlog measures a **constant ~5 ms and ~10 FPS drop on a relatively simple sublevel**, with the developer's fix being to suppress dirty-tile marking entirely. Epic's own optimisation documentation names the cost drivers and contains no empirical timings whatsoever. `P-87` is measured at ~1 ms working on the field directly and is one representation-conversion shorter.

**7. Perceptual studies as evidence for LOD popping.** The best perceptual-visibility dataset of 2024 (`10.1038/s41598-024-78254-0`, 16 observers, 720 binary masks, a validated CNN metric) tested texture resolution, shadow-map resolution and anti-aliasing complexity — **LOD change, geometry simplification and popping were explicitly out of scope**. There is no perceptual instrument to borrow. `M-121`'s 3.14 cells and the pixels-of-pop rule are this project's own instrument and remain the better evidence; `P-91` measures against them and does not attempt a study.

### One thing that is not foreclosed and is deliberately left alone

**Multi-material as a public API change.** The PhaseTree proposal (multi-material as *phases* of the field, <25% overhead) is flagged in the existing docs as a scope question rather than a ticket, and `P-75` touches it — per-vertex material weights are the shading half of the same question. `P-75` is scoped to weights carried as a vertex attribute, which is additive and breaks nothing. **Whether `Sdf` grows a phase concept is the owner's decision and this document does not make it.**

---

## Part 5 — Acquisitions

In priority order. DOIs resolved through Crossref, OpenAlex or a publisher page; the unverified ones say so. Given `M-371`, **check any of these that home-still reports present against the three discriminators before believing it**: markdown length, `chunks_indexed`, and whether `pdf_path` ends in `.pdf`.

1. **Macklin, Erleben, Müller, Chentanez, Jeschke & Corse, *Local Optimization for Robust Signed Distance Field Collision*, PACMCGIT 3(1), 2020** — `10.1145/3384538`. Open PDF at `mmacklin.com/sdfcontact.pdf`. Blocks `P-81`; also the correct citation to replace the docs' current SDF-collision reference.
2. **Müller, Chentanez & Kim, *Real Time Dynamic Fracture with Volumetric Approximate Convex Decompositions*, TOG 2013** — `10.1145/2461912.2461934`. Open PDF at `matthias-research.github.io`. Blocks `P-84` and is the missing citation on `A-027`.
3. **Massonnat & Verbrugge, *Efficient Octree-based 3D Pathfinding*, IEEE CoG 2024** — `10.1109/CoG60054.2024.10645669`. Open PDF at `sable.mcgill.ca`. Blocks `P-87` and `P-88`.
4. **Jin, Lewis & West, *A comparison of algorithms for vertex normal computation*, The Visual Computer 21(1–2), 2005** — `10.1007/s00371-004-0271-1`. Open PDF at `users.tricity.wsu.edu`. Blocks `P-73`.
5. **Pelletier-Guénette, Mercier-Aubin & Andrews, *Real-Time Triangle-SDF Continuous Collision Detection*, PACMCGIT 8(4) / SCA 2025** — `10.1145/3747862`. Open PDF at `profs.etsmtl.ca`. Blocks `P-82`.
6. **Hartmann & Ewougsi Tekeu, *Gauss divergence theorem for the calculation of the mass and area moment of inertia tensors*, Acta Mechanica 2025** — `10.1007/s00707-025-04419-1`. Blocks `P-83`.
7. **Tan, Chua, Koh & Bhojan, *RTSDF*, GRAPP 2022** — `10.5220/0010996200003124`, `arXiv:2210.06160` (open). Blocks `P-74`'s baseline.
8. **Mäkitalo, Jääskeläinen, Ikkala & Lauttia, *k-DOP Clipping*, SIGGRAPH Asia 2024 Tech Comms** — `10.1145/3681758.3697996`. Blocks `P-77`'s C2. GPU not named in the abstract; flag when citing.
9. **Ladeuil, Trabucato, Vaisse & Faraj, *Construction of clustered HLOD with As-Simplified-As-Possible boundaries*, CGF 2026** — `10.1111/cgf.70380`. Publisher 403 during the sweep; **numbers unverified**. Read for `P-91`'s boundary mechanism.
10. **Evans, *Learning from Failure*, SIGGRAPH 2015 Advances in Real-Time Rendering** — slides at `advances.realtimerendering.com`. Industry, no DOI. The reference design behind `P-90`, and the source of the >99% culling figure.

**Two industry sources with numbers, cite as testimony not measurement:** Gustafsson's *The unlikely story of Teardown Multiplayer* (`blog.voxagon.se`, 13 Mar 2026 — the ~1 Mbit/client figure and the fixed-point-then-back-to-float account, behind `P-94` and `P-97`), and Godot Voxel Tools' performance documentation (the "collider is 3–5× meshing" claim, no absolute timings, corroborating `P-85`).

**Already in the corpus and still uncited, worth pulling forward:** `10.1109/tvcg.2009.194` (Etiene et al., verifiable visualization — convergence-order verification of the extracted surface *and its normals*, which is the quantitative version of what `P-73` measures) and `10.1007/s10851-017-0769-6` (well-composedness, behind `P-41`'s census and `✗33`'s falsified repair).

---

## Part 6 — Ordering

Scored on *what a player notices* over effort, which is a different ranking from Phase 23's. The first six are all small.

| # | id | what | effort | why here |
|---|---|---|---|---|
| 1 | **P-85** | attribute the collider's 45% | **S** | `✗52`'s instrument on the largest unexamined number in the pipeline. Decides `P-81` vs `P-84` and costs a day |
| 2 | **P-89** | granularity below 2³ | **S** | `M-377`'s own open question, an hour, and it either confirms the 51× or reopens it |
| 3 | **P-73** | gradient normals over the pseudonormal | **S** | Stops a planned ~40-line change that a 140-citation paper says makes shading worse on exactly this output |
| 4 | **P-96** | how far apart smooth union's answers are | **S** | Fixture exists; changes a protocol decision `M-38` appeared to have settled |
| 5 | **P-98** | the fused certificate at 0.0658 | **S** | `M-378` computed the number and filed nothing; turns a debug gate into a shippable capability |
| 6 | **P-102** | ✗43's rate, with an artefact | **S** | The audit item that was downgraded; `P-63` already built the machinery |
| 7 | **P-93** | what 4³ costs the upload | **M** | `M-377` and `✗54` together make a computable trade nobody has computed |
| 8 | **P-74** | AO and shadows from the resident field | **M** | The rendering gap's cheapest entry; the construction cost every paper pays is already sunk here |
| 9 | **P-81** | capsule against the field | **M** | Attacks the 45% from the query side; deletes the ghost-contact class rather than mitigating it |
| 10 | **P-83** | mass properties by surface integral | **S** | Core-crate eligible, one pass over triangles already emitted, a day of work |
| 11 | **P-99** | the witness on the right metric | **M** | `✗53` said what the next registration should be and did not make it |
| 12 | **P-87** | octree navigation with local repair | **M** | ~1 ms measured on a laptop, against this project's own 0.36 s Morse–Smale route |
| 13 | **P-92** | regenerate, do not transmit | **M** | Turns `M-31`'s determinism from a regression gate into a streaming protocol; no published prior art |
| 14 | **P-75** | material weights at the vertex | **M** | `M-138`'s 2.33× is per sample; a shading pass pays it per fragment per frame |
| 15 | **P-77** | what a dig costs the temporal budget | **M** | Gates `P-76` and `P-91`; nobody has measured history rejection under destruction |
| 16 | **P-94** | how big is a dig, in bytes | **M** | The only public datum in the world is one blog post's ~1 Mbit/client |
| 17 | **P-91** | geomorph against dither, in pixels | **M** | Two existing doc items meet here and neither names the head-to-head |
| 18 | **P-100** | an octahedrally invariant decomposition | **M** | `M-372` found the obstruction and said no placement rule reaches it — this asks whether a different decomposition does |
| 19 | **P-101** | the duals, by canonical accumulation | **M** | 3 of 16 rows already reach 48, which a structural obstruction would not permit |
| 20 | **P-90** | per-brick edit-list culling at 4³ | **M** | `P-39` and `P-72` create this question together |
| 21 | **P-86** | how many stops are slivers | **M** | Connects a recorded metric to a gameplay symptom; a null closes it cheaply |
| 22 | **P-95** | undo, and the checkpoint cadence | **M** | The feature `game_editor` is unusable without |
| 23 | **P-80** | the LOD residual as a normal map | **M** | `M-72`'s aliasing exploited rather than authored around; no measured study exists |
| 24 | **P-78** | probes invalidated per dig | **M** | The first question GI asks, answerable with `M-311`'s existing instrument |
| 25 | **P-88** | clearance from the octree | **M** | Rides `P-87`; a route to the dossier's best feature that skips its hardest sub-problem |
| 26 | **P-82** | tunnelling through a thin wall | **M** | CPU-only, measured, and the exact configuration a digger creates |
| 27 | **P-97** | determinism at replay scale | **M** | The value is entirely in the failure case, and the failure case is the worst bug available |
| 28 | **P-76** | one triplanar plane instead of three | **M** | Gated on `P-77`; may be free or may be unaffordable and the answer is the same experiment |
| 29 | **P-79** | shadow pages per edit | **M** | Epic stopped caching and published no number; this is that number |
| 30 | **P-84** | convexity preserved, not recovered | **L** | The only route past `M-116`'s 241 ms, and the 2026 state of the art is 17 seconds on a GPU |

**Six are expected to return nulls, and that is registered rather than hoped:** `P-85`'s C1 if the collider cost is diffuse, `P-86`'s C1 if slivers do not stop anyone, `P-89`'s C1 (which is the point — confirming a boundary), `P-90`'s C3 if the chunk-level cull already found everything, `P-91`'s C2 if neither transition works on cave topology, and `P-101`'s C1 if `M-177`'s obstruction turns out to cover the octahedral case too. Phase 23's own most useful rows were `✗51` and `✗54`, both of which said *do not build this*, and both of which cost a day each to learn.

**One dependency worth stating:** `P-85` should run before `P-81` and `P-84`, because it decides whether the collider's cost is on the query side or the construction side, and those two rows attack different halves. Everything else in the list is independent.
