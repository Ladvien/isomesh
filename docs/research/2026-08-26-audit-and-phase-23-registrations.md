# Audit of Phase 21, and twelve registrations for Phase 23

**Written:** 2026-08-26 · **Reads:** `FINDINGS.md` (459 entries), `BACKLOG.md`, `CLAUDE.md`, the three 2026-08-23 sweep docs, the four Phase 21 CSVs under `docs/experiments/`, the home-still corpus (9,518 documents, 290,612 embedded chunks), and roughly 280 external fetches.

**What this document is for.** Two things, in this order. First, a skeptical audit of the four Phase 21 experiments and the ✗43 entry — not to re-litigate their verdicts, which mostly stand, but because five of the defects found are the `✗35` failure mode recurring and one of them means three of four Phase 21 datasets correspond to no commit in this repository. Second, and the reason the work was commissioned: **twelve new pre-registrations**, drawn from mathematics, formal logic and systems results that the 2026-08-23 sweeps did not reach.

**What this document is careful not to be.** The 2026-08-23 sweeps produced thirteen proposals (`A1`–`A13`), nine corpus finds (`B1`–`B9`), and a foreclosure list of twenty-five items. Nothing below re-proposes a foreclosed row, and where a proposal overlaps an already-open `A`/`B` row it says so and defers to the original. Part 5 adds nine new foreclosures with reasons, which is the cheap half of this result and the half most likely to save time later.

---

## Part 1 — The audit

Working tree at `cfecdde`, 496 commits, tracked files clean. Registration ordering verifies: `5872abe` ("P-57 to P-60 registered, before any harness exists") is an ancestor of all four harness commits. The protocol held. What follows is everything downstream of that.

### 1.1 The two provenance defects, which are the serious ones

**A-1. Three of the four Phase 21 CSVs name a commit that does not exist in this repository.** `docs/experiments/p-57.csv`, `p-58.csv` and `p-60.csv` all carry `# commit 0d81ae6 (WORKING TREE DIRTY)`. `git log -1 0d81ae6` returns *unknown revision*. It was rewritten away by a rebase or an amend. `p-59.csv`'s `f55416c` does resolve. **Nothing in the repository identifies the code that produced three of these four datasets.** The stamp exists precisely so that this cannot happen, and it did not fail — the rebase happened after it was written, which is a hole in the mechanism rather than in its use.

The fix is small and belongs in `preflight.sh`: a gate that every `# commit <sha>` line in `docs/experiments/*.csv` resolves. It would have caught this the next time preflight ran.

**A-2. All four runs were made against a dirty working tree.** Every provenance line says `WORKING TREE DIRTY`. Combined with A-1, the traceability guarantee these headers are supposed to provide does not hold for Phase 21 at all — not for three of four, for four of four. The honest reading is that the Phase 21 numbers are reproducible in the sense that the harnesses are committed and re-runnable, and not reproducible in the sense the header claims.

**A-3. ✗43's prevalence sweep has no artefact anywhere in the repository.** `FINDINGS.md:9640` and `:9652` report *"2 of 8,064 not closed, both at 6³"* before the fix and *"0 of 8,064"* after. `grep -rn "8064\|8,064\|1,152\|1152"` over the whole tree returns nothing but unrelated hits in `p-21.csv` and `Cargo.lock`. There is no test, no bench, and no CSV. The entry's whole claim to be *"a rate and not an anecdote"* — and its only evidence that the per-ring apex fix generalises — is unreproducible. This is `✗35` again, and it is a direct breach of `CLAUDE.md` hard rule 4.

Related and smaller: `FINDINGS.md:9627` says *"Clean at every other size: 7³–18³ all closed"*, but the committed test `the_cell_that_fanned_two_contours_from_one_apex_is_now_manifold` asserts at 6³ and then loops over `[9, 13, 17, 25, 33]`. Sizes 7, 8, 10, 11, 12, 14, 15, 16 and 18 are covered by nothing.

### 1.2 Defects in the entries themselves

Ordered by how much they change what the entry says.

**A-4. ✗41 / M-358's headline `20×` is not like-for-like, and the second cause is unmeasured.** Two separate problems.

The fixture (`benches/experiment_p59.rs:386,393`) starts from a **base sphere of radius 6**, so 51 of 64 chunks have `necessary = 0`, and 44 of those are non-empty and mesh bit-identically from an empty brush list — their geometry is entirely the base field. Restricted to the 13 chunks where any brush matters at all, the ratio is **690 survivors / 73 necessary = 9.45×**, not 20.6×. Separately, P-39's `64 → 19` is a per-chunk *median* and `1507 → 73` is a world-wide *total*; chaining them as "a further 20.6× on top of" mixes two statistics, and there is no clean median counterpart because the median `necessary` is 0.

Then `FINDINGS.md:9494` states *"the second route is domination rather than distance"* as a result. The CSV has a `dominant_adds` column and it reads **0 on all 64 rows**, and its definition in the bench measures a different notion than the one the prose invokes. Nothing in the file measures the invoked one. The 216 near-surface unnecessary survivors currently have an untested explanation stated as a finding.

**A-5. ✗42 / M-359's stated intermediate ratio is algebraically wrong.** `FINDINGS.md:9540` gives `|σ(σ−1) + τ(1−τ)| / (σ(1−σ))` and says it reduces to `|σ − 2τ|/σ`. It does not: at `τ = 1/5`, `σ = 0.5` the stated expression gives `0.36` against the correct `0.20`. The error is a frame mix — the numerator's error formula `s(1−s) + τ² − τ` is in the **shifted-segment** fraction `u = σ − τ` while the denominator is in the **standard-cell** fraction `σ`, and the entry substitutes `σ` into both. Substituting `u = σ − τ` gives `(1−σ)(σ − 2τ)` and the conclusion follows.

**The conclusion is correct** — independently re-derived here, including the `σ < τ` branch `(σ + 1 − 2τ)/(1 − σ)`, both matching the entry's `τ = 1/5` specialisations exactly, and the `82%` measure figure reproducing to `0.823529`. The paper's ~8 dB also cross-checks: the RMS-over-cell gain at `τ_opt` is `√6 = 2.449`, i.e. 7.78 dB. But anyone re-deriving from line 9540 gets a different function, and in a ledger whose value is that its algebra can be re-checked, that is the defect.

**A-6. ✗42's closed form is confirmed only on the synthetic control, and it mispredicts both real non-degenerate fields.** `predicted_error_ratio` is NaN on 5 of 8 rows (`σ = 0` exactly). On the two computable non-degenerate rows: `gyroid` predicted `0.372965` against measured `1.972134` (5.3× off); `fbm_terrain` predicted `0.758948` against measured `0.204694` (3.7× off). The entry explains away `gyroid` as an inflection point and **is silent about `fbm_terrain`'s equally large disagreement**. The claim at `:9597` that the gain is "exact, bounded by 1" is contradicted by the only real measurement in the file that is not at the bisection floor.

**A-7. ✗39 / M-356's closing claim is contradicted by six rows its own table omits.** The "who survives" table at `:9301-9307` accounts for 36 of the 42 `fixture_can_fail` dual rows. The six missing rows are `surface_nets` at `elements_vertex_exact = 6, pure_permutation_exact = 6` on `sphere` 33/25, `torus` 25, `csg_difference` 33/25 and `thin_plate` 33 — all fixture-can-fail. That falsifies the consumer conclusion at `:9344` that equivariance is available *"for the six pure axis permutations on `marching_cubes` and for nothing else"*. Surface Nets gets all six pure permutations on 6 of its 14 fixture-can-fail rows.

**A-8. ✗39's `worst_component_ulp` still carries the artefact the entry says was discarded.** `:9291` attaches a range of "(0–2,432)" to the tetrahedral rows; the actual CSV values there are ~9.2e18, and **55 of 112 rows carry a value above 1e6**. The sorted-merge multiset difference evidently still pairs a `+x` leftover against a `−x` leftover once the vertex sets differ. The consumer advice at `:9345` — "expect vertex positions to move by up to `worst_component_ulp`" — is meaningless on those 55 rows.

**A-9. Two clauses in Phase 21 had no falsifier that could discriminate, which is Part 5 rule 12.** P-58's C1 registered the paper's ordering-independence claim and then tested a tie-break this crate invented; reversing `(value, linear_index)` changes `g′` itself and therefore the lower-star *partition*, not the traversal order inside a part — which is what §3.1's sentence is about. The verdict at `:9375-9384` says this correctly, but it repairs the registration after the fact rather than the registration having been right. P-58's C2 had no theorem behind it at all, which `:9414` concedes. Both produced useful nulls; neither was a valid test as registered.

**A-10. Two overclaims worth a one-line correction each.** `:9425`'s *"on 18 of 24 rows the ratio is exactly 1.0000 … the two sets are essentially disjoint"* is vacuous on 15 of those 18 rows (`ambiguous_cells = 0`); where ambiguity actually exists the ratio is 0.483–0.992. And `:9247`'s C3 verdict row names only `torus`, but `noise_cavity` also misses its clause under the consecutive reading (3.758, 2.289 — neither above 4); the prose at `:9435` says so, the table does not.

**A-11. ✗43's vertex-identification argument is incomplete.** `:9621-9623` argues *"1 of 3 coordinates on the sample lattice — an edge vertex has 2 of 3, so this is the cell-local interior vertex"*. That rules out edge vertices but not **face** vertices, which also have exactly 1 of 3 on the lattice; the pinch vertex's `y = −1.2000000000000` is exactly `−2 + 1·0.8`, i.e. on a cell-face plane. The conclusion is almost certainly right — the eight-faces-in-two-groups evidence carries it — but the stated argument does not establish it.

### 1.3 What survives the audit intact

P-58's instrument is the strongest of the four and was independently re-verified: Proposition 4's identity `2·pairs + critical_total == (2n−1)³` holds on all 24 forward rows, `χ = c₀ − c₁ + c₂ − c₃ = 1` on all 24 forward and all 24 reverse censuses, `max_lower_star_cells ∈ {8, 18, 27}` with 27 occurring on exactly the three fields with `critical_3 > 0`, and every headline ratio recomputes. Its vacuity control is real, not assumed: `tied_voxels` is enormous on 7 of 8 fields, so the 19 C1 passes are not vacuous, and the entry flags the three `fbm_terrain` rows that are.

P-59's aggregates all reproduce exactly, `necessary_only_hash_unchanged` is true on 64/64, and the `unnecessary_far_by_hi = 0` dead half-predicate is flagged by the entry itself. ✗43's mechanism, fix, seed and Euler arithmetic all verify (`χ = 27 − 72 + 48 = 3` for two closed sheets sharing a vertex; `28 − 72 + 48 = 4` after). ✗39's C3 arithmetic checks: 799 across exactly 103 of 112 rows, `worst_differing_vertices == order_sensitive_edges` on all 16 Marching Cubes rows, and `5,488 = 112 × 48 + 112`.

### 1.4 Recommended corrections

Six edits, none of which reverses a verdict:

1. A `preflight.sh` gate that every `# commit` in `docs/experiments/*.csv` resolves. **(A-1)**
2. Commit ✗43's 8,064-case sweep as a bench with a CSV, or restate the entry's claim to what the committed test actually covers. **(A-3)**
3. Restate ✗41's headline as `9.45×` on the 13 chunks where brushes matter, with the 20.6× kept as the world-wide figure it is; and either measure domination or move that sentence into an open question. **(A-4)**
4. Fix the derivation at `:9540` and add `fbm_terrain`'s 3.7× disagreement to ✗42's own text. **(A-5, A-6)**
5. Add the six missing `surface_nets` rows to ✗39's table and correct `:9344`; scope `:9291`'s ulp range to the 57 rows where it means something. **(A-7, A-8)**
6. Three one-line corrections: `:9425`'s vacuity, `:9247`'s missing `noise_cavity`, `:9621`'s face-vertex gap. **(A-10, A-11)**

---

## Part 2 — Three findings to file before anything is built

Each is verified from a primary source, each contradicts something currently written down, and each changes a decision.

### V-tier — `wgpu` 29 exposes subgroup operations as **native-only**, and the browser backend does not get them

`SUBGROUP`, `SUBGROUP_VERTEX` and `SUBGROUP_BARRIER` are all in `FeaturesWGPU`, the native half of the split `Features` struct. `FeaturesWebGPU` carries `TIMESTAMP_QUERY`, `SHADER_F16`, `CLIP_DISTANCES`, `IMMEDIATES` and `DUAL_SOURCE_BLENDING` — and **no subgroup entry**. So subgroups reach native Metal, Vulkan and DX12, and do not reach `wgpu`'s browser backend, *even though Chrome 134 ships the feature unflagged*. The WGSL builtins have existed in naga since v0.20.0.

**Consequence:** any subgroup path needs a fallback, and the nine wasm demos are on the fallback. Two further constraints that will bite: `subgroup_invocation_id` is validated only in **1-D workgroup** compute shaders (`gfx-rs/wgpu#5555`), so a `@workgroup_size(x,y,z)` dispatch with `y` or `z > 1` is a hard blocker today; and `subgroup_min_size`/`subgroup_max_size` moved `Limits → AdapterInfo` in v28, so in 29 they are read from `adapter.get_info()`.

*Verified:* `docs.rs/wgpu-types/29.0.0`, `docs.rs/wgpu/29.0.0`, the `wgpu` trunk CHANGELOG, `gfx-rs/wgpu#5555`, `developer.chrome.com/blog/new-in-webgpu-134`.

### V-tier — a correction to `V-23` / `GPU-007`: WGSL mesh shaders on Metal need passthrough MSL at `wgpu` 29, and **full WGSL support landed in `wgpu` v30.0.0**

`CLAUDE.md`'s characterisation is right for the pinned version and about to stop being right. From the v28.0.0 major-changes section, verbatim: *"They are now fully supported on Vulkan, and supported on Metal and DX12 with passthrough shaders."* From v30.0.0's Metal section, verbatim: *"Added full support for mesh shaders, including in WGSL shaders. By @inner-daemons in [#8739]."* Trunk's own `docs/api-specs/mesh_shading.md` still lists the naga MSL backend as "❌ Planned" and is **stale against its own changelog**.

**Consequence:** the mesh-shader draw is worth 6.7% at 129³ (`M-149`). Paying for that at `wgpu` 29 means hand-written MSL, a forked shader pipeline and a second source of truth, in a rule-5 environment. 6.7% does not buy that. When Bevy moves to `wgpu` 30 the same 6.7% arrives for a feature flag and a WGSL entry point. **The decision is "wait", and it should be recorded as a decision rather than an omission.** `O-5` — the Metal probe — is still unrun and is now cheap to interpret, because you know what answer `wgpu` 29 should give.

### V-tier — CGAL shipped an Isosurfacing package in 6.1, and its octree domain openly leaves holes at resolution transitions

Announced 2025-10-10 by Rouxel-Labbé, Stahl, Zint and Alliez. It ships `marching_cubes()`, `marching_cubes(use_topologically_correct_marching_cubes = true)` (which is Grosso 2016's TMC, `10.1111/cgf.12975`) and `dual_contouring()`. Only the TMC row claims 2-manifold, watertight and topologically correct, and it claims it *"as long as the isosurface does not intersect the domain boundaries"* — the same caveat this crate's `is_closed()`/`is_manifold()` split already encodes. Two things matter here.

**The architecture is worth copying and the code is not.** The domain concept splits three ways — a **partition** (`Cartesian_grid_3`, `Octree_partition`), a **value field**, and a **gradient field** as a *separately supplied capability* rather than a finite-difference assumption baked into the extractor. That is a sharper factoring than `Sdf` + a `NormalStrategy`, and it is free to adopt. The package is **GPL**, so the papers are fair game and the implementation is not.

**The seam result is the finding.** CGAL's own octree example admits: *"the surface can enter and leave a cell without involving the cell's vertex. In practice, that means a hole if at nearby adjacent cells the voxels did get refined."* A mature C++ library, shipped in 2025, does not have adaptive seam handling. Combined with a sweep that found **no peer-reviewed 2023–2026 work on chunk seams, LOD stitching or crack-free transitions for the DC/MC family at all** — the line still terminates at Schaefer/Ju/Warren 2007 and Lengyel's Transvoxel — `M-128`, `M-132`, `M-133` and `M-106` are not this crate catching up to the literature. There is no literature to catch up to. That is worth saying out loud in the README, and it is worth a paper.

*Also unbenchmarked in CGAL's manual, and worth not inheriting:* "DC generates fewer faces and higher quality faces than MC, in general" and "MC is substantially faster than Delaunay refinement" appear with no numbers anywhere in the package.

---

## Part 3 — Twelve registrations

Numbering continues from `P-60`. Each follows Phase 15's protocol: the `P-` id goes into `crates/isomesh/src/experiment.rs` **in its own commit, before any harness**, with a `falsified_by` that is required by the type; numbers land in `docs/experiments/p-n.csv` with a resolving SHA (see A-1); and an `E×n` row is owed either way.

Three groups, in the order the request named them: **mathematics and logic**, then **new capability**, then **speed**. Within each group they are ordered by measured-prior-evidence over effort, which is not the same as ordered by interest.

### Group A — mathematics and logic

#### P-61 — the crossing is stored as a signed offset from the edge midpoint, and the mesh becomes bit-exactly equivariant under all 48 octahedral elements

**The hook is ✗39.** `M-356` found that a mesh is bit-exactly equivariant under axis *relabelling* and not under *reflection*, and attributed it to `a/(a−b)` and `b/(b−a)` being two different divisions of the same values. That attribution is right and incomplete, and completing it produces a fix.

The subtraction is not the culprit: IEEE round-to-nearest is sign-symmetric, so `fl(b−a) = −fl(a−b)` **exactly**, and the two denominators are exact negations. The asymmetry enters in two places — the division rounding (`fl(a/d)` and `−fl(b/d)` are independent roundings that do not sum to 1), and the **anchor**: `cube.rs:186` returns a parameter measured *from the lower corner*, and `surface_nets.rs:119` and `trilinear.rs:1024` then place the vertex at `lower + t`. A reflection swaps which corner is "lower", so the correct reflected parameter is `1 − t`, and `1 − t` is not representable as `b/(b−a)`.

**The construction.** Store the crossing as a signed offset from the edge *midpoint* rather than a parameter from the lower corner:

```
d = ((a + b) / 2) / (a − b)          # in [−1/2, +1/2]
position = edge_midpoint + h · d
```

This is exactly antisymmetric under the simultaneous swap, and the proof is three lines: `fl(a+b) = fl(b+a)`; halving is exact; `fl(b−a) = −fl(a−b)`; and `fl(S/(−D)) = −fl(S/D)` because round-to-nearest is odd. Every step is an IEEE 754 guarantee, not an empirical observation. The `[0,1]` parameter frame cannot host this because reflection acts there as `0 ↔ 1`, an *affine* map; the centred frame makes it a sign flip, which is what floating point respects exactly.

**Pre-measured, off-repo, and to be re-measured by the harness** (2,000,000 random straddling pairs per row, `f64`, script in the appendix):

| form | reflection mismatches, cell-local | world frame, `h = 0.125` | `h = 0.1` | `h = 3/32` |
|---|---|---|---|---|
| current `L + h·(a/(a−b))` | 1,125,580 / 2,000,000 (56.3%) | 25,554 / 300,000 (8.52%) | 74,961 / 100,000 (75.0%) | 9,118 / 100,000 (9.1%) |
| proposed `mid + h·(((a+b)/2)/(a−b))` | **0 / 2,000,000** | **0 / 300,000** | **0 / 100,000** | **0 / 100,000** |

Accuracy against exact rational arithmetic, world positions, 300,000 samples: current mean **0.086** ulp / worst **422** ulp; proposed mean **0.052** ulp / worst **757** ulp. It is a wash, slightly better in the mean and slightly worse in the tail. `0 / 300,000` out-of-cell. At endpoint magnitudes of `1e−300` the current form mismatches 129,164 / 200,000 and the proposed form 0 / 200,000.

**C1.** With the centred form, the ✗39 harness reports `elements_vertex_exact = 48` on **every** `fixture_can_fail` row of both the primal and dual families, at both resolutions. *Falsified by:* any row below 48 that `fixture_can_fail` marks true.

**C2.** No reference field's golden hash is unchanged. This is a **cost clause, not a benefit clause** — the vertex positions genuinely move, `T-007`'s 216 golden hashes must be rebaselined in the same commit, and a claim that they do not move would mean the change did nothing. *Falsified by:* any field whose hash survives, which would mean the centred form is not on the path that produces vertices.

**C3.** Symmetric Hausdorff on all eight reference fields at 33³/65³ changes by less than 1% in either direction, and the self-intersection-per-1,000 counts change by less than 1%. *Falsified by:* any field moving more than 1% on either metric — in which case the ~3-ulp tail is buying a real geometric cost and the trade is a decision, not a fix.

**C4, registered as a question with a stated risk of a null.** `M-32` says chunk seams are bit-exact only when the cell size is a power of two. The mechanism `M-32` names is `world_of_sample`, not the crossing, so the honest prediction is that C4 changes nothing. But the `h = 0.1` row above moves from 75.0% to 0%, so the seam sweep is re-run at `h = 0.1` and `h = 3/32` and the answer recorded either way. *Falsified by:* nothing — this clause is registered as measurement, and a null is the expected outcome.

**Why it is worth doing even if C4 is null.** ✗39's own consumer note says equivariance is currently available for six of 48 elements on one extractor. At 48 of 48 across the family, the octahedral group becomes a **metamorphic test oracle with 48 independent relations per fixture** rather than an interesting negative result — which is the `B8`/Villar "equivariance by construction" line arriving by a different door than the vertex rule.

#### P-62 — the Plantinga–Vegter gradient predicate as a *positive* per-cell certificate, checked for soundness against the tunnel classifier

**The gap.** `P-48` gave the crate a certificate that a cell is **empty** (`M-347`: zero unsound over 1.07 × 10⁹ evaluations). `P-54` gave it a tighter one via affine arithmetic (`M-354`: 3.85× more rejections on `gyroid`). There is no certificate in the other direction: nothing in the crate can say *"this cell's surface patch has no hidden topology"*. That is the difference between a mesher that is correct and a mesher that can **state where it is correct**, and it is the CAD half of the crate's mandate.

**The construction, and it is already half-owned.** Plantinga & Vegter's `C1` predicate (`10.1145/1057432.1057465`, in corpus, already used at `T-015`) is

```
0 ∉ □f(C)   ∨   ⟨□∇f(C), □∇f(C)⟩ > 0
```

The second clause says the gradient direction varies by less than π/2 across the cell, so `f` is strictly monotone in some coordinate direction and the zero set inside `C` is a **graph over a coordinate plane** — no cavity, no tunnel, no second sheet. Cost: one interval evaluation of `∇f` over the box (which a dual method already computes for QEF normals), then three interval multiplies and two adds. Constant, branchless, no allocation, no iteration, chunk-local, `libm`-only.

**The soundness check is the part that makes this a real experiment rather than a port.** This crate already owns a ground truth the PV literature does not: `A-020`'s classifier counts tunnels and twelve-vertex contours from the trilinear itself — `M-214` recorded 2,053 and 173 in 396,000 cells, and `M-222` established that χ falls by exactly two per tunnel. A cell containing a tunnel is a cell whose patch is *not* a graph. So:

**C1 (soundness, one-sided, the kill-shot).** Over eight reference fields at 17³/33³/65³, **zero** cells are `C1`-certified while the `A-020` classifier reports a tunnel or a twelve-vertex contour in them. *Falsified by:* one such cell. A single unsound certificate kills the direction, and `M-214`'s counts prove the fixture can produce the configuration — this is not a `M-44`-style pass over an unreached case.

**C2 (yield).** The certified fraction of *surface* cells is above 50% on `sphere`, `torus` and `box_exact` at 33³, and rises monotonically with resolution on all eight fields. *Falsified by:* a fraction below 50% on any of those three, or a non-monotone sequence on any field — the latter would mean the predicate is measuring the interval arithmetic's slack rather than the field's geometry.

**C3 (cost).** The predicate costs under 5% of extraction wall time on `marching_cubes` at 65³. *Falsified by:* above 5%, in which case it is a debug gate rather than a shippable capability.

**Caveat, registered rather than discovered.** `C1` guarantees the patch is a *graph*, not that its planar domain is connected — a graph over a disconnected planar region still has several components. PV closes that globally with a *balanced* octree, which this crate does not have. So the honest claim is "no hidden topology in this cell", not "exactly one component", and the entry must say so. Lin & Yap document the same gap (`10.1007/s00454-011-9345-9`).

#### P-63 — the vertex-link question in `O-12` is finite at 2¹⁸, and one sweep settles it for Marching Cubes

**`O-12` is the oldest open question in the ledger** — *"is Marching Cubes unconditionally manifold now?"* — and its own text says what would settle it: *"an exhaustive search over configurations spanning more than two cells … or a proof that a cell-local cycle triangulation plus shared face segments cannot produce a non-manifold vertex."* ✗43 found the first counterexample and it was inside one cell. The question stands because *"a vertex whose two face groups sit in different cells would be"* a third mechanism.

**The search space is much smaller than it looks.** In this crate every Marching Cubes vertex sits on a grid **edge** (`cube.rs`'s `EDGE_CORNERS` is lower-corner-first, and the interior vertex is cell-local after ✗43's per-ring apex). Every face incident to an edge vertex therefore comes from one of the **four cells sharing that grid edge**. Those four cells span a 3 × 3 × 2 block of grid corners — **18 corners, 2¹⁸ = 262,144 sign patterns**. That is not a sample. It is the whole space, and it runs in seconds.

This is a **proof by exhaustion of the vertex-link case for Marching Cubes**, and it is the case Chernikov & Xu's Coq work does not cover: their 2013 IMR proof (`10.1007/978-3-319-02335-9_28`) enumerates all 2⁸ single-cube sign configurations "disregarding any perceived symmetry" and proves cohesion and water-tightness, then composes to a grid via face-local consistency. Face-local consistency is exactly the argument that does not reach a vertex link.

**C1.** Over all 262,144 patterns, meshing the four cells, welding, and walking the connected components of the incident-face link of the shared edge vertex yields **zero** non-manifold vertices, with the interior-ambiguity rule both off and on. *Falsified by:* one pattern that does not, which is `✗49` and the third mechanism `O-12` asks about.

**C2 (the fixture can fail).** Injecting ✗43's pre-fix single-apex fan into the same sweep produces a **non-zero** count. *Falsified by:* zero — which would mean the link walk cannot see the one defect known to exist, and the sweep proves nothing.

**C3 (the dual family is where it is interesting).** The same sweep over `surface_nets`, `dual_contouring` and `manifold_dual_contouring` produces a non-zero count, and that count is a **function of the well-composedness census** in the sense of `M-338`'s bijection — the number of link-defective patterns equals the number of critical sign configurations in the block. *Falsified by:* a count that is non-zero and not equal, which is the more interesting outcome and would mean the bijection is cell-local and does not extend to a block.

**Scope note.** C1 and C2 are complete for Marching Cubes. For the dual family a vertex lives at a cell centre and its link involves the cell's 26 neighbours, which is 4³ = 64 corners and out of reach; C3 is therefore a **necessary-condition sweep on the same 18-corner block**, and the entry must not claim more. The full dual sweep at 2²⁷ over a 3 × 3 × 3 corner block (134,217,728 patterns, an estimated 4–45 minutes single-threaded) is a **nightly** gate and a separate ticket, deliberately not registered here.

#### P-64 — bounded model checking proves the combinatorics that property tests only sample

**The split that makes this tractable.** Bit-blasting IEEE 754 to SAT is the adversarial case for a model checker, and a harness over eight nondeterministic `f32` corner values is 256 bits of unconstrained float — precisely the shape Kani is worst at. But this crate's correctness risk is not in the arithmetic; `CLAUDE.md` rule 5 names it exactly: *"wrong case tables produce meshes that look fine and are subtly non-manifold."* That is **combinatorics over eight sign bits**, which is 256 states and trivial for BMC.

So: verify the combinatorics, keep testing the arithmetic. **Kani** (`arXiv:2607.01504`) ranges over the eight sign bits and proves, for all 256 patterns, that no case-table index goes out of range, no emitted index is ≥ the vertex count, no triangle carries two equal indices, and nothing panics. **Flux** (`10.1145/3591283`) carries length refinements through `MeshBuffer` so index safety in the buffer-writing layer is a *type* property checked on every build with no harness at all — and its published result is that it cut Prusti's annotation burden from up to 24% of code size to nothing by inferring loop invariants.

Both are dev tools with no runtime footprint, so hard rule 3 is not engaged. **No published use of either on geometry or graphics code was found**, which makes this novel as well as useful.

**C1.** Kani proves all four properties over all 256 sign patterns for `marching_cubes` with the interior rule off, in under 10 minutes on the M5. *Falsified by:* a timeout, or a property that cannot be expressed against the sign abstraction — the second outcome is the more informative one and means the abstraction is wrong, not the tool.

**C2.** It finds nothing the existing suite does not already cover. This is registered as the **expected** outcome and it is still worth the run: a proof and a passing property test are different objects, and `M-208`–`M-213` is five pre-registered claims that were true on seven fields and false on the eighth. *Falsified by:* Kani finding a reachable violation, which is a `✗` entry and the most valuable outcome available here.

**C3.** Turning the interior-ambiguity rule on keeps C1 under 30 minutes. *Falsified by:* a blow-up, which localises the state explosion to the interior rule and is itself a finding about that rule's branching.

**Scope note.** Neither tool touches vertex placement, and the registration must say so: placement stays under proptest and golden hashes. The honest scope is *"the table cannot be indexed wrongly"*, not *"the mesh is correct"*.

#### P-65 — MCPro's procedural construction against the one configuration this crate refuses to mesh

**The hook is `Error::UnresolvedSixSaddle`.** The README states it plainly: a cell whose contours run past Grosso's Corollary 6 bound has no published triangulation, so `extract` returns an error rather than emitting a hole (`A-002b`, `A-020`, `M-228`). `M-231` found that the `[9,3]` cell is not a topological subcase but a **singular face** the strict interior test lets through, and `M-233` recorded the blocker: a singular face needs a third routing and the resolution mask has two.

**MCPro is the paper that says it built the third routing.** Stahl & Grosso, *MCPro: A Procedural Method for Topologically Correct Isosurface Extraction Based on Marching Cubes*, GRAPP 2025, `10.5220/0013309800003912`. **No lookup table at all**: it classifies each face's bilinear restriction via the asymptotic decider into hyperbola / singular cross / single line / degenerate plane, divides the face into quadrants around the asymptote centre, assembles face segments into a per-cell halfedge structure, solves for the interior saddles, and builds the inner hexagon that decides tunnel-versus-disk. What is new against MC33 is exactly **singular faces**, and the authors exhibit configurations with three adjacent singular faces where MC33 produces topologically incorrect non-manifold triangulations. Validation: passes all 20,000 test cases of Etiene et al. (2012) on Betti numbers and Euler characteristic. Stahl is also a co-author of the CGAL package.

**One honest disclosure from the paper that this crate must carry if it adopts the method:** a trilinear isosurface **can genuinely be non-manifold** — singular vertices and edges are real features of the interpolant — so MCPro guarantees a topologically correct triangulation *of a possibly non-manifold surface*. On a field that produces a singular face, `is_manifold()` failing is **correct behaviour**, not a defect. That is a sharper statement than CGAL's TMC row and it changes what the validity gate should assert.

**C1.** The procedural construction produces a triangulation for **every** configuration that currently returns `UnresolvedSixSaddle`, and for the `[9,3]` singular-face cells `M-231` identified. *Falsified by:* any configuration it also cannot resolve, which would move the blocker rather than remove it.

**C2.** On the eight reference fields at 33³/65³, MCPro's output has the same Euler characteristic as the current extractor's on every cell where the current extractor succeeds, and χ differs by exactly `−2` per tunnel where it does not — i.e. `M-222` survives the change. *Falsified by:* any χ disagreement outside the tunnel accounting.

**C3 (cost).** The paper reports ~10% more vertices and triangles than MC33 on a skull volume and gives **no timings at all**. This crate's `M-42` measured the asymptotic decider as free to within a few percent, and `M-223` measured the interior rule at 1.95% at 33³. Predicted: MCPro costs under 25% more than the current interior-rule path on `noise_cavity`, the only field that exercises it. *Falsified by:* above 25%, in which case it is a correctness mode rather than a default.

**Effort warning, stated at registration.** This is the largest item in this document. The determinism risk is real and specific: asymptotic-decider comparisons near zero are exactly where `f32` and `f64` diverge, and `M-43` already established that the decider needs no division and no epsilon — that property must be shown to survive the procedural construction, or the golden hashes become scalar-dependent.

### Group B — new capability

#### P-66 — the monotone-edge witness, held to a ground truth this crate already computes

**The line this replaces died twice, and the third attempt has to be a different mechanism.** `P-43`/✗29 tried one evaluation at the cell centre as an under-sampling witness and was falsified on both clauses. `P-44`/✗31 tried the mean residual instead and was falsified out of sample. Both were *value* witnesses: they asked whether the trilinear's value at a point disagrees with the field's. The failure mode they were chasing is not a value disagreement, it is a **missed root** — an edge whose two endpoints have the same sign while the field crosses zero twice between them, or opposite signs while it crosses three times.

**The witness is a derivative sign test, and it comes from a two-week-old paper.** Finken, Li, Wang, Guo & Levine, *Topology-Preserving Meshing of Implicit Scalar Fields via Monotonicity Constraints*, `arXiv:2608.12142`, IEEE Vis 2026 short paper. Their central statement: if every edge of a PL mesh is monotonic with respect to `f`, the PL approximation is topologically consistent with `f`'s critical points. The test itself is one line — sample `∇f · ê` at `max(2, ⌈‖e‖/w⌉ + 1)` points along the edge and declare it non-monotonic when any two sampled projections disagree in sign.

**Do not port the paper.** It is explicitly 2D and the authors say so; its sampling-density argument, its Theorem 1 case analysis (3D Morse theory has four critical-point types and a spherical link, not three and a circle) and its separatrix refinement all fail to generalise, and it wants a Hessian this crate's field trait does not expose. **Take the edge test alone**, as a diagnostic rather than an extraction rule.

**What makes this registrable here rather than anywhere else is the oracle.** The subgrid extractor already finds *all* roots along each tet edge — `M-94` resolves a slab at 1/1000 of the edge length, `M-168` gave each crossing an identity, and `M-169` established that identity-based sharing is complete exactly when no root lands on a grid sample point. So the number of roots per edge is a **known quantity in this repository**, and the monotone-edge test can be scored against it rather than against a hunch.

**C1 (soundness, one-sided).** On eight reference fields at 17³/33³/65³, with `k = 5` samples per edge, every edge the subgrid root finder reports with more than one root is flagged non-monotonic. **Zero false negatives.** *Falsified by:* one multi-root edge the test calls monotonic.

**C2 (yield).** The false-positive rate — edges flagged non-monotonic that carry exactly one root — is below 20% at `k = 5` on the six smooth fields, and the rate falls as `k` rises. *Falsified by:* above 20% on any of the six, or a rate that does not fall with `k`, which would mean the test is measuring sampling noise.

**C3 (it is a resolution witness, which is the point).** The per-chunk non-monotonic edge fraction falls monotonically with resolution on all eight fields, and is highest on `thin_plate` — the field whose sub-cell features Marching Cubes structurally cannot see (`M-100`). *Falsified by:* a non-monotone sequence, or `thin_plate` not ranking first.

**Consequence if it holds.** A chunk can report *"this grid under-resolves this field, here"* as a number, per chunk, cheaply — which is the missing input to an LOD decision that `M-121`'s 3.14-cell surface pop and `M-72`'s aliasing both want and neither has.

#### P-67 — reduced affine arithmetic: the same rejection rate in a fixed-size struct

**`P-54` held and left a structural question open.** `M-354` measured affine arithmetic rejecting 3.85× more cells on `gyroid` than intervals, and exactly zero more where `min`/`max` destroys the correlation. The form it used grows a noise symbol per non-affine operation, so its size is a function of the tape — which for a twelve-brush edit log means an allocation whose size depends on the scene, inside a `no_std` crate whose whole pitch is that it does not do that.

**Reduced affine arithmetic fixes the symbol count.** All condensed error is folded into a single accumulated term, so the form is a fixed-size struct — `[R; N+1]` — with no allocation and a fixed operation sequence, which is the determinism property that matters here. The published cost/tightness comparison is on exactly this workload, not on ODEs: Knoll, Hijazi, Kensler, Schott, Hansen & Hagen, *Fast Ray Tracing of Arbitrary Implicit Surfaces with Interval and Affine Arithmetic*, CGF 28(1):26–40, 2009, `10.1111/j.1467-8659.2008.01189.x`. Verbatim: *"For typical functions with fairly low-order coefficients and moderate cross-multiplication of terms, RAA is generally 1.5–2× faster than IA. For functions with high bound overestimation… RAA is frequently three to four times faster."* And the counterexample they record: *"IA remains far more efficient for superquadrics."*

**C1.** RAA with 3 retained symbols reaches at least 90% of the full affine form's cell-rejection rate on `gyroid`, `noise_cavity` and a twelve-brush CSG tape. *Falsified by:* below 90% on any of the three, which means the condensation eats the correlation the form exists to preserve.

**C2.** RAA allocates zero bytes and its operation count is independent of tape length, measured by a counter rather than asserted. *Falsified by:* any allocation, or an op count that varies with the tape.

**C3.** RAA is between 1.5× and 4× faster than the interval form on the CSG tape, reproducing Knoll's band; and it is **slower** than intervals on `box_exact`, reproducing his superquadric exception. *Falsified by:* outside the band, or `box_exact` not showing the inversion — the second is the more interesting failure, because Knoll's exception is the part of his result that is a mechanism rather than a measurement.

#### P-68 — a running error bound turns a vertex position into a vertex position *and a certified interval*

**The gap is stated in the crate's own mandate.** A CAD tool wants to know how much to trust a coordinate. `M-142` found GPU and CPU agree on every triangle and disagree on 6% of vertices by exactly one ULP; `M-144` found bit-identity is a property of the cell size, not the port; `M-30` found an unclamped solve can fling a vertex 3.18 cells out of its own cell. Every one of those is a statement about error that had to be *measured after the fact*. Nothing in the crate reports, per vertex, at run time, how wide the interval containing the true crossing is.

**The construction is two extra flops per crossing and one per solve.** A running error bound propagates a first-order error term alongside the value — the same machinery as Shewchuk's adaptive predicates (`10.1007/PL00009321`, already mined in five files) and the static/semi-static/dynamic filter hierarchy of Bartels, Fisikopoulos & Weiser (`10.1007/s10543-023-00975-x`, in corpus, uncited), but used to *report* rather than to *branch*. P-61's centred form makes this cheap in a way the parameter form does not: the offset `d` is a single quotient with no cancellation, so its error bound is `|d| · (2 + |a−b|_err/|a−b|) · u` rather than the compounded bound a subtract-then-lerp accumulates.

**C1.** The reported bound is **sound**: on eight reference fields at 33³, comparing against exact rational arithmetic on 10⁶ crossings, the true error never exceeds the reported bound. Zero violations. *Falsified by:* one violation.

**C2.** It is not vacuous: the median reported bound is under 4 ulp and the 99th percentile under 64 ulp on the six smooth fields. *Falsified by:* a median above 4 ulp, which means the bound is a formality rather than information.

**C3.** Carrying it costs under 3% of extraction wall time when enabled and exactly zero when the feature is off, verified by a golden-hash comparison against the current build. *Falsified by:* above 3%, or any hash movement with the feature off.

**C4 (the reason a game wants it too).** On `csg_difference`, the reported bound is largest exactly at the seam cells `M-350` bounded — the entry proved a central difference across a CSG seam is wrong by at most half the crease angle, and this is the same locality showing up in the position rather than the normal. *Falsified by:* no correlation between the bound and seam proximity, which would mean the bound tracks magnitude rather than conditioning.

### Group C — speed

#### P-69 — autovectorising the sample loop, with bit-identity as the gate rather than the hope

**This is the cheapest large speedup available and it needs no dependency, no `unsafe` and no new crate.** `core::simd` is nightly and is staying nightly — the `LaneCount<N>: SupportedLaneCount` bound and the mask-element-type mismatch are unresolved (`rust-lang/portable-simd#364`), and the maintainers' own 2025 summary is *"nightly-only and will remain such for the foreseeable future."* So the lever is autovectorisation, and the measured prior says autovectorisation is enough: Wilcox's AArch64/NEON study on 100k `f32` samples measured scalar 77.67 µs, hand-written intrinsics 25.78 µs, and **autovectorised safe Rust 25.54 µs** — safe code matched intrinsics.

The patterns that decide it are specific and all of them are shape, not machinery: struct-of-fields rather than index arithmetic (`dst[i*2+0]` **defeated** vectorisation where a two-field struct did not), pre-slicing once outside the loop so LLVM can prove the bound, and `chunks_exact` / `zip` iterators, which vectorise.

**The float caveat, stated precisely, because it cuts the right way here.** The blanket claim that autovectorisation fails on floats is about **reductions** — LLVM will not reassociate float adds without fast-math, which stable Rust does not expose. Elementwise float map and zip vectorise fine. `isomesh`'s field evaluation is elementwise (`sin·cos`, `sqrt`, `floor` over independent samples); its accumulations — active-cell popcounts, vertex counts — are **integer**. Both halves land on the good side of the caveat.

**C1.** Restructuring the sample loop to `chunks_exact` over `[R; 8]` with the bound hoisted gives at least **2×** on the marginal `f32` cost — `M-20`'s 4.75 ns/sample falls below 2.4 ns — on `sphere` and `gyroid` at 129³ on the M5. *Falsified by:* under 2×, in which case the ceiling is NEON's 4-wide `f32` and the honest number is smaller than the prior suggests.

**C2 (the gate, and the clause most likely to kill this).** All 216 golden hashes are **unchanged**. Vectorisation must not move one bit; IEEE elementwise operations are exact per lane, so a hash movement means LLVM reassociated something, and the change is rejected rather than rebaselined. *Falsified by:* any hash movement — and unlike P-61's C2, here a movement is a defect, not a cost.

**C3.** The `f64` gain is at most half the `f32` gain, because NEON is 2-wide at `f64`. *Falsified by:* an `f64` gain above half, which would mean the `f32` path was not the vector path and C1's number came from something else.

**Verification requirement, stated at registration.** `cargo-show-asm` output for the monomorphised `f32` instance goes in the ticket. The crate is generic over `Real`, and LLVM vectorises the monomorphised instance or does not — that must be inspected per instantiation, not assumed, and a Criterion delta alone cannot distinguish a vectorised loop from a lucky one.

#### P-70 — subgroup ballot compaction, native only, with the current scan retained

**The measured prior is on this exact hardware class and this exact API.** Smith, Levien & Owens, *Decoupled Fallback: A Portable Single-Pass GPU Scan*, SPAA '25, `10.1145/3694906.3743326`. Inclusive prefix sum over 2²⁵ elements, all measured on WebGPU/Dawn: Apple M1 Max 36.85 vs 25.75 G ele/s (**1.43×**), Apple M3 10.87 vs 7.46 (**1.46×**), RTX 2080 Super 1.49×, RX 7900 XT 1.33×, Mali-G78 1.35×, Intel HD 620 1.43×.

**Two of their findings foreclose the obvious plan and should be recorded whatever this experiment returns.** First, verbatim: *"data collected by Sorensen et al. shows that ARM and Apple GPUs do not [provide forward-progress guarantees]"* — plain decoupled look-back **times out on M1 Max and M3**. Single-pass chained scan is off the table on this crate's primary target unless it carries a fallback. Second, the ceiling is structural: reduce-then-scan moves O(3n), chained scan O(2n), so **50% is the theoretical maximum** and 30–50% is what they measure. Scan is memory-bound and subgroups do not change that.

**Which is why the target is not a faster scan.** `M-150` already took the prefix sum from 15.01 to 9.65 ms; the literature says a further scan improvement is capped near 1.5× and much of it is banked. The larger win is **not scanning globally at all**: `subgroupBallot` plus `subgroupExclusiveAdd` plus one wave-level `atomicAdd` collapses per-workgroup compaction into registers and reduces the global scan to a much shorter one over per-workgroup totals. On Apple Silicon the subgroup size is 32 with ballot capped at 64 bits (`gpuweb#3950`), so one ballot covers half a word of `P-40`'s active-cell bitmap — a clean mapping rather than a fight.

**C1.** Subgroup compaction takes the 129³ extraction below 7.0 ms on the Zen 3 / RTX 3090 rig, from `M-150`'s 9.65. *Falsified by:* above 7.0 ms, which would mean the compaction was not the residue.

**C2.** Output is bit-identical to the current path on all eight fields — same triangles, same order. *Falsified by:* any difference; a compaction that reorders is a determinism regression and is rejected regardless of speed.

**C3 (the constraint clause).** The fallback path is exercised, not merely present: the web build takes it, and a native run with `SUBGROUP` forced off produces the same output as the subgroup run. *Falsified by:* a fallback that is not tested, which is `M-44`'s zero-that-could-not-have-been-non-zero in a new place.

**Blocker to resolve before writing the shader:** naga validates `subgroup_invocation_id` only in **1-D workgroup** compute shaders. If the bitmap pass dispatches with `y` or `z > 1` it must be flattened first, and that flattening is its own change with its own hash risk.

#### P-71 — the 83% is a blocking round-trip, and both targets can avoid it

**`M-167` is the largest single number this project owns about its own GPU path:** synchronisation was 83% of an extraction. `M-159` localised it — the last four bytes cost 0.033 ms to move and 0.375 ms to wait for, because `poll(Wait)` drains every dispatch queued before it. `M-160` showed what removing it buys: CPU time flat at ~0.17 ms from 33³ to 129³.

**What `wgpu` 29 actually gives you, and it splits the two targets.** `PollType::Poll` is *"check the device for a single time without blocking"*; `PollType::Wait` is *"block until the given submission has completed execution"*. And, verbatim from the docs: **"On WebGPU, this has no effect. Callbacks are invoked from the window event loop."** So native Bevy has a real CPU stall to design away and the web build has no blocking primitive at all — meaning any code shaped around `Wait` is native-only scaffolding, and the restructuring must not become a `#[cfg]` fork, for the same one-path reason the `libm` justification already gives.

`TIMESTAMP_QUERY` is in `FeaturesWebGPU` and supported on Vulkan, DX12, **Metal**, OpenGL and WebGPU — the one feature here that behaves identically on both targets.

**C1 (instrument first).** Timestamp queries attribute `M-167`'s 83% into submit / execute / map / copy, and the largest single component is **map-wait**, not execute. *Falsified by:* execute being the largest, which would mean the arithmetic did move after all and `M-167`'s "the arithmetic never moved and was never the point" needs re-tiering.

**C2.** Feeding the vertex and index counts into `draw_indirect` from a GPU buffer — so the CPU never learns the count — removes at least **60%** of the measured synchronisation at 129³. *Falsified by:* under 60%.

**C3.** An N-frame-delayed double-buffered staging ring for the paths that genuinely need CPU-side data (collider generation) holds the amortised per-frame cost within one chunk of the budget across a 320× range, i.e. `M-124`'s property survives the added latency. *Falsified by:* the amortised cost drifting outside one chunk, which means the ring traded a stall for a queue.

**A design question that is the owner's, not the harness's.** C3 costs one to two frames of latency on collision. For a voxel game that is invisible; for a CAD tool it is a decision. The registration records the question rather than assuming an answer, per `CLAUDE.md`'s rule about design decisions.

#### P-72 — the granularity of the active-cell structure is a first-class parameter, and it has never been swept

**The measured prior is a 256× spread from one knob.** Hoetzlein, *GVDB: Raytracing Sparse Voxel Database Structures on the GPU*, HPG 2016, all timings on a Quadro M6000. Tree build at 2048³: ⟨3,3,3,3⟩ 616,444 bricks in **461 ms**; ⟨3,3,3,4⟩ 83,218 bricks in 69.8 ms; ⟨3,3,3,6⟩ 2,036 bricks in **1.8 ms**. Same resolution, same data, 256× apart, and the paper's own conclusion is *"larger brick sizes produce a fewer number of bricks resulting in faster tree changes."* Also measured there: octrees were 30–40% slower on node insertion than N³-trees.

**`P-40` chose 64 cells per word and never asked whether 64 was right.** `M-337` measured the stage at 5.5× and 12/12 bit-identical, which settles that the bitmap works and says nothing about its granularity. `P-39`'s Lipschitz brush pruning (`M-341`, 3.36× median) is the direct analogue of GVDB's topology cull, and GVDB's result says its yield should likewise be granularity-dependent.

**C1.** Sweeping the bitmap block granularity across 8³, 16³, 32³ and 64³ cells per pruning unit on a live edit trace produces a **pronounced optimum** — the best and worst edit-plus-remesh times differ by at least 2×. *Falsified by:* a flat curve, which is a genuine null and means `M-337`'s 64 was already at a plateau. That null is worth having and is the expected outcome on the smooth fields.

**C2.** The optimum is **field-dependent**, and specifically differs between `gyroid` (surface everywhere) and `fbm_terrain` (surface on a sheet). *Falsified by:* one granularity winning on both, which would make it a constant rather than a parameter.

**C3.** The spread is smaller than GVDB's 256×, and predicted below 4×. GVDB's figure is a tree rebuild on a 2048³ SPH volume, far larger than a chunk here, and its level-set numbers (5–6× over CPU) are consistently weaker than its volume numbers (60×). *Falsified by:* above 4×, which would be the more valuable outcome and would say the chunk-size regime is not the damping factor it looks like.

---

## Part 4 — Foreclosed, and why

The 2026-08-23 dossier's framing applies unchanged: *"each of these looks transferable under its lens and is not. Recording them stops the next sweep paying for them again."* Nine new rows, each chased far enough to be sure.

**1. Topological degree computation as a "this cell contains a surface" certificate.** It has a type error, not a cost problem. `deg(F, B, 0)` is defined for `F: ℝⁿ → ℝⁿ` and certifies *isolated points*; this crate has `f: ℝ³ → ℝ` whose zero set is a **surface**, so `deg(f, B, 0)` is not defined. Worse, every degree method requires `f ≠ 0` on `∂B`, and for a cell the isosurface crosses the boundary by construction — the hypothesis fails before the cost matters. Degree theory is usable here only against `∇f = 0` (three equations, three unknowns), i.e. to certify **critical points**, and that is a different experiment. Kearfott's method (`10.1007/BF01404868`) is finite only for a given refinement; Franek & Ratschan (`10.1090/S0025-5718-2014-02877-9`) prove termination *by induction on refinement* and report 5.6 s to 175 s on dimension-8 examples, which is the unbounded-iteration rule's exact target.

**2. Sperner's lemma / constructive fixed-point labelling.** The lemma is constructive but *locating* the panchromatic simplex is a path-following walk with worst-case exponential path length, and the black-box problem is PPAD-hard. Unbounded iteration with no useful bound. Out on structure.

**3. Cheng & Wen's uniqueness criterion.** It does give existence and uniqueness together in a finite per-box computation, and it is in the corpus. Their own stated limitation kills it here: *"our method is invalid for singular roots"* — and singular and degenerate configurations are exactly where isosurface case tables go wrong. It also remains an arXiv preprint with no journal publication found. Note the dossier already foreclosed the surrounding interval-Newton/Krawczyk family at `2026-08-23-discovery-dossier.md:294-298`; this is the one member of it worth having chased, and it still fails.

**4. Taylor models and higher-order inclusion functions.** Two reasons, either sufficient. First, for the **trilinear** they buy nothing: a multi-affine function is a convex combination of its eight corner values, so `[min corner, max corner]` is the *exact* range — measured here at 0 violations in 1,200,000 samples — and affine arithmetic, Taylor models and SOS all sit above an exact answer. (Worth knowing separately: naive interval arithmetic *on the trilinear expression* over-estimates that exact range by a median of **2.56×**, p95 3.34×, max 12.57× over 200,000 random cells. If any cell test in this crate evaluates the expression rather than comparing the corners, that is free tightness on the floor.) Second, for the **CSG tape**, Taylor models need dynamic-degree multivariate coefficient arrays — heap allocation proportional to tape complexity, operation counts that vary with the input, and the wrapping-effect failure mode Neumaier documents (`10.1023/A:1023061927787`). None of that survives golden hashes. **P-67's reduced affine form is the version of this idea that fits**, which is why it is registered and this is not.

**5. Sum-of-squares / Positivstellensatz certificates.** Finding the certificate needs an SDP solver: a large dependency, an iterative interior-point method with no termination bound, and floating-point output that is not bit-reproducible. Checking is cheap only with exact rational arithmetic, which `libm` does not provide. And it is the wrong hierarchy anyway — a box is a *polytope*, so the natural Positivstellensatz is **Handelman's**, which is an LP rather than an SDP, and the **Bernstein expansion over a box is already such a certificate in closed form with no solver**. For a multilinear `f` it collapses to the eight corner values, i.e. to row 4 above. If a stronger-than-interval emptiness certificate is ever wanted, it is Bernstein, not SOS — and `2026-08-23-unmined-mathematics-for-meshing.md:197-224` already foreclosed the Bernstein–Bézier form for the trilinear, with an explicit reopening condition (*"revisit the entire hypothesis if the crate ever adopts a tricubic reconstruction filter"*).

**6. Reproducible / order-independent floating-point summation.** Correct technique, and this crate does not have the bug. Ahrens, Demmel & Nguyen (`10.1145/3389360`) give bitwise-identical results regardless of summation order at 1.2–1.6× overhead — for the case where a *parallel reduction schedules nondeterministically*. A chunk-local extractor has a fixed operation order by construction, which is why `M-31`'s 216 golden hashes already match across macOS/arm64 and Linux/x86-64. Spending 1.2–1.6× on a hot-path sum to fix a problem that does not exist is a regression. **The reopening condition is specific:** if a rayon-style reduction with a nondeterministic combine order ever enters normal averaging or cross-chunk QEF accumulation, the cheap fix is a **fixed reduction tree** — deterministic by shape, zero arithmetic overhead — and only if the tree shape genuinely cannot be fixed does the accumulator earn its cost. Worth recording that it would fit if needed: six words is a plain `[R; 6]`, no `unsafe`, no allocation.

**7. Cache-oblivious and van Emde Boas layouts for the sample grid.** Three independent reasons. Bender, Kuszmaul, Teng & Wang (`10.1007/s00224-009-9242-2`, `arXiv:0705.1033`, in corpus and uncited) is **purely theoretical** — no experimental section, no runtime, no cache-miss count, no hardware, and no discussion of regular grids at all. Morin's careful array-layout study finds vEB only intermediate on a modern core and concludes *"the vEB layout doesn't become effective until the data gets really big — the index calculations are just too complicated"*; a 64³ `f32` chunk is 1 MB and largely resident. And Marching Cubes is a **linear sweep**, not a random-access traversal, so the hardware prefetcher is already winning — Nocentino & Rhodes measured Morton cutting memory transactions up to 10× at small blocks and **identical to linear at 16×16**, because linear already coalesces. This supersedes `B4` from the 2026-08-23 corpus sweep, which listed it as open on the strength of the optimality proof. The proof is real and is not evidence about this workload. *The cheap alternative that is indicated instead:* tile the sweep so a two-slice stencil window stays in L2, and check for a bend in `resolution_sweep`'s existing `t = a + b·n³` fit. A clean fit across 16³–256³ means there is no cache cliff to fix, and that null is worth more than the experiment.

**8. Temporal / incremental remeshing — reusing the previous frame's mesh across an edit.** There is **no published prior art with numbers**, which is the finding. The nearest hit, HVOFusion (`arXiv:2404.17974`), is RGB-D SLAM where the "edit" is a new sensor frame at the reconstruction frontier, and its own text never compares incremental extraction against full re-extraction — its baselines are neural SLAM systems. Building this would be original research with no baseline, in a repo whose rules forbid inventing algorithm details. And most of what it would buy is already banked by a cheaper route: `P-39`'s brush pruning skips chunks an edit cannot touch and `P-40`'s bitmap skips cells within a chunk, both *spatial* incrementality. True *temporal* incrementality additionally needs stable vertex identity across edits, which Marching Cubes does not give — a corner crossing zero re-cases the cell and changes its triangle count. *The better-evidenced target in the same direction* is Flying Edges' four-pass edge-centric structure, which computes exact output sizes before writing and eliminates the dedup hash; Kitware reports 10–100× over classic MC on a 16-thread laptop, a range too wide to plan against but structurally sound, and it needs no cross-frame state.

**9. The FlexiCubes / TetWeave line, and adopting a VDB-family structure.** TetWeave (`10.1145/3730851`) gets its watertight, 2-manifold, intersection-free guarantee entirely from Marching Tetrahedra on a **Delaunay complex** — a global triangulation, which is the one thing chunk-locality forbids. Strip the learning and the differentiability and what remains is "MT on a good tet grid", which is 1991. FlexiCubes' own contribution is per-cell weights learned by gradient descent; with no loss there is nothing to weight. Separately: **NanoVDB is disqualified by its own documentation** — *"values can be modified in a NanoVDB grid, its tree topology cannot"*, and a brush that carves new empty space outside the current narrow band **is** a topology change. GVDB is CUDA-only, dead against the Metal target; fVDB is a deep-learning framework; Aokana is a raymarched renderer whose DAG deduplication makes edits expensive by construction. **P-72 takes GVDB's brick-size finding without taking any of the libraries**, which is the only part that transfers.

**And one item foreclosed by timing rather than by structure:** mesh shaders on Metal, until Bevy moves to `wgpu` 30. See Part 2.

### Not foreclosed, not re-proposed — the prior sweeps' open rows

These stay live and this document defers to them rather than duplicating them. From `2026-08-23-unmined-mathematics-for-meshing.md`, ranked by that memo: **the margin certificate `f₍₂₎ − f₍₁₎`** (its own #1, one comparison in the existing fold, and still the best value-per-effort item anyone has proposed for this crate), the multi-affine `b₀ ≤ 4` enumeration, the pseudo-Boolean / Lovász restatement, the power-diagram pruner, PV certified subdivision for body saddles, the **box-constrained QEF** (the cell clamp is projection-after-solve, not the constrained minimiser, and nobody has measured the gap), and optimal rectangle partition as greedy meshing's baseline. From the corpus sweep: high-order seeding (Saye), the narrowed O(3)-equivariant vertex rule — which **P-61 partly subsumes**, since equivariance-by-construction in the *crossing* is a weaker version of the same idea and should be run first — and cellular sheaves, still flagged *"read before proposing"*. From the dossier's reserve: Nielson's `DisC[T]`/`G[T]` invariants, deterministic simplicial collapse, the A15 acute-tetrahedra stencil, and per-vertex edit provenance riding `P-39`.

---

## Part 5 — Acquisitions

Papers this document leans on that are **not** in the corpus, in priority order. Never guess these DOIs; each below was resolved through Crossref, OpenAlex or an authoritative publisher page, and the two unresolved ones say so.

1. **Stahl & Grosso, *MCPro*, GRAPP 2025** — `10.5220/0013309800003912`. Blocks P-65 entirely. SciTePress, open.
2. **Knoll, Hijazi, Kensler, Schott, Hansen & Hagen, CGF 28(1):26–40, 2009** — `10.1111/j.1467-8659.2008.01189.x`. The published IA-vs-RAA comparison P-67 is measured against.
3. **Finken, Li, Wang, Guo & Levine, `arXiv:2608.12142`** — free. P-66's edge test.
4. **Chernikov & Xu, IMR 2013** — `10.1007/978-3-319-02335-9_28`, and the 2015 *Computational Geometry* follow-up (**DOI unverified**). P-63's prior art and the reason its scope note is written the way it is.
5. **Ahrens, Demmel & Nguyen, ACM TOMS 46(3):22, 2020** — `10.1145/3389360`. Foreclosure 6's reopening condition; acquire so the condition is checkable rather than remembered.
6. **Smith, Levien & Owens, SPAA '25** — `10.1145/3694906.3743326`, open PDF at eScholarship. P-70's prior, and the Apple forward-progress finding.
7. **Baktash, Gillespie & Crane, TOG 2026** — `10.1145/3811358`, `arXiv:2606.00454`. **Check before acquiring**: the §3.1/§3.2.1/§3.2.2/§3.2.3 structure, Theorems B.4 and B.6, the Figure-13 subdivision stencil, Steiner fans, the contractible spanning disk and normal coordinates that `A-014`, `M-80`–`M-101` and `M-161`–`M-205` all reference make it very likely this **is** the paper the subgrid extractor already implements, in which case the acquisition is the published version rather than a new source.
8. **Etiene et al., *Topology Verification for Isosurface Extraction*, IEEE TVCG 2012** (**DOI unverified**; PubMed 21690649) — the 20,000-case suite MCPro validates against, and the closest published relative of this crate's own validity gate. Its companion `10.1109/tvcg.2009.194` **is** in the corpus and is cited nowhere.
9. **Nishidate & Fujishiro, PACMCGIT 7(1), I3D 2024** — `10.1145/3651285`. Mesh shaders for MC-style reconstruction; **its results table could not be retrieved from any source**, so treat every claim in it as unverified until read.

Two corpus documents worth pulling forward without acquiring anything: `10.1109/tvcg.2007.70429` (Entezari, Van De Ville & Möller, box splines on the BCC lattice — the claim that BCC reconstruction is *twice as efficient as tensor-product B-splines at the same sampling density* is an R-tier claim in this corpus until measured, and note the catalog misspells the first author as "Enterazi"), and `10.1007/s10851-017-0769-6` (Najman, Boutry & Géraud, well-composedness — the theory behind `P-41`'s census and `P-46`'s falsified repair, still uncited).

---

## Part 6 — Ordering

If the whole set is not run, this is the order. It is by evidence-over-effort, and it deliberately front-loads the two that are proofs rather than measurements.

| # | id | what | effort | why here |
|---|---|---|---|---|
| 1 | **P-61** | centred crossing offset | S | Pre-measured at 0 / 2,000,000 with a three-line proof, and it turns the octahedral group into a 48-relation test oracle |
| 2 | **P-63** | 2¹⁸ vertex-link sweep | S | Settles `O-12` for Marching Cubes by exhaustion, in seconds, and the fixture-can-fail control already exists as ✗43's pre-fix code |
| 3 | **P-69** | autovectorise field eval | S | 3.0× measured prior, no dependency, no `unsafe`, and C2 makes it self-rejecting if it is wrong |
| 4 | **P-71** | kill `poll(Wait)` | M | The largest number this project owns about itself (83%), and C1 costs one feature flag |
| 5 | **P-62** | PV `C1` certificate | M | The crate's first *positive* per-cell certificate, with a soundness oracle it already computes |
| 6 | **P-72** | granularity sweep | S | One bench, no dependency, and a null is a real answer |
| 7 | **P-64** | Kani + Flux | M | Dev-tools only; novel in graphics; the expected outcome is a null and it is still worth the proof |
| 8 | **P-66** | monotone-edge witness | M | Third attempt at a line that died twice, but a different mechanism with an oracle the first two lacked |
| 9 | **P-67** | reduced affine form | M | Turns `P-54`'s held result into something shippable in `no_std` |
| 10 | **P-70** | subgroup compaction | M | Real prior, real 1.5× ceiling, and a 1-D-workgroup blocker to clear first |
| 11 | **P-68** | per-vertex error bar | M | Depends on P-61 landing; a CAD feature with no competitor |
| 12 | **P-65** | MCPro | **L** | The largest item here and the one that would close `Error::UnresolvedSixSaddle`; do not start it before P-63 has said whether the manifoldness question is settled |

**Four of the twelve are expected to return nulls** — P-64's C2, P-72's C1 on smooth fields, P-61's C4, and P-67's C1 if the condensation is worse than Knoll's band suggests. That is registered on purpose. `M-249`'s directional-Lipschitz null and `M-329`'s modal kill-shot are two of the more useful entries in the ledger, and a phase where every clause holds is a phase whose clauses were too easy.

---

## Appendix — the P-61 pre-measurement

Run off-repo in CPython `f64`, to be re-run inside the harness in `Real` for both scalars. Reproduced here so the numbers in P-61 are checkable before any Rust exists.

```python
import random, struct
from fractions import Fraction as F

def ulps(x, y):
    ix = struct.unpack('<q', struct.pack('<d', x))[0]
    iy = struct.unpack('<q', struct.pack('<d', y))[0]
    if ix < 0: ix = -0x8000000000000000 - ix
    if iy < 0: iy = -0x8000000000000000 - iy
    return abs(ix - iy)

h, origin = 0.125, -2.0
random.seed(2026)
bad_cur = bad_mid = 0
for _ in range(300_000):
    i = random.randrange(0, 32)
    L, U = origin + h*i, origin + h*(i+1)
    a, b = random.uniform(0, 1), -random.uniform(0, 1)   # a > 0 > b: the edge is cut

    # current: parameter from the lower corner, then placed
    t  = a/(a-b);   x  = L + h*t
    t2 = b/(b-a);   x2 = (-U) + h*t2        # the mirrored cell: corner roles swap
    if x2 != -x: bad_cur += 1

    # proposed: signed offset from the edge midpoint, never via t
    d  = ((a+b)*0.5)/(a-b);   m  = (L+U)*0.5;        y  = m  + h*d
    d2 = ((b+a)*0.5)/(b-a);   m2 = ((-U)+(-L))*0.5;  y2 = m2 + h*d2
    if y2 != -y: bad_mid += 1
# h=0.125 -> current 25554/300000 (8.52%), proposed 0/300000
# h=0.1   -> current 74961/100000 (75.0%), proposed 0/100000
# h=3/32  -> current  9118/100000 (9.12%), proposed 0/100000
```

The exactness argument does not depend on the measurement, and is the reason to expect the zeros to survive translation into `Real`. Under the simultaneous swap `(p₀, a) ↔ (p₁, b)` composed with the coordinate sign flip: `fl(a+b) = fl(b+a)` because addition is commutative in IEEE 754; halving is exact because 2 is a power of two; `fl(b−a) = −fl(a−b)` because round-to-nearest is an odd function; and `fl(S / −D) = −fl(S / D)` for the same reason. Every step is a guarantee, not an observation. What the parameter form loses is not any of these — it is that reflection acts on `[0, 1]` as `0 ↔ 1`, an **affine** map, and floating point respects sign flips exactly and affine maps only approximately.

Verified independently against exact rational arithmetic, 300,000 world positions: current mean 0.086 ulp / worst 422 ulp, proposed mean 0.052 ulp / worst 757 ulp. Zero of 300,000 proposed offsets fall outside the cell.
