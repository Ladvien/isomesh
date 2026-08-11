# isomesh — findings ledger

**Started:** 2026-08-11 · **Append-only.** Entries are never deleted, only re-tiered when new evidence
arrives — with the old verdict left visible.

This is the project's epistemic state: what we believe, how strongly, and **on what evidence**. It
exists because this project has already been wrong six times in ways that would have silently
propagated into code, and because the research corpus contains several published figures that failed
verification. A belief with no recorded falsification method is not a finding, it's a preference.

## How to use this file

- **Before acting on a "known" fact, look for it here.** If it's not here, it hasn't been checked.
- **When a measurement contradicts something written down, that's an entry** — the contradiction is
  the finding, not an inconvenience.
- **Falsified entries stay.** They're the most valuable rows, because they tell you which *sources* to
  distrust, not just which facts.
- Every entry names how it could be shown wrong. If you can't write that line, you have an opinion.

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

**Consequence:** SN's case must rest on quad connectivity and inner-loop cost, not output size. In
M-001 the count columns are a **checksum with a predicted value**, not a result.
**Would be shown wrong by:** any closed-manifold pair on the same grid where the difference ≠ 2χ.
**Legitimately breaks at:** boundary-clipped meshes (**incoming at G-001 chunking** — expect the
assertion to fail there and do not "fix" it), A-013 welding, MC vs MC33 differing in χ.

### ✗2 — "You can have a manifold mesh or an intersection-free one, not both"

**Believed because:** folklore, repeated in several secondary sources.
**Falsified by:** literature review round 1. Manson & Schaefer 2010 achieved both. ODC (2024) measured
Manifold DC at **100% of models self-intersecting** against ODC at **0 of 1500**.
**Consequence:** guaranteed intersection-free extraction is on the table, which is the premise under
A-009 and the runtime-convex-decomposition opportunity.

### ✗3 — "Every interior Surface Nets vertex has four neighbours"

**Believed because:** written into isomesh's own module docs before measuring.
**Falsified by:** A-004. Measured max degree **10** — higher than MC's **9**.
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

### ✗9 — "MC's cost inside the volumetric loop was never measured" / "navmesh rebuild cost was never measured"

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

### ✗11 — "Plain MC has ambiguous faces and produces holes"

**Believed because:** stated in this repo's own implementation brief (Stage 2, "Plain MC has ambiguous
faces and produces holes"), carried into `BACKLOG.md`'s A-002 acceptance criterion, and near-universal
folklore about Marching Cubes.
**Tested by:** `validate_table()` (`mc/mod.rs:319`), which checks all 256 cases structurally, and the
assertion `assert_eq!(report.face_disagreements, 0)` at `mc/tests.rs:30`.
**Result:** zero face disagreements, across all 256 cases.
**What's true instead:** holes require two cells sharing a face to *disagree* about how the surface
crosses it. In this implementation a face's segments are a function of that face's own four corner
signs and nothing else — the two cells meeting on a face read the same four corners, so they cannot
disagree. The property is structural, not empirical, and it falls out of the table being **derived at
compile time by walking each face counter-clockwise** rather than transcribed from a diagram.

The folklore is not wrong about Marching Cubes in general; it is wrong about *this* Marching Cubes.
Lorensen & Cline's original table was transcribed per-case and its ambiguous cases were resolved
inconsistently between complementary configurations, which is where the holes came from.

**Consequence:** A-002's acceptance criterion was unsatisfiable and has been re-scoped. MC33's
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

Worse, the audit's diagnosis of *why Dual Contouring pops* is the branch itself: DC's hard SVD
truncation at σ < 0.1 is a discontinuous branch, and over 20,000 trials seeded at the threshold in f32
the rank branch disagreed after a rotation in **454 cases**, with `‖f(Rx) − Rf(x)‖` median **2.13** and
max **9.10** — a several-cell vertex pop from an infinitesimal rotation. A triple-product threshold is
the same construction with a different discriminant, so the split would have reintroduced the exact
failure the rule exists to remove.

The measured equivariance residual (f32, coordinates in [0,256], 4000 random cells) also shows the
"fast path" is not the accurate one:

| rule | median | p99 | max |
|---|---:|---:|---:|
| DC normal equations | 6.80e−05 | 2.48e−01 | 5.6e+02 |
| dual basis (Cramer) | 1.61e−05 | 7.23e−04 | 3.6e−01 |
| **Tikhonov adjugate** | **1.59e−05** | **1.81e−04** | **6.4e−04** |

Tikhonov dominates Cramer on both tail columns, so nothing is traded away by dropping the three-plane
form. The two paths also do not agree to within noise, which means the branch would have been
*observable* in the output.

**Consequence:** A-007 and A-008 merged into one ticket with one unconditional path. Two requirements
the audit states and no ticket had recorded are now in it: **magnitude-sorted 3-term dot products**
(4328/9600 equivariance failures unsorted, **0/9600** sorted — the guarantee does not hold in f32
without this), and the derivation of **λ = 0.01** as the value that reproduces DC's σ = 0.1 truncation
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

So the 3×3 lives in `dc/solve.rs` as a six-entry symmetric matrix over `[R; 3]`, about 40 lines, and
**the crate stays at one dependency.** The "as light as possible" pitch survives A-007 intact.
**Would be shown wrong by:** glam gaining a generic scalar parameter, or this crate dropping `f64`.
**Note this is ✗10's second correction.** ✗10 moved glam from "day one" to "A-007"; the deferral target
was wrong too. The recurring error is reasoning about glam from its reputation rather than its API.

### ✗15 — "Marching Cubes is unconditionally manifold"

**Believed because:** every measurement in this repo said so. M-4 contrasts Surface Nets'
non-manifoldness against *"Marching Cubes' zero at every resolution"*, the README says *"Marching Cubes
stays manifold"*, and the mechanism looked airtight — MC places vertices on grid **edges** rather than
one per cell, so the multi-sheet argument that sinks SN does not apply. `SurfaceGate::Closed`'s own doc
comment asserted it.
**Falsified by:** T-005b's `marching_cubes_meshes_sphere_unions`, on its first run against a fresh
proptest seed — during T-006, which is a nice demonstration that the property tests keep working after
the ticket that wrote them.
**Result:** a union of three spheres at `h = 2/3` gives **2 non-manifold edges and 3 non-manifold
vertices** on a mesh that is otherwise perfect: closed, `χ = 2`, one component, consistently oriented,
zero boundary edges.

**What's true instead:** MC is manifold when the grid **resolves** the surface. Where the surface
*pinches* inside a single cell — two lobes of a union meeting at sub-cell scale — the shared grid edge
ends up carrying four faces. Refinement fixes it, and sharply:

| n | 7 | 9 | 13 | 17 | 25 | 33 | 49 | 65 |
|---|---|---|---|---|---|---|---|---|
| non-manifold edges | **2** | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

**Consequence:** the property suite's gate is renamed from `ClosedAllowingMultiSheet` to
**`ClosedAllowingUnresolvedTopology`** and now covers MC on generated fields too, because the condition
was never about the algorithm — it is about whether the grid resolves the field. The strict gate is
still asserted where it is actually true, on the seven reference fields in `mc/tests.rs`. The exact
counts are pinned in both directions by
`an_under_resolved_pinch_makes_marching_cubes_non_manifold`, following M-4's precedent: the defect is
an assertion, not an exclusion, so it fails if it spreads *and* if it silently disappears.
**Would be shown wrong by:** the same field at `h = 2/3` coming back manifold, or any *resolved* field
coming back non-manifold under MC.
**Worth noting for G-001:** a chunk boundary is a place where a surface can be under-resolved relative
to the chunk it lands in. This is a plausible source of seam defects later.

### ✗14 — "Surface Nets is the cheapest thing in the family and the natural default"

**Believed because:** this repo's own algorithm catalog states it as the engine verdict —
`docs/research/2026-08-10-meshing-algorithm-catalog-v2.md:163`, *"cheapest thing in the family and the
natural default"* — reinforced by the same folklore ✗1 already corrected once.
**Tested by:** T-006's resolution sweep, `cargo bench --bench resolution_sweep`. Sphere, `f32`, single
thread, median of 5 timed runs after 2 warm-ups, identical grid and reused output buffers for both
algorithms. Raw data committed at `docs/measurements/resolution_sweep.csv`.
**Result:** Surface Nets is cheaper only below about 48³. The crossover sits between 48³ and 64³, and
past it Surface Nets loses steadily and then sharply:

| n | MC ms | SN ms | SN/MC | MC ns/sample | SN ns/sample |
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
would localise this to one cache hierarchy rather than to the algorithm. **This is a one-machine
result** — Apple M5, single thread — and nothing here has been run elsewhere. See O-11.

**Narrowed by A-007.** Dual Contouring, which shares the dual topology and differs only in vertex
placement, shows the *same* curve: `218.9 ms` at 256³ against Surface Nets' `212.5`, and the same
negative fitted intercept. So the superlinearity is a property of the **shared dual engine** — its
sampling and its strided quad walk — and not of either vertex rule. That is a real narrowing of O-11's
search: the suspect is the gather in `emit_quads`, whose `z`-stride is `n²` cells apart.

---

## Part 2 — Measured here (tier M)

| # | Finding | Evidence |
|---|---|---|
| M-1 | **surface cells = crossed edges + χ** | 4 topologies × 3 resolutions, table in ✗1 |
| M-2 | `V_sn = V_mc + χ`, `F_sn = F_mc + 2χ` | A-004 tests, all four clean fields |
| M-3 | Surface Nets max vertex degree **10**; MC **9** | A-004 |
| M-4 | SN is non-manifold where one cell carries two sheets: **48** non-manifold edges on capped gyroid, **15** on fbm_terrain at 33³ | A-004; pinned as non-zero assertions, not excluded silently |
| M-5 | On `box_exact`, SN's nearest vertex to the corner (1,1,1) is **1.15 cells** away | A-004 — this gap is what E-104 exists to show |
| M-6 | `libm::sqrtf` lowers to hardware `fsqrt` (aarch64+neon) / `sqrtss` (x86-64+sse2) | libm 0.2.16 source: `src/math/arch/aarch64.rs` raw asm, dispatched by `select_implementation!` on `target_feature` |
| M-7 | dev-dependencies do not propagate: consumer resolves **3 packages**, the crate's own lockfile has **137** | Experiment, cloud container |
| M-8 | Cargo silently co-resolves two wgpu majors — **317 packages, both 29.0.4 and 30.0.0**, no resolution error; fails later as `expected TextureFormat, found a different TextureFormat` | Experiment |
| M-9 | Workspace feature unification leaks: `-p isomesh` alone gives glam `libm`; whole-workspace gives it `std`, `serde`, `bytemuck`, `encase`, `rand` | Experiment — the reason `bevy_isomesh` is excluded |
| M-10 | **Unit sphere at 64³ (`h = 4/63`), symmetric Hausdorff: MC `1.380e-3`, SN `2.288e-3`.** Mean absolute error MC `6.50e-4`, SN `1.367e-3`. SN is **1.66×** worse than MC on both | T-003, `a_unit_sphere_at_64_cubed_is_within_one_cell_diagonal` |
| M-11 | **T-003's own acceptance criterion is loose by ~80×.** One cell diagonal is `0.10997`; MC measures `0.00138`. A harness returning a constant `0.01` would pass it | T-003 — which is why the ticket also ships a convergence-order test and closed-form fixtures |
| M-12 | **MC's error falls like `h²`, measured.** Mean error `2.7168e-3` at 32³ against `6.5015e-4` at 64³ — a ratio of **4.179**, against the ideal `((4/31)/(4/63))² = 4.13` | T-003, `the_error_falls_like_h_squared` |
| M-13 | **Surface cells ≈ `1.5·A/h²`, not `A/h²`.** Measured `1.450` (25³), `1.442` (33³), `1.517` (64³) on the unit sphere. The constant is derivable: a plane of unit normal `n` crosses `(\|nₓ\|+\|n_y\|+\|n_z\|)/h²` cells per unit area, and `E[\|nₓ\|] = ½` over the sphere, so an isotropic surface gives `E[Σ\|nᵢ\|] = 3/2` | T-003. Predicted 6,430 triangles at 64³ from `A/h²` and measured **9,452** — a 1.47× miss, which is this factor |
| M-14 | **The reverse direction finds defects the forward direction structurally cannot.** `box_exact` at 33³: forward `0.0833`, reverse `0.1443` — the reverse number is MC's rounding of the sharp corner. `thin_plate` at 33³: forward `0.0083`, reverse `0.0893` — an under-resolved plate | T-003. Deleting one face of an octahedron leaves `mesh_to_field` bit-identical and moves `field_to_mesh` to `√(3/2 − 2/√3)` |
| M-15 | **Surface Nets' non-manifoldness is a resolution effect, not a topology effect.** M-4 measured it on `gyroid` (48 edges) and `fbm_terrain` (15) and read it as a high-genus / open-field property. T-005b finds it on a randomly generated **convex body** — 1–2 non-manifold edges, 3–4 non-manifold vertices, zero boundary edges. Any feature thinner than one cell forces two sheets through it | T-005b, `surface_nets_meshes_convex_bodies`. This is why the sweep has a named `SurfaceGate` rather than a per-field exception |
| M-16 | **The even-`χ` parity check is not independent of manifoldness — it is a corollary of it.** `χ = 2 − 2g` holds for a closed *orientable manifold*, so a gate that waives manifoldness and keeps parity is incoherent. Measured: SN on a generated convex body gives **`χ = 1`** with one non-manifold edge and zero boundary edges | T-005b. Cost one wrong gate before it was noticed; `SurfaceGate::ClosedAllowingMultiSheet` now documents the omission rather than leaving it to be rediscovered |
| M-17 | **A case-table entry naming an *uncut* edge is caught inside the crate, before any mesh exists** — `edge_crossing`'s `is_inside(a) != is_inside(b)` precondition fires. Real defence, but it is a `debug_assert`, so it is absent from a release build | T-005b. The mutation check therefore confines its wrong-edge corruption to *cut* edges, which is both the plausible transcription error and the one that actually reaches the validity gate |
| M-18 | *(refined by T-008 — the arithmetic below is about **adjacent cells**; the effect on a real mesh is **gradual**, because a mesh whose vertices span many lattice cells keeps resolving them for a while after neighbours have merged. Measured on a 1158-vertex sphere at `h = 0.125`: no collapse at all at `1e4`, and **1158 → 918 buckets**, a 21% loss, at `1e6`. Fixed by anchoring the lattice to the mesh's own minimum corner, so the scale depends on the mesh's extent rather than its position)* **`quantise`'s weld lattice collapses beyond ~105 world units, and it is a performance cliff rather than a correctness one.** It scales absolute coordinates by `1/(h·1e-4)` — 160,000 at `h = 0.0625` — so it passes `f32`'s exact-integer range at `2²⁴·weld_epsilon ≈ 104.86`. Measured: at `p = 104` two cells one apart stay distinct; at `p = 105` they collapse; by `p = 1000` a whole region is one bucket. **Correctness survives** — coarsening only *merges* buckets, and the 27-neighbour probe plus exact distance test still finds every duplicate — but the scan degrades toward quadratic, silently, at exactly the coordinates G-001 chunking and G-007 streaming produce. `TriangleGrid` is immune because it quantises *relative* to its own AABB origin, which is the fix pattern | T-005b follow-up, ✗13. Ticketed as T-008 |
| M-19 | **There is no meaningful fixed cost on the CPU extraction path, and the prediction saying so was written down before the run.** MC's fitted `a` is `0.49 ms` against a largest measured run of `80.3 ms` — **0.61%**. V-6's "73% of a published 64³ figure was fixed launch overhead" is a *GPU dispatch* property and does not transfer; the "stop trusting single-grid numbers" rule belongs to Phase 6, not here. **Caveat that matters:** `a` is 543% of the *smallest* measured run, so the fit must not be extrapolated below 16³ — down there the `O(n²)` surface term dominates | T-006, `benches/resolution_sweep.rs`. The prediction is in that file's module docs, committed before the first measurement |
| M-20 | **Marching Cubes' marginal cost is `4.75 ns/sample` — `211 M samples/s`, single-threaded, `f32`, Apple M5.** Per-sample cost is flat within 2% from 128³ upward, `r² = 0.99986` | T-006. Against V-3's `5.42 G voxel/s` on an RTX 2080 Ti that is a **~26× gap**, which is the number the Phase 6 GPU decision should be argued from rather than from folklore |
| M-21 | **Surface Nets is not `O(n³)` over this range; Marching Cubes is.** SN's fitted intercept is **negative** — `−3.13 ms` full sweep, `−7.32 ms` on the tail — which is physically impossible and is the signature of a curve convex in `n³`. `r² = 0.9899` against MC's `0.99986`. Per-sample cost rises `9.0 → 13.19 ns` while MC's falls and flattens | T-006. Cause unmeasured — see O-11. This is why ✗14's gap widens rather than staying constant |
| M-22 | **✗1's identity holds at every resolution to 256³**: `V_sn − V_mc = 2` and `F_sn − F_mc = 4` exactly, nine resolutions, `χ = 2`. The original table topped out at 49³, so this is corroboration at **5× the resolution** and 16.8 M samples | T-006's sweep, which records vertex and triangle counts alongside the timings |
| M-23 | **`f64` costs 8–10% on extraction paths with no matrix solve in them.** At 65³ on a sphere: MC `1.3928 ms` (f32) against `1.5083 ms` (f64), **+8.3%**; SN `2.3625` against `2.6036`, **+10.2%**. Not the 2× a naive "twice the bytes" guess suggests, because the work is dominated by field evaluation and branchy table lookup rather than by memory bandwidth | T-006, `benches/extract.rs`, the `precision` group. **Partially answers O-8** for the non-QEF paths; A-007's solve is where `AᵀA` squares the condition number and the answer may differ |
| M-24 | **Bit-exact lattice equivariance needs magnitude-ordered *products*, not just sums.** The audit prescribes "magnitude-sorted 3-term dot products", which is necessary and **not sufficient**: a cofactor expansion of `det(M+λI)` along a fixed row selects three of the six entries *by position*, so relabelling the axes evaluates a different expression. Measured **19 ULP** disagreement under a cyclic permutation, on all three fixtures, with the dots already sorted. Fixed by the symmetric determinant form with magnitude-ordered 3-factor products — FP multiplication is commutative but not associative, so `(a·b)·c ≠ (b·c)·a`. Now **72/72** rotation×fixture cases are bit-identical | A-007, `the_vertex_is_bit_exactly_equivariant_under_lattice_rotations`, which failed before the fix |
| M-25 | **The sharp-feature solve is nearly free: Dual Contouring costs 3% more than Surface Nets.** At 256³ on a sphere, `218.9 ms` against `212.5 ms`; marginal `78.1` against `80.4 M samples/s`. A full 3×3 regularized solve per surface cell, and it barely registers — because both methods are dominated by the *shared* dual topology (sampling and the quad walk), not by vertex placement | A-007, `benches/resolution_sweep.rs`, now sweeping three algorithms |
| M-26 | **Dual Contouring reaches a box corner to `0.01` cells where Surface Nets stops at `0.58`** — measured at 27³ on `box_exact`, `0.0009` against `0.0888` in world units. The resolution is deliberately **not** grid-aligned; on an aligned grid this measures the zero-classification rule instead (E-103's trap) | A-007, `the_corner_is_sharper_than_surface_nets`. This is E-104's money shot, measured before the example exists |
| M-27 | **The two dual methods differ *only* at features, with a 14-order-of-magnitude gap and nothing in between.** On `box_exact` at 27³: **864** of 1016 vertices agree with Surface Nets to within `2e-15` cells, **152** move by `0.35`–`0.57` cells, **0** land between. Exact reason: on a planar patch every crossing lies in the plane, so `pᵢ − c ⊥ n`, every `dᵢ` is exactly zero, `g` is exactly zero and the solve returns the centroid | A-007, `dual_contouring_moves_only_the_feature_vertices`. Consequence: E-104's side-by-side measures the feature and nothing else. Note the agreement is to *rounding*, not to the bit — the two centroids are computed by different expressions |
| M-28 | **The cell clamp eliminates placement-caused self-intersections entirely, and costs nothing in sharpness.** λ (pairs per 1,000 triangles) at 33³, clamp off → on: `torus` **2.66 → 0**, `gyroid` **71.43 → 3.12**, `fbm_terrain` **189.46 → 13.84**; `sphere`, `box_exact`, `csg_difference` and `thin_plate` were already 0. Corner gap on `box_exact` at 27³: **0.0057 cells either way** — identical, because a convex corner's solution is interior to its own cell, so the constraint never binds where the feature is | A-009, `the_clamp_measured_on_every_reference_field`. Default is now `Clamp::ToCell`, chosen by this measurement rather than by preference |
| M-29 | **The literature's two branches both fire, on disjoint fields — which is a sharper answer than either alone.** The review states the rule in advance: λ→0 means placement was the cause, λ unchanged means the defect is *connectivity* and needs A-010. Measured: λ → **exactly 0** on five of seven fields, and drops **23×** and **13.7×** on `gyroid` and `fbm_terrain` without reaching it. Those two are precisely the fields with multi-sheet cells (M-4, M-15). So the clamp removes the placement failure completely and the residue is exactly A-010's problem, with nothing left unaccounted for | A-009 |
| M-30 | **An unclamped solve can fling a vertex 3.18 cells out of its own cell** — measured max displacement on `gyroid` at 33³, with 618 of 5240 vertices outside; `fbm_terrain` 2.17 cells and 1097 of 1958. On the smooth closed fields it never leaves at all: `sphere`, `box_exact` and `thin_plate` have **zero** vertices outside | A-009. This is the failure mode the clamp exists for, quantified rather than asserted |

---

## Part 3 — Verified from primary sources (tier V)

| # | Finding | Source |
|---|---|---|
| V-1 | Bevy 0.19 pins **wgpu / wgpu-types / naga 29.0.3, glam 0.32.0, encase 0.12** | `bevy_render/Cargo.toml` @ v0.19.0 |
| V-2 | Bevy 0.19 **removed `RenderGraph`**; passes are systems in ECS schedules; non-camera work targets the `RenderGraph` schedule | 0.18→0.19 migration guide |
| V-3 | MC peak: **5.42 G voxel/s, 330 M tri/s** (RTX 2080 Ti). DMC costs 1.52–3.50×; FlexiCubes 2.77–3.92× | Grosso & Zint, `10.1007/s00371-021-02139-w` |
| V-4 | **Contouring 68 ms vs halfedge construction 58 ms** — extraction is 54% of a usable mesh | same, Table 5 |
| V-5 | On unstructured grids, Delaunay/MT ratio **15.3×–81.5×** — contouring is 1–2% of the pipeline | TetWeave Table 3 |
| V-6 | **73% of FlexiCubes' 64³ MC timing is fixed launch overhead** (fitted a ≈ 1.88 ms) | fit over FlexiCubes' own resolution series |
| V-7 | Cross-paper reproducibility floor is ~1.5× **in opposite directions**: TetWeave re-measured FlexiCubes at 128³ as 9.63/15.25 ms vs FlexiCubes' own 14.06/9.53 | both papers |
| V-8 | GPU MC throughput has not tracked hardware: **10.7× more bandwidth bought ~1.7× more throughput** (GTS 450 → 2080 Ti) | speed analysis |
| V-9 | Same MC, compute shader → mesh shader: **114.2 → 2679.4 fps (23.4×)** | Elliott MSc, Waikato 2022 |
| V-10 | CBT sum-reduction, atomics → LDS staging: **5.78 → 0.40 ms** | Unity SIGGRAPH 2021 (see ✗7) |
| V-11 | Meshlet compression: **15.5 M tri in 0.59 ms** (RX 7900 XTX) | `10.2312/vmv20241204` |
| V-12 | Work graphs: 79,710 instances in 3.74 ms — **but 2.8–3.4× slower** on classification workloads | `10.1145/3675376` + independent profile |
| V-13 | nvblox: meshing is the least GPU-accelerable stage, **×3–13 vs fusion's ×174–177** | nvblox |
| V-14 | Aokana renders **10¹⁰ voxels at 6 ms**, 5% resident, RTX 3060 Ti — **explicitly not editable** | Aokana |
| V-15 | CoACD vs V-HACD: **49% → 80%** downstream manipulation success | CoACD |
| V-16 | Dimforge migrated parry (0.26.0) and rapier (0.32) **off nalgebra onto glam**, citing rust-gpu support; performance *"nothing changed, at all"* | dimforge.com, 2026-01-09 |
| V-17 | **No paper since 2020 benchmarks MC vs Surface Nets vs DC against each other.** Surface Nets has no credible published timings at all | literature review round 1 |
| V-18 | **DC's own paper quantifies the f32 QEF failure.** At 256³, `bᵀb` reaches ~10⁶; f32 carries six decimal digits, so `E[x]` evaluated on a flat region — where it should be zero — has error **on the order of 1**. The paper's own remedy is double precision | Ju, Losasso, Schaefer & Warren 2002, `10.1145/566570.566586`, §2.3, read this session |
| V-19 | **DC's topology is Surface Nets' topology.** The paper's algorithm is literally: vertex at the QEF minimizer for each sign-changing cube, quad joining the four cubes of each sign-changing edge. Only vertex *placement* differs | same, §2.2 |
| V-20 | A QEF is stored as `AᵀA` (symmetric 3×3), `Aᵀb` (3-vector) and `bᵀb` (scalar) — 10 floats — rather than as `A` and `b` | same, §2.3 |

---

## Part 4 — Open questions

Each has the test that would settle it. **An open question with no proposed test is a wish.**

| # | Question | Settled by | Why it matters |
|---|---|---|---|
| O-1 | What fraction of cells actually change per brush stroke? | G-002 instrumentation; hash cell slabs, log per stroke | **Unpublished.** Ceiling on every incremental-repair idea in the opportunities doc |
| O-2 | Does clamping the QEF vertex to (1−ε) inside its cell eliminate self-intersections? | A-009: measure per 1,000 triangles, clamp on vs off, all seven fields | Decides whether guaranteed intersection-free extraction is free → whether runtime convex decomposition can stop failing |
| O-3 | MC vs SN vs DC vs MT — actual relative speed on one machine? | M-001 | The comparison does not exist (V-17). We'd have the only apples-to-apples measurement |
| O-4 | Do brush operations commute? | G-003: 8 ops × 40,320 orderings, count distinct results. Expect 1 | If not 1, the coordination-free multiplayer story dies — cheaply, before anything is built on it |
| O-5 | Do mesh shaders work on macOS/Metal? | GPU-007 capability probe | **Sources contradict:** wgpu's spec table lists MSL as *planned*; the tracking issue says the Metal HAL backend merged. Neither is trustworthy until probed |
| O-6 | What is amortized meshing cost per frame under continuous editing? | E-206 under a deliberately overloaded queue | The only number a game cares about, and no paper reports it |
| O-7 | What fraction of *our* pipeline is contouring vs everything else? | M-003 | V-4 says 54% for someone else's code with no physics. Ours is probably worse |
| O-8 | Does DC's vertex placement need f64 in practice, or is f32 enough? | E-112, with the QEF condition number in the HUD | `M = AᵀA` squares the condition number. **Half answered by V-18**: the original paper measures f32 error ~1 on flat regions at 256³, and recommends f64. **Partially answered by M-23**: on extraction paths with no solve, `f64` costs only 8–10% of wall time, so precision is cheap where there is no QEF. Still open for the vertex solve itself, and for *our* fields at *our* resolutions — which sidesteps `AᵀA` entirely and may not degrade the same way |
| O-9 | How much does T-003's gradient-flow chord **over**-estimate distance at a concave seam? | A comparison against nearest-point search over a dense surface point cloud, or E-104 once DC lands | The chord follows `∇f` to the zero set, which near `csg_difference`'s seam can land further away than the true nearest point. The bias direction is known and safe for a "below X" gate; the *magnitude* is not measured, and M-001's shootout column would inherit it. `csg_difference` measured forward `0.0833` at 33³ — how much of that is seam bias is unknown |
| O-10 | What is Surface Nets' non-manifold **rate** as a function of feature thickness over `h`? | A-010, which must drive it to zero; a sweep over a slab of shrinking thickness would answer it sooner | M-15 established it is a resolution effect rather than a topology one, and M-4 has counts at two resolutions on two fields. Nobody has the curve. It decides whether SN is usable at game resolutions or needs A-010 first |
| O-11 | **Why does the dual topology go superlinear in `n³` while Marching Cubes does not?** *(narrowed by A-007: Dual Contouring shows the same curve, so it is the shared engine rather than either vertex rule — see ✗14)* | A profile or cache-miss counter at 192³ vs 256³, or simply re-running T-006's sweep on another machine — `big` (Ryzen 9 5900X) would separate algorithm from cache hierarchy for the cost of one `cargo bench` | M-21 measures the effect (`9.0 → 13.19 ns/sample`, with the sharp step between 192³ and 256³) without explaining it. The hypothesis is working-set: SN gathers the four cells around each crossed edge, one of those strides being `n²` apart, which at 256³ is 65,536 cells. **Untested.** ✗14's verdict rests on this being a property of the algorithm rather than of one machine, so it is the highest-value cheap experiment currently open |

---

## Part 5 — Method rules, and the failure that earned each

Rules with no incident behind them get ignored. These all have one.

| Rule | Earned from |
|---|---|
| A typed error at the call site is louder than an abort — make the invalid state unrepresentable where you can, report it where you can't, and never substitute a default | The no-panic rule, reconciled with "fail loudly": `ValidateConfig` has private fields and one checked constructor, so the validator needs no runtime guard at all |
| Corpus presence is decided by `catalog_read`, never by `distill_search` | ✗4 — 342 documents readable but unsearchable |
| **Never guess a DOI or arXiv ID.** Look it up or stop | A subagent guessed an ID from memory and downloaded an unrelated condensed-matter physics paper under a meshing DOI |
| Verify the *source* separately from the *number* | ✗7 — right figure, wrong attribution |
| "Nobody measured X" needs the same evidence as any other claim | ✗9 — asserted twice, false twice |
| No performance number without the benchmark that produced it, in the repo | The corpus contains several published figures that failed verification |
| A doc comment the test suite disproves is worse than no doc comment | ✗3 |
| Assert the identity, not the inequality — a weak assertion hides a strong fact | ✗1 |
| **Verify that a property test can actually fail.** Corrupt an input and confirm red | A test that cannot fail is decoration, not evidence |
| Record an assertion's break conditions *next to it* | ✗1 — G-001 chunking will break it correctly, and it will look like a regression |
| Pin known defects as non-zero assertions rather than excluding them | M-4 — the numbers only move when someone means them to |
| Single-grid timings measure dispatch latency; sweep resolution and report the fixed cost | V-6 |
| Treat any published cross-paper ratio below ~2× as noise | V-7 |
| A green local run on one platform is not a green build. CI is the first real test of anything platform-shaped, and it will find things a local pass structurally cannot | First push: every job passed except `bevy_isomesh` on Linux, where Bevy 0.19's default Wayland backend needs `libwayland-dev` / `libxkbcommon-dev`. No such package exists on macOS, so no amount of local verification could have caught it |
| **A ticket's acceptance criterion is itself a claim about the code. Check it against the code before starting the ticket, not after.** | ✗11 — A-002 carried an `L`-sized acceptance criterion that the existing test suite had already made unsatisfiable. Nothing flagged it, because acceptance criteria are read as instructions rather than as assertions to verify |
| A property that falls out of *how a table was constructed* outranks folklore about the algorithm the table implements | ✗11 — "MC produces holes" is true of a transcribed table and false of a derived one; the distinction is invisible if you reason about "Marching Cubes" rather than about this code |
| **When a ticket paraphrases a research doc, re-read the doc.** A paraphrase can invert the property that made the technique worth adopting | ✗12 — "branch-free, handles all degeneracies" became "falls through when the triple product is near zero" across three documents, turning the rule's central guarantee into its opposite |
| **When an acceptance criterion passes by two orders of magnitude, it is not the test — find the one that fails.** Ship it anyway, and ship the real one beside it | M-11 — T-003's stated criterion passes with 80× margin, so a constant-returning harness satisfies it. The convergence-order test and the closed-form fixtures are what actually constrain the code |
| Estimate a count from the geometry, then **measure it before writing it down** — the tidy formula is usually missing a constant | M-13 — `A/h²` under-predicted the triangle count by 1.47×, because a surface crosses `3/2` cells per unit area, not one |
| **A test double is evidence only if it is pinned to the thing it stands in for.** Write the equivalence test first, and make it bit-exact | T-005b — the case-table mutation check runs a corrupted table through a local marcher. Without `the_double_reproduces_marching_cubes` comparing the two bit-for-bit on the *uncorrupted* table, a corrupted-table failure would be indistinguishable from the double having drifted |
| **Before waiving a property in a gate, check what else was resting on it.** Derived checks fail with their premise | M-16 — a gate that waived manifoldness kept the even-`χ` assertion, which is a *corollary* of manifoldness. It failed on the first non-manifold mesh it saw, and the assertion was the bug, not the mesh |
| When a mutation test passes for the wrong reason, the message tells you: **check where the panic came from, not just that there was one** | M-17 — the wrong-edge corruption first tripped `edge_crossing`'s precondition deep inside the crate, so the validity gate under test was never reached. `should_panic` alone would have called that a pass |
| **Write the prediction into the benchmark before the first run**, in the file, committed. Then the result cannot be rationalised afterwards | M-19 — T-006 predicted `a ≈ 0` on the CPU path from the fact that V-6's 73% figure is GPU dispatch overhead. It came out at 0.61%. Had it come out large, the prediction being on record is what would have forced the awkward question instead of a tidy story |
| **A fitted coefficient means nothing until it is compared to the data's own range.** Report it against both ends | M-19 — MC's `a` is 0.61% of the largest run and 543% of the smallest. Either number alone tells a different and misleading story; the pair says "negligible at scale, do not extrapolate below the range" |
| **A physically impossible fitted parameter is the model telling you it is wrong.** Do not report it as a value | M-21 — Surface Nets' fitted fixed cost is *negative*. Reported as "there is no fixed cost" it would be nonsense; read correctly it says the cost grows faster than `n³` and the whole two-term model does not apply |
| **A property that has held in every measurement so far is still a hypothesis, not a mechanism.** Say which condition it depends on | ✗15 — "MC is manifold" held on seven reference fields at every resolution ever tried, and the mechanism offered for it (vertices on edges, not one per cell) was real but insufficient. The true condition is "the grid resolves the surface", which nothing had stated |
| **A workspace that is excluded from the root is excluded from the root's CI commands too.** Check each one separately | E-111 — the lint job runs `cargo fmt --all --check` from the root, which excludes `bevy_isomesh`, so that crate's formatting had never been checked in 20 tickets and an example was committed unformatted. `cargo check`/`clippy`/`test` had their own steps in the bevy job; `fmt` was the one nobody noticed was missing |
| **Count what the claim is about, not what is easy to count.** A bitwise "did it move" is not a measure of "was it wrong" | A-009 — the clamp box is the cell shrunk by `ε`, so a vertex sitting exactly *on* a cell face is nudged by `5e-5` cells and a bitwise count calls it clamped. On a grid-aligned `box_exact` that read as **1176 of 1352 vertices displaced**; the honest count, against a `1e-3`-cell threshold, is **zero**. The first number would have gone into a commit message and been wrong |
| **A remedy stated for one operation does not cover the pipeline.** If a property is claimed end-to-end, check every reduction in it | M-24 — the audit's "magnitude-sorted dot products" is real and insufficient; the determinant needed the same treatment, and nothing said so. The equivariance test caught it because it asserted bit-equality rather than a tolerance |
| **Read a dependency's API before believing what it is for.** Reputation is not a type signature | ✗16 — glam is "the" Rust math library and was written into four documents as A-007's dependency. It has no generic scalar, so it cannot serve a crate generic over `f32` and `f64` |
| Before believing a performance verdict, ask **how many machines it has run on.** One is a hypothesis | ✗14 — Surface Nets loses to Marching Cubes by 2.76× at 256³ on an Apple M5, and the mechanism is probably cache. That is a strong result and a weak generalisation until it runs somewhere else (O-11) |

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

Two things that make this file worth keeping rather than a chore:

**Re-tier rather than rewrite.** When an R becomes an M, leave the R text and add the measurement
below it. The gap between what was reported and what was measured is itself data.

**Record the ones we got right for the wrong reason.** ✗10 is in the falsified section even though
the outcome was fine, because the reasoning was wrong and the reasoning is what generalizes.
