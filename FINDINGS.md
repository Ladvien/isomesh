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
| M-21 | **Surface Nets is not `O(n³)` over this range; Marching Cubes is.** Surface Nets' fitted intercept is **negative** — `−3.13 ms` full sweep, `−7.32 ms` on the tail — which is physically impossible and is the signature of a curve convex in `n³`. `r² = 0.9899` against Marching Cubes' `0.99986`. Per-sample cost rises `9.0 → 13.19 ns` while Marching Cubes' falls and flattens | T-006. Cause unmeasured — see O-11. This is why ✗14's gap widens rather than staying constant | **Amended at O-11: the intercept-sign diagnostic is machine-specific; the per-sample rise it stands in for is not.** On Zen 3 *both* fitted intercepts come back negative and both are numerically negligible — −0.04% of the largest run either way, `r² = 0.99999` — so the sign diagnoses nothing there. What reproduces is the underlying effect: Surface Nets' per-sample cost rises **31%** across the sweep (37.4 → 49.1 ns) while Marching Cubes' *falls* 13% (15.2 → 13.2). Report the per-sample curve, not the intercept
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
| M-45 | **✗14 reproduces on a second machine and gets worse; its crossover does not reproduce at all.** Same sweep, same field, same commit, AMD Ryzen 9 5900X (Zen 3, x86-64, single thread) against the Apple M5. Per-sample cost, `16³ → 256³`: **Marching Cubes M5 24.99 → 4.78 ns, Marching Cubes Zen 3 15.18 → 13.19; Surface Nets M5 8.40 → 12.66, Surface Nets Zen 3 37.38 → 49.08.** So Surface Nets degrades on both — the effect is the algorithm's memory pattern, not one cache hierarchy — and `Surface Nets/Marching Cubes` at 256³ is **3.72× on Zen 3 against 2.65× on the M5**. But Surface Nets never wins on Zen 3, at any resolution: `2.46×` behind even at 16³. The M5 crossover exists only because Marching Cubes starts expensive there and converges, which Zen 3's Marching Cubes does not do because it is flat from the start | O-11, `cargo bench --bench resolution_sweep` on `big` at commit `d2ab82a`. Raw data committed as `docs/measurements/resolution_sweep-ryzen9-5900x.csv` beside the M5's. Also: the M5 is **2.76× faster than the Ryzen on Marching Cubes at 256³** (80.2 vs 221.4 ms) single-threaded, while the Ryzen is faster below about 32³ |
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
| M-66 | **On a sharp field the geometry and the field disagree by an angle that does not fall with resolution.** Mean angle between area-weighted face normals and the analytic gradient: `sphere` **3.25° → 2.16° → 1.08°** and `torus` **11.65° → 6.07° → 2.45°** at 17/33/65³, both falling; `box_exact`'s mean falls **13.55° → 6.73° → 3.34°** but its **worst is 35.796° at all three resolutions, identical to six figures**. Refining a grid does not soften a corner — the disagreement there is geometry, not discretisation, where on a smooth field it is discretisation and shrinks | A-012. So "which way does the surface face" and "which way does the field increase" are different questions wherever the surface has a crease, and that is the whole reason the strategy is selectable rather than fixed. Asserted as the mechanism (box worst resolution-invariant, sphere worst falling) rather than as a pinned constant, since the constant is a property of the box's corner angle and would move with a different field |
| M-67 | **A sign test cannot distinguish 95.6% of the configurations a tetrahedron can actually be in.** Over every edge-coordinate vector with counts 0–3 on a tet's six edges, **181** satisfy normal surface theory's two conditions (even sum, triangle inequality) and only **8** are *classic* — every edge carrying at most one crossing. The other **173** put two or more crossings on some edge, where a sign test reads the parity alone and returns the classic configuration with the same parities. And classic Marching Tetrahedra is exactly the 0/1 corner of this encoding: taking `eᵢⱼ = 1` where corner signs differ reproduces A-003's triangle count on **all 96 (tet, configuration) pairs**, as 48 corner cuts and 36 diagonal cuts | A-014a, from Baktash, Gillespie & Crane `10.48550/arXiv.2606.00454` §2. The paper's own framing is that marching tetrahedra *"reinvented a small piece of this story"* — this puts a number on how small. It is also the quantitative form of A-005's `thin_plate` result: a feature thinner than a cell does not exist to a sign-based method, and 173 of 181 is how much else does not either |
| M-68 | **`parry3d`'s constructor is not a validity check: the only mesh it refuses is one with no triangles.** `TriMesh::new` returns `Result`, and its documented failure is *"the index buffer is empty (at least one triangle is required)."* Measured: it accepts a single zero-area triangle (three collinear points) and it accepts a two-chunk mesh with an unwelded seam, both without complaint. What *does* check is `set_flags(TriMeshFlags::HALF_EDGE_TOPOLOGY)`, which builds the half-edge adjacency and returns `TopologyError` when the mesh cannot support one, and `TriMeshFlags::ORIENTED`, which needs a closed consistently-oriented surface to compute pseudo-normals from | G-005, parry3d 0.30.2. So a caller who treats "the constructor took it" as "the mesh is fine" has checked nothing, which is the gap `collider::ColliderReadiness` exists to fill. The carved acceptance case — `csg_difference` at 41³, 4,484 triangles — passes both flags |
| M-69 | **A chunk seam costs 72 boundary edges, and welding removes exactly those and nothing else.** Two adjacent 16-cell chunks of a torus, meshed independently and concatenated: **36 duplicate vertices and 180 boundary edges**. After a weld at `1e-4` cells: **0 duplicate vertices and 108 boundary edges.** The 108 that remain are the two-chunk slab's own outer border, which is legitimately open — the surface leaves through the sides — so the weld closed the 72 that were the seam and left the rest alone | G-005. A renderer draws the unwelded version correctly and a physics engine reads those 72 as a hole, which is the concrete form of G-005's ticket note *"a chunked collider must be welded first or parry sees a seam of unshared vertices."* M-46 measured the same seam at A-013 as 80 boundary edges and 40 duplicated vertices on a different chunk pair; the mechanism is the same and the count is a property of the pair |
| M-70 | **Field-derived LOD is exact, not approximate: a coarse sample position is bit-identical to the fine one it sits on.** Level `k` doubles the spacing `k` times, so a level-`k` sample at index `s` and a level-0 sample at index `2^k·s` must land on the same world point — and they do, **bit for bit**, over cell sizes `0.125`, `4/35`, `0.1` and `1/3` and levels 0–3. Doubling is exact in IEEE and so is doubling a small integer, so `(h·2^k)·s` and `h·(2^k·s)` are the same real rounded the same way | G-004. Asserted rather than argued, because **M-32 and M-49 both caught this crate assuming an algebraic identity IEEE did not honour** — `cell_of` round-trips 3 of 3 cell corners at `h = 0.125` and 1 of 3 at `h = 4/35`. This one holds at every spacing tried, including those two. It is the precondition A-011b rests on: no coordinate drift can open a crack at an LOD boundary before transition cells get a chance to be wrong |
| M-71 | **Cells fall by 8 per LOD level and triangles by 4 — and the 4 degrades exactly where the grid stops resolving the surface.** Unit sphere over a fixed world extent, levels 0–3: **262,144 / 32,768 / 4,096 / 512 cells** (exactly `8×` each, by construction) against **9,512 / 2,312 / 536 / 104 triangles**, ratios **4.114, 4.313, 5.154**. A surface is two-dimensional so its triangle count tracks `area / h²`, but that is a *continuum* claim: by level 3 the sphere is four cells across and a 104-triangle staircase is not approximating anything smoothly | G-004. So each level buys back `8×` the sampling work and only `4×` the rendering, which is the whole economics of LOD and the reason the ticket's own acceptance figure is about cells. The tight `3.8–4.6` bound is asserted only on the two steps where the premise holds, and the upward drift is asserted as a *direction* rather than waived as tolerance |
| M-72 | **A sub-cell feature does not vanish under coarsening — it aliases, which is worse.** `thin_plate` across LOD 0–3: **4,088 → 1,016 → 248 → 56** triangles, still 56 at `h = 0.5` where the plate is a fraction of a cell thick. The test was written asserting it would be *gone* by the coarsest level and that assertion is what failed. Marching Cubes samples **corners** and cuts **edges**, so whichever edges happen to straddle a thin slab still register a sign change and what comes back is a partial, holey remnant | G-004. **The contrast is the mechanism:** A-005 measured the same field returning **zero** triangles under greedy quads, which asks one question per cell *centre* and therefore misses it cleanly. For a streamed world the aliasing is the worse behaviour — a feature that vanishes at a known distance can be faded, one that disintegrates into a resolution-dependent scatter pops. It is also the cost M-67 quantified from the other side |

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

---

## Part 4 — Open questions

Each has the test that would settle it. **An open question with no proposed test is a wish.**

| # | Question | Settled by | Why it matters |
|---|---|---|---|
| O-1 | What fraction of cells actually change per brush stroke? | G-002 instrumentation; hash cell slabs, log per stroke | **Unpublished.** Ceiling on every incremental-repair idea in the opportunities doc |
| O-2 | Does clamping the QEF vertex to (1−ε) inside its cell eliminate self-intersections? | A-009: measure per 1,000 triangles, clamp on vs off, all seven fields | Decides whether guaranteed intersection-free extraction is free → whether runtime convex decomposition can stop failing |
| O-3 | Marching Cubes vs Surface Nets vs Dual Contouring vs MT — actual relative speed on one machine? | M-001 | The comparison does not exist (V-17). We'd have the only apples-to-apples measurement |
| O-4 | Do brush operations commute? | G-003: 8 ops × 40,320 orderings, count distinct results. Expect 1 | If not 1, the coordination-free multiplayer story dies — cheaply, before anything is built on it |
| O-5 | Do mesh shaders work on macOS/Metal? | GPU-007 capability probe | **Sources contradict:** wgpu's spec table lists MSL as *planned*; the tracking issue says the Metal HAL backend merged. Neither is trustworthy until probed |
| O-6 | What is amortized meshing cost per frame under continuous editing? | E-206 under a deliberately overloaded queue | The only number a game cares about, and no paper reports it |
| O-7 | What fraction of *our* pipeline is contouring vs everything else? | M-003 | V-4 says 54% for someone else's code with no physics. Ours is probably worse |
| O-8 | Does Dual Contouring's vertex placement need f64 in practice, or is f32 enough? | E-112, with the QEF condition number in the HUD | `M = AᵀA` squares the condition number. **Half answered by V-18**: the original paper measures f32 error ~1 on flat regions at 256³, and recommends f64. **Partially answered by M-23**: on extraction paths with no solve, `f64` costs only 8–10% of wall time, so precision is cheap where there is no QEF. Still open for the vertex solve itself, and for *our* fields at *our* resolutions — which sidesteps `AᵀA` entirely and may not degrade the same way |
| O-9 | How much does T-003's gradient-flow chord **over**-estimate distance at a concave seam? | A comparison against nearest-point search over a dense surface point cloud, or E-104 once Dual Contouring lands | The chord follows `∇f` to the zero set, which near `csg_difference`'s seam can land further away than the true nearest point. The bias direction is known and safe for a "below X" gate; the *magnitude* is not measured, and M-001's shootout column would inherit it. `csg_difference` measured forward `0.0833` at 33³ — how much of that is seam bias is unknown |
| O-10 | ~~What is Surface Nets' non-manifold **rate** as a function of feature thickness over `h`?~~ **First curve measured at A-010 (M-60)**, as the multi-sheet-cell rate: `gyroid` 3.13% → 2.05% → 0.53% and `fbm_terrain` 1.70% → 0.84% → 0.77% at 17³/25³/33³, and exactly zero on the other five fields at every resolution. Still open only as the *slab* sweep, which would give thickness-over-`h` directly rather than resolution-at-fixed-field | A-010 drove it to zero, which was the ticket's job; a sweep over a slab of shrinking thickness would give the parametrised form | M-15 established it is a resolution effect rather than a topology one, and M-4 has counts at two resolutions on two fields. It decides whether Surface Nets is usable at game resolutions or needs A-010 first |
| O-11 | **Why does the dual topology go superlinear in `n³` while Marching Cubes does not?** *(Half-answered at M-45: it is not one machine's cache hierarchy. Surface Nets degrades on Zen 3 too — 37.4 → 49.1 ns/sample — and the `Surface Nets/Marching Cubes` ratio is worse there than on the M5. What remains open is the mechanism, not whether the effect is real.)* | A profile or cache-miss counter at 192³ vs 256³. The cross-machine experiment is **done**; a second one would not add anything | The working-set hypothesis survives and is strengthened: Surface Nets gathers the four cells around each crossed edge with one stride `n²` cells apart, and that stride is architecture-independent, which is exactly the kind of cost that would reproduce across microarchitectures. Note both machines show a per-sample **spike at 128³** (M5 Surface Nets 9.35, Zen 3 Surface Nets 53.84 against 45.6 at 96³ and 47.3 at 192³) — a working-set effect at one specific grid size on two unrelated cache hierarchies, which is itself a clue nobody has followed |
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

**P-6 — splitting the vertex reduces self-intersection on the two multi-sheet fields.** Registered
before running, reasoning from M-29 (*"the residue is exactly A-010's problem"*): if the clamp removed
every placement-caused intersection and the remainder is caused by two sheets sharing a vertex, giving
them separate vertices should remove the remainder. **Falsified at A-010 (M-61): it goes up**, `gyroid`
3.118 → 5.669 and `fbm_terrain` 13.837 → 15.434. The premise was right and the inference was wrong —
two vertices in one cell is exactly what breaks the within-cell partition the clamp's guarantee rests
on. ✗2's competing figure (ODC measuring Manifold Dual Contouring at 100% self-intersecting) was on
the record the whole time and should have been weighted against M-29 before registering this.

---

## Part 5 — Method rules, and the failure that earned each

Rules with no incident behind them get ignored. These all have one.

| Rule | Earned from |
|---|---|
| A typed error at the call site is louder than an abort — make the invalid state unrepresentable where you can, report it where you can't, and never substitute a default | The no-panic rule, reconciled with "fail loudly": `ValidateConfig` has private fields and one checked constructor, so the validator needs no runtime guard at all |
| Corpus presence is decided by `catalog_read`, never by `distill_search` | ✗4 — 342 documents readable but unsearchable |
| **Search home-still before believing a doc that says a paper is missing** | M-63 — `docs/research/` lists Manifold Dual Contouring and the Transvoxel dissertation as "genuinely absent, blocking"; both are in the corpus and both were read in one session. A-010's ticket named the wrong paper for its own algorithm because nobody had opened the right one |
| **Assert the property you believe, not the one that is easy — and when the assertion fails, the counterexample is usually the finding** | M-64, and M-59 before it. "A lateral link always crosses the resolution boundary" was false, and the case that broke it is precisely what transition cells exist to do. "Manifold Dual Contouring is manifold" was too strong, and the case that broke it is a second mechanism nobody had named |
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

Two things that make this file worth keeping rather than a chore:

**Re-tier rather than rewrite.** When an R becomes an M, leave the R text and add the measurement
below it. The gap between what was reported and what was measured is itself data.

**Record the ones we got right for the wrong reason.** ✗10 is in the falsified section even though
the outcome was fine, because the reasoning was wrong and the reasoning is what generalizes.
