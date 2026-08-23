# SOTA speed and feature frontier — and the harnesses that would falsify it

**Date:** 2026-08-16
**Scope:** speed advances and capability advances for isosurface meshing in games, from the local
corpus and from 2023–2026 external literature, each stated as something a testing harness could
refute. Plus the harnesses themselves — which to adopt, and which do not exist anywhere.

**How to read this.** Everything here is **tier R until you check it**. Four parallel sweeps produced
it (local corpus × speed, local corpus × quality, external 2023–2026, methodology), and each sweep
marked what it *read* versus what it saw in an abstract. Those marks are preserved. Where a claim
touches something `FINDINGS.md` already measured, the M-row is cited and **the M-row wins** — three
candidates below are demoted for exactly that reason, and one external claim is flat contradicted by
a measurement this repo already owns (§4.1).

Every DOI and arXiv ID below was returned by a tool during the sweeps. None was reconstructed from
memory. Where a source has no identifier, it has a URL and says so.

---

# Part 0 — The three cross-checks that reorder everything

Before any candidate: three numbers already in `FINDINGS.md` decide which half of the literature is
worth reading, and they point somewhere other than the extractor.

## 0.1 The extractor is 29% of the job, and the field is most of the rest

**M-135** measured the stage split on a usable mesh: **contour 29.0%, weld 25.5%, collider check
45.0%, normals 0.4%**. **M-136** measured the contour's own share ranging **13.1% → 74.3%** across
seven fields, with `fbm_terrain` at 65–74% *because fBm noise is expensive per sample*.

Compose those and the picture is stark:

| Where the time goes | Share | Who is optimising it |
|---|---|---|
| Field evaluation | dominant term on procedural fields (M-136) | **nobody, in this repo** |
| Collider readiness check | **45%** of a usable mesh (M-135) | nobody |
| Weld | 25.5% | nobody |
| The extraction algorithm itself | **≤29%**, ceiling **1.41×** if driven to zero | every paper in Part 1 |

**Consequence.** The 2003 SIMD-restructuring result (§1.5), the HistoPyramid-vs-scan question (§1.7),
the table-driven DC reformulation — all of them are optimising a term with a 1.41× ceiling on this
repo's own measurement. They are not wrong, they are **capped**. The candidates that attack field
evaluation (§1.1, §1.2) and the ones that attack the 45% collider line (§1.3) have no such cap.

**M-155 says the same thing from the GPU side**: cumulative 15.01 ms → 0.54 ms at 129³ across three
tickets, and `count + emit` — the actual Marching Cubes — went 0.11 ms → 0.04 ms. **None of the 28×
came from the extractor.** M-167 states it outright: synchronisation was 83% of an extraction and the
payload 7%.

## 0.2 Published constants are scene properties — this repo has proven it five times

✗14 (Surface Nets' crossover), M-51 (MT's 2–3× is really 2.87–3.91×), M-56 (greedy meshing's 2.76×
is really 1.70×–256×), M-60 (Nielson's 1.3% is really 0–3.13%), M-40 (the ambiguous face is absent
from five of seven fields). Five for five.

Every quality figure in Part 2 comes from **medical volumes or scanned CAD**. Two of the sweeps
flagged this independently and unprompted: the sliver rates (10.7–13.5%), DC's off-surface deviation
(39.09 vs 0.00), CoACD's 4.6× component reduction, the power-diagram Chamfer table — all measured on
data whose corner-value distribution has nothing to do with a procedural SDF. **Do not port a single
one as a baseline.** Every prediction in Part 2 is therefore written as *"this constant will not
reproduce, and here is what will"*.

## 0.3 Nobody publishes a comparable benchmark, which is a hole this repo is already standing in

Three independent findings, from three sweeps:

- The 2024 isocontouring survey (`10.3390/a17020083`) and MCPro (`10.5220/0013309800003912`)
  publish **no comparative benchmark**.
- **CGAL 6.1's Isosurfacing package** (shipped Oct 2025, TMC + DC, Cartesian + octree) ships with
  **zero timing numbers** in both the announcement and the manual.
- **No public voxel/SDF benchmark corpus exists.** The neural line (FlexiCubes, TetWeave, NDC) each
  rolls its own from Thingi10K/ABC via `SDFGen`, and FlexiCubes states in its own paper that its
  Table 5 numbers cannot be compared against its Section 5 numbers.

`docs/measurements/resolution_sweep.csv` plus the golden hashes plus M-31's cross-platform
bit-identity is, as far as four sweeps could find, **the only artifact of its kind in public**. That
is Part 3's last section, and it is the one item here that is a paper rather than a ticket.

---

# Part 1 — Speed candidates, ranked

Ranked by (measured evidence × transferability) ÷ cost, **after** applying §0.1's cap.

---

## 1.1 — Conservative-SDF empty-space skipping: cut evaluations, not instructions

**CLAIM.** A true signed distance field proves its own emptiness: `|f(p)| = d` means no surface within
radius `d` of `p`, so a cell whose centre sample exceeds its own half-diagonal can be skipped without
evaluating a single corner — and the output is **bit-identical**, not approximate.

**SOURCE.** The idea is stated as a by-product in Kohlbrenner & Alexa, *Isosurface Extraction for
Signed Distance Functions using Power Diagrams*, CGF 2025, DOI `10.1111/cgf.70037` (ABSTRACT only),
whose power diagram is a spherical-cap envelope of exactly this. Also the premise of Schott et al.,
*Sphere Carving: Bounding Volumes for Signed Distance Fields*, TOG/SIGGRAPH 2025, DOI
`10.1145/3730845` (ABSTRACT only), which claims acceleration of "Sphere Tracing **and
polygonization**" with no figure in the abstract. The complementary hierarchical form — a per-brick
`(min, max)` span test — is measured in Schmitz et al., CGF 2010, DOI
`10.1111/j.1467-8659.2010.01825.x` (READ), whose §5 states verbatim: *"the bottleneck of both the
HistoPyramids and Span Space data structures is the number of triangles in the isosurface, and not the
size of the volume."* NOISE's `O(√n + k)` search bound is Livnat, Shen & Johnson, DOI
`10.1109/2945.489388` (abstract READ).

**WHY IT OUTRANKS EVERYTHING ELSE HERE.** It attacks the **dominant** term (§0.1), it is bounded by a
*correctness invariant* rather than a heuristic, and its acceptance test is one this repo already
owns: **the golden hashes must not change by one bit**. There is no accuracy/speed trade to argue
about. It also composes with `SurfaceGate` and the `conservative` property some fields have and others
(`fbm_terrain`, `gyroid`) do not — which is itself the interesting measurement.

**COST.** A `Sdf::conservative() -> bool` (or a `Lipschitz` marker) on the field trait, plus a
hierarchical rejection pass. No dependency. The honest complication is that **two of the seven
reference fields are not conservative**: `fbm_terrain` is fBm and `gyroid` is trig — neither is
1-Lipschitz, so the skip is unsound for them without a stated Lipschitz constant. That is not a
blocker, it is the API question: a field that can supply `L` gets the skip, one that cannot gets the
brick span test instead.

**PRE-REGISTERED PREDICTION.** On `sphere` at 128³ a conservative skip reduces field evaluations by
**> 60%** with **zero** golden-hash changes. On `fbm_terrain` at 128³ the conservative route is
unavailable and an 8³ min–max brick pass rejects **> 60% of bricks** and cuts total extraction time by
**≥ 1.8×**. On `gyroid` at 128³ the brick pass costs **≥ 5% and gains < 10%** — a dense triply-periodic
surface has no empty space to skip.

**FALSIFIED BY.** Any golden hash changing (the skip is unsound). Or `fbm_terrain`'s gain landing under
1.3×, which would mean the brick size is wrong rather than the idea.

**HARNESS.** Instrument the existing per-field benches to emit `bricks_total`, `bricks_active`,
`cells_visited`, `field_evaluations`, `triangles_emitted`. Then plot **time against
`triangles_emitted` across all seven fields**: Schmitz's claim predicts a straight line through the
origin with volume size falling out as a parameter. That single plot is the whole experiment and it is
a `FINDINGS.md` row either way.

---

## 1.2 — Accelerating the field itself: in-tree optimisation nodes

**CLAIM.** Embedding proxy nodes and continuous-LOD nodes *inside* an implicit-surface construction
tree, while preserving the Lipschitz property, speeds up SDF evaluation by up to three orders of
magnitude — and the speedup applies to polygonisation, not only to sphere tracing.

**SOURCE.** Hubert-Brierre, Guérin, Peytavie, Galin (LIRIS), *Accelerating Signed Distance Functions*,
CGF 44, Oct 2025, DOI `10.1111/cgf.70258` — **ABSTRACT only**. "Up to three orders of magnitude" is
the authors' phrasing with no hardware attached; treat as a best case on deep trees until read.
Same group's newest: *The PhaseTree: Multiphase Signed Distance Fields*, TOG/SIGGRAPH 2026, DOI
`10.1145/3811379` (ABSTRACT READ) — multiphase sphere tracing at *"less than a 25% runtime overhead
compared to single-phase Sphere Tracing"*. Adaptive-octree SDF with proved `C⁰`/`C¹` continuity across
nodes: Pujol & Chica, C&G 2023, DOI `10.1016/j.cag.2023.06.020`, code at `github.com/UPC-ViRVIG/SdfLib`
— reports **100,000 particles at 5.5 ms/frame vs 1,000 particles at 10.9 ms/frame** against a sparse
voxel hierarchy, i.e. a ~200× query-throughput claim on one scene.

**TRANSFER.** This is the same object as `BrushStack`. **M-50** measured the walk directly: cost per
re-meshed chunk **0.158 / 0.354 / 0.525 / 0.589 ms** for edit logs of 1–15 / 16–30 / 31–45 / 46–60 —
3.7× for 7× the log, and flattening. **M-138** measured the same shape on `PaintStack`: 2.33× for 40×
the log. So the walk is real, sub-linear, and nothing has tried to make a *node* out of it. An
in-tree proxy node is precisely "collapse a prefix of the edit log into a cheaper conservative
approximation over a region."

**PhaseTree is separately a scope question, not a ticket.** A multiphase SDF *is* a multi-material
voxel world (stone / dirt / ore). The public API is currently one scalar field; multi-material is an
API change and belongs to the user, not to a research doc.

**PRE-REGISTERED PREDICTION.** Two, and the first is the gate on the second. **(a)** On `fbm_terrain`
at 64³, **over 70% of Marching Cubes wall-clock is inside the field closure** — measurable in one
afternoon by running `cargo bench --bench extract` twice, once with the real field and once with an
`#[inline(never)]` constant-returning stub of identical signature. M-136's 65–74% is the same quantity
measured a different way, so this should confirm rather than surprise. **(b)** If (a) holds, a single
conservative proxy node collapsing a 60-brush log to one bounded primitive over its own AABB cuts
per-chunk re-mesh cost by **≥ 2×** for logs past 30 brushes, with mesh differences confined to the
brushes' own support.

**FALSIFIED BY.** (a) coming back under 40%, which would mean M-136 is `fbm_terrain`-specific in a way
nobody has stated. Or the proxy changing geometry outside its own AABB, which is a soundness failure
rather than a slow result.

**HARNESS.** The stub-field ratio bench for (a). For (b), reuse E-202's `ISOMESH_AUTOCARVE` trace and
report per-chunk cost against log length with and without collapsing, plus a byte-identity assertion
outside the collapsed region (§3.7's edit-locality instrument does this exactly).

---

## 1.3 — The 45% nobody is looking at: collider readiness

**CLAIM.** Validating that a mesh is collider-ready costs more than extracting it, and no paper in
four sweeps measures this at all.

**SOURCE.** **M-135**, this repo. `collider::readiness` is 45.0% of a usable mesh, against contouring's
29.0%. The external anchor is **M-68**: `parry3d`'s `TriMesh::new` refuses only an *empty* index
buffer, so the check has to exist; and **M-69**: a chunk seam costs 72 boundary edges, which a
renderer draws correctly and a physics engine reads as a hole.

**WHY IT IS HERE AT ALL.** Because §0.1 says the extractor has a 1.41× ceiling and this line is bigger
than the extractor. The literature is silent — every paper measures extraction, and the one
methodology paper that measures a *pipeline* (Grosso & Zint, V-3/`10.1007/s00371-021-02139-w`) measures
halfedge construction, not physics readiness.

**PRE-REGISTERED PREDICTION.** The 45% is dominated by the self-intersection census, not by the
manifold walk, because self-intersection is the only super-linear check in the set. Specifically:
`self_intersections` will be **> 60%** of `readiness` on `gyroid` at 33³ and the manifold/orientation
walks together **< 25%**. And **most chunks do not need the check at all** — a chunk whose extractor
carries a structural guarantee (Marching Cubes: M-53 measured 0 non-manifold and 0 self-intersections
across every field and both grids) can skip it, so the readiness check is a *dual-method* tax rather
than a pipeline tax.

**FALSIFIED BY.** The manifold walk dominating, which would mean the check is cheap to make incremental
and the whole framing is wrong.

**HARNESS.** Break `readiness` into per-check timings in `stage_breakdown.csv`, same shape as M-135's
existing rows. Then a second column: `readiness` cost with the extractor's own guarantees consulted
first. If Marching Cubes chunks can skip the census outright, the 45% collapses for the default
extractor and stays for dual methods — which is a finding about the *pipeline*, not about a check.

---

## 1.4 — Decoupled Fallback: the prefix scan that does not crash on Apple

**CLAIM.** Single-pass decoupled-lookback scan needs forward-progress guarantees that Apple M-series
and ARM Mali do not provide, and on those devices it TDR-timeouts or runs *slower* than the two-pass
method. A work-stealing "fallback" variant restores single-pass performance with no 64-bit atomics and
no explicit memory barriers, at 98–104% of memcpy across six GPUs.

**SOURCE.** Smith, Levien & Owens, *Decoupled Fallback: A Portable Single-Pass GPU Scan*, SPAA '25,
DOI `10.1145/3694906.3743326` — **READ in full**, open access at
`escholarship.org/content/qt0bk9z4bt/qt0bk9z4bt.pdf`. Reference implementations across CUDA / D3D12 /
Unity / **WGPU** at `github.com/b0nes164/GPUPrefixSums`.

**NUMBERS** (2²⁵-element 32-bit inclusive prefix sum, 10⁹ elements/s):

| Device | forward-progress guarantee | memcpy | Decoupled Fallback | Reduce-then-Scan | DF/RTS |
|---|---|---|---|---|---|
| Apple M1 Max 32-core | **no** | 37.42 | 36.85 | 25.75 | **1.43×** |
| Apple M3 10-core | **no** | 10.82 | 10.87 | 7.46 | 1.46× |
| ARM Mali-G78 MP20 | **no** | 2.626 | 2.725 | 2.025 | 1.35× |
| RTX 2080 Super | yes | 50.95 | 51.25 | 34.36 | 1.49× |
| RX 7900 XT | yes | 73.68 | 83.27 | 62.57 | 1.33× |

On M1 Max, DF completes in ~1 ms against a 2000 ms TDR threshold, where **decoupled lookback
frequently triggers timeouts or exceeds 100 ms**. The paper also documents that WebGPU and Metal lack
both 64-bit atomics and explicit memory barriers, which is *why* lookback cannot be ported.

**TRANSFER — and the demotion.** GPU-010a already moved the scan onto the GPU and measured
**5.24 ms → 0.37 ms**, a factor of 14, and **M-155** then took the whole 129³ path to **0.54 ms**. So
the scan is no longer a meaningful share of anything: a further 1.4× on a stage that is already a small
part of 0.54 ms is worth microseconds. **This candidate is therefore ranked here for the trap, not for
the win.** The load-bearing content is: *if anyone ever ports a CUDA lookback scan into `isomesh-gpu`,
it will hang on the M5.* Two of this project's machines are an M5 Mac and a Linux/RTX 3090 box, and
the crash is device-class-specific, so it would pass CI on the Ryzen and die on the laptop.

**PRE-REGISTERED PREDICTION.** A naive decoupled-lookback scan in WGSL through wgpu will, on the M5,
either fail naga validation (no 64-bit atomics) or hang at input sizes ≥ 2²⁰ — while completing
normally on the RTX 3090. **Do not run this to find out**; the paper documents it as a device crash.

**HARNESS.** None needed. This is a note in `CLAUDE.md`'s GPU section and a link in the scan module's
doc comment. If the current scan is ever benched against DF, the falsifier is DF/RTS landing outside
1.3–1.5× on the 3090, which would mean the port is wrong.

---

## 1.5 — Stage-streamed extraction (Flying Edges' ancestor)

**CLAIM.** Restructuring Marching Cubes from "loop over cells, do everything per cell" into four
streaming stages — sign classify → edge gather+interpolate → topology lookup → triangle assembly —
gives ~4× on CPU, and the cell-by-cell formulation **cannot be vectorised at all**.

**SOURCE.** Newman, Byrd, Emani, Narayanan, Dastmalchi, *High performance SIMD marching cubes isosurface
extraction on commodity computers*, C&G 2003, DOI `10.1016/j.cag.2003.12.008` — **READ in full**.
Measured on Pentium III 450 / Athlon 750 against compiler-optimised C: comparison stage **3.0× / 5.0×**
(~2.33–2.67 clocks/element, ~50% of peak memory bandwidth); edge interpolation **2–5× / 2.5–8×**;
topology **2.5–3.0× / 3.25–4.0×**; triangle construction only **1.03–1.53× / 1.09–2.32×**, and higher
*only when the active-cell fraction is low*. Prefetch worth >25%.

**THE GAP THAT MATTERS MORE THAN THE PAPER.** **Flying Edges** (Schroeder, Maynard, Geveci, LDAV 2015)
— the shared-memory-parallel Marching Cubes that VTK ships, and the conclusion of exactly this
restructuring — is **not in the corpus**. Presence was checked against the full 9,354-row
`catalog_list` dump, not by a failed search (✗4's rule). It is the single highest-value acquisition in
this document.

**TRANSFER, capped.** §0.1 caps this at 1.41× on the whole job, and `core::simd` is unstable so an
explicit-SIMD version needs a dependency (rule 3). What survives the cap and needs no dependency:
**(a)** autovectorisation of the classify pass over a `&[R]` slice, which is where the measured 3–5×
lives; **(b)** an edge-indexed intersection array *is* the vertex dedup, with no hash map — and A-014g
already established (M-168) that giving a crossing an identity rather than a position removes 5.01× of
the subgrid extractor's vertices with **zero** triangle-count changes, which is the same idea one
algorithm over; **(c)** it is the decomposition the GPU path already has, so CPU and GPU stop being
two algorithms.

**PRE-REGISTERED PREDICTION.** On `sphere` at 64³, stage-splitting makes total wall time **≥ 1.6×**
faster than the per-cell loop in release on arm64, and the classify pass alone runs at **< 4 cycles per
grid point**. If it comes out *slower*, the cause is the intermediate arrays' memory traffic — each of
the paper's stages sat at ~50% of peak bandwidth — which is the measurable thing to check rather than
a mystery.

**HARNESS.** A `stage_streamed` variant behind a bench-only feature; `cargo bench --bench extract` on
all seven fields; refit `t = a + b·n³` and compare the `b` coefficient, which the sweep already
produces (M-62).

---

## 1.6 — Workgroup-shared staging for GPU reductions

**CLAIM.** In a GPU tree/scan pipeline the *reduction*, not the geometry, is the bottleneck, and staging
its writes through workgroup shared memory plus `countbits` over packed 32-bit words gives an order of
magnitude on that pass.

**SOURCE.** Deliot, Yao, Dupuy, Rijnen (Unity), *Experimenting with Concurrent Binary Trees: Large Scale
Terrain Rendering*, SIGGRAPH 2021 Advances — speaker notes **READ**. GTX 1080 Ti, CBT depth 28, 2 km²
terrain, 1,409,792 triangles: **subdivision 0.03 ms; first four sum-reduction passes 5.78 ms → 0.40 ms
after optimisation, "14.5x faster"**. Nsight attribution: Long Scoreboard (buffer-access) stalls before,
Short Scoreboard (shared memory) after. Note the underlying CBT paper, Dupuy 2020 DOI `10.1145/3406186`,
**is in the corpus as a HAL landing page only** — abstract and metadata, no body, no tables.

**Ledger note.** This is ✗7/V-10 already: the ledger has the 5.78 → 0.40 figure with the attribution
corrected to the Unity talk. What is *new* here is the mechanism — shared-memory write batching (one
`u32` per 16 threads) plus `countOneBits` over the 1-bit leaf level to skip reduction levels 2–5 — all
three expressible in WGSL today, none needing mesh shaders.

**Demoted for the same reason as §1.4.** GPU-010a's scan already landed; M-150 measured the stage at
0.37 ms and M-155 took the path to 0.54 ms total. The remaining value is a checklist for the *next*
GPU reduction anyone writes here, not a ticket.

---

## 1.7 — Output-driven (HistoPyramid) vs scatter (scan): a crossover, not a winner

**CLAIM.** For sparse volumes, output-driven pyramid traversal beats scatter-based prefix scan despite
its `O(log n)` per-output gather, because cache locality pays for the depth; for dense volumes the scan
wins.

**SOURCE.** Dyken, Ziegler, Theobalt, Seidel, *High-speed Marching Cubes using HistoPyramids*, CGF 27
(2008) 2028–2039, DOI `10.1111/j.1467-8659.2008.01182.x` — **body READ, Table 1 unusable**: the OCR
duplicates column headers and corrupts dataset dimensions ("255 × 255 × 25"). Not one number from that
table is trustworthy. What is clean is the analysis: *"Scan's output extraction iterates over all input
elements and scatters the relevant ones to output, while HP iterates on the output elements instead… if
a lot of the input elements are to be culled, which is the case with MC for larger and sparse volumes,
the HP algorithms can play out their strengths, despite the deep gathering traversal."* Also: the OpenGL
vertex-shader variant beat the geometry-shader variant even on DX10 hardware.

**PRE-REGISTERED PREDICTION.** On a 128³ chunk with **< 3% active cells**, output-driven traversal beats
scan-and-scatter by **≥ 1.5×**; at **> 25% active cells** scan wins by ≥ 1.2×. There is a crossover and
locating it is the deliverable. Caveat stated up front: this is 2008 hardware and the table was
unreadable, so tier R at best.

**HARNESS.** Two kernels behind one host API; sweep the isovalue on `gyroid` and `fbm_terrain` to move
active-cell density from 1% to 40%; plot the crossover.

---

## 1.8 — Layout: Morton across chunks, never Hilbert, row-major inside

**CLAIM.** Morton (Z-order) improves memory-bound throughput over row-major once the working set exceeds
cache; Hilbert does **not** pay for itself — its locality advantage is ~1.6% while its index computation
is an order of magnitude more expensive.

**SOURCE.** Reissmann, Meyer, Jahre, *A Study of Energy and Locality Effects using Space-filling Curves*,
arXiv `1606.06133` — chunks READ. Dual-socket Intel, 40 MB aggregate LLC: Hilbert LLC data read misses
**16.78×10⁶** vs Morton **17.06×10⁶** (~1.6% better) while Hilbert's *"absolute execution time is an
order of magnitude higher than RM and MO"*. Row-major wins outright for in-cache sizes: *"the order of
data element access is practically insignificant."* Dyken independently chose Morton for the same reason.

**PRE-REGISTERED PREDICTION.** Switching a 32³ chunk's internal sample array to Morton changes extraction
time by **< ±5%** (it fits in L2), and switching *inter-chunk* iteration to Morton improves a 16-chunk
re-mesh by **≥ 10%**. **The null result at 32³ is the finding worth recording**, because "Morton is
always better" is exactly the folklore this ledger exists to kill.

**HARNESS.** A/B both orderings in the existing sweep at 32³ and 256³.

---

## 1.9 — Fewer triangles at equal error is a speed result wearing a quality hat

**CLAIM.** Subgrid Marching Tetrahedra emits fewer polygons than classic marching **at equal
reconstruction error**, so at matched quality it reduces every downstream cost: upload, rasterise,
simplify, collide.

**SOURCE.** Baktash, Gillespie & Crane, *Subgrid Marching Tetrahedra*, TOG 45(4) art. 57 (SIGGRAPH 2026),
DOI `10.1145/3811358`, arXiv `2606.00454`. Fig. 21 caption **READ**: *"Our subgrid method generally
produces fewer polygons than classic marching for equal reconstruction error… on the DORA dataset."*
Project pages state **output at 180³ had error comparable to classical marching tets at 500³**.
**No runtime in ms and no hardware are published on either project page** — the "similar compute time"
claim is unverified anywhere public, which matters because **M-98** measured this implementation at
**70× classic MT and 196× Marching Cubes**.

**THE ACTIONABLE PART IS MEASUREMENT, NOT IMPLEMENTATION.** The algorithm is already in the crate. What
does not exist is the **error-vs-triangle-count curve**, which is the form that converts a quality result
into a speed result for a Bevy consumer. And **M-51** already built half of the instrument: it reproduced
M-10's `1.380e-3` exactly on the same harness.

**And the 70× has a named fix whose precondition is now met.** The 2026-08-14 round-2 doc identified the
grid-edge root cache and recorded that it was deferred because the correctness precondition was
unstated. **A-014g stated it and then went further**: M-168 gave every crossing a global identity keyed
on `(the tet edge's two grid points, the root's ordinal)`, and M-184 completed it by naming endpoint
roots after the grid point they sit on. **That key is the cache key.** The reason the cache was deferred
is gone — 576 field evaluations per cell against Marching Cubes' 8 shared corner samples is the largest
single measured cost in the crate, and it is now redundancy against a key that already exists.

**PRE-REGISTERED PREDICTION.** (a) At matched two-sided Hausdorff against analytic `torus` and `gyroid`,
subgrid emits **≥ 20% fewer** triangles than plain Marching Tetrahedra at the same error while costing
**< 2×** the extraction time — net positive end-to-end wherever the mesh is uploaded or simplified.
(b) A grid-edge root cache keyed on M-168's identity cuts subgrid's field evaluations by **≥ 5×** with
**zero golden-hash changes** — because a cached root and a re-found root are the same root through a
deterministic finder on identical endpoints (M-95's "same arguments, same output" boundary).

**FALSIFIED BY.** (b) changing any golden hash — which would mean the identity is not actually a function
of the grid edge and M-168's key is weaker than measured.

**HARNESS.** An `error_vs_triangles` bench: sweep resolution for MT and subgrid, compute Hausdorff against
the *analytic* field (no reference mesh needed — the fields are analytic), plot triangles against error.
Cheapest high-value experiment in this document, and it needs no new algorithm.

---

## 1.10 — Calibration anchors worth benchmarking against

Not techniques — targets. Rule 4 forbids a performance claim without an in-repo benchmark; these are
published numbers to be embarrassed by.

| Anchor | Number | Hardware | Source |
|---|---|---|---|
| GPU dual MC vs GPU MC | Skull 512²×641: **DMC 68 ms / MC 31 ms**; Pelvic 512×512×1047 → 3.36M verts in **92 ms** | RTX 2080 Ti, CUDA | Grosso & Zint, DOI `10.1007/s00371-021-02139-w` (tables READ) — this repo's V-3 |
| Adaptive tessellation, **WebGPU** | per 1M tris: **0.36 ms** single lit pass vs OpenSubdiv 0.50; **0.81 ms** vs 2.51 for lit + 4 shadow passes | RTX 3070 | Hable (Meta), SIGGRAPH 2026 Advances, slides READ — no DOI |
| Work-graph procedural geometry | **3.13 ms** to generate + G-buffer unique trees; 79,710 instances in **3.74 ms** | RX 7900 XTX | Kuth et al., HPG 2025 (no DOI found) and HPG 2024 DOI `10.1145/3675376` |
| Voxel raymarching at 10¹⁰ voxels | **6 ms/frame**, 2–4× faster than HashDAG above 32K res, **24 MB VRAM (2%)** | RTX 3060 Ti, Vulkan | Aokana, arXiv `2505.02017` / DOI `10.1145/3728299` — this repo's V-14 |
| Shipped triangle budgets | **4.7M desktop / 2.4M tablet / 1.1M phone** per frame after 92–98% reduction | — | Roblox SLIM, SIGGRAPH 2026 Advances slides READ |

**The Hable row is the most transferable number in this document** and nothing in the ledger cites it:
it is **WebGPU**, i.e. the same API surface class as wgpu, so it transfers far better than any D3D12
work-graph result. Its architectural argument — *tessellate once into flat buffers, consume across all
passes* — applies verbatim to chunk extraction: the advantage grows with pass count, which is the
0.50 → 2.51 ms blowup. Its memory warning is worth as much: 16M triangles + 16M vertices = **~704 MB**
of flat buffers, which is why `MeshBuffer` reuse is doing real work.

---

# Part 2 — Feature and quality candidates

Ranked by what they let a game *do* that it cannot now. Every published constant here is flagged as
scene-specific per §0.2.

---

## 2.1 — Intersection-free contouring by quad splitting — the named fix for M-61

**CLAIM.** Dual-method self-intersections can be **eliminated**, not merely reduced, by splitting the
dual quads under a cell-local rule — and the extension to *multiple vertices per cell* (Manifold Dual
Contouring's case) is published.

**SOURCE.** The technique is **Ju & Udeshi, "Intersection-free contouring on an octree grid", Pacific
Graphics 2006**, which is **cited three times in the corpus and present zero times** (checked against
the full stem index). The extension is Hwang & Sung, *Occupancy-Based Dual Contouring*, 2024, arXiv
`2409.13418`, **which is in the corpus and indexed**, stating: *"We incorporate the quad splitting
technique from Intersection-Free Contouring (IC) [Ju and Udeshi 2006] to resolve this issue. IC was
originally designed to handle a unique 3D point within each cell; we extend this idea to accommodate
multiple 3D points per cell."* The detail is in their appendix §A.1.5. Independent industry
corroboration: a SIGGRAPH 2015 *Advances in Real-Time Rendering* deck (corpus stem
`s2015-advances-pdf-230-mb`, indexed) whose speaker notes read *"ALSO! Oh no, there are self
intersections! This makes the lighting look glitched - fix em: … 'Intersection-free Contouring on An
Octree Grid' Tao Ju, Tushar Udeshi."*

**WHY THIS IS THE STRONGEST FEATURE CANDIDATE HERE.** ✗2 already cites ODC's measurement — Manifold Dual
Contouring at **100% of models self-intersecting** against ODC at **0 of 1500** — and this repo then
reproduced the direction independently at **M-61**: splitting the vertex makes self-intersection
*worse*, `gyroid` 3.118 → 5.669 and `fbm_terrain` 13.837 → 15.434, because A-009's cell-clamp partition
argument assumes one vertex per cell. **So the ledger already contains both the problem and the
citation of the paper that solves it, and nobody has extracted the mechanism.** M-29's branch rule said
the residue after clamping is a *connectivity* problem; quad splitting is a connectivity fix.

The industry note is the sharpest part: the shipping symptom was **broken lighting**, not broken
geometry — which means it is invisible to every check this repo currently runs except the
self-intersection census.

**COST.** Splitting quads raises triangle count and breaks the "one vertex per cell, one quad per
sign-change edge" invariant that greedy quad meshing (A-005) depends on.

**PRE-REGISTERED PREDICTION.** Applying a quad-split rule to Dual Contouring on `gyroid` at 33³ takes
self-intersections from **3.118 per 1k to exactly 0**, and on `fbm_terrain` from 13.837 to exactly 0,
while raising triangle count by **< 15%** and leaving symmetric Hausdorff unchanged to within 1%. On
Manifold Dual Contouring the same rule takes M-61's 5.669 and 15.434 to 0 as well — and if it does
*not*, the residue is M-225's four-quads-on-one-dual-edge mechanism, which is a different defect and
needs a different fix.

**FALSIFIED BY.** Anything other than exactly zero on a field with no multi-sheet cells. Or the triangle
count rising more than 15%, which would make it a cost question rather than a free guarantee.

**HARNESS.** Mostly built. Run the existing census across seven fields × {MC, DC, MDC} × three
resolutions and **publish the baseline table first** — M-53 has most of it already. Then compare against
the two published four-property tables (TetWeave's 2-manifold / sharp / intersection-free / fair;
FlexiCubes' five-column taxonomy) so the numbers are comparable rather than isolated. **Warning
flagged by the sweep: FlexiCubes' Table 1 taxonomy is column-scrambled in the corpus markdown — do not
read values off it.**

**ACQUIRE.** Ju & Udeshi 2006. It is the load-bearing citation for a guarantee this crate does not hold
and the one paper the corpus can only see second-hand.

---

## 2.2 — Method of Manufactured Solutions: a test *class* this repo does not have

**CLAIM.** Isosurface extractors have a *formal order of accuracy* per property, and measuring observed
order against it catches implementation bugs that manifoldness, orientation and χ checks pass silently.

**SOURCE.** Etiene, Scheidegger, Nonato, Kirby, Silva, *Verifiable Visualization for Isosurface
Extraction*, TVCG 2009, DOI `10.1109/tvcg.2009.194` — **in corpus, READ in full**.

| Property | Error definition | Formal order |
|---|---|---|
| Vertex position | `max_j \|λ − f(v_j)\|` — algebraic distance, **L∞**, measured at vertices | **O(h²)** |
| Normals | `max_j \|θ\|` between face normal and true normal at the point nearest the **centroid**, L∞ | **O(h)** |
| Area | `\|A(S) − A(S̃)\|` | **no formal order exists**; ≈2.0 observed empirically |
| Curvature | angle-deficit Gaussian, L∞ | **O(1) — predicted NOT to converge** |

Observed (sphere on `[−4,4]³`): VTK MC **1.94 / 0.93 / 2.00 / −3.35**; **buggy** DC **1.02 / −0.11 /
0.69 / −2.08**; the same code fixed **1.96 / 0.96 / 1.89 / −0.15**. The bug was *"a hard-coded limit in
the number of steps in a root-finding procedure that was being triggered by the high resolution of the
volume"* — invisible at coarse `h`, and the authors state plainly that *"both faulty implementations
performed appropriately for large values of h."*

**WHAT THIS REPO HAS AND WHAT IT LACKS.** M-12 measured Marching Cubes' `h²` position convergence
(ratio 4.179 against an ideal 4.13) and M-65 measured `h²` normal-direction convergence for central
differences. What is missing is four specific things the paper says are load-bearing: **the L∞ norm**
(L2/mean wash out off-by-one grid sizes and node-vs-cell-centric confusion); **measurement-location
discipline** (vertices for algebraic distance, centroids for normals — a harness that mixes them reports
nonsense); **a wide `h` range** (they swept `h ∈ (0.001, 1.0)`; three or four resolutions is not enough);
and **negative controls** (curvature predicted `O(1)`, a linear field predicted `p = 0` — without which
a passing suite proves less than it looks like).

**And the fixture warning is this repo's own rule arriving from outside.** A sphere gave the *buggy* DC
a passing `p = 0` for a linear field; Afront's spline reproduced a sphere to numerical error, giving a
spurious `p ≈ 0`. They had to switch to `x²+y²+z² + cos(Ax)² + cos(Ay)² + cos(Az)²` before the real
deficiency appeared. **`gyroid` and `fbm_terrain` are the right instrument; `sphere` alone is a
known-weak oracle** — the fixture trap, published, eight years before this repo hit it eight times.

**PRE-REGISTERED PREDICTION, and it is a sharp one.** On `sphere` over `h ∈ [1/256, 1/16]`, Marching
Cubes fits vertex slope **1.8–2.1** and normal slope **0.8–1.1**. **Dual Contouring with the Tikhonov
regulariser at fixed λ = 0.01 fits vertex slope < 1.5 and normal slope < 0.5** — because a fixed λ is an
`O(1)` bias that does not vanish under refinement, so DC cannot be second-order unless λ is annealed
with `h`. ✗12 derived λ = 0.01 as the value reproducing DC's σ = 0.1 truncation smoothly; **nothing has
ever asked whether that constant is allowed to be constant.**

**FALSIFIED BY.** DC fitting 1.8–2.1 at fixed λ, which would mean the regulariser's bias is below the
discretisation error over the whole range and the concern is empty. Either outcome is a `FINDINGS.md`
row, and the second is more interesting.

**HARNESS.** `benches/convergence.rs`: per analytic field with a known closest-point map, extract at
`n ∈ {16, 32, 64, 128, 256}`, report `max|f(v)|`, `max` centroid normal angle, `|A_exact − A_mesh|`,
fit three log-log slopes per (algorithm, field), report R². ~200 LOC, no dependencies, belongs in a
bench tier not `cargo test`. **The expensive part is not the code — it is deriving or citing the formal
order per (field, property, algorithm), which the authors say took them longer than the tests.**

---

## 2.3 — An analytic χ oracle for fields where no χ is derivable

**CLAIM.** Stratified Morse theory gives the Euler characteristic the *correct* isosurface must have,
computed analytically from a random trilinear field's critical points — so χ stops being something you
record and becomes something you predict.

**SOURCE.** Etiene, Nonato, Scheidegger, Tierny, Peters, Pascucci, Kirby, Silva, *Topology Verification
for Isosurface Extraction*, DOI `10.1109/tvcg.2011.109` (identifier from `paper_search`; tech-report
version at `julien-tierny.github.io/stuff/papers/etiene_techrep10.pdf`). Three instruments: a
**consistency test** (random fields, interior vertex links must be circles, boundary links lines — the
link walk this repo already does, driven by fuzzing); **MMS–SMT**, which partitions each cell into
strata and sums critical-point contributions to get expected χ, and **works for surfaces with
boundary**, hence on individual chunks; and **MMS–DS**, a digital-topology Betti oracle needing closed
surfaces and grid refinement to unambiguity.

**What these found in real code**, useful as regression targets: MACET failed to traverse boundary cells
→ non-manifold output. **MC33 (Chernyaev): a coding error in configuration 13, and an orientation
problem in case 13.5.2 requiring a criterion not documented in the original publications.** DelIso: 28
cases with holes, 15 with missing triangles, 7 with duplicated vertices. SnapMC with non-zero snap:
**> 50% failure rate**. MATLAB `isosurface`: non-manifold edges. **MCFLOW: 0% mismatch — the only
passing code, and the only one that used these procedures during development.**

**Why this lands here specifically.** This repo pins observed χ in golden fixtures for `gyroid` and
`fbm_terrain` because no χ is derivable a priori. SMT derives it a priori for **random trilinear
fields** — a class not currently tested at all, and strictly harder than seven analytic fields. It is
also the exact instrument the A-002 series wanted: M-208 measured that **0 of 68,385 reference-field
surface cells reach six body saddles**, so every gate the series could lean on was structurally
incapable of exercising the code it exists to write.

**PRE-REGISTERED PREDICTION.** Over 1,000 random trilinear cells, this repo's Marching Cubes agrees with
the SMT-derived χ in **100%** of cases, and the MC33 decider agrees in **100%** — but the trilinear
interior rule (A-002b) disagrees on the `[9,3]` case-13 configurations M-228/M-229 already found, which
are precisely Chernyaev's configuration 13 and precisely what Etiene's team found broken.

**FALSIFIED BY.** Any disagreement outside case 13, which would be a new defect rather than a known one.

**COST AND RISK.** ~400–600 LOC. **This is the highest transcription risk in this entire document.**
Eq. 4.7 must be read from the paper, never reconstructed — rule 5, and M-219's warning that even a
reference implementation can carry a one-line typo.

---

## 2.4 — Isotopy, not just manifoldness: interval arithmetic as a topology gate

**CLAIM.** For an analytically-evaluable field, an octree whose subdivision criterion is an *interval*
test on `∇F` yields a mesh **regularly isotopic** to the true surface — same topology *and* same
embedding, no self-intersections — with no critical-point computation.

**SOURCE.** Plantinga & Vegter, *Isotopic approximation of implicit curves and surfaces*, DOI
`10.1145/1057432.1057465` — in corpus, indexed. Companion: Boissonnat, Cohen-Steiner, Vegter, DOI
`10.1007/s00454-007-9011-4` — *"the first certified algorithm for the more difficult problem of isotopic
implicit surface polygonization… our output can be continuously deformed into the actual implicit
surface without introducing self-intersections."*

**The authors' own honesty is what makes this useful.** *"The current implementation does not try to
create a good approximation in terms of Hausdorff distance. An open problem remains how to improve the
mesh quality without loosing isotopy."* And it gives **no sharp features** — a direct conflict with the
DC path.

**SO DO NOT BUILD THE OCTREE. BUILD THE MEASUREMENT.** Implement `Interval<f64>` for add/mul/sin/cos/sqrt
(~200 lines, no dependency), extend `sphere` and `gyroid`, and on a *fixed uniform grid* report the
fraction of surface-crossing cells where the interval test says "not parametrisable" — i.e. **the
percentage of cells this crate's output resolves by fiat**. That is a number nobody has for a game
mesher.

**PRE-REGISTERED PREDICTION.** On `sphere` the unresolved fraction at 33³ is **< 0.1%**; on `gyroid` it is
**2–5%**; on `fbm_terrain` it is **> 5%** and does not fall as fast as `1/n` under refinement, because
fBm has structure at every scale (M-108 measured exactly this scale-invariance for the self-intersection
residue: `fbm_terrain` wandered 25.40 → 20.43 across nine resolutions while `gyroid` fell 6.4×).

**FALSIFIED BY.** The unresolved fraction being near zero on all seven fields, which would mean the idea
buys nothing here and should be closed with a reason.

**Second, sharper harness — the genus-stability sweep.** Sweep the grid origin over 20 sub-voxel offsets
on `gyroid` at 64³ and count distinct observed χ. Uniform MC should produce **≥ 3** distinct values; an
isotopy-certified method exactly **1**. This needs no interval arithmetic at all and is a day's work.

---

## 2.5 — Off-surface deviation: two lines that measure the Tikhonov bias

**CLAIM.** Dual Contouring vertices deviate from the isovalue by one to two orders of magnitude more than
MC/TMC/DMC vertices, which sit exactly on the trilinear interpolant.

**SOURCE.** Grosso & Zint, DOI `10.1007/s00371-021-02139-w`, Table 6 average `Δf(v)`: **Skull — DMC 0.00,
DMCS 4.34, DC 39.09, TMC 0.00, MC 0.00**; Body DC 5.26; Skeleton 4.41; iWP DC 3.44 vs DMC 0.97. Their
explanation is structural: DC uses a quadratic error function, so its vertices are *in general not on
the isosurface*. **And the paper contains its own scene-dependence control**: *"For more smooth surfaces
like gen2 all methods behave very similar"* — gen2 shows DC at 0.68–0.70 against MC 0.69–0.70, i.e. **no
difference at all**.

**PRE-REGISTERED PREDICTION** (§0.2 applied): on `sphere`, `torus`, `gyroid`, mean `|f(v)|` for DC and MC
differ by **< 2×** — the 39.09-vs-0.00 gap **will not reproduce**. On `box_exact` and `csg_difference`,
where the features are sharp, DC's mean exceeds MC's by **> 5×**. So the published constant is a property
of the skull's feature content, not of DC.

**WHY IT IS WORTH TWO LINES OF CODE ANYWAY.** `mean |f(v) − isovalue|` and `max |f(v) − isovalue|` per
algorithm per field is **the only direct measurement of the Tikhonov regulariser's bias**, and nothing
currently measures it. It is also the cheapest possible cross-check on §2.2's convergence prediction: if
DC's mean `|f(v)|` does not fall like `h²`, the convergence result is already decided before the sweep
runs.

---

## 2.6 — A thickness-parameterised reference field (the eighth field, again)

**CLAIM.** Nothing in the current seven fields probes the *representational floor* — the thickness at
which a feature stops being representable — and that floor is where every game-visible failure lives.

**SOURCE.** Splietker & Behnke, *Directional TSDF*, arXiv `1908.05146`, states the failure precisely:
*"it cannot represent anything thinner than the voxel size. While loss of fine details would be
acceptable, the extracted surfaces can become **completely wrong**."* Their fix — six direction-bucketed
fields combined by a weighted vote at the level of MC indices, so *"the intersection of the indices is
equivalent to the binary `and` operation"* — is cheap at the meshing layer and **6× the field memory**,
which is probably disqualifying for a chunked world at full resolution.

**Take the field, not the technique.** `thin_slab(t)` parameterised by thickness, swept over
`t/h ∈ [0.1, 3]`, reporting (component count, surface area, χ, is_closed) for every extractor. This
repo has three measurements that all live on this axis and none that parameterises it: **M-95** (subgrid
resolves a slab 1/20 of a cell thick where a sign test returns nothing), **M-72** (Marching Cubes
*aliases* a sub-cell feature rather than losing it — `thin_plate` 4,088 → 1,016 → 248 → **56** triangles
and still 56 at `h = 0.5`), **M-56/M-57** (greedy quads return **zero** on the same field, cleanly,
because they ask one question per cell centre).

**PRE-REGISTERED PREDICTION.** Surface area against `t/h` shows three distinct signatures: greedy quads
**discontinuously collapse to zero** at `t/h ≈ 1`; Marching Cubes **degrades continuously into a
resolution-dependent scatter** with no clean cutoff; subgrid MT stays within 5% of the true area down to
`t/h ≈ 0.05`. **The middle one is the worst behaviour for a game** and the round-2 doc already names why:
a feature that vanishes at a known distance can be faded; one that disintegrates pops.

**HARNESS.** One new field, one sweep, four counters. Cheapest addition in Part 2, and it makes M-72,
M-95 and M-56 rows on a single curve instead of three anecdotes.

---

## 2.7 — Quality metrics that make this repo's numbers comparable to published tables

**CLAIM.** The validity report currently records quantities that no paper reports, and omits three that
every recent paper reports — so its numbers cannot be placed beside anyone else's.

| Adopt | Definition | Replaces / adds to |
|---|---|---|
| **SA<10°** | % of triangles whose smallest internal angle is below 10° | replaces "degenerate (near-zero-area)" — scale-invariant, and directly comparable to FlexiCubes and TetWeave |
| **IN>5°** | % of nearest-point pairs whose face-normal angle exceeds 5° | new; the metric with the strongest claim to tracking what a player sees, because normals drive shading |
| **NV / NE / SI (%)** | non-manifold vertices / edges / self-intersecting triangles, as percentages | this repo reports counts; percentages are what published tables carry |
| **F-score at τ** | precision and recall of points within τ | separates "made stuff up" from "missed stuff", which Hausdorff and Chamfer conflate |

**Three warnings, each of which would corrupt a number if ignored.**

1. **Chamfer distance is not portable and the papers say so.** FlexiCubes samples 100,000 points;
   TetWeave samples 1,000,000; FlexiCubes states outright that its Table 5 CD numbers "cannot be
   directly compared" against its own Section 5 numbers. Tatarchenko et al., arXiv `1905.03678`,
   demonstrate CD is dominated by *where* the wrong geometry sits. **Do not adopt CD.**
2. **The F-score threshold IS the metric.** FlexiCubes τ = 0.003, TetWeave τ = 0.001, Tatarchenko
   recommends ≤ 1% of side length. Quoting F1 without τ is meaningless.
3. **"Radius ratio" has at least three incompatible definitions in current use, and one published
   inconsistency.** Thingi10K: aspect ratio = circumradius / incircle diameter, lower better. TetWeave:
   AR = longest edge / shortest altitude, lower better; **RR = inradius / circumradius, stated "lower
   values indicat[e] better triangle regularity"** — which is backwards, since an equilateral triangle
   maximises `r/R` at 1/2 and degenerate triangles drive it to 0. Both FlexiCubes and TetWeave present
   RR as a headline number. **Pin the definition or do not use the name.**

**A published result worth knowing before adopting NV(%).** FlexiCubes and TetWeave report MC and MCSeg
at **NE = 0.0, SI = 0.0 but NV = 47–52%**. If standard Marching Cubes output really is ~50%
non-manifold *at vertices* while being edge-manifold and intersection-free, then any harness reporting
0 there is measuring something else — which is direct external support for this repo's insistence on a
link walk rather than an incidence count. **M-53 reports 0 non-manifold vertices for Marching Cubes
across every field.** Those two statements cannot both be about the same quantity, and reconciling them
is a half-day that would sharpen O-12.

**And the tension nobody has resolved.** Etiene's Table 2, for Macet: no fix → *good quality, bad
accuracy*; either fix → *bad quality, good accuracy*. They raise, without resolving, whether this is a
theorem or an implementation accident. **Any project optimising triangle quality in a dual method should
treat that as a live hypothesis and test it, not assume it away.**

---

## 2.8 — Smaller items, each with the reason it is small

| # | Claim | Source | Verdict |
|---|---|---|---|
| a | Gradient-free sharp-feature DC: a per-cell QP from *sampled SDF values alone*, no Hermite data | *Dual Contouring of Signed Distance Data*, arXiv `2604.00157`, SIGGRAPH 2026 (ABSTRACT only, no numbers) | **Read it.** This is the destructible case exactly — after a player carves, the analytic field is gone and only samples remain. Stays a 3×3-ish solve, so ✗16's no-glam argument survives |
| b | DC over 2×2×2 "expanded cubes" extracts manifold **and non-manifold** surfaces from *unsigned* fields | *DCx*, TOG/SIGGRAPH 2026, DOI `10.1145/3811388`, code `github.com/jjjkkyz/DCx` | The cheapest published route to thin sheets and open surfaces, and it stays table-driven. **Transcribe from the code, never reconstruct** |
| c | Guaranteed plane angles: isosurface stuffing bounds min angle at **16.43°** and max at **125°**, competitive with Chew's 120° "at a fraction of the effort" | Labelle & Shewchuk, DOI `10.1145/1275808.1276448`; improved: DOI `10.1145/2504459.2504507` | Needs a BCC lattice — new topology, new table, new LOD story — **and it rounds off sharp corners**, the opposite of the DC path. **Measure first**: add a min-plane-angle histogram; if MC's observed minimum is already 5–10°, a 16.4° guarantee buys little |
| d | Warping the *grid* (not the mesh) removes degenerate triangles at ≤ 2× cost | Dietrich et al., DOI `10.1109/tvcg.2008.60` — **corpus copy is a repository landing page, not the paper**; corroborating rates from Custodio et al. DOI `10.1186/s13173-019-0086-6` (aneurysm **13.48%** of triangles with q < 0.15; bonsai **10.7%**) | Both rate figures are medical volumes. **Measure this repo's own first.** Also: it moves vertices, so every golden hash changes |
| e | Attribute-aware simplification needs **less than half the faces** for equal RMS error; λ = 0.02 for normals, 0.1 for colour | Hoppe, DOI `10.1109/visual.1999.809869` | Predicted **< 20%** win here, because this repo's attributes are piecewise-constant material IDs whose error concentrates entirely on boundaries, where a boundary constraint is cheaper and exact |
| f | Booleans on *arbitrary* oriented meshes via generalized winding number thresholded at ½ | Jacobson, arXiv `1601.07953`; foundation DOI `10.1145/2461912.2461916` | **Wrong layer.** Field-level CSG (which `BrushStack` already does) beats mesh-level booleans by orders of magnitude. **But steal the metric**: sample `w` at points whose inside/outside is known from the analytic field and report the disagreement rate — a solidity check *stronger than manifoldness* that this repo does not have |
| g | Angle-weighted pseudonormal is provably the **only** common analytic pseudonormal that correctly signs inside/outside, including on the medial axis | Bærentzen & Aanæs, DOI `10.1109/tvcg.2005.49` | Area-weighted vertex normals — the common default — are *not* the correct extension of face normals under this analysis. Predicted difference **< 2°** on well-shaped triangles, **> 15°** on triangles with radius ratio < 0.15, i.e. it matters exactly where the slivers are. ~40 lines to measure |
| h | Manifold-guaranteed adaptive simplification via a vertex tree *separate from* the octree; LOD by **re-tagging, not re-meshing** | DOI `10.1109/tvcg.2007.1012` | The re-tag property is the interesting half: changing an LOD threshold should cost **< 5%** of a full re-extract. Measure the premise first — sweep this repo's existing LOD threshold and count non-manifold edges; if flat zero, the second hierarchy is unjustified |
| i | Runtime convex decomposition: **4.6× fewer components** at lower concavity than HACD, and volume-based concavity wrongly fills slots and spouts | CoACD, arXiv `2205.02961` (this repo's V-15) | **~200 s per shape, single-threaded.** M-116 already measured this repo's own version at 241–272 ms per fragment = 14–22 whole frames. Offline only, and predicted **±20%** of V-HACD on near-convex terrain chunks. **Measure the premise**: histogram hull-to-mesh volume ratio per chunk; if most are near-convex, one hull per chunk is the right primitive and the line is moot |
| j | Fracture modes: `k` natural ways of breaking from a sparsity-regularised eigenproblem, impact response by linear projection, *"at no runtime cost"* | *Breaking Good*, DOI `10.1145/3549540` | Needs a **tet volume mesh** and per-shape offline precomputation, so it cannot apply to freshly re-meshed chunks. Value is zero unless a tet-extraction path exists. Prerequisite experiment: can this repo produce a usable tet mesh from a chunk at all — report inverted-tet count and min dihedral |
| k | SDF↔SDF collision without meshing: penetration depth, contact points and normals, real-time, accepting continuous *and* discrete SDFs | DOI `10.1016/j.cagd.2024.102305` (Tencent MoreFun co-authored) | Directly relevant to §1.3's 45%: a game that collides against the *field* skips collider generation entirely for field-defined geometry. Different product, but it reframes the 45% as optional rather than inherent |
| l | Pop-free LOD by geomorphing — and its unsolved half, stated by Lengyel himself | `transvoxel_dissertation_lengyel2010`, §6.1 and ch. 6 | Positions morph cleanly; *"there is no analogous solution for the morph between the primary and secondary normal vectors."* His scalar shortcut is offered with the guarantee withheld: *"(It is unclear whether this point always exists.)"* **That is a cheap, settleable open question — see §3.6** |

---

# Part 3 — Harnesses

Two kinds: instruments to adopt, and instruments that do not exist anywhere. The second list is the one
worth caring about.

---

## 3.1 — Metamorphic relations: the highest value per line in this document

No ground truth required at all. Each relation tests a property of the *transformation*, not of the
output. **The methodology sweep found no published application of metamorphic testing to isosurface
extraction or mesh generation** — the graphics precedent is compiler testing (Lascu's thesis, Donaldson's
MET work), not geometry.

| Relation | Statement | Exactness | What it catches that nothing here catches |
|---|---|---|---|
| Isovalue shift | `extract(f + c, λ + c) == extract(f, λ)` | **bit-identical** | isovalue plumbing, sign-rule asymmetry |
| Positive scaling | `extract(k·f, k·λ) == extract(f, λ)`, `k > 0` | exact in exact arithmetic; **the f32 drift magnitude is itself a conditioning measurement** | numerical conditioning of the vertex solve |
| **Sign flip** | `extract(−f, −λ)` == the same surface with **every triangle reversed** | exact | **strictly stronger than the edge-traversal orientation check** — it exercises the whole case table's sign handling rather than internal consistency. A-019 drove flipped edges to zero; this asks whether the winding is *right*, not merely consistent |
| Integer-cell translation | output translates identically, bit-identical | exact | indexing off-by-ones |
| Fractional-cell translation | output does **not** match, but Hausdorff to the analytic surface is statistically invariant | statistical | **a bimodal distribution indicates a case-table asymmetry** — and this is §2.4's genus-stability sweep from the other side |
| 90° rotation | output rotates identically, bit-identical | exact | case-table symmetry group. Complements ✗12/M-24's lattice equivariance, which covers the *solve* and not the *table* |
| **Chunk decomposition** | `union(extract(chunk_i)) ≡ extract(whole)` on shared interiors | exact **only at power-of-two cell sizes** | the core product requirement. Per-chunk manifoldness cannot see it |
| Resolution doubling | Hausdorff obeys the MMS order | — | this is §2.2 restated |

**The chunk-decomposition row needs M-32's condition attached or it will fail for the wrong reason.**
M-32 measured that two adjacent chunks agree bit-for-bit on **16 of 16** shared-plane vertices at
`h = 0.125` and **0 of 14** at `h = 4/35`, with 22% of 200,000 random `(origin, h, cells, chunk)`
combinations disagreeing by one or two ulp. So the relation is *exact at power-of-two `h`* and
*within-tolerance elsewhere*, and stating it that way turns M-32 from a caveat into an asserted law.

**Cost: ~300 LOC total for all eight, zero dependencies, seconds of runtime.**

---

## 3.2 — Exhaustive case-table sweep and the group-orbit oracle

The Marching Cubes configuration space is **finite and tiny**. Everything about the table can be checked
exhaustively rather than sampled:

- every configuration's boundary edges match the face-boundary pattern of its neighbours — the crack-free
  condition, over `256 × 6` face-adjacent pairs;
- every configuration is consistently oriented;
- the table is **equivariant under the 48-element cube symmetry group** and, for MC33-style tables, under
  complement;
- all `256² × 6` face-adjacency pairs produce compatible boundaries — ~400k checks, milliseconds.

**The orbit oracle.** *"Counting Cases in Marching Cubes: Toward a Generic Algorithm for Producing
Substitopes"* (`cs.upc.edu/~pere/PapersWeb/SGI/MarchingCubes.pdf`, **no DOI returned by any tool**) uses
computational group theory: a coloring group acts on the set of all colorings, **each orbit is one case,
so the number of orbits is the number of cases**. Their own framing — *"a tool for computational algebra
independently confirms results that cannot reasonably be checked by hand"* — is exactly rule 5's problem.

**Partly redundant here, and precisely where it is not.** ✗11 established that this repo's Marching Cubes
table is **derived at compile time by walking each face counter-clockwise**, which is why
`face_disagreements` is structurally zero. So the orbit oracle has nothing to catch there. **Where it has
teeth is everything transcribed**: the MC33 decider, Grosso's inner-hexagon ring (V-31 found the paper's
own listing corrupt), the tunnel triangulation (M-228 found an undefined case reachable only on
configuration 13), and the interior-vertex branch selection taken from a reference implementation that
M-219 found carries a typo. **Every one of those is a transcription, and none has an independent oracle.**

**Cost: ~200–400 LOC, zero dependencies, sub-second, runs as a `#[test]`.**

---

## 3.3 — Per-chunk latency tail, because a mean cannot see a hitch

Criterion reports a mean with a confidence interval. **A voxel editor's failure mode is not "mean re-mesh
is 3 ms", it is "one chunk in a thousand takes 40 ms and drops a frame."** M-124 already measured the
amortised side beautifully — budget swept over 320×, mean landing within 5% under the budget asked for,
and a control showing the unbudgeted queue costing 20.62 ms in a single frame against a budgeted peak of
2.10 ms, a **9.8× lower peak for identical total work**. What is missing is the distribution.

**NVIDIA's percentile spacing is `100 − 2^n`** — 99.0, 99.5, 99.8, 99.9, 99.95 — concentrating resolution
where stutter lives, with the explicit caveat that percentiles near 99.95% on small datasets may
represent a single frame and will not reproduce. Their argument against summary statistics is worth
quoting in the harness's doc comment: *"a metric like '4% of frame times were at least 20 ms above
average' directly correlates to subjective perception, whereas a standard deviation figure remains
abstract."* And **the consumer-side "1% low" has two incompatible definitions in circulation** — the 99th
percentile frame time, and the mean of the worst 1% of frames. If this repo ever reports one, it must say
which.

**Cost: ~200 LOC plus a trace format. The hard part is defining a *representative* edit trace, which is a
design decision rather than a coding one** — and E-202's `ISOMESH_AUTOCARVE=60` is already most of it.

---

## 3.4 — Three statistical rules that would change existing conclusions

**(a) Harmonic mean, not geometric mean, for aggregating speedups.** The Ghent IISWC20 methodology paper
(`users.elis.ugent.be/~leeckhou/papers/iiswc20.pdf`, fetched in full) measured PyPy at **4.4× by geometric
mean and 0.97× by harmonic mean** against CPython — the two aggregations disagree about whether the
system is faster *at all*. The harmonic mean has a physical meaning; the geometric mean assumes
log-normality, which they disproved on their own data via Shapiro-Wilk, D'Agostino and Anderson-Darling.
**Directly actionable here**: any future "algorithm X is 1.3× faster across the seven fields" is
true or false depending on the aggregation, and this repo has already been bitten by field-dependent
constants five times (§0.2). Report the harmonic mean *and* the per-field numbers.

**(b) Randomise order and environment.** Mytkowicz, Diwan, Hauswirth & Sweeney, *"Producing Wrong Data
Without Doing Anything Obviously Wrong!"* (ASPLOS 2009, no DOI returned; PDF at
`users.cs.northwestern.edu/~robby/courses/322-2013-spring/mytkowicz-wrong-data.pdf`) is the canonical
demonstration that **UNIX environment variable size and object-file link order** shift measured
performance enough to reverse a study's conclusion. **M-197 hit exactly this in this repo** — batching
read-backs appeared to make an untouched path 75% slower, four runs consistently, and swapping the order
of the two measurements moved the slowness to the other path. The rule earned there generalises: the
prescribed defence is randomising the biasing factors so bias becomes noise.

**(c) Separate start-up from steady-state, and never average them.** Their concrete protocol: multiple
process invocations for start-up, continuing until the CI is within 5% of the mean or 30 invocations;
steady state declared when the coefficient of variation of the last **k = 4** iterations is **< 2%**;
report 95% CIs with the Student t. **M-145 is what this rule looks like when it is violated** — the
GPU-vs-CPU conclusion was drawn from cold numbers, one of which absorbed 10.76 ms of shader compilation.

**One thing this repo does that the sources do not.** No published precedent for **pre-registration** in
graphics performance work was found in any sweep. P-1 through P-7 appear to be unpublished practice, and
that is itself a finding.

---

## 3.5 — GPU timing: build the guard before the instrument

**gfx-rs/wgpu issue #9414: on macOS 26 with Metal 4, `MTLCounterSampleBuffer`-based timestamp queries
return all zeros.** The buffer creates successfully and reports support, but the driver no longer
populates it — `supportsCounterSampling` now returns false for draw, dispatch and blit boundaries.
Reported against **wgpu 29.0.1** on macOS 26.3.1. Apple's replacement is `MTL4CounterHeap`; wgpu-hal has
not migrated.

**This repo is pinned to wgpu 29 and one of its two machines is an M5 Mac.** M-197 already records that
`ExtractTimings` are CPU wall-clock around each submission and that GPU timestamp queries "need a device
feature this crate does not request" — that decision now has a second, stronger reason. **If timestamps
are ever adopted, the first thing to write is a five-line assertion that the deltas are non-zero**,
because the failure mode is silence, not an error.

**Bottleneck attribution without timestamps** is the pre-timestamp method and it still works: vary the
workload or clock of one stage and see whether total time moves. *Real-Time Rendering 4th* §18.2–18.3 and
*GPU Gems 1* ch. 28 (both in corpus) agree on it, with two warnings that matter — **the bottleneck moves
within a single frame** (*"there is rarely only one bottleneck"*), and **lengthening a shader to test it
can be optimised away by the driver**, so the added work must touch constant registers whose values are
not known at compile time.

---

## 3.6 — Settle Lengyel's own open question, cheaply

Lengyel's dissertation offers a one-scalar geomorph (store a distance along the normal instead of a full
second position, saving 12 bytes per vertex) **with the guarantee explicitly withheld**: *"(It is unclear
whether this point always exists.)"*

**PRE-REGISTERED PREDICTION.** The shortcut fails on `gyroid`: for **> 1%** of high-detail vertices, the
ray `p + t·n` has no intersection with the low-detail mesh in the containing cell whose interpolated
normal has positive dot product with `n`. On `fbm_terrain` the failure rate is **≈ 0%**, because a
heightfield is star-shaped along the up axis and a triply-periodic minimal surface is not.

**HARNESS.** Extract the same region at LOD `n` and `n+1`; for every fine vertex, cast along its normal,
intersect the coarse mesh restricted to the containing cell, and count (a) no intersection, (b)
intersection with negative normal dot product, (c) success. **One field of storage is saved if and only
if the failure rate is zero**, and a non-zero rate turns a twelve-year-old open question into a measured
row. This repo already has both halves: M-121 extracts a block at both its old and new level at the
instant of the switch, and A-011b/c own the transition geometry.

---

## 3.7 — Novel instrument #1: the redundant re-mesh factor

**Nothing like this is published.** The time-varying isosurface literature is entirely about *accelerating*
re-extraction — temporal index trees, out-of-core streaming, compression. Four sweeps found **no paper
that measures whether re-extraction output is locally stable**. The chunked-voxel-game literature is not
academic and publishes no metrics at all.

**The instrument.** After a field edit confined to a ball of radius `r`, the output outside a known
dilation of that ball must be **byte-identical** to the previous output, and inside it must differ only
where the field differs. Measure (a) the symmetric difference of the vertex and triangle sets, (b) the
*spatial support* of that difference, (c) the ratio of `actually-changed triangles` to
`necessarily-changed triangles`, the latter bounded by the edit support dilated by one cell. Call (c) the
**redundant re-mesh factor**; a minimal mesher scores 1.0.

**This repo is unusually close to it already.** E1 measures the *input* side — M-33 (a brush changes
15–36% of cells in its bounding box), M-34 (counting value changes overstates by 2.8–3.7×), M-50 (both
reproduce live under a mouse). What is missing is the *output* side, and it is the half that decides
whether the incremental path is correct: **any change outside the dilated support is a bug**, and a
factor near the chunk/edit volume ratio means the incremental path buys nothing.

**Why it is the single most decision-relevant instrument in this document.** It answers two questions at
once — "is the incremental path correct?" and "is it worth having?" — and the second one is currently
answered by an inference from M-33 rather than by a measurement.

**Cost: ~150 LOC, zero dependencies, reuses the golden-hash machinery for the byte-identity half.**

---

## 3.8 — Novel instrument #2: a field-space popping metric in projected pixels

**The only published popping metrics are image-space.** StopThePop (arXiv `2402.00525`) warps frame `F_i`
to `F_{i+t}` with RAFT optical flow and computes FLIP (DOI `10.1145/3406183`) between warped and true
frames, at `t = 1` and `t = 7`, arguing `t = 7` is more reliable because error accumulates. They show
explicitly that **MSE does not detect popping where FLIP does**. The alternative is arXiv `2208.12674`, a
ResNet18 trained to classify LOD transitions from image pairs — also image-space, and requiring a trained
model. Both need a renderer, an optical-flow network and a perceptual model, i.e. a Python stack that
cannot run in a Rust test.

**A geometric equivalent does not exist.** The proposal: given `M_n` and `M_{n+1}`, compute the
**per-vertex displacement distribution** — not just Hausdorff, because a single spike and a uniform shift
feel completely different — then project through a camera model: at distance `d`, vertical FOV `θ`,
vertical resolution `H`, a world displacement `δ` subtends `≈ δ·H / (2·d·tan(θ/2))` pixels. Report the
**pixels-of-pop distribution and its p99**, optionally weighted by whether the displacement is along or
across the view direction, since across-view motion is far more visible.

**This repo has already taken the hard first step.** **M-121 measured the pop** — worst vertex-to-nearest-
vertex distance **3.136 cells**, typically 0.6–1.6, measured at the instant of the switch by extracting
the block twice on that frame — and noted that *"its size is what decides whether it can be hidden by a
fade, a morph, or nothing at all, and no figure for it exists in the literature review."* The camera
projection is the missing half, and it converts a number into a **policy**: *switch LOD at the distance
where p99 pixels-of-pop < 1*. No published instrument produces that rule.

**The honest caveat, which must be in the doc comment.** Without validation against `FLIP_t` on rendered
sequences or a StopThePop-style user study (18 participants, forced choice, Wilcoxon signed-rank), this
is a **plausible** metric, not a perceptual one, and should be labelled as such.

**Cost: ~200 LOC for the geometric half, zero dependencies.**

---

## 3.9 — Novel instrument #3: a mutation-coverage map of geometry instruments

**Every verification paper in this area demonstrates that *its* instrument found *some* bugs. None
publishes the coverage matrix.** Etiene 2012 comes closest — eight codes against three instruments — but
that measures the codes, not the instruments.

**The instrument.** Hand-curate a mutant set (flip a winding; swap two case-table indices; change `<` to
`<=` at an exact zero; perturb one interpolation coefficient; drop an axis from an exclusion list) and
record **which harness catches which mutant class**: Euler characteristic, orientation traversal, vertex
link walk, MMS position order, MMS normal order, Hausdorff, self-intersection census, golden hash,
metamorphic sign-flip.

**This repo has been building that matrix by accident for months and has never written it down.** M-119
(the obvious repair would have gutted the assertion it repaired, proven by mutation), M-200 (inverting
the `1 − t` branch passes all 32 tests in its module), M-203 (the detector moved and had to be
re-verified rather than assumed), M-83 (the self-intersection counter is structurally blind inside a
Steiner fan — 190 of 190 pairs skipped), M-204 (a passing test could not catch the defect it asserted
against), M-195 (an instrument reading 0 on a mesh with 56 holes). **Six data points for a matrix nobody
has published.**

A matrix saying *"orientation traversal catches 100% of winding mutants and 0% of case-table-index
mutants; the link walk is the only thing that catches bowtie mutants; the golden hash catches everything
and tells you nothing"* is directly useful to every future implementer, and it is exactly the kind of
finding this ledger is built to hold.

**Cost:** `cargo-mutants` as a dev tool plus a curated set — automated mutants are dominated by
uninteresting ones. ~1 day of compute for a full run. **The value is the matrix, produced once.**

---

## 3.10 — Novel contribution #4: the benchmark corpus that does not exist

§0.3 established the absence. What a fair voxel/SDF corpus would need, synthesising Thingi10K's design
argument (arXiv `1605.04797`) with Berger's three-phase pipeline (DOI `10.1145/2451236.2451246`, code at
`github.com/fwilliams/surface-reconstruction-benchmark`) and Etiene's negative results:

1. **Analytic fields with derivable formal orders per property**, spanning: exactly-representable (linear
   — the null control), smooth low-curvature (sphere), high-curvature (thin-tube torus),
   transcendental/high-frequency (gyroid — the field class that broke Afront's spline shortcut),
   boundary-crossing (fbm heightfield), and thin features approaching one voxel (§2.6's `thin_slab`).
2. **Ambiguous-cell density as a controlled axis**, not an accident. The MC33 case-13 and SnapMC failures
   are both ambiguity failures, and M-40 already measured this repo's own density: 0 on five of seven
   fields, 0.515% on `gyroid`, 1.532% on `fbm_terrain`.
3. **Near-degeneracy as a controlled axis** — corner values approaching exactly `λ`. Where SnapMC failed
   > 50%, where MC's slivers come from, and where M-231 found that quantisation *creates* the singular
   configuration.
4. Real CT volumes for performance and topology-invariant claims only, **explicitly marked as having no
   accuracy ground truth**.
5. Mesh-derived SDFs from Thingi10K/ABC subsets, with conversion error **bounded and reported separately**
   from meshing error.
6. Per-field published metadata so the harness never branches on field identity.

**Why this repo is unusually well placed, and it is not a matter of effort.** The three hardest parts of
making such a corpus reproducible are already done: the `ReferenceField` metadata design; the rule that
test code never branches on field identity; and **M-31 — 63 golden hashes generated on macOS/arm64 pass
unchanged on Linux/x86-64, bit-for-bit, because the crate is `libm`-only.** Every existing corpus depends
on platform libm and therefore cannot ship reproducible expected values. That property is the thing a
distributable benchmark needs and the thing nobody else has.

**Cost: this is a paper, not a ticket.**

---

# Part 4 — Corrections, contradictions and resolved open questions

## 4.1 — One external claim is contradicted by a measurement in this repo

The external sweep concluded that *"the crate pins wgpu 29.0.3 to Bevy 0.19, which predates even the
experimental flag"* and predicted a capability probe would report `EXPERIMENTAL_MESH_SHADER` **absent in
wgpu 29 entirely**.

**That is false, and M-146 measured it false.** GPU-007's probe on wgpu 29.0.4 reported
`EXPERIMENTAL_MESH_SHADER` **advertised** on the RTX 3090/Vulkan adapter, and M-147 then measured
`mesh_shader=true` on Bevy 0.19's own `RenderDevice` with no configuration and no `unsafe` in this
repository. The reasoning behind the external claim — reading a feature table rather than probing —
is the same error M-146 itself made in the other direction and M-147 corrected.

**What the external sweep does add, correctly, is the wgpu v30 status**: the mesh-shading spec lists
SPIR-V and HLSL as in progress and **MSL as "Planned"**, with `Features::EXPERIMENTAL_MESH_SHADER`,
`_MULTIVIEW` and `_POINTS`, plus a documented "queries are unsupported" limitation. That is consistent
with **V-23**, which resolved the same contradiction from the source: the *feature* reaches Metal, the
*WGSL compiler* does not, so on Metal a caller supplies pre-compiled MSL. **O-5 can be closed**: mesh
shaders exist in Metal-the-API since Metal 3 with hardware support on M3/A17 and later; wgpu's Metal
backend does not expose the WGSL path; both halves of the original contradiction were true.

**The general lesson is one this repo keeps earning.** A capability table is a claim about a backend; a
probe is a measurement of a device. M-147's amendment already says the runtime probe is load-bearing
rather than belt-and-braces, and this is a second instance of an outside reader reaching the wrong
conclusion from the table alone.

## 4.2 — Three candidates demoted by measurements already in the ledger

| Candidate | Why the literature ranks it high | Why it ranks lower here |
|---|---|---|
| Decoupled-lookback / Decoupled Fallback scan (§1.4) | 1.4× on the scan, and a documented crash class | GPU-010a already moved the scan on-GPU (5.24 → 0.37 ms) and M-155 took the whole path to 0.54 ms. Kept as a **trap warning**, not a ticket |
| Workgroup-shared reduction staging (§1.6) | **14.5×**, the largest single multiplier in the corpus | Same reason. Kept as a checklist for the next GPU reduction written here |
| Stage-streamed / SIMD restructuring (§1.5) | ~4× overall, measured, read in full | §0.1 caps it at **1.41×** on the whole job. The dedup and CPU/GPU-unification side effects are worth more than the speed |

## 4.3 — Corpus integrity notes that would corrupt a number if ignored

Collected from all four sweeps. Each is a place where reading the corpus markdown gives a wrong figure.

- **Landing pages masquerading as papers.** `10.1109_tvcg.2008.60` (Dietrich, edge transformations) is a
  UFRGS repository page with Portuguese navigation chrome — the abstract is quotable, the tables and
  error bounds are not in the corpus. `appearance-preserving-simplification-4ispnnapzy` is a
  citation-graph stub. **Dupuy's CBT paper `10.1145/3406186` is a HAL landing page** — the 5.78 → 0.40
  figure comes from the Unity slides, not from it.
- **Column-scrambled tables.** FlexiCubes' Table 1 property taxonomy, CoACD's Table 5, Aokana's Table 1,
  Dyken's Table 1 (unusable outright), Schmitz's Figure 5 (image-only), Thingi10K's Figure 7 (renders
  every quality percentage as 100%, plainly wrong — the body text says 45% of models self-intersect).
- **`abstract_search` is broken for this domain.** A query for "adaptive octree dual contouring
  error-guided simplification" returned three unrelated Nature/NeurIPS papers with `null` titles, while
  `distill_search` on the same query returned Ju et al. 2002 at 0.743. Use `distill_search` plus a
  catalog grep.
- **Year fields in `distill_search` metadata are unreliable** — `labelle_shewchuk_isosurface_stuffing_2007`
  → 1976; `sig2024_Multi-Material…` → 2006. Take years from body text.
- **Many high-value stems have empty `title` and `doi`** — `labelle_shewchuk…`, `dualsimp_tvcg`,
  `10.1145_566570.566586`, `SchaeferWarren2`, every `sig20xx_*`. `catalog_backfill_title` would pay for
  itself.
- **Substantial duplication** — `10.48550_arXiv.2205.02961`/lowercase twin, `dualsimp_tvcg` /
  `10.1109_tvcg.2007.1012`, every `GameAIPro2_Chapter*` / `gameaipro2-ch*` pair.

## 4.4 — Corpus absences, established against the full stem index

Per ✗4's rule, these were checked against a complete `catalog_list` dump (9,354 stems), not by a failed
search.

| Missing | Why it matters |
|---|---|
| **Flying Edges** (Schroeder, Maynard, Geveci, LDAV 2015) | The SOTA shared-memory-parallel Marching Cubes, what VTK ships, and the conclusion of §1.5's restructuring. Highest-value acquisition in this document |
| **Ju & Udeshi 2006, "Intersection-free contouring on an octree grid"** | Cited 3× in the corpus, present 0×. The load-bearing citation for §2.1, a guarantee this crate does not hold |
| **Any GPU work-graphs paper** | The `s2024-advances` decks are the closest and cover cluster/visibility-buffer rendering instead |
| **Any mesh-shader capability or performance paper** | The only mesh-shader content is a forward-looking CBT slide with no measurement. This repo's M-146/M-147/V-23 are, locally, better evidence than anything in the corpus |
| **Incremental / temporally-coherent re-meshing of an edited SDF** | The actual workload. Nearest neighbours are SLAM reconstruction (`10.1145/2508363.2508374` voxel hashing, `arXiv 2311.00626` nvblox, `10.1109/icra.2018.8463157`), none read. **This is the gap behind §3.7 being a novel instrument** |
| **Temporal coherence metrics for isosurfaces** | Searching returns video retargeting and flow-field topology. **This is the gap behind §3.8** |
| **Triplanar / UV alternatives for procedural surfaces** | All texturing hits are mesh parameterisation or virtual texturing |

## 4.5 — Open questions this document can move

| # | Status |
|---|---|
| **O-5** (mesh shaders on Metal) | **Closeable.** V-23 resolved it from wgpu's source; §4.1 adds the wgpu v30 spec status. Feature reaches Metal, WGSL compiler does not |
| **O-7** (what fraction of *our* pipeline is contouring) | **Answered by M-135/M-136 and sharper than V-4's 54%**: 29.0% mean, 13.1–74.3% range, with the collider check at 45%. §0.1 is that answer written as a decision |
| **O-8** (does DC need f64) | §2.2's convergence prediction is a second, independent probe: if fixed-λ Tikhonov is not second-order, the question is about the *regulariser*, not the precision |
| **O-10** (Surface Nets' non-manifold rate vs thickness/`h`) | §2.6's `thin_slab(t)` sweep gives the parametrised form the row explicitly says is still open |
| **O-11** (why the dual topology goes superlinear) | §1.8's Morton experiment is the cheapest remaining probe. The 128³ per-sample spike on two unrelated cache hierarchies is still unfollowed |
| **O-12** (is Marching Cubes unconditionally manifold) | §2.7's NV(%) reconciliation is a new angle: two published tables report MC at 47–52% non-manifold *vertices*; M-53 reports 0. Both cannot be measuring the same quantity |

---

# Part 5 — Acquisition list

Ranked. **Only identifiers returned by a tool.**

| # | Work | Identifier | Why |
|---|---|---|---|
| 1 | **Flying Edges** — Schroeder, Maynard, Geveci, LDAV 2015 | none returned; search by title | Absent from the corpus; the SOTA parallel MC and the conclusion of §1.5 |
| 2 | **Ju & Udeshi, Intersection-free contouring on an octree grid**, PG 2006 | none returned; PDF at `cs.wustl.edu/~taoju/research/interfree_paper_final.pdf` | §2.1. Cited 3× in the corpus, present 0× |
| 3 | **Subgrid Marching Tetrahedra** — Baktash, Gillespie, Crane | DOI `10.1145/3811358`, arXiv `2606.00454` | Already implemented here; **acquire for the missing timing table** that no public page carries |
| 4 | **Accelerating Signed Distance Functions** — Hubert-Brierre et al., CGF 2025 | DOI `10.1111/cgf.70258` | §1.2. If field evaluation is >70% of runtime this is where the budget belongs |
| 5 | **Dual Contouring of Signed Distance Data** — Carrera, Wang, Batty, Stein, Sellán | arXiv `2604.00157` | §2.8a. Gradient-free sharp features — the post-carve case |
| 6 | **Topology Verification for Isosurface Extraction** — Etiene et al. | DOI `10.1109/tvcg.2011.109` | §2.3's Eq. 4.7. Must be read, never reconstructed |
| 7 | **Occupancy-Based Dual Contouring** appendix §A.1.5 — Hwang & Sung | arXiv `2409.13418` (**already in corpus, indexed**) | The multi-vertex extension of §2.1. Read the appendix rather than downloading anything |
| 8 | **Decoupled Fallback** — Smith, Levien, Owens, SPAA '25 | DOI `10.1145/3694906.3743326` | Acquire the **code** (`github.com/b0nes164/GPUPrefixSums`), not the paper |
| 9 | **DCx** — Bao et al., TOG 2026 | DOI `10.1145/3811388`, code `github.com/jjjkkyz/DCx` | §2.8b. Thin sheets and open surfaces, table-driven. **Transcribe from code** |
| 10 | **Sphere Carving** — Schott et al., TOG 2025 | DOI `10.1145/3730845` | §1.1's hierarchical form, and it claims polygonisation specifically |
| 11 | **Real-Time GPU Tree Generation** — Kuth et al., HPG 2025 | none found; `coburggraphicslab.github.io/publication/Kuth25RTG.html` | Acquire for the **baseline comparison the SIGGRAPH 2025 course slides omit** — nobody has published how much of the win is work graphs versus merely being GPU-resident |
| 12 | **The PhaseTree** — Galin et al., TOG 2026 | DOI `10.1145/3811379` | Multi-material at <25% overhead. **An API design question to raise, not decide** |

**Explicitly rejected, with reasons, so they are not rediscovered.** Power-diagram adaptive extraction
(DOI `10.1111/cgf.70037`, arXiv `2506.09579`) — sequential incremental Delaunay insertion fights both
chunked editing and the GPU path, and its win is in *evaluation count*, which is nearly free for
procedural fields. TetWeave (DOI `10.1145/3730851`) — rebuilds a Delaunay triangulation per optimisation
step; structurally offline; read once for the guarantee vocabulary. Fracture modes (DOI `10.1145/3549540`)
and CoACD (arXiv `2205.02961`) — offline per-shape costs, incompatible with per-edit re-meshing, and
M-116 already measured the local version at 14–22 frames per fragment. Directional TSDF (arXiv
`1908.05146`) — 6× field memory; take the `thin_slab` field, not the technique. Mesh-level booleans
(arXiv `1601.07953`) — wrong layer; take the winding-number *validity probe* instead.

**And a standing check worth one tool call a year.** No 2024–2026 Gaussian-splatting paper found by three
keyword sweeps touches editing, destruction, collision or CSG — the SIGGRAPH 2026 3DGS cohort is entirely
rendering-side (quality, scale, convergence, displays). The prediction: searching SIGGRAPH 2027 and I3D
2027 titles for "destruction", "fracture", "carve", "CSG" or "collision" alongside 3DGS returns zero.
Splatting owns captured static appearance; meshing still owns authored geometry that must be collided
against and carved.

---

# Part 6 — What I would do first, and why in this order

| | Action | Cost | What it decides |
|---|---|---|---|
| 1 | **The stub-field ratio bench** (§1.2a) — run `cargo bench --bench extract` twice, once against an `#[inline(never)]` constant field | one afternoon, zero deps | Whether *any* further extractor work is worth doing. If field evaluation is >70%, the entire Part 1 ranking inverts and §1.1/§1.2 are the only speed tickets that matter |
| 2 | **The metamorphic relation suite** (§3.1) | ~300 LOC | Whether the case table's sign handling, symmetry and chunk seams are correct — with no ground truth at all. Highest value per line in this document, and the sign-flip relation is strictly stronger than the orientation check that exists |
| 3 | **The self-intersection baseline table, then quad splitting** (§2.1) | table is nearly free; the fix is M | Whether a *guarantee* replaces a recorded metric. ✗2 and M-61 already frame the problem; the mechanism is published and cited three times in the corpus |
| 4 | **The redundant re-mesh factor** (§3.7) | ~150 LOC | Whether the incremental path is correct *and* whether it is worth having. Currently the second half is an inference from M-33 rather than a measurement. Nobody has published this instrument |
| 5 | **The MMS convergence gate with negative controls** (§2.2) | ~200 LOC + the research | Whether fixed λ = 0.01 is allowed to be fixed. ✗12 derived the constant and nothing has asked whether it can survive refinement |
| 6 | **The subgrid grid-edge root cache** (§1.9b) | M | The largest single measured cost in the crate (M-98, 70×/196×). The precondition that deferred it is now met by M-168's identity key, and the acceptance test is zero golden-hash changes |
| 7 | **Pixels-of-pop** (§3.8) | ~200 LOC | Turns M-121's 3.136 cells into an LOD-distance *policy*. Genuinely novel, and the geometric half is cheap |

**The honest ordering argument.** Items 1–4 are all bounded, all falsifiable, and none needs a paper this
repo does not have. Item 1 comes first because it can invalidate the ranking of everything below it, which
is the most useful thing a cheap experiment can do. Item 5 is the one most likely to produce an
uncomfortable result, which is why it is worth doing. Items 6 and 7 are the two places where this
project's existing measurements have already done the hard part and nobody has taken the last step.

**What this document deliberately does not do.** It does not propose a new extractor. Four sweeps across
a 9,354-document corpus and the 2023–2026 literature found no algorithm that would beat what is already
implemented here on the axis that matters, and this repo's own M-135/M-136/M-155 say the extractor is not
where the time is. The frontier for this project is **field evaluation, data movement, and instruments** —
and two of the three novel instruments in Part 3 are eighty percent built already.
