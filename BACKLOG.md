# isomesh — BACKLOG

**Updated:** 2026-08-11
**Companions:** `CLAUDE.md` (rules), `docs/2026-08-11-implementation-brief.md` (the how),
`docs/2026-08-11-bevy-examples-catalog.md` (example detail), `docs/research/` (the why).

---

## How to work this backlog

1. Take the **topmost unblocked, unchecked ticket**. Don't cherry-pick interesting ones — the order
   encodes dependencies and the test harness exists so later work is cheap.
2. One ticket = one commit (or a short stack). Commit message starts with the ticket ID.
3. **Check the box in this file as part of that same commit.** This file is the state.
4. If a ticket can't be finished, leave it unchecked, add a `> BLOCKED:` line under it saying exactly
   what's in the way, and move to the next unblocked ticket. Do not half-finish and check the box.
5. If a ticket turns out to be wrong or to need splitting, edit it and say so in the commit.

### Definition of done — applies to every ticket

- Code compiles with no warnings. `cargo clippy -- -D warnings` clean.
- `cargo test -p isomesh` green. `cargo fmt` applied.
- **`grep -ri bevy crates/` returns nothing.** Non-negotiable — see `CLAUDE.md` rule 2.
- Any algorithm ticket also passes the T-001 validity suite. No exceptions, no "I'll add tests later."
- Any perf claim has a committed benchmark that produced it.
- Public items have doc comments. Anything with a sign convention, a coordinate order, or a winding
  order says so **in the doc comment**, not in a code comment.

**Size key:** `S` ≈ one sitting · `M` ≈ a day · `L` ≈ multi-day, consider splitting.

---

## Phase 0 — Foundation and the test harness

Everything downstream is cheap if this is right and expensive if it isn't. **Do not start Phase 1
before T-001 passes** — an algorithm without a validity harness is an algorithm you can't trust.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **T-003** | **Accuracy harness.** One-sided and symmetric Hausdorff distance from the mesh to an analytic SDF, sampled on triangle centroids and vertices. Also mean absolute error. **Acceptance:** a unit sphere meshed at 64³ has max error below one cell diagonal. | M | I-003 |
| ☐ | **T-005b** | **Property tests over extraction.** Wire T-005a's generators × grid sizes × algorithms through `assert_extracted_mesh_is_valid`, and verify the suite fails when a case-table entry is corrupted. **Acceptance:** the mutation check — a property test that can't fail is decoration. | S | T-005a, A-001 |
| ☐ | **T-006** | **Benchmark harness.** `criterion`, plus a resolution sweep that fits `t = a + b·n³` and **reports `a` separately**. CSV out to `docs/measurements/`. Rationale in the speed analysis: 73% of a published 64³ figure was fixed launch overhead. **Acceptance:** sweep runs 16³→256³ and prints the fitted fixed cost. | M | I-003 |
| ☐ | **T-007** | **Golden-hash regression.** Stable hash of (positions, normals, indices) for each (algorithm, field, resolution). Committed as a JSON fixture. **Acceptance:** a deliberate one-bit change to a case table fails the test with a useful message naming which combination drifted. | S | T-004 |

---

## Phase 1 — The usual suspects

Each algorithm ticket is done when: it passes T-001 with zero violations on all seven test fields, at
three resolutions; T-004 determinism passes; T-005 covers it; and a benchmark exists.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **A-002** | **MC33 / asymptotic decider.** Face saddle `S = (f00·f11 − f10·f01) / (f00 + f11 − f10 − f01)`, guard the denominator. Interior ambiguity via the trilinear body saddle. **Acceptance:** a field that holes under A-001 is closed under A-002, and the test asserts *both* — one that passes with plain MC isn't testing this. | L | A-001 |
| ☐ | **A-003** | **Marching Tetrahedra.** 6-tet decomposition of the cube (document which one — the choice affects the output). Unambiguous by construction, so this is the topological reference the others get compared against. Expect a large triangle count; record it. | M | A-001 |
| ☐ | **A-005** | **Greedy quads / blocky.** Binary occupancy, per-face masks, greedy merge. The budget end of the tradeoff table and the comparison baseline for triangle counts. | M | T-001 |
| ☐ | **A-006** | **Hermite data extraction.** Per sign-changing edge: crossing position (bisection or analytic) **and surface normal there**. `HermiteCell` holding up to 12 crossings. This is DC's real input and the reason DC can do sharp features at all. **Acceptance:** `hermite_debug` gizmo data is correct on `box_exact` — normals on a flat face are identical, normals across an edge are not. | M | I-003 |
| ☐ | **A-007** | **Dual Contouring — closed-form vertex.** The three-plane rotation-equivariant rule from the audit doc. Exactly equivariant, no iterative solve, no condition-number squaring. Falls through to A-008 when the triple product is near zero. | L | A-006 |
| ☐ | **A-008** | **Dual Contouring — regularized general solve.** `M = Σnᵢnᵢᵀ`, `g = Σdᵢnᵢ`, `x = c + adj(M+λI)·g / det(M+λI)`, λ≈0.01. For >3 planes and degenerate cells. **Document that `M = AᵀA` squares the condition number** — that's why QR/Givens exists in the literature, and why f64 matters here. | L | A-007 |
| ☐ | **A-009** | **The cell clamp + measurement.** Clamp the solved vertex to (1−ε) inside its own cell, ε≈1e-4, behind a flag. **Measure self-intersections per 1,000 triangles with and without, on all seven fields, and put both numbers in the commit message.** Cheapest high-value experiment in the project — it decides whether guaranteed intersection-free extraction is free. | M | A-008, T-002 |
| ☐ | **A-010** | **Manifold DC.** Vertex splitting so each cell can emit >1 vertex where topology demands it. **Acceptance:** `non_manifold_edges == 0` on `gyroid` and `csg_difference`, where plain DC will not manage it. | L | A-008 |
| ☐ | **A-011** | **Transvoxel transition cells.** Half-resolution transition cells at LOD boundaries. **Acceptance:** two adjacent chunks at differing LOD produce zero boundary gaps — assert on the geometry, then confirm visually in E-107. | L | A-001 |
| ☐ | **A-012** | **Normal estimation strategies.** Analytic gradient / central differences / area-weighted face normals, selectable. **Acceptance:** all three produce unit-length normals; analytic and central-difference agree within tolerance on `sphere`. | S | A-001 |
| ☐ | **A-013** | **Vertex welding and dedup.** Spatial-hash weld with configurable epsilon, index remap. Deterministic ordering — this is the classic determinism leak, so T-004 must cover it explicitly. | M | A-001, T-004 |
| ☐ | **A-014** | **Subgrid Marching Tetrahedra.** Integer edge-intersection counts instead of vertex signs — resolves features finer than one cell. The differentiator; do it after the usual suspects are solid. | L | A-003 |

---

## Phase 2 — Game-shaped infrastructure

Still zero Bevy. This is the machinery a game needs, living in the core crate where CAD can use it too.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **G-001** | **Chunk abstraction.** `Chunk<T>` with a **1-cell positive-face overlap** so neighbours agree on shared cells. `ChunkId`, world↔chunk↔local coordinate conversion, all three round-trip tested. **Acceptance:** two adjacent chunks meshed independently produce coincident vertices on the shared plane — assert on coordinates, don't eyeball it. | M | A-004 |
| ☐ | **G-002** | **Dirty-set re-meshing.** Track which cells changed; re-mesh only affected chunks. `mesh_dirty(&mut self, out: &mut _)`. Instrument **fraction of cells actually changed per edit** and expose it — this is E1 from the research and it's the ceiling on every incremental idea in the opportunities doc. Unpublished number; log it. | M | G-001 |
| ☐ | **G-003** | **Brush operations.** Add/subtract sphere, box, capsule against the field. `min`/`max` for union/difference, smooth-min variant. Fixed-point (`i32`) storage option for bit-exact determinism. **Acceptance:** a determinism test — 8 brush ops applied in all 40,320 orderings; count distinct results. Expect 1; if it's not 1, you've lost commutativity and the multiplayer story dies here. | M | G-001 |
| ☐ | **G-004** | **Field-derived LOD.** Mip the field (not the mesh), mesh at level N. **Acceptance:** LOD 0..3 all mesh cleanly; LOD *k* has roughly 1/8^k the cells. Pairs with A-011 for the seams. | M | A-011, G-001 |
| ☐ | **G-005** | **Collider export.** `MeshBuffer` → `parry3d::TriMesh` (`Vec<[u32;3]>` indices — parry takes plain arrays). Behind an optional `parry` feature. Optional convex decomposition path. **Acceptance:** a carved shape builds a `TriMesh` without error and passes parry's own validity check. | M | A-013 |
| ☐ | **G-006** | **Frame-budget scheduler.** `mesh_within_budget(ms)` — process the dirty queue until a time budget is exhausted, resume next call. Priority by camera distance. This is the constraint a real game actually operates under and the reason "how fast is the algorithm" is the wrong question. | M | G-002 |
| ☐ | **G-007** | **Chunk streaming.** Load/unload by camera distance with hysteresis so chunks at the boundary don't thrash. | M | G-004, G-006 |

---

## Phase 3 — `bevy_isomesh`

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **B-003** | **Plugin and component API.** `IsomeshPlugin`, a `VoxelVolume` component, a `NeedsRemesh` marker, systems that consume the frame budget from G-006 and drive `AsyncComputeTaskPool` so meshing is off the main thread. **Acceptance:** meshing a large volume does not stall the render loop — show it in the frame-time graph. | L | B-002, G-006 |

---

## Phase 4 — Examples

Two groups. The algorithm demos are quick and prove correctness visually. **The game-shaped ones are
the point** — they're how someone decides whether this crate is usable.

### 4a — Algorithm demos

| | ID | Example | Blocked by |
|---|---|---|---|
| ☐ | **E-102** | `mc33_ambiguity` — split screen, holes vs closed, Euler χ in the HUD for both | A-002 |
| ☐ | **E-104** | `dual_contouring_cube` — **the money shot.** `box_exact` + `csg_difference`: SN rounds the corners, DC holds them. README image. | A-007 |
| ☐ | **E-105** | `marching_tetrahedra` — same field, much higher triangle count | A-003 |
| ☐ | **E-106** | `greedy_quads` — the blocky path, quads before/after merge | A-005 |
| ☐ | **E-107** | `transvoxel_seams` — two LODs adjacent, toggle transition cells. **Commit both screenshots.** | A-011 |
| ☐ | **E-108** | `subgrid_features` — letters carved thinner than a voxel; toggle and watch them vanish | A-014 |
| ☐ | **E-109** | `sharp_features` — live slider on the normal-deviation threshold, through to over-sharpening | A-008 |
| ☐ | **E-110** | `qef_clamp` — clamp toggle, live self-intersections/1k, offending triangles in red | A-009 |
| ☐ | **E-111** | `manifold_check` — non-manifold edges as thick red lines, across all fields × algorithms | T-001 |
| ☐ | **E-112** | `precision_f32_vs_f64` — same field at ~1e6 offsets; f32 cracks, f64 doesn't. Condition number in the HUD. | A-008 |
| ☐ | **E-113** | `normal_estimation` — three panels, lit; differences live in the speculars, not the wireframe | A-012 |
| ☐ | **E-114** | `hermite_debug` — crossings, normals, solved vertex, cell box as gizmos. The view you debug A-007 in. | A-006 |

### 4b — Game-shaped

These use the algorithms the way a game does: chunked, edited, budgeted, collided against.

| | ID | Example | What it has to prove | Blocked by |
|---|---|---|---|---|
| ☐ | **E-201** | `game_terrain_stream` — walk a large fBm world, chunks stream by distance | Sustained 60 fps while streaming. HUD: chunks resident, meshing ms/frame, MB. | G-007, B-003 |
| ☐ | **E-202** | `game_dig` — first-person, click to carve tunnels | The core Minecraft/Deep Rock loop. Re-mesh is imperceptible. Chunks-touched count on screen. | G-002, G-003 |
| ☐ | **E-203** | `game_walk` — character controller on meshed terrain, parry3d colliders | **The acid test.** Walk every chunk seam. No falling through, no invisible walls. If this fails, G-001's overlap is wrong. | G-005, E-201 |
| ☐ | **E-204** | `game_destruction` — shoot a wall, it craters, debris becomes rigid bodies | Runtime fragments are correct physics bodies, not pre-fractured props. Carve a spiral and a hollow shell — that's where decomposition fails. | G-005, G-003 |
| ☐ | **E-205** | `game_lod_flyover` — fly out and back across LOD transitions | No popping, no cracks, no hitching. Transvoxel doing its job at speed. | G-004, A-011 |
| ☐ | **E-206** | `game_budget` — a deliberately overloaded edit queue under a frame budget | Frame time stays flat while the backlog drains. **Amortized cost per frame is the number no paper measures and the only one a game cares about.** | G-006 |
| ☐ | **E-207** | `game_editor` — sculpt with brushes, undo/redo over an op log | The CAD/editor use case. Undo is a re-fold of the log, not a snapshot. | G-003 |
| ☐ | **E-208** | `game_paint` — spray colour on a wall, then blow a hole through it | Paint on the remaining wall is exactly where you sprayed it. Row 4 of the opportunities table. | G-003, B-002 |
| ☐ | **E-209** | `game_csg_props` — place and boolean CSG primitives into the world live | Re-mesh per frame under moving primitives; concave sharp edges hold up. | G-003, A-008 |

---

## Phase 5 — Measurement

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **M-001** | **`bench_shootout`** — MC / MC33 / SN / DC / MT, identical fields and grids, one process, one run. Table: ms, verts, tris, non-manifold edges, self-int/1k, Hausdorff. CSV to `docs/measurements/`. **This comparison does not exist in the literature for post-2020 hardware.** | M | A-001..A-004, T-006 |
| ☐ | **M-002** | **`bench_resolution_sweep`** — 16³→256³, live plot, fits `t = a + b·n³`, **prints `a`**. Expect a large fixed cost at small grids and stop trusting single-grid numbers thereafter. | S | T-006 |
| ☐ | **M-003** | **`bench_stage_breakdown`** — stacked bar: contour / normals / weld / collider / upload. Published comparison: contouring 68 ms vs halfedge construction 58 ms — **the contour was 54% of a usable mesh.** Find your ratio before optimizing anything. | M | G-005, M-001 |
| ☐ | **M-004** | **Write up M-001..M-003** as `docs/research/YYYY-MM-DD-measured-comparison.md`. Numbers, method, hardware, and what surprised you. This is publishable on its own. | S | M-003 |

---

## Phase 6 — GPU (do not start before Phase 5)

The speed analysis is explicit that stage placement dominates the extraction algorithm by roughly an
order of magnitude. Which means GPU work is worth doing — and worth doing *after* you know your own
numbers, or you won't be able to tell what the port bought you.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **GPU-001** | `isomesh-gpu` skeleton, `wgpu 29.0.3`. Public API takes `&wgpu::Device` / `&Queue` / `&mut CommandEncoder`. **Never a Bevy type.** | M | M-001 |
| ☐ | **GPU-002** | Shader composition: `include_str!` + ~40-line `#include`/`#ifdef` preprocessor. **Not `naga_oil`** — see `CLAUDE.md`. | S | GPU-001 |
| ☐ | **GPU-003** | naga CI validation of every shader permutation. No GPU required. ~30 lines, highest value-per-line in the repo. | S | GPU-002 |
| ☐ | **GPU-004** | Compute-shader Marching Cubes + readback. **Headless harness first, no Bevy in the room** — if it can't run against raw wgpu, the abstraction leaked. | L | GPU-003 |
| ☐ | **GPU-005** | `E-301 gpu_compute_mc` — asserts **bit-identical** to CPU, or documents the exact divergence. "Looks the same" is not an acceptance criterion. | M | GPU-004, B-003 |
| ☐ | **GPU-006** | `E-302 gpu_vs_cpu` — both live, timing HUD, resolution slider. Watch the gap **close** at small grids: launch overhead made visible. | M | GPU-005 |
| ☐ | **GPU-007** | **Mesh shader capability probe.** Print what this adapter reports for `EXPERIMENTAL_MESH_SHADER` and stop. **macOS/Metal is the unverified case** — wgpu's spec table says MSL is *planned*, the tracking issue says the Metal backend merged. Report the truth before writing a line of shader. | S | GPU-004 |
| ☐ | **GPU-008** | `E-303 gpu_mesh_shader` — feature-gated, off by default, graceful fallback, never panics on an unsupported adapter. | L | GPU-007 |

---

## Deliberately not in scope yet

Recorded so they don't get picked up early, and so it's clear they weren't forgotten.

- Nanite-style mesh-space cluster simplification — the research concludes it can't be repaired
  edit-proportionally (no local validity certificate). Field-derived LOD is the bet instead.
- Networked/concurrent editing — depends on G-003's commutativity result landing first.
- Neural / differentiable extraction (FlexiCubes, TetWeave) — different problem, different crate.
- Publishing to crates.io — **but reserve the name with a 0.0.0 placeholder early.** `megamesh` was
  taken 48 hours before we checked.
