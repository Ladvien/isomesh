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
| ☑ | **I-001** | **Workspace skeleton.** Root `Cargo.toml` with `members = ["crates/*"]`, `exclude = ["bevy_isomesh"]`. `crates/isomesh` with `libm 0.2` as its only dependency, edition 2024, `#![no_std]` + `extern crate alloc`, unconditionally. **Acceptance:** `cargo tree -p isomesh -e normal` prints exactly two packages. | S | — |
| | | *Amended during implementation.* Three changes, each with its reason. **(a)** `-e normal` added: cargo's default edges are `normal,build,dev`, so the bare command starts showing proptest's and criterion's trees the moment T-005 and T-006 land — it does not express the intended claim today either. **(b)** `glam` → `libm` as the sole dependency: `core` has no `sqrt`/`floor`/`sin`/`cos` on stable and `Real` needs all four, whereas nothing in Phase 0 needs vector types — fields are scalar and `validate.rs` wants one hand-written cross product. glam is still the sanctioned internal math library and lands at A-007. **(c)** The `std` default feature dropped: `libm` is used unconditionally so there is one float backend and one set of results, which also makes T-007's golden hashes portable across macOS and Linux. With nothing left for `std` to gate, keeping it would be a forward-declared stub. Justifications recorded in `CLAUDE.md`. | | |
| ☑ | **I-002** | **Core traits.** `Real` (sealed, f32/f64), `Sdf`, `Shape3`, `MeshSink`, `MeshBuffer`. Signatures per the brief — **plain arrays only, no glam in any public item.** `MeshBuffer::reset()` clears without shrinking. **Acceptance:** a test asserts capacity survives `reset()`; a test asserts `Sdf` is object-safe if you intend `dyn` use, or documents that it isn't. | M | I-001 |
| | | *One deviation.* `MeshSink` carries `type Scalar: Real` and `MeshBuffer` is `MeshBuffer<R: Real = f32>`, where the brief writes `[f32; 3]` literally. The brief also says twelve lines earlier that `fast-surface-nets`' `SignedDistance: Into<f32>` bound "forecloses f64" for CAD and not to repeat it — an f32-only sink repeats exactly that foreclosure at the output end, throwing an f64 field and an f64 QEF solve through f32 on the way out. The common spelling is unchanged (`MeshBuffer`), `dyn` use becomes `dyn MeshSink<Scalar = f32>`, and narrowing moves downstream to `bevy_isomesh` where GPU buffers are f32 anyway. Cheap now; it would break every public signature later. | | |
| ☑ | **I-003** | **Shared test fields** in `crates/isomesh/src/fields.rs`: `sphere`, `torus`, `box_exact`, `csg_difference` (box − sphere), `fbm_terrain`, `thin_plate`, `gyroid`. Each with an analytic gradient where one exists. **These are shared by tests, benches, and every Bevy example** — one definition, no drift. **Acceptance:** each field has a test asserting sign at a known inside and outside point. | M | I-002 |
| | | *Scope added, deliberately.* **(a)** A `ReferenceField` trait (`NAME`, `domain`, `closed_in_domain`, `expected_euler`, `is_exact_distance`). Without it the harness needs a per-field `if` ladder, because `gyroid` is periodic (a finite box cuts it, so it has boundary) and `fbm_terrain` is a heightfield (open by construction) — neither has an Euler characteristic assertable a priori, and inventing one would violate rule 5. **(b)** The canonical gyroid is capped to a sphere so it is closed at all, which needs an `Intersection` combinator alongside `Difference`. **(c)** `for_each_reference_field!` sweeps all seven with concrete types; `Vec<Box<dyn Sdf>>` would put a vtable on the innermost loop of every benchmark. | | |
| | | *Two numbers worth recording.* No document in this repo contains a single field formula, so all seven originate here, each cited in its doc comment. `fbm_terrain`'s noise is hand-rolled (rule 3 forbids a noise crate) and evaluates **no transcendental at all** — integer hash plus a polynomial — so it is bit-identical across platforms, which is what T-007's committed goldens need. Its amplitude bound is the **provable** `\|n\| ≤ 2` rather than the widely quoted `√3/2 ≈ 0.866`; measured max is **0.907716** over 43,200 samples, i.e. the folklore figure is *wrong for this implementation* and hard-coding it would have been a silent bug. | | |
| ☑ | **T-001** | **Mesh validity harness** — `crates/isomesh/src/validate.rs`, compiled in normal builds (consumers want it too). Reports, doesn't just assert: `euler_characteristic`, `non_manifold_edges` (edge used by ≠2 faces), `non_manifold_vertices`, `boundary_edges`, `degenerate_triangles`, `out_of_range_indices`, `duplicate_vertices`. Returns a `MeshReport` struct that `Display`s as one block. **Acceptance:** unit-tested against a hand-built tetrahedron (χ=2, 0 violations) and a deliberately broken mesh (each violation type detected). | M | I-002 |
| | | *Definitional choices, each a judgment call rather than a fact.* **(a)** `non_manifold_edges` counts ≥3 faces and `boundary_edges` exactly 1, splitting the ticket's literal "≠2" — lumping them together double-counts and makes zero unachievable for any open mesh, i.e. for `fbm_terrain` and for every individual chunk. **(b)** `V` in χ is *referenced* vertices, with `unreferenced_vertices` reported separately; otherwise a stale vertex pool reports χ=3 for a sphere and reads as a topology bug when it's an allocation bug. **(c)** `duplicate_vertices` is "has an earlier vertex within ε" — epsilon-closeness isn't transitive, so equivalence classes aren't well defined but this is. **(d)** `genus` is `Some` only for a single, consistently oriented, manifold component with no skipped faces; the orientation precondition is what makes `χ = 2−2g−b` sound. **(e)** No `Result`, no panic, no short-circuit: a `Result` forces a branch at every call site and `?` discards the partial report exactly when it's most useful. | | |
| | | *Three counters added beyond the list, each because its absence is a hole.* `inconsistently_oriented_edges` — a case table with one flipped triangle passes χ, edge-manifoldness *and* vertex-manifoldness while being inside out (fixture 10 demonstrates it). `boundary_loops` — `χ=1` alone can't distinguish a hole from a handle, and it's `b` not the edge count that enters `χ = 2−2g−b`. `non_finite_positions` — a NaN vertex quantises into one bucket, never matches as a duplicate, and `NaN <= threshold` is false so it isn't degenerate either; without this it passes silently. Also: `degenerate_triangles` is excluded from `violations()` and is a recorded metric, since MC emits slivers for ordinary reasons. | | |
| ☐ | **T-002** | **Self-intersection counter.** Triangle–triangle intersection test, BVH or uniform grid broadphase. Reports intersections per 1,000 triangles. Doesn't need to be fast; needs to be right. **Acceptance:** two hand-placed crossing triangles detected; two adjacent triangles sharing an edge **not** counted (this is the trap). | M | T-001 |
| ☐ | **T-003** | **Accuracy harness.** One-sided and symmetric Hausdorff distance from the mesh to an analytic SDF, sampled on triangle centroids and vertices. Also mean absolute error. **Acceptance:** a unit sphere meshed at 64³ has max error below one cell diagonal. | M | I-003 |
| ☐ | **T-004** | **Determinism harness.** Run any extractor twice, assert byte-identical output buffers. **Acceptance:** wired into the test suite as a helper that every algorithm ticket calls. Catches `HashMap` iteration order leaking into vertex order — which it will. | S | T-001 |
| ☐ | **T-005** | **Property-test scaffolding.** `proptest` over randomized fields (random sphere unions, random planes) × grid sizes × algorithms. Every generated mesh must pass T-001 with zero violations. **Acceptance:** 1,000 cases green; a deliberately corrupted case table makes it fail. Verify that second part — a property test that can't fail is decoration. | M | T-001, I-003 |
| ☐ | **T-006** | **Benchmark harness.** `criterion`, plus a resolution sweep that fits `t = a + b·n³` and **reports `a` separately**. CSV out to `docs/measurements/`. Rationale in the speed analysis: 73% of a published 64³ figure was fixed launch overhead. **Acceptance:** sweep runs 16³→256³ and prints the fitted fixed cost. | M | I-003 |
| ☐ | **T-007** | **Golden-hash regression.** Stable hash of (positions, normals, indices) for each (algorithm, field, resolution). Committed as a JSON fixture. **Acceptance:** a deliberate one-bit change to a case table fails the test with a useful message naming which combination drifted. | S | T-004 |
| ☐ | **I-004** | **CI.** GitHub Actions: fmt, clippy `-D warnings`, `cargo test` in the root workspace, **and `cd bevy_isomesh && cargo build --examples`**. That last one is the whole reason examples live where they do — `block-mesh`'s rotted to bevy 0.13 because nothing ever compiled them. **Acceptance:** CI green on a clean checkout. | S | I-001 |

---

## Phase 1 — The usual suspects

Each algorithm ticket is done when: it passes T-001 with zero violations on all seven test fields, at
three resolutions; T-004 determinism passes; T-005 covers it; and a benchmark exists.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **A-001** | **Marching Cubes.** Standard 256-case table. **Also write the structural validator** described in the brief: for all 256 configurations, generated triangles' boundary edges must lie on cube faces and be consistent with that face's corner signs. Cite the table's source in a comment. Document the sign convention (**negative = inside**) and winding order on the public fn. | L | T-001, T-005 |
| ☐ | **A-002** | **MC33 / asymptotic decider.** Face saddle `S = (f00·f11 − f10·f01) / (f00 + f11 − f10 − f01)`, guard the denominator. Interior ambiguity via the trilinear body saddle. **Acceptance:** a field that holes under A-001 is closed under A-002, and the test asserts *both* — one that passes with plain MC isn't testing this. | L | A-001 |
| ☐ | **A-003** | **Marching Tetrahedra.** 6-tet decomposition of the cube (document which one — the choice affects the output). Unambiguous by construction, so this is the topological reference the others get compared against. Expect a large triangle count; record it. | M | A-001 |
| ☐ | **A-004** | **Surface Nets.** One vertex per sign-changing cell at the centroid of edge crossings, quad or triangle output. Optional Laplacian smoothing passes as a parameter. **Note for the record:** this has no credible published timings anywhere despite being what engines ship — your benchmark will be the reference. | M | T-001 |
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
| ☐ | **B-001** | **Crate skeleton.** Own workspace (excluded from root). Depends on **leaf crates** — `bevy_app`, `bevy_ecs`, `bevy_asset`, `bevy_mesh` — with `bevy_render` **optional** behind a `gpu` feature, so a CPU-only consumer never compiles the renderer. Pin `bevy 0.19`. | S | I-002 |
| ☐ | **B-002** | **`MeshSink` → Bevy `Mesh`.** Write straight into attribute arrays, no intermediate copy. Correct `Indices::U32`, `PrimitiveTopology::TriangleList`, positions/normals/UVs. | M | B-001 |
| ☐ | **B-003** | **Plugin and component API.** `IsomeshPlugin`, a `VoxelVolume` component, a `NeedsRemesh` marker, systems that consume the frame budget from G-006 and drive `AsyncComputeTaskPool` so meshing is off the main thread. **Acceptance:** meshing a large volume does not stall the render loop — show it in the frame-time graph. | L | B-002, G-006 |
| ☐ | **B-004** | **`examples/common/`.** Orbit camera, HUD (rolling 30-frame medians, not instantaneous), universal keybindings (`W` wire, `N` normals, `G` grid, `Space` pause, `R` remesh, `F12` screenshot, `Esc` quit), field picker. Build this **before** any example — it's what makes 27 examples cheap. | M | B-002 |

---

## Phase 4 — Examples

Two groups. The algorithm demos are quick and prove correctness visually. **The game-shaped ones are
the point** — they're how someone decides whether this crate is usable.

### 4a — Algorithm demos

| | ID | Example | Blocked by |
|---|---|---|---|
| ☐ | **E-101** | `mc_sphere` — MC baseline, wireframe toggle | A-001, B-004 |
| ☐ | **E-102** | `mc33_ambiguity` — split screen, holes vs closed, Euler χ in the HUD for both | A-002 |
| ☐ | **E-103** | `surface_nets_sphere` — SN next to MC, triangle counts | A-004 |
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
