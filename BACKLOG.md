# isomesh — BACKLOG

**Updated:** 2026-08-12
**Companions:** `CLAUDE.md` (rules), `FINDINGS.md` (what we know and how well),
`BACKLOG_ARCHIVE.md` (completed tickets + why they changed),
`docs/2026-08-11-implementation-brief.md` (the how),
`docs/2026-08-11-bevy-examples-catalog.md` (example detail), `docs/research/` (the why).

**36 tickets archived, 40 open.** Completed rows move to `BACKLOG_ARCHIVE.md` with their amendments
attached — read that before re-litigating a decision this project already made.

---

## How to work this backlog

1. Take the **topmost unblocked, unchecked ticket**. Don't cherry-pick interesting ones — the order
   encodes dependencies and the test harness exists so later work is cheap.
2. One ticket = one commit (or a short stack). Commit message starts with the ticket ID.
3. **Check the box in this file as part of that same commit.** This file is the state.
4. If a ticket can't be finished, leave it unchecked, add a `> BLOCKED:` line under it saying exactly
   what's in the way, and move to the next unblocked ticket. Do not half-finish and check the box.
5. If a ticket turns out to be wrong or to need splitting, edit it and say so in the commit.
6. **On completion, move the row to `BACKLOG_ARCHIVE.md`** with an indented annotation recording any
   amendment, deviation or falsified premise. The annotation is the point; the checkmark is not.
7. **New tickets** slot in by dependency, not by number. A ticket split after the fact takes a letter
   suffix (`T-005a`/`T-005b`, `A-002`/`A-002b`); a genuinely new one takes the next free number in its
   series even if that puts it out of numeric order (`A-015`).

### Definition of done — applies to every ticket

- Code compiles with no warnings. `cargo clippy -- -D warnings` clean.
- **`RUSTDOCFLAGS=-D warnings cargo doc --workspace --no-deps` clean.** A third of the lint job, and
  the third that clippy and fmt do not cover — it is what catches a doc link pointing at nothing, or
  at a private item. Added at A-002, which shipped one and found out from CI.
- `cargo test -p isomesh` green. `cargo fmt` applied.
- **`grep -ri bevy crates/` returns nothing.** Non-negotiable — see `CLAUDE.md` rule 2.
- Any algorithm ticket also passes the T-001 validity suite. No exceptions, no "I'll add tests later."
- Any perf claim has a committed benchmark that produced it.
- Public items have doc comments. Anything with a sign convention, a coordinate order, or a winding
  order says so **in the doc comment**, not in a code comment.
- **`FINDINGS.md` updated in the same commit** if the ticket measured something, contradicted
  something written down, or earned a method rule. A measurement that only exists in a commit message
  is not retrievable six weeks later.

**Size key:** `S` ≈ one sitting · `M` ≈ a day · `L` ≈ multi-day, consider splitting.

---

## Phase 0 — Foundation and the test harness ✅ complete

All eleven tickets archived (I-001..I-004, T-001..T-008). The bet paid: every algorithm since has
been cheap to validate because the harness predated it.

---

## Phase 1 — The usual suspects

Each algorithm ticket is done when: T-001 reports **no unexplained violations** on all seven test
fields at three resolutions; T-004 determinism passes; T-005 covers it; and a benchmark exists.

> **Amended 2026-08-12.** This originally said "zero violations on all seven fields." M-4 falsified
> that as a universal gate: Surface Nets is *legitimately* non-manifold where one cell carries two
> sheets (48 edges on capped gyroid, 15 on fbm_terrain), and those counts are pinned as **non-zero
> assertions** precisely so they can't drift silently. A known defect with a pinned number and a
> ticket that owns it satisfies this gate. An unexplained one does not.
>
> **A-010 has since landed and closed that owner.** Surface Nets' and Dual Contouring's counts stay
> pinned as non-zero — they are properties of one-vertex-per-cell, not bugs — and
> `manifold_dual_contouring` is the entry that takes the zero. Its own residue is M-59's parallel-edge
> collapse, pinned at one edge on the ✗15 fixture and zero everywhere else, owned by O-16.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **A-011** | **Transvoxel transition cells.** Half-resolution transition cells at LOD boundaries. **Acceptance:** two adjacent chunks at differing LOD produce zero boundary gaps — assert on the geometry, then confirm visually in E-107. | L | A-001 |
| ☐ | **A-012** | **Normal estimation strategies.** Analytic gradient / central differences / area-weighted face normals, selectable. **Acceptance:** all three produce unit-length normals; analytic and central-difference agree within tolerance on `sphere`. | S | A-001 |
| ☐ | **A-014** | **Subgrid Marching Tetrahedra.** Integer edge-intersection counts instead of vertex signs — resolves features finer than one cell. The differentiator; do it after the usual suspects are solid. | L | A-003 |
| ☐ | **A-002b** | **Marching Cubes 33 interior ambiguity — the trilinear body saddle.** Deliberately deferred at A-002, on evidence, not forgotten. Three reasons. (1) `catalog-v2.md:107` is explicit: *"Skip the interior test; spend the budget on chunk seams"* — a game needs topological *consistency*, which A-001 already has, over *correctness*. (2) **There is no correct published table to transcribe.** Custodio et al. 2013 (`10.1016/j.cag.2013.04.004` §5.1) prove Chernyaev's interior test tracks a quadratic where the true saddle trajectory is hyperbolic with an asymptote, so case 13.5.2 is misread as 13.5.1 — counterexample values in their Appendix A — and Lewiner's reference implementation omits disambiguation for cases 10 and 12 entirely. Rule 5 forbids inventing the missing one. (3) The v1 catalog prices it: the decider is *"~free"*, the guaranteed version is **730 subcases in the LUT**. Also needs cell-interior vertices for tunnels, which the grid-edge-keyed vertex cache has no slot for. **Acceptance:** a cell where the body saddle says "tunnel", meshed as a tunnel, with the sign tracked by Custodio's correction rather than Chernyaev's `F(t)`. | L | A-002 |

---

## Phase 2 — Game-shaped infrastructure

Still zero Bevy. This is the machinery a game needs, living in the core crate where CAD can use it too.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **G-004** | **Field-derived LOD.** Mip the field (not the mesh), mesh at level N. **Acceptance:** LOD 0..3 all mesh cleanly; LOD *k* has roughly 1/8^k the cells. Pairs with A-011 for the seams. | M | A-011, G-001 |
| ☐ | **G-005** | **Collider export.** `MeshBuffer` → `parry3d::TriMesh` (`Vec<[u32;3]>` indices — parry takes plain arrays). Behind an optional `parry` feature. Optional convex decomposition path. **Acceptance:** a carved shape builds a `TriMesh` without error and passes parry's own validity check. **A-013 landed, so this is unblocked** — and note a chunked collider must be welded first or parry sees a seam of unshared vertices. | M | A-013 |
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
| ☐ | **E-105** | `marching_tetrahedra` — same field, much higher triangle count | A-003 |
| ☐ | **E-106** | `greedy_quads` — the blocky path, quads before/after merge | A-005 |
| ☐ | **E-107** | `transvoxel_seams` — two LODs adjacent, toggle transition cells. **Commit both screenshots.** | A-011 |
| ☐ | **E-108** | `subgrid_features` — letters carved thinner than a voxel; toggle and watch them vanish | A-014 |
| ☐ | **E-109** | `sharp_features` — live slider on the normal-deviation threshold, through to over-sharpening | A-007 |
| ☐ | **E-110** | `qef_clamp` — clamp toggle, live self-intersections/1k, offending triangles in red | A-009 |
| ☐ | **E-112** | `precision_f32_vs_f64` — same field at ~1e6 offsets; f32 cracks, f64 doesn't. Condition number in the HUD. | A-007 |
| ☐ | **E-113** | `normal_estimation` — three panels, lit; differences live in the speculars, not the wireframe | A-012 |
| ☐ | **E-114** | `hermite_debug` — crossings, normals, solved vertex, cell box as gizmos. The view you debug A-007 in. | A-006 |

### 4b — Game-shaped

These use the algorithms the way a game does: chunked, edited, budgeted, collided against.

| | ID | Example | What it has to prove | Blocked by |
|---|---|---|---|---|
| ☐ | **E-201** | `game_terrain_stream` — walk a large fBm world, chunks stream by distance | Sustained 60 fps while streaming. HUD: chunks resident, meshing ms/frame, MB. | G-007, B-003 |
| ☐ | **E-203** | `game_walk` — character controller on meshed terrain, parry3d colliders | **The acid test.** Walk every chunk seam. No falling through, no invisible walls. If this fails, G-001's overlap is wrong. | G-005, E-201 |
| ☐ | **E-204** | `game_destruction` — shoot a wall, it craters, debris becomes rigid bodies | Runtime fragments are correct physics bodies, not pre-fractured props. Carve a spiral and a hollow shell — that's where decomposition fails. | G-005, G-003 |
| ☐ | **E-205** | `game_lod_flyover` — fly out and back across LOD transitions | No popping, no cracks, no hitching. Transvoxel doing its job at speed. | G-004, A-011 |
| ☐ | **E-206** | `game_budget` — a deliberately overloaded edit queue under a frame budget | Frame time stays flat while the backlog drains. **Amortized cost per frame is the number no paper measures and the only one a game cares about.** | G-006 |
| ☐ | **E-207** | `game_editor` — sculpt with brushes, undo/redo over an op log | The CAD/editor use case. Undo is a re-fold of the log, not a snapshot. | G-003 |
| ☐ | **E-208** | `game_paint` — spray colour on a wall, then blow a hole through it | Paint on the remaining wall is exactly where you sprayed it. Row 4 of the opportunities table. | G-003, B-002 |
| ☐ | **E-209** | `game_csg_props` — place and boolean CSG primitives into the world live | Re-mesh per frame under moving primitives; concave sharp edges hold up. | G-003, A-007 |

---

## Phase 5 — Measurement

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **M-002** | **`bench_resolution_sweep` — the live plot only.** ~~fits `t = a + b·n³`, **prints `a`**~~ **The numeric half is done and its premise is falsified (M-62).** The bench had been printing `NaN` since the algorithm names were spelled out; fixed, and the answer is the opposite of what this ticket expected — Marching Cubes' `a` is **0.5118 ms, 0.64% of the largest run**, so there is *no* large fixed cost, and Surface Nets and Dual Contouring fit `a < 0`, meaning the model does not describe them. What is left is the **live plot**, which is a Bevy example and belongs in Phase 4. | S | T-006 |
| ☐ | **M-003** | **`bench_stage_breakdown`** — stacked bar: contour / normals / weld / collider / upload. Published comparison: contouring 68 ms vs halfedge construction 58 ms — **the contour was 54% of a usable mesh.** Find your ratio before optimizing anything. | M | G-005, M-001a |
| ☐ | **M-004** | **Write up M-001a..M-003** as `docs/research/YYYY-MM-DD-measured-comparison.md`. Numbers, method, hardware, and what surprised you. This is publishable on its own. | S | M-003 |

---

## Phase 6 — GPU (do not start before Phase 5)

The speed analysis is explicit that stage placement dominates the extraction algorithm by roughly an
order of magnitude. Which means GPU work is worth doing — and worth doing *after* you know your own
numbers, or you won't be able to tell what the port bought you.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **GPU-001** | `isomesh-gpu` skeleton, `wgpu 29.0.3`. Public API takes `&wgpu::Device` / `&Queue` / `&mut CommandEncoder`. **Never a Bevy type.** | M | M-001a |
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
- Networked/concurrent editing — **precondition resolved.** G-003 landed with commutativity measured
  over all 40,320 orderings; record the verdict here and promote this to a real ticket or close it out
  explicitly. Leaving it as a pending "depends on" after the dependency has landed is how a decision
  gets silently dropped.
- Neural / differentiable extraction (FlexiCubes, TetWeave) — different problem, different crate.
- Publishing to crates.io. **But `I-005 — reserve the name` is now overdue, not deferred.** Publish a
  `0.0.0` placeholder. `megamesh` was taken 48 hours before we checked it; `isomesh` has been sitting
  unreserved for a day with a public repo pointing at it. Ten minutes, unbounded downside.
