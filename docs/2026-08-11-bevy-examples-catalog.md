# isomesh — Bevy examples catalog

**Date:** 2026-08-11
**Where they live:** `bevy_isomesh/examples/` — a workspace excluded from the root, so none of this
touches the core crate's dependency graph.
**Run:** `cd bevy_isomesh && cargo run --example <name> --release` (always `--release`).

26 examples in five tiers. Build them alongside the stage that enables them, not in a batch at the
end — an algorithm without its example is not done.

---

## Shared conventions — build this first

`bevy_isomesh/examples/common/mod.rs`, used by every example. Getting this right once is what makes 26
examples cheap instead of 26× expensive.

**Camera.** One orbit controller. LMB-drag orbit, RMB-drag pan, scroll zoom. Every example gets the
same one; nobody re-learns navigation.

**HUD.** A single top-left text block, one line per metric:
```
alg: dual_contouring   grid: 64³   verts: 12,443   tris: 24,102
extract: 3.21 ms       normals: 0.44 ms            total: 3.65 ms
non-manifold edges: 0  self-int/1k: 0.00
```
Metrics are pushed by the example; the HUD just renders whatever it's given. **Timings are rolling
medians over 30 frames, not instantaneous** — instantaneous frame times are noise and will make you
chase ghosts.

**Universal keybindings.** Same in every example, no exceptions:

| Key | |
|---|---|
| `W` | wireframe toggle |
| `N` | normals-as-gizmos toggle |
| `G` | grid/cell overlay toggle |
| `Space` | pause / resume any animation |
| `R` | re-mesh now |
| `F12` | screenshot to `screenshots/<example>-<n>.png` |
| `Esc` | quit |

**Test fields.** One shared `fields.rs` so comparisons are apples-to-apples: `sphere`, `torus`,
`box_exact` (sharp edges — the Dual Contouring discriminator), `csg_difference` (box minus sphere, concave sharp
edges), `fbm_terrain`, `thin_plate` (sub-voxel feature), `gyroid` (high genus, stresses topology).

**Screenshot on `F12`** is not a nicety. Several examples below are *visual* acceptance tests, and the
before/after pair belongs in the commit.

---

## Tier 1 — One per algorithm (stages 1–4, 6)

| # | Example | Demonstrates | What you see |
|---|---|---|---|
| 1 | `marching_cubes_sphere` | Marching Cubes baseline | A sphere. `W` shows the characteristic Marching Cubes triangle pattern with varying triangle sizes per cell. The reference every other example is judged against. |
| 2 | `marching_cubes_ambiguity` | Ambiguous face resolution | **Built, and re-specified — the row below is falsified by ✗11 and is kept for the record.** Neither side holes and neither can: the case table is derived by walking each face counter-clockwise, so two cells sharing a face cannot disagree about it. What ships instead is a split screen of `FaceAmbiguity::Separate` against `::AsymptoticDecider`, with every ambiguous cell boxed — amber where the decider agreed, magenta where it joined — and the HUD reporting χ for each side plus whether the two meshes are byte-identical. It has to run on `gyroid` or `fbm_terrain`; on the other five reference fields an ambiguous face never occurs and the split screen would show two identical meshes (M-40). |
| | | *Original spec, falsified:* Split screen. Left: plain Marching Cubes, with **visible holes** on the ambiguous configuration. Right: Marching Cubes 33, closed. HUD shows the Euler characteristic for each — 2 on the right, not 2 on the left. |
| 15 | `chunk_seam_weld` | Chunk seam, and A-013 welding | Two chunks of one torus, meshed independently, with the boundary edges on the shared plane drawn red and the chunks' own outer faces amber — the distinction `validate` makes between `boundary_edges` and `non_manifold_edges`, made visible. `V` welds. The spacing selector is the point: `h = 0.125` is bit-exact at a seam and `h = 4/35` is not (M-32), and one epsilon covers both. **This is the only visual for Phase 2**, which is otherwise three tickets of chunking with nothing to look at. |
| 202 | `game_dig` | The core carve loop | First person, click to carve. Re-meshes only the chunks G-002's `mark_edit` marks dirty, and outlines them. **The HUD carries E1 live** — cells in the brush box that actually re-mesh against cells whose sample merely moved — which is the number the incremental story rests on and had only ever been measured offline. `ISOMESH_AUTOCARVE=n` drives the loop without a mouse, because a screenshot cannot click. |
| 3 | `surface_nets_sphere` | Surface Nets vs Marching Cubes | Same field, side by side, with triangle counts. Surface Nets' dual mesh is visibly more regular under `W`. |
| 4 | `dual_contouring_cube` | **Sharp features** | `box_exact` and `csg_difference`. Surface Nets rounds the corners into mush; Dual Contouring holds them crisp. This is the single most persuasive image the project will produce — put it in the README. |
| 5 | `marching_tetrahedra` | MT and its cost | Same sphere. HUD triangle count is dramatically higher than Marching Cubes on the identical field. Shows why MT is chosen for guarantees, not for budget. |
| 6 | `subgrid_mt` | **Sub-voxel features** | A stone wall with letters carved into it, where the stroke width is *thinner than one voxel*. Classic marching cannot resolve this at all — toggle between the two and the text vanishes. |
| 7 | `transvoxel_seams` | LOD transitions | Two adjacent chunks at different resolutions. Toggle transition cells: cracks appear, cracks vanish. Screenshot both. |
| 8 | `greedy_quads` | The blocky path | Minecraft-style output on the same field, for the budget axis of the tradeoff table. Quad count before/after merging. |

---

## Tier 2 — The tradeoff axes (stage 4)

These are why someone picks one algorithm over another. Each is a live A/B, not a static render.

| # | Example | Demonstrates | What you see |
|---|---|---|---|
| 9 | `sharp_features` | The sharpness threshold | One model, a slider for the normal-deviation threshold. Watch edges snap from rounded to sharp, and watch it over-sharpen into spikes past the useful range. Teaches the parameter better than any doc. |
| 10 | `qef_clamp` | **The cell clamp** | Toggle the (1−ε) clamp. HUD shows **self-intersections per 1,000 triangles** updating live, and offending triangles highlight in red. Settles whether guaranteed intersection-free extraction is available for free. |
| 11 | `manifold_check` | Topological soundness | Non-manifold edges drawn as thick red lines, non-manifold vertices as red spheres. Run it across all six test fields and all four algorithms. |
| 12 | `precision_f32_vs_f64` | Why CAD needs f64 | The same field evaluated at large coordinate offsets (e.g. 1e6). f32 visibly cracks and jitters; f64 doesn't. HUD shows the QEF matrix condition number. |
| 13 | `normal_estimation` | Three normal strategies | Analytic gradient vs central differences vs area-weighted face normals, three panels, lit. Differences show up in specular highlights on curved regions, not in wireframe. |
| 14 | `hermite_debug` | What Dual Contouring actually operates on | Gizmo view: edge crossings as dots, surface normals as lines, the solved cell vertex as a larger dot, the cell as a wire box. The debug view you'll live in while stage 4 is broken. |

---

## Tier 3 — The engine-shaped workload (stage 5)

The examples that prove the crate is usable for a real game, not just a paper.

| # | Example | Demonstrates | What you see |
|---|---|---|---|
| 15 | `chunked_terrain` | Chunking, seams | fBm terrain across many chunks. `G` overlays chunk bounds. **Zero visible cracks** — hunt for them at grazing angles, that's where they hide. |
| 16 | `edit_brush` | Carving | Click to add/subtract material. HUD: chunks re-meshed this stroke, ms spent, chunks touched vs chunks in world. The core interaction loop of the game. |
| 17 | `dirty_set_metrics` | **E1 from the research** | Same as 16, plus a histogram of *what fraction of cells actually changed* per stroke. This number is the ceiling on every incremental-repair idea in the opportunities doc, and it is unpublished. Log it to CSV. |
| 18 | `lod_field_derived` | Field-derived LOD | Mip the field, mesh at 4 LOD levels, fly out and back. No popping, no cracks. Demonstrates the architectural bet — derive LOD from the field, not the mesh. |
| 19 | `collider_roundtrip` | Physics-grade output | Carve an arbitrary shape → mesh → `parry3d::TriMesh` → drop 200 rigid bodies on it. Do they rest, or fall through? Carve a spiral and a hollow shell; those are where naive decomposition fails. |
| 20 | `game_paint` | **Paint survives destruction** | Spray colour on a wall, then blow a hole through it. The paint on the remaining wall must be exactly where you sprayed it — not smeared, not reset. Row 4 of the opportunities table, and the cheapest genuinely-new thing in it. **Shipped at E-208 as `game_paint`**, not `attribute_transfer`: there is no transfer in it. Paint lives in the edit log, so a carve cannot move it and the L²-nearest machinery this row was priced on is unnecessary rather than cheap (M-137). |
| 21 | `csg_live` | Boolean ops | Union / difference / intersection of moving primitives, re-meshed per frame. The CAD-facing demo, and a stress test for concave sharp edges. |

---

## Tier 4 — Measurement (stage 7)

| # | Example | Demonstrates | What you see |
|---|---|---|---|
| 22 | `bench_shootout` | **The comparison that doesn't exist** | Marching Cubes / Marching Cubes 33 / Surface Nets / Dual Contouring / MT on identical fields and grids. On-screen table: ms, verts, tris, non-manifold edges, self-int per 1k, Hausdorff error vs analytic. Writes CSV to `docs/measurements/`. |
| 23 | `bench_resolution_sweep` | Fixed cost vs marginal cost | Time vs grid size, 16³ → 256³, plotted live. Fits `t = a + b·n³` and **prints `a`**. In the published FlexiCubes numbers 73% of the 64³ time was fixed overhead — expect to find the same and to stop trusting single-grid comparisons. |
| 24 | `bench_stage_breakdown` | Where the time actually goes | Stacked bar: contour / normals / weld / collider build / GPU upload. Grosso & Zint measured contouring at 68 ms against 58 ms of halfedge construction — **the contour is 54% of a usable mesh**. Find out what yours is before optimizing anything. |

---

## Tier 5 — GPU (stages 8–9)

| # | Example | Demonstrates | What you see |
|---|---|---|---|
| 25 | `gpu_compute_mc` | Compute-shader extraction | Same field, GPU path, readback. HUD asserts **bit-identical to CPU** or reports the exact divergence. Not "looks the same." |
| 26 | `gpu_vs_cpu` | The real speedup | Both paths live, side by side, timing HUD, resolution slider. Watch the gap open as the grid grows — and watch it *close* at small grids, which is the launch-overhead finding made visible. |
| 27 | `gpu_mesh_shader` | Feature-gated, exploratory | **First job is the capability probe**: print what the adapter reports for `EXPERIMENTAL_MESH_SHADER` and stop. On macOS/Metal this is unverified upstream. If unsupported, the example must print why and fall back cleanly — never panic. |

---

## Build order

Ship each example with the stage that enables it:

| Stage | Examples |
|---|---|
| 0 | `common/` module, `fields.rs`, orbit camera, HUD |
| 1 | 1 |
| 2 | 2 |
| 3 | 3 |
| 4 | 4, 5, 9, 10, 11, 12, 13, 14 |
| 5 | 15, 16, 17, 19, 20, 21 |
| 6 | 6, 7, 8, 18 |
| 7 | 22, 23, 24 |
| 8 | 25, 26 |
| 9 | 27 |

---

## Rules for examples

1. **Every example is self-contained and runnable with one command.** No setup steps, no asset
   downloads, no "first run the other example." Fields are generated in code.
2. **Every example displays its own metrics.** An example that renders something pretty but reports
   nothing is a screensaver.
3. **Comparison examples run both paths simultaneously**, in the same process, on the same frame.
   Sequential A/B across two runs is not a comparison — hardware state drifts.
4. **No example may `unwrap()` on GPU capability.** Probe, report, degrade.
5. **Examples are load-bearing tests.** `cargo build --examples` runs in CI. This is deliberate:
   `block-mesh`'s examples rotted to bevy 0.13 and `fast-surface-nets`' to bevy 0.7 precisely because
   nothing in CI ever compiled them.
6. **A visual acceptance test needs committed screenshots.** #2, #4, #6, #7, #15 and #20 are all
   "toggle a thing, look at the difference" — commit the before/after pair so a regression is visible
   in a diff rather than requiring someone to run the example.
