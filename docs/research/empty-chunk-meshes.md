# Empty chunk meshes: published to the game, withheld from the renderer

**Issue #10** (Sanctuary, 2026-10-28): 102 of the 856 chunks of Sanctuary's destructible wall mesh to nothing, and each one costs a `WARN bevy_render::mesh::allocator: Mesh … contains no vertices.` line at every world build. Written 2026-10-30, loop iteration 1, before any code.

## Sources

| Kind | Stem | Where |
|---|---|---|
| Algorithm | `transvoxel_dissertation_lengyel2010` (E. Lengyel, *Voxel-Based Terrain for Real-Time Virtual Simulations*, UC Davis, 2010) | §3 p. 28, corpus lines 559-569: a block is a renderable node only when the surface crosses it |
| Shipping game | `transvoxel_dissertation_lengyel2010` | ch. 7 p. 99, corpus lines 3239-3252: the whole method shipped in the C4 Engine |
| Technical | `book_real_time_rendering_4th_akenine` (Akenine-Möller et al., *Real-Time Rendering*, 4th ed.) | pp. 817-818, corpus lines 23737-23775: per-object and per-call overhead the application should not issue |
| Technical, primary | Bevy 0.20 source in the cargo registry | `bevy_render-0.20.0/src/render_asset.rs:346-367`, `bevy_render-0.20.0/src/mesh/allocator.rs:431-436`, `bevy_pbr-0.20.0/src/material.rs:1097`, `bevy_pbr-0.20.0/src/render/mesh.rs:4705-4710`, `bevy_mesh-0.20.0/src/mesh.rs:44,244,703-711,985-1008` |

## Craft principle

Lengyel, p. 28:

> "When generating a triangle mesh for voxel-based terrain, we operate on each 16 x 16 x 16 block in isolation, but we must also have access to voxels belonging to neighboring blocks. If any part of the entire terrain mesh intersects a particular block, then that part becomes an individually renderable and individually cullable geometry node in the scene. Whenever a region of the voxel map is edited, all blocks intersecting the region's bounding box have their triangle meshes regenerated."

Two facts are separate in that passage, and this crate had them fused. *Every* block in an edit's bounding box is regenerated, surface or not: being meshed is a property of the block. Only a block the surface *intersects* becomes a renderable node: being drawn is a property of the result. `bevy_isomesh` publishes the first correctly (every finished chunk gets a `ChunkMesh`) and then asks the renderer to take the second for every chunk too, empty or not.

## How shipping games do it

The same dissertation, p. 99:

> "The full content of this dissertation has been implemented and deployed in the C4 Engine, a multi-platform computer game engine [Lengyel 2005]. The robustness of the techniques we have developed has been field tested and proven by a user base consisting of thousands of engineers, designers, and artists working on many types of applications."

So the p. 28 rule (a node only where the surface is) is shipped practice, not a proposal. *Real-Time Rendering* 4th ed., pp. 817-818, is the cost side of it:

> "Galeano [511] shows how ignoring such filtering and issuing unneeded state calls would have cost their WebGL application up to nearly 2 ms/frame."

> "The cost savings from consolidation are not just from avoiding API draw calls. There are also savings from the application itself handling fewer objects."

An empty chunk handed to the render world is exactly an object the application should not have handed over: Bevy extracts it, clones its (empty) attribute map, prepares a `RenderMesh`, and then the allocator rejects it with a warning.

## What Bevy 0.20 does with it, read from source (V)

- `render_asset.rs:346-367`: an asset is extracted only if `asset_usage.contains(RenderAssetUsages::RENDER_WORLD)`. Without that bit nothing crosses to the render world: no clone, no `RenderMesh`, no allocator visit.
- `mesh/allocator.rs:431-436`: for each extracted mesh, `if vertex_buffer_size == 0 { warn!("Mesh {:?} contains no vertices.", mesh_id); continue; }`. The warning's condition is **zero vertices**, not zero triangles.
- `bevy_pbr` `material.rs:1097` (`let Some(mesh) = render_meshes.get(..) else { continue }`) and `render/mesh.rs:4705-4710` (`RenderCommandResult::Skip` when the mesh or its vertex slice is absent): an entity whose `Mesh3d` names a mesh that was never extracted is skipped quietly. This answers the issue's open edge case ("a consumer that attaches `Mesh3d` to every `ChunkMesh` regardless") from source; the smoke run below confirms it on Metal.
- `bevy_mesh` `mesh.rs:44,703-711,985-1008`: `count_vertices()` and `try_attributes()` read the main-world attribute map, which a mesh keeps exactly when `MAIN_WORLD` is set. `mesh.rs:244`: `asset_usage` is a public field.

## Model adopted

One rule, at the point where this crate constructs a `Mesh`, in both constructors (`MeshBuilder::into_mesh` and `to_bevy_mesh`):

- zero vertices → `RenderAssetUsages::MAIN_WORLD`;
- otherwise → `RenderAssetUsages::default()` (`MAIN_WORLD | RENDER_WORLD`), unchanged.

The plugin is unchanged in shape: every finished chunk, empty or not, still gets a `ChunkMesh`, because "finished" is the block-was-meshed fact and a game reads it (occupancy, colliders, "is this chunk done"). A chunk that gains or loses its surface after an edit gets a new asset from the same rule, and the old asset is freed when its handle drops, as today.

**Assumed decision (owner may override)**, recorded on issue #10:

1. *The rule in both mesh constructors* (picked): the issue's own proposal; one place, every caller (plugin, examples, consumers calling `into_mesh` directly) inherits it; no API change.
2. A plugin-only rule (set the usage in `apply_finished_meshes`): fixes Sanctuary but leaves every direct `into_mesh` caller warning, which the issue names as an edge case to cover.
3. No `ChunkMesh` for an empty chunk: breaks the "every finished chunk has a `ChunkMesh`" contract the crate's own test states (`plugin/tests.rs`, "the chunk never finished").
4. A plugin option: API surface for a choice nobody has asked to make the other way.

## Parameters

None. The condition is Bevy's own (`vertex_buffer_size == 0`), so there is no threshold to choose. A mesh with vertices and no indices is not special-cased: Bevy neither warns nor fails on it, and no extractor in this crate emits one for a non-empty surface.

## Validation targets

1. After a volume with at least one empty chunk has meshed, no `ChunkMesh` asset has zero vertices and `RENDER_WORLD` (issue acceptance 1), and every non-empty one keeps `MAIN_WORLD | RENDER_WORLD`.
2. Every finished chunk still has a `ChunkMesh`; `an_edit_during_extraction_is_requeued_rather_than_swallowed` passes unchanged (acceptance 2).
3. An empty mesh reads: `count_vertices() == 0` without panicking, `try_attributes()` is `Ok`, `from_bevy_mesh` returns two empty vectors (acceptance 3).
4. A chunk that gains its surface gets a mesh with vertices and `RENDER_WORLD`; one that loses it gets an empty `MAIN_WORLD`-only mesh; the old asset is gone from `Assets<Mesh>` once its handle drops (acceptance 4).
5. On Metal, an example that streams chunks with empty ones among them and attaches `Mesh3d` to every `ChunkMesh` logs **0** `contains no vertices` lines after the change against a counted number before, and its frame shows the same terrain.
6. `crates/` is untouched, so the 378 golden hashes, the extract medians and the shootout are unchanged by construction; measured anyway.

## What great looks like

A world build's log carries no allocator warning caused by this crate, so the warnings a game does get are readable. The rendered frame is the same terrain as before, with no chunk missing and no stale geometry left behind when a wall is breached to nothing. An empty chunk is still an entity with a readable `ChunkMesh`, so a game's occupancy, collider or "is this chunk done" code sees "finished, nothing here" rather than "never finished". The rule is documented where a caller meets it, including the one way to defeat it (filling an empty mesh in place).

## Deviations

- **Empty meshes are `MAIN_WORLD` only**, against `RenderAssetUsages::default()` everywhere else. Justified by the allocator's condition (`bevy_render-0.20.0/src/mesh/allocator.rs:433`) and Lengyel p. 28. Non-empty meshes are unchanged.
- **Filling an empty mesh in place** (`Assets::get_mut` then `insert_attribute`) leaves it `MAIN_WORLD` only, so it will not draw. Replacing the whole `Mesh` with a new `into_mesh()` or `to_bevy_mesh()` result re-derives the usage. This is stated in both constructors' rustdoc rather than worked around: the crate cannot see a later in-place edit. Checked: no example adds vertices to an empty mesh from either constructor. `game_destruction` overwrites a fragment's positions with re-centred ones of the same length and returns early on an empty fragment (`examples/game_destruction.rs:497`), and the two examples that replace a mesh in place (`tape_pruning`, `game_edit_tape_trim`) assign a whole new `Mesh`.

## Pending index

None. All three source kinds are indexed in home-still (distill healthy, 2026-10-30); the Bevy source is read from the local cargo registry, not the corpus.
