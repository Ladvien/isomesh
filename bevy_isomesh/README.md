# bevy_isomesh

> ⚠️ **Vibe Coded.** This crate was written by an AI agent working from a research corpus and a ticket queue. Every number below is produced by a test in this repository and every algorithm cites its source, but it has not been through human code review. Read the tests before trusting it with anything that matters.

[![crates.io](https://img.shields.io/crates/v/bevy_isomesh.svg)](https://crates.io/crates/bevy_isomesh) [![docs.rs](https://img.shields.io/docsrs/bevy_isomesh)](https://docs.rs/bevy_isomesh) [![CI](https://github.com/ladvien/isomesh/actions/workflows/ci.yml/badge.svg)](https://github.com/ladvien/isomesh/actions/workflows/ci.yml) [![license](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](https://github.com/ladvien/isomesh/blob/main/LICENSE-MIT)

**Bevy 0.19 integration for [`isomesh`](https://crates.io/crates/isomesh): signed distance field in, `Mesh` asset out.** The core crates stay engine-agnostic — this is the one place Bevy types appear.

![Flying through caves and arches meshed from a signed distance field](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/screenshots/e210-showcase-hero.png)

*Caves, arches and a roof over your head, from a nine-line field — the whole world is `max(heightfield, |gyroid| − thickness)`, streamed and meshed off the main thread while the camera flies through it. `cargo run --example game_showcase --release`.*

## Compatibility

| `bevy_isomesh` | `bevy` | `wgpu` |
|---|---|---|
| 0.0.x | 0.19 | 29.x |

Bevy 0.19 pins `wgpu` 29.0.3 and `glam` 0.32, and those move together or not at all: Cargo resolves two `wgpu` majors side by side with no error, and the failure surfaces much later as `expected TextureFormat, found a different TextureFormat`.

## Quickstart

Spawn a volume, mark its chunks, and attach the render component when the mesh arrives:

```rust,no_run
use bevy::prelude::*;
use bevy_isomesh::{ChunkMesh, IsomeshPlugin, NeedsRemesh, VoxelChunk, VoxelVolume};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::Sphere;

fn main() -> Result<(), isomesh::Error> {
    let mut app = App::new();
    // AssetPlugin must come first -- IsomeshPlugin asserts it and names the fix.
    app.add_plugins(DefaultPlugins).add_plugins(IsomeshPlugin);

    // 8 samples per chunk axis, 0.25 world units per cell, origin at zero.
    let layout = ChunkLayout::new(8, 0.25, [0.0; 3])?;
    let volume = app
        .world_mut()
        .spawn(VoxelVolume::new(layout, Sphere::<f32>::canonical()))
        .id();

    for x in -2..2 {
        app.world_mut().spawn((
            VoxelChunk { id: ChunkId::new([x, 0, 0]), volume },
            NeedsRemesh,
        ));
    }

    app.add_systems(Update, attach);
    app.run();
    Ok(())
}

// The plugin stops at a `Handle<Mesh>` on purpose. Attaching `Mesh3d` is the
// application's line to write, and that boundary is what keeps a CPU-only
// consumer from ever compiling the renderer.
fn attach(meshes: Query<(Entity, &ChunkMesh), Without<Mesh3d>>, mut commands: Commands) {
    for (entity, chunk_mesh) in &meshes {
        commands.entity(entity).insert(Mesh3d(chunk_mesh.0.clone()));
    }
}
```

Extraction runs on the `AsyncComputeTaskPool`, never in a system; what the frame budget bounds is turning finished extractions into assets. Editing a chunk mid-extraction re-queues it rather than swallowing the edit — the headless behavioral tests state these as claims (`spawning_the_work_does_not_do_the_work`, `an_edit_during_extraction_is_requeued_rather_than_swallowed`).

## Two layers — use either

**`MeshBuilder`** is a `MeshSink` whose buffers are exactly the arrays a Bevy `Mesh` wants, so `into_mesh()` hands the `Vec`s over by move — no copy on the way into the asset. This one runs, and does, on every `cargo test`:

```rust
use bevy_isomesh::MeshBuilder;
use isomesh::fields::Sphere;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::RuntimeShape3;

let mut builder = MeshBuilder::new();
let mut mc = MarchingCubes::<f32>::new();
let shape = RuntimeShape3::new([33; 3]).expect("valid shape");
mc.extract(&Sphere::<f32>::canonical(), &shape, [-2.0; 3], 0.125, &mut builder).expect("extraction");

let mesh = builder.into_mesh();
assert!(mesh.count_vertices() > 0);
```

**`IsomeshPlugin`** is the chunked, frame-budgeted layer above it, shown in the quickstart. `to_bevy_mesh` is the copying alternative for an `isomesh::MeshBuffer` you are reusing across chunks.

## What the plugin exposes

The contract, stated rather than implied:

- **Components** — `VoxelVolume` (a `ChunkLayout` + an `Arc<dyn VolumeField>` + an `Extractor`, default `MarchingCubes`), `VoxelChunk { id, volume }`, `NeedsRemesh` (marker; insert it to request work), `ChunkMesh(pub Handle<Mesh>)` (the result; yours to render).
- **Resources** — `MeshBudget` (default: 4 ms of asset-application per frame, at most 2× the task-pool parallelism in flight), `MeshStats` (spawned / applied / in-flight / waiting counters).
- **Systems** — `(spawn_meshing_tasks, apply_finished_meshes).chain()` in `Update`, both private. There is no public `SystemSet` yet; if you need to order against the plugin, that is a missing feature to request, not an omission to work around.
- **Traits** — `VolumeField`: any `isomesh::Sdf<Scalar = f32> + Send + Sync + 'static`, blanket-implemented.
- **Cargo features** — none.

`IsomeshPlugin` auto-adds `bevy_mesh::MeshPlugin` if absent, and asserts `AssetPlugin` was added first with a message that names the fix.

## Choosing an extractor for a chunked world

`Extractor::chunk_seams()` reports what was measured, per variant: `MarchingCubes` and `Subgrid` tile across chunk seams with **zero** boundary edges; `SurfaceNets` and `DualContouring` are **`Gapped`** — up to 5 open edges measured on a single seam, a structural property of one-vertex-per-cell placement, not a bug a future release fixes. For a chunked world, use `MarchingCubes` or `Subgrid { samples }`; the dual methods are for single-volume extraction where sharp features matter more than tiling.

## Why a separate workspace

This directory is deliberately excluded from the repository's root workspace. In a shared workspace, Cargo's feature unification would hand `glam` the `std`, `serde`, `bytemuck` and `encase` features Bevy asks for, and `cargo test` at the root would stop testing what a consumer of `isomesh` actually gets. The cost is running the two workspaces separately; the benefit is a root lockfile that proves the core crate's `no_std + libm` claim on every push.

Dependencies are Bevy's leaf crates (`bevy_asset`, `bevy_mesh`, `bevy_app`, `bevy_ecs`, `bevy_tasks`), never the `bevy` umbrella, so a CPU-only consumer never compiles the renderer.

Every Bevy example in the project also lives here, and CI builds them all on every commit. That is not tidiness — it is the one thing that stops them rotting. `block-mesh`'s examples are pinned to Bevy 0.13 and `fast-surface-nets`' to Bevy 0.7, both because they lived somewhere nothing ever built.

## Examples

31 of them, each a measured experiment with its claim in the header. A taste, with what each one proves:

| Example | The claim it makes | Evidence |
|---|---|---|
| `game_showcase` | caves and arches a heightfield cannot represent, flown through at 60 fps | [screenshot](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/screenshots/e210-showcase-hero.png) |
| `game_walk` | 495 chunk-seam crossings walked, 0 holes — the ball stands on the triangles, not the field | [screenshot](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/screenshots/e203-game-walk.png) |
| `game_dig` | carving tunnels re-meshes only the chunks the brush touched | [screenshot](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/screenshots/e202-game-dig.png) |
| `game_destruction` | the debris *is* the removed geometry — one boolean, two meshes | [screenshot](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/screenshots/e204-destruction-wall.png) |
| `game_lod_flyover` | an LOD ladder with 0 open seam edges, and 71–111 with transitions off | [screenshot](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/screenshots/e205-lod-flyover-stitched.png) |
| `game_paint` | graffiti survives the wall being blown open, drift exactly zero | [screenshot](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/screenshots/e208-game-paint.png) |
| `game_budget` | every chunk dirtied at once, drained under a frame budget that holds | [screenshot](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/screenshots/e206-budget.png) |
| `gpu_compute_mc` | GPU Marching Cubes agrees with the CPU vertex-for-vertex, worst case one ULP | [screenshot](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/screenshots/e301-gpu-compute-mc.png) |

```bash
cd bevy_isomesh
cargo run --example game_showcase --release
```

**Always `--release`.** A debug build meshes 37–62× slower — both ends measured in this repository ([FINDINGS](https://github.com/ladvien/isomesh/blob/main/FINDINGS.md) M-152) — and will convince you something is wrong with the algorithm.

The full inventory is [`bevy_isomesh/examples/`](https://github.com/ladvien/isomesh/tree/main/bevy_isomesh/examples); the measured figures live on the [gameplay](https://github.com/ladvien/isomesh/blob/main/docs/demos/gameplay.md), [algorithms](https://github.com/ladvien/isomesh/blob/main/docs/demos/algorithms.md) and [correctness](https://github.com/ladvien/isomesh/blob/main/docs/demos/correctness.md) pages.

## Requirements

Rust 1.95 (Bevy 0.19's own floor), edition 2024.

## License

MIT OR Apache-2.0, at your option.
