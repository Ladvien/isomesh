# bevy_isomesh

> ⚠️ **Vibe Coded.** This crate was written by an AI agent working from a research corpus and a ticket
> queue. Every number it quotes is produced by a test or benchmark in this repository and every
> algorithm cites its source, but it has not been through human code review. Read the tests before
> trusting it with anything that matters.

**Bevy 0.19 integration for [`isomesh`](../crates/isomesh): signed distance field in, `Mesh` asset
out.** The core crates stay engine-agnostic — this is the one place Bevy types appear.

Two layers, use either:

- **`MeshBuilder` / `to_bevy_mesh`** — a `MeshSink` that writes extraction output straight into the
  `Vec`s a Bevy `Mesh` will own, so conversion is a move, not a copy.
- **`IsomeshPlugin`** — chunked, frame-budgeted re-meshing on the `AsyncComputeTaskPool`: mark a
  `VoxelVolume`'s chunks `NeedsRemesh`, get `ChunkMesh` components back without blocking a frame.

## Why a separate workspace

This directory is deliberately excluded from the repository's root workspace. In a shared workspace,
Cargo's feature unification would hand `glam` the `std`, `serde`, `bytemuck` and `encase` features
Bevy asks for, and `cargo test` at the root would stop testing what a consumer of `isomesh` actually
gets. The cost is running the two workspaces separately; the benefit is a root lockfile that proves
the core crate's `no_std + libm` claim on every push.

Dependencies are Bevy's leaf crates (`bevy_asset`, `bevy_mesh`, `bevy_app`, …), never the `bevy`
umbrella, so a CPU-only consumer never compiles the renderer.

## Examples

Every example in the project lives here, so CI compiles all of them (the root crate has none by
design). The catalog with per-example detail is
[`docs/2026-08-11-bevy-examples-catalog.md`](../docs/2026-08-11-bevy-examples-catalog.md); the
rendered demo pages are linked from the [root README](../README.md).

```bash
cd bevy_isomesh
cargo run --example marching_cubes_sphere --release
cargo run --example game_showcase --release        # the kitchen sink
cargo run --example gpu_vs_cpu --release           # bevy's device driving isomesh-gpu
```

**Always `--release`.** A debug build meshes 37–62× slower — both ends measured in this repository
(FINDINGS M-152; the resolution-sweep incident in `crates/isomesh/Cargo.toml`) — and will convince
you something is wrong with the algorithm.

## License

MIT OR Apache-2.0, matching the workspace.
