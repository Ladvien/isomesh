// The README is this crate's front page — `docs/bevy_plugins.md` rule 1, the
// header Bevy's own crates use. Its fences compile as doctests on every
// `cargo test`, which is why the quickstart cannot rot; everything the old
// crate-level docs said (the separate-workspace rationale, the version pins,
// why every example lives here) moved into it.
#![doc = include_str!("../README.md")]

/// Writing an extracted mesh into a Bevy [`Mesh`](bevy_mesh::Mesh).
pub mod mesh;
/// The plugin: chunked, asynchronous meshing under a frame budget.
pub mod plugin;

pub use mesh::{MeshBuilder, SoupError, from_bevy_mesh, to_bevy_mesh};
pub use plugin::{
    ChunkMesh, ChunkSeams, Extractor, IsomeshPlugin, MeshBudget, MeshStats, NeedsRemesh,
    VolumeField, VoxelChunk, VoxelVolume,
};
