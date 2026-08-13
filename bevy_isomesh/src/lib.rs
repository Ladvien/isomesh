//! Bevy integration for [`isomesh`].
//!
//! `isomesh` itself knows nothing about Bevy — its public API is `[f32; 3]` and
//! `[u32; 3]` and it has one dependency. This crate is the whole of the Bevy
//! knowledge, and it lives in its own workspace so that Cargo's feature
//! unification cannot leak Bevy's feature choices into the core crate's
//! lockfile.
//!
//! # Why the examples live here
//!
//! Every Bevy example in this project is under `bevy_isomesh/examples/`, and CI
//! builds them on every commit. That is not organisational tidiness — it is the
//! one thing that stops them rotting. `block-mesh`'s examples are pinned to
//! Bevy 0.13 and `fast-surface-nets`' to Bevy 0.7, both because they lived in a
//! separate crate that nothing ever compiled.
//!
//! # Versions
//!
//! Bevy 0.19, which pins `wgpu` 29.0.3 and `glam` 0.32. Those move together or
//! not at all: Cargo will happily resolve two majors of `wgpu` side by side with
//! no error, and the failure surfaces much later as
//! `expected TextureFormat, found a different TextureFormat`.

/// Writing an extracted mesh into a Bevy [`Mesh`](bevy_mesh::Mesh).
pub mod mesh;
/// The plugin: chunked, asynchronous meshing under a frame budget.
pub mod plugin;

pub use mesh::{MeshBuilder, to_bevy_mesh};
pub use plugin::{
    ChunkMesh, Extractor, IsomeshPlugin, MeshBudget, MeshStats, NeedsRemesh, VolumeField,
    VoxelChunk, VoxelVolume,
};
