//! GPU isosurface extraction for [`isomesh`], on raw `wgpu`.
//!
//! ```no_run
//! use isomesh_gpu::{MarchingCubesGpu, headless};
//!
//! // Any wgpu device will do -- Bevy's via `RenderDevice::wgpu_device()`, a CAD
//! // tool's own, or this one. The API never names an engine type.
//! let gpu = headless::Gpu::new()?;
//! let mut mc = MarchingCubesGpu::new(gpu.device(), gpu.queue())?;
//! # let _ = &mut mc;
//! # Ok::<(), isomesh_gpu::Error>(())
//! ```
//!
//! # The one rule this crate exists to keep
//!
//! **Every public entry point takes `&wgpu::Device`, `&wgpu::Queue` or
//! `&mut wgpu::CommandEncoder`, and never an engine type.** A Bevy consumer
//! reaches the raw device through `RenderDevice::wgpu_device()` and hands it in;
//! a CAD tool passes its own; a test passes one from [`headless`]. There is no
//! second entry point that takes a `RenderDevice`, because that would make the
//! engine a dependency of the algorithm rather than a caller of it.
//!
//! That is also why [`headless`] exists in the library and not only in tests: a
//! CAD tool with no renderer at all still needs a device, and a GPU path that
//! can only be reached through a game engine has already leaked.
//!
//! # The wgpu version is a hard pin, and getting it wrong does not fail loudly
//!
//! `wgpu 29.0.3` exactly, matching Bevy 0.19. Cargo resolves two wgpu majors
//! **side by side with no resolution error** — verified in
//! `docs/research/2026-08-11-meshing-crate-architecture.md` §6, where adding
//! `wgpu 30` alongside `bevy 0.19` locks 317 packages containing both. The
//! failure arrives much later and reads as a type error about `TextureFormat`
//! not matching itself.
//!
//! | `isomesh-gpu` | `wgpu` | `bevy` |
//! |---|---|---|
//! | 0.0.x | 29.0.3 | 0.19 |
//!
//! # What is here now
//!
//! The substrate an extraction pipeline sits on, and it is complete rather than
//! sketched: a validated [`GridParams`] with the std140 packing a shader will
//! read, a [`FieldBuffer`] that puts a sampled `isomesh` field into GPU memory,
//! and [`read_buffer`] to get results back. Shaders, their composition and
//! Marching Cubes itself are GPU-002 through GPU-004; nothing here is a
//! placeholder waiting for them.
//!
//! ```no_run
//! use isomesh::fields::Sphere;
//! use isomesh_gpu::{FieldBuffer, GridParams, headless};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let gpu = headless::Gpu::new()?;
//! let grid = GridParams::new([33; 3], [-2.0; 3], 0.125)?;
//!
//! // Takes &Device and &Queue -- never an engine type.
//! let field = FieldBuffer::sampled(
//!     gpu.device(),
//!     gpu.queue(),
//!     grid,
//!     &Sphere::<f32>::canonical(),
//! )?;
//! assert_eq!(field.params().sample_count(), 33 * 33 * 33);
//! # Ok(())
//! # }
//! ```

mod block_on;
mod buffers;
mod error;
mod field_sampler;
mod grid;
mod marching_cubes;
mod mesh_render;
mod mesh_shader;
mod scan;
mod shader;

pub mod headless;
pub mod jump_flood;

pub use buffers::{FieldBuffer, read_buffer, read_buffer_u32, read_bytes, read_bytes_many};
pub use error::{Error, Result};
pub use field_sampler::{FieldSampler, GpuBrush, GpuField, GpuOp, GpuShape};
pub use grid::GridParams;
pub use jump_flood::JumpFlood;
pub use marching_cubes::{
    ExtractTimings, GpuGeometry, GpuMesh, IndirectGeometry, MarchingCubesGpu, case_table_bytes,
};
pub use mesh_render::MeshShaderRenderer;
pub use mesh_shader::{MeshShaderReport, probe_mesh_shaders};
pub use scan::{DeferredScan, PrefixScan, ScanOutput, cpu_prefix_sum};
pub use shader::{
    Composer, FEATURES, FIELD_WGSL, GRID_WGSL, JUMP_FLOOD_WGSL, MARCHING_CUBES_WGSL,
    MESH_RENDER_WGSL, SCAN_WGSL,
};

/// Compiles the README's example as a doctest, without putting the README into
/// these docs — the same pattern `isomesh` uses. The README's fence is
/// `rust,no_run`: it must compile everywhere, including a CI runner with no
/// adapter, and execute nowhere, because [`headless::Gpu::new`] refuses to fall
/// back to a software device.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeExample;
