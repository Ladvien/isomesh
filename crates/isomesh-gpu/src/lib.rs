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
//! # What is here
//!
//! A whole extraction pipeline, not a substrate waiting for one.
//!
//! | | |
//! |---|---|
//! | **Extraction** | [`MarchingCubesGpu`] — the case table is uploaded rather than transcribed, so there is one table in the repository and the GPU reads the same bytes the CPU does ([`case_table_bytes`]) |
//! | **Fields on the GPU** | [`GpuField`], [`GpuShape`], [`GpuOp`], [`GpuBrush`], [`FieldSampler`] — brushes folded device-side, so an edit does not cross the bus |
//! | **The plumbing extraction needs** | [`PrefixScan`] and [`DeferredScan`] for compaction, [`GridParams`] for std140 packing, [`FieldBuffer`], [`read_buffer`] |
//! | **Drawing without a readback** | [`MeshShaderRenderer`], [`IndirectGeometry`], and [`probe_mesh_shaders`] to ask first |
//! | **Distance transforms** | [`JumpFlood`] |
//! | **Shader sources** | [`MARCHING_CUBES_WGSL`], [`FIELD_WGSL`], [`GRID_WGSL`], [`SCAN_WGSL`], [`JUMP_FLOOD_WGSL`], [`MESH_RENDER_WGSL`], composed by [`Composer`] |
//! | **No device of your own** | [`headless`] |
//!
//! # Is it faster than the CPU? Above about 33³, yes — by 37× at 129³
//!
//! Sphere, warmed, median of three, RTX 3090 over Vulkan, against a
//! single-threaded CPU extraction. `docs/measurements/gpu_vs_cpu.csv`.
//!
//! | samples/axis | CPU | GPU, field evaluated on the GPU | |
//! |---|---|---|---|
//! | 17³ | **0.06 ms** | 0.22 ms | CPU ahead 3.7× |
//! | 33³ | 0.34 ms | **0.23 ms** | GPU ahead 1.5× |
//! | 65³ | 2.44 ms | **0.27 ms** | GPU ahead 9× |
//! | 129³ | 20.14 ms | **0.54 ms** | GPU ahead **37×** |
//!
//! **The shape is the finding.** That GPU column is nearly flat across a 420×
//! rise in cell count, because extraction was never the cost: `count + emit` is
//! **0.045 ms at 129³** and does not move with resolution. Below ~33³ a fixed
//! ~0.22 ms of setup is larger than the whole job, and the CPU wins.
//!
//! # Where you evaluate the field decides the rest
//!
//! Sample on the CPU and hand over a [`FieldBuffer`] and the **upload is 87% of
//! the path** — 8.37 ms at 129³, so 2.4× ahead instead of 37×. Worse, that
//! design does not take field evaluation off the CPU's budget; it adds a copy to
//! it, and field evaluation is 65–74% of the whole job on a noise field.
//!
//! Evaluate it in the shader instead — [`GpuField`], [`GpuShape`], [`GpuOp`],
//! [`GpuBrush`] — and the upload stops existing, because the samples are
//! produced where they are read.
//!
//! A base the shader has no name for — an arbitrary analytic SDF, which is half
//! this crate's audience — takes the third option:
//! [`FieldSampler::fold_into`] reads that base from a buffer the caller filled
//! and folds the edit log over it on the device. **One upload per grid instead
//! of none, and none per edit**, which is what makes it the right shape for an
//! editor: `bevy_isomesh`'s `game_dig` samples its terrain once per chunk and
//! thereafter moves 64 bytes per surviving brush and no samples at all.
//!
//! Three tickets took this path from 15.01 ms to 0.54 ms at 129³ and **none of
//! them made the extractor faster.** Every gain was data movement removed: a GPU
//! prefix scan so 8.4 MB of per-cell counts never come home (M-150), then
//! device-side field evaluation (M-155). What is left after that is the geometry
//! read-back, and [`MeshShaderRenderer`] removes even that by drawing straight
//! out of the compute output.
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

pub use buffers::{
    FieldBuffer, Readback, read_buffer, read_buffer_u32, read_bytes, read_bytes_many,
    read_bytes_many_deferred,
};
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
