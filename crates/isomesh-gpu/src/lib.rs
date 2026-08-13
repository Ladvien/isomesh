//! GPU isosurface extraction for [`isomesh`], on raw `wgpu`.
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
mod grid;
mod shader;

pub mod headless;

pub use buffers::{FieldBuffer, read_buffer};
pub use error::{Error, Result};
pub use grid::GridParams;
pub use shader::{Composer, FEATURES, GRID_WGSL};
