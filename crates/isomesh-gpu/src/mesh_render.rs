//! A mesh-shader pipeline that draws the compute output where it already is.
//!
//! # What this is for
//!
//! M-145 measured the whole GPU route at 15 ms on a 129³ sphere, of which the
//! *extraction* was 0.13 ms. Everything else was moving data. This pipeline
//! consumes the position and normal buffers `marching_cubes.wgsl` wrote, in
//! place, so **drawing** builds no vertex or index buffer and reads nothing
//! back.
//!
//! # How much that is worth, stated rather than implied
//!
//! An earlier version of this page said it removed "the largest remaining
//! piece", and that was wrong (M-149). GPU-010a has since moved the prefix sum
//! onto the GPU, so the breakdown at 129³ is now:
//!
//! | | ms | share of the GPU path |
//! |---|---:|---:|
//! | upload | 7.04 | **86%** |
//! | **geometry read-back — what this removes** | **0.63** | **7.7%** |
//! | prefix scan + 4-byte total | 0.37 | 4.5% |
//! | count + emit | 0.11 | 1.3% |
//!
//! So this is worth about **7.7%** at 129³ — a larger share than before only
//! because the path around it got shorter — 1.56× at GPU-010a, another 1.18×
//! at GPU-012. **The dominant cost is now the
//! upload, at 87%**, and that is field evaluation on the CPU rather than
//! anything a renderer can fix: `FieldBuffer::sampled` evaluates the SDF host-
//! side and copies the samples over. Evaluating the field *in* the shader is
//! the next real lever, and it is not ticketed yet.
//!
//! # It refuses rather than falls back
//!
//! [`new`](MeshShaderRenderer::new) returns [`Error::MeshShadersUnavailable`]
//! on a device without `EXPERIMENTAL_MESH_SHADER`. It does **not** quietly
//! build a vertex-buffer pipeline instead — that would be a second execution
//! path for one feature, and a caller told "drawing" while a different pipeline
//! ran has been misinformed about the only thing they asked.
//!
//! A caller that wants to degrade gracefully asks [`is_supported`](
//! MeshShaderRenderer::is_supported) first and chooses, visibly, in its own
//! code. That is one path selected by a measurement, which is a different thing
//! from two paths selected silently.
//!
//! # Availability is narrower than "the adapter supports it"
//!
//! Three gates, and all three have to pass (M-146, M-147, V-23):
//!
//! - The **adapter** advertises `EXPERIMENTAL_MESH_SHADER`.
//! - The **device** was created with it. That needs an `ExperimentalFeatures`
//!   token whose constructor is `unsafe`, so this crate cannot make one — but
//!   Bevy does, and its default `Functionality` priority requests every feature
//!   the adapter has. A device from [`headless::Gpu`](crate::headless::Gpu)
//!   therefore does **not** qualify, which is what GPU-009 records.
//! - The backend is **Vulkan**. WGSL mesh shaders go through naga, and wgpu's
//!   own source says naga supports them on Vulkan only; Metal and DX12 need
//!   pre-compiled passthrough shaders, which this crate does not produce.

use crate::{Composer, Error, Result};

/// Triangles one mesh workgroup emits. Must match `BATCH` in the shader.
const BATCH: u32 = 32;

/// A compiled mesh-shader pipeline for isomesh's triangle soup.
#[derive(Debug)]
pub struct MeshShaderRenderer {
    pipeline: wgpu::RenderPipeline,
    layout: wgpu::BindGroupLayout,
}

impl MeshShaderRenderer {
    /// Whether `device` can run a mesh-shader pipeline at all.
    ///
    /// Checks the **device**'s features rather than the adapter's, because a
    /// device only has what was requested at creation — an adapter that
    /// advertises mesh shaders says nothing about a device that did not ask
    /// for them.
    ///
    /// Also checks the **backend** (V-23): this pipeline's shader is WGSL,
    /// naga translates the mesh and task stages to SPIR-V only, and its MSL
    /// writer hits `unimplemented!()` on them (naga 29.0.4,
    /// `back/msl/writer.rs:6937`). Metal genuinely *advertises*
    /// `EXPERIMENTAL_MESH_SHADER` — the feature is real for callers shipping
    /// pre-compiled passthrough MSL, which this crate does not — so the
    /// feature bit alone would send [`new`](Self::new) into an abort inside
    /// pipeline creation instead of an error value. The positive gate
    /// (`== Vulkan` rather than `!= Metal`) is deliberate: naga's HLSL
    /// backend lacks the stages too, so allow-listing the one backend that
    /// compiles them is the single path.
    #[must_use]
    pub fn is_supported(device: &wgpu::Device) -> bool {
        device
            .features()
            .contains(wgpu::Features::EXPERIMENTAL_MESH_SHADER)
            && device.adapter_info().backend == wgpu::Backend::Vulkan
    }

    /// Compile the pipeline for a colour target of `format`.
    ///
    /// # Errors
    ///
    /// [`Error::MeshShadersUnavailable`] if the device lacks the feature or
    /// the backend cannot compile WGSL mesh stages (only Vulkan can — V-23).
    /// This is the whole "never panics on an unsupported adapter" requirement:
    /// the answer is an error value, checked at the call site, not an abort
    /// inside a driver.
    /// `samples` is the view's MSAA sample count and **must** match the pass
    /// this pipeline is used in. It is a parameter rather than a default
    /// because getting it wrong is not a quality difference — wgpu refuses the
    /// draw outright with *"the RenderPass uses textures with sample count 4 but
    /// the RenderPipeline uses attachments with format 1"*, which is how this
    /// signature came to have it. Bevy defaults to 4.
    pub fn new(
        device: &wgpu::Device,
        format: wgpu::TextureFormat,
        depth: Option<wgpu::DepthStencilState>,
        samples: u32,
    ) -> Result<Self> {
        if !Self::is_supported(device) {
            return Err(Error::MeshShadersUnavailable);
        }

        let source = Composer::with_builtins().compose("mesh_render", &[])?;
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("isomesh mesh render"),
            source: wgpu::ShaderSource::Wgsl(source.into()),
        });

        let uniform = |binding: u32| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::MESH,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Uniform,
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let storage = |binding: u32| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::MESH,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Storage { read_only: true },
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("isomesh mesh render bindings"),
            entries: &[uniform(0), uniform(1), storage(2), storage(3)],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("isomesh mesh render layout"),
            bind_group_layouts: &[Some(&layout)],
            immediate_size: 0,
        });

        let pipeline = device.create_mesh_pipeline(&wgpu::MeshPipelineDescriptor {
            label: Some("isomesh mesh render"),
            layout: Some(&pipeline_layout),
            // No task stage. The dispatch count is known on the CPU -- it is
            // the triangle count divided by the batch -- so there is nothing
            // for a task shader to decide, and adding one to look modern would
            // be a stage that computes a constant.
            task: None,
            mesh: wgpu::MeshState {
                module: &module,
                entry_point: Some("draw_mesh"),
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            },
            primitive: wgpu::PrimitiveState {
                // Counter-clockwise from outside, which is this project's
                // winding convention, stated once in `isomesh`'s crate docs.
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: Some(wgpu::Face::Back),
                ..wgpu::PrimitiveState::default()
            },
            depth_stencil: depth,
            multisample: wgpu::MultisampleState {
                count: samples,
                ..wgpu::MultisampleState::default()
            },
            fragment: Some(wgpu::FragmentState {
                module: &module,
                entry_point: Some("shade"),
                compilation_options: wgpu::PipelineCompilationOptions::default(),
                targets: &[Some(wgpu::ColorTargetState {
                    format,
                    blend: None,
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            multiview: None,
            cache: None,
        });

        Ok(Self { pipeline, layout })
    }

    /// The bind group layout, for building a group against.
    #[must_use]
    pub const fn bind_group_layout(&self) -> &wgpu::BindGroupLayout {
        &self.layout
    }

    /// Bind the camera, the draw parameters and the two geometry buffers.
    ///
    /// `positions` and `normals` are the buffers
    /// [`MarchingCubesGpu`](crate::MarchingCubesGpu) wrote. They are passed as
    /// `&wgpu::Buffer` rather than a `GpuMesh`, because the entire point is
    /// that they were never read back into one.
    #[must_use]
    pub fn bind_group(
        &self,
        device: &wgpu::Device,
        camera: &wgpu::Buffer,
        draw: &wgpu::Buffer,
        positions: &wgpu::Buffer,
        normals: &wgpu::Buffer,
    ) -> wgpu::BindGroup {
        fn entry(binding: u32, buffer: &wgpu::Buffer) -> wgpu::BindGroupEntry<'_> {
            wgpu::BindGroupEntry {
                binding,
                resource: buffer.as_entire_binding(),
            }
        }
        device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("isomesh mesh render"),
            layout: &self.layout,
            entries: &[
                entry(0, camera),
                entry(1, draw),
                entry(2, positions),
                entry(3, normals),
            ],
        })
    }

    /// Record a draw whose workgroup count the **GPU** supplies.
    ///
    /// The count and the triangle total were written by
    /// [`MarchingCubesGpu::extract_indirect`](crate::MarchingCubesGpu::extract_indirect)
    /// and never came home, so nothing on this path waits for the extraction to
    /// finish. Bind `draw_params` as binding 1 rather than a CPU-written
    /// uniform — the two are the same 16 bytes, one filled by a kernel.
    pub fn draw_indirect(
        &self,
        pass: &mut wgpu::RenderPass<'_>,
        bind_group: &wgpu::BindGroup,
        indirect: &wgpu::Buffer,
    ) {
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, bind_group, &[]);
        pass.draw_mesh_tasks_indirect(indirect, 0);
    }

    /// Record the draw. `triangles` is the soup's triangle count.
    pub fn draw(
        &self,
        pass: &mut wgpu::RenderPass<'_>,
        bind_group: &wgpu::BindGroup,
        triangles: u32,
    ) {
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, bind_group, &[]);
        pass.draw_mesh_tasks(triangles.div_ceil(BATCH), 1, 1);
    }

    /// The camera uniform's bytes: a column-major 4×4, as WGSL reads `mat4x4`.
    ///
    /// `[[f32; 4]; 4]` rather than a matrix type, because rule 1 keeps math
    /// libraries out of public signatures — a consumer on glam 0.32 and one on
    /// 0.33 both hand over the same array.
    #[must_use]
    pub fn camera_bytes(view_proj: [[f32; 4]; 4]) -> [u8; 64] {
        let mut out = [0u8; 64];
        for (column, chunk) in view_proj.iter().zip(out.chunks_exact_mut(16)) {
            for (value, slot) in column.iter().zip(chunk.chunks_exact_mut(4)) {
                slot.copy_from_slice(&value.to_le_bytes());
            }
        }
        out
    }

    /// The draw uniform's bytes. One `u32`, padded to a 16-byte binding.
    #[must_use]
    pub fn draw_bytes(triangles: u32) -> [u8; 16] {
        let mut out = [0u8; 16];
        out[..4].copy_from_slice(&triangles.to_le_bytes());
        out
    }
}

#[cfg(test)]
mod tests;
