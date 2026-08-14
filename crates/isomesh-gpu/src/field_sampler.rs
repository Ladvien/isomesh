//! Evaluating the field on the GPU, so nothing is uploaded.
//!
//! # What this is worth, and what it costs
//!
//! After GPU-010a and GPU-012 the upload is **86% of the extraction path** —
//! 7.04 ms of 8.15 at 129³ — and what is left of it is `field.sample()` on the
//! CPU plus an irreducible memcpy. Producing the samples here removes both: they
//! are written where they are read.
//!
//! **The cost is the sample-identity that made GPU-005's comparison mean
//! something.** With [`FieldBuffer::sampled`](crate::FieldBuffer::sampled) both
//! sides read *the same `f32`s*, because the CPU produced them — which is why
//! M-142 could report vertex positions agreeing to one ULP and attribute the
//! residue to a fused multiply-add. Evaluate independently on each side and the
//! agreement becomes a property of each platform's arithmetic. See
//! [`GpuField`] for the measured per-field deviation.
//!
//! The CPU path keeps its guarantee: `libm` is bit-reproducible across machines
//! (M-31). **The GPU path acquires no such guarantee**, and that is a property
//! of the approach rather than of this implementation.
//!
//! # This is the cheapest version on purpose
//!
//! Four reference fields selected by a uniform, enough to answer whether the
//! 86% is recoverable and to measure the drift. Interpreting an edit log, or
//! composing a consumer's own WGSL, is GPU-011b — and the deviation numbers
//! here are what that design has to be built around.

use crate::{Composer, FieldBuffer, GridParams, Result};

/// Threads per workgroup, matching `@workgroup_size(64)` in `field.wgsl`.
const WORKGROUP: u32 = 64;

/// A reference field this crate can evaluate on the GPU.
///
/// # Measured deviation from `isomesh`'s CPU evaluation (M-154)
///
/// Same grid, same expressions in the same order, at `h = 0.125` from origin
/// `-2.0` — both powers of two, so the sample *positions* are bit-identical and
/// everything below is the field arithmetic alone.
///
/// | field | bit-exact | worst absolute |
/// |---|---:|---:|
/// | `sphere` | 26,009 / 35,937 | 2.38e-7 |
/// | `torus` | 25,449 / 35,937 | 4.77e-7 |
/// | `box_exact` | 30,637 / 35,937 | 1.19e-7 |
/// | `gyroid` | 3,873 / 35,937 | 6.56e-7 |
///
/// **Nothing agrees bit-for-bit, including the fields with no transcendentals**
/// — a GPU may contract `x*x + y*y + z*z` into fused multiply-adds, rounding
/// once where `libm` rounds twice. `box_exact` comes closest because its
/// expression is `abs`/`max`/`min`, exact in IEEE, over a `sqrt` whose argument
/// is identically zero inside the box.
///
/// Every deviation is under `7e-7`, about `5e-6` of a cell at this spacing —
/// four orders below anything that could move a crossing. **The shape of the
/// result is what GPU-011b has to design around**: the drift is a property of
/// the *expression*, not of the GPU, so a field built from exact operations
/// crosses unchanged and one built from products does not.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GpuField {
    /// Centre origin, radius 1.
    Sphere,
    /// Centre origin, major 1, minor 0.3.
    Torus,
    /// Centre origin, half-extents `[1, 1, 1]`.
    BoxExact,
    /// Scale 1, iso 0. The only one with `sin`/`cos`.
    Gyroid,
}

impl GpuField {
    /// Every field, in shader-id order.
    pub const ALL: [Self; 4] = [Self::Sphere, Self::Torus, Self::BoxExact, Self::Gyroid];

    /// The name `isomesh`'s `ReferenceField` uses.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Sphere => "sphere",
            Self::Torus => "torus",
            Self::BoxExact => "box_exact",
            Self::Gyroid => "gyroid",
        }
    }

    /// The `id` the shader switches on.
    #[must_use]
    pub const fn id(self) -> u32 {
        match self {
            Self::Sphere => 0,
            Self::Torus => 1,
            Self::BoxExact => 2,
            Self::Gyroid => 3,
        }
    }

    /// The same field, evaluated by `isomesh` on the CPU.
    ///
    /// The reference the GPU is compared against. It calls the crate's own
    /// canonical constructors rather than re-stating the constants, so the two
    /// cannot drift apart in the comparison that is supposed to detect drift.
    #[must_use]
    pub fn sample_on_cpu(self, p: [f32; 3]) -> f32 {
        use isomesh::Sdf;
        use isomesh::fields::{BoxExact, Gyroid, Sphere, Torus};
        match self {
            Self::Sphere => Sphere::<f32>::canonical().sample(p),
            Self::Torus => Torus::<f32>::canonical().sample(p),
            Self::BoxExact => BoxExact::<f32>::canonical().sample(p),
            Self::Gyroid => Gyroid::<f32>::canonical().sample(p),
        }
    }
}

/// A compiled pipeline that fills a [`FieldBuffer`] on the GPU.
#[derive(Debug)]
pub struct FieldSampler {
    pipeline: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
}

impl FieldSampler {
    /// Compile the kernel.
    ///
    /// # Errors
    ///
    /// Propagates composition failure, which would mean this crate's own
    /// shaders are inconsistent.
    pub fn new(device: &wgpu::Device) -> Result<Self> {
        let source = Composer::with_builtins().compose("field", &[])?;
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("isomesh field sampler"),
            source: wgpu::ShaderSource::Wgsl(source.into()),
        });

        let entry = |binding: u32, ty: wgpu::BufferBindingType| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty,
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("isomesh field sampler bindings"),
            entries: &[
                entry(0, wgpu::BufferBindingType::Uniform),
                entry(1, wgpu::BufferBindingType::Uniform),
                entry(2, wgpu::BufferBindingType::Storage { read_only: false }),
            ],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("isomesh field sampler layout"),
            bind_group_layouts: &[Some(&layout)],
            immediate_size: 0,
        });

        Ok(Self {
            pipeline: device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("isomesh field sampler"),
                layout: Some(&pipeline_layout),
                module: &module,
                entry_point: Some("sample_field"),
                compilation_options: wgpu::PipelineCompilationOptions::default(),
                cache: None,
            }),
            layout,
        })
    }

    /// Allocate a [`FieldBuffer`] and fill it on the GPU.
    ///
    /// Nothing is uploaded but two uniforms — the grid and the field id, 48
    /// bytes between them — where
    /// [`FieldBuffer::sampled`](crate::FieldBuffer::sampled) copies four bytes
    /// per sample.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooLarge`](crate::Error::GridTooLarge) if the sample count
    /// needs more workgroups than the adapter will dispatch.
    pub fn sample(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        params: GridParams,
        field: GpuField,
    ) -> Result<FieldBuffer> {
        let samples = params.sample_count();
        let groups = samples.div_ceil(u64::from(WORKGROUP));
        if groups > u64::from(device.limits().max_compute_workgroups_per_dimension) {
            return Err(crate::Error::GridTooLarge {
                samples: params.samples(),
            });
        }

        let grid = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("isomesh grid params"),
            size: GridParams::UNIFORM_SIZE,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&grid, 0, &params.to_std140());

        let select = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("isomesh field select"),
            size: 16,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut bytes = [0u8; 16];
        bytes[..4].copy_from_slice(&field.id().to_le_bytes());
        queue.write_buffer(&select, 0, &bytes);

        let out = FieldBuffer::new(device, params);

        fn entry(binding: u32, buffer: &wgpu::Buffer) -> wgpu::BindGroupEntry<'_> {
            wgpu::BindGroupEntry {
                binding,
                resource: buffer.as_entire_binding(),
            }
        }
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("isomesh field sampler"),
            layout: &self.layout,
            entries: &[entry(0, &grid), entry(1, &select), entry(2, out.buffer())],
        });

        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("isomesh field sampler"),
        });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("isomesh field sampler"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            // Bounded against the adapter's limit above.
            pass.dispatch_workgroups(groups as u32, 1, 1);
        }
        queue.submit(Some(encoder.finish()));

        Ok(out)
    }
}

#[cfg(test)]
mod tests;
