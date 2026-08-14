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
//! # The general mechanism is an interpreted edit log, and here is the argument
//!
//! GPU-011a proved the win with four fields hard-coded in WGSL, which is not
//! shippable — a consumer's own field cannot be one of four. GPU-011b had two
//! candidates and chose **interpreting an edit log**
//! ([`sample_stack`](FieldSampler::sample_stack)) over **composing the
//! consumer's WGSL** through [`Composer`](crate::Composer). Four reasons, in
//! the order they mattered:
//!
//! 1. **It is what this crate's world already is.** A scene here is a base
//!    field plus an ordered op list — `BrushStack`, and `paint::Edit` after
//!    E-208. E-207's undo is a re-fold of that log. Consuming the same log puts
//!    the GPU path on the architecture that exists rather than beside it.
//! 2. **It survives editing, and the alternative does not.** A carve pushes an
//!    op and the GPU re-reads a buffer. Composing WGSL would **recompile a
//!    shader per edit**, which is tolerable for a static CAD model and fatal
//!    for `game_dig` or `game_editor` — the demos where a GPU path would
//!    actually be used.
//! 3. **It is exactly testable.** `brush::apply` and `brush::smooth_min` define
//!    the fold, so the interpreter is asserted against them sample-for-sample.
//!    Composed WGSL can only be checked against a *second* statement of the
//!    field, which is the drift rule 5 exists to prevent — testable, but only
//!    after the duplication is already there.
//! 4. **M-154 says the drift is in the expression, not the GPU.** A field of
//!    `abs`/`max`/`min` crosses within 1 ULP; one of sums-of-products does not.
//!    An interpreter evaluating the same primitives has the same drift as
//!    compiled WGSL doing the same maths, so option (b) buys no accuracy.
//!
//! **What it costs, stated plainly:** the log is bounded by the primitives the
//! shader knows — sphere, box, capsule, and the three brush ops. A consumer
//! with an arbitrary analytic SDF is **not** served, and that is the CAD half
//! of this crate's audience. Composing consumer WGSL remains the right answer
//! for them and is additive rather than contradictory: a second *feature*, not
//! a second path to this one.

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

/// One brush shape, matching `isomesh`'s three.
///
/// Implements [`Sdf`](isomesh::Sdf) by delegating to `isomesh`'s own types, so
/// the CPU side of every comparison is the real definition rather than a second
/// copy of it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GpuShape {
    /// A sphere.
    Sphere {
        /// Centre.
        center: [f32; 3],
        /// Radius.
        radius: f32,
    },
    /// An axis-aligned box, exact distance.
    BoxExact {
        /// Centre.
        center: [f32; 3],
        /// Half-extents.
        half_extents: [f32; 3],
    },
    /// The set of points within `radius` of a segment.
    Capsule {
        /// One end.
        a: [f32; 3],
        /// The other.
        b: [f32; 3],
        /// Distance from the segment that counts as inside.
        radius: f32,
    },
}

impl isomesh::Sdf for GpuShape {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        use isomesh::brush::Capsule;
        use isomesh::fields::{BoxExact, Sphere};
        match *self {
            Self::Sphere { center, radius } => Sphere { center, radius }.sample(p),
            Self::BoxExact {
                center,
                half_extents,
            } => BoxExact {
                center,
                half_extents,
            }
            .sample(p),
            Self::Capsule { a, b, radius } => Capsule { a, b, radius }.sample(p),
        }
    }
}

/// What a brush does to the field it is applied to, matching `isomesh`'s
/// [`BrushOp`](isomesh::brush::BrushOp).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GpuOp {
    /// Union: `min(field, shape)`.
    Add,
    /// Difference: `max(field, -shape)`.
    Subtract,
    /// Union with a rounded join of width `k`.
    SmoothAdd {
        /// Join width, in world units.
        k: f32,
    },
}

/// One entry in the edit log.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GpuBrush {
    /// The shape being applied.
    pub shape: GpuShape,
    /// What to do with it.
    pub op: GpuOp,
}

impl GpuBrush {
    /// Bytes this occupies in the shader's brush array.
    pub const STRIDE: u64 = 64;

    /// The same brush as `isomesh` sees it.
    #[must_use]
    pub fn to_cpu(self) -> isomesh::brush::Brush<GpuShape> {
        use isomesh::brush::{Brush, BrushOp};
        Brush {
            shape: self.shape,
            op: match self.op {
                GpuOp::Add => BrushOp::Add,
                GpuOp::Subtract => BrushOp::Subtract,
                GpuOp::SmoothAdd { k } => BrushOp::SmoothAdd { k: f64::from(k) },
            },
        }
    }

    /// The 64 bytes the shader reads.
    ///
    /// Four `vec4`s, so std140 and std430 agree. The slot meanings are written
    /// out in `field.wgsl` next to the matching struct — a disagreement between
    /// the two produces a brush of the right kind in the wrong place, which
    /// looks like a meshing bug.
    #[must_use]
    pub fn to_std140(self) -> [u8; Self::STRIDE as usize] {
        let (kind, a, a_w, b, b_w) = match self.shape {
            GpuShape::Sphere { center, radius } => (0u32, center, radius, [0.0f32; 3], 0.0f32),
            GpuShape::BoxExact {
                center,
                half_extents,
            } => (1, center, 0.0, half_extents, 0.0),
            GpuShape::Capsule { a, b, radius } => (2, a, radius, b, 0.0),
        };
        let (op, k) = match self.op {
            GpuOp::Add => (0u32, 0.0f32),
            GpuOp::Subtract => (1, 0.0),
            GpuOp::SmoothAdd { k } => (2, k),
        };

        let mut out = [0u8; Self::STRIDE as usize];
        let words: [[u8; 4]; 16] = [
            kind.to_le_bytes(),
            op.to_le_bytes(),
            0u32.to_le_bytes(),
            0u32.to_le_bytes(),
            a[0].to_le_bytes(),
            a[1].to_le_bytes(),
            a[2].to_le_bytes(),
            a_w.to_le_bytes(),
            b[0].to_le_bytes(),
            b[1].to_le_bytes(),
            b[2].to_le_bytes(),
            b_w.to_le_bytes(),
            k.to_le_bytes(),
            0f32.to_le_bytes(),
            0f32.to_le_bytes(),
            0f32.to_le_bytes(),
        ];
        for (slot, word) in out.chunks_exact_mut(4).zip(words) {
            slot.copy_from_slice(&word);
        }
        out
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
                entry(3, wgpu::BufferBindingType::Storage { read_only: true }),
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
        self.sample_stack(device, queue, params, field, &[])
    }

    /// Allocate a [`FieldBuffer`] and fill it from a base field **and an edit
    /// log**, on the GPU.
    ///
    /// The general mechanism GPU-011b chose, and the argument for choosing it
    /// over composing a consumer's WGSL is on [`FieldSampler`].
    ///
    /// `brushes` is folded first to last over the base, exactly as
    /// [`BrushStack`](isomesh::brush::BrushStack) does. **That order is part of
    /// the value**: a log mixing adds and subtracts does not commute, and one
    /// containing a smooth add is not even associative (M-36..M-38), so a
    /// shader that reordered for parallelism would compute a different solid.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooLarge`](crate::Error::GridTooLarge) if the sample count
    /// needs more workgroups than the adapter will dispatch.
    pub fn sample_stack(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        params: GridParams,
        field: GpuField,
        brushes: &[GpuBrush],
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
        let count = u32::try_from(brushes.len()).map_err(|_| crate::Error::GridTooLarge {
            samples: params.samples(),
        })?;
        bytes[4..8].copy_from_slice(&count.to_le_bytes());
        queue.write_buffer(&select, 0, &bytes);

        // At least one element: wgpu rejects a zero-sized binding, and an empty
        // log is the ordinary case rather than an error.
        let log = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("isomesh brush log"),
            size: GpuBrush::STRIDE * brushes.len().max(1) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        if !brushes.is_empty() {
            let mut packed = Vec::with_capacity(brushes.len() * GpuBrush::STRIDE as usize);
            for brush in brushes {
                packed.extend_from_slice(&brush.to_std140());
            }
            queue.write_buffer(&log, 0, &packed);
        }

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
            entries: &[
                entry(0, &grid),
                entry(1, &select),
                entry(2, out.buffer()),
                entry(3, &log),
            ],
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
