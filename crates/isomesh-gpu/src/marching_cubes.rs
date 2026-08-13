//! Marching Cubes on a compute shader, and the mesh it reads back.
//!
//! # The case table is shared, not transcribed
//!
//! `CLAUDE.md` rule 5 forbids guessing a case table, and the usual way a GPU
//! port breaks it is by pasting a second copy of the 256 cases into WGSL. There
//! is no second copy here. [`case_table_bytes`] packs
//! [`isomesh::marching_cubes::table::CASES`] — which is itself *derived* by a
//! `const fn` rather than copied from a paper — and uploads it, so the shader
//! runs the same table the CPU does and the two cannot disagree.
//!
//! # Two passes, because the output has to be dense *and* ordered
//!
//! 1. `count_cells` classifies every cell and writes its triangle count.
//! 2. The CPU prefix-sums those counts and sizes the output from the total.
//! 3. `emit_cells` classifies again and writes at this cell's offset.
//!
//! Classifying twice is the cheaper mistake. An atomic bump allocator needs one
//! pass, but then the output order depends on which workgroup arrived first,
//! and a mesh whose vertex order changes between runs cannot be compared with
//! anything — including itself, which is what the determinism rule asks for.
//! Cell order is the CPU extractor's order, so the two line up by construction.
//!
//! # What differs from the CPU path, stated rather than discovered
//!
//! **Normals.** The CPU calls `Sdf::gradient`, which every reference field
//! overrides analytically. A shader has only the uploaded samples, so it uses
//! central differences at the cell size over a trilinear read of the grid —
//! the same quantity `isomesh`'s `CentralDifference { step: h }` computes.
//! M-65 measured that difference at **0.460° worst, 0.299° mean at 17³,
//! converging at `h²`**.
//!
//! **Vertices are not shared.** This emits a triangle soup: three vertices per
//! triangle, no cache. The CPU path keys a cache on the grid edge and emits one
//! vertex per crossing. So vertex *counts* differ by construction and the
//! surfaces are compared as geometry. Welding is `isomesh::weld`'s job and runs
//! on either.
//!
//! **Plain Marching Cubes.** The uploaded table is `CASES`, the all-separate
//! resolution. A-002's asymptotic decider is a per-cell run-time choice and is
//! not ported here.

use isomesh::marching_cubes::table::{self, MAX_TRIANGLES};

use crate::{Composer, Error, FieldBuffer, GridParams, Result, read_buffer, read_buffer_u32};

/// Words per case in the uploaded table: a header, then one per triangle.
const CASE_STRIDE: usize = 1 + MAX_TRIANGLES;

/// Threads per workgroup, matching `@workgroup_size(64)` in the shader.
const WORKGROUP: u32 = 64;

/// The 256-case table, packed for the shader.
///
/// 13 words per case, little-endian:
///
/// ```text
/// header    = count | (centroids << 8)
/// triangle  = a | (b << 8) | (c << 16)
/// ```
///
/// Corner codes are one byte each and `NO_EDGE` is `0xff`, which never appears
/// in the first `count` triangles — the only ones either side reads.
#[must_use]
pub fn case_table_bytes() -> Vec<u8> {
    let mut out = Vec::with_capacity(256 * CASE_STRIDE * 4);
    for entry in &table::CASES {
        let header = u32::from(entry.count) | (u32::from(entry.centroids) << 8);
        out.extend_from_slice(&header.to_le_bytes());
        for tri in &entry.triangles {
            let word = u32::from(tri[0]) | (u32::from(tri[1]) << 8) | (u32::from(tri[2]) << 16);
            out.extend_from_slice(&word.to_le_bytes());
        }
    }
    out
}

/// Where the time went inside one [`extract`](MarchingCubesGpu::extract).
///
/// Always filled in, rather than sitting behind a flag or a second `_timed`
/// entry point. Five `Instant::now()` calls against a dispatch is nothing, and
/// a timing path that has to be opted into is a timing path nobody has run —
/// which is how "the GPU is slower" gets reported without anyone knowing
/// *which part* of it is slower.
///
/// Wall-clock from the CPU's side. `poll(Wait)` inside the read-backs is what
/// makes these meaningful at all: without a wait the submit returns
/// immediately and every millisecond lands on whichever call happens to block
/// first. Timestamp queries would attribute GPU-side time more precisely and
/// need a device feature this crate does not request.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ExtractTimings {
    /// Classifying every cell and writing per-cell triangle counts.
    pub count_ms: f64,
    /// Reading those counts back — the barrier the two-pass design costs.
    pub counts_readback_ms: f64,
    /// The CPU-side exclusive prefix sum and the offsets upload.
    pub prefix_ms: f64,
    /// Classifying again and writing the triangles.
    pub emit_ms: f64,
    /// Reading positions and normals back.
    pub geometry_readback_ms: f64,
}

impl ExtractTimings {
    /// Everything above, summed.
    #[must_use]
    pub fn total_ms(&self) -> f64 {
        self.count_ms
            + self.counts_readback_ms
            + self.prefix_ms
            + self.emit_ms
            + self.geometry_readback_ms
    }

    /// The share spent moving data back to the CPU rather than computing.
    ///
    /// The number that decides whether a GPU path is worth it at all: a
    /// consumer rendering straight from GPU memory never pays it, and one that
    /// needs a collider always does.
    #[must_use]
    pub fn readback_share(&self) -> f64 {
        let total = self.total_ms();
        if total > 0.0 {
            (self.counts_readback_ms + self.geometry_readback_ms) / total
        } else {
            0.0
        }
    }
}

/// A triangle soup read back from the GPU.
///
/// Positions and normals are parallel and three per triangle, in cell order.
/// There is no index buffer because there are no shared vertices — see the
/// module docs.
/// # Not `PartialEq`, deliberately
///
/// [`timings`](Self::timings) is wall-clock and never repeats, so a derived
/// equality would make two runs of the same input compare **unequal** — the
/// exact opposite of what anyone reaching for `==` on a mesh wants. Comparing
/// geometry means comparing `positions` and `indices`; saying so at the call
/// site is clearer than a hand-written `PartialEq` that silently drops a field.
#[derive(Clone, Debug, Default)]
pub struct GpuMesh {
    /// One per vertex.
    pub positions: Vec<[f32; 3]>,
    /// One per vertex, parallel to `positions`. Unit length, or zero where the
    /// field's gradient vanishes.
    pub normals: Vec<[f32; 3]>,
    /// Where the time went producing this.
    pub timings: ExtractTimings,
}

impl GpuMesh {
    /// Triangles in the soup.
    #[must_use]
    pub fn triangle_count(&self) -> usize {
        self.positions.len() / 3
    }
}

/// A compiled Marching Cubes pipeline, and the table it dispatches against.
///
/// Built once and reused: shader compilation and the table upload are the
/// expensive parts, and neither depends on the field or the grid.
#[derive(Debug)]
pub struct MarchingCubesGpu {
    count: wgpu::ComputePipeline,
    emit: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
    cases: wgpu::Buffer,
}

impl MarchingCubesGpu {
    /// Compose the shader, compile both entry points, and upload the table.
    ///
    /// # Errors
    ///
    /// [`Error::ShaderModuleMissing`] and friends if composition fails, which
    /// would mean this crate's own shaders are inconsistent.
    pub fn new(device: &wgpu::Device, queue: &wgpu::Queue) -> Result<Self> {
        // Composed through the same registry GPU-003's sweep validates, so
        // the source compiled here is the source that was checked.
        let source = Composer::with_builtins().compose("marching_cubes", &[])?;

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("isomesh marching cubes"),
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
        let read = wgpu::BufferBindingType::Storage { read_only: true };
        let write = wgpu::BufferBindingType::Storage { read_only: false };
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("isomesh marching cubes bindings"),
            entries: &[
                entry(0, wgpu::BufferBindingType::Uniform),
                entry(1, read),
                entry(2, read),
                entry(3, write),
                entry(4, write),
                entry(5, write),
            ],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("isomesh marching cubes layout"),
            bind_group_layouts: &[Some(&layout)],
            // wgpu 29 renamed push constants to "immediates" and sizes them
            // rather than taking ranges. Nothing here uses them.
            immediate_size: 0,
        });

        let compile = |entry_point: &str| {
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(entry_point),
                layout: Some(&pipeline_layout),
                module: &module,
                entry_point: Some(entry_point),
                compilation_options: wgpu::PipelineCompilationOptions::default(),
                cache: None,
            })
        };

        let bytes = case_table_bytes();
        let cases = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("isomesh marching cubes case table"),
            size: bytes.len() as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&cases, 0, &bytes);

        Ok(Self {
            count: compile("count_cells"),
            emit: compile("emit_cells"),
            layout,
            cases,
        })
    }

    /// Extract `field` and read the triangles back.
    ///
    /// Blocks twice: once for the counts, once for the geometry. That is what
    /// the two-pass design costs and it is the price of a deterministic,
    /// densely packed output.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooLarge`] if the cell count exceeds what one dispatch can
    /// cover on this adapter, plus anything read-back can report.
    pub fn extract(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        field: &FieldBuffer,
    ) -> Result<GpuMesh> {
        let params = field.params();
        let cells = params.cell_count();
        let groups = cells.div_ceil(u64::from(WORKGROUP));
        if groups > u64::from(device.limits().max_compute_workgroups_per_dimension) {
            return Err(Error::GridTooLarge {
                samples: params.samples(),
            });
        }
        // Every buffer below is indexed by a u32 in the shader, so the cell
        // count has to fit one. Checked here rather than trusted, because the
        // symptom of not checking is a silently wrapped index writing over
        // another cell's triangles.
        let cell_words = u32::try_from(cells).map_err(|_| Error::GridTooLarge {
            samples: params.samples(),
        })?;

        let uniform = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("isomesh grid params"),
            size: GridParams::UNIFORM_SIZE,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&uniform, 0, &params.to_std140());

        let counts_bytes = u64::from(cell_words) * 4;
        let counts = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("isomesh cell triangle counts"),
            size: counts_bytes,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_DST
                | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        // Pass one needs the output bindings to exist but writes nothing to
        // them, so they start at the minimum wgpu will bind rather than at the
        // size the mesh will need -- which is not known yet.
        let placeholder = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("isomesh unused output"),
            size: 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        let mut timings = ExtractTimings::default();

        let started = std::time::Instant::now();
        self.dispatch(
            device,
            queue,
            &uniform,
            field,
            &counts,
            &placeholder,
            &placeholder,
            &self.count,
            groups,
        );
        timings.count_ms = started.elapsed().as_secs_f64() * 1000.0;

        let started = std::time::Instant::now();
        let per_cell = read_buffer_u32(device, queue, &counts, counts_bytes)?;
        timings.counts_readback_ms = started.elapsed().as_secs_f64() * 1000.0;
        let started = std::time::Instant::now();

        // Exclusive prefix sum on the CPU. A scan kernel is the obvious
        // optimisation and is not this ticket -- the readback is already the
        // dominant cost and a wrong scan is a class of bug this does not have.
        let mut running = 0u32;
        let mut offsets = Vec::with_capacity(per_cell.len());
        for n in &per_cell {
            offsets.push(running);
            running = running.checked_add(*n).ok_or(Error::GridTooLarge {
                samples: params.samples(),
            })?;
        }
        let triangles = running;
        if triangles == 0 {
            timings.prefix_ms = started.elapsed().as_secs_f64() * 1000.0;
            return Ok(GpuMesh {
                timings,
                ..GpuMesh::default()
            });
        }

        let mut offset_bytes = Vec::with_capacity(offsets.len() * 4);
        for o in &offsets {
            offset_bytes.extend_from_slice(&o.to_le_bytes());
        }
        queue.write_buffer(&counts, 0, &offset_bytes);
        timings.prefix_ms = started.elapsed().as_secs_f64() * 1000.0;

        let vertex_floats = u64::from(triangles) * 3 * 3;
        let geometry = |label| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: vertex_floats * 4,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            })
        };
        let positions = geometry("isomesh positions");
        let normals = geometry("isomesh normals");

        let started = std::time::Instant::now();
        self.dispatch(
            device, queue, &uniform, field, &counts, &positions, &normals, &self.emit, groups,
        );
        timings.emit_ms = started.elapsed().as_secs_f64() * 1000.0;

        let started = std::time::Instant::now();
        let flat_positions = read_buffer(device, queue, &positions, vertex_floats * 4)?;
        let flat_normals = read_buffer(device, queue, &normals, vertex_floats * 4)?;
        timings.geometry_readback_ms = started.elapsed().as_secs_f64() * 1000.0;
        let triple = |flat: &[f32]| -> Vec<[f32; 3]> {
            flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
        };
        Ok(GpuMesh {
            positions: triple(&flat_positions),
            normals: triple(&flat_normals),
            timings,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn dispatch(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        uniform: &wgpu::Buffer,
        field: &FieldBuffer,
        counts: &wgpu::Buffer,
        positions: &wgpu::Buffer,
        normals: &wgpu::Buffer,
        pipeline: &wgpu::ComputePipeline,
        groups: u64,
    ) {
        fn binding(index: u32, buffer: &wgpu::Buffer) -> wgpu::BindGroupEntry<'_> {
            wgpu::BindGroupEntry {
                binding: index,
                resource: buffer.as_entire_binding(),
            }
        }
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("isomesh marching cubes"),
            layout: &self.layout,
            entries: &[
                binding(0, uniform),
                binding(1, field.buffer()),
                binding(2, &self.cases),
                binding(3, counts),
                binding(4, positions),
                binding(5, normals),
            ],
        });

        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("isomesh marching cubes"),
        });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("isomesh marching cubes"),
                timestamp_writes: None,
            });
            pass.set_pipeline(pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            // `groups` was bounded against the adapter's limit by the caller.
            pass.dispatch_workgroups(groups as u32, 1, 1);
        }
        queue.submit(Some(encoder.finish()));
    }
}

#[cfg(test)]
mod tests;
