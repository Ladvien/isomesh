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

use crate::{Composer, Error, FieldBuffer, GridParams, PrefixScan, Result, read_buffer};

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
    /// The GPU prefix scan, including the four-byte total read-back it ends
    /// with.
    ///
    /// One field where there were two. Before GPU-010a this stage was a
    /// **read-back of every per-cell count** followed by a CPU prefix sum —
    /// 1.97 ms and 3.28 ms of a 15.03 ms run at 129³, and 8 MB copied home to
    /// add up (M-149). Both are now one dispatch chain and four bytes, so
    /// reporting them separately would be reporting a stage that no longer
    /// exists.
    pub scan_ms: f64,
    /// Classifying again and writing the triangles.
    pub emit_ms: f64,
    /// Reading positions and normals back.
    pub geometry_readback_ms: f64,
}

impl ExtractTimings {
    /// Everything above, summed.
    #[must_use]
    pub fn total_ms(&self) -> f64 {
        self.count_ms + self.scan_ms + self.emit_ms + self.geometry_readback_ms
    }

    /// The share spent moving data back to the CPU rather than computing.
    ///
    /// Only the geometry read-back counts. The scan ends with four bytes, which
    /// is a latency cost rather than a bandwidth one, and folding it in here
    /// would make a `u32` look like data movement.
    #[must_use]
    pub fn readback_share(&self) -> f64 {
        let total = self.total_ms();
        if total > 0.0 {
            self.geometry_readback_ms / total
        } else {
            0.0
        }
    }
}

/// Extraction output left in GPU memory.
///
/// What [`MarchingCubesGpu::extract_buffers`] returns and what
/// [`MeshShaderRenderer`](crate::MeshShaderRenderer) draws. The buffers carry
/// `STORAGE | COPY_SRC`, so they can be bound to a shader *or* copied home —
/// the choice is the caller's and is made after extraction rather than during
/// it.
#[derive(Debug)]
pub struct GpuGeometry {
    /// Three positions per triangle, flat `f32` triples, in cell order.
    pub positions: wgpu::Buffer,
    /// Parallel to `positions`.
    pub normals: wgpu::Buffer,
    /// Triangles in the soup. Zero means the buffers are placeholders.
    pub triangles: u32,
    /// Where the time went. `geometry_readback_ms` is zero here by definition.
    pub timings: ExtractTimings,
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
    scan: PrefixScan,
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
            scan: PrefixScan::new(device)?,
        })
    }

    /// Extract `field` and read the triangles back to the CPU.
    ///
    /// Blocks twice: once for the counts, once for the geometry. A caller that
    /// only wants to *draw* the result should use
    /// [`extract_buffers`](Self::extract_buffers) instead and skip the second
    /// wait entirely — M-149 measures that at 6.7% of the path at 129³.
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
        let geometry = self.extract_buffers(device, queue, field)?;
        if geometry.triangles == 0 {
            return Ok(GpuMesh {
                timings: geometry.timings,
                ..GpuMesh::default()
            });
        }

        let mut timings = geometry.timings;
        let floats = u64::from(geometry.triangles) * 9;
        let started = std::time::Instant::now();
        let flat_positions = read_buffer(device, queue, &geometry.positions, floats * 4)?;
        let flat_normals = read_buffer(device, queue, &geometry.normals, floats * 4)?;
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

    /// Extract `field` and leave the geometry **on the GPU**.
    ///
    /// The half of [`extract`](Self::extract) a renderer wants: it still waits
    /// once, for the per-cell counts the prefix sum needs, but the positions
    /// and normals are never copied home. Hand the buffers to
    /// [`MeshShaderRenderer`](crate::MeshShaderRenderer) and they are drawn
    /// where they were written.
    ///
    /// The remaining wait is not a detail to gloss: M-149 measures the counts
    /// read-back at **1.97 ms of 15.03** at 129³, larger than the geometry
    /// read-back this removes. Losing it needs a GPU scan (GPU-010).
    ///
    /// # Errors
    ///
    /// As [`extract`](Self::extract).
    pub fn extract_buffers(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        field: &FieldBuffer,
    ) -> Result<GpuGeometry> {
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

        // Prefix-sum on the GPU. Before GPU-010a this read every per-cell count
        // home -- 8 MB at 129^3 -- and added them up on the CPU, which M-149
        // measured at 35% of the whole path. What comes back now is the total,
        // four bytes, needed to size the geometry buffers.
        let started = std::time::Instant::now();
        let scanned = self.scan.scan(device, queue, &counts, cell_words)?;
        let triangles = scanned.total;
        timings.scan_ms = started.elapsed().as_secs_f64() * 1000.0;

        if triangles == 0 {
            // An empty surface still needs bindable buffers: wgpu rejects a
            // zero-sized binding, and a caller looping over `triangles` will
            // never read them.
            let empty = |label| {
                device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some(label),
                    size: 4,
                    usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
                    mapped_at_creation: false,
                })
            };
            return Ok(GpuGeometry {
                positions: empty("isomesh positions (empty)"),
                normals: empty("isomesh normals (empty)"),
                triangles: 0,
                timings,
            });
        }

        let vertex_floats = u64::from(triangles) * 3 * 3;
        // The CPU's `checked_add` used to catch a count that overflowed a u32;
        // a GPU scan wraps silently instead, so the guard moves to the thing
        // that actually breaks -- a binding larger than the device accepts.
        // Without this the failure is a driver-side validation error rather
        // than a named one.
        if vertex_floats * 4 > device.limits().max_storage_buffer_binding_size {
            return Err(Error::GridTooLarge {
                samples: params.samples(),
            });
        }
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
            device,
            queue,
            &uniform,
            field,
            // The scan's output, not the raw counts: `emit_cells` reads binding
            // 3 as "where does this cell's first triangle go".
            &scanned.offsets,
            &positions,
            &normals,
            &self.emit,
            groups,
        );
        timings.emit_ms = started.elapsed().as_secs_f64() * 1000.0;

        Ok(GpuGeometry {
            positions,
            normals,
            triangles,
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
