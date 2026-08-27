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
// The standard library's `Instant::now()` compiles for
// `wasm32-unknown-unknown` and then **panics at run time** -- the bug
// `bevy_isomesh/Cargo.toml` already records being bitten by, and part of the
// reason this crate was kept out of the web graph. `web-time-1.1.0/src/lib.rs`
// re-exports the standard library's `time` module verbatim off-wasm, so this is
// a pure textual swap with bit-identical native behaviour: one clock, not two,
// and no `cfg` fork through the timing code below. Grepping this crate for the
// std path is the gate that keeps a second clock from creeping back in, which
// is why the path is not spelled out here.
use web_time::Instant;

use crate::{Composer, Error, FieldBuffer, GridParams, PrefixScan, Result, read_bytes_many};

/// Words per case in the uploaded table: a header, then one per triangle.
const CASE_STRIDE: usize = 1 + MAX_TRIANGLES;

/// Threads per workgroup, matching `@workgroup_size(64)` in the shader.
const WORKGROUP: u32 = 64;

/// Triangles a mesh-shader workgroup emits. Must match `BATCH` in
/// `mesh_render.wgsl` and `MeshShaderRenderer`'s own constant.
const MESH_BATCH: u32 = 32;

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

/// What both extraction entry points build before their first dispatch.
///
/// Extracted at GPU-013 because the two copies had **already drifted**: `counts`
/// carried `STORAGE | COPY_SRC` in one and `STORAGE | COPY_DST | COPY_SRC` in
/// the other, `counts_bytes` was a named local in one and an inline expression
/// in the other, and both explanatory comments existed in only one — so a reader
/// of `extract_indirect` could not see why either buffer was shaped that way.
///
/// The unified `counts` drops `COPY_DST`. Nothing copies or writes into that
/// buffer anywhere in the crate — the shader is its only author — so the flag
/// was dead in the copy that had it rather than missing from the one that did
/// not.
///
/// The ticket named the wrong pair. `extract` has no prologue at all: it
/// delegates to `extract_buffers` and reads two buffers home. The duplication
/// was between `extract_indirect` and `extract_buffers`, and the ticket's stated
/// justification — that the indirect-draw budget-clamp fix touched this region —
/// is also wrong: those hunks land 15 and 74 lines past its end (M-197).
struct Prologue {
    /// The grid, as both entry points go on to report it.
    params: GridParams,
    /// Workgroups the per-cell passes dispatch.
    groups: u64,
    /// Cells, as the `u32` the shader indexes with.
    cell_words: u32,
    /// Grid parameters, already written.
    uniform: wgpu::Buffer,
    /// Per-cell triangle counts, written by pass one.
    counts: wgpu::Buffer,
    /// A minimal output binding for the pass that writes no geometry.
    placeholder: wgpu::Buffer,
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

/// An extraction that never synchronised: the triangle count is still on the
/// GPU, and the draw reads it from there.
///
/// The zero-read-back form. Its cost is that the geometry buffers were sized
/// from a **budget** rather than from the answer, so a surface with more
/// triangles than the budget is truncated — see
/// [`MarchingCubesGpu::extract_indirect`] for why that is a contract rather
/// than a silent failure.
#[derive(Debug)]
pub struct IndirectGeometry {
    /// Three positions per triangle, sized to the budget.
    pub positions: wgpu::Buffer,
    /// Parallel to `positions`.
    pub normals: wgpu::Buffer,
    /// `[group_count_x, y, z]` for `draw_mesh_tasks_indirect`, clamped to the
    /// budget — the draw covers what the emit pass actually wrote.
    pub indirect: wgpu::Buffer,
    /// The mesh shader's own uniform: the triangle count, clamped to the
    /// budget for the same reason. The un-clamped total stays in the scan's
    /// own buffer for the caller who asks whether truncation happened.
    pub draw_params: wgpu::Buffer,
    /// A one-element buffer holding the triangle count.
    ///
    /// Reading it is the only way to learn whether the budget was exceeded, and
    /// it is the caller's choice when — or whether — to pay for that.
    pub total: wgpu::Buffer,
    /// Triangles the buffers were sized for.
    pub budget: u32,
    /// Where the time went. Every read-back field is zero here by definition.
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
    /// GPU-side spans for the two compute passes, when the extractor was built
    /// by [`with_timestamps`](Self::with_timestamps) (R-069, P-71).
    ///
    /// **`None` is not a second path.** `timestamp_writes` is a field of
    /// `ComputePassDescriptor`, so passing `None` or `Some(..)` is passing data
    /// through the same dispatch code. Nothing reads this and branches.
    timestamps: Option<crate::StageTimestamps>,
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
            timestamps: None,
        })
    }

    /// The same extractor, carrying a query set so its compute passes report
    /// **GPU-side** spans.
    ///
    /// Ticket: R-069 (P-71). `ExtractTimings` is wall-clock from the CPU's side,
    /// and its own documentation says why that cannot be split finer: without a
    /// wait the submit returns immediately and every millisecond lands on
    /// whichever call happens to block first. This is the device feature that
    /// documentation names.
    ///
    /// Read the spans with [`take_timestamps`](Self::take_timestamps) after an
    /// extraction. The extraction path is unchanged — `timestamp_writes` is data
    /// in a descriptor, not a branch.
    ///
    /// # Errors
    ///
    /// [`Error::TimestampsUnsupported`] if the device was not created with
    /// [`wgpu::Features::TIMESTAMP_QUERY`] or the queue's timestamp period is
    /// zero. **It refuses rather than degrading**: an extractor that quietly
    /// reported `Instant::now()` deltas under a GPU-side name would be reporting
    /// a column it did not measure.
    pub fn with_timestamps(device: &wgpu::Device, queue: &wgpu::Queue) -> Result<Self> {
        let mut me = Self::new(device, queue)?;
        me.timestamps = Some(crate::StageTimestamps::new(device, queue)?);
        Ok(me)
    }

    /// Resolve and clear the GPU-side spans recorded since the last call.
    ///
    /// `None` if this extractor was not built with
    /// [`with_timestamps`](Self::with_timestamps) — which is a different answer
    /// from "no passes ran" and is why it is an `Option` rather than an empty
    /// [`Spans`](crate::Spans).
    ///
    /// # Errors
    ///
    /// Anything [`StageTimestamps::resolve`](crate::StageTimestamps::resolve)
    /// can report.
    pub fn take_timestamps(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
    ) -> Result<Option<crate::Spans>> {
        match &self.timestamps {
            Some(t) => Ok(Some(t.resolve(device, queue)?)),
            None => Ok(None),
        }
    }

    /// Build [`Prologue`] — the buffers both entry points need before their
    /// first dispatch.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooLarge`] if the dispatch would exceed the device's
    /// workgroup limit, or if the cell count does not fit a `u32`.
    fn prologue(
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        field: &FieldBuffer,
    ) -> Result<Prologue> {
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

        let counts = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("isomesh cell triangle counts"),
            size: u64::from(cell_words) * 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
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

        Ok(Prologue {
            params,
            groups,
            cell_words,
            uniform,
            counts,
            placeholder,
        })
    }

    /// Extract `field` and read the triangles back to the CPU.
    ///
    /// Blocks twice: once for the counts, once for the geometry. A caller that
    /// only wants to *draw* the result should use
    /// [`extract_buffers`](Self::extract_buffers) instead and skip the second
    /// wait entirely — M-149 measures that at 6.7% of the path at 129³.
    ///
    /// **"Twice" became true at GPU-013.** It said so while the geometry
    /// read-back was two separate calls, one for positions and one for normals,
    /// each with its own `poll(Wait)` — so the count was three, and the doc had
    /// been describing the intent rather than the code. Positions and normals
    /// now travel in one submission and one wait (M-197).
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
        let bytes = u64::from(geometry.triangles) * 9 * 4;
        let started = Instant::now();
        // **One submission, one device wait, for both.** Positions and normals
        // are independent and nothing between them consumes a result, so reading
        // them separately bought two full `poll(Wait)` queue drains where one
        // does -- and M-167 puts synchronisation at 82.6% of this path's time,
        // so the second drain was the expensive half of the operation.
        let flat = read_bytes_many(
            device,
            queue,
            &[(&geometry.positions, bytes), (&geometry.normals, bytes)],
        )?;
        timings.geometry_readback_ms = started.elapsed().as_secs_f64() * 1000.0;

        let triple = |raw: &[u8]| -> Vec<[f32; 3]> {
            raw.as_chunks::<12>()
                .0
                .iter()
                .map(|c| {
                    let at = |k: usize| f32::from_le_bytes([c[k], c[k + 1], c[k + 2], c[k + 3]]);
                    [at(0), at(4), at(8)]
                })
                .collect()
        };
        let (Some(positions), Some(normals)) = (flat.first(), flat.get(1)) else {
            // `read_bytes_many` returns one result per request and it was given two.
            return Err(Error::DeviceLost);
        };
        Ok(GpuMesh {
            positions: triple(positions),
            normals: triple(normals),
            timings,
        })
    }

    /// Extract `field` with **no read-back at all**, sizing the output from a
    /// budget.
    ///
    /// [`extract_buffers`](Self::extract_buffers) still waits once, for the four
    /// bytes of the triangle count, and that wait also drains every dispatch
    /// queued before it — measured at **0.375 ms of a 0.454 ms extraction** at
    /// 129³ (M-159). This form removes the synchronisation entirely: the total
    /// stays in GPU memory, a kernel turns it into
    /// `draw_mesh_tasks_indirect` arguments, and the CPU carries on recording.
    ///
    /// # The budget is a contract, and here is its cost
    ///
    /// Without the total, the geometry buffers cannot be sized from the answer.
    /// Three strategies were considered and the budget is the one that ships:
    ///
    /// - **Worst case.** `MAX_TRIANGLES` is 12 per cell, so 129³ would need
    ///   **906 MB per buffer**, 1.8 GB for both, against a measured 38,456
    ///   triangles — a factor of **190×** wasted. Viable only on a large card and
    ///   never a good idea.
    /// - **Grow and retry.** Needs to know it overflowed, which needs the total,
    ///   which is the wait this exists to remove.
    /// - **A budget the caller sets.** What this does.
    ///
    /// A surface exceeding `budget` is **truncated**, and the shader stops at
    /// the buffer's own `arrayLength` rather than running off the end. That is
    /// not a silent fallback: [`total`](IndirectGeometry::total) holds the real
    /// count, and a caller who wants certainty reads it — explicitly, when they
    /// choose, exactly as `collider::readiness` is a check the caller runs
    /// rather than a guarantee the extractor pretends to.
    ///
    /// # Errors
    ///
    /// As [`extract_buffers`](Self::extract_buffers).
    pub fn extract_indirect(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        field: &FieldBuffer,
        budget: u32,
    ) -> Result<IndirectGeometry> {
        let Prologue {
            params,
            groups,
            cell_words,
            uniform,
            counts,
            placeholder,
        } = Self::prologue(device, queue, field)?;

        let mut timings = ExtractTimings::default();
        let started = Instant::now();
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

        let started = Instant::now();
        let scanned = self
            .scan
            .scan_deferred(device, queue, &counts, cell_words)?;
        timings.scan_ms = started.elapsed().as_secs_f64() * 1000.0;

        let floats = u64::from(budget) * 9;
        if floats * 4 > device.limits().max_storage_buffer_binding_size {
            return Err(Error::GridTooLarge {
                samples: params.samples(),
            });
        }
        let geometry = |label| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: (floats * 4).max(4),
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            })
        };
        let positions = geometry("isomesh positions (budgeted)");
        let normals = geometry("isomesh normals (budgeted)");

        let started = Instant::now();
        self.dispatch(
            device,
            queue,
            &uniform,
            field,
            &scanned.offsets,
            &positions,
            &normals,
            &self.emit,
            groups,
        );
        timings.emit_ms = started.elapsed().as_secs_f64() * 1000.0;

        let indirect = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("isomesh mesh draw args"),
            size: 12,
            usage: wgpu::BufferUsages::INDIRECT
                | wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let draw_params = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("isomesh draw params"),
            size: 16,
            usage: wgpu::BufferUsages::UNIFORM
                | wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        self.scan.write_draw_args(
            device,
            queue,
            &scanned.total,
            &indirect,
            &draw_params,
            MESH_BATCH,
            budget,
        );

        Ok(IndirectGeometry {
            positions,
            normals,
            indirect,
            draw_params,
            total: scanned.total,
            budget,
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
        let Prologue {
            params,
            groups,
            cell_words,
            uniform,
            counts,
            placeholder,
        } = Self::prologue(device, queue, field)?;

        let mut timings = ExtractTimings::default();

        let started = Instant::now();
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
        let started = Instant::now();
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

        let started = Instant::now();
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
                // Data, not a branch: `None` when this extractor has no query
                // set, `Some` when it has (R-069). The dispatch below is the
                // same code either way.
                timestamp_writes: self.timestamps.as_ref().and_then(|t| {
                    t.writes(if std::ptr::eq(pipeline, &self.count) {
                        "count"
                    } else {
                        "emit"
                    })
                }),
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
