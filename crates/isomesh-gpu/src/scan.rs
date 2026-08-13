//! Exclusive prefix sum on the GPU, so the counts never come home.
//!
//! # What this replaces, and what it is worth
//!
//! Before this, extraction did: count pass → **read every per-cell count back**
//! → prefix-sum on the CPU → upload the offsets → emit pass. M-149 measured the
//! two CPU-side stages at **3.28 ms and 1.97 ms of a 15.03 ms** run at 129³ —
//! 35% of the whole GPU path — and the read-back alone moves **4 bytes per
//! cell**, which is 8 MB at that size, purely to add the numbers up.
//!
//! With the scan here the CPU needs **four bytes**: the grand total, to size the
//! output buffer. That last read-back disappears too once the emit pass is
//! dispatched indirectly and the geometry buffer is sized without it, which is
//! GPU-010b.
//!
//! **Measured (M-150):** the two stages it replaces were `5.24 ms` at 129³ and
//! the scan is `0.37 ms` — **14× faster** — taking the whole GPU path from
//! `15.01` to `9.65 ms` and from `1.27×` to `1.95×` ahead of a single-threaded
//! CPU. It is a small *loss* below about 25³, where two extra dispatches cost
//! more than adding four thousand numbers up on the CPU.
//!
//! # The shape of it
//!
//! Hierarchical, with the level loop on the CPU because the depth depends on the
//! cell count — three levels at 129³, one below 256 cells:
//!
//! 1. Each workgroup scans its own block of [`PrefixScan::BLOCK`] elements and
//!    publishes the
//!    block's total.
//! 2. The block totals are themselves scanned, one level up, until a level fits
//!    in a single workgroup.
//! 3. Each level then adds its parent's exclusive scan back down.
//!
//! The grand total falls out for free: the deepest level has exactly one block,
//! so its single published total is the sum of everything below it.

use crate::{Composer, GridParams, Result, read_buffer_u32};

/// Elements one workgroup scans. Must match `BLOCK` and `@workgroup_size` in
/// `scan.wgsl`.
const BLOCK: u32 = 256;

/// One level of the hierarchy: the elements it scans and where it puts them.
#[derive(Debug)]
struct Level {
    /// Workgroups, i.e. blocks, i.e. how many totals this level publishes.
    ///
    /// The element count it came from lives only in [`Level::params`], because
    /// that is the copy the shader reads — keeping a second one in Rust would
    /// be two numbers that must agree and one place for them to stop agreeing.
    blocks: u32,
    /// This level's exclusive scan, `n` elements.
    out: wgpu::Buffer,
    /// Per-block totals, `blocks` elements. Scanned by the level above.
    sums: wgpu::Buffer,
    /// `n`, as the shader's uniform.
    params: wgpu::Buffer,
}

/// What a scan produced.
#[derive(Debug)]
pub struct ScanOutput {
    /// The exclusive prefix sum, one `u32` per input element.
    pub offsets: wgpu::Buffer,
    /// The sum of every input element.
    pub total: u32,
    /// Levels the hierarchy needed. One below [`PrefixScan::BLOCK`] elements,
    /// three at 129³.
    ///
    /// Exposed because a test that does not know the depth cannot tell a
    /// single-level scan from a multi-level one, and the multi-level path is
    /// the one with a cross-block bug to have.
    pub levels: usize,
}

/// The compiled scan pipelines.
#[derive(Debug)]
pub struct PrefixScan {
    scan: wgpu::ComputePipeline,
    add: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
}

impl PrefixScan {
    /// Elements one workgroup scans.
    pub const BLOCK: u32 = BLOCK;

    /// Compile both entry points.
    ///
    /// # Errors
    ///
    /// Propagates composition failure, which would mean this crate's own
    /// shaders are inconsistent.
    pub fn new(device: &wgpu::Device) -> Result<Self> {
        let source = Composer::with_builtins().compose("scan", &[])?;
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("isomesh prefix scan"),
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
            label: Some("isomesh prefix scan bindings"),
            entries: &[
                entry(0, wgpu::BufferBindingType::Uniform),
                entry(1, read),
                entry(2, write),
                entry(3, write),
            ],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("isomesh prefix scan layout"),
            bind_group_layouts: &[Some(&layout)],
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

        Ok(Self {
            scan: compile("scan_blocks"),
            add: compile("add_block_offsets"),
            layout,
        })
    }

    /// Exclusive prefix sum of the first `n` `u32`s in `counts`.
    ///
    /// `counts` needs `STORAGE`. The returned `offsets` buffer carries
    /// `STORAGE | COPY_SRC`, so it can be bound to the emit pass or read back
    /// for comparison.
    ///
    /// Blocks once, for four bytes: the grand total, which the caller needs to
    /// size whatever the scan is feeding.
    ///
    /// # Errors
    ///
    /// Anything the four-byte read-back can report.
    pub fn scan(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        counts: &wgpu::Buffer,
        n: u32,
    ) -> Result<ScanOutput> {
        // Sizes first, buffers second, dispatches third. Working out the whole
        // hierarchy before allocating any of it keeps the borrow of the level
        // below out of the loop that builds the level above.
        let mut sizes = Vec::new();
        let mut at = n.max(1);
        loop {
            let blocks = at.div_ceil(BLOCK);
            sizes.push((at, blocks));
            if blocks <= 1 {
                break;
            }
            at = blocks;
        }

        let levels: Vec<Level> = sizes
            .iter()
            .map(|&(n, blocks)| {
                let storage = |label: &str, elements: u32| {
                    device.create_buffer(&wgpu::BufferDescriptor {
                        label: Some(label),
                        size: u64::from(elements.max(1)) * 4,
                        usage: wgpu::BufferUsages::STORAGE
                            | wgpu::BufferUsages::COPY_SRC
                            | wgpu::BufferUsages::COPY_DST,
                        mapped_at_creation: false,
                    })
                };
                let params = device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some("isomesh scan params"),
                    size: 16,
                    usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
                    mapped_at_creation: false,
                });
                // One `u32`, padded to a 16-byte uniform binding.
                let mut bytes = [0u8; 16];
                bytes[..4].copy_from_slice(&n.to_le_bytes());
                queue.write_buffer(&params, 0, &bytes);
                Level {
                    blocks,
                    out: storage("isomesh scan output", n),
                    sums: storage("isomesh scan block sums", blocks),
                    params,
                }
            })
            .collect();

        // Up: scan each level, publishing block totals to the level above.
        for (index, level) in levels.iter().enumerate() {
            let input = if index == 0 {
                counts
            } else {
                &levels[index - 1].sums
            };
            self.dispatch(
                device,
                queue,
                &self.scan,
                level,
                input,
                &level.out,
                &level.sums,
            );
        }

        // Down: add each parent's exclusive scan into its children's blocks.
        // Skipped entirely when there is one level, because a single block's
        // exclusive scan is already the answer.
        for index in (0..levels.len().saturating_sub(1)).rev() {
            let level = &levels[index];
            self.dispatch(
                device,
                queue,
                &self.add,
                level,
                &level.sums,
                &level.out,
                &levels[index + 1].out,
            );
        }

        // The deepest level has exactly one block, so the single total it
        // published is the sum of everything. Four bytes, against 4 bytes per
        // cell before.
        let deepest = levels.last().expect("the loop pushes at least one level");
        let total = read_buffer_u32(device, queue, &deepest.sums, 4)?
            .first()
            .copied()
            .unwrap_or(0);

        let offsets = levels
            .into_iter()
            .next()
            .expect("the loop pushes at least one level")
            .out;
        Ok(ScanOutput {
            offsets,
            total,
            levels: sizes.len(),
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn dispatch(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        pipeline: &wgpu::ComputePipeline,
        level: &Level,
        input: &wgpu::Buffer,
        output: &wgpu::Buffer,
        block_sums: &wgpu::Buffer,
    ) {
        fn entry(binding: u32, buffer: &wgpu::Buffer) -> wgpu::BindGroupEntry<'_> {
            wgpu::BindGroupEntry {
                binding,
                resource: buffer.as_entire_binding(),
            }
        }
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("isomesh prefix scan"),
            layout: &self.layout,
            entries: &[
                entry(0, &level.params),
                entry(1, input),
                entry(2, output),
                entry(3, block_sums),
            ],
        });
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("isomesh prefix scan"),
        });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("isomesh prefix scan"),
                timestamp_writes: None,
            });
            pass.set_pipeline(pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(level.blocks, 1, 1);
        }
        queue.submit(Some(encoder.finish()));
    }
}

/// The exclusive prefix sum, on the CPU.
///
/// The reference the GPU scan is asserted against, and the implementation the
/// extraction used before GPU-010a. Kept because a parallel scan that is subtly
/// wrong produces a mesh that looks entirely plausible, so the comparison has to
/// be against something obviously right rather than against a second clever
/// implementation.
#[must_use]
pub fn cpu_prefix_sum(counts: &[u32]) -> (Vec<u32>, u32) {
    let mut running = 0u32;
    let mut offsets = Vec::with_capacity(counts.len());
    for n in counts {
        offsets.push(running);
        running = running.wrapping_add(*n);
    }
    (offsets, running)
}

/// Cells a grid scans, as the scan sees them.
#[must_use]
pub fn scan_length(params: GridParams) -> u32 {
    u32::try_from(params.cell_count()).unwrap_or(u32::MAX)
}

#[cfg(test)]
mod tests;
