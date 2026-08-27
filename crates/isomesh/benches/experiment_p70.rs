//! **P-70 - subgroup ballot compaction, measured and deliberately not landed.**
//!
//! Ticket: R-068. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p70
//! ```
//!
//! Writes `docs/experiments/p-70.csv`.
//!
//! # C1 was unreachable at registration, and this harness says so with numbers
//!
//! `docs/measurements/gpu_vs_cpu.csv` at 129³ reads `gpu_total_ms` **8.3694**,
//! of which `scan_ms` is **0.3657** - **4.37%**. C1 asks for below 7.0 ms, i.e.
//! 1.3694 ms removed. A **free** scan leaves 8.0037. The residue is `upload_ms`
//! at 7.3240 ms, **87.5%**. Those five numbers are read from the committed CSV
//! by this bench rather than quoted, so the row cannot drift from the artefact
//! it describes.
//!
//! # Why the shader is here and not in `isomesh-gpu`
//!
//! A second WGSL path in the shipped crate is what `CLAUDE.md`'s one-path rule
//! forbids, and 4.37% capped at the literature's 1.5x - a **1.4%** improvement
//! to `gpu_total_ms` - does not buy an exception. So the subgroup scan is
//! compiled from inline WGSL **in this bench**, measured against the shipped
//! path, and not landed. `P-71`'s `E×7` is the same shape.
//!
//! # The three-way correctness anchor
//!
//! The bench carries a Hillis-Steele entry point as well, because an A/B needs
//! both arms in one module and one submission. **A copy of a shipped shader
//! inside a bench is exactly the second definition this ledger keeps finding**,
//! so it is not trusted: both arms must agree with the shipped
//! `PrefixScan::scan` **and** with `cpu_prefix_sum`, on every input size. A
//! drifted transcription is then caught by the crate's own oracle rather than by
//! a human reading two shaders side by side.

#![allow(clippy::cast_precision_loss, clippy::too_many_lines)]

mod common;

use std::time::Instant;

use isomesh_gpu::wgpu;
use isomesh_gpu::{PrefixScan, cpu_prefix_sum};

/// Elements scanned. `PrefixScan::BLOCK` is 256, so these straddle the
/// one-level / two-level boundary and the 129³ cell count.
const SIZES: [u32; 4] = [65_536, 262_144, 1_048_576, 2_097_152];

/// Repeats per arm, median taken.
const REPS: usize = 7;

/// Both arms in one module, so the A/B is one pipeline layout and one
/// submission shape.
///
/// `scan_hillis` is a transcription of the shipped `scan_blocks` entry point,
/// anchored by the three-way agreement check rather than by inspection.
/// `scan_subgroup` replaces the Hillis-Steele ladder with one
/// `subgroupExclusiveAdd` per subgroup plus a short cross-subgroup scan: at
/// `BLOCK` 256 and subgroup 32 that is 8 partials, so 16 workgroup barriers
/// become 2.
const SCAN_AB_WGSL: &str = r#"

const BLOCK: u32 = 256u;

struct Params { n: u32, pad: u32 }

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> input: array<u32>;
@group(0) @binding(2) var<storage, read_write> output: array<u32>;
@group(0) @binding(3) var<storage, read_write> block_sums: array<u32>;

var<workgroup> temp: array<u32, BLOCK>;
// One slot per subgroup at the smallest subgroup size WGSL permits (4), which
// bounds the count at BLOCK / 4 = 64. Sized for the bound rather than for this
// adapter, so the module validates everywhere.
var<workgroup> partials: array<u32, 64>;

@compute @workgroup_size(256)
fn scan_hillis(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_index) tid: u32,
    @builtin(workgroup_id) wg: vec3<u32>,
) {
    let i = gid.x;
    var value = 0u;
    if (i < params.n) { value = input[i]; }
    temp[tid] = value;
    workgroupBarrier();

    for (var offset = 1u; offset < BLOCK; offset = offset << 1u) {
        var addend = 0u;
        if (tid >= offset) { addend = temp[tid - offset]; }
        workgroupBarrier();
        if (tid >= offset) { temp[tid] = temp[tid] + addend; }
        workgroupBarrier();
    }

    if (i < params.n) { output[i] = temp[tid] - value; }
    if (tid == BLOCK - 1u) { block_sums[wg.x] = temp[tid]; }
}

@compute @workgroup_size(256)
fn scan_subgroup(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_index) tid: u32,
    @builtin(workgroup_id) wg: vec3<u32>,
    @builtin(subgroup_size) sg_size: u32,
    @builtin(subgroup_invocation_id) sg_lane: u32,
) {
    let i = gid.x;
    var value = 0u;
    if (i < params.n) { value = input[i]; }

    // Lane-local exclusive scan, in registers, no barrier.
    let lane_prefix = subgroupExclusiveAdd(value);
    // The subgroup's total. `subgroupAdd` is a reduction that lands the same
    // value in every lane, so no broadcast is needed -- and naga requires a
    // `subgroupBroadcast` lane index to be a CONST-EXPRESSION, which
    // `sg_size - 1u` is not. That constraint is why this is two wave ops rather
    // than one plus a broadcast.
    let sg_total = subgroupAdd(value);

    let sg_id = tid / sg_size;
    if (sg_lane == 0u) { partials[sg_id] = sg_total; }
    workgroupBarrier();

    // Serial scan over the partials by one lane. At BLOCK 256 there are at most
    // 64 of them and in practice 8; a ladder here would cost more barriers than
    // it saves adds, which is the whole reason the ladder was worth removing
    // from the 256-element case.
    if (tid == 0u) {
        var running = 0u;
        let count = BLOCK / sg_size;
        for (var s = 0u; s < count; s = s + 1u) {
            let t = partials[s];
            partials[s] = running;
            running = running + t;
        }
        temp[0] = running;
    }
    workgroupBarrier();

    if (i < params.n) { output[i] = partials[sg_id] + lane_prefix; }
    if (tid == 0u) { block_sums[wg.x] = temp[0]; }
}
"#;

/// One arm at one size.
struct Row {
    arm: &'static str,
    elements: u32,
    scan_ms: f64,
    matches_shipped: bool,
    matches_cpu: bool,
}

/// Deterministic counts, in the shape the extractor produces: mostly zero, a
/// few small positives.
///
/// A uniform random fill would exercise the scan's arithmetic and not its data:
/// `M-337` measured 97% of cells producing nothing on a sphere at 128³, so a
/// fixture where every element is non-zero is a fixture the extractor never
/// hands it.
fn counts(n: u32) -> Vec<u32> {
    let mut state = 0x2026_u64 ^ 0x5EED_1234;
    (0..n)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            // 3% non-zero, values 1..=5 -- the active-cell density M-337
            // measured, and small enough that the total cannot overflow u32.
            if (state >> 40) % 100 < 3 {
                1 + ((state >> 20) % 5) as u32
            } else {
                0
            }
        })
        .collect()
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-70");

    // ── the committed artefact C1 is denominated in ───────────────────────────
    //
    // Read rather than quoted. A number in a comment is a number that drifts.
    let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let measured = std::fs::read_to_string(root.join("docs/measurements/gpu_vs_cpu.csv"))
        .expect("docs/measurements/gpu_vs_cpu.csv");
    let mut header: Vec<&str> = Vec::new();
    let mut gpu_total_129 = f64::NAN;
    let mut scan_129 = f64::NAN;
    let mut upload_129 = f64::NAN;
    for line in measured.lines() {
        let cells: Vec<&str> = line.split(',').collect();
        if header.is_empty() {
            header = cells;
            continue;
        }
        let get = |name: &str| -> Option<f64> {
            header
                .iter()
                .position(|h| *h == name)
                .and_then(|i| cells.get(i))
                .and_then(|v| v.parse().ok())
        };
        if get("samples") == Some(129.0) {
            gpu_total_129 = get("gpu_total_ms").unwrap_or(f64::NAN);
            scan_129 = get("scan_ms").unwrap_or(f64::NAN);
            upload_129 = get("upload_ms").unwrap_or(f64::NAN);
        }
    }
    assert!(
        gpu_total_129.is_finite() && scan_129.is_finite(),
        "no 129 row in gpu_vs_cpu.csv, so C1 has no denominator"
    );
    let scan_share = scan_129 / gpu_total_129;
    let free_scan_total = gpu_total_129 - scan_129;

    println!("C1's denominator, read from docs/measurements/gpu_vs_cpu.csv at 129³:");
    println!("  gpu_total_ms {gpu_total_129:.4}");
    println!(
        "  upload_ms    {upload_129:.4}  ({:.2}% -- the residue)",
        100.0 * upload_129 / gpu_total_129
    );
    println!("  scan_ms      {scan_129:.4}  ({:.2}%)", 100.0 * scan_share);
    println!("  a FREE scan leaves {free_scan_total:.4} ms; C1 asks for < 7.0\n");

    // ── the device, and C3's own subject ─────────────────────────────────────
    let gpu = isomesh_gpu::headless::Gpu::with_subgroups()
        .expect("a device with SUBGROUP; P-70 is VOID without it");
    let device = gpu.device();
    let queue = gpu.queue();
    // The adapter's own numbers, not a guess. `min == max` on every desktop
    // adapter this crate has met; both are printed so a machine where they
    // differ is visible rather than averaged.
    let report = gpu.report();
    let subgroup_size = if report.subgroup_min_size == report.subgroup_max_size {
        report.subgroup_min_size.to_string()
    } else {
        format!("{}-{}", report.subgroup_min_size, report.subgroup_max_size)
    };
    let has_subgroup = device.features().contains(wgpu::Features::SUBGROUP);
    assert!(
        has_subgroup,
        "VOID: the device does not advertise SUBGROUP, so every clause would be about code that \
         did not run"
    );
    println!("adapter subgroup size {subgroup_size}, SUBGROUP advertised {has_subgroup}\n");

    let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("p70 scan a/b"),
        source: wgpu::ShaderSource::Wgsl(SCAN_AB_WGSL.into()),
    });
    let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("p70"),
        entries: &[
            wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 1,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: true },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 2,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: false },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 3,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: false },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
        ],
    });
    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("p70"),
        bind_group_layouts: &[Some(&layout)],
        immediate_size: 0,
    });
    let make = |entry: &str| {
        device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(entry),
            layout: Some(&pipeline_layout),
            module: &module,
            entry_point: Some(entry),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        })
    };
    let hillis = make("scan_hillis");
    let subgroup = make("scan_subgroup");

    let shipped = PrefixScan::new(device).expect("shipped scan");
    let mut rows: Vec<Row> = Vec::new();

    for n in SIZES {
        let data = counts(n);
        let (cpu_offsets, _cpu_total) = cpu_prefix_sum(&data);

        // The shipped path's answer on the same input, for the three-way anchor.
        let input = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("counts"),
            size: u64::from(n) * 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&input, 0, &bytemuck_cast(&data));
        let shipped_out = shipped
            .scan(device, queue, &input, n)
            .expect("shipped scan");
        let shipped_offsets =
            isomesh_gpu::read_buffer_u32(device, queue, &shipped_out.offsets, u64::from(n) * 4)
                .expect("read shipped offsets");

        let groups = n.div_ceil(256);
        let params = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("params"),
            size: 8,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&params, 0, &[n.to_le_bytes(), [0u8; 4]].concat());
        let out = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("out"),
            size: u64::from(n) * 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let sums = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("sums"),
            size: u64::from(groups.max(1)) * 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let bind = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("p70"),
            layout: &layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: input.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: out.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: sums.as_entire_binding(),
                },
            ],
        });

        for (arm, pipeline) in [("hillis", &hillis), ("subgroup", &subgroup)] {
            // Warm, then median of REPS. One dispatch per timing so the number
            // is a block scan and not a pipeline creation.
            let mut samples: Vec<f64> = Vec::with_capacity(REPS);
            for rep in 0..=REPS {
                let started = Instant::now();
                let mut enc =
                    device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
                {
                    let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                        label: Some(arm),
                        timestamp_writes: None,
                    });
                    pass.set_pipeline(pipeline);
                    pass.set_bind_group(0, &bind, &[]);
                    pass.dispatch_workgroups(groups, 1, 1);
                }
                queue.submit(Some(enc.finish()));
                let _ = device.poll(wgpu::PollType::Wait {
                    submission_index: None,
                    timeout: None,
                });
                if rep > 0 {
                    samples.push(started.elapsed().as_nanos() as f64 / 1e6);
                }
            }
            samples.sort_unstable_by(|a, b| a.partial_cmp(b).expect("finite"));
            let scan_ms = samples[REPS / 2];

            // Per-block exclusive scan, which is what both entry points write.
            // The shipped `scan` writes a *global* exclusive scan, so the
            // comparison is per block: offset within the block.
            let got = isomesh_gpu::read_buffer_u32(device, queue, &out, u64::from(n) * 4)
                .expect("read out");
            let mut matches_cpu = true;
            let mut matches_shipped = true;
            for b in 0..groups as usize {
                let base = b * 256;
                let end = (base + 256).min(n as usize);
                let block_base = cpu_offsets[base];
                for i in base..end {
                    if got[i] != cpu_offsets[i] - block_base {
                        matches_cpu = false;
                    }
                    if got[i] + shipped_offsets[base] != shipped_offsets[i] {
                        matches_shipped = false;
                    }
                }
            }

            println!(
                "{arm:>9} n={n:>8} {scan_ms:>8.4} ms  cpu={matches_cpu} shipped={matches_shipped}"
            );
            rows.push(Row {
                arm,
                elements: n,
                scan_ms,
                matches_shipped,
                matches_cpu,
            });
        }
    }

    // ── controls ─────────────────────────────────────────────────────────────
    for r in &rows {
        assert!(
            r.matches_cpu,
            "{} at n={} disagrees with cpu_prefix_sum, so its time measures the wrong \
             computation",
            r.arm, r.elements
        );
        assert!(
            r.matches_shipped,
            "{} at n={} disagrees with the shipped PrefixScan, so the bench's transcription has \
             drifted from the crate",
            r.arm, r.elements
        );
    }

    // ── verdict ──────────────────────────────────────────────────────────────
    let mut speedups: Vec<(u32, f64)> = Vec::new();
    for n in SIZES {
        let h = rows
            .iter()
            .find(|r| r.arm == "hillis" && r.elements == n)
            .expect("hillis row")
            .scan_ms;
        let s = rows
            .iter()
            .find(|r| r.arm == "subgroup" && r.elements == n)
            .expect("subgroup row")
            .scan_ms;
        speedups.push((n, h / s));
        println!("speedup at n={n}: {:.4}x", h / s);
    }
    let best = speedups.iter().map(|(_, s)| *s).fold(0.0f64, f64::max);

    // C1: unreachable by construction, and the harness computes it rather than
    // asserting it. Even a free scan leaves `free_scan_total`.
    let c1 = free_scan_total < 7.0;
    // C2: both arms agree with the shipped path and the CPU oracle at every
    // size, which is bit-identity on the quantity being changed. The registered
    // "all eight fields" is NOT what this measures, and the entry says so.
    let c2 = rows.iter().all(|r| r.matches_shipped && r.matches_cpu);
    // C3: the fallback is the shipped path, and it is exercised by every GPU
    // test in the crate. What this harness adds is that a device WITH subgroups
    // still runs the non-subgroup arm and gets the same answer.
    let c3 = rows
        .iter()
        .any(|r| r.arm == "hillis" && r.matches_shipped && r.matches_cpu);

    println!(
        "\nC1 129³ below 7.0 ms: a FREE scan leaves {free_scan_total:.4} ms -> {}",
        if c1 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C2 bit-identical on the changed quantity, every size -> {}",
        if c2 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C3 the non-subgroup arm runs and agrees on a subgroup-capable device -> {}",
        if c3 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "\nbest subgroup speedup on the scan itself: {best:.4}x. Applied to scan_ms \
         {scan_129:.4}, gpu_total_ms goes {gpu_total_129:.4} -> {:.4}, i.e. {:.2}%.",
        gpu_total_129 - scan_129 + scan_129 / best,
        100.0 * (scan_129 - scan_129 / best) / gpu_total_129
    );
    println!(
        "NOT LANDED: a second WGSL path in the shipped crate is what CLAUDE.md's one-path rule \
         forbids, and that percentage does not buy an exception."
    );

    common::experiment::run(prereg, |run| {
        for r in &rows {
            let speedup = if r.arm == "subgroup" {
                speedups
                    .iter()
                    .find(|(n, _)| *n == r.elements)
                    .map_or(f64::NAN, |(_, s)| *s)
            } else {
                1.0
            };
            run.record(&[
                ("arm", r.arm.to_string()),
                ("elements", r.elements.to_string()),
                ("scan_ms", format!("{:.6}", r.scan_ms)),
                ("speedup_vs_hillis", format!("{speedup:.6}")),
                ("matches_shipped_scan", r.matches_shipped.to_string()),
                ("matches_cpu_oracle", r.matches_cpu.to_string()),
                ("subgroup_size", subgroup_size.clone()),
                ("subgroup_feature", has_subgroup.to_string()),
                ("gpu_total_ms_129", format!("{gpu_total_129:.6}")),
                ("scan_share_of_total", format!("{scan_share:.6}")),
                (
                    "reachable_total_if_scan_free",
                    format!("{free_scan_total:.6}"),
                ),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
            ]);
        }
    });
}

/// `&[u32]` as bytes, little-endian, without pulling in `bytemuck`.
///
/// `isomesh` has two normal dependencies and rule 3 counts them; a cast this
/// short does not earn a third. Correct on every target this crate supports,
/// all of which are little-endian, and `to_le_bytes` makes that explicit rather
/// than assumed.
fn bytemuck_cast(v: &[u32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(v.len() * 4);
    for x in v {
        out.extend_from_slice(&x.to_le_bytes());
    }
    out
}
