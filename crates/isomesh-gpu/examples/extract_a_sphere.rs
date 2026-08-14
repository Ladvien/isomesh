//! The whole crate in one file: field on the GPU, triangles out, checked against the CPU.
//!
//! ```bash
//! cargo run -p isomesh-gpu --example extract_a_sphere --release
//! ```
//!
//! No Bevy and no window, for the same reason [`mesh_shader_probe`] has none:
//! this crate's entire API takes `&wgpu::Device`, `&wgpu::Queue` and
//! `&mut wgpu::CommandEncoder`, so a demo that needed an engine to run would be
//! demonstrating the engine. A CAD tool with no renderer is a first-class
//! consumer here, and this is what its call sequence looks like.
//!
//! # The three steps
//!
//! 1. **Evaluate the field on the GPU.** [`FieldSampler`] fills a
//!    [`FieldBuffer`] without the samples ever crossing the bus. That is the
//!    single largest win in this crate's history — GPU-011a took a 129³
//!    extraction from 8.37 ms to 0.54 ms, because what was 86% of the run became
//!    nothing when the samples were produced where they are read.
//! 2. **Extract.** [`MarchingCubesGpu`] runs the same 256-case table `isomesh`
//!    derives at compile time, uploaded rather than transcribed — so there is no
//!    second copy of the table to disagree with the first.
//! 3. **Read it back, or don't.** `extract` returns vertices in host memory.
//!    A caller that only wants to *draw* the result should use
//!    `extract_buffers` and never wait at all.
//!
//! # Why it prints a comparison rather than a triangle count
//!
//! A GPU extractor that returns *some* triangles is easy. This one has to return
//! the *same* triangles as the CPU path, and the interesting part is where it
//! does not: normals are computed by central differences on the uploaded samples
//! rather than from the field's analytic gradient, which is a documented
//! divergence of 0.46° worst case at 17³ (M-65), not a defect. The run below
//! reports the agreement it actually achieves instead of asserting one.

use isomesh::{MeshBuffer, RuntimeShape3, Sdf};
use isomesh_gpu::headless::Gpu;
use isomesh_gpu::{FieldSampler, GpuField, GridParams, MarchingCubesGpu};

/// A [`GpuField`] seen as an `Sdf`, so the CPU comparison samples the same
/// function the shader does rather than a re-derivation of it.
struct FieldOf(GpuField);

impl Sdf for FieldOf {
    type Scalar = f32;
    fn sample(&self, p: [f32; 3]) -> f32 {
        self.0.sample_on_cpu(p)
    }
}

fn main() {
    // No software fallback, by design: a benchmark that silently ran on
    // lavapipe would report numbers three orders of magnitude off and look
    // merely slow.
    let gpu = match Gpu::new() {
        Ok(gpu) => gpu,
        Err(e) => {
            eprintln!("no GPU adapter: {e}");
            eprintln!("this crate does not fall back to a software rasteriser -- see Gpu::new");
            return;
        }
    };
    let report = gpu.report();
    println!(
        "adapter: {} ({:?}, {:?})",
        report.name, report.backend, report.device_type
    );
    println!();

    let field = GpuField::Sphere;
    let n = 65u32;
    // A sphere of radius 1 at the origin, so a grid spanning [-1.5, 1.5] holds
    // it with room for the boundary. `h` is deliberately not a power of two --
    // M-143 caught a `mul_add` here that was invisible at 0.125 and wrong at 0.1.
    let cell = 3.0 / (n as f32 - 1.0);
    let grid = match GridParams::new([n; 3], [-1.5; 3], cell) {
        Ok(grid) => grid,
        Err(e) => {
            eprintln!("grid: {e}");
            return;
        }
    };

    let sampler = match FieldSampler::new(gpu.device()) {
        Ok(sampler) => sampler,
        Err(e) => {
            eprintln!("field sampler: {e}");
            return;
        }
    };
    let buffer = match sampler.sample(gpu.device(), gpu.queue(), grid, field) {
        Ok(buffer) => buffer,
        Err(e) => {
            eprintln!("sampling: {e}");
            return;
        }
    };

    let mc = match MarchingCubesGpu::new(gpu.device(), gpu.queue()) {
        Ok(mc) => mc,
        Err(e) => {
            eprintln!("pipeline: {e}");
            return;
        }
    };
    // Warm once. The first call compiles pipelines and allocates, and reporting
    // that as the extraction time is how a GPU path acquires a reputation for
    // being slow -- M-145 is this repo's own instance of that mistake.
    if let Err(e) = mc.extract(gpu.device(), gpu.queue(), &buffer) {
        eprintln!("warm-up: {e}");
        return;
    }
    let mesh = match mc.extract(gpu.device(), gpu.queue(), &buffer) {
        Ok(mesh) => mesh,
        Err(e) => {
            eprintln!("extraction: {e}");
            return;
        }
    };

    println!("field {} at {n}^3, cell {cell:.6}", field.name());
    println!("  GPU triangles      {:>10}", mesh.triangle_count());
    println!("  GPU vertices       {:>10}", mesh.positions.len());
    println!();
    println!("  count pass         {:>10.3} ms", mesh.timings.count_ms);
    println!("  prefix scan        {:>10.3} ms", mesh.timings.scan_ms);
    println!("  emit pass          {:>10.3} ms", mesh.timings.emit_ms);
    println!(
        "  geometry read-back {:>10.3} ms   <- {:.1}% of the total, and avoidable",
        mesh.timings.geometry_readback_ms,
        mesh.timings.readback_share() * 100.0
    );
    println!("  total              {:>10.3} ms", mesh.timings.total_ms());

    // The same extraction on the CPU, through `isomesh` proper.
    let mut cpu = MeshBuffer::<f32>::new();
    let Ok(shape) = RuntimeShape3::new(grid.samples()) else {
        eprintln!("shape");
        return;
    };
    if let Err(e) = isomesh::marching_cubes::MarchingCubes::<f32>::new().extract(
        &FieldOf(field),
        &shape,
        grid.origin(),
        grid.cell_size(),
        &mut cpu,
    ) {
        eprintln!("cpu extraction: {e}");
        return;
    }

    println!();
    println!("  CPU triangles      {:>10}", cpu.triangle_count());
    println!(
        "  CPU vertices       {:>10}   <- fewer, because the CPU path shares them between cells",
        cpu.positions.len()
    );
    if cpu.triangle_count() == mesh.triangle_count() {
        println!("  -> identical triangle count: the uploaded case table is isomesh's own");
    } else {
        println!("  -> DIFFERENT triangle count, which would be a defect rather than a tolerance");
        return;
    }

    // **The two vertex buffers are not comparable index for index, and assuming
    // they were is the first thing this example got wrong.** The CPU path keys a
    // cache on the grid edge and emits each vertex once; the GPU path writes
    // three per triangle at a scanned offset and shares nothing. So the
    // comparison has to go through the CPU's index buffer, which is what puts
    // the two in the same order -- both walk cells x-fastest and both use the
    // same case table, so triangle `i` is triangle `i` on either side.
    let mut worst_position = 0.0f32;
    let mut worst_angle = 0.0f32;
    let mut sum_angle = 0.0f64;
    for corner in 0..cpu.indices.len() {
        let Some(&index) = cpu.indices.get(corner) else {
            break;
        };
        let (Some(a), Some(b)) = (
            cpu.positions.get(index as usize),
            mesh.positions.get(corner),
        ) else {
            break;
        };
        let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        worst_position = worst_position.max((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt());

        // Normals are the documented divergence: the shader has only the
        // uploaded samples and takes central differences at the cell size, where
        // the CPU asks the field for its analytic gradient. M-65 measured that
        // at 0.460 deg worst and 0.299 deg mean at 17^3, converging at h^2.
        if let (Some(na), Some(nb)) = (cpu.normals.get(index as usize), mesh.normals.get(corner)) {
            let dot = (na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2]).clamp(-1.0, 1.0);
            let angle = dot.acos().to_degrees();
            worst_angle = worst_angle.max(angle);
            sum_angle += f64::from(angle);
        }
    }
    let mean_angle = sum_angle / cpu.indices.len().max(1) as f64;

    println!();
    println!(
        "  worst vertex disagreement  {worst_position:.3e} world units, against {cell:.6} of spacing"
    );
    println!("  worst normal disagreement  {worst_angle:.3} deg");
    println!("  mean  normal disagreement  {mean_angle:.3} deg");
    println!("  -> central differences against an analytic gradient (M-65), not a defect");
}
