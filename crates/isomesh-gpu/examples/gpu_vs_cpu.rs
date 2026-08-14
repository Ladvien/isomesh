//! Where the GPU starts winning, and by how much — measured, not asserted.
//!
//! ```bash
//! cargo run -p isomesh-gpu --example gpu_vs_cpu --release
//! ```
//!
//! **Run it with `--release` or the answer is meaningless.** A debug-build CPU
//! extraction is 20–50× slower, which would flatter the GPU by roughly the
//! factor this example exists to measure.
//!
//! # What is being compared
//!
//! Three paths over the same field and the same grid:
//!
//! | | |
//! |---|---|
//! | **CPU** | `isomesh`'s Marching Cubes, single-threaded, evaluating the field as it goes |
//! | **GPU + read-back** | field evaluated on the GPU, extracted, vertices copied home |
//! | **GPU, no read-back** | the same without the copy — what a caller that only *draws* the result pays |
//!
//! The third column is the honest one for a renderer and the second is the
//! honest one for a CAD tool that wants the triangles in host memory. Quoting
//! either as "the GPU speedup" without saying which is how a benchmark becomes
//! folklore.
//!
//! # What this repo already measured, so you can see whether it reproduces
//!
//! GPU-011a took the 129³ path from **8.37 ms to 0.54 ms** by evaluating the
//! field on the GPU instead of uploading it — what had been 86% of the run
//! became nothing, because the samples are produced where they are read (M-155).
//! Two things about the shape of that result matter more than the ratio:
//!
//! - **It is nearly flat.** 0.22 → 0.54 ms across a **420×** increase in cells.
//!   The extractor is not remotely saturated at 129³.
//! - **The crossover is around 25³.** Below it the GPU loses on fixed overhead,
//!   and a demo that only showed 129³ would be hiding that.
//!
//! Your numbers will differ — that is the point of running it rather than
//! reading it.

use std::time::Instant;

use isomesh::{MeshBuffer, RuntimeShape3, Sdf};
use isomesh_gpu::headless::Gpu;
use isomesh_gpu::{FieldSampler, GpuField, GridParams, MarchingCubesGpu};

/// A [`GpuField`] seen as an `Sdf`, so both sides evaluate the same function.
struct FieldOf(GpuField);

impl Sdf for FieldOf {
    type Scalar = f32;
    fn sample(&self, p: [f32; 3]) -> f32 {
        self.0.sample_on_cpu(p)
    }
}

/// Median of three, warmed twice.
///
/// M-145 is this repo's own record of getting this wrong: a first call that
/// compiles pipelines and allocates was reported as the extraction time, and a
/// 10.76 ms figure went into a finding before anyone noticed it was a cold
/// start.
fn median_ms(mut run: impl FnMut()) -> f64 {
    run();
    run();
    let mut samples = [0.0f64; 3];
    for slot in &mut samples {
        let started = Instant::now();
        run();
        *slot = started.elapsed().as_secs_f64() * 1000.0;
    }
    samples.sort_by(f64::total_cmp);
    samples[1]
}

fn main() {
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
        "adapter: {} ({:?}, {:?}) driver {}",
        report.name, report.backend, report.device_type, report.driver
    );
    if cfg!(debug_assertions) {
        println!();
        println!("*** DEBUG BUILD -- these numbers are not worth reading. Use --release. ***");
    }
    println!();

    let field = GpuField::Sphere;
    let Ok(sampler) = FieldSampler::new(gpu.device()) else {
        eprintln!("field sampler");
        return;
    };
    let Ok(mc) = MarchingCubesGpu::new(gpu.device(), gpu.queue()) else {
        eprintln!("pipeline");
        return;
    };

    println!(
        "{:>6}  {:>10}  {:>12}  {:>12}  {:>12}  {:>9}  {:>9}",
        "n", "triangles", "cpu ms", "gpu+read ms", "gpu ms", "vs cpu", "no-read"
    );

    for n in [17u32, 33, 49, 65, 97, 129] {
        let cell = 3.0 / (n as f32 - 1.0);
        let Ok(grid) = GridParams::new([n; 3], [-1.5; 3], cell) else {
            eprintln!("{n}: grid");
            continue;
        };

        // The GPU path, field and all. Sampling is inside the timed region
        // because the CPU pays for its field evaluation too -- timing only the
        // extraction would be comparing the GPU's easy half against the CPU's
        // whole job.
        let mut triangles = 0usize;
        let gpu_read = median_ms(|| {
            let Ok(buffer) = sampler.sample(gpu.device(), gpu.queue(), grid, field) else {
                return;
            };
            if let Ok(mesh) = mc.extract(gpu.device(), gpu.queue(), &buffer) {
                triangles = mesh.triangle_count();
            }
        });

        let gpu_only = median_ms(|| {
            let Ok(buffer) = sampler.sample(gpu.device(), gpu.queue(), grid, field) else {
                return;
            };
            let _ = mc.extract_buffers(gpu.device(), gpu.queue(), &buffer);
        });

        let Ok(shape) = RuntimeShape3::new(grid.samples()) else {
            eprintln!("{n}: shape");
            continue;
        };
        let mut out = MeshBuffer::<f32>::new();
        // Constructed once, outside the timed closure: the extractor owns its
        // scratch precisely so re-meshing does not re-allocate ("Construct
        // once, call extract as often as you like" -- its own docs). A fresh
        // extractor per timed run would charge ~34 MB of cold allocation to
        // the CPU column and flatter the ratios this example exists to
        // measure. The warmup runs grow the scratch; the timed runs reuse it.
        let mut cpu_mc = isomesh::marching_cubes::MarchingCubes::<f32>::new();
        let cpu = median_ms(|| {
            out.reset();
            let _ = cpu_mc.extract(
                &FieldOf(field),
                &shape,
                grid.origin(),
                grid.cell_size(),
                &mut out,
            );
        });

        // Reported as a ratio in the direction that is true, rather than always
        // as "GPU speedup" -- below the crossover the GPU loses, and a table
        // that only ever printed a number greater than one would be lying at
        // the small end.
        let ratio = |a: f64, b: f64| {
            if b <= a {
                format!("{:.1}x", a / b.max(f64::EPSILON))
            } else {
                format!("{:.1}x slower", b / a.max(f64::EPSILON))
            }
        };

        println!(
            "{n:>6}  {triangles:>10}  {cpu:>12.3}  {gpu_read:>12.3}  {gpu_only:>12.3}  {:>9}  {:>9}",
            ratio(cpu, gpu_read),
            ratio(cpu, gpu_only),
        );

        if out.triangle_count() != triangles {
            println!(
                "        ^ CPU produced {} triangles and the GPU {triangles} -- that is a defect, not a tolerance",
                out.triangle_count()
            );
        }
    }

    println!();
    println!("`vs cpu` is against the read-back path, `no-read` against leaving the");
    println!("geometry on the device. A renderer pays the second and a CAD tool the first.");
}
