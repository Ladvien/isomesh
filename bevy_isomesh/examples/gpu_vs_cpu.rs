//! E-302 — both extractors live, and where the GPU's time actually goes.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example gpu_vs_cpu --release
//! ```
//!
//! `[` `]` resolution, `1`–`5` field, `S` run the sweep, `G`/`C` show the GPU
//! or CPU mesh. Drag to orbit, scroll to zoom.
//!
//! # The ticket expected the gap to close at small grids. It does not close —
//! it never opens
//!
//! GPU-006 asks to *"watch the gap close at small grids: launch overhead made
//! visible"*, on the reasonable assumption that the GPU wins at large grids and
//! loses at small ones. Measured, **the GPU is behind at every resolution this
//! runs at**, and the reason is not launch overhead.
//!
//! It is **read-back**. The HUD breaks one extraction into its five parts, and
//! the two that copy memory back to the CPU dominate. That distinction is the
//! whole point of the demo, because it decides who should use this at all:
//!
//! - A consumer that **renders from GPU memory** never pays read-back. For them
//!   the number that matters is `count + emit`.
//! - A consumer that needs a **collider**, or a validity check, or anything the
//!   CPU has to look at, pays all of it. M-003 measured that the collider check
//!   is already 45% of a usable mesh; a GPU path that has to come home first
//!   adds to that rather than replacing it.
//!
//! # And the field never moves to the GPU, which caps this from the start
//!
//! `FieldBuffer::sampled` evaluates the SDF **on the CPU** and uploads the
//! samples. So this path does not remove field evaluation from the CPU's
//! budget — it adds an upload to it. M-136 measured field evaluation at 65–74%
//! of the whole job on `fbm_terrain` and 13% on `sphere`, so on exactly the
//! workload where a GPU would help most, this design helps least. Evaluating
//! the field in the shader is what would change that, and it is not this
//! ticket.
//!
//! The upload is timed separately for that reason: it is the line item that
//! says so.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy::render::renderer::{RenderDevice, RenderQueue};
use common::{CommonPlugin, DemoStats, ViewFlags};
use isomesh::fields::{BoxExact, Sphere, ThinPlate, Torus, csg_difference};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};
use isomesh_gpu::{
    ExtractTimings, FieldBuffer, FieldSampler, GpuField, GridParams, MarchingCubesGpu,
};

/// Samples per axis, unless `ISOMESH_SAMPLES` says otherwise.
const DEFAULT_SAMPLES: u32 = 65;

/// Resolutions the sweep visits. Powers-of-two spacing where possible, so the
/// comparison is not also measuring M-143's rounding difference.
const SWEEP: [u32; 6] = [17, 33, 49, 65, 97, 129];

/// The GPU-side equivalent of the example's field index, where one exists.
fn gpu_field_at(index: usize) -> Option<GpuField> {
    match index % 5 {
        0 => Some(GpuField::Sphere),
        1 => Some(GpuField::Torus),
        2 => Some(GpuField::BoxExact),
        // csg_difference and thin_plate have no GPU implementation yet.
        _ => None,
    }
}

/// The fields `1`-`5` select.
fn field_at(index: usize) -> (&'static str, Box<dyn Sdf<Scalar = f32> + Send + Sync>) {
    match index % 5 {
        0 => ("sphere", Box::new(Sphere::<f32>::canonical())),
        1 => ("torus", Box::new(Torus::<f32>::canonical())),
        2 => ("box_exact", Box::new(BoxExact::<f32>::canonical())),
        3 => ("csg_difference", Box::new(csg_difference::<f32>())),
        _ => ("thin_plate", Box::new(ThinPlate::<f32>::canonical())),
    }
}

/// One resolution, both extractors.
#[derive(Clone, Copy, Default)]
struct Point {
    samples: u32,
    triangles: usize,
    cpu_ms: f64,
    upload_ms: f64,
    /// The whole path with the field **produced on the GPU** instead of
    /// uploaded: the sampling pass plus the extraction, end to end.
    ///
    /// One number rather than a breakdown, because the sampling pass cannot be
    /// timed on its own from the CPU — `sample()` submits and returns, so a
    /// clock around it reads zero regardless of what the GPU does. What can be
    /// timed is a span that ends in a wait, and `extract_buffers` supplies one.
    /// `None` when the field has no GPU implementation — GPU-011a covers four
    /// of the seven.
    gpu_field_ms: Option<f64>,
    gpu: ExtractTimings,
}

impl Point {
    /// Everything a caller pays to get triangles onto the CPU.
    fn gpu_total_ms(&self) -> f64 {
        self.upload_ms + self.gpu.total_ms()
    }

    /// What a caller pays if the mesh stays on the GPU.
    fn gpu_compute_ms(&self) -> f64 {
        self.gpu.count_ms + self.gpu.emit_ms
    }

    /// The same, with the field produced on the GPU instead of uploaded.
    fn gpu_field_total_ms(&self) -> Option<f64> {
        self.gpu_field_ms
    }
}

/// Run the sweep at startup and print it, from `ISOMESH_SWEEP=1`.
///
/// The same reason `game_dig` has `ISOMESH_AUTOCARVE`: the product of this
/// example is a table of numbers, and a table that can only be produced by
/// pressing a key cannot be regenerated from a command line or checked in CI.
#[derive(Resource)]
struct AutoSweep(bool);

#[derive(Resource, Default)]
struct State {
    field: &'static str,
    samples: u32,
    current: Point,
    sweep: Vec<Point>,
    show_cpu: bool,
}

#[derive(Resource)]
struct Gpu(MarchingCubesGpu, FieldSampler);

#[derive(Component)]
struct Surface;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-302 gpu vs cpu".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_resource::<State>()
        .insert_resource(AutoSweep(
            std::env::var("ISOMESH_SWEEP").is_ok_and(|v| v != "0"),
        ))
        .add_systems(Startup, setup)
        .add_systems(Update, (keys, run, report).chain())
        .run();
}

fn setup(
    mut commands: Commands,
    device: Res<RenderDevice>,
    queue: Res<RenderQueue>,
    mut state: ResMut<State>,
) {
    let sampler = FieldSampler::new(device.wgpu_device())
        .expect("the field sampler is this crate's own shader");
    match MarchingCubesGpu::new(device.wgpu_device(), &queue) {
        Ok(pipeline) => commands.insert_resource(Gpu(pipeline, sampler)),
        // No CPU fallback: a comparison demo that silently ran one side twice
        // would report a ratio of 1.0 and look like a result.
        Err(why) => panic!("could not build the GPU pipeline on Bevy's device: {why}"),
    }
    state.samples = common::samples_override().unwrap_or(DEFAULT_SAMPLES);
}

/// Run both extractors once at `samples`, timing each part.
#[allow(clippy::too_many_arguments)]
fn measure(
    gpu: &MarchingCubesGpu,
    sampler: &FieldSampler,
    gpu_field: Option<GpuField>,
    device: &RenderDevice,
    queue: &RenderQueue,
    field: &dyn Sdf<Scalar = f32>,
    samples: u32,
) -> Option<Point> {
    let cell = 4.0f32 / (samples - 1) as f32;
    let grid = GridParams::new([samples; 3], [-2.0; 3], cell).ok()?;

    // Timed apart from extraction because it is the line item that says the
    // field never moved to the GPU: this is a CPU evaluation plus an upload.
    let started = Instant::now();
    let buffer = FieldBuffer::sampled(device.wgpu_device(), queue, grid, &field).ok()?;
    let upload_ms = started.elapsed().as_secs_f64() * 1000.0;

    // The same extraction with the field produced by a compute pass rather than
    // uploaded. Timed as one span ending in `extract_buffers`, which waits.
    let gpu_field_ms = gpu_field.map(|which| {
        let started = Instant::now();
        if let Ok(produced) = sampler.sample(device.wgpu_device(), queue, grid, which) {
            let _ = gpu.extract_buffers(device.wgpu_device(), queue, &produced);
        }
        started.elapsed().as_secs_f64() * 1000.0
    });

    let mesh = gpu.extract(device.wgpu_device(), queue, &buffer).ok()?;

    let started = Instant::now();
    let mut cpu = MeshBuffer::<f32>::new();
    let shape = RuntimeShape3::new([samples; 3]).ok()?;
    MarchingCubes::<f32>::new()
        .extract(&field, &shape, [-2.0; 3], cell, &mut cpu)
        .ok()?;
    let cpu_ms = started.elapsed().as_secs_f64() * 1000.0;

    Some(Point {
        samples,
        triangles: mesh.triangle_count(),
        cpu_ms,
        upload_ms,
        gpu_field_ms,
        gpu: mesh.timings,
    })
}

/// [`measure`] three times, component-wise median.
///
/// Three because a single shot of a sub-millisecond dispatch is mostly noise,
/// and the median of three is what `bench_stage_breakdown` uses for the same
/// reason. Component-wise rather than picking the run with the median total:
/// the breakdown is the product here, and a total can be median while every
/// part of it is an outlier.
#[allow(clippy::too_many_arguments)]
fn measure_median(
    gpu: &MarchingCubesGpu,
    sampler: &FieldSampler,
    gpu_field: Option<GpuField>,
    device: &RenderDevice,
    queue: &RenderQueue,
    field: &dyn Sdf<Scalar = f32>,
    samples: u32,
) -> Option<Point> {
    let runs: Vec<Point> = (0..3)
        .filter_map(|_| measure(gpu, sampler, gpu_field, device, queue, field, samples))
        .collect();
    if runs.len() < 3 {
        return None;
    }
    let median = |mut xs: [f64; 3]| {
        xs.sort_by(f64::total_cmp);
        xs[1]
    };
    let pick = |f: fn(&Point) -> f64| median([f(&runs[0]), f(&runs[1]), f(&runs[2])]);
    Some(Point {
        samples,
        triangles: runs[0].triangles,
        cpu_ms: pick(|p| p.cpu_ms),
        upload_ms: pick(|p| p.upload_ms),
        gpu_field_ms: runs[0]
            .gpu_field_ms
            .map(|_| median([0, 1, 2].map(|i| runs[i].gpu_field_ms.unwrap_or(0.0)))),
        gpu: ExtractTimings {
            count_ms: pick(|p| p.gpu.count_ms),
            scan_ms: pick(|p| p.gpu.scan_ms),
            emit_ms: pick(|p| p.gpu.emit_ms),
            geometry_readback_ms: pick(|p| p.gpu.geometry_readback_ms),
        },
    })
}

fn keys(keyboard: Res<ButtonInput<KeyCode>>, mut state: ResMut<State>) {
    if keyboard.just_pressed(KeyCode::BracketLeft) {
        state.samples = (state.samples - 16).max(17);
    }
    if keyboard.just_pressed(KeyCode::BracketRight) {
        state.samples = (state.samples + 16).min(161);
    }
    if keyboard.just_pressed(KeyCode::KeyC) {
        state.show_cpu = true;
    }
    if keyboard.just_pressed(KeyCode::KeyG) {
        state.show_cpu = false;
    }
}

#[allow(clippy::too_many_arguments)]
fn run(
    mut commands: Commands,
    gpu: Option<Res<Gpu>>,
    device: Res<RenderDevice>,
    queue: Res<RenderQueue>,
    keyboard: Res<ButtonInput<KeyCode>>,
    flags: Res<ViewFlags>,
    mut state: ResMut<State>,
    mut auto: ResMut<AutoSweep>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    existing: Query<Entity, With<Surface>>,
    mut last: Local<Option<(usize, u32, bool)>>,
) {
    let Some(gpu) = gpu else {
        return;
    };
    let (name, field) = field_at(flags.field);
    let gpu_field = gpu_field_at(flags.field);

    if keyboard.just_pressed(KeyCode::KeyS) || auto.0 {
        auto.0 = false;
        state.sweep.clear();
        // Warm up before timing anything. The first extraction in a process
        // absorbs shader compilation, the first submit and driver
        // initialisation -- measured at **10.76 ms** landing on one
        // `counts_readback`, which is 20x the steady-state figure for the same
        // call. A sweep that starts cold reports that as the smallest grid
        // being the slowest, which is exactly backwards (M-145).
        for _ in 0..2 {
            measure(
                &gpu.0,
                &gpu.1,
                gpu_field,
                &device,
                &queue,
                field.as_ref(),
                SWEEP[0],
            );
        }
        for samples in SWEEP {
            if let Some(point) = measure_median(
                &gpu.0,
                &gpu.1,
                gpu_field,
                &device,
                &queue,
                field.as_ref(),
                samples,
            ) {
                info!(
                    "{name} {:>3}^3: {:>7} tris | cpu {:>8.2} ms | gpu total {:>8.2} = upload {:>7.2} + count {:>6.2} + scan {:>6.2} + emit {:>6.2} + geom-rb {:>7.2} | compute-only {:>6.2} | readback {:>4.0}% | field-on-gpu {:>7.2} ms",
                    point.samples,
                    point.triangles,
                    point.cpu_ms,
                    point.gpu_total_ms(),
                    point.upload_ms,
                    point.gpu.count_ms,
                    point.gpu.scan_ms,
                    point.gpu.emit_ms,
                    point.gpu.geometry_readback_ms,
                    point.gpu_compute_ms(),
                    100.0 * point.gpu.readback_share(),
                    point.gpu_field_total_ms().unwrap_or(f64::NAN),
                );
                state.sweep.push(point);
            }
        }
        // Committed alongside the other measurement CSVs, so the table can be
        // read without running a GPU -- and so a later run can be diffed
        // against it rather than remembered.
        write_sweep_csv(name, &state.sweep);
    }

    let want = (flags.field, state.samples, state.show_cpu);
    if *last == Some(want) && !flags.remesh_requested {
        return;
    }
    *last = Some(want);

    let Some(point) = measure(
        &gpu.0,
        &gpu.1,
        gpu_field,
        &device,
        &queue,
        field.as_ref(),
        state.samples,
    ) else {
        return;
    };
    state.field = name;
    state.current = point;

    // Rebuild whichever mesh is on show. Both are extracted every time
    // regardless -- the timings are the product, and only measuring the one
    // being drawn would make the comparison depend on which key was last
    // pressed.
    let cell = 4.0f32 / (state.samples - 1) as f32;
    let Ok(grid) = GridParams::new([state.samples; 3], [-2.0; 3], cell) else {
        return;
    };
    let mesh = if state.show_cpu {
        let mut cpu = MeshBuffer::<f32>::new();
        let Ok(shape) = RuntimeShape3::new([state.samples; 3]) else {
            return;
        };
        if MarchingCubes::<f32>::new()
            .extract(&field.as_ref(), &shape, [-2.0; 3], cell, &mut cpu)
            .is_err()
        {
            return;
        }
        bevy_isomesh::to_bevy_mesh(&cpu)
    } else {
        let Ok(buffer) = FieldBuffer::sampled(device.wgpu_device(), &queue, grid, &field.as_ref())
        else {
            return;
        };
        let Ok(gpu_mesh) = gpu.0.extract(device.wgpu_device(), &queue, &buffer) else {
            return;
        };
        let mut mesh = Mesh::new(
            bevy::render::mesh::PrimitiveTopology::TriangleList,
            bevy::asset::RenderAssetUsages::default(),
        );
        mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, gpu_mesh.positions.clone());
        mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, gpu_mesh.normals.clone());
        let indices: Vec<u32> = (0..gpu_mesh.positions.len() as u32).collect();
        mesh.insert_indices(bevy::render::mesh::Indices::U32(indices));
        mesh
    };

    for entity in &existing {
        commands.entity(entity).despawn();
    }
    commands.spawn((
        Mesh3d(meshes.add(mesh)),
        MeshMaterial3d(materials.add(StandardMaterial {
            base_color: if state.show_cpu {
                Color::srgb(0.55, 0.62, 0.72)
            } else {
                Color::srgb(0.72, 0.62, 0.50)
            },
            perceptual_roughness: 0.8,
            ..default()
        })),
        Surface,
    ));
}

/// Write the sweep to `docs/measurements/gpu_vs_cpu.csv`.
///
/// Failure is reported and not fatal: a demo that refuses to run because a
/// directory is missing is worse than one that says so and carries on. The
/// numbers are already on screen and in the log either way.
fn write_sweep_csv(field: &str, sweep: &[Point]) {
    let mut out = String::from(
        "field,samples,triangles,cpu_ms,upload_ms,count_ms,scan_ms,emit_ms,geometry_readback_ms,gpu_total_ms,gpu_compute_ms,gpu_field_total_ms\n",
    );
    for p in sweep {
        out.push_str(&format!(
            "{field},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
            p.samples,
            p.triangles,
            p.cpu_ms,
            p.upload_ms,
            p.gpu.count_ms,
            p.gpu.scan_ms,
            p.gpu.emit_ms,
            p.gpu.geometry_readback_ms,
            p.gpu_total_ms(),
            p.gpu_compute_ms(),
            p.gpu_field_total_ms().unwrap_or(f64::NAN),
        ));
    }
    let path = "../docs/measurements/gpu_vs_cpu.csv";
    match std::fs::write(path, out) {
        Ok(()) => info!("wrote {path}"),
        Err(why) => warn!("could not write {path}: {why}"),
    }
}

fn report(state: Res<State>, mut stats: ResMut<DemoStats>) {
    let p = &state.current;
    stats.title = format!(
        "E-302 gpu vs cpu - {} at {}^3 ({})",
        state.field,
        p.samples,
        if state.show_cpu {
            "showing cpu"
        } else {
            "showing gpu"
        }
    );
    stats.triangles = p.triangles;
    // Three per triangle: this path emits a soup, not an indexed mesh.
    stats.vertices = p.triangles * 3;
    stats.extract_ms = p.gpu_total_ms();

    let bar = |ms: f64| {
        let total = p.gpu_total_ms().max(1e-9);
        let width = ((ms / total) * 40.0).round() as usize;
        "#".repeat(width.min(40))
    };
    let ratio = if p.cpu_ms > 0.0 {
        p.gpu_total_ms() / p.cpu_ms
    } else {
        0.0
    };
    let compute_ratio = if p.cpu_ms > 0.0 {
        p.gpu_compute_ms() / p.cpu_ms
    } else {
        0.0
    };

    let mut extra = vec![
        format!(
            "cpu  single-threaded, samples + extract   {:>8.2} ms",
            p.cpu_ms
        ),
        String::new(),
        format!(
            "gpu  everything, triangles on the CPU     {:>8.2} ms   {ratio:>5.2}x cpu",
            p.gpu_total_ms()
        ),
        format!(
            "     upload (cpu samples the field!)     {:>8.2} ms  {}",
            p.upload_ms,
            bar(p.upload_ms)
        ),
        format!(
            "     count pass                          {:>8.2} ms  {}",
            p.gpu.count_ms,
            bar(p.gpu.count_ms)
        ),
        format!(
            "     prefix scan (gpu) + 4-byte total    {:>8.2} ms  {}",
            p.gpu.scan_ms,
            bar(p.gpu.scan_ms)
        ),
        format!(
            "     emit pass                           {:>8.2} ms  {}",
            p.gpu.emit_ms,
            bar(p.gpu.emit_ms)
        ),
        format!(
            "     geometry read-back                  {:>8.2} ms  {}",
            p.gpu.geometry_readback_ms,
            bar(p.gpu.geometry_readback_ms)
        ),
        String::new(),
        format!(
            "read-back is {:>4.0}% of the gpu path",
            100.0 * p.gpu.readback_share()
        ),
        format!(
            "compute only (mesh stays on the gpu)     {:>8.2} ms   {compute_ratio:>5.2}x cpu",
            p.gpu_compute_ms()
        ),
        match p.gpu_field_total_ms() {
            Some(ms) => format!(
                "field evaluated on the gpu, end to end   {:>8.2} ms   {:>5.2}x cpu   (no upload at all)",
                ms,
                if p.cpu_ms > 0.0 { ms / p.cpu_ms } else { 0.0 }
            ),
            None => "field evaluated on the gpu: not implemented for this field".to_string(),
        },
    ];

    if !state.sweep.is_empty() {
        extra.push(String::new());
        extra.push("sweep    tris     cpu ms   gpu ms   ratio   compute-only   ratio".to_string());
        for point in &state.sweep {
            let r = if point.cpu_ms > 0.0 {
                point.gpu_total_ms() / point.cpu_ms
            } else {
                0.0
            };
            let cr = if point.cpu_ms > 0.0 {
                point.gpu_compute_ms() / point.cpu_ms
            } else {
                0.0
            };
            extra.push(format!(
                "{:>4}^3 {:>7}  {:>8.2} {:>8.2}  {:>5.2}x   {:>8.2} ms  {:>5.2}x",
                point.samples,
                point.triangles,
                point.cpu_ms,
                point.gpu_total_ms(),
                r,
                point.gpu_compute_ms(),
                cr
            ));
        }
    }

    extra.push(String::new());
    extra.push("[S] sweep   [ ] resolution   [1-5] field   [G]/[C] gpu/cpu mesh".to_string());
    stats.extra = extra;
}
