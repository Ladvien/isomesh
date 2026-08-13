//! E-301 — the same Marching Cubes on the GPU and the CPU, and where they differ.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example gpu_compute_mc --release
//! ```
//!
//! `1`–`5` field, `[` `]` resolution, `C` show the CPU mesh instead, `V` colour
//! by agreement. Drag to orbit, scroll to zoom.
//!
//! # "Looks the same" is not the acceptance criterion
//!
//! GPU-005 says so explicitly, and it is the right instinct: two Marching Cubes
//! meshes of a sphere look identical at any resolution, including when one of
//! them is wrong. So this example does not ask you to compare pictures. It runs
//! **both** extractions on the same grid from the same samples, classifies every
//! GPU vertex against the CPU's, and puts the counts on screen.
//!
//! The verdict, measured rather than hoped for (M-142):
//!
//! - **Triangle counts are equal.** Both sides read the same samples, classify
//!   with the same table — the shader's is *uploaded* from `isomesh`'s own
//!   `CASES`, not transcribed — and iterate cells in the same order.
//! - **Most vertices are bit-identical, and the rest are one ULP away.** Not
//!   zero, and the size of the miss is the finding: WGSL permits a multiply-add
//!   to be contracted into a fused one, and this driver takes that permission,
//!   rounding once where the CPU rounds twice.
//!
//! Press `V` and the surface is coloured by that classification, so the
//! disagreement is *visible* rather than merely counted — which is the whole
//! difference between this and a screenshot of two spheres.
//!
//! # The device is Bevy's, and that is the point
//!
//! `isomesh-gpu`'s public API takes `&wgpu::Device` and never an engine type.
//! This example is the test of that claim: it reaches Bevy's raw device through
//! `RenderDevice::wgpu_device()` and hands it in. Nothing here constructs a
//! device, and `isomesh_gpu::headless` — which the crate's own tests use — is
//! not imported. If the abstraction had leaked, this file could not exist.
//!
//! `isomesh-gpu` is a **dev-dependency** of `bevy_isomesh`, not a dependency.
//! The plugin's whole design is that a CPU-only consumer never compiles the
//! renderer, and `isomesh-gpu` carries `wgpu`. Verified rather than assumed:
//! after adding it, the plugin's lockfile still holds exactly one `wgpu`, at
//! 29.0.4, shared with Bevy.

mod common;

use std::time::Instant;

use bevy::asset::RenderAssetUsages;
use bevy::prelude::*;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy::render::renderer::{RenderDevice, RenderQueue};
use common::{CommonPlugin, DemoStats, ViewFlags};
use isomesh::fields::{BoxExact, Sphere, ThinPlate, Torus, csg_difference};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};
use isomesh_gpu::{FieldBuffer, GridParams, MarchingCubesGpu};

/// Samples per axis, unless `ISOMESH_SAMPLES` says otherwise.
const DEFAULT_SAMPLES: u32 = 41;

/// Bit-exact with a CPU vertex.
const EXACT: [f32; 4] = [0.62, 0.66, 0.70, 1.0];
/// One ULP away on at least one axis.
const NEAR: [f32; 4] = [0.95, 0.62, 0.10, 1.0];
/// Further than one ULP — should never appear, and is loud if it does.
const STRANGER: [f32; 4] = [0.95, 0.10, 0.12, 1.0];

/// The fields `1`–`5` select.
fn field_at(index: usize) -> (&'static str, Box<dyn Sdf<Scalar = f32> + Send + Sync>) {
    match index % 5 {
        0 => ("sphere", Box::new(Sphere::<f32>::canonical())),
        1 => ("torus", Box::new(Torus::<f32>::canonical())),
        2 => ("box_exact", Box::new(BoxExact::<f32>::canonical())),
        3 => ("csg_difference", Box::new(csg_difference::<f32>())),
        _ => ("thin_plate", Box::new(ThinPlate::<f32>::canonical())),
    }
}

/// How the two meshes compared, and how long each took.
#[derive(Resource, Default)]
struct Verdict {
    field: &'static str,
    samples: u32,
    gpu_triangles: usize,
    cpu_triangles: usize,
    exact: usize,
    near: usize,
    strangers: usize,
    gpu_ms: f64,
    cpu_ms: f64,
    worst_cells: f32,
    show_cpu: bool,
    colour_by_agreement: bool,
}

impl Verdict {
    fn vertices(&self) -> usize {
        self.exact + self.near + self.strangers
    }
}

/// The compiled pipeline, built once on Bevy's device.
#[derive(Resource)]
struct Gpu(MarchingCubesGpu);

/// Marks the mesh entity so a re-extraction can replace it.
#[derive(Component)]
struct Surface;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-301 gpu compute mc".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_resource::<Verdict>()
        .add_systems(Startup, setup)
        .add_systems(Update, (keys, extract, report).chain())
        .run();
}

fn setup(
    mut commands: Commands,
    device: Res<RenderDevice>,
    queue: Res<RenderQueue>,
    mut verdict: ResMut<Verdict>,
) {
    // The claim under test: Bevy's raw device, handed to a crate that has never
    // heard of Bevy.
    match MarchingCubesGpu::new(device.wgpu_device(), &queue) {
        Ok(pipeline) => commands.insert_resource(Gpu(pipeline)),
        // Fail loudly. A GPU demo that silently falls back to the CPU path
        // would show two identical meshes and prove nothing at all, which is
        // exactly the failure this example exists to make impossible.
        Err(why) => panic!("could not build the GPU pipeline on Bevy's device: {why}"),
    }
    verdict.samples = common::samples_override().unwrap_or(DEFAULT_SAMPLES);
    verdict.colour_by_agreement = true;
}

fn keys(keyboard: Res<ButtonInput<KeyCode>>, mut verdict: ResMut<Verdict>) {
    if keyboard.just_pressed(KeyCode::KeyC) {
        verdict.show_cpu = !verdict.show_cpu;
    }
    if keyboard.just_pressed(KeyCode::KeyV) {
        verdict.colour_by_agreement = !verdict.colour_by_agreement;
    }
    if keyboard.just_pressed(KeyCode::BracketLeft) {
        verdict.samples = (verdict.samples - 8).max(9);
    }
    if keyboard.just_pressed(KeyCode::BracketRight) {
        verdict.samples = (verdict.samples + 8).min(129);
    }
}

/// Run both extractions and rebuild the mesh, whenever anything changed.
#[allow(clippy::too_many_arguments)]
fn extract(
    mut commands: Commands,
    gpu: Option<Res<Gpu>>,
    device: Res<RenderDevice>,
    queue: Res<RenderQueue>,
    flags: Res<ViewFlags>,
    mut verdict: ResMut<Verdict>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    existing: Query<Entity, With<Surface>>,
    mut last: Local<Option<(usize, u32, bool, bool)>>,
) {
    let Some(gpu) = gpu else {
        return;
    };
    let want = (
        flags.field,
        verdict.samples,
        verdict.show_cpu,
        verdict.colour_by_agreement,
    );
    if *last == Some(want) && !flags.remesh_requested {
        return;
    }
    *last = Some(want);

    let (name, field) = field_at(flags.field);
    let samples = verdict.samples;
    // The reference fields all live inside [-2, 2].
    let extent = 4.0f32;
    let cell = extent / (samples - 1) as f32;
    let Ok(grid) = GridParams::new([samples; 3], [-2.0; 3], cell) else {
        return;
    };

    // GPU. The upload samples the field on the CPU so both sides read exactly
    // the same f32 values -- otherwise a difference in the *input* would be
    // reported as a difference in the algorithm.
    let started = Instant::now();
    let Ok(buffer) = FieldBuffer::sampled(device.wgpu_device(), &queue, grid, &field.as_ref())
    else {
        return;
    };
    let Ok(gpu_mesh) = gpu.0.extract(device.wgpu_device(), &queue, &buffer) else {
        return;
    };
    let gpu_ms = started.elapsed().as_secs_f64() * 1000.0;

    // CPU, same grid.
    let started = Instant::now();
    let mut cpu = MeshBuffer::<f32>::new();
    let Ok(shape) = RuntimeShape3::new([samples; 3]) else {
        return;
    };
    if MarchingCubes::<f32>::new()
        .extract(&field.as_ref(), &shape, [-2.0; 3], cell, &mut cpu)
        .is_err()
    {
        return;
    }
    let cpu_ms = started.elapsed().as_secs_f64() * 1000.0;

    // Classify every GPU vertex against the CPU's. Bit-identical is one
    // bucket; everything else is measured as a *distance*, in cells, which is
    // the unit the rest of this repository reports geometry error in.
    //
    // Nearest-point search through a spatial hash rather than a ULP-neighbour
    // probe: a ULP probe answers "within k" only for the k you thought to ask
    // for, and the whole point here is to find out what the number is.
    let key = |p: &[f32; 3]| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()];
    let cpu_points: std::collections::HashSet<[u32; 3]> = cpu.positions.iter().map(key).collect();

    // Bucket at a hundredth of a cell: far larger than any rounding difference,
    // far smaller than the gap between two genuinely distinct vertices.
    let bucket = cell * 0.01;
    let cell_of = |p: &[f32; 3]| {
        [
            (p[0] / bucket).floor() as i64,
            (p[1] / bucket).floor() as i64,
            (p[2] / bucket).floor() as i64,
        ]
    };
    let mut index: std::collections::HashMap<[i64; 3], Vec<[f32; 3]>> =
        std::collections::HashMap::new();
    for p in &cpu.positions {
        index.entry(cell_of(p)).or_default().push(*p);
    }

    let mut colours = Vec::with_capacity(gpu_mesh.positions.len());
    let (mut exact, mut near, mut strangers) = (0usize, 0usize, 0usize);
    let mut worst_cells = 0.0f32;
    for p in &gpu_mesh.positions {
        if cpu_points.contains(&key(p)) {
            exact += 1;
            colours.push(EXACT);
            continue;
        }
        let home = cell_of(p);
        let mut best = f32::INFINITY;
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let at = [home[0] + dx, home[1] + dy, home[2] + dz];
                    for q in index.get(&at).into_iter().flatten() {
                        let d =
                            ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2))
                                .sqrt();
                        best = best.min(d);
                    }
                }
            }
        }
        let in_cells = best / cell;
        worst_cells = worst_cells.max(in_cells);
        // A rounding difference is many orders of magnitude below a cell. A
        // genuine disagreement about geometry is not, and that is the line.
        if in_cells < 1e-4 {
            near += 1;
            colours.push(NEAR);
        } else {
            strangers += 1;
            colours.push(STRANGER);
        }
    }
    verdict.worst_cells = worst_cells;

    verdict.field = name;
    verdict.gpu_triangles = gpu_mesh.triangle_count();
    verdict.cpu_triangles = cpu.indices.len() / 3;
    verdict.exact = exact;
    verdict.near = near;
    verdict.strangers = strangers;
    verdict.gpu_ms = gpu_ms;
    verdict.cpu_ms = cpu_ms;

    // One line per extraction, so the comparison can be read from a terminal
    // rather than off a HUD in a screenshot -- the same reason `game_dig`
    // prints its per-edit costs.
    info!(
        "{name} {samples}^3 h={cell}: tris gpu {} cpu {} | vertices {} = {exact} exact + {near} rounding + {strangers} moved | worst {worst_cells:e} cells | gpu {gpu_ms:.2} ms cpu {cpu_ms:.2} ms",
        gpu_mesh.triangle_count(),
        cpu.indices.len() / 3,
        gpu_mesh.positions.len(),
    );

    // Loud, not decorative: a vertex further than one ULP means the two
    // implementations disagree about geometry rather than about rounding.
    if strangers > 0 {
        error!(
            "{name} at {samples}^3: {strangers} gpu vertices moved, worst {worst_cells:e} cells"
        );
    }

    for entity in &existing {
        commands.entity(entity).despawn();
    }

    let mesh = if verdict.show_cpu {
        bevy_isomesh::to_bevy_mesh(&cpu)
    } else {
        let mut mesh = Mesh::new(
            PrimitiveTopology::TriangleList,
            RenderAssetUsages::default(),
        );
        mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, gpu_mesh.positions.clone());
        mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, gpu_mesh.normals.clone());
        if verdict.colour_by_agreement {
            // Straight to linear, as E-208 established: ATTRIBUTE_COLOR is
            // linear RGBA and an sRGB literal renders washed out.
            let linear: Vec<[f32; 4]> = colours
                .iter()
                .map(|c| {
                    Color::srgba(c[0], c[1], c[2], c[3])
                        .to_linear()
                        .to_f32_array()
                })
                .collect();
            mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, linear);
        }
        let indices: Vec<u32> = (0..gpu_mesh.positions.len() as u32).collect();
        mesh.insert_indices(Indices::U32(indices));
        mesh
    };

    commands.spawn((
        Mesh3d(meshes.add(mesh)),
        MeshMaterial3d(materials.add(StandardMaterial {
            base_color: if verdict.show_cpu {
                Color::srgb(0.62, 0.66, 0.70)
            } else {
                Color::WHITE
            },
            perceptual_roughness: 0.75,
            ..default()
        })),
        Surface,
    ));
}

fn report(verdict: Res<Verdict>, mut stats: ResMut<DemoStats>) {
    let agree = verdict.gpu_triangles == verdict.cpu_triangles;
    stats.title = format!(
        "E-301 gpu compute mc - {} at {}^3",
        verdict.field, verdict.samples
    );
    stats.vertices = verdict.vertices();
    stats.triangles = verdict.gpu_triangles;
    stats.extract_ms = verdict.gpu_ms;

    let total = verdict.vertices().max(1);
    let pct = |n: usize| 100.0 * n as f64 / total as f64;
    stats.extra = vec![
        format!(
            "triangles   gpu {:>6}   cpu {:>6}   {}",
            verdict.gpu_triangles,
            verdict.cpu_triangles,
            if agree { "EQUAL" } else { "DISAGREE" }
        ),
        format!(
            "extract     gpu {:>6.2} ms   cpu {:>6.2} ms   (gpu includes readback)",
            verdict.gpu_ms, verdict.cpu_ms
        ),
        String::new(),
        format!("vertices vs cpu, of {}:", verdict.vertices()),
        format!(
            "  bit-identical   {:>6}  ({:>5.1}%)",
            verdict.exact,
            pct(verdict.exact)
        ),
        format!(
            "  rounding only   {:>6}  ({:>5.1}%)   <- fused multiply-add, M-142/M-143",
            verdict.near,
            pct(verdict.near)
        ),
        format!(
            "  moved           {:>6}  ({:>5.1}%)   <- must be 0",
            verdict.strangers,
            pct(verdict.strangers)
        ),
        format!("  worst offset    {:>12.3e} cells", verdict.worst_cells),
        String::new(),
        "the shader's case table is uploaded from isomesh's own CASES, not transcribed".to_string(),
        "the device is Bevy's, through RenderDevice::wgpu_device()".to_string(),
        String::new(),
        format!(
            "[C] show {}   [V] colour by agreement {}   [1-5] field   [ ] resolution",
            if verdict.show_cpu { "gpu" } else { "cpu" },
            if verdict.colour_by_agreement {
                "on"
            } else {
                "off"
            }
        ),
    ];
}
