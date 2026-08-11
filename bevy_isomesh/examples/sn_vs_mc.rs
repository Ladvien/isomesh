//! E-103 — Surface Nets beside Marching Cubes, on the same field.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example sn_vs_mc --release
//! ```
//!
//! Left is Marching Cubes, right is Surface Nets. `1`–`4` switch field, `[` and
//! `]` change resolution, `S` cycles smoothing passes, `W` shows the wireframe.
//!
//! # What this is actually for
//!
//! The examples catalog bills this as "SN next to MC, triangle counts", on the
//! expectation that the counts differ substantially. **They do not**, and that
//! is the more interesting result, so the HUD puts it front and centre:
//!
//! ```text
//! F_sn − F_mc = 2χ        exactly, on every closed field
//! ```
//!
//! Marching Cubes places one vertex per crossed grid edge; Surface Nets emits
//! two triangles per crossed grid edge; and every closed triangulated surface
//! obeys `F = 2V − 2χ`. The counts are pinned to each other by the number of
//! crossed edges. Surface Nets is not the cheaper method by triangle count, and
//! the HUD recomputes the identity live so you can watch it hold as resolution
//! and field change.
//!
//! What *does* differ is visible rather than numeric, and the wireframe is where
//! to look:
//!
//! - **Corners, and a trap.** Neither method reaches the box corner, and
//!   Marching Cubes lands *further* from it than Surface Nets does — 1.41 cells
//!   against 1.15 at 25³. The reason is that `box_exact` is exactly zero across
//!   its whole boundary, so a grid plane lying on a box face classifies entirely
//!   as outside and the sign change moves a cell inward. Over this domain, 25
//!   and 33 samples are grid-aligned and 27 is not; press `[`/`]` and watch the
//!   corner snap in and out. Closing that gap properly is dual contouring's job,
//!   and E-104 is where it happens — on a box that is *not* grid-aligned, since
//!   benchmarking sharp features on an aligned one measures the classification
//!   rule rather than the algorithm.
//! - **Connectivity.** Marching Cubes produces irregular triangle fans; Surface
//!   Nets produces quads, each split into two triangles along a shared diagonal.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{BoxExact, ReferenceField, Sphere, Torus};
use isomesh::mc::MarchingCubes;
use isomesh::sn::SurfaceNets;
use isomesh::validate::{MeshReport, ValidateConfig, validate_indexed};
use isomesh::{RuntimeShape3, Sdf};

#[derive(Resource)]
struct Demo {
    field: usize,
    samples: u32,
    smoothing: u32,
}

#[derive(Resource)]
struct Materials {
    marching_cubes: Handle<StandardMaterial>,
    surface_nets: Handle<StandardMaterial>,
}

/// Which side of the screen a mesh sits on.
#[derive(Component, Clone, Copy, PartialEq, Eq)]
enum Side {
    MarchingCubes,
    SurfaceNets,
}

const FIELDS: [&str; 4] = ["sphere", "box_exact", "torus", "csg_difference"];
const MIN_SAMPLES: u32 = 5;
const MAX_SAMPLES: u32 = 97;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-103 surface nets vs marching cubes".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Demo {
            field: 0,
            samples: 25,
            smoothing: 0,
        })
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, remesh))
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    for mut orbit in &mut camera {
        // Look straight down -z at the pair, so neither side is nearer than the
        // other and perspective does not make one look bigger than it is. A
        // side-by-side comparison viewed from an angle is not a comparison.
        orbit.yaw = std::f32::consts::FRAC_PI_2;
        orbit.pitch = 0.10;
        orbit.radius = 11.0;
    }
    commands.insert_resource(Materials {
        marching_cubes: materials.add(StandardMaterial {
            base_color: Color::srgb(0.78, 0.80, 0.86),
            perceptual_roughness: 0.45,
            ..default()
        }),
        surface_nets: materials.add(StandardMaterial {
            base_color: Color::srgb(0.86, 0.74, 0.55),
            perceptual_roughness: 0.45,
            ..default()
        }),
    });
}

fn controls(keys: Res<ButtonInput<KeyCode>>, flags: Res<ViewFlags>, mut demo: ResMut<Demo>) {
    if keys.just_pressed(KeyCode::BracketRight) {
        demo.samples = (demo.samples + 2).min(MAX_SAMPLES);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        demo.samples = demo.samples.saturating_sub(2).max(MIN_SAMPLES);
    }
    if keys.just_pressed(KeyCode::KeyS) {
        demo.smoothing = (demo.smoothing + 1) % 5;
    }
    if flags.field < FIELDS.len() {
        demo.field = flags.field;
    }
}

/// One algorithm's output plus what the validity harness says about it.
struct Extraction {
    builder: MeshBuilder,
    report: MeshReport,
    millis: f64,
}

fn extract_pair<F: Sdf<Scalar = f32> + ReferenceField>(
    field: &F,
    samples: u32,
    smoothing: u32,
) -> (Extraction, Extraction, f32) {
    let (min, max) = field.domain();
    let cell_size = (max[0] - min[0]) / (samples - 1) as f32;
    let shape = RuntimeShape3::new([samples; 3]);
    let cfg = ValidateConfig::from_cell_size(f64::from(cell_size));

    let mut mc_builder = MeshBuilder::new();
    let started = Instant::now();
    MarchingCubes::<f32>::new().extract(field, &shape, min, cell_size, &mut mc_builder);
    let mc_millis = started.elapsed().as_secs_f64() * 1000.0;
    let mc_report = validate_indexed(mc_builder.positions(), mc_builder.indices(), &cfg);

    let mut sn = SurfaceNets::<f32>::new();
    sn.set_smoothing_passes(smoothing);
    let mut sn_builder = MeshBuilder::new();
    let started = Instant::now();
    sn.extract(field, &shape, min, cell_size, &mut sn_builder);
    let sn_millis = started.elapsed().as_secs_f64() * 1000.0;
    let sn_report = validate_indexed(sn_builder.positions(), sn_builder.indices(), &cfg);

    (
        Extraction {
            builder: mc_builder,
            report: mc_report,
            millis: mc_millis,
        },
        Extraction {
            builder: sn_builder,
            report: sn_report,
            millis: sn_millis,
        },
        max[0] - min[0],
    )
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    demo: Res<Demo>,
    mut stats: ResMut<DemoStats>,
    mut meshes: ResMut<Assets<Mesh>>,
    materials: Res<Materials>,
    mut commands: Commands,
    mut query: Query<(&mut Mesh3d, &mut Transform, &Side)>,
    mut last: Local<Option<(usize, u32, u32)>>,
) {
    let key = (demo.field, demo.samples, demo.smoothing);
    if *last == Some(key) {
        return;
    }
    *last = Some(key);

    let (mc, sn, width) = match demo.field {
        1 => extract_pair(&BoxExact::<f32>::canonical(), demo.samples, demo.smoothing),
        2 => extract_pair(&Torus::<f32>::canonical(), demo.samples, demo.smoothing),
        3 => extract_pair(
            &isomesh::fields::csg_difference::<f32>(),
            demo.samples,
            demo.smoothing,
        ),
        _ => extract_pair(&Sphere::<f32>::canonical(), demo.samples, demo.smoothing),
    };

    let difference = sn.builder.triangle_count() as i64 - mc.builder.triangle_count() as i64;
    let two_chi = 2 * sn.report.euler_characteristic;

    stats.title = format!(
        "E-103  marching cubes (grey, left)  vs  surface nets (tan, right)\n       field: {}   [1-4] to switch",
        FIELDS[demo.field.min(FIELDS.len() - 1)]
    );
    stats.vertices = mc.builder.vertex_count() + sn.builder.vertex_count();
    stats.triangles = mc.builder.triangle_count() + sn.builder.triangle_count();
    stats.extract_ms = mc.millis + sn.millis;
    stats.extra = vec![
        format!("{:>9} samples/axis   [ and ] to change", demo.samples),
        format!(
            "{:>9} smoothing passes on SN   [S] to cycle",
            demo.smoothing
        ),
        String::new(),
        format!("            {:>10}  {:>10}", "march.cubes", "surf.nets"),
        format!(
            "triangles   {:>10}  {:>10}",
            mc.builder.triangle_count(),
            sn.builder.triangle_count()
        ),
        format!(
            "vertices    {:>10}  {:>10}",
            mc.builder.vertex_count(),
            sn.builder.vertex_count()
        ),
        format!(
            "chi         {:>10}  {:>10}",
            mc.report.euler_characteristic, sn.report.euler_characteristic
        ),
        format!(
            "non-manif.  {:>10}  {:>10}",
            mc.report.non_manifold_edges, sn.report.non_manifold_edges
        ),
        format!("extract ms  {:>10.3}  {:>10.3}", mc.millis, sn.millis),
        String::new(),
        format!(
            "F_sn - F_mc = {difference:<5}   2*chi = {two_chi:<5}   {}",
            if difference == two_chi {
                "EQUAL  <- pinned by Euler, not a coincidence"
            } else {
                "!! differs"
            }
        ),
    ];

    // Side by side, each shifted half a domain apart.
    let offset = width * 0.62;
    let pairs = [
        (Side::MarchingCubes, mc.builder, -offset),
        (Side::SurfaceNets, sn.builder, offset),
    ];

    if query.is_empty() {
        for (side, builder, x) in pairs {
            let material = match side {
                Side::MarchingCubes => materials.marching_cubes.clone(),
                Side::SurfaceNets => materials.surface_nets.clone(),
            };
            commands.spawn((
                Mesh3d(meshes.add(builder.into_mesh())),
                MeshMaterial3d(material),
                Transform::from_xyz(x, 0.0, 0.0),
                DemoMesh,
                side,
            ));
        }
    } else {
        let mut built: Vec<(Side, Handle<Mesh>, f32)> = Vec::new();
        for (side, builder, x) in pairs {
            built.push((side, meshes.add(builder.into_mesh()), x));
        }
        for (mut mesh, mut transform, side) in &mut query {
            if let Some((_, handle, x)) = built.iter().find(|(s, _, _)| s == side) {
                mesh.0 = handle.clone();
                transform.translation.x = *x;
            }
        }
    }
}
