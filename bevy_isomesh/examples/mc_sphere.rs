//! E-101 — Marching Cubes, the baseline.
//!
//! A sphere meshed with [`isomesh::mc::MarchingCubes`], with a wireframe toggle
//! and a live resolution slider on the `[` and `]` keys. The first thing in this
//! project you can actually look at.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example mc_sphere --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 20-50x slower and will convince
//! you something is wrong with the algorithm.
//!
//! The HUD reports what the validity harness says about the mesh on screen —
//! Euler characteristic and manifoldness, recomputed on every re-mesh. That is
//! the point of the example: not that a sphere appears, but that the sphere is
//! demonstrably a closed 2-manifold with `χ = 2`.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{ReferenceField, Sphere};
use isomesh::mc::MarchingCubes;
use isomesh::validate::{ValidateConfig, validate_indexed};
use isomesh::{RuntimeShape3, Shape3};

/// Samples per axis. `[` and `]` step it.
#[derive(Resource)]
struct Resolution(u32);

const MIN_SAMPLES: u32 = 5;
const MAX_SAMPLES: u32 = 129;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-101 marching cubes".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Resolution(33))
        .add_systems(Startup, setup)
        .add_systems(Update, (change_resolution, remesh))
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    let (min, max) = Sphere::<f32>::canonical().domain();
    commands.spawn(DemoDomain {
        min: Vec3::from(min),
        max: Vec3::from(max),
    });

    for mut orbit in &mut camera {
        orbit.radius = 6.5;
    }

    // No placeholder mesh. A `Mesh` with no attributes and no indices makes
    // Bevy's slab allocator report a use-after-free, because there is nothing to
    // allocate for it. The entity is spawned by the first `remesh`, which runs
    // on the first update.
    commands.insert_resource(SurfaceMaterial(common::surface_material(&mut materials)));
}

/// Held so `remesh` can spawn the mesh entity the first time it runs.
#[derive(Resource)]
struct SurfaceMaterial(Handle<StandardMaterial>);

fn change_resolution(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut resolution: ResMut<Resolution>,
) {
    // While recording, sweep the resolution in step with the captured frames
    // rather than with wall-clock time, so the sequence is reproducible.
    if capture.is_active() {
        const LOW: u32 = 9;
        const HIGH: u32 = 81;
        let steps = (HIGH - LOW) / 2 + 1;
        let phase = capture.taken % (steps * 2);
        let step = if phase < steps { phase } else { steps * 2 - phase - 1 };
        resolution.0 = LOW + step * 2;
        return;
    }

    // A cell count that is a multiple of 4 puts grid samples exactly on the unit
    // sphere, which is a real and visible effect -- it produces zero-area
    // slivers. Stepping by 2 samples walks through both kinds so you can see it.
    if keys.just_pressed(KeyCode::BracketRight) {
        resolution.0 = (resolution.0 + 2).min(MAX_SAMPLES);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        resolution.0 = resolution.0.saturating_sub(2).max(MIN_SAMPLES);
    }
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    resolution: Res<Resolution>,
    flags: Res<ViewFlags>,
    mut stats: ResMut<DemoStats>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut query: Query<&mut Mesh3d, With<DemoMesh>>,
    mut commands: Commands,
    material: Res<SurfaceMaterial>,
    mut mc: Local<MarchingCubes<f32>>,
    mut last: Local<u32>,
) {
    let changed = resolution.0 != *last;
    if !changed && !flags.remesh_requested {
        return;
    }
    *last = resolution.0;

    let field = Sphere::<f32>::canonical();
    let (min, max) = field.domain();
    let samples = resolution.0;
    let cell_size = (max[0] - min[0]) / (samples - 1) as f32;
    // Every fallible call is handled rather than unwrapped: a demo that aborts
    // on a bad parameter teaches nothing about the parameter.
    let shape = match RuntimeShape3::new([samples; 3]) {
        Ok(shape) => shape,
        Err(error) => {
            error!("grid {samples}^3 rejected: {error}");
            return;
        }
    };

    let mut builder = MeshBuilder::new();
    let started = Instant::now();
    if let Err(error) = mc.extract(&field, &shape, min, cell_size, &mut builder) {
        error!("extraction failed at {samples}^3: {error}");
        stats.extra = vec![format!("extraction failed: {error}")];
        return;
    }
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;

    // The claim the example exists to make, checked on the mesh being shown
    // rather than asserted in prose.
    let Ok(cfg) = ValidateConfig::from_cell_size(f64::from(cell_size)) else {
        error!("cell size {cell_size} is not a usable spacing");
        return;
    };
    let report = validate_indexed(builder.positions(), builder.indices(), &cfg);

    stats.title = "E-101  marching cubes on a sphere".into();
    stats.vertices = builder.vertex_count();
    stats.triangles = builder.triangle_count();
    stats.extract_ms = extract_ms;
    stats.extra = vec![
        format!("{:>9} samples/axis   [ and ] to change", samples),
        format!(
            "{:>9} cells  ({} total)",
            samples - 1,
            shape.element_count()
        ),
        String::new(),
        format!(
            "{:>9} euler characteristic (2 = closed sphere)",
            report.euler_characteristic
        ),
        format!("{:>9} non-manifold edges", report.non_manifold_edges),
        format!("{:>9} boundary edges", report.boundary_edges),
        format!("{:>9} degenerate triangles", report.degenerate_triangles),
        format!(
            "          {}",
            if report.is_closed() {
                "MANIFOLD, CLOSED"
            } else {
                "!! NOT CLOSED"
            }
        ),
    ];

    let handle = meshes.add(builder.into_mesh());
    if query.is_empty() {
        commands.spawn((Mesh3d(handle), MeshMaterial3d(material.0.clone()), DemoMesh));
    } else {
        for mut mesh in &mut query {
            mesh.0 = handle.clone();
        }
    }
}
