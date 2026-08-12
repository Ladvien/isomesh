//! E-106 — the blocky path: face culling against greedy merging.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example greedy_quads --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 20-50x slower.
//!
//! `1`-`6` field, `[` `]` resolution, `W` wireframe (**on by default here**),
//! `F12` screenshot.
//!
//! # What this shows
//!
//! Both panels are the same occupancy grid and the same visible faces. The left
//! emits one quad per face; the right merges coplanar runs. So the wireframe is
//! the whole demo — the surfaces are identical and only the *quads* differ, which
//! is why this example starts wireframed like E-105 does.
//!
//! # The number in the HUD is the finding
//!
//! Greedy meshing is usually quoted as **`2.76x` fewer triangles than face
//! culling**, from one UE5 benchmark. A-005 measured it across seven fields at
//! the same resolution and it is not a constant (M-56):
//!
//! | field | saving |
//! |---|---|
//! | `gyroid` | 1.70x |
//! | `sphere` | 1.94x |
//! | `torus` | 2.69x |
//! | `fbm_terrain` | 4.60x |
//! | `csg_difference` | 10.64x |
//! | `box_exact` | **256x** |
//!
//! Merging pays for **flat runs**, so a grid-aligned box collapses to six quads
//! at every resolution — twelve triangles at 17³, 33³ and 65³ alike — while a
//! sphere's staircase surface is short runs and barely merges. The published
//! figure happens to land beside `torus`.
//!
//! Press `2` for `box_exact` and then `]` a few times to watch the right panel
//! stay at twelve triangles while the left one grows with the grid.
//!
//! # Two honest limitations on display
//!
//! - **`thin_plate` comes back empty.** It is 0.4 cells thick, and this algorithm
//!   asks one question per cell — is its *centre* inside — so a feature thinner
//!   than a cell does not exist to it. That is not a bug to fix; it is the
//!   premise, and it is what A-014's subgrid work exists to lift.
//! - **The mesh is open, on purpose.** A cube corner needs three normals, so
//!   vertices are split and the index buffer describes a surface with boundary.
//!   Welding closes it wherever quads meet corner to corner and **cannot** close
//!   a T-junction, where a long quad butts against several short ones (M-57).

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{
    BoxExact, FbmTerrain, ReferenceField, Sphere, ThinPlate, Torus, csg_difference,
};
use isomesh::greedy_quads::{GreedyQuads, Merge};
use isomesh::{RuntimeShape3, Sdf};

#[derive(Resource)]
struct Demo {
    field: usize,
    samples: u32,
}

const FIELDS: [&str; 6] = [
    "sphere",
    "box_exact",
    "torus",
    "csg_difference",
    "fbm_terrain",
    "thin_plate (empty — thinner than a cell)",
];
const MIN_SAMPLES: u32 = 5;
const MAX_SAMPLES: u32 = 65;

#[derive(Component, Clone, Copy, PartialEq, Eq)]
enum Side {
    Culled,
    Merged,
}

#[derive(Resource)]
struct Materials {
    culled: Handle<StandardMaterial>,
    merged: Handle<StandardMaterial>,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-106 greedy quads".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Demo {
            field: 4,
            samples: common::samples_override().unwrap_or(33),
        })
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, remesh))
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
    mut flags: ResMut<ViewFlags>,
) {
    // Both panels are the same surface. Only the quads differ, and quads are
    // only visible in the wireframe.
    flags.wireframe = true;

    for mut orbit in &mut camera {
        orbit.yaw = std::f32::consts::FRAC_PI_2 + 0.5;
        orbit.pitch = 0.38;
        orbit.radius = 14.0;
    }
    commands.insert_resource(Materials {
        culled: materials.add(StandardMaterial {
            base_color: Color::srgb(0.78, 0.70, 0.60),
            perceptual_roughness: 0.72,
            ..default()
        }),
        merged: materials.add(StandardMaterial {
            base_color: Color::srgb(0.88, 0.66, 0.42),
            perceptual_roughness: 0.62,
            ..default()
        }),
    });
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    mut demo: ResMut<Demo>,
    capture: Res<Capture>,
    mut flags: ResMut<ViewFlags>,
) {
    if capture.is_active() {
        return;
    }
    for (key, index) in [
        (KeyCode::Digit1, 0),
        (KeyCode::Digit2, 1),
        (KeyCode::Digit3, 2),
        (KeyCode::Digit4, 3),
        (KeyCode::Digit5, 4),
        (KeyCode::Digit6, 5),
    ] {
        if keys.just_pressed(key) {
            demo.field = index;
        }
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        demo.samples = (demo.samples + 4).min(MAX_SAMPLES);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        demo.samples = demo.samples.saturating_sub(4).max(MIN_SAMPLES);
    }
    if keys.just_pressed(KeyCode::KeyR) {
        flags.remesh_requested = true;
    }
}

struct Extraction {
    builder: MeshBuilder,
    millis: f64,
}

fn extract_one<F: Sdf<Scalar = f32> + ReferenceField>(
    field: &F,
    shape: &RuntimeShape3,
    origin: [f32; 3],
    cell_size: f32,
    merge: Merge,
) -> Option<Extraction> {
    let mut builder = MeshBuilder::new();
    let mut mesher = GreedyQuads::<f32>::new();
    mesher.set_merge(merge);
    let started = Instant::now();
    mesher.extract(field, shape, origin, cell_size, &mut builder).ok()?;
    Some(Extraction {
        builder,
        millis: started.elapsed().as_secs_f64() * 1000.0,
    })
}

fn extract_pair<F: Sdf<Scalar = f32> + ReferenceField>(
    field: &F,
    samples: u32,
) -> Option<(Extraction, Extraction, f32)> {
    let (min, max) = field.domain();
    let cell_size = (max[0] - min[0]) / (samples - 1) as f32;
    let shape = RuntimeShape3::new([samples; 3]).ok()?;
    let culled = extract_one(field, &shape, min, cell_size, Merge::Off)?;
    let merged = extract_one(field, &shape, min, cell_size, Merge::Greedy)?;
    Some((culled, merged, max[0] - min[0]))
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    demo: Res<Demo>,
    mut stats: ResMut<DemoStats>,
    mut meshes: ResMut<Assets<Mesh>>,
    materials: Res<Materials>,
    mut commands: Commands,
    mut query: Query<(&mut Mesh3d, &mut Transform, &Side)>,
    mut camera: Query<&mut OrbitCamera>,
    mut flags: ResMut<ViewFlags>,
    mut last: Local<Option<(usize, u32)>>,
) {
    let key = (demo.field, demo.samples);
    if *last == Some(key) && !flags.remesh_requested {
        return;
    }
    *last = Some(key);
    flags.remesh_requested = false;

    let built = match demo.field {
        0 => extract_pair(&Sphere::<f32>::canonical(), demo.samples),
        1 => extract_pair(&BoxExact::<f32>::canonical(), demo.samples),
        2 => extract_pair(&Torus::<f32>::canonical(), demo.samples),
        3 => extract_pair(&csg_difference::<f32>(), demo.samples),
        4 => extract_pair(&FbmTerrain::<f32>::canonical(), demo.samples),
        _ => extract_pair(&ThinPlate::<f32>::canonical(), demo.samples),
    };
    let Some((culled, merged, width)) = built else {
        return;
    };

    let saving = if merged.builder.triangle_count() == 0 {
        0.0
    } else {
        culled.builder.triangle_count() as f64 / merged.builder.triangle_count() as f64
    };

    stats.title = format!(
        "E-106  face culling | greedy merge   field {} ({})   [1-6] field, [ ] resolution",
        demo.field + 1,
        FIELDS[demo.field]
    );
    stats.vertices = culled.builder.vertex_count() + merged.builder.vertex_count();
    stats.triangles = culled.builder.triangle_count() + merged.builder.triangle_count();

    let mut lines = vec![
        format!("{} samples/axis   [ and ] to change", demo.samples),
        String::new(),
        format!("{:<12} {:>12} {:>12}", "", "culled", "merged"),
        format!(
            "{:<12} {:>12} {:>12}",
            "quads",
            culled.builder.triangle_count() / 2,
            merged.builder.triangle_count() / 2
        ),
        format!(
            "{:<12} {:>12} {:>12}",
            "triangles",
            culled.builder.triangle_count(),
            merged.builder.triangle_count()
        ),
        format!(
            "{:<12} {:>12.3} {:>12.3}",
            "extract ms", culled.millis, merged.millis
        ),
        String::new(),
        format!("merging saves {saving:.2}x here"),
    ];
    if merged.builder.triangle_count() == 0 {
        lines.push("this field is thinner than one cell, and this algorithm asks".into());
        lines.push("one question per cell -- is its centre inside. so: nothing.".into());
    } else {
        lines.push("the published figure is 2.76x, from one scene. measured across".into());
        lines.push("seven fields it ranges 1.70x to 256x (M-56) -- merging pays for".into());
        lines.push("flat runs, so a grid-aligned box collapses to six quads at any".into());
        lines.push("resolution and a sphere's staircase barely merges at all.".into());
    }
    stats.extra = lines;

    // A blocky field fills its whole domain, so the two panels need more room
    // between them than a sphere does — and more distance, or they run off the
    // bottom of the frame.
    let offset = width * 0.66;
    for mut orbit in &mut camera {
        orbit.radius = width * 3.35;
    }

    let pairs = [
        (Side::Culled, culled.builder, -offset),
        (Side::Merged, merged.builder, offset),
    ];

    if query.is_empty() {
        for (side, builder, x) in pairs {
            let material = match side {
                Side::Culled => materials.culled.clone(),
                Side::Merged => materials.merged.clone(),
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
        let built: Vec<(Side, Handle<Mesh>, f32)> = pairs
            .into_iter()
            .map(|(side, builder, x)| (side, meshes.add(builder.into_mesh()), x))
            .collect();
        for (mut mesh, mut transform, side) in &mut query {
            if let Some((_, handle, x)) = built.iter().find(|(s, _, _)| s == side) {
                mesh.0 = handle.clone();
                transform.translation.x = *x;
            }
        }
    }
}
