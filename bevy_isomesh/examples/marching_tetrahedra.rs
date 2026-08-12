//! E-105 — Marching Cubes against Marching Tetrahedra, same field, same grid.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example marching_tetrahedra --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 20-50x slower.
//!
//! `1`-`5` field, `[` `]` resolution, `W` wireframe (**on by default here**),
//! `F12` screenshot.
//!
//! # What this shows, and why it is not a recommendation
//!
//! The two surfaces are nearly the same shape. The difference is entirely in how
//! densely they are cut up, so this example starts with the wireframe **on** — a
//! shaded view of these two panels is two pictures of the same object, and the
//! whole finding is in the triangle count.
//!
//! The ratio in the HUD is the interesting number, and A-003 measured that it is
//! not the constant the literature quotes:
//!
//! - The tier-R figure of *"2–3× more vertices"* measures **2.87× to 3.91×** here,
//!   and covers only the two roughest fields (M-51).
//! - The spread is not noise. It is `4.0` when the surface normal lies inside one
//!   octant and `2.0` when it changes sign, and the isotropic average of those is
//!   `2.992` (M-52) — which is why a grid-aligned box reads 3.91 and a sphere,
//!   which samples every octant, reads 3.04.
//! - What the extra triangles buy is **4.3%** better geometry on a smooth field,
//!   and they buy it in the wrong direction from what was predicted: Marching
//!   Tetrahedra is *more* accurate than Marching Cubes on sharp fields, because
//!   its extra edge families sample a corner from more directions (M-55).
//!
//! So the value on show is not "use this algorithm". It is that the crate lets
//! you put two of them on the same grid and read the trade off a HUD.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{BoxExact, ReferenceField, Sphere, ThinPlate, Torus, csg_difference};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_tetrahedra::MarchingTetrahedra;
use isomesh::{RuntimeShape3, Sdf};

#[derive(Resource)]
struct Demo {
    field: usize,
    samples: u32,
}

const FIELDS: [&str; 5] = ["sphere", "box_exact", "torus", "csg_difference", "thin_plate"];
const MIN_SAMPLES: u32 = 5;
const MAX_SAMPLES: u32 = 65;

#[derive(Component, Clone, Copy, PartialEq, Eq)]
enum Side {
    MarchingCubes,
    MarchingTetrahedra,
}

#[derive(Resource)]
struct Materials {
    marching_cubes: Handle<StandardMaterial>,
    marching_tetrahedra: Handle<StandardMaterial>,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-105 marching tetrahedra".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Demo {
            field: 0,
            samples: common::samples_override().unwrap_or(21),
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
    // The difference *is* the tessellation, so start where it is visible. Every
    // other example in this repo starts shaded; this one would be showing two
    // identical pictures if it did.
    flags.wireframe = true;

    for mut orbit in &mut camera {
        orbit.yaw = std::f32::consts::FRAC_PI_2;
        orbit.pitch = 0.20;
        orbit.radius = 11.0;
    }
    commands.insert_resource(Materials {
        // Two warm tones rather than one, so a glance tells you which panel is
        // which without reading the labels.
        marching_cubes: materials.add(StandardMaterial {
            base_color: Color::srgb(0.86, 0.74, 0.58),
            perceptual_roughness: 0.55,
            ..default()
        }),
        marching_tetrahedra: materials.add(StandardMaterial {
            base_color: Color::srgb(0.72, 0.80, 0.86),
            perceptual_roughness: 0.55,
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

fn extract_pair<F: Sdf<Scalar = f32> + ReferenceField>(
    field: &F,
    samples: u32,
) -> Option<(Extraction, Extraction, f32)> {
    let (min, max) = field.domain();
    let cell_size = (max[0] - min[0]) / (samples - 1) as f32;
    let shape = RuntimeShape3::new([samples; 3]).ok()?;

    let mut mc = MeshBuilder::new();
    let started = Instant::now();
    MarchingCubes::<f32>::new()
        .extract(field, &shape, min, cell_size, &mut mc)
        .ok()?;
    let mc_millis = started.elapsed().as_secs_f64() * 1000.0;

    let mut mt = MeshBuilder::new();
    let started = Instant::now();
    MarchingTetrahedra::<f32>::new()
        .extract(field, &shape, min, cell_size, &mut mt)
        .ok()?;
    let mt_millis = started.elapsed().as_secs_f64() * 1000.0;

    Some((
        Extraction {
            builder: mc,
            millis: mc_millis,
        },
        Extraction {
            builder: mt,
            millis: mt_millis,
        },
        max[0] - min[0],
    ))
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
        _ => extract_pair(&ThinPlate::<f32>::canonical(), demo.samples),
    };
    let Some((mc, mt, width)) = built else {
        return;
    };

    let ratio = |a: usize, b: usize| {
        if a == 0 {
            0.0
        } else {
            b as f64 / a as f64
        }
    };
    let vertex_ratio = ratio(mc.builder.vertex_count(), mt.builder.vertex_count());
    let triangle_ratio = ratio(mc.builder.triangle_count(), mt.builder.triangle_count());

    stats.title = format!(
        "E-105  marching cubes | marching tetrahedra   field {} ({})   [1-5] field, [ ] resolution",
        demo.field + 1,
        FIELDS[demo.field]
    );
    stats.vertices = mc.builder.vertex_count() + mt.builder.vertex_count();
    stats.triangles = mc.builder.triangle_count() + mt.builder.triangle_count();
    stats.extra = vec![
        format!("{} samples/axis   [ and ] to change", demo.samples),
        String::new(),
        format!("{:<12} {:>12} {:>12} {:>9}", "", "cubes", "tetrahedra", "ratio"),
        format!(
            "{:<12} {:>12} {:>12} {:>8.3}x",
            "vertices",
            mc.builder.vertex_count(),
            mt.builder.vertex_count(),
            vertex_ratio
        ),
        format!(
            "{:<12} {:>12} {:>12} {:>8.3}x",
            "triangles",
            mc.builder.triangle_count(),
            mt.builder.triangle_count(),
            triangle_ratio
        ),
        format!(
            "{:<12} {:>12.3} {:>12.3}",
            "extract ms", mc.millis, mt.millis
        ),
        String::new(),
        "the ratio is 4.0 where the normal stays in one octant and 2.0 across a".into(),
        "sign change (M-52), so a grid-aligned box reads ~3.9 and a sphere ~3.0.".into(),
        "what it buys is 4.3% better geometry on a sphere -- and *better* than".into(),
        "marching cubes on a sharp field, which is not what was predicted.".into(),
    ];

    // Wide enough that the pair clears the HUD on the left and leaves the frame
    // with a margin. A wireframe comparison is unreadable when a keybinding line
    // is drawn across one of the two meshes.
    let offset = width * 0.58;
    for mut orbit in &mut camera {
        orbit.radius = width * 2.55;
    }

    let pairs = [
        (Side::MarchingCubes, mc.builder, -offset),
        (Side::MarchingTetrahedra, mt.builder, offset),
    ];

    if query.is_empty() {
        for (side, builder, x) in pairs {
            let material = match side {
                Side::MarchingCubes => materials.marching_cubes.clone(),
                Side::MarchingTetrahedra => materials.marching_tetrahedra.clone(),
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
