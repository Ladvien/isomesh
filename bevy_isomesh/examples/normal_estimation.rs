//! E-113 — three normal-estimation strategies, lit, side by side.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example normal_estimation --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 20-50x slower.
//!
//! `1`-`5` field, `[` `]` resolution, `W` wireframe, `F12` screenshot.
//!
//! # Why this example is lit and the others are not
//!
//! **The three meshes are geometrically identical.** Same positions, same
//! indices, same triangle count — `normals::recompute` overwrites the normal
//! buffer and touches nothing else. A wireframe view of this comparison shows
//! three identical pictures, which is exactly why the wireframe is off by default
//! and the material is glossy.
//!
//! A normal is only ever observed through shading. So the honest way to show what
//! a normal-estimation strategy costs is to light it and look at the speculars,
//! and the honest way to *quantify* it is the angle readout in the HUD — which is
//! the same number `area_weighted_normals_track_the_field_on_smooth_geometry_and_not_on_sharp`
//! asserts in the test suite.
//!
//! # What to look for
//!
//! - **`sphere`** — all three panels look the same, and the numbers say so:
//!   central differences sit under a tenth of a degree from the analytic gradient
//!   at ordinary resolutions.
//! - **`box_exact`** — the right panel's corners and edges are visibly rounded in
//!   the shading while the silhouette stays sharp. That is M-66: the area-weighted
//!   normal at a corner is the average of three face normals, and the field's
//!   gradient there is one of them.
//! - **Press `]` on `box_exact`.** The mean disagreement falls with resolution and
//!   the **worst does not** — 35.796 degrees at 17³, 33³ and 65³ alike. Refining a
//!   grid does not soften a corner.

mod common;

use bevy::prelude::*;
use bevy_isomesh::to_bevy_mesh;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{BoxExact, ReferenceField, Sphere, ThinPlate, Torus, csg_difference};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::normals::{NormalStrategy, recompute};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

#[derive(Resource)]
struct Demo {
    field: usize,
    samples: u32,
}

const FIELDS: [&str; 5] = [
    "sphere",
    "box_exact",
    "torus",
    "csg_difference",
    "thin_plate",
];
const MIN_SAMPLES: u32 = 5;
const MAX_SAMPLES: u32 = 97;

/// Which panel a mesh is, left to right.
#[derive(Component, Clone, Copy, PartialEq, Eq)]
enum Panel {
    Analytic,
    CentralDifference,
    AreaWeighted,
}

impl Panel {
    const ALL: [Self; 3] = [Self::Analytic, Self::CentralDifference, Self::AreaWeighted];

    fn label(self) -> &'static str {
        match self {
            Self::Analytic => "analytic gradient",
            Self::CentralDifference => "central diff (step = h)",
            Self::AreaWeighted => "area-weighted faces",
        }
    }
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-113 normal estimation".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Demo {
            field: 3,
            samples: common::samples_override().unwrap_or(25),
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
        // Straight on, so no panel is nearer than another and perspective does
        // not make one look bigger than it is.
        orbit.yaw = std::f32::consts::FRAC_PI_2;
        orbit.pitch = 0.30;
        orbit.radius = 9.5;
    }

    // Glossy on purpose. A matte surface hides exactly the difference this
    // example exists to show — a normal is only ever observed through shading,
    // and a rough material integrates the very variation being compared.
    commands.insert_resource(PanelMaterial(materials.add(StandardMaterial {
        // Warm clay, and glossy. The colour is not decoration: a cold grey reads
        // as a diagnostic dump, and this image has to make someone *want* to look
        // at the difference before it can show them one.
        base_color: Color::srgb(0.87, 0.72, 0.56),
        perceptual_roughness: 0.20,
        metallic: 0.02,
        reflectance: 0.35,
        ..default()
    })));
}

#[derive(Resource)]
struct PanelMaterial(Handle<StandardMaterial>);

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

/// One panel's mesh and how far its normals sit from the analytic answer.
struct Panelled {
    buffer: MeshBuffer<f32>,
    worst: f32,
    mean: f32,
}

fn build<F: Sdf<Scalar = f32> + ReferenceField>(
    field: &F,
    samples: u32,
) -> Option<(Vec<(Panel, Panelled)>, f32)> {
    let (min, max) = field.domain();
    let cell_size = (max[0] - min[0]) / (samples - 1) as f32;
    let shape = RuntimeShape3::new([samples; 3]).ok()?;

    let mut base = MeshBuffer::<f32>::new();
    MarchingCubes::<f32>::new()
        .extract(field, &shape, min, cell_size, &mut base)
        .ok()?;
    if base.is_empty() {
        return None;
    }

    // The analytic normals are the reference every panel is measured against, so
    // take them once before anything overwrites them.
    let reference = base.normals.clone();

    let mut out = Vec::new();
    for panel in Panel::ALL {
        let strategy = match panel {
            Panel::Analytic => NormalStrategy::AnalyticGradient,
            Panel::CentralDifference => NormalStrategy::CentralDifference { step: cell_size },
            Panel::AreaWeighted => NormalStrategy::AreaWeightedFaces,
        };
        let mut buffer = base.clone();
        // A degenerate normal is a real answer, not a reason to draw nothing —
        // report it and skip the panel rather than substituting something.
        if let Err(error) = recompute(&mut buffer, field, strategy) {
            error!("{}: {error}", panel.label());
            return None;
        }

        let (worst, mean) = deviation(&reference, &buffer.normals);
        out.push((
            panel,
            Panelled {
                buffer,
                worst,
                mean,
            },
        ));
    }
    Some((out, max[0] - min[0]))
}

/// Worst and mean angle between two normal sets, in degrees.
fn deviation(a: &[[f32; 3]], b: &[[f32; 3]]) -> (f32, f32) {
    let mut worst = 0.0f32;
    let mut total = 0.0f32;
    for (x, y) in a.iter().zip(b) {
        let dot = (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]).clamp(-1.0, 1.0);
        let angle = dot.acos().to_degrees();
        worst = worst.max(angle);
        total += angle;
    }
    let count = a.len().max(1) as f32;
    (worst, total / count)
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    demo: Res<Demo>,
    mut stats: ResMut<DemoStats>,
    mut meshes: ResMut<Assets<Mesh>>,
    material: Res<PanelMaterial>,
    mut commands: Commands,
    mut query: Query<(&mut Mesh3d, &mut Transform, &Panel)>,
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
        0 => build(&Sphere::<f32>::canonical(), demo.samples),
        1 => build(&BoxExact::<f32>::canonical(), demo.samples),
        2 => build(&Torus::<f32>::canonical(), demo.samples),
        3 => build(&csg_difference::<f32>(), demo.samples),
        _ => build(&ThinPlate::<f32>::canonical(), demo.samples),
    };
    let Some((panels, width)) = built else {
        return;
    };

    let vertices = panels.first().map_or(0, |(_, p)| p.buffer.vertex_count());
    let triangles = panels.first().map_or(0, |(_, p)| p.buffer.triangle_count());

    stats.title = format!(
        "E-113  normal estimation   field {} ({})   [1-5] field, [ ] resolution",
        demo.field + 1,
        FIELDS[demo.field]
    );
    stats.vertices = vertices;
    stats.triangles = triangles;
    let mut lines = vec![
        format!("{} samples/axis   [ and ] to change", demo.samples),
        "geometry is identical in all three panels; only the shading differs".into(),
        String::new(),
        format!("{:<26} {:>9} {:>9}", "strategy", "worst", "mean"),
    ];
    for (panel, p) in &panels {
        lines.push(format!(
            "{:<26} {:>8.3}d {:>8.3}d",
            panel.label(),
            p.worst,
            p.mean
        ));
    }
    lines.push(String::new());
    lines.push("degrees from the analytic gradient, left panel to each".into());
    stats.extra = lines;

    // Tight enough that the three fill the frame. A comparison the reader has to
    // squint at is not a comparison — and the difference this example is about
    // lives in a specular highlight a few pixels across.
    let offset = width * 0.60;
    for mut orbit in &mut camera {
        orbit.radius = width * 2.35;
    }

    let placed: Vec<(Panel, Handle<Mesh>, f32)> = panels
        .iter()
        .enumerate()
        .map(|(i, (panel, p))| {
            let x = (i as f32 - 1.0) * offset;
            (*panel, meshes.add(to_bevy_mesh(&p.buffer)), x)
        })
        .collect();

    if query.is_empty() {
        for (panel, handle, x) in placed {
            commands.spawn((
                Mesh3d(handle),
                MeshMaterial3d(material.0.clone()),
                Transform::from_xyz(x, 0.0, 0.0),
                DemoMesh,
                panel,
            ));
        }
    } else {
        for (mut mesh, mut transform, panel) in &mut query {
            if let Some((_, handle, x)) = placed.iter().find(|(p, _, _)| p == panel) {
                mesh.0 = handle.clone();
                transform.translation.x = *x;
            }
        }
    }
}
