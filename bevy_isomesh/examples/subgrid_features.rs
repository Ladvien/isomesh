//! E-108 — letters thinner than a voxel, and the extractor that can see them.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example subgrid_features --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 20-50x slower, and this example
//! runs the most expensive extractor in the crate.
//!
//! `[` `]` resolution, `-` `=` letter thickness, `W` wireframe, `F12` screenshot.
//!
//! # What this shows
//!
//! Two panels, one field, one grid. On the left, Marching Cubes. On the right,
//! subgrid Marching Tetrahedra. The letters are extruded to a thickness measured
//! **in voxels**, and `-` drives that below one.
//!
//! Somewhere just under a voxel the left panel stops being a picture of the
//! letters and starts being a picture of the sampling: first a scatter, then
//! nothing at all. The right panel does not change. Same grid, same field, same
//! triangle budget — the difference is only what question each extractor asks of
//! an edge.
//!
//! A sign test asks *what sign is this endpoint* and gets one bit. Subgrid
//! marching asks *where are all the zeros along this edge* and gets a list. M-67
//! puts a number on the gap: a sign test cannot distinguish **95.6%** of the
//! configurations a tetrahedron can actually be in.
//!
//! # The two failure modes on the left are different, and both are on show
//!
//! A-005 measured `thin_plate` returning **zero** triangles from greedy quads,
//! because that method asks one question per cell *centre* and no centre is
//! inside a plate 0.4 cells thick. That is the clean failure: the feature is
//! simply absent.
//!
//! Marching Cubes fails the other way, and M-72 measured it: the feature does
//! not vanish, it **aliases**. It samples corners and cuts edges, so whichever
//! edges happen to straddle the slab still register a sign change, and what
//! comes back is a partial, holey remnant that changes shape with the grid. For
//! a streamed world that is the worse behaviour — a feature that disappears at a
//! known distance can be faded; one that disintegrates into a resolution-
//! dependent scatter pops.
//!
//! Push `-` slowly and both are visible in order: solid, then holey, then gone.
//!
//! Getting that to be *true* took a correction. The sheet is offset off the grid
//! lattice and tilted slightly, because centred on `z = 0` it always contains a
//! whole plane of grid nodes and Marching Cubes never loses it — the left panel
//! reported the same 576 triangles at every thickness (M-100).
//!
//! # What it costs
//!
//! The HUD reports both extraction times, and the ratio is not small — M-98
//! measured **70×** classic Marching Tetrahedra on `sphere`, and the constant is
//! field evaluations rather than anything algorithmic: at 16 samples per edge a
//! cell costs `6 tets × 6 edges × 16` evaluations before refinement, against
//! Marching Cubes' 8 shared corner samples.
//!
//! The comparison worth making is not the one on the HUD. It is "subgrid at this
//! resolution against Marching Cubes at whatever resolution resolves the same
//! feature" — and for a letter 0.4 voxels thick, there is no such resolution.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::subgrid::extract::SubgridMarchingTetrahedra;
use isomesh::{RuntimeShape3, Sdf};

/// Half the extent of the sampling box, in world units.
const HALF_EXTENT: f32 = 1.6;
/// 1D samples per tetrahedron edge. Enough to bracket a tenth of a voxel.
const EDGE_SAMPLES: u32 = 16;

const MIN_SAMPLES: u32 = 9;
const MAX_SAMPLES: u32 = 41;
/// Letter thickness, in voxels. Below 1.0 is where the left panel starts losing.
const MIN_THICKNESS: f32 = 0.10;
const MAX_THICKNESS: f32 = 3.0;

#[derive(Resource)]
struct Demo {
    samples: u32,
    /// Letter thickness in **voxels**, so the interesting number is the one the
    /// user is turning.
    thickness_voxels: f32,
}

#[derive(Component, Clone, Copy, PartialEq, Eq)]
enum Side {
    MarchingCubes,
    Subgrid,
}

#[derive(Resource)]
struct Materials {
    marching_cubes: Handle<StandardMaterial>,
    subgrid: Handle<StandardMaterial>,
}

/// The word `ISO`, as thin extruded strokes.
///
/// Strokes are ~1.6 voxels **wide** and a settable number of voxels **thick**,
/// so only one dimension is sub-voxel. That separation is the whole point: a
/// feature that is small in every direction is just a small feature, and any
/// method loses it. A feature that is large in two directions and thin in the
/// third is a *sheet*, it is what real geometry looks like at a wall or a fin,
/// and it is exactly what a sign test cannot represent.
struct Letters {
    half_thickness: f32,
    /// Where the sheet's mid-plane sits, in world units.
    ///
    /// **Deliberately not zero, and deliberately not a multiple of the cell.**
    /// Centred on `z = 0` the demo is broken in a way that looks fine: `z = 0`
    /// is a grid plane at every odd resolution, so the sheet always contains a
    /// whole plane of nodes, Marching Cubes always finds a sign change on the
    /// vertical edges through them, and the left panel reports the *same 576
    /// triangles* at 0.45, 0.30 and 0.15 voxels. It never loses the letters at
    /// all. See M-100 — the fifth time a fixture in this repo has been placed
    /// exactly where the property under test cannot fail.
    centre: f32,
    /// Slope of the mid-plane in `y`, so the sheet crosses grid planes at a
    /// shallow angle instead of lying between two of them everywhere.
    ///
    /// Without it the failure is binary — every node inside, or none — and the
    /// left panel jumps from whole letters to nothing. The tilt is what makes
    /// M-72's middle phase visible: the sheet is inside the lattice in some
    /// places and between planes in others, so what comes back is the holey
    /// remnant, which is the failure mode a streamed world actually suffers.
    tilt: f32,
}

/// Letter thickness in voxels, requested through `ISOMESH_THICKNESS`.
fn thickness_override() -> Option<f32> {
    std::env::var("ISOMESH_THICKNESS")
        .ok()
        .and_then(|raw| raw.parse::<f32>().ok())
        .map(|value| value.clamp(MIN_THICKNESS, MAX_THICKNESS))
}

/// Distance from `p` to the segment `a`-`b`, in the plane.
fn segment_distance(p: [f32; 2], a: [f32; 2], b: [f32; 2]) -> f32 {
    let pa = [p[0] - a[0], p[1] - a[1]];
    let ba = [b[0] - a[0], b[1] - a[1]];
    let denominator = ba[0] * ba[0] + ba[1] * ba[1];
    let t = if denominator > 0.0 {
        ((pa[0] * ba[0] + pa[1] * ba[1]) / denominator).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let d = [pa[0] - ba[0] * t, pa[1] - ba[1] * t];
    (d[0] * d[0] + d[1] * d[1]).sqrt()
}

/// Every stroke of `ISO`, as 2D segments in a box roughly `[-1.4, 1.4] × [-0.7, 0.7]`.
const STROKES: [[[f32; 2]; 2]; 10] = [
    // I
    [[-0.95, -0.6], [-0.95, 0.6]],
    // S, as five bars -- a seven-segment S, which reads clearly and needs no
    // curve primitive.
    [[-0.5, 0.6], [0.1, 0.6]],
    [[-0.5, 0.6], [-0.5, 0.0]],
    [[-0.5, 0.0], [0.1, 0.0]],
    [[0.1, 0.0], [0.1, -0.6]],
    [[-0.5, -0.6], [0.1, -0.6]],
    // O, as a rectangular ring.
    [[0.5, -0.6], [0.5, 0.6]],
    [[1.1, -0.6], [1.1, 0.6]],
    [[0.5, 0.6], [1.1, 0.6]],
    [[0.5, -0.6], [1.1, -0.6]],
];

/// Half the stroke width, in world units.
const STROKE_HALF_WIDTH: f32 = 0.08;

impl Sdf for Letters {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        let flat = [p[0], p[1]];
        let mut nearest = f32::INFINITY;
        for [a, b] in STROKES {
            nearest = nearest.min(segment_distance(flat, a, b));
        }
        // A 2D stroke, extruded. The `max` is the exact SDF of an extrusion only
        // outside the solid; inside it under-estimates, which is harmless here
        // because every consumer of this field either tests the sign or finds a
        // root, and both only care about where it vanishes.
        let in_plane = nearest - STROKE_HALF_WIDTH;
        let mid_plane = self.centre + self.tilt * p[1];
        let out_of_plane = (p[2] - mid_plane).abs() - self.half_thickness;
        in_plane.max(out_of_plane)
    }
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-108 subgrid features".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Demo {
            samples: common::samples_override().unwrap_or(25),
            // Below the threshold where Marching Cubes still gets it, so the
            // first frame already shows the difference rather than needing a
            // keypress. `ISOMESH_THICKNESS` overrides it, which is what lets a
            // capture pick its frame without a keyboard -- the same reason the
            // harness has `ISOMESH_SAMPLES`.
            thickness_voxels: thickness_override().unwrap_or(0.35),
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
        // Looking slightly down the letters' plane, so a sheet reads as a sheet
        // rather than as a line.
        orbit.yaw = std::f32::consts::FRAC_PI_2;
        orbit.pitch = 0.55;
        orbit.radius = 7.5;
    }
    commands.insert_resource(Materials {
        marching_cubes: materials.add(StandardMaterial {
            base_color: Color::srgb(0.86, 0.74, 0.58),
            perceptual_roughness: 0.55,
            double_sided: true,
            cull_mode: None,
            ..default()
        }),
        subgrid: materials.add(StandardMaterial {
            base_color: Color::srgb(0.72, 0.80, 0.86),
            perceptual_roughness: 0.55,
            double_sided: true,
            cull_mode: None,
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
    if keys.just_pressed(KeyCode::BracketRight) {
        demo.samples = (demo.samples + 4).min(MAX_SAMPLES);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        demo.samples = demo.samples.saturating_sub(4).max(MIN_SAMPLES);
    }
    if keys.just_pressed(KeyCode::Equal) {
        demo.thickness_voxels = (demo.thickness_voxels + 0.15).min(MAX_THICKNESS);
    }
    if keys.just_pressed(KeyCode::Minus) {
        demo.thickness_voxels = (demo.thickness_voxels - 0.15).max(MIN_THICKNESS);
    }
    if keys.just_pressed(KeyCode::KeyR) {
        flags.remesh_requested = true;
    }
}

struct Extraction {
    builder: MeshBuilder,
    millis: f64,
}

fn extract_pair(samples: u32, thickness_voxels: f32) -> Option<(Extraction, Extraction)> {
    let min = [-HALF_EXTENT; 3];
    let cell_size = (2.0 * HALF_EXTENT) / (samples - 1) as f32;
    let shape = RuntimeShape3::new([samples; 3]).ok()?;
    // The thickness the user is turning is in voxels, so it tracks the grid: a
    // resolution change alone must not make the letters thicker or thinner in
    // the units this demo is about.
    let field = Letters {
        half_thickness: 0.5 * thickness_voxels * cell_size,
        // 0.41 of a cell off the lattice: far enough from a node plane that a
        // sheet under ~0.8 voxels thick fits between two of them, and not a
        // round fraction, so no resolution lands it back on one.
        centre: 0.41 * cell_size,
        tilt: 0.03,
    };

    let mut mc = MeshBuilder::new();
    let started = Instant::now();
    MarchingCubes::<f32>::new()
        .extract(&field, &shape, min, cell_size, &mut mc)
        .ok()?;
    let mc_millis = started.elapsed().as_secs_f64() * 1000.0;

    let mut sub = MeshBuilder::new();
    let started = Instant::now();
    SubgridMarchingTetrahedra::<f32>::new(EDGE_SAMPLES)
        .ok()?
        .extract(&field, &shape, min, cell_size, &mut sub)
        .ok()?;
    let sub_millis = started.elapsed().as_secs_f64() * 1000.0;

    Some((
        Extraction {
            builder: mc,
            millis: mc_millis,
        },
        Extraction {
            builder: sub,
            millis: sub_millis,
        },
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
    mut flags: ResMut<ViewFlags>,
    mut last: Local<Option<(u32, u32)>>,
) {
    // Thickness is quantised into the key so a float comparison never re-meshes
    // every frame; the step is 0.15 voxels, so hundredths are plenty.
    let key = (demo.samples, (demo.thickness_voxels * 100.0) as u32);
    if *last == Some(key) && !flags.remesh_requested {
        return;
    }
    *last = Some(key);
    flags.remesh_requested = false;

    let Some((mc, sub)) = extract_pair(demo.samples, demo.thickness_voxels) else {
        return;
    };

    stats.title = format!(
        "E-108  marching cubes | subgrid marching tetrahedra   letters {:.2} voxels thick   \
         [-/= thickness, [ ] resolution]",
        demo.thickness_voxels
    );
    stats.vertices = mc.builder.vertex_count() + sub.builder.vertex_count();
    stats.triangles = mc.builder.triangle_count() + sub.builder.triangle_count();

    let verdict = if mc.builder.triangle_count() == 0 {
        "gone -- no grid edge straddles the sheet"
    } else if sub.builder.triangle_count() > 4 * mc.builder.triangle_count().max(1) {
        "aliasing -- a holey remnant, not the letters (M-72)"
    } else {
        "both resolve it; press - to thin the letters"
    };

    stats.extra = vec![
        format!(
            "{} samples/axis, {EDGE_SAMPLES} per tet edge   [ and ] to change",
            demo.samples
        ),
        String::new(),
        format!("{:<12} {:>12} {:>12}", "", "cubes", "subgrid"),
        format!(
            "{:<12} {:>12} {:>12}",
            "vertices",
            mc.builder.vertex_count(),
            sub.builder.vertex_count()
        ),
        format!(
            "{:<12} {:>12} {:>12}",
            "triangles",
            mc.builder.triangle_count(),
            sub.builder.triangle_count()
        ),
        format!(
            "{:<12} {:>12.2} {:>12.2}",
            "extract ms", mc.millis, sub.millis
        ),
        String::new(),
        format!("left panel: {verdict}"),
        String::new(),
        "a sign test gets one bit per edge; subgrid marching gets every zero".into(),
        "along it. M-67: a sign test cannot tell apart 95.6% of the".into(),
        "configurations a tetrahedron can be in.".into(),
        String::new(),
        "the cost is real -- M-98 measured 70x classic marching tetrahedra --".into(),
        "but the honest comparison is against whatever resolution would".into(),
        "resolve the same feature, and below one voxel there is none.".into(),
    ];

    let offset = HALF_EXTENT * 1.25;
    let pairs = [
        (Side::MarchingCubes, mc.builder, -offset),
        (Side::Subgrid, sub.builder, offset),
    ];

    if query.is_empty() {
        for (side, builder, x) in pairs {
            let material = match side {
                Side::MarchingCubes => materials.marching_cubes.clone(),
                Side::Subgrid => materials.subgrid.clone(),
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
