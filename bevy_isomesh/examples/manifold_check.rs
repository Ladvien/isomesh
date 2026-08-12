//! E-111 — where a mesh stops being a manifold, drawn on the mesh.
//!
//! Every reference field, both extractors, with the offending topology drawn in
//! place: **non-manifold edges as thick red lines, non-manifold vertices as red
//! spheres, boundary edges in amber.**
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example manifold_check --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 20-50x slower.
//!
//! `1`-`7` switch field, `A` switches algorithm, `B` toggles the boundary
//! overlay. The rest of the keys are the shared ones — `W` wireframe, `[` `]`
//! resolution, `F12` screenshot.
//!
//! # What this is for
//!
//! The counts have been in the HUD since E-101, and a count tells you a mesh is
//! broken without telling you *where*. Two findings in this project were about
//! **where**, and neither is visible in a number:
//!
//! - Surface Nets goes non-manifold wherever two sheets of surface share a cell.
//!   On the capped gyroid that is 48 edges, scattered through the tunnels. M-15
//!   then found it on a plain convex body, because the real condition is
//!   resolution, not topology.
//! - Marching Cubes does it too, which nobody expected — see ✗15. It needs the
//!   surface to pinch inside a single cell, so it is rare and it is *local*, and
//!   a single red edge on an otherwise clean sphere-union is the entire evidence.
//!
//! Both are one red mark somewhere on a mesh with thousands of triangles.
//!
//! # Why the marks come from the validator
//!
//! The overlay is drawn from [`isomesh::validate::validate_features`], which
//! returns the offending edges and vertices from the *same pass* that produces
//! the counts beside them in the HUD. Recomputing "which edges look wrong" here
//! would let the picture and the caption drift apart, and there would be no way
//! to tell which of the two was lying.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::surface_nets::SurfaceNets;
use isomesh::validate::{ValidateConfig, validate_features};
use isomesh::{RuntimeShape3, Sdf};

/// Defect gizmos get their own config group so they can be thick and drawn in
/// front without dragging the shared wireframe along with them.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct DefectGizmos;

/// Samples per axis.
#[derive(Resource)]
struct Resolution(u32);

const MIN_SAMPLES: u32 = 5;
const MAX_SAMPLES: u32 = 97;

/// Which extractor is running.
#[derive(Resource, Clone, Copy, PartialEq, Eq)]
enum Algorithm {
    MarchingCubes,
    SurfaceNets,
}

impl Algorithm {
    fn name(self) -> &'static str {
        match self {
            Self::MarchingCubes => "marching cubes",
            Self::SurfaceNets => "surface nets",
        }
    }

    /// `ISOMESH_ALGORITHM=mc|sn`, so a capture needs no keyboard — the same
    /// reason `ISOMESH_FIELD` and `ISOMESH_VIEW` exist. Defaults to surface
    /// nets, which is the one with something to show.
    fn from_env() -> Self {
        match std::env::var("ISOMESH_ALGORITHM")
            .unwrap_or_default()
            .as_str()
        {
            "mc" | "marching_cubes" => Self::MarchingCubes,
            _ => Self::SurfaceNets,
        }
    }

    fn toggled(self) -> Self {
        match self {
            Self::MarchingCubes => Self::SurfaceNets,
            Self::SurfaceNets => Self::MarchingCubes,
        }
    }
}

/// Show the boundary overlay. Off by default: on an open field every edge of the
/// clip boundary is a boundary edge, and thousands of amber lines bury the two
/// red ones that matter.
#[derive(Resource)]
struct ShowBoundary(bool);

/// World-space geometry of everything the validator objected to.
///
/// Resolved to positions once per re-mesh rather than per frame — the draw
/// system runs every frame and the mesh changes rarely.
#[derive(Resource, Default)]
struct Overlay {
    non_manifold_edges: Vec<[Vec3; 2]>,
    non_manifold_vertices: Vec<Vec3>,
    inconsistent_edges: Vec<[Vec3; 2]>,
    boundary_edges: Vec<[Vec3; 2]>,
    /// Longest cell dimension, for sizing the vertex spheres.
    cell_size: f32,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-111 manifold check".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<DefectGizmos>()
        .insert_resource(Resolution(33))
        .insert_resource(Algorithm::from_env())
        .insert_resource(ShowBoundary(false))
        .init_resource::<Overlay>()
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, remesh, draw_defects))
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut gizmo_config: ResMut<GizmoConfigStore>,
) {
    // Thick, and biased towards the camera. A non-manifold edge lies exactly on
    // the surface, so at the default bias it z-fights with the triangles sharing
    // it and flickers in and out — which is indistinguishable from the defect
    // being intermittent.
    let (config, _) = gizmo_config.config_mut::<DefectGizmos>();
    config.line.width = 5.0;
    config.depth_bias = -0.4;

    commands.insert_resource(SurfaceMaterial(common::surface_material(&mut materials)));
}

#[derive(Resource)]
struct SurfaceMaterial(Handle<StandardMaterial>);

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut resolution: ResMut<Resolution>,
    mut algorithm: ResMut<Algorithm>,
    mut boundary: ResMut<ShowBoundary>,
) {
    if capture.is_active() {
        // Captured sequences sweep resolution in step with frames, so the GIF is
        // reproducible rather than wall-clock dependent.
        const LOW: u32 = 17;
        const HIGH: u32 = 49;
        let steps = (HIGH - LOW) / 2 + 1;
        let phase = capture.taken % (steps * 2);
        let step = if phase < steps {
            phase
        } else {
            steps * 2 - phase - 1
        };
        resolution.0 = LOW + step * 2;
        return;
    }

    if keys.just_pressed(KeyCode::BracketRight) {
        resolution.0 = (resolution.0 + 2).min(MAX_SAMPLES);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        resolution.0 = resolution.0.saturating_sub(2).max(MIN_SAMPLES);
    }
    if keys.just_pressed(KeyCode::KeyA) {
        *algorithm = algorithm.toggled();
    }
    if keys.just_pressed(KeyCode::KeyB) {
        boundary.0 = !boundary.0;
    }
}

/// Extract, validate, and rebuild the overlay.
///
/// The `for_each_reference_field!` sweep would be the obvious way to reach all
/// seven fields, but it expands to seven concrete blocks and this needs to pick
/// *one* at runtime — so the field index is matched here instead, which is the
/// same shape `surface_nets_vs_marching_cubes` uses.
#[allow(clippy::too_many_arguments)]
fn remesh(
    resolution: Res<Resolution>,
    algorithm: Res<Algorithm>,
    flags: Res<ViewFlags>,
    mut stats: ResMut<DemoStats>,
    mut overlay: ResMut<Overlay>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut query: Query<&mut Mesh3d, With<DemoMesh>>,
    mut commands: Commands,
    material: Res<SurfaceMaterial>,
    mut domain: Query<&mut DemoDomain>,
    mut camera: Query<&mut OrbitCamera>,
    mut last: Local<Option<(u32, usize, Algorithm)>>,
) {
    let key = (resolution.0, flags.field, *algorithm);
    if *last == Some(key) && !flags.remesh_requested {
        return;
    }
    let field_changed = last.map(|(_, f, _)| f) != Some(flags.field) || last.is_none();
    *last = Some(key);

    let field = flags.field.min(FIELD_COUNT - 1);
    let Some(built) = build(field, resolution.0, *algorithm) else {
        return;
    };

    for mut d in &mut domain {
        d.min = Vec3::from(built.domain_min);
        d.max = Vec3::from(built.domain_max);
    }

    // Frame the field rather than assuming one size. The seven domains differ by
    // 4x -- the compact five are half-extent 2, the capped gyroid is 7 and
    // fbm_terrain is 8 -- so a fixed radius puts the camera comfortably inside
    // the gyroid and the demo shows the inner wall of a tunnel.
    if field_changed {
        let half = built.domain_max[0] - built.domain_min[0];
        for mut orbit in &mut camera {
            orbit.radius = half * 1.6;
        }
    }

    *overlay = built.overlay;
    stats.title = format!(
        "E-111  manifold check - {}   field {} ({})   [1-7] field, A algorithm, B boundary",
        algorithm.name(),
        flags.field + 1,
        built.field_name,
    );
    stats.vertices = built.vertices;
    stats.triangles = built.triangles;
    stats.extract_ms = built.extract_ms;
    stats.extra = built.lines;

    let handle = meshes.add(built.mesh);
    if query.is_empty() {
        commands.spawn((Mesh3d(handle), MeshMaterial3d(material.0.clone()), DemoMesh));
    } else {
        for mut mesh in &mut query {
            mesh.0 = handle.clone();
        }
    }
}

const FIELD_COUNT: usize = 7;

struct Built {
    mesh: Mesh,
    overlay: Overlay,
    lines: Vec<String>,
    field_name: &'static str,
    domain_min: [f32; 3],
    domain_max: [f32; 3],
    vertices: usize,
    triangles: usize,
    extract_ms: f64,
}

/// Dispatch on the field index, then do the work once in [`extract_and_check`].
fn build(field: usize, samples: u32, algorithm: Algorithm) -> Option<Built> {
    use isomesh::fields::{
        BoxExact, FbmTerrain, Sphere, ThinPlate, Torus, capped_gyroid, csg_difference,
    };
    match field {
        0 => extract_and_check(&Sphere::<f32>::canonical(), samples, algorithm),
        1 => extract_and_check(&Torus::<f32>::canonical(), samples, algorithm),
        2 => extract_and_check(&BoxExact::<f32>::canonical(), samples, algorithm),
        3 => extract_and_check(&csg_difference::<f32>(), samples, algorithm),
        4 => extract_and_check(&ThinPlate::<f32>::canonical(), samples, algorithm),
        5 => extract_and_check(&capped_gyroid::<f32>(), samples, algorithm),
        _ => extract_and_check(&FbmTerrain::<f32>::canonical(), samples, algorithm),
    }
}

fn extract_and_check<F>(field: &F, samples: u32, algorithm: Algorithm) -> Option<Built>
where
    F: Sdf<Scalar = f32> + ReferenceField,
{
    let (min, max) = field.domain();
    let cell_size = (max[0] - min[0]) / (samples - 1) as f32;
    let shape = match RuntimeShape3::new([samples; 3]) {
        Ok(shape) => shape,
        Err(error) => {
            error!("grid {samples}^3 rejected: {error}");
            return None;
        }
    };

    let mut builder = MeshBuilder::new();
    let started = Instant::now();
    let extracted = match algorithm {
        Algorithm::MarchingCubes => {
            MarchingCubes::<f32>::new().extract(field, &shape, min, cell_size, &mut builder)
        }
        Algorithm::SurfaceNets => {
            SurfaceNets::<f32>::new().extract(field, &shape, min, cell_size, &mut builder)
        }
    };
    if let Err(error) = extracted {
        error!("extraction failed at {samples}^3: {error}");
        return None;
    }
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;

    let cfg = match ValidateConfig::from_cell_size(f64::from(cell_size)) {
        Ok(cfg) => cfg,
        Err(error) => {
            error!("cell size {cell_size} is not a usable spacing: {error}");
            return None;
        }
    };
    let (report, features) = validate_features(builder.positions(), builder.indices(), &cfg);

    // The overlay is the validator's own output, resolved to positions. Nothing
    // here decides what counts as a defect.
    let positions = builder.positions();
    let at = |i: u32| Vec3::from(positions[i as usize]);
    let edges = |list: &[[u32; 2]]| -> Vec<[Vec3; 2]> {
        list.iter().map(|e| [at(e[0]), at(e[1])]).collect()
    };
    let overlay = Overlay {
        non_manifold_edges: edges(&features.edges),
        non_manifold_vertices: features.vertices.iter().map(|v| at(*v)).collect(),
        inconsistent_edges: edges(&features.inconsistently_oriented_edges),
        boundary_edges: edges(&features.boundary_edges),
        cell_size,
    };

    let verdict = if report.is_closed() {
        "MANIFOLD, CLOSED"
    } else if report.is_manifold() {
        "MANIFOLD, WITH BOUNDARY"
    } else {
        "!! NON-MANIFOLD - see the red marks"
    };

    let lines = vec![
        format!("{:>9} samples/axis   [ and ] to change", samples),
        format!("{:>9.4} cell size", cell_size),
        String::new(),
        format!(
            "{:>9} non-manifold edges      (red lines)",
            report.non_manifold_edges
        ),
        format!(
            "{:>9} non-manifold vertices   (red spheres)",
            report.non_manifold_vertices
        ),
        format!(
            "{:>9} inconsistently oriented (magenta)",
            report.inconsistently_oriented_edges
        ),
        format!(
            "{:>9} boundary edges          (amber, B to show)",
            report.boundary_edges
        ),
        String::new(),
        format!("{:>9} euler characteristic", report.euler_characteristic),
        format!("{:>9} components", report.components),
        format!(
            "{:>9} closed in domain (from the field, not a guess)",
            field.closed_in_domain()
        ),
        format!("          {verdict}"),
    ];

    Some(Built {
        overlay,
        lines,
        field_name: F::NAME,
        domain_min: min,
        domain_max: max,
        vertices: builder.vertex_count(),
        triangles: builder.triangle_count(),
        extract_ms,
        mesh: builder.into_mesh(),
    })
}

/// Draw the overlay. Runs every frame; the overlay itself changes only on
/// re-mesh.
fn draw_defects(
    overlay: Res<Overlay>,
    boundary: Res<ShowBoundary>,
    mut gizmos: Gizmos<DefectGizmos>,
) {
    const RED: Color = Color::srgb(1.0, 0.13, 0.13);
    const MAGENTA: Color = Color::srgb(1.0, 0.0, 0.85);
    const AMBER: Color = Color::srgb(1.0, 0.65, 0.1);

    if boundary.0 {
        for [a, b] in &overlay.boundary_edges {
            gizmos.line(*a, *b, AMBER);
        }
    }
    for [a, b] in &overlay.inconsistent_edges {
        gizmos.line(*a, *b, MAGENTA);
    }
    // Drawn after the others so a defect is never hidden under the boundary.
    for [a, b] in &overlay.non_manifold_edges {
        gizmos.line(*a, *b, RED);
    }
    for v in &overlay.non_manifold_vertices {
        gizmos
            .sphere(
                Isometry3d::from_translation(*v),
                overlay.cell_size * 0.22,
                RED,
            )
            .resolution(8);
    }
}
