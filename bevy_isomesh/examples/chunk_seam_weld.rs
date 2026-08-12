//! E-115 — the crack between two chunks, and welding it shut.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example chunk_seam_weld --release
//! ```
//!
//! Two adjacent chunks, meshed **independently** — as a real game does, because
//! an edit only dirties the chunks it touches. `V` welds and unwelds, `E`
//! explodes the two chunks apart so you can see they are separate meshes, `[`
//! and `]` change the chunk resolution, `1`/`2` switch spacing.
//!
//! Red lines are **boundary edges on the seam plane**: mesh edges used by
//! exactly one triangle, where the surface should be continuous. That is the
//! crack. Amber lines are the chunks' own outer boundary, which is not a defect
//! — a chunk is an open patch by construction, and telling those two apart is
//! the whole reason `validate` counts `boundary_edges` separately from
//! `non_manifold_edges`.
//!
//! # Why the spacing selector matters more than it looks
//!
//! `1` is `h = 0.125` and `2` is `h = 4/35`. Both look arbitrary and only one
//! is: **M-32** measured that two chunks agree on their shared sample plane
//! bit-for-bit only when the cell size is a power of two, because one computes
//! `(o + h·cn) + h·n` and the other `o + h·(c+1)n` — equal by algebra, not by
//! IEEE — and 22% of random `(origin, h, cells, chunk)` combinations disagree by
//! an ulp. `4/35` comes from that search.
//!
//! So a weld keyed on exact equality would close the seam under `1` and leave it
//! open under `2`. This crate's weld is an epsilon weld for exactly that reason,
//! and pressing `1`/`2` with `V` on is the demonstration: the seam closes at both
//! spacings, and the HUD shows the vertices being merged at both.
//!
//! # What welding costs, and what it is not
//!
//! Welding is a **post-pass over one buffer**, so the two chunks have to be
//! concatenated first (`MeshBuffer::append`) and the indices shifted. It is not
//! a meshing mode and not something the extractor can do for you: neither chunk
//! can know the other exists.
//!
//! The rule is first fit against the vertices already kept, in input index
//! order, lowest index wins — stated on `isomesh::weld` because epsilon-closeness
//! is not transitive and "weld everything within ε" does not define equivalence
//! classes.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::Torus;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{ValidateConfig, validate_features};
use isomesh::weld::Welder;
use isomesh::{MeshBuffer, MeshSink, Sdf};

/// Gizmos for the seam overlay, so they can be drawn thick and on top.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct SeamGizmos;

#[derive(Resource)]
struct Demo {
    welded: bool,
    exploded: bool,
    cells: u32,
    /// Index into [`SPACINGS`].
    spacing: usize,
}

/// `0.125` is a power of two and bit-exact at a seam; `4/35` is not, and comes
/// from M-32's search over the 22% of combinations that disagree.
const SPACINGS: [(f32, &str); 2] = [(0.125, "0.125 (power of two)"), (4.0 / 35.0, "4/35")];

/// The seam overlay, rebuilt only when the mesh is.
#[derive(Resource, Default)]
struct Overlay {
    seam_edges: Vec<[Vec3; 2]>,
    outer_edges: Vec<[Vec3; 2]>,
    cell_size: f32,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-115 chunk seam weld".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<SeamGizmos>()
        .insert_resource(Demo {
            // `ISOMESH_WELD=1` starts welded, so the before/after pair can be
            // captured without a human pressing `V` -- the same reason the
            // harness has `ISOMESH_VIEW`.
            welded: std::env::var("ISOMESH_WELD").is_ok(),
            exploded: false,
            cells: common::samples_override().unwrap_or(28),
            spacing: 1,
        })
        .init_resource::<Overlay>()
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, remesh, draw_seam))
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
    mut config: ResMut<GizmoConfigStore>,
) {
    for mut orbit in &mut camera {
        // Down the seam plane's normal is useless — the crack is edge-on there.
        // Three quarters round and slightly above shows both chunks and the
        // plane between them.
        orbit.yaw = 0.9;
        orbit.pitch = 0.35;
        orbit.radius = 4.6;
    }
    let (seam, _) = config.config_mut::<SeamGizmos>();
    seam.line.width = 3.0;
    seam.depth_bias = -0.2;

    commands.insert_resource(Materials {
        left: materials.add(StandardMaterial {
            base_color: Color::srgb(0.78, 0.80, 0.86),
            perceptual_roughness: 0.45,
            ..default()
        }),
        right: materials.add(StandardMaterial {
            base_color: Color::srgb(0.86, 0.74, 0.55),
            perceptual_roughness: 0.45,
            ..default()
        }),
    });
}

#[derive(Resource)]
struct Materials {
    left: Handle<StandardMaterial>,
    right: Handle<StandardMaterial>,
}

#[derive(Component, Clone, Copy, PartialEq, Eq)]
enum Half {
    Left,
    Right,
}

fn controls(keys: Res<ButtonInput<KeyCode>>, flags: Res<ViewFlags>, mut demo: ResMut<Demo>) {
    if keys.just_pressed(KeyCode::KeyV) {
        demo.welded = !demo.welded;
    }
    if keys.just_pressed(KeyCode::KeyE) {
        demo.exploded = !demo.exploded;
    }
    if keys.just_pressed(KeyCode::Digit1) {
        demo.spacing = 0;
    }
    if keys.just_pressed(KeyCode::Digit2) {
        demo.spacing = 1;
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        demo.cells = (demo.cells - 4).max(24);
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        demo.cells = (demo.cells + 4).min(44);
    }
    if flags.remesh_requested {
        demo.welded = !demo.welded;
    }
}

/// Mesh one chunk on its own, exactly as a dirty-set re-mesh would.
fn mesh_chunk<F: Sdf<Scalar = f32>>(
    layout: &ChunkLayout<f32>,
    field: &F,
    id: ChunkId,
) -> Option<MeshBuffer<f32>> {
    let shape = layout.sample_shape().ok()?;
    let mut out = MeshBuffer::<f32>::new();
    MarchingCubes::<f32>::new()
        .extract(
            field,
            &shape,
            layout.sample_origin(id),
            layout.cell_size(),
            &mut out,
        )
        .ok()?;
    Some(out)
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    demo: Res<Demo>,
    mut stats: ResMut<DemoStats>,
    mut overlay: ResMut<Overlay>,
    mut meshes: ResMut<Assets<Mesh>>,
    materials: Res<Materials>,
    mut commands: Commands,
    mut query: Query<(&mut Mesh3d, &mut Transform, &Half)>,
    mut last: Local<Option<(bool, bool, u32, usize)>>,
) {
    let key = (demo.welded, demo.exploded, demo.cells, demo.spacing);
    if *last == Some(key) {
        return;
    }
    *last = Some(key);

    let (h, spacing_name) = SPACINGS[demo.spacing];
    // Place the seam at x = 0, straight through the middle of the torus, with
    // the pair spanning the whole field in y and z. Anchoring the layout at a
    // round number instead would put the seam wherever it happened to land and
    // clip a corner of the surface, which shows the machinery and not the point.
    // The lower bound on `cells` is what keeps a chunk wide enough to contain
    // the torus at the coarser spacing.
    let width = demo.cells as f32 * h;
    let Ok(layout) = ChunkLayout::<f32>::new(demo.cells, h, [-width, -width * 0.5, -width * 0.5])
    else {
        return;
    };
    // A torus, because it crosses the seam plane in two separate places -- one
    // crossing could be a coincidence.
    let field = Torus::<f32>::canonical();
    let left_id = ChunkId::new([0, 0, 0]);
    let right_id = left_id.neighbour(0, 1);
    let started = Instant::now();
    let (Some(left), Some(right)) = (
        mesh_chunk(&layout, &field, left_id),
        mesh_chunk(&layout, &field, right_id),
    ) else {
        return;
    };
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;
    let seam_x = layout.sample_origin(right_id)[0];

    let mut joined = left.clone();
    joined.append(&right);
    let before = joined.vertex_count();

    let welded = if demo.welded {
        Welder::<f32>::new()
            .weld(&mut joined, h * ValidateConfig::WELD_EPSILON_REL as f32)
            .ok()
    } else {
        None
    };

    let Ok(cfg) = ValidateConfig::from_cell_size(f64::from(h)) else {
        return;
    };
    let (report, features) = validate_features(&joined.positions, &joined.indices, &cfg);

    // Split the boundary edges into the seam and the chunks' own outer faces.
    // Both are "boundary"; only one of them is a defect, and an overlay that
    // could not tell them apart would be showing alarm rather than information.
    let near_seam = |p: [f32; 3]| (p[0] - seam_x).abs() <= h * 0.25;
    overlay.seam_edges.clear();
    overlay.outer_edges.clear();
    overlay.cell_size = h;
    for [a, b] in &features.boundary_edges {
        let pa = joined.positions[*a as usize];
        let pb = joined.positions[*b as usize];
        let line = [Vec3::from(pa), Vec3::from(pb)];
        if near_seam(pa) && near_seam(pb) {
            overlay.seam_edges.push(line);
        } else {
            overlay.outer_edges.push(line);
        }
    }

    stats.title = format!(
        "E-115 chunk seam - {} cells, h = {spacing_name}",
        demo.cells
    );
    stats.vertices = joined.vertex_count();
    stats.triangles = joined.triangle_count();
    stats.extract_ms = extract_ms;
    stats.extra = vec![
        format!(
            "two chunks, meshed independently: {} + {} verts",
            left.vertex_count(),
            right.vertex_count()
        ),
        String::new(),
        match welded {
            Some(w) => format!(
                "WELDED   {before} -> {} verts  ({} merged, {} tris collapsed)",
                w.vertices_after,
                w.vertices_removed(),
                w.triangles_collapsed
            ),
            None => format!("UNWELDED {before} verts        [V] to weld"),
        },
        String::new(),
        format!(
            "seam boundary edges  {:>5}   <- red, the crack",
            overlay.seam_edges.len()
        ),
        format!(
            "outer boundary edges {:>5}   <- amber, the chunk's own faces",
            overlay.outer_edges.len()
        ),
        format!(
            "duplicate vertices   {:>5}   chi {}",
            report.duplicate_vertices, report.euler_characteristic
        ),
        String::new(),
        if overlay.seam_edges.is_empty() {
            "the seam carries no boundary: the two chunks are one surface".to_string()
        } else {
            "the seam is open: every red line is a triangle with no neighbour".to_string()
        },
        String::new(),
        "[V] weld  [E] explode  [1]/[2] spacing  [ ] resolution".to_string(),
    ];

    // Exploding shifts the halves apart. It changes nothing about the meshing --
    // it is there to make "these are two independent meshes" visible, since a
    // closed seam looks exactly like a single mesh, which is the point.
    let shift = if demo.exploded { h * 6.0 } else { 0.0 };
    let halves = [(Half::Left, -shift), (Half::Right, shift)];

    // One buffer holds both chunks, so the render mesh is built once and both
    // entities show it; exploding is the only reason there are two entities.
    let mut builder = MeshBuilder::new();
    for i in 0..joined.vertex_count() {
        builder.vertex(joined.positions[i], joined.normals[i]);
    }
    for t in joined.indices.chunks_exact(3) {
        builder.triangle(t[0], t[1], t[2]);
    }
    let handle = meshes.add(builder.into_mesh());

    if query.is_empty() {
        for (half, x) in halves {
            let material = match half {
                Half::Left => materials.left.clone(),
                Half::Right => materials.right.clone(),
            };
            commands.spawn((
                Mesh3d(handle.clone()),
                MeshMaterial3d(material),
                Transform::from_xyz(x, 0.0, 0.0),
                DemoMesh,
                half,
            ));
        }
    } else {
        for (mut mesh, mut transform, half) in &mut query {
            mesh.0 = handle.clone();
            transform.translation.x = halves
                .iter()
                .find(|(h, _)| h == half)
                .map_or(0.0, |(_, x)| *x);
        }
    }
}

fn draw_seam(overlay: Res<Overlay>, flags: Res<ViewFlags>, mut gizmos: Gizmos<SeamGizmos>) {
    const RED: Color = Color::srgb(1.0, 0.13, 0.13);
    const AMBER: Color = Color::srgb(1.0, 0.65, 0.1);

    if flags.grid {
        for [a, b] in &overlay.outer_edges {
            gizmos.line(*a, *b, AMBER);
        }
    }
    for [a, b] in &overlay.seam_edges {
        gizmos.line(*a, *b, RED);
    }
}
