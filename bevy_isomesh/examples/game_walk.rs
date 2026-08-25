//! E-203 — walk every chunk seam, and count what happens.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_walk --release
//! ```
//!
//! **Always `--release`.**
//!
//! `Space` pause the walk · `[` `]` view distance · `W` wireframe · `F12` shot.
//!
//! # This one is designed to fail
//!
//! Every other example here demonstrates something. This one tries to falsify
//! a foundational assumption, and the backlog says so: *"the acid test. Walk
//! every chunk seam. No falling through, no invisible walls. If this fails,
//! G-001's overlap is wrong."*
//!
//! So it does not just let you walk around and look. A walker crosses seam after
//! seam, and every step casts a ray straight down against the **meshed
//! triangles** — through `parry3d`, the same library a consumer would hand this
//! geometry to. Three numbers come back:
//!
//! - **fall-throughs**: steps where the ray hit nothing. A gap between two
//!   chunks is invisible from above until something tries to stand on it.
//! - **worst step**: the largest vertical discontinuity between consecutive
//!   ground samples, in cells. A seam that does not line up is a lip, and a lip
//!   is an invisible wall.
//! - **seams crossed**: how much of the test has actually been exercised, so a
//!   clean run cannot be a run that never reached a boundary (M-44).
//!
//! # What makes the answer trustworthy
//!
//! The ray is cast against the mesh, not against the field. Asking the field
//! would test the field, which is not in doubt — G-001's overlap is what decides
//! whether two independently meshed chunks *meet*, and only the triangles know
//! that.
//!
//! And the worst step is compared against the terrain's own local slope rather
//! than against zero: real terrain has real steps, and a threshold that flagged
//! them would be measuring the landscape instead of the seams.

mod common;

use std::collections::HashMap;

use bevy::prelude::*;
use bevy_isomesh::{
    ChunkMesh, Extractor, IsomeshPlugin, MeshBudget, NeedsRemesh, VoxelChunk, VoxelVolume,
};
use common::{Capture, CommonPlugin, DemoStats, OrbitCamera, ViewFlags};
use isomesh::chunk::stream::{ChunkStream, StreamConfig, StreamUpdate};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::FbmTerrain;
use parry3d::math::{Pose, Vec3 as PVec3};
use parry3d::query::{Ray, RayCast};
use parry3d::shape::TriMesh;

const CHUNK_CELLS: u32 = 16;
const CELL_SIZE: f32 = 0.5;
/// Chunk layers that can hold the surface. See E-201 and M-104.
const VERTICAL_LAYERS: std::ops::RangeInclusive<i32> = -1..=0;

const WALK_SPEED: f32 = 7.0;
const HYSTERESIS: f32 = 0.25;

/// How high above the terrain the ray starts, and how far it reaches.
const RAY_HEIGHT: f32 = 40.0;
const RAY_LENGTH: f32 = 120.0;

#[derive(Resource)]
struct Walk {
    volume: Entity,
    stream: ChunkStream,
    update: StreamUpdate,
    view: f32,
    walking: bool,
    travelled: f32,
    /// Where the walker last stood, once the ground has been found.
    stood: Option<Vec3>,
    /// Adjacent probe pairs tested that straddle a chunk boundary.
    seam_pairs: usize,
    /// Probes that hit nothing at all.
    holes: usize,
    /// Largest vertical discontinuity across a chunk boundary, in cells.
    worst_lip_cells: f32,
    /// Largest discontinuity *not* across a boundary, for comparison -- this is
    /// the terrain's own roughness, and the seam figure is only meaningful
    /// against it.
    worst_interior_cells: f32,
}

#[derive(Resource, Default)]
struct Loaded(HashMap<ChunkId, Entity>);

/// One `parry3d` collider per chunk, built from the mesh the plugin produced.
#[derive(Resource, Default)]
struct Colliders(HashMap<ChunkId, TriMesh>);

#[derive(Resource)]
struct Look {
    terrain: Handle<StandardMaterial>,
    walker: Handle<StandardMaterial>,
    sphere: Handle<Mesh>,
}

#[derive(Component)]
struct Walker;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-203 game walk".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .add_plugins(IsomeshPlugin)
        .init_resource::<Loaded>()
        .init_resource::<Colliders>()
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (controls, stream, attach_meshes, build_colliders, step, hud).chain(),
        )
        .run();
}

fn layout() -> ChunkLayout<f32> {
    ChunkLayout::new(CHUNK_CELLS, CELL_SIZE, [0.0; 3]).expect("a valid layout")
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    for mut orbit in &mut camera {
        orbit.pitch = 0.30;
        orbit.yaw = 0.0;
        orbit.radius = 22.0;
    }

    let volume = commands
        .spawn(
            VoxelVolume::new(layout(), FbmTerrain::<f32>::canonical())
                .with_extractor(Extractor::MarchingCubes),
        )
        .id();

    commands.insert_resource(Walk {
        volume,
        stream: ChunkStream::new(),
        update: StreamUpdate::new(),
        view: 32.0,
        walking: true,
        travelled: 0.0,
        stood: None,
        seam_pairs: 0,
        holes: 0,
        worst_lip_cells: 0.0,
        worst_interior_cells: 0.0,
    });
    commands.insert_resource(Look {
        terrain: materials.add(StandardMaterial {
            base_color: Color::srgb(0.62, 0.66, 0.52),
            perceptual_roughness: 0.85,
            ..default()
        }),
        walker: materials.add(StandardMaterial {
            base_color: Color::srgb(0.92, 0.42, 0.30),
            perceptual_roughness: 0.4,
            ..default()
        }),
        sphere: meshes.add(Sphere::new(0.45).mesh().uv(24, 12)),
    });
    commands.insert_resource(MeshBudget {
        per_frame: std::time::Duration::from_micros(3_000),
        max_in_flight: 8,
    });
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut walk: ResMut<Walk>,
    mut flags: ResMut<ViewFlags>,
) {
    if capture.is_active() {
        return;
    }
    if keys.just_pressed(KeyCode::Space) {
        walk.walking = !walk.walking;
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        walk.view = (walk.view + 8.0).min(80.0);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        walk.view = (walk.view - 8.0).max(16.0);
    }
    if keys.just_pressed(KeyCode::KeyR) {
        flags.remesh_requested = true;
    }
}

/// Where the walker is, horizontally, after `travelled` metres.
///
/// A diagonal with a slow wobble, deliberately: an axis-aligned path crosses
/// seams on one axis only, and the seams that matter most are the ones where
/// three chunks meet.
fn path(travelled: f32) -> Vec2 {
    Vec2::new(
        travelled * 0.82,
        travelled * 0.57 + 9.0 * (travelled * 0.05).sin(),
    )
}

fn stream(mut commands: Commands, mut walk: ResMut<Walk>, mut loaded: ResMut<Loaded>) {
    let here = path(walk.travelled);
    let focus = [here.x, 0.0, here.y];
    let view = walk.view;
    let Ok(config) = StreamConfig::new(view, view * (1.0 + HYSTERESIS)) else {
        return;
    };

    let layout = layout();
    let volume = walk.volume;
    let Walk { stream, update, .. } = &mut *walk;
    if stream.update(&layout, focus, &config, update).is_err() {
        return;
    }

    for id in update.unloaded.drain(..) {
        if let Some(entity) = loaded.0.remove(&id) {
            commands.entity(entity).despawn();
        }
    }
    for id in update.loaded.drain(..) {
        if !VERTICAL_LAYERS.contains(&id.coords[1]) {
            continue;
        }
        let entity = commands
            .spawn((VoxelChunk { id, volume }, NeedsRemesh))
            .id();
        loaded.0.insert(id, entity);
    }
}

fn attach_meshes(
    mut commands: Commands,
    look: Res<Look>,
    fresh: Query<(Entity, &ChunkMesh), Without<Mesh3d>>,
) {
    for (entity, mesh) in &fresh {
        // `try_insert`, not `insert`, and this is not a swallowed error. The
        // streamer despawns a chunk the moment it leaves residency, and that
        // can happen in the same frame this system queues the attach -- both
        // are commands, and whichever order the queue flushes, the insert may
        // land on an entity that is already gone. Attaching a mesh to a chunk
        // that no longer exists is *vacuous*, not degraded: there is nothing
        // else this could correctly do, and `insert` panics instead.
        // `game_terrain_stream` is where it actually fired -- "Entity
        // despawned: the entity with ID 404v0 is invalid" -- and the other
        // three carried the same race latently.
        commands.entity(entity).try_insert((
            Mesh3d(mesh.0.clone()),
            MeshMaterial3d(look.terrain.clone()),
            Transform::default(),
        ));
    }
}

/// Turn each chunk's mesh into the collider a consumer would use.
fn build_colliders(
    meshes: Res<Assets<Mesh>>,
    mut colliders: ResMut<Colliders>,
    loaded: Res<Loaded>,
    chunks: Query<(&VoxelChunk, &ChunkMesh)>,
) {
    for (chunk, mesh) in &chunks {
        if colliders.0.contains_key(&chunk.id) {
            continue;
        }
        let Some(asset) = meshes.get(&mesh.0) else {
            continue;
        };
        let Some(collider) = to_trimesh(asset) else {
            continue;
        };
        colliders.0.insert(chunk.id, collider);
    }
    // Drop colliders for chunks that are no longer resident, or the map grows
    // without bound as the walk goes on.
    colliders.0.retain(|id, _| loaded.0.contains_key(id));
}

fn to_trimesh(mesh: &Mesh) -> Option<TriMesh> {
    let positions = mesh.attribute(Mesh::ATTRIBUTE_POSITION)?.as_float3()?;
    let indices = mesh.indices()?;
    let vertices: Vec<PVec3> = positions
        .iter()
        .map(|p| PVec3::new(p[0], p[1], p[2]))
        .collect();
    let mut faces: Vec<[u32; 3]> = Vec::with_capacity(indices.len() / 3);
    let flat: Vec<u32> = indices.iter().map(|i| i as u32).collect();
    for tri in flat.as_chunks::<3>().0 {
        faces.push([tri[0], tri[1], tri[2]]);
    }
    if vertices.is_empty() || faces.is_empty() {
        return None;
    }
    TriMesh::new(vertices, faces).ok()
}

/// Advance the walker, then probe the seams ahead of it.
///
/// The walker is the demo; the **probe sweep** is the test. Waiting for one
/// walker to reach a seam tests one seam per few seconds, which is not a
/// measurement — so every frame also samples a dense transect across the loaded
/// region and asks the same question at every point of it.
fn step(
    time: Res<Time>,
    colliders: Res<Colliders>,
    mut walk: ResMut<Walk>,
    mut commands: Commands,
    look: Res<Look>,
    mut walker: Query<&mut Transform, With<Walker>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    if walk.walking {
        walk.travelled += time.delta_secs() * WALK_SPEED;
    }
    let here = path(walk.travelled);
    let layout = layout();

    // Move the visible walker.
    if let Some(y) = ground_at(&colliders.0, &layout, here) {
        let standing = Vec3::new(here.x, y, here.y);
        walk.stood = Some(standing);
        if walker.is_empty() {
            commands.spawn((
                Mesh3d(look.sphere.clone()),
                MeshMaterial3d(look.walker.clone()),
                Transform::from_translation(standing + Vec3::Y * 0.45),
                Walker,
            ));
        } else {
            for mut transform in &mut walker {
                transform.translation = standing + Vec3::Y * 0.45;
            }
        }
        for mut orbit in &mut camera {
            orbit.focus = standing;
        }
    }

    // The transect. Fine enough that consecutive probes are a fraction of a
    // cell apart, so a lip at a seam cannot hide between two samples.
    const PROBES: usize = 400;
    let span = walk.view * 1.4;
    let step_len = span / PROBES as f32;
    let direction = Vec2::new(0.83, 0.56).normalize();
    let start = here - direction * (span * 0.5);

    let mut previous: Option<(Vec2, f32, ChunkId)> = None;
    for i in 0..PROBES {
        let at = start + direction * (i as f32 * step_len);
        let chunk = layout.chunk_of([at.x, 0.0, at.y]);
        let Some(y) = ground_at(&colliders.0, &layout, at) else {
            // A miss only counts as a hole once **every** layer that could
            // hold the surface has been meshed. `||` here is wrong and reads as
            // a catastrophe: the surface sits in whichever layer it sits in, and
            // with one layer still in the queue a perfectly sound column reports
            // a hole. That is what the first version did -- 439 of them, and a
            // verdict declaring G-001's overlap broken. See M-105.
            if VERTICAL_LAYERS.clone().all(|layer| {
                colliders
                    .0
                    .contains_key(&ChunkId::new([chunk.coords[0], layer, chunk.coords[2]]))
            }) {
                walk.holes += 1;
            }
            previous = None;
            continue;
        };

        if let Some((was_at, was_y, was_chunk)) = previous {
            let travelled = (at - was_at).length();
            if travelled > 1.0e-5 {
                let jump_cells = (y - was_y).abs() / CELL_SIZE;
                if was_chunk != chunk {
                    walk.seam_pairs += 1;
                    if jump_cells > walk.worst_lip_cells {
                        walk.worst_lip_cells = jump_cells;
                    }
                } else if jump_cells > walk.worst_interior_cells {
                    walk.worst_interior_cells = jump_cells;
                }
            }
        }
        previous = Some((at, y, chunk));
    }
}

/// The ground height at a point, from the colliders of the chunks above and
/// below it.
///
/// Only the two vertical layers that can contain the surface are consulted,
/// rather than every resident collider: a transect of 400 probes against 200
/// colliders would be 80,000 ray casts a frame, and the layout already says
/// which chunk a point is in.
fn ground_at(
    colliders: &HashMap<ChunkId, TriMesh>,
    layout: &ChunkLayout<f32>,
    at: Vec2,
) -> Option<f32> {
    let column = layout.chunk_of([at.x, 0.0, at.y]);
    let origin = PVec3::new(at.x, RAY_HEIGHT, at.y);
    let ray = Ray::new(origin, PVec3::new(0.0, -1.0, 0.0));
    let mut best: Option<f32> = None;
    for layer in VERTICAL_LAYERS {
        let id = ChunkId::new([column.coords[0], layer, column.coords[2]]);
        let Some(collider) = colliders.get(&id) else {
            continue;
        };
        if let Some(toi) = collider.cast_ray(&Pose::IDENTITY, &ray, RAY_LENGTH, false) {
            let y = RAY_HEIGHT - toi;
            best = Some(best.map_or(y, |current: f32| current.max(y)));
        }
    }
    best
}

fn hud(
    walk: Res<Walk>,
    loaded: Res<Loaded>,
    colliders: Res<Colliders>,
    mut demo: ResMut<DemoStats>,
) {
    demo.title = format!(
        "E-203  walking the seams   {} crossings tested   {}   [space] walk, [ ] view",
        walk.seam_pairs,
        if walk.walking { "walking" } else { "paused" }
    );

    let verdict = if walk.seam_pairs < 200 {
        "not enough seam crossings tested yet to mean anything"
    } else if walk.holes > 0 {
        "HOLE at a seam -- G-001's overlap is wrong"
    } else if walk.worst_lip_cells > walk.worst_interior_cells * 1.5 {
        "a lip at the seams that the terrain itself does not explain"
    } else {
        "no holes, and seams no rougher than the terrain. they hold"
    };

    demo.extra = vec![
        format!(
            "{:<18} {:>8}   adjacent probes straddling a chunk boundary",
            "seam crossings", walk.seam_pairs
        ),
        format!(
            "{:<18} {:>8}   probes that hit nothing over a resident chunk",
            "holes", walk.holes
        ),
        format!(
            "{:<18} {:>8.3}   cells, worst jump ACROSS a seam",
            "seam lip", walk.worst_lip_cells
        ),
        format!(
            "{:<18} {:>8.3}   cells, worst jump WITHIN one chunk",
            "terrain roughness", walk.worst_interior_cells
        ),
        String::new(),
        format!("{:<18} {:>8}", "chunks resident", loaded.0.len()),
        format!("{:<18} {:>8}", "colliders", colliders.0.len()),
        String::new(),
        format!("verdict: {verdict}"),
        String::new(),
        "the ray is cast against the meshed triangles through parry3d, not".into(),
        "against the field. asking the field would test the field, which is".into(),
        "not in doubt -- whether two independently meshed chunks *meet* is".into(),
        "decided by G-001's overlap, and only the triangles know it.".into(),
    ];
}
