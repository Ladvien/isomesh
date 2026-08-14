//! E-206b — a body that slides, not a ray that hits.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_capsule_walk --release
//! ```
//!
//! **Always `--release`.**
//!
//! `Space` pause · `[` `]` view distance · `R` reset the walk.
//!
//! # What E-203 could not test
//!
//! E-203 casts 400 rays straight down every frame and asks *"is there a triangle
//! under this point"*. It answers that well — M-106: 495 seam crossings, zero
//! probes hitting nothing, worst vertical step **0.412 cells at a seam** against
//! **0.539 cells inside one chunk**, so the joins are measurably smoother than the
//! terrain they join.
//!
//! It cannot answer *"can a body move through here"*, because nothing in it
//! slides. A ray hits a surface or misses it; it never gets **caught** on one. A
//! seam lip that a ray reports as a 0.2-unit step is, to a moving capsule, either
//! nothing at all or a wall — and which one it is depends on the capsule, its
//! speed, and the direction it crosses in. That distinction is invisible to every
//! measurement E-203 makes.
//!
//! So this drives an actual rigid body across the same terrain and measures what
//! it is *prevented* from doing.
//!
//! # The measurement: commanded distance against travelled distance
//!
//! The capsule is asked to move at a fixed speed along a fixed path. Every frame
//! the demo records how far it was **asked** to go and how far it **actually
//! went**, horizontally. The difference is a stall, and a stall is the thing a
//! ray sweep structurally cannot see.
//!
//! A stall is not automatically a defect — terrain has slopes, and walking uphill
//! costs horizontal progress. So stalls are bucketed the same way M-106 bucketed
//! steps, and for the same reason: **the number that means something is the seam
//! figure compared against the interior figure**, not either alone. A fixed stall
//! threshold would be measuring the landscape.
//!
//! # Why a physics engine here, and where it is allowed to live
//!
//! V-22 settled this: a raycast cannot provide simulation, and Avian is built on
//! parry, so the geometry contract is the one G-005 already encodes — welded,
//! manifold, correctly wound. It is a **dev-dependency of `bevy_isomesh` only**.
//! The plugin stops at `Handle<Mesh>` precisely so a consumer picks their own
//! renderer, and picking their physics engine for them would be the same mistake.
//! Nothing a consumer of this crate depends on gains an edge to Avian.
//!
//! **It does cost a duplicate `glam` in this crate's dev tree** — Avian 0.7 pulls
//! `parry3d` 0.27, which brings `glam` 0.33 alongside Bevy 0.19's pinned 0.32.
//! `CLAUDE.md` warns that a `glam` mismatch compiles two copies silently, and that
//! is exactly what happens here. It is confined to `[dev-dependencies]`, so it
//! reaches examples and nothing else, and the whole workspace compiles clean —
//! but it is a real cost and it is written down rather than discovered later.

mod common;

use std::collections::HashMap;

use avian3d::prelude::*;
use bevy::prelude::*;
use bevy_isomesh::{
    ChunkMesh, Extractor, IsomeshPlugin, MeshBudget, NeedsRemesh, VoxelChunk, VoxelVolume,
};
use common::{Capture, CommonPlugin, DemoStats, OrbitCamera, ViewFlags};
use isomesh::chunk::stream::{ChunkStream, StreamConfig, StreamUpdate};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::FbmTerrain;

/// The same world E-203 walks, so the two are comparable.
const CHUNK_CELLS: u32 = 16;
const CELL_SIZE: f32 = 0.5;
const VERTICAL_LAYERS: std::ops::RangeInclusive<i32> = -1..=0;
const HYSTERESIS: f32 = 0.25;

const WALK_SPEED: f32 = 7.0;
const CAPSULE_RADIUS: f32 = 0.4;
const CAPSULE_LENGTH: f32 = 0.9;

/// How close to a chunk boundary counts as "at a seam", in world units.
///
/// One cell. A capsule caught *on* a seam is touching geometry from both sides of
/// it, and its centre is within roughly a cell of the plane when that happens.
const SEAM_BAND: f32 = CELL_SIZE;

/// Ignore the first moments after a chunk appears.
///
/// A body resting on terrain that has not finished streaming is not being
/// blocked, it is falling. Counting that as a stall would make the headline
/// number a measure of load latency.
const SETTLE_SECONDS: f32 = 0.75;

/// Dropped from just above the terrain rather than from the sky.
///
/// `fbm_terrain` peaks a little under 8 world units here, so this is a short
/// fall. The first draft spawned at 24 and spent most of the run falling, which
/// made every number taken before ~2 seconds a measurement of gravity.
const SPAWN_HEIGHT: f32 = 12.0;

/// How often the run prints its accumulated numbers.
///
/// A still is taken at frame 90 -- about 1.5 seconds -- which is far too early
/// for a walk statistic to mean anything, and no screenshot can carry a figure
/// that needs a minute to accumulate. So the measurement leaves through the log
/// and the screenshot is only the picture.
const REPORT_SECONDS: f32 = 5.0;

#[derive(Resource)]
struct Walk {
    volume: Entity,
    stream: ChunkStream,
    update: StreamUpdate,
    view: f32,
    walking: bool,
    /// Distance along the path used to steer. Advances from the first frame,
    /// including while the body is still falling.
    commanded: f32,
    /// Distance the body was asked to cover *while being measured*, and what it
    /// actually covered. Both accumulate in the same place for the same frames --
    /// the first draft advanced `commanded` from frame zero and `travelled` only
    /// after settling, and reported 12.8% progress on a body that was walking
    /// perfectly well.
    asked: f32,
    travelled: f32,
    elapsed: f32,
    /// Frames in which the body advanced less than it was asked to.
    stalls_at_seam: u32,
    stalls_inside: u32,
    /// Worst single-frame shortfall, as a fraction of what was commanded.
    worst_at_seam: f32,
    worst_inside: f32,
    /// Where the worst seam stall happened, for the HUD.
    worst_seam_at: Option<Vec3>,
    last_position: Option<Vec3>,
}

#[derive(Resource, Default)]
struct Loaded(HashMap<ChunkId, Entity>);

#[derive(Resource)]
struct Look {
    terrain: Handle<StandardMaterial>,
}

#[derive(Component)]
struct Capsule;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-206b game capsule walk".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .add_plugins(IsomeshPlugin)
        .add_plugins(PhysicsPlugins::default())
        .init_resource::<Loaded>()
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (
                controls,
                stream,
                attach_meshes,
                drive,
                measure,
                follow,
                report,
                hud,
            )
                .chain(),
        )
        .run();
}

fn layout() -> ChunkLayout<f32> {
    ChunkLayout::new(CHUNK_CELLS, CELL_SIZE, [0.0; 3]).expect("a valid layout")
}

/// Where the capsule is asked to be, horizontally, after `d` metres.
///
/// The same diagonal-with-a-wobble E-203 uses, and for the same reason: an
/// axis-aligned path crosses seams on one axis only, and the seams that matter
/// most are where three chunks meet.
fn path(d: f32) -> Vec2 {
    Vec2::new(d * 0.82, d * 0.57 + 9.0 * (d * 0.05).sin())
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    for mut orbit in &mut camera {
        orbit.pitch = 0.32;
        orbit.yaw = 0.0;
        orbit.radius = 24.0;
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
        commanded: 0.0,
        asked: 0.0,
        travelled: 0.0,
        elapsed: 0.0,
        stalls_at_seam: 0,
        stalls_inside: 0,
        worst_at_seam: 0.0,
        worst_inside: 0.0,
        worst_seam_at: None,
        last_position: None,
    });
    commands.insert_resource(Look {
        terrain: materials.add(StandardMaterial {
            base_color: Color::srgb(0.62, 0.66, 0.52),
            perceptual_roughness: 0.85,
            ..default()
        }),
    });
    commands.insert_resource(MeshBudget {
        per_frame: std::time::Duration::from_micros(3_000),
        max_in_flight: 8,
    });

    // The body. Dynamic rather than kinematic on purpose: a kinematic capsule is
    // moved by writing its transform, which means *it* decides where it ends up
    // and the terrain never gets to refuse. Only a dynamic body can be stopped by
    // geometry, and being stopped is the entire measurement.
    commands.spawn((
        RigidBody::Dynamic,
        Collider::capsule(CAPSULE_RADIUS, CAPSULE_LENGTH),
        // Without this the capsule tips over on the first slope and the run
        // becomes a demo of a falling cylinder.
        LockedAxes::ROTATION_LOCKED,
        Friction::new(0.4),
        Mesh3d(meshes.add(bevy::prelude::Capsule3d::new(CAPSULE_RADIUS, CAPSULE_LENGTH).mesh())),
        MeshMaterial3d(materials.add(StandardMaterial {
            base_color: Color::srgb(0.92, 0.42, 0.30),
            perceptual_roughness: 0.4,
            ..default()
        })),
        Transform::from_xyz(0.0, SPAWN_HEIGHT, 0.0),
        Capsule,
    ));
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut walk: ResMut<Walk>,
    mut capsule: Query<&mut Transform, With<Capsule>>,
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
        walk.commanded = 0.0;
        walk.asked = 0.0;
        walk.travelled = 0.0;
        walk.elapsed = 0.0;
        walk.stalls_at_seam = 0;
        walk.stalls_inside = 0;
        walk.worst_at_seam = 0.0;
        walk.worst_inside = 0.0;
        walk.worst_seam_at = None;
        walk.last_position = None;
        for mut t in &mut capsule {
            t.translation = Vec3::new(0.0, SPAWN_HEIGHT, 0.0);
        }
    }
}

fn stream(
    mut commands: Commands,
    mut walk: ResMut<Walk>,
    mut loaded: ResMut<Loaded>,
    capsule: Query<&Transform, With<Capsule>>,
) {
    // Residency follows the *body*, not the path. If the capsule is stuck, the
    // world must stay under it -- streaming ahead of a body that never arrives is
    // how a demo hides the very stall it is looking for.
    let Ok(here) = capsule.single() else {
        return;
    };
    let focus = [here.translation.x, 0.0, here.translation.z];
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

/// Give each freshly meshed chunk a renderer and a static collider.
///
/// The collider is built from the same `Handle<Mesh>` the renderer draws, which
/// is the point: a consumer hands the physics engine exactly what is on screen,
/// and any disagreement between the two would be a defect in this crate rather
/// than in the demo.
fn attach_meshes(
    mut commands: Commands,
    look: Res<Look>,
    meshes: Res<Assets<Mesh>>,
    fresh: Query<(Entity, &ChunkMesh), Without<Mesh3d>>,
) {
    for (entity, chunk) in &fresh {
        let Some(mesh) = meshes.get(&chunk.0) else {
            continue;
        };
        // An empty chunk -- all air or all rock -- has no triangles and no
        // collider. Skipping it is correct; trying to build one is not.
        let Some(collider) = Collider::trimesh_from_mesh(mesh) else {
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
                Mesh3d(chunk.0.clone()),
                MeshMaterial3d(look.terrain.clone()),
                Transform::default(),
            ));
            continue;
        };
        commands.entity(entity).try_insert((
            Mesh3d(chunk.0.clone()),
            MeshMaterial3d(look.terrain.clone()),
            Transform::default(),
            RigidBody::Static,
            collider,
        ));
    }
}

/// Command the horizontal velocity, and let gravity and contacts do the rest.
///
/// Only the horizontal components are written. Writing `y` too would mean
/// overriding gravity every frame, and a body that cannot fall cannot be caught
/// on anything either.
fn drive(
    time: Res<Time>,
    mut walk: ResMut<Walk>,
    mut capsule: Query<(&Transform, &mut LinearVelocity), With<Capsule>>,
) {
    let dt = time.delta_secs();
    walk.elapsed += dt;
    let Ok((transform, mut velocity)) = capsule.single_mut() else {
        return;
    };
    if !walk.walking {
        velocity.0.x = 0.0;
        velocity.0.z = 0.0;
        return;
    }

    // Steer toward the point on the path a little ahead of the commanded
    // distance, so the body follows the curve rather than a fixed heading.
    walk.commanded += WALK_SPEED * dt;
    let target = path(walk.commanded);
    let to =
        Vec2::new(target.x, target.y) - Vec2::new(transform.translation.x, transform.translation.z);
    let heading = to.normalize_or_zero();
    velocity.0.x = heading.x * WALK_SPEED;
    velocity.0.z = heading.y * WALK_SPEED;
}

/// Compare what the body was asked to do with what it did.
fn measure(time: Res<Time>, mut walk: ResMut<Walk>, capsule: Query<&Transform, With<Capsule>>) {
    let Ok(transform) = capsule.single() else {
        return;
    };
    let here = transform.translation;
    let Some(last) = walk.last_position.replace(here) else {
        return;
    };
    if !walk.walking || walk.elapsed < SETTLE_SECONDS {
        return;
    }

    let moved = Vec2::new(here.x - last.x, here.z - last.z).length();
    let asked = WALK_SPEED * time.delta_secs();
    if asked <= f32::EPSILON {
        return;
    }
    walk.travelled += moved;
    walk.asked += asked;
    let shortfall = ((asked - moved) / asked).clamp(0.0, 1.0);
    // A body on a slope loses horizontal progress legitimately, so only a
    // substantial shortfall counts as a stall at all. The threshold is not the
    // measurement -- the seam/interior *comparison* is.
    if shortfall < 0.5 {
        return;
    }

    // How close is the capsule to a chunk boundary? A chunk spans
    // CHUNK_CELLS * CELL_SIZE, so distance to the nearest plane on each axis is
    // the distance to the nearest multiple of that.
    let span = CHUNK_CELLS as f32 * CELL_SIZE;
    let to_plane = |v: f32| {
        let m = (v / span).round() * span;
        (v - m).abs()
    };
    let at_seam = to_plane(here.x).min(to_plane(here.z)) < SEAM_BAND;

    if at_seam {
        walk.stalls_at_seam += 1;
        if shortfall > walk.worst_at_seam {
            walk.worst_at_seam = shortfall;
            walk.worst_seam_at = Some(here);
        }
    } else {
        walk.stalls_inside += 1;
        if shortfall > walk.worst_inside {
            walk.worst_inside = shortfall;
        }
    }
}

/// Keep the camera on the body.
///
/// Without this the orbit stays at the origin and the capsule walks out of frame
/// within a few seconds -- so the one thing the demo exists to show is off screen
/// for every capture after the first.
fn follow(capsule: Query<&Transform, With<Capsule>>, mut camera: Query<&mut OrbitCamera>) {
    let Ok(here) = capsule.single() else {
        return;
    };
    for mut orbit in &mut camera {
        orbit.focus = here.translation;
    }
}

/// One CSV line every few seconds, so the walk can be measured over a minute
/// instead of over the frame a screenshot happens to land on.
fn report(time: Res<Time>, walk: Res<Walk>, mut next: Local<f32>) {
    if walk.elapsed < SETTLE_SECONDS || walk.elapsed < *next {
        return;
    }
    *next = walk.elapsed + REPORT_SECONDS;
    let _ = time;
    info!(
        "capsule,{:.1},{:.2},{:.2},{:.1},{},{},{:.3},{:.3}",
        walk.elapsed,
        walk.asked,
        walk.travelled,
        if walk.asked > 0.0 {
            100.0 * walk.travelled / walk.asked
        } else {
            100.0
        },
        walk.stalls_at_seam,
        walk.stalls_inside,
        walk.worst_at_seam,
        walk.worst_inside,
    );
}

fn hud(
    walk: Res<Walk>,
    loaded: Res<Loaded>,
    meshes: Res<Assets<Mesh>>,
    chunks: Query<&ChunkMesh>,
    mut stats: ResMut<DemoStats>,
    mut flags: ResMut<ViewFlags>,
) {
    flags.grid = false;
    // What the body is actually standing on, summed over resident chunks.
    let (mut vertices, mut triangles) = (0usize, 0usize);
    for chunk in &chunks {
        if let Some(mesh) = meshes.get(&chunk.0) {
            vertices += mesh.count_vertices();
            triangles += mesh.indices().map_or(0, |i| i.len() / 3);
        }
    }
    stats.vertices = vertices;
    stats.triangles = triangles;

    let efficiency = if walk.asked > 0.0 {
        100.0 * walk.travelled / walk.asked
    } else {
        100.0
    };
    let verdict = if walk.elapsed < SETTLE_SECONDS {
        "settling -- the body is still falling to the ground".to_string()
    } else if walk.stalls_at_seam == 0 {
        "NO SEAM STALL -- a moving body crosses every join it has met".to_string()
    } else if walk.worst_at_seam <= walk.worst_inside {
        "seams are no worse than the terrain -- which is the comparison that counts".to_string()
    } else {
        format!(
            "!! SEAM WORSE THAN TERRAIN -- {:.0}% against {:.0}%",
            100.0 * walk.worst_at_seam,
            100.0 * walk.worst_inside
        )
    };

    stats.title = format!(
        "E-206b  capsule walk   {} chunks resident   {}",
        loaded.0.len(),
        if walk.walking { "walking" } else { "paused" }
    );
    stats.extra = vec![
        format!(
            "{:<26} {:>9.1}   metres asked for, once settled",
            "asked", walk.asked
        ),
        format!(
            "{:<26} {:>9.1}   metres actually covered",
            "travelled", walk.travelled
        ),
        format!("{:<26} {:>9.1}%  of what was asked", "progress", efficiency),
        String::new(),
        "a stall is a frame where the body advanced less than half what it was".into(),
        "asked to. bucketed by where it happened, because a slope costs progress".into(),
        "legitimately and a fixed threshold would just measure the landscape:".into(),
        String::new(),
        format!(
            "{:<26} {:>9}   frames, worst shortfall {:.0}%",
            "stalls AT a seam",
            walk.stalls_at_seam,
            100.0 * walk.worst_at_seam
        ),
        format!(
            "{:<26} {:>9}   frames, worst shortfall {:.0}%",
            "stalls INSIDE a chunk",
            walk.stalls_inside,
            100.0 * walk.worst_inside
        ),
        match walk.worst_seam_at {
            Some(p) => format!(
                "{:<26} {:>9}",
                "worst seam stall at",
                format!("{:.1},{:.1},{:.1}", p.x, p.y, p.z)
            ),
            None => format!("{:<26} {:>9}", "worst seam stall at", "--"),
        },
        String::new(),
        verdict,
        String::new(),
        "E-203 casts 400 rays a frame and finds no holes (M-106). a ray cannot be".into(),
        "CAUGHT on anything, so it cannot find a lip that stops a body. this is".into(),
        "the same terrain with a rigid body on it, and the number that means".into(),
        "something is the seam row against the interior row -- not either alone.".into(),
        String::new(),
        format!(
            "view {:.0}m   [ and ] to change   Space pauses   R resets",
            walk.view
        ),
    ];
}
