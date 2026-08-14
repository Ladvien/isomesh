//! E-201 — a world larger than memory, streamed past a moving camera.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_terrain_stream --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 20-50x slower and the whole
//! point here is that meshing keeps up.
//!
//! `Space` pause the flight · `[` `]` view distance · `W` wireframe · `F12` shot.
//!
//! # What this puts together
//!
//! Everything in Phase 2 at once, and it is the first example where none of the
//! pieces are visible on their own:
//!
//! - **G-007** decides which chunks exist, with hysteresis, so a camera drifting
//!   across a threshold does not re-mesh the same chunk every frame.
//! - **B-003** meshes them on [`AsyncComputeTaskPool`], never in a system, and
//!   applies finished meshes under a frame budget.
//! - **G-001**'s layout is what makes a chunk's world position exact rather than
//!   merely close, which is what stops seams (M-32).
//!
//! The terrain is `fbm_terrain` sampled without bound — it is a function, so
//! there is no edge to reach. Fly for long enough and every chunk you see was
//! meshed while you were watching.
//!
//! # The number to watch
//!
//! Not the triangle count. **`ms/frame`, while chunks are landing.** A streaming
//! world that hitches is a streaming world that does its meshing on the main
//! thread, and the HUD reports the median of the last 30 frames precisely so a
//! single unlucky sample cannot be read as a stall — or a single lucky one as
//! proof of its absence.
//!
//! The budget is deliberately visible: `applied/frame` is what the frame budget
//! let through, and `waiting` is what the in-flight cap held back. Both being
//! non-zero while the frame time stays flat is the property this example exists
//! to show.

mod common;

use bevy::prelude::*;
use bevy_isomesh::{
    ChunkMesh, Extractor, IsomeshPlugin, MeshBudget, MeshStats, NeedsRemesh, VoxelChunk,
    VoxelVolume,
};
use common::{Capture, CommonPlugin, DemoStats, OrbitCamera, ViewFlags};
use isomesh::chunk::stream::{ChunkStream, StreamConfig, StreamUpdate};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::FbmTerrain;

/// Cells per chunk edge, and the world size of one chunk.
const CHUNK_CELLS: u32 = 16;
const CELL_SIZE: f32 = 0.5;

/// How far the camera flies per second.
const FLIGHT_SPEED: f32 = 9.0;

const MIN_VIEW: f32 = 16.0;
const MAX_VIEW: f32 = 96.0;

#[derive(Resource)]
struct World {
    volume: Entity,
    stream: ChunkStream,
    update: StreamUpdate,
    /// Load radius. The unload radius trails it by [`HYSTERESIS`].
    view: f32,
    flying: bool,
    travelled: f32,
}

/// Which chunk layers can contain the surface.
///
/// `fbm_terrain`'s amplitude is bounded, and a chunk is `CHUNK_CELLS * CELL_SIZE`
/// tall, so two layers straddle it with room to spare. Everything above is air
/// and everything below is rock; both mesh to nothing and cost the same as
/// terrain to find that out.
const VERTICAL_LAYERS: std::ops::RangeInclusive<i32> = -1..=0;

/// How much further than the load radius a chunk must go before it is dropped.
///
/// A quarter of the view distance. Too small and the camera's own jitter
/// re-meshes the boundary every frame; too large and memory holds chunks nobody
/// will look at again. G-007 refuses to construct a zero-width band at all.
const HYSTERESIS: f32 = 0.25;

/// Chunk entities keyed by id, so streaming can despawn by id without a scan.
#[derive(Resource, Default)]
struct Loaded(std::collections::HashMap<ChunkId, Entity>);

#[derive(Resource)]
struct Terrain(Handle<StandardMaterial>);

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-201 game terrain stream".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .add_plugins(IsomeshPlugin)
        .init_resource::<Loaded>()
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, fly, stream, attach_meshes, hud).chain())
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    for mut orbit in &mut camera {
        orbit.pitch = 0.42;
        orbit.yaw = 0.0;
        orbit.radius = 26.0;
    }

    let layout = ChunkLayout::new(CHUNK_CELLS, CELL_SIZE, [0.0; 3]).expect("a valid layout");
    // Marching Cubes rather than the subgrid extractor: M-98 measured that one
    // at 70x, which is the right trade for a sub-voxel feature and the wrong one
    // for a heightfield with nothing thinner than a cell in it. E-108 is where
    // the other choice earns its cost.
    let volume = commands
        .spawn(
            VoxelVolume::new(layout, FbmTerrain::<f32>::canonical())
                .with_extractor(Extractor::MarchingCubes),
        )
        .id();

    commands.insert_resource(World {
        volume,
        stream: ChunkStream::new(),
        update: StreamUpdate::new(),
        view: 40.0,
        flying: true,
        travelled: 0.0,
    });
    commands.insert_resource(Terrain(materials.add(StandardMaterial {
        base_color: Color::srgb(0.62, 0.66, 0.52),
        perceptual_roughness: 0.85,
        ..default()
    })));
    // A budget tight enough that the cap and the queue are both visible in the
    // HUD rather than being drained in the frame they appear.
    commands.insert_resource(MeshBudget {
        per_frame: std::time::Duration::from_micros(2_500),
        max_in_flight: 8,
    });
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut world: ResMut<World>,
    mut flags: ResMut<ViewFlags>,
) {
    if capture.is_active() {
        return;
    }
    if keys.just_pressed(KeyCode::Space) {
        world.flying = !world.flying;
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        world.view = (world.view + 8.0).min(MAX_VIEW);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        world.view = (world.view - 8.0).max(MIN_VIEW);
    }
    if keys.just_pressed(KeyCode::KeyR) {
        flags.remesh_requested = true;
    }
}

/// Fly the camera along a gentle arc, so the streamed set turns over
/// continuously rather than sweeping one axis and stopping.
fn fly(time: Res<Time>, mut world: ResMut<World>, mut camera: Query<&mut OrbitCamera>) {
    if !world.flying {
        return;
    }
    world.travelled += time.delta_secs() * FLIGHT_SPEED;
    let t = world.travelled;
    for mut orbit in &mut camera {
        orbit.focus = Vec3::new(t, 0.0, 18.0 * (t * 0.03).sin());
    }
}

/// Settle residency against the camera, and spawn or despawn to match.
fn stream(mut commands: Commands, mut world: ResMut<World>, mut loaded: ResMut<Loaded>) {
    let focus = Vec3::new(world.travelled, 0.0, 18.0 * (world.travelled * 0.03).sin());
    let view = world.view;
    let config = match StreamConfig::new(view, view * (1.0 + HYSTERESIS)) {
        Ok(config) => config,
        Err(_) => return,
    };

    let layout = ChunkLayout::new(CHUNK_CELLS, CELL_SIZE, [0.0; 3]).expect("a valid layout");
    let volume = world.volume;
    let World { stream, update, .. } = &mut *world;
    if stream
        .update(&layout, focus.to_array(), &config, update)
        .is_err()
    {
        return;
    }

    for id in update.unloaded.drain(..) {
        if let Some(entity) = loaded.0.remove(&id) {
            // Dropping the entity drops its `ChunkMesh`, which drops the last
            // handle to the asset -- so the mesh is freed without this example
            // touching `Assets` at all.
            commands.entity(entity).despawn();
        }
    }
    for id in update.loaded.drain(..) {
        // A heightfield wants a **slab**, not a ball. `ChunkStream` is
        // radius-based because that is what a general residency rule is, and
        // filtering it here is what a real game does: `fbm_terrain` has a
        // bounded amplitude, so every chunk outside these layers is either
        // empty air or solid rock and meshes to nothing at full cost.
        //
        // Without this the demo loads 952 chunks and leaves 606 permanently
        // waiting -- a queue that never drains, which shows up as holes in the
        // terrain that never fill. See M-104.
        if !VERTICAL_LAYERS.contains(&id.coords[1]) {
            continue;
        }
        let entity = commands
            .spawn((VoxelChunk { id, volume }, NeedsRemesh))
            .id();
        loaded.0.insert(id, entity);
    }
}

/// Give a finished chunk something to render with.
///
/// This is the line `bevy_isomesh` deliberately does not write: the plugin stops
/// at a `Handle<Mesh>` so that a consumer meshing for collision, or on a server,
/// never compiles `bevy_render`.
fn attach_meshes(
    mut commands: Commands,
    material: Res<Terrain>,
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
            MeshMaterial3d(material.0.clone()),
            // Chunk meshes are extracted in world space, so the transform is
            // identity -- the layout already put every vertex where it belongs.
            Transform::default(),
        ));
    }
}

#[allow(clippy::too_many_arguments)]
fn hud(
    world: Res<World>,
    stats: Res<MeshStats>,
    loaded: Res<Loaded>,
    meshes: Res<Assets<Mesh>>,
    chunks: Query<&ChunkMesh>,
    mut demo: ResMut<DemoStats>,
) {
    let mut vertices = 0usize;
    let mut triangles = 0usize;
    for chunk in &chunks {
        if let Some(mesh) = meshes.get(&chunk.0) {
            vertices += mesh.count_vertices();
            triangles += mesh.indices().map_or(0, |i| i.len() / 3);
        }
    }
    // Positions and normals are 3 f32 each, indices one u32. Close enough to be
    // the number a memory budget would use, and labelled as an estimate because
    // it ignores whatever the renderer does downstream.
    let bytes = vertices * (3 + 3) * 4 + triangles * 3 * 4;

    demo.title = format!(
        "E-201  streaming fbm terrain   view {:.0}   {}   [space] fly, [ ] view",
        world.view,
        if world.flying { "flying" } else { "paused" }
    );
    demo.vertices = vertices;
    demo.triangles = triangles;
    demo.extra = vec![
        format!(
            "{:<16} {:>8}   (hysteresis band {:.0} to {:.0})",
            "chunks resident",
            loaded.0.len(),
            world.view,
            world.view * (1.0 + HYSTERESIS)
        ),
        format!("{:<16} {:>8}", "meshing now", stats.in_flight),
        format!(
            "{:<16} {:>8}   held back by the in-flight cap",
            "waiting", stats.waiting
        ),
        format!(
            "{:<16} {:>8}   this frame, under the budget",
            "applied", stats.applied
        ),
        format!(
            "{:<16} {:>8.1} MB   estimated",
            "geometry",
            bytes as f64 / 1.0e6
        ),
        format!("{:<16} {:>8.0}", "travelled", world.travelled),
        String::new(),
        "the number to watch is ms/frame while chunks are landing. meshing".into(),
        "runs on the task pool (B-003); what the budget bounds is turning".into(),
        "finished extractions into assets, which is all that touches this".into(),
        "thread. hysteresis (G-007) is why drifting across the boundary does".into(),
        "not re-mesh the same chunk every frame.".into(),
    ];
}
