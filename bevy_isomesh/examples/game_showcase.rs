//! E-210 — a world with a roof over your head.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_showcase --release
//! ```
//!
//! **Always `--release`.**
//!
//! `Space` pause · `[` `]` fly speed · `1`–`3` how much of the world is cave.
//!
//! # Why this example exists at all
//!
//! Every other gameplay demo in this repository runs on `fbm_terrain`, which is a
//! **heightfield**: one height per column, no overhangs, no caves, nothing above
//! anything else. That is a perfectly good test field and it is the wrong
//! advertisement, because a heightfield is exactly the case you do *not* need a
//! voxel mesher for. A quadtree of grids meshes one in a fraction of the time.
//!
//! So the demos proved the crate correct and quietly undersold what it is for.
//! This one flies a camera **through** the terrain — under arches, into tunnels,
//! out the far side of a wall — which a heightfield cannot represent at all.
//!
//! # The field, and why it is this one
//!
//! ```text
//! solid(p) = max( p.y − height(x, z) ,  |gyroid(p)| − thickness )
//! ```
//!
//! A `max` is an intersection, so material exists only where a point is **below
//! the terrain surface** *and* **inside a thickened gyroid**. The gyroid is
//! triply periodic — it tunnels in all three axes by construction — so the
//! intersection is a landscape shot through with caves that connect, arches that
//! carry rock over open air, and ceilings. None of those three things has a
//! height.
//!
//! It is also honest about where the geometry comes from: the gyroid is the
//! crate's own reference field, the terrain is a sum of sines, and the whole
//! composition is nine lines. Nothing here is authored, sculpted or baked.
//!
//! # What is still measured, with the HUD off
//!
//! Capture this with `ISOMESH_VIEW=nohud` and it is a picture; leave the HUD on
//! and it is still a demo. Chunks resident, triangles, frame time and the
//! streaming budget are all on screen, because a showcase that hides its cost is
//! the thing this repository exists not to be.
//!
//! # Why Marching Cubes, and not the sharper extractors
//!
//! Because this world is **chunked**, and the dual methods do not tile. Marching
//! Cubes places vertices on grid *edges*, which two neighbouring chunks compute
//! from identical corner values and agree on exactly. Surface Nets and Dual
//! Contouring place one vertex per cell *interior*, and a boundary quad needs the
//! neighbour's vertex — which the chunk does not have, so it stops short.
//!
//! Measured on two adjacent chunks, boundary edges lying in the shared plane:
//! **Marching Cubes 0, Surface Nets 5, Dual Contouring 4** (M-128). This example
//! shipped on Dual Contouring for one commit, having been switched during an
//! unrelated aliasing investigation, and the README's lead image was a cracked
//! world until someone looked at it.
//!
//! # Scope
//!
//! Streaming, budgeting and the field. **Carving is deliberately not here** —
//! `game_dig` and `game_destruction` already show it, and the acceptance for this
//! ticket is that a viewer who reads no numbers can see something a heightfield
//! cannot do. Adding a brush would not have moved that.

mod common;

use std::collections::HashMap;

use bevy::prelude::*;
use bevy_isomesh::{
    ChunkMesh, Extractor, IsomeshPlugin, MeshBudget, NeedsRemesh, VoxelChunk, VoxelVolume,
};
use common::{Capture, CommonPlugin, DemoStats, OrbitCamera, ViewFlags};
use isomesh::Sdf;
use isomesh::brush::smooth_min;
use isomesh::chunk::stream::{ChunkStream, StreamConfig, StreamUpdate};
use isomesh::chunk::{ChunkId, ChunkLayout};

const CHUNK_CELLS: u32 = 16;
const CELL_SIZE: f32 = 0.5;
/// Four layers, not two. A cave world needs room *under* the ground for the
/// caves to be in; `game_walk`'s two layers are enough for a heightfield and
/// would clip the roofs off this one.
const VERTICAL_LAYERS: std::ops::RangeInclusive<i32> = -3..=1;
const HYSTERESIS: f32 = 0.25;

/// A landscape intersected with a thickened gyroid.
///
/// Negative is solid, as everywhere in this crate. The gyroid term is what makes
/// this impossible to store as a height: it is triply periodic, so it tunnels in
/// `x`, `y` *and* `z`, and intersecting it with a half-space under the terrain
/// leaves caves that connect and rock that hangs over open ground.
#[derive(Clone, Copy)]
struct CaveWorld {
    /// Iso-level on the gyroid: rock is `g < thickness`. Larger is more rock
    /// and narrower tunnels.
    thickness: f32,
}

impl Sdf for CaveWorld {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        // A rolling surface. Sines rather than the crate's fBm because this one
        // has to read as landscape at a glance, and three octaves of sine is
        // both enough and legible in the source.
        let height = 2.6 * (p[0] * 0.11).sin() * (p[2] * 0.09).cos()
            + 1.1 * (p[0] * 0.23 + 1.7).sin()
            + 0.8 * (p[2] * 0.27 - 0.6).cos();
        let ground = p[1] - height;

        // The gyroid's period is `2*pi/K`, and that number is the whole look.
        // At K = 0.19 the period is 33 m, so a single tunnel is wider than the
        // view and the world reads as a cracked plain rather than as caves. At
        // 0.5 the period is 12.6 m and several arches are in frame at once.
        const K: f32 = 0.5;
        let g = (K * p[0]).sin() * (K * p[1]).cos()
            + (K * p[1]).sin() * (K * p[2]).cos()
            + (K * p[2]).sin() * (K * p[0]).cos();

        // Intersection: below the surface AND inside one of the two labyrinths
        // the gyroid separates.
        //
        // **`g - level`, not `|g| - thickness`.** The absolute value gives the
        // *shell* around the gyroid surface, and a shell's thickness is
        // `2t/|grad g|`, which varies: wherever the gradient is steep the sheet
        // pinches below a cell and the mesher tears it into stair-stepped
        // slashes. Thickening the shell moves the pinch without removing it,
        // because the pinch is where the shell's own rim is. One side of the
        // labyrinth has no rim and no thin sheets -- it is bulk rock with
        // tunnels through it, which is what a cave is.
        //
        // **Smoothed, and not for looks.** A hard `max` puts a knife edge
        // wherever the two surfaces cross, and that rim is thinner than a cell
        // over much of its length, so it aliases into a sawtooth -- M-72's
        // failure mode, arriving from the field rather than from the extractor.
        // Switching Surface Nets for Dual Contouring did not touch it, which is
        // what ruled the extractor out. A smooth intersection rounds the rim to
        // `BLEND` wide, comfortably above the 0.5 m cell, and real rock does not
        // have knife edges either.
        //
        // `smooth_max(a, b) = -smooth_min(-a, -b)`.
        const BLEND: f32 = 0.7;
        -smooth_min(-ground, -(g - self.thickness), BLEND)
    }
}

#[derive(Resource)]
struct Show {
    volume: Entity,
    stream: ChunkStream,
    update: StreamUpdate,
    view: f32,
    flying: bool,
    speed: f32,
    travelled: f32,
    thickness: f32,
}

#[derive(Resource, Default)]
struct Loaded(HashMap<ChunkId, Entity>);

#[derive(Resource)]
struct Look(Handle<StandardMaterial>);

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-210 game showcase".into(),
                ..default()
            }),
            ..default()
        }))
        // A sky, so the picture has somewhere to end. The harness leaves this at
        // the default charcoal, which is right for a diagram and wrong for a
        // landscape.
        .insert_resource(ClearColor(Color::srgb(0.42, 0.55, 0.68)))
        .add_plugins(CommonPlugin)
        .add_plugins(IsomeshPlugin)
        .init_resource::<Loaded>()
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, stream, attach, fly, hud).chain())
        .run();
}

fn layout() -> ChunkLayout<f32> {
    ChunkLayout::new(CHUNK_CELLS, CELL_SIZE, [0.0; 3]).expect("a valid layout")
}

/// How much of the labyrinth is rock.
///
/// The gyroid runs about `-1.5..1.5`, and `g < level` selects the solid side, so
/// `0.0` is an even split and larger is more rock with narrower tunnels.
fn thickness_for(step: u8) -> f32 {
    match step {
        0 => -0.35,
        1 => 0.25,
        _ => 0.85,
    }
}

/// The shell the demo opens on. See [`thickness_for`].
const DEFAULT_STEP: u8 = 1;

/// An `f32` from the environment, or `fallback`.
fn number(key: &str, fallback: f32) -> f32 {
    std::env::var(key)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(fallback)
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
    lights: Query<Entity, With<DirectionalLight>>,
) {
    for mut orbit in &mut camera {
        orbit.pitch = 0.10;
        orbit.yaw = 0.4;
        // Far enough back that arches read as arches. Flying *through* the rock
        // is what the motion shows; a still has to show the form, and at close
        // range every frame is one wall.
        orbit.radius = 20.0;
    }

    // The harness's single flat key light is built for reading a HUD off a mesh.
    // Caves need contrast or every interior reads as one grey mass, so it is
    // replaced with a warm key and a cool fill.
    for light in &lights {
        commands.entity(light).despawn();
    }
    commands.spawn((
        DirectionalLight {
            illuminance: 11_000.0,
            shadow_maps_enabled: true,
            ..default()
        },
        Transform::from_xyz(30.0, 40.0, 18.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
    commands.spawn((
        DirectionalLight {
            illuminance: 2_600.0,
            shadow_maps_enabled: false,
            ..default()
        },
        Transform::from_xyz(-24.0, 8.0, -30.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));

    let thickness = thickness_for(DEFAULT_STEP);
    let volume = commands
        .spawn(
            VoxelVolume::new(layout(), CaveWorld { thickness })
                .with_extractor(Extractor::MarchingCubes),
        )
        .id();

    commands.insert_resource(Show {
        volume,
        stream: ChunkStream::new(),
        update: StreamUpdate::new(),
        // Both settable from the environment, because both are keyboard-only
        // otherwise and a capture ignores the keyboard. A recording wants a
        // slower fly and a wider stream radius than a person exploring does:
        // chunks then arrive well before they are prominent, so the frontier
        // pops in the far distance instead of in the middle of the frame.
        view: number("ISOMESH_STREAM_VIEW", 34.0),
        flying: true,
        speed: number("ISOMESH_SPEED", 5.0),
        travelled: 0.0,
        thickness,
    });
    commands.insert_resource(Look(materials.add(StandardMaterial {
        base_color: Color::srgb(0.50, 0.45, 0.40),
        perceptual_roughness: 0.92,
        // Caves are lit from one side and the far wall of a tunnel is a
        // back-face; without this they render as holes in the world.
        double_sided: true,
        cull_mode: None,
        ..default()
    })));
    commands.insert_resource(MeshBudget {
        per_frame: std::time::Duration::from_micros(3_000),
        max_in_flight: 8,
    });
}

/// Where the camera is after `d` metres.
///
/// Held near the terrain's own height rather than above it, so the path goes
/// *through* the rock instead of over it. That is the entire shot.
fn path(d: f32) -> Vec3 {
    let x = d * 0.78;
    let z = d * 0.44 + 11.0 * (d * 0.035).sin();
    // A little above the surface, so the camera skims the tops of the arches
    // rather than boring through solid rock for the whole flight.
    let y = 2.6 * (x * 0.11).sin() * (z * 0.09).cos() + 1.4 * (d * 0.06).cos() + 3.0;
    Vec3::new(x, y, z)
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut show: ResMut<Show>,
    mut commands: Commands,
    mut loaded: ResMut<Loaded>,
    volumes: Query<Entity, With<VoxelVolume>>,
) {
    if capture.is_active() {
        return;
    }
    if keys.just_pressed(KeyCode::Space) {
        show.flying = !show.flying;
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        show.speed = (show.speed + 1.0).min(16.0);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        show.speed = (show.speed - 1.0).max(1.0);
    }
    for (key, step) in [
        (KeyCode::Digit1, 0u8),
        (KeyCode::Digit2, 1),
        (KeyCode::Digit3, 2),
    ] {
        if !keys.just_pressed(key) {
            continue;
        }
        // The field is owned by the volume, so changing it means a new volume
        // and a cleared world. Cheap here, and it keeps `VoxelVolume` immutable
        // rather than adding a mutation path the plugin does not have.
        show.thickness = thickness_for(step);
        for entity in &volumes {
            commands.entity(entity).despawn();
        }
        for (_, entity) in loaded.0.drain() {
            commands.entity(entity).despawn();
        }
        let Show { stream, update, .. } = &mut *show;
        stream.clear(update);
        update.reset();
        show.volume = commands
            .spawn(
                VoxelVolume::new(
                    layout(),
                    CaveWorld {
                        thickness: show.thickness,
                    },
                )
                .with_extractor(Extractor::MarchingCubes),
            )
            .id();
    }
}

fn stream(mut commands: Commands, mut show: ResMut<Show>, mut loaded: ResMut<Loaded>) {
    let here = path(show.travelled);
    let view = show.view;
    let Ok(config) = StreamConfig::new(view, view * (1.0 + HYSTERESIS)) else {
        return;
    };
    let layout = layout();
    let volume = show.volume;
    let Show { stream, update, .. } = &mut *show;
    if stream
        .update(&layout, [here.x, here.y, here.z], &config, update)
        .is_err()
    {
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

fn attach(
    mut commands: Commands,
    look: Res<Look>,
    fresh: Query<(Entity, &ChunkMesh), Without<Mesh3d>>,
) {
    for (entity, chunk) in &fresh {
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
            MeshMaterial3d(look.0.clone()),
            Transform::default(),
        ));
    }
}

fn fly(
    time: Res<Time>,
    mut show: ResMut<Show>,
    mut camera: Query<&mut OrbitCamera>,
    mut flags: ResMut<ViewFlags>,
) {
    flags.grid = false;
    if show.flying {
        show.travelled += show.speed * time.delta_secs();
    }
    let here = path(show.travelled);
    // Aim a little ahead, so the camera looks along the flight rather than at
    // the wall it is about to pass through.
    let ahead = path(show.travelled + 6.0);
    for mut orbit in &mut camera {
        orbit.focus = ahead;
        let to = ahead - here;
        orbit.yaw = to.x.atan2(to.z);
        orbit.pitch = 0.22;
    }
}

fn hud(
    show: Res<Show>,
    loaded: Res<Loaded>,
    meshes: Res<Assets<Mesh>>,
    chunks: Query<&ChunkMesh>,
    mut stats: ResMut<DemoStats>,
) {
    let (mut vertices, mut triangles) = (0usize, 0usize);
    for chunk in &chunks {
        if let Some(mesh) = meshes.get(&chunk.0) {
            vertices += mesh.count_vertices();
            triangles += mesh.indices().map_or(0, |i| i.len() / 3);
        }
    }
    stats.title = format!(
        "E-210  showcase   {} chunks resident   {:.0} m flown",
        loaded.0.len(),
        show.travelled
    );
    stats.vertices = vertices;
    stats.triangles = triangles;
    stats.extra = vec![
        format!(
            "{:<24} {:>8.2}   gyroid iso-level, rock is g < this  [1] [2] [3]",
            "cave density", show.thickness
        ),
        format!(
            "{:<24} {:>8.0}   metres of view",
            "stream radius", show.view
        ),
        String::new(),
        "solid(p) = max( p.y - height(x,z) , gyroid(p) - level )".into(),
        String::new(),
        "a max is an intersection: rock exists where a point is BELOW the".into(),
        "surface and INSIDE the gyroid labyrinth. it is triply periodic, so it".into(),
        "tunnels in x, y and z -- which is why there are ceilings.".into(),
        String::new(),
        "none of this has a height. a heightfield stores one number per column".into(),
        "and cannot represent an arch, a cave, or the rock you are under right".into(),
        "now. that is the whole reason to reach for a voxel mesher.".into(),
        String::new(),
        format!("speed {:.0}   [ and ] to change   Space pauses", show.speed),
    ];
}
