//! **The shortest path from an SDF to a mesh on screen. Copy this file.**
//!
//! ```bash
//! cargo run --example quickstart --release
//! ```
//!
//! Every other example in this directory teaches the library's *internals* — a
//! decider firing, a seam holding, a budget draining, a tunnel meshing as a
//! tunnel. Useful once you are using the crate, and the wrong first thing to
//! read. This one teaches nothing. It puts a sphere on screen and stops, so the
//! shape of a working app is visible in one file with nothing to filter out.
//!
//! # Its relationship to the README
//!
//! `README.md` carries the same sequence as a `no_run` doctest, so the API calls
//! below are compile-checked on every `cargo test` and cannot rot silently. What
//! a doctest cannot have is a camera, a light and a window — without those the
//! snippet is correct and shows you nothing. **This file is that snippet made
//! runnable**, and the difference between them is exactly the three lines in
//! `setup` that Bevy needs and `bevy_isomesh` does not provide.
//!
//! # The one boundary worth noticing
//!
//! The plugin stops at a [`ChunkMesh`], which holds a `Handle<Mesh>`. Attaching
//! `Mesh3d` is the application's line to write, and `attach` below is it. That
//! boundary is deliberate: it is what lets a headless consumer — a server
//! building colliders, a bake step — depend on this crate without ever compiling
//! a renderer.

use bevy::prelude::*;
use bevy_isomesh::{ChunkMesh, IsomeshPlugin, NeedsRemesh, VoxelChunk, VoxelVolume};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::Sphere;

fn main() -> Result<(), isomesh::Error> {
    let mut app = App::new();
    app.add_plugins(DefaultPlugins).add_plugins(IsomeshPlugin);

    // 16 *cells* per chunk axis at 1/16 of a world unit each, so one chunk spans
    // exactly 1.0 and the eight below tile the cube from -1 to 1. `Sphere` is
    // radius 1 at the origin, so it fills that cube and its surface crosses every
    // one of the eight — which is the point: they are meshed independently and
    // the result is still watertight.
    let layout = ChunkLayout::new(16, 1.0 / 16.0, [0.0; 3])?;
    let volume = app
        .world_mut()
        .spawn(VoxelVolume::new(layout, Sphere::<f32>::canonical()))
        .id();

    for z in -1..1 {
        for y in -1..1 {
            for x in -1..1 {
                app.world_mut().spawn((
                    VoxelChunk {
                        id: ChunkId::new([x, y, z]),
                        volume,
                    },
                    NeedsRemesh,
                ));
            }
        }
    }

    app.add_systems(Startup, setup);
    app.add_systems(Update, attach);
    app.run();
    Ok(())
}

/// The camera and light Bevy needs. Nothing here is isomesh's.
fn setup(mut commands: Commands) {
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(3.0, 2.5, 3.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
    commands.spawn((
        DirectionalLight::default(),
        Transform::from_xyz(4.0, 8.0, 4.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
}

/// Give each finished chunk mesh something to render with.
///
/// `Without<Mesh3d>` is what keeps this from re-inserting every frame; a chunk
/// that is re-meshed after an edit gets a new [`ChunkMesh`] and comes back
/// through here once.
fn attach(
    meshes: Query<(Entity, &ChunkMesh), Without<Mesh3d>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut commands: Commands,
) {
    for (entity, chunk_mesh) in &meshes {
        commands.entity(entity).insert((
            Mesh3d(chunk_mesh.0.clone()),
            MeshMaterial3d(materials.add(Color::srgb(0.8, 0.7, 0.6))),
        ));
    }
}
