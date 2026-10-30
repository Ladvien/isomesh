//! Throwaway smoke for issue #10, iteration 2. Never committed.
//!
//! Headless: the camera renders to an image, so the result does not depend on a
//! window being visible. Three ways a game meets an empty mesh, each with a
//! `Mesh3d` attached:
//! 1. the plugin: 16 chunks over a sphere at x = -2, eight of them empty;
//! 2. a direct `meshes.add(to_bevy_mesh(&empty))`;
//! 3. a drawn sphere at x = +2 replaced IN PLACE (`*meshes.get_mut(h) = ...`)
//!    by an empty mesh -- it must vanish, not linger.
//!
//! Phases are wall-clock, not frame counts: on a throttled host pipelines take
//! a long time to compile, and a frame-count screenshot photographs nothing.

use std::time::{Duration, Instant};

use bevy::app::ScheduleRunnerPlugin;
use bevy::camera::RenderTarget;
use bevy::prelude::*;
use bevy::render::render_resource::{TextureFormat, TextureUsages};
use bevy::render::view::screenshot::{Screenshot, save_to_disk};
use bevy::window::ExitCondition;
use bevy::winit::WinitPlugin;
use bevy_isomesh::{ChunkMesh, IsomeshPlugin, NeedsRemesh, VoxelChunk, VoxelVolume, to_bevy_mesh};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::Sphere;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3};

#[derive(Resource)]
struct InPlace(Handle<Mesh>);

#[derive(Resource)]
struct Target(Handle<Image>);

#[derive(Resource)]
struct Smoke {
    dir: String,
    started: Instant,
    frame: u32,
    phase: u32,
    mark: Instant,
    mark_frame: u32,
    settle: Duration,
}

#[derive(Resource)]
struct Mat(Handle<StandardMaterial>);

fn main() {
    let dir = std::env::var("SMOKE_DIR").unwrap_or_else(|_| "/tmp/i2/smoke".into());
    let settle = std::env::var("SMOKE_SETTLE_S")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(60);
    let _ = std::fs::create_dir_all(&dir);
    let now = Instant::now();
    App::new()
        .add_plugins(
            DefaultPlugins
                .set(WindowPlugin {
                    primary_window: None,
                    exit_condition: ExitCondition::DontExit,
                    ..default()
                })
                .disable::<WinitPlugin>(),
        )
        .add_plugins(ScheduleRunnerPlugin::run_loop(Duration::from_secs_f64(
            1.0 / 60.0,
        )))
        .add_plugins(IsomeshPlugin)
        .insert_resource(Smoke {
            dir,
            started: now,
            frame: 0,
            phase: 0,
            mark: now,
            mark_frame: 0,
            settle: Duration::from_secs(settle),
        })
        .add_systems(Startup, setup)
        .add_systems(Update, (attach, drive).chain())
        .run();
}

fn sphere_buffer(center: [f32; 3]) -> MeshBuffer<f32> {
    let mut buffer = MeshBuffer::<f32>::new();
    let shape = RuntimeShape3::new([33; 3]).expect("shape");
    let origin = [center[0] - 1.0, center[1] - 1.0, center[2] - 1.0];
    MarchingCubes::<f32>::new()
        .extract(
            &Sphere::<f32> { center, radius: 0.8 },
            &shape,
            origin,
            2.0 / 32.0,
            &mut buffer,
        )
        .expect("extract");
    buffer
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut images: ResMut<Assets<Image>>,
) {
    let mut target = Image::new_target_texture(960, 540, TextureFormat::Rgba8UnormSrgb, None);
    target.texture_descriptor.usage |= TextureUsages::COPY_SRC;
    let target = images.add(target);
    commands.insert_resource(Target(target.clone()));

    commands.spawn((
        Camera3d::default(),
        RenderTarget::Image(target.into()),
        Transform::from_xyz(0.0, 2.5, 6.5).looking_at(Vec3::ZERO, Vec3::Y),
    ));
    commands.spawn((
        DirectionalLight::default(),
        Transform::from_xyz(4.0, 8.0, 4.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
    let material = materials.add(Color::srgb(0.8, 0.7, 0.6));
    commands.insert_resource(Mat(material.clone()));

    // 1. The plugin: chunk x in -3..=0 spans world x in [-3, 1); the sphere
    //    covers x in [-2.8, -1.2], so chunks x = -1 and x = 0 are empty.
    let layout = ChunkLayout::new(16, 1.0 / 16.0, [0.0; 3]).expect("layout");
    let volume = commands
        .spawn(VoxelVolume::new(
            layout,
            Sphere::<f32> {
                center: [-2.0, 0.0, 0.0],
                radius: 0.8,
            },
        ))
        .id();
    for z in -1..1 {
        for y in -1..1 {
            for x in -3..1 {
                commands.spawn((
                    VoxelChunk {
                        id: ChunkId::new([x, y, z]),
                        volume,
                    },
                    NeedsRemesh,
                ));
            }
        }
    }

    // 2. A direct empty mesh, drawn by a Mesh3d.
    commands.spawn((
        Mesh3d(meshes.add(to_bevy_mesh(&MeshBuffer::<f32>::new()))),
        MeshMaterial3d(material.clone()),
    ));

    // 3. A drawn sphere at x = +2, to be emptied in place.
    let handle = meshes.add(to_bevy_mesh(&sphere_buffer([2.0, 0.0, 0.0])));
    commands.spawn((Mesh3d(handle.clone()), MeshMaterial3d(material)));
    commands.insert_resource(InPlace(handle));
}

fn attach(
    mut commands: Commands,
    material: Res<Mat>,
    fresh: Query<(Entity, &ChunkMesh), Without<Mesh3d>>,
) {
    for (entity, mesh) in &fresh {
        commands
            .entity(entity)
            .try_insert((Mesh3d(mesh.0.clone()), MeshMaterial3d(material.0.clone())));
    }
}

#[allow(clippy::too_many_arguments)]
fn drive(
    mut smoke: ResMut<Smoke>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    in_place: Res<InPlace>,
    target: Res<Target>,
    chunks: Query<&ChunkMesh>,
    mut exit: MessageWriter<AppExit>,
) {
    smoke.frame += 1;
    let since = smoke.mark.elapsed();
    let frames = smoke.frame - smoke.mark_frame;
    let ready = since >= smoke.settle && frames >= 30;
    let shot = |commands: &mut Commands, path: String| {
        commands
            .spawn(Screenshot::image(target.0.clone()))
            .observe(save_to_disk(path));
    };
    match smoke.phase {
        0 if ready && chunks.iter().count() == 16 => {
            shot(&mut commands, format!("{}/a_before_swap.png", smoke.dir));
            info!(
                "smoke: shot A at frame {} after {:.1}s",
                smoke.frame,
                smoke.started.elapsed().as_secs_f64()
            );
            smoke.phase = 1;
            smoke.mark = Instant::now();
            smoke.mark_frame = smoke.frame;
        }
        1 if since >= Duration::from_secs(5) && frames >= 10 => {
            if let Some(mut mesh) = meshes.get_mut(&in_place.0) {
                *mesh = to_bevy_mesh(&MeshBuffer::<f32>::new());
            }
            info!("smoke: in-place swap to an empty mesh at frame {}", smoke.frame);
            smoke.phase = 2;
            smoke.mark = Instant::now();
            smoke.mark_frame = smoke.frame;
        }
        2 if since >= Duration::from_secs(20) && frames >= 30 => {
            shot(&mut commands, format!("{}/b_after_swap.png", smoke.dir));
            info!("smoke: shot B at frame {}", smoke.frame);
            smoke.phase = 3;
            smoke.mark = Instant::now();
            smoke.mark_frame = smoke.frame;
        }
        3 if since >= Duration::from_secs(10) && frames >= 10 => {
            let mut total = 0;
            let mut empty = 0;
            let mut empty_usages = Vec::new();
            for chunk in &chunks {
                total += 1;
                if let Some(mesh) = meshes.get(&chunk.0)
                    && mesh.count_vertices() == 0
                {
                    empty += 1;
                    empty_usages.push(format!("{:?}", mesh.asset_usage));
                }
            }
            empty_usages.sort();
            empty_usages.dedup();
            let swapped = meshes
                .get(&in_place.0)
                .map(|m| (m.count_vertices(), format!("{:?}", m.asset_usage)));
            info!(
                "smoke: {total} chunk meshes, {empty} empty, empty usages {empty_usages:?}; in-place mesh now {swapped:?}; {} frames in {:.1}s",
                smoke.frame,
                smoke.started.elapsed().as_secs_f64()
            );
            exit.write(AppExit::Success);
        }
        _ => {}
    }
}
