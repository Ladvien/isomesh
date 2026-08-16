//! E-204 — the debris is the geometry that was removed, not a prop.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_destruction --release
//! ```
//!
//! **Always `--release`.**
//!
//! `Space` fire · `[` `]` charge radius · `X` reset · `1`–`3` target.
//!
//! # What "runtime fragments" has to mean
//!
//! The cheap version of this demo pre-fractures a wall into props at build time
//! and hides them until you shoot. That proves nothing about a meshing crate: the
//! fragments were authored, and the field never had to produce them.
//!
//! Here a shot does two things with the *same* signed distance field. It appends
//! `Brush::subtract(sphere)` to the wall's edit log, which craters the wall; and
//! it meshes the **intersection** of the solid-before-the-shot with that same
//! sphere, which is exactly the volume that was removed. That intersection
//! becomes the debris. Nothing is authored, and the crater and the fragment are
//! two views of one boolean.
//!
//! # The hard cases are the point, and the ticket names them
//!
//! *"Carve a spiral and a hollow shell — that's where decomposition fails."* Both
//! are here as targets, because both defeat the naive answer:
//!
//! - A **hollow shell**'s convex hull is a solid ball. Hand a physics engine the
//!   hull and every shot passes through a wall that is visibly there.
//! - A **spiral**'s hull is a fat cylinder that swallows the gaps. Anything
//!   dropped near it rests on air.
//!
//! So fragments get a **convex decomposition**, not a hull, and the HUD reports
//! how many parts each one cost. That number is the price of the correctness, and
//! it is worth seeing next to the shape that caused it.
//!
//! # What is measured
//!
//! A fragment is a *correct* physics body if it comes to rest on the world rather
//! than falling through it. So the demo counts fragments that settle against
//! fragments that end up below the floor plane, which is what tunnelling through
//! a too-thin or inside-out collider looks like from the outside.
//!
//! Decomposition time is reported per fragment, because it is the cost that
//! decides whether this is usable at all — and it is charged on the frame the
//! shot lands, which is the worst possible frame for it.

mod common;

use std::time::Instant;

use avian3d::prelude::*;
use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoStats, OrbitCamera};
use isomesh::brush::{Brush, BrushStack};
use isomesh::fields::{BoxExact, Difference, Intersection, Sphere};
use isomesh::{RuntimeShape3, Sdf};

/// Grid the *targets* are meshed on.
const TARGET_SAMPLES: u32 = 65;
/// Grid a *fragment* is meshed on. Smaller than the target's, because the box is
/// smaller — but not *much* smaller, because the debris is what this example is
/// named for. At 21 the fragments' silhouettes are visibly faceted and the demo
/// reads as low-poly rubble rather than as the boolean it is.
const FRAGMENT_SAMPLES: u32 = 41;
const HALF_EXTENT: f32 = 4.0;

/// Anything below this has fallen out of the world.
const FLOOR_Y: f32 = -8.0;

const MIN_RADIUS: f32 = 0.35;
/// How fast a fragment leaves the crater it was cut from, in world units per
/// second, and how fast it tumbles.
///
/// Chosen so the fragment separates within a frame or two at 60 Hz. Slower and
/// the coincident surfaces stay visible long enough to be captured; much faster
/// and the debris leaves the frame before it can be seen to be the right shape.
const EJECT_SPEED: f32 = 9.0;
const EJECT_SPIN: f32 = 3.0;

const MAX_RADIUS: f32 = 1.2;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Target {
    Wall,
    HollowShell,
    Spiral,
}

impl Target {
    fn name(self) -> &'static str {
        match self {
            Self::Wall => "wall",
            Self::HollowShell => "hollow shell",
            Self::Spiral => "spiral",
        }
    }

    /// Why this one is here.
    fn note(self) -> &'static str {
        match self {
            Self::Wall => "the ordinary case -- a slab, and a crater in it",
            Self::HollowShell => "hull is a solid ball, so a hull collider is a lie",
            Self::Spiral => "hull is a fat cylinder that swallows every gap",
        }
    }

    /// `ISOMESH_TARGET=wall|shell|spiral`, so the two shapes the ticket names as
    /// decomposition's failure cases can be captured and measured without a
    /// keyboard.
    fn from_env() -> Self {
        match std::env::var("ISOMESH_TARGET").unwrap_or_default().as_str() {
            "shell" | "hollow_shell" => Self::HollowShell,
            "spiral" => Self::Spiral,
            _ => Self::Wall,
        }
    }

    fn next(self) -> Self {
        match self {
            Self::Wall => Self::HollowShell,
            Self::HollowShell => Self::Spiral,
            Self::Spiral => Self::Wall,
        }
    }
}

/// A vertical slab — [`BoxExact`], which is the same distance function this
/// example used to spell out, plus the analytic gradient it did not have.
fn wall() -> BoxExact<f32> {
    BoxExact {
        center: [0.0; 3],
        half_extents: [2.6, 2.2, 0.45],
    }
}

/// A closed spherical shell: a ball with a smaller ball removed from inside.
///
/// Its convex hull is the outer ball, so a hull collider fills the cavity that
/// makes it a shell. This is the shape the ticket names first.
///
/// `Difference` is `max(fa, −fb)`, and `−(r − 1.7)` is `1.7 − r`, so this is
/// exactly the `(r - 2.3).max(1.7 - r)` it replaces — evaluated by the crate,
/// with the gradient of whichever operand is active rather than six extra
/// samples.
fn hollow_shell() -> Difference<Sphere<f32>, Sphere<f32>> {
    Difference {
        a: Sphere {
            center: [0.0; 3],
            radius: 2.3,
        },
        b: Sphere {
            center: [0.0; 3],
            radius: 1.7,
        },
    }
}

/// A helical tube.
///
/// Distance to the nearest turn of the helix, minus the tube radius. The nearest
/// turn is found by rounding rather than searched: `y` determines which turn is
/// closest, and the two neighbours either side are checked because a point
/// between turns is nearer one of them than the rounded one.
#[derive(Clone, Copy)]
struct Spiral;

impl Sdf for Spiral {
    type Scalar = f32;
    fn sample(&self, p: [f32; 3]) -> f32 {
        const COIL: f32 = 1.7;
        const PITCH: f32 = 1.1;
        const TUBE: f32 = 0.42;
        let theta = p[2].atan2(p[0]);
        let radial = (p[0] * p[0] + p[2] * p[2]).sqrt() - COIL;
        let mut best = f32::MAX;
        // `y` on the helix is PITCH * (theta + 2*pi*k) / (2*pi).
        let base = (p[1] / PITCH) - theta / std::f32::consts::TAU;
        for k in [base.floor(), base.ceil()] {
            let y = PITCH * (theta / std::f32::consts::TAU + k);
            let dy = p[1] - y;
            let d = (radial * radial + dy * dy).sqrt() - TUBE;
            if d < best {
                best = d;
            }
        }
        best
    }
}

/// One target, as a field, so the three share one code path.
#[derive(Clone, Copy)]
enum Solid {
    Wall(BoxExact<f32>),
    Shell(Difference<Sphere<f32>, Sphere<f32>>),
    Spiral(Spiral),
}

impl Sdf for Solid {
    type Scalar = f32;
    fn sample(&self, p: [f32; 3]) -> f32 {
        match self {
            Self::Wall(f) => f.sample(p),
            Self::Shell(f) => f.sample(p),
            Self::Spiral(f) => f.sample(p),
        }
    }

    /// **Forwarded, and that is the point of consuming `fields::` at all.**
    /// [`Sdf::gradient`]'s default is central differences — six extra `sample`
    /// calls per normal — so a dispatch layer that implements only `sample`
    /// throws away the analytic gradients the crate's fields carry. `Spiral` has
    /// no analytic gradient to forward and falls back to the same default it
    /// always used.
    fn gradient(&self, p: [f32; 3]) -> [f32; 3] {
        match self {
            Self::Wall(f) => f.gradient(p),
            Self::Shell(f) => f.gradient(p),
            Self::Spiral(f) => f.gradient(p),
        }
    }
}

impl From<Target> for Solid {
    fn from(t: Target) -> Self {
        match t {
            Target::Wall => Self::Wall(wall()),
            Target::HollowShell => Self::Shell(hollow_shell()),
            Target::Spiral => Self::Spiral(Spiral),
        }
    }
}

#[derive(Resource)]
struct World {
    target: Target,
    shots: Vec<Brush<Sphere<f32>>>,
    radius: f32,
    /// Fragments spawned, and what they cost.
    fragments: u32,
    convex_parts: u32,
    worst_decompose_ms: f64,
    total_decompose_ms: f64,
    /// Fragments that could not be given a collider at all.
    without_collider: u32,
    /// Frames since the last shot, for the scripted sequence.
    since_shot: f32,
    fired: u32,
}

#[derive(Resource)]
struct Look {
    solid: Handle<StandardMaterial>,
    debris: Handle<StandardMaterial>,
}

#[derive(Component)]
struct TargetMesh;

#[derive(Component)]
struct Debris;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-204 game destruction".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .add_plugins(PhysicsPlugins::default())
        .insert_resource(World {
            target: Target::from_env(),
            shots: Vec::new(),
            radius: 0.7,
            fragments: 0,
            convex_parts: 0,
            worst_decompose_ms: 0.0,
            total_decompose_ms: 0.0,
            without_collider: 0,
            since_shot: 0.0,
            fired: 0,
        })
        .add_systems(Startup, setup)
        .init_resource::<Lost>()
        .add_systems(
            Update,
            (controls, fire, remesh_target, retire, report, hud).chain(),
        )
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    for mut orbit in &mut camera {
        orbit.yaw = 0.6;
        orbit.pitch = 0.25;
        orbit.radius = 12.0;
    }
    // A floor, without which "lost" measures nothing. The first version had no
    // ground at all, so every fragment eventually fell past FLOOR_Y and the HUD
    // reported 15 of 23 as having "fallen through something" -- an accusation
    // against colliders that were working perfectly. A metric for tunnelling
    // needs a surface to tunnel through.
    commands.spawn((
        RigidBody::Static,
        Collider::cuboid(40.0, 1.0, 40.0),
        Mesh3d(meshes.add(Cuboid::new(40.0, 1.0, 40.0).mesh())),
        MeshMaterial3d(materials.add(StandardMaterial {
            base_color: Color::srgb(0.30, 0.32, 0.36),
            perceptual_roughness: 0.95,
            ..default()
        })),
        Transform::from_xyz(0.0, -4.0, 0.0),
    ));

    commands.insert_resource(Look {
        solid: materials.add(StandardMaterial {
            base_color: Color::srgb(0.66, 0.68, 0.72),
            perceptual_roughness: 0.7,
            ..default()
        }),
        debris: materials.add(StandardMaterial {
            base_color: Color::srgb(0.90, 0.48, 0.28),
            perceptual_roughness: 0.55,
            ..default()
        }),
    });
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut world: ResMut<World>,
    mut commands: Commands,
    debris: Query<Entity, With<Debris>>,
) {
    if capture.is_active() {
        return;
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        world.radius = (world.radius + 0.1).min(MAX_RADIUS);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        world.radius = (world.radius - 0.1).max(MIN_RADIUS);
    }
    for (key, target) in [
        (KeyCode::Digit1, Target::Wall),
        (KeyCode::Digit2, Target::HollowShell),
        (KeyCode::Digit3, Target::Spiral),
    ] {
        if keys.just_pressed(key) {
            world.target = target;
            reset(&mut world, &mut commands, &debris);
        }
    }
    if keys.just_pressed(KeyCode::KeyT) {
        world.target = world.target.next();
        reset(&mut world, &mut commands, &debris);
    }
    if keys.just_pressed(KeyCode::KeyX) {
        reset(&mut world, &mut commands, &debris);
    }
}

fn reset(world: &mut World, commands: &mut Commands, debris: &Query<Entity, With<Debris>>) {
    world.shots.clear();
    world.fragments = 0;
    world.convex_parts = 0;
    world.worst_decompose_ms = 0.0;
    world.total_decompose_ms = 0.0;
    world.without_collider = 0;
    world.fired = 0;
    for entity in debris {
        commands.entity(entity).despawn();
    }
}

/// Where shot `n` lands.
///
/// Scripted rather than aimed, so a capture and a measured run are the same run.
/// The points walk a lissajous across the target's face, which spreads impacts
/// over the shape instead of drilling one hole.
fn impact_point(n: u32, target: Target) -> Vec3 {
    let t = n as f32;
    let (u, v) = ((t * 1.1).sin(), (t * 0.7).cos());
    match target {
        Target::Wall => Vec3::new(u * 1.9, v * 1.5, 0.0),
        Target::HollowShell => Vec3::new(u * 1.6, v * 1.6, (t * 0.4).sin() * 1.2),
        Target::Spiral => Vec3::new(u * 1.7, v * 1.4, (t * 0.9).cos() * 1.7),
    }
}

/// Fire on a timer, carve the solid, and turn what was removed into a body.
#[allow(clippy::too_many_arguments)]
fn fire(
    time: Res<Time>,
    keys: Res<ButtonInput<KeyCode>>,
    mut world: ResMut<World>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    look: Res<Look>,
    camera: Query<&GlobalTransform, With<Camera3d>>,
) {
    world.since_shot += time.delta_secs();
    let manual = keys.just_pressed(KeyCode::Space);
    // A shot every 0.9 s, so the demo is doing something whenever it is looked
    // at and a capture never lands on an untouched wall.
    if !manual && (world.since_shot < 0.9 || world.fired >= 24) {
        return;
    }
    world.since_shot = 0.0;

    let centre = impact_point(world.fired, world.target);
    world.fired += 1;

    // The direction the charge was travelling, and the reason the fragment has
    // to be given one.
    //
    // The fragment is `solid ∩ charge` and the crater is `solid − charge`, so
    // the two share a surface *exactly*: the same sphere, extracted twice from
    // two different fields. A fragment left at rest sits in the hole it came
    // from with its faces coincident against the crater's, which renders as a
    // chunk fused into the wall rather than knocked out of it.
    //
    // **The obvious direction is the wrong one, measured.** The first attempt
    // used the solid's own gradient at the impact, which is "away from the
    // material" and reads correctly. It returns exactly `[0, 0, 0]`: the wall
    // is a slab and `impact_point` puts every charge at `z = 0`, which is the
    // slab's **medial plane** -- equidistant from both faces, where a distance
    // field's gradient vanishes by definition. Where it was non-zero it pointed
    // `+Y`, *along* the wall rather than out of it, because the nearest surface
    // to that charge was the top edge.
    //
    // The charge's own line of travel has neither problem. It is defined
    // wherever the camera is, it never lands on a medial axis, and it is what
    // actually determines where debris goes.
    let Ok(eye) = camera.single() else {
        return;
    };
    let shot = (centre - eye.translation()).normalize_or_zero();
    // Debris comes back out the way the charge went in. That is spall, it is
    // what an impact actually throws, and here it is also the only direction
    // that keeps the fragment visible rather than hidden behind the target.
    let spall = -shot;
    let radius = world.radius;
    let shape = Sphere {
        center: [centre.x, centre.y, centre.z],
        radius,
    };

    // The solid as it stands *before* this shot. The fragment is the part of it
    // inside the charge, so it has to be sampled before the log grows.
    let base = Solid::from(world.target);
    let before = BrushStack {
        base,
        brushes: &world.shots,
    };

    // Mesh the intersection over a box that just contains the charge. One cell
    // of padding so a crossing exactly on the boundary is not clipped.
    let pad = radius * 1.25;
    let min = [centre.x - pad, centre.y - pad, centre.z - pad];
    let cell = (pad * 2.0) / (FRAGMENT_SAMPLES - 1) as f32;
    let mut fragment = MeshBuilder::new();
    let carved = Intersection {
        a: before,
        b: shape,
    };
    let Ok(shape3) = RuntimeShape3::new([FRAGMENT_SAMPLES; 3]) else {
        return;
    };
    if isomesh::marching_cubes::MarchingCubes::<f32>::new()
        .extract(&carved, &shape3, min, cell, &mut fragment)
        .is_err()
    {
        return;
    }

    world.shots.push(Brush::subtract(shape));

    // An empty intersection means the charge missed the solid. That is a miss,
    // not a failure, and it must not count as a fragment.
    if fragment.triangle_count() == 0 {
        return;
    }

    // Re-centre on the fragment's own centroid so the rigid body spins about
    // itself rather than about the world origin.
    let positions = fragment.positions();
    let mut centroid = Vec3::ZERO;
    for p in positions {
        centroid += Vec3::from(*p);
    }
    centroid /= positions.len() as f32;
    let local: Vec<[f32; 3]> = positions
        .iter()
        .map(|p| [p[0] - centroid.x, p[1] - centroid.y, p[2] - centroid.z])
        .collect();

    let mut mesh = fragment.into_mesh();
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, local);

    // The decomposition, and its cost. Charged on the frame the shot lands,
    // which is the worst frame it could be charged on.
    let started = Instant::now();
    let collider = Collider::convex_decomposition_from_mesh(&mesh);
    let ms = started.elapsed().as_secs_f64() * 1000.0;
    world.total_decompose_ms += ms;
    if ms > world.worst_decompose_ms {
        world.worst_decompose_ms = ms;
    }

    world.fragments += 1;
    let handle = meshes.add(mesh);
    match collider {
        Some(collider) => {
            world.convex_parts += convex_parts(&collider);
            commands.spawn((
                RigidBody::Dynamic,
                collider,
                Mesh3d(handle),
                MeshMaterial3d(look.debris.clone()),
                // Placed at the **mouth** of the crater rather than inside
                // it, and this is the part that actually fixes the artefact.
                //
                // Velocity alone was not enough: measured, some fragments leave
                // at 9 m/s and others sit at the impact with a velocity of
                // 0.004 and no response to gravity -- mass is a healthy 1.19
                // and sleeping is off, so the cause is inside the solver and
                // not somewhere this example can reach. A fragment that fails
                // to move must therefore fail *clear of the wall*, because
                // sitting still inside the crater is the one place its faces
                // are coincident with the crater's by construction.
                Transform::from_translation(centroid + spall * radius),
                // The spin is not decoration: a fragment that translates
                // without rotating reads as a prop being slid, and this example
                // exists to say the debris is the geometry that was removed.
                LinearVelocity(spall * EJECT_SPEED),
                AngularVelocity(spall.cross(Vec3::Y) * EJECT_SPIN),
                // Never let a fragment fall asleep.
                //
                // A sleeping body in avian ignores gravity, and the target has
                // no collider at all -- it is geometry, not a body -- so a
                // fragment whose velocity is killed at spawn (by contact with
                // an older fragment that has not cleared yet) freezes exactly
                // where it was cut and stays there. Measured: debris resting at
                // y = 1.498 with velocity 0.004, which nothing but sleep can
                // explain. That frozen fragment is the artefact -- its faces
                // are coincident with the crater's by construction, so a
                // fragment that never leaves reads as a chunk fused into the
                // wall.
                SleepingDisabled,
                Debris,
            ));
        }
        None => {
            // Reported rather than substituted. A convex hull here would make
            // the demo look like it worked and quietly be the wrong shape --
            // which is exactly the failure this example exists to show.
            world.without_collider += 1;
            commands.spawn((
                Mesh3d(handle),
                MeshMaterial3d(look.debris.clone()),
                Transform::from_translation(centroid + spall * radius),
                Debris,
            ));
            // No RigidBody on this arm, so it has nothing to give a velocity
            // to -- a fragment with no collider stays where it was carved. That
            // is the honest rendering of the failure, and `without_collider` on
            // the HUD is where it is counted.
        }
    }
}

/// How many convex parts a compound collider ended up with.
fn convex_parts(collider: &Collider) -> u32 {
    collider
        .shape()
        .as_compound()
        .map_or(1, |compound| compound.shapes().len() as u32)
}

/// Re-mesh the target whenever the edit log changes.
#[allow(clippy::too_many_arguments)]
fn remesh_target(
    world: Res<World>,
    look: Res<Look>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    existing: Query<Entity, With<TargetMesh>>,
    mut last: Local<Option<(usize, u8)>>,
) {
    let key = (
        world.shots.len(),
        match world.target {
            Target::Wall => 0u8,
            Target::HollowShell => 1,
            Target::Spiral => 2,
        },
    );
    if *last == Some(key) {
        return;
    }
    *last = Some(key);

    let base = Solid::from(world.target);
    let field = BrushStack {
        base,
        brushes: &world.shots,
    };
    let Ok(shape) = RuntimeShape3::new([TARGET_SAMPLES; 3]) else {
        return;
    };
    let cell = (HALF_EXTENT * 2.0) / (TARGET_SAMPLES - 1) as f32;
    let mut builder = MeshBuilder::new();
    if isomesh::marching_cubes::MarchingCubes::<f32>::new()
        .extract(&field, &shape, [-HALF_EXTENT; 3], cell, &mut builder)
        .is_err()
    {
        return;
    }
    if builder.triangle_count() == 0 {
        return;
    }

    let mesh = builder.into_mesh();
    // The target is static, and it gets a decomposition too -- a hull would let
    // debris rest inside the shell's cavity.
    let collider = Collider::convex_decomposition_from_mesh(&mesh);
    let handle = meshes.add(mesh);
    for entity in &existing {
        commands.entity(entity).despawn();
    }
    let mut entity = commands.spawn((
        Mesh3d(handle),
        MeshMaterial3d(look.solid.clone()),
        Transform::default(),
        TargetMesh,
    ));
    if let Some(collider) = collider {
        entity.insert((RigidBody::Static, collider));
    }
}

/// How many fragments have left the world.
#[derive(Resource, Default)]
struct Lost(u32);

/// Count what fell out of the world, and stop simulating it.
fn retire(
    mut commands: Commands,
    mut lost: ResMut<Lost>,
    debris: Query<(Entity, &Transform), With<Debris>>,
) {
    for (entity, transform) in &debris {
        if transform.translation.y < FLOOR_Y {
            lost.0 += 1;
            commands.entity(entity).despawn();
        }
    }
}

/// One CSV row per fragment, so the decomposition cost can be read over a long
/// run rather than off the frame a screenshot lands on.
fn report(world: Res<World>, lost: Res<Lost>, mut last: Local<u32>) {
    if world.fragments == *last {
        return;
    }
    *last = world.fragments;
    info!(
        "destruction,{},{},{},{},{:.2},{:.2},{},{}",
        world.target.name(),
        world.fragments,
        world.convex_parts,
        world.without_collider,
        world.total_decompose_ms / f64::from(world.fragments.max(1)),
        world.worst_decompose_ms,
        lost.0,
        world.shots.len(),
    );
}

fn hud(
    world: Res<World>,
    lost: Res<Lost>,
    meshes: Res<Assets<Mesh>>,
    target: Query<&Mesh3d, With<TargetMesh>>,
    mut stats: ResMut<DemoStats>,
    debris: Query<&Transform, With<Debris>>,
) {
    let alive = debris.iter().count();
    let fell = lost.0;
    // The target mesh, which is what "vertices" and "triangles" mean everywhere
    // else in this harness. The first draft put the debris count in those rows
    // and the HUD read "2 vertices" beside a wall with thousands.
    let (mut vertices, mut triangles) = (0usize, 0usize);
    for handle in &target {
        if let Some(mesh) = meshes.get(&handle.0) {
            vertices += mesh.count_vertices();
            triangles += mesh.indices().map_or(0, |i| i.len() / 3);
        }
    }
    let mean_parts = if world.fragments > 0 {
        f64::from(world.convex_parts) / f64::from(world.fragments)
    } else {
        0.0
    };
    let mean_ms = if world.fragments > 0 {
        world.total_decompose_ms / f64::from(world.fragments)
    } else {
        0.0
    };

    let verdict = if world.without_collider > 0 {
        format!(
            "!! {} fragment(s) got NO collider -- reported, not replaced by a hull",
            world.without_collider
        )
    } else if fell > 0 {
        format!("!! {fell} fragment(s) went through the floor -- a collider is wrong")
    } else if world.fragments == 0 {
        "no fragments yet".to_string()
    } else {
        "every fragment is a body, and every one of them is still on the world".to_string()
    };

    stats.title = format!(
        "E-204  destruction   target: {} ({})",
        world.target.name(),
        world.target.note()
    );
    stats.vertices = vertices;
    stats.triangles = triangles;
    stats.extra = vec![
        format!(
            "{:<26} {:>8}   shots that hit the solid",
            "fragments", world.fragments
        ),
        format!("{:<26} {:>8}   still simulating", "debris alive", alive),
        format!("{:<26} {:>8}   fell out of the world", "lost", fell),
        String::new(),
        "the fragment is the INTERSECTION of the solid with the charge -- the same".into(),
        "boolean that made the crater, meshed. nothing here is pre-fractured.".into(),
        String::new(),
        format!(
            "{:<26} {:>8.1}   convex parts per fragment, mean",
            "decomposition", mean_parts
        ),
        format!(
            "{:<26} {:>8}   convex parts in total",
            "", world.convex_parts
        ),
        format!("{:<26} {:>8.2}   ms per fragment, mean", "cost", mean_ms),
        format!(
            "{:<26} {:>8.2}   ms worst single fragment",
            "", world.worst_decompose_ms
        ),
        String::new(),
        "a convex HULL would be cheaper and wrong: a hollow shell's hull is a".into(),
        "solid ball, and a spiral's is a fat cylinder that swallows every gap.".into(),
        "press 2 and 3 -- those are the shapes the ticket names.".into(),
        String::new(),
        verdict,
        String::new(),
        format!(
            "charge radius {:.2}   [ and ] to change   Space fires   T target   X resets",
            world.radius
        ),
    ];
}
