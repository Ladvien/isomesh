//! E-208 — paint that survives destruction.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_paint --release
//! ```
//!
//! `WASD`/`QE` to move, mouse to look, **left click sprays**, **right click
//! blows a hole**. `1`–`5` pick a colour, `[` and `]` size the nozzle, `Tab`
//! releases the cursor, `X` clears the log.
//!
//! # The claim being demonstrated, and why it is an equality
//!
//! Row 4 of `docs/research/2026-08-11-novel-gameplay-opportunities.md`: *"you
//! spray graffiti on a wall, then blow a hole through it, and the paint on the
//! remaining wall is still exactly where you sprayed it — not smeared, not
//! reset."*
//!
//! That row prices the feature on **L²-nearest attribute transfer over a common
//! subdivision** — machinery for carrying per-vertex data from an old mesh to a
//! new one. This example does not contain any, because this crate does not need
//! any. A world here is a base field plus an ordered log of edits, so paint goes
//! *in the log* ([`isomesh::paint`]) and is a function of world position. The
//! carve moves the surface; the paint was never on the surface.
//!
//! So the drift is not small, it is **zero**, and the HUD says so continuously
//! rather than in a commit message. A transfer-based implementation could only
//! ever report a tolerance. See M-137.
//!
//! # A spray dirties chunks; it does not dirty the *field*
//!
//! This is the one place the example needed something `game_dig` did not.
//! [`mark_edit`] finds chunks to re-mesh by comparing two fields, and a spray
//! leaves the field bit-identical — [`PaintStack::sample`] skips sprays without
//! evaluating them. Ask `mark_edit` about a spray and it correctly answers
//! "nothing moved."
//!
//! So the two edits take different routes to the same dirty set, and the split
//! is real rather than incidental:
//!
//! | edit | what changed | how its chunks are found |
//! |---|---|---|
//! | carve | geometry | [`mark_edit`], exactly as `game_dig` does |
//! | spray | vertex colour only | the nozzle's bounding box, re-shaded in place |
//!
//! # The spacing is a power of two, for M-32's reason
//!
//! `h = 0.125`. Two chunks agree on their shared sample plane bit-for-bit only
//! at a power-of-two cell size, and this example keeps each chunk as its own
//! `Mesh3d` without a weld, so it uses the spacing where the seam is exact.

mod common;

use std::time::Instant;

use bevy::input::mouse::AccumulatedMouseMotion;
use bevy::prelude::*;
use bevy::window::{CursorGrabMode, CursorOptions, PrimaryWindow};
use bevy_isomesh::MeshBuilder;
use common::{CommonPlugin, DemoStats, OrbitCamera};
use isomesh::Sdf;
use isomesh::brush::Brush;
use isomesh::chunk::dirty::{DirtySet, mark_edit};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::{BoxExact, Sphere};
use isomesh::paint::{Edit, PaintStack, Splat};

/// Chunk edge, in cells.
const CHUNK_CELLS: u32 = 16;
/// A power of two — see the module docs.
const CELL_SIZE: f32 = 0.125;
/// Chunks along x, y, z. One deep: this is a wall, not a landscape.
const EXTENT: [i32; 3] = [3, 2, 1];

/// Unpainted concrete, as authored — see [`linear`].
const BACKGROUND_SRGB: [f32; 4] = [0.58, 0.56, 0.53, 1.0];

/// The nozzle palette, `1`–`5`, as authored.
const PALETTE: [(&str, [f32; 4]); 5] = [
    ("red", [0.85, 0.12, 0.14, 1.0]),
    ("yellow", [0.95, 0.78, 0.10, 1.0]),
    ("cyan", [0.10, 0.72, 0.85, 1.0]),
    ("violet", [0.55, 0.22, 0.80, 1.0]),
    ("black", [0.05, 0.05, 0.06, 1.0]),
];

/// The shape of a splat's alias, written once because it appears in the
/// resource, the field and every helper.
type WorldEdit = Edit<Sphere<f32>, Sphere<f32>, f32>;

/// Colours as a human picks them (sRGB) into the values the renderer wants.
///
/// [`Mesh::ATTRIBUTE_COLOR`] is **linear** RGBA, and the difference is not
/// subtle: `[0.85, 0.12, 0.14]` fed in raw renders as pale pink rather than
/// red. Converting here rather than in [`isomesh::paint`] is deliberate — the
/// core module interpolates numbers and has no business knowing what colour
/// space they are in, and blending in linear is the correct thing to do anyway.
fn linear(srgb: [f32; 4]) -> [f32; 4] {
    Color::srgba(srgb[0], srgb[1], srgb[2], srgb[3])
        .to_linear()
        .to_f32_array()
}

/// The unpainted wall colour, ready for the attribute.
fn background() -> [f32; 4] {
    linear(BACKGROUND_SRGB)
}

/// The wall, before anything happened to it.
fn wall() -> BoxExact<f32> {
    BoxExact {
        center: [0.0, 0.0, 0.0],
        half_extents: [2.6, 1.6, 0.35],
    }
}

/// A point whose colour was recorded at spray time, and what it was.
///
/// The instrument behind the HUD's drift readout. These are world positions,
/// not mesh vertices, which is the whole point: a mesh vertex is a thing the
/// carve can delete, and a world position is not.
struct Probe {
    at: Vec3,
    color: [f32; 4],
}

#[derive(Resource)]
struct World {
    layout: ChunkLayout<f32>,
    edits: Vec<WorldEdit>,
    radius: f32,
    ink: usize,
    probes: Vec<Probe>,
    /// Largest per-channel colour change at any probe, over the whole session.
    worst_drift: f32,
    /// Probes still within a nozzle-width of solid surface.
    probes_on_surface: usize,
    last_ms: f64,
    last_chunks: usize,
    last_action: &'static str,
    grabbed: bool,
}

impl World {
    /// The field and its colours: everything sprayed and carved so far.
    fn field(&self) -> PaintStack<'_, BoxExact<f32>, Sphere<f32>, Sphere<f32>> {
        PaintStack {
            base: wall(),
            edits: &self.edits,
            background: background(),
        }
    }
}

#[derive(Component)]
struct Chunk(ChunkId);

#[derive(Resource)]
struct Look {
    yaw: f32,
    pitch: f32,
}

#[derive(Resource)]
struct SurfaceMaterial(Handle<StandardMaterial>);

/// A scripted spray-then-carve, driven by `ISOMESH_AUTOPAINT`.
///
/// The acceptance is about what a carve does to existing paint, and a
/// screenshot cannot click. Without this the example could be committed
/// compiling, rendering, and never having sprayed anything — so the sequence
/// runs itself through exactly the code path a click takes.
#[derive(Resource, Default)]
struct AutoPaint {
    remaining: u32,
    total: u32,
    step: u32,
}

impl AutoPaint {
    fn from_env() -> Self {
        let total = std::env::var("ISOMESH_AUTOPAINT")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0);
        Self {
            remaining: total,
            total,
            step: 0,
        }
    }

    /// Where spray `n` lands: a figure that fills the wall rather than a line,
    /// so a long run paints different chunks instead of the same one over and
    /// over — which is what makes the growing-log timings mean anything.
    fn spray_at(n: u32) -> Vec3 {
        let t = n as f32;
        Vec3::new(2.2 * (t * 0.9).sin(), 1.2 * (t * 0.55).cos(), 0.35)
    }

    /// Step `n` of a script of `total` steps: sprays, then two carves.
    ///
    /// The carves are **centred on patches that were definitely painted**,
    /// which is the only arrangement that tests anything — a hole that misses
    /// the paint proves nothing about whether the paint moved. They take their
    /// centres from [`spray_at`](Self::spray_at) rather than hard-coded
    /// coordinates, so changing the spray pattern cannot silently move the
    /// holes off the paint.
    ///
    /// Returns the world point, whether it is a spray, and the nozzle radius.
    fn action(n: u32, total: u32) -> (Vec3, bool, f32) {
        let sprays = total.saturating_sub(2);
        if n < sprays {
            (Self::spray_at(n), true, 0.38)
        } else {
            let target = if n == sprays { 0 } else { sprays / 2 };
            let at = Self::spray_at(target);
            // Radius 0.5 from just inside the front face reaches through the
            // wall's 0.7 of thickness and takes the probe ring with it.
            (Vec3::new(at.x, at.y, 0.1), false, 0.5)
        }
    }
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-208 game paint".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Look {
            yaw: 0.0,
            pitch: 0.0,
        })
        .insert_resource(AutoPaint::from_env())
        .add_systems(Startup, setup)
        .add_systems(Update, (grab, fly, edit, report))
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
    camera: Query<Entity, With<OrbitCamera>>,
) {
    // Same move as `game_dig`: take `OrbitCamera` off rather than despawning,
    // so the orbit system skips this camera and the rest of the harness still
    // finds one.
    for entity in &camera {
        commands
            .entity(entity)
            .remove::<OrbitCamera>()
            .insert(Transform::from_xyz(0.0, 0.0, 6.4));
    }

    let layout =
        ChunkLayout::<f32>::new(CHUNK_CELLS, CELL_SIZE, [-3.0, -2.0, -1.0]).expect("valid layout");

    // White base colour, because the vertex colour is multiplied into it and a
    // tinted material would silently recolour the paint.
    let material = materials.add(StandardMaterial {
        base_color: Color::WHITE,
        perceptual_roughness: 0.9,
        ..default()
    });

    let world = World {
        layout,
        edits: Vec::new(),
        radius: 0.38,
        ink: 0,
        probes: Vec::new(),
        worst_drift: 0.0,
        probes_on_surface: 0,
        last_ms: 0.0,
        last_chunks: 0,
        last_action: "nothing yet",
        grabbed: false,
    };

    let mut dirty = DirtySet::new();
    for z in 0..EXTENT[2] {
        for y in 0..EXTENT[1] {
            for x in 0..EXTENT[0] {
                dirty.insert(ChunkId::new([x, y, z]));
            }
        }
    }
    let field = world.field();
    let layout = world.layout;
    let mut scratch = Vec::new();
    dirty.mesh_dirty(&layout, |id, origin| {
        if let Some(mesh) = mesh_chunk(&layout, &field, origin, &mut scratch) {
            commands.spawn((
                Mesh3d(meshes.add(mesh)),
                MeshMaterial3d(material.clone()),
                Chunk(id),
            ));
        }
    });

    commands.insert_resource(SurfaceMaterial(material));
    commands.insert_resource(world);
}

/// Extract one chunk and colour it.
///
/// The two passes are separate on purpose and both write into arrays the `Mesh`
/// ends up owning: [`isomesh::paint::shade`] fills `scratch`, which is then
/// swapped into the builder rather than copied. `scratch` is the caller's, so
/// meshing every chunk in the volume allocates once.
fn mesh_chunk(
    layout: &ChunkLayout<f32>,
    field: &PaintStack<'_, BoxExact<f32>, Sphere<f32>, Sphere<f32>>,
    origin: [f32; 3],
    scratch: &mut Vec<[f32; 4]>,
) -> Option<Mesh> {
    let shape = layout.sample_shape().ok()?;
    let mut builder = MeshBuilder::new();
    isomesh::marching_cubes::MarchingCubes::<f32>::new()
        .extract(field, &shape, origin, layout.cell_size(), &mut builder)
        .ok()?;
    if builder.indices().is_empty() {
        return None;
    }
    isomesh::paint::shade(builder.positions(), field, scratch);
    core::mem::swap(builder.colors_mut(), scratch);
    Some(builder.into_mesh())
}

/// Sphere-trace the current field to find what the player is pointing at.
///
/// The field is an [`Sdf`] and the wall is an exact box, so marching by the
/// field value converges in a handful of steps and lands *on* the surface —
/// which is what spraying needs. Returns `None` when the ray leaves the volume.
fn aim<F: Sdf<Scalar = f32>>(field: &F, origin: Vec3, direction: Vec3) -> Option<Vec3> {
    let mut t = 0.0f32;
    for _ in 0..96 {
        let p = origin + direction * t;
        let d = field.sample([p.x, p.y, p.z]);
        if d < 0.002 {
            return Some(p);
        }
        t += d.max(0.01);
        if t > 24.0 {
            return None;
        }
    }
    None
}

/// Every chunk whose box overlaps a world-space box, inserted into `dirty`.
///
/// A spray needs this because [`mark_edit`] cannot see it: the field is
/// unchanged, so the honest answer to "what moved" is "nothing", and what needs
/// re-shading has to be named geometrically instead.
fn dirty_box(layout: &ChunkLayout<f32>, min: Vec3, max: Vec3, dirty: &mut DirtySet) {
    let lo = layout.chunk_of([min.x, min.y, min.z]).coords;
    let hi = layout.chunk_of([max.x, max.y, max.z]).coords;
    for z in lo[2]..=hi[2] {
        for y in lo[1]..=hi[1] {
            for x in lo[0]..=hi[0] {
                dirty.insert(ChunkId::new([x, y, z]));
            }
        }
    }
}

/// Spray, carve, and re-mesh what either one touched.
#[allow(clippy::too_many_arguments)]
fn edit(
    buttons: Res<ButtonInput<MouseButton>>,
    keys: Res<ButtonInput<KeyCode>>,
    mut world: ResMut<World>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    material: Res<SurfaceMaterial>,
    camera: Query<&Transform, With<Camera3d>>,
    chunks: Query<(Entity, &Chunk)>,
    mut auto: ResMut<AutoPaint>,
) {
    for (index, key) in [
        KeyCode::Digit1,
        KeyCode::Digit2,
        KeyCode::Digit3,
        KeyCode::Digit4,
        KeyCode::Digit5,
    ]
    .into_iter()
    .enumerate()
    {
        if keys.just_pressed(key) {
            world.ink = index;
        }
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        world.radius = (world.radius - 0.06).max(0.12);
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        world.radius = (world.radius + 0.06).min(1.0);
    }

    let scripted = if auto.remaining > 0 {
        auto.remaining -= 1;
        let action = AutoPaint::action(auto.step, auto.total);
        // The script drives its own nozzle through the same fields the keys
        // set, so a scripted run and a played one take one code path. The
        // colour cycles: a run that sprays one colour cannot tell paint staying
        // put apart from paint being repainted the same shade, so the drift
        // number would read zero either way.
        world.radius = action.2;
        world.ink = (auto.step as usize) % PALETTE.len();
        auto.step += 1;
        Some(action)
    } else {
        None
    };
    let clear = keys.just_pressed(KeyCode::KeyX);
    let spray = matches!(scripted, Some((_, true, _)))
        || (world.grabbed && buttons.just_pressed(MouseButton::Left));
    let carve = matches!(scripted, Some((_, false, _)))
        || (world.grabbed && buttons.just_pressed(MouseButton::Right));
    if !(spray || carve || clear) {
        return;
    }

    let Ok(view) = camera.single() else {
        return;
    };

    let started = Instant::now();
    let layout = world.layout;
    let radius = world.radius;

    let mut dirty = DirtySet::new();

    if clear {
        world.edits.clear();
        world.probes.clear();
        world.worst_drift = 0.0;
        world.last_action = "cleared the log";
        dirty_box(
            &layout,
            Vec3::new(-3.0, -2.0, -1.0),
            Vec3::new(3.0, 2.0, 1.0),
            &mut dirty,
        );
    } else {
        // Where the edit lands: the scripted point, or whatever the player is
        // looking at.
        let target = match scripted {
            Some((point, _, _)) => Some(point),
            None => {
                let field = world.field();
                aim(&field, view.translation, *view.forward())
            }
        };
        let Some(hit) = target else {
            return;
        };
        let centre = [hit.x, hit.y, hit.z];

        if spray {
            let color = linear(PALETTE[world.ink].1);
            world.edits.push(Edit::Spray(Splat {
                shape: Sphere {
                    center: centre,
                    radius,
                },
                color,
                softness: radius * 0.4,
                depth: 0.12,
            }));
            world.last_action = "sprayed";

            // Record probes on what was just painted, so a later carve has
            // something to be measured against.
            let fresh: Vec<Probe> = {
                let field = world.field();
                (0..8)
                    .map(|i| {
                        let angle = i as f32 * core::f32::consts::TAU / 8.0;
                        let offset = Vec3::new(angle.cos(), angle.sin(), 0.0) * radius * 0.45;
                        let at = hit + offset;
                        Probe {
                            color: field.color_at([at.x, at.y, at.z]),
                            at,
                        }
                    })
                    .collect()
            };
            world.probes.extend(fresh);

            // Only colour moved, so `mark_edit` would report nothing. Name the
            // nozzle's reach directly instead.
            let reach = radius + layout.cell_size();
            dirty_box(
                &layout,
                hit - Vec3::splat(reach),
                hit + Vec3::splat(reach),
                &mut dirty,
            );
        } else {
            world.edits.push(Edit::Carve(Brush::subtract(Sphere {
                center: centre,
                radius,
            })));
            world.last_action = "carved";

            // Geometry moved, so this is `game_dig`'s route: compare the field
            // either side of the one edit just pushed. Padded by a cell for
            // M-32's reason -- `cell_of` inverts `world_of_sample` reliably in a
            // cell's interior and not on its corner.
            let split = world.edits.len() - 1;
            let before = PaintStack {
                base: wall(),
                edits: &world.edits[..split],
                background: background(),
            };
            let after = PaintStack {
                base: wall(),
                edits: &world.edits,
                background: background(),
            };
            let reach = radius + layout.cell_size();
            let min_cell = layout.cell_of([hit.x - reach, hit.y - reach, hit.z - reach]);
            let max_cell = layout.cell_of([hit.x + reach, hit.y + reach, hit.z + reach]);
            mark_edit(&layout, &before, &after, min_cell, max_cell, &mut dirty)
                .expect("a paint brush spans a few cells, far inside the u32 sample space");
        }
    }

    let field = world.field();
    let mut scratch = Vec::new();
    let mut rebuilt = 0usize;
    dirty.mesh_dirty(&layout, |id, origin| {
        rebuilt += 1;
        let mesh = mesh_chunk(&layout, &field, origin, &mut scratch).map(|m| meshes.add(m));
        let existing = chunks.iter().find(|(_, c)| c.0 == id).map(|(e, _)| e);
        match (existing, mesh) {
            (Some(entity), Some(handle)) => {
                commands.entity(entity).insert(Mesh3d(handle));
            }
            (Some(entity), None) => {
                commands.entity(entity).despawn();
            }
            (None, Some(handle)) => {
                commands.spawn((
                    Mesh3d(handle),
                    MeshMaterial3d(material.0.clone()),
                    Chunk(id),
                ));
            }
            (None, None) => {}
        }
    });

    // The measurement this example exists for, and the reason it is scoped to
    // carves: **a spray is supposed to change the colour**, so comparing across
    // one would report a large drift for paint behaving perfectly. The claim is
    // narrower and sharper than "colour never changes" — it is that *removing
    // geometry* does not disturb the colour of what remains.
    //
    // So a spray re-baselines every probe (the wall now legitimately looks
    // different) and a carve is measured against that baseline without
    // re-baselining. A run of carves therefore accumulates.
    let (drift, on_surface) = {
        let field = world.field();
        let mut worst = 0.0f32;
        let mut on_surface = 0usize;
        for probe in &world.probes {
            let at = [probe.at.x, probe.at.y, probe.at.z];
            let now = field.color_at(at);
            for (channel, was) in now.iter().zip(&probe.color) {
                worst = worst.max((channel - was).abs());
            }
            if field.sample(at).abs() < layout.cell_size() {
                on_surface += 1;
            }
        }
        (worst, on_surface)
    };
    if spray || clear {
        let baselines: Vec<[f32; 4]> = {
            let field = world.field();
            world
                .probes
                .iter()
                .map(|p| field.color_at([p.at.x, p.at.y, p.at.z]))
                .collect()
        };
        for (probe, fresh) in world.probes.iter_mut().zip(baselines) {
            probe.color = fresh;
        }
    } else {
        world.worst_drift = world.worst_drift.max(drift);
    }
    world.probes_on_surface = on_surface;
    world.last_ms = started.elapsed().as_secs_f64() * 1000.0;
    world.last_chunks = rebuilt;

    if scripted.is_some() {
        info!(
            "{:>8} {:>2}: {} chunks in {:.3} ms, drift {:.6} over {} probes ({} still on surface)",
            world.last_action,
            world.edits.len(),
            rebuilt,
            world.last_ms,
            drift,
            world.probes.len(),
            on_surface
        );
    }
}

fn report(
    world: Res<World>,
    mut stats: ResMut<DemoStats>,
    meshes: Res<Assets<Mesh>>,
    chunks: Query<(&Chunk, &Mesh3d)>,
) {
    let mut resident = 0usize;
    let mut vertices = 0usize;
    let mut triangles = 0usize;
    for (_, handle) in &chunks {
        resident += 1;
        if let Some(mesh) = meshes.get(&handle.0) {
            vertices += mesh.count_vertices();
            triangles += mesh.indices().map_or(0, |i| i.len() / 3);
        }
    }

    let sprays = world
        .edits
        .iter()
        .filter(|e| matches!(e, Edit::Spray(_)))
        .count();
    let carves = world.edits.len() - sprays;

    stats.title = format!("E-208 game paint - {resident} chunks resident");
    stats.vertices = vertices;
    stats.triangles = triangles;
    stats.extract_ms = world.last_ms;
    stats.extra = vec![
        format!(
            "log {:>3} edits = {sprays} sprays + {carves} carves    nozzle {:.2} {}",
            world.edits.len(),
            world.radius,
            PALETTE[world.ink].0
        ),
        String::new(),
        format!(
            "last: {} - {} chunks re-meshed in {:.2} ms",
            world.last_action, world.last_chunks, world.last_ms
        ),
        String::new(),
        format!(
            "PAINT DRIFT {:.6}   over {} probes, {} still on surface",
            world.worst_drift,
            world.probes.len(),
            world.probes_on_surface
        ),
        "  worst per-channel change across a carve; sprays re-baseline, as they must".to_string(),
        "  zero by construction: paint is in the field, so a carve cannot move it".to_string(),
        String::new(),
        "[LMB] spray  [RMB] blow a hole  [1-5] colour  [ ] nozzle".to_string(),
        "[WASD/QE] move  [Shift] fast  [Tab] cursor  [X] clear log".to_string(),
    ];
}

fn grab(
    keys: Res<ButtonInput<KeyCode>>,
    buttons: Res<ButtonInput<MouseButton>>,
    mut cursor: Query<&mut CursorOptions, With<PrimaryWindow>>,
    mut world: ResMut<World>,
) {
    let Ok(mut cursor) = cursor.single_mut() else {
        return;
    };
    if buttons.just_pressed(MouseButton::Left) && !world.grabbed {
        world.grabbed = true;
    }
    if keys.just_pressed(KeyCode::Tab) {
        world.grabbed = !world.grabbed;
    }
    let (mode, visible) = if world.grabbed {
        (CursorGrabMode::Locked, false)
    } else {
        (CursorGrabMode::None, true)
    };
    cursor.grab_mode = mode;
    cursor.visible = visible;
}

fn fly(
    keys: Res<ButtonInput<KeyCode>>,
    motion: Res<AccumulatedMouseMotion>,
    time: Res<Time>,
    world: Res<World>,
    mut look: ResMut<Look>,
    mut camera: Query<&mut Transform, With<Camera3d>>,
) {
    let Ok(mut transform) = camera.single_mut() else {
        return;
    };
    if world.grabbed {
        let sensitivity = 0.0022;
        look.yaw -= motion.delta.x * sensitivity;
        look.pitch = (look.pitch - motion.delta.y * sensitivity).clamp(-1.5, 1.5);
    }
    transform.rotation = Quat::from_euler(EulerRot::YXZ, look.yaw, look.pitch, 0.0);

    let mut direction = Vec3::ZERO;
    let forward = *transform.forward();
    let right = *transform.right();
    for (key, delta) in [
        (KeyCode::KeyW, forward),
        (KeyCode::KeyS, -forward),
        (KeyCode::KeyA, -right),
        (KeyCode::KeyD, right),
        (KeyCode::KeyQ, -Vec3::Y),
        (KeyCode::KeyE, Vec3::Y),
    ] {
        if keys.pressed(key) {
            direction += delta;
        }
    }
    let speed = if keys.pressed(KeyCode::ShiftLeft) {
        6.0
    } else {
        2.5
    };
    if direction != Vec3::ZERO {
        transform.translation += direction.normalize() * speed * time.delta_secs();
    }
}
