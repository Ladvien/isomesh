//! E-202 — carving tunnels, the way a game does it.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_dig --release
//! ```
//!
//! `WASD` to move, mouse to look, **left click to carve**, right click to fill.
//! `Tab` releases the cursor. `[` and `]` change the brush radius, `X` clears
//! the edit log, `C` outlines the chunks that were re-meshed by the last edit.
//!
//! # What this is actually testing
//!
//! Not "can it mesh a field" — three examples already do that. This is the first
//! one where the mesh is **re-built while someone is holding the mouse down**,
//! and it exists to put two numbers on screen that nothing else in the repo
//! measures under load:
//!
//! - **Chunks touched per edit.** A brush changes the field everywhere, because
//!   an SDF is global; what it changes *visibly* is a shell. G-002's `mark_edit`
//!   compares the field either side of one edit and marks only the chunks whose
//!   cells actually moved.
//! - **E1 — the fraction of the brush's own bounding box that changed.** M-33
//!   measured 15–36% offline. This shows it live, per edit, and it is the number
//!   the entire incremental story rests on: if it were 100%, re-meshing the
//!   bounding box would be as cheap as being clever about it.
//!
//! # The spacing is a power of two, deliberately
//!
//! `h = 0.125`. **M-32** measured that two chunks agree on their shared sample
//! plane bit-for-bit only at a power-of-two cell size; anywhere else they differ
//! by an ulp and the seam needs A-013's weld to close. This example does not
//! weld — each chunk is its own `Mesh3d`, exactly as an engine would keep them —
//! so it uses the spacing where the seam is exact and the surface is continuous
//! without one. `chunk_seam_weld` is the example that shows the other case.
//!
//! # The edit log grows, and the cost grows with it — sub-linearly
//!
//! Edits compose rather than mutate: the field is a `BrushStack` over the base
//! terrain, and carving pushes a brush. That is what makes undo a re-fold of the
//! log rather than a snapshot (E-207's premise), and it means **every field
//! sample walks every brush**.
//!
//! So the cost grows, and it is worth being precise about how much rather than
//! waving at it. Measured over a 60-carve scripted run (`ISOMESH_AUTOCARVE=60`,
//! which prints one line per edit), median milliseconds per re-meshed chunk:
//!
//! | edits in the log | 1–15 | 16–30 | 31–45 | 46–60 |
//! |---|---|---|---|---|
//! | ms per chunk | 0.158 | 0.354 | 0.525 | 0.589 |
//!
//! **3.7× for 7× the log, and flattening** — not proportional, even though every
//! sample really does walk every brush. So the stack walk is a real cost and not
//! the dominant one at these lengths; what else is in there has not been
//! measured and is not asserted here. Press `X` to clear the log and watch it
//! drop back.

mod common;

use bevy::input::mouse::AccumulatedMouseMotion;
use bevy::platform::time::Instant;
use bevy::prelude::*;
use bevy::window::{CursorGrabMode, CursorOptions, PrimaryWindow};
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoStats, OrbitCamera};
use isomesh::Sdf;
use isomesh::brush::{Brush, BrushStack};
use isomesh::chunk::dirty::{DirtySet, EditReport, mark_edit};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::Sphere;

/// Chunk edge, in cells.
const CHUNK_CELLS: u32 = 16;
/// See the module docs: a power of two, so the seam is bit-exact without a weld.
const CELL_SIZE: f32 = 0.125;
/// Chunks along x and z, and up in y. Small enough to mesh at startup in one go.
const EXTENT: [i32; 3] = [3, 2, 3];

/// Gizmos for the re-meshed-chunk outline.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct ChunkGizmos;

/// The terrain before any edit: a slab with a rolling top.
///
/// Hand-rolled rather than `FbmTerrain`, because this needs a floor a player can
/// stand on and a ceiling to dig into, and it must be cheap — it is sampled
/// inside the edit loop.
#[derive(Clone, Copy)]
struct Ground;

impl Sdf for Ground {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        // Distance to a wavy height field, negative below it.
        let height = 0.35 * (p[0] * 0.9).sin() * (p[2] * 0.7).cos() + 0.15 * (p[0] * 2.1).sin();
        p[1] - height
    }
}

#[derive(Resource)]
struct World {
    layout: ChunkLayout<f32>,
    brushes: Vec<Brush<Sphere<f32>>>,
    dirty: DirtySet,
    radius: f32,
    /// Chunks re-meshed by the most recent edit, for the outline.
    last_touched: Vec<ChunkId>,
    last_edit: Option<EditReport>,
    last_edit_ms: f64,
    last_chunks: usize,
    show_chunks: bool,
    grabbed: bool,
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

/// A scripted sequence of carves, one per frame, driven by `ISOMESH_AUTOCARVE`.
///
/// The acceptance criterion for this ticket is about what happens *while
/// someone is clicking*, and a screenshot cannot click. Without this the example
/// could be committed compiling, rendering, and silently not carving at all —
/// so the loop runs itself, through exactly the same code path a click takes,
/// and the committed screenshot is of a tunnel that was actually dug.
#[derive(Resource, Default)]
struct AutoCarve {
    remaining: u32,
    step: u32,
    every: u32,
}

impl AutoCarve {
    fn from_env() -> Self {
        Self {
            remaining: std::env::var("ISOMESH_AUTOCARVE")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(0),
            // Captured frames per carve. One carve per frame digs the whole
            // visible tunnel in the first half-second of a clip and leaves the
            // rest of it static, which reads as a jump cut rather than as
            // digging.
            every: std::env::var("ISOMESH_AUTOCARVE_EVERY")
                .ok()
                .and_then(|v| v.parse().ok())
                .filter(|n| *n >= 1)
                .unwrap_or(1),
            step: 0,
        }
    }

    /// Where the `n`th scripted carve goes: a tunnel boring into the hill.
    fn centre(n: u32) -> Vec3 {
        let t = n as f32;
        Vec3::new(-0.9 + t * 0.30, 0.55 - t * 0.045, 2.2 - t * 0.34)
    }
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-202 game dig".into(),
                // Web only, inert on native: bind to the 1280x720 canvas the
                // page supplies rather than letting Bevy append its own. The HUD
                // panels are laid out in pixels for that size, so the canvas is
                // fixed and CSS scales it -- `fit_canvas_to_parent` stays at its
                // `false` default for the same reason.
                canvas: Some("#isomesh-canvas".into()),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<ChunkGizmos>()
        .insert_resource(Look {
            yaw: 0.0,
            pitch: -0.15,
        })
        .insert_resource(AutoCarve::from_env())
        .add_systems(Startup, setup)
        .add_systems(Update, (grab, fly, dig, report, outline_chunks))
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut config: ResMut<GizmoConfigStore>,
    camera: Query<Entity, With<OrbitCamera>>,
) {
    // The shared harness spawns an orbit camera. Take its `OrbitCamera` off
    // rather than despawning the entity: the orbit system then skips it, this
    // example drives the same camera directly, and everything else in the
    // harness that expects a camera to exist still finds one.
    for entity in &camera {
        commands
            .entity(entity)
            .remove::<OrbitCamera>()
            .insert(Transform::from_xyz(0.0, 1.6, 6.0));
    }
    let (chunk_gizmos, _) = config.config_mut::<ChunkGizmos>();
    chunk_gizmos.line.width = 2.0;

    let layout = ChunkLayout::<f32>::new(
        CHUNK_CELLS,
        CELL_SIZE,
        [
            -(EXTENT[0] as f32) * CHUNK_CELLS as f32 * CELL_SIZE * 0.5,
            -1.4,
            -(EXTENT[2] as f32) * CHUNK_CELLS as f32 * CELL_SIZE * 0.5,
        ],
    )
    .expect("valid layout");

    let material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.62, 0.58, 0.52),
        perceptual_roughness: 0.85,
        ..default()
    });

    let mut world = World {
        layout,
        brushes: Vec::new(),
        dirty: DirtySet::new(),
        radius: 0.55,
        last_touched: Vec::new(),
        last_edit: None,
        last_edit_ms: 0.0,
        last_chunks: 0,
        show_chunks: true,
        grabbed: false,
    };

    // Mesh every chunk once. After this, only edited chunks are ever re-meshed.
    for z in 0..EXTENT[2] {
        for y in 0..EXTENT[1] {
            for x in 0..EXTENT[0] {
                world.dirty.insert(ChunkId::new([x, y, z]));
            }
        }
    }
    let field = BrushStack {
        base: Ground,
        brushes: &world.brushes,
    };
    let layout = world.layout;
    world.dirty.mesh_dirty(&layout, |id, origin| {
        if let Some(mesh) = mesh_chunk(&layout, &field, origin) {
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

/// Extract one chunk. Returns `None` when the chunk holds no surface, so empty
/// air costs an entity with an empty mesh rather than a draw call over nothing.
fn mesh_chunk<F: Sdf<Scalar = f32>>(
    layout: &ChunkLayout<f32>,
    field: &F,
    origin: [f32; 3],
) -> Option<Mesh> {
    let shape = layout.sample_shape().ok()?;
    let mut builder = MeshBuilder::new();
    isomesh::marching_cubes::MarchingCubes::<f32>::new()
        .extract(field, &shape, origin, layout.cell_size(), &mut builder)
        .ok()?;
    if builder.indices().is_empty() {
        return None;
    }
    Some(builder.into_mesh())
}

fn grab(
    keys: Res<ButtonInput<KeyCode>>,
    buttons: Res<ButtonInput<MouseButton>>,
    // `CursorOptions` is its own component on the window entity in Bevy 0.19,
    // not a field of `Window`.
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

/// The loop this example exists for: one click, one brush, one incremental
/// re-mesh.
#[allow(clippy::too_many_arguments)]
fn dig(
    buttons: Res<ButtonInput<MouseButton>>,
    keys: Res<ButtonInput<KeyCode>>,
    mut world: ResMut<World>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    material: Res<SurfaceMaterial>,
    camera: Query<&Transform, With<Camera3d>>,
    chunks: Query<(Entity, &Chunk)>,
    mut auto: ResMut<AutoCarve>,
    capture: Res<Capture>,
) {
    if keys.just_pressed(KeyCode::BracketLeft) {
        world.radius = (world.radius - 0.1).max(0.2);
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        world.radius = (world.radius + 0.1).min(1.6);
    }
    if keys.just_pressed(KeyCode::KeyC) {
        world.show_chunks = !world.show_chunks;
    }

    // **One carve per *captured* frame, not per rendered frame.** The tunnel
    // bores away from the camera at 0.34 units a step, so at 60 Hz it is out of
    // frame in a quarter of a second -- and a recording cannot begin until the
    // window has stopped resizing, which is 30 frames in. A clip made that way
    // photographs the aftermath while all the action happens off-screen. Every
    // other capture-driven example in this directory paces itself off
    // `capture.taken` for the same reason.
    let ready = !capture.is_active() || capture.taken > 0;
    let due = !capture.is_active() || auto.step <= capture.taken / auto.every.max(1);
    let scripted = if auto.remaining > 0 && ready && due {
        auto.remaining -= 1;
        let centre = AutoCarve::centre(auto.step);
        auto.step += 1;
        Some(centre)
    } else {
        None
    };
    let clear = keys.just_pressed(KeyCode::KeyX);
    let carve = scripted.is_some() || (world.grabbed && buttons.just_pressed(MouseButton::Left));
    let fill = world.grabbed && buttons.just_pressed(MouseButton::Right);
    if !(carve || fill || clear) {
        return;
    }

    let Ok(view) = camera.single() else {
        return;
    };

    let started = Instant::now();
    let layout = world.layout;

    // The region to re-check. For an edit it is the brush's own bounding box,
    // padded by a cell so a crossing exactly on the boundary is included; for a
    // clear it is every chunk, because undoing the whole log can change anything
    // any brush ever touched.
    let (min_cell, max_cell) = if clear {
        world.brushes.clear();
        let cells = i64::from(CHUNK_CELLS);
        (
            [0, 0, 0],
            [
                i64::from(EXTENT[0]) * cells,
                i64::from(EXTENT[1]) * cells,
                i64::from(EXTENT[2]) * cells,
            ],
        )
    } else {
        let centre =
            scripted.unwrap_or_else(|| view.translation + *view.forward() * (world.radius + 1.2));
        let shape = Sphere {
            center: [centre.x, centre.y, centre.z],
            radius: world.radius,
        };
        let brush = if carve {
            Brush::subtract(shape)
        } else {
            Brush::add(shape)
        };
        world.brushes.push(brush);

        // Padded by a cell, which is not tidiness: `cell_of` inverts
        // `world_of_sample` in a cell's interior and not reliably on its corner,
        // for the same power-of-two reason as M-32. A padded range cannot lose a
        // crossing to that; an exact one can.
        let reach = world.radius + layout.cell_size();
        (
            layout.cell_of([centre.x - reach, centre.y - reach, centre.z - reach]),
            layout.cell_of([centre.x + reach, centre.y + reach, centre.z + reach]),
        )
    };

    // `before` is the field without the brush just pushed, `after` is with it.
    // Splitting the log rather than keeping two copies is what makes this
    // exact: the two fields differ by precisely one term.
    let split = if clear { 0 } else { world.brushes.len() - 1 };
    let before = BrushStack {
        base: Ground,
        brushes: &world.brushes[..split],
    };
    let after = BrushStack {
        base: Ground,
        brushes: &world.brushes,
    };
    let mut dirty = DirtySet::new();
    let report = mark_edit(&layout, &before, &after, min_cell, max_cell, &mut dirty)
        .expect("a dig brush spans a few cells, far inside the u32 sample space");

    let touched: Vec<ChunkId> = dirty.iter().collect();
    let field = BrushStack {
        base: Ground,
        brushes: &world.brushes,
    };
    let mut rebuilt = 0usize;
    dirty.mesh_dirty(&layout, |id, origin| {
        rebuilt += 1;
        let mesh = mesh_chunk(&layout, &field, origin).map(|m| meshes.add(m));
        let existing = chunks.iter().find(|(_, c)| c.0 == id).map(|(e, _)| e);
        match (existing, mesh) {
            (Some(entity), Some(handle)) => {
                commands.entity(entity).insert(Mesh3d(handle));
            }
            (Some(entity), None) => {
                // Carved away entirely. Despawning beats keeping an empty mesh.
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

    world.last_edit_ms = started.elapsed().as_secs_f64() * 1000.0;
    // One line per scripted edit, so the log-growth claim in the module docs can
    // be checked from a terminal rather than by reading a HUD off a screenshot.
    if scripted.is_some() {
        info!(
            "edit {:>3}: {} chunks in {:.3} ms, E1 {:.1}%",
            world.brushes.len(),
            rebuilt,
            world.last_edit_ms,
            100.0 * report.changed_fraction()
        );
    }
    world.last_edit = Some(report);
    world.last_chunks = rebuilt;
    world.last_touched = touched;
}

fn report(
    world: Res<World>,
    mut stats: ResMut<DemoStats>,
    meshes: Res<Assets<Mesh>>,
    chunks: Query<(&Chunk, &Mesh3d)>,
) {
    // Totals across every resident chunk, read back from the assets rather than
    // tracked alongside them -- a running counter would be one more thing that
    // can disagree with what is actually on screen.
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
    stats.title = format!("E-202 game dig - {resident} chunks resident");
    stats.vertices = vertices;
    stats.triangles = triangles;
    stats.extract_ms = world.last_edit_ms;
    stats.extra = vec![
        format!(
            "edit log {:>4} brushes    brush radius {:.2}",
            world.brushes.len(),
            world.radius
        ),
        String::new(),
        format!(
            "last edit: {:>3} chunks re-meshed in {:.2} ms",
            world.last_chunks, world.last_edit_ms
        ),
        match world.last_edit {
            // `output_changed_cells`, not `value_changed_cells`. M-34: counting
            // cells whose *samples* moved reads 100% and says incremental
            // meshing is pointless; counting cells whose *triangles* move is
            // 15-36% and says the opposite. E1 is the second one.
            Some(r) => format!(
                "           {} of {} cells in the box re-mesh = E1 {:.1}%  ({} moved a sample)",
                r.output_changed_cells,
                r.region_cells,
                100.0 * r.changed_fraction(),
                r.value_changed_cells
            ),
            None => "           (click to carve)".to_string(),
        },
        match world.last_edit {
            Some(r) => format!(
                "           {} of {} chunks in the box were dirty = {:.1}%",
                r.dirty_chunks,
                r.region_chunks,
                100.0 * r.dirty_chunk_fraction()
            ),
            None => String::new(),
        },
        String::new(),
        "every field sample walks the log: measured 3.7x ms/chunk for 7x the log".to_string(),
        String::new(),
        "[LMB] carve  [RMB] fill  [WASD/QE] move  [Shift] fast".to_string(),
        "[Tab] cursor  [ ] radius  [X] clear log  [C] chunk outlines".to_string(),
    ];
}

/// Outline the chunks the last edit re-meshed.
///
/// This is the "chunks-touched count on screen" the ticket asks for, made
/// spatial: the count says how much work an edit cost, the boxes say *where*,
/// and only the second one shows you that a brush straddling a corner costs
/// eight chunks rather than one.
fn outline_chunks(world: Res<World>, mut gizmos: Gizmos<ChunkGizmos>) {
    if !world.show_chunks {
        return;
    }
    let span = world.layout.cell_size() * CHUNK_CELLS as f32;
    for id in &world.last_touched {
        let origin = world.layout.sample_origin(*id);
        let centre = Vec3::new(
            origin[0] + span * 0.5,
            origin[1] + span * 0.5,
            origin[2] + span * 0.5,
        );
        gizmos.cube(
            Transform::from_translation(centre).with_scale(Vec3::splat(span)),
            Color::srgb(0.20, 0.85, 1.0),
        );
    }
}
