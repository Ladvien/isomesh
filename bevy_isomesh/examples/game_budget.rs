//! E-206 — the same work, spread. And what the guarantee costs.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_budget --release
//! ```
//!
//! **Always `--release`.**
//!
//! `Space` dirty the whole world · `U` drain it **unbudgeted** · `[` `]` budget ·
//! `R` reset the statistics.
//!
//! # The number a game cares about, and the one papers report
//!
//! A meshing paper reports throughput: chunks per second, or milliseconds for a
//! grid. A game cannot spend that number, because it does not get a second — it
//! gets 16.7 milliseconds and then it has to draw something. The question is not
//! *how fast does the queue drain* but **what does it cost the frame it lands
//! on**, and that figure is amortized, not total.
//!
//! So this demo overloads an edit queue on purpose — every chunk in the world
//! dirty at once — and drains it under [`DirtySet::mesh_within_budget`], with the
//! frame time on screen the whole way.
//!
//! # Flat frame time means nothing without the version that spikes
//!
//! `U` drains the identical queue through [`DirtySet::mesh_dirty`], which does
//! all of it now. That is the control, and it is the only thing that makes the
//! budgeted run's flatness worth reporting: **the total work is the same either
//! way.** A budget does not make meshing cheaper. It decides which frame pays,
//! and the headline is that the sum is unchanged while the maximum is not.
//!
//! # What the never-livelock guarantee actually costs
//!
//! `mesh_within_budget` consults its predicate **after** each chunk, never
//! before, so a budget too small for one chunk still meshes one. The doc calls
//! overshooting by at most one chunk "the price", and prices the alternative — a
//! queue that cannot make progress while it grows — as worse. Both halves of that
//! are correct and neither had a number.
//!
//! This demo measures the overshoot: how far past the budget a frame actually
//! ran, in milliseconds and as a multiple of the budget asked for. Set the budget
//! below one chunk's cost with `[` and the overshoot *is* the whole frame, which
//! is the guarantee doing exactly what it promises and is worth seeing rather
//! than trusting.

mod common;

use std::time::{Duration, Instant};

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoStats, OrbitCamera, ViewFlags};
use isomesh::Sdf;
use isomesh::chunk::dirty::DirtySet;
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::marching_cubes::MarchingCubes;

const CHUNK_CELLS: u32 = 16;
const CELL_SIZE: f32 = 0.25;
/// Chunks along x/z, and layers in y. 12 x 2 x 12 = 288 chunks.
///
/// Sized so the **unbudgeted** control misses a frame outright: at roughly
/// 0.066 ms a chunk, 288 of them is about 19 ms against a 16.7 ms frame. With
/// the 128 the first version used the control came in at 8.5 ms, which is a
/// hitch nobody would notice and would have made the comparison look like a
/// smaller win than it is.
const SPAN: i32 = 12;
const LAYERS: std::ops::RangeInclusive<i32> = 0..=1;

/// Budgets the `[` and `]` keys step through, in microseconds.
///
/// The low end is deliberately below one chunk's cost, because that is where the
/// never-livelock guarantee starts paying and the overshoot becomes visible.
const BUDGETS_US: [u64; 8] = [25, 50, 200, 500, 1_000, 2_000, 4_000, 8_000];

/// `ISOMESH_BUDGET_US` picks one without a keyboard, so the sweep below is a
/// shell loop rather than eight screenshots.
fn budget_from_env() -> usize {
    let Ok(us) = std::env::var("ISOMESH_BUDGET_US") else {
        return DEFAULT_BUDGET;
    };
    let Ok(us) = us.parse::<u64>() else {
        return DEFAULT_BUDGET;
    };
    BUDGETS_US
        .iter()
        .position(|&b| b == us)
        .unwrap_or(DEFAULT_BUDGET)
}

/// 2 ms -- an eighth of a frame, the sort of slice a game gives background work.
const DEFAULT_BUDGET: usize = 5;

/// A blobby solid, cheap enough that a frame's budget buys several chunks.
#[derive(Clone, Copy)]
struct Blobs;

impl Sdf for Blobs {
    type Scalar = f32;
    fn sample(&self, p: [f32; 3]) -> f32 {
        let ground = p[1] - 1.4 * (p[0] * 0.35).sin() * (p[2] * 0.31).cos() - 2.0;
        let r = (p[0] * p[0] + (p[1] - 2.0) * (p[1] - 2.0) + p[2] * p[2]).sqrt();
        ground.min(r - 6.0)
    }
}

#[derive(Resource)]
struct Bench {
    layout: ChunkLayout<f32>,
    dirty: DirtySet,
    budget: usize,
    /// Set while a budgeted drain is in progress.
    draining: bool,
    /// Statistics for the drain in progress, or the last one.
    frames: u32,
    chunks: usize,
    total_ms: f64,
    worst_frame_ms: f64,
    /// Worst amount by which one frame ran past the budget it was given.
    worst_overshoot_ms: f64,
    /// What the same queue cost when drained all at once.
    unbudgeted_ms: Option<f64>,
    unbudgeted_chunks: usize,
    /// Milliseconds the most recent frame spent meshing.
    last_frame_ms: f64,
    /// The unbudgeted control runs itself once, so a still is self-contained.
    control_done: bool,
}

impl Bench {
    fn budget(&self) -> Duration {
        Duration::from_micros(BUDGETS_US[self.budget])
    }
    fn budget_ms(&self) -> f64 {
        self.budget().as_secs_f64() * 1000.0
    }
    fn reset_stats(&mut self) {
        self.frames = 0;
        self.chunks = 0;
        self.total_ms = 0.0;
        self.worst_frame_ms = 0.0;
        self.worst_overshoot_ms = 0.0;
        self.last_frame_ms = 0.0;
    }
}

#[derive(Component)]
struct Chunk(ChunkId);

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-206 game budget".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, drain, report, hud).chain())
        .run();
}

fn all_chunks() -> impl Iterator<Item = ChunkId> {
    (0..SPAN)
        .flat_map(|x| LAYERS.flat_map(move |y| (0..SPAN).map(move |z| ChunkId::new([x, y, z]))))
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    let layout = ChunkLayout::new(CHUNK_CELLS, CELL_SIZE, [0.0; 3]).expect("a valid layout");
    let span = SPAN as f32 * CHUNK_CELLS as f32 * CELL_SIZE;
    for mut orbit in &mut camera {
        orbit.focus = Vec3::new(span * 0.5, 0.0, span * 0.5);
        orbit.radius = span * 1.7;
        orbit.yaw = 0.7;
        // Looking down at the whole grid: the subject is 128 chunks re-meshing,
        // and from ground level all you see is the nearest one.
        orbit.pitch = 0.75;
    }

    let material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.68, 0.66, 0.60),
        perceptual_roughness: 0.8,
        ..default()
    });
    for id in all_chunks() {
        let Some(mesh) = mesh_chunk(&layout, id) else {
            continue;
        };
        commands.spawn((
            Mesh3d(meshes.add(mesh)),
            MeshMaterial3d(material.clone()),
            Transform::default(),
            Chunk(id),
        ));
    }

    commands.insert_resource(Bench {
        layout,
        dirty: DirtySet::new(),
        budget: budget_from_env(),
        draining: false,
        frames: 0,
        chunks: 0,
        total_ms: 0.0,
        worst_frame_ms: 0.0,
        worst_overshoot_ms: 0.0,
        unbudgeted_ms: None,
        unbudgeted_chunks: 0,
        last_frame_ms: 0.0,
        control_done: false,
    });
}

/// Extract one chunk.
fn mesh_chunk(layout: &ChunkLayout<f32>, id: ChunkId) -> Option<Mesh> {
    let shape = layout.sample_shape().ok()?;
    let origin = layout.sample_origin(id);
    let mut builder = MeshBuilder::new();
    MarchingCubes::<f32>::new()
        .extract(&Blobs, &shape, origin, layout.cell_size(), &mut builder)
        .ok()?;
    Some(builder.into_mesh())
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut bench: ResMut<Bench>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut chunks: Query<(&mut Mesh3d, &Chunk)>,
    mut flags: ResMut<ViewFlags>,
) {
    flags.grid = false;

    // **The queue re-fills the moment it empties.** "Deliberately overloaded"
    // means the backlog never reaches zero, which is the state a game is in
    // while a player is carving -- and it also means any capture, at any frame,
    // is mid-drain. The first version only re-dirtied during a capture sequence,
    // so a single screenshot photographed an idle world and every number on it
    // was 0.00.
    let auto = !bench.draining;
    if auto || (!capture.is_active() && keys.just_pressed(KeyCode::Space)) {
        for id in all_chunks() {
            bench.dirty.insert(id);
        }
        // Statistics accumulate across refills rather than resetting, so the
        // amortized figure is over many passes instead of one.
        bench.draining = true;
    }
    if capture.is_active() {
        return;
    }

    if keys.just_pressed(KeyCode::BracketRight) {
        bench.budget = (bench.budget + 1).min(BUDGETS_US.len() - 1);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        bench.budget = bench.budget.saturating_sub(1);
    }
    if keys.just_pressed(KeyCode::KeyR) {
        bench.reset_stats();
        bench.unbudgeted_ms = None;
    }

    // The control runs itself once, before any budgeted drain, so a screenshot
    // carries the comparison without anyone pressing a key.
    let wanted = !bench.control_done || (!capture.is_active() && keys.just_pressed(KeyCode::KeyU));
    if wanted {
        bench.control_done = true;
        for id in all_chunks() {
            bench.dirty.insert(id);
        }
        let layout = bench.layout;
        let started = Instant::now();
        let mut built: Vec<(ChunkId, Mesh)> = Vec::new();
        let done = bench.dirty.mesh_dirty(&layout, |id, _| {
            if let Some(mesh) = mesh_chunk(&layout, id) {
                built.push((id, mesh));
            }
        });
        let elapsed = started.elapsed().as_secs_f64() * 1000.0;
        for (id, mesh) in built {
            let handle = meshes.add(mesh);
            for (mut slot, chunk) in &mut chunks {
                if chunk.0 == id {
                    slot.0 = handle.clone();
                }
            }
        }
        bench.unbudgeted_ms = Some(elapsed);
        bench.unbudgeted_chunks = done;
        bench.draining = false;
    }
}

fn drain(
    mut bench: ResMut<Bench>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut chunks: Query<(&mut Mesh3d, &Chunk)>,
    camera: Query<&OrbitCamera>,
) {
    if !bench.draining || bench.dirty.is_empty() {
        bench.draining = false;
        return;
    }
    let focus = camera.iter().next().map_or(Vec3::ZERO, |orbit| orbit.focus);

    let layout = bench.layout;
    let budget = bench.budget();
    let started = Instant::now();
    let mut built: Vec<(ChunkId, Mesh)> = Vec::new();
    let report = bench.dirty.mesh_within_budget(
        &layout,
        [focus.x, focus.y, focus.z],
        |id, _| {
            if let Some(mesh) = mesh_chunk(&layout, id) {
                built.push((id, mesh));
            }
        },
        // The predicate the crate asks for. `no_std` cannot read a clock, so the
        // caller owns it -- G-006's whole design, and the reason this is a
        // closure rather than a `Duration` argument.
        || started.elapsed() < budget,
    );
    let spent = started.elapsed().as_secs_f64() * 1000.0;

    // Uploading is the caller's cost, not the budget's, and is deliberately
    // outside the timed region: `mesh_within_budget` bounds *extraction*, and
    // charging it for asset creation would measure Bevy rather than G-006.
    for (id, mesh) in built {
        let handle = meshes.add(mesh);
        for (mut slot, chunk) in &mut chunks {
            if chunk.0 == id {
                slot.0 = handle.clone();
            }
        }
    }

    bench.frames += 1;
    bench.chunks += report.meshed;
    bench.total_ms += spent;
    bench.last_frame_ms = spent;
    if spent > bench.worst_frame_ms {
        bench.worst_frame_ms = spent;
    }
    let over = spent - bench.budget_ms();
    if over > bench.worst_overshoot_ms {
        bench.worst_overshoot_ms = over;
    }
    if report.is_drained() {
        bench.draining = false;
    }
}

/// One CSV row per drained pass, so the budget can be swept from a shell.
fn report(bench: Res<Bench>, mut last: Local<u32>) {
    if bench.frames == *last || !bench.frames.is_multiple_of(40) || bench.frames == 0 {
        return;
    }
    *last = bench.frames;
    info!(
        "budget,{:.3},{},{},{:.3},{:.3},{:.3},{:.2}",
        bench.budget_ms(),
        bench.frames,
        bench.chunks,
        bench.total_ms / f64::from(bench.frames),
        bench.worst_frame_ms,
        bench.worst_overshoot_ms,
        bench.chunks as f64 / f64::from(bench.frames),
    );
}

fn hud(bench: Res<Bench>, mut stats: ResMut<DemoStats>) {
    let remaining = bench.dirty.len();
    let mean = if bench.frames > 0 {
        bench.total_ms / f64::from(bench.frames)
    } else {
        0.0
    };
    let per_frame = if bench.frames > 0 {
        bench.chunks as f64 / f64::from(bench.frames)
    } else {
        0.0
    };
    let budget_ms = bench.budget_ms();

    stats.title = format!(
        "E-206  frame budget   {:.1} ms/frame   {}",
        budget_ms,
        if bench.draining {
            format!("draining, {remaining} left")
        } else {
            "idle -- Space to overload".to_string()
        }
    );
    stats.vertices = bench.chunks;
    stats.triangles = remaining;

    let mut lines = vec![
        format!(
            "{:<26} {:>9.2}   ms asked for per frame  [ and ]",
            "budget", budget_ms
        ),
        format!("{:<26} {:>9}   chunks still queued", "backlog", remaining),
        String::new(),
        format!(
            "{:<26} {:>9}   frames spent draining",
            "frames", bench.frames
        ),
        format!("{:<26} {:>9}   chunks meshed", "chunks", bench.chunks),
        format!("{:<26} {:>9.2}   chunks per frame", "rate", per_frame),
        format!(
            "{:<26} {:>9.2}   ms per frame, AMORTIZED",
            "mean cost", mean
        ),
        format!(
            "{:<26} {:>9.2}   ms, worst single frame",
            "peak cost", bench.worst_frame_ms
        ),
        String::new(),
        format!(
            "{:<26} {:>9.2}   ms past the budget, worst frame",
            "overshoot", bench.worst_overshoot_ms
        ),
        "                                      the predicate is consulted AFTER each".into(),
        "                                      chunk, so a budget too small for one".into(),
        "                                      still meshes one. that is the price of".into(),
        "                                      never livelocking, and this is it.".into(),
        String::new(),
    ];

    match bench.unbudgeted_ms {
        Some(ms) => {
            lines.push(format!(
                "{:<26} {:>9.2}   ms in ONE frame, {} chunks  [U]",
                "unbudgeted", ms, bench.unbudgeted_chunks
            ));
            if bench.worst_frame_ms > 0.0 {
                lines.push(format!(
                    "{:<26} {:>9.1}x  lower peak, for the same total work",
                    "budgeted peak is",
                    ms / bench.worst_frame_ms
                ));
            }
        }
        None => lines.push(format!(
            "{:<26} {:>9}   press U to drain the same queue all at once",
            "unbudgeted", "--"
        )),
    }

    lines.extend([
        String::new(),
        "a budget does not make meshing cheaper. the total is the same either".into(),
        "way; what changes is which frame pays. that is the only reason the".into(),
        "flat line above is worth anything, and why U is on the keyboard.".into(),
    ]);
    stats.extra = lines;
}
