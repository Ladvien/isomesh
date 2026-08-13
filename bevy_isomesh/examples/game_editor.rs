//! E-207 — undo is a re-fold, and the log's order is load-bearing.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_editor --release
//! ```
//!
//! **Always `--release`.**
//!
//! `Z` undo · `Y` redo · `E` one more edit · `S` swap the last two ops ·
//! `X` clear.
//!
//! # Undo without a snapshot
//!
//! The edits are a **log**, and the field is a fold of that log over a base:
//! `BrushStack { base, brushes: &log[..cursor] }`. Undo moves the cursor back by
//! one. Nothing is stored, nothing is copied, and there is no "before" image
//! anywhere — the previous state is *recomputed* from the ops that made it.
//!
//! That is the property a CAD tool wants, because a snapshot of a voxel world is
//! the whole world and a log entry is a shape and an enum. It costs the re-fold,
//! which is measured here: every field sample walks the whole log, so the cost of
//! an undo grows with how much history is in front of it.
//!
//! # The order of the log is not free to change, and the crate says so
//!
//! `BrushOp::commutes_with` returns *"the honest answer rather than an optimistic
//! one: only identical hard operations commute."* M-36 measured a run of
//! same-kind hard edits reordering **bit-for-bit** — one result from all 40,320
//! orderings — and M-37 measured an add/subtract boundary giving **11 distinct
//! results**, a difference that is semantic and that no storage format repairs.
//!
//! `S` swaps the last two entries in the log, re-folds, and compares the mesh
//! hash before and after. So the predicate is not taken on trust: when it says
//! two ops commute the hash must be **identical**, and when it says they do not
//! the hash must **differ**. A predicate that were merely conservative — refusing
//! to promise commutation that actually holds — would show up as "said no, and
//! nothing changed", and the HUD names that case explicitly.
//!
//! # What is measured
//!
//! - **Round trip.** Undo then redo must return a bit-identical mesh. The hash is
//!   over every position and index in the world, so a single vertex moving by one
//!   bit anywhere fails it.
//! - **Re-fold cost against log length**, which is the price of not snapshotting.
//! - **Whether `commutes_with` predicts reality**, per swap.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::to_bevy_mesh;
use common::{Capture, CommonPlugin, DemoStats, OrbitCamera, ViewFlags};
use isomesh::brush::{Brush, BrushOp, BrushStack};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::Sphere;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, Sdf};

const CHUNK_CELLS: u32 = 16;
const CELL_SIZE: f32 = 0.125;
const SPAN: i32 = 4;
const LAYERS: std::ops::RangeInclusive<i32> = 0..=1;

/// The block being carved: a slab with a gently rolling top.
#[derive(Clone, Copy)]
struct Stock;

impl Sdf for Stock {
    type Scalar = f32;
    fn sample(&self, p: [f32; 3]) -> f32 {
        let top = 2.4 + 0.35 * (p[0] * 1.1).sin() * (p[2] * 0.9).cos();
        let d = [p[0] - 4.0, p[1] - top * 0.5, p[2] - 4.0];
        let half = [3.4, top * 0.5, 3.4];
        let q = [
            d[0].abs() - half[0],
            d[1].abs() - half[1],
            d[2].abs() - half[2],
        ];
        let outside = [q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)];
        (outside[0] * outside[0] + outside[1] * outside[1] + outside[2] * outside[2]).sqrt()
            + q[0].max(q[1]).max(q[2]).min(0.0)
    }
}

/// The scripted edit sequence.
///
/// **Deliberately mixed.** A run of same-kind edits would commute at every swap
/// (M-36) and the `S` key would never show a difference; alternating them puts an
/// add/subtract boundary at most positions, which is the case M-37 measured as
/// *not* commuting. A sequence that could only ever agree would demonstrate
/// nothing.
fn scripted(n: usize) -> Brush<Sphere<f32>> {
    let t = n as f32;
    // Consecutive brushes have to **overlap**, or an add/subtract boundary
    // commutes anyway and the interesting case never fires. The first version
    // moved the centre by ~3 units between ops against radii of ~1, so every
    // pair was disjoint: the audit reported 0 of 7 pairs where order mattered,
    // on a fixture built to show that order matters. Slow frequencies keep
    // consecutive centres about 1.5 apart, inside the sum of the radii.
    let shape = Sphere {
        center: [
            4.0 + 2.2 * (t * 0.42).sin(),
            1.7 + 0.7 * (t * 0.31).cos(),
            4.0 + 2.2 * (t * 0.37).cos(),
        ],
        radius: 1.15 + 0.25 * (t * 0.55).sin().abs(),
    };
    // Runs of three, not strict alternation. Alternating puts an add/subtract
    // boundary at *every* adjacent pair, so `commutes_with` would answer `false`
    // seven times out of seven and its `true` branch would never be exercised.
    // Runs give both kinds of pair.
    if n % 5 < 3 {
        Brush::subtract(shape)
    } else {
        Brush::add(shape)
    }
}

/// What one swap proved.
#[derive(Clone, Copy)]
struct Swap {
    predicted_commute: bool,
    hash_changed: bool,
}

#[derive(Resource)]
struct Editor {
    layout: ChunkLayout<f32>,
    log: Vec<Brush<Sphere<f32>>>,
    /// How much of the log is applied. Undo moves this back; the tail is kept so
    /// redo can move it forward again.
    cursor: usize,
    hash: u64,
    /// Hash of the state before the last undo, to check the round trip.
    round_trip: Option<bool>,
    refold_ms: f64,
    triangles: usize,
    last_swap: Option<Swap>,
    swaps_checked: u32,
    swaps_disagreeing: u32,
    /// Adjacent pairs whose swap left the mesh bit-identical, split by what
    /// `commutes_with` predicted. Filled by the startup audit.
    agreed_commute: u32,
    agreed_differ: u32,
    conservative: u32,
    /// Of the conservative pairs, how many had **disjoint** shapes -- which is
    /// the mechanism rather than a defect: the predicate answers on operations
    /// alone and cannot see that two brushes never touch.
    conservative_disjoint: u32,
    audited: bool,
}

#[derive(Component)]
struct Chunk(ChunkId);

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-207 game editor".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, hud).chain())
        .run();
}

fn all_chunks() -> impl Iterator<Item = ChunkId> {
    (0..SPAN)
        .flat_map(|x| LAYERS.flat_map(move |y| (0..SPAN).map(move |z| ChunkId::new([x, y, z]))))
}

/// Re-mesh the whole world from `log[..cursor]`, and hash it.
///
/// Everything is re-meshed rather than only the dirty region, because the thing
/// being measured is the **fold**, and a hash over part of the world could not
/// establish that undo and redo agree everywhere. `game_dig` already measures the
/// dirty-set path; this one measures history.
fn refold(
    layout: &ChunkLayout<f32>,
    log: &[Brush<Sphere<f32>>],
    cursor: usize,
    meshes: &mut Assets<Mesh>,
    chunks: &mut Query<(&mut Mesh3d, &Chunk)>,
) -> (u64, usize, f64) {
    let field = BrushStack {
        base: Stock,
        brushes: &log[..cursor],
    };
    let started = Instant::now();
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    let mut triangles = 0usize;
    let mut built: Vec<(ChunkId, Mesh)> = Vec::new();

    for id in all_chunks() {
        let Ok(shape) = layout.sample_shape() else {
            continue;
        };
        let mut buffer = MeshBuffer::<f32>::new();
        if MarchingCubes::<f32>::new()
            .extract(
                &field,
                &shape,
                layout.sample_origin(id),
                layout.cell_size(),
                &mut buffer,
            )
            .is_err()
        {
            continue;
        }
        // FNV-1a over the bits of every position and every index. A vertex
        // moving by one bit anywhere in the world changes this, which is what
        // makes "bit-identical" a claim rather than a hope.
        for p in &buffer.positions {
            for v in p {
                hash ^= u64::from(v.to_bits());
                hash = hash.wrapping_mul(0x100_0000_01b3);
            }
        }
        for i in &buffer.indices {
            hash ^= u64::from(*i);
            hash = hash.wrapping_mul(0x100_0000_01b3);
        }
        triangles += buffer.indices.len() / 3;

        built.push((id, to_bevy_mesh(&buffer)));
    }
    let refold_ms = started.elapsed().as_secs_f64() * 1000.0;

    for (id, mesh) in built {
        let handle = meshes.add(mesh);
        for (mut slot, chunk) in chunks.iter_mut() {
            if chunk.0 == id {
                slot.0 = handle.clone();
            }
        }
    }
    (hash, triangles, refold_ms)
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    let layout = ChunkLayout::new(CHUNK_CELLS, CELL_SIZE, [0.0; 3]).expect("a valid layout");
    for mut orbit in &mut camera {
        orbit.focus = Vec3::new(4.0, 1.2, 4.0);
        orbit.radius = 11.0;
        orbit.yaw = 0.8;
        orbit.pitch = 0.42;
    }
    let material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.72, 0.68, 0.62),
        perceptual_roughness: 0.75,
        double_sided: true,
        cull_mode: None,
        ..default()
    });
    for id in all_chunks() {
        commands.spawn((
            Mesh3d(meshes.add(Mesh::from(Cuboid::new(0.0, 0.0, 0.0)))),
            MeshMaterial3d(material.clone()),
            Transform::default(),
            Chunk(id),
        ));
    }
    commands.insert_resource(Editor {
        layout,
        log: (0..8).map(scripted).collect(),
        cursor: 8,
        hash: 0,
        round_trip: None,
        refold_ms: 0.0,
        triangles: 0,
        last_swap: None,
        swaps_checked: 0,
        swaps_disagreeing: 0,
        agreed_commute: 0,
        agreed_differ: 0,
        conservative: 0,
        conservative_disjoint: 0,
        audited: false,
    });
}

#[allow(clippy::too_many_arguments)]
fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut editor: ResMut<Editor>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut chunks: Query<(&mut Mesh3d, &Chunk)>,
    mut flags: ResMut<ViewFlags>,
    mut started: Local<bool>,
) {
    flags.grid = false;

    let act = |editor: &mut Editor, meshes: &mut Assets<Mesh>, chunks: &mut Query<_, _>| {
        let (hash, triangles, ms) =
            refold(&editor.layout, &editor.log, editor.cursor, meshes, chunks);
        editor.hash = hash;
        editor.triangles = triangles;
        editor.refold_ms = ms;
    };

    if !*started {
        *started = true;
        act(&mut editor, &mut meshes, &mut chunks);
        return;
    }

    // The audit, once, so a still carries both claims without a keypress.
    //
    // Every adjacent pair is swapped, re-folded and restored, and the mesh hash
    // is compared against what `commutes_with` promised. One keypress would
    // check one pair; this checks all of them, which is the difference between
    // an anecdote and a measurement.
    if !editor.audited {
        editor.audited = true;
        let base = editor.hash;

        // Round trip first: undo, redo, and require the world back bit-for-bit
        // *and* require the undone state to have actually differed, or the
        // check passes on a no-op.
        if editor.cursor > 0 {
            editor.cursor -= 1;
            act(&mut editor, &mut meshes, &mut chunks);
            let undone = editor.hash;
            editor.cursor += 1;
            act(&mut editor, &mut meshes, &mut chunks);
            editor.round_trip = Some(editor.hash == base && undone != base);
        }

        for i in 0..editor.cursor.saturating_sub(1) {
            let (a, b) = (editor.log[i], editor.log[i + 1]);
            let predicted = a.op.commutes_with(b.op);
            editor.log.swap(i, i + 1);
            act(&mut editor, &mut meshes, &mut chunks);
            let changed = editor.hash != base;
            editor.log.swap(i, i + 1);

            editor.swaps_checked += 1;
            match (predicted, changed) {
                (true, false) => editor.agreed_commute += 1,
                (false, true) => editor.agreed_differ += 1,
                (false, false) => {
                    editor.conservative += 1;
                    // Two spheres that never touch commute whatever their ops
                    // are. `commutes_with` takes only the ops, so it cannot
                    // know that, and saying so is sound rather than wrong.
                    let d = [
                        a.shape.center[0] - b.shape.center[0],
                        a.shape.center[1] - b.shape.center[1],
                        a.shape.center[2] - b.shape.center[2],
                    ];
                    let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                    if dist >= a.shape.radius + b.shape.radius {
                        editor.conservative_disjoint += 1;
                    }
                }
                (true, true) => editor.swaps_disagreeing += 1,
            }
        }
        act(&mut editor, &mut meshes, &mut chunks);
        info!(
            "editor,{},{},{},{},{},{},{:.2}",
            editor.swaps_checked,
            editor.agreed_commute,
            editor.agreed_differ,
            editor.conservative,
            editor.conservative_disjoint,
            editor.swaps_disagreeing,
            editor.refold_ms,
        );
        return;
    }

    // A capture cycles undo and redo so a recorded sequence shows history
    // moving, and so the round-trip check is exercised without a keyboard.
    if capture.is_active() {
        if capture.taken % 16 < 8 {
            if editor.cursor > 0 {
                editor.cursor -= 1;
                act(&mut editor, &mut meshes, &mut chunks);
            }
        } else if editor.cursor < editor.log.len() {
            editor.cursor += 1;
            act(&mut editor, &mut meshes, &mut chunks);
        }
        return;
    }

    if keys.just_pressed(KeyCode::KeyZ) && editor.cursor > 0 {
        // The round trip: remember this state, undo, redo, and require the hash
        // to come back. Checked here rather than asserted in a test because the
        // claim is about the *editor*, and an editor whose undo does not round
        // trip is broken in a way no unit test of the fold would catch.
        let before = editor.hash;
        editor.cursor -= 1;
        act(&mut editor, &mut meshes, &mut chunks);
        let undone = editor.hash;
        editor.cursor += 1;
        act(&mut editor, &mut meshes, &mut chunks);
        editor.round_trip = Some(editor.hash == before && undone != before);
        // Leave the editor in the undone state, which is what Z is for.
        editor.cursor -= 1;
        act(&mut editor, &mut meshes, &mut chunks);
    }
    if keys.just_pressed(KeyCode::KeyY) && editor.cursor < editor.log.len() {
        editor.cursor += 1;
        act(&mut editor, &mut meshes, &mut chunks);
    }
    if keys.just_pressed(KeyCode::KeyE) {
        // A new edit truncates the redo tail, as every editor does.
        let cursor = editor.cursor;
        editor.log.truncate(cursor);
        let next = scripted(editor.log.len());
        editor.log.push(next);
        editor.cursor += 1;
        act(&mut editor, &mut meshes, &mut chunks);
    }
    if keys.just_pressed(KeyCode::KeyX) {
        editor.log.clear();
        editor.cursor = 0;
        editor.round_trip = None;
        editor.last_swap = None;
        act(&mut editor, &mut meshes, &mut chunks);
    }

    // Swap the last two applied ops and find out whether the crate's own
    // prediction holds.
    if keys.just_pressed(KeyCode::KeyS) && editor.cursor >= 2 {
        let (a, b) = (editor.log[editor.cursor - 2], editor.log[editor.cursor - 1]);
        let predicted = a.op.commutes_with(b.op);
        let before = editor.hash;
        let cursor = editor.cursor;
        editor.log.swap(cursor - 2, cursor - 1);
        act(&mut editor, &mut meshes, &mut chunks);
        let changed = editor.hash != before;

        editor.swaps_checked += 1;
        // Disagreement means the predicate promised commutation and the geometry
        // moved -- an unsound `true`, which is the only genuinely wrong answer.
        if predicted && changed {
            editor.swaps_disagreeing += 1;
        }
        editor.last_swap = Some(Swap {
            predicted_commute: predicted,
            hash_changed: changed,
        });
    }
}

fn hud(editor: Res<Editor>, mut stats: ResMut<DemoStats>) {
    stats.title = format!(
        "E-207  editor   log {} of {}   {} triangles",
        editor.cursor,
        editor.log.len(),
        editor.triangles
    );
    stats.vertices = editor.cursor;
    stats.triangles = editor.triangles;

    let ops: String = editor.log[..editor.cursor]
        .iter()
        .map(|b| match b.op {
            BrushOp::Add => '+',
            BrushOp::Subtract => '-',
            _ => '~',
        })
        .collect();

    let mut lines = vec![
        format!(
            "{:<24} {:>8}   ops applied  [Z] undo  [Y] redo  [E] edit",
            "log cursor", editor.cursor
        ),
        format!(
            "{:<24} {:>8}   entries kept for redo",
            "log length",
            editor.log.len()
        ),
        format!("{:<24} {:>8}", "applied", ops),
        format!(
            "{:<24} {:>8.2}   ms to re-fold the world from the log",
            "re-fold", editor.refold_ms
        ),
        String::new(),
        "undo is a re-fold, not a snapshot: the field IS the log, and the".into(),
        "previous state is recomputed rather than restored. a snapshot of a".into(),
        "voxel world is the whole world; a log entry is a shape and an enum.".into(),
        String::new(),
        match editor.round_trip {
            Some(true) => {
                "ROUND TRIP OK -- undo then redo gave a bit-identical world, and the".to_string()
            }
            Some(false) => {
                "!! ROUND TRIP FAILED -- undo/redo did not restore the world".to_string()
            }
            None => "press Z to check the undo/redo round trip".to_string(),
        },
    ];
    if editor.round_trip == Some(true) {
        lines.push("undone state really was different, so the check is not vacuous.".into());
    }

    lines.extend([
        String::new(),
        format!(
            "{:<24} {:>8}   swaps checked against BrushOp::commutes_with  [S]",
            "order", editor.swaps_checked
        ),
        format!(
            "{:<24} {:>8}   said COMMUTE, mesh identical -- correct",
            "  agreed", editor.agreed_commute
        ),
        format!(
            "{:<24} {:>8}   said NO, mesh moved -- correct, and not merely cautious",
            "  agreed", editor.agreed_differ
        ),
        format!(
            "{:<24} {:>8}   said NO, mesh unchanged ({} of them disjoint shapes)",
            "  cautious", editor.conservative, editor.conservative_disjoint
        ),
        format!(
            "{:<24} {:>8}   said COMMUTE and the mesh MOVED -- unsound",
            "  unsound", editor.swaps_disagreeing
        ),
    ]);
    match editor.last_swap {
        Some(swap) => {
            let verdict = match (swap.predicted_commute, swap.hash_changed) {
                (true, false) => "predicted commute, mesh identical -- correct",
                (true, true) => "!! predicted commute and the mesh MOVED -- unsound",
                (false, true) => {
                    "predicted no commute, mesh moved -- correct, and not merely cautious"
                }
                (false, false) => "predicted no commute, mesh unchanged -- conservative here",
            };
            lines.push(format!("{:<24} {:>8}   {verdict}", "last swap", ""));
        }
        None => lines.push(format!(
            "{:<24} {:>8}   press S to swap the last two ops",
            "last swap", "--"
        )),
    }

    lines.extend([
        String::new(),
        "M-36 measured a run of same-kind hard edits reordering bit-for-bit --".into(),
        "one result from all 40,320 orderings. M-37 measured an add/subtract".into(),
        "boundary giving 11 distinct results, a difference that is semantic and".into(),
        "that no storage format repairs. so the log's order is load-bearing and".into(),
        "commutes_with returns the honest answer rather than the optimistic one.".into(),
    ]);
    stats.extra = lines;
}
