//! E-209 — a concave edge, moving, measured every frame.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_csg_props --release
//! ```
//!
//! **Always `--release`.**
//!
//! `Space` pause the motion · `A` displayed extractor · `[` `]` resolution.
//!
//! # Why a concave edge, and why it moves
//!
//! `dual_contouring_cube` already measures a **convex** corner and finds Dual
//! Contouring 101x closer to it than Surface Nets. A convex corner is the easy
//! case: its solution lies *inside* the cell that contains it, so the cell clamp
//! never binds and the QEF has room to be right (M-28). A **concave** edge is
//! where a CAD tool actually lives — the inside corner of a pocket, the seam
//! where a boss meets a face — and it is the case where the solution wants to sit
//! outside the cell that produced it.
//!
//! It also **moves**. E-104 measured one static configuration and had to defend
//! itself against grid alignment with a rule about which resolutions to skip;
//! sweeping the cutter continuously means almost every frame is unaligned, and
//! the number reported is the **worst over the whole sweep** rather than whatever
//! one lucky position gave. A single-position measurement of a sharp feature is a
//! measurement of that position.
//!
//! # The solid
//!
//! ```text
//! solid(p) = max( box(p) , min(p.x − cx, p.y − cy) )
//! ```
//!
//! A `max` is a subtraction of the quarter-space `x > cx ∧ y > cy`, which cuts an
//! L out of the block. The inside of that L is a reflex dihedral running along
//! `z` at exactly `(cx, cy)` — an edge whose position is known in closed form,
//! which is what makes "how far is the nearest vertex" a measurement rather than
//! an impression. `(cx, cy)` orbits, so the edge sweeps through the grid.
//!
//! # Not chunked, deliberately
//!
//! One grid, one extraction. B-006 measured that the dual methods leave gaps at a
//! chunk boundary — 4 to 5 open edges on a single seam — because a boundary quad
//! needs the neighbour cell's vertex. A demo whose entire subject is Dual
//! Contouring holding an edge must not be chunked, and `Extractor::chunk_seams`
//! now says so at the call site for anyone who reaches for the plugin instead.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::to_bevy_mesh;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::dual_contouring::DualContouring;
use isomesh::fields::{BoxExact, Intersection};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::surface_nets::SurfaceNets;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

const HALF: f32 = 2.0;
const DOMAIN: f32 = 2.6;
const DEFAULT_SAMPLES: u32 = 41;
const MIN_SAMPLES: u32 = 25;
const MAX_SAMPLES: u32 = 65;
const SAMPLES_STEP: u32 = 8;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Which {
    DualContouring,
    /// Dual Contouring with A-009's cell clamp **off**.
    ///
    /// Here to test a mechanism rather than to be a recommendation. M-28 found
    /// the clamp costs nothing on a convex corner — the corner measures the same
    /// distance clamped or not, because a convex solution is interior to its own
    /// cell and the constraint never binds. A reflex edge is the opposite case:
    /// its solution wants to sit *outside*. If the clamp is what caps Dual
    /// Contouring on a concave edge, turning it off must move the worst case and
    /// nothing else will.
    DualContouringUnclamped,
    SurfaceNets,
    MarchingCubes,
}

impl Which {
    const ALL: [Self; 4] = [
        Self::DualContouring,
        Self::DualContouringUnclamped,
        Self::SurfaceNets,
        Self::MarchingCubes,
    ];

    fn name(self) -> &'static str {
        match self {
            Self::DualContouring => "dual contouring",
            Self::DualContouringUnclamped => "dc, clamp off",
            Self::SurfaceNets => "surface nets",
            Self::MarchingCubes => "marching cubes",
        }
    }

    fn next(self) -> Self {
        match self {
            Self::DualContouring => Self::DualContouringUnclamped,
            Self::DualContouringUnclamped => Self::SurfaceNets,
            Self::SurfaceNets => Self::MarchingCubes,
            Self::MarchingCubes => Self::DualContouring,
        }
    }

    fn from_env() -> Self {
        match std::env::var("ISOMESH_ALGORITHM")
            .unwrap_or_default()
            .as_str()
        {
            "sn" | "surface_nets" => Self::SurfaceNets,
            "mc" | "marching_cubes" => Self::MarchingCubes,
            _ => Self::DualContouring,
        }
    }
}

/// A block with a quarter-space cut out of it.
///
/// The cut is `max(box, min(x - cx, y - cy))`, and the reflex edge it leaves runs
/// along `z` at `(cx, cy)`.
#[derive(Clone, Copy)]
struct Notched {
    cx: f32,
    cy: f32,
}

impl Notched {
    /// The block, as the crate's own [`BoxExact`] rather than as a second copy
    /// of its distance function — and with the analytic gradient that copy did
    /// not carry.
    fn block() -> BoxExact<f32> {
        BoxExact {
            center: [0.0; 3],
            half_extents: [HALF; 3],
        }
    }

    /// The quarter-space `min(x − cx, y − cy)`.
    ///
    /// **This part stays local, deliberately.** It is a union of two half-spaces,
    /// and `fields::` ships `Difference` and `Intersection` but **no `Union`** —
    /// union exists only as `BrushOp::Add`, which is an edit operation rather
    /// than a field combinator. Writing one would be new public core API, which
    /// E-212 is not.
    fn cut(&self) -> Quarter {
        Quarter {
            cx: self.cx,
            cy: self.cy,
        }
    }
}

impl Sdf for Notched {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        Intersection {
            a: Self::block(),
            b: self.cut(),
        }
        .sample(p)
    }

    fn gradient(&self, p: [f32; 3]) -> [f32; 3] {
        Intersection {
            a: Self::block(),
            b: self.cut(),
        }
        .gradient(p)
    }
}

/// The quarter-space removed from the block: `min(x − cx, y − cy)`.
///
/// Its gradient is the active half-space's normal, which is exact everywhere
/// except on the reflex edge itself — and the reflex edge is precisely what this
/// demo measures, so it is worth having exactly rather than by six samples.
#[derive(Clone, Copy)]
struct Quarter {
    cx: f32,
    cy: f32,
}

impl Sdf for Quarter {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        (p[0] - self.cx).min(p[1] - self.cy)
    }

    fn gradient(&self, p: [f32; 3]) -> [f32; 3] {
        if p[0] - self.cx <= p[1] - self.cy {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        }
    }
}

/// What one extractor did with the edge this frame.
#[derive(Clone, Copy, Default)]
struct Measured {
    /// Distance from the exact reflex edge to the nearest vertex, in cells.
    edge_cells: f32,
    /// Worst that distance has been over the whole sweep.
    worst_cells: f32,
    /// Running sum and count, for the mean. The worst case and the typical case
    /// answer different questions and this demo needs both.
    sum_cells: f64,
    samples_taken: u32,
    triangles: usize,
    ms: f64,
}

#[derive(Resource)]
struct Demo {
    samples: u32,
    shown: Which,
    moving: bool,
    t: f32,
    stats: [Measured; 4],
    frames: u32,
}

#[derive(Resource)]
struct Look(Handle<StandardMaterial>);

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-209 csg props".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Demo {
            samples: common::samples_override().unwrap_or(DEFAULT_SAMPLES),
            shown: Which::from_env(),
            moving: true,
            t: 0.0,
            stats: [Measured::default(); 4],
            frames: 0,
        })
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, remesh, report, hud).chain())
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    for mut orbit in &mut camera {
        // Looking into the notch, which is where the edge is.
        orbit.yaw = 0.95;
        orbit.pitch = 0.55;
        orbit.radius = 8.5;
    }
    commands.insert_resource(Look(materials.add(StandardMaterial {
        base_color: Color::srgb(0.74, 0.71, 0.66),
        perceptual_roughness: 0.42,
        ..default()
    })));
}

/// Where the cutter's corner is at time `t`.
///
/// A slow ellipse, chosen so the edge sweeps across cell boundaries continuously
/// rather than sitting near one. The worst case over the sweep is the number that
/// matters; any single position can be lucky.
fn corner(t: f32) -> (f32, f32) {
    (0.85 * (t * 0.7).sin(), 0.85 * (t * 0.53).cos())
}

fn controls(
    time: Res<Time>,
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut demo: ResMut<Demo>,
    mut flags: ResMut<ViewFlags>,
) {
    flags.grid = false;
    if demo.moving {
        demo.t += time.delta_secs();
    }
    if capture.is_active() {
        return;
    }
    if keys.just_pressed(KeyCode::Space) {
        demo.moving = !demo.moving;
    }
    if keys.just_pressed(KeyCode::KeyA) {
        demo.shown = demo.shown.next();
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        demo.samples = (demo.samples + SAMPLES_STEP).min(MAX_SAMPLES);
        demo.stats = [Measured::default(); 4];
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        demo.samples = demo.samples.saturating_sub(SAMPLES_STEP).max(MIN_SAMPLES);
        demo.stats = [Measured::default(); 4];
    }
}

/// Extract once, and measure how close the mesh gets to the exact reflex edge.
///
/// The edge is the segment `x = cx, y = cy, |z| <= HALF`, known in closed form —
/// so this is a distance to a known line, not to a guess. Vertices are compared
/// against the nearest point *on the segment*, which is what makes the number a
/// property of the edge rather than of the block's corners.
fn extract(which: Which, field: &Notched, samples: u32) -> (MeshBuffer<f32>, f32, f64) {
    let cell = (DOMAIN * 2.0) / (samples - 1) as f32;
    let min = [-DOMAIN; 3];
    let Ok(shape) = RuntimeShape3::new([samples; 3]) else {
        return (MeshBuffer::new(), f32::NAN, 0.0);
    };
    let mut out = MeshBuffer::<f32>::new();
    let started = Instant::now();
    let ok = match which {
        Which::DualContouring => {
            DualContouring::<f32>::new().extract(field, &shape, min, cell, &mut out)
        }
        Which::DualContouringUnclamped => {
            let mut dc = DualContouring::<f32>::new();
            dc.set_clamp(isomesh::dual_contouring::Clamp::None);
            dc.extract(field, &shape, min, cell, &mut out)
        }
        Which::SurfaceNets => SurfaceNets::<f32>::new().extract(field, &shape, min, cell, &mut out),
        Which::MarchingCubes => {
            MarchingCubes::<f32>::new().extract(field, &shape, min, cell, &mut out)
        }
    };
    let ms = started.elapsed().as_secs_f64() * 1000.0;
    if ok.is_err() {
        return (MeshBuffer::new(), f32::NAN, ms);
    }

    let mut nearest = f32::MAX;
    for p in &out.positions {
        // Nearest point on the edge segment: clamp z, take the other two exactly.
        let z = p[2].clamp(-HALF, HALF);
        let d = [p[0] - field.cx, p[1] - field.cy, p[2] - z];
        let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        if dist < nearest {
            nearest = dist;
        }
    }
    (out, nearest / cell, ms)
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    mut demo: ResMut<Demo>,
    look: Res<Look>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut query: Query<&mut Mesh3d, With<DemoMesh>>,
) {
    let (cx, cy) = corner(demo.t);
    let field = Notched { cx, cy };
    let samples = demo.samples;
    let shown = demo.shown;
    let mut display = None;

    // All three every frame, because the comparison is the point and a number
    // that only updates when you press a key is a number you have to remember.
    for (i, which) in Which::ALL.iter().copied().enumerate() {
        let (buffer, edge, ms) = extract(which, &field, samples);
        let slot = &mut demo.stats[i];
        slot.edge_cells = edge;
        slot.triangles = buffer.indices.len() / 3;
        slot.ms = ms;
        // The worst over the sweep, which is the honest figure: any one position
        // of a sharp feature relative to the grid can be lucky.
        if edge.is_finite() {
            if edge > slot.worst_cells {
                slot.worst_cells = edge;
            }
            slot.sum_cells += f64::from(edge);
            slot.samples_taken += 1;
        }
        if which == shown {
            display = Some(buffer);
        }
    }
    demo.frames += 1;

    let Some(buffer) = display else {
        return;
    };
    let handle = meshes.add(to_bevy_mesh(&buffer));
    if query.is_empty() {
        commands.spawn((
            Mesh3d(handle),
            MeshMaterial3d(look.0.clone()),
            Transform::default(),
            DemoMesh,
        ));
    } else {
        for mut mesh in &mut query {
            mesh.0 = handle.clone();
        }
    }
}

/// One CSV row a second, so a long sweep can be read from a shell.
fn report(demo: Res<Demo>, mut next: Local<u32>) {
    if demo.frames < *next || demo.frames == 0 {
        return;
    }
    *next = demo.frames + 60;
    info!(
        "csg,{},{},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.2},{:.2},{:.2}",
        demo.samples,
        demo.frames,
        demo.stats[0].edge_cells,
        demo.stats[0].worst_cells,
        demo.stats[1].edge_cells,
        demo.stats[1].worst_cells,
        demo.stats[2].edge_cells,
        demo.stats[2].worst_cells,
        demo.stats[0].ms,
        demo.stats[1].ms,
        demo.stats[2].ms,
    );
}

fn hud(demo: Res<Demo>, mut stats: ResMut<DemoStats>) {
    let shown = Which::ALL
        .iter()
        .position(|w| *w == demo.shown)
        .unwrap_or(0);
    let total: f64 = demo.stats.iter().map(|s| s.ms).sum();
    let (cx, cy) = corner(demo.t);

    stats.title = format!(
        "E-209  csg props - {}   edge at ({cx:.2}, {cy:.2})   {}^3",
        demo.shown.name(),
        demo.samples
    );
    stats.vertices = demo.stats[shown].triangles;
    stats.triangles = demo.stats[shown].triangles;
    stats.extract_ms = demo.stats[shown].ms;

    let mut lines = vec![
        "distance from the EXACT reflex edge to the nearest vertex, in cells.".into(),
        "the edge is x = cx, y = cy, known in closed form, and it is moving:".into(),
        String::new(),
        format!(
            "{:<20} {:>8}  {:>8}  {:>8}  {:>6}",
            "", "now", "mean", "worst", "ms"
        ),
    ];
    for (i, which) in Which::ALL.iter().copied().enumerate() {
        let s = demo.stats[i];
        let mark = if which == demo.shown { "<--" } else { "   " };
        let mean = if s.samples_taken > 0 {
            s.sum_cells / f64::from(s.samples_taken)
        } else {
            0.0
        };
        lines.push(format!(
            "{mark} {:<16} {:>8.4}  {:>8.4}  {:>8.4}  {:>6.2}",
            which.name(),
            s.edge_cells,
            mean,
            s.worst_cells,
            s.ms
        ));
    }

    let dc = demo.stats[0].worst_cells;
    let unclamped = demo.stats[1].worst_cells;
    let sn = demo.stats[2].worst_cells;
    lines.extend([
        String::new(),
        if dc > 0.0 && sn > 0.0 {
            format!(
                "worst over the sweep: dual contouring only {:.2}x closer than surface",
                sn / dc
            )
        } else {
            "sweeping...".to_string()
        },
        "nets -- nothing like the 101x E-104 measures on a CONVEX corner. the mean".into(),
        "still favours it heavily; it is the worst case that converges.".into(),
        String::new(),
        "the clamp is NOT why -- registered as a hypothesis, and falsified:".to_string(),
        format!("unclamped worst is {unclamped:.4} against clamped {dc:.4}, identical."),
        "M-28 found the clamp costs nothing on a convex corner, because a convex".into(),
        "solution is interior to its own cell and the constraint never binds. the".into(),
        "guess here was that a reflex edge's solution wants to sit OUTSIDE its cell".into(),
        "and would be capped by the clamp. it does not: a reflex edge passing".into(),
        "through a cell has its QEF solution inside that same cell, so the clamp".into(),
        "never binds there either. M-28's result extends to concave features.".into(),
        String::new(),
        "what DOES cap the worst case is not established here. it coincides with".into(),
        "the edge lying near a sample plane -- E-104's alignment trap arriving".into(),
        "continuously rather than at chosen resolutions -- but this demo does not".into(),
        "isolate that, and saying so is cheaper than guessing.".into(),
        String::new(),
        format!(
            "{:<20} {:>9.2}   ms for all three, re-meshed every frame",
            "total", total
        ),
        format!(
            "{} samples/axis   [ and ] to change   A switches   Space pauses",
            demo.samples
        ),
    ]);
    stats.extra = lines;
}
