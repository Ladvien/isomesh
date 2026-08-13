//! E-114 — what Dual Contouring actually operates on.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example hermite_debug --release
//! ```
//!
//! **Always `--release`.**
//!
//! Arrows move the cell in x/y · `,` `.` in z · `-` `=` λ · `[` `]` resolution ·
//! `1`–`3` field · `H` hide the mesh.
//!
//! # One cell, everything the solve knows about it
//!
//! Dual Contouring places one vertex per cell, and the input to that decision is
//! not the field — it is **Hermite data**: the points where the surface crosses
//! the cell's twelve edges, and the surface normal at each of those points. This
//! example draws exactly that, for one cell at a time.
//!
//! - **Grey spheres** — the eight corners. Filled means inside the solid.
//! - **Amber dots** — edge crossings, one per cut edge.
//! - **Amber lines** — the surface normal at each crossing, the *only* directional
//!   information the solve has.
//! - **Green dot** — where the QEF put the vertex.
//! - **Red dot and line** — where it put it *before* the cell clamp, when the
//!   clamp had to drag it back.
//! - **White box** — the cell.
//!
//! Everything else on screen is context: the extracted mesh, drawn faintly, with
//! `H` to remove it.
//!
//! # Why the interesting cell is chosen for you
//!
//! On a 13³ grid there are 1,728 cells and almost all of them are boring — a flat
//! patch of surface whose normals all point the same way, where the QEF has one
//! obvious answer. The cells worth looking at are the ones where the normals
//! *disagree*, because that disagreement is the entire reason Dual Contouring
//! exists and the entire reason its solve can go wrong.
//!
//! So on every field or resolution change this jumps to the cell whose unclamped
//! solution sits furthest from its own centre, measured in cells. That is
//! precisely M-30's quantity, and it lands you on a sharp corner or a
//! near-degenerate cell without hunting. The arrows then walk from there, and
//! `ISOMESH_CELL=x,y,z` pins one for a capture.
//!
//! # The corner order is duplicated here, and checked rather than trusted
//!
//! [`HermiteCell::from_corners`] takes eight corner samples *"in this crate's
//! corner order"* — and that order lives in a private module, so no consumer can
//! read it. The layout is the obvious one, corner `c` at
//! `[c & 1, (c >> 1) & 1, (c >> 2) & 1]`, and `examples/common` already relies on
//! it for the domain box.
//!
//! Relying on it *silently* is what this file will not do. [`check_corner_order`]
//! runs at startup: it builds a cell from a plane whose crossing position is known
//! in closed form, and reports loudly if the crate disagrees. A demo that draws
//! the wrong corner as inside would look entirely plausible and teach the reader
//! something false.

mod common;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::dual_contouring::{CLAMP_EPSILON, DualContouring, solve};
use isomesh::fields::{BoxExact, ReferenceField, Sphere, csg_difference};
use isomesh::hermite::HermiteCell;
use isomesh::{RuntimeShape3, Sdf};

const FIELDS: [&str; 3] = ["box_exact", "csg_difference", "sphere"];

/// Coarse on purpose. This is a demo of *one cell*, and at 33³ a cell is four
/// pixels across.
///
/// **Every value in this range avoids E-104's alignment trap, and the step is
/// what enforces it.** `box_exact` is exactly zero across its whole boundary, so
/// on a grid whose planes land on the box faces the sign convention decides the
/// answer rather than the algorithm. Over the ±2 domain that happens exactly when
/// `n − 1` is a multiple of 4 — and the first draft of this file defaulted to
/// **13**, which is one of them: corner 7 sampled `-0.0000` and the demo opened
/// on the degenerate case it exists to explain. Stepping by 4 from 11 gives
/// `11, 15, 19, 23`, where `n − 1 ≡ 2 (mod 4)` throughout, so no reachable
/// resolution is aligned.
const DEFAULT_SAMPLES: u32 = 11;
const MIN_SAMPLES: u32 = 7;
const MAX_SAMPLES: u32 = 23;
const SAMPLES_STEP: u32 = 4;

const MIN_LAMBDA: f64 = 1.0e-6;
const MAX_LAMBDA: f64 = 1.0;
const LAMBDA_STEP: f64 = 4.0;

/// Debug marks get their own group so they can be thick and drawn in front. An
/// unbiased line lying on the surface z-fights and reads as intermittent —
/// `manifold_check` earned this first.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct HermiteGizmos;

#[derive(Resource)]
struct Demo {
    field: usize,
    samples: u32,
    lambda: f64,
    /// Selected cell, in grid coordinates.
    cell: [u32; 3],
    /// Cleared once the automatic pick has run for this (field, resolution).
    auto_pick: bool,
    show_mesh: bool,
}

#[derive(Resource)]
struct Surface(Handle<StandardMaterial>);

/// Everything drawn for the selected cell, resolved to world space once per
/// change rather than per frame.
#[derive(Resource, Default)]
struct Overlay {
    size: f32,
    /// The eight corners and whether each is inside the solid.
    corners: Vec<(Vec3, bool)>,
    /// Crossing position and its unit normal.
    crossings: Vec<(Vec3, Vec3)>,
    solved: Option<Vec3>,
    /// Present only when the clamp actually moved the vertex.
    unclamped: Option<Vec3>,
}

/// The local offset of corner `c`, as grid steps.
///
/// Duplicated from the core crate's private `cube::corner_offset`, and verified
/// against it at startup by [`check_corner_order`] rather than assumed.
fn corner_offset(c: u8) -> [u32; 3] {
    [
        u32::from(c & 1),
        u32::from((c >> 1) & 1),
        u32::from((c >> 2) & 1),
    ]
}

/// The crate's own sign convention: negative is inside.
fn is_inside(v: f32) -> bool {
    v < 0.0
}

/// Confirm the duplicated corner order still matches the crate's.
///
/// A plane `f(p) = p.x - 0.25` cut across the unit cell crosses the four
/// x-aligned edges at `x = 0.25` and nothing else. If the corner order this file
/// assumes were a permutation of the crate's, the corner *values* handed to
/// [`HermiteCell::from_corners`] would be permuted too, and the crossings would
/// land on different edges at different places. Checking the count and the
/// x-coordinate catches that.
///
/// **Mutation-tested, and the count alone would not have caught it.** Swapping x
/// and y in [`corner_offset`] — the obvious transcription slip — still yields
/// *four* crossings, and is caught only by the position: worst x error `7.5e-1`.
/// A check that counted crossings and stopped would have passed on a demo whose
/// every corner marker was wrong.
fn check_corner_order() {
    struct Plane;
    impl Sdf for Plane {
        type Scalar = f32;
        fn sample(&self, p: [f32; 3]) -> f32 {
            p[0] - 0.25
        }
    }

    let mut values = [0.0f32; 8];
    for (c, slot) in values.iter_mut().enumerate() {
        let o = corner_offset(c as u8);
        *slot = Plane.sample([o[0] as f32, o[1] as f32, o[2] as f32]);
    }
    let cell = HermiteCell::from_corners(&Plane, &values, [0.0; 3], 1.0);

    let count = cell.iter().count();
    let worst = cell
        .iter()
        .map(|c| (c.position[0] - 0.25).abs())
        .fold(0.0f32, f32::max);
    if count != 4 || worst > 1e-5 {
        error!(
            "corner order check FAILED: expected 4 crossings at x = 0.25, got {count} \
             with worst x error {worst:e}. This file's corner_offset no longer matches \
             the crate's, and every corner marker it draws is wrong."
        );
    }
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-114 hermite debug".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<HermiteGizmos>()
        .insert_resource(Demo {
            field: 0,
            samples: common::samples_override().unwrap_or(DEFAULT_SAMPLES),
            lambda: solve::LAMBDA,
            cell: cell_from_env().unwrap_or([0; 3]),
            auto_pick: cell_from_env().is_none(),
            show_mesh: true,
        })
        .init_resource::<Overlay>()
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, rebuild, draw))
        .run();
}

/// `ISOMESH_CELL=x,y,z` pins the cell, so a capture does not depend on the
/// automatic pick and can show a cell chosen for the picture.
fn cell_from_env() -> Option<[u32; 3]> {
    let raw = std::env::var("ISOMESH_CELL").ok()?;
    let mut parts = raw.split(',').map(|p| p.trim().parse::<u32>());
    match (parts.next(), parts.next(), parts.next(), parts.next()) {
        (Some(Ok(x)), Some(Ok(y)), Some(Ok(z)), None) => Some([x, y, z]),
        _ => {
            warn!("ISOMESH_CELL={raw:?} is not x,y,z; ignoring");
            None
        }
    }
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut gizmo_config: ResMut<GizmoConfigStore>,
    mut camera: Query<&mut OrbitCamera>,
    mut demo: ResMut<Demo>,
    flags: Res<ViewFlags>,
) {
    check_corner_order();
    demo.field = flags.field.min(FIELDS.len() - 1);

    let (config, _) = gizmo_config.config_mut::<HermiteGizmos>();
    config.line.width = 3.0;
    config.depth_bias = -0.6;

    for mut orbit in &mut camera {
        orbit.yaw = 0.7;
        orbit.pitch = 0.35;
        orbit.radius = 6.0;
    }
    commands.insert_resource(Surface(materials.add(StandardMaterial {
        base_color: Color::srgba(0.72, 0.76, 0.82, 0.35),
        perceptual_roughness: 0.5,
        // The cell under inspection is usually inside the mesh, so the context
        // surface is translucent and double-sided or it hides its own subject.
        alpha_mode: AlphaMode::Blend,
        double_sided: true,
        cull_mode: None,
        ..default()
    })));
}

fn controls(keys: Res<ButtonInput<KeyCode>>, mut demo: ResMut<Demo>, mut flags: ResMut<ViewFlags>) {
    let limit = demo.samples.saturating_sub(2);
    let mut step = |axis: usize, by: i32| {
        let v = demo.cell[axis] as i32 + by;
        demo.cell[axis] = v.clamp(0, limit as i32) as u32;
        demo.auto_pick = false;
    };
    for (key, axis, by) in [
        (KeyCode::ArrowRight, 0, 1),
        (KeyCode::ArrowLeft, 0, -1),
        (KeyCode::ArrowUp, 1, 1),
        (KeyCode::ArrowDown, 1, -1),
        (KeyCode::Period, 2, 1),
        (KeyCode::Comma, 2, -1),
    ] {
        if keys.just_pressed(key) {
            step(axis, by);
        }
    }
    for (key, index) in [
        (KeyCode::Digit1, 0),
        (KeyCode::Digit2, 1),
        (KeyCode::Digit3, 2),
    ] {
        if keys.just_pressed(key) {
            demo.field = index;
            demo.auto_pick = true;
        }
    }
    if keys.just_pressed(KeyCode::Equal) {
        demo.lambda = (demo.lambda * LAMBDA_STEP).min(MAX_LAMBDA);
    }
    if keys.just_pressed(KeyCode::Minus) {
        demo.lambda = (demo.lambda / LAMBDA_STEP).max(MIN_LAMBDA);
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        demo.samples = (demo.samples + SAMPLES_STEP).min(MAX_SAMPLES);
        demo.auto_pick = true;
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        demo.samples = demo.samples.saturating_sub(SAMPLES_STEP).max(MIN_SAMPLES);
        demo.auto_pick = true;
    }
    if keys.just_pressed(KeyCode::KeyH) {
        demo.show_mesh = !demo.show_mesh;
    }
    if keys.just_pressed(KeyCode::KeyR) {
        flags.remesh_requested = true;
    }
}

/// Reproduces the crate's cell clamp, which is `pub(crate)`.
///
/// Uses the crate's own [`CLAMP_EPSILON`] rather than a number of its own, so
/// the two cannot drift in the part that matters.
fn clamp_to_cell(x: [f32; 3], origin: [f32; 3], size: f32) -> [f32; 3] {
    let half = size * 0.5;
    let inset = half * (1.0 - CLAMP_EPSILON as f32);
    let mut out = x;
    for (axis, slot) in out.iter_mut().enumerate() {
        let centre = origin[axis] + half;
        *slot = slot.clamp(centre - inset, centre + inset);
    }
    out
}

/// The eight corner samples of one cell, in the crate's corner order.
fn corner_values<F: Sdf<Scalar = f32>>(
    field: &F,
    min: [f32; 3],
    cell: [u32; 3],
    h: f32,
) -> ([f32; 8], [f32; 3]) {
    let origin = [
        min[0] + cell[0] as f32 * h,
        min[1] + cell[1] as f32 * h,
        min[2] + cell[2] as f32 * h,
    ];
    let mut values = [0.0f32; 8];
    for (c, slot) in values.iter_mut().enumerate() {
        let o = corner_offset(c as u8);
        *slot = field.sample([
            origin[0] + o[0] as f32 * h,
            origin[1] + o[1] as f32 * h,
            origin[2] + o[2] as f32 * h,
        ]);
    }
    (values, origin)
}

/// The cell whose crossing normals disagree most.
///
/// Scored as `1 - |mean(normals)|`: unit normals that all point the same way
/// average to length 1 and score 0, while the three mutually perpendicular
/// normals of a box corner average to `1/√3` and score `0.42`. So the maximum is
/// a sharp feature, which is the only kind of cell worth drawing — a flat patch
/// has one obvious answer and demonstrates nothing.
///
/// **The obvious metric was tried first and is nearly degenerate here.** Picking
/// the cell whose unclamped solve sits furthest from its own centre is M-30's
/// quantity, but M-30 also records that `box_exact` has *zero* vertices outside
/// their cells — so on this demo's default field the score is ~0.006 everywhere
/// and the winner is chosen by rounding. It happened to land on a corner, which
/// is exactly how a broken heuristic survives review.
fn most_interesting<F: Sdf<Scalar = f32>>(
    field: &F,
    min: [f32; 3],
    samples: u32,
    h: f32,
) -> Option<[u32; 3]> {
    let last = samples.saturating_sub(2);
    let mut best = None;
    let mut worst = -1.0f32;
    for z in 0..=last {
        for y in 0..=last {
            for x in 0..=last {
                let (values, origin) = corner_values(field, min, [x, y, z], h);
                let cell = HermiteCell::from_corners(field, &values, origin, h);
                // Two crossings cannot show a corner, and one cannot disagree
                // with anything.
                if cell.len() < 3 {
                    continue;
                }
                let mut mean = [0.0f32; 3];
                for c in cell.iter() {
                    for (axis, slot) in mean.iter_mut().enumerate() {
                        *slot += c.normal[axis];
                    }
                }
                let n = cell.len() as f32;
                let len = (mean[0] * mean[0] + mean[1] * mean[1] + mean[2] * mean[2]).sqrt() / n;
                let spread = 1.0 - len;
                if spread.is_finite() && spread > worst {
                    worst = spread;
                    best = Some([x, y, z]);
                }
            }
        }
    }
    best
}

/// (field, samples, λ bits, cell, whether the pick is still automatic).
type RebuildKey = (usize, u32, u64, [u32; 3], bool);

struct Built {
    mesh: Mesh,
    vertices: usize,
    triangles: usize,
    overlay: Overlay,
    lines: Vec<String>,
    field_name: &'static str,
    cell: [u32; 3],
}

fn build(demo: &Demo) -> Option<Built> {
    match demo.field {
        0 => inspect(&BoxExact::<f32>::canonical(), demo),
        1 => inspect(&csg_difference::<f32>(), demo),
        _ => inspect(&Sphere::<f32>::canonical(), demo),
    }
}

fn inspect<F: Sdf<Scalar = f32> + ReferenceField>(field: &F, demo: &Demo) -> Option<Built> {
    let (min, max) = field.domain();
    let h = (max[0] - min[0]) / (demo.samples - 1) as f32;
    let shape = match RuntimeShape3::new([demo.samples; 3]) {
        Ok(shape) => shape,
        Err(error) => {
            error!("grid {}^3 rejected: {error}", demo.samples);
            return None;
        }
    };

    let mut builder = MeshBuilder::new();
    let mut dc = DualContouring::<f32>::new();
    dc.set_lambda(Some(demo.lambda));
    if let Err(error) = dc.extract(field, &shape, min, h, &mut builder) {
        error!("extraction failed: {error}");
        return None;
    }

    let picked = if demo.auto_pick {
        most_interesting(field, min, demo.samples, h).unwrap_or(demo.cell)
    } else {
        let limit = demo.samples.saturating_sub(2);
        [
            demo.cell[0].min(limit),
            demo.cell[1].min(limit),
            demo.cell[2].min(limit),
        ]
    };

    let (values, origin) = corner_values(field, min, picked, h);
    let cell = HermiteCell::from_corners(field, &values, origin, h);

    let corners = (0..8u8)
        .map(|c| {
            let o = corner_offset(c);
            (
                Vec3::new(
                    origin[0] + o[0] as f32 * h,
                    origin[1] + o[1] as f32 * h,
                    origin[2] + o[2] as f32 * h,
                ),
                is_inside(values[c as usize]),
            )
        })
        .collect();
    let crossings: Vec<(Vec3, Vec3)> = cell
        .iter()
        .map(|c| (Vec3::from(c.position), Vec3::from(c.normal)))
        .collect();

    let raw = solve::solve_with(&cell, demo.lambda as f32);
    let held = raw.map(|v| clamp_to_cell(v, origin, h));
    // Only worth drawing when the clamp actually did something. A vertex nudged
    // by 5e-5 of a cell is not a runaway, and drawing it as one would be the
    // bitwise-count mistake A-009 already made once.
    let moved = match (raw, held) {
        (Some(r), Some(c)) => {
            let d = [r[0] - c[0], r[1] - c[1], r[2] - c[2]];
            let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() / h;
            (dist > 1e-3).then_some((Vec3::from(r), dist))
        }
        _ => None,
    };

    let inside = values.iter().filter(|v| is_inside(**v)).count();
    let centre = [
        origin[0] + h * 0.5,
        origin[1] + h * 0.5,
        origin[2] + h * 0.5,
    ];
    let away = held.map_or(f32::NAN, |v| {
        let d = [v[0] - centre[0], v[1] - centre[1], v[2] - centre[2]];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() / h
    });

    let mut lines = vec![
        format!(
            "{:<22} {:>8}   arrows move x/y, , and . move z",
            "cell",
            format!("{},{},{}", picked[0], picked[1], picked[2])
        ),
        format!("{:<22} {:>8.4}   cell size", "h", h),
        format!(
            "{:<22} {:>8}   of 8, so this cell is {}",
            "corners inside",
            inside,
            match inside {
                0 | 8 => "EMPTY -- the surface misses it",
                _ => "cut by the surface",
            }
        ),
        format!(
            "{:<22} {:>8}   of 12 edges (the amber dots)",
            "crossings",
            crossings.len()
        ),
        String::new(),
    ];
    match held {
        Some(v) => lines.push(format!(
            "{:<22} {:>8}   {:.3} cells from the cell centre",
            "solved vertex",
            format!("{:.2},{:.2},{:.2}", v[0], v[1], v[2]),
            away
        )),
        None => lines.push(format!(
            "{:<22} {:>8}   no crossings, so nothing to solve",
            "solved vertex", "--"
        )),
    }
    match moved {
        Some((_, dist)) => lines.push(format!(
            "{:<22} {:>8.3}   cells the CLAMP had to drag it back (red)",
            "unclamped overshoot", dist
        )),
        None => lines.push(format!(
            "{:<22} {:>8}   the clamp did not bind here",
            "unclamped overshoot", "--"
        )),
    }
    lines.extend([
        format!(
            "{:<22} {:>8.0e}   [-] softer, [=] sharper",
            "lambda", demo.lambda
        ),
        String::new(),
        "the eight corner samples, in the crate's corner order:".into(),
    ]);
    for c in 0..8u8 {
        let o = corner_offset(c);
        lines.push(format!(
            "   corner {c} ({},{},{})   {:>9.4}   {}",
            o[0],
            o[1],
            o[2],
            values[c as usize],
            if is_inside(values[c as usize]) {
                "inside"
            } else {
                "outside"
            }
        ));
    }
    lines.extend([
        String::new(),
        "this is the whole input to the vertex decision: crossings, and the".into(),
        "normal at each one. the QEF finds the point minimising squared distance".into(),
        "to every crossing's tangent plane, which is why disagreeing normals".into(),
        "produce a sharp corner and agreeing ones produce a flat patch.".into(),
        String::new(),
        format!(
            "grid alignment           {:>8}   (n-1) mod 4 = {}; on box_exact a grid",
            if (demo.samples - 1) % 4 == 0 {
                "ALIGNED"
            } else {
                "safe"
            },
            (demo.samples - 1) % 4
        ),
        "                                    aligned to the box faces lets the sign".into(),
        "                                    convention decide instead of the solve (E-104)."
            .into(),
        String::new(),
        format!(
            "{} samples/axis   [ and ] to change   1-3 field   H hides the mesh",
            demo.samples
        ),
    ]);

    Some(Built {
        overlay: Overlay {
            size: h,
            corners,
            crossings,
            solved: held.map(Vec3::from),
            unclamped: moved.map(|(v, _)| v),
        },
        lines,
        field_name: F::NAME,
        cell: picked,
        vertices: builder.vertex_count(),
        triangles: builder.triangle_count(),
        mesh: builder.into_mesh(),
    })
}

#[allow(clippy::too_many_arguments)]
fn rebuild(
    mut demo: ResMut<Demo>,
    mut stats: ResMut<DemoStats>,
    mut overlay: ResMut<Overlay>,
    mut meshes: ResMut<Assets<Mesh>>,
    surface: Res<Surface>,
    mut commands: Commands,
    mut query: Query<(&mut Mesh3d, &mut Visibility), With<DemoMesh>>,
    mut flags: ResMut<ViewFlags>,
    mut last: Local<Option<RebuildKey>>,
) {
    let key = (
        demo.field,
        demo.samples,
        demo.lambda.to_bits(),
        demo.cell,
        demo.auto_pick,
    );
    if *last == Some(key) && !flags.remesh_requested {
        // The mesh does not change when only its visibility does.
        for (_, mut visible) in &mut query {
            *visible = if demo.show_mesh {
                Visibility::Inherited
            } else {
                Visibility::Hidden
            };
        }
        return;
    }
    *last = Some(key);
    flags.remesh_requested = false;

    let Some(built) = build(&demo) else {
        return;
    };
    // Adopt the automatic pick so the arrows walk from it rather than from
    // wherever the cursor happened to be.
    demo.cell = built.cell;
    demo.auto_pick = false;

    stats.title = format!(
        "E-114  hermite debug   cell {},{},{}   field {} ({})   {}^3",
        built.cell[0],
        built.cell[1],
        built.cell[2],
        demo.field + 1,
        built.field_name,
        demo.samples,
    );
    stats.vertices = built.vertices;
    stats.triangles = built.triangles;
    stats.extra = built.lines;
    *overlay = built.overlay;

    let handle = meshes.add(built.mesh);
    if query.is_empty() {
        commands.spawn((
            Mesh3d(handle),
            MeshMaterial3d(surface.0.clone()),
            Transform::default(),
            DemoMesh,
        ));
    } else {
        for (mut mesh, _) in &mut query {
            mesh.0 = handle.clone();
        }
    }
}

fn draw(overlay: Res<Overlay>, mut gizmos: Gizmos<HermiteGizmos>) {
    const WHITE: Color = Color::srgb(0.95, 0.95, 0.95);
    const AMBER: Color = Color::srgb(1.0, 0.70, 0.15);
    const GREEN: Color = Color::srgb(0.30, 0.95, 0.45);
    const RED: Color = Color::srgb(1.0, 0.13, 0.13);
    const GREY: Color = Color::srgb(0.55, 0.58, 0.65);

    let h = overlay.size;
    if h <= 0.0 {
        return;
    }

    // The cell, drawn from its own corners rather than as a primitive: two
    // corners share an edge exactly when their indices differ in one bit, which
    // is the same xyz bit layout `corner_offset` uses. Deriving the box from the
    // corners keeps one convention on screen instead of two.
    for a in 0..overlay.corners.len() {
        for bit in 0..3 {
            let b = a ^ (1 << bit);
            if b > a && b < overlay.corners.len() {
                gizmos.line(overlay.corners[a].0, overlay.corners[b].0, WHITE);
            }
        }
    }

    for (p, inside) in &overlay.corners {
        // Inside corners get a filled-looking sphere, outside ones a sparse
        // outline, so the sign pattern reads at a glance without the HUD.
        let (colour, resolution) = if *inside { (WHITE, 8) } else { (GREY, 3) };
        gizmos
            .sphere(Isometry3d::from_translation(*p), h * 0.06, colour)
            .resolution(resolution);
    }

    for (p, n) in &overlay.crossings {
        gizmos
            .sphere(Isometry3d::from_translation(*p), h * 0.05, AMBER)
            .resolution(6);
        // Scaled to the cell, not to the world: the normal is a direction and
        // its length here means nothing except legibility.
        gizmos.line(*p, *p + *n * h * 0.45, AMBER);
    }

    if let Some(v) = overlay.solved {
        gizmos
            .sphere(Isometry3d::from_translation(v), h * 0.10, GREEN)
            .resolution(10);
        if let Some(raw) = overlay.unclamped {
            gizmos
                .sphere(Isometry3d::from_translation(raw), h * 0.08, RED)
                .resolution(8);
            gizmos.line(raw, v, RED);
        }
    }
}
