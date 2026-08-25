//! E-109 — the sharpness knob, and what it costs at both ends.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example sharp_features --release
//! ```
//!
//! **Always `--release`.**
//!
//! `-` `=` λ · `C` clamp · `1`–`4` field · `[` `]` resolution · `W` wireframe.
//!
//! # The parameter the ticket asked for, under its real name
//!
//! E-109 was written asking for a slider on a *normal-deviation threshold*, and
//! this implementation of Dual Contouring has no such parameter. What actually
//! trades sharpness against stability here is **λ, the Tikhonov regularizer** in
//! the vertex solve — a number that was a compile-time constant until this
//! example needed to turn it.
//!
//! Toward **zero**, the solve is the unregularized plane intersection. A corner
//! where three planes meet comes out exactly, which is the entire reason to run
//! Dual Contouring rather than Surface Nets. But a *flat* cell has a rank-1
//! system, nothing determines its vertex along two directions, and with no
//! regularizer pulling it back it leaves — M-30 measured an unclamped solve
//! flinging a vertex **3.18 cells** out of its own cell on `gyroid`.
//!
//! Toward **large**, every vertex is pulled to the centroid of its cell's
//! crossings. Nothing flies anywhere and every sharp edge rounds over, which is
//! Surface Nets with extra arithmetic.
//!
//! The default, `0.01`, sits near the bottom of the usable range on purpose.
//!
//! # Why the clamp toggle is here too
//!
//! Because with the clamp **on** you cannot see the failure. A-009's cell clamp
//! confines each vertex to its own cell, so the runaway that small λ causes is
//! caught before it becomes geometry — the mesh degrades quietly instead of
//! spiking. Turn the clamp off to see what λ is actually protecting you from,
//! and turn it back on to see why the crate defaults to it.
//!
//! That is the honest version of "over-sharpening into spikes": the spikes are
//! real, and one line of defence already stops them.
//!
//! # Two numbers, because the two failures are not the same failure
//!
//! **Worst |f| / h** measures the *rounding* end: how far λ has pulled the
//! furthest vertex off the surface, in cells.
//!
//! **Worst clamp move** measures the *runaway* end: how far the clamp has to drag
//! a vertex back into its own cell.
//!
//! The second one exists because the first is blind to the runaway, which is
//! worth knowing before trusting either. A flat cell is rank 1 and its
//! unconstrained directions lie **within the surface** — so a vertex with nothing
//! holding it slides along the plane and stays exactly on it. At λ = 1e-6 on
//! `box_exact`, `|f| / h` reads **0.000**, which looks like a perfect mesh and is
//! a vertex several cells from where it belongs.

mod common;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::dual_contouring::{Clamp, DualContouring};
use isomesh::fields::{BoxExact, FbmTerrain, ReferenceField, capped_gyroid, csg_difference};
use isomesh::{RuntimeShape3, Sdf};

/// The field choice is not decoration: **the two failure modes live on different
/// fields**, and M-30 says which. The runaway was measured at 3.18 cells on
/// `gyroid` and 2.17 on `fbm_terrain`, and explicitly *not* on the smooth closed
/// fields — *"sphere, box_exact and thin_plate have zero vertices outside"*.
///
/// So `box_exact` is where you watch λ round a corner over, and `gyroid` is where
/// you watch it let a vertex leave. Opening on `box_exact` and sweeping λ shows
/// only half the story, which is exactly the mistake this list is arranged to
/// stop.
const FIELDS: [&str; 4] = ["box_exact", "csg_difference", "gyroid", "fbm_terrain"];

/// λ is swept geometrically — the interesting range spans four decades, and a
/// linear slider would spend all of its travel in the region where nothing
/// changes.
const MIN_LAMBDA: f64 = 1.0e-6;
const MAX_LAMBDA: f64 = 1.0;
const LAMBDA_STEP: f64 = 2.0;

#[derive(Resource)]
struct Demo {
    field: usize,
    samples: u32,
    lambda: f64,
    clamp: bool,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-109 sharp features".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Demo {
            field: 0,
            samples: common::samples_override().unwrap_or(25),
            lambda: 0.01,
            // Off by default here, and *only* here: this example exists to show
            // the failure the clamp hides.
            clamp: false,
        })
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, remesh))
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
    mut demo: ResMut<Demo>,
    flags: Res<ViewFlags>,
) {
    // `ISOMESH_FIELD` is the harness contract for choosing a field without a
    // keyboard, and it is what makes a capture of a *particular* field
    // reproducible. Honoured here rather than only by the digit keys.
    demo.field = flags.field.min(FIELDS.len() - 1);

    for mut orbit in &mut camera {
        orbit.yaw = 0.7;
        orbit.pitch = 0.35;
        orbit.radius = 7.0;
    }
    commands.insert_resource(Surface(materials.add(StandardMaterial {
        base_color: Color::srgb(0.80, 0.78, 0.72),
        perceptual_roughness: 0.45,
        double_sided: true,
        cull_mode: None,
        ..default()
    })));
}

#[derive(Resource)]
struct Surface(Handle<StandardMaterial>);

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut demo: ResMut<Demo>,
    mut flags: ResMut<ViewFlags>,
) {
    // A capture sweeps λ itself, in step with the frames rather than the clock.
    if capture.is_active() {
        const STEPS: u32 = 24;
        let phase = capture.taken % (STEPS * 2);
        let step = if phase < STEPS {
            phase
        } else {
            STEPS * 2 - phase - 1
        };
        let t = f64::from(step) / f64::from(STEPS - 1);
        demo.lambda = MIN_LAMBDA * (MAX_LAMBDA / MIN_LAMBDA).powf(t);
        return;
    }
    for (key, index) in [
        (KeyCode::Digit1, 0),
        (KeyCode::Digit2, 1),
        (KeyCode::Digit3, 2),
        (KeyCode::Digit4, 3),
    ] {
        if keys.just_pressed(key) {
            demo.field = index;
        }
    }
    if keys.just_pressed(KeyCode::Equal) {
        demo.lambda = (demo.lambda * LAMBDA_STEP).min(MAX_LAMBDA);
    }
    if keys.just_pressed(KeyCode::Minus) {
        demo.lambda = (demo.lambda / LAMBDA_STEP).max(MIN_LAMBDA);
    }
    if keys.just_pressed(KeyCode::KeyC) {
        demo.clamp = !demo.clamp;
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        demo.samples = (demo.samples + 4).min(49);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        demo.samples = demo.samples.saturating_sub(4).max(9);
    }
    if keys.just_pressed(KeyCode::KeyR) {
        flags.remesh_requested = true;
    }
}

/// Extract, and measure both failure modes.
///
/// Returns the mesh, the worst **rounding** (a vertex's `|f| / h`, which grows as
/// λ pulls vertices off the surface toward the centroid) and the worst
/// **runaway** (how far the clamp had to move a vertex, in cells, which grows as
/// λ → 0 lets an under-determined cell's vertex leave).
///
/// Two numbers because the two failures are not the same failure and one metric
/// cannot see both. `|f| / h` is blind to the runaway: a flat cell is rank 1 and
/// its *unconstrained directions lie within the surface*, so a vertex with
/// nothing holding it slides along the plane and stays exactly on it — measured
/// **0.000** at λ = 1e-6 on `box_exact`, which reads as "perfect" and is not.
/// M-30's quantity is distance out of the cell, so that is what the second
/// number is: the same extraction run with the clamp on, differenced.
fn extract<F: Sdf<Scalar = f32> + ReferenceField>(
    field: &F,
    demo: &Demo,
) -> Option<(MeshBuilder, f32, f32)> {
    let (min, max) = field.domain();
    let cell = (max[0] - min[0]) / (demo.samples - 1) as f32;
    let shape = RuntimeShape3::new([demo.samples; 3]).ok()?;

    let mut dc = DualContouring::<f32>::new();
    dc.set_lambda(Some(demo.lambda));
    dc.set_clamp(if demo.clamp {
        Clamp::ToCell
    } else {
        Clamp::None
    });

    let mut out = MeshBuilder::new();
    dc.extract(field, &shape, min, cell, &mut out).ok()?;

    // How far off the surface the worst vertex is, in cells. A vertex that is
    // still on the surface reads ~0 whatever λ is; one that has flown off reads
    // how far it flew, and it moves long before the mesh looks wrong.
    let mut rounding = 0.0f32;
    for p in out.positions() {
        let off = (field.sample(*p) / cell).abs();
        if off.is_finite() && off > rounding {
            rounding = off;
        }
    }

    // The same extraction with the clamp on. Clamping moves positions and never
    // changes connectivity, so the vertices correspond index for index and the
    // difference is exactly how far each one had left its cell.
    let mut clamped = DualContouring::<f32>::new();
    clamped.set_lambda(Some(demo.lambda));
    clamped.set_clamp(Clamp::ToCell);
    let mut reference = MeshBuilder::new();
    clamped
        .extract(field, &shape, min, cell, &mut reference)
        .ok()?;

    let mut runaway = 0.0f32;
    if reference.vertex_count() == out.vertex_count() {
        for (loose, held) in out.positions().iter().zip(reference.positions()) {
            let d = [loose[0] - held[0], loose[1] - held[1], loose[2] - held[2]];
            let moved = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() / cell;
            if moved.is_finite() && moved > runaway {
                runaway = moved;
            }
        }
    }
    Some((out, rounding, runaway))
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    demo: Res<Demo>,
    mut stats: ResMut<DemoStats>,
    mut meshes: ResMut<Assets<Mesh>>,
    surface: Res<Surface>,
    mut commands: Commands,
    mut query: Query<&mut Mesh3d, With<DemoMesh>>,
    mut flags: ResMut<ViewFlags>,
    mut last: Local<Option<(usize, u32, u64, bool)>>,
) {
    let key = (demo.field, demo.samples, demo.lambda.to_bits(), demo.clamp);
    if *last == Some(key) && !flags.remesh_requested {
        return;
    }
    *last = Some(key);
    flags.remesh_requested = false;

    let built = match demo.field {
        0 => extract(&BoxExact::<f32>::canonical(), &demo),
        1 => extract(&csg_difference::<f32>(), &demo),
        2 => extract(&capped_gyroid::<f32>(), &demo),
        _ => extract(&FbmTerrain::<f32>::canonical(), &demo),
    };
    let Some((builder, rounding, runaway)) = built else {
        return;
    };

    let verdict = if runaway > 1.0 {
        "RUNAWAY -- a vertex has left its own cell entirely"
    } else if rounding > 0.15 {
        "rounded -- lambda is pulling vertices off the surface"
    } else if runaway > 0.1 {
        "sharp, but vertices are starting to wander"
    } else {
        "sharp and stable -- the useful range"
    };

    stats.title = format!(
        "E-109  sharp features   lambda {:.0e}   clamp {}   field {} ({})",
        demo.lambda,
        if demo.clamp { "on" } else { "OFF" },
        demo.field + 1,
        FIELDS[demo.field]
    );
    stats.vertices = builder.vertex_count();
    stats.triangles = builder.triangle_count();
    stats.extra = vec![
        format!(
            "{:<18} {:>10.0e}   [-] softer, [=] sharper",
            "lambda", demo.lambda
        ),
        format!(
            "{:<18} {:>10}   [C] toggles -- the crate defaults to on",
            "cell clamp",
            if demo.clamp { "on" } else { "OFF" }
        ),
        format!(
            "{:<18} {:>10.3}   cells off the surface -- the ROUNDING end",
            "worst |f| / h", rounding
        ),
        format!(
            "{:<18} {:>10.3}   cells out of its cell -- the RUNAWAY end",
            "worst clamp move", runaway
        ),
        format!("{:<18} {:>10}", "samples/axis", demo.samples),
        String::new(),
        verdict.to_string(),
        String::new(),
        "lambda is the Tikhonov regularizer in the vertex solve, and the whole".into(),
        "sharpness trade in one number. the two failures need two numbers: |f|/h".into(),
        "sees the rounding, and it is BLIND to the runaway, because a flat cell's".into(),
        "unconstrained directions lie *within* the surface -- a vertex with".into(),
        "nothing holding it slides along the plane and stays on it. so the second".into(),
        "number is M-30's: how far the clamp had to drag it back.".into(),
        String::new(),
        "and the two ends live on different fields. M-30 measured the runaway at".into(),
        "3.18 cells on gyroid and 2.17 on fbm_terrain, and *zero* on box_exact --".into(),
        "press 3 or 4 to see it. box_exact is where the rounding shows.".into(),
        String::new(),
        "with the clamp ON you cannot see any of this: A-009 confines each".into(),
        "vertex to its own cell, so the failure never becomes geometry. that is".into(),
        "why it is the default, and why it is off here.".into(),
    ];

    let mesh = meshes.add(builder.into_mesh());
    if query.is_empty() {
        commands.spawn((
            Mesh3d(mesh),
            MeshMaterial3d(surface.0.clone()),
            Transform::default(),
            DemoMesh,
        ));
    } else {
        for mut handle in &mut query {
            handle.0 = mesh.clone();
        }
    }
}
