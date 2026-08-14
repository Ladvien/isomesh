//! E-116 — the MC33 interior decider, and the 12.6% where the classic test is wrong.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example marching_cubes_interior --release
//! ```
//!
//! `Space` pauses the sweep, `[` and `]` step it by hand, `R` restarts, `1` and
//! `2` switch between the counterexample and a configuration where the two tests
//! agree.
//!
//! # What is on screen
//!
//! One cell, and a plane sweeping through it from the bottom face to the top.
//! Marching Cubes 33 decides whether a cell's same-signed corners are joined
//! *through its interior* by asking whether any such cutting plane has a
//! **positive saddle** — Custodio §3.1: *"if there is a plane cutting the cube
//! such that its saddle point is positive … the positive vertices are connected
//! inside the cube."*
//!
//! The white dot is that saddle, and the trail behind it is the path Custodio's
//! Figure 4 draws. It is not a straight line and it does not stay on the face:
//! the saddle's position is a linear function divided by `Δ(t)`, so it is a
//! **hyperbola**, and where `Δ` has a root the saddle runs off to infinity and
//! comes back from the other side. The magenta plane is that pole.
//!
//! # The disagreement this example exists to show
//!
//! The saddle's *value* is `F(t) / Δ(t)` with `F` quadratic and `Δ` linear.
//! Chernyaev's test tracks the sign of `F` alone — and a quadratic has two sign
//! changes to spend where the quotient needs three. So on configuration `1` the
//! numerator is **negative for the whole sweep** while the saddle is
//! **positive past the pole**: the classic test reports the corners separated
//! and they are joined.
//!
//! The HUD prints both verdicts every frame. They differ on `1` and agree on
//! `2`, and the second one is there because an example that only ever showed the
//! failure would misrepresent how often the failure arrives.
//!
//! # What this is not
//!
//! It is not a tunnel being meshed. Knowing a cell has one is a different thing
//! from building it, and MC33's tunnel cases need vertices in the cell interior
//! that this crate's grid-edge-keyed vertex cache has no slot for — that is
//! A-002b. This example shows the decider, which is what A-002c delivered.
//!
//! # The rate, measured
//!
//! Among face pairs opposed in sign — the only structure that can put a pole
//! inside the sweep, and the structure Custodio's own Appendix A counterexample
//! has — the numerator-only test is wrong in **1,966 of 15,625, or 12.6%**
//! (M-165). That is *not* the same as Custodio's "once in 10,000 random 5×5×5
//! fields": theirs is a rate over fields, and it is small because the opposed
//! family is itself rare. The HUD says which number it is showing.

mod common;

use bevy::color::palettes::css;
use bevy::prelude::*;
use common::{CommonPlugin, DemoStats, OrbitCamera};
use isomesh::marching_cubes::interior::{Interior, SweptFaces};

/// One labelled configuration to sweep.
struct Configuration {
    name: &'static str,
    /// `[A, B, C, D]` on the bottom face and the top, in Custodio's cyclic order.
    lo: [f64; 4],
    hi: [f64; 4],
    note: &'static str,
}

/// Both faces ambiguous with their diagonals *opposed* — `A`/`C` positive below
/// and negative above — which is what puts `Δ`'s root inside the sweep.
///
/// `F` is convex here and negative at both ends, and a convex function's maximum
/// on a closed interval is at an endpoint, so `F < 0` throughout and the
/// numerator-only test finds nothing. Derived rather than transcribed: Appendix
/// A's own numbers did not survive the paper's conversion legibly (✗22), and
/// fitting values until they matched would have been the exact failure rule 5
/// exists to prevent. It has the same sign *structure* as Appendix A's, which
/// was not aimed at.
const CONFIGURATIONS: [Configuration; 2] = [
    Configuration {
        name: "opposed — the tests disagree",
        lo: [0.1, -2.0, 10.0, -2.0],
        hi: [-10.0, 2.0, -0.1, 2.0],
        note: "F < 0 for the whole sweep, saddle > 0 past the pole",
    },
    Configuration {
        name: "aligned — the tests agree",
        lo: [1.0, -1.0, 1.0, -1.0],
        hi: [2.0, -3.0, 0.5, -1.0],
        note: "no pole: Δ keeps one sign, so the quotient's sign is F's",
    },
];

#[derive(Resource)]
struct Sweep {
    which: usize,
    t: f64,
    running: bool,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-116 marching cubes interior decider".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Sweep {
            which: 0,
            t: 0.0,
            running: true,
        })
        .add_systems(Startup, aim_camera)
        .add_systems(Update, (controls, advance, draw, report).chain())
        .run();
}

fn aim_camera(mut camera: Query<(&mut Transform, &mut OrbitCamera)>) {
    // The cell is the unit cube centred on the origin, so the camera only needs
    // to be far enough out to see the saddle leave it.
    for (mut transform, mut orbit) in &mut camera {
        orbit.focus = Vec3::ZERO;
        orbit.radius = 3.4;
        *transform =
            Transform::from_translation(Vec3::new(2.2, 1.6, 2.2)).looking_at(Vec3::ZERO, Vec3::Y);
    }
}

fn controls(keys: Res<ButtonInput<KeyCode>>, mut sweep: ResMut<Sweep>) {
    if keys.just_pressed(KeyCode::Space) {
        sweep.running = !sweep.running;
    }
    if keys.just_pressed(KeyCode::KeyR) {
        sweep.t = 0.0;
        sweep.running = true;
    }
    for (key, at) in [(KeyCode::Digit1, 0), (KeyCode::Digit2, 1)] {
        if keys.just_pressed(key) {
            sweep.which = at;
            sweep.t = 0.0;
        }
    }
    // Stepping by hand pauses, because scrubbing while it runs fights itself.
    let step = 0.004;
    if keys.pressed(KeyCode::BracketRight) {
        sweep.running = false;
        sweep.t = (sweep.t + step).min(1.0);
    }
    if keys.pressed(KeyCode::BracketLeft) {
        sweep.running = false;
        sweep.t = (sweep.t - step).max(0.0);
    }
}

fn advance(time: Res<Time>, mut sweep: ResMut<Sweep>) {
    if !sweep.running {
        return;
    }
    sweep.t += f64::from(time.delta_secs()) * 0.12;
    if sweep.t > 1.0 {
        sweep.t = 0.0;
    }
}

/// The cell, the sweeping plane, the saddle and its trail.
fn draw(sweep: Res<Sweep>, mut gizmos: Gizmos) {
    let Some(config) = CONFIGURATIONS.get(sweep.which) else {
        return;
    };
    let Ok(faces) = SweptFaces::new(config.lo, config.hi) else {
        return;
    };

    // The unit cube, centred. `A` at (0,0), `B` at (1,0), `C` at (1,1),
    // `D` at (0,1) in the plane's own coordinates -- the layout
    // `saddle_position` documents.
    let at = |u: f64, v: f64, t: f64| {
        Vec3::new((u - 0.5) as f32, (t - 0.5) as f32 * 1.2, (v - 0.5) as f32)
    };
    let corner_uv = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];

    for t in [0.0, 1.0] {
        for k in 0..4 {
            let (u0, v0) = corner_uv[k];
            let (u1, v1) = corner_uv[(k + 1) % 4];
            gizmos.line(at(u0, v0, t), at(u1, v1, t), css::DIM_GRAY);
        }
    }
    for (u, v) in corner_uv {
        gizmos.line(at(u, v, 0.0), at(u, v, 1.0), css::DIM_GRAY);
    }

    // Corner signs. Red is positive and blue negative, matching the paper's own
    // "black dots are positive" only in spirit -- a screenshot needs hue.
    for (t, values) in [(0.0, config.lo), (1.0, config.hi)] {
        for (k, (u, v)) in corner_uv.iter().enumerate() {
            let colour = if values[k] > 0.0 {
                css::ORANGE_RED
            } else {
                css::DODGER_BLUE
            };
            gizmos.sphere(Isometry3d::from_translation(at(*u, *v, t)), 0.045, colour);
        }
    }

    // The pole: the one height where the saddle does not exist.
    if let Some(pole) = faces.pole() {
        for k in 0..4 {
            let (u0, v0) = corner_uv[k];
            let (u1, v1) = corner_uv[(k + 1) % 4];
            gizmos.line(at(u0, v0, pole), at(u1, v1, pole), css::MAGENTA);
        }
    }

    // The cutting plane at the current height.
    for k in 0..4 {
        let (u0, v0) = corner_uv[k];
        let (u1, v1) = corner_uv[(k + 1) % 4];
        gizmos.line(at(u0, v0, sweep.t), at(u1, v1, sweep.t), css::WHITE);
    }

    // The saddle's trajectory. Drawn as points rather than a polyline because it
    // is a hyperbola with a pole in the middle -- joining the last sample before
    // the pole to the first after it would draw a chord across a discontinuity
    // and hide the very thing this example is about.
    let steps = 240;
    for step in 0..=steps {
        let t = f64::from(step) / f64::from(steps);
        if faces.pole() == Some(t) {
            continue;
        }
        let [u, v] = faces.saddle_position(t);
        if !(-0.6..=1.6).contains(&u) || !(-0.6..=1.6).contains(&v) {
            continue;
        }
        let positive = faces.saddle(t) > 0.0;
        let colour = if positive { css::GOLD } else { css::STEEL_BLUE };
        gizmos.sphere(Isometry3d::from_translation(at(u, v, t)), 0.008, colour);
    }

    // Where it is now.
    let [u, v] = faces.saddle_position(sweep.t);
    if (-0.6..=1.6).contains(&u) && (-0.6..=1.6).contains(&v) {
        let here = at(u, v, sweep.t);
        gizmos.sphere(Isometry3d::from_translation(here), 0.05, css::WHITE);
        // A tick toward the plane's centre, so the dot reads as being *on* the
        // plane rather than floating near it.
        gizmos.line(here, at(0.5, 0.5, sweep.t), css::DIM_GRAY);
    }
}

/// Both verdicts, every frame, with the numbers behind them.
fn report(sweep: Res<Sweep>, mut stats: ResMut<DemoStats>) {
    let Some(config) = CONFIGURATIONS.get(sweep.which) else {
        return;
    };
    let Ok(faces) = SweptFaces::new(config.lo, config.hi) else {
        return;
    };

    // Chernyaev's test tracks the sign of the numerator alone. Sampled here
    // rather than solved, because this is a demonstration of what it sees and a
    // dense sample is the most generous reading of it -- if it still misses the
    // positive saddle, it is not missing it for want of resolution.
    let steps = 2000;
    let numerator_positive =
        (0..=steps).any(|k| faces.numerator(f64::from(k) / f64::from(steps)) > 0.0);
    let classic = if numerator_positive {
        Interior::Joined
    } else {
        Interior::Separated
    };
    let corrected = faces.test();

    let saddle = faces.saddle(sweep.t);
    let pole = faces
        .pole()
        .map_or_else(|| "none in (0, 1)".to_string(), |t| format!("t = {t:.4}"));

    stats.title = format!("E-116  interior decider — {}", config.name);
    stats.extra = vec![
        format!(
            "[{}]  t {:.4}",
            if sweep.running { "running" } else { "paused" },
            sweep.t
        ),
        format!(
            "F(t) {:>12.5}   D(t) {:>12.5}   saddle {:>12.5}",
            faces.numerator(sweep.t),
            faces.denominator(sweep.t),
            if saddle.is_finite() { saddle } else { f64::NAN },
        ),
        format!("pole: {pole}"),
        String::new(),
        format!("Chernyaev  (numerator only): {classic:?}"),
        format!("Custodio   (corrected)     : {corrected:?}"),
        if classic == corrected {
            "the two agree here".to_string()
        } else {
            "*** they DISAGREE — the classic test is wrong ***".to_string()
        },
        String::new(),
        config.note.to_string(),
        "12.6% of opposed face pairs are decided wrongly by the numerator (M-165)".to_string(),
        "space pause | [ ] step | R restart | 1 2 configuration".to_string(),
    ];
}
