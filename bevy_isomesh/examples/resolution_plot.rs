//! M-002 — the fit, drawn, and the point where the model stops describing it.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example resolution_plot --release
//! ```
//!
//! **Always `--release`.**
//!
//! `R` re-run the sweep · `1`–`3` which curve is highlighted.
//!
//! # What this is for, given the numbers already exist
//!
//! `cargo bench --bench resolution_sweep` already fits `t = a + b·n³` and writes
//! the coefficients out. This is the same fit with the residuals visible, and the
//! residuals are the entire point: the two-term model **describes Marching Cubes
//! and does not describe Surface Nets**, and no pair of coefficients can tell you
//! that. A plot can.
//!
//! # The ticket asked for the wrong thing, and its own bench said so
//!
//! M-002 wanted the fixed cost `a` printed, expecting it to be large — the
//! premise being that a big constant overhead per extraction is what a chunked
//! game pays. It is not. `marching_cubes` fits **a = 0.5118 ms, 0.64% of the
//! largest run** (M-62): there is no meaningful fixed cost to find.
//!
//! Over the range a *live* plot can afford — 17³ to 89³, because 256³ is seven
//! seconds for one point — every fitted `a` is under **0.18 ms** in absolute
//! value. There is no fixed cost here worth the name, which is M-62's conclusion
//! reached from this machine instead of the committed CSV.
//!
//! The **sign** is mostly stable and it agrees with M-21. Across eight runs:
//! Surface Nets negative in **7 of 8**, Marching Cubes positive in 7 of 8, Dual
//! Contouring positive in 8 of 8. So the *direction* reproduces and the
//! *magnitude* does not — M-21's `−3.13 ms` comes from the full sweep to 256³,
//! and the same coefficient here is about `−0.07`, forty times smaller. `cargo
//! bench --bench resolution_sweep` is the authority for the numbers; this is the
//! shape, and the shape is what shows the model fitting one curve better than
//! another.
//!
//! # Measured live, not read from the CSV
//!
//! One resolution per frame, timed on this machine, plotted as it arrives — so
//! the curve is this hardware's rather than the committed one's, and the shape
//! can be compared against `docs/measurements/resolution_sweep.csv` by eye.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use common::{Capture, CommonPlugin, DemoStats, OrbitCamera, ViewFlags};
use isomesh::dual_contouring::DualContouring;
use isomesh::fields::ReferenceField;
use isomesh::fields::Sphere;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::surface_nets::SurfaceNets;
use isomesh::{MeshBuffer, RuntimeShape3};

/// Samples per axis at each step of the sweep.
///
/// Starts well above the smallest useful grid: M-19 records that Marching Cubes'
/// fitted `a` is 0.61% of the largest run and **543% of the smallest**, so a
/// sweep that begins too low measures dispatch rather than extraction and makes
/// the intercept look meaningful.
const STEPS: [u32; 10] = [17, 25, 33, 41, 49, 57, 65, 73, 81, 89];

/// Plot rectangle, in world units. The camera looks straight at it.
const PLOT_W: f32 = 9.0;
const PLOT_H: f32 = 6.0;

#[derive(Default, Reflect, GizmoConfigGroup)]
struct PlotGizmos;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Curve {
    MarchingCubes,
    SurfaceNets,
    DualContouring,
}

impl Curve {
    const ALL: [Self; 3] = [Self::MarchingCubes, Self::SurfaceNets, Self::DualContouring];

    fn name(self) -> &'static str {
        match self {
            Self::MarchingCubes => "marching cubes",
            Self::SurfaceNets => "surface nets",
            Self::DualContouring => "dual contouring",
        }
    }

    fn colour(self) -> Color {
        match self {
            Self::MarchingCubes => Color::srgb(0.45, 0.80, 1.00),
            Self::SurfaceNets => Color::srgb(1.00, 0.70, 0.25),
            Self::DualContouring => Color::srgb(0.55, 0.95, 0.55),
        }
    }
}

/// `t = a + b·x` by least squares, with the fit quality.
///
/// `x` is `n³` rather than `n`, which is the model the ticket names. Returned as
/// it comes out: **a negative `a` is not clamped**, because that value is the
/// finding (M-21) and hiding it behind a `max(0.0)` would turn "the model does
/// not apply" into "the fixed cost is zero", which is a different and false
/// statement.
#[derive(Clone, Copy, Default)]
struct Fit {
    a: f64,
    b: f64,
    r2: f64,
}

fn fit(points: &[(f64, f64)]) -> Fit {
    if points.len() < 2 {
        return Fit::default();
    }
    let n = points.len() as f64;
    let mx = points.iter().map(|p| p.0).sum::<f64>() / n;
    let my = points.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx: f64 = points.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum();
    let sxy: f64 = points.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    if sxx <= 0.0 {
        return Fit::default();
    }
    let b = sxy / sxx;
    let a = my - b * mx;
    let ss_tot: f64 = points.iter().map(|p| (p.1 - my) * (p.1 - my)).sum();
    let ss_res: f64 = points
        .iter()
        .map(|p| {
            let e = p.1 - (a + b * p.0);
            e * e
        })
        .sum();
    Fit {
        a,
        b,
        r2: if ss_tot > 0.0 {
            1.0 - ss_res / ss_tot
        } else {
            0.0
        },
    }
}

#[derive(Resource)]
struct Sweep {
    /// Next step to measure, or `STEPS.len()` when finished.
    next: usize,
    /// `(n³, ms)` per curve.
    points: [Vec<(f64, f64)>; 3],
    fits: [Fit; 3],
    highlight: usize,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — M-002 resolution plot".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<PlotGizmos>()
        .insert_resource(Sweep {
            next: 0,
            points: [Vec::new(), Vec::new(), Vec::new()],
            fits: [Fit::default(); 3],
            highlight: 0,
        })
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, measure, draw, report, hud).chain())
        .run();
}

fn setup(mut gizmo_config: ResMut<GizmoConfigStore>, mut camera: Query<&mut OrbitCamera>) {
    let (config, _) = gizmo_config.config_mut::<PlotGizmos>();
    config.line.width = 2.5;
    for mut orbit in &mut camera {
        // Straight on. This is a chart, and any pitch or yaw is a chart drawn in
        // perspective, which is a chart you cannot read values off.
        //
        // `yaw = PI/2`, not 0: the harness builds its direction as
        // `(cos yaw cos pitch, sin pitch, sin yaw cos pitch)`, so `yaw = 0` looks
        // along **+X** and renders an XY chart edge-on -- which is exactly what
        // the first capture showed, three coloured vertical lines and no plot.
        // Focus raised so the chart sits *below* the HUD rather than behind it.
        // This example's HUD is twenty lines and the axes' origin -- where a
        // negative intercept shows -- is the one part that must not be covered.
        orbit.focus = Vec3::new(PLOT_W * 0.5, PLOT_H * 0.5 + 3.6, 0.0);
        orbit.yaw = std::f32::consts::FRAC_PI_2;
        orbit.pitch = 0.0;
        orbit.radius = 20.0;
    }
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut sweep: ResMut<Sweep>,
    mut flags: ResMut<ViewFlags>,
) {
    flags.grid = false;
    if capture.is_active() {
        return;
    }
    if keys.just_pressed(KeyCode::KeyR) {
        sweep.next = 0;
        sweep.points = [Vec::new(), Vec::new(), Vec::new()];
        sweep.fits = [Fit::default(); 3];
    }
    for (key, i) in [
        (KeyCode::Digit1, 0),
        (KeyCode::Digit2, 1),
        (KeyCode::Digit3, 2),
    ] {
        if keys.just_pressed(key) {
            sweep.highlight = i;
        }
    }
}

/// One resolution per frame, all three curves.
///
/// Deliberately one step per frame rather than the whole sweep at once: the
/// largest grids take tens of milliseconds each and doing them in a single frame
/// would freeze the window for a second with nothing on screen, which is a worse
/// demo than watching it fill in.
fn measure(mut sweep: ResMut<Sweep>) {
    if sweep.next >= STEPS.len() {
        return;
    }
    let n = STEPS[sweep.next];
    sweep.next += 1;

    let field = Sphere::<f32>::canonical();
    let (lo, hi) = field.domain();
    let cell = (hi[0] - lo[0]) / (n - 1) as f32;
    let Ok(shape) = RuntimeShape3::new([n; 3]) else {
        return;
    };
    let x = f64::from(n) * f64::from(n) * f64::from(n);

    for (i, curve) in Curve::ALL.iter().copied().enumerate() {
        // Median of three, not one shot. A single timing per resolution is
        // noisy enough that the fitted intercept swung 0.25 / 0.20 / 0.18 ms
        // across consecutive runs of the first version -- the *sign* was stable,
        // which is the finding, but the magnitude was not, and a plot that
        // reports a different number every time it is opened invites the reader
        // to believe whichever one they saw. `cargo bench --bench
        // resolution_sweep` remains the authority; this is the shape.
        let mut runs = [0.0f64; 3];
        let mut failed = false;
        for slot in &mut runs {
            let mut out = MeshBuffer::<f32>::new();
            let started = Instant::now();
            let ok = match curve {
                Curve::MarchingCubes => {
                    MarchingCubes::<f32>::new().extract(&field, &shape, lo, cell, &mut out)
                }
                Curve::SurfaceNets => {
                    SurfaceNets::<f32>::new().extract(&field, &shape, lo, cell, &mut out)
                }
                Curve::DualContouring => {
                    DualContouring::<f32>::new().extract(&field, &shape, lo, cell, &mut out)
                }
            };
            if ok.is_err() {
                failed = true;
                break;
            }
            *slot = started.elapsed().as_secs_f64() * 1000.0;
        }
        if failed {
            continue;
        }
        runs.sort_by(f64::total_cmp);
        let ms = runs[1];
        sweep.points[i].push((x, ms));
        sweep.fits[i] = fit(&sweep.points[i]);
    }
}

fn draw(sweep: Res<Sweep>, mut gizmos: Gizmos<PlotGizmos>) {
    const AXIS: Color = Color::srgb(0.55, 0.58, 0.65);
    const ZERO: Color = Color::srgb(0.85, 0.35, 0.35);

    let max_x = STEPS
        .last()
        .map_or(1.0, |n| f64::from(*n) * f64::from(*n) * f64::from(*n));
    let max_y = sweep
        .points
        .iter()
        .flatten()
        .map(|p| p.1)
        .fold(1.0f64, f64::max)
        * 1.15;

    // The axes. `y = 0` is drawn separately and in red, because a fitted line
    // crossing below it is the whole finding and it has to be visible where.
    gizmos.line(Vec3::ZERO, Vec3::new(PLOT_W, 0.0, 0.0), ZERO);
    gizmos.line(Vec3::ZERO, Vec3::new(0.0, PLOT_H, 0.0), AXIS);

    let to_world = |x: f64, y: f64| {
        Vec3::new(
            (x / max_x) as f32 * PLOT_W,
            (y / max_y) as f32 * PLOT_H,
            0.0,
        )
    };

    for (i, curve) in Curve::ALL.iter().copied().enumerate() {
        let points = &sweep.points[i];
        if points.is_empty() {
            continue;
        }
        let colour = curve.colour();
        let dim = colour.with_alpha(if i == sweep.highlight { 1.0 } else { 0.35 });

        // The measurements, joined.
        for pair in points.windows(2) {
            gizmos.line(
                to_world(pair[0].0, pair[0].1),
                to_world(pair[1].0, pair[1].1),
                dim,
            );
        }
        for p in points {
            let at = to_world(p.0, p.1);
            gizmos.circle(Isometry3d::from_translation(at), 0.055, dim);
        }

        // The fit, drawn from x = 0 so its intercept is on screen. A line whose
        // `a` is negative starts *below* the red axis, which is the picture the
        // coefficient alone cannot give.
        let f = sweep.fits[i];
        if points.len() >= 2 {
            let y0 = f.a;
            let y1 = f.a + f.b * max_x;
            gizmos.line(to_world(0.0, y0), to_world(max_x, y1), dim);
        }
    }
}

/// One CSV row when the sweep finishes, so the intercept's stability across runs
/// can be counted rather than asserted.
fn report(sweep: Res<Sweep>, mut said: Local<bool>) {
    if *said || sweep.next < STEPS.len() {
        return;
    }
    *said = true;
    info!(
        "fit,{:.4},{:.4},{:.4},{:.5},{:.5},{:.5}",
        sweep.fits[0].a,
        sweep.fits[1].a,
        sweep.fits[2].a,
        sweep.fits[0].r2,
        sweep.fits[1].r2,
        sweep.fits[2].r2,
    );
}

fn hud(sweep: Res<Sweep>, mut stats: ResMut<DemoStats>) {
    let done = sweep.next >= STEPS.len();
    stats.title = format!(
        "M-002  resolution sweep   {} of {} steps   sphere, f32",
        sweep.next.min(STEPS.len()),
        STEPS.len()
    );
    stats.vertices = 0;
    stats.triangles = 0;

    let mut lines = vec![
        "t = a + b n^3, fitted live on this machine. a is the FIXED cost the".into(),
        "ticket asked for, and the plot is here because two coefficients cannot".into(),
        "show you that the model fits one curve and not another.".into(),
        String::new(),
        format!(
            "{:<20} {:>11}  {:>11}  {:>8}",
            "", "a (ms)", "b (ms/n^3)", "r^2"
        ),
    ];
    for (i, curve) in Curve::ALL.iter().copied().enumerate() {
        let f = sweep.fits[i];
        let mark = if i == sweep.highlight { "<--" } else { "   " };
        lines.push(format!(
            "{mark} {:<16} {:>11.4} {:>12.3e} {:>8.5}",
            curve.name(),
            f.a,
            f.b,
            f.r2
        ));
    }

    let largest = sweep
        .points
        .first()
        .and_then(|p| p.last())
        .map_or(0.0, |p| p.1);
    lines.extend([
        String::new(),
        if done && largest > 0.0 {
            format!(
                "marching cubes' a is {:.2}% of its largest run -- there is no",
                100.0 * sweep.fits[0].a / largest
            )
        } else {
            "sweeping...".to_string()
        },
        "meaningful fixed cost to find, which is the opposite of what M-002".into(),
        "expected (M-62). and M-19's rule is why the range is on screen: a".into(),
        "coefficient means nothing until it is compared against the data's own".into(),
        "ends -- the same 0.61% is 543% of the SMALLEST run.".into(),
        String::new(),
        "over THIS range every a is under 0.18 ms in absolute value -- there is".into(),
        "no fixed cost here worth the name. the SIGN is mostly stable and it".into(),
        "agrees with M-21: across eight runs surface nets came out negative in".into(),
        "7, marching cubes positive in 7, dual contouring positive in 8.".into(),
        String::new(),
        "so the direction reproduces and the MAGNITUDE does not: M-21's -3.13".into(),
        "ms for surface nets comes from the committed sweep out to 256^3, and".into(),
        "the same coefficient here is about -0.07, forty times smaller. this".into(),
        "stops at 89^3 because a live plot cannot spend seven seconds on one".into(),
        "point. cargo bench --bench resolution_sweep is the authority for the".into(),
        "numbers; this is the shape, and the shape is what shows the model".into(),
        "fitting one curve better than another.".into(),
        String::new(),
        "[1] [2] [3] highlight a curve   [R] re-run the sweep".into(),
    ]);
    stats.extra = lines;
}
