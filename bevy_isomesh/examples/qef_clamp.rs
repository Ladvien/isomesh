//! E-110 — the cell clamp, and the half of the problem it cannot reach.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example qef_clamp --release
//! ```
//!
//! **Always `--release`.**
//!
//! `C` clamp · `A` algorithm · `I` red · `1`–`4` field · `[` `]` resolution · `W` wireframe.
//!
//! # The catalog says this settles a question. It doesn't — A-009 already did.
//!
//! The examples catalog describes this demo as settling *"whether guaranteed
//! intersection-free extraction is available for free."* That question is O-2,
//! and A-009 answered it. What is left is to make the answer legible **on the
//! mesh**, because the answer has two halves and they live on different fields:
//!
//! - **Free for placement.** Five of the seven reference fields go to *exactly*
//!   zero when the clamp is on, and `box_exact`'s corner measures the same
//!   distance clamped or not — a convex corner's solution is interior to its own
//!   cell, so the constraint never binds where the feature is (M-28).
//! - **Not sufficient overall.** `gyroid` and `fbm_terrain` keep 3.12 and 13.84
//!   pairs per 1,000 triangles, and those two are precisely the fields with
//!   multi-sheet cells (M-4, M-15). What survives is a **connectivity** failure,
//!   not a placement one, and no amount of clamping reaches it (M-29).
//!
//! And the obvious fix makes it worse. Manifold Dual Contouring splits the shared
//! vertex, which is exactly the assumption the clamp's partition argument rests
//! on — `gyroid` goes 3.118 → 5.669 and `fbm_terrain` 13.837 → 15.434 (M-61).
//! Press `A`.
//!
//! # Why every configuration is measured on every re-mesh
//!
//! The headline quantity is a **ratio between two settings**, and a HUD that
//! shows one number at a time asks the reader to hold the other in their head
//! across a keypress. So all four combinations of {clamp off, clamp on} ×
//! {dual contouring, manifold dual contouring} are extracted and counted every
//! re-mesh, and all four are on screen at once. `C` and `A` then change only
//! *which mesh you are looking at and where the red comes from* — no number on
//! screen moves when you press them.
//!
//! `sharp_features` already runs a second reference extraction for the same
//! reason. Here it is four, which is affordable precisely because it happens on
//! re-mesh and never per frame.
//!
//! # λ is deliberately not a knob here
//!
//! `sharp_features` (E-109) owns λ, the Tikhonov regularizer. Turning both it and
//! the clamp in one demo would confound the two, and the numbers would stop being
//! comparable to M-28 — which was measured at the default λ = 0.01. This example
//! leaves the solve alone and moves one thing.
//!
//! # Why the red marks come from the report and not from a second opinion
//!
//! Nothing here decides what counts as an intersection; the marks are
//! [`isomesh::validate::self_intersections`]'s own `pairs`, resolved to
//! positions. Recomputing "which triangles look wrong" for the picture would let
//! the picture and the caption drift with no way to tell which one was lying —
//! the same argument `manifold_check` makes.
//!
//! That has one sharp edge, and it is guarded rather than assumed. `pairs`
//! indexes the validator's **filtered** triangle list, not `indices` — a face is
//! kept only if its indices are in range and distinct. On well-formed output the
//! two coincide, which is exactly what makes indexing `indices` directly
//! dangerous: it would work on every mesh anyone tested it on, and point at
//! arbitrary triangles the first time it didn't. So [`red_triangles`] rebuilds
//! the filtered list with the same predicate and refuses to draw anything unless
//! its length matches the count the crate computed.
//!
//! # What this demo cannot see
//!
//! The counter skips any pair of triangles sharing a vertex index, and dual
//! contouring's quads share vertices with their neighbours across every cell
//! face. So a fold that pinches exactly at a shared vertex draws no red and is
//! counted in nothing (M-83). The number of skipped pairs is on screen next to
//! the number of tested ones, because otherwise this demo would be claiming
//! *"these are the offending triangles"* when it means *"these are the ones the
//! counter can see."*
//!
//! # A prediction, registered before the sweep was run
//!
//! Every self-intersection figure in `FINDINGS.md` (M-28, M-53, M-61) is at a
//! single resolution, 33³. `[` and `]` sweep it, so:
//!
//! > The clamped residue on `gyroid` and `fbm_terrain` **falls as resolution
//! > rises.** M-29 attributes it to multi-sheet cells, and M-15 establishes that
//! > multi-sheet is a *resolution* effect rather than a topological one — "any
//! > feature thinner than one cell forces two sheets through it." Finer cells
//! > should mean fewer such cells and a smaller residue.
//!
//! If it comes out flat or rising, then M-29's attribution or M-15's mechanism is
//! wrong, and that contradiction is the finding rather than a bug in this demo.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::dual_contouring::{CLAMP_EPSILON, Clamp, DualContouring};
use isomesh::fields::{FbmTerrain, ReferenceField, Torus, capped_gyroid, csg_difference};
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::validate::{SelfIntersectionReport, self_intersections};
use isomesh::{RuntimeShape3, Sdf};

/// Ordered so the **default** field shows the effect, which is the mistake M-107
/// exists to stop: `sphere`, `box_exact`, `thin_plate` and `csg_difference` read
/// **zero either way** (M-28), and opening on one of them shows a `C` key that
/// does nothing.
///
/// `gyroid` first — the largest removal in the crate, 71.43 → 3.118 per 1,000
/// triangles, *and* a visible residue, so one field carries both halves of O-2's
/// answer at once. `torus` second — the only reference field that reaches
/// **exactly zero**, which is the "free for placement" half on its own.
/// `fbm_terrain` third — the worst case in both directions, 189.46 → 13.837.
///
/// `csg_difference` last, and it is the **control**: zero either way, present so
/// a reader can confirm the meter reads zero when it should. The HUD says so, so
/// that it does not read as a broken demo.
const FIELDS: [&str; 4] = ["gyroid", "torus", "fbm_terrain", "csg_difference"];

/// M-28's grid, and the default for that reason alone: every rate in this HUD is
/// comparable to `FINDINGS.md` and to `docs/measurements/shootout.csv` without a
/// conversion only at this resolution.
const DEFAULT_SAMPLES: u32 = 33;
const MIN_SAMPLES: u32 = 17;
/// Well below the shootout's 65³. `self_intersections` is a broadphase-accelerated
/// all-pairs test whose own docs say it *"measures a mesh, it does not produce
/// one"*, and this example runs it four times per re-mesh. 49³ on `gyroid` is
/// about 24,000 triangles per run.
const MAX_SAMPLES: u32 = 49;
/// Keeps 33 on the lattice: 17, 21, 25, 29, **33**, 37, 41, 45, 49.
const SAMPLES_STEP: u32 = 4;

/// Frames a capture dwells on one clamp state before flipping it.
const CAPTURE_DWELL: u32 = 8;

/// Offending triangles get their own config group so they can be thick and drawn
/// in front. An unbiased line lying exactly on the surface z-fights with the
/// triangles it outlines and flickers, which is indistinguishable from the defect
/// being intermittent — `manifold_check` learned this first.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct IntersectionGizmos;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Algorithm {
    DualContouring,
    ManifoldDualContouring,
}

impl Algorithm {
    const ALL: [Self; 2] = [Self::DualContouring, Self::ManifoldDualContouring];

    fn name(self) -> &'static str {
        match self {
            Self::DualContouring => "dual contouring",
            Self::ManifoldDualContouring => "manifold dual contouring",
        }
    }

    /// Short enough to keep the four HUD rows aligned.
    fn short(self) -> &'static str {
        match self {
            Self::DualContouring => "dual contouring",
            Self::ManifoldDualContouring => "manifold dc",
        }
    }

    /// `ISOMESH_ALGORITHM=dc|mdc`, so a capture needs no keyboard — the same
    /// reason `ISOMESH_FIELD` and `ISOMESH_VIEW` exist.
    fn from_env() -> Self {
        match std::env::var("ISOMESH_ALGORITHM")
            .unwrap_or_default()
            .as_str()
        {
            "mdc" | "manifold_dual_contouring" => Self::ManifoldDualContouring,
            _ => Self::DualContouring,
        }
    }

    fn toggled(self) -> Self {
        match self {
            Self::DualContouring => Self::ManifoldDualContouring,
            Self::ManifoldDualContouring => Self::DualContouring,
        }
    }

    fn index(self) -> usize {
        match self {
            Self::DualContouring => 0,
            Self::ManifoldDualContouring => 1,
        }
    }
}

/// `ISOMESH_CLAMP=on|off`. A screenshot cannot press `C`, and the before/after
/// pair this example exists to produce has to come off a command line.
fn clamp_from_env() -> bool {
    matches!(
        std::env::var("ISOMESH_CLAMP").unwrap_or_default().as_str(),
        "on" | "1" | "true"
    )
}

#[derive(Resource)]
struct Demo {
    field: usize,
    samples: u32,
    /// Off by default: this example exists to show what the clamp removes, and
    /// with it already on there is nothing to remove.
    clamp: bool,
    algorithm: Algorithm,
    show_red: bool,
}

#[derive(Resource)]
struct SurfaceMaterial(Handle<StandardMaterial>);

/// What one (algorithm, clamp) configuration measured.
#[derive(Clone)]
struct Measured {
    per_1k: f64,
    pairs: u64,
    adjacent_skipped: u64,
    degenerate: u64,
    /// Faces the validator's filter dropped before indexing. Reported rather than
    /// asserted away: it is what tells a reader whether the naive index mapping
    /// would have been equivalent on *this* mesh.
    dropped: usize,
    /// Distinct triangles carrying a pair — not two per pair. A triangle in five
    /// pairs is one red triangle.
    red_triangles: usize,
    extract_ms: f64,
    check_ms: f64,
    /// Set when the counter refused the mesh, and printed instead of a rate.
    refused: Option<String>,
}

/// World-space outlines for the displayed configuration. Resolved once per
/// re-mesh; the draw system only iterates them.
#[derive(Resource, Default)]
struct Overlay {
    red: Vec<[Vec3; 3]>,
    /// Set when the report could not be mapped back to triangles. Nothing is
    /// drawn while this is `Some`, and the HUD prints it instead of a count.
    disabled: Option<String>,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-110 qef clamp".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<IntersectionGizmos>()
        .insert_resource(Demo {
            field: 0,
            samples: common::samples_override().unwrap_or(DEFAULT_SAMPLES),
            clamp: clamp_from_env(),
            algorithm: Algorithm::from_env(),
            show_red: true,
        })
        .init_resource::<Overlay>()
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, remesh, draw_intersections))
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut gizmo_config: ResMut<GizmoConfigStore>,
    mut demo: ResMut<Demo>,
    flags: Res<ViewFlags>,
) {
    // `ISOMESH_FIELD` is the harness contract for choosing a field without a
    // keyboard, and a capture of the wrong field is how E-109 nearly shipped a
    // sweep of a field that could not show the effect.
    demo.field = flags.field.min(FIELDS.len() - 1);

    let (config, _) = gizmo_config.config_mut::<IntersectionGizmos>();
    config.line.width = 5.0;
    config.depth_bias = -0.4;

    commands.insert_resource(SurfaceMaterial(common::surface_material(&mut materials)));
    // Spawned so `G` is not a dead key: on the capped gyroid, knowing where the
    // clip boundary is changes how you read the marks near it.
    commands.spawn(DemoDomain {
        min: Vec3::splat(-1.0),
        max: Vec3::splat(1.0),
    });
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut demo: ResMut<Demo>,
    mut flags: ResMut<ViewFlags>,
) {
    // The parameter this demo sweeps is the clamp itself, so a capture flips it
    // in step with the frames rather than with the clock — the GIF is then the
    // same GIF on any machine.
    if capture.is_active() {
        demo.clamp = (capture.taken / CAPTURE_DWELL) % 2 == 1;
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
    if keys.just_pressed(KeyCode::KeyC) {
        demo.clamp = !demo.clamp;
    }
    if keys.just_pressed(KeyCode::KeyA) {
        demo.algorithm = demo.algorithm.toggled();
    }
    if keys.just_pressed(KeyCode::KeyI) {
        demo.show_red = !demo.show_red;
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        demo.samples = (demo.samples + SAMPLES_STEP).min(MAX_SAMPLES);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        demo.samples = demo.samples.saturating_sub(SAMPLES_STEP).max(MIN_SAMPLES);
    }
    if keys.just_pressed(KeyCode::KeyR) {
        flags.remesh_requested = true;
    }
}

/// Resolve the report's pairs to triangle outlines, or refuse.
///
/// [`SelfIntersectionReport::pairs`] indexes the validator's *filtered* triangle
/// list, so this rebuilds that list with the same predicate rather than indexing
/// `indices` directly. `report.triangles` is *defined* as that list's length, so
/// the equality below is precisely the statement "my copy of the filter agrees
/// with yours" — it cannot pass by accident. If the crate ever changes the
/// predicate, this returns `Err` and the overlay goes dark with a message, rather
/// than painting red on triangles that were never in a pair.
///
/// Returns the outlines and the number of faces the filter dropped, which the HUD
/// reports so a reader can see whether the naive mapping would have been
/// equivalent on this particular mesh.
fn red_triangles(
    positions: &[[f32; 3]],
    indices: &[u32],
    report: &SelfIntersectionReport,
) -> Result<(Vec<[Vec3; 3]>, usize), String> {
    let whole = indices.len() - indices.len() % 3;
    let mut tris: Vec<[u32; 3]> = Vec::with_capacity(whole / 3);
    for tri in indices[..whole].as_chunks::<3>().0 {
        let in_range = tri.iter().all(|&i| (i as usize) < positions.len());
        let distinct = tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2];
        if in_range && distinct {
            tris.push([tri[0], tri[1], tri[2]]);
        }
    }

    if tris.len() as u64 != report.triangles {
        return Err(format!(
            "overlay off: filtered {} faces, report counted {} -- the filter has drifted",
            tris.len(),
            report.triangles
        ));
    }

    // One outline per distinct triangle, not two per pair.
    let mut ids: Vec<u32> = report
        .pairs
        .iter()
        .flat_map(|p| p.iter().copied())
        .collect();
    ids.sort_unstable();
    ids.dedup();

    // A stale report against a newer mesh would index past the end. The list is
    // already sorted, so this is one comparison.
    if ids.last().is_some_and(|&i| i as usize >= tris.len()) {
        return Err("overlay off: a pair indexes past the triangle list".into());
    }

    let at = |i: u32| Vec3::from(positions[i as usize]);
    let red = ids
        .iter()
        .map(|&t| {
            let f = tris[t as usize];
            [at(f[0]), at(f[1]), at(f[2])]
        })
        .collect();
    Ok((red, whole / 3 - tris.len()))
}

/// One extraction, counted.
struct Run {
    builder: MeshBuilder,
    measured: Measured,
    red: Vec<[Vec3; 3]>,
    overlay_error: Option<String>,
}

fn run_one<F>(field: &F, samples: u32, algorithm: Algorithm, clamp: bool) -> Option<Run>
where
    F: Sdf<Scalar = f32> + ReferenceField,
{
    let (min, max) = field.domain();
    let cell_size = (max[0] - min[0]) / (samples - 1) as f32;
    let shape = match RuntimeShape3::new([samples; 3]) {
        Ok(shape) => shape,
        Err(error) => {
            error!("grid {samples}^3 rejected: {error}");
            return None;
        }
    };
    let setting = if clamp { Clamp::ToCell } else { Clamp::None };

    let mut builder = MeshBuilder::new();
    let started = Instant::now();
    // λ is left at the crate default on both paths. That is what makes these
    // numbers M-28's numbers rather than new ones.
    let extracted = match algorithm {
        Algorithm::DualContouring => {
            let mut dc = DualContouring::<f32>::new();
            dc.set_clamp(setting);
            dc.extract(field, &shape, min, cell_size, &mut builder)
        }
        Algorithm::ManifoldDualContouring => {
            let mut mdc = ManifoldDualContouring::<f32>::new();
            mdc.set_clamp(setting);
            mdc.extract(field, &shape, min, cell_size, &mut builder)
        }
    };
    if let Err(error) = extracted {
        error!("extraction failed at {samples}^3: {error}");
        return None;
    }
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;

    let counting = Instant::now();
    let counted = self_intersections(builder.positions(), builder.indices(), f64::from(cell_size));
    let check_ms = counting.elapsed().as_secs_f64() * 1000.0;

    // A refusal is reachable here and nowhere else in the repo. With the clamp
    // off, M-30 measured a vertex 3.18 cells outside its own cell on `gyroid`,
    // and a triangle joining a runaway vertex to a normal one can span more
    // broadphase cells than the guard allows. Report it; do not retry at a
    // coarser spacing, because a rate measured at one spacing is not comparable
    // to a rate measured at another, and comparability is the whole point.
    let report = match counted {
        Ok(report) => report,
        Err(error) => {
            return Some(Run {
                builder,
                measured: Measured {
                    per_1k: 0.0,
                    pairs: 0,
                    adjacent_skipped: 0,
                    degenerate: 0,
                    dropped: 0,
                    red_triangles: 0,
                    extract_ms,
                    check_ms,
                    refused: Some(format!("{error}")),
                },
                red: Vec::new(),
                overlay_error: None,
            });
        }
    };

    let (red, dropped, overlay_error) =
        match red_triangles(builder.positions(), builder.indices(), &report) {
            Ok((red, dropped)) => (red, dropped, None),
            Err(message) => {
                error!("{message}");
                (Vec::new(), 0, Some(message))
            }
        };

    Some(Run {
        measured: Measured {
            per_1k: report.per_thousand_triangles(),
            pairs: report.count(),
            adjacent_skipped: report.adjacent_pairs_skipped,
            degenerate: report.degenerate_triangles,
            dropped,
            red_triangles: red.len(),
            extract_ms,
            check_ms,
            refused: None,
        },
        red,
        overlay_error,
        builder,
    })
}

struct Built {
    mesh: Mesh,
    overlay: Overlay,
    lines: Vec<String>,
    field_name: &'static str,
    domain_min: [f32; 3],
    domain_max: [f32; 3],
    vertices: usize,
    triangles: usize,
    extract_ms: f64,
}

/// Index into the four-configuration table.
fn slot(algorithm: Algorithm, clamp: bool) -> usize {
    algorithm.index() * 2 + usize::from(clamp)
}

/// Dispatch on the field index, then do the work once in [`measure_all`].
fn build(demo: &Demo) -> Option<Built> {
    match demo.field {
        0 => measure_all(&capped_gyroid::<f32>(), demo),
        1 => measure_all(&Torus::<f32>::canonical(), demo),
        2 => measure_all(&FbmTerrain::<f32>::canonical(), demo),
        _ => measure_all(&csg_difference::<f32>(), demo),
    }
}

fn measure_all<F>(field: &F, demo: &Demo) -> Option<Built>
where
    F: Sdf<Scalar = f32> + ReferenceField,
{
    let (min, max) = field.domain();

    // All four, every re-mesh. `C` and `A` then move only which one is shown.
    let mut runs: Vec<Option<Run>> = Vec::with_capacity(4);
    for algorithm in Algorithm::ALL {
        for clamp in [false, true] {
            runs.push(run_one(field, demo.samples, algorithm, clamp));
        }
    }

    let shown = slot(demo.algorithm, demo.clamp);
    // The displayed configuration must exist; the other three may legitimately
    // be missing and are reported as such rather than silently blanked.
    runs.get(shown)?.as_ref()?;

    let at = |algorithm: Algorithm, clamp: bool| -> Option<&Measured> {
        runs.get(slot(algorithm, clamp))
            .and_then(Option::as_ref)
            .map(|r| &r.measured)
    };
    let rate = |algorithm: Algorithm, clamp: bool| -> Option<f64> {
        at(algorithm, clamp)
            .filter(|m| m.refused.is_none())
            .map(|m| m.per_1k)
    };

    let mut lines = vec![
        "self-intersections per 1,000 triangles. all four are extracted and".into(),
        "counted every re-mesh, so no number here moves when you press C or A:".into(),
        String::new(),
    ];
    for algorithm in Algorithm::ALL {
        for clamp in [false, true] {
            let marker = if slot(algorithm, clamp) == shown {
                "<--"
            } else {
                "   "
            };
            let body = match at(algorithm, clamp) {
                None => "          did not extract".to_string(),
                Some(m) => match &m.refused {
                    Some(why) => format!("          refused: {why}"),
                    None => format!("{:>10.3} {:>7} pairs", m.per_1k, m.pairs),
                },
            };
            lines.push(format!(
                "  {marker} {:<16} clamp {:<3}{body}",
                algorithm.short(),
                if clamp { "ON" } else { "off" },
            ));
        }
    }

    let off = rate(demo.algorithm, false);
    let on = rate(demo.algorithm, true);
    lines.push(String::new());
    match (off, on) {
        (Some(off), Some(on)) if off > 0.0 => {
            lines.push(format!(
                "  the clamp removed {:.1}% on this field, {:.1}x lower",
                100.0 * (off - on) / off,
                off / on.max(f64::MIN_POSITIVE),
            ));
        }
        // A percentage of zero is not a number, and both "NaN%" and "0.0%" would
        // be lies about a field that never had any to remove.
        (Some(_), Some(_)) => {
            lines.push("  removed: n/a -- nothing to remove on this field".into())
        }
        _ => lines.push("  removed: n/a -- a configuration did not report".into()),
    }
    if let (Some(dc), Some(mdc)) = (
        rate(Algorithm::DualContouring, true),
        rate(Algorithm::ManifoldDualContouring, true),
    ) && dc > 0.0
    {
        lines.push(format!(
            "  splitting the vertex instead: {:.2}x {} than dual contouring (M-61)",
            mdc / dc,
            if mdc > dc { "WORSE" } else { "better" },
        ));
    }

    let m = at(demo.algorithm, demo.clamp)?;
    let check_ms: f64 = runs.iter().flatten().map(|r| r.measured.check_ms).sum();
    lines.extend([
        String::new(),
        format!(
            "{:>9} triangles in red (distinct, not two per pair)   [I] hides them",
            m.red_triangles
        ),
        format!(
            "{:>9} pairs skipped for sharing a vertex -- a fold pinching exactly",
            m.adjacent_skipped
        ),
        "          at a shared vertex is NOT counted (M-83)".into(),
        format!("{:>9} degenerate triangles", m.degenerate),
        format!(
            "{:>9} faces dropped before indexing -- red is 1:1 with the report",
            m.dropped
        ),
        format!("{check_ms:>9.1} ms counting all four -- per re-mesh, never per frame"),
        String::new(),
        format!(
            "clamp inset {:<12} of the half cell      {} samples/axis",
            format!("(1 - {CLAMP_EPSILON:e})"),
            demo.samples
        ),
        "[C] clamp  [A] algorithm  [I] red  [1-4] field  [ and ] resolution".into(),
        String::new(),
    ]);

    let verdict = match (off, on) {
        (Some(off), Some(on)) if off == 0.0 && on == 0.0 => {
            "NOTHING TO REMOVE HERE -- this field has none either way. that is the".to_string()
        }
        (Some(_), Some(0.0)) => {
            "FREE, HERE -- the clamp took this field to EXACTLY zero (M-28), and the".to_string()
        }
        (Some(off), Some(on)) => format!(
            "NOT FREE OVERALL -- the clamp removed {:.1}% and left {on:.3} per 1k.",
            100.0 * (off - on) / off
        ),
        _ => "INCOMPLETE -- a configuration refused; see the rows above.".to_string(),
    };
    lines.push(verdict);
    lines.extend(match (off, on) {
        (Some(off), Some(on)) if off == 0.0 && on == 0.0 => vec![
            "CONTROL, not the demo: it shows the meter reads zero when it should.".into(),
            "press 1 or 3 for a field where the clamp has work to do.".into(),
        ],
        (Some(_), Some(0.0)) => vec![
            "sharp corner it was protecting measures the same distance either way,".into(),
            "because a convex corner's solution is interior to its own cell.".into(),
        ],
        _ => vec![
            "the clamp fixes PLACEMENT and cannot reach CONNECTIVITY: what is left".into(),
            "is multi-sheet cells, where two sheets of the surface are forced".into(),
            "through one vertex (M-4, M-15, M-29). press A -- splitting that vertex".into(),
            "is the obvious fix and it makes the count WORSE, not better (M-61).".into(),
        ],
    });

    // One CSV line per re-mesh, so a resolution sweep is a shell loop over
    // `ISOMESH_SAMPLES` rather than nine screenshots read by eye. The prediction
    // in this file's header is checked against exactly this output.
    info!(
        "rates,{},{},{},{},{},{}",
        F::NAME,
        demo.samples,
        rate(Algorithm::DualContouring, false).unwrap_or(f64::NAN),
        rate(Algorithm::DualContouring, true).unwrap_or(f64::NAN),
        rate(Algorithm::ManifoldDualContouring, false).unwrap_or(f64::NAN),
        rate(Algorithm::ManifoldDualContouring, true).unwrap_or(f64::NAN),
    );

    // Taken only now. Everything above reads all four runs uniformly through
    // `at`, so the displayed one needs no special case while the HUD is built.
    let run = runs.get_mut(shown).and_then(Option::take)?;
    let overlay = Overlay {
        red: run.red,
        disabled: run.overlay_error,
    };

    Some(Built {
        overlay,
        lines,
        field_name: F::NAME,
        domain_min: min,
        domain_max: max,
        vertices: run.builder.vertex_count(),
        triangles: run.builder.triangle_count(),
        extract_ms: run.measured.extract_ms,
        mesh: run.builder.into_mesh(),
    })
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    demo: Res<Demo>,
    mut stats: ResMut<DemoStats>,
    mut overlay: ResMut<Overlay>,
    mut meshes: ResMut<Assets<Mesh>>,
    material: Res<SurfaceMaterial>,
    mut commands: Commands,
    mut query: Query<&mut Mesh3d, With<DemoMesh>>,
    mut domain: Query<&mut DemoDomain>,
    mut camera: Query<&mut OrbitCamera>,
    mut flags: ResMut<ViewFlags>,
    mut last: Local<Option<(usize, u32, bool, Algorithm)>>,
) {
    // `clamp` and `algorithm` are in the key even though every configuration is
    // measured, because they decide which report the overlay is built from.
    // Leaving them out would leave the previous setting's red marks drawn over
    // the new mesh -- the silent-wrong-marks failure arriving through the back
    // door.
    let key = (demo.field, demo.samples, demo.clamp, demo.algorithm);
    if *last == Some(key) && !flags.remesh_requested {
        return;
    }
    let field_changed = last.map(|(f, ..)| f) != Some(demo.field) || last.is_none();
    *last = Some(key);
    flags.remesh_requested = false;

    let Some(built) = build(&demo) else {
        return;
    };

    for mut d in &mut domain {
        d.min = Vec3::from(built.domain_min);
        d.max = Vec3::from(built.domain_max);
    }
    // Frame the field rather than assuming one size: the capped gyroid's extent
    // is 3.5x the compact fields', and a fixed radius puts the camera inside a
    // tunnel looking at its inner wall.
    if field_changed {
        let extent = built.domain_max[0] - built.domain_min[0];
        for mut orbit in &mut camera {
            orbit.radius = extent * 1.6;
        }
    }

    *overlay = built.overlay;
    stats.title = format!(
        "E-110  qef clamp - {}   clamp {}   field {} ({})   {}^3",
        demo.algorithm.name(),
        if demo.clamp { "on" } else { "OFF" },
        demo.field + 1,
        built.field_name,
        demo.samples,
    );
    stats.vertices = built.vertices;
    stats.triangles = built.triangles;
    stats.extract_ms = built.extract_ms;
    stats.extra = built.lines;

    let handle = meshes.add(built.mesh);
    if query.is_empty() {
        commands.spawn((
            Mesh3d(handle),
            MeshMaterial3d(material.0.clone()),
            Transform::default(),
            DemoMesh,
        ));
    } else {
        for mut mesh in &mut query {
            mesh.0 = handle.clone();
        }
    }
}

fn draw_intersections(
    overlay: Res<Overlay>,
    demo: Res<Demo>,
    mut gizmos: Gizmos<IntersectionGizmos>,
) {
    const RED: Color = Color::srgb(1.0, 0.13, 0.13);

    if !demo.show_red || overlay.disabled.is_some() {
        return;
    }
    for [a, b, c] in &overlay.red {
        gizmos.line(*a, *b, RED);
        gizmos.line(*b, *c, RED);
        gizmos.line(*c, *a, RED);
    }
}
