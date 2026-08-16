//! E-112 — two failures, two laws, and neither one is where the ticket said.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example precision_f32_vs_f64 --release
//! ```
//!
//! **Always `--release`.**
//!
//! `-` `=` offset · `[` `]` resolution · `1`–`3` field · `W` wireframe.
//!
//! # The premise was measured before this was written, and it was wrong
//!
//! E-112 asked for *"the same field at ~1e6 offsets; f32 cracks, f64 doesn't"*.
//! **At 1e6 `f32` does not crack.** It returns the same 1,160 vertices and 2,316
//! triangles as `f64`, `χ = 2`, zero boundary edges — a topologically perfect
//! sphere. What moves is accuracy: 1.38 cells against `f64`'s 0.0362.
//!
//! The crack is an order of magnitude further out, and the reason it is *exactly*
//! there is the point of this demo.
//!
//! # Two knobs, and they are not the same knob
//!
//! **Accuracy is relative.** The worst distance from a vertex to the true surface
//! grows with `ulp(offset) / h` — one representable step measured in cells. From
//! this example at power-of-two offsets, worst `|f|/h`:
//!
//! | offset | `ulp/h` @33³ | worst @33³ | `ulp/h` @65³ | worst @65³ |
//! |---|---:|---:|---:|---:|
//! | `2²⁰` | 1 | 1.3808 | 2 | 3.1938 |
//! | `2²¹` | 2 | 3.1010 | 4 | 6.2020 |
//! | `2²²` | 4 | 4.0000 | 8 | 8.0000 |
//! | `2²³` | 8 | 5.8564 | 16 | 11.7128 |
//!
//! The clean law is in the *rows*, not the columns: **at the same offset, halving
//! `h` exactly doubles the error measured in cells** — `3.1010 → 6.2020`,
//! `4.0000 → 8.0000`, `5.8564 → 11.7128`. The absolute error is set by
//! `ulp(offset)` alone; expressing it in cells is what makes a finer grid look
//! worse. Press `[` and watch the left column degrade while the offset stays put.
//!
//! There is deliberately no fitted constant here. An earlier draft claimed
//! `≈1.4 · ulp/h`, which held at `ulp/h = 1` and `2` and then broke: the constant
//! also depends on how the offset happens to sit on the `f32` lattice, so a
//! non-power-of-two offset with the *same* `ulp/h = 4` measures `5.8564` rather
//! than `4.0000`. The table is measured; the formula was tidier than the truth.
//!
//! **Topology is absolute.** The mesh tears when `ulp(f32) ≥ 1`, which is
//! `offset ≥ 2²³ = 8,388,608`, and that threshold **does not move with the cell
//! size**. `χ` drops 2 → 1, the vertex count collapses, and boundary edges — real
//! holes — appear.
//!
//! The two are independent, and one fixture proves it: at 65³ and `ulp/h = 8`
//! (offset `2²²`) the mesh is topologically **clean**, `χ = 2` with zero boundary
//! edges, while at 33³ and the *same* `ulp/h = 8` (offset `2²³`) it is torn,
//! `χ = 1` with 42 holes. So the crack cannot be read off `ulp/h`, and the
//! blurring cannot be read off the offset alone. A demo showing one number would
//! have shown half the finding.
//!
//! # Why there is no condition number here, though the ticket asked for one
//!
//! The QEF matrix's condition number describes a *cell's normals* — how
//! ill-posed that vertex is to place. Translating the whole field leaves the
//! normals alone, so κ is very nearly unchanged by the thing this demo varies. It
//! would have sat in the HUD looking relevant and explaining nothing, which is
//! precisely what `|f|/h` did to E-109. The instruments that do predict the two
//! failures are `ulp/h` and `ulp ≥ 1`, so those are the ones on screen.
//!
//! # Two suspects that were ruled out before the demo was built
//!
//! Both were plausible and both were wrong, which is why they are recorded rather
//! than left for someone to re-derive (M-112).
//!
//! *The validator, not the mesh.* `validate.rs` quantises with `as_f32() as iN`
//! and T-008 anchors that lattice to the mesh's own minimum. Re-validating the
//! same `f32` vertices after moving them back to the origin **in `f64`** gives
//! bit-identical reports — `χ = 1` and 54 boundary edges either way. The holes are
//! in the mesh.
//!
//! *The fixture's gradient.* [`Sdf::gradient`]'s default is a central difference
//! whose step scales with `|p|`, and its own doc warns that a field whose
//! characteristic length is far from 1 should override it. Re-run with an
//! analytic sphere normal, every number is bit-identical. The gradient is not
//! involved.
//!
//! # What the `f64` side is also demonstrating
//!
//! It keeps `f64` through the whole extraction and narrows to `f32` only at the
//! end, *after* subtracting the offset — which is what `Real::as_f32`'s doc tells
//! a CAD consumer to do. The `f32` side subtracts its offset in `f32`, because
//! that is what a consumer who chose `f32` actually gets. Neither side is being
//! handicapped to make a point.

mod common;

use std::time::Instant;

use bevy::asset::RenderAssetUsages;
use bevy::mesh::{Indices, PrimitiveTopology};
use bevy::prelude::*;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::dual_contouring::DualContouring;
use isomesh::fields::{BoxExact, FieldBound, ReferenceField, Sphere, Torus};
use isomesh::validate::{ValidateConfig, validate_features};
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

const FIELDS: [&str; 3] = ["sphere", "torus", "box_exact"];

/// Offsets are swept in **powers of two**, not decades, so the sweep lands
/// exactly on `2²³` — the offset at which `ulp(f32)` reaches 1 and the mesh
/// tears. A decade sweep steps straight over the only interesting point.
const MIN_LOG2_OFFSET: i32 = 0;
const MAX_LOG2_OFFSET: i32 = 26;
/// `2²³ = 8,388,608`, the first offset that tears at any cell size. The demo
/// opens on the failure rather than on the healthy case.
const DEFAULT_LOG2_OFFSET: i32 = 23;

/// `ISOMESH_OFFSET=<k>` selects `2^k`, the same way `ISOMESH_SAMPLES` selects a
/// grid. The offset *is* the parameter this demo varies, so without it a capture
/// could only ever photograph the default — and the finding is the difference
/// between two offsets, which no single still can carry.
fn log2_offset_from_env() -> i32 {
    std::env::var("ISOMESH_OFFSET")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(DEFAULT_LOG2_OFFSET)
        .clamp(MIN_LOG2_OFFSET, MAX_LOG2_OFFSET)
}

const DEFAULT_SAMPLES: u32 = 33;
const MIN_SAMPLES: u32 = 17;
const MAX_SAMPLES: u32 = 65;
const SAMPLES_STEP: u32 = 8;

/// Which precision an entity is showing.
#[derive(Component, Clone, Copy, PartialEq, Eq)]
enum Side {
    F32,
    F64,
}

#[derive(Resource)]
struct Demo {
    field: usize,
    samples: u32,
    log2_offset: i32,
}

impl Demo {
    fn offset(&self) -> f64 {
        (2.0f64).powi(self.log2_offset)
    }
}

#[derive(Resource)]
struct Materials {
    f32_side: Handle<StandardMaterial>,
    f64_side: Handle<StandardMaterial>,
}

/// What one precision produced at this offset.
struct Measured {
    vertices: usize,
    triangles: usize,
    euler: i64,
    boundary_edges: u64,
    non_manifold_edges: u64,
    /// Worst `|f(p)| / h` over the vertices, **evaluated in `f64`** whichever
    /// precision produced them. Asking the `f32` field how far off its own
    /// vertices are would ask the damaged arithmetic to grade itself.
    worst_off_surface: f64,
    extract_ms: f64,
    mesh: Mesh,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-112 precision f32 vs f64".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Demo {
            field: 0,
            samples: common::samples_override().unwrap_or(DEFAULT_SAMPLES),
            log2_offset: log2_offset_from_env(),
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
    // keyboard, and it is what makes a capture of a particular field reproducible.
    demo.field = flags.field.min(FIELDS.len() - 1);

    for mut orbit in &mut camera {
        orbit.yaw = 0.6;
        orbit.pitch = 0.3;
        // Framing a *pair* takes more than framing one mesh: the two are offset
        // along x by 0.6 of a domain each way, so the thing to fit is roughly
        // 2.2 domains of width, not one.
        orbit.radius = 14.0;
        // And the pair is pushed below the HUD rather than centred behind it.
        // This example's HUD is thirty lines of two-column table, so it covers
        // the upper left -- which is exactly where the f32 mesh sits, and the
        // f32 mesh is the subject. Centring the pair photographs the argument
        // with its evidence hidden.
        orbit.focus = Vec3::new(0.0, 3.0, 0.0);
    }
    let mut side = |r: f32, g: f32, b: f32| {
        materials.add(StandardMaterial {
            base_color: Color::srgb(r, g, b),
            perceptual_roughness: 0.45,
            // The torn `f32` mesh is full of holes, and a single-sided material
            // renders their far walls invisible -- which reads as a smaller
            // defect than it is.
            double_sided: true,
            cull_mode: None,
            ..default()
        })
    };
    commands.insert_resource(Materials {
        f32_side: side(0.86, 0.55, 0.42),
        f64_side: side(0.62, 0.72, 0.84),
    });
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut demo: ResMut<Demo>,
    mut flags: ResMut<ViewFlags>,
) {
    // A capture walks the offset in step with the frames rather than the clock,
    // and ping-pongs so the sequence loops.
    if capture.is_active() {
        const LOW: i32 = 18;
        let steps = (MAX_LOG2_OFFSET - LOW + 1) as u32;
        let phase = capture.taken % (steps * 2);
        let step = if phase < steps {
            phase
        } else {
            steps * 2 - phase - 1
        };
        demo.log2_offset = LOW + step as i32;
        return;
    }
    for (key, index) in [
        (KeyCode::Digit1, 0),
        (KeyCode::Digit2, 1),
        (KeyCode::Digit3, 2),
    ] {
        if keys.just_pressed(key) {
            demo.field = index;
        }
    }
    if keys.just_pressed(KeyCode::Equal) {
        demo.log2_offset = (demo.log2_offset + 1).min(MAX_LOG2_OFFSET);
    }
    if keys.just_pressed(KeyCode::Minus) {
        demo.log2_offset = (demo.log2_offset - 1).max(MIN_LOG2_OFFSET);
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

/// Any field, translated to `by` on every axis.
///
/// Forwards [`Sdf::gradient`] as well as [`Sdf::sample`]. It must: forwarding
/// only `sample` would silently substitute the central-difference default for an
/// analytic gradient, which is the trap the `&S` impl in the core crate already
/// documents — and here it would have quietly changed what the demo measures.
struct Offset<F: Sdf> {
    inner: F,
    by: F::Scalar,
}

impl<F: Sdf> Sdf for Offset<F> {
    type Scalar = F::Scalar;

    fn sample(&self, p: [Self::Scalar; 3]) -> Self::Scalar {
        self.inner
            .sample([p[0] - self.by, p[1] - self.by, p[2] - self.by])
    }

    fn gradient(&self, p: [Self::Scalar; 3]) -> [Self::Scalar; 3] {
        self.inner
            .gradient([p[0] - self.by, p[1] - self.by, p[2] - self.by])
    }
}

/// A translated field is still the same field: same topology, same Euler
/// characteristic, same distance property. Only the box it lives in moves, and
/// **that translation is the one piece of arithmetic this demo is about** — at
/// large `by` in `f32` the shifted domain corners are themselves inexact, which
/// is exactly the error under study rather than a flaw in the wrapper.
impl<F: Sdf + ReferenceField> ReferenceField for Offset<F> {
    const NAME: &'static str = F::NAME;

    fn domain(&self) -> ([Self::Scalar; 3], [Self::Scalar; 3]) {
        let (min, max) = self.inner.domain();
        (
            [min[0] + self.by, min[1] + self.by, min[2] + self.by],
            [max[0] + self.by, max[1] + self.by, max[2] + self.by],
        )
    }

    fn closed_in_domain(&self) -> bool {
        self.inner.closed_in_domain()
    }

    fn expected_euler(&self) -> Option<i64> {
        self.inner.expected_euler()
    }

    /// Forwarded, because translation preserves it exactly.
    ///
    /// `f(p − by)` has the same Lipschitz constant as `f`, and if `|f|` was the
    /// distance to the zero set then so is `|f(p − by)|` to the translated one.
    /// **This replaced a forward of `is_exact_distance`, which F-001 removed
    /// from the trait** — a `bool` could not say *how* inexact, and this example
    /// went on calling it for 58 commits because the root workspace excludes
    /// `bevy_isomesh` and nothing local compiles it (M-293, I-008).
    fn bound(&self) -> FieldBound {
        self.inner.bound()
    }
}

/// One extraction at one precision, recentred and measured.
///
/// `R` is the precision under test. `reference` is the *same field* built in
/// `f64` and is what grades the result — the whole question is how far the `R`
/// arithmetic drifted from the truth, and only `f64` can say.
///
/// `widen` is passed in rather than taken from [`Real`], and that is not a
/// stylistic choice. `Real`'s only downward conversion is
/// [`as_f32`](Real::as_f32), which is *lossy for `f64`* — grading through it
/// quantises the `f64` vertices to `f32` before measuring them, so both columns
/// report the same number and the demo silently compares `f32` with itself. That
/// is exactly what the first version of this file did: it printed `5.8564` for
/// both sides at 2²³ and looked entirely plausible. The caller knows the concrete
/// type and can widen without loss, so it does.
fn extract<R, F, G>(
    field: &F,
    reference: &G,
    widen: impl Fn(R) -> f64,
    offset: f64,
    samples: u32,
) -> Option<Measured>
where
    R: Real,
    F: Sdf<Scalar = R> + ReferenceField,
    G: Sdf<Scalar = f64>,
{
    let (min, max) = field.domain();
    let cell = (max[0] - min[0]) / R::from_f64(f64::from(samples - 1));
    let shape = match RuntimeShape3::new([samples; 3]) {
        Ok(shape) => shape,
        Err(error) => {
            error!("grid {samples}^3 rejected: {error}");
            return None;
        }
    };

    let mut buffer: MeshBuffer<R> = MeshBuffer::new();
    let started = Instant::now();
    if let Err(error) = DualContouring::<R>::new().extract(field, &shape, min, cell, &mut buffer) {
        error!("extraction failed at {samples}^3, offset {offset:e}: {error}");
        return None;
    }
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;

    let h = cell.as_f32() as f64;
    let cfg = match ValidateConfig::from_cell_size(h) {
        Ok(cfg) => cfg,
        Err(error) => {
            error!("cell size {h} is not a usable spacing: {error}");
            return None;
        }
    };
    let (report, _) = validate_features(&buffer.positions, &buffer.indices, &cfg);

    // Grade in f64, against the same field, widening losslessly -- see `widen`.
    let mut worst_off_surface = 0.0f64;
    for p in &buffer.positions {
        let q = [widen(p[0]), widen(p[1]), widen(p[2])];
        let off = (reference.sample(q) / h).abs();
        if off.is_finite() && off > worst_off_surface {
            worst_off_surface = off;
        }
    }

    // Recentre in the precision that produced the mesh. The f32 side subtracts
    // its offset in f32 because that is what a consumer who chose f32 gets; the
    // f64 side subtracts in f64 and narrows afterwards, which is what
    // `Real::as_f32`'s doc tells a CAD consumer to do.
    let shift = R::from_f64(offset);
    let positions: Vec<[f32; 3]> = buffer
        .positions
        .iter()
        .map(|p| {
            [
                (p[0] - shift).as_f32(),
                (p[1] - shift).as_f32(),
                (p[2] - shift).as_f32(),
            ]
        })
        .collect();
    let normals: Vec<[f32; 3]> = buffer
        .normals
        .iter()
        .map(|n| [n[0].as_f32(), n[1].as_f32(), n[2].as_f32()])
        .collect();

    let mut mesh = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::default(),
    );
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
    mesh.insert_indices(Indices::U32(buffer.indices.clone()));

    Some(Measured {
        vertices: buffer.positions.len(),
        triangles: buffer.indices.len() / 3,
        euler: report.euler_characteristic,
        boundary_edges: report.boundary_edges,
        non_manifold_edges: report.non_manifold_edges,
        worst_off_surface,
        extract_ms,
        mesh,
    })
}

struct Built {
    narrow: Measured,
    wide: Measured,
    field_name: &'static str,
    /// Cell size of the `f64` side, which is the undamaged one.
    cell: f64,
    lines: Vec<String>,
}

/// Dispatch on the field index, building the *same* field twice — once at the
/// precision under test and once in `f64` to grade it.
fn build(demo: &Demo) -> Option<Built> {
    let offset = demo.offset();
    match demo.field {
        0 => pair(
            Sphere::<f32>::canonical(),
            Sphere::<f64>::canonical(),
            demo,
            offset,
        ),
        1 => pair(
            Torus::<f32>::canonical(),
            Torus::<f64>::canonical(),
            demo,
            offset,
        ),
        _ => pair(
            BoxExact::<f32>::canonical(),
            BoxExact::<f64>::canonical(),
            demo,
            offset,
        ),
    }
}

fn pair<A, B>(narrow_field: A, wide_field: B, demo: &Demo, offset: f64) -> Option<Built>
where
    A: Sdf<Scalar = f32> + ReferenceField,
    B: Sdf<Scalar = f64> + ReferenceField + Clone,
{
    let grader = Offset {
        inner: wide_field.clone(),
        by: offset,
    };
    let narrow = Offset {
        inner: narrow_field,
        by: offset as f32,
    };
    let wide = Offset {
        inner: wide_field,
        by: offset,
    };

    // Read the cell size off the f64 side before it is consumed: it is the
    // undamaged one, and every ratio in the HUD is against it.
    let (min, max) = wide.domain();
    let cell = (max[0] - min[0]) / f64::from(demo.samples - 1);

    let narrow = extract::<f32, _, _>(&narrow, &grader, f64::from, offset, demo.samples)?;
    let wide = extract::<f64, _, _>(&wide, &grader, |x| x, offset, demo.samples)?;

    Some(Built {
        lines: hud(demo, &narrow, &wide, cell),
        narrow,
        wide,
        field_name: A::NAME,
        cell,
    })
}

/// One representable step of `f32` at `x`.
///
/// `f32` carries a 24-bit significand, so the gap between neighbours at `x` is
/// `2^(⌊log₂x⌋ − 23)`. This is the whole mechanism: it is the smallest distance
/// the format can tell apart *at that offset*, and every number below is it
/// compared against something.
fn ulp_f32(x: f64) -> f64 {
    if x <= 0.0 || !x.is_finite() {
        return 0.0;
    }
    (2.0f64).powi(x.abs().log2().floor() as i32 - 23)
}

fn hud(demo: &Demo, narrow: &Measured, wide: &Measured, cell: f64) -> Vec<String> {
    let offset = demo.offset();
    let ulp = ulp_f32(offset);
    let ratio = ulp / cell;
    let torn = narrow.boundary_edges > 0 || narrow.euler != wide.euler;

    let verdict = if torn {
        "CRACKED -- f32 has real holes here, and f64 does not."
    } else if narrow.worst_off_surface > wide.worst_off_surface * 4.0 {
        "BLURRED -- topology still perfect, accuracy already gone."
    } else {
        "both intact -- f32 has precision to spare at this offset."
    };

    let mut lines = vec![
        format!(
            "{:<20} {:>14}   [-] nearer, [=] further",
            "offset",
            format!("2^{} = {:.0}", demo.log2_offset, offset)
        ),
        format!("{:<20} {:>14.6}   cell size", "h", cell),
        format!(
            "{:<20} {:>14.6}   one representable f32 step at this offset",
            "ulp(f32)", ulp
        ),
        String::new(),
        "two failures, two laws, and they are NOT the same knob:".into(),
        format!(
            "{:<20} {:>14.3}   ACCURACY is relative -- worst |f|/h tracks this,",
            "ulp / h", ratio
        ),
        "                                        press [ and this doubles while the".into(),
        "                                        offset stays put -- and so does |f|/h.".into(),
        format!(
            "{:<20} {:>14}   TOPOLOGY is absolute -- it tears when ulp >= 1,",
            "ulp >= 1 ?",
            if ulp >= 1.0 { "YES" } else { "no" }
        ),
        "                                        i.e. offset >= 2^23, and that does".into(),
        "                                        NOT move when h does.".into(),
        String::new(),
        "                            f32              f64".into(),
        format!(
            "{:<20} {:>14} {:>16}",
            "vertices", narrow.vertices, wide.vertices
        ),
        format!(
            "{:<20} {:>14} {:>16}",
            "triangles", narrow.triangles, wide.triangles
        ),
        format!(
            "{:<20} {:>14} {:>16}",
            "euler characteristic", narrow.euler, wide.euler
        ),
        format!(
            "{:<20} {:>14} {:>16}",
            "boundary edges (holes)", narrow.boundary_edges, wide.boundary_edges
        ),
        format!(
            "{:<20} {:>14} {:>16}",
            "non-manifold edges", narrow.non_manifold_edges, wide.non_manifold_edges
        ),
        format!(
            "{:<20} {:>14.4} {:>16.4}   cells off the true surface",
            "worst |f| / h", narrow.worst_off_surface, wide.worst_off_surface
        ),
        format!(
            "{:<20} {:>14.3} {:>16.3}   ms",
            "extract", narrow.extract_ms, wide.extract_ms
        ),
        String::new(),
        verdict.into(),
        String::new(),
    ];

    lines.extend([
        "the ticket said f32 cracks at ~1e6. it does not -- at 2^20 the f32 mesh".into(),
        "is topologically perfect, chi = 2 with zero holes, and 1.38 cells off.".into(),
        "press - down to 2^20 and watch the boundary count. the tear is at 2^23,".into(),
        "where one representable step reaches a whole world unit (M-112).".into(),
        String::new(),
        "no condition number here, though the ticket asked for one: it describes a".into(),
        "cell's normals, and translating the field leaves those alone. it would".into(),
        "have looked relevant and explained nothing.".into(),
        String::new(),
        format!(
            "{} samples/axis   [ and ] to change   1-3 field",
            demo.samples
        ),
    ]);
    lines
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    demo: Res<Demo>,
    mut stats: ResMut<DemoStats>,
    mut meshes: ResMut<Assets<Mesh>>,
    materials: Res<Materials>,
    mut commands: Commands,
    mut query: Query<(&mut Mesh3d, &Side), With<DemoMesh>>,
    mut flags: ResMut<ViewFlags>,
    mut last: Local<Option<(usize, u32, i32)>>,
) {
    let key = (demo.field, demo.samples, demo.log2_offset);
    if *last == Some(key) && !flags.remesh_requested {
        return;
    }
    *last = Some(key);
    flags.remesh_requested = false;

    let Some(built) = build(&demo) else {
        return;
    };

    // One CSV line per re-mesh, so sweeping offset and resolution is a shell
    // loop rather than a stack of screenshots read by eye. The two laws in the
    // module header are checked against exactly this output.
    info!(
        "precision,{},{},{},{:.6},{},{},{},{},{:.4},{:.4}",
        built.field_name,
        demo.samples,
        demo.log2_offset,
        ulp_f32(demo.offset()) / built.cell,
        built.narrow.euler,
        built.wide.euler,
        built.narrow.boundary_edges,
        built.wide.boundary_edges,
        built.narrow.worst_off_surface,
        built.wide.worst_off_surface,
    );

    stats.title = format!(
        "E-112  precision   offset 2^{} = {:.0}   field {} ({})   {}^3",
        demo.log2_offset,
        demo.offset(),
        demo.field + 1,
        built.field_name,
        demo.samples,
    );
    // The f32 side is the subject; its counts are the ones in the standard rows.
    stats.vertices = built.narrow.vertices;
    stats.triangles = built.narrow.triangles;
    stats.extract_ms = built.narrow.extract_ms;
    stats.extra = built.lines;

    // Offset along x so the two read as a comparison rather than as one object.
    let span = (built.cell * f64::from(demo.samples - 1)) as f32;
    let shift = span * 0.62;
    let narrow = meshes.add(built.narrow.mesh);
    let wide = meshes.add(built.wide.mesh);

    if query.is_empty() {
        commands.spawn((
            Mesh3d(narrow.clone()),
            MeshMaterial3d(materials.f32_side.clone()),
            Transform::from_xyz(-shift, 0.0, 0.0),
            DemoMesh,
            Side::F32,
        ));
        commands.spawn((
            Mesh3d(wide.clone()),
            MeshMaterial3d(materials.f64_side.clone()),
            Transform::from_xyz(shift, 0.0, 0.0),
            DemoMesh,
            Side::F64,
        ));
    } else {
        for (mut handle, side) in &mut query {
            handle.0 = match side {
                Side::F32 => narrow.clone(),
                Side::F64 => wide.clone(),
            };
        }
    }
}
