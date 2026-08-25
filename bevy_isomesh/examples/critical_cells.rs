//! E-304 — the defect is in the sign lattice, before the mesh exists.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example critical_cells --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 20-50x slower, and this one
//! samples the field twice.
//!
//! `1`-`5` switch field, `F` flies to the densest cluster and back, `H` hides
//! the surface. The rest are the shared keys — `W` wireframe, `G` domain box,
//! `F12` screenshot.
//!
//! Under `ISOMESH_CAPTURE` it needs no keyboard: the field advances every 16
//! captured frames and the camera dives into the cluster within each field, so
//! `record_gif.sh`'s default 80 frames is exactly one pass through all five
//! fields. `ISOMESH_SPIN=0.004` adds yaw on top, `ISOMESH_FIELD` pins one field
//! and `ISOMESH_SURFACE=off` starts with the surface hidden.
//!
//! ```bash
//! # COLORS=64 because the dive and the spin move the whole frame, which defeats
//! # the inter-frame compression a GIF relies on: 4.73 MB against 6.34 MB at the
//! # default 256, and the palette here is a dark grey, three saturated gizmo
//! # colours and white text, so 64 costs nothing visible.
//! ISOMESH_SPIN=0.004 COLORS=64 scripts/record_gif.sh critical_cells docs/gifs/e304.gif
//! ```
//!
//! # What is on screen
//!
//! - **Translucent grey** — the Dual Contouring surface. Context, not subject;
//!   `H` removes it.
//! - **Cyan cages** — cells whose `2x2x2` corner signs host a **2D-critical**
//!   configuration: some `2x2` face of the cell is a checkerboard, so its two
//!   inside corners are one diagonal of that face and its two outside corners
//!   are the other. The inside pair shares only a cell *edge*.
//! - **Magenta cages** — cells hosting a **3D-critical** configuration: the
//!   inside set is exactly two corners differing in all three coordinates, so
//!   they share only a cell *vertex*; or the complementary six-corner case.
//! - **Yellow dots** — the mesh's actual **non-manifold vertices**, from
//!   `isomesh::validate::validate_features`.
//!
//! Each cage is drawn at the exact bounds of its cell, so a dot inside a cage is
//! a dot inside that cell and nothing has been fudged to make them line up.
//!
//! # The demo is the coincidence
//!
//! M-338 / P-41 (`docs/experiments/p-41.csv`) measured this at 65 samples per
//! axis on the eight reference fields. Every yellow dot sits inside a cage,
//! every cage holds exactly one dot, and the two counts on the HUD are the same
//! number:
//!
//! | field | 2D-critical | 3D-critical | total | non-manifold vertices |
//! |---|---:|---:|---:|---:|
//! | `noise_cavity` | 567 | 35 | **602** | **602** |
//! | `gyroid` | 132 | 9 | **141** | **141** |
//! | `fbm_terrain` | 58 | 0 | **58** | **58** |
//! | `sphere`, `csg_difference` | 0 | 0 | **0** | **0** |
//!
//! Co-location is 2442/2442 = **100%** pooled over both dual extractors, against
//! a chance baseline — the share of *vertex-hosting* cells that are critical — of
//! 0.66% to 2.1%. So the relation is not a correlation and not a detector that
//! fires near trouble. The census does not predict where the non-manifold
//! vertices are; it **counts** them, from the sign bytes alone, before a single
//! triangle is emitted.
//!
//! The two clean fields are half of the demo rather than filler. A census that
//! only ever came back non-zero would say nothing about whether it is measuring
//! the defect or the field's complexity; `sphere` and `csg_difference` read zero
//! across the whole row, and their extractors produce zero non-manifold output,
//! so the implication runs both ways.
//!
//! # The 256-entry table is enumerated here, not transcribed
//!
//! [`sign_table`] walks all 256 sign bytes and decides each one combinatorially
//! from the definitions above — face sets built from the axis they fix,
//! diagonality from an XOR of corner indices. `CLAUDE.md`'s rule 5 forbids
//! writing a case table from memory, and a 256-entry table of critical
//! configurations is exactly the kind of thing that would be wrong in one entry
//! and unfalsifiable on screen.
//!
//! The enumeration is checked at startup and the three numbers are logged: **120
//! 2D-critical, 8 3D-critical, 0 in both** — the 8 being four main diagonals
//! times two complements. Disjointness is then re-checked against the live
//! census on every field, because a table that is disjoint over abstract bytes
//! and a census that never sees a cell in both classes are two different claims.
//!
//! # Inside is `value < 0.0`, and that is not the IEEE sign bit
//!
//! Copied from `cube.rs::is_inside` rather than reinvented: the census has to
//! partition the samples the *extractor's* way or it is a census of a different
//! lattice. `-0.0` has its sign bit set and `-0.0 < 0.0` is **false**, so a
//! negative zero is outside. That is reachable rather than theoretical —
//! `box_exact` is exactly zero across its whole boundary — and a census built on
//! `is_sign_negative` would flag cells the extractor never treats as mixed.
//!
//! # `f64`, and Dual Contouring only
//!
//! M-338 was measured in `f64`, so the numbers on the HUD are reproducible only
//! in `f64`; the surface is cast to `f32` on its way into the [`Mesh`] asset and
//! nothing but the picture depends on that.
//!
//! Dual Contouring rather than Surface Nets because its `Clamp::ToCell` insets
//! by `(1 - eps)` about the cell centre, so mapping a vertex back to its cell by
//! `floor` is exact. M-338 measured the difference and it is real: Surface Nets'
//! centroid rule places 5,400 of 5,768 vertices *exactly on* a cell boundary on
//! `box_exact`, where `floor` names the neighbour. It changed no attribution
//! there, and this example does not need to rely on that.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::dual_contouring::DualContouring;
use isomesh::fields::{
    FbmTerrain, ReferenceField, Sphere, capped_gyroid, csg_difference, noise_cavity,
};
use isomesh::validate::{ValidateConfig, validate_features};
use isomesh::{MeshBuffer, MeshSink, RuntimeShape3, Sdf};

/// The registered resolution. `65` samples span `64` cells per axis, and it is
/// the grid every number in this file's table was measured on.
const DEFAULT_SAMPLES: u32 = 65;

/// Below this a 2x2x2 block does not exist.
const MIN_SAMPLES: u32 = 3;

/// Above this the census allocates more than the demo is worth.
const MAX_SAMPLES: u32 = 129;

/// The fields offered, in the order the digit keys select them.
///
/// `noise_cavity` leads because it is the densest in the sweep, and the last two
/// are the clean controls — a census that never returns zero would be
/// unfalsifiable on screen.
const FIELD_COUNT: usize = 5;

// ─── the 256-byte classification, enumerated ────────────────────────────────

/// Which of the 256 possible cell sign bytes host each critical configuration.
///
/// Corner bit layout: bit `i` is corner `(x, y, z)` with `i = x + 2y + 4z`, so
/// two corners are cell-diagonal exactly when `i ^ j == 0b111`, and two corners
/// of the face fixing axis `a` are face-diagonal exactly when
/// `i ^ j == 0b111 ^ (1 << a)`.
#[derive(Resource)]
struct SignTable {
    /// A checkerboard `2x2` face: the inside pair shares only a cell edge.
    two_d: [bool; 256],
    /// A main-diagonal inside pair, or its complement: the pair shares only a
    /// cell vertex.
    three_d: [bool; 256],
}

/// Inside, the way the extractors decide it.
///
/// Copied from `cube.rs::is_inside`: strictly negative. See the module docs for
/// why this is not `is_sign_negative`.
fn is_inside(value: f64) -> bool {
    value < 0.0
}

/// Exactly two inside corners, differing in all three coordinates.
fn is_vertex_diagonal_pair(byte: u32) -> bool {
    let mut found = [0u32; 8];
    let mut n = 0usize;
    for corner in 0..8u32 {
        if (byte >> corner) & 1 == 1 {
            found[n] = corner;
            n += 1;
        }
    }
    n == 2 && (found[0] ^ found[1]) == 0b111
}

/// Some `2x2` face of the cell is a checkerboard.
///
/// The face fixing axis `a` at side `s` is the four corners with bit `a` equal
/// to `s`; its two diagonals are the corner pairs with
/// `i ^ j == 0b111 ^ (1 << a)`. A checkerboard is two inside corners forming one
/// of those diagonals, which forces the other two — the outside pair — to be the
/// other.
fn has_checkerboard_face(byte: u32) -> bool {
    for axis in 0..3u32 {
        let diagonal = 0b111 ^ (1 << axis);
        for side in 0..2u32 {
            let mut inside = [0u32; 4];
            let mut n = 0usize;
            for corner in 0..8u32 {
                if (corner >> axis) & 1 == side && (byte >> corner) & 1 == 1 {
                    inside[n] = corner;
                    n += 1;
                }
            }
            if n == 2 && (inside[0] ^ inside[1]) == diagonal {
                return true;
            }
        }
    }
    false
}

/// Decide all 256 sign bytes from the definitions, once, and say what came out.
///
/// The three numbers are logged rather than asserted because this is a demo a
/// stranger runs: a wrong table must be loud and must not take the window down
/// with it. M-338 gives 120 / 8 / 0, and anything else means the enumeration
/// above no longer matches Latecki's definitions.
fn sign_table() -> SignTable {
    let mut two_d = [false; 256];
    let mut three_d = [false; 256];
    for byte in 0..256u32 {
        two_d[byte as usize] = has_checkerboard_face(byte);
        // The complementary case is the same configuration seen from the other
        // side, and Latecki lists both: six inside corners whose two outside
        // corners share only a vertex is just as non-well-composed.
        three_d[byte as usize] =
            is_vertex_diagonal_pair(byte) || is_vertex_diagonal_pair(!byte & 0xFF);
    }

    let flat = two_d.iter().filter(|c| **c).count();
    let solid = three_d.iter().filter(|c| **c).count();
    let both = (0..256).filter(|b| two_d[*b] && three_d[*b]).count();
    if (flat, solid, both) == (120, 8, 0) {
        info!(
            "sign-byte classification, enumerated from the definitions: \
             {flat} 2D-critical, {solid} 3D-critical, {both} in both -- matches M-338"
        );
    } else {
        error!(
            "sign-byte classification is {flat} 2D-critical, {solid} 3D-critical, \
             {both} in both. M-338 measured 120 / 8 / 0, so this file's \
             has_checkerboard_face or is_vertex_diagonal_pair no longer matches \
             Latecki's definitions and every cage on screen is suspect."
        );
    }

    SignTable { two_d, three_d }
}

// ─── the grid ───────────────────────────────────────────────────────────────

/// The grid the field is sampled and meshed on.
///
/// One definition, so the census and the extraction cannot disagree about which
/// lattice they looked at.
struct Grid {
    /// World position of sample `[0, 0, 0]`.
    origin: [f64; 3],
    cell_size: f64,
    samples: u32,
    /// Cells per axis: `samples - 1`.
    cells: u32,
}

impl Grid {
    /// `i = x + y*sx + z*sx*sy`, the crate's order.
    fn sample_index(&self, x: usize, y: usize, z: usize) -> usize {
        let n = self.samples as usize;
        x + y * n + z * n * n
    }

    /// The same order over cells rather than samples.
    fn cell_index(&self, cell: [u32; 3]) -> usize {
        let c = self.cells as usize;
        cell[0] as usize + cell[1] as usize * c + cell[2] as usize * c * c
    }

    fn point(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        [
            self.origin[0] + self.cell_size * x as f64,
            self.origin[1] + self.cell_size * y as f64,
            self.origin[2] + self.cell_size * z as f64,
        ]
    }

    /// The cell a dual vertex belongs to.
    ///
    /// Exact rather than a guess: Dual Contouring's `Clamp::ToCell` insets by
    /// `(1 - eps)` about the cell centre, so no vertex lands on a face. The
    /// clamp to `0..cells-1` is belt and braces for a vertex that escaped
    /// anyway — it can only mis-attribute, never index out of bounds.
    fn cell_of(&self, p: [f64; 3]) -> usize {
        let mut cell = [0u32; 3];
        let last = self.cells.saturating_sub(1);
        for (axis, slot) in cell.iter_mut().enumerate() {
            let t = ((p[axis] - self.origin[axis]) / self.cell_size).floor();
            *slot = if t < 0.0 {
                0
            } else if t > f64::from(last) {
                last
            } else {
                t as u32
            };
        }
        self.cell_index(cell)
    }
}

// ─── the overlay ────────────────────────────────────────────────────────────

/// Everything drawn on top of the surface, resolved to world space once per
/// rebuild rather than per frame.
#[derive(Resource, Default)]
struct Overlay {
    /// Cell size in world units, which is also the cage size.
    cell_size: f32,
    /// Minimum corner of each 2D-critical cell.
    two_d: Vec<Vec3>,
    /// Minimum corner of each 3D-critical cell.
    three_d: Vec<Vec3>,
    /// The mesh's non-manifold vertices.
    markers: Vec<Vec3>,
}

/// Cages get their own group so they can be drawn in front of the translucent
/// surface without dragging the shared wireframe along with them.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct CageGizmos;

/// Markers get a third group, biased harder still, so a yellow dot is never lost
/// behind the cyan cage it sits inside. `manifold_check` earned the first of
/// these; a line lying on the surface z-fights and reads as intermittent, which
/// is indistinguishable from the defect being intermittent.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct MarkerGizmos;

// ─── camera framing ─────────────────────────────────────────────────────────

/// How many bins per axis the cluster search buckets critical cells into.
///
/// A bin is `1/8` of the domain, which at 65 samples is 8 cells across — the
/// scale a cluster of critical cells actually has, and coarse enough that the
/// pick is stable under a resolution change.
const CLUSTER_BINS: usize = 8;

/// Wide view radius, in domain extents.
///
/// From the field's own `domain()` rather than a fixed number, because the five
/// domains differ by 4x — the compact three are half-extent 2, the capped gyroid
/// is 7 and `fbm_terrain` is 8 — and a hardcoded radius puts the camera
/// comfortably *inside* the gyroid. E-110 found E-109's committed screenshot was
/// a picture of an inner wall for exactly that reason.
const WIDE_RADIUS_EXTENTS: f32 = 1.6;

/// Where the dolly ends on a field with no critical cells, in domain extents.
const CLEAN_RADIUS_EXTENTS: f32 = 1.05;

/// Cluster view radius, in bin widths. Wide enough to keep some surface in
/// frame, close enough that a cage is tens of pixels rather than five. Measured
/// on a 640x360 capture: at 3.0 the cells nearest the camera fill a third of the
/// frame each and the cluster reads as clutter, at 4.0 about thirty cells span
/// the width and one dot per cage is legible.
const CLUSTER_VIEW_BINS: f32 = 4.0;

/// Where the subject sits in frame, as a fraction of the orbit radius, right and
/// down from the centre.
///
/// **The HUD is fourteen lines in the upper left and the cluster is the
/// subject.** Centring it photographs the argument with its evidence hidden,
/// which is E-112's lesson and E-109's committed screenshot. Applied in the
/// camera's own basis rather than as a world offset, so it holds while
/// `ISOMESH_SPIN` yaws.
const SUBJECT_OFFSET: Vec2 = Vec2::new(0.20, 0.10);

/// The two camera framings a rebuild computes.
#[derive(Resource, Default)]
struct Framing {
    wide_focus: Vec3,
    wide_radius: f32,
    /// Equal to the wide framing when the census is empty — there is no cluster
    /// to fly to on a clean field, and inventing one would point the camera at
    /// an arbitrary corner.
    cluster_focus: Vec3,
    cluster_radius: f32,
}

/// Interactive view state. Under capture the blend is driven by frame number
/// instead, so a GIF is not a still frame.
#[derive(Resource, Default)]
struct Zoomed(bool);

/// Whether the translucent surface is drawn.
#[derive(Resource)]
struct ShowSurface(bool);

impl Default for ShowSurface {
    /// `ISOMESH_SURFACE=off` starts with the surface hidden, so the cage cloud
    /// on its own can be captured without a keyboard.
    ///
    /// The harness's own contract — `ISOMESH_VIEW`, `ISOMESH_FIELD`,
    /// `ISOMESH_SAMPLES` all exist for this reason — is that anything a capture
    /// depends on is reachable from the environment. A view that can only be
    /// reached by pressing `H` is a view no committed image can be regenerated
    /// from, which is the hole E-115 found in `ISOMESH_SAMPLES`.
    fn default() -> Self {
        Self(!matches!(
            std::env::var("ISOMESH_SURFACE")
                .unwrap_or_default()
                .trim()
                .to_ascii_lowercase()
                .as_str(),
            "off" | "0" | "no" | "hide"
        ))
    }
}

/// Captured frames spent on each field.
///
/// Five fields at 16 frames each is 80, which is `record_gif.sh`'s default
/// `ISOMESH_CAPTURE_FRAMES` — so the default capture is exactly one pass through
/// the census, ending where it started.
const CAPTURE_FRAMES_PER_FIELD: u32 = 16;

/// Which field is showing.
#[derive(Resource)]
struct Field(usize);

/// Samples per axis.
#[derive(Resource)]
struct Resolution(u32);

#[derive(Resource)]
struct SurfaceMaterial(Handle<StandardMaterial>);

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-304 critical cells".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<CageGizmos>()
        .init_gizmo_group::<MarkerGizmos>()
        // Enumerated here rather than in `setup`, so that every system can take
        // it as a plain `Res` -- a resource inserted by a `Startup` command is
        // present in `Update`, but only because of a schedule ordering rule, and
        // an `Option<Res<_>>` that guards against it reads like a fallback for a
        // state that cannot happen. `DefaultPlugins` has already installed the
        // log subscriber by this point, so the self-check is heard.
        .insert_resource(sign_table())
        .insert_resource(Resolution(
            common::samples_override()
                .unwrap_or(DEFAULT_SAMPLES)
                .clamp(MIN_SAMPLES, MAX_SAMPLES),
        ))
        .insert_resource(Field(0))
        .init_resource::<ShowSurface>()
        .init_resource::<Overlay>()
        .init_resource::<Framing>()
        .init_resource::<Zoomed>()
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (controls, rebuild, frame_camera, draw_overlay).chain(),
        )
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut gizmo_config: ResMut<GizmoConfigStore>,
    mut camera: Query<&mut OrbitCamera>,
) {
    // Off both axes on purpose. A cage seen down an axis is a square, and a
    // square hides the dot inside it -- three quarters of what makes the
    // coincidence readable is that the cages read as boxes.
    for mut orbit in &mut camera {
        orbit.yaw = 0.62;
        orbit.pitch = 0.32;
    }

    let (cages, _) = gizmo_config.config_mut::<CageGizmos>();
    cages.line.width = 1.6;
    cages.depth_bias = -0.4;

    let (marks, _) = gizmo_config.config_mut::<MarkerGizmos>();
    marks.line.width = 3.2;
    marks.depth_bias = -0.8;

    // The cages and the dots are inside the surface, so the surface is
    // translucent and double-sided or it hides its own subject — the same reason
    // `hermite_debug` is.
    commands.insert_resource(SurfaceMaterial(materials.add(StandardMaterial {
        base_color: Color::srgba(0.70, 0.74, 0.80, 0.22),
        perceptual_roughness: 0.5,
        alpha_mode: AlphaMode::Blend,
        double_sided: true,
        cull_mode: None,
        ..default()
    })));

    // Spawned rather than assumed: `draw_domain` queries for it, and without one
    // the `G` toggle silently does nothing. Filled in by the first rebuild.
    commands.spawn(DemoDomain {
        min: Vec3::splat(-1.0),
        max: Vec3::splat(1.0),
    });
}

/// Field selection and the two view toggles.
///
/// Under capture the field advances on frame count and the zoom is driven by
/// [`frame_camera`], because an example whose subject only changes on a keypress
/// captures as a still frame.
fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    mut field: ResMut<Field>,
    mut zoomed: ResMut<Zoomed>,
    mut surface: ResMut<ShowSurface>,
) {
    if capture.is_active() {
        field.0 = (capture.taken / CAPTURE_FRAMES_PER_FIELD) as usize % FIELD_COUNT;
    } else {
        field.0 = flags.field.min(FIELD_COUNT - 1);
        if keys.just_pressed(KeyCode::KeyF) {
            zoomed.0 = !zoomed.0;
        }
    }
    if keys.just_pressed(KeyCode::KeyH) {
        surface.0 = !surface.0;
    }
}

/// Census the sign lattice, mesh it, validate the mesh, and cross-tabulate.
#[allow(clippy::too_many_arguments)]
fn rebuild(
    field: Res<Field>,
    resolution: Res<Resolution>,
    flags: Res<ViewFlags>,
    table: Res<SignTable>,
    mut stats: ResMut<DemoStats>,
    mut overlay: ResMut<Overlay>,
    mut framing: ResMut<Framing>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut query: Query<&mut Mesh3d, With<DemoMesh>>,
    mut commands: Commands,
    material: Res<SurfaceMaterial>,
    mut domain: Query<&mut DemoDomain>,
    mut last: Local<Option<(usize, u32)>>,
) {
    let key = (field.0, resolution.0);
    if *last == Some(key) && !flags.remesh_requested {
        return;
    }
    *last = Some(key);

    let Some(built) = build(field.0, resolution.0, &table) else {
        return;
    };

    for mut d in &mut domain {
        d.min = built.domain_min;
        d.max = built.domain_max;
    }
    *framing = built.framing;
    *overlay = built.overlay;

    stats.title = format!(
        "E-304  critical cells - {}   {}^3   [1-5] field, F cluster, H surface",
        built.field_name, resolution.0,
    );
    stats.vertices = built.vertices;
    stats.triangles = built.triangles;
    stats.extract_ms = built.extract_ms;
    stats.extra = built.lines;

    let handle = meshes.add(built.mesh);
    if query.is_empty() {
        commands.spawn((Mesh3d(handle), MeshMaterial3d(material.0.clone()), DemoMesh));
    } else {
        for mut mesh in &mut query {
            mesh.0 = handle.clone();
        }
    }
}

/// Everything one field produced.
struct Built {
    mesh: Mesh,
    overlay: Overlay,
    framing: Framing,
    lines: Vec<String>,
    field_name: &'static str,
    domain_min: Vec3,
    domain_max: Vec3,
    vertices: usize,
    triangles: usize,
    extract_ms: f64,
}

/// Dispatch on the field index, then do the work once in [`census_and_mesh`].
///
/// The eight reference fields are eight different types, so
/// `for_each_reference_field!` cannot serve a runtime choice — the index is
/// matched here instead, the same shape `manifold_check` uses.
fn build(field: usize, samples: u32, table: &SignTable) -> Option<Built> {
    match field {
        0 => census_and_mesh(&noise_cavity::<f64>(), samples, table),
        1 => census_and_mesh(&capped_gyroid::<f64>(), samples, table),
        2 => census_and_mesh(&FbmTerrain::<f64>::canonical(), samples, table),
        3 => census_and_mesh(&Sphere::<f64>::canonical(), samples, table),
        _ => census_and_mesh(&csg_difference::<f64>(), samples, table),
    }
}

fn census_and_mesh<F>(field: &F, samples: u32, table: &SignTable) -> Option<Built>
where
    F: Sdf<Scalar = f64> + ReferenceField,
{
    let (min, max) = field.domain();
    let cells = samples.saturating_sub(1);
    if cells == 0 {
        error!("{} samples per axis leaves no cells to census", samples);
        return None;
    }
    let grid = Grid {
        origin: min,
        cell_size: (max[0] - min[0]) / f64::from(cells),
        samples,
        cells,
    };

    // ── the census: sample, take signs, classify every cell ─────────────────
    //
    // Before the extractor runs, and reading nothing the extractor produced.
    // That ordering is the claim this example exists to make legible.
    let census_started = Instant::now();
    let n = samples as usize;
    let mut inside = vec![false; n * n * n];
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                inside[grid.sample_index(x, y, z)] = is_inside(field.sample(grid.point(x, y, z)));
            }
        }
    }

    let c = cells as usize;
    let mut critical = vec![false; c * c * c];
    let mut two_d = Vec::new();
    let mut three_d = Vec::new();
    let mut both_classes = 0usize;
    let mut active_cells = 0usize;
    for cz in 0..c {
        for cy in 0..c {
            for cx in 0..c {
                let mut byte = 0u32;
                for corner in 0..8usize {
                    let dx = corner & 1;
                    let dy = (corner >> 1) & 1;
                    let dz = (corner >> 2) & 1;
                    if inside[grid.sample_index(cx + dx, cy + dy, cz + dz)] {
                        byte |= 1 << corner;
                    }
                }
                if byte != 0x00 && byte != 0xFF {
                    active_cells += 1;
                }
                let flat = table.two_d[byte as usize];
                let solid = table.three_d[byte as usize];
                if !flat && !solid {
                    continue;
                }
                if flat && solid {
                    both_classes += 1;
                }
                critical[grid.cell_index([cx as u32, cy as u32, cz as u32])] = true;
                let corner = grid.point(cx, cy, cz);
                let at = Vec3::new(corner[0] as f32, corner[1] as f32, corner[2] as f32);
                // The classes are disjoint over the 256 bytes, so this `if` is a
                // partition rather than a precedence. `both_classes` is what
                // would say otherwise, on real cells rather than abstract bytes.
                if solid {
                    three_d.push(at);
                } else {
                    two_d.push(at);
                }
            }
        }
    }
    let census_ms = census_started.elapsed().as_secs_f64() * 1000.0;
    if both_classes > 0 {
        error!(
            "{} cells were classified BOTH 2D- and 3D-critical. M-338 measures the \
             two classes as disjoint over all 256 sign bytes, so either the \
             enumeration or that finding is wrong.",
            both_classes
        );
    }

    // ── the mesh, and where it is non-manifold ──────────────────────────────
    let shape = match RuntimeShape3::new([samples; 3]) {
        Ok(shape) => shape,
        Err(error) => {
            error!("grid {samples}^3 rejected: {error}");
            return None;
        }
    };
    let mut buffer = MeshBuffer::<f64>::new();
    let extract_started = Instant::now();
    if let Err(error) = DualContouring::<f64>::new().extract(
        field,
        &shape,
        grid.origin,
        grid.cell_size,
        &mut buffer,
    ) {
        error!("dual contouring failed at {samples}^3: {error}");
        return None;
    }
    let extract_ms = extract_started.elapsed().as_secs_f64() * 1000.0;

    let cfg = match ValidateConfig::from_cell_size(grid.cell_size) {
        Ok(cfg) => cfg,
        Err(error) => {
            error!(
                "cell size {} is not a usable spacing: {error}",
                grid.cell_size
            );
            return None;
        }
    };
    let (_, features) = validate_features(&buffer.positions, &buffer.indices, &cfg);

    // ── co-location: which incidents happened in a critical cell ────────────
    //
    // M-338's registered rule. A non-manifold vertex names one cell, so there is
    // nothing to decide. A non-manifold edge joins two dual vertices and
    // therefore two cells, and it counts as co-located when *either* endpoint's
    // cell is flagged — that being the reading under which "this incident
    // occurred in a critical cell" is true.
    let cell_of_vertex: Vec<usize> = buffer.positions.iter().map(|p| grid.cell_of(*p)).collect();
    let at_vertex = |v: u32| -> Option<usize> { cell_of_vertex.get(v as usize).copied() };
    let vertex_hits = features
        .vertices
        .iter()
        .filter(|v| at_vertex(**v).is_some_and(|cell| critical[cell]))
        .count();
    let edge_hits = features
        .edges
        .iter()
        .filter(|e| {
            e.iter()
                .any(|v| at_vertex(*v).is_some_and(|cell| critical[cell]))
        })
        .count();
    let incidents = features.edges.len() + features.vertices.len();
    let hits = vertex_hits + edge_hits;

    // Distinct critical cells hosting a non-manifold vertex. Equal to both the
    // critical count and the vertex count is a bijection rather than a
    // coincidence, and that is the sentence this demo is a picture of.
    let mut hosting: Vec<usize> = features
        .vertices
        .iter()
        .filter_map(|v| at_vertex(*v))
        .filter(|cell| critical[*cell])
        .collect();
    hosting.sort_unstable();
    hosting.dedup();

    let critical_count = two_d.len() + three_d.len();
    let markers: Vec<Vec3> = features
        .vertices
        .iter()
        .filter_map(|v| buffer.positions.get(*v as usize))
        .map(|p| Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32))
        .collect();

    // ── framing ────────────────────────────────────────────────────────────
    let domain_min = Vec3::new(min[0] as f32, min[1] as f32, min[2] as f32);
    let domain_max = Vec3::new(max[0] as f32, max[1] as f32, max[2] as f32);
    let extent = domain_max.x - domain_min.x;
    let centre = (domain_min + domain_max) * 0.5;
    let bin = extent / CLUSTER_BINS as f32;
    let framing = match cluster(&two_d, &three_d, domain_min, extent) {
        Some(focus) => Framing {
            wide_focus: centre,
            wide_radius: extent * WIDE_RADIUS_EXTENTS,
            cluster_focus: focus,
            cluster_radius: bin * CLUSTER_VIEW_BINS,
        },
        // There is no cluster to fly to on a clean field, and inventing one
        // would point the camera at an arbitrary corner. It still has to move:
        // two of the five fields are clean, so a framing that held still would
        // spend two fifths of a capture on a frozen frame and read as a hung
        // demo rather than as a zero. So it dollies in on the same lerp, one
        // code path, and the HUD says why nothing is marked.
        None => Framing {
            wide_focus: centre,
            wide_radius: extent * WIDE_RADIUS_EXTENTS,
            cluster_focus: centre,
            cluster_radius: extent * CLEAN_RADIUS_EXTENTS,
        },
    };

    let colocation = if incidents == 0 {
        // Never 100% on an empty set. An empty numerator over an empty
        // denominator manufactures agreement out of silence, which is the
        // failure M-338's own harness records `n/a` to avoid.
        String::from("      n/a co-location (no incidents to place)")
    } else {
        format!(
            "{:>8.1}% co-location ({hits} of {incidents} incidents in critical cells)",
            100.0 * hits as f64 / incidents as f64
        )
    };

    let lines = vec![
        format!(
            "{:>9} samples/axis   {} cells   {census_ms:.1} ms census",
            samples,
            c * c * c
        ),
        String::new(),
        format!("{:>9} 2D-critical cells      (cyan cages)", two_d.len()),
        format!(
            "{:>9} 3D-critical cells      (magenta cages)",
            three_d.len()
        ),
        format!("{critical_count:>9} critical cells         total"),
        String::new(),
        format!(
            "{:>9} non-manifold vertices  (yellow dots)",
            features.vertices.len()
        ),
        format!("{:>9} non-manifold edges", features.edges.len()),
        format!("{:>9} critical cells hosting one", hosting.len()),
        colocation,
        String::new(),
        format!(
            "{:>9} sign-active cells        {:.3}% of them critical",
            active_cells,
            if active_cells == 0 {
                0.0
            } else {
                100.0 * critical_count as f64 / active_cells as f64
            }
        ),
        format!(
            "          sign table: {} 2D / {} 3D / {} both, enumerated at startup",
            table.two_d.iter().filter(|c| **c).count(),
            table.three_d.iter().filter(|c| **c).count(),
            (0..256usize)
                .filter(|b| table.two_d[*b] && table.three_d[*b])
                .count(),
        ),
    ];

    // The HUD is the evidence and a headless capture has no HUD to read. One
    // line per rebuild, so `ISOMESH_CAPTURE` leaves the census in the log where
    // a script can check it against M-338 -- E-203 learned this the hard way:
    // a measurement that only exists on screen cannot be verified from a
    // terminal.
    info!(
        "{} at {}^3: critical {} = {} 2D + {} 3D; non-manifold {} vertices, {} edges; \
         {} critical cells host one; co-location {}/{}; census {:.1} ms, extract {:.1} ms",
        F::NAME,
        samples,
        critical_count,
        two_d.len(),
        three_d.len(),
        features.vertices.len(),
        features.edges.len(),
        hosting.len(),
        hits,
        incidents,
        census_ms,
        extract_ms,
    );

    Some(Built {
        mesh: to_mesh(&buffer),
        overlay: Overlay {
            cell_size: grid.cell_size as f32,
            two_d,
            three_d,
            markers,
        },
        framing,
        lines,
        field_name: F::NAME,
        domain_min,
        domain_max,
        vertices: buffer.vertex_count(),
        triangles: buffer.triangle_count(),
        extract_ms,
    })
}

/// Where the critical cells are densest, or `None` when there are none.
///
/// Bucketed into `CLUSTER_BINS³` bins rather than searched pairwise: a quadratic
/// scan over every critical cell is fine at 602 and is not at 129 samples per
/// axis, and a camera target does not need more resolution than a bin. The focus
/// is the mean of the cells in the winning bin's neighbourhood, so a cluster
/// straddling a bin boundary is still centred rather than clipped.
fn cluster(two_d: &[Vec3], three_d: &[Vec3], origin: Vec3, extent: f32) -> Option<Vec3> {
    if two_d.is_empty() && three_d.is_empty() || extent <= 0.0 {
        return None;
    }
    let bin = extent / CLUSTER_BINS as f32;
    let all = || two_d.iter().chain(three_d.iter());
    let index = |p: &Vec3| -> usize {
        let mut out = 0usize;
        for axis in 0..3usize {
            let t = ((p[axis] - origin[axis]) / bin).max(0.0) as usize;
            out = out * CLUSTER_BINS + t.min(CLUSTER_BINS - 1);
        }
        out
    };

    let mut counts = vec![0u32; CLUSTER_BINS * CLUSTER_BINS * CLUSTER_BINS];
    for p in all() {
        counts[index(p)] += 1;
    }
    let winner = counts
        .iter()
        .enumerate()
        .max_by_key(|(i, count)| (**count, std::cmp::Reverse(*i)))
        .map(|(i, _)| i)?;
    let winner_centre = {
        let x = winner / (CLUSTER_BINS * CLUSTER_BINS);
        let y = (winner / CLUSTER_BINS) % CLUSTER_BINS;
        let z = winner % CLUSTER_BINS;
        origin + Vec3::new(x as f32 + 0.5, y as f32 + 0.5, z as f32 + 0.5) * bin
    };

    // One bin of slack on every side, so a knot split across a boundary is
    // averaged whole.
    let reach = bin * 1.5;
    let mut sum = Vec3::ZERO;
    let mut count = 0u32;
    for p in all() {
        if (*p - winner_centre).abs().max_element() <= reach {
            sum += *p;
            count += 1;
        }
    }
    if count == 0 {
        return Some(winner_centre);
    }
    Some(sum / count as f32)
}

/// The `f64` extraction as a Bevy mesh.
///
/// Cast rather than re-extracted in `f32`: the numbers on the HUD are M-338's
/// and they are `f64` numbers, so the mesh the picture is drawn from has to be
/// the one they were computed on.
fn to_mesh(buffer: &MeshBuffer<f64>) -> Mesh {
    let mut builder = MeshBuilder::new();
    for i in 0..buffer.positions.len() {
        let (Some(p), Some(n)) = (buffer.positions.get(i), buffer.normals.get(i)) else {
            continue;
        };
        builder.vertex(
            [p[0] as f32, p[1] as f32, p[2] as f32],
            [n[0] as f32, n[1] as f32, n[2] as f32],
        );
    }
    for t in buffer.indices.as_chunks::<3>().0 {
        let (Some(a), Some(b), Some(c)) = (t.first(), t.get(1), t.get(2)) else {
            continue;
        };
        builder.triangle(*a, *b, *c);
    }
    builder.into_mesh()
}

/// Point the orbit camera at the wide view or at the cluster.
///
/// Under capture the blend follows the frame counter, so the GIF dives into the
/// cluster once per field rather than sitting still. `paused` freezes it, which
/// is what makes a paused capture inspectable.
fn frame_camera(
    capture: Res<Capture>,
    framing: Res<Framing>,
    zoomed: Res<Zoomed>,
    flags: Res<ViewFlags>,
    mut camera: Query<&mut OrbitCamera>,
) {
    let t = if capture.is_active() {
        let phase = (capture.taken % CAPTURE_FRAMES_PER_FIELD) as f32
            / (CAPTURE_FRAMES_PER_FIELD.saturating_sub(1).max(1)) as f32;
        // Hold wide, ease in, arrive before the cut to the next field.
        let raw = ((phase - 0.15) / 0.7).clamp(0.0, 1.0);
        raw * raw * (3.0 - 2.0 * raw)
    } else if zoomed.0 {
        1.0
    } else {
        0.0
    };

    for mut orbit in &mut camera {
        if flags.paused {
            continue;
        }
        let target = framing.wide_focus.lerp(framing.cluster_focus, t);
        let radius = framing.wide_radius + (framing.cluster_radius - framing.wide_radius) * t;

        // The camera's own basis, from the same yaw/pitch the harness's
        // `orbit_camera` builds its transform from -- so the offset is exactly
        // one screen-space nudge however far the spin has turned. `orbit_camera`
        // places the eye at `focus + dir * radius`, so the view direction is
        // `-dir` and a focus moved along `-right` puts the target right of
        // centre.
        let dir = Vec3::new(
            orbit.yaw.cos() * orbit.pitch.cos(),
            orbit.pitch.sin(),
            orbit.yaw.sin() * orbit.pitch.cos(),
        );
        let forward = -dir;
        let right = forward.cross(Vec3::Y).normalize_or_zero();
        let up = right.cross(forward).normalize_or_zero();

        orbit.focus =
            target - right * (SUBJECT_OFFSET.x * radius) + up * (SUBJECT_OFFSET.y * radius);
        orbit.radius = radius;
    }
}

/// Draw the cages and the markers. Runs every frame; the overlay changes only on
/// rebuild.
fn draw_overlay(
    overlay: Res<Overlay>,
    surface: Res<ShowSurface>,
    mut visibility: Query<&mut Visibility, With<DemoMesh>>,
    mut cages: Gizmos<CageGizmos>,
    mut marks: Gizmos<MarkerGizmos>,
) {
    const CYAN: Color = Color::srgb(0.15, 0.85, 1.0);
    const MAGENTA: Color = Color::srgb(1.0, 0.25, 0.90);
    const YELLOW: Color = Color::srgb(1.0, 0.95, 0.20);

    // Written only when it differs. A `*visible = ...` every frame marks the
    // component changed every frame, and Bevy's visibility propagation is
    // change-driven — so an unconditional write turns a toggle nobody pressed
    // into per-frame work on every descendant.
    let wanted = if surface.0 {
        Visibility::Inherited
    } else {
        Visibility::Hidden
    };
    for mut visible in &mut visibility {
        if *visible != wanted {
            *visible = wanted;
        }
    }

    for corner in &overlay.two_d {
        cage(&mut cages, *corner, overlay.cell_size, CYAN);
    }
    for corner in &overlay.three_d {
        cage(&mut cages, *corner, overlay.cell_size, MAGENTA);
    }
    for marker in &overlay.markers {
        marks
            .sphere(
                Isometry3d::from_translation(*marker),
                overlay.cell_size * 0.24,
                YELLOW,
            )
            .resolution(6);
    }
}

/// The twelve edges of one cell, at its exact bounds.
///
/// Exact rather than inflated: a cage that is larger than its cell would make
/// every dot look contained whether or not it was, which is the one thing this
/// picture must not do. Corner indexing matches the extractor's — bit `i` of the
/// corner index is axis `i`.
fn cage(gizmos: &mut Gizmos<CageGizmos>, min: Vec3, size: f32, colour: Color) {
    let corner = |i: usize| {
        min + Vec3::new(
            if i & 1 == 0 { 0.0 } else { size },
            if i & 2 == 0 { 0.0 } else { size },
            if i & 4 == 0 { 0.0 } else { size },
        )
    };
    for i in 0..8usize {
        for axis in 0..3usize {
            let bit = 1 << axis;
            if i & bit == 0 {
                gizmos.line(corner(i), corner(i | bit), colour);
            }
        }
    }
}
