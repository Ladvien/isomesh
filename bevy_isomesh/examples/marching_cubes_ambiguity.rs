//! E-102 — the asymptotic decider, and how rarely it fires.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example marching_cubes_ambiguity --release
//! ```
//!
//! Left is plain Marching Cubes, right is the same extraction with
//! `FaceAmbiguity::AsymptoticDecider`. `1`–`3` switch field, `[` and `]` change
//! resolution, `A` toggles the cell markers, `W` shows the wireframe.
//!
//! # The framing this example was written with was wrong
//!
//! The examples catalog asks for "visible holes on the left, closed on the
//! right". **That cannot be shown, because it does not happen.** This crate's
//! case table is derived at compile time by walking each face
//! counter-clockwise, so a face's segments are a function of that face's own
//! four corner signs and two cells sharing a face cannot disagree — `✗11`, with
//! `validate_table()` checking all 256 cases. Neither side holes, at any
//! resolution, on any field.
//!
//! What the decider actually changes is *which* surface is built on an
//! ambiguous face: whether the two diagonally opposite inside corners are joined
//! across it or cut off separately. That is a topology change, and the HUD reads
//! it off as a **difference in Euler characteristic**.
//!
//! # The measurement this example exists to make visible
//!
//! Every box drawn is a cell with at least one ambiguous face — **amber where
//! the decider agreed with plain Marching Cubes and separated the corners,
//! magenta where it disagreed and joined them.** Magenta is the only place the
//! two meshes can differ, and there is very little of it.
//!
//! That is the result, and it is easier to believe when you can count the boxes.
//! Measured at 33³ over the seven reference fields (M-40): an ambiguous face
//! occurs on **0.515%** of the gyroid's surface cells, **1.532%** of
//! `fbm_terrain`'s, and on `sphere`, `torus`, `box_exact`, `csg_difference` and
//! `thin_plate` it **never occurs at all**. Press `3` for the sphere and the two
//! sides are not merely similar — the HUD reports them **byte-identical**, which
//! the committed golden fixture also pins.
//!
//! So the field picker deliberately offers a field where the rule does nothing.
//! An example that only ever showed the interesting case would misrepresent how
//! often the interesting case arrives.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{CappedGyroid, FbmTerrain, ReferenceField, Sphere, capped_gyroid};
use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::table::{AMBIGUOUS_FACES, is_inside};
use isomesh::marching_cubes::{FaceAmbiguity, MarchingCubes};
use isomesh::validate::{MeshReport, ValidateConfig, validate_indexed};
use isomesh::{RuntimeShape3, Sdf};

#[derive(Resource)]
struct Demo {
    field: usize,
    samples: u32,
    markers: bool,
}

#[derive(Resource)]
struct Materials {
    separate: Handle<StandardMaterial>,
    decider: Handle<StandardMaterial>,
}

#[derive(Component, Clone, Copy, PartialEq, Eq, Debug)]
enum Side {
    Separate,
    Decider,
}

/// One cell with an ambiguous face, in world space.
#[derive(Clone, Copy)]
struct AmbiguousCell {
    centre: Vec3,
    /// True where the decider joined at least one face, i.e. where the two
    /// meshes can actually differ.
    joined: bool,
}

/// Where the markers live between frames, so the gizmo pass does not re-sample
/// the field every frame.
#[derive(Resource, Default)]
struct Markers {
    cells: Vec<AmbiguousCell>,
    cell_size: f32,
    offset: f32,
}

struct Extraction {
    builder: MeshBuilder,
    report: MeshReport,
    millis: f64,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-102 marching cubes ambiguity".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Demo {
            field: 0,
            samples: common::samples_override().unwrap_or(33),
            markers: true,
        })
        .init_resource::<Markers>()
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, remesh, draw_markers))
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    for mut orbit in &mut camera {
        // Straight down -z at the pair, so neither side is nearer than the other
        // and perspective cannot make one look bigger than it is.
        orbit.yaw = std::f32::consts::FRAC_PI_2;
        orbit.pitch = 0.12;
        orbit.radius = 34.0;
    }
    commands.insert_resource(Materials {
        separate: materials.add(StandardMaterial {
            base_color: Color::srgb(0.78, 0.80, 0.86),
            perceptual_roughness: 0.45,
            ..default()
        }),
        decider: materials.add(StandardMaterial {
            base_color: Color::srgb(0.62, 0.80, 0.72),
            perceptual_roughness: 0.45,
            ..default()
        }),
    });
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    flags: Res<ViewFlags>,
    capture: Res<Capture>,
    mut demo: ResMut<Demo>,
) {
    if capture.is_active() {
        // While recording, sweep resolution in step with the captured frames so
        // the sequence shows the ambiguous-cell count changing with the grid
        // rather than with wall-clock time.
        let steps = [17u32, 21, 25, 29, 33, 41, 49];
        demo.samples = steps[(capture.taken as usize / 4) % steps.len()];
        return;
    }
    if keys.just_pressed(KeyCode::Digit1) {
        demo.field = 0;
    }
    if keys.just_pressed(KeyCode::Digit2) {
        demo.field = 1;
    }
    if keys.just_pressed(KeyCode::Digit3) {
        demo.field = 2;
    }
    if keys.just_pressed(KeyCode::KeyA) {
        demo.markers = !demo.markers;
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        demo.samples = (demo.samples - 4).max(9);
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        demo.samples = (demo.samples + 4).min(65);
    }
    if flags.field != demo.field {
        demo.field = flags.field;
    }
}

/// Extract the same field twice, once per rule, and census the ambiguous cells.
///
/// The census walks the same grid `extract` walks and asks the same two
/// questions the mesher asks — `AMBIGUOUS_FACES[case]`, then
/// [`joined_mask`] — rather than inferring anything from the meshes. Reading the
/// rule directly is the only way the marker can be trusted to mean what the
/// caption says it means.
fn extract_pair<F: Sdf<Scalar = f32> + ReferenceField>(
    field: &F,
    samples: u32,
) -> Option<(Extraction, Extraction, Vec<AmbiguousCell>, f32, f32)> {
    let (min, max) = field.domain();
    let cell_size = (max[0] - min[0]) / (samples - 1) as f32;
    let shape = RuntimeShape3::new([samples; 3]).ok()?;
    let cfg = ValidateConfig::from_cell_size(f64::from(cell_size)).ok()?;

    let extract_with = |rule: FaceAmbiguity| -> Option<Extraction> {
        let mut mesher = MarchingCubes::<f32>::new();
        mesher.set_face_ambiguity(rule);
        let mut builder = MeshBuilder::new();
        let started = Instant::now();
        mesher
            .extract(field, &shape, min, cell_size, &mut builder)
            .ok()?;
        let millis = started.elapsed().as_secs_f64() * 1000.0;
        let report = validate_indexed(builder.positions(), builder.indices(), &cfg);
        Some(Extraction {
            builder,
            report,
            millis,
        })
    };
    let separate = extract_with(FaceAmbiguity::Separate)?;
    let decider = extract_with(FaceAmbiguity::AsymptoticDecider)?;

    // ── the census ──────────────────────────────────────────────────────────
    let corner_offset = |c: u8| {
        [
            u32::from(c & 1),
            u32::from((c >> 1) & 1),
            u32::from((c >> 2) & 1),
        ]
    };
    let mut cells = Vec::new();
    let mut surface_cells = 0u32;
    for z in 0..samples - 1 {
        for y in 0..samples - 1 {
            for x in 0..samples - 1 {
                let mut corner_value = [0.0f32; 8];
                let mut case = 0u8;
                for (c, slot) in corner_value.iter_mut().enumerate() {
                    let o = corner_offset(c as u8);
                    *slot = field.sample([
                        min[0] + cell_size * (x + o[0]) as f32,
                        min[1] + cell_size * (y + o[1]) as f32,
                        min[2] + cell_size * (z + o[2]) as f32,
                    ]);
                    if is_inside(*slot) {
                        case |= 1 << c;
                    }
                }
                if case == 0 || case == 255 {
                    continue;
                }
                surface_cells += 1;
                let ambiguous = AMBIGUOUS_FACES[case as usize];
                if ambiguous == 0 {
                    continue;
                }
                cells.push(AmbiguousCell {
                    centre: Vec3::new(
                        min[0] + cell_size * (x as f32 + 0.5),
                        min[1] + cell_size * (y as f32 + 0.5),
                        min[2] + cell_size * (z as f32 + 0.5),
                    ),
                    joined: joined_mask(&corner_value, ambiguous) != 0,
                });
            }
        }
    }

    Some((separate, decider, cells, cell_size, surface_cells as f32))
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    demo: Res<Demo>,
    mut stats: ResMut<DemoStats>,
    mut markers: ResMut<Markers>,
    mut meshes: ResMut<Assets<Mesh>>,
    materials: Res<Materials>,
    mut commands: Commands,
    mut query: Query<(&mut Mesh3d, &mut Transform, &Side)>,
    mut camera: Query<&mut OrbitCamera>,
    mut last: Local<Option<(usize, u32)>>,
) {
    let key = (demo.field, demo.samples);
    if *last == Some(key) {
        return;
    }
    *last = Some(key);

    let (name, extracted, width) = match demo.field {
        0 => {
            let field: CappedGyroid<f32> = capped_gyroid();
            let (lo, hi) = field.domain();
            ("gyroid", extract_pair(&field, demo.samples), hi[0] - lo[0])
        }
        1 => {
            let field = FbmTerrain::<f32>::canonical();
            let (lo, hi) = field.domain();
            (
                "fbm_terrain",
                extract_pair(&field, demo.samples),
                hi[0] - lo[0],
            )
        }
        _ => {
            let field = Sphere::<f32>::canonical();
            let (lo, hi) = field.domain();
            ("sphere", extract_pair(&field, demo.samples), hi[0] - lo[0])
        }
    };
    let Some((separate, decider, cells, cell_size, surface_cells)) = extracted else {
        return;
    };

    let joined = cells.iter().filter(|c| c.joined).count();
    // Byte-identical is the honest headline on a field with no ambiguous face,
    // and it is exactly what the golden fixture pins for five of the seven.
    let identical = separate.builder.positions() == decider.builder.positions()
        && separate.builder.indices() == decider.builder.indices();

    // ASCII only: the HUD font has no em-dash and no superscript, and a
    // previous example shipped with both rendering as empty boxes.
    stats.title = format!("E-102 ambiguity - {name} at {}^3", demo.samples);
    stats.vertices = decider.builder.positions().len();
    stats.triangles = decider.builder.indices().len() / 3;
    stats.extract_ms = decider.millis;
    stats.extra = vec![
        "                separate     decider".to_string(),
        format!(
            "triangles   {:>10}  {:>10}",
            separate.builder.indices().len() / 3,
            decider.builder.indices().len() / 3
        ),
        format!(
            "chi         {:>10}  {:>10}",
            separate.report.euler_characteristic, decider.report.euler_characteristic
        ),
        format!(
            "non-manif.  {:>10}  {:>10}",
            separate.report.non_manifold_edges, decider.report.non_manifold_edges
        ),
        format!(
            "extract ms  {:>10.3}  {:>10.3}",
            separate.millis, decider.millis
        ),
        String::new(),
        format!(
            "surface cells {surface_cells:>8.0}   ambiguous {:>5} ({:.3}%)",
            cells.len(),
            100.0 * cells.len() as f32 / surface_cells.max(1.0)
        ),
        format!(
            "magenta: decider JOINED {joined:>5}   amber: agreed {:>5}",
            cells.len() - joined
        ),
        String::new(),
        // Three outcomes, not two. A mesh can differ without chi differing --
        // the gyroid at 33^3 does exactly that -- and reporting "chi differs by
        // 0" alongside "different surface" is a contradiction on screen.
        match (
            identical,
            decider.report.euler_characteristic - separate.report.euler_characteristic,
        ) {
            (true, _) => "the two meshes are BYTE-IDENTICAL - the rule never fired".to_string(),
            (false, 0) => "same chi, different mesh <- the decider re-paired a face without \
                 changing the topology"
                .to_string(),
            (false, d) => format!(
                "chi differs by {} <- a different surface, not a different fan",
                d.abs()
            ),
        },
    ];

    let offset = width * 0.62;
    for mut orbit in &mut camera {
        // Tight enough that the markers are findable; the pair is what is
        // being looked at, not the empty space round it.
        orbit.radius = width * 2.05;
    }
    markers.cells = cells;
    markers.cell_size = cell_size;
    markers.offset = offset;

    let pairs = [
        (Side::Separate, separate.builder, -offset),
        (Side::Decider, decider.builder, offset),
    ];
    if query.is_empty() {
        for (side, builder, x) in pairs {
            let material = match side {
                Side::Separate => materials.separate.clone(),
                Side::Decider => materials.decider.clone(),
            };
            commands.spawn((
                Mesh3d(meshes.add(builder.into_mesh())),
                MeshMaterial3d(material),
                Transform::from_xyz(x, 0.0, 0.0),
                DemoMesh,
                side,
            ));
        }
    } else {
        let mut built: Vec<(Side, Handle<Mesh>, f32)> = Vec::new();
        for (side, builder, x) in pairs {
            built.push((side, meshes.add(builder.into_mesh()), x));
        }
        for (mut mesh, mut transform, side) in &mut query {
            if let Some((_, handle, x)) = built.iter().find(|(s, _, _)| s == side) {
                mesh.0 = handle.clone();
                transform.translation.x = *x;
            }
        }
    }
}

/// A box round every cell with an ambiguous face, on both halves.
///
/// Drawn on both so the eye can compare the same cell under the two rules —
/// a marker on one side only would leave you hunting for its partner.
fn draw_markers(demo: Res<Demo>, markers: Res<Markers>, mut gizmos: Gizmos) {
    if !demo.markers || markers.cells.is_empty() {
        return;
    }
    // A touch larger than the cell, so the box reads as a marker rather than
    // disappearing into the surface it sits on.
    let size = Vec3::splat(markers.cell_size * 1.6);
    for cell in &markers.cells {
        let colour = if cell.joined {
            Color::srgb(0.95, 0.25, 0.85)
        } else {
            Color::srgb(0.95, 0.72, 0.20)
        };
        for x in [-markers.offset, markers.offset] {
            gizmos.cube(
                Transform::from_translation(cell.centre + Vec3::X * x).with_scale(size),
                colour,
            );
        }
    }
}
