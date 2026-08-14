//! E-107 — two levels of detail meeting, with and without transition cells.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example transvoxel_seams --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 20-50x slower.
//!
//! `T` transition cells on/off · `1`-`4` field · `[` `]` resolution ·
//! `W` wireframe · `F12` screenshot.
//!
//! The two committed screenshots are the toggle's two states, and they are
//! reproducible from a command line rather than from a keypress:
//!
//! ```bash
//! ISOMESH_FIELD=3 ISOMESH_SAMPLES=6 ISOMESH_TRANSITIONS=0 \
//!   ISOMESH_VIEW=nogrid ISOMESH_CAPTURE=<dir> cargo run --example transvoxel_seams --release
//! ```
//!
//! The gyroid at six coarse cells is the field the crack is *visible* on. On a
//! sphere at the same settings the gap is real and about `0.03` world units
//! wide — the HUD counts it, and you cannot see it. High curvature across the
//! seam is what turns a countable crack into a visible one.
//!
//! # What this shows
//!
//! Two blocks of the same field meshed at different resolutions, meeting on one
//! plane. The left is full resolution, the right is half. Meshed independently
//! they **do not meet**: the fine side ends on a contour of `2n x 2n` sub-squares
//! and the coarse side on one of `n x n` squares, and the difference is a ring of
//! unmatched boundary edges — a crack you can see the background through.
//!
//! The HUD counts them, so the toggle is a number and not an impression. At the
//! default settings it reads **88 without, 0 with** (M-76), and the count is taken
//! the same way the test takes it: boundary edges lying **wholly in the seam
//! plane**, because both blocks are legitimately open at their outer borders and a
//! global count would drown the signal.
//!
//! # Why the transition cells are not flat
//!
//! They have a **width** — Lengyel's `w(k) = 2^(k-2)`, half a coarse cell — and
//! the coarse block's boundary cells are scaled inward by the same amount to make
//! room (his Equation 4.2, M-77).
//!
//! A zero width also closes the crack, and closing it is *all* it does: M-74
//! measured a zero-width patch as **exactly** perpendicular to the surface, so it
//! shades as a hard crease. That is §4.3's *"severe shading problems"*, and it is
//! why this example ships with a real width rather than the simpler thing.

mod common;

use bevy::prelude::*;
use bevy_isomesh::to_bevy_mesh;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{ReferenceField, Sphere, Torus, capped_gyroid, csg_difference};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::transvoxel::cell::TransitionCell;
use isomesh::transvoxel::inset::{face_bit, inset_boundary};
use isomesh::validate::{ValidateConfig, validate_features};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

#[derive(Resource)]
struct Demo {
    field: usize,
    /// Coarse cells per axis on the half-resolution block.
    coarse_cells: u32,
    transitions: bool,
}

const FIELDS: [&str; 4] = ["sphere", "torus", "csg_difference", "gyroid"];
const MIN_CELLS: u32 = 2;
const MAX_CELLS: u32 = 24;

#[derive(Component, Clone, Copy, PartialEq, Eq)]
enum Part {
    Fine,
    Coarse,
    Transition,
}

#[derive(Resource)]
struct Materials {
    fine: Handle<StandardMaterial>,
    coarse: Handle<StandardMaterial>,
    transition: Handle<StandardMaterial>,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-107 transvoxel seams".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Demo {
            field: std::env::var("ISOMESH_FIELD")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .map_or(0, |v| v.min(FIELDS.len() - 1)),
            coarse_cells: common::samples_override().unwrap_or(4),
            // Settable for a capture, the same way ISOMESH_FIELD and
            // ISOMESH_VIEW are -- the two screenshots this example owes are the
            // toggle's two states, and pressing a key is not reproducible.
            transitions: std::env::var("ISOMESH_TRANSITIONS").map_or(true, |v| v != "0"),
        })
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, remesh))
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    for mut orbit in &mut camera {
        // Across the seam, not along it: `x` runs left-right so both blocks are
        // in frame, and slightly off-axis so a crack shows as background rather
        // than as a shading change.
        orbit.yaw = std::f32::consts::FRAC_PI_2 + 0.30;
        orbit.pitch = 0.26;
        orbit.radius = 6.2;
    }
    commands.insert_resource(Materials {
        fine: materials.add(StandardMaterial {
            base_color: Color::srgb(0.86, 0.74, 0.58),
            perceptual_roughness: 0.45,
            ..default()
        }),
        coarse: materials.add(StandardMaterial {
            base_color: Color::srgb(0.72, 0.78, 0.86),
            perceptual_roughness: 0.45,
            ..default()
        }),
        // The stitch itself, so you can see what is doing the work.
        transition: materials.add(StandardMaterial {
            base_color: Color::srgb(0.92, 0.52, 0.32),
            perceptual_roughness: 0.35,
            ..default()
        }),
    });
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    mut demo: ResMut<Demo>,
    capture: Res<Capture>,
    mut flags: ResMut<ViewFlags>,
) {
    if capture.is_active() {
        return;
    }
    if keys.just_pressed(KeyCode::KeyT) {
        demo.transitions = !demo.transitions;
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
    if keys.just_pressed(KeyCode::BracketRight) {
        demo.coarse_cells = (demo.coarse_cells + 2).min(MAX_CELLS);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        demo.coarse_cells = demo.coarse_cells.saturating_sub(2).max(MIN_CELLS);
    }
    if keys.just_pressed(KeyCode::KeyR) {
        flags.remesh_requested = true;
    }
}

struct Built {
    /// The field's extent, so the camera can size itself — the gyroid's domain is
    /// several times the others' and a fixed radius puts the camera inside it.
    extent: f32,
    fine: MeshBuffer<f32>,
    coarse: MeshBuffer<f32>,
    transition: MeshBuffer<f32>,
    gaps: usize,
    fine_h: f32,
    width: f32,
}

/// Mesh both blocks, optionally stitch them, and count what is left open in the
/// seam plane.
fn build<F: Sdf<Scalar = f32> + ReferenceField>(
    field: &F,
    coarse_cells: u32,
    transitions: bool,
) -> Option<Built> {
    let (lo, hi) = field.domain();
    let extent = hi[0] - lo[0];
    // The half-resolution block owns the upper half of the domain on x.
    let coarse_h = (extent * 0.5) / coarse_cells as f32;
    let fine_h = coarse_h * 0.5;
    let fine_cells = coarse_cells * 2;
    let seam_x = lo[0] + extent * 0.5;
    let width = if transitions { fine_h } else { 0.0 };

    // In-plane extent: the whole domain, at each block's own spacing.
    let fine_shape =
        RuntimeShape3::new([fine_cells + 1, 2 * fine_cells + 1, 2 * fine_cells + 1]).ok()?;
    let mut fine = MeshBuffer::<f32>::new();
    MarchingCubes::<f32>::new()
        .extract(field, &fine_shape, lo, fine_h, &mut fine)
        .ok()?;

    let coarse_shape =
        RuntimeShape3::new([coarse_cells + 1, 2 * coarse_cells + 1, 2 * coarse_cells + 1]).ok()?;
    let coarse_origin = [seam_x, lo[1], lo[2]];
    let mut coarse = MeshBuffer::<f32>::new();
    MarchingCubes::<f32>::new()
        .extract(field, &coarse_shape, coarse_origin, coarse_h, &mut coarse)
        .ok()?;

    let mut transition = MeshBuffer::<f32>::new();
    if transitions {
        // Equation 4.2 on the coarse block's low-x face -- the only one with a
        // finer neighbour.
        inset_boundary(
            &mut coarse,
            coarse_origin,
            coarse_cells,
            coarse_h,
            width,
            face_bit(0, 0),
        )
        .ok()?;

        for jz in 0..2 * coarse_cells as i64 {
            for jy in 0..2 * coarse_cells as i64 {
                let cell = TransitionCell::sample(
                    field,
                    lo,
                    fine_h,
                    [i64::from(fine_cells), 2 * jy, 2 * jz],
                    1,
                    2,
                    width,
                );
                cell.emit(field, 0, &mut transition);
            }
        }
    }

    // What is still open in the seam plane, counted on the assembled mesh.
    let mut assembled = MeshBuffer::<f32>::new();
    assembled
        .append(&fine)
        .expect("the meshes fit the u32 index space");
    assembled
        .append(&coarse)
        .expect("the meshes fit the u32 index space");
    assembled
        .append(&transition)
        .expect("the meshes fit the u32 index space");
    isomesh::weld::Welder::<f32>::new()
        .weld(&mut assembled, fine_h * 1e-5)
        .ok()?;
    let cfg = ValidateConfig::from_cell_size(f64::from(fine_h)).ok()?;
    let (_report, features) = validate_features(&assembled.positions, &assembled.indices, &cfg);
    let planes = [seam_x, seam_x + width];
    let gaps = features
        .boundary_edges
        .iter()
        .filter(|[a, b]| {
            let (pa, pb) = (
                assembled.positions[*a as usize][0],
                assembled.positions[*b as usize][0],
            );
            planes.iter().any(|plane| {
                (pa - plane).abs() < fine_h * 1e-4 && (pb - plane).abs() < fine_h * 1e-4
            })
        })
        .count();

    Some(Built {
        extent,
        fine,
        coarse,
        transition,
        gaps,
        fine_h,
        width,
    })
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    demo: Res<Demo>,
    mut stats: ResMut<DemoStats>,
    mut meshes: ResMut<Assets<Mesh>>,
    materials: Res<Materials>,
    mut commands: Commands,
    mut query: Query<(&mut Mesh3d, &mut Visibility, &Part)>,
    mut camera: Query<&mut OrbitCamera>,
    mut flags: ResMut<ViewFlags>,
    mut last: Local<Option<(usize, u32, bool)>>,
) {
    let key = (demo.field, demo.coarse_cells, demo.transitions);
    if *last == Some(key) && !flags.remesh_requested {
        return;
    }
    *last = Some(key);
    flags.remesh_requested = false;

    let built = match demo.field {
        0 => build(
            &Sphere::<f32>::canonical(),
            demo.coarse_cells,
            demo.transitions,
        ),
        1 => build(
            &Torus::<f32>::canonical(),
            demo.coarse_cells,
            demo.transitions,
        ),
        2 => build(
            &csg_difference::<f32>(),
            demo.coarse_cells,
            demo.transitions,
        ),
        _ => build(&capped_gyroid::<f32>(), demo.coarse_cells, demo.transitions),
    };
    let Some(built) = built else {
        return;
    };

    for mut orbit in &mut camera {
        orbit.radius = built.extent * 1.55;
    }

    stats.title = format!(
        "E-107  transvoxel seams   field {} ({})   [T] transitions {}   [1-4] field, [ ] resolution",
        demo.field + 1,
        FIELDS[demo.field],
        if demo.transitions { "ON" } else { "OFF" }
    );
    stats.vertices =
        built.fine.vertex_count() + built.coarse.vertex_count() + built.transition.vertex_count();
    stats.triangles = built.fine.triangle_count()
        + built.coarse.triangle_count()
        + built.transition.triangle_count();
    stats.extra = vec![
        format!(
            "fine h {:.4}   coarse h {:.4}   transition width {:.4}",
            built.fine_h,
            built.fine_h * 2.0,
            built.width
        ),
        String::new(),
        format!(
            "{:>6} unmatched boundary edges in the seam plane",
            built.gaps
        ),
        if built.gaps == 0 {
            "       the two resolutions meet".into()
        } else {
            "       ^^ the crack. press T.".into()
        },
        String::new(),
        format!(
            "{:>6} transition triangles (orange)",
            built.transition.triangle_count()
        ),
        String::new(),
        "counted in the seam plane only: both blocks are legitimately open".into(),
        "at their outer borders, so a global boundary count says nothing.".into(),
    ];

    let parts = [
        (Part::Fine, &built.fine, materials.fine.clone()),
        (Part::Coarse, &built.coarse, materials.coarse.clone()),
        (
            Part::Transition,
            &built.transition,
            materials.transition.clone(),
        ),
    ];

    // An **empty** buffer is never handed to Bevy. With the toggle off there are
    // no transition triangles at all, and uploading a zero-vertex mesh makes
    // `bevy_render`'s slab allocator report a use-after-free and render nothing
    // -- the whole frame, not just that entity. So an empty part keeps whatever
    // mesh it had and is hidden instead.
    if query.is_empty() {
        for (part, buffer, material) in parts {
            if buffer.is_empty() {
                continue;
            }
            commands.spawn((
                Mesh3d(meshes.add(to_bevy_mesh(buffer))),
                MeshMaterial3d(material),
                Transform::default(),
                DemoMesh,
                part,
            ));
        }
        return;
    }

    let rebuilt: Vec<(Part, Option<Handle<Mesh>>)> = parts
        .into_iter()
        .map(|(part, buffer, _)| {
            let handle = (!buffer.is_empty()).then(|| meshes.add(to_bevy_mesh(buffer)));
            (part, handle)
        })
        .collect();
    for (mut mesh, mut visibility, part) in &mut query {
        let Some((_, handle)) = rebuilt.iter().find(|(p, _)| p == part) else {
            continue;
        };
        match handle {
            Some(handle) => {
                mesh.0 = handle.clone();
                *visibility = Visibility::Inherited;
            }
            None => *visibility = Visibility::Hidden,
        }
    }
}
