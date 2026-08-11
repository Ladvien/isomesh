//! E-104 — the corner, side by side.
//!
//! Surface Nets (tan, left) against Dual Contouring (grey, right) on the same
//! field, the same grid and the same crossings. The only difference between the
//! two meshes on screen is one function: where a cell's vertex goes.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example dual_contouring_cube --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 20-50x slower.
//!
//! `1` and `2` switch field, `[` `]` resolution, `C` toggles the cell clamp,
//! `W` wireframe, `F12` screenshot.
//!
//! # The resolution is not a free parameter
//!
//! **`box_exact` is exactly zero across its entire boundary.** `f(1,0,0)`,
//! `f(1,1,0)` and `f(1,1,1)` are all `+0`, and the crate's convention is that
//! zero is *outside* — so a grid plane lying on a box face classifies wholly
//! outside and the sign change moves a cell inward. On such a grid this example
//! would be measuring the zero-classification rule rather than either algorithm.
//!
//! E-103 measured what that does: on an aligned grid Marching Cubes lands
//! *further* from the corner than Surface Nets, which is the opposite of the
//! intuition that edge-placed vertices should reach corners better.
//!
//! Over the ±2 domain a grid is aligned exactly when `n − 1` is a multiple of 4.
//! So this example **steps resolution by 2 from an odd base and skips the
//! aligned ones**, and the HUD says which regime you are in. 27³ is the default
//! and is not aligned; 25³ and 33³ are.
//!
//! # What the numbers underneath mean
//!
//! The HUD reports each side's distance from the true corner `(1,1,1)` to its
//! nearest vertex, in cells, recomputed on every re-mesh. Measured at 27³:
//! Surface Nets **0.58 cells**, Dual Contouring **0.01**.
//!
//! That gap is the entire case for dual contouring, and it is why this is the
//! image the README leads with.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::dc::{Clamp, DualContouring};
use isomesh::fields::{BoxExact, ReferenceField, csg_difference};
use isomesh::sn::SurfaceNets;
use isomesh::{RuntimeShape3, Sdf};

/// Samples per axis.
#[derive(Resource)]
struct Resolution(u32);

const MIN_SAMPLES: u32 = 11;
const MAX_SAMPLES: u32 = 81;

/// Whether the solved vertex is confined to its cell. A-009 measured this as
/// free — no sharpness cost — so it is on, and the toggle is here to show that
/// rather than to offer a choice.
#[derive(Resource)]
struct ClampOn(bool);

/// Which side of the screen a mesh sits on.
#[derive(Component, Clone, Copy, PartialEq, Eq)]
enum Side {
    SurfaceNets,
    DualContouring,
}

#[derive(Resource)]
struct Materials {
    surface_nets: Handle<StandardMaterial>,
    dual_contouring: Handle<StandardMaterial>,
}

const FIELDS: [&str; 2] = ["box_exact", "csg_difference"];

/// Which corner each field is measured at, and why it differs.
///
/// `box_exact` keeps all eight, so `(1,1,1)` is as good as any. **`csg_difference`
/// does not**: it is that box minus a radius-`0.75` sphere centred on
/// `(0.6, 0.6, 0.6)`, and `|(1,1,1) − (0.6,0.6,0.6)| = 0.693 < 0.75`, so that
/// corner is carved away entirely. Measuring it there compares two meshes'
/// distance to a point that is not on either of them, which reads as the two
/// algorithms being nearly equal when the picture plainly shows otherwise.
///
/// `(−1,−1,−1)` is the far corner, well outside the subtracted sphere, and it
/// survives.
fn measured_corner(field: usize) -> [f32; 3] {
    if field == 0 {
        [1.0, 1.0, 1.0]
    } else {
        [-1.0, -1.0, -1.0]
    }
}

/// A grid is aligned to `box_exact`'s faces when `n - 1` is a multiple of 4 over
/// the ±2 domain, which makes the comparison measure the sign convention instead
/// of the algorithm.
fn is_grid_aligned(samples: u32) -> bool {
    (samples - 1) % 4 == 0
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-104 dual contouring".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Resolution(27))
        .insert_resource(ClampOn(true))
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
        // Two constraints pulling against each other, and both matter.
        //
        // The pair is offset along world `x`, so `x` has to stay roughly
        // screen-horizontal or the two meshes drift diagonally apart and stop
        // reading as a comparison. That wants `yaw = π/2`, looking down `-z`.
        //
        // But a corner is only visibly a corner when two of the faces meeting at
        // it are in view; face-on shows a silhouette, which is the one angle at
        // which both algorithms look identical. That wants an oblique view.
        //
        // So: `π/2` plus a small rotation, and enough pitch to bring the top
        // face in. Three faces visible, pair still level.
        orbit.yaw = std::f32::consts::FRAC_PI_2 + 0.42;
        orbit.pitch = 0.42;
        orbit.radius = 7.2;
        // Look slightly up and to the left of the pair, which pushes it down and
        // right on screen and out from under the HUD. The readout is part of the
        // evidence here -- 0.577 against 0.006 cells -- so it must not sit on top
        // of the thing it is describing.
        orbit.focus = Vec3::new(-0.55, 0.55, 0.0);
    }
    commands.insert_resource(Materials {
        surface_nets: materials.add(StandardMaterial {
            base_color: Color::srgb(0.86, 0.74, 0.55),
            perceptual_roughness: 0.45,
            ..default()
        }),
        dual_contouring: materials.add(StandardMaterial {
            base_color: Color::srgb(0.78, 0.80, 0.86),
            perceptual_roughness: 0.45,
            ..default()
        }),
    });
}

fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut resolution: ResMut<Resolution>,
    mut clamp: ResMut<ClampOn>,
) {
    if capture.is_active() {
        // Sweep in step with captured frames so a recording is reproducible,
        // and step by 4 from an odd base so every frame stays off the aligned
        // grid — a GIF that wandered onto an aligned resolution would show the
        // comparison inverting for a reason no viewer could guess.
        const LOW: u32 = 15;
        const HIGH: u32 = 63;
        let steps = (HIGH - LOW) / 4 + 1;
        let phase = capture.taken % (steps * 2);
        let step = if phase < steps {
            phase
        } else {
            steps * 2 - phase - 1
        };
        resolution.0 = LOW + step * 4;
        return;
    }

    // Step by 2 and skip the aligned resolutions, so the comparison stays a
    // comparison. See the module docs.
    if keys.just_pressed(KeyCode::BracketRight) {
        let mut next = (resolution.0 + 2).min(MAX_SAMPLES);
        if is_grid_aligned(next) {
            next = (next + 2).min(MAX_SAMPLES);
        }
        resolution.0 = next;
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        let mut next = resolution.0.saturating_sub(2).max(MIN_SAMPLES);
        if is_grid_aligned(next) {
            next = next.saturating_sub(2).max(MIN_SAMPLES);
        }
        resolution.0 = next;
    }
    if keys.just_pressed(KeyCode::KeyC) {
        clamp.0 = !clamp.0;
    }
}

/// Distance from the true corner to the nearest vertex, in cells.
fn corner_gap(positions: &[[f32; 3]], corner: [f32; 3], cell_size: f32) -> f32 {
    positions
        .iter()
        .map(|p| {
            ((p[0] - corner[0]).powi(2) + (p[1] - corner[1]).powi(2) + (p[2] - corner[2]).powi(2))
                .sqrt()
        })
        .fold(f32::INFINITY, f32::min)
        / cell_size
}

struct Built {
    builder: MeshBuilder,
    gap: f32,
    extract_ms: f64,
}

fn build<F: Sdf<Scalar = f32> + ReferenceField>(
    field: &F,
    samples: u32,
    side: Side,
    clamp: bool,
    corner: [f32; 3],
) -> Option<(Built, f32, [f32; 3], [f32; 3])> {
    let (min, max) = field.domain();
    let cell_size = (max[0] - min[0]) / (samples - 1) as f32;
    let shape = match RuntimeShape3::new([samples; 3]) {
        Ok(shape) => shape,
        Err(error) => {
            error!("grid {samples}^3 rejected: {error}");
            return None;
        }
    };

    let mut builder = MeshBuilder::new();
    let started = Instant::now();
    let extracted = match side {
        Side::SurfaceNets => {
            SurfaceNets::<f32>::new().extract(field, &shape, min, cell_size, &mut builder)
        }
        Side::DualContouring => {
            let mut dc = DualContouring::<f32>::new();
            dc.set_clamp(if clamp { Clamp::ToCell } else { Clamp::None });
            dc.extract(field, &shape, min, cell_size, &mut builder)
        }
    };
    if let Err(error) = extracted {
        error!("extraction failed at {samples}^3: {error}");
        return None;
    }
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;
    let gap = corner_gap(builder.positions(), corner, cell_size);

    Some((
        Built {
            builder,
            gap,
            extract_ms,
        },
        cell_size,
        min,
        max,
    ))
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    resolution: Res<Resolution>,
    clamp: Res<ClampOn>,
    flags: Res<ViewFlags>,
    mut stats: ResMut<DemoStats>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut query: Query<(&mut Mesh3d, &mut Transform, &Side)>,
    mut commands: Commands,
    materials: Res<Materials>,
    mut domain: Query<&mut DemoDomain>,
    mut last: Local<Option<(u32, usize, bool)>>,
) {
    let field_index = flags.field.min(FIELDS.len() - 1);
    let key = (resolution.0, field_index, clamp.0);
    if *last == Some(key) && !flags.remesh_requested {
        return;
    }
    *last = Some(key);

    let samples = resolution.0;
    let corner = measured_corner(field_index);
    let built = if field_index == 0 {
        let f = BoxExact::<f32>::canonical();
        (
            build(&f, samples, Side::SurfaceNets, clamp.0, corner),
            build(&f, samples, Side::DualContouring, clamp.0, corner),
        )
    } else {
        let f = csg_difference::<f32>();
        (
            build(&f, samples, Side::SurfaceNets, clamp.0, corner),
            build(&f, samples, Side::DualContouring, clamp.0, corner),
        )
    };
    let (Some((sn, cell_size, min, max)), Some((dc, _, _, _))) = built else {
        return;
    };

    for mut d in &mut domain {
        d.min = Vec3::from(min);
        d.max = Vec3::from(max);
    }

    let aligned = is_grid_aligned(samples);
    stats.title = format!(
        "E-104  surface nets (tan)  vs  dual contouring (grey)   field: {}   [1-2] field, [ ] resolution, C clamp",
        FIELDS[field_index]
    );
    stats.vertices = sn.builder.vertex_count() + dc.builder.vertex_count();
    stats.triangles = sn.builder.triangle_count() + dc.builder.triangle_count();
    stats.extract_ms = sn.extract_ms + dc.extract_ms;
    stats.extra = vec![
        format!("{:>9} samples/axis   h = {cell_size:.4}", samples),
        format!(
            "          {}",
            if aligned {
                "!! GRID-ALIGNED - measuring the sign rule, not the algorithm"
            } else {
                "not grid-aligned, so this measures the algorithm"
            }
        ),
        String::new(),
        format!(
            "          distance from the true corner ({}, {}, {}), in cells",
            corner[0], corner[1], corner[2]
        ),
        format!("{:>9.3} surface nets", sn.gap),
        format!("{:>9.3} dual contouring", dc.gap),
        format!("{:>9.1}x closer", sn.gap / dc.gap.max(1e-6)),
        String::new(),
        format!(
            "          cell clamp {} (C) - A-009 measured it free",
            if clamp.0 { "ON" } else { "OFF" }
        ),
        format!(
            "{:>9} triangles each side (identical: same topology)",
            dc.builder.triangle_count()
        ),
    ];

    // Two meshes, offset either side of the origin, so neither is nearer the
    // camera than the other.
    let width = max[0] - min[0];
    // Close enough that the eye compares them rather than scanning between them.
    let offset = width * 0.34;
    let sides = [
        (Side::SurfaceNets, sn.builder, -offset),
        (Side::DualContouring, dc.builder, offset),
    ];

    if query.is_empty() {
        for (side, builder, x) in sides {
            let material = match side {
                Side::SurfaceNets => materials.surface_nets.clone(),
                Side::DualContouring => materials.dual_contouring.clone(),
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
        let mut handles = [None, None];
        for (side, builder, _) in sides {
            let slot = match side {
                Side::SurfaceNets => 0,
                Side::DualContouring => 1,
            };
            handles[slot] = Some(meshes.add(builder.into_mesh()));
        }
        for (mut mesh, mut transform, side) in &mut query {
            let slot = match side {
                Side::SurfaceNets => 0,
                Side::DualContouring => 1,
            };
            if let Some(handle) = handles[slot].clone() {
                mesh.0 = handle;
            }
            transform.translation.x = if *side == Side::SurfaceNets {
                -offset
            } else {
                offset
            };
        }
    }
}
