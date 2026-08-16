//! **What a position-only weld costs you, and the one line that buys it back.**
//!
//! ```bash
//! cargo run --example weld_creases --release
//! ```
//!
//! Ticket: B-014. Two cubes, same input, same tolerance. The left one is welded
//! on position alone and its eight corners swallow all 24 vertices — the faces
//! now share normals, so the flat shading that made it read as a cube is gone and
//! it renders like a rounded blob. The right one is welded with a key built from
//! the vertex normals, keeps all 24, and still looks like a cube.
//!
//! **The point is that both welds are correct.** Neither is a bug. The left one
//! did exactly what a position weld means; it is just that "merge coincident
//! vertices" and "preserve the crease" are different requests, and only the
//! caller knows which one it wants.
//!
//! # Why the key is a quantum and not a smoothing angle
//!
//! The conventional test is "merge if the normals are within 30°", and this crate
//! deliberately does not offer it. **An angle threshold is not transitive**: `a`
//! within 30° of `b` and `b` within 30° of `c` does not put `a` within 30° of
//! `c`. So it is not an equivalence relation, and applied to a `k`-way
//! coincidence it merges some members and refuses others — leaving the leftover
//! representative a bowtie. This repository measured that exact shape adding **up
//! to 791 non-manifold vertices** (E×4 in `FINDINGS.md`).
//!
//! Quantising to a lattice *is* transitive, so a class always splits into
//! complete sub-classes. Its failure mode is a **missed merge** at a bucket
//! boundary — a visible seam, harmless topologically — which is the right failure
//! to prefer over a manufactured bowtie.
//!
//! # What the readout is telling you
//!
//! The vertex counts, and one number that surprises people: the key-welded mesh
//! has **24 boundary edges** where the position-welded one has none. That is not
//! damage. Keeping a cube's six faces apart necessarily opens the edges between
//! them, so the split mesh is a **surface**, not a solid — and its open edges are
//! a recorded number rather than a failure. `isomesh::validate::SurfaceGate` is
//! how you say which of those two things you meant (M-305, ✗22).

use bevy::prelude::*;
use bevy_isomesh::{WeldKeyConfig, from_bevy_mesh, weld_keys};
use isomesh::MeshBuffer;
use isomesh::validate::{SurfaceGate, ValidateConfig, validate_indexed};
use isomesh::weld::{Welder, epsilon_for};

/// Edge length of each cube.
const SIDE: f32 = 1.6;

/// How far apart the two are placed.
const GAP: f32 = 1.3;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_systems(Startup, setup)
        .run();
}

/// One cube, welded with the given keys, as a Bevy mesh plus its readout.
///
/// `keys` empty means the unconditional position weld.
fn welded(keys_from: Option<WeldKeyConfig>) -> (Mesh, String) {
    let source = Mesh::from(bevy_math::primitives::Cuboid::new(SIDE, SIDE, SIDE));
    let (positions, indices) = match from_bevy_mesh(&source) {
        Ok(pair) => pair,
        // A cuboid is a triangle list with positions and indices, so this arm is
        // unreachable -- but an example that unwraps teaches unwrapping.
        Err(err) => return (source, format!("could not read the cuboid: {err}")),
    };
    let normals: Vec<[f32; 3]> = source
        .attribute(Mesh::ATTRIBUTE_NORMAL)
        .and_then(|values| values.as_float3())
        .map(<[[f32; 3]]>::to_vec)
        .unwrap_or_default();

    let mut buffer = MeshBuffer::<f32>::new();
    buffer.positions.extend_from_slice(&positions);
    buffer.normals.extend_from_slice(&normals);
    buffer.indices.extend_from_slice(&indices);
    let before = buffer.positions.len();

    let keys = keys_from
        .map(|config| weld_keys(&source, config))
        .unwrap_or_default();
    let epsilon = epsilon_for(SIDE);
    if let Err(err) = Welder::default().weld_split_by(&mut buffer, epsilon, &keys) {
        return (source, format!("weld failed: {err}"));
    }

    // The gate this artefact earns. A position-welded cube is a closed solid; a
    // key-split one is a surface, and asserting `Closed` on it would report
    // correct output as broken.
    let gate = if keys.is_empty() {
        SurfaceGate::Closed
    } else {
        SurfaceGate::Manifold
    };
    let report = match ValidateConfig::from_cell_size(f64::from(SIDE)) {
        Ok(cfg) => {
            let p: Vec<[f64; 3]> = buffer
                .positions
                .iter()
                .map(|v| [f64::from(v[0]), f64::from(v[1]), f64::from(v[2])])
                .collect();
            Some(validate_indexed(&p, &buffer.indices, &cfg))
        }
        Err(_) => None,
    };

    let readout = match report {
        Some(r) => format!(
            "{before} -> {} vertices\n{} boundary edges\n{:?}: {}",
            buffer.positions.len(),
            r.boundary_edges,
            gate,
            if r.satisfies(gate) { "met" } else { "NOT met" }
        ),
        None => format!("{before} -> {} vertices", buffer.positions.len()),
    };

    (bevy_isomesh::to_bevy_mesh(&buffer), readout)
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.75, 0.76, 0.80),
        perceptual_roughness: 0.65,
        ..default()
    });

    let (flat, flat_text) = welded(None);
    let (creased, creased_text) = welded(Some(WeldKeyConfig::default()));

    for (mesh, x) in [(flat, -GAP), (creased, GAP)] {
        commands.spawn((
            Mesh3d(meshes.add(mesh)),
            MeshMaterial3d(material.clone()),
            Transform::from_xyz(x, 0.0, 0.0),
        ));
    }

    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(0.0, 2.2, 4.6).looking_at(Vec3::ZERO, Vec3::Y),
    ));
    commands.spawn((
        DirectionalLight {
            illuminance: 6_000.0,
            shadow_maps_enabled: false,
            ..default()
        },
        Transform::from_xyz(3.0, 5.0, 4.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));

    commands.spawn((
        Text::new(format!(
            "weld_creases  —  B-014\n\n\
             left   position only\n{flat_text}\n\n\
             right  split on the normal key\n{creased_text}\n\n\
             Both welds are correct. The left one did exactly what a position\n\
             weld means; only the caller knows which was wanted.\n\
             The 24 boundary edges on the right are the split seen from the\n\
             edge column, not damage — so the right mesh is a surface, not a solid."
        )),
        TextFont {
            font_size: bevy::text::FontSize::Px(14.0),
            ..default()
        },
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(12.0),
            left: Val::Px(12.0),
            ..default()
        },
    ));
}
