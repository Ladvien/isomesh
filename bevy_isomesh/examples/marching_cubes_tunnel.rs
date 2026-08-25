//! E-213 — the tunnel, meshed as a tunnel.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example marching_cubes_tunnel --release
//! ```
//!
//! `1` and `2` switch configuration, `H` toggles the inner hexagon, `C` toggles
//! the contours, `W` toggles wireframe-only.
//!
//! # What is on screen
//!
//! One cell, meshed twice. **Left** is the face rule alone — Marching Cubes 33's
//! asymptotic decider, which resolves every ambiguous *face* and then treats the
//! interior as if the faces settled it. **Right** is the same cell with the
//! interior rule as well. Left gives two separate discs; right gives one
//! cylinder passing through the cell, which is what the trilinear interpolant
//! actually does there.
//!
//! The third sibling of `marching_cubes_ambiguity` (an ambiguous *face*, E-102)
//! and `marching_cubes_interior` (the *decider* that detects an ambiguous cell,
//! E-116). Those two show the decision; this one shows the surface it decides
//! about, and is the first of the three where the mesh changes.
//!
//! # The inner hexagon is the whole construction
//!
//! Six points, drawn as a closed ring in gold. They are the trilinear
//! interpolant's **body saddles** — where two of the straight lines that lie on
//! the level set cross — and Grosso's construction is built entirely out of
//! them: each contour vertex is assigned to its nearest hexagon vertex, and each
//! contour edge is then closed with one, two or three triangles depending on how
//! far apart around the ring its two ends landed.
//!
//! That is also what buys manifoldness. Chernyaev's tunnel triangulation lays
//! part of the tunnel *on* an ambiguous face, so the two cells sharing that face
//! both claim it and the shared edge ends up with four triangles. Grosso's stays
//! strictly inside the cell — every triangle drawn on the right touches the cell
//! boundary only along a contour edge, which its neighbour shares exactly once.
//!
//! # The number that grades it
//!
//! A tunnel is a handle, and a handle costs a closed surface exactly two of its
//! Euler characteristic. So the claim this example is illustrating is arithmetic:
//! **χ falls by two per tunnel and by nothing else** (M-222). Everything else the
//! interior rule does — giving an ambiguous contour an interior vertex and
//! fanning from it — adds one vertex, three edges and two faces, and
//! `1 − 3 + 2 = 0`.
//!
//! The HUD prints both χ values and their difference. On a single cell the
//! surface is open, so χ is not `2 − 2g` and the *difference* is what to read.
//!
//! # Why the configurations are hand-written
//!
//! Because almost nothing reaches this. Interior ambiguity occurs in **0 of
//! 68,385** surface cells across the seven original reference fields (M-208) —
//! that is what `noise_cavity` was added for, and even there it is three or four
//! cells in a grid. A demo that sampled a field would be showing an empty cell
//! almost every time.
//!
//! So these eight corner values are searched rather than invented: swept for a
//! cell whose body saddles number six and whose contours form a tunnel, then
//! rounded to one decimal so the constant is readable.

mod common;

use bevy::color::palettes::css;
use bevy::prelude::*;
use common::{CommonPlugin, DemoStats, OrbitCamera};
use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::table::{AMBIGUOUS_FACES, EDGE_CORNERS};
use isomesh::marching_cubes::trilinear::{BodySaddles, Contours, Topology, local_crossing};
use isomesh::marching_cubes::{FaceAmbiguity, InteriorAmbiguity, MarchingCubes};
use isomesh::validate::{ValidateConfig, validate_indexed};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// One labelled cell.
struct Configuration {
    name: &'static str,
    /// The eight corner values, in `isomesh::cube`'s numbering — which is also
    /// Grosso's, a coincidence pinned by `grosso_corner_numbering_is_ours`.
    corners: [f64; 8],
    note: &'static str,
}

const CONFIGURATIONS: [Configuration; 2] = [
    Configuration {
        name: "tunnel — contours of six and three",
        corners: [-0.2, -0.9, 0.7, -0.3, 0.8, -0.4, -0.9, 0.6],
        note: "Corollary 6's shape: one contour of at most six, one of three",
    },
    Configuration {
        name: "tunnel — contours of four and three",
        corners: [0.2, -0.6, 0.2, -0.2, 0.4, 0.1, -0.9, 1.0],
        note: "the same topology through a smaller pair of contours",
    },
];

#[derive(Resource)]
struct Show {
    which: usize,
    hexagon: bool,
    contours: bool,
    wireframe: bool,
}

/// The trilinear interpolant of eight corner values over the unit cell.
///
/// Sampling this on a `2×2×2` grid reproduces the corner values exactly, which
/// is the point: the configurations above are cell data, not a field, and this
/// is the smallest honest way to hand them to an extractor.
struct Cell {
    values: [f64; 8],
}

impl Sdf for Cell {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let mut sum = 0.0;
        for (i, &v) in self.values.iter().enumerate() {
            let mut w = v;
            for (axis, s) in p.iter().enumerate() {
                let bit = (i >> axis) & 1;
                w *= if bit == 1 { *s } else { 1.0 - *s };
            }
            sum += w;
        }
        sum
    }
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-213 the tunnel, meshed as a tunnel".into(),
                // Web only, inert on native: bind to the 1280x720 canvas the
                // page supplies rather than letting Bevy append its own. The HUD
                // panels are laid out in pixels for that size, so the canvas is
                // fixed and CSS scales it -- `fit_canvas_to_parent` stays at its
                // `false` default for the same reason.
                canvas: Some("#isomesh-canvas".into()),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Show {
            which: 0,
            hexagon: true,
            contours: true,
            wireframe: false,
        })
        .add_systems(Startup, aim_camera)
        .add_systems(Update, (controls, fill_patches, draw, report).chain())
        .run();
}

fn aim_camera(mut camera: Query<(&mut Transform, &mut OrbitCamera)>) {
    for (mut transform, mut orbit) in &mut camera {
        orbit.focus = Vec3::new(0.0, 0.5, 0.5);
        orbit.radius = 4.2;
        *transform = Transform::from_translation(Vec3::new(2.4, 2.2, 3.4))
            .looking_at(Vec3::new(0.0, 0.5, 0.5), Vec3::Y);
    }
}

fn controls(keys: Res<ButtonInput<KeyCode>>, mut show: ResMut<Show>) {
    if keys.just_pressed(KeyCode::Digit1) {
        show.which = 0;
    }
    if keys.just_pressed(KeyCode::Digit2) {
        show.which = 1;
    }
    if keys.just_pressed(KeyCode::KeyH) {
        show.hexagon = !show.hexagon;
    }
    if keys.just_pressed(KeyCode::KeyC) {
        show.contours = !show.contours;
    }
    if keys.just_pressed(KeyCode::KeyW) {
        show.wireframe = !show.wireframe;
    }
}

/// Mesh the cell under one interior rule, as a `2×2×2` grid.
fn mesh_with(values: [f64; 8], interior: InteriorAmbiguity) -> MeshBuffer<f64> {
    let field = Cell { values };
    let shape = RuntimeShape3::new([2; 3]).expect("a 2^3 grid");
    let mut mc = MarchingCubes::<f64>::new();
    mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
    mc.set_interior_ambiguity(interior);
    let mut out = MeshBuffer::<f64>::new();
    // Refuses only on the case A-020 owns, which neither configuration reaches.
    mc.extract(&field, &shape, [0.0; 3], 1.0, &mut out)
        .expect("neither configuration is A-020's undefined case");
    out
}

/// Marks the two filled patches, so they can be replaced when the configuration
/// changes.
#[derive(Component)]
struct Patch;

/// Spawn the two surfaces as **solid** meshes, not only as gizmo outlines.
///
/// # Why this exists
///
/// The example drew every triangle as three lines and nothing else, with a
/// per-triangle normal tick added *"so the two sides read as surfaces rather
/// than as a cage of lines"* — which concedes the problem. Interactively it is a
/// good diagram: you can see through it to the hexagon and the contours. As a
/// picture it fails, because the claim being made is about **surfaces** — two
/// discs against one cylinder — and a cage of lines does not show a surface at
/// all. Recorded to a GIF it read as scaffolding.
///
/// So the patches are filled, and the gizmo overlay stays on top of them. The
/// mesh comes from the same `mesh_with` call the gizmos use, so the outline and
/// the fill cannot disagree about what was extracted.
fn fill_patches(
    mut commands: Commands,
    show: Res<Show>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    existing: Query<Entity, With<Patch>>,
    mut last: Local<Option<usize>>,
) {
    if *last == Some(show.which) {
        return;
    }
    *last = Some(show.which);
    for entity in &existing {
        commands.entity(entity).despawn();
    }
    let Some(config) = CONFIGURATIONS.get(show.which) else {
        return;
    };

    for (shift, interior, colour) in [
        (
            -1.3f32,
            InteriorAmbiguity::Ignore,
            Color::srgb(0.62, 0.66, 0.70),
        ),
        (
            0.9,
            InteriorAmbiguity::Trilinear,
            Color::srgb(0.20, 0.68, 0.92),
        ),
    ] {
        let buffer = mesh_with(config.corners, interior);
        if buffer.triangle_count() == 0 {
            continue;
        }
        let positions: Vec<[f32; 3]> = buffer
            .positions
            .iter()
            .map(|p| [p[0] as f32 + shift, p[1] as f32, p[2] as f32])
            .collect();
        let normals: Vec<[f32; 3]> = buffer
            .normals
            .iter()
            .map(|n| [n[0] as f32, n[1] as f32, n[2] as f32])
            .collect();
        let mut mesh = Mesh::new(
            bevy::render::mesh::PrimitiveTopology::TriangleList,
            bevy::asset::RenderAssetUsages::default(),
        );
        mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
        mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
        mesh.insert_indices(bevy::render::mesh::Indices::U32(buffer.indices.clone()));
        commands.spawn((
            Mesh3d(meshes.add(mesh)),
            MeshMaterial3d(materials.add(StandardMaterial {
                base_color: colour,
                // Both sides matter: a disc seen from behind is still the claim.
                cull_mode: None,
                perceptual_roughness: 0.75,
                ..default()
            })),
            Patch,
        ));
    }
}

fn to_vec3(p: [f64; 3], shift: f32) -> Vec3 {
    Vec3::new(p[0] as f32 + shift, p[1] as f32, p[2] as f32)
}

fn draw(show: Res<Show>, mut gizmos: Gizmos) {
    let Some(config) = CONFIGURATIONS.get(show.which) else {
        return;
    };
    let values = config.corners;

    for (shift, interior, tint) in [
        (-1.3f32, InteriorAmbiguity::Ignore, css::SLATE_GRAY),
        (0.9, InteriorAmbiguity::Trilinear, css::DEEP_SKY_BLUE),
    ] {
        // The cell, so the surface is read against something.
        for (a, b) in EDGE_CORNERS.iter().map(|[a, b]| (*a, *b)) {
            let p = |c: u8| {
                Vec3::new(
                    (c & 1) as f32 + shift,
                    ((c >> 1) & 1) as f32,
                    ((c >> 2) & 1) as f32,
                )
            };
            gizmos.line(p(a), p(b), css::DARK_SLATE_GRAY);
        }
        // Corner signs: filled for inside, hollow for outside.
        for (c, &v) in values.iter().enumerate() {
            let at = Vec3::new(
                (c & 1) as f32 + shift,
                ((c >> 1) & 1) as f32,
                ((c >> 2) & 1) as f32,
            );
            let colour = if v < 0.0 { css::ORANGE_RED } else { css::WHITE };
            gizmos.sphere(Isometry3d::from_translation(at), 0.045, colour);
        }

        let mesh = mesh_with(values, interior);
        for t in mesh.indices.chunks_exact(3) {
            let v = |i: u32| to_vec3(mesh.positions[i as usize], shift);
            let (a, b, c) = (v(t[0]), v(t[1]), v(t[2]));
            gizmos.line(a, b, tint);
            gizmos.line(b, c, tint);
            gizmos.line(c, a, tint);
            if !show.wireframe {
                // A short normal tick per triangle, so the two sides read as
                // surfaces rather than as a cage of lines.
                let centre = (a + b + c) / 3.0;
                let n = (b - a).cross(c - a).normalize_or_zero();
                gizmos.line(centre, centre + n * 0.06, tint.with_alpha(0.55));
            }
        }
    }

    // The construction itself, on the right-hand copy only.
    let saddles = BodySaddles::of(&values);
    if let Some(ring) = saddles.inner_hexagon().filter(|_| show.hexagon) {
        {
            for k in 0..6 {
                let a = to_vec3(ring[k], 0.9);
                let b = to_vec3(ring[(k + 1) % 6], 0.9);
                gizmos.line(a, b, css::GOLD);
                gizmos.sphere(Isometry3d::from_translation(a), 0.03, css::GOLD);
            }
        }
    }
    if show.contours {
        let mut case = 0u8;
        for (c, &v) in values.iter().enumerate() {
            if v < 0.0 {
                case |= 1 << c;
            }
        }
        let mask = joined_mask(&values, AMBIGUOUS_FACES[case as usize]);
        let contours = Contours::of(case, mask);
        for r in 0..contours.count() {
            let ring = contours.ring(r);
            let colour = if r == 0 {
                css::SPRING_GREEN
            } else {
                css::MAGENTA
            };
            for k in 0..ring.len() {
                let a = to_vec3(local_crossing(ring[k], &values), 0.9);
                let b = to_vec3(local_crossing(ring[(k + 1) % ring.len()], &values), 0.9);
                gizmos.line(a, b, colour);
            }
        }
    }
}

/// Both meshes' topology, and the difference that is the claim.
fn report(show: Res<Show>, mut stats: ResMut<DemoStats>) {
    let Some(config) = CONFIGURATIONS.get(show.which) else {
        return;
    };
    let values = config.corners;

    let cfg = ValidateConfig::from_cell_size(1.0).expect("a valid spacing");
    let chi_of = |interior| {
        let mesh = mesh_with(values, interior);
        let report = validate_indexed(&mesh.positions, &mesh.indices, &cfg);
        (
            report.euler_characteristic,
            report.components,
            mesh.indices.len() / 3,
        )
    };
    let (chi_face, parts_face, tris_face) = chi_of(InteriorAmbiguity::Ignore);
    let (chi_full, parts_full, tris_full) = chi_of(InteriorAmbiguity::Trilinear);

    let mut case = 0u8;
    for (c, &v) in values.iter().enumerate() {
        if v < 0.0 {
            case |= 1 << c;
        }
    }
    let saddles = BodySaddles::of(&values);
    let mask = joined_mask(&values, AMBIGUOUS_FACES[case as usize]);
    let contours = Contours::of(case, mask);
    let sizes: Vec<usize> = (0..contours.count())
        .map(|r| contours.ring(r).len())
        .collect();
    let topology = contours.topology(&saddles);

    stats.title = format!("E-213  the tunnel — {}", config.name);
    stats.extra = vec![
        format!("case {case:#010b}   face mask {mask:#08b}"),
        format!(
            "body saddles {} of 6   topology {topology:?}   contours {sizes:?}",
            saddles.inside_count()
        ),
        String::new(),
        format!(
            "left   face rule only : chi {chi_face:>4}   parts {parts_face}   {tris_face} tris"
        ),
        format!(
            "right  + interior rule: chi {chi_full:>4}   parts {parts_full}   {tris_full} tris"
        ),
        format!(
            "difference: {}   {}",
            chi_full - chi_face,
            if topology == Topology::Tunnel {
                "a tunnel is a handle, and a handle costs exactly two"
            } else {
                "no tunnel here, so chi must not move"
            }
        ),
        String::new(),
        config.note.to_string(),
        "interior ambiguity: 0 of 68,385 cells on the seven original fields (M-208)".to_string(),
        "1 2 configuration | H hexagon | C contours | W wireframe".to_string(),
    ];
}
