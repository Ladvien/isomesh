//! E-216 — **building** a field, not meshing one.
//!
//! ```bash
//! cargo run --example sdf_authoring --release
//! ```
//!
//! **Always `--release`.** A debug build meshes 37–62× slower (FINDINGS M-152)
//! and the blend sweep will feel broken.
//!
//! # The gap this fills
//!
//! Every other example here starts from a field that already exists — a
//! reference field, or one a brush has edited. So a reader learns that isomesh
//! *meshes* fields and never learns how to **write** the field a game would
//! ship. This one has no meshing content at all: the extractor is the default,
//! the resolution is fixed, and the only thing that changes is the expression.
//!
//! # What is on screen
//!
//! Left, the four primitives, each labelled with the expression that produced
//! it. Right, one asset assembled from them — a mushroom, because it needs a
//! union, a smooth union and a difference to look like anything, which is three
//! of the four operators doing visible work:
//!
//! ```text
//! stem   = Capsule { a: (0,-0.75,0), b: (0,0.15,0), radius: 0.18 }
//! cap    = Sphere  { center: (0,0.34,0), radius: 0.52 }
//! flat   = BoxExact half-extent (1, 0.34, 1) centred at (0,0.62,0)
//! gills  = Torus   { major: 0.34, minor: 0.07 } at (0,0.16,0)
//!
//! mushroom = SmoothUnion { a: stem, b: Difference { a: cap, b: flat }, k }
//!            unioned with gills
//! ```
//!
//! # The knob
//!
//! `[` and `]` sweep the smooth-union blend radius `k` and re-mesh. **That is
//! the parameter a level designer reaches for** and nothing else in this
//! repository puts it on screen: at `k = 0` the stem meets the cap in a crease,
//! and by `k = 0.25` it is a fillet. The HUD reports `k`, the triangle count and
//! the extraction time, so the cost of the blend is visible rather than
//! asserted.
//!
//! `1`–`4` isolate a single primitive; `0` returns to the assembly.

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use isomesh::brush::Capsule;
use isomesh::fields::{BoxExact, Difference, SmoothUnion, Sphere, Torus, Union};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{RuntimeShape3, Sdf};

mod common;
use common::{CommonPlugin, DemoStats};

/// Samples per axis for every field on screen. Fixed, because this example is
/// about the expression and not about resolution.
const SAMPLES: u32 = 65;

/// The domain every field is meshed over.
const HALF_EXTENT: f32 = 1.25;

/// Where each gallery primitive sits, and where the assembly sits.
///
/// Hand-placed rather than computed from a spacing constant, because the HUD
/// occupies the top-left and a row centred on the origin runs straight through
/// it — which the first recording showed and no amount of even spacing fixes.
const GALLERY_X: [f32; 4] = [-3.4, -1.9, -0.4, 1.1];

/// Gallery row height. Below the assembly, and below the text.
const GALLERY_Y: f32 = -1.45;

/// Where the finished asset sits.
const ASSEMBLY: [f32; 2] = [3.3, 0.15];

/// What the viewer is currently looking at.
#[derive(Resource)]
struct Authoring {
    /// Smooth-union blend radius, in world units.
    blend: f32,
    /// `0` for the assembly, `1..=4` for one primitive on its own.
    isolate: usize,
    dirty: bool,
}

impl Default for Authoring {
    fn default() -> Self {
        Self {
            blend: 0.12,
            isolate: 0,
            dirty: true,
        }
    }
}

/// The stem: a capsule from below the origin to just inside the cap.
fn stem() -> Capsule<f32> {
    Capsule {
        a: [0.0, -0.75, 0.0],
        b: [0.0, 0.15, 0.0],
        radius: 0.18,
    }
}

/// The cap, before it is flattened.
fn cap_sphere() -> Sphere<f32> {
    Sphere {
        center: [0.0, 0.34, 0.0],
        radius: 0.52,
    }
}

/// The half-space that flattens the cap's top, as a box far larger than the
/// domain in `x` and `z`.
fn cap_cutter() -> BoxExact<f32> {
    BoxExact {
        center: [0.0, 0.62, 0.0],
        half_extents: [1.0, 0.34, 1.0],
    }
}

/// A plain box for the gallery.
///
/// **Not the same box as [`cap_cutter`]**, deliberately. The cutter is a
/// half-space — a slab far wider than the domain, which is what makes it flatten
/// the cap rather than dent it — and putting that on a shelf of primitives shows
/// a wall, not a box.
fn gallery_box() -> BoxExact<f32> {
    BoxExact {
        center: [0.0, 0.0, 0.0],
        half_extents: [0.42, 0.42, 0.42],
    }
}

/// The ring under the cap.
fn gills() -> Torus<f32> {
    Torus {
        center: [0.0, 0.16, 0.0],
        major: 0.34,
        minor: 0.07,
    }
}

/// The whole asset, as one expression.
///
/// Read outward: the cap is a sphere with its top cut off, the stem is
/// *smoothly* unioned to it so the join is a fillet rather than a crease, and
/// the gills are unioned hard because a ring should stay a ring.
fn mushroom(blend: f32) -> impl Sdf<Scalar = f32> {
    Union {
        a: SmoothUnion {
            a: stem(),
            b: Difference {
                a: cap_sphere(),
                b: cap_cutter(),
            },
            k: blend,
        },
        b: gills(),
    }
}

/// Mesh one field into a Bevy mesh, and say what it cost.
fn build(field: &impl Sdf<Scalar = f32>) -> (Mesh, usize, f64) {
    let shape = RuntimeShape3::new([SAMPLES; 3]).expect("valid shape");
    let h = 2.0 * HALF_EXTENT / (SAMPLES - 1) as f32;
    let mut builder = MeshBuilder::new();
    let started = std::time::Instant::now();
    MarchingCubes::<f32>::new()
        .extract(field, &shape, [-HALF_EXTENT; 3], h, &mut builder)
        .expect("extraction");
    let ms = started.elapsed().as_secs_f64() * 1e3;
    let triangles = builder.triangle_count();
    (builder.into_mesh(), triangles, ms)
}

#[derive(Component)]
struct Piece;

/// One gallery entry: its label, and how to mesh it.
///
/// Named because the tuple is genuinely hard to read inline, and because the
/// boxed closure is the point — each primitive is a different concrete `Sdf`
/// type, so they cannot sit in an array without erasing to a common return.
type GalleryEntry = (&'static str, Box<dyn Fn() -> (Mesh, usize, f64)>);

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-216 building a field, not meshing one".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_resource::<Authoring>()
        .add_systems(Startup, aim_camera)
        .add_systems(Update, (controls, rebuild, report).chain())
        .run();
}

fn aim_camera(mut cameras: Query<&mut Transform, With<Camera3d>>) {
    for mut transform in &mut cameras {
        *transform =
            Transform::from_xyz(0.6, 1.2, 8.2).looking_at(Vec3::new(0.2, -0.55, 0.0), Vec3::Y);
    }
}

fn controls(keys: Res<ButtonInput<KeyCode>>, mut state: ResMut<Authoring>) {
    let mut changed = false;
    if keys.just_pressed(KeyCode::BracketRight) {
        state.blend = (state.blend + 0.02).min(0.45);
        changed = true;
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        state.blend = (state.blend - 0.02).max(0.0);
        changed = true;
    }
    for (key, which) in [
        (KeyCode::Digit0, 0),
        (KeyCode::Digit1, 1),
        (KeyCode::Digit2, 2),
        (KeyCode::Digit3, 3),
        (KeyCode::Digit4, 4),
    ] {
        if keys.just_pressed(key) {
            state.isolate = which;
            changed = true;
        }
    }
    if changed {
        state.dirty = true;
    }
}

#[allow(clippy::too_many_lines)]
fn rebuild(
    mut commands: Commands,
    mut state: ResMut<Authoring>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut stats: ResMut<DemoStats>,
    existing: Query<Entity, With<Piece>>,
) {
    if !state.dirty {
        return;
    }
    state.dirty = false;
    for entity in &existing {
        commands.entity(entity).despawn();
    }

    let stone = materials.add(StandardMaterial {
        base_color: Color::srgb(0.78, 0.74, 0.68),
        perceptual_roughness: 0.85,
        ..default()
    });
    let accent = materials.add(StandardMaterial {
        base_color: Color::srgb(0.85, 0.62, 0.42),
        perceptual_roughness: 0.7,
        ..default()
    });

    let mut total_triangles = 0usize;
    let mut total_ms = 0.0f64;
    let mut spawn =
        |mesh: Mesh, tris: usize, ms: f64, at: [f32; 2], material: Handle<StandardMaterial>| {
            total_triangles += tris;
            total_ms += ms;
            commands.spawn((
                Mesh3d(meshes.add(mesh)),
                MeshMaterial3d(material),
                Transform::from_xyz(at[0], at[1], 0.0),
                Piece,
            ));
        };

    // The gallery: one primitive per column, so the expression and the shape sit
    // next to each other.
    let gallery: [GalleryEntry; 4] = [
        ("Capsule", Box::new(|| build(&stem()))),
        ("Sphere", Box::new(|| build(&cap_sphere()))),
        ("BoxExact", Box::new(|| build(&gallery_box()))),
        ("Torus", Box::new(|| build(&gills()))),
    ];

    if state.isolate == 0 {
        for (index, (_, make)) in gallery.iter().enumerate() {
            let (mesh, tris, ms) = make();
            spawn(mesh, tris, ms, [GALLERY_X[index], GALLERY_Y], stone.clone());
        }
        let (mesh, tris, ms) = build(&mushroom(state.blend));
        spawn(mesh, tris, ms, ASSEMBLY, accent);
    } else {
        let (_, make) = &gallery[state.isolate - 1];
        let (mesh, tris, ms) = make();
        spawn(mesh, tris, ms, [0.0, 0.0], stone);
    }

    let what = match state.isolate {
        0 => "assembly + gallery",
        n => gallery[n - 1].0,
    };
    stats.title = format!("E-216  authoring — {what}");
    stats.triangles = total_triangles;
    stats.extract_ms = total_ms;
    stats.extra = vec![
        format!(
            "blend radius k = {:.2}          [ and ] to sweep",
            state.blend
        ),
        String::new(),
        "mushroom = Union {".into(),
        "    a: SmoothUnion { stem, Difference { cap, flat }, k },".into(),
        "    b: gills,".into(),
        "}".into(),
        String::new(),
        "k = 0 leaves a crease where the stem meets the cap;".into(),
        "k = 0.25 is a fillet. That is the whole knob.".into(),
        String::new(),
        "[1-4] isolate a primitive   [0] the assembly".into(),
    ];
}

fn report(state: Res<Authoring>, mut last: Local<Option<f32>>) {
    if *last != Some(state.blend) {
        *last = Some(state.blend);
        info!("blend radius k = {:.2}", state.blend);
    }
}
