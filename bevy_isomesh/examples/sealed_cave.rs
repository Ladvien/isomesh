//! **Did I just seal the cave?** — connectivity across a chunk seam.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example sealed_cave --release
//! ```
//!
//! `F` plugs the tunnel, `G` re-opens it, `C` outlines the three chunks.
//!
//! Ticket: R-028. Two chambers in **different chunks**, joined by a tunnel
//! through the chunk between them. Plugging the tunnel is one edit inside one
//! chunk, and the question it answers — *are these two rooms still connected?* —
//! is one no single chunk can answer alone.
//!
//! # Why this example is not a meshing demo
//!
//! Thirty-odd examples here mesh a field. This one is about the query beside the
//! mesh. `isomesh::validate` can tell you a mesh is closed, manifold, correctly
//! wound and Hausdorff-close to the field, and **none of that tells you whether
//! the water in the left chamber can reach the right one.** That is a question
//! about the connected components of the air region, and it is asked after every
//! edit rather than once at build time.
//!
//! # What the readout is showing
//!
//! - **`components`** — air components across the whole world, seams included.
//!   It goes `1 → 2` on the plug and back on the dig.
//! - **`A↔B`** — [`AirWorld::connected`] between one sample in each chamber, in
//!   global sample coordinates. This is the query the mechanic is built on.
//! - **`visited`** — samples the replacement search touched, and the number
//!   R-028 exists to bound. **Compare it to `chunk`, not to `world`**: the search
//!   runs inside one [`Air`] and cannot exceed it, however many chunks are
//!   loaded (M-322).
//! - **`stitch`** — what the restitch cost. Its `nodes` are *components*, not
//!   samples, which is what makes rebuilding the whole global graph after every
//!   edit affordable.
//!
//! # The number worth watching, and what actually sets it
//!
//! Severing a passage is the shape that beats the unchunked structure: lockstep
//! search stops when all but one frontier exhausts, so a fill that splits two
//! pieces of *similar* size walks both. M-321 measured exactly this edit at
//! **1.1× a full rebuild** on one large unchunked grid, and M-322 at **0.97× a
//! chunk** rebuild when both halves lived inside the edited chunk.
//!
//! **This example reports about 138 — roughly 0.03× a chunk — and the gap is the
//! lesson.** The chambers are in chunks 0 and 2; the plug is in chunk 1. So the
//! pieces the *local* search separates are two short tunnel stubs, not two
//! caverns, and the split between the large components is resolved by the
//! **boundary graph**, whose nodes are components rather than samples.
//!
//! So the bound and the cost are different statements. **Chunking bounds the
//! search by the chunk** — that part is structural (M-322). **What it costs is
//! the edited chunk's share of the severed component**, which is geometry, and
//! here that share is a pair of stubs. Move the chambers inside one chunk and
//! the same edit approaches the bound instead.
//!
//! # The plug drives the mesh and the connectivity from one flag
//!
//! `Cave::plugged` is read by the extractor and by the sampling that feeds
//! [`AirWorld`], so the picture and the answer cannot disagree. A real engine
//! would apply one edit to one field and route it to both; this keeps them in
//! step by construction, which is the point rather than a shortcut.
//!
//! # Spacing is a power of two, deliberately
//!
//! `h = 0.125`. **M-32** measured that two chunks agree on their shared sample
//! plane bit-for-bit only at a power-of-two cell size. The seam stitch here does
//! not depend on that — it joins components by **integer sample identity**, not
//! by position — but the *mesh* does, and a visible crack would distract from
//! what this example is about.

mod common;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{CommonPlugin, DemoStats, OrbitCamera};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::connectivity::AirWorld;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, MeshSink, Sdf};

/// A [`MeshBuffer`] as a Bevy mesh.
fn to_mesh(buffer: &MeshBuffer<f32>) -> Mesh {
    let mut builder = MeshBuilder::new();
    for i in 0..buffer.positions.len() {
        let (Some(p), Some(n)) = (buffer.positions.get(i), buffer.normals.get(i)) else {
            continue;
        };
        builder.vertex(*p, *n);
    }
    for t in buffer.indices.chunks_exact(3) {
        let (Some(a), Some(b), Some(c)) = (t.first(), t.get(1), t.get(2)) else {
            continue;
        };
        builder.triangle(*a, *b, *c);
    }
    builder.into_mesh()
}

/// Cells per chunk. A chunk is `17³` samples.
const CELLS: u32 = 16;

/// Sample spacing. A power of two — see the module docs.
const H: f32 = 0.125;

/// Chunks along x. The chambers live in the first and last.
const CHUNKS: i32 = 3;

/// Chamber radius, world units.
const CHAMBER: f32 = 0.72;

/// Tunnel radius. Wide enough to be air at this spacing — about two samples.
const TUNNEL: f32 = 0.25;

/// Plug radius. Comfortably wider than the tunnel, so it severs rather than
/// pinches.
const PLUG: f32 = 0.45;

/// World-space centre of each chamber.
const CHAMBER_X: f32 = 2.0;

/// The cave: two chambers joined by a tunnel, optionally plugged.
///
/// Negative inside the rock, positive inside the cave — the crate's convention,
/// and what makes `value >= 0` mean *air* to [`AirWorld`].
#[derive(Clone, Copy)]
struct Cave {
    plugged: bool,
}

fn sphere(p: [f32; 3], c: [f32; 3], r: f32) -> f32 {
    let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - r
}

impl Sdf for Cave {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        // The cave as a solid: union of two chambers and the tunnel between.
        let a = sphere(p, [-CHAMBER_X, 0.0, 0.0], CHAMBER);
        let b = sphere(p, [CHAMBER_X, 0.0, 0.0], CHAMBER);
        // A capsule along x, which is a cylinder with rounded caps -- the caps
        // sit inside the chambers, so they are never seen.
        let x = p[0].clamp(-CHAMBER_X, CHAMBER_X);
        let tunnel = sphere(p, [x, 0.0, 0.0], TUNNEL);
        let mut cave = a.min(b).min(tunnel);

        if self.plugged {
            // Subtract the plug: difference is `max(shape, -cutter)`.
            let plug = sphere(p, [0.0, 0.0, 0.0], PLUG);
            cave = cave.max(-plug);
        }
        // Rock outside the cave, air inside.
        -cave
    }
}

#[derive(Resource)]
struct Demo {
    plugged: bool,
    outline: bool,
    /// Rebuilt when `plugged` changes.
    world: AirWorld<f32>,
    layout: ChunkLayout<f32>,
    /// Last reported search cost.
    visited: u64,
    splits: u64,
    seeds: u64,
    /// Samples in one chunk, and in the whole world, for scale.
    per_chunk: u64,
    world_samples: u64,
}

#[derive(Component)]
struct ChunkMesh(ChunkId);

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — R-028 sealed cave".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, remesh, report, outline))
        .run();
}

/// The chunk ids this demo loads, west to east.
fn ids() -> impl Iterator<Item = ChunkId> {
    (0..CHUNKS).map(|c| ChunkId::new([c, 0, 0]))
}

/// Sample one chunk of `field` into the array [`AirWorld::load`] wants.
fn samples_of(layout: &ChunkLayout<f32>, id: ChunkId, field: &Cave) -> Vec<f32> {
    let n = CELLS + 1;
    let mut values = Vec::with_capacity((n * n * n) as usize);
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let g = layout.global_sample(id, [x, y, z]);
                values.push(field.sample(layout.world_of_sample(g)));
            }
        }
    }
    values
}

/// The chunk-local samples the plug removes, and the chunk holding them.
///
/// The plug sits at the world origin, which is inside the middle chunk — so this
/// is one edit to one chunk, and the search it triggers cannot leave that chunk.
fn plug_samples(layout: &ChunkLayout<f32>) -> (ChunkId, Vec<[u32; 3]>) {
    let id = ChunkId::new([CHUNKS / 2, 0, 0]);
    let mut out = Vec::new();
    let n = CELLS + 1;
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let g = layout.global_sample(id, [x, y, z]);
                let w = layout.world_of_sample(g);
                if sphere(w, [0.0, 0.0, 0.0], PLUG) <= 0.0 {
                    out.push([x, y, z]);
                }
            }
        }
    }
    (id, out)
}

/// A global sample inside each chamber, for the `A↔B` query.
fn chamber_samples(layout: &ChunkLayout<f32>) -> ([i64; 3], [i64; 3]) {
    let mid = CELLS / 2;
    (
        layout.global_sample(ChunkId::new([0, 0, 0]), [mid, mid, mid]),
        layout.global_sample(ChunkId::new([CHUNKS - 1, 0, 0]), [mid, mid, mid]),
    )
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    for mut orbit in &mut camera {
        orbit.yaw = 0.7;
        orbit.pitch = 0.45;
        orbit.radius = 9.0;
    }

    // Centre the world on the origin, so the plug is at the middle chunk's
    // middle and the camera has something symmetric to look at.
    let span = f32::from(CHUNKS as i16) * f32::from(CELLS as i16) * H;
    let origin = [
        -span * 0.5,
        -f32::from(CELLS as i16) * H * 0.5,
        -f32::from(CELLS as i16) * H * 0.5,
    ];
    let Ok(layout) = ChunkLayout::<f32>::new(CELLS, H, origin) else {
        error!("chunk layout rejected {CELLS} cells at {H}");
        return;
    };

    let field = Cave { plugged: false };
    let mut world = AirWorld::new(layout);
    for id in ids() {
        let values = samples_of(&layout, id, &field);
        if let Err(e) = world.load(id, &values) {
            error!("loading {id:?}: {e}");
            return;
        }
    }

    let n = u64::from(CELLS + 1);
    let material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.80, 0.78, 0.72),
        perceptual_roughness: 0.55,
        ..default()
    });
    for id in ids() {
        let mesh = mesh_chunk(&layout, &field, id).unwrap_or_default();
        commands.spawn((
            Mesh3d(meshes.add(to_mesh(&mesh))),
            MeshMaterial3d(material.clone()),
            ChunkMesh(id),
        ));
    }

    commands.insert_resource(Demo {
        plugged: false,
        outline: false,
        world,
        layout,
        visited: 0,
        splits: 0,
        seeds: 0,
        per_chunk: n * n * n,
        world_samples: (u64::from(CHUNKS as u32) * u64::from(CELLS) + 1) * n * n,
    });
}

/// Mesh one chunk on its own, exactly as a dirty-set re-mesh would.
fn mesh_chunk(layout: &ChunkLayout<f32>, field: &Cave, id: ChunkId) -> Option<MeshBuffer<f32>> {
    let shape = layout.sample_shape().ok()?;
    let mut out = MeshBuffer::<f32>::new();
    MarchingCubes::<f32>::new()
        .extract(
            field,
            &shape,
            layout.sample_origin(id),
            layout.cell_size(),
            &mut out,
        )
        .ok()?;
    Some(out)
}

fn controls(keys: Res<ButtonInput<KeyCode>>, demo: Option<ResMut<Demo>>) {
    let Some(mut demo) = demo else { return };
    if keys.just_pressed(KeyCode::KeyF) && !demo.plugged {
        demo.plugged = true;
        let (id, samples) = plug_samples(&demo.layout);
        // The edit, and the whole point of the example: one fill, in one chunk.
        if let Some(f) = demo.world.fill(id, &samples, || true) {
            demo.visited = f.visited;
            demo.splits = f.splits;
            demo.seeds = f.seeds;
        }
    }
    if keys.just_pressed(KeyCode::KeyG) && demo.plugged {
        demo.plugged = false;
        let (id, samples) = plug_samples(&demo.layout);
        if let Some(r) = demo.world.dig(id, &samples, || true) {
            demo.visited = r.relabels;
            demo.splits = 0;
            demo.seeds = 0;
        }
    }
    if keys.just_pressed(KeyCode::KeyC) {
        demo.outline = !demo.outline;
    }
}

fn remesh(
    demo: Option<Res<Demo>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut query: Query<(&mut Mesh3d, &ChunkMesh)>,
    mut last: Local<Option<bool>>,
) {
    let Some(demo) = demo else { return };
    if *last == Some(demo.plugged) {
        return;
    }
    *last = Some(demo.plugged);

    let field = Cave {
        plugged: demo.plugged,
    };
    for (mut handle, chunk) in &mut query {
        // Only the middle chunk's geometry can have moved, but re-meshing all
        // three keeps this loop honest about what it did rather than relying on
        // the plug staying where it is.
        let Some(buffer) = mesh_chunk(&demo.layout, &field, chunk.0) else {
            continue;
        };
        handle.0 = meshes.add(to_mesh(&buffer));
    }
}

fn report(demo: Option<Res<Demo>>, mut stats: ResMut<DemoStats>) {
    let Some(demo) = demo else { return };
    let (a, b) = chamber_samples(&demo.layout);
    let linked = demo.world.connected(a, b);

    stats.title = "R-028 — sealed cave".into();
    stats.extra = alloc_lines(&demo, linked);
}

fn alloc_lines(demo: &Demo, linked: bool) -> Vec<String> {
    let seams = demo.world.last_seams();
    vec![
        format!(
            "chunks {}   samples  chunk {}   world {}",
            demo.world.loaded(),
            demo.per_chunk,
            demo.world_samples
        ),
        format!(
            "components {}      A<->B {}",
            demo.world.components(),
            if linked { "CONNECTED" } else { "SEALED" }
        ),
        format!(
            "last repair: visited {}  seeds {}  splits {}",
            demo.visited, demo.seeds, demo.splits
        ),
        format!(
            "  = {:.3} x one chunk   ({:.3} x the world)",
            demo.visited as f64 / demo.per_chunk as f64,
            demo.visited as f64 / demo.world_samples as f64
        ),
        format!(
            "stitch: {} seams rescanned, {} graph nodes, {} pairs",
            seams.rescanned, seams.nodes, seams.pairs
        ),
        String::new(),
        "F plug the tunnel    G re-open it    C outline chunks".into(),
    ]
}

/// The twelve edges of a unit box, as corner sign pairs.
const EDGES: [([f32; 3], [f32; 3]); 12] = [
    ([-1.0, -1.0, -1.0], [1.0, -1.0, -1.0]),
    ([-1.0, 1.0, -1.0], [1.0, 1.0, -1.0]),
    ([-1.0, -1.0, 1.0], [1.0, -1.0, 1.0]),
    ([-1.0, 1.0, 1.0], [1.0, 1.0, 1.0]),
    ([-1.0, -1.0, -1.0], [-1.0, 1.0, -1.0]),
    ([1.0, -1.0, -1.0], [1.0, 1.0, -1.0]),
    ([-1.0, -1.0, 1.0], [-1.0, 1.0, 1.0]),
    ([1.0, -1.0, 1.0], [1.0, 1.0, 1.0]),
    ([-1.0, -1.0, -1.0], [-1.0, -1.0, 1.0]),
    ([1.0, -1.0, -1.0], [1.0, -1.0, 1.0]),
    ([-1.0, 1.0, -1.0], [-1.0, 1.0, 1.0]),
    ([1.0, 1.0, -1.0], [1.0, 1.0, 1.0]),
];

fn outline(demo: Option<Res<Demo>>, mut gizmos: Gizmos) {
    let Some(demo) = demo else { return };
    if !demo.outline {
        return;
    }
    let side = f32::from(CELLS as i16) * H;
    for id in ids() {
        let o = demo.layout.sample_origin(id);
        let centre = Vec3::new(o[0] + side * 0.5, o[1] + side * 0.5, o[2] + side * 0.5);
        let colour = if id.coords[0] == CHUNKS / 2 {
            Color::srgb(0.95, 0.65, 0.25)
        } else {
            Color::srgb(0.35, 0.45, 0.60)
        };
        let h = side * 0.5;
        // The twelve edges of the box, as pairs of corner sign patterns.
        for (u, v) in EDGES {
            let corner = |m: [f32; 3]| centre + Vec3::new(m[0] * h, m[1] * h, m[2] * h);
            gizmos.line(corner(u), corner(v), colour);
        }
    }
}
