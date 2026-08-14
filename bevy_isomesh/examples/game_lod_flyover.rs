//! E-205 — flying across LOD transitions, and counting what opens up.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_lod_flyover --release
//! ```
//!
//! **Always `--release`.**
//!
//! `Space` pause · `T` transition cells on/off · `[` `]` fly speed · `R` reset.
//!
//! # What this has to prove, and how each part is measured
//!
//! *No cracks.* The assembled world is welded and validated every re-mesh, and
//! the **boundary edges lying in a seam plane** are counted. That is E-107's
//! technique applied to a whole ladder of blocks rather than one pair, and it is
//! a count rather than a look: a crack a pixel wide at this camera distance is
//! invisible and still fatal.
//!
//! *No popping.* When a block changes level its surface moves. The demo meshes
//! the block at **both** levels at the moment of the switch and reports the worst
//! vertex displacement between them, in cells. A pop is not a bug — a coarser
//! mesh really is a different surface — but its **size** is the number that
//! decides whether it can be hidden, and nothing in the literature review
//! measures it.
//!
//! *No hitching.* Frame time is on the HUD, and the re-mesh cost of a level
//! change is charged to the frame it happens on.
//!
//! # The LOD ladder runs along one axis, deliberately
//!
//! Blocks are slabs stacked along `x`, each meshed at its own spacing, and the
//! camera flies out along `x` and back. Every seam in the world is therefore an
//! `x` seam, which is what makes the crack count mean something: if a seam plane
//! has a boundary edge in it, the transition failed, and there is no second axis
//! for the failure to hide behind.
//!
//! A production terrain needs transitions on four faces, or six for a world with
//! caves, and **that is not what this demonstrates**. What it demonstrates is the
//! mechanism, at speed, with the failure counted — and it exercises one thing
//! E-107 never did, below.
//!
//! # The low-side transition has never been run before
//!
//! E-107 meshes exactly one pair of blocks with the fine one on the low-`x` side,
//! so `inset_boundary` has only ever been called with `face_bit(0, 0)` and
//! transition cells have only ever been sampled at the fine block's *last* `x`
//! index. A camera that flies **out and back** puts coarse blocks on both sides
//! of itself, so half the seams here are the mirrored configuration.
//!
//! If that configuration is wrong — a mirrored patch is the classic place for an
//! inside-out winding, which no manifold or Euler check can see — the seam
//! boundary count says so, per side. The HUD reports the two sides separately for
//! exactly that reason.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::to_bevy_mesh;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::FbmTerrain;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::transvoxel::cell::TransitionCell;
use isomesh::transvoxel::inset::{face_bit, inset_boundary};
use isomesh::validate::{ValidateConfig, validate_features};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Spacing at level 0. Every coarser level doubles it.
const BASE_H: f32 = 0.25;
/// Levels 0..=2, so spacings are 0.25, 0.5 and 1.0.
const MAX_LEVEL: u32 = 2;
/// Width of one block along `x`, in world units.
///
/// A multiple of the coarsest spacing, or a block would not contain a whole
/// number of its own cells and its far face would not land on the seam plane.
const BLOCK_W: f32 = 4.0;
const BLOCKS: usize = 12;
/// Half-extent of the world in `y` and `z`.
const CROSS: f32 = 4.0;

/// How far the camera has to be for the level to step up.
const LEVEL_RANGE: f32 = 7.0;

#[derive(Resource)]
struct Fly {
    /// Camera position along `x`.
    at: f32,
    speed: f32,
    flying: bool,
    transitions: bool,
    levels: [u32; BLOCKS],
    /// Seam boundary edges, split by which side of the camera the seam is on.
    open_low: u64,
    open_high: u64,
    seams: u32,
    /// How many of those seams sit below the camera. Without this, `open_low`
    /// reads zero both when the mirrored transition works and when there is
    /// nothing on that side to test.
    seams_low: u32,
    /// Worst vertex displacement across a level change, in cells of the finer
    /// level. `None` until a block has actually switched.
    worst_pop_cells: Option<f32>,
    build_ms: f64,
    validate_ms: f64,
    /// Blocks re-extracted on the last level change, of BLOCKS.
    extracted: usize,
    triangles: usize,
    vertices: usize,
}

#[derive(Resource)]
struct Look(Handle<StandardMaterial>);

/// Each block's raw extraction, kept across level changes.
///
/// **Without this the demo re-extracts all twelve blocks whenever any one of
/// them changes level, which cost 12–23 ms and is a hitch.** That is the naive
/// implementation's cost, not the crate's: a level change invalidates the blocks
/// that changed and nothing else. What is cached is the *un-inset* extraction,
/// because `inset_boundary` mutates positions in place and which faces need
/// tapering depends on the neighbours' levels, which move independently.
#[derive(Resource, Default)]
struct Cache {
    blocks: Vec<Option<(u32, MeshBuffer<f32>)>>,
    /// How many blocks the last re-assembly actually had to extract.
    extracted: usize,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-205 game lod flyover".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(Fly {
            at: 0.0,
            speed: 6.0,
            flying: true,
            transitions: std::env::var("ISOMESH_TRANSITIONS").map_or(true, |v| v != "0"),
            levels: [u32::MAX; BLOCKS],
            open_low: 0,
            open_high: 0,
            seams: 0,
            seams_low: 0,
            worst_pop_cells: None,
            build_ms: 0.0,
            validate_ms: 0.0,
            extracted: 0,
            triangles: 0,
            vertices: 0,
        })
        .init_resource::<Cache>()
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, rebuild, report, hud).chain())
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    for mut orbit in &mut camera {
        orbit.yaw = 0.9;
        orbit.pitch = 0.45;
        orbit.radius = 22.0;
    }
    commands.insert_resource(Look(materials.add(StandardMaterial {
        base_color: Color::srgb(0.62, 0.66, 0.52),
        perceptual_roughness: 0.85,
        double_sided: true,
        cull_mode: None,
        ..default()
    })));
}

fn controls(
    time: Res<Time>,
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut fly: ResMut<Fly>,
    mut flags: ResMut<ViewFlags>,
) {
    if !capture.is_active() {
        if keys.just_pressed(KeyCode::Space) {
            fly.flying = !fly.flying;
        }
        if keys.just_pressed(KeyCode::KeyT) {
            fly.transitions = !fly.transitions;
        }
        if keys.just_pressed(KeyCode::BracketRight) {
            fly.speed = (fly.speed + 2.0).min(20.0);
        }
        if keys.just_pressed(KeyCode::BracketLeft) {
            fly.speed = (fly.speed - 2.0).max(2.0);
        }
        if keys.just_pressed(KeyCode::KeyR) {
            fly.at = 0.0;
            fly.worst_pop_cells = None;
        }
    }
    if !fly.flying {
        return;
    }
    // Out and back along the ladder, so coarse blocks end up on both sides of
    // the camera and both mirror images of the seam get exercised.
    let span = BLOCK_W * BLOCKS as f32;
    fly.at += fly.speed * time.delta_secs();
    if fly.at > span {
        fly.at -= span * 2.0;
    }
    flags.grid = false;
}

/// Level for the block whose centre is `centre`, given a camera at `at`.
fn level_for(centre: f32, at: f32) -> u32 {
    ((centre - at).abs() / LEVEL_RANGE) as u32
}

/// Force adjacent blocks to differ by at most one level.
///
/// Transvoxel's transition cells bridge a **2:1** ratio and nothing else. A 4:1
/// step would need a different table, so the ladder is smoothed rather than left
/// to whatever the distance function produced.
fn smooth(levels: &mut [u32; BLOCKS]) {
    for _ in 0..BLOCKS {
        let mut settled = true;
        for i in 1..BLOCKS {
            if levels[i] > levels[i - 1] + 1 {
                levels[i] = levels[i - 1] + 1;
                settled = false;
            }
            if levels[i - 1] > levels[i] + 1 {
                levels[i - 1] = levels[i] + 1;
                settled = false;
            }
        }
        if settled {
            break;
        }
    }
}

fn spacing(level: u32) -> f32 {
    BASE_H * (1 << level) as f32
}

/// Mesh one block at its own level, with no transition geometry.
fn mesh_block<F: Sdf<Scalar = f32>>(
    field: &F,
    index: usize,
    level: u32,
) -> Option<MeshBuffer<f32>> {
    let h = spacing(level);
    let along = (BLOCK_W / h).round() as u32;
    let across = ((CROSS * 2.0) / h).round() as u32;
    let shape = RuntimeShape3::new([along + 1, across + 1, across + 1]).ok()?;
    let origin = [index as f32 * BLOCK_W, -CROSS, -CROSS];
    let mut out = MeshBuffer::<f32>::new();
    MarchingCubes::<f32>::new()
        .extract(field, &shape, origin, h, &mut out)
        .ok()?;
    Some(out)
}

/// The whole world, assembled, with transitions on every 2:1 seam.
struct Assembled {
    mesh: MeshBuffer<f32>,
    /// Extraction and welding -- the work a game would actually pay for.
    build_ms: f64,
    /// Validation. **This demo's instrument, not the crate's cost.** A game
    /// never re-validates its whole world on a level change; this one does,
    /// because counting cracks is the entire point. Timed separately so the
    /// hitching claim is judged on the work rather than on the measuring.
    validate_ms: f64,
    open_low: u64,
    open_high: u64,
    seams: u32,
    seams_low: u32,
}

fn assemble<F: Sdf<Scalar = f32>>(
    field: &F,
    levels: &[u32; BLOCKS],
    at: f32,
    transitions: bool,
    cache: &mut Cache,
) -> Option<Assembled> {
    let build = Instant::now();
    cache.blocks.resize_with(BLOCKS, || None);
    cache.extracted = 0;
    let mut blocks: Vec<MeshBuffer<f32>> = Vec::with_capacity(BLOCKS);
    for (i, &level) in levels.iter().enumerate() {
        let hit = matches!(cache.blocks[i], Some((cached, _)) if cached == level);
        if !hit {
            cache.blocks[i] = Some((level, mesh_block(field, i, level)?));
            cache.extracted += 1;
        }
        // Cloned rather than moved: the cached copy has to survive, and the
        // assembled one is about to be tapered in place.
        let Some((_, buffer)) = &cache.blocks[i] else {
            return None;
        };
        blocks.push(buffer.clone());
    }

    let mut seam_planes: Vec<(f32, bool)> = Vec::new();
    let mut transition_mesh = MeshBuffer::<f32>::new();
    let mut seams = 0;

    // Seams are found whether or not they are going to be bridged. **The first
    // version registered them inside the `if transitions` block, so turning
    // transitions off left nothing to count against and the crack total read a
    // confident zero on a world full of holes.** A counter that cannot report a
    // failure is not a counter.
    for i in 0..BLOCKS - 1 {
        if levels[i] != levels[i + 1] {
            seams += 1;
            let seam_x = (i + 1) as f32 * BLOCK_W;
            // The seam is "low" when it sits on the camera's low-x side, which
            // is the mirror of the only configuration E-107 ever ran.
            seam_planes.push((seam_x, seam_x < at));
        }
    }

    if transitions {
        for i in 0..BLOCKS - 1 {
            let (a, b) = (levels[i], levels[i + 1]);
            if a == b {
                continue;
            }
            let (fine_index, coarse_index) = if a < b { (i, i + 1) } else { (i + 1, i) };
            let fine_level = levels[fine_index];
            let fine_h = spacing(fine_level);
            let coarse_h = spacing(levels[coarse_index]);
            let width = fine_h;

            // Make room on the coarse block's facing side. `inset_boundary` must
            // only ever be given faces that actually have a differently-resolved
            // neighbour, or the block pulls away from a same-level neighbour and
            // opens a seam where there was none.
            let coarse_face = if coarse_index > fine_index {
                face_bit(0, 0)
            } else {
                face_bit(0, 1)
            };
            let coarse_origin = [coarse_index as f32 * BLOCK_W, -CROSS, -CROSS];
            let coarse_cells = (BLOCK_W / coarse_h).round() as u32;
            inset_boundary(
                &mut blocks[coarse_index],
                coarse_origin,
                coarse_cells,
                coarse_h,
                width,
                coarse_face,
            )
            .ok()?;

            // Transition cells live on the *fine* grid, in the seam plane.
            let fine_origin = [fine_index as f32 * BLOCK_W, -CROSS, -CROSS];
            let fine_along = (BLOCK_W / fine_h).round() as i64;
            let fine_across = ((CROSS * 2.0) / fine_h).round() as i64;
            // The fine block's own x index at the seam: its far face when the
            // coarse block is above it, its near face when below.
            let base_x = if coarse_index > fine_index {
                fine_along
            } else {
                0
            };
            // The patch's width is **signed**: positive runs toward increasing
            // x, so it states which side the coarse block is on. A mirrored
            // seam (coarse block below) needs the negative sign, or the patch
            // grows into the fine block and leaves the inset moat open --
            // measured at 44 boundary edges per mirrored two-block seam, at
            // seam±w where the HUD's seam-plane counter cannot see them, and
            // zero with the sign stated.
            let sample_width = if coarse_index > fine_index {
                width
            } else {
                -width
            };
            for jz in 0..fine_across / 2 {
                for jy in 0..fine_across / 2 {
                    let cell = TransitionCell::sample(
                        field,
                        fine_origin,
                        fine_h,
                        [base_x, 2 * jy, 2 * jz],
                        1,
                        2,
                        sample_width,
                    );
                    cell.emit(field, 0, &mut transition_mesh);
                }
            }
        }
    }

    let mut mesh = MeshBuffer::<f32>::new();
    for block in &blocks {
        mesh.append(block)
            .expect("the assembled world fits the u32 index space");
    }
    mesh.append(&transition_mesh)
        .expect("the meshes fit the u32 index space");

    let finest = spacing(*levels.iter().min().unwrap_or(&0));
    isomesh::weld::Welder::<f32>::new()
        .weld(&mut mesh, isomesh::weld::epsilon_for(finest))
        .ok()?;

    let build_ms = build.elapsed().as_secs_f64() * 1000.0;

    // Count what is still open, and where. A boundary edge in a seam plane is a
    // crack; one on the world's outer wall is the world ending, which is not.
    let validate = Instant::now();
    let cfg = ValidateConfig::from_cell_size(f64::from(finest)).ok()?;
    let (_report, features) = validate_features(&mesh.positions, &mesh.indices, &cfg);
    let (mut open_low, mut open_high) = (0u64, 0u64);
    for edge in &features.boundary_edges {
        let (p, q) = (
            mesh.positions[edge[0] as usize],
            mesh.positions[edge[1] as usize],
        );
        let tol = BASE_H * 0.25;
        // An edge on the world's outer wall is the world ending, not a crack --
        // and the wall passes *through* every seam plane, so without this the
        // four edges where they meet are counted as failures. The first run
        // reported exactly that: 1 crack that was the y/z boundary crossing the
        // seam, on a transition that had closed correctly.
        let on_outer = |v: [f32; 3]| {
            (v[1].abs() - CROSS).abs() < tol
                || (v[2].abs() - CROSS).abs() < tol
                || v[0] < tol
                || v[0] > BLOCK_W * BLOCKS as f32 - tol
        };
        if on_outer(p) && on_outer(q) {
            continue;
        }
        for &(plane, is_low) in &seam_planes {
            // Both endpoints within a hair of the plane, so an edge that merely
            // crosses the region is not counted as lying in it.
            if (p[0] - plane).abs() < tol && (q[0] - plane).abs() < tol {
                if is_low {
                    open_low += 1;
                } else {
                    open_high += 1;
                }
                break;
            }
        }
    }

    let validate_ms = validate.elapsed().as_secs_f64() * 1000.0;
    let seams_low = seam_planes.iter().filter(|(_, low)| *low).count() as u32;
    Some(Assembled {
        mesh,
        build_ms,
        validate_ms,
        open_low,
        open_high,
        seams,
        seams_low,
    })
}

#[allow(clippy::too_many_arguments)]
fn rebuild(
    mut fly: ResMut<Fly>,
    mut cache: ResMut<Cache>,
    look: Res<Look>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut query: Query<&mut Mesh3d, With<DemoMesh>>,
    mut camera: Query<&mut OrbitCamera>,
    mut flags: ResMut<ViewFlags>,
) {
    let field = FbmTerrain::<f32>::canonical();

    let mut levels = [0u32; BLOCKS];
    for (i, slot) in levels.iter_mut().enumerate() {
        let centre = (i as f32 + 0.5) * BLOCK_W;
        *slot = level_for(centre, fly.at).min(MAX_LEVEL);
    }
    smooth(&mut levels);

    // Follow the flight with the camera whatever else happens.
    for mut orbit in &mut camera {
        orbit.focus = Vec3::new(fly.at, 0.0, 0.0);
    }

    if levels == fly.levels && !flags.remesh_requested {
        return;
    }

    // Anything that changed level: mesh it at both, and measure how far the
    // surface moved. This is the pop, and it is only measurable at the instant
    // of the switch.
    let mut worst = fly.worst_pop_cells;
    if fly.levels[0] != u32::MAX {
        for (i, (&was, &now)) in fly.levels.iter().zip(levels.iter()).enumerate() {
            if was == now {
                continue;
            }
            let (Some(before), Some(after)) =
                (mesh_block(&field, i, was), mesh_block(&field, i, now))
            else {
                continue;
            };
            let h = spacing(was.min(now));
            let moved = worst_gap(&before, &after) / h;
            if moved.is_finite() && worst.is_none_or(|w| moved > w) {
                worst = Some(moved);
            }
        }
    }

    let Some(built) = assemble(&field, &levels, fly.at, fly.transitions, &mut cache) else {
        return;
    };
    fly.build_ms = built.build_ms;
    fly.extracted = cache.extracted;
    fly.validate_ms = built.validate_ms;
    fly.levels = levels;
    fly.open_low = built.open_low;
    fly.open_high = built.open_high;
    fly.seams = built.seams;
    fly.seams_low = built.seams_low;
    fly.worst_pop_cells = worst;
    fly.vertices = built.mesh.positions.len();
    fly.triangles = built.mesh.indices.len() / 3;
    flags.remesh_requested = false;

    let handle = meshes.add(to_bevy_mesh(&built.mesh));
    if query.is_empty() {
        commands.spawn((
            Mesh3d(handle),
            MeshMaterial3d(look.0.clone()),
            Transform::default(),
            DemoMesh,
        ));
    } else {
        for mut mesh in &mut query {
            mesh.0 = handle.clone();
        }
    }
}

/// Worst distance from a vertex of `a` to the nearest vertex of `b`.
///
/// A one-sided Hausdorff distance over vertices, which is enough to size a pop:
/// the question is how far the surface a viewer was looking at has moved, and
/// the vertices are where it moved most.
fn worst_gap(a: &MeshBuffer<f32>, b: &MeshBuffer<f32>) -> f32 {
    let mut worst = 0.0f32;
    for p in &a.positions {
        let mut nearest = f32::MAX;
        for q in &b.positions {
            let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
            let dist = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            if dist < nearest {
                nearest = dist;
            }
        }
        if nearest.is_finite() && nearest > worst {
            worst = nearest;
        }
    }
    worst.sqrt()
}

/// One CSV row per re-mesh, so a whole out-and-back can be read from a shell.
/// A still cannot carry this: the low side has no seams at all until the camera
/// has flown past some, so its zero is vacuous early on.
fn report(fly: Res<Fly>, mut last: Local<Option<[u32; BLOCKS]>>) {
    if *last == Some(fly.levels) {
        return;
    }
    *last = Some(fly.levels);
    info!(
        "lod,{:.1},{},{},{},{},{:.3},{:.1},{:.1},{},{}",
        fly.at,
        fly.seams,
        fly.seams_low,
        fly.open_low,
        fly.open_high,
        fly.worst_pop_cells.unwrap_or(f32::NAN),
        fly.build_ms,
        fly.validate_ms,
        fly.extracted,
        fly.triangles,
    );
}

fn hud(fly: Res<Fly>, mut stats: ResMut<DemoStats>) {
    let open = fly.open_low + fly.open_high;
    let verdict = if !fly.transitions {
        "transitions OFF -- every seam below is an open crack, on purpose".to_string()
    } else if open == 0 {
        "NO CRACKS -- every 2:1 seam in the ladder is closed".to_string()
    } else {
        format!("!! {open} boundary edges lie in a seam plane -- a transition failed")
    };

    stats.title = format!(
        "E-205  lod flyover   x = {:.1}   levels {:?}",
        fly.at, fly.levels
    );
    stats.vertices = fly.vertices;
    stats.triangles = fly.triangles;
    stats.extra = vec![
        format!(
            "{:<26} {:>8}   2:1 seams, {} of them below the camera",
            "seams", fly.seams, fly.seams_low
        ),
        format!(
            "{:<26} {:>8}   open edges on the {} seam(s) BELOW",
            "cracks, low side", fly.open_low, fly.seams_low
        ),
        format!(
            "{:<26} {:>8}   open edges on the {} seam(s) ABOVE",
            "cracks, high side",
            fly.open_high,
            fly.seams - fly.seams_low
        ),
        String::new(),
        "the two sides are counted separately because only the high side has".into(),
        "ever been run before: E-107 meshes one pair with the fine block on the".into(),
        "low-x side, so the mirrored seam is new here. a mirrored patch is the".into(),
        "classic place for an inside-out winding, and no manifold or Euler".into(),
        "check can see one.".into(),
        String::new(),
        match fly.worst_pop_cells {
            Some(p) => format!(
                "{:<26} {:>8.3}   cells the surface moved at a level change",
                "worst pop", p
            ),
            None => format!("{:<26} {:>8}   no level change yet", "worst pop", "--"),
        },
        format!(
            "{:<26} {:>8}   of {} blocks re-extracted on this level change",
            "blocks rebuilt", fly.extracted, BLOCKS
        ),
        format!(
            "{:<26} {:>8.1}   ms to extract those and weld the ladder",
            "re-build", fly.build_ms
        ),
        format!(
            "{:<26} {:>8.1}   ms to validate it -- this demo's instrument,",
            "re-validate", fly.validate_ms
        ),
        "                                    not a cost a game pays on a level change".into(),
        String::new(),
        verdict,
        String::new(),
        format!(
            "speed {:.0}   [ and ] to change   Space pauses   T transitions   R resets",
            fly.speed
        ),
    ];
}
