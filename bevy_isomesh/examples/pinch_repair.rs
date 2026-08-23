//! E-310 — the repair is exact and it is still not the thing to ship.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example pinch_repair --release
//! ```
//!
//! **Always `--release`.** `bonsai` is 256³: a debug build spends minutes in the
//! march alone, and this example runs it twice before the window opens.
//!
//! `1` `fuel`, `2` `bonsai`, `B` switches between the baseline and the repaired
//! arm, and the mouse wheel zooms. The rest are the shared keys — `G` domain
//! box, `F12` screenshot, `Esc` quit. **`W` is a trap here**: the shared
//! wireframe draws every triangle edge as a gizmo line and `bonsai` has a
//! million triangles.
//!
//! Under `ISOMESH_CAPTURE` it needs no keyboard: the capture is four equal
//! stages — `fuel` baseline, `fuel` repaired, `bonsai` baseline, `bonsai`
//! repaired — so a clip of any length carries the whole comparison.
//! `ISOMESH_FIELD=0`/`=1` pins one volume and the two arms then split the clip
//! in half.
//!
//! ```bash
//! # 1024x576 rather than the script's default 1280x720, and that is measured
//! # rather than a preference. `record_gif.sh` scales to 900px wide, so 1280
//! # shrinks the harness's 13px HUD to about 9px; at 1024 it lands near 11.4px
//! # and the topology table is comfortably readable in the clip. This example's
//! # whole argument is the numbers in the corner, so the HUD sets the
//! # capture size rather than the other way round.
//! ISOMESH_WINDOW=1024x576 ISOMESH_SPIN=0.003 \
//!   ./scripts/record_gif.sh pinch_repair docs/gifs/e310.gif
//! ```
//!
//! Nothing but the camera moves between frames — the meshes are prebuilt and the
//! arm switch is a handle swap — so a GIF of this compresses well at the
//! script's default palette and dither.
//!
//! Both arms of both volumes are meshed **once, at startup**, and the toggles
//! swap a `Handle<Mesh>`. `bonsai`'s march is ~500 ms and there are two of them
//! plus two validator passes, so this is the one thing that must not be
//! per-frame.
//!
//! Demonstrates **M-352 / P-53**.
//!
//! # What the repair does, and what it costs
//!
//! This crate labels a corner with `value < 0`, so a sample sitting *exactly* on
//! the isosurface is outside. On a cut edge whose outside endpoint is that
//! exactly-zero corner, `t = a / (a − b)` is `0 / (0 − b)` = **exactly 0**, so
//! the crossing lands **on the corner**; every cut edge meeting there does the
//! same thing, and the cell emits a triangle whose corners are all one point.
//! Custodio, Pesco & Silva (`10.1186/s13173-019-0086-6`) name the fix: a **third
//! corner label**, so the equal case is not folded into one of the two sides.
//!
//! P-53 replayed the march for provenance, applied that label as a pre-pass, and
//! measured both arms. The repair works, exactly and completely:
//!
//! | | `fuel` 64³ | `bonsai` 256³ |
//! |---|---:|---:|
//! | degenerate triangles | **164 → 0** | **58,097 → 0** |
//! | `max_snap_distance` | **0** | **0** |
//! | Euler characteristic | 19 → 19 | **517 → 585** |
//! | non-manifold edges | 0 → 0 | **0 → 561** |
//! | boundary edges | 24 → 24 | **4366 → 3716** |
//! | collapse groups | 50 | 17,201 |
//! | **pinch groups** | **0** | **516** |
//! | components welded | 0 | **520** |
//!
//! `max_snap_distance` is exactly zero on both, and that is the load-bearing
//! number: with `origin = 0`, `h = 1` and integer samples, `t` is exactly 0 or 1
//! on an edge incident to an `=` corner, so `lo + (hi − lo)·t` is bit-exactly a
//! corner position. **Nothing moves.** The label's entire effect is which
//! vertices are declared to be the same vertex, which makes this a pure
//! connectivity decision — and a connectivity decision is exactly the kind that
//! can rewire a surface without displacing a single point.
//!
//! # The pinch is the whole finding
//!
//! Two vertices snapped to one `=` corner are in one of two situations:
//!
//! - They **already share a triangle**. That triangle is one of the degenerate
//!   ones; merging them flattens a fold. No edge, boundary or component can
//!   move. Drawn **amber**.
//! - They **share no triangle**. They are on different pieces of the surface
//!   that happen to meet at that sample — the isosurface genuinely touches
//!   itself there — and identifying the point is a change of topology no
//!   relabelling can avoid. Drawn **magenta**, and marked with a cross.
//!
//! `fuel` has **0 pinches of 50 groups**. Every collapse there is a fold, the
//! three topology counters do not budge, and the repair is free. `bonsai` has
//! **516 of 17,201**, and those 516 weld **520** previously separate pieces
//! together: χ moves 517 → 585, 561 non-manifold edges appear where there were
//! none, and 650 boundary edges vanish into the welds.
//!
//! So the shippable result is the **precondition**, not the repair. A pipeline
//! that guarantees no sample sits exactly on the isovalue — contour at a
//! half-integer, which `u8` data cannot attain — gets `fuel`'s row on every
//! volume. A pipeline that applies the label to data that *does* have equal
//! corners gets a mesh with no degenerate triangles and a different topology
//! than the one it extracted.
//!
//! # What is on screen
//!
//! The surface, at isovalue 32, painted per vertex:
//!
//! - **Dark grey** — ordinary surface.
//! - **Amber** — a vertex in a collapse group that already shared a triangle: a
//!   fold about to be flattened. Present in the baseline arm, gone in the
//!   repaired one, because there is nothing left to see there.
//! - **Magenta** — a vertex in a **pinch** group. Present in *both* arms: before,
//!   two sheets about to be welded; after, the weld. `fuel` never shows one.
//! - **Red** — a vertex of a degenerate triangle that the `=` label does not
//!   reach: the third corner of a sliver whose other two coincide. Reported on
//!   the HUD rather than claimed to be zero.
//! - **Magenta crosses** — the 516 pinch sites, depth-biased so they read
//!   through the tree. There are none on `fuel`, and that empty overlay is half
//!   the demonstration.
//!
//! A degenerate triangle has zero area and draws nothing, so it cannot be shaded;
//! what is highlighted is the site, from the same exact census the ledger used.
//!
//! # `f64`, Marching Cubes, and the ledger
//!
//! P-53 was measured in `f64` under `isomesh::marching_cubes`, so the mesh is
//! extracted that way and cast to `f32` on its way into the [`Mesh`] asset.
//! Everything the HUD reports is re-derived here from the raw bytes and then
//! **compared column by column against `docs/experiments/p-53.csv`**, including
//! `mesh_hash` — a 64-bit hash of both arms' vertex and index buffers. If this
//! file's replay were off by one triangle, the hash would say so.
//!
//! The march is replayed for provenance only — cell iteration order, the
//! eight-corner case index, the vertex cache keyed on `(lower sample, axis)` —
//! and its index buffer is compared against the crate's element for element.
//! `replay_matches_crate` is on the ledger and is checked here too; the tagging
//! is only licensed because that comparison holds.

#![allow(
    clippy::float_cmp,
    reason = "exact equality is the phenomenon: a corner on the isosurface, a \
              vertex on that corner, a triangle of exactly zero area"
)]

mod common;

use std::path::{Path, PathBuf};
use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::construct::SampledField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_cubes::table::{self, CASES, EDGE_AXIS, EDGE_CORNERS, is_inside};
use isomesh::validate::{ValidateConfig, mesh_hash, validate_indexed};
use isomesh::{MeshBuffer, MeshSink, RuntimeShape3, Sdf, Shape3};

// ─── the registered fixture ─────────────────────────────────────────────────

/// One world unit per voxel: both datasets declare `spacing 1x1x1`.
const CELL: f64 = 1.0;

/// Sample `[0, 0, 0]` sits at the world origin, so a corner's world position is
/// its integer grid coordinate exactly. That is what makes `max_snap_distance`
/// meaningful as an absolute number rather than a relative one.
const ORIGIN: [f64; 3] = [0.0; 3];

/// P-53's integer isovalue. Integer is the point: a `u8` sample **can** equal
/// it, which is the input class the whole finding is about.
const ISOVALUE: f64 = 32.0;

/// How the ledger spells [`ISOVALUE`] in its `isovalue` column.
const ISOVALUE_LABEL: &str = "32";

/// The raw byte that sits exactly on [`ISOVALUE`].
const EQUAL: u8 = 32;

/// A volume to read.
struct Volume {
    /// File name under `docs/measurements/volumes`.
    file: &'static str,
    /// How the ledger spells it in its `volume` column.
    short: &'static str,
    /// Samples per axis; the files are cubes.
    n: u32,
}

/// The two CT volumes P-53 measured, in the order the digit keys select them.
///
/// `fuel` leads because it is the clean row: 50 collapse groups and not one
/// pinch. An example that only ever showed `bonsai` would leave "the repair
/// rewires the surface" looking like a property of the repair rather than of the
/// data it is given.
const VOLUMES: [Volume; 2] = [
    Volume {
        file: "fuel_64x64x64_uint8.raw",
        short: "fuel",
        n: 64,
    },
    Volume {
        file: "bonsai_256x256x256_uint8.raw",
        short: "bonsai",
        n: 256,
    },
];

/// How many volumes the digit keys offer.
const VOLUME_COUNT: usize = VOLUMES.len();

/// Captured frames are split into this many stages: each volume's two arms.
const STAGES: u32 = 4;

// ─── the vertex classification, and its colours ─────────────────────────────

/// An ordinary vertex: no collapse group, no degenerate triangle.
const CLASS_PLAIN: u8 = 0;

/// A vertex of a degenerate triangle that the `=` label does not reach — the
/// third corner of a sliver whose other two coincide.
const CLASS_DEGENERATE: u8 = 1;

/// A member of a collapse group whose vertices already shared a triangle: a fold
/// about to be flattened, and nothing else.
const CLASS_SAFE: u8 = 2;

/// A member of a collapse group whose vertices shared **no** triangle: two
/// sheets about to be welded.
const CLASS_PINCH: u8 = 3;

/// Ordinary surface.
///
/// **Darker than this repo's usual surface grey, measured rather than chosen.**
/// Both volumes fill the frame, so the whole HUD sits on top of the
/// surface, and Bevy's tonemapper compresses hard: measured on a 1280x720
/// capture at the harness's ambient 220 lux, a `0.29` grey lights to `168/255`,
/// `0.14` to `140` and this to about `100`, against HUD text at `235`. The
/// numbers are the evidence, so the rock loses — and the three accents below
/// have to be the brightest thing in frame anyway.
const COLOUR_PLAIN: [f32; 4] = [0.05, 0.055, 0.085, 1.0];

/// [`CLASS_DEGENERATE`].
const COLOUR_DEGENERATE: [f32; 4] = [1.0, 0.30, 0.18, 1.0];

/// [`CLASS_SAFE`].
const COLOUR_SAFE: [f32; 4] = [1.0, 0.71, 0.10, 1.0];

/// [`CLASS_PINCH`].
const COLOUR_PINCH: [f32; 4] = [1.0, 0.16, 0.86, 1.0];

/// The colour a class paints its vertex.
fn class_colour(class: u8) -> [f32; 4] {
    match class {
        CLASS_DEGENERATE => COLOUR_DEGENERATE,
        CLASS_SAFE => COLOUR_SAFE,
        CLASS_PINCH => COLOUR_PINCH,
        _ => COLOUR_PLAIN,
    }
}

// ─── camera framing ─────────────────────────────────────────────────────────

/// Wide view distance, as a multiple of the mesh's **largest half-extent**.
///
/// The largest half-extent, not the bounding-sphere radius. Both volumes throw
/// off isolated specks at isovalue 32 — a `u8` CT scan has noise below the
/// isosurface everywhere — and those specks inflate the bounding sphere without
/// contributing anything to look at, so framing on it put the subject at 43% of
/// the frame width. Measured at 1280x720 and re-checked at the 900px the GIFs
/// are published at.
///
/// An exact vertical fit at Bevy's default 45° field of view is
/// `1 / tan(22.5°)` = 2.414; this is a hair tighter, so a corner-on view crops
/// the box edges rather than leaving a margin.
const WIDE_RADIUS: f32 = 2.30;

/// Half-length of a pinch cross, as a fraction of the volume's largest extent.
///
/// Relative rather than absolute because the two volumes differ by 4x, and a
/// marker sized for `bonsai` would swallow `fuel`.
const MARKER_FRACTION: f32 = 1.0 / 58.0;

// ─── reading the volume ─────────────────────────────────────────────────────

/// Where the `.raw` files live.
fn volume_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../docs/measurements/volumes")
}

/// Where P-53's committed rows live.
fn ledger_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../docs/experiments/p-53.csv")
}

/// Read a `uint8` volume, raw. The length is checked against the dimensions
/// rather than trusted from the filename, as `benches/volumes.rs` does.
fn read_u8(path: &Path, n: u32) -> Result<Vec<u8>, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let want = (n as usize).pow(3);
    if bytes.len() != want {
        return Err(format!(
            "{}: {} bytes, expected {want} for {n}^3",
            path.display(),
            bytes.len()
        ));
    }
    Ok(bytes)
}

// ─── the march, replayed for provenance ─────────────────────────────────────

/// The local offset of cube corner `i`, as grid steps: bit 0 is x, bit 1 is y,
/// bit 2 is z.
///
/// `isomesh::cube::corner_offset` is `pub(crate)`. This is its documented
/// definition rather than a second convention.
const fn corner_offset(corner: u8) -> [u32; 3] {
    [
        (corner & 1) as u32,
        ((corner >> 1) & 1) as u32,
        ((corner >> 2) & 1) as u32,
    ]
}

/// Marching Cubes' own vertex numbering and cell attribution, recovered without
/// recomputing a single position.
struct Provenance {
    /// The index buffer in emission order, to be compared against the crate's.
    indices: Vec<u32>,
    /// Per triangle: the linear sample index of its cell's base corner.
    cell: Vec<u32>,
    /// Per vertex: the lower sample of its grid edge, or `u32::MAX` for a
    /// cell-local cycle centroid.
    edge_lo: Vec<u32>,
    /// Per vertex: the axis of that grid edge, or `3` for a centroid.
    edge_axis: Vec<u8>,
    /// Cells with at least one corner inside and at least one not.
    surface_cells: u64,
    /// Cells whose triangulation needed a cycle centroid. Plain Marching Cubes
    /// tops out at cycle length 7 and `safe_apex` covers 3..=7, so this is zero;
    /// it is counted rather than assumed, because a centroid vertex is allocated
    /// on a path the replay has to follow to stay in step.
    centroid_cells: u64,
    /// Distinct samples that are a corner of at least one surface cell.
    surface_cell_corners: u64,
    /// How many of those are exactly on the isosurface.
    equal_corners: u64,
    /// Surface cells whose **base** corner is exactly on the isosurface —
    /// M-316's own narrower census, reproduced to check this loader against the
    /// record.
    m316_equal_base_corners: u64,
}

/// Replay the march for provenance, reading the same values the crate read.
fn replay(values: &[f64], raw: &[u8], shape: &RuntimeShape3) -> Provenance {
    let size = shape.size();
    let samples = shape.element_count();
    let mut p = Provenance {
        indices: Vec::new(),
        cell: Vec::new(),
        edge_lo: Vec::new(),
        edge_axis: Vec::new(),
        surface_cells: 0,
        centroid_cells: 0,
        surface_cell_corners: 0,
        equal_corners: 0,
        m316_equal_base_corners: 0,
    };
    // The crate's cache: one slot per grid edge, keyed on the lower sample plus
    // the axis, so the key is the same whichever cell arrives first.
    let mut edge_vertex = vec![u32::MAX; samples * 3];
    let mut corner_seen = vec![0u64; samples.div_ceil(64)];
    let mut next = 0u32;

    for z in 0..size[2] - 1 {
        for y in 0..size[1] - 1 {
            for x in 0..size[0] - 1 {
                let base = [x, y, z];
                let mut case = 0u8;
                let mut sample = [0u32; 8];
                for (c, slot) in sample.iter_mut().enumerate() {
                    let o = corner_offset(c as u8);
                    let s = shape.linearize([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
                    *slot = s;
                    if is_inside(values[s as usize]) {
                        case |= 1 << c;
                    }
                }
                if case == 0 || case == u8::MAX {
                    continue;
                }
                p.surface_cells += 1;
                if raw[sample[0] as usize] == EQUAL {
                    p.m316_equal_base_corners += 1;
                }
                for s in sample {
                    corner_seen[s as usize / 64] |= 1 << (s % 64);
                }

                let entry = CASES[case as usize];
                if entry.count == 0 {
                    continue;
                }
                // Cycle centroids are allocated before any triangle of the cell
                // and are cell-local, so no grid edge names them.
                let mut centroid = [u32::MAX; table::MAX_CENTROIDS];
                if entry.centroids > 0 {
                    p.centroid_cells += 1;
                }
                for slot in centroid.iter_mut().take(entry.centroids as usize) {
                    *slot = next;
                    next += 1;
                    p.edge_lo.push(u32::MAX);
                    p.edge_axis.push(3);
                }

                let cell = shape.linearize(base);
                for tri in &entry.triangles[..entry.count as usize] {
                    for &code in tri {
                        let index = if table::is_centroid(code) {
                            centroid[(code - table::CENTROID_BASE) as usize]
                        } else {
                            let axis = EDGE_AXIS[code as usize];
                            let lo = sample[EDGE_CORNERS[code as usize][0] as usize];
                            let key = lo as usize * 3 + axis as usize;
                            if edge_vertex[key] == u32::MAX {
                                edge_vertex[key] = next;
                                next += 1;
                                p.edge_lo.push(lo);
                                p.edge_axis.push(axis);
                            }
                            edge_vertex[key]
                        };
                        p.indices.push(index);
                    }
                    p.cell.push(cell);
                }
            }
        }
    }

    for (w, word) in corner_seen.iter().enumerate() {
        let mut bits = *word;
        while bits != 0 {
            let s = w * 64 + bits.trailing_zeros() as usize;
            bits &= bits - 1;
            p.surface_cell_corners += 1;
            if raw[s] == EQUAL {
                p.equal_corners += 1;
            }
        }
    }

    p
}

// ─── the degeneracy census, exactly as the ledger counted it ────────────────

/// Bit pattern of a position, so coincidence is bitwise rather than approximate.
fn bits(p: [f64; 3]) -> [u64; 3] {
    [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()]
}

/// Doubled area, as a vector: `(b − a) × (c − a)`.
fn doubled_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

/// The crate's normal rule, restated over the same public [`Sdf::gradient`]:
/// `marching_cubes::unit_gradient` is private, and a snapped vertex has to get
/// the normal the crate would have given it or the arms differ in two things.
fn unit_gradient<S: Sdf<Scalar = f64>>(field: &S, p: [f64; 3]) -> [f64; 3] {
    let g = field.gradient(p);
    let inv = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt().recip();
    [g[0] * inv, g[1] * inv, g[2] * inv]
}

/// Is this triangle degenerate?
///
/// **No epsilon.** Two of the three indices are the same vertex, or two of the
/// three positions are bit-identical, or the doubled-area cross product is
/// exactly the zero vector. All three are exact — integers, bit patterns, and a
/// float against zero rather than a tuned threshold — so the count is the same
/// on every machine, which a relative-area threshold is not obliged to be.
fn is_degenerate(idx: [u32; 3], p: [[f64; 3]; 3]) -> bool {
    idx[0] == idx[1]
        || idx[1] == idx[2]
        || idx[2] == idx[0]
        || bits(p[0]) == bits(p[1])
        || bits(p[1]) == bits(p[2])
        || bits(p[2]) == bits(p[0])
        || doubled_area(p[0], p[1], p[2]) == [0.0; 3]
}

/// The three positions of a triangle.
fn triangle_positions(mesh: &MeshBuffer<f64>, idx: [u32; 3]) -> [[f64; 3]; 3] {
    [
        mesh.positions[idx[0] as usize],
        mesh.positions[idx[1] as usize],
        mesh.positions[idx[2] as usize],
    ]
}

/// Degenerate triangles, and how many are attributable to an `=` corner.
#[derive(Clone, Copy, Default)]
struct Degeneracy {
    /// Triangles [`is_degenerate`] accepts.
    total: u64,
    /// Of those, how many sit in a cell with a corner exactly on the isosurface.
    from_equal_corners: u64,
    /// Two of three indices the same vertex.
    repeated_index: u64,
    /// Two of three positions bit-identical.
    coincident_position: u64,
    /// Exactly zero area with neither of the above.
    zero_area_only: u64,
    /// What `isomesh::validate` counts at its own `area <= 1e-6 * h^2`, recorded
    /// beside the exact count so the two definitions can be read against each
    /// other rather than confused.
    validator_epsilon: u64,
}

/// Does this cell have a corner whose raw value is exactly the isovalue?
fn cell_has_equal_corner(shape: &RuntimeShape3, raw: &[u8], cell: u32) -> bool {
    let base = shape.delinearize(cell);
    (0..8u8).any(|c| {
        let o = corner_offset(c);
        let s = shape.linearize([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
        raw[s as usize] == EQUAL
    })
}

/// Count degenerate triangles and tag each one with its own cell's corners.
///
/// `cell[t]` is the base sample of the cell that emitted triangle `t`, so the
/// tag is read off that cell rather than inferred from a correlation.
fn census(mesh: &MeshBuffer<f64>, cell: &[u32], raw: &[u8], shape: &RuntimeShape3) -> Degeneracy {
    let mut d = Degeneracy::default();
    for (t, tri) in mesh.indices.chunks_exact(3).enumerate() {
        let idx = [tri[0], tri[1], tri[2]];
        let p = triangle_positions(mesh, idx);
        let repeated = idx[0] == idx[1] || idx[1] == idx[2] || idx[2] == idx[0];
        let coincident =
            bits(p[0]) == bits(p[1]) || bits(p[1]) == bits(p[2]) || bits(p[2]) == bits(p[0]);
        let flat = doubled_area(p[0], p[1], p[2]) == [0.0; 3];
        if !(repeated || coincident || flat) {
            continue;
        }
        d.total += 1;
        if repeated {
            d.repeated_index += 1;
        }
        if coincident {
            d.coincident_position += 1;
        }
        if flat && !repeated && !coincident {
            d.zero_area_only += 1;
        }
        if cell_has_equal_corner(shape, raw, cell[t]) {
            d.from_equal_corners += 1;
        }
    }
    d
}

/// Mark every vertex of a degenerate triangle that no collapse group claims.
///
/// Returns how many vertices that is. A degenerate triangle has zero area and
/// draws nothing, so the *site* is what can be highlighted; two of its three
/// corners coincide and are therefore in a group, and this is the third one.
fn mark_degenerate_vertices(mesh: &MeshBuffer<f64>, class: &mut [u8]) -> u64 {
    let mut unreached = 0;
    for tri in mesh.indices.chunks_exact(3) {
        let idx = [tri[0], tri[1], tri[2]];
        if !is_degenerate(idx, triangle_positions(mesh, idx)) {
            continue;
        }
        for v in idx {
            if class[v as usize] == CLASS_PLAIN {
                class[v as usize] = CLASS_DEGENERATE;
                unreached += 1;
            }
        }
    }
    unreached
}

// ─── the collapse ───────────────────────────────────────────────────────────

/// Disjoint-set over vertex indices, path halving, unioned to the **lower**
/// root.
///
/// `validate::Dsu` is private and this needs the same property it has: the
/// result depends only on which unions were requested, never on the order they
/// arrived in, so the component count below is a pure function of the mesh.
struct Dsu(Vec<u32>);

impl Dsu {
    /// `n` singletons.
    fn new(n: usize) -> Self {
        Self((0..n as u32).collect())
    }

    /// The root of `x`, halving the path on the way.
    fn find(&mut self, mut x: u32) -> u32 {
        while self.0[x as usize] != x {
            let parent = self.0[x as usize];
            self.0[x as usize] = self.0[parent as usize];
            x = self.0[x as usize];
        }
        x
    }

    /// Join two sets, keeping the lower root.
    fn union(&mut self, a: u32, b: u32) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
        self.0[hi as usize] = lo;
    }
}

/// What the ternary label did, and the mesh it produced.
struct Ternary {
    /// The repaired mesh.
    mesh: MeshBuffer<f64>,
    /// Per surviving triangle: its cell's base sample, as in [`Provenance`].
    cell: Vec<u32>,
    /// Vertices moved onto an `=` corner.
    snapped_vertices: u64,
    /// Groups with more than one member — the ones that actually merge.
    collapsed_groups: u64,
    /// Vertices the merge removed.
    vertices_removed: u64,
    /// Triangles that then named one vertex twice.
    triangles_dropped: u64,
    /// The largest move the snap actually made. Expected to be exactly zero.
    max_snap: f64,
    /// Vertices whose position is bit-identical to an earlier vertex's and which
    /// the `=` label does not reach.
    unexplained_coincident: u64,
    /// Collapsed groups whose members did **not** all already share triangles.
    pinch_groups: u64,
    /// Summed `components − 1` over those groups: how many separate pieces the
    /// collapse welded together in total.
    pinch_excess_components: u64,
    /// Per baseline vertex: [`CLASS_PLAIN`], [`CLASS_SAFE`] or [`CLASS_PINCH`].
    class: Vec<u8>,
    /// Per baseline vertex: its index in [`Self::mesh`], or `u32::MAX` when it
    /// was merged into another.
    new_index: Vec<u32>,
    /// World positions of the corners where the collapse welds pieces that
    /// shared no triangle.
    pinch_sites: Vec<[f64; 3]>,
}

/// The world position of a grid corner. With [`ORIGIN`] zero and [`CELL`] one it
/// is the integer grid coordinate exactly.
fn corner_position(shape: &RuntimeShape3, s: u32) -> [f64; 3] {
    let c = shape.delinearize(s);
    [
        ORIGIN[0] + CELL * f64::from(c[0]),
        ORIGIN[1] + CELL * f64::from(c[1]),
        ORIGIN[2] + CELL * f64::from(c[2]),
    ]
}

/// Label, snap, collapse — P-53's treatment arm, replayed.
///
/// The label is the pre-pass: a ternary sign per corner, read off the raw bytes.
/// Every vertex on an edge incident to an `=` corner snaps to that corner, the
/// vertices sharing one become the lowest of their indices, and a triangle that
/// then names a vertex twice is dropped. Nothing else moves.
fn ternary(
    field: &SampledField<'_, f64, RuntimeShape3>,
    base: &MeshBuffer<f64>,
    prov: &Provenance,
    raw: &[u8],
    shape: &RuntimeShape3,
) -> Ternary {
    let size = shape.size();
    let stride = [1, size[0], size[0] * size[1]];
    let n = base.positions.len();

    // ── the label, and what each vertex snaps to ────────────────────────────
    let mut target = vec![u32::MAX; n];
    let mut pairs: Vec<(u32, u32)> = Vec::new();
    for (v, slot) in target.iter_mut().enumerate() {
        let lo = prov.edge_lo[v];
        if lo == u32::MAX {
            continue;
        }
        let hi = lo + stride[prov.edge_axis[v] as usize];
        // A cut edge has one endpoint strictly inside, and `=` is on the outside
        // of `is_inside`, so at most one of the two can be `=`.
        let s = if raw[lo as usize] == EQUAL {
            lo
        } else if raw[hi as usize] == EQUAL {
            hi
        } else {
            continue;
        };
        *slot = s;
        pairs.push((s, v as u32));
    }
    // Sorted, then scanned in runs: no map, no iteration order to leak. The key
    // ends in the vertex index, so no two entries compare equal.
    pairs.sort_unstable();

    let mut remap: Vec<u32> = (0..n as u32).collect();
    let mut collapsed_groups = 0u64;
    let mut groups = 0u64;
    let mut i = 0;
    while i < pairs.len() {
        let mut j = i + 1;
        while j < pairs.len() && pairs[j].0 == pairs[i].0 {
            j += 1;
        }
        groups += 1;
        if j - i > 1 {
            collapsed_groups += 1;
        }
        // The lowest-indexed member represents the group, which is the tie-break
        // `isomesh::weld` documents.
        let rep = pairs[i].1;
        for &(_, v) in &pairs[i..j] {
            remap[v as usize] = rep;
        }
        i = j;
    }

    // ── the snap, measured ──────────────────────────────────────────────────
    let mut max_snap = 0.0f64;
    for (v, &s) in target.iter().enumerate() {
        if s == u32::MAX {
            continue;
        }
        let to = corner_position(shape, s);
        let from = base.positions[v];
        for a in 0..3 {
            max_snap = max_snap.max((to[a] - from[a]).abs());
        }
    }

    // ── the collapsed mesh ──────────────────────────────────────────────────
    let mut mesh = MeshBuffer::<f64>::new();
    let mut new_index = vec![u32::MAX; n];
    for (v, slot) in new_index.iter_mut().enumerate() {
        if remap[v] as usize != v {
            continue;
        }
        let (position, normal) = if target[v] == u32::MAX {
            (base.positions[v], base.normals[v])
        } else {
            let p = corner_position(shape, target[v]);
            (p, unit_gradient(field, p))
        };
        *slot = mesh.vertex(position, normal);
    }

    let mut cell = Vec::new();
    let mut triangles_dropped = 0u64;
    for (t, tri) in base.indices.chunks_exact(3).enumerate() {
        let a = new_index[remap[tri[0] as usize] as usize];
        let b = new_index[remap[tri[1] as usize] as usize];
        let c = new_index[remap[tri[2] as usize] as usize];
        if a == b || b == c || c == a {
            triangles_dropped += 1;
            continue;
        }
        mesh.triangle(a, b, c);
        cell.push(prov.cell[t]);
    }

    // ── coincidence the label does not explain ──────────────────────────────
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_unstable_by(|&a, &b| {
        let (p, q) = (base.positions[a as usize], base.positions[b as usize]);
        // `total_cmp`, so a NaN coordinate sorts into view instead of vanishing
        // through a partial comparison.
        p[0].total_cmp(&q[0])
            .then(p[1].total_cmp(&q[1]))
            .then(p[2].total_cmp(&q[2]))
            .then(a.cmp(&b))
    });
    let mut unexplained_coincident = 0u64;
    let mut i = 0;
    while i < order.len() {
        let key = bits(base.positions[order[i] as usize]);
        let mut j = i + 1;
        while j < order.len() && bits(base.positions[order[j] as usize]) == key {
            j += 1;
        }
        if j - i > 1 {
            for &v in &order[i..j] {
                if target[v as usize] == u32::MAX {
                    unexplained_coincident += 1;
                }
            }
        }
        i = j;
    }

    // ── which collapses join pieces, rather than flatten folds ──────────────
    //
    // Two vertices of one group that already share a triangle are two corners of
    // a triangle about to be dropped: the fold flattens and nothing moves. Two
    // that share no triangle are on different pieces of the surface, and
    // identifying their point is a pinch. This is the difference between the
    // topology holding and the topology moving, so it is counted rather than
    // argued.
    let mut dsu = Dsu::new(n);
    for tri in base.indices.chunks_exact(3) {
        for (a, b) in [(0, 1), (1, 2), (2, 0)] {
            let (u, w) = (tri[a], tri[b]);
            if u != w && remap[u as usize] == remap[w as usize] {
                dsu.union(u, w);
            }
        }
    }
    let mut class = vec![CLASS_PLAIN; n];
    let mut pinch_groups = 0u64;
    let mut pinch_excess_components = 0u64;
    let mut pinch_sites = Vec::new();
    let mut g = 0;
    while g < pairs.len() {
        let mut j = g + 1;
        while j < pairs.len() && pairs[j].0 == pairs[g].0 {
            j += 1;
        }
        if j - g > 1 {
            let mut roots: Vec<u32> = pairs[g..j].iter().map(|&(_, v)| dsu.find(v)).collect();
            roots.sort_unstable();
            roots.dedup();
            let pinch = roots.len() > 1;
            if pinch {
                pinch_groups += 1;
                pinch_excess_components += roots.len() as u64 - 1;
                pinch_sites.push(corner_position(shape, pairs[g].0));
            }
            let paint = if pinch { CLASS_PINCH } else { CLASS_SAFE };
            for &(_, v) in &pairs[g..j] {
                class[v as usize] = paint;
            }
        }
        g = j;
    }

    Ternary {
        mesh,
        cell,
        snapped_vertices: pairs.len() as u64,
        collapsed_groups,
        vertices_removed: pairs.len() as u64 - groups,
        triangles_dropped,
        max_snap,
        unexplained_coincident,
        pinch_groups,
        pinch_excess_components,
        class,
        new_index,
        pinch_sites,
    }
}

// ─── the committed ledger ───────────────────────────────────────────────────

/// P-53's committed rows, read once at startup.
struct Ledger {
    /// The header line, split on commas.
    header: Vec<String>,
    /// Every data row, split on commas.
    rows: Vec<Vec<String>>,
}

impl Ledger {
    /// Read `docs/experiments/p-53.csv`. Comment lines start with `#`; the first
    /// line that does not is the header.
    fn read(path: &Path) -> Result<Self, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let mut header: Option<Vec<String>> = None;
        let mut rows = Vec::new();
        for line in text.lines() {
            let line = line.trim_end();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let cells: Vec<String> = line.split(',').map(str::to_string).collect();
            match header {
                None => header = Some(cells),
                Some(_) => rows.push(cells),
            }
        }
        let header = header.ok_or_else(|| format!("{}: no header line", path.display()))?;
        if rows.is_empty() {
            return Err(format!("{}: no data rows", path.display()));
        }
        Ok(Self { header, rows })
    }

    /// Where a named column sits.
    fn column(&self, name: &str) -> Option<usize> {
        self.header.iter().position(|h| h == name)
    }

    /// The row for one volume and one label rule at [`ISOVALUE_LABEL`].
    fn row(&self, volume: &str, rule: &str) -> Option<&[String]> {
        let v = self.column("volume")?;
        let i = self.column("isovalue")?;
        let r = self.column("label_rule")?;
        self.rows
            .iter()
            .find(|row| {
                row.get(v).is_some_and(|c| c == volume)
                    && row.get(i).is_some_and(|c| c == ISOVALUE_LABEL)
                    && row.get(r).is_some_and(|c| c == rule)
            })
            .map(Vec::as_slice)
    }

    /// One cell of a row.
    fn cell<'a>(&self, row: &'a [String], name: &str) -> Option<&'a str> {
        Some(row.get(self.column(name)?)?.as_str())
    }
}

/// How many ledger columns were compared, and how many reproduced.
#[derive(Default, Clone, Copy)]
struct Tally {
    /// Columns compared.
    checked: u32,
    /// Columns whose live value equals the committed one.
    matched: u32,
}

impl Tally {
    /// Compare one column, parsing the committed cell the same way the live
    /// value is typed.
    ///
    /// A mismatch is **loud and does not take the window down with it**: this is
    /// a demo a stranger runs, and a disagreement with the ledger is a finding
    /// that has to be readable rather than an abort.
    fn check<T>(&mut self, ledger: &Ledger, row: &[String], who: &str, name: &str, live: T)
    where
        T: std::str::FromStr + PartialEq + std::fmt::Display,
        <T as std::str::FromStr>::Err: std::fmt::Display,
    {
        self.checked += 1;
        let Some(cell) = ledger.cell(row, name) else {
            error!("p-53.csv has no column {name}, so {who} cannot be checked against it");
            return;
        };
        match cell.parse::<T>() {
            Ok(want) if want == live => self.matched += 1,
            Ok(want) => error!(
                "{who} {name}: this run says {live}, p-53.csv says {want}. \
                 The ledger is the record; this run is the disagreement."
            ),
            Err(e) => error!("{who} {name}: ledger cell {cell:?} does not parse ({e})"),
        }
    }

    /// Whether every compared column reproduced.
    fn holds(self) -> bool {
        self.checked > 0 && self.checked == self.matched
    }
}

// ─── one volume, both arms ──────────────────────────────────────────────────

/// The numbers one arm of one volume produced.
#[derive(Clone, Copy, Default)]
struct Arm {
    /// Vertices in the mesh.
    vertices: usize,
    /// Triangles in the mesh.
    triangles: usize,
    /// Degenerate triangles, by the exact test.
    degenerate: u64,
    /// Of those, how many came from a cell with an `=` corner.
    degenerate_from_equal: u64,
    /// `V − E + F` over the valid subset.
    chi: i64,
    /// Edges used by three or more faces.
    non_manifold_edges: u64,
    /// Edges used by exactly one face.
    boundary_edges: u64,
}

/// Everything one volume contributes, computed once before the window opens.
struct VolumeData {
    /// The ledger's name for it.
    name: &'static str,
    /// Samples per axis.
    samples: u32,
    /// What the crate's own Marching Cubes took, for scale. Gates nothing.
    extract_ms: f64,
    /// Samples exactly on the isosurface, among surface-cell corners.
    equal_corners: u64,
    /// The baseline arm: the crate's own march.
    base: Arm,
    /// The repaired arm: the same mesh under the ternary label.
    repaired: Arm,
    /// Groups with more than one member.
    collapsed_groups: u64,
    /// Of those, how many weld pieces that shared no triangle.
    pinch_groups: u64,
    /// How many separate pieces those welds joined.
    pinch_excess: u64,
    /// The largest move the snap made. Exactly zero, which is the point.
    max_snap: f64,
    /// Vertices of a degenerate triangle that no collapse group claims.
    unreached_degenerate_vertices: u64,
    /// Where the pinches are, in world space.
    pinch_sites: Vec<Vec3>,
    /// Half-length of a pinch cross on this volume.
    marker_half: f32,
    /// Wide camera target and distance.
    wide: (Vec3, f32),
    /// The sampled grid, for the `G` box.
    domain: (Vec3, Vec3),
    /// How the ledger comparison came out.
    tally: Tally,
}

impl VolumeData {
    /// One line, for the log and the HUD, carrying this volume's headline row.
    fn headline(&self) -> String {
        format!(
            "{:<6} {} -> {}  chi {} -> {}  {}/{} pinch, {} welded",
            self.name,
            self.base.degenerate,
            self.repaired.degenerate,
            self.base.chi,
            self.repaired.chi,
            self.pinch_groups,
            self.collapsed_groups,
            self.pinch_excess,
        )
    }
}

/// Read one volume, march it, replay it, repair it, and check every number
/// against the ledger.
///
/// Returns the numbers and the two display meshes. Both arms are built here
/// because `bonsai`'s march alone is ~500 ms and there are two validator passes
/// on top; doing any of it per frame would make the demo a slideshow.
fn build_volume(dir: &Path, v: &Volume, ledger: &Ledger) -> (VolumeData, [Mesh; 2]) {
    let raw = match read_u8(&dir.join(v.file), v.n) {
        Ok(raw) => raw,
        Err(e) => {
            error!("{e}");
            std::process::exit(1);
        }
    };
    let shape = match RuntimeShape3::new([v.n; 3]) {
        Ok(shape) => shape,
        Err(e) => {
            error!("{}: {e}", v.file);
            std::process::exit(1);
        }
    };
    // `iso - value`, so a dense voxel is negative and the crate's sign
    // convention holds unchanged. `benches/volumes.rs` forms it the same way.
    let values: Vec<f64> = raw.iter().map(|b| ISOVALUE - f64::from(*b)).collect();
    let field = match SampledField::new(&values, &shape, ORIGIN, CELL) {
        Ok(field) => field,
        Err(e) => {
            error!("{}: {e}", v.file);
            std::process::exit(1);
        }
    };
    let cfg = match ValidateConfig::from_cell_size(CELL) {
        Ok(cfg) => cfg,
        Err(e) => {
            error!("{e}");
            std::process::exit(1);
        }
    };

    // ── the baseline: the crate's own Marching Cubes ────────────────────────
    let mut base = MeshBuffer::<f64>::new();
    let t0 = Instant::now();
    if let Err(e) = MarchingCubes::<f64>::new().extract(&field, &shape, ORIGIN, CELL, &mut base) {
        error!("{} at iso {ISOVALUE_LABEL}: {e}", v.file);
        std::process::exit(1);
    }
    let extract_ms = t0.elapsed().as_secs_f64() * 1e3;

    let prov = replay(&values, &raw, &shape);
    let replay_matches = prov.indices == base.indices && prov.edge_lo.len() == base.positions.len();
    if !replay_matches {
        error!(
            "{}: the replay is not the crate's march ({} vs {} indices, {} vs {} vertices) \
             -- every attribution below would be a guess",
            v.short,
            prov.indices.len(),
            base.indices.len(),
            prov.edge_lo.len(),
            base.positions.len()
        );
    }

    let mut base_deg = census(&base, &prov.cell, &raw, &shape);
    let base_report = validate_indexed(&base.positions, &base.indices, &cfg);
    base_deg.validator_epsilon = base_report.degenerate_triangles;
    let base_hash = mesh_hash(&base);

    // ── the repaired arm ───────────────────────────────────────────────────
    let tern = ternary(&field, &base, &prov, &raw, &shape);
    let mut tern_deg = census(&tern.mesh, &tern.cell, &raw, &shape);
    let tern_report = validate_indexed(&tern.mesh.positions, &tern.mesh.indices, &cfg);
    tern_deg.validator_epsilon = tern_report.degenerate_triangles;
    let tern_hash = mesh_hash(&tern.mesh);

    // ── the paint ──────────────────────────────────────────────────────────
    let mut base_class = tern.class.clone();
    let unreached = mark_degenerate_vertices(&base, &mut base_class);
    let mut repaired_class = vec![CLASS_PLAIN; tern.mesh.positions.len()];
    for (vertex, &new) in tern.new_index.iter().enumerate() {
        if new != u32::MAX && tern.class[vertex] == CLASS_PINCH {
            repaired_class[new as usize] = CLASS_PINCH;
        }
    }

    // ── the ledger comparison ──────────────────────────────────────────────
    let mut tally = Tally::default();
    let shared = |t: &mut Tally, row: &[String], who: &str| {
        t.check(ledger, row, who, "samples_per_axis", u64::from(v.n));
        t.check(ledger, row, who, "cells", u64::from(v.n - 1).pow(3));
        t.check(ledger, row, who, "surface_cells", prov.surface_cells);
        t.check(
            ledger,
            row,
            who,
            "surface_cell_corners",
            prov.surface_cell_corners,
        );
        t.check(ledger, row, who, "equal_corners", prov.equal_corners);
        t.check(
            ledger,
            row,
            who,
            "m316_equal_base_corners",
            prov.m316_equal_base_corners,
        );
        t.check(ledger, row, who, "centroid_cells", prov.centroid_cells);
        t.check(ledger, row, who, "replay_matches_crate", replay_matches);
        t.check(
            ledger,
            row,
            who,
            "unexplained_coincident_vertices",
            tern.unexplained_coincident,
        );
    };
    let arm_columns = |t: &mut Tally,
                       row: &[String],
                       who: &str,
                       mesh: &MeshBuffer<f64>,
                       deg: &Degeneracy,
                       report: &isomesh::validate::MeshReport,
                       hash: u64| {
        t.check(ledger, row, who, "triangles", mesh.triangle_count() as u64);
        t.check(ledger, row, who, "vertices", mesh.vertex_count() as u64);
        t.check(ledger, row, who, "degenerate_triangles", deg.total);
        t.check(
            ledger,
            row,
            who,
            "degenerate_from_equal_corners",
            deg.from_equal_corners,
        );
        t.check(
            ledger,
            row,
            who,
            "degenerate_repeated_index",
            deg.repeated_index,
        );
        t.check(
            ledger,
            row,
            who,
            "degenerate_coincident_position",
            deg.coincident_position,
        );
        t.check(
            ledger,
            row,
            who,
            "degenerate_zero_area_only",
            deg.zero_area_only,
        );
        t.check(
            ledger,
            row,
            who,
            "degenerate_validator_epsilon",
            deg.validator_epsilon,
        );
        t.check(
            ledger,
            row,
            who,
            "euler_characteristic",
            report.euler_characteristic,
        );
        t.check(
            ledger,
            row,
            who,
            "non_manifold_edges",
            report.non_manifold_edges,
        );
        t.check(ledger, row, who, "boundary_edges", report.boundary_edges);
        t.check(ledger, row, who, "mesh_hash", hash);
    };

    match ledger.row(v.short, "binary") {
        Some(row) => {
            let who = format!("{}/binary", v.short);
            shared(&mut tally, row, &who);
            arm_columns(
                &mut tally,
                row,
                &who,
                &base,
                &base_deg,
                &base_report,
                base_hash,
            );
        }
        None => error!("p-53.csv has no {}/{ISOVALUE_LABEL}/binary row", v.short),
    }
    match ledger.row(v.short, "ternary") {
        Some(row) => {
            let who = format!("{}/ternary", v.short);
            shared(&mut tally, row, &who);
            arm_columns(
                &mut tally,
                row,
                &who,
                &tern.mesh,
                &tern_deg,
                &tern_report,
                tern_hash,
            );
            tally.check(ledger, row, &who, "snapped_vertices", tern.snapped_vertices);
            tally.check(ledger, row, &who, "collapsed_groups", tern.collapsed_groups);
            tally.check(ledger, row, &who, "vertices_removed", tern.vertices_removed);
            tally.check(
                ledger,
                row,
                &who,
                "triangles_dropped",
                tern.triangles_dropped,
            );
            tally.check(ledger, row, &who, "max_snap_distance", tern.max_snap);
            tally.check(ledger, row, &who, "pinch_groups", tern.pinch_groups);
            tally.check(
                ledger,
                row,
                &who,
                "pinch_excess_components",
                tern.pinch_excess_components,
            );
        }
        None => error!("p-53.csv has no {}/{ISOVALUE_LABEL}/ternary row", v.short),
    }

    // ── framing ────────────────────────────────────────────────────────────
    let (aabb_min, aabb_max) = aabb(&base.positions);
    let centre = (aabb_min + aabb_max) * 0.5;
    let half = (aabb_max - aabb_min) * 0.5;
    let wide = (centre, WIDE_RADIUS * half.max_element().max(1.0));
    let sites: Vec<Vec3> = tern
        .pinch_sites
        .iter()
        .map(|p| Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32))
        .collect();
    let extent = (aabb_max - aabb_min).max_element().max(1.0);

    let data = VolumeData {
        name: v.short,
        samples: v.n,
        extract_ms,
        equal_corners: prov.equal_corners,
        base: Arm {
            vertices: base.vertex_count(),
            triangles: base.triangle_count(),
            degenerate: base_deg.total,
            degenerate_from_equal: base_deg.from_equal_corners,
            chi: base_report.euler_characteristic,
            non_manifold_edges: base_report.non_manifold_edges,
            boundary_edges: base_report.boundary_edges,
        },
        repaired: Arm {
            vertices: tern.mesh.vertex_count(),
            triangles: tern.mesh.triangle_count(),
            degenerate: tern_deg.total,
            degenerate_from_equal: tern_deg.from_equal_corners,
            chi: tern_report.euler_characteristic,
            non_manifold_edges: tern_report.non_manifold_edges,
            boundary_edges: tern_report.boundary_edges,
        },
        collapsed_groups: tern.collapsed_groups,
        pinch_groups: tern.pinch_groups,
        pinch_excess: tern.pinch_excess_components,
        max_snap: tern.max_snap,
        unreached_degenerate_vertices: unreached,
        pinch_sites: sites,
        marker_half: extent * MARKER_FRACTION,
        wide,
        domain: (Vec3::ZERO, Vec3::splat((v.n - 1) as f32)),
        tally,
    };

    info!(
        "E-310 self-check, {} at {}^3 iso {ISOVALUE_LABEL}: {}",
        data.name,
        data.samples,
        data.headline()
    );
    info!(
        "    degenerate {} -> {} ({} of them from an = corner), nm_e {} -> {}, bnd {} -> {}, \
         max snap {:e}, snapped {} vertices, dropped {} triangles, extract {extract_ms:.3} ms",
        data.base.degenerate,
        data.repaired.degenerate,
        data.base.degenerate_from_equal,
        data.base.non_manifold_edges,
        data.repaired.non_manifold_edges,
        data.base.boundary_edges,
        data.repaired.boundary_edges,
        data.max_snap,
        tern.snapped_vertices,
        tern.triangles_dropped,
    );
    info!(
        "    mesh spans {:?} to {:?}, so the wide view sits {:.1} back from {:?}; \
         {} pinch sites",
        aabb_min,
        aabb_max,
        data.wide.1,
        data.wide.0,
        data.pinch_sites.len(),
    );
    if tally.holds() {
        info!(
            "    p-53.csv: all {} columns reproduce, mesh_hash {} / {} included",
            tally.checked, base_hash, tern_hash
        );
    } else {
        error!(
            "    p-53.csv: {} of {} columns reproduce -- SEE THE ERRORS ABOVE. The committed \
             ledger is the record; this run is the disagreement, and that is the finding.",
            tally.matched, tally.checked
        );
    }

    let meshes = [
        to_mesh(&base, &base_class),
        to_mesh(&tern.mesh, &repaired_class),
    ];
    (data, meshes)
}

/// The axis-aligned bounds of a position array.
fn aabb(positions: &[[f64; 3]]) -> (Vec3, Vec3) {
    let mut min = Vec3::splat(f32::INFINITY);
    let mut max = Vec3::splat(f32::NEG_INFINITY);
    for p in positions {
        let q = Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32);
        min = min.min(q);
        max = max.max(q);
    }
    if positions.is_empty() {
        (Vec3::ZERO, Vec3::ONE)
    } else {
        (min, max)
    }
}

/// The `f64` extraction as a Bevy mesh, painted per vertex.
///
/// Cast rather than re-extracted in `f32`: the numbers on the HUD are P-53's and
/// they are `f64` numbers, so the mesh the picture is drawn from has to be the
/// one they were computed on. Positions are integers below 2^24, so the cast is
/// exact.
fn to_mesh(buffer: &MeshBuffer<f64>, class: &[u8]) -> Mesh {
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
    let colours = builder.colors_mut();
    colours.reserve(class.len());
    colours.extend(class.iter().copied().map(class_colour));
    for t in buffer.indices.chunks_exact(3) {
        let (Some(a), Some(b), Some(c)) = (t.first(), t.get(1), t.get(2)) else {
            continue;
        };
        builder.triangle(*a, *b, *c);
    }
    builder.into_mesh()
}

// ─── state ──────────────────────────────────────────────────────────────────

/// Everything both volumes produced, keyed by the digit that selects them.
#[derive(Resource)]
struct Prepared(Vec<VolumeData>);

/// The display meshes, waiting for [`Assets<Mesh>`] to exist.
#[derive(Resource)]
struct Pending(Vec<[Mesh; 2]>);

/// Per volume: the baseline and the repaired mesh asset.
#[derive(Resource)]
struct Arms(Vec<[Handle<Mesh>; 2]>);

/// A volume pinned by `ISOMESH_FIELD`, which overrides the capture's stepping.
#[derive(Resource)]
struct Pinned(Option<usize>);

/// What is on screen this frame.
#[derive(Resource, Default)]
struct View {
    /// Index into [`VOLUMES`].
    volume: usize,
    /// Whether the repaired arm is showing.
    repaired: bool,
}

/// The pinch crosses draw in front of the surface, so they need their own depth
/// bias: a marker lying on a surface z-fights and reads as intermittent, which
/// is indistinguishable from the defect being intermittent.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct PinchGizmos;

// ─── app ────────────────────────────────────────────────────────────────────

fn main() {
    // Checked before a single Bevy plugin is built, and with `println!` rather
    // than `info!`, because there is no log subscriber yet and because a clean
    // clone with no network is not an error -- M-006. The volumes are not
    // committed; the fetch script is.
    let dir = volume_dir();
    let missing: Vec<&str> = VOLUMES
        .iter()
        .filter(|v| !dir.join(v.file).is_file())
        .map(|v| v.file)
        .collect();
    if !missing.is_empty() {
        println!(
            "E-310 needs the CT volumes P-53 was measured on, and {} of {} \
             {} not in {}:\n  {}\n\nFetch them, then run this again:\n  \
             ./scripts/fetch_volumes.sh\n\nThe files are ~17 MB total and are \
             verified against the publisher's own SHA-512; see \
             docs/measurements/volumes/PROVENANCE.md.",
            missing.len(),
            VOLUMES.len(),
            if missing.len() == 1 { "is" } else { "are" },
            dir.display(),
            missing.join("\n  "),
        );
        return;
    }

    let mut app = App::new();
    app.add_plugins(DefaultPlugins.set(WindowPlugin {
        primary_window: Some(Window {
            title: "isomesh - E-310 pinch repair".into(),
            ..default()
        }),
        ..default()
    }))
    .add_plugins(CommonPlugin)
    .init_gizmo_group::<PinchGizmos>();

    // Here rather than in `setup`, so the self-check is logged before the window
    // opens and after `add_plugins` has installed the log subscriber -- and so
    // every system below can take the result as a plain `Res`.
    let ledger = match Ledger::read(&ledger_path()) {
        Ok(ledger) => ledger,
        Err(e) => {
            error!(
                "{e} -- this example's whole claim is that its numbers are the committed \
                 ones, so there is nothing to show without the ledger to check them against"
            );
            std::process::exit(1);
        }
    };
    let built: Vec<(VolumeData, [Mesh; 2])> = VOLUMES
        .iter()
        .map(|v| build_volume(&dir, v, &ledger))
        .collect();
    let (prepared, pending): (Vec<VolumeData>, Vec<[Mesh; 2]>) = built.into_iter().unzip();
    let pinned = pinned_volume();
    info!(
        "E-310: the repair is exact on both volumes and free on only one. {} | {}",
        prepared[0].headline(),
        prepared[VOLUMES.len() - 1].headline(),
    );

    app.insert_resource(Prepared(prepared))
        .insert_resource(Pending(pending))
        .insert_resource(Pinned(pinned))
        // `ISOMESH_FIELD` picks the starting volume whether or not a capture is
        // running. The harness's contract is that anything a capture depends on
        // is reachable from the environment, and a still is a one-frame capture.
        .insert_resource(View {
            volume: pinned.unwrap_or(0),
            ..default()
        })
        .add_systems(Startup, setup)
        // `PreUpdate`, for the reason `aperture_gate` documents at length: the
        // harness's `update_hud` and `capture_sequence` both run in `Update` and
        // give no ordering, so a HUD written there renders a frame-old
        // `DemoStats` beside a current mesh. Here the arm, the camera, the
        // crosses and the HUD in any one frame are all the same volume.
        .add_systems(
            PreUpdate,
            (controls, apply_view, draw_pinches, report)
                .chain()
                .after(bevy::input::InputSystems),
        )
        .run();
}

/// The volume `ISOMESH_FIELD` asks for, if it asks for one.
fn pinned_volume() -> Option<usize> {
    let raw = std::env::var("ISOMESH_FIELD").ok()?;
    match raw.trim().parse::<usize>() {
        Ok(index) if index < VOLUME_COUNT => Some(index),
        _ => {
            error!("ISOMESH_FIELD={raw} is not one of 0..{VOLUME_COUNT}");
            None
        }
    }
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut config: ResMut<GizmoConfigStore>,
    mut pending: ResMut<Pending>,
) {
    let (pinch, _) = config.config_mut::<PinchGizmos>();
    pinch.line.width = 2.6;
    // Negative, so a cross reads through the tree it is inside. All 516 of them
    // showing at once is the picture; the ones facing away are not decoration.
    pinch.depth_bias = -1.0;

    let arms: Vec<[Handle<Mesh>; 2]> = std::mem::take(&mut pending.0)
        .into_iter()
        .map(|[base, repaired]| [meshes.add(base), meshes.add(repaired)])
        .collect();
    commands.insert_resource(Arms(arms));
    commands.remove_resource::<Pending>();

    // White base colour, because `StandardMaterial` multiplies it by the vertex
    // colour and the vertex colour is the entire classification. Double-sided:
    // a pinch weld is a non-manifold junction and the two sheets meeting there
    // are seen from both faces.
    let material = materials.add(StandardMaterial {
        base_color: Color::WHITE,
        perceptual_roughness: 0.75,
        metallic: 0.03,
        double_sided: true,
        cull_mode: None,
        ..default()
    });

    // `Mesh3d::default()` names no asset, so nothing is uploaded until
    // `apply_view` picks an arm on the first frame.
    commands.spawn((Mesh3d::default(), MeshMaterial3d(material), DemoMesh));

    // Spawned rather than assumed: `draw_domain` queries for it, and without one
    // the `G` toggle silently does nothing. Filled in by the first view.
    commands.spawn(DemoDomain {
        min: Vec3::ZERO,
        max: Vec3::ONE,
    });
}

/// Frames a capture runs for.
///
/// Read from the environment rather than from [`Capture`], which keeps its
/// length private, because pacing the stages off the capture is what stops a
/// six-frame smoke test and an eighty-frame clip from both being a still.
fn capture_frames() -> u32 {
    std::env::var("ISOMESH_CAPTURE_FRAMES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(60)
        .max(1)
}

/// Decide the volume and the arm for this frame.
///
/// Under capture both come off the captured-frame counter — four equal stages,
/// each volume's two arms — so a clip of any length carries the whole
/// comparison. With `ISOMESH_FIELD` pinning a volume the two arms split the clip
/// in half instead.
fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    pinned: Res<Pinned>,
    mut view: ResMut<View>,
) {
    if capture.is_active() {
        let total = capture_frames();
        let stage = (u64::from(capture.taken) * u64::from(STAGES) / u64::from(total)) as u32;
        let stage = stage.min(STAGES - 1);
        view.repaired = stage % 2 == 1;
        view.volume = pinned
            .0
            .unwrap_or((stage / 2) as usize)
            .min(VOLUME_COUNT - 1);
        return;
    }

    if keys.just_pressed(KeyCode::Digit1) {
        view.volume = 0;
    }
    if keys.just_pressed(KeyCode::Digit2) {
        view.volume = (VOLUME_COUNT - 1).min(1);
    }
    if keys.just_pressed(KeyCode::KeyB) {
        view.repaired = !view.repaired;
    }
}

/// Swap the arm, the domain box and the camera — only when the view changed.
fn apply_view(
    view: Res<View>,
    prepared: Res<Prepared>,
    arms: Res<Arms>,
    mut surface: Query<&mut Mesh3d, With<DemoMesh>>,
    mut domain: Query<&mut DemoDomain>,
    mut camera: Query<&mut OrbitCamera>,
    mut last: Local<Option<(usize, bool)>>,
) {
    let key = (view.volume, view.repaired);
    if *last == Some(key) {
        return;
    }
    let Some(data) = prepared.0.get(view.volume) else {
        return;
    };
    let Some(pair) = arms.0.get(view.volume) else {
        return;
    };
    *last = Some(key);

    for mut mesh in &mut surface {
        mesh.0 = pair[usize::from(view.repaired)].clone();
    }

    for mut d in &mut domain {
        d.min = data.domain.0;
        d.max = data.domain.1;
    }

    for mut orbit in &mut camera {
        orbit.focus = data.wide.0;
        orbit.radius = data.wide.1;
    }
}

/// The pinch crosses: one per group whose members shared no triangle.
///
/// Drawn in **both** arms, at the same places, because the positions are the
/// whole point — `max_snap_distance` is exactly zero, so the repair moved
/// nothing and these are where two sheets meet before and after.
fn draw_pinches(
    view: Res<View>,
    prepared: Res<Prepared>,
    flags: Res<ViewFlags>,
    mut gizmos: Gizmos<PinchGizmos>,
) {
    if !flags.grid && !flags.hud {
        // `nogrid,nohud` is the "just the geometry" view, and a cross is
        // annotation.
        return;
    }
    let Some(data) = prepared.0.get(view.volume) else {
        return;
    };
    let colour = Color::srgb(COLOUR_PINCH[0], COLOUR_PINCH[1], COLOUR_PINCH[2]);
    let h = data.marker_half;
    for &site in &data.pinch_sites {
        gizmos.line(site - Vec3::X * h, site + Vec3::X * h, colour);
        gizmos.line(site - Vec3::Y * h, site + Vec3::Y * h, colour);
        gizmos.line(site - Vec3::Z * h, site + Vec3::Z * h, colour);
    }
}

/// The HUD. The numbers are the demo.
fn report(view: Res<View>, prepared: Res<Prepared>, mut stats: ResMut<DemoStats>) {
    let Some(data) = prepared.0.get(view.volume) else {
        return;
    };
    let arm = if view.repaired {
        &data.repaired
    } else {
        &data.base
    };

    // Every line below is kept inside 72 characters, the longest being the
    // `bonsai` check row. At the harness's 13px font that is about 570 logical
    // pixels, so nothing wraps in the 1024-wide capture the module docs
    // recommend -- and a wrapped line in a GIF reads as a bug.
    stats.title = format!(
        "E-310  pinch repair - {} {}^3  iso {ISOVALUE_LABEL}  [1-2] volume  [B] arm",
        data.name, data.samples,
    );
    stats.vertices = arm.vertices;
    stats.triangles = arm.triangles;
    stats.extract_ms = data.extract_ms;

    let tally = data.tally;
    stats.extra = vec![
        format!(
            "arm       {}",
            if view.repaired {
                "REPAIRED - the ternary label; [B] shows the baseline"
            } else {
                "BASELINE - the crate as it ships; [B] shows the repair"
            },
        ),
        format!(
            "degen     {} -> {} triangles, all {} from an = corner",
            data.base.degenerate, data.repaired.degenerate, data.base.degenerate_from_equal,
        ),
        format!(
            "groups    {} collapsed, {} pinches, {} components welded",
            data.collapsed_groups, data.pinch_groups, data.pinch_excess,
        ),
        format!(
            "snap      max_snap_distance {:e} - exactly zero, nothing moved",
            data.max_snap
        ),
        String::from("topology            chi   nm_edges   boundary"),
        format!(
            "  baseline  {:>11} {:>10} {:>10}",
            data.base.chi, data.base.non_manifold_edges, data.base.boundary_edges,
        ),
        format!(
            "  repaired  {:>11} {:>10} {:>10}",
            data.repaired.chi, data.repaired.non_manifold_edges, data.repaired.boundary_edges,
        ),
        format!(
            "colour    magenta pinch, amber fold, red sliver ({} of them)",
            data.unreached_degenerate_vertices
        ),
        format!(
            "corners   {} samples exactly on {ISOVALUE_LABEL}   ledger {}",
            data.equal_corners,
            if tally.holds() {
                format!("p-53.csv {}/{} ok", tally.matched, tally.checked)
            } else {
                format!(
                    "p-53.csv {}/{} -- SEE THE LOG",
                    tally.matched, tally.checked
                )
            },
        ),
        format!("check     {}", prepared.0[0].headline()),
        format!("check     {}", prepared.0[VOLUME_COUNT - 1].headline()),
        String::from("finding   ship the precondition, not the repair"),
    ];
}
