//! **P-20 — does a weld key move anything the key does not name?**
//!
//! Ticket: R-010. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p20
//! ```
//!
//! Writes `docs/experiments/p-20.csv`.
//!
//! # What E×4 leaves
//!
//! The link-condition gate was measured on this exact fixture and reverted as
//! strictly worse: over 56 configurations it removed **at most 4** non-manifold
//! edges and added **up to 791** non-manifold vertices — `noise_cavity` +
//! subgrid went 301 → 1,092, and `sphere` + Marching Cubes went 0 → 96. The
//! mechanism is that a `k`-way coincidence is manifold only if **all `k`** merge,
//! and a *pairwise* test refuses one member of the set, leaving its
//! representative a bowtie. That is why the damage sat in the vertex column while
//! the edge column barely moved.
//!
//! Equality on a key is an **equivalence relation**, so it partitions each class
//! into complete sub-classes and cannot leave a lone representative behind. This
//! measures whether that reasoning survives contact with eight fields and every
//! extractor.
//!
//! # Three arms, and the first two are controls
//!
//! - **`none`** — the unconditional weld, `weld()`. The baseline every other row
//!   is read against.
//! - **`constant`** — every vertex keyed `0`. Must be **identical** to `none` on
//!   every metric. This is a control on the *plumbing*: if a constant key moves
//!   anything, the hook itself is doing something and no result from the third
//!   arm can be trusted. Registered as an arm rather than assumed, per M-279 —
//!   "the hook perturbs the weld" is exactly the rival that a varying-key-only
//!   experiment could not tell apart from "the key does what it says."
//! - **`normal`** — the key a Bevy consumer would actually build (B-014's case):
//!   the vertex normal quantised to a coarse lattice. This is the arm with
//!   something to find.
//!
//! # Why the normal is quantised before hashing
//!
//! Two normals that differ in the last bit must land on the same key or they will
//! not merge, and an exact float comparison makes that a coin toss decided by
//! rounding. Quantising to `1/16` of a unit is coarse enough that a seam's two
//! copies of one vertex agree, and fine enough that a cube corner's three face
//! normals do not.

mod common;

use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::extractor::Extractor;
use isomesh::validate::{ValidateConfig, validate_indexed};
use isomesh::weld::{Welder, epsilon_for};
use isomesh::{MeshBuffer, Sdf};

/// Cells per chunk. P-8's, and M-274's correction with it: 18 centres the 2×2×2
/// block on the `[-2, 2]` domain so the chunk seams cross the surface.
const CELLS: u32 = 18;

/// Cell size. `4/35` rather than a power of two, deliberately — M-32 measured
/// that a seam is bit-exact *only* at a power of two, so a nicer spacing would
/// hand every arm a mesh with no duplicates to reason about.
const CELL_SIZE: f64 = 4.0 / 35.0;

/// Block origin.
const ORIGIN: f64 = -2.0;

/// Quantisation of the normal key, in units of a unit normal's components.
///
/// Coarse enough that a seam's two copies of one vertex agree; fine enough that
/// a crease does not collapse. See the module docs.
const NORMAL_QUANTUM: f64 = 16.0;

/// The 2×2×2 block of chunks, joined into one buffer with its seam duplicates
/// intact — which is the point, since those duplicates are what a weld merges.
fn eight_chunks<E: Extractor<f64>>(
    field: &impl Sdf<Scalar = f64>,
    layout: &ChunkLayout<f64>,
    extractor: &mut E,
) -> MeshBuffer<f64> {
    let shape = layout.sample_shape().expect("valid shape");
    let mut joined = MeshBuffer::<f64>::new();
    for z in 0..2 {
        for y in 0..2 {
            for x in 0..2 {
                let id = ChunkId::new([x, y, z]);
                let mut piece = MeshBuffer::<f64>::new();
                extractor
                    .extract_into(
                        field,
                        &shape,
                        layout.sample_origin(id),
                        layout.cell_size(),
                        &mut piece,
                    )
                    .expect("extraction");
                joined.append(&piece).expect("the meshes fit u32");
            }
        }
    }
    joined
}

/// One `u64` per vertex, from the quantised normal.
fn normal_keys(mesh: &MeshBuffer<f64>) -> Vec<u64> {
    mesh.normals
        .iter()
        .map(|n| {
            // Round rather than truncate, so a component either side of zero does
            // not land in two different buckets for the same direction.
            let q = |v: f64| (v * NORMAL_QUANTUM).round() as i64 as u64;
            // Three 21-bit fields is ample: a quantised unit component is in
            // [-16, 16], so it fits far inside.
            q(n[0]) & 0x1F_FFFF | (q(n[1]) & 0x1F_FFFF) << 21 | (q(n[2]) & 0x1F_FFFF) << 42
        })
        .collect()
}

/// A unit cube as 6 quads, **24 vertices**, four per face, each carrying its own
/// face normal.
///
/// **The field sweep cannot test the interesting half of H and this is why.** At
/// a chunk seam the two copies of a vertex come from the same field evaluated at
/// the same point, so their normals agree to far inside the quantum and the key
/// splits nothing — measured, 0 splits over all 50 field × extractor pairs. A
/// key that never fires is a control that cannot discriminate, and it reads
/// exactly like a real negative (M-279).
///
/// Here every one of the 8 corners is shared by **three** vertices with three
/// different normals, so the unconditional weld collapses 24 → 8 and the normal
/// key must keep all 24. That is also B-014's acceptance case, in the crate that
/// owns the mechanism.
fn creased_cube() -> MeshBuffer<f64> {
    let mut mesh = MeshBuffer::<f64>::new();
    // (origin, edge u, edge v, face normal), per face of the unit cube.
    type Face = ([f64; 3], [f64; 3], [f64; 3], [f64; 3]);
    let faces: [Face; 6] = [
        (
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0],
        ),
        (
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ),
        (
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
        ),
        (
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
        ),
        (
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
        ),
        (
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ),
    ];
    for (o, u, v, n) in faces {
        let base = mesh.positions.len() as u32;
        for (su, sv) in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
            mesh.positions.push([
                o[0] + u[0] * su + v[0] * sv,
                o[1] + u[1] * su + v[1] * sv,
                o[2] + u[2] * su + v[2] * sv,
            ]);
            mesh.normals.push(n);
        }
        mesh.indices
            .extend_from_slice(&[base, base + 1, base + 2, base, base + 2, base + 3]);
    }
    mesh
}

struct Arm {
    vertices_after: usize,
    non_manifold_edges: u64,
    non_manifold_vertices: u64,
    boundary_edges: u64,
}

/// Weld a copy of `mesh` with the given keys and validate the result.
fn arm(mesh: &MeshBuffer<f64>, keys: &[u64], cfg: &ValidateConfig) -> Arm {
    let mut work = mesh.clone();
    Welder::default()
        .weld_split_by(&mut work, epsilon_for(CELL_SIZE), keys)
        .expect("weld");
    let report = validate_indexed(&work.positions, &work.indices, cfg);
    Arm {
        vertices_after: work.positions.len(),
        non_manifold_edges: report.non_manifold_edges,
        non_manifold_vertices: report.non_manifold_vertices,
        boundary_edges: report.boundary_edges,
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-20");
    common::experiment::run(prereg, |run| {
        let layout = ChunkLayout::<f64>::new(CELLS, CELL_SIZE, [ORIGIN; 3]).expect("valid layout");
        let cfg = ValidateConfig::from_cell_size(CELL_SIZE).expect("valid cell size");

        println!(
            "{:<16} {:<22} {:<9} {:>9} {:>7} {:>7} {:>7} {:>7}",
            "field", "extractor", "key", "verts", "splits", "nm_e", "nm_v", "bnd_e"
        );

        // Violations of the two falsifiers, counted rather than eyeballed.
        let mut plumbing_moved = 0usize;
        let mut vertices_over_splits = 0usize;

        isomesh::for_each_reference_field!(f64, |name, field| {
            // Inline blocks, so no `return` in either of these (M-253).
            isomesh::for_each_extractor!(f64, |ename, extractor| {
                let mesh = eight_chunks(&field, &layout, &mut extractor);
                if !mesh.indices.is_empty() {
                    let constant = vec![0u64; mesh.positions.len()];
                    let varying = normal_keys(&mesh);

                    let none = arm(&mesh, &[], &cfg);
                    let flat = arm(&mesh, &constant, &cfg);
                    let split = arm(&mesh, &varying, &cfg);

                    // Falsifier one: a constant key must move nothing at all.
                    if flat.vertices_after != none.vertices_after
                        || flat.non_manifold_edges != none.non_manifold_edges
                        || flat.non_manifold_vertices != none.non_manifold_vertices
                        || flat.boundary_edges != none.boundary_edges
                    {
                        plumbing_moved += 1;
                    }

                    // Falsifier two: E*4's shape reappearing. The vertex count
                    // rising is what the key NAMES; non-manifold vertices rising
                    // by more than that would be a partial refusal within a class.
                    let splits = split.vertices_after as i64 - none.vertices_after as i64;
                    let nm_v_rise =
                        split.non_manifold_vertices as i64 - none.non_manifold_vertices as i64;
                    if nm_v_rise > splits.max(0) {
                        vertices_over_splits += 1;
                    }

                    for (label, a, s) in [
                        ("none", &none, 0i64),
                        (
                            "constant",
                            &flat,
                            flat.vertices_after as i64 - none.vertices_after as i64,
                        ),
                        ("normal", &split, splits),
                    ] {
                        println!(
                            "{name:<16} {ename:<22} {label:<9} {:>9} {s:>7} {:>7} {:>7} {:>7}",
                            a.vertices_after,
                            a.non_manifold_edges,
                            a.non_manifold_vertices,
                            a.boundary_edges
                        );
                        run.record(&[
                            ("field", name.to_string()),
                            ("extractor", ename.to_string()),
                            ("key", label.to_string()),
                            ("vertices_after", a.vertices_after.to_string()),
                            ("splits", s.to_string()),
                            ("non_manifold_edges", a.non_manifold_edges.to_string()),
                            ("non_manifold_vertices", a.non_manifold_vertices.to_string()),
                            ("boundary_edges", a.boundary_edges.to_string()),
                        ]);
                    }
                }
            });
        });

        // The crease fixture: the only rows where the key has anything to do.
        {
            let mesh = creased_cube();
            let constant = vec![0u64; mesh.positions.len()];
            let varying = normal_keys(&mesh);
            let none = arm(&mesh, &[], &cfg);
            let flat = arm(&mesh, &constant, &cfg);
            let split = arm(&mesh, &varying, &cfg);

            if flat.vertices_after != none.vertices_after
                || flat.non_manifold_edges != none.non_manifold_edges
                || flat.non_manifold_vertices != none.non_manifold_vertices
                || flat.boundary_edges != none.boundary_edges
            {
                plumbing_moved += 1;
            }
            let splits = split.vertices_after as i64 - none.vertices_after as i64;
            let nm_v_rise = split.non_manifold_vertices as i64 - none.non_manifold_vertices as i64;
            if nm_v_rise > splits.max(0) {
                vertices_over_splits += 1;
            }

            println!(
                "splits on the crease fixture: {splits} \
                 (0 would mean the varying arm discriminates nowhere at all)"
            );

            for (label, a, s) in [
                ("none", &none, 0i64),
                (
                    "constant",
                    &flat,
                    flat.vertices_after as i64 - none.vertices_after as i64,
                ),
                ("normal", &split, splits),
            ] {
                println!(
                    "{:<16} {:<22} {label:<9} {:>9} {s:>7} {:>7} {:>7} {:>7}",
                    "creased_cube",
                    "hand-built",
                    a.vertices_after,
                    a.non_manifold_edges,
                    a.non_manifold_vertices,
                    a.boundary_edges
                );
                run.record(&[
                    ("field", "creased_cube".to_string()),
                    ("extractor", "hand-built".to_string()),
                    ("key", label.to_string()),
                    ("vertices_after", a.vertices_after.to_string()),
                    ("splits", s.to_string()),
                    ("non_manifold_edges", a.non_manifold_edges.to_string()),
                    ("non_manifold_vertices", a.non_manifold_vertices.to_string()),
                    ("boundary_edges", a.boundary_edges.to_string()),
                ]);
            }
        }

        println!();
        println!("constant-key rows that moved a metric: {plumbing_moved} (falsifier: any)");
        println!(
            "rows where non-manifold vertices rose by more than the split count: \
             {vertices_over_splits} (falsifier: any)"
        );
    });
}
