//! **P-8 — does gating the weld on the link condition remove non-manifoldness?**
//!
//! Ticket: R-001. Pre-registered at R-000 before a line of this file existed;
//! `isomesh::experiment!("P-8")` is a compile error otherwise.
//!
//! ```bash
//! cargo bench --bench experiment_p8
//! ```
//!
//! Writes `docs/experiments/p-8.csv`.
//!
//! # The condition, and why it is not vacuous here
//!
//! Merging two vertices `u` and `v` that are **not** adjacent is safe when their
//! links are disjoint: if they share a neighbour `w`, the merged vertex and `w`
//! are joined by an edge carrying faces from both sides, which is how an edge
//! acquires three faces.
//!
//! At a plain chunk seam the condition is trivially satisfied — the two copies
//! live in different chunks and share nothing. **It stops being trivial where
//! more than two chunks meet**, and where a grid sample lands exactly on the
//! isosurface (M-48): several vertices coincide, they are merged one at a time,
//! and each merge changes the links the next one is tested against. Dey, Fan &
//! Wang's decomposition of a `k`-way merge into `k − 1` pairwise merges *in the
//! intermediate complex* is the reason the order matters, and is what P-9 is
//! about.
//!
//! # Both arms share the candidate pairing, so only the gate varies
//!
//! The ungated arm is `Welder` itself. The gated arm **replays `Welder`'s own
//! `remap`** — so the two arms agree on which vertices are candidates to merge
//! and differ only in whether a candidate is accepted. Re-deriving the
//! candidates with a second lattice would have made the comparison a comparison
//! of broadphases.

mod common;

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::extractor::Extractor;
use isomesh::validate::{ValidateConfig, validate_indexed};
use isomesh::weld::{Welder, epsilon_for};
use isomesh::{MeshBuffer, Sdf};

/// Cells per chunk.
///
/// **18, not 8, and the first version's 8 was a fixture bug (M-274).** A 2×2×2
/// block of 8-cell chunks at `h = 4/35` spans 1.83 units from `-2.0`, which does
/// not reach the origin — so it clipped a corner off every field and no chunk
/// edge crossed the surface anywhere interesting. 36 cells spans 4.11, which
/// covers the `[-2, 2]` domain with the block **centred**, so the four-chunk
/// edges and the eight-chunk corner sit where the surface is.
const CELLS: u32 = 18;

/// Cell size. `4/35` rather than a power of two, deliberately: M-32 measured
/// that a seam is bit-exact **only** at a power of two, so a nicer spacing would
/// hand both arms a mesh with no duplicates to reason about.
const CELL_SIZE: f64 = 4.0 / 35.0;

/// Block origin, placed so the 2×2×2 block is **centred on the field**.
///
/// The eight-chunk corner then sits at the origin and the four-chunk edges cross
/// the surface, which is where a bucket can hold more than two vertices at all.
const ORIGIN: f64 = -(2.0 * CELLS as f64) * CELL_SIZE / 2.0;

/// Vertices adjacent to each vertex, from the triangle list.
///
/// The link of a vertex in a triangle mesh is the cycle of opposite edges; two
/// links intersect exactly when the vertices share a neighbour, so the adjacency
/// sets are what the test needs and the cycles are not.
fn adjacency(indices: &[u32], n: usize) -> Vec<BTreeSet<u32>> {
    let mut adj = vec![BTreeSet::new(); n];
    for t in indices.chunks_exact(3) {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            if a != b {
                adj[a as usize].insert(b);
                adj[b as usize].insert(a);
            }
        }
    }
    adj
}

/// What one arm produced.
struct Arm {
    non_manifold_edges: u64,
    non_manifold_vertices: u64,
    vertices: usize,
    ms: f64,
}

/// The ungated weld: `Welder`, unchanged.
fn ungated(mesh: &MeshBuffer<f64>, cfg: &ValidateConfig) -> (Arm, Vec<u32>) {
    let mut work = mesh.clone();
    let mut welder = Welder::<f64>::new();
    let start = Instant::now();
    welder
        .weld(&mut work, epsilon_for(CELL_SIZE))
        .expect("valid epsilon");
    let ms = start.elapsed().as_secs_f64() * 1e3;
    let remap = welder.remap().to_vec();
    let report = validate_indexed(&work.positions, &work.indices, cfg);
    (
        Arm {
            non_manifold_edges: report.non_manifold_edges,
            non_manifold_vertices: report.non_manifold_vertices,
            vertices: work.positions.len(),
            ms,
        },
        remap,
    )
}

/// The gated weld: the same candidates, each accepted only if the links are
/// disjoint **in the complex as it stands when that merge is considered**.
///
/// Returns the arm and how many candidate merges were refused.
fn gated(mesh: &MeshBuffer<f64>, remap: &[u32], cfg: &ValidateConfig) -> (Arm, usize) {
    let start = Instant::now();
    let n = mesh.positions.len();

    // `Welder`'s remap sends every vertex to its group's output index, so
    // grouping by that value recovers the candidate sets without re-running a
    // broadphase.
    let mut groups: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for (v, &to) in remap.iter().enumerate() {
        groups.entry(to).or_default().push(v as u32);
    }

    // Start with every vertex its own class, then accept merges one at a time.
    let mut class: Vec<u32> = (0..n as u32).collect();
    let mut adj = adjacency(&mesh.indices, n);
    let mut rejected = 0usize;

    for members in groups.values() {
        if members.len() < 2 {
            continue;
        }
        // Lowest index is the representative, matching `Welder`'s own rule.
        let rep = members[0];
        for &v in &members[1..] {
            let (a, b) = (class[rep as usize], class[v as usize]);
            if a == b {
                continue;
            }
            // The link condition, on the current complex. `a` and `b` are not
            // adjacent (coincident vertices from different chunks never are),
            // so sharing any neighbour is exactly what makes the merge unsafe.
            let shared = adj[a as usize].intersection(&adj[b as usize]).count();
            if shared > 0 {
                rejected += 1;
                continue;
            }
            // Accept: fold `b` into `a`, in the intermediate complex, so the
            // next merge in this group is tested against the result.
            let moved: Vec<u32> = adj[b as usize].iter().copied().collect();
            for w in moved {
                adj[w as usize].remove(&b);
                adj[w as usize].insert(a);
                adj[a as usize].insert(w);
            }
            adj[b as usize].clear();
            adj[a as usize].remove(&b);
            for c in &mut class {
                if *c == b {
                    *c = a;
                }
            }
        }
    }

    // Realise the classes as a mesh, compacting indices.
    let mut out_of: BTreeMap<u32, u32> = BTreeMap::new();
    let mut out = MeshBuffer::<f64>::new();
    for v in 0..n as u32 {
        let c = class[v as usize];
        out_of.entry(c).or_insert_with(|| {
            let at = out.positions.len() as u32;
            out.positions.push(mesh.positions[c as usize]);
            out.normals.push(mesh.normals[c as usize]);
            at
        });
    }
    for t in mesh.indices.chunks_exact(3) {
        let a = out_of[&class[t[0] as usize]];
        let b = out_of[&class[t[1] as usize]];
        let c = out_of[&class[t[2] as usize]];
        // A triangle whose corners collapsed together has no area; dropping it
        // is what `Welder` does too, for the same reason.
        if a != b && b != c && a != c {
            out.indices.extend_from_slice(&[a, b, c]);
        }
    }
    let ms = start.elapsed().as_secs_f64() * 1e3;

    let report = validate_indexed(&out.positions, &out.indices, cfg);
    (
        Arm {
            non_manifold_edges: report.non_manifold_edges,
            non_manifold_vertices: report.non_manifold_vertices,
            vertices: out.positions.len(),
            ms,
        },
        rejected,
    )
}

/// Eight chunks in a 2×2×2 block, meshed independently and appended.
///
/// **A block, not a pair.** Two chunks share a plane and their seam vertices
/// have disjoint links by construction, which would make the gate vacuous.
/// Eight chunks share edges and a corner, which is where several vertices
/// coincide at once and where the order of a `k`-way merge can matter.
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

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-8");
    common::experiment::run(prereg, |run| {
        let layout = ChunkLayout::<f64>::new(CELLS, CELL_SIZE, [ORIGIN; 3]).expect("valid layout");
        let cfg = ValidateConfig::from_cell_size(CELL_SIZE).expect("valid cell size");

        println!(
            "{:<16} {:<22} {:>8} {:>9} {:>8} {:>9} {:>9} {:>7}",
            "field",
            "extractor",
            "nm_e_un",
            "nm_e_gate",
            "nm_v_un",
            "nm_v_gate",
            "rejected",
            "Δverts"
        );

        isomesh::for_each_reference_field!(f64, |name, field| {
            // Inline blocks, so no `return` in either of these (M-253).
            isomesh::for_each_extractor!(f64, |ename, extractor| {
                let mesh = eight_chunks(&field, &layout, &mut extractor);
                if !mesh.indices.is_empty() {
                    let (un, remap) = ungated(&mesh, &cfg);
                    let (ga, rejected) = gated(&mesh, &remap, &cfg);

                    println!(
                        "{name:<16} {ename:<22} {:>8} {:>9} {:>8} {:>9} {rejected:>9} {:>7}",
                        un.non_manifold_edges,
                        ga.non_manifold_edges,
                        un.non_manifold_vertices,
                        ga.non_manifold_vertices,
                        ga.vertices as i64 - un.vertices as i64
                    );

                    run.record(&[
                        ("field", name.to_string()),
                        ("extractor", ename.to_string()),
                        (
                            "non_manifold_edges_ungated",
                            un.non_manifold_edges.to_string(),
                        ),
                        (
                            "non_manifold_edges_gated",
                            ga.non_manifold_edges.to_string(),
                        ),
                        (
                            "non_manifold_vertices_ungated",
                            un.non_manifold_vertices.to_string(),
                        ),
                        (
                            "non_manifold_vertices_gated",
                            ga.non_manifold_vertices.to_string(),
                        ),
                        ("rejected_merges", rejected.to_string()),
                        (
                            "vertex_delta",
                            (ga.vertices as i64 - un.vertices as i64).to_string(),
                        ),
                        ("weld_ms_ungated", format!("{:.4}", un.ms)),
                        ("weld_ms_gated", format!("{:.4}", ga.ms)),
                    ]);
                }
            });
        });
    });
}
