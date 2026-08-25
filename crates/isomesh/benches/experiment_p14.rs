//! **P-14 — where does the dual extractors' residual non-manifoldness live?**
//!
//! Ticket: R-003, re-scoped by R-001. Pre-registered at R-000 in its own commit,
//! before this file existed.
//!
//! ```bash
//! cargo bench --bench experiment_p14
//! ```
//!
//! Writes `docs/experiments/p-14.csv`.
//!
//! # What R-001 left behind
//!
//! P-8 measured `noise_cavity` under `surface_nets` at **276 non-manifold edges
//! and 536 non-manifold vertices**, and showed the weld neither caused it nor
//! could remove it — the link-gated weld rejects one merge there and changes the
//! edge count by zero. So it is the extractor's own output. A count is not a
//! diagnosis, and this is the diagnosis.
//!
//! # The discriminator
//!
//! Surface Nets and Dual Contouring place **one vertex per cell**. Manifold Dual
//! Contouring exists precisely because a cell can contain more than one surface
//! component, and it splits those cells into several vertices. So MDC's own
//! output is a *census of multi-component cells*, measured by the crate rather
//! than asserted — and cross-referencing it against where the non-manifold
//! vertices actually are settles whether the one-vertex-per-cell rule is the
//! cause.
//!
//! # A single chunk, and no weld anywhere
//!
//! P-8 and P-9 used eight chunks because they were about welding. This is about
//! the extractor, so a chunk seam would only add a second explanation. One grid,
//! one extraction, straight into `validate_features`.

mod common;

use std::collections::BTreeMap;

use isomesh::extractor::Extractor;
use isomesh::fields::ReferenceField;
use isomesh::validate::{ValidateConfig, validate_features};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Samples per axis. Large enough that `noise_cavity` shows the defect at the
/// scale P-8 saw, small enough for a bench.
const SAMPLES: u32 = 49;

/// Which cell a dual vertex belongs to.
///
/// Both rules place their vertex strictly inside the cell it belongs to, so
/// flooring the offset recovers the cell without the extractor having to report
/// it. A vertex exactly on a cell boundary would be attributed to the higher
/// cell; that is a tie-break, not an approximation, and it is the same one on
/// both arms so it cannot bias the cross-reference.
fn cell_of(p: [f64; 3], origin: [f64; 3], h: f64) -> [i64; 3] {
    [
        ((p[0] - origin[0]) / h).floor() as i64,
        ((p[1] - origin[1]) / h).floor() as i64,
        ((p[2] - origin[2]) / h).floor() as i64,
    ]
}

/// Mesh one grid with one extractor.
fn mesh<E: Extractor<f64>>(
    field: &impl Sdf<Scalar = f64>,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    h: f64,
    extractor: &mut E,
) -> MeshBuffer<f64> {
    let mut out = MeshBuffer::<f64>::new();
    extractor
        .extract_into(field, shape, origin, h, &mut out)
        .expect("extraction");
    out
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-14");
    common::experiment::run(prereg, |run| {
        println!(
            "{:<16} {:<26} {:>7} {:>7} {:>10} {:>14} {:>7} {:>7}",
            "field",
            "extractor",
            "nm_v",
            "nm_e",
            "multi_cell",
            "in_single_cell",
            "worst_k",
            "worst_deg"
        );

        isomesh::for_each_reference_field!(f64, |name, field| {
            // Inline block, so no `return` in here (M-253).
            let (lo, hi) = field.domain();
            let h = (hi[0] - lo[0]) / f64::from(SAMPLES - 1);
            let shape = RuntimeShape3::new([SAMPLES; 3]).expect("valid shape");
            let cfg = ValidateConfig::from_cell_size(h).expect("valid cell size");

            // The census of multi-component cells, taken from Manifold Dual
            // Contouring's own output rather than re-derived.
            let mdc = mesh(
                &field,
                &shape,
                lo,
                h,
                &mut isomesh::manifold_dual_contouring::ManifoldDualContouring::<f64>::new(),
            );
            let mut per_cell: BTreeMap<[i64; 3], u32> = BTreeMap::new();
            for p in &mdc.positions {
                *per_cell.entry(cell_of(*p, lo, h)).or_insert(0) += 1;
            }
            let multi: std::collections::BTreeSet<[i64; 3]> = per_cell
                .iter()
                .filter(|&(_, &n)| n > 1)
                .map(|(c, _)| *c)
                .collect();

            for (ename, out) in [
                (
                    "surface_nets",
                    mesh(
                        &field,
                        &shape,
                        lo,
                        h,
                        &mut isomesh::surface_nets::SurfaceNets::<f64>::new(),
                    ),
                ),
                (
                    "dual_contouring",
                    mesh(
                        &field,
                        &shape,
                        lo,
                        h,
                        &mut isomesh::dual_contouring::DualContouring::<f64>::new(),
                    ),
                ),
                ("manifold_dual_contouring", mdc.clone()),
            ] {
                let (report, features) = validate_features(&out.positions, &out.indices, &cfg);

                // Of the offending vertices, how many are in a cell MDC did
                // **not** split? Those are the ones the hypothesis cannot
                // explain.
                let in_single = features
                    .vertices
                    .iter()
                    .filter(|&&v| !multi.contains(&cell_of(out.positions[v as usize], lo, h)))
                    .count();

                // How badly split the worst link is: an incident-face link with
                // `k` components is `k` cones sharing an apex, and `k = 2` is the
                // classic bowtie. Recomputed here rather than taken from the
                // report, which counts vertices and not components.
                let worst = features
                    .vertices
                    .iter()
                    .map(|&v| link_components(&out.indices, v))
                    .max()
                    .unwrap_or(0);

                // **A second shape of non-manifold vertex, and the first run
                // could not see it.** `worst_link_components == 1` on a vertex
                // `validate` flagged is not a contradiction: a link can be
                // *connected* and still not a simple cycle, where one link
                // vertex is reached by four link edges instead of two. That is a
                // **pinch**, not a bowtie, and it is a different defect with a
                // different cause. An extra column rather than a changed
                // registration -- `records` lists metrics that must be
                // reported, not the whole schema (M-273).
                let worst_degree = features
                    .vertices
                    .iter()
                    .map(|&v| worst_link_degree(&out.indices, v))
                    .max()
                    .unwrap_or(0);

                println!(
                    "{name:<16} {ename:<26} {:>7} {:>7} {:>10} {:>14} {worst:>7} \
                     {worst_degree:>7}",
                    report.non_manifold_vertices,
                    report.non_manifold_edges,
                    multi.len(),
                    in_single
                );
                run.record(&[
                    ("field", name.to_string()),
                    ("extractor", ename.to_string()),
                    (
                        "non_manifold_vertices",
                        report.non_manifold_vertices.to_string(),
                    ),
                    ("non_manifold_edges", report.non_manifold_edges.to_string()),
                    ("multi_vertex_cells", multi.len().to_string()),
                    ("nm_vertices_in_single_vertex_cells", in_single.to_string()),
                    ("worst_link_components", worst.to_string()),
                    ("worst_link_vertex_degree", worst_degree.to_string()),
                ]);
            }
        });
    });
}

/// Connected components of the incident-face link at `v`.
///
/// Two faces are in the same component when they share an edge **through `v`**.
/// Sharing only `v` itself is exactly the bowtie, which is why the adjacency is
/// on the opposite-edge endpoints rather than on the faces' vertex sets.
fn link_components(indices: &[u32], v: u32) -> usize {
    // The link is the set of opposite edges; components are its connected parts.
    let mut edges: Vec<[u32; 2]> = Vec::new();
    for t in indices.as_chunks::<3>().0 {
        if let Some(i) = t.iter().position(|&x| x == v) {
            let a = t[(i + 1) % 3];
            let b = t[(i + 2) % 3];
            edges.push([a, b]);
        }
    }
    if edges.is_empty() {
        return 0;
    }

    let mut seen = vec![false; edges.len()];
    let mut components = 0;
    for start in 0..edges.len() {
        if seen[start] {
            continue;
        }
        components += 1;
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(i) = stack.pop() {
            for j in 0..edges.len() {
                if seen[j] {
                    continue;
                }
                let touching = edges[i].iter().any(|a| edges[j].contains(a));
                if touching {
                    seen[j] = true;
                    stack.push(j);
                }
            }
        }
    }
    components
}

/// Highest degree of any vertex in the link at `v`.
///
/// A manifold vertex's link is a simple cycle, so every link vertex has degree
/// exactly two. Degree four is a **pinch**: two fans meeting at a single link
/// vertex rather than being wholly disjoint, which leaves the link connected —
/// so [`link_components`] reads one and misses it entirely.
fn worst_link_degree(indices: &[u32], v: u32) -> usize {
    let mut degree: BTreeMap<u32, usize> = BTreeMap::new();
    for t in indices.as_chunks::<3>().0 {
        if let Some(i) = t.iter().position(|&x| x == v) {
            *degree.entry(t[(i + 1) % 3]).or_insert(0) += 1;
            *degree.entry(t[(i + 2) % 3]).or_insert(0) += 1;
        }
    }
    degree.into_values().max().unwrap_or(0)
}
