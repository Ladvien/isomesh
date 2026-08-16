//! **A-022 — the dual path's ambiguous face, with the source finally read.**
//!
//! ```bash
//! cargo bench --bench a022_decider
//! ```
//!
//! Writes `docs/measurements/a022-decider.csv`.
//!
//! # The paper was not paywalled, it was on the author's own page
//!
//! A-022 sat blocked on *"the source is paywalled and not in the corpus"* —
//! Schaefer, Ju & Warren, **Manifold Dual Contouring**, IEEE TVCG 13(3) 2007,
//! `10.1109/TVCG.2007.1012`, with `paper_download` reporting *"No open-access
//! PDF found"*. It is at `cs.wustl.edu/~taoju/research/dualsimp_tvcg.pdf`,
//! which is Tao Ju's own publications page and is the exact filename the ticket
//! already named. V-31's rule, in a new costume.
//!
//! # What it says, and it settles the question the ticket could not
//!
//! §3, *Contouring on a Uniform Grid*, in full on the point:
//!
//! > *"One of the limitations of DC is that it allows no more than one vertex
//! > within each grid cell. On a uniform grid, **DC leads to nonmanifold
//! > vertices and edges for all of the ambiguous sign configurations** in the
//! > original MC algorithm.*
//! >
//! > *To combat this effect, Nielson's modification allows multiple vertices to
//! > be placed in a single cell. In particular, **Nielson associates one vertex
//! > with each cycle of a modified MC table [26]**. Since each cycle consists of
//! > a list of edges on the cubic cell, each vertex is associated with a set of
//! > edges, and **each edge is associated with exactly one vertex**. … this
//! > surface is always a manifold because the original MC algorithm always
//! > constructs a manifold and the dual preserves the topology of the surface."*
//!
//! **Reference [26] is Nielson & Hamann, *The Asymptotic Decider*.** So the
//! criterion is neither "face-based" nor "component-based" as A-022 framed the
//! choice: it is **one vertex per cycle of a table whose ambiguous faces have
//! already been resolved by the decider**. The face ambiguity is settled
//! upstream, inside the table the cycles are read from, and the dual walk needs
//! no rule of its own.
//!
//! Two consequences, and the second is what this bench measures.
//!
//! **A-022's acceptance was unreachable as written.** It asked for M-276's 314
//! non-manifold edges to go to zero *under Surface Nets and Dual Contouring*.
//! Those are one vertex per cell, and the paper says in as many words that this
//! is nonmanifold for **all** ambiguous configurations. The 314 is the
//! literature's own prediction, not a defect to remove.
//!
//! **This crate's `ManifoldDualContouring` defaults to `FaceAmbiguity::Separate`,
//! which is not the table the paper specifies.** It reads its cycles from
//! `segment_links` on the *unmodified* table. M-276's residual **53** was
//! therefore measured against a construction the source does not describe. So:
//! census both rules, on the fields M-276 used, with Surface Nets and Dual
//! Contouring beside them as the control that says what one vertex per cell
//! costs.
//!
//! **Marching Cubes itself is in the census, because the paper's argument rests
//! on it.** Its manifoldness claim is *"because the original MC algorithm always
//! constructs a manifold and the dual preserves the topology of the surface"* —
//! and ✗15 in this repository is *"Marching Cubes is unconditionally manifold",
//! falsified*. If the premise fails here then the conclusion is expected to, and
//! the two counts should be related rather than independent.
//!
//! No chunking and no weld, matching A-021 — the weld can create a non-manifold
//! edge (M-226) and would put a second mechanism in the count.

use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use isomesh::dual_contouring::DualContouring;
use isomesh::extractor::Extractor;
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::marching_cubes::{FaceAmbiguity, MarchingCubes};
use isomesh::surface_nets::SurfaceNets;
use isomesh::validate::{ValidateConfig, validate_indexed};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf, for_each_reference_field};

/// `f64`, matching the validity suite.
type Scalar = f64;

/// Samples per axis. 49 is the resolution M-276 counted its 314 at.
const RESOLUTIONS: [u32; 3] = [33, 49, 65];

/// One measured configuration.
struct Row {
    field: &'static str,
    extractor: &'static str,
    samples: u32,
    non_manifold_edges: u64,
    non_manifold_vertices: u64,
    boundary_edges: u64,
    triangles: usize,
}

fn measure<E, F>(
    extractor: &mut E,
    name: &'static str,
    field: &F,
    field_name: &'static str,
    samples: u32,
) -> Row
where
    E: Extractor<Scalar>,
    F: Sdf<Scalar = Scalar> + isomesh::fields::ReferenceField,
{
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("the fixture fits u32");
    let mut mesh = MeshBuffer::<Scalar>::new();
    extractor
        .extract_into(field, &shape, lo, cell_size, &mut mesh)
        .expect("extraction");
    let cfg = ValidateConfig::from_cell_size(cell_size).expect("a positive cell size");
    let report = validate_indexed(&mesh.positions, &mesh.indices, &cfg);
    Row {
        field: field_name,
        extractor: name,
        samples,
        non_manifold_edges: report.non_manifold_edges,
        non_manifold_vertices: report.non_manifold_vertices,
        boundary_edges: report.boundary_edges,
        triangles: mesh.triangle_count(),
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!(
        "{:<16} {:<28} {:>5} {:>10} {:>10} {:>10} {:>10}",
        "field", "extractor", "n", "nm edges", "nm verts", "boundary", "triangles"
    );

    let mut rows: Vec<Row> = Vec::new();
    for samples in RESOLUTIONS {
        for_each_reference_field!(Scalar, |name, field| {
            let mut separate = ManifoldDualContouring::<Scalar>::new();
            let mut decider = ManifoldDualContouring::<Scalar>::new();
            decider.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
            // Marching Cubes itself, because the paper's argument rests on it.
            let mut mc_decider = MarchingCubes::<Scalar>::new();
            mc_decider.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
            let measured = [
                measure(
                    &mut MarchingCubes::<Scalar>::new(),
                    "marching_cubes",
                    &field,
                    name,
                    samples,
                ),
                measure(&mut mc_decider, "mc+decider", &field, name, samples),
                measure(&mut separate, "mdc+separate", &field, name, samples),
                measure(&mut decider, "mdc+decider", &field, name, samples),
                measure(
                    &mut SurfaceNets::<Scalar>::new(),
                    "surface_nets",
                    &field,
                    name,
                    samples,
                ),
                measure(
                    &mut DualContouring::<Scalar>::new(),
                    "dual_contouring",
                    &field,
                    name,
                    samples,
                ),
            ];
            for r in measured {
                println!(
                    "{:<16} {:<28} {:>5} {:>10} {:>10} {:>10} {:>10}",
                    r.field,
                    r.extractor,
                    r.samples,
                    r.non_manifold_edges,
                    r.non_manifold_vertices,
                    r.boundary_edges,
                    r.triangles
                );
                rows.push(r);
            }
        });
    }

    let total = |name: &str| -> (u64, u64) {
        rows.iter()
            .filter(|r| r.extractor == name)
            .fold((0, 0), |(e, v), r| {
                (e + r.non_manifold_edges, v + r.non_manifold_vertices)
            })
    };
    println!("\ntotals over every field and resolution — non-manifold edges / vertices:");
    for name in [
        "marching_cubes",
        "mc+decider",
        "mdc+separate",
        "mdc+decider",
        "surface_nets",
        "dual_contouring",
    ] {
        let (e, v) = total(name);
        println!("  {name:<20} {e:>8} / {v:<8}");
    }

    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/measurements");
    fs::create_dir_all(&dir).expect("create docs/measurements");
    let mut csv = String::from("# A-022: does the decider-modified table remove MDC's residue?\n");
    let _ = writeln!(
        csv,
        "field,extractor,samples,non_manifold_edges,non_manifold_vertices,boundary_edges,triangles"
    );
    for r in &rows {
        let _ = writeln!(
            csv,
            "{},{},{},{},{},{},{}",
            r.field,
            r.extractor,
            r.samples,
            r.non_manifold_edges,
            r.non_manifold_vertices,
            r.boundary_edges,
            r.triangles
        );
    }
    let path = dir.join("a022-decider.csv");
    fs::write(&path, csv).expect("write csv");
    println!("\n{} rows → {}", rows.len(), path.display());
}
