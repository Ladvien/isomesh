//! **P-21 — does a freshly extracted mesh seal what the field seals?**
//!
//! Ticket: R-024. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p21
//! ```
//!
//! Writes `docs/experiments/p-21.csv`.
//!
//! # The question, and why nothing else in the crate answers it
//!
//! Every other validity metric here judges a mesh against **itself** —
//! manifoldness, orientation, Euler characteristic — or against the field's
//! **geometry**, in `validate::accuracy`. None asks whether the mesh partitions
//! *space* the way the field's sign does, and neither claim implies the other: a
//! mesh can be closed, manifold, correctly wound and Hausdorff-close while
//! sealing a passage the field leaves open, or opening one it seals.
//!
//! The instrument is Wojtan, Thürey, Gross & Turk's **complex edge test**
//! (`10.1145/1778765.1778787`) and is theirs, not ours — V-37. What is new is
//! running it as a correctness audit of *extraction*, across a family of
//! extractors on one fixture. See `validate::sealing` for the mechanics and for
//! the two degeneracies that had to be resolved before any number here means
//! anything.
//!
//! # Reading the columns
//!
//! `unsealed_walls` is a **hole**: the field separates two samples and the mesh
//! does not. `spurious_walls` is a **membrane**: the reverse. `mixed_regions`
//! is the same defect as the first at component scale, and is what survives when
//! several holes conspire. The two component counts are R-024's own wording —
//! *"compute connected components of the air sublevel set; compute connected
//! components of the mesh complement; assert they agree."*
//!
//! # Resolutions
//!
//! `[17, 25, 33]`, the suite's standard three — the same grids `golden.rs` pins
//! and the subgrid tests sweep, so a row here is comparable with every other
//! per-field number in the repo.

mod common;

use isomesh::MeshBuffer;
use isomesh::extractor::Extractor;
use isomesh::validate::sealing;

/// Samples per axis. `golden.rs`'s `RESOLUTIONS`, unchanged.
const RESOLUTIONS: [u32; 3] = [17, 25, 33];

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-21");
    common::experiment::run(prereg, |run| {
        println!(
            "{:<16} {:<28} {:>4} {:>7} {:>7} {:>8} {:>8} {:>5} {:>5} {:>5} {:>6}",
            "field",
            "extractor",
            "n",
            "probes",
            "walls",
            "unsealed",
            "spurious",
            "mixed",
            "f_air",
            "m_air",
            "degen"
        );

        // The two falsifiers, counted rather than eyeballed.
        let mut disagreements = 0usize;
        let mut primal_disagreements = 0usize;
        let mut rows = 0usize;

        isomesh::for_each_reference_field!(f64, |name, field| {
            for samples in RESOLUTIONS {
                let (shape, origin, cell_size) = common::grid::<f64, _>(&field, samples);
                // Inline block, so no `return` in here (M-253).
                isomesh::for_each_extractor!(f64, |ename, extractor| {
                    let mut out = MeshBuffer::<f64>::new();
                    let meshed = extractor
                        .extract_into(&field, &shape, origin, cell_size, &mut out)
                        .is_ok();
                    if meshed && !out.indices.is_empty() {
                        let r = sealing(
                            &field,
                            &shape,
                            origin,
                            cell_size,
                            &out.positions,
                            &out.indices,
                        );
                        rows += 1;
                        if !r.agrees() {
                            disagreements += 1;
                            // "Primal" is placement on the grid edge, which is
                            // what H's mechanism turns on -- everything whose
                            // name starts `marching`.
                            if ename.starts_with("marching") || ename.starts_with("subgrid") {
                                primal_disagreements += 1;
                            }
                        }
                        println!(
                            "{name:<16} {ename:<28} {samples:>4} {:>7} {:>7} {:>8} {:>8} \
                             {:>5} {:>5} {:>5} {:>6}",
                            r.probes,
                            r.field_walls,
                            r.unsealed_walls,
                            r.spurious_walls,
                            r.mixed_regions,
                            r.field_air_components,
                            r.mesh_air_components,
                            r.degenerate_probes
                        );
                        run.record(&[
                            ("field", name.to_string()),
                            ("extractor", ename.to_string()),
                            ("samples_per_axis", samples.to_string()),
                            ("field_air_components", r.field_air_components.to_string()),
                            ("mesh_air_components", r.mesh_air_components.to_string()),
                            ("unsealed_walls", r.unsealed_walls.to_string()),
                            (
                                "unsealed_on_domain_face",
                                r.unsealed_on_domain_face.to_string(),
                            ),
                            ("spurious_walls", r.spurious_walls.to_string()),
                            ("mixed_regions", r.mixed_regions.to_string()),
                            ("probes", r.probes.to_string()),
                            ("field_walls", r.field_walls.to_string()),
                            ("mesh_walls", r.mesh_walls.to_string()),
                            ("merged_crossings", r.merged_crossings.to_string()),
                            ("coplanar_probes", r.coplanar_probes.to_string()),
                            ("degenerate_triangles", r.degenerate_triangles.to_string()),
                            ("boundary_samples", r.boundary_samples.to_string()),
                            ("endpoint_crossings", r.endpoint_crossings.to_string()),
                            ("degenerate_probes", r.degenerate_probes.to_string()),
                            ("triangles", (out.indices.len() / 3).to_string()),
                            ("agrees", r.agrees().to_string()),
                        ]);
                    }
                });
            }
        });

        println!();
        println!("{rows} rows, {disagreements} disagreeing, {primal_disagreements} of them primal");
        if disagreements == 0 {
            println!(
                "P-21 FALSIFIED by its own stated falsifier: universal agreement. That is a \
                 stronger correctness statement than this crate has made before."
            );
        }
    });
}
