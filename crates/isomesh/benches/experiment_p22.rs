//! **P-22 — mean-ratio triangle quality, against a baseline borrowed from another
//! implementation.**
//!
//! Ticket: T-026. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p22
//! ```
//!
//! Writes `docs/experiments/p-22.csv`.
//!
//! # The metric, and why it is this one
//!
//! Grosso & Zint's mean ratio (`10.1007/s00371-021-02139-w` §5),
//! `q = 4√3·A / Σᵢlᵢ²` — **1 for equilateral, 0 for degenerate**. Chosen over
//! the `AR > 4` figures the neural-isosurfacing line reports because those come
//! from a differentiable-rendering loop on **learned** fields and are not
//! comparable with meshing an analytic field on a uniform grid (V-38). Grosso &
//! Zint mesh uniform grids and report Marching Cubes, topologically correct
//! Marching Cubes and Dual Contouring **by name** (V-39).
//!
//! # Two clauses, and the second is the risk
//!
//! **Clause 1 is the paper's own observation.** Its MC and TMC columns agree to
//! two decimals on all seven rows, because both *"place their vertices also on
//! the trilinear interpolant but along the voxel edges not within the voxel
//! cells."* This crate has that pair — `marching_cubes` and
//! `marching_cubes+decider` — so a face rule should change which crossings are
//! joined and not where they are.
//!
//! **Clause 2 borrows a band from somebody else's code.** Their MC sits at
//! 0.65–0.71, and `gen2` — their synthetic volume — is flat at 0.71 across 64³,
//! 128³ and 256³. Every cross-source comparison this project has attempted has
//! needed an amendment, so this is registered as a real risk rather than a
//! formality.
//!
//! # Irregular vertices carry a caveat, stated before the numbers
//!
//! Valence ≠ 6 is their definition and is written for **closed** medical
//! volumes. Every boundary vertex on an open field is irregular by construction,
//! so `gyroid` and `fbm_terrain` read high for a reason that has nothing to do
//! with triangle quality. The column is raw to match the published definition.

mod common;

use isomesh::MeshBuffer;
use isomesh::extractor::Extractor;
use isomesh::validate::{ValidateConfig, validate_indexed};

/// Samples per axis. The suite's standard three.
const RESOLUTIONS: [u32; 3] = [17, 25, 33];

/// The band Grosso & Zint's Marching Cubes occupies, Table 7.
const THEIR_MC: (f64, f64) = (0.65, 0.71);

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-22");
    common::experiment::run(prereg, |run| {
        println!(
            "{:<16} {:<28} {:>4} {:>10} {:>9} {:>7} {:>8}",
            "field", "extractor", "n", "mean_ratio", "irreg", "verts", "tris"
        );

        // Clause 1: the two Marching Cubes entries, per (field, resolution).
        let mut mc_pairs = 0usize;
        let mut mc_pairs_differing = 0usize;
        // Clause 2: plain Marching Cubes on the smooth analytic fields.
        let mut smooth_rows = 0usize;
        let mut smooth_outside = 0usize;

        isomesh::for_each_reference_field!(f64, |name, field| {
            // The smooth analytic fields, which is what their `gen2` row is.
            let smooth = matches!(name, "sphere" | "torus");
            for samples in RESOLUTIONS {
                let (shape, origin, cell_size) = common::grid::<f64, _>(&field, samples);
                let cfg = ValidateConfig::from_cell_size(cell_size).expect("valid cell size");
                let mut plain = f64::NAN;
                let mut decided = f64::NAN;

                isomesh::for_each_extractor!(f64, |ename, extractor| {
                    let mut out = MeshBuffer::<f64>::new();
                    let ok = extractor
                        .extract_into(&field, &shape, origin, cell_size, &mut out)
                        .is_ok();
                    if ok && !out.indices.is_empty() {
                        let r = validate_indexed(&out.positions, &out.indices, &cfg);
                        if ename == "marching_cubes" {
                            plain = r.mean_ratio;
                            if smooth {
                                smooth_rows += 1;
                                if r.mean_ratio < THEIR_MC.0 || r.mean_ratio > THEIR_MC.1 {
                                    smooth_outside += 1;
                                }
                            }
                        }
                        if ename == "marching_cubes+decider" {
                            decided = r.mean_ratio;
                        }
                        println!(
                            "{name:<16} {ename:<28} {samples:>4} {:>10.4} {:>9} {:>7} {:>8}",
                            r.mean_ratio, r.irregular_vertices, r.referenced_vertices, r.faces
                        );
                        run.record(&[
                            ("field", name.to_string()),
                            ("extractor", ename.to_string()),
                            ("samples_per_axis", samples.to_string()),
                            ("mean_ratio", format!("{:.9}", r.mean_ratio)),
                            ("irregular_vertices", r.irregular_vertices.to_string()),
                            ("referenced_vertices", r.referenced_vertices.to_string()),
                            ("triangles", r.faces.to_string()),
                            ("degenerate_triangles", r.degenerate_triangles.to_string()),
                            ("boundary_edges", r.boundary_edges.to_string()),
                        ]);
                    }
                });

                // Bit-identical, not merely close: the claim is that the face
                // rule moves no geometry at all.
                if plain.is_finite() && decided.is_finite() {
                    mc_pairs += 1;
                    #[allow(clippy::float_cmp, reason = "identity is the hypothesis")]
                    let same = plain == decided;
                    if !same {
                        mc_pairs_differing += 1;
                        println!(
                            "  DIFFER {name} {samples}: plain {plain:.12} decider {decided:.12}"
                        );
                    }
                }
            }
        });

        println!();
        println!(
            "clause 1: {mc_pairs_differing} of {mc_pairs} Marching Cubes pairs differ \
             (H says 0)"
        );
        println!(
            "clause 2: {smooth_outside} of {smooth_rows} smooth-field rows outside \
             {:.2}-{:.2} (H says 0)",
            THEIR_MC.0, THEIR_MC.1
        );
    });
}
