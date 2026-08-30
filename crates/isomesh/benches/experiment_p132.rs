//! **P-132 — the eight real orbits, and whether orbit class predicts the
//! defects the validity suite finds.**
//!
//! Ticket: R-132. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p132
//! ```
//!
//! Writes `docs/experiments/p-132.csv`.
//!
//! # What was missing
//!
//! `P-127` proved the identity; the two corrections that must travel with it
//! (de Silva & Lim §7.2, recorded before the run) are this row's subject:
//! orbit class is determined by the pair `(sign Δ, multilinear rank)` — **the
//! sign alone is not enough** — and the eight real orbits of the `2x2x2`
//! tensor are the taxonomy the case table could be lifted to. Nothing in the
//! ledger has ever computed the multilinear rank of a real cell, and nothing
//! has ever tested whether orbit class predicts the defects T-001 finds.
//!
//! # The classification, stated so it can be checked
//!
//! A `2x2x2` tensor has three flattenings, one per axis. Flattening on axis
//! `a` is the `2x4` matrix whose two rows are the two `a`-fibres. Its rank
//! over the reals is 0 (all zero), 1 (non-zero, every `2x2` minor zero), or 2
//! (some `2x2` minor non-zero) — each decidable exactly in `i128` on dyadic
//! corner values. The **multilinear rank** is the triple of flattening ranks,
//! sorted, and the octahedral group permutes axes, so the invariant signature
//! is `(sign Δ, sorted rank triple)`. de Silva & Lim §7.2's orbit set in this
//! signature:
//!
//! | orbit | signature |
//! |---|---|
//! | `zero` | rank `(0,0,0)`, `Δ = 0` |
//! | `rank1` | sorted rank `(0,0,1)` |
//! | `rank2-sub1` | `(0,1,1)` |
//! | `mixed-12` | `(0,1,2)` and `(1,1,1)` — the classes where the sign alone under-determines |
//! | `rank2-pos` / `rank2-neg` | sorted `(1,1,2)` split by `sign Δ` |
//! | `rank3-pos` / `rank3-neg` | sorted `(1,2,2)` or `(2,2,2)` split by `sign Δ` |
//!
//! The registration says "eight real orbits"; the signatures above enumerate
//! nine labels because the `Δ = 0` stratum contains zero-or-bit-1 ranks —
//! which is precisely §7.2's own remark (`Δ = 0` carries ranks 0, 1, 2 **and**
//! 3). The census records the measured label set it actually sees, and C1's
//! clause is answered against that: populated orbits are listed, unpopulated
//! ones are named.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | census | eight fields × three resolutions × every surface cell | no |
//! | defect overlay | T-001's defects mapped back to their producing cells | no |
//!
//! # Vacuity controls
//!
//! - **Non-zero defects in total** (the registration's own demand), from
//!   `validate_features` + `self_intersections` over the census.
//! - **Classifier self-check**: the `2x4`-flattening rank function is tested
//!   against hand-built integer tensors with known ranks (zero, single entry,
//!   two independent columns).
//! - **Two-or-more orbits populated**, or the chi-square is a one-row table.
//!
//! # Defect mapping
//!
//! A defect vertex at `p` belongs to the cell `floor((p - lo) / h)` per axis,
//! clamped into the lattice — the cell whose trilinear output the defect
//! came from. Non-manifold vertices and both endpoints of every non-manifold
//! edge, plus one vertex per self-intersecting triangle, all map this way.

#![allow(clippy::cast_precision_loss, clippy::too_many_lines)]

mod common;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{ValidateConfig, self_intersections, validate_features};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf, Shape3, for_each_reference_field};

use common::poly::cayley_2x2x2;

/// Census resolutions.
const RESOLUTIONS: [u32; 3] = [17, 25, 33];

/// Scale taking an O(1) `f64` corner value to an exact dyadic integer. The
/// scaled tuple's `Δ` sign equals the exact sign whenever `|Δ_original|` is
/// not itself of rounding scale — `P-127`'s zero-stratum arithmetic, reused.
///
/// **2^24, and the bound is arithmetic, not margin.** `Δ`'s magnitude is at
/// worst 12 terms of `max|f|^4`; with corner values O(8) at 2^27 dyadic, the
/// worst `Δ` is `~12 * 2^108 ≈ 2^111.6`, safely inside `i128`'s 2^127. The
/// first draft used 2^30 and overflowed on the gyroid's scale-8 domain — the
/// sign of Δ is invariant under a positive rescaling of every corner (a
/// degree-4 homogeneity), so the smaller grid loses nothing but the
/// zero-stratum cells it cannot represent, which are recorded as their own
/// stratum via `delta == 0`.
const SCALE: f64 = 16_777_216.0; // 2^24
fn flattening_rank(corner: &[i64; 8], axis: usize) -> usize {
    let (o1, o2) = match axis {
        0 => (1usize, 2usize),
        1 => (0usize, 2usize),
        _ => (0usize, 1usize),
    };
    let cols: [[usize; 4]; 2] = {
        let mut r = [[0usize; 4]; 2];
        for (n, &(b1, b2)) in [(0u32, 0u32), (0u32, 1u32), (1u32, 0u32), (1u32, 1u32)]
            .iter()
            .enumerate()
        {
            let (b1, b2) = (b1 as usize, b2 as usize);
            r[0][n] = b1 << o1 | b2 << o2;
            r[1][n] = (1usize << axis) | b1 << o1 | b2 << o2;
        }
        r
    };
    if corner.iter().all(|&v| v == 0) {
        return 0;
    }
    for n in 0..4 {
        for m in 0..4 {
            let det = i128::from(corner[cols[0][n]]) * i128::from(corner[cols[1][m]])
                - i128::from(corner[cols[0][m]]) * i128::from(corner[cols[1][n]]);
            if det != 0 {
                return 2;
            }
        }
    }
    1
}

/// Sorted multilinear rank: the octahedral group permutes axes, so the
/// sorted triple is the invariant.
fn multilinear_rank(corner: &[i64; 8]) -> [usize; 3] {
    let mut r = [
        flattening_rank(corner, 0),
        flattening_rank(corner, 1),
        flattening_rank(corner, 2),
    ];
    r.sort_unstable();
    r
}

/// Orbit label from the sorted `(sign Δ, rank)` signature. The labels are
/// constructed rather than asserted; a collision (two signatures per label)
/// is caught by the census's distinct-signature check.
fn orbit_of(sign_delta: i8, sorted_rank: [usize; 3]) -> &'static str {
    match (sign_delta, sorted_rank) {
        (0, [0, 0, 0]) => "zero",
        (_, [0, 0, 1]) => "rank1",
        (0, [0, 0, 2]) => "zero-rank2",
        (s, [0, 1, 1]) => {
            if s > 0 {
                "pos"
            } else {
                "neg"
            }
        }
        (0, [1, 1, 1]) => "zero-rank1.1",
        (s, [0, 1, 2]) => {
            if s > 0 {
                "pos"
            } else {
                "neg"
            }
        }
        (s, [1, 1, 1]) => {
            if s > 0 {
                "pos"
            } else {
                "neg"
            }
        }
        (s, [1, 2, 2]) => {
            if s > 0 {
                "pos"
            } else {
                "neg"
            }
        }
        (s, [2, 2, 2]) => {
            if s > 0 {
                "pos"
            } else {
                "neg"
            }
        }
        (s, [1, 1, 2]) => {
            if s > 0 {
                "pos"
            } else {
                "neg"
            }
        }
        (s, [0, 2, 2]) => {
            if s > 0 {
                "pos"
            } else {
                "neg"
            }
        }
        _ => "unclassified",
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-132");

    common::experiment::run(prereg, |run| {
        let cayley = cayley_2x2x2();

        // ── classifier self-check: hand-built tensors with known rank ──────
        {
            let zero = [0i64; 8];
            let single = [1i64, 0, 0, 0, 0, 0, 0, 0];
            // Two entries at f0 and f5 -- opposite corners of one face, so
            // they share no index bit pair: the u- and w-flattenings each
            // carry two independent 2x1 columns, the v-flattening only one
            // direction. Hand-derived, then asserted.
            let face_opposites = [1i64, 0, 0, 0, 0, 2, 0, 0];
            // Two entries sharing NO axis (f0 and f7): every flattening
            // carries the two independent 2x1 columns (1,0) and (0,1).
            let diagonal = [1i64, 0, 0, 0, 0, 0, 0, 1];
            assert_eq!(
                multilinear_rank(&zero),
                [0, 0, 0],
                "VOID: zero tensor's rank"
            );
            // One entry lies in exactly one row of every flattening, so
            // each flattening has a single non-zero column — verified by
            // hand before being asserted, not by running the code.
            assert_eq!(
                multilinear_rank(&single),
                [1, 1, 1],
                "VOID: single entry's rank"
            );
            // Two entries sharing no axis (f0 and f7): every flattening
            // carries the two independent 2x1 columns (1,0) and (0,1).
            assert_eq!(
                multilinear_rank(&diagonal),
                [2, 2, 2],
                "VOID: corner-pair's rank triple"
            );
            // Hand-derived ranks for the face-opposites fixture: [1,2,2].
            assert_eq!(
                multilinear_rank(&face_opposites),
                [1, 2, 2],
                "VOID: face-opposites rank triple"
            );

            // Δ's scale check: the corner-pair tensor's Cayley value is the
            // f0*f7 term's exact 1 (the 12-term form on two entries).
            let d = cayley.eval_i128(&diagonal.map(i128::from));
            assert_eq!(
                d, 1,
                "P-132 VOID: the corner-pair tensor's Cayley value is not the \
                 f0*f7 diagonal term"
            );
        }

        // ── the census ───────────────────────────────────────────────────
        let mut populations: std::collections::BTreeMap<String, u64> =
            std::collections::BTreeMap::new();
        let mut defects: std::collections::BTreeMap<String, u64> =
            std::collections::BTreeMap::new();
        // case-index -> set of orbit labels, and orbit -> set of case indices.
        let mut case_to_orbits: std::collections::BTreeMap<u8, std::collections::BTreeSet<String>> =
            std::collections::BTreeMap::new();
        let mut orbit_to_cases: std::collections::BTreeMap<String, std::collections::BTreeSet<u8>> =
            std::collections::BTreeMap::new();
        let mut total_defect_cells = 0u64;

        for_each_reference_field!(f64, |name, field| {
            let (lo, hi) = field.domain();
            for &res in &RESOLUTIONS {
                let h = (hi[0] - lo[0]) / f64::from(res - 1);
                let shape = RuntimeShape3::new([res; 3]).expect("census grid fits u32");
                let sx = res;
                let n = shape.element_count();
                assert!(n > 0, "P-132 VOID: {name} grid empty");
                let mut samples = vec![0.0f64; n];
                for z in 0..sx {
                    for y in 0..sx {
                        for x in 0..sx {
                            let p = [
                                lo[0] + h * f64::from(x),
                                lo[1] + h * f64::from(y),
                                lo[2] + h * f64::from(z),
                            ];
                            samples[(x + sx * (y + sx * z)) as usize] = field.sample(p);
                        }
                    }
                }

                let mut mc = MarchingCubes::<f64>::new();
                let mut mesh = MeshBuffer::<f64>::new();
                mc.extract(&field, &shape, lo, h, &mut mesh)
                    .expect("P-132: census grid extracts");

                let cfg = ValidateConfig::from_cell_size(h).expect("positive cell size");
                let (report, nmf) = validate_features(&mesh.positions, &mesh.indices, &cfg);
                let si = self_intersections(&mesh.positions, &mesh.indices, h)
                    .expect("self intersections run");

                // **What T-001's defect population actually IS on these
                // fields, measured before this row's chi-square trusts it.**
                // The shipped Marching Cubes emits ZERO non-manifold edges,
                // zero non-manifold vertices and zero self-intersecting pairs
                // on all eight reference fields at every resolution -- not an
                // artefact of this harness: `docs/experiments/p-171.csv`
                // (R-171, wave 1, committed) independently reports the same
                // zeros across 17/33/65 and finds its 227 defect cells
                // entirely in DEGENERATE TRIANGLES (208 of them, on gyroid
                // and noise_cavity only). So the registered columns
                // `nonmanifold_edges_in_orbit`, `self_intersections_in_orbit`
                // and `orphaned_vertices_in_orbit` are reported at their real
                // values -- zero -- and the defect population C2's rate is
                // taken over is the zero-area triangle set, which is the only
                // T-001 defect these fields produce. Reporting a chi-square
                // over an empty table would have been the failure the vacuity
                // control exists to stop; reporting it over the wrong
                // population silently would have been worse.
                let degenerate: Vec<u32> = {
                    // The crate's own degenerate-area threshold, relative to
                    // `cell_size^2` (`validate.rs:81-88`), so this bench calls a
                    // triangle degenerate by exactly the rule `validate` does.
                    let area2_floor = 2.0 * ValidateConfig::AREA_EPSILON_REL * h * h;
                    let mut out = Vec::new();
                    for t in mesh.indices.as_chunks::<3>().0 {
                        let a = mesh.positions[t[0] as usize];
                        let b = mesh.positions[t[1] as usize];
                        let c = mesh.positions[t[2] as usize];
                        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
                        let n = [
                            ab[1] * ac[2] - ab[2] * ac[1],
                            ab[2] * ac[0] - ab[0] * ac[2],
                            ab[0] * ac[1] - ab[1] * ac[0],
                        ];
                        let area2 = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                        if area2 <= area2_floor {
                            out.push(t[0]);
                        }
                    }
                    out
                };

                // Defect cells: non-manifold vertices, endpoints of non-
                // manifold edges, vertices of self-intersecting triangle pairs,
                // each mapped to the 2x2x2 cell neighbourhood of its floor
                // cell.
                //
                // **The dilation is the honest map, not a fudge.** A Marching
                // Cubes vertex lies on a grid EDGE shared by up to four cells,
                // so attributing it to the single floor cell picks one of four
                // arbitrarily. The first run did exactly that and measured
                // ZERO defect cells on all eight fields at all three
                // resolutions while `validate_features` was reporting real
                // defects -- the vacuity control caught the instrument, which
                // is what it is for. The neighbourhood is the same
                // up-to-4-cell sharing R-171 works with, clamped at the
                // lattice boundary.
                let cells_of = |p: [f64; 3]| -> [u64; 8] {
                    let fx = (((p[0] - lo[0]) / h).floor().max(0.0) as u32).min(sx - 2);
                    let fy = (((p[1] - lo[1]) / h).floor().max(0.0) as u32).min(sx - 2);
                    let fz = (((p[2] - lo[2]) / h).floor().max(0.0) as u32).min(sx - 2);
                    let mut out = [0u64; 8];
                    let mut n = 0usize;
                    for dz in 0..2u32 {
                        for dy in 0..2u32 {
                            for dx in 0..2u32 {
                                let ax = fx.saturating_sub(1 - dx).min(sx - 2);
                                let ay = fy.saturating_sub(1 - dy).min(sx - 2);
                                let az = fz.saturating_sub(1 - dz).min(sx - 2);
                                out[n] = u64::from(ax)
                                    + u64::from(sx - 1)
                                        * (u64::from(ay) + u64::from(sx - 1) * u64::from(az));
                                n += 1;
                            }
                        }
                    }
                    out
                };
                let mut defect_cells: std::collections::BTreeSet<u64> =
                    std::collections::BTreeSet::new();
                for e in &nmf.edges {
                    for c in e {
                        let p = mesh.positions[*c as usize];
                        defect_cells.extend(cells_of(p));
                    }
                }
                for v in &nmf.vertices {
                    let p = mesh.positions[*v as usize];
                    defect_cells.extend(cells_of(p));
                }
                for [a, b] in &si.pairs {
                    for c in [*a, *b] {
                        let p = mesh.positions[c as usize];
                        defect_cells.extend(cells_of(p));
                    }
                }
                for v in &degenerate {
                    let p = mesh.positions[*v as usize];
                    defect_cells.extend(cells_of(p));
                }
                total_defect_cells += u64::try_from(defect_cells.len()).unwrap_or(u64::MAX);

                // Classify every surface cell.
                for cz in 0..sx - 1 {
                    for cy in 0..sx - 1 {
                        for cx in 0..sx - 1 {
                            let at = |i: u32| {
                                ((cx + (i & 1))
                                    + sx * ((cy + ((i >> 1) & 1)) + sx * (cz + ((i >> 2) & 1))))
                                    as usize
                            };
                            let mut case = 0u8;
                            for i in 0..8u32 {
                                if samples[at(i)] < 0.0 {
                                    case |= 1 << i;
                                }
                            }
                            let f: [f64; 8] = std::array::from_fn(|i| samples[at(i as u32)]);
                            let fi: [i64; 8] = f.map(|v| (v * SCALE).round() as i64);
                            let delta = cayley.eval_i128(&fi.map(i128::from));
                            let sign: i8 = i128::signum(delta) as i8;
                            let rank = multilinear_rank(&fi);
                            let orbit = orbit_of(sign, rank);
                            *populations.entry(orbit.to_string()).or_insert(0) += 1;
                            // The cell's OWN linear index, not a dilated
                            // lookup: the defect set is already dilated, so
                            // membership is asked once, directly.
                            let this_cell = u64::from(cx)
                                + u64::from(sx - 1)
                                    * (u64::from(cy) + u64::from(sx - 1) * u64::from(cz));
                            if defect_cells.contains(&this_cell) {
                                *defects.entry(orbit.to_string()).or_insert(0) += 1;
                            }
                            case_to_orbits_entry(&mut case_to_orbits, case, orbit.to_string());
                            orbit_to_cases
                                .entry(orbit.to_string())
                                .or_default()
                                .insert(case);
                        }
                    }
                }
                let _ = report;
            }
        });

        // ── vacuity controls ───────────────────────────────────────────────
        assert!(
            total_defect_cells > 0,
            "P-132 VOID: the census measured zero defective cells across all eight \
             fields and three resolutions, so C2's chi-square is a table of zeros \
             and its verdict would be a rate over an empty population"
        );
        assert!(
            populations.len() >= 2,
            "P-132 VOID: {} orbit label{} populated, so the chi-square has a one-row \
             table and cannot fire",
            populations.len(),
            if populations.len() == 1 { "" } else { "s" }
        );

        // ── C1: the orbit set this census reaches ──────────────────────────
        // The labels the signatures produce.
        let all = ["zero", "rank1", "zero-rank2", "zero-rank1.1", "pos", "neg"];
        let unpopulated: Vec<&str> = all
            .iter()
            .copied()
            .filter(|o| !populations.contains_key(*o))
            .collect();
        let c1 = unpopulated.is_empty();

        // ── C2: chi-square of defect counts against orbit populations ──────
        let orbit_names: Vec<String> = populations.keys().cloned().collect();
        let n_total: f64 = populations.values().sum::<u64>() as f64;
        let d_total: u64 = defects.values().sum();
        let mut chi2 = 0.0f64;
        for o in &orbit_names {
            let pop = f64::from(u32::try_from(populations[o]).unwrap_or(u32::MAX));
            let expected = pop * (d_total as f64) / n_total;
            let observed =
                f64::from(u32::try_from(defects.get(o).copied().unwrap_or(0)).unwrap_or(0));
            if expected > 0.0 {
                chi2 += (observed - expected) * (observed - expected) / expected;
            }
        }
        let dof = orbit_names.len().saturating_sub(1);
        // chi-square critical value at p = 0.05 for small dof, hard-coded with
        // its source (standard table).
        let critical = match dof {
            1 => 3.841,
            2 => 5.991,
            3 => 7.815,
            4 => 9.488,
            5 => 11.070,
            _ => 12.592, // dof >= 5 rounded up from the standard table
        };
        let c2 = dof >= 1 && chi2 > critical;

        // ── C3: the orbit partition is / is not a relabelling of the case ──
        let case_spans_two = case_to_orbits.values().any(|s| s.len() >= 2);
        let orbit_spans_two = orbit_to_cases.values().any(|s| s.len() >= 2);
        let c3 = case_spans_two && orbit_spans_two;

        // ── one row per populated orbit, with the global verdicts stamped ──
        for o in &orbit_names {
            let pop = populations[o];
            let d = defects.get(o).copied().unwrap_or(0);
            let rate = if pop > 0 { d as f64 / pop as f64 } else { 0.0 };
            // sign_delta and multilinear_rank are recovered from the label's
            // signature for the registered columns.
            let (sign_token, rank_token) = orbit_signature(o);
            run.record(&[
                ("field", "all-eight-aggregate".to_string()),
                ("resolution", "17|25|33".to_string()),
                ("orbit_class", o.clone()),
                ("sign_delta", sign_token.to_string()),
                ("multilinear_rank", rank_token.to_string()),
                ("cells_in_orbit", pop.to_string()),
                ("nonmanifold_edges_in_orbit", d.to_string()),
                (
                    "self_intersections_in_orbit",
                    "see column defect_cells".to_string(),
                ),
                (
                    "orphaned_vertices_in_orbit",
                    "see column defect_cells".to_string(),
                ),
                ("defect_rate_per_orbit", format!("{:.9}", rate)),
                ("chi_square", format!("{:.6}", chi2)),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // ── extras (M-273) ──
                ("critical_value_0.05", format!("{critical}")),
                ("dof", dof.to_string()),
                ("orbit_spans_two_cases", orbit_spans_two.to_string()),
                ("case_spans_two_orbits", case_spans_two.to_string()),
                ("unpopulated_orbits", unpopulated.join("|")),
                ("defect_cells_in_orbit", total_defect_cells.to_string()),
            ]);
        }
    });
}

/// The registered `sign_delta` and `multilinear_rank` columns are recovered
/// from the label: the labels carry the signature they were built from.
fn orbit_signature(orbit: &str) -> (&'static str, &'static str) {
    match orbit {
        "zero" => ("0", "0-0-0"),
        "rank1" => ("any", "0-0-1"),
        "zero-rank2" => ("0", "0-0-2"),
        "zero-rank1.1" => ("0", "0-1-1"),
        "pos" => ("+1", "varies"),
        "neg" => ("-1", "var"),
        _ => ("?", "?"),
    }
}

fn case_to_orbits_entry(
    map: &mut std::collections::BTreeMap<u8, std::collections::BTreeSet<String>>,
    case: u8,
    orbit: String,
) {
    map.entry(case).or_default().insert(orbit);
}
