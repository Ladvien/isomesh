//! **P-66 - the monotone-edge witness, scored against the root finder.**
//!
//! Ticket: R-064. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p66
//! ```
//!
//! Writes `docs/experiments/p-66.csv`.
//!
//! # The two things being compared, and neither is an approximation of the other
//!
//! - **The oracle**: `subgrid::roots::all_roots`, which divides the edge into
//!   [`ORACLE_SAMPLES`] intervals and bisects every interval carrying a sign
//!   change. Its answer is a **root count**.
//! - **The witness**: `k` samples of `∇f · ê` along the edge. Its answer is a
//!   **boolean**, and the paper's claim is that `false` implies at most one root.
//!
//! So C1 is an implication in one direction and nothing else: multi-root implies
//! flagged. The converse is C2's false-positive rate, and it is expected to be
//! non-zero — a monotonicity violation with only one root is a real feature of
//! the field, not an error.
//!
//! # The oracle's resolution is a limit, and it is recorded rather than assumed
//!
//! `all_roots` cannot see two roots closer together than one of its intervals.
//! At [`ORACLE_SAMPLES`] per edge that is 1/128 of a cell. A root pair inside one
//! interval is invisible to **both** instruments, so it cannot produce a false
//! negative here — but it also means C1's "zero false negatives" is scoped to
//! pairs the oracle can resolve, and the entry has to say so. `oracle_samples`
//! is a column for exactly that reason.
//!
//! # Which edges
//!
//! Every axis-aligned **grid edge** of the sample lattice, which is the edge set
//! Marching Cubes actually interpolates along — `3 · (n-1) · n · n` of them. Not
//! tetrahedron edges: the paper is about the mesh's own edges, and this crate's
//! primary extractor is cell-based.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]

mod common;

use isomesh::Sdf;
use isomesh::fields::ReferenceField;
use isomesh::for_each_reference_field;
use isomesh::subgrid::roots::all_roots;

/// Samples per axis. The registered three.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// Directional-derivative samples per edge. The registered `k = 5` plus the
/// sweep C2's second half needs.
const KS: [u32; 5] = [2, 3, 5, 9, 17];

/// The `k` C1 and C3 are stated at.
const K_REGISTERED: u32 = 5;

/// Intervals the oracle divides each edge into.
///
/// 128 rather than `M-94`'s 1000: the sweep walks every grid edge at 65³, which
/// is 811,200 edges, and 1000 intervals each would be 811 million field
/// evaluations per field per resolution. 128 resolves a root pair separated by
/// more than 1/128 of a cell, which is finer than any feature Marching Cubes can
/// represent on that grid, and the number is a recorded column so the scope of
/// C1 is on the artefact rather than in this comment.
const ORACLE_SAMPLES: u32 = 128;

/// One `(field, resolution, k)` row.
struct Row {
    field: &'static str,
    samples: u32,
    k: u32,
    edges: u64,
    single_root: u64,
    multi_root: u64,
    flagged: u64,
    /// Multi-root edges the witness called monotonic. C1's whole subject.
    false_negatives: u64,
    /// Single-root edges the witness flagged.
    false_positives: u64,
}

/// Is this edge monotonic under the witness, at `k` samples?
///
/// `k` points **inclusive of both endpoints**, so `k = 2` is the endpoints alone
/// and `k = 5` adds three interior points. The paper's own count is
/// `max(2, ⌈‖e‖/w⌉ + 1)`; here `k` is swept directly, because C2's second half is
/// about how the rate moves with it.
///
/// A sampled projection of **exactly zero** is not a sign disagreement. Zero is
/// the boundary between the two signs and treating it as either would make the
/// test's answer depend on which endpoint it was compared against first; the
/// crate's own convention is that zero is not inside, and the analogue here is
/// that zero is not evidence.
fn is_monotonic<S: Sdf<Scalar = f64>>(sdf: &S, from: [f64; 3], to: [f64; 3], k: u32) -> bool {
    let dir = [to[0] - from[0], to[1] - from[1], to[2] - from[2]];
    let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
    if len == 0.0 {
        return true;
    }
    let unit = [dir[0] / len, dir[1] / len, dir[2] / len];

    let mut seen_positive = false;
    let mut seen_negative = false;
    for i in 0..k {
        let t = f64::from(i) / f64::from(k - 1);
        let p = [
            from[0] + dir[0] * t,
            from[1] + dir[1] * t,
            from[2] + dir[2] * t,
        ];
        let g = sdf.gradient(p);
        let proj = g[0] * unit[0] + g[1] * unit[1] + g[2] * unit[2];
        if proj > 0.0 {
            seen_positive = true;
        } else if proj < 0.0 {
            seen_negative = true;
        }
        if seen_positive && seen_negative {
            return false;
        }
    }
    true
}

/// Every grid edge of one field at one resolution, scored at every `k`.
///
/// The oracle runs **once per edge** and its answer is reused across the `k`
/// sweep. Re-running it per `k` would multiply the cost by five and could not
/// change an answer, and two oracle calls on one edge that disagreed would be a
/// determinism bug this harness is not the place to find.
fn measure<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    field_name: &'static str,
    samples: u32,
) -> Vec<Row> {
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);
    let at = |x: u32, y: u32, z: u32| -> [f64; 3] {
        [
            lo[0] + h * f64::from(x),
            lo[1] + h * f64::from(y),
            lo[2] + h * f64::from(z),
        ]
    };

    let mut rows: Vec<Row> = KS
        .iter()
        .map(|k| Row {
            field: field_name,
            samples,
            k: *k,
            edges: 0,
            single_root: 0,
            multi_root: 0,
            flagged: 0,
            false_negatives: 0,
            false_positives: 0,
        })
        .collect();

    let mut roots: Vec<f64> = Vec::with_capacity(8);
    let n = samples;
    for axis in 0..3usize {
        // An axis-aligned edge runs from sample `p` to `p + ê`, so that axis
        // stops one short and the other two run the full extent.
        let mut extent = [n, n, n];
        extent[axis] = n - 1;
        for z in 0..extent[2] {
            for y in 0..extent[1] {
                for x in 0..extent[0] {
                    let a = at(x, y, z);
                    let mut step = [x, y, z];
                    step[axis] += 1;
                    let b = at(step[0], step[1], step[2]);

                    roots.clear();
                    all_roots(a, b, field, ORACLE_SAMPLES, &mut roots);
                    let count = roots.len();

                    for row in &mut rows {
                        row.edges += 1;
                        let monotonic = is_monotonic(field, a, b, row.k);
                        if !monotonic {
                            row.flagged += 1;
                        }
                        if count > 1 {
                            row.multi_root += 1;
                            if monotonic {
                                row.false_negatives += 1;
                            }
                        } else if count == 1 {
                            row.single_root += 1;
                            if !monotonic {
                                row.false_positives += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    rows
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-66");
    let mut rows: Vec<Row> = Vec::new();

    for samples in RESOLUTIONS {
        for_each_reference_field!(f64, |name, field| {
            rows.extend(measure(&field, name, samples));
        });
    }

    println!(
        "{:>15} {:>5} {:>4} {:>9} {:>8} {:>8} {:>9} {:>6} {:>8} {:>9}",
        "field", "n", "k", "edges", "1-root", "n-root", "flagged", "FN", "FP", "FP rate"
    );
    for r in &rows {
        let fp_rate = if r.single_root == 0 {
            f64::NAN
        } else {
            r.false_positives as f64 / r.single_root as f64
        };
        println!(
            "{:>15} {:>5} {:>4} {:>9} {:>8} {:>8} {:>9} {:>6} {:>8} {:>9.4}",
            r.field,
            r.samples,
            r.k,
            r.edges,
            r.single_root,
            r.multi_root,
            r.flagged,
            r.false_negatives,
            r.false_positives,
            fp_rate
        );
    }

    // ── controls ─────────────────────────────────────────────────────────────
    let multi_total: u64 = rows
        .iter()
        .filter(|r| r.k == K_REGISTERED)
        .map(|r| r.multi_root)
        .sum();
    assert!(
        multi_total > 0,
        "VOID: the oracle found no multi-root edge anywhere in the sweep at {ORACLE_SAMPLES} \
         intervals, so C1's zero false negatives is a zero over an empty population"
    );
    let flagged_total: u64 = rows
        .iter()
        .filter(|r| r.k == K_REGISTERED)
        .map(|r| r.flagged)
        .sum();
    assert!(
        flagged_total > 0,
        "VOID: the witness flagged nothing at k = {K_REGISTERED}, so it cannot report the bad news"
    );
    let unflagged: u64 = rows
        .iter()
        .filter(|r| r.k == K_REGISTERED)
        .map(|r| r.edges - r.flagged)
        .sum();
    assert!(
        unflagged > 0,
        "VOID: the witness flagged EVERY edge at k = {K_REGISTERED}, so C1 passes trivially"
    );

    // ── verdict ──────────────────────────────────────────────────────────────
    let fields: Vec<&'static str> = {
        let mut v: Vec<&'static str> = rows.iter().map(|r| r.field).collect();
        v.sort_unstable();
        v.dedup();
        v
    };

    // C1: zero false negatives at the registered k, over the whole sweep.
    let false_negatives: u64 = rows
        .iter()
        .filter(|r| r.k == K_REGISTERED)
        .map(|r| r.false_negatives)
        .sum();
    let c1 = false_negatives == 0;

    // C2, first half: rate below 20% at k = 5 on the six smooth fields. The two
    // that are not smooth are named here rather than guessed: `fbm_terrain` and
    // `noise_cavity` are the noise fields, and the registration's "six smooth"
    // is the other six.
    let rough = ["fbm_terrain", "noise_cavity"];
    let mut c2_rate = true;
    for f in &fields {
        if rough.contains(f) {
            continue;
        }
        for r in rows
            .iter()
            .filter(|r| r.field == *f && r.k == K_REGISTERED && r.single_root > 0)
        {
            let rate = r.false_positives as f64 / r.single_root as f64;
            if rate >= 0.20 {
                c2_rate = false;
                println!("C2 rate {f} at {}³: {rate:.4} >= 0.20", r.samples);
            }
        }
    }

    // C2, second half: the rate falls as k rises. Measured at the finest
    // resolution, where there is the most data.
    let mut falls_with_k: Vec<(&'static str, bool, Vec<f64>)> = Vec::new();
    for f in &fields {
        let seq: Vec<f64> = KS
            .iter()
            .filter_map(|k| {
                rows.iter()
                    .find(|r| r.field == *f && r.samples == 65 && r.k == *k)
                    .filter(|r| r.single_root > 0)
                    .map(|r| r.false_positives as f64 / r.single_root as f64)
            })
            .collect();
        let falls = seq.windows(2).all(|w| w[1] <= w[0]);
        println!(
            "C2 falls with k, {f} at 65³: {} -> {}",
            seq.iter()
                .map(|v| format!("{v:.4}"))
                .collect::<Vec<_>>()
                .join(" "),
            if falls { "falling" } else { "NOT FALLING" }
        );
        falls_with_k.push((f, falls, seq));
    }
    let c2 = c2_rate && falls_with_k.iter().all(|(_, ok, _)| *ok);

    // C3: the non-monotonic fraction falls with resolution, and `thin_plate`
    // ranks first at the registered k.
    let mut falls_with_res: Vec<(&'static str, bool, Vec<f64>)> = Vec::new();
    for f in &fields {
        let seq: Vec<f64> = RESOLUTIONS
            .iter()
            .filter_map(|n| {
                rows.iter()
                    .find(|r| r.field == *f && r.samples == *n && r.k == K_REGISTERED)
                    .map(|r| r.flagged as f64 / r.edges as f64)
            })
            .collect();
        let falls = seq.windows(2).all(|w| w[1] <= w[0]);
        println!(
            "C3 falls with resolution, {f}: {} -> {}",
            seq.iter()
                .map(|v| format!("{v:.4}"))
                .collect::<Vec<_>>()
                .join(" "),
            if falls { "falling" } else { "NOT FALLING" }
        );
        falls_with_res.push((f, falls, seq));
    }
    let top = falls_with_res
        .iter()
        .max_by(|a, b| {
            a.2.first()
                .unwrap_or(&0.0)
                .partial_cmp(b.2.first().unwrap_or(&0.0))
                .expect("finite")
        })
        .map(|(f, _, _)| *f)
        .unwrap_or("none");
    let thin_first = top == "thin_plate";
    println!("C3 highest non-monotonic fraction at 17³: {top}");
    let c3 = falls_with_res.iter().all(|(_, ok, _)| *ok) && thin_first;

    println!(
        "\nC1 zero false negatives at k = {K_REGISTERED} over {multi_total} multi-root edges: \
         {false_negatives} -> {}",
        if c1 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C2 rate < 20% on the six smooth AND falling with k -> {}",
        if c2 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C3 falling with resolution on all eight AND thin_plate first -> {}",
        if c3 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "\nSCOPE: the oracle divides each edge into {ORACLE_SAMPLES} intervals and cannot resolve \
         a root pair closer together than one of them, so C1's zero is scoped to pairs it can \
         see. That is a recorded column, not a footnote."
    );

    common::experiment::run(prereg, |run| {
        for r in &rows {
            let fp_rate = if r.single_root == 0 {
                f64::NAN
            } else {
                r.false_positives as f64 / r.single_root as f64
            };
            let fk = falls_with_k
                .iter()
                .find(|(f, _, _)| *f == r.field)
                .is_some_and(|(_, ok, _)| *ok);
            let fr = falls_with_res
                .iter()
                .find(|(f, _, _)| *f == r.field)
                .is_some_and(|(_, ok, _)| *ok);
            run.record(&[
                ("field", r.field.to_string()),
                ("samples_per_axis", r.samples.to_string()),
                ("k", r.k.to_string()),
                ("oracle_samples", ORACLE_SAMPLES.to_string()),
                ("edges", r.edges.to_string()),
                ("single_root_edges", r.single_root.to_string()),
                ("multi_root_edges", r.multi_root.to_string()),
                ("flagged_non_monotonic", r.flagged.to_string()),
                ("false_negatives", r.false_negatives.to_string()),
                ("false_positives", r.false_positives.to_string()),
                ("false_positive_rate", format!("{fp_rate:.6}")),
                (
                    "non_monotonic_fraction",
                    format!("{:.6}", r.flagged as f64 / r.edges as f64),
                ),
                ("falls_with_k", fk.to_string()),
                ("falls_with_resolution", fr.to_string()),
                ("thin_plate_ranks_first", thin_first.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
            ]);
        }
    });
}
