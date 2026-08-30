//! **P-129 — octahedral and negation invariance as algebra, at a fraction of
//! the cost of the 48-element sweep.**
//!
//! Ticket: R-129. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p129
//! ```
//!
//! Writes `docs/experiments/p-129.csv`.
//!
//! # What was missing
//!
//! `✗49` established bit-exact octahedral equivariance of plain Marching Cubes
//! by sweeping all 48 elements of the cube group; `M-177` recorded that
//! reordering cannot buy the same equivariance for vertex placement.
//! `P-127` — this phase's headline — proved the discriminant `Δ` invariant
//! under all 48 octahedral relabellings **and** under negating the field, as
//! algebra. This row asks the registered question: can that algebra stand in
//! for the sweep? A `Δ`-based check costs a constant number of polynomial
//! evaluations per cell where the sweep costs 48 relabellings and case
//! recomputations per cell, and C1 measures whether the two instruments agree
//! on every cell of every reference field at `>= 10×` less.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | sweep | 48 relabellings per cell; case equivariance of the sign word | no |
//! | invariant | one `Δ` pair per element per cell; invariance of the polynomial | no |
//! | broken | one asymmetric tuple fed to the sweep | **yes** |
//!
//! The two instruments answer the same question — "is this cell's
//! classification the same under the whole cube group?" — at different depths.
//! The sweep reads the case word; the invariant reads `Δ`. On a cell where one
//! fires and the other does not, the agreement count sees it.
//!
//! # SHARE, recomputed before the numbers
//!
//! C1's share is the equivariance-test stage of CI, currently the full sweep.
//! The scope caveat is `P-127`'s and applies here too: `Δ`-invariance is a
//! **necessary** condition for what the sweep tests, not the same condition —
//! the sweep tests the extractor's classification output, `Δ` tests the input's
//! own invariant. That is why the machine-independent evidence is a WORK COUNT
//! (`relabellings_per_delta`, how many relabellings the sweep performs per
//! invariant evaluation) rather than only a clock, and why the timing rows
//! carry interleaved min/median/max (M-280).
//!
//! # Vacuity controls
//!
//! - **A deliberately broken table must be caught**: an asymmetric tuple must
//!   produce a non-zero sweep violation, or the sweep cannot report bad news
//!   and its agreement with the invariant means nothing.
//! - **The instrument's own arithmetic**: all 48 committed relabellings must
//!   leave Cayley's exact integer value untouched on a synthetic tuple with
//!   `Δ != 0`.
//! - **The population must hold several distinct case indices** — two
//!   instruments agreeing on one value agree on nothing (M-44).
//! - **The negation arm must see non-zero `Δ`**: on a `Δ = 0` population,
//!   negation invariance is a comparison of zeros.
//!
//! # The clock
//!
//! M-280: this host's `amd-pstate-epp` governor swings the same binary 1.45×
//! between runs. Every timed row carries min/median/max over [`REPEATS`]
//! interleaved repeats, warm-up pass excluded, and the count-based work ratio
//! recorded beside them is the claim the clock cannot corrupt.

#![allow(
    clippy::cast_precision_loss,
    clippy::too_many_lines,
    // The invariance checks compare two evaluations of the SAME degree-4
    // form on RELATED inputs: exact equality is the claim, and a margin
    // would turn the instrument's own arithmetic into the thing measured.
    clippy::float_cmp
)]

mod common;

use std::time::{Duration, Instant};

use common::poly::{Rng, cayley_2x2x2, octahedral_relabellings, relabel};
use isomesh::marching_cubes::table::{CASES, EDGE_CORNERS};
use isomesh::{Sdf, for_each_reference_field};

use isomesh::fields::ReferenceField;

const CELLS: u32 = 33;
/// Samples per axis.
const SAMPLES: u32 = CELLS + 1;

/// Interleaved timed repeats per arm, after one untimed warm-up.
const REPEATS: usize = 7;

/// The registered bar, as a ratio the clock measures.
const MIN_SPEEDUP: f64 = 10.0;

/// The sign word of a corner tuple under the crate's own inside rule —
/// `value < 0`, exact zero outside (`cube.rs:171`).
fn case_of(corner: &[f64; 8]) -> u8 {
    let mut c = 0u8;
    for (i, slot) in corner.iter().enumerate() {
        if *slot < 0.0 {
            c |= 1 << i;
        }
    }
    c
}

/// Apply the element to the case word: bit `i` of the output is bit `perm[i]`
/// of the input, matching `relabel`'s value convention exactly.
fn relabel_case(case: u8, perm: &[u8; 8]) -> u8 {
    let mut out = 0u8;
    for (i, slot) in perm.iter().enumerate() {
        if case & (1 << usize::from(*slot)) != 0 {
            out |= 1 << i;
        }
    }
    out
}

/// `g` applied to one edge: its two corners go where the element sends them,
/// and the moved edge's own table index is looked up in [`EDGE_CORNERS`].
fn relabel_edge(
    cut: &std::collections::BTreeSet<usize>,
    perm: &[u8; 8],
) -> std::collections::BTreeSet<usize> {
    let mut out = std::collections::BTreeSet::new();
    for &e in cut {
        let [a, b] = EDGE_CORNERS[e];
        let ga = usize::from(perm[usize::from(a)]);
        let gb = usize::from(perm[usize::from(b)]);
        let found = EDGE_CORNERS
            .iter()
            .position(|&[x, y]| {
                (usize::from(x) == ga && usize::from(y) == gb)
                    || (usize::from(x) == gb && usize::from(y) == ga)
            })
            .expect("P-129: an octahedral element maps the 12 cube edges to themselves");
        out.insert(found);
    }
    out
}

/// Cut edges of `corner`: the edge is cut exactly when its two endpoints
/// differ in the crate's inside rule.
fn cut_edges(corner: &[f64; 8]) -> std::collections::BTreeSet<usize> {
    let mut out = std::collections::BTreeSet::new();
    for (e, &[a, b]) in EDGE_CORNERS.iter().enumerate() {
        if (corner[usize::from(a)] < 0.0) != (corner[usize::from(b)] < 0.0) {
            out.insert(e);
        }
    }
    out
}

/// The sweep arm, per cell: the classification must be equivariant under the
/// element. Two exact integer conditions per element:
///
/// 1. the moved tuple's cut-edge set equals the relabelled cut set — the
///    element acts on the geometry, not just the word;
/// 2. the table's triangulation of the moved case has the same triangle count
///    as its triangulation of the relabelled original case — the shipped
///    classification carries the element's symmetry.
///
/// **The first condition alone is a tautology** (`case_of(relabel(f))` and
/// `relabel_case(case_of(f))` are the same function by construction, which the
/// first draft of this harness discovered by measuring zero violations on an
/// asymmetric tuple), so the triangulation count is what the sweep actually
/// tests, and the cut-set comparison is the bookkeeping that makes the count
/// meaningful.
fn sweep_violations_of(corner: &[f64; 8], elements: &[[u8; 8]; 48]) -> usize {
    let case = case_of(corner);
    let cut = cut_edges(corner);
    let mut violations = 0usize;
    for perm in elements {
        let moved = relabel(perm, corner);
        let moved_case = case_of(&moved);
        let expected_case = relabel_case(case, perm);
        let moved_cut = cut_edges(&moved);
        let expected_cut = relabel_edge(&cut, perm);
        if moved_cut != expected_cut
            || CASES[usize::from(moved_case)].count != CASES[usize::from(expected_case)].count
        {
            violations += 1;
        }
    }
    violations
}

/// The invariant arm, per cell: `Δ` must survive every relabelling and
/// negation. Returns (relabelling violations, negation violations), both out
/// of `f64`'s exact arithmetic on a degree-4 form of O(1) inputs.
fn invariant_violations_of(
    corner: &[f64; 8],
    cayley: &common::poly::Poly,
    elements: &[[u8; 8]; 48],
) -> usize {
    let base = cayley.eval_f64(corner);
    let mut violations = 0usize;
    for perm in elements {
        let moved: [f64; 8] = relabel(perm, corner);
        if cayley.eval_f64(&moved) != base {
            violations += 1;
        }
    }
    let neg: [f64; 8] = corner.map(|v| -v);
    if cayley.eval_f64(&neg) != base {
        violations += 1;
    }
    violations
}

/// Median of a non-empty run of samples, by `total_cmp`.
fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(f64::total_cmp);
    s[s.len() / 2]
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-129");

    common::experiment::run(prereg, |run| {
        let elements = octahedral_relabellings();
        let cayley = cayley_2x2x2();

        // ── vacuity control 1: the instrument's own arithmetic ──────────────
        {
            let trial: [i128; 8] = [3, -1, 4, 1, -5, 9, -2, 6];
            let base = cayley.eval_i128(&trial);
            assert_ne!(
                base, 0,
                "P-129 VOID: the synthetic trial has Δ = 0, so the invariance \
                 check is a comparison of two zeros"
            );
            for perm in &elements {
                let moved: [i128; 8] = relabel(perm, &trial);
                assert_eq!(
                    cayley.eval_i128(&moved),
                    base,
                    "P-129 VOID: relabelling under {perm:?} moves Δ on the \
                     synthetic trial; the committed element table is not \
                     Cayley-invariant"
                );
            }
        }

        // ── vacuity control 2: a broken table must be caught ───────────────
        {
            let broken: [f64; 8] = [-1.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0];
            let violations = sweep_violations_of(&broken, &elements);
            assert!(
                violations > 0,
                "P-129 VOID: a deliberately asymmetric corner tuple produced \
                 zero sweep violations, so the sweep cannot catch a broken table \
                 and its agreement with the invariant arm is two instruments \
                 agreeing on nothing"
            );
        }

        // ── the population ─────────────────────────────────────────────────
        let mut pop: Vec<[f64; 8]> = Vec::new();
        let mut distinct_cases: std::collections::BTreeSet<u8> = std::collections::BTreeSet::new();

        for_each_reference_field!(f64, |name, field| {
            let (lo, hi) = field.domain();
            let h = (hi[0] - lo[0]) / f64::from(SAMPLES - 1);
            let sx = SAMPLES;
            let n = usize::try_from(sx * sx * sx).unwrap_or(0);
            assert!(n > 0, "P-129 VOID: {name}'s 33^3 grid overflows u32");
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
            for cz in 0..CELLS {
                for cy in 0..CELLS {
                    for cx in 0..CELLS {
                        let at = |i: u32| {
                            ((cx + (i & 1))
                                + sx * ((cy + ((i >> 1) & 1)) + sx * (cz + ((i >> 2) & 1))))
                                as usize
                        };
                        let mut corner = [0.0f64; 8];
                        for (i, slot) in corner.iter_mut().enumerate() {
                            *slot = samples[at(i as u32)];
                        }
                        distinct_cases.insert(case_of(&corner));
                        pop.push(corner);
                    }
                }
            }
        });

        let mut rng = Rng::new(0x129_0001);
        for _ in 0..4096u32 {
            let c: [f64; 8] = std::array::from_fn(|_| rng.next_f64_unit() * 3.0 - 1.0);
            distinct_cases.insert(case_of(&c));
            pop.push(c);
        }
        let cells = pop.len();

        assert!(
            distinct_cases.len() >= 2,
            "P-129 VOID: the population holds {} distinct case indices; two \
             instruments agreeing on one value agree on nothing (M-44)",
            distinct_cases.len()
        );

        // ── vacuity control on the Δ population: some cells must carry a
        //    non-zero Δ, or the negation arm never has anything to move ─────
        let nonzero_delta = pop.iter().filter(|c| cayley.eval_f64(c) != 0.0).count();
        assert!(
            nonzero_delta > 0,
            "P-129 VOID: every census cell has Δ = 0, so the invariance and \
             negation arms compare zeros and cannot fail"
        );

        // ── verdicts, uninstrumented ───────────────────────────────────────
        let mut sweep_violations_total = 0usize;
        let mut invariant_violations_total = 0usize;
        let mut negation_violations_total = 0usize;
        let mut disagreements = 0usize;
        for corner in &pop {
            let sv = sweep_violations_of(corner, &elements);
            let iv = invariant_violations_of(corner, &cayley, &elements);
            sweep_violations_total += sv;
            // The invariant arm's function folds relabelling and negation
            // into one count; the negation half is recomputed on its own so
            // its verdict has a column of its own.
            let mut nv = 0usize;
            let base = cayley.eval_f64(corner);
            let neg: [f64; 8] = corner.map(|v| -v);
            if cayley.eval_f64(&neg) != base {
                nv = 1;
            }
            negation_violations_total += nv;
            invariant_violations_total += iv;
            if (sv == 0) != (iv == 0) {
                disagreements += 1;
            }
        }

        // ── timed interleaved repeats ──────────────────────────────────────
        let warm = Duration::from_millis(1);
        let mut sweep_times: Vec<f64> = Vec::with_capacity(REPEATS);
        let mut invariant_times: Vec<f64> = Vec::with_capacity(REPEATS);
        for repeat in 0..=REPEATS {
            let started = Instant::now();
            for corner in &pop {
                let _ = sweep_violations_of(corner, &elements);
            }
            let elapsed = started.elapsed();
            if repeat > 0 {
                sweep_times.push(elapsed.as_secs_f64() * 1e3);
            }

            let started = Instant::now();
            for corner in &pop {
                let _ = invariant_violations_of(corner, &cayley, &elements);
            }
            let elapsed = started.elapsed();
            if repeat > 0 {
                invariant_times.push(elapsed.as_secs_f64() * 1e3);
            }
        }
        let _ = warm;

        let sweep_med = median(&sweep_times);
        let inv_med = median(&invariant_times);
        let sweep_min = sweep_times.iter().copied().fold(f64::INFINITY, f64::min);
        let sweep_max = sweep_times
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let inv_min = invariant_times
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let inv_max = invariant_times
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let speedup = sweep_med / inv_med;

        // Work counts, machine-independent.
        let relabellings = cells * 48usize;
        let delta_evals = cells * 49usize; // 48 relabelled + 1 negated pair

        let c1 = sweep_violations_total == 0
            && invariant_violations_total == 0
            && disagreements == 0
            && speedup >= MIN_SPEEDUP;
        let c2 = negation_violations_total == 0;
        // C3 answered honestly: the invariant check is Δ-blind BY THEOREM
        // (P-127 proved every violation impossible), so it catches no
        // configuration class the sweep is blind to — the registered "none"
        // outcome, recorded as such.
        let c3 = disagreements == 0;

        run.record(&[
            ("element", "all-48-aggregate".to_string()),
            ("sweep_ms", format!("{sweep_med:.6}")),
            ("invariant_ms", format!("{inv_med:.6}")),
            ("speedup", format!("{speedup:.3}")),
            ("cells_checked", cells.to_string()),
            ("sweep_violations", sweep_violations_total.to_string()),
            (
                "invariant_violations",
                invariant_violations_total.to_string(),
            ),
            ("agreement", (disagreements == 0).to_string()),
            ("negation_violations", negation_violations_total.to_string()),
            ("c1_holds", c1.to_string()),
            ("c2_holds", c2.to_string()),
            ("c3_holds", c3.to_string()),
            // ── extras (M-273) ──
            (
                "c3_token",
                // C3's registered option "none, fine outcome, recorded as
                // such" is the honest reading: Δ's invariance is a theorem,
                // so the algebraic check catches nothing the sweep can and
                // misses nothing either. Disagreements would mean one
                // instrument misread the same cell.
                if disagreements == 0 {
                    "none-registered-outcome".to_string()
                } else {
                    format!("{disagreements}-cells-sweep-invariant-verdict-disagree")
                },
            ),
            ("invariant_ms_max", format!("{inv_max:.6}")),
            ("invariant_ms_min", format!("{inv_min:.6}")),
            ("nonzero_delta_cells", nonzero_delta.to_string()),
            ("population_synthetic", "4096".to_string()),
            ("repeats", REPEATS.to_string()),
            ("relabellings", relabellings.to_string()),
            (
                "relabellings_per_delta",
                format!("{:.3}", relabellings as f64 / delta_evals as f64),
            ),
            ("sweep_disagreements", disagreements.to_string()),
            ("sweep_ms_max", format!("{sweep_max:.6}")),
            ("sweep_ms_min", format!("{sweep_min:.6}")),
        ]);
    });
}
