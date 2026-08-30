//! **P-99 — the under-resolution witness on the metric `✗53`'s own data supports.**
//!
//! Ticket: R-099. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p99
//! ```
//!
//! Writes `docs/experiments/p-99.csv`.
//!
//! # What is inherited and what is new
//!
//! Everything about the *measurement* is `P-66`'s and is copied rather than
//! re-invented: the predicate (`k` samples of `∇f · ê` along a grid edge, flagged
//! non-monotonic when two sampled projections disagree in sign), the eight
//! reference fields, the three resolutions, the five `k`, and the oracle —
//! `subgrid::roots::all_roots` at [`ORACLE_SAMPLES`] intervals per edge. The
//! sample points are computed by the same expression in the same order, so every
//! column this shares with `p-66.csv` must reproduce it bit for bit; that
//! agreement is the reproduction control and it is checked by eye against the
//! committed file rather than asserted here, because `p-66.csv` is another
//! experiment's artefact.
//!
//! What is new is the **denominator**. `✗53` scored the witness over *all* grid
//! edges, and that is why `thin_plate` — a plate in a mostly-empty box — did not
//! top the ranking: almost every edge in that box is far from the plate, where
//! the field is smooth. The all-edges rate measures how much of the **volume**
//! has turning gradient. This row measures the same numerator against the edges
//! that carry the **surface**, and reports both rates on every row, because the
//! denominator was the defect.
//!
//! # Convergence, defined, because an unstated criterion is not a clause
//!
//! A field **converges** at a fixed `k` when its `rate_single_root` sequence over
//! the three resolutions 17³ → 33³ → 65³ is **non-increasing**: each refinement
//! weakly reduces the rate. It **fails to converge** when the sequence rises at
//! either step. `converges_strictly` is the stronger form — strictly decreasing
//! at both steps — and it is recorded separately so that a field whose rate is
//! identically zero is not silently counted as evidence of convergence.
//! `rate_ever_nonzero` is on the row for the same reason.
//!
//! # Ranking, defined, because ties would make "first" ambiguous
//!
//! `rank = 1 + (number of fields at this resolution and `k` with a strictly
//! greater rate)`, so tied fields share the better rank. That makes "ranks
//! first" mean "no field scores strictly higher", which would be satisfied by a
//! tie — so `rank1_strict_single_root` records whether the top rate is strictly
//! above the second, and it is a control rather than a decoration.
//!
//! # SHARE
//!
//! No clause here is a ratio of a runtime total, so `✗51`'s arithmetic does not
//! apply and there is no share to recompute against a speedup bar. The
//! registration says so itself: *"this is a rate; it moves nothing"*. All three
//! clauses are comparisons of exactly counted integer populations — a rank, a
//! sequence direction, and two monotonicities — and each is reachable given that
//! the sweep has at least two resolutions and two `k`, which it has three and
//! five of.
//!
//! One reachability caveat that is **not** arithmetic and is worth stating
//! before the run: the metric C1 endorses is not oracle-free. Its denominator is
//! *"edges the root finder reports with exactly one root"*, so computing it needs
//! `all_roots`, exactly like the false-negative count `✗53` set aside as
//! *"offline"*. The oracle-free proxy is the **sign-changing** edge — an edge
//! whose endpoints straddle the isosurface, which Marching Cubes already
//! classifies for free — and it is recorded alongside as `rate_sign_change` so
//! the deployability of the LOD signal is on the artefact rather than assumed.
//!
//! # The vacuity control, and its caveat
//!
//! `oracle_samples` is `✗53`'s column, inherited at **128 intervals**, and so is
//! the caveat that comes with it: `all_roots` cannot see two roots closer
//! together than one interval, i.e. 1/128 of an edge. **A root pair inside one
//! interval is invisible to both instruments — it cannot produce a false
//! negative here, and it also cannot be ruled out.** What the control asserts is
//! the thing that makes C3's numbers measurements rather than empty sets: at 128
//! intervals the oracle finds a non-empty multi-root population, and the
//! false-negative count at the coarsest `k` is non-zero — so a false-negative
//! count of zero at `k = 17` is a zero that could have been non-zero.
//!
//! # Which edges
//!
//! Every axis-aligned **grid edge** of the sample lattice, `3 · (n-1) · n · n` of
//! them, which is the edge set Marching Cubes interpolates along. Because every
//! edge is axis-aligned, `ê` is a basis vector and `∇f · ê` is the single
//! gradient component `g[axis]` — the same number `P-66` obtained by normalising
//! `(h, 0, 0)` and taking a three-term dot product, without the division or the
//! zero-length branch.

mod common;

use isomesh::Sdf;
use isomesh::fields::ReferenceField;
use isomesh::for_each_reference_field;
use isomesh::subgrid::roots::all_roots;

/// Samples per axis. `P-66`'s three, inherited.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// Directional-derivative samples per edge. `P-66`'s five, inherited, because
/// C3 is a statement about how two counts move across exactly this sweep.
const KS: [u32; 5] = [2, 3, 5, 9, 17];

/// The `k` C1 and C2 are stated at.
///
/// `✗53`'s quoted `thin_plate` sequence — 0.3889, 0.4412, 0.4697 — is the
/// `k = 5` column of `p-66.csv`, so that is the `k` the two clauses inherit.
/// Ranks and convergence are computed at **every** `k` and recorded on every
/// row, so a reader can see whether either clause depends on the choice.
const K_REGISTERED: u32 = 5;

/// Intervals the oracle divides each edge into. The registered vacuity control.
///
/// 128 rather than `M-94`'s 1000, for `P-66`'s reason: the sweep walks 811,200
/// edges at 65³ and 1000 intervals each would be 811 million field evaluations
/// per field per resolution. It resolves a root pair separated by more than
/// 1/128 of a cell and nothing finer, and that limit is a recorded column rather
/// than a footnote.
const ORACLE_SAMPLES: u32 = 128;

/// One `(field, resolution, k)` row.
struct Row {
    field: &'static str,
    resolution: u32,
    k: u32,
    /// All axis-aligned grid edges. `rate_all_edges`' denominator.
    edges: u64,
    /// Edges the oracle reports with exactly one root. `rate_single_root`'s
    /// denominator, and the whole subject of C1.
    single_root: u64,
    /// Edges the oracle reports with more than one root. C3's population.
    multi_root: u64,
    /// Edges the witness flagged, over the whole grid. `rate_all_edges`'
    /// numerator — `✗53`'s `flagged_non_monotonic`.
    flagged: u64,
    /// Single-root edges the witness flagged.
    ///
    /// This is `✗53`'s `false_positives` column under P-99's name. Both are
    /// registered and both are written, from this one field: they are the same
    /// count, and pretending otherwise by measuring it twice would invent an
    /// agreement instead of recording an identity.
    non_monotonic_single_root: u64,
    /// Multi-root edges the witness called monotonic.
    false_negatives: u64,
    /// Edges whose endpoints straddle the isosurface. Oracle-free.
    sign_change: u64,
    /// Sign-changing edges the witness flagged.
    non_monotonic_sign_change: u64,
    /// `non_monotonic_single_root` split by the edge's axis. The mechanism
    /// instrument: it says *where* the flagged edges are, not just how many.
    nm_axis: [u64; 3],
    /// Per `(field, resolution)`: edges flagged at one `k` and **not** flagged
    /// at the next larger `k`.
    ///
    /// C3's second half as a set inclusion rather than a count comparison. Two
    /// counts can both be non-decreasing while individual edges swap in and out,
    /// and that would be a predicate violation the count check cannot see.
    flag_dropped: u64,
    /// Per `(field, resolution)`: edges where the endpoint sign change disagrees
    /// with the parity of the oracle's root count.
    ///
    /// The oracle counts `inside`-transitions over its intervals, so its count
    /// is odd exactly when the endpoints straddle. This cross-checks that the
    /// two instruments are being read with the same convention.
    parity_violations: u64,
}

impl Row {
    fn new(field: &'static str, resolution: u32, k: u32) -> Self {
        Self {
            field,
            resolution,
            k,
            edges: 0,
            single_root: 0,
            multi_root: 0,
            flagged: 0,
            non_monotonic_single_root: 0,
            false_negatives: 0,
            sign_change: 0,
            non_monotonic_sign_change: 0,
            nm_axis: [0; 3],
            flag_dropped: 0,
            parity_violations: 0,
        }
    }
}

/// Is this axis-aligned edge monotonic under the witness, at `k` samples?
///
/// `k` points **inclusive of both endpoints**, so `k = 2` is the endpoints alone.
/// The five `k` are recomputed independently rather than derived from the
/// finest — every coarser sample set here happens to be a bit-exact subset of
/// the `k = 17` set, so subsetting would be 2.1× cheaper *and* would make C3's
/// second half true by construction of this function instead of by construction
/// of the predicate.
///
/// A sampled projection of **exactly zero** is not a sign disagreement, which is
/// `P-66`'s rule and the crate's convention that zero is not inside: treating it
/// as either sign would make the answer depend on comparison order.
fn is_monotonic<S: Sdf<Scalar = f64>>(
    sdf: &S,
    from: [f64; 3],
    to: [f64; 3],
    axis: usize,
    k: u32,
) -> bool {
    let mut seen_positive = false;
    let mut seen_negative = false;
    for i in 0..k {
        let t = f64::from(i) / f64::from(k - 1);
        let p = [
            from[0] + (to[0] - from[0]) * t,
            from[1] + (to[1] - from[1]) * t,
            from[2] + (to[2] - from[2]) * t,
        ];
        // `ê` is a basis vector, so the dot product is one component.
        let proj = sdf.gradient(p)[axis];
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
/// sweep: re-running it per `k` could not change an answer, and two oracle calls
/// on one edge that disagreed would be a determinism bug this harness is not the
/// place to find.
fn measure<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    field_name: &'static str,
    resolution: u32,
) -> Vec<Row> {
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(resolution - 1);
    let at = |x: u32, y: u32, z: u32| -> [f64; 3] {
        [
            lo[0] + h * f64::from(x),
            lo[1] + h * f64::from(y),
            lo[2] + h * f64::from(z),
        ]
    };

    let mut rows: Vec<Row> = KS
        .iter()
        .map(|k| Row::new(field_name, resolution, *k))
        .collect();

    let mut roots: Vec<f64> = Vec::with_capacity(8);
    let mut flags = [false; KS.len()];
    let mut flag_dropped = 0u64;
    let mut parity_violations = 0u64;
    let n = resolution;
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

                    // The endpoints are sampled through the oracle's own
                    // parameterisation — `from + (to - from) * t` at `t = 0` and
                    // `t = 1`, the same expression in the same order — so the
                    // parity control below compares two instruments and not two
                    // roundings.
                    let ends = [0.0f64, 1.0].map(|t| {
                        field.sample([
                            a[0] + (b[0] - a[0]) * t,
                            a[1] + (b[1] - a[1]) * t,
                            a[2] + (b[2] - a[2]) * t,
                        ])
                    });
                    let sign_change = (ends[0] < 0.0) != (ends[1] < 0.0);
                    if sign_change != (count % 2 == 1) {
                        parity_violations += 1;
                    }

                    for (i, row) in rows.iter_mut().enumerate() {
                        let monotonic = is_monotonic(field, a, b, axis, row.k);
                        flags[i] = !monotonic;
                        row.edges += 1;
                        if !monotonic {
                            row.flagged += 1;
                        }
                        if sign_change {
                            row.sign_change += 1;
                            if !monotonic {
                                row.non_monotonic_sign_change += 1;
                            }
                        }
                        if count > 1 {
                            row.multi_root += 1;
                            if monotonic {
                                row.false_negatives += 1;
                            }
                        } else if count == 1 {
                            row.single_root += 1;
                            if !monotonic {
                                row.non_monotonic_single_root += 1;
                                row.nm_axis[axis] += 1;
                            }
                        }
                    }
                    for i in 1..KS.len() {
                        if flags[i - 1] && !flags[i] {
                            flag_dropped += 1;
                        }
                    }
                }
            }
        }
    }

    for row in &mut rows {
        row.flag_dropped = flag_dropped;
        row.parity_violations = parity_violations;
    }
    rows
}

/// A field's `rate_single_root` sequence over the three resolutions at one `k`.
struct Sequence {
    field: &'static str,
    k: u32,
    rate: [f64; 3],
    /// Non-increasing at both steps. The registered definition of convergence.
    converges: bool,
    /// Strictly decreasing at both steps.
    strictly: bool,
    /// Non-zero somewhere, so `converges` is not a statement about a flat zero.
    ever_nonzero: bool,
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-99");
    let mut rows: Vec<Row> = Vec::new();

    for resolution in RESOLUTIONS {
        for_each_reference_field!(f64, |name, field| {
            rows.extend(measure(&field, name, resolution));
        });
    }

    let mut fields: Vec<&'static str> = Vec::new();
    for r in &rows {
        if !fields.contains(&r.field) {
            fields.push(r.field);
        }
    }

    let find = |field: &str, resolution: u32, k: u32| -> &Row {
        rows.iter()
            .find(|r| r.field == field && r.resolution == resolution && r.k == k)
            .expect("every (field, resolution, k) was swept")
    };
    let rate_single_root = |r: &Row| r.non_monotonic_single_root as f64 / r.single_root as f64;
    let rate_all_edges = |r: &Row| r.flagged as f64 / r.edges as f64;
    let rate_sign_change = |r: &Row| r.non_monotonic_sign_change as f64 / r.sign_change as f64;

    // ── the two rankings, side by side, per (resolution, k) ──────────────────
    //
    // `rank = 1 + #{strictly greater}`, so ties share the better rank.
    let cohort =
        |resolution: u32, k: u32, rate: &dyn Fn(&Row) -> f64| -> Vec<(&'static str, f64)> {
            fields
                .iter()
                .map(|f| (*f, rate(find(f, resolution, k))))
                .collect()
        };
    let rank_in = |c: &[(&'static str, f64)], field: &str| -> usize {
        let v = c
            .iter()
            .find(|(g, _)| *g == field)
            .map(|(_, v)| *v)
            .expect("field in cohort");
        1 + c.iter().filter(|(_, w)| *w > v).count()
    };
    // Is the top of this cohort strictly above the second? Without it, "ranks
    // first" could be satisfied by a tie.
    let rank1_strict = |c: &[(&'static str, f64)]| -> bool {
        let mut v: Vec<f64> = c.iter().map(|(_, w)| *w).collect();
        v.sort_unstable_by(|a, b| b.partial_cmp(a).expect("finite rates"));
        v[0] > v[1]
    };

    // ── convergence, per (field, k) ──────────────────────────────────────────
    let mut sequences: Vec<Sequence> = Vec::new();
    for f in &fields {
        for k in KS {
            let rate = RESOLUTIONS.map(|n| rate_single_root(find(f, n, k)));
            sequences.push(Sequence {
                field: f,
                k,
                converges: rate[1] <= rate[0] && rate[2] <= rate[1],
                strictly: rate[1] < rate[0] && rate[2] < rate[1],
                ever_nonzero: rate.iter().any(|v| *v > 0.0),
                rate,
            });
        }
    }
    let seq_of = |field: &str, k: u32| -> &Sequence {
        sequences
            .iter()
            .find(|s| s.field == field && s.k == k)
            .expect("every (field, k) has a sequence")
    };

    // ── C3's two monotonicities, aggregate and per (field, resolution) ───────
    let fn_total: Vec<u64> = KS
        .iter()
        .map(|k| {
            rows.iter()
                .filter(|r| r.k == *k)
                .map(|r| r.false_negatives)
                .sum()
        })
        .collect();
    let fp_total: Vec<u64> = KS
        .iter()
        .map(|k| {
            rows.iter()
                .filter(|r| r.k == *k)
                .map(|r| r.non_monotonic_single_root)
                .sum()
        })
        .collect();
    let fn_falls_aggregate = fn_total.windows(2).all(|w| w[1] <= w[0]);
    let fn_falls_strictly = fn_total.windows(2).all(|w| w[1] < w[0]);
    let fp_rises_aggregate = fp_total.windows(2).all(|w| w[1] >= w[0]);

    let fp_seq = |field: &str, resolution: u32| -> Vec<u64> {
        KS.iter()
            .map(|k| find(field, resolution, *k).non_monotonic_single_root)
            .collect()
    };
    let fn_seq = |field: &str, resolution: u32| -> Vec<u64> {
        KS.iter()
            .map(|k| find(field, resolution, *k).false_negatives)
            .collect()
    };
    let fp_non_decreasing = |field: &str, n: u32| fp_seq(field, n).windows(2).all(|w| w[1] >= w[0]);
    let fn_non_increasing = |field: &str, n: u32| fn_seq(field, n).windows(2).all(|w| w[1] <= w[0]);
    // Did the count move at all across the sweep? A constant sequence satisfies
    // "non-decreasing" without ever having been able to fail it.
    let fp_moves = |field: &str, n: u32| fp_seq(field, n).windows(2).any(|w| w[1] > w[0]);

    // ── the table: both denominators on every row ───────────────────────────
    println!(
        "{:>15} {:>4} {:>3} {:>9} {:>7} {:>7} {:>7} {:>9} {:>8} {:>4} {:>8} {:>4} {:>8} {:>4} \
         {:>4} {:>4}",
        "field",
        "n",
        "k",
        "edges",
        "1-root",
        "n-root",
        "signch",
        "nm-1root",
        "rate-1r",
        "rk",
        "rate-all",
        "rk",
        "rate-sc",
        "rk",
        "FN",
        "conv"
    );
    for r in &rows {
        let c_sr = cohort(r.resolution, r.k, &rate_single_root);
        let c_all = cohort(r.resolution, r.k, &rate_all_edges);
        let c_sc = cohort(r.resolution, r.k, &rate_sign_change);
        println!(
            "{:>15} {:>4} {:>3} {:>9} {:>7} {:>7} {:>7} {:>9} {:>8.4} {:>4} {:>8.4} {:>4} \
             {:>8.4} {:>4} {:>4} {:>4}",
            r.field,
            r.resolution,
            r.k,
            r.edges,
            r.single_root,
            r.multi_root,
            r.sign_change,
            r.non_monotonic_single_root,
            rate_single_root(r),
            rank_in(&c_sr, r.field),
            rate_all_edges(r),
            rank_in(&c_all, r.field),
            rate_sign_change(r),
            rank_in(&c_sc, r.field),
            r.false_negatives,
            if seq_of(r.field, r.k).converges {
                "y"
            } else {
                "N"
            }
        );
    }

    // ── controls, all before the artefact is written ─────────────────────────
    //
    // The registered vacuity control is `oracle_samples`. Pinning the constant
    // is the cheap half; the half that matters is that at 128 intervals the
    // oracle finds a population C3 can be wrong about.
    assert_eq!(
        ORACLE_SAMPLES, 128,
        "the registration names 128 intervals as the vacuity control, and the recorded column \
         has to be that number"
    );
    let multi_total: u64 = rows
        .iter()
        .filter(|r| r.k == K_REGISTERED)
        .map(|r| r.multi_root)
        .sum();
    assert!(
        multi_total > 0,
        "VACUOUS: the oracle found no multi-root edge anywhere at {ORACLE_SAMPLES} intervals, so \
         C3's false-negative counts are a sequence of zeros over an empty population"
    );
    assert!(
        fn_total[0] > 0,
        "VACUOUS: the false-negative count is already zero at k = {}, so C3's fall with k is a \
         zero that could not have been non-zero",
        KS[0]
    );
    for r in &rows {
        assert!(
            r.single_root > 0,
            "VACUOUS: {} at {}³ has no single-root edge, so C1's rate has no denominator",
            r.field,
            r.resolution
        );
        assert!(
            r.sign_change >= r.single_root,
            "{} at {}³: every single-root edge straddles the surface, so the oracle-free \
             denominator cannot be smaller than the oracle's",
            r.field,
            r.resolution
        );
    }
    let flagged_total: u64 = rows
        .iter()
        .filter(|r| r.k == K_REGISTERED)
        .map(|r| r.flagged)
        .sum();
    let unflagged_total: u64 = rows
        .iter()
        .filter(|r| r.k == K_REGISTERED)
        .map(|r| r.edges - r.flagged)
        .sum();
    assert!(
        flagged_total > 0,
        "VACUOUS: the witness flagged nothing at k = {K_REGISTERED}, so no rate can rank anything"
    );
    assert!(
        unflagged_total > 0,
        "VACUOUS: the witness flagged EVERY edge at k = {K_REGISTERED}, so every rate is 1 and \
         every rank is a tie"
    );
    for n in [33, 65] {
        let c = cohort(n, K_REGISTERED, &rate_single_root);
        assert!(
            rank1_strict(&c),
            "VACUOUS: the top two single-root rates at {n}³ are equal, so C1's \"ranks first\" is \
             a tie rather than a result"
        );
    }
    assert!(
        rows.iter()
            .any(|r| r.k == K_REGISTERED && fp_moves(r.field, r.resolution)),
        "VACUOUS: no field's flagged single-root count moves across the k sweep, so C3's \
         non-decreasing half is satisfied by a constant"
    );
    assert!(
        sequences
            .iter()
            .any(|s| s.k == K_REGISTERED && s.converges && s.ever_nonzero),
        "VACUOUS: every converging field converges at a flat zero, so C2's \"converges on the \
         other seven\" is a statement about empty numerators"
    );

    // ── verdicts ─────────────────────────────────────────────────────────────
    //
    // C1: the single-root rate ranks `thin_plate` first at 33³ and 65³, at the
    // k the registration's quoted numbers come from.
    let mut c1 = true;
    for n in [33, 65] {
        let c_sr = cohort(n, K_REGISTERED, &rate_single_root);
        let c_all = cohort(n, K_REGISTERED, &rate_all_edges);
        let r_sr = rank_in(&c_sr, "thin_plate");
        let r_all = rank_in(&c_all, "thin_plate");
        println!(
            "C1 thin_plate at {n}³, k = {K_REGISTERED}: rank_single_root = {r_sr}, rank_all_edges = {r_all}"
        );
        c1 &= r_sr == 1;
    }

    // C2: fails to converge on `thin_plate`, converges on the other seven.
    let mut c2 = true;
    for f in &fields {
        let s = seq_of(f, K_REGISTERED);
        let want = *f != "thin_plate";
        println!(
            "C2 {f} at k = {K_REGISTERED}: {} -> {}{}",
            s.rate
                .iter()
                .map(|v| format!("{v:.4}"))
                .collect::<Vec<_>>()
                .join(" "),
            if s.converges {
                if s.strictly {
                    "converges (strictly)"
                } else {
                    "converges (not strictly)"
                }
            } else {
                "DOES NOT CONVERGE"
            },
            if s.ever_nonzero {
                ""
            } else {
                " [rate is a flat zero]"
            }
        );
        c2 &= s.converges == want;
    }

    // C3: false negatives fall with k, false positives do not.
    let mut c3 = fn_falls_aggregate && fp_rises_aggregate;
    for f in &fields {
        for n in RESOLUTIONS {
            c3 &= fp_non_decreasing(f, n) && fn_non_increasing(f, n);
        }
    }
    println!(
        "C3 false negatives over k = {KS:?}: {fn_total:?} -> {}",
        if fn_falls_aggregate {
            if fn_falls_strictly {
                "non-increasing (and strictly falling)"
            } else {
                "non-increasing"
            }
        } else {
            "RISES"
        }
    );
    println!(
        "C3 flagged single-root edges over k = {KS:?}: {fp_total:?} -> {}",
        if fp_rises_aggregate {
            "non-decreasing"
        } else {
            "FALLS"
        }
    );

    println!(
        "\nC1 single-root rate ranks thin_plate first at 33³ and 65³ -> {}",
        if c1 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C2 fails to converge on thin_plate and converges on the other seven -> {}",
        if c2 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C3 false negatives non-increasing in k, flagged single-root non-decreasing -> {}",
        if c3 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "\nSCOPE: the oracle divides each edge into {ORACLE_SAMPLES} intervals and cannot resolve \
         a root pair closer together than one of them. Such a pair is invisible to BOTH \
         instruments: it cannot produce a false negative here, and it cannot be ruled out."
    );

    common::experiment::run(prereg, |run| {
        for r in &rows {
            let c_sr = cohort(r.resolution, r.k, &rate_single_root);
            let c_all = cohort(r.resolution, r.k, &rate_all_edges);
            let c_sc = cohort(r.resolution, r.k, &rate_sign_change);
            let s = seq_of(r.field, r.k);
            run.record(&[
                ("field", r.field.to_string()),
                ("resolution", r.resolution.to_string()),
                ("k", r.k.to_string()),
                ("single_root_edges", r.single_root.to_string()),
                ("multi_root_edges", r.multi_root.to_string()),
                (
                    "non_monotonic_single_root",
                    r.non_monotonic_single_root.to_string(),
                ),
                ("rate_single_root", format!("{:.6}", rate_single_root(r))),
                ("rate_all_edges", format!("{:.6}", rate_all_edges(r))),
                ("rank_single_root", rank_in(&c_sr, r.field).to_string()),
                ("rank_all_edges", rank_in(&c_all, r.field).to_string()),
                ("converges", s.converges.to_string()),
                ("false_negatives", r.false_negatives.to_string()),
                ("false_positives", r.non_monotonic_single_root.to_string()),
                ("oracle_samples", ORACLE_SAMPLES.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // Extras. The denominators, the oracle-free proxy, the
                // mechanism, and the controls.
                ("edges", r.edges.to_string()),
                ("flagged_non_monotonic", r.flagged.to_string()),
                ("sign_change_edges", r.sign_change.to_string()),
                (
                    "non_monotonic_sign_change",
                    r.non_monotonic_sign_change.to_string(),
                ),
                ("rate_sign_change", format!("{:.6}", rate_sign_change(r))),
                ("rank_sign_change", rank_in(&c_sc, r.field).to_string()),
                ("converges_strictly", s.strictly.to_string()),
                ("rate_ever_nonzero", s.ever_nonzero.to_string()),
                ("rank1_strict_single_root", rank1_strict(&c_sr).to_string()),
                ("nm_single_root_axis_x", r.nm_axis[0].to_string()),
                ("nm_single_root_axis_y", r.nm_axis[1].to_string()),
                ("nm_single_root_axis_z", r.nm_axis[2].to_string()),
                ("flag_dropped_at_higher_k", r.flag_dropped.to_string()),
                ("oracle_parity_violations", r.parity_violations.to_string()),
                (
                    "fp_non_decreasing_with_k",
                    fp_non_decreasing(r.field, r.resolution).to_string(),
                ),
                (
                    "fn_non_increasing_with_k",
                    fn_non_increasing(r.field, r.resolution).to_string(),
                ),
                (
                    "fp_moves_with_k",
                    fp_moves(r.field, r.resolution).to_string(),
                ),
            ]);
        }
    });

    // ── after the artefact: the two integrity checks ─────────────────────────
    //
    // These are predicate results rather than vacuity gates, so they run after
    // the CSV is on disk: a violation is a finding someone has to be able to
    // read the numbers for.
    let parity: u64 = rows
        .iter()
        .filter(|r| r.k == K_REGISTERED)
        .map(|r| r.parity_violations)
        .sum();
    assert_eq!(
        parity, 0,
        "the endpoint sign change disagrees with the parity of the oracle's root count on \
         {parity} edges, so the two instruments are not being read with the same convention"
    );
    let dropped: u64 = rows
        .iter()
        .filter(|r| r.k == K_REGISTERED)
        .map(|r| r.flag_dropped)
        .sum();
    assert_eq!(
        dropped, 0,
        "{dropped} edges are flagged at one k and unflagged at the next larger one, which the \
         predicate cannot do: adding a sample point can only add an opportunity to disagree"
    );
}
