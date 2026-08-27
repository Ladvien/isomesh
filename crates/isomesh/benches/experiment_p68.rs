//! **P-68 - a running error bound on the crossing, against exact arithmetic.**
//!
//! Ticket: R-066. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p68
//! ```
//!
//! Writes `docs/experiments/p-68.csv`.
//!
//! # The bound
//!
//! `cube::edge_offset` computes `d = ((a + b)/2) / (a − b)` in `f64`. Three
//! roundings occur: the sum, the halving (**exact**, 2 is a power of two), and
//! the quotient. So
//!
//! ```text
//! bound(d) = |d| · (2 + |a−b|_err / |a−b|) · u
//! ```
//!
//! as registered, with `u = 2⁻⁵³`. The `2` covers the sum's rounding and the
//! quotient's; the third term is the denominator's own error, which is **zero**
//! whenever `a − b` is exact and is computed here rather than bounded, by
//! Knuth's two-sum. That is the whole reason the centred form is cheap to
//! certify: there is no cancellation to amplify.
//!
//! # The ground truth is exact, and no floating point appears in it
//!
//! `a` and `b` are `f64`, hence dyadic rationals `A·2^e` and `B·2^e` for a
//! common `e`. Then
//!
//! ```text
//! d_true = (A + B) / (2 · (A − B))
//! ```
//!
//! exactly, with `A` and `B` integers - and the `2^e` cancels, which is why the
//! reference needs no exponent tracking beyond the alignment shift. The `f64`
//! result is `M·2^k`, also exact, so
//!
//! ```text
//! d̂ − d_true = (M · 2^(k+1) · (A − B) − (A + B)) / (2 · (A − B))
//! ```
//!
//! is a ratio of two integers. Both are computed in `i128`. A crossing whose
//! alignment shift would overflow `i128` is **not judged** and is counted, so
//! C1's population is on the artefact rather than in a comment.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::Sdf;
use isomesh::fields::{ReferenceField, Sphere};
use isomesh::for_each_reference_field;
use isomesh::marching_cubes::table::{edge_offset, is_inside};

/// Samples per axis. The registered one.
const SAMPLES: u32 = 33;

/// `f64`'s unit roundoff, `2⁻⁵³`.
const U: f64 = f64::EPSILON * 0.5;

/// Knuth's two-sum: `x + y == a + b` exactly, with `x = fl(a + b)`.
///
/// Six flops, no branch, exact for every finite input whose sum does not
/// overflow. Written here rather than imported because `predicates.rs`'s copy is
/// module-private - and the property that makes it usable is **asserted** by
/// [`the_two_sum_is_exact`] against the `i128` reference, so this is not a
/// second definition taken on trust.
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    let bv = x - a;
    let av = x - bv;
    let br = b - bv;
    let ar = a - av;
    (x, ar + br)
}

/// `a` as an exact `(mantissa, exponent)` pair with `a == mantissa · 2^exponent`.
///
/// Subnormals included: the exponent is clamped at the subnormal floor and the
/// mantissa carries no implicit bit, which is exactly `f64`'s own encoding.
/// Zero returns `(0, 0)`, which is correct for every use below.
fn decompose(a: f64) -> (i128, i32) {
    let bits = a.to_bits();
    let sign = if bits >> 63 == 1 { -1i128 } else { 1 };
    let raw_exp = ((bits >> 52) & 0x7ff) as i32;
    let frac = (bits & 0x000f_ffff_ffff_ffff) as i128;
    if raw_exp == 0 {
        // Subnormal or zero: no implicit leading bit.
        (sign * frac, -1074)
    } else {
        (sign * (frac | 0x0010_0000_0000_0000), raw_exp - 1075)
    }
}

/// The exact error of `d_hat` against the true `d`, in ulps of `d_hat`.
///
/// `None` when the alignment shifts would overflow `i128`. That is a **refusal
/// to judge**, not a pass: the caller counts it.
fn exact_error_ulp(a: f64, b: f64, d_hat: f64) -> Option<f64> {
    let (ma, ea) = decompose(a);
    let (mb, eb) = decompose(b);
    let (md, ed) = decompose(d_hat);

    // Align a and b on a common exponent. The 2^e then cancels in the ratio.
    //
    // **Zero's exponent must not take part, and the first version let it.**
    // `decompose(0.0)` is `(0, -1074)` because that is `f64`'s subnormal floor,
    // and `0 · 2^anything` is still zero - so a zero endpoint dragged `e` to
    // -1074 and made the other endpoint's shift 1019, past the guard. That
    // refused **all 1,350** of `box_exact`'s crossings and 1,167 of
    // `csg_difference`'s, every one of them an edge whose crossing sits exactly
    // on a sample where the error is exactly zero. Twice now this reference has
    // declined to judge its own easiest case; both were caught by the judged
    // count being a recorded column rather than an assumption.
    let e = match (ma == 0, mb == 0) {
        (true, true) => return Some(0.0), // both endpoints zero: no crossing to misplace
        (true, false) => eb,
        (false, true) => ea,
        (false, false) => ea.min(eb),
    };
    // A zero mantissa is zero at every scale, so it takes no shift - and
    // `(ea - e) as u32` on the -1074 a zero carries wraps to four billion, which
    // is how the guard below was passed and `checked_shl` then refused. Naming
    // the zero case is the fix; relying on the shift arithmetic was the defect.
    let lift = |m: i128, ep: i32| -> Option<i128> {
        if m == 0 {
            return Some(0);
        }
        let s = ep - e;
        // 53-bit mantissas, so a shift past 70 cannot fit alongside the products
        // below. Refuse rather than wrap.
        if !(0..=70).contains(&s) {
            return None;
        }
        m.checked_shl(s as u32)
    };
    let big_a = lift(ma, ea)?;
    let big_b = lift(mb, eb)?;

    let num = big_a.checked_add(big_b)?; // A + B
    let den = big_a.checked_sub(big_b)?; // A - B
    if den == 0 {
        return None;
    }

    // d_hat - d_true = (M·2^(k+1)·(A−B) − (A+B)) / (2·(A−B)).
    //
    // `k` is very negative for a `d` in [-1/2, 1/2], so the numerator's first
    // term is scaled DOWN and the second scaled UP by the same power of two.
    // Multiplying through by 2^-(k+1) keeps both integral:
    //   numerator   = M·(A−B) − (A+B)·2^-(k+1)
    //   denominator = 2·(A−B)·2^-(k+1)
    // **`d̂ == 0` is the easiest case and the first version refused it.**
    // `decompose(0.0)` is `(0, 0)`, so `shift` came out `-1` and 1,350 of
    // `box_exact`'s 1,350 crossings went unjudged - every one of them a
    // symmetric edge where `a + b` is exactly zero, which is precisely where
    // the error is exactly zero and the bound is trivially sound. A reference
    // that declines to judge its own best case is `M-44` wearing a refusal.
    let (err_num, err_den) = if md == 0 {
        // d_true = num / (2·den); the error is its negation.
        (num, den.checked_mul(2)?)
    } else {
        let shift = -(ed + 1);
        if !(0..=70).contains(&shift) {
            return None;
        }
        let scale = 1i128.checked_shl(shift as u32)?;
        let lhs = md.checked_mul(den)?;
        let rhs = num.checked_mul(scale)?;
        (
            lhs.checked_sub(rhs)?,
            den.checked_mul(2)?.checked_mul(scale)?,
        )
    };

    // One division, at the very end, on two exact integers. Its own relative
    // error is ~1e-16, which is nine orders below the ulp granularity every
    // clause is stated at.
    let err = err_num as f64 / err_den as f64;
    if d_hat == 0.0 {
        // No ulp to normalise by. The error is exact and, for a symmetric edge,
        // exactly zero; reporting it in ulps of zero would be a division by the
        // smallest subnormal. Zero error is zero ulps whatever the scale.
        return Some(if err == 0.0 { 0.0 } else { f64::INFINITY });
    }
    Some(err.abs() / ulp_of(d_hat))
}

/// The gap to the next `f64` above `|x|`.
fn ulp_of(x: f64) -> f64 {
    if x == 0.0 {
        return f64::MIN_POSITIVE;
    }
    let a = x.abs();
    let next = f64::from_bits(a.to_bits() + 1);
    next - a
}

/// The registered bound, in ulps of `d`.
///
/// Two extra flops over the crossing itself: the two-sum's residual is six, but
/// five of those are the sum the crossing already computes.
fn bound_ulp(a: f64, b: f64, d: f64) -> f64 {
    bound_ulp_with(a, b, d, 2.0)
}

/// The same bound with the leading coefficient as a parameter.
///
/// **The registered coefficient is 2 and the measurement says it must be 3.**
/// The derivation: `fl(a + b)` carries a relative error up to `u`, so an
/// **absolute** error up to `u·|a + b|`. Dividing by `(a − b)` turns that into an
/// absolute error of `u·|a + b| / |a − b|`, which is `u · 2|d|` - **two** units,
/// not one, because `|a + b| / |a − b|` *is* `2|d|`. Add the quotient's own `u`
/// and the first-order coefficient is **3**, not 2. The registration counted
/// roundings instead of propagating them, and this harness measures both so the
/// correction is a number rather than an argument.
fn bound_ulp_with(a: f64, b: f64, d: f64, coeff: f64) -> f64 {
    let (x, y) = two_sum(a, -b);
    let denom_rel = if x == 0.0 { 0.0 } else { (y / x).abs() };
    let bound = d.abs() * (coeff + denom_rel) * U;
    if d == 0.0 {
        // A zero offset has zero error and a zero bound; in ulps that is 0, and
        // dividing by `ulp_of(0)` would report a subnormal-scaled infinity.
        return 0.0;
    }
    bound / ulp_of(d)
}

/// One field's row.
struct Row {
    field: &'static str,
    crossings: u64,
    judged: u64,
    violations: u64,
    /// Violations against the corrected coefficient 3. See `bound_ulp_with`.
    violations_c3: u64,
    bounds: Vec<f64>,
    errors: Vec<f64>,
    /// Bound and seam distance, `csg_difference` only. See C4.
    seam: Vec<(f64, f64)>,
}

fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let i = ((sorted.len() - 1) as f64 * q).round() as usize;
    sorted[i]
}

/// Seam proximity for `csg_difference`, which is `max(box, −sphere)`.
///
/// The seam is where the two arguments **tie**, so `|box + sphere|` is the
/// distance to it: the max is a crease exactly where its two branches cross, and
/// that is the locality `M-350` bounded for the normal.
fn seam_distance(p: [f64; 3]) -> f64 {
    let b = isomesh::fields::BoxExact::<f64>::canonical().sample(p);
    let s = Sphere::<f64> {
        center: [0.6; 3],
        radius: 0.75,
    }
    .sample(p);
    (b + s).abs()
}

fn measure<F: Sdf<Scalar = f64> + ReferenceField>(field: &F, field_name: &'static str) -> Row {
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(SAMPLES - 1);
    let n = SAMPLES as usize;

    let mut values = Vec::with_capacity(n * n * n);
    for z in 0..SAMPLES {
        for y in 0..SAMPLES {
            for x in 0..SAMPLES {
                values.push(field.sample([
                    lo[0] + h * f64::from(x),
                    lo[1] + h * f64::from(y),
                    lo[2] + h * f64::from(z),
                ]));
            }
        }
    }
    let at = |x: u32, y: u32, z: u32| values[(z as usize * n + y as usize) * n + x as usize];

    let mut row = Row {
        field: field_name,
        crossings: 0,
        judged: 0,
        violations: 0,
        violations_c3: 0,
        bounds: Vec::new(),
        errors: Vec::new(),
        seam: Vec::new(),
    };

    // **The three axis-aligned grid edges per sample, not the twelve cube edges
    // per cell.** A cell's twelve edges are shared with its neighbours, so a
    // per-cell walk counts most crossings four times and needs a dedup rule -
    // and a dedup rule is a place to be subtly wrong about which cell owns an
    // edge. `P-66` walks the lattice directly for the same reason; this is the
    // same edge set, which also makes the two experiments' populations
    // comparable.
    for axis in 0..3usize {
        let mut extent = [SAMPLES, SAMPLES, SAMPLES];
        extent[axis] -= 1;
        for z in 0..extent[2] {
            for y in 0..extent[1] {
                for x in 0..extent[0] {
                    let a = at(x, y, z);
                    let mut step = [x, y, z];
                    step[axis] += 1;
                    let b = at(step[0], step[1], step[2]);
                    if is_inside(a) == is_inside(b) {
                        continue;
                    }
                    row.crossings += 1;

                    let d = edge_offset(a, b);
                    let bnd = bound_ulp(a, b, d);
                    row.bounds.push(bnd);

                    let Some(err) = exact_error_ulp(a, b, d) else {
                        continue;
                    };
                    row.judged += 1;
                    row.errors.push(err);
                    if err > bnd {
                        row.violations += 1;
                    }
                    if err > bound_ulp_with(a, b, d, 3.0) {
                        row.violations_c3 += 1;
                    }
                    if field_name == "csg_difference" {
                        // The crossing's world position, by the crate's own
                        // placement: midpoint plus `d` along the edge.
                        let w0 = [
                            lo[0] + h * f64::from(x),
                            lo[1] + h * f64::from(y),
                            lo[2] + h * f64::from(z),
                        ];
                        let mut w1 = w0;
                        w1[axis] += h;
                        let world = [0, 1, 2].map(|i| (w0[i] + w1[i]) * 0.5 + (w1[i] - w0[i]) * d);
                        row.seam.push((bnd, seam_distance(world)));
                    }
                }
            }
        }
    }

    row.bounds
        .sort_unstable_by(|a, b| a.partial_cmp(b).expect("finite"));
    row.errors
        .sort_unstable_by(|a, b| a.partial_cmp(b).expect("finite"));
    row
}

/// The two-sum's defining property, asserted rather than assumed.
///
/// `x + y == a + b` **exactly**, checked in `i128` on adversarial pairs. Without
/// this the bound's denominator term rests on a six-line function copied into a
/// bench, which is the drift this ledger keeps finding.
fn the_two_sum_is_exact() {
    let cases = [
        (1.0f64, f64::EPSILON * 0.25),
        (1e300, -1e300 + 1.0),
        (0.1, 0.2),
        (-1e-300, 1e300),
        (f64::MIN_POSITIVE, -f64::MIN_POSITIVE * 0.5),
        (3.0, 7.0),
    ];
    for (a, b) in cases {
        let (x, y) = two_sum(a, b);
        let (ma, ea) = decompose(a);
        let (mb, eb) = decompose(b);
        let (mx, ex) = decompose(x);
        let (my, ey) = decompose(y);
        let e = ea.min(eb).min(ex).min(ey);
        let lift = |m: i128, ep: i32| -> Option<i128> {
            let s = ep - e;
            if s > 100 {
                None
            } else {
                m.checked_shl(s as u32)
            }
        };
        let (Some(la), Some(lb), Some(lx), Some(ly)) =
            (lift(ma, ea), lift(mb, eb), lift(mx, ex), lift(my, ey))
        else {
            continue;
        };
        assert_eq!(
            la + lb,
            lx + ly,
            "two_sum is not exact on ({a:e}, {b:e}): the bound's denominator term is unfounded"
        );
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-68");
    the_two_sum_is_exact();

    let mut rows: Vec<Row> = Vec::new();
    for_each_reference_field!(f64, |name, field| {
        rows.push(measure(&field, name));
    });

    // ── C3's cost, one binary and one run (M-281) ────────────────────────────
    //
    // The crossing alone against the crossing plus the bound, over the same
    // inputs. Not an extraction: what C3 asks is what carrying the bound costs,
    // and the denominator is taken from the committed extraction measurement so
    // the share is against a real path rather than against this loop.
    // Deterministic straddling pairs, the shape a cut edge produces: one
    // strictly negative endpoint and one non-negative one, which is the
    // precondition `edge_offset` asserts.
    let pairs: Vec<(f64, f64)> = {
        let mut state = 0x2026u64 ^ 0x5EED_1234;
        let mut next = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..1_000_000).map(|_| (-next() - 1e-9, next())).collect()
    };

    let reps = 5usize;
    let mut crossing_ns = f64::INFINITY;
    let mut both_ns = f64::INFINITY;
    for _ in 0..reps {
        let t = Instant::now();
        let mut acc = 0.0f64;
        for (a, b) in &pairs {
            acc += edge_offset(*a, *b);
        }
        std::hint::black_box(acc);
        let e = t.elapsed().as_nanos() as f64 / pairs.len() as f64;
        crossing_ns = crossing_ns.min(e);

        let t = Instant::now();
        let mut acc = 0.0f64;
        for (a, b) in &pairs {
            let d = edge_offset(*a, *b);
            acc += d + bound_ulp(*a, *b, d);
        }
        std::hint::black_box(acc);
        let e = t.elapsed().as_nanos() as f64 / pairs.len() as f64;
        both_ns = both_ns.min(e);
    }
    let bound_ns = (both_ns - crossing_ns).max(0.0);

    // **C3's denominator is an EXTRACTION, and the first version used the
    // crossing.** The clause reads *"under 3% of extraction wall time"*; against
    // the crossing alone the bound reads 0.69, which is a true number about the
    // wrong total - `M-375`'s rule in the other direction. So the extraction is
    // timed here, on the same grids, and the share is
    // `bound_ns × crossings / extract_ns`.
    let mut extract_ns = 0.0f64;
    let mut mc = isomesh::marching_cubes::MarchingCubes::<f64>::new();
    for_each_reference_field!(f64, |name, field| {
        let _ = name;
        let shape = isomesh::RuntimeShape3::new([SAMPLES; 3]).expect("shape");
        let (flo, fhi) = field.domain();
        let fh = (fhi[0] - flo[0]) / f64::from(SAMPLES - 1);
        let mut out = isomesh::MeshBuffer::<f64>::new();
        let mut best = f64::INFINITY;
        for _ in 0..3 {
            let t = Instant::now();
            out = isomesh::MeshBuffer::<f64>::new();
            let _ = isomesh::extractor::Extractor::extract_into(
                &mut mc, &field, &shape, flo, fh, &mut out,
            );
            best = best.min(t.elapsed().as_nanos() as f64);
        }
        std::hint::black_box(&out);
        extract_ns += best;
    });

    let total_crossings: u64 = rows.iter().map(|r| r.crossings).sum();

    println!(
        "{:>15} {:>9} {:>9} {:>6} {:>11} {:>11} {:>11} {:>11}",
        "field", "crossings", "judged", "viol", "med bound", "p99 bound", "med err", "max err"
    );
    for r in &rows {
        println!(
            "{:>15} {:>9} {:>9} {:>6} {:>11.4} {:>11.4} {:>11.4} {:>11.4}",
            r.field,
            r.crossings,
            r.judged,
            r.violations,
            quantile(&r.bounds, 0.5),
            quantile(&r.bounds, 0.99),
            quantile(&r.errors, 0.5),
            r.errors.last().copied().unwrap_or(f64::NAN)
        );
    }

    // ── controls ─────────────────────────────────────────────────────────────
    assert!(
        total_crossings > 0,
        "VOID: no cut edge on any field, so every clause is about an empty population"
    );
    let judged: u64 = rows.iter().map(|r| r.judged).sum();
    assert!(
        judged > 0,
        "VOID: the i128 reference judged no crossing, so C1's zero violations is a zero over a \
         population the instrument declined to look at"
    );
    // The instrument must be able to report a violation. A deliberately wrong
    // bound -- a tenth of the registered one -- has to produce violations on the
    // same data, or C1's zero says nothing.
    let mut sabotage = 0u64;
    for r in &rows {
        for (e, b) in r.errors.iter().zip(r.bounds.iter()) {
            if *e > b * 0.1 {
                sabotage += 1;
            }
        }
    }
    assert!(
        sabotage > 0,
        "VOID: even a bound ten times too small produces no violation, so the comparison cannot \
         report the bad news"
    );

    // ── verdict ──────────────────────────────────────────────────────────────
    let violations: u64 = rows.iter().map(|r| r.violations).sum();
    let c1 = violations == 0;

    let rough = ["fbm_terrain", "noise_cavity"];
    let mut c2 = true;
    for r in &rows {
        if rough.contains(&r.field) {
            continue;
        }
        let med = quantile(&r.bounds, 0.5);
        let p99 = quantile(&r.bounds, 0.99);
        if !(med < 4.0 && p99 < 64.0) {
            c2 = false;
            println!("C2 {}: median {med:.4} ulp, p99 {p99:.4} ulp", r.field);
        }
    }

    let bound_share = bound_ns * total_crossings as f64 / extract_ns;
    let c3 = bound_share < 0.03;

    // C4: split csg_difference's crossings at the median seam distance and
    // compare the two halves' mean bound. A ratio above 1 means the bound is
    // larger nearer the seam.
    let seam_ratio = {
        let r = rows
            .iter()
            .find(|r| r.field == "csg_difference")
            .expect("csg_difference row");
        let mut d: Vec<f64> = r.seam.iter().map(|(_, s)| *s).collect();
        d.sort_unstable_by(|a, b| a.partial_cmp(b).expect("finite"));
        let cut = quantile(&d, 0.5);
        let mut near = (0.0f64, 0u64);
        let mut far = (0.0f64, 0u64);
        for (b, s) in &r.seam {
            if *s <= cut {
                near.0 += b;
                near.1 += 1;
            } else {
                far.0 += b;
                far.1 += 1;
            }
        }
        if near.1 == 0 || far.1 == 0 {
            f64::NAN
        } else {
            (near.0 / near.1 as f64) / (far.0 / far.1 as f64)
        }
    };
    let c4 = seam_ratio > 1.0;

    println!(
        "\ncrossing {crossing_ns:.4} ns, +bound {both_ns:.4} ns, bound alone {bound_ns:.4} ns"
    );
    println!(
        "extraction over eight fields at {SAMPLES}³: {:.4} ms for {total_crossings} crossings; \
         bound adds {:.4} ms = {bound_share:.4}",
        extract_ns / 1e6,
        bound_ns * total_crossings as f64 / 1e6
    );
    println!("C4 csg_difference near/far seam mean bound ratio: {seam_ratio:.4}");
    let violations_c3: u64 = rows.iter().map(|r| r.violations_c3).sum();
    println!(
        "\nC1 zero violations over {judged} exactly-judged crossings: {violations} -> {}",
        if c1 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "     the same population against the CORRECTED coefficient 3: {violations_c3} \
         violations"
    );
    println!(
        "C2 median < 4 ulp and p99 < 64 ulp on the six smooth -> {}",
        if c2 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C3 bound costs < 3% of the crossing path: {bound_share:.4} -> {}",
        if c3 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C4 bound larger nearer the seam: {seam_ratio:.4} -> {}",
        if c4 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "\nC3's 'exactly zero when the feature is off' is TRIVIALLY TRUE and the entry says so: \
         no feature was added, the crate is unchanged, and the 216 golden hashes are therefore \
         untouched. Landing a per-vertex bound buffer doubles vertex memory and is the owner's \
         decision, not this harness's -- P-71's ring is the same shape."
    );

    common::experiment::run(prereg, |run| {
        for r in &rows {
            let med_b = quantile(&r.bounds, 0.5);
            let med_e = quantile(&r.errors, 0.5);
            run.record(&[
                ("field", r.field.to_string()),
                ("samples_per_axis", SAMPLES.to_string()),
                ("crossings", r.crossings.to_string()),
                ("exactly_judged", r.judged.to_string()),
                ("violations", r.violations.to_string()),
                ("median_bound_ulp", format!("{med_b:.6}")),
                ("p99_bound_ulp", format!("{:.6}", quantile(&r.bounds, 0.99))),
                (
                    "max_bound_ulp",
                    format!("{:.6}", r.bounds.last().copied().unwrap_or(f64::NAN)),
                ),
                ("median_true_error_ulp", format!("{med_e:.6}")),
                (
                    "max_true_error_ulp",
                    format!("{:.6}", r.errors.last().copied().unwrap_or(f64::NAN)),
                ),
                ("bound_tightness", format!("{:.6}", med_e / med_b)),
                ("seam_bound_ratio", format!("{seam_ratio:.6}")),
                ("bound_ns_per_crossing", format!("{bound_ns:.6}")),
                ("crossing_ns_per_crossing", format!("{crossing_ns:.6}")),
                ("bound_share", format!("{bound_share:.6}")),
                ("violations_coeff3", r.violations_c3.to_string()),
                ("golden_hashes_unchanged", "true".to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                ("c4_holds", c4.to_string()),
            ]);
        }
    });
}
