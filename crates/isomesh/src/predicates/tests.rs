//! T-024a and T-024b. The oracle is exact `i128` arithmetic over integer-valued
//! inputs,
//! deliberately **not** this module's own expansion code — checking an exact
//! predicate against itself proves only that it is self-consistent.

// Exact equality is the entire subject here. `two_product` returning a value
// that is *approximately* the rounded product would be a defect, not a rounding
// difference, so every comparison in this file is deliberately strict.
#![allow(clippy::float_cmp)]

use super::*;

/// The determinant, computed exactly in `i128`.
///
/// Valid only for integer-valued inputs small enough that the products cannot
/// overflow. Callers here draw from `±2^31` at the widest, where a difference is
/// at most `2^32` and a product at most `2^64` — comfortably inside `i128`, and
/// every coordinate stays below `2^53` so its `f64` image is exact.
fn oracle(a: [i64; 2], b: [i64; 2], c: [i64; 2]) -> i128 {
    let ax = i128::from(a[0] - c[0]);
    let ay = i128::from(a[1] - c[1]);
    let bx = i128::from(b[0] - c[0]);
    let by = i128::from(b[1] - c[1]);
    ax * by - ay * bx
}

fn as_f64(p: [i64; 2]) -> [f64; 2] {
    [p[0] as f64, p[1] as f64]
}

/// A deterministic integer generator. No `rand` dependency, and the sequence is
/// identical on every platform, so a failure reproduces exactly.
struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        // Numerical Recipes' constants; the low bits are poor, so the callers
        // take from the high half.
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    /// Uniform in `-range..=range`.
    fn coord(&mut self, range: i64) -> i64 {
        let span = (range as u64) * 2 + 1;
        ((self.next_u64() >> 16) % span) as i64 - range
    }
}

// ── the primitives ──────────────────────────────────────────────────────────

/// An integer-valued `f64` as an exact `i128`. Sound only because every value
/// this module's tests feed it is an integer below `2^53`, or a sum/product of
/// such integers that a `two_*` routine returned — all of which are integers
/// exactly representable in both types.
fn exact(v: f64) -> i128 {
    assert_eq!(v.fract(), 0.0, "oracle needs an integer-valued float: {v}");
    v as i128
}

/// `two_product` recovers the exact roundoff, which is what fails silently if a
/// compiler ever contracts the multiply into an FMA.
///
/// The module's correctness rests on this, so it is asserted rather than
/// assumed. See the module docs.
#[test]
fn two_product_is_exact() {
    let mut lcg = Lcg(24_001);
    for _ in 0..50_000 {
        // Factors near 2^31, so the product needs ~62 bits and genuinely rounds.
        // Both factors stay exactly representable, which is what lets the i128
        // oracle stand outside this module's own arithmetic.
        let a = lcg.coord(1 << 31);
        let b = lcg.coord(1 << 31);
        let (x, y) = two_product(a as f64, b as f64);
        assert_eq!(
            x,
            (a as f64) * (b as f64),
            "high word is the rounded product"
        );
        assert_eq!(
            exact(x) + exact(y),
            i128::from(a) * i128::from(b),
            "two_product({a}, {b}) must be exact -- an FMA contraction breaks this"
        );
    }
}

#[test]
fn two_sum_and_two_diff_are_exact() {
    let mut lcg = Lcg(24_002);
    for _ in 0..50_000 {
        // One large and one small, so the sum needs more bits than it has.
        let a = lcg.coord(1 << 52) * 16;
        let b = lcg.coord(1 << 20);
        let (x, y) = two_sum(a as f64, b as f64);
        assert_eq!(x, (a as f64) + (b as f64));
        assert_eq!(
            exact(x) + exact(y),
            i128::from(a) + i128::from(b),
            "two_sum({a}, {b}) must be exact"
        );

        let (dx, dy) = two_diff(a as f64, b as f64);
        assert_eq!(dx, (a as f64) - (b as f64));
        assert_eq!(
            exact(dx) + exact(dy),
            i128::from(a) - i128::from(b),
            "two_diff({a}, {b}) must be exact"
        );
    }
}

#[test]
fn two_two_diff_is_exact() {
    let mut lcg = Lcg(24_003);
    for _ in 0..20_000 {
        let a = lcg.coord(1 << 31);
        let b = lcg.coord(1 << 31);
        let c = lcg.coord(1 << 31);
        let d = lcg.coord(1 << 31);
        let (a1, a0) = two_product(a as f64, b as f64);
        let (b1, b0) = two_product(c as f64, d as f64);
        let want = i128::from(a) * i128::from(b) - i128::from(c) * i128::from(d);
        let got: i128 = two_two_diff(a1, a0, b1, b0)
            .iter()
            .copied()
            .map(exact)
            .sum();
        assert_eq!(got, want, "({a}*{b}) - ({c}*{d})");
    }
}

#[test]
fn fast_two_sum_agrees_with_two_sum_when_ordered() {
    let mut lcg = Lcg(24_004);
    for _ in 0..10_000 {
        let a = lcg.coord(1 << 20) as f64 * 1.5;
        let b = lcg.coord(1 << 10) as f64 * 0.25;
        let (a, b) = if a.abs() >= b.abs() { (a, b) } else { (b, a) };
        assert_eq!(fast_two_sum(a, b), two_sum(a, b));
    }
}

#[test]
fn splitter_matches_the_significand_width() {
    // 2^ceil(p/2) + 1 for each type. Wrong values here would not fail loudly --
    // two_product would just stop being exact -- so they are pinned.
    assert_eq!(f64::SPLITTER, 134_217_729.0);
    assert_eq!(f32::SPLITTER, 4_097.0);
    assert_eq!(f64::UNIT_ROUNDOFF, f64::EPSILON / 2.0);
    assert_eq!(f32::UNIT_ROUNDOFF, f32::EPSILON / 2.0);
    // The paper's epsilon, stated explicitly: 2^-53 and 2^-24.
    assert_eq!(f64::UNIT_ROUNDOFF, 2.0_f64.powi(-53));
    assert_eq!(f32::UNIT_ROUNDOFF, 2.0_f32.powi(-24));
}

// ── expansion addition ──────────────────────────────────────────────────────

#[test]
fn expansion_sum_totals_are_preserved() {
    let mut lcg = Lcg(24_005);
    for _ in 0..5_000 {
        // Two expansions built the way the predicate builds them.
        let e = two_two_diff(
            lcg.coord(1 << 20) as f64,
            lcg.coord(1 << 4) as f64 * f64::EPSILON,
            lcg.coord(1 << 20) as f64,
            lcg.coord(1 << 4) as f64 * f64::EPSILON,
        );
        let f = two_two_diff(
            lcg.coord(1 << 20) as f64,
            lcg.coord(1 << 4) as f64 * f64::EPSILON,
            lcg.coord(1 << 20) as f64,
            lcg.coord(1 << 4) as f64 * f64::EPSILON,
        );
        let mut h = [0.0; 8];
        let len = fast_expansion_sum(&e, &f, &mut h).expect("buffer is exactly sized");
        assert!(len <= 8);
        let want: f64 = e.iter().chain(f.iter()).sum();
        let got: f64 = h.iter().take(len).sum();
        assert_eq!(
            got, want,
            "the expansion's total must equal the inputs' total"
        );
        // Zero elimination: no interior zero components.
        assert!(
            h.iter().take(len.saturating_sub(1)).all(|&c| c != 0.0),
            "zeros must be eliminated"
        );
    }
}

#[test]
fn expansion_sum_refuses_a_short_buffer_rather_than_truncating() {
    let e = [1.0_f64, 2.0];
    let f = [3.0_f64, 4.0];
    let mut too_small = [0.0; 3];
    assert_eq!(fast_expansion_sum(&e, &f, &mut too_small), None);
    let mut exactly_sized = [0.0; 4];
    assert!(fast_expansion_sum(&e, &f, &mut exactly_sized).is_some());
}

// ── orient2d ────────────────────────────────────────────────────────────────

#[test]
fn orient2d_matches_the_exact_oracle_on_random_input() {
    let mut lcg = Lcg(24_006);
    for _ in 0..200_000 {
        let a = [lcg.coord(1 << 20), lcg.coord(1 << 20)];
        let b = [lcg.coord(1 << 20), lcg.coord(1 << 20)];
        let c = [lcg.coord(1 << 20), lcg.coord(1 << 20)];
        let want = oracle(a, b, c).signum() as i32;
        let got = orient2d(as_f64(a), as_f64(b), as_f64(c));
        let got_sign = if got > 0.0 {
            1
        } else if got < 0.0 {
            -1
        } else {
            0
        };
        assert_eq!(got_sign, want, "orient2d({a:?}, {b:?}, {c:?}) = {got}");
    }
}

/// The adversarial case: points drawn from a tiny range so that near-degenerate
/// and exactly-degenerate configurations are common rather than vanishingly
/// rare. Random input over a wide range almost never exercises the exact path.
#[test]
fn orient2d_matches_the_oracle_where_degeneracy_is_common() {
    let mut lcg = Lcg(24_007);
    let mut collinear = 0u32;
    for _ in 0..200_000 {
        let a = [lcg.coord(3), lcg.coord(3)];
        let b = [lcg.coord(3), lcg.coord(3)];
        let c = [lcg.coord(3), lcg.coord(3)];
        let want = oracle(a, b, c).signum() as i32;
        if want == 0 {
            collinear += 1;
        }
        let got = orient2d(as_f64(a), as_f64(b), as_f64(c));
        let got_sign = if got > 0.0 {
            1
        } else if got < 0.0 {
            -1
        } else {
            0
        };
        assert_eq!(got_sign, want, "orient2d({a:?}, {b:?}, {c:?}) = {got}");
    }
    assert!(
        collinear > 10_000,
        "fixture should be degenerate often; saw {collinear}"
    );
}

/// T-024a's acceptance fixture: a configuration the naive float determinant gets
/// **wrong**, not merely imprecise, and that the exact path gets right.
///
/// The points are collinear by construction — `c` is on the line through `a` and
/// `b` — but the coordinates are chosen so that the two products in the naive
/// determinant are each rounded, and their difference is a pure rounding
/// artefact with a nonzero sign.
#[test]
fn the_float_path_misclassifies_a_fixture_the_exact_path_gets_right() {
    // Found by exhaustive search over Bezout pairs (see the commit that added
    // this test). `a` and `b` are chosen so that the EXACT determinant is 1 --
    // the smallest nonzero value it can take on integer input -- while the two
    // products land near 2^61, where one ulp is 512. The naive difference is
    // then entirely rounding, and it collapses to exactly 0.0: the float path
    // does not merely lose precision, it reports three points as COLLINEAR that
    // are not. Every coordinate is below 2^53, so all three are exactly
    // representable and the i128 oracle is exact.
    let ai = [2_147_483_647_i64, 2_147_483_645];
    let bi = [-1_073_741_823_i64, -1_073_741_822];
    let ci = [0_i64, 0];

    let a = as_f64(ai);
    let b = as_f64(bi);
    let c = as_f64(ci);

    // Naive, unfiltered: exactly the expression `orient2d`'s stage A computes,
    // but returned without the error-bound check that makes it safe.
    let naive = (a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0]);

    assert_eq!(oracle(ai, bi, ci), 1, "the true determinant is 1");
    assert_eq!(
        naive, 0.0,
        "the fixture exists because the naive form says collinear -- if a \
         compiler change makes it exact, replace the fixture rather than \
         deleting the test"
    );
    assert!(
        orient2d(a, b, c) > 0.0,
        "the exact path must report counterclockwise where the float path said collinear"
    );
}

#[test]
fn orient2d_is_antisymmetric_and_sign_correct() {
    let ccw = ([0.0_f64, 0.0], [1.0_f64, 0.0], [0.0_f64, 1.0]);
    assert!(orient2d(ccw.0, ccw.1, ccw.2) > 0.0, "counterclockwise");
    assert!(orient2d(ccw.0, ccw.2, ccw.1) < 0.0, "clockwise");
    // Swapping any two arguments flips the sign, exactly.
    let d = orient2d(ccw.0, ccw.1, ccw.2);
    assert_eq!(orient2d(ccw.1, ccw.0, ccw.2), -d);
    assert_eq!(orient2d(ccw.0, ccw.2, ccw.1), -d);
}

#[test]
fn orient2d_is_exactly_zero_on_collinear_and_coincident_input() {
    // Collinear along each axis and the diagonal.
    assert_eq!(orient2d([0.0_f64, 0.0], [1.0, 0.0], [2.0, 0.0]), 0.0);
    assert_eq!(orient2d([0.0_f64, 0.0], [0.0, 1.0], [0.0, 2.0]), 0.0);
    assert_eq!(orient2d([0.0_f64, 0.0], [1.0, 1.0], [2.0, 2.0]), 0.0);
    // Two coincident points make the triangle degenerate whatever the third is.
    assert_eq!(orient2d([1.5_f64, -2.25], [1.5, -2.25], [9.0, 4.0]), 0.0);
    assert_eq!(orient2d([9.0_f64, 4.0], [1.5, -2.25], [1.5, -2.25]), 0.0);
}

#[test]
fn orient2d_works_in_f32_too() {
    let mut lcg = Lcg(24_008);
    for _ in 0..50_000 {
        // f32 has a 24-bit significand, so integers to 2^11 keep the oracle's
        // products exactly representable as f32 inputs.
        let a = [lcg.coord(1 << 10), lcg.coord(1 << 10)];
        let b = [lcg.coord(1 << 10), lcg.coord(1 << 10)];
        let c = [lcg.coord(1 << 10), lcg.coord(1 << 10)];
        let want = oracle(a, b, c).signum() as i32;
        let got = orient2d(
            [a[0] as f32, a[1] as f32],
            [b[0] as f32, b[1] as f32],
            [c[0] as f32, c[1] as f32],
        );
        let got_sign = if got > 0.0 {
            1
        } else if got < 0.0 {
            -1
        } else {
            0
        };
        assert_eq!(got_sign, want, "f32 orient2d({a:?}, {b:?}, {c:?}) = {got}");
    }
}

#[test]
fn orient2d_is_deterministic() {
    let a = [0.1_f64, 0.2];
    let b = [0.3, 0.4];
    let c = [0.5, 0.600_000_000_000_000_1];
    let first = orient2d(a, b, c);
    for _ in 0..1_000 {
        assert_eq!(orient2d(a, b, c), first);
    }
}

// ── incircle ────────────────────────────────────────────────────────────────

/// The lifted 4x4 determinant, computed exactly in `i128`.
///
/// `z = x^2 + y^2`, so at coordinates up to `2^20` a lift is `2^41`, a 2x2 cross
/// product `2^41`, and their product `2^82` -- well inside `i128`.
fn incircle_oracle(a: [i64; 2], b: [i64; 2], c: [i64; 2], d: [i64; 2]) -> i128 {
    let cross = |u: [i64; 2], v: [i64; 2]| -> i128 {
        i128::from(u[0]) * i128::from(v[1]) - i128::from(u[1]) * i128::from(v[0])
    };
    let lift = |p: [i64; 2]| -> i128 {
        i128::from(p[0]) * i128::from(p[0]) + i128::from(p[1]) * i128::from(p[1])
    };
    let (ab, ac, ad) = (cross(a, b), cross(a, c), cross(a, d));
    let (bc, bd, cd) = (cross(b, c), cross(b, d), cross(c, d));
    lift(a) * (bc - bd + cd) - lift(b) * (ac - ad + cd) + lift(c) * (ab - ad + bd)
        - lift(d) * (ab - ac + bc)
}

fn sign_of(v: f64) -> i32 {
    if v > 0.0 {
        1
    } else if v < 0.0 {
        -1
    } else {
        0
    }
}

#[test]
fn incircle_matches_a_hand_computed_circle() {
    // The unit circle through (0,0), (1,0), (0,1) is centred (0.5, 0.5). These
    // three are counterclockwise -- orient2d says so -- which fixes the sign.
    let a = [0.0_f64, 0.0];
    let b = [1.0_f64, 0.0];
    let c = [0.0_f64, 1.0];
    assert!(orient2d(a, b, c) > 0.0, "fixture must be counterclockwise");

    // The centre is inside.
    assert!(incircle(a, b, c, [0.5, 0.5]) > 0.0);
    // (1, 1) is exactly on the circle.
    assert_eq!(incircle(a, b, c, [1.0, 1.0]), 0.0);
    // Far away is outside.
    assert!(incircle(a, b, c, [5.0, 5.0]) < 0.0);
    // Each defining point is on its own circle, exactly.
    assert_eq!(incircle(a, b, c, a), 0.0);
    assert_eq!(incircle(a, b, c, b), 0.0);
    assert_eq!(incircle(a, b, c, c), 0.0);
}

#[test]
fn incircle_inverts_with_the_winding_of_abc() {
    let a = [0.0_f64, 0.0];
    let b = [1.0_f64, 0.0];
    let c = [0.0_f64, 1.0];
    let inside = [0.5_f64, 0.5];
    // Swapping two of the first three reverses the winding, hence the sign.
    assert_eq!(
        sign_of(incircle(a, b, c, inside)),
        -sign_of(incircle(b, a, c, inside)),
        "the documented sign convention is stated for counterclockwise abc"
    );
}

#[test]
fn incircle_matches_the_exact_oracle_on_random_input() {
    let mut lcg = Lcg(24_009);
    for _ in 0..100_000 {
        let a = [lcg.coord(1 << 20), lcg.coord(1 << 20)];
        let b = [lcg.coord(1 << 20), lcg.coord(1 << 20)];
        let c = [lcg.coord(1 << 20), lcg.coord(1 << 20)];
        let d = [lcg.coord(1 << 20), lcg.coord(1 << 20)];
        let want = incircle_oracle(a, b, c, d).signum() as i32;
        let got = incircle(as_f64(a), as_f64(b), as_f64(c), as_f64(d));
        assert_eq!(
            sign_of(got),
            want,
            "incircle({a:?}, {b:?}, {c:?}, {d:?}) = {got}"
        );
    }
}

/// Small coordinates so cocircular and near-cocircular configurations are common
/// and the exact path actually runs.
#[test]
fn incircle_matches_the_oracle_where_degeneracy_is_common() {
    let mut lcg = Lcg(24_010);
    let mut cocircular = 0u32;
    for _ in 0..100_000 {
        let a = [lcg.coord(3), lcg.coord(3)];
        let b = [lcg.coord(3), lcg.coord(3)];
        let c = [lcg.coord(3), lcg.coord(3)];
        let d = [lcg.coord(3), lcg.coord(3)];
        let want = incircle_oracle(a, b, c, d).signum() as i32;
        if want == 0 {
            cocircular += 1;
        }
        let got = incircle(as_f64(a), as_f64(b), as_f64(c), as_f64(d));
        assert_eq!(
            sign_of(got),
            want,
            "incircle({a:?}, {b:?}, {c:?}, {d:?}) = {got}"
        );
    }
    assert!(
        cocircular > 5_000,
        "fixture should be degenerate often; saw {cocircular}"
    );
}

/// Large coordinates, where the lifted terms reach `2^80` and the filtered
/// estimate has no chance -- this is the case that exercises `incircle_exact`
/// and its 384-component accumulator hardest.
#[test]
fn incircle_matches_the_oracle_at_large_coordinates() {
    let mut lcg = Lcg(24_011);
    for _ in 0..50_000 {
        let a = [lcg.coord(1 << 25), lcg.coord(1 << 25)];
        let b = [lcg.coord(1 << 25), lcg.coord(1 << 25)];
        let c = [lcg.coord(1 << 25), lcg.coord(1 << 25)];
        // Put d very close to a, so the determinant is small against the lifts.
        let d = [a[0] + lcg.coord(2), a[1] + lcg.coord(2)];
        let want = incircle_oracle(a, b, c, d).signum() as i32;
        let got = incircle(as_f64(a), as_f64(b), as_f64(c), as_f64(d));
        assert_eq!(
            sign_of(got),
            want,
            "incircle({a:?}, {b:?}, {c:?}, {d:?}) = {got}"
        );
    }
}

#[test]
fn incircle_works_in_f32_too() {
    let mut lcg = Lcg(24_012);
    for _ in 0..50_000 {
        let a = [lcg.coord(1 << 8), lcg.coord(1 << 8)];
        let b = [lcg.coord(1 << 8), lcg.coord(1 << 8)];
        let c = [lcg.coord(1 << 8), lcg.coord(1 << 8)];
        let d = [lcg.coord(1 << 8), lcg.coord(1 << 8)];
        let want = incircle_oracle(a, b, c, d).signum() as i32;
        let f = |p: [i64; 2]| [p[0] as f32, p[1] as f32];
        let got = incircle(f(a), f(b), f(c), f(d));
        assert_eq!(
            sign_of(f64::from(got)),
            want,
            "f32 incircle({a:?}, {b:?}, {c:?}, {d:?}) = {got}"
        );
    }
}

#[test]
fn incircle_is_deterministic() {
    let a = [0.1_f64, 0.2];
    let b = [0.3, 0.4];
    let c = [0.5, 0.600_000_000_000_000_1];
    let d = [0.7, 0.8];
    let first = incircle(a, b, c, d);
    for _ in 0..1_000 {
        assert_eq!(incircle(a, b, c, d), first);
    }
}

#[test]
fn scale_expansion_is_exact() {
    let mut lcg = Lcg(24_013);
    for _ in 0..20_000 {
        let p = lcg.coord(1 << 31);
        let q = lcg.coord(1 << 31);
        let r = lcg.coord(1 << 31);
        let s = lcg.coord(1 << 31);
        let (hi, lo) = two_product(p as f64, q as f64);
        let (hi2, lo2) = two_product(r as f64, s as f64);
        let e = two_two_diff(hi, lo, hi2, lo2);
        let scale = lcg.coord(1 << 20);
        let mut out = [0.0; 8];
        let len = scale_expansion(&e, scale as f64, &mut out).expect("sized 2n");
        let want =
            (i128::from(p) * i128::from(q) - i128::from(r) * i128::from(s)) * i128::from(scale);
        let got: i128 = out.iter().take(len).copied().map(exact).sum();
        assert_eq!(got, want, "({p}*{q} - {r}*{s}) * {scale}");
    }
}

/// The buffer chain in `incircle_exact` is sized at its exact worst case, and
/// every `None` arm there is unreachable only because of it. If a buffer were
/// ever too small the function would return `ZERO`, which reads as "exactly
/// cocircular" -- a wrong answer that looks like a legitimate degeneracy rather
/// than a failure. So the lengths are pinned here directly, on the helpers.
#[test]
fn expansion_buffers_are_exactly_the_worst_case() {
    let mut lcg = Lcg(24_014);
    let mut worst_minor = 0usize;
    let mut worst_term = 0usize;
    for _ in 0..20_000 {
        // Full-width coordinates, which is what makes the expansions longest.
        let p = |lcg: &mut Lcg| [lcg.coord(1 << 26) as f64, lcg.coord(1 << 26) as f64];
        let (a, b, c, d) = (p(&mut lcg), p(&mut lcg), p(&mut lcg), p(&mut lcg));

        let bc = cross_exact(b, c);
        let bd = cross_exact(b, d);
        let cd = cross_exact(c, d);
        assert_eq!(bc.len(), 4, "cross_exact is always four components");

        let mut minor = [0.0; 12];
        let minor_len = sum3(&bc, &negated(bd), &cd, &mut minor).expect("12 is the bound");
        worst_minor = worst_minor.max(minor_len);

        let mut term = [0.0; 96];
        let term_len = lift(a, &minor[..minor_len], &mut term).expect("96 is the bound");
        worst_term = worst_term.max(term_len);
    }
    // Observed lengths must sit inside the declared buffers. These are ceilings,
    // not equalities -- zero elimination usually keeps them well under.
    assert!(worst_minor <= 12, "minor reached {worst_minor}");
    assert!(worst_term <= 96, "term reached {worst_term}");
    // The accumulator's bound follows arithmetically: four terms of at most 96
    // is 384, which is what `incircle_exact` declares. Nothing to assert -- the
    // load-bearing checks are the two above, plus the 300,000 oracle-matched
    // cases, any of which would report a spurious zero if a buffer overflowed.
}
