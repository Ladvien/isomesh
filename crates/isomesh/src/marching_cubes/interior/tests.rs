//! A-002c's tests.
//!
//! The one the ticket exists for is
//! `the_numerator_alone_reads_a_joined_cell_as_separated`: it is Custodio's
//! Figure 6 phenomenon, on a configuration derived here rather than transcribed.

use super::{Interior, SweptFaces, chernyaev_numerator_test};
use crate::Error;
use crate::marching_cubes::ambiguity::face_is_joined;

/// The sign structure that admits a pole, and the one Appendix A's counterexample
/// has: the `A`/`C` diagonal positive on the low face and negative on the high
/// one, with `B`/`D` the other way round.
///
/// Both faces are ambiguous — one diagonal strictly negative and the other
/// strictly positive — which is the precondition the interior test is stated
/// under, so this is a configuration MC33 can actually be asked about.
fn opposed() -> SweptFaces<f64> {
    SweptFaces::new([0.1, -2.0, 10.0, -2.0], [-10.0, 2.0, -0.1, 2.0]).expect("ambiguous faces")
}

#[test]
fn the_derived_coefficients_reproduce_the_numerator() {
    // The module's docs derive `F(t) = a t^2 + b t + c` and note that the three
    // coefficients agree term for term with the ones Custodio prints. That
    // agreement is only worth anything if the derivation is also what the code
    // computes, so this checks the expansion against the direct evaluation --
    // `numerator` interpolates the corners and multiplies, `numerator_roots`
    // uses the expanded form, and they must be the same polynomial.
    let faces = opposed();
    let d = |k: usize| faces.hi[k] - faces.lo[k];
    let a = d(0) * d(2) - d(1) * d(3);
    let b = faces.lo[2] * d(0) + faces.lo[0] * d(2) - faces.lo[3] * d(1) - faces.lo[1] * d(3);
    let c = faces.lo[0] * faces.lo[2] - faces.lo[1] * faces.lo[3];

    for step in 0..=20 {
        let t = f64::from(step) / 20.0;
        let expanded = a * t * t + b * t + c;
        let direct = faces.numerator(t);
        assert!(
            (expanded - direct).abs() <= 1e-12 * direct.abs().max(1.0),
            "t = {t}: expanded {expanded}, direct {direct}"
        );
    }

    // And the roots the sign walk uses really are roots.
    for root in faces.numerator_roots() {
        assert!(
            faces.numerator(root).abs() <= 1e-12 * a.abs().max(1.0),
            "root {root} is not a root: F = {}",
            faces.numerator(root)
        );
    }
}

#[test]
fn the_numerator_alone_reads_a_joined_cell_as_separated() {
    // **A-002c's acceptance, and Custodio §5.1's phenomenon.**
    //
    // `F` is convex here -- a = (A1-A0)(C1-C0) - (B1-B0)(D1-D0) = 86.01 > 0 --
    // and negative at both ends, and the maximum of a convex function on a
    // closed interval is at an endpoint. So `F < 0` on the whole sweep and
    // Chernyaev's numerator test finds no positive saddle anywhere.
    //
    // The saddle value is `F / D`, and `D` runs from +14.1 to -14.1, so beyond
    // the pole at t = 0.5 a negative numerator over a negative denominator is a
    // **positive saddle**. The cell is joined and the quadratic cannot say so:
    // it has two sign changes to spend and the quotient needs three.
    let faces = opposed();

    // The denominator runs from positive to negative, which is what puts a pole
    // inside the sweep. The signs are the claim; the magnitudes are incidental
    // and asserting them exactly would be asserting a rounding.
    assert!(faces.denominator(0.0) > 0.0);
    assert!(faces.denominator(1.0) < 0.0);
    assert_eq!(faces.pole(), Some(0.5));

    // The numerator never turns positive.
    for step in 0..=100 {
        let t = f64::from(step) / 100.0;
        assert!(faces.numerator(t) < 0.0, "F({t}) = {}", faces.numerator(t));
    }
    // But the saddle does, past the pole.
    assert!(
        faces.saddle(0.75) > 0.0,
        "saddle(0.75) = {}",
        faces.saddle(0.75)
    );

    assert_eq!(chernyaev_numerator_test(&faces), Interior::Separated);
    assert_eq!(faces.test(), Interior::Joined);
}

#[test]
fn the_saddle_position_is_where_the_bilinear_gradient_vanishes() {
    // A saddle is a critical point, so the bilinear function's two partial
    // derivatives must both be zero there -- checked rather than trusted,
    // because Equation (1) is transcribed geometry and a swapped coordinate
    // would still look plausible.
    //
    // On the unit square with A at (0,0), B at (1,0), C at (1,1), D at (0,1):
    //   f(u,v) = A(1-u)(1-v) + Bu(1-v) + Cuv + D(1-u)v
    let faces = opposed();
    for step in 0..=20 {
        let t = f64::from(step) / 20.0;
        if faces.pole() == Some(t) {
            continue;
        }
        let s = 1.0 - t;
        let c = |k: usize| faces.lo[k] * s + faces.hi[k] * t;
        let (a, b, cc, d) = (c(0), c(1), c(2), c(3));
        let [u, v] = faces.saddle_position(t);

        let df_du = -a * (1.0 - v) + b * (1.0 - v) + cc * v - d * v;
        let df_dv = -a * (1.0 - u) - b * u + cc * u + d * (1.0 - u);
        let scale = (a.abs() + b.abs() + cc.abs() + d.abs()).max(1.0);
        assert!(
            df_du.abs() <= 1e-9 * scale && df_dv.abs() <= 1e-9 * scale,
            "t = {t}: gradient at the saddle is ({df_du}, {df_dv}), not zero"
        );

        // And the value there is the one `saddle` reports.
        let f = a * (1.0 - u) * (1.0 - v) + b * u * (1.0 - v) + cc * u * v + d * (1.0 - u) * v;
        assert!(
            (f - faces.saddle(t)).abs() <= 1e-9 * scale,
            "t = {t}: bilinear at the saddle is {f}, saddle() says {}",
            faces.saddle(t)
        );
    }
}

#[test]
fn without_a_pole_the_correction_changes_nothing() {
    // Where the denominator keeps one sign the quotient's sign is the
    // numerator's, up to that constant sign -- so the correction can only
    // matter when the sweep crosses a pole. Both faces ambiguous the *same* way
    // round is exactly that case: the denominator is a sum of a positive
    // diagonal and a negated negative one at both ends, so it never vanishes.
    let aligned = [
        ([1.0, -1.0, 1.0, -1.0], [2.0, -3.0, 0.5, -1.0]),
        ([0.25, -4.0, 4.0, -0.25], [3.0, -0.5, 0.5, -3.0]),
        ([5.0, -1.0, 0.2, -2.0], [0.1, -0.1, 9.0, -9.0]),
    ];
    for (lo, hi) in aligned {
        let faces = SweptFaces::new(lo, hi).expect("ambiguous faces");
        assert_eq!(faces.pole(), None, "{lo:?} -> {hi:?}");
        assert!(faces.denominator(0.0) > 0.0 && faces.denominator(1.0) > 0.0);
        assert_eq!(
            faces.test(),
            chernyaev_numerator_test(&faces),
            "{lo:?} -> {hi:?}: the two tests should agree with no pole"
        );
    }
}

#[test]
fn a_face_with_no_saddle_is_rejected_rather_than_defaulted() {
    // A + C - B - D == 0 means the bilinear function has no saddle, so the
    // criterion has nothing to evaluate. It cannot arise on an ambiguous face,
    // which is why reaching it is reported as the caller's error rather than
    // absorbed into a topology.
    assert_eq!(
        SweptFaces::new([1.0, 1.0, 1.0, 1.0], [1.0, -1.0, 1.0, -1.0]),
        Err(Error::DegenerateSweep)
    );
    assert_eq!(
        SweptFaces::new([1.0, -1.0, 1.0, -1.0], [2.0, 2.0, 2.0, 2.0]),
        Err(Error::DegenerateSweep)
    );
    // The pole itself is not an error: it is interior to the sweep, not on a
    // face, and it is the whole subject of this module.
    assert!(SweptFaces::new([0.1, -2.0, 10.0, -2.0], [-10.0, 2.0, -0.1, 2.0]).is_ok());
}

#[test]
fn how_often_the_correction_changes_the_answer() {
    // **How much this matters, measured rather than asserted.** Custodio reports
    // the misread occurring once in 10,000 random 5x5x5 fields and six times
    // across 50 isosurfaces of the Skull dataset, which is rare enough that the
    // catalog's advice to skip the interior test for a game is defensible. That
    // is a rate over *fields*; this is the rate over the configurations that can
    // exhibit it at all, which is a different and much larger number, and the
    // two should not be confused.
    //
    // The sweep is over face pairs opposed in sign -- the structure Appendix A's
    // counterexample has and the only one that puts a pole inside the sweep.
    let magnitudes = [0.1, 0.25, 1.0, 4.0, 10.0];
    let mut total = 0u32;
    let mut disagree = 0u32;
    let mut with_pole = 0u32;
    for a0 in magnitudes {
        for c0 in magnitudes {
            for b0 in magnitudes {
                for a1 in magnitudes {
                    for c1 in magnitudes {
                        for b1 in magnitudes {
                            // Low face: A, C positive, B, D negative.
                            // High face: the same diagonals, signs reversed.
                            let lo = [a0, -b0, c0, -b0];
                            let hi = [-a1, b1, -c1, b1];
                            let Ok(faces) = SweptFaces::new(lo, hi) else {
                                continue;
                            };
                            total += 1;
                            if faces.pole().is_some() {
                                with_pole += 1;
                            }
                            if faces.test() != chernyaev_numerator_test(&faces) {
                                disagree += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    std::println!(
        "opposed face pairs: {total} configurations, {with_pole} with a pole in the sweep, \
         {disagree} where the numerator-only test gives the wrong answer"
    );
    // Pinned, so the instrument cannot quietly stop finding anything.
    //
    // **12.6%.** Every opposed pair has a pole, by construction -- the
    // denominator is a sum of positives at one end and of negatives at the
    // other -- and in one in eight of them the pole falls where it changes the
    // answer. So the correction is not a rounding detail on this family; it is
    // rare in *fields* because the family itself is rare, which is a different
    // claim and the one worth quoting.
    assert_eq!((total, with_pole, disagree), (15625, 15625, 1966));
}

#[test]
fn the_test_is_deterministic_and_independent_of_diagonal_order() {
    // The two diagonals are `{A, C}` and `{B, D}`, and rotating the face by one
    // corner swaps them -- which negates the denominator and leaves the
    // numerator alone, so it negates the saddle. That is not a symmetry of the
    // *answer*: it is the question "are the A/C regions joined" becoming "are
    // the B/D regions joined". What must hold is that rotating by *two* corners,
    // which preserves both diagonals, changes nothing at all.
    let faces = opposed();
    let rotate2 = |v: [f64; 4]| [v[2], v[3], v[0], v[1]];
    let rotated =
        SweptFaces::new(rotate2(faces.lo), rotate2(faces.hi)).expect("still ambiguous faces");

    assert_eq!(faces.test(), rotated.test());
    assert_eq!(faces.pole(), rotated.pole());

    // **The saddle values are asserted equal to a tolerance rather than
    // bit-for-bit, and the difference is A-002b's problem rather than a
    // weaker test.** `ambiguity` can promise bit-identical agreement between
    // two cells because its decider reduces to comparing two *products*, and
    // IEEE multiplication is commutative. This module's denominator is
    // `((A + C) - B) - D`, a fixed subtraction order that a rotation permutes,
    // and floating-point addition is not associative. So two cells meeting on a
    // face could evaluate the same sweep to different bits and disagree about a
    // tunnel -- which is M-32's caveat in a new place, and has to be settled
    // before this is wired into an extractor.
    //
    // The pole is skipped rather than tolerated: the saddle is genuinely
    // undefined there, so it comes back infinite, and `inf - inf` is a NaN that
    // no tolerance admits. That is the module's contract appearing in a test,
    // not an inconvenience worked around.
    let pole = faces.pole();
    for step in 0..=20 {
        let t = f64::from(step) / 20.0;
        if pole == Some(t) {
            assert!(
                !faces.saddle(t).is_finite(),
                "the pole should not be finite"
            );
            continue;
        }
        let (a, b) = (faces.saddle(t), rotated.saddle(t));
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "t = {t}: {a} vs {b} -- two cells reading one face must agree \
             bit-for-bit, or they can disagree about a tunnel"
        );
    }
    // Twice, identical: nothing here reads a hash map or an address.
    assert_eq!(faces.test(), faces.test());
}

#[test]
fn catastrophic_cancellation_does_not_hide_a_positive_saddle() {
    // **The stable quadratic, witnessed.** `a = αγ − βδ` is a difference of
    // near-equal products and nothing keeps it away from zero: this fixture
    // has `a ≈ 1e-11` against `b ≈ 2.79`, and its only positive-saddle
    // interval is the 4.8e-6-wide gap between `F`'s small root (0.3240821…)
    // and the pole (0.3240869…). The textbook `(−b + √disc)/2a` computes
    // `√disc − b` — two O(1) values agreeing to eleven digits — and lands the
    // root 8.0e-6 low; the walk then brackets an interval whose midpoint sits
    // below the true root, where the saddle is still negative, and answers
    // Separated. Kahan's `q = −(b + signum(b)·√disc)/2`, roots `q/a` and
    // `c/q`, places the small root to full precision (|F| there: 1e-17
    // against the textbook root's 2e-5) and the walk finds the interval.
    //
    // Both faces are ambiguous — one diagonal strictly negative, the other
    // strictly positive — so this is a configuration MC33 can actually ask
    // about, not a synthetic corner.
    let faces = SweptFaces::new(
        [0.3, -1.0, 0.3, -0.9926951163344483],
        [-1.7, 1.0, -1.7, 1.0073048836605518],
    )
    .expect("ambiguous faces");

    // The positive interval is real whatever the root solver does: just below
    // the pole the saddle is far above zero, while both endpoints read
    // negative — so only the interval walk can find it, and it must.
    let pole = faces.pole().expect("Δ changes sign inside this sweep");
    assert!(
        faces.saddle(pole - 1e-7) > 0.0,
        "saddle(pole − 1e-7) = {}",
        faces.saddle(pole - 1e-7)
    );
    assert!(faces.saddle(0.0) < 0.0 && faces.saddle(1.0) < 0.0);

    assert_eq!(faces.test(), Interior::Joined);
}

#[test]
#[allow(clippy::float_cmp)] // The exact-zero comparison is the guard's point.
fn joined_names_the_positive_regions_and_the_face_rule_names_the_negative() {
    // **The polarity relation between the two deciders, pinned.** Custodio
    // states the interior criterion for the *positive* vertices and this
    // module keeps that sign: positive saddle → `Interior::Joined`, and in
    // this crate positive is *outside*. `ambiguity::face_is_joined` speaks the
    // other way round — its `true` joins the *inside* (negative) corners, and
    // its `d_in > d_out` comparison is algebraically `S < 0`. So on a shared
    // face the two agree about the geometry and use "joined" for opposite
    // sign classes: `face_is_joined(face)` must equal `saddle < 0` there,
    // never `saddle > 0`. A-002b has to translate between them, not equate
    // them; this is the translation, as an assertion.
    let faces = opposed();
    for (v, t) in [(faces.lo, 0.0), (faces.hi, 1.0)] {
        let s = faces.saddle(t);
        assert_ne!(s, 0.0, "a tied saddle at {v:?} would pin nothing");
        assert_eq!(
            face_is_joined(v),
            s < 0.0,
            "face {v:?}: saddle {s} disagrees with the face rule's polarity"
        );
    }

    // The endpoint shortcut in `test` fires on this very fixture: the high
    // face's saddle is positive, so the sweep answers Joined — while the face
    // rule reads the same four values as *not* joining its inside corners.
    // Same numbers, opposite words; that is the trap this test exists for.
    assert!(faces.saddle(1.0) > 0.0);
    assert!(!face_is_joined(faces.hi));
    assert_eq!(faces.test(), Interior::Joined);
}

/// **M-166, closed.** Two cells reading one shared face agree bit-for-bit.
///
/// The face decider can promise this for free because it compares two
/// *products* and IEEE multiplication is commutative. This module's denominator
/// could not: `((A + C) − B) − D` is a fixed subtraction order that a rotation
/// permutes, and floating-point addition is not associative, so two cells
/// meeting on a face could disagree about a tunnel.
///
/// # The fixture is searched for, not chosen
///
/// `the_test_is_deterministic_and_independent_of_diagonal_order` asserts the
/// same property and **passes on its own fixture even with the defect present**,
/// because `opposed()`'s corners are all of similar magnitude and the two
/// orders round identically. That is the fixture trap this project has now hit
/// five times (M-32, M-38, M-44, G-003, and here), so the fixture below comes
/// from a search over corner magnitudes: `(1, 1, 1, 10⁻⁸)` is the first
/// quadruple where the orders disagree, `0.99999999` against
/// `0.9999999900000001`. Over 20,166 quadruples the old order disagrees on
/// **2,764** and the grouped one on **none**.
///
/// The first assertion is that the fixture still distinguishes the two
/// anchorings — without it this could pass by having stopped exercising the
/// case, which is exactly how the older test passes.
#[test]
fn two_cells_reading_one_face_agree_bit_for_bit() {
    // The two orders, spelled out, so the fixture can be shown to separate them.
    let old = |v: [f64; 4]| ((v[0] + v[2]) - v[1]) - v[3];
    let rotate2 = |v: [f64; 4]| [v[2], v[3], v[0], v[1]];
    let corners = [1.0, 1.0, 1.0, 1e-8];
    assert_ne!(
        old(corners).to_bits(),
        old(rotate2(corners)).to_bits(),
        "the fixture no longer distinguishes the two subtraction orders, so it \
         cannot show that grouping the diagonals is what fixes them"
    );

    let faces = SweptFaces::new(corners, corners).expect("a non-degenerate sweep");
    let rotated =
        SweptFaces::new(rotate2(corners), rotate2(corners)).expect("a non-degenerate sweep");

    // A two-corner rotation preserves both diagonals, so nothing may move.
    for step in 0..=20 {
        let t = f64::from(step) / 20.0;
        assert_eq!(
            faces.denominator(t).to_bits(),
            rotated.denominator(t).to_bits(),
            "t = {t}: the denominator moved under a rotation that preserves \
             both diagonals"
        );
        if faces.pole() == Some(t) {
            continue;
        }
        assert_eq!(
            faces.saddle(t).to_bits(),
            rotated.saddle(t).to_bits(),
            "t = {t}: the saddle moved, so two cells could disagree about a tunnel"
        );
    }
    assert_eq!(faces.test(), rotated.test());
    assert_eq!(faces.pole(), rotated.pole());
}
