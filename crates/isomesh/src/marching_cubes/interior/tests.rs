//! A-002c's tests.
//!
//! The one the ticket exists for is
//! `the_numerator_alone_reads_a_joined_cell_as_separated`: it is Custodio's
//! Figure 6 phenomenon, on a configuration derived here rather than transcribed.

use super::{Interior, SweptFaces, chernyaev_numerator_test};
use crate::Error;

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
        assert!(
            (a - b).abs() <= 1e-12 * a.abs().max(1.0),
            "t = {t}: {a} vs {b}"
        );
    }
    // Twice, identical: nothing here reads a hash map or an address.
    assert_eq!(faces.test(), faces.test());
}
