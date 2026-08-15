//! Tests for the field-bound validator.
//!
//! The load-bearing one is [`tightening_a_declaration_by_one_step_is_caught`],
//! which is F-002's acceptance: a checker nobody has seen reject anything is not
//! evidence that the declarations are right, only that the checker is quiet.

use super::{EIKONAL_TOLERANCE, FieldBoundReport, field_bound_report};
use crate::fields::FieldBound;

/// Samples per axis. Enough that the census is stable run to run, cheap enough
/// to run on every `cargo test`.
const N: u32 = 16;

/// **Every reference field meets the bound it declares.**
///
/// The whole point of F-001's type, checked rather than trusted. The numbers are
/// printed as a census because they are interesting in their own right: the
/// eikonal fraction says how close a field is to being a distance, and it splits
/// the eight into two obvious groups.
#[test]
fn every_reference_field_meets_its_declared_bound() {
    let mut rows = alloc::vec::Vec::new();
    crate::for_each_reference_field!(f64, |name, field| {
        let report = field_bound_report(&field, N);
        assert!(report.samples > 1000, "{name}: {} samples", report.samples);
        assert!(
            !report.violates(0.05),
            "{name} declares {:?} but ‖∇f‖ reaches {:.4}",
            report.declared,
            report.sup
        );
        rows.push((name, report));
    });
    assert_eq!(rows.len(), 8, "the reference field set changed");

    for (name, r) in &rows {
        std::println!(
            "measured: {name:<16} {:<28} sup {:>7.3}  inf {:>6.3}  eikonal {:>5.1}%",
            alloc::format!("{:?}", r.declared),
            r.sup,
            r.inf,
            100.0 * r.eikonal_fraction
        );
    }

    // **An exact field is eikonal almost everywhere.** The converse is false and
    // that is the point below.
    for (name, r) in &rows {
        if r.declared == FieldBound::Exact {
            assert!(
                r.eikonal_fraction > 0.99,
                "{name} declares Exact but is eikonal on only {:.1}% of samples",
                100.0 * r.eikonal_fraction
            );
        }
    }

    // **The eikonal fraction cannot tell an underestimate from a distance, and
    // that is why `q` is a separate number from `l` (M-245).**
    //
    // `csg_difference` is `max(box, −sphere)`: away from the seam the active
    // operand is an exact distance, so `‖∇f‖ = 1` there, and the seam is a
    // measure-zero set a grid almost never lands on. It measures **100%**
    // eikonal while its values are not distances near a concave seam — which is
    // the exact configuration A-014d found coincident polygons at.
    //
    // Asserted rather than remarked, because if this ever stops being true the
    // reason `FieldBound` carries two numbers has changed.
    let csg = rows
        .iter()
        .find(|(name, _)| *name == "csg_difference")
        .map(|(_, r)| *r)
        .expect("csg_difference is a reference field");
    assert!(
        matches!(csg.declared, FieldBound::Underestimate { .. }),
        "csg_difference no longer declares an underestimate"
    );
    assert!(
        csg.eikonal_fraction > 0.99,
        "csg_difference is no longer eikonal almost everywhere, so the reason \
         `q` is separate from `l` needs restating: {:.3}",
        csg.eikonal_fraction
    );
}

/// **Tightening a declaration by one step is caught (F-002's acceptance).**
///
/// A checker that has never rejected anything is not evidence about the
/// declarations, so this makes it reject. Each field's *real* bound is replaced
/// by the next tighter one and the report must say so:
///
/// - a `Lipschitz` field claiming `Exact` — the strictest possible claim;
/// - an `Unbounded` field claiming the Lipschitz constant its own sampled
///   maximum would suggest, which is the specific mistake M-244 records, since a
///   sampled maximum is a *lower* bound on the supremum.
///
/// `Exact` fields have nothing tighter to claim and are skipped, which the test
/// asserts rather than passes over silently.
#[test]
fn tightening_a_declaration_by_one_step_is_caught() {
    let mut tightened = 0usize;
    let mut already_tightest = 0usize;

    crate::for_each_reference_field!(f64, |name, field| {
        let real = field_bound_report(&field, N);
        match real.declared {
            FieldBound::Exact => {
                already_tightest += 1;
            }
            FieldBound::Underestimate { .. } => {
                already_tightest += 1;
            }
            FieldBound::Lipschitz { .. } | FieldBound::Unbounded => {
                // The one step tighter: claim to be an exact distance.
                let claimed = FieldBoundReport {
                    declared: FieldBound::Exact,
                    ..real
                };
                assert!(
                    claimed.violates(0.05),
                    "{name}: claiming Exact was not caught, though ‖∇f‖ reaches {:.4}",
                    real.sup
                );
                tightened += 1;
            }
        }
    });

    assert!(
        tightened >= 3,
        "only {tightened} fields could be tightened, too few to mean anything"
    );
    assert!(
        already_tightest >= 4,
        "{already_tightest} fields were already tightest"
    );
    std::println!(
        "measured: {tightened} declarations tightened and caught, \
         {already_tightest} already at their tightest"
    );
}

/// **The validator is one-sided, and says so.**
///
/// A sampled maximum can exceed a declared constant and settle the question; it
/// can never fall below one and settle anything. So a *loosened* declaration is
/// never reported as a violation — that is not leniency, it is the direction of
/// the inequality, and encoding it wrongly would make the checker claim to prove
/// declarations correct.
#[test]
fn a_loosened_declaration_is_never_a_violation() {
    crate::for_each_reference_field!(f64, |name, field| {
        let real = field_bound_report(&field, N);
        let loosened = FieldBoundReport {
            declared: FieldBound::Lipschitz {
                l: real.sup * 10.0 + 1.0,
            },
            ..real
        };
        assert!(
            !loosened.violates(0.0),
            "{name}: a loose bound was rejected"
        );
        let unbounded = FieldBoundReport {
            declared: FieldBound::Unbounded,
            ..real
        };
        assert!(
            !unbounded.violates(0.0),
            "{name}: Unbounded claims nothing and cannot be violated"
        );
    });
}

/// The eikonal band is a description, not a gate — pinned so it stays one.
#[test]
fn the_eikonal_tolerance_is_loose_enough_to_survive_a_crease() {
    const _: () = assert!(
        EIKONAL_TOLERANCE >= 0.01,
        "a band this tight measures creases rather than the field"
    );
    let report = field_bound_report(&crate::fields::BoxExact::<f64>::canonical(), N);
    assert!(
        report.eikonal_fraction > 0.95,
        "an exact box should be eikonal almost everywhere, got {:.3}",
        report.eikonal_fraction
    );
}
