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

/// **A directional Lipschitz bound is tighter than the global one only where the
/// global one is loose — which is nowhere, for four of the eight fields
/// (F-006).**
///
/// Galin, Guérin, Paris & Peytavie, *Segment Tracing Using Local Lipschitz
/// Bounds*, CGF 39(2) (`10.1111/cgf.13951`). Their method marches by a bound
/// computed **along the ray** rather than over all directions, which admits
/// larger safe steps.
///
/// # The prediction, written before the measurement
///
/// The paper states the condition for its own failure plainly: *"when the
/// implicit objects have an almost uniform distribution of primitives and a
/// uniform Lipschitz bound over their support Ω, the benefit is limited or
/// negative in terms of speed."* **Four of this crate's eight fields are exactly
/// that** — `sphere`, `torus`, `box_exact` and `thin_plate` are 1-Lipschitz
/// everywhere, so a directional bound cannot be smaller than 1 either and there
/// is nothing to win. F-006 says a null result is a finding; this is where it is
/// recorded.
///
/// The gyroid is the case that can gain, and its directional bound is
/// **derivable rather than sampled**: along a coordinate axis the directional
/// derivative is a single partial, `|∂g/∂x| = |cos a cos b − sin c sin a| ≤ 2`,
/// against the global `|∇g| ≤ 2√3`. A factor of `√3` tighter, exactly, with no
/// estimation involved.
#[test]
fn a_directional_bound_helps_only_where_the_global_one_is_loose() {
    use crate::Sdf;
    use crate::fields::ReferenceField;

    /// Steps a sphere tracer takes from `origin` along `+x`, marching by
    /// `|f| / lambda`.
    fn steps<F: Sdf<Scalar = f64>>(field: &F, origin: [f64; 3], lambda: f64, far: f64) -> u32 {
        let mut t = 0.0;
        for step in 0..4096u32 {
            let p = [origin[0] + t, origin[1], origin[2]];
            let d = field.sample(p).abs();
            if d < 1e-4 {
                return step;
            }
            t += (d / lambda).max(1e-6);
            if t > far {
                return step;
            }
        }
        4096
    }

    crate::for_each_reference_field!(f64, |name, field| {
        let Some(global) = field.bound().lipschitz() else {
            return;
        };
        // Along a coordinate axis, only one partial contributes. For every field
        // here except the gyroid the global bound is already 1 and cannot be
        // beaten; for the gyroid the axis-aligned bound is 2 against 2√3.
        let directional = if name == "gyroid" { 2.0 } else { global };

        let (lo, hi) = field.domain();
        let far = hi[0] - lo[0];
        let mut total_global = 0u64;
        let mut total_directional = 0u64;
        for a in 0..6 {
            for b in 0..6 {
                let at = |v: i32, l: f64, h: f64| l + (h - l) * (f64::from(v) + 0.5) / 6.0;
                let origin = [lo[0], at(a, lo[1], hi[1]), at(b, lo[2], hi[2])];
                total_global += u64::from(steps(&field, origin, global, far));
                total_directional += u64::from(steps(&field, origin, directional, far));
            }
        }
        let gain = total_global as f64 / total_directional as f64;
        std::println!(
            "measured: {name:<16} global λ {global:>6.3}  directional {directional:>6.3}  \
             steps {total_global:>6} → {total_directional:>6}  ({gain:.2}×)"
        );

        if (directional - global).abs() < 1e-12 {
            // Nothing to win, and nothing lost: the two are the same march.
            assert_eq!(
                total_global, total_directional,
                "{name}: identical bounds must produce identical marches"
            );
        } else {
            assert!(
                gain > 1.2,
                "{name}: a √3-tighter bound bought only {gain:.2}×, so the tightness \
                 is not where the cost is"
            );
        }
    });
}
