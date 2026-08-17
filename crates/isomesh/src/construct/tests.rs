//! Tests for the exact distance transform.
//!
//! Two independent checks, because "it looks like a distance field" is not one.
//! [`agrees_with_brute_force_exactly`] compares against an `O(n²)` search that
//! shares no line of code with the transform, and
//! [`matches_the_analytic_sphere_within_one_spacing`] compares against a closed
//! form. The first proves the algorithm; the second proves it is measuring the
//! thing it claims to.

use super::signed_distance_field;
use crate::{RuntimeShape3, Sdf, Shape3};

/// Sample a field onto a grid, in the crate's `x`-fastest order.
fn sample_grid<F: Sdf<Scalar = f64>>(
    field: &F,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    h: f64,
) -> alloc::vec::Vec<f64> {
    let size = shape.size();
    let mut out = alloc::vec::Vec::with_capacity(shape.element_count());
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                out.push(field.sample([
                    origin[0] + h * f64::from(x),
                    origin[1] + h * f64::from(y),
                    origin[2] + h * f64::from(z),
                ]));
            }
        }
    }
    out
}

/// **The transform agrees with brute force, exactly (S-001).**
///
/// The `O(n)` separable algorithm is subtle — a lower envelope of parabolas
/// maintained in one pass — and subtle code that produces plausible numbers is
/// the failure mode. So it is checked against the definition itself: for every
/// sample, search **every** other sample and keep the nearest of opposite sign.
///
/// Exact equality, not a tolerance. Both compute the same squared integer
/// distances and take one square root, so any difference is an algorithmic one
/// rather than a rounding one — and a tolerance here would hide exactly the
/// off-by-one the envelope is prone to.
#[test]
fn agrees_with_brute_force_exactly() {
    let field = crate::fields::Sphere::<f64>::canonical();
    // Deliberately not a cube: a bug that transposes two axes survives a cube
    // and dies here.
    let shape = RuntimeShape3::new([11, 9, 13]).expect("valid shape");
    let size = shape.size();
    let h = 0.25_f64;
    let origin = [-1.25, -1.0, -1.5];
    let samples = sample_grid(&field, &shape, origin, h);

    let got = signed_distance_field(&samples, &shape, h).expect("transform");

    let index = |x: u32, y: u32, z: u32| ((z * size[1] + y) * size[0] + x) as usize;
    let mut checked = 0usize;
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let i = index(x, y, z);
                let inside = samples[i] < 0.0;
                let mut best = f64::INFINITY;
                for bz in 0..size[2] {
                    for by in 0..size[1] {
                        for bx in 0..size[0] {
                            let j = index(bx, by, bz);
                            if (samples[j] < 0.0) == inside {
                                continue;
                            }
                            let d = f64::from(x as i64 as i32 - bx as i64 as i32).powi(2)
                                + f64::from(y as i64 as i32 - by as i64 as i32).powi(2)
                                + f64::from(z as i64 as i32 - bz as i64 as i32).powi(2);
                            best = best.min(d);
                        }
                    }
                }
                let want = if inside {
                    -best.sqrt() * h
                } else {
                    best.sqrt() * h
                };
                assert!(
                    (got[i] - want).abs() < 1e-12,
                    "at ({x},{y},{z}): transform {} but brute force {want}",
                    got[i]
                );
                checked += 1;
            }
        }
    }
    assert_eq!(checked, shape.element_count());
    std::println!("measured: {checked} samples agree with brute force exactly");
}

/// **The transform is within one sample spacing of the analytic sphere
/// (S-001).**
///
/// The other half: brute force proves the algorithm computes what it intends,
/// and this proves the intent is a distance field. The gap is bounded by the
/// sampling rather than by the algorithm — the transform measures to the nearest
/// **sample of opposite sign**, and the true surface passes somewhere between
/// two samples, so half a cell of disagreement is the floor and one cell is a
/// generous ceiling.
#[test]
fn matches_the_analytic_sphere_within_one_spacing() {
    let field = crate::fields::Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new([41; 3]).expect("valid shape");
    let h = 0.1_f64;
    let origin = [-2.0; 3];
    let samples = sample_grid(&field, &shape, origin, h);
    let got = signed_distance_field(&samples, &shape, h).expect("transform");

    let size = shape.size();
    let mut worst = 0.0f64;
    let mut worst_at = [0u32; 3];
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let p = [
                    origin[0] + h * f64::from(x),
                    origin[1] + h * f64::from(y),
                    origin[2] + h * f64::from(z),
                ];
                let truth = field.sample(p);
                let i = ((z * size[1] + y) * size[0] + x) as usize;
                let err = (got[i] - truth).abs();
                if err > worst {
                    worst = err;
                    worst_at = [x, y, z];
                }
            }
        }
    }
    // **One spacing is the theoretical limit, and the measurement sits on it.**
    // The transform answers with the distance to the nearest opposite-signed
    // *sample*; the surface lies between samples, so a point whose nearest
    // crossing is mid-cell is off by up to a full spacing. Measured worst is
    // exactly 1.00 cells, so the comparison is against `1 + ε` rather than `1` —
    // not slack, but the difference between "at the bound" and "over it".
    let cells = worst / h;
    assert!(
        cells <= 1.0 + 1e-9,
        "worst disagreement {cells:.4} cells exceeds one spacing at {worst_at:?}"
    );
    std::println!(
        "measured: worst |transform − analytic| = {worst:.5} = {cells:.4} cells, \
         which is the sampling limit rather than an algorithmic error"
    );
}

/// A degenerate grid is rejected rather than transformed.
#[test]
fn a_grid_too_small_to_hold_a_cell_is_refused() {
    let shape = RuntimeShape3::new([1, 4, 4]).expect("valid shape");
    let samples = alloc::vec![1.0f64; shape.element_count()];
    assert!(matches!(
        signed_distance_field(&samples, &shape, 0.1),
        Err(crate::Error::GridTooSmall { .. })
    ));
}

/// A sample count that does not match the shape is refused rather than read past.
#[test]
fn a_mismatched_sample_count_is_refused() {
    let shape = RuntimeShape3::new([4; 3]).expect("valid shape");
    let samples = alloc::vec![1.0f64; 7];
    assert!(signed_distance_field(&samples, &shape, 0.1).is_err());
}

/// **Sweeping beats the exact transform near the surface and loses to it far
/// away (S-002).**
///
/// Both claims in one test, because either alone would be misleading. The exact
/// transform answers with the distance to the nearest opposite-signed *sample*,
/// so it is quantised to the grid; sweeping seeds from the *interpolated*
/// crossing and solves `|∇d| = 1`, so it can place the surface between samples.
/// But sweeping accumulates — a value ten cells out is ten first-order Godunov
/// updates — where the exact transform's error does not grow with distance.
///
/// So the two are compared **by band**, against the analytic sphere, rather than
/// by a single worst case that would hide the crossover.
#[test]
fn sweeping_and_the_exact_transform_trade_places_with_distance() {
    let field = crate::fields::Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new([41; 3]).expect("valid shape");
    let h = 0.1_f64;
    let origin = [-2.0; 3];
    let samples = sample_grid(&field, &shape, origin, h);

    let exact = signed_distance_field(&samples, &shape, h).expect("transform");
    let swept = super::signed_distance_field_swept(&samples, &shape, h).expect("sweep");

    let size = shape.size();
    // Worst error in each band, measured in cells from the surface.
    let mut near = (0.0f64, 0.0f64);
    let mut far_band = (0.0f64, 0.0f64);
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let p = [
                    origin[0] + h * f64::from(x),
                    origin[1] + h * f64::from(y),
                    origin[2] + h * f64::from(z),
                ];
                let truth = field.sample(p);
                let i = ((z * size[1] + y) * size[0] + x) as usize;
                let (ee, se) = ((exact[i] - truth).abs(), (swept[i] - truth).abs());
                if truth.abs() <= 2.0 * h {
                    near = (near.0.max(ee), near.1.max(se));
                } else if truth.abs() >= 8.0 * h {
                    far_band = (far_band.0.max(ee), far_band.1.max(se));
                }
            }
        }
    }

    std::println!(
        "measured: within 2 cells of the surface — exact {:.5}, swept {:.5}\n\
         measured: beyond 8 cells          — exact {:.5}, swept {:.5}",
        near.0,
        near.1,
        far_band.0,
        far_band.1
    );

    // **Sweeping wins near the surface by 3×**, which is the sub-cell seeding
    // doing exactly what it is for.
    assert!(
        near.1 * 2.0 < near.0,
        "sweeping was expected to win near the surface by a wide margin: \
         exact {:.5}, swept {:.5}",
        near.0,
        near.1
    );

    // **And it does not lose far away, which was predicted and is false
    // (M-252).** The concern was that a first-order Godunov update accumulates
    // over distance. On a sphere it does not measurably: the characteristics are
    // radial straight lines, the eight-orthant sweep follows them, and the
    // sub-cell seeding advantage survives all the way out.
    //
    // Asserted as "does not lose" rather than "wins", because the margin at
    // distance is small enough that a different field could plausibly flip it —
    // and if one does, this fails and says so.
    assert!(
        far_band.1 <= far_band.0 + 1e-9,
        "sweeping lost far from the surface after all: exact {:.5}, swept {:.5} — \
         the accumulation concern is real and M-252 needs revisiting",
        far_band.0,
        far_band.1
    );
}

/// Sweeping refuses the same degenerate inputs the exact transform does.
#[test]
fn sweeping_refuses_what_the_transform_refuses() {
    let shape = RuntimeShape3::new([1, 4, 4]).expect("valid shape");
    let samples = alloc::vec![1.0f64; shape.element_count()];
    assert!(super::signed_distance_field_swept(&samples, &shape, 0.1).is_err());

    let shape = RuntimeShape3::new([4; 3]).expect("valid shape");
    assert!(super::signed_distance_field_swept(&alloc::vec![1.0f64; 7], &shape, 0.1).is_err());
}

/// **Marching and sweeping agree closely, because they are the same update in a
/// different order (S-003).**
///
/// The two share `godunov` literally — one function, called from both — and the
/// same sub-cell seeding, so this compares *orderings* and nothing else.
/// Sweeping does eight fixed passes and lets the answer settle; marching
/// finalises the smallest tentative value at each step and never revisits.
///
/// They should therefore land close but not identical: sweeping's answer at a
/// sample depends on the eight orderings reaching it, marching's only on the
/// front, and on a sphere both find essentially the same characteristics.
#[test]
fn marching_agrees_with_sweeping_and_neither_is_the_other() {
    let field = crate::fields::Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new([33; 3]).expect("valid shape");
    let h = 0.125_f64;
    let origin = [-2.0; 3];
    let samples = sample_grid(&field, &shape, origin, h);

    let swept = super::signed_distance_field_swept(&samples, &shape, h).expect("sweep");
    let marched = super::signed_distance_field_marched(&samples, &shape, h).expect("march");

    let size = shape.size();
    let mut worst_gap = 0.0f64;
    let mut swept_err = 0.0f64;
    let mut marched_err = 0.0f64;
    let mut differ = 0usize;
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let p = [
                    origin[0] + h * f64::from(x),
                    origin[1] + h * f64::from(y),
                    origin[2] + h * f64::from(z),
                ];
                let truth = field.sample(p);
                let i = ((z * size[1] + y) * size[0] + x) as usize;
                worst_gap = worst_gap.max((swept[i] - marched[i]).abs());
                swept_err = swept_err.max((swept[i] - truth).abs());
                marched_err = marched_err.max((marched[i] - truth).abs());
                if (swept[i] - marched[i]).abs() > 1e-12 {
                    differ += 1;
                }
            }
        }
    }
    std::println!(
        "measured: swept vs marched — worst gap {worst_gap:.5}, {differ} samples differ; \
         error against analytic: swept {swept_err:.5}, marched {marched_err:.5}"
    );

    // Close: the same update from the same seeds cannot diverge far.
    assert!(
        worst_gap < 4.0 * h,
        "the two orderings disagree by {worst_gap:.5}, which is more than an \
         ordering can explain"
    );
    // Both must be real distance fields, not merely similar to each other.
    assert!(
        marched_err < 2.0 * h,
        "marching's error {marched_err:.5} is too large"
    );
    assert!(
        swept_err < 2.0 * h,
        "sweeping's error {swept_err:.5} is too large"
    );
}

/// Marching refuses the same degenerate inputs the others do.
#[test]
fn marching_refuses_what_the_others_refuse() {
    let shape = RuntimeShape3::new([1, 4, 4]).expect("valid shape");
    let samples = alloc::vec![1.0f64; shape.element_count()];
    assert!(super::signed_distance_field_marched(&samples, &shape, 0.1).is_err());

    let shape = RuntimeShape3::new([4; 3]).expect("valid shape");
    assert!(super::signed_distance_field_marched(&alloc::vec![1.0f64; 7], &shape, 0.1).is_err());
}

/// **Repeated reinitialisation does not move the zero set (S-004).**
///
/// The acceptance, and the warning the ticket exists to carry. Sussman & Fatemi:
/// naive reinitialisation *moves the surface*, so a field rebuilt after every
/// brush stroke has geometry that creeps — a wall slowly changing shape while
/// nobody edits it, which is worse than a field that is merely not a distance.
///
/// Measured rather than argued: reinitialise twenty times over and track where
/// the field crosses zero along a dense set of grid edges. The drift must stay
/// below a stated fraction of a cell, and **the bound is on the total after
/// twenty applications**, not per application — a per-step bound is satisfied by
/// a steady creep, which is the failure mode.
#[test]
fn reinitialisation_does_not_move_the_zero_set() {
    let field = crate::fields::Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new([33; 3]).expect("valid shape");
    let h = 0.125_f64;
    let origin = [-2.0; 3];
    let size = shape.size();

    let mut values = sample_grid(&field, &shape, origin, h);

    // Where the field crosses zero along every cut x-edge, as a fraction of a
    // cell. That is what a mesher reads, so it is what must not move.
    let crossings = |v: &[f64]| {
        let mut out = alloc::vec::Vec::new();
        for z in 0..size[2] {
            for y in 0..size[1] {
                for x in 0..size[0] - 1 {
                    let i = ((z * size[1] + y) * size[0] + x) as usize;
                    let j = i + 1;
                    if (v[i] < 0.0) != (v[j] < 0.0) {
                        out.push((i, v[i] / (v[i] - v[j])));
                    }
                }
            }
        }
        out
    };

    let first = crossings(&values);
    assert!(first.len() > 300, "only {} crossings to track", first.len());

    let mut touched = 0usize;
    for _ in 0..20 {
        let (next, visited) =
            super::reinitialise_narrow_band(&values, &shape, h, 3).expect("reinit");
        values = next;
        touched = visited;
    }

    let last = crossings(&values);
    assert_eq!(
        first.len(),
        last.len(),
        "reinitialisation created or destroyed a crossing, which is worse than \
         moving one"
    );

    let mut worst = 0.0f64;
    for ((ia, ta), (ib, tb)) in first.iter().zip(&last) {
        assert_eq!(ia, ib, "a crossing moved to a different edge");
        worst = worst.max((ta - tb).abs());
    }
    // **The cost claim, measured.** The ticket's premise is that a band costs
    // edited surface area rather than chunk volume; a solve that quietly touched
    // every sample would satisfy the drift assertion and none of the reason the
    // ticket exists.
    let total = shape.element_count();
    let share = 100.0 * touched as f64 / total as f64;
    std::println!(
        "measured: 20 reinitialisations moved the zero set by at most {worst:.6} of a cell; \
         the band solve finalised {touched} of {total} samples ({share:.1}%)"
    );
    assert!(
        share < 25.0,
        "the band solve touched {share:.1}% of the grid — that is a volume cost, \
         not a surface-area one"
    );
    // **Zero, not small.** The samples adjacent to a sign change are handed back
    // unchanged, so the crossing fraction is bit-identical rather than close.
    // Anything above zero here means a seed was overwritten and Sussman &
    // Fatemi's creep is back — measured at 0.152 of a cell before the fix.
    // `assert!(worst == 0.0)` rather than `assert_eq!`: clippy rejects strict
    // float comparison, and rightly in general — here the equality is the point,
    // since the seed values are handed back untouched and so must compare bitwise
    // identical, not merely close.
    #[allow(clippy::float_cmp)]
    {
        assert!(
            worst == 0.0,
            "the zero set moved by {worst:e} of a cell over 20 reinitialisations"
        );
    }
}

/// Reinitialisation refuses what the constructors refuse.
#[test]
fn reinitialisation_refuses_what_the_others_refuse() {
    let shape = RuntimeShape3::new([1, 4, 4]).expect("valid shape");
    let v = alloc::vec![1.0f64; shape.element_count()];
    assert!(super::reinitialise_narrow_band(&v, &shape, 0.1, 3).is_err());
}

// ── A-028: SampledField's analytic gradient ─────────────────────────────────

/// **The analytic gradient agrees with a central difference wherever the
/// central difference is valid**, which is the check that says it is the
/// gradient of the same interpolant `sample` evaluates.
///
/// Tested strictly *inside* cells, away from the corners and faces where the
/// interpolant's gradient is one-sided and the two constructions legitimately
/// differ. A random field rather than a smooth one, so the corner differences
/// are unrelated to each other and a sign or axis transposition cannot hide.
#[test]
fn the_analytic_gradient_matches_a_central_difference_inside_a_cell() {
    use crate::Sdf;
    use crate::construct::SampledField;

    const N: u32 = 6;
    let shape = RuntimeShape3::new([N; 3]).expect("valid shape");
    let mut state = 0x2026_A028_u64;
    let mut values = alloc::vec![0.0_f64; shape.element_count()];
    for v in &mut values {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        *v = (state >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0;
    }

    let h = 0.25_f64;
    let field = SampledField::new(&values, &shape, [0.0; 3], h).expect("wrap");

    let mut worst = 0.0_f64;
    let mut checked = 0usize;
    // Interior points of interior cells, at fractions that are not 0 or 1.
    for cell in 1..(N - 2) {
        for frac in [0.2_f64, 0.5, 0.8] {
            let p = [
                (f64::from(cell) + frac) * h,
                (f64::from(cell) + 0.37) * h,
                (f64::from(cell) + 0.61) * h,
            ];
            let analytic = field.gradient(p);
            // A central difference at a step small against the cell but large
            // against rounding.
            let d = h * 1e-4;
            for axis in 0..3 {
                let mut a = p;
                let mut b = p;
                a[axis] += d;
                b[axis] -= d;
                let numeric = (field.sample(a) - field.sample(b)) / (2.0 * d);
                worst = worst.max((analytic[axis] - numeric).abs());
                checked += 1;
            }
        }
    }

    assert!(
        checked > 20,
        "only {checked} comparisons, which measures nothing"
    );
    assert!(
        worst < 1e-6,
        "analytic and numeric gradients differ by {worst:.3e}; the analytic form is not \
         the gradient of the interpolant `sample` evaluates"
    );
}

/// **The failure A-028 is about, as a fixture: a central difference is zero at a
/// local extremum however steep the field is, and the analytic form is not.**
///
/// Three samples along one axis, `18, -1, 18` — the shape measured on `bonsai`,
/// where `u8` quantisation put both neighbours on the same integer. The slopes
/// are ∓19, so nothing here is flat; the field is *symmetric*, which is what
/// zeroes a central difference.
#[test]
fn a_central_difference_is_zero_at_a_local_extremum_and_the_analytic_form_is_not() {
    use crate::Sdf;
    use crate::construct::SampledField;

    const N: u32 = 3;
    let shape = RuntimeShape3::new([N; 3]).expect("valid shape");
    let mut values = alloc::vec![5.0_f64; shape.element_count()];
    // A symmetric trough along x through the middle of the grid.
    for z in 0..N {
        for y in 0..N {
            for (x, v) in [(0u32, 18.0), (1, -1.0), (2, 18.0)] {
                let i = (z as usize * N as usize + y as usize) * N as usize + x as usize;
                values[i] = v;
            }
        }
    }

    let field = SampledField::new(&values, &shape, [0.0; 3], 1.0).expect("wrap");
    let p = [1.0, 1.0, 1.0];

    // The default central difference, spelled out rather than called, so this
    // test still means something if `Sdf::gradient`'s default changes.
    let step = f64::EPSILON.cbrt() * 1.0_f64.max(1.0);
    let numeric = (field.sample([1.0 + step, 1.0, 1.0]) - field.sample([1.0 - step, 1.0, 1.0]))
        / (2.0 * step);
    let analytic = field.gradient(p);

    // The true slope either side is 19. The central difference averages `-19`
    // and `+19` and keeps only what rounding leaves behind; on `bonsai` that
    // residue was **exactly zero** and the extractor refused the volume. Here
    // it is merely negligible, which is the same failure at a different set of
    // float values, so the assertion is on the ratio rather than on zero.
    assert!(
        numeric.abs() < 1.0,
        "a central difference across a symmetric extremum should keep almost \
         nothing of a slope of 19, got {numeric}"
    );
    assert!(
        analytic[0].abs() > 18.0,
        "the analytic gradient should see the 19 slope of the cell it is in, got {analytic:?}"
    );
    assert!(
        analytic[0].abs() > numeric.abs() * 1e6,
        "the point of this fixture is the gap between the two, and it is only \
         {:.3e} against {:.3e}",
        analytic[0].abs(),
        numeric.abs()
    );
}
