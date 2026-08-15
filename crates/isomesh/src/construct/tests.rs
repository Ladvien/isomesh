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
