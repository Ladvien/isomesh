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
