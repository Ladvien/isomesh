//! S-005's acceptance: *"error against S-001 quantified rather than assumed;
//! 'approximate' is a measurement, not an adjective."*

use isomesh::construct::{signed_distance_field, signed_distance_field_swept};
use isomesh::fields::ReferenceField;
use isomesh::{RuntimeShape3, Sdf, Shape3};

use super::JumpFlood;
use crate::{FieldBuffer, GridParams, headless};

/// Strides halve from `n/2` to `1`, largest first.
#[test]
fn strides_halve_from_half_the_longest_axis() {
    let p = GridParams::new([65, 33, 17], [-1.0; 3], 0.03125).expect("valid grid");
    assert_eq!(JumpFlood::strides(p), std::vec![32, 16, 8, 4, 2, 1]);

    // Non-power-of-two grids still terminate at 1, which is the property that
    // matters -- a stride of 1 makes the last pass a plain 27-neighbour
    // propagation, so nothing is left unreachable by rounding.
    let odd = GridParams::new([50; 3], [-1.0; 3], 0.04).expect("valid grid");
    let s = JumpFlood::strides(odd);
    assert_eq!(s.first(), Some(&25));
    assert_eq!(s.last(), Some(&1));
}

/// **The measurement (M-257).** Jump flooding against the CPU constructors, on
/// every reference field whose analytic value is a distance.
///
/// # Two references, because S-001 alone gave the wrong answer
///
/// The ticket says *"error against S-001"*, and taking that literally produced a
/// gate the flood failed for the wrong reason: it disagreed with
/// [`signed_distance_field`] by a full cell near the surface **while being three
/// times closer to the analytic truth than the transform was** (0.082 against
/// 0.250 on `sphere` at 17³). The exact transform is exact to the nearest
/// *sample*, so near the surface it is the coarser of the two, and asserting
/// agreement with it asserts that the flood reproduce a quantisation error.
///
/// So the near-surface reference is [`signed_distance_field_swept`], which seeds
/// from the same sub-cell crossings the flood does — that comparison isolates
/// the flood's own approximation, which is what S-005 is asking about.
///
/// Columns:
///
/// - **vs_exact** / **vs_swept** — worst disagreement anywhere.
/// - **near_swept** — worst disagreement within two cells of the surface,
///   against the constructor that shares the flood's seeding. **The gate.**
/// - **flood_err** / **swept_err** / **xform_err** — each against the analytic
///   field. Ground truth, so that "disagrees" and "is wrong" stay separable.
#[test]
fn jump_flood_error_against_the_exact_transform() {
    let Ok(gpu) = headless::Gpu::new() else {
        std::eprintln!("no GPU adapter; skipping");
        return;
    };
    let flood = JumpFlood::new(gpu.device()).expect("compiles");

    let mut csv = std::string::String::from(
        "field,samples,vs_exact,vs_swept,near_swept,flood_err,swept_err,xform_err\n",
    );
    std::println!(
        "{:<16} {:>7} {:>9} {:>9} {:>10} {:>10} {:>10} {:>10}",
        "field",
        "samples",
        "vs_exact",
        "vs_swept",
        "near_swept",
        "flood_err",
        "swept_err",
        "xform_err"
    );

    isomesh::for_each_reference_field!(f32, |name, field| {
        // `if` rather than `let ... else { return }`: the macro expands inline
        // blocks, so a `return` would exit this test at the first skipped field
        // and pass while covering nothing (M-253).
        if field.bound().is_exact() {
            for samples in [17u32, 33, 65] {
                let (lo, hi) = field.domain();
                let h = (hi[0] - lo[0]) / (samples - 1) as f32;
                let params = GridParams::new([samples; 3], lo, h).expect("valid grid");
                let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

                let mut grid = std::vec::Vec::with_capacity(shape.element_count());
                for z in 0..samples {
                    for y in 0..samples {
                        for x in 0..samples {
                            grid.push(field.sample([
                                lo[0] + h * x as f32,
                                lo[1] + h * y as f32,
                                lo[2] + h * z as f32,
                            ]));
                        }
                    }
                }

                let buffer = FieldBuffer::uploaded(gpu.device(), gpu.queue(), params, &grid)
                    .expect("upload");
                let flooded = flood
                    .build(gpu.device(), gpu.queue(), &buffer)
                    .expect("flood");
                let exact = signed_distance_field(&grid, &shape, h).expect("transform");
                let swept = signed_distance_field_swept(&grid, &shape, h).expect("sweep");

                let mut vs_exact = 0.0f32;
                let mut vs_swept = 0.0f32;
                let mut near_swept = 0.0f32;
                let mut flood_err = 0.0f32;
                let mut swept_err = 0.0f32;
                let mut xform_err = 0.0f32;
                for z in 0..samples {
                    for y in 0..samples {
                        for x in 0..samples {
                            let i = ((z * samples + y) * samples + x) as usize;
                            let truth = field.sample([
                                lo[0] + h * x as f32,
                                lo[1] + h * y as f32,
                                lo[2] + h * z as f32,
                            ]);
                            vs_exact = vs_exact.max((flooded[i] - exact[i]).abs());
                            let gap = (flooded[i] - swept[i]).abs();
                            vs_swept = vs_swept.max(gap);
                            if truth.abs() <= 2.0 * h {
                                near_swept = near_swept.max(gap);
                            }
                            flood_err = flood_err.max((flooded[i] - truth).abs());
                            swept_err = swept_err.max((swept[i] - truth).abs());
                            xform_err = xform_err.max((exact[i] - truth).abs());
                        }
                    }
                }

                let total = shape.element_count();
                std::println!(
                    "{name:<16} {samples:>7} {vs_exact:>9.5} {vs_swept:>9.5} \
                     {near_swept:>10.5} {flood_err:>10.5} {swept_err:>10.5} {xform_err:>10.5}"
                );
                csv.push_str(&std::format!(
                    "{name},{samples},{vs_exact:.8},{vs_swept:.8},{near_swept:.8},\
                     {flood_err:.8},{swept_err:.8},{xform_err:.8}\n"
                ));

                // **The sign must never disagree.** Both take it from the field
                // itself, so a mismatch is a bug in this crate rather than the
                // flood's approximation, and it would invert a chunk.
                for i in 0..total {
                    assert_eq!(
                        flooded[i] < 0.0,
                        exact[i] < 0.0,
                        "{name} at {samples}³ sample {i}: flood {} vs transform {}",
                        flooded[i],
                        exact[i]
                    );
                }

                // **The gate is against ground truth, not against agreement.**
                //
                // Two gates were written before this one and both were wrong in
                // the same way: they asserted the flood *agree* with a CPU
                // constructor, which asserts it reproduce that constructor's
                // error. Against `signed_distance_field` that is a quantisation
                // to the nearest sample; against `signed_distance_field_swept` it
                // is the Godunov update's first-order truncation. `thin_plate` at
                // 17³ failed a half-cell agreement bound at 0.523 cells while
                // being **the most accurate of the three** against the analytic
                // field (0.250, versus the sweep's 0.407 and the transform's
                // 0.354).
                //
                // So the assertion is the one that matters to a consumer: an
                // approximate GPU method has no reason to exist if it loses to
                // either exact CPU one. `near_swept` stays a recorded metric.
                assert!(
                    flood_err <= xform_err,
                    "{name} at {samples}³: flood is {flood_err} from truth against \
                     the transform's {xform_err}"
                );
                assert!(
                    flood_err <= swept_err,
                    "{name} at {samples}³: flood is {flood_err} from truth against \
                     the sweep's {swept_err}"
                );
            }
        }
    });

    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements/jump_flood.csv");
    if let Some(dir) = path.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let _ = std::fs::write(&path, csv);
    std::println!("\nwrote {}", path.display());
}

/// A field with no sign change at all leaves every sample unreached.
///
/// The interesting half is that it must not produce a NaN: `resolve` writes
/// `1e30` rather than infinity for exactly the reason `construct::far` does on
/// the CPU, so the first consumer to subtract two of them gets a finite answer.
#[test]
fn a_field_with_no_surface_stays_finite() {
    let Ok(gpu) = headless::Gpu::new() else {
        std::eprintln!("no GPU adapter; skipping");
        return;
    };
    let flood = JumpFlood::new(gpu.device()).expect("compiles");
    let params = GridParams::new([9; 3], [0.0; 3], 0.1).expect("valid grid");
    let grid = std::vec![1.0f32; params.sample_count() as usize];
    let buffer = FieldBuffer::uploaded(gpu.device(), gpu.queue(), params, &grid).expect("upload");

    let out = flood
        .build(gpu.device(), gpu.queue(), &buffer)
        .expect("flood");
    assert!(out.iter().all(|v| v.is_finite() && *v > 0.0));
}
