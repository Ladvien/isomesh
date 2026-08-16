//! T-016's acceptance: the head-to-head, and whether M-72's aliasing is the
//! predicted failure of an operator the literature already rejects.

extern crate std;

use crate::fields::ReferenceField;
use crate::marching_cubes::MarchingCubes;
use crate::mesh::MeshBuffer;
use crate::{RuntimeShape3, Sdf, Shape3};

use super::{Downsample, downsample};

/// Sample a field onto a grid, x fastest.
fn sample_grid<F: Sdf<Scalar = f64>>(
    field: &F,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    h: f64,
) -> std::vec::Vec<f64> {
    let size = shape.size();
    let mut out = std::vec::Vec::with_capacity(shape.element_count());
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

/// Mesh a grid of samples with Marching Cubes.
fn mesh(values: &[f64], shape: &RuntimeShape3, origin: [f64; 3], h: f64) -> MeshBuffer<f64> {
    let mut out = MeshBuffer::new();
    let mut mc = MarchingCubes::<f64>::new();
    let field = crate::construct::SampledField::new(values, shape, origin, h).expect("wrap");
    mc.extract(&field, shape, origin, h, &mut out)
        .expect("extraction");
    out
}

/// The shapes and shifts are right: coarse sample `i` sits on fine sample `2i`.
///
/// The property the whole module rests on, because a half-sample drift per level
/// is systematic and invisible until three levels down.
#[test]
fn coarsening_halves_the_grid_without_shifting_it() {
    let shape = RuntimeShape3::new([9; 3]).expect("valid shape");
    let size = shape.size();

    // A field linear in x. Every centred kernel reproduces a linear function
    // exactly, so any shift shows up as a constant offset.
    let mut fine = std::vec::Vec::with_capacity(shape.element_count());
    for _z in 0..size[2] {
        for _y in 0..size[1] {
            for x in 0..size[0] {
                fine.push(f64::from(x));
            }
        }
    }

    for op in Downsample::ALL {
        let (coarse, cshape) = downsample(&fine, &shape, op).expect("downsample");
        assert_eq!(cshape.size(), [5; 3], "{}", op.name());
        assert_eq!(coarse.len(), cshape.element_count());

        // Interior coarse samples only: the edge ones see clamped neighbours,
        // and `Min` legitimately reads low there.
        for z in 1..4u32 {
            for y in 1..4u32 {
                for x in 1..4u32 {
                    let got = coarse[((z * 5 + y) * 5 + x) as usize];
                    let want = f64::from(x * 2) - if op == Downsample::Min { 1.0 } else { 0.0 };
                    assert!(
                        (got - want).abs() < 1e-12,
                        "{} at ({x},{y},{z}): {got} against {want}",
                        op.name()
                    );
                }
            }
        }
    }
}

/// Non-nesting grids are refused rather than silently truncated.
#[test]
fn it_refuses_a_grid_that_does_not_halve() {
    let shape = RuntimeShape3::new([10, 9, 9]).expect("valid shape");
    let v = std::vec![0.0f64; shape.element_count()];
    assert!(downsample(&v, &shape, Downsample::Mean).is_err());

    let tiny = RuntimeShape3::new([2; 3]).expect("valid shape");
    let v = std::vec![0.0f64; tiny.element_count()];
    assert!(downsample(&v, &tiny, Downsample::Mean).is_err());
}

/// **The head-to-head (M-265).** Every operator against re-sampling, on every
/// reference field, across four levels.
///
/// Two columns carry the argument:
///
/// - **triangles** — what survives. The literature's claim is that a sub-cell
///   feature should *disappear* under re-sampling and *alias* under box
///   averaging.
/// - **worst** — the largest `|analytic field|` over the output vertices, so a
///   surface that survives but sits in the wrong place is not scored as a win.
#[test]
fn every_operator_against_resampling() {
    // 65³ down to 9³: four levels, each `2ᵏ + 1`.
    const LEVELS: usize = 4;
    let mut csv = std::string::String::from("field,operator,level,samples,triangles,worst\n");

    for op in Downsample::ALL {
        std::println!("\n=== {} ===", op.name());
        std::println!(
            "{:<16} {:>6} {:>10} {:>12} {:>10} {:>12}",
            "field",
            "level",
            "triangles",
            "resampled",
            "worst",
            "resampled"
        );

        crate::for_each_reference_field!(f64, |name, field| {
            // Inline block, so no `return` anywhere in here (M-253).
            let (lo, hi) = field.domain();

            let mut samples = 65u32;
            let mut h = (hi[0] - lo[0]) / f64::from(samples - 1);
            let mut shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
            let mut values = sample_grid(&field, &shape, lo, h);

            for level in 0..LEVELS {
                samples = shape.size()[0];
                // Re-sampling: evaluate the field at this level's spacing.
                let resampled = sample_grid(&field, &shape, lo, h);
                let a = mesh(&values, &shape, lo, h);
                let b = mesh(&resampled, &shape, lo, h);

                let worst = |m: &MeshBuffer<f64>| {
                    m.positions
                        .iter()
                        .map(|p| field.sample(*p).abs())
                        .fold(0.0f64, f64::max)
                };
                let (wa, wb) = (worst(&a), worst(&b));
                std::println!(
                    "{name:<16} {level:>6} {:>10} {:>12} {wa:>10.5} {wb:>12.5}",
                    a.indices.len() / 3,
                    b.indices.len() / 3
                );
                csv.push_str(&std::format!(
                    "{name},{},{level},{samples},{},{wa:.8}\n",
                    op.name(),
                    a.indices.len() / 3
                ));
                if op == Downsample::Decimate {
                    // Written once, not once per operator.
                    csv.push_str(&std::format!(
                        "{name},resample,{level},{samples},{},{wb:.8}\n",
                        b.indices.len() / 3
                    ));

                    // **Decimation and re-sampling are the same operator here,
                    // and this pins it (M-265).** Decimating a grid that was
                    // itself sampled from the field keeps a *subset* of those
                    // same points, and decimation composes — every level of the
                    // chain is still a subset. So the literature's "re-sample,
                    // do not downsample" is, for this operator on a nested
                    // grid, a distinction with no difference. It stops being
                    // one the moment the fine level is edited rather than
                    // sampled, which is why this is an equality and not an
                    // approximation.
                    assert_eq!(
                        a.indices, b.indices,
                        "{name} at level {level}: decimation diverged from \
                         re-sampling"
                    );
                    assert_eq!(a.positions, b.positions);
                }

                if level + 1 < LEVELS {
                    let (next, next_shape) = downsample(&values, &shape, op).expect("downsample");
                    values = next;
                    shape = next_shape;
                    h *= 2.0;
                }
            }
        });
    }

    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements/downsample.csv");
    if let Some(dir) = path.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let _ = std::fs::write(&path, csv);
    std::println!("\nwrote {}", path.display());
}

/// **M-72's prediction, tested head-on.**
///
/// The ticket asserts that `thin_plate`'s aliasing is the predicted failure of
/// box-filter averaging, and that under re-sampling the plate would *"correctly
/// disappear"* once thinner than a cell. Both halves are checked here, on the
/// same field at the same four spacings, because a prediction that names one
/// operator as the culprit is falsified the moment another operator does the
/// same thing.
#[test]
fn the_thin_plate_across_every_operator() {
    let field = crate::fields::ThinPlate::<f64>::canonical();
    let (lo, hi) = field.domain();

    std::println!(
        "{:<10} {:>8} {:>8} {:>8} {:>8} {:>10}",
        "level",
        "h",
        "decimate",
        "mean",
        "tent",
        "min"
    );

    let mut counts: std::vec::Vec<(Downsample, std::vec::Vec<usize>)> = std::vec::Vec::new();
    for op in Downsample::ALL {
        let samples = 65u32;
        let mut h = (hi[0] - lo[0]) / f64::from(samples - 1);
        let mut shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
        let mut values = sample_grid(&field, &shape, lo, h);
        let mut row = std::vec::Vec::new();
        for level in 0..4 {
            row.push(mesh(&values, &shape, lo, h).indices.len() / 3);
            if level + 1 < 4 {
                let (next, next_shape) = downsample(&values, &shape, op).expect("downsample");
                values = next;
                shape = next_shape;
                h *= 2.0;
            }
        }
        counts.push((op, row));
    }

    // Re-sampling, for the same four spacings.
    let mut resampled = std::vec::Vec::new();
    for level in 0..4u32 {
        let samples = 64 / 2u32.pow(level) + 1;
        let h = (hi[0] - lo[0]) / f64::from(samples - 1);
        let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
        let values = sample_grid(&field, &shape, lo, h);
        resampled.push(mesh(&values, &shape, lo, h).indices.len() / 3);
    }

    for level in 0..4 {
        let samples = 64 / 2u32.pow(level as u32) + 1;
        let h = (hi[0] - lo[0]) / f64::from(samples - 1);
        std::println!(
            "{level:<10} {h:>8.4} {:>8} {:>8} {:>8} {:>10}",
            counts[0].1[level],
            counts[1].1[level],
            counts[2].1[level],
            counts[3].1[level],
        );
    }
    std::println!("resample:            {resampled:?}");
    for (op, row) in &counts {
        std::println!("{:<10} {row:?}", op.name());
    }

    // Every operator must at least agree at level 0, where nothing has been
    // coarsened yet -- otherwise the comparison is between different fields.
    for (op, row) in &counts {
        assert_eq!(
            row[0],
            resampled[0],
            "{} disagrees with re-sampling before any coarsening",
            op.name()
        );
    }

    // `Min` is the conservative operator and must never lose the plate: taking
    // the minimum can only pull values negative, so a solid region can grow and
    // never shrink.
    let min_row = &counts[3].1;
    assert!(
        min_row[3] > 0,
        "the conservative `min` operator lost the plate entirely: {min_row:?}"
    );
}

/// **The mechanism behind M-72, isolated (M-266).**
///
/// `ThinPlate` is `0.4 × h` thick and centred at `y = 0`. Every grid here has an
/// *odd* sample count, so `y = 0` is always a sample plane — the plate is
/// perfectly aligned with the sampling at every level, and that is why
/// re-sampling keeps finding it however coarse the grid gets.
///
/// M-72 read the survival as *"whichever edges happen to straddle a thin slab"*,
/// which suggests chance. It is not chance. Shifting the plate by half a cell
/// puts it strictly between two sample planes, and the surface should then
/// vanish completely at every level — including the finest, where it is still
/// only 0.4 cells thick.
#[test]
fn the_aliasing_is_alignment_and_a_half_cell_shift_removes_it() {
    let base = crate::fields::ThinPlate::<f64>::canonical();
    let (lo, hi) = base.domain();

    std::println!(
        "{:<10} {:>8} {:>10} {:>10}",
        "level",
        "h",
        "centred",
        "shifted"
    );
    for level in 0..4u32 {
        let samples = 64 / 2u32.pow(level) + 1;
        let h = (hi[0] - lo[0]) / f64::from(samples - 1);
        let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

        let shifted = crate::fields::ThinPlate::<f64> {
            center: [0.0, h * 0.5, 0.0],
            half_extents: base.half_extents,
        };

        let a = mesh(&sample_grid(&base, &shape, lo, h), &shape, lo, h);
        let b = mesh(&sample_grid(&shifted, &shape, lo, h), &shape, lo, h);
        std::println!(
            "{level:<10} {h:>8.4} {:>10} {:>10}",
            a.indices.len() / 3,
            b.indices.len() / 3
        );

        assert!(
            !a.indices.is_empty(),
            "level {level}: the centred plate vanished, so the alignment claim \
             has nothing to explain"
        );
        assert!(
            b.indices.is_empty(),
            "level {level}: the plate shifted half a cell off the sample plane \
             still produced {} triangles — the survival is not alignment",
            b.indices.len() / 3
        );
    }
}
