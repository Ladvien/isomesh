//! A-012's acceptance, and the measurement the acceptance does not ask for.

use alloc::vec::Vec;

use super::*;
use crate::fields::{BoxExact, ReferenceField, Sphere, Torus};
use crate::marching_cubes::MarchingCubes;
use crate::{MeshBuffer, RuntimeShape3};

/// Mesh a reference field with Marching Cubes, returning the mesh and the cell
/// size. Marching Cubes because it puts vertices *on grid edges*, which is where
/// the three strategies are most likely to disagree.
fn mesh<F: Sdf<Scalar = f64> + ReferenceField>(field: &F, samples: u32) -> (MeshBuffer<f64>, f64) {
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let mut out = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(field, &shape, lo, h, &mut out)
        .expect("extraction");
    (out, h)
}

/// Worst and mean angle between two normal sets, in degrees.
fn deviation(a: &[[f64; 3]], b: &[[f64; 3]]) -> (f64, f64) {
    assert_eq!(a.len(), b.len());
    let mut worst = 0.0f64;
    let mut total = 0.0f64;
    for (x, y) in a.iter().zip(b) {
        let d = vec3::dot(*x, *y).clamp(-1.0, 1.0);
        let angle = d.acos().to_degrees();
        worst = worst.max(angle);
        total += angle;
    }
    (worst, total / a.len() as f64)
}

// ─── the acceptance criterion ───────────────────────────────────────────────

/// **A-012, first half:** all three produce unit-length normals.
///
/// Checked on every strategy against three fields at two resolutions, and at
/// `f64`'s own tolerance rather than a generous one — a normal that is only
/// approximately unit length makes every lighting calculation downstream
/// approximately wrong, and there is no reason for it to be.
#[test]
fn every_strategy_produces_unit_length_normals() {
    let mut checked = 0;
    for samples in [17u32, 33] {
        for (name, mut buffer, h, analytic, central, area) in [
            {
                let (m, h) = mesh(&Sphere::<f64>::canonical(), samples);
                ("sphere", m, h, true, true, true)
            },
            {
                let (m, h) = mesh(&Torus::<f64>::canonical(), samples);
                ("torus", m, h, true, true, true)
            },
            {
                let (m, h) = mesh(&BoxExact::<f64>::canonical(), samples);
                ("box_exact", m, h, true, true, true)
            },
        ] {
            let field_of = |n: &str| -> alloc::boxed::Box<dyn Sdf<Scalar = f64>> {
                match n {
                    "sphere" => alloc::boxed::Box::new(Sphere::<f64>::canonical()),
                    "torus" => alloc::boxed::Box::new(Torus::<f64>::canonical()),
                    _ => alloc::boxed::Box::new(BoxExact::<f64>::canonical()),
                }
            };
            let field = field_of(name);

            let mut strategies: Vec<NormalStrategy<f64>> = Vec::new();
            if analytic {
                strategies.push(NormalStrategy::AnalyticGradient);
            }
            if central {
                strategies.push(NormalStrategy::CentralDifference { step: h });
            }
            if area {
                strategies.push(NormalStrategy::AreaWeightedFaces);
            }

            for strategy in strategies {
                recompute(&mut buffer, &field, strategy).expect("normals");
                for (v, n) in buffer.normals.iter().enumerate() {
                    let len = vec3::length(*n);
                    assert!(
                        (len - 1.0).abs() < 1e-12,
                        "{name} at {samples}^3, {strategy:?}: vertex {v} normal length {len}"
                    );
                }
                checked += 1;
            }
        }
    }
    assert_eq!(checked, 2 * 3 * 3);
}

/// **A-012, second half:** analytic and central-difference agree within
/// tolerance on `sphere`.
///
/// The tolerance is stated as an angle rather than a component difference,
/// because a normal is a direction and the number a renderer cares about is how
/// far off that direction is. Both are reported so the number is on the record
/// rather than only the pass.
#[test]
fn the_analytic_and_differenced_gradients_agree_on_a_sphere() {
    for samples in [17u32, 33, 65] {
        let (mut buffer, h) = mesh(&Sphere::<f64>::canonical(), samples);
        let field = Sphere::<f64>::canonical();

        recompute(&mut buffer, &field, NormalStrategy::AnalyticGradient).expect("normals");
        let analytic = buffer.normals.clone();

        recompute(
            &mut buffer,
            &field,
            NormalStrategy::CentralDifference { step: h },
        )
        .expect("normals");

        let (worst, mean) = deviation(&analytic, &buffer.normals);
        std::println!(
            "sphere {samples}^3, step = h = {h:.4}: worst {worst:.6} deg, mean {mean:.6} deg"
        );
        assert!(
            worst < 1.0,
            "sphere at {samples}^3: worst deviation {worst} degrees"
        );
    }
}

// ─── what the acceptance does not ask for ───────────────────────────────────

/// The post-pass must be the same path extraction already took.
///
/// [`NormalStrategy::AnalyticGradient`] is exactly what every extractor does
/// inline, so recomputing with it has to give back **bit-identical** normals. If
/// this ever fails, the two have diverged and one of them is now producing
/// something nobody chose.
#[test]
fn recomputing_with_the_analytic_gradient_reproduces_extraction() {
    // Exact comparison on purpose: this pins that two code paths are one.
    #![allow(clippy::float_cmp)]
    for samples in [17u32, 33] {
        let (mut buffer, _) = mesh(&Sphere::<f64>::canonical(), samples);
        let extracted = buffer.normals.clone();
        recompute(
            &mut buffer,
            &Sphere::<f64>::canonical(),
            NormalStrategy::AnalyticGradient,
        )
        .expect("normals");
        assert_eq!(
            extracted, buffer.normals,
            "the post-pass and the extractor disagree at {samples}^3"
        );
    }
}

/// Where the mesh and the field disagree, and by how much.
///
/// Area-weighted normals come from the geometry, so on a smooth field they track
/// the field closely and on a **sharp** one they cannot: a box's corner vertex
/// gets the average of three face normals where the field's gradient gives one
/// of them. That is not a defect in either — it is the difference between "which
/// way does the surface face" and "which way does the field increase", and it is
/// the whole reason this is selectable.
///
/// Recorded as a measurement, gated only where the direction of the effect is
/// certain.
#[test]
fn area_weighted_normals_track_the_field_on_smooth_geometry_and_not_on_sharp() {
    let mut rows: Vec<(&str, u32, f64, f64)> = Vec::new();
    for samples in [17u32, 33, 65] {
        for name in ["sphere", "torus", "box_exact"] {
            let (mut buffer, _) = match name {
                "sphere" => mesh(&Sphere::<f64>::canonical(), samples),
                "torus" => mesh(&Torus::<f64>::canonical(), samples),
                _ => mesh(&BoxExact::<f64>::canonical(), samples),
            };
            let field: alloc::boxed::Box<dyn Sdf<Scalar = f64>> = match name {
                "sphere" => alloc::boxed::Box::new(Sphere::<f64>::canonical()),
                "torus" => alloc::boxed::Box::new(Torus::<f64>::canonical()),
                _ => alloc::boxed::Box::new(BoxExact::<f64>::canonical()),
            };

            recompute(&mut buffer, &field, NormalStrategy::AnalyticGradient).expect("normals");
            let analytic = buffer.normals.clone();
            recompute(&mut buffer, &field, NormalStrategy::AreaWeightedFaces).expect("normals");
            let (worst, mean) = deviation(&analytic, &buffer.normals);
            rows.push((name, samples, worst, mean));
        }
    }

    for (name, samples, worst, mean) in &rows {
        std::println!("{name} {samples}^3: worst {worst:.3} deg, mean {mean:.3} deg from analytic");
    }

    // The direction of the effect is the assertion; the magnitudes are the
    // measurement. A box is all corners and edges, a sphere has neither.
    let mean_of = |want: &str, n: u32| {
        rows.iter()
            .find(|r| r.0 == want && r.1 == n)
            .map(|r| r.3)
            .expect("row")
    };
    for samples in [17u32, 33, 65] {
        assert!(
            mean_of("box_exact", samples) > mean_of("sphere", samples),
            "at {samples}^3 the box does not disagree more than the sphere"
        );
    }

    // **The mechanism, as an assertion.** Refining a grid does not soften a
    // corner, so the *worst* disagreement on a sharp field is a property of the
    // feature and not of the resolution — measured identical to six figures at
    // 17^3, 33^3 and 65^3. On a smooth field it falls with resolution, because
    // there the disagreement is discretisation rather than geometry.
    let worst_of = |want: &str, n: u32| {
        rows.iter()
            .find(|r| r.0 == want && r.1 == n)
            .map(|r| r.2)
            .expect("row")
    };
    let box_worst: Vec<f64> = [17u32, 33, 65]
        .iter()
        .map(|n| worst_of("box_exact", *n))
        .collect();
    for w in &box_worst {
        assert!(
            (w - box_worst[0]).abs() < 1e-6,
            "the box's worst disagreement moved with resolution: {box_worst:?}"
        );
    }
    let sphere_worst: Vec<f64> = [17u32, 33, 65]
        .iter()
        .map(|n| worst_of("sphere", *n))
        .collect();
    assert!(
        sphere_worst[2] < sphere_worst[0] * 0.75,
        "the sphere's worst disagreement did not fall with resolution: {sphere_worst:?}"
    );
}

/// The step is the measurement, so a meaningless one is refused at the door
/// rather than producing a black mesh.
#[test]
fn a_meaningless_differencing_step_is_rejected() {
    let (mut buffer, _) = mesh(&Sphere::<f64>::canonical(), 17);
    let field = Sphere::<f64>::canonical();
    for step in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let result = recompute(
            &mut buffer,
            &field,
            NormalStrategy::CentralDifference { step },
        );
        assert!(
            matches!(result, Err(crate::Error::InvalidCellSize { .. })),
            "step {step} was accepted"
        );
    }
}

/// A vertex with no incident area has no normal, and that is reported rather
/// than substituted.
#[test]
fn a_vertex_with_no_incident_area_is_an_error() {
    let mut buffer = MeshBuffer::<f64>::new();
    // Two triangles' worth of vertices, but only one degenerate triangle: every
    // cross product is zero, so nothing accumulates.
    for p in [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]] {
        buffer.positions.push(p);
        buffer.normals.push([0.0, 0.0, 1.0]);
    }
    buffer.indices.extend_from_slice(&[0, 1, 2]);

    let result = recompute(
        &mut buffer,
        &Sphere::<f64>::canonical(),
        NormalStrategy::AreaWeightedFaces,
    );
    assert!(
        matches!(result, Err(crate::Error::DegenerateNormal { vertex: 0 })),
        "collinear triangle accepted: {result:?}"
    );
}

/// Same input twice, byte-identical output — T-004's rule, applied to the pass
/// that now decides what a normal is.
#[test]
fn recomputation_is_deterministic() {
    // Exact comparison on purpose.
    #![allow(clippy::float_cmp)]
    let field = Torus::<f64>::canonical();
    for strategy in [
        NormalStrategy::AnalyticGradient,
        NormalStrategy::CentralDifference { step: 0.125 },
        NormalStrategy::AreaWeightedFaces,
    ] {
        let (mut a, _) = mesh(&field, 25);
        let (mut b, _) = mesh(&field, 25);
        recompute(&mut a, &field, strategy).expect("normals");
        recompute(&mut b, &field, strategy).expect("normals");
        assert_eq!(a.normals, b.normals, "{strategy:?} is not deterministic");

        // And running it twice on the same buffer must not drift.
        let once = a.normals.clone();
        recompute(&mut a, &field, strategy).expect("normals");
        assert_eq!(once, a.normals, "{strategy:?} is not idempotent");
    }
}

/// Differencing at the cell size is the voxel-grid case, and it is measurably
/// coarser than differencing at the field's own step.
///
/// Worth pinning because it is the number a game actually gets: a sampled field
/// has nothing finer than its spacing to difference over, so this is the ceiling
/// on normal quality for anyone without an analytic field.
#[test]
fn differencing_at_the_cell_size_costs_something_and_the_cost_is_recorded() {
    let mut means: Vec<f64> = Vec::new();
    for samples in [17u32, 33, 65] {
        let (mut buffer, h) = mesh(&Sphere::<f64>::canonical(), samples);
        let field = Sphere::<f64>::canonical();

        recompute(&mut buffer, &field, NormalStrategy::AnalyticGradient).expect("normals");
        let analytic = buffer.normals.clone();

        recompute(
            &mut buffer,
            &field,
            NormalStrategy::CentralDifference { step: h },
        )
        .expect("normals");
        let (coarse_worst, coarse_mean) = deviation(&analytic, &buffer.normals);

        recompute(
            &mut buffer,
            &field,
            NormalStrategy::CentralDifference { step: h * 1e-4 },
        )
        .expect("normals");
        let (fine_worst, fine_mean) = deviation(&analytic, &buffer.normals);

        std::println!(
            "sphere {samples}^3: step h -> worst {coarse_worst:.6} deg mean {coarse_mean:.6}; \
             step h/1e4 -> worst {fine_worst:.6} deg mean {fine_mean:.6}"
        );
        assert!(
            coarse_mean >= fine_mean,
            "a coarser step was not worse at {samples}^3"
        );
        means.push(coarse_mean);
    }

    // Central differences are O(h^2), so halving the spacing must cut the error
    // by roughly four. Measured 3.76x and 3.92x — the same convergence order
    // M-12 found for *position* error, now for direction.
    for pair in means.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.0..5.0).contains(&ratio),
            "convergence ratio {ratio:.3} is not h^2-like; means {means:?}"
        );
    }
}
