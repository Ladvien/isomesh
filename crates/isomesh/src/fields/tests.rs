//! Tests for the reference fields.
//!
//! Three sweeps run over all seven fields through
//! [`for_each_reference_field!`](crate::for_each_reference_field), plus one
//! explicit sign test per field with hand-computed values.

// The sign tests compare against values computed by hand, exactly. An
// approximate comparison would not be a weaker test, it would be a different one.
#![allow(clippy::float_cmp)]

use super::*;

/// Sample points in normalised coordinates, scaled by each field's domain.
///
/// Hand-chosen so that within each point the three coordinate magnitudes are
/// distinct and none is near zero. That matters: the box's gradient is
/// non-unique on its interior medial axis, where two components of `q` tie, and
/// the right response is to keep sample points off it **by construction** rather
/// than to filter failures out of a loop.
const SAMPLES: [[f64; 3]; 12] = [
    [0.137, 0.523, -0.311],
    [-0.412, 0.229, 0.671],
    [0.733, -0.158, 0.394],
    [-0.259, -0.617, 0.142],
    [0.481, 0.352, -0.706],
    [-0.694, 0.443, -0.217],
    [0.312, -0.771, 0.556],
    [-0.148, 0.685, 0.423],
    [0.607, -0.294, -0.469],
    [-0.535, -0.371, -0.628],
    [0.226, 0.114, 0.759],
    [-0.783, 0.596, 0.183],
];

/// Two radii: one well inside the compact fields' surfaces, one well outside, so
/// both branches of the box gradient get exercised.
const RADII: [f64; 2] = [0.4, 0.9];

fn scaled_samples<S: ReferenceField<Scalar = f64>>(field: &S) -> impl Iterator<Item = [f64; 3]> {
    let (_, hi) = field.domain();
    let half = hi[0];
    RADII.into_iter().flat_map(move |r| {
        SAMPLES
            .into_iter()
            .map(move |p| [p[0] * half * r, p[1] * half * r, p[2] * half * r])
    })
}

// ─── sweeps over all seven fields ───────────────────────────────────────────

/// The analytic gradient against a central difference, componentwise.
///
/// Componentwise absolute agreement is the right check rather than a comparison
/// of magnitudes: a single flipped component barely moves `|∇f|` but shows up
/// here as a difference of twice the component. That is the actual bug class in
/// a hand-derived gradient.
///
/// Run in `f64`. In `f32` a central difference has a noise floor near
/// `cbrt(ε) ≈ 6e-3`, far too loose to catch a sign error in one component.
#[test]
fn analytic_gradients_match_central_differences() {
    fn check<S: ReferenceField<Scalar = f64>>(name: &str, field: &S) {
        let h = 1e-5f64;
        for p in scaled_samples(field) {
            let analytic = field.gradient(p);
            for axis in 0..3 {
                let mut lo = p;
                let mut hi = p;
                lo[axis] -= h;
                hi[axis] += h;
                let numeric = (field.sample(hi) - field.sample(lo)) / (2.0 * h);
                assert!(
                    (analytic[axis] - numeric).abs() < 1e-6,
                    "{name} at {p:?} axis {axis}: analytic {} vs numeric {numeric}",
                    analytic[axis],
                );
            }
        }
    }

    for_each_reference_field!(f64, |name, field| {
        check(name, &field);
    });
}

/// Where a field claims to be a signed distance, `|∇f|` must actually be one.
///
/// This is the test that catches a missing normalisation, and it is why
/// `is_exact_distance` lives on the trait instead of in a comment.
#[test]
fn exact_distance_fields_have_unit_gradients() {
    fn check<S: ReferenceField<Scalar = f64>>(name: &str, field: &S) {
        if !field.is_exact_distance() {
            return;
        }
        for p in scaled_samples(field) {
            let g = field.gradient(p);
            let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
            assert!(
                (len - 1.0).abs() < 1e-12,
                "{name} at {p:?}: |grad| = {len}, grad = {g:?}"
            );
        }
    }

    for_each_reference_field!(f64, |name, field| {
        check(name, &field);
    });
}

/// Every closed field's surface must lie strictly inside its declared domain.
///
/// Mechanically catches "someone tightened a domain and the gyroid now touches
/// the wall", which would silently invalidate every Euler-characteristic
/// expectation downstream.
#[test]
fn closed_fields_have_constant_sign_on_their_domain_boundary() {
    fn check<S: ReferenceField<Scalar = f64>>(name: &str, field: &S) {
        if !field.closed_in_domain() {
            return;
        }
        let (lo, hi) = field.domain();
        let steps = 8;
        for axis in 0..3 {
            for &wall in &[lo[axis], hi[axis]] {
                for i in 0..=steps {
                    for j in 0..=steps {
                        let a = axis;
                        let b = (axis + 1) % 3;
                        let c = (axis + 2) % 3;
                        let mut p = [0.0f64; 3];
                        p[a] = wall;
                        p[b] = lo[b] + (hi[b] - lo[b]) * f64::from(i) / f64::from(steps);
                        p[c] = lo[c] + (hi[c] - lo[c]) * f64::from(j) / f64::from(steps);
                        let v = field.sample(p);
                        assert!(v > 0.0, "{name}: sample{p:?} = {v}, not outside the solid");
                    }
                }
            }
        }
    }

    for_each_reference_field!(f64, |name, field| {
        check(name, &field);
    });
}

/// Both scalar monomorphisations exist and behave. The `f64` sweeps above would
/// pass even if the `f32` instantiation did not compile.
#[test]
fn every_field_instantiates_in_f32() {
    for_each_reference_field!(f32, |name, field| {
        let (lo, hi) = field.domain();
        assert!(hi[0] > lo[0], "{name}");
        let v = field.sample([0.1f32, 0.2, 0.3]);
        assert!(v.is_finite(), "{name}: {v}");
        let g = field.gradient([0.1f32, 0.2, 0.3]);
        assert!(g.iter().all(|c| c.is_finite()), "{name}: {g:?}");
    });
}

// ─── per-field sign tests, with hand-computed values ────────────────────────

/// The sign-convention guard. Negative is inside; half of all inside-out bugs
/// are this flipping across a module boundary.
#[test]
fn negative_is_inside() {
    for_each_reference_field!(f64, |name, field| {
        let _ = name;
        let _ = &field;
    });
    assert!(Sphere::<f64>::canonical().sample([0.0, 0.0, 0.0]) < 0.0);
    assert!(BoxExact::<f64>::canonical().sample([0.0, 0.0, 0.0]) < 0.0);
}

#[test]
fn sphere_signs() {
    let f = Sphere::<f64>::canonical();
    assert_eq!(f.sample([0.0, 0.0, 0.0]), -1.0); // centre
    assert_eq!(f.sample([2.0, 0.0, 0.0]), 1.0); // outside
    assert_eq!(f.sample([1.0, 0.0, 0.0]), 0.0); // on the surface
}

#[test]
fn torus_signs() {
    let f = Torus::<f64>::canonical();
    // Ring lies in the xz-plane: (major, 0, 0) is the tube's core circle.
    assert_eq!(f.sample([1.0, 0.0, 0.0]), -0.3);
    assert_eq!(f.sample([0.0, 0.0, 1.0]), -0.3); // and on the z axis too
    assert_eq!(f.sample([0.0, 0.0, 0.0]), 0.7); // through the hole
    // `major + minor` is where the surface is, but 1.3 has no exact binary
    // representation, so this point is a few ulps off the surface rather than on
    // it. Asserted to within rounding, not exactly -- the alternative would be
    // asserting a number that is precise about the wrong thing.
    assert!(f.sample([1.3, 0.0, 0.0]).abs() < 8.0 * f64::EPSILON);
    // The axis of revolution is +y, so a point straight up is outside.
    assert!(f.sample([0.0, 1.0, 0.0]) > 0.0);
}

#[test]
fn box_signs() {
    let f = BoxExact::<f64>::canonical();
    assert_eq!(f.sample([0.0, 0.0, 0.0]), -1.0);
    assert_eq!(f.sample([2.0, 0.0, 0.0]), 1.0); // face
    assert_eq!(f.sample([1.0, 0.0, 0.0]), 0.0); // on the face
    assert_eq!(f.sample([2.0, 2.0, 2.0]), 3.0f64.sqrt()); // corner: the exact
    // form, where the `max(q)` bound would wrongly report 1.0.
}

#[test]
fn thin_plate_signs() {
    let f = ThinPlate::<f64>::canonical();
    let h = ThinPlate::<f64>::CANONICAL_CELL_SIZE;
    let half_thickness = h * ThinPlate::<f64>::THICKNESS_IN_CELLS * 0.5;

    assert_eq!(f.sample([0.0, 0.0, 0.0]), -half_thickness);
    assert_eq!(f.sample([0.0, half_thickness, 0.0]), 0.0);
    // One whole cell above the mid-plane is comfortably outside -- which is the
    // property that makes marching cubes miss the plate entirely.
    assert_eq!(f.sample([0.0, h, 0.0]), h - half_thickness);
    assert!(half_thickness * 2.0 < h, "the plate must be sub-voxel");
}

#[test]
fn csg_difference_signs() {
    let f = csg_difference::<f64>();
    // Material: near the -y corner, far from the sphere.
    assert!(f.sample([0.9, -0.9, 0.0]) < 0.0);
    // Scooped away: inside the box but inside the sphere too.
    assert!(f.sample([0.9, 0.9, 0.9]) > 0.0);
    // Outside the box entirely.
    assert_eq!(f.sample([2.0, 0.0, 0.0]), 1.0);
    // The box corner the sphere reaches: |(1,1,1) - (0.6,0.6,0.6)| = 0.6928 < 0.75.
    assert!(f.sample([1.0, 1.0, 1.0]) > 0.0);
}

#[test]
fn gyroid_signs() {
    let f = capped_gyroid::<f64>();
    // g(0) = 0 exactly: every term is sin(0)*cos(0).
    assert_eq!(f.sample([0.0, 0.0, 0.0]), 0.0);
    // g(-pi/2, 0, 0) = sin(-pi/2)*cos(0) = -1.
    let inside = f.sample([-core::f64::consts::FRAC_PI_2, 0.0, 0.0]);
    assert!((inside + 1.0).abs() < 1e-12, "{inside}");
    let outside = f.sample([core::f64::consts::FRAC_PI_2, 0.0, 0.0]);
    assert!((outside - 1.0).abs() < 1e-12, "{outside}");
    // Beyond the spherical cap, the intersection is positive whatever g does.
    assert!(f.sample([6.5, 0.0, 0.0]) > 0.0);

    // The uncapped surface really is not a distance field -- this is the reason
    // `is_exact_distance` exists.
    assert!(!f.is_exact_distance());
}

#[test]
fn fbm_terrain_signs() {
    let f = FbmTerrain::<f64>::canonical();
    let bound = f.height_bound();
    // Proved from the bound, never from a sampled height: below every possible
    // surface is inside, above every possible surface is outside.
    assert!(f.sample([0.0, -(bound + 1.0), 0.0]) < 0.0);
    assert!(f.sample([0.0, bound + 1.0, 0.0]) > 0.0);
    // 2.0 amplitude * 2.0 per-octave bound * (1 + 1/2 + 1/4 + 1/8)
    assert_eq!(bound, 7.5);
    assert!(!f.closed_in_domain());
    assert!(f.expected_euler().is_none());
}

// ─── metadata ───────────────────────────────────────────────────────────────

#[test]
fn reference_field_names_match_the_sweep() {
    let mut count = 0;
    for_each_reference_field!(f64, |name, field| {
        let _ = &field;
        count += 1;
        assert!(!name.is_empty());
    });
    assert_eq!(count, 8);

    assert_eq!(Sphere::<f64>::NAME, "sphere");
    assert_eq!(Torus::<f64>::NAME, "torus");
    assert_eq!(BoxExact::<f64>::NAME, "box_exact");
    assert_eq!(CsgDifference::<f64>::NAME, "csg_difference");
    assert_eq!(ThinPlate::<f64>::NAME, "thin_plate");
    assert_eq!(CappedGyroid::<f64>::NAME, "gyroid");
    assert_eq!(FbmTerrain::<f64>::NAME, "fbm_terrain");
    assert_eq!(<NoiseCavity<f64> as ReferenceField>::NAME, "noise_cavity");
}

/// Two fields have no analytically known Euler characteristic, and saying so is
/// the point -- inventing one would be exactly the kind of guess the project
/// rules forbid.
#[test]
fn only_the_analytically_known_euler_characteristics_are_declared() {
    assert_eq!(Sphere::<f64>::canonical().expected_euler(), Some(2));
    assert_eq!(Torus::<f64>::canonical().expected_euler(), Some(0)); // genus 1
    assert_eq!(BoxExact::<f64>::canonical().expected_euler(), Some(2));
    assert_eq!(ThinPlate::<f64>::canonical().expected_euler(), Some(2));
    assert_eq!(csg_difference::<f64>().expected_euler(), Some(2));
    assert_eq!(capped_gyroid::<f64>().expected_euler(), None);
    assert_eq!(FbmTerrain::<f64>::canonical().expected_euler(), None);
}

/// A field defined by `canonical()` and nothing else is a field that cannot
/// drift, which is the entire reason the constructor exists.
#[test]
fn default_is_canonical() {
    assert_eq!(Sphere::<f32>::default(), Sphere::canonical());
    assert_eq!(Torus::<f32>::default(), Torus::canonical());
    assert_eq!(BoxExact::<f32>::default(), BoxExact::canonical());
    assert_eq!(ThinPlate::<f32>::default(), ThinPlate::canonical());
    assert_eq!(Gyroid::<f32>::default(), Gyroid::canonical());
    assert_eq!(FbmTerrain::<f32>::default(), FbmTerrain::canonical());
}

/// The combinators are generic over any pair of fields, not just the reference
/// ones, and the difference's gradient really does switch operands.
#[test]
fn difference_gradient_follows_the_active_operand() {
    let f = csg_difference::<f64>();
    // Deep in the box, far from the sphere: the box is active.
    let g = f.gradient([-0.9, -0.5, 0.0]);
    assert_eq!(g, BoxExact::<f64>::canonical().gradient([-0.9, -0.5, 0.0]));

    // Inside the scoop: the negated sphere is active.
    let p = [0.85, 0.85, 0.85];
    let g = f.gradient(p);
    let sphere_grad = Sphere::<f64> {
        center: [0.6; 3],
        radius: 0.75,
    }
    .gradient(p);
    assert_eq!(g, [-sphere_grad[0], -sphere_grad[1], -sphere_grad[2]]);
}
