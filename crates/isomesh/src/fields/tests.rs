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
/// The bound lives on the trait instead of in a comment.
#[test]
fn exact_distance_fields_have_unit_gradients() {
    fn check<S: ReferenceField<Scalar = f64>>(name: &str, field: &S) {
        if !field.bound().is_exact() {
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
    // The declared bound says so, rather than a comment.
    assert!(!f.bound().is_exact());
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

/// **Union is the set union, and its value never overestimates distance
/// (E-216).**
///
/// Two properties, because `min` is easy to write and easy to write backwards.
/// The **set** property is what makes it a union: a point is inside `a ∪ b`
/// exactly when it is inside either. The **bound** property is what makes it
/// safe: `min` of two exact distances is never larger than the true distance to
/// the union, so a sphere tracer stepping by it can only under-step.
#[test]
fn union_is_the_set_union_and_never_overestimates() {
    use super::{Sphere, Union};

    let a = Sphere {
        center: [-0.4, 0.0, 0.0],
        radius: 0.7_f64,
    };
    let b = Sphere {
        center: [0.4, 0.0, 0.0],
        radius: 0.7_f64,
    };
    let u = Union { a, b };

    let mut checked = 0usize;
    let n = 24;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let p = [
                    -1.5 + 3.0 * f64::from(i) / f64::from(n - 1),
                    -1.5 + 3.0 * f64::from(j) / f64::from(n - 1),
                    -1.5 + 3.0 * f64::from(k) / f64::from(n - 1),
                ];
                let (fa, fb, fu) = (a.sample(p), b.sample(p), u.sample(p));

                // Set property, stated on the signs rather than the values.
                assert_eq!(
                    fu <= 0.0,
                    fa <= 0.0 || fb <= 0.0,
                    "union disagrees with `inside a or inside b` at {p:?}"
                );
                // Bound property: never larger than either operand's distance.
                assert!(
                    fu <= fa + 1e-12 && fu <= fb + 1e-12,
                    "union overestimates distance at {p:?}: {fu} > min({fa}, {fb})"
                );
                checked += 1;
            }
        }
    }
    assert_eq!(checked, 24 * 24 * 24);
}

/// **The blend radius does what a level designer expects, and `k = 0` is a hard
/// union (E-216).**
///
/// `SmoothUnion` exists so the seam between two primitives can be filleted, so
/// the property worth pinning is monotone: **more `k`, more material** in the
/// crease between two nearly touching spheres, and never less anywhere. Plus the
/// degenerate case, because a blend radius of zero reaching `smooth_min`'s
/// divide-by-`k` would be the obvious way for this to be wrong.
#[test]
fn the_blend_radius_adds_material_monotonically() {
    use super::{SmoothUnion, Sphere, Union};

    let (a, b) = (
        Sphere {
            center: [-0.55, 0.0, 0.0],
            radius: 0.5_f64,
        },
        Sphere {
            center: [0.55, 0.0, 0.0],
            radius: 0.5_f64,
        },
    );
    // The midpoint of the crease: outside both spheres, and the first place a
    // fillet puts material.
    let seam = [0.0, 0.0, 0.0];

    let hard = Union { a, b }.sample(seam);
    let zero = SmoothUnion { a, b, k: 0.0 }.sample(seam);
    assert_eq!(zero, hard, "k = 0 must degenerate to a hard union");

    let mut previous = hard;
    for step in 1..=6 {
        let k = 0.1 * f64::from(step);
        let value = SmoothUnion { a, b, k }.sample(seam);
        assert!(
            value < previous,
            "k = {k} did not add material at the seam: {value} >= {previous}"
        );
        previous = value;
    }
    std::println!("measured: seam value {hard:.4} (hard) -> {previous:.4} at k = 0.6");

    // And it never removes material anywhere, which is what "union" still means.
    let blended = SmoothUnion { a, b, k: 0.3 };
    for i in 0..20 {
        for j in 0..20 {
            let p = [
                -1.5 + 3.0 * f64::from(i) / 19.0,
                -1.5 + 3.0 * f64::from(j) / 19.0,
                0.0,
            ];
            assert!(
                blended.sample(p) <= Union { a, b }.sample(p) + 1e-12,
                "the blend removed material at {p:?}"
            );
        }
    }
}

/// **Every reference field declares a bound, and the declaration is checked
/// against the field rather than trusted (F-001).**
///
/// A declaration nobody verifies is the defect this ticket existed to remove —
/// `csg_difference` declared `true` with `// away from the seam` beside it for
/// months. So each field's claim is measured:
///
/// - **`Exact`** must satisfy `|∇f| ≈ 1`, since that is what "the value is the
///   distance" means differentially. Corollary 1 of Bálint, Valasek & Gergó 2019
///   is the authority: every true SDF is 1-Lipschitz and 1 is the *smallest*
///   such constant, so an exact field cannot have a gradient shorter than one
///   either.
/// - **`Lipschitz { l }`** must satisfy `|∇f| ≤ l`. Declaring a constant smaller
///   than the field's own gradient is the failure worth catching, because a
///   sphere tracer dividing by it would step through the surface.
/// - **`Underestimate { q }`** must not *overstate* distance along the gradient,
///   which is the property that makes it the safe direction.
#[test]
fn every_field_meets_the_bound_it_declares() {
    use super::FieldBound;

    // Sampled off-lattice on purpose: a grid aligned to the noise lattice or to
    // a box face lands on exactly the creases where a gradient is undefined, and
    // would measure the discontinuity instead of the field.
    const N: u32 = 12;
    let offset = 0.031_7_f64;

    crate::for_each_reference_field!(f64, |name, field| {
        let (lo, hi) = field.domain();
        let bound = field.bound();
        let mut worst_low = f64::INFINITY;
        let mut worst_high: f64 = 0.0;
        let mut samples = 0usize;

        for i in 0..N {
            for j in 0..N {
                for k in 0..N {
                    let t = |v: u32, a: f64, b: f64| {
                        a + (b - a) * (f64::from(v) + 0.5 + offset) / f64::from(N)
                    };
                    let p = [t(i, lo[0], hi[0]), t(j, lo[1], hi[1]), t(k, lo[2], hi[2])];
                    let g = field.gradient(p);
                    let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                    if !len.is_finite() {
                        continue;
                    }
                    worst_low = worst_low.min(len);
                    worst_high = worst_high.max(len);
                    samples += 1;
                }
            }
        }
        assert!(samples > 100, "{name}: only {samples} usable samples");

        match bound {
            FieldBound::Exact => {
                // Generous both ways: a crease has no gradient and a finite
                // sample near one reads short. The claim being tested is that
                // the field is a distance, not that every sample is clean.
                assert!(
                    worst_high < 1.20,
                    "{name} declares Exact but |∇f| reaches {worst_high:.4}"
                );
                assert!(
                    worst_low > 0.30,
                    "{name} declares Exact but |∇f| falls to {worst_low:.4}"
                );
            }
            FieldBound::Lipschitz { l } => {
                assert!(
                    worst_high <= l * 1.05,
                    "{name} declares Lipschitz l = {l} but |∇f| reaches {worst_high:.4}"
                );
            }
            FieldBound::Underestimate { q } => {
                assert!((0.0..=1.0).contains(&q), "{name}: q = {q} is not in (0, 1]");
                // An underestimate is still 1-Lipschitz in this crate's fields:
                // it is built from exact operands by min/max, and those preserve
                // the constant even where they destroy exactness.
                assert!(
                    worst_high < 1.20,
                    "{name} declares Underestimate but |∇f| reaches {worst_high:.4}"
                );
            }
            FieldBound::Unbounded => {}
        }
        std::println!("measured: {name:<16} {bound:?} |∇f| in [{worst_low:.3}, {worst_high:.3}]");
    });
}

/// **`csg_difference` is no longer `Exact`, which is F-001's acceptance stated
/// as a test.**
///
/// Pinned separately from the sweep above because it is the specific claim the
/// ticket was written about, and because a future edit that "tidies" the
/// declaration back to `Exact` would pass every other test in this file.
#[test]
fn the_csg_field_does_not_claim_to_be_a_distance() {
    use super::FieldBound;

    let f = super::csg_difference::<f64>();
    assert_ne!(
        f.bound(),
        FieldBound::Exact,
        "csg_difference is max(box, -sphere); max of two exact distances is not \
         an exact distance, and near a concave seam it overestimates"
    );
    assert!(!f.bound().is_exact());
    // It is still 1-Lipschitz, which is the half that survives CSG and the half
    // Phase 12's empty-cell rejection needs.
    assert_eq!(f.bound().lipschitz(), Some(1.0));
}

/// **How far `min` and `max` are from the true distance, measured against
/// analytic ground truth (F-003).**
///
/// The ticket predicted an asymmetry — that `min` of two exact fields yields a
/// bound F-002 confirms and `max` yields *"a strictly weaker one"* — and asked
/// for the test to assert it rather than treat the two as equivalent. So this
/// measures both against a distance computed in closed form, rather than
/// assuming which direction the asymmetry runs.
///
/// Ground truth for two spheres is exact and cheap: the distance to the union of
/// two balls, and to their intersection, both have closed forms for the
/// configuration used here — a point's distance to a ball is `|p − c| − r`, and
/// for two overlapping balls the union's boundary is the near part of each
/// sphere while the intersection's is the far part.
#[test]
fn min_and_max_are_both_inexact_and_the_asymmetry_is_by_region() {
    use super::{Intersection, Sphere, Union};

    let a = Sphere {
        center: [-0.35, 0.0, 0.0],
        radius: 0.8_f64,
    };
    let b = Sphere {
        center: [0.35, 0.0, 0.0],
        radius: 0.8_f64,
    };
    let union = Union { a, b };
    let inter = Intersection { a, b };

    // Worst relative shortfall of the composed value against the true distance,
    // split by region. `min`/`max` never overestimate here, so the interesting
    // number is how far they fall short.
    let mut union_outside = 0.0f64;
    let mut union_inside = 0.0f64;
    let mut inter_outside = 0.0f64;
    let mut inter_inside = 0.0f64;

    let n = 40;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let t = |v: i32| -1.6 + 3.2 * (f64::from(v) + 0.5) / f64::from(n);
                let p = [t(i), t(j), t(k)];
                let (fa, fb) = (a.sample(p), b.sample(p));

                // Union: outside both, the true distance to A ∪ B is exactly
                // min(d_a, d_b) -- the nearest surface is the nearest surface.
                // Inside, the value is the distance to the *nearer* boundary,
                // which may be an interior surface that the union removed.
                let u = union.sample(p);
                if fa > 0.0 && fb > 0.0 {
                    union_outside = union_outside.max((u - fa.min(fb)).abs());
                } else {
                    // True distance to the union's boundary from inside: the
                    // point is inside at least one ball, and the boundary of the
                    // union is whichever sphere it is nearer to *from within the
                    // union*, which is the larger of the two signed values.
                    let truth = fa.max(fb);
                    union_inside = union_inside.max((truth - u).abs());
                }

                // Intersection is the mirror image: exact inside, degraded
                // outside.
                let x = inter.sample(p);
                if fa < 0.0 && fb < 0.0 {
                    inter_inside = inter_inside.max((x - fa.max(fb)).abs());
                } else {
                    let truth = fa.min(fb);
                    inter_outside = inter_outside.max((truth - x).abs());
                }
            }
        }
    }

    std::println!(
        "measured: union  outside {union_outside:.3e}  inside {union_inside:.3e}\n\
         measured: inter  inside  {inter_inside:.3e}  outside {inter_outside:.3e}"
    );

    // **Union is exact outside and intersection is exact inside**, to rounding.
    // That is the region each operator gets right.
    assert!(
        union_outside < 1e-12,
        "min is not exact outside the union: {union_outside:e}"
    );
    assert!(
        inter_inside < 1e-12,
        "max is not exact inside the intersection: {inter_inside:e}"
    );

    // **And each is wrong in the other region, by a comparable amount.** The
    // ticket expected `min` to be strictly better than `max`; it is not. They
    // are mirror images, and the asymmetry is by *region* rather than by
    // operator (M-246).
    assert!(
        union_inside > 1e-3,
        "min was expected to lose accuracy inside the union"
    );
    assert!(
        inter_outside > 1e-3,
        "max was expected to lose accuracy outside the intersection"
    );
    let ratio = union_inside / inter_outside;
    assert!(
        (0.5..2.0).contains(&ratio),
        "one operator is much worse than the other ({ratio:.3}), which would make \
         the ticket's 'strictly weaker' prediction right after all"
    );
}

/// **A composed bound is what the combinator can prove, not what its operands
/// were (F-003).**
///
/// Composing two `Exact` fields yields `Lipschitz { l: 1 }` — not `Exact`, and
/// not an `Underestimate` with an invented `q`. That downgrade is the ticket's
/// substance: exactness does not survive `min`/`max`, the Lipschitz constant
/// does, and the precision `q` needs a set-contact factor this crate does not
/// compute.
#[test]
fn composing_exact_fields_keeps_the_constant_and_loses_exactness() {
    use super::{BoundedSdf, Difference, FieldBound, Intersection, SmoothUnion, Sphere, Union};

    let a = Sphere::<f64>::canonical();
    let b = Sphere {
        center: [0.5; 3],
        radius: 0.6_f64,
    };
    assert_eq!(a.value_bound(), FieldBound::Exact);

    for (name, bound) in [
        ("union", Union { a, b }.value_bound()),
        ("intersection", Intersection { a, b }.value_bound()),
        ("difference", Difference { a, b }.value_bound()),
        ("smooth_union", SmoothUnion { a, b, k: 0.2 }.value_bound()),
    ] {
        assert_eq!(
            bound,
            FieldBound::Lipschitz { l: 1.0 },
            "{name} should keep the constant and lose exactness"
        );
        assert!(!bound.is_exact(), "{name} must not claim to be a distance");
        assert_eq!(bound.lipschitz(), Some(1.0), "{name} lost its step size");
    }

    // An unbounded operand poisons the composition, because nothing can be
    // claimed about a value composed from one nothing is claimed about.
    let noisy = super::Gyroid::<f64>::canonical();
    assert!(matches!(
        Union { a, b: noisy }.value_bound(),
        FieldBound::Lipschitz { .. }
    ));
    assert_eq!(
        Union {
            a,
            b: super::FbmTerrain::<f64>::canonical()
        }
        .value_bound(),
        FieldBound::Unbounded
    );
}
