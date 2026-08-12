//! Property-test scaffolding: randomized fields, randomized grids, and the
//! assertion bundle every algorithm ticket runs on a mesh it produced.
//!
//! Test-only. These generators exist so that an extraction algorithm can be
//! attacked with a thousand fields nobody chose by hand, which is the only way
//! to reach the configurations a hand-written fixture never thinks of.
//!
//! # The property the random fields are checked against
//!
//! Randomly generated fields have no known surface, so there is nothing to
//! compare a sample against. What *is* exactly true, and strong enough to catch
//! a broken combinator, is that they are **1-Lipschitz**: `min` and `max` of
//! 1-Lipschitz functions are 1-Lipschitz, and an exact sphere or half-space
//! distance is 1-Lipschitz to begin with. So
//!
//! ```text
//! |f(a) − f(b)| ≤ |a − b|    for all a, b
//! ```
//!
//! holds everywhere, including at the seams where the field is not
//! differentiable and a gradient check would have to be skipped.
//! `a_field_that_is_not_one_lipschitz_is_caught` is the negative control: a
//! property test that cannot fail is decoration.

use alloc::vec;
use alloc::vec::Vec;

use proptest::prelude::*;

use crate::fields::Sphere;
use crate::validate::{
    SelfIntersectionReport, ValidateConfig, self_intersections, validate_features,
};
use crate::vec3;
use crate::{Real, RuntimeShape3, Sdf, Shape3};

pub(crate) mod extraction;

/// Half-extent of the box every generated field lives inside.
pub(crate) const DOMAIN: f64 = 2.0;

// ─── randomized fields ──────────────────────────────────────────────────────

/// The union of several spheres, as `min` over them.
///
/// Not an exact distance field inside the union, but 1-Lipschitz everywhere,
/// which is the property being checked.
#[derive(Clone, Debug)]
pub(crate) struct SphereUnion {
    pub spheres: Vec<Sphere<f64>>,
}

impl Sdf for SphereUnion {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let mut best = f64::INFINITY;
        for s in &self.spheres {
            best = best.min(s.sample(p));
        }
        best
    }

    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let mut best = f64::INFINITY;
        let mut grad = [0.0; 3];
        for s in &self.spheres {
            let v = s.sample(p);
            if v < best {
                best = v;
                grad = s.gradient(p);
            }
        }
        grad
    }
}

/// `f(p) = n·p − offset`, with `n` a unit vector. An exact distance to a plane.
#[derive(Clone, Copy, Debug)]
pub(crate) struct HalfSpace {
    pub normal: [f64; 3],
    pub offset: f64,
}

/// The intersection of several half-spaces with a bounding sphere: a convex body
/// with sharp edges, as `max` over the operands.
///
/// The bound is what keeps it closed — an intersection of half-spaces alone is
/// generally unbounded.
#[derive(Clone, Debug)]
pub(crate) struct ConvexBody {
    pub planes: Vec<HalfSpace>,
    pub bound: Sphere<f64>,
}

impl Sdf for ConvexBody {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let mut worst = self.bound.sample(p);
        for h in &self.planes {
            worst = worst.max(vec3::dot(h.normal, p) - h.offset);
        }
        worst
    }
}

fn sphere() -> impl Strategy<Value = Sphere<f64>> {
    (-1.0f64..=1.0, -1.0f64..=1.0, -1.0f64..=1.0, 0.25f64..=0.9).prop_map(|(x, y, z, radius)| {
        Sphere {
            center: [x, y, z],
            radius,
        }
    })
}

pub(crate) fn sphere_union() -> impl Strategy<Value = SphereUnion> {
    proptest::collection::vec(sphere(), 1..=4).prop_map(|spheres| SphereUnion { spheres })
}

/// Uniform on the sphere, and exactly unit by construction.
///
/// Sampling `cos φ` uniformly rather than `φ` avoids both the polar clustering
/// of naive spherical angles and the near-zero-length guard a normalised random
/// vector would need.
fn unit_vector() -> impl Strategy<Value = [f64; 3]> {
    (0.0f64..core::f64::consts::TAU, -1.0f64..=1.0).prop_map(|(theta, cos_phi)| {
        let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
        [
            sin_phi * Real::cos(theta),
            sin_phi * Real::sin(theta),
            cos_phi,
        ]
    })
}

pub(crate) fn convex_body() -> impl Strategy<Value = ConvexBody> {
    proptest::collection::vec(
        (unit_vector(), -0.8f64..=0.8).prop_map(|(normal, offset)| HalfSpace { normal, offset }),
        1..=6,
    )
    .prop_map(|planes| ConvexBody {
        planes,
        bound: Sphere {
            center: [0.0; 3],
            radius: 1.5,
        },
    })
}

/// A point anywhere in the shared domain.
pub(crate) fn point() -> impl Strategy<Value = [f64; 3]> {
    (-DOMAIN..=DOMAIN, -DOMAIN..=DOMAIN, -DOMAIN..=DOMAIN).prop_map(|(x, y, z)| [x, y, z])
}

/// Grid dimensions, deliberately allowed to be non-cubic — a cubic shape hides
/// every stride bug.
pub(crate) fn resolution() -> impl Strategy<Value = [u32; 3]> {
    (2u32..=20, 2u32..=20, 2u32..=20).prop_map(|(x, y, z)| [x, y, z])
}

// ─── the assertion bundle ───────────────────────────────────────────────────

/// Which validity gate a mesh is held to.
///
/// Deliberately an enum rather than a `bool`, and deliberately three cases: a
/// blanket gate is unsatisfiable for at least one field *and* at least one
/// algorithm, so the caller has to name which one applies and why.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SurfaceGate {
    /// A closed, oriented 2-manifold. **The gate for a closed field on a grid
    /// that resolves it** — the seven reference fields at their own resolutions,
    /// meshed by Marching Cubes.
    ///
    /// Note the condition is on the *pair*, not on the algorithm. Marching Cubes
    /// was believed to earn this unconditionally, by placing vertices on grid
    /// edges rather than one per cell. ✗15 shows it does not: refine far enough
    /// and it does, but a surface that pinches inside one cell defeats it.
    Closed,
    /// A 2-manifold, possibly with boundary. **The gate for an open field**, and
    /// for any single chunk once G-001 lands.
    Manifold,
    /// Closed, correctly oriented and wholly inside the grid, but *permitted* to
    /// be non-manifold.
    ///
    /// **The gate for a grid that may not resolve the field's topology.** Note
    /// what this is *not* keyed on: it started life as "the gate for
    /// one-vertex-per-cell methods", and that was wrong twice over.
    ///
    /// Two distinct mechanisms land here, which is why the name is about the
    /// grid rather than the algorithm:
    ///
    /// - **Surface Nets and plain Dual Contouring** place one vertex per cell, so
    ///   two sheets passing through one cell must share it. The literature calls
    ///   this DC's *"actual structural defect"*; A-010 fixes it architecturally
    ///   by vertex splitting. M-4 measured it on `gyroid` and `fbm_terrain` and
    ///   read it as a high-genus/open-field effect; M-15 corrected that — a
    ///   generated **convex body** does it too, so it is about resolution.
    /// - **Marching Cubes**, which was believed unconditionally manifold, does it
    ///   too where the surface *pinches* inside a single cell: the shared grid
    ///   edge ends up carrying four faces. See ✗15 and
    ///   `an_under_resolved_pinch_makes_marching_cubes_non_manifold`, which pins
    ///   the exact counts at `h = 2/3` and their disappearance by `h = 1/2`.
    ///
    /// The strict [`Closed`](Self::Closed) claim is still tested where it is
    /// actually true — the seven reference fields at their own resolutions, in
    /// `mc/tests.rs`. It is only the *generated* fields, which are adversarial by
    /// construction and go as coarse as `h = 2/3`, that need this.
    ///
    /// What is still asserted is everything unrelated to unresolved topology: no
    /// structural errors, no boundary (the surface did not leave the grid), and
    /// consistent winding.
    ///
    /// **The even-`χ` parity check is deliberately *not* asserted here**, and
    /// that is not an oversight. `χ = 2 − 2g` — hence `χ` even — holds for a
    /// closed *orientable manifold*, so parity is a corollary of manifoldness
    /// rather than an independent check. Waiving manifoldness and keeping the
    /// parity check is incoherent, and measurably so: Surface Nets on a generated
    /// convex body produces `χ = 1` with one non-manifold edge and zero boundary
    /// edges. A-010 is where this becomes assertable again.
    ClosedAllowingUnresolvedTopology,
}

/// Everything an algorithm ticket must assert about a mesh it produced.
///
/// Returns the self-intersection report rather than gating on it: dual methods
/// are *expected* to be non-zero there, and A-009 is the ticket that measures
/// what the cell clamp changes. Callers assert whatever they know about their
/// own algorithm.
///
/// `gate` comes from the field and the algorithm, never from the caller's
/// intuition — see `ReferenceField::closed_in_domain` and [`SurfaceGate`].
pub(crate) fn assert_extracted_mesh_is_valid(
    label: &str,
    positions: &[[f64; 3]],
    indices: &[u32],
    cell_size: f64,
    gate: SurfaceGate,
) -> SelfIntersectionReport {
    let (report, features) = validate_features(
        positions,
        indices,
        &ValidateConfig::from_cell_size(cell_size).expect("valid cell size"),
    );
    assert!(
        !report.has_structural_errors(),
        "{label}: malformed mesh\n{report}"
    );

    // The lists and the counts come from one pass, so they cannot disagree --
    // but E-111 draws the lists beside the counts, so a drift here would put a
    // wrong picture next to a right number. Checked on every case the bundle
    // sees, which is over a thousand generated meshes.
    assert_eq!(
        features.edges.len() as u64,
        report.non_manifold_edges,
        "{label}: non-manifold edge list disagrees with the count\n{report}"
    );
    assert_eq!(
        features.vertices.len() as u64,
        report.non_manifold_vertices,
        "{label}: non-manifold vertex list disagrees with the count\n{report}"
    );
    assert_eq!(
        features.boundary_edges.len() as u64,
        report.boundary_edges,
        "{label}: boundary edge list disagrees with the count\n{report}"
    );
    assert_eq!(
        features.inconsistently_oriented_edges.len() as u64,
        report.inconsistently_oriented_edges,
        "{label}: orientation list disagrees with the count\n{report}"
    );
    match gate {
        SurfaceGate::Closed => assert!(
            report.is_closed(),
            "{label}: expected a closed surface\n{report}"
        ),
        SurfaceGate::Manifold => assert!(
            report.is_manifold(),
            "{label}: expected a manifold with boundary\n{report}"
        ),
        SurfaceGate::ClosedAllowingUnresolvedTopology => {
            assert_eq!(
                report.boundary_edges, 0,
                "{label}: the surface left the grid\n{report}"
            );
            assert_eq!(
                report.inconsistently_oriented_edges, 0,
                "{label}: inconsistent winding\n{report}"
            );
        }
    }
    self_intersections(positions, indices, cell_size).expect("self-intersection scan")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lipschitz_violation<S: Sdf<Scalar = f64>>(field: &S, a: [f64; 3], b: [f64; 3]) -> bool {
        let distance = vec3::length(vec3::sub(a, b));
        let change = (field.sample(a) - field.sample(b)).abs();
        // A relative slack plus an absolute floor, so the check is neither
        // fooled by rounding at large coordinates nor vacuous near zero.
        change > distance * (1.0 + 1e-9) + 1e-9
    }

    proptest! {
        #![proptest_config(ProptestConfig { cases: 1000, ..ProptestConfig::default() })]

        #[test]
        fn sphere_unions_are_one_lipschitz(
            field in sphere_union(),
            a in point(),
            b in point(),
        ) {
            prop_assert!(
                !lipschitz_violation(&field, a, b),
                "|f(a) - f(b)| exceeded |a - b| for {field:?}"
            );
        }

        #[test]
        fn convex_bodies_are_one_lipschitz(
            field in convex_body(),
            a in point(),
            b in point(),
        ) {
            prop_assert!(
                !lipschitz_violation(&field, a, b),
                "|f(a) - f(b)| exceeded |a - b| for {field:?}"
            );
        }

        #[test]
        fn generated_fields_are_finite(field in sphere_union(), p in point()) {
            prop_assert!(field.sample(p).is_finite());
            prop_assert!(field.gradient(p).iter().all(|c| c.is_finite()));
        }

        /// A sphere's own centre is inside the union that contains it, whatever
        /// the other spheres do — `min` can only make the value smaller.
        #[test]
        fn a_sphere_centre_is_inside_its_union(field in sphere_union()) {
            for s in &field.spheres {
                prop_assert!(field.sample(s.center) <= -s.radius);
            }
        }

        /// Grid sizes, deliberately non-cubic. The same round-trip the
        /// hand-written fixture checks at [3, 5, 7], over a thousand shapes.
        #[test]
        fn grid_indices_round_trip(size in resolution()) {
            let shape = RuntimeShape3::new(size).expect("generated sizes fit u32");
            for i in 0..shape.element_count() as u32 {
                prop_assert_eq!(shape.linearize(shape.delinearize(i)), i);
            }
        }
    }

    /// The negative control for the Lipschitz property. A field scaled by two
    /// changes twice as fast as the distance between samples, so the check must
    /// reject it — otherwise the thousand cases above prove nothing.
    #[test]
    fn a_field_that_is_not_one_lipschitz_is_caught() {
        struct Steep(Sphere<f64>);
        impl Sdf for Steep {
            type Scalar = f64;
            fn sample(&self, p: [f64; 3]) -> f64 {
                2.0 * self.0.sample(p)
            }
        }
        let field = Steep(Sphere::canonical());
        assert!(
            lipschitz_violation(&field, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
            "the property test cannot fail, so it is not testing anything"
        );
    }

    /// The negative control for the mesh bundle. A non-manifold mesh must stop
    /// it, or every algorithm ticket that calls it is asserting nothing.
    #[test]
    #[should_panic(expected = "expected a closed surface")]
    fn the_mesh_bundle_rejects_a_non_manifold_mesh() {
        // Three triangles sharing one edge.
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
        ];
        let indices = [0, 1, 2, 0, 1, 3, 0, 1, 4];
        assert_extracted_mesh_is_valid(
            "negative control",
            &positions,
            &indices,
            1.0,
            SurfaceGate::Closed,
        );
    }

    /// And that it rejects malformed input before it reaches the topology.
    #[test]
    #[should_panic(expected = "malformed mesh")]
    fn the_mesh_bundle_rejects_a_bad_index() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0f64]];
        let indices = [0, 1, 9];
        assert_extracted_mesh_is_valid(
            "negative control",
            &positions,
            &indices,
            1.0,
            SurfaceGate::Manifold,
        );
    }

    /// A well-formed closed mesh passes, so the bundle is not merely a
    /// panic generator.
    #[test]
    fn the_mesh_bundle_accepts_a_tetrahedron() {
        let positions = vec![
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ];
        let indices = [0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2];
        let si = assert_extracted_mesh_is_valid(
            "tetrahedron",
            &positions,
            &indices,
            1.0,
            SurfaceGate::Closed,
        );
        assert!(si.is_intersection_free(), "{si}");
    }
}
