//! S-007's acceptance: *"classifies correctly on a deliberately hole-punched
//! mesh where S-006 fails."*

extern crate std;

use crate::fields::Sphere;
use crate::marching_cubes::MarchingCubes;
use crate::mesh::MeshBuffer;
use crate::{RuntimeShape3, Sdf, Shape3};

use super::{signed_distance_from_mesh_winding, winding_numbers};
use crate::construct::from_mesh::signed_distance_from_mesh;

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

/// A closed sphere: the winding number is 1 inside and 0 outside.
///
/// The calibration. With no boundary edges the cone is empty, the correction
/// term is exactly zero, and this reduces to plain ray-parity — so a failure
/// here is in the intersection routine and not in the construction.
#[test]
fn a_closed_sphere_winds_once_inside_and_not_at_all_outside() {
    let field = Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new([17; 3]).expect("valid shape");
    let h = 0.25_f64;
    let origin = [-2.0; 3];

    let analytic = sample_grid(&field, &shape, origin, h);
    let m = mesh(&analytic, &shape, origin, h);
    let w = winding_numbers(&m.positions, &m.indices, &shape, origin, h).expect("winding");

    let size = shape.size();
    let mut wrong = 0usize;
    let mut worst_inside = 1.0f64;
    let mut worst_outside = 0.0f64;
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let i = ((z * size[1] + y) * size[0] + x) as usize;
                let truth = analytic[i];
                // Skip the band where the meshed surface and the analytic one
                // genuinely disagree about which side a sample is on.
                if truth.abs() <= h {
                    continue;
                }
                if truth < 0.0 {
                    worst_inside = worst_inside.min(w[i]);
                    if w[i] < 0.5 {
                        wrong += 1;
                    }
                } else {
                    worst_outside = worst_outside.max(w[i]);
                    if w[i] > 0.5 {
                        wrong += 1;
                    }
                }
            }
        }
    }
    std::println!(
        "measured: closed sphere — least winding inside {worst_inside:.9}, \
         greatest outside {worst_outside:.9}"
    );
    assert_eq!(wrong, 0, "the closed mesh misclassified {wrong} samples");
}

/// **The acceptance (M-262).** Hole-punched spheres, where the pseudonormal
/// fails and the winding number does not.
///
/// Holes are cut by removing every triangle whose centroid has `x > cut`, which
/// takes a cap off the sphere and leaves a genuine open boundary. **The cut is
/// swept** rather than fixed at one value: a single small hole separates the two
/// methods by only five samples out of 4,491, which is not evidence of anything.
/// Four hole sizes turn that into a trend — the gap widens from 5-versus-0 to
/// 1,435-versus-88 as the cap grows.
///
/// **The comparison is the ticket.** S-006's sign is Bærentzen & Aanæs's
/// theorem, which is about *closed* meshes — on an open one it is not
/// inaccurate, it is answering a different question, and the point of this test
/// is to show that with numbers rather than assert it in prose.
#[test]
fn the_winding_number_survives_holes_that_defeat_the_pseudonormal() {
    let field = Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new([17; 3]).expect("valid shape");
    let h = 0.25_f64;
    let origin = [-2.0; 3];

    let analytic = sample_grid(&field, &shape, origin, h);
    let whole = mesh(&analytic, &shape, origin, h);
    let total = whole.indices.len() / 3;

    std::println!(
        "{:>6} {:>8} {:>9} {:>13} {:>13} {:>12}",
        "cut",
        "removed",
        "boundary",
        "pseudonormal",
        "winding",
        "mean |w-½|"
    );

    let mut worst_normal = 0usize;
    let mut first = true;
    for cut in [0.6_f64, 0.3, 0.0, -0.5] {
        let mut indices = std::vec::Vec::new();
        for tri in whole.indices.as_chunks::<3>().0 {
            let cx = (whole.positions[tri[0] as usize][0]
                + whole.positions[tri[1] as usize][0]
                + whole.positions[tri[2] as usize][0])
                / 3.0;
            if cx <= cut {
                indices.extend_from_slice(tri);
            }
        }
        let removed = total - indices.len() / 3;
        let boundary = super::boundary_edges(&indices);
        assert!(removed > 0, "cut {cut} removed no triangles");
        assert!(!boundary.is_empty(), "cut {cut} left no boundary");

        let by_normal = signed_distance_from_mesh(&whole.positions, &indices, &shape, origin, h)
            .expect("pseudonormal");
        let by_winding =
            signed_distance_from_mesh_winding(&whole.positions, &indices, &shape, origin, h, 0.5)
                .expect("winding");

        let raw = winding_numbers(&whole.positions, &indices, &shape, origin, h).expect("winding");

        let size = shape.size();
        let mut normal_wrong = 0usize;
        let mut winding_wrong = 0usize;
        let mut counted = 0usize;
        let mut ambiguity = 0.0f64;
        for z in 0..size[2] {
            for y in 0..size[1] {
                for x in 0..size[0] {
                    let i = ((z * size[1] + y) * size[0] + x) as usize;
                    let truth = analytic[i];
                    // Skip the band where the meshed surface and the analytic
                    // one genuinely disagree about which side a sample is on.
                    if truth.abs() <= h {
                        continue;
                    }
                    counted += 1;
                    if (truth < 0.0) != (by_normal[i] < 0.0) {
                        normal_wrong += 1;
                    }
                    if (truth < 0.0) != (by_winding[i] < 0.0) {
                        winding_wrong += 1;
                        // How far from the ½ threshold the winding number was.
                        // Near zero means the *measure* was ambiguous rather
                        // than wrong -- which is what a GWN is for.
                        ambiguity += (raw[i] - 0.5).abs();
                    }
                }
            }
        }

        let mean_ambiguity = if winding_wrong == 0 {
            0.0
        } else {
            ambiguity / winding_wrong as f64
        };
        std::println!(
            "{cut:>6.1} {removed:>8} {:>9} {normal_wrong:>7}/{counted:<5} \
             {winding_wrong:>7}/{counted:<5} {mean_ambiguity:>12.4}",
            boundary.len()
        );
        worst_normal = worst_normal.max(normal_wrong);

        // **A small hole must be classified perfectly.** The correction term is
        // exact, so with the geometry still mostly enclosing a point there is no
        // room for the answer to be anything else.
        if first {
            assert_eq!(
                winding_wrong, 0,
                "cut {cut}: the winding number misclassified {winding_wrong} \
                 samples on the smallest hole, where it has no excuse"
            );
            first = false;
        }

        // **At every hole size it must beat the pseudonormal outright.** Not
        // "be perfect", and the `mean |w-½|` column is why: it *grows* with the
        // hole, from 0.0 to 0.27. That is not the winding number becoming
        // confidently wrong -- it is the question becoming wrong. Once half the
        // sphere is gone the mesh no longer encloses the points under the hole,
        // so calling them outside is the correct answer for the surface that
        // actually exists, and scoring against the analytic sphere is scoring
        // against geometry that was deleted. The papers are explicit that a GWN
        // measures *"how confident we can be that a point is inside"*, and the
        // column is that confidence at exactly the samples this counts as
        // wrong.
        assert!(
            winding_wrong * 3 < normal_wrong.max(1),
            "cut {cut}: winding {winding_wrong} against pseudonormal \
             {normal_wrong} -- not the decisive margin the ticket claims"
        );
    }

    // And the pseudonormal must actually fail, or the test proves nothing about
    // the winding number.
    assert!(
        worst_normal > 0,
        "the pseudonormal classified every holed mesh correctly, so this test \
         asserts nothing -- cut a bigger hole"
    );
}

/// A single triangle's solid angle is a hemisphere's worth when it should be.
///
/// An equilateral triangle seen from directly above its centroid at a distance
/// tending to zero subtends `2π`; the check here is the weaker but exact one
/// that a triangle and its reversal subtend equal and opposite angles, which is
/// what the orientation sign has to mean.
#[test]
fn solid_angle_flips_sign_with_orientation() {
    let q = [0.0f64, 0.0, 0.0];
    let a = [1.0f64, 0.0, 1.0];
    let b = [0.0f64, 1.0, 1.0];
    let c = [-1.0f64, -1.0, 1.0];
    let forward = super::solid_angle(q, a, b, c);
    let reverse = super::solid_angle(q, a, c, b);
    assert!(forward.abs() > 0.1, "degenerate test triangle: {forward}");
    assert!(
        (forward + reverse).abs() < 1e-12,
        "{forward} and {reverse} do not cancel"
    );
}

/// The boundary of a closed mesh is empty, and of an open one is not.
#[test]
fn boundary_edges_are_the_directed_ones_with_no_partner() {
    // Two triangles sharing an edge: four boundary edges, not six.
    let closed_pair = std::vec![0u32, 1, 2, 2, 1, 3];
    assert_eq!(super::boundary_edges(&closed_pair).len(), 4);

    // A lone triangle: three.
    assert_eq!(super::boundary_edges(&[0u32, 1, 2]).len(), 3);

    // A tetrahedron, consistently oriented outward: none.
    let tet = std::vec![0u32, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
    assert!(super::boundary_edges(&tet).is_empty());
}

/// T-025. The two functions compute one magnitude, not two.
///
/// On a **closed** mesh the pseudonormal sign is valid and the winding sign
/// agrees with it, so the two functions must agree outright — every sample, bit
/// for bit. That is a stronger statement than "the magnitudes match" and it is
/// the one that fails if a second distance implementation ever reappears here.
// Bit-identity is the claim, not approximate agreement: one implementation
// produces both numbers, so a tolerance would test nothing.
#[allow(clippy::float_cmp)]
#[test]
fn the_magnitude_is_the_pseudonormal_paths_magnitude() {
    use crate::construct::from_mesh::signed_distance_from_mesh;

    // A closed, consistently oriented tetrahedron -- the smallest mesh on which
    // both signs are meaningful.
    let positions = [
        [0.0_f64, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let indices = [0u32, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
    let shape = RuntimeShape3::new([9; 3]).expect("valid shape");
    let origin = [-0.35_f64, -0.35, -0.35];
    let h = 0.2_f64;

    let pseudonormal =
        signed_distance_from_mesh(&positions, &indices, &shape, origin, h).expect("grid is valid");
    let winding = signed_distance_from_mesh_winding(&positions, &indices, &shape, origin, h, 0.5)
        .expect("grid is valid");

    assert_eq!(pseudonormal.len(), winding.len());
    for (i, (&a, &b)) in pseudonormal.iter().zip(winding.iter()).enumerate() {
        assert_eq!(
            a.abs(),
            b.abs(),
            "sample {i}: magnitudes must come from one implementation"
        );
        assert_eq!(
            a, b,
            "sample {i}: on a closed mesh both signs are valid and must agree"
        );
    }
    // The fixture is only meaningful if it straddles the surface.
    assert!(winding.iter().any(|&v| v < 0.0), "some sample is inside");
    assert!(winding.iter().any(|&v| v > 0.0), "some sample is outside");
}
