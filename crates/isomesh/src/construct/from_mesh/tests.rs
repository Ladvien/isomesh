//! S-006's acceptance: *"round-trip — mesh a sphere, convert back to a field,
//! re-mesh, and compare against the original."*

extern crate std;

use crate::fields::Sphere;
use crate::marching_cubes::MarchingCubes;
use crate::mesh::MeshBuffer;
use crate::validate;
use crate::{RuntimeShape3, Sdf, Shape3};

use super::signed_distance_from_mesh;

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

/// **The round trip (M-259).** Sphere → mesh → field → mesh, compared against
/// the original at every stage.
///
/// This is the end-to-end test the crate did not have. It exercises the
/// extractor, the pseudonormal sign, the closest-point routine and the grid
/// binning against each other, and a fault in any of them shows up as geometry
/// that moved.
#[test]
fn a_sphere_survives_the_round_trip() {
    let field = Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new([33; 3]).expect("valid shape");
    let h = 0.125_f64;
    let origin = [-2.0; 3];

    let analytic = sample_grid(&field, &shape, origin, h);
    let first = mesh(&analytic, &shape, origin, h);
    assert!(!first.indices.is_empty(), "nothing meshed");

    let rebuilt = signed_distance_from_mesh(&first.positions, &first.indices, &shape, origin, h)
        .expect("mesh to field");
    let second = mesh(&rebuilt, &shape, origin, h);

    // **The sign must agree with the analytic field at every sample.** This is
    // Theorem 1's whole claim, and it is the assertion that would catch a
    // pseudonormal taken from the wrong feature -- a face normal where a vertex
    // one was needed is right for most samples and wrong near a crease.
    let mut sign_disagreements = 0usize;
    for (i, (&a, &b)) in analytic.iter().zip(&rebuilt).enumerate() {
        // Samples within half a cell of the surface are excluded: the meshed
        // surface is a chord of the analytic one there, so the two genuinely
        // disagree about which side a sample is on, and that is discretisation
        // rather than a wrong sign.
        if a.abs() > 0.5 * h && (a < 0.0) != (b < 0.0) {
            sign_disagreements += 1;
            if sign_disagreements < 4 {
                std::println!("  sample {i}: analytic {a:.5}, rebuilt {b:.5}");
            }
        }
    }
    assert_eq!(
        sign_disagreements, 0,
        "the pseudonormal sign disagrees with the analytic field"
    );

    // Geometry: the re-meshed surface must sit on the original.
    let mut worst = 0.0f64;
    for p in &second.positions {
        worst = worst.max(field.sample(*p).abs());
    }
    let mut worst_first = 0.0f64;
    for p in &first.positions {
        worst_first = worst_first.max(field.sample(*p).abs());
    }

    let cfg = validate::ValidateConfig::from_cell_size(h).expect("valid cell size");
    let a = validate::validate(&first, &cfg);
    let b = validate::validate(&second, &cfg);
    std::println!(
        "measured: round trip — {} → {} vertices, {} → {} triangles; \
         worst |field| on vertices {worst_first:.6} → {worst:.6}",
        first.positions.len(),
        second.positions.len(),
        first.indices.len() / 3,
        second.indices.len() / 3,
    );
    std::println!(
        "measured: χ {} → {}, boundary edges {} → {}, non-manifold edges {} → {}",
        a.euler_characteristic,
        b.euler_characteristic,
        a.boundary_edges,
        b.boundary_edges,
        a.non_manifold_edges,
        b.non_manifold_edges
    );

    // The round trip must not degrade the topology. A closed sphere stays a
    // closed sphere; anything else means the field picked up a sign error that
    // survived the excursion above.
    assert!(b.is_closed(), "the round trip opened the surface");
    assert_eq!(
        a.euler_characteristic, b.euler_characteristic,
        "the round trip changed the genus"
    );

    // **Half a cell.** The round trip resamples a discretised surface onto the
    // same grid, so a vertex can move by at most the difference between the two
    // fields at that grid, which is bounded by the sampling.
    assert!(
        worst < 0.5 * h,
        "re-meshed vertices sit {worst} from the analytic sphere, which is {} cells",
        worst / h
    );
}

/// A single triangle's pseudonormals collapse to the paper's stated cases.
///
/// The face case must equal the face normal exactly, and an edge with one
/// incident face must equal it too — both are what `Σ αᵢ nᵢ` reduces to when
/// there is one term, and getting either wrong would be invisible on a closed
/// mesh's interior while wrong at every boundary.
#[test]
fn pseudonormals_collapse_to_the_face_normal_on_a_lone_triangle() {
    let positions = std::vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let indices = std::vec![0u32, 1, 2];
    let n = super::Pseudonormals::build(&positions, &indices);

    // Exact equality on purpose: the cross product of two axis-aligned unit
    // vectors is +z with no rounding at all, and a tolerance here would accept a
    // normalisation that quietly went wrong.
    #[allow(clippy::float_cmp)]
    {
        assert!(n.face[0] == [0.0, 0.0, 1.0]);
        for edge in 0..3u8 {
            assert!(n.at(&indices, 0, super::Feature::Edge(edge)) == [0.0, 0.0, 1.0]);
        }
    }

    // The vertex sums are the face normal scaled by that corner's angle, so
    // they point the same way and have length equal to the angle. π/2 at the
    // right-angled corner, π/4 at each of the others -- and they sum to π,
    // which is the triangle's angle sum and a check that the angles are the
    // interior ones.
    let total: f64 = (0..3).map(|i| n.vertex[i][2]).sum();
    assert!(
        (total - core::f64::consts::PI).abs() < 1e-12,
        "interior angles sum to {total}, not π"
    );
    assert!((n.vertex[0][2] - core::f64::consts::FRAC_PI_2).abs() < 1e-12);
}

/// A closed cube: the sign is correct at the eight places a face normal is not.
///
/// The paper's motivating case. A point diagonally outside a cube corner is
/// closest to that **vertex**, and *"may be outside the mesh but behind a face of
/// the mesh"* — so any single face normal gives the wrong sign there while the
/// angle-weighted sum gives the right one.
#[test]
fn the_sign_is_right_diagonally_outside_a_cube_corner() {
    // Unit cube from [0,0,0] to [1,1,1], outward-oriented.
    let positions = std::vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    let indices: std::vec::Vec<u32> = std::vec![
        0, 2, 1, 0, 3, 2, // -z
        4, 5, 6, 4, 6, 7, // +z
        0, 1, 5, 0, 5, 4, // -y
        3, 7, 6, 3, 6, 2, // +y
        0, 4, 7, 0, 7, 3, // -x
        1, 2, 6, 1, 6, 5, // +x
    ];

    let shape = RuntimeShape3::new([9; 3]).expect("valid shape");
    let h = 0.25_f64;
    let origin = [-0.5; 3];
    let out =
        signed_distance_from_mesh(&positions, &indices, &shape, origin, h).expect("mesh to field");

    let size = shape.size();
    let at = |x: u32, y: u32, z: u32| out[((z * size[1] + y) * size[0] + x) as usize];

    // The eight corners of the sample grid are diagonally outside the cube's
    // own corners, which is exactly the configuration the paper singles out.
    for (x, y, z) in [
        (0, 0, 0),
        (8, 0, 0),
        (0, 8, 0),
        (8, 8, 0),
        (0, 0, 8),
        (8, 0, 8),
        (0, 8, 8),
        (8, 8, 8),
    ] {
        let v = at(x, y, z);
        assert!(v > 0.0, "sample ({x},{y},{z}) is outside but read {v}");
        // Distance from (-0.5,-0.5,-0.5) to (0,0,0) is √3/2.
        assert!(
            (v - (3.0f64).sqrt() / 2.0).abs() < 1e-12,
            "corner distance {v}"
        );
    }

    // The centre is inside, at distance ½ from the nearest face.
    assert!(
        (at(4, 4, 4) + 0.5).abs() < 1e-12,
        "centre read {}",
        at(4, 4, 4)
    );
}

/// Malformed input is refused rather than sampled into nonsense.
#[test]
fn it_refuses_what_it_cannot_mesh() {
    let positions = std::vec![[0.0f64; 3]; 3];
    let shape = RuntimeShape3::new([4; 3]).expect("valid shape");

    // Not a multiple of three.
    assert!(signed_distance_from_mesh(&positions, &[0, 1], &shape, [0.0; 3], 0.1).is_err());
    // Names a vertex that does not exist.
    assert!(signed_distance_from_mesh(&positions, &[0, 1, 9], &shape, [0.0; 3], 0.1).is_err());
    // Grid too small.
    let tiny = RuntimeShape3::new([1, 4, 4]).expect("valid shape");
    assert!(signed_distance_from_mesh(&positions, &[0, 1, 2], &tiny, [0.0; 3], 0.1).is_err());
}

/// **The box reject is worth what the doc claims (M-260).**
///
/// The module header claims the two-level bounding-box reject beats an
/// unaccelerated scan. That is a performance claim, so the benchmark that
/// produced it lives in the repo — rule 4.
///
/// The brute-force comparison here is a **test-local reference**, not a second
/// path in the library: it exists to prove the accelerated one returns
/// *bit-identical* answers, which is the half that makes the speedup mean
/// anything. An earlier uniform-grid version was 3.9× **slower** and this test
/// is what said so.
#[test]
fn the_reject_agrees_with_brute_force_and_is_faster() {
    use std::time::Instant;

    let field = Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new([21; 3]).expect("valid shape");
    let h = 0.2_f64;
    let origin = [-2.0; 3];

    let analytic = sample_grid(&field, &shape, origin, h);
    let m = mesh(&analytic, &shape, origin, h);
    let triangles = m.indices.len() / 3;

    let start = Instant::now();
    let fast = signed_distance_from_mesh(&m.positions, &m.indices, &shape, origin, h)
        .expect("mesh to field");
    let fast_ms = start.elapsed().as_secs_f64() * 1e3;

    // Brute force: every sample against every triangle, same sign rule.
    let start = Instant::now();
    let normals = super::Pseudonormals::build(&m.positions, &m.indices);
    let size = shape.size();
    let mut slow = std::vec::Vec::with_capacity(shape.element_count());
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let p = [
                    origin[0] + h * f64::from(x),
                    origin[1] + h * f64::from(y),
                    origin[2] + h * f64::from(z),
                ];
                let mut best = f64::INFINITY;
                let mut sign = 1.0;
                for t in 0..triangles {
                    let tri = &m.indices[t * 3..t * 3 + 3];
                    let (c, feature) = super::closest_on_triangle(
                        p,
                        m.positions[tri[0] as usize],
                        m.positions[tri[1] as usize],
                        m.positions[tri[2] as usize],
                    );
                    let r = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
                    let d = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
                    if d < best {
                        best = d;
                        let n = normals.at(&m.indices, t, feature);
                        sign = if n[0] * r[0] + n[1] * r[1] + n[2] * r[2] < 0.0 {
                            -1.0
                        } else {
                            1.0
                        };
                    }
                }
                slow.push(best * sign);
            }
        }
    }
    let brute_ms = start.elapsed().as_secs_f64() * 1e3;

    // **Bit-identical, not merely close.** Both walk the same triangles with the
    // same routine; the binning only changes the *order*, and a `<` comparison
    // is order-independent except across exact ties, where both keep the first
    // seen. A tolerance here would hide a real divergence.
    #[allow(clippy::float_cmp)]
    {
        for (i, (a, b)) in fast.iter().zip(&slow).enumerate() {
            assert!(*a == *b, "sample {i}: reject {a} vs brute force {b}");
        }
    }

    std::println!(
        "measured: {} samples × {triangles} triangles — box reject {fast_ms:.1} ms, \
         brute force {brute_ms:.1} ms ({:.1}×)",
        shape.element_count(),
        brute_ms / fast_ms
    );
    assert!(
        fast_ms < brute_ms,
        "the box reject lost to brute force: {fast_ms} vs {brute_ms} ms"
    );
}
