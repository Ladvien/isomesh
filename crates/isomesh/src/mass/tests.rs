//! Closed forms, because a mass-properties routine that is only checked against
//! another integrator is checked against a shared misunderstanding.
//!
//! Every fixture here has an exact answer derived from
//! `∫_simplex x^a y^b z^c = a!b!c!/(a+b+c+3)!` or from elementary geometry, so
//! the four deferred leading factors — `1/18`, `1/96`, `1/400`, `1/240`, none of
//! which the paper writes out in one place — have nowhere to hide.

use super::{MassProperties, mass_properties};
use crate::Error;

/// `[0,1]³` as twelve outward-wound triangles.
const CUBE_POSITIONS: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];

const CUBE_TRIANGLES: [[u32; 3]; 12] = [
    [0, 2, 1],
    [0, 3, 2],
    [4, 5, 6],
    [4, 6, 7],
    [0, 1, 5],
    [0, 5, 4],
    [1, 2, 6],
    [1, 6, 5],
    [2, 3, 7],
    [2, 7, 6],
    [3, 0, 4],
    [3, 4, 7],
];

/// The corner simplex: origin plus the three unit axis points.
const TET_POSITIONS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

/// Outward from the solid: the three faces on the coordinate planes wind so
/// their normals point along the negative axes, the slanted face along `+1`.
const TET_TRIANGLES: [[u32; 3]; 4] = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * b.abs().max(1.0)
}

#[test]
fn the_unit_cube_matches_its_closed_form() {
    let p = mass_properties(&CUBE_POSITIONS, &CUBE_TRIANGLES).expect("cube is a solid");

    assert!(close(p.volume, 1.0, 1e-15), "volume {}", p.volume);
    for axis in 0..3 {
        assert!(
            close(p.center_of_mass[axis], 0.5, 1e-15),
            "com {:?}",
            p.center_of_mass
        );
        // About the centre: I = m(a² + b²)/12 = 1/6 for a unit cube.
        assert!(
            close(p.inertia[axis][axis], 1.0 / 6.0, 1e-14),
            "inertia {:?}",
            p.inertia
        );
        // About the origin corner: ∫(y² + z²) dV = 2/3.
        assert!(
            close(p.inertia_about_origin[axis][axis], 2.0 / 3.0, 1e-14),
            "inertia about origin {:?}",
            p.inertia_about_origin
        );
    }
    // Products about the corner: −∫xy dV = −1/4. Centred, they vanish.
    assert!(close(p.inertia_about_origin[0][1], -0.25, 1e-14));
    assert!(p.inertia[0][1].abs() < 1e-15, "{:?}", p.inertia);
}

/// The tetrahedron is the fixture that separates a wrong *cubic* factor from a
/// wrong *quadratic* one: the cube's symmetry makes several sign errors
/// invisible and the simplex's does not.
#[test]
fn the_corner_simplex_matches_its_closed_form() {
    let p = mass_properties(&TET_POSITIONS, &TET_TRIANGLES).expect("tet is a solid");

    assert!(close(p.volume, 1.0 / 6.0, 1e-15), "volume {}", p.volume);
    for axis in 0..3 {
        assert!(close(p.center_of_mass[axis], 0.25, 1e-14));
        // ∫x² dV = 2!/5! = 1/60, so Iₓₓ about the origin is 2/60 = 1/30.
        assert!(
            close(p.inertia_about_origin[axis][axis], 1.0 / 30.0, 1e-14),
            "{:?}",
            p.inertia_about_origin
        );
    }
    // ∫xy dV = 1/5! = 1/120, negated by the physics convention.
    assert!(close(p.inertia_about_origin[0][1], -1.0 / 120.0, 1e-13));
    // Parallel axis: Iₓᵧ,com = −1/120 + V·cₓcᵧ = −1/120 + (1/6)(1/16) = +1/480,
    // and it changes SIGN, which is what makes it the assertion worth writing —
    // a shift applied the wrong way round lands on −1/80 rather than near zero.
    assert!(close(p.inertia[0][1], 1.0 / 480.0, 1e-13), "{:?}", p.inertia);
    // Iₓₓ,com = 1/30 − V·(cᵧ² + c𝓏²) = 1/30 − (1/6)(2/16) = 1/80.
    assert!(close(p.inertia[0][0], 1.0 / 80.0, 1e-13), "{:?}", p.inertia);
}

/// Translating the mesh must move the centroid by exactly that much and leave
/// the centred tensor alone. This is the property that catches a parallel-axis
/// shift applied in the wrong direction, which the origin-centred fixtures
/// above cannot see.
#[test]
fn the_centred_tensor_is_translation_invariant() {
    let shift = [3.5, -2.25, 11.0];
    let moved: [[f64; 3]; 4] = core::array::from_fn(|v| {
        core::array::from_fn(|axis| TET_POSITIONS[v][axis] + shift[axis])
    });

    let here = mass_properties(&TET_POSITIONS, &TET_TRIANGLES).expect("tet is a solid");
    let there = mass_properties(&moved, &TET_TRIANGLES).expect("shifted tet is a solid");

    assert!(close(there.volume, here.volume, 1e-13));
    for (axis, offset) in shift.iter().enumerate() {
        assert!(close(
            there.center_of_mass[axis],
            here.center_of_mass[axis] + offset,
            1e-13
        ));
        for other in 0..3 {
            assert!(
                close(there.inertia[axis][other], here.inertia[axis][other], 1e-11),
                "{:?} vs {:?}",
                there.inertia,
                here.inertia
            );
        }
    }
}

/// Both tensors are symmetric *exactly*, not to a tolerance — the mirrored
/// entry is the same expression, stored twice.
#[test]
fn both_tensors_are_bit_exactly_symmetric() {
    let p = mass_properties(&TET_POSITIONS, &TET_TRIANGLES).expect("tet is a solid");
    for i in 0..3 {
        for j in 0..3 {
            assert_eq!(p.inertia[i][j].to_bits(), p.inertia[j][i].to_bits());
            assert_eq!(
                p.inertia_about_origin[i][j].to_bits(),
                p.inertia_about_origin[j][i].to_bits()
            );
        }
    }
}

/// The leak detector, both ways round.
///
/// On a closed surface the two divergence forms of `∫xᵢxⱼ dV` differ only by
/// summation order, so the residual is round-off. Remove one triangle and the
/// theorem no longer applies: the residual jumps to the scale of the hole. A
/// harness that only ever saw the closed case could not tell the difference
/// between a working detector and a constant zero.
#[test]
fn asymmetry_separates_a_closed_mesh_from_a_leaking_one() {
    let closed = mass_properties(&CUBE_POSITIONS, &CUBE_TRIANGLES).expect("cube is a solid");
    // Shifted off the origin, where the two forms have something to disagree
    // about: at the origin every cubic moment of a symmetric fixture is tiny.
    let shifted: [[f64; 3]; 8] = core::array::from_fn(|v| {
        core::array::from_fn(|axis| CUBE_POSITIONS[v][axis] + [7.0, 5.0, 3.0][axis])
    });
    let far = mass_properties(&shifted, &CUBE_TRIANGLES).expect("cube is a solid");

    // Round-off floor, relative to the tensor it sits beside.
    assert!(
        closed.asymmetry <= 1e-14 * closed.inertia_about_origin[0][0].abs().max(1.0),
        "closed cube asymmetry {}",
        closed.asymmetry
    );
    assert!(
        far.asymmetry <= 1e-13 * far.inertia_about_origin[0][0].abs(),
        "shifted cube asymmetry {} against {}",
        far.asymmetry,
        far.inertia_about_origin[0][0]
    );

    let mut leaking = CUBE_TRIANGLES.to_vec();
    leaking.pop();
    let leaky = mass_properties(&shifted, &leaking).expect("still encloses a positive volume");
    assert!(
        leaky.asymmetry > 1e-3 * far.inertia_about_origin[0][0].abs(),
        "a missing triangle left asymmetry at {}",
        leaky.asymmetry
    );
}

#[test]
fn an_inward_wound_mesh_is_rejected_rather_than_flipped() {
    let inward: [[u32; 3]; 12] =
        core::array::from_fn(|t| [CUBE_TRIANGLES[t][0], CUBE_TRIANGLES[t][2], CUBE_TRIANGLES[t][1]]);
    let err = mass_properties(&CUBE_POSITIONS, &inward).expect_err("negative volume");
    assert!(matches!(err, Error::MassPropertiesUndefined { .. }), "{err:?}");
}

#[test]
fn an_empty_surface_encloses_nothing_and_says_so() {
    let err = mass_properties::<f64>(&[], &[]).expect_err("no solid");
    assert!(matches!(
        err,
        Error::MassPropertiesUndefined {
            volume: 0.0,
            largest_moment: 0.0
        }
    ));
}

#[test]
fn an_index_past_the_end_names_itself() {
    let err = mass_properties(&TET_POSITIONS, &[[0, 1, 9]]).expect_err("bad index");
    assert_eq!(
        err,
        Error::IndexOutOfRange {
            at: 2,
            index: 9,
            vertices: 4
        }
    );
}

/// `f32` is the default scalar everywhere else in the crate, so the module has
/// to work there too — at `f32`'s own accuracy, not at `f64`'s.
#[test]
fn f32_gets_the_same_answer_at_its_own_precision() {
    let positions: [[f32; 3]; 8] =
        core::array::from_fn(|v| core::array::from_fn(|axis| CUBE_POSITIONS[v][axis] as f32));
    let p: MassProperties<f32> =
        mass_properties(&positions, &CUBE_TRIANGLES).expect("cube is a solid");
    assert!((p.volume - 1.0).abs() <= 1e-6, "{}", p.volume);
    assert!((p.inertia[0][0] - 1.0 / 6.0).abs() <= 1e-6, "{:?}", p.inertia);
}
