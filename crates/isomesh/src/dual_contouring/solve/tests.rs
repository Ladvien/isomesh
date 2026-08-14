//! The load-bearing tests here are the corner (which has a closed form), the
//! lattice-rotation equivariance, and the negative control that shows the
//! magnitude sort is what buys it.

use alloc::vec;
use alloc::vec::Vec;

use super::{LAMBDA, dot_equivariant, solve};
use crate::hermite::{HermiteCell, HermiteCrossing};

fn crossing(position: [f64; 3], normal: [f64; 3]) -> HermiteCrossing<f64> {
    let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    HermiteCrossing {
        position,
        normal: [normal[0] / len, normal[1] / len, normal[2] / len],
    }
}

fn cell(pairs: &[(u8, HermiteCrossing<f64>)]) -> HermiteCell<f64> {
    HermiteCell::from_crossings(pairs)
}

fn close(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
    (0..3).all(|i| (a[i] - b[i]).abs() <= tol)
}

/// Three orthogonal planes through a common point: the case the whole rule
/// exists for, and the one with an exact answer.
///
/// With unit axis normals `M = I`, so `A = (1+λ)I`, `adj(A) = (1+λ)²I` and
/// `det(A) = (1+λ)³`. Every plane passing through `q` makes `dᵢ = nᵢ·(q − c)`,
/// hence `g = M(q − c) = q − c`, and the whole solve collapses to
///
/// ```text
/// x = c + (q − c)/(1 + λ)
/// ```
///
/// So the vertex lands on the segment from centroid to corner, short of the
/// corner by exactly `λ/(1+λ)` of the way — **0.990099…** of it at `λ = 0.01`.
/// That is the regularizer's cost in sharpness, in closed form, and it is why
/// `λ` is small.
#[test]
fn a_perfect_corner_lands_at_the_closed_form() {
    let q = [1.0, 1.0, 1.0];
    let c = cell(&[
        (0, crossing([1.0, 0.4, 0.7], [1.0, 0.0, 0.0])),
        (1, crossing([0.3, 1.0, 0.9], [0.0, 1.0, 0.0])),
        (2, crossing([0.6, 0.2, 1.0], [0.0, 0.0, 1.0])),
    ]);

    let centroid = c.centroid().expect("three crossings");
    let x = solve(&c).expect("solvable");

    let expected = [
        centroid[0] + (q[0] - centroid[0]) / (1.0 + LAMBDA),
        centroid[1] + (q[1] - centroid[1]) / (1.0 + LAMBDA),
        centroid[2] + (q[2] - centroid[2]) / (1.0 + LAMBDA),
    ];
    assert!(close(x, expected, 1e-12), "got {x:?}, want {expected:?}");

    // And it is genuinely near the corner rather than near the centroid, which
    // is the difference between dual contouring and surface nets.
    let to_corner = (0..3).map(|i| (x[i] - q[i]).powi(2)).sum::<f64>().sqrt();
    let centroid_to_corner = (0..3)
        .map(|i| (centroid[i] - q[i]).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(
        to_corner < 0.02 * centroid_to_corner,
        "vertex should sit at the corner: {to_corner} vs {centroid_to_corner}"
    );
}

/// A flat region gives a rank-1 `M`, which is exactly the case an unregularized
/// solve cannot handle — it is one equation in three unknowns.
///
/// The vertex must stay put laterally and only move onto the plane. With all
/// normals equal, `A n = (k+λ)n`, so the offset is `D/(k+λ)` along `n` and
/// nothing across it — bounded, not a blow-up.
#[test]
fn a_flat_region_does_not_fly_off() {
    let n = [0.0, 1.0, 0.0];
    let c = cell(&[
        (0, crossing([0.1, 0.5, 0.2], n)),
        (1, crossing([0.9, 0.5, 0.3], n)),
        (2, crossing([0.4, 0.5, 0.8], n)),
        (3, crossing([0.7, 0.5, 0.6], n)),
    ]);

    let centroid = c.centroid().expect("four crossings");
    let x = solve(&c).expect("solvable");

    // Lateral position unchanged: the solve has no information across the plane
    // and must not invent any.
    assert!((x[0] - centroid[0]).abs() < 1e-12, "moved in x: {x:?}");
    assert!((x[2] - centroid[2]).abs() < 1e-12, "moved in z: {x:?}");
    // And it stays in the neighbourhood of the plane rather than diverging.
    assert!((x[1] - 0.5).abs() < 0.5, "left the plane: {x:?}");
}

/// A single crossing is the most under-determined case there is. It must still
/// produce a finite answer near the crossing.
#[test]
fn one_crossing_is_finite_and_local() {
    let c = cell(&[(5, crossing([0.5, 0.5, 0.5], [1.0, 1.0, 0.0]))]);
    let x = solve(&c).expect("solvable");
    assert!(x.iter().all(|v| v.is_finite()), "{x:?}");
    assert!(close(x, [0.5, 0.5, 0.5], 0.5), "{x:?}");
}

/// Two parallel, opposed planes — a slab thinner than the cell. Degenerate for a
/// three-plane rule, ordinary for this one.
#[test]
fn opposed_planes_are_handled_by_the_same_arithmetic() {
    let c = cell(&[
        (0, crossing([0.4, 0.2, 0.1], [0.0, 1.0, 0.0])),
        (1, crossing([0.6, 0.6, 0.9], [0.0, -1.0, 0.0])),
    ]);
    let x = solve(&c).expect("solvable");
    assert!(x.iter().all(|v| v.is_finite()), "{x:?}");
}

#[test]
fn an_empty_cell_has_no_vertex() {
    assert!(solve(&cell(&[])).is_none());
}

// ─── equivariance ───────────────────────────────────────────────────────────

/// The 24 rotations of the cube, as signed axis permutations.
///
/// Generated rather than tabulated: every permutation of the three axes crossed
/// with every sign pattern, keeping the 24 with determinant `+1`. A transcribed
/// table would be one more thing to get wrong.
fn octahedral_rotations() -> Vec<[[f64; 3]; 3]> {
    let mut out = Vec::new();
    for perm in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        for signs in 0..8u8 {
            let s = [
                if signs & 1 == 0 { 1.0 } else { -1.0 },
                if signs & 2 == 0 { 1.0 } else { -1.0 },
                if signs & 4 == 0 { 1.0 } else { -1.0 },
            ];
            let mut m = [[0.0f64; 3]; 3];
            for row in 0..3 {
                m[row][perm[row]] = s[row];
            }
            // det = +1 keeps the proper rotations and drops the reflections.
            let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
            if det > 0.0 {
                out.push(m);
            }
        }
    }
    out
}

fn apply(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// A representative spread of cells, deliberately including asymmetric ones —
/// a symmetric fixture would be equivariant under any rule at all.
fn fixtures() -> Vec<Vec<(u8, HermiteCrossing<f64>)>> {
    vec![
        vec![
            (0, crossing([0.83, 0.21, 0.44], [1.0, 0.0, 0.0])),
            (1, crossing([0.17, 0.92, 0.36], [0.0, 1.0, 0.0])),
            (2, crossing([0.55, 0.63, 0.71], [0.0, 0.0, 1.0])),
        ],
        vec![
            (3, crossing([0.31, 0.77, 0.12], [0.7, 0.3, -0.2])),
            (4, crossing([0.68, 0.24, 0.95], [-0.1, 0.9, 0.4])),
            (5, crossing([0.49, 0.58, 0.37], [0.2, -0.4, 0.8])),
            (6, crossing([0.13, 0.86, 0.61], [0.5, 0.5, 0.5])),
        ],
        vec![
            (7, crossing([0.29, 0.41, 0.88], [0.33, -0.81, 0.19])),
            (8, crossing([0.74, 0.15, 0.52], [-0.62, 0.27, 0.74])),
        ],
    ]
}

/// **The headline property.** Rotating a cell by any of the cube's 24 rotations
/// must rotate its vertex by the same rotation — *exactly*, to the bit.
///
/// This is what the vertex rule was chosen for. It is why a brush dragged across
/// a lattice does not pop: the answer depends on the geometry, not on which axis
/// happens to be called `x`.
#[test]
fn the_vertex_is_bit_exactly_equivariant_under_lattice_rotations() {
    let rotations = octahedral_rotations();
    assert_eq!(rotations.len(), 24, "the octahedral group has 24 rotations");

    let mut checked = 0u32;
    for fixture in fixtures() {
        let base = solve(&cell(&fixture)).expect("solvable");
        for m in &rotations {
            let rotated: Vec<(u8, HermiteCrossing<f64>)> = fixture
                .iter()
                .map(|&(e, c)| {
                    (
                        e,
                        HermiteCrossing {
                            position: apply(*m, c.position),
                            normal: apply(*m, c.normal),
                        },
                    )
                })
                .collect();
            let got = solve(&cell(&rotated)).expect("solvable");
            let want = apply(*m, base);
            for axis in 0..3 {
                assert_eq!(
                    got[axis].to_bits(),
                    want[axis].to_bits(),
                    "rotation {m:?} axis {axis}: got {got:?}, want {want:?}"
                );
            }
            checked += 1;
        }
    }
    assert_eq!(checked, 24 * 3);
}

/// **A-016's property, and the half the test above structurally cannot see.**
///
/// The vertex must be a function of the *set* of crossings, never of the edge
/// labels they arrived under. A lattice rotation permutes those labels — and
/// `the_vertex_is_bit_exactly_equivariant_under_lattice_rotations` carries `e`
/// through unchanged at its own line 216, so both sides of its comparison visit
/// the crossings in the same order and any accumulation-order defect cancels out
/// of the difference it measures.
///
/// Relabelling is the isolated form of the same permutation: nothing geometric
/// changes, because `solve_with` reads only `position` and `normal`. Anything
/// but bit-identity means the vertex is a function of the labelling.
#[test]
fn the_vertex_does_not_depend_on_which_edges_the_crossings_arrived_on() {
    let mut checked = 0u32;
    for fixture in fixtures() {
        let base = solve(&cell(&fixture)).expect("solvable");

        // Every relabelling keeps the crossings and moves only their names,
        // including ones that reverse the visit order outright.
        let labellings: [&[u8]; 4] = [
            &[11, 10, 9, 8],
            &[0, 5, 7, 11],
            &[9, 2, 6, 1],
            &[4, 3, 1, 0],
        ];
        for labels in labellings {
            let relabelled: Vec<(u8, HermiteCrossing<f64>)> = fixture
                .iter()
                .zip(labels)
                .map(|(&(_, c), &e)| (e, c))
                .collect();
            assert_eq!(relabelled.len(), fixture.len());

            let got = solve(&cell(&relabelled)).expect("solvable");
            for axis in 0..3 {
                assert_eq!(
                    got[axis].to_bits(),
                    base[axis].to_bits(),
                    "labels {labels:?} axis {axis}: got {got:?}, want {base:?}"
                );
            }
            checked += 1;
        }
    }
    assert_eq!(checked, 4 * 3);
}

/// The negative control for the test above, and the reason it is not decoration.
///
/// The accumulation this replaced summed in visit order. Reproduced here on one
/// fixture: the same crossings under two labellings give two different centroids
/// once the sum is taken in label order, and the QEF is solved *relative to that
/// centroid*, so the vertex moves with it.
///
/// If this ever stops finding a disagreement, the fixture has gone too tame to
/// detect the defect and the test above is proving nothing.
#[test]
fn an_in_visit_order_sum_really_does_depend_on_the_labelling() {
    let fixture = &fixtures()[1];
    let in_visit_order = |labels: &[u8]| {
        let mut pairs: Vec<(u8, [f64; 3])> = fixture
            .iter()
            .zip(labels)
            .map(|(&(_, c), &e)| (e, c.position))
            .collect();
        pairs.sort_by_key(|&(e, _)| e);
        let mut sum = [0.0f64; 3];
        for (_, p) in &pairs {
            for axis in 0..3 {
                sum[axis] += p[axis];
            }
        }
        sum
    };
    let a = in_visit_order(&[0, 1, 2, 3]);
    let b = in_visit_order(&[11, 10, 9, 8]);
    assert!(
        (0..3).any(|axis| a[axis].to_bits() != b[axis].to_bits()),
        "fixture cannot detect an accumulation-order defect: {a:?} vs {b:?}"
    );
}

/// The negative control for the rotation test.
///
/// A dot product summed in index order is not permutation invariant, so under a
/// rotation that relabels the axes the three products are added in a different
/// order and the last bits differ. The audit measured **4328 of 9600** lattice
/// trials failing that way, and **0 of 9600** with the magnitude sort.
///
/// If this test ever stops finding a disagreement, `dot_equivariant`'s sort has
/// stopped being load-bearing and the equivariance test above is proving
/// nothing.
#[test]
fn an_unsorted_dot_product_would_break_equivariance() {
    fn dot_naive(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    // Magnitudes spread far enough apart that the summation order decides the
    // last bits -- which is exactly the situation a real cell is in when one
    // component of a normal dominates.
    let a = [1.0, 1e-16, 1e-8];
    let b = [1.0, 1.0, 1.0];
    let swapped_a = [a[1], a[2], a[0]];
    let swapped_b = [b[1], b[2], b[0]];

    assert_ne!(
        dot_naive(a, b).to_bits(),
        dot_naive(swapped_a, swapped_b).to_bits(),
        "the naive dot must be order-dependent, or this control proves nothing"
    );
    assert_eq!(
        dot_equivariant(a, b).to_bits(),
        dot_equivariant(swapped_a, swapped_b).to_bits(),
        "the sorted dot must not be"
    );
}

/// Translation invariance: the rule is stated relative to the centroid, so
/// moving a cell must move its vertex by exactly the same amount.
#[test]
fn the_vertex_is_translation_equivariant() {
    for fixture in fixtures() {
        let base = solve(&cell(&fixture)).expect("solvable");
        for shift in [[3.0, 0.0, 0.0], [-7.0, 2.0, 11.0], [100.0, 100.0, 100.0]] {
            let moved: Vec<(u8, HermiteCrossing<f64>)> = fixture
                .iter()
                .map(|&(e, c)| {
                    (
                        e,
                        HermiteCrossing {
                            position: [
                                c.position[0] + shift[0],
                                c.position[1] + shift[1],
                                c.position[2] + shift[2],
                            ],
                            normal: c.normal,
                        },
                    )
                })
                .collect();
            let got = solve(&cell(&moved)).expect("solvable");
            let want = [base[0] + shift[0], base[1] + shift[1], base[2] + shift[2]];
            assert!(
                close(got, want, 1e-9),
                "shift {shift:?}: got {got:?}, want {want:?}"
            );
        }
    }
}

/// `λ` is added to a dimensionless matrix, so the rule must not care what unit
/// the cell is measured in. A cell scaled by `s` must give a vertex scaled by
/// `s`.
#[test]
fn the_rule_is_scale_invariant() {
    for fixture in fixtures() {
        let base = solve(&cell(&fixture)).expect("solvable");
        for scale in [1e-3, 0.5, 2.0, 1e3] {
            let scaled: Vec<(u8, HermiteCrossing<f64>)> = fixture
                .iter()
                .map(|&(e, c)| {
                    (
                        e,
                        HermiteCrossing {
                            position: [
                                c.position[0] * scale,
                                c.position[1] * scale,
                                c.position[2] * scale,
                            ],
                            // Normals are directions: unchanged by a scale.
                            normal: c.normal,
                        },
                    )
                })
                .collect();
            let got = solve(&cell(&scaled)).expect("solvable");
            let want = [base[0] * scale, base[1] * scale, base[2] * scale];
            let tol = 1e-9 * scale.max(1.0);
            assert!(
                close(got, want, tol),
                "scale {scale}: got {got:?}, want {want:?}"
            );
        }
    }
}

/// Determinism, in the sense T-004 means it: same input, bit-identical output.
#[test]
fn the_solve_is_deterministic() {
    for fixture in fixtures() {
        let a = solve(&cell(&fixture)).expect("solvable");
        let b = solve(&cell(&fixture)).expect("solvable");
        for axis in 0..3 {
            assert_eq!(a[axis].to_bits(), b[axis].to_bits());
        }
    }
}

/// The same fixture in `f32`, because the crate promises both widths and the
/// module docs make a specific claim about `f32` being the risky one.
#[test]
fn f32_solves_the_corner_too() {
    let c = HermiteCell::<f32>::from_crossings(&[
        (
            0,
            HermiteCrossing {
                position: [1.0, 0.4, 0.7],
                normal: [1.0, 0.0, 0.0],
            },
        ),
        (
            1,
            HermiteCrossing {
                position: [0.3, 1.0, 0.9],
                normal: [0.0, 1.0, 0.0],
            },
        ),
        (
            2,
            HermiteCrossing {
                position: [0.6, 0.2, 1.0],
                normal: [0.0, 0.0, 1.0],
            },
        ),
    ]);
    let x = solve(&c).expect("solvable");
    let centroid = c.centroid().expect("three crossings");
    let lambda = LAMBDA as f32;
    for axis in 0..3 {
        let want = centroid[axis] + (1.0 - centroid[axis]) / (1.0 + lambda);
        assert!(
            (x[axis] - want).abs() < 1e-6,
            "axis {axis}: got {x:?}, want {want}"
        );
    }
}
