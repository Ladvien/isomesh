//! Tests for the asymptotic decider.
//!
//! Three lines of defence, in increasing order of how much they would have cost
//! to get wrong:
//!
//! 1. The two anchors — a strongly negative diagonal joins, a strongly positive
//!    one separates. These would catch an inverted comparison.
//! 2. Invariance under rotating and reflecting the corner order. This is the
//!    crack-free argument, and nothing else in the suite sees a violation of it
//!    directly: a decider that disagreed between two cells would show up only as
//!    boundary edges on some field, eventually, on some resolution.
//! 3. Agreement with the bilinear interpolant itself, by dense sampling and
//!    flood fill. This is the one that says the answer is *right* rather than
//!    merely consistent, and it is what A-002's acceptance criterion means by
//!    "matching the trilinear connectivity" on a face.

use alloc::vec;
use alloc::vec::Vec;

use super::{FaceAmbiguity, face_is_joined, joined_mask};
use crate::cube::{face_corners, is_inside};
use crate::mc::table::{AMBIGUOUS_FACES, face_bit};

// ─── anchors ────────────────────────────────────────────────────────────────

/// A deep inside diagonal against a shallow outside one: the inside region owns
/// the middle of the face, so the two inside corners are joined.
#[test]
fn a_dominant_inside_diagonal_joins() {
    assert!(face_is_joined([-2.0f64, 1.0, -2.0, 1.0]));
    assert!(face_is_joined([1.0f64, -2.0, 1.0, -2.0]));
}

/// And the reverse: a shallow inside diagonal loses the middle to the outside,
/// so the inside corners are cut off separately.
#[test]
fn a_dominant_outside_diagonal_separates() {
    assert!(!face_is_joined([-1.0f64, 2.0, -1.0, 2.0]));
    assert!(!face_is_joined([2.0f64, -1.0, 2.0, -1.0]));
}

/// `d_in == d_out` is `S == 0` exactly — the hyperbola degenerates to two
/// crossing lines and the regions touch at a point. Separated is the documented
/// tie-break, and it is the one that agrees with plain Marching Cubes.
#[test]
fn a_degenerate_saddle_separates() {
    assert!(!face_is_joined([-1.0f64, 1.0, -1.0, 1.0]));
    assert!(!face_is_joined([-2.0f64, 1.0, -0.5, 1.0]));
}

/// A sample of exactly zero is outside (`crate::cube::is_inside`), so an outside
/// diagonal can have a zero product while the inside one cannot. Joined, and no
/// division happens that could turn that into a NaN.
#[test]
fn a_zero_on_the_outside_diagonal_joins() {
    assert!(face_is_joined([-1.0f64, 0.0, -1.0, 3.0]));
    assert!(face_is_joined([-1.0f64, 0.0, -1.0, 0.0]));
}

// ─── the crack-free property ────────────────────────────────────────────────

/// Two cells meeting on a face list its corners in orders that differ by a
/// rotation and a reflection. If the answer moved under either, the two cells
/// would build different surfaces on a shared face — a crack.
///
/// Bit-identical, not merely equal: the products are the same two IEEE
/// multiplications in either order, because multiplication is commutative and
/// correctly rounded.
#[test]
fn the_decider_is_invariant_under_rotation_and_reflection() {
    for v in ambiguous_quadruples() {
        let expected = face_is_joined(v);
        for r in 0..4usize {
            let rotated = [v[r], v[(r + 1) % 4], v[(r + 2) % 4], v[(r + 3) % 4]];
            assert_eq!(face_is_joined(rotated), expected, "rotation {r} of {v:?}");

            let reflected = [rotated[0], rotated[3], rotated[2], rotated[1]];
            assert_eq!(
                face_is_joined(reflected),
                expected,
                "reflected rotation {r} of {v:?}"
            );
        }
    }
}

// ─── against the interpolant itself ─────────────────────────────────────────

/// The decider must agree with the bilinear interpolant it claims to read.
///
/// Ground truth is computed the slow, obvious way: sample `B(s, t)` on a dense
/// grid over the face, flood-fill the negative region with 4-connectivity, and
/// ask whether the two inside corners land in one component. No formula is
/// shared with the implementation, so an algebra error in the derivation cannot
/// hide in both.
///
/// Near-degenerate saddles are skipped rather than asserted, because the neck
/// between the two hyperbola branches is about `2·sqrt(|S/a|)` wide and a finite
/// grid cannot see it below its own spacing. The skip count is recorded so the
/// test cannot quietly become vacuous.
#[test]
fn the_decider_agrees_with_the_bilinear_interpolant() {
    const N: usize = 256;
    // Comfortably above the (2/N)² the grid can resolve.
    const TOL: f64 = 2e-3;

    let mut checked = 0usize;
    let mut skipped = 0usize;
    let mut joined_seen = 0usize;

    for v in ambiguous_quadruples() {
        // B(s,t) = a·s·t + b·s + c·t + d over the unit square, with the corners
        // in the order face_corners produces: (0,0), (1,0), (1,1), (0,1).
        let a = v[0] - v[1] + v[2] - v[3];
        let saddle = (v[0] * v[2] - v[1] * v[3]) / (v[0] + v[2] - v[1] - v[3]);
        if saddle.abs() <= TOL * a.abs() {
            skipped += 1;
            continue;
        }

        let truth = negative_corners_are_connected::<N>(v);
        let decided = face_is_joined(v);
        assert_eq!(
            decided, truth,
            "values {v:?}: decider says joined={decided}, dense sampling says {truth} \
             (saddle {saddle})"
        );
        checked += 1;
        joined_seen += usize::from(truth);
    }

    // Both answers have to occur, or the agreement is trivial.
    assert!(checked > 100, "only {checked} quadruples resolved");
    assert!(
        joined_seen > 0 && joined_seen < checked,
        "joined={joined_seen} of {checked}"
    );
    assert!(
        skipped * 4 < checked,
        "{skipped} skipped against {checked} checked — the tolerance is eating the test"
    );
}

// ─── the per-cell mask ──────────────────────────────────────────────────────

/// `joined_mask` must consult the faces it is handed and no others, and it must
/// place each answer in that face's own bit.
#[test]
fn joined_mask_reads_only_the_faces_it_is_told_to() {
    // Corner values chosen so that every face the case marks ambiguous is
    // decided "joined": the inside corners are deep, the outside ones shallow.
    let corner_value = [-4.0f64, 0.5, 0.5, -4.0, 0.5, -4.0, -4.0, 0.5];
    let case = case_of(&corner_value);
    let ambiguous = AMBIGUOUS_FACES[case as usize];
    assert_ne!(ambiguous, 0, "fixture has no ambiguous face");

    assert_eq!(joined_mask(&corner_value, ambiguous), ambiguous);
    assert_eq!(joined_mask(&corner_value, 0), 0);

    // One face at a time: the answer for a face never lands in another's bit.
    for axis in 0..3usize {
        for side in 0..2u8 {
            let bit = face_bit(axis, side);
            if ambiguous & bit == 0 {
                continue;
            }
            assert_eq!(joined_mask(&corner_value, bit), bit);
        }
    }
}

/// The default is Marching Cubes proper. Changing it would silently restate
/// every measurement taken against A-001.
#[test]
fn the_default_is_separate() {
    assert_eq!(FaceAmbiguity::default(), FaceAmbiguity::Separate);
}

// ─── helpers ────────────────────────────────────────────────────────────────

/// Every ambiguous four-value face over a fixed ladder of magnitudes.
///
/// Systematic, not chosen: picking a fixture by eye has landed in the
/// degenerate region twice in this repo (M-32, M-38). The ladder spans two
/// orders of magnitude either side of zero and includes exactly zero, which is
/// the boundary case the sign rule turns on.
fn ambiguous_quadruples() -> Vec<[f64; 4]> {
    const LADDER: [f64; 7] = [-4.0, -1.0, -0.25, 0.0, 0.25, 1.0, 4.0];
    let mut out = Vec::new();
    for &v0 in &LADDER {
        for &v1 in &LADDER {
            for &v2 in &LADDER {
                for &v3 in &LADDER {
                    let v = [v0, v1, v2, v3];
                    if is_inside(v0) == is_inside(v2)
                        && is_inside(v1) == is_inside(v3)
                        && is_inside(v0) != is_inside(v1)
                    {
                        out.push(v);
                    }
                }
            }
        }
    }
    out
}

/// Are the face's two inside corners in one connected component of `{B < 0}`?
///
/// `N + 1` samples per axis so the four corners are hit exactly. 4-connectivity:
/// away from the degenerate band the neck at the saddle is a region of positive
/// width, so diagonal-only links do not arise.
fn negative_corners_are_connected<const N: usize>(v: [f64; 4]) -> bool {
    let stride = N + 1;
    let at = |i: usize, j: usize| -> f64 {
        let s = i as f64 / N as f64;
        let t = j as f64 / N as f64;
        v[0] * (1.0 - s) * (1.0 - t) + v[1] * s * (1.0 - t) + v[2] * s * t + v[3] * (1.0 - s) * t
    };
    let index = |i: usize, j: usize| j * stride + i;

    // The two inside corners, in (i, j) grid coordinates.
    let (from, to) = if is_inside(v[0]) {
        ((0, 0), (N, N)) // corners c0 and c2
    } else {
        ((N, 0), (0, N)) // corners c1 and c3
    };

    let mut seen = vec![false; stride * stride];
    let mut stack = vec![from];
    seen[index(from.0, from.1)] = true;
    while let Some((i, j)) = stack.pop() {
        if (i, j) == to {
            return true;
        }
        let push = |i: usize, j: usize, stack: &mut Vec<(usize, usize)>, seen: &mut [bool]| {
            let k = index(i, j);
            if !seen[k] && is_inside(at(i, j)) {
                seen[k] = true;
                stack.push((i, j));
            }
        };
        if i > 0 {
            push(i - 1, j, &mut stack, &mut seen);
        }
        if i < N {
            push(i + 1, j, &mut stack, &mut seen);
        }
        if j > 0 {
            push(i, j - 1, &mut stack, &mut seen);
        }
        if j < N {
            push(i, j + 1, &mut stack, &mut seen);
        }
    }
    false
}

/// The 8-bit case index for a set of corner values.
fn case_of(corner_value: &[f64; 8]) -> u8 {
    let mut case = 0u8;
    for (c, &v) in corner_value.iter().enumerate() {
        if is_inside(v) {
            case |= 1 << c;
        }
    }
    case
}

/// Unused-import guard: `face_corners` is what defines the corner order the
/// decider is documented against, and this pins that the order is a 4-cycle of
/// cube-adjacent corners so the diagonals really are `{0, 2}` and `{1, 3}`.
#[test]
fn the_diagonals_are_corners_zero_two_and_one_three() {
    for axis in 0..3usize {
        for side in 0..2u8 {
            let c = face_corners(axis, side);
            // Adjacent corners differ in one bit; diagonal ones differ in two.
            assert_eq!((c[0] ^ c[2]).count_ones(), 2, "axis {axis} side {side}");
            assert_eq!((c[1] ^ c[3]).count_ones(), 2, "axis {axis} side {side}");
        }
    }
}
