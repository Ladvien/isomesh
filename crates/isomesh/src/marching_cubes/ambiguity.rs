//! MC33's face rule: Nielson & Hamann's asymptotic decider.
//!
//! Nielson, G. M. and Hamann, B., *The Asymptotic Decider: Resolving the
//! Ambiguity in Marching Cubes*, Proceedings of Visualization '91, pp. 83–91
//! (`10.1109/visual.1991.175782`).
//!
//! # What it decides
//!
//! An **ambiguous face** has alternating corner signs, so all four of its edges
//! are cut and the two inside corners can be either joined across the face or
//! separated. Marching Cubes proper always separates them (see
//! [`super::table`]). That choice is crack-free, because it is a function of the
//! face's own four corner signs, but it is arbitrary: it agrees with the field's
//! own bilinear interpolant on the face only by luck.
//!
//! The decider replaces it with the interpolant's own answer. On the unit square
//! the bilinear function has a single saddle point, and the sign of the function
//! *at that saddle* says which pair of diagonally opposite corners the level set
//! leaves connected. For corner values `v0..v3` counter-clockwise, with the
//! diagonals `{v0, v2}` and `{v1, v3}`:
//!
//! ```text
//! S = (v0·v2 − v1·v3) / (v0 + v2 − v1 − v3)
//! ```
//!
//! # Why there is no division and no epsilon
//!
//! On an ambiguous face one diagonal is strictly negative — this crate's inside
//! — and the other is non-negative, because a sample of exactly zero counts as
//! outside (`crate::cube::is_inside`). So the denominator is a strictly negative
//! sum minus a non-negative one, or the reverse: **it can never be zero**, and
//! its sign is known from the signs alone. The implementation brief's "guard the
//! denominator near zero" is unnecessary here, for the same reason
//! `crate::cube::edge_crossing` needs no epsilon.
//!
//! Only `sign(S)` is wanted, and with the denominator's sign already known that
//! reduces to comparing the two diagonal products. Writing `d_in` for the
//! product of the inside diagonal and `d_out` for the outside one, both branches
//! of the derivation collapse to the same test:
//!
//! ```text
//! joined  ⟺  d_in > d_out
//! ```
//!
//! Ties resolve to *separated*, matching Marching Cubes. `d_in == d_out` is
//! `S == 0` exactly: the two branches of the hyperbola degenerate into crossing
//! straight lines and the regions meet at a single point, so either answer is
//! defensible and the one that changes nothing is preferred.
//!
//! # Why this is still crack-free
//!
//! Two cells meeting on a face read the *same four sample values*, in orders
//! that differ by a rotation and a reflection — the shared face is
//! counter-clockwise-from-outside for each of them, and their outsides are
//! opposite. Neither transformation changes which diagonal is which, and IEEE
//! multiplication is commutative and correctly rounded, so `v0·v2` and `v2·v0`
//! are the same bits. The two cells therefore cannot disagree.
//! `the_decider_is_invariant_under_rotation_and_reflection` is that argument as
//! a test.

#[cfg(test)]
mod tests;

use crate::Real;
use crate::cube::{face_corners, is_inside};

use super::table::face_bit;

/// How a cell resolves an ambiguous face.
///
/// Both settings are crack-free. They differ in which surface they produce
/// there, and therefore in the topology of the result — see A-002's archive
/// entry for the measured difference.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum FaceAmbiguity {
    /// Separate the two inside corners. Marching Cubes proper, and the default.
    ///
    /// The default because it is what A-001 shipped and what the committed
    /// golden hashes pin; switching the default would silently restate every
    /// measurement taken against plain Marching Cubes.
    #[default]
    Separate,
    /// Ask the bilinear interpolant, via the asymptotic decider. MC33's face
    /// half.
    AsymptoticDecider,
}

/// Does the surface leave the face's two inside corners joined?
///
/// `v` is the face's four corner values, counter-clockwise as seen from outside
/// the cube — the order [`crate::cube::face_corners`] produces. The result is
/// invariant under rotating or reflecting that order, which is what lets two
/// adjacent cells agree.
///
/// # Panics
///
/// In debug builds, if the face is not ambiguous. There is nothing to decide
/// then, and asking would mean the caller has confused a face's identity.
#[must_use]
pub fn face_is_joined<R: Real>(v: [R; 4]) -> bool {
    debug_assert!(
        is_inside(v[0]) == is_inside(v[2])
            && is_inside(v[1]) == is_inside(v[3])
            && is_inside(v[0]) != is_inside(v[1]),
        "the asymptotic decider only applies to an ambiguous face"
    );

    let d02 = v[0] * v[2];
    let d13 = v[1] * v[3];
    // Whichever diagonal is inside supplies `d_in`. Both arms are the same
    // comparison; they differ only in which product plays which role.
    if is_inside(v[0]) {
        d02 > d13
    } else {
        d13 > d02
    }
}

/// The per-face resolution mask for one cell, over the faces `ambiguous` marks.
///
/// Bit layout is [`face_bit`]'s. Faces outside `ambiguous` are left clear, since
/// [`super::table::segment_links`] ignores them anyway and leaving them clear
/// keeps the mask canonical for the validator.
pub(crate) fn joined_mask<R: Real>(corner_value: &[R; 8], ambiguous: u8) -> u8 {
    let mut mask = 0u8;
    for axis in 0..3usize {
        for side in 0..2u8 {
            let bit = face_bit(axis, side);
            if ambiguous & bit == 0 {
                continue;
            }
            let c = face_corners(axis, side);
            let v = [
                corner_value[c[0] as usize],
                corner_value[c[1] as usize],
                corner_value[c[2] as usize],
                corner_value[c[3] as usize],
            ];
            if face_is_joined(v) {
                mask |= bit;
            }
        }
    }
    mask
}
