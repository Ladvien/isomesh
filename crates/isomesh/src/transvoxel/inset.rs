//! Making room for the transition cells, which is Lengyel's Equation 4.2.
//!
//! A transition cell occupies a slab of thickness `w` on the **low-resolution**
//! side of a block boundary. Something has to vacate that slab, and it is the
//! coarse block's own boundary cells — §4.4:
//!
//! > If adjacent blocks are rendered at different levels of detail … then regular
//! > cells in the block for the lower level of detail must be **scaled in one or
//! > more directions to make space for transition cells** … we store two positions
//! > for each vertex belonging to regular cells on the boundary of a
//! > low-resolution block, a **primary** position used when transition cells are
//! > not rendered, and a **secondary** position used when transition cells are
//! > rendered.
//!
//! # The formula, and why it is simpler than it looks
//!
//! Equation 4.2 is written in level-0 cell units, which is what makes the `2^−k`
//! appear:
//!
//! ```text
//! Δx = (1 − 2^−k·x)·w(k)          if x < 2^k
//!      0                          if 2^k ≤ x ≤ 2^k(s − 1)
//!      (s − 1 − 2^−k·x)·w(k)      if x > 2^k(s − 1)
//! ```
//!
//! Substituting `x = c·2^k`, where `c` is the coordinate in **this block's own
//! cells**, the level index cancels completely:
//!
//! ```text
//! Δ = (1 − c)·w        if c < 1
//!     0                if 1 ≤ c ≤ s − 1
//!     (s − 1 − c)·w     if c > s − 1
//! ```
//!
//! So it is a linear taper across the first and last cell of the block, and
//! nothing in between moves. A vertex exactly on the boundary plane (`c = 0`)
//! moves by exactly `w`, which is precisely the displacement
//! [`TransitionCell`](super::cell::TransitionCell) gives its half-resolution
//! face — and that coincidence is what keeps the seam closed at a non-zero width.
//! `the_seam_stays_closed_at_a_real_width` is what checks it rather than trusting
//! the algebra.
//!
//! # This is a post-pass, deliberately
//!
//! Lengyel stores two positions per boundary vertex because his vertex program
//! picks between them per frame. That is a renderer's concern: which of a block's
//! neighbours are coarser is not known when the block is meshed, and can change
//! while it is resident. So the core crate meshes the block once, and
//! [`inset_boundary`] applies the taper to a copy when a transition is actually
//! being rendered on some face — the primary mesh is the un-inset one, and the
//! secondary is this.

use crate::{MeshBuffer, Real};

/// Which faces of a block have a transition rendered on them.
///
/// Bit `axis * 2 + side`, where side `0` is the low face and `1` the high one —
/// the same layout [`marching_cubes::table::face_bit`](crate::marching_cubes::table::face_bit)
/// uses, stated here because a transposition between the two would be invisible.
#[inline]
#[must_use]
pub const fn face_bit(axis: usize, side: u8) -> u8 {
    1 << (axis * 2 + side as usize)
}

/// Every face. The case where all six neighbours are coarser.
pub const ALL_FACES: u8 = 0b0011_1111;

/// Apply Equation 4.2's taper to a block's mesh, in place.
///
/// `origin` is the world position of the block's minimum corner, `cells` is `s` —
/// its size in **its own** cells — and `cell_size` is its own spacing. `width` is
/// the transition width, in world units. `faces` selects which faces to make room
/// on; a face with no coarser neighbour must not be tapered, or the block will
/// pull away from a same-resolution neighbour and open a seam where there was
/// none.
///
/// Only positions move. Normals are left alone: the taper is a shear of at most
/// `w` over one cell and re-deriving normals from it would report the *shear's*
/// geometry rather than the field's. Anything wanting normals from the displaced
/// geometry should say so explicitly through
/// [`normals::recompute`](crate::normals::recompute).
///
/// # Errors
///
/// [`Error::InvalidCellSize`](crate::Error::InvalidCellSize) if `cell_size` is
/// not finite and positive, or if `width` is negative or not finite. A width of
/// zero is allowed and is a no-op — that is the configuration A-011b shipped, and
/// M-74 records what it costs.
///
/// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if `cells` is less than
/// two. The taper needs a first cell and a last cell to be different cells; with
/// one, the two tapers would fight over the same vertices.
pub fn inset_boundary<R: Real>(
    mesh: &mut MeshBuffer<R>,
    origin: [R; 3],
    cells: u32,
    cell_size: R,
    width: R,
    faces: u8,
) -> crate::Result<()> {
    if !cell_size.is_finite() || cell_size <= R::ZERO {
        return Err(crate::Error::InvalidCellSize {
            value: f64::from(cell_size.as_f32()),
        });
    }
    if !width.is_finite() || width < R::ZERO {
        return Err(crate::Error::InvalidCellSize {
            value: f64::from(width.as_f32()),
        });
    }
    if cells < 2 {
        return Err(crate::Error::GridTooSmall {
            size: [cells + 1; 3],
        });
    }
    if width == R::ZERO || faces == 0 {
        return Ok(());
    }

    let span = R::from_f64(f64::from(cells));
    let one = R::ONE;
    let inv = cell_size.recip();

    for position in &mut mesh.positions {
        for axis in 0..3usize {
            // The coordinate in this block's own cells.
            let c = (position[axis] - origin[axis]) * inv;

            if c < one && faces & face_bit(axis, 0) != 0 {
                position[axis] += (one - c) * width;
            } else if c > span - one && faces & face_bit(axis, 1) != 0 {
                position[axis] += (span - one - c) * width;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
