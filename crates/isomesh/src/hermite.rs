//! Hermite data: where the surface crosses each cell edge, **and which way it
//! faces there**.
//!
//! Marching Cubes and Surface Nets need only the crossing positions. Dual
//! contouring needs the normals too, and that difference is the entire reason it
//! can reproduce a sharp corner: a position says the surface passes through
//! here, a position *and* a normal says it passes through here *going that way*,
//! and three of those intersect at a point that a corner actually is.
//!
//! # Source
//!
//! Ju, Losasso, Schaefer & Warren, *"Dual Contouring of Hermite Data"*, SIGGRAPH
//! 2002, `10.1145/566570.566586`. §2.2 defines the input as "the pairs
//! `pᵢ, nᵢ` which "correspond to the intersections (and unit normals) of the
//! contour with the edges of the cube" — which is exactly
//! [`HermiteCrossing`].
//!
//! # Why the position comes from the same interpolation as everywhere else
//!
//! The crossing is placed by the same linear interpolation the crate's shared
//! `cube::edge_crossing` gives Marching Cubes and Surface Nets, rather than by a
//! bisection or Newton refinement. That is deliberate and it is
//! about measurement, not cost: if dual contouring used *better* crossings than
//! Surface Nets, then E-104's side-by-side would be comparing two changes at
//! once — different crossings *and* different vertex placement — and could not
//! attribute the corner it recovers to either. Same crossings in, so the only
//! variable left is what the algorithm does with them.
//!
//! A refined root is a legitimate accuracy improvement later, but it belongs to
//! every extractor at once or to none.

use crate::cube::{EDGE_CORNERS, EDGE_COUNT, corner_offset, edge_crossing, is_inside};
use crate::vec3;
use crate::{Real, Sdf};

/// One edge crossing: where the surface is, and which way it faces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HermiteCrossing<R: Real> {
    /// World position of the crossing.
    pub position: [R; 3],
    /// Unit surface normal there, pointing away from the solid.
    ///
    /// Taken from the field's own gradient rather than estimated from geometry.
    /// For an exact signed distance field this is the true normal; for a field
    /// that is not a distance function it is still the direction of steepest
    /// increase, which is what the vertex solve wants.
    pub normal: [R; 3],
}

/// The crossings on one cell's twelve edges.
///
/// Up to twelve, in the crate's own edge numbering, with a bitmask saying which
/// are present. Fixed-size rather than a `Vec`, because a dual contourer builds
/// one of these per surface cell, and allocating per cell is exactly the cost
/// the reusable-buffer rule exists to avoid.
#[derive(Clone, Copy, Debug)]
pub struct HermiteCell<R: Real> {
    crossings: [HermiteCrossing<R>; EDGE_COUNT],
    /// Bit `e` set when edge `e` is cut.
    mask: u16,
}

impl<R: Real> HermiteCell<R> {
    /// Sample the crossings of one cell.
    ///
    /// `corner_values` are the eight corner samples in this crate's corner
    /// order, `cell_origin` is the world position of corner 0, and `cell_size`
    /// is the edge length. The caller supplies the corner values because it has
    /// already sampled them — re-sampling here would double the field
    /// evaluations, which dominate extraction.
    ///
    /// Costs one [`Sdf::gradient`] per crossing and nothing else.
    pub fn from_corners<S: Sdf<Scalar = R>>(
        sdf: &S,
        corner_values: &[R; 8],
        cell_origin: [R; 3],
        cell_size: R,
    ) -> Self {
        let empty = HermiteCrossing {
            position: [R::ZERO; 3],
            normal: [R::ZERO; 3],
        };
        let mut crossings = [empty; EDGE_COUNT];
        let mut mask = 0u16;

        for (edge, [lo, hi]) in EDGE_CORNERS.iter().copied().enumerate() {
            let (a, b) = (corner_values[lo as usize], corner_values[hi as usize]);
            if is_inside(a) == is_inside(b) {
                continue;
            }

            let t = edge_crossing(a, b);
            let (lo_offset, hi_offset) = (corner_offset(lo), corner_offset(hi));
            let mut position = [R::ZERO; 3];
            for (axis, slot) in position.iter_mut().enumerate() {
                let from = R::from_f64(f64::from(lo_offset[axis]));
                let to = R::from_f64(f64::from(hi_offset[axis]));
                *slot = cell_origin[axis] + cell_size * (from + (to - from) * t);
            }

            let gradient = sdf.gradient(position);
            let length = vec3::length(gradient);
            // A zero gradient at a crossing means the field is degenerate there.
            // It cannot happen for an exact distance field, where |grad| is 1.
            debug_assert!(length > R::ZERO, "zero gradient at an edge crossing");
            let normal = vec3::scale(gradient, length.recip());

            crossings[edge] = HermiteCrossing { position, normal };
            mask |= 1 << edge;
        }

        Self { crossings, mask }
    }

    /// Build a cell from explicit `(edge, crossing)` pairs.
    ///
    /// Internal, and it exists for one reason: the dual-contouring vertex rule
    /// is a piece of linear algebra whose interesting cases — three orthogonal
    /// planes, a rank-1 flat region, a lattice rotation applied to both input
    /// and expected output — are stated as *planes*, not as a field and a grid.
    /// Reaching them through [`from_corners`](Self::from_corners) would mean
    /// hand-designing an SDF for each, and the test would then be exercising the
    /// SDF as much as the solve.
    ///
    /// Edges outside `0..EDGE_COUNT` are ignored rather than rejected: this
    /// takes a fixed test fixture, not user input.
    #[cfg(test)]
    pub(crate) fn from_crossings(pairs: &[(u8, HermiteCrossing<R>)]) -> Self {
        let mut cell = Self {
            crossings: [HermiteCrossing {
                position: [R::ZERO; 3],
                normal: [R::ZERO; 3],
            }; EDGE_COUNT],
            mask: 0,
        };
        for &(edge, crossing) in pairs {
            if (edge as usize) < EDGE_COUNT {
                cell.crossings[edge as usize] = crossing;
                cell.mask |= 1 << edge;
            }
        }
        cell
    }

    /// How many edges the surface crosses.
    #[must_use]
    pub fn len(&self) -> usize {
        self.mask.count_ones() as usize
    }

    /// `true` when the surface misses this cell entirely.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.mask == 0
    }

    /// Whether edge `edge` is cut.
    #[must_use]
    pub fn contains(&self, edge: u8) -> bool {
        edge < EDGE_COUNT as u8 && self.mask & (1 << edge) != 0
    }

    /// The crossing on one edge, if there is one.
    #[must_use]
    pub fn get(&self, edge: u8) -> Option<&HermiteCrossing<R>> {
        if self.contains(edge) {
            Some(&self.crossings[edge as usize])
        } else {
            None
        }
    }

    /// Every crossing, in edge order.
    ///
    /// Edge order rather than insertion order, so the sequence is a function of
    /// the cell's corner signs alone — which is what keeps a vertex solve built
    /// on it deterministic.
    pub fn iter(&self) -> impl Iterator<Item = &HermiteCrossing<R>> {
        (0..EDGE_COUNT as u8).filter_map(|edge| self.get(edge))
    }

    /// The average of the crossing positions.
    ///
    /// This is precisely the vertex Surface Nets would place, and dual
    /// contouring's starting point: the QEF is minimised subject to staying near
    /// it, and it is the answer when the normals carry no directional
    /// information at all. Returns `None` for an empty cell.
    #[must_use]
    pub fn centroid(&self) -> Option<[R; 3]> {
        if self.is_empty() {
            return None;
        }
        let mut sum = [R::ZERO; 3];
        for crossing in self.iter() {
            for (axis, slot) in sum.iter_mut().enumerate() {
                *slot += crossing.position[axis];
            }
        }
        let inverse = R::from_f64(self.len() as f64).recip();
        Some([sum[0] * inverse, sum[1] * inverse, sum[2] * inverse])
    }
}

#[cfg(test)]
mod tests;
