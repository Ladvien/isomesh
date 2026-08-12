//! Where a transition cell's vertices go, and the one identity the seam rests on.
//!
//! [`table`](super::table) says *which* edges are cut and how they link. This
//! says *where* the crossings are — and the whole of A-011b depends on one
//! property of that placement:
//!
//! > A crossing on a **half-resolution** edge must land exactly where the coarse
//! > neighbour's own Marching Cubes pass puts its vertex.
//!
//! Not "within a tolerance". Exactly. Two vertices a float apart leave a crack
//! the width of a rounding error, and A-013 measured what that costs: a seam of
//! unshared vertices that a renderer draws correctly and a collider reads as a
//! hole (M-69). The weld can close a seam it can *see*; it cannot invent a shared
//! vertex where the two sides disagree in the last bit.
//!
//! Three facts make the identity hold, and each was measured rather than assumed:
//!
//! - **The endpoints coincide.** A half-resolution edge joins two *corner*
//!   samples of the transition face, and M-70 measured that a level-`k` sample
//!   position is bit-identical to the level-0 position of the sample it sits on,
//!   at every spacing tried including `4/35`.
//! - **The values coincide**, because they are the same field evaluated at the
//!   same point.
//! - **The interpolation coincides**, because [`TransitionCell::crossing`] is
//!   `lo + (hi − lo)·t` with `t = a/(a − b)` — character for character what
//!   `marching_cubes`' `edge_position` computes, and the reason it is written out
//!   here rather than shared is that the two take their operands from different
//!   places. `the_half_resolution_crossings_are_the_coarse_neighbours_vertices`
//!   is what keeps them equal.

use crate::cube::{edge_crossing, is_inside};
use crate::vec3;
use crate::{Real, Sdf};

use super::table::{EDGE_COUNT, EDGE_SAMPLES, SAMPLE_COUNT};

/// One transition cell's nine sample values and their world positions.
///
/// The four half-resolution corners are *not* stored twice: they are samples
/// 0, 2, 6 and 8 of the same nine, which is exactly why the seam closes. See the
/// module docs.
#[derive(Clone, Copy, Debug)]
pub struct TransitionCell<R: Real> {
    /// Field value at each of the nine samples.
    pub value: [R; SAMPLE_COUNT],
    /// World position of each of the nine samples.
    pub position: [[R; 3]; SAMPLE_COUNT],
}

impl<R: Real> TransitionCell<R> {
    /// Sample a transition face.
    ///
    /// `origin` is the **grid's** origin — the world position of global sample
    /// `[0, 0, 0]` — and `base` is the global **fine** sample index of this
    /// face's sample 0. `step` is the fine spacing, so the face spans `2·step`
    /// on each in-plane axis, and `u`/`v` are those axes. The nine samples are
    /// laid out `u`-fastest, matching [`table`](super::table)'s numbering.
    ///
    /// # Why an index and not a face origin
    ///
    /// The first version of this took the face's own world origin and added
    /// local offsets to it. That is wrong, and the test caught it: at
    /// `h = 4/14` a half-resolution crossing came out at `y = -1.11e-16` where
    /// the coarse mesh had it at exactly `0`, because
    ///
    /// ```text
    /// (origin + h·i) + h    ≠    origin + h·(i + 1)
    /// ```
    ///
    /// in IEEE at a spacing that is not a power of two. A hairline difference in
    /// the last bit is not a rounding curiosity here — it is a **crack**, and one
    /// no weld can close, because the weld can only merge vertices it can see are
    /// the same and these two are not.
    ///
    /// Indexing from the grid origin is the same expression
    /// [`ChunkLayout::world_of_sample`](crate::chunk::ChunkLayout::world_of_sample)
    /// uses, so the coarse grid's `origin + (2h)·c` and this cell's
    /// `origin + h·(2c)` are bit-identical by M-70. See M-73.
    pub fn sample<S: Sdf<Scalar = R>>(
        sdf: &S,
        origin: [R; 3],
        step: R,
        base: [i64; 3],
        u: usize,
        v: usize,
    ) -> Self {
        let mut value = [R::ZERO; SAMPLE_COUNT];
        let mut position = [[R::ZERO; 3]; SAMPLE_COUNT];
        for s in 0..SAMPLE_COUNT {
            let mut index = base;
            index[u] += (s % 3) as i64;
            index[v] += (s / 3) as i64;
            let mut p = [R::ZERO; 3];
            for (axis, slot) in p.iter_mut().enumerate() {
                *slot = origin[axis] + step * R::from_f64(index[axis] as f64);
            }
            position[s] = p;
            value[s] = sdf.sample(p);
        }
        Self { value, position }
    }

    /// Whether the surface crosses this edge.
    #[must_use]
    pub fn is_cut(&self, edge: u8) -> bool {
        let [a, b] = EDGE_SAMPLES[edge as usize];
        is_inside(self.value[a as usize]) != is_inside(self.value[b as usize])
    }

    /// Where the surface crosses `edge`, or `None` if it does not.
    ///
    /// `lo + (hi − lo)·t` with `t = a/(a − b)`, which is Marching Cubes' own
    /// placement. On a cut edge exactly one endpoint is strictly negative and the
    /// other is `>= 0`, so `a − b` is never zero and no epsilon guard is needed —
    /// and an epsilon here would snap resolvable crossings to the midpoint, which
    /// is precisely the sub-voxel detail A-014 exists to keep.
    #[must_use]
    pub fn crossing(&self, edge: u8) -> Option<[R; 3]> {
        if !self.is_cut(edge) {
            return None;
        }
        let [lo, hi] = EDGE_SAMPLES[edge as usize];
        let (a, b) = (self.value[lo as usize], self.value[hi as usize]);
        let t = edge_crossing(a, b);
        let (lo_pos, hi_pos) = (self.position[lo as usize], self.position[hi as usize]);
        Some([
            lo_pos[0] + (hi_pos[0] - lo_pos[0]) * t,
            lo_pos[1] + (hi_pos[1] - lo_pos[1]) * t,
            lo_pos[2] + (hi_pos[2] - lo_pos[2]) * t,
        ])
    }

    /// Every crossing, in edge order, with its edge index.
    pub fn crossings(&self) -> impl Iterator<Item = (u8, [R; 3])> + '_ {
        (0..EDGE_COUNT as u8).filter_map(|e| self.crossing(e).map(|p| (e, p)))
    }

    /// The case index this cell presents to [`table`](super::table).
    #[must_use]
    pub fn case(&self) -> u16 {
        let mut case = 0u16;
        for (s, v) in self.value.iter().enumerate() {
            if is_inside(*v) {
                case |= 1 << s;
            }
        }
        case
    }

    /// The centroid of this cell's crossings, or `None` when the surface misses
    /// it.
    ///
    /// Not a vertex rule — a transition cell's vertices are its crossings. This
    /// is here for A-011c, which needs a point inside the cell to fan a long
    /// cycle from, the same way `marching_cubes` fans one from a centroid rather
    /// than from one of its own vertices (A-015).
    #[must_use]
    pub fn centroid(&self) -> Option<[R; 3]> {
        let mut sum = [R::ZERO; 3];
        let mut n = 0u32;
        for (_, p) in self.crossings() {
            for (axis, slot) in sum.iter_mut().enumerate() {
                *slot += p[axis];
            }
            n += 1;
        }
        if n == 0 {
            return None;
        }
        Some(vec3::scale(sum, R::from_f64(f64::from(n)).recip()))
    }
}

#[cfg(test)]
mod tests;
