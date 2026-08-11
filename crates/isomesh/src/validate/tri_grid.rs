//! A uniform spatial index over triangles, and the point-triangle distance it
//! accelerates.
//!
//! # Layout
//!
//! Compressed-sparse-row: `offsets` holds one entry per cell plus a terminator,
//! and `items` holds the triangle indices of every cell back to back. Cell `c`
//! owns `items[offsets[c]..offsets[c + 1]]`.
//!
//! There is no hash map, for the reason [the parent module](super) gives, and
//! there is no sort either — which makes the determinism argument stronger here
//! than it is there:
//!
//! 1. `offsets` is a prefix sum over a count vector, a pure function of the
//!    triangle geometry.
//! 2. The fill pass walks triangles in ascending index order, so every cell's
//!    item list comes out ascending. No sort, no tie-break, and no
//!    `sort_unstable` instability to reason about.
//! 3. The query visits cells in a fixed order and accumulates a **minimum**. A
//!    minimum is a value rather than an index, so the answer does not depend on
//!    visit order at all.
//!
//! # Cost
//!
//! Build is `O(F · cells-per-triangle)` with two passes and no comparison sort.
//! Memory is `4·(cells + 1) + 4·entries` bytes. A query is `O(1)` in the
//! occupancy of the shells it visits, and terminates as soon as the best
//! distance found is provably closer than anything unvisited.

use alloc::vec;
use alloc::vec::Vec;

use crate::Real;
use crate::vec3;

/// Broadphase blow-up guard, in grid cells touched per triangle.
///
/// With a `cell_size` matching the extraction grid, a triangle spans at most a
/// couple of cells per axis. Wildly more means the caller passed a spacing that
/// does not describe this mesh, and the grid would grow until it exhausted
/// memory. Better to say so.
///
/// Shared with [`self_intersections`](super::self_intersections) so that "this
/// spacing does not describe this mesh" has one definition rather than two that
/// can drift apart.
pub(crate) const MAX_CELLS_PER_TRIANGLE: usize = 512;

/// Total cells a grid may allocate.
///
/// Two jobs. It bounds memory against a caller who passes a spacing far finer
/// than the mesh, and it keeps every axis index well below `2^24`, which is
/// where [`cell_of`]'s narrowing through `f32` would start to lose integers.
pub(crate) const MAX_TOTAL_CELLS: u64 = 1 << 22;

/// Floor to a grid coordinate. A non-finite input lands in cell zero and is
/// filtered out by the caller's finiteness check before it can matter.
#[inline]
pub(crate) fn cell_of<R: Real>(scaled: R) -> i64 {
    let f = scaled.floor();
    if f.is_finite() { f.as_f32() as i64 } else { 0 }
}

/// Squared distance from a point to a triangle.
///
/// The region/Voronoi form: classify the point against the triangle's seven
/// features (three vertices, three edges, the face interior) and project onto
/// whichever one owns it.
///
/// Squared throughout, so the inner loop takes no square root — the caller takes
/// one per query rather than one per triangle.
///
/// Reference: Christer Ericson, *Real-Time Collision Detection*, Morgan
/// Kaufmann 2004, §5.1.5 "Closest Point on Triangle to Point".
///
/// # Robustness
///
/// The face branch is reached only when all three of `va`, `vb`, `vc` are
/// strictly positive, which means every edge test failed. So `denom > 0` there
/// and the barycentric coordinates cannot blow up, even on a needle triangle.
/// Exactly collinear input gives `va == vb == vc == 0` and fires a vertex or
/// edge branch first, so the face branch never sees it.
pub(crate) fn point_triangle_distance_squared<R: Real>(
    p: [R; 3],
    a: [R; 3],
    b: [R; 3],
    c: [R; 3],
) -> R {
    let ab = vec3::sub(b, a);
    let ac = vec3::sub(c, a);
    let ap = vec3::sub(p, a);

    let d1 = vec3::dot(ab, ap);
    let d2 = vec3::dot(ac, ap);
    if d1 <= R::ZERO && d2 <= R::ZERO {
        return vec3::length_squared(ap); // vertex region A
    }

    let bp = vec3::sub(p, b);
    let d3 = vec3::dot(ab, bp);
    let d4 = vec3::dot(ac, bp);
    if d3 >= R::ZERO && d4 <= d3 {
        return vec3::length_squared(bp); // vertex region B
    }

    let vc_ = d1 * d4 - d3 * d2;
    if vc_ <= R::ZERO && d1 >= R::ZERO && d3 <= R::ZERO {
        let v = d1 / (d1 - d3); // edge region AB
        return vec3::length_squared(vec3::sub(ap, vec3::scale(ab, v)));
    }

    let cp = vec3::sub(p, c);
    let d5 = vec3::dot(ab, cp);
    let d6 = vec3::dot(ac, cp);
    if d6 >= R::ZERO && d5 <= d6 {
        return vec3::length_squared(cp); // vertex region C
    }

    let vb_ = d5 * d2 - d1 * d6;
    if vb_ <= R::ZERO && d2 >= R::ZERO && d6 <= R::ZERO {
        let w = d2 / (d2 - d6); // edge region AC
        return vec3::length_squared(vec3::sub(ap, vec3::scale(ac, w)));
    }

    let va_ = d3 * d6 - d5 * d4;
    if va_ <= R::ZERO && (d4 - d3) >= R::ZERO && (d5 - d6) >= R::ZERO {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6)); // edge region BC
        let bc = vec3::sub(c, b);
        return vec3::length_squared(vec3::sub(bp, vec3::scale(bc, w)));
    }

    // Face region.
    let denom = (va_ + vb_ + vc_).recip();
    let v = vb_ * denom;
    let w = vc_ * denom;
    let closest = vec3::sub(vec3::sub(ap, vec3::scale(ab, v)), vec3::scale(ac, w));
    vec3::length_squared(closest)
}

/// A uniform grid binning triangles by their axis-aligned bounds.
#[derive(Clone, Debug)]
pub(crate) struct TriangleGrid<R: Real> {
    /// World position of the `[0, 0, 0]` cell's minimum corner.
    lo: [R; 3],
    cell_size: R,
    inv_cell: R,
    dims: [u32; 3],
    /// One entry per cell, plus a terminator. `offsets[c]..offsets[c + 1]`.
    offsets: Vec<u32>,
    /// Triangle indices, grouped by cell, ascending within each cell.
    items: Vec<u32>,
}

impl<R: Real> TriangleGrid<R> {
    /// Bin `tris` by the cells their bounds touch.
    ///
    /// `tris` must already be filtered: every index in range, every coordinate
    /// finite. The grid does not re-check, because the caller has to walk the
    /// triangles to filter them anyway and doing it twice invites the two checks
    /// to disagree.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidCellSize`](crate::Error::InvalidCellSize) if `cell_size`
    /// is not finite and positive.
    ///
    /// [`Error::CellSizeMismatch`](crate::Error::CellSizeMismatch) if one
    /// triangle spans more than [`MAX_CELLS_PER_TRIANGLE`] cells, or if the grid
    /// as a whole would exceed [`MAX_TOTAL_CELLS`]. Both mean the spacing does
    /// not describe this mesh.
    pub(crate) fn build(
        positions: &[[R; 3]],
        tris: &[[u32; 3]],
        cell_size: f64,
    ) -> crate::Result<Self> {
        if !cell_size.is_finite() || cell_size <= 0.0 {
            return Err(crate::Error::InvalidCellSize { value: cell_size });
        }
        let h = R::from_f64(cell_size);
        let inv_cell = h.recip();

        // An empty grid is a legal grid: it answers every query with infinity.
        if tris.is_empty() {
            return Ok(Self {
                lo: [R::ZERO; 3],
                cell_size: h,
                inv_cell,
                dims: [1; 3],
                offsets: vec![0, 0],
                items: Vec::new(),
            });
        }

        let mut lo = [R::INFINITY; 3];
        let mut hi = [R::NEG_INFINITY; 3];
        for t in tris {
            for &i in t {
                let v = positions[i as usize];
                for a in 0..3 {
                    lo[a] = lo[a].min(v[a]);
                    hi[a] = hi[a].max(v[a]);
                }
            }
        }

        let mut dims = [1u32; 3];
        let mut total: u64 = 1;
        for a in 0..3 {
            let span = cell_of((hi[a] - lo[a]) * inv_cell) + 1;
            let d = span.clamp(1, i64::from(u32::MAX)) as u32;
            dims[a] = d;
            total = total.saturating_mul(u64::from(d));
        }
        if total > MAX_TOTAL_CELLS {
            return Err(crate::Error::CellSizeMismatch {
                triangle: u64::MAX,
                cells: u128::from(total),
                cell_size,
            });
        }

        let cells = (dims[0] as usize) * (dims[1] as usize) * (dims[2] as usize);
        let mut grid = Self {
            lo,
            cell_size: h,
            inv_cell,
            dims,
            offsets: vec![0; cells + 1],
            items: Vec::new(),
        };

        // Pass 1 — count. `offsets[c + 1]` accumulates the population of `c`, so
        // the prefix sum in pass 2 needs no shift.
        for (ti, t) in tris.iter().enumerate() {
            let (blo, bhi) = grid.bounds_cells(positions, *t);
            let span = span_cells(blo, bhi);
            if span > MAX_CELLS_PER_TRIANGLE as u128 {
                return Err(crate::Error::CellSizeMismatch {
                    triangle: ti as u64,
                    cells: span,
                    cell_size,
                });
            }
            for z in blo[2]..=bhi[2] {
                for y in blo[1]..=bhi[1] {
                    for x in blo[0]..=bhi[0] {
                        grid.offsets[grid_index(dims, x, y, z) + 1] += 1;
                    }
                }
            }
        }

        // Pass 2 — prefix sum.
        for c in 0..cells {
            grid.offsets[c + 1] += grid.offsets[c];
        }

        // Pass 3 — fill. Ascending triangle order gives ascending item lists.
        let entries = grid.offsets[cells] as usize;
        grid.items = vec![0; entries];
        let mut cursor = grid.offsets.clone();
        for (ti, t) in tris.iter().enumerate() {
            let (blo, bhi) = grid.bounds_cells(positions, *t);
            for z in blo[2]..=bhi[2] {
                for y in blo[1]..=bhi[1] {
                    for x in blo[0]..=bhi[0] {
                        let c = grid_index(dims, x, y, z);
                        grid.items[cursor[c] as usize] = ti as u32;
                        cursor[c] += 1;
                    }
                }
            }
        }

        Ok(grid)
    }

    /// The triangles registered in one cell.
    #[inline]
    pub(crate) fn cell(&self, x: u32, y: u32, z: u32) -> &[u32] {
        let c = grid_index(self.dims, x, y, z);
        &self.items[self.offsets[c] as usize..self.offsets[c + 1] as usize]
    }

    /// Squared distance from `q` to the nearest triangle, or
    /// [`Real::INFINITY`](Real::INFINITY) when the grid is empty.
    ///
    /// Visits cells in shells of increasing Chebyshev radius around `q`'s own
    /// cell and stops as soon as the best distance found is no greater than the
    /// distance from `q` to the nearest unexamined region.
    ///
    /// That bound is exact rather than heuristic. Every cell inside the examined
    /// box has been visited, and a triangle is registered in every cell its
    /// bounds touch — including the one containing its own closest point. So any
    /// triangle not yet seen has its closest point outside the box, and every
    /// point outside the box is at least `reach` from `q`. A box face that sits
    /// on the grid boundary contributes no bound at all, because there is
    /// nothing beyond it.
    pub(crate) fn nearest_distance_squared(
        &self,
        q: [R; 3],
        positions: &[[R; 3]],
        tris: &[[u32; 3]],
    ) -> R {
        if self.items.is_empty() {
            return R::INFINITY;
        }
        let c = [
            self.axis_cell(q[0], 0),
            self.axis_cell(q[1], 1),
            self.axis_cell(q[2], 2),
        ];

        let mut best = R::INFINITY;
        let mut r: u32 = 0;
        loop {
            let blo = [
                c[0].saturating_sub(r),
                c[1].saturating_sub(r),
                c[2].saturating_sub(r),
            ];
            let bhi = [
                (c[0].saturating_add(r)).min(self.dims[0] - 1),
                (c[1].saturating_add(r)).min(self.dims[1] - 1),
                (c[2].saturating_add(r)).min(self.dims[2] - 1),
            ];

            for z in blo[2]..=bhi[2] {
                for y in blo[1]..=bhi[1] {
                    for x in blo[0]..=bhi[0] {
                        // Shell only: anything closer than `r` was visited on an
                        // earlier pass.
                        if r > 0 && chebyshev(c, [x, y, z]) < r {
                            continue;
                        }
                        for &ti in self.cell(x, y, z) {
                            let t = tris[ti as usize];
                            let d = point_triangle_distance_squared(
                                q,
                                positions[t[0] as usize],
                                positions[t[1] as usize],
                                positions[t[2] as usize],
                            );
                            if d < best {
                                best = d;
                            }
                        }
                    }
                }
            }

            // Distance from `q` to the exterior of the examined box, ignoring
            // faces that lie on the grid boundary.
            let mut reach = R::INFINITY;
            let mut bounded = false;
            for a in 0..3 {
                if blo[a] > 0 {
                    let wall = self.lo[a] + self.cell_size * R::from_f64(f64::from(blo[a]));
                    reach = reach.min(q[a] - wall);
                    bounded = true;
                }
                if bhi[a] < self.dims[a] - 1 {
                    let wall = self.lo[a] + self.cell_size * R::from_f64(f64::from(bhi[a] + 1));
                    reach = reach.min(wall - q[a]);
                    bounded = true;
                }
            }
            if !bounded {
                return best; // the box covers the whole grid
            }
            let reach = reach.max(R::ZERO);
            if best <= reach * reach {
                return best;
            }
            r += 1;
        }
    }

    /// Cell index along one axis, clamped into the grid.
    #[inline]
    fn axis_cell(&self, v: R, a: usize) -> u32 {
        let i = cell_of((v - self.lo[a]) * self.inv_cell);
        i.clamp(0, i64::from(self.dims[a] - 1)) as u32
    }

    /// Inclusive cell bounds of one triangle.
    fn bounds_cells(&self, positions: &[[R; 3]], t: [u32; 3]) -> ([u32; 3], [u32; 3]) {
        let a = positions[t[0] as usize];
        let b = positions[t[1] as usize];
        let c = positions[t[2] as usize];
        let mut lo = [0u32; 3];
        let mut hi = [0u32; 3];
        for ax in 0..3 {
            lo[ax] = self.axis_cell(a[ax].min(b[ax]).min(c[ax]), ax);
            hi[ax] = self.axis_cell(a[ax].max(b[ax]).max(c[ax]), ax);
        }
        (lo, hi)
    }
}

/// Cells spanned by an inclusive cell-coordinate box.
#[inline]
fn span_cells(lo: [u32; 3], hi: [u32; 3]) -> u128 {
    u128::from(hi[0] - lo[0] + 1) * u128::from(hi[1] - lo[1] + 1) * u128::from(hi[2] - lo[2] + 1)
}

/// `x` varies fastest, matching the crate's index convention.
#[inline]
fn grid_index(dims: [u32; 3], x: u32, y: u32, z: u32) -> usize {
    (x as usize)
        + (y as usize) * (dims[0] as usize)
        + (z as usize) * (dims[0] as usize) * (dims[1] as usize)
}

/// Chebyshev (chessboard) distance between two cell coordinates.
#[inline]
fn chebyshev(a: [u32; 3], b: [u32; 3]) -> u32 {
    a[0].abs_diff(b[0])
        .max(a[1].abs_diff(b[1]))
        .max(a[2].abs_diff(b[2]))
}
