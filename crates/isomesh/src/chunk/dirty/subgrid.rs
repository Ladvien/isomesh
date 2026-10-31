//! [`mark_edit_subgrid`]: which chunks an edit changed, for the extractor that
//! reads the field between the lattice corners.
//!
//! [`mark_edit`](super::mark_edit) compares two fields at the eight corners of
//! each cell and calls a cell a surface cell when its corners disagree. For
//! Marching Cubes that is the extractor's whole input. Subgrid Marching
//! Tetrahedra reads more, and Baktash, Gillespie & Crane say what
//! (`10.48550/arXiv.2606.00454` §1):
//!
//! > We replace 0-dimensional sampling (evaluate `f` at each grid node), with
//! > 1-dimensional root finding (find all zeros of `f` along each grid edge).
//!
//! A feature thinner than a cell leaves all eight corners on one side and
//! crosses the surface twice along an edge, so corner signs cannot tell whether
//! the cell's mesh changed (issue #11). This module asks the question of the
//! points the extractor reads.

use alloc::vec::Vec;

use super::{DirtySet, EditReport, tally};
use crate::chunk::ChunkLayout;
use crate::cube::{corner_offset, is_inside};
use crate::marching_tetrahedra::table::{TET_COUNT, TET_EDGE_COUNT, tet_edge_corners};
use crate::subgrid::roots::{edge_point, edge_step, sample_parameter};
use crate::{Real, Sdf};

/// A probe point's flags: the two fields differ there, and which side of the
/// surface each puts it on.
const CHANGED: u8 = 1;
const INSIDE_BEFORE: u8 = 2;
const INSIDE_AFTER: u8 = 4;

/// What one edge, or one cell, showed over all its probe points.
const SHOWED_CHANGE: u8 = 1;
const SHOWED_SIGN: u8 = 2;
const CROSSED_BEFORE: u8 = 4;
const CROSSED_AFTER: u8 = 8;

/// A tetrahedron edge of a cell, as the layered sampler addresses it.
///
/// The edge runs from cube corner `lo` to cube corner `hi`, where `lo`'s bits
/// are a subset of `hi`'s (`TETS` orders every tetrahedron by inclusion), so the
/// step is `0` or `1` on each axis.
#[derive(Clone, Copy, Default)]
struct EdgeSpec {
    lo: [usize; 3],
    hi: [usize; 3],
    /// Which of the plane's three in-plane directions (x, y, xy), or of the
    /// slab's four through-plane ones (z, xz, yz, xyz).
    direction: usize,
    /// Whether the edge leaves its plane: those live in the slab between two
    /// planes, the others in a plane.
    through_z: bool,
}

/// Every distinct tetrahedron edge of a cell, and how many there are.
///
/// Read off the extractor's own tables rather than written out, so a change to
/// the decomposition changes this with it. Kuhn's split gives 19: the 12 cube
/// edges, one diagonal per face (6) and the main diagonal.
fn cell_edges() -> ([EdgeSpec; TET_COUNT * TET_EDGE_COUNT], usize) {
    let mut specs = [EdgeSpec::default(); TET_COUNT * TET_EDGE_COUNT];
    let mut seen = [[false; 8]; 8];
    let mut count = 0usize;
    for tet in 0..TET_COUNT {
        for edge in 0..TET_EDGE_COUNT {
            let [lo, hi] = tet_edge_corners(tet, edge);
            let (l, h) = (usize::from(lo), usize::from(hi));
            if seen[l][h] {
                continue;
            }
            seen[l][h] = true;
            let a = corner_offset(lo);
            let b = corner_offset(hi);
            // The bits that differ are the axes the edge runs along. Bit 2 is
            // z, so a difference of 4 or more leaves the plane.
            let step = lo ^ hi;
            let through_z = step >= 4;
            specs[count] = EdgeSpec {
                lo: [a[0] as usize, a[1] as usize, a[2] as usize],
                hi: [b[0] as usize, b[1] as usize, b[2] as usize],
                direction: usize::from(if through_z { step - 4 } else { step - 1 }),
                through_z,
            };
            count += 1;
        }
    }
    (specs, count)
}

/// The two lattice planes bounding one layer of cells, and the interior points
/// of the edges that stay in a plane.
struct Plane {
    /// One flag byte per lattice point.
    corner: Vec<u8>,
    /// Flags for the `samples - 1` interior points of each in-plane edge,
    /// indexed by the edge's lower lattice point: direction x, y, xy.
    along: [Vec<u8>; 3],
}

/// The interior points of the edges that run between two adjacent planes, by
/// their lower lattice point: direction z, xz, yz, xyz.
struct Slab {
    along: [Vec<u8>; 4],
}

/// The two fields, the region, and the lattice, bundled so the plane and slab
/// fills read as what they do.
struct Probes<'a, R: Real, A, B> {
    layout: &'a ChunkLayout<R>,
    before: &'a A,
    after: &'a B,
    min_cell: [i64; 3],
    /// Lattice points per axis in the region.
    extent: [usize; 3],
    samples: u32,
    step: R,
}

impl<R, A, B> Probes<'_, R, A, B>
where
    R: Real,
    A: Sdf<Scalar = R>,
    B: Sdf<Scalar = R>,
{
    /// The interior points per edge.
    fn interior(&self) -> usize {
        self.samples as usize - 1
    }

    /// World position of region lattice point `(x, y, z)`.
    fn position(&self, x: usize, y: usize, z: usize) -> [R; 3] {
        self.layout.world_of_sample([
            self.min_cell[0] + x as i64,
            self.min_cell[1] + y as i64,
            self.min_cell[2] + z as i64,
        ])
    }

    /// Both fields at `p`, as flags.
    ///
    /// Bit comparison, not `!=`, as in [`mark_edit`](super::mark_edit): a sign
    /// flip on a zero is a real change, since zero is outside.
    fn flags(&self, p: [R; 3]) -> u8 {
        let a = self.before.sample(p);
        let b = self.after.sample(p);
        let mut f = 0;
        if a.total_cmp(&b) != core::cmp::Ordering::Equal {
            f |= CHANGED;
        }
        if is_inside(a) {
            f |= INSIDE_BEFORE;
        }
        if is_inside(b) {
            f |= INSIDE_AFTER;
        }
        f
    }

    fn new_plane(&self) -> Plane {
        let n = self.extent[0] * self.extent[1];
        let m = n * self.interior();
        Plane {
            corner: alloc::vec![0; n],
            along: [
                alloc::vec![0; m],
                alloc::vec![0; m],
                alloc::vec![0; m],
            ],
        }
    }

    fn new_slab(&self) -> Slab {
        let m = self.extent[0] * self.extent[1] * self.interior();
        Slab {
            along: [
                alloc::vec![0; m],
                alloc::vec![0; m],
                alloc::vec![0; m],
                alloc::vec![0; m],
            ],
        }
    }

    /// Fill `out` with the points of region lattice plane `z`: every lattice
    /// point, and the interior points of the x, y and xy edges inside the plane.
    fn fill_plane(&self, out: &mut Plane, z: usize) {
        let [ex, ey, _] = self.extent;
        let interior = self.interior();
        for y in 0..ey {
            for x in 0..ex {
                out.corner[y * ex + x] = self.flags(self.position(x, y, z));
            }
        }
        // (dx, dy) for direction x, y, xy.
        for (direction, (dx, dy)) in [(1, 0), (0, 1), (1, 1)].into_iter().enumerate() {
            for y in 0..ey.saturating_sub(dy) {
                for x in 0..ex.saturating_sub(dx) {
                    let from = self.position(x, y, z);
                    let to = self.position(x + dx, y + dy, z);
                    let slot = (y * ex + x) * interior;
                    for k in 1..self.samples {
                        let t = sample_parameter(k, self.samples, self.step);
                        out.along[direction][slot + k as usize - 1] =
                            self.flags(edge_point(from, to, t));
                    }
                }
            }
        }
    }

    /// Fill `out` with the interior points of the edges between plane `z` and
    /// plane `z + 1`: z, xz, yz and xyz.
    fn fill_slab(&self, out: &mut Slab, z: usize) {
        let [ex, ey, _] = self.extent;
        let interior = self.interior();
        for (direction, (dx, dy)) in [(0, 0), (1, 0), (0, 1), (1, 1)].into_iter().enumerate() {
            for y in 0..ey.saturating_sub(dy) {
                for x in 0..ex.saturating_sub(dx) {
                    let from = self.position(x, y, z);
                    let to = self.position(x + dx, y + dy, z + 1);
                    let slot = (y * ex + x) * interior;
                    for k in 1..self.samples {
                        let t = sample_parameter(k, self.samples, self.step);
                        out.along[direction][slot + k as usize - 1] =
                            self.flags(edge_point(from, to, t));
                    }
                }
            }
        }
    }
}

/// What the points of one edge show: any change, any classification change,
/// and whether either field's classification flips between neighbouring points
/// (the extractor finds a root there).
#[inline]
fn scan_edge(first: u8, inner: &[u8], last: u8) -> u8 {
    let point = |f: u8| (f & CHANGED) | ((((f >> 1) ^ (f >> 2)) & 1) << 1);
    let mut shown = point(first);
    let mut prev = first;
    for &f in inner.iter().chain(core::iter::once(&last)) {
        shown |= point(f);
        // INSIDE_BEFORE is bit 1 and INSIDE_AFTER bit 2; the crossings are bits
        // 2 and 3.
        shown |= ((prev ^ f) & (INSIDE_BEFORE | INSIDE_AFTER)) << 1;
        prev = f;
    }
    shown
}

/// Compare two fields over a region at the points **Subgrid Marching
/// Tetrahedra reads**, mark the chunks that changed, and report how much of the
/// region moved. The counterpart of [`mark_edit`](super::mark_edit) for that
/// extractor.
///
/// `samples` is the extractor's 1D resolution
/// ([`SubgridMarchingTetrahedra::new`](crate::subgrid::extract::SubgridMarchingTetrahedra::new)):
/// every one of a cell's 19 tetrahedron edges is cut into that many intervals,
/// and the field is compared at each cut. A cell is then:
///
/// - **value changed** when any such point differs between the fields;
/// - a **surface cell** before (after) when the field classifies two
///   neighbouring points of some edge on opposite sides, which is exactly when
///   the extractor finds a root and emits triangles;
/// - **output changed**, and its owner chunk dirty, when it is value-changed and
///   a surface cell before or after.
///
/// The counters of the returned [`EditReport`] mean what they mean in
/// [`mark_edit`](super::mark_edit), with "corner" read as "probe point".
/// **`samples == 1` is `mark_edit`'s rule exactly**: the cube edges connect the
/// eight corners, so some edge changes class if and only if the corners
/// disagree, and the report is equal field by field.
///
/// # What it is exact for, and what it is not
///
/// A feature thinner than a cell leaves every corner outside and still puts two
/// roots on an edge; the corner rule cannot see it and this one does. It is
/// exact with respect to the **probe points**. Roots are then refined by
/// bisection, and vertex normals come from the field's gradient at each root,
/// and neither reads a probe. An edit that lies strictly between two probes and
/// moves no probe value is invisible here, as it is to the extractor's own root
/// search (`roots`: a pair of roots closer together than the sample spacing is
/// invisible), and Baktash et al. (§1.3) price it the same way: "we are no worse
/// off than classic marching".
///
/// Probe positions use the extractor's own expression along each edge, between
/// the lattice points [`ChunkLayout::world_of_sample`] gives, as
/// [`mark_edit`](super::mark_edit) does. The extractor computes lattice points as
/// `origin + cell_size * local`, which is bit-equal for a power-of-two cell size
/// and can differ in the last place otherwise (M-32); only a probe within an ulp
/// of a zero could land on the other side.
///
/// # Cost
///
/// `1 + 7 (samples - 1)` evaluations of each field per cell, amortised: the
/// lattice point, and the `samples - 1` interior points of one edge in each of
/// seven directions (x, y, z, xy, xz, yz, xyz), each point evaluated once
/// however many cells share it. At `samples == 1` that is the one evaluation per
/// cell `mark_edit` makes. The region is the caller's, so the cost is theirs to
/// bound: pass the edit's bounding box rather than the volume's.
///
/// Memory is one layer of cells at a time (two lattice planes and the slab
/// between them, one byte per point), not a grid the size of the region.
///
/// # Errors
///
/// [`Error::InvalidCellSize`](crate::Error::InvalidCellSize) if `samples` is
/// `0`: the extractor cannot be built, so no rule can match it, and this is its
/// constructor's error. Checked before the region, so an inverted region does
/// not hide a bad argument.
///
/// [`Error::IndexSpaceExhausted`](crate::Error::IndexSpaceExhausted) if the
/// probe grid, at most `1 + 7 (samples - 1)` points per lattice point of the
/// region, exceeds the `u32` sample space: reported rather than attempted, as in
/// [`mark_edit`](super::mark_edit).
pub fn mark_edit_subgrid<R, A, B>(
    layout: &ChunkLayout<R>,
    before: &A,
    after: &B,
    min_cell: [i64; 3],
    max_cell: [i64; 3],
    samples: u32,
    dirty: &mut DirtySet,
) -> crate::Result<EditReport>
where
    R: Real,
    A: Sdf<Scalar = R>,
    B: Sdf<Scalar = R>,
{
    if samples == 0 {
        return Err(crate::Error::InvalidCellSize { value: 0.0 });
    }
    let mut report = EditReport::default();
    for axis in 0..3 {
        if max_cell[axis] < min_cell[axis] {
            return Ok(report);
        }
    }

    // Extents in u64, where even `i64::MAX - i64::MIN` fits: `max >= min` held
    // above, so the wrapping difference is the true difference.
    let mut lattice = 1u64;
    let mut extent = [0usize; 3];
    for axis in 0..3 {
        let planes = (max_cell[axis].wrapping_sub(min_cell[axis]) as u64).saturating_add(2);
        lattice = lattice.saturating_mul(planes);
        // Sound to narrow only because the gate below refuses anything past the
        // u32 sample space before `extent` is ever read.
        extent[axis] = planes as usize;
    }
    let probes = lattice.saturating_mul(1 + 7 * u64::from(samples - 1));
    if probes > u64::from(u32::MAX) {
        return Err(crate::Error::IndexSpaceExhausted { needed: probes });
    }

    let grid = Probes {
        layout,
        before,
        after,
        min_cell,
        extent,
        samples,
        step: edge_step::<R>(samples),
    };
    let interior = grid.interior();
    let [ex, ey, ez] = extent;
    let (edges, edge_count) = cell_edges();
    let edges = &edges[..edge_count];

    let mut low = grid.new_plane();
    let mut high = grid.new_plane();
    let mut slab = grid.new_slab();
    grid.fill_plane(&mut low, 0);

    // Distinct owners this call dirtied, so the count is this edit's and not the
    // caller's whole queue.
    let mut dirtied = DirtySet::new();
    let mut region_chunks = DirtySet::new();
    for cz in 0..ez - 1 {
        grid.fill_plane(&mut high, cz + 1);
        grid.fill_slab(&mut slab, cz);

        for cy in 0..ey - 1 {
            for cx in 0..ex - 1 {
                let mut shown = 0u8;
                for e in edges {
                    let a = (cy + e.lo[1]) * ex + cx + e.lo[0];
                    let b = (cy + e.hi[1]) * ex + cx + e.hi[0];
                    let slot = a * interior;
                    shown |= if e.through_z {
                        scan_edge(
                            low.corner[a],
                            &slab.along[e.direction][slot..slot + interior],
                            high.corner[b],
                        )
                    } else {
                        // Both ends are in one plane: the layer's low plane for
                        // an edge on the cell's bottom face, its high plane for
                        // one on the top face.
                        let plane = if e.lo[2] == 0 { &low } else { &high };
                        scan_edge(
                            plane.corner[a],
                            &plane.along[e.direction][slot..slot + interior],
                            plane.corner[b],
                        )
                    };
                }

                let cell = [
                    min_cell[0] + cx as i64,
                    min_cell[1] + cy as i64,
                    min_cell[2] + cz as i64,
                ];
                // A cell is meshed by the chunk that owns it, and cell `k` is
                // owned by the chunk owning sample `k`.
                let (owner, _) = layout.local_sample(cell);
                region_chunks.insert(owner);

                if tally(
                    &mut report,
                    shown & SHOWED_CHANGE != 0,
                    shown & SHOWED_SIGN != 0,
                    shown & CROSSED_BEFORE != 0,
                    shown & CROSSED_AFTER != 0,
                ) {
                    dirtied.insert(owner);
                }
            }
        }
        core::mem::swap(&mut low, &mut high);
    }

    report.region_chunks = region_chunks.len() as u64;
    report.dirty_chunks = dirtied.len() as u64;
    for id in dirtied.iter() {
        dirty.insert(id);
    }
    Ok(report)
}

#[cfg(test)]
mod tests;
