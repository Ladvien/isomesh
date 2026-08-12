//! Surface Nets — the dual method.
//!
//! Marching Cubes puts vertices *on grid edges* and triangulates within each
//! cell. Surface Nets does the opposite: one vertex per cell that the surface
//! passes through, and one quad per grid edge that the surface crosses, joining
//! the four cells around that edge. Hence "dual" — vertices and faces swap roles.
//!
//! The consequences are worth stating precisely, because two of the three
//! things usually claimed for Surface Nets turn out to be measurably false:
//!
//! - **It does *not* produce fewer triangles.** Measured against Marching Cubes
//!   on four fields, the counts differ by exactly `2χ` — four triangles on a
//!   sphere, zero on a torus. That is forced by Euler, not luck: Marching Cubes
//!   places one vertex per crossed grid edge and Surface Nets emits two
//!   triangles per crossed grid edge, and any closed triangulated surface has
//!   `F = 2V − 2χ`. See `triangle_counts_track_marching_cubes_up_to_two_chi`.
//! - **Quads, so the connectivity is regular** — but not uniformly degree four
//!   once triangulated. Measured max degree on a sphere: 10, against Marching
//!   Cubes' 9.
//! - **Rounded corners, genuinely.** A vertex placed at the centroid of its
//!   cell's edge crossings cannot sit on a sharp feature, because an average of
//!   points on a corner's two faces lands between them. Measured on `box_exact`
//!   at 33³: the nearest vertex to the corner `(1,1,1)` is 1.15 cells away. That
//!   is not a defect to fix here — it is exactly what dual contouring changes,
//!   and it is the whole point of E-104.
//!
//!   Worth knowing before drawing conclusions from that number: Marching Cubes
//!   does *worse* on the same measurement (1.41 cells), because `box_exact` is
//!   exactly zero across its entire boundary and a grid plane on a box face
//!   therefore classifies as outside. See
//!   `neither_method_reaches_a_box_corner_and_the_reason_is_the_grid`.
//!
//! # The structural limit
//!
//! One vertex per cell means that where two sheets of the surface pass through
//! the same cell, they are forced to share a vertex and the result is
//! **non-manifold**. Measured: 48 non-manifold edges on the capped gyroid at
//! 49³, 15 on `fbm_terrain` at 33³. The literature review calls this out as "DC's
//! actual structural defect", fixed architecturally by vertex splitting rather
//! than by patching — which is A-010.
//!
//! # On the absence of published timings
//!
//! `docs/research/2026-08-11-meshing-speed-analysis.md` records that Surface
//! Nets and greedy meshing **have no credible published timings at all**,
//! despite being the two things game engines actually ship. Every number this
//! crate measures for it is therefore a reference rather than a comparison, and
//! should be reported with its hardware and grid attached.

use crate::cube::{EDGE_CORNERS, corner_offset, edge_crossing, is_inside};
use crate::dual::{CellVertices, DualMesher, VertexRule};
use crate::{MeshSink, Real, Sdf, Shape3};

/// Surface Nets over a sampled grid.
///
/// Owns its scratch buffers, for the same reason
/// [`MarchingCubes`](crate::marching_cubes::MarchingCubes) does: the real workload re-meshes
/// thousands of chunks and allocation dominates.
///
/// The topology — one vertex per crossed cell, one quad per crossed edge — is
/// shared with [`DualContouring`](crate::dual_contouring::DualContouring) and lives in
/// `crate::dual`. This type is that engine plus [`Centroid`], and the centroid
/// is the whole of what makes it Surface Nets rather than dual contouring.
///
/// # Example
///
/// ```
/// use isomesh::{MeshBuffer, RuntimeShape3};
/// use isomesh::fields::Sphere;
/// use isomesh::surface_nets::SurfaceNets;
///
/// let mut sn = SurfaceNets::<f32>::new();
/// let mut out = MeshBuffer::<f32>::new();
/// let shape = RuntimeShape3::new([33; 3])?;
/// sn.extract(&Sphere::<f32>::canonical(), &shape, [-2.0; 3], 0.125, &mut out)?;
///
/// assert!(out.triangle_count() > 0);
/// # Ok::<(), isomesh::Error>(())
/// ```
#[derive(Debug)]
pub struct SurfaceNets<R: Real> {
    mesher: DualMesher<R>,
}

/// The Surface Nets vertex rule: the centroid of the cell's edge crossings.
///
/// Cheap, and **always inside the cell** — so the mesh cannot self-intersect
/// through bad placement, which is a guarantee dual contouring gives up and
/// A-009 is the ticket that measures the cost of getting back.
///
/// It cannot represent a sharp feature, and not by accident: an average of
/// points lying on a corner's two faces lands between them, never on the corner.
/// Measured on `box_exact` at 33³, the nearest vertex to `(1,1,1)` is 1.15 cells
/// away. That is the gap [`DualContouring`](crate::dual_contouring::DualContouring) closes.
#[derive(Clone, Copy, Debug, Default)]
pub struct Centroid;

impl<R: Real> VertexRule<R> for Centroid {
    /// One vertex owning the whole cell — Surface Nets does not separate sheets,
    /// which is exactly the limitation A-010 measured and
    /// [`ManifoldDualContouring`](crate::manifold_dual_contouring::ManifoldDualContouring)
    /// lifts.
    fn place<S: Sdf<Scalar = R>>(
        &self,
        _sdf: &S,
        corner: &[R; 8],
        base: [u32; 3],
        origin: [R; 3],
        cell_size: R,
        out: &mut CellVertices<R>,
    ) {
        let mut sum = [R::ZERO; 3];
        let mut crossings = 0u32;
        for [lo, hi] in EDGE_CORNERS {
            let (a, b) = (corner[lo as usize], corner[hi as usize]);
            if is_inside(a) == is_inside(b) {
                continue;
            }
            let t = edge_crossing(a, b);
            let (lo_o, hi_o) = (corner_offset(lo), corner_offset(hi));
            for axis in 0..3 {
                let from = R::from_f64(f64::from(lo_o[axis]));
                let to = R::from_f64(f64::from(hi_o[axis]));
                sum[axis] += from + (to - from) * t;
            }
            crossings += 1;
        }
        // A cell with both signs present always has at least one cut edge: the
        // cube graph is connected, so some edge joins an inside corner to an
        // outside one.
        debug_assert!(crossings > 0);
        if crossings == 0 {
            return;
        }

        let inv = R::from_f64(f64::from(crossings)).recip();
        out.push_whole_cell([
            origin[0] + cell_size * (R::from_f64(f64::from(base[0])) + sum[0] * inv),
            origin[1] + cell_size * (R::from_f64(f64::from(base[1])) + sum[1] * inv),
            origin[2] + cell_size * (R::from_f64(f64::from(base[2])) + sum[2] * inv),
        ]);
    }
}

impl<R: Real> SurfaceNets<R> {
    /// A mesher that has allocated nothing yet.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            mesher: DualMesher::new(),
        }
    }

    /// How many Laplacian smoothing passes to run. Default zero.
    ///
    /// Each pass replaces every vertex with the average of itself and the
    /// vertices of the face-adjacent active cells. It visibly relaxes the
    /// staircase that a coarse grid produces.
    ///
    /// **It can move a vertex outside its own cell**, and vertices that leave
    /// their cells are how a dual method starts self-intersecting. The
    /// self-intersection counter is the way to see that happening; A-009 is the
    /// ticket that measures the clamp which prevents it.
    pub fn set_smoothing_passes(&mut self, passes: u32) {
        self.mesher.smoothing_passes = passes;
    }

    /// Extract the zero level set into `out`.
    ///
    /// `shape` counts **samples**, so `[n; 3]` spans `n - 1` cells per axis.
    /// `origin` is the world position of sample `[0, 0, 0]`.
    ///
    /// # Conventions
    ///
    /// Identical to Marching Cubes, deliberately, so the two can be compared on
    /// the same field without an adapter: negative is inside, zero is outside,
    /// winding is counter-clockwise seen from outside the solid, and normals come
    /// from the field's gradient.
    ///
    /// Output is triangles, two per quad, because [`MeshSink`] takes triangles.
    /// The quad structure is still visible in the output — each pair shares a
    /// diagonal.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if any axis has fewer
    /// than two samples. Surface Nets places at most one vertex per cell, and a
    /// grid with no cells has nothing to place.
    pub fn extract<S, M>(
        &mut self,
        sdf: &S,
        shape: &impl Shape3,
        origin: [R; 3],
        cell_size: R,
        out: &mut M,
    ) -> crate::Result<()>
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        self.mesher
            .extract(&Centroid, sdf, shape, origin, cell_size, out)
    }
}

impl<R: Real> Default for SurfaceNets<R> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests;
