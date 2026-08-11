//! Dual Contouring — the sharp-feature method.
//!
//! Same topology as [`SurfaceNets`](crate::sn::SurfaceNets), different vertex.
//! That sentence is the whole algorithm, and it is not a simplification: V-19
//! records that the original paper's method *is* one vertex per sign-changing
//! cube and one quad per sign-changing edge, with only placement differing. Both
//! methods therefore run the same engine in `crate::dual` and supply a
//! vertex rule each.
//!
//! Surface Nets averages the cell's edge crossings, which cannot land on a
//! corner. Dual Contouring solves for the point where the crossings' **tangent
//! planes** agree, which can — that is what A-006's Hermite data is for, and why
//! this method is the one a CAD tool needs.
//!
//! ```
//! use isomesh::{MeshBuffer, RuntimeShape3};
//! use isomesh::dc::DualContouring;
//! use isomesh::fields::BoxExact;
//!
//! let mut dc = DualContouring::<f32>::new();
//! let mut out = MeshBuffer::<f32>::new();
//! let shape = RuntimeShape3::new([33; 3])?;
//! dc.extract(&BoxExact::<f32>::canonical(), &shape, [-2.0; 3], 0.125, &mut out)?;
//!
//! assert!(out.triangle_count() > 0);
//! # Ok::<(), isomesh::Error>(())
//! ```
//!
//! # What it costs
//!
//! Two things, both structural rather than incidental.
//!
//! - **The vertex can leave its cell.** Surface Nets' centroid is inside by
//!   construction; a solved vertex is not, and a vertex outside its own cell is
//!   how a dual method self-intersects. A-009 adds the `(1−ε)` clamp and measures
//!   what it costs in sharpness against what it buys in self-intersections.
//!   Until then this is the unclamped solve, so that A-009 has a baseline.
//! - **One vertex per cell.** Where two sheets of surface share a cell they share
//!   its vertex, and the mesh is non-manifold there. Inherited from the shared
//!   topology, measured on Surface Nets as M-4 and M-15, and fixed
//!   architecturally by A-010's vertex splitting rather than by patching.
//!
//! Neither is a bug to be fixed here. Both are the reason A-009 and A-010 exist.

pub mod solve;

use crate::dual::{DualMesher, VertexRule};
use crate::hermite::HermiteCell;
use crate::{MeshSink, Real, Sdf, Shape3};

/// Dual Contouring over a sampled grid.
///
/// Owns its scratch buffers; reuse it across chunks.
///
/// Deliberately exposes no smoothing control, unlike
/// [`SurfaceNets`](crate::sn::SurfaceNets): a Laplacian pass averages each vertex
/// with its neighbours, which is exactly the operation that destroys the sharp
/// feature the solve just recovered.
#[derive(Debug)]
pub struct DualContouring<R: Real> {
    mesher: DualMesher<R>,
}

/// The Dual Contouring vertex rule: the regularized plane-intersection solve.
///
/// Builds the cell's [`HermiteCell`] — a position and a surface normal per edge
/// crossing — and hands it to [`solve::solve`]. All of the mathematics, and all
/// of the justification for this particular form of it, is in that module.
#[derive(Clone, Copy, Debug, Default)]
pub struct Qef;

impl<R: Real> VertexRule<R> for Qef {
    fn place<S: Sdf<Scalar = R>>(
        &self,
        sdf: &S,
        corner: &[R; 8],
        base: [u32; 3],
        origin: [R; 3],
        cell_size: R,
    ) -> Option<[R; 3]> {
        let cell_origin = [
            origin[0] + cell_size * R::from_f64(f64::from(base[0])),
            origin[1] + cell_size * R::from_f64(f64::from(base[1])),
            origin[2] + cell_size * R::from_f64(f64::from(base[2])),
        ];
        // `from_corners` costs one gradient evaluation per crossing and re-uses
        // the corner samples the engine already took, so this adds field
        // evaluations only where the surface actually is.
        let cell = HermiteCell::from_corners(sdf, corner, cell_origin, cell_size);
        solve::solve(&cell)
    }
}

impl<R: Real> DualContouring<R> {
    /// A mesher that has allocated nothing yet.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            mesher: DualMesher::new(),
        }
    }

    /// Extract the zero level set into `out`.
    ///
    /// `shape` counts **samples**, so `[n; 3]` spans `n - 1` cells per axis.
    /// `origin` is the world position of sample `[0, 0, 0]`.
    ///
    /// # Conventions
    ///
    /// Identical to Marching Cubes and Surface Nets, deliberately, so the three
    /// can be compared on the same field without an adapter: negative is inside,
    /// zero is outside, winding is counter-clockwise seen from outside the solid,
    /// and normals come from the field's gradient.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if any axis has fewer
    /// than two samples.
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
            .extract(&Qef, sdf, shape, origin, cell_size, out)
    }
}

impl<R: Real> Default for DualContouring<R> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests;
