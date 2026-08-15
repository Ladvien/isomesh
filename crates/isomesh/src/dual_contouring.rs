//! Dual Contouring — the sharp-feature method.
//!
//! Same topology as [`SurfaceNets`](crate::surface_nets::SurfaceNets), different vertex.
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
//! use isomesh::dual_contouring::DualContouring;
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
//!   how a dual method self-intersects. A-009 added the `(1−ε)` clamp and
//!   measured what it costs in sharpness against what it buys in
//!   self-intersections — [`Clamp::ToCell`] is the default (M-28, M-29), and
//!   [`Clamp::None`] remains only as the measured baseline configuration.
//! - **One vertex per cell.** Where two sheets of surface share a cell they share
//!   its vertex, and the mesh is non-manifold there. Inherited from the shared
//!   topology, measured on Surface Nets as M-4 and M-15, and fixed
//!   architecturally by A-010's vertex splitting rather than by patching.
//!
//! Neither is a bug to be fixed here. Both are the reason A-009 and A-010 exist.

pub mod solve;

use crate::dual::{CellVertices, DualMesher, VertexRule};
use crate::hermite::HermiteCell;
use crate::{MeshSink, Real, Sdf, Shape3};

/// Dual Contouring over a sampled grid.
///
/// Owns its scratch buffers; reuse it across chunks.
///
/// Deliberately exposes no smoothing control, unlike
/// [`SurfaceNets`](crate::surface_nets::SurfaceNets): a Laplacian pass averages each vertex
/// with its neighbours, which is exactly the operation that destroys the sharp
/// feature the solve just recovered.
/// # The rule is a type parameter, and that is X-002's whole point
///
/// `V` defaults to [`Qef`], so `DualContouring::<f32>::new()` means exactly what
/// it always meant and production monomorphises to one path with no runtime
/// branch anywhere. An *experiment* names a different rule —
/// `DualContouring::<f32, Centroid>::with_rule(Centroid)` — and gets a second,
/// separately compiled path that shares every other line of the algorithm.
///
/// This resolves a tension the crate had recorded and left alone.
/// `property/extraction.rs` argues, correctly, that adding a swappable table
/// *parameter* to the public API *"would be a second execution path in
/// production code, and the crate's rule is one"*. A type parameter is not that:
/// there is no value to branch on, so the compiler emits one rule per
/// instantiation and the binary a consumer ships contains only the rule they
/// named.
#[derive(Debug)]
pub struct DualContouring<R: Real, V: VertexRule<R> = Qef> {
    mesher: DualMesher<R>,
    rule: V,
}

/// Whether a solved vertex is confined to the cell that produced it.
///
/// Surface Nets' centroid is inside its cell by construction. A solved vertex is
/// not, and **a vertex outside its own cell is how a dual method
/// self-intersects** — the cells no longer partition space, so two patches can
/// overlap.
///
/// Manson & Schaefer's guarantee is exactly this constraint:
///
/// > As long as the dual vertices lie within their corresponding cells, this
/// > process cannot produce any inverted tetrahedra and creates a partition of
/// > space… Since our decomposition is one-to-one, our surface cannot intersect
/// > itself.
///
/// A-009 measured what it costs here; the numbers are on that ticket's archive
/// entry and in M-28.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Clamp {
    /// Leave the vertex where the solve put it.
    ///
    /// The baseline A-007 shipped, and what the clamp is measured against.
    None,
    /// Confine it to the cell, shrunk about its centre by `(1 − ε)`.
    #[default]
    ToCell,
}

/// How far inside the cell the clamp keeps a vertex, as a fraction of the cell.
///
/// The literature's `(1 − ε)` with `ε = 1e-4`: the cell is scaled about its own
/// centre, so the vertex stays *strictly* interior rather than landing on a face
/// it shares with a neighbour. Strictness is what rules out the degenerate
/// coincidences the partition argument would otherwise have to handle.
pub const CLAMP_EPSILON: f64 = 1e-4;

/// The Dual Contouring vertex rule: the regularized plane-intersection solve.
///
/// Builds the cell's [`HermiteCell`] — a position and a surface normal per edge
/// crossing — and hands it to [`solve::solve`]. All of the mathematics, and all
/// of the justification for this particular form of it, is in that module.
#[derive(Clone, Copy, Debug, Default)]
pub struct Qef {
    /// Whether to confine the result to its cell.
    pub clamp: Clamp,
    /// The Tikhonov regularizer, or [`None`] for
    /// [`solve::LAMBDA`].
    ///
    /// The sharpness knob. Toward zero the solve reproduces a three-plane corner
    /// exactly and lets a flat cell's vertex fly off; large, it pulls every
    /// vertex to the centroid of its crossings and rounds every edge over. See
    /// [`solve_with`](crate::dual_contouring::solve::solve_with).
    ///
    /// `Option` rather than a bare value so that [`Default`] means "the default
    /// λ" without this type having to restate it — one definition of the
    /// number, in the module that explains it.
    pub lambda: Option<f64>,
}

impl<R: Real> VertexRule<R> for Qef {
    /// One vertex owning the whole cell.
    ///
    /// Where two sheets share a cell they share this vertex, and the mesh pinches
    /// — the structural defect A-010 lifts by solving the QEF per surface
    /// component instead. See
    /// [`ManifoldDualContouring`](crate::manifold_dual_contouring::ManifoldDualContouring).
    fn place<S: Sdf<Scalar = R>>(
        &self,
        sdf: &S,
        corner: &[R; 8],
        base: [u32; 3],
        origin: [R; 3],
        cell_size: R,
        out: &mut CellVertices<R>,
    ) {
        let cell_origin = [
            origin[0] + cell_size * R::from_f64(f64::from(base[0])),
            origin[1] + cell_size * R::from_f64(f64::from(base[1])),
            origin[2] + cell_size * R::from_f64(f64::from(base[2])),
        ];
        // `from_corners` costs one gradient evaluation per crossing and re-uses
        // the corner samples the engine already took, so this adds field
        // evaluations only where the surface actually is.
        let cell = HermiteCell::from_corners(sdf, corner, cell_origin, cell_size);
        let lambda = R::from_f64(self.lambda.unwrap_or(solve::LAMBDA));
        let Some(x) = solve::solve_with(&cell, lambda) else {
            return;
        };
        out.push_whole_cell(apply_clamp(self.clamp, x, cell_origin, cell_size));
    }
}

/// Confine a solved vertex to its own cell, if the rule says to.
///
/// Shared by [`Qef`] and
/// [`CycleQef`](crate::manifold_dual_contouring::CycleQef) so the two cannot
/// disagree about what "inside the cell" means — a disagreement that would show
/// up only as a self-intersection count nobody could account for.
///
/// The cell is scaled about its own centre by `(1 − ε)`. `min`/`max` is a
/// *continuous* operation, which is the reason this does not reintroduce the
/// failure ✗12 is about: a vertex pushed against a wall slides along it as the
/// field moves, where an SVD rank branch jumps. Clamping costs sharpness, not
/// stability.
pub(crate) fn apply_clamp<R: Real>(
    clamp: Clamp,
    x: [R; 3],
    cell_origin: [R; 3],
    cell_size: R,
) -> [R; 3] {
    match clamp {
        Clamp::None => x,
        Clamp::ToCell => {
            let half = cell_size * R::HALF;
            let inset = half * R::from_f64(1.0 - CLAMP_EPSILON);
            let mut out = x;
            for (axis, slot) in out.iter_mut().enumerate() {
                let centre = cell_origin[axis] + half;
                *slot = slot.clamp(centre - inset, centre + inset);
            }
            out
        }
    }
}

impl<R: Real> DualContouring<R> {
    /// A mesher that has allocated nothing yet.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            mesher: DualMesher::new(),
            rule: Qef {
                clamp: Clamp::ToCell,
                lambda: None,
            },
        }
    }

    /// The Tikhonov regularizer, the knob that trades sharpness for stability.
    ///
    /// `None` restores [`solve::LAMBDA`]. Toward zero, a three-plane corner comes
    /// out exactly and a flat cell's vertex flies off; large, every edge rounds
    /// over. `sharp_features` is the example that makes that visible.
    pub fn set_lambda(&mut self, lambda: Option<f64>) {
        self.rule.lambda = lambda;
    }

    /// Whether solved vertices are confined to their own cells.
    ///
    /// Defaults to [`Clamp::ToCell`]. See [`Clamp`] for why, and A-009's archive
    /// entry for what it measured.
    pub fn set_clamp(&mut self, clamp: Clamp) {
        self.rule.clamp = clamp;
    }
}

impl<R: Real, V: VertexRule<R>> DualContouring<R, V> {
    /// A mesher that places its vertices with `rule`.
    ///
    /// The experiment's entry point, and the reason the rule is a type
    /// parameter. [`new`](DualContouring::new) is this with [`Qef`], and the two
    /// share every line of the algorithm below the placement — which is what
    /// makes a comparison between two rules a measurement of the *rule* rather
    /// than of two implementations that happen to differ.
    ///
    /// ```
    /// use isomesh::dual_contouring::DualContouring;
    /// use isomesh::surface_nets::Centroid;
    ///
    /// // Dual Contouring's topology with Surface Nets' placement.
    /// let mesher = DualContouring::<f32, _>::with_rule(Centroid);
    /// # let _ = mesher;
    /// ```
    #[must_use]
    pub const fn with_rule(rule: V) -> Self {
        Self {
            mesher: DualMesher::new(),
            rule,
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
            .extract(&self.rule, sdf, shape, origin, cell_size, out)
    }
}

impl<R: Real> Default for DualContouring<R> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests;
