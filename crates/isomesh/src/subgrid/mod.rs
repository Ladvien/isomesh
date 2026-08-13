//! Subgrid Marching Tetrahedra — features finer than one cell.
//!
//! Every algorithm in this crate so far shares one limitation, and Baktash,
//! Gillespie & Crane state it exactly (`10.48550/arXiv.2606.00454`, SIGGRAPH
//! 2026, §1, read this session):
//!
//! > Output frequency is limited by grid resolution, since the sign test used by
//! > these methods produces at most one vertex per grid edge. Hence, when `f` has
//! > multiple zeros along an edge, one gets **aliasing of fine geometric
//! > features, thin sheets, etc.** One common remedy is to simply increase grid
//! > resolution — yielding over-tessellated output which must be aggressively
//! > decimated.
//!
//! A-005 measured that limitation from the other side: `thin_plate` is 0.4 cells
//! thick and greedy quads returns **zero triangles** for it, because no cell
//! centre is inside. A feature thinner than a cell does not exist to a
//! sign-based method.
//!
//! The fix is to stop asking *what sign is this vertex* and start asking **how
//! many times does the surface cross this edge**:
//!
//! > We replace 0-dimensional sampling (evaluate `f` at each grid node), with
//! > 1-dimensional root finding (find all zeros of `f` along each grid edge)…
//! > Rather than rely on a finite lookup table of output configurations, we
//! > develop a deterministic algorithm that reconstructs a local polygonal
//! > approximation given **any number of intersections** along the edges of a
//! > tetrahedron.
//!
//! > This encoding sidesteps the usual Nyquist–Shannon limit, **putting no lower
//! > bound on the size of features that can be resolved on a fixed grid.**
//!
//! # Status
//!
//! A-014a — [`coordinates`], the encoding and its algebra — is here. The
//! reconstruction of boundary curves (§3.1), the surface fill with Steiner points
//! (§3.2–3.3) and all-roots edge finding (§4.3.2) are A-014b onward. Nothing here
//! places a vertex in world space yet.
//!
//! # Why the encoding is worth having on its own
//!
//! It is the part everything else is defined against, and it is the part that can
//! be checked against something this repo already trusts. The paper notes that
//!
//! > The marching tetrahedra algorithm **reinvented a small piece of this
//! > story**, but the isosurfacing literature makes no reference to the broader
//! > theory.
//!
//! So classic Marching Tetrahedra must be the special case of this encoding where
//! every edge coordinate is 0 or 1 — and
//! `classic_marching_tetrahedra_is_the_zero_one_case_of_this_encoding` asserts
//! exactly that against A-003's own table, on all sixteen tet configurations.

pub mod coordinates;
pub mod curves;
pub mod extract;
pub mod roots;
pub mod surface;
