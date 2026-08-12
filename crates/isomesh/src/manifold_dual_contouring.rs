//! Manifold Dual Contouring — one vertex per surface component, not per cell.
//!
//! [`DualContouring`](crate::dual_contouring::DualContouring) places one vertex
//! per cell. Where two sheets of surface pass through the same cell they are
//! forced to share it, and the mesh pinches: `non_manifold_edges` is 128 on the
//! shootout's fields (M-53), and M-29 established that this residue is *all* of
//! what the A-009 cell clamp could not remove. This module removes it.
//!
//! # The construction, and where it comes from
//!
//! Schaefer, Ju & Warren, *Manifold Dual Contouring*, TVCG 13(3) 2007,
//! `10.1109/TVCG.2007.1012`, §2.2 and §3 — read this session:
//!
//! > A problem of DC that has been of common interest in almost all subsequent
//! > work is the restriction of DC in maintaining no more than one contour vertex
//! > within each grid cell. To relax this restriction on a uniform grid, multiple
//! > contour components in a cell can be detected either by identifying
//! > edge-connected components of positive (or negative) cell corners \[9],\[12]
//! > or by utilizing the cycles in the MC lookup table \[7],\[13]. In this paper,
//! > we follow the Dual MC approach of Nielson \[13].
//!
//! > Nielson associates one vertex with each cycle of a modified MC table \[26].
//! > Since each cycle consists of a list of edges on the cubic cell, each vertex
//! > is associated with a set of edges, and each edge is associated with exactly
//! > one vertex. To create polygons, the algorithm constructs one polygon
//! > connecting the vertices associated with that edge in the four adjacent
//! > cells. **This algorithm creates a quadrilateral surface that is the dual of
//! > the surface created using MC** … **this surface is always a manifold because
//! > the original MC algorithm always constructs a manifold and the dual
//! > preserves the topology of the surface.**
//!
//! > To incorporate Hermite data into Nielson's Dual MC algorithm, we simply
//! > construct a QEF for each vertex using the Hermite data on the edges
//! > associated with that vertex. We place this vertex at the location that
//! > minimizes that error function.
//!
//! So the uniform-grid algorithm is Nielson's Dual MC (IEEE Vis 2004,
//! `10.1109/VISUAL.2004.28`) with the QEF from A-007 solved per cycle. **The
//! manifold guarantee is Nielson's, not the 2007 paper's** — the 2007
//! contribution is the octree vertex-clustering criterion, which needs an octree
//! this crate does not have and which G-004 does not want (LOD here mips the
//! *field*, not the mesh). Recorded as a scoping correction on A-010.
//!
//! # Why the cycles are derived, not transcribed
//!
//! Nielson's paper is paywalled, exactly as Doi & Koide was at A-003, and hard
//! rule 5 forbids inventing the table it would have supplied. Nothing is
//! invented: the cycles come from
//! [`crate::marching_cubes::table::segment_links`], the same
//! function the 256-case Marching Cubes table is *built* from. `next[e]` is the
//! cut edge the segment leaving `e` arrives at, so walking it from every unvisited
//! cut edge enumerates precisely the cycles the paper means, for precisely the
//! configuration this crate's Marching Cubes would triangulate.
//!
//! That identity is the whole guarantee. The output is the dual of *this crate's*
//! Marching Cubes — the one measured at zero non-manifold edges and zero
//! self-intersections in M-53 — rather than the dual of some other paper's. And
//! because the dual of a surface has `V' = F`, `E' = E`, `F' = V`, the Euler
//! characteristic is carried across unchanged, which
//! `euler_characteristic_matches_marching_cubes` asserts on all seven fields.
//!
//! # What it does not fix
//!
//! Self-intersection. ✗2 records ODC (2024) measuring Manifold Dual Contouring at
//! **100% of models self-intersecting**, and Manson & Schaefer's within-cell
//! partition argument — the one [`Clamp::ToCell`] rests on — assumes *one* vertex
//! per cell, which is the assumption this module drops. The counts are measured
//! rather than assumed; see A-010's archive entry.

use crate::cube::{EDGE_COUNT, is_inside};
use crate::dual::{CellVertices, DualMesher, VertexRule};
use crate::dual_contouring::{Clamp, apply_clamp, solve};
use crate::hermite::HermiteCell;
use crate::marching_cubes::FaceAmbiguity;
use crate::marching_cubes::ambiguity::joined_mask;
use crate::marching_cubes::table::{AMBIGUOUS_FACES, NO_EDGE, segment_links};
use crate::{MeshSink, Real, Sdf, Shape3};

/// Manifold Dual Contouring over a sampled grid.
///
/// Sharp features from the same QEF as
/// [`DualContouring`](crate::dual_contouring::DualContouring), and manifold
/// connectivity from splitting a cell's vertex once per Marching Cubes cycle.
///
/// Costs one extra vertex per multi-sheet cell and nothing anywhere else — the
/// cycle walk is over a 12-entry array the table already computes, and cells with
/// one cycle take the identical path they took before.
///
/// # Example
///
/// ```
/// use isomesh::{MeshBuffer, RuntimeShape3};
/// use isomesh::fields::capped_gyroid;
/// use isomesh::manifold_dual_contouring::ManifoldDualContouring;
/// use isomesh::validate::{ValidateConfig, validate};
///
/// let h = 4.0 / 32.0;
/// let mut mdc = ManifoldDualContouring::<f64>::new();
/// let mut out = MeshBuffer::<f64>::new();
/// let shape = RuntimeShape3::new([33; 3])?;
/// mdc.extract(&capped_gyroid::<f64>(), &shape, [-2.0; 3], h, &mut out)?;
///
/// // The gyroid is the field one-vertex-per-cell cannot mesh cleanly.
/// let report = validate(&out, &ValidateConfig::from_cell_size(h)?);
/// assert_eq!(report.non_manifold_edges, 0);
/// # Ok::<(), isomesh::Error>(())
/// ```
#[derive(Debug)]
pub struct ManifoldDualContouring<R: Real> {
    mesher: DualMesher<R>,
    rule: CycleQef,
}

/// The Manifold Dual Contouring vertex rule: one QEF solve per Marching Cubes
/// cycle.
///
/// Identical to [`Qef`](crate::dual_contouring::Qef) on any cell the surface
/// passes through once — the single cycle owns every cut edge, the restricted
/// cell is the whole cell, and the solve sees the same crossings. The two differ
/// only where a cell carries more than one sheet, which is the only place the
/// pinch was.
#[derive(Clone, Copy, Debug, Default)]
pub struct CycleQef {
    /// Whether to confine each solved vertex to its own cell.
    pub clamp: Clamp,
    /// How an ambiguous face is resolved.
    ///
    /// This must be read as choosing *which* Marching Cubes surface to take the
    /// dual of. Both settings are manifold; they disagree about the topology
    /// being dualised, exactly as they do in
    /// [`MarchingCubes`](crate::marching_cubes::MarchingCubes).
    pub face_ambiguity: FaceAmbiguity,
}

impl<R: Real> VertexRule<R> for CycleQef {
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

        let mut case = 0u8;
        for (c, value) in corner.iter().enumerate() {
            if is_inside(*value) {
                case |= 1 << c;
            }
        }
        let ambiguous = match self.face_ambiguity {
            FaceAmbiguity::Separate => 0,
            FaceAmbiguity::AsymptoticDecider => AMBIGUOUS_FACES[case as usize],
        };
        let next = segment_links(case, joined_mask(corner, ambiguous));

        // One Hermite sample of the cell, shared by every component: the
        // crossings do not depend on which cycle owns them, and re-deriving them
        // per component would repeat a gradient evaluation per crossing.
        let cell = HermiteCell::from_corners(sdf, corner, cell_origin, cell_size);

        let mut visited = 0u16;
        for start in 0..EDGE_COUNT as u8 {
            if next[start as usize] == NO_EDGE || visited & (1 << start) != 0 {
                continue;
            }

            // Walk this cut edge's cycle. Every cut edge is an entry on exactly
            // one of its two faces and an exit on the other, so the links form
            // disjoint cycles covering all of them.
            let mut edges = 0u16;
            let mut current = start;
            while visited & (1 << current) == 0 {
                visited |= 1 << current;
                edges |= 1 << current;
                current = next[current as usize];
            }

            // A cycle closes across faces, so it has at least three edges and the
            // restricted cell is never empty.
            debug_assert!(
                edges.count_ones() >= 3,
                "case {case:#010b} cycle {edges:#014b}"
            );
            let Some(x) = solve::solve(&cell.restricted(edges)) else {
                // Unreachable for a finite field: `solve` declines only on an
                // empty cell or a non-finite normal. Leaving the component
                // unplaced is caught by the engine, which then treats the whole
                // cell as inactive rather than emitting a quad corner that does
                // not exist.
                debug_assert!(false, "a non-empty component always solves");
                continue;
            };
            out.push_component(apply_clamp(self.clamp, x, cell_origin, cell_size), edges);
        }
    }
}

impl<R: Real> ManifoldDualContouring<R> {
    /// A mesher that has allocated nothing yet.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            mesher: DualMesher::new(),
            rule: CycleQef {
                clamp: Clamp::ToCell,
                face_ambiguity: FaceAmbiguity::Separate,
            },
        }
    }

    /// Whether solved vertices are confined to their own cells.
    ///
    /// Defaults to [`Clamp::ToCell`], as
    /// [`DualContouring`](crate::dual_contouring::DualContouring) does. Note the
    /// clamp's partition argument assumes one vertex per cell and this algorithm
    /// does not honour that assumption — see the module docs.
    pub fn set_clamp(&mut self, clamp: Clamp) {
        self.rule.clamp = clamp;
    }

    /// Which Marching Cubes surface to take the dual of on an ambiguous face.
    ///
    /// Defaults to [`FaceAmbiguity::Separate`] — Marching Cubes proper — for the
    /// same reason Marching Cubes does.
    pub fn set_face_ambiguity(&mut self, face_ambiguity: FaceAmbiguity) {
        self.rule.face_ambiguity = face_ambiguity;
    }

    /// Extract the zero level set into `out`.
    ///
    /// `shape` counts **samples**, so `[n; 3]` spans `n - 1` cells per axis.
    /// `origin` is the world position of sample `[0, 0, 0]`.
    ///
    /// # Conventions
    ///
    /// Identical to every other extractor here, deliberately, so they can be
    /// compared on the same field without an adapter: negative is inside, zero is
    /// outside, winding is counter-clockwise seen from outside the solid, and
    /// normals come from the field's gradient.
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

impl<R: Real> Default for ManifoldDualContouring<R> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests;
