//! Handing a mesh to a physics engine, and knowing first whether it will take it.
//!
//! A collider is not a render mesh with a different material. A renderer will
//! draw a zero-area triangle, an unwelded seam and an open surface without
//! complaint; a physics engine will variously refuse them, build a broken
//! adjacency structure from them, or answer "is this point inside" wrongly. This
//! module says which of those a given [`MeshBuffer`] would hit, **before** it is
//! handed over.
//!
//! # What this deliberately does not do
//!
//! Depend on a physics engine. `parry3d` is a **dev-dependency** here, not an
//! optional feature, and the ticket that asked for the feature is amended rather
//! than followed. Three reasons:
//!
//! - **`parry3d` takes `Vec<[u32; 3]>` for its index buffer** — the architecture
//!   research calls it *"the most math-opinionated crate in the ecosystem"* and
//!   notes it still uses plain arrays there. So the only thing a conversion adds
//!   is mapping `[f32; 3]` to parry's vector type, which is one line at the call
//!   site and is shown below.
//! - **The core crate stays at one dependency in every configuration.** An
//!   optional dependency is still a version this manifest has to track, and
//!   `parry3d` is pre-1.0 and moving fast — it migrated off nalgebra onto
//!   glam/glamx in 0.26.0, which is why it wants glam `^0.33` where the engine
//!   version this workspace pins against wants `^0.32`. That skew is documented
//!   in `CLAUDE.md`'s hard version pins and is the reason every public signature
//!   here is an array.
//! - **The acceptance criterion is still met, and by parry itself.**
//!   `a_carved_shape_builds_a_parry_trimesh` builds a real `TriMesh` from a
//!   carved field and runs parry's own topology check on it. A dev-dependency
//!   does not propagate, so `cargo tree -p isomesh -e normal` is still two
//!   packages.
//!
//! The conversion, in full, for a consumer who wants it:
//!
//! ```ignore
//! let indices = isomesh::collider::triangle_indices(&mesh);
//! let vertices = mesh.positions.iter().map(|p| Vector::new(p[0], p[1], p[2])).collect();
//! let trimesh = TriMesh::new(vertices, indices)?;
//! ```
//!
//! # Weld first. This is the part that bites.
//!
//! Meshing a world in chunks produces, at every seam, two vertices at the same
//! position that no index shares. A renderer cannot tell. A physics engine sees a
//! **boundary edge** — a hole — and a character walks through it. A-013 measured
//! the seam: 80 boundary edges and 40 duplicated vertices between two chunks of a
//! torus, closing to zero after a weld (M-46). So
//! [`readiness`] reports `duplicate_vertices` first, and
//! `an_unwelded_seam_is_reported_as_not_ready` is the test that keeps it honest.

use alloc::vec::Vec;

use crate::validate::{MeshReport, ValidateConfig, validate_indexed};
use crate::{MeshBuffer, Real};

/// The index buffer in the shape a physics engine wants it.
///
/// [`MeshBuffer::indices`](MeshBuffer) is flat, because that is what a GPU takes.
/// `parry3d` — and every other engine checked — takes triples. A trailing partial
/// triangle is dropped, exactly as the validator counts it in `trailing_indices`
/// rather than rounding it up into a triangle nobody asked for.
#[must_use]
pub fn triangle_indices<R: Real>(mesh: &MeshBuffer<R>) -> Vec<[u32; 3]> {
    mesh.indices
        .chunks_exact(3)
        .map(|t| [t[0], t[1], t[2]])
        .collect()
}

/// What a physics engine would make of this mesh.
///
/// Every count is read from the T-001 validator rather than recomputed, so this
/// cannot drift from what the rest of the suite reports about the same mesh. What
/// is added is the *interpretation*: which of those numbers a collider actually
/// cares about, and which it does not.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ColliderReadiness {
    /// Triangles the engine would receive — the validator's `faces`, which
    /// excludes any it had to skip for a bad index.
    pub triangles: u64,
    /// Triangles the validator skipped: an out-of-range or a repeated index.
    pub triangles_skipped: u64,
    /// Vertices at the same position that no index shares.
    ///
    /// **The chunk-seam failure.** Invisible to a renderer, a hole to a physics
    /// engine. Weld before exporting.
    pub duplicate_vertices: u64,
    /// Triangles with no area, so no normal.
    ///
    /// Marching Cubes emits these whenever a grid corner sits near zero — it is
    /// the algorithm, not a bug (T-001 records them rather than gating on them).
    /// A contact normal derived from one is undefined, so a collider is the one
    /// consumer that has to care.
    pub degenerate_triangles: u64,
    /// Indices naming a vertex that does not exist.
    pub out_of_range_indices: u64,
    /// Indices past the last complete triangle.
    pub trailing_indices: u64,
    /// Positions that are not finite.
    pub non_finite_positions: u64,
    /// Edges used by exactly one triangle — the surface is open here.
    pub boundary_edges: u64,
    /// Edges used by three or more triangles.
    pub non_manifold_edges: u64,
    /// Vertices whose incident faces form more than one fan.
    ///
    /// **The one manifold fault no edge counter can see.** A bowtie — two cones
    /// meeting at a shared apex — gives every edge exactly two faces, so
    /// [`non_manifold_edges`](Self::non_manifold_edges) is zero, `χ` can come out
    /// right, and the surface still is not a 2-manifold. The apex has no single
    /// normal, so a pseudo-normal query there is arbitrary rather than wrong.
    ///
    /// Read from the validator's link walk, not recomputed. Forwarded at T-021
    /// after M-300 found it missing while auditing what a convex decomposer
    /// requires of its input — CoACD's plane cutting is stated over manifold
    /// meshes, and this is the half of "manifold" the collider view was
    /// silently dropping.
    pub non_manifold_vertices: u64,
    /// Edges whose two triangles traverse them the same way.
    pub inconsistently_oriented_edges: u64,
}

impl ColliderReadiness {
    /// Whether this mesh can be handed over at all.
    ///
    /// The four structural faults: a non-finite position, an index naming a
    /// vertex that does not exist, a trailing partial triangle, and no triangles
    /// whatsoever. `parry3d`'s own constructor rejects the last of those —
    /// *"the index buffer is empty (at least one triangle is required)"* — and
    /// will happily accept the rest and misbehave later, which is why they are
    /// checked here instead.
    ///
    /// **A duplicate vertex is not a structural fault and does not fail this.**
    /// It is a *correctness* fault for a chunked world and a non-issue for a
    /// single mesh, and only the caller knows which it has. See
    /// [`is_seam_free`](Self::is_seam_free).
    #[must_use]
    pub fn is_usable(&self) -> bool {
        self.triangles > 0
            && self.triangles_skipped == 0
            && self.out_of_range_indices == 0
            && self.trailing_indices == 0
            && self.non_finite_positions == 0
    }

    /// Whether the mesh has no unwelded seam.
    ///
    /// Check this on anything assembled from more than one chunk.
    #[must_use]
    pub fn is_seam_free(&self) -> bool {
        self.duplicate_vertices == 0
    }

    /// Whether an inside/outside test on this mesh would mean anything.
    ///
    /// `parry3d` computes pseudo-normals for that under
    /// `TriMeshFlags::ORIENTED`, and a pseudo-normal is only defined on a closed,
    /// consistently oriented surface. An open surface has no inside, and a
    /// non-manifold edge has no single normal to speak of — so a query against
    /// one returns an answer that is arbitrary rather than wrong, which is worse.
    ///
    /// **A bowtie vertex fails this too, since T-021.** It has to be checked
    /// separately because no edge counter can see it: two cones sharing an apex
    /// give every edge exactly two faces (M-300). Before T-021 this returned
    /// `true` for such a mesh, which is the arbitrary-answer case above, at the
    /// one point on the surface where it is guaranteed.
    #[must_use]
    pub fn supports_inside_outside(&self) -> bool {
        self.is_usable()
            && self.boundary_edges == 0
            && self.non_manifold_edges == 0
            && self.non_manifold_vertices == 0
            && self.inconsistently_oriented_edges == 0
    }
}

/// Read a mesh through a collider's eyes.
///
/// Runs the T-001 validator and reinterprets it. `cfg` decides the weld epsilon
/// and the degenerate-area threshold, both relative to the grid spacing — see
/// [`ValidateConfig::from_cell_size`].
#[must_use]
pub fn readiness<R: Real>(mesh: &MeshBuffer<R>, cfg: &ValidateConfig) -> ColliderReadiness {
    from_report(&validate_indexed(&mesh.positions, &mesh.indices, cfg))
}

/// The same reinterpretation, from a report the caller already has.
///
/// Exists so a caller that has just validated a mesh does not validate it twice;
/// the validator is a multi-pass walk over every triangle and is not free.
#[must_use]
pub fn from_report(report: &MeshReport) -> ColliderReadiness {
    ColliderReadiness {
        triangles: report.faces,
        triangles_skipped: report.faces_skipped,
        duplicate_vertices: report.duplicate_vertices,
        degenerate_triangles: report.degenerate_triangles,
        out_of_range_indices: report.out_of_range_indices,
        trailing_indices: report.trailing_indices,
        non_finite_positions: report.non_finite_positions,
        boundary_edges: report.boundary_edges,
        non_manifold_edges: report.non_manifold_edges,
        non_manifold_vertices: report.non_manifold_vertices,
        inconsistently_oriented_edges: report.inconsistently_oriented_edges,
    }
}

#[cfg(test)]
mod tests;
