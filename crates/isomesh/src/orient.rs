//! Give a mesh a **coherent** orientation, one connected component at a time.
//!
//! A locally correct winding is not a globally consistent one, and the
//! difference is invisible to every other check this crate runs. A mesh can pass
//! the Euler characteristic, edge-manifoldness *and* vertex-manifoldness with a
//! patch wound inside out; only the orientation test sees it, and by then the
//! mesh is already wrong.
//!
//! # Why this is a pass and not a rule inside an extractor
//!
//! Because orientation is a property of a *surface*, and an extractor working
//! one cell at a time does not have one. Subgrid Marching Tetrahedra is the case
//! that forced the point (A-014e): §3.2 fixes each polygon's vertex order from
//! its own boundary curve, which carries no relation to which side the field
//! calls inside, so each triangle is flipped to agree with the field's gradient
//! at its own centroid. That vote is deliberately **per triangle** rather than
//! per patch, because a sheet thinner than a cell puts two oppositely-facing
//! surfaces inside one tetrahedron and a per-patch vote would flip one of them
//! wrongly.
//!
//! The vote is right at every triangle and still not coherent across the mesh.
//! `gyroid` is where the two come apart: a triply periodic surface passes close
//! to itself, so the nearest sheet to a triangle's centroid is not always the
//! sheet that triangle is on. Measured at A-014f: of the 186 triangles standing
//! on its 138 inconsistently-oriented edges, **171 have a decisive vote on both
//! sides of the edge and the two answers disagree** (M-164).
//!
//! Propagation needs connectivity, and a per-tetrahedron soup has none until it
//! is welded — which is exactly why this is a pass the caller composes after
//! [`weld`](crate::weld) rather than something an extractor could do on its own.
//!
//! # The seed, and why it is the most confident triangle
//!
//! Propagation fixes *relative* orientation: it can make a component consistent,
//! but not decide which way out is. That still comes from the field, through the
//! per-vertex normals the extractor already wrote — so each component is seeded
//! from the triangle whose own geometric normal agrees most strongly with the
//! field's, and the rest of the component is made to agree with the seed.
//!
//! Taking the *most confident* triangle rather than the first one is the whole
//! point. A triangle whose gradient lies in its own plane has a vote of
//! approximately zero and no information in it; seeding from one would flip a
//! component on a coin toss. Sorting by confidence puts those triangles at the
//! end of the propagation, where they inherit an answer instead of inventing one.
//!
//! # What it does not do
//!
//! It does not make a non-orientable surface orientable — nothing can. On a
//! Möbius-like patch, propagation returns to its start with the opposite sign and
//! the pass reports the component as non-orientable rather than silently
//! choosing one of the two answers.

use alloc::collections::BTreeMap;
use alloc::vec;
use alloc::vec::Vec;

use crate::mesh::MeshBuffer;
use crate::real::Real;
use crate::vec3;

/// What [`orient`] did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OrientReport {
    /// Connected components found, counting by face adjacency across shared
    /// edges. A mesh of isolated triangles has one component per triangle.
    pub components: u64,
    /// Triangles whose winding was reversed.
    pub triangles_flipped: u64,
    /// Components on which propagation returned to a triangle it had already
    /// oriented and disagreed with itself.
    ///
    /// Non-zero means that component is **not orientable** — a Möbius-like
    /// patch — and its triangles are left as propagation last set them. It is
    /// not an error: a mesh can legitimately contain one, and reporting it is
    /// more useful than refusing the whole call.
    pub non_orientable_components: u64,
    /// Triangles whose own normal and the field's disagreed by less than a
    /// right angle's worth of confidence — recorded because a component made
    /// entirely of these was oriented on very little evidence.
    pub low_confidence_seeds: u64,
}

impl OrientReport {
    /// Whether every component came out orientable.
    #[must_use]
    pub const fn is_orientable(&self) -> bool {
        self.non_orientable_components == 0
    }

    /// Whether nothing needed changing.
    #[must_use]
    pub const fn is_noop(&self) -> bool {
        self.triangles_flipped == 0
    }
}

/// An undirected edge, as its two vertex indices with the smaller first.
///
/// Undirected on purpose: the *directed* edge is what says which way a triangle
/// is wound, and two correctly-wound neighbours traverse their shared edge in
/// opposite directions. Keying on the undirected form is what lets those two be
/// found as neighbours in the first place.
fn undirected(a: u32, b: u32) -> (u32, u32) {
    if a <= b { (a, b) } else { (b, a) }
}

/// Make every connected component of `mesh` coherently oriented.
///
/// Run this **after** [`weld`](crate::weld::Welder::weld). On an unwelded
/// triangle soup every edge belongs to exactly one triangle, so there are no
/// neighbours, every triangle is its own component, and this does nothing at
/// considerable expense.
///
/// Winding is left as counter-clockwise seen from outside the solid, the
/// convention every extractor here uses — the seed triangle of each component is
/// oriented against the field's own normals and the rest follows it.
///
/// # Errors
///
/// [`Error::IndexOutOfRange`](crate::Error::IndexOutOfRange) if any index does
/// not address a vertex. Checked at the door rather than during the walk,
/// because a partially reoriented mesh is worse than a rejected one.
///
/// # Example
///
/// ```
/// use isomesh::MeshBuffer;
/// use isomesh::orient::orient;
///
/// // Two triangles sharing edge 1-2, the second wound the wrong way round.
/// let mut mesh = MeshBuffer::<f64>::new();
/// mesh.positions = vec![
///     [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0],
/// ];
/// mesh.normals = vec![[0.0, 0.0, 1.0]; 4];
/// // Coherent would be `1, 3, 2`: two neighbours must traverse their shared
/// // edge in opposite directions. This one goes the same way round as its
/// // neighbour, so it is inside out relative to it.
/// mesh.indices = vec![0, 1, 2, /* reversed: */ 1, 2, 3];
///
/// let report = orient(&mut mesh)?;
///
/// assert_eq!(report.components, 1);
/// assert_eq!(report.triangles_flipped, 1);
/// assert!(report.is_orientable());
/// # Ok::<(), isomesh::Error>(())
/// ```
pub fn orient<R: Real>(mesh: &mut MeshBuffer<R>) -> crate::Result<OrientReport> {
    let vertices = mesh.positions.len();
    let faces = mesh.indices.len() / 3;

    for (at, &index) in mesh.indices.iter().enumerate() {
        if index as usize >= vertices {
            return Err(crate::Error::IndexOutOfRange {
                at: at as u64,
                index,
                vertices: vertices as u64,
            });
        }
    }
    if faces == 0 {
        return Ok(OrientReport::default());
    }

    // Edge -> the faces on it. A manifold edge has two; a boundary edge one; a
    // non-manifold edge more, and those are traversed too rather than skipped,
    // since refusing to orient a mesh because it has a bad edge somewhere else
    // would make this pass useless on exactly the meshes that need it.
    let mut incident: BTreeMap<(u32, u32), Vec<u32>> = BTreeMap::new();
    for face in 0..faces {
        let tri = &mesh.indices[face * 3..face * 3 + 3];
        for k in 0..3 {
            incident
                .entry(undirected(tri[k], tri[(k + 1) % 3]))
                .or_default()
                .push(face as u32);
        }
    }

    // Confidence per face, computed once and up front: the sort needs it, the
    // seeding needs it, and a closure holding `mesh` borrowed would stop the
    // walk from reorienting anything.
    let confidence: Vec<R> = (0..faces).map(|face| confidence_of(mesh, face)).collect();

    let mut visited = vec![false; faces];
    let mut report = OrientReport::default();
    let mut queue: Vec<u32> = Vec::new();

    // Seeds in descending confidence, so a component is always entered through
    // its best-evidenced triangle. Stable by face index on a tie, so the result
    // is a function of the mesh and not of the map's iteration order.
    let mut order: Vec<usize> = (0..faces).collect();
    order.sort_by(|&a, &b| {
        confidence[b]
            .total_cmp(&confidence[a])
            .then_with(|| a.cmp(&b))
    });

    for seed in order {
        if visited[seed] {
            continue;
        }
        report.components += 1;

        // Orient the seed against the field, then let the component follow it.
        if confidence[seed] <= R::ZERO {
            report.low_confidence_seeds += 1;
        }
        if faces_outward(mesh, seed) == Some(false) {
            flip(mesh, seed);
            report.triangles_flipped += 1;
        }

        visited[seed] = true;
        queue.clear();
        queue.push(seed as u32);
        let mut disagreed = false;

        while let Some(face) = queue.pop() {
            let tri = face_indices(mesh, face as usize);
            for k in 0..3 {
                let (from, to) = (tri[k], tri[(k + 1) % 3]);
                let Some(neighbours) = incident.get(&undirected(from, to)) else {
                    continue;
                };
                // **Propagation stops at a non-manifold edge (A-019).** With
                // three or more faces on one edge there is no "the neighbour" to
                // agree with: the walk would take all of them, commit to
                // whichever it reached first, and carry a *consistent* winding
                // across a whole patch that is consistent with the wrong side.
                // That is how orientation used to make the count **worse** —
                // 1,580 → 2,422 on `noise_cavity` at 25³ (M-213). Stopping
                // instead leaves each sheet to be seeded and oriented on its own
                // evidence, which is the only evidence there is.
                if neighbours.len() > 2 {
                    continue;
                }
                for &other in neighbours {
                    if other == face {
                        continue;
                    }
                    // Two coherently-wound faces traverse a shared edge in
                    // *opposite* directions. Same direction means one of them is
                    // inside out relative to the other.
                    let agrees = !traverses(mesh, other, from, to);
                    if visited[other as usize] {
                        // Already fixed. If it disagrees now, propagation has
                        // come back around with the opposite sign, which is what
                        // a non-orientable component looks like from the inside.
                        disagreed |= !agrees;
                        continue;
                    }
                    visited[other as usize] = true;
                    if !agrees {
                        flip(mesh, other as usize);
                        report.triangles_flipped += 1;
                    }
                    queue.push(other);
                }
            }
        }

        if disagreed {
            report.non_orientable_components += 1;
        }
    }

    Ok(report)
}

/// Whether face `face` traverses the directed edge `from -> to`.
fn traverses<R: Real>(mesh: &MeshBuffer<R>, face: u32, from: u32, to: u32) -> bool {
    let tri = face_indices(mesh, face as usize);
    (0..3).any(|k| tri[k] == from && tri[(k + 1) % 3] == to)
}

/// Whether this triangle's winding already points away from the solid, or
/// `None` when there is no direction to compare it with.
fn faces_outward<R: Real>(mesh: &MeshBuffer<R>, face: usize) -> Option<bool> {
    let tri = face_indices(mesh, face);
    let geometric = geometric_normal(mesh, tri);
    let field = field_direction(mesh, tri);
    let vote = vec3::dot(geometric, field);
    if vote > R::ZERO {
        Some(true)
    } else if vote < R::ZERO {
        Some(false)
    } else {
        None
    }
}

/// Reverse a triangle's winding, by swapping two of its three indices.
fn flip<R: Real>(mesh: &mut MeshBuffer<R>, face: usize) {
    mesh.indices.swap(face * 3 + 1, face * 3 + 2);
}

/// The field direction a triangle should agree with: the sum of its vertices'
/// own normals.
///
/// Summed rather than averaged because only the *sign* of a dot product with it
/// is ever used, and the divide would cost an operation to change nothing.
///
/// Returns the zero vector when the mesh carries no normal for a vertex, which
/// makes every vote on that triangle indecisive rather than wrong — it then
/// inherits its orientation from a neighbour, which is the correct outcome for a
/// triangle with no evidence of its own.
fn field_direction<R: Real>(mesh: &MeshBuffer<R>, tri: [u32; 3]) -> [R; 3] {
    let mut sum = [R::ZERO; 3];
    for index in tri {
        let Some(n) = mesh.normals.get(index as usize) else {
            return [R::ZERO; 3];
        };
        for axis in 0..3 {
            sum[axis] += n[axis];
        }
    }
    sum
}

/// How strongly a triangle's own geometric normal agrees with the field normals
/// its vertices carry, as `|cos|` between them.
///
/// The absolute value, because the question this answers is *how much
/// information the vote has*, not which way it points — a triangle whose
/// gradient lies in its own plane scores zero either way round, and that is
/// exactly the triangle that must not be trusted to seed a component.
///
/// Zero for anything indecisive: a degenerate triangle, a vertex with no normal,
/// or a non-finite one.
fn confidence_of<R: Real>(mesh: &MeshBuffer<R>, face: usize) -> R {
    let tri = face_indices(mesh, face);
    let geometric = geometric_normal(mesh, tri);
    let length = vec3::length(geometric);
    let field = field_direction(mesh, tri);
    let field_length = vec3::length(field);
    // `is_finite` first, so the comparison that follows is total by the time it
    // runs and NaN cannot slip through as "not greater than zero".
    if !length.is_finite() || !field_length.is_finite() {
        return R::ZERO;
    }
    if length <= R::ZERO || field_length <= R::ZERO {
        return R::ZERO;
    }
    vec3::dot(
        vec3::scale(geometric, length.recip()),
        vec3::scale(field, field_length.recip()),
    )
    .abs()
}

/// The three indices of one face.
fn face_indices<R: Real>(mesh: &MeshBuffer<R>, face: usize) -> [u32; 3] {
    [0, 1, 2].map(|k| mesh.indices[face * 3 + k])
}

/// A triangle's unnormalised geometric normal, from its winding.
fn geometric_normal<R: Real>(mesh: &MeshBuffer<R>, tri: [u32; 3]) -> [R; 3] {
    let p = tri.map(|i| mesh.positions[i as usize]);
    vec3::cross(vec3::sub(p[1], p[0]), vec3::sub(p[2], p[0]))
}

#[cfg(test)]
mod tests;
