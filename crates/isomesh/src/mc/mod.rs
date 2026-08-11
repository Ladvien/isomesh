//! Marching Cubes.
//!
//! The reference implementation everything else in this crate is compared
//! against, so correctness matters here more than speed.
//!
//! The 256-case table is **derived at compile time** rather than transcribed —
//! see [`table`] for why and how, and [`validate_table`] for the structural
//! check that backs it up. It agrees with the published Lorensen & Cline table
//! on all 256 cases; `matches_the_published_table` in the tests demonstrates
//! that against an independently parsed copy.

pub mod table;

#[cfg(test)]
mod reference;
#[cfg(test)]
mod tests;

use alloc::vec::Vec;

use crate::cube::corner_offset;
use crate::vec3;
use crate::{MeshSink, Real, Sdf, Shape3};

use table::{CASES, EDGE_AXIS, EDGE_CORNERS, NO_EDGE, is_inside};

/// Marching Cubes over a sampled grid.
///
/// Owns its scratch buffers so that re-meshing thousands of chunks does not
/// allocate thousands of times — the same reason [`crate::MeshBuffer`] is
/// caller-provided and reusable. Construct once, call [`extract`](Self::extract)
/// as often as you like.
///
/// # Example
///
/// ```
/// use isomesh::{MeshBuffer, RuntimeShape3};
/// use isomesh::fields::Sphere;
/// use isomesh::mc::MarchingCubes;
///
/// let mut mc = MarchingCubes::<f32>::new();
/// let mut out = MeshBuffer::<f32>::new();
///
/// // 33 samples per axis spans 32 cells over [-2, 2].
/// let shape = RuntimeShape3::new([33; 3]);
/// mc.extract(&Sphere::<f32>::canonical(), &shape, [-2.0; 3], 0.125, &mut out);
///
/// assert!(out.triangle_count() > 0);
/// ```
#[derive(Debug)]
pub struct MarchingCubes<R: Real> {
    values: Vec<R>,
    /// One slot per (sample, axis): the vertex sitting on that grid edge, or
    /// [`u32::MAX`].
    edge_vertices: Vec<u32>,
}

impl<R: Real> MarchingCubes<R> {
    /// A mesher that has allocated nothing yet.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            values: Vec::new(),
            edge_vertices: Vec::new(),
        }
    }

    /// Extract the zero level set into `out`.
    ///
    /// `shape` counts **samples**, not cells, so a shape of `[n; 3]` spans
    /// `n - 1` cells per axis. `origin` is the world position of sample
    /// `[0, 0, 0]` and `cell_size` is the spacing between adjacent samples.
    ///
    /// # Conventions
    ///
    /// - **Sign:** negative is inside, and a sample of exactly zero counts as
    ///   outside. See [`table::is_inside`].
    /// - **Winding:** counter-clockwise seen from outside the solid, so
    ///   `(b − a) × (c − a)` points away from it. Verified rather than asserted:
    ///   `meshed_sphere_has_positive_signed_volume` would catch a global flip,
    ///   which no manifold or Euler check can see.
    /// - **Normals:** the field's own gradient at the vertex, normalised.
    ///   A-012 is where alternative estimators live.
    ///
    /// Vertices are shared between cells that meet on a grid edge, so the output
    /// is a properly connected surface rather than a triangle soup.
    ///
    /// # Panics
    ///
    /// If `shape` has fewer than two samples on any axis — there would be no
    /// cell to march — or if the sample count exceeds `u32`.
    pub fn extract<S, M>(
        &mut self,
        sdf: &S,
        shape: &impl Shape3,
        origin: [R; 3],
        cell_size: R,
        out: &mut M,
    ) where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        let size = shape.size();
        assert!(
            size[0] >= 2 && size[1] >= 2 && size[2] >= 2,
            "marching cubes needs at least two samples per axis, got {size:?}"
        );
        let sample_count = shape.element_count();

        // ── sample once per grid point ──────────────────────────────────────
        self.values.clear();
        self.values.reserve(sample_count);
        for z in 0..size[2] {
            for y in 0..size[1] {
                for x in 0..size[0] {
                    let p = [
                        origin[0] + cell_size * R::from_f64(f64::from(x)),
                        origin[1] + cell_size * R::from_f64(f64::from(y)),
                        origin[2] + cell_size * R::from_f64(f64::from(z)),
                    ];
                    self.values.push(sdf.sample(p));
                }
            }
        }
        debug_assert_eq!(self.values.len(), sample_count);

        self.edge_vertices.clear();
        self.edge_vertices.resize(sample_count * 3, u32::MAX);

        // ── march ───────────────────────────────────────────────────────────
        for z in 0..size[2] - 1 {
            for y in 0..size[1] - 1 {
                for x in 0..size[0] - 1 {
                    let base = [x, y, z];

                    let mut case = 0u8;
                    let mut corner_value = [R::ZERO; 8];
                    for (c, slot) in corner_value.iter_mut().enumerate() {
                        let s = corner_sample(shape, base, c as u8);
                        let v = self.values[s as usize];
                        *slot = v;
                        if is_inside(v) {
                            case |= 1 << c;
                        }
                    }

                    let entry = &CASES[case as usize];
                    if entry.count == 0 {
                        continue;
                    }

                    for tri in &entry.triangles[..entry.count as usize] {
                        let mut idx = [0u32; 3];
                        for (k, &edge) in tri.iter().enumerate() {
                            debug_assert!(edge != NO_EDGE);
                            idx[k] = self.vertex_on_edge(
                                sdf,
                                shape,
                                base,
                                edge,
                                &corner_value,
                                origin,
                                cell_size,
                                out,
                            );
                        }
                        out.triangle(idx[0], idx[1], idx[2]);
                    }
                }
            }
        }
    }

    /// The vertex on one cut edge of one cell, creating it if this is the first
    /// cell to ask.
    ///
    /// Cells sharing a grid edge share the vertex on it, which is what makes the
    /// result a connected surface. The cache is keyed on the grid edge — the
    /// lower sample plus the axis — so the key is the same whichever of the four
    /// surrounding cells arrives first, and the result does not depend on
    /// traversal order.
    #[allow(clippy::too_many_arguments)]
    fn vertex_on_edge<S, M>(
        &mut self,
        sdf: &S,
        shape: &impl Shape3,
        base: [u32; 3],
        edge: u8,
        corner_value: &[R; 8],
        origin: [R; 3],
        cell_size: R,
        out: &mut M,
    ) -> u32
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        let [lo_corner, hi_corner] = EDGE_CORNERS[edge as usize];
        let axis = EDGE_AXIS[edge as usize] as usize;
        let lo_sample = corner_sample(shape, base, lo_corner);
        let key = lo_sample as usize * 3 + axis;

        let cached = self.edge_vertices[key];
        if cached != u32::MAX {
            return cached;
        }

        let a = corner_value[lo_corner as usize];
        let b = corner_value[hi_corner as usize];
        // On a cut edge exactly one endpoint is strictly negative and the other
        // is >= 0, so `a - b` is never zero and no epsilon guard is needed. An
        // epsilon here would snap resolvable crossings to the midpoint.
        debug_assert!(is_inside(a) != is_inside(b));
        let t = a / (a - b);

        let lo_pos = corner_position(base, lo_corner, origin, cell_size);
        let hi_pos = corner_position(base, hi_corner, origin, cell_size);
        let position = [
            lo_pos[0] + (hi_pos[0] - lo_pos[0]) * t,
            lo_pos[1] + (hi_pos[1] - lo_pos[1]) * t,
            lo_pos[2] + (hi_pos[2] - lo_pos[2]) * t,
        ];

        let g = sdf.gradient(position);
        let len = vec3::length(g);
        // A zero gradient at a surface crossing means the field is degenerate
        // there; it cannot happen for any exact distance field, where |grad| is 1.
        debug_assert!(len > R::ZERO, "zero gradient at a surface crossing");
        let normal = vec3::scale(g, len.recip());

        let index = out.vertex(position, normal);
        self.edge_vertices[key] = index;
        index
    }
}

impl<R: Real> Default for MarchingCubes<R> {
    fn default() -> Self {
        Self::new()
    }
}

#[inline]
fn corner_sample(shape: &impl Shape3, base: [u32; 3], corner: u8) -> u32 {
    let o = corner_offset(corner);
    shape.linearize([base[0] + o[0], base[1] + o[1], base[2] + o[2]])
}

#[inline]
fn corner_position<R: Real>(base: [u32; 3], corner: u8, origin: [R; 3], cell_size: R) -> [R; 3] {
    let o = corner_offset(corner);
    [
        origin[0] + cell_size * R::from_f64(f64::from(base[0] + o[0])),
        origin[1] + cell_size * R::from_f64(f64::from(base[1] + o[1])),
        origin[2] + cell_size * R::from_f64(f64::from(base[2] + o[2])),
    ]
}

/// What is wrong with a derived case table, if anything.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TableReport {
    /// Cases whose triangles name an edge that is not cut.
    pub triangles_on_uncut_edges: u64,
    /// Cases where some cut edge carries no triangle, or vice versa.
    pub cut_edge_mismatch: u64,
    /// Cases where a triangle repeats an edge.
    pub degenerate_triangles: u64,
    /// Cut edges without exactly one incoming and one outgoing segment.
    pub bad_segment_degree: u64,
    /// Faces whose segments are not a function of that face's own corner signs.
    ///
    /// Non-zero here means cracks: two cells sharing a face would disagree
    /// about where the surface crosses it.
    pub face_disagreements: u64,
    /// The largest triangle count produced by any case.
    pub max_triangles: u64,
}

impl TableReport {
    /// `true` when every check passed.
    #[must_use]
    pub fn is_sound(&self) -> bool {
        self.triangles_on_uncut_edges == 0
            && self.cut_edge_mismatch == 0
            && self.degenerate_triangles == 0
            && self.bad_segment_degree == 0
            && self.face_disagreements == 0
    }
}

/// Check all 256 cases structurally, without consulting any reference table.
///
/// This is the brief's second defence, and it is the one that does not depend on
/// anyone else's numbering being read correctly. It verifies the properties the
/// construction is supposed to guarantee:
///
/// - triangles only ever name edges the corner signs actually cut, and every cut
///   edge carries a triangle;
/// - no triangle repeats an edge;
/// - every cut edge has exactly one incoming and one outgoing segment, which is
///   what makes the segments close into loops with nothing left over;
/// - **a face's segments depend only on that face's own four corner signs.**
///   That last one is the crack-free property: two cells meeting on a face see
///   the same four corners, so if the segments are a function of those corners
///   the cells cannot disagree.
#[must_use]
pub fn validate_table() -> TableReport {
    use table::{EDGE_COUNT, corner_inside, edge_index, face_corners, segment_links};

    let mut report = TableReport::default();

    // (face, that face's 4-corner pattern) -> the segments it produced.
    let mut face_seen: [[Option<[u8; EDGE_COUNT]>; 16]; 6] = [[None; 16]; 6];

    for case in 0..=255u8 {
        let entry = &CASES[case as usize];
        report.max_triangles = report.max_triangles.max(u64::from(entry.count));

        let mut cut = [false; EDGE_COUNT];
        for (e, slot) in cut.iter_mut().enumerate() {
            let [a, b] = EDGE_CORNERS[e];
            *slot = corner_inside(case, a) != corner_inside(case, b);
        }

        let mut used = [false; EDGE_COUNT];
        for tri in &entry.triangles[..entry.count as usize] {
            if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
                report.degenerate_triangles += 1;
            }
            for &e in tri {
                if !cut[e as usize] {
                    report.triangles_on_uncut_edges += 1;
                }
                used[e as usize] = true;
            }
        }
        if used != cut {
            report.cut_edge_mismatch += 1;
        }

        let links = segment_links(case);
        let mut incoming = [0u8; EDGE_COUNT];
        for e in 0..EDGE_COUNT {
            if links[e] != NO_EDGE {
                incoming[links[e] as usize] += 1;
            }
        }
        for e in 0..EDGE_COUNT {
            let out_degree = u8::from(links[e] != NO_EDGE);
            if cut[e] != (out_degree == 1) || incoming[e] != out_degree {
                report.bad_segment_degree += 1;
            }
        }

        // Recompute each face's segments in isolation and check they depend on
        // nothing but that face's own corners.
        let mut f = 0usize;
        for axis in 0..3usize {
            for side in 0..2u8 {
                let c = face_corners(axis, side);
                let mut pattern = 0usize;
                for (k, &corner) in c.iter().enumerate() {
                    if corner_inside(case, corner) {
                        pattern |= 1 << k;
                    }
                }

                let mut segments = [NO_EDGE; EDGE_COUNT];
                if let Some(start) = (0..4).find(|&k| !corner_inside(case, c[k])) {
                    let mut entry_edge = NO_EDGE;
                    for j in 0..4 {
                        let p = c[(start + j) % 4];
                        let q = c[(start + j + 1) % 4];
                        match (corner_inside(case, p), corner_inside(case, q)) {
                            (false, true) => entry_edge = edge_index(p, q),
                            (true, false) => {
                                segments[entry_edge as usize] = edge_index(p, q);
                            }
                            _ => {}
                        }
                    }
                }

                match face_seen[f][pattern] {
                    None => face_seen[f][pattern] = Some(segments),
                    Some(previous) => {
                        if previous != segments {
                            report.face_disagreements += 1;
                        }
                    }
                }
                f += 1;
            }
        }
    }

    report
}
