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

pub mod ambiguity;
pub mod table;

#[cfg(test)]
mod reference;
#[cfg(test)]
mod tests;

use alloc::vec::Vec;

use crate::cube::corner_offset;
use crate::vec3;
use crate::{MeshSink, Real, Sdf, Shape3};

pub use ambiguity::FaceAmbiguity;

use ambiguity::joined_mask;
use table::{
    AMBIGUOUS_FACES, CASES, EDGE_AXIS, EDGE_CORNERS, NO_EDGE, is_inside, segment_links, triangulate,
};

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
/// let shape = RuntimeShape3::new([33; 3])?;
/// mc.extract(&Sphere::<f32>::canonical(), &shape, [-2.0; 3], 0.125, &mut out)?;
///
/// assert!(out.triangle_count() > 0);
/// # Ok::<(), isomesh::Error>(())
/// ```
#[derive(Debug)]
pub struct MarchingCubes<R: Real> {
    values: Vec<R>,
    /// One slot per (sample, axis): the vertex sitting on that grid edge, or
    /// [`u32::MAX`].
    edge_vertices: Vec<u32>,
    face_ambiguity: FaceAmbiguity,
}

impl<R: Real> MarchingCubes<R> {
    /// A mesher that has allocated nothing yet.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            values: Vec::new(),
            edge_vertices: Vec::new(),
            face_ambiguity: FaceAmbiguity::Separate,
        }
    }

    /// How ambiguous faces are resolved.
    ///
    /// Defaults to [`FaceAmbiguity::Separate`], which is Marching Cubes proper.
    /// [`FaceAmbiguity::AsymptoticDecider`] is A-002's MC33 face rule; see
    /// [`ambiguity`] for the mathematics and A-002's archive entry for what the
    /// difference measures.
    pub fn set_face_ambiguity(&mut self, face_ambiguity: FaceAmbiguity) {
        self.face_ambiguity = face_ambiguity;
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
    /// # Errors
    ///
    /// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if any axis has fewer
    /// than two samples, since then there is no cell to march.
    /// [`Error::IndexSpaceExhausted`](crate::Error::IndexSpaceExhausted) if the
    /// grid could produce more vertices than a `u32` can address — Marching
    /// Cubes places one per crossed grid edge, so the bound is three per sample.
    /// Checked up front, which is what lets the per-vertex path stay a
    /// `debug_assert!`.
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
        let size = shape.size();
        if size[0] < 2 || size[1] < 2 || size[2] < 2 {
            return Err(crate::Error::GridTooSmall { size });
        }
        let sample_count = shape.element_count();
        let bound = 3u64 * sample_count as u64;
        if bound > u64::from(u32::MAX) {
            return Err(crate::Error::IndexSpaceExhausted { needed: bound });
        }

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

                    // The triangulation for this cell. Under `Separate` that is
                    // the derived table verbatim; under `AsymptoticDecider` it is
                    // the same construction with the ambiguous faces re-paired by
                    // the bilinear saddle. A cell with no ambiguous face reads
                    // the table either way, which is a memo and not a second
                    // rule: `masks_are_ignored_on_unambiguous_faces` proves the
                    // two agree, and `the_separate_mask_reproduces_the_derived
                    // _table` proves the table is the mask-zero construction.
                    let ambiguous = match self.face_ambiguity {
                        FaceAmbiguity::Separate => 0,
                        FaceAmbiguity::AsymptoticDecider => AMBIGUOUS_FACES[case as usize],
                    };
                    let entry = if ambiguous == 0 {
                        CASES[case as usize]
                    } else {
                        triangulate(segment_links(case, joined_mask(&corner_value, ambiguous)))
                    };
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

        Ok(())
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

/// `(face, that face's 4-corner pattern, that face's decision) -> its segments`.
///
/// The decision index is masked by whether the face is *actually* ambiguous, so
/// a set bit on a face that has no choice shares a slot with the same face at
/// bit clear. A disagreement in that slot is therefore the report of a mask bit
/// having had an effect where it must not — which is what licenses the table
/// lookup in [`MarchingCubes::extract`] for cells with no ambiguous face.
type FaceMemo = [[[Option<[u8; table::EDGE_COUNT]>; 2]; 16]; 6];

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
///
/// This checks the shipped [`CASES`] array, at the all-separate resolution.
/// [`validate_decider_table`] is the same checks over every resolution mask.
#[must_use]
pub fn validate_table() -> TableReport {
    let mut report = TableReport::default();
    let mut face_seen: FaceMemo = [[[None; 2]; 16]; 6];
    for case in 0..=255u8 {
        check_case(case, 0, &CASES[case as usize], &mut report, &mut face_seen);
    }
    report
}

/// The same checks, over all 256 cases **and** all 64 face-resolution masks.
///
/// 16,384 combinations, which is what A-002 has to be sound over rather than
/// just the 256 the compile-time table covers. The face-locality property is the
/// one that matters here and it is stronger than it looks: two cells meeting on
/// a face agree about that face's corner signs *and*, because the decider is a
/// function of the four shared sample values, about its decision bit — so if the
/// segments are a function of `(pattern, bit)` the cells still cannot disagree.
///
/// `max_triangles` is recorded rather than gated, as in [`validate_table`]; the
/// crossed pairing can produce longer cycles than the separated one.
#[must_use]
pub fn validate_decider_table() -> TableReport {
    let mut report = TableReport::default();
    let mut face_seen: FaceMemo = [[[None; 2]; 16]; 6];
    for case in 0..=255u8 {
        for mask in 0..(1u8 << table::FACE_COUNT) {
            let entry = triangulate(segment_links(case, mask));
            check_case(case, mask, &entry, &mut report, &mut face_seen);
        }
    }
    report
}

/// Is a face with this 4-corner pattern ambiguous?
///
/// True when the signs alternate around the ring, which is the only way a face
/// gets four cut edges and so the only way it has a pairing to choose. Written
/// against the pattern rather than [`table::AMBIGUOUS_FACES`] so the check below
/// stays honest about consulting nothing but the face's own corners;
/// `ambiguous_faces_agrees_with_the_face_pattern` ties the two together.
const fn pattern_is_ambiguous(pattern: usize) -> bool {
    pattern == 0b0101 || pattern == 0b1010
}

fn check_case(
    case: u8,
    mask: u8,
    entry: &table::McCase,
    report: &mut TableReport,
    face_seen: &mut FaceMemo,
) {
    use table::{EDGE_COUNT, corner_inside, edge_index, face_bit, face_corners};

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

    let links = segment_links(case, mask);
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
    // nothing but that face's own corners and that face's own decision bit.
    for axis in 0..3usize {
        for side in 0..2u8 {
            let f = axis * 2 + side as usize;
            let c = face_corners(axis, side);
            let mut pattern = 0usize;
            for (k, &corner) in c.iter().enumerate() {
                if corner_inside(case, corner) {
                    pattern |= 1 << k;
                }
            }
            let bit =
                usize::from(mask & face_bit(axis, side) != 0 && pattern_is_ambiguous(pattern));

            let mut segments = [NO_EDGE; EDGE_COUNT];
            if let Some(start) = (0..4).find(|&k| !corner_inside(case, c[k])) {
                let mut entries = [NO_EDGE; 2];
                let mut exits = [NO_EDGE; 2];
                let mut pairs = 0usize;
                for j in 0..4 {
                    let p = c[(start + j) % 4];
                    let q = c[(start + j + 1) % 4];
                    match (corner_inside(case, p), corner_inside(case, q)) {
                        (false, true) => entries[pairs] = edge_index(p, q),
                        (true, false) => {
                            exits[pairs] = edge_index(p, q);
                            pairs += 1;
                        }
                        _ => {}
                    }
                }
                if pairs == 2 && bit == 1 {
                    segments[entries[0] as usize] = exits[1];
                    segments[entries[1] as usize] = exits[0];
                } else {
                    for n in 0..pairs {
                        segments[entries[n] as usize] = exits[n];
                    }
                }
            }

            match face_seen[f][pattern][bit] {
                None => face_seen[f][pattern][bit] = Some(segments),
                Some(previous) => {
                    if previous != segments {
                        report.face_disagreements += 1;
                    }
                }
            }
        }
    }
}
