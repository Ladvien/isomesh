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
pub mod interior;
pub mod table;
pub mod trilinear;

#[cfg(test)]
mod reference;
#[cfg(test)]
mod tests;

use alloc::vec::Vec;

use crate::cube::corner_offset;
use crate::vec3;
use crate::{MeshSink, Real, Sdf, Shape3};

pub use ambiguity::FaceAmbiguity;

/// How a cell resolves an ambiguous **interior** — the trilinear body saddle.
///
/// [`FaceAmbiguity`] decides what happens on a *face*. This decides what happens
/// inside the cell, which is a different and larger question: two faces can each
/// be resolved and the cell's topology still be undetermined, because the regions
/// they carry may or may not be joined through the interior.
///
/// Both settings are crack-free, for the same structural reason: neither touches
/// face connectivity. See [`trilinear::Contours`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum InteriorAmbiguity {
    /// Mesh each cell from its face connectivity alone. Marching Cubes and MC33's
    /// face half both do this, and it is the default because it is what every
    /// committed golden hash pins.
    #[default]
    Ignore,
    /// Ask the trilinear interpolant, via Grosso's construction: contours from the
    /// cut edges, topology from one quadratic, and a tunnel meshed **as** a tunnel.
    ///
    /// Only cells with an ambiguous face take this path; every other cell reads
    /// the same table it always did, so the cost is confined to the roughly one
    /// cell in two hundred that can differ.
    Trilinear,
}

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
/// use isomesh::marching_cubes::MarchingCubes;
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
    interior_ambiguity: InteriorAmbiguity,
    crossing_refinement: u32,
}

impl<R: Real> MarchingCubes<R> {
    /// A mesher that has allocated nothing yet.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            values: Vec::new(),
            edge_vertices: Vec::new(),
            face_ambiguity: FaceAmbiguity::Separate,
            interior_ambiguity: InteriorAmbiguity::Ignore,
            crossing_refinement: 0,
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

    /// How the cell's **interior** ambiguity is resolved.
    ///
    /// Defaults to [`InteriorAmbiguity::Ignore`]. [`InteriorAmbiguity::Trilinear`]
    /// is A-002b's, and it only does anything on a cell that also has an ambiguous
    /// face — so it is normally set together with
    /// [`FaceAmbiguity::AsymptoticDecider`], whose answer it builds on.
    pub fn set_interior_ambiguity(&mut self, interior_ambiguity: InteriorAmbiguity) {
        self.interior_ambiguity = interior_ambiguity;
    }

    /// Bisection steps spent locating each edge crossing on the real field.
    ///
    /// Ticket: F-007. Defaults to **0** — plain linear interpolation, which is
    /// what every committed golden hash pins.
    ///
    /// # When it changes anything
    ///
    /// `t = a / (a − b)` is exact for a field that is linear along the edge, and
    /// an analytic primitive is close enough over one cell that refinement moves
    /// the vertex by rounding. **A CSG field is not linear along an edge that
    /// crosses a seam**: `min`/`max` select an operand pointwise, so the field is
    /// two straight pieces meeting at an angle and a line through the endpoints
    /// misses the root.
    ///
    /// The sign is untouched — `{min(f,g) ≤ 0}` *is* the union — so the case
    /// classification and the topology are already correct, and this moves
    /// vertices without moving triangles. That is why it is a much narrower
    /// repair than redistancing the field.
    ///
    /// # Cost
    ///
    /// One field evaluation per step per **cut** edge, and only cut edges, so it
    /// scales with surface area rather than volume.
    pub fn set_crossing_refinement(&mut self, steps: u32) {
        self.crossing_refinement = steps;
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
    /// grid could produce more vertices than a `u32` can address — one per
    /// crossed grid edge, so three per sample, plus up to
    /// [`table::MAX_CENTROIDS`] cell-local cycle-centroid vertices per cell,
    /// which only the decider's joined cycles can add (A-015). Checked up
    /// front, which is what lets the per-vertex path stay a `debug_assert!`.
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
        // Three edge vertices per sample, plus up to `MAX_CENTROIDS` cell-local
        // cycle-centroid vertices per cell — A-015's fan centroids are created
        // per cell and never cached, so the edge count does not cover them.
        // Both shape types guarantee the sample product fits in `u32` and cells
        // are fewer than samples, so this sum stays far below `u64::MAX`.
        let cells = u64::from(size[0] - 1) * u64::from(size[1] - 1) * u64::from(size[2] - 1);
        // A cell's interior vertices are cell-local and uncached, so they are
        // budgeted per cell rather than per edge. A-015's cycle centroids need
        // three; A-002h's tunnel names all six vertices of the inner hexagon, so
        // the larger of the two is what has to be covered (M-218).
        let per_cell = if table::MAX_CENTROIDS > trilinear::MAX_INTERIOR_VERTICES {
            table::MAX_CENTROIDS
        } else {
            trilinear::MAX_INTERIOR_VERTICES
        };
        let bound = 3u64 * sample_count as u64 + per_cell as u64 * cells;
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
                    let mask = if ambiguous == 0 {
                        0
                    } else {
                        joined_mask(&corner_value, ambiguous)
                    };

                    // The trilinear path, and it is deliberately *only* reachable
                    // on a cell that already has an ambiguous face: everything
                    // else reads the same table it always did, byte for byte,
                    // which is what keeps every existing golden hash intact.
                    if ambiguous != 0 && self.interior_ambiguity == InteriorAmbiguity::Trilinear {
                        self.emit_trilinear(
                            sdf,
                            shape,
                            base,
                            case,
                            mask,
                            &corner_value,
                            origin,
                            cell_size,
                            out,
                        )?;
                        continue;
                    }

                    let entry = if ambiguous == 0 {
                        CASES[case as usize]
                    } else {
                        triangulate(segment_links(case, mask))
                    };
                    if entry.count == 0 {
                        continue;
                    }

                    // Cycle centroids first, because a triangle that names one
                    // needs every edge vertex of that cycle averaged before it
                    // can be emitted. They are cell-local by construction and so
                    // never cached — that locality is the whole point (A-015).
                    let mut centroid = [0u32; table::MAX_CENTROIDS];
                    for (c, slot) in centroid
                        .iter_mut()
                        .enumerate()
                        .take(entry.centroids as usize)
                    {
                        let code = table::CENTROID_BASE + c as u8;
                        let mut sum = [R::ZERO; 3];
                        let mut n = 0u32;
                        // A cycle's edges are exactly the non-centroid corners of
                        // the triangles naming it, each appearing twice.
                        for tri in &entry.triangles[..entry.count as usize] {
                            if tri[0] != code {
                                continue;
                            }
                            let position = edge_position(
                                sdf,
                                base,
                                tri[1],
                                &corner_value,
                                origin,
                                cell_size,
                                self.crossing_refinement,
                            );
                            sum = [
                                sum[0] + position[0],
                                sum[1] + position[1],
                                sum[2] + position[2],
                            ];
                            n += 1;
                        }
                        debug_assert!(n >= 4, "a centroid stands for a cycle of four or more");
                        let scale = R::from_f64(f64::from(n)).recip();
                        let position = [sum[0] * scale, sum[1] * scale, sum[2] * scale];
                        *slot = out.vertex(position, unit_gradient(sdf, position));
                    }

                    for tri in &entry.triangles[..entry.count as usize] {
                        let mut idx = [0u32; 3];
                        for (k, &code) in tri.iter().enumerate() {
                            debug_assert!(code != NO_EDGE);
                            idx[k] = if table::is_centroid(code) {
                                centroid[(code - table::CENTROID_BASE) as usize]
                            } else {
                                self.vertex_on_edge(
                                    sdf,
                                    shape,
                                    base,
                                    code,
                                    &corner_value,
                                    origin,
                                    cell_size,
                                    out,
                                )
                            };
                        }
                        out.triangle(idx[0], idx[1], idx[2]);
                    }
                }
            }
        }

        Ok(())
    }

    /// Mesh one ambiguous cell by Grosso's construction.
    ///
    /// Contours from the cut edges, topology from the body saddles, and a tunnel
    /// meshed **as** a tunnel. Interior vertices are cell-local and never cached —
    /// the same rule A-015's cycle centroids follow, and for the same reason: no
    /// other cell can name them.
    ///
    /// Normals come from the field's own gradient at the vertex, as everywhere
    /// else here. The reference implementation interpolates the eight corner
    /// normals instead; this crate has the field on hand and `unit_gradient` is
    /// its one rule for what a normal is.
    #[allow(clippy::too_many_arguments)]
    fn emit_trilinear<S, M>(
        &mut self,
        sdf: &S,
        shape: &impl Shape3,
        base: [u32; 3],
        case: u8,
        mask: u8,
        corner_value: &[R; 8],
        origin: [R; 3],
        cell_size: R,
        out: &mut M,
    ) -> crate::Result<()>
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        use trilinear::{BodySaddles, Contours, INTERIOR, MAX_INTERIOR_VERTICES, Topology};

        let contours = Contours::of(case, mask);
        if contours.count() == 0 {
            return Ok(());
        }
        let saddles = BodySaddles::of(corner_value);

        // **Refused before a single vertex is emitted.** Six body saddles with a
        // contour past Corollary 6's bound is not a tunnel — its inside region is
        // two blobs rather than one cylinder (M-229) — and neither §5.1's tunnel
        // rule nor §5.2's twelve-vertex rule applies to it. Sending it to
        // `fan_tunnel` anyway is what produced the hole M-228 found; the ring
        // count admitted a cell the corollary excludes, and the missing three-step
        // rule was the symptom rather than the defect. A-020b owns the
        // triangulation this needs.
        if contours.topology(&saddles) == Topology::SeparateDisks {
            return Err(crate::Error::UnresolvedSixSaddle {
                case,
                mask,
                longest: contours.longest(),
            });
        }

        // Cell-local coordinates to world. The cell is a unit cube in `(u,v,w)`.
        let to_world = |p: [R; 3]| {
            [
                origin[0] + cell_size * (R::from_f64(f64::from(base[0])) + p[0]),
                origin[1] + cell_size * (R::from_f64(f64::from(base[1])) + p[1]),
                origin[2] + cell_size * (R::from_f64(f64::from(base[2])) + p[2]),
            ]
        };

        // Interior vertices first, for the reason the centroids go first: a
        // triangle that names one cannot be emitted until it exists.
        let mut interior = [u32::MAX; MAX_INTERIOR_VERTICES];
        let hexagon = saddles.inner_hexagon();
        if let Some(ring) = hexagon {
            for (slot, local) in interior.iter_mut().zip(ring) {
                let position = to_world(local);
                *slot = out.vertex(position, unit_gradient(sdf, position));
            }
        } else if let Some(local) = saddles.interior_vertex() {
            debug_assert!(
                local.iter().all(|&c| c >= R::ZERO && c <= R::ONE),
                "interior vertex outside the cell: {local:?} mask {:#08b}",
                saddles.inside_mask()
            );
            let position = to_world(local);
            interior[0] = out.vertex(position, unit_gradient(sdf, position));
        }

        // `fan`/`fan_tunnel` hand back codes; resolve them here so the two paths
        // share one resolver and cannot disagree about what a code means.
        let mut triangles = [[0u8; 3]; trilinear::MAX_PATCH_TRIANGLES];
        let mut count = 0usize;
        if hexagon.is_some() {
            let unresolved = contours.fan_tunnel(&saddles, corner_value, |t| {
                triangles[count] = t;
                count += 1;
            });
            // **Loud rather than holed.** A contour edge whose endpoints land
            // three steps apart on the inner hexagon has no rule in Grosso's
            // construction and none in the authors' implementation, which simply
            // emits nothing there. Emitting nothing is a hole in the surface, and
            // a hole that only appears on Marching Cubes' case 13 with particular
            // face resolutions is exactly the kind of defect that reaches a
            // consumer's collider before anyone notices (M-228). A-020 owns
            // deriving the missing triangulation; until then this refuses.
            if unresolved != 0 {
                return Err(crate::Error::UnresolvedTunnel {
                    case,
                    mask,
                    edges: unresolved,
                });
            }
        } else {
            let fanned = saddles.interior_vertex().is_some();
            contours.fan(fanned, |t| {
                triangles[count] = t;
                count += 1;
            });
        }

        for tri in &triangles[..count] {
            let mut idx = [0u32; 3];
            for (k, &code) in tri.iter().enumerate() {
                idx[k] = if code >= INTERIOR {
                    let slot = interior[(code - INTERIOR) as usize];
                    debug_assert!(
                        slot != u32::MAX,
                        "a triangle named an interior vertex that was never created"
                    );
                    slot
                } else {
                    self.vertex_on_edge(
                        sdf,
                        shape,
                        base,
                        code,
                        corner_value,
                        origin,
                        cell_size,
                        out,
                    )
                };
            }
            out.triangle(idx[0], idx[1], idx[2]);
        }
        Ok(())
    }

    /// The vertex on one cut edge of one cell and where it sits, creating it if
    /// this is the first cell to ask.
    ///
    /// Cells sharing a grid edge share the vertex on it, which is what makes the
    /// result a connected surface. The cache is keyed on the grid edge — the
    /// lower sample plus the axis — so the key is the same whichever of the four
    /// surrounding cells arrives first, and the result does not depend on
    /// traversal order.
    ///
    /// The cache hit returns before any arithmetic, which is the point: a grid
    /// edge is asked for by up to four cells and only the first does work.
    /// Cycle centroids need positions rather than indices and take
    /// [`edge_position`] directly, so the position formula still lives in one
    /// place without putting it on this path.
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
        let axis = EDGE_AXIS[edge as usize] as usize;
        let lo_sample = corner_sample(shape, base, EDGE_CORNERS[edge as usize][0]);
        let key = lo_sample as usize * 3 + axis;

        let cached = self.edge_vertices[key];
        if cached != u32::MAX {
            return cached;
        }

        let position = edge_position(
            sdf,
            base,
            edge,
            corner_value,
            origin,
            cell_size,
            self.crossing_refinement,
        );
        let index = out.vertex(position, unit_gradient(sdf, position));
        self.edge_vertices[key] = index;
        index
    }
}

/// Where the surface crosses one cut edge of one cell.
///
/// The single definition of an edge vertex's position: [`MarchingCubes::extract`]
/// reaches it through the vertex cache, and a cycle centroid averages it over its
/// own cycle.
fn edge_position<S, R: Real>(
    sdf: &S,
    base: [u32; 3],
    edge: u8,
    corner_value: &[R; 8],
    origin: [R; 3],
    cell_size: R,
    refinement: u32,
) -> [R; 3]
where
    S: Sdf<Scalar = R>,
{
    let [lo_corner, hi_corner] = EDGE_CORNERS[edge as usize];
    let a = corner_value[lo_corner as usize];
    let b = corner_value[hi_corner as usize];
    // On a cut edge exactly one endpoint is strictly negative and the other is
    // >= 0, so `a - b` is never zero and no epsilon guard is needed. An epsilon
    // here would snap resolvable crossings to the midpoint.
    debug_assert!(is_inside(a) != is_inside(b));
    let t = a / (a - b);

    let lo_pos = corner_position(base, lo_corner, origin, cell_size);
    let hi_pos = corner_position(base, hi_corner, origin, cell_size);
    let t = refine_crossing(sdf, lo_pos, hi_pos, a, t, refinement);
    [
        lo_pos[0] + (hi_pos[0] - lo_pos[0]) * t,
        lo_pos[1] + (hi_pos[1] - lo_pos[1]) * t,
        lo_pos[2] + (hi_pos[2] - lo_pos[2]) * t,
    ]
}

/// Refine a linearly interpolated crossing by bisecting the *actual* field.
///
/// Ticket: F-007. Source: Pujol & Chica, *Adaptive approximation of signed
/// distance fields* (`10.1016/j.cag.2023.06.020`).
///
/// # What linear interpolation assumes, and where CSG breaks it
///
/// `t = a / (a − b)` is exact when `f` is linear along the edge, and every
/// analytic primitive here is close enough to linear over one cell that it does
/// not matter. **A CSG field is not.** `min`/`max` select an operand pointwise,
/// so along an edge that crosses a seam the field is a *kinked* piecewise
/// function: two straight pieces meeting at an angle. A line through its
/// endpoints misses the root by an amount that does not shrink with the cell
/// size the way a smooth field's error does — it shrinks only as the seam
/// occupies less of the edge.
///
/// **The sign is untouched by this.** `{min(f,g) ≤ 0}` *is* the union, exactly,
/// so the case classification and therefore the topology are already right.
/// Only the crossing's *position* is wrong, which is why this is a much narrower
/// repair than redistancing the field.
///
/// # Bisection rather than the CSG tree
///
/// The tree would say exactly where the kink is, and this crate's `Sdf` does not
/// expose one — a field is a closure as far as the extractor is concerned. So
/// the kink is found rather than known: bisect on the sign of the real field,
/// which converges on the true root whatever shape the field has between the
/// endpoints, at the cost of `steps` extra evaluations per cut edge.
///
/// Returns the parameter in `[0, 1]` along the edge.
fn refine_crossing<S, R>(sdf: &S, lo: [R; 3], hi: [R; 3], a: R, t0: R, steps: u32) -> R
where
    S: Sdf<Scalar = R>,
    R: Real,
{
    if steps == 0 {
        return t0;
    }
    // The bracket is the whole edge: the endpoints differ in sign by
    // construction, which is what makes bisection safe here without a guard.
    let (mut low, mut high) = (R::ZERO, R::ONE);
    let a_inside = is_inside(a);
    let mut t = t0;
    for _ in 0..steps {
        let p = [
            lo[0] + (hi[0] - lo[0]) * t,
            lo[1] + (hi[1] - lo[1]) * t,
            lo[2] + (hi[2] - lo[2]) * t,
        ];
        if is_inside(sdf.sample(p)) == a_inside {
            low = t;
        } else {
            high = t;
        }
        t = (low + high) * R::HALF;
    }
    t
}

/// The field's own gradient at a point, normalised — this crate's normal rule,
/// stated once so the edge vertices and the cycle centroids cannot use two.
///
/// # Panics
///
/// In debug builds, if the gradient vanishes. That means the field is degenerate
/// there; it cannot happen for any exact distance field, where `|grad|` is 1.
#[inline]
fn unit_gradient<R: Real, S: Sdf<Scalar = R>>(sdf: &S, position: [R; 3]) -> [R; 3] {
    let g = sdf.gradient(position);
    let len = vec3::length(g);
    debug_assert!(len > R::ZERO, "zero gradient at a surface vertex");
    vec3::scale(g, len.recip())
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
    /// Triangles naming a centroid this case does not declare, or a declared
    /// centroid standing for fewer than four triangles.
    ///
    /// A cycle centroid is cell-local, so a stale or spurious reference would
    /// silently attach geometry to the wrong vertex. See
    /// [`table::CENTROID_BASE`].
    pub bad_centroid_reference: u64,
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
            && self.bad_centroid_reference == 0
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
    let mut centroid_uses = [0u32; table::MAX_CENTROIDS];
    for tri in &entry.triangles[..entry.count as usize] {
        if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
            report.degenerate_triangles += 1;
        }
        for &code in tri {
            if table::is_centroid(code) {
                let c = (code - table::CENTROID_BASE) as usize;
                if c >= entry.centroids as usize {
                    report.bad_centroid_reference += 1;
                } else {
                    centroid_uses[c] += 1;
                }
                continue;
            }
            if !cut[code as usize] {
                report.triangles_on_uncut_edges += 1;
            }
            used[code as usize] = true;
        }
    }
    if used != cut {
        report.cut_edge_mismatch += 1;
    }
    // A centroid stands for a cycle of four or more, and a cycle of `k` gives it
    // exactly `k` triangles. Fewer than four means it was emitted for a cycle
    // that did not need one.
    for &uses in centroid_uses.iter().take(entry.centroids as usize) {
        if uses < 4 {
            report.bad_centroid_reference += 1;
        }
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
