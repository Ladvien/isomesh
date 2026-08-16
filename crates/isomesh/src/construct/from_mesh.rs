//! A distance field from a triangle mesh, signed by the angle-weighted
//! pseudonormal.
//!
//! Ticket: S-006. Bærentzen & Aanæs, *Signed Distance Computation Using the
//! Angle Weighted Pseudonormal*, IEEE TVCG 11(3), pp. 243–253 (2005)
//! (`10.1109/TVCG.2005.49`).
//!
//! # This is a proof, not a heuristic
//!
//! The paper's Theorem 1, verified in the corpus rather than recalled. With `c`
//! a closest point on the mesh to `p` and `r = p − c`:
//!
//! ```text
//! n_α · r > 0   if p is outside
//! n_α · r < 0   if p is inside
//! n_α · r = 0   if p is on the surface
//! ```
//!
//! where the **angle-weighted pseudonormal** at a point `x` on the mesh is
//!
//! ```text
//! n_α = Σ αᵢ nᵢ / ‖Σ αᵢ nᵢ‖
//! ```
//!
//! over the faces incident on `x`, with `αᵢ` the incident angle. The paper is
//! explicit about what that collapses to, which is why the three cases here are
//! transcribed and not invented:
//!
//! - **face** — one incident face, `α = 2π`, so `n_α` *is* the face normal;
//! - **edge** — *"both face normals have weight π and the result is the same as
//!   when computing the unweighted average"*;
//! - **vertex** — the incident angles at that vertex, which is Thürmer &
//!   Wüthrich's original construction.
//!
//! And, quoting: *"we do not assume that the closest point is unique. The proof
//! requires only that `c` is a closest point"* — so this is correct on the medial
//! axis too, which is exactly where a ray-parity or plain-normal test fails.
//!
//! The paper is equally explicit that the obvious alternatives do **not** work:
//! the unweighted mean of normals (Gouraud) and the least-squares plane fit
//! (Glassner) *"cannot be used for sign computation in general"*. Substituting
//! either would produce a field that is right almost everywhere, which is the
//! worst failure mode available.
//!
//! # Cost, measured rather than asserted
//!
//! Every sample is tested against every triangle, with two exact rejects layered
//! over that: triangles are grouped into blocks of 64 with a block
//! bounding box, and each triangle carries its own. A box whose nearest corner
//! is further than the best distance so far cannot contain a closer point, so
//! the test is exact — it changes what is *visited*, never what is *found*, and
//! `the_reject_agrees_with_brute_force_and_is_faster` asserts bit-identical output
//! against an unaccelerated scan.
//!
//! **A uniform grid over the sample cells was implemented first and was 3.9×
//! slower (M-260).** Expanding shells around a sample cost `O(k³)` bins to reach
//! radius `k`, and a grid whose corners are twelve cells from the surface makes
//! every one of those samples walk most of the grid. The block reject wins
//! because it is the far samples — the majority — that reject fastest.
//!
//! It is **not** a BVH. The paper uses one, and a BVH is what makes this
//! `O(log m)` rather than `O(m)` per sample. The blocks work here because
//! Marching Cubes emits triangles in grid order, so consecutive ones are
//! spatially close and their block box is tight; a mesh whose triangles arrive
//! in random order degrades to the unaccelerated scan and stays correct.
//!

#[cfg(test)]
mod tests;

use alloc::collections::BTreeMap;
use alloc::vec;
use alloc::vec::Vec;

use crate::real::Real;
use crate::sdf::Sdf;
use crate::shape::Shape3;

/// Which feature of a triangle a closest point landed on.
///
/// The whole reason this is returned rather than discarded: the pseudonormal
/// differs per feature, and picking the face normal for a point closest to a
/// vertex is precisely the error the paper exists to correct.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Feature {
    /// Interior of the face.
    Face,
    /// Edge between local vertices `i` and `(i + 1) % 3`.
    Edge(u8),
    /// Local vertex `i`.
    Vertex(u8),
}

/// `a − b`.
pub(super) fn sub<R: Real>(a: [R; 3], b: [R; 3]) -> [R; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// `a + b`.
pub(super) fn add<R: Real>(a: [R; 3], b: [R; 3]) -> [R; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// `a · s`.
pub(super) fn scale<R: Real>(a: [R; 3], s: R) -> [R; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// `a · b`.
pub(super) fn dot<R: Real>(a: [R; 3], b: [R; 3]) -> R {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `a × b`.
pub(super) fn cross<R: Real>(a: [R; 3], b: [R; 3]) -> [R; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// `‖a‖`.
pub(super) fn norm<R: Real>(a: [R; 3]) -> R {
    dot(a, a).sqrt()
}

/// Closest point on triangle `abc` to `p`, and which feature it landed on.
///
/// Ericson, *Real-Time Collision Detection*, §5.1.5 — the barycentric region
/// test, which decides the feature as a by-product of finding the point instead
/// of needing a second classification pass that could disagree with the first.
pub(super) fn closest_on_triangle<R: Real>(
    p: [R; 3],
    a: [R; 3],
    b: [R; 3],
    c: [R; 3],
) -> ([R; 3], Feature) {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= R::ZERO && d2 <= R::ZERO {
        return (a, Feature::Vertex(0));
    }

    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= R::ZERO && d4 <= d3 {
        return (b, Feature::Vertex(1));
    }

    let vc = d1 * d4 - d3 * d2;
    if vc <= R::ZERO && d1 >= R::ZERO && d3 <= R::ZERO {
        let denom = d1 - d3;
        let v = if denom == R::ZERO {
            R::ZERO
        } else {
            d1 / denom
        };
        return (add(a, scale(ab, v)), Feature::Edge(0));
    }

    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= R::ZERO && d5 <= d6 {
        return (c, Feature::Vertex(2));
    }

    let vb = d5 * d2 - d1 * d6;
    if vb <= R::ZERO && d2 >= R::ZERO && d6 <= R::ZERO {
        let denom = d2 - d6;
        let w = if denom == R::ZERO {
            R::ZERO
        } else {
            d2 / denom
        };
        return (add(a, scale(ac, w)), Feature::Edge(2));
    }

    let va = d3 * d6 - d5 * d4;
    if va <= R::ZERO && (d4 - d3) >= R::ZERO && (d5 - d6) >= R::ZERO {
        let denom = (d4 - d3) + (d5 - d6);
        let w = if denom == R::ZERO {
            R::ZERO
        } else {
            (d4 - d3) / denom
        };
        return (add(b, scale(sub(c, b), w)), Feature::Edge(1));
    }

    let denom = va + vb + vc;
    if denom == R::ZERO {
        // Degenerate triangle: every barycentric region collapsed. Its vertices
        // are still on the surface, so `a` is a correct closest point even
        // though it is not the nearest one -- and a degenerate triangle's
        // pseudonormal contributes nothing to the sums anyway.
        return (a, Feature::Vertex(0));
    }
    let inv = R::ONE / denom;
    let v = vb * inv;
    let w = vc * inv;
    (add(add(a, scale(ab, v)), scale(ac, w)), Feature::Face)
}

/// Angle at `x` in the triangle `x, y, z`, in radians.
///
/// The `αᵢ` of the paper. Returns zero for a degenerate corner rather than a
/// NaN, so a sliver triangle contributes nothing instead of poisoning the sum.
fn corner_angle<R: Real>(x: [R; 3], y: [R; 3], z: [R; 3]) -> R {
    let u = sub(y, x);
    let v = sub(z, x);
    let lu = norm(u);
    let lv = norm(v);
    if lu == R::ZERO || lv == R::ZERO {
        return R::ZERO;
    }
    let mut c = dot(u, v) / (lu * lv);
    if c > R::ONE {
        c = R::ONE;
    }
    if c < -R::ONE {
        c = -R::ONE;
    }
    c.acos()
}

/// The three pseudonormal tables the sign test reads.
///
/// Unnormalised sums throughout: only the **sign** of `n_α · r` is used, and
/// dividing by `‖Σ αᵢ nᵢ‖` cannot change it. Skipping the normalisation removes
/// the one place a zero-length sum could produce a NaN.
#[derive(Debug)]
struct Pseudonormals<R: Real> {
    /// Per triangle, the unit face normal.
    face: Vec<[R; 3]>,
    /// Per undirected edge `(min, max)`, the sum of its incident face normals.
    edge: BTreeMap<(u32, u32), [R; 3]>,
    /// Per vertex, `Σ αᵢ nᵢ` over incident faces.
    vertex: Vec<[R; 3]>,
}

impl<R: Real> Pseudonormals<R> {
    /// Build the tables in one pass over the triangles.
    fn build(positions: &[[R; 3]], indices: &[u32]) -> Self {
        let mut face = Vec::with_capacity(indices.len() / 3);
        let mut edge: BTreeMap<(u32, u32), [R; 3]> = BTreeMap::new();
        let mut vertex = vec![[R::ZERO; 3]; positions.len()];

        for tri in indices.chunks_exact(3) {
            let (ia, ib, ic) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
            let (a, b, c) = (positions[ia], positions[ib], positions[ic]);
            let n = cross(sub(b, a), sub(c, a));
            let len = norm(n);
            // A degenerate triangle has no normal and no incident angles worth
            // the name. It contributes zero to every sum, which is the only
            // answer that leaves the other faces' contributions intact.
            let unit = if len == R::ZERO {
                [R::ZERO; 3]
            } else {
                scale(n, R::ONE / len)
            };
            face.push(unit);

            for (x, y, z) in [(ia, ib, ic), (ib, ic, ia), (ic, ia, ib)] {
                let angle = corner_angle(positions[x], positions[y], positions[z]);
                vertex[x] = add(vertex[x], scale(unit, angle));
            }

            for (u, v) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if u < v { (u, v) } else { (v, u) };
                let slot = edge.entry(key).or_insert([R::ZERO; 3]);
                *slot = add(*slot, unit);
            }
        }

        Self { face, edge, vertex }
    }

    /// The pseudonormal at the closest point of triangle `t`, per feature.
    fn at(&self, indices: &[u32], t: usize, feature: Feature) -> [R; 3] {
        let tri = &indices[t * 3..t * 3 + 3];
        match feature {
            Feature::Face => self.face[t],
            Feature::Vertex(i) => self.vertex[tri[i as usize] as usize],
            Feature::Edge(i) => {
                let u = tri[i as usize];
                let v = tri[(i as usize + 1) % 3];
                let key = if u < v { (u, v) } else { (v, u) };
                // An edge missing from the table cannot happen -- it was built
                // from these same indices -- but falling back to the face normal
                // would be a second path. The face normal *is* the sum when only
                // one face is incident, which is what a boundary edge means, and
                // that is what the table already holds.
                self.edge.get(&key).copied().unwrap_or(self.face[t])
            }
        }
    }
}

/// Triangles per block. One cache line of `f32` indices, and small enough that a
/// block's bounding box stays tight for grid-ordered input.
const BLOCK: usize = 64;

/// Axis-aligned bounds, per triangle and per block of [`BLOCK`] triangles.
#[derive(Debug)]
struct Bounds<R: Real> {
    /// `[lo, hi]` per triangle.
    tri: Vec<[[R; 3]; 2]>,
    /// `[lo, hi]` per block.
    block: Vec<[[R; 3]; 2]>,
}

/// Squared distance from `p` to the box `[lo, hi]`. Zero inside.
///
/// Squared because it is only ever compared against another squared distance,
/// and the square root would be a per-triangle cost paid to learn nothing.
fn box_distance_sq<R: Real>(p: [R; 3], b: &[[R; 3]; 2]) -> R {
    let mut total = R::ZERO;
    for axis in 0..3 {
        let d = if p[axis] < b[0][axis] {
            b[0][axis] - p[axis]
        } else if p[axis] > b[1][axis] {
            p[axis] - b[1][axis]
        } else {
            R::ZERO
        };
        total += d * d;
    }
    total
}

impl<R: Real> Bounds<R> {
    fn build(positions: &[[R; 3]], indices: &[u32]) -> Self {
        let triangles = indices.len() / 3;
        let mut tri = Vec::with_capacity(triangles);
        for t in 0..triangles {
            let v = &indices[t * 3..t * 3 + 3];
            let mut lo = positions[v[0] as usize];
            let mut hi = lo;
            for &i in &v[1..] {
                let p = positions[i as usize];
                for axis in 0..3 {
                    if p[axis] < lo[axis] {
                        lo[axis] = p[axis];
                    }
                    if p[axis] > hi[axis] {
                        hi[axis] = p[axis];
                    }
                }
            }
            tri.push([lo, hi]);
        }

        let mut block = Vec::with_capacity(tri.len().div_ceil(BLOCK));
        for chunk in tri.chunks(BLOCK) {
            let mut lo = chunk[0][0];
            let mut hi = chunk[0][1];
            for b in &chunk[1..] {
                for axis in 0..3 {
                    if b[0][axis] < lo[axis] {
                        lo[axis] = b[0][axis];
                    }
                    if b[1][axis] > hi[axis] {
                        hi[axis] = b[1][axis];
                    }
                }
            }
            block.push([lo, hi]);
        }

        Self { tri, block }
    }
}

/// Sample a signed distance field from a closed triangle mesh.
///
/// `positions` and `indices` are the mesh; `shape`, `origin` and `cell_size`
/// describe the grid to sample onto. The result is one value per sample, in the
/// same x-fastest order everything else in this crate uses.
///
/// # The mesh must be closed and consistently oriented
///
/// Not a soft requirement. The sign is `sign(n_α · r)`, and an open mesh has
/// boundary edges whose pseudonormal is a single face normal — which is a
/// perfectly good vector that answers a question with no answer, because there
/// is no inside. A mesh with an inconsistently oriented triangle produces a
/// field that is wrong in a region and right everywhere else.
///
/// This is why S-006 is *"the right tool for geometry isomesh produced itself"*:
/// that geometry already passes this crate's manifold and orientation gates. For
/// imported or damaged input the generalized winding number is the tool, which
/// is S-007.
///
/// **Nothing here checks the mesh**, deliberately: [`validate`](crate::validate)
/// already does, better, and duplicating a weaker check here would be a second
/// path to the same answer.
///
/// # Errors
///
/// [`Error::GridTooSmall`](crate::Error::GridTooSmall) for a grid under 2×2×2,
/// [`Error::ShapeOverflow`](crate::Error::ShapeOverflow) if `indices` is not a
/// multiple of three or names a vertex that does not exist.
pub fn signed_distance_from_mesh<R: Real>(
    positions: &[[R; 3]],
    indices: &[u32],
    shape: &impl Shape3,
    origin: [R; 3],
    cell_size: R,
) -> crate::Result<Vec<R>> {
    let size = shape.size();
    if size[0] < 2 || size[1] < 2 || size[2] < 2 {
        return Err(crate::Error::GridTooSmall { size });
    }
    if !indices.len().is_multiple_of(3) {
        return Err(crate::Error::ShapeOverflow {
            size,
            product: indices.len() as u64,
        });
    }
    if let Some(&bad) = indices.iter().find(|&&i| i as usize >= positions.len()) {
        return Err(crate::Error::ShapeOverflow {
            size,
            product: u64::from(bad),
        });
    }

    let field = MeshField::new(positions, indices)?;
    let count = shape.element_count();
    let (nx, ny, nz) = (size[0], size[1], size[2]);
    let mut out = Vec::with_capacity(count);

    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                out.push(field.sample([
                    origin[0] + cell_size * R::from_f64(f64::from(x)),
                    origin[1] + cell_size * R::from_f64(f64::from(y)),
                    origin[2] + cell_size * R::from_f64(f64::from(z)),
                ]));
            }
        }
    }

    Ok(out)
}

/// A mesh's signed distance field, evaluated on demand.
///
/// Ticket: S-008. [`signed_distance_from_mesh`] answers *"sample this mesh onto
/// that grid"*. This answers *"what is the distance **here**"*, which is what a
/// genuinely sparse consumer wants — a point probe, a collision query, an
/// empty-cell rejection by sphere tracing.
///
/// **Not Manifold Dual Contouring**, which S-008 named as the motivating
/// consumer and which is **not one** (D-011). Its `extract` calls
/// `self.sample(sdf, shape, origin, cell_size)` at `dual.rs:257`, looping every
/// one of the N³ grid points into a buffer before anything else runs, so it
/// reads a grid like every other extractor here and the dense path is the right
/// price for it. The claim came from a summary of the paper rather than from
/// this codebase; whether an on-demand field beats the batch one for a *truly*
/// sparse consumer is open, and S-009 owns it with a crossover to pre-register.
///
/// The acceleration structures — the angle-weighted pseudonormals and the
/// blocked bounding boxes — are built once by [`new`](Self::new) and reused by
/// every [`sample`](Sdf::sample). The mesh itself is borrowed, never copied.
///
/// # There is one implementation of the query, and this is it
///
/// [`signed_distance_from_mesh`] is this type in a loop, literally, since S-008.
/// Two implementations of *"where is the nearest triangle and which side of it
/// am I on"* would be free to disagree at the eighth digit and nothing would
/// notice; `the_grid_path_is_this_field_in_a_loop` asserts they are bit-equal.
/// What batching still buys is building the pseudonormals and the block boxes
/// once for a whole grid rather than once per caller.
///
/// # The same requirement as the batch path, for the same reason
///
/// **The mesh must be closed and consistently oriented.** The sign is
/// `sign(n_α · r)`, and an open mesh has boundary edges whose pseudonormal
/// answers a question that has no answer, because there is no inside. See
/// [`signed_distance_from_mesh`] for the full statement, and
/// [`winding`](super::winding) for the tool that handles imported or damaged
/// input instead.
///
/// # Not available for the winding backend, and that is measured rather than
/// unfinished
///
/// S-008 was scoped for both backends behind one type. The generalized winding
/// number does not fit: [`winding_numbers`](super::winding::winding_numbers)
/// casts **one ray per grid row** and shares it across every sample in that row,
/// which is Martens & Bessmeltsev's *"to compute voxelizations of resolution
/// N³, we only need to shoot N² rays."* A per-point query cannot share a cast
/// with points it has never seen, so an on-demand winding field would cast N³
/// rays for the same grid — a factor of N, not a constant. The batch function
/// keeps the row amortisation and there is no on-demand twin of it.
#[derive(Debug)]
pub struct MeshField<'a, R: Real> {
    positions: &'a [[R; 3]],
    indices: &'a [u32],
    normals: Pseudonormals<R>,
    bounds: Bounds<R>,
}

impl<'a, R: Real> MeshField<'a, R> {
    /// Build the query structures over `positions` and `indices`.
    ///
    /// A trailing partial triangle is dropped rather than rejected — the same
    /// convention [`validate`](crate::validate) counts as `trailing_indices` and
    /// [`collider::triangle_indices`](crate::collider::triangle_indices) applies
    /// with `chunks_exact`. One answer to a ragged index buffer, in all three
    /// places.
    ///
    /// **Nothing here checks that the mesh is closed**, deliberately:
    /// [`validate`](crate::validate) already does it better, and a weaker copy
    /// of that check here would be a second path to the same verdict.
    ///
    /// # Errors
    ///
    /// [`Error::IndexOutOfRange`](crate::Error::IndexOutOfRange) if an index
    /// names a vertex the buffer does not have.
    pub fn new(positions: &'a [[R; 3]], indices: &'a [u32]) -> crate::Result<Self> {
        if let Some((at, &index)) = indices
            .iter()
            .enumerate()
            .find(|&(_, &i)| i as usize >= positions.len())
        {
            return Err(crate::Error::IndexOutOfRange {
                at: at as u64,
                index,
                vertices: positions.len() as u64,
            });
        }

        let whole = indices.len() - indices.len() % 3;
        let indices = &indices[..whole];
        Ok(Self {
            positions,
            indices,
            normals: Pseudonormals::build(positions, indices),
            bounds: Bounds::build(positions, indices),
        })
    }

    /// Triangles the field is built over, after any ragged tail is dropped.
    #[must_use]
    pub const fn triangles(&self) -> usize {
        self.indices.len() / 3
    }
}

impl<R: Real> Sdf for MeshField<'_, R> {
    type Scalar = R;

    fn sample(&self, p: [R; 3]) -> R {
        let triangles = self.triangles();

        // Squared throughout, rooted once at the end. The comparison that drives
        // the reject is monotone in the square, so the square root is pure cost
        // inside the loop.
        let mut best_sq = super::far::<R>();
        let mut best_sign = R::ONE;

        for (b, bbox) in self.bounds.block.iter().enumerate() {
            if box_distance_sq(p, bbox) >= best_sq {
                continue;
            }
            let lo = b * BLOCK;
            let hi = (lo + BLOCK).min(triangles);
            for t in lo..hi {
                if box_distance_sq(p, &self.bounds.tri[t]) >= best_sq {
                    continue;
                }
                let tri = &self.indices[t * 3..t * 3 + 3];
                let (c, feature) = closest_on_triangle(
                    p,
                    self.positions[tri[0] as usize],
                    self.positions[tri[1] as usize],
                    self.positions[tri[2] as usize],
                );
                let r = sub(p, c);
                let d_sq = dot(r, r);
                if d_sq < best_sq {
                    best_sq = d_sq;
                    // Theorem 1: the sign is that of `n_α · r`, with `n_α` the
                    // pseudonormal of the *feature* the closest point landed on.
                    let n = self.normals.at(self.indices, t, feature);
                    best_sign = if dot(n, r) < R::ZERO { -R::ONE } else { R::ONE };
                }
            }
        }

        best_sq.sqrt() * best_sign
    }
}
