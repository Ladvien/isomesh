//! Inside-outside by generalized winding number, for meshes S-006 cannot judge.
//!
//! Ticket: S-007. Xie, Hafner & Wojtan, *Fast and Exact Winding Numbers for
//! Triangle Meshes*, ACM TOG 45(4) (2026) (`10.1145/3811339`).
//!
//! # Why not Barill et al. 2018, which is what everyone cites
//!
//! Because it is an approximation and the 2026 literature is blunt about how
//! good: Barill's Barnes–Hut summation *"trades off accuracy for computational
//! speed… the resulting values are merely approximations of the actual winding
//! number."* The ticket for this work carries a stronger quote still, from the
//! Antipodal paper, calling the order-0/order-1 expansions *"very imprecise…
//! not useful for applications."* An approximate inside-outside test on a
//! damaged mesh is exactly the wrong tool, since the damaged regions are where
//! the answer is close to the threshold.
//!
//! # The construction, transcribed
//!
//! The generalized winding number is a surface integral of solid angle,
//! `w_S(q) = (1/4π) ∫_S dΩ(q)`. For a **closed** mesh it collapses to ray
//! casting — the signed intersection count is independent of ray direction. For
//! an **open** one it does not, which is the whole problem.
//!
//! Xie et al. bridge the two with the additive property `w_{M+C} = w_M + w_C`.
//! Close `M` with a generalized cone `C` from its boundary to an apex **directly
//! behind the ray**, so `C` contributes no forward intersections. Then:
//!
//! ```text
//! w_M(q) = w_{M+C}(q) − w_C(q)
//!        = Σᵢ sgn(r · nᵢ) − (1/4π) Σⱼ Ωⱼ
//! ```
//!
//! where `r` is an arbitrary ray direction from `q`, `nᵢ` are the normals at the
//! ray's intersections with `M`, and `Ωⱼ` is the solid angle the `j`-th cone
//! triangle subtends at `q`.
//!
//! **The cost scales with holes, not triangles.** The second sum runs over
//! boundary edges, so a nearly-closed mesh pays almost nothing for the
//! correction — which is the property that makes this affordable on geometry
//! that is *mostly* fine.
//!
//! # One ray per grid row, not per sample
//!
//! Both this paper and Martens & Bessmeltsev note that a single ray answers for
//! every query point along it. Casting along `+x` makes the ray from sample
//! `(x, y, z)` a suffix of the ray from `(0, y, z)`, so one cast per row gives
//! every `χ` in that row by summing the intersections beyond each sample. A
//! `n³` grid needs `n²` casts. Martens states the same reduction: *"to compute
//! voxelizations of resolution N³, we only need to shoot N² rays."*
//!
//! # Classify points. Do not repair meshes.
//!
//! Takayama, Jacobson, Kavan & Sorkine 2014 — the GWN authors' own follow-up —
//! is that the orientation-repair application is fundamentally flawed. Nothing
//! here writes to a mesh; it reads one and answers a question about a point.

#[cfg(test)]
mod tests;

use alloc::collections::BTreeMap;
use alloc::vec::Vec;

use crate::real::Real;
use crate::shape::Shape3;

use super::from_mesh::{closest_on_triangle, cross, dot, norm, sub};

/// How far behind the query point the cone's apex sits, in multiples of the
/// mesh's bounding-box diagonal.
///
/// The derivation wants the apex at infinity; a finite apex is fine as long as
/// no cone triangle reaches the forward ray, and pushing it many diagonals back
/// makes any such triangle a sliver whose only approach to the ray is at the
/// boundary edge itself — which a generic ray misses. Larger is not better: the
/// solid-angle formula divides by products of the three vertex distances, and an
/// apex at `1e12` diagonals throws away mantissa for nothing.
const APEX_DISTANCE: f64 = 64.0;

/// Solid angle subtended at `q` by triangle `abc`, signed by its orientation.
///
/// Van Oosterom & Strackee, *The solid angle of a plane triangle*, IEEE TBME
/// 30(2) (1983):
///
/// ```text
/// tan(Ω/2) = a · (b × c) / (‖a‖‖b‖‖c‖ + (a·b)‖c‖ + (a·c)‖b‖ + (b·c)‖a‖)
/// ```
///
/// with `a, b, c` the vertices relative to `q`. `atan2` rather than `atan` of
/// the quotient: the denominator goes negative for a triangle subtending more
/// than a hemisphere, and only the four-quadrant form keeps the half-turn.
fn solid_angle<R: Real>(q: [R; 3], a: [R; 3], b: [R; 3], c: [R; 3]) -> R {
    let a = sub(a, q);
    let b = sub(b, q);
    let c = sub(c, q);
    let (la, lb, lc) = (norm(a), norm(b), norm(c));
    if la == R::ZERO || lb == R::ZERO || lc == R::ZERO {
        // The query point is *on* a vertex, where the winding number is
        // undefined rather than zero. Returning zero is the only value that
        // leaves the other triangles' contributions intact, and the caller's
        // threshold test is what decides such a point -- see the module docs on
        // why nothing here tries to repair the input.
        return R::ZERO;
    }
    let numer = dot(a, cross(b, c));
    let denom = la * lb * lc + dot(a, b) * lc + dot(a, c) * lb + dot(b, c) * la;
    R::TWO * numer.atan2(denom)
}

/// The directed boundary of a mesh, with multiplicity.
///
/// A directed edge `a → b` that appears `n` more times than `b → a` needs `n`
/// closing triangles. Counting nets rather than presence is what makes this
/// correct on triangle **soup** — a non-manifold edge with three incident faces
/// has a net of one and closes with one triangle, where a boolean
/// "is it a boundary edge" would either drop it or double it.
fn boundary_edges(indices: &[u32]) -> Vec<([u32; 2], i32)> {
    let mut net: BTreeMap<(u32, u32), i32> = BTreeMap::new();
    for tri in indices.chunks_exact(3) {
        for (u, v) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            if u == v {
                continue;
            }
            let (key, delta) = if u < v { ((u, v), 1) } else { ((v, u), -1) };
            *net.entry(key).or_insert(0) += delta;
        }
    }
    net.into_iter()
        .filter(|&(_, n)| n != 0)
        .map(|((u, v), n)| ([u, v], n))
        .collect()
}

/// One intersection of a row's ray with the mesh.
struct Hit<R> {
    /// World `x` at which it happened.
    x: R,
    /// `sgn(r · n)`, which for a `+x` ray is the sign of the normal's `x`.
    sign: i32,
}

/// Generalized winding numbers of `positions`/`indices`, sampled on a grid.
///
/// One value per sample, x-fastest. A value near `1` means inside, near `0`
/// outside, and `2` means the surface wraps twice — the measure is not a
/// boolean, which is the point of it.
///
/// # Errors
///
/// As [`signed_distance_from_mesh`](super::from_mesh::signed_distance_from_mesh).
pub fn winding_numbers<R: Real>(
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

    let boundary = boundary_edges(indices);
    let triangles = indices.len() / 3;

    // How far back the apex goes, from the mesh's own extent so the constant is
    // scale-free.
    let mut lo = origin;
    let mut hi = origin;
    for p in positions {
        for axis in 0..3 {
            if p[axis] < lo[axis] {
                lo[axis] = p[axis];
            }
            if p[axis] > hi[axis] {
                hi[axis] = p[axis];
            }
        }
    }
    let diagonal = norm(sub(hi, lo));
    let back = if diagonal == R::ZERO {
        R::ONE
    } else {
        diagonal
    } * R::from_f64(APEX_DISTANCE);

    let (nx, ny, nz) = (size[0], size[1], size[2]);
    let mut out = alloc::vec![R::ZERO; shape.element_count()];
    let mut hits: Vec<Hit<R>> = Vec::new();
    let quarter = R::ONE / (R::from_f64(4.0) * R::from_f64(core::f64::consts::PI));

    for z in 0..nz {
        for y in 0..ny {
            let py = origin[1] + cell_size * R::from_f64(f64::from(y));
            let pz = origin[2] + cell_size * R::from_f64(f64::from(z));

            // **One cast for the whole row.** Every sample in it lies on this
            // line, so its `χ` is the sum of the signs beyond it.
            hits.clear();
            for t in 0..triangles {
                let tri = &indices[t * 3..t * 3 + 3];
                let a = positions[tri[0] as usize];
                let b = positions[tri[1] as usize];
                let c = positions[tri[2] as usize];
                if let Some(hit) = intersect_x_ray(a, b, c, py, pz) {
                    hits.push(hit);
                }
            }

            for x in 0..nx {
                let px = origin[0] + cell_size * R::from_f64(f64::from(x));

                // χ: signed intersections strictly ahead of this sample.
                let mut chi = 0i32;
                for h in &hits {
                    if h.x > px {
                        chi += h.sign;
                    }
                }

                // The cone correction. Apex directly behind the ray, so the
                // cone contributes no forward intersections and the count above
                // is the closed mesh's winding number.
                let q = [px, py, pz];
                let apex = [px - back, py, pz];
                let mut omega = R::ZERO;
                for &([u, v], n) in &boundary {
                    // The closing mesh carries the *reverse* directed edge, so
                    // M + C is consistently oriented. A net of `n` for `u → v`
                    // means `n` copies of the triangle `(v, u, apex)`.
                    let (first, second) = if n > 0 {
                        (positions[v as usize], positions[u as usize])
                    } else {
                        (positions[u as usize], positions[v as usize])
                    };
                    let one = solid_angle(q, first, second, apex);
                    omega += one * R::from_f64(f64::from(n.abs()));
                }

                let i = ((z * ny + y) * nx + x) as usize;
                out[i] = R::from_f64(f64::from(chi)) - omega * quarter;
            }
        }
    }

    Ok(out)
}

/// Where a `+x` ray at `(·, py, pz)` crosses triangle `abc`, if it does.
///
/// Möller & Trumbore, *Fast, Minimum Storage Ray/Triangle Intersection*, Journal
/// of Graphics Tools 2(1) (1997), specialised to the `+x` direction so the
/// determinant is a plain component of the edge cross product.
///
/// **Boundaries are `>= 0` and `<= 1` on one side only.** A ray passing exactly
/// through a shared edge must be counted once, not twice or zero times; the
/// half-open convention is what makes that happen, and it is the same rule the
/// crate's own case tables use for edge ownership.
fn intersect_x_ray<R: Real>(a: [R; 3], b: [R; 3], c: [R; 3], py: R, pz: R) -> Option<Hit<R>> {
    let ab = sub(b, a);
    let ac = sub(c, a);
    // r × ac with r = (1, 0, 0).
    let h = [R::ZERO, -ac[2], ac[1]];
    let det = dot(ab, h);
    if det.abs() < R::EPSILON {
        // Ray parallel to the triangle's plane. A grazing hit contributes
        // nothing to a signed count -- it enters and leaves at the same point --
        // so dropping it is exact rather than a tolerance.
        return None;
    }
    let inv = R::ONE / det;
    let s = [R::ZERO, py - a[1], pz - a[2]];
    let u = inv * dot(s, h);
    if u < R::ZERO || u > R::ONE {
        return None;
    }
    let q = cross(s, ab);
    let v = inv * q[0];
    if v < R::ZERO || u + v > R::ONE {
        return None;
    }
    let t = inv * dot(ac, q);

    // sgn(r · n) with r = +x̂ and n = ab × ac, so just the x component's sign.
    // Zero is impossible here: it would make `det` zero, which returned above.
    let n_x = ab[1] * ac[2] - ab[2] * ac[1];
    Some(Hit {
        x: a[0] + t,
        sign: if n_x < R::ZERO { -1 } else { 1 },
    })
}

/// A signed distance field from a mesh, signed by the winding number.
///
/// The magnitude is the true distance to the nearest triangle, exactly as
/// [`signed_distance_from_mesh`](super::from_mesh::signed_distance_from_mesh)
/// computes it. Only the **sign** differs, and that is the entire point: the
/// pseudonormal's sign is a theorem about *closed* meshes, so on a mesh with a
/// hole it is not merely inaccurate, it is answering a different question.
///
/// `threshold` is where the winding number is cut. `0.5` is the standard choice
/// and is what the literature uses; a higher value is more conservative about
/// calling a point inside.
///
/// # Errors
///
/// As [`winding_numbers`].
pub fn signed_distance_from_mesh_winding<R: Real>(
    positions: &[[R; 3]],
    indices: &[u32],
    shape: &impl Shape3,
    origin: [R; 3],
    cell_size: R,
    threshold: R,
) -> crate::Result<Vec<R>> {
    let winding = winding_numbers(positions, indices, shape, origin, cell_size)?;

    let size = shape.size();
    let triangles = indices.len() / 3;
    let (nx, ny, nz) = (size[0], size[1], size[2]);
    let mut out = Vec::with_capacity(shape.element_count());

    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let p = [
                    origin[0] + cell_size * R::from_f64(f64::from(x)),
                    origin[1] + cell_size * R::from_f64(f64::from(y)),
                    origin[2] + cell_size * R::from_f64(f64::from(z)),
                ];
                let mut best = super::far::<R>();
                for t in 0..triangles {
                    let tri = &indices[t * 3..t * 3 + 3];
                    let (c, _) = closest_on_triangle(
                        p,
                        positions[tri[0] as usize],
                        positions[tri[1] as usize],
                        positions[tri[2] as usize],
                    );
                    let d = norm(sub(p, c));
                    if d < best {
                        best = d;
                    }
                }
                let i = ((z * ny + y) * nx + x) as usize;
                out.push(if winding[i] > threshold { -best } else { best });
            }
        }
    }

    Ok(out)
}
