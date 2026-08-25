//! Where a vertex normal comes from, made a decision rather than an assumption.
//!
//! Every extractor here has always produced one normal per vertex, from the
//! field's gradient at the vertex's final position. That is one of three
//! defensible answers and it was never selectable, which meant the other two
//! could not be compared against it.
//!
//! | Strategy | Costs | Right when |
//! |---|---|---|
//! | [`AnalyticGradient`](NormalStrategy::AnalyticGradient) | one [`Sdf::gradient`] per vertex | the field implements `gradient`, so the answer is exact |
//! | [`CentralDifference`](NormalStrategy::CentralDifference) | six [`Sdf::sample`] per vertex | the field is *sampled* — a voxel buffer has no analytic gradient to ask for, and this is what a game actually has |
//! | [`AreaWeightedFaces`](NormalStrategy::AreaWeightedFaces) | no field evaluations at all | the mesh is all you have, or you want the normal to agree with the geometry rather than with the field |
//!
//! # Why a post-pass and not a setting on six extractors
//!
//! One path. A `set_normal_strategy` on each of the six extractors would be six
//! copies of the same branch, six places for the default to drift, and it would
//! still not serve a caller who wants to re-derive normals on a mesh that has
//! already been welded or merged. [`recompute`] takes the finished
//! [`MeshBuffer`] and re-derives every normal in it, which works for all six and
//! for anything downstream of them.
//!
//! It also makes the no-op assertable:
//! `recomputing_with_the_analytic_gradient_reproduces_extraction` requires
//! [`AnalyticGradient`](NormalStrategy::AnalyticGradient) to give back exactly
//! the normals extraction already produced, bit for bit. If that ever fails, the
//! post-pass and the extractors have diverged.
//!
//! # The one thing this does not do
//!
//! Smooth. There is no averaging over a neighbourhood and no crease angle. A
//! crease threshold needs *split* vertices to be any use — a cube corner wants
//! three normals and a `MeshBuffer` vertex carries one — and splitting is
//! [`crate::greedy_quads`]' business, not this module's.

use crate::{Error, MeshBuffer, Real, Result, Sdf, vec3};

/// Where each vertex normal comes from.
///
/// Non-exhaustive in spirit but not in fact: three is the whole set that A-012
/// names, and adding a fourth should be a decision with a measurement attached
/// rather than a convenience.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NormalStrategy<R: Real> {
    /// Ask the field, via [`Sdf::gradient`].
    ///
    /// **The default every extractor already uses.** Exact wherever the field
    /// overrides `gradient` — all seven reference fields do — and silently the
    /// central-difference default where it does not, because that is what the
    /// trait's default implementation is.
    AnalyticGradient,

    /// Difference the field, ignoring any analytic gradient it offers.
    ///
    /// Six [`Sdf::sample`] calls with an isotropic step, which is the same
    /// stencil [`Sdf::gradient`]'s default uses — *the same on all three axes, so
    /// the returned direction is unbiased.*
    ///
    /// `step` is absolute. Passing the grid's cell size is the case worth
    /// measuring: it is what a voxel buffer can actually offer, since a sampled
    /// field has nothing finer than its own spacing to difference over.
    CentralDifference {
        /// The differencing step. Must be finite and positive.
        step: R,
    },

    /// Accumulate the incident triangles' cross products, then normalise.
    ///
    /// **Area weighting comes free from not normalising the face normal**: for a
    /// triangle `(a, b, c)` the cross product `(b − a) × (c − a)` has magnitude
    /// twice the area, so summing raw cross products *is* the area-weighted sum
    /// and dividing each by its own length would be the thing that throws the
    /// weighting away.
    ///
    /// Uses no field evaluations, so it is the only strategy available to a
    /// caller holding a mesh and nothing else.
    AreaWeightedFaces,
}

/// Re-derive every normal in `mesh`.
///
/// Overwrites `mesh.normals` in place; positions and indices are untouched. The
/// buffer keeps its allocation, per rule 6.
///
/// # Errors
///
/// - [`Error::InvalidCellSize`] if a
///   [`CentralDifference`](NormalStrategy::CentralDifference) step is not finite
///   and positive. A zero step makes every difference zero and every normal
///   undefined, which is worth refusing at the door rather than discovering as a
///   black mesh.
/// - [`Error::DegenerateNormal`] if a vertex ends up with nothing to normalise: a
///   zero gradient under either field strategy, or, under
///   [`AreaWeightedFaces`](NormalStrategy::AreaWeightedFaces), a vertex whose
///   incident triangles have no area or which no triangle references at all.
///   **Reported rather than papered over** — substituting `[0, 0, 1]` there
///   would give a mesh that renders and is wrong in a way nothing downstream
///   could attribute.
pub fn recompute<R, S>(mesh: &mut MeshBuffer<R>, sdf: &S, strategy: NormalStrategy<R>) -> Result<()>
where
    R: Real,
    S: Sdf<Scalar = R>,
{
    match strategy {
        NormalStrategy::AnalyticGradient => from_field(mesh, |p| sdf.gradient(p)),
        NormalStrategy::CentralDifference { step } => {
            if !step.is_finite() || step <= R::ZERO {
                return Err(Error::InvalidCellSize {
                    value: f64::from(step.as_f32()),
                });
            }
            from_field(mesh, |p| central_difference(sdf, p, step))
        }
        NormalStrategy::AreaWeightedFaces => area_weighted(mesh),
    }
}

/// Central differences at an explicit step, ignoring the field's own gradient.
///
/// Split out because it is the whole of what
/// [`CentralDifference`](NormalStrategy::CentralDifference) means, and a caller
/// comparing an analytic gradient against a numerical one wants to call both
/// directly rather than through a mesh.
///
/// The step is the same on all three axes. A per-axis step would bias every
/// direction it returned, which is the reason [`Sdf::gradient`]'s default says so
/// too.
#[must_use]
pub fn central_difference<R, S>(sdf: &S, p: [R; 3], step: R) -> [R; 3]
where
    R: Real,
    S: Sdf<Scalar = R>,
{
    let inv = (R::TWO * step).recip();
    [
        (sdf.sample([p[0] + step, p[1], p[2]]) - sdf.sample([p[0] - step, p[1], p[2]])) * inv,
        (sdf.sample([p[0], p[1] + step, p[2]]) - sdf.sample([p[0], p[1] - step, p[2]])) * inv,
        (sdf.sample([p[0], p[1], p[2] + step]) - sdf.sample([p[0], p[1], p[2] - step])) * inv,
    ]
}

fn from_field<R, F>(mesh: &mut MeshBuffer<R>, mut gradient: F) -> Result<()>
where
    R: Real,
    F: FnMut([R; 3]) -> [R; 3],
{
    for (vertex, position) in mesh.positions.iter().enumerate() {
        let g = gradient(*position);
        let length = vec3::length(g);
        // `length > ZERO` is written positively and then negated by the `else`
        // rather than as `!(length > ZERO)`: a NaN length is not greater than
        // zero, so it falls to the error arm, which is what we want.
        if length > R::ZERO && length.is_finite() {
            mesh.normals[vertex] = vec3::scale(g, length.recip());
        } else {
            return Err(Error::DegenerateNormal {
                vertex: vertex as u64,
            });
        }
    }
    Ok(())
}

fn area_weighted<R: Real>(mesh: &mut MeshBuffer<R>) -> Result<()> {
    mesh.normals.fill([R::ZERO; 3]);

    for tri in mesh.indices.as_chunks::<3>().0 {
        let (i, j, k) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
        // A triangle naming a vertex that does not exist is the validator's
        // finding, not this pass's; skipping it here keeps the two from
        // disagreeing about the same mesh.
        if i >= mesh.positions.len() || j >= mesh.positions.len() || k >= mesh.positions.len() {
            continue;
        }
        let a = mesh.positions[i];
        let cross = vec3::cross(
            vec3::sub(mesh.positions[j], a),
            vec3::sub(mesh.positions[k], a),
        );
        for vertex in [i, j, k] {
            for (axis, slot) in mesh.normals[vertex].iter_mut().enumerate() {
                *slot += cross[axis];
            }
        }
    }

    for (vertex, normal) in mesh.normals.iter_mut().enumerate() {
        let length = vec3::length(*normal);
        if length > R::ZERO && length.is_finite() {
            *normal = vec3::scale(*normal, length.recip());
        } else {
            return Err(Error::DegenerateNormal {
                vertex: vertex as u64,
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
