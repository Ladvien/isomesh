//! Volume, centre of mass and the inertia tensor of the solid the triangles
//! bound — computed on the surface, with no volume mesh and no round trip
//! through a physics crate.
//!
//! Ticket: R-083, pre-registered as `P-83`.
//!
//! # The construction
//!
//! Hartmann & Ewougsi Tekeu, *Gauss divergence theorem for the calculation of
//! the mass and area moment of inertia tensors*, Acta Mechanica 236 (2025),
//! [`10.1007/s00707-025-04419-1`](https://doi.org/10.1007/s00707-025-04419-1),
//! open access. Their equation numbers are cited below; this module is their
//! §5.1 with nothing added.
//!
//! The mass properties of a solid `Ω` are volume integrals, and the obvious way
//! to get them is to fill `Ω` with tetrahedra or voxels and sum. This crate has
//! already produced `∂Ω`, so filling it again is paying twice for geometry it is
//! holding. Substituting the position vector into the divergence theorem removes
//! the volume entirely. With `div x = 3` (their Eq. 5–6),
//!
//! ```text
//! V = ⅓ ∮ x·n dA
//! ```
//!
//! with `div(x ⊗ x) = 4x` (Eq. 8, 12),
//!
//! ```text
//! r = 1/(4V) ∮ (x·n) x dA
//! ```
//!
//! and — the part they call novel — with `div((x·x)x) = 5(x·x)` and
//! `grad((x·x)x) = 2 x⊗x + (x·x)I` (Eq. 17–21), the whole inertia tensor about
//! the origin:
//!
//! ```text
//! Θ₀ = 3/10 (∮ (x·x)(x·n) dA) I − ½ ∮ (x·x) x ⊗ n dA
//! ```
//!
//! # The trap, which is the reason [`MassProperties::asymmetry`] is public
//!
//! `Θ₀` is symmetric by definition, so the second surface integral must satisfy
//! their Eq. (22), `∮(x·x) x⊗n dA = ∮(x·x) n⊗x dA`. The paper's own words:
//! *"which holds in the exact case but might be violated by discretization
//! schemes"*, and their §6.1 measures it violated on a triangulated cylinder —
//! `Θˣᶻ ≠ Θᶻˣ` — with the remedy `Θ = ½(Θ + Θᵀ)` (Eq. 66).
//!
//! The mechanism is visible in the sum: each triangle contributes the rank-one
//! `G ⊗ N`, which is antisymmetric-in-part and only cancels against the rest of
//! the surface. So the asymmetry is a *global cancellation residual*, and it
//! reports on the surface rather than on the arithmetic — see the field's docs.
//!
//! # One quadrature kernel, not ten
//!
//! On a triangle `(X₁, X₂, X₃)` with `N = (X₂ − X₁) × (X₃ − X₁)`, the outward
//! `n dA` is the constant `N` times the parametric area element, so each
//! integral above is a polynomial in the linear shape functions
//! `n₁ = 1 − ξ − η`, `n₂ = ξ`, `n₃ = η` over the unit triangle (their Eq. 44),
//! and `∫ n₁^a n₂^b n₃^c dξ dη = a!b!c!/(a+b+c+2)!`. Their Table 1 is the cubic
//! case of that, and it collects — writing `S = ΣXₐ` and `Q = Σ Xₐ·Xₐ` — into
//!
//! ```text
//! 120·∫ (x·x) x dξdη = (S·S)S + Q·S + 2·Σₐ(S·Xₐ)Xₐ + 2·Σₐ(Xₐ·Xₐ)Xₐ
//! ```
//!
//! which is this module's one quadrature kernel. `‖N‖` cancels between `dA` and
//! the unit normal, so
//! **no square root is evaluated anywhere in this module** and the triangle's
//! area is never formed.
//!
//! # Cost
//!
//! One pass, no allocation, no scratch buffer — the accumulator is fourteen
//! scalars. There is nothing for a caller to provide and nothing to reuse, which
//! is why this module has no `&mut` output parameter while the extractors all
//! do.

use crate::vec3::{cross, dot, sub};
use crate::{Error, Real, Result};

#[cfg(test)]
mod tests;

/// What a closed triangle mesh weighs, where it balances, and how it spins.
///
/// Unit density throughout: multiply [`volume`](Self::volume) by a density to
/// get mass, and either tensor by the same density to get the physical inertia.
/// Density is not a parameter because a uniform one is a scalar multiple of
/// every field here, and taking it would let a caller believe the module
/// supports a varying one.
///
/// # Units
///
/// With positions in metres: `volume` is m³, `center_of_mass` is m, both tensors
/// are m⁵ — inertia per unit density — and [`asymmetry`](Self::asymmetry) is m⁵
/// so that it is directly comparable against the tensor entries beside it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassProperties<R: Real> {
    /// Enclosed volume. Always finite and strictly positive — a surface that
    /// bounds anything else is an [`Error::MassPropertiesUndefined`] instead.
    pub volume: R,

    /// Centroid of the enclosed solid, in the same frame as the positions.
    pub center_of_mass: [R; 3],

    /// Inertia tensor about [`center_of_mass`](Self::center_of_mass),
    /// row-major, at unit density. `Θ_M` in the paper.
    ///
    /// `inertia[0][0] = ∫(y² + z²) dV` and `inertia[0][1] = −∫xy dV`, i.e. the
    /// physics convention where the off-diagonals are the *negated* products of
    /// inertia and `I·ω` is the angular momentum. Exactly symmetric: each
    /// off-diagonal entry is computed once and mirrored.
    pub inertia: [[R; 3]; 3],

    /// Inertia tensor about the coordinate origin, row-major, at unit density.
    /// `Θ₀` in the paper.
    ///
    /// Kept alongside [`inertia`](Self::inertia) because both are already
    /// computed — the centred one is this minus the Steiner contribution
    /// `V[(r·r)I − r⊗r]` — and because performing that shift is the single most
    /// common place a caller loses a sign. A reference integrator accumulating
    /// moments about the origin compares against this one directly.
    pub inertia_about_origin: [[R; 3]; 3],

    /// How far the discretised tensor was from symmetric before it was
    /// symmetrised, in tensor units: `maxᵢ<ⱼ |Θᵢⱼ − Θⱼᵢ|`.
    ///
    /// The paper's Eq. (22) says `∮(x·x) x⊗n dA` must equal its own transpose,
    /// and adds that this *"might be violated by discretization schemes"*.
    /// This is the violation, measured rather than assumed away.
    ///
    /// **It is a leak detector, and that is why it is public.** Each triangle
    /// contributes a rank-one `G ⊗ N` whose antisymmetric part is cancelled
    /// only by the rest of a *closed* surface. On a watertight mesh the
    /// cancellation is exact in exact arithmetic, so this sits at the round-off
    /// floor of that cancellation — small, and not zero. On a mesh with a
    /// boundary edge, a T-junction or a flipped triangle there is nothing to
    /// cancel against, and it rises to the scale of the hole. Nothing else this
    /// module returns can tell a caller that, and the caller cannot recompute it
    /// from the symmetrised tensor.
    pub asymmetry: R,
}

/// Mass properties of the solid bounded by `triangles`.
///
/// `triangles` indexes `positions`, and the winding is the crate's:
/// counter-clockwise viewed from outside the solid, so `(X₂ − X₁) × (X₃ − X₁)`
/// points away from the material. The paper makes the same demand — *"a
/// fundamental problem with triangulation is the numbering of the nodes, which
/// results in a normal pointing outward from the body"* (§5.1). A consistently
/// *inward*-wound mesh yields a negative volume and is rejected rather than
/// silently flipped: the sign is the caller's contract, and repairing it here
/// would make an inside-out mesh indistinguishable from a correct one.
///
/// The mesh must be closed. Nothing here can check that cheaply — closedness is
/// an edge-adjacency property and this is a single pass over triangles — so the
/// module reports [`MassProperties::asymmetry`] instead, which rises off the
/// round-off floor exactly when the surface leaks. Use
/// [`validate_indexed`](crate::validate::validate_indexed) when a definite
/// answer is wanted.
///
/// Triangles are consumed in slice order and nothing is sorted, so the result is
/// a deterministic function of the input order. Float addition is not
/// associative; a caller that reorders its triangle list will get a different
/// last bit.
///
/// # Errors
///
/// [`Error::IndexOutOfRange`] if a triangle names a vertex `positions` does not
/// have. [`Error::MassPropertiesUndefined`] if the surface encloses no positive
/// finite volume, or if a moment overflowed to a non-finite value — in either
/// case the centre of mass is a division by something that is not a volume, and
/// the crate reports that rather than returning one.
///
/// # Example
///
/// ```
/// use isomesh::mass::mass_properties;
///
/// // The unit cube [0,1]³, twelve triangles, wound outward.
/// let positions = [
///     [0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
///     [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
/// ];
/// let triangles = [
///     [0u32, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
///     [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
///     [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7],
/// ];
///
/// let props = mass_properties(&positions, &triangles)?;
/// assert!((props.volume - 1.0).abs() < 1e-15);
/// assert!((props.center_of_mass[0] - 0.5).abs() < 1e-15);
/// // A unit cube about its centre: I = 1/6 on the diagonal, 0 off it.
/// assert!((props.inertia[0][0] - 1.0 / 6.0).abs() < 1e-15);
/// assert!(props.inertia[0][1].abs() < 1e-15);
/// # Ok::<(), isomesh::Error>(())
/// ```
pub fn mass_properties<R: Real>(
    positions: &[[R; 3]],
    triangles: &[[u32; 3]],
) -> Result<MassProperties<R>> {
    // Fourteen accumulators. `outer` is Σ G ⊗ N *unsymmetrised* and in the order
    // the triangles arrive, which is what makes `asymmetry` a measurement rather
    // than a restatement of the symmetrisation that follows it.
    let mut volume_sum = R::ZERO;
    let mut centroid_sum = [R::ZERO; 3];
    let mut trace_sum = R::ZERO;
    let mut outer_sum = [[R::ZERO; 3]; 3];

    for (t, triangle) in triangles.iter().enumerate() {
        let mut x = [[R::ZERO; 3]; 3];
        for (k, &index) in triangle.iter().enumerate() {
            let Some(&corner) = positions.get(index as usize) else {
                return Err(Error::IndexOutOfRange {
                    at: (t as u64) * 3 + k as u64,
                    index,
                    vertices: positions.len() as u64,
                });
            };
            x[k] = corner;
        }

        // `n dA` for the whole triangle: constant over it, magnitude twice the
        // area, pointing out of the solid under the crate's winding. The paper's
        // Eq. (46), (47), (49) — with `‖N‖` left in, because it cancels against
        // the area element in every integral below.
        let n = cross(sub(x[1], x[0]), sub(x[2], x[0]));
        let s = [
            x[0][0] + x[1][0] + x[2][0],
            x[0][1] + x[1][1] + x[2][1],
            x[0][2] + x[1][2] + x[2][2],
        ];

        // The three per-corner contractions every integral below is built from.
        // Hoisted, because each of them was otherwise recomputed once per axis:
        // nine dot products where three will do, on the innermost loop of the
        // whole module.
        let mut xn = [R::ZERO; 3];
        for (slot, corner) in xn.iter_mut().zip(&x) {
            *slot = dot(*corner, n);
        }
        let sn = dot(s, n);

        // V = ⅓∮x·n dA, and ∫nₐ dξdη = 1/6, so the triangle gives (S·N)/6.
        volume_sum += sn;

        // r = 1/(4V) ∮(x·n)x dA, and ∫nₐn_b dξdη is 1/12 on the diagonal and
        // 1/24 off it, which collects into (S·N)S + Σₐ(Xₐ·N)Xₐ over 24.
        for axis in 0..3 {
            let mut own = R::ZERO;
            for (weight, corner) in xn.iter().zip(&x) {
                own += *weight * corner[axis];
            }
            centroid_sum[axis] += sn * s[axis] + own;
        }

        // Θ₀ = 3/10 (∮(x·x)(x·n) dA) I − ½ ∮(x·x) x⊗n dA, with the vector
        // integral shared between the two terms exactly as in the paper's
        // Eq. (41): one `G` per triangle, contracted for the trace and outer-
        // multiplied for the tensor.
        let g = cubic_moment(&x, s);
        trace_sum += dot(g, n);
        for i in 0..3 {
            for j in 0..3 {
                outer_sum[i][j] += g[i] * n[j];
            }
        }
    }

    // Leading factors, applied once. One correctly rounded division per number
    // rather than a rounding per triangle.
    //
    //   volume   ⅓ · 1/6            = 1/18
    //   centroid 1/(4V) · 1/24      = 1/96 then /V
    //   trace    3/10 · 1/120       = 1/400
    //   outer    ½ · 1/120          = 1/240
    let volume = volume_sum / R::from_f64(18.0);
    let trace = trace_sum / R::from_f64(400.0);

    // Θ₀ before symmetrisation, exactly as the discretisation produces it.
    let mut raw = [[R::ZERO; 3]; 3];
    let mut largest = R::ZERO;
    for i in 0..3 {
        for j in 0..3 {
            let diagonal = if i == j { trace } else { R::ZERO };
            raw[i][j] = diagonal - outer_sum[i][j] / R::from_f64(240.0);
            largest = largest.max(raw[i][j].abs());
        }
    }

    if !(volume.is_finite() && volume > R::ZERO) || !largest.is_finite() {
        return Err(Error::MassPropertiesUndefined {
            volume: volume.as_f64(),
            largest_moment: largest.as_f64(),
        });
    }

    // The paper's Eq. (66), and the residual it hides.
    let mut asymmetry = R::ZERO;
    let mut inertia_about_origin = [[R::ZERO; 3]; 3];
    for i in 0..3 {
        inertia_about_origin[i][i] = raw[i][i];
        for j in (i + 1)..3 {
            asymmetry = asymmetry.max((raw[i][j] - raw[j][i]).abs());
            let mean = (raw[i][j] + raw[j][i]) * R::HALF;
            inertia_about_origin[i][j] = mean;
            inertia_about_origin[j][i] = mean;
        }
    }

    let scale = volume * R::from_f64(96.0);
    let c = [
        centroid_sum[0] / scale,
        centroid_sum[1] / scale,
        centroid_sum[2] / scale,
    ];

    // Steiner, Eq. (25): Θ_M = Θ₀ − V[(r·r)I − r⊗r], unit density so m = V.
    // Written once per *pair* rather than once per entry, so the result stays
    // exactly symmetric instead of symmetric-to-an-ulp.
    let mut inertia = [[R::ZERO; 3]; 3];
    for i in 0..3 {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        inertia[i][i] = inertia_about_origin[i][i] - volume * (c[j] * c[j] + c[k] * c[k]);
        let off = inertia_about_origin[i][j] + volume * c[i] * c[j];
        inertia[i][j] = off;
        inertia[j][i] = off;
    }

    Ok(MassProperties {
        volume,
        center_of_mass: c,
        inertia,
        inertia_about_origin,
        asymmetry,
    })
}

/// `120 · ∫ (x·x) x dξ dη` over the unit triangle, for `x` interpolated
/// linearly from the three corners.
///
/// The paper's Eq. (48) with its Table 1 folded in: a cubic product of
/// barycentrics integrates to `1/20` when all three indices coincide, `1/60`
/// when exactly two do and `1/120` otherwise, which is
/// `(1 + [a=b] + [b=c] + [a=c] + 2[a=b=c])/120`. Summing that against
/// `(Xₐ·X_b)X_c` gives the five terms below, so the triple sum over 27 shape
/// products never has to be written out.
///
/// `s` is `X₁ + X₂ + X₃`, passed in because the caller already formed it for the
/// volume.
///
/// The `1/120` is deferred to the caller, so that the trace term's `3/10` and
/// the outer term's `½` are each applied in a single division at the end.
///
/// The two per-corner contractions are formed once rather than once per axis —
/// this is the innermost arithmetic of the module and the naive nesting does
/// eighteen dot products where six suffice.
#[inline]
fn cubic_moment<R: Real>(x: &[[R; 3]; 3], s: [R; 3]) -> [R; 3] {
    let mut weight = [R::ZERO; 3];
    let mut q = R::ZERO;
    for (slot, corner) in weight.iter_mut().zip(x) {
        let own = dot(*corner, *corner);
        q += own;
        // The `Σₐ(S·Xₐ)Xₐ` and `Σₐ(Xₐ·Xₐ)Xₐ` terms share their vector factor, so
        // they share one multiply per axis instead of two.
        *slot = dot(s, *corner) + own;
    }
    let k = dot(s, s) + q;
    let mut out = [R::ZERO; 3];
    for axis in 0..3 {
        let mut acc = R::ZERO;
        for (w, corner) in weight.iter().zip(x) {
            acc += *w * corner[axis];
        }
        out[axis] = k * s[axis] + acc * R::TWO;
    }
    out
}
