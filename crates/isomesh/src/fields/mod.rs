//! The eight reference fields.
//!
//! One definition, shared by tests, benchmarks and every example, so that
//! comparisons between algorithms are actually comparisons between algorithms.
//! These are public API rather than `#[cfg(test)]` items precisely because the
//! engine wrapper lives in a separate workspace and consumes this crate
//! normally — a field defined twice is a field that drifts.
//!
//! Every field is negative inside, positive outside, and every one of them
//! overrides [`Sdf::gradient`] with an analytic gradient. The central-difference
//! default is never used by a reference field.
//!
//! # Canonical instances
//!
//! Each field has a `canonical()` constructor. **That constructor is the
//! anti-drift mechanism**: nothing anywhere hard-codes a radius or a half-extent,
//! so changing a parameter changes it everywhere at once. Use
//! [`for_each_reference_field!`](crate::for_each_reference_field) to sweep all
//! eight without dynamic dispatch.
//!
//! # Not all of them are closed, and not all of them are distances
//!
//! `gyroid` is triply periodic, so any finite sampling box cuts it and the
//! result has boundary; `fbm_terrain` is a heightfield and leaves through the
//! sides of its domain by construction. Neither has an Euler characteristic that
//! can be asserted a priori. That is what [`ReferenceField`] is for — a harness
//! asks the field what to expect rather than hard-coding `χ == 2` and
//! discovering the problem later.

pub(crate) mod noise;

use crate::vec3::{length, scale, sub};
use crate::{Real, Sdf};

/// What a field's *value* is, as distinct from where its zero set is.
///
/// Ticket: F-001. Sources: Bálint, Valasek & Gergó, *Operations on signed
/// distance functions*, Acta Cybernetica 24(1) (`10.14232/actacyb.24.1.2019.3`),
/// Corollary 1; and their follow-up *Operations on Signed Distance Function
/// Estimates*, CAD&A 20(6) (`10.14733/cadaps.2023.1154-1174`), for `q`.
///
/// # Why a `bool` was the wrong type
///
/// This replaced `is_exact_distance() -> bool`, and the thing that made it worth
/// replacing was sitting in the source: [`CsgDifference`] returned
///
/// ```text
/// fn is_exact_distance(&self) -> bool {
///     true // away from the seam
/// }
/// ```
///
/// — a comment admitting the invariant is false, on a function returning true.
/// The comment is right and the value is wrong, and no `bool` can hold both.
///
/// # The two numbers, and what they are for
///
/// **`l`, the Lipschitz constant.** Corollary 1 of the 2019 paper is exact about
/// this: *"Every signed distance function is Lipschitz continuous and their
/// smallest Lipschitz constant is 1."* So `l` is what bounds how fast the field
/// can change — a sphere tracer may step by `|f| / l` and be guaranteed not to
/// pass through the surface. `min`/`max` of 1-Lipschitz functions is
/// 1-Lipschitz, so **the Lipschitz bound survives arbitrary CSG**, which is what
/// makes Phase 12's empty-cell rejection correct under unlimited player editing.
///
/// **`q`, the underestimate ratio.** The 2023 paper's precision: `q · d(p) ≤
/// |f(p)|` for the true distance `d`, so the field never overstates distance and
/// understates it by at most a factor `q`. Exactness does *not* survive CSG —
/// the paper is blunt that *"the minimum of two exact signed distance functions…
/// is not an exact SDF, and the error can be arbitrarily large under certain
/// conditions"* — but a `q` does, degrading in a way their Theorem 6 quantifies.
///
/// # What reads this
///
/// [`accuracy`](crate::validate::accuracy) and the shootout, which must not
/// report Hausdorff against a field whose values are not distances: the harness
/// measures geometry against `|sample|`, so a field that understates distance
/// makes an extractor look better than it is.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FieldBound {
    /// `|f(p)|` **is** the distance to the zero set, everywhere in the domain.
    ///
    /// Equivalently `Lipschitz { l: 1 }` with `q = 1`. Only this admits a
    /// Hausdorff measurement against the field's own values.
    Exact,
    /// `|f(x) − f(y)| ≤ l · |x − y|`, but the value is not the distance.
    ///
    /// A sphere tracer may still step by `|f| / l`. `l > 1` is the interesting
    /// case: an implicit function that crosses zero at the right place while
    /// changing faster than distance does.
    Lipschitz {
        /// The smallest constant known to hold. `1` for a true distance.
        l: f64,
    },
    /// A conservative lower bound: `q · d(p) ≤ |f(p)| ≤ d(p)`, with `q ∈ (0, 1]`.
    ///
    /// **The safe direction**, and the one CSG produces. A tracer stepping by
    /// `|f|` only ever under-steps, which costs iterations rather than
    /// correctness.
    Underestimate {
        /// The precision. `1` would be exact; smaller is looser.
        q: f64,
    },
    /// Nothing is claimed beyond the sign.
    ///
    /// The zero set is still the surface — extraction only reads signs and
    /// interpolates — but no step size can be derived and no accuracy figure
    /// against `|sample|` means anything.
    Unbounded,
}

impl FieldBound {
    /// Whether `|f|` may be compared against a true distance.
    ///
    /// The one question the accuracy harness asks. Deliberately not a
    /// `matches!` at each call site: there is one definition of "exact enough to
    /// measure against", and it is here.
    #[must_use]
    pub const fn is_exact(self) -> bool {
        matches!(self, Self::Exact)
    }

    /// The Lipschitz constant a sphere tracer may divide by, if one is known.
    ///
    /// `Exact` and `Underestimate` are both 1-Lipschitz — the underestimate
    /// ratio bounds the *value*, not the rate of change — so both answer `1`.
    #[must_use]
    pub const fn lipschitz(self) -> Option<f64> {
        match self {
            Self::Exact | Self::Underestimate { .. } => Some(1.0),
            Self::Lipschitz { l } => Some(l),
            Self::Unbounded => None,
        }
    }
}

/// What a validity or accuracy harness needs to know about a field before it can
/// decide what a correct extraction looks like.
///
/// Without this, a test suite ends up with a per-field `if` ladder. With it there
/// is one rule: closed fields must produce a closed manifold, open fields a
/// manifold with boundary, and the Euler characteristic is asserted only where it
/// is analytically known.
pub trait ReferenceField: Sdf {
    /// Stable identifier, used as the key in golden-hash fixtures.
    const NAME: &'static str;

    /// The axis-aligned box this field is meant to be sampled over, as
    /// `(min, max)`.
    ///
    /// Chosen so the zero set never touches the wall — for the closed fields,
    /// `sample` has a constant sign on the whole boundary of this box.
    fn domain(&self) -> ([Self::Scalar; 3], [Self::Scalar; 3]);

    /// `true` when the surface inside [`domain`](ReferenceField::domain) is
    /// closed, i.e. bounded by geometry lying strictly inside it.
    ///
    /// **False for [`FbmTerrain`]**, which exits through the sides. A caller must
    /// not require zero boundary edges when this is false.
    fn closed_in_domain(&self) -> bool;

    /// The Euler characteristic a correct extraction must produce, when it is
    /// known analytically.
    ///
    /// `None` means *not derivable a priori* — record the observed value in a
    /// golden fixture rather than asserting an invented one. `None` for
    /// [`CappedGyroid`], whose genus depends on how many tunnels the cap
    /// encloses at that scale, and for [`FbmTerrain`], which is not closed.
    fn expected_euler(&self) -> Option<i64>;

    /// `true` when `|∇f| == 1` almost everywhere, i.e. the field really is a
    /// signed distance.
    ///
    /// **False for [`Gyroid`] and [`FbmTerrain`].** An accuracy harness must not
    /// treat `|sample(v)|` as a distance to the surface when this is false.
    fn bound(&self) -> FieldBound;
}

/// Runs a block once per reference field, with the field bound to a **concrete
/// type**.
///
/// A `Vec<Box<dyn Sdf>>` would be shorter and would put a virtual call on the
/// innermost loop of every benchmark, so the numbers would measure dispatch.
/// This expands instead.
///
/// ```
/// # use isomesh::{for_each_reference_field, fields::ReferenceField, Sdf};
/// let mut n = 0;
/// for_each_reference_field!(f32, |name, field| {
///     assert!(!name.is_empty());
///     let (lo, _hi) = field.domain();
///     let _ = field.sample(lo);
///     n += 1;
/// });
/// assert_eq!(n, 8);
/// ```
///
/// # It looks like a closure and it is not
///
/// The `|name, field|` is syntax, not a closure: the body is **inlined once per
/// field**, because the eight fields are eight different types and no single
/// closure can take all of them. So a `return` in the body returns from the
/// **enclosing function**, not from one iteration — a test that skips fields
/// with `if name != "…" { return; }` exits on `sphere` and silently stops,
/// passing while asserting nothing (M-199). Select with `if name == "…" { … }`
/// instead.
///
/// `continue` and `break` are safe only inside a loop the body itself opens; a
/// bare one does not compile, which is why `return` is the only shape of this
/// mistake that survives to run.
#[macro_export]
macro_rules! for_each_reference_field {
    ($scalar:ty, |$name:ident, $field:ident| $body:block) => {{
        {
            let $name = "sphere";
            let $field = $crate::fields::Sphere::<$scalar>::canonical();
            $body
        }
        {
            let $name = "torus";
            let $field = $crate::fields::Torus::<$scalar>::canonical();
            $body
        }
        {
            let $name = "box_exact";
            let $field = $crate::fields::BoxExact::<$scalar>::canonical();
            $body
        }
        {
            let $name = "csg_difference";
            let $field = $crate::fields::csg_difference::<$scalar>();
            $body
        }
        {
            let $name = "thin_plate";
            let $field = $crate::fields::ThinPlate::<$scalar>::canonical();
            $body
        }
        {
            let $name = "gyroid";
            let $field = $crate::fields::capped_gyroid::<$scalar>();
            $body
        }
        {
            let $name = "fbm_terrain";
            let $field = $crate::fields::FbmTerrain::<$scalar>::canonical();
            $body
        }
        {
            let $name = "noise_cavity";
            let $field = $crate::fields::noise_cavity::<$scalar>();
            $body
        }
    }};
}

/// The domain half-extent shared by the five compact fields.
const COMPACT_DOMAIN: f64 = 2.0;

#[inline]
fn cube_domain<R: Real>(half: f64) -> ([R; 3], [R; 3]) {
    let h = R::from_f64(half);
    ([-h, -h, -h], [h, h, h])
}

// ─── sphere ─────────────────────────────────────────────────────────────────

/// Exact signed distance to a sphere: `f(p) = |p − c| − r`.
///
/// The simplest closed surface, and the baseline every algorithm is checked
/// against first. `|∇f| == 1` everywhere except the centre, where the gradient is
/// undefined; that is a single point and no grid corner in the canonical domain
/// lands on it unless the grid is both odd-sized and centred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sphere<R: Real> {
    /// Centre.
    pub center: [R; 3],
    /// Radius.
    pub radius: R,
}

impl<R: Real> Sphere<R> {
    /// Unit sphere at the origin.
    #[must_use]
    pub fn canonical() -> Self {
        Self {
            center: [R::ZERO; 3],
            radius: R::ONE,
        }
    }
}

impl<R: Real> Default for Sphere<R> {
    fn default() -> Self {
        Self::canonical()
    }
}

impl<R: Real> Sdf for Sphere<R> {
    type Scalar = R;

    #[inline]
    fn sample(&self, p: [R; 3]) -> R {
        let d = sub(p, self.center);
        length(d) - self.radius
    }

    #[inline]
    fn gradient(&self, p: [R; 3]) -> [R; 3] {
        let d = sub(p, self.center);
        let len = length(d);
        scale(d, len.recip())
    }
}

impl<R: Real> ReferenceField for Sphere<R> {
    const NAME: &'static str = "sphere";
    fn domain(&self) -> ([R; 3], [R; 3]) {
        cube_domain(COMPACT_DOMAIN)
    }
    fn closed_in_domain(&self) -> bool {
        true
    }
    fn expected_euler(&self) -> Option<i64> {
        Some(2)
    }
    /// Exact. `|p − c| − r` is the distance to the sphere, inside and out.
    fn bound(&self) -> FieldBound {
        FieldBound::Exact
    }
}

// ─── torus ──────────────────────────────────────────────────────────────────

/// Exact signed distance to a torus.
///
/// # Orientation
///
/// **The ring lies in the xz-plane and the axis of revolution is +y.** Stated
/// here because it is a convention, not a fact, and every visual comparison
/// depends on it.
///
/// ```text
/// s = |(p.x, p.z)|              distance from the axis
/// q = (s − major, p.y)
/// f = |q| − minor
/// ```
///
/// `|∇f| == 1` exactly. Undefined on the y-axis (`s == 0`) and on the tube's core
/// circle (`|q| == 0`); both are measure zero.
///
/// The only reference field with a non-trivial genus known in closed form, which
/// makes it the one that proves a validity harness computes `χ` rather than
/// hard-coding `2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Torus<R: Real> {
    /// Centre.
    pub center: [R; 3],
    /// Distance from the centre to the tube's core circle.
    pub major: R,
    /// Tube radius.
    pub minor: R,
}

impl<R: Real> Torus<R> {
    /// Major radius 1, minor radius 0.3, at the origin.
    #[must_use]
    pub fn canonical() -> Self {
        Self {
            center: [R::ZERO; 3],
            major: R::ONE,
            minor: R::from_f64(0.3),
        }
    }
}

impl<R: Real> Default for Torus<R> {
    fn default() -> Self {
        Self::canonical()
    }
}

impl<R: Real> Sdf for Torus<R> {
    type Scalar = R;

    #[inline]
    fn sample(&self, p: [R; 3]) -> R {
        let d = sub(p, self.center);
        let s = (d[0] * d[0] + d[2] * d[2]).sqrt();
        let q = [s - self.major, d[1]];
        (q[0] * q[0] + q[1] * q[1]).sqrt() - self.minor
    }

    #[inline]
    fn gradient(&self, p: [R; 3]) -> [R; 3] {
        let d = sub(p, self.center);
        let s = (d[0] * d[0] + d[2] * d[2]).sqrt();
        let q = [s - self.major, d[1]];
        let qlen = (q[0] * q[0] + q[1] * q[1]).sqrt();
        let radial = q[0] / qlen; // d|q| / ds
        [radial * (d[0] / s), q[1] / qlen, radial * (d[2] / s)]
    }
}

impl<R: Real> ReferenceField for Torus<R> {
    const NAME: &'static str = "torus";
    fn domain(&self) -> ([R; 3], [R; 3]) {
        cube_domain(COMPACT_DOMAIN)
    }
    fn closed_in_domain(&self) -> bool {
        true
    }
    fn expected_euler(&self) -> Option<i64> {
        Some(0) // genus 1
    }
    /// Exact. The revolved-circle form is a true distance in both regions.
    fn bound(&self) -> FieldBound {
        FieldBound::Exact
    }
}

// ─── box ────────────────────────────────────────────────────────────────────

/// `q = |p − c| − b`, shared by [`BoxExact`] and [`ThinPlate`] so the two cannot
/// drift apart.
#[inline]
fn box_q<R: Real>(p: [R; 3], center: [R; 3], half: [R; 3]) -> ([R; 3], [R; 3]) {
    let d = sub(p, center);
    let q = [
        d[0].abs() - half[0],
        d[1].abs() - half[1],
        d[2].abs() - half[2],
    ];
    (d, q)
}

#[inline]
fn box_sample<R: Real>(p: [R; 3], center: [R; 3], half: [R; 3]) -> R {
    let (_, q) = box_q(p, center, half);
    let outside = [q[0].max(R::ZERO), q[1].max(R::ZERO), q[2].max(R::ZERO)];
    let inside = q[0].max(q[1]).max(q[2]).min(R::ZERO);
    length(outside) + inside
}

#[inline]
fn box_gradient<R: Real>(p: [R; 3], center: [R; 3], half: [R; 3]) -> [R; 3] {
    let (d, q) = box_q(p, center, half);
    let max_q = q[0].max(q[1]).max(q[2]);

    if max_q > R::ZERO {
        // Exterior. `max(q, 0)` is non-zero here, so the division is safe --
        // which is exactly why the branch is `> 0` and not `>= 0`.
        let m = [q[0].max(R::ZERO), q[1].max(R::ZERO), q[2].max(R::ZERO)];
        let inv = length(m).recip();
        [
            d[0].signum() * m[0] * inv,
            d[1].signum() * m[1] * inv,
            d[2].signum() * m[2] * inv,
        ]
    } else {
        // Interior and on the surface. `|max(q, 0)|` is identically zero in a
        // neighbourhood, so its gradient is exactly zero and `f` reduces to
        // `q_j` for `j = argmax q`: the outward normal of the nearest face.
        let j = if q[0] >= q[1] && q[0] >= q[2] {
            0
        } else if q[1] >= q[2] {
            1
        } else {
            2
        };
        let mut g = [R::ZERO; 3];
        g[j] = d[j].signum();
        g
    }
}

/// Exact signed distance to an axis-aligned box.
///
/// ```text
/// q = |p − c| − b
/// f = |max(q, 0)| + min(max(q.x, q.y, q.z), 0)
/// ```
///
/// This is the **exact** form, not the cheaper `max(q.x, q.y, q.z)` bound. The
/// bound underestimates distance outside near an edge or corner, and dual
/// contouring's Hermite normals would inherit that error exactly where sharp
/// features are supposed to be recovered.
///
/// The discriminator between surface nets and dual contouring: surface nets
/// rounds these corners, dual contouring holds them.
///
/// # Where the gradient is not unique
///
/// Two loci, and the code makes a *deterministic selection from the
/// subdifferential* at each rather than inventing a value.
///
/// - **The interior medial axis**, where two or more components of `q` tie for
///   the maximum, i.e. `p` is equidistant from two faces. The true
///   subdifferential is the convex hull of the tied face normals; this returns
///   the tied axis of lowest index, x before y before z.
/// - **Box edges and corners**, where the exterior limit depends on the approach
///   direction. The gradient is therefore *discontinuous* across a box edge —
///   which is correct for an exact box distance, since any single-valued
///   selection must be, and is precisely the behaviour Hermite extraction needs
///   to see in order to recover the edge.
///
/// `p_i == 0` looks like it should be a third case, because `|p_i|` is not
/// differentiable there, and is not: at `p_i == 0` the component `q_i` is at its
/// minimum `−b_i`, so it is never the argmax for a non-degenerate box, and in the
/// exterior branch `max(q, 0)_i` is zero from both sides.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoxExact<R: Real> {
    /// Centre.
    pub center: [R; 3],
    /// Half-extents along x, y, z.
    pub half_extents: [R; 3],
}

impl<R: Real> BoxExact<R> {
    /// The `[-1, 1]³` cube.
    #[must_use]
    pub fn canonical() -> Self {
        Self {
            center: [R::ZERO; 3],
            half_extents: [R::ONE; 3],
        }
    }
}

impl<R: Real> Default for BoxExact<R> {
    fn default() -> Self {
        Self::canonical()
    }
}

impl<R: Real> Sdf for BoxExact<R> {
    type Scalar = R;

    #[inline]
    fn sample(&self, p: [R; 3]) -> R {
        box_sample(p, self.center, self.half_extents)
    }

    #[inline]
    fn gradient(&self, p: [R; 3]) -> [R; 3] {
        box_gradient(p, self.center, self.half_extents)
    }
}

impl<R: Real> ReferenceField for BoxExact<R> {
    const NAME: &'static str = "box_exact";
    fn domain(&self) -> ([R; 3], [R; 3]) {
        cube_domain(COMPACT_DOMAIN)
    }
    fn closed_in_domain(&self) -> bool {
        true
    }
    fn expected_euler(&self) -> Option<i64> {
        Some(2)
    }
    /// Exact, and that is what distinguishes this from a cheaper box: the
    /// `max(d, 0).length() + min(max(d), 0)` form is a true distance outside
    /// *and* inside, where the common `max(d)` form is only correct outside.
    fn bound(&self) -> FieldBound {
        FieldBound::Exact
    }
}

// ─── thin plate ─────────────────────────────────────────────────────────────

/// A slab deliberately thinner than one grid cell — the sub-voxel feature case.
///
/// Geometrically a box, and it shares [`BoxExact`]'s formula and gradient
/// exactly, with no second implementation. It is a distinct type so that its
/// intent and its cell-size-relative constructor have somewhere to live.
///
/// # Why it exists
///
/// At the cell size it was built for, both endpoints of every vertical grid edge
/// fall outside the plate, so classic marching cubes sees no sign change and
/// emits **nothing at all**. An extractor that counts edge *intersections*
/// instead of corner *signs* sees two crossings on that edge and recovers the
/// plate. That contrast is the whole demonstration.
///
/// # Resolution caveat
///
/// [`expected_euler`](ReferenceField::expected_euler) reports `Some(2)`, which is
/// what a *resolving* extraction must produce. At the cell size the plate was
/// constructed for, a corner-sign extractor correctly produces an empty mesh;
/// that is a separate, explicit expectation and not a tolerated failure. Sample
/// at roughly four times the constructing resolution to resolve it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThinPlate<R: Real> {
    /// Centre.
    pub center: [R; 3],
    /// Half-extents. The y component is the half-thickness.
    pub half_extents: [R; 3],
}

impl<R: Real> ThinPlate<R> {
    /// Plate thickness as a fraction of the cell size it is built for.
    ///
    /// Below `1.0` so that no grid phase can ever put a corner inside the plate,
    /// with margin; far enough above zero that a four-fold refinement resolves
    /// it comfortably.
    pub const THICKNESS_IN_CELLS: f64 = 0.4;

    /// A plate `THICKNESS_IN_CELLS × h` thick, spanning `[-1, 1]` in x and z.
    #[must_use]
    pub fn for_cell_size(h: R) -> Self {
        let half_thickness = h * R::from_f64(Self::THICKNESS_IN_CELLS * 0.5);
        Self {
            center: [R::ZERO; 3],
            half_extents: [R::ONE, half_thickness, R::ONE],
        }
    }

    /// The cell size [`canonical`](Self::canonical) is built for: a 64³ grid
    /// spanning this field's domain.
    pub const CANONICAL_CELL_SIZE: f64 = 2.0 * COMPACT_DOMAIN / 64.0;

    /// A plate sub-voxel at 64³ over the `[-2, 2]³` domain.
    #[must_use]
    pub fn canonical() -> Self {
        Self::for_cell_size(R::from_f64(Self::CANONICAL_CELL_SIZE))
    }
}

impl<R: Real> Default for ThinPlate<R> {
    fn default() -> Self {
        Self::canonical()
    }
}

impl<R: Real> Sdf for ThinPlate<R> {
    type Scalar = R;

    #[inline]
    fn sample(&self, p: [R; 3]) -> R {
        box_sample(p, self.center, self.half_extents)
    }

    #[inline]
    fn gradient(&self, p: [R; 3]) -> [R; 3] {
        box_gradient(p, self.center, self.half_extents)
    }
}

impl<R: Real> ReferenceField for ThinPlate<R> {
    const NAME: &'static str = "thin_plate";
    fn domain(&self) -> ([R; 3], [R; 3]) {
        cube_domain(COMPACT_DOMAIN)
    }
    fn closed_in_domain(&self) -> bool {
        true
    }
    fn expected_euler(&self) -> Option<i64> {
        Some(2)
    }
    /// Exact. A slab is a box, and the same exact form is used.
    fn bound(&self) -> FieldBound {
        FieldBound::Exact
    }
}

// ─── combinators ────────────────────────────────────────────────────────────

/// Set difference `a − b`, as `max(f_a, −f_b)`.
///
/// Set-exact for any pair of implicit functions, whether or not either is a true
/// distance field.
///
/// # Gradient
///
/// The gradient of whichever operand is active, and it is **discontinuous at the
/// seam** `f_a == −f_b`. That discontinuity is the concave sharp edge, and is the
/// thing a sharp-feature extractor is supposed to reproduce. On the seam itself
/// the first operand wins — again a deterministic selection from the
/// subdifferential rather than an invented value.
///
/// # Distance
///
/// If both operands are exact distances then so is this, away from the seam. The
/// *value* near a concave seam is only a lower bound on true distance, so an
/// accuracy harness should measure against geometry rather than against
/// `|sample|`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Difference<A, B> {
    /// The solid being cut.
    pub a: A,
    /// The tool being subtracted.
    pub b: B,
}

impl<A: Sdf, B: Sdf<Scalar = A::Scalar>> Sdf for Difference<A, B> {
    type Scalar = A::Scalar;

    #[inline]
    fn sample(&self, p: [Self::Scalar; 3]) -> Self::Scalar {
        let (fa, fb) = (self.a.sample(p), self.b.sample(p));
        if fa >= -fb { fa } else { -fb }
    }

    #[inline]
    fn gradient(&self, p: [Self::Scalar; 3]) -> [Self::Scalar; 3] {
        let (fa, fb) = (self.a.sample(p), self.b.sample(p));
        if fa >= -fb {
            self.a.gradient(p)
        } else {
            let g = self.b.gradient(p);
            [-g[0], -g[1], -g[2]]
        }
    }
}

/// Set intersection `a ∩ b`, as `max(f_a, f_b)`.
///
/// Set-exact for implicit functions even when neither operand is a distance
/// field, which is what makes it usable to cap [`Gyroid`]. Gradient is that of
/// the active operand; the first wins on the seam.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Intersection<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
}

impl<A: Sdf, B: Sdf<Scalar = A::Scalar>> Sdf for Intersection<A, B> {
    type Scalar = A::Scalar;

    #[inline]
    fn sample(&self, p: [Self::Scalar; 3]) -> Self::Scalar {
        let (fa, fb) = (self.a.sample(p), self.b.sample(p));
        if fa >= fb { fa } else { fb }
    }

    #[inline]
    fn gradient(&self, p: [Self::Scalar; 3]) -> [Self::Scalar; 3] {
        let (fa, fb) = (self.a.sample(p), self.b.sample(p));
        if fa >= fb {
            self.a.gradient(p)
        } else {
            self.b.gradient(p)
        }
    }
}

/// Set union `a ∪ b`, as `min(f_a, f_b)`.
///
/// The combinator the crate went without until E-216 needed it, and the absence
/// is worth a sentence: `Difference` caps `csg_difference` and `Intersection`
/// caps `Gyroid`, so both had a reference field asking for them. **Nothing in
/// the suite unions anything**, so the most basic CSG operation was the one
/// missing — a property of the fixtures rather than of the design (M-240).
///
/// # Distance, and why this one is the safe direction
///
/// `min` of two 1-Lipschitz functions is 1-Lipschitz, and it **never
/// overestimates** the true distance: away from the seam it is exact, and near
/// one it is a conservative lower bound. That is the direction a sphere tracer
/// can survive — a step that is too short only costs iterations. `max`, which
/// [`Intersection`] and [`Difference`] use, overestimates near concave seams and
/// is the direction that lets a tracer step *through* a surface. Phase 11's
/// `F-001` is where that distinction gets a type.
///
/// Gradient is that of the active operand; the first wins on the seam, which is
/// the same deterministic selection from the subdifferential the other two make.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Union<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
}

impl<A: Sdf, B: Sdf<Scalar = A::Scalar>> Sdf for Union<A, B> {
    type Scalar = A::Scalar;

    #[inline]
    fn sample(&self, p: [Self::Scalar; 3]) -> Self::Scalar {
        let (fa, fb) = (self.a.sample(p), self.b.sample(p));
        if fa <= fb { fa } else { fb }
    }

    #[inline]
    fn gradient(&self, p: [Self::Scalar; 3]) -> [Self::Scalar; 3] {
        let (fa, fb) = (self.a.sample(p), self.b.sample(p));
        if fa <= fb {
            self.a.gradient(p)
        } else {
            self.b.gradient(p)
        }
    }
}

/// Union with a rounded seam of radius `k`, as
/// [`smooth_min`](crate::brush::smooth_min).
///
/// **The parameter a level designer actually reaches for.** A hard [`Union`]
/// leaves a crease where two primitives meet; `k` is how wide the fillet is, in
/// world units, and sweeping it from zero is the difference between two spheres
/// touching and one blob.
///
/// # It is not a distance field, and the blend is where it stops being one
///
/// `smooth_min` is smaller than `min` by up to `k/4`, so the value understates
/// distance inside the blend region. It stays a conservative *lower* bound —
/// the safe direction, as for [`Union`] — but an accuracy harness must measure
/// against geometry rather than against `|sample|`. The gradient is the exact
/// derivative of the blend rather than either operand's, because on this seam
/// there is no active operand to pick: that is what "smooth" means.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmoothUnion<A, B, R> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
    /// Blend radius, in world units. Zero degenerates to a hard [`Union`].
    pub k: R,
}

impl<A: Sdf, B: Sdf<Scalar = A::Scalar>> Sdf for SmoothUnion<A, B, A::Scalar> {
    type Scalar = A::Scalar;

    #[inline]
    fn sample(&self, p: [Self::Scalar; 3]) -> Self::Scalar {
        crate::brush::smooth_min(self.a.sample(p), self.b.sample(p), self.k)
    }

    /// Central differences, because the blend has no active operand.
    ///
    /// [`Union`] can hand back whichever operand is active; here both contribute
    /// everywhere inside the blend, so the analytic gradient would be a chain
    /// rule through `smooth_min`'s own `h`. Differencing the composed field is
    /// shorter, is exact to `O(h²)`, and cannot disagree with `sample` — which
    /// an independently written analytic form could.
    #[inline]
    fn gradient(&self, p: [Self::Scalar; 3]) -> [Self::Scalar; 3] {
        let e = Self::Scalar::from_f64(1e-4);
        let at = |q: [Self::Scalar; 3]| self.sample(q);
        let two = Self::Scalar::ONE + Self::Scalar::ONE;
        [
            (at([p[0] + e, p[1], p[2]]) - at([p[0] - e, p[1], p[2]])) / (two * e),
            (at([p[0], p[1] + e, p[2]]) - at([p[0], p[1] - e, p[2]])) / (two * e),
            (at([p[0], p[1], p[2] + e]) - at([p[0], p[1], p[2] - e])) / (two * e),
        ]
    }
}

/// A box with a sphere subtracted from one corner.
pub type CsgDifference<R> = Difference<BoxExact<R>, Sphere<R>>;

/// The `[-1, 1]³` cube minus a sphere of radius `0.75` centred at
/// `(0.6, 0.6, 0.6)`.
///
/// The sphere reaches the `+x`, `+y` and `+z` faces (`0.4 < 0.75`) and none of
/// the opposite three (`1.6 > 0.75`), so it scoops the `+++` corner and leaves
/// **one closed surface with `χ = 2`**: convex sharp edges from the box, and
/// three concave circular seams where the sphere cuts each face.
///
/// Concave sharp edges are the harder half of sharp-feature extraction, and this
/// is the field that shows whether they survive.
#[must_use]
pub fn csg_difference<R: Real>() -> CsgDifference<R> {
    Difference {
        a: BoxExact::canonical(),
        b: Sphere {
            center: [R::from_f64(0.6); 3],
            radius: R::from_f64(0.75),
        },
    }
}

impl<R: Real> ReferenceField for CsgDifference<R> {
    const NAME: &'static str = "csg_difference";
    fn domain(&self) -> ([R; 3], [R; 3]) {
        cube_domain(COMPACT_DOMAIN)
    }
    fn closed_in_domain(&self) -> bool {
        true
    }
    fn expected_euler(&self) -> Option<i64> {
        Some(2)
    }
    /// **Not exact, and this declaration is why F-001 existed.**
    ///
    /// It returned `true` with `// away from the seam` beside it — a comment
    /// admitting the invariant is false, on a function returning true. The field
    /// is `max(box, −sphere)`, and `max` of two exact distances is not an exact
    /// distance: near a concave seam it **overestimates**, which is the direction
    /// that lets a sphere tracer step through the surface. That seam is exactly
    /// where A-014d located this field's coincident polygons.
    ///
    /// **`q` is declared conservatively rather than derived.** The 2023 paper's
    /// Theorem 2 gives `min(q_f, q_g)` for `max` on the interior, which is 1 for
    /// two exact operands, but its Theorem 6 shows the whole-space bound degrades
    /// by a geometry-dependent factor this crate does not compute. F-003 is the
    /// ticket that computes it from the operands; until then a loose bound that
    /// is true beats a tight one that is guessed.
    fn bound(&self) -> FieldBound {
        FieldBound::Underestimate { q: 0.5 }
    }
}

// ─── gyroid ─────────────────────────────────────────────────────────────────

/// Schoen's gyroid, in its standard trigonometric approximation.
///
/// ```text
/// g(p) = sin(sx)·cos(sy) + sin(sy)·cos(sz) + sin(sz)·cos(sx) − iso
/// ```
///
/// # This is not a distance field
///
/// `g` is an implicit function whose zero set is the surface, and nothing more.
/// `|∇g|` varies with position, so its [`bound`](ReferenceField::bound) is
/// [`Lipschitz`](FieldBound::Lipschitz) rather than
/// [`Exact`](FieldBound::Exact), and `|g(v)|` must not be used as a distance.
///
/// It is deliberately **not** normalised by `|∇g|`. Doing so would buy nothing
/// and cost something: marching cubes' linear edge interpolation `t = a/(a−b)` is
/// a first-order root find whose error goes as `O(h²·f″/f′)` — a *ratio*, which
/// scaling the field leaves unchanged to leading order — while dividing by
/// `|∇g|` introduces a second failure mode wherever `|∇g| → 0`. The extracted
/// vertex is first-order accurate either way.
///
/// Triply periodic with period `2π/scale`, and high genus at any useful size,
/// which is what makes it the field that stresses topology.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Gyroid<R: Real> {
    /// Spatial frequency. The period is `2π / scale`.
    pub scale: R,
    /// Level set to extract. Zero gives the balanced surface.
    pub iso: R,
}

impl<R: Real> Gyroid<R> {
    /// Unit scale, zero iso value: period `2π`.
    #[must_use]
    pub fn canonical() -> Self {
        Self {
            scale: R::ONE,
            iso: R::ZERO,
        }
    }
}

impl<R: Real> Default for Gyroid<R> {
    fn default() -> Self {
        Self::canonical()
    }
}

impl<R: Real> Sdf for Gyroid<R> {
    type Scalar = R;

    #[inline]
    fn sample(&self, p: [R; 3]) -> R {
        let (a, b, c) = (self.scale * p[0], self.scale * p[1], self.scale * p[2]);
        a.sin() * b.cos() + b.sin() * c.cos() + c.sin() * a.cos() - self.iso
    }

    #[inline]
    fn gradient(&self, p: [R; 3]) -> [R; 3] {
        let (a, b, c) = (self.scale * p[0], self.scale * p[1], self.scale * p[2]);
        let (sa, ca) = (a.sin(), a.cos());
        let (sb, cb) = (b.sin(), b.cos());
        let (sc, cc) = (c.sin(), c.cos());
        [
            self.scale * (ca * cb - sc * sa),
            self.scale * (cb * cc - sa * sb),
            self.scale * (cc * ca - sb * sc),
        ]
    }
}

/// A gyroid capped to a sphere, so that it is a closed surface at all.
pub type CappedGyroid<R> = Intersection<Gyroid<R>, Sphere<R>>;

/// The canonical gyroid entry: [`Gyroid::canonical`] intersected with a sphere of
/// radius 6, over the `[-7, 7]³` domain.
///
/// The cap exists because an uncapped triply periodic surface has **boundary**
/// wherever the sampling box cuts it, so its Euler characteristic is neither `2`
/// nor predictable. Capping makes it closed and therefore checkable; its genus
/// still is not known in closed form, which is why
/// [`expected_euler`](ReferenceField::expected_euler) returns `None` and the
/// observed value belongs in a golden fixture instead.
///
/// The cap sits a full unit inside the domain wall, so the surface never touches
/// it. `max`-style intersection is set-exact even though `Gyroid` is not a
/// distance field; the one caveat is that the two operands have very different
/// magnitudes near the seam, so an edge straddling it interpolates poorly. That
/// is inherent to CSG on mismatched-scale fields rather than a defect here.
#[must_use]
pub fn capped_gyroid<R: Real>() -> CappedGyroid<R> {
    Intersection {
        a: Gyroid::canonical(),
        b: Sphere {
            center: [R::ZERO; 3],
            radius: R::from_f64(6.0),
        },
    }
}

impl<R: Real> ReferenceField for CappedGyroid<R> {
    const NAME: &'static str = "gyroid";
    fn domain(&self) -> ([R; 3], [R; 3]) {
        cube_domain(7.0)
    }
    fn closed_in_domain(&self) -> bool {
        true
    }
    fn expected_euler(&self) -> Option<i64> {
        None // genus depends on how many tunnels the cap encloses
    }
    /// Lipschitz, and not a distance — with the constant **derived** rather than
    /// guessed, because the first attempt at this guessed `√3` and was wrong.
    ///
    /// The field is `sin(a)cos(b) + sin(b)cos(c) + sin(c)cos(a)` at unit scale,
    /// whose partials are
    ///
    /// ```text
    /// ∂/∂x = cos(a)cos(b) − sin(c)sin(a)
    /// ```
    ///
    /// and permutations. Each is a difference of two products of quantities
    /// bounded by one, so `|∂/∂x| ≤ 2` and `|∇g| ≤ 2√3 ≈ 3.464`. Loose — the
    /// sampled maximum is about 1.7 — but **sound**, which is the property that
    /// matters: a sphere tracer dividing by a constant smaller than the true one
    /// steps through the surface, so an under-declared Lipschitz bound is worse
    /// than none.
    ///
    /// Capped by an exact sphere, and `max` cannot raise the constant above the
    /// larger of its operands', so the cap does not change this.
    fn bound(&self) -> FieldBound {
        FieldBound::Lipschitz {
            l: 3.464_101_615_137_754_6,
        }
    }
}

// ─── volumetric noise ───────────────────────────────────────────────────────

/// A single octave of 3D gradient noise, as a volume rather than a heightfield.
///
/// `f(p) = perlin(frequency · p) − iso`, so the surface is the boundary between
/// the regions where the noise is below and above `iso` — the blobby, branching
/// shape a voxel game gets from carving caves out of solid rock, and the shape
/// [`FbmTerrain`] cannot produce because it only ever samples a horizontal plane.
///
/// # Why this field exists, and why noise rather than another analytic solid
///
/// Every other field here has an *interior ambiguity* rate of exactly zero.
/// Measured over all seven at 17³, 33³ and 65³: **not one of 68,385 surface cells
/// has six body saddles**, and only five cells in the whole sweep reach even five
/// (M-208). So the trilinear interpolant's tunnel case — the thing MC33's interior
/// rule exists for — was unreachable by this crate's own test suite, and a
/// per-cell proof was all A-002's series could ever have got.
///
/// Smooth analytic solids do not produce tunnels because they are too smooth: a
/// tunnel needs the field to reverse twice across one cell, and a sphere or a box
/// never does. Both papers behind MC33's corrections used **randomly generated
/// scalar fields** for exactly this reason — Custodio et al. count their
/// non-manifold case *"once in 10000"* random 5×5×5 grids, and Grosso's own tunnel
/// statistics come from CT data and random volumes. Gradient noise is the smooth,
/// deterministic, analytically differentiable version of the same thing.
///
/// # The parameters are searched, and `iso` may not be zero
///
/// `frequency` and `iso` were **searched** rather than chosen — over 610
/// combinations, keeping only those whose mesh is *closed* at 17³, 25³ and 33³
/// **and** which reach six body saddles at all three rather than at a lucky one.
/// 97 combinations qualify, so this is a plateau rather than a knife edge; these
/// are the best of them (M-209).
///
/// **`iso` is deliberately not zero, and that is not a tuning choice.** Perlin
/// noise is *exactly* zero at every point of its own integer lattice — measured,
/// `0.000e0` at all six lattice points probed, against `1.5e-1`-ish just off it.
/// So the zero level set of gradient noise **contains the whole lattice**, which
/// makes the surface pass through a regular array of its own critical points; the
/// first version of this field used `iso = 0` and tripped
/// [`Sdf::gradient`](crate::Sdf)'s zero-gradient assertion during extraction
/// (M-210). Any non-zero level avoids the lattice entirely.
///
/// # Undersampling is the point, not a defect
///
/// The features here are about `1 / frequency ≈ 0.29` across while the coarsest
/// golden grid has cells `0.25` across, so this field is **deliberately sampled
/// near its own feature size**. That is not sloppiness: a tunnel requires the
/// field to reverse twice across a single cell, which cannot happen when the cell
/// is small relative to the features. Grosso reports the same relationship from
/// the other end — refining a volume with 16 tunnels twice leaves **one** — and
/// it reproduces here, six-saddle cells thinning to **zero** by 65³.
///
/// The frequency is also the **gentlest** that still reaches the configuration at
/// all three golden resolutions, and that is deliberate too. Volumetric noise
/// sampled near its feature size is hard on every extractor, not just on the
/// interior rule: at frequency `4.9` this field reaches four times as many
/// six-saddle cells and roughly **quadruples** Subgrid Marching Tetrahedra's
/// non-manifold and flipped-edge counts, which is collateral rather than coverage
/// (M-209). The ticket asked for a field that reaches the configuration, not one
/// that maximises it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NoiseVolume<R: Real> {
    /// Which noise field. Any `u32` gives a different but equally valid one.
    pub seed: u32,
    /// Spatial frequency. The features are about `1 / frequency` across.
    pub frequency: R,
    /// Level set to extract. Zero gives the noise's own balanced surface.
    pub iso: R,
}

impl<R: Real> NoiseVolume<R> {
    /// The searched parameters: frequency `3.45`, iso `0.25`.
    ///
    /// The seed is [`FbmTerrain`]'s, deliberately — the two fields drive the same
    /// generator through different domains (a horizontal plane against the whole
    /// volume) and sharing it means one fewer arbitrary constant, not one more
    /// coincidence.
    #[must_use]
    pub fn canonical() -> Self {
        Self {
            seed: 0x5EED_1234,
            frequency: R::from_f64(3.45),
            iso: R::from_f64(0.25),
        }
    }
}

impl<R: Real> Default for NoiseVolume<R> {
    fn default() -> Self {
        Self::canonical()
    }
}

impl<R: Real> Sdf for NoiseVolume<R> {
    type Scalar = R;

    #[inline]
    fn sample(&self, p: [R; 3]) -> R {
        let (n, _) = noise::perlin(scale(p, self.frequency), self.seed);
        n - self.iso
    }

    #[inline]
    fn gradient(&self, p: [R; 3]) -> [R; 3] {
        // Chain rule: the noise is evaluated at `frequency · p`, so its gradient
        // scales by `frequency`. Analytic, from the same evaluation as the value —
        // `perlin` returns both, which is why this field can be `ReferenceField`
        // at all (`fields/tests.rs` checks the analytic gradient against a central
        // difference, and a finite-difference stand-in would be comparing a
        // difference with itself).
        let (_, g) = noise::perlin(scale(p, self.frequency), self.seed);
        scale(g, self.frequency)
    }
}

/// [`NoiseVolume`] capped to a sphere, so that it is a closed surface.
pub type NoiseCavity<R> = Intersection<NoiseVolume<R>, Sphere<R>>;

/// The canonical volumetric-noise entry: [`NoiseVolume::canonical`] intersected
/// with a sphere of radius `1.5`, over the `[-2, 2]³` domain.
///
/// Capped for [`CappedGyroid`]'s reason: an uncapped noise level set has
/// **boundary** wherever the sampling box cuts it, so its Euler characteristic is
/// neither `2` nor predictable. The cap sits half a unit inside the domain wall.
/// Its genus is not known in closed form — it depends on how many noise blobs the
/// sphere happens to enclose — so
/// [`expected_euler`](ReferenceField::expected_euler) is `None` and the observed
/// value belongs in a golden fixture, exactly as for the gyroid.
#[must_use]
pub fn noise_cavity<R: Real>() -> NoiseCavity<R> {
    Intersection {
        a: NoiseVolume::canonical(),
        b: Sphere {
            center: [R::ZERO; 3],
            radius: R::from_f64(1.5),
        },
    }
}

impl<R: Real> ReferenceField for NoiseCavity<R> {
    const NAME: &'static str = "noise_cavity";
    fn domain(&self) -> ([R; 3], [R; 3]) {
        cube_domain(2.0)
    }
    fn closed_in_domain(&self) -> bool {
        true
    }
    fn expected_euler(&self) -> Option<i64> {
        None // genus depends on how many noise blobs the cap encloses
    }
    /// **Unbounded, because the constant has not been derived and a guess was
    /// caught (M-244).**
    ///
    /// This first declared `Lipschitz { l: 2.598 }`, reasoning loosely from
    /// Perlin noise's gradient bound. `every_field_meets_the_bound_it_declares`
    /// measured `|∇f|` reaching **7.73** and failed it on the first run — which
    /// is the whole reason a declared bound is checked rather than trusted, and
    /// the same defect as the `// away from the seam` comment this ticket
    /// replaced, made by the person replacing it.
    ///
    /// A real constant is derivable from the noise's amplitude, frequency and
    /// octave count, and **F-002 is where that gets established** — a sampled
    /// maximum is a *lower* bound on the supremum, so measuring one and
    /// declaring it would be unsound in exactly the dangerous direction. Until
    /// then this claims nothing beyond the sign, which is all extraction reads.
    fn bound(&self) -> FieldBound {
        FieldBound::Unbounded
    }
}

// ─── terrain ────────────────────────────────────────────────────────────────

/// A fractal heightfield: `f(p) = p.y − (base + amplitude · fbm(p.x, 0, p.z))`.
///
/// # Not a distance field, and not closed
///
/// The value is a vertical distance only, and its Lipschitz constant exceeds one
/// wherever the terrain is steep, so
/// its [`bound`](ReferenceField::bound) is [`Unbounded`](FieldBound::Unbounded). The
/// surface also leaves through the sides of the domain, so
/// [`closed_in_domain`](ReferenceField::closed_in_domain) is `false` and a caller
/// must not require zero boundary edges. It is the reference field for the case
/// a game actually meshes.
///
/// # Determinism
///
/// The noise underneath evaluates no transcendental at all — see the module it
/// comes from. Every default below is arbitrary but **committed**: changing any
/// of them invalidates every golden hash that mentions this field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FbmTerrain<R: Real> {
    /// Lattice hash seed.
    pub seed: u32,
    /// Number of octaves summed.
    pub octaves: u32,
    /// Frequency multiplier per octave. `2.0` is an exact binary scale, so no
    /// drift accumulates across octaves.
    pub lacunarity: R,
    /// Amplitude multiplier per octave.
    pub gain: R,
    /// Frequency of the first octave.
    pub frequency: R,
    /// Overall height scale.
    pub amplitude: R,
    /// Height the terrain varies about.
    pub base_height: R,
}

impl<R: Real> FbmTerrain<R> {
    /// The committed parameter set.
    #[must_use]
    pub fn canonical() -> Self {
        Self {
            seed: 0x5EED_1234,
            octaves: 4,
            lacunarity: R::from_f64(2.0),
            gain: R::from_f64(0.5),
            frequency: R::from_f64(0.25),
            amplitude: R::from_f64(2.0),
            base_height: R::ZERO,
        }
    }

    /// A provable bound on how far the surface can sit from
    /// [`base_height`](Self::base_height).
    ///
    /// Used by the sign tests, so that "this point is above the terrain" is a
    /// proof rather than a sampled guess that happens to hold today.
    #[must_use]
    pub fn height_bound(&self) -> R {
        self.amplitude.abs() * noise::fbm_bound::<R>(self.octaves, self.gain)
    }
}

impl<R: Real> Default for FbmTerrain<R> {
    fn default() -> Self {
        Self::canonical()
    }
}

impl<R: Real> Sdf for FbmTerrain<R> {
    type Scalar = R;

    #[inline]
    fn sample(&self, p: [R; 3]) -> R {
        let (n, _) = noise::fbm(
            [p[0], R::ZERO, p[2]],
            self.seed,
            self.octaves,
            self.lacunarity,
            self.gain,
            self.frequency,
        );
        p[1] - (self.base_height + self.amplitude * n)
    }

    #[inline]
    fn gradient(&self, p: [R; 3]) -> [R; 3] {
        let (_, g) = noise::fbm(
            [p[0], R::ZERO, p[2]],
            self.seed,
            self.octaves,
            self.lacunarity,
            self.gain,
            self.frequency,
        );
        [-self.amplitude * g[0], R::ONE, -self.amplitude * g[2]]
    }
}

impl<R: Real> ReferenceField for FbmTerrain<R> {
    const NAME: &'static str = "fbm_terrain";
    fn domain(&self) -> ([R; 3], [R; 3]) {
        cube_domain(8.0)
    }
    fn closed_in_domain(&self) -> bool {
        false // a heightfield exits through the sides
    }
    fn expected_euler(&self) -> Option<i64> {
        None // not closed, so there is nothing to assert
    }
    /// Unbounded. A heightfield's `y − h(x, z)` is a distance only where the
    /// surface is level; on a slope it understates by `cos θ`, and fbm has no
    /// slope bound worth asserting. The sign is right everywhere, which is all
    /// extraction needs.
    fn bound(&self) -> FieldBound {
        FieldBound::Unbounded
    }
}

#[cfg(test)]
mod tests;
