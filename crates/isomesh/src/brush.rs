//! Brush operations: the edits a sculpting tool makes to a field.
//!
//! A brush is a shape plus an operation. Applied to a field it produces another
//! field, so a stack of edits is a composition rather than a mutation of stored
//! data — which is what lets [`chunk::dirty`](crate::chunk::dirty) express an
//! edit as *two fields* and measure exactly what moved.
//!
//! # Commutativity, which is the whole reason this ticket has an acceptance test
//!
//! Whether edits can be reordered decides whether concurrent editing is possible
//! at all. G-003 measured it over all `8! = 40,320` orderings of eight brushes,
//! and the answer is three different answers:
//!
//! | edits | orderings agree? | why |
//! |---|---|---|
//! | all [`Add`](BrushOp::Add), or all [`Subtract`](BrushOp::Subtract) | **yes, bit-exactly** | `min` and `max` are commutative *and* associative in IEEE, with no rounding at all |
//! | mixed add and subtract | **no** | carving a hole and then filling it is a different solid from filling and then carving. Semantic, not numerical |
//! | any [`SmoothAdd`](BrushOp::SmoothAdd) | **no** | smooth-min is **not associative**, and not bit-commutative either |
//!
//! See M-36, M-37 and M-38. The practical consequence is in
//! [`BrushOp::commutes_with`].

use crate::{Real, Sdf, vec3};

/// A capsule: the set of points within `radius` of a segment.
///
/// The third brush shape, alongside [`Sphere`](crate::fields::Sphere) and
/// [`BoxExact`](crate::fields::BoxExact). An exact distance field — the classic
/// point-to-segment distance, which is what makes it usable as a brush without
/// distorting the neighbourhood it is combined into.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Capsule<R: Real> {
    /// One end of the segment.
    pub a: [R; 3],
    /// The other end.
    pub b: [R; 3],
    /// Distance from the segment that counts as inside.
    pub radius: R,
}

impl<R: Real> Sdf for Capsule<R> {
    type Scalar = R;

    fn sample(&self, p: [R; 3]) -> R {
        let ab = vec3::sub(self.b, self.a);
        let ap = vec3::sub(p, self.a);
        let denom = vec3::dot(ab, ab);
        // A zero-length capsule is a sphere, which is the right answer rather
        // than a degenerate case to reject.
        let t = if denom > R::ZERO {
            (vec3::dot(ap, ab) / denom).clamp(R::ZERO, R::ONE)
        } else {
            R::ZERO
        };
        vec3::length(vec3::sub(ap, vec3::scale(ab, t))) - self.radius
    }
}

/// What a brush does to the field it is applied to.
///
/// `PartialEq` but not `Eq`, because [`SmoothAdd`](Self::SmoothAdd) carries a
/// join width and floats are not totally ordered by equality.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BrushOp {
    /// Union: `min(field, shape)`. Adds material.
    ///
    /// `min` is commutative *and* associative in IEEE and introduces no rounding
    /// — it selects one of its inputs rather than computing a new value — so any
    /// number of `Add`s reorder bit-exactly.
    Add,
    /// Difference: `max(field, −shape)`. Removes material.
    ///
    /// Same argument as [`Add`](Self::Add): `max` is exact and associative, so
    /// `Subtract`s reorder among themselves bit-exactly. They do **not** commute
    /// with `Add`s, and that is geometry rather than arithmetic.
    Subtract,
    /// Union with a rounded join of width `k`, as a polynomial smooth minimum.
    ///
    /// Looks better and **breaks reordering twice over**. Smooth-min is not
    /// associative even in exact arithmetic — `smin(smin(a,b),c) ≠
    /// smin(a,smin(b,c))` — and it is not *bit*-commutative either, because
    /// swapping the arguments evaluates a different expression that agrees only
    /// to rounding. A stack containing one of these has an order that is part of
    /// its meaning: measured, all eight over 40,320 orderings give **40,317**
    /// distinct results.
    SmoothAdd {
        /// Join width, in world units. Larger is rounder.
        k: f64,
    },
}

impl BrushOp {
    /// Whether two operations can be swapped without changing the result.
    ///
    /// The honest answer rather than an optimistic one: only identical hard
    /// operations commute. A caller building a concurrent-editing protocol can
    /// reorder freely within a run of `Add`s or a run of `Subtract`s, and must
    /// preserve order across a boundary between them or around any smooth
    /// operation.
    #[must_use]
    pub fn commutes_with(self, other: Self) -> bool {
        matches!(
            (self, other),
            (Self::Add, Self::Add) | (Self::Subtract, Self::Subtract)
        )
    }
}

/// One edit: a shape and what to do with it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Brush<S> {
    /// The shape being applied.
    pub shape: S,
    /// What to do with it.
    pub op: BrushOp,
}

impl<S> Brush<S> {
    /// Add this shape.
    pub const fn add(shape: S) -> Self {
        Self {
            shape,
            op: BrushOp::Add,
        }
    }

    /// Subtract this shape.
    pub const fn subtract(shape: S) -> Self {
        Self {
            shape,
            op: BrushOp::Subtract,
        }
    }

    /// Add this shape with a rounded join.
    pub const fn smooth_add(shape: S, k: f64) -> Self {
        Self {
            shape,
            op: BrushOp::SmoothAdd { k },
        }
    }
}

/// Apply one brush to a field value.
///
/// `field` is what the field says at this point already; `shape` is what the
/// brush's shape says there.
#[must_use]
pub fn apply<R: Real>(op: BrushOp, field: R, shape: R) -> R {
    match op {
        BrushOp::Add => field.min(shape),
        BrushOp::Subtract => field.max(-shape),
        BrushOp::SmoothAdd { k } => smooth_min(field, shape, R::from_f64(k)),
    }
}

/// Polynomial smooth minimum.
///
/// ```text
/// h = clamp(0.5 + 0.5·(b − a)/k, 0, 1)
/// smin = mix(b, a, h) − k·h·(1 − h)
/// ```
///
/// # Two separate failures of reordering, both measured
///
/// **Not associative**, in exact arithmetic — this is the one that matters, and
/// it is why a stack of smooth unions depends on the order it was folded in.
///
/// **Commutative algebraically but not bit-exactly.** Swapping `a` and `b` maps
/// `h` to `1 − h` and leaves both terms invariant on paper, but it evaluates a
/// different expression, and the two agree only to rounding: measured 1 ULP
/// apart at `smooth_min(-0.75, 0.4, 0.1)`. That was written here as an
/// unqualified "commutative" until a test disagreed.
///
/// A `k` of zero degenerates to an ordinary `min` rather than dividing by zero.
#[must_use]
pub fn smooth_min<R: Real>(a: R, b: R, k: R) -> R {
    if k <= R::ZERO {
        return a.min(b);
    }
    let h = (R::HALF + R::HALF * (b - a) / k).clamp(R::ZERO, R::ONE);
    (b + (a - b) * h) - k * h * (R::ONE - h)
}

/// A base field with a stack of brushes applied in order.
///
/// The order is part of the value: see the module docs for exactly when it can
/// be permuted and when it cannot.
#[derive(Clone, Debug)]
pub struct BrushStack<'a, F, S> {
    /// The field the brushes are applied to.
    pub base: F,
    /// The edits, applied first to last.
    pub brushes: &'a [Brush<S>],
}

impl<F, S> Sdf for BrushStack<'_, F, S>
where
    F: Sdf,
    S: Sdf<Scalar = F::Scalar>,
{
    type Scalar = F::Scalar;

    fn sample(&self, p: [Self::Scalar; 3]) -> Self::Scalar {
        let mut value = self.base.sample(p);
        for brush in self.brushes {
            value = apply(brush.op, value, brush.shape.sample(p));
        }
        value
    }
}

#[cfg(test)]
mod tests;
