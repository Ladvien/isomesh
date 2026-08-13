//! Surface attributes that survive destruction: paint that stays where it was
//! sprayed when the wall under it is carved away.
//!
//! Row 4 of `docs/research/2026-08-11-novel-gameplay-opportunities.md` — *"you
//! spray graffiti on a wall, then blow a hole through it, and the paint on the
//! remaining wall is still exactly where you sprayed it — not smeared, not
//! reset."*
//!
//! # Why there is no attribute *transfer* here
//!
//! The research prices that row on **common subdivision plus L²-nearest
//! attribute transfer** (Integer Coordinates §3.8, §4.4), on the reasoning that
//! the shared tetrahedral grid is the common refinement so the expensive half
//! comes free. That machinery exists to move attributes from an old mesh to a
//! new one, and it is needed exactly when the attributes live **on the mesh**.
//!
//! They do not have to. This crate's world is not stored voxels; it is a base
//! field plus an ordered log of edits (the [`brush`] module), which is what
//! lets [`chunk::dirty`](crate::chunk::dirty) express an edit as *two fields*
//! and what lets an editor treat undo as a re-fold of the log rather than a
//! snapshot. Put the paint in the same log and it is a function of world
//! position, so re-meshing cannot move it: the carve changes the *surface*, and
//! the paint is not on the surface, it is in the space the surface passes
//! through.
//!
//! So the transfer is not cheap here. It is **unnecessary**, and the result is
//! exact rather than L²-nearest. See M-137.
//!
//! # One log, both kinds of edit
//!
//! Sprays and carves interleave in a single [`Edit`] list, and that ordering is
//! load-bearing rather than tidy. Paint is confined to a shell around the
//! surface **as it stood when the spray happened** — that is what makes it
//! graffiti rather than dyed material, and it is what leaves the inside of a
//! fresh hole bare. Walking one log in order means the field value the walk is
//! already carrying *is* that surface, at no extra cost and with no index into
//! another list to get wrong.
//!
//! ```
//! use isomesh::brush::Brush;
//! use isomesh::fields::{BoxExact, Sphere};
//! use isomesh::paint::{Edit, PaintStack, Splat};
//! use isomesh::Sdf;
//!
//! // A wall.
//! let wall = BoxExact::<f64> { center: [0.0; 3], half_extents: [2.0, 2.0, 0.25] };
//!
//! let log = [
//!     // Spray red on the front face...
//!     Edit::Spray(Splat {
//!         shape: Sphere { center: [0.0, 0.0, -0.25], radius: 0.75 },
//!         color: [1.0, 0.0, 0.0, 1.0],
//!         softness: 0.1,
//!         depth: 0.05,
//!     }),
//!     // ...then blow a hole through it, well away from the paint.
//!     Edit::Carve(Brush::subtract(Sphere { center: [1.2, 0.0, 0.0], radius: 0.5 })),
//! ];
//!
//! let world = PaintStack { base: wall, edits: &log, background: [0.5, 0.5, 0.5, 1.0] };
//!
//! // Still red where it was sprayed, after the carve.
//! assert!(world.color_at([0.0, 0.0, -0.25])[0] > 0.9);
//! // And the field is a field: extractors take this directly.
//! assert!(world.sample([1.2, 0.0, 0.0]) > 0.0);
//! ```

use alloc::vec::Vec;

use crate::brush::{self, Brush};
use crate::{Real, Sdf};

/// One spray: a shape, a colour, and how far the paint reaches.
///
/// The shape is any [`Sdf`], so a spray can be a sphere, a capsule dragged
/// along a stroke, or a stencil — this type does not care.
///
/// # The two widths, which do different jobs
///
/// [`softness`](Self::softness) is measured **outward from the shape**, in
/// world units: coverage is full where the shape is solid and fades to nothing
/// `softness` beyond its boundary. It is the edge of the spray cone.
///
/// [`depth`](Self::depth) is measured **from the surface**, and it is the one
/// that makes this paint rather than dye. Coverage is full at the surface as it
/// stood when this splat was sprayed and falls to nothing `depth` away from it,
/// on both sides. Carve into painted material and the newly exposed interior is
/// bare, because it was never within `depth` of the old surface.
///
/// # The shape decides reach, and it is the *only* thing that does
///
/// There is no spray direction here, and its absence is a decision rather than
/// an omission. The shape is an arbitrary [`Sdf`], so it already carries every
/// bit of directionality a caller might want — a cone, a capsule swept along a
/// stroke, a stencil intersected with a half-space — and a `direction` field
/// would duplicate what the shape already expresses.
///
/// The consequence to know about: **a shape that passes through a thin wall
/// paints the far face too**, because both faces are surface and both are
/// inside the shape. A spray that must not reach through gets a shape that does
/// not reach through. There is a test for both halves of that.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Splat<P, R: Real> {
    /// Where the spray lands. Negative inside, like every field here.
    pub shape: P,
    /// Straight (non-premultiplied) RGBA.
    ///
    /// The alpha channel is **coverage**: it scales how much of this colour
    /// reaches the accumulator. It is not written to the output — see
    /// [`PaintStack::color_at`].
    pub color: [R; 4],
    /// Falloff width outside the shape, in world units. Zero is a hard edge.
    pub softness: R,
    /// How far from the spray-time surface the paint reaches, in world units.
    pub depth: R,
}

/// One entry in the world's edit log: a carve or a spray.
///
/// The order is part of the value. [`Brush`] entries already have measured
/// reordering rules — see [`BrushOp::commutes_with`](crate::brush::BrushOp::commutes_with)
/// and M-36..M-38 — and [`Spray`](Self::Spray) entries are stricter still: a
/// spray reads the surface as it stands *at its position in the log*, so moving
/// one across a carve changes what it paints. None of them commute.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Edit<S, P, R: Real> {
    /// Change the shape of the solid.
    Carve(Brush<S>),
    /// Change the colour of the surface.
    Spray(Splat<P, R>),
}

/// A base field with an ordered log of carves and sprays applied to it.
///
/// This is an [`Sdf`], so it goes straight into any extractor in the crate;
/// [`color_at`](Self::color_at) is the second output, evaluated at whatever
/// positions the extractor produced.
///
/// # Relationship to [`BrushStack`](crate::brush::BrushStack)
///
/// `BrushStack` is the unpainted path and stays the right choice when there are
/// no attributes in play. A `PaintStack` whose log contains no
/// [`Spray`](Edit::Spray) samples **bit-identically** to the `BrushStack` over
/// the same brushes in the same order — both fold through [`brush::apply`], so
/// they cannot drift, and there is a test that says so.
#[derive(Clone, Debug)]
pub struct PaintStack<'a, F, S, P>
where
    F: Sdf,
    S: Sdf<Scalar = F::Scalar>,
    P: Sdf<Scalar = F::Scalar>,
{
    /// The field the log is applied to.
    pub base: F,
    /// The edits, applied first to last.
    pub edits: &'a [Edit<S, P, F::Scalar>],
    /// The colour of unpainted surface.
    pub background: [F::Scalar; 4],
}

/// A linear ramp from 1 at `x <= 0` to 0 at `x >= width`.
///
/// A zero or negative `width` is the limit of that ramp — a hard step — rather
/// than a division by zero, on the same reasoning as
/// [`smooth_min`](crate::brush::smooth_min) treating `k <= 0` as an ordinary
/// `min`: the degenerate case has a correct answer, so it gets that answer
/// instead of a rejection.
fn ramp<R: Real>(x: R, width: R) -> R {
    if width <= R::ZERO {
        if x <= R::ZERO { R::ONE } else { R::ZERO }
    } else {
        (R::ONE - x / width).clamp(R::ZERO, R::ONE)
    }
}

impl<F, S, P> PaintStack<'_, F, S, P>
where
    F: Sdf,
    S: Sdf<Scalar = F::Scalar>,
    P: Sdf<Scalar = F::Scalar>,
{
    /// The colour at `p`, from one walk of the log.
    ///
    /// Splats composite in log order, later over earlier. Each contributes
    ///
    /// ```text
    /// coverage = ramp(shape(p), softness) · ramp(|f(p)|, depth) · color.a
    /// ```
    ///
    /// where `f` is the field **as it stood at that point in the log** — which
    /// is the value this walk is already carrying, and the reason paint does
    /// not chase a later carve.
    ///
    /// # Alpha
    ///
    /// The returned alpha is [`background`](Self::background)'s, unchanged.
    /// Each splat's alpha is consumed as coverage in the expression above; a
    /// surface does not become more or less transparent by being painted.
    pub fn color_at(&self, p: [F::Scalar; 3]) -> [F::Scalar; 4] {
        let mut field = self.base.sample(p);
        let mut color = self.background;
        for edit in self.edits {
            match edit {
                Edit::Carve(b) => field = brush::apply(b.op, field, b.shape.sample(p)),
                Edit::Spray(s) => {
                    let coverage = ramp(s.shape.sample(p), s.softness)
                        * ramp(field.abs(), s.depth)
                        * s.color[3];
                    // RGB only: alpha is coverage, consumed above, not written.
                    for (channel, &target) in color.iter_mut().zip(&s.color).take(3) {
                        *channel = *channel + (target - *channel) * coverage;
                    }
                }
            }
        }
        color
    }
}

impl<F, S, P> Sdf for PaintStack<'_, F, S, P>
where
    F: Sdf,
    S: Sdf<Scalar = F::Scalar>,
    P: Sdf<Scalar = F::Scalar>,
{
    type Scalar = F::Scalar;

    /// The field, ignoring paint.
    ///
    /// Sprays are skipped without evaluating their shapes, so meshing a painted
    /// world costs one `match` arm more than meshing an unpainted one and not a
    /// field evaluation more.
    fn sample(&self, p: [Self::Scalar; 3]) -> Self::Scalar {
        let mut value = self.base.sample(p);
        for edit in self.edits {
            if let Edit::Carve(b) = edit {
                value = brush::apply(b.op, value, b.shape.sample(p));
            }
        }
        value
    }
}

/// Colour every position, into a caller-owned buffer.
///
/// The vertex-attribute counterpart to
/// [`normals::recompute`](crate::normals::recompute): hand it the positions an
/// extractor produced and it fills `out` with one RGBA per vertex, parallel to
/// them.
///
/// `out` is truncated and refilled without releasing its capacity, so the same
/// buffer serves every chunk of a re-mesh — rule 6, the same contract as
/// [`MeshBuffer::reset`](crate::MeshBuffer::reset).
pub fn shade<F, S, P>(
    positions: &[[F::Scalar; 3]],
    stack: &PaintStack<'_, F, S, P>,
    out: &mut Vec<[F::Scalar; 4]>,
) where
    F: Sdf,
    S: Sdf<Scalar = F::Scalar>,
    P: Sdf<Scalar = F::Scalar>,
{
    out.clear();
    out.reserve(positions.len());
    out.extend(positions.iter().map(|&p| stack.color_at(p)));
}

#[cfg(test)]
mod tests;
