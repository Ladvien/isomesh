//! §4.3.2 — every zero along an edge, not just the first.
//!
//! This is the half of subgrid Marching Tetrahedra that
//! [`surface`](super::surface) is defined against. Baktash, Gillespie & Crane,
//! `10.48550/arXiv.2606.00454` §1:
//!
//! > We replace 0-dimensional sampling (evaluate `f` at each grid node), with
//! > 1-dimensional root finding (find all zeros of `f` along each grid edge).
//!
//! A sign test asks one question per edge and gets one bit back. This asks *how
//! many times* and *where*, which is what lets a feature thinner than a cell
//! survive — M-67 measured the cost of the other approach: a sign test cannot
//! distinguish **95.6%** of the configurations a tetrahedron can be in.
//!
//! # Sampling is the primary path, not a fallback
//!
//! The 1D search below samples and refines rather than solving analytically, and
//! that is a deliberate single path rather than a degraded one. §1.3 prices it:
//!
//! > 1D marching can of course miss intersections, \[but\] we are no worse off
//! > than classic marching.
//!
//! An analytic root solver would need each field to expose its polynomial form,
//! which [`Sdf`] deliberately does not — the trait's whole point is
//! that a field is a black box. Offering both would mean two root finders and
//! two sets of results for one input, which is exactly what `CLAUDE.md`'s
//! one-path rule exists to prevent. **What sampling loses is bounded and
//! stated**: a pair of roots closer together than the sample spacing is invisible,
//! exactly as a pair closer than the *grid* spacing is invisible to classic
//! marching.
//!
//! # The sign convention, and why zero is not a root
//!
//! Inside is `f < 0`, and a sample of exactly zero counts as **outside** — the
//! same convention as [`marching_cubes`](crate::marching_cubes) and
//! [`marching_tetrahedra`](crate::marching_tetrahedra), stated here because a
//! root finder is precisely where an inconsistency would show up as a mesh that
//! is subtly not the field's. A crossing is a change in `f < 0` between two
//! samples, so a field that merely touches zero without passing through it
//! contributes no root, and neither does one that sits at zero along a stretch
//! of edge.

use alloc::vec::Vec;

use crate::real::Real;
use crate::sdf::Sdf;

/// How many refinement steps a bracketed root gets before it is accepted.
///
/// The bisection stops early once the midpoint stops differing from an endpoint,
/// which is machine precision for whichever [`Real`] is in use — 24 bits of
/// mantissa for `f32`, 53 for `f64`. This is the backstop that guarantees
/// termination regardless, and 80 is comfortably past `f64`'s worst case of 53.
const REFINEMENTS: u32 = 80;

/// Every crossing along one edge, as parameters in `[0, 1]` from `from`.
///
/// `samples` is the 1D marching resolution: the edge is divided into that many
/// intervals, and each interval carrying a sign change yields one root. It must
/// be at least 1.
///
/// Results are appended to `out` in **ascending order**, which is
/// [`TetCrossings`](super::surface::TetCrossings)' contract. Two roots can never
/// coincide, because each comes from a distinct bracketing interval.
///
/// # Determinism
///
/// The parameters are computed as `i / samples` and refined by bisection on
/// those, so a given `(from, to, sdf, samples)` always produces bit-identical
/// output. That matters more here than anywhere else in the crate: two
/// tetrahedra sharing an edge call this independently and must agree exactly, or
/// the conformity §3.1 guarantees combinatorially is lost geometrically.
///
/// **The caller is responsible for calling this with identical endpoints** from
/// both sides. M-32 is the precedent: `origin + h·i` and `(origin + h·c·n) + h·i`
/// are equal by algebra and not by IEEE.
pub fn all_roots<R: Real, S: Sdf<Scalar = R>>(
    from: [R; 3],
    to: [R; 3],
    sdf: &S,
    samples: u32,
    out: &mut Vec<R>,
) {
    if samples == 0 {
        return;
    }

    let at = |t: R| -> R {
        sdf.sample([
            from[0] + (to[0] - from[0]) * t,
            from[1] + (to[1] - from[1]) * t,
            from[2] + (to[2] - from[2]) * t,
        ])
    };
    // Inside is `f < 0`; a sample of exactly zero is outside.
    let inside = |v: R| v < R::ZERO;

    let step = R::ONE / R::from_f64(f64::from(samples));
    let mut low = R::ZERO;
    let mut low_inside = inside(at(low));

    for i in 1..=samples {
        // `i * step` rather than an accumulating add: the running sum would
        // drift, and the last value would not be exactly 1.
        let high = if i == samples {
            R::ONE
        } else {
            R::from_f64(f64::from(i)) * step
        };
        let high_inside = inside(at(high));

        if high_inside != low_inside {
            out.push(refine(low, high, low_inside, &at, &inside));
        }

        low = high;
        low_inside = high_inside;
    }
}

/// Bisect a bracketed sign change down to machine precision.
///
/// The bracket always contains a sign change on entry, and every step preserves
/// that, so the returned parameter is strictly inside the original interval —
/// which is what keeps roots from different intervals from colliding.
fn refine<R: Real>(
    mut low: R,
    mut high: R,
    low_inside: bool,
    at: &impl Fn(R) -> R,
    inside: &impl Fn(R) -> bool,
) -> R {
    let half = R::from_f64(0.5);
    for _ in 0..REFINEMENTS {
        let mid = (low + high) * half;
        // Once the midpoint is one of the endpoints, the bracket is as tight as
        // this `Real` can express and further steps would not move.
        if mid <= low || mid >= high {
            break;
        }
        if inside(at(mid)) == low_inside {
            low = mid;
        } else {
            high = mid;
        }
    }
    // Either end is within one ulp; the upper one is chosen so the result is
    // strictly greater than the interval's start, keeping the ascending-and-
    // distinct contract even when a root sits exactly on a sample.
    high
}

#[cfg(test)]
mod tests;
