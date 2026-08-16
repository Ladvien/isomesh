//! Does a field meet the bound it declares?
//!
//! Ticket: F-002. The other half of F-001, and the half that makes the
//! declaration mean something.
//!
//! # Why a declaration needs a checker
//!
//! F-001 replaced `is_exact_distance() -> bool` because `csg_difference` had
//! declared `true` with `// away from the seam` beside it for months and nothing
//! looked. Replacing the *type* does not fix that: a `FieldBound` nobody checks
//! is the same defect with more cases. It is not a hypothetical either — the
//! first draft of F-001 declared `noise_cavity` at `l = 2.598` when its gradient
//! reaches **7.73**, and only a checker caught it (M-244).
//!
//! # What is measured, and what the numbers mean
//!
//! `‖∇f‖`, sampled over the field's own domain.
//!
//! - **`sup`** is the largest gradient seen. Against
//!   [`Lipschitz`](crate::fields::FieldBound::Lipschitz) this is the number that
//!   must not exceed `l`.
//! - **`eikonal_fraction`** is the share of samples with `‖∇f‖ ≈ 1`, which is the
//!   differential form of *"the value is the distance"*. An
//!   [`Exact`](crate::fields::FieldBound::Exact) field should be near 1.0 and a
//!   gyroid nowhere near it.
//!
//! # The direction of the error is the whole point
//!
//! **A sampled maximum is a lower bound on a supremum.** So this can prove a
//! declaration *wrong* and can never prove one right: finding `‖∇f‖ = 7.73`
//! against a declared `2.598` is conclusive, while finding nothing above `1.7`
//! against a declared `3.46` says only that the sampling missed. That asymmetry
//! is why [`FieldBoundReport::violates`] answers a one-sided question, and why
//! `noise_cavity` is declared [`Unbounded`](crate::fields::FieldBound::Unbounded)
//! rather than given the measured figure — declaring a sampled maximum as a
//! Lipschitz constant is unsound in exactly the direction that lets a sphere
//! tracer step through a surface.

#[cfg(test)]
mod tests;

use crate::fields::{FieldBound, ReferenceField};
use crate::real::Real;

/// How close to `1` a gradient must be to count as eikonal.
///
/// Loose on purpose. An exact field's gradient is one everywhere it is defined,
/// but a finite sample lands near creases — a box's edges, a CSG seam — where it
/// is not defined at all and a central difference reads short. The fraction is a
/// description of the field, not a gate, so the band only has to be wide enough
/// that a crease does not dominate it.
pub const EIKONAL_TOLERANCE: f64 = 0.05;

/// What sampling `‖∇f‖` over a field's domain found.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldBoundReport {
    /// The bound the field declares.
    pub declared: FieldBound,
    /// Largest `‖∇f‖` seen. A **lower** bound on the true supremum.
    pub sup: f64,
    /// Smallest `‖∇f‖` seen.
    pub inf: f64,
    /// Fraction of samples with `‖∇f‖` within [`EIKONAL_TOLERANCE`] of one.
    pub eikonal_fraction: f64,
    /// How many samples produced a finite gradient.
    pub samples: usize,
}

impl FieldBoundReport {
    /// The declared bound is **provably** wrong, given what was sampled.
    ///
    /// One-sided, per this module's header: a sampled maximum can exceed a
    /// declared constant and settle the question, but never falling below one
    /// settles nothing.
    ///
    /// - [`Lipschitz`](FieldBound::Lipschitz): violated when `sup > l`.
    /// - [`Exact`](FieldBound::Exact): violated when `sup` exceeds one, since
    ///   Corollary 1 of Bálint, Valasek & Gergó makes `1` the *smallest*
    ///   Lipschitz constant of a true distance — an exact field cannot exceed it.
    /// - [`Underestimate`](FieldBound::Underestimate): the same, because this
    ///   crate's underestimates are built from exact operands by `min`/`max`,
    ///   which preserve the constant while destroying exactness.
    /// - [`Unbounded`](FieldBound::Unbounded): never violated. It claims nothing.
    #[must_use]
    pub fn violates(&self, slack: f64) -> bool {
        match self.declared {
            FieldBound::Unbounded => false,
            FieldBound::Lipschitz { l } => self.sup > l * (1.0 + slack),
            FieldBound::Exact | FieldBound::Underestimate { .. } => self.sup > 1.0 * (1.0 + slack),
        }
    }
}

/// Sample `‖∇f‖` over `field`'s domain and report what was found.
///
/// `n` is samples per axis, so cost is `n³` gradient evaluations.
///
/// # The grid is offset, deliberately
///
/// Samples land at cell *centres* plus an irrational-ish nudge. A grid aligned
/// to the domain lands exactly on a box's faces and on the noise lattice, which
/// are precisely the points where the gradient is undefined or degenerate — so
/// an aligned sweep measures the discontinuities instead of the field, and its
/// `inf` is an artefact.
#[must_use]
pub fn field_bound_report<F>(field: &F, n: u32) -> FieldBoundReport
where
    F: ReferenceField,
    F::Scalar: Real,
{
    let (lo, hi) = field.domain();
    let mut sup = 0.0f64;
    let mut inf = f64::INFINITY;
    let mut eikonal = 0usize;
    let mut samples = 0usize;

    // Chosen to avoid landing on lattice points and face planes; the exact value
    // does not matter, only that it is not a simple fraction of the domain.
    let nudge = 0.031_7_f64;

    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let at = |v: u32, a: F::Scalar, b: F::Scalar| -> F::Scalar {
                    let t = (f64::from(v) + 0.5 + nudge) / f64::from(n);
                    a + (b - a) * F::Scalar::from_f64(t)
                };
                let p = [
                    at(i, lo[0], hi[0]),
                    at(j, lo[1], hi[1]),
                    at(k, lo[2], hi[2]),
                ];
                let g = field.gradient(p);
                let len2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
                let len = len2.sqrt().as_f64();
                if !len.is_finite() {
                    continue;
                }
                samples += 1;
                sup = sup.max(len);
                inf = inf.min(len);
                if (len - 1.0).abs() <= EIKONAL_TOLERANCE {
                    eikonal += 1;
                }
            }
        }
    }

    FieldBoundReport {
        declared: field.bound(),
        sup,
        inf: if samples == 0 { 0.0 } else { inf },
        eikonal_fraction: if samples == 0 {
            0.0
        } else {
            eikonal as f64 / samples as f64
        },
        samples,
    }
}
