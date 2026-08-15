//! Speculative algorithms, off by default and exempt from semver.
//!
//! Ticket: X-003. Enable with `features = ["experimental"]`.
//!
//! # What this module is for
//!
//! Phases 11 through 15 are experiments: the same algorithm with one rule
//! swapped, measured against the rule it replaces. X-002 built the seam that
//! makes that possible without a second execution path; this is where the other
//! side of the seam lives, so a speculative rule does not have to enter the
//! stable public API to be measured.
//!
//! # Exempt from semver, and from nothing else
//!
//! Anything here may change or disappear in a patch release. It is **not**
//! exempt from correctness: T-001's validity suite, T-004's determinism check
//! and the property tests apply exactly as they do to the stable extractors, and
//! an experimental rule that cannot pass them is a wrong answer rather than an
//! early one.
//!
//! It is also not exempt from the dependency rule. The feature adds no
//! dependencies — `cargo tree -p isomesh -e normal` is `isomesh + libm` with the
//! feature on or off, which `the_experimental_feature_adds_no_dependencies`
//! checks rather than assumes.

use crate::Sdf;
use crate::dual::{CellVertices, VertexRule};
use crate::dual_contouring::{Clamp, apply_clamp, solve};
use crate::hermite::HermiteCell;
use crate::real::Real;

/// Trettner & Kobbelt's probabilistic plane quadric, which for this crate's
/// formulation is the existing solve with a **crossing-count-scaled**
/// regularizer.
///
/// Ticket: X-004. Source: *Fast and Robust QEF Minimization using Probabilistic
/// Quadrics*, Trettner & Kobbelt, CGF 39(2) (`10.1111/cgf.13933`), §3.1.
///
/// # The paper's rule reduces to ours, and the derivation is the finding
///
/// X-004 was written expecting this to *supersede* the Tikhonov regularizer.
/// It does not — it **is** that regularizer, with one number changed, and the
/// reason is that this crate already solves in centroid-relative coordinates
/// (M-238).
///
/// Equations (6) and (7) give the probabilistic plane quadric under Gaussian
/// noise `Σₙ` on the normal:
///
/// ```text
/// A = Σ nᵢnᵢᵀ + N·Σₙ        b = Σ nᵢnᵢᵀqᵢ + Σₙ·Σ qᵢ
/// ```
///
/// Write `x = c + Δ` and `qᵢ = c + rᵢ` for the crossings' centroid `c`, take
/// `Σₙ = σ²I`, and the system `Ax = b` becomes
///
/// ```text
/// (M + Nσ²I) Δ = Σ nᵢdᵢ + σ² Σ rᵢ
/// ```
///
/// where `dᵢ = nᵢ·(qᵢ − c)`. **`Σ rᵢ` is identically zero**, because `c` is
/// defined as the arithmetic mean of the `qᵢ` — so the extra term the paper adds
/// to `b` vanishes, and what remains is exactly
/// [`solve_with`](crate::dual_contouring::solve::solve_with) at `λ = Nσ²`.
///
/// `the_probabilistic_quadric_is_the_existing_solve` pins that as a numeric
/// identity rather than leaving it as algebra.
///
/// # So what is actually different
///
/// The **scaling**. [`Qef`](crate::dual_contouring::Qef) applies one fixed `λ`
/// to every cell; this applies `N·σ²`, so a cell with twelve crossings is
/// regularized four times as hard as one with three. That is the whole
/// experimental content of this rule, and it is measurable — which is why it
/// exists as a rule rather than as a paragraph.
///
/// # What is *not* implemented, and why that is the honest boundary
///
/// The paper's actual novelty is **anisotropic** `Σₙ`: a per-plane normal
/// covariance, which is what lets it distinguish noise from a sharp feature.
/// That needs a noise model, and this crate's fields are analytic with exact
/// gradients — there is no measurement error to model. Inventing a covariance so
/// the formula has somewhere to put one would be fitting the data to the method.
/// A consumer meshing scanned or quantised data has such a model and would be
/// the reason to add it.
#[derive(Clone, Copy, Debug)]
pub struct ProbabilisticQuadric {
    /// Whether the result is confined to its cell. Same meaning as
    /// [`Qef::clamp`](crate::dual_contouring::Qef::clamp).
    pub clamp: Clamp,
    /// The isotropic normal variance `σ²`, in units of the regularizer.
    ///
    /// The effective regularizer is `N·σ²` for a cell with `N` crossings, so
    /// this is *per crossing* where [`Qef`](crate::dual_contouring::Qef)'s λ is
    /// per cell. [`DEFAULT_SIGMA_SQUARED`] is chosen so a typical cell lands
    /// near the λ the stable path uses.
    pub sigma_squared: f64,
}

/// The default per-crossing variance.
///
/// Chosen so that a cell with the **median** crossing count reproduces the
/// stable path's regularizer rather than by fitting anything: the reference
/// fields' dual cells carry four crossings far more often than any other count,
/// so `σ² = LAMBDA / 4` makes the common cell agree exactly and lets the tails
/// differ, which is the effect being measured.
pub const DEFAULT_SIGMA_SQUARED: f64 = solve::LAMBDA / 4.0;

impl Default for ProbabilisticQuadric {
    fn default() -> Self {
        Self {
            clamp: Clamp::ToCell,
            sigma_squared: DEFAULT_SIGMA_SQUARED,
        }
    }
}

impl<R: Real> VertexRule<R> for ProbabilisticQuadric {
    fn place<S: Sdf<Scalar = R>>(
        &self,
        sdf: &S,
        corner: &[R; 8],
        base: [u32; 3],
        origin: [R; 3],
        cell_size: R,
        out: &mut CellVertices<R>,
    ) {
        let cell_origin = [
            origin[0] + cell_size * R::from_f64(f64::from(base[0])),
            origin[1] + cell_size * R::from_f64(f64::from(base[1])),
            origin[2] + cell_size * R::from_f64(f64::from(base[2])),
        ];
        let cell = HermiteCell::from_corners(sdf, corner, cell_origin, cell_size);
        // The one line that differs from `Qef`: the regularizer scales with the
        // number of planes, per the paper's `N·Σₙ`.
        let lambda = R::from_f64(cell.len() as f64 * self.sigma_squared);
        let Some(x) = solve::solve_with(&cell, lambda) else {
            return;
        };
        out.push_whole_cell(apply_clamp(self.clamp, x, cell_origin, cell_size));
    }
}

#[cfg(test)]
mod tests;
