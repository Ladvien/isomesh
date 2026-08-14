//! MC33's interior rule: Chernyaev's sweeping plane, with Custodio's correction.
//!
//! Chernyaev, E. V., *Marching Cubes 33: Construction of Topologically Correct
//! Isosurfaces*, Technical Report CN/95-17, Institute for High Energy Physics
//! (1995), and Custodio, L., Etiene, T., Pesco, S. and Silva, C., *Practical
//! considerations on Marching Cubes 33 topological correctness*, Computers &
//! Graphics 37(7), pp. 840–850 (`10.1016/j.cag.2013.04.004`), §3.1 and §5.1.
//!
//! # What it decides, and how it differs from the face rule
//!
//! [`super::ambiguity`] resolves an ambiguous *face* by evaluating the bilinear
//! interpolant at that face's saddle point. This resolves an ambiguous *cell*:
//! two opposite faces can each be ambiguous and separately resolved, and the
//! cell's own topology still not be determined, because the regions they carry
//! may or may not be joined through the cell's **interior**. Custodio §3.1
//! states the criterion:
//!
//! > if there is a plane cutting the cube such that its saddle point is
//! > positive, it means that there is a positive area crossing the cube, i.e.
//! > the positive vertices are connected inside the cube.
//!
//! So sweep a plane between the two faces. At height `t` the trilinear
//! interpolant restricts to a bilinear function whose four corner values are
//! linear in `t`, and that function has one saddle, whose value is
//!
//! ```text
//! f(x_c(t)) = (A_t·C_t − B_t·D_t) / (A_t + C_t − B_t − D_t) = F(t) / Δ(t)
//! ```
//!
//! # The correction, which is the whole point of this module
//!
//! `F` is a **quadratic** in `t` and `Δ` is **linear**, and Chernyaev's test
//! tracks the sign of `F` alone. Custodio §5.1:
//!
//! > the polynomial F(t) … used by Chernyaev's MC33 algorithm for tracking the
//! > sign of the saddle point, is a second order equation in t and thus can only
//! > allow for two sign changes. Therefore, the sign tracked by the MC33
//! > algorithm will not match the expected one at some point.
//!
//! The quantity that actually matters is the **quotient**, and a quotient of a
//! quadratic by a linear is a hyperbola: where `Δ` has a root the saddle value
//! has a pole, and the sign jumps across it. That admits **three** sign changes
//! on `(0, 1)`, which no quadratic can track. Their Figure 6 is a case 13.5.2
//! that Chernyaev's test reads as 13.5.1 for exactly this reason.
//!
//! This module therefore evaluates the criterion **exactly**, rather than
//! through any polynomial proxy: `F/Δ` has constant sign on each interval
//! between its own roots and poles, so the breakpoints are enumerated and one
//! interior point of each resulting subinterval is tested. There is no
//! tolerance, no sampling and no iteration count — the answer is a finite sign
//! computation.
//!
//! # The coefficients, derived here rather than transcribed
//!
//! Writing `α = A₁ − A₀` and so on, `A_t = A₀ + tα`, so
//!
//! ```text
//! A_t·C_t = A₀C₀ + t(A₀γ + C₀α) + t²αγ
//! B_t·D_t = B₀D₀ + t(B₀δ + D₀β) + t²βδ
//! ```
//!
//! and subtracting gives `F(t) = a t² + b t + c` with
//!
//! ```text
//! a = αγ − βδ
//! b = C₀α + A₀γ − D₀β − B₀δ
//! c = A₀C₀ − B₀D₀
//! ```
//!
//! Those agree term for term with the three the paper prints, which is why they
//! are written out above: the agreement is the check (V-25). Rule 5 forbids
//! guessing a published formula, and re-deriving one is how you avoid having to.
//!
//! # What is *not* here
//!
//! The triangulation. Knowing that a cell has a tunnel is not the same as
//! meshing one, and MC33's tunnel cases need vertices in the cell **interior**
//! that this crate's grid-edge-keyed vertex cache has no slot for. That is
//! A-002b, and it is a separate ticket for that reason.

#[cfg(test)]
mod tests;

use crate::real::Real;
use crate::{Error, Result};

/// The two faces a sweep runs between, as the bilinear corner values on each.
///
/// Each array is `[A, B, C, D]` in the cyclic order Custodio's Figure 2 uses, so
/// that `A` and `C` are one diagonal and `B` and `D` are the other. `lo` is the
/// face at `t = 0` and `hi` the face at `t = 1`; the corners must correspond,
/// meaning `lo[k]` and `hi[k]` are the two ends of a cell edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SweptFaces<R> {
    lo: [R; 4],
    hi: [R; 4],
}

/// What the sweep found.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Interior {
    /// Some cutting plane has a positive saddle: the same-signed regions the two
    /// faces carry are connected through the cell's interior.
    Joined,
    /// No cutting plane has a positive saddle: they are separated.
    Separated,
}

impl<R: Real> SweptFaces<R> {
    /// The two faces of a sweep.
    ///
    /// # Errors
    ///
    /// [`Error::DegenerateSweep`] if either face's bilinear denominator
    /// `A + C − B − D` is zero, which means that face's bilinear function has no
    /// saddle and the criterion this module evaluates is undefined there.
    ///
    /// **This is a precondition rather than a case to handle.** On an ambiguous
    /// face one diagonal is strictly negative and the other non-negative — the
    /// same argument [`super::ambiguity`] makes — so the denominator is a
    /// strictly negative sum minus a non-negative one, or the reverse, and
    /// cannot be zero. A caller reaching this error has applied the interior
    /// test to a face that is not ambiguous, and is asking a question that has
    /// no answer rather than one whose answer is hard.
    pub fn new(lo: [R; 4], hi: [R; 4]) -> Result<Self> {
        let faces = Self { lo, hi };
        if faces.denominator(R::ZERO) == R::ZERO || faces.denominator(R::ONE) == R::ZERO {
            return Err(Error::DegenerateSweep);
        }
        Ok(faces)
    }

    /// `Δ(t) = A_t + C_t − B_t − D_t`, the bilinear denominator at height `t`.
    ///
    /// Linear in `t`, and its root is the pole Chernyaev's quadratic cannot see.
    #[must_use]
    pub fn denominator(&self, t: R) -> R {
        let at = |v: &[R; 4]| v[0] + v[2] - v[1] - v[3];
        at(&self.lo) + (at(&self.hi) - at(&self.lo)) * t
    }

    /// `F(t) = A_t·C_t − B_t·D_t`, the bilinear numerator at height `t`.
    ///
    /// Evaluated from the interpolated corner values rather than from the
    /// expanded `a t² + b t + c`, so that the definition and the implementation
    /// are the same expression. The expansion is derived in this module's docs
    /// because the *degree* is what the correction turns on, not because the
    /// coefficients are needed to evaluate it.
    #[must_use]
    pub fn numerator(&self, t: R) -> R {
        let s = R::ONE - t;
        let v = |k: usize| self.lo[k] * s + self.hi[k] * t;
        v(0) * v(2) - v(1) * v(3)
    }

    /// The saddle value `f(x_c(t))` of the cutting plane at height `t`.
    ///
    /// Undefined at [`Self::pole`], which is why the decision below never
    /// evaluates it there.
    #[must_use]
    pub fn saddle(&self, t: R) -> R {
        self.numerator(t) / self.denominator(t)
    }

    /// Where the cutting plane's saddle *sits*, in the plane's own `[0, 1]²`
    /// coordinates.
    ///
    /// Custodio's Equation (1), with `A` at `(0, 0)`, `B` at `(1, 0)`, `C` at
    /// `(1, 1)` and `D` at `(0, 1)` — so `A`/`C` are one diagonal and `B`/`D`
    /// the other.
    ///
    /// **This is what makes the correction visible rather than merely provable.**
    /// Their Figure 4 plots this path and calls it hyperbolic: it is a linear
    /// function over `Δ(t)`, so as the sweep approaches [`pole`](Self::pole) the
    /// saddle runs off to infinity and returns from the other side. A point that
    /// leaves `[0, 1]²` has left the face, which is why the value there can
    /// change sign without the numerator doing anything.
    #[must_use]
    pub fn saddle_position(&self, t: R) -> [R; 2] {
        let s = R::ONE - t;
        let v = |k: usize| self.lo[k] * s + self.hi[k] * t;
        let d = self.denominator(t);
        [(v(0) - v(3)) / d, (v(0) - v(1)) / d]
    }

    /// Where the saddle value has its pole, if that is inside the sweep.
    ///
    /// `Δ` is linear and non-zero at both ends, so it has a root in `(0, 1)`
    /// exactly when its two end values differ in sign. **This is the term
    /// Chernyaev's test drops**, and Custodio's Figure 6 counterexample is a
    /// configuration where it lies inside the sweep.
    #[must_use]
    pub fn pole(&self) -> Option<R> {
        let (d0, d1) = (self.denominator(R::ZERO), self.denominator(R::ONE));
        if (d0 < R::ZERO) == (d1 < R::ZERO) {
            return None;
        }
        // d0 − d1 cannot be zero: the two have strictly opposite signs.
        let t = d0 / (d0 - d1);
        if t > R::ZERO && t < R::ONE {
            Some(t)
        } else {
            None
        }
    }

    /// Custodio's corrected interior test.
    ///
    /// Answers the criterion exactly. `F/Δ` is continuous and non-zero on each
    /// open interval between consecutive breakpoints — the real roots of `F` and
    /// the pole of `Δ` — so its sign there is the sign at any interior point, and
    /// testing one midpoint per subinterval decides the whole sweep.
    ///
    /// The endpoints `0` and `1` are included as candidates: the faces
    /// themselves are cutting planes, and a positive saddle on one of them is
    /// the face rule's own answer.
    #[must_use]
    pub fn test(&self) -> Interior {
        let mut breaks = [R::ZERO; 3];
        let mut count = 0;
        for t in self.numerator_roots().chain(self.pole()) {
            if t > R::ZERO && t < R::ONE {
                breaks[count] = t;
                count += 1;
            }
        }
        let breaks = &mut breaks[..count];
        // Insertion sort: at most three entries, and it avoids depending on a
        // sort that `no_std` slices only offer for `Ord`, which `R` is not.
        for i in 1..breaks.len() {
            let mut j = i;
            while j > 0 && breaks[j - 1] > breaks[j] {
                breaks.swap(j - 1, j);
                j -= 1;
            }
        }

        // The faces, then one interior point of every subinterval between them.
        let half = R::ONE / (R::ONE + R::ONE);
        if self.saddle(R::ZERO) > R::ZERO || self.saddle(R::ONE) > R::ZERO {
            return Interior::Joined;
        }
        let mut previous = R::ZERO;
        for k in 0..=breaks.len() {
            let next = breaks.get(k).copied().unwrap_or(R::ONE);
            let mid = (previous + next) * half;
            // A midpoint coincides with a breakpoint only if the subinterval is
            // empty, which happens when two breakpoints are equal; the saddle is
            // then zero or undefined there and the neighbouring subintervals
            // carry the answer.
            if mid > previous && mid < next && self.saddle(mid) > R::ZERO {
                return Interior::Joined;
            }
            previous = next;
        }
        Interior::Separated
    }

    /// The real roots of `F` in `(0, 1)`, as the breakpoints the sign walk needs.
    ///
    /// Returns at most two. The quadratic is solved from the derived
    /// coefficients; when `a` is zero it is linear, which is not a special case
    /// so much as a smaller polynomial, and is solved as one.
    fn numerator_roots(&self) -> impl Iterator<Item = R> {
        let d = |k: usize| self.hi[k] - self.lo[k];
        let a = d(0) * d(2) - d(1) * d(3);
        let b = self.lo[2] * d(0) + self.lo[0] * d(2) - self.lo[3] * d(1) - self.lo[1] * d(3);
        let c = self.lo[0] * self.lo[2] - self.lo[1] * self.lo[3];

        let mut roots = [R::ZERO; 2];
        let mut count = 0;
        if a == R::ZERO {
            if b != R::ZERO {
                roots[0] = -c / b;
                count = 1;
            }
        } else {
            let two = R::ONE + R::ONE;
            let disc = b * b - two * two * a * c;
            if disc >= R::ZERO {
                let root = disc.sqrt();
                roots[0] = (-b - root) / (two * a);
                roots[1] = (-b + root) / (two * a);
                count = 2;
            }
        }
        roots.into_iter().take(count)
    }
}

/// Chernyaev's numerator-only test, kept purely as a cross-check.
///
/// **This is not the test the crate uses**, in the same way and for the same
/// reason [`super::reference`] is not the table it uses: it exists so that
/// [`SweptFaces::test`] can be checked against the construction it corrects, and
/// so that the disagreement Custodio §5.1 predicts can be *demonstrated* rather
/// than described.
///
/// It reports whether `F` is positive anywhere in `[0, 1]`, which is the sign
/// Chernyaev's three conditions track. Where `Δ > 0` throughout the sweep the
/// two agree, because dividing by a positive number does not change a sign.
#[cfg(test)]
pub(super) fn chernyaev_numerator_test<R: Real>(faces: &SweptFaces<R>) -> Interior {
    // The maximum of a quadratic over a closed interval is at an endpoint or at
    // its vertex, so three evaluations decide it — no sweep required.
    let d = |k: usize| faces.hi[k] - faces.lo[k];
    let a = d(0) * d(2) - d(1) * d(3);
    let b = faces.lo[2] * d(0) + faces.lo[0] * d(2) - faces.lo[3] * d(1) - faces.lo[1] * d(3);
    let two = R::ONE + R::ONE;

    let mut positive = faces.numerator(R::ZERO) > R::ZERO || faces.numerator(R::ONE) > R::ZERO;
    if a < R::ZERO {
        let vertex = -b / (two * a);
        if vertex > R::ZERO && vertex < R::ONE {
            positive = positive || faces.numerator(vertex) > R::ZERO;
        }
    }
    if positive {
        Interior::Joined
    } else {
        Interior::Separated
    }
}
