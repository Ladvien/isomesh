//! The dual-contouring vertex rule: one unconditional path.
//!
//! Given the crossings on a cell's edges — a position and a surface normal per
//! crossing, which is what A-006's [`HermiteCell`]
//! holds — place the cell's vertex where those tangent planes best agree.
//!
//! ```text
//! c = centroid of the crossing positions        dᵢ = nᵢ·(pᵢ − c)
//! M = Σ nᵢnᵢᵀ        g = Σ dᵢnᵢ        λ = 0.01
//! x = c + adj(M + λI)·g / det(M + λI)
//! ```
//!
//! Source: `docs/research/2026-08-10-adjacent-math-transfer-audit.md` §3, which
//! derives it from classical invariant theory — see Villar et al., *Scalars are
//! universal*, `10.48550/arxiv.2106.06610`, and Fels & Olver on moving frames,
//! `10.1023/a:1005878210297`. The QEF it approximates is Ju, Losasso, Schaefer &
//! Warren, *Dual Contouring of Hermite Data*, SIGGRAPH 2002,
//! `10.1145/566570.566586` §2.3.
//!
//! # Why there is no fast path for three planes
//!
//! The audit gives an exact Cramer rule for the 3-crossing case, and the backlog
//! originally split it into two tickets with a triple-product test between them.
//! That split was a misreading and is recorded as ✗12. Two reasons it is gone:
//!
//! - The audit's diagnosis of *why* dual contouring pops is the branch itself.
//!   Over 20,000 trials seeded at DC's `σ = 0.1` SVD threshold in `f32`, the rank
//!   branch disagreed after a rotation in **454 cases**, and when it flipped the
//!   vertex moved a median of **2.13 cells**, at worst **9.10**. A triple-product
//!   threshold is the same construction with a different discriminant.
//! - The dropped path is not the accurate one. Measured equivariance residual,
//!   `f32`, 4000 random cells: Cramer p99 `7.23e−04` and max `3.6e−01`, against
//!   this form's `1.81e−04` and `6.4e−04`. It is better in the tail by three
//!   orders of magnitude, and the tail is what a user sees.
//!
//! So: no eigendecomposition, no SVD, no iteration, and no data-dependent
//! branch. Around 90 flops.
//!
//! # Why `λ = 0.01`, and why it is unitless
//!
//! `λ` is the Tikhonov regularizer that stops an under-determined cell — a flat
//! region, where `M` is rank 1 — from flying off to where the two constraints it
//! does have happen to meet. The audit's value **reproduces DC's `σ = 0.1`
//! truncation smoothly**, which is the whole point: the same regularisation
//! without the discontinuity.
//!
//! It needs no length scale. The normals are unit, so `M` is dimensionless and
//! `λ` adds to a dimensionless quantity; `g` carries the length units and comes
//! back out through `adj·g/det`. A cell measured in metres and the same cell in
//! millimetres get the same answer, which is not true of a `λ` with units.
//!
//! Note the corpus circulates three constants — this `0.01`, Subgrid MT's QEF
//! `λ = 0.1`, and DC's SVD threshold `σ = 0.1`. They are not interchangeable.
//!
//! # `M = AᵀA` squares the condition number
//!
//! Forming `M = Σnᵢnᵢᵀ` is forming `AᵀA`, and that squares the condition number
//! of the underlying system. This is the exact reason the QR/Givens formulation
//! exists in the literature: DC's own paper measures `bᵀb` reaching `~10⁶` on a
//! 256³ grid, so in `f32` — six decimal digits — `E[x]` evaluated on a flat
//! region has error *on the order of 1*, and recommends `f64` (V-18).
//!
//! Two things make that less alarming here than it sounds. The regularizer keeps
//! `A = M + λI` away from singular by construction, and the adjugate form never
//! divides by a small pivot — it divides once, by `det(A) ≥ λ³ > 0`. It is still
//! the reason `Real` spans `f64` at all, and E-112 is the ticket that measures
//! where `f32` gives out.

use crate::Real;
use crate::hermite::HermiteCell;
use crate::vec3;

/// The Tikhonov regularizer. Unitless — see the module docs.
pub const LAMBDA: f64 = 0.01;

/// A dot product summed smallest-magnitude-first.
///
/// Floating-point addition is not associative, so `a·b` depends on the order the
/// three products are summed, and the order a naive loop uses depends on which
/// axis is which. Under a 90° rotation the components permute, the summation
/// order changes with them, and the result differs in the last bits — so the
/// same cell rotated onto itself does not give the same vertex.
///
/// The audit measured that directly: **4328 of 9600** lattice-symmetry trials
/// disagreed with an unsorted dot product, and **0 of 9600** with a sorted one.
/// Sorting by magnitude is a permutation-invariant order, so the sum is a
/// function of the *set* of products rather than of the axis labelling.
///
/// Three elements, so this is a sorting network rather than a sort: five
/// comparisons, no allocation, no branch on data beyond the swaps.
#[inline]
pub fn dot_equivariant<R: Real>(a: [R; 3], b: [R; 3]) -> R {
    sum_equivariant([a[0] * b[0], a[1] * b[1], a[2] * b[2]])
}

/// Sort ascending by magnitude, in place.
///
/// `total_cmp` on the absolute values is a total order even across NaN and
/// signed zero, so the network cannot produce an order that depends on which
/// comparison ran first. Ties need no tie-break: IEEE addition and
/// multiplication are both *commutative* — only associativity fails — so two
/// terms of equal magnitude give the same answer either way round.
#[inline]
fn sort_by_magnitude<R: Real, const N: usize>(t: &mut [R; N]) {
    // Insertion sort: N is 3 or 5 here, and it is branch-predictable and
    // allocation-free, which a comparison sort over a slice would not be in
    // `no_std`.
    let mut i = 1;
    while i < N {
        let mut j = i;
        while j > 0 && t[j - 1].abs().total_cmp(&t[j].abs()) == core::cmp::Ordering::Greater {
            t.swap(j - 1, j);
            j -= 1;
        }
        i += 1;
    }
}

/// Sum smallest-magnitude-first, so the result is a function of the *set* of
/// terms rather than of the order they were written in.
#[inline]
fn sum_equivariant<R: Real, const N: usize>(mut t: [R; N]) -> R {
    sort_by_magnitude(&mut t);
    let mut acc = R::ZERO;
    let mut i = 0;
    while i < N {
        acc += t[i];
        i += 1;
    }
    acc
}

/// Multiply smallest-magnitude-first, for the same reason.
///
/// Needed because floating-point multiplication is commutative but **not
/// associative**: `(a·b)·c` and `(b·c)·a` differ in the last bits, and a lattice
/// rotation permutes exactly which of the matrix entries play the part of `a`,
/// `b` and `c` in the determinant's terms.
#[inline]
fn mul_equivariant<R: Real>(mut t: [R; 3]) -> R {
    sort_by_magnitude(&mut t);
    (t[0] * t[1]) * t[2]
}

/// A symmetric 3×3 matrix, stored as its six distinct entries.
///
/// ```text
/// [ xx xy xz ]
/// [ xy yy yz ]
/// [ xz yz zz ]
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
struct Symmetric3<R: Real> {
    xx: R,
    xy: R,
    xz: R,
    yy: R,
    yz: R,
    zz: R,
}

impl<R: Real> Symmetric3<R> {
    const ZERO: Self = Self {
        xx: R::ZERO,
        xy: R::ZERO,
        xz: R::ZERO,
        yy: R::ZERO,
        yz: R::ZERO,
        zz: R::ZERO,
    };

    /// Accumulate the outer product `n nᵀ`.
    #[inline]
    fn add_outer(&mut self, n: [R; 3]) {
        self.xx += n[0] * n[0];
        self.xy += n[0] * n[1];
        self.xz += n[0] * n[2];
        self.yy += n[1] * n[1];
        self.yz += n[1] * n[2];
        self.zz += n[2] * n[2];
    }

    /// Add `λ` to the diagonal.
    #[inline]
    fn regularized(mut self, lambda: R) -> Self {
        self.xx += lambda;
        self.yy += lambda;
        self.zz += lambda;
        self
    }

    /// The adjugate — the transpose of the cofactor matrix — which is symmetric
    /// when the input is, and `det(A) · A⁻¹` when `A` is invertible.
    ///
    /// Used instead of an inverse so there is exactly one division in the whole
    /// solve, and it is by a quantity bounded below by `λ³`.
    #[inline]
    fn adjugate(self) -> Self {
        Self {
            xx: self.yy * self.zz - self.yz * self.yz,
            xy: self.xz * self.yz - self.xy * self.zz,
            xz: self.xy * self.yz - self.xz * self.yy,
            yy: self.xx * self.zz - self.xz * self.xz,
            yz: self.xy * self.xz - self.xx * self.yz,
            zz: self.xx * self.yy - self.xy * self.xy,
        }
    }

    /// `det(self)`, in the symmetric form.
    ///
    /// **Not** a cofactor expansion along a row. A row expansion picks three of
    /// the six entries by *position*, and a lattice rotation relabels the axes —
    /// so the rotated matrix gets expanded along a different set of entries, the
    /// arithmetic associates differently, and the result differs in the last
    /// bits. That is a real bug and it cost the equivariance test a failure
    /// before this form replaced it.
    ///
    /// ```text
    /// det = xx·yy·zz + 2·xy·yz·xz − xx·yz² − yy·xz² − zz·xy²
    /// ```
    ///
    /// Under a permutation of the axes the first two terms map to themselves and
    /// the last three permute among each other, so the *multiset* of five terms
    /// is invariant. Summing them magnitude-first, with each three-factor product
    /// also ordered by magnitude, makes the evaluation invariant too.
    #[inline]
    fn determinant(self) -> R {
        let two = R::TWO;
        sum_equivariant([
            mul_equivariant([self.xx, self.yy, self.zz]),
            two * mul_equivariant([self.xy, self.yz, self.xz]),
            -mul_equivariant([self.xx, self.yz, self.yz]),
            -mul_equivariant([self.yy, self.xz, self.xz]),
            -mul_equivariant([self.zz, self.xy, self.xy]),
        ])
    }

    /// `self · v`.
    #[inline]
    fn mul_vec(self, v: [R; 3]) -> [R; 3] {
        [
            dot_equivariant([self.xx, self.xy, self.xz], v),
            dot_equivariant([self.xy, self.yy, self.yz], v),
            dot_equivariant([self.xz, self.yz, self.zz], v),
        ]
    }
}

/// Where a cell's vertex goes.
///
/// Returns `None` only when the cell has no crossings at all, which the caller
/// has usually already ruled out by checking corner signs. Every other case —
/// one crossing, twelve, all coplanar, all parallel — is handled by the same
/// arithmetic, which is the property the whole module is built around.
///
/// The result is **not** clamped to the cell. A-009 is the ticket that adds the
/// clamp and measures what it costs in sharpness and buys in
/// self-intersections; until then this is the unmodified solve, so that A-009
/// has an unclamped baseline to measure against.
#[must_use]
pub fn solve<R: Real>(cell: &HermiteCell<R>) -> Option<[R; 3]> {
    let centroid = cell.centroid()?;

    let mut m = Symmetric3::<R>::ZERO;
    let mut g = [R::ZERO; 3];
    for crossing in cell.iter() {
        let n = crossing.normal;
        // Plane distance from the centroid, along this crossing's own normal.
        let d = dot_equivariant(n, vec3::sub(crossing.position, centroid));
        m.add_outer(n);
        g = [g[0] + n[0] * d, g[1] + n[1] * d, g[2] + n[2] * d];
    }

    let a = m.regularized(R::from_f64(LAMBDA));
    let adj = a.adjugate();
    let det = a.determinant();

    // `det >= λ³` for any set of unit normals, since `M` is positive
    // semi-definite and `λ > 0`, so this division is safe by construction rather
    // than by a guard. The finiteness check below is for a caller who supplied a
    // non-finite normal, not for the algebra.
    let offset = vec3::scale(adj.mul_vec(g), det.recip());
    let x = [
        centroid[0] + offset[0],
        centroid[1] + offset[1],
        centroid[2] + offset[2],
    ];
    if x[0].is_finite() && x[1].is_finite() && x[2].is_finite() {
        Some(x)
    } else {
        None
    }
}

#[cfg(test)]
mod tests;
