//! Per-cell isotopy certificate, from normal variation.
//!
//! Ticket: T-015. Plantinga & Vegter, *Isotopic approximation of implicit curves
//! and surfaces*, SGP 2004 (`10.1145/1057432.1057465`).
//!
//! # Hausdorff distance does not certify topology, and cannot
//!
//! Two surfaces can be arbitrarily Hausdorff-close and not homeomorphic — a
//! sphere and a sphere with a hairline handle differ by as little as you like in
//! distance and are different objects. So every geometric error this crate
//! reports, however small, says nothing about whether the mesh has the right
//! topology. Every real theorem in this space adds a second hypothesis, and the
//! isosurface-specific ones make it a bound on how much the **normal** turns.
//!
//! # The condition, transcribed
//!
//! Plantinga & Vegter, verbatim: a grid `G` such that for each cell `C`
//!
//! ```text
//! 0 ∉ □F(C)  ∨  ⟨□∇F(C), □∇F(C)⟩ > 0
//! ```
//!
//! where `□` is an interval bound over the cell. Their reasoning for why the
//! second clause is what matters: *"Due to the inner product constraint, S is
//! parametrizable in the direction of one of the axes. Therefore we cannot have
//! alternating signs of F at the vertices of C, since F would have to increase
//! along one edge, and decrease along the other parallel edge."* A cell that
//! passes has at most one surface component in it, isotopic to a disc, and the
//! single facet the mesher emits is the right answer *topologically* and not
//! merely nearby.
//!
//! # Why this is exact here rather than an interval-arithmetic approximation
//!
//! The general form needs interval arithmetic over an arbitrary `F`, which this
//! crate has no way to do — an [`Sdf`](crate::Sdf) hands back point values. A
//! sampled hull of the gradient would be a *lower* bound on its variation, so
//! the predicate could pass where the truth fails, which is the one direction a
//! certificate must never err in.
//!
//! But the surface Marching Cubes actually approximates is not `F`. It is the
//! **trilinear interpolant** of the eight corner values — that is the whole
//! subject of the `trilinear` module and of Grosso's papers — and for a
//! trilinear function the interval bounds are exact and closed-form:
//!
//! - `∂F/∂x` is *bilinear* in `(y, z)`, so over the cell it is a convex
//!   combination of the four `x`-edge differences. Its exact range is their min
//!   and max. Likewise for `y` and `z`.
//! - `F` itself is a convex combination of the eight corner values, so
//!   `0 ∉ □F(C)` is exactly *"all eight corners share a sign"* — an inactive
//!   cell. **The first clause is therefore free**, and for active cells the
//!   predicate is the inner product alone.
//!
//! So this is a genuine certificate against the surface the extractor is
//! approximating, with no interval library and no sampling. What it does not
//! certify is the analytic field against its trilinear interpolant — that is a
//! separate question, and it is the one the sampling error already measures.
//!
//! # The cell size cancels
//!
//! `∂F/∂x` is a corner difference divided by `h`, and so are the other two. The
//! predicate tests the **sign** of a sum of three such squares, so `h²` factors
//! out and the arithmetic runs on raw corner differences. That is only true for
//! **isotropic** cells; an anisotropic grid would need each axis divided by its
//! own spacing before the sum.

#[cfg(test)]
mod tests;

use crate::real::Real;
use crate::shape::Shape3;

/// The four corner index pairs whose differences give `∂F/∂x` over a cell.
///
/// Corner `i` sits at `(i&1, (i>>1)&1, (i>>2)&1)` — `cube.rs`'s numbering, bit 0
/// being `x`. So an `x` difference pairs corners differing only in bit 0.
const X_PAIRS: [[usize; 2]; 4] = [[0, 1], [2, 3], [4, 5], [6, 7]];
/// Likewise for `∂F/∂y`: corners differing only in bit 1.
const Y_PAIRS: [[usize; 2]; 4] = [[0, 2], [1, 3], [4, 6], [5, 7]];
/// Likewise for `∂F/∂z`: corners differing only in bit 2.
const Z_PAIRS: [[usize; 2]; 4] = [[0, 4], [1, 5], [2, 6], [3, 7]];

/// Exact range of one partial derivative over the cell, in corner-difference
/// units.
fn partial_range<R: Real>(corner: &[R; 8], pairs: &[[usize; 2]; 4]) -> [R; 2] {
    let mut lo = corner[pairs[0][1]] - corner[pairs[0][0]];
    let mut hi = lo;
    for p in &pairs[1..] {
        let d = corner[p[1]] - corner[p[0]];
        if d < lo {
            lo = d;
        }
        if d > hi {
            hi = d;
        }
    }
    [lo, hi]
}

/// Lower bound of the interval `[a, b] · [a, b]`.
///
/// Interval multiplication takes the minimum of the four corner products. For a
/// self-product those are `a²`, `ab` (twice) and `b²`, so the minimum is `a²` or
/// `b²` when the interval keeps a sign and `ab ≤ 0` when it straddles zero.
/// **That straddle is the whole test**: an axis whose derivative changes sign
/// inside the cell contributes a negative term and can only be rescued by
/// another axis dominating it.
fn self_product_low<R: Real>([a, b]: [R; 2]) -> R {
    if a > R::ZERO {
        a * a
    } else if b < R::ZERO {
        b * b
    } else {
        a * b
    }
}

/// Does this cell satisfy Plantinga & Vegter's condition?
///
/// `corner` is the eight corner values in `cube.rs`'s numbering. `true` means
/// the trilinear surface inside the cell is isotopic to what a single-facet
/// approximation produces — one component, no handle, no second sheet.
///
/// Inactive cells pass trivially by the first clause, since `0 ∉ □F(C)` is
/// exactly "all eight corners share a sign" for a trilinear interpolant.
#[must_use]
pub fn cell_is_certified<R: Real>(corner: &[R; 8]) -> bool {
    // Clause one: `0 ∉ □F(C)`. Exact, because `F` is a convex combination of
    // the corners.
    let inside = corner[0] < R::ZERO;
    if corner.iter().all(|v| (*v < R::ZERO) == inside) {
        return true;
    }

    // Clause two: `⟨□∇F, □∇F⟩ > 0`, the interval inner product of the gradient
    // bound with itself.
    let x = self_product_low(partial_range(corner, &X_PAIRS));
    let y = self_product_low(partial_range(corner, &Y_PAIRS));
    let z = self_product_low(partial_range(corner, &Z_PAIRS));
    x + y + z > R::ZERO
}

/// What the certificate found over a whole grid.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct IsotopyReport {
    /// Cells examined.
    pub cells: u64,
    /// Cells the surface passes through — those with a sign change.
    pub active_cells: u64,
    /// Active cells satisfying the condition.
    pub certified: u64,
    /// Active cells that do not. **Not a count of errors**: an uncertified cell
    /// may still be meshed correctly; what is lost is the guarantee.
    pub uncertified: u64,
    /// The first uncertified cell, by grid index, for a caller that wants to
    /// look at it.
    pub first_failure: Option<[u32; 3]>,
}

impl IsotopyReport {
    /// Share of active cells that carry the guarantee, in `[0, 1]`.
    ///
    /// One means the whole extracted surface is isotopic to the trilinear
    /// isosurface. Anything less means some region is not certified — which is
    /// weaker than saying it is wrong.
    #[must_use]
    pub fn certified_fraction(&self) -> f64 {
        if self.active_cells == 0 {
            1.0
        } else {
            self.certified as f64 / self.active_cells as f64
        }
    }

    /// Every active cell is certified.
    #[must_use]
    pub const fn is_certified(&self) -> bool {
        self.uncertified == 0
    }
}

/// Evaluate the certificate over a sampled grid.
///
/// # Errors
///
/// [`Error::GridTooSmall`](crate::Error::GridTooSmall) for a grid under 2×2×2,
/// [`Error::ShapeOverflow`](crate::Error::ShapeOverflow) if `values` is not one
/// entry per sample.
pub fn isotopy_report<R: Real>(values: &[R], shape: &impl Shape3) -> crate::Result<IsotopyReport> {
    let size = shape.size();
    if size[0] < 2 || size[1] < 2 || size[2] < 2 {
        return Err(crate::Error::GridTooSmall { size });
    }
    if values.len() != shape.element_count() {
        return Err(crate::Error::ShapeOverflow {
            size,
            product: values.len() as u64,
        });
    }

    let (nx, ny) = (size[0], size[1]);
    let mut report = IsotopyReport::default();

    for z in 0..size[2] - 1 {
        for y in 0..ny - 1 {
            for x in 0..nx - 1 {
                let at = |dx: u32, dy: u32, dz: u32| {
                    values[(((z + dz) * ny + (y + dy)) * nx + (x + dx)) as usize]
                };
                let mut corner = [R::ZERO; 8];
                for (i, slot) in corner.iter_mut().enumerate() {
                    let i = i as u32;
                    *slot = at(i & 1, (i >> 1) & 1, (i >> 2) & 1);
                }

                report.cells += 1;
                let inside = corner[0] < R::ZERO;
                if corner.iter().all(|v| (*v < R::ZERO) == inside) {
                    continue;
                }
                report.active_cells += 1;
                if cell_is_certified(&corner) {
                    report.certified += 1;
                } else {
                    report.uncertified += 1;
                    if report.first_failure.is_none() {
                        report.first_failure = Some([x, y, z]);
                    }
                }
            }
        }
    }

    Ok(report)
}
