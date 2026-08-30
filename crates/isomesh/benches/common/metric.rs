//! Metric-based anisotropic mesh adaptation — the machinery of Loseille &
//! Alauzet, driven entirely through the `isomesh` public API.
//!
//! Ticket: R-146, the `L` row and the longest pole of Group D. Consumed
//! unchanged by R-147 (the rate-versus-constant decomposition), R-148 (metric
//! interpolation across a chunk seam), R-149 (the per-feature mechanism),
//! R-150 (the within-a-norm methodology guard) and R-151 (goal-oriented
//! refinement). Five benches read one eigensolver on purpose: a second copy of
//! a Jacobi sweep is two answers to one question, and the question here is
//! whether a metric's smallest eigenvalue is a measurement or an artifact.
//!
//! # Where `M_Lp` comes from
//!
//! Interpolating a function `u` by a piecewise-linear reconstruction over an
//! element `K` leaves an error controlled entirely by the Hessian: the classical
//! bound is `|u − Π_h u|_{∞,K} ≲ max_{x,y ∈ K} ⟨y − x, |H_u| (y − x)⟩`, where
//! `|H_u|` is the Hessian with the absolute value taken on its spectrum — the
//! sign of a curvature does not make an element cheaper. A Riemannian metric
//! `M` calls an element *unit* when that same quadratic form equals one, so
//! declaring `M ∝ |H_u|` equidistributes the error and makes the element's
//! anisotropy follow the curvature's anisotropy, direction for direction.
//!
//! That fixes only the *ratios* of `M`'s eigenvalues. The local *scale* is still
//! free, and it is pinned by minimising the global `L^p` error subject to a
//! fixed budget. The budget is the **complexity**
//! `C(M) = ∫_Ω √det M`, which stands in for the vertex count because a
//! unit-metric simplex occupies Euclidean volume `1/√det M`, so `√det M` is a
//! point density. The stationary point of that constrained minimisation is
//!
//! ```text
//!     M_Lp = D_Lp · det(|H_u|)^(−1/(2p + d)) · |H_u|,        d = 3,
//! ```
//!
//! with `D_Lp` the single global Lagrange multiplier that pins `C(M)` to the
//! requested budget. Source: **NASA NTRS 20200003084**, which restates Loseille
//! & Alauzet verbatim; the two SIAM originals are paywalled and the corpus holds
//! only their landing pages, so the restatement is the primary this module was
//! written from and is the one a finding must cite.
//!
//! The Hessian itself needs nothing new from the crate. `isomesh` already
//! samples the trilinear and already differences at cell size (`M-65`), and
//! [`hessian`] is the same stencil at the same step.
//!
//! # The two constants this module chooses, and why
//!
//! **`D_LP = 1`.** `D_Lp` multiplies the entire field by one scalar. Every
//! quantity Group D reports either divides it out or is blind to it: a
//! complexity *ratio* between two arms cancels `D_Lp^{3/2}`, [`aspect_ratio`] is
//! scale-free by construction, [`am_gm_gap`] is a ratio of two norms over the
//! same population, and [`density_from_metric`] is normalised to unit mean. So
//! `D_Lp` is folded to 1 and a budget is matched by moving the sample count,
//! never by fitting the multiplier. Recording `D_Lp` as a number would be
//! recording a choice of units.
//!
//! **`H_FLOOR = 1e-9`.** A zero Hessian eigenvalue is not an edge case to be
//! papered over — it is *exactly* the flat direction the whole of Group D is
//! about, and it is the common case: a cylinder, a ridge, a CSG box face and a
//! plane all have one. But `det(|H|) = 0` makes `det^(−1/(2p+3))` infinite, and
//! an infinity that silently propagates into a metric is a second path through
//! the same code. So [`metric_lp`] floors each `|eigenvalue|` at `H_FLOOR`
//! before forming the determinant. The value is absolute and in the units of a
//! second derivative (reciprocal length): the reference fields' genuine
//! curvatures run from about `1e-2` to about `1e2`, so `1e-9` sits seven orders
//! below the smallest real one and cannot mask a curvature, while capping
//! `det^(−1/(2p+3))` at a factor of about twenty for `p = 2`.
//!
//! **The floor is visible in the output and a consumer must say so.** Where a
//! direction is genuinely flat, [`aspect_ratio`] returns
//! `|λ|max / H_FLOOR` — a number of order `1e11` that is *the floor talking*,
//! not a measured anisotropy. `H_FLOOR` is public for precisely this reason: a
//! bench reporting `aspect_ratio_max` should also report how many cells sit at
//! the floor, or its maximum is a restatement of this constant. R-146's vacuity
//! control ("maximum aspect ratio above 3 on at least one field") is only
//! informative when it is met by cells that are *not* at the floor.
//!
//! # The eigensolver
//!
//! Cyclic symmetric Jacobi, sweeping `(0,1) → (0,2) → (1,2)`, at most
//! [`JACOBI_SWEEPS`] `= 12` sweeps, exiting early when the off-diagonal
//! Frobenius norm falls to [`JACOBI_TOLERANCE`] `= 1e-14` times the **full**
//! Frobenius norm of the matrix.
//!
//! The scale is the Frobenius norm and deliberately **not** the trace, even
//! though the trace is the more usual choice: a Hessian is indefinite, and
//! `diag(1, −1, 0)` is a perfectly ordinary saddle whose trace is zero. Against
//! a trace scale its tolerance would be identically zero and the early exit
//! could never fire. `‖M‖_F = √(Σ λ²)` is invariant under Jacobi rotations,
//! never zero for a non-zero matrix, and is the one scale the off-diagonal norm
//! is directly comparable to. The `12`-sweep cap makes convergence
//! unconditional either way — for `3 × 3`, three or four sweeps is normal — so
//! the criterion only decides how early the loop stops, never whether it
//! terminates.
//!
//! Eigenvalues come back **ascending**, tie-broken lexicographically on the
//! eigenvector column and then on the original axis index, so a degenerate
//! spectrum still produces one ordering rather than whichever one the sweep
//! happened to leave. Eigenvectors are orthonormal to round-off and each column
//! is signed so that its first component of magnitude above `1e-12` is
//! positive. No `HashMap`, no `partial_cmp`, no `unwrap`.
//!
//! # Why this module is `f64` throughout and not generic over `Real`
//!
//! The metric algebra needs `exp`, `ln` and `powf`: the log-Euclidean scheme is
//! a matrix logarithm and `M_Lp` is a fractional power of a determinant.
//! `isomesh::Real` has `sqrt`, `sin`, `cos`, `atan2` and `acos` and none of
//! those three (real.rs), and it is sealed, so there is nowhere to add them.
//! Every signature here is concrete `f64`. A caller working in `f32` converts at
//! the boundary — which is the honest place for that loss, since a second
//! difference divided by `h²` is the one quantity in the crate that cannot
//! afford `f32`.

use isomesh::Sdf;
use std::cmp::Ordering;

/// The global Lagrange constant of the optimal `L^p` metric, folded to 1.
///
/// See the module header: it cancels in every quantity this phase reports.
pub(crate) const D_LP: f64 = 1.0;

/// The absolute floor applied to each `|Hessian eigenvalue|` before the
/// `det(|H|)^(−1/(2p+d))` factor is formed.
///
/// See the module header for the derivation of the value and for the warning
/// that an `aspect_ratio` near `|λ|max / H_FLOOR` is reporting this constant
/// rather than a measurement.
pub(crate) const H_FLOOR: f64 = 1e-9;

/// Maximum cyclic Jacobi sweeps. Three or four is normal for `3 × 3`; twelve is
/// a cap that makes termination unconditional, not an expectation.
pub(crate) const JACOBI_SWEEPS: usize = 12;

/// Early-exit tolerance for the Jacobi sweep, relative to `‖M‖_F`.
pub(crate) const JACOBI_TOLERANCE: f64 = 1e-14;

/// Below this gradient magnitude the level set has no tangent plane and
/// [`principal_curvatures`] returns `None`.
///
/// Absolute, because the reference fields are eikonal or near-eikonal: away
/// from a critical point `‖∇f‖ ≈ 1`, so `1e-12` names critical points and
/// nothing else.
pub(crate) const GRAD_FLOOR: f64 = 1e-12;

/// A component of an eigenvector counts as non-zero for the sign convention
/// when its magnitude exceeds this. Columns are unit length, so the test is
/// already relative.
const SIGN_EPS: f64 = 1e-12;

/// `(i, j)` → index into [`Sym3`]'s six-entry array.
const IDX: [[usize; 3]; 3] = [[0, 1, 2], [1, 3, 4], [2, 4, 5]];

/// A symmetric 3x3 matrix, stored as the six upper entries
/// `[xx, xy, xz, yy, yz, zz]`.
///
/// Six `f64`s rather than nine: a Hessian and a metric are symmetric by
/// construction, and storing the lower triangle separately invites the two
/// halves to disagree after a rotation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Sym3(pub(crate) [f64; 6]);

impl Sym3 {
    /// The identity metric — unit-isotropic, the "uniform grid" arm's metric.
    pub(crate) fn identity() -> Self {
        Self([1.0, 0.0, 0.0, 1.0, 0.0, 1.0])
    }

    /// The zero matrix. A valid Hessian (a linear field has one) and never a
    /// valid metric.
    pub(crate) fn zero() -> Self {
        Self([0.0; 6])
    }

    /// Entry `(i, j)`, either triangle. Panics for an index outside `0..3`.
    pub(crate) fn get(&self, i: usize, j: usize) -> f64 {
        self.0[IDX[i][j]]
    }

    /// Determinant, by cofactor expansion on the first row.
    ///
    /// For a metric out of [`metric_lp`] this is bounded below by
    /// `(D_LP · H_FLOOR)³` times the two larger floored eigenvalues, so the
    /// cancellation in the expansion — of order `ε · |λ|max³` — stays five
    /// orders below the answer even at the floor.
    pub(crate) fn det(&self) -> f64 {
        let [xx, xy, xz, yy, yz, zz] = self.0;
        xx * (yy * zz - yz * yz) - xy * (xy * zz - yz * xz) + xz * (xy * yz - yy * xz)
    }

    /// Trace. For `|H|` (see [`Sym3::abs`]) this is `Σ|λ|`, which is `d` times
    /// the arithmetic mean of the curvature magnitudes — the denominator of
    /// [`am_gm_gap`].
    pub(crate) fn trace(&self) -> f64 {
        self.0[0] + self.0[3] + self.0[5]
    }

    /// Jacobi eigendecomposition. Returns (eigenvalues ascending, eigenvectors
    /// as columns).
    ///
    /// `vectors[row][col]` is component `row` of eigenvector `col`, so
    /// `values[c]` pairs with the column `c`. Convergence criterion, sweep cap,
    /// ordering and sign convention are all stated in the module header.
    ///
    /// Panics on a non-finite entry: a NaN pivot would flow silently into every
    /// eigenvalue and every consumer, and a field that sampled NaN is a defect
    /// upstream of this module rather than a case to be handled here.
    pub(crate) fn eigen(&self) -> ([f64; 3], [[f64; 3]; 3]) {
        assert!(
            self.0.iter().all(|entry| entry.is_finite()),
            "Sym3::eigen on a non-finite matrix: {self:?}"
        );

        let [xx, xy, xz, yy, yz, zz] = self.0;
        let mut work = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]];
        let mut basis = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

        // ‖M‖_F is invariant under the sweep, so it is fixed once here.
        let frobenius = (xx * xx + yy * yy + zz * zz + 2.0 * (xy * xy + xz * xz + yz * yz)).sqrt();
        let tolerance = JACOBI_TOLERANCE * frobenius;

        for _ in 0..JACOBI_SWEEPS {
            let off = (2.0
                * (work[0][1] * work[0][1] + work[0][2] * work[0][2] + work[1][2] * work[1][2]))
                .sqrt();
            if off <= tolerance {
                break;
            }
            for (row, col) in [(0usize, 1usize), (0, 2), (1, 2)] {
                let pivot = work[row][col];
                if pivot.abs() <= 0.0 {
                    continue;
                }

                // The rotation that annihilates `work[row][col]`:
                // θ = (a_qq − a_pp) / 2a_pq, t = sgn θ / (|θ| + √(θ² + 1)).
                // Taking the smaller root keeps |t| ≤ 1, which is what makes
                // the sweep unconditionally stable. A huge θ (a pivot
                // vanishingly small beside the diagonal gap) sends `t` to zero,
                // i.e. no rotation — the correct limit, reached without a
                // branch.
                let diff = work[col][col] - work[row][row];
                let theta = diff / (2.0 * pivot);
                let sign = if theta < 0.0 { -1.0 } else { 1.0 };
                let tan = sign / (theta.abs() + (theta * theta + 1.0).sqrt());
                let cos = 1.0 / (tan * tan + 1.0).sqrt();
                let sin = tan * cos;
                // s/(1+c) form of the update: algebraically `a'_rp = c·a_rp −
                // s·a_rq` but written so the correction is added to the old
                // value rather than the value being rebuilt from scratch.
                let half_tan = sin / (1.0 + cos);

                work[row][row] -= tan * pivot;
                work[col][col] += tan * pivot;
                work[row][col] = 0.0;
                work[col][row] = 0.0;

                // The third index: (0,1)→2, (0,2)→1, (1,2)→0.
                let other = 3 - row - col;
                let from_row = work[other][row];
                let from_col = work[other][col];
                let to_row = from_row - sin * (from_col + half_tan * from_row);
                let to_col = from_col + sin * (from_row - half_tan * from_col);
                work[other][row] = to_row;
                work[row][other] = to_row;
                work[other][col] = to_col;
                work[col][other] = to_col;

                for basis_row in &mut basis {
                    let old_row = basis_row[row];
                    let old_col = basis_row[col];
                    basis_row[row] = cos * old_row - sin * old_col;
                    basis_row[col] = sin * old_row + cos * old_col;
                }
            }
        }

        let raw = [work[0][0], work[1][1], work[2][2]];

        // Ascending, tie-broken lexicographically on the column and finally on
        // the axis index, so a degenerate spectrum has exactly one ordering.
        let mut order = [0usize, 1, 2];
        order.sort_by(|&left, &right| {
            raw[left].total_cmp(&raw[right]).then_with(|| {
                for basis_row in &basis {
                    let cmp = basis_row[left].total_cmp(&basis_row[right]);
                    if cmp != Ordering::Equal {
                        return cmp;
                    }
                }
                left.cmp(&right)
            })
        });

        let values = [raw[order[0]], raw[order[1]], raw[order[2]]];
        let mut vectors = [[0.0f64; 3]; 3];
        for (col, &src) in order.iter().enumerate() {
            for row in 0..3 {
                vectors[row][col] = basis[row][src];
            }
        }

        // Sign convention: first component of magnitude above SIGN_EPS is
        // positive. An eigenvector is only defined up to sign, so without this
        // two calls on numerically identical matrices can disagree.
        for col in 0..3 {
            let mut lead = 0usize;
            while lead < 2 && vectors[lead][col].abs() <= SIGN_EPS {
                lead += 1;
            }
            if vectors[lead][col] < 0.0 {
                for row in &mut vectors {
                    row[col] = -row[col];
                }
            }
        }

        (values, vectors)
    }

    /// Rebuild from eigenvalues and eigenvectors: `M = Σ_c λ_c v_c v_cᵀ`.
    ///
    /// `vectors[row][col]` is component `row` of eigenvector `col`, matching
    /// [`Sym3::eigen`]'s output exactly. The sum is symmetric for any input, but
    /// it reconstructs the intended matrix only for an orthonormal basis.
    pub(crate) fn from_eigen(values: [f64; 3], vectors: [[f64; 3]; 3]) -> Self {
        let mut out = [0.0f64; 6];
        for (col, &value) in values.iter().enumerate() {
            for i in 0..3 {
                for j in i..3 {
                    out[IDX[i][j]] += value * vectors[i][col] * vectors[j][col];
                }
            }
        }
        Self(out)
    }

    /// `|M|`: the same eigenvectors with `|eigenvalue|` in place of eigenvalue.
    ///
    /// This is the `|H_u|` of the interpolation-error bound. A negative
    /// curvature does not make an element cheaper, so only the magnitude
    /// survives into the metric.
    pub(crate) fn abs(&self) -> Self {
        self.spectral_map(f64::abs)
    }

    /// Matrix logarithm via the eigendecomposition. Requires strictly positive
    /// eigenvalues.
    ///
    /// Panics otherwise, and that is the intended behaviour: `log` of a matrix
    /// with a non-positive eigenvalue has no real value, and a metric out of
    /// [`metric_lp`] is floored at `D_LP · H_FLOOR > 0` and cannot reach the
    /// panic. Reaching it means the caller passed something that is not a
    /// metric.
    pub(crate) fn log(&self) -> Self {
        let (values, vectors) = self.eigen();
        assert!(
            values[0] > 0.0,
            "Sym3::log needs strictly positive eigenvalues; the smallest is {}, from {self:?}. \
             A metric from metric_lp is floored at H_FLOOR = {H_FLOOR:e} and cannot land here",
            values[0]
        );
        Self::from_eigen([values[0].ln(), values[1].ln(), values[2].ln()], vectors)
    }

    /// Matrix exponential via the eigendecomposition.
    ///
    /// Defined for every symmetric matrix and always positive-definite, which
    /// is what makes log-Euclidean interpolation stay inside the metric cone.
    pub(crate) fn exp(&self) -> Self {
        self.spectral_map(f64::exp)
    }

    /// `k · M`.
    pub(crate) fn scale(&self, k: f64) -> Self {
        let mut out = self.0;
        for entry in &mut out {
            *entry *= k;
        }
        Self(out)
    }

    /// `M + N`.
    pub(crate) fn add(&self, other: &Self) -> Self {
        let mut out = self.0;
        for (entry, &addend) in out.iter_mut().zip(other.0.iter()) {
            *entry += addend;
        }
        Self(out)
    }

    /// max eigenvalue / min eigenvalue, by magnitude; `INFINITY` if the min is
    /// zero.
    ///
    /// For a metric this is the anisotropy of the element it prescribes: the
    /// element's edge lengths go as `1/√λ`, so the length ratio is the square
    /// root of this. The eigenvalues come back ascending by *value*, not by
    /// magnitude, so both extremes are taken over the absolute values here.
    ///
    /// The zero matrix prescribes nothing and reports `INFINITY` rather than
    /// `NaN` — a degenerate metric is infinitely anisotropic, not undefined.
    pub(crate) fn aspect_ratio(&self) -> f64 {
        let (values, _) = self.eigen();
        let mut lo = f64::INFINITY;
        let mut hi = 0.0f64;
        for value in values {
            let magnitude = value.abs();
            if magnitude < lo {
                lo = magnitude;
            }
            if magnitude > hi {
                hi = magnitude;
            }
        }
        if lo > 0.0 { hi / lo } else { f64::INFINITY }
    }

    /// Apply a scalar function to the spectrum, keeping the eigenvectors.
    fn spectral_map(&self, mut map: impl FnMut(f64) -> f64) -> Self {
        let (values, vectors) = self.eigen();
        Self::from_eigen([map(values[0]), map(values[1]), map(values[2])], vectors)
    }
}

/// Central-difference Hessian of `sdf` at `p`, with step `h` (the cell size).
///
/// The diagonal uses the 7-point stencil
/// `(f(p + h eᵢ) − 2 f(p) + f(p − h eᵢ)) / h²`; the mixed entries use the
/// 4-point stencil
/// `(f⁺⁺ − f⁺⁻ − f⁻⁺ + f⁻⁻) / 4h²`. Nineteen samples in all.
///
/// # Exactness on a quadratic
///
/// Both stencils have **identically zero truncation error** for any quadratic:
/// the Taylor remainder starts at the fourth derivative for the diagonal and at
/// the third mixed derivative for the off-diagonal, and a quadratic has neither.
/// So on a quadratic the only error is the floating-point cancellation in the
/// subtraction, of relative size `≈ ε · |f| / (h² · |f''|)` — the numerator is
/// `O(f'' h²)` assembled from terms of size `O(f)`. That is the sense in which
/// this is exact to `f64` rounding, and it is why `h` should be the cell size
/// rather than something much smaller: shrinking `h` does not reduce a
/// truncation error that is already zero, it only amplifies the cancellation.
pub(crate) fn hessian<S: Sdf<Scalar = f64>>(sdf: &S, p: [f64; 3], h: f64) -> Sym3 {
    assert!(
        h > 0.0 && h.is_finite(),
        "hessian needs a positive finite step; got {h}"
    );

    let at = |offset: [f64; 3]| sdf.sample([p[0] + offset[0], p[1] + offset[1], p[2] + offset[2]]);
    let centre = at([0.0; 3]);
    let inv_h2 = 1.0 / (h * h);

    let mut out = [0.0f64; 6];

    for (axis, slot) in [(0usize, 0usize), (1, 3), (2, 5)] {
        let mut plus = [0.0; 3];
        plus[axis] = h;
        let mut minus = [0.0; 3];
        minus[axis] = -h;
        out[slot] = (at(plus) - 2.0 * centre + at(minus)) * inv_h2;
    }

    for (i, j, slot) in [(0usize, 1usize, 1usize), (0, 2, 2), (1, 2, 4)] {
        let corner = |si: f64, sj: f64| {
            let mut offset = [0.0; 3];
            offset[i] = si * h;
            offset[j] = sj * h;
            at(offset)
        };
        out[slot] = (corner(1.0, 1.0) - corner(1.0, -1.0) - corner(-1.0, 1.0) + corner(-1.0, -1.0))
            * (0.25 * inv_h2);
    }

    Sym3(out)
}

/// The optimal `L^p` metric: `D_Lp · det(|H|)^(−1/(2p + d)) · |H|`, `d = 3`.
///
/// `D_Lp` is [`D_LP`] `= 1` — a global constant that cancels in every ratio this
/// phase reports (module header). Each `|eigenvalue|` is floored at [`H_FLOOR`]
/// **before** the determinant is formed, so a flat direction produces a very
/// anisotropic metric rather than an infinity.
///
/// `p` must be positive. `p = f64::INFINITY` is allowed and is the `L^∞`
/// metric: the exponent goes to `−0`, `det^(−0) = 1`, and the metric is `|H|`
/// itself — the answer the theory gives in that limit, reached by the same
/// arithmetic rather than by a special case.
pub(crate) fn metric_lp(hessian: &Sym3, p: f64) -> Sym3 {
    assert!(
        p > 0.0,
        "metric_lp needs p > 0 (p = INFINITY is the L-infinity metric and is allowed); got {p}"
    );

    let (values, vectors) = hessian.eigen();
    let floored = [
        values[0].abs().max(H_FLOOR),
        values[1].abs().max(H_FLOOR),
        values[2].abs().max(H_FLOOR),
    ];
    let determinant = floored[0] * floored[1] * floored[2];
    let exponent = -1.0 / (2.0 * p + 3.0);
    let factor = D_LP * determinant.powf(exponent);

    Sym3::from_eigen(
        [
            floored[0] * factor,
            floored[1] * factor,
            floored[2] * factor,
        ],
        vectors,
    )
}

/// Complexity `C(M) = ∫ √det M` over the domain, as a Riemann sum on the sample
/// grid. This stands in for the vertex count.
///
/// `√det M` is a point density — a unit-metric simplex occupies Euclidean
/// volume `1/√det M` — so the integral is the number of unit elements the metric
/// asks for, and it is the budget `D_Lp` would be fitted to pin.
///
/// Panics on a non-positive determinant. That is not defensive coding: a metric
/// is positive-definite, [`metric_lp`] guarantees it through [`H_FLOOR`], and a
/// non-positive determinant here means the caller assembled something that is
/// not a metric field. Clamping at zero and carrying on would report a
/// complexity for a mesh that cannot exist.
pub(crate) fn complexity(metrics: &[Sym3], cell_volume: f64) -> f64 {
    assert!(
        cell_volume > 0.0 && cell_volume.is_finite(),
        "complexity needs a positive finite cell volume; got {cell_volume}"
    );

    let mut sum = 0.0f64;
    for (index, metric) in metrics.iter().enumerate() {
        let determinant = metric.det();
        assert!(
            determinant > 0.0,
            "complexity: metric {index} has det {determinant} <= 0, so it is not a metric; \
             metric_lp floors at H_FLOOR = {H_FLOOR:e} and cannot produce this"
        );
        sum += determinant.sqrt();
    }
    sum * cell_volume
}

/// Component-wise linear interpolation of two metrics.
///
/// `(1 − t) A + t B`, entry by entry. Positive-definite for `t ∈ [0, 1]` because
/// the cone is convex, but **not** determinant-monotone: interpolating between
/// two metrics whose anisotropies point in different directions swells
/// `det`, which in mesh terms means the seam quietly asks for coarser elements
/// than either side did. Measuring that swell is P-148 C1.
pub(crate) fn interp_componentwise(a: &Sym3, b: &Sym3, t: f64) -> Sym3 {
    let mut out = [0.0f64; 6];
    for (slot, (&left, &right)) in out.iter_mut().zip(a.0.iter().zip(b.0.iter())) {
        *slot = (1.0 - t) * left + t * right;
    }
    Sym3(out)
}

/// Log-Euclidean interpolation, `exp((1 − t) log A + t log B)`.
/// Determinant-monotone.
///
/// Monotone because `log det` is linear along this path:
/// `log det exp(X) = tr X`, so
/// `log det M(t) = (1 − t) log det A + t log det B` exactly, and the determinant
/// is the geometric interpolant of the two endpoints' determinants. It can
/// therefore never exceed both, which is the swell [`interp_componentwise`]
/// admits.
///
/// Both arguments must be positive-definite; see [`Sym3::log`].
pub(crate) fn interp_log_euclidean(a: &Sym3, b: &Sym3, t: f64) -> Sym3 {
    a.log().scale(1.0 - t).add(&b.log().scale(t)).exp()
}

/// max/min `|eigenvalue|` of the metric, i.e. the anisotropy of the element it
/// prescribes.
///
/// Free-function form of [`Sym3::aspect_ratio`], for the consumers that map it
/// over a slice. Read the module header's warning about the floor before
/// reporting a maximum of this.
pub(crate) fn aspect_ratio(m: &Sym3) -> f64 {
    m.aspect_ratio()
}

/// The AM–GM gap: `‖√|det H|‖ / ‖tr|H| / d‖` over a population of Hessians,
/// both in the `l^τ` norm named by `tau`. This is the quantity P-147 C2
/// correlates against.
///
/// Per Hessian the numerator's term is `√|det H| = √(Π|λᵢ|)` and the
/// denominator's is `tr|H| / d = (Σ|λᵢ|) / 3`, the arithmetic mean of the
/// curvature magnitudes. Both come from one eigendecomposition, so `tr|H|` here
/// is `Σ|λ|` and not `|tr H|` — those differ on every saddle and the bound is
/// about the former.
///
/// The ratio is the factor by which the anisotropic error constant beats the
/// isotropic one: by AM–GM the geometric mean never exceeds the arithmetic
/// mean, and the gap is largest exactly where one curvature is small.
///
/// # A fidelity choice, stated
///
/// `√|det H|` is `λ^{3/2}` in three dimensions while `tr|H|/d` is `λ^1`, so the
/// ratio is **not** dimensionless here — it carries one power of a curvature.
/// The literal `‖√|det H|‖ / ‖tr|H|/d‖` is what P-147's registration names, and
/// it is homogeneous in the two-dimensional form of the bound it was quoted
/// from. It is implemented exactly as registered rather than silently repaired
/// to `det^{1/d}`, because P-147 C2 asks only for a **rank correlation** across
/// fields, which is invariant under any monotone reparametrisation of the gap,
/// and because amending a registration once a run exists is forbidden. A
/// consumer must not compare this number across different length units.
///
/// `tau = f64::INFINITY` gives the max norm. Panics on an empty population —
/// a norm over nothing is not zero, it is unasked.
pub(crate) fn am_gm_gap(hessians: &[Sym3], tau: f64) -> f64 {
    assert!(
        !hessians.is_empty(),
        "VOID: am_gm_gap over an empty population of Hessians"
    );
    assert!(tau > 0.0, "am_gm_gap needs tau > 0; got {tau}");

    let mut geometric = Vec::with_capacity(hessians.len());
    let mut arithmetic = Vec::with_capacity(hessians.len());
    for hess in hessians {
        let (values, _) = hess.eigen();
        let determinant = values[0] * values[1] * values[2];
        geometric.push(determinant.abs().sqrt());
        arithmetic.push((values[0].abs() + values[1].abs() + values[2].abs()) / 3.0);
    }

    let numerator = l_tau_norm(&geometric, tau);
    let denominator = l_tau_norm(&arithmetic, tau);
    if denominator > 0.0 {
        numerator / denominator
    } else {
        f64::INFINITY
    }
}

/// `(Σ |xᵢ|^τ)^{1/τ}`, or the max norm for `τ = ∞`. Summed in slice order, so
/// the answer is reproducible bit for bit.
fn l_tau_norm(xs: &[f64], tau: f64) -> f64 {
    if tau.is_infinite() {
        return xs.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));
    }
    let mut sum = 0.0f64;
    for x in xs {
        sum += x.abs().powf(tau);
    }
    sum.powf(1.0 / tau)
}

/// The two principal curvatures of the level set of `sdf` at `p`, ascending by
/// magnitude, from the shape operator restricted to the tangent plane. Returns
/// `None` where the gradient is degenerate.
///
/// # How it is computed
///
/// With `g = ∇f` and `n = g/‖g‖` the outward unit normal, the second
/// fundamental form of the level set is `II(X, Y) = ⟨∇_X n, Y⟩ = Xᵀ H Y / ‖g‖`
/// for tangent `X`, `Y` — the projector `I − n nᵀ` drops out because `X` and `Y`
/// are already tangent. So the method builds an orthonormal tangent pair
/// `(t₁, t₂)`, forms the `2 × 2` matrix `Bab = tₐᵀ H t_b / ‖g‖`, and takes its
/// two eigenvalues in closed form:
/// `κ = mean ∓ √(((B₀₀ − B₁₁)/2)² + B₀₁²)` with `mean = (B₀₀ + B₁₁)/2`.
///
/// Sign convention: positive where the surface curves **away** from the outward
/// normal, so the unit sphere `‖x‖ − 1` gives `κ = (1, 1)` and a saddle gives
/// one of each sign. A plane gives `(0, 0)`, a cylinder of radius `R` gives
/// `(0, 1/R)` — that zero is the flat direction Group D exists to exploit, and
/// P-149's `principal_curvature_ratio` is `|κ₀| / |κ₁|`.
///
/// The gradient is differenced at the **same** step `h` as the Hessian rather
/// than through `Sdf::gradient`, whose default step is
/// `Real::DIFF_STEP · max(|pᵢ|, 1)` and therefore a different discrete object.
/// `κ` is a ratio of the two; measuring numerator and denominator at two
/// different steps would put a step mismatch into the curvature ratio, which is
/// exactly the column P-149 correlates against.
///
/// The tangent basis is chosen deterministically — cross `n` with the axis in
/// which `n` is smallest, first index on a tie. That axis satisfies
/// `|n_axis| ≤ 1/√3`, so the cross product has length at least `√(2/3)` and
/// cannot degenerate; there is no second branch to take.
pub(crate) fn principal_curvatures<S: Sdf<Scalar = f64>>(
    sdf: &S,
    p: [f64; 3],
    h: f64,
) -> Option<[f64; 2]> {
    assert!(
        h > 0.0 && h.is_finite(),
        "principal_curvatures needs a positive finite step; got {h}"
    );

    let mut gradient = [0.0f64; 3];
    for (axis, slot) in [(0usize, 0usize), (1, 1), (2, 2)] {
        let mut plus = p;
        plus[axis] += h;
        let mut minus = p;
        minus[axis] -= h;
        gradient[slot] = (sdf.sample(plus) - sdf.sample(minus)) / (2.0 * h);
    }

    let norm =
        (gradient[0] * gradient[0] + gradient[1] * gradient[1] + gradient[2] * gradient[2]).sqrt();
    if !norm.is_finite() || norm <= GRAD_FLOOR {
        return None;
    }
    let normal = [gradient[0] / norm, gradient[1] / norm, gradient[2] / norm];

    let mut smallest = 0usize;
    for axis in [1usize, 2] {
        if normal[axis].abs() < normal[smallest].abs() {
            smallest = axis;
        }
    }
    let mut axis_vector = [0.0f64; 3];
    axis_vector[smallest] = 1.0;

    let raw_t1 = cross(normal, axis_vector);
    let t1_norm = (raw_t1[0] * raw_t1[0] + raw_t1[1] * raw_t1[1] + raw_t1[2] * raw_t1[2]).sqrt();
    let t1 = [
        raw_t1[0] / t1_norm,
        raw_t1[1] / t1_norm,
        raw_t1[2] / t1_norm,
    ];
    let t2 = cross(normal, t1);

    let hess = hessian(sdf, p, h);
    let b00 = quadratic_form(&hess, t1, t1) / norm;
    let b01 = quadratic_form(&hess, t1, t2) / norm;
    let b11 = quadratic_form(&hess, t2, t2) / norm;

    let mean = 0.5 * (b00 + b11);
    let half_gap = 0.5 * (b00 - b11);
    let radius = (half_gap * half_gap + b01 * b01).sqrt();
    let mut curvatures = [mean - radius, mean + radius];
    if curvatures[1].abs() < curvatures[0].abs() {
        curvatures.swap(0, 1);
    }
    Some(curvatures)
}

/// A metric-driven sampling density field: per cell, the number of samples the
/// metric prescribes relative to a uniform grid of the same total budget.
/// Deterministic. Returns one multiplier per cell in x-fastest order.
///
/// The metric's point density is `√det M` (see [`complexity`]); a uniform grid
/// spends the same budget evenly, so it spends the population mean everywhere.
/// The multiplier is therefore `√det Mᵢ / mean(√det M)`, whose mean is exactly
/// one — the two arms are matched on total budget by construction rather than by
/// a fitted normalisation, which is what makes a triangle-count comparison
/// between them mean anything.
///
/// Input and output are in the caller's own cell order and the sum is taken in
/// slice order, so the result is reproducible bit for bit. An empty input gives
/// an empty output. A non-positive determinant panics, for the same reason as in
/// [`complexity`].
pub(crate) fn density_from_metric(metrics: &[Sym3]) -> Vec<f64> {
    let mut densities = Vec::with_capacity(metrics.len());
    let mut sum = 0.0f64;
    for (index, metric) in metrics.iter().enumerate() {
        let determinant = metric.det();
        assert!(
            determinant > 0.0,
            "density_from_metric: cell {index} has det {determinant} <= 0, so it is not a \
             metric; metric_lp floors at H_FLOOR = {H_FLOOR:e} and cannot produce this"
        );
        let density = determinant.sqrt();
        densities.push(density);
        sum += density;
    }
    if densities.is_empty() {
        return densities;
    }

    let mean = sum / densities.len() as f64;
    for density in &mut densities {
        *density /= mean;
    }
    densities
}

/// `a × b`.
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// `uᵀ M v`.
fn quadratic_form(m: &Sym3, u: [f64; 3], v: [f64; 3]) -> f64 {
    let mut total = 0.0f64;
    for (i, &ui) in u.iter().enumerate() {
        for (j, &vj) in v.iter().enumerate() {
            total += ui * m.get(i, j) * vj;
        }
    }
    total
}
