//! **P-54 — a revised affine form against Hart's Lipschitz ball.**
//!
//! Ticket: R-049. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p54
//! ```
//!
//! Writes `docs/experiments/p-54.csv`.
//!
//! # The two rejection tests, side by side
//!
//! **Hart's**, which the crate already has, transcribed from
//! `subgrid/extract.rs::cell_is_provably_empty`: one sample at the cell centre,
//! rejected when `|f(centre)| > l · r` with `r = (√3/2)·h` the circumradius of
//! the cell and `l` the field's declared Lipschitz constant. Strict, so a value
//! exactly on the bound subdivides. It reasons over the **ball** that
//! circumscribes the cell, and it knows nothing about `f` except one value and
//! one constant.
//!
//! **The revised affine form** (Fryazinov, Pasko & Comninos,
//! `10.1016/j.cag.2010.07.003`) reasons over the **box** and knows the whole
//! expression. `x̂ = x0 + x1·e1 + x2·e2 + x3·e3 + ex·[-1,1]` with `ex ≥ 0`: five
//! stored `f64`, fixed size, **never growing**, no heap, because every non-affine
//! error accumulates into the trailing `ex` rather than opening a new symbol. A
//! cell box maps in as `centre_k + half_extent_k · e_k` with `ex = 0`, so the
//! three noise symbols *are* the three axes and stay that way through the whole
//! expression. The range is `[x0 − rad, x0 + rad]` with
//! `rad = |x1| + |x2| + |x3| + ex`, and the cell is provably empty exactly when
//! that range excludes zero.
//!
//! The mechanism is correlation. `sin(x)cos(y) + sin(y)cos(z) + sin(z)cos(x)`
//! cannot have all three terms extremal at once; a per-term interval bound throws
//! that away, and so does a ball of radius `(√3/2)·h` around one sample.
//!
//! # `min`/`max` are not in the source paper, and that is the point
//!
//! §4.4 of Fryazinov et al. concludes that general conditionals are an open
//! problem, and the paper gives **no** affine rule for `min`, `max` or `abs`. The
//! only sound treatment available is therefore to collapse both operands to their
//! ranges, take the interval `min`/`max` — which is exact on *independent*
//! intervals and merely sound on correlated ones — and hand back a **degenerate**
//! form carrying the whole width in `ex`. Every noise symbol is gone at that
//! point, so downstream operations have nothing left to cancel.
//!
//! That is not a fallback path bolted onto a primary one; it is the only rule
//! there is, applied unconditionally, and it is what makes the prediction
//! non-uniform: a field built from smooth trig keeps its correlation and a field
//! built from `min`/`max` loses it at the first combinator.
//!
//! `abs` is the one place a shortcut is taken, and it is exact rather than
//! opportunistic: when the form's range does not straddle zero, `|x̂|` **is**
//! `±x̂`, with the sign decided by a range that has already been computed. Only a
//! straddling range collapses. `box_exact` reaches the exact branch on every cell
//! that does not meet a coordinate plane, so this is not where its correlation
//! goes — the `max(q, 0)` and `min(max q, 0)` that follow are.
//!
//! # Chebyshev, not min-range, and why
//!
//! `sin`, `cos` and `sqrt` are not affine, so each is replaced by
//! `α·x̂ + ζ ± δ`. Two standard choices for `α`:
//!
//! - **min-range** takes `α` from an endpoint derivative, which makes the
//!   *resulting range* exactly the function's true range on a monotone branch,
//!   and pays for it with a larger `δ`.
//! - **Chebyshev** takes `α` as the secant slope, the best `L∞` linear fit, which
//!   minimises `δ` and pays with a resulting range slightly wider than the truth.
//!
//! This uses **Chebyshev**, deliberately. `δ` lands in `ex`, and `ex` is exactly
//! the part that **cannot cancel** against anything: two operands' `e_k` terms
//! subtract, two operands' `ex` terms add. In a one-operation expression
//! min-range wins and the choice would be the other way; in a sum of products of
//! shared arguments — which is every field here that the prediction is about —
//! the quantity to minimise is the incompressible one.
//!
//! `δ` is computed rather than bounded away: the residual `f(t) − α·t` is
//! extremal at an endpoint or where `f'(t) = α`, so those critical points are
//! enumerated (`cos t = α` for `sin`, `sin t = −α` for `cos`, one interior point
//! for `sqrt`) and the residual is evaluated at each. The containment test uses a
//! guard that errs toward *including* a critical point, which can only widen.
//!
//! # Soundness is checked, not asserted
//!
//! An unsound bound rejects a cell that contains surface, which is a **hole** —
//! and a hole is invisible to every validity gate this crate has, since the mesh
//! is simply missing a piece and remains perfectly manifold. So soundness is
//! measured twice, and both numbers are on every CSV row:
//!
//! 1. [`verify_substitutions`] checks the affine substitution itself, per noise
//!    symbol rather than merely per range: for a single-symbol input `x̂`, every
//!    `ε ∈ [-1, 1]` must satisfy `|f(x0 + x1·ε) − (y0 + y1·ε)| ≤ ex_y`. That is
//!    the property composition depends on, and a range check would not see it
//!    fail.
//! 2. `soundness_violations` checks the whole pipeline against the crate's own
//!    `Sdf::sample` at `9³` points inside every cell of every field. The affine
//!    range must contain all of them; `soundness_min_margin` is the tightest it
//!    ever came, in field units.
//!
//! Both slacks are `256 · f64::EPSILON · max(1, |value|)` per operation, P-48's
//! constant for P-48's reason: `core` has no directed rounding, and the enclosure
//! has to cover the *field's* own `f64` rounding as well as its own. Every
//! transcendental below is the same `libm` call `real.rs` binds for `f64`.
//!
//! # `gyroid` is a `max`, and the registration's baseline is that `max`
//!
//! The registration predicts C1 on "gyroid … built only from smooth trig with
//! shared arguments, NO min/max", against a Lipschitz baseline of 688 of 4,096.
//! Those two are different fields. The crate's reference field named `gyroid` is
//! [`capped_gyroid`] — `Intersection<Gyroid, Sphere>`, a `max` — over `[-7, 7]³`,
//! and it is the one that rejects 688. The bare [`Gyroid`], the pure-trig field
//! the mechanism is about, rejects **zero** cells on the same grid, because
//! `l · r = 2√3 · (√3/2) · 0.875 = 2.625` exceeds the gyroid's entire range of
//! `±1.5`: Hart's test cannot reject a single cell of it at any resolution where
//! the cell is that large.
//!
//! Both are therefore rows. `gyroid` is the registered baseline, unaltered.
//! `gyroid_uncapped` is the same `Gyroid::canonical()` on the same domain with
//! the same declared constant, minus the cap, and it is where the correlation
//! claim can actually be read. Neither is a substitute for the other and the CSV
//! carries both.
//!
//! ## Is it the constant or the ball? `rejected_at_gradient_sup` decides it
//!
//! The registration's falsification clause reads C1 failing as evidence that
//! "M-267's 2× gap is genuinely attainable by the gradient rather than an
//! artefact of the ball". That is answerable directly rather than by inference,
//! so it is a column: Hart's test re-run with the **measured** supremum `1.731`
//! in place of the declared `2√3`.
//!
//! It is a counterfactual and not an option. A sampled maximum is a *lower*
//! bound on the true supremum, so declaring it would be unsound in precisely the
//! direction `validate/field_bound.rs` refuses to err in — the column is never
//! gated against a mesh and nothing here suggests using it.
//!
//! # What the mesh columns mean
//!
//! C3 is the gate. Two independent checks, both on every row:
//!
//! - **`baseline_vertices_in_rejected_cells`** is the whole-grid witness. Every
//!   vertex the subgrid extractor emits from a cell lies in that cell's closed
//!   box, so if no vertex of the **ungated** whole-grid mesh lies in any rejected
//!   cell's closed box, skipping those cells removes nothing. It must be zero.
//! - **`mesh_hash_rig_gated`** vs **`mesh_hash_rig_ungated`** is the byte
//!   comparison. The extractor exposes no per-cell hook, so the rig extracts each
//!   cell on its own `2³` grid and appends; gating is then simply not appending.
//!   `mesh_identical` is the conjunction of the two checks.
//!
//! `mesh_hash_lipschitz_gated` and `mesh_hash_ungated` are the crate's own
//! whole-grid hashes with and without `set_lipschitz`, recorded so the rig can be
//! seen to agree with the path the crate actually ships: they are equal, and
//! `rig_triangles` equals `wholegrid_triangles` on every field.
//!
//! ## And the instrument is checked against a known-bad rejection
//!
//! `mesh_identical` reading `true` on twelve rows proves nothing unless a wrong
//! rejection would have made it read `false`. So every field is gated a third
//! way, by Hart's test with the declared constant divided by
//! [`CONTROL_DIVISOR`] — M-244's incident, where a hand-reasoned constant was
//! wrong by 3× on the first try, made deliberate. `control_detected` must be
//! `true` on every row; a `false` there would mean the rig cannot see a hole and
//! would void C3 rather than confirm it.
//!
//! # Counted, not timed
//!
//! `rejected_cells` is an integer and identical on every machine.
//! `bound_ns_per_cell` sits beside it and **gates nothing** — M-348 is the
//! incident where a discovery was demoted for resting on a wall clock.

mod common;

use std::ops::{Add, Mul, Neg, Sub};
use std::time::Instant;

use common::experiment::Run;
use isomesh::fields::{
    BoundedSdf, BoxExact, Difference, Gyroid, Intersection, ReferenceField, Sphere, Torus,
    capped_gyroid, csg_difference,
};
use isomesh::subgrid::extract::SubgridMarchingTetrahedra;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// The registered resolution: 17 samples, so 16 cells per axis and 4,096 cells.
const SAMPLES: u32 = 17;

/// Cells per axis.
const CELLS_PER_AXIS: u32 = SAMPLES - 1;

/// Cells in the grid. The registration's denominator.
const CELLS: u32 = CELLS_PER_AXIS * CELLS_PER_AXIS * CELLS_PER_AXIS;

/// Samples per tetrahedron edge for the subgrid extractor, matching the crate's
/// own `rejection_does_not_change_the_mesh`.
const EDGE_SAMPLES: u32 = 6;

/// Dense-sampling resolution per cell axis for the soundness check: `9³ = 729`
/// points in every cell, endpoints exactly on the cell faces.
const DENSE: u32 = 9;

/// Samples per interval for [`verify_substitutions`].
const PROBE_SAMPLES: u32 = 4097;

/// ULP-scale units of slack per operation, for `core`'s lack of directed
/// rounding *and* for the field implementation's own `f64` error. P-48's
/// constant, for P-48's reason.
const SLACK_ULPS: f64 = 256.0;

/// Absolute guard for the critical-point containment test, in the argument's
/// units. Erring toward inclusion only ever widens.
const CRITICAL_GUARD: f64 = 1e-12;

/// Half the space diagonal of a unit cube: the circumradius factor in
/// `cell_is_provably_empty`, transcribed to the digit.
const CIRCUMRADIUS: f64 = 1.732_050_807_568_877_2;

/// The negative control divides the declared Lipschitz constant by this.
///
/// **The instrument has to be able to report the bad news.** `mesh_identical`
/// reading `true` on every row is worthless unless a wrong rejection would make
/// it read `false`, so every field is also gated by Hart's test with a constant
/// four times too small — M-244's incident exactly, where a hand-reasoned
/// constant was wrong by 3× on the first try. `control_mesh_identical` is
/// expected `false`, and a `true` there voids C3 on that field rather than
/// confirming it.
const CONTROL_DIVISOR: f64 = 4.0;

const PI: f64 = std::f64::consts::PI;
const TAU: f64 = std::f64::consts::TAU;

// ─── the revised affine form ────────────────────────────────────────────────

/// `x0 + x1·e1 + x2·e2 + x3·e3 + ex·[-1, 1]`, with `ex ≥ 0`.
///
/// Five `f64`. The size is the invariant: no operation here opens a fourth noise
/// symbol, so a hundred-operation expression costs the same forty bytes as one.
#[derive(Clone, Copy)]
struct Af {
    /// The centre, `x0`.
    c: f64,
    /// Coefficients of the three noise symbols, one per cell axis.
    e: [f64; 3],
    /// The accumulator. Every non-affine error ever committed, summed.
    ex: f64,
}

impl Af {
    /// A constant: no symbols, no error.
    const fn point(v: f64) -> Self {
        Self {
            c: v,
            e: [0.0; 3],
            ex: 0.0,
        }
    }

    /// One cell axis: `centre + half·e_axis`, the exact image of the cell box.
    fn symbol(centre: f64, half: f64, axis: usize) -> Self {
        let mut e = [0.0; 3];
        e[axis] = half;
        Self {
            c: centre,
            e,
            ex: 0.0,
        }
    }

    /// A degenerate form covering `[lo, hi]`, with the width in `ex`.
    ///
    /// What `min`, `max` and a straddling `abs` return, and the only shape in
    /// which correlation is lost.
    fn interval(lo: f64, hi: f64) -> Self {
        let c = (lo + hi) * 0.5;
        let ex = (hi - lo) * 0.5;
        Self {
            c,
            e: [0.0; 3],
            ex: ex.max(0.0),
        }
        .widen(lo.abs().max(hi.abs()))
    }

    /// `|x1| + |x2| + |x3|`: the width the noise symbols alone can reach.
    fn linear(&self) -> f64 {
        self.e[0].abs() + self.e[1].abs() + self.e[2].abs()
    }

    /// `|x1| + |x2| + |x3| + ex`: half the width of the range.
    fn radius(&self) -> f64 {
        self.linear() + self.ex
    }

    fn lo(&self) -> f64 {
        self.c - self.radius()
    }

    fn hi(&self) -> f64 {
        self.c + self.radius()
    }

    /// The largest magnitude any value of the form can reach. The scale the
    /// per-operation slack is measured against.
    fn magnitude(&self) -> f64 {
        self.c.abs() + self.radius()
    }

    /// `0 ∉ [lo, hi]`. The certificate: the cell is provably empty.
    fn excludes_zero(&self) -> bool {
        self.lo() > 0.0 || self.hi() < 0.0
    }

    /// Slack for one rounded operation on quantities of magnitude `m`.
    fn slack(m: f64) -> f64 {
        SLACK_ULPS * f64::EPSILON * m.max(1.0)
    }

    /// Absorb one operation's rounding into the accumulator.
    fn widen(mut self, m: f64) -> Self {
        self.ex += Self::slack(m);
        self
    }

    /// Multiply by a constant. Exact in the symbols, so correlation survives.
    fn scaled(self, k: f64) -> Self {
        Self {
            c: self.c * k,
            e: [self.e[0] * k, self.e[1] * k, self.e[2] * k],
            ex: self.ex * k.abs(),
        }
    }

    /// Add a constant. The one operation with no error at all.
    fn shifted(self, k: f64) -> Self {
        Self {
            c: self.c + k,
            ..self
        }
    }

    /// The square, tighter than `self * self` and for a reason that matters.
    ///
    /// Writing `x̂ = x0 + L + E`, the square is
    /// `x0² + 2x0·L + 2x0·E + (L + E)²`, and `(L + E)²` lies in `[0, rad²]`
    /// **exactly** — it is a square. So the last term contributes `rad²/2` to the
    /// centre and `rad²/2` to `ex`, where [`Mul`] would contribute `rad²` to `ex`
    /// alone and leave the centre at `x0²`, which is twice as wide *and*
    /// mis-centred. Every distance field here is a sum of squares under a root,
    /// so this is not a micro-optimisation.
    fn sqr(self) -> Self {
        let r = self.radius();
        let half = r * r * 0.5;
        let out = Self {
            c: self.c * self.c + half,
            e: [
                2.0 * self.c * self.e[0],
                2.0 * self.c * self.e[1],
                2.0 * self.c * self.e[2],
            ],
            ex: 2.0 * self.c.abs() * self.ex + half,
        };
        let m = out.magnitude();
        out.widen(m)
    }

    /// Absolute value.
    ///
    /// Exact and correlation-preserving unless the range straddles zero, in which
    /// case there is nothing to preserve: `|x̂|` is then genuinely not affine in
    /// the symbols and the exact interval `[0, max(−lo, hi)]` is the most that can
    /// be said. Not in the source paper, and treated as `min`/`max` is for the
    /// same reason.
    fn abs(self) -> Self {
        let (lo, hi) = (self.lo(), self.hi());
        if lo >= 0.0 {
            self
        } else if hi <= 0.0 {
            -self
        } else {
            Self::interval(0.0, (-lo).max(hi))
        }
    }

    /// `min`, by collapsing both operands to their ranges.
    ///
    /// Sound, exact on independent operands, and the end of every noise symbol in
    /// the expression. See the module docs: the paper offers nothing else.
    fn amin(self, o: Self) -> Self {
        Self::interval(self.lo().min(o.lo()), self.hi().min(o.hi()))
    }

    /// `max`, by the same collapse.
    fn amax(self, o: Self) -> Self {
        Self::interval(self.lo().max(o.lo()), self.hi().max(o.hi()))
    }

    /// Replace `f(x̂)` by `α·x̂ + ζ ± δ`, with `δ` bounding the residual
    /// `f(t) − α·t` over the form's own range.
    ///
    /// Sound because the form's values sweep exactly `[lo, hi]` — it is a
    /// continuous linear image of a box — so bounding the residual there bounds it
    /// at every `(ε, η)`. The new `δ` merges into `ex` rather than opening a
    /// symbol, which over-approximates by treating it as independent of the old
    /// `ex`, and over-approximating is the safe direction.
    fn substitute(self, alpha: f64, f: impl Fn(f64) -> f64, criticals: &[f64]) -> Self {
        let (lo, hi) = (self.lo(), self.hi());
        let mut rlo = f(lo) - alpha * lo;
        let mut rhi = rlo;
        let end = f(hi) - alpha * hi;
        rlo = rlo.min(end);
        rhi = rhi.max(end);
        for &t in criticals {
            // Guarded toward inclusion: a critical point just outside the range
            // contributes a residual that can only be larger, which widens.
            if t > lo - CRITICAL_GUARD && t < hi + CRITICAL_GUARD {
                let r = f(t) - alpha * t;
                rlo = rlo.min(r);
                rhi = rhi.max(r);
            }
        }
        let zeta = (rlo + rhi) * 0.5;
        let delta = ((rhi - rlo) * 0.5).max(0.0);
        let mut out = self.scaled(alpha);
        out.c += zeta;
        out.ex += delta;
        let m = out.magnitude().max(rhi.abs()).max(rlo.abs());
        out.widen(m)
    }

    /// `sin` by Chebyshev substitution.
    ///
    /// The residual `sin(t) − α·t` is extremal at an endpoint or where
    /// `cos t = α`, which is `±acos(α) + 2kπ`. Four periods either side of the
    /// range's own are enumerated, which is three more than a range narrower than
    /// `2π` can need.
    fn sin(self) -> Self {
        let (lo, hi) = (self.lo(), self.hi());
        if hi - lo >= TAU {
            return Self::interval(-1.0, 1.0);
        }
        let alpha = if hi > lo {
            ((libm::sin(hi) - libm::sin(lo)) / (hi - lo)).clamp(-1.0, 1.0)
        } else {
            libm::cos(lo)
        };
        let ac = libm::acos(alpha);
        let k0 = libm::floor(lo / TAU);
        let mut criticals = [0.0; 8];
        for (n, dk) in [-1.0, 0.0, 1.0, 2.0].into_iter().enumerate() {
            let base = (k0 + dk) * TAU;
            criticals[2 * n] = ac + base;
            criticals[2 * n + 1] = -ac + base;
        }
        self.substitute(alpha, libm::sin, &criticals)
    }

    /// `cos` by Chebyshev substitution. Critical points solve `sin t = −α`.
    fn cos(self) -> Self {
        let (lo, hi) = (self.lo(), self.hi());
        if hi - lo >= TAU {
            return Self::interval(-1.0, 1.0);
        }
        let alpha = if hi > lo {
            ((libm::cos(hi) - libm::cos(lo)) / (hi - lo)).clamp(-1.0, 1.0)
        } else {
            -libm::sin(lo)
        };
        let asn = libm::asin((-alpha).clamp(-1.0, 1.0));
        let k0 = libm::floor(lo / TAU);
        let mut criticals = [0.0; 8];
        for (n, dk) in [-1.0, 0.0, 1.0, 2.0].into_iter().enumerate() {
            let base = (k0 + dk) * TAU;
            criticals[2 * n] = asn + base;
            criticals[2 * n + 1] = PI - asn + base;
        }
        self.substitute(alpha, libm::cos, &criticals)
    }

    /// `sqrt` by Chebyshev substitution over the non-negative part.
    ///
    /// Clamping the lower end at zero is exact rather than a fallback: every
    /// radicand here is a sum of squares, so the true value set is contained in
    /// `[max(lo, 0), hi]` whatever the enclosure's arithmetic did to the lower
    /// end. `√` is concave, so its single interior critical point for the secant
    /// slope is `((√lo + √hi)/2)²`.
    fn sqrt(self) -> Self {
        let hi = self.hi();
        if hi <= 0.0 {
            return Self::point(0.0);
        }
        let lo = self.lo().max(0.0);
        let (sl, sh) = (libm::sqrt(lo), libm::sqrt(hi));
        let alpha = 1.0 / (sl + sh);
        let mid = (sl + sh) * 0.5;
        // Correlation survives only when the enclosure never went negative; when
        // it did, the clamp has already discarded part of the form's own range
        // and the symbols no longer describe what is left.
        let base = if self.lo() >= 0.0 {
            self
        } else {
            Self::interval(lo, hi)
        };
        base.substitute(alpha, libm::sqrt, &[mid * mid])
    }
}

impl Add for Af {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        let out = Self {
            c: self.c + o.c,
            e: [self.e[0] + o.e[0], self.e[1] + o.e[1], self.e[2] + o.e[2]],
            ex: self.ex + o.ex,
        };
        let m = out.magnitude();
        out.widen(m)
    }
}

impl Sub for Af {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        let out = Self {
            c: self.c - o.c,
            e: [self.e[0] - o.e[0], self.e[1] - o.e[1], self.e[2] - o.e[2]],
            ex: self.ex + o.ex,
        };
        let m = out.magnitude();
        out.widen(m)
    }
}

impl Neg for Af {
    type Output = Self;
    /// Exact, and the only operation that adds nothing to `ex`.
    fn neg(self) -> Self {
        Self {
            c: -self.c,
            e: [-self.e[0], -self.e[1], -self.e[2]],
            ex: self.ex,
        }
    }
}

impl Mul for Af {
    type Output = Self;
    /// The tight rule: the linear part is the bilinear cross terms and every
    /// non-linear term goes into `ex`.
    ///
    /// With `x̂ = x0 + Lx + Ex` and `ŷ = y0 + Ly + Ey`, the product is
    /// `x0·y0 + (x0·Ly + y0·Lx) + x0·Ey + y0·Ex + (Lx + Ex)(Ly + Ey)`. The middle
    /// bracket is affine and keeps its symbols — that is where correlation
    /// survives a product. The rest is bounded by
    /// `|x0|·ey + |y0|·ex + (rx + ex)(ry + ey)`.
    fn mul(self, o: Self) -> Self {
        let (rx, ry) = (self.linear(), o.linear());
        let out = Self {
            c: self.c * o.c,
            e: [
                self.c * o.e[0] + o.c * self.e[0],
                self.c * o.e[1] + o.c * self.e[1],
                self.c * o.e[2] + o.c * self.e[2],
            ],
            ex: self.c.abs() * o.ex + o.c.abs() * self.ex + (rx + self.ex) * (ry + o.ex),
        };
        let m = out.magnitude();
        out.widen(m)
    }
}

// ─── the fields, over the affine form ───────────────────────────────────────

/// A field whose formula has been reimplemented over [`Af`].
///
/// Bench-local by ticket: `crates/isomesh/src/**` is read-only for R-049, so this
/// mirrors what P-48 did for interval arithmetic. Nothing checks the
/// transcription at compile time, which is why the soundness columns exist.
trait AffineForm {
    /// Whether the expression contains a `min`, a `max` or a straddling `abs` —
    /// the operations that collapse the form. C2's discriminator.
    const HAS_MIN_MAX: bool;

    /// The enclosure over the box the three symbols describe.
    fn affine(&self, p: &[Af; 3]) -> Af;
}

impl AffineForm for Sphere<f64> {
    const HAS_MIN_MAX: bool = false;
    fn affine(&self, p: &[Af; 3]) -> Af {
        let d = [
            p[0].shifted(-self.center[0]),
            p[1].shifted(-self.center[1]),
            p[2].shifted(-self.center[2]),
        ];
        ((d[0].sqr() + d[1].sqr()) + d[2].sqr())
            .sqrt()
            .shifted(-self.radius)
    }
}

impl AffineForm for Torus<f64> {
    const HAS_MIN_MAX: bool = false;
    fn affine(&self, p: &[Af; 3]) -> Af {
        let d = [
            p[0].shifted(-self.center[0]),
            p[1].shifted(-self.center[1]),
            p[2].shifted(-self.center[2]),
        ];
        let s = (d[0].sqr() + d[2].sqr()).sqrt();
        let q0 = s.shifted(-self.major);
        (q0.sqr() + d[1].sqr()).sqrt().shifted(-self.minor)
    }
}

impl AffineForm for BoxExact<f64> {
    const HAS_MIN_MAX: bool = true;
    fn affine(&self, p: &[Af; 3]) -> Af {
        // `box_sample` verbatim: q = |p − c| − b, then ‖max(q, 0)‖ + min(max q, 0).
        let q = [
            p[0].shifted(-self.center[0])
                .abs()
                .shifted(-self.half_extents[0]),
            p[1].shifted(-self.center[1])
                .abs()
                .shifted(-self.half_extents[1]),
            p[2].shifted(-self.center[2])
                .abs()
                .shifted(-self.half_extents[2]),
        ];
        let zero = Af::point(0.0);
        let outside =
            ((q[0].amax(zero).sqr() + q[1].amax(zero).sqr()) + q[2].amax(zero).sqr()).sqrt();
        let inside = q[0].amax(q[1]).amax(q[2]).amin(zero);
        outside + inside
    }
}

impl AffineForm for Gyroid<f64> {
    const HAS_MIN_MAX: bool = false;
    fn affine(&self, p: &[Af; 3]) -> Af {
        let a = p[0].scaled(self.scale);
        let b = p[1].scaled(self.scale);
        let c = p[2].scaled(self.scale);
        (((a.sin() * b.cos()) + (b.sin() * c.cos())) + (c.sin() * a.cos())).shifted(-self.iso)
    }
}

impl<A: AffineForm, B: AffineForm> AffineForm for Intersection<A, B> {
    const HAS_MIN_MAX: bool = true;
    fn affine(&self, p: &[Af; 3]) -> Af {
        self.a.affine(p).amax(self.b.affine(p))
    }
}

impl<A: AffineForm, B: AffineForm> AffineForm for Difference<A, B> {
    const HAS_MIN_MAX: bool = true;
    fn affine(&self, p: &[Af; 3]) -> Af {
        self.a.affine(p).amax(-self.b.affine(p))
    }
}

// ─── soundness of the substitutions ─────────────────────────────────────────

/// What [`verify_substitutions`] found. Carried onto every CSV row.
struct Probe {
    samples: u64,
    violations: u64,
    /// The tightest the slack ever came, in the function's own units. Negative
    /// would mean an unsound substitution and would void the whole experiment.
    min_margin: f64,
}

/// **Is the affine substitution conservative? Checked per noise symbol.**
///
/// The property composition rests on is stronger than range containment: for a
/// single-symbol input `x̂ = x0 + x1·e`, the output `ŷ = y0 + y1·e + ey·[-1,1]`
/// must satisfy `|f(x0 + x1·ε) − (y0 + y1·ε)| ≤ ey` for **every** `ε ∈ [-1, 1]`,
/// not merely `f(x) ∈ [ŷ.lo, ŷ.hi]`. A substitution can pass the range check and
/// fail this one, and it is this one that a later `mul` relies on.
///
/// Centres and radii are a fixed deterministic sweep, chosen off any multiple of
/// `π` and spanning from far narrower to far wider than a cell of this
/// experiment: `sin` and `cos` over radii up to `7`, `sqrt` over radicands whose
/// enclosure stays non-negative, so the clamp is not what is being measured.
fn verify_substitutions() -> Probe {
    let mut samples = 0u64;
    let mut violations = 0u64;
    let mut min_margin = f64::INFINITY;

    let radii = [0.001, 0.01, 0.1, 0.437_5, 0.875, 1.5, 3.0, 6.0, 7.0];

    let mut check = |x: Af, y: Af, f: &dyn Fn(f64) -> f64| {
        for i in 0..PROBE_SAMPLES {
            let eps = 2.0 * f64::from(i) / f64::from(PROBE_SAMPLES - 1) - 1.0;
            let t = x.c + x.e[0] * eps;
            let predicted = y.c + y.e[0] * eps;
            let margin = y.ex - (f(t) - predicted).abs();
            if margin < 0.0 {
                violations += 1;
            }
            min_margin = min_margin.min(margin);
            samples += 1;
        }
    };

    for i in 0..137 {
        // 0.37 is incommensurate with π to well past f64's reach, so no centre
        // lands on a critical point and the enumeration is exercised rather than
        // sidestepped.
        let centre = -25.0 + 0.37 * f64::from(i);
        for &r in &radii {
            let x = Af::symbol(centre, r, 0);
            check(x, x.sin(), &libm::sin);
            check(x, x.cos(), &libm::cos);
            // Only radicands whose enclosure is entirely non-negative: the clamp
            // is a separate, exact argument and would otherwise be measured here.
            if centre - r > 0.0 {
                check(x, x.sqrt(), &libm::sqrt);
            }
        }
    }

    Probe {
        samples,
        violations,
        min_margin,
    }
}

// ─── the grid ───────────────────────────────────────────────────────────────

/// Flat index of a cell, `x` fastest, matching the extractor's own loop order.
fn cell_index(i: u32, j: u32, k: u32) -> usize {
    ((k * CELLS_PER_AXIS + j) * CELLS_PER_AXIS + i) as usize
}

/// The two rejected-cell sets, and the soundness of the affine one.
struct Verdicts {
    lipschitz: Vec<bool>,
    affine: Vec<bool>,
    /// The negative control: Hart's test with the constant divided by
    /// [`CONTROL_DIVISOR`], which is unsound by construction.
    control: Vec<bool>,
    rejected_lipschitz: u32,
    rejected_affine: u32,
    rejected_control: u32,
    only_lipschitz: u32,
    only_affine: u32,
    dense_samples: u64,
    violations: u64,
    min_margin: f64,
}

/// Hart's predicate, transcribed from `cell_is_provably_empty`.
///
/// One centre sample against `l · r` with `r` the circumradius, **strict**, so a
/// value exactly on the bound subdivides.
fn lipschitz_rejects<F: Sdf<Scalar = f64>>(
    field: &F,
    lo: [f64; 3],
    h: f64,
    cell: [u32; 3],
    l: f64,
) -> bool {
    let half = h * 0.5;
    let centre = [
        lo[0] + h * f64::from(cell[0]) + half,
        lo[1] + h * f64::from(cell[1]) + half,
        lo[2] + h * f64::from(cell[2]) + half,
    ];
    let radius = half * CIRCUMRADIUS;
    field.sample(centre).abs() > l * radius
}

/// The affine form over one cell's box: the three symbols are the three axes.
fn cell_box(lo: [f64; 3], h: f64, cell: [u32; 3]) -> [Af; 3] {
    let half = h * 0.5;
    [
        Af::symbol(lo[0] + h * f64::from(cell[0]) + half, half, 0),
        Af::symbol(lo[1] + h * f64::from(cell[1]) + half, half, 1),
        Af::symbol(lo[2] + h * f64::from(cell[2]) + half, half, 2),
    ]
}

/// Both rejection tests over the whole grid, plus the dense soundness check.
fn decide<F>(field: &F, lo: [f64; 3], h: f64, l: f64) -> Verdicts
where
    F: Sdf<Scalar = f64> + AffineForm,
{
    let n = CELLS as usize;
    let mut out = Verdicts {
        lipschitz: vec![false; n],
        affine: vec![false; n],
        control: vec![false; n],
        rejected_lipschitz: 0,
        rejected_affine: 0,
        rejected_control: 0,
        only_lipschitz: 0,
        only_affine: 0,
        dense_samples: 0,
        violations: 0,
        min_margin: f64::INFINITY,
    };

    for k in 0..CELLS_PER_AXIS {
        for j in 0..CELLS_PER_AXIS {
            for i in 0..CELLS_PER_AXIS {
                let cell = [i, j, k];
                let idx = cell_index(i, j, k);
                let lip = lipschitz_rejects(field, lo, h, cell, l);
                let form = field.affine(&cell_box(lo, h, cell));
                let aff = form.excludes_zero();
                out.lipschitz[idx] = lip;
                out.affine[idx] = aff;
                out.rejected_lipschitz += u32::from(lip);
                out.rejected_affine += u32::from(aff);
                out.only_lipschitz += u32::from(lip && !aff);
                out.only_affine += u32::from(aff && !lip);
                let ctl = lipschitz_rejects(field, lo, h, cell, l / CONTROL_DIVISOR);
                out.control[idx] = ctl;
                out.rejected_control += u32::from(ctl);

                // The enclosure must contain every value the crate's own `sample`
                // returns inside the closed cell. `t = i / (DENSE − 1)` is exactly
                // 0 and 1 at the ends, so every probe is in the box the form was
                // asked about.
                let (flo, fhi) = (form.lo(), form.hi());
                let corner = [
                    lo[0] + h * f64::from(i),
                    lo[1] + h * f64::from(j),
                    lo[2] + h * f64::from(k),
                ];
                for a in 0..DENSE {
                    for b in 0..DENSE {
                        for c in 0..DENSE {
                            let step = |m: u32| h * f64::from(m) / f64::from(DENSE - 1);
                            let p = [
                                corner[0] + step(a),
                                corner[1] + step(b),
                                corner[2] + step(c),
                            ];
                            let v = field.sample(p);
                            let margin = (v - flo).min(fhi - v);
                            if margin < 0.0 {
                                out.violations += 1;
                            }
                            out.min_margin = out.min_margin.min(margin);
                            out.dense_samples += 1;
                        }
                    }
                }
            }
        }
    }
    out
}

// ─── the mesh ───────────────────────────────────────────────────────────────

/// What one field's meshes came to.
struct Meshes {
    /// The crate's whole-grid mesh with no rejection at all.
    ungated: u64,
    /// The crate's whole-grid mesh with `set_lipschitz`.
    lipschitz_gated: u64,
    triangles: usize,
    skipped_tetrahedra: u64,
    /// Per-cell rig, nothing gated.
    rig_ungated: u64,
    rig_triangles: usize,
    /// Per-cell rig, this field's Lipschitz-rejected cells skipped.
    rig_lipschitz: u64,
    /// Per-cell rig, this field's affine-rejected cells skipped.
    rig_affine: u64,
    /// Per-cell rig, gated by the deliberately unsound control.
    rig_control: u64,
    /// Vertices of the ungated whole-grid mesh lying in a rejected cell's closed
    /// box, per test. Must be zero: a vertex there is a hole.
    witness_lipschitz: u64,
    witness_affine: u64,
    witness_control: u64,
}

/// Every vertex the extractor emits from a cell lies in that cell's closed box,
/// so counting the ungated mesh's vertices against a rejected set decides whether
/// gating could have removed anything.
fn witness(mesh: &MeshBuffer<f64>, lo: [f64; 3], h: f64, rejected: &[bool]) -> u64 {
    let eps = h * 1e-9;
    let last = (CELLS_PER_AXIS - 1) as i64;
    let mut hits = 0u64;
    for p in &mesh.positions {
        let mut span = [[0i64; 2]; 3];
        for (axis, s) in span.iter_mut().enumerate() {
            let t = (p[axis] - lo[axis]) / h;
            let a = libm::floor(t - eps) as i64;
            let b = libm::floor(t + eps) as i64;
            *s = [a.clamp(0, last), b.clamp(0, last)];
        }
        let mut inside = false;
        for k in span[2][0]..=span[2][1] {
            for j in span[1][0]..=span[1][1] {
                for i in span[0][0]..=span[0][1] {
                    if rejected[cell_index(i as u32, j as u32, k as u32)] {
                        inside = true;
                    }
                }
            }
        }
        hits += u64::from(inside);
    }
    hits
}

/// Extract every way, hash every result.
fn meshes<F>(field: &F, lo: [f64; 3], h: f64, l: f64, v: &Verdicts) -> Meshes
where
    F: Sdf<Scalar = f64>,
{
    let shape = RuntimeShape3::new([SAMPLES; 3]).expect("experiment grid fits u32");

    let mut plain = MeshBuffer::<f64>::new();
    let mut ext = SubgridMarchingTetrahedra::<f64>::new(EDGE_SAMPLES).expect("valid resolution");
    ext.extract(field, &shape, lo, h, &mut plain)
        .expect("whole-grid extraction");
    let skipped = ext.report().skipped_tetrahedra;

    let mut gated = MeshBuffer::<f64>::new();
    let mut fast = SubgridMarchingTetrahedra::<f64>::new(EDGE_SAMPLES).expect("valid resolution");
    fast.set_lipschitz(Some(l));
    fast.extract(field, &shape, lo, h, &mut gated)
        .expect("whole-grid gated extraction");

    // The rig. The extractor exposes no per-cell hook, so each cell is extracted
    // on its own 2³ grid and appended; gating is then not appending. Cell
    // geometry is `origin + cell_size` either way, so the triangles are the same
    // ones — only vertex identity at shared edges differs, and that difference is
    // present in both the gated and the ungated rig mesh.
    let unit = RuntimeShape3::new([2; 3]).expect("unit grid fits u32");
    let mut cellwise =
        SubgridMarchingTetrahedra::<f64>::new(EDGE_SAMPLES).expect("valid resolution");
    let mut scratch = MeshBuffer::<f64>::new();
    let mut rig_ungated = MeshBuffer::<f64>::new();
    let mut rig_lip = MeshBuffer::<f64>::new();
    let mut rig_aff = MeshBuffer::<f64>::new();
    let mut rig_ctl = MeshBuffer::<f64>::new();
    for k in 0..CELLS_PER_AXIS {
        for j in 0..CELLS_PER_AXIS {
            for i in 0..CELLS_PER_AXIS {
                let idx = cell_index(i, j, k);
                let corner = [
                    lo[0] + h * f64::from(i),
                    lo[1] + h * f64::from(j),
                    lo[2] + h * f64::from(k),
                ];
                scratch.reset();
                cellwise
                    .extract(field, &unit, corner, h, &mut scratch)
                    .expect("per-cell extraction");
                rig_ungated.append(&scratch).expect("rig mesh fits u32");
                if !v.lipschitz[idx] {
                    rig_lip.append(&scratch).expect("rig mesh fits u32");
                }
                if !v.affine[idx] {
                    rig_aff.append(&scratch).expect("rig mesh fits u32");
                }
                if !v.control[idx] {
                    rig_ctl.append(&scratch).expect("rig mesh fits u32");
                }
            }
        }
    }

    Meshes {
        ungated: mesh_hash(&plain),
        lipschitz_gated: mesh_hash(&gated),
        triangles: plain.triangle_count(),
        skipped_tetrahedra: skipped,
        rig_ungated: mesh_hash(&rig_ungated),
        rig_triangles: rig_ungated.triangle_count(),
        rig_lipschitz: mesh_hash(&rig_lip),
        rig_affine: mesh_hash(&rig_aff),
        rig_control: mesh_hash(&rig_ctl),
        witness_lipschitz: witness(&plain, lo, h, &v.lipschitz),
        witness_affine: witness(&plain, lo, h, &v.affine),
        witness_control: witness(&plain, lo, h, &v.control),
    }
}

// ─── cost, recorded and gating nothing ──────────────────────────────────────

/// Best of three passes, in nanoseconds per cell. Beside the verdict, never in it.
fn ns_per_cell(mut pass: impl FnMut() -> u64) -> f64 {
    let mut best = f64::INFINITY;
    for _ in 0..3 {
        let start = Instant::now();
        let acc = pass();
        let elapsed = start.elapsed().as_nanos() as f64;
        std::hint::black_box(acc);
        best = best.min(elapsed / f64::from(CELLS));
    }
    best
}

// ─── the experiment ─────────────────────────────────────────────────────────

/// One field: decide every cell both ways, extract every mesh, record two rows.
///
/// `sup` is the field's **measured** gradient supremum, which is not the same
/// number as `l` and is not usable as one: a sampled maximum is a lower bound on
/// the true supremum, so declaring it would be unsound in exactly the direction
/// `validate/field_bound.rs` refuses. It is here only to answer the question the
/// registration's own falsification clause asks — whether M-267's 2× gap between
/// `gyroid`'s declared `2√3` and its measured `1.731` is what costs Hart's test,
/// or whether the ball is. `rejected_at_gradient_sup` is therefore a
/// **counterfactual** column, never gated against the mesh and never a bound
/// anyone may use.
fn sweep<F>(
    run: &mut Run,
    name: &'static str,
    field: &F,
    half_domain: f64,
    l: f64,
    sup: f64,
    probe: &Probe,
) where
    F: Sdf<Scalar = f64> + AffineForm,
{
    let lo = [-half_domain; 3];
    let h = 2.0 * half_domain / f64::from(CELLS_PER_AXIS);

    let v = decide(field, lo, h, l);
    let m = meshes(field, lo, h, l, &v);

    let mut at_sup = 0u32;
    for k in 0..CELLS_PER_AXIS {
        for j in 0..CELLS_PER_AXIS {
            for i in 0..CELLS_PER_AXIS {
                at_sup += u32::from(lipschitz_rejects(field, lo, h, [i, j, k], sup));
            }
        }
    }

    let lip_ns = ns_per_cell(|| {
        let mut acc = 0u64;
        for k in 0..CELLS_PER_AXIS {
            for j in 0..CELLS_PER_AXIS {
                for i in 0..CELLS_PER_AXIS {
                    acc += u64::from(lipschitz_rejects(field, lo, h, [i, j, k], l));
                }
            }
        }
        acc
    });
    let aff_ns = ns_per_cell(|| {
        let mut acc = 0u64;
        for k in 0..CELLS_PER_AXIS {
            for j in 0..CELLS_PER_AXIS {
                for i in 0..CELLS_PER_AXIS {
                    acc += u64::from(field.affine(&cell_box(lo, h, [i, j, k])).excludes_zero());
                }
            }
        }
        acc
    });

    let ratio = |rejected: u32| -> f64 {
        if v.rejected_lipschitz == 0 {
            if rejected == 0 { 1.0 } else { f64::INFINITY }
        } else {
            f64::from(rejected) / f64::from(v.rejected_lipschitz)
        }
    };

    let control_identical = m.rig_control == m.rig_ungated && m.witness_control == 0;

    println!(
        "  {name:<16} lipschitz {:>5}  affine {:>5}  ratio {:>7.3}  \
         min margin {:>10.3e}  violations {}  control ({} cells) detected: {}",
        v.rejected_lipschitz,
        v.rejected_affine,
        ratio(v.rejected_affine),
        v.min_margin,
        v.violations,
        v.rejected_control,
        !control_identical
    );

    for (bound, rejected, gated_hash, wit, ns) in [
        (
            "hart_lipschitz",
            v.rejected_lipschitz,
            m.rig_lipschitz,
            m.witness_lipschitz,
            lip_ns,
        ),
        (
            "revised_affine",
            v.rejected_affine,
            m.rig_affine,
            m.witness_affine,
            aff_ns,
        ),
    ] {
        let identical = gated_hash == m.rig_ungated && wit == 0;
        run.record(&[
            ("field", name.to_string()),
            ("samples_per_axis", SAMPLES.to_string()),
            ("bound", bound.to_string()),
            ("cells", CELLS.to_string()),
            ("rejected_cells", rejected.to_string()),
            (
                "rejected_fraction",
                format!("{:.6}", f64::from(rejected) / f64::from(CELLS)),
            ),
            (
                "rejected_ratio_vs_lipschitz",
                format!("{:.6}", ratio(rejected)),
            ),
            ("has_min_max", F::HAS_MIN_MAX.to_string()),
            ("mesh_identical", identical.to_string()),
            ("mesh_hash", format!("{gated_hash:016x}")),
            // One bound evaluation per cell, both tests. The integer the
            // registration asked for; the nanoseconds beside it gate nothing.
            ("bound_evals", CELLS.to_string()),
            ("bound_ns_per_cell", format!("{ns:.1}")),
            // ── beyond the registration ──
            ("lipschitz_constant", format!("{l:.9}")),
            ("rejected_lipschitz", v.rejected_lipschitz.to_string()),
            ("rejected_affine", v.rejected_affine.to_string()),
            ("rejected_only_lipschitz", v.only_lipschitz.to_string()),
            ("rejected_only_affine", v.only_affine.to_string()),
            ("mesh_hash_rig_ungated", format!("{:016x}", m.rig_ungated)),
            ("mesh_hash_rig_gated", format!("{gated_hash:016x}")),
            ("mesh_hash_ungated", format!("{:016x}", m.ungated)),
            (
                "mesh_hash_lipschitz_gated",
                format!("{:016x}", m.lipschitz_gated),
            ),
            ("baseline_vertices_in_rejected_cells", wit.to_string()),
            ("wholegrid_triangles", m.triangles.to_string()),
            ("rig_triangles", m.rig_triangles.to_string()),
            ("skipped_tetrahedra", m.skipped_tetrahedra.to_string()),
            ("edge_samples", EDGE_SAMPLES.to_string()),
            ("dense_per_axis", DENSE.to_string()),
            ("soundness_samples", v.dense_samples.to_string()),
            ("soundness_violations", v.violations.to_string()),
            ("soundness_min_margin", format!("{:.6e}", v.min_margin)),
            ("substitution_samples", probe.samples.to_string()),
            ("substitution_violations", probe.violations.to_string()),
            (
                "substitution_min_margin",
                format!("{:.6e}", probe.min_margin),
            ),
            // The negative control, so `mesh_identical` is a gate and not a
            // tautology. `control_detected` must be `true`.
            ("control_divisor", format!("{CONTROL_DIVISOR:.1}")),
            ("control_rejected_cells", v.rejected_control.to_string()),
            (
                "control_vertices_in_rejected_cells",
                m.witness_control.to_string(),
            ),
            ("control_mesh_hash", format!("{:016x}", m.rig_control)),
            ("control_detected", (!control_identical).to_string()),
            // Counterfactual, and unsound as a bound: see `sweep`'s docs.
            ("gradient_sup", format!("{sup:.9}")),
            ("rejected_at_gradient_sup", at_sup.to_string()),
        ]);
    }
}

fn main() {
    let prereg = isomesh::experiment!("P-54");

    let probe = verify_substitutions();
    println!(
        "substitution soundness: {} samples, {} violations, min margin {:.6e}",
        probe.samples, probe.violations, probe.min_margin
    );
    println!();

    // The measured gradient supremum per field, for the counterfactual column.
    // The five compact fields are exact distances or built from them by min/max,
    // so `‖∇f‖ = 1` and the counterfactual coincides with the declared constant.
    // `gyroid` is the one field where the two differ, and M-267 is where 1.731
    // comes from.
    const EIKONAL_SUP: f64 = 1.0;
    const GYROID_SUP: f64 = 1.731;

    common::experiment::run(prereg, |run| {
        let sphere = Sphere::<f64>::canonical();
        sweep(run, "sphere", &sphere, 2.0, 1.0, EIKONAL_SUP, &probe);

        let torus = Torus::<f64>::canonical();
        sweep(run, "torus", &torus, 2.0, 1.0, EIKONAL_SUP, &probe);

        let boxed = BoxExact::<f64>::canonical();
        sweep(run, "box_exact", &boxed, 2.0, 1.0, EIKONAL_SUP, &probe);

        let csg = csg_difference::<f64>();
        sweep(run, "csg_difference", &csg, 2.0, 1.0, EIKONAL_SUP, &probe);

        // The registered baseline: the crate's reference field named `gyroid`,
        // which is `Intersection<Gyroid, Sphere>` and rejects 688 of 4,096.
        let capped = capped_gyroid::<f64>();
        let capped_l = capped
            .bound()
            .lipschitz()
            .expect("capped gyroid declares a constant");
        sweep(run, "gyroid", &capped, 7.0, capped_l, GYROID_SUP, &probe);

        // The field the mechanism is actually about: the same gyroid on the same
        // domain with the same declared constant, minus the cap.
        let bare = Gyroid::<f64>::canonical();
        let bare_l = bare
            .value_bound()
            .lipschitz()
            .expect("gyroid declares a constant");
        sweep(
            run,
            "gyroid_uncapped",
            &bare,
            7.0,
            bare_l,
            GYROID_SUP,
            &probe,
        );
    });
}
