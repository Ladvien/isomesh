//! E-308 — a bound that reads the *expression* rejects 3.85x more empty cells.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example affine_rejection --release
//! ```
//!
//! **Always `--release`.** A debug build evaluates the affine form 4,096 times
//! per field at startup and the window takes tens of seconds to appear.
//!
//! `1` the capped gyroid, `2` `box_exact`, `3` the uncapped gyroid. `A` draws
//! every cell at once instead of one slab, `Space` freezes the slab sweep. The
//! rest are the shared keys — `W` wireframe, `G` domain box, `F12` screenshot.
//!
//! Under `ISOMESH_CAPTURE` it needs no keyboard: the slab sweeps the domain once
//! per field and the field advances on the captured-frame counter, so a capture
//! of *any* length shows all three. `ISOMESH_FIELD` pins one field even under
//! capture.
//!
//! ```bash
//! ISOMESH_CAPTURE_FRAMES=90 ISOMESH_CAPTURE_EVERY=2 ISOMESH_WINDOW=1280x720 \
//!   FPS=15 ./scripts/record_gif.sh affine_rejection docs/gifs/e308.gif
//! ```
//!
//! Demonstrates **M-354 / P-54** (`docs/experiments/p-54.csv`).
//!
//! # The two rejection tests, side by side
//!
//! Both answer one question — *is this cell provably empty, so extraction may
//! skip it* — and both are sound. They differ in what they are allowed to look
//! at.
//!
//! **Hart's**, which the crate already ships in
//! `subgrid/extract.rs::cell_is_provably_empty`: one sample at the cell centre,
//! rejected when `|f(centre)| > l · r` with `r = (√3/2)·h` the circumradius of
//! the cell and `l` the field's declared Lipschitz constant. It reasons over the
//! **ball** that circumscribes the cell, and it knows nothing about `f` except
//! one value and one constant.
//!
//! **The revised affine form** (Fryazinov, Pasko & Comninos,
//! `10.1016/j.cag.2010.07.003`): `x̂ = x0 + x1·e1 + x2·e2 + x3·e3 + ex·[-1,1]`,
//! five stored `f64` that never grow, because every non-affine error accumulates
//! into the trailing `ex` rather than opening a new symbol. A cell **box** maps
//! in as `centre_k + half_k · e_k`, so the three noise symbols *are* the three
//! axes; the range is `[x0 − rad, x0 + rad]` with `rad = |x1|+|x2|+|x3| + ex`,
//! and the cell is provably empty exactly when that range excludes zero.
//!
//! [`Af`] below is that form, transcribed from
//! `crates/isomesh/benches/experiment_p54.rs`, which is where its soundness is
//! measured: `9³` probes inside every cell of every field — `2,985,984` per
//! field — with **zero** violations on all six, the enclosure never coming
//! closer to failing than `2.54e-13`. Nothing is re-tuned here; this file
//! reproduces the bench's integers and says so.
//!
//! # The mechanism is correlation, not a better constant
//!
//! `sin(x)cos(y) + sin(y)cos(z) + sin(z)cos(x)` cannot have all three terms
//! extremal at once. A per-term interval bound throws that away, and so does a
//! ball of radius `(√3/2)·h` around one sample. The affine form keeps it,
//! because two operands' `e_k` coefficients **subtract**.
//!
//! That is why the win is not uniform, and the field cycle is the whole demo:
//!
//! | field | Hart | affine | ratio | only affine |
//! |---|---:|---:|---:|---:|
//! | `gyroid` (capped) | 688 | **2647** | **3.847384x** | **1959** |
//! | `box_exact` | 3312 | 3312 | 1.000000x | **0** |
//! | `gyroid_uncapped` | **0** | **1098** | inf | **1098** |
//!
//! ## Why `box_exact` gains exactly nothing
//!
//! §4.4 of Fryazinov et al. concludes that general conditionals are an open
//! problem, and the paper gives **no** affine rule for `min`, `max` or `abs`.
//! The only sound treatment is to collapse both operands to their ranges, take
//! the interval `min`/`max`, and hand back a **degenerate** form carrying the
//! whole width in `ex`. Every noise symbol is gone at that point.
//!
//! `box_exact` is `‖max(q, 0)‖ + min(max q, 0)` with `q = |p − c| − b`, and its
//! `abs` and `max(q, 0)` are at the **leaves**. So the form collapses on the
//! first combinator and everything after it is plain interval arithmetic — which
//! on an exact distance field is exactly what Hart's test already was. The two
//! rejected sets are then not merely the same size: they are the same set, cell
//! for cell, `only_affine` **and** `only_lipschitz` both zero.
//!
//! ## And why the capped gyroid gains anyway, despite `has_min_max`
//!
//! `gyroid` is `max(gyroid, sphere(6))`, so it carries a `max` too — and still
//! rejects 3.85x more. The `max` is the **last** operation. The tightness the
//! correlated trig earned is already in the operand's range when `max` reads it,
//! and there is nothing downstream left to cancel, so the collapse costs
//! nothing. Position in the expression is what matters, not presence.
//!
//! ## `gyroid_uncapped` is where Hart's test has nothing to say at all
//!
//! Same field, same domain, same declared constant, minus the cap.
//! `l · r = 2√3 · (√3/2) · 0.875 = 2.625` exceeds the gyroid's entire range of
//! `±1.5`, so Hart's test cannot reject **a single cell** at this resolution —
//! `0` of 4,096 — while the affine form rejects `1098`. The ratio column in the
//! ledger reads `inf`, and that is not a division artefact: it is a bound with
//! a true statement to make where the other one is silent.
//!
//! # This is not "affine always wins"
//!
//! It is not, and the ledger says so on rows this demo does not show: on `torus`
//! the affine form rejects **88 fewer** cells than Hart's (3680 against 3768),
//! and on `csg_difference` **17 fewer**. Both are exact distance fields, where
//! Hart's constant is already tight and the box the affine form reasons over
//! buys nothing back. The magenta colour below is `only_lipschitz` — cells Hart
//! rejects and the affine form does not — and it is drawn whenever it is
//! non-empty precisely so that the picture is capable of showing the loss. On
//! these three fields it is empty, which is the claim that the affine set here
//! **contains** Hart's rather than merely outnumbering it.
//!
//! # What is on screen
//!
//! A 17³ grid — 16 cells per axis, 4,096 cells — one slab at a time, sweeping
//! along `x`, with the extracted surface inside it.
//!
//! - **Amber cages** — `only_affine`. Cells the revised affine form proves empty
//!   and Hart's test cannot. **This is the finding.**
//! - **Dim slate cages** — rejected by both. The bulk that was always cheap.
//! - **Magenta cages** — `only_lipschitz`, rejected by Hart's and not by the
//!   affine form. Empty on all three fields here; see above.
//! - **No cage** — rejected by neither. The surface band, left clear so the mesh
//!   is visible through it.
//!
//! A cage is drawn at the exact bounds of its cell. `A` drops the slab and draws
//! the whole volume, which is the same information with the count made visceral
//! and the geometry buried.
//!
//! # The self-check runs before the window opens
//!
//! `docs/experiments/p-54.csv` is **compiled in** with [`include_str!`], not
//! quoted from memory: [`ledger_row`] parses the committed file and
//! [`SelfCheck`] compares every integer this example computes live against it.
//! Three rows, six integers and a ratio each. A disagreement logs at `error!`
//! and says so on the HUD rather than panicking, because a demo a stranger runs
//! is not the place for an assertion — but it is also not a thing to hide.
//!
//! # `f64`, and the ledger's own extractor
//!
//! M-354 was measured in `f64`, so the numbers on the HUD are reproducible only
//! in `f64`; the surface is cast to `f32` on its way into the [`Mesh`] asset.
//! The mesh is [`SubgridMarchingTetrahedra`] at 6 edge samples, ungated, which
//! is the ledger's `wholegrid_triangles` column — so the triangle count on the
//! HUD is a fourth number checked against the CSV.
//!
//! Every transcendental goes through [`isomesh::Real`], which binds the same
//! `libm` entry points `real.rs` gives `f64` and the bench calls directly.
//! `f64`'s **inherent** `sin` from `std` is a different function, so every call
//! below is spelled out through the trait rather than as a method.

mod common;

use std::ops::{Add, Mul, Neg, Sub};
use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{
    BoundedSdf, BoxExact, Gyroid, Intersection, ReferenceField, Sphere, capped_gyroid,
};
use isomesh::subgrid::extract::SubgridMarchingTetrahedra;
use isomesh::{MeshBuffer, MeshSink, Real, RuntimeShape3, Sdf};

// ─── the registered fixture ─────────────────────────────────────────────────

/// P-54's resolution. 17 samples span 16 cells per axis, and every integer in
/// this file was measured on that grid.
const SAMPLES: u32 = 17;

/// Cells per axis.
const CELLS_PER_AXIS: u32 = SAMPLES - 1;

/// The registration's denominator.
const CELLS: u32 = CELLS_PER_AXIS * CELLS_PER_AXIS * CELLS_PER_AXIS;

/// Samples per tetrahedron edge, matching the bench's whole-grid extraction.
const EDGE_SAMPLES: u32 = 6;

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

const PI: f64 = std::f64::consts::PI;
const FRAC_PI_2: f64 = std::f64::consts::FRAC_PI_2;
const TAU: f64 = std::f64::consts::TAU;

// ─── libm, through the trait the crate binds it with ────────────────────────
//
// `f64` has inherent `sin`, `cos`, `sqrt`, `floor` and `acos` from `std`, and an
// inherent method wins method resolution over a trait one. Those are the
// platform's libm; `isomesh::Real` binds the `libm` **crate**, which is what the
// bench calls and therefore what the ledger's integers were produced by. So each
// is spelled out once here and never called as a method again.

/// `libm::sin`, the way `real.rs` binds it for `f64`.
fn sin(x: f64) -> f64 {
    <f64 as Real>::sin(x)
}

/// `libm::cos`.
fn cos(x: f64) -> f64 {
    <f64 as Real>::cos(x)
}

/// `libm::sqrt`. Correctly rounded by IEEE-754, so this one could not differ.
fn sqrt(x: f64) -> f64 {
    <f64 as Real>::sqrt(x)
}

/// `libm::acos`.
fn acos(x: f64) -> f64 {
    <f64 as Real>::acos(x)
}

/// `libm::floor`. Exact.
fn floor(x: f64) -> f64 {
    <f64 as Real>::floor(x)
}

/// `asin`, via the identity `asin x = π/2 − acos x`.
///
/// `Real` binds no `asin`, and reaching for `std`'s would swap libm
/// implementations mid-expression. The identity costs one subtraction, and it is
/// used only to place `cos`'s critical points — where the residual being bounded
/// is stationary, so a last-bit difference in the location moves the bound by the
/// *square* of it.
fn asin(x: f64) -> f64 {
    FRAC_PI_2 - acos(x)
}

// ─── the revised affine form ────────────────────────────────────────────────

/// `x0 + x1·e1 + x2·e2 + x3·e3 + ex·[-1, 1]`, with `ex ≥ 0`.
///
/// Five `f64`. The size is the invariant: no operation here opens a fourth noise
/// symbol, so a hundred-operation expression costs the same forty bytes as one.
///
/// Transcribed from `crates/isomesh/benches/experiment_p54.rs`, where the
/// soundness of every rule below is measured rather than argued.
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
    /// Writing `x̂ = x0 + L + E`, the square is `x0² + 2x0·L + 2x0·E + (L + E)²`,
    /// and `(L + E)²` lies in `[0, rad²]` **exactly** — it is a square. So the
    /// last term contributes `rad²/2` to the centre and `rad²/2` to `ex`, where
    /// [`Mul`] would contribute `rad²` to `ex` alone and leave the centre at
    /// `x0²`, which is twice as wide *and* mis-centred.
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
    /// Exact and correlation-preserving unless the range straddles zero, in
    /// which case there is nothing to preserve: `|x̂|` is then genuinely not
    /// affine in the symbols and the exact interval `[0, max(−lo, hi)]` is the
    /// most that can be said.
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
    /// Sound, exact on independent operands, and the end of every noise symbol
    /// in the expression. See the module docs: the paper offers nothing else.
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
    /// continuous linear image of a box — so bounding the residual there bounds
    /// it at every `(ε, η)`.
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

    /// `sin` by **Chebyshev** substitution: `α` is the secant slope, the best
    /// `L∞` linear fit, which minimises `δ`.
    ///
    /// Chebyshev rather than min-range deliberately. `δ` lands in `ex`, and `ex`
    /// is the part that **cannot cancel**: two operands' `e_k` terms subtract,
    /// two operands' `ex` terms add. In a sum of products of shared arguments —
    /// which is the field this whole demo is about — the quantity to minimise is
    /// the incompressible one.
    ///
    /// The residual `sin(t) − α·t` is extremal at an endpoint or where
    /// `cos t = α`, which is `±acos(α) + 2kπ`.
    fn sin(self) -> Self {
        let (lo, hi) = (self.lo(), self.hi());
        if hi - lo >= TAU {
            return Self::interval(-1.0, 1.0);
        }
        let alpha = if hi > lo {
            ((sin(hi) - sin(lo)) / (hi - lo)).clamp(-1.0, 1.0)
        } else {
            cos(lo)
        };
        let ac = acos(alpha);
        let k0 = floor(lo / TAU);
        let mut criticals = [0.0; 8];
        for (n, dk) in [-1.0, 0.0, 1.0, 2.0].into_iter().enumerate() {
            let base = (k0 + dk) * TAU;
            criticals[2 * n] = ac + base;
            criticals[2 * n + 1] = -ac + base;
        }
        self.substitute(alpha, sin, &criticals)
    }

    /// `cos` by Chebyshev substitution. Critical points solve `sin t = −α`.
    fn cos(self) -> Self {
        let (lo, hi) = (self.lo(), self.hi());
        if hi - lo >= TAU {
            return Self::interval(-1.0, 1.0);
        }
        let alpha = if hi > lo {
            ((cos(hi) - cos(lo)) / (hi - lo)).clamp(-1.0, 1.0)
        } else {
            -sin(lo)
        };
        let asn = asin((-alpha).clamp(-1.0, 1.0));
        let k0 = floor(lo / TAU);
        let mut criticals = [0.0; 8];
        for (n, dk) in [-1.0, 0.0, 1.0, 2.0].into_iter().enumerate() {
            let base = (k0 + dk) * TAU;
            criticals[2 * n] = asn + base;
            criticals[2 * n + 1] = PI - asn + base;
        }
        self.substitute(alpha, cos, &criticals)
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
        let (sl, sh) = (sqrt(lo), sqrt(hi));
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
        base.substitute(alpha, sqrt, &[mid * mid])
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
    /// bracket is affine and keeps its symbols — **that is where correlation
    /// survives a product**, and it is the whole reason `sin(x)cos(y)` summed
    /// three ways is tighter than three intervals added.
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
/// Example-local, mirroring the bench: `crates/isomesh/src/**` carries no affine
/// arithmetic, and nothing checks this transcription at compile time. What
/// checks it is that the counts below have to come out equal to a committed CSV.
trait AffineForm {
    /// Whether the expression contains a `min`, a `max` or a straddling `abs` —
    /// the operations that collapse the form. The CSV's `has_min_max`.
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

// ─── the two predicates over the grid ───────────────────────────────────────

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

/// Which bounds proved this cell empty. What the cage colour is.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Verdict {
    /// Neither bound rejects it: the surface band. Left clear.
    Neither,
    /// Both bounds reject it.
    Both,
    /// The affine form alone. **The finding.**
    OnlyAffine,
    /// Hart's alone. Empty on all three fields here, drawn if it ever is not.
    OnlyLipschitz,
}

// ─── the committed ledger, compiled in ──────────────────────────────────────

/// `docs/experiments/p-54.csv`, byte for byte, at compile time.
///
/// Compiled in rather than read at runtime: a path relative to the working
/// directory is a thing that can be missing, and a self-check that silently
/// skips when its evidence is absent is not a check.
const LEDGER_CSV: &str = include_str!("../../docs/experiments/p-54.csv");

/// The columns of one `revised_affine` row of P-54.
///
/// One row per field carries both bounds' counts, so the `hart_lipschitz` row is
/// redundant here and only `revised_affine` is read.
struct LedgerRow {
    cells: u32,
    rejected_lipschitz: u32,
    rejected_affine: u32,
    only_affine: u32,
    only_lipschitz: u32,
    /// Kept as text: the `gyroid_uncapped` ratio is `inf`, and comparing the
    /// formatted live value against the committed token handles that without a
    /// special case.
    ratio: String,
    has_min_max: bool,
    triangles: u32,
}

/// Parse one field's `revised_affine` row out of [`LEDGER_CSV`].
///
/// Columns are located by **name** from the header line rather than by position.
/// P-54 writes 40 columns and a column added between runs would silently shift
/// every index; a name cannot.
fn ledger_row(field: &str) -> Option<LedgerRow> {
    let mut lines = LEDGER_CSV.lines().filter(|line| !line.starts_with('#'));
    let header: Vec<&str> = lines.next()?.split(',').collect();
    let col = |name: &str| header.iter().position(|h| *h == name);
    let (c_field, c_bound) = (col("field")?, col("bound")?);
    let c_cells = col("cells")?;
    let c_lip = col("rejected_lipschitz")?;
    let c_aff = col("rejected_affine")?;
    let c_only_aff = col("rejected_only_affine")?;
    let c_only_lip = col("rejected_only_lipschitz")?;
    let c_ratio = col("rejected_ratio_vs_lipschitz")?;
    let c_minmax = col("has_min_max")?;
    let c_tris = col("wholegrid_triangles")?;

    for line in lines {
        let v: Vec<&str> = line.split(',').collect();
        if v.get(c_field) != Some(&field) || v.get(c_bound) != Some(&"revised_affine") {
            continue;
        }
        let num = |i: usize| v.get(i).and_then(|s| s.parse::<u32>().ok());
        return Some(LedgerRow {
            cells: num(c_cells)?,
            rejected_lipschitz: num(c_lip)?,
            rejected_affine: num(c_aff)?,
            only_affine: num(c_only_aff)?,
            only_lipschitz: num(c_only_lip)?,
            ratio: (*v.get(c_ratio)?).to_owned(),
            has_min_max: *v.get(c_minmax)? == "true",
            triangles: num(c_tris)?,
        });
    }
    None
}

// ─── one field, decided and meshed ──────────────────────────────────────────

/// The fields, in the order the digit keys select them.
///
/// `gyroid` first because it is the registered baseline and the headline ratio;
/// `box_exact` second because the bright set vanishing is the control; the
/// uncapped gyroid last because it is where Hart's test rejects nothing at all,
/// which is the frame to end a clip on.
const FIELD_COUNT: usize = 3;

/// CSV `field` values, and the keys into [`ledger_row`].
const FIELD_NAMES: [&str; FIELD_COUNT] = ["gyroid", "box_exact", "gyroid_uncapped"];

/// What each field is, in one HUD line.
///
/// **Every HUD body line in this file is budgeted to 62 characters**, prefix
/// included; the title, which carries the longest field name, is the one line
/// that reaches 68. At the harness's 13px font 62 characters is about 490
/// logical pixels, and the clip is captured at 1280x720 and published at 900
/// wide, so the text has to clear the geometry at 0.7x — a line that runs under
/// the subject is a line nobody reads. These three, the two `why` lines and the
/// self-check rows are the ones long enough to have needed measuring.
const FIELD_FORMULA: [&str; FIELD_COUNT] = [
    "max(gyroid, sphere r6) over [-7,7]^3",
    "|p|-b, |max(q,0)| + min(max q,0), over [-2,2]^3",
    "sin x cos y + sin y cos z + sin z cos x, [-7,7]^3",
];

/// Why the affine form gains what it gains here. Two HUD lines each.
const FIELD_WHY: [[&str; 2]; FIELD_COUNT] = [
    [
        "the max is the LAST op, so its collapse costs",
        "nothing: the operand's range is already tight",
    ],
    [
        "min/max at the LEAVES: the form collapses to an",
        "interval at once -- the same set, cell for cell",
    ],
    [
        "pure trig, shared args: no two terms peak at once,",
        "and l*r = 2.625 exceeds the whole range of +/-1.5",
    ],
];

/// Everything one field produced, decided once at startup.
struct Analysis {
    /// Index into [`FIELD_NAMES`].
    index: usize,
    origin: [f64; 3],
    cell_size: f64,
    half_domain: f32,
    lipschitz: f64,
    has_min_max: bool,
    /// One entry per cell, `x` fastest, matching the extractor's loop order.
    verdict: Vec<Verdict>,
    rejected_lipschitz: u32,
    rejected_affine: u32,
    only_affine: u32,
    only_lipschitz: u32,
    ratio: f64,
    decide_ms: f64,
    extract_ms: f64,
    buffer: MeshBuffer<f64>,
}

impl Analysis {
    fn name(&self) -> &'static str {
        FIELD_NAMES[self.index]
    }

    fn triangles(&self) -> usize {
        self.buffer.triangle_count()
    }

    /// Flat index of a cell, `x` fastest.
    fn cell_index(i: u32, j: u32, k: u32) -> usize {
        ((k * CELLS_PER_AXIS + j) * CELLS_PER_AXIS + i) as usize
    }

    /// The world-space minimum corner of a cell.
    fn corner(&self, i: u32, j: u32, k: u32) -> Vec3 {
        Vec3::new(
            (self.origin[0] + self.cell_size * f64::from(i)) as f32,
            (self.origin[1] + self.cell_size * f64::from(j)) as f32,
            (self.origin[2] + self.cell_size * f64::from(k)) as f32,
        )
    }
}

/// Decide all 4,096 cells both ways and extract the ungated whole-grid mesh.
fn analyse<F>(index: usize, field: &F, half_domain: f64, l: f64) -> Analysis
where
    F: Sdf<Scalar = f64> + AffineForm,
{
    let lo = [-half_domain; 3];
    let h = 2.0 * half_domain / f64::from(CELLS_PER_AXIS);

    let started = Instant::now();
    let mut verdict = vec![Verdict::Neither; CELLS as usize];
    let (mut rejected_lipschitz, mut rejected_affine) = (0u32, 0u32);
    let (mut only_affine, mut only_lipschitz) = (0u32, 0u32);
    for k in 0..CELLS_PER_AXIS {
        for j in 0..CELLS_PER_AXIS {
            for i in 0..CELLS_PER_AXIS {
                let cell = [i, j, k];
                let lip = lipschitz_rejects(field, lo, h, cell, l);
                let aff = field.affine(&cell_box(lo, h, cell)).excludes_zero();
                rejected_lipschitz += u32::from(lip);
                rejected_affine += u32::from(aff);
                only_affine += u32::from(aff && !lip);
                only_lipschitz += u32::from(lip && !aff);
                verdict[Analysis::cell_index(i, j, k)] = match (lip, aff) {
                    (true, true) => Verdict::Both,
                    (false, true) => Verdict::OnlyAffine,
                    (true, false) => Verdict::OnlyLipschitz,
                    (false, false) => Verdict::Neither,
                };
            }
        }
    }
    let decide_ms = started.elapsed().as_secs_f64() * 1e3;

    // The bench's `ratio` exactly, `inf` case included: a bound that rejects
    // where the other rejects nothing is not a division by zero to be papered
    // over, it is the reading.
    let ratio = if rejected_lipschitz == 0 {
        if rejected_affine == 0 {
            1.0
        } else {
            f64::INFINITY
        }
    } else {
        f64::from(rejected_affine) / f64::from(rejected_lipschitz)
    };

    let mut buffer = MeshBuffer::<f64>::new();
    let mut extract_ms = 0.0;
    match (
        RuntimeShape3::new([SAMPLES; 3]),
        SubgridMarchingTetrahedra::<f64>::new(EDGE_SAMPLES),
    ) {
        (Ok(shape), Ok(mut ext)) => {
            let started = Instant::now();
            if let Err(error) = ext.extract(field, &shape, lo, h, &mut buffer) {
                error!("{} extraction failed: {error}", FIELD_NAMES[index]);
            }
            extract_ms = started.elapsed().as_secs_f64() * 1e3;
        }
        (shape, ext) => error!(
            "{} rig rejected: shape {:?}, extractor {:?}",
            FIELD_NAMES[index],
            shape.err(),
            ext.err()
        ),
    }

    Analysis {
        index,
        origin: lo,
        cell_size: h,
        half_domain: half_domain as f32,
        lipschitz: l,
        has_min_max: F::HAS_MIN_MAX,
        verdict,
        rejected_lipschitz,
        rejected_affine,
        only_affine,
        only_lipschitz,
        ratio,
        decide_ms,
        extract_ms,
        buffer,
    }
}

/// Build one field by index. Three concrete types, so the index is matched here.
fn field_analysis(index: usize) -> Option<Analysis> {
    match index {
        0 => {
            let field = capped_gyroid::<f64>();
            let l = field.bound().lipschitz()?;
            Some(analyse(index, &field, 7.0, l))
        }
        1 => {
            let field = BoxExact::<f64>::canonical();
            let l = field.bound().lipschitz()?;
            Some(analyse(index, &field, 2.0, l))
        }
        2 => {
            // Not a `ReferenceField`: the bare gyroid is not closed in any
            // domain. Same declared constant, same `[-7,7]³`, minus the cap.
            let field = Gyroid::<f64>::canonical();
            let l = field.value_bound().lipschitz()?;
            Some(analyse(index, &field, 7.0, l))
        }
        _ => None,
    }
}

// ─── the startup self-check ─────────────────────────────────────────────────

/// Every live integer held against the committed CSV, before the window opens.
///
/// A disagreement is **loud and must not take the window down with it**: it logs
/// at `error!` and the HUD says which row failed. That is the house rule for a
/// demo a stranger runs, and it is also the only outcome worth reporting — a
/// live run that disagrees with the ledger is a finding, not something to tune
/// away.
#[derive(Resource)]
struct SelfCheck {
    /// One HUD line per field: the live counts and whether the committed row
    /// agrees. A table rather than a sentence, because three rows compared
    /// column by column is what the check actually is.
    table: Vec<String>,
    /// Whether all three rows matched every column, to the digit.
    all_match: bool,
}

impl SelfCheck {
    fn run(fields: &[Analysis]) -> Self {
        info!(
            "E-308 self-check: P-54 at {SAMPLES}^3 = {CELLS} cells, against the \
             compiled-in docs/experiments/p-54.csv"
        );

        let mut matched = [false; FIELD_COUNT];
        let mut table = Vec::with_capacity(FIELD_COUNT);
        let mut headline = Vec::with_capacity(FIELD_COUNT);
        for (slot, a) in matched.iter_mut().zip(fields) {
            let name = a.name();
            let live_ratio = format!("{:.6}", a.ratio);
            headline.push(format!(
                "{name} {}->{} {live_ratio}x",
                a.rejected_lipschitz, a.rejected_affine
            ));
            let Some(row) = ledger_row(name) else {
                error!(
                    "p-54.csv has no revised_affine row for {name}, so this example cannot \
                     check itself against the ledger it claims to reproduce"
                );
                table.push(format!("{name:<16} NO LEDGER ROW"));
                continue;
            };
            let ok = row.cells == CELLS
                && row.rejected_lipschitz == a.rejected_lipschitz
                && row.rejected_affine == a.rejected_affine
                && row.only_affine == a.only_affine
                && row.only_lipschitz == a.only_lipschitz
                && row.has_min_max == a.has_min_max
                && row.ratio == live_ratio
                && row.triangles as usize == a.triangles();
            *slot = ok;
            // The live numbers, with the ledger's verdict beside them. Sized to
            // 59 characters so the widest HUD line clears the geometry at 900
            // wide, which is the size these clips are published at.
            table.push(format!(
                "{name:<15}{:>5} -> {:>4}  {:>9}  {}",
                a.rejected_lipschitz,
                a.rejected_affine,
                format!("{live_ratio}x"),
                if ok { "MATCHES" } else { "DIFFERS" },
            ));
            if ok {
                info!(
                    "  {name:<16} lipschitz {:>4} -> affine {:>4}  ratio {live_ratio:>8}x  \
                     only_affine {:>4}  only_lipschitz {:>3}  has_min_max {:<5}  \
                     triangles {:>5}  MATCHES p-54.csv",
                    a.rejected_lipschitz,
                    a.rejected_affine,
                    a.only_affine,
                    a.only_lipschitz,
                    a.has_min_max,
                    a.triangles(),
                );
            } else {
                error!(
                    "  {name:<16} DISAGREES WITH p-54.csv. live: cells {CELLS} lipschitz {} \
                     affine {} ratio {live_ratio} only_affine {} only_lipschitz {} \
                     has_min_max {} triangles {}. ledger: cells {} lipschitz {} affine {} \
                     ratio {} only_affine {} only_lipschitz {} has_min_max {} triangles {}. \
                     That is a finding about this build, not a number to adjust.",
                    a.rejected_lipschitz,
                    a.rejected_affine,
                    a.only_affine,
                    a.only_lipschitz,
                    a.has_min_max,
                    a.triangles(),
                    row.cells,
                    row.rejected_lipschitz,
                    row.rejected_affine,
                    row.ratio,
                    row.only_affine,
                    row.only_lipschitz,
                    row.has_min_max,
                    row.triangles,
                );
            }
        }

        let all_match = matched.iter().all(|ok| *ok) && table.len() == FIELD_COUNT;
        if all_match {
            info!(
                "E-308 self-check: 3 of 3 rows reproduce p-54.csv to the digit -- {}",
                headline.join(" | ")
            );
        } else {
            error!(
                "E-308 self-check: {} of {FIELD_COUNT} rows reproduce p-54.csv -- {}",
                matched.iter().filter(|ok| **ok).count(),
                headline.join(" | ")
            );
        }

        Self { table, all_match }
    }
}

// ─── state ──────────────────────────────────────────────────────────────────

/// Every field, decided and meshed once at startup.
///
/// All three cost about a fifth of a second together, so nothing is deferred: a
/// field switch is a handle swap, and a capture never photographs a rebuild.
#[derive(Resource)]
struct Study(Vec<Analysis>);

/// What is being drawn this frame.
#[derive(Resource)]
struct Shown {
    field: usize,
    /// Which `x` slab of cells is cageed, `0..CELLS_PER_AXIS`.
    slab: u32,
    /// Draw the whole volume instead of one slab.
    all: bool,
}

/// A field pinned by `ISOMESH_FIELD`, which overrides the capture's stepping.
#[derive(Resource)]
struct Pinned(Option<usize>);

/// One `Handle<Mesh>` per field, uploaded once.
#[derive(Resource)]
struct Surfaces(Vec<Handle<Mesh>>);

/// Cages get their own group so they can be drawn in front of the surface
/// without dragging the shared wireframe along with them.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct CellGizmos;

// ─── framing ────────────────────────────────────────────────────────────────

/// Orbit radius, in domain half-extents.
///
/// The three domains differ by 3.5x — `box_exact` is half-extent 2 and both
/// gyroids are 7 — so a fixed world radius would put the camera inside one of
/// them.
///
/// `3.5` is the number that just fits the domain vertically at the harness's
/// 45° field of view: the tallest silhouette a half-extent-`H` cube presents at
/// this yaw and pitch is `1.413·H`, and `tan(22.5°)·3.5·H = 1.449·H`. At `3.1`
/// the box was measurably clipped at the bottom of a 1280x720 capture, which on
/// the uncapped gyroid — the one field that fills its whole domain — reads as a
/// surface that carries on past the frame rather than a grid that was measured.
const RADIUS_EXTENTS: f32 = 3.5;

/// How far right of centre the subject sits, as a fraction of the orbit radius.
///
/// **The HUD fills the upper left and the slab is the subject.**
/// Centring it photographs the argument with its evidence hidden. Applied in
/// the camera's own basis rather than as a world offset, so it holds while
/// `ISOMESH_SPIN` yaws.
///
/// Right only, and no downward nudge: the HUD is entirely to the *left* of the
/// subject at this radius, so vertical room spent clearing it would only crop
/// the domain, which is the thing [`RADIUS_EXTENTS`] was widened to stop.
const SUBJECT_SHIFT: f32 = 0.32;

/// Seconds for the slab to cross the domain once, when nobody is capturing.
const SWEEP_SECONDS: f32 = 7.0;

// ─── colours ────────────────────────────────────────────────────────────────

/// Rejected by both. Subordinate to the amber, because it is the part that was
/// always cheap — but not *invisible*: measured on a 1280x720 capture, at
/// `0.30, 0.34, 0.43` a slate cage was indistinguishable from both the dark
/// background and the harness's own domain box, so "both bounds agree here"
/// read as "nothing here", which is the wrong claim.
const BOTH: Color = Color::srgb(0.38, 0.52, 0.70);

/// Rejected by the affine form alone. **The finding.**
const ONLY_AFFINE: Color = Color::srgb(1.0, 0.70, 0.10);

/// Rejected by Hart's alone. Empty on these three fields; drawn if it is not.
const ONLY_LIPSCHITZ: Color = Color::srgb(1.0, 0.24, 0.85);

// ─── app ────────────────────────────────────────────────────────────────────

fn main() {
    let mut app = App::new();

    // **`LogPlugin` alone first, then the self-check, then everything else.**
    //
    // Not a flourish: `WinitPlugin` builds the event loop inside `build`, so
    // `add_plugins(DefaultPlugins)` *panics* on a machine with no display —
    // before a single line of this example runs. A self-check that only happens
    // once a window opens cannot be run on the machine where somebody doubts
    // the numbers, and a demo whose claim is "these integers are the committed
    // ones" should be able to say so from a terminal.
    //
    // Adding the subscriber here and disabling the copy inside `DefaultPlugins`
    // is what makes that possible: `LogPlugin` installs a *global* subscriber
    // and would panic if it were installed twice.
    app.add_plugins(bevy::log::LogPlugin::default());
    let study = build_study();
    let check = SelfCheck::run(&study.0);

    app.add_plugins(
        DefaultPlugins
            .set(WindowPlugin {
                primary_window: Some(Window {
                    title: "isomesh - E-308 affine rejection".into(),
                    ..default()
                }),
                ..default()
            })
            .disable::<bevy::log::LogPlugin>(),
    )
    .add_plugins(CommonPlugin)
    .init_gizmo_group::<CellGizmos>();

    app.insert_resource(study)
        .insert_resource(check)
        .insert_resource(Pinned(pinned_field()))
        .insert_resource(Shown {
            field: pinned_field().unwrap_or(0),
            slab: 0,
            all: false,
        })
        .add_systems(Startup, setup)
        // `PreUpdate` for E-306's reason: the harness's `update_hud` renders
        // `DemoStats` and its `capture_sequence` advances `Capture::taken`, both
        // in `Update` with no ordering against an example's own systems. In
        // `Update` the HUD would render a frame-old readout beside a current
        // slab, which for a demo whose whole claim is "these numbers describe
        // this picture" is the one defect that matters.
        .add_systems(
            PreUpdate,
            (advance, apply, frame_camera, draw_cells, report)
                .chain()
                .after(bevy::input::InputSystems),
        )
        .run();
}

/// Decide and mesh all three fields, then hold them against the ledger.
fn build_study() -> Study {
    let fields: Vec<Analysis> = (0..FIELD_COUNT).filter_map(field_analysis).collect();
    if fields.len() != FIELD_COUNT {
        error!(
            "only {} of {FIELD_COUNT} fields could be built; the missing ones declare no \
             Lipschitz constant, so Hart's test has nothing to run",
            fields.len()
        );
    }
    Study(fields)
}

/// The field `ISOMESH_FIELD` asks for, if it asks for one.
fn pinned_field() -> Option<usize> {
    let raw = std::env::var("ISOMESH_FIELD").ok()?;
    match raw.trim().parse::<usize>() {
        Ok(index) if index < FIELD_COUNT => Some(index),
        _ => {
            error!("ISOMESH_FIELD={raw} is not one of 0..{FIELD_COUNT}");
            None
        }
    }
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut config: ResMut<GizmoConfigStore>,
    mut camera: Query<&mut OrbitCamera>,
    study: Res<Study>,
) {
    let (cells, _) = config.config_mut::<CellGizmos>();
    cells.line.width = 1.7;
    // Negative, so a cage is never lost behind the surface that passes through
    // it. A cell's verdict is a statement about the whole box, including the
    // half of it facing away from the camera.
    cells.depth_bias = -0.5;

    // Darker than the harness's usual surface grey. The amber cages are the
    // subject and the light HUD text sits over the upper left; at the default
    // 0.72 both wash out against a gyroid that fills the frame.
    let surface = materials.add(StandardMaterial {
        base_color: Color::srgb(0.30, 0.33, 0.39),
        perceptual_roughness: 0.62,
        metallic: 0.04,
        double_sided: true,
        cull_mode: None,
        ..default()
    });

    let handles: Vec<Handle<Mesh>> = study.0.iter().map(|a| meshes.add(to_mesh(a))).collect();

    // `Mesh3d::default()` names no asset, so nothing is uploaded until `apply`
    // picks a field. An empty mesh would be worse: `MeshAllocator` skips a
    // zero-byte vertex buffer and then copies into it anyway, once per frame.
    commands.spawn((Mesh3d::default(), MeshMaterial3d(surface), DemoMesh));

    // Spawned rather than assumed: `draw_domain` queries for it, and without one
    // the `G` toggle silently does nothing. Filled in by `apply`.
    commands.spawn(DemoDomain {
        min: Vec3::splat(-1.0),
        max: Vec3::splat(1.0),
    });

    for mut orbit in &mut camera {
        orbit.yaw = 0.72;
        orbit.pitch = 0.34;
    }

    commands.insert_resource(Surfaces(handles));
}

/// The `f64` extraction as a Bevy mesh.
///
/// Cast rather than re-extracted in `f32`: the triangle count on the HUD is
/// checked against the ledger, and the ledger's is an `f64` number.
fn to_mesh(analysis: &Analysis) -> Mesh {
    let buffer = &analysis.buffer;
    let mut builder = MeshBuilder::new();
    for i in 0..buffer.positions.len() {
        let (Some(p), Some(n)) = (buffer.positions.get(i), buffer.normals.get(i)) else {
            continue;
        };
        builder.vertex(
            [p[0] as f32, p[1] as f32, p[2] as f32],
            [n[0] as f32, n[1] as f32, n[2] as f32],
        );
    }
    for t in buffer.indices.chunks_exact(3) {
        let (Some(a), Some(b), Some(c)) = (t.first(), t.get(1), t.get(2)) else {
            continue;
        };
        builder.triangle(*a, *b, *c);
    }
    builder.into_mesh()
}

/// Frames a capture runs for.
///
/// Read from the environment rather than from [`Capture`], which keeps its
/// length private, because pacing the sweep off the capture is what stops a
/// six-frame smoke test and a ninety-frame clip from both being a still.
fn capture_frames() -> u32 {
    std::env::var("ISOMESH_CAPTURE_FRAMES")
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .unwrap_or(60)
        .max(1)
}

/// Choose the field and the slab for this frame.
///
/// Under capture both come off the captured-frame counter, so a clip of any
/// length shows all three fields and one full sweep of each. Interactively the
/// digits pick the field and the slab sweeps on a loop, so the picture moves
/// without anybody pressing anything.
fn advance(
    time: Res<Time>,
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    pinned: Res<Pinned>,
    mut shown: ResMut<Shown>,
    mut clock: Local<f32>,
) {
    if keys.just_pressed(KeyCode::KeyA) {
        shown.all = !shown.all;
    }
    for (n, key) in [KeyCode::Digit1, KeyCode::Digit2, KeyCode::Digit3]
        .into_iter()
        .enumerate()
    {
        if keys.just_pressed(key) {
            shown.field = n;
        }
    }

    if capture.is_active() {
        // A pinned field owns the whole clip, so its sweep is the whole clip.
        // Dividing by `FIELD_COUNT` anyway would spend two thirds of a pinned
        // capture re-photographing the same five slabs.
        let per_field = match pinned.0 {
            Some(_) => capture_frames(),
            // **Rounded up, not down.** `80 / 3` is 26, which leaves frames 78
            // and 79 past the end of the third field's share: measured, and the
            // clip's last two frames restarted the uncapped gyroid's sweep at
            // slab 1 after it had reached slab 16. `div_ceil` spends the
            // remainder on the last field instead, which is one frame short of
            // a full sweep rather than a visible rewind.
            None => (capture_frames().div_ceil(FIELD_COUNT as u32)).max(1),
        };
        let step = capture.taken / per_field.max(1);
        shown.field = pinned
            .0
            .unwrap_or_else(|| (step as usize).min(FIELD_COUNT - 1));
        // The slab sweeps once within the field's share of the clip, and the
        // phase is taken from the captured-frame counter rather than the clock
        // so the same command produces the same frames.
        let phase = f32::from(u16::try_from(capture.taken % per_field.max(1)).unwrap_or(0))
            / f32::from(u16::try_from(per_field).unwrap_or(1).max(1));
        shown.slab = ((phase * CELLS_PER_AXIS as f32) as u32).min(CELLS_PER_AXIS - 1);
        return;
    }

    if let Some(index) = pinned.0 {
        shown.field = index;
    }
    if !flags.paused {
        *clock = (*clock + time.delta_secs() / SWEEP_SECONDS).fract();
    }
    shown.slab = ((*clock * CELLS_PER_AXIS as f32) as u32).min(CELLS_PER_AXIS - 1);
}

/// Swap the surface and the domain box when the field changes.
fn apply(
    shown: Res<Shown>,
    study: Res<Study>,
    surfaces: Res<Surfaces>,
    mut last: Local<Option<usize>>,
    mut mesh: Query<&mut Mesh3d, With<DemoMesh>>,
    mut domain: Query<&mut DemoDomain>,
) {
    if *last == Some(shown.field) {
        return;
    }
    let (Some(analysis), Some(handle)) = (
        study.0.get(shown.field),
        surfaces.0.get(shown.field).cloned(),
    ) else {
        return;
    };
    *last = Some(shown.field);

    for mut slot in &mut mesh {
        slot.0 = handle.clone();
    }
    let half = analysis.half_domain;
    for mut d in &mut domain {
        d.min = Vec3::splat(-half);
        d.max = Vec3::splat(half);
    }
}

/// Keep the domain filling the frame, offset clear of the HUD.
fn frame_camera(shown: Res<Shown>, study: Res<Study>, mut camera: Query<&mut OrbitCamera>) {
    let Some(analysis) = study.0.get(shown.field) else {
        return;
    };
    let radius = analysis.half_domain * RADIUS_EXTENTS;
    for mut orbit in &mut camera {
        // The camera's own basis, from the same yaw/pitch the harness's
        // `orbit_camera` builds its transform from, so the offset is one
        // screen-space nudge however far `ISOMESH_SPIN` has turned.
        let dir = Vec3::new(
            orbit.yaw.cos() * orbit.pitch.cos(),
            orbit.pitch.sin(),
            orbit.yaw.sin() * orbit.pitch.cos(),
        );
        // `orbit_camera` places the eye at `focus + dir * radius`, so the view
        // direction is `-dir` and a focus moved along `-right` puts the subject
        // right of centre.
        let forward = -dir;
        let right = forward.cross(Vec3::Y).normalize_or_zero();
        orbit.focus = -right * (SUBJECT_SHIFT * radius);
        orbit.radius = radius;
    }
}

/// Cage every rejected cell in view, coloured by which bound rejected it.
///
/// One `x` slab at a time by default. The whole volume is the same information —
/// `A` shows it — but 2,647 cages on `gyroid` is 31,764 line segments in front of
/// the surface they are about, and a picture in which nothing can be seen is not
/// a stronger claim for being denser.
fn draw_cells(shown: Res<Shown>, study: Res<Study>, mut gizmos: Gizmos<CellGizmos>) {
    let Some(analysis) = study.0.get(shown.field) else {
        return;
    };
    let size = analysis.cell_size as f32;
    let slab = shown.slab.min(CELLS_PER_AXIS - 1);

    for k in 0..CELLS_PER_AXIS {
        for j in 0..CELLS_PER_AXIS {
            for i in 0..CELLS_PER_AXIS {
                if !shown.all && i != slab {
                    continue;
                }
                let colour = match analysis.verdict[Analysis::cell_index(i, j, k)] {
                    Verdict::Neither => continue,
                    Verdict::Both => BOTH,
                    Verdict::OnlyAffine => ONLY_AFFINE,
                    Verdict::OnlyLipschitz => ONLY_LIPSCHITZ,
                };
                cage(&mut gizmos, analysis.corner(i, j, k), size, colour);
            }
        }
    }
}

/// The twelve edges of one cell, at its exact bounds.
///
/// Exact rather than inflated: a cage larger than its cell would make the
/// rejected set look like it covers the surface band, which is the one thing
/// this picture must not say. Corner indexing matches the extractor's — bit `i`
/// of the corner index is axis `i`.
fn cage(gizmos: &mut Gizmos<CellGizmos>, min: Vec3, size: f32, colour: Color) {
    let corner = |i: usize| {
        min + Vec3::new(
            if i & 1 == 0 { 0.0 } else { size },
            if i & 2 == 0 { 0.0 } else { size },
            if i & 4 == 0 { 0.0 } else { size },
        )
    };
    for i in 0..8usize {
        for axis in 0..3usize {
            let bit = 1 << axis;
            if i & bit == 0 {
                gizmos.line(corner(i), corner(i | bit), colour);
            }
        }
    }
}

/// The HUD. The numbers are the demo.
fn report(
    shown: Res<Shown>,
    study: Res<Study>,
    check: Res<SelfCheck>,
    mut stats: ResMut<DemoStats>,
) {
    let Some(a) = study.0.get(shown.field) else {
        return;
    };
    let index = a.index;

    // The title reaches 68 characters and every line below it stops at 62 --
    // see [`FIELD_FORMULA`] for why those are the budgets.
    stats.title = format!(
        "E-308  affine rejection - {}  {SAMPLES}^3  [1-3] field [A] all",
        a.name(),
    );
    stats.vertices = a.buffer.positions.len();
    stats.triangles = a.triangles();
    stats.extract_ms = a.extract_ms;

    let mut extra = vec![
        format!("field     {}", FIELD_FORMULA[index]),
        format!(
            "          l {:.6}  h {:.6}  min/max {}",
            a.lipschitz, a.cell_size, a.has_min_max,
        ),
        String::new(),
        format!("cells                {CELLS}"),
        format!(
            "rejected_lipschitz  {:>4}   Hart ball, |f(centre)| > l*r",
            a.rejected_lipschitz,
        ),
        format!(
            "rejected_affine     {:>4}   affine range over the cell box",
            a.rejected_affine,
        ),
        format!("ratio          {:>10}", format!("{:.6}x", a.ratio)),
        format!(
            "only_affine         {:>4}   AMBER cages -- the finding",
            a.only_affine,
        ),
        format!(
            "only_lipschitz      {:>4}   {}",
            a.only_lipschitz,
            if a.only_lipschitz == 0 {
                "none: affine contains Hart's"
            } else {
                "MAGENTA -- affine loses these"
            },
        ),
        String::new(),
        format!("why       {}", FIELD_WHY[index][0]),
        format!("          {}", FIELD_WHY[index][1]),
        String::new(),
    ];

    // The self-check table, with the row for the field on screen marked. The
    // three rows are the three the ticket names, and they are the live counts
    // -- the verdict column is what holds them against the committed CSV.
    for (n, row) in check.table.iter().enumerate() {
        let lead = if n == 0 { "selfcheck" } else { "         " };
        let here = if n == index { ">" } else { " " };
        extra.push(format!("{lead} {here}{row}"));
    }
    extra.push(format!(
        "          {} vs docs/experiments/p-54.csv",
        if check.all_match {
            "3/3 rows MATCH"
        } else {
            "MISMATCH -- SEE THE LOG"
        },
    ));
    extra.push(format!(
        "view      {}  decide {:.1} ms / {CELLS} cells",
        if shown.all {
            format!("all {CELLS_PER_AXIS} slabs")
        } else {
            format!("x slab {:>2}/{CELLS_PER_AXIS}", shown.slab + 1)
        },
        a.decide_ms,
    ));

    stats.extra = extra;
}
