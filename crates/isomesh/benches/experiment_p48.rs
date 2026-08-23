//! **P-48 — a sound inclusion function, so the certificate is against the field.**
//!
//! Ticket: R-044. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p48
//! ```
//!
//! Writes `docs/experiments/p-48.csv`.
//!
//! # The gap being closed, in `isotopy.rs`'s own words
//!
//! `validate/isotopy.rs` states it: *"The general form needs interval arithmetic
//! over an arbitrary `F`, which this crate has no way to do — an `Sdf` hands back
//! point values. A sampled hull of the gradient would be a lower bound on its
//! variation, so the predicate could pass where the truth fails, which is the one
//! direction a certificate must never err in."* It therefore certifies the
//! **trilinear interpolant**, where `0 ∉ □F(C)` is exactly "all eight corners
//! share a sign", and says plainly that it does not certify the analytic field.
//!
//! This harness supplies the missing half: a compositional inclusion function
//! over the crate's own field types. It lives here, not in the crate, and it adds
//! no dependency — `libm` is already `isomesh`'s single dependency, so a bench
//! target has it in scope, and every transcendental below is the *same* `libm`
//! call the field itself makes.
//!
//! **This is not a speedup.** It is strictly more work than sampling: it is a
//! proof obligation discharged per cell. The only thing it buys is that a cell it
//! calls surface-free has no surface in it, as opposed to no surface at the eight
//! places anyone looked.
//!
//! # The arithmetic: `(lo, hi)` and one ULP per operation
//!
//! `core` has no directed rounding, so an interval computed with
//! round-to-nearest can be up to half an ULP too narrow at each end — and too
//! narrow is precisely the unsound direction. Every rounding operation therefore
//! ends with [`Iv::widen`], which is one `next_down` on the lower end and one
//! `next_up` on the upper. That is at least twice the error it needs to absorb.
//!
//! Six operations round and are widened: **add, subtract, multiply, square root,
//! sine and cosine**. Four do not round at all and are exact: **negate**,
//! **absolute value**, and interval **min**/**max** — the last two matter,
//! because `min` and `max` of intervals are exactly `min` and `max` of the
//! endpoints, which is what makes the CSG combinators lossless.
//!
//! ## A second, separate slack for the *field's* own rounding
//!
//! One-ULP widening makes the interval a sound enclosure of the **exact real**
//! function. The thing being certified is not the exact real function: it is the
//! `f64` code in `fields/mod.rs`, which has its own rounding error. So each
//! primitive enclosure ends with [`Iv::impl_slack`], widening by
//! `256 · f64::EPSILON · max(1, |lo|, |hi|)`. Two hundred and fifty-six is
//! generous by construction: the longest primitive here (`perlin`) performs well
//! under a hundred rounded operations on quantities bounded by two, so its
//! accumulated error cannot reach a tenth of that. In absolute terms the slack is
//! about `6e-14` against a per-cell field variation of order `0.1`, so it costs
//! nothing measurable in reach and removes the last way a certificate could be
//! wrong.
//!
//! # Each enclosure, and why it is sound
//!
//! - **`Sphere`** — `√(Σ dᵢ²) − r`. `Iv::sqr` is the exact interval square
//!   (`[0, max(a², b²)]` when the interval straddles zero, which `mul` would get
//!   merely sound and loose), so the radicand's lower end is never negative.
//! - **`BoxExact`, `ThinPlate`** — `box_sample` verbatim:
//!   `|d| − b`, then `‖max(q, 0)‖ + min(max q, 0)`. `abs`, `min` and `max` are
//!   all exact on intervals, so the only rounding is the norm.
//! - **`Torus`** — `s = √(dx² + dz²)`, then `√((s − R)² + dy²) − a`.
//! - **`Union`, `Intersection`, `Difference`** — `min(A, B)`, `max(A, B)` and
//!   `max(A, −B)`. Exact, and the reason a composite is no looser than its worst
//!   operand.
//! - **`Gyroid`** — `sin a cos b + sin b cos c + sin c cos a − iso`, over
//!   intervals. The trig is a real monotonic-branch analysis and not an
//!   approximation: `sin` is monotone between consecutive critical points, so on
//!   an interval narrower than `2π` the range is given by the two endpoints
//!   **unless** a critical point `π/2 + 2kπ` (maximum `+1`) or `−π/2 + 2kπ`
//!   (minimum `−1`) lies inside, in which case that end is pinned to `±1`. The
//!   containment test is decided with a guard so that it errs toward *including*
//!   the extremum, which widens and cannot narrow. `cos` is the same with
//!   critical points at `2kπ` and `π + 2kπ`. Endpoints are widened three times
//!   rather than once, because `libm`'s `sin`/`cos` are not correctly rounded:
//!   `libm(x) ≤ true(x) + 1ulp ≤ true(hi) + 1ulp ≤ libm(hi) + 2ulp` on a
//!   monotone branch, so two would do and three is the margin.
//! - **`NoiseVolume`, `FbmTerrain`** — the interesting case, and it is derived
//!   rather than bounded away. `perlin` is
//!   `Σ_{c ∈ {0,1}³} w_c(t) · (g_c · d_c)` with `w_c ≥ 0`, `Σ w_c = 1`,
//!   `d_c = t − c` and `t` the fractional part. Every factor is interval-
//!   evaluable: the quintic fade `u = t³(6t² − 15t + 10)` has derivative
//!   `30t²(t − 1)² ≥ 0` **everywhere**, so it is monotone and its endpoints give
//!   its range exactly; each `g_c · d_c` is `±dᵢ ± dⱼ` because every `GRAD12`
//!   entry has exactly two `±1` components and one zero, so the dot product is
//!   exact sign flips and two additions. A box spanning several lattice cells is
//!   handled by enumerating the cells it overlaps and taking the hull, which is
//!   sound because the cells cover the box. The `t` interval is clamped to
//!   `[0, 1]`, which is exact rather than a fallback: the true fractional part is
//!   always in `[0, 1)`.
//!   `fbm` is `Σᵢ gainⁱ · perlin(lacⁱ · freq · p + offsetᵢ)`, so it is the sum of
//!   per-octave enclosures. `FbmTerrain` passes `y = 0` into the noise, so the
//!   noise is two-dimensional and the `y` interval enters only through the
//!   outer `p[1] − (base + amplitude · n)`.
//!
//! **Term-wise rather than by the convex-combination bound**, deliberately. The
//! bound `min_c inf(g_c·d_c) ≤ value ≤ max_c sup(g_c·d_c)` follows from the
//! weights being a partition of unity and is a one-liner, but it is badly loose
//! near a lattice corner: it admits the far corner's `g·d ≈ −2` while the actual
//! weight on that corner is under a thousandth. Evaluating the sum term by term
//! keeps that correlation and is what makes the noise fields decide anything at
//! all.
//!
//! # The transcription this depends on, stated as a liability
//!
//! `hash3`, `GRAD12` and `OCTAVE_OFFSET` are `pub(crate)`-invisible — private to
//! `fields/noise.rs` — so the noise enclosure has to **transcribe** them, and a
//! change to any of them in the crate would silently invalidate it. Rather than
//! trust the transcription, [`verify_noise_transcription`] checks it: a
//! reimplemented `perlin`/`fbm` is compared **bit for bit** against
//! `NoiseVolume::sample` and `FbmTerrain::sample` at 68,921 points before any
//! certificate is issued, and the result is on every CSV row as
//! `noise_transcription_verified`. If that column is ever `false` the noise rows
//! are not evidence of anything. The clean remedy is crate-side visibility, which
//! is a source change this ticket does not authorise.
//!
//! # What the three-way verdict means
//!
//! `undecided_fraction` is the ambiguity in the registration, and it is resolved
//! toward the reading that carries information. The certificate has three
//! outcomes per cell, not two:
//!
//! - **certified empty** — `0 ∉ □F(C)`, proven surface-free.
//! - **definitely active** — the eight corners disagree in sign, so by continuity
//!   the surface really is in there. Not a failure of the certificate.
//! - **undecided** — neither. This is the enclosure's slack, and it is the only
//!   number that measures how good the inclusion function is.
//!
//! So `undecided_fraction` is the third of those. `undecided_fraction_naive`
//! (`1 − certified_fraction`, the other reading) is recorded beside it so the
//! artefact carries both.
//!
//! # The comparison that gives it meaning
//!
//! `certified_vs_trilinear` is `IsotopyReport::certified_fraction()` on the same
//! grid: the share of **active** cells that `isotopy.rs` certifies against the
//! interpolant. Note the denominators differ — that one is over active cells,
//! `certified_fraction` is over all cells — so `trilinear_inactive_fraction` is
//! recorded as the apples-to-apples partner: the share of **all** cells whose
//! eight corners agree in sign, which is exactly `isotopy.rs`'s first clause and
//! exactly "the interpolant is provably surface-free here".
//!
//! Those two are nested, always: if `0 ∉ □F(C)` over the whole cell then in
//! particular the eight corners agree, so **certified ⊆ inactive**. The
//! difference, `trust_gap_fraction`, is the set of cells where the corners agree
//! but the analytic field is not proven free — the cells `isotopy.rs` passes over
//! in silence and this certificate refuses. That gap is the quantity this
//! experiment exists to produce.

mod common;

use std::ops::{Add, AddAssign, Mul, Neg, Sub};

use common::experiment::Run;
use isomesh::fields::{
    BoxExact, Difference, FbmTerrain, Gyroid, Intersection, NoiseVolume, ReferenceField, Sphere,
    ThinPlate, Torus, Union,
};
use isomesh::validate::isotopy_report;
use isomesh::{RuntimeShape3, Sdf};

/// The registered resolution: 33 samples, so 32 cells per axis.
const SAMPLES: u32 = 33;

/// Cells per axis.
const CELLS_PER_AXIS: u32 = SAMPLES - 1;

/// Dense-sampling resolution per cell axis. `16³ = 4096`, as registered.
const DENSE: u32 = 16;

const PI: f64 = std::f64::consts::PI;
const TAU: f64 = std::f64::consts::TAU;

/// ULP-scale units of slack for the *field implementation's* own rounding.
///
/// See the module docs: the longest primitive here runs well under a hundred
/// rounded operations on quantities bounded by two.
const IMPL_SLACK_ULPS: f64 = 256.0;

/// Absolute guard for the trig critical-point containment test, in radians.
///
/// Every argument in this experiment is under `2⁴` in magnitude, so this is
/// ~10⁴ ULPs — far more than the test needs, and it only ever widens.
const TRIG_GUARD: f64 = 1e-12;

// ─── the interval ───────────────────────────────────────────────────────────

/// A sound enclosure. Every value the field returns on the box is in `[lo, hi]`.
#[derive(Clone, Copy)]
struct Iv {
    lo: f64,
    hi: f64,
}

impl Iv {
    /// A degenerate interval.
    const fn point(v: f64) -> Self {
        Self { lo: v, hi: v }
    }

    /// An interval from ordered endpoints.
    fn of(lo: f64, hi: f64) -> Self {
        debug_assert!(lo <= hi, "interval endpoints out of order: {lo} > {hi}");
        Self { lo, hi }
    }

    /// One ULP outward at each end. The only thing standing in for directed
    /// rounding.
    fn widen(self) -> Self {
        Self {
            lo: self.lo.next_down(),
            hi: self.hi.next_up(),
        }
    }

    /// Slack for the rounding error of the *field's* own `f64` evaluation, which
    /// outward-rounded interval arithmetic over the exact real function does not
    /// cover.
    fn impl_slack(self) -> Self {
        let m = self.lo.abs().max(self.hi.abs()).max(1.0);
        let e = IMPL_SLACK_ULPS * f64::EPSILON * m;
        Self {
            lo: self.lo - e,
            hi: self.hi + e,
        }
    }

    /// The exact interval square.
    ///
    /// Not `self * self`: for an interval straddling zero that would give a
    /// negative lower bound — sound, but it would then poison [`Iv::sqrt`]'s
    /// radicand and lose the whole point of a tight distance enclosure.
    fn sqr(self) -> Self {
        let (a, b) = (self.lo * self.lo, self.hi * self.hi);
        let out = if self.lo > 0.0 {
            Self { lo: a, hi: b }
        } else if self.hi < 0.0 {
            Self { lo: b, hi: a }
        } else {
            Self {
                lo: 0.0,
                hi: a.max(b),
            }
        };
        out.widen()
    }

    /// Square root over the non-negative part.
    ///
    /// Clamping the lower endpoint at zero is exact rather than a fallback: the
    /// radicand here is a sum of squares, so its true lower bound is never
    /// negative, and only the one-ULP widening can push the endpoint under.
    fn sqrt(self) -> Self {
        let lo = if self.lo > 0.0 { libm::sqrt(self.lo) } else { 0.0 };
        let hi = if self.hi > 0.0 { libm::sqrt(self.hi) } else { 0.0 };
        Self { lo, hi }.widen()
    }

    /// Absolute value. Exact — negation introduces no rounding.
    fn abs(self) -> Self {
        if self.lo >= 0.0 {
            self
        } else if self.hi <= 0.0 {
            Self {
                lo: -self.hi,
                hi: -self.lo,
            }
        } else {
            Self {
                lo: 0.0,
                hi: (-self.lo).max(self.hi),
            }
        }
    }

    /// Interval minimum. Exact: `min` over a box is `min` of the endpoints.
    fn imin(self, o: Self) -> Self {
        Self {
            lo: self.lo.min(o.lo),
            hi: self.hi.min(o.hi),
        }
    }

    /// Interval maximum. Exact, for the same reason.
    fn imax(self, o: Self) -> Self {
        Self {
            lo: self.lo.max(o.lo),
            hi: self.hi.max(o.hi),
        }
    }

    /// Smallest interval containing both.
    fn hull(self, o: Self) -> Self {
        Self {
            lo: self.lo.min(o.lo),
            hi: self.hi.max(o.hi),
        }
    }

    /// `0 ∈ [lo, hi]`. The negation of the certificate.
    fn straddles_zero(self) -> bool {
        self.lo <= 0.0 && self.hi >= 0.0
    }

    fn width(self) -> f64 {
        self.hi - self.lo
    }

    /// `sin` over the interval, by monotonic-branch analysis.
    fn sin(self) -> Self {
        if self.width() >= TAU {
            return Self { lo: -1.0, hi: 1.0 };
        }
        let (a, b) = (libm::sin(self.lo), libm::sin(self.hi));
        let mut out = Self {
            lo: a.min(b),
            hi: a.max(b),
        };
        if self.holds_critical(PI / 2.0) {
            out.hi = 1.0;
        }
        if self.holds_critical(-PI / 2.0) {
            out.lo = -1.0;
        }
        // Three ULPs: `libm`'s sin is not correctly rounded, so an interior point
        // can exceed the endpoint value by up to two of them.
        out.widen().widen().widen()
    }

    /// `cos` over the interval, by monotonic-branch analysis.
    fn cos(self) -> Self {
        if self.width() >= TAU {
            return Self { lo: -1.0, hi: 1.0 };
        }
        let (a, b) = (libm::cos(self.lo), libm::cos(self.hi));
        let mut out = Self {
            lo: a.min(b),
            hi: a.max(b),
        };
        if self.holds_critical(0.0) {
            out.hi = 1.0;
        }
        if self.holds_critical(PI) {
            out.lo = -1.0;
        }
        out.widen().widen().widen()
    }

    /// Is some `c + 2kπ` inside the interval?
    ///
    /// Four candidate `k` around the obvious one, and a guard that errs toward
    /// answering yes — a spurious yes pins an end to `±1`, which widens.
    fn holds_critical(self, c: f64) -> bool {
        let k0 = libm::floor((self.lo - c) / TAU);
        let range = (self.lo - TRIG_GUARD)..=(self.hi + TRIG_GUARD);
        [-1.0, 0.0, 1.0, 2.0]
            .into_iter()
            .any(|dk| range.contains(&(c + (k0 + dk) * TAU)))
    }
}

impl Add for Iv {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self {
            lo: self.lo + o.lo,
            hi: self.hi + o.hi,
        }
        .widen()
    }
}

impl AddAssign for Iv {
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}

impl Sub for Iv {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self {
            lo: self.lo - o.hi,
            hi: self.hi - o.lo,
        }
        .widen()
    }
}

impl Mul for Iv {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let p = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        let mut lo = p[0];
        let mut hi = p[0];
        for &v in &p[1..] {
            lo = lo.min(v);
            hi = hi.max(v);
        }
        Self { lo, hi }.widen()
    }
}

impl Neg for Iv {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
}

/// Euclidean norm of a three-interval vector.
fn norm_iv(v: [Iv; 3]) -> Iv {
    (v[0].sqr() + v[1].sqr() + v[2].sqr()).sqrt()
}

// ─── the inclusion function ─────────────────────────────────────────────────

/// A sound inclusion function for a field.
///
/// `enclose(b)` must contain every value `Sdf::sample` returns for a point in
/// `b`. Everything else in this file depends on that and nothing checks it at
/// compile time, which is why clause one exists.
trait Enclose {
    /// Short tag for the CSV, so a row says which kind of enclosure produced it.
    const KIND: &'static str;

    fn enclose(&self, b: &[Iv; 3]) -> Iv;
}

impl Enclose for Sphere<f64> {
    const KIND: &'static str = "exact_distance";

    fn enclose(&self, b: &[Iv; 3]) -> Iv {
        let d = [
            b[0] - Iv::point(self.center[0]),
            b[1] - Iv::point(self.center[1]),
            b[2] - Iv::point(self.center[2]),
        ];
        (norm_iv(d) - Iv::point(self.radius)).impl_slack()
    }
}

/// `box_sample`'s enclosure, shared by [`BoxExact`] and [`ThinPlate`] exactly as
/// `box_q` is shared by the fields themselves.
fn box_enclose(b: &[Iv; 3], center: [f64; 3], half: [f64; 3]) -> Iv {
    let mut q = [Iv::point(0.0); 3];
    for (k, slot) in q.iter_mut().enumerate() {
        *slot = (b[k] - Iv::point(center[k])).abs() - Iv::point(half[k]);
    }
    let zero = Iv::point(0.0);
    let outside = [q[0].imax(zero), q[1].imax(zero), q[2].imax(zero)];
    let inside = q[0].imax(q[1]).imax(q[2]).imin(zero);
    (norm_iv(outside) + inside).impl_slack()
}

impl Enclose for BoxExact<f64> {
    const KIND: &'static str = "exact_distance";

    fn enclose(&self, b: &[Iv; 3]) -> Iv {
        box_enclose(b, self.center, self.half_extents)
    }
}

impl Enclose for ThinPlate<f64> {
    const KIND: &'static str = "exact_distance";

    fn enclose(&self, b: &[Iv; 3]) -> Iv {
        box_enclose(b, self.center, self.half_extents)
    }
}

impl Enclose for Torus<f64> {
    const KIND: &'static str = "exact_distance";

    fn enclose(&self, b: &[Iv; 3]) -> Iv {
        let d = [
            b[0] - Iv::point(self.center[0]),
            b[1] - Iv::point(self.center[1]),
            b[2] - Iv::point(self.center[2]),
        ];
        let s = (d[0].sqr() + d[2].sqr()).sqrt();
        let q0 = s - Iv::point(self.major);
        let q1 = d[1];
        ((q0.sqr() + q1.sqr()).sqrt() - Iv::point(self.minor)).impl_slack()
    }
}

impl<A: Enclose, B: Enclose> Enclose for Union<A, B> {
    const KIND: &'static str = "csg_min";

    fn enclose(&self, b: &[Iv; 3]) -> Iv {
        self.a.enclose(b).imin(self.b.enclose(b))
    }
}

impl<A: Enclose, B: Enclose> Enclose for Intersection<A, B> {
    const KIND: &'static str = "csg_max";

    fn enclose(&self, b: &[Iv; 3]) -> Iv {
        self.a.enclose(b).imax(self.b.enclose(b))
    }
}

impl<A: Enclose, B: Enclose> Enclose for Difference<A, B> {
    const KIND: &'static str = "csg_max";

    fn enclose(&self, b: &[Iv; 3]) -> Iv {
        self.a.enclose(b).imax(-self.b.enclose(b))
    }
}

impl Enclose for Gyroid<f64> {
    const KIND: &'static str = "trig";

    fn enclose(&self, b: &[Iv; 3]) -> Iv {
        let s = Iv::point(self.scale);
        let (a, bb, c) = (b[0] * s, b[1] * s, b[2] * s);
        let v = a.sin() * bb.cos() + bb.sin() * c.cos() + c.sin() * a.cos();
        (v - Iv::point(self.iso)).impl_slack()
    }
}

impl Enclose for NoiseVolume<f64> {
    const KIND: &'static str = "lattice_noise";

    fn enclose(&self, b: &[Iv; 3]) -> Iv {
        let f = Iv::point(self.frequency);
        let q = [b[0] * f, b[1] * f, b[2] * f];
        (perlin_enclose(&q, self.seed) - Iv::point(self.iso)).impl_slack()
    }
}

impl Enclose for FbmTerrain<f64> {
    const KIND: &'static str = "lattice_noise";

    fn enclose(&self, b: &[Iv; 3]) -> Iv {
        let n = fbm_enclose(
            &[b[0], Iv::point(0.0), b[2]],
            self.seed,
            self.octaves,
            self.lacunarity,
            self.gain,
            self.frequency,
        );
        let height = Iv::point(self.base_height) + Iv::point(self.amplitude) * n;
        (b[1] - height).impl_slack()
    }
}

// ─── the noise, transcribed from `fields/noise.rs` ──────────────────────────

/// Transcribed. A change to `fields/noise.rs::GRAD12` silently invalidates every
/// noise enclosure here, which is why [`verify_noise_transcription`] runs first.
const GRAD12: [[i8; 3]; 12] = [
    [1, 1, 0],
    [-1, 1, 0],
    [1, -1, 0],
    [-1, -1, 0],
    [1, 0, 1],
    [-1, 0, 1],
    [1, 0, -1],
    [-1, 0, -1],
    [0, 1, 1],
    [0, -1, 1],
    [0, 1, -1],
    [0, -1, -1],
];

/// Transcribed: `octave × (1/φ, 1/φ², 1/φ³)`.
const OCTAVE_OFFSET: [f64; 3] = [
    0.618_033_988_749_894_8,
    0.381_966_011_250_105_15,
    0.236_067_977_499_789_7,
];

/// Transcribed lattice hash.
const fn hash3(ix: i32, iy: i32, iz: i32, seed: u32) -> u32 {
    let mut h = seed;
    h = h.wrapping_add((ix as u32).wrapping_mul(0x9E37_79B1));
    h = h.wrapping_add((iy as u32).wrapping_mul(0x85EB_CA6B));
    h = h.wrapping_add((iz as u32).wrapping_mul(0xC2B2_AE35));
    h ^= h >> 16;
    h = h.wrapping_mul(0x85EB_CA6B);
    h ^= h >> 13;
    h = h.wrapping_mul(0xC2B2_AE35);
    h ^= h >> 16;
    h
}

/// Transcribed: the crate narrows through `f32` before the integer cast, and the
/// hash input has to match bit for bit.
fn lattice_index(floored: f64) -> i32 {
    floored as f32 as i32
}

/// Perlin's quintic fade, in the crate's operation order.
fn fade_point(t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    t3 * (t * (t * 6.0 - 15.0) + 10.0)
}

/// The fade over an interval.
///
/// `u' = 30t²(t − 1)² ≥ 0` for **every** real `t`, so the quintic is monotone
/// increasing and its endpoints give its range with no case analysis. The result
/// is clamped to `[0, 1]`, which is exact because the true argument is a
/// fractional part.
fn fade_iv(t: Iv) -> Iv {
    let lo = fade_point(t.lo);
    let hi = fade_point(t.hi);
    let w = Iv::of(lo.min(hi), lo.max(hi)).widen().widen();
    let clamped_lo = w.lo.max(0.0);
    Iv::of(clamped_lo, w.hi.min(1.0).max(clamped_lo))
}

/// `perlin`'s value, reimplemented so the transcription can be checked.
fn perlin_point(p: [f64; 3], seed: u32) -> f64 {
    let base = [libm::floor(p[0]), libm::floor(p[1]), libm::floor(p[2])];
    let i = [
        lattice_index(base[0]),
        lattice_index(base[1]),
        lattice_index(base[2]),
    ];
    let t = [p[0] - base[0], p[1] - base[1], p[2] - base[2]];
    let (ux, uy, uz) = (fade_point(t[0]), fade_point(t[1]), fade_point(t[2]));

    let mut value = 0.0_f64;
    for cz in 0..2i32 {
        for cy in 0..2i32 {
            for cx in 0..2i32 {
                let g = GRAD12[(hash3(i[0] + cx, i[1] + cy, i[2] + cz, seed) % 12) as usize];
                let gv = [f64::from(g[0]), f64::from(g[1]), f64::from(g[2])];
                let d = [
                    t[0] - f64::from(cx),
                    t[1] - f64::from(cy),
                    t[2] - f64::from(cz),
                ];
                let dot = gv[0] * d[0] + gv[1] * d[1] + gv[2] * d[2];
                let wx = if cx == 1 { ux } else { 1.0 - ux };
                let wy = if cy == 1 { uy } else { 1.0 - uy };
                let wz = if cz == 1 { uz } else { 1.0 - uz };
                value += wx * wy * wz * dot;
            }
        }
    }
    value
}

/// `fbm`'s value, reimplemented for the same reason.
fn fbm_point(p: [f64; 3], seed: u32, octaves: u32, lacunarity: f64, gain: f64, frequency: f64) -> f64 {
    let mut value = 0.0_f64;
    let mut freq = frequency;
    let mut amp = 1.0_f64;
    for octave in 0..octaves {
        let k = f64::from(octave);
        let q = [
            p[0] * freq + k * OCTAVE_OFFSET[0],
            p[1] * freq + k * OCTAVE_OFFSET[1],
            p[2] * freq + k * OCTAVE_OFFSET[2],
        ];
        value += amp * perlin_point(q, seed);
        freq *= lacunarity;
        amp *= gain;
    }
    value
}

/// `g · d` for one axis. Exact: every `GRAD12` component is `0` or `±1`.
fn axis_term(g: i8, d: Iv) -> Iv {
    match g {
        1 => d,
        -1 => -d,
        _ => Iv::point(0.0),
    }
}

/// `perlin` over a box lying inside one lattice cell, term by term.
fn perlin_cell_enclose(q: &[Iv; 3], i: [i32; 3], seed: u32) -> Iv {
    let mut t = [Iv::point(0.0); 3];
    for (k, slot) in t.iter_mut().enumerate() {
        let base = f64::from(i[k]);
        let lo = (q[k].lo - base).max(0.0);
        let hi = (q[k].hi - base).min(1.0);
        assert!(
            lo <= hi,
            "lattice cell {base} does not overlap [{}, {}]",
            q[k].lo,
            q[k].hi
        );
        *slot = Iv::of(lo, hi).widen();
    }
    let u = [fade_iv(t[0]), fade_iv(t[1]), fade_iv(t[2])];
    let one = Iv::point(1.0);

    let mut value = Iv::point(0.0);
    for cz in 0..2i32 {
        for cy in 0..2i32 {
            for cx in 0..2i32 {
                let g = GRAD12[(hash3(i[0] + cx, i[1] + cy, i[2] + cz, seed) % 12) as usize];
                let d = [
                    t[0] - Iv::point(f64::from(cx)),
                    t[1] - Iv::point(f64::from(cy)),
                    t[2] - Iv::point(f64::from(cz)),
                ];
                let dot = axis_term(g[0], d[0]) + axis_term(g[1], d[1]) + axis_term(g[2], d[2]);
                let wx = if cx == 1 { u[0] } else { one - u[0] };
                let wy = if cy == 1 { u[1] } else { one - u[1] };
                let wz = if cz == 1 { u[2] } else { one - u[2] };
                value += wx * wy * wz * dot;
            }
        }
    }
    value
}

/// `perlin` over an arbitrary box: the hull over every lattice cell it meets.
fn perlin_enclose(q: &[Iv; 3], seed: u32) -> Iv {
    let ilo = [
        lattice_index(libm::floor(q[0].lo)),
        lattice_index(libm::floor(q[1].lo)),
        lattice_index(libm::floor(q[2].lo)),
    ];
    let ihi = [
        lattice_index(libm::floor(q[0].hi)),
        lattice_index(libm::floor(q[1].hi)),
        lattice_index(libm::floor(q[2].hi)),
    ];
    let mut out: Option<Iv> = None;
    for iz in ilo[2]..=ihi[2] {
        for iy in ilo[1]..=ihi[1] {
            for ix in ilo[0]..=ihi[0] {
                let cell = perlin_cell_enclose(q, [ix, iy, iz], seed);
                out = Some(match out {
                    None => cell,
                    Some(prev) => prev.hull(cell),
                });
            }
        }
    }
    out.expect("every box meets at least one lattice cell")
}

/// `fbm` over a box: the sum of the per-octave enclosures.
fn fbm_enclose(
    b: &[Iv; 3],
    seed: u32,
    octaves: u32,
    lacunarity: f64,
    gain: f64,
    frequency: f64,
) -> Iv {
    let mut value = Iv::point(0.0);
    let mut freq = frequency;
    let mut amp = 1.0_f64;
    for octave in 0..octaves {
        let k = f64::from(octave);
        let f = Iv::point(freq);
        let q = [
            b[0] * f + Iv::point(k * OCTAVE_OFFSET[0]),
            b[1] * f + Iv::point(k * OCTAVE_OFFSET[1]),
            b[2] * f + Iv::point(k * OCTAVE_OFFSET[2]),
        ];
        value += Iv::point(amp) * perlin_enclose(&q, seed);
        freq *= lacunarity;
        amp *= gain;
    }
    value
}

/// Bit-for-bit check of the transcription against the crate's own fields.
///
/// Returns `(points, matches, max_abs_diff)`. Anything but a full match makes the
/// two noise rows worthless, and the columns say so on every row rather than in a
/// comment.
fn verify_noise_transcription() -> (u64, u64, f64) {
    let nv = NoiseVolume::<f64>::canonical();
    let terrain = FbmTerrain::<f64>::canonical();
    let mut points = 0u64;
    let mut matches = 0u64;
    let mut worst = 0.0_f64;

    // A deterministic 41³ lattice over each field's own domain, deliberately not
    // aligned with the 33-sample grid so the check does not only see the
    // coordinates the certificate happens to use.
    for (lo, hi, which) in [(-2.0_f64, 2.0_f64, 0u8), (-8.0, 8.0, 1)] {
        let step = (hi - lo) / 40.0;
        for iz in 0..41u32 {
            for iy in 0..41u32 {
                for ix in 0..41u32 {
                    let p = [
                        lo + f64::from(ix) * step,
                        lo + f64::from(iy) * step,
                        lo + f64::from(iz) * step,
                    ];
                    let (theirs, mine) = if which == 0 {
                        let q = [p[0] * nv.frequency, p[1] * nv.frequency, p[2] * nv.frequency];
                        (nv.sample(p), perlin_point(q, nv.seed) - nv.iso)
                    } else {
                        let n = fbm_point(
                            [p[0], 0.0, p[2]],
                            terrain.seed,
                            terrain.octaves,
                            terrain.lacunarity,
                            terrain.gain,
                            terrain.frequency,
                        );
                        (
                            terrain.sample(p),
                            p[1] - (terrain.base_height + terrain.amplitude * n),
                        )
                    };
                    points += 1;
                    if theirs.to_bits() == mine.to_bits() {
                        matches += 1;
                    }
                    worst = worst.max((theirs - mine).abs());
                }
            }
        }
    }
    (points, matches, worst)
}

// ─── the experiment ─────────────────────────────────────────────────────────

/// What one field's grid produced.
struct Tally {
    cells: u64,
    certified: u64,
    definitely_active: u64,
    undecided: u64,
    dense_free: u64,
    dense_free_and_certified: u64,
    unsound: u64,
    width_sum: f64,
    first_unsound: Option<[u32; 3]>,
}

/// Min and max of `DENSE³` samples over the closed cell.
///
/// `t = i / (DENSE − 1)` is exactly `0` and `1` at the ends, so every sample lies
/// in the closed cell the enclosure was asked about — which is what makes an
/// unsound certification a genuine contradiction rather than an off-by-one.
fn dense_min_max<F: Sdf<Scalar = f64>>(field: &F, origin: [f64; 3], h: f64) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    let last = f64::from(DENSE - 1);
    for iz in 0..DENSE {
        let tz = f64::from(iz) / last;
        for iy in 0..DENSE {
            let ty = f64::from(iy) / last;
            for ix in 0..DENSE {
                let tx = f64::from(ix) / last;
                let v = field.sample([
                    origin[0] + tx * h,
                    origin[1] + ty * h,
                    origin[2] + tz * h,
                ]);
                lo = lo.min(v);
                hi = hi.max(v);
            }
        }
    }
    (lo, hi)
}

/// The whole grid for one field.
fn tally<F>(field: &F, lo: [f64; 3], h: f64) -> Tally
where
    F: Sdf<Scalar = f64> + Enclose,
{
    let mut out = Tally {
        cells: 0,
        certified: 0,
        definitely_active: 0,
        undecided: 0,
        dense_free: 0,
        dense_free_and_certified: 0,
        unsound: 0,
        width_sum: 0.0,
        first_unsound: None,
    };

    for cz in 0..CELLS_PER_AXIS {
        for cy in 0..CELLS_PER_AXIS {
            for cx in 0..CELLS_PER_AXIS {
                let origin = [
                    lo[0] + f64::from(cx) * h,
                    lo[1] + f64::from(cy) * h,
                    lo[2] + f64::from(cz) * h,
                ];
                let b = [
                    Iv::of(origin[0], origin[0] + h).widen(),
                    Iv::of(origin[1], origin[1] + h).widen(),
                    Iv::of(origin[2], origin[2] + h).widen(),
                ];
                let e = field.enclose(&b);
                let certified = !e.straddles_zero();

                // Provably active: the eight corners disagree, so continuity puts
                // the surface inside. Needs no enclosure at all.
                let mut neg = false;
                let mut non_neg = false;
                for corner in 0..8u32 {
                    let v = field.sample([
                        origin[0] + f64::from(corner & 1) * h,
                        origin[1] + f64::from((corner >> 1) & 1) * h,
                        origin[2] + f64::from((corner >> 2) & 1) * h,
                    ]);
                    if v < 0.0 {
                        neg = true;
                    } else {
                        non_neg = true;
                    }
                }
                let definitely_active = neg && non_neg;

                let (dmin, dmax) = dense_min_max(field, origin, h);
                let dense_free = dmin > 0.0 || dmax < 0.0;

                out.cells += 1;
                out.width_sum += e.width();
                if certified {
                    out.certified += 1;
                    if dense_free {
                        out.dense_free_and_certified += 1;
                    } else {
                        out.unsound += 1;
                        if out.first_unsound.is_none() {
                            out.first_unsound = Some([cx, cy, cz]);
                        }
                    }
                } else if definitely_active {
                    out.definitely_active += 1;
                } else {
                    out.undecided += 1;
                }
                if dense_free {
                    out.dense_free += 1;
                }
            }
        }
    }
    out
}

/// Provenance of the noise transcription, carried onto every row.
struct Provenance {
    points: u64,
    matches: u64,
    worst: f64,
}

/// One field: sample, certify, densely check, record.
fn sweep<F>(run: &mut Run, name: &str, field: &F, prov: &Provenance)
where
    F: ReferenceField + Sdf<Scalar = f64> + Enclose,
{
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(SAMPLES - 1);
    let shape = RuntimeShape3::new([SAMPLES; 3]).expect("33³ fits u32");

    let n = SAMPLES as usize;
    let mut values = Vec::with_capacity(n * n * n);
    for iz in 0..SAMPLES {
        for iy in 0..SAMPLES {
            for ix in 0..SAMPLES {
                values.push(field.sample([
                    lo[0] + f64::from(ix) * h,
                    lo[1] + f64::from(iy) * h,
                    lo[2] + f64::from(iz) * h,
                ]));
            }
        }
    }
    let tri = isotopy_report(&values, &shape).expect("33³ grid is large enough");
    let t = tally(field, lo, h);

    let cells = t.cells as f64;
    let certified_fraction = t.certified as f64 / cells;
    let undecided_fraction = t.undecided as f64 / cells;
    let trilinear_inactive = (t.cells - tri.active_cells) as f64 / cells;
    let of_sampled_free = if t.dense_free == 0 {
        1.0
    } else {
        t.dense_free_and_certified as f64 / t.dense_free as f64
    };

    println!(
        "{name:>14}  certified {:>5}/{:<5} ({certified_fraction:6.4})  \
         dense-free {:>5} → {of_sampled_free:6.4} of them  undecided {undecided_fraction:6.4}  \
         active {:>5}  UNSOUND {}  trilinear inactive {trilinear_inactive:6.4} \
         (gap {:+.4})  mean width {:.4}",
        t.certified,
        t.cells,
        t.dense_free,
        t.definitely_active,
        t.unsound,
        trilinear_inactive - certified_fraction,
        t.width_sum / cells,
    );

    run.record(&[
        ("field", name.to_string()),
        ("samples_per_axis", SAMPLES.to_string()),
        ("cells", t.cells.to_string()),
        (
            "cells_surface_free_sampled",
            t.dense_free.to_string(),
        ),
        ("cells_certified_empty", t.certified.to_string()),
        ("certified_fraction", format!("{certified_fraction:.9}")),
        ("unsound_certifications", t.unsound.to_string()),
        ("undecided_fraction", format!("{undecided_fraction:.9}")),
        (
            "certified_vs_trilinear",
            format!("{:.9}", tri.certified_fraction()),
        ),
        // ── clause two's actual quantity ────────────────────────────────────
        (
            "certified_of_sampled_free",
            format!("{of_sampled_free:.9}"),
        ),
        (
            "certified_and_sampled_free",
            t.dense_free_and_certified.to_string(),
        ),
        // ── the three-way verdict, and the other reading of "undecided" ─────
        (
            "cells_definitely_active",
            t.definitely_active.to_string(),
        ),
        ("cells_undecided", t.undecided.to_string()),
        (
            "undecided_fraction_naive",
            format!("{:.9}", 1.0 - certified_fraction),
        ),
        // ── the apples-to-apples trilinear comparison ──────────────────────
        (
            "trilinear_inactive_fraction",
            format!("{trilinear_inactive:.9}"),
        ),
        (
            "trust_gap_fraction",
            format!("{:.9}", trilinear_inactive - certified_fraction),
        ),
        ("trilinear_active_cells", tri.active_cells.to_string()),
        ("trilinear_certified", tri.certified.to_string()),
        ("trilinear_uncertified", tri.uncertified.to_string()),
        // ── the enclosure itself ───────────────────────────────────────────
        ("enclosure_kind", F::KIND.to_string()),
        (
            "mean_enclosure_width",
            format!("{:.9}", t.width_sum / cells),
        ),
        (
            "mean_enclosure_width_over_h",
            format!("{:.9}", t.width_sum / cells / h),
        ),
        ("cell_size", format!("{h:.9}")),
        ("domain_half_extent", format!("{:.6}", hi[0])),
        (
            "dense_points_per_cell",
            (DENSE * DENSE * DENSE).to_string(),
        ),
        (
            "first_unsound_cell",
            match t.first_unsound {
                Some([x, y, z]) => format!("{x}:{y}:{z}"),
                None => "none".to_string(),
            },
        ),
        // ── the transcription this row's soundness may rest on ─────────────
        (
            "noise_transcription_verified",
            (prov.matches == prov.points).to_string(),
        ),
        (
            "noise_transcription_points",
            prov.points.to_string(),
        ),
        (
            "noise_transcription_max_abs_diff",
            format!("{:.3e}", prov.worst),
        ),
    ]);
}

fn main() {
    let prereg = isomesh::experiment!("P-48");

    let (points, matches, worst) = verify_noise_transcription();
    println!(
        "noise transcription: {matches}/{points} bit-identical, max |Δ| = {worst:.3e}\n"
    );
    let prov = Provenance {
        points,
        matches,
        worst,
    };

    common::experiment::run(prereg, |run| {
        isomesh::for_each_reference_field!(f64, |name, field| {
            sweep(run, name, &field, &prov);
        });
    });
}
