//! The scalar type the crate is generic over.

use core::cmp::Ordering;
use core::fmt::Debug;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

mod sealed {
    pub trait Sealed {}
}

/// The scalar an [`Sdf`](crate::Sdf) samples in and the vertex solvers work in.
///
/// Implemented for `f32` and `f64`, and **sealed** — downstream crates cannot
/// add implementations. Fixed-point and interval scalars are deliberately out of
/// scope: the extraction algorithms assume IEEE-754 semantics (signed zero,
/// infinities, a total order) in enough places that a third implementation would
/// be a second execution path rather than a new type.
///
/// # Precision policy
///
/// `f32` is the default everywhere. `f64` exists because `M = AᵀA` squares the
/// condition number in the dual-contouring solve, and because CAD consumers work
/// at coordinates where `f32` has already lost the surface. The crate never
/// narrows a scalar internally; narrowing happens only where a caller asks for
/// it, through [`Real::as_f32`].
///
/// Note the direction of the conversions: [`from_f64`](Real::from_f64) builds a
/// scalar at the target's own precision, and [`as_f32`](Real::as_f32) narrows on
/// the way out to a caller who asked for it. There is deliberately **no**
/// `Self: Into<f32>` bound — that is what disqualifies `fast-surface-nets` for
/// CAD, and having it would foreclose `f64` here for exactly the same reason.
///
/// # Float backend
///
/// [`sqrt`](Real::sqrt), [`floor`](Real::floor), [`sin`](Real::sin) and
/// [`cos`](Real::cos) come from `libm` unconditionally — see the crate-level
/// docs for why. `sqrt` and `floor` are correctly rounded under IEEE-754, so
/// they are bit-identical everywhere regardless. `sin` and `cos` are not
/// specified to be correctly rounded by anyone, which is precisely why routing
/// them through a pure-Rust implementation rather than the platform's libm is
/// what makes committed golden hashes portable.
///
/// # Not here, on purpose
///
/// `mul_add` is absent. A fused multiply-add gives different results from a
/// separate multiply and add, and GPU-005 requires the compute path to be
/// bit-identical to this one. If a call site later proves it needs FMA, add it
/// *and* amend GPU-005's acceptance criterion in the same commit.
#[allow(private_bounds)] // The seal is real: `sealed::Sealed` is unnameable downstream.
pub trait Real:
    sealed::Sealed
    + Copy
    + Debug
    + Default
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Mul<Output = Self>
    + MulAssign
    + Div<Output = Self>
    + DivAssign
    + Neg<Output = Self>
    + Send
    + Sync
    + 'static
{
    /// `0.0`.
    const ZERO: Self;
    /// `1.0`.
    const ONE: Self;
    /// `2.0`. Present because the central-difference denominator `2h` is on a
    /// hot path and `ONE + ONE` reads worse than it computes.
    const TWO: Self;
    /// `0.5`.
    const HALF: Self;
    /// Machine epsilon: the gap between `1.0` and the next representable value.
    const EPSILON: Self;
    /// Positive infinity. Initialises bounding boxes and minimum searches.
    const INFINITY: Self;
    /// Negative infinity.
    const NEG_INFINITY: Self;

    /// Step size for the central-difference gradient in
    /// [`Sdf::gradient`](crate::Sdf::gradient), for a field of unit scale.
    ///
    /// Equal to `cbrt(EPSILON)`. A central difference carries truncation error
    /// `≈ h²·f‴/6` and round-off error `≈ EPSILON·|f|/h`; their sum is minimised
    /// at `h = cbrt(EPSILON)` when `f` and `f‴` are `O(1)` in the units of `p`.
    ///
    /// The default gradient scales this by coordinate magnitude — see
    /// [`Sdf::gradient`](crate::Sdf::gradient). A field whose characteristic
    /// length is far from `1` should override `gradient` rather than expect this
    /// constant to suit it.
    const DIFF_STEP: Self;

    // ── from libm ───────────────────────────────────────────────────────────

    /// Square root. Correctly rounded under IEEE-754, hence identical on every
    /// conforming platform.
    #[must_use]
    fn sqrt(self) -> Self;

    /// Largest integer not greater than `self`. Exact under IEEE-754.
    #[must_use]
    fn floor(self) -> Self;

    /// Sine, in radians.
    ///
    /// **Not** correctly rounded — no implementation guarantees that. This one
    /// is `libm`'s, which is pure Rust and therefore identical across platforms;
    /// the platform's own libm is not.
    #[must_use]
    fn sin(self) -> Self;

    /// Cosine, in radians. Same caveat as [`Real::sin`].
    #[must_use]
    fn cos(self) -> Self;

    // ── from core ───────────────────────────────────────────────────────────

    /// Absolute value.
    #[must_use]
    fn abs(self) -> Self;

    /// `1.0`, `-1.0`, or `NaN`. Note that `signum(-0.0)` is `-1.0`.
    #[must_use]
    fn signum(self) -> Self;

    /// `1.0 / self`.
    #[must_use]
    fn recip(self) -> Self;

    /// IEEE-754 `minNum`: returns the other operand if one of them is `NaN`.
    #[must_use]
    fn min(self, other: Self) -> Self;

    /// IEEE-754 `maxNum`: returns the other operand if one of them is `NaN`.
    #[must_use]
    fn max(self, other: Self) -> Self;

    /// Restrict to `[low, high]`.
    ///
    /// # Panics
    ///
    /// If `low > high`, or if either bound is `NaN`.
    #[must_use]
    fn clamp(self, low: Self, high: Self) -> Self;

    /// `true` if neither infinite nor `NaN`.
    fn is_finite(self) -> bool;

    /// The IEEE-754 `totalOrder` predicate.
    ///
    /// **This is the only ordering the crate may sort by.** `partial_cmp`
    /// returns `None` on `NaN`, so every `unwrap` of it is a panic path; more
    /// importantly a non-total order makes sort results depend on input order,
    /// which is exactly the determinism leak T-004 exists to catch. Distinguishes
    /// `-0.0 < 0.0`.
    fn total_cmp(&self, other: &Self) -> Ordering;

    // ── conversion ──────────────────────────────────────────────────────────

    /// Build a scalar from a literal, **at this type's own precision**.
    ///
    /// This is how every constant in the crate is written. The obvious
    /// alternative — a `From<f32>` supertrait — silently stores `f32`-accurate
    /// constants in an `f64` field, so a torus asked for in `f64` would come back
    /// with a minor radius of `0.30000001192092896`. That is a CAD consumer
    /// getting `f32` geometry after explicitly asking not to.
    ///
    /// For `f64` this is the identity. For `f32` it rounds the literal once, to
    /// the nearest `f32` — which is the correct constant for that precision, not
    /// a loss of one.
    #[must_use]
    fn from_f64(value: f64) -> Self;

    /// Narrow to `f32`.
    ///
    /// **Lossy for `f64`**, and infallible: values outside `f32`'s range become
    /// infinite and roughly 29 bits of mantissa are discarded. This is the only
    /// narrowing operation in the crate, and the crate itself never calls it —
    /// it exists for consumers writing into an `f32` vertex buffer. A CAD
    /// consumer that needs `f64` output should use a
    /// [`MeshSink`](crate::MeshSink) with `Scalar = f64` and never call this.
    #[must_use]
    fn as_f32(self) -> f32;
}

macro_rules! impl_real {
    (
        $ty:ty,
        diff_step = $diff_step:expr,
        sqrt = $sqrt:path, floor = $floor:path, sin = $sin:path, cos = $cos:path $(,)?
    ) => {
        impl sealed::Sealed for $ty {}

        impl Real for $ty {
            const ZERO: Self = 0.0;
            const ONE: Self = 1.0;
            const TWO: Self = 2.0;
            const HALF: Self = 0.5;
            const EPSILON: Self = <$ty>::EPSILON;
            const INFINITY: Self = <$ty>::INFINITY;
            const NEG_INFINITY: Self = <$ty>::NEG_INFINITY;
            const DIFF_STEP: Self = $diff_step;

            #[inline]
            fn sqrt(self) -> Self {
                $sqrt(self)
            }
            #[inline]
            fn floor(self) -> Self {
                $floor(self)
            }
            #[inline]
            fn sin(self) -> Self {
                $sin(self)
            }
            #[inline]
            fn cos(self) -> Self {
                $cos(self)
            }

            #[inline]
            fn abs(self) -> Self {
                <$ty>::abs(self)
            }
            #[inline]
            fn signum(self) -> Self {
                <$ty>::signum(self)
            }
            #[inline]
            fn recip(self) -> Self {
                <$ty>::recip(self)
            }
            #[inline]
            fn min(self, other: Self) -> Self {
                <$ty>::min(self, other)
            }
            #[inline]
            fn max(self, other: Self) -> Self {
                <$ty>::max(self, other)
            }
            #[inline]
            fn clamp(self, low: Self, high: Self) -> Self {
                <$ty>::clamp(self, low, high)
            }
            #[inline]
            fn is_finite(self) -> bool {
                <$ty>::is_finite(self)
            }
            #[inline]
            fn total_cmp(&self, other: &Self) -> Ordering {
                <$ty>::total_cmp(self, other)
            }

            #[inline]
            #[allow(clippy::unnecessary_cast)] // Identity for f64; rounds once for f32.
            fn from_f64(value: f64) -> Self {
                value as $ty
            }

            #[inline]
            #[allow(clippy::unnecessary_cast)] // A no-op for f32; the narrowing point for f64.
            fn as_f32(self) -> f32 {
                self as f32
            }
        }
    };
}

impl_real!(
    f32,
    diff_step = 4.921_567e-3,
    sqrt = libm::sqrtf,
    floor = libm::floorf,
    sin = libm::sinf,
    cos = libm::cosf,
);

impl_real!(
    f64,
    diff_step = 6.055_454_452_393_343e-6,
    sqrt = libm::sqrt,
    floor = libm::floor,
    sin = libm::sin,
    cos = libm::cos,
);

#[cfg(test)]
mod tests {
    // Several of these tests assert bit-exactness on purpose -- that IEEE-754
    // makes `sqrt` correctly rounded, that the `Real` constants are the literals
    // they claim to be, that `as_f32` narrows the way it says. An approximate
    // comparison would not be a weaker test, it would be a different one.
    #![allow(clippy::float_cmp)]

    use super::*;

    #[test]
    fn constants_are_exact() {
        assert_eq!(f32::ZERO, 0.0);
        assert_eq!(f32::ONE, 1.0);
        assert_eq!(f32::TWO, 2.0);
        assert_eq!(f32::HALF, 0.5);
        assert_eq!(<f32 as Real>::EPSILON, f32::EPSILON);
        assert!(<f32 as Real>::INFINITY.is_infinite());

        assert_eq!(f64::ZERO, 0.0);
        assert_eq!(f64::ONE, 1.0);
        assert_eq!(f64::TWO, 2.0);
        assert_eq!(f64::HALF, 0.5);
        assert_eq!(<f64 as Real>::EPSILON, f64::EPSILON);
        assert!(<f64 as Real>::NEG_INFINITY.is_infinite());
    }

    /// `DIFF_STEP` is a hand-written literal, so a transcription error would
    /// otherwise be invisible: the gradient would still look roughly right and
    /// just be less accurate. Pin it to within a factor of two of `cbrt(EPSILON)`.
    #[test]
    fn diff_step_is_the_cube_root_of_epsilon() {
        let cube = f64::from(f32::DIFF_STEP).powi(3);
        let eps = f64::from(f32::EPSILON);
        assert!(
            cube > 0.5 * eps && cube < 2.0 * eps,
            "f32: {cube:e} vs {eps:e}"
        );

        let cube = f64::DIFF_STEP.powi(3);
        let eps = <f64 as Real>::EPSILON;
        assert!(
            cube > 0.5 * eps && cube < 2.0 * eps,
            "f64: {cube:e} vs {eps:e}"
        );
    }

    /// Reference values come from `core::consts` and from literals — never from
    /// `f32::sqrt`. Comparing `Real::sqrt` against `f32::sqrt` would be a
    /// tautology that passes even if the trait dispatched somewhere else
    /// entirely; `core`'s constants are an independent source.
    ///
    /// `sqrt` is asserted exactly because IEEE-754 requires it to be correctly
    /// rounded, so any conforming implementation returns this bit pattern.
    #[test]
    fn backend_matches_reference_values() {
        assert_eq!(Real::sqrt(2.0f32), core::f32::consts::SQRT_2);
        assert_eq!(Real::sqrt(2.0f64), core::f64::consts::SQRT_2);
        assert_eq!(Real::floor(-1.5f32), -2.0);
        assert_eq!(Real::floor(-1.5f64), -2.0);
        assert_eq!(Real::abs(-0.0f32), 0.0);

        // sin/cos are not correctly rounded by anyone; 4 ULP is generous and
        // still tight enough to catch a swapped or mis-scaled implementation.
        assert!((Real::sin(1.0f32) - 0.841_470_98_f32).abs() <= 4.0 * f32::EPSILON);
        assert!((Real::cos(1.0f32) - 0.540_302_3_f32).abs() <= 4.0 * f32::EPSILON);
        assert!((Real::sin(1.0f64) - 0.841_470_984_807_896_5).abs() <= 4.0 * f64::EPSILON);
        assert!((Real::cos(1.0f64) - 0.540_302_305_868_139_7).abs() <= 4.0 * f64::EPSILON);
    }

    /// The reason `from_f64` exists rather than a `From<f32>` supertrait: a
    /// literal must land at the *target's* precision. With `From<f32>` the `f64`
    /// column below would read `0.30000001192092896`, so an `f64` field would
    /// carry `f32`-accurate geometry — precisely what a CAD consumer asked not to
    /// get by choosing `f64`.
    #[test]
    fn from_f64_builds_at_the_targets_precision() {
        fn literal<R: Real>() -> R {
            R::from_f64(0.3)
        }
        assert_eq!(literal::<f64>(), 0.3f64);
        assert_eq!(literal::<f32>(), 0.3f32);
        assert_ne!(literal::<f64>(), f64::from(0.3f32));

        // Exact values are exact in both.
        assert_eq!(<f64 as Real>::from_f64(0.75), 0.75f64);
        assert_eq!(<f32 as Real>::from_f64(0.75), 0.75f32);
    }

    /// Pins the documented lossy, infallible behaviour of the one narrowing
    /// operation in the crate.
    #[test]
    fn as_f32_narrows_and_saturates() {
        assert!(Real::as_f32(1e300f64).is_infinite());
        assert_eq!(Real::as_f32(1.5f64), 1.5f32);
        assert_eq!(Real::as_f32(1.5f32), 1.5f32);
    }

    /// The ordering used for every deterministic sort in the crate must be
    /// total, which `partial_cmp` is not.
    #[test]
    fn total_cmp_separates_signed_zero() {
        assert_eq!(Real::total_cmp(&-0.0f32, &0.0f32), Ordering::Less);
        assert_eq!(Real::total_cmp(&-0.0f64, &0.0f64), Ordering::Less);
        assert_eq!((-0.0f32).partial_cmp(&0.0f32), Some(Ordering::Equal));
    }
}
