//! The scalar field an extraction algorithm sees.

use alloc::boxed::Box;

use crate::Real;

/// A scalar field sampled at a point. This is the only thing an extraction
/// algorithm needs to see.
///
/// # Sign convention
///
/// **Negative is inside.** [`sample`](Sdf::sample) returns a negative value
/// strictly inside the solid, a positive value strictly outside, and the surface
/// is the zero level set. Half of all "my mesh is inside out" bugs are this
/// convention flipping across a module boundary, so it is stated here, in the
/// crate-level docs, and on every reference field that implements this trait.
///
/// # Exactness is not required
///
/// Nothing here demands `|∇f| == 1`. A true signed distance field satisfies it,
/// but `gyroid` is an implicit function and `fbm_terrain` is a heightfield, and
/// both are legitimate inputs. Algorithms that need a genuine distance must say
/// so, and the reference fields advertise whether they provide one.
///
/// # `dyn` use
///
/// `Sdf` is dyn-compatible once the scalar is named: `&dyn Sdf<Scalar = f32>`
/// and `Box<dyn Sdf<Scalar = f32>>` are valid types. The blanket
/// `impl Sdf for &S` below means such a reference *also* satisfies the generic
/// `S: Sdf` bound that every extraction function uses, so a runtime-selected
/// field needs no special path.
///
/// Prefer the generic form regardless. [`sample`](Sdf::sample) is the innermost
/// call in the crate — a 64³ grid is a quarter of a million evaluations before
/// any refinement — and the generic form inlines it where `dyn` cannot.
pub trait Sdf {
    /// The scalar this field is sampled in.
    type Scalar: Real;

    /// The field value at `p`. Negative inside, positive outside.
    fn sample(&self, p: [Self::Scalar; 3]) -> Self::Scalar;

    /// The **raw, unnormalised** gradient `∇f` at `p`.
    ///
    /// # Direction and magnitude
    ///
    /// Points along increasing [`sample`](Sdf::sample), i.e. **away from the
    /// solid**. For an exact signed distance field `|∇f| == 1` and this is the
    /// outward unit normal, but the crate does not require exactness, and for
    /// `gyroid` and `fbm_terrain` the magnitude is the local Lipschitz constant.
    ///
    /// It is returned raw because raw is the one from which the other is
    /// recoverable: a caller wanting a normal divides by the length, while a
    /// caller refining a root along an edge needs the length itself. Normalising
    /// here would destroy information and add a second failure mode wherever
    /// `|∇f| → 0`.
    ///
    /// # Default implementation
    ///
    /// Central differences — six calls to [`sample`](Sdf::sample) — with
    ///
    /// ```text
    /// h = Scalar::DIFF_STEP * max(|pₓ|, |p_y|, |p_z|, 1)
    /// ```
    ///
    /// The same `h` on all three axes, so the stencil is isotropic and the
    /// returned *direction* is unbiased; a per-axis `h` would bias every normal
    /// and every dual-contouring plane.
    ///
    /// The magnitude scaling is load-bearing rather than decorative. At
    /// `p ≈ [1e6, 0, 0]` in `f32` one ULP is `0.0625`, so a fixed
    /// `h = 4.9e-3` would give `p + h == p` exactly and the gradient would
    /// collapse to `[0, 0, 0]` — a bug in the differencing rather than a
    /// demonstration of anything about `f32`.
    ///
    /// Override this whenever an analytic gradient exists. It is six samples
    /// against one evaluation, and it is exact rather than `O(h²)`.
    #[inline]
    fn gradient(&self, p: [Self::Scalar; 3]) -> [Self::Scalar; 3] {
        let scale = p[0]
            .abs()
            .max(p[1].abs())
            .max(p[2].abs())
            .max(Self::Scalar::ONE);
        let h = Self::Scalar::DIFF_STEP * scale;
        let inv = (Self::Scalar::TWO * h).recip();
        [
            (self.sample([p[0] + h, p[1], p[2]]) - self.sample([p[0] - h, p[1], p[2]])) * inv,
            (self.sample([p[0], p[1] + h, p[2]]) - self.sample([p[0], p[1] - h, p[2]])) * inv,
            (self.sample([p[0], p[1], p[2] + h]) - self.sample([p[0], p[1], p[2] - h])) * inv,
        ]
    }
}

/// Forwards to the referent, so `&dyn Sdf<Scalar = f32>` and `&ConcreteField`
/// both satisfy the generic `S: Sdf` bound.
///
/// Note that this forwards [`gradient`](Sdf::gradient) as well as
/// [`sample`](Sdf::sample). It must: forwarding only `sample` would silently
/// substitute the central-difference default for every analytic gradient the
/// moment a field crossed a generic boundary by reference — six times the cost,
/// `O(h²)` instead of exact, and nothing anywhere would fail.
impl<S: Sdf + ?Sized> Sdf for &S {
    type Scalar = S::Scalar;

    #[inline]
    fn sample(&self, p: [Self::Scalar; 3]) -> Self::Scalar {
        S::sample(self, p)
    }

    #[inline]
    fn gradient(&self, p: [Self::Scalar; 3]) -> [Self::Scalar; 3] {
        S::gradient(self, p)
    }
}

/// Forwards to the boxed field, including [`gradient`](Sdf::gradient) — see the
/// note on the `&S` implementation for why that matters.
impl<S: Sdf + ?Sized> Sdf for Box<S> {
    type Scalar = S::Scalar;

    #[inline]
    fn sample(&self, p: [Self::Scalar; 3]) -> Self::Scalar {
        S::sample(self, p)
    }

    #[inline]
    fn gradient(&self, p: [Self::Scalar; 3]) -> [Self::Scalar; 3] {
        S::gradient(self, p)
    }
}

/// Sample a field on a grid, one **pre-sliced contiguous write per row**.
///
/// Ticket: R-067 (P-69). The single definition of "evaluate this field at every
/// grid point", and it exists because there were **three** — `marching_cubes`,
/// `marching_tetrahedra` and the dual path each had their own triple loop
/// pushing into a `Vec`, differing only in whether rows carry padding.
///
/// # Why the shape, and not just the sharing
///
/// The three copies all had the same defect. `Vec::push` in the innermost loop
/// re-proves the capacity bound on **every element**, so the store is not a
/// contiguous write LLVM can widen, and the `y`/`z` coordinates were recomputed
/// per element inside it. Here the buffer is sized once, each row is sliced once
/// so its bound is proved once, the two outer coordinates are hoisted, and the
/// inner loop is an `iter_mut().enumerate()` over a slice of known length —
/// which is the shape the autovectorisation literature names (`chunks_exact` and
/// `zip`, pre-slicing outside the loop, struct-of-fields over index arithmetic).
///
/// **Whether it actually vectorises is a property of the field, not of this
/// loop**, and it is decided by `libm`: `sqrtf` carries a `select_implementation`
/// on `target_feature = "sse2"` and on `aarch64`+`neon`, so a field built from
/// `sqrt` reduces to hardware instructions, while `sinf` and `cosf` carry no arch
/// selection at all and are software with argument-reduction branches. A loop
/// containing one cannot be widened at any shape.
///
/// # Bit-identity
///
/// Every value is the same expression it was: `origin[k] + cell_size · k`, with
/// the `y` and `z` terms hoisted rather than re-associated. Hoisting a loop
/// invariant does not change a rounding, so the output is bit-identical to the
/// three loops it replaces — which is `T-007`'s 216 golden hashes, and P-69's C2
/// registers that as a **veto**: a moved hash means LLVM reassociated something
/// and the change is rejected rather than rebaselined.
///
/// `row_stride` is the distance between row starts and must be at least
/// `size[0]`; the excess slots are left at [`Real::ZERO`], which the dual path
/// relies on for its odd-stride padding. `out` is cleared first, so a shorter
/// grid cannot leave a longer one's values behind.
pub(crate) fn sample_grid<R: Real, S: Sdf<Scalar = R> + ?Sized>(
    sdf: &S,
    size: [u32; 3],
    origin: [R; 3],
    cell_size: R,
    row_stride: usize,
    out: &mut alloc::vec::Vec<R>,
) {
    let nx = size[0] as usize;
    debug_assert!(row_stride >= nx, "a row cannot be shorter than the grid");
    let rows = size[1] as usize * size[2] as usize;
    out.clear();
    // Sized once. The pad slots are zero from here and are never written again,
    // which is what lets the inner loop be exactly `nx` wide.
    out.resize(row_stride * rows, R::ZERO);

    for z in 0..size[2] {
        let pz = origin[2] + cell_size * R::from_f64(f64::from(z));
        for y in 0..size[1] {
            let py = origin[1] + cell_size * R::from_f64(f64::from(y));
            let start = row_stride * (y as usize + size[1] as usize * z as usize);
            // One slice, one bound check, `nx` stores.
            let row = &mut out[start..start + nx];
            for (x, slot) in row.iter_mut().enumerate() {
                *slot = sdf.sample([origin[0] + cell_size * R::from_f64(x as f64), py, pz]);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    // The gradient assertions here compare exact values because the fields
    // under test are linear or analytic and the expected result is a specific
    // number, not an approximation of one.
    #![allow(clippy::float_cmp)]

    use super::*;

    /// `f(p) = n·p - d`. Analytic gradient is exactly `n`.
    struct Plane<R: Real> {
        n: [R; 3],
        d: R,
    }

    impl<R: Real> Sdf for Plane<R> {
        type Scalar = R;
        fn sample(&self, p: [R; 3]) -> R {
            self.n[0] * p[0] + self.n[1] * p[1] + self.n[2] * p[2] - self.d
        }
        // Deliberately does NOT override `gradient` -- these tests exercise the
        // default central-difference implementation.
    }

    /// `f(p) = |p| - r`, with an analytic gradient that returns a recognisable
    /// sentinel so forwarding can be distinguished from the default.
    struct Sentinel;

    const SENTINEL: [f32; 3] = [42.0, -17.0, 3.5];

    impl Sdf for Sentinel {
        type Scalar = f32;
        fn sample(&self, p: [f32; 3]) -> f32 {
            p[0]
        }
        fn gradient(&self, _p: [f32; 3]) -> [f32; 3] {
            SENTINEL
        }
    }

    fn norm(v: [f32; 3]) -> f32 {
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }

    /// I-002's stated acceptance criterion: `Sdf` is object-safe, and we intend
    /// `dyn` use.
    #[test]
    fn sdf_is_dyn_compatible() {
        const _: Option<&dyn Sdf<Scalar = f32>> = None;
        let field = Sentinel;
        let dynamic: &dyn Sdf<Scalar = f32> = &field;
        assert_eq!(dynamic.sample([2.0, 0.0, 0.0]), 2.0);
        assert_eq!(dynamic.gradient([0.0; 3]), SENTINEL);

        let boxed: Box<dyn Sdf<Scalar = f32>> = Box::new(Sentinel);
        assert_eq!(boxed.gradient([0.0; 3]), SENTINEL);
    }

    /// A `dyn` field drops into the generic bound every extractor uses, so a
    /// runtime-selected field needs no separate code path.
    #[test]
    fn dyn_sdf_satisfies_the_generic_bound() {
        fn takes_generic<S: Sdf<Scalar = f32>>(s: S) -> f32 {
            s.sample([3.0, 0.0, 0.0])
        }
        let field = Sentinel;
        let dynamic: &dyn Sdf<Scalar = f32> = &field;
        assert_eq!(takes_generic(dynamic), 3.0);
    }

    /// The forwarding trap. If `impl Sdf for &S` forwarded only `sample`, this
    /// would silently return a central difference and every analytic gradient in
    /// `fields` would be discarded whenever a field crossed a generic boundary
    /// by reference.
    #[test]
    fn ref_impl_forwards_analytic_gradient() {
        let field = Sentinel;

        // Named explicitly. `(&field).gradient(..)` would auto-deref straight to
        // the inherent impl and never touch `impl Sdf for &S` at all -- it looks
        // like this test and is not.
        let by_ref: &Sentinel = &field;
        assert_eq!(
            <&Sentinel as Sdf>::gradient(&by_ref, [1.0, 2.0, 3.0]),
            SENTINEL
        );

        let boxed: Box<Sentinel> = Box::new(Sentinel);
        assert_eq!(
            <Box<Sentinel> as Sdf>::gradient(&boxed, [1.0, 2.0, 3.0]),
            SENTINEL
        );

        // And through the generic bound, where `S` really is `&Sentinel`.
        fn through_generic<S: Sdf<Scalar = f32>>(s: S) -> [f32; 3] {
            s.gradient([1.0, 2.0, 3.0])
        }
        assert_eq!(through_generic(&field), SENTINEL);
        assert_eq!(through_generic(Box::new(Sentinel)), SENTINEL);
    }

    /// A plane is linear, so the central difference is exact up to rounding.
    #[test]
    fn default_gradient_is_exact_on_a_plane() {
        let n = [0.267_261_24_f32, 0.534_522_5, 0.801_783_7]; // (1,2,3)/sqrt(14)
        let plane = Plane { n, d: 0.5 };
        let g = plane.gradient([0.25, -0.5, 0.75]);
        for i in 0..3 {
            assert!(
                (g[i] - n[i]).abs() <= 16.0 * f32::EPSILON,
                "axis {i}: {g:?}"
            );
        }
    }

    #[test]
    fn default_gradient_matches_analytic_on_a_sphere() {
        struct Sphere;
        impl Sdf for Sphere {
            type Scalar = f64;
            fn sample(&self, p: [f64; 3]) -> f64 {
                (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() - 1.0
            }
        }
        let p = [0.3, -0.7, 0.5];
        let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        let analytic = [p[0] / len, p[1] / len, p[2] / len];
        let g = Sphere.gradient(p);
        for i in 0..3 {
            assert!((g[i] - analytic[i]).abs() < 1e-9, "axis {i}: {g:?}");
        }
    }

    /// Proves the `max(|p|, 1)` scaling. With a fixed step this returns
    /// `[0, 0, 0]`, because at this magnitude `p + h == p` in `f32`.
    #[test]
    fn default_gradient_survives_large_coordinates() {
        let plane = Plane {
            n: [1.0f32, 0.0, 0.0],
            d: 0.0,
        };
        let p = [1.0e6f32, 0.0, 0.0];

        // The failure this guards against, made explicit.
        assert_eq!(
            p[0] + f32::DIFF_STEP,
            p[0],
            "premise: a fixed step vanishes here"
        );

        let g = plane.gradient(p);
        assert!(g.iter().all(|c| c.is_finite()));
        assert!((g[0] - 1.0).abs() < 1e-3, "{g:?}");
    }

    /// Pins the raw-not-normalised contract: a field scaled by 3 has a gradient
    /// of length 3.
    #[test]
    fn default_gradient_is_not_normalised() {
        let plane = Plane {
            n: [3.0f32, 0.0, 0.0],
            d: 0.0,
        };
        let g = plane.gradient([0.5, 0.5, 0.5]);
        assert!((norm(g) - 3.0).abs() < 1e-3, "{g:?}");
    }
}
