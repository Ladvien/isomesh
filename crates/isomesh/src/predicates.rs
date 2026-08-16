//! Exact geometric predicates: sign tests that are never wrong.
//!
//! A floating-point `orient2d` does not merely lose accuracy near degeneracy —
//! it returns the **wrong sign**, which turns a triangulation's combinatorial
//! decisions into contradictions: a point that is left of `ab` and left of `ba`,
//! a "convex" hull that is not, an edge flip that never terminates. The remedy
//! is not a tolerance. It is to compute the sign exactly.
//!
//! Jonathan Richard Shewchuk, *Adaptive Precision Floating-Point Arithmetic and
//! Fast Robust Geometric Predicates*, Discrete & Computational Geometry 18(3),
//! 1997. [`10.1007/pl00009321`](https://doi.org/10.1007/pl00009321).
//!
//! # This is one path, not a fast path with a fallback
//!
//! [`orient2d`] and [`incircle`] each compute a floating-point estimate, and
//! return it **only when a proven error bound shows its sign cannot be wrong**.
//! Otherwise they compute the determinant exactly. Both branches return the same
//! sign for the same
//! input — the filter is an early exit from one algorithm, not a cheaper
//! substitute for it, and there is no input for which the two disagree. That is
//! the distinction between this and the degraded-fallback pattern the crate
//! forbids: a fallback answers a *different, worse* question when the primary
//! path fails; a certified filter answers the *same* question and declines to
//! answer when it cannot prove it.
//!
//! # What it assumes
//!
//! Radix-2 arithmetic with exact rounding (IEEE 754), and **no fused
//! multiply-add**. If a compiler contracts `a * b - c * d` into an FMA, the
//! roundoff term that `two_product` recovers is no longer the roundoff of the
//! product that was actually computed, and every downstream expansion is
//! silently wrong. Rust does not contract by default, and [`Real`] deliberately
//! omits `mul_add` for a related reason — but this is asserted rather than
//! assumed: see `two_product_is_exact` in the test module, which fails loudly if
//! that ever changes.
//!
//! No allocation: every expansion here has a compile-time length bound, so the
//! whole module is `no_std` with fixed-size arrays and no `alloc`.

use crate::Real;

// ── Shewchuk's primitives (§2.3, §2.5) ──────────────────────────────────────

/// `a ⊕ b` and its exact roundoff, **requiring `|a| ≥ |b|`**.
///
/// Shewchuk Theorem 6. Three operations instead of [`two_sum`]'s six, at the
/// cost of a precondition the caller must guarantee. Where the precondition is
/// merely *likely*, use [`two_sum`]: a violation here is silent, not a panic.
#[inline]
fn fast_two_sum<R: Real>(a: R, b: R) -> (R, R) {
    let x = a + b;
    let b_virtual = x - a;
    (x, b - b_virtual)
}

/// `a ⊕ b` and its exact roundoff, for any `a` and `b`.
///
/// Shewchuk Theorem 7. `x + y = a + b` exactly, with `x` the rounded sum.
#[inline]
fn two_sum<R: Real>(a: R, b: R) -> (R, R) {
    let x = a + b;
    let b_virtual = x - a;
    let a_virtual = x - b_virtual;
    let b_roundoff = b - b_virtual;
    let a_roundoff = a - a_virtual;
    (x, a_roundoff + b_roundoff)
}

/// `a ⊖ b` and its exact roundoff. [`two_sum`] with the sign folded in.
#[inline]
fn two_diff<R: Real>(a: R, b: R) -> (R, R) {
    let x = a - b;
    let b_virtual = a - x;
    let a_virtual = x + b_virtual;
    let b_roundoff = b_virtual - b;
    let a_roundoff = a - a_virtual;
    (x, a_roundoff + b_roundoff)
}

/// Split `a` into two halves of at most `⌈p/2⌉` significand bits each.
///
/// Dekker–Veltkamp, via Shewchuk Theorem 17. The halves are what make the four
/// cross products in [`two_product`] individually exact.
#[inline]
fn split<R: Real>(a: R) -> (R, R) {
    let c = R::SPLITTER * a;
    let a_big = c - a;
    let a_hi = c - a_big;
    (a_hi, a - a_hi)
}

/// `a ⊗ b` and its exact roundoff, so that `x + y = a·b` exactly.
///
/// Shewchuk Theorem 18, which requires `p > 6` — true of `f32` and `f64` by a
/// wide margin.
#[inline]
fn two_product<R: Real>(a: R, b: R) -> (R, R) {
    let x = a * b;
    let (a_hi, a_lo) = split(a);
    let (b_hi, b_lo) = split(b);
    let err1 = x - (a_hi * b_hi);
    let err2 = err1 - (a_lo * b_hi);
    let err3 = err2 - (a_hi * b_lo);
    (x, (a_lo * b_lo) - err3)
}

/// `(a1 + a0) − (b1 + b0)` as a four-component expansion, least significant
/// first.
///
/// Shewchuk's `Two_Two_Diff`, expressed as two nested one-component
/// subtractions. Used to turn a difference of two exact products into an exact
/// expansion.
#[inline]
fn two_two_diff<R: Real>(a1: R, a0: R, b1: R, b0: R) -> [R; 4] {
    // Subtract the less significant half of b, then the more significant half.
    let (i, x0) = two_diff(a0, b0);
    let (j, k) = two_sum(a1, i);
    let (i2, x1) = two_diff(k, b1);
    let (x3, x2) = two_sum(j, i2);
    [x0, x1, x2, x3]
}

// ── Expansion addition (§2.4) ───────────────────────────────────────────────

/// `FAST-EXPANSION-SUM` with zero elimination.
///
/// Adds two nonoverlapping expansions, each in order of increasing magnitude,
/// and writes their exact sum to `h` in the same order. Returns the number of
/// components written; components that come out exactly zero are dropped, so the
/// result is usually far shorter than `e.len() + f.len()`.
///
/// Shewchuk §2.4. The paper's form is: merge `e` and `f` into one sequence `g`
/// of nondecreasing magnitude, then sweep it maintaining an approximate running
/// total `Q` — *"the `Qᵢ` terms maintain an approximate running total"* (Fig. 8)
/// — where each step emits one component and `"Qᵢ and hᵢ₋₁ are produced by a
/// TWO-SUM or FAST-TWO-SUM operation"`. Correctness under round-to-even is
/// Theorem 13; that both inputs are *strongly* nonoverlapping is the hypothesis,
/// and it holds for everything this module feeds it.
///
/// Returns `None` rather than truncating if `h` cannot hold `e.len() + f.len()`
/// components. Every caller here passes a buffer sized at compile time, so the
/// `None` arm is unreachable in this crate — it exists so that a too-small
/// buffer can never silently produce a shortened expansion, which would read as
/// a legitimate answer with the wrong sign.
fn fast_expansion_sum<R: Real>(e: &[R], f: &[R], h: &mut [R]) -> Option<usize> {
    let total = e.len().checked_add(f.len())?;
    if h.len() < total {
        return None;
    }

    // Line 1: merge by nondecreasing magnitude. Both inputs are already sorted
    // that way, so one linear pass suffices and no scratch buffer is needed.
    // Indices rather than iterators so the borrow checker permits the interleave
    // without a self-referential closure; `get` keeps it free of panics.
    let mut ei = 0usize;
    let mut fi = 0usize;
    fn take<R: Real>(e: &[R], ei: &mut usize, f: &[R], fi: &mut usize) -> Option<R> {
        match (e.get(*ei), f.get(*fi)) {
            (Some(&ev), Some(&fv)) => {
                if ev.abs() <= fv.abs() {
                    *ei += 1;
                    Some(ev)
                } else {
                    *fi += 1;
                    Some(fv)
                }
            }
            (Some(&ev), None) => {
                *ei += 1;
                Some(ev)
            }
            (None, Some(&fv)) => {
                *fi += 1;
                Some(fv)
            }
            (None, None) => None,
        }
    }

    let Some(g0) = take(e, &mut ei, f, &mut fi) else {
        return Some(0);
    };
    let Some(g1) = take(e, &mut ei, f, &mut fi) else {
        // A single component is already its own expansion. Zero-eliminate it.
        return Some(if g0 == R::ZERO {
            0
        } else {
            *h.first_mut()? = g0;
            1
        });
    };

    // Line 2: the first pair. FAST-TWO-SUM is justified because g1 is the later
    // element of a nondecreasing-magnitude merge, so |g1| >= |g0|.
    let (mut q, first) = fast_two_sum(g1, g0);
    let mut written = 0usize;
    let emit = |value: R, h: &mut [R], written: &mut usize| {
        if value != R::ZERO {
            match h.get_mut(*written) {
                Some(slot) => {
                    *slot = value;
                    *written += 1;
                }
                None => return false,
            }
        }
        true
    };
    if !emit(first, h, &mut written) {
        return None;
    }

    // Lines 3-4: sweep the rest, emitting one roundoff term per step.
    while let Some(g) = take(e, &mut ei, f, &mut fi) {
        let (q_new, roundoff) = two_sum(q, g);
        q = q_new;
        if !emit(roundoff, h, &mut written) {
            return None;
        }
    }

    // Line 5: the running total is the most significant component. It is kept
    // even when zero if nothing else survived, so that the empty expansion is
    // represented as a single zero rather than as nothing at all.
    if q != R::ZERO || written == 0 {
        *h.get_mut(written)? = q;
        written += 1;
    }

    Some(written)
}

/// `SCALE-EXPANSION` with zero elimination: the exact product of an expansion
/// and a single scalar.
///
/// Shewchuk §2.6, Theorem 19. The paper's line structure is fixed by its own
/// proof, which names the invariant
/// `Q₂ᵢ + Σⱼ₌₁^(2i−1) hⱼ = Σⱼ₌₁^i eⱼb`, states that it *"holds for i = 1 after
/// Line 1 is executed"*, inducts *"on Lines 3, 4, and 5"*, notes *"the use of
/// FAST-TWO-SUM in Line 5"*, and concludes that *"after Line 6 is executed,
/// Σⱼ₌₁^(2m) hⱼ = b Σⱼ₌₁^m eⱼ"*.
///
/// `h` must have room for `2 · e.len()`; returns `None` rather than truncating,
/// for the reason [`fast_expansion_sum`] does.
fn scale_expansion<R: Real>(e: &[R], b: R, h: &mut [R]) -> Option<usize> {
    if h.len() < e.len().checked_mul(2)? {
        return None;
    }
    let mut written = 0usize;
    let push = |value: R, h: &mut [R], written: &mut usize| -> bool {
        if value == R::ZERO {
            return true;
        }
        match h.get_mut(*written) {
            Some(slot) => {
                *slot = value;
                *written += 1;
                true
            }
            None => false,
        }
    };

    let Some((&first, rest)) = e.split_first() else {
        return Some(0);
    };
    // Line 1.
    let (mut q, h1) = two_product(first, b);
    if !push(h1, h, &mut written) {
        return None;
    }
    // Lines 2-5.
    for &e_i in rest {
        let (t_hi, t_lo) = two_product(e_i, b);
        let (q_sum, h_even) = two_sum(q, t_lo);
        if !push(h_even, h, &mut written) {
            return None;
        }
        // FAST-TWO-SUM is sound here because `t_hi` is a product's high word and
        // `q_sum` is the accumulated remainder, so |t_hi| >= |q_sum|.
        let (h_odd, q_next) = fast_two_sum(t_hi, q_sum);
        if !push(h_odd, h, &mut written) {
            return None;
        }
        q = q_next;
    }
    // Line 6.
    if (q != R::ZERO || written == 0) && !push_final(q, h, &mut written) {
        return None;
    }
    Some(written)
}

/// The `q != 0 || nothing written` tail shared by the two expansion routines.
fn push_final<R: Real>(value: R, h: &mut [R], written: &mut usize) -> bool {
    match h.get_mut(*written) {
        Some(slot) => {
            *slot = value;
            *written += 1;
            true
        }
        None => false,
    }
}

// ── ORIENT2D (§4.3) ─────────────────────────────────────────────────────────

/// Error bound for the unfiltered floating-point estimate: `(3ε + 16ε²)`.
///
/// Shewchuk Table 1, stage A. The paper notes this *"should be exactly computed
/// once at program initialization"*; here it is a `const fn` of [`Real`]'s
/// [`UNIT_ROUNDOFF`](Real::UNIT_ROUNDOFF), which the compiler folds.
#[inline]
fn ccw_error_bound_a<R: Real>() -> R {
    let e = R::UNIT_ROUNDOFF;
    (R::from_f64(3.0) + R::from_f64(16.0) * e) * e
}

/// Twice the signed area of triangle `(a, b, c)`; **the sign is always correct**.
///
/// Positive when `a`, `b`, `c` are counterclockwise, negative when clockwise, and
/// exactly zero **only** when the three points are exactly collinear. The
/// magnitude is a good approximation to the determinant but is not itself exact;
/// only the sign carries a guarantee, which is what every combinatorial use
/// needs.
///
/// Equivalently: positive when `c` lies to the left of the directed line `ab`.
///
/// # Coordinate order
///
/// Each point is `[x, y]`. The determinant computed is
/// `(aₓ−cₓ)(b_y−c_y) − (a_y−c_y)(bₓ−cₓ)`.
///
/// # Cost
///
/// One multiply-and-subtract plus a comparison in the overwhelmingly common
/// case. The exact path runs only when the estimate's proven error bound does
/// not separate it from zero, which for well-separated input essentially never
/// happens. See the module docs on why this is one algorithm rather than two.
#[must_use]
pub fn orient2d<R: Real>(a: [R; 2], b: [R; 2], c: [R; 2]) -> R {
    let det_left = (a[0] - c[0]) * (b[1] - c[1]);
    let det_right = (a[1] - c[1]) * (b[0] - c[0]);
    let det = det_left - det_right;

    // When the two products have strictly opposite signs -- or either is zero --
    // the subtraction cannot cancel, so the estimate's sign is already the
    // determinant's. Shewchuk notes this is why the bound is zero "in the common
    // case that all three input points lie on a horizontal or vertical line".
    let det_sum = if det_left > R::ZERO {
        if det_right <= R::ZERO {
            return det;
        }
        det_left + det_right
    } else if det_left < R::ZERO {
        if det_right >= R::ZERO {
            return det;
        }
        -det_left - det_right
    } else {
        return det;
    };

    let error_bound = ccw_error_bound_a::<R>() * det_sum;
    if det >= error_bound || -det >= error_bound {
        return det;
    }

    orient2d_exact(a, b, c)
}

/// The determinant as an exact expansion, reduced to its most significant
/// component.
///
/// Shewchuk's `orient2dexact`. Expands
/// `aₓb_y − aₓc_y + bₓc_y − bₓa_y + cₓa_y − cₓb_y` into three four-component
/// terms and sums them. The most significant component of a nonoverlapping
/// expansion carries its sign, which is the only thing promised.
fn orient2d_exact<R: Real>(a: [R; 2], b: [R; 2], c: [R; 2]) -> R {
    let (axby1, axby0) = two_product(a[0], b[1]);
    let (axcy1, axcy0) = two_product(a[0], c[1]);
    let a_term = two_two_diff(axby1, axby0, axcy1, axcy0);

    let (bxcy1, bxcy0) = two_product(b[0], c[1]);
    let (bxay1, bxay0) = two_product(b[0], a[1]);
    let b_term = two_two_diff(bxcy1, bxcy0, bxay1, bxay0);

    let (cxay1, cxay0) = two_product(c[0], a[1]);
    let (cxby1, cxby0) = two_product(c[0], b[1]);
    let c_term = two_two_diff(cxay1, cxay0, cxby1, cxby0);

    let mut ab = [R::ZERO; 8];
    let mut total = [R::ZERO; 12];

    // Both sums are sized exactly at their inputs' worst case, so neither can
    // report `None`. Treating that arm as a zero determinant would be a silent
    // wrong answer, so it is reached only if this module is edited wrongly, and
    // the tests cover the arities.
    let Some(ab_len) = fast_expansion_sum(&a_term, &b_term, &mut ab) else {
        return R::ZERO;
    };
    let Some(total_len) = fast_expansion_sum(ab.get(..ab_len).unwrap_or(&[]), &c_term, &mut total)
    else {
        return R::ZERO;
    };

    // The most significant component carries the sign. `total_len` is at least
    // one: `fast_expansion_sum` keeps the running total even when it is zero.
    total_len
        .checked_sub(1)
        .and_then(|last| total.get(last).copied())
        .unwrap_or(R::ZERO)
}

// ── INCIRCLE (§4.4) ─────────────────────────────────────────────────────────

/// Error bound for the unfiltered floating-point estimate: `(10ε + 96ε²)`.
///
/// Shewchuk Table 5, stage A — *"Error bounds for the expansions calculated by
/// INCIRCLE"*. Note it is **not** [`ccw_error_bound_a`]'s constant: the incircle
/// determinant is 3×3 with squared entries, so it accumulates more roundoff.
#[inline]
fn incircle_error_bound_a<R: Real>() -> R {
    let e = R::UNIT_ROUNDOFF;
    (R::from_f64(10.0) + R::from_f64(96.0) * e) * e
}

/// `u ⨯ v` as an exact four-component expansion: `uₓv_y − u_yvₓ`.
#[inline]
fn cross_exact<R: Real>(u: [R; 2], v: [R; 2]) -> [R; 4] {
    let (p1, p0) = two_product(u[0], v[1]);
    let (q1, q0) = two_product(u[1], v[0]);
    two_two_diff(p1, p0, q1, q0)
}

/// Negate every component. Exact: negation of a float is exact, and it preserves
/// both the nonoverlapping property and the magnitude ordering.
#[inline]
fn negated<R: Real>(e: [R; 4]) -> [R; 4] {
    [-e[0], -e[1], -e[2], -e[3]]
}

/// Whether `d` lies inside the circle through `a`, `b` and `c`; **the sign is
/// always correct**.
///
/// Positive when `d` is strictly inside, negative when strictly outside, and
/// exactly zero **only** when the four points are exactly cocircular. As with
/// [`orient2d`], only the sign carries a guarantee.
///
/// # Orientation matters
///
/// The sign is stated for `a`, `b`, `c` in **counterclockwise** order — that is,
/// `orient2d(a, b, c) > 0`. If they are clockwise the sense inverts, because the
/// determinant is antisymmetric under swapping two points. Callers that cannot
/// guarantee the winding should test it with [`orient2d`] first; this function
/// does not, because doing so silently would cost every caller a predicate they
/// usually already know the answer to.
///
/// # Coordinate order
///
/// Each point is `[x, y]`. The determinant is the standard lifted one,
/// `|aₓ a_y aₓ²+a_y² 1; …|`, expanded by cofactors so that no coordinate
/// difference is ever formed — differences round, and an exact predicate cannot
/// afford that.
///
/// # Cost
///
/// The filtered estimate plus a comparison in the common case; the exact
/// expansion only where the proven bound does not separate the estimate from
/// zero. See the module docs on why this is one algorithm rather than two.
#[must_use]
pub fn incircle<R: Real>(a: [R; 2], b: [R; 2], c: [R; 2], d: [R; 2]) -> R {
    let adx = a[0] - d[0];
    let ady = a[1] - d[1];
    let bdx = b[0] - d[0];
    let bdy = b[1] - d[1];
    let cdx = c[0] - d[0];
    let cdy = c[1] - d[1];

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let alift = adx * adx + ady * ady;

    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let blift = bdx * bdx + bdy * bdy;

    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;
    let clift = cdx * cdx + cdy * cdy;

    let det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);

    // The permanent: the same expression with every subtraction made an addition
    // and every factor its magnitude, which bounds the roundoff the subtractions
    // can hide. This is orient2d's `det_sum` construction one dimension up.
    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * alift
        + (cdxady.abs() + adxcdy.abs()) * blift
        + (adxbdy.abs() + bdxady.abs()) * clift;

    let error_bound = incircle_error_bound_a::<R>() * permanent;
    if det > error_bound || -det > error_bound {
        return det;
    }

    incircle_exact(a, b, c, d)
}

/// The lifted determinant as an exact expansion, reduced to its most significant
/// component.
///
/// Expanded by cofactors along the `x² + y²` column, which gives
/// `z_a·M_a − z_b·M_b + z_c·M_c − z_d·M_d` where each `M` is a signed sum of
/// three 2×2 cross products and `z = x² + y²`. Every `z·M` is formed as
/// `x·(x·M) + y·(y·M)` so that only expansion-times-scalar is ever needed —
/// [`scale_expansion`] — rather than expansion-times-expansion.
///
/// **Nothing here forms a coordinate difference.** The filtered estimate above
/// subtracts `d` from each point first, which is what makes it cheap and what
/// makes it inexact; the exact path cannot, because those differences round.
fn incircle_exact<R: Real>(a: [R; 2], b: [R; 2], c: [R; 2], d: [R; 2]) -> R {
    let ab = cross_exact(a, b);
    let ac = cross_exact(a, c);
    let ad = cross_exact(a, d);
    let bc = cross_exact(b, c);
    let bd = cross_exact(b, d);
    let cd = cross_exact(c, d);

    // Buffer sizes are the exact worst case, derived rather than guessed, and
    // `expansion_buffers_are_exactly_the_worst_case` pins the chain. Every
    // `None` arm below is therefore unreachable -- which matters, because the
    // only thing this function could return on overflow is `ZERO`, and `ZERO`
    // means "exactly cocircular". A silently shortened expansion here would not
    // look like a failure, it would look like a degenerate configuration.
    //
    //   cross_exact                     -> 4          (two_two_diff is always 4)
    //   sum3      = 4+4 -> 8, 8+4       -> 12
    //   lift      = scale(12) -> 24, scale(24) -> 48, per coordinate
    //               48 + 48             -> 96
    //   total     = 96, 192, 288        -> 384        (four terms accumulated)
    //
    // `fast_expansion_sum` needs `e.len() + f.len()`, and the final accumulation
    // is 288 + 96 = 384, so `scratch` is exactly filled at the worst case.
    let mut minor = [R::ZERO; 12];
    let mut term = [R::ZERO; 96];
    let mut total = [R::ZERO; 384];
    let mut total_len = 0usize;
    let mut scratch = [R::ZERO; 384];

    // (point, its minor's three signed cross products, the term's sign)
    let contributions: [([R; 2], [[R; 4]; 3], bool); 4] = [
        (a, [bc, negated(bd), cd], false),
        (b, [ac, negated(ad), cd], true),
        (c, [ab, negated(ad), bd], false),
        (d, [ab, negated(ac), bc], true),
    ];

    for (point, parts, negate_term) in contributions {
        let Some(minor_len) = sum3(&parts[0], &parts[1], &parts[2], &mut minor) else {
            return R::ZERO;
        };
        let Some(term_len) = lift(point, minor.get(..minor_len).unwrap_or(&[]), &mut term) else {
            return R::ZERO;
        };
        if negate_term {
            for slot in term.iter_mut().take(term_len) {
                *slot = -*slot;
            }
        }
        let Some(next) = fast_expansion_sum(
            total.get(..total_len).unwrap_or(&[]),
            term.get(..term_len).unwrap_or(&[]),
            &mut scratch,
        ) else {
            return R::ZERO;
        };
        match total.get_mut(..next) {
            Some(dst) => match scratch.get(..next) {
                Some(src) => dst.copy_from_slice(src),
                None => return R::ZERO,
            },
            None => return R::ZERO,
        }
        total_len = next;
    }

    total_len
        .checked_sub(1)
        .and_then(|last| total.get(last).copied())
        .unwrap_or(R::ZERO)
}

/// `e + f + g` for three four-component expansions, into a buffer of 12.
fn sum3<R: Real>(e: &[R; 4], f: &[R; 4], g: &[R; 4], out: &mut [R; 12]) -> Option<usize> {
    let mut ef = [R::ZERO; 8];
    let ef_len = fast_expansion_sum(e, f, &mut ef)?;
    fast_expansion_sum(ef.get(..ef_len).unwrap_or(&[]), g, out)
}

/// `(x² + y²) · minor`, as `x·(x·minor) + y·(y·minor)`.
///
/// Squaring by scaling twice keeps every multiplication expansion-by-scalar. The
/// intermediate lengths are `12 → 24 → 48` per coordinate, so the sum fits 96.
fn lift<R: Real>(point: [R; 2], minor: &[R], out: &mut [R; 96]) -> Option<usize> {
    let mut once = [R::ZERO; 24];
    let mut x_squared = [R::ZERO; 48];
    let mut y_squared = [R::ZERO; 48];

    let once_len = scale_expansion(minor, point[0], &mut once)?;
    let x_len = scale_expansion(once.get(..once_len)?, point[0], &mut x_squared)?;

    let once_len = scale_expansion(minor, point[1], &mut once)?;
    let y_len = scale_expansion(once.get(..once_len)?, point[1], &mut y_squared)?;

    fast_expansion_sum(x_squared.get(..x_len)?, y_squared.get(..y_len)?, out)
}

#[cfg(test)]
mod tests;
