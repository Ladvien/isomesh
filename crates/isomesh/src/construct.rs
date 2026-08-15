//! Building a distance field, rather than consuming one.
//!
//! Ticket: S-001, the first of Phase 13. Everything else in this crate takes an
//! [`Sdf`](crate::Sdf) and produces triangles; this takes samples and produces a
//! distance.
//!
//! # Why the exact transform comes first
//!
//! It is the **ground truth** every other constructor is measured against —
//! fast sweeping, fast marching and jump flooding are all approximations of this
//! — and it is the cheapest of them to get exactly right, because it is exact by
//! construction rather than by convergence.
//!
//! Felzenszwalb & Huttenlocher, *Distance Transforms of Sampled Functions*,
//! Theory of Computing 8(19), pp. 415–428 (2012). The algorithm is separable and
//! `O(n)` per dimension: the 3D transform is the 1D transform run along each
//! axis in turn, which is what makes an exact answer affordable at all.
//!
//! # The 1D step is a lower envelope of parabolas
//!
//! For a sampled function `f`, the transform is
//!
//! ```text
//! D(q) = min over p of ( (q − p)² + f(p) )
//! ```
//!
//! Each `p` contributes a parabola of the same shape rooted at `(p, f(p))`, so
//! `D` is their lower envelope. Two parabolas cross at exactly one point, which
//! is what makes the envelope computable in one forward pass: maintain the
//! parabolas currently on it and where each takes over, and a new parabola can
//! only ever displace them from the right.
//!
//! **Squared distances throughout, and the square root taken once at the end.**
//! The recurrence is exact in squared space and would accumulate error in
//! rooted space, and the intersection formula below divides by `2(q − p)`, which
//! only has that form because the terms are squares.

#[cfg(test)]
mod tests;

use alloc::vec;
use alloc::vec::Vec;

use crate::real::Real;
use crate::shape::Shape3;

/// A value standing in for "no seed anywhere near here".
///
/// Not infinity: the recurrence adds `(q − p)²` to it, and `∞ + x` is `∞` but
/// `∞ − ∞` in the intersection formula is a NaN that propagates silently through
/// the envelope. A large finite number keeps every comparison meaningful.
fn far<R: Real>() -> R {
    R::from_f64(1e30)
}

/// One pass of the exact squared distance transform along a single axis.
///
/// `f` is the input row and `d` receives the transform. `v` and `z` are scratch
/// of length `n` and `n + 1` — passed in rather than allocated so a 3D transform
/// allocates once rather than once per row.
///
/// This is Felzenszwalb & Huttenlocher's Figure 1, with their variable names
/// kept so the source and the paper can be read side by side: `k` indexes the
/// rightmost parabola on the envelope, `v[k]` is which sample it belongs to, and
/// `z[k]` is where it takes over from its predecessor.
fn transform_row<R: Real>(f: &[R], d: &mut [R], v: &mut [usize], z: &mut [R]) {
    let n = f.len();
    if n == 0 {
        return;
    }
    let mut k = 0usize;
    v[0] = 0;
    z[0] = -far::<R>();
    z[1] = far::<R>();

    for q in 1..n {
        let fq = f[q];
        let qr = R::from_f64(q as f64);
        loop {
            let p = v[k];
            let pr = R::from_f64(p as f64);
            // Where the parabola from `q` meets the one from `p`. The `2(q − p)`
            // denominator is never zero because `q > p` holds by construction.
            let s = ((fq + qr * qr) - (f[p] + pr * pr)) / ((qr - pr) + (qr - pr));
            if s > z[k] {
                k += 1;
                v[k] = q;
                z[k] = s;
                z[k + 1] = far::<R>();
                break;
            }
            // `q`'s parabola hides the last one on the envelope; drop it and
            // reconsider against the one before.
            if k == 0 {
                v[0] = q;
                z[0] = -far::<R>();
                z[1] = far::<R>();
                break;
            }
            k -= 1;
        }
    }

    let mut k = 0usize;
    for (q, slot) in d.iter_mut().enumerate() {
        let qr = R::from_f64(q as f64);
        while z[k + 1] < qr {
            k += 1;
        }
        let p = v[k];
        let pr = R::from_f64(p as f64);
        *slot = (qr - pr) * (qr - pr) + f[p];
    }
}

/// The exact squared Euclidean distance transform of a 3D grid, in place.
///
/// `grid` holds the seed cost per sample — zero where the feature is, [`far`]
/// where it is not — in `x`-fastest order, matching [`Shape3`].
fn squared_edt<R: Real>(grid: &mut [R], size: [u32; 3]) {
    let (nx, ny, nz) = (size[0] as usize, size[1] as usize, size[2] as usize);
    let longest = nx.max(ny).max(nz);
    let mut f = vec![R::ZERO; longest];
    let mut d = vec![R::ZERO; longest];
    let mut v = vec![0usize; longest];
    let mut z = vec![R::ZERO; longest + 1];

    let at = |x: usize, y: usize, z: usize| (z * ny + y) * nx + x;

    for zi in 0..nz {
        for yi in 0..ny {
            for xi in 0..nx {
                f[xi] = grid[at(xi, yi, zi)];
            }
            transform_row(&f[..nx], &mut d[..nx], &mut v[..nx], &mut z[..nx + 1]);
            for xi in 0..nx {
                grid[at(xi, yi, zi)] = d[xi];
            }
        }
    }
    for zi in 0..nz {
        for xi in 0..nx {
            for yi in 0..ny {
                f[yi] = grid[at(xi, yi, zi)];
            }
            transform_row(&f[..ny], &mut d[..ny], &mut v[..ny], &mut z[..ny + 1]);
            for yi in 0..ny {
                grid[at(xi, yi, zi)] = d[yi];
            }
        }
    }
    for yi in 0..ny {
        for xi in 0..nx {
            for zi in 0..nz {
                f[zi] = grid[at(xi, yi, zi)];
            }
            transform_row(&f[..nz], &mut d[..nz], &mut v[..nz], &mut z[..nz + 1]);
            for zi in 0..nz {
                grid[at(xi, yi, zi)] = d[zi];
            }
        }
    }
}

/// Build a signed distance field from sampled values.
///
/// `samples` are read only for their **sign**: negative is inside, matching the
/// crate's convention everywhere else. The result is the Euclidean distance to
/// the nearest sample of opposite sign, negated inside, in world units.
///
/// # Two transforms, not one
///
/// A single transform gives the distance to the inside, which is zero
/// everywhere inside and useless there. So the inside and the outside are
/// transformed separately and subtracted: `d_out − d_in`. That is the standard
/// construction and it is why the result is signed rather than merely a
/// distance.
///
/// # What the answer is a distance *to*
///
/// The nearest **sample of opposite sign**, not the nearest point of the true
/// surface. Those differ by up to half a cell, because the surface passes
/// between two samples and this cannot see where. That is the resolution limit
/// of a sampled transform rather than an error in it, and
/// `matches_the_analytic_sphere_within_one_spacing` pins the size of the gap.
///
/// # Errors
///
/// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if any axis has fewer
/// than two samples, and
/// [`Error::ShapeOverflow`](crate::Error::ShapeOverflow) if `samples` does not
/// match `shape`.
pub fn signed_distance_field<R: Real>(
    samples: &[R],
    shape: &impl Shape3,
    cell_size: R,
) -> crate::Result<Vec<R>> {
    let size = shape.size();
    if size[0] < 2 || size[1] < 2 || size[2] < 2 {
        return Err(crate::Error::GridTooSmall { size });
    }
    let count = shape.element_count();
    if samples.len() != count {
        return Err(crate::Error::ShapeOverflow {
            size,
            product: samples.len() as u64,
        });
    }

    // Seeds are the samples on the *other* side, so each transform measures the
    // distance to the boundary rather than to itself.
    let mut inside = vec![far::<R>(); count];
    let mut outside = vec![far::<R>(); count];
    for (i, &s) in samples.iter().enumerate() {
        if s < R::ZERO {
            inside[i] = R::ZERO;
        } else {
            outside[i] = R::ZERO;
        }
    }
    squared_edt(&mut inside, size);
    squared_edt(&mut outside, size);

    let mut out = vec![R::ZERO; count];
    for (i, slot) in out.iter_mut().enumerate() {
        // Distance to the nearest opposite-signed sample, in samples, then
        // scaled to world units once.
        let d = if samples[i] < R::ZERO {
            -outside[i].sqrt()
        } else {
            inside[i].sqrt()
        };
        *slot = d * cell_size;
    }
    Ok(out)
}
