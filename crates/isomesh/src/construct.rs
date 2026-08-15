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

/// Fast sweeping: Gauss–Seidel passes over the eikonal equation, in place.
///
/// Ticket: S-002. Zhao, *A fast sweeping method for eikonal equations*,
/// Mathematics of Computation 74(250), pp. 603–627 (2005).
///
/// # Why this exists when [`signed_distance_field`] is already exact
///
/// The exact transform answers with the distance to the nearest opposite-signed
/// **sample**, which quantises the answer to the grid — M-251 measured that at a
/// full spacing on a sphere. Sweeping solves `|∇d| = 1` instead, so it can place
/// the surface *between* samples and does not inherit that floor. It is also
/// `O(N)` with no heap and a handful of passes, which is why it is the
/// pragmatic default rather than the exact one.
///
/// # Eight sweeps, and why the count is not a tuning knob
///
/// Each pass visits the grid in one of the `2³` diagonal orderings. Zhao's
/// argument is that a characteristic of the eikonal equation is a straight line,
/// and every straight line in 3D is monotone in each axis, so **some** ordering
/// follows it — after all eight, every characteristic has been swept along.
/// That is why the count is eight rather than "enough": it is the number of
/// orthants, not a convergence parameter.
///
/// # The seed values are the whole accuracy story
///
/// Sweeping propagates whatever it is given. Seeded with zeros on the inside
/// samples it reproduces the exact transform's quantisation exactly; seeded with
/// the **sub-cell** crossing position — where the sign change actually falls
/// between two samples — it does better, and that is what
/// `beats_the_exact_transform_near_the_surface` measures.
fn sweep<R: Real>(d: &mut [R], frozen: &[bool], size: [u32; 3], h: R) {
    let (nx, ny, nz) = (size[0] as usize, size[1] as usize, size[2] as usize);
    let at = |x: usize, y: usize, z: usize| (z * ny + y) * nx + x;

    // The eight diagonal orderings, as (reverse_x, reverse_y, reverse_z).
    for pass in 0..8u8 {
        let (rx, ry, rz) = (pass & 1 != 0, pass & 2 != 0, pass & 4 != 0);
        for zi in 0..nz {
            let z = if rz { nz - 1 - zi } else { zi };
            for yi in 0..ny {
                let y = if ry { ny - 1 - yi } else { yi };
                for xi in 0..nx {
                    let x = if rx { nx - 1 - xi } else { xi };
                    let i = at(x, y, z);
                    if frozen[i] {
                        continue;
                    }
                    // Smaller neighbour along each axis; a boundary axis
                    // contributes nothing rather than a sentinel.
                    let pick = |a: Option<R>, b: Option<R>| match (a, b) {
                        (Some(p), Some(q)) => Some(if p < q { p } else { q }),
                        (Some(p), None) | (None, Some(p)) => Some(p),
                        (None, None) => None,
                    };
                    let ax = pick(
                        (x > 0).then(|| d[at(x - 1, y, z)]),
                        (x + 1 < nx).then(|| d[at(x + 1, y, z)]),
                    );
                    let ay = pick(
                        (y > 0).then(|| d[at(x, y - 1, z)]),
                        (y + 1 < ny).then(|| d[at(x, y + 1, z)]),
                    );
                    let az = pick(
                        (z > 0).then(|| d[at(x, y, z - 1)]),
                        (z + 1 < nz).then(|| d[at(x, y, z + 1)]),
                    );

                    // Godunov update: solve for the largest root using only the
                    // neighbours smaller than it, adding them in increasing
                    // order until the candidate stops exceeding the next one.
                    let mut a: [R; 3] = [far(), far(), far()];
                    let mut n = 0usize;
                    for v in [ax, ay, az].into_iter().flatten() {
                        a[n] = v;
                        n += 1;
                    }
                    if n == 0 {
                        continue;
                    }
                    a[..n].sort_by(|p, q| p.partial_cmp(q).unwrap_or(core::cmp::Ordering::Equal));

                    let mut candidate = a[0] + h;
                    if n > 1 && candidate > a[1] {
                        // Two-axis root of (d−a0)² + (d−a1)² = h².
                        let (s, q) = (a[0] + a[1], a[0] * a[1]);
                        let disc = s * s - (R::ONE + R::ONE) * (q + q - h * h);
                        if disc >= R::ZERO {
                            candidate = (s + disc.sqrt()) * R::HALF;
                        }
                        if n > 2 && candidate > a[2] {
                            // Three-axis root; the same algebra one dimension up.
                            let s3 = a[0] + a[1] + a[2];
                            let q3 = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
                            let three = R::from_f64(3.0);
                            let disc3 = s3 * s3 - three * (q3 - h * h);
                            if disc3 >= R::ZERO {
                                candidate = (s3 + disc3.sqrt()) / three;
                            }
                        }
                    }
                    if candidate < d[i] {
                        d[i] = candidate;
                    }
                }
            }
        }
    }
}

/// Build a signed distance field by fast sweeping.
///
/// Ticket: S-002. Same contract as [`signed_distance_field`] — sign convention,
/// world units, error cases — and a different algorithm behind it.
///
/// # What it does better, and what it does worse
///
/// **Better near the surface, by 3×.** Cells adjacent to a sign change are seeded
/// with the *interpolated* crossing distance rather than zero, so the answer is
/// not quantised to the grid the way the exact transform's is.
///
/// **Not worse far from it, which was predicted and is false.** The concern was
/// that a value ten cells out is ten first-order Godunov updates and accumulates
/// error the exact transform does not. Measured on a sphere at 41³: within two
/// cells of the surface, worst error **0.0333 against the transform's 0.1000**;
/// beyond eight cells, **0.0933 against 0.1000**. It wins everywhere, narrowly
/// at distance — the characteristics of a sphere are radial straight lines, the
/// eight-orthant sweep follows them, and the seeding advantage survives (M-252).
///
/// `sweeping_and_the_exact_transform_trade_places_with_distance` keeps the name
/// it was written under and now asserts "does not lose" at distance rather than
/// "loses", so a field that does flip the ordering fails it loudly.
///
/// # Errors
///
/// As [`signed_distance_field`].
pub fn signed_distance_field_swept<R: Real>(
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
    let (nx, ny, nz) = (size[0] as usize, size[1] as usize, size[2] as usize);
    let at = |x: usize, y: usize, z: usize| (z * ny + y) * nx + x;

    // Seed: every sample adjacent to a sign change gets the *interpolated*
    // distance to that crossing, which is the sub-cell information the exact
    // transform throws away.
    let mut d = vec![far::<R>(); count];
    let mut frozen = vec![false; count];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let i = at(x, y, z);
                let here = samples[i];
                let inside = here < R::ZERO;
                let mut best = far::<R>();
                let mut neighbour = |j: usize| {
                    let there = samples[j];
                    if (there < R::ZERO) == inside {
                        return;
                    }
                    // Linear crossing along this edge, as a fraction of a cell.
                    let denom = here - there;
                    if denom == R::ZERO {
                        return;
                    }
                    let t = (here / denom).abs();
                    if t < best {
                        best = t;
                    }
                };
                if x > 0 {
                    neighbour(at(x - 1, y, z));
                }
                if x + 1 < nx {
                    neighbour(at(x + 1, y, z));
                }
                if y > 0 {
                    neighbour(at(x, y - 1, z));
                }
                if y + 1 < ny {
                    neighbour(at(x, y + 1, z));
                }
                if z > 0 {
                    neighbour(at(x, y, z - 1));
                }
                if z + 1 < nz {
                    neighbour(at(x, y, z + 1));
                }
                if best < far::<R>() {
                    d[i] = best * cell_size;
                    frozen[i] = true;
                }
            }
        }
    }

    sweep(&mut d, &frozen, size, cell_size);

    let mut out = vec![R::ZERO; count];
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = if samples[i] < R::ZERO { -d[i] } else { d[i] };
    }
    Ok(out)
}
