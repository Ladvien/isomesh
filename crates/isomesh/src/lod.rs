//! Coarsening a sampled field, and the four ways people do it.
//!
//! Ticket: T-016. The head-to-head that does not exist in the literature.
//!
//! # There are two families here, and only one of them is in this module
//!
//! **Re-sampling** builds each level by evaluating the field at that spacing.
//! Frisken's adaptive distance fields and Koschier's hp-adaptive work both do
//! this, and it needs no API: the caller already knows how to sample a field, so
//! adding a `Downsample::Reevaluate` variant would be a second path to something
//! that is not downsampling at all.
//!
//! **Downsampling** builds each level from the level below. That is what this
//! module does, and the literature's position is that it is the wrong thing —
//! T-016 exists to measure by how much.
//!
//! # The grids here are sample-centred, which decides the kernels
//!
//! A level of `2ᵏ + 1` samples coarsens to `2ᵏ⁻¹ + 1`, so coarse sample `i` sits
//! exactly on fine sample `2i` and the two levels share their corners. That is
//! the nesting every chunked LOD scheme wants, and it rules out the textbook
//! two-tap box filter: averaging fine samples `2i` and `2i+1` produces a value
//! belonging at `2i + ½`, which is not a coarse sample position. **A half-sample
//! shift per level is a systematic drift, not a smoothing choice**, so every
//! kernel here is odd-length and centred.
//!
//! | operator | kernel | who uses it |
//! |---|---|---|
//! | [`Decimate`](Downsample::Decimate) | `[1]` | naive LOD; the free option |
//! | [`Mean`](Downsample::Mean) | `[1,1,1]/3` | box filter, the classic |
//! | [`Tent`](Downsample::Tent) | `[1,2,1]/4` | linear B-spline; the Haar's successor |
//! | [`Min`](Downsample::Min) | `min` over 3×3×3 | conservative for solids |
//!
//! **The Haar scaling function *is* the box average**, so "wavelet" and "mean"
//! are the same operator in the Haar basis and there is no fifth column to fill.
//! [`Tent`](Downsample::Tent) is the next scaling function up, which is what
//! "use a wavelet instead" actually means in practice.
//!
//! [`Min`](Downsample::Min) is the one that is not a filter. It is what a voxel
//! engine reaches for when inside is negative and losing solid matter is worse
//! than gaining it — the coarse field is then a **conservative under-estimate**
//! of the surface, and nothing thin ever disappears. Whether that is a feature
//! is exactly what the measurement is for.

#[cfg(test)]
mod tests;

use alloc::vec::Vec;

use crate::real::Real;
use crate::shape::{RuntimeShape3, Shape3};

/// How to build a coarse level from a fine one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Downsample {
    /// Take every other sample. No filtering at all.
    Decimate,
    /// Mean over the 3×3×3 neighbourhood — the separable box filter.
    Mean,
    /// Weighted mean with a `[1,2,1]³` kernel — the linear B-spline.
    Tent,
    /// Minimum over the 3×3×3 neighbourhood.
    Min,
}

/// Half-width of every kernel here, in fine samples.
const RADIUS: i64 = 1;

impl Downsample {
    /// Separable weights over `[-RADIUS, RADIUS]`, or `None` for [`Min`], which
    /// is not a linear filter.
    const fn weights(self) -> Option<[f64; 3]> {
        match self {
            Self::Decimate => Some([0.0, 1.0, 0.0]),
            Self::Mean => Some([1.0, 1.0, 1.0]),
            Self::Tent => Some([1.0, 2.0, 1.0]),
            Self::Min => None,
        }
    }

    /// Short name, for tables.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Decimate => "decimate",
            Self::Mean => "mean",
            Self::Tent => "tent",
            Self::Min => "min",
        }
    }

    /// Every operator, for a sweep.
    pub const ALL: [Self; 4] = [Self::Decimate, Self::Mean, Self::Tent, Self::Min];
}

/// Halve a sampled grid.
///
/// `size` must be `2ᵏ + 1` on every axis, with `k ≥ 1`; the result is
/// `2ᵏ⁻¹ + 1`. Returns the coarse values and their shape.
///
/// Samples outside the fine grid are **clamped to the edge**, not treated as
/// zero or as background: a zero would introduce a sign change at every boundary
/// face and mesh a wall around the chunk.
///
/// # Errors
///
/// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if an axis is not
/// `2ᵏ + 1` with `k ≥ 1`, and
/// [`Error::ShapeOverflow`](crate::Error::ShapeOverflow) if `fine` is not one
/// entry per sample.
pub fn downsample<R: Real>(
    fine: &[R],
    shape: &impl Shape3,
    op: Downsample,
) -> crate::Result<(Vec<R>, RuntimeShape3)> {
    let size = shape.size();
    for n in size {
        // `2ᵏ + 1` for `k ≥ 1` means `n - 1` is a power of two and at least 2.
        if n < 3 || !(n - 1).is_power_of_two() {
            return Err(crate::Error::GridTooSmall { size });
        }
    }
    if fine.len() != shape.element_count() {
        return Err(crate::Error::ShapeOverflow {
            size,
            product: fine.len() as u64,
        });
    }

    let coarse_size = [size[0] / 2 + 1, size[1] / 2 + 1, size[2] / 2 + 1];
    let coarse_shape = RuntimeShape3::new(coarse_size)?;
    let (nx, ny, nz) = (size[0] as i64, size[1] as i64, size[2] as i64);
    let at = |x: i64, y: i64, z: i64| -> R {
        let c = |v: i64, n: i64| v.clamp(0, n - 1);
        fine[((c(z, nz) * ny + c(y, ny)) * nx + c(x, nx)) as usize]
    };

    let weights = op.weights();
    let mut out = Vec::with_capacity(coarse_shape.element_count());
    for cz in 0..coarse_size[2] as i64 {
        for cy in 0..coarse_size[1] as i64 {
            for cx in 0..coarse_size[0] as i64 {
                let (fx, fy, fz) = (cx * 2, cy * 2, cz * 2);
                let value = match weights {
                    Some(w) => {
                        let mut sum = R::ZERO;
                        let mut total = R::ZERO;
                        for dz in -RADIUS..=RADIUS {
                            for dy in -RADIUS..=RADIUS {
                                for dx in -RADIUS..=RADIUS {
                                    let k = w[(dx + RADIUS) as usize]
                                        * w[(dy + RADIUS) as usize]
                                        * w[(dz + RADIUS) as usize];
                                    if k == 0.0 {
                                        continue;
                                    }
                                    let k = R::from_f64(k);
                                    sum += at(fx + dx, fy + dy, fz + dz) * k;
                                    total += k;
                                }
                            }
                        }
                        sum / total
                    }
                    None => {
                        let mut best = at(fx, fy, fz);
                        for dz in -RADIUS..=RADIUS {
                            for dy in -RADIUS..=RADIUS {
                                for dx in -RADIUS..=RADIUS {
                                    let v = at(fx + dx, fy + dy, fz + dz);
                                    if v < best {
                                        best = v;
                                    }
                                }
                            }
                        }
                        best
                    }
                };
                out.push(value);
            }
        }
    }

    Ok((out, coarse_shape))
}
