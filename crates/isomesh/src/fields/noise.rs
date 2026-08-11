//! Deterministic gradient noise, hand-rolled.
//!
//! A noise crate would be a second dependency, which the crate's dependency
//! policy does not allow for a test field. More importantly it would be a
//! dependency whose bit-level output this crate does not control, and
//! `fbm_terrain` feeds committed golden hashes.
//!
//! # Platform determinism
//!
//! This module evaluates **no transcendental function at all**. The only
//! floating-point operations are `+ − * /` and `floor`, every one of which
//! IEEE-754 requires to be correctly rounded, so the result is bit-identical on
//! every conforming platform. The lattice hash is pure integer arithmetic, and
//! the gradients have components of exactly `0` and `±1`, so gradient selection
//! introduces no rounding at all.
//!
//! Rust does not contract `a * b + c` into a fused multiply-add on its own, and
//! nothing here asks for one. Do not introduce `mul_add` in this module: it
//! rounds once where the written expression rounds twice, and the difference is
//! visible in a golden hash.

use crate::Real;

/// Lattice hash.
///
/// `0x9E37_79B1` is 2³²/φ (Knuth's multiplicative constant). `0x85EB_CA6B` and
/// `0xC2B2_AE35` are the finalisation constants from MurmurHash3's `fmix32`
/// (Appleby, 2011). Negative lattice coordinates wrap through two's complement,
/// which Rust defines and which is therefore platform-independent.
#[inline]
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

/// Perlin's twelve improved-noise gradients: the midpoints of a cube's edges.
///
/// Components are exactly `0` or `±1`, so every `g · d` is exact additions and
/// subtractions and no rounding enters gradient selection.
///
/// Twelve edge midpoints rather than the original sixteen-entry table, which has
/// a known directional bias. Selection is `hash % 12`; the modulo bias over 2³²
/// is about 1e-10 and irrelevant here.
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

/// The largest coordinate magnitude for which the `f32` round-trip in
/// [`lattice_index`] is exact: `f32` represents every integer below 2²⁴.
const LATTICE_LIMIT: f64 = 16_777_216.0;

/// Integer lattice coordinate from an already-floored scalar.
///
/// Exact for `|floored| < 2²⁴`, which every reference field stays far inside.
/// Beyond that the `f32` round-trip would start skipping integers, so it is a
/// debug assertion rather than a silent wrong answer.
#[inline]
fn lattice_index<R: Real>(floored: R) -> i32 {
    debug_assert!(
        floored.abs() < R::from_f64(LATTICE_LIMIT),
        "noise coordinate beyond the exactly-representable lattice range"
    );
    floored.as_f32() as i32
}

/// Perlin's quintic fade `u = t³(6t² − 15t + 10)` and its derivative
/// `u' = 30t²(t − 1)²`, which is exact rather than approximated.
#[inline]
fn fade<R: Real>(t: R) -> (R, R) {
    let six = R::from_f64(6.0);
    let fifteen = R::from_f64(15.0);
    let ten = R::from_f64(10.0);
    let thirty = R::from_f64(30.0);

    let t2 = t * t;
    let t3 = t2 * t;
    let u = t3 * (t * (t * six - fifteen) + ten);
    let d = t - R::ONE;
    (u, thirty * t2 * d * d)
}

/// Three-dimensional Perlin gradient noise, with its **exact analytic
/// gradient**.
///
/// Returns `(value, ∇value)`. Both come out of one traversal of the eight cell
/// corners, so the gradient cannot drift away from the value it belongs to.
/// Callers that need only the value get the gradient arithmetic removed by dead
/// code elimination — it is pure arithmetic on locals.
///
/// The value is exactly zero at every integer lattice point, which is a property
/// of the construction rather than a coincidence: at a corner all the weight
/// falls on that corner and its offset vector is zero.
#[inline]
pub(crate) fn perlin<R: Real>(p: [R; 3], seed: u32) -> (R, [R; 3]) {
    let base = [p[0].floor(), p[1].floor(), p[2].floor()];
    let i = [
        lattice_index(base[0]),
        lattice_index(base[1]),
        lattice_index(base[2]),
    ];
    let t = [p[0] - base[0], p[1] - base[1], p[2] - base[2]];

    let (ux, dux) = fade(t[0]);
    let (uy, duy) = fade(t[1]);
    let (uz, duz) = fade(t[2]);

    let mut value = R::ZERO;
    let mut grad = [R::ZERO; 3];

    for cz in 0..2i32 {
        for cy in 0..2i32 {
            for cx in 0..2i32 {
                let g = GRAD12[(hash3(i[0] + cx, i[1] + cy, i[2] + cz, seed) % 12) as usize];
                let gv = [
                    R::from_f64(f64::from(g[0])),
                    R::from_f64(f64::from(g[1])),
                    R::from_f64(f64::from(g[2])),
                ];

                // Offset from this corner to the sample point.
                let d = [
                    t[0] - R::from_f64(f64::from(cx)),
                    t[1] - R::from_f64(f64::from(cy)),
                    t[2] - R::from_f64(f64::from(cz)),
                ];
                let dot = gv[0] * d[0] + gv[1] * d[1] + gv[2] * d[2];

                // Trilinear weight for this corner, and its three derivatives.
                let (wx, dwx) = if cx == 1 {
                    (ux, dux)
                } else {
                    (R::ONE - ux, -dux)
                };
                let (wy, dwy) = if cy == 1 {
                    (uy, duy)
                } else {
                    (R::ONE - uy, -duy)
                };
                let (wz, dwz) = if cz == 1 {
                    (uz, duz)
                } else {
                    (R::ONE - uz, -duz)
                };
                let w = wx * wy * wz;

                value += w * dot;
                // Product rule: the weight varies with position and so does the
                // corner's linear ramp.
                grad[0] += dwx * wy * wz * dot + w * gv[0];
                grad[1] += wx * dwy * wz * dot + w * gv[1];
                grad[2] += wx * wy * dwz * dot + w * gv[2];
            }
        }
    }

    (value, grad)
}

/// A conservative and *provable* bound on `|perlin(..)|`.
///
/// Each gradient has exactly two `±1` components and each corner offset lies in
/// `[-1, 1]³`, so `|g · d| ≤ 2`; the eight weights are non-negative and sum to
/// one, so the interpolant cannot exceed the largest of them.
///
/// The figure usually quoted for three-dimensional Perlin noise is `√3/2 ≈
/// 0.866`. No source in this repository derives it, so it is not used here — and
/// it is worth noting that **this implementation measurably exceeds it**:
/// `noise_stays_within_the_provable_bound` observes a maximum of **0.907716**
/// over 43,200 samples and prints it on every run. Whether that reflects the
/// twelve-gradient variant, the sampling, or the quoted figure simply being
/// wrong is not something this crate needs to settle; it is a reason not to have
/// hard-coded `0.866` and called it a bound.
pub(crate) const PERLIN_BOUND: f64 = 2.0;

/// Per-octave lattice offsets: `octave × (1/φ, 1/φ², 1/φ³)`.
///
/// Without them every octave would vanish at the same points. Perlin noise is
/// exactly zero on its integer lattice, and with `lacunarity = 2` every later
/// octave's lattice contains octave zero's, so the sum would return to the base
/// height on a visible regular grid. The golden ratio's powers are irrational,
/// so no octave's lattice aligns with another's, and the three axes shift by
/// different amounts.
///
/// This does not disturb the analytic gradient: `d/dp n(f·p + k) = f·(∇n)(f·p + k)`
/// for any constant `k`.
const OCTAVE_OFFSET: [f64; 3] = [
    0.618_033_988_749_894_8,  // 1/φ
    0.381_966_011_250_105_15, // 1/φ²
    0.236_067_977_499_789_7,  // 1/φ³
];

/// Fractional Brownian motion over [`perlin`], with its exact analytic gradient.
///
/// ```text
/// fbm(p)  = Σᵢ gainⁱ · n(lacunarityⁱ · frequency · p + offsetᵢ)
/// ∇fbm(p) = Σᵢ gainⁱ · lacunarityⁱ · frequency · ∇n(…)
/// ```
///
/// `lacunarity = 2` is an exact binary scale, so no drift accumulates across
/// octaves.
#[inline]
pub(crate) fn fbm<R: Real>(
    p: [R; 3],
    seed: u32,
    octaves: u32,
    lacunarity: R,
    gain: R,
    frequency: R,
) -> (R, [R; 3]) {
    let mut value = R::ZERO;
    let mut grad = [R::ZERO; 3];
    let mut freq = frequency;
    let mut amp = R::ONE;

    for octave in 0..octaves {
        let k = R::from_f64(f64::from(octave));
        let q = [
            p[0] * freq + k * R::from_f64(OCTAVE_OFFSET[0]),
            p[1] * freq + k * R::from_f64(OCTAVE_OFFSET[1]),
            p[2] * freq + k * R::from_f64(OCTAVE_OFFSET[2]),
        ];
        let (v, g) = perlin(q, seed);

        value += amp * v;
        let chain = amp * freq;
        grad[0] += chain * g[0];
        grad[1] += chain * g[1];
        grad[2] += chain * g[2];

        freq *= lacunarity;
        amp *= gain;
    }

    (value, grad)
}

/// Bound on `|fbm(..)|`: the geometric sum of the per-octave bounds.
#[inline]
pub(crate) fn fbm_bound<R: Real>(octaves: u32, gain: R) -> R {
    let mut total = R::ZERO;
    let mut amp = R::from_f64(PERLIN_BOUND);
    for _ in 0..octaves {
        total += amp;
        amp *= gain;
    }
    total
}

#[cfg(test)]
mod tests {
    // The reduction and bound tests assert exact equality deliberately: one
    // octave at unit gain must be *the same computation* as a single noise call,
    // not merely close to it.
    #![allow(clippy::float_cmp)]

    use super::*;
    use alloc::vec::Vec;

    /// The platform-portability canary. These are committed outputs of pure
    /// integer arithmetic, including negative coordinates, which exercise the
    /// two's-complement wrap. If this ever fails, every `fbm_terrain` golden
    /// hash downstream is invalid.
    #[test]
    fn hash3_matches_committed_values() {
        let cases = [
            (0, 0, 0, 0u32),
            (1, 0, 0, 0u32),
            (0, 1, 0, 0u32),
            (0, 0, 1, 0u32),
            (-1, -1, -1, 0u32),
            (-7, 13, -29, 0x5EED_1234u32),
            (i32::MIN, i32::MAX, 0, 0x5EED_1234u32),
        ];
        let got: Vec<u32> = cases
            .iter()
            .map(|&(x, y, z, s)| hash3(x, y, z, s))
            .collect();
        assert_eq!(
            got,
            [
                0x0000_0000,
                0x11FD_02EB,
                0xCB72_770F,
                0x0C66_C024,
                0x6AB1_7FCE,
                0x738F_6E3F,
                0x43BC_BA29,
            ]
        );
    }

    /// A structural property of the construction, not a coincidence: at a
    /// lattice point one corner takes all the weight and its offset vector is
    /// zero. Cheap, and it catches an off-by-one in the corner loop.
    #[test]
    fn perlin_is_zero_at_lattice_points() {
        for z in -3..=3 {
            for y in -3..=3 {
                for x in -3..=3 {
                    let p = [f64::from(x), f64::from(y), f64::from(z)];
                    let (v, _) = perlin(p, 0x5EED_1234);
                    assert!(v.abs() < 1e-15, "perlin{p:?} = {v}");
                }
            }
        }
    }

    /// The analytic gradient against a central difference, in `f64`. This is the
    /// bug class that matters: a single flipped component barely moves the
    /// magnitude, so direction is checked separately from it.
    #[test]
    fn perlin_gradient_matches_central_difference() {
        let h = 1e-6f64;
        for i in 0..40 {
            let f = f64::from(i);
            // Deliberately off-lattice and off-axis, so no sample lands on a
            // point where the fade derivative is degenerate.
            let p = [0.137 + f * 0.313, -0.241 + f * 0.517, 0.359 - f * 0.211_f64];
            let (_, analytic) = perlin(p, 0x5EED_1234);
            for axis in 0..3 {
                let mut lo = p;
                let mut hi = p;
                lo[axis] -= h;
                hi[axis] += h;
                let numeric = (perlin(hi, 0x5EED_1234).0 - perlin(lo, 0x5EED_1234).0) / (2.0 * h);
                assert!(
                    (analytic[axis] - numeric).abs() < 1e-6,
                    "axis {axis} at {p:?}: analytic {} vs numeric {numeric}",
                    analytic[axis]
                );
            }
        }
    }

    /// Asserts the provable bound and *records* the observed maximum, which is
    /// the number anyone actually wants and which no doc in this repo derives.
    #[test]
    fn noise_stays_within_the_provable_bound() {
        let mut observed = 0.0f64;
        for i in 0..60 {
            for j in 0..60 {
                for k in 0..12 {
                    let p = [
                        -6.0 + f64::from(i) * 0.211,
                        -6.0 + f64::from(j) * 0.211,
                        -1.0 + f64::from(k) * 0.211,
                    ];
                    let (v, _) = perlin(p, 0x5EED_1234);
                    observed = observed.max(v.abs());
                }
            }
        }
        assert!(observed <= PERLIN_BOUND, "observed {observed}");
        std::println!("measured: max |perlin| over 43200 samples = {observed:.6}");
    }

    /// With one octave and unit gain, fBm is a single noise evaluation. Catches
    /// an accumulator initialised wrong.
    #[test]
    fn single_octave_fbm_reduces_to_one_noise_call() {
        let p = [0.31, -0.72, 0.44f64];
        let (v, g) = fbm(p, 0x5EED_1234, 1, 2.0, 1.0, 1.0);
        let (v0, g0) = perlin(p, 0x5EED_1234);
        assert_eq!(v, v0);
        assert_eq!(g, g0);
    }

    #[test]
    fn fbm_gradient_matches_central_difference() {
        let h = 1e-6f64;
        for i in 0..20 {
            let f = f64::from(i);
            let p = [0.137 + f * 0.41, 0.0, 0.359 - f * 0.23f64];
            let (_, analytic) = fbm(p, 0x5EED_1234, 4, 2.0, 0.5, 0.25);
            for axis in [0usize, 2] {
                let mut lo = p;
                let mut hi = p;
                lo[axis] -= h;
                hi[axis] += h;
                let numeric = (fbm(hi, 0x5EED_1234, 4, 2.0, 0.5, 0.25).0
                    - fbm(lo, 0x5EED_1234, 4, 2.0, 0.5, 0.25).0)
                    / (2.0 * h);
                assert!(
                    (analytic[axis] - numeric).abs() < 1e-6,
                    "axis {axis} at {p:?}: {} vs {numeric}",
                    analytic[axis]
                );
            }
        }
    }

    #[test]
    fn fbm_bound_is_the_geometric_sum() {
        // 2 * (1 + 1/2 + 1/4 + 1/8) = 3.75
        assert_eq!(fbm_bound::<f64>(4, 0.5), 3.75);
    }
}
