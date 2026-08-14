//! Internal vector helpers over `[R; 3]`.
//!
//! Not a math library and not public API — [`Real`] arithmetic on plain arrays,
//! nothing more. It exists so that the two modules needing a cross product do
//! not each carry their own copy of one.
//!
//! `glam` was once expected to land with the dual-contouring solve; ✗16
//! falsified that. glam has no scalar abstraction — `Mat3` is `f32`, `DMat3` is
//! `f64`, no generic `Mat3<T>` spans them — so it cannot serve a crate generic
//! over [`Real`], and the solve shipped as a six-entry symmetric matrix in
//! `dual_contouring/solve.rs` instead. These helpers *are* the crate's math
//! library. Do not add glam here without reopening ✗16 and CLAUDE.md rule 3.

use crate::Real;

#[inline]
pub(crate) fn sub<R: Real>(a: [R; 3], b: [R; 3]) -> [R; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
pub(crate) fn scale<R: Real>(a: [R; 3], s: R) -> [R; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
pub(crate) fn dot<R: Real>(a: [R; 3], b: [R; 3]) -> R {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub(crate) fn cross<R: Real>(a: [R; 3], b: [R; 3]) -> [R; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
pub(crate) fn length_squared<R: Real>(a: [R; 3]) -> R {
    dot(a, a)
}

#[inline]
pub(crate) fn length<R: Real>(a: [R; 3]) -> R {
    length_squared(a).sqrt()
}

/// Index of the component with the largest magnitude.
///
/// Ties resolve to the lowest index, deterministically — the caller uses this to
/// pick a projection axis, and any consistent choice is correct.
#[inline]
pub(crate) fn dominant_axis<R: Real>(a: [R; 3]) -> usize {
    let m = [a[0].abs(), a[1].abs(), a[2].abs()];
    if m[0] >= m[1] && m[0] >= m[2] {
        0
    } else if m[1] >= m[2] {
        1
    } else {
        2
    }
}
