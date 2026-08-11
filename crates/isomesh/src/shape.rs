//! Linearisation of a 3D index space.

/// Maps a 3D grid coordinate to a flat buffer index and back.
///
/// Deliberately array-based: `[u32; 3]` in, `u32` out. No math library appears
/// here, for the same reason it appears nowhere else in the public API.
///
/// # Index order
///
/// **`x` varies fastest.** For a shape of size `[sx, sy, sz]`:
///
/// ```text
/// linearize([x, y, z]) == x + y * sx + z * sx * sy
/// strides              == [1, sx, sx * sy]
/// ```
///
/// So a scan over increasing `x` is contiguous in memory. Every extraction inner
/// loop in this crate walks `x` innermost and relies on that. Reversing the
/// convention does not produce wrong meshes — it produces slow ones, which is
/// considerably harder to notice.
///
/// This is stated as strides rather than as "row-major", which is ambiguous in
/// three dimensions and would otherwise be re-litigated in every review.
///
/// # Bounds
///
/// [`linearize`](Shape3::linearize) does **not** bounds-check in release builds;
/// it is on the hottest path in the crate. It `debug_assert!`s instead.
/// Out-of-range input in release produces a meaningless index, not a panic.
pub trait Shape3 {
    /// `[sx, sy, sz]`.
    fn size(&self) -> [u32; 3];

    /// Flat index for a grid coordinate. See the trait docs for the formula.
    fn linearize(&self, p: [u32; 3]) -> u32;

    /// Inverse of [`linearize`](Shape3::linearize).
    fn delinearize(&self, i: u32) -> [u32; 3];

    /// Total number of elements, `sx * sy * sz`.
    ///
    /// `usize` because it is what you hand to `Vec::with_capacity`. Named
    /// `element_count` rather than `len` because this type also has a `size`,
    /// and "length" of a 3D shape is not a thing anyone should have to guess at.
    #[inline]
    fn element_count(&self) -> usize {
        let [x, y, z] = self.size();
        x as usize * y as usize * z as usize
    }
}

/// A shape whose dimensions are known at run time.
///
/// Required by anything that sweeps resolution — the benchmark harness meshes
/// the same field from 16³ to 256³, so the resolution cannot live in the type.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct RuntimeShape3 {
    size: [u32; 3],
}

impl RuntimeShape3 {
    /// # Errors
    ///
    /// [`Error::ShapeOverflow`] if `sx * sy * sz` does not fit in `u32`. This is
    /// reported rather than wrapped because a silent wrap aliases distinct cells
    /// onto the same index and produces a mesh no validity test can catch — the
    /// topology comes out self-consistent, it is just not the field's.
    pub fn new(size: [u32; 3]) -> crate::Result<Self> {
        let fits = size[0]
            .checked_mul(size[1])
            .and_then(|xy| xy.checked_mul(size[2]))
            .is_some();
        if fits {
            Ok(Self { size })
        } else {
            Err(crate::Error::ShapeOverflow {
                size,
                product: u64::from(size[0]) * u64::from(size[1]) * u64::from(size[2]),
            })
        }
    }
}

impl Shape3 for RuntimeShape3 {
    #[inline]
    fn size(&self) -> [u32; 3] {
        self.size
    }

    #[inline]
    fn linearize(&self, p: [u32; 3]) -> u32 {
        debug_assert!(
            p[0] < self.size[0] && p[1] < self.size[1] && p[2] < self.size[2],
            "coordinate {p:?} out of range for shape {:?}",
            self.size
        );
        p[0] + self.size[0] * (p[1] + self.size[1] * p[2])
    }

    #[inline]
    fn delinearize(&self, i: u32) -> [u32; 3] {
        debug_assert!(
            (i as usize) < self.element_count(),
            "index {i} out of range for shape {:?}",
            self.size
        );
        let [sx, sy, _] = self.size;
        [i % sx, (i / sx) % sy, i / (sx * sy)]
    }
}

/// A shape whose dimensions are known at compile time.
///
/// Turns [`linearize`](Shape3::linearize) into multiplication by literals, which
/// is the entire reason this exists alongside [`RuntimeShape3`]. Chunk meshing
/// knows its chunk size statically; the benchmark sweep does not.
///
/// These are two implementations of one trait, not two paths for one feature:
/// a caller picks by whether its dimensions are static, and both compute
/// identical indices — asserted in the tests.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ConstShape3<const X: u32, const Y: u32, const Z: u32>;

impl<const X: u32, const Y: u32, const Z: u32> ConstShape3<X, Y, Z> {
    /// # Compile errors
    ///
    /// Fails to compile if `X * Y * Z` exceeds `u32::MAX`. An inline `const`
    /// block inherits the enclosing item's generic parameters, so the check
    /// fires at the instantiation site — the loudest failure available.
    #[must_use]
    pub const fn new() -> Self {
        const { assert!(X as u64 * Y as u64 * Z as u64 <= u32::MAX as u64) }
        Self
    }
}

impl<const X: u32, const Y: u32, const Z: u32> Default for ConstShape3<X, Y, Z> {
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl<const X: u32, const Y: u32, const Z: u32> Shape3 for ConstShape3<X, Y, Z> {
    #[inline]
    fn size(&self) -> [u32; 3] {
        [X, Y, Z]
    }

    #[inline]
    fn linearize(&self, p: [u32; 3]) -> u32 {
        debug_assert!(
            p[0] < X && p[1] < Y && p[2] < Z,
            "coordinate {p:?} out of range for shape [{X}, {Y}, {Z}]"
        );
        p[0] + X * (p[1] + Y * p[2])
    }

    #[inline]
    fn delinearize(&self, i: u32) -> [u32; 3] {
        debug_assert!(
            i < X * Y * Z,
            "index {i} out of range for shape [{X}, {Y}, {Z}]"
        );
        [i % X, (i / X) % Y, i / (X * Y)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec::Vec;

    // [3, 5, 7], never a cube. A cubic shape hides every stride bug, because
    // swapping sx and sy is invisible when they are equal.
    const SX: u32 = 3;
    const SY: u32 = 5;
    const SZ: u32 = 7;
    const N: u32 = SX * SY * SZ; // 105

    fn all_coords() -> Vec<[u32; 3]> {
        let mut v = Vec::new();
        for z in 0..SZ {
            for y in 0..SY {
                for x in 0..SX {
                    v.push([x, y, z]);
                }
            }
        }
        v
    }

    #[test]
    fn runtime_shape_round_trips_non_cubic() {
        let shape = RuntimeShape3::new([SX, SY, SZ]).expect("valid shape");
        for i in 0..N {
            assert_eq!(shape.linearize(shape.delinearize(i)), i);
        }
        for p in all_coords() {
            assert_eq!(shape.delinearize(shape.linearize(p)), p);
        }
    }

    #[test]
    fn const_shape_round_trips_non_cubic() {
        let shape = ConstShape3::<SX, SY, SZ>::new();
        for i in 0..N {
            assert_eq!(shape.linearize(shape.delinearize(i)), i);
        }
        for p in all_coords() {
            assert_eq!(shape.delinearize(shape.linearize(p)), p);
        }
    }

    /// The index-order convention, as an executable statement rather than a
    /// sentence in a doc comment.
    #[test]
    fn x_is_the_contiguous_axis() {
        let shape = RuntimeShape3::new([SX, SY, SZ]).expect("valid shape");
        assert_eq!(shape.linearize([1, 0, 0]), 1);
        assert_eq!(shape.linearize([0, 1, 0]), SX);
        assert_eq!(shape.linearize([0, 0, 1]), SX * SY);
    }

    /// Keeps the two implementations from drifting apart.
    #[test]
    fn runtime_and_const_shapes_agree() {
        let runtime = RuntimeShape3::new([SX, SY, SZ]).expect("valid shape");
        let constant = ConstShape3::<SX, SY, SZ>::new();
        assert_eq!(runtime.size(), constant.size());
        assert_eq!(runtime.element_count(), constant.element_count());
        for p in all_coords() {
            assert_eq!(runtime.linearize(p), constant.linearize(p));
        }
        for i in 0..N {
            assert_eq!(runtime.delinearize(i), constant.delinearize(i));
        }
    }

    #[test]
    fn element_count_is_the_product() {
        assert_eq!(
            RuntimeShape3::new([SX, SY, SZ])
                .expect("valid shape")
                .element_count(),
            105
        );
        assert_eq!(ConstShape3::<SX, SY, SZ>::new().element_count(), 105);
    }

    #[test]
    fn runtime_shape_rejects_overflow() {
        let error = RuntimeShape3::new([u32::MAX, 2, 2]).expect_err("should not fit in u32");
        assert!(
            matches!(error, crate::Error::ShapeOverflow { .. }),
            "{error}"
        );
        // The message has to name the numbers, or it is a worse panic.
        let text = alloc::format!("{error}");
        assert!(text.contains("4294967295"), "{text}");
    }

    #[test]
    fn shape3_is_dyn_compatible() {
        const _: Option<&dyn Shape3> = None;
        let shape = RuntimeShape3::new([SX, SY, SZ]).expect("valid shape");
        let dynamic: &dyn Shape3 = &shape;
        assert_eq!(dynamic.linearize([0, 1, 0]), SX);
    }
}
