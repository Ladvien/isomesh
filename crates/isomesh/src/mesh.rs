//! Where triangles go.

use alloc::vec::Vec;

use crate::Real;

/// Receives the triangles an extraction algorithm produces.
///
/// Implemented by [`MeshBuffer`], by an engine wrapper writing straight into a
/// renderer's attribute arrays, and by anything a CAD consumer wants. The point
/// is that an algorithm never allocates the mesh it produces — see
/// [`MeshBuffer::reset`].
///
/// # Winding
///
/// `a, b, c` in [`triangle`](MeshSink::triangle) are **counter-clockwise when
/// viewed from outside the solid**, in a right-handed coordinate system.
/// Equivalently `(b − a) × (c − a)` points the same way as the surface normal,
/// which points away from the solid — the direction of increasing
/// [`Sdf::sample`](crate::Sdf::sample), since negative is inside.
///
/// # Index contract
///
/// [`vertex`](MeshSink::vertex) returns the index by which that vertex may be
/// referenced from [`triangle`](MeshSink::triangle). **Callers must not assume
/// the returned index equals a running counter.** A sink is permitted to weld:
/// to recognise a vertex it has already seen and return that earlier index
/// instead of appending. A welding sink must document its epsilon and its
/// tie-breaking rule, because weld ordering is the classic determinism leak.
///
/// Indices are valid only on the sink that returned them, and only until that
/// sink is reset. Forward references are not permitted: every index passed to
/// `triangle` must have been returned by an earlier `vertex` call on the same
/// sink.
///
/// [`MeshBuffer`] never welds. It always appends, and its returned index is
/// always `vertex_count() - 1`.
pub trait MeshSink {
    /// The scalar positions and normals are written in.
    ///
    /// `f32` for anything destined for a GPU; `f64` for CAD, where the whole
    /// reason the crate is generic is that narrowing the output would throw away
    /// the precision the solve was done in.
    type Scalar: Real;

    /// Append a vertex and return its index. See the trait's index contract.
    fn vertex(&mut self, position: [Self::Scalar; 3], normal: [Self::Scalar; 3]) -> u32;

    /// Append a triangle. See the trait's winding convention.
    fn triangle(&mut self, a: u32, b: u32, c: u32);

    /// A capacity hint.
    ///
    /// Implementors must never *shrink* in response, and must remain correct if
    /// this is never called.
    fn reserve(&mut self, _vertices: usize, _triangles: usize) {}
}

/// The default sink: three parallel `Vec`s, reusable across chunks.
///
/// Reuse is the entire point. A brush stroke re-meshes thousands of chunks and
/// allocation dominates, so the intended lifecycle is one buffer that is
/// [`reset`](MeshBuffer::reset) between chunks rather than one buffer per chunk.
///
/// # Invariants
///
/// The fields are public so a consumer can write straight into a GPU buffer or a
/// mesh attribute array with no intermediate copy. That also means a caller
/// *can* break the invariants: `positions.len() == normals.len()`,
/// `indices.len() % 3 == 0`, and no index `>= positions.len()`. The validity
/// harness checks all three.
///
/// # Example
///
/// ```
/// use isomesh::{MeshBuffer, MeshSink};
///
/// let mut out = MeshBuffer::<f32>::new();
/// for _chunk in 0..2 {
///     out.reset(); // keeps every allocation made so far
///     let a = out.vertex([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
///     let b = out.vertex([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
///     let c = out.vertex([0.0, 0.0, 1.0], [0.0, 1.0, 0.0]);
///     out.triangle(a, b, c);
/// }
/// assert_eq!(out.triangle_count(), 1);
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct MeshBuffer<R: Real = f32> {
    /// One entry per vertex.
    pub positions: Vec<[R; 3]>,
    /// One entry per vertex, parallel to `positions`. Unit length.
    pub normals: Vec<[R; 3]>,
    /// Flat triples: `[a0, b0, c0, a1, b1, c1, …]`.
    pub indices: Vec<u32>,
}

impl<R: Real> MeshBuffer<R> {
    /// An empty buffer that has allocated nothing.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            positions: Vec::new(),
            normals: Vec::new(),
            indices: Vec::new(),
        }
    }

    /// An empty buffer with room for `vertices` vertices and `triangles`
    /// triangles.
    #[must_use]
    pub fn with_capacity(vertices: usize, triangles: usize) -> Self {
        Self {
            positions: Vec::with_capacity(vertices),
            normals: Vec::with_capacity(vertices),
            indices: Vec::with_capacity(triangles * 3),
        }
    }

    /// Truncate to zero length **without releasing capacity**.
    ///
    /// Call this between chunks. Afterwards [`vertex_count`](Self::vertex_count)
    /// is zero and every allocation this buffer has ever made is still held: it
    /// does not reallocate, and it does not touch the heap at all.
    pub fn reset(&mut self) {
        self.positions.clear();
        self.normals.clear();
        self.indices.clear();
    }

    /// Number of vertices.
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }

    /// Number of triangles.
    #[must_use]
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// `true` if no triangles have been written.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }

    /// Release unused capacity.
    ///
    /// **The only method on this type that frees memory.** Everything else is
    /// deliberately non-shrinking, so this is the explicit opt-out rather than
    /// something that happens on its own.
    pub fn shrink_to_fit(&mut self) {
        self.positions.shrink_to_fit();
        self.normals.shrink_to_fit();
        self.indices.shrink_to_fit();
    }

    /// Append another mesh, shifting its indices to follow this one's vertices.
    ///
    /// The reason this exists is [`crate::weld`]: a chunk seam only becomes
    /// weldable once both chunks' vertices are in one buffer, and the index
    /// shift is exactly the step a caller would get wrong. Nothing is welded
    /// here — the vertices are concatenated as they are, duplicates and all.
    pub fn append(&mut self, other: &Self) {
        let base = self.positions.len() as u32;
        self.positions.extend_from_slice(&other.positions);
        self.normals.extend_from_slice(&other.normals);
        self.indices.extend(other.indices.iter().map(|&i| i + base));
    }
}

// Hand-written rather than derived: `derive(Default)` would add a spurious
// `R: Default` bound to the impl.
impl<R: Real> Default for MeshBuffer<R> {
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl<R: Real> MeshSink for MeshBuffer<R> {
    type Scalar = R;

    #[inline]
    fn vertex(&mut self, position: [R; 3], normal: [R; 3]) -> u32 {
        let index = self.positions.len();
        // A `debug_assert!` rather than a check on the hot path: the extractors
        // bound their vertex count against `u32` before they start, and return
        // `Error::IndexSpaceExhausted` if it could not fit. This catches a sink
        // driven directly past the limit, in the builds where that is worth
        // paying for.
        debug_assert!(
            index < u32::MAX as usize,
            "MeshBuffer exceeded the u32 index space at {index} vertices"
        );
        self.positions.push(position);
        self.normals.push(normal);
        index as u32
    }

    #[inline]
    fn triangle(&mut self, a: u32, b: u32, c: u32) {
        self.indices.push(a);
        self.indices.push(b);
        self.indices.push(c);
    }

    #[inline]
    fn reserve(&mut self, vertices: usize, triangles: usize) {
        // `Vec::reserve` never shrinks, so this inherits that property.
        self.positions.reserve(vertices);
        self.normals.reserve(vertices);
        self.indices.reserve(triangles * 3);
    }
}

/// Forwards to the referent, so a sink passed by mutable reference still
/// satisfies the generic `M: MeshSink` bound extraction functions use.
///
/// Note that this forwards [`reserve`](MeshSink::reserve) as well. It must:
/// omitting it would silently substitute the no-op default, so every reservation
/// an algorithm makes would be discarded whenever the sink arrived by reference.
impl<T: MeshSink + ?Sized> MeshSink for &mut T {
    type Scalar = T::Scalar;

    #[inline]
    fn vertex(&mut self, position: [Self::Scalar; 3], normal: [Self::Scalar; 3]) -> u32 {
        T::vertex(self, position, normal)
    }

    #[inline]
    fn triangle(&mut self, a: u32, b: u32, c: u32) {
        T::triangle(self, a, b, c);
    }

    #[inline]
    fn reserve(&mut self, vertices: usize, triangles: usize) {
        T::reserve(self, vertices, triangles);
    }
}

#[cfg(test)]
mod tests {
    // Position and index round-trips are asserted exactly: a sink that stores
    // a vertex must give back that vertex, bit for bit.
    #![allow(clippy::float_cmp)]

    use super::*;

    fn fill(out: &mut MeshBuffer<f32>, vertices: u32) {
        for i in 0..vertices {
            let f = i as f32;
            out.vertex([f, 0.0, 0.0], [0.0, 1.0, 0.0]);
        }
        for i in 0..vertices / 3 {
            out.triangle(3 * i, 3 * i + 1, 3 * i + 2);
        }
    }

    /// I-002's stated acceptance criterion.
    ///
    /// The pointer check is what actually proves "resizes without reallocating":
    /// capacity alone would survive a same-size reallocation to a different
    /// address, which is exactly the thing the reuse contract forbids.
    #[test]
    fn reset_preserves_capacity() {
        let mut out = MeshBuffer::<f32>::new();
        fill(&mut out, 96);

        let caps = (
            out.positions.capacity(),
            out.normals.capacity(),
            out.indices.capacity(),
        );
        let ptrs = (
            out.positions.as_ptr(),
            out.normals.as_ptr(),
            out.indices.as_ptr(),
        );
        // Guard against passing vacuously on a buffer that never allocated.
        assert!(caps.0 > 0 && caps.1 > 0 && caps.2 > 0);

        out.reset();

        assert_eq!(out.vertex_count(), 0);
        assert_eq!(out.triangle_count(), 0);
        assert!(out.is_empty());
        assert_eq!(
            (
                out.positions.capacity(),
                out.normals.capacity(),
                out.indices.capacity()
            ),
            caps
        );
        assert_eq!(
            (
                out.positions.as_ptr(),
                out.normals.as_ptr(),
                out.indices.as_ptr()
            ),
            ptrs
        );
    }

    #[test]
    fn reset_is_idempotent() {
        let mut out = MeshBuffer::<f64>::new();
        for i in 0..12 {
            out.vertex([f64::from(i), 0.0, 0.0], [0.0, 1.0, 0.0]);
        }
        out.reset();
        let cap = out.positions.capacity();
        out.reset();
        assert_eq!(out.positions.capacity(), cap);
        assert_eq!(out.vertex_count(), 0);
    }

    #[test]
    fn reserve_never_shrinks() {
        let mut out = MeshBuffer::<f32>::new();
        out.reserve(1000, 1000);
        let caps = (out.positions.capacity(), out.indices.capacity());
        assert!(caps.0 >= 1000 && caps.1 >= 3000);
        out.reserve(1, 1);
        assert_eq!((out.positions.capacity(), out.indices.capacity()), caps);
    }

    /// Pins the documented non-welding behaviour, which the index contract lets
    /// other sinks vary.
    #[test]
    fn mesh_buffer_never_welds() {
        let mut out = MeshBuffer::<f32>::new();
        let p = [1.0, 2.0, 3.0];
        let n = [0.0, 1.0, 0.0];
        assert_eq!(out.vertex(p, n), 0);
        assert_eq!(out.vertex(p, n), 1);
        assert_eq!(out.vertex_count(), 2);
    }

    #[test]
    fn mesh_sink_is_dyn_compatible() {
        const _: Option<&mut dyn MeshSink<Scalar = f32>> = None;
        let mut out = MeshBuffer::<f32>::new();
        let dynamic: &mut dyn MeshSink<Scalar = f32> = &mut out;
        dynamic.reserve(4, 1);
        let a = dynamic.vertex([0.0; 3], [0.0, 1.0, 0.0]);
        assert_eq!(a, 0);
    }

    /// The second forwarding trap. If `impl MeshSink for &mut T` omitted
    /// `reserve`, this would silently reach the no-op default.
    #[test]
    fn mut_ref_impl_forwards_reserve() {
        struct Counting {
            inner: MeshBuffer<f32>,
            reserves: usize,
        }
        impl MeshSink for Counting {
            type Scalar = f32;
            fn vertex(&mut self, p: [f32; 3], n: [f32; 3]) -> u32 {
                self.inner.vertex(p, n)
            }
            fn triangle(&mut self, a: u32, b: u32, c: u32) {
                self.inner.triangle(a, b, c);
            }
            fn reserve(&mut self, v: usize, t: usize) {
                self.reserves += 1;
                self.inner.reserve(v, t);
            }
        }

        fn through_generic<M: MeshSink>(mut sink: M) {
            sink.reserve(64, 32);
        }

        let mut sink = Counting {
            inner: MeshBuffer::new(),
            reserves: 0,
        };
        through_generic(&mut sink);
        assert_eq!(sink.reserves, 1);
        assert!(sink.inner.positions.capacity() >= 64);
    }

    #[test]
    fn indices_are_flat_triples() {
        let mut out = MeshBuffer::<f32>::new();
        out.triangle(0, 1, 2);
        assert_eq!(out.indices, [0, 1, 2]);
        assert_eq!(out.triangle_count(), 1);
    }

    #[test]
    fn shrink_to_fit_is_the_only_release() {
        let mut out = MeshBuffer::<f32>::new();
        fill(&mut out, 96);
        out.reset();
        assert!(out.positions.capacity() > 0);
        out.shrink_to_fit();
        assert_eq!(out.positions.capacity(), 0);
        assert_eq!(out.normals.capacity(), 0);
        assert_eq!(out.indices.capacity(), 0);
    }
}
