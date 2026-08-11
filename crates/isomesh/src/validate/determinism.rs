//! Checking that meshing the same thing twice produces the same bytes.
//!
//! # Why `==` is the wrong comparison
//!
//! "Byte-identical" is not what float equality means, and it is wrong in *both*
//! directions: `+0.0 == -0.0` is true although the bit patterns differ, and
//! `NaN == NaN` is false although they may be the same bits. Both cases are
//! reachable — a sign flip on a zero coordinate is exactly what a reordered
//! summation produces, and it is exactly the sort of thing a golden hash would
//! then catch far downstream instead of here.
//!
//! Comparison therefore goes through [`Real::total_cmp`], the IEEE `totalOrder`
//! predicate, which compares the full bit pattern. `Equal` means identical bits
//! and nothing else.
//!
//! # What is checked
//!
//! Two things, because the crate's whole API is built on reusing output buffers:
//!
//! 1. **The same call twice** produces identical output. This is what catches
//!    iteration order leaking into vertex order.
//! 2. **A reused buffer** produces the same output as a fresh one. Nothing else
//!    in the test suite checks that an extractor is correct after
//!    [`MeshBuffer::reset`], and every algorithm here is meant to be driven that
//!    way.

use core::cmp::Ordering;
use core::fmt;

use crate::{MeshBuffer, Real};

/// Which pair of runs disagreed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RunPair {
    /// The same call, made twice into two fresh buffers.
    ///
    /// A failure here is non-determinism in the extractor itself — most often
    /// an unordered container's iteration order reaching the output.
    RepeatedCall,
    /// A fresh buffer against one that was used, [`MeshBuffer::reset`], and used
    /// again.
    ///
    /// A failure here means the result depends on the output buffer's prior
    /// state, which breaks the reuse contract the whole API is built on.
    ReusedBuffer,
}

impl fmt::Display for RunPair {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RepeatedCall => write!(f, "the same call made twice"),
            Self::ReusedBuffer => write!(f, "a reused buffer against a fresh one"),
        }
    }
}

/// The first place two runs differed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Divergence<R: Real> {
    /// Different numbers of vertices.
    VertexCount {
        /// From the first run.
        first: u64,
        /// From the second.
        second: u64,
    },
    /// Different numbers of indices.
    IndexCount {
        /// From the first run.
        first: u64,
        /// From the second.
        second: u64,
    },
    /// A position differs in its bits.
    Position {
        /// Vertex it belongs to.
        vertex: u64,
        /// Component, `0..3`.
        axis: usize,
        /// From the first run.
        first: R,
        /// From the second.
        second: R,
    },
    /// A normal differs in its bits.
    Normal {
        /// Vertex it belongs to.
        vertex: u64,
        /// Component, `0..3`.
        axis: usize,
        /// From the first run.
        first: R,
        /// From the second.
        second: R,
    },
    /// An index differs.
    Index {
        /// Position in the flat index buffer.
        at: u64,
        /// From the first run.
        first: u32,
        /// From the second.
        second: u32,
    },
}

impl<R: Real> fmt::Display for Divergence<R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::VertexCount { first, second } => {
                write!(f, "vertex count: {first} then {second}")
            }
            Self::IndexCount { first, second } => {
                write!(f, "index count: {first} then {second}")
            }
            Self::Position {
                vertex,
                axis,
                first,
                second,
            } => write!(f, "position[{vertex}][{axis}]: {first:?} then {second:?}"),
            Self::Normal {
                vertex,
                axis,
                first,
                second,
            } => write!(f, "normal[{vertex}][{axis}]: {first:?} then {second:?}"),
            Self::Index { at, first, second } => {
                write!(f, "index[{at}]: {first} then {second}")
            }
        }
    }
}

/// The result of running an extractor several times.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeterminismReport<R: Real> {
    /// Vertices the first run produced.
    pub vertices: u64,
    /// Triangles the first run produced.
    pub triangles: u64,
    /// The first disagreement found, if any.
    pub divergence: Option<(RunPair, Divergence<R>)>,
}

impl<R: Real> DeterminismReport<R> {
    /// `true` when every run produced identical bytes.
    #[must_use]
    pub fn is_deterministic(&self) -> bool {
        self.divergence.is_none()
    }

    /// The loud path, for tests that want to stop.
    ///
    /// # Panics
    ///
    /// If any run diverged, naming which pair of runs and exactly where.
    pub fn panic_if_divergent(&self) {
        if let Some((pair, divergence)) = self.divergence {
            panic!("non-deterministic output between {pair}: {divergence}");
        }
    }
}

impl<R: Real> fmt::Display for DeterminismReport<R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "determinism report")?;
        writeln!(f, "  vertices                 {:8}", self.vertices)?;
        writeln!(f, "  triangles                {:8}", self.triangles)?;
        match self.divergence {
            None => write!(f, "  DETERMINISTIC"),
            Some((pair, d)) => write!(f, "  DIVERGED between {pair}\n    {d}"),
        }
    }
}

/// Bit-level inequality.
///
/// [`Real::total_cmp`] is the IEEE `totalOrder` predicate and compares full bit
/// patterns, so `Equal` means identical bits. See the module docs for why `==`
/// cannot be used here.
#[inline]
fn differs<R: Real>(a: R, b: R) -> bool {
    a.total_cmp(&b) != Ordering::Equal
}

fn compare<R: Real>(first: &MeshBuffer<R>, second: &MeshBuffer<R>) -> Option<Divergence<R>> {
    if first.positions.len() != second.positions.len() {
        return Some(Divergence::VertexCount {
            first: first.positions.len() as u64,
            second: second.positions.len() as u64,
        });
    }
    if first.indices.len() != second.indices.len() {
        return Some(Divergence::IndexCount {
            first: first.indices.len() as u64,
            second: second.indices.len() as u64,
        });
    }

    for (v, (a, b)) in first.positions.iter().zip(&second.positions).enumerate() {
        for axis in 0..3 {
            if differs(a[axis], b[axis]) {
                return Some(Divergence::Position {
                    vertex: v as u64,
                    axis,
                    first: a[axis],
                    second: b[axis],
                });
            }
        }
    }

    // Normals are compared only where both runs produced them; a length
    // mismatch against `positions` is the validity harness's business, not this
    // one's.
    for (v, (a, b)) in first.normals.iter().zip(&second.normals).enumerate() {
        for axis in 0..3 {
            if differs(a[axis], b[axis]) {
                return Some(Divergence::Normal {
                    vertex: v as u64,
                    axis,
                    first: a[axis],
                    second: b[axis],
                });
            }
        }
    }
    if first.normals.len() != second.normals.len() {
        return Some(Divergence::VertexCount {
            first: first.normals.len() as u64,
            second: second.normals.len() as u64,
        });
    }

    for (i, (a, b)) in first.indices.iter().zip(&second.indices).enumerate() {
        if a != b {
            return Some(Divergence::Index {
                at: i as u64,
                first: *a,
                second: *b,
            });
        }
    }

    None
}

/// Run `extract` three times and report the first bit-level disagreement.
///
/// The intended call shape is one line per algorithm ticket:
///
/// ```
/// use isomesh::{MeshBuffer, MeshSink};
/// use isomesh::validate::check_determinism;
///
/// let report = check_determinism(|out: &mut MeshBuffer<f32>| {
///     let a = out.vertex([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
///     let b = out.vertex([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
///     let c = out.vertex([0.0, 0.0, 1.0], [0.0, 1.0, 0.0]);
///     out.triangle(a, b, c);
/// });
///
/// assert!(report.is_deterministic());
/// assert_eq!(report.triangles, 1);
/// ```
///
/// Three runs rather than two: two fresh buffers to catch non-determinism in the
/// extractor, and one reused buffer to catch output that depends on the buffer's
/// prior state. The second is not paranoia — every algorithm in this crate is
/// meant to be driven by resetting one buffer across thousands of chunks, and
/// nothing else checks that it survives being driven that way.
pub fn check_determinism<R, F>(mut extract: F) -> DeterminismReport<R>
where
    R: Real,
    F: FnMut(&mut MeshBuffer<R>),
{
    let mut first = MeshBuffer::<R>::new();
    extract(&mut first);

    let mut report = DeterminismReport {
        vertices: first.positions.len() as u64,
        triangles: (first.indices.len() / 3) as u64,
        divergence: None,
    };

    let mut second = MeshBuffer::<R>::new();
    extract(&mut second);
    if let Some(d) = compare(&first, &second) {
        report.divergence = Some((RunPair::RepeatedCall, d));
        return report;
    }

    let mut reused = MeshBuffer::<R>::new();
    extract(&mut reused);
    reused.reset();
    extract(&mut reused);
    if let Some(d) = compare(&first, &reused) {
        report.divergence = Some((RunPair::ReusedBuffer, d));
    }

    report
}

#[cfg(test)]
mod tests;
