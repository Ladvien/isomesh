//! The one error type.
//!
//! Every public entry point that can reject its input returns
//! [`Result`](crate::Result) rather than panicking. Internal invariants — things
//! that are this crate's own bugs rather than a caller's — stay as
//! `debug_assert!`, which compiles out in release and so cannot abort a shipped
//! build.
//!
//! "Fail loudly" is still the rule; a typed error at the call site is louder
//! than an abort, because the caller can print it, log it, or attach it to the
//! chunk that produced it. What is *not* acceptable is a degraded substitute —
//! no extractor here silently clamps a bad grid or invents a spacing.
//!
//! # Where errors are made impossible instead
//!
//! [`ValidateConfig`](crate::validate::ValidateConfig) has private fields and
//! one checked constructor, so a meaningless threshold cannot be constructed at
//! all and the validator needs no runtime check. Making a state unrepresentable
//! beats reporting it.

use core::fmt;

/// What went wrong.
///
/// Every variant carries the numbers that produced it, so the message is
/// actionable without reaching for a debugger.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum Error {
    /// A grid whose sample count does not fit in `u32`.
    ///
    /// A silent wrap here would alias distinct cells onto the same index and
    /// produce a mesh that is self-consistent and simply not the field's — which
    /// no validity test can catch.
    ShapeOverflow {
        /// The requested size.
        size: [u32; 3],
        /// What `sx * sy * sz` actually comes to.
        product: u64,
    },

    /// Fewer than two samples on some axis, so there is no cell to extract from.
    GridTooSmall {
        /// The requested size.
        size: [u32; 3],
    },

    /// The output would need more vertices than a `u32` index can address.
    IndexSpaceExhausted {
        /// The vertex demand that failed to fit: an a-priori upper bound where
        /// one exists, or the count at which subgrid's running guard tripped.
        needed: u64,
    },

    /// A triangle index referring to a vertex the buffer does not have.
    ///
    /// [`MeshBuffer`](crate::MeshBuffer)'s fields are public, so a caller can
    /// construct this. It is rejected at the door rather than part-way through
    /// an operation, which would leave the buffer half-rewritten.
    IndexOutOfRange {
        /// Offset into `indices` of the offending entry.
        at: u64,
        /// The index found there.
        index: u32,
        /// How many vertices the buffer has.
        vertices: u64,
    },

    /// A cell size that is not finite and positive.
    ///
    /// Thresholds in this crate are relative to the grid spacing, so a
    /// meaningless spacing makes every threshold meaningless with it.
    InvalidCellSize {
        /// The value given.
        value: f64,
    },

    /// A weld tolerance that is not finite and positive.
    ///
    /// Zero is rejected rather than treated as "weld only exact matches",
    /// because M-32 measured that seam vertices are bit-identical only for
    /// power-of-two cell sizes — an exact-match weld is precisely the thing that
    /// works on the fixture and fails in the field.
    InvalidWeldEpsilon {
        /// The value given.
        value: f64,
    },

    /// A spacing that does not describe the mesh it was given.
    ///
    /// Reported rather than absorbed: a broadphase grid finer than the mesh
    /// would grow until it exhausted memory, and guessing a better spacing would
    /// silently change what the measurement means.
    CellSizeMismatch {
        /// Index of the triangle that triggered it.
        triangle: u64,
        /// How many grid cells that triangle's bounding box spans.
        cells: u128,
        /// The spacing given.
        cell_size: f64,
    },

    /// A sweep whose bilinear denominator vanishes on one of its two faces.
    ///
    /// `A + C - B - D` being zero means that face's bilinear function has no
    /// saddle point, so
    /// [`SweptFaces::test`](crate::marching_cubes::interior::SweptFaces::test)'s
    /// criterion -- *is there a cutting plane whose saddle is positive* -- has
    /// no answer there rather than a hard one.
    ///
    /// It cannot happen on an ambiguous face, where one diagonal is strictly
    /// negative and the other non-negative. Reported rather than defaulted for
    /// that reason: reaching it means the interior test was applied to a face
    /// that is not ambiguous, and a default would hide the caller's mistake
    /// behind a plausible topology.
    DegenerateSweep,

    /// A vertex whose normal cannot be derived, so there is nothing to normalise.
    ///
    /// A zero or non-finite field gradient, or — under
    /// [`AreaWeightedFaces`](crate::normals::NormalStrategy::AreaWeightedFaces) —
    /// a vertex whose incident triangles have no area, or which no triangle
    /// references at all.
    ///
    /// Reported rather than substituted. Writing `[0, 0, 1]` there gives a mesh
    /// that renders and is wrong in a way nothing downstream can attribute, which
    /// is exactly the class of bug the one-path rule exists to prevent.
    DegenerateNormal {
        /// Index of the offending vertex.
        vertex: u64,
    },

    /// Subgrid Marching Tetrahedra met a tetrahedron it could not triangulate.
    ///
    /// Every case §3.2 defines is implemented, so this is not "an unsupported
    /// configuration" — it is a defect, and it names the cell and tetrahedron so
    /// the offending edge coordinates can be recovered. Reported rather than
    /// skipped: a silently dropped tetrahedron is a hole in the mesh, and a hole
    /// is indistinguishable from the thin feature this extractor exists to
    /// resolve.
    SubgridUnfilled {
        /// The cell whose tetrahedron failed, in grid coordinates.
        cell: [u32; 3],
        /// Which of the six tetrahedra of that cell.
        tet: u8,
        /// What [`fill`](crate::subgrid::surface::fill) reported, as its
        /// `Debug` form.
        reason: &'static str,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ShapeOverflow { size, product } => write!(
                f,
                "grid {}x{}x{} has {product} samples, which does not fit in u32",
                size[0], size[1], size[2]
            ),
            Self::GridTooSmall { size } => write!(
                f,
                "grid {}x{}x{} has an axis with fewer than two samples, so it contains no cell",
                size[0], size[1], size[2]
            ),
            Self::IndexSpaceExhausted { needed } => write!(
                f,
                "extraction could need {needed} vertices, beyond the u32 index space"
            ),
            Self::IndexOutOfRange {
                at,
                index,
                vertices,
            } => write!(
                f,
                "indices[{at}] is {index}, but the buffer has {vertices} vertices"
            ),
            Self::InvalidCellSize { value } => {
                write!(f, "cell size must be finite and positive, got {value}")
            }
            Self::InvalidWeldEpsilon { value } => {
                write!(f, "weld epsilon must be finite and positive, got {value}")
            }
            Self::CellSizeMismatch {
                triangle,
                cells,
                cell_size,
            } => write!(
                f,
                "triangle {triangle} spans {cells} grid cells at cell size {cell_size}; \
                 that spacing does not describe this mesh"
            ),
            Self::DegenerateSweep => write!(
                f,
                "a swept face has a zero bilinear denominator, so it has no saddle point \
                 and the interior test has nothing to evaluate"
            ),
            Self::DegenerateNormal { vertex } => write!(
                f,
                "vertex {vertex} has no normal to derive: a zero gradient, or no incident area"
            ),
            Self::SubgridUnfilled { cell, tet, reason } => write!(
                f,
                "subgrid marching tetrahedra could not fill tetrahedron {tet} of cell \
                 [{}, {}, {}]: {reason}",
                cell[0], cell[1], cell[2]
            ),
        }
    }
}

impl core::error::Error for Error {}

/// This crate's result type.
pub type Result<T> = core::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::format;

    /// Every variant has to say what happened *and* with which numbers, or the
    /// error is a worse version of a panic.
    #[test]
    fn messages_carry_their_numbers() {
        let cases = [
            (
                Error::ShapeOverflow {
                    size: [70000, 70000, 2],
                    product: 9_800_000_000,
                },
                "9800000000",
            ),
            (Error::GridTooSmall { size: [1, 4, 4] }, "1x4x4"),
            (
                Error::IndexSpaceExhausted {
                    needed: 5_000_000_000,
                },
                "5000000000",
            ),
            (Error::InvalidCellSize { value: 0.0 }, "got 0"),
            (Error::InvalidWeldEpsilon { value: -1.0 }, "got -1"),
            (
                Error::CellSizeMismatch {
                    triangle: 7,
                    cells: 900,
                    cell_size: 0.001,
                },
                "900",
            ),
        ];
        for (error, expected) in cases {
            let text = format!("{error}");
            assert!(text.contains(expected), "{error:?} rendered as {text:?}");
        }
    }

    #[test]
    fn errors_compare_by_value() {
        assert_eq!(
            Error::GridTooSmall { size: [1, 2, 2] },
            Error::GridTooSmall { size: [1, 2, 2] }
        );
        assert_ne!(
            Error::GridTooSmall { size: [1, 2, 2] },
            Error::GridTooSmall { size: [2, 1, 2] }
        );
    }
}
