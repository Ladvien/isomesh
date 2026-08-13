//! What can go wrong, reported rather than panicked.

use core::fmt;

/// The result of a fallible `isomesh-gpu` operation.
pub type Result<T> = core::result::Result<T, Error>;

/// Everything this crate can fail at.
///
/// One variant per genuinely distinct cause. There is deliberately no
/// catch-all `Other(String)`: a caller that cannot tell "no adapter" from
/// "wrong buffer length" cannot react to either, and a variant that carries a
/// message is a way of not deciding what the error is.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Error {
    /// A grid axis had fewer than two samples, so it contains no cells.
    DegenerateGrid {
        /// The offending sample counts.
        samples: [u32; 3],
    },
    /// The cell size was not a finite positive number.
    InvalidCellSize,
    /// The grid origin was not finite.
    InvalidOrigin,
    /// The grid describes more samples than a `u64` byte count can address.
    GridTooLarge {
        /// The offending sample counts.
        samples: [u32; 3],
    },
    /// An upload's slice length disagreed with the grid it is being written to.
    SampleCountMismatch {
        /// What the grid says it holds.
        expected: u64,
        /// What the caller passed.
        got: u64,
    },
    /// No adapter matched the request. There is no software fallback here on
    /// purpose: a caller told "GPU extraction is running" while it is on a
    /// CPU reference driver has been misinformed about the only thing that
    /// matters.
    NoAdapter,
    /// The adapter refused to create a device with the requested limits.
    DeviceUnavailable,
    /// Mapping a buffer for read-back failed.
    MapFailed,
    /// The device disconnected or a submission never completed.
    DeviceLost,
    /// A shader module was included, or composed, but never registered.
    ShaderModuleMissing {
        /// The name that was asked for.
        name: String,
    },
    /// A shader module includes itself, through some chain.
    ShaderCircularInclude {
        /// The module the cycle closes on.
        name: String,
    },
    /// A preprocessor directive was malformed or unbalanced.
    ShaderDirective {
        /// Module the directive is in.
        module: String,
        /// One-based line number within that module.
        line: usize,
    },
    /// A read-back range was not a whole number of elements.
    UnalignedReadback {
        /// Bytes asked for.
        bytes: u64,
        /// Bytes per element.
        stride: u64,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DegenerateGrid { samples } => write!(
                f,
                "every grid axis needs at least 2 samples to hold a cell, got {samples:?}"
            ),
            Self::InvalidCellSize => f.write_str("cell size must be finite and greater than zero"),
            Self::InvalidOrigin => f.write_str("grid origin must be finite"),
            Self::GridTooLarge { samples } => {
                write!(f, "grid {samples:?} does not fit in an addressable buffer")
            }
            Self::SampleCountMismatch { expected, got } => {
                write!(f, "grid holds {expected} samples, got {got}")
            }
            Self::NoAdapter => f.write_str("no GPU adapter matched the request"),
            Self::DeviceUnavailable => f.write_str("the adapter would not create a device"),
            Self::MapFailed => f.write_str("mapping a buffer for read-back failed"),
            Self::DeviceLost => f.write_str("the device was lost or a submission never completed"),
            Self::ShaderModuleMissing { name } => {
                write!(f, "no shader module registered as `{name}`")
            }
            Self::ShaderCircularInclude { name } => {
                write!(f, "shader module `{name}` includes itself")
            }
            Self::ShaderDirective { module, line } => {
                write!(f, "malformed or unbalanced directive at {module}:{line}")
            }
            Self::UnalignedReadback { bytes, stride } => write!(
                f,
                "read-back of {bytes} bytes is not a whole number of {stride}-byte elements"
            ),
        }
    }
}

impl core::error::Error for Error {}
