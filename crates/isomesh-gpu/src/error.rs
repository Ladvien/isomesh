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
    /// Timestamp queries were asked for and are not usable on this device.
    ///
    /// Ticket: R-069 (P-71). Raised when the device was not created with
    /// [`wgpu::Features::TIMESTAMP_QUERY`], when the queue reports a timestamp
    /// period of zero — a driver that advertises the feature and does not
    /// implement it — or when a resolved pair ends before it begins.
    ///
    /// **It refuses rather than degrading, on purpose.** Falling back to
    /// `Instant::now()` under a GPU-side column name would report a number that
    /// was named and not measured, which is the failure the experiment harness's
    /// missing-column panic exists to prevent, arriving by a different door.
    TimestampsUnsupported,
    /// The adapter does not advertise a feature the caller asked for.
    ///
    /// R-068. `Gpu::open` previously returned [`Self::TimestampsUnsupported`]
    /// for **any** missing feature, which was correct while timestamps were its
    /// only caller and became a lie the moment `with_subgroups` existed: a
    /// device without `SUBGROUP` would have reported that timestamps were
    /// unavailable. Carrying the missing set means the error names what is
    /// actually absent, and adding a third capability cannot reintroduce the
    /// bug.
    FeaturesUnsupported {
        /// The features the adapter lacks.
        missing: wgpu::Features,
    },
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
    /// The device cannot run mesh shaders.
    ///
    /// Returned rather than substituting a vertex-buffer pipeline: a caller
    /// told "drawing" while a different pipeline ran has been misinformed
    /// about the one thing they asked for.
    MeshShadersUnavailable,
    /// A read-back range was not a whole number of elements.
    UnalignedReadback {
        /// Bytes asked for.
        bytes: u64,
        /// Bytes per element.
        stride: u64,
    },
    /// A scan was asked for more elements than one dispatch can cover.
    ///
    /// Level 0 dispatches one workgroup per `PrefixScan::BLOCK` elements, and
    /// `max_compute_workgroups_per_dimension` caps how many that can be.
    ScanTooLong {
        /// Elements the caller asked to scan.
        elements: u32,
        /// The most this device will scan in one call.
        max: u32,
    },
    /// [`crate::DeferredGeometry`] was asked to hold one more read-back than it
    /// has room for.
    ///
    /// R-071 (P-71 C3). Its own variant rather than [`Self::DeviceLost`] or a
    /// silent drop, because it is the one error here that is **not** a failure:
    /// it is the queue telling a frame-budget scheduler that this frame's
    /// submissions have caught up with last frame's collections, which is the
    /// signal a caller uses to stop submitting. A dropped read-back would be a
    /// chunk of geometry that was extracted, paid for and never delivered, and
    /// nothing downstream could tell that from a chunk with no surface in it.
    DeferredQueueFull {
        /// In-flight read-backs the queue was built to hold.
        capacity: usize,
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
            Self::TimestampsUnsupported => f.write_str(
                "timestamp queries are not usable on this device: the feature was \
                 not enabled, the timestamp period is zero, or a span ended before \
                 it began",
            ),
            Self::FeaturesUnsupported { missing } => {
                write!(f, "the adapter does not advertise {missing:?}")
            }
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
            Self::MeshShadersUnavailable => f.write_str(
                "mesh shaders need EXPERIMENTAL_MESH_SHADER on a Vulkan device; \
                 naga compiles WGSL mesh stages for SPIR-V only",
            ),
            Self::UnalignedReadback { bytes, stride } => write!(
                f,
                "read-back of {bytes} bytes is not a whole number of {stride}-byte elements"
            ),
            Self::ScanTooLong { elements, max } => write!(
                f,
                "cannot scan {elements} elements; this device dispatches at most {max}"
            ),
            Self::DeferredQueueFull { capacity } => {
                write!(f, "the deferred read-back queue is full at {capacity}")
            }
        }
    }
}

impl core::error::Error for Error {}
