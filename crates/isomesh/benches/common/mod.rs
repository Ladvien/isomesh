//! Setup shared by the benchmark targets.
//!
//! One definition of "this field at N³", so that a criterion benchmark and the
//! resolution sweep cannot silently disagree about which grid they measured.
//!
//! # Why the blanket `dead_code` allow
//!
//! Cargo compiles `mod common;` **once per bench target**, so every helper is
//! unused from the point of view of every bench that does not happen to call it
//! — `grid` is dead in `experiment_p8`, `experiment` is dead in `shootout`, and
//! so on. The alternative is one `#[allow]` per item, which grows with the file
//! and says the same thing eight times.

#![allow(
    dead_code,
    reason = "compiled once per bench target, so every helper is unused in most of them"
)]

pub(crate) mod beta;
pub(crate) mod boolean;
pub(crate) mod experiment;
pub(crate) mod heat;
pub(crate) mod lattice;
pub(crate) mod metric;
pub(crate) mod poly;
pub(crate) mod tpms;
pub(crate) mod wedge;

// Hardware counters come from `perf_event_open`, a Linux system call with no
// macOS equivalent a bench can reach. Callers on other platforms either refuse
// (`experiment_p12`) or record the columns as `unavailable` (`family`).
#[cfg(target_os = "linux")]
pub(crate) mod counters;

use isomesh::fields::ReferenceField;
use isomesh::{Real, RuntimeShape3, Sdf};

/// The grid a reference field is meant to be sampled on at `samples` per axis.
///
/// `shape` counts **samples**, so `n` samples span `n − 1` cells. This is the
/// same convention `mc/tests.rs` uses, kept identical on purpose: a benchmark
/// and a test that both say "64³" have to mean the same grid, or the accuracy
/// numbers and the timing numbers describe different things.
pub(crate) fn grid<R, F>(field: &F, samples: u32) -> (RuntimeShape3, [R; 3], R)
where
    R: Real,
    F: ReferenceField + Sdf<Scalar = R>,
{
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / R::from_f64(f64::from(samples - 1));
    let shape = RuntimeShape3::new([samples; 3]).expect("benchmark grid fits u32");
    (shape, lo, cell_size)
}
