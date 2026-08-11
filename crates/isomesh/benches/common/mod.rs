//! Setup shared by the benchmark targets.
//!
//! One definition of "this field at N³", so that a criterion benchmark and the
//! resolution sweep cannot silently disagree about which grid they measured.

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
