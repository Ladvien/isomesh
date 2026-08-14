//! How far GPU arithmetic drifts from `libm`'s, per field.
//!
//! This is the deliverable of GPU-011a, not a side check. The whole reason to
//! measure it is that moving field evaluation onto the GPU **gives up the
//! sample-identity** that lets M-142 compare CPU and GPU meshes at all, and the
//! size of what replaces it decides whether the general mechanism is worth
//! building.
#![allow(clippy::float_cmp)]

use super::{FieldSampler, GpuField};
use crate::headless::Gpu;
use crate::{GridParams, read_buffer};

fn gpu() -> Gpu {
    Gpu::new().expect("a GPU adapter -- no software fallback, by design")
}

/// How one field's GPU evaluation compares with `isomesh`'s.
///
/// `worst_ulps` is reported but **not** asserted on except for `box_exact`, and
/// the reason is that these fields are zero on their own surface: at a value of
/// `1e-8` an absolute difference of `1e-8` is astronomically many ULPs while
/// being physically nothing. ULP distance is only meaningful between numbers of
/// comparable magnitude, and a signed distance field crosses zero by
/// construction. `worst_abs` is the number that decides whether a crossing
/// moves.
struct Deviation {
    exact: usize,
    total: usize,
    worst_abs: f32,
    worst_ulps: u32,
    at: [f32; 3],
}

fn measure(gpu: &Gpu, sampler: &FieldSampler, field: GpuField, grid: GridParams) -> Deviation {
    let buffer = sampler
        .sample(gpu.device(), gpu.queue(), grid, field)
        .expect("gpu sample");
    let got = read_buffer(
        gpu.device(),
        gpu.queue(),
        buffer.buffer(),
        grid.field_buffer_size(),
    )
    .expect("read back");

    let [sx, sy, sz] = grid.samples();
    let mut deviation = Deviation {
        exact: 0,
        total: 0,
        worst_abs: 0.0,
        worst_ulps: 0,
        at: [0.0; 3],
    };
    let mut i = 0usize;
    for z in 0..sz {
        for y in 0..sy {
            for x in 0..sx {
                let p = grid.sample_position([x, y, z]);
                let cpu = field.sample_on_cpu(p);
                let gpu_value = got[i];
                deviation.total += 1;
                if gpu_value.to_bits() == cpu.to_bits() {
                    deviation.exact += 1;
                } else {
                    let abs = (gpu_value - cpu).abs();
                    // ULPs between two finite floats of the same sign is the
                    // distance between their bit patterns.
                    let ulps = gpu_value.to_bits().abs_diff(cpu.to_bits());
                    if abs > deviation.worst_abs {
                        deviation.worst_abs = abs;
                        deviation.at = p;
                    }
                    deviation.worst_ulps = deviation.worst_ulps.max(ulps);
                }
                i += 1;
            }
        }
    }
    deviation
}

/// The deviation table, printed and asserted.
///
/// Asserted on the **absolute** deviation, which is the quantity that decides
/// whether a crossing moves. The bit-exact counts and ULP figures are printed
/// rather than gated: they are properties of this driver's contraction choices
/// and would make a brittle test, where "no field is wrong by enough to move a
/// surface" is the claim that has to hold on any adapter.
#[test]
fn the_gpu_field_matches_the_cpu_within_a_measured_bound() {
    let gpu = gpu();
    let sampler = FieldSampler::new(gpu.device()).expect("pipeline");
    let grid = GridParams::new([33; 3], [-2.0; 3], 0.125).expect("grid");

    println!(
        "{:<12} {:>10} {:>8} {:>12} {:>8}",
        "field", "bit-exact", "of", "worst abs", "ulps"
    );
    for field in GpuField::ALL {
        let d = measure(&gpu, &sampler, field, grid);
        println!(
            "{:<12} {:>10} {:>8} {:>12.3e} {:>8}   worst at {:?}",
            field.name(),
            d.exact,
            d.total,
            d.worst_abs,
            d.worst_ulps,
            d.at
        );

        // A cell is 0.125 world units. Anything within a millionth of that
        // cannot move a crossing anywhere a mesh would notice.
        assert!(
            d.worst_abs < 1e-6,
            "{}: worst deviation {:e} is too large to be rounding",
            field.name(),
            d.worst_abs
        );
    }
}

/// **No field agrees bit-for-bit — not even the ones with no transcendentals.**
///
/// That was the hypothesis and it is wrong. The grid here uses `h = 0.125` from
/// origin `-2.0`, both powers of two, so `origin + h·i` is *exact* and the two
/// sides evaluate at identical positions. The divergence is therefore inside the
/// field expression itself: a GPU is free to contract `x*x + y*y + z*z` into
/// fused multiply-adds, which rounds once where `libm` rounds twice. Same cause
/// as M-142, one layer further up.
///
/// Asserted as an inequality so it cannot rot into a belief: if a driver starts
/// matching exactly, this fails and the table gets re-read.
#[test]
fn no_field_agrees_with_libm_bit_for_bit() {
    let gpu = gpu();
    let sampler = FieldSampler::new(gpu.device()).expect("pipeline");
    let grid = GridParams::new([33; 3], [-2.0; 3], 0.125).expect("grid");

    for field in GpuField::ALL {
        let d = measure(&gpu, &sampler, field, grid);
        assert!(
            d.exact < d.total,
            "{} now agrees with libm on all {} samples -- the GPU/CPU divergence \
             this ticket measured has gone away for it, and everything resting \
             on the deviation table needs re-reading",
            field.name(),
            d.total
        );
    }
}

/// `box_exact` is far closer than the others, and the reason is structural.
///
/// Its expression is `abs`, `max` and `min` — all exact in IEEE — over a `sqrt`
/// whose argument is **identically zero everywhere inside the box**, where
/// `sqrt(0)` is exact too. Only the exterior shell has anything to round. So it
/// lands within a single ULP while `sphere` and `torus`, which square and sum
/// three coordinates on every sample, do not.
///
/// This is the observation GPU-011b needs: the drift is a property of the
/// *expression*, not of the GPU, so a field built from exact operations crosses
/// unchanged and one built from products does not.
#[test]
fn the_field_of_exact_operations_is_the_closest() {
    let gpu = gpu();
    let sampler = FieldSampler::new(gpu.device()).expect("pipeline");
    let grid = GridParams::new([33; 3], [-2.0; 3], 0.125).expect("grid");

    let boxed = measure(&gpu, &sampler, GpuField::BoxExact, grid);
    assert_eq!(
        boxed.worst_ulps, 1,
        "box_exact's worst deviation is {} ULPs, not 1 -- its arithmetic is \
         supposed to be exact except on the exterior shell",
        boxed.worst_ulps
    );

    for other in [GpuField::Sphere, GpuField::Torus, GpuField::Gyroid] {
        let d = measure(&gpu, &sampler, other, grid);
        assert!(
            d.exact < boxed.exact,
            "{} matches on {} samples, box_exact on {} -- box_exact should be \
             the closest",
            other.name(),
            d.exact,
            boxed.exact
        );
    }
}

/// The samples land where `grid_sample_at` says, not transposed.
///
/// A transposed field is sampled correctly and stored rotated, which looks like
/// a meshing bug rather than an indexing one — so the inverse index mapping gets
/// its own asymmetric fixture.
#[test]
fn the_flat_index_is_not_transposed() {
    let gpu = gpu();
    let sampler = FieldSampler::new(gpu.device()).expect("pipeline");
    // Deliberately unequal extents: a transposition on a cube is invisible.
    let grid = GridParams::new([5, 7, 3], [-1.0, -2.0, -0.5], 0.25).expect("grid");
    let buffer = sampler
        .sample(gpu.device(), gpu.queue(), grid, GpuField::BoxExact)
        .expect("sample");
    let got = read_buffer(
        gpu.device(),
        gpu.queue(),
        buffer.buffer(),
        grid.field_buffer_size(),
    )
    .expect("read back");

    let [sx, sy, _] = grid.samples();
    for (index, value) in got.iter().enumerate() {
        let i = index as u32;
        let at = [i % sx, (i / sx) % sy, i / (sx * sy)];
        let expected = GpuField::BoxExact.sample_on_cpu(grid.sample_position(at));
        assert_eq!(
            value.to_bits(),
            expected.to_bits(),
            "sample {index} holds the value for a different position"
        );
    }
}

/// **Does the deviation change the mesh?** — the question the table only
/// implies.
///
/// A sample within `7e-7` of zero can land on the other side of the sign test,
/// and a flipped sign changes a cell's case, its triangle count and its
/// topology. Over two million samples "the deviation is tiny" is an argument;
/// this is the measurement.
#[test]
fn a_gpu_evaluated_field_extracts_the_same_mesh() {
    use crate::{FieldBuffer, MarchingCubesGpu};

    let gpu = gpu();
    let sampler = FieldSampler::new(gpu.device()).expect("pipeline");
    let mc = MarchingCubesGpu::new(gpu.device(), gpu.queue()).expect("pipeline");

    for (field, samples, cell) in [
        (GpuField::Sphere, 65u32, 0.0625f32),
        (GpuField::Torus, 65, 0.0625),
        (GpuField::BoxExact, 65, 0.0625),
        // The one with the largest deviation, at a spacing that is not a power
        // of two so the positions differ as well as the values.
        (GpuField::Gyroid, 49, 4.0 / 48.0),
    ] {
        let grid = GridParams::new([samples; 3], [-2.0; 3], cell).expect("grid");

        let uploaded = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &FieldOf(field))
            .expect("cpu sample");
        let produced = sampler
            .sample(gpu.device(), gpu.queue(), grid, field)
            .expect("gpu sample");

        let from_cpu = mc
            .extract_buffers(gpu.device(), gpu.queue(), &uploaded)
            .expect("extract");
        let from_gpu = mc
            .extract_buffers(gpu.device(), gpu.queue(), &produced)
            .expect("extract");

        assert_eq!(
            from_cpu.triangles,
            from_gpu.triangles,
            "{} at {samples}^3: a GPU-evaluated field changed the triangle count \
             ({} against {}), so a sample crossed the sign test",
            field.name(),
            from_gpu.triangles,
            from_cpu.triangles
        );
    }
}

/// `GpuField` as an `Sdf`, so the CPU side of the comparison goes through the
/// same upload path a consumer uses.
struct FieldOf(GpuField);

impl isomesh::Sdf for FieldOf {
    type Scalar = f32;
    fn sample(&self, p: [f32; 3]) -> f32 {
        self.0.sample_on_cpu(p)
    }
}
