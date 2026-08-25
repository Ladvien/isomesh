//! These need a real device, and they say so if they cannot get one.
//!
//! There is no CPU fallback: `headless::Gpu::new` refuses a software adapter on
//! purpose, so a machine with no GPU fails these rather than quietly measuring
//! lavapipe. That is the intended behaviour and it is why CI for this crate has
//! to name its runner.
//!
//! Exact float comparison throughout: upload and download are copies, and a
//! tolerance would hide a layout bug as readily as it would absorb one.
#![allow(clippy::float_cmp)]

use isomesh::fields::Sphere;
use isomesh::{Real, Sdf};

use super::{FieldBuffer, read_buffer, read_bytes_many_deferred};
use crate::headless::Gpu;
use crate::{Error, GridParams};

fn gpu() -> &'static Gpu {
    crate::headless::shared()
}

/// The round trip that everything above this depends on: CPU field in, GPU
/// memory, CPU floats out, **bit for bit**.
///
/// Bit-exact rather than approximate because nothing here is arithmetic --
/// upload and download are copies, and a tolerance would hide a layout bug as
/// easily as it would absorb one.
#[test]
fn a_sampled_field_survives_the_round_trip_bit_for_bit() {
    let gpu = gpu();
    let grid = GridParams::new([33; 3], [-2.0; 3], 0.125).expect("valid grid");
    let sphere = Sphere::<f32>::canonical();

    let field = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &sphere)
        .expect("upload a sampled field");
    let back = read_buffer(
        gpu.device(),
        gpu.queue(),
        field.buffer(),
        grid.field_buffer_size(),
    )
    .expect("read the samples back");

    assert_eq!(back.len() as u64, grid.sample_count());

    // Re-derive on the CPU in the documented index order and compare.
    let [sx, sy, sz] = grid.samples();
    let mut i = 0usize;
    for z in 0..sz {
        for y in 0..sy {
            for x in 0..sx {
                let expected = sphere.sample(grid.sample_position([x, y, z]));
                assert_eq!(
                    back[i].to_bits(),
                    expected.to_bits(),
                    "sample [{x}, {y}, {z}] came back changed"
                );
                i += 1;
            }
        }
    }
}

/// The index order is a convention, so it gets a test that would fail if the
/// loops were transposed -- which comparing against the same loops would not.
#[test]
fn x_varies_fastest() {
    let gpu = gpu();
    let grid = GridParams::new([4, 3, 2], [0.0; 3], 1.0).expect("valid grid");

    /// `f(p) = p.x + 10·p.y + 100·p.z`, so every sample names its own index.
    #[derive(Clone, Copy)]
    struct Ramp;
    impl Sdf for Ramp {
        type Scalar = f32;
        fn sample(&self, p: [f32; 3]) -> f32 {
            p[0] + 10.0 * p[1] + 100.0 * p[2]
        }
    }

    let field =
        FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &Ramp).expect("upload the ramp");
    let back = read_buffer(
        gpu.device(),
        gpu.queue(),
        field.buffer(),
        grid.field_buffer_size(),
    )
    .expect("read back");

    assert_eq!(back[0], 0.0);
    assert_eq!(back[1], 1.0, "index 1 must be x = 1");
    assert_eq!(back[4], 10.0, "index sx must be y = 1");
    assert_eq!(back[12], 100.0, "index sx*sy must be z = 1");
}

/// A wrong-length upload is rejected at the door rather than padded or
/// truncated: both repairs produce a buffer that looks meshable and describes a
/// different surface.
#[test]
fn a_mismatched_sample_count_is_refused_rather_than_repaired() {
    let gpu = gpu();
    let grid = GridParams::new([4, 4, 4], [0.0; 3], 1.0).expect("valid grid");
    let expected = grid.sample_count();

    for got in [expected - 1, expected + 1, 0] {
        let samples = vec![0.0f32; got as usize];
        assert_eq!(
            FieldBuffer::uploaded(gpu.device(), gpu.queue(), grid, &samples).err(),
            Some(Error::SampleCountMismatch { expected, got }),
            "accepted {got} samples for a {expected}-sample grid"
        );
    }
}

#[test]
fn an_unaligned_readback_is_refused() {
    let gpu = gpu();
    let grid = GridParams::new([4, 4, 4], [0.0; 3], 1.0).expect("valid grid");
    let field = FieldBuffer::new(gpu.device(), grid);

    assert_eq!(
        read_buffer(gpu.device(), gpu.queue(), field.buffer(), 6).err(),
        Some(Error::UnalignedReadback {
            bytes: 6,
            stride: 4
        })
    );
}

/// The two crates have to agree on the scalar as well as the layout.
#[test]
fn f32_is_what_crosses_the_boundary() {
    let gpu = gpu();
    let grid = GridParams::new([8; 3], [-1.0; 3], 0.25).expect("valid grid");
    let samples: Vec<f32> = (0..grid.sample_count())
        .map(|i| f32::from_f64(i as f64) * 0.5 - 3.0)
        .collect();

    let field = FieldBuffer::uploaded(gpu.device(), gpu.queue(), grid, &samples).expect("upload");
    let back = read_buffer(
        gpu.device(),
        gpu.queue(),
        field.buffer(),
        grid.field_buffer_size(),
    )
    .expect("read back");

    for (a, b) in samples.iter().zip(&back) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

/// The deferred read-back has to return exactly what the blocking one would,
/// **and** it has to split a batch back into the buffers that were asked for.
///
/// Two requests of different lengths, because one request cannot tell a correct
/// span split from a `to_vec()` of the whole mapping.
#[test]
fn a_deferred_readback_returns_the_bytes_it_was_asked_for() {
    let gpu = gpu();

    // Distinguishable patterns rather than counters from the same sequence: a
    // swapped pair of requests would still pass against `0, 1, 2, ...`.
    let first: Vec<u32> = (0..16u32).map(|i| 0xa000_0000 | i).collect();
    let second: Vec<u32> = (0..4u32).map(|i| 0x5b00_0000 | i).collect();

    let upload = |words: &[u32]| {
        let bytes: Vec<u8> = words.iter().flat_map(|w| w.to_le_bytes()).collect();
        let buffer = gpu.device().create_buffer(&wgpu::BufferDescriptor {
            label: Some("deferred readback source"),
            size: bytes.len() as u64,
            usage: wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        gpu.queue().write_buffer(&buffer, 0, &bytes);
        (buffer, bytes)
    };
    let (buffer_a, expected_a) = upload(&first);
    let (buffer_b, expected_b) = upload(&second);

    let readback = read_bytes_many_deferred(
        gpu.device(),
        gpu.queue(),
        &[
            (&buffer_a, expected_a.len() as u64),
            (&buffer_b, expected_b.len() as u64),
        ],
    )
    .expect("start the deferred readback");

    // The polling loop a frame-driven caller runs, bounded so a mapping that
    // never completes fails the test instead of hanging the suite.
    let mut spins = 0;
    while !readback.ready(gpu.device()) {
        spins += 1;
        assert!(spins < 5_000, "the mapping never completed");
        std::thread::sleep(std::time::Duration::from_millis(1));
    }

    let out = readback.take().expect("take the mapped bytes");
    assert_eq!(out.len(), 2, "one Vec per request, in request order");
    assert_eq!(out[0], expected_a);
    assert_eq!(out[1], expected_b);
}

/// Nothing requested is not an error, and it is not a zero-size buffer either --
/// wgpu rejects those, so the empty case must be represented without one.
#[test]
fn a_deferred_readback_of_nothing_is_ready_at_once() {
    let gpu = gpu();

    let readback =
        read_bytes_many_deferred(gpu.device(), gpu.queue(), &[]).expect("an empty request list");

    assert!(
        readback.ready(gpu.device()),
        "an empty readback has nothing to wait for"
    );
    assert!(readback.take().expect("take nothing").is_empty());
}

/// The same door check as the blocking sibling, refused before a buffer is
/// allocated rather than rounded up to an aligned size.
#[test]
fn an_unaligned_deferred_readback_is_refused() {
    let gpu = gpu();
    let grid = GridParams::new([4, 4, 4], [0.0; 3], 1.0).expect("valid grid");
    let field = FieldBuffer::new(gpu.device(), grid);

    assert_eq!(
        read_bytes_many_deferred(gpu.device(), gpu.queue(), &[(field.buffer(), 6)]).err(),
        Some(Error::UnalignedReadback {
            bytes: 6,
            stride: 4
        })
    );
}
