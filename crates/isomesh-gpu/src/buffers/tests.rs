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

use super::{DeferredGeometry, FieldBuffer, read_buffer, read_bytes_many_deferred};
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

/// A zero-capacity queue is refused at construction rather than becoming a
/// scheduler that meshes nothing.
///
/// Every `submit` into it would fail, so the only thing a `new(0)` can produce
/// is a caller who thinks they have a queue.
#[test]
fn a_deferred_queue_of_no_capacity_is_refused() {
    assert_eq!(
        DeferredGeometry::<u32>::new(0).err(),
        Some(Error::DeferredQueueFull { capacity: 0 })
    );
}

/// The queue is a **budget**, so it has to refuse rather than grow, and the
/// refusal must not consume a slot.
///
/// `in_flight` before and after the refused submit is the assertion that
/// matters: a queue that counted the failure would report itself fuller than it
/// is and starve the scheduler by one chunk a frame, forever.
#[test]
fn a_full_deferred_queue_refuses_and_stays_where_it_was() {
    let gpu = gpu();
    let mut queue = DeferredGeometry::new(2).expect("a queue of two");

    for key in 0u32..2 {
        assert!(queue.has_room(), "room for {key}");
        let readback =
            read_bytes_many_deferred(gpu.device(), gpu.queue(), &[]).expect("an empty readback");
        queue.submit(key, readback).expect("submit within capacity");
    }

    assert!(!queue.has_room());
    assert_eq!(queue.in_flight(), 2);
    assert_eq!(queue.capacity(), 2);

    let extra =
        read_bytes_many_deferred(gpu.device(), gpu.queue(), &[]).expect("an empty readback");
    assert_eq!(
        queue.submit(99, extra).err(),
        Some(Error::DeferredQueueFull { capacity: 2 })
    );
    assert_eq!(queue.in_flight(), 2, "a refused submit consumes no slot");
}

/// The whole point: bytes out, **under the key they went in with**, and the slot
/// freed.
///
/// Geometry that comes back without saying which chunk it belongs to is geometry
/// a caller cannot install, so the key travelling with the bytes is the contract
/// and not a convenience. Eight bytes of a recognisable pattern rather than
/// zeros, because a queue that returned a fresh allocation would pass against
/// zeros.
#[test]
fn a_deferred_queue_returns_the_bytes_under_their_key() {
    let gpu = gpu();
    let expected: Vec<u8> = 0xdead_beef_0bad_f00du64.to_le_bytes().to_vec();
    let source = gpu.device().create_buffer(&wgpu::BufferDescriptor {
        label: Some("deferred queue source"),
        size: expected.len() as u64,
        usage: wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    gpu.queue().write_buffer(&source, 0, &expected);

    let mut queue = DeferredGeometry::new(1).expect("a queue of one");
    let readback = read_bytes_many_deferred(
        gpu.device(),
        gpu.queue(),
        &[(&source, expected.len() as u64)],
    )
    .expect("start the deferred readback");
    queue.submit("chunk", readback).expect("submit");
    assert_eq!(queue.in_flight(), 1);

    // The frame loop a scheduler runs, bounded so a mapping that never completes
    // fails rather than hangs.
    let mut collected = Vec::new();
    let mut frames = 0;
    while collected.is_empty() {
        collected = queue.drain_ready(gpu.device()).expect("drain");
        frames += 1;
        assert!(frames < 5_000, "the mapping never completed");
        std::thread::sleep(std::time::Duration::from_millis(1));
    }

    assert_eq!(collected.len(), 1);
    assert_eq!(collected[0].0, "chunk", "the key it was submitted under");
    assert_eq!(collected[0].1, vec![expected], "one Vec per request");
    assert_eq!(queue.in_flight(), 0, "the slot is free again");
    assert!(queue.has_room());
}

/// Draining an empty queue is a no-op that returns an empty harvest, which is
/// what makes `drain_ready` safe to call unconditionally once a frame.
///
/// It never reaches [`super::Readback::ready`] because there is no slot to ask,
/// so it does not poll the device either — a scheduler with nothing in flight
/// pays nothing for asking. Called twice: an empty `Vec` must mean "nothing yet"
/// on every call, never "nothing ever".
#[test]
fn draining_an_empty_deferred_queue_harvests_nothing() {
    let gpu = gpu();
    let mut queue = DeferredGeometry::<u32>::new(4).expect("a queue of four");

    for _ in 0..2 {
        assert!(queue.drain_ready(gpu.device()).expect("drain").is_empty());
        assert_eq!(queue.in_flight(), 0);
        assert!(queue.has_room());
    }
}
