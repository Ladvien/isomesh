//! The scan against the CPU prefix sum it replaces, element for element.
//!
//! That comparison is the whole acceptance. A parallel scan that is subtly
//! wrong — one block's offset missing, a partial tail mis-summed — produces
//! triangles at the wrong offsets, which produces a mesh that renders and is
//! silently corrupt. There is no picture that catches it and no aggregate that
//! catches it; only the elements do.

use super::{PrefixScan, cpu_prefix_sum};
use crate::headless::Gpu;
use crate::read_buffer_u32;

fn gpu() -> &'static Gpu {
    crate::headless::shared()
}

/// Upload `counts` as a storage buffer.
fn upload(device: &wgpu::Device, queue: &wgpu::Queue, counts: &[u32]) -> wgpu::Buffer {
    let buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("scan test input"),
        size: (counts.len().max(1) * 4) as u64,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let mut bytes = Vec::with_capacity(counts.len() * 4);
    for n in counts {
        bytes.extend_from_slice(&n.to_le_bytes());
    }
    queue.write_buffer(&buffer, 0, &bytes);
    buffer
}

/// Scan `counts` on the GPU and assert it equals the CPU's answer exactly.
///
/// Returns the level count so a caller can assert which path it exercised.
fn check(gpu: &Gpu, scan: &PrefixScan, counts: &[u32], what: &str) -> usize {
    let input = upload(gpu.device(), gpu.queue(), counts);
    let out = scan
        .scan(gpu.device(), gpu.queue(), &input, counts.len() as u32)
        .expect("scan");

    let got = read_buffer_u32(
        gpu.device(),
        gpu.queue(),
        &out.offsets,
        (counts.len().max(1) * 4) as u64,
    )
    .expect("read offsets");

    let (expected, total) = cpu_prefix_sum(counts);
    assert_eq!(out.total, total, "{what}: grand total differs");
    for (i, (a, b)) in got.iter().zip(&expected).enumerate() {
        assert_eq!(
            a, b,
            "{what}: offset[{i}] is {a}, cpu says {b} (levels {})",
            out.levels
        );
    }
    out.levels
}

/// One block, no hierarchy at all.
#[test]
fn a_single_block_matches_the_cpu() {
    let gpu = gpu();
    let scan = PrefixScan::new(gpu.device()).expect("pipelines");
    let counts: Vec<u32> = (0..200).map(|i| (i % 7) as u32).collect();
    assert_eq!(check(gpu, &scan, &counts, "200 elements"), 1);
}

/// Exactly one full block, where the partial-tail handling is not exercised and
/// an off-by-one in the block total would show.
#[test]
fn exactly_one_full_block_matches_the_cpu() {
    let gpu = gpu();
    let scan = PrefixScan::new(gpu.device()).expect("pipelines");
    let counts: Vec<u32> = (0..PrefixScan::BLOCK).map(|i| i % 5).collect();
    assert_eq!(check(gpu, &scan, &counts, "256 elements"), 1);
}

/// Two levels, with a partial tail block — the ordinary case.
#[test]
fn two_levels_with_a_partial_tail_match_the_cpu() {
    let gpu = gpu();
    let scan = PrefixScan::new(gpu.device()).expect("pipelines");
    let counts: Vec<u32> = (0..1000).map(|i| (i % 11) as u32).collect();
    assert_eq!(check(gpu, &scan, &counts, "1000 elements"), 2);
}

/// Three levels, which is what a 129³ extraction actually uses.
///
/// The level that only exists here is where a cross-block offset can go
/// missing, so a suite that stopped at two levels would not be testing the
/// path that ships.
#[test]
fn three_levels_match_the_cpu() {
    let gpu = gpu();
    let scan = PrefixScan::new(gpu.device()).expect("pipelines");
    // 300_000 -> 1172 blocks -> 5 blocks -> 1.
    let counts: Vec<u32> = (0..300_000).map(|i| (i % 13) as u32).collect();
    assert_eq!(check(gpu, &scan, &counts, "300k elements"), 3);
}

/// The degenerate inputs, which is where a scan usually breaks.
#[test]
fn degenerate_inputs_match_the_cpu() {
    let gpu = gpu();
    let scan = PrefixScan::new(gpu.device()).expect("pipelines");

    check(gpu, &scan, &[0], "one zero");
    check(gpu, &scan, &[7], "one element");
    check(gpu, &scan, &vec![0u32; 1000], "all zeros");
    check(gpu, &scan, &vec![1u32; 1000], "all ones");
    // A block boundary landing exactly on the end, one past, and one short.
    for n in [
        PrefixScan::BLOCK - 1,
        PrefixScan::BLOCK + 1,
        PrefixScan::BLOCK * 2,
        PrefixScan::BLOCK * PrefixScan::BLOCK,
    ] {
        let counts: Vec<u32> = (0..n).map(|i| i % 3).collect();
        check(gpu, &scan, &counts, &format!("{n} elements"));
    }
}

/// Marching Cubes counts are sparse and clustered — most cells are empty and
/// the surface cells are adjacent. A scan tested only on smooth data has not
/// met its actual input.
#[test]
fn a_realistic_sparse_distribution_matches_the_cpu() {
    let gpu = gpu();
    let scan = PrefixScan::new(gpu.device()).expect("pipelines");

    // A shell: zero almost everywhere, a run of non-zero where a surface would
    // cross, exactly as `count_cells` produces.
    let counts: Vec<u32> = (0..50_000u32)
        .map(|i| {
            let on_surface = (i % 997) < 12;
            if on_surface { (i % 5) + 1 } else { 0 }
        })
        .collect();
    check(gpu, &scan, &counts, "sparse shell");
}

/// The CPU reference is only worth comparing against if it is obviously right.
#[test]
fn the_cpu_reference_is_an_exclusive_scan() {
    let (offsets, total) = cpu_prefix_sum(&[3, 0, 4, 1]);
    assert_eq!(offsets, vec![0, 3, 3, 7]);
    assert_eq!(total, 8);

    let (empty, zero) = cpu_prefix_sum(&[]);
    assert!(empty.is_empty());
    assert_eq!(zero, 0);
}

/// A scan of zero elements is the empty scan, not a scan of one stale word.
///
/// The nonzero word is planted deliberately: scanned with `n = 0`, it must not
/// come back as the total the way `n.max(1)` once made it.
#[test]
fn scanning_zero_elements_matches_the_cpu_and_reads_nothing() {
    let gpu = gpu();
    let scan = PrefixScan::new(gpu.device()).expect("pipelines");
    let input = upload(gpu.device(), gpu.queue(), &[0xdead_beef]);
    let out = scan
        .scan(gpu.device(), gpu.queue(), &input, 0)
        .expect("scan");
    let (expected, total) = cpu_prefix_sum(&[]);
    assert!(expected.is_empty());
    assert_eq!(out.total, total, "empty scan: grand total differs");
    assert_eq!(out.levels, 0, "no elements need no hierarchy");
}

/// More elements than one dispatch covers is a named error, not a wgpu panic.
#[test]
fn a_scan_wider_than_the_dispatch_limit_is_refused() {
    let gpu = gpu();
    let scan = PrefixScan::new(gpu.device()).expect("pipelines");
    let limit = gpu.device().limits().max_compute_workgroups_per_dimension;
    let Some(n) = limit
        .checked_mul(PrefixScan::BLOCK)
        .and_then(|x| x.checked_add(1))
    else {
        return; // a device this wide cannot express the failing input in u32
    };
    // The guard fires before `counts` is ever bound, so a one-word buffer with
    // a lying `n` is the honest probe -- without the guard this exact call is
    // a wgpu validation panic, not an error value.
    let input = upload(gpu.device(), gpu.queue(), &[0]);
    let err = scan.scan(gpu.device(), gpu.queue(), &input, n);
    assert!(
        matches!(err, Err(crate::Error::ScanTooLong { .. })),
        "expected ScanTooLong, got {err:?}"
    );
}
