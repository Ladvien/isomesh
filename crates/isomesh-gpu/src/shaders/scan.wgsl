// Exclusive prefix sum over the per-cell triangle counts, on the GPU.
//
// This replaces the largest avoidable cost in the extraction path. M-149 broke
// a 129^3 run into upload 8.61 ms, CPU prefix sum 3.28, counts read-back 1.97,
// geometry read-back 1.04 -- and the middle two exist only because the scan was
// done on the CPU, which meant copying 4 bytes per cell (8 MB at 129^3) home to
// add them up. Scanning here leaves the CPU needing four bytes: the grand
// total, to size the output buffer.
//
// # The algorithm, and why this one
//
// A hierarchical scan: each workgroup scans its own block in workgroup memory
// and publishes the block's total; the block totals are themselves scanned, one
// level up, until a level fits in a single workgroup; then each level adds its
// parent's exclusive scan back down. The level loop lives on the CPU because
// the number of levels depends on the cell count -- three for 129^3.
//
// The in-block scan is Hillis-Steele. It does O(n log n) work where a Blelloch
// up-down sweep does O(n), and it is chosen anyway: n here is 256, the whole
// thing is in workgroup memory, and the correctness argument is one loop rather
// than two sweeps with a mid-point swap. The measured cost is what decides
// whether that was right -- see the acceptance in GPU-010a -- and a scan that
// is subtly wrong produces a mesh that looks entirely plausible, which is the
// failure this trades against.

// Elements scanned per workgroup. Must equal the `@workgroup_size` below and
// `PrefixScan::BLOCK` on the CPU side.
const BLOCK: u32 = 256u;

struct ScanParams {
    // Elements at this level. The last block is partial whenever this is not a
    // multiple of BLOCK, which is the usual case.
    n: u32,
}

@group(0) @binding(0) var<uniform> params: ScanParams;
@group(0) @binding(1) var<storage, read> input: array<u32>;
@group(0) @binding(2) var<storage, read_write> output: array<u32>;
@group(0) @binding(3) var<storage, read_write> block_sums: array<u32>;

var<workgroup> temp: array<u32, BLOCK>;

// Scan one block, write its exclusive scan to `output`, its total to
// `block_sums[workgroup]`.
@compute @workgroup_size(256)
fn scan_blocks(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_index) tid: u32,
    @builtin(workgroup_id) wg: vec3<u32>,
) {
    let i = gid.x;
    // Threads past the end contribute zero, which keeps the block total right
    // without a separate partial-block path.
    var value = 0u;
    if (i < params.n) {
        value = input[i];
    }
    temp[tid] = value;
    workgroupBarrier();

    // Hillis-Steele inclusive scan. The two barriers are not decoration: every
    // read of `temp[tid - offset]` must complete before any write of
    // `temp[tid]` in the same round, or a thread reads a value from the round
    // it is currently in. Both barriers sit in uniform control flow -- the loop
    // bound is uniform and the `if`s contain no barrier -- which WGSL requires.
    for (var offset = 1u; offset < BLOCK; offset = offset << 1u) {
        var addend = 0u;
        if (tid >= offset) {
            addend = temp[tid - offset];
        }
        workgroupBarrier();
        if (tid >= offset) {
            temp[tid] = temp[tid] + addend;
        }
        workgroupBarrier();
    }

    // Exclusive from inclusive: subtract this element's own contribution.
    if (i < params.n) {
        output[i] = temp[tid] - value;
    }
    // The last lane holds the inclusive total of the whole block, including the
    // zeros any out-of-range lanes contributed.
    if (tid == BLOCK - 1u) {
        block_sums[wg.x] = temp[tid];
    }
}

// Add the parent level's exclusive scan back into this level's blocks.
//
// `block_sums` here is the *scanned* totals, so block 0 adds zero and the rest
// add everything before them.
@compute @workgroup_size(256)
fn add_block_offsets(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(workgroup_id) wg: vec3<u32>,
) {
    let i = gid.x;
    if (i >= params.n) {
        return;
    }
    output[i] = output[i] + block_sums[wg.x];
}
