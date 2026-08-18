use super::*;
use crate::RuntimeShape3;

const CELLS: u32 = 4;

fn layout() -> ChunkLayout<f64> {
    ChunkLayout::new(CELLS, 1.0, [0.0; 3]).expect("valid layout")
}

/// A deterministic pseudo-random phase field over **global** samples, so the
/// world and the single grid are looking at the same thing by construction.
///
/// Derived from the coordinate rather than from an iteration order: the two
/// fixtures have to agree bit for bit, and a generator advanced in loop order
/// would not.
fn phase(g: [i64; 3]) -> f64 {
    let mut h = 0x9e37_79b9_7f4a_7c15u64;
    for c in g {
        h ^= (c as u64).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        h = h.wrapping_mul(0x94d0_49bb_1331_11eb);
        h ^= h >> 31;
    }
    // Roughly half air, which is where a lattice has the most components.
    if h % 100 < 45 { 1.0 } else { -1.0 }
}

/// `chunks` chunks side by side in x, plus the single grid spanning all of them.
fn world_and_grid(chunks: i32) -> (AirWorld<f64>, Air, [i64; 3]) {
    let l = layout();
    let mut world = AirWorld::new(l);
    let n = CELLS as usize + 1;

    for c in 0..chunks {
        let id = ChunkId::new([c, 0, 0]);
        let mut values = Vec::with_capacity(n * n * n);
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let g = l.global_sample(id, [x as u32, y as u32, z as u32]);
                    values.push(phase(g));
                }
            }
        }
        world.load(id, &values).expect("load");
    }

    // The same samples as one grid. Chunks overlap by a plane, so the span is
    // `chunks * cells + 1` rather than `chunks * (cells + 1)`.
    let span = chunks as usize * CELLS as usize + 1;
    let shape = RuntimeShape3::new([span as u32, n as u32, n as u32]).expect("shape");
    let mut values = Vec::with_capacity(span * n * n);
    for z in 0..n {
        for y in 0..n {
            for x in 0..span {
                values.push(phase([x as i64, y as i64, z as i64]));
            }
        }
    }
    let (grid, _) = Air::build(&values, &shape).expect("build");
    // The addressable span, which is NOT the sample span. `local_sample` gives
    // the overlap plane to the *next* chunk, so the plane at `chunks * cells`
    // belongs to a chunk nobody loaded — and likewise the top plane in y and z.
    // The chunks' arrays still cover those samples and label them; they simply
    // cannot be named from the world. Comparing them would test the fixture's
    // arithmetic rather than the stitch.
    let owned = [chunks as i64 * CELLS as i64, CELLS as i64, CELLS as i64];
    (world, grid, owned)
}

/// **R-028's falsifier: a stitched world must answer exactly as one grid does.**
///
/// A cost measurement cannot see a structure that is fast and wrong — the same
/// argument P-26 made, and M-321 proved it the hard way when an `O(n³)` defect
/// passed all fifteen correctness tests. So membership is compared **pairwise**,
/// not merely by component count: two wrong stitchings can agree on a total.
#[test]
fn a_stitched_world_agrees_with_one_grid() {
    let (world, grid, dims) = world_and_grid(3);
    assert_eq!(world.loaded(), 3);

    assert_eq!(
        world.components(),
        grid.components(),
        "world and single grid disagree on how many components exist"
    );

    let mut air_seen = 0;
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let g = [x, y, z];
                let in_world = world.component_at(g).is_some();
                let in_grid = grid.label_of([x as u32, y as u32, z as u32]).is_some();
                assert_eq!(in_world, in_grid, "air disagrees at {g:?}");
                if in_world {
                    air_seen += 1;
                }
            }
        }
    }
    assert!(air_seen > 0, "the fixture has no air to compare");

    // Pairwise membership. This is what catches a seam joined in the wrong
    // direction, which every count-based check passes.
    for z in 0..dims[2] {
        for x in 0..dims[0] {
            for x2 in 0..dims[0] {
                let (p, q) = ([x, 1, z], [x2, 2, z]);
                let w = world.connected(p, q);
                let s = grid.connected(
                    [p[0] as u32, p[1] as u32, p[2] as u32],
                    [q[0] as u32, q[1] as u32, q[2] as u32],
                );
                assert_eq!(w, s, "connected({p:?}, {q:?}) diverged");
            }
        }
    }
}

/// **A tunnel crossing a seam is one component, and the seam is what makes it
/// so.**
#[test]
fn a_tunnel_across_a_seam_is_one_component() {
    let l = layout();
    let mut world = AirWorld::new(l);
    let n = CELLS as usize + 1;
    let solid = alloc::vec![-1.0_f64; n * n * n];
    for c in 0..2 {
        world.load(ChunkId::new([c, 0, 0]), &solid).expect("load");
    }
    assert_eq!(world.components(), 0);

    // A line along x through both chunks, meeting on the shared plane: chunk 0's
    // local `cells` is chunk 1's local 0.
    let line: Vec<[u32; 3]> = (0..=CELLS).map(|x| [x, 2, 2]).collect();
    world.dig(ChunkId::new([0, 0, 0]), &line, || true);
    assert_eq!(world.components(), 1, "one chunk dug");
    world.dig(ChunkId::new([1, 0, 0]), &line, || true);

    assert_eq!(world.components(), 1, "the seam joins them into one");
    assert!(
        world.connected([0, 2, 2], [2 * i64::from(CELLS) - 1, 2, 2]),
        "opposite ends of the world are connected through the seam"
    );
}

/// **Two chunks with air that does not meet at the seam stay separate.**
///
/// The complement of the test above, and the one that fails if the stitch joins
/// labels that merely appear on a shared plane rather than at the same sample.
#[test]
fn air_that_does_not_meet_at_the_seam_stays_separate() {
    let l = layout();
    let mut world = AirWorld::new(l);
    let n = CELLS as usize + 1;
    let solid = alloc::vec![-1.0_f64; n * n * n];
    for c in 0..2 {
        world.load(ChunkId::new([c, 0, 0]), &solid).expect("load");
    }

    // Both tunnels reach the shared plane, at different heights.
    let left: Vec<[u32; 3]> = (0..=CELLS).map(|x| [x, 1, 1]).collect();
    let right: Vec<[u32; 3]> = (0..=CELLS).map(|x| [x, 3, 3]).collect();
    world.dig(ChunkId::new([0, 0, 0]), &left, || true);
    world.dig(ChunkId::new([1, 0, 0]), &right, || true);

    assert_eq!(
        world.components(),
        2,
        "they touch the plane, not each other"
    );
    assert!(!world.connected([0, 1, 1], [2 * i64::from(CELLS) - 1, 3, 3]));
}

/// **Filling at a seam severs a component that spans two chunks.**
///
/// The sealed-volume mechanic across a chunk boundary — the case a per-chunk
/// structure with no stitching cannot even represent.
#[test]
fn filling_at_a_seam_severs_across_chunks() {
    let l = layout();
    let mut world = AirWorld::new(l);
    let n = CELLS as usize + 1;
    let solid = alloc::vec![-1.0_f64; n * n * n];
    for c in 0..2 {
        world.load(ChunkId::new([c, 0, 0]), &solid).expect("load");
    }
    let line: Vec<[u32; 3]> = (0..=CELLS).map(|x| [x, 2, 2]).collect();
    world.dig(ChunkId::new([0, 0, 0]), &line, || true);
    world.dig(ChunkId::new([1, 0, 0]), &line, || true);
    assert_eq!(world.components(), 1);

    // Cut the shared plane. It is owned by chunk 1 at local 0 and chunk 0 holds
    // a copy at local `cells`; both have to go, or the seam still sees air on
    // one side. That duplication is the price of the overlap and is why this
    // test exists rather than being obvious.
    world.fill(ChunkId::new([0, 0, 0]), &[[CELLS, 2, 2]], || true);
    world.fill(ChunkId::new([1, 0, 0]), &[[0, 2, 2]], || true);

    assert_eq!(world.components(), 2, "the passage is sealed");
    assert!(!world.connected([0, 2, 2], [2 * i64::from(CELLS) - 1, 2, 2]));
}

/// **The bisect search is bounded by the chunk, not by the world.**
///
/// R-028's hypothesis, and the reason the ticket exists. M-321 measured a tunnel
/// bisection at 1.1× a full rebuild on one large grid; here the same edit cannot
/// visit more than one chunk holds, whatever the world's size.
#[test]
fn a_bisect_visits_no_more_than_one_chunk() {
    let l = layout();
    let n = CELLS as usize + 1;
    let per_chunk = (n * n * n) as u64;

    let mut visited_by_width = Vec::new();
    for width in [2i32, 4, 8] {
        let mut world = AirWorld::new(l);
        let solid = alloc::vec![-1.0_f64; n * n * n];
        for c in 0..width {
            world.load(ChunkId::new([c, 0, 0]), &solid).expect("load");
        }
        // One tunnel spanning every chunk.
        let line: Vec<[u32; 3]> = (0..=CELLS).map(|x| [x, 2, 2]).collect();
        for c in 0..width {
            world.dig(ChunkId::new([c, 0, 0]), &line, || true);
        }
        assert_eq!(world.components(), 1, "width {width}: one tunnel");

        // Bisect inside the middle chunk.
        let mid = ChunkId::new([width / 2, 0, 0]);
        let f = world.fill(mid, &[[2, 2, 2]], || true).expect("loaded");
        assert!(
            f.visited <= per_chunk,
            "width {width}: visited {} exceeds one chunk's {per_chunk} samples",
            f.visited
        );
        visited_by_width.push(f.visited);
    }

    // The bound is the point: doubling the world twice must not move it.
    let first = visited_by_width.first().copied().unwrap_or(0);
    assert!(
        visited_by_width.iter().all(|v| *v == first),
        "the search cost tracked the world rather than the chunk: {visited_by_width:?}"
    );
}
