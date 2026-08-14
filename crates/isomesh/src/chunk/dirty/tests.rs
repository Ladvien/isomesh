//! The load-bearing output here is E1 — the fraction of a brush's bounding box
//! that genuinely changes — which the research records as unpublished.

use alloc::vec::Vec;

use super::{BudgetReport, DirtySet, EditReport, mark_edit};
use crate::chunk::{ChunkId, ChunkLayout};
use crate::fields::{BoxExact, Difference, Sphere};

fn layout() -> ChunkLayout<f64> {
    // Power-of-two spacing, per M-32: chunk seams are only bit-exact there.
    ChunkLayout::new(16, 0.0625, [-4.0; 3]).expect("valid layout")
}

/// A brush: the solid with a sphere subtracted from it.
fn carved(radius: f64, centre: [f64; 3]) -> Difference<BoxExact<f64>, Sphere<f64>> {
    Difference {
        a: BoxExact::canonical(),
        b: Sphere {
            center: centre,
            radius,
        },
    }
}

/// The cell range a world-space sphere touches — what a brush marks, knowing
/// nothing about what actually changed inside it.
fn brush_region(l: &ChunkLayout<f64>, centre: [f64; 3], radius: f64) -> ([i64; 3], [i64; 3]) {
    let h = l.cell_size();
    let mut min = [0i64; 3];
    let mut max = [0i64; 3];
    for axis in 0..3 {
        min[axis] = ((centre[axis] - radius + 4.0) / h).floor() as i64;
        max[axis] = ((centre[axis] + radius + 4.0) / h).ceil() as i64;
    }
    (min, max)
}

// ─── the dirty set ──────────────────────────────────────────────────────────

#[test]
fn the_set_is_sorted_deduplicated_and_order_independent() {
    let mut a = DirtySet::new();
    let mut b = DirtySet::new();
    let ids = [
        ChunkId::new([2, 0, 0]),
        ChunkId::new([-1, 3, 1]),
        ChunkId::new([0, 0, 0]),
        ChunkId::new([2, 0, 0]),
    ];
    for id in ids {
        a.insert(id);
    }
    for id in ids.iter().rev() {
        b.insert(*id);
    }
    assert_eq!(a, b, "insertion order must not survive into the queue");
    assert_eq!(a.len(), 3, "duplicates must collapse");

    let seen: Vec<ChunkId> = a.iter().collect();
    let mut sorted = seen.clone();
    sorted.sort();
    assert_eq!(seen, sorted, "iteration must be ascending");
}

#[test]
fn mesh_dirty_visits_each_chunk_once_and_clears() {
    let l = layout();
    let mut dirty = DirtySet::new();
    for id in [ChunkId::new([0, 0, 0]), ChunkId::new([1, 0, 0])] {
        dirty.insert(id);
        dirty.insert(id);
    }

    let mut visited = Vec::new();
    let done = dirty.mesh_dirty(&l, |id, origin| {
        // The origin is handed in rather than computed by the caller, so a
        // consumer cannot get the seam arithmetic wrong.
        let want = l.sample_origin(id);
        assert!(
            (0..3).all(|a| origin[a].to_bits() == want[a].to_bits()),
            "origin {origin:?} is not the layout\'s {want:?}"
        );
        visited.push(id);
    });

    assert_eq!(done, 2);
    assert_eq!(visited.len(), 2);
    assert!(dirty.is_empty(), "the queue must be cleared after meshing");
}

#[test]
fn an_inverted_region_changes_nothing() {
    let l = layout();
    let mut dirty = DirtySet::new();
    let field = BoxExact::<f64>::canonical();
    let report = mark_edit(&l, &field, &field, [5, 5, 5], [4, 5, 5], &mut dirty)
        .expect("an empty region is in bounds");
    assert_eq!(report, EditReport::default());
    assert!(dirty.is_empty());
}

/// An edit that changes nothing must dirty nothing, or every incremental scheme
/// built on this degenerates into re-meshing the world.
#[test]
fn an_identical_field_dirties_no_chunks() {
    let l = layout();
    let mut dirty = DirtySet::new();
    let field = carved(0.5, [0.0, 0.0, 0.0]);
    let (min, max) = brush_region(&l, [0.0, 0.0, 0.0], 0.6);
    let report = mark_edit(&l, &field, &field, min, max, &mut dirty).expect("a region in bounds");

    assert!(report.region_cells > 0, "the region must be non-trivial");
    assert_eq!(report.value_changed_cells, 0);
    assert_eq!(report.output_changed_cells, 0);
    assert_eq!(report.sign_changed_cells, 0);
    assert_eq!(report.swept_cells, 0);
    assert!(dirty.is_empty());
}

// ─── E1 ─────────────────────────────────────────────────────────────────────

/// **E1, the number the research records as unpublished.**
///
/// A sphere brush is carved two cells deeper, and this reports how much of the
/// brush's own bounding box actually changed. Reported rather than gated: the
/// answer depends on the brush and the grid, and the point is to know it.
///
/// Two numbers, because they answer different questions. **Value-changed** cells
/// genuinely need re-meshing, since a vertex interpolates the samples and moves
/// when they move. **Sign-changed** cells are where the triangles differ rather
/// than merely shifting.
#[test]
fn e1_the_fraction_of_a_brush_that_actually_changes() {
    let l = layout();
    std::println!(
        "measured: G-002 E1 -- fraction of a brush's bounding box that actually changes (h = {}, {} cells/chunk)",
        l.cell_size(),
        l.cells()
    );
    std::println!(
        "  {:>7} {:>9} {:>11} {:>10} {:>10} {:>8} {:>14}",
        "radius",
        "region",
        "value-chg",
        "output-chg",
        "sign-chg",
        "E1",
        "chunks dirty"
    );

    let centre = [0.1, 0.2, -0.15];
    for radius in [0.25f64, 0.5, 0.75, 1.0] {
        let grown = radius + l.cell_size() * 2.0;
        let before = carved(radius, centre);
        let after = carved(grown, centre);
        let (min, max) = brush_region(&l, centre, grown);

        let mut dirty = DirtySet::new();
        let r = mark_edit(&l, &before, &after, min, max, &mut dirty).expect("a region in bounds");

        std::println!(
            "  {radius:>7.2} {:>9} {:>11} {:>10} {:>9} {:>7} {:>7.1}% {:>5}/{:<7}",
            r.region_cells,
            r.value_changed_cells,
            r.output_changed_cells,
            r.sign_changed_cells,
            r.swept_cells,
            r.changed_fraction() * 100.0,
            r.dirty_chunks,
            r.region_chunks
        );

        assert!(r.output_changed_cells > 0, "the edit changed nothing");
        assert!(
            r.output_changed_cells <= r.value_changed_cells,
            "an output change implies a value change"
        );
        // Deliberately NOT `sign <= output`: a cell the surface swept entirely
        // through has every corner flipped and no triangles either side.
        assert!(
            r.swept_cells <= r.sign_changed_cells,
            "a swept cell is a sign-changed cell"
        );
        assert!(
            r.value_changed_cells <= r.region_cells,
            "more cells changed than the region contains"
        );
        assert_eq!(
            dirty.len() as u64,
            r.dirty_chunks,
            "the queue and the report disagree about how many chunks are dirty"
        );
    }
}

/// The same edit at three chunk sizes.
///
/// One changed cell dirties its whole chunk, so the *cell* fraction is a property
/// of the edit while the *chunk* fraction is a property of the chunk size. The
/// gap between them is what a finer dirty granularity would buy, measured rather
/// than argued.
#[test]
fn chunk_size_decides_how_much_the_dirty_set_over_marks() {
    std::println!("measured: G-002 granularity -- one edit, three chunk sizes");
    let centre = [0.1, 0.2, -0.15];
    for cells in [8u32, 16, 32] {
        let l = ChunkLayout::new(cells, 0.0625, [-4.0; 3]).expect("valid layout");
        let before = carved(0.5, centre);
        let after = carved(0.5 + l.cell_size() * 2.0, centre);
        let (min, max) = brush_region(&l, centre, 0.63);

        let mut dirty = DirtySet::new();
        let r = mark_edit(&l, &before, &after, min, max, &mut dirty).expect("a region in bounds");
        std::println!(
            "  {cells:>3} cells/chunk -> E1 {:>5.1}% of cells, {:>5.1}% of chunks ({}/{})",
            r.changed_fraction() * 100.0,
            r.dirty_chunk_fraction() * 100.0,
            r.dirty_chunks,
            r.region_chunks
        );
        assert!(
            r.dirty_chunk_fraction() >= r.changed_fraction() - 1e-9,
            "chunk granularity can only over-mark, never under-mark"
        );
    }
}

/// The dirty set must name chunks that really overlap the edit, or `mesh_dirty`
/// re-meshes the wrong ones and the seams drift apart.
#[test]
fn dirtied_chunks_overlap_the_region() {
    let l = layout();
    let centre = [0.1, 0.2, -0.15];
    let before = carved(0.5, centre);
    let after = carved(0.5 + l.cell_size() * 2.0, centre);
    let (min, max) = brush_region(&l, centre, 0.63);

    let mut dirty = DirtySet::new();
    let report = mark_edit(&l, &before, &after, min, max, &mut dirty).expect("a region in bounds");
    assert!(!dirty.is_empty());

    for id in dirty.iter() {
        let base = l.base_sample(id);
        let overlaps = (0..3).all(|a| base[a] + i64::from(l.cells()) > min[a] && base[a] <= max[a]);
        assert!(
            overlaps,
            "chunk {id:?} was dirtied but does not overlap the region"
        );
    }
    assert_eq!(dirty.len() as u64, report.dirty_chunks);
}

// ─── G-006: the frame budget ────────────────────────────────────────────────

/// Nearest first, and the order is a pure function of the contents and the
/// camera — never of the order edits arrived.
#[test]
fn a_budgeted_pass_meshes_nearest_first() {
    let layout = ChunkLayout::<f64>::new(8, 0.125, [0.0; 3]).expect("valid layout");
    let mut set = DirtySet::new();
    // Inserted deliberately out of distance order.
    for coords in [[3, 0, 0], [0, 0, 0], [2, 0, 0], [1, 0, 0]] {
        set.insert(ChunkId::new(coords));
    }

    let mut visited = Vec::new();
    let report = set.mesh_within_budget(
        &layout,
        [0.0; 3],
        |id, _| visited.push(id.coords[0]),
        || true,
    );
    assert_eq!(visited, alloc::vec![0, 1, 2, 3]);
    assert_eq!(report.meshed, 4);
    assert!(report.is_drained());

    // From the far side, the order reverses.
    let mut set = DirtySet::new();
    for coords in [[3, 0, 0], [0, 0, 0], [2, 0, 0], [1, 0, 0]] {
        set.insert(ChunkId::new(coords));
    }
    let mut visited = Vec::new();
    set.mesh_within_budget(
        &layout,
        [10.0, 0.0, 0.0],
        |id, _| visited.push(id.coords[0]),
        || true,
    );
    assert_eq!(visited, alloc::vec![3, 2, 1, 0]);
}

/// The budget stops the pass and the remainder survives to the next call, with
/// nothing lost and nothing done twice.
#[test]
fn an_exhausted_budget_leaves_the_rest_for_next_time() {
    let layout = ChunkLayout::<f64>::new(8, 0.125, [0.0; 3]).expect("valid layout");
    let mut set = DirtySet::new();
    for x in 0..10i32 {
        set.insert(ChunkId::new([x, 0, 0]));
    }

    let mut all = Vec::new();
    let mut passes = 0;
    while !set.is_empty() {
        let mut this_pass = 0;
        let report = set.mesh_within_budget(
            &layout,
            [0.0; 3],
            |id, _| all.push(id.coords[0]),
            || {
                this_pass += 1;
                this_pass < 3
            },
        );
        assert!(report.meshed > 0, "a pass must always make progress");
        passes += 1;
        assert!(passes <= 10, "the queue is not draining");
    }

    // Every chunk exactly once, in distance order across the passes.
    assert_eq!(all, (0..10).collect::<Vec<i32>>());
}

/// **A budget too small for one chunk must still make progress.**
///
/// Consulting `spend` *before* a chunk instead of after would drain nothing,
/// forever, while the queue grew — a livelock that presents as a memory leak.
/// Overshooting by at most one chunk is the price and it is the right way round.
#[test]
fn a_budget_that_is_already_gone_still_meshes_one_chunk() {
    let layout = ChunkLayout::<f64>::new(8, 0.125, [0.0; 3]).expect("valid layout");
    let mut set = DirtySet::new();
    for x in 0..4i32 {
        set.insert(ChunkId::new([x, 0, 0]));
    }

    let mut passes = 0;
    while !set.is_empty() {
        let report = set.mesh_within_budget(&layout, [0.0; 3], |_, _| {}, || false);
        assert_eq!(report.meshed, 1, "exactly one, never zero");
        passes += 1;
        assert!(passes <= 4, "the queue is not draining");
    }
    assert_eq!(passes, 4);
}

/// Ties break by `ChunkId`, so equidistant chunks have a defined order rather
/// than whatever the sort happened to do.
#[test]
fn equidistant_chunks_break_ties_deterministically() {
    let layout = ChunkLayout::<f64>::new(8, 0.125, [0.0; 3]).expect("valid layout");
    // Six chunks the same distance from the origin-centred camera.
    let ring = [
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1],
    ];

    let mut orders = Vec::new();
    for start in 0..ring.len() {
        let mut set = DirtySet::new();
        // Insert in a rotated order each time.
        for k in 0..ring.len() {
            set.insert(ChunkId::new(ring[(start + k) % ring.len()]));
        }
        let centre = layout.cell_size() * 4.0;
        let mut visited = Vec::new();
        set.mesh_within_budget(
            &layout,
            [centre, centre, centre],
            |id, _| visited.push(id.coords),
            || true,
        );
        orders.push(visited);
    }
    for window in orders.windows(2) {
        assert_eq!(
            window[0], window[1],
            "insertion order leaked into the schedule"
        );
    }
}

/// The origin handed to the mesh callback is the layout's, not one the caller
/// recomputed — the same guarantee `mesh_dirty` gives, and for M-32's reason.
#[test]
fn the_budgeted_pass_supplies_the_layouts_own_origin() {
    // Exact: the whole point is that the caller never does this arithmetic.
    #![allow(clippy::float_cmp)]
    let layout = ChunkLayout::<f64>::new(8, 4.0 / 35.0, [-1.0; 3]).expect("valid layout");
    let mut set = DirtySet::new();
    let id = ChunkId::new([2, -1, 3]);
    set.insert(id);

    let mut seen = None;
    set.mesh_within_budget(&layout, [0.0; 3], |_, origin| seen = Some(origin), || true);
    assert_eq!(seen, Some(layout.sample_origin(id)));
}

/// An empty queue is not an error and does no work.
#[test]
fn an_empty_queue_meshes_nothing() {
    let layout = ChunkLayout::<f64>::new(8, 0.125, [0.0; 3]).expect("valid layout");
    let mut set = DirtySet::new();
    let mut called = 0;
    let report = set.mesh_within_budget(&layout, [0.0; 3], |_, _| called += 1, || true);
    assert_eq!(called, 0);
    assert_eq!(report, BudgetReport::default());
    assert!(report.is_drained());
}
