//! `mark_edit_subgrid` exists to predict one thing: whether a chunk's Subgrid
//! mesh changed. So the tests that matter extract every chunk before and after
//! an edit and hold the prediction to the mesh.

use alloc::boxed::Box;
use alloc::vec;
use alloc::vec::Vec;
use core::cell::Cell;

use super::mark_edit_subgrid;
use crate::chunk::dirty::{DirtySet, EditReport, mark_edit};
use crate::chunk::{ChunkId, ChunkLayout};
use crate::fields::{Difference, Sphere, Union};
use crate::subgrid::extract::SubgridMarchingTetrahedra;
use crate::{Error, MeshBuffer, Sdf};

const SAMPLES: u32 = 4;

/// Four chunks of 8 cells along x, 0.25 apart: chunk `[x, 0, 0]` covers world x
/// in `[2x, 2x + 2)`. A power-of-two spacing, per M-32.
fn layout() -> ChunkLayout<f64> {
    ChunkLayout::new(8, 0.25, [0.0; 3]).expect("a valid layout")
}

const ROW: ([i64; 3], [i64; 3]) = ([0, 0, 0], [31, 7, 7]);

fn row_chunks() -> impl Iterator<Item = ChunkId> {
    (0..4).map(|x| ChunkId::new([x, 0, 0]))
}

/// A slab `2 * half` thick, perpendicular to x, centred at `at`. Negative inside.
struct Plate {
    at: f64,
    half: f64,
}

impl Sdf for Plate {
    type Scalar = f64;
    fn sample(&self, p: [f64; 3]) -> f64 {
        (p[0] - self.at).abs() - self.half
    }
}

/// Solid everywhere, a hair inside.
struct Solid;

impl Sdf for Solid {
    type Scalar = f64;
    fn sample(&self, _: [f64; 3]) -> f64 {
        -0.01
    }
}

/// [`Solid`] with a crack `2 * half` wide at `at`, which is outside.
///
/// At the corners of the cell the crack sits in, `half - |x - at|` is below
/// `-0.01`, so **no corner value differs from `Solid`'s**, and the crack is
/// still there between them.
struct Cracked {
    at: f64,
    half: f64,
}

impl Sdf for Cracked {
    type Scalar = f64;
    fn sample(&self, p: [f64; 3]) -> f64 {
        (-0.01f64).max(self.half - (p[0] - self.at).abs())
    }
}

type Field = Box<dyn Sdf<Scalar = f64>>;

fn sphere(center: [f64; 3], radius: f64) -> Sphere<f64> {
    Sphere { center, radius }
}

fn mesh(field: &Field, l: &ChunkLayout<f64>, id: ChunkId) -> MeshBuffer<f64> {
    let mut out = MeshBuffer::new();
    let shape = l.sample_shape().expect("a sample shape");
    SubgridMarchingTetrahedra::<f64>::new(SAMPLES)
        .expect("a valid resolution")
        .extract(field, &shape, l.sample_origin(id), l.cell_size(), &mut out)
        .expect("the chunk extracts");
    out
}

/// The chunks of the row whose Subgrid mesh differs between the two fields.
fn chunks_that_changed(l: &ChunkLayout<f64>, before: &Field, after: &Field) -> Vec<ChunkId> {
    row_chunks()
        .filter(|id| mesh(before, l, *id) != mesh(after, l, *id))
        .collect()
}

fn marked_by_subgrid(
    l: &ChunkLayout<f64>,
    before: &Field,
    after: &Field,
    samples: u32,
) -> (Vec<ChunkId>, EditReport) {
    let mut dirty = DirtySet::new();
    let report = mark_edit_subgrid(l, before, after, ROW.0, ROW.1, samples, &mut dirty)
        .expect("a region in bounds");
    (dirty.iter().collect(), report)
}

fn marked_by_corners(
    l: &ChunkLayout<f64>,
    before: &Field,
    after: &Field,
) -> (Vec<ChunkId>, EditReport) {
    let mut dirty = DirtySet::new();
    let report =
        mark_edit(l, before, after, ROW.0, ROW.1, &mut dirty).expect("a region in bounds");
    (dirty.iter().collect(), report)
}

/// Edits of every kind a caller makes, each of which changes some chunk's mesh.
fn edits() -> Vec<(&'static str, Field, Field)> {
    let plate = || Plate {
        at: 4.125,
        half: 0.075,
    };
    let with_plate = || -> Field {
        Box::new(Union {
            a: sphere([1.0, 1.0, 1.0], 0.5),
            b: plate(),
        })
    };
    let alone = || -> Field { Box::new(sphere([1.0, 1.0, 1.0], 0.5)) };
    let crack = || Cracked {
        at: 1.125,
        half: 0.075,
    };
    vec![
        (
            "a sphere grows in place",
            alone(),
            Box::new(sphere([1.0, 1.0, 1.0], 0.6)),
        ),
        (
            "a sphere moves from chunk 0 into chunk 2",
            alone(),
            Box::new(sphere([5.0, 1.0, 1.0], 0.5)),
        ),
        ("a plate thinner than a cell is added", alone(), with_plate()),
        ("that plate is removed", with_plate(), alone()),
        (
            "a crack no corner sees is cut into solid",
            Box::new(Solid),
            Box::new(crack()),
        ),
        (
            "that crack is closed",
            Box::new(crack()),
            Box::new(Solid),
        ),
        (
            "a plate is added beside a sphere's surface",
            Box::new(sphere([1.0, 1.0, 1.0], 1.1)),
            Box::new(Union {
                a: sphere([1.0, 1.0, 1.0], 1.1),
                b: Plate {
                    at: 2.125,
                    half: 0.075,
                },
            }),
        ),
        (
            "a brush carves across the seam between chunks 0 and 1",
            Box::new(sphere([1.5, 1.0, 1.0], 0.9)),
            Box::new(Difference {
                a: sphere([1.5, 1.0, 1.0], 0.9),
                b: sphere([2.0, 1.0, 1.0], 0.3),
            }),
        ),
    ]
}

/// **The soundness oracle.** Whatever an edit does, every chunk whose Subgrid
/// mesh differs afterwards must have been marked, or the caller keeps drawing
/// and colliding with a mesh the field no longer describes. Over-marking costs a
/// re-mesh; under-marking is a stale world.
///
/// Reports how many chunks were marked without changing, which is the price.
#[test]
fn every_chunk_whose_subgrid_mesh_changed_is_marked() {
    let l = layout();
    let (mut changed_total, mut marked_total, mut over) = (0, 0, 0);
    for (name, before, after) in edits() {
        let changed = chunks_that_changed(&l, &before, &after);
        assert!(
            !changed.is_empty(),
            "{name}: the edit changes no chunk's mesh, so it proves nothing"
        );
        let (marked, _) = marked_by_subgrid(&l, &before, &after, SAMPLES);
        for id in &changed {
            assert!(
                marked.contains(id),
                "{name}: chunk {id:?} changed and was not marked (marked: {marked:?})"
            );
        }
        let wasted = marked.iter().filter(|id| !changed.contains(id)).count();
        std::println!(
            "{name}: changed {}, marked {}, marked without changing {wasted}",
            changed.len(),
            marked.len()
        );
        changed_total += changed.len();
        marked_total += marked.len();
        over += wasted;
    }
    std::println!("total: changed {changed_total}, marked {marked_total}, over-marked {over}");
}

/// Issue #11's witness: every corner of the cell the plate sits in is outside
/// before and after, and the chunk's mesh changes.
#[test]
fn a_plate_thinner_than_a_cell_dirties_the_chunk_that_meshes_it() {
    let l = layout();
    let all = edits();
    for name in ["a plate thinner than a cell is added", "that plate is removed"] {
        let (_, before, after) = all.iter().find(|(n, _, _)| *n == name).expect("an edit");
        let (marked, report) = marked_by_subgrid(&l, before, after, SAMPLES);
        assert!(marked.contains(&ChunkId::new([2, 0, 0])), "{name}: {marked:?}");
        assert!(report.output_changed_cells > 0, "{name}");
    }
}

/// The contrast that makes the new rule necessary: a crack between the corners
/// changes no corner value, so the corner rule sees nothing and the mesh changes.
#[test]
fn a_crack_no_corner_sees_is_missed_by_mark_edit_and_found_by_mark_edit_subgrid() {
    let l = layout();
    let all = edits();
    let (_, before, after) = all
        .iter()
        .find(|(n, _, _)| n.starts_with("a crack no corner sees"))
        .expect("an edit");
    assert_eq!(
        chunks_that_changed(&l, before, after),
        vec![ChunkId::new([0, 0, 0])],
        "the crack must really change the mesh"
    );

    let (corner_marked, corner_report) = marked_by_corners(&l, before, after);
    assert_eq!(corner_report.value_changed_cells, 0, "no corner value moved");
    assert!(corner_marked.is_empty());

    let (subgrid_marked, _) = marked_by_subgrid(&l, before, after, SAMPLES);
    assert_eq!(subgrid_marked, vec![ChunkId::new([0, 0, 0])]);
}

/// Marking is not "everything". An edit that stays in one chunk leaves the others
/// to keep their meshes.
#[test]
fn an_edit_inside_one_chunk_leaves_the_others_alone() {
    let l = layout();
    let before: Field = Box::new(sphere([1.0, 1.0, 1.0], 0.5));
    let after: Field = Box::new(sphere([1.0, 1.0, 1.0], 0.6));
    let (marked, _) = marked_by_subgrid(&l, &before, &after, SAMPLES);
    assert_eq!(marked, vec![ChunkId::new([0, 0, 0])]);
}

/// At one sample per edge the extractor reads the corners and nothing else, and
/// the rule must collapse to `mark_edit`'s, counter for counter.
#[test]
fn one_sample_per_edge_is_the_corner_rule() {
    let l = layout();
    for (name, before, after) in edits() {
        let (corner_marked, corner_report) = marked_by_corners(&l, &before, &after);
        let (marked, report) = marked_by_subgrid(&l, &before, &after, 1);
        assert_eq!(report, corner_report, "{name}: the reports differ");
        assert_eq!(marked, corner_marked, "{name}: the marked chunks differ");
    }
}

/// Counts the evaluations made of it.
struct Counting<'a>(&'a Cell<u64>);

impl Sdf for Counting<'_> {
    type Scalar = f64;
    fn sample(&self, _: [f64; 3]) -> f64 {
        self.0.set(self.0.get() + 1);
        1.0
    }
}

/// The cost is a claim in the docs, so it is counted: every distinct point once,
/// however many cells share it.
#[test]
fn each_probe_point_is_evaluated_once() {
    let l = layout();
    let [cx, cy, cz] = [3u64, 2, 4];
    for samples in [1u32, 2, 4, 5] {
        let s = u64::from(samples) - 1;
        let corners = (cx + 1) * (cy + 1) * (cz + 1);
        let edges = cx * (cy + 1) * (cz + 1)
            + (cx + 1) * cy * (cz + 1)
            + (cx + 1) * (cy + 1) * cz
            + cx * cy * (cz + 1)
            + cx * (cy + 1) * cz
            + (cx + 1) * cy * cz
            + cx * cy * cz;
        let expected = corners + edges * s;

        let (n_before, n_after) = (Cell::new(0), Cell::new(0));
        let mut dirty = DirtySet::new();
        mark_edit_subgrid(
            &l,
            &Counting(&n_before),
            &Counting(&n_after),
            [0, 0, 0],
            [cx as i64 - 1, cy as i64 - 1, cz as i64 - 1],
            samples,
            &mut dirty,
        )
        .expect("a region in bounds");
        assert_eq!(n_before.get(), expected, "samples {samples}: before");
        assert_eq!(n_after.get(), expected, "samples {samples}: after");
    }
}

#[test]
fn an_identical_field_dirties_nothing() {
    let l = layout();
    let field: Field = Box::new(sphere([1.0, 1.0, 1.0], 0.5));
    let (marked, report) = marked_by_subgrid(&l, &field, &field, SAMPLES);
    assert!(marked.is_empty());
    assert_eq!(report.value_changed_cells, 0);
    assert_eq!(report.dirty_chunks, 0);
    assert!(report.region_cells > 0);
}

/// `dirty_chunks` is this edit's count, however long the caller's queue already
/// is: a backlog is not something the edit did.
#[test]
fn dirty_chunks_counts_this_edit_and_not_the_queue() {
    let l = layout();
    let field: Field = Box::new(sphere([1.0, 1.0, 1.0], 0.5));
    let mut dirty = DirtySet::new();
    dirty.insert(ChunkId::new([0, 0, 0]));
    let report = mark_edit_subgrid(&l, &field, &field, ROW.0, ROW.1, SAMPLES, &mut dirty)
        .expect("a region in bounds");
    assert_eq!(report.dirty_chunks, 0);
    assert_eq!(dirty.len(), 1, "the caller's queue is kept");
}

#[test]
fn bad_arguments_are_refused_rather_than_panicked_on() {
    let l = layout();
    let field: Field = Box::new(sphere([1.0, 1.0, 1.0], 0.5));
    let mut dirty = DirtySet::new();

    // The extractor cannot be built with no samples, whatever the region.
    for (min, max) in [(ROW.0, ROW.1), ([5, 5, 5], [4, 5, 5])] {
        let r = mark_edit_subgrid(&l, &field, &field, min, max, 0, &mut dirty);
        assert!(matches!(r, Err(Error::InvalidCellSize { .. })), "{r:?}");
    }

    // An inverted region is nothing to do.
    let r = mark_edit_subgrid(&l, &field, &field, [5, 5, 5], [4, 5, 5], SAMPLES, &mut dirty)
        .expect("an inverted region is empty, not an error");
    assert_eq!(r, EditReport::default());

    // An enormous one is refused before anything is allocated.
    let r = mark_edit_subgrid(
        &l,
        &field,
        &field,
        [i64::MIN, i64::MIN, i64::MIN],
        [i64::MAX, i64::MAX, i64::MAX],
        SAMPLES,
        &mut dirty,
    );
    assert!(matches!(r, Err(Error::IndexSpaceExhausted { .. })), "{r:?}");
    // A region the corner rule accepts, but whose probe grid does not fit.
    let r = mark_edit_subgrid(
        &l,
        &field,
        &field,
        [0, 0, 0],
        [1000, 1000, 1000],
        SAMPLES,
        &mut dirty,
    );
    assert!(matches!(r, Err(Error::IndexSpaceExhausted { .. })), "{r:?}");
    assert!(dirty.is_empty());
}
