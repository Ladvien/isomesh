use alloc::vec;
use alloc::vec::Vec;

use super::*;
use crate::extractor::Extractor;
use crate::fields::{ReferenceField, Sphere};
use crate::marching_cubes::MarchingCubes;
use crate::weld::{Welder, epsilon_for};
use crate::{MeshBuffer, RuntimeShape3};

/// One world unit per cell, so `weld_epsilon` is `1e-4` absolute and every
/// fixture below can be written in whole units without thinking about scale.
fn cfg() -> ValidateConfig {
    ValidateConfig::from_cell_size(1.0).expect("a positive cell size")
}

fn census(positions: &[[f64; 3]], indices: &[u32]) -> (PinchReport, PinchGroups) {
    pinch_features(positions, indices, &cfg())
}

/// **A fold, which is the case the repair exists to remove.**
///
/// Two vertices at one point, both corners of the same triangle, with a second
/// triangle sharing that triangle's edge. Merging them flattens the fold: the
/// degenerate face goes and nothing else can move.
#[test]
fn a_fold_is_not_a_pinch() {
    let positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ];
    let indices = [0, 1, 2, 0, 2, 3];
    let (report, groups) = census(&positions, &indices);

    assert_eq!(report.collapse_groups, 1);
    assert_eq!(report.collapsing_vertices, 2);
    assert_eq!(report.pinch_groups, 0);
    assert_eq!(report.pieces_joined, 0);
    // Only the one degenerate face is dropped, and only its own edge unions.
    assert_eq!(report.folding_faces, 1);
    assert_eq!(report.sharing_edges, 1);
    assert!(report.is_pinch_free());
    assert!(report.moves_no_geometry());

    assert_eq!(groups.group_count(), 1);
    assert_eq!(groups.members(0), [0, 1]);
    assert_eq!(groups.clusters_in(0), 1);
    assert!(!groups.clusters_are_split(0));
    // The label is the least member of the cluster, whichever order the unions
    // arrived in.
    assert_eq!(groups.clusters_of(0), [0, 0]);
}

/// **A pinch, which is the case that welds two pieces.**
///
/// Two triangles sharing no index and exactly one position. Identifying that
/// point joins two pieces the mesh kept apart, and no face is degenerate.
#[test]
fn two_sheets_touching_at_a_point_are_a_pinch() {
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    ];
    let indices = [0, 1, 2, 3, 4, 5];
    let (report, groups) = census(&positions, &indices);

    assert_eq!(report.collapse_groups, 1);
    assert_eq!(report.pinch_groups, 1);
    assert_eq!(report.pieces_joined, 1);
    assert_eq!(report.folding_faces, 0);
    // The instrument found no sharing edge to find, and says so rather than
    // leaving the pinch to be believed on faith.
    assert_eq!(report.sharing_edges, 0);
    assert!(!report.is_pinch_free());

    assert_eq!(groups.members(0), [0, 3]);
    assert_eq!(groups.clusters_of(0), [0, 3]);
    assert_eq!(groups.clusters_in(0), 2);
    assert!(groups.clusters_are_split(0));
}

/// **Three sheets at one point join two pieces, not three.**
///
/// This is the difference between `M-352`'s 516 pinch groups and its 520 welded
/// pieces, in the smallest fixture that shows it: a group can span more clusters
/// than two, and the piece count is the cluster count less one.
#[test]
fn a_three_way_touch_joins_two_pieces() {
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
    ];
    let indices = [0, 1, 2, 3, 4, 5, 6, 7, 8];
    let (report, groups) = census(&positions, &indices);

    assert_eq!(report.collapse_groups, 1);
    assert_eq!(report.pinch_groups, 1);
    assert_eq!(report.pieces_joined, 2);
    assert_eq!(groups.clusters_in(0), 3);
    assert_eq!(groups.members(0), [0, 3, 6]);
    assert_eq!(groups.clusters_of(0), [0, 3, 6]);
}

/// **`P-125`'s C2, asked of the unit test as well as of the harness.**
///
/// `✗26` was a face-iteration-order leak found after the fact. The union-find
/// here consumes the faces in whatever order they arrive, so the report *and*
/// the feature lists have to be identical under any permutation of them —
/// including the cluster labels, which is why the union keeps the lower root
/// rather than the larger set.
#[test]
fn the_census_does_not_depend_on_face_order() {
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        [5.0, 0.0, 0.0],
        [5.0, 0.0, 0.0],
        [5.0, 1.0, 0.0],
        [6.0, 1.0, 0.0],
    ];
    let faces: [[u32; 3]; 4] = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [6, 8, 9]];
    let flat: Vec<u32> = faces.iter().flatten().copied().collect();
    let expected = census(&positions, &flat);

    // One pinch and one fold, so both branches of the predicate are exercised
    // by every permutation rather than only one of them.
    assert_eq!(expected.0.collapse_groups, 2);
    assert_eq!(expected.0.pinch_groups, 1);
    assert_eq!(expected.0.folding_faces, 1);

    // All 24 orderings of four faces, enumerated rather than sampled.
    let mut seen = 0;
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let order = [a, b, c, d];
                    let mut sorted = order;
                    sorted.sort_unstable();
                    if sorted != [0, 1, 2, 3] {
                        continue;
                    }
                    seen += 1;
                    let permuted: Vec<u32> =
                        order.iter().flat_map(|&f| faces[f]).collect::<Vec<u32>>();
                    assert_eq!(census(&positions, &permuted), expected, "order {order:?}");
                }
            }
        }
    }
    assert_eq!(seen, 24, "every permutation of four faces");
}

/// **The instrument has to be able to read zero, and to read it for the right
/// reason.**
///
/// No two positions within `weld_epsilon`, so there is nothing to collapse. The
/// zero is a property of the mesh rather than of the lattice: the fixture's
/// vertices are a whole unit apart and the epsilon is `1e-4`.
#[test]
fn a_mesh_with_no_coincidence_reports_nothing() {
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let indices = [0, 1, 2];
    let (report, groups) = census(&positions, &indices);

    assert_eq!(report.collapse_groups, 0);
    assert_eq!(report.collapsing_vertices, 0);
    assert_eq!(report.pinch_groups, 0);
    assert_eq!(report.folding_faces, 0);
    assert_eq!(report.sharing_edges, 0);
    assert_eq!(report.triangles, 1);
    assert!(report.is_pinch_free());
    assert_eq!(groups.group_count(), 0);
    assert_eq!(groups, PinchGroups::default());
}

/// **A bucket is not the same claim as a coincidence, and the report says which
/// it found.**
///
/// Two vertices a micron apart land in one `1e-4` cell without being the same
/// point. `M-352`'s repair moved nothing at all, so a caller reading a pinch
/// count needs to know whether the collapse it is being told about would also
/// move geometry.
#[test]
fn a_group_whose_members_merely_share_a_cell_is_reported_as_moving_geometry() {
    let positions = [
        [0.0, 0.0, 0.0],
        [1e-6, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [2.0, 0.0, 0.0],
    ];
    let indices = [0, 2, 3, 1, 3, 4];
    let (report, groups) = census(&positions, &indices);

    assert_eq!(report.collapse_groups, 1);
    assert_eq!(groups.members(0), [0, 1]);
    assert_eq!(report.pinch_groups, 1);
    assert_eq!(report.groups_moving_geometry, 1);
    assert!(!report.moves_no_geometry());
}

/// **A malformed mesh is reported, not dereferenced.**
///
/// The same intake rule as `validate_features`: an index past the buffer cannot
/// be looked up and a repeated index has no edge set, so the face is skipped and
/// counted, and a partial triangle at the end is a column rather than a panic.
#[test]
fn a_malformed_face_is_skipped_rather_than_dereferenced() {
    let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    // One unusable face per reason, then two trailing indices.
    let indices = [0, 1, 99, 0, 0, 2, 0, 1];
    let report = pinch_census(&positions, &indices, &cfg());

    assert_eq!(report.triangles, 2);
    assert_eq!(report.trailing_indices, 2);
    assert_eq!(report.faces_skipped, 2);
    // The coincidence is still found — the groups come from the positions — but
    // no face was usable, so nothing could share a triangle and the collapse
    // reads as a pinch.
    assert_eq!(report.collapse_groups, 1);
    assert_eq!(report.sharing_edges, 0);
    assert_eq!(report.pinch_groups, 1);
}

/// **The census against the weld it is a census of, on a real extraction.**
///
/// `M-48` measured this exact case: on `sphere` at 25³ a grid sample lands on
/// the isosurface, `t` is 0 or 1 there, and every cut edge meeting the sample
/// places its own vertex at the same point — **48 vertices and 96 collapsed
/// triangles**. The welder's rule (lowest-indexed representative within `ε`,
/// probing 27 cells) and this one (ε-connected components over the same probe)
/// are not the same rule — the welder's classes refine these — so their
/// agreement is measured rather than assumed, exactly as `validate`'s duplicate
/// count is measured against the weld it bounds.
///
/// **This is also the measurement that rejected the cell-key reading of a
/// group.** Grouping by lattice cell instead of by ε-closure gives 17 groups
/// over 47 vertices here — 30 of the welder's 48 merges — because a coincidence
/// class straddling a cell face is split between two buckets. The closure sees
/// all 48.
#[test]
fn the_census_predicts_the_weld_on_a_real_extraction() {
    const N: u32 = 25;
    let field = Sphere::<f64>::canonical();
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(N - 1);
    let shape = RuntimeShape3::new([N; 3]).expect("valid shape");
    let mut mesh = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract_into(&field, &shape, lo, h, &mut mesh)
        .expect("extraction");

    let cfg = ValidateConfig::from_cell_size(h).expect("a positive cell size");
    let report = pinch_census(&mesh.positions, &mesh.indices, &cfg);

    // The population is non-empty before anything is concluded from it.
    assert!(
        report.collapse_groups > 0,
        "no coincidence on the one fixture M-48 measured one on:\n{report}"
    );
    assert_eq!(report.collapse_groups, 24);
    assert_eq!(report.collapsing_vertices - report.collapse_groups, 48);
    assert_eq!(report.folding_faces, 96);
    // Every coincidence here is a fold: the sphere is one sheet and it does not
    // touch itself. The 144 sharing edges are what licenses that zero.
    assert_eq!(report.sharing_edges, 144);
    assert!(report.is_pinch_free());

    let mut welded = mesh.clone();
    let weld = Welder::<f64>::new()
        .weld(&mut welded, epsilon_for(h))
        .expect("a positive epsilon");
    assert_eq!(
        u64::try_from(weld.vertices_removed()).expect("a count"),
        report.collapsing_vertices - report.collapse_groups,
        "the census and the weld disagree about how many vertices coincide"
    );
    assert_eq!(
        u64::try_from(weld.triangles_collapsed).expect("a count"),
        report.folding_faces,
        "the census and the weld disagree about which faces fold"
    );
}

/// **Compressed-row storage is only reserved once, and `len == capacity` is what
/// says so.**
///
/// `P-125`'s C2 is falsified by the report allocating per group. The workspace
/// forbids `unsafe_code`, so a counting allocator is unrepresentable here; a
/// `Vec` whose length equals its capacity is one that never grew, which is the
/// observable half of the same claim. Asserted on a fixture with several groups
/// so a single-group accident cannot pass it.
#[test]
fn the_returned_buffers_are_reserved_exactly_once() {
    let mut positions: Vec<[f64; 3]> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    for k in 0..7 {
        let x = f64::from(k) * 10.0;
        let base = u32::try_from(positions.len()).expect("a small fixture");
        positions.extend(vec![
            [x, 0.0, 0.0],
            [x + 1.0, 0.0, 0.0],
            [x, 1.0, 0.0],
            [x, 0.0, 0.0],
            [x - 1.0, 0.0, 0.0],
            [x, -1.0, 0.0],
        ]);
        indices.extend([base, base + 1, base + 2, base + 3, base + 4, base + 5]);
    }
    let (report, groups) = census(&positions, &indices);

    assert_eq!(report.collapse_groups, 7);
    assert_eq!(report.pinch_groups, 7);
    assert_eq!(groups.vertices.len(), groups.vertices.capacity());
    assert_eq!(groups.clusters.len(), groups.clusters.capacity());
    assert_eq!(groups.starts.len(), groups.starts.capacity());
    assert_eq!(groups.starts.len(), 8);
}
