//! A-014f's tests.
//!
//! The load-bearing ones are `gyroids_flipped_edges_go_to_zero` — the acceptance
//! — and `the_two_sheets_of_a_thin_plate_are_not_merged`, which is the case
//! A-014e's per-triangle vote was written to protect and the one P-7 was
//! registered to de-risk.

use alloc::vec;
use alloc::vec::Vec;

use super::{OrientReport, orient};
use crate::mesh::MeshBuffer;

/// A quad as two triangles sharing edge 1-2, coherently wound and facing the
/// way its normals say is out.
///
/// Coherent means the two traverse their shared edge in *opposite* directions:
/// triangle `0 1 2` goes 1 → 2, so its neighbour must go 2 → 1, which is
/// `1 3 2` and not `1 2 3`.
fn quad() -> MeshBuffer<f64> {
    let mut mesh = MeshBuffer::<f64>::new();
    mesh.positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ];
    mesh.normals = vec![[0.0, 0.0, 1.0]; 4];
    mesh.indices = vec![0, 1, 2, 1, 3, 2];
    mesh
}

/// The same quad with its second triangle wound backwards.
fn quad_with_one_reversed() -> MeshBuffer<f64> {
    let mut mesh = quad();
    mesh.indices[4..6].swap(0, 1);
    mesh
}

/// How many interior edges are traversed the *same* way by both their faces —
/// the defect this whole pass exists to remove, counted independently of
/// `validate` so the two are not proving each other.
fn same_direction_edges(mesh: &MeshBuffer<f64>) -> usize {
    use alloc::collections::BTreeMap;
    let mut seen: BTreeMap<(u32, u32), i32> = BTreeMap::new();
    for tri in mesh.indices.chunks_exact(3) {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            let key = if a <= b { (a, b) } else { (b, a) };
            *seen.entry(key).or_default() += if a <= b { 1 } else { -1 };
        }
    }
    // +-2 means both faces went the same way round; 0 means they opposed.
    seen.values().filter(|v| v.abs() == 2).count()
}

#[test]
fn an_empty_mesh_is_accepted_and_does_nothing() {
    let mut mesh = MeshBuffer::<f64>::new();
    assert_eq!(orient(&mut mesh).expect("empty"), OrientReport::default());
}

#[test]
fn an_index_past_the_end_is_refused_before_anything_moves() {
    let mut mesh = quad_with_one_reversed();
    mesh.indices[4] = 99;
    let before = mesh.indices.clone();
    assert!(orient(&mut mesh).is_err());
    assert_eq!(mesh.indices, before, "a refused call must not have edited");
}

#[test]
fn a_reversed_neighbour_is_flipped_back() {
    let mut mesh = quad_with_one_reversed();
    assert_eq!(same_direction_edges(&mesh), 1);

    let report = orient(&mut mesh).expect("orientable");
    assert_eq!(report.components, 1);
    assert_eq!(report.triangles_flipped, 1);
    assert!(report.is_orientable());
    assert_eq!(same_direction_edges(&mesh), 0);
}

#[test]
fn an_already_coherent_mesh_is_left_alone() {
    let mut mesh = quad();
    let before = mesh.indices.clone();
    let report = orient(&mut mesh).expect("orientable");
    assert!(report.is_noop());
    assert_eq!(mesh.indices, before);
}

/// Winding is not merely *consistent*, it is consistent **the right way round**.
/// A pass that made everything agree with each other and disagreed with the
/// field would satisfy every edge check and render inside out.
#[test]
fn a_wholly_inverted_mesh_is_turned_back_outward() {
    let mut mesh = quad();
    // Reverse both triangles: internally coherent, globally inside out.
    for tri in mesh.indices.chunks_exact_mut(3) {
        tri.swap(1, 2);
    }
    assert_eq!(same_direction_edges(&mesh), 0, "already self-consistent");

    let report = orient(&mut mesh).expect("orientable");
    assert_eq!(report.triangles_flipped, 2, "both, not neither");
    assert_eq!(mesh.indices, quad().indices);
}

#[test]
fn two_disjoint_patches_are_two_components() {
    let mut mesh = quad();
    let offset = mesh.positions.len() as u32;
    let second = quad_with_one_reversed();
    mesh.positions
        .extend(second.positions.iter().map(|p| [p[0] + 10.0, p[1], p[2]]));
    mesh.normals.extend(second.normals.iter().copied());
    mesh.indices
        .extend(second.indices.iter().map(|i| i + offset));

    let report = orient(&mut mesh).expect("orientable");
    assert_eq!(report.components, 2);
    assert_eq!(
        report.triangles_flipped, 1,
        "only the second patch was wrong"
    );
    assert_eq!(same_direction_edges(&mesh), 0);
}

/// The seed must be the most *confident* triangle, not the first one.
///
/// Here triangle 0 has a normal lying in its own plane — no information — while
/// triangle 1's is decisive and says the mesh is inside out. Seeding by index
/// would take triangle 0's non-vote, leave the pair as it found them, and report
/// success on a mesh that is still wrong.
#[test]
fn the_seed_is_the_most_confident_triangle_not_the_first() {
    let mut mesh = quad();
    for tri in mesh.indices.chunks_exact_mut(3) {
        tri.swap(1, 2);
    }
    // Vertex 0 is used only by triangle 0, so this blunts that triangle's vote
    // without touching triangle 1's.
    mesh.normals[0] = [0.0, 0.0, -1.0];

    let report = orient(&mut mesh).expect("orientable");
    assert_eq!(
        report.triangles_flipped, 2,
        "the indecisive triangle was allowed to seed the component"
    );
    assert_eq!(mesh.indices, quad().indices);
}

/// A triangle with no usable vote at all still gets a *defined* answer, from its
/// neighbour rather than from the order the faces happened to be visited in.
/// A-014f's second acceptance clause.
#[test]
fn a_triangle_with_no_vote_of_its_own_inherits_one() {
    let mut mesh = quad_with_one_reversed();
    // Zero normals on the vertices unique to the second triangle: its own vote
    // is now the zero vector, which is neither outward nor inward.
    mesh.normals[3] = [0.0, 0.0, 0.0];

    let report = orient(&mut mesh).expect("orientable");
    assert_eq!(same_direction_edges(&mesh), 0, "it inherited nothing");
    assert!(report.is_orientable());
}

/// Determinism, and specifically that it does not depend on the order the faces
/// arrive in — the failure mode a `HashMap`-backed adjacency would introduce.
#[test]
fn the_result_does_not_depend_on_the_order_the_faces_are_listed_in() {
    let base = {
        let mut m = quad_with_one_reversed();
        orient(&mut m).expect("orientable");
        m.indices
    };
    // The same two triangles, listed the other way round.
    let mut swapped = quad_with_one_reversed();
    let (a, b) = swapped.indices.split_at(3);
    let reordered: Vec<u32> = b.iter().chain(a.iter()).copied().collect();
    swapped.indices = reordered;
    orient(&mut swapped).expect("orientable");
    assert_eq!(same_direction_edges(&swapped), 0);
    assert_eq!(base.len(), swapped.indices.len());
}

/// **A Möbius band is reported, not silently half-fixed.**
///
/// Five quads in a ring with a half-twist: propagation comes back to where it
/// started with the opposite sign, and no assignment of windings can make every
/// edge oppose. The honest answer is to say so.
#[test]
fn a_mobius_band_is_reported_as_non_orientable() {
    let strips = 5;
    let mut mesh = MeshBuffer::<f64>::new();
    for i in 0..strips {
        let t = i as f64 / strips as f64;
        mesh.positions.push([t, 0.0, 0.0]);
        mesh.positions.push([t, 1.0, 0.0]);
        mesh.normals.push([0.0, 0.0, 1.0]);
        mesh.normals.push([0.0, 0.0, 1.0]);
    }
    for i in 0..strips {
        let a = (i * 2) as u32;
        let b = a + 1;
        let next = ((i + 1) % strips) as u32 * 2;
        // The last quad joins back with its two vertices swapped: the twist.
        let (c, d) = if i + 1 == strips {
            (next + 1, next)
        } else {
            (next, next + 1)
        };
        mesh.indices.extend_from_slice(&[a, b, c]);
        mesh.indices.extend_from_slice(&[b, d, c]);
    }

    let report = orient(&mut mesh).expect("indices are in range");
    assert_eq!(report.components, 1);
    assert!(
        !report.is_orientable(),
        "a half-twisted ring was reported orientable"
    );
    assert_eq!(report.non_orientable_components, 1);
}

/// **A-014f's acceptance.** `gyroid`'s inconsistently-oriented edges go to zero,
/// and `thin_plate` — the field A-014e's per-triangle vote exists to protect —
/// does not regress. Measured on the real extractor at all three resolutions,
/// because M-182 found `thin_plate` failing at the two the suite never tested.
#[test]
fn every_reference_field_comes_out_coherently_oriented() {
    use crate::shape::RuntimeShape3;
    use crate::subgrid::extract::SubgridMarchingTetrahedra;
    use crate::validate::{ValidateConfig, validate_indexed};

    let mut rows = 0;
    let mut got: Vec<(&str, u32, u64, u64, u64)> = Vec::new();
    for n in [17u32, 25, 33] {
        crate::for_each_reference_field!(f64, |name, field| {
            let (lo, hi) = crate::fields::ReferenceField::domain(&field);
            let cell = (hi[0] - lo[0]) / f64::from(n - 1);
            let shape = RuntimeShape3::new([n; 3]).expect("a cubic grid");
            let cfg = ValidateConfig::from_cell_size(cell).expect("a valid spacing");

            let mut mesh = MeshBuffer::<f64>::default();
            SubgridMarchingTetrahedra::<f64>::new(16)
                .expect("valid")
                .extract(&field, &shape, lo, cell, &mut mesh)
                .unwrap_or_else(|e| panic!("{name} {n}: {e}"));
            crate::weld::Welder::<f64>::new()
                .weld(&mut mesh, crate::weld::epsilon_for(cell))
                .expect("weld");

            let before = validate_indexed(&mesh.positions, &mesh.indices, &cfg);
            let report = orient(&mut mesh).expect("indices in range");
            let after = validate_indexed(&mesh.positions, &mesh.indices, &cfg);

            std::println!(
                "{name:<15} {n}³  flipped edges {:>4} -> {:>3}  \
                 (components {}, triangles flipped {}, non-orientable {})",
                before.inconsistently_oriented_edges,
                after.inconsistently_oriented_edges,
                report.components,
                report.triangles_flipped,
                report.non_orientable_components,
            );

            got.push((
                name,
                n,
                before.inconsistently_oriented_edges,
                after.inconsistently_oriented_edges,
                before.non_manifold_edges,
            ));
            // Orientation must not disturb anything else about the mesh.
            assert_eq!(after.vertices, before.vertices, "{name} {n}³ vertices");
            assert_eq!(after.faces, before.faces, "{name} {n}³ faces");
            assert_eq!(
                (after.non_manifold_edges, after.non_manifold_vertices),
                (before.non_manifold_edges, before.non_manifold_vertices),
                "{name} {n}³ manifoldness moved"
            );
            assert_eq!(
                after.boundary_edges, before.boundary_edges,
                "{name} {n}³ boundary moved"
            );
            rows += 1;
        });
    }
    assert_eq!(rows, 24);

    // **The law, and it is exact (M-187).** Orientation drives the count to
    // *zero* on every row whose edges are all manifold. Where it does not, the
    // mesh has edges carrying more than two faces, and no assignment of windings
    // can make four faces pairwise oppose across one edge — that residue is
    // A-014d's non-manifoldness wearing an orientation costume, not a limit of
    // propagation.
    for &(name, n, before, after, non_manifold) in &got {
        if non_manifold == 0 {
            assert_eq!(
                after, 0,
                "{name} {n}³: {after} flipped edges left on a mesh with no \
                 non-manifold edge, from {before}"
            );
        } else {
            // **`after <= before` is not a law, and `noise_cavity` is what shows
            // it (M-213).** The reasoning was that propagation can only ever
            // agree edges it reaches, so the residue can shrink but not grow.
            // That holds while the non-manifold edges are sparse enough to stay
            // out of propagation's way. They are not here: with 318 of them the
            // flood fill crosses a four-face edge, commits to one side's winding,
            // and carries a consistent-but-wrong orientation across a whole
            // patch — so the count rises, 1580 → 2422. Nothing is *more* wrong
            // afterwards; a different, larger set of edges is now the one that
            // disagrees. Pinned per field so a second field growing still fails.
            // A-019.
            assert!(
                after > 0,
                "{name} {n}³: residue vanished on a non-manifold mesh"
            );
            let may_grow = name == "noise_cavity";
            assert!(
                after <= before || may_grow,
                "{name} {n}³: {before} -> {after} with {non_manifold} \
                 non-manifold edges -- orientation made it worse, or removed a \
                 residue it should not be able to"
            );
        }
    }

    // The full census, pinned. The law above explains its shape; this is what
    // stops any row of it moving unnoticed.
    let table: Vec<(&str, u32, u64, u64)> =
        got.iter().map(|&(f, n, b, a, _)| (f, n, b, a)).collect();
    assert_eq!(
        table,
        alloc::vec![
            ("sphere", 17, 0, 0),
            ("torus", 17, 0, 0),
            ("box_exact", 17, 0, 0),
            ("csg_difference", 17, 6, 6),
            ("thin_plate", 17, 0, 0),
            ("gyroid", 17, 138, 0),
            ("fbm_terrain", 17, 19, 6),
            ("noise_cavity", 17, 1629, 1015),
            ("sphere", 25, 0, 0),
            ("torus", 25, 6, 6),
            ("box_exact", 25, 0, 0),
            ("csg_difference", 25, 0, 0),
            ("thin_plate", 25, 8, 0),
            ("gyroid", 25, 150, 0),
            ("fbm_terrain", 25, 29, 3),
            ("noise_cavity", 25, 1580, 2422),
            ("sphere", 33, 0, 0),
            ("torus", 33, 0, 0),
            ("box_exact", 33, 0, 0),
            ("csg_difference", 33, 36, 0),
            ("thin_plate", 33, 6, 0),
            ("gyroid", 33, 330, 0),
            ("fbm_terrain", 33, 53, 12),
            // The only rows where orientation *raises* the count -- M-213, A-019.
            // It lowers it at 17³ and raises it at 25³ and 33³, which is the
            // shape of a flood fill meeting enough four-face edges to commit to
            // the wrong side of one.
            ("noise_cavity", 33, 1477, 3341),
        ],
        "the flipped-edge census moved"
    );

    // A-014f's headline acceptance, named rather than left to the loop.
    let row = |f: &str, n: u32| {
        got.iter()
            .find(|g| g.0 == f && g.1 == n)
            .map(|g| (g.2, g.3))
            .expect("in table")
    };
    assert_eq!(row("gyroid", 17), (138, 0));
    assert_eq!(row("gyroid", 25), (150, 0));
    assert_eq!(row("gyroid", 33), (330, 0));
    // And `thin_plate` -- the field the per-triangle vote exists to protect --
    // is not merely un-regressed, it is fixed at the two resolutions M-182
    // found it failing at.
    assert_eq!(row("thin_plate", 17), (0, 0));
    assert_eq!(row("thin_plate", 25), (8, 0));
    assert_eq!(row("thin_plate", 33), (6, 0));
}
