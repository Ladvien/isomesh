//! Tests for A-013.
//!
//! The load-bearing one is `a_chunk_seam_welds_at_a_non_power_of_two_spacing`.
//! Everything else here is a unit test of the rule; that one is the reason the
//! module exists, and it is deliberately run at a spacing M-32 identified as
//! *not* bit-exact, because a seam test at `h = 0.125` passes without welding
//! anything and proves nothing.

use alloc::vec;
use alloc::vec::Vec;

use super::{WeldReport, Welder, epsilon_for};
use crate::chunk::{ChunkId, ChunkLayout};
use crate::fields::{ReferenceField, Sphere};
use crate::marching_cubes::MarchingCubes;
use crate::validate::{ValidateConfig, check_determinism, validate_indexed};
use crate::{Error, MeshBuffer, Sdf};

/// A buffer with the given positions, one normal each, and no triangles.
fn points(positions: &[[f64; 3]]) -> MeshBuffer<f64> {
    MeshBuffer::<f64> {
        positions: positions.to_vec(),
        normals: vec![[0.0, 0.0, 1.0]; positions.len()],
        indices: Vec::new(),
    }
}

// ─── the rule ───────────────────────────────────────────────────────────────

#[test]
fn coincident_vertices_collapse_to_the_lowest_index() {
    let mut mesh = points(&[[0.0; 3], [1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]]);
    mesh.indices = vec![0, 1, 2];
    let report = Welder::<f64>::new()
        .weld(&mut mesh, 1e-9)
        .expect("valid epsilon");

    assert_eq!(report.vertices_before, 4);
    assert_eq!(report.vertices_after, 2);
    assert_eq!(report.vertices_removed(), 2);
    // Vertices 0, 2 and 3 all became output vertex 0.
    assert_eq!(
        Welder::<f64>::new().remap().len(),
        0,
        "fresh welder is empty"
    );
    // And the triangle that used two of them collapsed.
    assert_eq!(report.triangles_collapsed, 1);
    assert_eq!(mesh.triangle_count(), 0);
}

/// The tie-break is the *lowest* index, not the first the 27-cell probe reaches
/// — otherwise the answer would depend on the order those cells are visited in,
/// which is an implementation detail and exactly the kind of thing that changes
/// under a refactor without anyone noticing.
#[test]
fn the_representative_is_the_lowest_index_not_the_nearest() {
    // Three points in a line, all within epsilon of each other.
    let mut mesh = points(&[[0.0; 3], [0.4e-9, 0.0, 0.0], [0.8e-9, 0.0, 0.0]]);
    let mut welder = Welder::<f64>::new();
    welder.weld(&mut mesh, 1e-9).expect("valid epsilon");

    // 2 is nearer to 1 than to 0, and still joins 0.
    assert_eq!(welder.remap(), &[0, 0, 0]);
    assert_eq!(mesh.positions, vec![[0.0; 3]]);
}

/// Epsilon-closeness is not transitive, and the rule does not pretend otherwise:
/// comparing against **kept** vertices only is what stops a chain of near-misses
/// dragging distant vertices together.
#[test]
fn a_chain_of_near_misses_does_not_merge_its_ends() {
    // 0 ~ 1 and 1 ~ 2, but 0 and 2 are two epsilons apart.
    let mut mesh = points(&[[0.0; 3], [0.9, 0.0, 0.0], [1.8, 0.0, 0.0]]);
    let mut welder = Welder::<f64>::new();
    let report = welder.weld(&mut mesh, 1.0).expect("valid epsilon");

    // 1 joins 0. 2 is then compared against kept vertices — which is 0 alone,
    // since 1 was welded away — and 0 is 1.8 away, so 2 is kept.
    assert_eq!(welder.remap(), &[0, 0, 1]);
    assert_eq!(report.vertices_after, 2);
    assert_eq!(mesh.positions, vec![[0.0; 3], [1.8, 0.0, 0.0]]);
}

#[test]
fn a_vertex_exactly_epsilon_away_welds_and_one_beyond_does_not() {
    let mut inside = points(&[[0.0; 3], [1.0, 0.0, 0.0]]);
    let report = Welder::<f64>::new()
        .weld(&mut inside, 1.0)
        .expect("valid epsilon");
    assert_eq!(
        report.vertices_after, 1,
        "exactly epsilon is within epsilon"
    );

    let mut outside = points(&[[0.0; 3], [1.0 + 1e-12, 0.0, 0.0]]);
    let report = Welder::<f64>::new()
        .weld(&mut outside, 1.0)
        .expect("valid epsilon");
    assert!(report.is_noop(), "one ulp beyond epsilon must not weld");
}

#[test]
fn a_mesh_with_nothing_to_weld_is_left_alone() {
    let mut mesh = points(&[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
    mesh.indices = vec![0, 1, 2];
    let before = mesh.clone();
    let report = Welder::<f64>::new()
        .weld(&mut mesh, 1e-9)
        .expect("valid epsilon");

    assert!(report.is_noop());
    assert_eq!(mesh.positions, before.positions);
    assert_eq!(mesh.normals, before.normals);
    assert_eq!(mesh.indices, before.indices);
}

#[test]
fn an_empty_mesh_is_not_an_error() {
    let mut mesh = MeshBuffer::<f64>::new();
    let report = Welder::<f64>::new()
        .weld(&mut mesh, 1e-9)
        .expect("valid epsilon");
    assert_eq!(report, WeldReport::default());
}

/// Zero is rejected rather than read as "weld exact matches only", which is the
/// case M-32 says fails at a seam.
#[test]
fn a_meaningless_epsilon_is_rejected() {
    for bad in [0.0f64, -1.0, f64::NAN, f64::INFINITY] {
        let mut mesh = points(&[[0.0; 3], [0.0; 3]]);
        let error = Welder::<f64>::new()
            .weld(&mut mesh, bad)
            .expect_err("epsilon must be finite and positive");
        assert!(matches!(error, Error::InvalidWeldEpsilon { .. }), "{bad}");
        assert_eq!(mesh.vertex_count(), 2, "a rejected weld changes nothing");
    }
}

/// A NaN position is bucketed at the origin and then fails the distance test, so
/// it never welds — including to another NaN.
#[test]
fn a_non_finite_position_never_welds() {
    let mut mesh = points(&[[f64::NAN, 0.0, 0.0], [f64::NAN, 0.0, 0.0], [0.0; 3]]);
    let report = Welder::<f64>::new()
        .weld(&mut mesh, 1.0)
        .expect("valid epsilon");
    assert_eq!(report.vertices_after, 3);
}

/// The remap has to be usable for the caller's own per-vertex data, or a mesh
/// with colours cannot be welded at all.
#[test]
fn the_remap_covers_every_input_vertex() {
    let mut mesh = points(&[[0.0; 3], [5.0, 0.0, 0.0], [0.0; 3], [5.0, 0.0, 0.0]]);
    let mut welder = Welder::<f64>::new();
    welder.weld(&mut mesh, 1e-9).expect("valid epsilon");
    assert_eq!(welder.remap(), &[0, 1, 0, 1]);
    assert!(
        welder
            .remap()
            .iter()
            .all(|&o| (o as usize) < mesh.vertex_count())
    );
}

// ─── the reason it exists ───────────────────────────────────────────────────

fn mesh_chunk(
    layout: &ChunkLayout<f64>,
    field: &impl Sdf<Scalar = f64>,
    chunk: ChunkId,
) -> MeshBuffer<f64> {
    let shape = layout.sample_shape().expect("valid shape");
    let mut out = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(
            field,
            &shape,
            layout.sample_origin(chunk),
            layout.cell_size(),
            &mut out,
        )
        .expect("extraction");
    out
}

/// **Two chunks, meshed independently, welded into one watertight surface.**
///
/// The spacing is `4/35` on purpose. M-32 measured that a chunk seam is
/// bit-exact only when the cell size is a power of two, and that 22% of random
/// `(origin, h, cells, chunk)` combinations disagree by an ulp — `4/35` is from
/// that search. A weld keyed on equality would pass at `h = 0.125` and fail
/// here, so this is the fixture that can tell the two apart.
///
/// The assertion is on the *seam*, not on the totals: before welding the two
/// chunk meshes share a plane of duplicated vertices and every triangle on it
/// has a boundary edge; after, the boundary edges on that plane are gone.
#[test]
fn a_chunk_seam_welds_at_a_non_power_of_two_spacing() {
    let h = 4.0 / 35.0;
    let layout = ChunkLayout::<f64>::new(8, h, [-2.0; 3]).expect("valid layout");
    let field = Sphere::<f64>::canonical();

    let a = ChunkId::new([1, 1, 1]);
    let b = a.neighbour(0, 1);

    let mut joined = mesh_chunk(&layout, &field, a);
    let second = mesh_chunk(&layout, &field, b);
    assert!(
        !joined.is_empty() && !second.is_empty(),
        "both chunks must carry surface or this measures nothing"
    );
    joined
        .append(&second)
        .expect("the meshes fit the u32 index space");

    let cfg = ValidateConfig::from_cell_size(h).expect("valid cell size");
    let before = validate_indexed(&joined.positions, &joined.indices, &cfg);

    let mut welder = Welder::<f64>::new();
    let report = welder
        .weld(&mut joined, epsilon_for(h))
        .expect("valid epsilon");
    let after = validate_indexed(&joined.positions, &joined.indices, &cfg);

    std::println!(
        "measured: A-013 seam at h = 4/35 -- {} vertices welded to {}, boundary edges {} -> {}, \
         duplicate_vertices {} -> {}, chi {} -> {}",
        report.vertices_before,
        report.vertices_after,
        before.boundary_edges,
        after.boundary_edges,
        before.duplicate_vertices,
        after.duplicate_vertices,
        before.euler_characteristic,
        after.euler_characteristic,
    );

    assert!(
        report.vertices_removed() > 0,
        "the seam was not welded at all"
    );
    assert!(
        after.boundary_edges < before.boundary_edges,
        "welding did not close any boundary: {} -> {}",
        before.boundary_edges,
        after.boundary_edges
    );
    // The welder's whole job: nothing coincides afterwards.
    assert_eq!(after.duplicate_vertices, 0, "{after}");
    assert_eq!(after.inconsistently_oriented_edges, 0, "{after}");
}

/// The same seam at a power-of-two spacing, where M-32 says the two chunks agree
/// bit for bit. It must weld there too — the point being that one epsilon covers
/// both cases, so there is no exact-match path and no spacing-dependent
/// behaviour.
#[test]
fn a_chunk_seam_welds_at_a_power_of_two_spacing_too() {
    let h = 0.125;
    let layout = ChunkLayout::<f64>::new(8, h, [-2.0; 3]).expect("valid layout");
    let field = Sphere::<f64>::canonical();

    let a = ChunkId::new([1, 1, 1]);
    let mut joined = mesh_chunk(&layout, &field, a);
    let second = mesh_chunk(&layout, &field, a.neighbour(0, 1));
    assert!(!joined.is_empty() && !second.is_empty());
    joined
        .append(&second)
        .expect("the meshes fit the u32 index space");

    let cfg = ValidateConfig::from_cell_size(h).expect("valid cell size");
    let mut welder = Welder::<f64>::new();
    let report = welder
        .weld(&mut joined, epsilon_for(h))
        .expect("valid epsilon");
    let after = validate_indexed(&joined.positions, &joined.indices, &cfg);

    assert!(report.vertices_removed() > 0);
    assert_eq!(after.duplicate_vertices, 0, "{after}");
}

// ─── T-004 ──────────────────────────────────────────────────────────────────

/// Weld ordering is the classic determinism leak, so this is not optional.
///
/// The extraction and the weld run together inside the harness, three times
/// including once into a reused buffer, and the comparison is bit-level through
/// `total_cmp` rather than `==`.
#[test]
fn welding_is_deterministic() {
    let h = 4.0 / 35.0;
    let layout = ChunkLayout::<f64>::new(8, h, [-2.0; 3]).expect("valid layout");
    let field = Sphere::<f64>::canonical();
    let a = ChunkId::new([1, 1, 1]);
    let b = a.neighbour(0, 1);

    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        let first = mesh_chunk(&layout, &field, a);
        let second = mesh_chunk(&layout, &field, b);
        out.append(&first)
            .expect("the meshes fit the u32 index space");
        out.append(&second)
            .expect("the meshes fit the u32 index space");
        Welder::<f64>::new()
            .weld(out, epsilon_for(h))
            .expect("valid epsilon");
    });
    assert!(report.is_deterministic(), "{report}");
}

/// Scratch reuse must not change the answer — the whole reason [`Welder`] owns
/// buffers is to be called thousands of times.
#[test]
fn a_reused_welder_gives_the_same_answer_as_a_fresh_one() {
    let h = 4.0 / 35.0;
    let layout = ChunkLayout::<f64>::new(8, h, [-2.0; 3]).expect("valid layout");
    let field = Sphere::<f64>::canonical();
    let eps = epsilon_for(h);

    let build = || {
        let a = ChunkId::new([1, 1, 1]);
        let mut m = mesh_chunk(&layout, &field, a);
        m.append(&mesh_chunk(&layout, &field, a.neighbour(0, 1)))
            .expect("the meshes fit the u32 index space");
        m
    };

    let mut reused = Welder::<f64>::new();
    // Warm the scratch on a different, larger mesh first.
    let mut other = build();
    other
        .append(&build())
        .expect("the meshes fit the u32 index space");
    reused.weld(&mut other, eps).expect("valid epsilon");

    let mut from_reused = build();
    let a = reused.weld(&mut from_reused, eps).expect("valid epsilon");
    let mut from_fresh = build();
    let b = Welder::<f64>::new()
        .weld(&mut from_fresh, eps)
        .expect("valid epsilon");

    assert_eq!(a, b);
    assert_eq!(from_reused.positions, from_fresh.positions);
    assert_eq!(from_reused.indices, from_fresh.indices);
}

// ─── against the validator ──────────────────────────────────────────────────

/// **`duplicate_vertices` is an upper bound on what a weld removes, not the
/// count.**
///
/// The validator asks whether *any* earlier vertex is within `ε`; the welder
/// asks for the lowest-indexed *kept* one. They differ exactly where a chain of
/// near-misses exists, because the validator counts the middle of a chain as a
/// duplicate of its start and the welder does not weld the end onto anything.
///
/// Predicted before running: **equal on a real seam**, where duplicates come in
/// pairs an ulp apart with no chains, and **different on the constructed chain**
/// above. Recorded either way.
#[test]
fn the_validator_bounds_the_weld_rather_than_predicting_it() {
    // A real seam: pairs, no chains.
    let h = 4.0 / 35.0;
    let layout = ChunkLayout::<f64>::new(8, h, [-2.0; 3]).expect("valid layout");
    let field = Sphere::<f64>::canonical();
    let a = ChunkId::new([1, 1, 1]);
    let mut joined = mesh_chunk(&layout, &field, a);
    joined
        .append(&mesh_chunk(&layout, &field, a.neighbour(0, 1)))
        .expect("the meshes fit the u32 index space");

    let cfg = ValidateConfig::from_cell_size(h).expect("valid cell size");
    let counted = validate_indexed(&joined.positions, &joined.indices, &cfg).duplicate_vertices;
    let removed = Welder::<f64>::new()
        .weld(&mut joined, epsilon_for(h))
        .expect("valid epsilon")
        .vertices_removed();
    std::println!("measured: seam -- validator counted {counted}, weld removed {removed}");
    assert_eq!(counted as usize, removed, "no chains at a seam");

    // The constructed chain, where they must differ.
    let mut chain = points(&[[0.0; 3], [0.9, 0.0, 0.0], [1.8, 0.0, 0.0]]);
    let cfg = ValidateConfig::from_cell_size(1.0 / ValidateConfig::WELD_EPSILON_REL)
        .expect("valid cell size");
    let counted = validate_indexed(&chain.positions, &chain.indices, &cfg).duplicate_vertices;
    let removed = Welder::<f64>::new()
        .weld(&mut chain, 1.0)
        .expect("valid epsilon")
        .vertices_removed();
    std::println!("measured: chain -- validator counted {counted}, weld removed {removed}");
    assert_eq!(counted, 2, "the validator counts 1 and 2 as duplicates");
    assert_eq!(removed, 1, "the weld only removes 1");
}

/// **A whole-volume mesh is not always weld-free, and the exception is exactly
/// the sliver case A-001 already recorded (M-48).**
///
/// The edge-vertex cache shares a vertex between the cells that meet on a grid
/// *edge*, and that is all it can do. When a grid **sample** lands on the
/// isosurface, `t` is 0 or 1 and the crossing sits *at that sample* — so every
/// cut edge meeting there places its own vertex at the same point, and they are
/// on different edges, so nothing shares them. Welding merges them, and the
/// triangles that used two of them collapse.
///
/// On `sphere` at 25³ that removes **48 vertices and 96 triangles**, and the 96
/// is not a coincidence: it is exactly the degenerate-sliver count A-001
/// measured at that resolution from the 30 lattice points that sit exactly on
/// the unit sphere. **Welding is therefore a fix for that class of sliver**,
/// which was not predicted.
///
/// Pinned in both directions, following M-4. A field not listed here welds to
/// nothing, and that is asserted too.
#[test]
fn a_whole_volume_mesh_welds_only_where_a_sample_sits_on_the_surface() {
    let mut found: Vec<(&str, u32, usize, usize)> = Vec::new();
    crate::for_each_reference_field!(f64, |name, field| {
        for samples in [17u32, 25, 33] {
            let (lo, hi) = field.domain();
            let h = (hi[0] - lo[0]) / f64::from(samples - 1);
            let shape = crate::RuntimeShape3::new([samples; 3]).expect("valid shape");
            let mut mesh = MeshBuffer::<f64>::new();
            MarchingCubes::<f64>::new()
                .extract(&field, &shape, lo, h, &mut mesh)
                .expect("extraction");
            if mesh.is_empty() {
                continue;
            }
            let report = Welder::<f64>::new()
                .weld(&mut mesh, epsilon_for(h))
                .expect("valid epsilon");
            if !report.is_noop() {
                found.push((
                    name,
                    samples,
                    report.vertices_removed(),
                    report.triangles_collapsed,
                ));
            }
        }
    });

    std::println!("measured: whole-volume welds (field, n, vertices, triangles) = {found:?}");
    assert_eq!(
        found,
        vec![
            ("sphere", 25, 48, 96),
            ("gyroid", 17, 2, 4),
            ("gyroid", 25, 2, 4),
            ("gyroid", 33, 2, 4),
            ("fbm_terrain", 33, 1, 2),
            // The eighth field welds far more than the others, and for the reason
            // M-212 records: volumetric noise sampled near its feature size puts many
            // crossings on or beside a grid point. One of these merges is the pair
            // that fuses two sheets (A-018); the rest are ordinary.
            ("noise_cavity", 17, 30, 54),
            ("noise_cavity", 25, 136, 266),
            ("noise_cavity", 33, 35, 68),
        ]
    );
}

/// **T-009's invariant: the welder and the validator mean the same thing by
/// "coincident".**
///
/// [`epsilon_for`] and
/// [`ValidateConfig::weld_epsilon`](crate::validate::ValidateConfig::weld_epsilon)
/// are two expressions of one policy, and if they ever drift the validator's
/// duplicate count is describing a different mesh than the welder produced —
/// which is precisely what `the_validator_bounds_the_weld_rather_than_predicting_it`
/// assumes and cannot check.
#[test]
fn the_one_policy_is_one_number_however_it_is_reached() {
    use crate::validate::ValidateConfig;

    // Spacings including M-32's non-power-of-two case, where bit-exactness fails
    // and the tolerance is actually doing work.
    for h in [1.0, 0.5, 0.125, 0.1, 1.0 / 3.0, 0.0625, 7.5] {
        let cfg = ValidateConfig::from_cell_size(h).expect("a valid spacing");
        assert_eq!(
            epsilon_for(h).to_bits(),
            cfg.weld_epsilon().to_bits(),
            "h = {h}: the welder and the validator disagree about coincidence"
        );
    }

    // And it is genuinely relative -- doubling the grid doubles the tolerance,
    // which is the property an absolute epsilon does not have.
    assert_eq!(
        epsilon_for(2.0f64).to_bits(),
        (2.0 * epsilon_for(1.0f64)).to_bits()
    );

    // f32 reaches the same policy through the same constant.
    assert_eq!(
        epsilon_for(1.0f32).to_bits(),
        (ValidateConfig::WELD_EPSILON_REL as f32).to_bits()
    );
}
