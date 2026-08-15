//! The load-bearing tests here are the topology identity against Surface Nets —
//! which is V-19 stated as an assertion — and the corner measurement, which is
//! the entire reason this algorithm exists.

use super::DualContouring;
use crate::fields::{BoxExact, ReferenceField, Sphere};
use crate::surface_nets::SurfaceNets;
use crate::validate::{ValidateConfig, check_determinism, validate_indexed};
use crate::{MeshBuffer, RuntimeShape3, Sdf};

/// Mesh a reference field at `samples` per axis, the same convention the rest of
/// the suite uses: `shape` counts samples, so `n` spans `n - 1` cells.
///
/// Uses the shipped default, which includes the cell clamp.
fn mesh_dc<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    samples: u32,
) -> (MeshBuffer<f64>, f64) {
    mesh_dc_with(field, samples, crate::dual_contouring::Clamp::ToCell)
}

/// As [`mesh_dc`], with the clamp chosen explicitly.
///
/// Tests that characterise the *vertex rule* use [`Clamp::None`], because the
/// clamp box is the cell shrunk by `CLAMP_EPSILON` and therefore nudges any
/// vertex sitting exactly *on* a cell face by `5e-5` cells. That is nothing
/// geometrically and it is fatal to a bit-comparison — and on a **grid-aligned**
/// `box_exact` it is most of the mesh, since the box faces coincide with grid
/// planes. Over the ±2 domain a grid is aligned when `n - 1` is a multiple of 4,
/// which 17 and 33 are.
fn mesh_dc_with<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    samples: u32,
    clamp: crate::dual_contouring::Clamp,
) -> (MeshBuffer<f64>, f64) {
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let mut dc = DualContouring::<f64>::new();
    dc.set_clamp(clamp);
    let mut out = MeshBuffer::<f64>::new();
    dc.extract(field, &shape, lo, cell_size, &mut out)
        .expect("extraction");
    (out, cell_size)
}

fn mesh_sn<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    samples: u32,
) -> (MeshBuffer<f64>, f64) {
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let mut sn = SurfaceNets::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    sn.extract(field, &shape, lo, cell_size, &mut out)
        .expect("extraction");
    (out, cell_size)
}

/// Distance from `target` to the nearest vertex in the mesh.
fn nearest_vertex(mesh: &MeshBuffer<f64>, target: [f64; 3]) -> f64 {
    mesh.positions
        .iter()
        .map(|p| {
            ((p[0] - target[0]).powi(2) + (p[1] - target[1]).powi(2) + (p[2] - target[2]).powi(2))
                .sqrt()
        })
        .fold(f64::INFINITY, f64::min)
}

/// **V-19, as an assertion.** Dual Contouring's topology *is* Surface Nets'
/// topology — one vertex per crossed cell, one quad per crossed edge — so on the
/// same field and grid the two must agree on every index, and differ only in
/// where the vertices sit.
///
/// This is what the shared engine in `crate::dual` buys, and it is the reason
/// E-104 is an honest comparison rather than two algorithms shrugging at each
/// other: same crossings, same connectivity, same winding, one function swapped.
#[test]
fn topology_is_identical_to_surface_nets() {
    for samples in [17u32, 27, 33] {
        let field = BoxExact::<f64>::canonical();
        let (dc, _) = mesh_dc_with(&field, samples, crate::dual_contouring::Clamp::None);
        let (sn, _) = mesh_sn(&field, samples);

        assert_eq!(
            dc.vertex_count(),
            sn.vertex_count(),
            "{samples}^3: vertex counts differ"
        );
        assert_eq!(
            dc.indices, sn.indices,
            "{samples}^3: connectivity differs — the shared engine is not shared"
        );
        // And the vertices genuinely do move, or the rule is not being used --
        // but only some of them, which is the subject of the next test.
        let moved = dc
            .positions
            .iter()
            .zip(&sn.positions)
            .filter(|(a, b)| (0..3).any(|k| a[k].to_bits() != b[k].to_bits()))
            .count();
        assert!(moved > 0, "{samples}^3: no vertex moved at all");
        assert!(
            moved < dc.vertex_count(),
            "{samples}^3: every vertex moved, so the flat-region identity is broken"
        );
    }
}

/// **Dual Contouring agrees with Surface Nets everywhere the surface is flat,
/// and differs only at features** — and the split is not a matter of degree.
///
/// The reason is exact: on a planar patch every crossing lies in the plane, so
/// the centroid does too, so `pᵢ − c` lies *in* the plane and is perpendicular to
/// `n`. Every `dᵢ` is exactly zero, `g` is exactly zero, and the solve returns
/// `x = c + adj(A)·0/det(A) = c`. Mathematically it *is* the Surface Nets vertex.
///
/// It is not bit-identical, and that is worth knowing rather than glossing: the
/// two centroids are computed by different expressions — Surface Nets
/// accumulates offsets in cell-local units and scales once,
/// [`HermiteCell`](crate::hermite::HermiteCell) works in world coordinates — so
/// they agree to rounding, not to the last bit.
///
/// What makes this an assertion rather than an observation is the **gap**.
/// Measured on `box_exact` at 27³: 864 of 1016 vertices agree to within
/// `2e-15` cells, 152 move by `0.35`–`0.57` cells, and **nothing lands in
/// between** — fourteen orders of magnitude of empty space. A vertex is either
/// on a flat face or on a feature; there is no continuum.
///
/// The consequence is the one that matters for E-104: the two methods differ
/// *only* where dual contouring has something to contribute, so a side-by-side
/// comparison is measuring the feature and nothing else.
#[test]
fn dual_contouring_moves_only_the_feature_vertices() {
    let field = BoxExact::<f64>::canonical();
    let samples = 27;
    let (dc, cell_size) = mesh_dc_with(&field, samples, crate::dual_contouring::Clamp::None);
    let (sn, _) = mesh_sn(&field, samples);

    /// Below this (in cells) two vertices are the same point computed twice.
    const ROUNDING: f64 = 1e-12;
    /// Above this (in cells) a vertex has genuinely been moved by the solve.
    const FEATURE: f64 = 0.1;

    let mut flat = 0usize;
    let mut feature = 0usize;
    let mut between = 0usize;
    for (a, b) in dc.positions.iter().zip(&sn.positions) {
        let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
            / cell_size;
        if d <= ROUNDING {
            flat += 1;
        } else if d >= FEATURE {
            feature += 1;
        } else {
            between += 1;
        }
    }
    std::println!(
        "measured: box_exact at {samples}^3 -> {flat} vertices agree with surface nets to rounding, {feature} moved by >= {FEATURE} cells, {between} in between"
    );

    assert_eq!(
        between, 0,
        "the split must be clean: {between} vertices landed between rounding and a feature"
    );
    assert!(feature > 0, "no vertex was moved by the solve");
    // A box is mostly flat, so the flat population must dominate.
    assert!(
        flat > 4 * feature,
        "a box is mostly flat: {flat} flat vs {feature} feature"
    );
}

/// **The measurement the algorithm exists for.**
///
/// Surface Nets averages a cell's crossings, so its vertex lands between a
/// corner's two faces and never on it. Dual Contouring solves for where the
/// tangent planes meet, which is the corner.
///
/// The resolution is **27³ and not 25³ or 33³**, and that is load-bearing.
/// `box_exact` is exactly zero across its whole boundary, so on a grid whose
/// planes land on the box faces the sign classification decides the answer
/// rather than the algorithm — E-103 measured Marching Cubes coming out *worse*
/// than Surface Nets on an aligned grid for exactly that reason. Over the ±2
/// domain a grid is aligned when `n - 1` is a multiple of 4; 27 is not.
#[test]
fn the_corner_is_sharper_than_surface_nets() {
    let field = BoxExact::<f64>::canonical();
    let samples = 27;
    let (dc, cell_size) = mesh_dc(&field, samples);
    let (sn, _) = mesh_sn(&field, samples);

    let corner = [1.0, 1.0, 1.0];
    let dc_gap = nearest_vertex(&dc, corner);
    let sn_gap = nearest_vertex(&sn, corner);

    std::println!(
        "measured: box_exact at {samples}^3 (h = {cell_size:.4}, not grid-aligned) -> nearest vertex to (1,1,1): dual contouring {dc_gap:.4} ({:.2} cells), surface nets {sn_gap:.4} ({:.2} cells)",
        dc_gap / cell_size,
        sn_gap / cell_size
    );

    assert!(
        dc_gap < sn_gap,
        "dual contouring must reach the corner Surface Nets rounds: {dc_gap} vs {sn_gap}"
    );
    // Within a fifth of a cell of the true corner. Surface Nets measures 0.58
    // cells here, so this is a real gap rather than a rounding difference.
    assert!(
        dc_gap < 0.2 * cell_size,
        "expected the corner within 0.2 cells, got {:.3} cells",
        dc_gap / cell_size
    );
}

/// Every reference field meshes into something structurally sound.
///
/// The gate is **not** `is_closed`, and the reason is measured rather than
/// assumed: Dual Contouring inherits one-vertex-per-cell from the shared dual
/// topology, so where two sheets share a cell the mesh is non-manifold by
/// construction — M-4 and M-15 on Surface Nets, and A-010 is the ticket that
/// fixes it. What is asserted is everything that is *not* a consequence of that:
/// no structural errors, no boundary on a closed field, consistent winding.
#[test]
fn every_reference_field_meshes_soundly() {
    crate::for_each_reference_field!(f64, |name, field| {
        for samples in [17u32, 27, 33] {
            let (mesh, cell_size) = mesh_dc(&field, samples);
            let cfg = ValidateConfig::from_cell_size(cell_size).expect("valid spacing");
            let report = validate_indexed(&mesh.positions, &mesh.indices, &cfg);

            assert!(
                !report.has_structural_errors(),
                "{name} at {samples}^3: malformed\n{report}"
            );
            assert_eq!(
                report.inconsistently_oriented_edges, 0,
                "{name} at {samples}^3: winding\n{report}"
            );
            if field.closed_in_domain() {
                assert_eq!(
                    report.boundary_edges, 0,
                    "{name} at {samples}^3: the surface left the grid\n{report}"
                );
            }
            std::println!(
                "measured: dual contouring {name:16} at {samples}^3 -> {} tris, chi {}, {} non-manifold edges",
                mesh.triangle_count(),
                report.euler_characteristic,
                report.non_manifold_edges
            );
        }
    });
}

/// T-004's harness, on this algorithm. The solve sorts by magnitude in several
/// places, and a sort with an unstable tie-break would show up here.
#[test]
fn extraction_is_deterministic() {
    let field = Sphere::<f64>::canonical();
    let (lo, hi) = field.domain();
    let samples = 25u32;
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        let mut dc = DualContouring::<f64>::new();
        dc.extract(&field, &shape, lo, cell_size, out)
            .expect("extraction");
    });
    assert!(report.is_deterministic(), "{report}");
}

/// The whole grid being outside is not an error, and neither is a grid too small
/// to hold a cell being rejected.
#[test]
fn edge_cases_are_handled_rather_than_panicking() {
    let field = Sphere::<f64>::canonical();
    let mut dc = DualContouring::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();

    // A one-sample axis has no cells.
    let flat = RuntimeShape3::new([1, 9, 9]).expect("valid shape");
    assert!(dc.extract(&field, &flat, [-2.0; 3], 0.5, &mut out).is_err());

    // A grid far outside the sphere: valid, and empty.
    let away = RuntimeShape3::new([5; 3]).expect("valid shape");
    out.reset();
    dc.extract(&field, &away, [50.0; 3], 0.1, &mut out)
        .expect("extraction");
    assert_eq!(out.triangle_count(), 0);
}

// ─── A-009: the cell clamp ──────────────────────────────────────────────────

/// **The experiment A-009 exists for**, on all seven reference fields.
///
/// Reports λ — intersecting pairs per 1,000 triangles — with the clamp off and
/// on, plus what the clamp costs in sharpness. λ and not `p`, the
/// fraction-of-meshes-with-any: `p = 1 − e^{−λT}` saturates with chunk size, so
/// a large enough chunk reads `p ≈ 1` for any λ > 0 and the statistic stops
/// meaning anything.
///
/// The literature review states the decision rule in advance, which is why this
/// is an experiment rather than a demonstration:
///
/// - **λ → 0** — the failures were unclamped placement. Fixed, at some
///   sharp-feature cost, and guaranteed intersection-free extraction is
///   available.
/// - **λ unchanged** — the failures are in *connectivity*: multiple sheets
///   through one cell, which is the structural defect A-010 fixes architecturally
///   rather than by patching.
#[test]
fn the_clamp_measured_on_every_reference_field() {
    use crate::dual_contouring::Clamp;
    use crate::validate::self_intersections;

    let samples = 33u32;
    std::println!(
        "measured: A-009 cell clamp at {samples}^3 -- lambda = intersecting pairs per 1,000 triangles"
    );
    std::println!(
        "  {:16} {:>8} {:>12} {:>12} {:>14} {:>9}",
        "field",
        "tris",
        "lambda off",
        "lambda on",
        "outside/verts",
        "max cells"
    );

    crate::for_each_reference_field!(f64, |name, field| {
        let (lo, hi) = field.domain();
        let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
        let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

        let measure = |clamp: Clamp| {
            let mut dc = DualContouring::<f64>::new();
            dc.set_clamp(clamp);
            let mut out = MeshBuffer::<f64>::new();
            dc.extract(&field, &shape, lo, cell_size, &mut out)
                .expect("extraction");
            let si = self_intersections(&out.positions, &out.indices, cell_size)
                .expect("self-intersection scan");
            (out, si)
        };

        let (unclamped, si_off) = measure(Clamp::None);
        let (clamped, si_on) = measure(Clamp::ToCell);

        // How far the clamp moved each vertex, in cells.
        //
        // Counted against a threshold rather than by bit-inequality, and the
        // difference is not cosmetic: the clamp box is the cell shrunk by
        // `CLAMP_EPSILON`, so a vertex sitting exactly *on* a cell face gets nudged
        // by `5e-5` cells and a bitwise count calls that "clamped". On a
        // grid-aligned `box_exact` that is most of the mesh, and the column then
        // reads as though the solve were flinging vertices everywhere when it is
        // doing nothing of the kind. What the guarantee is about is vertices that
        // were genuinely *outside*.
        let displacement: alloc::vec::Vec<f64> = unclamped
            .positions
            .iter()
            .zip(&clamped.positions)
            .map(|(a, b)| {
                ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
                    / cell_size
            })
            .collect();
        let outside = displacement.iter().filter(|d| **d > 1e-3).count();
        let worst = displacement.iter().copied().fold(0.0f64, f64::max);

        std::println!(
            "  {name:16} {:>8} {:>12.3} {:>12.3} {:>7}/{:<6} {:>9.2}",
            clamped.triangle_count(),
            si_off.per_thousand_triangles(),
            si_on.per_thousand_triangles(),
            outside,
            clamped.vertex_count(),
            worst
        );

        // The clamp must never make things worse, on any field.
        assert!(
            si_on.per_thousand_triangles() <= si_off.per_thousand_triangles(),
            "{name}: the clamp increased lambda, {} -> {}",
            si_off.per_thousand_triangles(),
            si_on.per_thousand_triangles()
        );
        // Topology is untouched by where a vertex sits.
        assert_eq!(
            unclamped.indices, clamped.indices,
            "{name}: the clamp changed connectivity"
        );
    });
}

/// What the clamp costs at the feature it exists to preserve.
///
/// The corner of `box_exact` is where the unclamped solve does its best work, so
/// it is where a constraint is most likely to hurt. Measured at 27³, off the
/// aligned grid, for E-103's reason.
#[test]
fn the_clamp_cost_in_sharpness_is_measured() {
    use crate::dual_contouring::Clamp;

    let field = BoxExact::<f64>::canonical();
    let samples = 27u32;
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

    let gap = |clamp: Clamp| {
        let mut dc = DualContouring::<f64>::new();
        dc.set_clamp(clamp);
        let mut out = MeshBuffer::<f64>::new();
        dc.extract(&field, &shape, lo, cell_size, &mut out)
            .expect("extraction");
        nearest_vertex(&out, [1.0, 1.0, 1.0]) / cell_size
    };

    let off = gap(Clamp::None);
    let on = gap(Clamp::ToCell);
    std::println!(
        "measured: A-009 sharpness cost on box_exact at {samples}^3 -> corner gap {off:.4} cells unclamped, {on:.4} cells clamped"
    );

    // Whatever it costs, the clamped result must still beat Surface Nets' 0.58
    // cells, or the clamp has given the feature back.
    assert!(
        on < 0.2,
        "the clamp must not surrender the corner: {on:.4} cells"
    );
}

/// **The two ablation arms differ in vertex positions and in nothing else
/// (X-002).**
///
/// The seam's correctness property, and a stronger statement than "it compiles".
/// `DualContouring<Qef>` and `DualContouring<Centroid>` must produce:
///
/// - **the same index buffer**, because the cell classification and the quad
///   walk are the same lines of code and the rule cannot reach them; and
/// - **different positions**, because otherwise the type parameter is not
///   swapping what it claims to and every ablation measured through it is
///   measuring nothing.
///
/// Both halves matter. Identical indices with identical positions is a seam that
/// does not work; differing indices is a seam that changes more than the rule,
/// which would make X-004's comparison of two vertex rules a comparison of two
/// algorithms.
#[test]
fn the_ablation_arms_differ_only_in_position() {
    use crate::surface_nets::Centroid;

    let field = crate::fields::Torus::<f64>::canonical();
    let shape = crate::RuntimeShape3::new([25; 3]).expect("valid shape");

    let mut qef_mesh = crate::MeshBuffer::<f64>::new();
    DualContouring::<f64>::new()
        .extract(&field, &shape, [-2.0; 3], 0.16, &mut qef_mesh)
        .expect("extraction");

    let mut centroid_mesh = crate::MeshBuffer::<f64>::new();
    DualContouring::<f64, Centroid>::with_rule(Centroid)
        .extract(&field, &shape, [-2.0; 3], 0.16, &mut centroid_mesh)
        .expect("extraction");

    assert_eq!(
        qef_mesh.indices, centroid_mesh.indices,
        "the vertex rule changed the topology, so it is not only the rule that swapped"
    );
    assert_eq!(qef_mesh.positions.len(), centroid_mesh.positions.len());
    assert!(
        qef_mesh.positions.len() > 100,
        "too small a mesh for this to mean anything"
    );

    // Exact comparison on purpose: the question is whether the rule moved the
    // vertex at all, and a tolerance would answer a different one. Two rules that
    // agreed to within an epsilon would still be two rules.
    #[allow(clippy::float_cmp)]
    let moved = qef_mesh
        .positions
        .iter()
        .zip(&centroid_mesh.positions)
        .filter(|(a, b)| a != b)
        .count();
    assert!(
        moved > qef_mesh.positions.len() / 2,
        "only {moved} of {} vertices moved — the rule is not being applied",
        qef_mesh.positions.len()
    );
    std::println!(
        "measured: {} vertices, identical indices, {moved} positions differ",
        qef_mesh.positions.len()
    );
}

/// **The seam is static dispatch, so there is no branch to emit (X-002).**
///
/// The acceptance asks for "no `if` in the hot loop". That is a type-system
/// property rather than a measurement: a generic parameter bound by a trait is
/// monomorphised, so each arm is compiled separately and there is no value for
/// the placement to test. The only way to lose it is to reach for `dyn`, which
/// would put a vtable call in exactly the loop this exists to keep clean.
///
/// So what is checkable is that nobody has: the dual path names no trait object.
#[test]
fn the_ablation_arms_are_not_branches() {
    for (name, source) in [
        ("dual.rs", include_str!("../dual.rs")),
        ("dual_contouring.rs", include_str!("../dual_contouring.rs")),
        ("surface_nets.rs", include_str!("../surface_nets.rs")),
        (
            "manifold_dual_contouring.rs",
            include_str!("../manifold_dual_contouring.rs"),
        ),
    ] {
        assert!(
            !source.contains("dyn VertexRule"),
            "{name} dispatches the vertex rule dynamically, which puts a vtable call \
             in the per-cell loop and gives up what the type parameter buys"
        );
    }
}
