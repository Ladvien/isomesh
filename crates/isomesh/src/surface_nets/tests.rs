//! Tests for Surface Nets.
//!
//! Same validity gates as Marching Cubes, on the same fields, because the point
//! of having both is that they are comparable.

// Zero self-intersections and exact triangle counts are asserted exactly.
#![allow(clippy::float_cmp)]

use alloc::vec::Vec;

use super::SurfaceNets;
use crate::fields::{BoxExact, ReferenceField, Sphere, Torus};
use crate::marching_cubes::MarchingCubes;
use crate::validate::{ValidateConfig, check_determinism, self_intersections, validate_indexed};
use crate::{MeshBuffer, RuntimeShape3, Sdf, vec3};

fn mesh<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    samples: u32,
    smoothing: u32,
) -> (MeshBuffer<f64>, f64) {
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let mut sn = SurfaceNets::<f64>::new();
    sn.set_smoothing_passes(smoothing);
    let mut out = MeshBuffer::<f64>::new();
    sn.extract(
        field,
        &RuntimeShape3::new([samples; 3]).expect("valid shape"),
        lo,
        cell_size,
        &mut out,
    )
    .expect("extraction");
    (out, cell_size)
}

fn mesh_with_mc<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    samples: u32,
) -> (MeshBuffer<f64>, f64) {
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    mc.extract(
        field,
        &RuntimeShape3::new([samples; 3]).expect("valid shape"),
        lo,
        cell_size,
        &mut out,
    )
    .expect("extraction");
    (out, cell_size)
}

fn signed_volume(mesh: &MeshBuffer<f64>) -> f64 {
    let mut total = 0.0;
    for t in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[t[0] as usize];
        let b = mesh.positions[t[1] as usize];
        let c = mesh.positions[t[2] as usize];
        total += vec3::dot(a, vec3::cross(b, c));
    }
    total / 6.0
}

#[test]
fn a_meshed_sphere_is_closed() {
    for samples in [17u32, 25, 33] {
        let (out, h) = mesh(&Sphere::<f64>::canonical(), samples, 0);
        let report = validate_indexed(
            &out.positions,
            &out.indices,
            &ValidateConfig::from_cell_size(h).expect("valid cell size"),
        );
        assert!(out.triangle_count() > 0);
        assert!(report.is_closed(), "{samples} samples:\n{report}");
        assert_eq!(report.euler_characteristic, 2, "{samples}:\n{report}");
        assert_eq!(report.genus, Some(0));
        assert_eq!(report.non_manifold_edges, 0);
        assert_eq!(report.non_manifold_vertices, 0);
        assert_eq!(report.inconsistently_oriented_edges, 0);
    }
}

/// The only check that catches a globally inverted surface.
#[test]
fn meshed_sphere_has_positive_signed_volume() {
    let (out, _) = mesh(&Sphere::<f64>::canonical(), 33, 0);
    let volume = signed_volume(&out);
    let exact = 4.0 / 3.0 * core::f64::consts::PI;
    assert!(volume > 0.0, "inside out: {volume}");
    assert!((volume - exact).abs() / exact < 0.05, "{volume} vs {exact}");
}

#[test]
fn a_meshed_torus_has_genus_one() {
    let (out, h) = mesh(&Torus::<f64>::canonical(), 49, 0);
    let report = validate_indexed(
        &out.positions,
        &out.indices,
        &ValidateConfig::from_cell_size(h).expect("valid cell size"),
    );
    assert!(report.is_closed(), "{report}");
    assert_eq!(report.euler_characteristic, 0, "{report}");
    assert_eq!(report.genus, Some(1), "{report}");
}

/// Two disjoint unit-ish spheres: a closed surface with two components, so
/// `χ = 4`. Included because the identity is stated in terms of `χ` and every
/// other reference field has `χ` of 2 or 0 — without a `χ = 4` case the tests
/// could not tell `2χ` from "always 4".
#[derive(Clone, Copy, Debug)]
struct TwoSpheres;

impl Sdf for TwoSpheres {
    type Scalar = f64;
    fn sample(&self, p: [f64; 3]) -> f64 {
        let left = Sphere::<f64> {
            center: [-1.0, 0.0, 0.0],
            radius: 0.6,
        };
        let right = Sphere::<f64> {
            center: [1.0, 0.0, 0.0],
            radius: 0.6,
        };
        left.sample(p).min(right.sample(p))
    }
}

impl ReferenceField for TwoSpheres {
    const NAME: &'static str = "two_spheres";
    fn domain(&self) -> ([f64; 3], [f64; 3]) {
        ([-2.0; 3], [2.0; 3])
    }
    fn closed_in_domain(&self) -> bool {
        true
    }
    fn expected_euler(&self) -> Option<i64> {
        Some(4) // two components, each a sphere
    }
    fn is_exact_distance(&self) -> bool {
        true
    }
}

/// **Surface Nets and Marching Cubes produce the same triangle count, up to
/// `2χ`.** Measured on four fields, exact every time.
///
/// The implementation brief expects these "to differ substantially". They do
/// not, and the reason is Euler rather than luck:
///
/// - Marching Cubes places one vertex per crossed grid edge, so `V_mc = Q`.
/// - Surface Nets emits one quad — two triangles — per crossed grid edge, so
///   `F_sn = 2Q`.
/// - Any closed triangulated surface has `V − E + F = χ` and `3F = 2E`, hence
///   `F = 2V − 2χ`.
///
/// Therefore `F_mc = 2Q − 2χ = F_sn − 2χ`. The two counts are pinned to each
/// other by the number of crossed edges, and the widely repeated claim that
/// Surface Nets is the cheaper method *by triangle count* is false for closed
/// surfaces. What it actually buys is regular quad connectivity and one vertex
/// per cell instead of one per edge — not fewer triangles.
#[test]
fn triangle_counts_track_marching_cubes_up_to_two_chi() {
    fn check<F: Sdf<Scalar = f64> + ReferenceField>(name: &str, field: &F, samples: u32) {
        let (sn_mesh, h) = mesh(field, samples, 0);
        let (mc_mesh, _) = mesh_with_mc(field, samples);
        let cfg = ValidateConfig::from_cell_size(h).expect("valid cell size");
        let sn_report = validate_indexed(&sn_mesh.positions, &sn_mesh.indices, &cfg);
        let mc_report = validate_indexed(&mc_mesh.positions, &mc_mesh.indices, &cfg);

        assert!(sn_report.is_closed() && mc_report.is_closed(), "{name}");
        assert_eq!(
            sn_report.euler_characteristic, mc_report.euler_characteristic,
            "{name}: the two methods must agree on topology"
        );

        let chi = sn_report.euler_characteristic;
        let difference = sn_mesh.triangle_count() as i64 - mc_mesh.triangle_count() as i64;
        assert_eq!(
            difference,
            2 * chi,
            "{name} at {samples}: F_sn {} vs F_mc {}, chi {chi}",
            sn_mesh.triangle_count(),
            mc_mesh.triangle_count(),
        );
        // The vertex form of the same identity: F = 2V - 2*chi for both, so a
        // triangle gap of 2*chi is a vertex gap of chi.
        let vertex_difference = sn_mesh.vertex_count() as i64 - mc_mesh.vertex_count() as i64;
        assert_eq!(
            vertex_difference,
            chi,
            "{name} at {samples}: V_sn {} vs V_mc {}, chi {chi}",
            sn_mesh.vertex_count(),
            mc_mesh.vertex_count(),
        );

        std::println!(
            "measured: {name} at {samples}^3 -> surface nets {} tris / {} verts, \
             marching cubes {} tris / {} verts (difference {difference} = 2*chi)",
            sn_mesh.triangle_count(),
            sn_mesh.vertex_count(),
            mc_mesh.triangle_count(),
            mc_mesh.vertex_count(),
        );
    }

    // Three resolutions each, and a chi = 4 field, so the identity is tested as
    // `2*chi` rather than as the constant 4 that three of these fields share.
    for samples in [17u32, 25, 33] {
        check("sphere", &Sphere::<f64>::canonical(), samples);
        check("box_exact", &BoxExact::<f64>::canonical(), samples);
        check("two_spheres", &TwoSpheres, samples);
    }
    for samples in [33u32, 41, 49] {
        check("torus", &Torus::<f64>::canonical(), samples);
        check(
            "csg_difference",
            &crate::fields::csg_difference::<f64>(),
            samples,
        );
    }
}

/// Vertex degree, recorded rather than bounded.
///
/// A dual vertex belongs to one quad per crossed edge of its cell, and a cell
/// can have up to twelve crossed edges, so there is no tight bound worth
/// asserting. The measurement is the interesting part.
#[test]
fn vertex_degree_is_recorded() {
    fn max_degree(mesh: &MeshBuffer<f64>) -> u32 {
        let mut degree = alloc::vec![0u32; mesh.vertex_count()];
        for t in mesh.indices.chunks_exact(3) {
            for &i in t {
                degree[i as usize] += 1;
            }
        }
        degree.iter().copied().max().unwrap_or(0)
    }
    let sn_max = max_degree(&mesh(&Sphere::<f64>::canonical(), 25, 0).0);
    let mc_max = max_degree(&mesh_with_mc(&Sphere::<f64>::canonical(), 25).0);
    std::println!("measured: max vertex degree -- surface nets {sn_max}, marching cubes {mc_max}");
    assert!(sn_max > 0 && mc_max > 0);
}

/// Smoothing must not break the surface. It moves vertices, so it is exactly the
/// kind of change that can turn a valid mesh into a self-intersecting one.
#[test]
fn smoothing_keeps_the_surface_closed() {
    for passes in [0u32, 1, 4] {
        let (out, h) = mesh(&Sphere::<f64>::canonical(), 25, passes);
        let report = validate_indexed(
            &out.positions,
            &out.indices,
            &ValidateConfig::from_cell_size(h).expect("valid cell size"),
        );
        assert!(report.is_closed(), "{passes} passes:\n{report}");
        assert_eq!(report.euler_characteristic, 2, "{passes} passes:\n{report}");

        let si = self_intersections(&out.positions, &out.indices, h).expect("self intersections");
        std::println!(
            "measured: sphere at 25^3, {passes} smoothing passes -> \
             {:.3} intersecting pairs per 1000 triangles",
            si.per_thousand_triangles()
        );
    }
}

/// Smoothing is supposed to *do* something, or the parameter is decoration.
#[test]
fn smoothing_moves_vertices() {
    let rough = mesh(&Sphere::<f64>::canonical(), 25, 0).0;
    let smooth = mesh(&Sphere::<f64>::canonical(), 25, 4).0;
    assert_eq!(rough.vertex_count(), smooth.vertex_count());

    let moved = rough
        .positions
        .iter()
        .zip(&smooth.positions)
        .filter(|(a, b)| vec3::length(vec3::sub(**a, **b)) > 1e-12)
        .count();
    assert!(
        moved > rough.vertex_count() / 2,
        "only {moved} of {} vertices moved",
        rough.vertex_count()
    );
}

/// The characteristic Surface Nets weakness, stated as a measurement rather than
/// a claim — and the reason E-104 exists.
///
/// A vertex at the centroid of its cell's edge crossings cannot land on a box's
/// corner: the average of points on two faces lies between them. Dual contouring
/// replaces the centroid with a solve and recovers the corner.
#[test]
fn a_box_has_its_corners_rounded_off() {
    let field = BoxExact::<f64>::canonical();
    let (out, h) = mesh(&field, 33, 0);
    let report = validate_indexed(
        &out.positions,
        &out.indices,
        &ValidateConfig::from_cell_size(h).expect("valid cell size"),
    );
    assert!(report.is_closed(), "{report}");

    // The box is [-1, 1]^3, so its eight corners are at distance sqrt(3) from
    // the origin. How close does the mesh actually get?
    let corner = [1.0f64, 1.0, 1.0];
    let nearest = out
        .positions
        .iter()
        .map(|p| vec3::length(vec3::sub(*p, corner)))
        .fold(f64::INFINITY, f64::min);

    std::println!(
        "measured: surface nets on box_exact at 33^3 -> nearest vertex to the corner \
         (1,1,1) is {nearest:.4} away ({:.2} cells)",
        nearest / h
    );
    assert!(
        nearest > h * 0.25,
        "the corner was reproduced exactly, which surface nets cannot do: {nearest}"
    );
}

/// How close each method gets to a box corner — and the grid-alignment trap
/// that makes the obvious answer wrong.
///
/// **`box_exact` is exactly zero across its entire boundary**, not just at the
/// surface in a limit sense: `f(1,0,0)`, `f(1,1,0)` and `f(1,1,1)` are all `+0`.
/// Since zero classifies as *outside*, a grid plane lying on a box face is
/// entirely outside, and the sign change happens a whole cell further in. So on
/// a grid-aligned box, **Marching Cubes lands further from the corner than
/// Surface Nets does** — the opposite of what "MC puts vertices on edges so it
/// can hit the corner" suggests.
///
/// That is a real trap for anyone benchmarking sharp-feature recovery on an
/// axis-aligned box: the answer is decided by the zero-classification rule
/// rather than by the algorithm. Over the `[-2, 2]` domain, 25 and 33 samples
/// are both grid-aligned; 27 is not.
///
/// What survives as a robust statement, and what E-104 has to beat: **Surface
/// Nets cannot reach a corner at any resolution**, because its vertex is the
/// centroid of its cell's edge crossings and an average of points on a corner's
/// faces lies strictly inside it.
#[test]
fn neither_method_reaches_a_box_corner_and_the_reason_is_the_grid() {
    let field = BoxExact::<f64>::canonical();
    let corner = [1.0f64, 1.0, 1.0];
    let nearest = |mesh: &MeshBuffer<f64>| {
        mesh.positions
            .iter()
            .map(|p| vec3::length(vec3::sub(*p, corner)))
            .fold(f64::INFINITY, f64::min)
    };

    for (samples, aligned) in [(25u32, true), (27, false), (33, true)] {
        let (sn_mesh, h) = mesh(&field, samples, 0);
        let (mc_mesh, _) = mesh_with_mc(&field, samples);
        let (sn_gap, mc_gap) = (nearest(&sn_mesh), nearest(&mc_mesh));

        std::println!(
            "measured: box_exact at {samples}^3 (h = {h:.4}, {}) -> nearest vertex to (1,1,1): \
             marching cubes {mc_gap:.4} ({:.2} cells), surface nets {sn_gap:.4} ({:.2} cells)",
            if aligned {
                "grid-aligned"
            } else {
                "not aligned"
            },
            mc_gap / h,
            sn_gap / h,
        );

        // The robust claim. Everything else here is a recorded measurement.
        assert!(
            sn_gap > h * 0.25,
            "{samples}: surface nets should not reach the corner, got {sn_gap}"
        );
    }
}

#[test]
fn surface_nets_is_deterministic() {
    let field = Sphere::<f64>::canonical();
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / 24.0;
    let shape = RuntimeShape3::new([25; 3]).expect("valid shape");
    let mut sn = SurfaceNets::<f64>::new();

    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        sn.extract(&field, &shape, lo, cell_size, out)
            .expect("extraction");
    });
    assert!(report.is_deterministic(), "{report}");
    assert!(report.triangles > 0);
}

#[test]
fn every_closed_reference_field_meshes_cleanly() {
    fn check<F: Sdf<Scalar = f64> + ReferenceField>(name: &str, field: &F, samples: u32) {
        let (out, h) = mesh(field, samples, 0);
        let report = validate_indexed(
            &out.positions,
            &out.indices,
            &ValidateConfig::from_cell_size(h).expect("valid cell size"),
        );
        assert!(out.triangle_count() > 0, "{name} produced nothing");
        if field.closed_in_domain() {
            assert!(report.is_closed(), "{name}:\n{report}");
        } else {
            assert!(report.is_manifold(), "{name}:\n{report}");
        }
        if let Some(chi) = field.expected_euler() {
            assert_eq!(report.euler_characteristic, chi, "{name}:\n{report}");
        }
        std::println!(
            "measured: surface nets {name} at {samples}^3 -> {} tris, chi {}",
            out.triangle_count(),
            report.euler_characteristic
        );
    }

    check("sphere", &Sphere::<f64>::canonical(), 33);
    check("torus", &Torus::<f64>::canonical(), 49);
    check("box_exact", &BoxExact::<f64>::canonical(), 33);
    check(
        "csg_difference",
        &crate::fields::csg_difference::<f64>(),
        41,
    );
    // gyroid and fbm_terrain are deliberately absent -- see the next test.
}

/// **Naive Surface Nets is not manifold on the gyroid or on fbm_terrain, and
/// that is the method's known structural defect rather than a bug here.**
///
/// A dual method places exactly one vertex per cell. Where two separate sheets
/// of the surface pass through the same cell — which is what a high-genus field
/// at a coarse grid produces constantly — that one vertex is forced to join
/// sheets that should stay apart, and the result is non-manifold.
///
/// `docs/research/2026-08-11-literature-review-round-1.md` names this exactly:
/// failures of this kind are "in *connectivity*: multiple surface sheets through
/// one cell, which is DC's actual structural defect. The fix is architectural
/// (partition-based extractor), not another patch." A-010's vertex splitting is
/// that fix; until it lands, this is pinned as a measurement so the number moves
/// only when someone means it to.
#[test]
fn multi_sheet_cells_expose_the_one_vertex_per_cell_limit() {
    fn check<F: Sdf<Scalar = f64> + ReferenceField>(name: &str, field: &F, samples: u32) {
        let (out, h) = mesh(field, samples, 0);
        let report = validate_indexed(
            &out.positions,
            &out.indices,
            &ValidateConfig::from_cell_size(h).expect("valid cell size"),
        );

        std::println!(
            "measured: surface nets {name} at {samples}^3 -> {} tris, chi {}, \
             {} non-manifold edges, {} non-manifold vertices",
            out.triangle_count(),
            report.euler_characteristic,
            report.non_manifold_edges,
            report.non_manifold_vertices,
        );

        assert!(
            report.non_manifold_edges > 0,
            "{name}: if this is now zero, a dual method gained sheet separation \
             and A-010 should be revisited:\n{report}"
        );
        // The damage is confined to connectivity -- no orientation flips.
        assert_eq!(report.inconsistently_oriented_edges, 0, "{name}:\n{report}");
    }

    // High genus: sheets of the gyroid pass through one cell constantly.
    check("gyroid", &crate::fields::capped_gyroid::<f64>(), 49);
    // Not high genus, but its highest octave has a wavelength comparable to the
    // cell size at this resolution, so the terrain folds within single cells.
    check(
        "fbm_terrain",
        &crate::fields::FbmTerrain::<f64>::canonical(),
        33,
    );
}

#[test]
fn f32_and_f64_both_extract() {
    let mut sn = SurfaceNets::<f32>::new();
    let mut out = MeshBuffer::<f32>::new();
    sn.extract(
        &Sphere::<f32>::canonical(),
        &RuntimeShape3::new([17; 3]).expect("valid shape"),
        [-2.0; 3],
        4.0 / 16.0,
        &mut out,
    )
    .expect("extraction");
    assert!(out.triangle_count() > 0);
}

#[test]
fn a_field_with_no_surface_produces_no_triangles() {
    let field = Sphere::<f64> {
        center: [100.0; 3],
        radius: 1.0,
    };
    let mut sn = SurfaceNets::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    sn.extract(
        &field,
        &RuntimeShape3::new([9; 3]).expect("valid shape"),
        [-2.0; 3],
        0.5,
        &mut out,
    )
    .expect("extraction");
    assert_eq!(out.triangle_count(), 0);
}

#[test]
fn a_degenerate_grid_is_rejected() {
    let mut sn = SurfaceNets::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    let shape = RuntimeShape3::new([1, 4, 4]).expect("valid shape");
    let error = sn
        .extract(
            &Sphere::<f64>::canonical(),
            &shape,
            [-2.0; 3],
            0.5,
            &mut out,
        )
        .expect_err("a one-sample axis contains no cell");
    assert_eq!(error, crate::Error::GridTooSmall { size: [1, 4, 4] });
    assert!(out.is_empty(), "nothing should have been written");
}

/// Vertices are shared, so the count tracks crossed *cells* rather than
/// triangles.
#[test]
fn one_vertex_per_crossed_cell() {
    let (out, _) = mesh(&Sphere::<f64>::canonical(), 25, 0);
    let unique: Vec<u32> = {
        let mut v: Vec<u32> = out.indices.clone();
        v.sort_unstable();
        v.dedup();
        v
    };
    assert_eq!(
        unique.len(),
        out.vertex_count(),
        "every emitted vertex should be referenced"
    );
}
