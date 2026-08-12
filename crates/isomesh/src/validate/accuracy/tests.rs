//! The load-bearing tests here are the closed-form fixtures and the
//! brute-force cross-check.
//!
//! The ticket's own acceptance criterion — a unit sphere at 64³ within one cell
//! diagonal — passes with roughly sixty times to spare, so on its own it cannot
//! tell a correct harness from a badly broken one. What actually pins this
//! module down is: three polyhedra whose every quantity is computable by hand,
//! a spatial index checked bit-for-bit against exhaustive search, and a
//! convergence-order test that a constant-returning harness cannot fake.

use alloc::vec::Vec;

use super::super::tri_grid::{TriangleGrid, point_triangle_distance_squared};
use super::{AccuracyConfig, accuracy};
use crate::fields::{BoxExact, ReferenceField, Sphere};
use crate::marching_cubes::MarchingCubes;
use crate::surface_nets::SurfaceNets;
use crate::{MeshBuffer, RuntimeShape3, Sdf};

// ─── fixtures ───────────────────────────────────────────────────────────────

/// A regular octahedron of circumradius `s`, wound counter-clockwise from
/// outside. Vertex `2k` is `+s` on axis `k`, vertex `2k+1` is `−s`.
fn octahedron(s: f64) -> (Vec<[f64; 3]>, Vec<u32>) {
    let positions = alloc::vec![
        [s, 0.0, 0.0],
        [-s, 0.0, 0.0],
        [0.0, s, 0.0],
        [0.0, -s, 0.0],
        [0.0, 0.0, s],
        [0.0, 0.0, -s],
    ];
    let indices = alloc::vec![
        0, 2, 4, //
        0, 5, 2, //
        0, 4, 3, //
        0, 3, 5, //
        1, 4, 2, //
        1, 2, 5, //
        1, 3, 4, //
        1, 5, 3, //
    ];
    (positions, indices)
}

/// The lattice every polyhedron fixture is measured on: spacing 1, origin
/// `[-3; 3]`, 7³ samples. Chosen so the integer lattice contains `(±1, 0, 0)`
/// and `(±1, ±1, ±1)` — the directions in which each fixture's reverse distance
/// is extremal — which is what makes the assertions exact rather than
/// approximate.
fn unit_lattice() -> (RuntimeShape3, [f64; 3], AccuracyConfig) {
    (
        RuntimeShape3::new([7; 3]).expect("valid shape"),
        [-3.0; 3],
        AccuracyConfig::from_cell_size(1.0).expect("valid spacing"),
    )
}

/// Extract with Marching Cubes, using the same samples-vs-cells convention as
/// `mc/tests.rs`: `shape` counts samples, so `n` samples span `n − 1` cells.
fn mc_mesh<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    samples: u32,
) -> (MeshBuffer<f64>, f64, RuntimeShape3, [f64; 3]) {
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    mc.extract(field, &shape, lo, cell_size, &mut out)
        .expect("extraction");
    (out, cell_size, shape, lo)
}

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

// ─── closed-form fixtures ───────────────────────────────────────────────────

/// Octahedron inscribed in the unit sphere: the six vertices sit exactly on the
/// surface, so **the whole forward error comes from the centroids**. A harness
/// that sampled only vertices would report zero here.
#[test]
fn an_inscribed_octahedron_matches_the_closed_form() {
    let (positions, indices) = octahedron(1.0);
    let (shape, origin, cfg) = unit_lattice();
    let report = accuracy(
        &positions,
        &indices,
        &Sphere::<f64>::canonical(),
        &shape,
        origin,
        &cfg,
    )
    .expect("measurable");

    // Face centroid is at radius 1/√3, so its distance is 1 − 1/√3.
    let d = 1.0 - 1.0 / 3.0f64.sqrt();
    assert_eq!(report.vertex_samples, 6);
    assert_eq!(report.centroid_samples, 8);
    assert!(
        close(report.mesh_to_field.max, d, 1e-12),
        "max {:?} != {d:?}\n{report}",
        report.mesh_to_field.max
    );
    // Six zeros and eight copies of d, over fourteen samples.
    assert!(
        close(report.mean_absolute_error(), 8.0 * d / 14.0, 1e-12),
        "mean {:?}\n{report}",
        report.mean_absolute_error()
    );
    // Reverse extremum is the face-centre direction (1,1,1)/√3, the same value.
    assert!(
        close(report.field_to_mesh.max, d, 1e-9),
        "reverse {:?} != {d:?}\n{report}",
        report.field_to_mesh.max
    );
}

/// The same octahedron scaled so its faces are *tangent* to the unit sphere:
/// now the eight centroids sit exactly on the surface and **the whole forward
/// error comes from the vertices**. Together with the inscribed case this pins
/// both populations — a harness that drops either one passes exactly one of the
/// two tests.
///
/// It also separates the two directions: the forward max is `√3 − 1` while the
/// reverse max is `1 − 1/√3`, so `symmetric_hausdorff` demonstrably takes the
/// larger of two independently computed numbers rather than reporting one twice.
#[test]
fn a_circumscribed_octahedron_matches_the_closed_form() {
    let (positions, indices) = octahedron(3.0f64.sqrt());
    let (shape, origin, cfg) = unit_lattice();
    let report = accuracy(
        &positions,
        &indices,
        &Sphere::<f64>::canonical(),
        &shape,
        origin,
        &cfg,
    )
    .expect("measurable");

    let vertex = 3.0f64.sqrt() - 1.0;
    let reverse = 1.0 - 1.0 / 3.0f64.sqrt();
    assert!(
        close(report.mesh_to_field.max, vertex, 1e-12),
        "max {:?} != {vertex:?}\n{report}",
        report.mesh_to_field.max
    );
    assert!(
        close(report.mean_absolute_error(), 6.0 * vertex / 14.0, 1e-12),
        "mean {:?}\n{report}",
        report.mean_absolute_error()
    );
    assert!(
        close(report.field_to_mesh.max, reverse, 1e-9),
        "reverse {:?} != {reverse:?}\n{report}",
        report.field_to_mesh.max
    );
    // The forward direction is the larger, and symmetric picks it.
    assert!(
        close(report.symmetric_hausdorff(), vertex, 1e-12),
        "symmetric {:?}\n{report}",
        report.symmetric_hausdorff()
    );
}

/// Delete one face of the inscribed octahedron.
///
/// The forward direction cannot see this at all — every surviving vertex and
/// centroid is exactly where it was, so `mesh_to_field` is unchanged. Only the
/// reverse direction notices, and it is now the larger of the two. **This is the
/// test that justifies the reverse direction's cost**: misplaced geometry shows
/// up forwards, missing geometry only shows up backwards.
///
/// The exact value is the distance from the hole's centre direction
/// `(1,1,1)/√3` to the nearest surviving edge, whose closest point is the edge
/// midpoint `(½, ½, 0)`:
///
/// ```text
/// d² = 2·(½ − 1/√3)² + ⅓ = 3/2 − 2/√3
/// ```
#[test]
fn a_hole_is_found_only_by_the_reverse_direction() {
    let (positions, indices) = octahedron(1.0);
    let (shape, origin, cfg) = unit_lattice();
    let field = Sphere::<f64>::canonical();

    let whole = accuracy(&positions, &indices, &field, &shape, origin, &cfg).expect("measurable");
    // Drop the first face, the (+,+,+) octant.
    let holed_indices: Vec<u32> = indices[3..].to_vec();
    let holed =
        accuracy(&positions, &holed_indices, &field, &shape, origin, &cfg).expect("measurable");

    assert_eq!(holed.triangles, 7);
    assert!(
        close(holed.mesh_to_field.max, whole.mesh_to_field.max, 1e-15),
        "the forward direction should not notice a hole\n{holed}"
    );

    let expected = (1.5 - 2.0 / 3.0f64.sqrt()).sqrt();
    assert!(
        close(holed.field_to_mesh.max, expected, 1e-9),
        "reverse {:?} != {expected:?}\n{holed}",
        holed.field_to_mesh.max
    );
    assert!(
        holed.field_to_mesh.max > holed.mesh_to_field.max,
        "the hole should dominate\n{holed}"
    );
    assert!(close(holed.symmetric_hausdorff(), expected, 1e-9));
}

/// Two triangles lying exactly on a face of `BoxExact`'s `[-1,1]³` cube.
///
/// The field is exactly zero there and its gradient is a unit axis vector, so
/// the Newton step is exactly zero and the reported distance must be exactly
/// `0.0` — not "small". This is the fixture that catches a projector which
/// always moves the point, or a harness with a systematic bias.
#[test]
fn a_mesh_on_the_surface_reports_exactly_zero() {
    let positions = alloc::vec![
        [1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0],
        [1.0, 1.0, 1.0],
        [1.0, -1.0, 1.0],
    ];
    let indices = alloc::vec![0, 1, 2, 0, 2, 3];
    let (shape, origin, cfg) = unit_lattice();
    let report = accuracy(
        &positions,
        &indices,
        &BoxExact::<f64>::canonical(),
        &shape,
        origin,
        &cfg,
    )
    .expect("measurable");

    // Bit-compared rather than `== 0.0`, so that `-0.0` would fail too.
    assert_eq!(
        report.mesh_to_field.max.to_bits(),
        0.0f64.to_bits(),
        "a mesh on the surface must measure exactly zero\n{report}"
    );
    assert_eq!(
        report.mesh_to_field.mean.to_bits(),
        0.0f64.to_bits(),
        "{report}"
    );
    assert_eq!(report.unconverged_mesh_samples, 0, "{report}");
}

/// A max and a mean must respond differently to one bad vertex.
#[test]
fn a_perturbed_vertex_moves_the_max_but_barely_the_mean() {
    let field = Sphere::<f64>::canonical();
    let (mesh, cell_size, shape, origin) = mc_mesh(&field, 25);
    let cfg = AccuracyConfig::from_cell_size(cell_size).expect("valid spacing");

    let before =
        accuracy(&mesh.positions, &mesh.indices, &field, &shape, origin, &cfg).expect("measurable");

    // Push one vertex a long way out along its own radius.
    let mut moved = mesh.positions.clone();
    let p = moved[0];
    let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
    moved[0] = [
        p[0] * (1.0 + 0.5 / len),
        p[1] * (1.0 + 0.5 / len),
        p[2] * (1.0 + 0.5 / len),
    ];

    let after = accuracy(&moved, &mesh.indices, &field, &shape, origin, &cfg).expect("measurable");

    assert!(
        after.mesh_to_field.max > 3.0 * before.mesh_to_field.max,
        "the max should jump: {:?} -> {:?}",
        before.mesh_to_field.max,
        after.mesh_to_field.max
    );
    let mean_growth = after.mesh_to_field.mean / before.mesh_to_field.mean;
    assert!(
        mean_growth < 1.5,
        "one vertex in {} should barely move the mean, grew {mean_growth:?}x",
        before.mesh_to_field.samples
    );
}

// ─── the nearest-triangle primitive ─────────────────────────────────────────

/// One fixture per Voronoi region of the triangle.
///
/// Without this, a transposed subscript in the region test produces answers that
/// are correct for most points and wrong inside one wedge — which no end-to-end
/// test would localise.
#[test]
fn point_triangle_distance_hits_every_region() {
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    let cases: [([f64; 3], f64, &str); 7] = [
        ([0.25, 0.25, 5.0], 5.0, "face"),
        ([-1.0, -1.0, 0.0], 2.0f64.sqrt(), "vertex A"),
        ([2.0, 0.0, 0.0], 1.0, "vertex B"),
        ([0.0, 2.0, 0.0], 1.0, "vertex C"),
        ([0.5, -1.0, 0.0], 1.0, "edge AB"),
        ([-1.0, 0.5, 0.0], 1.0, "edge AC"),
        ([1.0, 1.0, 0.0], 0.5f64.sqrt(), "edge BC"),
    ];
    for (p, want, region) in cases {
        let got = point_triangle_distance_squared(p, a, b, c).sqrt();
        assert!(
            close(got, want, 1e-15),
            "{region}: got {got:?}, want {want:?}"
        );
    }
}

/// Exhaustive search over the same triangles, for the cross-check below.
fn brute_force(q: [f64; 3], positions: &[[f64; 3]], tris: &[[u32; 3]]) -> f64 {
    let mut best = f64::INFINITY;
    for t in tris {
        let d = point_triangle_distance_squared(
            q,
            positions[t[0] as usize],
            positions[t[1] as usize],
            positions[t[2] as usize],
        );
        if d < best {
            best = d;
        }
    }
    best
}

/// The grid is an accelerator, not a different algorithm.
///
/// Compared bit-for-bit rather than approximately: both routes compute the same
/// per-triangle values and take a minimum, and a minimum is a value rather than
/// an index, so exact equality is the correct expectation. An approximate
/// comparison here would hide the very shell-termination bug this exists to
/// catch.
///
/// Run at three spacings, including one far coarser than the mesh — that case
/// collapses the grid to a handful of cells and must still agree.
#[test]
fn the_grid_agrees_with_brute_force() {
    let field = Sphere::<f64>::canonical();
    let (mesh, _, _, _) = mc_mesh(&field, 17);
    let tris: Vec<[u32; 3]> = mesh
        .indices
        .chunks_exact(3)
        .map(|t| [t[0], t[1], t[2]])
        .collect();
    assert!(tris.len() > 100, "fixture should be non-trivial");

    for spacing in [0.05, 0.25, 4.0] {
        let grid = TriangleGrid::build(&mesh.positions, &tris, spacing).expect("grid");
        let mut checked = 0;
        for i in -4..=4 {
            for j in -4..=4 {
                for k in -4..=4 {
                    let q = [f64::from(i) * 0.6, f64::from(j) * 0.6, f64::from(k) * 0.6];
                    let fast = grid.nearest_distance_squared(q, &mesh.positions, &tris);
                    let slow = brute_force(q, &mesh.positions, &tris);
                    assert_eq!(
                        fast.to_bits(),
                        slow.to_bits(),
                        "spacing {spacing}: at {q:?} grid {fast:?} != brute {slow:?}"
                    );
                    checked += 1;
                }
            }
        }
        // Points exactly on the mesh, and far outside it.
        for q in [mesh.positions[0], [50.0, 50.0, 50.0], [0.0, 0.0, 0.0]] {
            let fast = grid.nearest_distance_squared(q, &mesh.positions, &tris);
            let slow = brute_force(q, &mesh.positions, &tris);
            assert_eq!(fast.to_bits(), slow.to_bits(), "spacing {spacing} at {q:?}");
            checked += 1;
        }
        assert_eq!(checked, 9 * 9 * 9 + 3);
    }
}

// ─── behaviour ──────────────────────────────────────────────────────────────

#[test]
fn results_are_deterministic() {
    let field = Sphere::<f64>::canonical();
    let (mesh, cell_size, shape, origin) = mc_mesh(&field, 21);
    let cfg = AccuracyConfig::from_cell_size(cell_size).expect("valid spacing");
    let run = || {
        accuracy(&mesh.positions, &mesh.indices, &field, &shape, origin, &cfg).expect("measurable")
    };
    let a = run();
    let b = run();
    assert_eq!(a, b);
    assert_eq!(a.mesh_to_field.max.to_bits(), b.mesh_to_field.max.to_bits());
    assert_eq!(a.field_to_mesh.max.to_bits(), b.field_to_mesh.max.to_bits());
}

#[test]
fn malformed_indices_are_counted_not_fatal() {
    let (positions, mut indices) = octahedron(1.0);
    indices.extend_from_slice(&[0, 0, 2]); // repeated index
    indices.extend_from_slice(&[0, 2, 99]); // out of range
    indices.extend_from_slice(&[0, 2]); // trailing partial
    let (shape, origin, cfg) = unit_lattice();
    let report = accuracy(
        &positions,
        &indices,
        &Sphere::<f64>::canonical(),
        &shape,
        origin,
        &cfg,
    )
    .expect("measurable");

    assert_eq!(report.triangles, 8);
    assert_eq!(report.faces_skipped, 3);
}

#[test]
fn a_meaningless_cell_size_is_rejected() {
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            AccuracyConfig::from_cell_size(bad).is_err(),
            "{bad:?} should be rejected"
        );
    }
}

#[test]
fn an_empty_mesh_reports_no_coverage() {
    let (shape, origin, cfg) = unit_lattice();
    let report =
        accuracy(&[], &[], &Sphere::<f64>::canonical(), &shape, origin, &cfg).expect("measurable");
    assert_eq!(report.triangles, 0);
    assert!(!report.has_coverage(), "{report}");
    assert_eq!(report.mesh_to_field.samples, 0);
    assert_eq!(report.field_to_mesh.samples, 0);
}

#[test]
fn a_grid_too_small_to_seed_is_rejected() {
    let cfg = AccuracyConfig::from_cell_size(1.0).expect("valid spacing");
    let shape = RuntimeShape3::new([1, 7, 7]).expect("valid shape");
    let (positions, indices) = octahedron(1.0);
    assert!(
        accuracy(
            &positions,
            &indices,
            &Sphere::<f64>::canonical(),
            &shape,
            [-3.0; 3],
            &cfg
        )
        .is_err()
    );
}

/// The sphere's gradient is `0/0` at its own centre, which the field module
/// documents. An odd, centred lattice lands on it exactly once.
///
/// Pinning this proves the non-finite case is *routed* rather than swallowed —
/// `Real::max` follows IEEE `maxNum`, which ignores NaN, so a dropped guard here
/// would silently vanish instead of poisoning the maximum.
#[test]
fn the_sphere_centre_is_a_reported_singularity() {
    let (positions, indices) = octahedron(1.0);
    let (shape, origin, cfg) = unit_lattice();
    let report = accuracy(
        &positions,
        &indices,
        &Sphere::<f64>::canonical(),
        &shape,
        origin,
        &cfg,
    )
    .expect("measurable");
    assert_eq!(report.seeds_non_finite, 1, "{report}");
}

/// Every lattice point lands in exactly one bucket. Structural, cheap, and it
/// makes a silently dropped seed impossible rather than merely unlikely.
#[test]
fn the_seed_counters_account_for_every_lattice_point() {
    crate::for_each_reference_field!(f64, |name, field| {
        let (mesh, cell_size, shape, origin) = mc_mesh(&field, 17);
        let cfg = AccuracyConfig::from_cell_size(cell_size).expect("valid spacing");
        let r = accuracy(&mesh.positions, &mesh.indices, &field, &shape, origin, &cfg)
            .expect("measurable");
        assert_eq!(
            r.seeds,
            r.seeds_out_of_band
                + r.seeds_non_finite
                + r.seeds_unconverged
                + r.seeds_outside_domain
                + r.field_to_mesh.samples,
            "{name}: a seed went missing\n{r}"
        );
        assert_eq!(r.seeds, 17 * 17 * 17, "{name}");
    });
}

/// Report, do not judge. Two of the seven fields have no error bound that can
/// be asserted a priori — `fbm_terrain` is open, so its reverse direction is
/// legitimately `O(h)` at the walls, and the capped gyroid's cap seam mixes
/// operands of different scales. These rows are the baseline M-001's shootout
/// table consumes.
#[test]
fn every_reference_field_is_measurable() {
    crate::for_each_reference_field!(f64, |name, field| {
        let (mesh, cell_size, shape, origin) = mc_mesh(&field, 33);
        let cfg = AccuracyConfig::from_cell_size(cell_size).expect("valid spacing");
        let r = accuracy(&mesh.positions, &mesh.indices, &field, &shape, origin, &cfg)
            .expect("measurable");
        std::println!(
            "measured: {name:16} h {:.6}  tris {:6}  fwd max {:.6} mean {:.6}  rev max {:.6}  seeds {}/{}",
            cell_size,
            r.triangles,
            r.mesh_to_field.max,
            r.mesh_to_field.mean,
            r.field_to_mesh.max,
            r.field_to_mesh.samples,
            r.seeds
        );
        assert!(r.has_coverage(), "{name} produced no coverage\n{r}");
    });
}

// ─── acceptance ─────────────────────────────────────────────────────────────

/// The ticket's criterion, for both extractors.
///
/// Also asserts a gate an order of magnitude tighter. The loose one is the
/// contract; the tight one is what actually notices when something breaks,
/// because the loose one passes with roughly 60× to spare.
#[test]
fn a_unit_sphere_at_64_cubed_is_within_one_cell_diagonal() {
    let field = Sphere::<f64>::canonical();
    let (lo, hi) = field.domain();
    let samples = 64u32;
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let cfg = AccuracyConfig::from_cell_size(cell_size).expect("valid spacing");

    let mut mc = MarchingCubes::<f64>::new();
    let mut sn = SurfaceNets::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();

    for which in ["mc", "sn"] {
        out.reset();
        if which == "mc" {
            mc.extract(&field, &shape, lo, cell_size, &mut out)
                .expect("extraction");
        } else {
            sn.extract(&field, &shape, lo, cell_size, &mut out)
                .expect("extraction");
        }
        let r =
            accuracy(&out.positions, &out.indices, &field, &shape, lo, &cfg).expect("measurable");
        std::println!(
            "measured: {which} 64^3 h {:.8} diagonal {:.8}  tris {}  symmetric {:.8}  mean {:.8}",
            cell_size,
            cfg.cell_diagonal(),
            r.triangles,
            r.symmetric_hausdorff(),
            r.mean_absolute_error()
        );
        r.panic_if_worse_than(cfg.cell_diagonal());
        assert!(
            r.symmetric_hausdorff() < 0.25 * cell_size,
            "{which}: within the contract but far worse than this extractor should be\n{r}"
        );
    }
}

/// The error must fall like `h²`.
///
/// This is the real answer to the acceptance criterion's slack: a harness that
/// returns a constant, or one with its units confused, passes "below one cell
/// diagonal" and fails here. The ideal ratio between 32 and 64 samples over the
/// same domain is `((4/31)/(4/63))² = 4.13`; 3.0 leaves room for grid phase.
///
/// Uses the mean rather than the max — a supremum over samples is noisy, a mean
/// is not.
#[test]
fn the_error_falls_like_h_squared() {
    let field = Sphere::<f64>::canonical();
    let mean_at = |samples: u32| {
        let (mesh, cell_size, shape, origin) = mc_mesh(&field, samples);
        let cfg = AccuracyConfig::from_cell_size(cell_size).expect("valid spacing");
        accuracy(&mesh.positions, &mesh.indices, &field, &shape, origin, &cfg)
            .expect("measurable")
            .mean_absolute_error()
    };
    let coarse = mean_at(32);
    let fine = mean_at(64);
    let ratio = coarse / fine;
    std::println!("measured: mean error 32^3 {coarse:.9}, 64^3 {fine:.9}, ratio {ratio:.3}");
    assert!(
        ratio >= 3.0,
        "error should fall like h^2; got a ratio of {ratio:?}"
    );
}
