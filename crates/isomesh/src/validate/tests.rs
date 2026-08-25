//! Fixtures for the validity harness.
//!
//! One clean closed surface with genus 0, one with genus 1, and one broken mesh
//! per violation class. The genus-1 fixture is not decoration: it is what proves
//! `χ` is computed rather than hard-coded, and that a hole is distinguishable
//! from a handle.

use super::*;
use alloc::vec;
use alloc::vec::Vec;
use alloc::{format, string::String};

fn cfg() -> ValidateConfig {
    ValidateConfig::from_cell_size(1.0).expect("valid cell size")
}

// ─── clean fixtures ─────────────────────────────────────────────────────────

/// A regular tetrahedron centred at the origin, wound counter-clockwise as seen
/// from outside.
///
/// The winding is verified rather than asserted: for each face, the cross
/// product of its two edge vectors dotted with the face centroid is `+4`, and
/// since the solid is centred at the origin a positive dot means the normal
/// points away from it. `outward_winding_is_verified` re-checks that in code so
/// the claim cannot rot.
fn tetrahedron() -> (Vec<[f64; 3]>, Vec<u32>) {
    let positions = vec![
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ];
    let indices = vec![0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2];
    (positions, indices)
}

/// An `m × n` grid wrapped in both directions and split into triangles: a torus
/// with `V = mn`, `E = 3mn`, `F = 2mn`, and therefore `χ = 0` for any
/// `m, n >= 3`.
fn torus_grid(m: u32, n: u32) -> (Vec<[f64; 3]>, Vec<u32>) {
    let tau = 2.0 * core::f64::consts::PI;
    let (major, minor) = (1.0f64, 0.3f64);
    let mut positions = Vec::new();
    for i in 0..m {
        let theta = tau * f64::from(i) / f64::from(m);
        for j in 0..n {
            let phi = tau * f64::from(j) / f64::from(n);
            let ring = major + minor * Real::cos(phi);
            positions.push([
                ring * Real::cos(theta),
                minor * Real::sin(phi),
                ring * Real::sin(theta),
            ]);
        }
    }

    let at = |i: u32, j: u32| (i % m) * n + (j % n);
    let mut indices = Vec::new();
    for i in 0..m {
        for j in 0..n {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            indices.extend_from_slice(&[a, b, c, a, c, d]);
        }
    }
    (positions, indices)
}

#[test]
fn outward_winding_is_verified() {
    let (p, idx) = tetrahedron();
    for tri in idx.as_chunks::<3>().0 {
        let (a, b, c) = (p[tri[0] as usize], p[tri[1] as usize], p[tri[2] as usize]);
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        // The solid is centred at the origin, so a face centroid is an outward
        // direction and the normal must agree with it.
        let centroid = [
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        ];
        let dot = n[0] * centroid[0] + n[1] * centroid[1] + n[2] * centroid[2];
        assert!(dot > 0.0, "face {tri:?} is wound inward: {dot}");
    }
}

#[test]
fn tetrahedron_is_a_clean_closed_surface() {
    let (p, idx) = tetrahedron();
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.vertices, 4);
    assert_eq!(r.referenced_vertices, 4);
    assert_eq!(r.edges, 6);
    assert_eq!(r.faces, 4);
    assert_eq!(r.euler_characteristic, 2);
    assert_eq!(r.components, 1);
    assert_eq!(r.boundary_loops, 0);
    assert_eq!(r.genus, Some(0));
    assert_eq!(r.violations(), 0);
    assert_eq!(r.degenerate_triangles, 0);
    assert_eq!(r.duplicate_vertices, 0);
    assert!(r.is_manifold() && r.is_closed() && !r.has_structural_errors());
}

/// The fixture that proves `χ` is measured, not assumed — and that `genus`
/// distinguishes a handle from a hole.
#[test]
fn torus_grid_has_genus_one() {
    let (p, idx) = torus_grid(8, 5);
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.vertices, 40);
    assert_eq!(r.edges, 120);
    assert_eq!(r.faces, 80);
    assert_eq!(r.euler_characteristic, 0);
    assert_eq!(r.components, 1);
    assert_eq!(r.boundary_loops, 0);
    assert_eq!(r.genus, Some(1)); // a handle, not a hole
    assert_eq!(r.violations(), 0);
    assert!(r.is_closed());
}

// ─── one fixture per violation ──────────────────────────────────────────────

/// A tetrahedron with one face removed: a disk. `χ = 1` *with one boundary
/// loop*, which is what makes it unambiguous — `χ = 1` alone could equally be
/// something with a handle.
#[test]
fn open_tetra_is_a_disk_not_a_handle() {
    let (p, mut idx) = tetrahedron();
    idx.truncate(9);
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.faces, 3);
    assert_eq!(r.edges, 6);
    assert_eq!(r.euler_characteristic, 1);
    assert_eq!(r.boundary_edges, 3);
    assert_eq!(r.boundary_loops, 1);
    assert_eq!(r.genus, Some(0));
    assert_eq!(r.non_manifold_edges, 0);
    assert!(
        r.is_manifold(),
        "a manifold with boundary is still a manifold"
    );
    assert!(!r.is_closed());
}

#[test]
fn non_manifold_edge_is_detected() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, -1.0, 0.0],
    ];
    // Three faces sharing the edge {0, 1}.
    let idx = vec![0, 1, 2, 0, 1, 3, 0, 1, 4];
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.non_manifold_edges, 1);
    assert_eq!(r.boundary_edges, 6);
    // Both endpoints of the branching edge have an umbrella rather than a fan.
    assert_eq!(r.non_manifold_vertices, 2);
    assert!(!r.is_manifold());
}

/// The fixture that pays for the connected-components walk.
///
/// Two triangles meeting at a single vertex. Every edge has exactly one face, no
/// edge is non-manifold, and the cheap "incident faces equals incident edges"
/// test reports the apex as clean — it has two faces and two wing edges. Only a
/// walk over the link finds it.
#[test]
fn bowtie_is_detected_only_by_the_link_walk() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [-1.0, -1.0, 0.0],
    ];
    let idx = vec![0, 1, 2, 0, 3, 4];
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.non_manifold_vertices, 1);
    assert_eq!(
        r.non_manifold_edges, 0,
        "the cheap edge test cannot see this"
    );
    assert_eq!(r.components, 2);
    assert_eq!(r.euler_characteristic, 1);
    assert!(!r.is_manifold());

    // The apex has as many incident faces as wing edges, which is exactly why
    // counting them would report this mesh clean.
    let (faces_at_apex, wing_edges_at_apex) = (2, 2);
    assert_eq!(faces_at_apex, wing_edges_at_apex);
}

#[test]
fn degenerate_triangle_is_detected() {
    let p = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let idx = vec![0, 1, 2];
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.degenerate_triangles, 1);
    assert_eq!(r.faces, 1);
    // A sliver is not a structural error: marching cubes emits them for
    // perfectly ordinary reasons.
    assert!(!r.has_structural_errors());
    assert_eq!(r.violations(), 0);
}

#[test]
fn repeated_index_triangle_is_skipped_and_counted() {
    let (p, mut idx) = tetrahedron();
    idx.extend_from_slice(&[0, 1, 1]);
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.repeated_index_triangles, 1);
    assert_eq!(r.faces_skipped, 1);
    assert_eq!(r.faces, 4, "the four good faces still count");
    assert_eq!(r.euler_characteristic, 2);
    assert!(r.has_structural_errors());
    assert_eq!(r.genus, None, "genus is undefined once a face was skipped");
}

#[test]
fn out_of_range_index_does_not_panic() {
    let (p, mut idx) = tetrahedron();
    let last = idx.len() - 1;
    idx[last] = 9; // only four vertices exist
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.out_of_range_indices, 1);
    assert_eq!(r.faces_skipped, 1);
    assert_eq!(r.faces, 3);
    assert!(r.has_structural_errors());
    // Every other pass still ran over the valid subset.
    assert_eq!(r.edges, 6);
    assert_eq!(r.boundary_edges, 3);
}

#[test]
fn trailing_index_is_counted_and_ignored() {
    let (p, idx) = tetrahedron();
    let mut broken = idx.clone();
    broken.push(0);
    let clean = validate_indexed(&p, &idx, &cfg());
    let r = validate_indexed(&p, &broken, &cfg());

    assert_eq!(r.trailing_indices, 1);
    assert!(r.has_structural_errors());
    // The partial triangle is ignored; everything else is untouched.
    assert_eq!(r.faces, clean.faces);
    assert_eq!(r.edges, clean.edges);
    assert_eq!(r.euler_characteristic, clean.euler_characteristic);
}

#[test]
fn exact_duplicate_vertex_is_detected() {
    let (mut p, idx) = tetrahedron();
    p.push(p[0]); // never referenced
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.duplicate_vertices, 1);
    assert_eq!(r.unreferenced_vertices, 1);
    assert_eq!(r.vertices, 5);
    assert_eq!(r.referenced_vertices, 4);
    assert_eq!(
        r.euler_characteristic, 2,
        "the stale vertex must not inflate chi"
    );
}

/// An epsilon that is never crossed is decoration, so this crosses it in both
/// directions.
#[test]
fn near_duplicate_is_detected_at_one_epsilon_and_not_at_a_quarter_of_it() {
    let (mut p, idx) = tetrahedron();
    let coarse = cfg();
    let half_eps = coarse.weld_epsilon * 0.5;
    p.push([p[0][0] + half_eps, p[0][1], p[0][2]]);

    assert_eq!(validate_indexed(&p, &idx, &coarse).duplicate_vertices, 1);

    // A quarter of the cell size gives a quarter of the weld epsilon, which the
    // same pair no longer falls inside.
    let fine = ValidateConfig::from_cell_size(coarse.cell_size * 0.25).expect("valid cell size");
    assert!(half_eps > fine.weld_epsilon);
    assert_eq!(validate_indexed(&p, &idx, &fine).duplicate_vertices, 0);
}

/// The reason `inconsistently_oriented_edges` exists.
///
/// One flipped face passes the Euler check, edge manifoldness *and* vertex
/// manifoldness. Without this counter the mesh reports clean and is inside out.
#[test]
fn flipped_winding_passes_every_other_check() {
    let (p, mut idx) = tetrahedron();
    idx.swap(1, 2); // face 0 becomes (0, 2, 1)
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.euler_characteristic, 2, "chi alone cannot see this");
    assert_eq!(r.non_manifold_edges, 0);
    assert_eq!(r.non_manifold_vertices, 0);
    assert_eq!(r.boundary_edges, 0);

    assert_eq!(r.inconsistently_oriented_edges, 3);
    assert!(!r.is_manifold());
    assert_eq!(r.genus, None);
}

/// Two nested shells: one solid with an enclosed void. `χ = 4` because there are
/// two boundary components, and `genus` is `None` because the formula applies to
/// one component at a time.
#[test]
fn enclosed_cavity_has_two_components() {
    let (outer, outer_idx) = tetrahedron();
    let mut p = outer.clone();
    let mut idx = outer_idx.clone();
    for v in &outer {
        p.push([v[0] * 0.3, v[1] * 0.3, v[2] * 0.3]);
    }
    // The inner shell is wound inward, as the boundary of a void must be.
    for tri in outer_idx.as_chunks::<3>().0 {
        idx.extend_from_slice(&[tri[0] + 4, tri[2] + 4, tri[1] + 4]);
    }
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.components, 2);
    assert_eq!(r.faces, 8);
    assert_eq!(r.edges, 12);
    assert_eq!(r.euler_characteristic, 4);
    assert_eq!(r.boundary_loops, 0);
    assert_eq!(r.inconsistently_oriented_edges, 0);
    assert_eq!(r.genus, None, "genus needs a single component");
    assert!(r.is_closed());
}

/// The acceptance criterion: every violation class non-zero at the same time,
/// with nothing panicking and every pass still completing.
#[test]
fn composite_broken_mesh_detects_every_violation_class() {
    let mut p = vec![
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ];
    p.push(p[0]); // exact duplicate, unreferenced
    p.push([0.0, 0.0, 0.0]);
    p.push([1.0, 0.0, 0.0]);
    p.push([2.0, 0.0, 0.0]); // collinear with the previous two
    p.push([f64::NAN, 0.0, 0.0]); // non-finite

    let mut idx = vec![
        0, 2, 1, // flipped winding
        0, 2, 3, //
        0, 3, 1, //
        1, 3, 2, //
        0, 1, 9, // out of range
        0, 1, 1, // repeated index
        5, 6, 7, // zero area
    ];
    idx.push(4); // trailing partial triangle

    let r = validate_indexed(&p, &idx, &cfg());

    assert!(r.inconsistently_oriented_edges > 0);
    assert!(r.degenerate_triangles > 0);
    assert!(r.repeated_index_triangles > 0);
    assert!(r.out_of_range_indices > 0);
    assert!(r.trailing_indices > 0);
    assert!(r.duplicate_vertices > 0);
    assert!(r.unreferenced_vertices > 0);
    assert!(r.non_finite_positions > 0);
    assert!(r.faces_skipped > 0);
    assert!(r.has_structural_errors());
    assert!(!r.is_manifold());
    assert_eq!(r.genus, None);

    let block = format!("{r}");
    assert!(block.contains("!! STRUCTURAL ERRORS"));
    assert!(block.contains("(over the valid subset only)"));
    assert!(block.contains("INVALID:"));
}

// ─── behaviour of the harness itself ────────────────────────────────────────

/// Guards the concern the project rules single out: no iteration order may leak
/// into the output. The sort keys end in a face or vertex index, so no two
/// entries ever compare equal and the result is a pure function of the values.
#[test]
fn validation_is_deterministic() {
    let (p, idx) = torus_grid(6, 4);
    let a = validate_indexed(&p, &idx, &cfg());
    let b = validate_indexed(&p, &idx, &cfg());
    assert_eq!(a, b);
    assert_eq!(format!("{a}"), format!("{b}"));
}

/// The case every validator gets wrong once.
#[test]
fn empty_mesh_does_not_divide_by_zero() {
    let p: Vec<[f64; 3]> = Vec::new();
    let idx: Vec<u32> = Vec::new();
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.vertices, 0);
    assert_eq!(r.faces, 0);
    assert_eq!(r.edges, 0);
    assert_eq!(r.euler_characteristic, 0);
    assert_eq!(r.components, 0);
    assert_eq!(r.violations(), 0);
    let _ = format!("{r}");
}

/// `validate` is a forwarder plus the one check that needs the normals.
#[test]
fn mesh_buffer_normal_count_is_checked() {
    let (p, idx) = tetrahedron();
    let mut mesh = MeshBuffer::<f64> {
        positions: p,
        normals: Vec::new(),
        indices: idx,
    };
    let r = validate(&mesh, &cfg());
    assert!(r.normal_count_mismatch);
    assert!(r.has_structural_errors());

    mesh.normals = mesh.positions.clone();
    let r = validate(&mesh, &cfg());
    assert!(!r.normal_count_mismatch);
    assert!(r.is_closed());
}

#[test]
fn config_rejects_a_meaningless_scale() {
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let error = ValidateConfig::from_cell_size(bad)
            .expect_err("a meaningless spacing makes every threshold meaningless");
        assert!(
            matches!(error, crate::Error::InvalidCellSize { .. }),
            "{error}"
        );
    }
    // And the invalid state is unrepresentable: the fields are private, so the
    // checked constructor is the only way to obtain one.
    let good = ValidateConfig::from_cell_size(0.5).expect("valid");
    assert!((good.cell_size() - 0.5).abs() < f64::EPSILON);
    assert!(good.weld_epsilon() > 0.0 && good.area_epsilon_rel() > 0.0);
}

#[test]
#[should_panic(expected = "MANIFOLD, WITH BOUNDARY")]
fn panic_if_invalid_reports_the_whole_block() {
    let (p, mut idx) = tetrahedron();
    idx.truncate(9);
    validate_indexed(&p, &idx, &cfg()).panic_if_invalid(true);
}

/// The rule the whole harness exists to support: **one gate, selected by the
/// field**, with no per-field branch in test code.
///
/// There is no extractor yet, so each field is paired with a stand-in mesh of
/// the shape a correct extraction would produce — closed for a closed field,
/// bounded for an open one. What this proves is that the selection compiles and
/// discriminates; A-001 is where a real extracted mesh first flows through this
/// exact path, and this test is the shape it will take.
#[test]
fn the_gate_is_chosen_by_the_field_not_by_the_test() {
    use crate::fields::ReferenceField;
    use crate::for_each_reference_field;

    let (closed_p, closed_i) = tetrahedron();
    let (open_p, mut open_i) = tetrahedron();
    open_i.truncate(9);

    let mut checked = 0;
    for_each_reference_field!(f64, |name, field| {
        // Stand in for what a correct extraction of this field would produce.
        let (p, i) = if field.closed_in_domain() {
            (&closed_p, &closed_i)
        } else {
            (&open_p, &open_i)
        };
        let report = validate_indexed(p, i, &cfg());

        if field.closed_in_domain() {
            assert!(report.is_closed(), "{name}\n{report}");
        } else {
            assert!(report.is_manifold(), "{name}\n{report}");
            assert!(
                report.boundary_edges > 0,
                "{name}: expected an open surface"
            );
        }
        if let Some(chi) = field.expected_euler() {
            // The stand-in is a sphere, so only the genus-0 fields can match it;
            // what matters is that the expectation comes from the field.
            assert!(
                chi == 2 || chi == 0,
                "{name}: unexpected declared chi {chi}"
            );
        }
        checked += 1;
    });
    assert_eq!(checked, 8);
}

/// Freezes the rendered block. It is the thing a golden-hash regression will
/// hash, so its shape is deliberately pinned here rather than left to drift.
#[test]
fn display_block_is_stable() {
    let (p, idx) = tetrahedron();
    let block: String = format!("{}", validate_indexed(&p, &idx, &cfg()));
    assert_eq!(block, EXPECTED_TETRAHEDRON_BLOCK);
}

const EXPECTED_TETRAHEDRON_BLOCK: &str = "\
mesh report
  vertices                        4  (4 referenced, 0 unreferenced)
  edges                           6
  faces                           4
  euler characteristic            2
  components                      1
  boundary loops                  0
  genus                           0
  ------------------------------------------------------------
  non-manifold edges              0
  non-manifold vertices           0
  boundary edges                  0
  inconsistently oriented         0
  degenerate triangles            0  (area <= 1e-6 * h^2, h = 1)
  repeated-index triangles        0
  duplicate vertices              0  (within 1e-4)
  out-of-range indices            0
  trailing indices                0
  non-finite positions            0
  ------------------------------------------------------------
  MANIFOLD, CLOSED";

/// T-008: the weld lattice must not care where in the world the mesh is.
///
/// `quantise` narrows through `as_f32`, which stops distinguishing consecutive
/// integers above `2²⁴`. Scaled by `1/weld_epsilon` — 160,000 at `h = 0.0625` —
/// an *absolute* coordinate crosses that at about 105 world units (M-18). Past
/// it, neighbouring cells merge and the 27-cell probe degrades toward comparing
/// everything with everything.
///
/// **The count alone cannot see this**, which is why the assertion is on the
/// occupancy. Coarsening only ever merges buckets, the probe still covers them
/// and the exact distance test still runs, so `duplicate_vertices` stays right
/// while the scan slides toward quadratic. A test that checked only the answer
/// would have passed against the broken version.
#[test]
fn the_weld_lattice_does_not_collapse_far_from_the_origin() {
    use crate::fields::ReferenceField as _;
    let field = crate::fields::Sphere::<f64>::canonical();
    let (lo, _) = field.domain();
    let cell_size = 4.0 / 32.0;
    let shape = crate::RuntimeShape3::new([33; 3]).expect("valid shape");
    let mut mc = crate::marching_cubes::MarchingCubes::<f64>::new();
    let mut mesh = MeshBuffer::<f64>::new();
    mc.extract(&field, &shape, lo, cell_size, &mut mesh)
        .expect("extraction");
    let cfg = ValidateConfig::from_cell_size(cell_size).expect("valid spacing");

    let at = |offset: f64| {
        let moved: Vec<[f64; 3]> = mesh
            .positions
            .iter()
            .map(|p| [p[0] + offset, p[1] + offset, p[2] + offset])
            .collect();
        validate_indexed(&moved, &mesh.indices, &cfg)
    };

    let home = at(0.0);
    let far = at(1.0e4);
    let further = at(1.0e6);

    std::println!(
        "measured: T-008 weld buckets for {} vertices -- origin {}, +1e4 {}, +1e6 {}",
        mesh.vertex_count(),
        home.weld_buckets,
        far.weld_buckets,
        further.weld_buckets
    );

    // The answer is unchanged, which is the part that was already true.
    assert_eq!(home.duplicate_vertices, far.duplicate_vertices);
    assert_eq!(home.duplicate_vertices, further.duplicate_vertices);

    // The occupancy is unchanged, which is the part T-008 fixed. Every vertex
    // of an extracted sphere is its own bucket, so a collapse is unmissable.
    assert_eq!(
        home.weld_buckets, far.weld_buckets,
        "the lattice collapsed at 1e4"
    );
    assert_eq!(
        home.weld_buckets, further.weld_buckets,
        "the lattice collapsed at 1e6"
    );
    assert!(
        home.weld_buckets as usize > mesh.vertex_count() * 9 / 10,
        "a sphere's vertices should be spread across buckets, got {} for {}",
        home.weld_buckets,
        mesh.vertex_count()
    );
}

// ─── T-020: fan folds, and which counter can see one ────────────────────────

/// Fan-triangulate a closed planar loop around its centroid, giving every
/// triangle the `+z` normal by flipping its winding when the slice faces the
/// other way.
///
/// This reproduces `bevy_autogib`'s `push_cap_tri` rule rather than depending
/// on it: that crate is outside this workspace, and a fixture that moves when
/// somebody else's cap builder moves is not a fixture. The loop lies in `z = 0`
/// and the outward normal is `+z`, so a slice is flipped exactly when its
/// signed area is negative.
///
/// Vertex 0 is the apex. Vertex `1 + i` is `loop_xy[i]`.
fn cap_fan(loop_xy: &[[f64; 2]]) -> (Vec<[f64; 3]>, Vec<u32>) {
    let n = loop_xy.len();
    let cx = loop_xy.iter().map(|p| p[0]).sum::<f64>() / n as f64;
    let cy = loop_xy.iter().map(|p| p[1]).sum::<f64>() / n as f64;

    let mut positions = vec![[cx, cy, 0.0]];
    positions.extend(loop_xy.iter().map(|p| [p[0], p[1], 0.0]));

    let mut indices = Vec::new();
    for i in 0..n {
        let (a, b) = (1 + i as u32, 1 + ((i + 1) % n) as u32);
        let (u, v) = (loop_xy[i], loop_xy[(i + 1) % n]);
        // z-component of (u − c) × (v − c): twice the slice's signed area.
        let twice_area = (u[0] - cx) * (v[1] - cy) - (u[1] - cy) * (v[0] - cx);
        if twice_area >= 0.0 {
            indices.extend_from_slice(&[0, a, b]);
        } else {
            indices.extend_from_slice(&[0, b, a]);
        }
    }
    (positions, indices)
}

/// How many times the loop winds around its own centroid.
///
/// The independent witness that a fan is broken, computed from the loop alone
/// and never from the counters under test. A cap whose boundary winds twice
/// covers its disk twice, whatever the validator makes of it.
fn winding_turns(loop_xy: &[[f64; 2]]) -> f64 {
    let n = loop_xy.len();
    let cx = loop_xy.iter().map(|p| p[0]).sum::<f64>() / n as f64;
    let cy = loop_xy.iter().map(|p| p[1]).sum::<f64>() / n as f64;
    let tau = 2.0 * core::f64::consts::PI;

    let mut total = 0.0;
    for i in 0..n {
        let (u, v) = (loop_xy[i], loop_xy[(i + 1) % n]);
        let a = libm::atan2(u[1] - cy, u[0] - cx);
        let b = libm::atan2(v[1] - cy, v[0] - cx);
        let mut d = b - a;
        while d > core::f64::consts::PI {
            d -= tau;
        }
        while d < -core::f64::consts::PI {
            d += tau;
        }
        total += d;
    }
    total / tau
}

/// The control: a convex loop never reverses its sweep, so nothing fires.
#[test]
fn a_convex_cap_fan_is_consistently_oriented() {
    let hexagon: Vec<[f64; 2]> = (0..6)
        .map(|k| {
            let t = 2.0 * core::f64::consts::PI * f64::from(k) / 6.0;
            [Real::cos(t), Real::sin(t)]
        })
        .collect();

    assert!((winding_turns(&hexagon) - 1.0).abs() < 1e-12, "winds once");

    let (p, idx) = cap_fan(&hexagon);
    let r = validate_indexed(&p, &idx, &cfg());

    assert_eq!(r.faces, 6);
    assert_eq!(r.inconsistently_oriented_edges, 0);
    assert_eq!(r.non_manifold_edges, 0);
    assert_eq!(r.degenerate_triangles, 0);
    assert_eq!(r.boundary_edges, 6, "the rim, and only the rim");
}

/// **T-020's acceptance.** A fan that reverses its sweep is invisible to
/// `self_intersections` and unmissable in `inconsistently_oriented_edges`.
///
/// The mechanism: adjacent slices share the spoke `{apex, p_i}` and traverse it
/// in *opposite* directions whenever they agree about which way is outward. The
/// per-triangle flip reverses one of them, so a slice that disagrees with its
/// neighbour makes them traverse that spoke the **same** way — which is what
/// this counter is. Equivalently, it counts flip-state *changes* around the fan.
///
/// Both numbers are asserted in the same test on purpose: the zero is vacuous
/// (M-83 — every pair shares the apex, so none is ever compared) and reading it
/// alone would say this mesh is clean.
#[test]
fn a_fan_that_reverses_its_sweep_is_caught_by_orientation_and_missed_by_self_intersection() {
    // Four points sweeping counter-clockwise, then one that backtracks.
    let folded = [
        [2.0, 0.0],
        [0.0, 2.0],
        [-2.0, 0.0],
        [0.0, -2.0],
        [-1.0, -1.0],
    ];

    assert!(
        (winding_turns(&folded) - 1.0).abs() < 1e-12,
        "the loop still winds exactly once -- the defect is the reversal, not a double cover"
    );

    let (p, idx) = cap_fan(&folded);
    let r = validate_indexed(&p, &idx, &cfg());
    let si = self_intersections(&p, &idx, 1.0).expect("valid cell size");

    std::println!(
        "measured: T-020 reversing fold -- inconsistently_oriented_edges {}, \
         self-intersections {} from {} tested, {} skipped for sharing the apex",
        r.inconsistently_oriented_edges,
        si.count(),
        si.tested_pairs,
        si.adjacent_pairs_skipped
    );

    // Exactly one slice disagrees with its neighbours, so exactly two spokes
    // change flip state: the one entering the reversal and the one leaving it.
    assert_eq!(r.inconsistently_oriented_edges, 2);

    // And the counter that is supposed to find folds finds nothing, because it
    // never compared anything. Asserting the skip count rather than only the
    // zero is what keeps the green tick from reading as coverage.
    assert_eq!(si.count(), 0);
    assert_eq!(si.tested_pairs, 0);
    assert_eq!(
        si.adjacent_pairs_skipped, 10,
        "5 choose 2, all sharing apex"
    );

    assert_eq!(r.non_manifold_edges, 0, "a fold is not a non-manifold edge");
    assert_eq!(r.degenerate_triangles, 0);
}

/// **The limit of the equivalence, and why it must be stated with its
/// construction attached.**
///
/// A star polygon folds over itself without ever reversing its sweep: every
/// slice of a pentagram keeps the same sign, so no triangle is flipped, so no
/// spoke changes flip state. The cap double-covers its disk and *both* counters
/// report zero. QEx (Ebke, Kobbelt, Bommes & Campen, `10.1145/2508363.2508372`)
/// records the analogous special case for their fold-over test — "this
/// assumption fails in one special case … the entire triangle fan spans an
/// absolute range of less than 180°" — which is what prompted looking for this
/// one rather than asserting the equivalence was total.
///
/// So `inconsistently_oriented_edges > 0` is an exact test for a fan that
/// *reverses*, not for a fan that *folds*.
#[test]
fn a_star_polygon_fan_double_covers_its_disk_and_both_counters_report_zero() {
    // Every second vertex of a regular pentagon: each slice spans +144°.
    let pentagram: Vec<[f64; 2]> = (0..5)
        .map(|k| {
            let t = 2.0 * core::f64::consts::PI * f64::from(2 * k) / 5.0;
            [Real::cos(t), Real::sin(t)]
        })
        .collect();

    let turns = winding_turns(&pentagram);
    assert!(
        (turns - 2.0).abs() < 1e-12,
        "the boundary winds twice, so the cap covers its disk twice: {turns}"
    );

    let (p, idx) = cap_fan(&pentagram);
    let r = validate_indexed(&p, &idx, &cfg());
    let si = self_intersections(&p, &idx, 1.0).expect("valid cell size");

    std::println!(
        "measured: T-020 star polygon -- winding {turns}, \
         inconsistently_oriented_edges {}, self-intersections {} from {} tested",
        r.inconsistently_oriented_edges,
        si.count(),
        si.tested_pairs
    );

    assert_eq!(
        r.inconsistently_oriented_edges, 0,
        "no slice is flipped, so no spoke changes flip state"
    );
    assert_eq!(si.count(), 0);
    assert_eq!(si.tested_pairs, 0);
}

// ── T-023: the gate rule, now that it ships ─────────────────────────────────

/// A report with every counter clean, to be spoiled one field at a time.
fn clean_report() -> MeshReport {
    MeshReport {
        vertices: 4,
        referenced_vertices: 4,
        edges: 6,
        faces: 4,
        faces_skipped: 0,
        euler_characteristic: 2,
        components: 1,
        boundary_loops: 0,
        genus: Some(0),
        non_manifold_edges: 0,
        non_manifold_vertices: 0,
        boundary_edges: 0,
        inconsistently_oriented_edges: 0,
        degenerate_triangles: 0,
        mean_ratio: 1.0,
        irregular_vertices: 0,
        repeated_index_triangles: 0,
        duplicate_vertices: 0,
        unreferenced_vertices: 0,
        weld_buckets: 4,
        out_of_range_indices: 0,
        trailing_indices: 0,
        non_finite_positions: 0,
        normal_count_mismatch: false,
        config: ValidateConfig::from_cell_size(1.0).expect("valid cell size"),
    }
}

#[test]
fn a_closed_solid_satisfies_every_gate() {
    let report = clean_report();
    assert!(report.satisfies(SurfaceGate::Closed));
    assert!(report.satisfies(SurfaceGate::Manifold));
    assert!(report.satisfies(SurfaceGate::ClosedAllowingUnresolvedTopology));
}

/// The case the whole ticket exists for: an open surface is **not** a broken
/// solid, and only the gate can tell those apart.
#[test]
fn an_open_surface_fails_the_solid_gate_and_passes_the_surface_one() {
    let mut report = clean_report();
    report.boundary_edges = 12;
    report.boundary_loops = 1;
    report.euler_characteristic = 1;

    assert!(
        !report.satisfies(SurfaceGate::Closed),
        "a surface with boundary is not a closed solid"
    );
    assert!(
        report.satisfies(SurfaceGate::Manifold),
        "and that is not a defect -- an open field is supposed to have boundary"
    );
    assert!(
        !report.satisfies(SurfaceGate::ClosedAllowingUnresolvedTopology),
        "this gate still requires the surface not to leave the grid"
    );
}

/// The loosest gate waives non-manifoldness and nothing else.
#[test]
fn the_unresolved_topology_gate_waives_manifoldness_only() {
    let mut report = clean_report();
    report.non_manifold_edges = 3;
    report.euler_characteristic = 1;

    assert!(!report.satisfies(SurfaceGate::Closed));
    assert!(!report.satisfies(SurfaceGate::Manifold));
    assert!(
        report.satisfies(SurfaceGate::ClosedAllowingUnresolvedTopology),
        "a grid too coarse to resolve its field may be non-manifold"
    );

    // But winding and boundary are still asserted under it.
    let mut flipped = report;
    flipped.inconsistently_oriented_edges = 1;
    assert!(!flipped.satisfies(SurfaceGate::ClosedAllowingUnresolvedTopology));

    let mut leaked = report;
    leaked.boundary_edges = 1;
    assert!(!leaked.satisfies(SurfaceGate::ClosedAllowingUnresolvedTopology));
}

/// A malformed mesh satisfies nothing, including the loosest gate. `is_closed`
/// and `is_manifold` already fold this in; the third arm would not without the
/// explicit check, which is why `satisfies` hoists it.
#[test]
fn a_structurally_broken_mesh_satisfies_no_gate() {
    let mut report = clean_report();
    report.out_of_range_indices = 1;
    assert!(report.has_structural_errors());
    for gate in [
        SurfaceGate::Closed,
        SurfaceGate::Manifold,
        SurfaceGate::ClosedAllowingUnresolvedTopology,
    ] {
        assert!(
            !report.satisfies(gate),
            "{gate:?} must reject malformed input"
        );
    }
}

// ── T-026: mean ratio and irregular vertices ────────────────────────────────

/// **The metric is pinned on shapes whose answer is known in closed form**, not
/// on a mesh whose value nobody can check.
///
/// `q = 4√3·A / Σlᵢ²` is **1** for an equilateral triangle by construction, and
/// for a right isosceles triangle with legs `1` it is
/// `4√3 · ½ / (1 + 1 + 2) = √3/2 ≈ 0.8660`. Both are computed here rather than
/// quoted, so the test states the formula twice in different forms and would
/// catch a transcription error in either.
#[test]
fn mean_ratio_is_one_for_equilateral_and_known_for_a_right_isosceles() {
    let cfg = ValidateConfig::from_cell_size(1.0).expect("valid cell size");

    let h = 3.0_f64.sqrt() / 2.0;
    let equilateral = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, h, 0.0]];
    let r = validate_indexed(&equilateral, &[0, 1, 2], &cfg);
    assert!(
        (r.mean_ratio - 1.0).abs() < 1e-12,
        "equilateral should be exactly 1, got {}",
        r.mean_ratio
    );

    let right = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let r = validate_indexed(&right, &[0, 1, 2], &cfg);
    let expected = 3.0_f64.sqrt() / 2.0;
    assert!(
        (r.mean_ratio - expected).abs() < 1e-12,
        "right isosceles should be √3/2 = {expected}, got {}",
        r.mean_ratio
    );
}

/// A sliver tends to zero, which is the other end of the scale and the end the
/// metric exists to see.
#[test]
fn mean_ratio_falls_toward_zero_as_a_triangle_collapses() {
    let cfg = ValidateConfig::from_cell_size(1.0).expect("valid cell size");
    let mut last = 1.0;
    for height in [1e-1, 1e-2, 1e-3, 1e-4] {
        let tri = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, height, 0.0]];
        let q = validate_indexed(&tri, &[0, 1, 2], &cfg).mean_ratio;
        assert!(q < last, "quality should fall as the triangle flattens");
        assert!(q > 0.0, "and stay positive while the area is positive");
        last = q;
    }
    assert!(
        last < 1e-3,
        "a 1e-4-tall sliver should be near zero, got {last}"
    );
}

/// **Valence 6 is the regular interior value, and the fixture makes that
/// visible**: a hexagonal fan has one interior vertex of valence 6 and six rim
/// vertices of valence 3 or 4, so exactly the rim is irregular.
///
/// This is also the boundary caveat in miniature — every rim vertex here is
/// irregular because the patch has a boundary, not because its triangles are
/// poor. See the field's doc comment.
#[test]
fn only_the_boundary_of_a_hexagonal_fan_is_irregular() {
    let cfg = ValidateConfig::from_cell_size(1.0).expect("valid cell size");
    let mut positions = alloc::vec![[0.0_f64, 0.0, 0.0]];
    for k in 0..6 {
        let a = core::f64::consts::TAU * f64::from(k) / 6.0;
        positions.push([a.cos(), a.sin(), 0.0]);
    }
    let mut indices = alloc::vec::Vec::new();
    for k in 0..6u32 {
        indices.extend_from_slice(&[0, 1 + k, 1 + (k + 1) % 6]);
    }

    let r = validate_indexed(&positions, &indices, &cfg);
    assert_eq!(r.referenced_vertices, 7);
    assert_eq!(
        r.irregular_vertices, 6,
        "the centre has valence 6 and the six rim vertices do not"
    );
    // Six equilateral triangles, so the quality is 1 to rounding.
    assert!((r.mean_ratio - 1.0).abs() < 1e-9, "got {}", r.mean_ratio);
}
