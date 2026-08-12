//! T-005b — the extractors, run over generated fields on generated grids.
//!
//! [T-005a](super) built the generators and the assertion bundle; this module
//! points them at [`MarchingCubes`] and [`SurfaceNets`]. A hand-written fixture
//! tests the configurations someone thought of. These reach the ones nobody
//! did.
//!
//! # The mutation check
//!
//! A property test that cannot fail is decoration, so the suite has to be shown
//! failing. The defect worth simulating is a **corrupted case table**, because
//! that is the failure mode the whole 256-entry construction exists to prevent
//! and the one that produces meshes which look fine.
//!
//! The extractor reads [`CASES`] directly and there is no injection seam, which
//! is correct — adding a table parameter to the public API so that a test can
//! swap it would be a second execution path in production code, and the crate's
//! rule is one.
//!
//! So the corrupted table runs through [`march_with_table`], a marcher local to
//! this module. Every geometric primitive it uses is the crate's own —
//! `EDGE_CORNERS`, `EDGE_AXIS`, `corner_offset`, `is_inside`, `edge_crossing`,
//! the same gradient normalisation — so the *only* thing that can differ from
//! [`MarchingCubes::extract`] is which table it reads.
//!
//! That claim is not asserted, it is checked:
//! [`the_double_reproduces_marching_cubes`] compares the two bit-for-bit on the
//! real table across several fields and grids. Without that guard the mutation
//! check would prove nothing, because a test double that has drifted from the
//! code under test can fail for its own reasons.

use alloc::vec;

use proptest::prelude::*;

use super::{DOMAIN, SurfaceGate, assert_extracted_mesh_is_valid, convex_body, sphere_union};
use crate::cube::{corner_offset, edge_crossing};
use crate::dual_contouring::DualContouring;
use crate::marching_cubes::table::{
    CASES, EDGE_AXIS, EDGE_CORNERS, McCase, corner_inside, is_inside,
};
use crate::marching_cubes::{FaceAmbiguity, MarchingCubes};
use crate::surface_nets::SurfaceNets;
use crate::validate::{ValidateConfig, validate_indexed};
use crate::{MeshBuffer, MeshSink, RuntimeShape3, Sdf, Shape3, vec3};

/// Grid sizes for extraction: still deliberately non-cubic, but never so coarse
/// that there is no surface left to test.
///
/// [`resolution`](super::resolution) bottoms out at 2 because it exists for
/// index round-trips, where a one-cell axis is the interesting case. Here a
/// one-cell axis means every corner on that axis lies outside the object and the
/// mesh comes back empty, which is valid and proves nothing.
fn extraction_resolution() -> impl Strategy<Value = [u32; 3]> {
    (6u32..=18, 6u32..=18, 6u32..=18).prop_map(|(x, y, z)| [x, y, z])
}

/// A grid that is guaranteed to contain the generated object with the surface
/// strictly inside it.
///
/// The spacing comes from the **shortest** axis, so every axis spans at least
/// the full `2·DOMAIN` the generators live in; longer axes simply extend further
/// out into empty space. That is what makes "closed" an assertion rather than a
/// hope: every generated field is strictly positive on the grid boundary — a
/// sphere union reaches at most `1.0 + 0.9` on any axis and a convex body is
/// bounded by a radius-1.5 sphere, both inside `DOMAIN = 2.0`.
fn grid_for(size: [u32; 3]) -> (RuntimeShape3, [f64; 3], f64) {
    let shortest = size.iter().copied().min().expect("three axes");
    let cell_size = 2.0 * DOMAIN / f64::from(shortest - 1);
    (
        RuntimeShape3::new(size).expect("generated sizes fit u32"),
        [-DOMAIN; 3],
        cell_size,
    )
}

/// Extract with Marching Cubes and run the bundle over the result. Returns the
/// triangle count, so the caller can also assert the mesh is not empty.
///
/// **The strict [`SurfaceGate::Closed`] gate, since A-015** — and that is a
/// measured decision in both directions. It was waived at ✗15, which exhibited a
/// generated sphere union giving 2 non-manifold edges at `h = 2/3` and
/// attributed it to the surface pinching inside one cell. A-015 showed the cause
/// was the fan chord instead, took that fixture to zero, and the gate went back
/// on: **8,000 generated cases pass it** (`PROPTEST_CASES=4000` on each of the
/// two properties), where before A-015 it failed on the first fresh seed.
///
/// That is evidence, not proof — nothing here shows Marching Cubes is
/// *unconditionally* manifold, only that the one mechanism ever exhibited is
/// gone. If this ever fails it is a finding and not a regression: it would mean
/// a second mechanism exists, and the failing case would be the first example of
/// it. See O-12.
fn check_mc<S: Sdf<Scalar = f64>>(label: &str, field: &S, size: [u32; 3]) -> usize {
    let (shape, origin, cell_size) = grid_for(size);
    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    mc.extract(field, &shape, origin, cell_size, &mut out)
        .expect("extraction");
    assert_extracted_mesh_is_valid(
        label,
        &out.positions,
        &out.indices,
        cell_size,
        SurfaceGate::Closed,
    );
    out.triangle_count()
}

/// Extract with the asymptotic decider and run the bundle over the result.
///
/// Same gate as [`check_mc`], and the reason is the same one: the decider does
/// not change where vertices go or which grid edges carry them, only how an
/// ambiguous face pairs its four cut edges, so every argument about MC's
/// manifoldness carries over unchanged — including A-015's, since a cycle with
/// no chord-safe apex gets a centroid whichever rule produced it. Also measured
/// at 8,000 generated cases. What the generators add over the seven
/// reference fields is *reach* — `mc/tests.rs`'s census finds an ambiguous face
/// on only two of them, so without this the rule would be exercised on a
/// handful of cells in the whole suite.
fn check_mc33<S: Sdf<Scalar = f64>>(label: &str, field: &S, size: [u32; 3]) -> usize {
    let (shape, origin, cell_size) = grid_for(size);
    let mut mc = MarchingCubes::<f64>::new();
    mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
    let mut out = MeshBuffer::<f64>::new();
    mc.extract(field, &shape, origin, cell_size, &mut out)
        .expect("extraction");
    assert_extracted_mesh_is_valid(
        label,
        &out.positions,
        &out.indices,
        cell_size,
        SurfaceGate::Closed,
    );
    out.triangle_count()
}

/// Extract with Surface Nets and run the bundle over the result.
///
/// A weaker gate than [`check_mc`]'s, and the reason is measured rather than
/// assumed — see [`SurfaceGate::ClosedAllowingUnresolvedTopology`]. Surface Nets places
/// one vertex per cell, so any feature thinner than a cell forces two sheets to
/// share a vertex and the mesh is non-manifold by construction. Everything that
/// is *not* a consequence of that is still asserted.
fn check_sn<S: Sdf<Scalar = f64>>(label: &str, field: &S, size: [u32; 3]) -> usize {
    let (shape, origin, cell_size) = grid_for(size);
    let mut sn = SurfaceNets::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    sn.extract(field, &shape, origin, cell_size, &mut out)
        .expect("extraction");
    assert_extracted_mesh_is_valid(
        label,
        &out.positions,
        &out.indices,
        cell_size,
        SurfaceGate::ClosedAllowingUnresolvedTopology,
    );
    out.triangle_count()
}

/// Extract with Dual Contouring and run the bundle over the result.
///
/// Same gate as Surface Nets, and for the same reason: the two share the dual
/// topology, so Dual Contouring inherits one-vertex-per-cell and is non-manifold
/// wherever two sheets meet in a cell. What is different is that its vertex can
/// also leave its own cell, which is A-009's subject rather than this gate's.
fn check_dc<S: Sdf<Scalar = f64>>(label: &str, field: &S, size: [u32; 3]) -> usize {
    let (shape, origin, cell_size) = grid_for(size);
    let mut dc = DualContouring::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    dc.extract(field, &shape, origin, cell_size, &mut out)
        .expect("extraction");
    assert_extracted_mesh_is_valid(
        label,
        &out.positions,
        &out.indices,
        cell_size,
        SurfaceGate::ClosedAllowingUnresolvedTopology,
    );
    out.triangle_count()
}

/// The radius of the largest sphere in a generated union.
///
/// An object comfortably larger than the grid spacing **must** be found. Without
/// this the properties above pass on an empty mesh, which every gate accepts:
/// a sub-voxel sphere legitimately produces nothing, so "valid" alone is
/// satisfied by an extractor that returns nothing at all.
fn largest_radius(field: &super::SphereUnion) -> f64 {
    field
        .spheres
        .iter()
        .map(|s| s.radius)
        .fold(0.0f64, f64::max)
}

/// The spacing `grid_for` will choose for this size.
fn spacing_for(size: [u32; 3]) -> f64 {
    2.0 * DOMAIN / f64::from(size.iter().copied().min().expect("three axes") - 1)
}

proptest! {
    // 256 rather than T-005a's 1000, because a case here extracts, validates and
    // self-intersection-scans a whole mesh instead of evaluating a field at two
    // points. Measured: all four properties together run in 0.87 s in a debug
    // build, so this is a considered number rather than a guess, and there is
    // room to raise it if a regression ever needs hunting.
    #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

    #[test]
    fn marching_cubes_meshes_sphere_unions(field in sphere_union(), size in extraction_resolution()) {
        let tris = check_mc("mc / sphere union", &field, size);
        if largest_radius(&field) >= 2.0 * spacing_for(size) {
            prop_assert!(tris > 0, "a resolvable sphere was missed entirely");
        }
    }

    #[test]
    fn marching_cubes_meshes_convex_bodies(field in convex_body(), size in extraction_resolution()) {
        check_mc("mc / convex body", &field, size);
    }

    #[test]
    fn the_decider_meshes_sphere_unions(field in sphere_union(), size in extraction_resolution()) {
        let tris = check_mc33("mc33 / sphere union", &field, size);
        if largest_radius(&field) >= 2.0 * spacing_for(size) {
            prop_assert!(tris > 0, "a resolvable sphere was missed entirely");
        }
    }

    /// Unions of spheres are where an ambiguous face actually turns up: two
    /// lobes approaching each other put diagonally opposite corners inside one
    /// face. The reference fields barely produce the configuration at all.
    #[test]
    fn the_decider_meshes_convex_bodies(field in convex_body(), size in extraction_resolution()) {
        check_mc33("mc33 / convex body", &field, size);
    }

    #[test]
    fn surface_nets_meshes_sphere_unions(field in sphere_union(), size in extraction_resolution()) {
        let tris = check_sn("sn / sphere union", &field, size);
        if largest_radius(&field) >= 2.0 * spacing_for(size) {
            prop_assert!(tris > 0, "a resolvable sphere was missed entirely");
        }
    }

    #[test]
    fn surface_nets_meshes_convex_bodies(field in convex_body(), size in extraction_resolution()) {
        check_sn("sn / convex body", &field, size);
    }

    #[test]
    fn dual_contouring_meshes_sphere_unions(field in sphere_union(), size in extraction_resolution()) {
        let tris = check_dc("dc / sphere union", &field, size);
        if largest_radius(&field) >= 2.0 * spacing_for(size) {
            prop_assert!(tris > 0, "a resolvable sphere was missed entirely");
        }
    }

    /// The generator that matters most for this algorithm: a convex body is an
    /// intersection of half-spaces, so it is *all* sharp edges and corners --
    /// exactly the geometry the solve exists for and the geometry most likely to
    /// send a vertex somewhere absurd.
    #[test]
    fn dual_contouring_meshes_convex_bodies(field in convex_body(), size in extraction_resolution()) {
        check_dc("dc / convex body", &field, size);
    }
}

// ─── the injectable marcher ─────────────────────────────────────────────────

/// [`MarchingCubes::extract`]'s loop, reading a caller-supplied case table.
///
/// Test-only, and it exists for exactly one reason: to run a *corrupted* table
/// through otherwise-identical machinery. See the module docs for why the real
/// extractor does not take a table parameter, and
/// [`the_double_reproduces_marching_cubes`] for the guard that keeps this
/// honest.
pub(crate) fn march_with_table<S: Sdf<Scalar = f64>>(
    sdf: &S,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    cell_size: f64,
    cases: &[McCase; 256],
) -> MeshBuffer<f64> {
    let size = shape.size();

    let mut values = vec![0.0f64; shape.element_count()];
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let i = shape.linearize([x, y, z]) as usize;
                values[i] = sdf.sample(corner_world([x, y, z], [0, 0, 0], origin, cell_size));
            }
        }
    }

    let mut out = MeshBuffer::<f64>::new();
    let mut edge_vertices = vec![u32::MAX; shape.element_count() * 3];

    for z in 0..size[2] - 1 {
        for y in 0..size[1] - 1 {
            for x in 0..size[0] - 1 {
                let base = [x, y, z];

                let mut case = 0u8;
                let mut corner_value = [0.0f64; 8];
                for c in 0..8u8 {
                    let s = corner_linear(shape, base, c);
                    let v = values[s as usize];
                    corner_value[c as usize] = v;
                    if is_inside(v) {
                        case |= 1 << c;
                    }
                }

                let entry = &cases[case as usize];
                if entry.count == 0 {
                    continue;
                }

                for tri in &entry.triangles[..entry.count as usize] {
                    let mut idx = [0u32; 3];
                    for (k, &edge) in tri.iter().enumerate() {
                        let [lo_corner, hi_corner] = EDGE_CORNERS[edge as usize];
                        let axis = EDGE_AXIS[edge as usize] as usize;
                        let lo_sample = corner_linear(shape, base, lo_corner);
                        let key = lo_sample as usize * 3 + axis;

                        if edge_vertices[key] != u32::MAX {
                            idx[k] = edge_vertices[key];
                            continue;
                        }

                        let a = corner_value[lo_corner as usize];
                        let b = corner_value[hi_corner as usize];
                        let t = edge_crossing(a, b);
                        let lo_pos =
                            corner_world(base, corner_offset(lo_corner), origin, cell_size);
                        let hi_pos =
                            corner_world(base, corner_offset(hi_corner), origin, cell_size);
                        let position = [
                            lo_pos[0] + (hi_pos[0] - lo_pos[0]) * t,
                            lo_pos[1] + (hi_pos[1] - lo_pos[1]) * t,
                            lo_pos[2] + (hi_pos[2] - lo_pos[2]) * t,
                        ];

                        let g = sdf.gradient(position);
                        let len = vec3::length(g);
                        let normal = vec3::scale(g, len.recip());

                        let index = out.vertex(position, normal);
                        edge_vertices[key] = index;
                        idx[k] = index;
                    }
                    out.triangle(idx[0], idx[1], idx[2]);
                }
            }
        }
    }

    out
}

/// World position of `base + offset`.
fn corner_world(base: [u32; 3], offset: [u32; 3], origin: [f64; 3], cell_size: f64) -> [f64; 3] {
    [
        origin[0] + cell_size * f64::from(base[0] + offset[0]),
        origin[1] + cell_size * f64::from(base[1] + offset[1]),
        origin[2] + cell_size * f64::from(base[2] + offset[2]),
    ]
}

/// Linear index of one corner of one cell.
fn corner_linear(shape: &RuntimeShape3, base: [u32; 3], corner: u8) -> u32 {
    let o = corner_offset(corner);
    shape.linearize([base[0] + o[0], base[1] + o[1], base[2] + o[2]])
}

/// The guard the mutation check rests on.
///
/// If this ever fails, [`march_with_table`] has drifted from
/// [`MarchingCubes::extract`] and every conclusion drawn from a corrupted table
/// is worthless — the mesh would differ for a reason that has nothing to do with
/// the corruption.
#[test]
fn the_double_reproduces_marching_cubes() {
    let field = crate::fields::Sphere::<f64>::canonical();
    for size in [[9, 9, 9], [12, 7, 15], [6, 18, 11]] {
        let (shape, origin, cell_size) = grid_for(size);
        let mut mc = MarchingCubes::<f64>::new();
        let mut real = MeshBuffer::<f64>::new();
        mc.extract(&field, &shape, origin, cell_size, &mut real)
            .expect("extraction");

        let double = march_with_table(&field, &shape, origin, cell_size, &CASES);

        assert_eq!(double.indices, real.indices, "{size:?}: indices differ");
        assert_eq!(
            double.positions.len(),
            real.positions.len(),
            "{size:?}: vertex count differs"
        );
        for (i, (d, r)) in double.positions.iter().zip(&real.positions).enumerate() {
            for a in 0..3 {
                assert_eq!(
                    d[a].to_bits(),
                    r[a].to_bits(),
                    "{size:?}: vertex {i} axis {a}"
                );
            }
        }
        for (i, (d, r)) in double.normals.iter().zip(&real.normals).enumerate() {
            for a in 0..3 {
                assert_eq!(
                    d[a].to_bits(),
                    r[a].to_bits(),
                    "{size:?}: normal {i} axis {a}"
                );
            }
        }
        assert!(
            !real.indices.is_empty(),
            "{size:?}: fixture produced nothing"
        );
    }
}

// ─── the mutation check ─────────────────────────────────────────────────────

/// A grid and field the mutation checks share. Fine enough that every low-count
/// case appears many times over.
fn mutation_fixture() -> (crate::fields::Sphere<f64>, RuntimeShape3, [f64; 3], f64) {
    let (shape, origin, cell_size) = grid_for([15, 15, 15]);
    (crate::fields::Sphere::canonical(), shape, origin, cell_size)
}

/// Which case indices the fixture actually exercises.
///
/// A corruption applied to a case the mesh never reaches would change nothing
/// and the mutation check would pass while testing nothing — the exact failure
/// mode this whole module is about.
fn cases_used(
    shape: &RuntimeShape3,
    field: &crate::fields::Sphere<f64>,
    origin: [f64; 3],
    h: f64,
) -> [u32; 256] {
    let size = shape.size();
    let mut counts = [0u32; 256];
    for z in 0..size[2] - 1 {
        for y in 0..size[1] - 1 {
            for x in 0..size[0] - 1 {
                let mut case = 0u8;
                for c in 0..8u8 {
                    let o = corner_offset(c);
                    let p = corner_world([x, y, z], o, origin, h);
                    if is_inside(field.sample(p)) {
                        case |= 1 << c;
                    }
                }
                counts[case as usize] += 1;
            }
        }
    }
    counts
}

/// Reversing one triangle's winding in one case is the classic transcription
/// error, and it is invisible to every check except edge orientation: the mesh
/// still has the right Euler characteristic, the right edge degrees and the
/// right vertex links.
#[test]
#[should_panic(expected = "expected a closed surface")]
fn a_flipped_triangle_in_the_case_table_is_caught() {
    let (field, shape, origin, cell_size) = mutation_fixture();
    let used = cases_used(&shape, &field, origin, cell_size);

    let mut cases = CASES;
    let victim = (1..255u16)
        .map(|c| c as usize)
        .find(|&c| used[c] > 0 && cases[c].count > 0)
        .expect("the fixture must exercise some non-empty case");
    cases[victim].triangles[0].swap(1, 2);

    let mesh = march_with_table(&field, &shape, origin, cell_size, &cases);
    assert_extracted_mesh_is_valid(
        "flipped triangle",
        &mesh.positions,
        &mesh.indices,
        cell_size,
        SurfaceGate::Closed,
    );
}

/// Dropping a triangle from one case leaves a hole, which shows up as boundary
/// edges rather than as an orientation error — a different detector from the
/// one above, so the two together show the bundle is not resting on a single
/// counter.
#[test]
#[should_panic(expected = "expected a closed surface")]
fn a_missing_triangle_in_the_case_table_is_caught() {
    let (field, shape, origin, cell_size) = mutation_fixture();
    let used = cases_used(&shape, &field, origin, cell_size);

    let mut cases = CASES;
    let victim = (1..255u16)
        .map(|c| c as usize)
        .find(|&c| used[c] > 0 && cases[c].count > 0)
        .expect("the fixture must exercise some non-empty case");
    cases[victim].count -= 1;

    let mesh = march_with_table(&field, &shape, origin, cell_size, &cases);
    assert_extracted_mesh_is_valid(
        "missing triangle",
        &mesh.positions,
        &mesh.indices,
        cell_size,
        SurfaceGate::Closed,
    );
}

/// Repointing one triangle at a **different but still-cut** edge of the same
/// cell: the transcription error where a table entry names the wrong edge.
///
/// The "still-cut" part is the whole difficulty, and it is worth stating because
/// the naive version of this test does not test what it looks like. Pointing at
/// an *uncut* edge trips `edge_crossing`'s own precondition
/// (`is_inside(a) != is_inside(b)`) and panics inside the crate before a mesh
/// exists — which is a real and welcome defence, but it means the assertion
/// bundle was never reached and the test would have proved nothing about it.
/// Worse, that guard is a `debug_assert` and vanishes in a release build.
///
/// So the mutation is confined to edges the case actually cuts, which is exactly
/// the plausible transcription error, and the resulting mesh is well formed
/// enough to reach the bundle and be rejected on topology.
#[test]
#[should_panic(expected = "expected a closed surface")]
fn a_wrong_edge_reference_in_the_case_table_is_caught() {
    let (field, shape, origin, cell_size) = mutation_fixture();
    let used = cases_used(&shape, &field, origin, cell_size);

    let mut cases = CASES;
    let (victim, replacement) = (1..255u16)
        .map(|c| c as u8)
        .filter(|&c| used[c as usize] > 0 && cases[c as usize].count > 0)
        .find_map(|c| {
            let tri = cases[c as usize].triangles[0];
            // A cut edge of this case that this triangle does not already name.
            (0..12u8)
                .find(|&e| {
                    let [lo, hi] = EDGE_CORNERS[e as usize];
                    corner_inside(c, lo) != corner_inside(c, hi) && !tri.contains(&e)
                })
                .map(|e| (c as usize, e))
        })
        .expect("some exercised case cuts an edge its first triangle does not name");
    cases[victim].triangles[0][0] = replacement;

    let mesh = march_with_table(&field, &shape, origin, cell_size, &cases);
    assert_extracted_mesh_is_valid(
        "wrong edge reference",
        &mesh.positions,
        &mesh.indices,
        cell_size,
        SurfaceGate::Closed,
    );
}

/// And the control: the *uncorrupted* table through the same path passes.
///
/// Without this the three tests above would also pass if `march_with_table`
/// produced garbage unconditionally.
#[test]
fn the_uncorrupted_table_passes_the_same_bundle() {
    let (field, shape, origin, cell_size) = mutation_fixture();
    let mesh = march_with_table(&field, &shape, origin, cell_size, &CASES);
    let report = validate_indexed(
        &mesh.positions,
        &mesh.indices,
        &ValidateConfig::from_cell_size(cell_size).expect("valid cell size"),
    );
    assert!(report.is_closed(), "the control must pass\n{report}");
    assert_extracted_mesh_is_valid(
        "uncorrupted control",
        &mesh.positions,
        &mesh.indices,
        cell_size,
        SurfaceGate::Closed,
    );
}

/// Marching Cubes is **not** unconditionally manifold, and this is the case
/// that showed it.
///
/// Found by `marching_cubes_meshes_sphere_unions` on a fresh proptest seed, then
/// reduced. Three spheres whose union pinches inside a single cell at
/// `h = 2/3`: the surface touches itself, the shared grid edge ends up used by
/// four faces, and the result is closed, correctly oriented, `χ = 2` — and
/// non-manifold.
///
/// Pinned in both directions, following M-4's precedent of asserting a known
/// defect as a non-zero count rather than excluding it: the counts at 7³ are
/// exact, and every refinement from 9³ up is exactly zero. So this fails if the
/// defect ever spreads *and* if it ever silently disappears.
#[test]
fn the_fixture_that_falsified_unconditional_manifoldness_is_now_manifold() {
    use crate::fields::Sphere;

    let field = super::SphereUnion {
        spheres: alloc::vec![
            Sphere {
                center: [
                    0.216_424_612_766_318_28,
                    0.529_307_710_262_215_2,
                    -0.804_663_039_989_917_6
                ],
                radius: 0.619_553_810_790_568_1,
            },
            Sphere {
                center: [0.514_601_202_644_422_7, 0.230_953_855_883_975_85, 0.0],
                radius: 0.495_969_042_463_108_13,
            },
            Sphere {
                center: [0.449_324_060_565_480_9, -0.870_601_428_657_975_9, 0.0],
                radius: 0.875_530_864_840_149_2,
            },
        ],
    };

    let report_at = |n: u32| {
        let (shape, origin, h) = grid_for([n; 3]);
        let mut mc = MarchingCubes::<f64>::new();
        let mut out = MeshBuffer::<f64>::new();
        mc.extract(&field, &shape, origin, h, &mut out)
            .expect("extraction");
        validate_indexed(
            &out.positions,
            &out.indices,
            &ValidateConfig::from_cell_size(h).expect("valid cell size"),
        )
    };

    // Was 2 non-manifold edges and 3 non-manifold vertices at n = 7, which is
    // what falsified "Marching Cubes is unconditionally manifold" (✗15). A-015
    // took it to zero, and that is the finding: a **triangulation** change
    // cannot repair two sheets genuinely meeting inside a cell, so the cause was
    // never the geometry ✗15 attributed it to. It was the fan chord.
    let coarse = report_at(7);
    assert_eq!(coarse.non_manifold_edges, 0, "{coarse}");
    assert_eq!(coarse.non_manifold_vertices, 0, "{coarse}");
    assert_eq!(coarse.boundary_edges, 0, "{coarse}");
    assert_eq!(coarse.inconsistently_oriented_edges, 0, "{coarse}");
    // And chi moves 2 -> 0, which is the same defect seen from the other side.
    // A collided mesh edge is counted **once** in E though it carries four
    // faces, so E was short by exactly the two collisions and chi was long by
    // two. The old value was not a topology reading at all. At this spacing the
    // three lobes genuinely meet in a handle: closed, orientable, genus 1.
    assert_eq!(coarse.euler_characteristic, 0, "{coarse}");
    assert_eq!(coarse.genus, Some(1), "{coarse}");
    assert!(coarse.is_closed(), "{coarse}");

    for n in [9u32, 13, 17, 25, 33, 49, 65] {
        let r = report_at(n);
        assert_eq!(r.non_manifold_edges, 0, "n={n}\n{r}");
        assert_eq!(r.non_manifold_vertices, 0, "n={n}\n{r}");
        assert!(r.is_closed(), "n={n}\n{r}");
    }
}
