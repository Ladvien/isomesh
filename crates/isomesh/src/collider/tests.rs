//! G-005's acceptance, run against parry itself, plus the seam failure it warns
//! about.

use alloc::vec::Vec;

use parry3d::math::Vector;
use parry3d::shape::{TriMesh, TriMeshFlags};

use super::*;
use crate::chunk::{ChunkId, ChunkLayout};
use crate::fields::{ReferenceField, Sphere, Torus, csg_difference};
use crate::marching_cubes::MarchingCubes;
use crate::validate::ValidateConfig;
use crate::weld::Welder;
use crate::{MeshBuffer, RuntimeShape3, Sdf};

/// Mesh a field with Marching Cubes at `samples` per axis.
fn mesh<F: Sdf<Scalar = f32> + ReferenceField>(field: &F, samples: u32) -> (MeshBuffer<f32>, f32) {
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / (samples - 1) as f32;
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let mut out = MeshBuffer::<f32>::new();
    MarchingCubes::<f32>::new()
        .extract(field, &shape, lo, h, &mut out)
        .expect("extraction");
    (out, h)
}

fn config(h: f32) -> ValidateConfig {
    ValidateConfig::from_cell_size(f64::from(h)).expect("valid cell size")
}

/// The conversion a consumer writes, in the one place it is written.
fn to_parry(mesh: &MeshBuffer<f32>) -> Result<TriMesh, parry3d::shape::TriMeshBuilderError> {
    let vertices: Vec<Vector> = mesh
        .positions
        .iter()
        .map(|p| Vector::new(p[0], p[1], p[2]))
        .collect();
    TriMesh::new(vertices, triangle_indices(mesh))
}

// ─── the acceptance criterion, checked by parry ─────────────────────────────

/// **G-005's acceptance, verbatim:** *"a carved shape builds a `TriMesh` without
/// error and passes parry's own validity check."*
///
/// The carve matters. A pristine sphere is the easy case; a brush-carved field
/// has a concave rim, and that rim is where a contouring bug would put a
/// degenerate triangle or a flipped winding.
///
/// Parry's own check is `set_flags(TriMeshFlags::HALF_EDGE_TOPOLOGY)`, which
/// builds the half-edge adjacency and returns `TopologyError` if the mesh will
/// not support one. That is a stronger statement than "the constructor accepted
/// it" — `TriMesh::new` only refuses an empty index buffer.
#[test]
fn a_carved_shape_builds_a_parry_trimesh() {
    // `csg_difference` is a box with a sphere subtracted from it — a carved
    // shape, and the concave rim it leaves is where a contouring bug would put a
    // degenerate triangle or a flipped winding. A pristine sphere would not
    // exercise either.
    let field = csg_difference::<f32>();

    let (lo, hi) = field.domain();
    let samples = 41u32;
    let h = (hi[0] - lo[0]) / (samples - 1) as f32;
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let mut buffer = MeshBuffer::<f32>::new();
    MarchingCubes::<f32>::new()
        .extract(&field, &shape, lo, h, &mut buffer)
        .expect("extraction");

    let ready = readiness(&buffer, &config(h));
    std::println!(
        "carved csg_difference at {samples}^3: {} triangles, {} degenerate, {} duplicate vertices, \
         usable {}, inside/outside {}",
        ready.triangles,
        ready.degenerate_triangles,
        ready.duplicate_vertices,
        ready.is_usable(),
        ready.supports_inside_outside(),
    );
    assert!(ready.is_usable(), "{ready:?}");

    let mut trimesh = to_parry(&buffer).expect("parry accepts the mesh");
    assert_eq!(trimesh.num_triangles() as u64, ready.triangles);

    // Parry's own validity check.
    trimesh
        .set_flags(TriMeshFlags::HALF_EDGE_TOPOLOGY)
        .expect("parry builds a half-edge topology");
    assert!(trimesh.topology().is_some());

    // And the orientation flag, which is what an inside/outside query needs.
    assert!(ready.supports_inside_outside(), "{ready:?}");
    trimesh
        .set_flags(TriMeshFlags::ORIENTED)
        .expect("parry orients the mesh");
    assert!(trimesh.pseudo_normals_if_oriented().is_some());
}

/// The prediction this module is built on, tested rather than asserted: parry
/// refuses a mesh with no triangles, and accepts far more than it should.
///
/// `TriMesh::new`'s only documented failure is an empty index buffer. So a caller
/// who trusts the constructor to vet the mesh gets a `TriMesh` built from a
/// seam-ridden, degenerate-riddled buffer and finds out during play. That gap is
/// exactly what [`ColliderReadiness`] fills.
#[test]
fn parry_only_refuses_an_empty_index_buffer() {
    let empty = MeshBuffer::<f32>::new();
    assert!(to_parry(&empty).is_err(), "parry accepted an empty mesh");

    // A single degenerate triangle: three collinear points, zero area.
    let mut degenerate = MeshBuffer::<f32>::new();
    for p in [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]] {
        degenerate.positions.push(p);
        degenerate.normals.push([0.0, 0.0, 1.0]);
    }
    degenerate.indices.extend_from_slice(&[0, 1, 2]);
    assert!(
        to_parry(&degenerate).is_ok(),
        "parry refused a degenerate triangle, so this module's premise is wrong"
    );

    let ready = readiness(&degenerate, &config(1.0));
    assert_eq!(ready.degenerate_triangles, 1, "{ready:?}");
    std::println!(
        "parry accepted a zero-area triangle; readiness reports {} degenerate",
        ready.degenerate_triangles
    );
}

// ─── the seam, which is the reason the ticket has a note on it ──────────────

/// **The failure G-005's ticket warns about**, reproduced end to end: two chunks
/// meshed independently, concatenated, and handed over unwelded.
///
/// A renderer draws this correctly. Parry sees the seam as boundary edges — a
/// hole a character falls through — and `supports_inside_outside` says so before
/// anything is handed over. Welding closes it, and then it does not.
#[test]
fn an_unwelded_seam_is_reported_and_welding_fixes_it() {
    let field = Torus::<f32>::canonical();
    let layout = ChunkLayout::<f32>::new(16, 0.125, [-1.0, -1.0, -1.0]).expect("valid layout");
    let shape = layout.sample_shape().expect("valid shape");

    let mut joined = MeshBuffer::<f32>::new();
    for id in [ChunkId::new([0, 0, 0]), ChunkId::new([1, 0, 0])] {
        let mut chunk = MeshBuffer::<f32>::new();
        MarchingCubes::<f32>::new()
            .extract(
                &field,
                &shape,
                layout.sample_origin(id),
                layout.cell_size(),
                &mut chunk,
            )
            .expect("extraction");
        joined.append(&chunk);
    }

    let cfg = config(layout.cell_size());
    let before = readiness(&joined, &cfg);
    std::println!(
        "two chunks unwelded: {} duplicate vertices, {} boundary edges, seam-free {}",
        before.duplicate_vertices,
        before.boundary_edges,
        before.is_seam_free()
    );

    // The whole point of the ticket's note.
    assert!(
        before.is_usable(),
        "an unwelded seam is still structurally fine"
    );
    assert!(
        !before.is_seam_free(),
        "two chunks met with no duplicate vertices at all — the fixture is not exercising a seam"
    );

    // Parry takes it anyway. That is the trap.
    let taken = to_parry(&joined);
    assert!(
        taken.is_ok(),
        "parry refused the unwelded mesh, which would make this a non-issue"
    );

    let mut welder = Welder::<f32>::new();
    welder
        .weld(&mut joined, layout.cell_size() * 1e-4)
        .expect("weld");
    let after = readiness(&joined, &cfg);
    std::println!(
        "after welding: {} duplicate vertices, {} boundary edges, seam-free {}",
        after.duplicate_vertices,
        after.boundary_edges,
        after.is_seam_free()
    );

    assert!(after.is_seam_free(), "{after:?}");
    assert!(
        after.boundary_edges < before.boundary_edges,
        "welding did not close any boundary: {before:?} -> {after:?}"
    );
    assert!(to_parry(&joined).is_ok());
}

// ─── the interpretation itself ──────────────────────────────────────────────

/// A closed sphere supports an inside/outside query; a single open chunk does
/// not, and neither reading is a defect.
#[test]
fn only_a_closed_surface_supports_an_inside_outside_query() {
    let (whole, h) = mesh(&Sphere::<f32>::canonical(), 33);
    let ready = readiness(&whole, &config(h));
    assert!(ready.is_usable(), "{ready:?}");
    assert!(ready.supports_inside_outside(), "{ready:?}");
    assert_eq!(ready.boundary_edges, 0);

    // One chunk of a larger field leaves the surface through its sides, which is
    // the normal state of every chunk in a streamed world.
    let field = Torus::<f32>::canonical();
    let layout = ChunkLayout::<f32>::new(16, 0.125, [-1.0, -1.0, -1.0]).expect("valid layout");
    let shape = layout.sample_shape().expect("valid shape");
    let mut chunk = MeshBuffer::<f32>::new();
    MarchingCubes::<f32>::new()
        .extract(
            &field,
            &shape,
            layout.sample_origin(ChunkId::new([0, 0, 0])),
            layout.cell_size(),
            &mut chunk,
        )
        .expect("extraction");

    let ready = readiness(&chunk, &config(layout.cell_size()));
    std::println!(
        "one chunk: {} triangles, {} boundary edges, usable {}, inside/outside {}",
        ready.triangles,
        ready.boundary_edges,
        ready.is_usable(),
        ready.supports_inside_outside()
    );
    assert!(ready.is_usable(), "a chunk is still a usable collider");
    assert!(
        ready.boundary_edges > 0 && !ready.supports_inside_outside(),
        "a chunk of a larger field must be open: {ready:?}"
    );
}

/// The index buffer is regrouped, not rewritten.
#[test]
fn triangle_indices_regroups_without_changing_anything() {
    let (buffer, _) = mesh(&Sphere::<f32>::canonical(), 17);
    let triples = triangle_indices(&buffer);
    assert_eq!(triples.len(), buffer.triangle_count());
    let flat: Vec<u32> = triples.iter().flat_map(|t| t.iter().copied()).collect();
    assert_eq!(flat, buffer.indices);
}

/// A trailing partial triangle is dropped rather than completed, and the
/// readiness report names it so the caller is not left wondering where a triangle
/// went.
#[test]
fn a_trailing_partial_triangle_is_dropped_and_reported() {
    let (mut buffer, h) = mesh(&Sphere::<f32>::canonical(), 17);
    let complete = buffer.triangle_count();
    buffer.indices.push(0);
    buffer.indices.push(1);

    assert_eq!(triangle_indices(&buffer).len(), complete);
    let ready = readiness(&buffer, &config(h));
    assert_eq!(ready.trailing_indices, 2, "{ready:?}");
    assert!(
        !ready.is_usable(),
        "a trailing index means the caller built the buffer wrong: {ready:?}"
    );
}

/// Reading a report the caller already has must give the same answer as running
/// the validator again — or the two would be a second path.
#[test]
fn from_report_agrees_with_readiness() {
    let (buffer, h) = mesh(&Torus::<f32>::canonical(), 25);
    let cfg = config(h);
    let direct = readiness(&buffer, &cfg);
    let report = crate::validate::validate_indexed(&buffer.positions, &buffer.indices, &cfg);
    assert_eq!(direct, from_report(&report));
}
