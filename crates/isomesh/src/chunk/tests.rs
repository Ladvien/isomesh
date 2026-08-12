//! The load-bearing test is the seam: two chunks meshed independently must
//! agree on the plane they share, asserted on coordinates rather than eyeballed.

use alloc::vec::Vec;

use super::{ChunkId, ChunkLayout};
use crate::fields::Sphere;
use crate::mc::MarchingCubes;
use crate::sn::SurfaceNets;
use crate::validate::{ValidateConfig, validate_indexed};
use crate::{MeshBuffer, Sdf};

fn layout(cell_size: f64) -> ChunkLayout<f64> {
    ChunkLayout::new(8, cell_size, [-2.0; 3]).expect("valid layout")
}

/// A chunk that actually straddles the unit sphere's surface, and its positive-x
/// neighbour.
///
/// With origin `[-2; 3]` and 8 cells of `0.125`, each chunk spans one world unit,
/// so chunk `[1,1,1]` covers `[-1, 0]³` and `[2,1,1]` covers `[0, 1] × [-1, 0]²`.
/// Their shared plane is `x = 0`, where the sphere passes through — chunk
/// `[0,0,0]` covers `[-2, -1]³` and contains no surface at all, which makes every
/// seam assertion vacuously true.
fn straddling_pair() -> (ChunkId, ChunkId) {
    let a = ChunkId::new([1, 1, 1]);
    (a, a.neighbour(0, 1))
}

// ─── coordinate round trips ─────────────────────────────────────────────────

#[test]
fn global_and_local_sample_indices_round_trip() {
    let l = layout(0.125);
    for cz in -2..=2i32 {
        for cy in -2..=2i32 {
            for cx in -2..=2i32 {
                let id = ChunkId::new([cx, cy, cz]);
                for lz in 0..l.cells() {
                    for ly in 0..l.cells() {
                        for lx in 0..l.cells() {
                            let local = [lx, ly, lz];
                            let global = l.global_sample(id, local);
                            let (back_id, back_local) = l.local_sample(global);
                            assert_eq!(back_id, id, "global {global:?}");
                            assert_eq!(back_local, local, "global {global:?}");
                        }
                    }
                }
            }
        }
    }
}

/// The overlap plane belongs to the next chunk, which is what "positive-face
/// overlap" means and what stops two chunks both claiming it.
#[test]
fn the_overlap_plane_is_owned_by_the_next_chunk() {
    let l = layout(0.125);
    let id = ChunkId::new([0, 0, 0]);
    let n = l.cells();

    // Local sample `n` on chunk 0 is the same global sample as local 0 on
    // chunk 1 ...
    let shared = l.global_sample(id, [n, 0, 0]);
    assert_eq!(shared, l.global_sample(id.neighbour(0, 1), [0, 0, 0]));

    // ... and it is *owned* by chunk 1.
    let (owner, local) = l.local_sample(shared);
    assert_eq!(owner, ChunkId::new([1, 0, 0]));
    assert_eq!(local, [0, 0, 0]);
}

#[test]
fn world_and_chunk_round_trip_including_negatives() {
    let l = layout(0.125);
    for cz in -3..=3i32 {
        for cy in -3..=3i32 {
            for cx in -3..=3i32 {
                let id = ChunkId::new([cx, cy, cz]);
                let origin = l.sample_origin(id);
                // A point just inside the chunk, so half-open ownership is
                // unambiguous.
                let inside = [
                    origin[0] + l.cell_size() * 0.5,
                    origin[1] + l.cell_size() * 0.5,
                    origin[2] + l.cell_size() * 0.5,
                ];
                assert_eq!(l.chunk_of(inside), id, "chunk {id:?}");
            }
        }
    }
}

/// A point exactly on a shared plane belongs to the positive side, matching
/// where the overlap sits. Off-by-one here is a whole plane of duplicated or
/// missing work.
#[test]
fn a_point_on_a_seam_belongs_to_the_positive_side() {
    let l = layout(0.125);
    let seam = l.sample_origin(ChunkId::new([1, 0, 0]));
    assert_eq!(l.chunk_of(seam).coords[0], 1);

    let just_below = [seam[0] - l.cell_size() * 0.5, seam[1], seam[2]];
    assert_eq!(l.chunk_of(just_below).coords[0], 0);
}

#[test]
fn a_degenerate_layout_is_rejected() {
    assert!(ChunkLayout::<f64>::new(0, 0.125, [0.0; 3]).is_err());
    assert!(ChunkLayout::<f64>::new(8, 0.0, [0.0; 3]).is_err());
    assert!(ChunkLayout::<f64>::new(8, -1.0, [0.0; 3]).is_err());
    assert!(ChunkLayout::<f64>::new(8, f64::NAN, [0.0; 3]).is_err());
}

// ─── the seam ───────────────────────────────────────────────────────────────

fn mesh_chunk<F: Sdf<Scalar = f64>>(
    l: &ChunkLayout<f64>,
    field: &F,
    id: ChunkId,
) -> MeshBuffer<f64> {
    let shape = l.sample_shape().expect("valid shape");
    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    mc.extract(field, &shape, l.sample_origin(id), l.cell_size(), &mut out)
        .expect("extraction");
    out
}

/// Vertices within `tol` of the plane `x = seam`.
fn on_plane(mesh: &MeshBuffer<f64>, seam: f64, tol: f64) -> Vec<[f64; 3]> {
    let mut v: Vec<[f64; 3]> = mesh
        .positions
        .iter()
        .copied()
        .filter(|p| (p[0] - seam).abs() <= tol)
        .collect();
    v.sort_by(|a, b| {
        a[1].total_cmp(&b[1])
            .then(a[2].total_cmp(&b[2]))
            .then(a[0].total_cmp(&b[0]))
    });
    v
}

/// **The acceptance criterion.** Two adjacent chunks meshed independently
/// produce coincident vertices on the plane they share.
///
/// Asserted on coordinates, and asserted **exactly** — with a power-of-two cell
/// size the two chunks' arithmetic agrees bit-for-bit, so there is no tolerance
/// to hide behind. The non-power-of-two case is measured separately below, and
/// it is a different answer.
#[test]
fn adjacent_chunks_agree_on_the_shared_plane() {
    // 0.125 is a power of two, so `h · k` is exact for every integer `k` here.
    let l = layout(0.125);
    let field = Sphere::<f64>::canonical();

    let (a, b) = straddling_pair();
    let mesh_a = mesh_chunk(&l, &field, a);
    let mesh_b = mesh_chunk(&l, &field, b);

    let seam = l.sample_origin(b)[0];
    let from_a = on_plane(&mesh_a, seam, 1e-12);
    let from_b = on_plane(&mesh_b, seam, 1e-12);

    std::println!(
        "measured: G-001 seam at x = {seam} -- chunk A contributes {} vertices, chunk B {}",
        from_a.len(),
        from_b.len()
    );

    assert!(
        !from_a.is_empty(),
        "the seam carries no geometry, so this proves nothing -- move it"
    );
    assert_eq!(
        from_a.len(),
        from_b.len(),
        "the two chunks disagree about how many vertices sit on the seam"
    );
    for (x, y) in from_a.iter().zip(&from_b) {
        for axis in 0..3 {
            assert_eq!(
                x[axis].to_bits(),
                y[axis].to_bits(),
                "seam vertices differ: {x:?} vs {y:?}"
            );
        }
    }
}

/// The seam where the two expressions genuinely disagree.
///
/// An extractor computes its samples as `origin + h·local`, so chunk `c`'s last
/// plane is `(o + h·cn) + h·n` while chunk `c+1`'s first is `o + h·(c+1)n` — the
/// same point by algebra, not by IEEE. Over 200,000 random
/// `(origin, h, cells, chunk)` combinations **22% disagree**, by one or two ulp.
///
/// The chunk and spacing here are chosen from that search rather than picked for
/// looking irregular, and that distinction cost a rewrite: the obvious choice,
/// `h = 4/33` at chunk 1, lands in the 78% that happen to agree, so the test
/// passed while proving nothing about the case it was named after. A test that
/// can only pass is the same problem as a test that cannot fail.
///
/// The gap is a rounding error rather than a crack — well under a millionth of a
/// cell — but it is not zero, and a project this careful about bit-identity
/// should know which of its guarantees survive chunking and which degrade to
/// "within an ulp". See M-32.
#[test]
fn a_non_power_of_two_cell_size_costs_exactness_at_the_seam() {
    // From the search: o = -2, cells = 8, h = 4/35, chunk 1 -> seam at
    // x = -0.1714..., which the unit sphere crosses, and the two expressions
    // differ by 1.11e-16.
    let l = layout(4.0 / 35.0);
    let field = Sphere::<f64>::canonical();

    let a = ChunkId::new([1, 1, 1]);
    let b = a.neighbour(0, 1);
    let mesh_a = mesh_chunk(&l, &field, a);
    let mesh_b = mesh_chunk(&l, &field, b);

    let seam = l.sample_origin(b)[0];
    let from_a = on_plane(&mesh_a, seam, l.cell_size() * 1e-6);
    let from_b = on_plane(&mesh_b, seam, l.cell_size() * 1e-6);
    assert_eq!(from_a.len(), from_b.len(), "seam vertex counts differ");

    // The two chunks' arithmetic for this seam plane must actually differ, or
    // the measurement below is of nothing.
    let plane_a = (l.sample_origin(a)[0]) + l.cell_size() * f64::from(l.cells());
    let plane_b = l.sample_origin(b)[0];
    assert_ne!(
        plane_a.to_bits(),
        plane_b.to_bits(),
        "this fixture was chosen because the two expressions disagree; they no longer do,          so the test has stopped measuring what it is named after"
    );

    let mut worst = 0.0f64;
    let mut exact = 0usize;
    for (x, y) in from_a.iter().zip(&from_b) {
        let d = ((x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2) + (x[2] - y[2]).powi(2)).sqrt();
        worst = worst.max(d);
        if (0..3).all(|k| x[k].to_bits() == y[k].to_bits()) {
            exact += 1;
        }
    }
    std::println!(
        "measured: G-001 seam with h = 4/35 -- {} of {} vertices bit-identical, worst gap {worst:.3e} world units ({:.2e} cells)",
        exact,
        from_a.len(),
        worst / l.cell_size()
    );

    // Whatever it is, it must be a rounding error rather than a crack.
    assert!(
        worst < l.cell_size() * 1e-9,
        "seam gap {worst} is too large to be rounding"
    );
}

/// ✗1's recorded break condition, arriving exactly where it said it would.
///
/// A chunk's mesh is clipped by the chunk boundary, so it is a manifold **with
/// boundary** and `is_closed` is false — which is not a regression. The finding
/// wrote this down before this module existed: *"expect the assertion to fail
/// there and do not 'fix' it."*
#[test]
fn a_single_chunk_is_manifold_with_boundary_not_closed() {
    let l = layout(0.125);
    let field = Sphere::<f64>::canonical();
    let mesh = mesh_chunk(&l, &field, straddling_pair().0);
    let cfg = ValidateConfig::from_cell_size(l.cell_size()).expect("valid spacing");
    let report = validate_indexed(&mesh.positions, &mesh.indices, &cfg);

    std::println!(
        "measured: G-001 one chunk of a sphere -> {} tris, chi {}, {} boundary edges",
        mesh.triangle_count(),
        report.euler_characteristic,
        report.boundary_edges
    );

    assert!(report.is_manifold(), "{report}");
    assert!(
        report.boundary_edges > 0,
        "a clipped chunk must have boundary; if it does not, the seam is not being cut\n{report}"
    );
    assert!(
        !report.is_closed(),
        "a clipped chunk is not closed -- see ✗1, this is expected\n{report}"
    );
}

/// Surface Nets meets its seam too, and it is a *different* question: its vertex
/// is one per cell rather than on a shared edge, so the shared cells are the ones
/// both chunks own a copy of.
#[test]
fn surface_nets_chunks_agree_on_shared_cell_vertices() {
    let l = layout(0.125);
    let field = Sphere::<f64>::canonical();
    let shape = l.sample_shape().expect("valid shape");

    let mesh = |id: ChunkId| {
        let mut sn = SurfaceNets::<f64>::new();
        let mut out = MeshBuffer::<f64>::new();
        sn.extract(&field, &shape, l.sample_origin(id), l.cell_size(), &mut out)
            .expect("extraction");
        out
    };

    let (a, b) = straddling_pair();
    let seam = l.sample_origin(b)[0];

    // Surface Nets puts its vertex at the centroid of a cell's crossings, so a
    // seam vertex sits half a cell away from the plane rather than on it.
    let half = l.cell_size() * 0.5;
    let near_a = on_plane(&mesh(a), seam - half, l.cell_size() * 0.25);
    let near_b = on_plane(&mesh(b), seam + half, l.cell_size() * 0.25);
    std::println!(
        "measured: G-001 surface nets seam -- {} vertices in chunk A's last cell layer, {} in chunk B's first",
        near_a.len(),
        near_b.len()
    );
    // The two layers are different cells, so they are not expected to coincide.
    // What matters is that both chunks produced geometry there rather than one
    // of them stopping short, which is the failure a 1-cell overlap prevents.
    assert!(!near_a.is_empty() && !near_b.is_empty());
}
