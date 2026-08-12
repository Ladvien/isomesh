//! The load-bearing test is the seam: two chunks meshed independently must
//! agree on the plane they share, asserted on coordinates rather than eyeballed.

use alloc::vec::Vec;

use super::{ChunkId, ChunkLayout};
use crate::fields::Sphere;
use crate::marching_cubes::MarchingCubes;
use crate::surface_nets::SurfaceNets;
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

/// **`cell_of` inverts `world_of_sample` in a cell's interior, and not reliably
/// on its corner — which is M-32 wearing a different hat.**
///
/// `world_of_sample` computes `origin + h·sample` and `cell_of` computes
/// `floor((p − origin) / h)`. Those are inverse by algebra and not by IEEE: at
/// `h = 4/35` the round trip through a corner of cell `[7, -3, 12]` comes back
/// `[6, -4, 12]`, because the division lands a hair under the integer and
/// `floor` takes it down. At `h = 0.125` it round-trips exactly, and that is the
/// same power-of-two/not divide M-32 measured at chunk seams.
///
/// No epsilon is added to paper over it. A point *exactly* on a cell boundary
/// belongs to either cell by convention, and at a non-power-of-two spacing
/// "exactly on" is not a decidable question — snapping would trade a visible
/// ambiguity for an invisible one. Callers that need a cell *range* pad it,
/// which is what E-202 does with the brush radius.
///
/// Recorded in both directions: the interior is asserted, the corner is
/// measured, so this fails if the interior ever breaks *and* if the corner
/// silently starts working.
#[test]
fn cell_of_inverts_world_of_sample_inside_a_cell() {
    let cells = [[0i64, 0, 0], [7, -3, 12], [-19, 40, -1]];
    // `power_of_two` rather than comparing the spacing back to a literal: the
    // property is about the divide, not about a particular number.
    for (spacing, name, power_of_two) in [(0.125f64, "0.125", true), (4.0 / 35.0, "4/35", false)] {
        let l = layout(spacing);
        let mut exact_corners = 0usize;
        for cell in cells {
            let corner = l.world_of_sample(cell);
            if l.cell_of(corner) == cell {
                exact_corners += 1;
            }
            // The interior is unambiguous at any spacing, and is what a caller
            // converting a bounding box actually needs.
            let inside = [
                corner[0] + l.cell_size() * 0.5,
                corner[1] + l.cell_size() * 0.5,
                corner[2] + l.cell_size() * 0.5,
            ];
            assert_eq!(l.cell_of(inside), cell, "h = {name}, interior of {cell:?}");
        }
        std::println!(
            "measured: h = {name} -- {exact_corners} of {} cell corners round-trip exactly",
            cells.len()
        );
        let expected = if power_of_two { cells.len() } else { 1 };
        assert_eq!(exact_corners, expected, "h = {name}");
    }
}

#[test]
fn cell_of_agrees_with_chunk_of() {
    let l = layout(0.125);
    let n = i64::from(l.cells());
    for point in [[0.0f64, 0.0, 0.0], [1.3, -0.7, 2.9], [-3.1, 0.2, -0.4]] {
        let cell = l.cell_of(point);
        let expected = l.chunk_of(point);
        let from_cell = crate::chunk::ChunkId::new([
            cell[0].div_euclid(n) as i32,
            cell[1].div_euclid(n) as i32,
            cell[2].div_euclid(n) as i32,
        ]);
        assert_eq!(from_cell, expected, "{point:?}");
    }
}

#[test]
fn cell_of_does_not_wrap_on_a_non_finite_point() {
    let l = layout(0.125);
    assert_eq!(l.cell_of([f64::NAN, f64::INFINITY, 0.0])[0], 0);
    assert_eq!(l.cell_of([f64::NAN, f64::INFINITY, 0.0])[1], 0);
}

// ─── G-004: field-derived LOD ───────────────────────────────────────────────

/// **The property that makes field-derived LOD free**, and the reason it is a
/// property rather than a hope.
///
/// A level-`k` sample must land **bit-identically** on the level-0 sample `2^k`
/// times its index. Doubling a float is exact in IEEE and so is doubling a small
/// integer, so `(h·2^k)·s` and `h·(2^k·s)` are the same real number rounded the
/// same way — but M-32 and M-49 both caught this crate assuming an algebraic
/// identity that IEEE did not honour, so it is asserted rather than argued.
///
/// If it failed, an LOD boundary would show a crack from coordinate drift alone,
/// before any transition-cell work (A-011b) had a chance to be wrong.
#[test]
fn every_coarse_sample_lands_exactly_on_a_fine_one() {
    // Exact comparison on purpose: this pins a bit-level coincidence.
    #![allow(clippy::float_cmp)]
    for cell_size in [0.125, 4.0 / 35.0, 0.1, 1.0 / 3.0] {
        let fine = layout(cell_size);
        for level in 0..4u32 {
            let coarse = fine.at_lod(level).expect("valid level");
            let step = 1i64 << level;
            for s in [-9i64, -1, 0, 1, 5, 33, 128] {
                let sample = [s, s + 2, s - 3];
                let at_coarse = coarse.world_of_sample(sample);
                let at_fine =
                    fine.world_of_sample([sample[0] * step, sample[1] * step, sample[2] * step]);
                assert_eq!(
                    at_coarse, at_fine,
                    "h={cell_size} level={level} sample={sample:?}"
                );
            }
        }
    }
}

/// Level `k` doubles the spacing `k` times and keeps the cell count, so one
/// chunk covers `8^k` times the volume — which is the acceptance criterion's
/// "LOD *k* has roughly 1/8^k the cells" seen from the other side: covering a
/// fixed world takes that many fewer cells.
#[test]
fn each_level_covers_eight_times_the_volume_of_the_last() {
    // Exact: doubling is exact, so this is not an approximate claim.
    #![allow(clippy::float_cmp)]
    let base = layout(0.125);
    for level in 0..4u32 {
        let lod = base.at_lod(level).expect("valid level");
        assert_eq!(lod.cells(), base.cells(), "cell count must not change");
        assert_eq!(
            lod.cell_size(),
            base.cell_size() * f64::from(1u32 << level),
            "level {level}"
        );

        // Cells needed to cover one fixed world extent, per axis.
        let extent = 8.0f64;
        let per_axis = extent / lod.cell_size();
        let expected = extent / base.cell_size() / f64::from(1u32 << level);
        assert_eq!(per_axis, expected, "level {level}");
    }
}

/// An unrepresentable level is refused rather than saturated.
#[test]
fn a_level_that_overflows_the_spacing_is_rejected() {
    let base = layout(1.0);
    assert!(base.at_lod(64).is_ok(), "2^64 is still finite in f64");
    assert!(
        matches!(base.at_lod(2048), Err(crate::Error::InvalidCellSize { .. })),
        "a spacing doubled to infinity must be reported"
    );
}

/// **G-004's acceptance:** LOD 0..3 all mesh cleanly, and the cell count falls
/// by `8^k`.
///
/// The measurement worth having is the one the acceptance does *not* ask for:
/// **triangles fall by `4^k`, not `8^k`.** Cells fill a volume and a surface does
/// not, so each level buys `8x` the sampling work back and only `4x` the
/// rendering. That gap is the whole economics of LOD and it is why the ticket's
/// own figure is about cells.
#[test]
fn every_level_meshes_cleanly_and_the_counts_fall_as_predicted() {
    use crate::fields::ReferenceField;

    let field = Sphere::<f64>::canonical();
    let (lo, _hi) = field.domain();
    // 64 cells of 0.0625 spans the whole 4-unit domain at level 0.
    let base = ChunkLayout::<f64>::new(64, 0.0625, lo).expect("valid layout");

    let mut rows: Vec<(u32, u64, usize, usize)> = Vec::new();
    for level in 0..4u32 {
        // Keep the *world* extent fixed: halve the cell count as the spacing
        // doubles, so every level meshes the same region.
        let cells = 64u32 >> level;
        let lod = ChunkLayout::<f64>::new(cells, base.cell_size(), lo)
            .expect("valid layout")
            .at_lod(level)
            .expect("valid level");
        let shape = lod.sample_shape().expect("valid shape");

        let mut out = MeshBuffer::<f64>::new();
        MarchingCubes::<f64>::new()
            .extract(&field, &shape, lo, lod.cell_size(), &mut out)
            .expect("extraction");

        let report = validate_indexed(
            &out.positions,
            &out.indices,
            &ValidateConfig::from_cell_size(lod.cell_size()).expect("valid cell size"),
        );
        assert!(report.is_closed(), "level {level}:\n{report}");
        assert_eq!(report.euler_characteristic, 2, "level {level}:\n{report}");
        assert_eq!(report.non_manifold_edges, 0, "level {level}:\n{report}");

        let cells_total = u64::from(cells).pow(3);
        rows.push((level, cells_total, out.vertex_count(), out.triangle_count()));
    }

    for (level, cells, vertices, triangles) in &rows {
        std::println!("lod {level}: {cells} cells, {vertices} vertices, {triangles} triangles");
    }

    // Cells: exactly 1/8 per level, by construction.
    for pair in rows.windows(2) {
        assert_eq!(pair[0].1, pair[1].1 * 8, "cells must fall by 8 per level");
    }

    // Triangles: about 1/4 per level, because the surface is two-dimensional.
    //
    // The ratios drift *upward* as the grid stops resolving the sphere --
    // measured 4.114, 4.313, 5.154 -- and that drift is the interesting part.
    // "Triangles scale with area over h squared" is a continuum claim, and by the
    // coarsest level the sphere is four cells across, where a staircase of 104
    // triangles is no longer approximating anything smoothly. So the tight claim
    // is asserted only where its premise holds, and the loose one everywhere.
    let ratios: Vec<f64> = rows
        .windows(2)
        .map(|pair| pair[0].3 as f64 / pair[1].3 as f64)
        .collect();
    for (step, ratio) in ratios.iter().enumerate() {
        std::println!("lod {step} -> {}: triangles fell {ratio:.3}x", step + 1);
    }

    for (step, ratio) in ratios.iter().enumerate() {
        assert!(
            (3.0..6.0).contains(ratio),
            "lod {step} -> {}: {ratio:.3}x is not even roughly a quartering",
            step + 1
        );
    }
    // Where the grid still resolves the surface -- 32 and 16 cells across a unit
    // sphere -- the quartering is tight.
    for (step, ratio) in ratios.iter().take(2).enumerate() {
        assert!(
            (3.8..4.6).contains(ratio),
            "lod {step} -> {}: {ratio:.3}x, expected a tight quartering while the grid resolves",
            step + 1
        );
    }
    // And it is looser at the coarsest step, which is the finding rather than a
    // tolerance. If this ever tightens, the continuum model reaches further than
    // measured and that is worth knowing.
    let coarsest = *ratios.last().expect("three steps");
    assert!(
        coarsest > ratios[0],
        "the quartering did not degrade as the grid stopped resolving: {ratios:?}"
    );
}

/// **What LOD costs, and it is not what the ticket assumed.**
///
/// Written first as `a_feature_thinner_than_the_coarse_spacing_disappears`,
/// asserting the plate is gone by the coarsest level. It is not: measured
/// **4,088 → 1,016 → 248 → 56** triangles across levels 0–3, still 56 at
/// `h = 0.5` where the plate is a fraction of a cell thick.
///
/// The correction is the finding. Marching Cubes samples **corners** and cuts
/// **edges**, so a slab thinner than a cell is not missed — it is *aliased*.
/// Whichever edges happen to straddle it still register a sign change, and what
/// comes back is a partial, holey remnant rather than nothing. A-005 measured the
/// same field returning **zero** triangles under greedy quads, which asks one
/// question per cell *centre* and therefore does miss it cleanly.
///
/// For a streamed world that is the worse of the two behaviours: a feature that
/// vanishes at a known distance can be faded, and one that disintegrates into a
/// resolution-dependent scatter pops. It is also the concrete cost A-014's
/// subgrid work exists to remove — see M-67, where 95.6% of the configurations a
/// tet can be in are invisible to a sign test.
#[test]
fn a_sub_cell_feature_aliases_under_coarsening_rather_than_vanishing() {
    use crate::fields::{ReferenceField, ThinPlate};

    let field = ThinPlate::<f64>::canonical();
    let (lo, hi) = field.domain();
    let base = ChunkLayout::<f64>::new(64, (hi[0] - lo[0]) / 64.0, lo).expect("valid layout");

    let mut counts = Vec::new();
    for level in 0..4u32 {
        let cells = 64u32 >> level;
        let lod = ChunkLayout::<f64>::new(cells, base.cell_size(), lo)
            .expect("valid layout")
            .at_lod(level)
            .expect("valid level");
        let shape = lod.sample_shape().expect("valid shape");
        let mut out = MeshBuffer::<f64>::new();
        MarchingCubes::<f64>::new()
            .extract(&field, &shape, lo, lod.cell_size(), &mut out)
            .expect("extraction");
        std::println!(
            "thin_plate at lod {level} (h = {:.4}): {} triangles",
            lod.cell_size(),
            out.triangle_count()
        );
        counts.push(out.triangle_count());
    }

    assert!(counts[0] > 0, "the plate must exist at full resolution");

    // It thins out sharply -- the feature is genuinely being lost.
    for pair in counts.windows(2) {
        let ratio = pair[0] as f64 / pair[1] as f64;
        assert!(
            ratio > 3.0,
            "the plate is not being lost at all: {counts:?}"
        );
    }

    // And yet it never reaches zero, which is the whole point of this test.
    assert!(
        counts.iter().all(|c| *c > 0),
        "the plate vanished cleanly, which would make this crate's LOD story \
         easier than it is: {counts:?}"
    );
}
