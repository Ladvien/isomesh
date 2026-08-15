//! A-010's acceptance, and the identity the whole construction rests on.

use alloc::vec::Vec;

use super::ManifoldDualContouring;
use crate::cube::{EDGE_CORNERS, EDGE_COUNT, corner_inside};
use crate::dual::MAX_CELL_VERTICES;
use crate::dual_contouring::DualContouring;
use crate::fields::{ReferenceField, Sphere, Torus, capped_gyroid, csg_difference};
use crate::marching_cubes::MarchingCubes;
use crate::marching_cubes::table::{AMBIGUOUS_FACES, NO_EDGE, segment_links};
use crate::validate::{ValidateConfig, check_determinism, self_intersections, validate_indexed};
use crate::{MeshBuffer, RuntimeShape3, Sdf};

/// Mesh a reference field at `samples` per axis, returning the mesh and the cell
/// size — the same convention the rest of the suite uses.
fn mesh_mdc<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    samples: u32,
) -> (MeshBuffer<f64>, f64) {
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let mut out = MeshBuffer::<f64>::new();
    ManifoldDualContouring::<f64>::new()
        .extract(field, &shape, lo, h, &mut out)
        .expect("extraction");
    (out, h)
}

/// The same field through plain Dual Contouring, for the comparisons.
fn mesh_dc<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    samples: u32,
) -> (MeshBuffer<f64>, f64) {
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let mut out = MeshBuffer::<f64>::new();
    DualContouring::<f64>::new()
        .extract(field, &shape, lo, h, &mut out)
        .expect("extraction");
    (out, h)
}

fn report_of(mesh: &MeshBuffer<f64>, h: f64) -> crate::validate::MeshReport {
    validate_indexed(
        &mesh.positions,
        &mesh.indices,
        &ValidateConfig::from_cell_size(h).expect("valid cell size"),
    )
}

// ─── the table properties the construction rests on ─────────────────────────

/// How many cycles a configuration produces, by the same walk the rule uses.
fn cycle_count(case: u8, joined: u8) -> usize {
    let next = segment_links(case, joined);
    let mut visited = 0u16;
    let mut cycles = 0;
    for start in 0..EDGE_COUNT as u8 {
        if next[start as usize] == NO_EDGE || visited & (1 << start) != 0 {
            continue;
        }
        cycles += 1;
        let mut current = start;
        while visited & (1 << current) == 0 {
            visited |= 1 << current;
            current = next[current as usize];
        }
    }
    cycles
}

/// Every canonical `(case, joined)` pair — the ones a real cell can present.
fn every_configuration() -> impl Iterator<Item = (u8, u8)> {
    (0..=255u8).flat_map(|case| {
        (0..64u8)
            .filter(move |joined| joined & !AMBIGUOUS_FACES[case as usize] == 0)
            .map(move |joined| (case, joined))
    })
}

/// [`MAX_CELL_VERTICES`] is the array size the rule writes into, so an
/// underestimate is an out-of-bounds write. Checked over every case and every
/// face-resolution mask, not just the ones some field happens to produce.
#[test]
fn every_case_fits_the_slot_budget() {
    let mut worst = 0;
    let mut worst_case = 0u8;
    for (case, joined) in every_configuration() {
        let n = cycle_count(case, joined);
        if n > worst {
            worst = n;
            worst_case = case;
        }
    }
    std::println!("worst cycle count {worst}, first at case {worst_case:#010b}");
    assert!(
        worst <= MAX_CELL_VERTICES,
        "{worst} cycles needed, budget is {MAX_CELL_VERTICES}"
    );
    // Four is not an accident: the four corners of one tetrahedron, each
    // isolated from the others. If this ever drops, the budget is over-sized,
    // which is worth knowing too.
    assert_eq!(worst, 4);
}

/// Every cut edge belongs to **exactly one** cycle.
///
/// This is Schaefer, Ju & Warren's *"each edge is associated with exactly one
/// vertex"* checked against this crate's own table, and it is what makes the
/// quad walk's `(cell, edge) → vertex` lookup total. If it failed, some crossed
/// grid edge would have a corner with no vertex to use.
#[test]
fn every_cut_edge_belongs_to_exactly_one_cycle() {
    for (case, joined) in every_configuration() {
        let next = segment_links(case, joined);

        let mut cut = 0u16;
        for (edge, [lo, hi]) in EDGE_CORNERS.into_iter().enumerate() {
            if corner_inside(case, lo) != corner_inside(case, hi) {
                cut |= 1 << edge;
            }
        }

        let mut covered = 0u16;
        for start in 0..EDGE_COUNT as u8 {
            if next[start as usize] == NO_EDGE || covered & (1 << start) != 0 {
                continue;
            }
            let mut current = start;
            while covered & (1 << current) == 0 {
                covered |= 1 << current;
                current = next[current as usize];
            }
        }

        assert_eq!(
            covered, cut,
            "case {case:#010b} joined {joined:#08b}: covered {covered:#014b}, cut {cut:#014b}"
        );
    }
}

// ─── A-010's acceptance criterion ───────────────────────────────────────────

/// The ticket, verbatim: `non_manifold_edges == 0` on `gyroid` and
/// `csg_difference`, where plain Dual Contouring will not manage it.
///
/// Asserted **in both directions**. A test that only checked the new algorithm
/// would still pass if those fields had stopped being hard, and then it would be
/// measuring nothing.
#[test]
fn the_fields_dual_contouring_pinches_come_out_manifold() {
    let mut any_pinched = false;

    for samples in [33u32, 49] {
        let field = capped_gyroid::<f64>();
        let (mdc, h) = mesh_mdc(&field, samples);
        let (dc, _) = mesh_dc(&field, samples);
        let mdc_report = report_of(&mdc, h);
        let dc_report = report_of(&dc, h);

        std::println!(
            "gyroid {samples}^3: dual contouring {} non-manifold edges / {} vertices, \
             manifold dual contouring {} / {}",
            dc_report.non_manifold_edges,
            dc_report.non_manifold_vertices,
            mdc_report.non_manifold_edges,
            mdc_report.non_manifold_vertices,
        );
        assert_eq!(
            mdc_report.non_manifold_edges, 0,
            "gyroid at {samples}^3:\n{mdc_report}"
        );
        assert_eq!(
            mdc_report.non_manifold_vertices, 0,
            "gyroid at {samples}^3:\n{mdc_report}"
        );
        any_pinched |= dc_report.non_manifold_edges > 0;
    }

    let field = csg_difference::<f64>();
    for samples in [33u32, 49] {
        let (mdc, h) = mesh_mdc(&field, samples);
        let (dc, _) = mesh_dc(&field, samples);
        let mdc_report = report_of(&mdc, h);
        let dc_report = report_of(&dc, h);
        std::println!(
            "csg_difference {samples}^3: dual contouring {} non-manifold edges, \
             manifold dual contouring {}",
            dc_report.non_manifold_edges,
            mdc_report.non_manifold_edges,
        );
        assert_eq!(
            mdc_report.non_manifold_edges, 0,
            "csg_difference at {samples}^3:\n{mdc_report}"
        );
    }

    assert!(
        any_pinched,
        "dual contouring no longer pinches on the gyroid — the comparison is vacuous \
         and A-010's premise needs re-checking"
    );
}

/// Every reference field, manifold, by the gate the field itself publishes.
///
/// This is the gate Surface Nets and Dual Contouring cannot pass, and the reason
/// `SurfaceGate::ClosedAllowingUnresolvedTopology` exists. A-010 is what makes it
/// assertable again.
///
/// # A-010's zero was conditional, and A-002e is what found the condition (M-211)
///
/// It held on seven fields for one reason nobody had checked: **not one of them
/// produces a cell with an interior ambiguity** — 0 of 68,385 surface cells across
/// all seven at three resolutions (M-208). `noise_cavity` was added precisely to
/// reach that configuration, and Manifold Dual Contouring stops being manifold on
/// it. Whether that is a defect in this crate's construction or a limit of the
/// published guarantee is **A-017's**, not this test's.
///
/// So the count is pinned rather than asserted zero, which is M-4's precedent for
/// Surface Nets and Dual Contouring: *a known defect with a pinned number and a
/// ticket that owns it satisfies this gate; an unexplained one does not.* The
/// census is compared **whole**, so it fails if the defect spreads to another
/// field *and* if it silently disappears — and no gate is selected by a field's
/// name, which `CLAUDE.md` forbids outright.
#[test]
fn every_reference_field_meshes_manifold() {
    let mut checked = 0;
    let mut census: Vec<(&str, u32, u64, u64)> = Vec::new();
    crate::for_each_reference_field!(f64, |name, field| {
        for samples in [17u32, 33] {
            let (mesh, h) = mesh_mdc(&field, samples);
            if mesh.triangle_count() == 0 {
                continue;
            }
            let report = report_of(&mesh, h);
            if report.non_manifold_edges != 0 || report.non_manifold_vertices != 0 {
                census.push((
                    name,
                    samples,
                    report.non_manifold_edges,
                    report.non_manifold_vertices,
                ));
            }
            // Orientation is *not* relaxed. A non-manifold edge is a count; an
            // inside-out triangle is a different failure and nothing here excuses
            // one.
            assert_eq!(
                report.inconsistently_oriented_edges, 0,
                "{name} at {samples}^3:\n{report}"
            );
            if let Some(chi) = field.expected_euler() {
                assert_eq!(
                    report.euler_characteristic, chi,
                    "{name} at {samples}^3:\n{report}"
                );
            }
            checked += 1;
        }
    });
    assert_eq!(
        census, MDC_NON_MANIFOLD_CENSUS,
        "the manifold census moved — see A-017, and do not re-bless this without \
         reading which field changed"
    );
    assert!(
        checked >= 12,
        "only {checked} field/resolution pairs meshed"
    );
}

/// Where Manifold Dual Contouring is *not* manifold: `(field, samples, edges,
/// vertices)`, and every row is owned by **A-017**.
///
/// Empty for the seven original reference fields at every resolution. Non-empty
/// only where a cell carries an interior ambiguity, which before A-002e no field
/// in this crate reached.
/// Note the exact factor of two in both rows: **every offending edge carries
/// exactly two offending vertices, its own endpoints**, so the defect is a
/// property of edges and the vertex count is a corollary rather than a second
/// mechanism. Recorded here because it is the kind of identity ✗1 says to assert
/// and let the counterexample explain itself.
const MDC_NON_MANIFOLD_CENSUS: &[(&str, u32, u64, u64)] =
    &[("noise_cavity", 17, 30, 60), ("noise_cavity", 33, 64, 128)];

/// **P-5, pre-registered before running:** the output is the dual of Marching
/// Cubes, and the dual of a surface has `V' = F`, `E' = E`, `F' = V`, so the
/// Euler characteristic carries across *unchanged*.
///
/// Far stronger than a manifold count: it says the algorithm reproduces Marching
/// Cubes' topology field by field, which is the entire reason to expect it to be
/// manifold at all. Only the closed fields are compared — on an open one the two
/// methods trim different amounts at the grid border, so χ is not shared and the
/// claim does not apply.
///
/// # Where it stops holding, and why that is the same finding as above
///
/// The prediction assumed the dual is *a* dual — one vertex per Marching Cubes
/// face, one face per vertex. On a cell whose interior is ambiguous that
/// correspondence breaks before any code runs: the two methods are not describing
/// the same surface there, because Marching Cubes separates a tunnel that the
/// trilinear interpolant joins. `noise_cavity` is the first field in this crate
/// to contain such a cell (M-208), and it is exactly where χ parts company.
/// Pinned whole, and owned by **A-017** — the component count is *not* relaxed
/// with it, because fusing pieces is the defect A-010 exists to have fixed.
#[test]
fn euler_characteristic_matches_marching_cubes_on_closed_fields() {
    let mut checked = 0;
    let mut chi_census: Vec<(&str, u32, i64, i64)> = Vec::new();
    crate::for_each_reference_field!(f64, |name, field| {
        if field.closed_in_domain() {
            for samples in [17u32, 25, 33] {
                let (lo, hi) = field.domain();
                let h = (hi[0] - lo[0]) / f64::from(samples - 1);
                let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

                let mut mc_out = MeshBuffer::<f64>::new();
                MarchingCubes::<f64>::new()
                    .extract(&field, &shape, lo, h, &mut mc_out)
                    .expect("extraction");
                if mc_out.triangle_count() == 0 {
                    continue;
                }

                let (mdc, _) = mesh_mdc(&field, samples);
                let mc_chi = report_of(&mc_out, h).euler_characteristic;
                let mdc_chi = report_of(&mdc, h).euler_characteristic;

                if mdc_chi != mc_chi {
                    chi_census.push((name, samples, mc_chi, mdc_chi));
                }

                // Components too, and this is the sharper half. Where plain Dual
                // Contouring pinches, the shared vertex **fuses pieces that
                // Marching Cubes keeps apart** — the gyroid at 19^3 reads one
                // component under Dual Contouring and seven under both Marching
                // Cubes and this. So the pinch was not only a topological defect
                // in the index buffer, it was reporting the wrong object.
                let mc_parts = report_of(&mc_out, h).components;
                let mdc_parts = report_of(&mdc, h).components;
                assert_eq!(
                    mdc_parts, mc_parts,
                    "{name} at {samples}^3: marching cubes {mc_parts} components, dual {mdc_parts}"
                );
                checked += 1;
            }
        }
    });
    assert_eq!(
        chi_census, MDC_CHI_CENSUS,
        "the chi census moved — see A-017, and read which field changed before \
         re-blessing it"
    );
    assert!(checked >= 12, "only {checked} comparisons ran");
}

/// Where Manifold Dual Contouring's Euler characteristic parts company with
/// Marching Cubes': `(field, samples, marching cubes, dual)`. Owned by **A-017**.
///
/// Empty for the seven original fields. The dual identity `V' = F, E' = E,
/// F' = V` needs the two methods to be describing the same surface, and on a cell
/// with an interior ambiguity they are not.
const MDC_CHI_CENSUS: &[(&str, u32, i64, i64)] = &[
    ("noise_cavity", 17, -30, 0),
    ("noise_cavity", 25, -78, -1),
    ("noise_cavity", 33, -96, -32),
];

// ─── the standard per-algorithm gate ────────────────────────────────────────

#[test]
fn a_meshed_sphere_is_closed() {
    let (mesh, h) = mesh_mdc(&Sphere::<f64>::canonical(), 33);
    let report = report_of(&mesh, h);
    assert!(report.is_closed(), "{report}");
    assert_eq!(report.euler_characteristic, 2);
    assert_eq!(report.non_manifold_edges, 0);
    assert_eq!(report.boundary_edges, 0);
    assert_eq!(report.inconsistently_oriented_edges, 0);
}

#[test]
fn a_meshed_torus_has_genus_one() {
    let (mesh, h) = mesh_mdc(&Torus::<f64>::canonical(), 33);
    let report = report_of(&mesh, h);
    assert!(report.is_closed(), "{report}");
    assert_eq!(report.genus, Some(1), "{report}");
}

/// No manifold or Euler check can see a globally inverted surface; only the
/// signed volume can.
#[test]
fn meshed_sphere_has_positive_signed_volume() {
    let (mesh, _) = mesh_mdc(&Sphere::<f64>::canonical(), 33);
    let mut volume = 0.0f64;
    for tri in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[tri[0] as usize];
        let b = mesh.positions[tri[1] as usize];
        let c = mesh.positions[tri[2] as usize];
        volume += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
    assert!(
        volume > 0.0,
        "signed volume {volume} — the mesh is inside out"
    );
}

#[test]
fn extraction_is_deterministic() {
    let field = capped_gyroid::<f64>();
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / 24.0;
    let shape = RuntimeShape3::new([25; 3]).expect("valid shape");

    let report = check_determinism(|out| {
        ManifoldDualContouring::<f64>::new()
            .extract(&field, &shape, lo, h, out)
            .expect("extraction");
    });
    assert!(report.is_deterministic(), "{report}");
}

// ─── what the split costs, and where ────────────────────────────────────────

/// A sphere passes through no cell twice, so every cell has exactly one cycle
/// and the output must be **identical** to Dual Contouring's, bit for bit.
///
/// This is the test that says splitting costs nothing where there is nothing to
/// split — and it is also the strongest possible check that the per-cycle solve
/// degenerates to the whole-cell solve when the cycle *is* the whole cell.
#[test]
fn a_single_sheet_field_reproduces_dual_contouring_exactly() {
    for samples in [17u32, 33] {
        let field = Sphere::<f64>::canonical();
        let (mdc, _) = mesh_mdc(&field, samples);
        let (dc, _) = mesh_dc(&field, samples);
        assert_eq!(
            mdc.indices, dc.indices,
            "connectivity differs at {samples}^3"
        );
        assert_eq!(
            mdc.positions, dc.positions,
            "vertex positions differ at {samples}^3"
        );
    }
}

/// Where the extra vertices are, and how many.
///
/// Nielson reports multi-vertex configurations at *"about 1.3% of all
/// configurations"*; this is the same quantity measured on real fields instead
/// of over the case table, which is the number that actually decides the memory
/// cost.
#[test]
fn the_split_census_is_reported() {
    let mut rows: Vec<(&str, u32, usize, usize)> = Vec::new();
    crate::for_each_reference_field!(f64, |name, field| {
        for samples in [17u32, 25, 33] {
            let (mdc, _) = mesh_mdc(&field, samples);
            let (dc, _) = mesh_dc(&field, samples);
            rows.push((name, samples, mdc.vertex_count(), dc.vertex_count()));
        }
    });

    for (name, samples, mdc, dc) in &rows {
        let extra = mdc - dc;
        let pct = if *dc == 0 {
            0.0
        } else {
            100.0 * extra as f64 / *dc as f64
        };
        std::println!(
            "{name} {samples}^3: {mdc} vertices against {dc} — {extra} split ({pct:.2}%)"
        );
    }

    let split = |want: &str| -> bool { rows.iter().filter(|r| r.0 == want).any(|r| r.2 > r.3) };
    assert!(split("gyroid"), "no cell split on the gyroid");
    assert!(
        !split("sphere"),
        "a sphere has no multi-sheet cell and must never split"
    );
}

/// **The one thing vertex splitting does not fix, isolated and pinned.**
///
/// Splitting per cycle removes the shared-vertex pinch entirely. What remains is
/// a different mechanism, and this is the case that exhibits it — the *same*
/// three-sphere union at `h = 2/3` that falsified unconditional manifoldness for
/// Marching Cubes at ✗15, found again independently by
/// `manifold_dual_contouring_meshes_sphere_unions` shrinking to it.
///
/// The dual of a manifold surface is a manifold **complex**. It is not
/// necessarily a manifold *indexed mesh*, because an index buffer cannot carry
/// two distinct edges between the same pair of vertices. Where two cells share a
/// face that carries **two** surface segments, and each cell puts both segments
/// in the same cycle, the two dual edges have the same two endpoints and collapse
/// into one edge with four faces.
///
/// The arithmetic identifies the mechanism rather than merely observing the
/// symptom: each collapse removes exactly one edge from `E` and nothing else, so
///
/// ```text
/// chi_dual - chi_marching_cubes == non_manifold_edges
/// ```
///
/// must hold exactly. Asserted at `n = 7`, where both sides are 1, and across
/// every refinement, where both sides are 0.
#[test]
fn the_parallel_dual_edge_collapse_is_the_only_residue() {
    use crate::fields::Sphere;
    use crate::property::SphereUnion;

    // Verbatim the ✗15 fixture, so the two findings are demonstrably about the
    // same configuration and not merely about similar ones.
    let field = SphereUnion {
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

    // The property suite's own grid, so this is the case it found, exactly.
    let reports_at = |n: u32| {
        let half = crate::property::DOMAIN;
        let h = 2.0 * half / f64::from(n - 1);
        let shape = RuntimeShape3::new([n; 3]).expect("valid shape");
        let origin = [-half; 3];

        let mut mc_out = MeshBuffer::<f64>::new();
        MarchingCubes::<f64>::new()
            .extract(&field, &shape, origin, h, &mut mc_out)
            .expect("extraction");

        let mut mdc_out = MeshBuffer::<f64>::new();
        ManifoldDualContouring::<f64>::new()
            .extract(&field, &shape, origin, h, &mut mdc_out)
            .expect("extraction");

        (report_of(&mc_out, h), report_of(&mdc_out, h))
    };

    let (mc, mdc) = reports_at(7);
    std::println!(
        "n=7: marching cubes chi {} manifold-edges {}; dual chi {} manifold-edges {} vertices {}",
        mc.euler_characteristic,
        mc.non_manifold_edges,
        mdc.euler_characteristic,
        mdc.non_manifold_edges,
        mdc.non_manifold_vertices,
    );

    // Marching Cubes is clean here — A-015 saw to that — so the dual's defect is
    // the dual's own, not one inherited from the surface being dualised.
    assert_eq!(mc.non_manifold_edges, 0, "{mc}");
    assert_eq!(mc.euler_characteristic, 0, "{mc}");

    // Exactly one collapse. Pinned in both directions, per M-4's precedent: this
    // fails if the defect spreads *and* if it silently disappears.
    assert_eq!(mdc.non_manifold_edges, 1, "{mdc}");
    assert_eq!(mdc.boundary_edges, 0, "{mdc}");
    assert_eq!(mdc.inconsistently_oriented_edges, 0, "{mdc}");

    // The mechanism, as arithmetic.
    assert_eq!(
        mdc.euler_characteristic - mc.euler_characteristic,
        mdc.non_manifold_edges as i64,
        "a collapse must cost exactly one edge\nmc:\n{mc}\ndual:\n{mdc}"
    );

    // And it is a coarse-grid effect, exactly as ✗15's was: refine and it goes.
    for n in [9u32, 13, 17, 25, 33, 49] {
        let (mc, mdc) = reports_at(n);
        assert_eq!(mdc.non_manifold_edges, 0, "n={n}\n{mdc}");
        assert_eq!(mdc.non_manifold_vertices, 0, "n={n}\n{mdc}");
        assert_eq!(
            mdc.euler_characteristic, mc.euler_characteristic,
            "n={n}: the dual must carry chi across unchanged\nmc:\n{mc}\ndual:\n{mdc}"
        );
    }
}

/// Self-intersection is what this does **not** fix. ✗2 records ODC (2024)
/// measuring Manifold Dual Contouring at 100% of models self-intersecting, and
/// Manson & Schaefer's within-cell partition argument — which the clamp rests on
/// — assumes one vertex per cell.
///
/// Recorded, not gated, exactly as Dual Contouring's is.
#[test]
fn the_self_intersection_count_is_recorded() {
    crate::for_each_reference_field!(f64, |name, field| {
        let samples = 33u32;
        let (mdc, h) = mesh_mdc(&field, samples);
        let (dc, _) = mesh_dc(&field, samples);
        if mdc.triangle_count() > 0 {
            let a = self_intersections(&mdc.positions, &mdc.indices, h).expect("valid cell size");
            let b = self_intersections(&dc.positions, &dc.indices, h).expect("valid cell size");
            std::println!(
                "{name} {samples}^3: manifold dual contouring {:.3} per 1k, dual contouring {:.3}",
                a.per_thousand_triangles(),
                b.per_thousand_triangles(),
            );
        }
    });
}
