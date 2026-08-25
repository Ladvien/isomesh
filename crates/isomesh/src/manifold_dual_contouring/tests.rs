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
    for tri in mesh.indices.as_chunks::<3>().0 {
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

/// **A-017's characterisation, and it rules out the two explanations the ticket
/// offered.**
///
/// Manifold Dual Contouring is not manifold on `noise_cavity`. The ticket asked
/// whether that is a defect in this crate's construction or the published
/// guarantee not covering a cell whose interior the trilinear interpolant joins.
/// Measured, it is **neither, because it is not about tunnels at all** (M-224):
///
/// - **Only one** offending edge in each of 30 and 64 lies within `1.5h` of a
///   tunnel cell, while **all** of them lie within `1.5h` of an *ambiguous* cell.
///   The field has 193 and 502 ambiguous cells against 3 and 2 tunnels.
/// - **Every** offending edge carries exactly **four faces, all four distinct** —
///   no threes, no fives, and no duplication. Four distinct triangles on one edge
///   is two sheets meeting along it, which is a genuine junction rather than the
///   double-emission ✗17 found behind Marching Cubes' fan chords.
/// - It survives the **correct** face rule. Under `AsymptoticDecider` the count
///   falls from 30 to 8 at 17³ but does not reach zero — and the same setting
///   *introduces* 3 offending edges on `gyroid` at 25³, where `Separate` gives
///   none.
///
/// So the remaining hypothesis is about the quad walk around a crossed grid edge,
/// not about interior topology and not about face pairing. Pinned here in both
/// directions so that whatever A-017 eventually does has to move these numbers.
#[test]
fn the_manifold_dual_contouring_defect_is_four_distinct_faces_on_one_edge() {
    use crate::validate::validate_features;
    use alloc::collections::{BTreeMap, BTreeSet};

    let field = crate::fields::noise_cavity::<f64>();
    let mut rows = Vec::new();
    for samples in [17u32, 33] {
        let (mesh, h) = mesh_mdc(&field, samples);
        let (_report, features) = validate_features(
            &mesh.positions,
            &mesh.indices,
            &ValidateConfig::from_cell_size(h).expect("valid cell size"),
        );
        let offending: BTreeSet<(u32, u32)> = features.edges.iter().map(|e| (e[0], e[1])).collect();

        let mut faces: BTreeMap<(u32, u32), usize> = BTreeMap::new();
        let mut distinct: BTreeMap<(u32, u32), BTreeSet<[u32; 3]>> = BTreeMap::new();
        for t in mesh.indices.as_chunks::<3>().0 {
            let mut key = [t[0], t[1], t[2]];
            key.sort_unstable();
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                let e = (a.min(b), a.max(b));
                if offending.contains(&e) {
                    *faces.entry(e).or_insert(0) += 1;
                    distinct.entry(e).or_default().insert(key);
                }
            }
        }

        for (e, n) in &faces {
            assert_eq!(*n, 4, "edge {e:?} carries {n} faces, not four");
            assert_eq!(
                distinct[e].len(),
                4,
                "edge {e:?} carries {n} faces but only {} distinct triangles — \
                 that would be duplication, which is a different defect",
                distinct[e].len()
            );
        }
        rows.push((samples, offending.len()));
    }
    assert_eq!(
        rows,
        alloc::vec![(17, 30), (33, 64)],
        "the A-017 census moved"
    );
    std::println!("measured: A-017 offending edges {rows:?}, all four distinct faces");
}

/// **A-017's mechanism, and it predicts the count exactly under both face rules.**
///
/// An ambiguous face has **all four** of its edges cut. Manifold Dual Contouring
/// puts one vertex per cycle per cell, so if all four of those edges belong to
/// *one* cycle in each of the two cells sharing the face, then all four dual
/// quads — one per crossed grid edge — connect the **same pair** of cell
/// vertices. Four quads on one dual edge, and a quad contributes exactly one of
/// its two triangles to each of its sides, which is precisely the *four distinct
/// faces* measured above.
///
/// So the count is not an estimate. It is
///
/// ```text
/// non_manifold_edges  ==  shared ambiguous faces whose four cut edges
///                         lie in one cycle on both sides
/// ```
///
/// and it holds for `Separate` (30 and 64) and for `AsymptoticDecider` (8 and 40)
/// alike — computed here from the *grid*, with no mesh involved, against a count
/// the validator takes from the *mesh*, with no grid involved (M-225).
///
/// **This is a limit of one-vertex-per-cycle, not a bug in this transcription.**
/// Schaefer, Ju & Warren's argument separates sheets *within* a cell, and nothing
/// in it stops two different crossed edges of one shared face resolving to the
/// same pair of cycles. A-017 owns what to do about it; this test owns knowing
/// exactly what it is.
#[test]
fn the_defect_count_is_predicted_from_the_grid_alone() {
    use crate::cube::{corner_offset, edge_on_face, is_inside};
    use crate::marching_cubes::FaceAmbiguity;
    use crate::marching_cubes::ambiguity::joined_mask;
    use crate::marching_cubes::table::AMBIGUOUS_FACES;
    use crate::marching_cubes::trilinear::Contours;

    let field = crate::fields::noise_cavity::<f64>();
    let (lo, _hi) = field.domain();
    let mut rows = Vec::new();

    for rule in [FaceAmbiguity::Separate, FaceAmbiguity::AsymptoticDecider] {
        for samples in [17u32, 33] {
            let h = 4.0 / f64::from(samples - 1);
            let cell_of = |x: u32, y: u32, z: u32| {
                let mut corner = [0.0f64; 8];
                let mut case = 0u8;
                for (c, slot) in corner.iter_mut().enumerate() {
                    let o = corner_offset(c as u8);
                    *slot = field.sample([
                        lo[0] + h * f64::from(x + o[0]),
                        lo[1] + h * f64::from(y + o[1]),
                        lo[2] + h * f64::from(z + o[2]),
                    ]);
                    if is_inside(*slot) {
                        case |= 1 << c;
                    }
                }
                (case, corner)
            };
            let rings_of = |case: u8, corner: &[f64; 8]| {
                let mask = match rule {
                    FaceAmbiguity::Separate => 0,
                    FaceAmbiguity::AsymptoticDecider => {
                        joined_mask(corner, AMBIGUOUS_FACES[case as usize])
                    }
                };
                let contours = Contours::of(case, mask);
                let mut owner = [255u8; 12];
                for r in 0..contours.count() {
                    for &e in contours.ring(r) {
                        owner[e as usize] = r as u8;
                    }
                }
                owner
            };

            let mut predicted = 0usize;
            for axis in 0..3usize {
                for z in 0..samples - 1 {
                    for y in 0..samples - 1 {
                        for x in 0..samples - 1 {
                            let mut n = [x, y, z];
                            n[axis] += 1;
                            if n[axis] >= samples - 1 {
                                continue;
                            }
                            let (ca, va) = cell_of(x, y, z);
                            if AMBIGUOUS_FACES[ca as usize] & (1u8 << (axis * 2 + 1)) == 0 {
                                continue;
                            }
                            let (cb, vb) = cell_of(n[0], n[1], n[2]);
                            let oa = rings_of(ca, &va);
                            let ob = rings_of(cb, &vb);
                            let cut = |owner: &[u8; 12], side: u8| -> Vec<u8> {
                                (0..12u8)
                                    .filter(|&e| {
                                        edge_on_face(e, axis, side) && owner[e as usize] != 255
                                    })
                                    .collect()
                            };
                            let (cut_a, cut_b) = (cut(&oa, 1), cut(&ob, 0));
                            if cut_a.len() != 4 || cut_b.len() != 4 {
                                continue;
                            }
                            let one = |owner: &[u8; 12], edges: &[u8]| {
                                edges
                                    .iter()
                                    .all(|&e| owner[e as usize] == owner[edges[0] as usize])
                            };
                            if one(&oa, &cut_a) && one(&ob, &cut_b) {
                                predicted += 1;
                            }
                        }
                    }
                }
            }

            let (lo_d, hi_d) = field.domain();
            let cell = (hi_d[0] - lo_d[0]) / f64::from(samples - 1);
            let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
            let mut mdc = ManifoldDualContouring::<f64>::new();
            mdc.set_face_ambiguity(rule);
            let mut out = MeshBuffer::<f64>::new();
            mdc.extract(&field, &shape, lo_d, cell, &mut out)
                .expect("extraction");
            let observed = report_of(&out, cell).non_manifold_edges as usize;

            assert_eq!(
                predicted, observed,
                "{rule:?} at {samples}^3: the grid predicts {predicted} and the mesh has {observed}"
            );
            rows.push((rule, samples, predicted));
        }
    }
    std::println!("measured: A-017 predicted == observed for {rows:?}");
}

// ─── A-025: the defect, constructed rather than sampled ─────────────────────

const NX: usize = 4;
const NY: usize = 4;
const NZ: usize = 3;

/// A field defined by an explicit lattice of sample values, trilinear between
/// them.
///
/// Trilinear is not an arbitrary choice of interpolant: it is the one the
/// asymptotic decider itself assumes, so the constructed fixture and the rule
/// under test agree about what lies between the samples.
struct Lattice {
    /// `[x][y][z]`, sampled at integer world coordinates with `h = 1`.
    values: [[[f64; NZ]; NY]; NX],
}

impl Lattice {
    /// Base index and fraction on one axis. The index is clamped so the cell is
    /// always a real one; the fraction is **not**, so just outside the lattice
    /// the field extends linearly rather than flattening — a flat extension
    /// would hand the gradient a zero everywhere past the boundary.
    fn split(p: f64, n: usize) -> (usize, f64) {
        let floor = libm::floor(p);
        let i = if floor < 0.0 {
            0
        } else if floor as usize > n - 2 {
            n - 2
        } else {
            floor as usize
        };
        (i, p - i as f64)
    }

    fn corners(&self, p: [f64; 3]) -> ([f64; 8], [f64; 3]) {
        let (i, u) = Self::split(p[0], NX);
        let (j, v) = Self::split(p[1], NY);
        let (k, w) = Self::split(p[2], NZ);
        let mut c = [0.0; 8];
        for (n, slot) in c.iter_mut().enumerate() {
            *slot = self.values[i + (n & 1)][j + ((n >> 1) & 1)][k + ((n >> 2) & 1)];
        }
        (c, [u, v, w])
    }
}

impl Sdf for Lattice {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let (c, t) = self.corners(p);
        let mut acc = 0.0;
        for (n, value) in c.iter().enumerate() {
            let mut weight = 1.0;
            for (axis, &frac) in t.iter().enumerate() {
                weight *= if (n >> axis) & 1 == 0 {
                    1.0 - frac
                } else {
                    frac
                };
            }
            acc += value * weight;
        }
        acc
    }

    /// Analytic, because the default central difference straddles a cell
    /// boundary wherever the crossing sits on a grid edge — which is every
    /// crossing there is — and trilinear is only C⁰ there.
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let (c, t) = self.corners(p);
        let mut g = [0.0; 3];
        for (axis, slot) in g.iter_mut().enumerate() {
            for (n, value) in c.iter().enumerate() {
                let mut factor = 1.0;
                for (other, &frac) in t.iter().enumerate() {
                    let bit = (n >> other) & 1;
                    factor *= if other == axis {
                        if bit == 0 { -1.0 } else { 1.0 }
                    } else if bit == 0 {
                        1.0 - frac
                    } else {
                        frac
                    };
                }
                *slot += value * factor;
            }
        }
        g
    }
}

/// The 3×3×2-cell block, with the two-cell column at `(1, 1, ·)` carrying
/// `pattern` — bit `z·4 + (y≪1) + x`, set meaning inside — and everything
/// around it outside.
///
/// **Why 3×3×2 and not 1×1×2.** The defect is four dual quads landing on one
/// dual edge, and a quad exists only where all four cells around its grid edge
/// do. The shared face's four grid edges reach out one cell in `x` and one in
/// `y` on either side, so the two cells alone would produce no quads at all and
/// the fixture would measure nothing. This is the smallest block in which every
/// quad the mechanism needs can form.
fn column_lattice(pattern: u16) -> Lattice {
    let mut values = [[[1.0f64; NZ]; NY]; NX];
    for low in 0..4usize {
        let column = &mut values[1 + (low & 1)][1 + ((low >> 1) & 1)];
        for (z, slot) in column.iter_mut().enumerate() {
            if pattern & (1 << (z * 4 + low)) != 0 {
                *slot = -1.0;
            }
        }
    }
    Lattice { values }
}

/// The case index a cell of the column carries, in the crate's corner order.
fn column_case(pattern: u16, plane: u32) -> u8 {
    let mut case = 0u8;
    for c in 0..8u8 {
        let o = crate::cube::corner_offset(c);
        let low = (o[1] << 1) | o[0];
        if pattern & (1 << ((plane + o[2]) * 4 + low)) != 0 {
            case |= 1 << c;
        }
    }
    case
}

/// Which cycle owns each cube edge, by the walk the rule uses.
fn cycle_owners(case: u8, joined: u8) -> [u8; 12] {
    let contours = crate::marching_cubes::trilinear::Contours::of(case, joined);
    let mut owner = [255u8; 12];
    for r in 0..contours.count() {
        for &e in contours.ring(r) {
            owner[e as usize] = r as u8;
        }
    }
    owner
}

/// Every non-manifold edge of a mesh, with the distinct triangles sitting on it.
fn offending_edges(
    mesh: &MeshBuffer<f64>,
) -> alloc::collections::BTreeMap<(u32, u32), alloc::collections::BTreeSet<[u32; 3]>> {
    use alloc::collections::{BTreeMap, BTreeSet};
    let (_report, features) = crate::validate::validate_features(
        &mesh.positions,
        &mesh.indices,
        &ValidateConfig::from_cell_size(1.0).expect("valid cell size"),
    );
    let bad: BTreeSet<(u32, u32)> = features.edges.iter().map(|e| (e[0], e[1])).collect();
    let mut distinct: BTreeMap<(u32, u32), BTreeSet<[u32; 3]>> = BTreeMap::new();
    for t in mesh.indices.as_chunks::<3>().0 {
        let mut key = [t[0], t[1], t[2]];
        key.sort_unstable();
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            let e = (a.min(b), a.max(b));
            if bad.contains(&e) {
                distinct.entry(e).or_default().insert(key);
            }
        }
    }
    distinct
}

/// One of M-292's eighteen. Bit `z·4 + (y≪1) + x` of the column, set = inside.
const OFFENDER: u16 = 0b0111_0110_0111;

/// **A-025 — one hand-built configuration that makes Manifold Dual Contouring
/// non-manifold, with Marching Cubes manifold on the very same samples.**
///
/// M-290 measured the defect over eight fields and M-292 bounded it over all
/// 4,096 two-cell sign patterns; both are censuses. What neither gives is a
/// fixture small enough to read, and A-021's lesson is that a census names the
/// *rate* while a constructed case names the *mechanism*. This is the
/// constructed case: **48 samples**, no field, no noise.
///
/// # The configuration
///
/// Twelve of those samples, `−1` inside and `+1` outside, in a `2×2×3` column
/// embedded in a `4×4×3` lattice whose every other sample is outside:
///
/// ```text
///      z = 0          z = 1          z = 2
///    y↑ in  out     y↑ in  out     y↑ in  out
///       in  in         out in         in  in
///          → x
/// ```
///
/// The middle plane is the whole of it: the shared face's four corners alternate
/// `out, in, in, out`, which is the face saddle — all four of its edges cut, and
/// the two diagonals equally entitled to be joined. Of the 4,096 patterns this
/// is one of the **18** that offend under the mask `FaceAmbiguity::Separate`
/// produces (M-292), and all 18 are, exactly, the ones where both cells resolve
/// to a **single** cycle.
///
/// # What it demonstrates
///
/// One cycle per cell is one vertex per cell, so all four quads — one for each
/// cut edge of the shared face — connect the same two vertices. Four quads, one
/// dual edge, and each quad puts one of its two triangles on each side of it:
/// **four faces on an edge that should carry two**.
///
/// Marching Cubes reads the same twelve samples and comes out manifold, which is
/// ✗19 in a single fixture: Schaefer, Ju & Warren's premise — *"the original MC
/// algorithm always constructs a manifold"* — holds here, and their conclusion
/// — *"the dual preserves the topology of the surface"* — does not.
///
/// Plain Dual Contouring is measured beside it because the paper predicts it
/// (§3, *"DC leads to nonmanifold vertices and edges for all of the ambiguous
/// sign configurations"*), and on this configuration the manifold construction
/// buys **nothing**: splitting a cell by cycle cannot split a cell that has one.
#[test]
fn a_constructed_ambiguous_face_makes_the_dual_non_manifold() {
    use crate::cube::edge_on_face;
    use crate::marching_cubes::trilinear::Contours;

    // ── first, from the tables alone: this pattern offends, and why ──
    let (case_a, case_b) = (column_case(OFFENDER, 0), column_case(OFFENDER, 1));
    assert_eq!(
        (case_a, case_b),
        (103, 118),
        "the fixture's two cases moved"
    );
    assert_ne!(
        AMBIGUOUS_FACES[case_a as usize] & (1 << 5),
        0,
        "case {case_a}'s +z face is not ambiguous, so there is nothing to demonstrate"
    );
    let (owners_a, owners_b) = (cycle_owners(case_a, 0), cycle_owners(case_b, 0));
    for (owner, side) in [(&owners_a, 1u8), (&owners_b, 0u8)] {
        let cut: Vec<u8> = (0..12u8)
            .filter(|&e| edge_on_face(e, 2, side) && owner[e as usize] != 255)
            .collect();
        assert_eq!(cut.len(), 4, "an ambiguous face has all four edges cut");
        assert!(
            cut.iter()
                .all(|&e| owner[e as usize] == owner[cut[0] as usize]),
            "the four cut edges must land in one cycle — that is the mechanism"
        );
    }
    assert_eq!(
        (
            Contours::of(case_a, 0).count(),
            Contours::of(case_b, 0).count()
        ),
        (1, 1),
        "M-292: the default's eighteen offenders are exactly the (1, 1) bucket"
    );

    // ── then the meshes ──
    let field = column_lattice(OFFENDER);
    let shape = RuntimeShape3::new([NX as u32, NY as u32, NZ as u32]).expect("valid shape");
    let mut mdc = MeshBuffer::<f64>::new();
    ManifoldDualContouring::<f64>::new()
        .extract(&field, &shape, [0.0; 3], 1.0, &mut mdc)
        .expect("extraction");
    let mut dc = MeshBuffer::<f64>::new();
    DualContouring::<f64>::new()
        .extract(&field, &shape, [0.0; 3], 1.0, &mut dc)
        .expect("extraction");
    let mut mc = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(&field, &shape, [0.0; 3], 1.0, &mut mc)
        .expect("extraction");

    // **Marching Cubes, on the same twelve samples, is manifold.** The premise
    // of the paper's argument holds, which is what makes the rest of this a
    // falsification of the conclusion rather than of the premise.
    assert!(
        offending_edges(&mc).is_empty(),
        "Marching Cubes is non-manifold here, so the fixture would say nothing about the dual"
    );

    // **The dual is not**, and by the predicted mechanism exactly.
    let bad = offending_edges(&mdc);
    assert_eq!(bad.len(), 1, "expected the one shared face to offend, once");
    let (edge, faces) = bad.iter().next().expect("just asserted one");
    assert_eq!(
        faces.len(),
        4,
        "the mechanism is four *distinct* quads on one dual edge; fewer distinct \
         triangles would be duplication, which is a different defect"
    );

    // Its two ends are the two cells' vertices, one on each side of the shared
    // face — which is the claim "all four quads connect the same pair".
    let (p, q) = (
        mdc.positions[edge.0 as usize],
        mdc.positions[edge.1 as usize],
    );
    let (below, above) = if p[2] < q[2] { (p, q) } else { (q, p) };
    for v in [below, above] {
        assert!(
            (1.0..=2.0).contains(&v[0]) && (1.0..=2.0).contains(&v[1]),
            "a dual vertex escaped the column in x or y: {v:?}"
        );
    }
    assert!(
        below[2] < 1.0 && above[2] > 1.0,
        "the two ends should straddle the shared face at z = 1: {below:?} {above:?}"
    );

    // **Plain Dual Contouring measures the same defect**, which is the paper's
    // own prediction for it and the sharpest statement of what the manifold
    // construction is worth on this configuration: nothing.
    assert_eq!(
        offending_edges(&dc).len(),
        1,
        "one vertex per cell should offend here too — the paper says so in as many words"
    );
    assert_eq!(
        (
            mdc.triangle_count(),
            dc.triangle_count(),
            mc.triangle_count()
        ),
        (20, 20, 40),
        "the constructed fixture moved"
    );
    std::println!(
        "measured: 48 samples — MDC and DC each 1 non-manifold edge with 4 distinct faces, MC 0"
    );
}

/// **The same signs, and the asymptotic decider gives two different answers —
/// which is why no face rule closes A-025.**
///
/// [`a_constructed_ambiguous_face_makes_the_dual_non_manifold`] fixes the twelve
/// *signs*. The decider does not read signs; it reads the ambiguous face's four
/// corner *magnitudes* and asks where the bilinear saddle sits,
///
/// ```text
/// s = (v₀₀·v₁₁ − v₁₀·v₀₁) / (v₀₀ + v₁₁ − v₁₀ − v₀₁)
/// ```
///
/// — Nielson & Hamann, *The Asymptotic Decider* (Visualization '91). Scaling the
/// face's two **inside** corners, while leaving every sign alone, walks `s`
/// across zero, and the defect appears and disappears with it:
///
/// | inside corners | saddle `s` | the face | non-manifold edges |
/// |---|---|---|---|
/// | `−0.25` | `+0.375` | separated | **1** |
/// | `−1` (the fixture) | `0` — an exact tie | separated | **1** |
/// | `−4` | `−1.5` | **joined** | **0** |
///
/// Three consequences, and the third is the one A-025 is about.
///
/// **The tie is resolved, not undefined.** At `s = 0` this crate answers
/// *separated* (`ambiguity::face_is_joined`), so the perfectly symmetric saddle
/// — the one a hand-built fixture reaches first — is exactly the case the
/// decider declines to fix.
///
/// **The triangle count never moves.** Twenty either way: joining the face
/// changes which vertices the quads connect, not how many there are.
///
/// **The offending set is not a set of sign configurations.** M-292 enumerated
/// all 4,096 and found none that offends under every consistent mask; this is
/// what that looks like from the other side — one configuration, two answers,
/// chosen by the field rather than by the rule. Which is why the decider still
/// leaves 25–49 offending faces per resolution on `noise_cavity` (M-291) while
/// being combinatorially capable of avoiding all of them.
#[test]
fn the_decider_fixes_the_constructed_case_only_when_the_magnitudes_break_the_tie() {
    use crate::cube::{corner_offset, is_inside};
    use crate::marching_cubes::FaceAmbiguity;
    use crate::marching_cubes::ambiguity::joined_mask;

    /// Bit 5 of a face mask is the `+z` face — the one the column shares.
    const SHARED_FACE: u8 = 1 << 5;

    let mut rows = Vec::new();
    for scale in [0.25f64, 1.0, 4.0] {
        // Scale the shared face's two inside corners. No sign is touched.
        let mut field = column_lattice(OFFENDER);
        field.values[1][2][1] *= scale;
        field.values[2][1][1] *= scale;

        // What the decider makes of cell A's `+z` face.
        let mut corner = [0.0f64; 8];
        let mut case = 0u8;
        for (c, slot) in corner.iter_mut().enumerate() {
            let o = corner_offset(c as u8);
            *slot = field.values[1 + o[0] as usize][1 + o[1] as usize][o[2] as usize];
            if is_inside(*slot) {
                case |= 1 << c;
            }
        }
        assert_eq!(case, 103, "scaling a magnitude must not move a sign");
        let joined = joined_mask(&corner, AMBIGUOUS_FACES[case as usize]) & SHARED_FACE != 0;

        let shape = RuntimeShape3::new([NX as u32, NY as u32, NZ as u32]).expect("valid shape");
        let mut out = MeshBuffer::<f64>::new();
        let mut extractor = ManifoldDualContouring::<f64>::new();
        extractor.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
        extractor
            .extract(&field, &shape, [0.0; 3], 1.0, &mut out)
            .expect("extraction");
        assert_eq!(
            out.triangle_count(),
            20,
            "the face rule changes connectivity, not triangle count"
        );
        rows.push((scale, joined, report_of(&out, 1.0).non_manifold_edges));
    }

    assert_eq!(
        rows,
        alloc::vec![(0.25, false, 1), (1.0, false, 1), (4.0, true, 0)],
        "the decider's answer, and the defect with it, follows the magnitudes"
    );
    std::println!("measured: same signs, saddle crossing zero — {rows:?}");
}
