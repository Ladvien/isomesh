use super::*;
use crate::extractor::Extractor;
use crate::fields::{ReferenceField, Sphere};
use crate::{MeshBuffer, RuntimeShape3};

/// Samples per axis. Odd, so the sphere's centre is a sample and the fixture is
/// symmetric about it.
const N: u32 = 17;

/// The `sphere` grid, spelled the way the benches spell it.
fn grid() -> (Sphere<f64>, RuntimeShape3, [f64; 3], f64) {
    let field = Sphere::<f64>::canonical();
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(N - 1);
    let shape = RuntimeShape3::new([N; 3]).expect("valid shape");
    (field, shape, lo, h)
}

fn mesh_with<E: Extractor<f64>>(extractor: &mut E) -> MeshBuffer<f64> {
    let (field, shape, lo, h) = grid();
    let mut out = MeshBuffer::<f64>::new();
    extractor
        .extract_into(&field, &shape, lo, h, &mut out)
        .expect("extraction");
    out
}

fn report_of(mesh: &MeshBuffer<f64>) -> SealingReport {
    let (field, shape, lo, h) = grid();
    sealing(&field, &shape, lo, h, &mesh.positions, &mesh.indices)
}

/// **R-024's headline, and it is the one that can be wrong.**
///
/// Every extractor in the registry, on the field where the answer is least
/// ambiguous. A disagreement here would be a defect in the extractor, in the
/// harness, or in the claim — and the controls below exist so those can be told
/// apart.
#[test]
fn every_extractor_seals_a_sphere() {
    crate::for_each_extractor!(f64, |name, extractor| {
        let mesh = mesh_with(&mut extractor);
        assert!(
            !mesh.indices.is_empty(),
            "{name} meshed a sphere to nothing"
        );
        let r = report_of(&mesh);
        assert!(
            r.agrees(),
            "{name} does not seal a sphere:\n{r}\nunsealed {} spurious {} mixed {}",
            r.unsealed_walls,
            r.spurious_walls,
            r.mixed_regions
        );
    });
}

/// **The instrument has to be able to report the bad news** — E-208's rule, and
/// the reason a zero above means anything.
///
/// With no triangles at all, every probe reports zero crossings, so the mesh
/// separates nothing and every wall the field asserts is unsealed.
#[test]
fn an_empty_mesh_leaves_every_wall_unsealed() {
    let empty = MeshBuffer::<f64>::new();
    let r = report_of(&empty);

    assert!(r.field_walls > 0, "the fixture has no walls to miss");
    assert_eq!(r.unsealed_walls, r.field_walls);
    assert_eq!(r.mesh_walls, 0);
    assert_eq!(r.spurious_walls, 0);
    // Seven regions, and the six extra ones are not the mesh's doing: the six
    // samples the surface passes exactly through have all their probes set
    // aside, so each is isolated in **both** graphs. What is left is the one big
    // region, and it holds air and solid alike because nothing cut it.
    assert_eq!(r.mesh_regions, 1 + r.boundary_samples);
    assert_eq!(r.mixed_regions, 1);
    assert!(!r.agrees());
}

/// **A mesh in the wrong place, which is the sharpest control available.**
///
/// The same sphere mesh displaced by a non-multiple of the cell size. Nothing
/// about it is malformed — it is closed, manifold and correctly wound — and it
/// is wrong about the field in both directions at once, which is exactly what
/// distinguishes this report from every other one in `validate`.
#[test]
fn a_displaced_mesh_disagrees_in_both_directions() {
    let mut mesh = mesh_with(&mut crate::marching_cubes::MarchingCubes::<f64>::new());
    for p in &mut mesh.positions {
        p[0] += 0.3;
    }
    let r = report_of(&mesh);

    assert!(
        r.unsealed_walls > 0,
        "a displaced sphere leaves the true boundary unsealed:\n{r}"
    );
    assert!(
        r.spurious_walls > 0,
        "a displaced sphere walls off open air:\n{r}"
    );
    assert!(!r.agrees());
}

/// **One triangle across open air is a membrane, and is reported as one.**
///
/// Placed so that it contains exactly one grid sample's outgoing `x` probe and
/// comes nowhere near another: a small triangle in the plane `x = −1.625`,
/// around the sample at `(−1.75, −1.75)` in `(y, z)`, which is 3.03 from the
/// origin and so deep in air.
#[test]
fn a_membrane_across_open_air_is_a_spurious_wall() {
    let mut mesh = mesh_with(&mut crate::marching_cubes::MarchingCubes::<f64>::new());
    let clean = report_of(&mesh);
    assert!(clean.agrees(), "the fixture must start sealed:\n{clean}");

    let base = mesh.positions.len() as u32;
    for p in [
        [-1.625, -1.85, -1.80],
        [-1.625, -1.65, -1.80],
        [-1.625, -1.75, -1.65],
    ] {
        mesh.positions.push(p);
        mesh.normals.push([1.0, 0.0, 0.0]);
    }
    mesh.indices.extend_from_slice(&[base, base + 1, base + 2]);

    let r = report_of(&mesh);
    assert_eq!(
        r.spurious_walls, 1,
        "one membrane should wall off exactly one probe:\n{r}"
    );
    assert_eq!(r.unsealed_walls, 0, "the sphere itself is untouched:\n{r}");
    assert!(!r.agrees());
}

/// **The merge is doing real work, and the amount is a property of the family.**
///
/// A primal method places its vertex *on* the probed edge and its whole fan
/// contains that point, so raw triangle hits far outnumber distinct crossings. A
/// dual method places its vertex in the cell interior and its quads cross
/// transversally, so almost nothing merges. If this ever inverted, the
/// deduplication would be papering over something instead of resolving a known
/// degeneracy — see the module docs.
#[test]
fn the_merge_separates_the_primal_family_from_the_dual_one() {
    let primal = report_of(&mesh_with(
        &mut crate::marching_cubes::MarchingCubes::<f64>::new(),
    ));
    let dual = report_of(&mesh_with(
        &mut crate::dual_contouring::DualContouring::<f64>::new(),
    ));

    assert_eq!(
        primal.field_walls, dual.field_walls,
        "same field, same walls"
    );
    assert!(
        primal.merged_crossings > primal.field_walls,
        "a primal method should merge more than one hit per wall, got {} for {} walls",
        primal.merged_crossings,
        primal.field_walls
    );
    assert!(
        dual.merged_crossings * 5 < dual.field_walls,
        "a dual method should barely merge at all, got {} for {} walls",
        dual.merged_crossings,
        dual.field_walls
    );
}

/// **A sample exactly on the surface makes every probe touching it undecidable,
/// and the count is derivable rather than observed.**
///
/// The unit sphere on a `0.25` lattice over `[−2, 2]` passes exactly through its
/// six axis intercepts and nowhere else on the grid: `(±1, 0, 0)` and its
/// permutations. `a² + b² + c² = 1` in multiples of `0.25` means
/// `A² + B² + C² = 16` in integers, whose only solutions with `|A|, |B|, |C| ≤ 8`
/// are the permutations of `(±4, 0, 0)`.
///
/// Each is interior, so each has six incident grid edges, and `6 × 6 = 36`
/// probes are set aside. **The exclusion is symmetric** — both graphs lose the
/// same edges — which is why the two component counts still agree while each
/// gains the six isolated samples.
#[test]
fn a_sample_on_the_surface_makes_its_probes_undecidable() {
    let (field, _, lo, h) = grid();
    let mut on_surface = alloc::vec::Vec::new();
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let p = [
                    lo[0] + h * f64::from(x),
                    lo[1] + h * f64::from(y),
                    lo[2] + h * f64::from(z),
                ];
                #[allow(clippy::float_cmp, reason = "the exact tie is the subject")]
                let on = field.sample(p) == 0.0;
                if on {
                    on_surface.push(p);
                }
            }
        }
    }
    assert_eq!(on_surface.len(), 6, "expected the six axis intercepts");
    for p in &on_surface {
        #[allow(clippy::float_cmp, reason = "an intercept's other two are exactly 0")]
        let off_axis = p.iter().filter(|c| **c != 0.0).count();
        assert_eq!(off_axis, 1, "{p:?} is not an axis intercept");
    }

    let r = report_of(&mesh_with(
        &mut crate::marching_cubes::MarchingCubes::<f64>::new(),
    ));
    assert_eq!(r.boundary_samples, 6);
    assert_eq!(
        r.degenerate_probes, 36,
        "six interior samples, six edges each"
    );
    assert_eq!(
        r.field_air_components, r.mesh_air_components,
        "the exclusion has to be symmetric or it manufactures a difference"
    );
}

/// **A zero-area triangle answers "parallel" to every probe, so it is excluded
/// before binning rather than judged.**
///
/// Its normal is a cancellation residue, so the parallel branch fires whatever
/// the probe is, and it lands in `coplanar_probes` once per probe it is binned
/// near. Marching Tetrahedra collapses a triangle at each of the six samples the
/// surface passes through — the same six as the test above, `6 × 6 = 36` — and
/// before the exclusion those 36 slivers produced **6,624** coplanar events
/// between them while every other extractor produced none.
#[test]
fn a_zero_area_triangle_is_excluded_before_it_can_be_judged() {
    let mt = report_of(&mesh_with(
        &mut crate::marching_tetrahedra::MarchingTetrahedra::<f64>::new(),
    ));
    assert_eq!(mt.degenerate_triangles, 36);
    assert_eq!(
        mt.coplanar_probes, 0,
        "the exclusion is what takes this to zero; 6624 without it"
    );

    let mc = report_of(&mesh_with(
        &mut crate::marching_cubes::MarchingCubes::<f64>::new(),
    ));
    assert_eq!(mc.degenerate_triangles, 0);
    assert_eq!(mc.coplanar_probes, 0);
}

/// Same mesh twice, same report. The union-find's representative is fixed by the
/// lattice scan and nothing here iterates a map.
#[test]
fn sealing_is_deterministic() {
    let mesh = mesh_with(&mut crate::surface_nets::SurfaceNets::<f64>::new());
    assert_eq!(report_of(&mesh), report_of(&mesh));
}

/// **A dual method leaves the domain boundary unsealed, and only there.**
///
/// R-024's sharpest result, pinned. A dual emits one quad per sign-changing grid
/// edge and that quad needs all **four** cells around the edge. On a face of the
/// sampled domain only one or two exist, so no quad is emitted and the wall
/// stays open — while a primal method, which emits per *cell*, meshes every cell
/// it has and seals the same edge.
///
/// `fbm_terrain` is the only reference field whose surface leaves through the
/// sides, so it is the only one where this is reachable. All three duals report
/// the identical count, which is what says the mechanism is one-vertex-per-cell
/// rather than any particular solve.
///
/// **For a chunked world this is the chunk seam**, and it is why a dual chunk's
/// collider is not watertight on its own.
#[test]
fn a_dual_leaves_the_domain_boundary_unsealed_and_only_there() {
    use crate::fields::FbmTerrain;

    const M: u32 = 17;
    let field = FbmTerrain::<f64>::canonical();
    let (flo, fhi) = field.domain();
    let h = (fhi[0] - flo[0]) / f64::from(M - 1);
    let shape = RuntimeShape3::new([M; 3]).expect("valid shape");

    let mut duals = 0;
    crate::for_each_extractor!(f64, |name, extractor| {
        let mut out = MeshBuffer::<f64>::new();
        extractor
            .extract_into(&field, &shape, flo, h, &mut out)
            .expect("extraction");
        let r = sealing(&field, &shape, flo, h, &out.positions, &out.indices);

        // A method that places its crossing *on* the probed grid edge seals
        // everything; the split is that property, not primal-versus-dual.
        if name.starts_with("marching") {
            assert_eq!(r.unsealed_walls, 0, "{name} left a wall open:\n{r}");
        }
        if name.ends_with("dual_contouring") || name == "surface_nets" {
            duals += 1;
            assert_eq!(r.unsealed_walls, 92, "{name}:\n{r}");
            assert_eq!(
                r.unsealed_on_domain_face, r.unsealed_walls,
                "{name} left a hole away from the domain face, which is a \
                 different defect:\n{r}"
            );
        }
    });
    assert_eq!(duals, 3, "expected three dual entries in the registry");
}
