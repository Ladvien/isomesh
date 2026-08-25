//! **P-16 — are the past-90° normals the two faces meeting inside one cell?**
//!
//! Ticket: R-008. Pre-registered at R-008's own first commit.
//!
//! ```bash
//! cargo bench --bench experiment_p16
//! ```
//!
//! Writes `docs/experiments/p-16.csv`.
//!
//! # What M-283 left
//!
//! On an exact convex wedge, **6,959 vertices under Marching Cubes and 4,868
//! under Dual Contouring** carry an area-weighted normal more than 90° from the
//! field gradient, worst 128.0°. That should not be possible: an area-weighted
//! normal is a convex combination of incident face normals, so it lies inside
//! the cone those faces span, and two planes meeting at `θ` span at most
//! `(180° − θ)/2 ≤ 75°`. The escape — that the vertex is not on the surface,
//! since Marching Cubes interpolates linearly and a wedge is not linear at its
//! apex — was measured and closed: the median `|f(v)|/h` there is `0.0000`.
//!
//! So some incident face is not on either plane. **Where does it come from?**
//!
//! # The classification, and why it is about cells rather than triangles
//!
//! A triangle's *own* normal being off the two planes is the observation, not an
//! explanation — a triangle bridging a crease is expected to face somewhere
//! between. The question is whether such triangles come from cells the crease
//! passes through, which would make them inherent to one vertex per crossed
//! edge, or from cells wholly on one side, which would make them a winding or
//! ordering defect with a fix.
//!
//! A cell **straddles** when its eight corners do not all have the same nearer
//! plane. That is read off the field rather than off the mesh, so the
//! classification cannot inherit a mistake from the thing being classified.
//!
//! A vertex's incident cell is found from each incident *triangle*: every
//! Marching Cubes triangle lies inside one cell, and a triangle's centroid is
//! interior to it — which is why the centroid is used and the vertices are not.
//! M-49 is about `cell_of` on a cell **corner**, and this never asks it that.
//!
//! # The registered definition is wrong, and this file reports both
//!
//! P-16's registration says *"straddles the crease — its eight corners do not
//! all have the same nearer plane"*. The clause after the dash does not
//! implement the claim before it: `d0 ≥ d1` splits by the **bisector** plane,
//! which runs from the apex through the interior of the solid, not by the
//! crease. A cell deep inside the wedge and nowhere near the surface straddles
//! the bisector; a cell holding the crease may not.
//!
//! A registration is not edited after its experiment has run, so the registered
//! test stays and is reported as `past90_straddling_share`. The corrected one —
//! does the cell hold the **crease line** — is reported beside it as
//! `past90_on_crease_share` and is explicitly **post-hoc**: it is a measurement,
//! not evidence for a pre-registered claim. M-288 says which is which.
//!
//! # The control
//!
//! Vertices **under** 90° are classified the same way in the same run. If nearly
//! all vertices are in straddling cells then the hypothesis is vacuous, and the
//! `control_straddling_share` column is what says so.

mod common;

use isomesh::extractor::Extractor;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::normals::{NormalStrategy, recompute};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

use common::wedge::{Scalar, Wedge};

/// Dihedrals swept. `180` is dropped: there is no crease and M-283 measured
/// exactly zero disagreement there, so it has nothing to classify.
const DIHEDRALS: [f64; 7] = [30.0, 45.0, 60.0, 90.0, 120.0, 150.0, 170.0];

/// Rotations about the crease. M-283's lesson: a fixture on a grid axis reports
/// the prediction exactly and is wrong.
const ROTATIONS: [f64; 3] = [0.0, 17.0, 37.0];

/// Samples per axis.
const RESOLUTIONS: [u32; 3] = [33, 65, 129];

/// Domain half-extent, matching the crate's compact reference fields.
const HALF_EXTENT: f64 = 2.0;

/// Apex offset in cells — irrational, so the crease never lands on a sample
/// plane at any resolution in the sweep.
const OFFSET_CELLS: f64 = 0.447_213_595_499_957_9;

/// Vertices within this many cells of the domain boundary are excluded; a
/// boundary-clipped vertex has fewer incident faces than it should.
const BOUNDARY_MARGIN_CELLS: f64 = 1.0;

/// Which of the wedge's two planes is nearer at `p`.
///
/// Read from the field, not the mesh. `d0 >= d1` matches the tie-break
/// `Wedge::gradient` uses, so the two cannot disagree about a corner.
fn nearer_plane(field: &Wedge, p: [Scalar; 3]) -> bool {
    let q = field.local(p);
    let [n0, n1] = field.plane_normals();
    q[0] * n0[0] + q[1] * n0[1] >= q[0] * n1[0] + q[1] * n1[1]
}

/// Does the cell whose minimum corner is `lo` hold the **crease line**?
///
/// The corrected test, and it is **post-hoc** — see the module docs and M-288.
/// The crease is the line through the apex along `z`, and rotation is *about*
/// that line, so it does not move: a cell holds it exactly when the cell's `x`
/// and `y` ranges both contain the apex. `within` widens that by a number of
/// cells, because a vertex's incident triangles can come from a neighbour.
fn holds_crease(field: &Wedge, lo: [Scalar; 3], h: Scalar, within: Scalar) -> bool {
    let pad = within * h;
    (lo[0] - pad..lo[0] + h + pad).contains(&field.apex[0])
        && (lo[1] - pad..lo[1] + h + pad).contains(&field.apex[1])
}

/// Does the cell whose minimum corner is `lo` have corners on both sides?
///
/// **This is P-16's registered operational definition and it does not implement
/// P-16's claim** — see the module docs. `d0 ≥ d1` splits by the *bisector*
/// plane, which cuts the solid's interior from the apex outward, not by the
/// crease. Kept, and reported, because a registration is not edited after its
/// experiment runs.
fn straddles(field: &Wedge, lo: [Scalar; 3], h: Scalar) -> bool {
    let first = nearer_plane(field, lo);
    for corner in 0..8u8 {
        let p = [
            lo[0] + if corner & 1 != 0 { h } else { 0.0 },
            lo[1] + if corner & 2 != 0 { h } else { 0.0 },
            lo[2] + if corner & 4 != 0 { h } else { 0.0 },
        ];
        if nearer_plane(field, p) != first {
            return true;
        }
    }
    false
}

/// What one configuration produced.
#[derive(Default)]
struct Tally {
    past90: usize,
    past90_straddling: usize,
    past90_on_crease: usize,
    control: usize,
    control_straddling: usize,
    control_on_crease: usize,
    offending_faces: usize,
    /// Distance from the crease line to each past-90° vertex, in cells.
    crease_distance: Vec<Scalar>,
    /// Distance from the nearest domain wall to each past-90° vertex, in cells.
    wall_distance: Vec<Scalar>,
}

/// One row.
fn measure(dihedral: f64, samples: u32, rotation: f64, run: &mut common::experiment::Run) {
    let cell_size = 2.0 * HALF_EXTENT / f64::from(samples - 1);
    let shift = OFFSET_CELLS * cell_size;
    let field = Wedge::new(dihedral, [-HALF_EXTENT / 2.0 + shift, shift], rotation);
    let origin = [-HALF_EXTENT; 3];

    let shape = RuntimeShape3::new([samples; 3]).expect("the fixture fits u32");
    let mut mesh = MeshBuffer::<Scalar>::new();
    MarchingCubes::<Scalar>::new()
        .extract_into(&field, &shape, origin, cell_size, &mut mesh)
        .expect("extraction");
    recompute(&mut mesh, &field, NormalStrategy::AreaWeightedFaces).expect("normals");

    // Incident triangles per vertex.
    let mut incident: Vec<Vec<u32>> = vec![Vec::new(); mesh.vertex_count()];
    for (t, tri) in mesh.indices.as_chunks::<3>().0.iter().enumerate() {
        for &v in tri {
            incident[v as usize].push(t as u32);
        }
    }

    // The cone the two planes span, as a half-angle about their bisector. A
    // face outside it is on neither plane.
    let cone = (180.0 - dihedral) / 2.0;
    let [n0, n1] = field.plane_normals();
    let planes = [field.unrotate(n0), field.unrotate(n1)];

    let mut tally = Tally::default();
    // Where the offenders are, asked directly. The classifications above say
    // where they are *not*, and two negatives do not locate anything.
    let limit = HALF_EXTENT - BOUNDARY_MARGIN_CELLS * cell_size;
    for (v, incident_here) in incident.iter().enumerate() {
        let position = mesh.positions[v];
        if position.iter().any(|c| c.abs() > limit) {
            continue;
        }
        let normal = mesh.normals[v];
        let g = field.gradient(position);
        let gl = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        let nl = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        if gl == 0.0 || nl == 0.0 {
            continue;
        }
        let dot = (normal[0] * g[0] + normal[1] * g[1] + normal[2] * g[2]) / (gl * nl);
        let angle = dot.clamp(-1.0, 1.0).acos().to_degrees();

        let mut any_straddling = false;
        let mut any_on_crease = false;
        let mut offending = 0usize;
        for &t in incident_here {
            let tri = &mesh.indices[t as usize * 3..t as usize * 3 + 3];
            let p = [
                mesh.positions[tri[0] as usize],
                mesh.positions[tri[1] as usize],
                mesh.positions[tri[2] as usize],
            ];
            let centre = [
                (p[0][0] + p[1][0] + p[2][0]) / 3.0,
                (p[0][1] + p[1][1] + p[2][1]) / 3.0,
                (p[0][2] + p[1][2] + p[2][2]) / 3.0,
            ];
            // A triangle's centroid is interior to its cell, which is why the
            // centroid is used here and a vertex is not — M-49 is about
            // `cell_of` on a cell corner and this never asks it that.
            let mut lo = [0.0; 3];
            for (axis, slot) in lo.iter_mut().enumerate() {
                let index = ((centre[axis] - origin[axis]) / cell_size).floor();
                *slot = origin[axis] + cell_size * index;
            }
            if straddles(&field, lo, cell_size) {
                any_straddling = true;
            }
            if holds_crease(&field, lo, cell_size, 1.0) {
                any_on_crease = true;
            }

            // Is this face's normal outside the cone the two planes span?
            let e1 = [p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]];
            let e2 = [p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]];
            let fnorm = [
                e1[1] * e2[2] - e1[2] * e2[1],
                e1[2] * e2[0] - e1[0] * e2[2],
                e1[0] * e2[1] - e1[1] * e2[0],
            ];
            let fl = (fnorm[0] * fnorm[0] + fnorm[1] * fnorm[1] + fnorm[2] * fnorm[2]).sqrt();
            if fl > 0.0 {
                let worst = planes
                    .iter()
                    .map(|n| {
                        let d = (fnorm[0] * n[0] + fnorm[1] * n[1] + fnorm[2] * n[2]) / fl;
                        d.clamp(-1.0, 1.0).acos().to_degrees()
                    })
                    .fold(f64::INFINITY, f64::min);
                // A degree of slack, so a face lying on a plane to within
                // rounding is not counted as leaving it.
                if worst > cone + 1.0 {
                    offending += 1;
                }
            }
        }

        if angle > 90.0 {
            tally.past90 += 1;
            tally.offending_faces += offending;
            if any_straddling {
                tally.past90_straddling += 1;
            }
            if any_on_crease {
                tally.past90_on_crease += 1;
            }
            // Distance to the crease line, which runs along `z` through the
            // apex, and to the nearest domain wall — both in cells.
            let dx = position[0] - field.apex[0];
            let dy = position[1] - field.apex[1];
            tally
                .crease_distance
                .push((dx * dx + dy * dy).sqrt() / cell_size);
            let wall = position
                .iter()
                .map(|c| HALF_EXTENT - c.abs())
                .fold(f64::INFINITY, f64::min);
            tally.wall_distance.push(wall / cell_size);
        } else {
            tally.control += 1;
            if any_straddling {
                tally.control_straddling += 1;
            }
            if any_on_crease {
                tally.control_on_crease += 1;
            }
        }
    }

    let share = |a: usize, b: usize| if b == 0 { 0.0 } else { a as f64 / b as f64 };
    let median = |v: &mut Vec<Scalar>| {
        if v.is_empty() {
            return 0.0;
        }
        v.sort_by(f64::total_cmp);
        v[v.len() / 2]
    };
    let crease_median = median(&mut tally.crease_distance);
    let wall_median = median(&mut tally.wall_distance);
    println!(
        "{dihedral:>8.0} {rotation:>4.0} {samples:>5} {:>8} {:>10.4} {:>10.4} {:>11.2} {:>10.2}",
        tally.past90,
        share(tally.past90_straddling, tally.past90),
        share(tally.past90_on_crease, tally.past90),
        crease_median,
        wall_median
    );

    run.record(&[
        ("dihedral_deg", format!("{dihedral:.1}")),
        ("samples", samples.to_string()),
        ("past90_vertices", tally.past90.to_string()),
        (
            "past90_in_straddling_cell",
            tally.past90_straddling.to_string(),
        ),
        (
            "offending_faces_per_past90_vertex",
            format!("{:.4}", share(tally.offending_faces, tally.past90)),
        ),
        ("rotation_deg", format!("{rotation:.0}")),
        (
            "past90_straddling_share",
            format!("{:.6}", share(tally.past90_straddling, tally.past90)),
        ),
        ("control_vertices", tally.control.to_string()),
        (
            "control_straddling_share",
            format!("{:.6}", share(tally.control_straddling, tally.control)),
        ),
        (
            "past90_on_crease_share",
            format!("{:.6}", share(tally.past90_on_crease, tally.past90)),
        ),
        (
            "control_on_crease_share",
            format!("{:.6}", share(tally.control_on_crease, tally.control)),
        ),
        (
            "past90_median_crease_distance_cells",
            format!("{crease_median:.4}"),
        ),
        (
            "past90_median_wall_distance_cells",
            format!("{wall_median:.4}"),
        ),
        ("triangles", mesh.triangle_count().to_string()),
    ]);
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-16");
    common::experiment::run(prereg, |run| {
        println!(
            "{:>8} {:>4} {:>5} {:>8} {:>10} {:>10} {:>9} {:>10}",
            "dihedral", "rot", "n", "past90", "registered", "on crease", "d(crease)", "d(wall)"
        );
        for dihedral in DIHEDRALS {
            for rotation in ROTATIONS {
                for samples in RESOLUTIONS {
                    measure(dihedral, samples, rotation, run);
                }
            }
        }
        println!(
            "\nA cell straddles when its eight corners do not all have the same nearer plane, read \
             from\nthe field rather than the mesh. `ctl share` is the same classification over the \
             vertices\n**under** 90°, so a hypothesis that is true of every vertex is visible as \
             such."
        );
    });
}
