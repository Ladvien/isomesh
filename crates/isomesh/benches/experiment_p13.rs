//! **P-13 — is M-66's non-convergent angle a property of the feature?**
//!
//! Ticket: R-006. Pre-registered at R-000.
//!
//! ```bash
//! cargo bench --bench experiment_p13
//! ```
//!
//! Writes `docs/experiments/p-13.csv`.
//!
//! # The one error in this crate that does not converge
//!
//! Everything else falls with `h`. M-12: position error at `h²`. M-65: normal
//! direction at `h²`, 0.460° → 0.031° over 17³…65³. Then **M-66**: on
//! `box_exact` the *mean* angle between the area-weighted face normal and the
//! analytic gradient falls 13.55° → 6.73° → 3.34°, and the **worst is 35.796°
//! at all three resolutions, identical to six figures**.
//!
//! An error that does not fall with resolution is either a real property of
//! sharp features or a defect, and the two are worth telling apart.
//!
//! # The fixture the ticket asks for, and why a box could not answer it
//!
//! `box_exact` has one dihedral — 90° — and also has corners, where *three*
//! faces meet and the geometry is not the same problem. A single number at a
//! single angle cannot show that the angle is *set by* the dihedral.
//!
//! So [`Wedge`] is an exact signed distance field for a convex wedge of
//! controllable dihedral: one crease, no corners, no second feature. Exact
//! rather than `max(d₁, d₂)` of two half-spaces, which is a Pseudo-SDF that
//! **underestimates outside the crease** (Phase 11) — exactly where the
//! extractor interpolates its crossings, so the error would land on the
//! measurement rather than on the field.
//!
//! # The prediction, made concrete
//!
//! Across a crease of dihedral `θ` the surface normal turns by `180° − θ`. A
//! vertex on the crease is given, by area weighting, something between the two
//! face normals, so its disagreement with either is at most **half that turn**:
//!
//! ```text
//! predicted worst angle = (180° − θ) / 2
//! ```
//!
//! For a box's 90° that is 45°, against M-66's measured 35.796° — so the
//! prediction is not expected to be exact, and what P-13 claims is that the
//! measurement **tracks** it and does not move with `h`.
//!
//! # The alignment control, because M-266 exists
//!
//! *"M-72's aliasing is alignment, not chance."* A crease sitting exactly on a
//! grid plane is resolved perfectly and would report an angle that is a property
//! of the fixture. Every row is therefore run twice — apex on a sample plane and
//! apex offset by an irrational fraction of a cell — and the CSV says which.
//!
//! **Rotation, because the first sweep confounded two things.** A wedge of
//! dihedral `θ` has its two faces at `±(90° − θ/2)` to the bisector, so **as `θ`
//! changes the faces change their angle to the grid**, and at `θ = 170°` both
//! are within 5° of a coordinate plane — the worst case for Marching Cubes,
//! which terraces an almost-axis-aligned plane and gives the terraces normals
//! nothing to do with the surface. The first sweep duly reported 75–93° where
//! the crease can only turn 5–30°, and that is the fixture, not the feature.
//!
//! So the whole wedge is also rotated about the crease by `0°`, `17°` and `37°`
//! — none of them special — and the CSV carries the rotation. A quantity that is
//! a property of the dihedral survives it; a staircase does not.
//!
//! # `worst` is not a robust statistic, so it is not the only one
//!
//! One degenerate vertex sets the maximum, and grid-scale artefacts produce
//! angles at `arctan 2 = 63.43°` and at 90° that have nothing to do with the
//! crease. `p99` and the median are recorded beside it.

mod common;

use isomesh::dual_contouring::DualContouring;
use isomesh::extractor::Extractor;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::normals::{NormalStrategy, recompute};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

use common::wedge::{Scalar, Wedge};

/// Dihedral angles swept, in degrees.
///
/// `180` is the control: the two half-planes coincide, there is no crease, and
/// whatever the harness reports there is the discretisation floor rather than
/// the feature.
const DIHEDRALS: [f64; 8] = [30.0, 45.0, 60.0, 90.0, 120.0, 150.0, 170.0, 180.0];

/// Samples per axis. A-012's three, plus one more so a trend has somewhere to
/// go.
const RESOLUTIONS: [u32; 4] = [17, 33, 65, 129];

/// Rotations of the whole wedge about its crease, in degrees.
///
/// None of them special, and the point of them is in the module docs: a wedge's
/// faces change their angle to the grid as its dihedral changes, so without this
/// the sweep measures orientation and dihedral at once.
const ROTATIONS: [f64; 3] = [0.0, 17.0, 37.0];

/// Domain half-extent, matching the crate's compact reference fields.
const HALF_EXTENT: f64 = 2.0;

/// How far the apex is moved off the grid in the offset arm, in cells.
///
/// `1/√5` — irrational, and not near a half or a quarter, so no resolution in
/// the sweep can accidentally land the crease back on a sample plane.
const OFFSET_CELLS: f64 = 0.447_213_595_499_957_9;

/// Vertices within this many cells of the domain boundary are excluded.
///
/// A boundary-clipped vertex has fewer incident faces than it should, so its
/// area-weighted normal is wrong for a reason that has nothing to do with the
/// crease. One cell is enough: the clipping is exactly one cell deep.
const BOUNDARY_MARGIN_CELLS: f64 = 1.0;

/// The distribution of the angle, in degrees, between each vertex's normal and
/// the field's own gradient there.
struct Disagreement {
    worst: f64,
    p99: f64,
    median: f64,
    mean: f64,
    vertices: usize,
    /// Vertices whose normal has turned past 90° from the gradient — pointing
    /// into the solid rather than merely in the wrong direction.
    inverted: usize,
    /// The largest `|f(v)| / h` over those past-90° vertices.
    ///
    /// **Because past 90° should be impossible and the reason matters.** An
    /// area-weighted normal is a convex combination of face normals, so it lies
    /// inside the cone they span — at most `(180° − θ)/2 ≤ 75°` from either
    /// plane normal. Exceeding 90° therefore means the comparison is not against
    /// a plane normal at all, and the likeliest reason is that the vertex is not
    /// on the surface: Marching Cubes interpolates linearly, and the wedge is
    /// **not** linear near its apex, so a crossing there can sit a fraction of a
    /// cell off the true surface where the gradient is something else again.
    /// This measures that fraction rather than assuming it.
    inverted_offsurface: f64,
}

fn disagreement(
    mesh: &MeshBuffer<Scalar>,
    field: &Wedge,
    margin: Scalar,
    cell_size: Scalar,
) -> Disagreement {
    let mut angles: Vec<f64> = Vec::with_capacity(mesh.positions.len());
    let mut inverted = 0usize;
    let mut inverted_offsurface = 0.0f64;
    let limit = HALF_EXTENT - margin;
    for (position, normal) in mesh.positions.iter().zip(&mesh.normals) {
        // A boundary-clipped vertex has fewer incident faces than it should, so
        // its area-weighted normal is wrong for a reason that is not the crease.
        if position.iter().any(|c| c.abs() > limit) {
            continue;
        }
        let g = field.gradient(*position);
        let gl = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        let nl = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        if gl == 0.0 || nl == 0.0 {
            continue;
        }
        let dot = (normal[0] * g[0] + normal[1] * g[1] + normal[2] * g[2]) / (gl * nl);
        let angle = dot.clamp(-1.0, 1.0).acos().to_degrees();
        if angle > 90.0 {
            inverted += 1;
            inverted_offsurface =
                inverted_offsurface.max((field.sample(*position) / cell_size).abs());
        }
        angles.push(angle);
    }
    if angles.is_empty() {
        return Disagreement {
            worst: 0.0,
            p99: 0.0,
            median: 0.0,
            mean: 0.0,
            vertices: 0,
            inverted: 0,
            inverted_offsurface: 0.0,
        };
    }
    angles.sort_by(f64::total_cmp);
    let at = |q: f64| angles[(((angles.len() - 1) as f64) * q).round() as usize];
    Disagreement {
        worst: angles[angles.len() - 1],
        p99: at(0.99),
        median: at(0.5),
        mean: angles.iter().sum::<f64>() / angles.len() as f64,
        vertices: angles.len(),
        inverted,
        inverted_offsurface,
    }
}

/// One row.
fn measure<E: Extractor<Scalar>>(
    extractor: &mut E,
    name: &str,
    dihedral: f64,
    samples: u32,
    aligned: bool,
    rotation: f64,
    run: &mut common::experiment::Run,
) {
    let cell_size = 2.0 * HALF_EXTENT / f64::from(samples - 1);
    // The apex sits on the `x` axis so both faces are symmetric about `y = 0`,
    // and — in the offset arm — off the sample lattice on both in-plane axes.
    let shift = if aligned {
        0.0
    } else {
        OFFSET_CELLS * cell_size
    };
    let field = Wedge::new(dihedral, [-HALF_EXTENT / 2.0 + shift, shift], rotation);

    let shape = RuntimeShape3::new([samples; 3]).expect("the fixture fits u32");
    let mut mesh = MeshBuffer::<Scalar>::new();
    extractor
        .extract_into(&field, &shape, [-HALF_EXTENT; 3], cell_size, &mut mesh)
        .expect("extraction");
    assert!(
        mesh.triangle_count() > 0,
        "θ = {dihedral} at {samples}³ meshed to nothing"
    );

    // M-66's strategy: the normal the *geometry* implies, not the field's.
    recompute(&mut mesh, &field, NormalStrategy::AreaWeightedFaces).expect("normals");
    let d = disagreement(&mesh, &field, BOUNDARY_MARGIN_CELLS * cell_size, cell_size);

    let predicted = (180.0 - dihedral) / 2.0;
    let arm = if aligned { "aligned" } else { "offset" };
    println!(
        "{dihedral:>8.0} {samples:>6} {arm:<8} {:>12.4} {:>12.4} {:>10.4} {:>9}",
        d.worst, predicted, d.mean, d.vertices
    );

    run.record(&[
        ("dihedral_deg", format!("{dihedral:.1}")),
        ("samples", samples.to_string()),
        ("measured_angle_deg", format!("{:.4}", d.worst)),
        ("predicted_angle_deg", format!("{predicted:.4}")),
        ("apex", arm.to_string()),
        ("extractor", name.to_string()),
        ("rotation_deg", format!("{rotation:.0}")),
        ("p99_angle_deg", format!("{:.4}", d.p99)),
        ("median_angle_deg", format!("{:.4}", d.median)),
        ("mean_angle_deg", format!("{:.4}", d.mean)),
        ("inverted_vertices", d.inverted.to_string()),
        (
            "inverted_offsurface_cells",
            format!("{:.4}", d.inverted_offsurface),
        ),
        ("vertices_counted", d.vertices.to_string()),
        ("triangles", mesh.triangle_count().to_string()),
        ("cell_size", format!("{cell_size:.17e}")),
    ]);
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-13");
    common::experiment::run(prereg, |run| {
        println!(
            "{:<16} {:>6} {:>4} {:>5} {:<8} {:>9} {:>9} {:>8} {:>8} {:>8} {:>6}",
            "extractor",
            "theta",
            "rot",
            "n",
            "apex",
            "worst",
            "predict",
            "p99",
            "median",
            "mean",
            "inv"
        );
        for dihedral in DIHEDRALS {
            for rotation in ROTATIONS {
                for aligned in [true, false] {
                    for samples in RESOLUTIONS {
                        measure(
                            &mut MarchingCubes::<Scalar>::new(),
                            "marching_cubes",
                            dihedral,
                            samples,
                            aligned,
                            rotation,
                            run,
                        );
                        measure(
                            &mut DualContouring::<Scalar>::new(),
                            "dual_contouring",
                            dihedral,
                            samples,
                            aligned,
                            rotation,
                            run,
                        );
                    }
                }
            }
        }
        println!(
            "\n`predicted` is (180° − θ)/2, half the turn the surface normal makes across the \
             crease.\n`worst` is the largest disagreement between an area-weighted vertex normal \
             and the field's\nown gradient there, over vertices more than one cell from the domain \
             boundary."
        );
    });
}
