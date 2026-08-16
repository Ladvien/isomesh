//! **R-009 — find the one-to-three cells per cross-section that face backwards.**
//!
//! ```bash
//! cargo bench --bench r009_locate
//! ```
//!
//! Writes `docs/measurements/r009-locate.csv`.
//!
//! # Not a sweep, and that is the point
//!
//! R-008 (M-288) failed to locate 80% of the past-90° normals and bounded them
//! tightly instead. Counts are exact multiples of `n − 2` at every resolution;
//! the wedge is extruded along `z`, so the sampling repeats identically in every
//! layer and any cross-section feature repeats once per layer. **There are one,
//! two or three offending locations in the whole two-dimensional cross-section**
//! — on a surface that is otherwise two exact planes, which Marching Cubes
//! reproduces to `0.0000°`.
//!
//! Two classifiers have now said where they are *not*. A-021 is the model for
//! what to do next: it found its answer by printing a face-count histogram for a
//! **plain half-space**, not by widening a census. So this dumps one
//! configuration and reads it.
//!
//! # What it prints, and why each column is there
//!
//! For every past-90° vertex in `θ = 150°`, rotation 17°, `n = 65` — 126 of
//! them, so 126 / (65 − 2) = **two per layer**:
//!
//! - **its cross-section cell**, so the repeat-per-layer claim is checked rather
//!   than assumed;
//! - **how many faces are incident**, because a vertex with a two-triangle fan
//!   has an area-weighted normal that is not an average of anything and is the
//!   first thing to rule out;
//! - **each incident face's angle to the nearer plane and its area**, because a
//!   sliver contributes a numerically meaningless normal at a weight that is
//!   supposed to make it harmless;
//! - **the cell's Marching Cubes case index**, which names the configuration if
//!   it is one.

mod common;

use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use isomesh::extractor::Extractor;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::normals::{NormalStrategy, recompute};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

use common::wedge::{Scalar, Wedge};

/// The configuration R-008 says has 126 offenders, two per layer.
const DIHEDRAL: f64 = 150.0;
const ROTATION: f64 = 17.0;
const SAMPLES: u32 = 65;
const HALF_EXTENT: f64 = 2.0;
const OFFSET_CELLS: f64 = 0.447_213_595_499_957_9;
const BOUNDARY_MARGIN_CELLS: f64 = 1.0;

/// One offending vertex.
struct Offender {
    position: [Scalar; 3],
    angle: Scalar,
    cell: [i64; 3],
    faces: usize,
    /// Worst angle from any incident face to the nearer plane.
    worst_face: Scalar,
    /// Smallest incident face area, in units of `h²`.
    min_area: Scalar,
    /// Largest incident face area, in units of `h²`.
    max_area: Scalar,
    /// Marching Cubes case index of the cell the smallest face came from.
    case: u8,
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let cell_size = 2.0 * HALF_EXTENT / f64::from(SAMPLES - 1);
    let shift = OFFSET_CELLS * cell_size;
    let field = Wedge::new(DIHEDRAL, [-HALF_EXTENT / 2.0 + shift, shift], ROTATION);
    let origin = [-HALF_EXTENT; 3];

    let shape = RuntimeShape3::new([SAMPLES; 3]).expect("the fixture fits u32");
    let mut mesh = MeshBuffer::<Scalar>::new();
    MarchingCubes::<Scalar>::new()
        .extract_into(&field, &shape, origin, cell_size, &mut mesh)
        .expect("extraction");
    recompute(&mut mesh, &field, NormalStrategy::AreaWeightedFaces).expect("normals");

    let mut incident: Vec<Vec<u32>> = vec![Vec::new(); mesh.vertex_count()];
    for (t, tri) in mesh.indices.chunks_exact(3).enumerate() {
        for &v in tri {
            incident[v as usize].push(t as u32);
        }
    }

    let [n0, n1] = field.plane_normals();
    let planes = [field.unrotate(n0), field.unrotate(n1)];
    let limit = HALF_EXTENT - BOUNDARY_MARGIN_CELLS * cell_size;
    let h2 = cell_size * cell_size;

    let mut offenders: Vec<Offender> = Vec::new();
    for (v, faces) in incident.iter().enumerate() {
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
        if angle <= 90.0 {
            continue;
        }

        let mut worst_face = 0.0f64;
        let mut min_area = f64::INFINITY;
        let mut max_area = 0.0f64;
        let mut case = 0u8;
        let mut smallest_cell = [0i64; 3];
        for &t in faces {
            let tri = &mesh.indices[t as usize * 3..t as usize * 3 + 3];
            let p = [
                mesh.positions[tri[0] as usize],
                mesh.positions[tri[1] as usize],
                mesh.positions[tri[2] as usize],
            ];
            let e1 = [p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]];
            let e2 = [p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]];
            let cross = [
                e1[1] * e2[2] - e1[2] * e2[1],
                e1[2] * e2[0] - e1[0] * e2[2],
                e1[0] * e2[1] - e1[1] * e2[0],
            ];
            let len = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
            let area = len / 2.0 / h2;
            if area < min_area {
                min_area = area;
                let centre = [
                    (p[0][0] + p[1][0] + p[2][0]) / 3.0,
                    (p[0][1] + p[1][1] + p[2][1]) / 3.0,
                    (p[0][2] + p[1][2] + p[2][2]) / 3.0,
                ];
                for (axis, slot) in smallest_cell.iter_mut().enumerate() {
                    *slot = ((centre[axis] - origin[axis]) / cell_size).floor() as i64;
                }
                case = 0;
                for corner in 0..8u8 {
                    let c = [
                        origin[0]
                            + cell_size
                                * (smallest_cell[0] as f64 + f64::from(u8::from(corner & 1 != 0))),
                        origin[1]
                            + cell_size
                                * (smallest_cell[1] as f64 + f64::from(u8::from(corner & 2 != 0))),
                        origin[2]
                            + cell_size
                                * (smallest_cell[2] as f64 + f64::from(u8::from(corner & 4 != 0))),
                    ];
                    if field.sample(c) < 0.0 {
                        case |= 1 << corner;
                    }
                }
            }
            max_area = max_area.max(area);
            if len > 0.0 {
                let worst = planes
                    .iter()
                    .map(|n| {
                        let d = (cross[0] * n[0] + cross[1] * n[1] + cross[2] * n[2]) / len;
                        d.clamp(-1.0, 1.0).acos().to_degrees()
                    })
                    .fold(f64::INFINITY, f64::min);
                worst_face = worst_face.max(worst);
            }
        }
        offenders.push(Offender {
            position,
            angle,
            cell: smallest_cell,
            faces: faces.len(),
            worst_face,
            min_area,
            max_area,
            case,
        });
    }

    println!(
        "θ = {DIHEDRAL}°, rotation {ROTATION}°, n = {SAMPLES}: {} offenders over {} layers",
        offenders.len(),
        SAMPLES - 2
    );

    // Distinct cross-section cells. The repeat-per-layer claim, checked.
    let mut columns: Vec<[i64; 2]> = offenders.iter().map(|o| [o.cell[0], o.cell[1]]).collect();
    columns.sort_unstable();
    columns.dedup();
    println!(
        "distinct (x, y) cells: {} → {:?}",
        columns.len(),
        &columns[..columns.len().min(8)]
    );

    let mut fans: Vec<usize> = offenders.iter().map(|o| o.faces).collect();
    fans.sort_unstable();
    fans.dedup();
    println!("incident face counts seen: {fans:?}");

    println!(
        "\n{:>7} {:>26} {:>8} {:>6} {:>10} {:>11} {:>11} {:>6}",
        "angle", "cell", "faces", "case", "worst face", "min area/h²", "max area/h²", "z"
    );
    for o in offenders.iter().take(6) {
        println!(
            "{:>7.2} {:>26} {:>8} {:>6} {:>10.2} {:>11.2e} {:>11.4} {:>6}",
            o.angle,
            format!("({}, {}, {})", o.cell[0], o.cell[1], o.cell[2]),
            o.faces,
            o.case,
            o.worst_face,
            o.min_area,
            o.max_area,
            o.cell[2]
        );
    }

    // **One offender, in full.** The summary above says six faces all lying
    // exactly on a plane and a vertex normal 91° from the gradient, and those
    // two cannot both be true: an area-weighted normal is a positive
    // combination of its faces' normals, so it lies in the cone they span, which
    // for a 150° wedge is 30° wide. One of the two measurements is wrong and
    // this prints enough to say which.
    if let Some(o) = offenders.first() {
        let key = |p: &[Scalar; 3]| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()];
        let v = mesh
            .positions
            .iter()
            .position(|p| key(p) == key(&o.position))
            .expect("the offender came from this buffer");
        println!("\n--- one offender in full: vertex {v} at {:?}", o.position);
        let g = field.gradient(o.position);
        println!("    field value      {:.3e}", field.sample(o.position));
        println!("    gradient         {g:?}");
        println!("    stored normal    {:?}", mesh.normals[v]);
        println!("    plane normals    {:?}  {:?}", planes[0], planes[1]);
        let mut sum = [0.0f64; 3];
        for &t in &incident[v] {
            let tri = &mesh.indices[t as usize * 3..t as usize * 3 + 3];
            let p = [
                mesh.positions[tri[0] as usize],
                mesh.positions[tri[1] as usize],
                mesh.positions[tri[2] as usize],
            ];
            let e1 = [p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]];
            let e2 = [p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]];
            let cross = [
                e1[1] * e2[2] - e1[2] * e2[1],
                e1[2] * e2[0] - e1[0] * e2[2],
                e1[0] * e2[1] - e1[1] * e2[0],
            ];
            for (axis, slot) in sum.iter_mut().enumerate() {
                *slot += cross[axis];
            }
            let len = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
            println!(
                "    face {t:>6}  area/h² {:>9.4}  unit normal [{:>7.4}, {:>7.4}, {:>7.4}]",
                len / 2.0 / h2,
                cross[0] / len,
                cross[1] / len,
                cross[2] / len
            );
        }
        let sl = (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]).sqrt();
        println!(
            "    sum of crosses, normalised: [{:>7.4}, {:>7.4}, {:>7.4}]   |sum|/h² {:.4}",
            sum[0] / sl,
            sum[1] / sl,
            sum[2] / sl,
            sl / 2.0 / h2
        );
        println!(
            "    vertices at this exact position in the buffer: {}",
            mesh.positions
                .iter()
                .filter(|p| key(p) == key(&o.position))
                .count()
        );
    }

    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let dir = root.join("docs/measurements");
    fs::create_dir_all(&dir).expect("create docs/measurements");
    let mut csv = String::from("# R-009: locating the past-90° normals\n");
    let _ = writeln!(
        csv,
        "angle_deg,x,y,z,cell_x,cell_y,cell_z,incident_faces,case,worst_face_deg,\
         min_area_over_h2,max_area_over_h2"
    );
    for o in &offenders {
        let _ = writeln!(
            csv,
            "{:.4},{:.6},{:.6},{:.6},{},{},{},{},{},{:.4},{:.6e},{:.6}",
            o.angle,
            o.position[0],
            o.position[1],
            o.position[2],
            o.cell[0],
            o.cell[1],
            o.cell[2],
            o.faces,
            o.case,
            o.worst_face,
            o.min_area,
            o.max_area
        );
    }
    let path = dir.join("r009-locate.csv");
    fs::write(&path, csv).expect("write csv");
    println!("\n{} rows → {}", offenders.len(), path.display());
}
