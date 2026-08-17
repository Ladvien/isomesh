//! **How much of the residual error can *any* vertex placement reach?**
//!
//! Ticket: R-025. **Exploratory. Nothing is registered against this run** — it
//! exists to check the arithmetic of a hypothesis *before* it is registered,
//! which is now the house style after P-23 clause 3 and P-24 both died on
//! contact and M-313 and M-314 both came out right by measuring first.
//!
//! ```bash
//! cargo bench --bench placement_ceiling
//! ```
//!
//! Writes `docs/measurements/placement_ceiling.csv`.
//!
//! # The question R-025 has to answer before it can be asked
//!
//! R-025 proposes curvature-aware vertex placement and predicts *">20% better
//! Hausdorff on smooth fields"*. That presumes the residual error on a smooth
//! field **is** placement error. It may not be: a flat triangle inscribed in a
//! curved surface deviates from it at its **interior** by the sagitta,
//! `≈ h²/(8R)` for edge `h` on a sphere of radius `R`, and **no vertex placement
//! rule can touch that**. If the sagitta already accounts for the measured
//! error, the hypothesis is dead before any code is written.
//!
//! # Why a sphere, and why no root finding
//!
//! On `sphere` the exact distance from any point to the surface is
//! `| |p − c| − R |` in closed form, so this measures the thing itself rather
//! than a projection that might not converge. Three quantities per resolution:
//!
//! - **vertex error** — how far the placed vertices are from the true surface.
//!   This is what a better placement rule would reduce.
//! - **centroid error** — how far triangle *interiors* are. Chord-versus-arc,
//!   plus whatever the vertices contribute.
//! - **the floor** — vertices projected exactly onto the sphere, then centroids
//!   recomputed. **This is the ceiling on any placement rule**: it is what
//!   remains when placement is perfect.

mod common;

use std::fmt::Write as _;

use isomesh::dual_contouring::DualContouring;
use isomesh::extractor::Extractor;
use isomesh::fields::{ReferenceField, Sphere, Torus};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3};

const RESOLUTIONS: [u32; 4] = [17, 33, 65, 129];

/// Exact signed distance, in closed form, for the two smooth reference fields
/// that publish one. **Both are exact rather than approximate**, which is what
/// lets this measure the thing itself instead of a projection that might not
/// converge.
#[derive(Clone, Copy)]
enum Smooth {
    /// Unit sphere at the origin.
    Sphere,
    /// Major radius 1, minor 0.3 — `Torus::canonical`.
    Torus,
}

impl Smooth {
    fn distance(self, p: [f64; 3]) -> f64 {
        match self {
            Self::Sphere => ((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() - 1.0).abs(),
            Self::Torus => {
                let s = (p[0] * p[0] + p[2] * p[2]).sqrt();
                let q = [s - 1.0, p[1]];
                ((q[0] * q[0] + q[1] * q[1]).sqrt() - 0.3).abs()
            }
        }
    }

    /// The nearest point on the surface. Degenerate only on the axis of
    /// symmetry, which no vertex of either mesh reaches; guarded anyway.
    fn project(self, p: [f64; 3]) -> [f64; 3] {
        match self {
            Self::Sphere => {
                let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                if len <= 0.0 {
                    return p;
                }
                [p[0] / len, p[1] / len, p[2] / len]
            }
            Self::Torus => {
                let s = (p[0] * p[0] + p[2] * p[2]).sqrt();
                if s <= 0.0 {
                    return p;
                }
                // The ring centre nearest `p`, then step `minor` toward `p`.
                let ring = [p[0] / s, 0.0, p[2] / s];
                let q = [p[0] - ring[0], p[1], p[2] - ring[2]];
                let ql = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
                if ql <= 0.0 {
                    return p;
                }
                [
                    ring[0] + q[0] / ql * 0.3,
                    q[1] / ql * 0.3,
                    ring[2] + q[2] / ql * 0.3,
                ]
            }
        }
    }
}

struct Stats {
    max: f64,
    mean: f64,
}

fn stats(values: impl Iterator<Item = f64>) -> Stats {
    let (mut max, mut sum, mut n) = (0.0_f64, 0.0_f64, 0u64);
    for v in values {
        max = max.max(v);
        sum += v;
        n += 1;
    }
    Stats {
        max,
        mean: if n == 0 { 0.0 } else { sum / n as f64 },
    }
}

fn centroids(positions: &[[f64; 3]], indices: &[u32]) -> Vec<[f64; 3]> {
    let mut out = Vec::with_capacity(indices.len() / 3);
    for t in indices.chunks_exact(3) {
        let (Some(a), Some(b), Some(c)) = (
            positions.get(t[0] as usize),
            positions.get(t[1] as usize),
            positions.get(t[2] as usize),
        ) else {
            continue;
        };
        out.push([
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        ]);
    }
    out
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!(
        "{:<8} {:<18} {:>5} {:>11} {:>11} {:>11} {:>11} {:>12}",
        "field", "extractor", "n", "vertex max", "centroid mx", "hausdorff", "FLOOR", "ceiling %"
    );
    let mut csv = String::from(
        "field,extractor,samples_per_axis,cell_size,vertices,triangles,vertex_max,vertex_mean,\
         centroid_max,centroid_mean,hausdorff,floor_max,floor_mean,sagitta_estimate,ceiling_share\n",
    );

    for (field_name, kind) in [("sphere", Smooth::Sphere), ("torus", Smooth::Torus)] {
        for n in RESOLUTIONS {
            let Ok(shape) = RuntimeShape3::new([n; 3]) else {
                continue;
            };
            let (lo, hi) = match kind {
                Smooth::Sphere => Sphere::<f64>::canonical().domain(),
                Smooth::Torus => Torus::<f64>::canonical().domain(),
            };
            let h = (hi[0] - lo[0]) / f64::from(n - 1);

            for (name, mesh) in [
                ("dual_contouring", {
                    let mut out = MeshBuffer::<f64>::new();
                    let ok = match kind {
                        Smooth::Sphere => DualContouring::<f64>::new()
                            .extract_into(&Sphere::<f64>::canonical(), &shape, lo, h, &mut out)
                            .is_ok(),
                        Smooth::Torus => DualContouring::<f64>::new()
                            .extract_into(&Torus::<f64>::canonical(), &shape, lo, h, &mut out)
                            .is_ok(),
                    };
                    if ok { Some(out) } else { None }
                }),
                ("marching_cubes", {
                    let mut out = MeshBuffer::<f64>::new();
                    let ok = match kind {
                        Smooth::Sphere => MarchingCubes::<f64>::new()
                            .extract_into(&Sphere::<f64>::canonical(), &shape, lo, h, &mut out)
                            .is_ok(),
                        Smooth::Torus => MarchingCubes::<f64>::new()
                            .extract_into(&Torus::<f64>::canonical(), &shape, lo, h, &mut out)
                            .is_ok(),
                    };
                    if ok { Some(out) } else { None }
                }),
            ] {
                let Some(mesh) = mesh else { continue };
                if mesh.indices.is_empty() {
                    continue;
                }

                let vertex = stats(mesh.positions.iter().map(|p| kind.distance(*p)));
                let centroid = stats(
                    centroids(&mesh.positions, &mesh.indices)
                        .iter()
                        .map(|p| kind.distance(*p)),
                );

                // The floor: place every vertex perfectly, then look at what the
                // triangles still do. No placement rule can beat this.
                let perfect: Vec<[f64; 3]> =
                    mesh.positions.iter().map(|p| kind.project(*p)).collect();
                let floor = stats(
                    centroids(&perfect, &mesh.indices)
                        .iter()
                        .map(|p| kind.distance(*p)),
                );

                // Chord-versus-arc for an edge of about one cell on the unit sphere.
                let sagitta = h * h / 8.0;
                // **The ceiling on any placement rule.** The reported Hausdorff is
                // the worst mesh sample, and the harness samples vertices *and*
                // centroids — so it is `max(vertex, centroid)`. Place every vertex
                // perfectly and the vertex term goes to zero, leaving the centroid
                // term computed on perfect vertices. That is the floor, and the
                // difference is the most any placement rule can win.
                let current = vertex.max.max(centroid.max);
                let share = if current > 0.0 {
                    (current - floor.max) / current
                } else {
                    0.0
                };

                println!(
                    "{field_name:<8} {name:<18} {n:>5} {:>11.4e} {:>11.4e} {:>11.4e} {:>11.4e} {:>11.1}%",
                    vertex.max,
                    centroid.max,
                    current,
                    floor.max,
                    share * 100.0
                );
                let _ = writeln!(
                    csv,
                    "{field_name},{name},{n},{h:.9},{},{},{:.9e},{:.9e},{:.9e},{:.9e},{current:.9e},{:.9e},{:.9e},{sagitta:.9e},{share:.6}",
                    mesh.positions.len(),
                    mesh.indices.len() / 3,
                    vertex.max,
                    vertex.mean,
                    centroid.max,
                    centroid.mean,
                    floor.max,
                    floor.mean
                );
            }
        }
    }

    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements/placement_ceiling.csv");
    match std::fs::write(&path, &csv) {
        Ok(()) => println!("\nwrote {}", path.display()),
        Err(e) => println!("\n::error:: {}: {e}", path.display()),
    }
}
