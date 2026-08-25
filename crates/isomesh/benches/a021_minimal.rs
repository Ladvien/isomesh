//! **A-021 — the smallest sign configuration that makes a dual method emit an
//! edge with three or more faces.**
//!
//! ```bash
//! cargo bench --bench a021_minimal
//! ```
//!
//! Writes `docs/measurements/a021_minimal.csv`.
//!
//! # Minimised from a real occurrence, not guessed
//!
//! Rule 5 forbids inventing a configuration, and the ticket says so explicitly.
//! Reasoning from `DualMesher::emit_quads` gives a *plausible* answer — the mesh
//! edge between two face-adjacent cells is a side of the quad of **every**
//! crossed boundary edge of their shared face, and a face with a sign change has
//! two of those — which would make every dual mesh non-manifold everywhere. It
//! measurably does not. So the reasoning is wrong somewhere, and the way to find
//! out is to take a defect that exists and shrink it.
//!
//! # Block minimisation was tried first and did not answer the question
//!
//! Lifting the 4³ sample block around each of P-14's 314 defects and re-meshing
//! it alone reproduces **314 of 314**, shrinks to 3³ for **none**, and yields
//! **298 distinct sign patterns** out of 314 — the surrounding signs are mostly
//! irrelevant and the method cannot tell which. A census beats a minimisation
//! when the defect is a local *predicate* rather than a local *pattern*.
//!
//! So this counts, for every non-manifold edge, how many boundary edges of the
//! grid face its two cells share carry a sign change — **with the same census
//! over the manifold edges as a control**, because a number that is only ever
//! measured on the defect cannot say the defect is unusual.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use isomesh::extractor::Extractor;
use isomesh::fields::ReferenceField;
use isomesh::surface_nets::SurfaceNets;
use isomesh::validate::{ValidateConfig, validate_features};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf, Shape3};

/// Samples per axis for the field the defects are drawn from.
const SAMPLES: u32 = 49;

/// Which cell a dual vertex sits in.
fn cell_of(p: [f64; 3], origin: [f64; 3], h: f64) -> [i64; 3] {
    [
        ((p[0] - origin[0]) / h).floor() as i64,
        ((p[1] - origin[1]) / h).floor() as i64,
        ((p[2] - origin[2]) / h).floor() as i64,
    ]
}

/// Mesh a grid of values with Surface Nets.
fn mesh_values(values: &[f64], shape: &RuntimeShape3, h: f64) -> MeshBuffer<f64> {
    let mut out = MeshBuffer::<f64>::new();
    let field = isomesh::construct::SampledField::new(values, shape, [0.0; 3], h).expect("wrap");
    SurfaceNets::<f64>::new()
        .extract_into(&field, shape, [0.0; 3], h, &mut out)
        .expect("extraction");
    out
}

/// **A controlled 4³ case, printed rather than reasoned about.**
///
/// Reading `emit_quads` says the mesh edge between two face-adjacent cells is a
/// side of the quad of every crossed boundary edge of their shared face, and
/// that a face carrying a sign change has two of those — which would make every
/// dual mesh non-manifold everywhere, and measurably does not. Rather than
/// argue with the measurement, this prints the face count of every mesh edge on
/// a grid whose signs are a plain half-space, where nothing subtle is going on.
fn probe_plane(h: f64) {
    let n = 4u32;
    let shape = RuntimeShape3::new([n; 3]).expect("valid shape");
    // A half-space: negative for y below the middle. Nothing ambiguous anywhere.
    let mut values = Vec::new();
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let _ = (x, z);
                values.push(f64::from(y) - 1.5);
            }
        }
    }
    let out = mesh_values(&values, &shape, h);
    let mut faces: BTreeMap<[u32; 2], u32> = BTreeMap::new();
    for t in out.indices.as_chunks::<3>().0 {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let key = if a < b { [a, b] } else { [b, a] };
            *faces.entry(key).or_insert(0) += 1;
        }
    }
    let mut hist: BTreeMap<u32, u32> = BTreeMap::new();
    for c in faces.values() {
        *hist.entry(*c).or_insert(0) += 1;
    }
    println!(
        "\nprobe: a plain half-space at 4³ — {} vertices, {} triangles",
        out.positions.len(),
        out.indices.len() / 3
    );
    println!("  mesh edges by incident-face count: {hist:?}");
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let field = isomesh::fields::noise_cavity::<f64>();
    probe_plane(0.25);
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(SAMPLES - 1);
    let shape = RuntimeShape3::new([SAMPLES; 3]).expect("valid shape");
    let size = shape.size();

    let mut full = Vec::with_capacity(shape.element_count());
    for z in 0..SAMPLES {
        for y in 0..SAMPLES {
            for x in 0..SAMPLES {
                full.push(field.sample([
                    lo[0] + h * f64::from(x),
                    lo[1] + h * f64::from(y),
                    lo[2] + h * f64::from(z),
                ]));
            }
        }
    }

    let mut out = MeshBuffer::<f64>::new();
    SurfaceNets::<f64>::new()
        .extract_into(&field, &shape, lo, h, &mut out)
        .expect("extraction");
    let cfg = ValidateConfig::from_cell_size(h).expect("valid cell size");
    let (report, features) = validate_features(&out.positions, &out.indices, &cfg);

    println!(
        "noise_cavity at {SAMPLES}³ under surface_nets: {} non-manifold edges",
        report.non_manifold_edges
    );

    // **The census.** For each non-manifold edge, find the grid face its two
    // cells share and count how many of that face's four boundary edges are
    // crossed. `emit_quads` puts a quad side in exactly one triangle, so a mesh
    // edge carries one face per crossed boundary edge — 2 is manifold and 4 is
    // not. 4 crossed boundary edges on a square means the corner signs
    // alternate, which is a **sign-ambiguous face**.
    let sign_at = |x: i64, y: i64, z: i64| -> bool {
        full[((z * i64::from(size[1]) + y) * i64::from(size[0]) + x) as usize] < 0.0
    };
    let mut crossed_hist: BTreeMap<usize, u32> = BTreeMap::new();
    let mut not_face_adjacent = 0u32;
    for e in &features.edges {
        let a = cell_of(out.positions[e[0] as usize], lo, h);
        let b = cell_of(out.positions[e[1] as usize], lo, h);
        let delta = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let moved: Vec<usize> = (0..3).filter(|&i| delta[i] != 0).collect();
        if moved.len() != 1 || delta[moved[0]].abs() != 1 {
            // A quad *diagonal* joins cells differing on two axes. It sits in
            // two triangles of its own quad and appears in no other, so it can
            // only be non-manifold by some other route -- counted rather than
            // assumed away.
            not_face_adjacent += 1;
            continue;
        }
        // The shared face: fixed on the moved axis at the higher cell's
        // coordinate, spanning one cell on each of the other two.
        let n = moved[0];
        let (p, q) = ((n + 1) % 3, (n + 2) % 3);
        let lo_cell = if delta[n] > 0 { a } else { b };
        let mut corner = [0i64; 3];
        corner[n] = lo_cell[n] + 1;
        corner[p] = lo_cell[p];
        corner[q] = lo_cell[q];
        // The face's four corners, and its four boundary edges as corner pairs.
        let at = |dp: i64, dq: i64| {
            let mut c = corner;
            c[p] += dp;
            c[q] += dq;
            sign_at(c[0], c[1], c[2])
        };
        let s = [at(0, 0), at(1, 0), at(1, 1), at(0, 1)];
        let crossed = (0..4).filter(|&i| s[i] != s[(i + 1) % 4]).count();
        *crossed_hist.entry(crossed).or_insert(0) += 1;
    }
    println!(
        "\ncrossed boundary edges on the shared face, per non-manifold edge: {crossed_hist:?}"
    );
    println!("  non-manifold edges whose cells are not face-adjacent: {not_face_adjacent}");

    // The control: the same census over **every** mesh edge, so the number above
    // means something. An instrument that only ever sees the defect cannot say
    // the defect is unusual.
    let mut all_faces: BTreeMap<[u32; 2], u32> = BTreeMap::new();
    for t in out.indices.as_chunks::<3>().0 {
        for (x, y) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let key = if x < y { [x, y] } else { [y, x] };
            *all_faces.entry(key).or_insert(0) += 1;
        }
    }
    let mut control: BTreeMap<usize, u32> = BTreeMap::new();
    for (edge, count) in &all_faces {
        if *count != 2 {
            continue;
        }
        let a = cell_of(out.positions[edge[0] as usize], lo, h);
        let b = cell_of(out.positions[edge[1] as usize], lo, h);
        let delta = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let moved: Vec<usize> = (0..3).filter(|&i| delta[i] != 0).collect();
        if moved.len() != 1 || delta[moved[0]].abs() != 1 {
            continue;
        }
        let n = moved[0];
        let (p, q) = ((n + 1) % 3, (n + 2) % 3);
        let lo_cell = if delta[n] > 0 { a } else { b };
        let mut corner = [0i64; 3];
        corner[n] = lo_cell[n] + 1;
        corner[p] = lo_cell[p];
        corner[q] = lo_cell[q];
        let at = |dp: i64, dq: i64| {
            let mut c = corner;
            c[p] += dp;
            c[q] += dq;
            sign_at(c[0], c[1], c[2])
        };
        let s = [at(0, 0), at(1, 0), at(1, 1), at(0, 1)];
        let crossed = (0..4).filter(|&i| s[i] != s[(i + 1) % 4]).count();
        *control.entry(crossed).or_insert(0) += 1;
    }
    println!("  control — the same census over manifold (2-face) edges: {control:?}");

    let mut csv = String::from("population,crossed_boundary_edges,count\n");
    for (k, v) in &crossed_hist {
        let _ = writeln!(csv, "non_manifold,{k},{v}");
    }
    for (k, v) in &control {
        let _ = writeln!(csv, "manifold,{k},{v}");
    }
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements")
        .join("a021_minimal.csv");
    if let Some(dir) = path.parent() {
        let _ = fs::create_dir_all(dir);
    }
    let _ = fs::write(&path, csv);
    println!("\nwrote {}", path.display());
}
