//! **Every extractor over a real scanned volume.**
//!
//! Ticket: M-006.
//!
//! ```bash
//! ./scripts/fetch_volumes.sh          # once; the data is not committed
//! cargo bench --bench volumes
//! ```
//!
//! Writes `docs/measurements/volumes.csv`.
//!
//! # Why this exists
//!
//! Every reference field in this crate is **analytic**, which is what makes the
//! accuracy harness exact and is exactly why its numbers cannot be set beside a
//! published isosurfacing benchmark — those are run on CT and simulation data.
//! `docs/measurements/volumes/PROVENANCE.md` carries the sources, the hashes and
//! the reason each file was chosen.
//!
//! There is a second reason and it is the sharper one: **quantised data is what
//! makes Grosso 2017's singular faces reachable.** A continuous `f64` field
//! produces **0 of 299,215** over 400,000 random cells (M-220); `u8` voxels
//! collide readily. Both volumes here are `uint8` deliberately.
//!
//! # Absent data is not a failure
//!
//! The volumes are gitignored, so a clean clone does not have them. This bench
//! prints what to run and exits 0 rather than failing — the alternative is a
//! benchmark suite that cannot be run without a network.
//!
//! # Sign convention
//!
//! The volumes are **densities**: `fuel`'s page says *"the higher the density
//! value, the less presence of air"*, so high is solid. This crate's convention
//! is negative inside, so the field is `f = iso − value` and nothing about the
//! extractors changes.

mod common;

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use isomesh::construct::SampledField;
use isomesh::extractor::Extractor;
use isomesh::validate::{ValidateConfig, validate_indexed};
use isomesh::{MeshBuffer, RuntimeShape3};

/// A volume to read, and where the surface is put.
///
/// The isovalue is the one the dataset page's viewer opens at, which is the
/// closest thing to a canonical choice these files have.
struct Volume {
    file: &'static str,
    size: [u32; 3],
    iso: f64,
}

const VOLUMES: [Volume; 2] = [
    Volume {
        file: "fuel_64x64x64_uint8.raw",
        size: [64, 64, 64],
        iso: 32.0,
    },
    Volume {
        file: "bonsai_256x256x256_uint8.raw",
        size: [256, 256, 256],
        iso: 32.0,
    },
];

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/measurements/volumes")
}

/// Read a `uint8` volume as a field, negative inside.
///
/// The file is one byte per sample, x fastest, little-endian — which for a
/// single byte means nothing, and is why `uint8` is the easy case. The length is
/// checked against the dimensions rather than trusted from the filename.
fn read_u8(path: &Path, size: [u32; 3], iso: f64) -> Result<Vec<f64>, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let want = size[0] as usize * size[1] as usize * size[2] as usize;
    if bytes.len() != want {
        return Err(format!(
            "{}: {} bytes, expected {want} for {size:?}",
            path.display(),
            bytes.len()
        ));
    }
    // `iso - v`, so a dense voxel is negative and the crate's convention holds.
    Ok(bytes.iter().map(|b| iso - f64::from(*b)).collect())
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let dir = dir();
    let missing: Vec<&str> = VOLUMES
        .iter()
        .map(|v| v.file)
        .filter(|f| !dir.join(f).exists())
        .collect();
    if !missing.is_empty() {
        println!(
            "volumes: {} of {} absent — {}",
            missing.len(),
            VOLUMES.len(),
            missing.join(", ")
        );
        println!("run ./scripts/fetch_volumes.sh to fetch them (M-006); skipping, not failing");
        return;
    }

    let mut csv = String::from(
        "volume,samples,algorithm,iso,median_ms,vertices,triangles,\
         non_manifold_edges,non_manifold_vertices,boundary_edges,\
         degenerate_triangles,mean_ratio,irregular_vertices,euler_characteristic\n",
    );

    println!(
        "{:<14} {:<28} {:>9} {:>9} {:>9} {:>6} {:>6} {:>8}",
        "volume", "algorithm", "ms", "verts", "tris", "nm_e", "nm_v", "quality"
    );

    for v in &VOLUMES {
        let values = match read_u8(&dir.join(v.file), v.size, v.iso) {
            Ok(values) => values,
            Err(e) => {
                println!("::error:: {e}");
                std::process::exit(1);
            }
        };
        let shape = match RuntimeShape3::new(v.size) {
            Ok(shape) => shape,
            Err(e) => {
                println!("::error:: {}: {e}", v.file);
                std::process::exit(1);
            }
        };
        // One world unit per voxel. The datasets all declare `spacing 1x1x1`.
        let cell_size = 1.0_f64;
        let field = match SampledField::new(&values, &shape, [0.0; 3], cell_size) {
            Ok(field) => field,
            Err(e) => {
                println!("::error:: {}: {e}", v.file);
                std::process::exit(1);
            }
        };
        let cfg = match ValidateConfig::from_cell_size(cell_size) {
            Ok(cfg) => cfg,
            Err(e) => {
                println!("::error:: {e}");
                std::process::exit(1);
            }
        };
        let short = v.file.split('_').next().unwrap_or(v.file);

        isomesh::for_each_extractor!(f64, |name, extractor| {
            let mut out = MeshBuffer::<f64>::new();
            let mut ms = f64::INFINITY;
            let mut ok = true;
            let mut why = String::new();
            for _ in 0..3 {
                out.reset();
                let t = std::time::Instant::now();
                if let Err(e) =
                    extractor.extract_into(&field, &shape, [0.0; 3], cell_size, &mut out)
                {
                    ok = false;
                    why = format!("{e}");
                    break;
                }
                ms = ms.min(t.elapsed().as_secs_f64() * 1e3);
            }
            if ok && !out.indices.is_empty() {
                let r = validate_indexed(&out.positions, &out.indices, &cfg);
                println!(
                    "{short:<14} {name:<28} {ms:>9.1} {:>9} {:>9} {:>6} {:>6} {:>8.4}",
                    out.positions.len(),
                    out.indices.len() / 3,
                    r.non_manifold_edges,
                    r.non_manifold_vertices,
                    r.mean_ratio
                );
                let _ = writeln!(
                    csv,
                    "{short},{},{name},{},{ms:.6},{},{},{},{},{},{},{:.9},{},{}",
                    v.size[0],
                    v.iso,
                    out.positions.len(),
                    out.indices.len() / 3,
                    r.non_manifold_edges,
                    r.non_manifold_vertices,
                    r.boundary_edges,
                    r.degenerate_triangles,
                    r.mean_ratio,
                    r.irregular_vertices,
                    r.euler_characteristic
                );
            } else {
                println!("{short:<14} {name:<28} {:>9}  {why}", "refused");
            }
        });
    }

    let out = dir.join("../volumes.csv");
    match std::fs::write(&out, &csv) {
        Ok(()) => println!("\nwrote {}", out.display()),
        Err(e) => println!("\n::error:: {}: {e}", out.display()),
    }
}
