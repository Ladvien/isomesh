//! M-001 — every extractor, one process, one run.
//!
//! ```bash
//! cargo bench --bench family
//! ```
//!
//! Writes `docs/measurements/family.csv`.
//!
//! # The ticket nineteen references pointed at, and which had no row
//!
//! `BACKLOG_ARCHIVE.md` says three separate times that a deferred re-measurement
//! *"belongs to M-001"*, and `FINDINGS.md` names it ten more. Until R-005 went
//! looking, M-001 **had no ticket row in either file** — so the designated home
//! for every deferred re-measurement did not exist and none of them could ever
//! be taken off the queue.
//!
//! # Why it is a new bench rather than a wider `resolution_sweep`
//!
//! `resolution_sweep` answers T-006's question — does `t = a + b·n³` describe
//! this path, and what is `a` — over three algorithms, and its CSV is quoted by
//! ✗14, M-19, M-20, M-21, M-22 and O-11. Its own header says a fourth row would
//! rewrite that evidence. This is a different question with a different answer
//! shape (every extractor, one number each, side by side), so it gets its own
//! file and leaves the fit sweep and both of its committed CSVs untouched.
//!
//! # Cycles, not milliseconds
//!
//! M-280: two runs of one binary reported the same measurement as 8.13 and
//! 14.66 ns/sample while cycles/sample held at ~34, because this host's
//! governor spans 1.96–5.62 GHz. So every row carries `cycles_per_sample` and
//! the `ghz` it was taken at, and the milliseconds are there for continuity with
//! what is already written down rather than as the number to compare.
//!
//! Counters are `perf_event_open` and therefore Linux. On any other platform
//! those three columns read **`unavailable`** — a word, not a blank and not a
//! zero, so a reader cannot mistake a missing instrument for a measurement.
//!
//! # No silent truncation
//!
//! Subgrid Marching Tetrahedra is **70× classic Marching Tetrahedra** (M-98),
//! which at 256³ is hours rather than seconds. Rather than hard-code which
//! extractor stops where — a list that drifts from the thing it names — each one
//! climbs until a run exceeds [`BUDGET_MS`] and then stops, and the resolution
//! it reached is printed and recorded. A row that is absent is absent because a
//! stated budget said so.

mod common;

use std::fmt::Write as _;
use std::fs;
use std::hint::black_box;
use std::path::PathBuf;
use std::time::Instant;

use isomesh::extractor::Extractor;
use isomesh::fields::Sphere;
use isomesh::{MeshBuffer, for_each_extractor};

/// `f32`, matching `resolution_sweep`, so the two can be compared at all.
type Scalar = f32;

/// Samples per axis. The resolution set `resolution_sweep` already uses.
const RESOLUTIONS: [u32; 9] = [16, 24, 32, 48, 64, 96, 128, 192, 256];

/// Untimed runs first, so this is steady-state re-meshing rather than
/// first-touch page faults on freshly grown scratch.
const WARMUP_RUNS: u32 = 2;

/// Timed runs per resolution; the median is taken.
const TIMED_RUNS: usize = 5;

/// An extractor stops climbing once one run costs more than this.
///
/// Two seconds is chosen so the whole sweep is minutes rather than hours and so
/// that the three algorithms `resolution_sweep` already covers still reach 256³
/// on the machines this runs on. It is a budget, not a judgement: an entry that
/// stops early has its last resolution recorded rather than being quietly
/// shorter than its neighbours.
const BUDGET_MS: f64 = 2000.0;

/// One measured configuration.
struct Row {
    algorithm: &'static str,
    samples: u32,
    median_ms: f64,
    ns_per_sample: f64,
    /// `None` where there are no hardware counters.
    cycles_per_sample: Option<f64>,
    instructions_per_sample: Option<f64>,
    vertices: usize,
    triangles: usize,
}

impl Row {
    fn ipc(&self) -> Option<f64> {
        match (self.instructions_per_sample, self.cycles_per_sample) {
            (Some(i), Some(c)) if c > 0.0 => Some(i / c),
            _ => None,
        }
    }

    fn ghz(&self) -> Option<f64> {
        self.cycles_per_sample.map(|c| c / self.ns_per_sample)
    }
}

/// A cell that is a number where it was measured and the word `unavailable`
/// where the instrument does not exist.
fn cell(value: Option<f64>, places: usize) -> String {
    value.map_or_else(|| String::from("unavailable"), |v| format!("{v:.places$}"))
}

/// Time and count one extraction, `TIMED_RUNS` times, and keep the median.
fn measure<E: Extractor<Scalar>>(extractor: &mut E, samples: u32) -> Row {
    let field = Sphere::<Scalar>::canonical();
    let (shape, origin, cell_size) = common::grid(&field, samples);
    let mut mesh = MeshBuffer::<Scalar>::new();

    for _ in 0..WARMUP_RUNS {
        // The caller owns the buffer and resets it (rule 6). Leaving this out
        // makes later runs pay reallocation the extraction did not cause, and
        // the tell is a triangle count that is not monotone in `n` — M-279 lost
        // a run to exactly that.
        mesh.reset();
        extractor
            .extract_into(&field, &shape, origin, cell_size, &mut mesh)
            .expect("extraction");
        black_box(&mesh);
    }

    #[cfg(target_os = "linux")]
    let mut probe = common::counters::Probe::open();
    // `(nanos, cycles, instructions)`, sorted by the first.
    let mut runs: Vec<(u128, Option<u64>, Option<u64>)> = Vec::with_capacity(TIMED_RUNS);
    for _ in 0..TIMED_RUNS {
        mesh.reset();
        #[cfg(target_os = "linux")]
        probe.reset_and_enable();
        let started = Instant::now();
        extractor
            .extract_into(&field, &shape, origin, cell_size, &mut mesh)
            .expect("extraction");
        let nanos = started.elapsed().as_nanos();
        #[cfg(target_os = "linux")]
        {
            probe.disable();
            let counts = probe.read();
            assert!(
                counts.worst_ratio() >= common::counters::MIN_TIME_RATIO,
                "a counter ran only {:.1}% of the time it was enabled, so its value is a scaled \
                 estimate rather than a measurement",
                counts.worst_ratio() * 100.0
            );
            runs.push((
                nanos,
                Some(counts.cycles.count),
                Some(counts.instructions.count),
            ));
        }
        #[cfg(not(target_os = "linux"))]
        runs.push((nanos, None, None));
        black_box(&mesh);
    }
    runs.sort_unstable();
    let (nanos, cycles, instructions) = runs[TIMED_RUNS / 2];

    let total = f64::from(samples).powi(3);
    Row {
        algorithm: "",
        samples,
        median_ms: nanos as f64 / 1e6,
        ns_per_sample: nanos as f64 / total,
        cycles_per_sample: cycles.map(|c| c as f64 / total),
        instructions_per_sample: instructions.map(|i| i as f64 / total),
        vertices: mesh.vertex_count(),
        triangles: mesh.triangle_count(),
    }
}

/// Short output of a command, or `"unknown"`.
fn ask(program: &str, args: &[&str]) -> String {
    std::process::Command::new(program)
        .args(args)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| String::from("unknown"))
}

fn write_csv(rows: &[Row]) -> PathBuf {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let dir = root.join("docs/measurements");
    fs::create_dir_all(&dir).expect("create docs/measurements");
    let path = dir.join("family.csv");

    let machine = std::process::Command::new(root.join("scripts/machine.sh"))
        .arg("--slug")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| String::from("unknown"));
    let dirty = ask("git", &["status", "--porcelain"]);

    let mut csv = String::new();
    // `scripts/regress.sh` skips `#` lines, so provenance costs nothing and a
    // file from the wrong machine stops being invisible — which is the hazard
    // `resolution_sweep.csv` still carries.
    let _ = writeln!(csv, "# family sweep (M-001), sphere, f32");
    let _ = writeln!(
        csv,
        "# commit {}{} on {machine} at {}",
        ask("git", &["rev-parse", "--short", "HEAD"]),
        if dirty == "unknown" || dirty.is_empty() {
            ""
        } else {
            " (WORKING TREE DIRTY)"
        },
        ask("date", &["-u", "+%Y-%m-%dT%H:%M:%SZ"])
    );
    let _ = writeln!(
        csv,
        "# cycles/ipc/ghz are `unavailable` where perf_event_open is not; \
         budget {BUDGET_MS:.0} ms per run"
    );
    let _ = writeln!(
        csv,
        "algorithm,scalar,field,samples,n_cubed,median_ms,ns_per_sample,cycles_per_sample,\
         instructions_per_sample,ipc,ghz,vertices,triangles"
    );
    for r in rows {
        let _ = writeln!(
            csv,
            "{},f32,sphere,{},{:.0},{:.6},{:.4},{},{},{},{},{},{}",
            r.algorithm,
            r.samples,
            f64::from(r.samples).powi(3),
            r.median_ms,
            r.ns_per_sample,
            cell(r.cycles_per_sample, 3),
            cell(r.instructions_per_sample, 3),
            cell(r.ipc(), 3),
            cell(r.ghz(), 3),
            r.vertices,
            r.triangles
        );
    }
    fs::write(&path, csv).expect("write csv");
    path.canonicalize().unwrap_or(path)
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!(
        "{:<28} {:>5} {:>11} {:>10} {:>10} {:>7} {:>6} {:>10}",
        "algorithm", "n", "median_ms", "ns/sample", "cyc/sample", "IPC", "GHz", "triangles"
    );

    let mut rows: Vec<Row> = Vec::new();
    let mut capped: Vec<(&'static str, u32)> = Vec::new();
    // **Resolution outermost, and the extractor rebuilt inside it.** Both halves
    // of that are load-bearing and both were measured before being chosen
    // (M-281). A single extractor carried across resolutions arrives at 256³
    // with buffers grown from 192³, and Marching Cubes then measures **131.6 ms
    // against `resolution_sweep`'s 152.7** on the same machine at the same
    // clock — a 16% disagreement between two of this repo's own benches, which
    // is the thing M-001 exists to remove. Rebuilding per resolution is what the
    // sweep does, so the two now answer the same question the same way.
    // Interleaving the extractors also means a clock excursion hits all of them
    // rather than whichever one happened to be running.
    for samples in RESOLUTIONS {
        for_each_extractor!(f32, |name, extractor| {
            if !capped.iter().any(|(n, _)| *n == name) {
                let mut row = measure(&mut extractor, samples);
                row.algorithm = name;
                println!(
                    "{name:<28} {samples:>5} {:>11.3} {:>10.4} {:>10} {:>7} {:>6} {:>10}",
                    row.median_ms,
                    row.ns_per_sample,
                    cell(row.cycles_per_sample, 2),
                    cell(row.ipc(), 2),
                    cell(row.ghz(), 2),
                    row.triangles
                );
                if row.median_ms > BUDGET_MS {
                    capped.push((name, samples));
                }
                rows.push(row);
            }
        });
    }
    rows.sort_by_key(|r| (r.algorithm, r.samples));

    let path = write_csv(&rows);
    println!("\n{} rows → {}", rows.len(), path.display());
    if capped.is_empty() {
        println!("no entry hit the {BUDGET_MS:.0} ms budget; every one reached 256³");
    } else {
        println!("stopped at the {BUDGET_MS:.0} ms budget, last resolution measured:");
        for (name, samples) in &capped {
            println!("  {name:<28} {samples}³");
        }
    }
}
