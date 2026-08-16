//! **The shootout for the input half of the pipeline.**
//!
//! Ticket: T-018. `benches/shootout.rs` compares extractors — field in,
//! triangles out. This compares *constructors* — something in, field out — and
//! until S-001 there was nothing to compare.
//!
//! ```bash
//! cargo bench --bench constructors
//! ```
//!
//! Writes `docs/measurements/constructors.csv`.
//!
//! # Three families, and mixing them would be the mistake
//!
//! **Signs → field.** S-001 exact transform, S-002 sweeping, S-003 marching.
//! All three read the input's *signs* and its sub-cell crossings and discard its
//! magnitudes, so handing them an already-exact field gives them no advantage
//! and their error is a clean number.
//!
//! **Field → field.** S-004 narrow-band reinitialisation. It **keeps** the input
//! outside the band, so feeding it the analytic field makes most of its output
//! ground truth by construction — the first version of this bench did exactly
//! that and ranked it best on both speed and accuracy, which was a measurement
//! of the harness (M-270). It is fed a **degraded** field here: the analytic one
//! scaled by two, which has the same zero set and twice the gradient, the shape
//! a CSG scale or a careless edit produces. Its row is reported beside the
//! degraded input's *own* error, so what the column shows is what
//! reinitialisation bought rather than what it inherited.
//!
//! **Mesh → field.** S-006 pseudonormal, S-007 winding number. These take
//! triangles, so measuring them against an analytic field necessarily includes
//! the **meshing** that produced those triangles. Their rows are reported
//! separately and their error is a round-trip error, not a constructor error.
//! Putting all six in one ranked table would make the mesh-based pair look worse
//! than they are for a reason that has nothing to do with them.
//!
//! **S-005, jump flooding, is absent from this process and not from the
//! comparison.** It lives in `isomesh-gpu`, and rule 2 keeps that dependency out
//! of this crate — so its numbers are produced by its own test and committed to
//! `docs/measurements/jump_flood.csv`. The recommendation below quotes them.
//!
//! # Memory: what is here, and the part that is not
//!
//! **`out_kib` is exact** — the output buffer's size, which is the allocation a
//! consumer has to budget for and hold. It is computed, not sampled, because
//! `len × size_of::<R>()` has nothing to estimate.
//!
//! **Peak working set is not measured, and is not estimated either.** The
//! obvious instrument is a counting global allocator; that needs
//! `unsafe impl GlobalAlloc`, and this workspace sets `unsafe_code = "forbid"`,
//! which is the basis of the crate's "100% safe Rust" claim (M-147). Writing a
//! working-set figure derived from reading the algorithm would be a performance
//! number with no benchmark behind it — rule 4 — and it would have been wrong
//! anyway: the marching constructor's `BTreeSet` is data-dependent and can
//! exceed its own grid.
//!
//! **T-019 fills the gap from outside.** `--only <constructor>` runs exactly one
//! constructor on one field at one resolution and nothing else, so the
//! process's own peak resident set attributes to it;
//! `scripts/constructor_memory.sh` drives one process per constructor under the
//! platform's `time` tool and subtracts a `baseline` run that builds the input
//! and stops. The instrument is then the operating system rather than a
//! `GlobalAlloc`, and the crate stays free of `unsafe` rather than the rule
//! being bent.

use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use isomesh::MeshBuffer;
use isomesh::construct::{
    from_mesh::signed_distance_from_mesh, reinitialise_narrow_band, signed_distance_field,
    signed_distance_field_marched, signed_distance_field_swept,
    winding::signed_distance_from_mesh_winding,
};
use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{RuntimeShape3, Sdf, Shape3};

/// Grid sizes. `O(n log n)` against `O(n)` only separates with size.
const RESOLUTIONS: [u32; 2] = [33, 65];

/// Timed runs; the best is taken, since the fastest run is the one with the
/// least interference.
const RUNS: u32 = 3;

/// Band half-width, in cells, for the narrow-band reinitialiser.
const BAND: u32 = 3;

struct Row {
    field: &'static str,
    constructor: &'static str,
    family: &'static str,
    samples: u32,
    ms: f64,
    out_kib: f64,
    worst: f64,
    worst_near: f64,
}

/// The process's peak resident set, in KiB, or `None` where the kernel does not
/// publish one.
///
/// `VmHWM` from `/proc/self/status` — a **file read**, which is ordinary safe
/// Rust. The instrument that was ruled out is a counting `GlobalAlloc`, because
/// that needs `unsafe impl` and the workspace forbids it; nothing about the
/// measurement itself required going outside the crate. What it does require is
/// going outside the *process*, since `VmHWM` is a high-water mark over a
/// process's whole life and would otherwise attribute every constructor's peak
/// to whichever ran first — hence `--only`.
///
/// Linux only. macOS has no `/proc`, and its `ru_maxrss` needs `libc`, which
/// this crate does not depend on. CI runs Linux, which is where the number is
/// checked, so the figure exists exactly where it is used.
fn peak_rss_kib() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    status
        .lines()
        .find(|l| l.starts_with("VmHWM:"))
        .and_then(|l| l.split_whitespace().nth(1))
        .and_then(|v| v.parse().ok())
}

/// Value of `--<flag>`, if present.
fn flag(name: &str) -> Option<String> {
    let args: Vec<String> = std::env::args().collect();
    let i = args.iter().position(|a| a == name)?;
    args.get(i + 1).cloned()
}

/// One constructor, one field, one resolution, and nothing else in the process.
///
/// The whole point is what is *absent*: no other constructor runs, no CSV is
/// written, and nothing is printed but a single line, so the peak resident set
/// the caller reads belongs to this constructor plus the input it needed.
/// `baseline` and `mesh-baseline` build that input and stop, which is what makes
/// the subtraction meaningful.
fn run_one(which: &str) {
    let field_name = flag("--field").unwrap_or_else(|| String::from("sphere"));
    let samples: u32 = flag("--samples").and_then(|s| s.parse().ok()).unwrap_or(65);

    isomesh::for_each_reference_field!(f64, |name, field| {
        if name == field_name {
            let (lo, hi) = field.domain();
            let h = (hi[0] - lo[0]) / f64::from(samples - 1);
            let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
            let size = shape.size();
            let mut grid = Vec::with_capacity(shape.element_count());
            for z in 0..size[2] {
                for y in 0..size[1] {
                    for x in 0..size[0] {
                        grid.push(field.sample([
                            lo[0] + h * f64::from(x),
                            lo[1] + h * f64::from(y),
                            lo[2] + h * f64::from(z),
                        ]));
                    }
                }
            }

            let built: Vec<f64> = match which {
                "baseline" => Vec::new(),
                "exact" => signed_distance_field(&grid, &shape, h).expect("exact"),
                "swept" => signed_distance_field_swept(&grid, &shape, h).expect("swept"),
                "marched" => signed_distance_field_marched(&grid, &shape, h).expect("marched"),
                "band" => {
                    let degraded: Vec<f64> = grid.iter().map(|v| v * 2.0).collect();
                    reinitialise_narrow_band(&degraded, &shape, h, BAND)
                        .expect("band")
                        .0
                }
                "mesh-baseline" | "pseudonormal" | "winding" => {
                    let mesh = mesh_of(&grid, &shape, lo, h);
                    match which {
                        "pseudonormal" => {
                            signed_distance_from_mesh(&mesh.positions, &mesh.indices, &shape, lo, h)
                                .expect("pseudonormal")
                        }
                        "winding" => signed_distance_from_mesh_winding(
                            &mesh.positions,
                            &mesh.indices,
                            &shape,
                            lo,
                            h,
                            0.5,
                        )
                        .expect("winding"),
                        _ => {
                            // `mesh-baseline`: the mesh is the thing being paid
                            // for, so it has to survive to the end of the
                            // process or the peak would not include it.
                            std::hint::black_box(&mesh);
                            Vec::new()
                        }
                    }
                }
                other => {
                    eprintln!("unknown constructor {other}");
                    std::process::exit(2);
                }
            };
            std::hint::black_box(&built);
            match peak_rss_kib() {
                Some(kib) => println!("{which} {field_name} {samples} {} {kib}", built.len()),
                None => println!("{which} {field_name} {samples} {} unavailable", built.len()),
            }
        }
    });
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    if let Some(which) = flag("--only") {
        run_one(&which);
        return;
    }

    println!("constructors — the input half of the pipeline\n");
    println!(
        "{:<16} {:<12} {:>7} {:>10} {:>11} {:>11} {:>12}",
        "field", "constructor", "samples", "ms", "out KiB", "worst", "worst_near"
    );

    let mut rows: Vec<Row> = Vec::new();

    isomesh::for_each_reference_field!(f64, |name, field| {
        // Inline block, so no `return` in here (M-253).
        //
        // Error is measured against `field.sample`, which is only a distance
        // where `FieldBound::Exact` says so. A gyroid row would be a column of
        // nonsense, so it is skipped rather than annotated.
        if field.bound().is_exact() {
            for samples in RESOLUTIONS {
                let (lo, hi) = field.domain();
                let h = (hi[0] - lo[0]) / f64::from(samples - 1);
                let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
                let size = shape.size();

                let mut grid = Vec::with_capacity(shape.element_count());
                for z in 0..size[2] {
                    for y in 0..size[1] {
                        for x in 0..size[0] {
                            grid.push(field.sample([
                                lo[0] + h * f64::from(x),
                                lo[1] + h * f64::from(y),
                                lo[2] + h * f64::from(z),
                            ]));
                        }
                    }
                }

                let truth = |x: u32, y: u32, z: u32| {
                    field.sample([
                        lo[0] + h * f64::from(x),
                        lo[1] + h * f64::from(y),
                        lo[2] + h * f64::from(z),
                    ])
                };

                // Signs → field.
                let field_to_field: [(&str, &dyn Fn() -> Vec<f64>); 3] = [
                    ("exact", &|| {
                        signed_distance_field(&grid, &shape, h).expect("exact")
                    }),
                    ("swept", &|| {
                        signed_distance_field_swept(&grid, &shape, h).expect("swept")
                    }),
                    ("marched", &|| {
                        signed_distance_field_marched(&grid, &shape, h).expect("marched")
                    }),
                ];
                for (label, run) in field_to_field {
                    let (ms, out_kib, built) = timed(run);
                    let (worst, worst_near) = error(&built, size, h, &truth);
                    emit(
                        &mut rows, name, label, "field", samples, ms, out_kib, worst, worst_near,
                    );
                }

                // Field → field, on a **degraded** input. Scaling by two keeps
                // the zero set and doubles the gradient, which is what a scale
                // or a careless composition does to a distance field.
                let degraded: Vec<f64> = grid.iter().map(|v| v * 2.0).collect();
                let (worst_in, near_in) = error(&degraded, size, h, &truth);
                emit(
                    &mut rows,
                    name,
                    "degraded-in",
                    "repair",
                    samples,
                    0.0,
                    0.0,
                    worst_in,
                    near_in,
                );
                let (ms, out_kib, built) = timed(&|| {
                    reinitialise_narrow_band(&degraded, &shape, h, BAND)
                        .expect("band")
                        .0
                });
                let (worst, worst_near) = error(&built, size, h, &truth);
                emit(
                    &mut rows, name, "band", "repair", samples, ms, out_kib, worst, worst_near,
                );

                // Mesh → field. The mesh is built once and shared, so the two
                // constructors are compared against each other and not against
                // two different meshings.
                let mesh = mesh_of(&grid, &shape, lo, h);
                let mesh_to_field: [(&str, &dyn Fn() -> Vec<f64>); 2] = [
                    ("pseudonormal", &|| {
                        signed_distance_from_mesh(&mesh.positions, &mesh.indices, &shape, lo, h)
                            .expect("pseudonormal")
                    }),
                    ("winding", &|| {
                        signed_distance_from_mesh_winding(
                            &mesh.positions,
                            &mesh.indices,
                            &shape,
                            lo,
                            h,
                            0.5,
                        )
                        .expect("winding")
                    }),
                ];
                for (label, run) in mesh_to_field {
                    let (ms, out_kib, built) = timed(run);
                    let (worst, worst_near) = error(&built, size, h, &truth);
                    emit(
                        &mut rows, name, label, "mesh", samples, ms, out_kib, worst, worst_near,
                    );
                }
            }
        }
    });

    let mut csv = String::from("field,constructor,family,samples,ms,out_kib,worst,worst_near\n");
    for r in &rows {
        let _ = writeln!(
            csv,
            "{},{},{},{},{:.6},{:.1},{:.8},{:.8}",
            r.field, r.constructor, r.family, r.samples, r.ms, r.out_kib, r.worst, r.worst_near
        );
    }
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements")
        .join("constructors.csv");
    if let Some(dir) = path.parent() {
        let _ = fs::create_dir_all(dir);
    }
    let _ = fs::write(&path, csv);
    println!("\nwrote {}", path.display());

    recommend(&rows);
}

/// Best of `RUNS`, plus the exact size of what came back.
///
/// The best rather than the mean, since the fastest run is the one with the
/// least interference from everything else on the machine.
fn timed(run: &dyn Fn() -> Vec<f64>) -> (f64, f64, Vec<f64>) {
    let mut ms = f64::INFINITY;
    let mut built = Vec::new();
    for _ in 0..RUNS {
        let start = Instant::now();
        built = run();
        ms = ms.min(start.elapsed().as_secs_f64() * 1e3);
        std::hint::black_box(&built);
    }
    let out_kib = (built.len() * core::mem::size_of::<f64>()) as f64 / 1024.0;
    (ms, out_kib, built)
}

/// Worst error anywhere, and worst within two cells of the surface.
fn error(
    built: &[f64],
    size: [u32; 3],
    h: f64,
    truth: &dyn Fn(u32, u32, u32) -> f64,
) -> (f64, f64) {
    let mut worst = 0.0f64;
    let mut worst_near = 0.0f64;
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let i = ((z * size[1] + y) * size[0] + x) as usize;
                let t = truth(x, y, z);
                let err = (built[i] - t).abs();
                worst = worst.max(err);
                if t.abs() <= 2.0 * h {
                    worst_near = worst_near.max(err);
                }
            }
        }
    }
    (worst, worst_near)
}

/// Marching Cubes over a sampled grid.
fn mesh_of(grid: &[f64], shape: &RuntimeShape3, origin: [f64; 3], h: f64) -> MeshBuffer<f64> {
    let mut out = MeshBuffer::new();
    let mut mc = MarchingCubes::<f64>::new();
    let sampled = isomesh::construct::SampledField::new(grid, shape, origin, h).expect("wrap");
    mc.extract(&sampled, shape, origin, h, &mut out)
        .expect("extraction");
    out
}

#[expect(
    clippy::too_many_arguments,
    reason = "one call site, and a struct literal here would be longer than the \
              call it replaces"
)]
fn emit(
    rows: &mut Vec<Row>,
    field: &'static str,
    constructor: &'static str,
    family: &'static str,
    samples: u32,
    ms: f64,
    out_kib: f64,
    worst: f64,
    worst_near: f64,
) {
    println!(
        "{field:<16} {constructor:<12} {samples:>7} {ms:>10.3} {out_kib:>11.1} \
         {worst:>11.5} {worst_near:>12.5}"
    );
    rows.push(Row {
        field,
        constructor,
        family,
        samples,
        ms,
        out_kib,
        worst,
        worst_near,
    });
}

/// The recommendation T-018 asks for, with the numbers behind it.
fn recommend(rows: &[Row]) {
    let mean = |c: &str, pick: &dyn Fn(&Row) -> f64| -> f64 {
        let v: Vec<f64> = rows
            .iter()
            .filter(|r| r.constructor == c)
            .map(pick)
            .collect();
        v.iter().sum::<f64>() / v.len().max(1) as f64
    };

    for (family, members) in [
        ("signs → field", &["exact", "swept", "marched"][..]),
        ("field → field (repair)", &["degraded-in", "band"][..]),
        ("mesh → field", &["pseudonormal", "winding"][..]),
    ] {
        println!("\nmean over every exact field and resolution — {family}:");
        println!(
            "{:<14} {:>10} {:>11} {:>11} {:>12}",
            "constructor", "ms", "out KiB", "worst", "worst_near"
        );
        for c in members {
            println!(
                "{c:<14} {:>10.3} {:>11.1} {:>11.5} {:>12.5}",
                mean(c, &|r| r.ms),
                mean(c, &|r| r.out_kib),
                mean(c, &|r| r.worst),
                mean(c, &|r| r.worst_near)
            );
        }
    }

    println!("\nrecommendation:");
    let fastest = ["exact", "swept", "marched"]
        .into_iter()
        .min_by(|a, b| {
            mean(a, &|r| r.ms)
                .partial_cmp(&mean(b, &|r| r.ms))
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .unwrap_or("exact");
    let nearest = ["exact", "swept", "marched"]
        .into_iter()
        .min_by(|a, b| {
            mean(a, &|r| r.worst_near)
                .partial_cmp(&mean(b, &|r| r.worst_near))
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .unwrap_or("exact");
    println!(
        "  fastest signs → field:      {fastest} ({:.3} ms mean)",
        mean(fastest, &|r| r.ms)
    );
    println!(
        "  reinitialisation buys:      near-surface {:.5} → {:.5} on a 2× degraded\n           field, and leaves the far field alone at {:.5} → {:.5} — which is what a\n           band is, not a shortfall.",
        mean("degraded-in", &|r| r.worst_near),
        mean("band", &|r| r.worst_near),
        mean("degraded-in", &|r| r.worst),
        mean("band", &|r| r.worst)
    );
    println!(
        "  most accurate near surface: {nearest} ({:.5} mean)",
        mean(nearest, &|r| r.worst_near)
    );
    println!(
        "  S-005 jump flooding is GPU-only and absent from this process. Its own\n  \
         measurement (M-257, docs/measurements/jump_flood.csv) beats every CPU\n  \
         constructor here against analytic truth on all twelve of its rows, so a\n  \
         consumer with a device should use it and read this table for the CPU\n  \
         fallback it does not have."
    );
}
