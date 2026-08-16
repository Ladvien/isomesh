//! **A-024 — is the 128³ penalty a 64 KiB plane stride, and would padding fix it?**
//!
//! ```bash
//! cargo bench --bench a024_aliasing
//! ```
//!
//! Writes `docs/measurements/a024-aliasing.csv`.
//!
//! # The thing being explained
//!
//! Once A-023 removed the store-to-load stall that was covering it (M-286),
//! Surface Nets costs **84.3 cycles per sample at 128³ against 33.2 at 127³ and
//! 33.0 at 129³** — working sets 2% apart, a 2.6× penalty, on a field with no
//! surface at all. `DualMesher::values` is `n³` floats indexed by
//! `shape.linearize`, so its plane stride is `n²·4` bytes, which at `n = 128` is
//! **exactly 64 KiB** — a cache-set aliasing period on this machine.
//!
//! # The measurement, which needs no library change
//!
//! The aliasing depends on `size[0]·size[1]·4`, and the caller chooses the
//! shape. So adding **one sample** to `x` or `y` moves the plane stride off the
//! power of two while changing the amount of work by under 1%. If the spike is
//! the stride, `128×129×128` costs what `127³` costs; if it survives, the stride
//! is not the mechanism and A-024's remedy (a) would have been built on a guess.
//!
//! **This is exactly what padding would do**, arranged by the caller instead of
//! the mesher, so it prices remedy (a)'s benefit before anything is implemented.
//! It does not price its *cost*, which is a multiply-add per index at every
//! resolution and is what the ticket's acceptance asks for separately.
//!
//! # The z axis is the control
//!
//! `128×128×129` adds a sample to the axis that is **not** part of the plane
//! stride. Its work grows by the same 0.8% and its stride does not move, so if
//! the spike is the stride that row keeps it and the other two lose it. Without
//! it, "adding a sample helped" would be indistinguishable from "the fixture got
//! slightly different".

mod common;

#[cfg(target_os = "linux")]
mod probe {
    use std::fmt::Write as _;
    use std::fs;
    use std::hint::black_box;
    use std::path::PathBuf;
    use std::time::Instant;

    use isomesh::extractor::Extractor;
    use isomesh::fields::Sphere;
    use isomesh::surface_nets::SurfaceNets;
    use isomesh::{MeshBuffer, RuntimeShape3};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    type Scalar = f32;
    const WARMUP_RUNS: u32 = 2;
    const TIMED_RUNS: usize = 5;
    /// Entirely inside the domain: no crossing, no vertex, only the scaffolding.
    const EMPTY_RADIUS: Scalar = 10.0;

    /// Shapes, grouped by what each is for.
    ///
    /// `label` names the group so the CSV can be read without reconstructing the
    /// argument: `baseline` is the spike and its two neighbours, `plane` breaks
    /// the plane stride, `depth` is the control that does not.
    const CASES: [(&str, [u32; 3]); 13] = [
        ("baseline", [127, 127, 127]),
        ("baseline", [128, 128, 128]),
        ("baseline", [129, 129, 129]),
        ("plane", [129, 128, 128]),
        ("plane", [128, 129, 128]),
        ("plane", [130, 128, 128]),
        ("plane", [128, 130, 128]),
        ("depth", [128, 128, 129]),
        ("depth", [128, 128, 130]),
        ("baseline", [255, 255, 255]),
        ("baseline", [256, 256, 256]),
        ("baseline", [257, 257, 257]),
        // Row stride alone: `size[0]·4` is 512 bytes here and the plane stride
        // is not a power of two, so this separates the two aliasing periods.
        ("row", [128, 131, 131]),
    ];

    struct Row {
        label: &'static str,
        size: [u32; 3],
        cycles_per_sample: f64,
        misses_per_sample: f64,
        ipc: f64,
        ns_per_sample: f64,
        plane_bytes: u64,
    }

    fn measure(label: &'static str, size: [u32; 3]) -> Row {
        let field = Sphere::<Scalar> {
            center: [0.0; 3],
            radius: EMPTY_RADIUS,
        };
        let shape = RuntimeShape3::new(size).expect("the fixture fits u32");
        let cell_size = 4.0 / (f64::from(size[0] - 1)) as Scalar;
        let mut extractor = SurfaceNets::<Scalar>::new();
        let mut mesh = MeshBuffer::<Scalar>::new();
        for _ in 0..WARMUP_RUNS {
            mesh.reset();
            extractor
                .extract_into(&field, &shape, [-2.0; 3], cell_size, &mut mesh)
                .expect("extraction");
            black_box(&mesh);
        }

        let mut probe = Probe::open();
        let mut runs: Vec<(u128, u64, u64, u64)> = Vec::with_capacity(TIMED_RUNS);
        for _ in 0..TIMED_RUNS {
            mesh.reset();
            probe.reset_and_enable();
            let started = Instant::now();
            extractor
                .extract_into(&field, &shape, [-2.0; 3], cell_size, &mut mesh)
                .expect("extraction");
            let nanos = started.elapsed().as_nanos();
            probe.disable();
            let counted = probe.read();
            assert!(
                counted.worst_ratio() >= MIN_TIME_RATIO,
                "a counter ran only {:.1}% of the time it was enabled",
                counted.worst_ratio() * 100.0
            );
            runs.push((
                nanos,
                counted.cycles.count,
                counted.instructions.count,
                counted.cache_misses.count,
            ));
            black_box(&mesh);
        }
        runs.sort_unstable();
        let (nanos, cycles, instructions, misses) = runs[TIMED_RUNS / 2];
        let samples = f64::from(size[0]) * f64::from(size[1]) * f64::from(size[2]);
        Row {
            label,
            size,
            cycles_per_sample: cycles as f64 / samples,
            misses_per_sample: misses as f64 / samples,
            ipc: instructions as f64 / cycles as f64,
            ns_per_sample: nanos as f64 / samples,
            plane_bytes: u64::from(size[0]) * u64::from(size[1]) * 4,
        }
    }

    pub(crate) fn run() {
        println!(
            "{:<10} {:>16} {:>13} {:>11} {:>10} {:>7} {:>10}",
            "group", "shape", "plane bytes", "cyc/sample", "miss/samp", "IPC", "ns/sample"
        );
        let mut rows = Vec::new();
        for (label, size) in CASES {
            let r = measure(label, size);
            println!(
                "{:<10} {:>16} {:>13} {:>11.2} {:>10.3} {:>7.2} {:>10.3}",
                r.label,
                format!("{}x{}x{}", r.size[0], r.size[1], r.size[2]),
                if r.plane_bytes.is_power_of_two() {
                    format!("{} (2^n)", r.plane_bytes)
                } else {
                    r.plane_bytes.to_string()
                },
                r.cycles_per_sample,
                r.misses_per_sample,
                r.ipc,
                r.ns_per_sample
            );
            rows.push(r);
        }

        let cyc = |size: [u32; 3]| {
            rows.iter()
                .find(|r| r.size == size)
                .map_or(f64::NAN, |r| r.cycles_per_sample)
        };
        let spike = cyc([128, 128, 128]);
        let clean = (cyc([127, 127, 127]) + cyc([129, 129, 129])) / 2.0;
        println!(
            "\n128³ against the mean of 127³ and 129³: {:.2}×",
            spike / clean
        );
        println!(
            "one sample added to x: {:.2}×   to y: {:.2}×   to z (control): {:.2}×",
            cyc([129, 128, 128]) / clean,
            cyc([128, 129, 128]) / clean,
            cyc([128, 128, 129]) / clean
        );
        let clean256 = (cyc([255, 255, 255]) + cyc([257, 257, 257])) / 2.0;
        println!(
            "256³ against the mean of 255³ and 257³: {:.2}×",
            cyc([256, 256, 256]) / clean256
        );
        println!(
            "row stride alone (128x131x131, 512-byte rows, non-power-of-two plane): {:.2}× the \
             mean of 127³/129³",
            cyc([128, 131, 131]) / clean
        );

        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
        let dir = root.join("docs/measurements");
        fs::create_dir_all(&dir).expect("create docs/measurements");
        let mut csv = String::from("# A-024: is the 128³ penalty a 64 KiB plane stride?\n");
        let _ = writeln!(
            csv,
            "group,shape,samples,plane_bytes,plane_is_power_of_two,cycles_per_sample,\
             cache_misses_per_sample,ipc,ns_per_sample"
        );
        for r in &rows {
            let samples = f64::from(r.size[0]) * f64::from(r.size[1]) * f64::from(r.size[2]);
            let _ = writeln!(
                csv,
                "{},{}x{}x{},{:.0},{},{},{:.4},{:.4},{:.4},{:.4}",
                r.label,
                r.size[0],
                r.size[1],
                r.size[2],
                samples,
                r.plane_bytes,
                u8::from(r.plane_bytes.is_power_of_two()),
                r.cycles_per_sample,
                r.misses_per_sample,
                r.ipc,
                r.ns_per_sample
            );
        }
        let path = dir.join("a024-aliasing.csv");
        fs::write(&path, csv).expect("write csv");
        println!("\n{} rows → {}", rows.len(), path.display());
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    #[cfg(target_os = "linux")]
    probe::run();
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "A-024's probe needs hardware counters, and this platform has no perf_event_open."
        );
        std::process::exit(1);
    }
}
