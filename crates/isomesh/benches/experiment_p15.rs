//! **P-15 — where does the dual mesher's work go?**
//!
//! Ticket: R-007. Pre-registered at R-005.
//!
//! ```bash
//! cargo bench --bench experiment_p15
//! ```
//!
//! Writes `docs/experiments/p-15.csv`. **Linux only**, for the same reason as
//! `experiment_p12`: the counters are `perf_event_open`.
//!
//! # The question R-005 left
//!
//! M-279: Surface Nets executes **207 instructions per sample at IPC 1.22**
//! where Marching Cubes executes **132 at 4.04** — 1.57× the instructions and
//! 5.24× the cycles — and none of the candidates survived. Not the crossed-edge
//! gather (a field with no surface costs the same to within 0.9%), not branches
//! (they fall), not allocation (zero page faults), not the TLB, and not the
//! misses (a 2.4× swing moves cycles by 0.4%). M-282 then found the same IPC
//! wall under every rule built on `DualMesher` — 1.20, 1.35, 1.42 — and nowhere
//! else in the family.
//!
//! # How this measures a private function without touching it
//!
//! R-007's ticket offered two ways in and asked which is acceptable: counter
//! windows inside `DualMesher::extract`, or an ablation seam. **Neither is
//! used.** The first is impossible before it is undesirable — `isomesh` is
//! `no_std` and cannot call a Linux system call — and the second is public API
//! for one experiment.
//!
//! There is a third way, and it needs nothing from the library. The stages have
//! **different iteration counts**, and the counts depend on the grid's *shape*
//! rather than only its size:
//!
//! | stage | iterations |
//! |---|---|
//! | `sample` | `S = ∏ size[i]` |
//! | `resize` ×2, `place_vertices` | `C = ∏ (size[i] − 1)` |
//! | `emit_quads` | `Q = Σ_axis (size[axis] − 1)·(cells[u] − 1)·(cells[v] − 1)` |
//!
//! On a cube `Q/C ≈ 3`; on a slab two samples deep it is **1**; on a rod two
//! samples deep on *both* minor axes the inner loops are empty and it is
//! **exactly 0**, while `S/C` runs from 1 to 4 over the same shapes. So a sweep
//! of shapes with the same code and no surface makes the design matrix
//! separable, and least squares reads off the per-iteration cost of each stage
//! without a seam, a syscall, or a line of library change.
//!
//! **Instructions are fitted as well as cycles, and that is the load-bearing
//! half.** An instruction count is deterministic and cache-independent, so its
//! fit is a decomposition of the work rather than of the machine; if the cycle
//! fit is much worse, the difference is where the stalls are.
//!
//! # The field has no surface, and that is not a simplification
//!
//! M-279 measured the dual costing the same to within 0.9% with 153,552
//! triangles and with none, so the scaffolding *is* the cost and the empty field
//! removes the vertex-count terms that would otherwise need their own column.
//! The cubic shapes are run on `sphere` too, so the size of what is being
//! ignored is on the record rather than asserted.

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::hint::black_box;
    use std::time::Instant;

    use isomesh::extractor::Extractor;
    use isomesh::fields::Sphere;
    use isomesh::surface_nets::SurfaceNets;
    use isomesh::{MeshBuffer, RuntimeShape3};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    type Scalar = f32;

    /// Untimed runs first; the counted run is steady state.
    const WARMUP_RUNS: u32 = 2;
    /// Counted runs; the median by wall time is the one reported.
    const TIMED_RUNS: usize = 5;

    /// A sphere of radius 10 over a domain of half-extent 2 is entirely inside:
    /// no edge is crossed and no vertex is placed. See the module docs.
    const EMPTY_RADIUS: Scalar = 10.0;

    /// The shapes, chosen so `Q/C` spans 0 to 3 and `S/C` spans 1 to 4.
    ///
    /// Sizes are held within a factor of four of each other in total cells, so
    /// no single row dominates the least squares, and the two-deep shapes appear
    /// in all three orientations because the same `Q/C` reached by a different
    /// axis is the check that the fit is reading iteration counts and not memory
    /// layout.
    const SHAPES: [[u32; 3]; 13] = [
        [129, 129, 129],
        [161, 161, 161],
        [193, 193, 193],
        [257, 257, 33],
        [513, 513, 9],
        [725, 725, 5],
        [1025, 1025, 3],
        [1449, 1449, 2],
        [2, 1449, 1449],
        [1449, 2, 1449],
        [500_001, 2, 2],
        [2, 2, 500_001],
        [1_000_001, 2, 2],
    ];

    /// The shape held back from the fit, to check it predicts rather than
    /// merely describes.
    const HELD_OUT: [u32; 3] = [385, 385, 17];

    /// Iteration counts for one grid, from the loop bounds in `dual.rs`.
    fn counts(size: [u32; 3]) -> (f64, f64, f64) {
        let cells = [size[0] - 1, size[1] - 1, size[2] - 1];
        let s = f64::from(size[0]) * f64::from(size[1]) * f64::from(size[2]);
        let c = f64::from(cells[0]) * f64::from(cells[1]) * f64::from(cells[2]);
        let mut q = 0.0;
        for (axis, extent) in size.iter().enumerate() {
            let u = (axis + 1) % 3;
            let v = (axis + 2) % 3;
            // `for a in 0..size[axis] - 1`, `for b in 1..cells[u]`,
            // `for c in 1..cells[v]`.
            let a = f64::from(extent - 1);
            let b = f64::from(cells[u].saturating_sub(1));
            let cc = f64::from(cells[v].saturating_sub(1));
            q += a * b * cc;
        }
        (s, c, q)
    }

    /// One measured shape.
    pub(crate) struct Row {
        pub(crate) size: [u32; 3],
        pub(crate) field: &'static str,
        pub(crate) s: f64,
        pub(crate) c: f64,
        pub(crate) q: f64,
        pub(crate) cycles: f64,
        pub(crate) instructions: f64,
        pub(crate) nanos: f64,
        pub(crate) triangles: usize,
    }

    fn measure(size: [u32; 3], radius: Scalar, field_name: &'static str) -> Row {
        let field = Sphere::<Scalar> {
            center: [0.0; 3],
            radius,
        };
        let shape = RuntimeShape3::new(size).expect("the fixture fits u32");
        let cell_size = 4.0 / f64::from(size[0].max(size[1]).max(size[2]) - 1) as Scalar;
        let origin = [-2.0; 3];
        let mut extractor = SurfaceNets::<Scalar>::new();
        let mut mesh = MeshBuffer::<Scalar>::new();

        for _ in 0..WARMUP_RUNS {
            mesh.reset();
            extractor
                .extract_into(&field, &shape, origin, cell_size, &mut mesh)
                .expect("extraction");
            black_box(&mesh);
        }

        let mut probe = Probe::open();
        let mut runs: Vec<(u128, u64, u64)> = Vec::with_capacity(TIMED_RUNS);
        for _ in 0..TIMED_RUNS {
            mesh.reset();
            probe.reset_and_enable();
            let started = Instant::now();
            extractor
                .extract_into(&field, &shape, origin, cell_size, &mut mesh)
                .expect("extraction");
            let nanos = started.elapsed().as_nanos();
            probe.disable();
            let counted = probe.read();
            assert!(
                counted.worst_ratio() >= MIN_TIME_RATIO,
                "a counter ran only {:.1}% of the time it was enabled",
                counted.worst_ratio() * 100.0
            );
            runs.push((nanos, counted.cycles.count, counted.instructions.count));
            black_box(&mesh);
        }
        runs.sort_unstable();
        let (nanos, cycles, instructions) = runs[TIMED_RUNS / 2];

        let (s, c, q) = counts(size);
        Row {
            size,
            field: field_name,
            s,
            c,
            q,
            cycles: cycles as f64,
            instructions: instructions as f64,
            nanos: nanos as f64,
            triangles: mesh.triangle_count(),
        }
    }

    /// Least squares for `y = a·S + b·C + e·Q`, no intercept.
    ///
    /// No intercept because every stage is a loop over one of the three counts
    /// and a constant term would have nothing to represent but the allocator —
    /// which is itself `O(C)` bytes. Solved by normal equations on a 3×3
    /// symmetric system; the shapes are chosen so it is well conditioned, and
    /// the held-out prediction is what says whether that worked.
    fn fit(rows: &[&Row], y: impl Fn(&Row) -> f64) -> [f64; 3] {
        let mut ata = [[0.0f64; 3]; 3];
        let mut atb = [0.0f64; 3];
        for r in rows {
            let x = [r.s, r.c, r.q];
            let yi = y(r);
            for i in 0..3 {
                atb[i] += x[i] * yi;
                for j in 0..3 {
                    ata[i][j] += x[i] * x[j];
                }
            }
        }
        // Gaussian elimination with partial pivoting on a 3×3.
        let mut m = [
            [ata[0][0], ata[0][1], ata[0][2], atb[0]],
            [ata[1][0], ata[1][1], ata[1][2], atb[1]],
            [ata[2][0], ata[2][1], ata[2][2], atb[2]],
        ];
        for col in 0..3 {
            let mut pivot = col;
            for row in col + 1..3 {
                if m[row][col].abs() > m[pivot][col].abs() {
                    pivot = row;
                }
            }
            m.swap(col, pivot);
            let d = m[col][col];
            if d == 0.0 {
                return [0.0; 3];
            }
            for value in m[col].iter_mut().skip(col) {
                *value /= d;
            }
            let pivot_row = m[col];
            for (row, cells) in m.iter_mut().enumerate() {
                if row != col {
                    let f = cells[col];
                    for (k, value) in cells.iter_mut().enumerate().skip(col) {
                        *value -= f * pivot_row[k];
                    }
                }
            }
        }
        [m[0][3], m[1][3], m[2][3]]
    }

    fn r_squared(rows: &[&Row], y: impl Fn(&Row) -> f64 + Copy, k: [f64; 3]) -> f64 {
        let mean = rows.iter().map(|r| y(r)).sum::<f64>() / rows.len() as f64;
        let mut ss_res = 0.0;
        let mut ss_tot = 0.0;
        for r in rows {
            let predicted = k[0] * r.s + k[1] * r.c + k[2] * r.q;
            ss_res += (y(r) - predicted).powi(2);
            ss_tot += (y(r) - mean).powi(2);
        }
        if ss_tot == 0.0 {
            1.0
        } else {
            1.0 - ss_res / ss_tot
        }
    }

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        println!(
            "{:>20} {:>10} {:>6} {:>6} {:>10} {:>10} {:>8}",
            "shape", "cells", "S/C", "Q/C", "cyc/cell", "ins/cell", "tris"
        );
        let mut rows: Vec<Row> = Vec::new();
        for size in SHAPES {
            rows.push(measure(size, EMPTY_RADIUS, "empty"));
        }
        rows.push(measure(HELD_OUT, EMPTY_RADIUS, "empty"));
        // The cubic shapes on a real surface, so the size of what the empty
        // field ignores is on the record.
        for size in [[129u32; 3], [193; 3]] {
            rows.push(measure(size, 1.0, "sphere"));
        }

        for r in &rows {
            println!(
                "{:>20} {:>10.0} {:>6.2} {:>6.2} {:>10.2} {:>10.2} {:>8}",
                format!("{}x{}x{}", r.size[0], r.size[1], r.size[2]),
                r.c,
                r.s / r.c,
                r.q / r.c,
                r.cycles / r.c,
                r.instructions / r.c,
                r.triangles
            );
        }

        let training: Vec<&Row> = rows
            .iter()
            .filter(|r| r.field == "empty" && r.size != HELD_OUT)
            .collect();
        let ins = fit(&training, |r| r.instructions);
        let cyc = fit(&training, |r| r.cycles);
        let ins_r2 = r_squared(&training, |r| r.instructions, ins);
        let cyc_r2 = r_squared(&training, |r| r.cycles, cyc);

        println!(
            "\nper-iteration fit over {} shapes, no intercept:\n  \
             instructions:  sample {:>8.3}   place {:>8.3}   emit_quads {:>8.3}   r² = {:.6}\n  \
             cycles:        sample {:>8.3}   place {:>8.3}   emit_quads {:>8.3}   r² = {:.6}",
            training.len(),
            ins[0],
            ins[1],
            ins[2],
            ins_r2,
            cyc[0],
            cyc[1],
            cyc[2],
            cyc_r2
        );

        // The answer to P-15, on the cube the rest of the repo measures.
        let cube = rows
            .iter()
            .find(|r| r.size == [193, 193, 193] && r.field == "empty")
            .expect("the 193³ row is in SHAPES");
        let share = |k: [f64; 3]| k[2] * cube.q / (k[0] * cube.s + k[1] * cube.c + k[2] * cube.q);
        let ins_share = share(ins) * 100.0;
        let cyc_share = share(cyc) * 100.0;

        let held = rows
            .iter()
            .find(|r| r.size == HELD_OUT)
            .expect("the held-out row was measured");
        let predict = |k: [f64; 3]| k[0] * held.s + k[1] * held.c + k[2] * held.q;
        let ins_err = (predict(ins) / held.instructions - 1.0) * 100.0;
        let cyc_err = (predict(cyc) / held.cycles - 1.0) * 100.0;

        println!(
            "\nheld-out {}x{}x{}: instructions predicted {:+.2}%, cycles {:+.2}%",
            HELD_OUT[0], HELD_OUT[1], HELD_OUT[2], ins_err, cyc_err
        );
        println!(
            "emit_quads' share at 193³: {ins_share:.1}% of instructions, {cyc_share:.1}% of cycles"
        );

        // **The same conclusion without the least squares.** `1025x1025x3` and
        // `1449x1449x2` have the same cell count to 0.02% and differ in `Q/C` by
        // exactly 1, so the difference in cost per cell *is* the cost of one
        // `emit_quads` iteration, with only the small `S/C` change to correct
        // for. A conclusion that needs a three-parameter fit to see is weaker
        // than one two rows can show.
        let pair = |a: [u32; 3], b: [u32; 3]| -> Option<(f64, f64)> {
            let x = rows.iter().find(|r| r.size == a && r.field == "empty")?;
            let y = rows.iter().find(|r| r.size == b && r.field == "empty")?;
            let dq = x.q / x.c - y.q / y.c;
            let ds = x.s / x.c - y.s / y.c;
            Some((
                (x.cycles / x.c - y.cycles / y.c - ds * cyc[0]) / dq,
                (x.instructions / x.c - y.instructions / y.c - ds * ins[0]) / dq,
            ))
        };
        if let Some((cyc_each, ins_each)) = pair([1025, 1025, 3], [1449, 1449, 2]) {
            println!(
                "fit-free, from two rows of equal cell count: emit_quads costs {cyc_each:.1} \
                 cycles and {ins_each:.1} instructions per iteration"
            );
            println!(
                "  (least squares says {:.1} and {:.1}; the rods, where Q is exactly 0, say the \
                 rest of the mesher costs {:.1} cycles per cell)",
                cyc[2],
                ins[2],
                cyc[0] + cyc[1]
            );
        }

        for r in &rows {
            let stage_share = |k: [f64; 3], i: usize| {
                let x = [r.s, r.c, r.q];
                k[i] * x[i] / (k[0] * r.s + k[1] * r.c + k[2] * r.q)
            };
            run.record(&[
                (
                    "stage",
                    format!("{}x{}x{}", r.size[0], r.size[1], r.size[2]),
                ),
                ("cycles_per_sample", format!("{:.4}", r.cycles / r.s)),
                (
                    "instructions_per_sample",
                    format!("{:.4}", r.instructions / r.s),
                ),
                ("ipc", format!("{:.4}", r.instructions / r.cycles)),
                ("samples", format!("{:.0}", r.s)),
                ("field", r.field.to_string()),
                ("cells", format!("{:.0}", r.c)),
                ("emit_quad_iterations", format!("{:.0}", r.q)),
                ("ns_per_cell", format!("{:.4}", r.nanos / r.c)),
                ("triangles", r.triangles.to_string()),
                (
                    "emit_share_of_cycles",
                    format!("{:.4}", stage_share(cyc, 2)),
                ),
                (
                    "emit_share_of_instructions",
                    format!("{:.4}", stage_share(ins, 2)),
                ),
                ("held_out", u8::from(r.size == HELD_OUT).to_string()),
            ]);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-15");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} needs hardware performance counters, and this platform has no `perf_event_open`.",
            prereg.id
        );
        std::process::exit(1);
    }
}
