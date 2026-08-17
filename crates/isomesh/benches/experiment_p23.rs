//! **P-23 — does repairing air connectivity after a dig cost the edit, or the
//! lattice?**
//!
//! Ticket: R-022a. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p23
//! ```
//!
//! Writes `docs/experiments/p-23.csv`.
//!
//! # The design in one line
//!
//! **Hold the brush fixed and grow the lattice.** If repair cost is a property
//! of the edit, the union count does not move as `n³` grows by 60×. If it is a
//! property of the lattice, it tracks `n³`. Nothing else needs to vary.
//!
//! # Two arms, and the second is the control
//!
//! - **`incremental`** — `Air::dig`, the thing under test.
//! - **`rebuild`** — `Air::build` over the whole lattice, which is what an
//!   engine does today when it has no incremental structure. It is here so the
//!   incremental column has something to be flat *against*; a flat curve with
//!   nothing beside it is not evidence.
//!
//! **The rebuild's union count is NOT `Θ(n³)`, and P-23 said it would be.** A
//! union-find build unions only **air-air** edges, and the air volume is the
//! brush, which is fixed — so both union columns come out flat and the
//! comparison H proposed does not exist. What grows with `n³` is the **scan**:
//! the rebuild must visit every sample to discover which ones changed. That is
//! the `lattice_samples` column, and it is the honest control.
//!
//! # Counted, not timed
//!
//! Union calls are integers and identical on every machine (✗24). The
//! milliseconds are printed beside them and gated on nothing.

mod common;

use std::fmt::Write as _;

use isomesh::connectivity::Air;
use isomesh::{RuntimeShape3, Shape3};

/// Samples per axis. `129³` is 2.1 M samples, which is R-022's stated 128³.
const RESOLUTIONS: [u32; 3] = [33, 65, 129];

/// Brush radius in samples. Fixed across resolutions — that is the experiment.
const BRUSH_RADIUS: f64 = 6.0;

/// The samples a spherical brush of `BRUSH_RADIUS` centred at `c` covers.
fn brush(centre: [u32; 3]) -> Vec<[u32; 3]> {
    let r = BRUSH_RADIUS.ceil() as i64;
    let mut out = Vec::new();
    for dz in -r..=r {
        for dy in -r..=r {
            for dx in -r..=r {
                let d2 = (dx * dx + dy * dy + dz * dz) as f64;
                if d2 > BRUSH_RADIUS * BRUSH_RADIUS {
                    continue;
                }
                let p = [
                    i64::from(centre[0]) + dx,
                    i64::from(centre[1]) + dy,
                    i64::from(centre[2]) + dz,
                ];
                if p.iter().all(|c| *c >= 0) {
                    out.push([p[0] as u32, p[1] as u32, p[2] as u32]);
                }
            }
        }
    }
    out
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-23");
    common::experiment::run(prereg, |run| {
        println!(
            "{:>5} {:>10} {:>8} {:>13} {:>12} {:>9} {:>10} {:>11}",
            "n",
            "scanned",
            "dirty",
            "incr unions",
            "rebuild u",
            "per dirty",
            "incr ms",
            "rebuild ms"
        );

        let mut incremental_counts = Vec::new();
        let mut rebuild_counts = Vec::new();
        let mut over_degree = 0usize;

        for n in RESOLUTIONS {
            let shape = match RuntimeShape3::new([n; 3]) {
                Ok(shape) => shape,
                Err(e) => {
                    println!("::error:: {n}: {e}");
                    continue;
                }
            };
            let count = shape.element_count();

            // A solid block. The brush is the only air there is, so the dirty
            // set is exactly the brush and nothing about the field confounds
            // the measurement.
            let values = vec![-1.0_f64; count];

            let rebuild_start = std::time::Instant::now();
            let (mut air, rebuild) = match Air::build(&values, &shape) {
                Ok(pair) => pair,
                Err(e) => {
                    println!("::error:: {n}: {e}");
                    continue;
                }
            };
            let _initial_build_ms = rebuild_start.elapsed().as_secs_f64() * 1e3;

            let centre = [n / 2; 3];
            let cells = brush(centre);
            let dig_start = std::time::Instant::now();
            let incr = air.dig(&cells);
            let incr_ms = dig_start.elapsed().as_secs_f64() * 1e3;

            // The rebuild an engine without this structure would pay: the whole
            // lattice again, now that the brush has changed it.
            let mut after = values.clone();
            for c in &cells {
                if c[0] < n && c[1] < n && c[2] < n {
                    let i =
                        (c[2] as usize * n as usize + c[1] as usize) * n as usize + c[0] as usize;
                    if let Some(v) = after.get_mut(i) {
                        *v = 1.0;
                    }
                }
            }
            let after_start = std::time::Instant::now();
            let rebuild_after = match Air::build(&after, &shape) {
                Ok((_, r)) => r,
                Err(e) => {
                    println!("::error:: {n}: {e}");
                    continue;
                }
            };
            let rebuild_after_ms = after_start.elapsed().as_secs_f64() * 1e3;
            let _ = rebuild;

            if incr.unions > 6 * incr.dirty {
                over_degree += 1;
            }
            incremental_counts.push(incr.unions);
            rebuild_counts.push(rebuild_after.unions);

            println!(
                "{n:>5} {count:>10} {:>8} {:>13} {:>12} {:>9.2} {incr_ms:>10.3} {rebuild_after_ms:>11.1}",
                incr.dirty,
                incr.unions,
                rebuild_after.unions,
                incr.unions_per_dirty()
            );

            run.record(&[
                ("samples_per_axis", n.to_string()),
                ("dirty_samples", incr.dirty.to_string()),
                ("incremental_unions", incr.unions.to_string()),
                ("rebuild_unions", rebuild_after.unions.to_string()),
                (
                    "unions_per_dirty",
                    format!("{:.6}", incr.unions_per_dirty()),
                ),
                ("lattice_samples", count.to_string()),
                ("incremental_merges", incr.merges.to_string()),
                ("incremental_ms", format!("{incr_ms:.6}")),
                ("rebuild_ms", format!("{rebuild_after_ms:.6}")),
            ]);
        }

        println!();
        let flat = incremental_counts.windows(2).all(|w| w[0] == w[1]);
        println!(
            "clause 1: incremental union counts {incremental_counts:?} across {RESOLUTIONS:?} \
             -- {} (H says identical)",
            if flat { "IDENTICAL" } else { "NOT identical" }
        );
        println!(
            "clause 2: {over_degree} resolutions exceeded 6 unions per dirty sample (H says 0)"
        );
        let rebuild_flat = rebuild_counts.windows(2).all(|w| w[0] == w[1]);
        println!(
            "clause 3: rebuild union counts {rebuild_counts:?} -- {} \
             (H said these would grow as n^3; they do not, and the n^3 is in the scan)",
            if rebuild_flat { "FLAT" } else { "growing" }
        );
        let mut note = String::new();
        let _ = write!(
            note,
            "brush radius {BRUSH_RADIUS} samples, fixed across resolutions"
        );
        println!("{note}");
    });
}
