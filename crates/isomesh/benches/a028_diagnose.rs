//! **A-028: which branch of `DegenerateNormal` fires on `bonsai`, and why.**
//!
//! Ticket: A-028. **Diagnostic, not an experiment** — it exists to answer the
//! ticket's acceptance question (*which* of the two causes the message names
//! actually fired, and at which cell) before anything is changed.
//!
//! ```bash
//! ./scripts/fetch_volumes.sh
//! cargo bench --bench a028_diagnose
//! ```
//!
//! # What the message cannot tell you
//!
//! `Error::DegenerateNormal` reads *"a zero gradient, or no incident area"* and
//! is raised from **two** places with **different** causes: `normals.rs`, where a
//! vertex has no incident triangle area, and `subgrid/extract.rs:548`, where
//! `sdf.gradient(p)` has zero length. Only one of them can be the `bonsai`
//! failure, because only the subgrid extractor refuses that volume — so the
//! branch is already decided by the call site, and the message is the thing that
//! obscured it.
//!
//! This locates **where** and **why**.

// Every float comparison here is an exact-equality test on purpose: the whole
// diagnostic is about values that are *bit-identical* where they should differ,
// so a tolerance would hide the thing being measured.
#![allow(clippy::float_cmp, reason = "exact equality is the phenomenon")]

mod common;

use std::path::PathBuf;

use isomesh::construct::SampledField;
use isomesh::{RuntimeShape3, Sdf};

const FILE: &str = "bonsai_256x256x256_uint8.raw";
const N: u32 = 256;
const ISO: f64 = 32.0;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements/volumes")
        .join(FILE);
    let Ok(bytes) = std::fs::read(&path) else {
        println!("{FILE} absent — run ./scripts/fetch_volumes.sh (M-006). Skipping.");
        return;
    };
    let want = N as usize * N as usize * N as usize;
    if bytes.len() != want {
        println!("::error:: {} bytes, expected {want}", bytes.len());
        return;
    }

    let values: Vec<f64> = bytes.iter().map(|b| ISO - f64::from(*b)).collect();
    let Ok(shape) = RuntimeShape3::new([N; 3]) else {
        return;
    };
    let Ok(field) = SampledField::new(&values, &shape, [0.0; 3], 1.0) else {
        return;
    };

    let at = |x: u32, y: u32, z: u32| {
        values[(z as usize * N as usize + y as usize) * N as usize + x as usize]
    };

    // Walk grid corners in scan order — the order the extractor visits cells —
    // and stop at the first one on a surface-crossing cell whose gradient has
    // zero length. That is the exact condition `subgrid/extract.rs:548` tests.
    let mut found = 0usize;
    let mut scanned = 0u64;
    let mut zero_total = 0u64;
    'outer: for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                // Only corners of cells the surface actually crosses matter.
                if x + 1 >= N || y + 1 >= N || z + 1 >= N {
                    continue;
                }
                let mut lo = f64::INFINITY;
                let mut hi = f64::NEG_INFINITY;
                for dz in 0..2 {
                    for dy in 0..2 {
                        for dx in 0..2 {
                            let v = at(x + dx, y + dy, z + dz);
                            lo = lo.min(v);
                            hi = hi.max(v);
                        }
                    }
                }
                if !(lo < 0.0 && hi >= 0.0) {
                    continue;
                }
                scanned += 1;

                let p = [f64::from(x), f64::from(y), f64::from(z)];
                let g = field.gradient(p);
                let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                if len.is_finite() && len > 0.0 {
                    continue;
                }
                zero_total += 1;
                if found < 3 {
                    found += 1;
                    println!("--- zero gradient #{found} at corner [{x}, {y}, {z}] ---");
                    println!("  gradient {g:?}  length {len:?}");
                    println!("  the eight corner values of this cell (iso - density):");
                    for dz in 0..2 {
                        for dy in 0..2 {
                            for dx in 0..2 {
                                let v = at(x + dx, y + dy, z + dz);
                                println!(
                                    "    [{}, {}, {}] = {v:>6.1}   density {:>3}",
                                    x + dx,
                                    y + dy,
                                    z + dz,
                                    ISO - v
                                );
                            }
                        }
                    }
                    // Which cell does `floor` put this corner in, and is that
                    // cell uniform? A corner sits on the boundary of eight
                    // cells and `sample`/`gradient` both read the one at
                    // `floor(t)` -- the +side cell. If the surface is in a
                    // -side cell and the +side one is uniform, the analytic
                    // gradient is legitimately zero there.
                    {
                        let mut uniform = true;
                        let first = at(x, y, z);
                        for dz in 0..2 {
                            for dy in 0..2 {
                                for dx in 0..2 {
                                    if at(x + dx, y + dy, z + dz) != first {
                                        uniform = false;
                                    }
                                }
                            }
                        }
                        println!(
                            "  the +side cell (the one floor picks) is {}",
                            if uniform {
                                "UNIFORM -- analytic gradient is zero here"
                            } else {
                                "not uniform"
                            }
                        );
                    }
                    // What the central difference actually sampled.
                    let scale = p[0].abs().max(p[1].abs()).max(p[2].abs()).max(1.0);
                    let h = f64::EPSILON.cbrt() * scale;
                    println!("  DIFF_STEP*scale = {h:.6e}   (voxel is 1.0, so h/voxel = {h:.3e})");
                    for axis in 0..3 {
                        let mut a = p;
                        let mut b = p;
                        a[axis] += h;
                        b[axis] -= h;
                        println!(
                            "    axis {axis}: sample(+h) = {:.17}  sample(-h) = {:.17}  equal = {}",
                            field.sample(a),
                            field.sample(b),
                            field.sample(a) == field.sample(b)
                        );
                    }
                    // **The reading to confirm: is this a local extremum?** A
                    // central difference is zero at one whatever the slopes
                    // either side are, and quantised data makes them common.
                    let here = at(x, y, z);
                    let names = ["x", "y", "z"];
                    for axis in 0..3 {
                        let mut lo_p = [x, y, z];
                        let mut hi_p = [x, y, z];
                        if lo_p[axis] == 0 || hi_p[axis] + 1 >= N {
                            println!("    {} neighbours: on the grid boundary", names[axis]);
                            continue;
                        }
                        lo_p[axis] -= 1;
                        hi_p[axis] += 1;
                        let before = at(lo_p[0], lo_p[1], lo_p[2]);
                        let after = at(hi_p[0], hi_p[1], hi_p[2]);
                        let extremum = (before - here).signum() == (after - here).signum()
                            && (before - here) != 0.0;
                        let verdict = if extremum {
                            "LOCAL EXTREMUM"
                        } else {
                            "monotone"
                        };
                        let axis_name = names[axis];
                        let rise = here - before;
                        let fall = after - here;
                        println!(
                            "    {axis_name} neighbours: {before:>6.1} [{here:>6.1}] {after:>6.1}   \
                             slopes {rise:+.1} / {fall:+.1}  -> {verdict}"
                        );
                    }
                    if found == 3 {
                        break 'outer;
                    }
                }
            }
        }
    }

    println!();
    println!(
        "{zero_total} zero-length gradients among {scanned} surface-cell corners scanned \
         (stopped after {found} reports)"
    );
    if zero_total == 0 {
        println!(
            "none found at cell corners — the failing position is not a corner, so the \
             extractor is evaluating somewhere else and this diagnostic is looking in the \
             wrong place"
        );
    }
}
