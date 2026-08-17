//! **P-37 — Sabine RT60 per air component: the cheapest audible thing in the dossier.**
//!
//! Ticket: R-036. Pre-registered in the commit before this one; the premise
//! correction (no area accumulator existed) is the registration's first line,
//! and the accumulator this measures landed with its own invariant tests.
//!
//! ```bash
//! cargo bench --bench experiment_p37
//! ```
//!
//! Writes `docs/experiments/p-37.csv`.
//!
//! # The two clauses
//!
//! - **RT60 on the breach frame.** Dig the tunnel voxel that merges two
//!   chambers; on the frame of the merge, `RT60 = 0.161·V/(α·S)` is two
//!   accumulator reads and a divide. Registered: **< 0.1 ms** — structurally,
//!   since it is O(1) arithmetic on maintained counters, which is the whole
//!   point of R-036's accumulator.
//! - **The Planeverb-style re-bake.** A 64×64 2-D FDTD slice (leapfrog
//!   pressure–velocity, CFL step, 1,000 steps ≈ half a second of audio —
//!   the length a decay measurement needs; damped edges). Registered:
//!   **< 30 ms** single-threaded. Cost is the claim; acoustic fidelity is
//!   not, and the public Planeverb C++ is the fidelity reference
//!   (`10.1111/cgf.14099`).
//!
//! # Recorded beside them, not registered
//!
//! Dig/fill wall costs with the accumulator in place, and the split rate of
//! a 200-fill sequence against M-319's one-in-six — a divergence THERE would
//! be news; the clauses holding is not, and this run expects no banner.
//!
//! # Scale
//!
//! One voxel = 0.25 m and α = 0.15 (rock-ish), registered constants for the
//! recorded RT60 value; the *cost* clauses do not depend on either.

mod common;

use isomesh::connectivity::Air;
use isomesh::{RuntimeShape3, Shape3};
use std::time::Instant;

const N: u32 = 64;
/// Game scale for the recorded RT60 value.
const VOXEL_M: f64 = 0.25;
const ABSORB: f64 = 0.15;
/// FDTD: 64×64 cells of `VOXEL_M`, c = 343 m/s, CFL step, 1,000 steps.
const FDTD_STEPS: usize = 1_000;

struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    fn pick(&mut self, n: u32) -> u32 {
        (self.next_u64() >> 33) as u32 % n
    }
}

fn ball(c: [u32; 3], r: u32) -> Vec<[u32; 3]> {
    let mut out = Vec::new();
    for x in c[0].saturating_sub(r)..(c[0] + r + 1).min(N) {
        for y in c[1].saturating_sub(r)..(c[1] + r + 1).min(N) {
            for z in c[2].saturating_sub(r)..(c[2] + r + 1).min(N) {
                let d = (i64::from(x) - i64::from(c[0])).pow(2)
                    + (i64::from(y) - i64::from(c[1])).pow(2)
                    + (i64::from(z) - i64::from(c[2])).pow(2);
                if d <= i64::from(r * r) {
                    out.push([x, y, z]);
                }
            }
        }
    }
    out
}

/// The Sabine estimate off the maintained counters — the thing C1 times.
fn rt60(air: &Air, at: [u32; 3]) -> Option<f64> {
    let l = air.label_of(at)?;
    let v = f64::from(air.component_size(l)) * VOXEL_M * VOXEL_M * VOXEL_M;
    let s = f64::from(air.component_area(l)) * VOXEL_M * VOXEL_M;
    if s <= 0.0 {
        return None;
    }
    Some(0.161 * v / (ABSORB * s))
}

/// One 64×64 FDTD bake: leapfrog pressure–velocity, damped edges, impulse at
/// a fixed source. Returns the terminal pressure L1 so the loop cannot be
/// optimised away.
fn fdtd_bake() -> f64 {
    const W: usize = 64;
    let c = 343.0f64;
    let dt = 0.99 * VOXEL_M / (c * 2.0f64.sqrt());
    let cdt = (c * dt / VOXEL_M).powi(2);
    let mut p = vec![0.0f64; W * W];
    let mut p_prev = vec![0.0f64; W * W];
    p[20 * W + 20] = 1.0;
    for _ in 0..FDTD_STEPS {
        let mut p_next = vec![0.0f64; W * W];
        for y in 1..W - 1 {
            for x in 1..W - 1 {
                let i = y * W + x;
                let lap = p[i - 1] + p[i + 1] + p[i - W] + p[i + W] - 4.0 * p[i];
                p_next[i] = (2.0 * p[i] - p_prev[i] + cdt * lap) * 0.9995;
            }
        }
        p_prev = p;
        p = p_next;
    }
    p.iter().map(|v| v.abs()).sum()
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_unstable_by(f64::total_cmp);
    v[v.len() / 2]
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-37");
    common::experiment::run(prereg, |run| {
        let shape = RuntimeShape3::new([N; 3]).expect("grid fits");
        let values = vec![-1.0f64; shape.element_count()];
        let (mut air, _) = Air::build(&values, &shape).expect("build");

        // Two chambers and a tunnel, all but the breach voxel.
        let mut dig = Vec::new();
        for x in 8..24u32 {
            for y in 24..40u32 {
                for z in 24..40u32 {
                    dig.push([x, y, z]);
                }
            }
        }
        for x in 40..56u32 {
            for y in 24..40u32 {
                for z in 24..40u32 {
                    dig.push([x, y, z]);
                }
            }
        }
        for x in 24..32u32 {
            dig.push([x, 32, 32]);
        }
        for x in 33..40u32 {
            dig.push([x, 32, 32]);
        }
        air.dig(&dig, || true);
        assert_eq!(air.components(), 2, "fixture: two chambers pre-breach");

        // The breach: one voxel merges them, on this frame.
        let breach = air.dig(&[[32, 32, 32]], || true);
        assert!(
            breach.merges >= 1,
            "reachability: the breach voxel merged nothing"
        );
        assert_eq!(air.components(), 1);
        let player = [12u32, 32, 32];
        assert!(air.component_area(air.label_of(player).expect("air")) > 0);

        // ---- C1: the RT60 read, timed hot ----------------------------------
        let mut acc = 0.0f64;
        let reps = 10_000u32;
        let t0 = Instant::now();
        for _ in 0..reps {
            acc += rt60(&air, player).unwrap_or(0.0);
        }
        let per_call_us = t0.elapsed().as_secs_f64() * 1e6 / f64::from(reps);
        let rt60_value = rt60(&air, player).unwrap_or(0.0);
        assert!(acc.is_finite());

        // ---- C2: the FDTD bake, 11 runs, median ----------------------------
        let mut bake_ms = Vec::new();
        let mut sink = 0.0f64;
        for _ in 0..11 {
            let t = Instant::now();
            sink += fdtd_bake();
            bake_ms.push(t.elapsed().as_secs_f64() * 1e3);
        }
        assert!(
            sink.is_finite() && sink > 0.0,
            "FDTD energy vanished or diverged — the bake did not run"
        );
        let fdtd_med = median(bake_ms);

        // ---- recorded: 200-fill sequence, split rate and op costs ----------
        let mut lcg = Lcg(0x5EED_5EED_5EED_5EED);
        let mut fills = 0u64;
        let mut splits = 0u64;
        let mut fill_ms = Vec::new();
        let mut dig_ms = Vec::new();
        for step in 0..200 {
            let c = [8 + lcg.pick(48), 24 + lcg.pick(16), 24 + lcg.pick(16)];
            let r = 1 + lcg.pick(2);
            let b = ball(c, r);
            if step % 2 == 0 {
                let t = Instant::now();
                let out = air.fill(&b, || true);
                fill_ms.push(t.elapsed().as_secs_f64() * 1e3);
                fills += 1;
                splits += out.splits;
            } else {
                let t = Instant::now();
                air.dig(&b, || true);
                dig_ms.push(t.elapsed().as_secs_f64() * 1e3);
            }
        }
        let split_rate = splits as f64 / fills as f64;

        // ---- emit -----------------------------------------------------------
        let rows: Vec<(&str, f64, &str, &str, &str)> = vec![
            (
                "rt60_compute",
                per_call_us,
                "us",
                "100",
                if per_call_us < 100.0 { "true" } else { "false" },
            ),
            ("rt60_value", rt60_value, "s", "-", "-"),
            (
                "fdtd_bake",
                fdtd_med,
                "ms",
                "30",
                if fdtd_med < 30.0 { "true" } else { "false" },
            ),
            ("splits_per_fill", split_rate, "ratio", "-", "-"),
            ("fill_ms_median", median(fill_ms.clone()), "ms", "-", "-"),
            ("dig_ms_median", median(dig_ms.clone()), "ms", "-", "-"),
        ];
        println!(
            "\n{:>16} {:>12} {:>6} {:>7} {:>6}",
            "quantity", "value", "unit", "bound", "held"
        );
        for (q, v, u, b, h) in &rows {
            println!("{q:>16} {v:>12.4} {u:>6} {b:>7} {h:>6}");
            run.record(&[
                ("quantity", (*q).to_string()),
                ("value", format!("{v:.5}")),
                ("unit", (*u).to_string()),
                ("bound", (*b).to_string()),
                ("held", (*h).to_string()),
            ]);
        }

        println!();
        let c1 = per_call_us < 100.0;
        let c2 = fdtd_med < 30.0;
        println!(
            "C1 (RT60 on the breach frame): {per_call_us:.3} µs/call -- {} (H says < 0.1 ms; \
             falsified at 0.3 ms)",
            if c1 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C2 (64×64 FDTD re-bake): {fdtd_med:.2} ms -- {} (H says < 30 ms; falsified at 90 ms)",
            if c2 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "recorded: split rate {split_rate:.3}/fill against M-319's ~0.167; fill median \
             {:.3} ms, dig median {:.3} ms with the accumulator in place",
            median(fill_ms),
            median(dig_ms)
        );
    });
}
