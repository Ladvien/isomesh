//! **P-27 — the exact-zero medial reading, offset by half a voxel.**
//!
//! Ticket: R-029. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p27
//! ```
//!
//! Writes `docs/experiments/p-27.csv`.
//!
//! # The design in one line
//!
//! **Move the slab, not the instrument.** M-172's exact `[0, 0, 0]` came from a
//! mid-plane the sampling lattice happened to sit on; if the zero belongs to
//! the lattice rather than the field, a half-voxel offset kills it, while the
//! band `‖∇ρ‖ < 0.1` — a property of the field — stays put.
//!
//! # Two populations, from the registration
//!
//! Under voxel-step central differencing the slab's gradient magnitude at
//! distance `d` from the mid-plane is exactly `d/h` for `|d| < h`, so the
//! sub-threshold set is the band `|d| < 0.1h`: thickness `0.2h`, thinner than
//! the lattice pitch. Counted on the voxel lattice, clause 2 would collapse
//! under ANY misalignment — vacuously false, for the instrument's reason
//! rather than the field's. So:
//!
//! - **Exact zeros** are counted on the voxel lattice with the *default*
//!   [`Sdf::gradient`] (`DIFF_STEP`-scaled — M-172's actual instrument).
//! - **The band** is counted on a regular probe lattice of pitch `h/200`
//!   spanning `|y| ≤ 2h`, gradient by [`central_difference`] at step `h` —
//!   the regime a game's voxel buffer actually has.
//!
//! # The probe phase is half a pitch, and that is load-bearing
//!
//! The registration fixes pitch, span and step; the *phase* is the harness's
//! choice. Probes at integer multiples of the pitch would sit exactly on the
//! band edge (`0.1h` is 20 pitches), where `< 0.1` is decided by rounding
//! crumbs ~1e-16 wide and the count would wobble between arms for floating
//! point's reasons. At half-pitch phase every probe-to-edge distance is a
//! half-integer multiple of the pitch: margin `0.5·pitch`, five orders above
//! rounding, and the registered arithmetic — 40 ± 1 in the window, worst rigid
//! offset 2.5% — becomes exact. The half-voxel offset is 100 pitches, so both
//! arms share a probe phase and the *predicted* clause-2 change is exactly 0.
//!
//! # Three predictions and one inversion
//!
//! - **aligned** — 4,225 exact zeros (65², the full plane), the reachability
//!   arm: the offset arm's zero is meaningful only because this one fires.
//! - **offset_half** — 0 exact zeros, band count unchanged.
//! - **offset_full** — **3,136** exact zeros, corrected pre-run in the
//!   registration: the cancellation needs `fl(y₀ + h_cd) − y₀ == h_cd`, which
//!   holds unconditionally only at `y₀ = 0.0`, so a full-voxel offset restores
//!   the zero only where `DIFF_STEP·scale` happens to survive the rounding —
//!   the exact zero is a *coordinate-origin* artifact, not merely a
//!   lattice-alignment one.
//! - **naive probe pitch** (the inversion): the band counted at voxel pitch
//!   must show the ≈100% collapse — the `> 5%` branch demonstrated reachable,
//!   and the vacuity trap the registration routes around shown as a row.
//!
//! # Counted, not timed
//!
//! Every column is an integer count, identical on every machine (✗24). No
//! interleaving: there is no timing A/B here for run order to bias — M-197
//! does not apply, and this note exists so the absence is not read as an
//! oversight.

mod common;

use isomesh::Sdf;
use isomesh::normals::central_difference;

/// Samples per axis. 65 over `[-2, 2]` gives `h = 2⁻⁴` — every lattice
/// coordinate and both plane positions representable, so an exact zero is a
/// property of the field rather than of a rounding accident.
const N: usize = 65;
/// Domain half-width, matching the reference fields' compact domain.
const DOMAIN: f64 = 2.0;
/// Voxel pitch: `4/64 = 2⁻⁴`.
const H: f64 = 2.0 * DOMAIN / 64.0;
/// Slab half-thickness: 8 voxels, so the mid-plane band is far from both faces.
const HALF_THICKNESS: f64 = 0.5;
/// The registered stability threshold on `‖∇ρ‖`.
const THRESHOLD: f64 = 0.1;
/// Probe columns per axis: every second lattice line, spanning the domain.
const COLUMNS: usize = 33;
/// Registered probe pitch, `h/200`.
const PROBE_PITCH: f64 = H / 200.0;
/// Probes per column at the registered pitch: `4h` of span at `h/200`.
const PROBES: usize = 800;

/// The slab, `|y − y₀| − t`. Deliberately not `BoxExact` or `ThinPlate`: their
/// edges and corners add gradient kinks that would contaminate the band count
/// with geometry the hypothesis is not about.
struct Slab {
    y0: f64,
}

impl Sdf for Slab {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        (p[1] - self.y0).abs() - HALF_THICKNESS
    }
}

/// Lattice coordinate `i` of 65 across `[-2, 2]`. Exact: `i·2⁻⁴` is
/// representable for every `i ≤ 64` and the sum shifts within one binade.
fn lattice(i: usize) -> f64 {
    -DOMAIN + (i as f64) * H
}

/// Exact-zero count over the full 65³ lattice, default gradient.
fn exact_zeros(slab: &Slab) -> u64 {
    let mut zeros = 0u64;
    for i in 0..N {
        for j in 0..N {
            for k in 0..N {
                let p = [lattice(i), lattice(j), lattice(k)];
                let g = slab.gradient(p);
                // Exact comparison is the measurement: M-172's claim is about
                // literal `[0.0, 0.0, 0.0]` returns, not about small ones.
                #[allow(
                    clippy::float_cmp,
                    reason = "the hypothesis is about exact zeros; a tolerance would measure a different claim"
                )]
                if g == [0.0, 0.0, 0.0] {
                    zeros += 1;
                }
            }
        }
    }
    zeros
}

/// Band count `‖∇ρ‖ < 0.1` at voxel-step central differences, on a regular
/// y-lattice of `probes` points at `pitch` and `phase` (in pitches) spanning
/// from `−2h`, replicated over `COLUMNS²` (x, z) columns.
///
/// Returns the total and asserts every column agrees — the field is
/// independent of x and z, so a column that disagrees is a harness bug, not a
/// result.
fn band_count(slab: &Slab, pitch: f64, phase: f64, probes: usize) -> u64 {
    let mut per_column: Option<u64> = None;
    let mut total = 0u64;
    for i in 0..COLUMNS {
        for j in 0..COLUMNS {
            let x = -DOMAIN + (i as f64) * (2.0 * H);
            let z = -DOMAIN + (j as f64) * (2.0 * H);
            let mut count = 0u64;
            for k in 0..probes {
                let y = -2.0 * H + ((k as f64) + phase) * pitch;
                let g = central_difference(slab, [x, y, z], H);
                let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                if mag < THRESHOLD {
                    count += 1;
                }
            }
            match per_column {
                None => per_column = Some(count),
                Some(first) => assert_eq!(
                    first, count,
                    "column ({i},{j}) disagrees — the field is x/z-invariant, \
                     so this is a harness bug"
                ),
            }
            total += count;
        }
    }
    total
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-27");
    common::experiment::run(prereg, |run| {
        let aligned = Slab { y0: 0.0 };
        let half = Slab { y0: H / 2.0 };
        let full = Slab { y0: H };

        // Clause 1: the voxel lattice, default gradient.
        let zeros_aligned = exact_zeros(&aligned);
        let zeros_half = exact_zeros(&half);
        let zeros_full = exact_zeros(&full);

        // Reachability: the offset arm's zero means something only because the
        // aligned arm demonstrates the counter fires, at exactly the
        // registered plane population.
        assert_eq!(
            zeros_aligned,
            (N * N) as u64,
            "reachability: the aligned mid-plane must return exactly one full \
             lattice plane of `[0,0,0]`s — if it does not, the instrument \
             cannot see the thing the hypothesis is about"
        );

        // Clause 2: the registered probe lattice, voxel-step differences.
        let band_aligned = band_count(&aligned, PROBE_PITCH, 0.5, PROBES);
        let band_half = band_count(&half, PROBE_PITCH, 0.5, PROBES);

        // Inversion: the same band at voxel pitch — the vacuity trap, shown.
        let naive_aligned = band_count(&aligned, H, 0.0, 5);
        let naive_half = band_count(&half, H, 0.0, 5);
        let naive_change =
            (naive_aligned.abs_diff(naive_half)) as f64 / (naive_aligned.max(1)) as f64;
        assert!(
            naive_change > 0.5,
            "inversion: at voxel probe pitch the band count must collapse \
             under the offset (predicted ≈100%); if it does not, the vacuity \
             analysis behind the two-population design is wrong and the \
             derived 2.5% bound is unfounded"
        );

        let band_change = (band_aligned.abs_diff(band_half)) as f64 / (band_aligned.max(1)) as f64;

        println!(
            "{:>14} {:>13} {:>12} {:>11} {:>18}",
            "arm", "offset_voxels", "exact_zeros", "band_count", "probe_pitch_voxels"
        );
        let rows: [(&str, f64, u64, u64, f64); 5] = [
            ("aligned", 0.0, zeros_aligned, band_aligned, 0.005),
            ("offset_half", 0.5, zeros_half, band_half, 0.005),
            (
                "offset_full",
                1.0,
                zeros_full,
                band_count(&full, PROBE_PITCH, 0.5, PROBES),
                0.005,
            ),
            ("naive_aligned", 0.0, zeros_aligned, naive_aligned, 1.0),
            ("naive_half", 0.5, zeros_half, naive_half, 1.0),
        ];
        for (arm, offset, zeros, band, pitch) in &rows {
            println!("{arm:>14} {offset:>13} {zeros:>12} {band:>11} {pitch:>18}");
            run.record(&[
                ("arm", (*arm).to_string()),
                ("offset_voxels", format!("{offset}")),
                ("exact_zeros", zeros.to_string()),
                ("band_count", band.to_string()),
                ("probe_pitch_voxels", format!("{pitch}")),
            ]);
        }

        println!();
        println!(
            "clause 1: exact zeros {} -> {} under a half-voxel offset -- {} (H says 4,225 -> 0)",
            zeros_aligned,
            zeros_half,
            if zeros_half == 0 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "clause 2: band count {} -> {} ({:.2}% change) -- {} (H says < 5%, instrument bound 2.5%, predicted 0%)",
            band_aligned,
            band_half,
            band_change * 100.0,
            if band_change < 0.05 {
                "HELD"
            } else {
                "FALSIFIED"
            }
        );
        println!(
            "corrected full-voxel prediction: {} exact zeros -- {} (pre-run arithmetic says exactly 3,136 of 4,225)",
            zeros_full,
            if zeros_full == 3_136 { "HIT" } else { "MISSED" }
        );
        println!(
            "inversion: naive voxel-pitch band {} -> {} ({:.0}% collapse) -- PREDICTED, and the \
             reason clause 2 needs its own population",
            naive_aligned,
            naive_half,
            naive_change * 100.0
        );
    });
}
