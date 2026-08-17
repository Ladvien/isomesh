//! **What the interior decider is the sign of, measured.**
//!
//! Ticket: R-023. **Exploratory, and deliberately not a registered experiment.**
//!
//! ```bash
//! cargo bench --bench interior_margin
//! ```
//!
//! Writes `docs/measurements/interior_margin.csv`.
//!
//! # Why there is no `P-` id on this
//!
//! P-24 was registered on `sign(F(s))` at a body saddle and is **degenerate**:
//! `F` is zero at every body saddle by construction, measured at 2.3e-15 worst
//! case (M-312). That is the second time in this phase a clause was registered
//! without checking its arithmetic, so this run registers nothing. It measures
//! the distribution of the quantity first, and a hypothesis follows the numbers
//! rather than preceding them into the same wall.
//!
//! # The quantity
//!
//! `SweptFaces::margin()` — the largest saddle value the sweep reaches, over the
//! candidate points `test()` walks. `test()` is now literally its sign, so the
//! interior decider is the `ε = 0` member of a one-parameter family by
//! construction rather than by agreement.
//!
//! **The open question this measures:** `saddle(t)` is a ratio whose denominator
//! vanishes at the pole, so the margin has no *a priori* bound. If it is
//! well-behaved, a threshold in the field's units is meaningful. If it has a
//! heavy tail near poled sweeps, a threshold is measuring the pole rather than
//! the tunnel, and the whole approach needs a different scalar.
//!
//! # It is a decision margin, not a persistence
//!
//! It answers *how far the field would have to move for the interior answer to
//! flip*. That is what a persistence threshold is wanted for, and it is not a
//! persistence pair (V-42). It also does **not**, on its own, retire A-002b,
//! A-002i or A-020b — those are about which topology to emit, and this is about
//! how confidently the cell was classified.

mod common;

use std::fmt::Write as _;

use isomesh::marching_cubes::interior::{Interior, SweptFaces};

/// Corner quadruples per face, in the cyclic order `SweptFaces` wants: `A`/`C`
/// one diagonal, `B`/`D` the other. `z = 0` is corners 0, 1, 3, 2.
const LO: [usize; 4] = [0, 1, 3, 2];
const HI: [usize; 4] = [4, 5, 7, 6];

/// The same generator the trilinear census uses, so the cell population is the
/// one M-220 and M-232 drew from.
struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    /// A corner value in `[-1, 1]`, quantised to `steps` levels when asked.
    fn value(&mut self, steps: Option<u32>) -> f64 {
        let raw = (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0;
        match steps {
            None => raw,
            Some(k) => (raw * f64::from(k)).round() / f64::from(k),
        }
    }

    fn corners(&mut self, steps: Option<u32>) -> [f64; 8] {
        let mut f = [0.0; 8];
        for v in &mut f {
            *v = self.value(steps);
        }
        f
    }
}

/// Is this face ambiguous — one diagonal strictly negative, the other not?
fn ambiguous(face: [f64; 4]) -> bool {
    let neg = |v: f64| v < 0.0;
    (neg(face[0]) == neg(face[2]))
        && (neg(face[1]) == neg(face[3]))
        && (neg(face[0]) != neg(face[1]))
}

struct Row {
    label: String,
    cells: u64,
    ambiguous: u64,
    joined: u64,
    poled: u64,
    finite: u64,
    max_abs: f64,
    /// Margins whose absolute value exceeds the largest corner magnitude — the
    /// sign that the ratio, not the field, is setting the scale.
    over_corner_scale: u64,
    /// Margins above **half** the corner scale. An ambiguous face forces the
    /// denominator's four terms to add rather than cancel, and AM-GM on the
    /// numerator then suggests `|saddle| <= max|corner| / 2`. Counted rather
    /// than asserted.
    over_half_scale: u64,
    tiny: u64,
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let mut rows = Vec::new();

    for (label, steps) in [
        ("continuous", None),
        ("quantum 1/255", Some(255u32)),
        ("quantum 1/16", Some(16)),
        ("quantum 1/4", Some(4)),
    ] {
        let mut rng = Lcg(0x0000_A023_D000_0007);
        let mut row = Row {
            label: label.to_string(),
            cells: 0,
            ambiguous: 0,
            joined: 0,
            poled: 0,
            finite: 0,
            max_abs: 0.0,
            over_corner_scale: 0,
            over_half_scale: 0,
            tiny: 0,
        };

        for _ in 0..400_000u32 {
            let f = rng.corners(steps);
            row.cells += 1;
            let lo = [f[LO[0]], f[LO[1]], f[LO[2]], f[LO[3]]];
            let hi = [f[HI[0]], f[HI[1]], f[HI[2]], f[HI[3]]];
            // The interior test is only defined on an ambiguous face.
            if !ambiguous(lo) || !ambiguous(hi) {
                continue;
            }
            let Ok(swept) = SweptFaces::new(lo, hi) else {
                continue;
            };
            row.ambiguous += 1;

            let margin = swept.margin();
            if swept.test() == Interior::Joined {
                row.joined += 1;
            }
            if swept.pole().is_some() {
                row.poled += 1;
            }
            if margin.is_finite() {
                row.finite += 1;
                row.max_abs = row.max_abs.max(margin.abs());
                let scale = f.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
                if margin.abs() > scale {
                    row.over_corner_scale += 1;
                }
                if margin.abs() > 0.5 * scale {
                    row.over_half_scale += 1;
                }
                if margin.abs() < 1e-6 {
                    row.tiny += 1;
                }
            }
        }
        rows.push(row);
    }

    println!(
        "{:<16} {:>10} {:>8} {:>7} {:>12} {:>11} {:>10} {:>7}",
        "population",
        "ambiguous",
        "joined",
        "poled",
        "max |margin|",
        "over scale",
        "over half",
        "tiny"
    );
    let mut csv = String::from(
        "population,cells,ambiguous,joined,poled,finite,max_abs_margin,over_corner_scale,over_half_scale,tiny\n",
    );
    for r in &rows {
        println!(
            "{:<16} {:>10} {:>8} {:>7} {:>12.4e} {:>11} {:>10} {:>7}",
            r.label,
            r.ambiguous,
            r.joined,
            r.poled,
            r.max_abs,
            r.over_corner_scale,
            r.over_half_scale,
            r.tiny
        );
        let _ = writeln!(
            csv,
            "{},{},{},{},{},{},{:.9e},{},{},{}",
            r.label,
            r.cells,
            r.ambiguous,
            r.joined,
            r.poled,
            r.finite,
            r.max_abs,
            r.over_corner_scale,
            r.over_half_scale,
            r.tiny
        );
    }

    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements/interior_margin.csv");
    match std::fs::write(&path, &csv) {
        Ok(()) => println!("\nwrote {}", path.display()),
        Err(e) => println!("\n::error:: {}: {e}", path.display()),
    }
}
