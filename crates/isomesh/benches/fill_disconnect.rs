//! **How often does filling actually disconnect the air region?**
//!
//! Ticket: R-022b. **Exploratory. Nothing is registered against this run** — it
//! checks the arithmetic of a hypothesis before one is written, which is now the
//! house style: this phase failed twice by registering first (P-23 clause 3,
//! P-24) and got it right five times by measuring first.
//!
//! ```bash
//! cargo bench --bench fill_disconnect
//! ```
//!
//! Writes `docs/measurements/fill_disconnect.csv`.
//!
//! # The question R-022b turns on
//!
//! V-41 split R-022: **digging** only inserts, so a union-find answers it and
//! M-311 measured the repair at `O(|edit|)`. **Filling** deletes, and a deletion
//! *"may split a tree into two"* and then needs a **replacement-edge search** —
//! which is where every constant factor in the batch-dynamic literature lives,
//! and which a union-find cannot do at any price.
//!
//! But the expensive machinery is only needed **when a deletion actually
//! disconnects something**. R-022's own framing says *"most digging does not
//! alter connectivity"*. If the same is true of filling, then the design space
//! is not "implement dynamic connectivity" but "**detect the rare split
//! cheaply, and recompute only then**" — a different and much smaller ticket.
//!
//! So: apply brush fills to a field with real cave structure and count how many
//! change the air component count. That number decides the shape.
//!
//! # Why the component count is recomputed rather than maintained
//!
//! Maintaining it is the thing under question. Rebuilding is `O(n)` and this is
//! a measurement, not a hot path — and a rebuilt count cannot inherit a bug from
//! the structure whose necessity it is being used to judge.

// Exact equality throughout: the question is whether a value changed at all.
#![allow(
    clippy::float_cmp,
    reason = "the question is whether a value changed at all"
)]

mod common;

use std::fmt::Write as _;

use isomesh::connectivity::Air;
use isomesh::fields::{ReferenceField, noise_cavity};
use isomesh::{RuntimeShape3, Sdf, Shape3};

/// Samples per axis.
const RESOLUTIONS: [u32; 3] = [33, 49, 65];

/// Brush radius in samples, fixed across resolutions.
const BRUSH_RADIUS: f64 = 4.0;

/// Fills applied per resolution.
const FILLS: usize = 200;

/// The generator, so the fill sequence is the same on every machine.
struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    fn below(&mut self, n: u32) -> u32 {
        (self.next_u64() >> 33) as u32 % n
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!(
        "{:>5} {:>9} {:>7} {:>8} {:>9} {:>10} {:>11} {:>12}",
        "n", "air", "fills", "changed", "split", "emptied", "changed %", "air left %"
    );
    let mut csv = String::from(
        "samples_per_axis,samples,air_samples,components_before,fills,\
         fills_changing_components,fills_splitting,fills_emptying,mean_dirty,air_left_pct\n",
    );

    for n in RESOLUTIONS {
        let Ok(shape) = RuntimeShape3::new([n; 3]) else {
            continue;
        };
        let field = noise_cavity::<f64>();
        let (lo, hi) = field.domain();
        let h = (hi[0] - lo[0]) / f64::from(n - 1);

        // `noise_cavity` is the reference field with real cave structure — it
        // exists because none of the others produces an interior ambiguity
        // (M-208), and that same roughness is what gives the air region a
        // topology worth disconnecting.
        let mut values = Vec::with_capacity(shape.element_count());
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    values.push(field.sample([
                        lo[0] + h * f64::from(x),
                        lo[1] + h * f64::from(y),
                        lo[2] + h * f64::from(z),
                    ]));
                }
            }
        }

        let Ok((mut air, _)) = Air::build(&values, &shape) else {
            continue;
        };
        let air_samples = air.air_samples();
        let components_before = air.components();

        let mut rng = Lcg(0x0000_A022_B000_0001);
        let (mut changed, mut split, mut emptied, mut dirty_total) = (0u64, 0u64, 0u64, 0u64);
        let mut previous = components_before;
        let mut filled = values.clone();

        for _ in 0..FILLS {
            // A brush centred on a random sample. Filling means forcing solid.
            let centre = [rng.below(n), rng.below(n), rng.below(n)];
            let r = BRUSH_RADIUS.ceil() as i64;
            let mut dirty = 0u64;
            for dz in -r..=r {
                for dy in -r..=r {
                    for dx in -r..=r {
                        if (dx * dx + dy * dy + dz * dz) as f64 > BRUSH_RADIUS * BRUSH_RADIUS {
                            continue;
                        }
                        let p = [
                            i64::from(centre[0]) + dx,
                            i64::from(centre[1]) + dy,
                            i64::from(centre[2]) + dz,
                        ];
                        if p.iter().any(|c| *c < 0 || *c >= i64::from(n)) {
                            continue;
                        }
                        let i = (p[2] as usize * n as usize + p[1] as usize) * n as usize
                            + p[0] as usize;
                        if let Some(v) = filled.get_mut(i) {
                            // Solid. Only counts as dirty if it was air.
                            if *v >= 0.0 {
                                *v = -1.0;
                                dirty += 1;
                            }
                        }
                    }
                }
            }
            dirty_total += dirty;

            // Rebuilt, not maintained — see the module docs.
            let Ok((mut rebuilt, _)) = Air::build(&filled, &shape) else {
                continue;
            };
            let now = rebuilt.components();
            if now != previous {
                changed += 1;
                if now > previous {
                    split += 1;
                } else {
                    emptied += 1;
                }
            }
            previous = now;
        }

        let pct = changed as f64 / FILLS as f64 * 100.0;
        // How much air the run consumed. The fills accumulate, so a run that
        // eats most of the domain is measuring a shrinking field rather than a
        // steady one -- and that confound has to be visible, not inferred.
        let air_left = match Air::build(&filled, &shape) {
            Ok((a, _)) => a.air_samples(),
            Err(_) => 0,
        };
        let left_pct = if air_samples == 0 {
            0.0
        } else {
            air_left as f64 / air_samples as f64 * 100.0
        };
        println!(
            "{n:>5} {air_samples:>9} {FILLS:>7} {changed:>8} {split:>9} {emptied:>10} \
             {pct:>10.1}% {left_pct:>11.1}%"
        );
        let _ = writeln!(
            csv,
            "{n},{},{air_samples},{components_before},{FILLS},{changed},{split},{emptied},{:.2},{left_pct:.2}",
            shape.element_count(),
            dirty_total as f64 / FILLS as f64
        );
    }

    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements/fill_disconnect.csv");
    match std::fs::write(&path, &csv) {
        Ok(()) => println!("\nwrote {}", path.display()),
        Err(e) => println!("\n::error:: {}: {e}", path.display()),
    }
}
