//! **P-26 — does repairing air connectivity after a FILL cost the shed volume,
//! or the component?**
//!
//! Ticket: R-022b. Pre-registered in a commit before this one, as P-25 and then
//! as P-26 once ✗26 killed P-25's mechanism clause.
//!
//! ```bash
//! cargo bench --bench experiment_p26
//! ```
//!
//! Writes `docs/experiments/p-26.csv`.
//!
//! # Two fixtures, and the prediction is different for each
//!
//! This is the part that distinguishes P-26 from P-25. Lockstep search stops
//! when all but one frontier exhausts, so its cost is the **second-largest**
//! piece. M-320 measured the smaller side of a split at **one voxel at the
//! median** — but that is a property of the **edit distribution**, not of the
//! structure, and a harness that only ran that distribution would be measuring
//! M-320's fixture rather than the thing built.
//!
//! - **`distribution`** — `noise_cavity`, 200 random radius-4 brush fills, the
//!   same field, brush, seed and sequence M-319 and M-320 used. **Predicted
//!   flat** as the lattice grows.
//!
//! - **`bisect`** — two equal caverns joined by a one-voxel tunnel; fill the
//!   tunnel's midpoint. Both frontiers are then huge and the search walks until
//!   one exhausts. **Predicted to grow with `n`**, and that growth is *not*
//!   falsifying. A structure that came out flat on both would mean this fixture
//!   is not adversarial and needs rebuilding — a fixture failure, not a result.
//!
//! That edit is not exotic. Sealing a passage between two spaces *is* the
//! sealed-volume mechanic, which is the thing this layer exists to support.
//!
//! # The control
//!
//! `rebuild_visited` is what a from-scratch `Air::build` scans: every sample,
//! `n³`. It is the cost the incremental path exists to avoid, and it is the
//! only thing the `visited` column is interesting *against*.

mod common;

use std::fmt::Write as _;

use isomesh::connectivity::Air;
use isomesh::fields::{ReferenceField, noise_cavity};
use isomesh::{RuntimeShape3, Sdf, Shape3};

/// Samples per axis. Matches M-319 and M-320 so the tables can be read together.
const RESOLUTIONS: [u32; 3] = [33, 49, 65];

/// Brush radius in samples, fixed across resolutions. Matches M-320.
const BRUSH_RADIUS: f64 = 4.0;

/// Fills applied per resolution on the distribution fixture. Matches M-320.
const FILLS: usize = 200;

/// The generator, so the fill sequence is the same on every machine — and the
/// same sequence M-319 and M-320 used.
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

/// One row of the table.
struct Row {
    fixture: &'static str,
    n: u32,
    fills: u64,
    dirty: u64,
    seeds: u64,
    visited: u64,
    splits: u64,
    shed: u64,
    vanished: u64,
    rebuild_visited: u64,
}

/// The measured distribution: random brush fills into a cave field.
fn distribution(n: u32) -> Option<Row> {
    let shape = RuntimeShape3::new([n; 3]).ok()?;
    let field = noise_cavity::<f64>();
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(n - 1);

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

    let (mut air, _) = Air::build(&values, &shape).ok()?;
    let mut rng = Lcg(0x0000_A022_B000_0001);
    let mut row = Row {
        fixture: "distribution",
        n,
        fills: FILLS as u64,
        dirty: 0,
        seeds: 0,
        visited: 0,
        splits: 0,
        shed: 0,
        vanished: 0,
        rebuild_visited: shape.element_count() as u64,
    };

    let r = BRUSH_RADIUS.ceil() as i64;
    let mut brush: Vec<[u32; 3]> = Vec::new();
    for _ in 0..FILLS {
        let centre = [rng.below(n), rng.below(n), rng.below(n)];
        brush.clear();
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
                    brush.push([p[0] as u32, p[1] as u32, p[2] as u32]);
                }
            }
        }
        // Synchronous: the budget is a separate concern from the cost, and a
        // partial repair would make `visited` mean two different things.
        let f = air.fill(&brush, || true);
        row.dirty += f.dirty;
        row.seeds += f.seeds;
        row.visited += f.visited;
        row.splits += f.splits;
        row.shed += f.shed;
        row.vanished += f.vanished;
    }
    Some(row)
}

/// The adversarial fixture: two equal caverns joined by one tunnel, cut in half.
///
/// This is what M-320's distribution does not contain, and it is the shape the
/// sealed-volume mechanic actually produces.
fn bisect(n: u32) -> Option<Row> {
    let shape = RuntimeShape3::new([n; 3]).ok()?;
    let count = shape.element_count();
    // Solid everywhere; the caverns are dug rather than sampled, so their sizes
    // are exact and equal by construction rather than by luck.
    let values = vec![-1.0_f64; count];
    let (mut air, _) = Air::build(&values, &shape).ok()?;

    let mid = n / 2;
    let mut cavern: Vec<[u32; 3]> = Vec::new();
    // Two slabs either side of the midplane, each a fixed fraction of the
    // lattice, so both frontiers grow with n and neither exhausts early.
    for z in 1..n - 1 {
        for y in 1..n - 1 {
            for x in 1..n - 1 {
                if x != mid {
                    cavern.push([x, y, z]);
                }
            }
        }
    }
    // The tunnel: a single voxel on the midplane joining the two slabs.
    cavern.push([mid, mid, mid]);
    air.dig(&cavern, || true);
    if air.components() != 1 {
        return None;
    }

    let mut row = Row {
        fixture: "bisect",
        n,
        fills: 1,
        dirty: 0,
        seeds: 0,
        visited: 0,
        splits: 0,
        shed: 0,
        vanished: 0,
        rebuild_visited: count as u64,
    };
    let f = air.fill(&[[mid, mid, mid]], || true);
    row.dirty += f.dirty;
    row.seeds += f.seeds;
    row.visited += f.visited;
    row.splits += f.splits;
    row.shed += f.shed;
    row.vanished += f.vanished;
    Some(row)
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-26");
    println!("{}\n", prereg.id);

    println!(
        "{:>13} {:>5} {:>7} {:>9} {:>8} {:>9} {:>7} {:>6} {:>9} {:>10}",
        "fixture",
        "n",
        "fills",
        "dirty",
        "seeds",
        "visited",
        "splits",
        "shed",
        "vanished",
        "rebuild"
    );
    let mut csv = String::from(
        "samples_per_axis,fixture,fills,dirty_samples,seeds,visited,splits,\
         shed_components,vanished_components,rebuild_visited\n",
    );

    let mut rows: Vec<Row> = Vec::new();
    for n in RESOLUTIONS {
        rows.extend(distribution(n));
    }
    for n in RESOLUTIONS {
        rows.extend(bisect(n));
    }

    for row in &rows {
        println!(
            "{:>13} {:>5} {:>7} {:>9} {:>8} {:>9} {:>7} {:>6} {:>9} {:>10}",
            row.fixture,
            row.n,
            row.fills,
            row.dirty,
            row.seeds,
            row.visited,
            row.splits,
            row.shed,
            row.vanished,
            row.rebuild_visited
        );
        let _ = writeln!(
            csv,
            "{},{},{},{},{},{},{},{},{},{}",
            row.n,
            row.fixture,
            row.fills,
            row.dirty,
            row.seeds,
            row.visited,
            row.splits,
            row.shed,
            row.vanished,
            row.rebuild_visited
        );
    }

    // The verdict, stated in the terms the registration used.
    let read = |fixture: &str| -> Vec<u64> {
        rows.iter()
            .filter(|r| r.fixture == fixture)
            .map(|r| r.visited)
            .collect()
    };
    let dist = read("distribution");
    let bis = read("bisect");
    println!("\ndistribution visited: {dist:?}");
    println!("bisect       visited: {bis:?}");
    if let (Some(&d0), Some(&dn)) = (dist.first(), dist.last()) {
        println!(
            "distribution grew {:.2}x while the lattice grew {:.1}x",
            if d0 == 0 { 0.0 } else { dn as f64 / d0 as f64 },
            (65.0_f64 / 33.0).powi(3)
        );
    }
    if let (Some(&b0), Some(&bn)) = (bis.first(), bis.last()) {
        println!(
            "bisect       grew {:.2}x -- PREDICTED, and not falsifying",
            if b0 == 0 { 0.0 } else { bn as f64 / b0 as f64 }
        );
    }

    // The raw column is not the claim. P-23 could say "flat" because it held the
    // brush fixed and the edit really was constant; here the same brush removes
    // MORE air as the lattice grows, because there is more air to remove. So the
    // quantity P-26 is about -- work per unit of edit -- is visited/seed, and
    // that is what has to be flat.
    println!(
        "\n{:>13} {:>5} {:>9} {:>8} {:>12} {:>14}",
        "fixture", "n", "visited", "seeds", "per seed", "vs 200 rebuilds"
    );
    for row in &rows {
        let per_seed = if row.seeds == 0 {
            0.0
        } else {
            row.visited as f64 / row.seeds as f64
        };
        let against = row.rebuild_visited.saturating_mul(row.fills);
        let ratio = if row.visited == 0 {
            0.0
        } else {
            against as f64 / row.visited as f64
        };
        println!(
            "{:>13} {:>5} {:>9} {:>8} {per_seed:>12.2} {ratio:>13.1}x",
            row.fixture, row.n, row.visited, row.seeds
        );
    }

    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/experiments/p-26.csv");
    match std::fs::write(&path, &csv) {
        Ok(()) => println!("\nwrote {}", path.display()),
        Err(e) => println!("\n::error:: {}: {e}", path.display()),
    }
}
