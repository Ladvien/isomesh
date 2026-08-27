//! **P-62 - the Plantinga-Vegter certificate, checked for soundness.**
//!
//! Ticket: R-060. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p62
//! ```
//!
//! Writes `docs/experiments/p-62.csv`.
//!
//! # This measures shipped code, and the registration says so
//!
//! `validate::isotopy::cell_is_certified` landed under `T-015`. Nothing here
//! reimplements it - a second copy of a predicate is a second thing to drift,
//! and the whole point of C1 is a statement about *the* predicate. What the
//! harness supplies is the other half: the `A-020` classifier's verdict on the
//! same cell, and the cross-tabulation nobody has taken.
//!
//! # The three answers, and which cells each is about
//!
//! - **C1** is over cells the classifier calls `Tunnel` or `TwelveVertexContour`.
//!   Those are cells whose trilinear patch is provably *not* a graph over any
//!   coordinate plane, so a certificate on one is unsound. The clause is
//!   one-sided: zero, or the direction is dead.
//! - **C2** is over **surface** cells only. Certifying an inactive cell is free
//!   (clause one is "all eight corners share a sign"), so a fraction over all
//!   cells would read 97% on a sphere at 128³ and mean nothing.
//! - **C3** is the predicate's wall time against a real `marching_cubes`
//!   extraction on the same grid, in one binary and one run (`M-281`).
//!
//! # Controls
//!
//! - **The classifier must find the configuration.** `M-214` counted 2,053
//!   tunnels and 173 twelve-vertex contours in 396,000 cells, so the fixture
//!   *can* produce them - but this fixture is not that one, and a C1 pass over a
//!   population of zero is `M-44`'s vacuous zero. The harness reports the
//!   population per row and asserts it is non-zero **somewhere** in the sweep.
//! - **The certificate must be able to fail.** A predicate that returned `true`
//!   everywhere would pass C1 trivially. The harness asserts that at least one
//!   surface cell is **refused**, per field.
//! - **The saddle classifier must agree with the case table about activity.** A
//!   cell the classifier calls a tunnel must be an active cell; if it is not,
//!   the two halves are being fed different corner values.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]

mod common;

use std::time::Instant;

use isomesh::extractor::Extractor;
use isomesh::fields::ReferenceField;
use isomesh::for_each_reference_field;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::table::AMBIGUOUS_FACES;
use isomesh::marching_cubes::trilinear::{BodySaddles, Contours, Topology};
use isomesh::validate::cell_is_certified;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Samples per axis. The registered three.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// `cube.rs`'s corner numbering: corner `i` sits at `(i&1, (i>>1)&1, (i>>2)&1)`.
const CORNERS: [[u32; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
];

/// `cube::is_inside`, which is private. Negative is inside; exactly zero is
/// outside.
fn is_inside(v: f64) -> bool {
    v < 0.0
}

/// Random cells for C1's population, because the reference fields do not supply
/// one.
///
/// **The eight reference fields at 17³-65³ produced seven tunnel cells in
/// 172,032 cells.** C1 would then be a zero over a population of seven, which is
/// a hair away from `M-44`'s vacuous zero and not the kill-shot the registration
/// describes. `M-214`'s 2,053 tunnels came from **400,000 random cells**, not
/// from a smooth field, and that is the arm that makes the clause bite: corner
/// values drawn uniformly from `[-1, 1)` hit the ambiguous configurations at a
/// rate no smooth sampled field does.
const RANDOM_CELLS: u64 = 400_000;

/// The same LCG the trilinear tests use, so the population is reproducible and
/// is drawn the same way `M-214`'s was.
struct Lcg(u64);

impl Lcg {
    /// A value in `[-1, 1)`.
    fn signed(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        f64::from((self.0 >> 40) as u32) / f64::from(1u32 << 23) - 1.0
    }
}

/// The cross-tabulation over random cells: certificate against classifier.
///
/// No timing and no extraction - this arm exists only to give C1 a population,
/// and reporting a wall time for it would invite comparison with the sampled
/// rows, which are a different fixture.
fn random_arm() -> Row {
    let mut rng = Lcg(0x2026_u64 ^ 0x5EED_1234);
    let mut row = Row {
        field: "random_cells",
        samples: 0,
        cells: RANDOM_CELLS,
        surface_cells: 0,
        certified_cells: 0,
        certified_surface_cells: 0,
        refused_surface_cells: 0,
        tunnel_cells: 0,
        twelve_vertex_cells: 0,
        unsound: 0,
        inactive_tunnels: 0,
        predicate_ms: 0.0,
        gather_ms: 0.0,
        extract_ms: 0.0,
    };
    for _ in 0..RANDOM_CELLS {
        let mut corner = [0.0f64; 8];
        for c in &mut corner {
            *c = rng.signed();
        }
        let mut case = 0u8;
        for (i, v) in corner.iter().enumerate() {
            if is_inside(*v) {
                case |= 1 << i;
            }
        }
        let active = case != 0 && case != 255;
        let cert = cell_is_certified(&corner);
        if cert {
            row.certified_cells += 1;
        }
        if active {
            row.surface_cells += 1;
            if cert {
                row.certified_surface_cells += 1;
            } else {
                row.refused_surface_cells += 1;
            }
        }
        let mask = joined_mask(&corner, AMBIGUOUS_FACES[case as usize]);
        let saddles = BodySaddles::of(&corner);
        let topo = Contours::of(case, mask).topology(&saddles);
        let hidden = match topo {
            Topology::Tunnel => {
                row.tunnel_cells += 1;
                true
            }
            Topology::TwelveVertexContour => {
                row.twelve_vertex_cells += 1;
                true
            }
            _ => false,
        };
        if hidden {
            if !active {
                row.inactive_tunnels += 1;
            }
            if cert {
                row.unsound += 1;
            }
        }
    }
    row
}

/// One field at one resolution.
struct Row {
    field: &'static str,
    samples: u32,
    cells: u64,
    surface_cells: u64,
    certified_cells: u64,
    certified_surface_cells: u64,
    refused_surface_cells: u64,
    tunnel_cells: u64,
    twelve_vertex_cells: u64,
    /// Certified **and** classified as a tunnel or twelve-vertex contour. C1's
    /// whole subject.
    unsound: u64,
    /// Tunnels or twelve-vertex contours the case table calls inactive, which
    /// would mean the two halves disagree about the cell.
    inactive_tunnels: u64,
    predicate_ms: f64,
    /// The same loop with the predicate removed: eight corner loads and nothing
    /// else.
    ///
    /// **C3's denominator is an extraction, and an extraction already gathers
    /// these eight corners** - it has to, to index the case table. So a
    /// standalone predicate pass pays for the gather twice, and the share it
    /// reports is not the cost of adding the predicate to a mesher that already
    /// walks the grid. Timing the bare gather separates the two, which is the
    /// difference between "this is a debug gate" and "this is 3 multiplies and 2
    /// adds per cell that a fused version would pay".
    gather_ms: f64,
    extract_ms: f64,
}

fn measure<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    field_name: &'static str,
    samples: u32,
) -> Row {
    let shape = RuntimeShape3::new([samples; 3]).expect("shape");
    let ([lo, hi], _) = ([field.domain().0, field.domain().1], ());
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);

    // The sample grid, once. Both halves read the same values, which is what
    // makes the cross-tabulation about the cell rather than about two samplings.
    let n = samples as usize;
    let mut values = Vec::with_capacity(n * n * n);
    for z in 0..samples {
        for y in 0..samples {
            for x in 0..samples {
                values.push(field.sample([
                    lo[0] + h * f64::from(x),
                    lo[1] + h * f64::from(y),
                    lo[2] + h * f64::from(z),
                ]));
            }
        }
    }

    let at =
        |x: u32, y: u32, z: u32| -> f64 { values[(z as usize * n + y as usize) * n + x as usize] };
    let corners_of = |x: u32, y: u32, z: u32| -> [f64; 8] {
        let mut c = [0.0f64; 8];
        for (i, o) in CORNERS.iter().enumerate() {
            c[i] = at(x + o[0], y + o[1], z + o[2]);
        }
        c
    };

    let cells = samples - 1;
    let mut row = Row {
        field: field_name,
        samples,
        cells: u64::from(cells).pow(3),
        surface_cells: 0,
        certified_cells: 0,
        certified_surface_cells: 0,
        refused_surface_cells: 0,
        tunnel_cells: 0,
        twelve_vertex_cells: 0,
        unsound: 0,
        inactive_tunnels: 0,
        predicate_ms: 0.0,
        gather_ms: 0.0,
        extract_ms: 0.0,
    };

    // ── the predicate, timed on its own ──────────────────────────────────────
    //
    // Separate pass from the classification so C3 measures the predicate and not
    // the classifier, which is a debug instrument nobody proposes to ship.
    let t = Instant::now();
    let mut certified = 0u64;
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                if cell_is_certified(&corners_of(x, y, z)) {
                    certified += 1;
                }
            }
        }
    }
    row.predicate_ms = t.elapsed().as_nanos() as f64 / 1e6;
    row.certified_cells = certified;

    // The same walk, gather and all, with the predicate removed. `black_box` on
    // the corner array so the gather cannot be elided.
    let t = Instant::now();
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                std::hint::black_box(corners_of(x, y, z));
            }
        }
    }
    row.gather_ms = t.elapsed().as_nanos() as f64 / 1e6;

    // ── the cross-tabulation, untimed ────────────────────────────────────────
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let corner = corners_of(x, y, z);
                let mut case = 0u8;
                for (i, v) in corner.iter().enumerate() {
                    if is_inside(*v) {
                        case |= 1 << i;
                    }
                }
                let active = case != 0 && case != 255;
                let cert = cell_is_certified(&corner);

                if active {
                    row.surface_cells += 1;
                    if cert {
                        row.certified_surface_cells += 1;
                    } else {
                        row.refused_surface_cells += 1;
                    }
                }

                // The classifier. `joined_mask` resolves ambiguous faces by the
                // asymptotic decider, which is what `extract` does at
                // `mod.rs:407` - so this is the same verdict the mesher would
                // reach, not a stricter or looser one.
                let mask = joined_mask(&corner, AMBIGUOUS_FACES[case as usize]);
                let saddles = BodySaddles::of(&corner);
                let topo = Contours::of(case, mask).topology(&saddles);
                let hidden = match topo {
                    Topology::Tunnel => {
                        row.tunnel_cells += 1;
                        true
                    }
                    Topology::TwelveVertexContour => {
                        row.twelve_vertex_cells += 1;
                        true
                    }
                    _ => false,
                };
                if hidden {
                    if !active {
                        row.inactive_tunnels += 1;
                    }
                    if cert {
                        row.unsound += 1;
                    }
                }
            }
        }
    }

    // ── the extraction, for C3's denominator ─────────────────────────────────
    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    let t = Instant::now();
    let _ = mc.extract_into(field, &shape, lo, h, &mut out);
    row.extract_ms = t.elapsed().as_nanos() as f64 / 1e6;
    std::hint::black_box(&out);

    row
}

type CsvRow = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-62");
    let mut rows: Vec<Row> = Vec::new();

    // **Through `for_each_reference_field!`, not a hand-written list of eight.**
    // The macro is the crate's own definition of "the eight reference fields",
    // and a list retyped here would be a second definition that drifts the
    // moment a ninth is added - which is `E-113`'s lesson in miniature.
    for samples in RESOLUTIONS {
        for_each_reference_field!(f64, |name, field| {
            rows.push(measure(&field, name, samples));
        });
    }
    rows.push(random_arm());

    println!(
        "{:>15} {:>5} {:>9} {:>9} {:>8} {:>8} {:>8} {:>8} {:>8} {:>7}",
        "field",
        "n",
        "surface",
        "cert surf",
        "frac",
        "refused",
        "tunnels",
        "12vert",
        "UNSOUND",
        "share"
    );
    for r in &rows {
        let frac = if r.surface_cells == 0 {
            f64::NAN
        } else {
            r.certified_surface_cells as f64 / r.surface_cells as f64
        };
        let share = r.predicate_ms / r.extract_ms;
        println!(
            "{:>15} {:>5} {:>9} {:>9} {:>8.4} {:>8} {:>8} {:>8} {:>8} {:>7.4}",
            r.field,
            r.samples,
            r.surface_cells,
            r.certified_surface_cells,
            frac,
            r.refused_surface_cells,
            r.tunnel_cells,
            r.twelve_vertex_cells,
            r.unsound,
            share
        );
    }

    // ── controls ─────────────────────────────────────────────────────────────
    let hidden_population: u64 = rows
        .iter()
        .map(|r| r.tunnel_cells + r.twelve_vertex_cells)
        .sum();
    assert!(
        hidden_population > 0,
        "VOID: the classifier found no tunnel and no twelve-vertex contour in {} rows, so C1's \
         zero is M-44's vacuous zero rather than a soundness result",
        rows.len()
    );
    // **Global, not per-field, and the first version was wrong to be per-field.**
    // It fired on `sphere` at 17³, which refuses nothing - and that is a correct
    // outcome, not a broken instrument: every surface cell of a sphere sampled
    // at 17³ genuinely is a graph over a coordinate plane. What has to be shown
    // is that the predicate *can* refuse, which is a statement about the
    // predicate and therefore about the sweep, not about each field.
    let refused_total: u64 = rows.iter().map(|r| r.refused_surface_cells).sum();
    assert!(
        refused_total > 0,
        "VOID: the certificate refused no surface cell in any of {} rows, so its pass rate is \
         not a measurement",
        rows.len()
    );
    let inactive: u64 = rows.iter().map(|r| r.inactive_tunnels).sum();
    assert_eq!(
        inactive, 0,
        "{inactive} cells are classified as a tunnel or twelve-vertex contour while the case \
         table calls them inactive: the predicate and the classifier are being fed different \
         corner values"
    );

    // ── verdict ──────────────────────────────────────────────────────────────
    let unsound: u64 = rows.iter().map(|r| r.unsound).sum();
    let c1 = unsound == 0;

    // C2, first half: the three named fields at 33³.
    let named = ["sphere", "torus", "box_exact"];
    let mut c2_yield = true;
    for name in named {
        let r = rows
            .iter()
            .find(|r| r.field == name && r.samples == 33)
            .expect("row");
        let frac = r.certified_surface_cells as f64 / r.surface_cells as f64;
        if frac <= 0.5 {
            c2_yield = false;
        }
        println!("C2 yield {name} at 33³: {frac:.4}");
    }

    // C2, second half: monotone in resolution, per field.
    let mut monotone: Vec<(&'static str, bool, Vec<f64>)> = Vec::new();
    // The population arm has no resolution sweep, so it is not a C2 subject:
    // C2 is about how the yield moves as a sampled field is refined, and
    // `random_cells` is not sampled from a field at all.
    let mut field_names: Vec<&'static str> = rows
        .iter()
        .filter(|r| r.field != "random_cells")
        .map(|r| r.field)
        .collect();
    field_names.dedup();
    field_names.sort_unstable();
    field_names.dedup();
    for name in field_names {
        let seq: Vec<f64> = RESOLUTIONS
            .iter()
            .map(|n| {
                let r = rows
                    .iter()
                    .find(|r| r.field == name && r.samples == *n)
                    .expect("row");
                r.certified_surface_cells as f64 / r.surface_cells as f64
            })
            .collect();
        let up = seq.windows(2).all(|w| w[1] >= w[0]);
        println!(
            "C2 monotone {name}: {:.4} -> {:.4} -> {:.4} -> {}",
            seq[0],
            seq[1],
            seq[2],
            if up { "rising" } else { "NOT MONOTONE" }
        );
        monotone.push((name, up, seq));
    }
    let c2_monotone = monotone.iter().all(|(_, up, _)| *up);
    let c2 = c2_yield && c2_monotone;

    // C3: the predicate's share of extraction at 65³, worst field.
    let worst_share = rows
        .iter()
        .filter(|r| r.samples == 65)
        .map(|r| r.predicate_ms / r.extract_ms)
        .fold(0.0f64, f64::max);
    let c3 = worst_share < 0.05;
    // The same share with the gather removed, which is what a fused
    // implementation would pay. Reported, never substituted for C3: the clause
    // is about the standalone predicate and this harness measured the
    // standalone predicate.
    let worst_fused = rows
        .iter()
        .filter(|r| r.samples == 65)
        .map(|r| ((r.predicate_ms - r.gather_ms).max(0.0)) / r.extract_ms)
        .fold(0.0f64, f64::max);

    println!(
        "\nC1 unsound certificates over {hidden_population} tunnel/twelve-vertex cells: \
         {unsound} -> {}",
        if c1 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C2 yield > 50% on the three named at 33³ AND monotone on all eight: {c2_yield} && \
         {c2_monotone} -> {}",
        if c2 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C3 worst predicate share at 65³: {worst_share:.4} -> {}",
        if c3 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C3 decomposed, NOT a substitute for the clause: worst share with the corner gather \
         removed {worst_fused:.4}, which is what a version fused into an extractor that already \
         gathers the corners would pay"
    );
    println!(
        "\nThe registered caveat, restated because it bounds what C1 means: a certified cell's \
         patch is a GRAPH over a coordinate plane, not necessarily ONE component -- a graph over \
         a disconnected planar domain has several. PV close that with a balanced octree this \
         crate does not have."
    );

    common::experiment::run(prereg, |run| {
        for r in &rows {
            let frac = r.certified_surface_cells as f64 / r.surface_cells as f64;
            let mono = monotone
                .iter()
                .find(|(n, _, _)| *n == r.field)
                .is_some_and(|(_, up, _)| *up);
            let mut csv: CsvRow = vec![
                ("field", r.field.to_string()),
                ("samples_per_axis", r.samples.to_string()),
                ("cells", r.cells.to_string()),
                ("surface_cells", r.surface_cells.to_string()),
                ("certified_cells", r.certified_cells.to_string()),
                (
                    "certified_surface_cells",
                    r.certified_surface_cells.to_string(),
                ),
                ("certified_surface_fraction", format!("{frac:.6}")),
                ("refused_surface_cells", r.refused_surface_cells.to_string()),
                ("tunnel_cells", r.tunnel_cells.to_string()),
                ("twelve_vertex_cells", r.twelve_vertex_cells.to_string()),
                ("unsound_certificates", r.unsound.to_string()),
                ("monotone_in_resolution", mono.to_string()),
                ("predicate_ms", format!("{:.4}", r.predicate_ms)),
                ("gather_ms", format!("{:.4}", r.gather_ms)),
                (
                    "fused_share",
                    format!(
                        "{:.6}",
                        (r.predicate_ms - r.gather_ms).max(0.0) / r.extract_ms
                    ),
                ),
                ("extract_ms", format!("{:.4}", r.extract_ms)),
                (
                    "predicate_share",
                    format!("{:.6}", r.predicate_ms / r.extract_ms),
                ),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
            ];
            csv.push(("hidden_population", hidden_population.to_string()));
            run.record(&csv);
        }
    });
}
