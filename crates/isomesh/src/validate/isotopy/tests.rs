//! T-015's acceptance: *"the predicate is evaluated per cell and its pass rate
//! reported per field; a field engineered to violate it is correctly flagged."*

extern crate std;

use crate::fields::{ReferenceField, capped_gyroid};
use crate::{RuntimeShape3, Sdf, Shape3};

use super::{cell_is_certified, isotopy_report};

/// **The engineered violator, and it is the paper's own Figure 1.**
///
/// Plantinga & Vegter: *"we cannot have alternating signs of F at the vertices
/// of C, since F would have to increase along one edge, and decrease along the
/// other parallel edge."* So a cell with alternating corner signs is exactly
/// what the condition forbids, and the predicate must flag it.
#[test]
fn alternating_corner_signs_are_flagged() {
    // Corner `i` at `(i&1, (i>>1)&1, (i>>2)&1)`; alternate on the parity of the
    // three bits, which puts opposite signs on every cube edge.
    let mut corner = [0.0f64; 8];
    for (i, slot) in corner.iter_mut().enumerate() {
        let parity = (i & 1) ^ ((i >> 1) & 1) ^ ((i >> 2) & 1);
        *slot = if parity == 0 { 1.0 } else { -1.0 };
    }
    assert!(
        !cell_is_certified(&corner),
        "the paper's own counterexample was certified"
    );

    // A single plane through the cell is the case that must pass: every x-edge
    // difference is the same, so the gradient does not turn at all.
    let mut plane = [0.0f64; 8];
    for (i, slot) in plane.iter_mut().enumerate() {
        let x = (i & 1) as f64;
        *slot = x - 0.5;
    }
    assert!(cell_is_certified(&plane), "a plane was not certified");

    // An inactive cell passes by the first clause, whatever its gradient does.
    let uniform = [1.0f64, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0];
    assert!(cell_is_certified(&uniform));
}

/// A face-ambiguous cell — the configuration the whole A-002 series exists for.
///
/// Two diagonally opposite corners inside on one face, and the ambiguity is
/// precisely that the trilinear surface has two possible topologies there. So
/// the certificate must **refuse** it: a cell where the answer depends on a
/// tie-break is the definition of one that is not topologically determined by
/// the corners alone.
#[test]
fn a_face_ambiguous_cell_is_not_certified() {
    // Corners 0 and 3 inside, 1 and 2 outside, on the z = 0 face; everything on
    // z = 1 outside.
    let corner = [-1.0f64, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0];
    assert!(
        !cell_is_certified(&corner),
        "an ambiguous face was certified"
    );
}

/// **The measurement (M-264).** Pass rate per field, across resolutions.
///
/// The number the ticket asks for, and the one that makes the certificate worth
/// having: a field whose rate climbs toward 1 with resolution is being resolved,
/// and one whose rate stalls has a feature the grid cannot certify at any
/// spacing.
#[test]
fn the_certified_fraction_per_field_and_resolution() {
    std::println!(
        "{:<16} {:>7} {:>9} {:>11} {:>13}",
        "field",
        "samples",
        "active",
        "uncertified",
        "certified %"
    );

    let mut csv = std::string::String::from("field,samples,active,uncertified,certified\n");
    let mut rows: std::vec::Vec<(&str, u32, u64, u64)> = std::vec::Vec::new();

    crate::for_each_reference_field!(f64, |name, field| {
        // `if`-free: every field is measured. `for_each_reference_field!`
        // expands inline blocks, so a `return` would end this test early and
        // pass having covered one field (M-253).
        for samples in [17u32, 33, 65] {
            let (lo, hi) = field.domain();
            let h = (hi[0] - lo[0]) / f64::from(samples - 1);
            let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

            let mut grid = std::vec::Vec::with_capacity(shape.element_count());
            for z in 0..samples {
                for y in 0..samples {
                    for x in 0..samples {
                        grid.push(field.sample([
                            lo[0] + h * f64::from(x),
                            lo[1] + h * f64::from(y),
                            lo[2] + h * f64::from(z),
                        ]));
                    }
                }
            }

            let report = isotopy_report(&grid, &shape).expect("report");
            let pct = 100.0 * report.certified_fraction();
            std::println!(
                "{name:<16} {samples:>7} {:>9} {:>11} {pct:>12.3}%",
                report.active_cells,
                report.uncertified
            );
            csv.push_str(&std::format!(
                "{name},{samples},{},{},{:.6}\n",
                report.active_cells,
                report.uncertified,
                report.certified_fraction()
            ));
            rows.push((name, samples, report.active_cells, report.uncertified));

            assert!(
                report.active_cells > 0,
                "{name} at {samples}³ has no active cells to certify"
            );
            assert_eq!(
                report.cells,
                u64::from(samples - 1).pow(3),
                "{name}: cells counted wrong"
            );
            assert_eq!(
                report.certified + report.uncertified,
                report.active_cells,
                "{name}: an active cell was neither certified nor not"
            );
        }
    });

    // **How the uncertified set scales, which says what it *is*.** Halving `h`
    // quadruples the count of cells meeting a surface, because a surface is
    // two-dimensional. If the uncertified cells merely doubled, they are meeting
    // a **curve** — a sharp edge, a seam, a crease — and no resolution will ever
    // certify them, because the feature is genuinely not smooth. If they
    // quadrupled, the field is simply unresolved and finer sampling will fix it.
    // Reported as an **effective dimension** rather than a band: two doublings
    // multiply a `d`-dimensional set's cell count by `2^(2d)`, so
    // `d = log2(growth) / 2`. That reads directly instead of being compared
    // against a threshold somebody has to remember the meaning of.
    std::println!("\n17³ → 65³ growth, as effective dimension:");
    std::println!(
        "{:<16} {:>9} {:>14} {:>24}",
        "field",
        "active d",
        "uncertified d",
        "reading"
    );
    let mut names: std::vec::Vec<&str> = rows.iter().map(|r| r.0).collect();
    names.dedup();
    for name in names {
        let mine: std::vec::Vec<&(&str, u32, u64, u64)> =
            rows.iter().filter(|r| r.0 == name).collect();
        if mine.len() < 3 || mine[2].3 == 0 {
            std::println!(
                "{name:<16} {:>9} {:>14} {:>24}",
                "-",
                "-",
                "fully certified"
            );
            continue;
        }
        let dim = |from: u64, to: u64| (to as f64 / from as f64).log2() / 2.0;
        let active_d = dim(mine[0].2, mine[2].2);
        let bad_d = dim(mine[0].3, mine[2].3);
        let reading = if bad_d < 1.6 {
            "a feature curve"
        } else {
            "an under-resolved area"
        };
        std::println!("{name:<16} {active_d:>9.2} {bad_d:>14.2} {reading:>24}");

        // **The assertion that makes this a finding rather than a table.** The
        // active set is two-dimensional because a surface is. If the uncertified
        // set were *also* two-dimensional, the certificate would be failing
        // everywhere and would say nothing about the field; that it is strictly
        // lower-dimensional is what identifies those cells as a **feature** —
        // an edge, a seam, a crease — that no spacing will smooth away.
        assert!(
            active_d > 1.8,
            "{name}: the active set grew as dimension {active_d:.2}, not ~2 — \
             the measurement is not measuring what it thinks"
        );
        assert!(
            bad_d < active_d - 0.3,
            "{name}: the uncertified set grew as dimension {bad_d:.2} against \
             the active set's {active_d:.2} — the certificate is failing across \
             the surface rather than along a feature"
        );
    }

    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements/isotopy.csv");
    if let Some(dir) = path.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let _ = std::fs::write(&path, csv);
    std::println!("\nwrote {}", path.display());
}

/// Every cell the certificate refuses has an axis whose derivative changes sign.
///
/// The mechanism, asserted rather than assumed — this is what distinguishes the
/// condition from a generic quality heuristic, and a refusal with no sign change
/// anywhere would mean the interval arithmetic is wrong rather than the cell.
#[test]
fn a_refusal_always_has_a_turning_derivative() {
    // `CappedGyroid` rather than `Gyroid`: it is the one with a `ReferenceField`
    // domain, and it is the reference field with the most turning in it, which
    // is what this test needs a supply of.
    let field = capped_gyroid::<f64>();
    let samples = 25u32;
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

    let mut checked = 0usize;
    let size = shape.size();
    let mut grid = std::vec::Vec::with_capacity(shape.element_count());
    for z in 0..samples {
        for y in 0..samples {
            for x in 0..samples {
                grid.push(field.sample([
                    lo[0] + h * f64::from(x),
                    lo[1] + h * f64::from(y),
                    lo[2] + h * f64::from(z),
                ]));
            }
        }
    }

    for z in 0..size[2] - 1 {
        for y in 0..size[1] - 1 {
            for x in 0..size[0] - 1 {
                let at = |dx: u32, dy: u32, dz: u32| {
                    grid[(((z + dz) * size[1] + (y + dy)) * size[0] + (x + dx)) as usize]
                };
                let mut corner = [0.0f64; 8];
                for (i, slot) in corner.iter_mut().enumerate() {
                    let i = i as u32;
                    *slot = at(i & 1, (i >> 1) & 1, (i >> 2) & 1);
                }
                if cell_is_certified(&corner) {
                    continue;
                }
                checked += 1;
                let turns = |pairs: &[[usize; 2]; 4]| {
                    let [lo, hi] = super::partial_range(&corner, pairs);
                    lo <= 0.0 && hi >= 0.0
                };
                assert!(
                    turns(&super::X_PAIRS) || turns(&super::Y_PAIRS) || turns(&super::Z_PAIRS),
                    "cell ({x},{y},{z}) was refused with no derivative changing sign"
                );
            }
        }
    }
    std::println!("measured: {checked} refusals, every one with a turning derivative");
    assert!(
        checked > 0,
        "the capped gyroid produced no refusals to examine"
    );
}

/// Malformed input is refused.
#[test]
fn it_refuses_what_it_cannot_certify() {
    let tiny = RuntimeShape3::new([1, 4, 4]).expect("valid shape");
    let v = std::vec![0.0f64; tiny.element_count()];
    assert!(isotopy_report(&v, &tiny).is_err());

    let shape = RuntimeShape3::new([4; 3]).expect("valid shape");
    assert!(isotopy_report(&[0.0f64; 3], &shape).is_err());
}
