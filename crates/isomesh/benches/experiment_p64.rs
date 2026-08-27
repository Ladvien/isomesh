//! **P-64 - bounded model checking over the case tables.**
//!
//! Ticket: R-062. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p64
//! ```
//!
//! Writes `docs/experiments/p-64.csv`.
//!
//! # This bench runs the solver rather than restating its answer
//!
//! The proofs live in `marching_cubes::proofs` behind `cfg(kani)`. This harness
//! shells out to `cargo kani` once per proof, parses the check count, the failure
//! count, the verdict and the solver time out of Kani's own output, and records
//! those. **A CSV that named numbers a human had read off a terminal would be a
//! second copy of the result**, and the whole idiom here is that a dataset names
//! a commit and a machine and comes from the instrument.
//!
//! # Why the check count is a column and not a footnote
//!
//! `VERIFICATION:- SUCCESSFUL` over **zero** checks is `M-44`'s vacuous zero in
//! formal clothing, and it prints the same word as a real proof. The
//! registration made the check count a required record for exactly that reason,
//! and `the_properties_can_fail` is the other half: a `#[kani::should_panic]`
//! harness that corrupts one table entry and requires the assertions to fire. If
//! a property ever drifts into a tautology, that harness fails while the three
//! real ones keep printing `SUCCESSFUL`.
//!
//! # Requires Kani
//!
//! `cargo install --locked kani-verifier && cargo kani setup`. If `cargo kani`
//! is absent the bench **fails loudly** rather than recording a skip - a row
//! saying "not run" in a file whose purpose is to say what the solver proved is
//! worse than no row.

#![allow(clippy::cast_precision_loss)]

mod common;

use std::process::Command;
use std::time::Instant;

/// The proofs, in the order the CSV lists them.
///
/// `interior_rule` is what the row is about: `CASES` is the mask-zero table
/// `extract` reads with [`isomesh::marching_cubes::InteriorAmbiguity::Ignore`],
/// and the 64-mask table is every resolution the interior rule can select.
const HARNESSES: [(&str, &str, &str, u64); 3] = [
    (
        "shipped_case_table_is_indexable_for_every_sign_pattern",
        "shape, nameable, edge-is-cut, non-degenerate",
        "off",
        256,
    ),
    (
        "every_case_and_mask_triangulation_is_indexable",
        "shape, nameable, edge-is-cut, non-degenerate",
        "on",
        256 * 64,
    ),
    (
        "the_properties_can_fail",
        "control: a corrupted entry must be caught",
        "off",
        256,
    ),
];

/// What one `cargo kani --harness` run reported.
struct Verdict {
    checks: u64,
    failed: u64,
    status: String,
    solver_seconds: f64,
    wall_seconds: f64,
}

fn run_harness(name: &str) -> Verdict {
    let started = Instant::now();
    let out = Command::new("cargo")
        .args(["kani", "--harness", name])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("cargo kani must be installed: cargo install --locked kani-verifier");
    let wall_seconds = started.elapsed().as_secs_f64();
    let text = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );

    // `** N of M failed (K unreachable)`
    let mut checks = 0u64;
    let mut failed = 0u64;
    for line in text.lines() {
        let t = line.trim();
        let Some(rest) = t.strip_prefix("** ") else {
            continue;
        };
        if let Some((n, rest)) = rest.split_once(" of ")
            && let Some((m, _)) = rest.split_once(" failed")
            && let (Ok(n), Ok(m)) = (n.parse::<u64>(), m.parse::<u64>())
        {
            failed = n;
            checks = m;
        }
    }

    let status = text
        .lines()
        .find_map(|l| l.trim().strip_prefix("VERIFICATION:- "))
        .unwrap_or("NO VERDICT")
        .trim()
        .to_string();

    let solver_seconds = text
        .lines()
        .find_map(|l| {
            l.trim()
                .strip_prefix("Verification Time: ")
                .and_then(|v| v.trim_end_matches('s').parse::<f64>().ok())
        })
        .unwrap_or(f64::NAN);

    Verdict {
        checks,
        failed,
        status,
        solver_seconds,
        wall_seconds,
    }
}

fn kani_version() -> String {
    let out = Command::new("cargo")
        .args(["kani", "--version"])
        .output()
        .expect("cargo kani --version");
    String::from_utf8_lossy(&out.stdout).trim().to_string()
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-64");
    let version = kani_version();
    println!("{version}\n");

    println!(
        "{:>10} {:>8} {:>8} {:>10} {:>10}  status / harness",
        "interior", "checks", "failed", "solver s", "wall s"
    );

    let mut results: Vec<(&str, &str, &str, u64, Verdict)> = Vec::new();
    for (name, property, interior, patterns) in HARNESSES {
        let v = run_harness(name);
        println!(
            "{:>10} {:>8} {:>8} {:>10.4} {:>10.2}  {} / {name}",
            interior, v.checks, v.failed, v.solver_seconds, v.wall_seconds, v.status
        );
        results.push((name, property, interior, patterns, v));
    }

    // ── controls ─────────────────────────────────────────────────────────────
    for (name, _, _, _, v) in &results {
        assert!(
            v.checks > 0,
            "VOID: {name} reported ZERO checks. 'VERIFICATION:- SUCCESSFUL' over an empty check \
             set is M-44's vacuous zero in formal clothing and prints the same word as a proof"
        );
    }
    let control = results
        .iter()
        .find(|(n, _, _, _, _)| *n == "the_properties_can_fail")
        .expect("the sabotage control must be in the run");
    assert!(
        control.4.failed > 0,
        "VOID: the sabotage control caught nothing, so the four properties cannot fail and the \
         three SUCCESSFUL verdicts beside it are tautologies"
    );
    assert!(
        control.4.status.starts_with("SUCCESSFUL"),
        "the sabotage control must end SUCCESSFUL under #[kani::should_panic]; it said {}",
        control.4.status
    );

    // ── verdict ──────────────────────────────────────────────────────────────
    let c1 = results
        .iter()
        .find(|(n, _, _, _, _)| *n == "shipped_case_table_is_indexable_for_every_sign_pattern")
        .is_some_and(|(_, _, _, _, v)| {
            v.failed == 0 && v.status == "SUCCESSFUL" && v.wall_seconds < 600.0
        });
    // C2 is "it finds nothing the existing suite does not cover", which is
    // exactly "no real harness reported a failure". The control is excluded: it
    // is *supposed* to fail, and counting it would make C2 unfalsifiable in the
    // wrong direction.
    let c2 = results
        .iter()
        .filter(|(n, _, _, _, _)| *n != "the_properties_can_fail")
        .all(|(_, _, _, _, v)| v.failed == 0);
    let c3 = results
        .iter()
        .find(|(n, _, _, _, _)| *n == "every_case_and_mask_triangulation_is_indexable")
        .is_some_and(|(_, _, _, _, v)| {
            v.failed == 0 && v.status == "SUCCESSFUL" && v.wall_seconds < 1800.0
        });

    println!(
        "\nC1 all four properties over 256 patterns, interior off, under 10 min -> {}",
        if c1 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C2 no reachable violation the suite does not already cover -> {}",
        if c2 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C3 interior rule on -- 256 x 64 masks -- under 30 min -> {}",
        if c3 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "\nSCOPE, per the registration: neither the proofs nor this bench touch VERTEX \
         PLACEMENT. Placement stays under proptest and the 216 golden hashes. What is proved is \
         'the table cannot be indexed wrongly', not 'the mesh is correct'."
    );

    common::experiment::run(prereg, |run| {
        for (name, property, interior, patterns, v) in &results {
            run.record(&[
                ("harness", (*name).to_string()),
                ("property", (*property).to_string()),
                ("interior_rule", (*interior).to_string()),
                ("patterns", patterns.to_string()),
                ("checks", v.checks.to_string()),
                ("failed_checks", v.failed.to_string()),
                ("status", v.status.clone()),
                ("solver_seconds", format!("{:.4}", v.solver_seconds)),
                ("wall_seconds", format!("{:.2}", v.wall_seconds)),
                ("kani_version", version.clone()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
            ]);
        }
    });
}
