//! Running a pre-registered experiment, and writing down what it did.
//!
//! Ticket: R-000. The compile-time half lives in `isomesh::experiment`, which
//! refuses to build an experiment whose `P-` id is not registered. This is the
//! run-time half: it emits one CSV per experiment stamped with the git SHA, the
//! machine and the time, and prints the `FINDINGS.md` stanza ready to paste.
//!
//! # Why the columns are checked against the registration
//!
//! A prediction that names its records and then quietly drops one is a
//! prediction that cannot be falsified on that metric. [`Run::record`] takes
//! exactly the keys `Preregistration::records` lists and panics on any other
//! set — so a metric abandoned mid-experiment is a failure rather than a
//! silence. That is the same rule as F-002's one-sidedness in a different
//! costume: the instrument has to be able to report the bad news.
//!
//! # Provenance is a comment, not a column
//!
//! The `#` header carries SHA, machine and timestamp. `scripts/regress.sh`
//! already skips `#` lines (`rows_of`), so an experiment CSV can be diffed by
//! the same reader without a second format.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::PathBuf;
use std::process::Command;

use isomesh::experiment::Preregistration;

/// Rows accumulated by one experiment.
pub(crate) struct Run {
    prereg: &'static Preregistration,
    rows: Vec<BTreeMap<&'static str, String>>,
}

impl Run {
    /// Add one row.
    ///
    /// # Panics
    ///
    /// If the keys are not exactly `Preregistration::records`. Deliberate: a
    /// missing column is a metric that was predicted and then not measured, and
    /// finding that out from a silent CSV is how a falsified hypothesis survives.
    pub(crate) fn record(&mut self, values: &[(&'static str, String)]) {
        let mut row = BTreeMap::new();
        for (k, v) in values {
            assert!(
                self.prereg.records.contains(k),
                "{}: `{k}` is not one of the records this experiment registered: {:?}",
                self.prereg.id,
                self.prereg.records
            );
            row.insert(*k, v.clone());
        }
        for expected in self.prereg.records {
            assert!(
                row.contains_key(expected),
                "{}: row is missing `{expected}`, which the pre-registration \
                 promised to record",
                self.prereg.id
            );
        }
        self.rows.push(row);
    }
}

/// Short output of a command, or `"unknown"`.
///
/// Never a hard error: an experiment that produced real numbers should not be
/// discarded because the machine has no `git`. The word `unknown` in a
/// provenance header is honest and visible.
fn ask(program: &str, args: &[&str]) -> String {
    Command::new(program)
        .args(args)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| String::from("unknown"))
}

/// Run a pre-registered experiment, write its CSV, print its stanza.
///
/// The `Preregistration` can only come from `isomesh::experiment!`, which is a
/// compile error for an unregistered id — so there is no way to reach this
/// function without having registered first.
pub(crate) fn run(prereg: &'static Preregistration, body: impl FnOnce(&mut Run)) {
    let mut run = Run {
        prereg,
        rows: Vec::new(),
    };

    println!("{} — {}", prereg.id, prereg.ticket);
    println!("  H:            {}", prereg.hypothesis);
    println!("  falsified by: {}", prereg.falsified_by);
    println!();

    body(&mut run);

    let sha = ask("git", &["rev-parse", "--short", "HEAD"]);
    let dirty = ask("git", &["status", "--porcelain"]);
    let when = ask("date", &["-u", "+%Y-%m-%dT%H:%M:%SZ"]);
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let machine = Command::new(root.join("scripts/machine.sh"))
        .arg("--slug")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| String::from("unknown"));

    let mut csv = String::new();
    let _ = writeln!(csv, "# experiment {} ({})", prereg.id, prereg.ticket);
    let _ = writeln!(csv, "# hypothesis: {}", prereg.hypothesis);
    let _ = writeln!(csv, "# falsified by: {}", prereg.falsified_by);
    let _ = writeln!(
        csv,
        "# commit {sha}{} on {machine} at {when}",
        // A dirty tree means the numbers do not correspond to any commit, and
        // that has to be on the artefact rather than in someone's memory.
        if dirty == "unknown" || dirty.is_empty() {
            ""
        } else {
            " (WORKING TREE DIRTY)"
        }
    );
    let _ = writeln!(csv, "{}", prereg.records.join(","));
    for row in &run.rows {
        let line: Vec<&str> = prereg
            .records
            .iter()
            .map(|k| row.get(k).map_or("", String::as_str))
            .collect();
        let _ = writeln!(csv, "{}", line.join(","));
    }

    let dir = root.join("docs/experiments");
    let _ = std::fs::create_dir_all(&dir);
    let path = dir.join(format!("{}.csv", prereg.id.to_lowercase()));
    let _ = std::fs::write(&path, &csv);

    println!("\n{} rows → {}", run.rows.len(), path.display());
    println!("\n--- FINDINGS.md stanza, ready to paste ---");
    println!(
        "| {} | **<one line: what happened>** | H: {} | <both arms, with numbers> | \
         **HELD** / **FALSIFIED** — <why> | {}, `docs/experiments/{}.csv` |",
        prereg.id,
        prereg.hypothesis,
        prereg.ticket,
        prereg.id.to_lowercase()
    );
    println!("--- end ---");
}
