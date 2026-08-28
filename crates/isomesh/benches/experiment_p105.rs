//! **P-105 — Harley–Seal carry-save popcount for a per-chunk active count the crate does not yet compute.**
//!
//! Ticket: R-105. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p105
//! ```
//!
//! Writes `docs/experiments/p-105.csv`.
//!
//! # What was missing
//!
//! TODO: filled in by the harness author.
//!
//! # SHARE
//!
//! TODO: filled in by the harness author.

mod common;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-105");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-105 harness not written yet");
    });
}
