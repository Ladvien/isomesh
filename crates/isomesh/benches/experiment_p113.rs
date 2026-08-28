//! **P-113 — Roaring's density thresholds as the chunk-representation decision rule.**
//!
//! Ticket: R-113. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p113
//! ```
//!
//! Writes `docs/experiments/p-113.csv`.
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
    let prereg = isomesh::experiment!("P-113");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-113 harness not written yet");
    });
}
