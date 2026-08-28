//! **P-114 — a hierarchical bitmap above the active-cell bitmap.**
//!
//! Ticket: R-114. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p114
//! ```
//!
//! Writes `docs/experiments/p-114.csv`.
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
    let prereg = isomesh::experiment!("P-114");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-114 harness not written yet");
    });
}
