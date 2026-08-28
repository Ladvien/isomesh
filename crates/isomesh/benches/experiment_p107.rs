//! **P-107 — a rank directory over the active-cell bitmap gives the output slot index in O(1).**
//!
//! Ticket: R-107. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p107
//! ```
//!
//! Writes `docs/experiments/p-107.csv`.
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
    let prereg = isomesh::experiment!("P-107");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-107 harness not written yet");
    });
}
