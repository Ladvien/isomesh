//! **P-120 — array-based union-find for per-chunk labelling, answering `✗26` rather than reopening it.**
//!
//! Ticket: R-120. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p120
//! ```
//!
//! Writes `docs/experiments/p-120.csv`.
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
    let prereg = isomesh::experiment!("P-120");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-120 harness not written yet");
    });
}
