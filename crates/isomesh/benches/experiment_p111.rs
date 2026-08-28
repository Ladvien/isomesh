//! **P-111 — table-driven scalar compaction, 8 cells per lookup, branchless.**
//!
//! Ticket: R-111. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p111
//! ```
//!
//! Writes `docs/experiments/p-111.csv`.
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
    let prereg = isomesh::experiment!("P-111");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-111 harness not written yet");
    });
}
