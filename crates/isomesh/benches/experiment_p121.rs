//! **P-121 — what fraction of extraction is bit work — runs first, gates Group A.**
//!
//! Ticket: R-121. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p121
//! ```
//!
//! Writes `docs/experiments/p-121.csv`.
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
    let prereg = isomesh::experiment!("P-121");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-121 harness not written yet");
    });
}
