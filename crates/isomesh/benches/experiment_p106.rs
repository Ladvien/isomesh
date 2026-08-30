//! **P-106 — SWAR sign extraction and edge-crossing masks, exhaustively over all 256 patterns.**
//!
//! Ticket: R-106. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p106
//! ```
//!
//! Writes `docs/experiments/p-106.csv`.
//!
//! # What was missing
//!
//! TODO: filled in by the harness author.
//!
//! # SHARE
//!
//! TODO: filled in by the harness author.

mod common;

#[allow(clippy::unimplemented)] // skeleton until the harness commit; the allow and the stub go together
fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-106");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-106 harness not written yet");
    });
}
