//! **P-112 — count, scan, scatter, and the argument against the middle phase.**
//!
//! Ticket: R-112. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p112
//! ```
//!
//! Writes `docs/experiments/p-112.csv`.
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
    let prereg = isomesh::experiment!("P-112");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-112 harness not written yet");
    });
}
