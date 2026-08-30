//! **P-109 — Elias–Fano for the edge→vertex structure, which is a dense flat vec and not a map.**
//!
//! Ticket: R-109. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p109
//! ```
//!
//! Writes `docs/experiments/p-109.csv`.
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
    let prereg = isomesh::experiment!("P-109");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-109 harness not written yet");
    });
}
