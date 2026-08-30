//! **P-103 — the case index computed 64 cells at a time from bit-sliced sign planes.**
//!
//! Ticket: R-103. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p103
//! ```
//!
//! Writes `docs/experiments/p-103.csv`.
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
    let prereg = isomesh::experiment!("P-103");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-103 harness not written yet");
    });
}
