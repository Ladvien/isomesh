//! **P-124 — the monotone-edge condition on the ambient complex.**
//!
//! Ticket: R-052. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p124
//! ```
//!
//! Writes `docs/experiments/p-124.csv`.
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
    let prereg = isomesh::experiment!("P-124");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-124 harness not written yet");
    });
}
