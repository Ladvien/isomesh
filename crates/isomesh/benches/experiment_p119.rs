//! **P-119 — double-buffering as the determinism mechanism.**
//!
//! Ticket: R-119. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p119
//! ```
//!
//! Writes `docs/experiments/p-119.csv`.
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
    let prereg = isomesh::experiment!("P-119");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-119 harness not written yet");
    });
}
