//! **P-110 — mutable rank/select for a structure written during extraction, scored against the determinism gate.**
//!
//! Ticket: R-110. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p110
//! ```
//!
//! Writes `docs/experiments/p-110.csv`.
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
    let prereg = isomesh::experiment!("P-110");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-110 harness not written yet");
    });
}
