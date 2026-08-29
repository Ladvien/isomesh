//! **P-123 — where `M-318`'s 45x goes, decomposed into three terms.**
//!
//! Ticket: R-027a. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p123
//! ```
//!
//! Writes `docs/experiments/p-123.csv`.
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
    let prereg = isomesh::experiment!("P-123");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-123 harness not written yet");
    });
}
