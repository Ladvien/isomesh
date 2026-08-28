//! **P-104 — the interleaved layout for the active-cell bitmap, which packs samples and not cells.**
//!
//! Ticket: R-104. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p104
//! ```
//!
//! Writes `docs/experiments/p-104.csv`.
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
    let prereg = isomesh::experiment!("P-104");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-104 harness not written yet");
    });
}
