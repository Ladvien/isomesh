//! **P-125 — the pinch predicate, shipped as a `validate` report.**
//!
//! Ticket: R-053. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p125
//! ```
//!
//! Writes `docs/experiments/p-125.csv`.
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
    let prereg = isomesh::experiment!("P-125");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-125 harness not written yet");
    });
}
