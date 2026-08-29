//! **P-126 — `O-12`'s remaining half — the dual vertex link at 2^27, as a nightly gate.**
//!
//! Ticket: R-072. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p126
//! ```
//!
//! Writes `docs/experiments/p-126.csv`.
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
    let prereg = isomesh::experiment!("P-126");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-126 harness not written yet");
    });
}
