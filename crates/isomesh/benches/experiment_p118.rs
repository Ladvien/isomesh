//! **P-118 — Neal's superaccumulator for cross-cell float accumulation, aimed at `M-177`.**
//!
//! Ticket: R-118. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p118
//! ```
//!
//! Writes `docs/experiments/p-118.csv`.
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
    let prereg = isomesh::experiment!("P-118");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-118 harness not written yet");
    });
}
