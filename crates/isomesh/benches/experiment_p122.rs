//! **P-122 — Stream VByte's control/data split, applied to the case stream.**
//!
//! Ticket: R-122. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p122
//! ```
//!
//! Writes `docs/experiments/p-122.csv`.
//!
//! # What was missing
//!
//! TODO: filled in by the harness author.
//!
//! # SHARE
//!
//! TODO: filled in by the harness author.

mod common;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-122");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-122 harness not written yet");
    });
}
