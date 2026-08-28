//! **P-117 — FMA contraction as a latent golden-hash divergence — a risk audit, not an optimisation.**
//!
//! Ticket: R-117. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p117
//! ```
//!
//! Writes `docs/experiments/p-117.csv`.
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
    let prereg = isomesh::experiment!("P-117");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-117 harness not written yet");
    });
}
