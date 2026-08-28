//! **P-115 — Tree-Encoded Bitmaps for a subblock-empty summary, because WAH and EWAH are foreclosed.**
//!
//! Ticket: R-115. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p115
//! ```
//!
//! Writes `docs/experiments/p-115.csv`.
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
    let prereg = isomesh::experiment!("P-115");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-115 harness not written yet");
    });
}
