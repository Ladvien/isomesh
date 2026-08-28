//! **P-108 — broadword select to walk the set bits, with no `PEXT` and no table.**
//!
//! Ticket: R-108. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p108
//! ```
//!
//! Writes `docs/experiments/p-108.csv`.
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
    let prereg = isomesh::experiment!("P-108");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-108 harness not written yet");
    });
}
