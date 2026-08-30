//! **P-116 — GRAPHGEN's decision-table pipeline against a case table that is already `const`-derived.**
//!
//! Ticket: R-116. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p116
//! ```
//!
//! Writes `docs/experiments/p-116.csv`.
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
    let prereg = isomesh::experiment!("P-116");
    common::experiment::run(prereg, |_run| {
        unimplemented!("P-116 harness not written yet");
    });
}
