**What this changes**

<!-- One concept per PR, matching the repo's one-concept-per-commit habit. -->

**The gates**

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
cargo test -p isomesh
./scripts/backlog_gate.sh
./scripts/readme_sync.sh
cd bevy_isomesh && cargo test        # if the Bevy side is touched
```

- [ ] The gates above pass locally.
- [ ] Any performance claim comes with the committed benchmark that produced it, and names the
      machine (see the Machines block in `FINDINGS.md`).
- [ ] If a measurement contradicted something written down, `FINDINGS.md` gained the entry — the
      contradiction is the finding.
- [ ] If this completes a `BACKLOG.md` ticket, the row moved to `BACKLOG_ARCHIVE.md` with an
      annotation, in this same PR.
