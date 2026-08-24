#!/usr/bin/env bash
#
# Everything to run before committing, in one command — because "run both
# workspaces" written in prose is not a gate, and M-293 is the proof.
#
# `bevy_isomesh` is excluded from the root workspace deliberately (feature
# unification, M-190), so `cargo check --workspace --all-targets` does not
# compile it. F-001 changed a `ReferenceField` method, broke an example there,
# and **58 commits went by with every local gate green**. CI had the right
# command all along and simply had not seen the branch — 75 commits unpushed.
# `CLAUDE.md` already said "run both" in two places, which is the evidence that
# saying it does not work.
#
#   scripts/preflight.sh          the fast set: seconds, run it before every commit
#   scripts/preflight.sh --full   adds the test suites and MSRV: minutes, before pushing
#
# The split is the point. `cargo test -p isomesh` alone is over four minutes, and
# a gate nobody waits for is a gate nobody runs — which is how this repository
# got into the state that produced M-293 in the first place. The fast set is
# everything that is cheap and would have caught it.
#
# # This duplicates ci.yml, and that is a known cost
#
# `backlog_gate.sh`'s header states the rule this breaks: CI runs that file
# rather than a copy of its logic, so there is one path. Here there are two —
# this script and `.github/workflows/ci.yml` — and they can drift. The right fix
# is for CI to call this, which is a change to the workflow and therefore the
# repository owner's; until then the duplication is stated here rather than
# discovered later.

set -uo pipefail

cd "$(dirname "$0")/.."

FULL=0
case "${1:-}" in
    --full) FULL=1 ;;
    --fast | "") FULL=0 ;;
    *)
        printf 'usage: %s [--fast|--full]\n' "$0" >&2
        exit 2
        ;;
esac

failed=()
step() {
    local name=$1
    shift
    printf '\n\033[1m── %s\033[0m\n' "$name"
    if "$@"; then
        printf '   ok\n'
    else
        printf '   \033[31mFAILED\033[0m\n'
        failed+=("$name")
    fi
}

# `cargo` in a subdirectory needs the directory, and `step` takes a command, so
# the two-workspace steps go through this.
in_bevy() { (cd bevy_isomesh && "$@"); }

# ── the fast set ──────────────────────────────────────────────────────────────
# Ordered cheapest-first, so a formatting slip fails in a second rather than
# after a clippy run.
step "root: fmt" cargo fmt --all --check
step "bevy: fmt" in_bevy cargo fmt --all --check
step "backlog gate" ./scripts/backlog_gate.sh
step "findings index" ./scripts/findings_index.sh --check
step "doc facts" ./scripts/doc_facts.sh
step "readme sync" ./scripts/readme_sync.sh
step "toolchain drift" ./scripts/toolchain_drift.sh
step "root: clippy" cargo clippy --workspace --all-targets -- -D warnings
# **The step whose absence is M-293.** A type-check of every target in the other
# workspace, which is exactly what CI's `bevy` job runs and what nothing local
# ran.
step "bevy: check --all-targets" in_bevy cargo check --all-targets
step "bevy: clippy" in_bevy cargo clippy --all-targets -- -D warnings
step "root: rustdoc" env RUSTDOCFLAGS=-D\ warnings cargo doc --workspace --no-deps
step "bevy: rustdoc" in_bevy env RUSTDOCFLAGS=-D\ warnings cargo doc --no-deps
# Rule 2 and rule 3, which are cheap and are the crate's whole pitch.
step "no bevy in the resolved graph" bash -c \
    '[ "$(cargo metadata --format-version 1 | grep -c "\"name\":\"bevy")" = 0 ]'
step "isomesh depends on libm and nothing else" bash -c \
    '[ "$(cargo tree -p isomesh -e normal | wc -l)" = 2 ]'

# ── the full set ──────────────────────────────────────────────────────────────
if [ "$FULL" = 1 ]; then
    step "root: test" cargo test -p isomesh
    step "gpu: test" cargo test -p isomesh-gpu
    # `--lib` and `--doc` rather than a bare `cargo test`, matching the `bevy` CI
    # job for the reason that job's comment gives: `cargo test` also links all 49
    # examples at ~2.0 GB each, which is ~98 GB of writes here and does not fit on
    # a CI runner at all. Coverage is identical -- this crate has no `tests/`
    # directory. Two invocations because cargo refuses `--lib --doc` together.
    step "bevy: test" in_bevy cargo test --lib
    step "bevy: doc tests" in_bevy cargo test --doc
    step "msrv 1.89" cargo +1.89 check --workspace --all-targets
fi

printf '\n'
if [ ${#failed[@]} -eq 0 ]; then
    if [ "$FULL" = 1 ]; then
        printf '\033[32mpreflight --full: all green\033[0m\n'
    else
        printf '\033[32mpreflight: all green\033[0m — run --full before pushing\n'
    fi
    exit 0
fi
printf '\033[31mpreflight FAILED\033[0m: %s\n' "${failed[*]}"
exit 1
