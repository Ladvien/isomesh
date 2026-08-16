#!/usr/bin/env bash
#
# CI installs whatever `stable` is on the day it runs. A dev machine installs
# `stable` once and keeps it until somebody remembers to update. Clippy gains
# lints between releases, so those two facts together mean **a green
# `preflight.sh` does not imply a green CI** -- and on 2026-08-16 that cost the
# 0.0.6 release: local clippy 0.1.96 was clean, CI's 0.1.97 raised
# `question_mark` twice in `predicates.rs`, the lint job failed, and the publish
# step was skipped. See M-304.
#
# This is that check. It is deliberately a *warning* and not a failure: being a
# patch release behind is normal and blocking every local run on it would be the
# gate-that-cries-wolf Part 5 already warns about. What it must not do is stay
# silent, because silence is what made the drift invisible.
#
# Exit 0 always. The signal is the message.

set -euo pipefail

cd "$(dirname "$0")/.."

# `rustup check` is the only thing that knows what the current stable is without
# guessing at a release calendar. If rustup is absent -- a CI container, a
# distro-packaged toolchain -- there is nothing to compare and nothing to say.
if ! command -v rustup >/dev/null 2>&1; then
    echo "   rustup absent; cannot compare against the current stable"
    exit 0
fi

check=$(rustup check 2>/dev/null || true)
if [ -z "$check" ]; then
    echo "   rustup check produced nothing; skipping"
    exit 0
fi

stable_line=$(printf '%s\n' "$check" | grep -E '^stable-' || true)
if [ -z "$stable_line" ]; then
    echo "   no stable toolchain line; skipping"
    exit 0
fi

if printf '%s' "$stable_line" | grep -q 'update available'; then
    have=$(printf '%s' "$stable_line" | sed -E 's/.*: ([0-9.]+).*->.*/\1/')
    want=$(printf '%s' "$stable_line" | sed -E 's/.*-> ([0-9.]+).*/\1/')
    printf '\n::warning::local stable is %s, current stable is %s\n' "$have" "$want" >&2
    cat >&2 <<EOF
    CI installs the current stable, so its clippy knows lints yours does not.
    A clean run here is NOT evidence that CI's lint job will pass.

        rustup update stable

    M-304: this exact drift failed the 0.0.6 release with the lint job red and
    the publish step skipped, after a fully green local preflight.
EOF
    exit 0
fi

echo "   local stable matches the current release"
