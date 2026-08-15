#!/usr/bin/env bash
# D-008: the root README's quickstart is a verbatim copy of the crate README's
# doctested one. The crate copy is compiled by every `cargo test` — the
# cfg(doctest) include_str! in crates/isomesh/src/lib.rs — but the root copy
# cannot be: the root README lives outside the package directory, and a
# published crate must never reference a file it does not ship. So the root
# copy is held identical by diff instead. CI's lint job runs this beside the
# backlog gate.
set -euo pipefail
cd "$(dirname "$0")/.."

extract() {
  awk '/^```rust$/{inside=1; next} /^```/{inside=0} inside' "$1"
}

crate="$(extract crates/isomesh/README.md)"
root="$(extract README.md)"

if [[ -z "$crate" ]]; then
  echo 'readme sync: no ```rust fence found in crates/isomesh/README.md' >&2
  exit 1
fi
if [[ -z "$root" ]]; then
  echo 'readme sync: no ```rust fence found in README.md' >&2
  exit 1
fi
if [[ "$root" != "$crate" ]]; then
  echo 'readme sync: the root README rust fence differs from the doctested one in crates/isomesh/README.md' >&2
  diff <(printf '%s\n' "$crate") <(printf '%s\n' "$root") >&2 || true
  exit 1
fi
echo "readme sync: root quickstart matches the doctested crate snippet"
