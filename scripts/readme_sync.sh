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

# The **first** ```rust fence, which is the quickstart this exists to pin.
#
# It used to be every such fence, concatenated, so adding a second Rust snippet
# anywhere in either file made this fail with a diff that pointed at the new
# snippet rather than at any drift in the quickstart. That is a gate failing for
# a reason it is not about, which is how gates get switched off.
extract() {
  awk '/^```rust$/ && !seen {inside=1; seen=1; next} /^```/{inside=0} inside' "$1"
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
