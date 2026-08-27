#!/usr/bin/env bash
#
# Six manifest sites carry this project's version and nothing checked they
# agree. On `0.0.x` a missed pin is not cosmetic and it is not recoverable.
#
# Cargo reads `version = "0.0.9"` as `>=0.0.9, <0.0.10`, because on a
# zero-zero version every patch bump is a breaking change. So a published
# `bevy_isomesh 0.0.10` carrying an unbumped `isomesh` pin does not warn and
# does not fail to build: it resolves, for every consumer, against the **old**
# `isomesh 0.0.9`. The consumer then gets a crate whose changelog describes
# arithmetic it is not running. And a published version can be yanked but never
# deleted, so the fix is a third release rather than an edit.
#
# The three workspaces are the reason this cannot be one `cargo metadata` call:
# `bevy_isomesh` and `isomesh_web` are excluded from the root workspace (feature
# unification, M-190), so the root's metadata cannot see either of them --
# which is why `publish.sh` special-cases `bevy_isomesh` in its `ORDER` array.
# This reads the manifests as files instead.
#
# CI runs this same file rather than a copy of its logic in the workflow YAML,
# so a developer and CI run the same thing -- `backlog_gate.sh`'s rule.
#
# Exit 0 having printed the agreed version; exit 1 naming every site that
# disagrees, with its file and line.

set -euo pipefail

cd "$(dirname "$0")/.."

fail=0

# Extract one manifest site's version.
#
#   site <file> <section-or-dash> <key-prefix>
#
# Prints `<line>\t<version>`, or nothing if the site is absent.
#
# **The key prefix is what makes this safe.** These manifests also pin `wgpu
# 29.0.3`, Bevy `0.19.x` and two MSRVs, so a bare `grep -o '0\.0\.[0-9]*'`
# would match a dependency's version and report a mismatch that is not one.
# Matching on the *key* -- `version = `, `isomesh = `, `isomesh-gpu = ` -- at
# the start of the line is unambiguous, and `isomesh = ` does not match
# `isomesh-gpu = ` because the space and `=` are part of the prefix.
#
# `section` scopes a bare `version = ` to its stanza, because the root manifest
# carries it under `[workspace.package]` and a package manifest under
# `[package]`; `-` means "anywhere in the file", which is what the two path
# dependency pins need.
site() {
    awk -v want="$2" -v key="$3" '
        /^\[/ { section = $0 }
        index($0, key) == 1 {
            if (want != "-" && section != want) next
            if (match($0, /version *= *"[^"]*"/)) {
                v = substr($0, RSTART, RLENGTH)
                sub(/^version *= *"/, "", v)
                sub(/"$/, "", v)
                printf "%d\t%s\n", NR, v
                exit
            }
        }
    ' "$1"
}

# file:section:key:what-it-is. The order is the order a release edits them.
SITES=(
    'Cargo.toml:[workspace.package]:version = :the workspace version, which drives isomesh and isomesh-gpu'
    'crates/isomesh-gpu/Cargo.toml:-:isomesh = :isomesh-gpu -> isomesh pin'
    'bevy_isomesh/Cargo.toml:[package]:version = :bevy_isomesh own version'
    'bevy_isomesh/Cargo.toml:-:isomesh = :bevy_isomesh -> isomesh pin'
    'bevy_isomesh/Cargo.toml:-:isomesh-gpu = :bevy_isomesh -> isomesh-gpu pin'
    'isomesh_web/Cargo.toml:[package]:version = :isomesh_web own version (publish = false, kept in step so the site names one version)'
)

# The root workspace version is the reference every other site is held against.
ref_site=${SITES[0]}
ref_found=$(site "${ref_site%%:*}" "$(cut -d: -f2 <<<"$ref_site")" 'version = ')
if [ -z "$ref_found" ]; then
    printf '::error::version gate: no `version = ` under `[workspace.package]` in Cargo.toml\n' >&2
    exit 1
fi
WANT=${ref_found#*$'\t'}

for entry in "${SITES[@]}"; do
    file=${entry%%:*}
    rest=${entry#*:}
    section=${rest%%:*}
    rest=${rest#*:}
    key=${rest%%:*}
    what=${rest#*:}

    found=$(site "$file" "$section" "$key")
    if [ -z "$found" ]; then
        printf '::error::version gate: %s has no `%s` site (%s)\n' "$file" "$key" "$what" >&2
        fail=1
        continue
    fi

    line=${found%%$'\t'*}
    got=${found#*$'\t'}
    if [ "$got" != "$WANT" ]; then
        printf '::error::version gate: %s:%s is %s, expected %s -- %s\n' \
            "$file" "$line" "$got" "$WANT" "$what" >&2
        fail=1
    fi
done

if [ "$fail" -ne 0 ]; then
    printf '\nversion gate: FAILED -- a partial bump publishes a crate pinned to the previous release\n' >&2
    exit 1
fi

printf 'version gate: all six sites at %s\n' "$WANT"
