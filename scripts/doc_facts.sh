#!/usr/bin/env bash
#
# The documentation states counts, and the counts rot.
#
# An audit on 2026-08-16 found the reference-field count wrong in six places
# ("seven reference fields" -- there are eight), the shootout's algorithm count
# wrong in two, and the golden-hash count wrong in the root README for the
# **second** time: CHANGELOG 0.0.4 records fixing that exact number from 147 to
# 168, and it is now 216. Every one was found by a person reading, which is the
# failure mode `backlog_gate.sh` exists for. `readme_sync.sh` was the only
# doc-content check in CI and it compares one code fence, so nothing watched the
# prose.
#
# # Only facts with one mechanical source, and only phrasings that assert a total
#
# The first version of this file checked any number-word before "fields" or
# "extractors" and was unusable: **"five of the seven reference fields go to
# exactly zero"** is correct prose about a subset, and so is "one field, one
# grid, two extractors". A gate that cries wolf is a gate somebody disables, so
# the patterns here are deliberately narrow -- each one is a phrase that can
# only mean the total.
#
# Note the subtlety in that same example: the *subset* "five" is fine and the
# *total* "seven reference fields" inside it is wrong, and this catches the
# second without touching the first, because `reference fields` is the term of
# art for the whole set and `fields` alone is not.
#
# The rule for a writer is: **say nothing, or say the true number.** No document
# is required to mention any of this.
#
# Exit 0 if every document agrees with the source; exit 1 naming the file, the
# line, and both numbers.

set -uo pipefail

cd "$(dirname "$0")/.."

fail=0

problem() {
    printf '\n::error::%s\n' "$1" >&2
    shift
    if [ "$#" -gt 0 ]; then
        printf '    %s\n' "$@" >&2
    fi
    fail=1
}

# Documents whose prose is checked.
#
# `FINDINGS.md`, `BACKLOG.md` and `BACKLOG_ARCHIVE.md` are **excluded on
# purpose**: all three are append-only records of what was true when each entry
# was written, and a finding that says "five of the seven fields" was correct on
# the day it was measured. Rewriting history to satisfy a gate is the opposite of
# what those files are for -- FINDINGS' own rule is "re-tier rather than
# rewrite". `docs/research/` is excluded for the same reason.
DOCS=(
    README.md
    CLAUDE.md
    crates/isomesh/README.md
    crates/isomesh-gpu/README.md
    bevy_isomesh/README.md
    bevy_isomesh/DEMOS.md
    docs/demos/algorithms.md
    docs/demos/correctness.md
    docs/demos/gameplay.md
    docs/experiments.md
    docs/measurements/README.md
)

WORDS='zero|one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve'
WORD=(zero one two three four five six seven eight nine ten eleven twelve)

# `$1` the true count, `$2` what it counts, `$3` the regex that must follow the
# number for this to be a claim about the total.
check() {
    local truth=$1 name=$2 phrase=$3
    local truth_word=${WORD[$truth]:-__none__}
    local doc line text claim

    for doc in "${DOCS[@]}"; do
        [ -f "$doc" ] || continue
        while IFS=: read -r line text; do
            [ -n "${line:-}" ] || continue
            claim=$(printf '%s' "$text" | grep -oiE "($WORDS|[0-9]+)[ -]+$phrase" | head -1)
            problem "$doc:$line claims \"$claim\" — there are $truth $name" \
                "$(printf '%s' "$text" | sed -E 's/^ +//' | cut -c1-108)"
        done < <(grep -inE "($WORDS|[0-9]+)[ -]+$phrase" "$doc" 2>/dev/null |
            grep -viE "($truth|$truth_word)[ -]+$phrase" || true)
    done
}

# --- the sources of truth -----------------------------------------------------

# Reference fields: the macro every bench and every validity test iterates.
FIELDS=$(sed -n '/macro_rules! for_each_reference_field/,/^}/p' \
    crates/isomesh/src/fields/mod.rs | grep -c 'let \$name = "' || true)

# Extractors: the one macro that forwards the inherent `extract` methods.
EXTRACTORS=$(sed -n '/^forward_extractor!(/,/^);/p' crates/isomesh/src/extractor.rs |
    grep -c '^\s*crate::' || true)

# Golden hashes: the committed fixture itself.
HASHES=$(python3 -c 'import json,sys; print(len(json.load(open(sys.argv[1]))))' \
    crates/isomesh/golden_hashes.json 2>/dev/null || echo 0)

# The shootout's shape, from the committed CSV.
SHOOTOUT_ALGS=$(tail -n +2 docs/measurements/shootout.csv 2>/dev/null |
    grep -v '^#' | cut -d, -f2 | sort -u | grep -c . || echo 0)

# Bevy examples on disk. `common/` is a shared module, not an example.
EXAMPLES=$(find bevy_isomesh/examples -maxdepth 1 -name '*.rs' | grep -c . || true)

for pair in "FIELDS:$FIELDS" "EXTRACTORS:$EXTRACTORS" "HASHES:$HASHES" \
    "SHOOTOUT_ALGS:$SHOOTOUT_ALGS" "EXAMPLES:$EXAMPLES"; do
    if [ "${pair#*:}" -lt 1 ]; then
        problem "cannot derive ${pair%%:*} — the source moved and this gate is blind"
    fi
done

# --- the checks ---------------------------------------------------------------
#
# Each phrase below can only be a claim about the whole set. Adding one is
# cheap; adding a loose one costs everybody who runs this.

check "$FIELDS" "reference fields" 'reference fields'
check "$EXTRACTORS" "extractors" 'extractors, side by side'
check "$EXTRACTORS" "extractors" 'isosurface extractors'
check "$HASHES" "golden hashes" 'golden hashes'
check "$SHOOTOUT_ALGS" "algorithms in the shootout" '[-]?algorithm shootout'

# --- media referenced but absent ----------------------------------------------
#
# A broken image is worse than a missing section: it renders as alt text and a
# torn icon, and nothing else in CI opens these files. Badge URLs are skipped --
# they are generated by a service and there is nothing on disk to find.
for doc in "${DOCS[@]}"; do
    [ -f "$doc" ] || continue
    while read -r ref; do
        [ -n "$ref" ] || continue
        case "$ref" in
        https://img.shields.io/* | https://github.com/*/badge.svg) continue ;;
        https://raw.githubusercontent.com/ladvien/isomesh/main/*)
            path=${ref#https://raw.githubusercontent.com/ladvien/isomesh/main/}
            ;;
        http*) continue ;;
        /*) path=$ref ;;
        *) path="$(dirname "$doc")/$ref" ;;
        esac
        [ -e "$path" ] || problem "$doc references media that is not on disk" "$ref → $path"
    done < <(grep -ohE '\]\([^) ]*\.(gif|png|jpg|jpeg|webp|svg)\)' "$doc" 2>/dev/null |
        sed -E 's/^\]\(//; s/\)$//' | sort -u || true)
done

# --- verdict ------------------------------------------------------------------
if [ "$fail" -ne 0 ]; then
    printf '\ndoc facts FAILED\n' >&2
    exit 1
fi

printf 'doc facts: %s fields, %s extractors, %s golden hashes, %s shootout algorithms, %s examples — no drift\n' \
    "$FIELDS" "$EXTRACTORS" "$HASHES" "$SHOOTOUT_ALGS" "$EXAMPLES"
