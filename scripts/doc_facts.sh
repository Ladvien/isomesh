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
    web/index.md
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

# Playable web demos: the array `scripts/build_web.sh` iterates is the one place
# that decides what is built for the browser, so it is the source and the three
# prose sites follow it. One name per line, no comments inside the array, which
# is what makes this a `grep -c` rather than a parser.
PLAYABLE=$(sed -n '/^DEMOS=($/,/^)$/p' scripts/build_web.sh | grep -cE '^ +[a-z_]+$' || true)

for pair in "FIELDS:$FIELDS" "EXTRACTORS:$EXTRACTORS" "HASHES:$HASHES" \
    "SHOOTOUT_ALGS:$SHOOTOUT_ALGS" "EXAMPLES:$EXAMPLES" "PLAYABLE:$PLAYABLE"; do
    if [ "${pair#*:}" -lt 1 ]; then
        problem "cannot derive ${pair%%:*} — the source moved and this gate is blind"
    fi
done

# --- the checks ---------------------------------------------------------------
#
# Each phrase below can only be a claim about the whole set. Adding one is
# cheap; adding a loose one costs everybody who runs this.

# `[a-z]* ?` absorbs one adjective, because "seven **shipped** reference
# fields" hid from the first version of this check while "seven reference
# fields" three lines away did not.
check "$FIELDS" "reference fields" '[a-z]* ?reference fields'
check "$EXTRACTORS" "extractors" 'extractors, side by side'
check "$EXTRACTORS" "extractors" 'isosurface extractors'
# `web/index.md`'s lead paragraph, which was ungated until that file joined
# `DOCS` -- and neither phrase above matches it, so adding the file was not
# enough on its own.
check "$EXTRACTORS" "extractors" 'extractors sit behind one trait'
check "$HASHES" "golden hashes" 'golden hashes'
check "$SHOOTOUT_ALGS" "algorithms in the shootout" '[-]?algorithm shootout'

# Derived since the first version of this gate and never checked against the
# prose, which is the failure this whole script exists to prevent -- caught when
# B-014's example made it 35 and both READMEs went on saying 34.
#
# The phrase is deliberately the long one. A bare `examples` also matches
# "Three examples" in `isomesh-gpu/README.md`, which is a true claim about that
# crate's own three and has nothing to do with the Bevy set -- the loose version
# was written first and failed on exactly that, which is this script's own header
# warning ("adding a loose one costs everybody who runs this") landing on its
# author.
check "$EXAMPLES" "examples" 'examples, each with an animated capture'

# The site's own claim, in three files: `web/index.md`, `bevy_isomesh/DEMOS.md`
# and -- once it says so -- `README.md`. It went from three to nine in one
# afternoon, which is exactly the kind of number this script exists for.
check "$PLAYABLE" "playable demos" 'of the demos are playable'

# A structural check the prose helper cannot express, because it compares two
# files rather than a file against a count. Both failure modes are invisible in a
# green CI run: a demo built and not allow-listed is a 36 MB module nothing can
# reach, and one allow-listed and not built is a link to a 404. The `site` job
# passes either way.
ALLOWLIST=$(sed -n '/const DEMOS = {/,/};/p' web/play.html | grep -cE '^ +[a-z_]+:' || true)
if [ "$ALLOWLIST" -ne "$PLAYABLE" ]; then
    problem "web/play.html allow-lists $ALLOWLIST demos, scripts/build_web.sh builds $PLAYABLE"
fi

# And one per demo in the other direction: every allow-listed demo needs the
# `#notes-<name>` block `play.html`'s script unhides, or the page throws on a
# null `hidden` assignment before the module ever loads.
while read -r demo; do
    [ -n "$demo" ] || continue
    grep -qF "id=\"notes-$demo\"" web/play.html ||
        problem "web/play.html allow-lists $demo with no #notes-$demo block"
done < <(sed -n '/^DEMOS=($/,/^)$/p' scripts/build_web.sh | grep -oE '^ +[a-z_]+$' |
    tr -d ' ' || true)

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
        # `web/index.md` is rendered as the site *root*, so its relative targets
        # resolve against the repository root rather than against `web/` -- the
        # same `link_base` column `scripts/build_site.py` carries, and the one
        # thing about that file an implementer is most likely to assume away.
        # Every other document here resolves against its own directory.
        *)
            case "$doc" in
            web/*) path=$ref ;;
            *) path="$(dirname "$doc")/$ref" ;;
            esac
            ;;
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

printf 'doc facts: %s fields, %s extractors, %s golden hashes, %s shootout algorithms, %s examples, %s playable demos — no drift\n' \
    "$FIELDS" "$EXTRACTORS" "$HASHES" "$SHOOTOUT_ALGS" "$EXAMPLES" "$PLAYABLE"
