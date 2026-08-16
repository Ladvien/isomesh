#!/usr/bin/env bash
#
# BACKLOG.md is the state, and twice it has drifted from the work it describes.
# E-113 shipped at a0859e8 and its row sat in BACKLOG.md for four commits; the
# header counts disagreed with the row counts on a separate occasion. Both came
# from editing the file several times in one turn, where a later write clobbered
# an earlier one, and both were found by an audit rather than by a check.
#
# This is that check. CI runs this same file rather than a copy of the logic in
# the workflow YAML, so there is one path and it is the one you can run locally.
#
# Exit 0 if the two files agree; exit 1 having printed exactly what disagrees.

set -euo pipefail

cd "$(dirname "$0")/.."

BACKLOG=BACKLOG.md
ARCHIVE=BACKLOG_ARCHIVE.md
FINDINGS="FINDINGS.md"
# Ticket IDs: A-014b, GPU-001, T-005a, E-201.
ID='[A-Za-z]+-[0-9]+[a-z]?'

fail=0

# First argument is the headline; the rest are the offending lines.
problem() {
    printf '\n::error::%s\n' "$1" >&2
    shift
    if [ "$#" -gt 0 ]; then
        printf '    %s\n' "$@" >&2
    fi
    fail=1
}

# A *ticket row* is what counts here, never a bare mention. Both files reference
# archived IDs in prose constantly -- "Blocked by | A-014a", "A-010 has since
# landed and closed that owner" -- and matching bare IDs would call every one of
# those a duplicate. The row shape is the discriminator: an open row opens with
# the empty ballot box, an archived row with the checked one.
rows() {
    grep -oE "^\| *$1 *\| *\*\*$ID\*\*" "$2" | grep -oE "$ID" | sort || true
}

# `comm` wants files; these keep the empty case from turning into a blank line
# that matches another blank line and reports a phantom duplicate.
lines() { printf '%s' "$1" | grep -v '^$' || true; }
count() { lines "$1" | grep -c . || true; }

open=$(rows '☐' "$BACKLOG")
archived=$(rows '☑' "$ARCHIVE")

# --- 1. A ticket lives in exactly one file ------------------------------------
both=$(comm -12 <(lines "$open") <(lines "$archived"))
if [ -n "$both" ]; then
    problem "ticket rows appear in BOTH $BACKLOG and $ARCHIVE" $both
fi

# --- 2. No row in the wrong file ----------------------------------------------
# E-113's exact failure: the work was done, the row was never moved. The mirror
# case -- an unchecked row already in the archive -- is the same mistake made in
# the other direction, so both are checked.
stranded=$(grep -nE "^\| *☑ *\| *\*\*$ID\*\*" "$BACKLOG" || true)
if [ -n "$stranded" ]; then
    problem "checked rows still in $BACKLOG -- they belong in $ARCHIVE" "$stranded"
fi
premature=$(grep -nE "^\| *☐ *\| *\*\*$ID\*\*" "$ARCHIVE" || true)
if [ -n "$premature" ]; then
    problem "unchecked rows in $ARCHIVE -- an unfinished ticket was archived" "$premature"
fi

# --- 3. No ID twice in one file -----------------------------------------------
for pair in "open:$BACKLOG" "archived:$ARCHIVE"; do
    case $pair in
        open:*) set -- "$open" "${pair#open:}" ;;
        *) set -- "$archived" "${pair#archived:}" ;;
    esac
    dupes=$(lines "$1" | uniq -d)
    if [ -n "$dupes" ]; then
        problem "duplicate ticket rows in $2" $dupes
    fi
done

# --- 4. The declared counts match the rows ------------------------------------
# "**52 tickets archived, 29 open.**" near the top of BACKLOG.md, and the
# archive's own "52 tickets." above its index. Three numbers, two sources, and
# nothing but this keeps them honest.
n_open=$(count "$open")
n_archived=$(count "$archived")

header=$(grep -oE '\*\*[0-9]+ tickets archived, [0-9]+ open\.\*\*' "$BACKLOG" || true)
if [ -z "$header" ]; then
    problem "no '**N tickets archived, M open.**' line found in $BACKLOG"
else
    said_archived=$(printf '%s' "$header" | grep -oE '[0-9]+' | head -1)
    said_open=$(printf '%s' "$header" | grep -oE '[0-9]+' | tail -1)
    if [ "$said_archived" != "$n_archived" ] || [ "$said_open" != "$n_open" ]; then
        problem "$BACKLOG header disagrees with its rows" \
            "header says: $said_archived archived, $said_open open" \
            "rows say:    $n_archived archived, $n_open open"
    fi
fi

index_count=$(grep -oE '^[0-9]+ tickets\.' "$ARCHIVE" | grep -oE '[0-9]+' | head -1 || true)
if [ -z "$index_count" ]; then
    problem "no 'N tickets.' line found above the index in $ARCHIVE"
elif [ "$index_count" != "$n_archived" ]; then
    problem "$ARCHIVE index count disagrees with its rows" \
        "index says: $index_count" \
        "rows say:   $n_archived"
fi

# --- 5. Every dependency resolves ---------------------------------------------
# The last cell of a ticket row is its "Blocked by" list, in every table shape
# the file uses (Phase 1 and 2 have six columns, Phase 4a five, Phase 4b six).
# A typo'd or deleted blocker makes a ticket permanently unreachable without
# anything looking wrong, which is the quiet version of the same failure.
known=$(printf '%s\n%s\n' "$open" "$archived" | sort -u)
refs=$(grep -E "^\| *[☐☑] *\| *\*\*$ID\*\*" "$BACKLOG" \
    | awk -F'|' '{ print $(NF - 1) }' \
    | grep -oE "$ID" | sort -u || true)
dangling=$(comm -13 <(lines "$known") <(lines "$refs"))
if [ -n "$dangling" ]; then
    problem "'Blocked by' names tickets that exist in neither file" $dangling
fi

# --- 6. Every pre-registration reaches FINDINGS.md ----------------------------
# R-000 moved the `P-` predictions into `crates/isomesh/src/experiment.rs`, where
# the compiler can refuse an unregistered id. That file is the source; the prose
# in FINDINGS.md elaborates on it. Two copies of a hypothesis drift, and the
# drift is undetectable -- a pre-registration that quietly acquired a clause is
# worse than none -- so this only checks that the prose *mentions* each id, which
# is enough to make an unrecorded experiment visible without pinning the wording.
EXPERIMENTS="crates/isomesh/src/experiment.rs"
if [ -f "$EXPERIMENTS" ]; then
    registered=$(grep -oE 'id: "P-[0-9]+"' "$EXPERIMENTS" | grep -oE 'P-[0-9]+' | sort -u || true)
    for p in $registered; do
        if ! grep -q "$p" "$FINDINGS" 2>/dev/null; then
            problem "$p is pre-registered in $EXPERIMENTS but appears nowhere in $FINDINGS"
        fi
    done
fi

# --- verdict ------------------------------------------------------------------
if [ "$fail" -ne 0 ]; then
    printf '\nbacklog gate FAILED\n' >&2
    exit 1
fi

printf 'backlog gate: %s open, %s archived, no drift\n' "$n_open" "$n_archived"
