#!/usr/bin/env bash
#
# An experiment CSV's `# commit <sha>` header exists so the dataset can be
# traced to the code that produced it. Nothing checked it, and it stopped being
# true without anything looking wrong.
#
# The 2026-08-26 audit found it: `docs/experiments/p-57.csv`, `p-58.csv` and
# `p-60.csv` all carry `# commit 0d81ae6`, and `0d81ae6` is **not in this
# branch's history** — a rebase rewrote it away after the header was written.
# An independent count while writing this gate found two more the audit missed,
# `p-41.csv` and `p-42.csv` at `cceb112`, so it is five datasets and not three.
# The stamp did not fail; the mechanism had a hole, because the rewrite happens
# *after* the number is written and nothing looks again.
#
# # Why ancestry rather than `git cat-file -e`
#
# `git cat-file -e 0d81ae6^{commit}` **succeeds** in the working clone: the
# object is still on disk, reachable from a reflog, unreferenced by any branch
# and one `gc` from gone. So an existence check is green here and red in a fresh
# clone of the same repository, which makes it a property of the checkout rather
# than of the artefact. `git merge-base --is-ancestor <sha> HEAD` asks the
# question the header is for — *can a reader check this out* — and answers it
# the same way everywhere.
#
# # The three inherited debts are pinned by name, and the pins must shrink
#
# 34 of the 50 CSVs in the tree were run against a dirty tree, five name a SHA
# that left the history, and `p-26.csv` predates the header entirely. Rewriting
# those headers would be a lie — the runs really were dirty — so each set is
# listed below and asserted **exactly**: a new violation fails, and a file that
# gets re-run clean also fails until its name is removed. That is deliberate.
# An allow-list checked with `⊆` rots into a list of files that no longer
# violate anything, which is `M-4`'s rule (pin a known defect as an assertion,
# not as an exclusion) in a second place.
#
# Exit 0 if every governed CSV resolves; exit 1 naming the file and the clause.

set -uo pipefail

cd "$(dirname "$0")/.."

fail=0

problem() {
    printf '\033[31mcsv provenance:\033[0m %s\n' "$1" >&2
    shift
    for line in "$@"; do printf '  %s\n' "$line" >&2; done
    fail=1
}

# --- the inherited debts, one list per clause ---------------------------------
#
# Recorded 2026-08-26 at D-017, from the plan for Phase 23 and the audit at
# `docs/research/2026-08-26-audit-and-phase-23-registrations.md` (A-1, A-2).

# SHAs a history rewrite removed. Nothing in the repository identifies the code
# that produced these five datasets; the harnesses are committed and re-runnable
# and that is a weaker guarantee than the header claims.
ANCESTRY_DEBT=(p-41.csv p-42.csv p-57.csv p-58.csv p-60.csv)

# Runs made against a modified tree. Honest and visible on the artefact, which
# is what the flag is for, and not traceable.
DIRTY_DEBT=(
    p-8.csv p-9.csv p-12.csv p-14.csv p-15.csv p-16.csv p-19.csv p-20.csv
    p-21.csv p-22.csv p-23.csv p-38.csv p-39.csv p-40.csv p-41.csv p-42.csv
    p-43.csv p-44.csv p-45.csv p-46.csv p-47.csv p-48.csv p-49.csv p-50.csv
    p-51.csv p-52.csv p-53.csv p-54.csv p-55.csv p-56.csv p-57.csv p-58.csv
    p-59.csv p-60.csv
)

# Written before `common::experiment` stamped anything.
NO_HEADER_DEBT=(p-26.csv)

contains() {
    local needle=$1
    shift
    for item in "$@"; do [ "$item" = "$needle" ] && return 0; done
    return 1
}

# --- the scan -----------------------------------------------------------------

observed_ancestry=()
observed_dirty=()
observed_no_header=()
governed=0

for path in docs/experiments/*.csv; do
    name=$(basename "$path")
    header=$(grep -m1 '^# commit ' "$path" || true)

    if [ -z "$header" ]; then
        observed_no_header+=("$name")
        contains "$name" "${NO_HEADER_DEBT[@]}" ||
            problem "$name has no '# commit' provenance header" \
                "every experiment CSV is written by common::experiment::run, which stamps one"
        continue
    fi

    governed=$((governed + 1))
    sha=$(printf '%s' "$header" | awk '{ print $3 }')

    case "$header" in
        *'WORKING TREE DIRTY'*)
            observed_dirty+=("$name")
            contains "$name" "${DIRTY_DEBT[@]}" ||
                problem "$name was run against a dirty working tree" \
                    "$header" \
                    "re-run the bench with a clean tree so the numbers name a commit"
            ;;
    esac

    # `unknown` is what `ask()` writes when the machine has no `git`, and its
    # own documentation says so: honest, visible, and not a traceability claim
    # that can be false.
    [ "$sha" = "unknown" ] && continue

    if ! git merge-base --is-ancestor "$sha" HEAD 2>/dev/null; then
        observed_ancestry+=("$name")
        detail="not in HEAD's history"
        git cat-file -e "$sha^{commit}" 2>/dev/null &&
            detail="on disk but unreferenced — one 'git gc' from unreadable"
        contains "$name" "${ANCESTRY_DEBT[@]}" ||
            problem "$name names commit $sha, which is $detail" \
                "$header" \
                "the header exists so a reader can check this out; re-run the bench"
    fi
done

# --- the pins must shrink, not rot -------------------------------------------

check_pin() {
    local clause=$1 pinned=$2 observed=$3
    local stale
    stale=$(comm -23 <(printf '%s\n' $pinned | sort -u) <(printf '%s\n' $observed | sort -u))
    if [ -n "$stale" ]; then
        problem "$clause: these names are pinned as violations and no longer violate it" \
            $stale "remove them from the list in $0"
    fi
}

check_pin "ancestry" "${ANCESTRY_DEBT[*]}" "${observed_ancestry[*]:-}"
check_pin "dirty tree" "${DIRTY_DEBT[*]}" "${observed_dirty[*]:-}"
check_pin "missing header" "${NO_HEADER_DEBT[*]}" "${observed_no_header[*]:-}"

# --- verdict ------------------------------------------------------------------
if [ "$fail" -ne 0 ]; then
    printf '\033[31mcsv provenance: FAILED\033[0m\n' >&2
    exit 1
fi

printf 'csv provenance: %s headers, %s resolve to a commit in history, %s inherited debts pinned\n' \
    "$governed" \
    "$((governed - ${#observed_ancestry[@]}))" \
    "$((${#ANCESTRY_DEBT[@]} + ${#DIRTY_DEBT[@]} + ${#NO_HEADER_DEBT[@]}))"
