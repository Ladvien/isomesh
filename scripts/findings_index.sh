#!/usr/bin/env bash
#
# FINDINGS.md is the epistemic state and it is 387 KB, 945 lines and 298 entries.
# Past a certain size a ledger stops being a lookup table and becomes a diary --
# nobody reads it end to end, so a fact that is in there gets re-derived anyway,
# which is the exact failure the file exists to prevent. It has already happened:
# V-29 and V-32 are the same correction made twice, three days apart, two rows
# apart in the same table.
#
# So the file gets an index, and the index gets a generator, because a
# hand-maintained index of 298 rows is a rotting artefact with extra steps.
# `BACKLOG_ARCHIVE.md`'s index was built by hand and carries rows reading
# "(title not auto-extracted -- grep the ID)", which is what that costs.
#
#   findings_index.sh            rewrite the index block in place
#   findings_index.sh --check    exit 1 if the block is stale, print nothing else
#
# The --check mode is what CI runs, next to backlog_gate.sh. An index that is
# allowed to drift is worse than no index, because it answers confidently and
# wrongly.

set -euo pipefail
cd "$(dirname "$0")/.."

MODE="${1:-write}"

python3 - "$MODE" <<'PY'
import re, sys, os

MODE = sys.argv[1]
PATH = "FINDINGS.md"
BEGIN = "<!-- BEGIN GENERATED INDEX -- scripts/findings_index.sh -->"
END = "<!-- END GENERATED INDEX -->"

src = open(PATH).read()
lines = src.split("\n")


def tidy(text):
    """One line of prose from a table cell or heading."""
    text = re.sub(r"`([^`]*)`", r"\1", text)
    text = text.replace("**", "").replace("~~", "").replace("*", "").strip()
    text = re.sub(r"\s+", " ", text)
    return text


def claim(cell):
    """The entry's claim: its bolded lead where that is a sentence, else the cell.

    Some rows open with a bolded *label* rather than a claim -- O-13 and O-14
    both start `**Pre-registered:**` -- so a bare "first bold run" rule returns
    the label and indexes two entries as identical. Anything under 30 characters
    is treated as a label and the whole cell is used instead.
    """
    lead = re.search(r"\*\*(.+?)\*\*", cell)
    if lead and len(tidy(lead.group(1))) >= 30:
        text = tidy(lead.group(1))
    else:
        text = tidy(cell)
    return text[:120].rstrip() + "…" if len(text) > 120 else text


rows = []
for number, line in enumerate(lines, 1):
    # Table entries: M-, V- and O- all lead their row with the id.
    entry = re.match(r"^\| (M-\d+|V-\d+|O-\d+) \| (.*)$", line)
    if entry:
        body = entry.group(2).split(" | ")[0]
        rows.append((entry.group(1), number, claim(body)))
        continue
    # Falsified entries are headings rather than rows.
    killed = re.match(r"^### (✗\d+) — (.*)$", line)
    if killed:
        rows.append((killed.group(1), number, claim(killed.group(2))))

if not rows:
    sys.exit("findings_index: no entries matched -- the file's shape changed")


def sort_key(row):
    kind = row[0][0]
    order = {"✗": 0, "M": 1, "V": 2, "O": 3}.get(kind, 9)
    return (order, int(re.sub(r"\D", "", row[0])))


rows.sort(key=sort_key)
counts = {}
for identifier, _, _ in rows:
    counts[identifier[0]] = counts.get(identifier[0], 0) + 1

block = [
    BEGIN,
    "",
    f"**{len(rows)} entries** — "
    f"{counts.get('✗', 0)} falsified, {counts.get('M', 0)} measured, "
    f"{counts.get('V', 0)} verified, {counts.get('O', 0)} open. "
    "Regenerate with `scripts/findings_index.sh`; CI fails if this is stale.",
    "",
    "| # | Claim |",
    "|---|---|",
]
for identifier, number, text in rows:
    text = text.replace("|", "\\|")
    block.append(f"| `{identifier}` | {text} |")
block += ["", END]
block = "\n".join(block)

if BEGIN in src and END in src:
    head = src.split(BEGIN)[0]
    tail = src.split(END, 1)[1]
    updated = head + block + tail
else:
    # First run: place the index directly after the how-to-use section, which is
    # where a reader who has just been told to look things up arrives.
    anchor = "## Confidence tiers"
    if anchor not in src:
        sys.exit("findings_index: could not find '## Confidence tiers' to anchor the index")
    updated = src.replace(anchor, "## Index\n\n" + block + "\n\n" + anchor, 1)

if MODE == "--check":
    if updated != src:
        sys.exit(
            "findings_index: FINDINGS.md's index is stale. "
            "Run scripts/findings_index.sh and commit the result."
        )
    print(f"findings index: {len(rows)} entries, current")
else:
    open(PATH, "w").write(updated)
    print(f"findings index: {len(rows)} entries written")
PY
