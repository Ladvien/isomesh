#!/usr/bin/env bash
#
# Identify the machine a measurement was taken on.
#
#   scripts/machine.sh --slug    one token, for a filename
#   scripts/machine.sh --spec    the full provenance, as CSV comment lines
#
# # Why both, and why the specs are not in the filename (T-014)
#
# `resolution_sweep-ryzen9-5900x.csv` was written by hand and shows the problem:
# the filename has to stay short enough to type, so it carries a nickname and
# loses everything that makes the numbers interpretable -- core count, memory,
# kernel, and above all which *commit* produced them. A timing figure whose
# compiler version is unknown is a rumour.
#
# So the filename gets a slug and the file gets a header. `--spec` emits `#`
# lines, which every reader in this repo skips and which spreadsheets treat as a
# leading text row rather than choking on.
#
# The commit is the field that matters most and the one most easily lost: two
# runs a week apart on the same box are not comparable if the extractor changed
# between them, and nothing else in the CSV records that it did.

set -euo pipefail
cd "$(dirname "$0")/.."

cpu_model() {
    if [ -r /proc/cpuinfo ]; then
        grep -m1 '^model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//' && return 0
    fi
    if command -v sysctl >/dev/null 2>&1; then
        sysctl -n machdep.cpu.brand_string 2>/dev/null && return 0
    fi
    echo "unknown"
}

slug() {
    echo "$(cpu_model)" \
        | tr '[:upper:]' '[:lower:]' \
        | sed -E 's/\(r\)|\(tm\)|cpu|processor|@.*//g; s/[^a-z0-9]+/-/g; s/^-+|-+$//g'
}

spec() {
    local cores mem kernel rustc commit dirty
    cores="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo '?')"
    if [ -r /proc/meminfo ]; then
        mem="$(awk '/^MemTotal:/ {printf "%.0f GB", $2/1048576}' /proc/meminfo)"
    else
        mem="?"
    fi
    kernel="$(uname -sr)"
    rustc="$(rustc --version 2>/dev/null || echo '?')"
    commit="$(git rev-parse --short HEAD 2>/dev/null || echo '?')"
    # Only code can invalidate a measurement. A dirty README or a script edit
    # cannot change a triangle count or a timing, and flagging those would train
    # everyone to ignore the warning that matters.
    dirty=""
    git diff --quiet -- crates bevy_isomesh Cargo.toml Cargo.lock 2>/dev/null \
        || dirty=" (code modified since this commit -- these numbers are not attributable)"

    echo "# machine: $(cpu_model), ${cores} logical cores, ${mem}"
    echo "# system:  ${kernel}"
    echo "# rustc:   ${rustc}"
    echo "# commit:  ${commit}${dirty}"
    echo "# taken:   $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
}

case "${1:---slug}" in
--slug) slug ;;
--spec) spec ;;
*)
    echo "usage: $0 [--slug|--spec]" >&2
    exit 2
    ;;
esac
