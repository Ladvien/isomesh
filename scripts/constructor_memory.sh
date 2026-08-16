#!/usr/bin/env bash
#
# Peak working set per constructor, measured out of process.
#
#   scripts/constructor_memory.sh                 sphere at 33^3 and 65^3
#   scripts/constructor_memory.sh torus 17,33,65  any reference field, any sizes
#
# Writes docs/measurements/constructor_memory.csv.
#
# # Why this is a script and not a column in the benchmark
#
# T-018 reports `out_kib` exactly and leaves peak memory unmeasured, because the
# obvious instrument is a counting `GlobalAlloc` and this workspace sets
# `unsafe_code = "forbid"` -- the basis of the crate's "100% safe Rust" claim
# (M-147). Estimating a working set from reading the algorithm would be a
# performance number with no benchmark behind it, which rule 4 forbids.
#
# So the instrument is the operating system. Each constructor runs in its own
# process via `--only`, which does that one thing and stops, and the kernel's own
# high-water mark is read afterwards. The crate stays free of `unsafe` rather
# than the rule being bent for a column.
#
# # The baseline is what makes the number mean something
#
# A process's peak RSS includes the Rust runtime, the binary's pages, and the
# input grid the constructor was handed. `baseline` builds that input and exits;
# `mesh-baseline` builds the input *and* the mesh, because the mesh-based
# constructors need one and it is not their cost. Every reported figure is a
# subtraction, and a subtraction that comes out at or below zero is reported as
# zero rather than as a negative. That is a real answer, not a failure: it means
# the constructor's own peak sits below the peak of building the input it was
# handed -- which is what the mesh-based pair measure, since meshing allocates
# more transiently than they do.

set -euo pipefail
cd "$(dirname "$0")/.."

FIELD="${1:-sphere}"
SAMPLE_LIST="${2:-33,65}"
OUT="docs/measurements/constructor_memory.csv"

# The bench binary. `--no-run` builds it and the JSON tells us where it landed,
# which beats globbing `target/release/deps` for a hash-suffixed name.
BIN="$(cargo bench -p isomesh --bench constructors --no-run --message-format=json 2>/dev/null \
  | python3 -c '
import json, sys
for line in sys.stdin:
    try:
        msg = json.loads(line)
    except ValueError:
        continue
    if msg.get("reason") == "compiler-artifact" and msg.get("executable"):
        if msg.get("target", {}).get("name") == "constructors":
            print(msg["executable"])
')"

if [ -z "${BIN:-}" ] || [ ! -x "$BIN" ]; then
  echo "::error::could not locate the constructors bench binary" >&2
  exit 1
fi

# Peak RSS in KiB for one run.
#
# The bench binary reads its own `VmHWM` from `/proc/self/status` and prints it
# as the last field. That is a plain file read -- ordinary safe Rust. What was
# ruled out is a counting `GlobalAlloc`, which needs `unsafe impl` and which this
# workspace forbids; nothing about reading a kernel counter required leaving the
# crate. What it does require is leaving the *process*, because `VmHWM` is a
# high-water mark over a whole process life and would otherwise attribute every
# constructor's peak to whichever ran first.
#
# Linux only. macOS publishes no `/proc`, and its `ru_maxrss` needs `libc`, which
# `isomesh` does not depend on. CI runs Linux, which is where the figure is
# checked.
peak_kib() {
  local out
  out="$("$@" 2>/dev/null | tail -1)"
  echo "${out##* }"
}

run() {
  peak_kib "$BIN" --bench --only "$1" --field "$FIELD" --samples "$SAMPLES"
}

mkdir -p "$(dirname "$OUT")"
printf 'field,samples,constructor,peak_kib,above_baseline_kib\n' > "$OUT"

for SAMPLES in ${SAMPLE_LIST//,/ }; do
  echo "constructor memory -- $FIELD at ${SAMPLES}^3, peak RSS above baseline"
  printf '%-16s %12s %14s\n' "constructor" "peak KiB" "above base"

  BASE="$(run baseline)"
  if [ "$BASE" = "unavailable" ]; then
    echo "::error::this platform publishes no /proc/self/status; see the header" >&2
    exit 1
  fi
  MESH_BASE="$(run mesh-baseline)"
  printf '%-16s %12s %14s\n' "baseline" "$BASE" "-"
  printf '%-16s %12s %14s\n' "mesh-baseline" "$MESH_BASE" \
    "$((MESH_BASE > BASE ? MESH_BASE - BASE : 0))"
  printf '%s,%s,baseline,%s,0\n' "$FIELD" "$SAMPLES" "$BASE" >> "$OUT"
  printf '%s,%s,mesh-baseline,%s,%s\n' "$FIELD" "$SAMPLES" "$MESH_BASE" \
    "$((MESH_BASE > BASE ? MESH_BASE - BASE : 0))" >> "$OUT"

  for c in exact swept marched band; do
    p="$(run "$c")"
    d=$((p > BASE ? p - BASE : 0))
    printf '%-16s %12s %14s\n' "$c" "$p" "$d"
    printf '%s,%s,%s,%s,%s\n' "$FIELD" "$SAMPLES" "$c" "$p" "$d" >> "$OUT"
  done

  for c in pseudonormal winding; do
    p="$(run "$c")"
    d=$((p > MESH_BASE ? p - MESH_BASE : 0))
    printf '%-16s %12s %14s\n' "$c" "$p" "$d"
    printf '%s,%s,%s,%s,%s\n' "$FIELD" "$SAMPLES" "$c" "$p" "$d" >> "$OUT"
  done
  echo
done

echo "wrote $OUT"
