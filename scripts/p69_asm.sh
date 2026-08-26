#!/usr/bin/env bash
#
# P-69's registered verification requirement: assembly for the **monomorphised
# `f32` instance** of the sample loop.
#
# `cargo-show-asm` is not installed here and needs no install for this:
# `cargo rustc -- --emit=asm` writes one `.s` per codegen unit, and with
# `codegen-units = 1` the bench's own unit holds every monomorphisation of
# `push_loop` and `row_loop` the harness instantiates. `experiment_p69.rs`
# carries four `#[inline(never)]` monomorphic probes for exactly this reason —
# with thin LTO both loops otherwise inline into their callers and leave no
# symbol to inspect, so a dump without them shows nothing and looks like
# evidence.
#
# # What the classification is, and why it is a count rather than a look
#
# A Criterion delta cannot distinguish a vectorised loop from a lucky one, which
# is the registration's stated reason for requiring this. The discriminator is
# **`%ymm`**: AVX2's 256-bit registers. A loop that widened four or eight `f32`
# lanes has them; a loop that did not has only `%xmm` scalars. `packed` counts
# `...ps`/`...pd` arithmetic and `scalar` counts `...ss`/`...sd`, and the packed
# count is *not* on its own evidence — `xorps` and `andps` are the ordinary way
# to manipulate a sign bit on one scalar.
#
# Prints one row per monomorphisation. Exit 0 always: this reports, it does not
# gate. The gate is `golden_hashes_are_unchanged`.

set -uo pipefail

cd "$(dirname "$0")/.."

echo "── emitting assembly for the bench's own codegen unit ──"
rm -f target/release/deps/experiment_p69-*.s
# `cargo rustc` is a no-op if the fingerprint is current, and then the `.s`
# files just deleted are not re-emitted and the check reads as "no probes".
# Touching the source is what makes `--emit=asm` actually run.
touch crates/isomesh/benches/experiment_p69.rs
cargo rustc -p isomesh --release --bench experiment_p69 -- --emit=asm -C debuginfo=0 \
    >/dev/null 2>&1 || {
    printf 'p69_asm: the bench did not build\n' >&2
    exit 1
}

UNIT=$(grep -l 'asm_probe' target/release/deps/experiment_p69-*.s 2>/dev/null | head -1)
if [ -z "$UNIT" ]; then
    printf 'p69_asm: no codegen unit carries the asm probes — did they get \n'  >&2
    printf '         eliminated? They are kept alive by black_box in main.\n' >&2
    exit 1
fi
echo "unit: $UNIT"

python3 - "$UNIT" <<'PY'
import collections, re, sys

text = open(sys.argv[1], errors="replace").read()
labels = [(m.start(), m.group(1)) for m in
          re.finditer(r"^(\S*(?:row_loop|push_loop)\S*):$", text, re.M)]
if not labels:
    sys.exit("p69_asm: the unit has no row_loop/push_loop symbol")


def field_of(sym):
    for needle, name in (("6Spheref", "sphere"),
                         ("8BoxExact", "box_exact"),
                         ("12Intersection", "gyroid")):
        if needle in sym:
            return name
    return "?"


print(f"\n{'field':10s} {'shape':5s} {'scalar':6s} {'lines':>6s} {'%ymm':>6s} "
      f"{'packed':>7s} {'scalar_op':>10s} {'calls':>6s}  transcendental")
total_ymm = 0
for pos, sym in labels:
    end = text.find(f"\n\t.size\t{sym}", pos)
    body = text[pos:end] if end > 0 else text[pos:]
    shape = "row" if "row_loop" in sym else "push"
    scalar = "f64" if re.search(r"(?:row|push)_loopd", sym) else "f32"
    ymm = len(re.findall(r"%ymm", body))
    total_ymm += ymm
    packed = len(re.findall(r"\bv?(?:mul|add|sub|sqrt|max|min|and|xor)p[sd]\b", body))
    scalar_op = len(re.findall(r"\bv?(?:mul|add|sub|sqrt|max|min)s[sd]\b", body))
    calls = collections.Counter(re.findall(r"call\w*\s+(\S+)", body))
    trans = sum(v for k, v in calls.items()
                if re.search(r"3sin|3cos|sinf|cosf", k))
    print(f"{field_of(sym):10s} {shape:5s} {scalar:6s} {body.count(chr(10)):>6d} "
          f"{ymm:>6d} {packed:>7d} {scalar_op:>10d} {sum(calls.values()):>6d}  {trans}")

print(f"\ntotal %ymm across every monomorphisation: {total_ymm}")
print("A zero here means LLVM widened nothing, at either loop shape, on any "
      "field — which is P-69's C1 answered from the machine code rather than "
      "from a stopwatch.")
PY
