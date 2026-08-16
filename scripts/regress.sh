#!/usr/bin/env bash
#
# Diff a measurement CSV against its committed baseline and fail on a regression.
#
#   scripts/regress.sh                      check every baseline that exists
#   scripts/regress.sh shootout             check one
#   scripts/regress.sh --accept shootout    adopt the current run as the baseline
#   scripts/regress.sh --self-test          prove the checker detects a regression
#
# # What this catches that nothing else does
#
# T-007's golden hashes catch a changed *mesh*, bit for bit. They are blind to a
# change that leaves every vertex identical and doubles the runtime, and blind to
# one that worsens Hausdorff on a single field while the hashes are not even
# taken. `docs/measurements/*.csv` is rewritten fresh on every run, so before
# this there was no artefact that answered "did that make anything worse".
#
# # Two classes of metric, and only one of them is noisy
#
# The structural columns -- vertices, triangles, non_manifold_edges -- are
# deterministic: the same code on the same field at the same resolution emits the
# same counts on any machine, which T-004 already guarantees. They are compared
# EXACTLY, because a tolerance there would only hide a real change.
#
# The measured columns are not. `median_ms` moves with the machine, the governor
# and what else is running; `hausdorff` and `self_intersections_per_1k` are
# deterministic in principle but accumulate float differences across
# architectures. Each gets a stated tolerance below, and the timing tolerance is
# deliberately loose -- this is a tripwire for a 2x regression, not a benchmark.
#
# # The baseline is per machine, and the file says which
#
# A timing baseline from another box is worse than none. Baselines live in
# `docs/measurements/baseline/<name>-<machine>.csv` and the machine is part of
# the filename, so a run on a different host finds no baseline and says so rather
# than comparing against numbers that never applied to it.

set -euo pipefail
cd "$(dirname "$0")/.."

MEASUREMENTS="docs/measurements"
BASELINES="$MEASUREMENTS/baseline"

# One definition of "which machine is this", shared with anything else that needs
# it -- a second copy would drift and two baselines would disagree about their own
# names (T-014).
MACHINE="${REGRESS_MACHINE:-$(./scripts/machine.sh --slug)}"
SPEC="$(./scripts/machine.sh --spec)"

# **Quoted heredoc, and the values come through the environment.** Unquoted, bash
# expands the Python body -- and a backtick in a docstring becomes a command
# substitution, so `median_ms` in a comment was executed and reported as
# "command not found" on every run. Nothing failed; it just printed noise, which
# is the kind of bug that survives for months.
export MEASUREMENTS BASELINES MACHINE SPEC
python3 - "$@" <<'PYEOF'
import csv, os, sys

MEASUREMENTS = os.environ["MEASUREMENTS"]
BASELINES = os.environ["BASELINES"]
MACHINE = os.environ["MACHINE"]
SPEC = os.environ["SPEC"]

# Columns that identify a row. Everything else is a metric.
KEYS = ("field", "algorithm", "rule", "samples")

# Per-metric tolerance as a relative fraction. `None` means exact.
#
# The structural three are exact because they are deterministic and a tolerance
# would only hide a change. The rest are stated with the reason they are not.
TOLERANCE = {
    "vertices": None,
    "triangles": None,
    "non_manifold_edges": None,
    # Deterministic in principle; float accumulation differs across
    # architectures, so a small band rather than exact.
    "hausdorff": 0.02,
    "self_intersections_per_1k": 0.05,
    # A tripwire for a doubling, not a benchmark. Wall clock moves with the
    # governor, thermal state and whatever else is running.
    "median_ms": 0.60,
    "n_cubed": None,
    # T-017's field-quality columns. **All exact.** Every one is a max, a min or
    # a ratio of counts over deterministic `libm` arithmetic -- never a sum, so
    # there is no accumulation order to differ across architectures, and a
    # tolerance would only hide the change these exist to catch.
    #
    # A field is not supposed to move at all. If one of these shifts, either a
    # field's definition changed or a constant it depends on did, and both are
    # things somebody should have to acknowledge.
    "sup_grad": None,
    "inf_grad": None,
    "eikonal_pct": None,
    "bound_gap": None,
    "certified_pct": None,
}

# Metrics where a *smaller* number is better, so only an increase is a
# regression. Everything else is compared in both directions.
LOWER_IS_BETTER = {
    "hausdorff",
    "self_intersections_per_1k",
    "median_ms",
    "non_manifold_edges",
    # A field whose gradient supremum grows is worse-conditioned than it was.
    # `eikonal_pct` and `certified_pct` are the other way round -- higher is
    # better -- and are absent here on purpose, but since every field-quality
    # column is compared exactly the direction never actually decides anything.
    # It is stated so the next person adding a tolerance does not have to guess.
    "sup_grad",
}


def rows_of(path):
    """Rows, ignoring the `#` provenance header a baseline carries (T-014).

    Every reader goes through here. The self-test did not, at first, and read the
    header as data -- so it found no row with a `median_ms` and died on a bare
    StopIteration rather than saying so.
    """
    with open(path, newline="") as handle:
        body = [line for line in handle if not line.startswith("#")]
    return list(csv.DictReader(body))


def read(path):
    rows = rows_of(path)
    keyed = {}
    for row in rows:
        key = tuple(row.get(k, "") for k in KEYS)
        keyed[key] = row
    return keyed


def compare(name, current_path, baseline_path):
    current, baseline = read(current_path), read(baseline_path)
    problems = []

    missing = set(baseline) - set(current)
    added = set(current) - set(baseline)
    for key in sorted(missing):
        problems.append(f"{name}: row vanished: {' '.join(x for x in key if x)}")
    for key in sorted(added):
        print(f"  note: new row not in baseline: {' '.join(x for x in key if x)}")

    for key in sorted(set(current) & set(baseline)):
        for metric, tolerance in TOLERANCE.items():
            have, want = current[key].get(metric), baseline[key].get(metric)
            if have in (None, "") or want in (None, ""):
                continue
            try:
                have_f, want_f = float(have), float(want)
            except ValueError:
                continue
            label = " ".join(x for x in key if x)
            if tolerance is None:
                if have_f != want_f:
                    problems.append(
                        f"{name}: {label}: {metric} {want_f:g} -> {have_f:g} (exact match required)"
                    )
                continue
            if want_f == 0.0:
                # A metric that was zero and is not is always worth failing on --
                # that is 0 non-manifold edges becoming some.
                if have_f > 0.0 and metric in LOWER_IS_BETTER:
                    problems.append(f"{name}: {label}: {metric} 0 -> {have_f:g}")
                continue
            change = (have_f - want_f) / abs(want_f)
            worse = change > tolerance if metric in LOWER_IS_BETTER else abs(change) > tolerance
            if worse:
                problems.append(
                    f"{name}: {label}: {metric} {want_f:g} -> {have_f:g} "
                    f"({change:+.1%}, tolerance {tolerance:.0%})"
                )
    return problems


def baseline_path(name):
    return os.path.join(BASELINES, f"{name}-{MACHINE}.csv")


args = [a for a in sys.argv[1:]]
accept = "--accept" in args
self_test = "--self-test" in args
args = [a for a in args if not a.startswith("--")]

if self_test:
    # Prove the checker fails when it should, by perturbing a copy rather than by
    # trusting that it would. A gate nobody has seen fail is a gate nobody has
    # tested.
    import tempfile

    source = None
    for candidate in sorted(os.listdir(BASELINES)) if os.path.isdir(BASELINES) else []:
        if candidate.endswith(".csv"):
            source = os.path.join(BASELINES, candidate)
            break
    if source is None:
        sys.exit("regress --self-test: no baseline committed yet")

    with tempfile.TemporaryDirectory() as tmp:
        rows = rows_of(source)
        fields = rows[0].keys()
        # A "deliberately slowed extractor": one row, doubled.
        target = next(r for r in rows if r.get("median_ms"))
        target_label = " ".join(
            x for x in (target.get(k, "") for k in KEYS) if x
        )
        target["median_ms"] = str(float(target["median_ms"]) * 2.0)
        slowed = os.path.join(tmp, "slowed.csv")
        with open(slowed, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(fields))
            writer.writeheader()
            writer.writerows(rows)

        found = compare("self-test", slowed, source)
        if not any("median_ms" in p for p in found):
            sys.exit("regress --self-test: a doubled median_ms was NOT detected")
        if not any(target_label in p for p in found):
            sys.exit("regress --self-test: the failure did not name the row")
        print(f"regress self-test: a doubled median_ms is detected and named -- {found[0]}")
    sys.exit(0)

names = args or [
    f[: f.rindex(f"-{MACHINE}.csv")]
    for f in (sorted(os.listdir(BASELINES)) if os.path.isdir(BASELINES) else [])
    if f.endswith(f"-{MACHINE}.csv")
]

if not names:
    print(f"regress: no baseline for machine '{MACHINE}'.")
    print(f"regress: run the benches, then: scripts/regress.sh --accept shootout")
    sys.exit(0)

failed = False
for name in names:
    current = os.path.join(MEASUREMENTS, f"{name}.csv")
    base = baseline_path(name)
    if not os.path.exists(current):
        print(f"regress: {name}: no current run at {current} -- run the bench first")
        continue
    if accept:
        os.makedirs(BASELINES, exist_ok=True)
        # The provenance goes *in* the file. A filename can hold a nickname and
        # not a compiler version, and a timing figure whose commit is unknown is
        # a rumour -- two runs a week apart on one box are not comparable if the
        # extractor changed between them (T-014).
        with open(base, "w", newline="") as out:
            out.write(SPEC.rstrip("\n") + "\n")
            out.write(open(current).read())
        print(f"regress: accepted {current} as {base}")
        continue
    if not os.path.exists(base):
        print(f"regress: {name}: no baseline for '{MACHINE}'; --accept to create one")
        continue
    problems = compare(name, current, base)
    if problems:
        failed = True
        for problem in problems:
            print(f"::error::{problem}")
    else:
        print(f"regress: {name}: no regression against {os.path.basename(base)}")

sys.exit(1 if failed else 0)
PYEOF
