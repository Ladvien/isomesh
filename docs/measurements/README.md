# Measurements

Every number this project publishes comes from a committed benchmark, per `CLAUDE.md` rule 4. This
directory is where those benchmarks write, and this file is the procedure for taking a measurement
that someone else can use — which is a different and stricter thing than taking one.

## What is here

| File | Written by | What it is |
|---|---|---|
| `shootout.csv` | `cargo bench --bench shootout` | Every extractor × eight fields × three resolutions: timing, counts, manifoldness, self-intersection, Hausdorff |
| `ablation.csv` | `cargo bench --bench ablation` | One algorithm, two vertex rules — X-002's seam |
| `resolution_sweep.csv` | `cargo bench --bench resolution_sweep` | 16³…256³, and the `t = a + b·n³` fit |
| `stage_breakdown.csv` | `cargo bench --bench stage_breakdown` | Where the time goes inside one extraction |
| `gpu_vs_cpu.csv` | `cargo bench --bench gpu_vs_cpu` | The GPU path against the CPU, part by part |
| `baseline/` | `scripts/regress.sh --accept` | Committed reference runs, one set per machine |

`resolution_sweep-ryzen9-5900x.csv` predates this convention. It is a hand-named second-machine run
and is left alone deliberately: **M-45, ✗14 and O-11 quote figures from it by name**, so renaming it
to match would silently break five references to buy tidiness.

## Taking a measurement

1. **Release, and nothing else running.** Debug meshing is 37–62× slower (M-152) and the numbers are
   not merely scaled, they are a different shape.
2. **Commit first.** `scripts/machine.sh` records the commit into the baseline header, and it marks
   the run *not attributable* if anything under `crates/`, `bevy_isomesh/` or the manifests has been
   modified since. A doc or script edit does not trip it, because neither can change a triangle count.
3. Run the bench. It writes its CSV here, overwriting the previous run — that file is scratch.
4. `scripts/regress.sh` to compare against the committed baseline for this machine.

## Checking for a regression

```bash
scripts/regress.sh                      # every baseline this machine has
scripts/regress.sh shootout             # one
scripts/regress.sh --accept shootout    # adopt the current run as the new baseline
scripts/regress.sh --self-test          # prove the checker still detects a regression
```

**Two classes of column, and only one is noisy.** `vertices`, `triangles` and `non_manifold_edges`
are deterministic — T-004 guarantees it — so they are compared **exactly**, and a tolerance there
could only hide a real change. `hausdorff` gets 2% and `self_intersections_per_1k` 5%, for float
accumulation across architectures. `median_ms` gets **60%**: a tripwire for a doubling, not a
benchmark.

`--self-test` is what runs in CI. The benches take minutes and their timings are meaningless on a
shared runner, so what CI protects is the *checker* — the half that can rot without anyone noticing,
because a gate that has stopped detecting anything still exits zero.

## Cross-machine

**A baseline from another box is worse than none**, so baselines are per machine:
`baseline/<bench>-<machine slug>.csv`. A run on an unrecognised host finds no baseline and says so
rather than comparing against numbers that never applied to it.

The slug is short enough to be a filename and therefore cannot carry provenance, so the provenance
goes **inside the file** as `#` lines that every reader here skips:

```
# machine: AMD Ryzen 9 5900X 12-Core Processor, 24 logical cores, 31 GB
# system:  Linux 7.1.3-1-cachyos
# rustc:   rustc 1.96.0 (ac68faa20 2026-05-25)
# commit:  9568b82
# taken:   2026-08-15T21:16:42Z
```

The commit line is the one that matters most and the one most easily lost: two runs a week apart on
the same box are not comparable if the extractor changed between them, and nothing else in a CSV
records that it did.

**A second machine is worth the trouble, and there is a finding to prove it.** M-45 exists only
because the resolution sweep was run on a second box: Surface Nets' superlinearity reproduced there
and got *worse*, which is what ruled out one cache hierarchy as the explanation. A measurement taken
on one machine describes that machine until a second one agrees.
