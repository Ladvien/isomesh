# Agent brief — Phase 24, thirty registrations for the game

*Paste this to Claude Code in `~/isomesh`, or point it at this file. Everything below is addressed to the agent.*

---

## Mission

Work Phase 24: thirty experiments aimed at **what a player would notice**. `P-73` … `P-102` are **already registered** in `crates/isomesh/src/experiment.rs`, their prose is in `FINDINGS.md`'s Phase 24 section, and their tickets `R-073` … `R-102` are rows in `BACKLOG.md`. For each: **build the smallest harness that could refute the registered clauses, run it, and record the result — whichever way it came out.**

Six are expected to return nulls and that is registered rather than hoped. Phase 23's two most useful rows were `✗51` and `✗54`, and both said *do not build this*. **A falsified prediction is a completed task with the more valuable outcome.** A null that cost a day is a sprint you did not spend.

---

## Read first, in this order

1. `CLAUDE.md` — the project's rules. They override anything in this brief.
2. `FINDINGS.md` — Part 5 (method rules) **in full**, then the Phase 24 section. You will be adding to both.
3. `BACKLOG.md` — Phase 24's header, and check what is in flight before starting.
4. `docs/research/2026-08-27-thirty-experiments-for-the-game.md` — the source. Part 1 is Phase 23 read back and names five defects you will be asked to fix; Part 3 is the elaboration of each registration; Part 4 is the foreclosure list; Part 5 is the acquisitions.
5. `docs/research/2026-08-26-audit-and-phase-23-registrations.md` — Phase 23's own source, for the ordering rules it introduced.

---

## The registrations exist. Do not touch them.

`experiment.rs` is the source and `FINDINGS.md` quotes it. **Amending a registration after its harness exists is a rewrite of the prediction**, and the only honest way to change one is to register a new id and say why the first was inadequate. `P-43 → P-44` and `M-374`'s vacuous C3 are both precedents: *record the vacuity, let a new id carry a fixture that can fail.*

If a clause turns out to be unmeasurable, **say so in the entry and score it VACUOUS**. Do not quietly answer a different question and score it HELD — that is `P-70`'s C3, the weakest thing in Phase 23.

---

## Non-negotiables, in the order they will bite

**1. Every clause stated as a ratio must name the share it can move.** This is `✗51`'s rule and it is why that experiment was unreachable: the sample loop was 11.5% of the quantity C1 was denominated in, so halving it gave 1.06× against a registered 2×. **Before writing a harness, recompute the SHARE line in the registration and check the clause is still arithmetically reachable.** If it is not, say so *before* running — `P-70` did exactly that and the entry is better for it.

**2. Every run needs a control that could have failed.** Each registration carries a VACUITY CONTROL naming the column. Emit it. `M-44` is a zero over an unreached configuration; `P-62` added a 400,000-cell random arm because eight reference fields gave only seven tunnel cells — *"a hair from `M-44`'s vacuous zero"*.

**3. The CSV must resolve and the tree must be clean.** `scripts/csv_provenance.sh` runs in `preflight.sh`. Every Phase-23 dataset resolves to an ancestor of HEAD against a clean tree; keep that streak. Commit the harness, *then* run it, *then* commit the CSV.

**4. `Run::record` refuses a comma, a quote or a newline.** It gained that after `P-64`'s first CSV silently shifted every later column by three places. If a value needs one of those characters, the column is wrong.

**5. Numbers in an entry must come from the committed CSV.** `✗35` is this failure twice and `✗52` is it again — its headline figures live only in a superseded commit while its `M.` line names a file that no longer contains them, and its "within 0.5%" reconciliation is refuted by its own two numbers. If you supersede a run, **re-quote from the new file or say plainly which commit the old figures came from and do not claim agreement you have not computed.** `F-009` is the ticket for gating this mechanically; it is still open.

**6. Corpus presence is decided by three discriminators, not by `catalog_read`.** `M-371`: MCPro reported `markdown_path` set, conversion complete, `chunks_indexed: 1` — and its markdown is **383 characters**. Check markdown length against a real paper's (Finken's is 37,165), `chunks_indexed` against ~12, and whether `pdf_path` ends in `.pdf`. The third is exact rather than statistical.

**7. `bevy_isomesh` is excluded from the root workspace and from your pre-commit loop.** `M-293` is 58 green commits over a broken example. Eight Phase-24 rows live there. Run `cd bevy_isomesh && cargo check --all-targets` explicitly.

**8. Machine discipline.** GPU rows (`R-074`, `R-079`, `R-093`) need the Zen 3 / RTX 3090 rig. `M-280`: on a governed CPU a nanosecond is not a unit — report cycles and put the clock on the row. `M-281`: compare within one build and one run, or compare ratios. `M-005` is still blocked on a quiet machine; do not add to that pile.

---

## Order

**Do the chores first.** They are Phase 23's debt, they are cheap, and two of them are gates that would have caught the defects they exist for.

- **C1.** `✗49:10550` says "9.2 million straddling pairs"; its own `pairs` column sums to **5,800,000**. Fix the number.
- **C2.** `✗49:10611` says "2,285 of 28,124 … 7.3% of them". Both counts reproduce; the ratio is **8.12%**. Fix it.
- **C3.** `✗52` — re-quote C1 and C2 from the committed `p-71.csv`, or state which commit the figures came from; and delete or correct the "within 0.5%" sentence, which its own cited numbers refute (−5.83%, +4.55%).
- **C4.** `M-377`'s `spread` column carries gyroid's `51.045955` on the `fbm_terrain` rows. Either record the per-field value or rename the column.
- **C5.** `p-61.csv`'s `c1_rows_at_48` is not the clause's column (`c1_can_fail_rows_at_48` is). Rename so a reader who grabs the obvious one is not misled.
- **C6.** `P-70`'s C3 has no instrument. Either add a forced-`SUBGROUP`-off arm and re-score it, or re-score it VACUOUS. Do not leave a HELD with no column behind it.

**Then the work, in this order.** The first six are `S`.

1. **`R-085`** — attribute the collider's 45% with `✗52`'s instrument. **This runs before `R-081` and `R-084`**: it decides whether that cost is query-side or construction-side, and those two attack different halves.
2. **`R-089`** — granularity below 2³. An hour. Settles whether `M-377`'s 51× is a shape or a sample.
3. **`R-073`** — gradient normals against the pseudonormal. Stops a planned change rather than making one.
4. **`R-096`** — how far apart smooth union's 40,317 answers are. `M-38`'s fixture already exists.
5. **`R-098`** — the fused certificate at 0.0658. `M-378` computed the number and filed nothing.
6. **`R-102`** — `✗43`'s rate, with an artefact this time. `P-63` already built the machinery.

Then `R-093`, `R-074`, `R-081`, `R-083`, `R-099`, `R-087`, `R-092`, `R-075`, `R-077`, `R-094`, `R-091`, `R-100`, `R-101`, `R-090`, `R-086`, `R-095`, `R-080`, `R-078`, `R-088`, `R-082`, `R-097`, `R-076`, `R-079`, and `R-084` last — it is the only `L` and it should not start before `R-085` has answered.

**Three hard dependencies:** `R-085` before `R-081`/`R-084`. `R-077` before `R-076`/`R-091` — both spend a temporal budget `R-077` measures. `R-088` rides `R-087`.

---

## Per-ticket loop

1. **Re-read the registration in `experiment.rs`** — not the dossier's paraphrase. Recompute its SHARE line and confirm the clause is reachable. If it is not, write that finding first.
2. **Check the acquisition.** Part 5 of the dossier lists which paper each row needs and where an open PDF is. Apply discriminator 6 above before believing the corpus has it. `R-063` and `R-065` are blocked on exactly this failure; do not add a third.
3. **Write the harness as a bench**, `crates/isomesh/benches/experiment_pNN.rs`, with the hypothesis and falsifier in the header comment. `crates/isomesh/src/**` is read-only apart from the registrations — **except `R-083` and `R-101`, which *are* source changes and say so at registration**.
4. **Run it against a clean tree.** Commit the harness first.
5. **Write the entry** in `FINDINGS.md` under the Phase 24 section, tiered, with the verdict per clause, the numbers, the mechanism, and what changed as a result. Use the tier glyphs (🔬 measured-and-held, 💥 falsified, 🧊 measured-and-cold, 📖 verified from a source). Falsified entries take the next `✗` number.
6. **Move the ticket row to `BACKLOG_ARCHIVE.md`** with an annotation recording any amendment, deviation or falsified premise — the annotation is the point, the checkmark is not.
7. **Owe an `E×n` row in Part 4b** if you built a variant and put it back. Next free is `E×9`.
8. **`./scripts/preflight.sh --full`** before pushing.

---

## What "done" looks like for a row

A committed bench, a committed CSV whose provenance line resolves against a clean tree, a `FINDINGS.md` entry whose every number appears in that CSV, an archived ticket row with an annotation, and — if a clause could not have failed — that said out loud in the entry rather than found by the next audit.

**And one thing that is not done:** landing the mechanism. `P-70` and `P-71` both measured a real effect and did not merge it, because a second WGSL path and a doubled vertex buffer are the owner's decisions and 1.33% does not buy an exception to the one-path rule. **Measure it, surface the decision, and stop.** If a result wants a source change beyond the two named in step 3, that is a separate ticket and a separate conversation.

---

## Foreclosed — do not re-propose, and do not quietly re-enter

Part 4 of the dossier adds seven with reasons, on top of the 2026-08-23 foreclosures, the novelty table's rejected list, and the mechanics dossier's Tier-3 losers. The ones most likely to look attractive mid-task: **cluster LOD / Nanite-style virtual geometry** (≈3.5 ms to re-DAG one edited chunk, on zeux's own measured throughput); **meshlet compression** (all three strong results decode in a mesh shader, which `V-49` establishes is `unimplemented!()` on Metal in naga 29); **any VDB-family storage layer** (NanoVDB's own docs: topology cannot be modified); **Breaking Good and the learned-fracture family for carved geometry** (non-convex fragments, per-shape precompute); **fluids and granular material** (everything credible is GPU-resident and global; the gradient-projection primitive needs no paper); **Recast/Detour** (measured at a constant ~5 ms and a ~10 FPS drop on a simple sublevel); and **perceptual studies for LOD popping** (the best 2024 dataset explicitly excludes it).

If one of these looks newly viable, that is a finding and it belongs in `FINDINGS.md` with the evidence — not a quiet reopening.
