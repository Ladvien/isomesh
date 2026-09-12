# Question ledger

The state of `docs/research/2026-09-12-discovery-loop-prompt.md`. One row per question; a row moves
`OPEN` → `READY` → `REGISTERED` → `RUN`, or leaves via `BLOCKED:` / `KILLED:`. The loop reads this file
first, every iteration.

```
research_moves: 1
fields_mined: [{geometric inference / reach estimation, 1, P-177}]
last_governs_edit: 0        # iteration 0 — file created from axes-v2 + FINDINGS through M-489
next_p: 178                 # P-177 registered at iteration 1
next_r: 183                 # R-182 opened at iteration 1
```

| id | question (the assumption, as a question) | raised_by | field | terms | corpus | source_stem | state | P | expected_information | notes |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | Does this crate's topology grading rest on a regularity hypothesis — positive reach — that has never been given a number for a single reference field? | ✗124 / M-484 | geometric inference / reach estimation | `reach of a manifold`, `bottleneck`, `medial axis`, `condition number` | present, 0.708 top score, relevant top hit, **uncited in code** | `10.1214_19-ejs1551` | REGISTERED | P-177 | 7 | picked. Aamari et al. Theorem 3.4 gives the two halves without a medial axis, which is the only route left after the falsified entry 27 |
| Q2 | Is "needs a second vertex" a non-local property — sheet separation — rather than the local Hermite inconsistency M-486 measured? | M-486 | same | `bottleneck`, `reach` | folded into P-177's C3 rather than registered separately | — | KILLED: merged into P-177's C3 | P-177 | — | merged, not killed for cause: it is the same instrument on the same two fields, so a second row would be a second harness for one number. It does **not** count toward §1's three-REGISTERED threshold, which counts distinct `P-` rows |
| Q3 | Does any consumer in this repository read mesh quality at all, or is ✗125's "0 of 8 consumers" a property of the survey rather than of the crate? | ✗125 / M-488 | software architecture (no owning field of mathematics) | — | not probed | — | OPEN | — | — | §2b: no field owns it, so it is not ready. It is a survey, and a survey is a `BACKLOG.md` chore |
