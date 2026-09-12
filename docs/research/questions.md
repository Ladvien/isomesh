# Question ledger

The state of `docs/research/2026-09-12-discovery-loop-prompt.md`. One row per question; a row moves
`OPEN` → `READY` → `REGISTERED` → `RUN`, or leaves via `BLOCKED:` / `KILLED:`. The loop reads this file
first, every iteration.

```
research_moves: 0
fields_mined: []            # [{field, iteration, P}]
last_governs_edit: 0        # iteration 0 — file created from axes-v2 + FINDINGS through M-489
next_p: 177                 # highest registered in crates/isomesh/src/experiment.rs is P-176
next_r: 182                 # highest in BACKLOG.md / BACKLOG_ARCHIVE.md is R-181
```

| id | question (the assumption, as a question) | raised_by | field | terms | corpus | source_stem | state | P | expected_information | notes |
|---|---|---|---|---|---|---|---|---|---|---|
