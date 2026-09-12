# Kill list — fields and mechanisms probed and found empty

A field lands here when a RESEARCH move probed it and found no hook into the ledger, or found the corpus
empty with no reachable source. **Never re-probe a field on this list** unless a new `FINDINGS.md` entry
names it in its `Raises:` line — that is the only reader in the loop with standing to overturn a kill.

Each row records the date, the probe's `distill_search` top score (calibration: present ≈ 0.65–0.72 with a
relevant top hit, absent ≈ 0.55–0.64 with an unrelated one), and the reason.

| field | date | probe score | reason | reopened by |
|---|---|---|---|---|

**Empty at creation, 2026-09-12.** Nothing has been killed yet. Two adjacent facts, so the first iterations
do not rediscover them:

- **Anisotropic mesh adaptation is not killed, it is *graded and negative*** — `✗107 / M-459`, where the
  metric-driven arm was never cheaper and 2.298935× dearer on `thin_plate`. A new hypothesis in that field
  must say what it does differently, not re-run it.
- **Intrinsic Delaunay is not killed either, it is *invisible*** — `✗125 / M-488`, `vertex_positions_moved`
  **0** and `c3_holds` **0 of 8** consumers. A hypothesis there must name a consumer that reads intrinsic
  connectivity, or it has no hook.
