# Question ledger

The state of `docs/research/2026-09-12-discovery-loop-prompt.md`. One row per question; a row moves
`OPEN` → `READY` → `REGISTERED` → `RUN`, or leaves via `BLOCKED:` / `KILLED:`. The loop reads this file
first, every iteration.

```
research_moves: 3
experiment_moves: 1
fields_mined: [{geometric inference / reach estimation, 1, P-177}, {random fields / integral geometry, 2, P-178}, {optimal transport, 3, P-179, FORCED-EXPLORE}]
last_governs_edit: 4        # iteration 4 — ✗126 / M-490 edited the Instruments section
next_p: 180
next_r: 185
forced_explore_due: 6       # §1 — every third RESEARCH move; move 3 was one
```

| id | question (the assumption, as a question) | raised_by | field | terms | corpus | source_stem | state | P | expected_information | notes |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | Does this crate's topology grading rest on a regularity hypothesis — positive reach — that has never been given a number for a single reference field? | ✗124 / M-484 | geometric inference / reach estimation | `reach of a manifold`, `bottleneck`, `medial axis`, `condition number` | present, 0.708 top score, relevant top hit, **uncited in code** | `10.1214_19-ejs1551` | REGISTERED | P-177 | 7 | picked. Aamari et al. Theorem 3.4 gives the two halves without a medial axis, which is the only route left after the falsified entry 27 |
| Q2 | Is "needs a second vertex" a non-local property — sheet separation — rather than the local Hermite inconsistency M-486 measured? | M-486 | same | `bottleneck`, `reach` | folded into P-177's C3 rather than registered separately | — | KILLED: merged into P-177's C3 | P-177 | — | merged, not killed for cause: it is the same instrument on the same two fields, so a second row would be a second harness for one number. It does **not** count toward §1's three-REGISTERED threshold, which counts distinct `P-` rows |
| Q3 | Does any consumer in this repository read mesh quality at all, or is ✗125's "0 of 8 consumers" a property of the survey rather than of the crate? | ✗125 / M-488 | software architecture (no owning field of mathematics) | — | not probed | — | OPEN | — | — | §2b: no field owns it, so it is not ready. It is a survey, and a survey is a `BACKLOG.md` chore |
| Q4 | Is `M-458`'s "first ground-truth χ" for `noise_cavity` a ground truth at all, or two samplings of a quantity that has not converged — 82 at 33³ against −152 at 65³? | M-458, M-489 | random fields / integral geometry | `Gaussian kinematic formula`, `excursion set`, `expected Euler characteristic`, `Lipschitz–Killing` | present, 0.687 top score, Taylor 2006 top hit, **uncited in code** (`Lipschitz–Killing`, `Hermite polynomial`, `spectral moment` all 0 in the repo) | `10.1214_009117905000000594` | RUN | P-178 | 8 | **answered: not a ground truth.** ✗126 / M-490 — 82 / −152 / −186 at 33³/65³/129³, still moving. Row is VACUOUS as registered because the `4λ` control is unsatisfiable at `u = ±1` |
| Q5 | Does curvature error measured against an analytic curvature assume normal convergence the mesh never proved — the Schwarz lantern's objection applied to `M-487`? | M-487 | geometric measure theory (varifold convergence) | `normal convergence`, `Schwarz lantern`, `varifold` | not probed this move | — | OPEN | — | — | candidate 2 of 3. Kept open deliberately: it is the next question if `P-178` lands, and it shares `M-487`'s instrument |
| Q6 | Is the isovalue sweep in `M-489` measuring the field or the sampler — does the percolation threshold move with resolution? | M-489 | percolation theory | `critical level`, `finite-size scaling` | not probed this move | — | OPEN | — | — | candidate 3 of 3. Closest to a tune (a re-run of an existing instrument at other resolutions), which §2d gate 4 rejects, so it is parked rather than registered |
| Q7 | Does "accuracy" in this repository assume an exact distance oracle, when the object being graded is a measure and a transport distance needs only samples? | ✗107 / M-459 | optimal transport (FORCED-EXPLORE: named before the hook) | `sliced Wasserstein`, `optimal transport`, `1-D projections` | present, 0.721 on the targeted query, Peyré & Cuturi §10.4 the top hit; `Wasserstein`, `earth mover` and `Sinkhorn` are **0** in the repo and optimal transport appears only in research prose as *"uncited, in corpus"* | `10.48550_arXiv.1803.00567` | REGISTERED | P-179 | 5 | the three registered rows now exceed §1's threshold, so the next move is an EXPERIMENT and it takes `P-178` at 8 |
| Q8 | Why does `fbm_terrain` land inside the Gaussian kinematic band at 0.25σ while `noise_cavity` misses by 6.77σ with the sign of χ wrong, when both are hash noise? | ✗126 / M-490 | random fields — **stationarity**, not Gaussianity | `stationary random field`, `non-stationary excursion set`, `local isotropy` | not probed | — | OPEN | — | — | raised by the EXPERIMENT move of iteration 4. First candidate: `noise_cavity` carries a spherical cap of radius 1.5, and a stationary-field theorem integrated over the box does not allow one |
