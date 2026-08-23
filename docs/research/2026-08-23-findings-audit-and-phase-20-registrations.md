# The ledger audited against its own artefacts, and eight registrations that follow from what the audit found

**Date:** 2026-08-23 · **Repo state:** `61c6201` ("Working"), 431 entries, 442 commits
**Scope, as chosen:** the load-bearing rows — the ones cited elsewhere as justification — plus every
Phase 19 entry, because Phase 19 is where the ledger and the repository have parted company.

---

## 0. What I did, and one provenance correction I owe you

I checked entries against three things, in this order: the **committed CSV the entry names**, the
**source it makes claims about**, and the **primary literature it cites**. Where I could re-run
something deterministic I did — the whole `isomesh` suite passes on Linux/x86-64, 707 + 19 + 1 tests, at `61c6201`,
golden hashes included, which is M-31 re-confirmed on a third architecture rather than assumed.

**The correction.** My first pass read the repository at `59c641c`, reported that Phase 19 was absent,
and asked you a scoping question on that basis. That was true of the snapshot and false thirty minutes
later: `61c6201` landed at 13:46 and brought P-38…P-47, M-336…M-342, ✗28…✗31, ten benches and ten
CSVs with it. The finding I built on a stale read is withdrawn. It is worth one sentence anyway,
because the same failure mode is the subject of §1: **an artefact and a claim disagreeing is not
always the claim's fault, and the first move is to check which one moved.**

---

## Part A — Verification

### A.1 The rows that reproduce exactly

These I checked digit by digit against the committed CSVs. No corrections.

| Row | Artefact | Result |
|---|---|---|
| **✗1 / M-2 / M-22** — `V_sn = V_mc + χ`, `F_sn = F_mc + 2χ` | `family.csv` | **Exact at all nine resolutions.** Vertices +2 and triangles +4 on every row from 16³ to 256³ (76,776→76,778 and 153,548→153,552 at the top). The identity is not approximately true; it is true. |
| **✗25** — three occupants of the good corner | `shootout.csv`, 112 rows | **Exact.** MC 0/0, MC+decider 0/0, MT 0/0; SN 747 nm-edges / 0 self-int; DC 747 / 29.745; MDC 90 / 45.232; subgrid 446 / 0. The 90 are all `noise_cavity` (64 + 26). Every non-zero self-intersection row is at 33³ and none at 65³, exactly as the entry says. |
| **M-308** — ✗14's cost half overturned | `family.csv` | **Exact.** Ratios 0.933 / 0.958 / 0.981 / 1.030 / 1.101 / 1.161 / 1.291 against the entry's 0.93 / 0.96 / 0.98 / 1.03 / 1.10 / 1.16 / 1.29. The crossover sits between 32³ and 48³, so "about 40³" is honest. |
| **M-338 / P-41** — the critical-cell bijection | `p-41.csv`, 16 rows | **Exact, including the null control.** gyroid 132+9 = 141 critical = 141 nm-vertices = 141 hosting; fbm_terrain 58 = 58 = 58; noise_cavity 567+35 = 602 = 602 = 602. Co-location 1.000000 on every affected row; the chance baseline 0.66% / 0.69% / 2.1%. `box_exact` Surface Nets really does put 5,400 of 5,768 vertices on a cell boundary. |
| **M-340 / ✗30 / P-42** — the curvature bound | `p-42.csv` | **Exact, to every digit I checked.** Residual 2.81e-13 → 4.95e-12; ratios 4.013 and 4.394; residual-per-vertex 2.424e-16 / 2.367e-16 / 2.573e-16 against `ε = 2.22e-16`; `ε/h` = 1.398 / 1.593 / 3.608; fatness 4.164e-2 → 1.210e-2; the excluded-degenerate bound identical to six places. The mean-curvature relative errors halve-squared cleanly (0.249, 0.250). |
| **M-341 / P-39** — pruning is free and bit-exact | `p-39.csv`, 64 rows | **Exact but for one word.** Median survivor 0.2969, speedup median 3.3648, world 2.4726, min 0.9923, max 22.4695, 56 of 64 over 1.25×, 3 all-survive, 48 smooth chunks differing, `dominant_adds` 0, bound cost 4.4e-5–1.3e-3, bound 540–1450 ns. **One error:** the entry says *"seventeen share 173"*; the CSV has **nineteen**. |
| **M-342 / ✗31 / P-44** | `p-44.csv` | **Exact.** Pearson 0.980 / 0.982 / 0.99994 / 0.997; gaps 1.008 / 0.888 / 0.019 / 0.127; cost 0.563–1.262; witness exponents 0.968–0.999; and the sphere maximum bit-identical at `2.558110331e-1` on all four rows, which I re-derived independently: `|√3/2 − (3+3√2+√3)/8| = 0.2558110331`. Seven digits. |

**M-31 re-confirmed as a side effect.** The 216 golden hashes were generated on macOS/arm64 and pass
unchanged on this Linux/x86-64 container. That is a third platform, not a second.

### A.2 The two rows whose numbers are not in the repository

This is the substantive finding, and it is one failure repeated twice.

#### M-336 / ✗28 (P-38) — the entry's table is not `p-38.csv`

| row | entry: before / neighbour mean / ratio | `p-38.csv`: before / mean / ratio |
|---|---|---|
| MC 128 | 8.643 / 8.850 / **0.977** | 8.7419 / 8.7539 / **0.9986** |
| MC 256 | 8.570 / 8.603 / **0.996** | 8.5094 / 8.4618 / **1.0056** |
| SN 128 | 11.493 / 11.879 / 0.968 | 10.9271 / 11.0748 / 0.9867 |
| SN 256 | 12.047 / 12.842 / 0.938 | 11.8755 / 12.7414 / 0.9320 |

The surface-free control and the `pad_z` control disagree too — the entry says 8.116 against 8.046 and
a `pad_z` ratio of 1.0075; the CSV says 8.1792 with `scaffold_ratio_before` 1.0138 and `pad_z_ratio`
1.0038.

**The verdict survives and the numbers do not.** Every reading in both the entry and the artefact sits
between 0.93 and 1.01 against a registered threshold of 1.5×, so ✗28 — *Marching Cubes does not alias
at 128³, and the access pattern is why* — is not in question. What is in question is that the entry
cites `docs/experiments/p-38.csv` and quotes a run that is not in it. The CSV's own header says
`commit 9abc62f`, the harness commit; the entry was written at `64bef8a`.

#### M-337 (P-40) — same shape, and here it changes a verdict

The `before` column is a **committed fixture** (`docs/measurements/p40-baseline.csv`), so it matches
the entry to five digits, as it must. The `after` column is timed live, and half of it does not match:

| sphere row | entry after / ratio | `p-40.csv` after / ratio |
|---|---|---|
| SN 64 | 7.227 / 1.372 | 7.2148 / 1.3745 ✓ |
| **SN 128** | **7.960 / 1.336** | **8.9173 / 1.1925** ✗ |
| SN 256 | 8.747 / 1.247 | 8.9700 / 1.2160 ✗ |
| DC 64 | 21.698 / 1.122 | 21.7488 / 1.1195 ✓ |
| DC 128 | 14.576 / 1.179 | 14.6225 / 1.1755 ✓ |
| DC 256 | 11.855 / 1.184 | 12.0548 / 1.1642 ✗ |

Stage ratios likewise: the entry reads 5.47 / 5.49 / 5.34 / 5.31 / 5.49 / 5.33; the CSV reads 5.3455 /
5.3655 / 5.2468 / 4.9076 / 5.1535 / 5.2102.

**M-337 records C2 as "HELD for Surface Nets" on 1.336× against a registered ≥ 1.25×. The committed
artefact says 1.1925×, which fails the bar — and so does SN at 256³, at 1.2160×.** By the CSV in the
repository, C2 failed for both extractors, not one.

I am not claiming the mechanism is worth less. C1 (stage ratio ≈ 5.3×, deterministic in shape) and C3
(12 of 12 mesh hashes identical, deterministic outright) both hold in the artefact, and the bitmap
prepass is plainly a real win. What failed is the *clause*, and it failed because of what it was:

> **✗24 is the rule that already covers this**, and Phase 19 did not apply it. *"A wall-clock ratio is
> not a gate. Gate the count the ratio samples."* That entry cost a release. P-40's C2, P-38's whole
> table, and P-47's C2 are all wall-clock ratios registered as pass/fail thresholds, on a machine whose
> governor M-280 already caught swinging 1.8×.

**P-54 below is the fix**, and it is the cheap one: the quantity the bitmap prepass actually changes is
*how many eight-corner gathers run*, which is `active_cells` against `cells` — two integers, identical
on every machine, already in the CSV.

### A.3 Three experiments that ran, produced results that falsify them, and have no ledger entry

`p-45.csv`, `p-46.csv` and `p-47.csv` are committed with data. `FINDINGS.md` carries P-45, P-46 and
P-47 only as registrations. The index agrees — 333 measured, no M-343. **The results are in the
repository and the conclusions are not**, which is the state Part 5 already has a rule about.

I read the CSVs. All three registrations are falsified, and two of them interestingly.

**P-45 (R-041a, curvature additivity) — C1 and C3 falsified, and the artefact names the cause.**

| field | global Gaussian | chunk sum | gap | gap / π | excess incidence | Borel gap |
|---|---:|---:|---:|---:|---:|---:|
| sphere | 12.566370614359 | 351.858377202056 | 339.29 | **108.000000** | **108** | 1.12e-12 |
| torus | 0.000000000000 | 647.168086639496 | 647.17 | **206.000000** | **206** | 1.00e-12 |
| box_exact | 12.566370614359 | 351.858377202057 | 339.29 | **108.000000** | **108** | 3.55e-15 |

`gaussian_gap_matches_excess` is `true` on all three. The boundary term as transcribed counts a vertex
once **per chunk it touches**, and each extra incidence costs exactly π. The **Borel** accounting —
each vertex to exactly one chunk — closes to 1e-12, so additivity itself is fine and the transcription
is not. C2 (mean measure) holds at 1.7e-13–3.0e-13. C3 fails asymmetrically and informatively:
`isolated_gaussian_bit_identical` is **true** with zero mismatched chunks, `isolated_mean_bit_identical`
is **false** on 56 of 64 chunks on sphere and 64 of 64 on torus, and **true** on `box_exact`. So the
**Gaussian** measure is chunk-local bit-for-bit and the **mean** measure is not — which is exactly what
you would predict, since a dihedral needs both faces of a seam edge and a chunk has one. `box_exact`
passes because its seam edges are flat (`flat_edges` 16,548 of 17,292) and a flat edge contributes zero
either way.

**P-46 (R-040a, well-composed repair) — falsified on all three clauses, on one field.**

| field | critical → | sweeps | nm edges → | nm vertices → | Hausdorff ratio | stuck cells |
|---|---|---:|---|---|---:|---:|
| gyroid | 141 → **0** | **5** | 69 → 0 | 141 → 0 | 1.000000 | 0 |
| fbm_terrain | 58 → **0** | 1 | 29 → 0 | 58 → 0 | 1.000000 | 0 |
| noise_cavity | 602 → **118** | **64** | 322 → **70** | 602 → **118** | **4.45** (DC) / 1.69 (SN) | **6,307** |

C1 fails (118 survive), C2 fails on two of three fields (5 and 64 sweeps against a registered ≤ 2), C3
fails at 4.45× against ≤ 1.10×. And the field it fails on is `noise_cavity`, which exists precisely
because it is the only reference field with an interior ambiguity (M-208).

**The registration's stated interpretation of its own falsifier does not apply, and this is the part
worth writing down.** P-46 says a survivor *"would mean well-composedness of the sign lattice is not
sufficient"*. That inference needs the repair to have **reached** a well-composed lattice, and it did
not — 118 critical configurations remain, `stuck_cells` is 6,307, and `residual_exhausted_cells` is
118. Sufficiency is untested. What *is* tested, and holds: `residual_nm_vertices_in_critical` = 118 =
the residual critical count, so **M-338's bijection survives the repair**, which is a second,
independent confirmation of it on perturbed data.

The other thing the artefact says plainly: on `noise_cavity` the repair moved 5,190 corners (1.9% of
samples), took triangles from 57,764 to 65,272 (+13%), and moved Hausdorff from 0.253 to 1.126.
"Smallest representable step" is not small in aggregate.

**P-47 (R-043, dual-number normals) — C1 falsified by three orders of magnitude; C2 and C3 hold.**

| fixture | mean angular error | max | speedup | control |
|---|---:|---:|---:|---:|
| brush_stack_64 | **7.598e-5°** | **4.365°** | **2.84×** | 2.83e-16 |
| brush_stack_64_smooth | 7.650e-5° | 4.365° | 2.49× | 2.83e-16 |
| capsule | 7.0e-10° | 1.7e-9° | **0.877×** | — |
| sphere (analytic override) | 0 | 0 | 0.863× | — |

Registered: mean **> 0.1°** and max **> 5°**. Measured 7.6e-5° and 4.365°. Both miss, and the
registration's own falsifier called it — *"a real possibility, since `DIFF_STEP` scales with `|p|` and
these shapes are smooth away from their seams."* C2 holds at 2.84×, C3 at 2.83e-16, and
`dual_value_bit_exact` is `true`, which is the control that says the dual arithmetic reproduces the
field.

**The residue is much sharper than the hypothesis was.** `bulk_mean_angular_error_deg` is **1.9e-8°**
while the reported mean is 7.6e-5° — so the mean is carried by **one vertex in 57,470**.
`vertices_over_1deg` = 1, `vertices_over_5deg` = 0, and `worst_stencil_straddles_seam` = **true** from
32 brushes upward. Central differences are effectively exact on a composed field *everywhere except
where the six-sample stencil straddles a CSG seam*, and there they are off by 4.4°. That is a finding;
it is not the finding that was registered. **P-55 below turns it into one that can fail.**

### A.4 One gate is red at HEAD

```
::error::README.md:263 claims "36 examples" — there are 38
::error::bevy_isomesh/README.md:149 claims "36 examples" — there are 38
doc facts FAILED
```

`61c6201` added `critical_cells` (E-304) and did not move the count. `scripts/doc_facts.sh` exists
because this number has rotted twice (M-295); it caught the third and the commit went in anyway.
`findings_index.sh --check` and `backlog_gate.sh` are both green.

### A.5 One small numeric slip

✗29 says the inadequate-cell fraction runs *"3.8–19.4% on every registered row"*. Measured from
`p-43.csv`: gyroid 2.93 / 5.90 / 5.28 / 3.76, noise_cavity 14.43 / 19.38 / 18.29 / 16.17. The low end
is **2.93%**, not 3.8%. The claim it supports — never zero — is unaffected.

---

## Part B — What the citations actually say

I checked thirty identifiers against Crossref, arXiv and publisher copies. Twenty-three are clean.
Seven need something, and four of those matter.

### B.1 Corrections that change what a row may claim

**1. "Cohen-Steiner & Morvan's Theorem 6" does not exist.** `10.1145/777792.777839` (SoCG 2003) has
**one** numbered theorem, Theorem 1, plus Proposition 3 and Lemmas 2 and 4–9. The `C_S·K` bound
numbered **Theorem 6** is in **Sun & Morvan, `10.5802/acirm.50`** — which the ledger already cites, and
already correctly notes is the in-corpus one. P-42 and M-340 attribute the theorem to the wrong paper
of the two they cite. The substance is untouched; the attribution should swap.

**2. The FlexiCubes pointer is wrong and the PyVista attribution is unsupported.** V-38 says the AR>4
figure is in *"Fig. 15"* — it is **Table 4**. The `MC + Reg.` labelling that V-38 exists to correct is
right (plain `MC` is 11.71 in the same column), so the finding stands. But V-38 also says *"FlexiCubes'
Fig. 15 names **PyVista** as the measurement tool, which wraps VTK's Verdict library, so the definition
to match is Verdict's"* — and **PyVista appears nowhere in the paper**, which credits no tool and says
only *"we compute triangle aspect ratios, radius ratios, and min and max angles."* V-38 then instructs
*"Confirm it from Verdict before implementing — rule 5."* That instruction was built on an attribution
that is not in the source. V-39 superseded the metric anyway, so nothing downstream is wrong; the row
needs the line struck.

**3. Boutry, Géraud & Najman `10.1007/s10851-017-0769-6` is *A Tutorial on Well-Composedness* — a
survey, not the repair method.** P-46 calls its mechanism *"Boutry's self-dual repair
(`10.1007/s10851-017-0769-6`)"*. The self-dual repair is **`10.1007/978-3-319-18720-4_47`**, *How to
Make nD Functions Digitally Well-Composed in a Self-dual Way*, ISMM 2015, same three authors. Same
family, different paper. Given P-46's results this matters more than a citation tidy: it is worth
reading the actual method before concluding the repair does not converge, because **the method being
cited is not the method that was implemented**.

**4. `‖∇ρ‖ = cos(θ/2)` is not in the Attali/Boissonnat/Edelsbrunner report.** ✗27 and P-27 present it
as read from `10.1007/b106657_6`. Two passes over the freely available PDF find eq. (1)
`∇ρ(x) = (x − c(x))/ρ(x)` and, separately, `ρ = δ/sin(θ/2)`. The identity **follows** from eq. (1) — it
is a correct derivation — but it is a derivation, presented as a quotation. That is precisely the
family ✗21's amendment catalogues: *"a property lifted from a summary, not from the thing itself."*
Nothing measured changes; ✗27 and M-324 stand on their own numbers.

### B.2 Softer ones, recorded so they are not re-checked

- **P-29's "verbatim" overclaims.** `k₁ = 4×10⁻¹¹` is verbatim in `10.5194/hess-23-1995-2019`.
  `c_eq = 10⁻⁶ mol cm⁻³` and `D = 10⁻⁵ cm² s⁻¹` were **not found** in the open-access full text over
  three passes — only the composite relation. They are Dreybrodt's habitual values elsewhere. The
  registration says *"Everything the simulator needs is in it, verbatim"*; two of three constants are
  not.
- **P-36's Dziuk & Elliott quote is unverifiable.** `10.4171/ifb/182` is the right paper on the right
  subject, but *"the matrices depend only on the evaluation of the gradient of the level set function"*
  is not in the abstract and the body is paywalled. It is cited as READ.
- **Custodio §5.1's "second order, two sign changes"** (V-24) could not be checked — paywalled, no open
  copy. The correction itself is corroborated by the same group's 2019 JBCS paper.
- **Lieutier & Wintraecken `10.1145/3564246.3585113` is STOC 2023**, not SoCG. Exponents ½ and ¼
  confirmed exactly via arXiv:2303.04014, Props. 6.3 and 6.5.
- **Barbier's "two orders of magnitude" is a sphere-tracing speedup** (629× on a 6,023-node tree), not
  a primitive-count reduction. P-39 attributes the phrase to **Keeter**, whose abstract does say it, so
  P-39 is clean — but do not let the Barbier number migrate into a meshing row.
- **The near-miss DOIs the ledger flags are both correct.** `10.1006/cviu.1995.1013` really is
  Brechbühler's SPHARM paper and `10.1016/j.dam.2015.01.006` really is Turán numbers and batch codes.

### B.3 Confirmed verbatim, in case you want the certainty

Whiting's `0.1075` and `15.84°`; Schaefer/Ju/Warren's *"this surface is always a manifold…"* (✗19's
target, word for word bar the MC abbreviation); Grosso & Zint's mean-ratio formula and the MC/TMC
two-decimal agreement on **the average-quality column** (they are *not* identical in irregular-vertex
%, e.g. Skull 59.41 vs 60.17); van Gelder & Wilhelms' ambiguous-face definition; VisACD's 16.97 s /
36.31 s and the 35% merging figure; CoACD's 194.4–253.4 s single-threaded; Barbier's
polygonization-as-future-work sentence; and **`arXiv:2606.00454` is real** — Baktash, Gillespie &
Crane, *Subgrid Marching Tetrahedra*, TOG/SIGGRAPH 2026, `10.1145/3811358`, with §3.1, §3.2.1–3.2.3
and Appendices A–E as A-014's series describes them.

**One thing that is now worth doing and was not possible before:** the paper is published with
**expanded lookup tables and reference C and JavaScript implementations in supplemental**. A-014
derived those tables under rule 5 with the paper unobtainable. There is now a primary source to diff
against, and CLAUDE.md's rule 5 argues for doing it.

---

## Part C — Phase 20, registered before any harness exists

Eight, ranked. Each states what it predicts, what would refute it, the arithmetic derived **before**
the run, and the control that makes a null mean something. Two are repairs of clauses this audit
broke; six come from the literature sweep, and every one of those names the specific reason its
headline number might not transfer.

Sources are marked **[corpus]** where `docs/research/` already carries them and **[acquire]** where
they are not in the library.

---

### P-48 — the empty ball every sample asserts, and whether this crate's extractors respect it

**The idea, and why it is not a metric nobody thought of.** A signed distance sample `(p, d)` does not
only say *the field is `d` here*. It says **the open ball `B(p, |d|)` contains no surface** — that is
what a distance means. Marching Cubes, Surface Nets and Dual Contouring all read `d` as a number to
interpolate and discard the ball. Four 2025–2026 papers are built on recovering it — power diagrams
weighted by `d²` (`10.1111/cgf.70037`), polyhedral empty spheres (`10.1145/3721238.3730748`), gradient
approximation (`10.1111/cgf.70373`), regular triangulations for UDFs (`10.1111/cgf.70524`) — and every
one needs a weighted Delaunay construction that breaks rule 3, is not `no_std`, and is not per-cell
parallel. **[acquire]**

**The cheap half needs none of that.** Whether the extractors *violate* the constraint is one pass over
the output, and it decides whether the expensive half is worth reading further.

> **H.** With `violation(v) = max over samples p within one cell of (|d(p)| − ‖v − p‖)`, normalised by
> cell size: **(C1)** on the fields that declare `FieldBound::Exact` — `sphere`, `torus`, `box_exact`,
> `thin_plate`, `csg_difference` — `marching_cubes` violates on **fewer than 1 vertex per 1,000** at
> 65³, because a Marching Cubes vertex is the root of the interpolant on a grid edge and cannot be
> deeply inside a neighbour's ball. **(C2)** `dual_contouring` violates on **at least 20 per 1,000** on
> the same fields, because the QEF minimises distance to tangent *planes* and the cell clamp is the
> only thing holding it, so it will sit inside a ball whenever a neighbouring sample is closer to the
> surface than the plane fit believes. **(C3)** `surface_nets` sits strictly between them.

**Falsified by** C2 coming in under 20 per 1,000, which would say the constraint is already respected
by construction and closes four papers for this crate in one pass. **Or** by C1 exceeding C2, which
inverts the mechanism. Neither outcome is a wasted ticket, and C2 is the one I expect to be
interesting: M-27 measured that on `box_exact` **864 of 1,016 Dual Contouring vertices agree with the
centroid to 2e-15 and 152 move by 0.35–0.57 cells** — those 152 are exactly the population that can
enter a ball.

**Registered arithmetic, so a null cannot be blamed on tolerance.** The floor is not zero: a vertex is
placed by interpolation and a sample's own `|d|` carries the field's discretisation, so
`|d(p)| − ‖v − p‖` can be positive by `O(h²)` on a curved surface with no rule violated. The threshold
is therefore **`0.05 · h`**, derived from M-12's measured `h²` convergence (`sphere` at 65³ measures a
mean error of 6.5e-4 against `h = 0.0635`, i.e. **1.0% of a cell**), which puts the gate five times
above the honest floor and far below the 0.35–0.57 cells M-27 measured for a moved vertex.

**The control that must run first.** On `box_exact`, whose surface is planar and axis-aligned and where
M-27 measured *every* Dual Contouring vertex on a planar patch landing exactly on the centroid, the
violation count for **all three** extractors must be zero. An instrument that finds violations there is
measuring its own tolerance.

**Records** `field`, `extractor`, `samples_per_axis`, `violations_per_1k`, `worst_violation_cells`,
`vertices`, `samples_probed_per_vertex`, `control_box_exact_zero`.

**Why it might not transfer, stated up front.** The empty-ball premise needs the field to *be* a
distance. A voxel game's field after boolean edits is not one — F-004/M-247 measured the worst-case
underestimate ratio falling from 0.577 to 0.004 over 256 strokes. So this is a claim about the five
`Exact` fields and about the CAD consumer, and it must never be run on `gyroid`, `fbm_terrain` or
`noise_cavity` and reported as a defect.

---

### P-49 — the tangency energy, on the vertex-rule seam that already exists

**The claim being tested.** *Dual Contouring of Signed Distance Data* (arXiv:2604.00157, SIGGRAPH 2026,
`10.1145/3799902.3811116`) replaces the QEF's plane-distance rows with a **tangency** residual: each
sample is a sphere of radius `|s_j|` the surface must touch, and the linearised energy is
`((t_j − q_j)·d_j)²`. Reported on ABC at 100³: Chamfer `0.779e-3` against MC's `2.62e-3` and
gradient-estimated DC's `1.589e-3`, and **edge Chamfer 0.0262 against MC 0.417 and DC 0.350 — 13× on
sharp features**. **[acquire]**

**Why this crate can test it in an afternoon and almost nobody else can.** X-002 built the seam
deliberately: `VertexRule::place(sdf, corner, base, origin, cell_size, out)` at `dual.rs:135`, with
`Qef` and `Centroid` already running through *identical* cell classification, quad walk and buffers.
M-237 measured that swap and pinned the property that makes it a controlled comparison — vertex,
triangle and non-manifold counts **identical between arms on every field**, 680 vertices with
byte-identical index buffers and all 680 positions different. A third rule inherits all of that.

> **H.** With `Tangency` as a third `VertexRule`, at a fixed **two** inner iterations: **(C1)** on
> `box_exact` and `thin_plate` at 65³ the symmetric Hausdorff improves by at least **1.5×** over `Qef`;
> **(C2)** on `sphere` and `torus` it is within **±10%** of `Qef` — the sharp-feature rule must not buy
> creases by losing smooth surfaces; **(C3)** vertex, triangle and non-manifold-edge counts are
> **identical** to `Qef` on every field, which is what says the rule changed placement and nothing else.

**Falsified by** C1 under 1.5×, which says the tangency energy does not survive being clamped to its
cell and iterated twice; **or** C3 failing, which means the arm is not controlled and no number in it
means anything.

**The two ceilings this must be read against, both measured here.** M-315 measured that projecting every
vertex exactly onto the true surface buys only **1.5%–21.5%** of symmetric Hausdorff on `sphere` and
`torus`, so **C2 is a no-harm clause and C1 is where any win must live** — and M-315 also measured that
Dual Contouring's Hausdorff is **vertex-dominated on 8 of 8 rows** while Marching Cubes' is
centroid-dominated, so a placement rule is aimed at the right extractor. Second: the paper's headline
costs *~10 s per 100³ mesh at ~100×100 iterations*. Two iterations is not that method; it is the
question of whether the first two iterations carry the sharp-feature gain. If C1 fails at two and the
paper is right at a hundred, the honest conclusion is that this is offline CAD and not a chunk budget,
and that belongs in the ledger as a null with a number.

**Records** `field`, `samples_per_axis`, `rule`, `iterations`, `symmetric_hausdorff`,
`hausdorff_ratio_vs_qef`, `edge_chamfer`, `self_intersections_per_1k`, `vertices`, `triangles`,
`non_manifold_edges`, `counts_identical_to_qef`, `ns_per_sample`.

---

### P-50 — the third corner label, on the data that actually hits it

**This one targets a defect the crate has already measured on real volumes.** Custódio et al.'s
*extended triangulation* (`10.1186/s13173-019-0086-6`, in the same series as V-24 and V-25) observes
that MC33 classifies a corner **equal** to the isovalue as inside, which marks all three incident edges
cut and emits triangles with coincident vertices. Their fix is a third vertex label — `+`, `−`, `=` —
which removes the degeneracy while preserving topological correctness across all 33 cases. **[acquire]**

**On analytic fields this is measure-zero. On this crate's actual target it is not, and the number is
already on the record.** M-316 measured `bonsai` — `u8`, 256³, integer isovalue — at **16,284 of
529,508 surface-cell corners exactly on the isosurface, 3%**, and M-232 measured singular faces at 20
per 400,000 cells at `u8` density against **0** in continuous data. M-317's subgrid extractor declines
483 tetrahedra around 33 of those points. The `=` label is a one-enum-variant change aimed squarely at
the input class M-006 exists to test.

> **H.** On `fuel` (64³) and `bonsai` (256³) at an **integer** isovalue: **(C1)** `degenerate_triangles`
> under `marching_cubes` is non-zero and at least **80%** of it is attributable to cells with a corner
> exactly equal to the isovalue — measured by tagging, not inferred; **(C2)** with the `=` label,
> `degenerate_triangles` falls by at least **10×** while `euler_characteristic`, `non_manifold_edges`
> and `boundary_edges` are **unchanged**; **(C3)** at a **half-offset** isovalue, where an integer
> sample cannot equal the isosurface, the two paths produce **byte-identical** meshes.

**Falsified by** C1 under 80%, which would mean the degenerate count is dominated by ordinary
near-tangency slivers — the thing CLAUDE.md correctly refuses to gate on — and the paper solves a
problem this crate does not have. **Or** by C2 changing χ, which would mean the label is not
topology-preserving as implemented and rule 5 applies.

**C3 is the control and it is not decoration.** M-317's own guidance is to contour at `127.5` rather
than `127` precisely because a half-integer is unattainable by integer data. If the two paths differ
there, the label is doing something beyond the exact-equality case and the whole result is suspect.

**Records** `volume`, `isovalue`, `label_rule`, `degenerate_triangles`,
`degenerate_from_equal_corners`, `equal_corners`, `surface_cell_corners`, `euler_characteristic`,
`non_manifold_edges`, `boundary_edges`, `mesh_hash`, `half_offset_identical`.

---

### P-51 — the rejection bound, tightened by affine arithmetic rather than by a better constant

**The measured baseline this improves on.** M-248 measured empty-cell rejection by Hart's bound at
**16.8%** of cells on `gyroid` against **80.6–95.1%** on every other field, and M-306 identified why:
`gyroid` declares `Lipschitz { l = 2√3 ≈ 3.464 }`, derived correctly at M-244, while M-267 measured its
actual gradient supremum converging to **1.731**. The bound is sound and loose by **2×**, and the
rejection radius scales with it, so the one field where rejection pays least is the field whose constant
is loosest. That is not a coincidence and it is not fixable by declaring a smaller constant — M-267 also
measured that a sampled maximum is not a stable lower bound on a supremum, which is why the declaration
stands.

**What affine arithmetic changes.** Sharp & Jacobson's *Spelunking the Deep* (arXiv:2202.02444, TOG
2022, `10.1145/3528223.3530155`) evaluates a conservative range over a **box** by carrying noise
symbols through the expression graph. Unlike a Lipschitz ball it is *correlation-aware*: `sin(a)·cos(b)`
over a box gets a tighter interval than `|∇| · r` can give, because the terms cannot all be extremal at
once — which is exactly the slack M-267 measured in `2√3`. **[acquire]** F-006/M-249 already established
the shape of the answer for the *directional* variant: exactly nothing on five fields and **1.80× on
the gyroid**, the one field with a loose global bound. This asks whether the box variant does better.

> **H.** With a reduced-affine range over each cell's box, against `cell_is_provably_empty`'s Lipschitz
> test at `subgrid/extract.rs:484`: **(C1)** the rejected-cell **count** on `gyroid` at 17³ rises from
> 688 of 4,096 to at least **1,400**; **(C2)** on `sphere`, `torus`, `box_exact`, `csg_difference` and
> `thin_plate` — all `Lipschitz { l = 1 }`, where a distance field's ball bound is already tight — the
> rejected count rises by **less than 5%**; **(C3)** the mesh is **byte-identical** on every field, which
> is the only property a rejection test must have, since a wrong rejection produces a hole and a hole is
> invisible to every validity gate this crate has.

**Falsified by** C1 under 1,400 — the correlation slack is not where the looseness lives, and the
2× gap is genuinely attainable by the gradient rather than an artefact of the ball. **C2 failing
upward** would be more interesting than C1 holding: a tighter bound on a field that already declares
`l = 1` would mean the *ball* geometry, not the constant, is what costs, and that generalises to every
field.

**Counted, not timed** — ✗24, applied before the fact. Rejected cells are integers and identical on
every machine; the evaluation cost is printed beside them and gates nothing. And the cost side must be
reported honestly: affine arithmetic is `O(k)` per operation in noise symbols and wants a symbolic
field, which the `Real`-generic closure API may not expose without `alloc` in a hot path — against
rule 6. **If C1 holds, the follow-up is an engineering decision about the field API, not another
experiment.**

**Records** `field`, `samples_per_axis`, `bound`, `rejected_cells`, `cells`, `rejected_fraction`,
`mesh_identical`, `mesh_hash`, `bound_evals`, `bound_ns_per_cell`, `extract_ns`.

---

### P-52 — a monotone-edge certificate, for the two fields where χ cannot be asserted

**The gap in the suite, stated precisely.** `validate` checks manifoldness, orientation, Euler
characteristic, self-intersection, isotopy (T-015) and Hausdorff accuracy. **Nothing checks that the
mesh's critical-point structure matches the field's.** On `gyroid` and `fbm_terrain` the crate cannot
even assert χ — CLAUDE.md says so and gives the reason — so those two fields have *no* topological gate
at all beyond manifoldness.

*Topology-Preserving Meshing of Implicit Scalar Fields via Monotonicity Constraints*
(arXiv:2608.12142, Aug 2026) supplies one: if every mesh edge is **monotone** with respect to the
field, the piecewise-linear approximation is consistent with the field's critical points. **[acquire]**
Note what the abstract actually claims — *"correct with regards to critical points **in our
experiments**"* — an empirical claim, not a theorem, and the row should say so.

> **H.** Sampling `k = 16` points along each mesh edge: **(C1)** `marching_cubes` on `sphere`, `torus`
> and `box_exact` at 65³ has **zero** non-monotone edges — the surface is well resolved and an edge
> joining two crossings on adjacent grid edges cannot turn around. **(C2)** `gyroid` and `fbm_terrain`
> have a **non-zero** count that **falls** with resolution across 17³/33³/65³/129³, which makes it a
> resolution witness rather than a defect. **(C3)** `noise_cavity` has the highest count of the eight,
> because it is the field with interior ambiguity.

**Falsified by** C1 non-zero, which would mean the gate is measuring the sampling of the edge rather
than the mesh and `k` is the problem; **or** C2 flat in resolution, which would make it a property of
the field rather than of the grid and therefore useless as a witness.

**The tolerance is a design decision and it is raised rather than assumed**, per CLAUDE.md. A strict
monotone test on `k` samples will fire on float noise wherever the field is flat along an edge. The
registration proposes **"strictly monotone after discarding steps under `1e-12 · (|f(a)| + |f(b)|)`"**
and records the count at three tolerances so the answer's sensitivity is visible — but if you would
rather fix the tolerance differently, that changes the row and should be settled before the harness
runs.

**Records** `field`, `extractor`, `samples_per_axis`, `edges`, `non_monotone_edges`,
`non_monotone_per_1k`, `k`, `tolerance`, `worst_reversal`.

---

### P-53 — are the 216 golden hashes pinned to a definition, or to a version of `libm`?

**This is an epistemic experiment, not a performance one, and it questions a claim the crate rests
on.** CLAUDE.md justifies `libm` unconditionally on determinism: *"`libm` is pure Rust and
bit-reproducible everywhere"*, and M-31 verified it — 216 hashes generated on macOS/arm64, passing on
Linux/x86-64, and now on a third container. That is **platform** independence. It is not
**implementation** independence.

`libm`'s `sin`/`cos` are not correctly rounded. CORE-MATH ships correctly-rounded elementary functions
and `rust-lang/libm#248` tracks adopting them. **If `libm::sinf` and a correctly-rounded `sinf` differ
on any input `gyroid` reaches, then the golden hashes are pinned to a `libm` version rather than to a
mathematical definition, and a future patch bump can invalidate them silently.** That is a real
maintenance hazard for a crate whose central determinism claim is a committed hash set.

> **H.** **(C1)** over the `gyroid` sample set at 33³ (35,937 points, `sin` and `cos` each), `libm`'s
> result and a correctly-rounded reference differ on **at least one** input; **(C2)** at least one
> such difference propagates to a **changed golden hash** among the 24 `gyroid` rows; **(C3)** the
> other seven fields, which use only `sqrt`, `abs`, `min` and `max` — all correctly rounded in IEEE by
> definition — are **unchanged**, which is the control that says the effect is transcendental and not
> the harness.

**Falsified by** C1 finding zero differences, which would be a genuinely reassuring null and worth an
entry: it would mean `libm` *is* correctly rounded on this crate's inputs and the hashes are pinned to
the definition after all. C2 failing while C1 holds is the middle case and is also worth recording —
the difference exists and does not reach a mesh, which bounds the exposure.

**The reference must not become a dependency.** `core-math` the crate is a C binding: not `no_std`, a
hard no under rule 3. The comparison is a **dev-only** oracle, run once, and the deliverable is a
finding about a risk — not a proposal to change a dependency. If C1 and C2 both hold, the honest
consequence is one sentence in `CLAUDE.md` next to the `libm` justification: *the hashes are pinned to
`libm 0.2.x`, and a bump is a re-baseline.*

**Records** `function`, `inputs_tested`, `differing_inputs`, `worst_ulp`, `field`,
`golden_rows_checked`, `golden_rows_changed`, `control_fields_unchanged`.

---

### P-54 — the bitmap prepass, gated on a count instead of a stopwatch

**Why this exists.** §A.2. M-337's C2 is a wall-clock ratio, registered as a threshold, and the
committed `p-40.csv` does not clear it (1.1925× and 1.2160× on Surface Nets against ≥ 1.25×) while the
entry records it as held on 1.336×. ✗24 already earned the rule and cost a release doing it: *gate the
count the ratio samples.*

**What the mechanism actually changes is countable, and it is already in the CSV.** The bitmap prepass
does not make a gather faster; it **removes** gathers. The deterministic quantity is *eight-corner
gathers performed*, which under the scalar path is `cells` and under the bitmap path is `active_cells`,
plus the per-word bitmap build. `p-40.csv` already carries both columns: `sphere` at 128³ is
19,010 active of 2,048,383 — **0.93%**.

> **H.** With `gathers_performed` instrumented on both paths: **(C1)** the scalar path performs exactly
> `cells` gathers and the bitmap path exactly `active_cells`, on every row, **as an equality** — not a
> ratio and not a tolerance; **(C2)** the bitmap build performs exactly `sample_count` comparisons and
> `ceil(size_x/64) · size_y · size_z · c` word operations for a fixed small `c` determined by reading
> `active_word`, so the whole prepass is `O(n³/64)` word ops against the `O(n³)` gathers it replaces;
> **(C3)** the mesh hash is unchanged on all 12 rows, which `p-40.csv` already reports as 12 of 12 and
> which C1 and C2 make a *consequence* rather than a hope: the set-bit walk visits cells in ascending
> `x`, so vertex creation order is the scalar order.

**Falsified by** C1 not being an exact equality — a gather performed on an inactive cell means the mask
is wrong, and a mask that is wrong in the *other* direction produces a hole, so this is a correctness
gate wearing a performance gate's clothes. **Timing stays in the CSV and gates nothing**, reported with
`ghz` beside it per M-280.

**And the sign-bit trap stays named**, because it is the shortcut a later reader will reach for: the
bit is `v < 0`, **not** the IEEE sign bit. `-0.0` has the sign bit set and `-0.0 < 0.0` is false, and
`box_exact` is exactly zero across its entire boundary — so a sign-bit build would be faster, would
pass every timing clause, and would change a reference field's mesh.

**Records** `field`, `samples_per_axis`, `extractor`, `cells`, `active_cells`, `gathers_scalar`,
`gathers_bitmap`, `bitmap_word_ops`, `bitmap_comparisons`, `gather_ratio`, `mesh_hash`,
`mesh_identical`, `ns_per_sample`, `ghz`.

---

### P-55 — the one vertex in 57,470, and whether it is a seam or a coincidence

**P-47's residue, turned into a claim that can fail.** P-47's registered accuracy clause died by three
orders of magnitude — mean 7.6e-5° against a predicted 0.1° — and the artefact says why: the
**bulk** mean is **1.9e-8°** and one vertex carries the rest, at 4.365°. `vertices_over_1deg` = 1,
`vertices_over_5deg` = 0, and `worst_stencil_straddles_seam` is `true` at 32 and 64 brushes and `false`
below. So the hypothesis that survives is far narrower and far more useful than the one registered:
central differences on a composed field are effectively exact **except** where the six-sample stencil
crosses a CSG seam, and there the error is bounded by the seam's dihedral rather than by `h`.

**The mechanism, derived rather than guessed.** At a `min`/`max` seam the field is `C⁰` and not `C¹`.
A central difference straddling it averages two different gradients, so the returned direction lies in
the cone the two branches span and the error is at most **half the angle between them** — which is
M-283's `(180° − θ)/2` in a second setting, and which does **not** shrink with `h` because the stencil
step is `DIFF_STEP · |p|`, independent of the grid. The measured 4.365° therefore predicts a seam
dihedral of about **171°**, a very shallow crease, which is consistent with the error being rare.

> **H.** Over a swept family of two-sphere `Subtract` fixtures with the seam dihedral `θ` controlled
> from 30° to 175°: **(C1)** every vertex whose central-difference stencil straddles the seam has
> angular error against the dual-number normal **bounded by `(180° − θ)/2`**, on every fixture, with
> no exceptions; **(C2)** the count of such vertices scales like the seam's **length in cells**, so it
> is `O(n)` on an `n³` grid rather than `O(n²)`, which is why one vertex in 57,470 is the expected
> order and not a fluke; **(C3)** vertices whose stencil does **not** straddle the seam have mean error
> under **1e-6°**, which is the control that says the effect is the seam and not the tape.

**Falsified by** C1 exceeding the bound, which would mean the error is not the two-branch average and
M-283's mechanism does not transfer; **or** C2 scaling as `n²`, which would make it a surface-wide
effect and change what a consumer should do about it.

**What holds from P-47 and should be written up as it stands:** C2's **2.84×** speedup and C3's
`2.83e-16` control, and the fact that `dual_value_bit_exact` is `true` — the dual number reproduces the
field's *value* bit-for-bit, which is the strongest available evidence that the derivative arm is
wired to the same expression. Note also the two rows where the dual is **slower** — `capsule` 0.877×
and `sphere` 0.863× — because a field with an analytic override never pays the six samples. **The
speedup is a property of composed tapes, not of dual numbers**, and the row should say so.

**Records** `dihedral_deg`, `samples_per_axis`, `seam_cells`, `straddling_vertices`,
`straddling_max_error_deg`, `predicted_bound_deg`, `within_bound`, `non_straddling_mean_error_deg`,
`vertices`, `scaling_exponent`.

---

## Part D — Owed re-registrations, and what I did not pick

### D.1 Two re-registrations the artefacts already earned

**P-45a** — the Gaussian clause with **Borel** accounting (each vertex to exactly one chunk) rather
than per-incidence, since `borel_gaussian_gap` is already 1e-12 on all three fields and
`gaussian_gap_over_pi` is exactly `excess_chunk_incidence`. And clause three split: the **Gaussian**
measure is chunk-local bit-for-bit (0 mismatched chunks); the **mean** measure needs a one-ring halo
and the honest prediction is that it is composable and not chunk-local. `box_exact` passing C3 because
16,548 of its 17,292 edges are flat is the control that proves the instrument works.

**P-46a** — the repair, re-asked with the cascade as the *subject* rather than as a clause. The
registered ≤ 2 sweeps was falsified at 5 on a field that converged and 64 on one that did not, with
6,307 stuck cells. The question worth asking: **is the stuck set exactly the interior-ambiguity cells?**
M-208 says `noise_cavity` is the only reference field that reaches six body saddles; M-224/M-225 say
Manifold Dual Contouring's residue on that field is predicted exactly from the grid. If the 118
survivors sit in cells with an interior ambiguity, then well-composedness of the *sign lattice* is the
wrong object and the ambiguity is in the *interpolant*, which is a different repair. And read
`10.1007/978-3-319-18720-4_47` first — §B.1 item 3 — because the method P-46 cites is not the method it
implemented.

### D.2 Considered and rejected, so nobody searches twice

| Candidate | Why not |
|---|---|
| **Decoupled-fallback GPU scan** (`10.1145/3694906.3743326`, SPAA 2025) — 30–50% over reduce-then-scan, 98–113% of memcpy, **measured on Apple M3 and M1 Max**, and Apple silicon is named as lacking forward-progress guarantees | Genuinely good and platform-matched. Rejected **for now** only because M-167 measured that after GPU-010a/011a the scan is no longer where the time goes — synchronisation was 83% and the arithmetic 3.7%. Worth reopening if the GPU path is revived; the FPG hazard is a latent correctness issue on your target platform independent of speed. |
| **SpUDD / unsigned power diagrams** (arXiv:2604.19568) — has a real convergence **theorem** | Needs exact predicates and a regular triangulation: breaks rule 3, not `no_std`, not per-cell parallel. Its capability (open, non-orientable surfaces) also conflicts with your manifoldness gates. |
| **Adaptive tetrahedral grids for implicit complexes** (`10.1145/3658215`) | Opposite of a chunked uniform grid, offline setting, and its predecessor leans on exact predicates. |
| **GPU work graphs** (`10.1145/3675376`) | DX12/Vulkan only. No wgpu surface, no Metal path. Blocked at the API layer, not the math layer. |
| **Transform-aware sparse voxel DAGs** (`10.1145/3728301`) | DAGs are for static scenes and must be rebuilt on edit — the exact inverse of your access pattern. |
| **Reproducible summation** (`10.1145/3389360`) | A-016 (M-175/M-176/M-177) already settled the order question for the QEF accumulation, and established that negation equivariance is not achievable by ordering at all. |
| **Rust `core::simd`** | Still unstable. A `no_std` core cannot use portable SIMD on stable. |
| **New Marching Cubes case tables since 2023** | None exist. The line ends at Custódio 2019 / Vega JCGT 2019 — which is why P-50 is a 2019 paper and not a 2026 one. |
| **Neural extractors** | Every headline is bought with a learned prior on a shape dataset. The two learning-free exceptions are P-48's and P-49's sources, which is why they are in Part C. |

### D.3 One process item, which is not an experiment

`doc_facts.sh` gates counts. Nothing gates **numbers quoted inside a `FINDINGS.md` entry against the
CSV that entry names** — which is how M-336 and M-337 got in, and how M-287 got in before them. The
check is mechanical for the subset that matters: an entry that names `docs/experiments/p-NN.csv` and
prints a markdown table of numbers can have those numbers matched against that file's columns. It
would have failed on Phase 19 twice. Whether that is worth building is your call; it is the cheapest
defence against the one systematic error this audit found.

---

## The short version

- **Seven load-bearing rows reproduce exactly** against their committed artefacts, including ✗1's
  identity at every resolution, ✗25's whole table, M-338's bijection, and M-340's curvature bound to
  every digit. The suite passes on a third architecture, which re-confirms M-31.
- **Two Phase 19 entries quote runs that are not in the repository.** ✗28's conclusion survives it;
  **M-337's C2 verdict does not** — the committed CSV says 1.19×, under the registered 1.25× bar. Both
  clauses were wall-clock ratios, which ✗24 already said should not be gates.
- **P-45, P-46 and P-47 have run and falsified themselves, and have no ledger entries.** Their CSVs
  carry the diagnoses: an excess-incidence count that is exactly `gap/π`, a repair that stalls on 6,307
  cells of the one field with an interior ambiguity, and a normal-error mean carried by a single
  seam-straddling vertex.
- **Four citation corrections**, one of which — Boutry's tutorial standing in for Boutry's method — is
  worth acting on before P-46 is re-attempted.
- **One gate is red at HEAD** (38 examples, 36 claimed).
- **Eight registrations**, of which P-48 and P-50 are the cheapest and most likely to pay: one pass
  over existing output, and one enum variant aimed at a defect this crate has already measured on real
  CT data.
