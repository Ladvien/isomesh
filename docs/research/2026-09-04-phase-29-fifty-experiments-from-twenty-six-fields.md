# Fifty registrations from twenty-six fields of mathematics, and the instrument that was in the case table all along

**Date:** 2026-09-04 · **Repo state:** `31ac0b8` on `mac-air` (`main` = `origin/main`; **confirm `big` has not
moved before registering**), 600 findings entries, 25 open tickets
**Corpus:** home-still, 9,602 documents, 293,344 chunks, both scribe instances healthy, distill on CUDA
**Companion:** `docs/research/2026-09-04-phase-29-axes-and-vocabulary-v3.md` (the words)

**Numbering.** Phase 27 closed at `P-176`; Phase 28's five tickets (`R-177`–`R-181`) are landings and carry
no `P-`. This phase is drafted as **`P-177`–`P-226`** against tickets **`R-182`–`R-231`**. Phase 27 was
renumbered twice on merge; if `big`'s `experiment.rs` has moved past `P-176`, shift the whole block and
say so in `FINDINGS.md`, as `M-440` did.

**Method.** Twenty-six areas were probed against the corpus and the repository at once, on the rule the
2026-08-23 sweep earned. This time the probes were chosen by reading Phase 27's *results* rather than its
registrations: every group below hangs on a measured number from `M-440`–`M-489` that the field in question
can explain, bound, or replace. Roughly 40 semantic queries and 60 word-boundary greps; **23 papers acquired** (converted or
converting at time of writing), **27 recorded as paywalled or unreachable by the resolver**, and **2
acquired by a wrong DOI and flagged for removal** (§1).

**The calibration, restated.** Present at **0.65–0.72**; absent at **0.54–0.62** with an unrelated top hit.
Two subjects this sweep needed were *already present and never cited*: Niyogi–Smale–Weinberger (0.69) and
Cuel–Lachaud–Thibert's multigrid-convergent Voronoi estimator (in the corpus with a proved theorem). Both
are Group E's and Group C's fault for not having been searched, not the corpus's.

**Protocol.** Phase 15's protocol applies in full. All fifty `P-` entries are registered in
`crates/isomesh/src/experiment.rs` **before** any harness commit; that registration is the first commit of
the phase and is the only write to `crates/isomesh/src/**` that precedes a measurement. `✗51`'s rule is
applied: every clause stated as a ratio of a total carries a SHARE line, and every registration carries a
VACUITY CONTROL naming the column that proves the fixture could have failed.

---

## 0. The headline, and it is not a paper either

Phase 27's headline was an identity already in the crate's source. This phase's is a *sum* already
available to it. The extractor indexes every surface cell by its `2×2×2` sign configuration — the 256-entry
case index. **Ziegel & Kiderlen (2010) prove that a weighted sum over exactly those configuration counts
estimates surface area with an asymptotic worst-case error under 4%, and that no local weighting can do
better** [verified, read after conversion: *"Using this optimal weight vector, (9) will imply that the
relative worst case error is asymptotically less than 4%"*; only **five** configuration types carry weight,
the weight family is `λ̄ = (sr, 1.358, 2.352 − sr, 0.960, s(1 − r))` with `r ∈ [0,1]`, `s ∈ [1.663, 1.745]`;
Lindblad's earlier weights reach **7.3%**; the space spanned by all `2×2×2` h-functions has dimension
**50**; the result is stated for the **Gauss digitization**, which is what `field.sample()`'s sign grid is,
and the authors say it *"only hold[s] for digitizations with good or very good resolution"*]. Legland, Kiêu & Devaux (2007) give the lattice Crofton
version with **13 directions and an explicit directional bound of `−7.33%`/`+2.27%`** [verified, read in
full], and their Table 2 reads **1600 against a true 2400** for an axis-aligned cube with 3 directions.

That matters here for one reason. Phase 27's largest structural finding was that **four of eight reference
fields are ungradeable**: `validate::accuracy` needs `field.bound() == Exact`, and `✗107`'s anisotropy
mechanism demonstrably worked on the one field (`fbm_terrain`) that could not be graded. A mesh-free area,
mean curvature and Euler characteristic from the sign grid alone grade every field, and their error is a
theorem rather than a bound the field must volunteer. Group A lands that instrument and Groups B, J and Q
use it.

Three further things were verified during the sweep and are recorded once here:

- **The lattice-point discrepancy of a sphere has exponent between `1/2` and `21/32`** in `x = R²`:
  `P₃(x) = O(x^{21/32+ε})` (Heath-Brown), `Ω±(x^{1/2}(log x)^{1/2})` (Szegő; Tsang), and Jarník's
  mean-square `∫₀ˣ P₃² = c₃X² log X` [verified, Ivić et al. §2]. In radius: RMS jitter `~ R (log R)^{1/2}`
  against a surface-cell count `~ R²`. That is a *prediction* for how much a translating sphere's triangle
  count flickers, and Group M registers it.
- **The expected number of critical points of an isotropic Gaussian field depends on the covariance only
  through `ρ'(0)` and `ρ''(0)`** [verified, Cheng & Schwartzman §1, §3]. Two numbers predict the ambiguous
  cell census. Group B's fixture is what makes that testable.
- **Weyl–Steiner: the volume of the `r`-offset of a set is a polynomial in `r` whose coefficients are the
  intrinsic volumes, for `r` below the reach** [verified, Taylor 2006 §1, eq. 1.1]. For an SDF the
  `r`-offset is `{f ≤ r}` for free, so the reach is the radius at which a cubic fit breaks. Group E.

---

## 1. What the twenty-six probes found, in one table

| Area | Corpus | Repo cites | Verdict |
|---|---|---|---|
| Integral geometry / stereology (intrinsic volumes, Crofton, configuration weights) | ABSENT (0.60) → **acquired** Legland, Schladitz, Ziegel–Kiderlen, Lindblad ×2, Schmalzing | 0 | **The headline. Group A** |
| Random-field geometry (GKF, Kac–Rice, EC densities) | ABSENT (0.60) → **acquired** Taylor 2006, Cheng–Schwartzman, Letendre | 0 | Group B |
| Digital geometry (digital planes, multigrid convergence) | ABSENT (0.60) for planes; **PRESENT** Cuel et al. for multigrid convergence | 0 | Group C |
| Higher-order output (PN triangles, superconvergence, Richardson) | marginal (0.65, a textbook); PN/Wang–Olano/ZZ **paywalled** | 0 | Group D — the lever `✗116` named |
| Reach, tubes, ε-samples | **PRESENT 0.69** NSW; Berenfeld already held; **acquired** Aamari | **0** | Group E — paid for, unmined |
| Besov / Hurst / tree approximation | ABSENT (0.61); DeVore, CDDD **paywalled** | `Hurst` in 4 files, no exponent measured | Group F |
| Entropy rate / source coding | ABSENT (0.57); Ziv–Lempel **paywalled** | 0 | Group G — Phase 25's denominator |
| Kinetic data structures | **PRESENT 0.64** (Acar ×2, Agarwal) | 3 mentions, no certificate | Group H |
| Optimal transport | ABSENT for surfaces (0.66 wrong docs) → **acquired** Peyré–Cuturi, Solomon | 2 mentions | Group I |
| Robust statistics; conformal prediction | ABSENT (0.66 / 0.60) → **acquired** Angelopoulos–Bates; Huber **paywalled** | `Huber` 1 file | Group J |
| Quantum information (3-tangle, SLOCC) | ABSENT → **acquired** CKW, Dür–Vidal–Cirac | W-state named in `P-131` | Group K |
| Discrete Morse on cubical grids; QTAIM/Bader | ABSENT (0.60) → **acquired** Robins–Wood–Sheppard, Yu–Trinkle | ranked #7 on 2026-08-23, never done | Group L |
| Analytic number theory (lattice points) | ABSENT (0.59) → **acquired** Ivić et al. | 0 | Group M |
| Crystallographic symmetry reduction | ABSENT (0.58) | `space group` named in `P-143` | Group N |
| Shape / topological derivative | ABSENT (0.65, topology-optimisation paper) → **acquired** Sokolowski–Zochowski | 0 | Group O |
| Stochastic sampling / dither | ABSENT (0.62); Cook, Wannamaker **paywalled** | `Poisson-disk` 1 file | Group P |
| Design of experiments | not probed in corpus; McKay–Beckman–Conover **paywalled** | 0 | Group Q |
| Backward error / conditioning | ABSENT (0.62) | 0 | Group R |
| Affine arithmetic | **PRESENT 0.68** Fryazinov et al. | 8 files | Group S — a different use |
| Digital topology / Alexander duality | ABSENT (0.58); Kong–Rosenfeld **paywalled** | well-composedness mined | Group T |
| Currents / flat norm | ABSENT → **acquired** Morgan–Vixie | 0 | Group U |
| Approximation constants (Blu–Unser) | held since Phase 27 | 0 for the constant | Group V |
| Chen–Stein / Poisson approximation | ABSENT; Arratia et al. **paywalled** | 0 | Group W |
| Coarea / divergence theorem | ABSENT (0.55) | 0 | Group X |
| Quasi-Monte Carlo / discrepancy | ABSENT (0.60); Dick–Kuo–Sloan **paywalled** | 0 | Group Y |
| Harmonic analysis of the error field | not a corpus subject | 0 | Group Z |
| Isogeometric analysis | ABSENT → **acquired** Hughes et al. (HTML) | 0 | vocabulary only |
| Extreme-value theory for frame times | ABSENT (0.65, wrong field) | 0 | **killed** — not a meshing question |
| Gaussian-process implicit surfaces | **PRESENT** (sig2026) | 0 | deferred — Pöthkow paywalled, and Group J's conformal row covers the use |

**Killed during the probe, so nobody re-searches them.** *Extreme-value theory for the worst frame* — a
statistics-of-the-harness question, not a meshing one. *Sum-of-squares certificates for empty cells* — the
2026-08-26 audit already showed the box's Positivstellensatz is Handelman's LP and the Bernstein form is
that certificate; nothing new. *Spectral graph theory of the cave graph* — gameplay, not extraction.
*Well-composed repair before extraction* — mined on 2026-08-23. *Cellular automata, sheaves, category theory,
Ramsey theory* — probed for a hook into a measured number and none was found.

**Corrections owed to earlier documents.** `2026-08-29-phase-27-fifty-experiments...md:117-124` lists nine
paywalled DOIs and the prose says six; Phase 27's `M-440` already recorded that. This document's own count is
twenty-seven and is in §1 below. `2026-08-12-axes-and-vocabulary.md` Axis 10 says chunk consistency has no
theory; Group N's mirror-plane seams are a theorem for the symmetric case.

**Two DOI traps, both this session's own.** `10.48550/arXiv.1504.01720` was downloaded as Letendre 2016
and is *Isoperimetric regions in ℝⁿ with density rᵖ*; Letendre is `arXiv:1408.2107` (acquired).
`10.1016/j.cviu.2014.04.005` was downloaded as Coeurjolly–Lachaud–Levallois and is an action-recognition
paper; the correct DOI is `10.1016/j.cviu.2014.04.013` (acquired, HTML). **Both wrong papers are now in the
corpus and should be removed with the `hs` CLI** — the MCP has no delete.

**Blocked acquisitions — paywalled or unreachable by the resolver this session.** Lindblad 2005
`10.1016/j.imavis.2004.06.012` (the 2002 and 2003 conference versions were acquired instead); Vlachos et al.
*Curved PN triangles* `10.1145/364338.364387`; Wang & Olano `10.1145/2699276.2721400`; Zienkiewicz–Zhu
`10.1002/nme.1620330702`; DeVore 1998 `10.1017/S0962492900002816`; Cohen–Dahmen–Daubechies–DeVore
`10.1006/acha.2000.0336`; Basch–Guibas–Hershberger `10.1006/jagm.1998.0988`; Huber 1964
`10.1214/aoms/1177703732`; Henkelman et al. `10.1016/j.commatsci.2005.04.010`; Arratia–Goldstein–Gordon
`10.1214/aop/1176991491` and `10.1214/ss/1177012015`; Christiansen–Petersen `10.1007/BF01932705`; Roache
`10.1115/1.2910291`; Ziv–Lempel `10.1109/TIT.1978.1055934`; Dick–Kuo–Sloan `10.1017/S0962492913000044`;
Boubekeur–Alexa `10.1145/1409060.1409094`; Cook `10.1145/7529.8927`; Coeurjolly et al. 2003
`10.1007/3-540-36586-9_7`; Edelsbrunner et al. time-varying Reeb `10.1016/j.comgeo.2007.11.001`; Raman &
Wenger `10.1111/j.1467-8659.2008.01209.x`; Kong & Rosenfeld `10.1016/0734-189X(89)90147-3`; Worsley 1995
`10.2307/1427930`; Taylor & Adler 2003 `10.1214/aop/1046294306`; Adler 2000 `10.1214/aoap/1019737664`;
Wannamaker et al. `10.1109/78.824859`; McKay–Beckman–Conover `10.1080/00401706.1979.10489755`. **Free
copies located but not fetchable by DOI:** Wiley et al. 2003 (escholarship `2fp9j33g`); Edelsbrunner &
Mücke is already in the corpus as `arXiv:math/9410209`.

**Sixteen rows are expected to return nulls, and that is registered rather than hoped:** `P-184`, `P-191`,
`P-197`, `P-200`, `P-203`, `P-204`, `P-207`, `P-209`, `P-211`, `P-213`, `P-216`, `P-219`, `P-220`, `P-221`,
`P-223`, `P-225`. Two more (`P-198`, `P-222`) are denominators whose most valuable outcome is that the
numerator is already close.

**Dependencies, stated once.** `P-177` runs before `P-180` and `P-216`. `P-181` runs before `P-182`,
`P-183`, `P-184` and `P-223`. `P-189` runs before `P-188`, `P-190` and `P-191`. `P-193` runs before `P-192`
and `P-194`. `P-198` runs before `P-199` and `P-200`. `P-212` runs before `P-213`. `P-180` runs before
`P-226`. Everything else is independent. Sizes: **23 `S`**, **22 `M`**, **5 `L`** (`P-181`, `P-188`,
`P-201`, `P-212`, `P-217`).

---

# Group A — mesh-free measurement (integral geometry, stereology)

#### P-177 — registered for R-182, before the harness: The case table already carries a ≤4% area estimator, and no field is ungradeable

**Ticket:** `R-182` (M). **Records:** `field`, `resolution`, `bound_kind`, `config_count_nonzero`, `area_mesh`, `area_zk_weights`, `area_crofton_13`, `area_crofton_3`, `area_true_or_none`, `rel_err_mesh`, `rel_err_zk`, `rel_err_crofton_13`, `zk_within_4pct`, `estimator_cost_share`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** THE EXTRACTOR INDEXES EVERY SURFACE CELL BY ITS `2×2×2` SIGN CONFIGURATION, AND A WEIGHTED SUM OVER THOSE COUNTS IS A PROVEN SURFACE-AREA ESTIMATOR. Ziegel & Kiderlen determine `2×2×2` configuration weights whose *"relative worst case error is asymptotically less than 4%"* and show that bound is best possible among local weightings [verified]; their weight family is `(sr, 1.358, 2.352 − sr, 0.960, s(1 − r))` over the five informative configuration types, and Lindblad's weights reach 7.3%; Legland et al. give the 13-direction lattice Crofton estimator with directional error in `[−7.33%, +2.27%]` [verified]. (C1) On the four `Exact` fields at 33³, 65³ and 129³, the Ziegel–Kiderlen sum is within **4%** of the true area on every row and the 13-direction Crofton sum within **7.33%**, both computed from the sign grid with no mesh. (C2) On the same rows the *mesh* area is within **2%** of the true area, so the mesh is the better instrument where a truth exists — and where none exists (`csg_difference`, `gyroid`, `fbm_terrain`, `noise_cavity`) the configuration estimator and the mesh area agree to within the estimator's own bound on **every** row, which makes the estimator a grade for the four fields `✗107` could not grade. (C3) The estimator costs under **2%** of an extraction, because the counts are the case-index histogram the extractor already walks. SHARE: C3 is the only clause with a cost; C1 and C2 are accuracy clauses and move nothing.

**Falsified by.** C1 by any `Exact` row outside the proven bound, which would mean either the weights were transcribed wrong, the fixture is not a Gauss digitisation in the paper's sense (`field.sample()` is one, and the row must say so), or the resolution is below the *"good or very good"* regime the authors scope the bound to — in which case the row reports the rung at which the bound starts to hold, which is itself a number the paper does not give. C2 by the mesh and the estimator disagreeing beyond the bound on an ungradeable field, which is the interesting outcome: one of them is wrong and only a third instrument (`P-224`) can say which. C3 by a cost above 2%, which would say the histogram is not free after all. VACUITY CONTROL: `config_count_nonzero` must show at least 30 distinct configurations populated on every field, or the sum is over a table the fixture never exercises.

#### P-178 — registered for R-183, before the harness: The 6/26 split is an adjacency-pair choice, and the theorem says which pairs are consistent

**Ticket:** `R-183` (S). **Records:** `field`, `resolution`, `adjacency_pair`, `chi_lattice`, `chi_p145_oracle`, `chi_mesh`, `disagreement_cells`, `split_blocks`, `complementary_pair_agrees`, `c1_holds`, `c2_holds`.

**Hypothesis.** `M-458`'s oracle reports `connectivity_split_blocks` **487** on `noise_cavity` at 33³ and *"must be a reported column rather than a hidden choice between foreground-6 and foreground-26"*. Ohser, Nagel & Schladitz's theory of the Euler number on lattices says the choice is not free: foreground and background adjacencies must be chosen as a **complementary pair** — `(6,26)`, `(26,6)`, `(18,6)`, `(6,18)`, and the `14.1`/`14.2` pairs — for χ to satisfy the Euler formula on the reconstruction [verified, Legland et al. §"Euler–Poincaré measures", citing Ohser et al. 2002]. (C1) Computing χ by the lattice reconstruction with each complementary pair reproduces `P-145`'s oracle exactly on every field where `P-145` reports `split_blocks = 0`, and on `noise_cavity` the **487** split blocks are exactly the cells where the two members of a pair *disagree with each other* — i.e. the split count is the number of cells whose χ contribution depends on which of the pair is foreground. (C2) The mesh χ (Marching Cubes) sits on one specific pair — the registration predicts `(6,26)`-foreground-inside — and the disagreement between mesh and oracle on `noise_cavity` (two extractors, two numbers, `M-458`) is removed by choosing that pair. SHARE: none; this is a classification of an existing count.

**Falsified by.** C1 by the split count not coinciding with the pair-disagreement count, which would mean `P-145`'s split is measuring something the adjacency theory does not name. C2 by no pair reproducing the mesh χ, which would say Marching Cubes' implied adjacency is not a lattice adjacency at all — a real finding about the case table. VACUITY CONTROL: the non-complementary pair `(26,26)` must be run as a control and must **fail** the Euler formula on at least one field, or the harness cannot tell a consistent pair from an inconsistent one.

#### P-179 — registered for R-184, before the harness: Integrated mean curvature from the sign grid, as a third instrument against normal cycles

**Ticket:** `R-184` (M). **Records:** `field`, `resolution`, `mean_breadth_crofton_13`, `mean_breadth_true_or_none`, `integrated_H_normal_cycles`, `integrated_H_varifold`, `rel_err_crofton`, `rel_err_normal_cycles`, `direction_count`, `c1_holds`, `c2_holds`.

**Hypothesis.** The mean breadth `b̄ = ∫χ(P) dP` over planes is the third intrinsic volume and Legland et al. compute it from the sign grid with 13 plane directions [verified]; on a sphere of radius 30 it reads **39.9 against 40.0**, on a torus **46.2 against 47.1** [verified, Table 3]. `P-174` benched normal cycles against varifolds for mean curvature and neither had a mesh-free reference. (C1) On `sphere` and `torus` the lattice mean breadth is within **3%** of the closed form at 65³ and 129³. (C2) On every field the lattice value and the normal-cycles integral of `H` over the mesh agree to within **5%**, so the normal-cycles instrument is graded from the field for the first time. SHARE: none — a validation instrument.

**Falsified by.** C1 by errors above 3% on the two closed forms, which would mean the plane-direction weights were transcribed wrong. C2 by disagreement on `box_exact` in particular — where Legland's own Table 3 shows the cube is the worst case (**27.2 against 30.0**) — in which case the row records the *field-dependent* error and the instrument is scoped to smooth fields. VACUITY CONTROL: `direction_count` must be run at both 3 and 13, and the 3-direction cube error must reproduce the published **−33%**, or the harness is not computing the estimator the paper describes.

#### P-180 — registered for R-185, before the harness: The error has a Wulff plot, and the plot is where BCC and cubic actually differ

**Ticket:** `R-185` (S). **Records:** `normal_direction`, `lattice`, `resolution`, `area_mesh_over_true`, `hausdorff_over_h`, `needle_ratio`, `octahedral_class`, `class_spread`, `anisotropy_index`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `M-309` measured a `0.250` mean triangle-quality ratio on `box_exact` against Marching Cubes' `0.858`; `P-162`/`P-163` measured BCC's *mean* gain against cubic in dB and found it small. Both are averages over a quantity that is a **function on the sphere of normals**: Legland's Table 2 shows the same cube reading **1600, 1989.3, 2489.3** for a true 2400 at three orientations [verified]. (C1) Sampling a plane's normal over a 48-fold-reduced spherical grid (`✗49` licenses the reduction), mesh area error, Hausdorff error and needle ratio are each a smooth function of direction with octahedral symmetry, and the area-error plot's extremes bracket Legland's directional bound `[−7.33%, +2.27%]` for the 13-direction Crofton arm. (C2) The **same plot on BCC** has a different shape: its extremes are closer to zero (the registration predicts the max-min spread falls by at least **30%**) even though its mean is within 0.3 dB of cubic, which is what `✗119`'s **+5.73 dB on `noise_cavity`** against a predicted 0.011 dB was without a direction axis. (C3) The needle-ratio plot predicts `M-309`'s `0.250` as its value at the axis-aligned direction. SHARE: none; an instrument.

**Falsified by.** C1 by the plot lacking octahedral symmetry, which would contradict `✗49` and be a bug. C2 by BCC's spread not falling, which would say the lattices differ in mean only and `P-162`'s null stands as a full description. C3 by the axis-aligned value not reproducing `M-309`. VACUITY CONTROL: `class_spread` must be non-zero within at least one octahedral class across resolutions, or the plot is a constant and the symmetry test is vacuous.

---

# Group B — the field as a random object (random-field geometry, probability)

#### P-181 — registered for R-186, before the harness: A Gaussian random field fixture, and the closed-form χ that comes with it

**Ticket:** `R-186` (L). **Records:** `seed`, `resolution`, `correlation_length`, `rho1`, `rho2`, `spectral_moment_check`, `chi_measured_oracle`, `chi_measured_mesh`, `chi_expected_gkf`, `chi_expected_boundary_corrected`, `seeds`, `chi_mean`, `chi_sd`, `z_score`, `gaussianity_test_p`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `P-176` established `fbm_terrain` is hash-based lattice noise and not a Gaussian field; `P-144` found χ varies across seeds by **10.24** with no prediction to compare against. Taylor's Gaussian kinematic formula gives `E[χ(M ∩ f⁻¹[u,∞))]` as `Σ_j L_j(M) ρ_j(u)` — Lipschitz–Killing curvatures of the domain times EC densities of the field — with the domain's boundary handled by the `L_j(M)` terms [verified, Taylor 2006 §1 and Thm 4.1]. (C1) A new bench-local field `gaussian_rf` — a sum of `K ≥ 64` random-phase cosines with an isotropic spectral density, unit variance, and known `ρ'(0)`, `ρ''(0)` — passes a multivariate normality test at every resolution (`gaussianity_test_p > 0.01`), and its measured spectral moments agree with the construction to **1%**. (C2) Over **32** seeds at 65³, the mean of `P-145`'s oracle χ agrees with the GKF prediction within **2 standard errors**, at `u = 0`. (C3) Over the same seeds the two value-noise fields' χ, rescaled to the same correlation length, sit **more than 3 standard errors** from the GKF prediction — the quantitative version of `P-176`'s *"not a Gaussian field"*. SHARE: none — a fixture and three predictions.

**Falsified by.** C1 by a failed normality test, which would mean `K` is too small or the phases are not independent, and the fixture is not the object the theorem is about. C2 by a `z_score` above 2 on the GRF itself, in which case either the EC densities were transcribed wrong or the boundary terms were dropped — the row must compute both `chi_expected_gkf` (interior only) and `chi_expected_boundary_corrected` so the harness can tell which. C3 by the noise fields sitting *inside* the GRF band, which would be a strictly better result: it would say the closed form applies to them after all. VACUITY CONTROL: `chi_sd` must be non-zero across seeds and the oracle's `split_blocks` must be reported per seed, or C2 is a comparison against one number.

#### P-182 — registered for R-187, before the harness: Kac–Rice predicts the ambiguous-cell census from two numbers, before any mesh exists

**Ticket:** `R-187` (M). **Records:** `seed`, `resolution`, `cells_per_correlation_length`, `rho1`, `rho2`, `ambiguous_cells_measured`, `body_saddle_cells_measured`, `critical_points_expected_kr`, `ratio_measured_over_expected`, `ratio_stability_across_h`, `c1_holds`, `c2_holds`.

**Hypothesis.** Cheng & Schwartzman evaluate the Kac–Rice formula for isotropic Gaussian fields and *"the obtained formulae depend on the covariance function only through its first and second derivatives at zero"* [verified]; for `N=2` the expected count per unit area is `1/(√(3π) η²)`, `η = √(−ρ')/√(ρ'')` [verified, Example 3.8]. Body saddles are critical points of the trilinear, and `P-137` proved there are at most 2 per cell; the trilinear's critical points track the field's as `h → 0`. (C1) On `gaussian_rf`, the number of cells with a body saddle, divided by the Kac–Rice expected number of critical points of the *field* in the box (computed via the `N=3` GOE expectation the paper gives), converges to a constant as `cells_per_correlation_length` grows — `ratio_stability_across_h` under **10%** between 65³ and 129³. (C2) That constant is the same across **8** seeds to within **15%**, so the ambiguous-cell fraction of a field is predictable from `ρ'(0)`, `ρ''(0)` and `h` alone. SHARE: none; a prediction.

**Falsified by.** C1 by the ratio drifting with `h`, which would mean the trilinear's saddles are not tracking the field's — plausible, and it would bound how far below the correlation length the `2×2×2` model can be trusted. C2 by seed-to-seed spread above 15%, which would say the count is not concentrated and `P-223`'s Poisson row inherits the question. VACUITY CONTROL: `body_saddle_cells_measured` must be non-zero at every rung and `P-137`'s degree-2 eliminant must be the saddle instrument, not a re-transcription.

#### P-183 — registered for R-188, before the harness: The expected area of a Gaussian level set is `2√λ₂/π` per unit volume, and `M-13`'s constant should be it

**Ticket:** `R-188` (S). **Records:** `seed`, `resolution`, `lambda2`, `area_expected_closed_form`, `area_mesh`, `area_zk_weights`, `area_coarea`, `ratio_mesh_over_expected`, `ratio_zk_over_expected`, `derivation_check_symbolic`, `c1_holds`, `c2_holds`.

**Hypothesis.** For a stationary unit-variance Gaussian field with `Var(∂ᵢX) = λ₂ = −2ρ'(0)`, the Kac–Rice formula for the level set gives `E[Area({X=0})] / Vol = E|∇X| · φ(0) = √λ₂ · 2√(2/π) · (2π)^{−1/2} = 2√λ₂/π` **[derived here, not quoted — the row must check it]**. `M-13`'s area law fits a constant per field with no prediction. (C1) On `gaussian_rf` over 32 seeds, the mesh area per unit volume agrees with `2√λ₂/π` within **2 standard errors**, and so do `P-177`'s two mesh-free estimators. (C2) The derivation is checked symbolically in the bench (a committed symbolic integration of `‖z‖` against the 3-D standard normal reproduces `2√(2/π)`, and the product reproduces `2√λ₂/π`), so the constant is a theorem in the ledger rather than a recalled formula. SHARE: none.

**Falsified by.** C1 by a systematic offset, which would say either the mesh area is biased (Group A grades that separately) or the fixture's `λ₂` is not what the construction claims. C2 by the symbolic check disagreeing with the closed form as written here, in which case this document is wrong and the row's value is that it says so. VACUITY CONTROL: `derivation_check_symbolic` is the vacuity control — it must be a computed `true`, not a constant.

#### P-184 — registered for R-189, before the harness: A null registered on purpose — the value-noise percolation threshold is not the Gaussian one

**Ticket:** `R-189` (M). **Records:** `field`, `seed`, `resolution`, `percolation_isovalue`, `percolation_isovalue_grf`, `bracket_low`, `bracket_high`, `threshold_negative`, `giant_component_share_at_zero`, `c1_holds`, `c2_holds`.

**Hypothesis.** `M-489` measured the percolation isovalue on `noise_cavity` at **0.029340**. For Gaussian fields in `d ≥ 3` the level-set percolation threshold is a property of the covariance class and is strictly **negative** for the standard smooth examples (the Bargmann–Fock and Gaussian-free-field results the 2026-08-29 sweep cited). (C1) On `gaussian_rf` at matched correlation length, `percolation_isovalue_grf` is **negative** on 8 of 8 seeds — the giant component of `{X ≥ 0}` already exists at the zero level. (C2) `noise_cavity`'s positive **0.0293** is therefore *not* the Gaussian value, and the row records the sign difference as the second quantitative non-Gaussianity of the value-noise fields after `P-181` C3. SHARE: none.

**Falsified by.** C1 by a positive threshold on the GRF, which would either contradict the literature (unlikely) or say the fixture's correlation length is so short relative to the box that finite-size effects dominate — in which case `resolution` must be pushed until the sign settles, and the row says at what size. C2 by `noise_cavity` agreeing with the GRF after all, a strictly better result. VACUITY CONTROL: the bracket columns must narrow with resolution on both fields, or the threshold is not being located.

---

# Group C — digital geometry (arithmetic planes, multigrid convergence)

#### P-185 — registered for R-190, before the harness: Greedy meshing is finished only for axis-aligned planes, and digital plane recognition finishes the rest

**Ticket:** `R-190` (M). **Records:** `rotation_deg`, `axis`, `resolution`, `triangles_mc`, `triangles_greedy_axis_aligned`, `planes_recognised`, `cells_covered_by_planes`, `coverage_share`, `triangles_plane_quads`, `max_plane_deviation_cells`, `recognition_ms_share`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `M-478` found greedy meshing within **1.12%** of Eppstein's optimum and declared the stage finished — for axis-aligned faces, which is all greedy meshing can see. A flat wall at any other orientation leaves a *digital plane* `{μ ≤ ax+by+cz < μ+ω}` in the sign grid, and recognising digital planes is a solved problem with a review of its own [Brimkov, Coeurjolly & Klette 2007 — landing page; body unread]. (C1) On `box_exact` rotated by **10°, 20°, 35°** about one axis at 65³, a greedy digital-plane recogniser covers at least **80%** of surface cells with recognised planes whose true-face deviation is under **1 cell**. (C2) Emitting one quad per recognised plane (clipped to the chunk) cuts triangle count below Marching Cubes' by at least **5×** on those rows, against an axis-aligned greedy arm that cannot apply. (C3) Recognition costs under **50%** of an extraction. SHARE: C2 moves the triangle count of flat-walled fields only; the share of surface cells on flat walls must be reported per field and it is near zero on `sphere`, `torus` and `gyroid`.

**Falsified by.** C1 by coverage under 80%, which would say the sign grid of a rotated box is not a clean digital plane at this `h` — `P-186` then measures why. C2 by the quads not saving triangles, which cannot happen if C1 holds and is a harness bug. C3 by cost above 50%, which prices the mechanism out of the hot path and leaves it for the CAD consumer. VACUITY CONTROL: at **0°** the recogniser must reproduce `M-478`'s greedy rectangle count to within 5%, or it is not recognising planes.

#### P-186 — registered for R-191, before the harness: The staircase has a period, and the period is a continued fraction

**Ticket:** `R-191` (S). **Records:** `rotation_deg`, `slope_rational_approx`, `continued_fraction_depth`, `staircase_period_predicted`, `staircase_period_measured`, `needle_ratio_mean`, `needle_ratio_predicted`, `c1_holds`, `c2_holds`.

**Hypothesis.** A naive digital plane's staircase repeats with a period set by the slopes' continued-fraction expansions (the 2-D case — digital straight segments — is classical). `M-309`'s `0.250` needle ratio on `box_exact` is the axis-aligned limit of a function of rotation angle that this structure predicts. (C1) On rotated `box_exact` the measured staircase period along the wall (run-length of identical case indices) equals the continued-fraction prediction on **every** angle tested, to ±1 cell. (C2) Marching Tetrahedra's mean needle ratio as a function of angle is predicted by the staircase period to within **0.05** — the registration expects the ratio to *improve* off-axis, since the needles are an axis-aligned pathology. SHARE: none; a model.

**Falsified by.** C1 by periods that do not match, which would mean the sign grid of a trilinear-sampled plane is not the naive plane of the definition (thickness `ω` differs), and the row records the measured `ω`. C2 by the needle ratio being angle-independent, which would say the slivers are not a lattice-alignment effect and Treece's randomised split (novelty table #3) is the wrong remedy. VACUITY CONTROL: at least two angles must produce different predicted periods, or C1 is a single-point match.

#### P-187 — registered for R-192, before the harness: The normal estimator has a multigrid-convergence order, and it is the order DC's QEF inherits

**Ticket:** `R-192` (M). **Records:** `field`, `resolution`, `normal_estimator`, `normal_error_l2`, `normal_error_max`, `fitted_order`, `qef_vertex_error`, `qef_vertex_order`, `c1_holds`, `c2_holds`.

**Hypothesis.** The crate's Hermite normals are gradients of the trilinear; digital geometry states normal-estimator accuracy as *multigrid convergence* with a proved order — Cuel, Lachaud & Thibert 2014 (in the corpus) give an explicit bound with an `h^{1/2}` term for a Voronoi-covariance estimator [verified, Theorem 1]. (C1) The trilinear-gradient normal converges at fitted order **1** on smooth fields (its truncation error is `O(h)`), and the Voronoi-covariance estimator at its published order, measured on the same rows. (C2) Dual Contouring's vertex error order equals the normal error order **plus one** on `sphere` and `torus`, i.e. the QEF is normal-limited, which says a higher-order normal is the lever for DC that `✗116` said no filter is for MC. SHARE: none; orders.

**Falsified by.** C1 by the gradient normal converging faster than 1, which would be a pleasant surprise about the trilinear at the crossing point (superconvergence, Group D). C2 by the vertex order not tracking the normal order, which would say the QEF error is dominated by the *positions* of the Hermite samples, not their normals. VACUITY CONTROL: a deliberately degraded normal (one-sided difference) must lower the fitted order, or the instrument is not measuring the estimator.

---

# Group D — output representation order (CAGD, higher-order FEM, numerical analysis)

#### P-189 — registered for R-193, before the harness: Vertex error and surface error have different orders, and nobody has measured them apart

**Ticket:** `R-193` (S). **Records:** `field`, `resolution`, `filter`, `vertex_error_rms`, `vertex_error_order`, `barycenter_error_rms`, `barycenter_error_order`, `edge_midpoint_error_order`, `order_gap`, `c1_holds`, `c2_holds`.

**Hypothesis.** `✗116` measured a tricubic filter reconstructing at order **3.84–4.04** while the mesh converged at **1.77–1.98**, and concluded the second order belongs to the PL mesh. That conclusion has an untested corollary: the *vertices* of that mesh — points on the tricubic's zero set — should carry the filter's order, and only the flat triangles between them the PL order. (C1) Under the tricubic filter, `vertex_error_order` is **≥ 3.5** on `sphere` and `torus` while `barycenter_error_order` is **≤ 2.2** — an `order_gap` above 1.3. (C2) Under the shipped trilinear filter both orders are **2 ± 0.2**, so the gap is created by the filter and consumed by the representation. SHARE: none; reuses `p-157.csv`'s arms.

**Falsified by.** C1 by the vertices converging at order 2 under the tricubic, which would locate `✗116`'s bottleneck in *vertex placement along the edge* (the crossing-position lottery of `✗42`) rather than in the representation, and would retire `P-188`, `P-190` and `P-191` before they are built — the cheapest possible outcome for Group D and the reason this row runs first. C2 by a gap under the trilinear, which would be a superconvergence result in its own right. VACUITY CONTROL: an infinite-order `exact_oracle` arm (as in `✗116`) must show the same vertex/barycenter split, or the split is a filter artefact.

#### P-188 — registered for R-194, before the harness: Curved PN triangles from the Hermite data DC already has, and the order they buy

**Ticket:** `R-194` (L). **Records:** `field`, `resolution`, `extractor`, `patch_kind`, `tessellation_level`, `hausdorff_pl`, `hausdorff_patch`, `fitted_order_pl`, `fitted_order_patch`, `triangles_pl`, `triangles_patch_tessellated`, `error_at_matched_triangles`, `normal_source`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** Vlachos et al.'s PN triangle needs three vertices and three normals; Dual Contouring's `HermiteCell` has both; Wang & Olano applied PN triangles to Marching Cubes in 2015 (paywalled — prior art, unread, cited as such). `✗116` says the PL representation caps the order at 2. (C1) With PN-triangle patches built from DC's vertices and Hermite normals, tessellated at level 4, the Hausdorff exponent on `sphere` and `torus` rises to **≥ 2.8** from PL's `~1.9`. (C2) At *matched output triangle count* — the tessellated patch against a PL mesh at the finer resolution that produces the same count — the patch arm's Hausdorff error is lower on **3 of 4** smooth fields. (C3) On `box_exact` and `csg_difference` the patch arm is **no worse** than PL: the crease normals from DC's two-vertex cells must not bulge the patches across the crease. SHARE: C2 is the claim that matters; C1 is the mechanism. The cost of tessellation is reported but not claimed against.

**Falsified by.** C1 by the exponent staying near 2, which — if `P-189` C1 held — means the patch construction is wrong, and if `P-189` C1 failed means the vertices were never better than order 2 and the row should not have run. C2 by PL winning at matched count, which is the honest likely outcome on fields with any curvature variation and would scope patches to the CAD consumer. C3 by crease bulging, which names the exact defect a future row must fix. VACUITY CONTROL: `normal_source` must be run as both `hermite` and `finite_difference`, and the two must differ in C1, or the Hermite data is not load-bearing.

#### P-190 — registered for R-195, before the harness: The triquadratic's isosurface is a rational-quadratic patch per cell, and the table is tabulable

**Ticket:** `R-195` (M). **Records:** `field`, `resolution`, `cells_quadratic_fit_ok`, `patch_rank_deficient`, `hausdorff_patch`, `fitted_order`, `patch_eval_cost_share`, `crease_cells_excluded`, `c1_holds`, `c2_holds`.

**Hypothesis.** `M-451` priced the triquadratic rung: mixed volume **29**, case space `2^27`, *"the only rung that has both order and a tabulable table"*. Wiley et al. 2003 show a quadratic element's isosurface *"is represented as a rational-quadratic patch in parameter space and a rational-quartic patch in physical space"* [abstract; free copy at escholarship, not fetchable by DOI]. (C1) Fitting a triquadratic per cell from the 27-point stencil and emitting Wiley's rational-quadratic patch, the Hausdorff exponent on `sphere` and `torus` is **≥ 2.8**. (C2) The patch evaluation costs under **3×** an extraction, and the rank-deficient patches (`patch_rank_deficient`) are confined to cells `P-171` already flags as non-transverse. SHARE: none claimed; the row exists to price a rung `M-451` said was tabulable.

**Falsified by.** C1 by order 2, which would put this row in the same bin as `P-188`'s failure and say the bottleneck is upstream. C2 by rank deficiency spreading beyond the non-transverse cells, which would say a `3×3×3` stencil is too small for a stable quadratic fit on real fields. VACUITY CONTROL: `cells_quadratic_fit_ok` must be under 100% on at least one field, or the fit is never stressed.

#### P-191 — registered for R-196, before the harness: A null registered on purpose — Richardson extrapolation of vertex positions across two resolutions

**Ticket:** `R-196` (S). **Records:** `field`, `resolution_pair`, `shared_edges`, `vertex_error_h`, `vertex_error_h2`, `vertex_error_extrapolated`, `observed_order_from_pair`, `extrapolation_gain`, `expansion_regular_share`, `c1_holds`, `c2_holds`.

**Hypothesis.** If the vertex error has a regular expansion `e(h) = C h² + O(h⁴)` then `(4v_{h/2} − v_h)/3` on the edges the two grids share is `O(h⁴)`, for free. `✗42` measured the crossing-position lottery — the error's *sign* depends on where in the cell the surface crosses — which is precisely a failure of regularity, so the registration expects the null. (C1) On `sphere`, `expansion_regular_share` (edges where the two-resolution error ratio is within 20% of 4) is **under 60%**, and the extrapolated vertex error is **not** better than `vertex_error_h2` on more than 60% of shared edges. (C2) On the minority of regular edges the extrapolation gain exceeds **3×**, so the mechanism is real where its hypothesis holds and a per-edge regularity test would be the gate. SHARE: none.

**Falsified by.** C1 by the extrapolation working on most edges, a strictly better result that would say `✗42`'s lottery averages out at the vertex level. C2 by no gain even on the regular edges, which would say the expansion's leading term is not `h²` at the vertex — consistent with `P-189` finding a higher vertex order. VACUITY CONTROL: `observed_order_from_pair` must be reported per edge and its distribution must have non-zero spread, or the regular/irregular split is a constant.

---

# Group E — reach, tubes and sampling (geometric measure theory, learning theory)

#### P-193 — registered for R-197, before the harness: The reach of an SDF is where its tube polynomial breaks, and that is a measurement not a lookup

**Ticket:** `R-197` (S). **Records:** `field`, `resolution`, `offset_radii`, `volume_of_offset`, `cubic_fit_residual`, `breakpoint_radius`, `reach_closed_form_or_none`, `reach_medial_axis_estimate`, `intrinsic_volumes_from_fit`, `c1_holds`, `c2_holds`.

**Hypothesis.** Weyl–Steiner: `Vol({x : d(x,M) ≤ r}) = Σ_j L_j(M) ω_{k−j} r^{k−j}` for `r` small enough — a cubic in `r` whose coefficients are the intrinsic volumes, valid below the reach [verified, Taylor 2006 eq. 1.1]. For an SDF, `{f ≤ r}` is that offset with no computation. (C1) On `sphere` (reach `R`) and `torus` (reach `r_minor`) the volume of `{f ≤ r}`, counted on the grid, is fitted by a cubic in `r` with residual under **0.5%** up to the closed-form reach and the residual grows past it, so `breakpoint_radius` is within **one cell** of the true reach at 65³. (C2) The cubic's fitted coefficients reproduce the closed-form area and integrated mean curvature of `sphere` within **3%**, so the tube fit is a *fourth* mesh-free instrument for Group A's quantities. SHARE: none.

**Falsified by.** C1 by no breakpoint, which on `torus` would mean the offset volume stays cubic past `r_minor` — impossible for the true torus, so it would be a sampling artefact and the row records the resolution at which it appears. C2 by the coefficients disagreeing with Group A's, which would say one of the two integral-geometric routes was transcribed wrong. VACUITY CONTROL: `box_exact` (reach 0, a polyhedron) must **fail** the cubic fit at every `r > 0`, or the fit is not sensitive to the reach.

#### P-192 — registered for R-198, before the harness: `graph_cube_g5`'s hole at 9³ is a sampling-density event, and the reach predicts it

**Ticket:** `R-198` (M). **Records:** `field`, `genus`, `reach_measured`, `min_feature_size`, `resolution`, `h_over_reach`, `chi_correct`, `first_correct_rung`, `predicted_rung_from_reach`, `nonmonotone_rungs`, `c1_holds`, `c2_holds`.

**Hypothesis.** `M-453`: *"`graph_cube_g5`, the highest-genus fixture, is perfect at 7³ and 11³ and meshes literally nothing at 9³"*. Niyogi, Smale & Weinberger (in the corpus at 0.69, cited nowhere) give a sampling-density condition in terms of the reach under which homology is recovered. (C1) With `reach_measured` from `P-193` on all six prescribed-genus fields, the first rung at which χ is correct is predicted by a single threshold on `h_over_reach` — the same threshold for all six fields, to within one rung. (C2) The **non-monotone** rungs (correct, wrong, correct) all lie at `h_over_reach` values where a thin feature's diameter is within **±0.25 cells** of an integer number of cells — i.e. the 9³ hole is a *commensurability* event of the kind Group M names, not a density one, and the reach threshold governs only the *first* correct rung. SHARE: none.

**Falsified by.** C1 by no common threshold, which would say the six fields' failures are not governed by one regularity number — plausible for the graph-based ones whose thinnest feature is a tube, where the relevant number is the tube radius, not the reach; the row reports `min_feature_size` beside `reach_measured` for that reason. C2 by the non-monotone rungs not being commensurability events, which would leave `M-453`'s sharpest finding unexplained and worth a row of its own. VACUITY CONTROL: `nonmonotone_rungs` must be non-empty on at least one field, or C2 is about nothing.

#### P-194 — registered for R-199, before the harness: A topology certificate computed before extraction, from the reach and the grid spacing

**Ticket:** `R-199` (M). **Records:** `field`, `resolution`, `reach_measured`, `h`, `nsw_condition_satisfied`, `chi_correct_measured`, `certificate_true_and_wrong`, `certificate_false_and_right`, `certificate_cost_share`, `c1_holds`, `c2_holds`.

**Hypothesis.** `CLAUDE.md` says χ is unassertable outside `sphere` and `torus`; `P-145` added an oracle that reads the field; NSW's theorem adds a *certificate* that needs only the reach and `h`. (C1) On all fourteen post-`R-177` reference fields at 17/33/65/129, `certificate_true_and_wrong` is **0**: whenever the NSW condition is satisfied, the extracted χ equals the prescribed or oracle χ. (C2) The certificate is conservative — `certificate_false_and_right` is non-zero on at least half the rows — and the row reports how conservative, because a certificate that never fires is a certificate that costs `P-193`'s reach measurement for nothing. SHARE: C1 is a soundness clause and moves nothing.

**Falsified by.** C1 by a single true-and-wrong row, which is the important outcome: it would mean the NSW hypotheses (a `C²` manifold, the sample within the reach) are not met by the field or the extractor, and the row must say which. C2 by the certificate being tight, a strictly better result. VACUITY CONTROL: `box_exact` and `csg_difference` (reach 0) must produce `nsw_condition_satisfied = false` on every row, or the certificate ignores the reach.

---

# Group F — regularity predicts adaptivity (nonlinear approximation, fractal geometry)

#### P-195 — registered for R-200, before the harness: The Besov index of the field predicts `M-474`'s `s`, from Haar coefficients alone

**Ticket:** `R-200` (M). **Records:** `field`, `resolution`, `haar_levels`, `coefficient_decay_exponent`, `besov_index_estimate`, `s_adaptive_p161`, `s_uniform_p161`, `s_difference_p161`, `rank_correlation`, `c1_holds`, `c2_holds`.

**Hypothesis.** DeVore's survey (paywalled; its abstract states the theorem) characterises the rate of nonlinear approximation by a Besov regularity that is *"significantly weaker than required in the linear theory"*, and Besov regularity is read off wavelet coefficient decay. `M-474` measured `s` for both families on eight fields over two decades of field evaluations — the hard way. (C1) A Haar decomposition of the sampled field at 129³ gives a `besov_index_estimate` per field that rank-correlates with `M-474`'s `s_adaptive` at **≥ 0.7** across the eight fields. (C2) The *difference* between the Besov and Sobolev indices — the regularity gap the nonlinear theory says adaptivity can exploit — rank-correlates with `s_difference` at **≥ 0.7**, so `lod_helps` is predictable from the field's coefficients before an octree is built. SHARE: none; a predictor.

**Falsified by.** C1 by correlation under 0.7, which with eight fields is easy to fail and would say the Haar estimate at this resolution is too coarse — the row reports `haar_levels` so a finer sweep can follow. C2 by the gap not predicting `s_difference`, which would say `M-474`'s field-dependence is not a regularity effect. VACUITY CONTROL: `fbm_terrain`'s estimated index must come out *below* `sphere`'s, or the estimator is not measuring roughness.

#### P-196 — registered for R-201, before the harness: The word *fractal* pays its exponent — `fbm_terrain`'s triangle count scales as `h^{−(3−H)}`

**Ticket:** `R-201` (S). **Records:** `field`, `hurst_H_configured`, `hurst_H_measured`, `resolution`, `surface_cells`, `triangles`, `fitted_count_exponent`, `predicted_exponent`, `area_law_exponent_m13`, `c1_holds`, `c2_holds`.

**Hypothesis.** `M-13`'s area law says surface-cell count scales as `h^{−2}` times a per-field constant. Fractional Brownian terrain with Hurst exponent `H` has a level set of box-counting dimension `3 − H`, so its cell count must scale as `h^{−(3−H)}` and the area law must **fail** on it by exactly `1 − H` in the exponent — the axes document v2 listed *fractal (no measured exponent)* as folklore; this is the exponent. (C1) On `fbm_terrain` across 17/33/65/129/257, `fitted_count_exponent` equals `3 − H` within **0.1**, where `H` is measured from the field's own structure function (not taken from the generator). (C2) On the seven non-fractal fields the exponent is **2 ± 0.1**, so the deviation is specific to the fractal field. SHARE: none.

**Falsified by.** C1 by an exponent of 2, which would say the hashed lattice noise is smooth below its finest octave and the "fractal" description is wrong at these resolutions — `P-176` already found the field is not a Gaussian fBm, so this is a live possibility and worth one afternoon. C2 by another field showing a non-2 exponent, which would be a finding about that field. VACUITY CONTROL: `hurst_H_measured` must be reported with its fit's confidence interval, and the generator's configured `H` must lie outside that interval on at least the finest rung, or the two numbers are the same number.

#### P-197 — registered for R-202, before the harness: A null registered on purpose — the octree's tree constraint costs a constant, not a rate

**Ticket:** `R-202` (M). **Records:** `field`, `budget`, `error_free_nterm`, `error_tree_nterm`, `error_octree_shipped`, `tree_over_free_ratio`, `octree_over_tree_ratio`, `ratio_stable_across_budget`, `c1_holds`, `c2_holds`.

**Hypothesis.** Cohen, Dahmen, Daubechies & DeVore (paywalled) prove tree-structured `N`-term approximation achieves the same rate as unconstrained `N`-term approximation, at the cost of a constant. `M-473` measured the octree's gain against uniform; nothing measured its loss against the *unconstrained* best. (C1) On the four `Exact` fields at matched sample budget, the error of the best free `N`-term Haar approximation and of the best tree-constrained one differ by a `tree_over_free_ratio` that is stable across budgets (**±20%**), i.e. a constant. (C2) The shipped greedy octree sits within a further **2×** of the tree optimum, so the remaining headroom for adaptivity is a constant of at most `tree_over_free × 2` — a number the ledger does not have. SHARE: none; a denominator for Axis 7.

**Falsified by.** C1 by the ratio growing with budget, which would contradict the theorem's hypotheses being met and would say these fields are outside its Besov class. C2 by the greedy octree being far from the tree optimum, which would be a live lever. VACUITY CONTROL: `error_free_nterm` must be strictly below `error_tree_nterm` on every row, or the free arm is not free.

---

# Group G — the denominator for Phase 25 (information theory)

#### P-198 — registered for R-203, before the harness: The case-index stream has an entropy rate, and no bit-packing beats it

**Ticket:** `R-203` (S). **Records:** `field`, `resolution`, `scan_order`, `symbols`, `h0_marginal_bits`, `h1_conditional_bits`, `h2_conditional_bits`, `lz78_bits_per_symbol`, `entropy_rate_estimate`, `phase25_best_bits_per_cell`, `gap_to_floor_ratio`, `c1_holds`, `c2_holds`.

**Hypothesis.** Phase 25 measured twenty bit-packing numerators (`R-103`–`R-122`) and no denominator. The case-index stream along a scanline is a stationary source; its entropy rate bounds every lossless encoding of it, and the Lempel–Ziv code length per symbol converges to that rate (Ziv & Lempel 1978, paywalled — the theorem is standard). (C1) On every field at 65³, the order-2 conditional entropy `H(X_n | X_{n−1}, X_{n−2})` and the LZ78 length per symbol agree within **15%**, and both are below the marginal entropy `H₀` by at least **30%** — the stream is compressible beyond a per-symbol code. (C2) Phase 25's best measured bits-per-cell is within **1.5×** of the entropy-rate estimate on the majority of fields, which would say the bit-packing line is near its floor and should stop — the most valuable outcome. SHARE: none — this row is the denominator, and it reports `gap_to_floor_ratio` per field.

**Falsified by.** C1 by the two estimates disagreeing, which would mean the stream has long-range structure order-2 conditioning misses (Group Z's spectral row is the follow-up). C2 by a gap above 1.5×, which would say Phase 25 left a measurable amount on the table and names how much. VACUITY CONTROL: a scrambled-order arm (random permutation of cells) must raise the LZ78 length toward `H₀`, or the estimator is insensitive to order.

#### P-199 — registered for R-204, before the harness: Adjacent cells share four corners, so the case-index transition is a walk on a 16-out graph and its entropy is at most 4 bits

**Ticket:** `R-204` (S). **Records:** `field`, `resolution`, `axis`, `transitions_observed`, `transitions_possible_16`, `h1_conditional_bits`, `bound_4_bits`, `predictive_coder_bits_per_cell`, `arithmetic_coder_bits_per_cell`, `c1_holds`, `c2_holds`.

**Hypothesis.** Along any axis, the next cell's case index shares the four corner signs of the shared face with the current one, so at most `2⁴ = 16` successors are possible and `H(X_n | X_{n−1}) ≤ 4` bits by a counting argument, against `H₀ ≤ 8` — a **de Bruijn-like** transition structure the bit-packing sweep never used. (C1) `transitions_observed` never exceeds 16 per predecessor on any field (a correctness check on the encoding), and `h1_conditional_bits` is under **4** on every row. (C2) A predictive coder that transmits only the four *new* bits per cell, arithmetic-coded against the empirical conditional distribution, achieves bits-per-cell within **10%** of `P-198`'s entropy-rate estimate. SHARE: C2's bits-per-cell is compared to Phase 25's best; the row claims nothing about speed.

**Falsified by.** C1 by more than 16 successors, which would be a bug in the case-index encoding and worth finding. C2 by the coder sitting far from the rate, which would mean the order-1 model is too weak and the residual structure is across rows, not along them. VACUITY CONTROL: the harness must show `transitions_observed` reaching 16 for at least one predecessor on `gyroid`, or the bound is not being approached.

#### P-200 — registered for R-205, before the harness: A null registered on purpose — the mesh is far above the metric-entropy floor of the surface class, and by how much

**Ticket:** `R-205` (M). **Records:** `field`, `resolution`, `epsilon`, `mesh_bits_position`, `mesh_bits_index`, `mesh_bits_total`, `sign_grid_bits`, `metric_entropy_bound_bits`, `mesh_over_floor_ratio`, `signgrid_over_floor_ratio`, `c1_holds`, `c2_holds`.

**Hypothesis.** Axis 14 named ε-entropy — the number of bits any representation of a surface to accuracy ε must carry — and nobody computed it. For a `W^{2,∞}`-bounded surface in a unit box the ε-entropy scales as `ε^{−1}` (dimension 2, smoothness 2: `ε^{−d/s} = ε^{−1}`), with a constant the field's Hessian bound sets. (C1) The emitted mesh (positions at `f32`, indices at `u32`) carries at least **50×** the ε-entropy bound at its own Hausdorff ε on every `Exact` field, and the sign grid alone at least **10×**. (C2) The ratio is stable across resolutions to ±30%, i.e. it is a constant of the representation, not of the resolution. SHARE: none — a denominator that says how much a *future* compressed representation could gain, and the registration expects the answer to be "a lot", which is why the row is cheap and registered as a null-shaped result.

**Falsified by.** C1 by the mesh being within 50× of the floor, a surprise that would say the mesh is already a near-optimal code for the surface — worth knowing. C2 by a ratio that grows with resolution, which would say the mesh's redundancy is resolution-dependent and the compressibility story is about coarse grids only. VACUITY CONTROL: `metric_entropy_bound_bits` must be computed from the field's measured Hessian bound and reported with it, not from a constant.

---

# Group H — time (kinetic data structures)

#### P-201 — registered for R-206, before the harness: A sliding brush is a kinetic problem, and the number of certificate failures per frame is the incremental cost

**Ticket:** `R-206` (L). **Records:** `field`, `brush_radius_cells`, `speed_cells_per_frame`, `frames`, `certificates_total`, `failures_per_frame_mean`, `failures_per_frame_max`, `cells_remeshed_per_frame`, `cells_full_remesh`, `event_cost_share`, `output_identical_to_full`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** v1's Axis 8 read *"full re-mesh"* in every column. A KDS maintains a structure under continuous motion by certificates with failure times (Basch, Guibas & Hershberger 1999, paywalled; Acar et al.'s kinetic convex hulls and kinetic mesh refinement are in the corpus). For a brush SDF `b(x − vt)` combined by `min` with a static field, the certificate for a corner is *"sign(f) unchanged"* and its failure time is the root of a one-dimensional function of `t`. (C1) Over 120 frames of a brush of radius 8 cells moving at 0.5 cells/frame through `noise_cavity` at 65³, the number of certificate failures per frame is under **2%** of the surface cells, and remeshing only the cells whose certificates failed produces output **byte-identical** to a full remesh on every frame. (C2) `failures_per_frame_max` is under **3×** the mean — the worst frame is bounded, which `M-322`'s adversarial bisect was not. (C3) Event processing plus local remesh costs under **20%** of a full remesh per frame. SHARE: C3 is the claim; the share of a frame that is remesh must be reported from `edit_trace`'s numbers.

**Falsified by.** C1 by any non-identical frame, which would mean a certificate was missed — the union-find of `R-022` is the obvious culprit and the row names it. C2 by a worst frame well above the mean, which is the adversarial case again and would say the KDS bounds *events*, not their cost. C3 by cost above 20%, which would price the mechanism against a full remesh and settle Axis 8 for this workload. VACUITY CONTROL: a brush that moves **0** cells/frame must produce 0 failures, and a brush that jumps its full diameter per frame must produce failures on essentially every cell it covers, or the certificate set is wrong.

#### P-202 — registered for R-207, before the harness: Topology changes at discrete times, and the sealing check need only run at them

**Ticket:** `R-207` (M). **Records:** `field`, `brush_path`, `frames`, `topology_events_predicted`, `topology_events_observed`, `frames_with_sealing_change`, `events_missed`, `events_spurious`, `sealing_checks_saved_share`, `c1_holds`, `c2_holds`.

**Hypothesis.** Edelsbrunner et al.'s time-varying Reeb graph (paywalled) computes the times at which the topology of `f(x,t)` changes. A sealing check (`validate::sealing`) that runs every frame is wasted on frames between events. (C1) Over the `P-201` brush path, the set of frames at which the extracted mesh's χ or component count changes is **exactly** the set of predicted topology events (`events_missed = 0`, `events_spurious ≤ 5%`). (C2) Running the sealing check only at events saves at least **80%** of the checks with no change in what is caught. SHARE: C2's saving is of a validation stage, not of extraction, and must be reported as such.

**Falsified by.** C1 by missed events, which would mean the event predictor (critical points of `f` in `(x,t)`) is missing the discrete events the PL extractor creates that the smooth field does not — a genuine discrepancy between the continuous and PL topology and worth a row. C2 by savings under 80%, which would say events are dense along a dig path and the check is not the cost. VACUITY CONTROL: a path that does not touch the surface must produce **0** predicted and 0 observed events.

---

# Group I — optimal transport

#### P-203 — registered for R-208, before the harness: A null registered on purpose — sliced Wasserstein between mesh and surface measures, against Hausdorff

**Ticket:** `R-208` (M). **Records:** `field`, `resolution`, `extractor`, `hausdorff`, `sliced_w2`, `w2_sinkhorn_coarse`, `rank_correlation_with_hausdorff`, `topology_defect_planted`, `hausdorff_response`, `sw2_response`, `c1_holds`, `c2_holds`.

**Hypothesis.** Hausdorff distance is not a metric on *measures*: it is blind to how much surface sits where, and a mesh with a missing thin sheet can score well on it. The 2-Wasserstein distance between the mesh's area measure and the true surface's is a metric that sees mass, and its sliced form is `O(n log n)` (Peyré & Cuturi, acquired). (C1) On the four `Exact` fields at three resolutions, `sliced_w2` rank-correlates with Hausdorff at **≥ 0.9** — the null: on well-behaved meshes the two metrics agree and the new one adds nothing. (C2) With a *planted* defect — a thin sheet of `thin_plate` deleted from the mesh — Hausdorff moves by under **10%** and `sliced_w2` by over **50%**, which is the one situation where the metric earns its cost. SHARE: none; an instrument.

**Falsified by.** C1 by low correlation on clean meshes, which would say the two metrics measure different things even in the benign case and one of them is the wrong grade. C2 by `sliced_w2` not responding to the planted defect, which would retire it. VACUITY CONTROL: `topology_defect_planted` must be a committed fixture whose triangle count differs from the clean mesh by a stated amount.

#### P-204 — registered for R-209, before the harness: A null registered on purpose — displacement interpolation between LOD levels

**Ticket:** `R-209` (M). **Records:** `field`, `lod_pair`, `t`, `error_linear_blend`, `error_displacement_interp`, `error_ratio_at_half`, `precompute_ms`, `per_frame_ms`, `c1_holds`, `c2_holds`.

**Hypothesis.** LOD transitions pop; the usual fix is a linear vertex blend, which is not a geodesic in any shape space. McCann's displacement interpolation — the `W₂` geodesic between the two LOD meshes' measures, computable by Solomon et al.'s convolutional method on the grid (acquired) — is the canonical path. (C1) At `t = 0.5` between 33³ and 65³ meshes of `sphere` and `torus`, the displacement-interpolated mesh's Hausdorff error to the true surface is at least **30%** lower than the linear blend's. (C2) The precompute costs more than **10** extractions, which prices it out of the hot path regardless of C1 — the registered null. SHARE: C2 is the price and is expected to kill the row for the game; the CAD consumer is offline.

**Falsified by.** C1 by no improvement, which would say the linear blend is already near-geodesic for meshes this close. C2 by a cheap precompute, a strictly better result that reopens the row. VACUITY CONTROL: at `t = 0` and `t = 1` both arms must reproduce the endpoint meshes to `1e-6`, or the interpolant is not an interpolant.

---

# Group J — robust statistics and distribution-free grading

#### P-205 — registered for R-210, before the harness: A Huber QEF and a 2-means split on normals, against the `λ` the curl residual could not replace

**Ticket:** `R-210` (M). **Records:** `field`, `resolution`, `cells_second_vertex_m60`, `huber_delta`, `irls_iterations`, `huber_vertex_moved_cells`, `two_means_split_cells`, `auc_two_means_vs_m60`, `auc_curl_p173`, `self_intersections_before`, `self_intersections_after`, `cost_share`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `M-486`: the curl residual is free and *"does not predict `M-60`'s second-vertex cells"* (AUC **0.27–0.91**). `M-60` says two of seven fields ever need a second vertex, with no rule. A cell whose Hermite normals come from two sheets is a *mixture*, and the robust-statistics answer to a mixture in a least-squares fit is a Huber loss solved by IRLS, or a 2-means split on the normal directions. (C1) A 2-means split on the normals (two centroids on `S²`, one iteration from antipodal seeds) predicts `M-60`'s second-vertex cells with AUC **≥ 0.85** on every field where `M-60` has any, beating `P-173`'s curl AUC on every row. (C2) The Huber-IRLS QEF with `δ` at the median residual moves the vertex by more than `0.1` cells on **at least 80%** of `M-60`'s cells and on **under 5%** of the others — it acts where the mixture is. (C3) On `noise_cavity` at 33³ the Huber arm cuts `M-486`'s **388** self-intersections by at least half, and costs under **15%** of an extraction. SHARE: C3's share is the QEF stage's; `M-486`'s per-cell `λ` already got 388 → 195 and the row must beat that comparand, not the shipped one.

**Falsified by.** C1 by AUC under 0.85, which would say the second-vertex cells are not a two-cluster mixture in normal space — then they are something else and the row names it. C2 by the Huber vertex moving on ordinary cells, which would say `δ` is mis-set and the robust fit is over-reacting to trilinear noise. C3 by fewer than half, or by cost above 15%. VACUITY CONTROL: `auc_curl_p173` must be recomputed on the same rows from `p-173.csv`, not quoted, or the comparison is against a memory.

#### P-206 — registered for R-211, before the harness: Conformal prediction grades the four fields `validate::accuracy` cannot

**Ticket:** `R-211` (M). **Records:** `calibration_fields`, `test_field`, `resolution`, `alpha`, `nonconformity_score`, `interval_width`, `coverage_measured_on_exact_holdout`, `predicted_hausdorff_interval_unbounded`, `coverage_target`, `c1_holds`, `c2_holds`.

**Hypothesis.** `✗107`: *"`validate::accuracy` is meaningless where `field.bound()` is not `Exact`"*, and 20 of 40 rows read `unmeasurable`. Conformal prediction produces a prediction interval with guaranteed marginal coverage `1 − α` from a calibration set, under exchangeability alone (Angelopoulos & Bates, acquired). Per-cell features the extractor already has — `|Δ|`-normalised, `‖∇f‖`, `P-173`'s curl residual, the case index — are the covariates; per-cell Hausdorff on `Exact` fields is the label. (C1) Calibrated on three `Exact` fields and tested on the held-out fourth, the `α = 0.1` intervals cover the true per-cell Hausdorff on **88–92%** of cells — the coverage guarantee holds on a field the calibration never saw. (C2) Applied to the four ungradeable fields, the intervals are **informative**: their width is under half the cell size on at least **70%** of surface cells, so the grade is a number and not a shrug. SHARE: none; a grade.

**Falsified by.** C1 by coverage outside `[0.85, 0.95]`, which would say the fields are not exchangeable — the honest likely failure on `box_exact` (creases are a different population), and the row reports coverage per held-out field so the exchangeable subset is visible. C2 by wide intervals, which would say the covariates carry no information about error — a finding about the covariates. VACUITY CONTROL: a random-label control must produce coverage at `1 − α` with intervals as wide as the label range, or the harness is not computing conformal intervals.

---

# Group K — quantum information (the 3-tangle is the hyperdeterminant)

#### P-207 — registered for R-212, before the harness: A null registered on purpose — a W-class entanglement monotone on the cells where `|Δ|` is blind

**Ticket:** `R-212` (M). **Records:** `field`, `resolution`, `cells_delta_zero`, `cells_w_class`, `w_monotone_value`, `tau_abc_value`, `rank_correlation_with_hausdorff`, `rank_correlation_delta_p134`, `defect_rate_by_monotone_quartile`, `c1_holds`, `c2_holds`.

**Hypothesis.** Coffman, Kundu & Wootters' three-way tangle is `4|Δ|` for the `2×2×2` tensor (their abstract defines *"the three-way tangle … invariant under permutations of the qubits"*; the hyperdeterminant identity is standard). Dür, Vidal & Cirac's SLOCC classification is de Silva & Lim's orbit classification under another name, and the **W class** is exactly `P-131`'s `Δ = 0`, rank 3, border rank 2 stratum — where `|Δ|` vanishes identically and `✗103` found it uninformative on 5 of 8 fields. The physics literature has monotones that are non-zero on W states. (C1) On the cells `P-131` classified as W-class, a W-class monotone (the registration names the Dür–Vidal–Cirac `|⟨ψ|ψ̃⟩|`-based measure, computed from the normalised corner tensor) rank-correlates with per-cell Hausdorff at **≥ 0.5** — on cells where `|Δ|` correlates at exactly 0 because it is exactly 0. (C2) Outside the W class the monotone adds nothing over `|Δ|` (correlation within 0.05 of `✗103`'s), which is the registered null and the reason this is one row and not a group. SHARE: none.

**Falsified by.** C1 by no correlation on the W cells, which retires the dictionary as a source of *measures* while leaving it as a source of *names*. C2 by the monotone beating `|Δ|` outside the W class, a strictly better result. VACUITY CONTROL: `cells_w_class` must be non-zero on `csg_difference` (where `M-444` counted 2,222 rank-3/border-rank-2 hits), or C1 is over an empty set.

---

# Group L — critical points on grids (discrete Morse theory, quantum chemistry)

#### P-208 — registered for R-213, before the harness: Discrete Morse on the cubical complex is a third χ oracle, and its critical cells should be the body saddles

**Ticket:** `R-213` (M). **Records:** `field`, `resolution`, `critical_cells_by_dim`, `chi_dms`, `chi_p145_oracle`, `chi_p142_analytic`, `chi_mesh`, `saddle_cells_p137`, `critical_2cells_in_saddle_cells_share`, `cost_share`, `c1_holds`, `c2_holds`.

**Hypothesis.** Robins, Wood & Sheppard 2011 (acquired) construct a discrete Morse complex on the cubical complex of a grayscale image in linear time, with critical cells that correspond to the field's critical points. The 2026-08-23 memo ranked it #7 and it was never run. (C1) `χ = Σ_k (−1)^k c_k` over the discrete Morse critical cells equals `P-145`'s oracle on every field where that oracle reports `split_blocks = 0`, and equals `P-142`'s analytic `−8N³` on the periodic gyroid — a third oracle, this one with a *gradient field* attached. (C2) At least **80%** of critical 2-cells and 3-cells lie in cells `P-137`'s eliminant flags as body-saddle cells, so the combinatorial and the algebraic saddle instruments agree on where the topology is decided. SHARE: C1's cost share is reported against `P-145`'s **0.47–1.06×**.

**Falsified by.** C1 by χ disagreeing with the oracle where the oracle is trusted, which would mean the cubical complex's adjacency and the oracle's are different pairs — `P-178` would then say which. C2 by the critical cells and the algebraic saddles living in different places, which would be a real discrepancy between the discrete and the trilinear model of the same samples. VACUITY CONTROL: `critical_cells_by_dim` must be non-zero in dimensions 1 and 2 on `gyroid`, or the complex has collapsed to its minima.

#### P-209 — registered for R-214, before the harness: A null registered on purpose — a regular isovalue by Sard, located with the chemists' grid algorithm

**Ticket:** `R-214` (S). **Records:** `field`, `resolution`, `critical_values_found`, `isovalue_distance_to_nearest_critical`, `non_transverse_cells_at_zero_p171`, `non_transverse_cells_at_shifted`, `shift_applied`, `hausdorff_delta_from_shift`, `c1_holds`, `c2_holds`.

**Hypothesis.** `M-484`: transversality fails on **19%** of `box_exact`'s cells at 17³. Sard's theorem says the set of critical values has measure zero — so a *tiny* shift of the isovalue away from the nearest critical value makes the level set transverse almost surely. Yu & Trinkle's weight method (acquired) locates a grid function's basins and critical points in linear time with quadratic convergence [abstract]. (C1) The weight method finds the critical values of each sampled field; the shipped isovalue 0 is within `1e-3` of a critical value on `box_exact` and `csg_difference` (their faces sit *on* the isovalue by construction) and not on the smooth fields. (C2) Shifting the isovalue by `1e-4` on those two fields reduces `non_transverse_cells` by at least **90%** while moving the Hausdorff error by under **1e-3** cells — which is the null: it works, and it is the same thing as `M-317`'s exactly-zero-corner problem, so the correct landing is an exact predicate (novelty table #6), not a shift. SHARE: none.

**Falsified by.** C1 by the weight method missing the polyhedral critical values, which are degenerate (a whole face is critical) and may defeat a method built for smooth densities — the row records that as the boundary of the transfer. C2 by the shift not removing non-transversality, which would say the failures are not isovalue-coincidence and `P-171` needs another explanation. VACUITY CONTROL: `sphere` must show `non_transverse_cells = 0` before and after, or the shift is not neutral where it should be.

---

# Group M — analytic number theory (lattice points)

#### P-210 — registered for R-215, before the harness: A translating sphere's triangle count flickers with the sphere-problem exponent

**Ticket:** `R-215` (S). **Records:** `radius_cells`, `translations`, `surface_cells_mean`, `surface_cells_sd`, `triangles_mean`, `triangles_sd`, `relative_jitter`, `fitted_jitter_exponent`, `predicted_exponent_range`, `interior_count_discrepancy_rms`, `c1_holds`, `c2_holds`.

**Hypothesis.** The number of lattice points in a ball of radius `R` is `(4/3)πR³ + P₃(R²)` with `P₃(x) = O(x^{21/32+ε})` and `Ω±(x^{1/2}(log x)^{1/2})`, and Jarník's mean square `∫₀ˣ P₃² dx = c₃ X² log X` [verified, Ivić et al. §2, eqs. 2.5, 2.9, 2.10]. The surface-cell count of a digitised sphere is a shell count — a difference of two such discrepancies — so its jitter under sub-cell translation is bounded by the same exponents: RMS `~ R (log R)^{1/2}` against a mean `~ R²`. (C1) Over **64** sub-cell translations of `sphere` at radii 8, 16, 32, 64 cells, `interior_count_discrepancy_rms` grows with `R` at a fitted exponent in `[1.0, 1.32]` — Jarník's `1` with the log, Heath-Brown's `21/16 = 1.3125` as the ceiling. (C2) `relative_jitter` of the *triangle* count falls with `R` at exponent `≈ −1` (i.e. as `R^{1−2} · (log R)^{1/2}`), so a sphere of radius 32 cells has a triangle-count jitter under **3%** under motion — the number a frame-budget predictor needs. SHARE: none; a prediction about output-size variance under motion.

**Falsified by.** C1 by an exponent above 1.32, which would contradict Heath-Brown and be a bug; by an exponent well below 1, which would say the shifted-centre problem has a smaller discrepancy than the origin-centred one — a genuine number-theoretic question the row would then have data on. C2 by the triangle jitter not tracking the cell jitter, which would say triangles per cell vary under translation more than cells do. VACUITY CONTROL: `translations` must include the origin and at least one irrational offset per axis, or the sample is the rational special case Landau's theorem separates.

#### P-211 — registered for R-216, before the harness: A null registered on purpose — the gyroid's χ converges at integer cells-per-period and jitters otherwise, and that is Landau's rational/irrational dichotomy

**Ticket:** `R-216` (S). **Records:** `cells_per_period`, `rational`, `resolution`, `chi_measured`, `chi_analytic_p142`, `chi_error`, `chi_error_rms_over_phase`, `phase_offsets`, `c1_holds`, `c2_holds`.

**Hypothesis.** `M-455` reproduced `−8N³` at integer `N`; `M-457` found the noise field's χ *"still moving at 129³"*. For lattice points in ellipsoids, Landau and Jarník prove different discrepancy orders for **rational** and **irrational** quadratic forms [verified, Ivić et al. §3]. A periodic field sampled at an integer number of cells per period is the rational case; at an irrational ratio the sampling phase drifts across the box. (C1) At `cells_per_period ∈ {16, 32, 64}` the gyroid's χ equals `−8N³` exactly at every phase offset (`chi_error_rms_over_phase = 0`). (C2) At `cells_per_period ∈ {16·√2, 32·φ}` the χ error over 16 phase offsets has non-zero RMS that **does not fall** with resolution — the registered null: the irrational case never converges pointwise, only in mean, and the ledger should say so before someone reports it as a bug. SHARE: none.

**Falsified by.** C1 by any error at integer ratios, a bug. C2 by the irrational RMS falling with resolution, which would say the periodic surface's χ is more forgiving of incommensurate sampling than lattice-point theory suggests — an interesting positive. VACUITY CONTROL: `phase_offsets` must include at least 8 distinct values and the integer arm must show `chi_error = 0` on all of them, or the phase sweep is not sweeping.

---

# Group N — symmetry reduction (crystallography)

#### P-212 — registered for R-217, before the harness: Mesh one asymmetric unit of a symmetric field and replicate it, bit-identically, at `1/|G|` of the cost

**Ticket:** `R-217` (L). **Records:** `field`, `space_group`, `group_order`, `resolution`, `fundamental_domain_cells`, `full_cells`, `cell_ratio`, `replicated_hash`, `full_hash`, `hashes_identical`, `vertex_max_deviation`, `extract_ms_full`, `extract_ms_reduced_plus_replicate`, `speedup`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `M-456` classified the gyroid's primitive lattice and space group; `✗49` proved bit-exact octahedral equivariance of Marching Cubes. Together they license a construction crystallography has used for a century: compute in the **asymmetric unit** and generate the rest by the group. (C1) On `gyroid` (`Ia-3d`, 96 operations per conventional cell) and `schwarz_p` (`Im-3m`), extracting the fundamental domain and replicating by the group produces a mesh whose **golden hash is identical** to the full extraction's, after a canonical vertex-order sort — the registration predicts identity for the 48 point-group operations (`✗49`'s guarantee) and records separately whether the translational and glide operations, which `✗49` did not test, preserve it. (C2) `speedup` is at least **20×** at 129³, against a theoretical `|G|` of 96 minus replication cost. (C3) `sphere` (`O(3)` — infinite group; reduce by the 48-element cube group only) shows a `48×` cell reduction and hash identity. SHARE: C2 is a cost claim on symmetric fields only; the share of a game world that is symmetric is zero and the row says so — this is for the periodic-tile and CAD consumers.

**Falsified by.** C1 by non-identical hashes under glide or translation operations, which is the finding: it would locate the exact operation under which the extractor is not equivariant, and `✗49` would gain a caveat. C2 by replication cost eating the gain, which prices the mechanism. C3 by `sphere` failing, which cannot happen if `✗49` holds and would be a bug in the replication. VACUITY CONTROL: a deliberately mis-specified group (one operation removed) must produce a hash mismatch, or the identity check is comparing a mesh to itself.

#### P-213 — registered for R-218, before the harness: A null registered on purpose — seams on mirror planes cannot crack

**Ticket:** `R-218` (S). **Records:** `field`, `seam_plane`, `is_mirror_plane`, `seam_vertices_a`, `seam_vertices_b`, `seam_pairs`, `injective`, `crack_count`, `max_gap`, `c1_holds`, `c2_holds`.

**Hypothesis.** `M-461` found the gyroid's seam correspondence non-injective (**295 / 294 / 295**), which is `M-48` again. A chunk boundary that lies on a **mirror plane** of the field's symmetry group has a seam whose two sides are each other's images under a reflection the extractor is equivariant to, so the correspondence is the identity and the seam cannot crack — a theorem, not a weld. (C1) On `gyroid` chunked along its mirror planes, every seam has `injective = true` and `crack_count = 0` with no weld pass. (C2) Chunked along a plane that is *not* a mirror plane, the same field reproduces `M-461`'s non-injective seams — the null half, which says the mechanism is the symmetry and not the chunking. SHARE: none.

**Falsified by.** C1 by a crack on a mirror seam, which would say the extractor is not equivariant under that reflection — the same finding `P-212` C1 would make from the other side. C2 by non-mirror seams also being clean, which would say `M-461`'s defect was resolution-specific. VACUITY CONTROL: `is_mirror_plane` must be computed from the space group's operations, not asserted per row.

---

# Group O — shape calculus

#### P-214 — registered for R-219, before the harness: The shape derivative of area predicts a CSG edit's triangle delta before the remesh

**Ticket:** `R-219` (S). **Records:** `field`, `edit_kind`, `edit_radius_cells`, `triangles_before`, `triangles_after`, `delta_measured`, `delta_predicted_first_order`, `delta_predicted_with_curvature_term`, `relative_error`, `prediction_cost_share`, `c1_holds`, `c2_holds`.

**Hypothesis.** Hadamard's shape derivative of the area functional under a normal displacement `φ` is `∫ H φ dA` (mean curvature times displacement); Sokolowski & Zochowski's topological derivative (acquired) handles the case where the edit *creates* a hole. `M-13` makes triangle count proportional to area, so the derivative predicts the triangle delta of an edit to first order from data on the *old* mesh. (C1) For spherical subtractions of radius 2–8 cells from `sphere`, `torus` and `noise_cavity`, the first-order prediction is within **15%** of the measured delta on **80%** of edits; the topological-derivative term is needed exactly when the edit punches through (`edit_radius_cells` exceeds local thickness). (C2) The prediction costs under **1%** of a remesh. SHARE: none — a budget predictor for the scheduler, which is a consumer, not this crate; the row exists because the derivative is a theorem about the field the crate owns.

**Falsified by.** C1 by errors above 15%, which would say the second-order (Gaussian curvature) term is not negligible at these radii — the row reports the curvature-corrected prediction beside the first-order one for that reason. C2 by cost, which cannot happen for a surface integral over the old mesh. VACUITY CONTROL: an edit that misses the surface entirely must predict and measure a delta of **0**.

---

# Group P — stochastic sampling and dither (signal processing)

#### P-215 — registered for R-220, before the harness: A jittered grid kills the axis-aligned grain, and the Wulff plot says by how much

**Ticket:** `R-220` (M). **Records:** `field`, `resolution`, `jitter_amplitude_cells`, `needle_ratio_mean`, `needle_ratio_axis_aligned`, `wulff_spread_regular`, `wulff_spread_jittered`, `hausdorff_regular`, `hausdorff_jittered`, `hashes_deterministic`, `cost_share`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** Cook 1986 (paywalled): stochastic sampling converts aliasing into noise. A per-cell deterministic hash offset of the sample positions (bounded to `0.25` cells) is a jittered grid, and the extractor on a jittered grid is Marching Tetrahedra on a perturbed lattice — the Kuhn triangulation `M-452` proved regular survives small perturbation. Treece's randomised tet *split* (novelty table #3) randomises topology; this randomises *position*, which is the other half. (C1) On `box_exact` at 65³, the mean needle ratio moves from `M-309`'s **0.250** to above **0.5** at jitter 0.25. (C2) `P-180`'s Wulff-plot spread falls by at least **50%** under jitter, at a Hausdorff cost under **10%**. (C3) Output is deterministic (`hashes_deterministic = true` across runs), because the jitter is a pure function of grid coordinates. SHARE: C2's Hausdorff cost is the price and is reported against the regular grid on every field, including the smooth ones where jitter has nothing to fix.

**Falsified by.** C1 by the needle ratio not moving, which would say the slivers are a property of the tet split and not of the lattice alignment, and `P-186` C2 inherits the question. C2 by Hausdorff cost above 10%, which prices jitter out on smooth fields. C3 by any non-determinism, a bug. VACUITY CONTROL: jitter amplitude **0** must reproduce every golden hash, or the jittered path is not the shipped path at zero.

#### P-216 — registered for R-221, before the harness: A null registered on purpose — dithering the isovalue does not unbias the mesh area

**Ticket:** `R-221` (S). **Records:** `field`, `resolution`, `dither_amplitude`, `area_mesh_undithered`, `area_mesh_dithered_mean`, `area_true_or_zk`, `bias_undithered`, `bias_dithered`, `variance_dithered`, `c1_holds`, `c2_holds`.

**Hypothesis.** Subtractive dither makes quantisation error independent of the signal (Schuchman's condition; Wannamaker et al., paywalled). A sign grid is a 1-bit quantiser of the field, and the mesh area's orientation bias (`P-180`) is quantisation error. Dithering the *isovalue* per cell by a zero-mean amount would, if the analogy held, make the expected area unbiased. (C1) The registered null: over 32 dither seeds, the dithered mean area's bias against `P-177`'s estimator is **not** smaller than the undithered bias on `box_exact`, because the bias is *directional* (a function of the normal) and a scalar dither does not see direction. (C2) Dither *does* reduce the bias on `sphere`, where the bias averages over all normals, by at least **30%** — the half of the mechanism that works, and the reason the analogy is tempting. SHARE: none.

**Falsified by.** C1 by dither removing the box bias, a strictly better result that would give a two-line area fix. C2 by no effect on the sphere, which would say the mesh area bias is not a quantisation effect at all. VACUITY CONTROL: `variance_dithered` must grow with `dither_amplitude`, or the dither is not reaching the sign decision.

---

# Group Q — the roster is a designed experiment (statistics)

#### P-217 — registered for R-222, before the harness: A space-filling roster of reference fields, so that no future mechanism is judged by a roster accident

**Ticket:** `R-222` (L). **Records:** `field`, `curvature_ratio_max`, `besov_index_p195`, `genus`, `bound_kind`, `flat_axis_aligned_fraction`, `reach_over_h`, `design_cell`, `design_coverage_before`, `design_coverage_after`, `anisotropy_rerun_c1_winners`, `c1_holds`, `c2_holds`.

**Hypothesis.** `M-459`: *"this crate's reference-field roster is adversarial to global anisotropy"*; `M-462`: 70 of 112 rows unmeasurable because five of eight fields make two arms one grid. A roster chosen for coverage of the *invariants that decide experiments* — curvature ratio, regularity, genus, bound kind, axis-aligned fraction, reach-to-cell ratio — is a designed experiment (McKay, Beckman & Conover's Latin hypercube, paywalled; the construction is elementary). (C1) The current eight fields occupy at most **40%** of a 6-factor, 3-level design's cells; six new bench-local fields chosen by a maximin criterion raise coverage above **75%**. (C2) Re-running `P-146`'s anisotropy comparison on the designed roster produces `c1_winners` **≥ 3** where `M-459` had **0 of 4** — not because the mechanism changed but because the roster now contains fields with a flat direction *and* an `Exact` bound. SHARE: none — this row buys future rows their denominators.

**Falsified by.** C1 by the current roster already covering the design, which would say `M-459`'s complaint was about one factor, not the roster. C2 by anisotropy still losing on a roster designed to let it win, which would be the definitive negative for Group D of Phase 27 and worth having. VACUITY CONTROL: `design_coverage_before` must be computed from measured invariants (`P-195`, `P-193`), not from the fields' names.

---

# Group R — backward error and conditioning (numerical analysis)

#### P-218 — registered for R-223, before the harness: The backward error of an extraction, and the condition number that turns it into Hausdorff

**Ticket:** `R-223` (S). **Records:** `field`, `resolution`, `extractor`, `backward_error_vertex_rms`, `backward_error_face_max`, `forward_hausdorff`, `gradient_norm_min`, `condition_number_estimate`, `forward_over_backward`, `ratio_matches_condition`, `c1_holds`, `c2_holds`.

**Hypothesis.** Every accuracy number in the ledger is a *forward* error. The backward error — the smallest `‖δf‖∞` such that the mesh lies exactly on `{f + δf = 0}` — is `max |f(x)|` over the mesh, a field evaluation per vertex and no true surface needed, so it grades **every** field. Forward error is backward error times the condition number `1/‖∇f‖` at the surface. (C1) On the four `Exact` fields, `forward_over_backward` equals `1/gradient_norm_min` within a factor of **2** on every row — the condition-number relation holds. (C2) On the four ungradeable fields the backward error is reported and is within the range the `Exact` fields span, so a second mesh-free grade (after Group A's) exists for them, and `thin_plate`'s difficulty in `✗107` is explained as a **conditioning** effect (`gradient_norm_min` small where the plate is thin) rather than a mechanism failure. SHARE: none; an instrument.

**Falsified by.** C1 by the ratio not tracking the gradient, which would say the mesh error is not first-order in the field perturbation — a real statement about the extractor. C2 by backward errors on the unbounded fields outside the `Exact` range, which would be a finding about those fields. VACUITY CONTROL: the backward error must be computed with the *field*, not the trilinear interpolant, and the row must report both so the interpolation error is visible as their difference.

---

# Group S — affine arithmetic, pointed at the empty-cell test (interval analysis)

#### P-219 — registered for R-224, before the harness: A null registered on purpose — extended revised affine arithmetic on `gyroid`'s rejected share

**Ticket:** `R-224` (M). **Records:** `field`, `resolution`, `coarse_level`, `rejected_share_corner_sampling`, `rejected_share_interval`, `rejected_share_affine`, `false_rejections`, `cost_ratio_affine_over_corner`, `c1_holds`, `c2_holds`.

**Hypothesis.** `M-306`: corner sampling rejects **16.8%** of coarse cells on `gyroid` against **95.1%** on `thin_plate`. Fryazinov, Pasko & Comninos' extended revised affine arithmetic (in the corpus, 0.68) bounds procedural implicit functions more tightly than interval arithmetic. The 2026-08-26 audit rejected affine arithmetic for the *trilinear*, where the corner range is exact; the coarse-cell test evaluates the *field*, where it is not. (C1) Affine arithmetic raises `gyroid`'s rejected share from 16.8% to at least **40%** with `false_rejections = 0`. (C2) The registered null: it costs more than **5×** corner sampling per coarse cell, and the net traversal time does not fall, because `gyroid`'s surface cells are dense enough that rejection buys little — the row prices the trade rather than assuming it. SHARE: C2 is the traversal stage's share and must be reported from `M-306`'s own numbers.

**Falsified by.** C1 by a share under 40%, which would say `gyroid`'s `sin·cos` products defeat affine forms as they defeat intervals. C2 by a net win, a strictly better result. VACUITY CONTROL: `false_rejections` must be checked against a full-resolution extraction on every row, or the tighter bound could be tighter than the truth.

---

# Group T — digital topology (Alexander duality)

#### P-220 — registered for R-225, before the harness: A null registered on purpose — solid and air have dual Betti numbers, and the split blocks are where the duality fails

**Ticket:** `R-225` (S). **Records:** `field`, `resolution`, `adjacency_pair`, `b0_solid`, `b1_solid`, `b2_solid`, `b0_air`, `b1_air`, `duality_holds`, `split_blocks_p145`, `duality_failures_in_split_blocks_share`, `c1_holds`, `c2_holds`.

**Hypothesis.** For a closed surface in a ball, Alexander duality gives `b₁(solid) = b₁(air)` and relates the component counts. `P-145`'s oracle chooses a foreground adjacency and reports `split_blocks`; `P-178` names the complementary pairs. The registered null: with a complementary pair the duality holds exactly, and every failure of it under a non-complementary pair lies in a split block. (C1) Under `(6,26)` and `(26,6)`, `duality_holds` on **every** field at 33³ and 65³. (C2) Under `(26,26)` — deliberately inconsistent — `duality_holds` fails, and **100%** of the failing cells are `P-145`'s split blocks. SHARE: none.

**Falsified by.** C1 by a duality failure under a complementary pair, which would say the field is not closed in the box (the gyroid and `fbm_terrain` leave through the sides, so they are excluded by `expected_euler`'s own rule and the row says so). C2 by failures outside the split blocks, which would say the split count under-reports. VACUITY CONTROL: `(26,26)` must fail on at least one field, or the control is not a control.

---

# Group U — currents and the flat norm (geometric measure theory)

#### P-221 — registered for R-226, before the harness: A null registered on purpose — the flat norm sees a missing sheet that Hausdorff sees too, and costs an LP

**Ticket:** `R-226` (S). **Records:** `field`, `resolution`, `flat_norm_scale`, `flat_norm_distance`, `hausdorff`, `planted_defect`, `flat_norm_response`, `hausdorff_response`, `lp_solve_ms`, `cost_share`, `c1_holds`, `c2_holds`.

**Hypothesis.** Morgan & Vixie (acquired) show the flat norm between boundaries is computed by an `L¹TV` minimisation — a linear program on the grid. The flat norm is a true metric on currents that sees orientation flips and mass, unlike Hausdorff. (C1) On a planted orientation flip (one triangle strip inverted on `sphere`) the flat norm responds by more than **10×** its clean-mesh value while Hausdorff responds by **0** — the one defect class Hausdorff is structurally blind to and the orientation test already catches combinatorially. (C2) The registered null: on every other planted defect (`P-203`'s deleted sheet, a displaced vertex) Hausdorff and the flat norm respond in the same direction, and the LP costs more than **10** extractions, so the instrument is not worth shipping. SHARE: none.

**Falsified by.** C1 by no response to a flip, which would mean the current's orientation was dropped in the discretisation. C2 by a cheap LP, which would reopen the row as a validation metric. VACUITY CONTROL: the clean-mesh flat norm must be non-zero and reported, so responses are ratios and not absolutes.

---

# Group V — approximation constants (the 7× has three factors)

#### P-222 — registered for R-227, before the harness: `M-472`'s constant gap decomposes into a filter constant, a representation constant and a placement constant, and the ledger should say which is 5 of the 7

**Ticket:** `R-227` (S). **Records:** `field`, `resolution`, `constant_gap_m472`, `filter_constant_blu_unser_closed_form`, `filter_constant_measured`, `representation_constant`, `placement_constant`, `product_of_three`, `product_matches_gap`, `largest_factor`, `c1_holds`, `c2_holds`.

**Hypothesis.** `M-472` bounded the constant gap to the `W_inf^2` minimal error at **6.26–7.25**. Three stages each contribute a constant: the trilinear filter (Blu & Unser give the asymptotic constant in closed form — Part I paywalled, but the Blu paper the corpus already held states the formula), the PL representation (the interpolation constant of a triangle of the mesh's shape), and vertex placement along the edge (`✗42`'s lottery, a measurable factor). (C1) The three constants, each measured in isolation by holding the other two stages exact (an `exact_oracle` field for the filter, exact vertices for the representation, exact surface for placement), multiply to within **20%** of `M-472`'s gap on `sphere` and `torus`. (C2) The representation constant is the largest of the three on both fields, which is `✗116`'s conclusion stated as a number and is what licenses Group D. SHARE: none — a decomposition.

**Falsified by.** C1 by a product far from the gap, which would say the three stages interact and the decomposition is wrong — worth knowing before Group D spends `L` on the representation. C2 by the filter or placement constant being largest, which would redirect Group D's effort to `✗42`'s lottery or to Axis 11. VACUITY CONTROL: `filter_constant_measured` must reproduce the Blu–Unser closed form within 10%, or the filter arm is not isolated.

---

# Group W — Poisson approximation (probability)

#### P-223 — registered for R-228, before the harness: A null registered on purpose — ambiguous cells are Poisson at large correlation length

**Ticket:** `R-228` (S). **Records:** `seed`, `resolution`, `cells_per_correlation_length`, `ambiguous_cells`, `mean_over_seeds`, `variance_over_seeds`, `variance_over_mean`, `chi_variance_p144_style`, `chi_variance_predicted_from_ambiguous`, `c1_holds`, `c2_holds`.

**Hypothesis.** Chen–Stein (Arratia, Goldstein & Gordon, paywalled): rare, weakly dependent events are Poisson to a bounded total-variation error. Ambiguous cells on `gaussian_rf` are rare when the correlation length is many cells, and dependent only within a correlation length. (C1) Over **32** seeds, `variance_over_mean` of the ambiguous-cell count is within **[0.8, 1.2]** once `cells_per_correlation_length ≥ 8`, and drifts away from 1 below that. (C2) `P-144`'s per-seed χ variance is predicted by the ambiguous-cell count's Poisson variance times a per-cell χ contribution measured on the same rows, within a factor of **2** — the null being that the topological variance across seeds is *explained* by the rare-cell count and needs no further mechanism. SHARE: none.

**Falsified by.** C1 by over-dispersion at large correlation length, which would say ambiguous cells cluster and the dependence is longer-range than the covariance suggests. C2 by the χ variance not being explained, which would leave `M-457`'s spread as a live question. VACUITY CONTROL: the count's mean must exceed 20 on every row, or the Poisson test is on a handful of events.

---

# Group X — coarea and the divergence theorem (measure theory)

#### P-224 — registered for R-229, before the harness: Area from the coarea integral and volume from the divergence theorem, as the tie-breaker for Group A

**Ticket:** `R-229` (S). **Records:** `field`, `resolution`, `smoothing_epsilon`, `area_coarea`, `area_zk_p177`, `area_mesh`, `volume_grid_count`, `volume_mesh_divergence`, `volume_true_or_none`, `area_disagreement_arbitrated`, `c1_holds`, `c2_holds`.

**Hypothesis.** `∫ |∇f| δ_ε(f) dx` is a third area estimate, from the *values* rather than the signs; `∮ x·n dA / 3` over the mesh is its volume by the divergence theorem, against the voxel count. `P-177` C2 needs a tie-breaker when the mesh and the configuration estimator disagree on an ungradeable field. (C1) On the `Exact` fields the coarea area is within **3%** of truth at 65³ for `ε = 1.5h`, and the mesh's divergence-theorem volume within **1%** of the voxel count plus the surface correction. (C2) On any row where `P-177` recorded a mesh/estimator disagreement beyond the bound, the coarea value sides with one of them by a margin larger than its own 3%, so the disagreement is arbitrated. SHARE: none.

**Falsified by.** C1 by errors above 3%, which would say `δ_ε` at this `ε` is too wide for the field's gradient variation — the row sweeps `ε`. C2 by the coarea value sitting between the two, which would leave the disagreement unarbitrated and would say three instruments are not enough. VACUITY CONTROL: the divergence-theorem volume must **fail** on a mesh with a planted inverted strip (`P-221`'s fixture), or the orientation dependence that makes it a check is not present.

---

# Group Y — quasi-Monte Carlo (discrepancy theory)

#### P-225 — registered for R-230, before the harness: A null registered on purpose — eight Sobol points per coarse cell against eight corners

**Ticket:** `R-230` (S). **Records:** `field`, `coarse_level`, `probe_set`, `probes_per_cell`, `rejected_share`, `missed_surface_cells`, `koksma_hlawka_bound`, `cost_ratio`, `c1_holds`, `c2_holds`.

**Hypothesis.** Corner sampling is a `2×2×2` regular probe of a coarse cell; the Koksma–Hlawka inequality bounds a quasi-Monte Carlo probe's miss rate by the field's variation times the point set's discrepancy, and a Sobol set of 8 points has lower discrepancy than the corners. (C1) With 8 Sobol probes per coarse cell, `missed_surface_cells` (cells rejected that contain surface at full resolution) on `gyroid` is at least **30%** lower than with 8 corners at the same cost. (C2) The registered null: the corners are *exact* for the trilinear and the Sobol set is not, so on the smooth fields the Sobol arm rejects **fewer** cells and the mechanism is a wash — the only field it can help is the one whose surface is dense relative to the coarse cell, and there `P-219`'s affine arm is the comparand. SHARE: none.

**Falsified by.** C1 by no reduction in misses, which retires the row. C2 by Sobol winning on smooth fields, a surprise about the trilinear's corner exactness. VACUITY CONTROL: `probes_per_cell` must be equal across arms, or the comparison is of budgets.

---

# Group Z — the error field has a spectrum (harmonic analysis, reciprocal lattices)

#### P-226 — registered for R-231, before the harness: The mesh error's power spectrum peaks at the reciprocal lattice, and BCC moves the peaks

**Ticket:** `R-231` (M). **Records:** `field`, `lattice`, `resolution`, `error_field_resolution`, `spectrum_peaks`, `peaks_at_reciprocal_vectors_share`, `peak_power_share`, `isotropic_power_share`, `bcc_peak_set`, `cubic_peak_set`, `peak_sets_differ`, `c1_holds`, `c2_holds`.

**Hypothesis.** The signed distance from the true surface to the mesh, sampled on a fine grid, is a field whose power spectrum should concentrate at the **reciprocal lattice vectors** of the sampling lattice — the aliasing signature, which `P-180` measures in the direction domain and this row measures in the frequency domain. (C1) On `sphere` at 65³ with the error sampled at 4× resolution, at least **60%** of the error power above the isotropic floor sits within one bin of a reciprocal lattice vector of `Z³`. (C2) On BCC the peaks sit at the reciprocal (FCC) vectors instead — `peak_sets_differ = true` — and the *peak* power share falls even where `P-162` found the *total* error unchanged, which is the frequency-domain version of `P-180` C2. SHARE: none; an instrument that makes `P-162`'s null legible.

**Falsified by.** C1 by diffuse power, which would say the error is not aliasing-dominated and the lattice choice is genuinely irrelevant — closing Axis 12 for this crate. C2 by BCC's peaks not moving, a bug in the BCC arm. VACUITY CONTROL: a mesh with white-noise vertex perturbation must show a flat spectrum, or the FFT is seeing the grid rather than the error.

---

## 12. What this phase deliberately does not claim

- **No row claims a speedup on the shipped path.** Every cost clause is a *price* (`P-201`, `P-204`, `P-212`, `P-219`, `P-221`), reported so that the landing decision has a number.
- **No row re-runs a Phase 27 mechanism** except `P-217`, which re-runs `P-146`'s anisotropy on a roster designed to let it win, because `M-459` said the roster was the confound and that is a testable sentence.
- **Nothing here needs a paywalled paper to run.** Where the source is paywalled the row's claim is derived (`P-183`), classical (`P-198`, `P-209`), or uses an acquired conference version (`P-177`'s Lindblad 2002/2003).
- **The two wrongly-acquired papers are not cited anywhere in this document**, and are listed once in §1 for removal.

**Acquired and in the corpus at time of writing (23).** Taylor 2006 `10.1214/009117905000000594`; Letendre `arXiv:1408.2107`; Cheng & Schwartzman `arXiv:1511.06835`; Legland et al. `10.5566/ias.v26.p83-92`; Schladitz–Ohser–Nagel `10.1007/11907350_21`; Ziegel & Kiderlen `10.1016/j.imavis.2009.04.013`; Lindblad 2003 `10.1007/978-3-540-39966-7_33`; Lindblad & Nyström 2002 `10.1007/3-540-45986-3_24`; Schmalzing & Buchert `10.1086/310680`; Brimkov–Coeurjolly–Klette `10.1016/j.dam.2006.08.004` (HTML landing); Coeurjolly–Lachaud–Levallois `10.1016/j.cviu.2014.04.013` (HTML landing); Aamari et al. `10.1214/19-EJS1551`; Angelopoulos & Bates `arXiv:2107.07511`; Ivić et al. `arXiv:math/0410522`; Peyré & Cuturi `arXiv:1803.00567`; Solomon et al. `10.1145/2766963` (HTML); Sokolowski & Zochowski `10.1137/S0363012997323230` (HTML); Coffman–Kundu–Wootters `arXiv:quant-ph/9907047`; Dür–Vidal–Cirac `arXiv:quant-ph/0005115`; Robins–Wood–Sheppard `10.1109/TPAMI.2011.95`; Yu & Trinkle `10.1063/1.3553716`; Morgan & Vixie `arXiv:math/0612287`; Hughes–Cottrell–Bazilevs `10.1016/j.cma.2004.10.008` (HTML). **Already held and newly cited (6):** Niyogi–Smale–Weinberger `10.1007/s00454-008-9053-2`; Berenfeld et al. `10.1007/s00454-021-00290-8`; Cuel–Lachaud–Thibert `10.1007/978-3-319-09955-2_12`; Fryazinov et al. `10.1016/j.cag.2010.07.003`; Acar et al. kinetic hulls `10.1007/978-3-540-87744-8_3`; Edelsbrunner & Mücke `arXiv:math/9410209`.
