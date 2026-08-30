//! **P-166 — a null registered on purpose: greedy meshing is not a matroid
//! problem, but LOD budgeting might be.**
//!
//! Ticket: R-166. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p166
//! ```
//!
//! Writes `docs/experiments/p-166.csv`.
//!
//! # What was missing
//!
//! Two halves of one question, and neither had ever been asked here.
//!
//! **The meshing half.** `M-56` measured greedy meshing's saving over face
//! culling at **`1.70×` to `256×`** across seven reference fields
//! (`FINDINGS.md:1172`) and called it *"a property of one scene, not of the
//! algorithm"*. `M-57` then showed the merged output carries T-junctions no weld
//! can remove. Neither asked the **structural** question underneath: *what kind
//! of combinatorial problem is greedy meshing?* The answer decides which
//! guarantee, if any, is available. Greedy is exactly optimal on a matroid and
//! within `(1 − 1/e)` for monotone submodular maximisation under a cardinality
//! constraint — and greedy **meshing** is neither of those, which is registered
//! here as a negative and proved rather than asserted. `R-165` measures the
//! *optimality gap* against Eppstein's `arXiv:0908.3916` formula; this row
//! measures the *structure*, and the two do not overlap: a system can fail the
//! exchange property while greedy still attains its minimum, which is exactly
//! what the four hand fixtures below turn out to show.
//!
//! **The budget half.** `M-124` measured the amortised cost of the re-mesh queue
//! at **0.085–6.269 ms** against budgets of 0.025–8 ms, tracking the budget *to
//! within one chunk* across a 320× range, and `M-125` priced that one chunk at
//! **≈0.072 ms**. So the crate knows what a chunk budget *costs*. It has never
//! measured what a chunk budget *buys* — which chunks to spend it on, and whether
//! the marginal screen-space error a refinement removes is **diminishing**. That
//! is the only property that licenses a guarantee, and it is the property this
//! row hunts a counterexample to, because *"finding one such chunk is cheaper
//! than proving none exists"*.
//!
//! # The source, because Nemhauser, Wolsey & Fisher is paywalled
//!
//! The registration requires the `(1 − 1/e)` statement to come from *"a freely
//! available restatement"*. The primary is paywalled twice over —
//! `10.1007/BF01588971` (*Mathematical Programming* **14**, 265–294, 1978, Part I)
//! and `10.1007/BFb0121195` (*Mathematical Programming Studies* **8**, 73–87, 1978,
//! Part II) — and the local corpus holds neither.
//!
//! **The restatement actually read for this harness**, in full, is:
//!
//! > Shamak Dutta, Bahman Gharesifard and Stephen L. Smith, *Submodular
//! > Optimization with Applications to Decision and Control*,
//! > **arXiv:2606.10192v1 [math.OC]**, 8 June 2026.
//!
//! It is open access, it restates every result this row needs with proofs, and
//! four of its statements are load-bearing below. Quoted rather than paraphrased:
//!
//! - **Definition 3.1 (independence systems and matroids).** `(X, I)` is an
//!   *independence system* if (1) `∅ ∈ I` and (2) *"for every `A ∈ I`, we have
//!   that `B ∈ I` for all `B ⊆ A`"* (downward-closed). It is a **matroid** if it
//!   additionally satisfies (3) *"if `A, B ∈ I` and `|A| > |B|`, then there exists
//!   `a ∈ A \ B` such that `{a} ∪ B ∈ I`"* — **the augmentation property**, which
//!   is the exchange property C1 hunts a counterexample to.
//! - **Example 3.5 (uniform matroid).** `I = {A ⊆ X : |A| ≤ k}`. *"This
//!   corresponds to a cardinality constraint … The augmentation property is
//!   immediate."* This is the constraint family the LOD budget lives in, and it
//!   is why `is_matroid` reads **true** on every `lod_budget` row.
//! - **Theorem 4.2 (Nemhauser et al., 1978).** *"Consider the
//!   cardinality-constrained submodular maximization problem with a normalized,
//!   monotone submodular function `f : 2^X → R`, subject to `|S| ≤ k` … Then the
//!   greedy algorithm returns a set `S_G` satisfying `f(S_G)/f(S*) ≥ 1 − 1/e`."*
//!   `1 − 1/e = 0.632120558828557`, and the survey's own gloss is *"at least
//!   (1−1/e) ≈ 63% of the optimal"*.
//! - **Proposition 2.3.** `f` is submodular **iff** `Δ(x|S) ≥ Δ(x|T)` for all
//!   `S ⊆ T` and `x ∉ T` — diminishing returns. Definition 2.1 is the union /
//!   intersection form, `f(S) + f(T) ≥ f(S∪T) + f(S∩T)`. This harness implements
//!   **both** and asserts they agree, because two instruments that must agree are
//!   a check on each other.
//!
//! Two further statements are used and named but are *not* what C1 or C2 turn on.
//! **Section 3.3** exhibits the canonical augmentation failure — a knapsack with
//! `X = {a, b, c}`, `w = (2, 1, 1)`, `B = 2`, where `{a}` and `{b, c}` are both
//! feasible, `|{a}| < |{b, c}|`, *"yet adding either `b` or `c` to `A` gives total
//! weight `3 > B`"* — and that instance is reproduced here verbatim as a
//! **negative control** for the exchange tester. **Section 3.4** defines a
//! `p`-system: *"for every `Y ⊆ X`, any two inclusion-wise maximal independent
//! subsets of `Y` have sizes within a factor `p` of each other"*, with a
//! `1/(p+1)` guarantee from Fisher et al. (1978). The rectangle-packing system's
//! basis-size spread is a **lower bound** on its `p`, recorded as
//! `p_system_lower_bound`.
//!
//! **What the restatement does not cover, said plainly.** It scopes itself to
//! *maximisation*: *"For submodular minimization, we refer the reader to Bach
//! (2019)"*. Greedy meshing is a **minimisation** of cardinality, so the
//! classical home for its greedy analysis is Wolsey's submodular set-covering
//! line (`10.1007/BF02579435`, *Mathematical Programming* **23**, 1982) — which
//! this harness **names and did not read**, and therefore quotes no constant
//! from. Nothing below depends on it.
//!
//! # C1 — greedy meshing as a set system, and the counterexample
//!
//! Greedy meshing is the per-slice procedure `greedy_quads.rs:15-18` states in
//! prose and `greedy_quads.rs:238-272` implements: *"per 2D slice walk +X for a
//! run of identical voxels, extend +Y holding that width, emit one quad."* The
//! set system it lives in is forced by that:
//!
//! - **Ground set `X`** — every axis-aligned rectangle lying entirely inside the
//!   slice's face mask.
//! - **`I`** — every family of **pairwise-disjoint** rectangles from `X`.
//!
//! `∅ ∈ I` trivially, and `I` is downward-closed because a subfamily of a
//! pairwise-disjoint family is pairwise disjoint — so `(X, I)` **is** an
//! independence system by Definition 3.1(1)(2). The two axioms are checked
//! exhaustively over all `2^|X|` subsets wherever `|X| ≤ 16`.
//!
//! **A maximal member of `I` is exactly a rectangle partition of the mask**, and
//! that single observation is the counterexample. If a cell were left uncovered,
//! the `1×1` rectangle on it is in `X` and is disjoint from everything chosen, so
//! the family was not maximal. Hence every basis covers the mask — and different
//! bases have **different cardinalities**, which Definition 3.1(3) forbids.
//!
//! The certificate is explicit. Take `B` a **minimum** partition and `A` the
//! all-unit-cells partition. Then `|A| > |B|`, and **no** `a ∈ A \ B` has
//! `{a} ∪ B ∈ I`, because `B` covers every cell so every rectangle in `X`
//! overlaps it. `disjoint_from_small` is that count, computed rather than argued,
//! and it must be `0`.
//!
//! The four hand-built regions, drawn. `#` is a cell of the region, `.` a hole;
//! the CSV writes the same art with `/` between rows because a value may not
//! contain a newline.
//!
//! ```text
//!   square_2x2        plus            staircase          ring
//!     # #            . # .             # . .           # # #
//!     # #            # # #             # # .           # . #
//!                    . # .             # # #           # # #
//!
//!   4 cells         5 cells           6 cells          8 cells
//!   min basis 1     min basis 3       min basis 3      min basis 4
//!   max basis 4     max basis 5       max basis 6      max basis 8
//!   |X| = 9         |X| = 11          |X| = 15         |X| = 20
//! ```
//!
//! `square_2x2` is the **minimal** counterexample: one `2×2` quad and four `1×1`
//! quads are both bases, `4 > 1`, and the single `2×2` blocks every one of the
//! nine rectangles in `X`. Nothing smaller than a `2×2` block can fail, because a
//! region of one cell has exactly one partition. `ring` carries a **hole**, which
//! is the case `R-165`'s formula models with its `h` term; its `|X| = 20` exceeds
//! the exhaustive gate so it is certified by witness rather than by scan, and the
//! CSV says which via `axioms`.
//!
//! **The shipped mesher is the fifth witness, and it needs no fixture.**
//! `Merge::Off` emits one quad per visible cell face — the all-unit-cells basis.
//! `Merge::Greedy` emits the merged basis. Both are maximal in the same system,
//! over the same occupancy, in the same binary (`greedy_quads.rs:98-104` says the
//! switch exists precisely so there is one occupancy and not two). So
//! `quads_merge_off > quads_merge_greedy` **is** a measured violation of
//! equicardinality — and it is `M-56`'s number re-read as an algebraic fact
//! rather than as a saving.
//!
//! # C2 — the LOD chunk budget
//!
//! **The model, and how it maps onto `isomesh::chunk` and `isomesh::lod`.** In
//! this crate LOD comes from the field: `ChunkLayout::at_lod(k)` keeps the chunk's
//! cell count and doubles the spacing (`chunk.rs:137-144`), and a level-`k` sample
//! position is **bit-identical** to the level-0 position of the sample `2^k` times
//! its index (`chunk.rs:152-157`). So *"chunk `c` is rendered at level `L`"* is
//! well defined as *"the field over `c`'s region is reconstructed from samples on
//! the level-`L` lattice"*, and `lod::downsample(.., Downsample::Decimate)`
//! produces exactly those samples — its kernel is `[0, 1, 0]` (`lod.rs:74`), so
//! decimation of a sample-centred nested grid **is** re-evaluation on the coarse
//! lattice. That equality is asserted bit-exactly on the first chunk measured,
//! along with `at_lod(L).cell_size() == h · 2^L`, because the whole level model
//! rests on it.
//!
//! With `CHUNK_CELLS = 8` a chunk carries `9³` samples and the ladder is
//! `9 → 5 → 3 → 2`: **four levels**, spacings `h, 2h, 4h, 8h`, and **three**
//! refinement steps per chunk — comfortably past the registered *"at least two
//! LOD levels"*, and enough for **two** second differences per chunk, which is
//! the smallest number for which *"diminishing"* means anything.
//!
//! *Named simplification.* Each level-0 chunk picks its own level independently.
//! A real quadtree cut would force eight siblings to share one. The relaxation is
//! stated rather than hidden; it is the per-chunk LOD selection a distance-driven
//! streamer performs, and it is what makes the objective a sum over a partition
//! of the ground set, which is what Theorem 4.2 needs.
//!
//! **The error metric, measured and not modelled.** For chunk `c` at level `L`:
//!
//! ```text
//! residual(c, L) = max over all 9³ fine samples p of |trilinear_L(p) − f(p)|
//! lipschitz(c)   = max |Δf| / h over every adjacent fine sample pair, all 3 axes
//! geometric(c,L) = residual(c, L) / lipschitz(c)                  [world units]
//! screen(c, L)   = geometric(c, L) · 1080 / (2 · d(c) · tan 30°)   [pixels]
//! ```
//!
//! `d(c)` is the distance from the camera to the chunk's axis-aligned box, floored
//! at one fine cell so a camera inside a chunk cannot divide by zero. The camera
//! sits at `1.5 ×` the domain's upper corner, so `d` spans a wide range and the
//! projection factor is not a constant across the world.
//!
//! `geometric(c, 0) = 0` exactly — level 0 *is* the fine grid — so the three gains
//! telescope to `screen(c, 3)` and `f` is normalised by construction.
//!
//! Dividing by the measured local Lipschitz constant converts a value residual to
//! a **length**; it is the standard proxy and is not a certified Hausdorff bound,
//! which is said here rather than implied. **It cannot change C2's verdict.** The
//! projection factor and `1/lipschitz(c)` are both per-chunk positive constants,
//! and monotonicity of a chunk's gain sequence is invariant under multiplying that
//! chunk's whole ladder by a positive constant. They move magnitudes and they move
//! *which* chunk greedy refines first; they cannot move the sign of a second
//! difference. So the camera and the normalisation are load-bearing for
//! `greedy_ratio` and inert for `marginal_returns_monotone`.
//!
//! **The set function.** Ground set `X = {(c, j) : c a surface chunk,
//! j ∈ 1..=3}`, one refinement token per chunk per level, `|X| = 3 · chunks`.
//!
//! ```text
//! f(S) = Σ over chunks c of  value(c, |S ∩ X_c|)
//! value(c, t) = screen(c, 3) − screen(c, 3 − t)      value(c, 0) = 0
//! ```
//!
//! `f` depends on the *count* of tokens spent in a chunk and not on which, so it
//! is defined on all of `2^X` and the constraint is a bare `|S| ≤ k` — Example
//! 3.5's uniform matroid, which is why `is_matroid` is **true** on these rows and
//! **false** on the meshing rows. That contrast is the registration's two halves
//! read straight off one column.
//!
//! `f` is a sum, over a partition of `X`, of functions of a count. It is
//! therefore submodular **iff** every chunk's gain sequence is non-increasing,
//! and that is the same statement as *"marginal returns are diminishing"*. The
//! harness does not take the equivalence on trust: it measures the two ends
//! separately and asserts they meet.
//!
//! **Two instruments, and only one of them can be trusted to fire.** Definition
//! 2.1 is sampled on random subset pairs, and a random pair is a **whole-world**
//! question: `f(A) + f(B) − f(A∪B) − f(A∩B)` sums one slack per chunk, so a
//! minority of chunks with increasing returns can be **hidden** by the concave
//! majority's positive slack and the scan reports nothing. That is expected
//! rather than surprising, and it is exactly why the assertion is not on the
//! scan. It is on a **constructed** extremal witness: for the chunk with the
//! largest `gain[j+1] − gain[j]`, the pair `A = {tokens 1..=j+1}` and
//! `B = {tokens 1..=j} ∪ {token j+2}` both have `j+1` elements, `A∪B` has `j+2`
//! and `A∩B` has `j`, so (2.1) there collapses to
//! `2·value(j+1) ≥ value(j+2) + value(j)` — which is Proposition 2.3's
//! `gain[j] ≥ gain[j+1]` and nothing else. The constructed pair is deterministic
//! and cannot be missed. Both readings are in the CSV as
//! `submodular_scan_saw_it` and `witness_sees_violation`, so the file settles
//! whether they agreed rather than leaving it to this paragraph.
//!
//! **The mechanism to expect, named before the numbers.** For a smooth field the
//! trilinear interpolation error at spacing `s` goes as `s²`, and `s` doubles per
//! level, so `geometric(c, L) ≈ C·4^L` and the three gains fall as
//! `48C : 12C : 4C` — strictly diminishing, by a factor of four per step. What
//! breaks that is **saturation**: at the coarsest levels a chunk carrying an
//! oscillation shorter than `8h` is undersampled, the max residual stops growing
//! and pins near the field's own amplitude, and `screen(c,3) − screen(c,2)`
//! collapses while `screen(c,2) − screen(c,1)` does not. A chunk in that regime
//! has **increasing** marginal returns, which is C2's registered falsifier. So
//! the fields are chosen to straddle it on purpose: `sphere` is smooth at every
//! level and should be clean, `capped_gyroid` oscillates on a fixed period and
//! `fbm_terrain` carries four octaves, and both should alias at `8h`. A
//! falsification here would therefore be **a statement about which fields admit
//! the guarantee**, not a statement that no field does — and
//! `violating_chunks` is recorded per field so the entry can say which.
//!
//! **The optimum is exact, not a bound.** `f(S*)` under `|S| ≤ k` is a bounded
//! knapsack over the per-chunk `value(c, ·)` tables — one `O(chunks · k · 4)` DP
//! pass per field yields the optimum for **every** budget at once. So
//! `greedy_ratio` is `f(greedy) / f(optimal)` with no relaxation in it, and
//! `greedy_ratio == 1` is a proof that greedy was optimal rather than evidence
//! that it was close.
//!
//! # Arms
//!
//! | arm | what it varies | `is_control` |
//! |---|---|---|
//! | `greedy_meshing/square_2x2` · `plus` · `staircase` | hand-built region, axioms scanned exhaustively over `2^|X|` | no |
//! | `greedy_meshing/ring` | a region **with a hole**; `|X| = 20` so certified by witness | no |
//! | `greedy_meshing/shipped@box_exact_33` · `@sphere_33` | the shipped `GreedyQuads`, `Merge::Off` against `Merge::Greedy` | no |
//! | `control/uniform_matroid_U3_8` | Example 3.5's `U(3,8)`; the tester **must** report the exchange property **holding** | **yes** |
//! | `control/knapsack_B2_w211` | Section 3.3's published failure, verbatim; the tester must reproduce its witness | **yes** |
//! | `lod_budget/sphere@k=…` | `12³` chunks over `Sphere`, six budgets | no |
//! | `lod_budget/capped_gyroid@k=…` | `7³` chunks over `CappedGyroid`, six budgets | no |
//! | `lod_budget/fbm_terrain@k=…` | `10³` chunks over `FbmTerrain`, six budgets | no |
//! | `lod_budget/<field>@worst_chunk` | the single chunk with the largest `gain[j+1] − gain[j]` | no |
//! | `control/lod_convex_synthetic` | two chunks, gains `(1,1,10)` and `(2,0,0)`; greedy **must** miss | **yes** |
//! | `control/lod_reversed_<field>` | the measured ladders with each chunk's gains reversed | **yes** |
//!
//! # The nine registered columns, and what each means on each arm
//!
//! `c1_holds` and `c2_holds` are **global** clause verdicts and carry the same
//! value on every row; the per-arm reading is in the extras as `arm_verdict`.
//! Every other column has an honest per-arm meaning and none of them is filler:
//!
//! | column | on a `greedy_meshing` / `control` set-system row | on a `lod_budget` row |
//! |---|---|---|
//! | `problem` | the arm | the arm and its budget |
//! | `is_matroid` | Definition 3.1(3), measured | **true** — `|S| ≤ k` is Example 3.5's `U(k,n)`, and the `U(3,8)` control row is where that is verified |
//! | `is_submodular` | **true** — the objective is `|S|`, which is *modular*, hence submodular. Submodularity is **not** what greedy meshing is missing | Definition 2.1 on random subset pairs **and** on the constructed extremal witness |
//! | `marginal_returns_monotone` | **true** — a modular objective has *constant* marginals, and constant is non-increasing | Proposition 2.3 per chunk: `gain[0] ≥ gain[1] ≥ gain[2]` |
//! | `greedy_ratio` | `min basis found / greedy's basis`, in `(0,1]`, `1` = greedy optimal | `f(greedy) / f(optimal)`, in `(0,1]`, `1` = greedy optimal |
//! | `bound_applies` | whether Theorem 4.2's hypotheses are met | `normalized ∧ monotone ∧ submodular ∧ cardinality-constrained` |
//! | `violating_chunks` | `0`, **and it is not a measurement here** — this arm has no chunks. `arm_has_chunks=false` says so on the row | chunks with `gain[j+1] > gain[j]` for some `j` |
//!
//! `greedy_ratio` is one definition in both directions: **what greedy achieved
//! divided by the best achievable**, so it lands in `(0,1]` and `1` always means
//! optimal — for a minimisation that is `optimum / greedy`, for a maximisation
//! `greedy / optimum`.
//!
//! `bound_applies` is **false** on every meshing row and the reason is structural,
//! not a technicality: Theorem 4.2 needs a monotone submodular objective
//! *maximised* under a cardinality constraint. Greedy meshing *minimises*
//! cardinality subject to *covering*, and its constraint system fails Definition
//! 3.1(3). Two of the theorem's hypotheses fail, so the `(1 − 1/e)` guarantee is
//! not weak here — it is **inapplicable**, which is what C1 was registered to say.
//!
//! **The extras are sparse and that is deliberate.** The two halves measure
//! different objects, so `fixture_art` and `min_basis_rects` are blank on a
//! `lod_budget` row and `curve_chunks` and `f_optimal` are blank on a
//! `greedy_meshing` one. The alternative is a column whose meaning changes with
//! the row, which is the defect `M-273` and `P-64` were both about; an empty cell
//! is honest where a re-used one is not. `arm_kind` names which group a row
//! belongs to, so the file can be split on one column.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"C2 governs the LOD selection stage, whose
//! amortised cost `M-124` tracks to within one chunk."* Discharged, and it is
//! **not** zero.
//!
//! `M-124`/`M-125` already fix the denominator in *time*: the queue spends a
//! budget and lands within one chunk (≈0.072 ms) of it. What C2 governs is not
//! how long that budget takes but **which chunks it is spent on**, so the share is
//! a fraction of *error* at fixed cost, not a fraction of a frame. The harness
//! measures it directly by running a third ordering over the same ladders:
//!
//! - `f_greedy` — refine whichever chunk currently has the largest marginal.
//! - `f_distance` — refine the **nearest** chunk all the way down, then the next:
//!   the distance-only heuristic a streamer's load radius already implies.
//! - `f_opt` — the exact knapsack optimum.
//!
//! `share_over_distance_order = (f_greedy − f_distance) / f_opt` is the pixels of
//! error the *same* chunk budget buys by ordering on marginal gain instead of on
//! distance, as a fraction of the most any budget of that size could buy. It is
//! recorded on every `lod_budget` row and it is the number the SHARE sentence
//! should be quoted from.
//!
//! # Vacuity controls
//!
//! `M-44`: a zero that could not have been non-zero is not a measurement. Every
//! control below is `assert!`ed with a `VOID: ` message **and** written to the
//! CSV, so a run that cannot fire aborts instead of recording a pass.
//!
//! | zero or verdict at risk | control, asserted | why it licenses the reading |
//! |---|---|---|
//! | `is_matroid = false` on the meshing rows | `control/uniform_matroid_U3_8` reports the augmentation property **holding** | Example 3.5 *proves* `U(3,8)` is a matroid, so a failure there would be a tester bug and every meshing `false` would be worthless |
//! | the same, from the other side | `control/knapsack_B2_w211` reports **failure**, with `{b,c}` against `{a}` | Section 3.3 publishes that exact witness; the tester reproducing it is the tester being right about a case it did not choose |
//! | `disjoint_from_small = 0` | `max_basis > min_basis` on every fixture | a region with one partition has nothing to exchange, and `0` there would mean *"no larger basis exists"* rather than *"no augmentation is possible"* |
//! | the shipped witness | `quads_merge_off > quads_merge_greedy` | if merging emitted the same count, `Merge` would be inert and the shipped arm would witness nothing |
//! | `violating_chunks = 0` | `control/lod_convex_synthetic` reports `violating_chunks = 1` and `greedy_ratio = 1/3` | the same second-difference test and the same DP, on a ladder built to have increasing returns. `1/3 < 1 − 1/e`, so the instrument is shown able to report *both* an increasing return **and** a `greedy_ratio` that breaks the bound |
//! | the same, on real data | `control/lod_reversed_<field>` reports `violating_chunks > 0` | the measured ladders themselves, reversed — so the detector is exercised on the real magnitudes and not only on hand-picked ones |
//! | `curve_chunks` | `≥ 50` on **every** measured field, asserted per field | the registration's own bar, verbatim |
//! | `levels_spanned` | `≥ 3` levels and `≥ 2` second differences at compile time, and `screen[3] > screen[1] > 0` on some chunk at run time | *"at least two LOD levels"* has to mean the levels **differ**, not merely that four were named |
//! | the two submodularity instruments | Definition 2.1 at the **constructed** extremal witness must agree with Proposition 2.3's per-chunk verdict, asserted whenever a violating chunk exists | Proposition 2.3 says they are the same statement; disagreement means one implementation is wrong. The *random* scan is **not** what this is asserted on, because a random pair over a sum can be blind to a minority of convex chunks — both readings are recorded as `submodular_scan_saw_it` and `witness_sees_violation` so a disagreement is in the file rather than in a footnote |
//!
//! # Determinism
//!
//! One thread, no map iteration, `f64` throughout. Every sweep order is fixed:
//! chunks in `z`-major `ChunkId` order, fixtures in source order, budgets
//! ascending. The only randomness is the Definition 2.1 pair sampler, which draws
//! from `common::poly::Rng` — SplitMix64 — seeded at `0x5EED_0166`, stated here
//! and recorded as `rng_seed`. Every ordering uses [`f64::total_cmp`] with an
//! index tiebreak, never `partial_cmp`, so a NaN sorts into view rather than being
//! dropped. No clause here is a wall clock, so `M-280`'s 1.45× governor swing
//! cannot reach any verdict; `wall_ns` is recorded because it is interesting and
//! is read by nothing.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::needless_range_loop,
    clippy::too_many_lines
)]

mod common;

use std::cmp::Ordering;
use std::time::Instant;

use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::{BoxExact, FbmTerrain, ReferenceField, Sphere, capped_gyroid};
use isomesh::greedy_quads::{GreedyQuads, Merge};
use isomesh::lod::{Downsample, downsample};
use isomesh::{MeshBuffer, Sdf, Shape3};

use common::poly::Rng;

// ─── constants ──────────────────────────────────────────────────────────────

/// Levels in a chunk's ladder: `0` is the fine grid, `3` the coarsest.
const LEVELS: usize = 4;

/// The coarsest level, where every chunk starts.
const TOP: usize = LEVELS - 1;

/// Refinement steps per chunk, so `STEPS - 1` second differences.
const STEPS: usize = LEVELS - 1;

/// The ladder's shape is a property of this file rather than of the data, so it
/// is checked where it can be: at compile time. *"Diminishing"* is not a
/// statement about fewer than two marginal gains, and the registration asks for
/// at least two LOD levels — a ladder that cannot form a second difference
/// should not build. The **measured** span is a separate runtime control, because
/// four levels that all read the same error would satisfy this and still be
/// vacuous.
const _LADDER_IS_LONG_ENOUGH: () = assert!(
    LEVELS >= 3 && STEPS >= 2,
    "P-166: the ladder must carry at least two refinement steps"
);

/// Cells per chunk axis. `9` samples, halving to `5`, `3`, `2`.
const CHUNK_CELLS: u32 = 8;

/// The registration's bar: *"at least 50 chunks"*.
const MIN_CURVE_CHUNKS: usize = 50;

/// Vertical resolution the screen-space error is quoted in.
const SCREEN_HEIGHT_PX: f64 = 1080.0;

/// `tan(30°)` — half of a 60° vertical field of view. A literal rather than a
/// call so the projection factor is identical on every machine and every libm.
const TAN_HALF_FOV: f64 = 0.577_350_269_189_625_7;

/// The camera sits at this multiple of the domain's upper corner.
const CAMERA_STANDOFF: f64 = 1.5;

/// `1 − 1/e`, Theorem 4.2's constant.
const NWF_BOUND: f64 = 0.632_120_558_828_557_7;

/// Random subset pairs per field for Definition 2.1.
const SUBMODULAR_PAIRS: usize = 20_000;

/// Relative slack allowed in Definition 2.1: these are sums over hundreds of
/// chunks, so a bit-level cancellation is not a violation of submodularity.
const SUBMODULAR_REL_TOL: f64 = 1e-9;

/// Stated in the header and recorded as a column.
const RNG_SEED: u64 = 0x5EED_0166;

/// Largest ground set the axiom scan will enumerate: `2^16` subsets.
const MAX_EXHAUSTIVE_GROUND: usize = 16;

/// Guard on the exhaustive scan's pair loop.
const MAX_INDEPENDENT_SETS: usize = 4096;

/// Guard on the partition enumeration.
const MAX_PARTITIONS: usize = 200_000;

/// Samples per axis for the shipped `GreedyQuads` witness.
const SHIPPED_SAMPLES: u32 = 33;

/// The hand-built regions. `#` is a cell of the region, `.` a hole, `/` a row
/// break — `/` rather than a newline because a CSV value may not contain one.
const FIXTURES: [(&str, &str); 4] = [
    ("square_2x2", "##/##"),
    ("plus", ".#./###/.#."),
    ("staircase", "#../##./###"),
    ("ring", "###/#.#/###"),
];

// ─── C1: rectilinear regions and their packing system ───────────────────────

/// An axis-aligned rectangle of cells, minimum corner `(x, y)`.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Rect {
    /// Minimum corner along the mask's first axis.
    x: usize,
    /// Minimum corner along the mask's second axis.
    y: usize,
    /// Extent along the first axis, in cells. Always at least one.
    w: usize,
    /// Extent along the second axis, in cells. Always at least one.
    h: usize,
}

impl Rect {
    /// CSV-safe rendering: `2x3@0+1` is a `2×3` rectangle at `(0, 1)`.
    fn tag(self) -> String {
        format!("{}x{}@{}+{}", self.w, self.h, self.x, self.y)
    }
}

/// A rectilinear region, as a rectangular mask of present and absent cells.
struct Mask {
    /// Extent along the first axis.
    w: usize,
    /// Extent along the second axis.
    h: usize,
    /// Row-major occupancy, first axis fastest.
    on: Vec<bool>,
}

impl Mask {
    /// Parse the `#`/`.` art used by [`FIXTURES`].
    fn parse(art: &str) -> Self {
        let rows: Vec<&str> = art.split('/').collect();
        let h = rows.len();
        let w = rows[0].len();
        assert!(h > 0 && w > 0, "P-166: an empty fixture is not a region");
        assert!(
            rows.iter().all(|r| r.len() == w),
            "P-166: fixture rows must all be {w} wide: {art}"
        );
        let mut on = Vec::with_capacity(w * h);
        for row in &rows {
            for ch in row.chars() {
                assert!(
                    ch == '#' || ch == '.',
                    "P-166: fixture art takes # and . only"
                );
                on.push(ch == '#');
            }
        }
        Self { w, h, on }
    }

    /// Whether `(x, y)` is in the region.
    fn at(&self, x: usize, y: usize) -> bool {
        self.on[x + self.w * y]
    }

    /// Cells in the region, which is also the size of its largest basis.
    fn cells(&self) -> usize {
        self.on.iter().filter(|b| **b).count()
    }

    /// Every axis-aligned rectangle lying entirely inside the region: the ground
    /// set `X` of Definition 3.1.
    fn rectangles(&self) -> Vec<Rect> {
        let mut out = Vec::new();
        for y in 0..self.h {
            for x in 0..self.w {
                if !self.at(x, y) {
                    continue;
                }
                for h in 1..=(self.h - y) {
                    for w in 1..=(self.w - x) {
                        let mut full = true;
                        for dy in 0..h {
                            for dx in 0..w {
                                if !self.at(x + dx, y + dy) {
                                    full = false;
                                }
                            }
                        }
                        if full {
                            out.push(Rect { x, y, w, h });
                        }
                    }
                }
            }
        }
        out
    }

    /// One bit per cell of the mask, so rectangle disjointness is a bit test.
    fn cell_bits(&self, r: Rect) -> u128 {
        let mut bits = 0u128;
        for dy in 0..r.h {
            for dx in 0..r.w {
                bits |= 1u128 << (r.x + dx + self.w * (r.y + dy));
            }
        }
        bits
    }

    /// Every cell of the region, as one bit mask.
    fn full_bits(&self) -> u128 {
        let mut bits = 0u128;
        for i in 0..self.on.len() {
            if self.on[i] {
                bits |= 1u128 << i;
            }
        }
        bits
    }

    /// The shipped greedy merge, run on one mask.
    ///
    /// Transcribed from `crates/isomesh/src/greedy_quads.rs:238-272` with `du`
    /// and `dv` renamed to this mask's `w` and `h`: widest run along the first
    /// axis, then as many full-width rows along the second as continue it. The
    /// source is `pub` only as a whole extractor, so the slice loop is copied
    /// with the line it came from rather than reached for.
    fn greedy_partition(&self) -> Vec<Rect> {
        let mut m = self.on.clone();
        let mut out = Vec::new();
        for b in 0..self.h {
            let mut a = 0usize;
            while a < self.w {
                if !m[a + self.w * b] {
                    a += 1;
                    continue;
                }
                let mut width = 1usize;
                let mut height = 1usize;
                while a + width < self.w && m[a + width + self.w * b] {
                    width += 1;
                }
                'grow: while b + height < self.h {
                    for k in 0..width {
                        if !m[a + k + self.w * (b + height)] {
                            break 'grow;
                        }
                    }
                    height += 1;
                }
                for row in 0..height {
                    for k in 0..width {
                        m[a + k + self.w * (b + row)] = false;
                    }
                }
                out.push(Rect {
                    x: a,
                    y: b,
                    w: width,
                    h: height,
                });
                a += width;
            }
        }
        out
    }
}

/// Every rectangle partition of a region, which is every basis of its packing
/// system.
struct Partitions {
    /// How many were found.
    count: usize,
    /// Whether the enumeration finished inside [`MAX_PARTITIONS`].
    complete: bool,
    /// Cardinality of the smallest.
    min: usize,
    /// Cardinality of the largest.
    max: usize,
    /// One partition of minimum cardinality: `B` in the exchange certificate.
    small: Vec<Rect>,
}

/// Enumerate every rectangle partition, each exactly once.
///
/// The recursion assigns the **lexicographically first uncovered cell** to a
/// rectangle whose minimum corner is that cell — which is forced, because a
/// rectangle containing that cell and starting earlier would contain a cell that
/// is already covered. So each partition is reached along exactly one path and
/// no deduplication structure is needed.
fn enumerate_partitions(mask: &Mask, rects: &[Rect]) -> Partitions {
    let full = mask.full_bits();
    let bits: Vec<u128> = rects.iter().map(|r| mask.cell_bits(*r)).collect();
    let mut found = Partitions {
        count: 0,
        complete: true,
        min: usize::MAX,
        max: 0,
        small: Vec::new(),
    };
    let mut chosen: Vec<Rect> = Vec::new();
    walk(mask, rects, &bits, full, 0, &mut chosen, &mut found);
    found
}

/// One node of [`enumerate_partitions`]'s recursion.
fn walk(
    mask: &Mask,
    rects: &[Rect],
    bits: &[u128],
    full: u128,
    covered: u128,
    chosen: &mut Vec<Rect>,
    found: &mut Partitions,
) {
    if covered == full {
        found.count += 1;
        if chosen.len() < found.min {
            found.min = chosen.len();
            found.small = chosen.clone();
        }
        if chosen.len() > found.max {
            found.max = chosen.len();
        }
        if found.count >= MAX_PARTITIONS {
            found.complete = false;
        }
        return;
    }
    if !found.complete {
        return;
    }
    let first = (full & !covered).trailing_zeros() as usize;
    let (x0, y0) = (first % mask.w, first / mask.w);
    for (i, r) in rects.iter().enumerate() {
        if r.x != x0 || r.y != y0 || bits[i] & covered != 0 {
            continue;
        }
        chosen.push(*r);
        walk(mask, rects, bits, full, covered | bits[i], chosen, found);
        chosen.pop();
        if !found.complete {
            return;
        }
    }
}

/// What Definition 3.1's three axioms report on one set system.
struct Axioms {
    /// Independent sets found, including the empty one.
    independent_sets: usize,
    /// Axiom (1).
    empty_independent: bool,
    /// Axiom (2).
    downward_closed: bool,
    /// Axiom (3), the augmentation — and therefore exchange — property.
    augmentation_holds: bool,
    /// Ordered pairs `(A, B)` with `|A| > |B|` admitting no augmentation.
    augmentation_failures: u64,
    /// The first such pair, as element-index lists.
    witness: Option<(u32, u32)>,
}

/// Check Definition 3.1(1)(2)(3) exhaustively over all `2^n` subsets.
///
/// Exhaustive on purpose. The augmentation property is a statement about *every*
/// ordered pair of independent sets, and sampling pairs could only ever report
/// *"no counterexample found"* — which is not what C1 claims in either direction.
fn check_axioms<F>(n: usize, independent: F) -> Axioms
where
    F: Fn(u32) -> bool,
{
    assert!(
        n <= MAX_EXHAUSTIVE_GROUND,
        "P-166: a ground set of {n} is past the exhaustive gate of {MAX_EXHAUSTIVE_GROUND}"
    );
    let mut sets: Vec<u32> = Vec::new();
    for m in 0..(1u32 << n) {
        if independent(m) {
            sets.push(m);
        }
    }
    assert!(
        sets.len() <= MAX_INDEPENDENT_SETS,
        "P-166: {} independent sets is past the pair-loop guard of {MAX_INDEPENDENT_SETS}",
        sets.len()
    );

    let mut downward_closed = true;
    for &m in &sets {
        let mut rest = m;
        while rest != 0 {
            // `rest & rest.wrapping_neg()`, not `isolate_lowest_one()`: the
            // intrinsic is stable only since 1.97 and this crate's MSRV is 1.89.
            let bit = rest & rest.wrapping_neg();
            rest ^= bit;
            if !independent(m ^ bit) {
                downward_closed = false;
            }
        }
    }

    let mut failures = 0u64;
    let mut witness = None;
    for &a in &sets {
        for &b in &sets {
            if a.count_ones() <= b.count_ones() {
                continue;
            }
            let mut augmented = false;
            let mut extra = a & !b;
            while extra != 0 {
                let bit = extra & extra.wrapping_neg();
                extra ^= bit;
                if independent(b | bit) {
                    augmented = true;
                    break;
                }
            }
            if !augmented {
                failures += 1;
                if witness.is_none() {
                    witness = Some((a, b));
                }
            }
        }
    }

    Axioms {
        independent_sets: sets.len(),
        empty_independent: independent(0),
        downward_closed,
        augmentation_holds: failures == 0,
        augmentation_failures: failures,
        witness,
    }
}

/// A subset mask as a `|`-joined list of element indices, CSV-safe.
fn elements(mask: u32) -> String {
    if mask == 0 {
        return String::from("empty");
    }
    let mut parts: Vec<String> = Vec::new();
    for i in 0..32 {
        if mask >> i & 1 == 1 {
            parts.push(i.to_string());
        }
    }
    parts.join("|")
}

/// Everything one hand-built region contributes to the CSV.
struct Fixture {
    /// The arm's name.
    name: &'static str,
    /// The art, with `/` between rows.
    art: &'static str,
    /// Cells in the region.
    cells: usize,
    /// `|X|`.
    ground: usize,
    /// Every basis, or the certificate's endpoints if enumeration was capped.
    parts: Partitions,
    /// `Merge::Greedy`'s cardinality on this region.
    greedy: usize,
    /// Rectangles of `X` disjoint from the minimum basis. Must be zero.
    disjoint_from_small: usize,
    /// Definition 3.1's scan, where `|X|` allowed one.
    axioms: Option<Axioms>,
}

/// Measure one hand-built region.
fn measure_fixture(name: &'static str, art: &'static str) -> Fixture {
    let mask = Mask::parse(art);
    let rects = mask.rectangles();
    let parts = enumerate_partitions(&mask, &rects);
    let greedy = mask.greedy_partition().len();

    // Every rectangle of `X` overlaps the minimum basis, because that basis
    // covers every cell. Computed rather than argued: this is the half of the
    // certificate that says no augmentation is *possible*, not merely that a
    // larger basis exists.
    let mut small_bits = 0u128;
    for r in &parts.small {
        small_bits |= mask.cell_bits(*r);
    }
    let disjoint_from_small = rects
        .iter()
        .filter(|r| mask.cell_bits(**r) & small_bits == 0)
        .count();

    let axioms = if rects.len() <= MAX_EXHAUSTIVE_GROUND {
        let bits: Vec<u128> = rects.iter().map(|r| mask.cell_bits(*r)).collect();
        Some(check_axioms(rects.len(), |m| {
            let mut used = 0u128;
            let mut rest = m;
            while rest != 0 {
                let i = rest.trailing_zeros() as usize;
                rest &= rest - 1;
                if used & bits[i] != 0 {
                    return false;
                }
                used |= bits[i];
            }
            true
        }))
    } else {
        None
    };

    Fixture {
        name,
        art,
        cells: mask.cells(),
        ground: rects.len(),
        parts,
        greedy,
        disjoint_from_small,
        axioms,
    }
}

/// Quads the shipped mesher emits under one merge setting.
fn shipped_quads<F>(field: &F, merge: Merge) -> usize
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, cell) = common::grid::<f64, _>(field, SHIPPED_SAMPLES);
    let mut mesher = GreedyQuads::<f64>::new();
    mesher.set_merge(merge);
    let mut out = MeshBuffer::<f64>::new();
    mesher
        .extract(field, &shape, origin, cell, &mut out)
        .expect("greedy quads extracts a reference field at 33 samples");
    out.triangle_count() / 2
}

// ─── C2: the LOD chunk budget ───────────────────────────────────────────────

/// One level-0 chunk's error ladder, in world units and in pixels.
struct ChunkCurve {
    /// Chunk coordinates, for the witness row.
    coords: [i32; 3],
    /// Camera-to-box distance, floored at one fine cell.
    distance: f64,
    /// `max |Δf| / h` over adjacent fine samples: the local Lipschitz constant.
    lipschitz: f64,
    /// Geometric error per level, world units. `geometric[0]` is exactly zero.
    geometric: [f64; LEVELS],
    /// Screen-space error per level, pixels.
    screen: [f64; LEVELS],
    /// `gain[j]` is the pixels refinement step `j + 1` removes.
    gain: [f64; STEPS],
}

impl ChunkCurve {
    /// `value(t)` is the pixels the first `t` refinement steps remove.
    fn value(&self) -> [f64; LEVELS] {
        let mut v = [0.0f64; LEVELS];
        for t in 1..LEVELS {
            v[t] = v[t - 1] + self.gain[t - 1];
        }
        v
    }

    /// The largest `gain[j + 1] − gain[j]`, and where. Positive means the
    /// marginal returns **increase**, which is C2's falsifier.
    fn worst_second_difference(&self) -> (f64, usize) {
        let mut worst = f64::NEG_INFINITY;
        let mut at = 0usize;
        for j in 0..STEPS - 1 {
            let d = self.gain[j + 1] - self.gain[j];
            if d.total_cmp(&worst) == Ordering::Greater {
                worst = d;
                at = j;
            }
        }
        (worst, at)
    }

    /// Whether some step's marginal exceeds the one before it.
    fn violates(&self) -> bool {
        self.worst_second_difference().0 > 0.0
    }

    /// The same, but only when the inversion is more than one part in a hundred
    /// — so a float artefact is separable from a real inversion.
    fn violates_materially(&self) -> bool {
        (0..STEPS - 1).any(|j| self.gain[j + 1] > self.gain[j] * 1.01 + f64::EPSILON)
    }
}

/// Trilinear interpolation of a coarse grid, evaluated at a **fine** sample
/// index.
///
/// `step` is the fine samples per coarse cell edge, `2^L`. At the far face the
/// fine index lands exactly on the last coarse sample, so the cell is the last
/// one and the parameter is `1` — rather than an out-of-range base index.
fn trilinear(coarse: &[f64], n: usize, step: usize, index: [usize; 3]) -> f64 {
    let mut base = [0usize; 3];
    let mut t = [0.0f64; 3];
    for a in 0..3 {
        let q = index[a] / step;
        if q + 1 < n {
            base[a] = q;
            t[a] = (index[a] % step) as f64 / step as f64;
        } else {
            base[a] = n - 2;
            t[a] = 1.0;
        }
    }
    let at = |dx: usize, dy: usize, dz: usize| -> f64 {
        coarse[(base[0] + dx) + n * ((base[1] + dy) + n * (base[2] + dz))]
    };
    let lerp = |a: f64, b: f64, s: f64| a + (b - a) * s;
    let y0 = lerp(
        lerp(at(0, 0, 0), at(1, 0, 0), t[0]),
        lerp(at(0, 1, 0), at(1, 1, 0), t[0]),
        t[1],
    );
    let y1 = lerp(
        lerp(at(0, 0, 1), at(1, 0, 1), t[0]),
        lerp(at(0, 1, 1), at(1, 1, 1), t[0]),
        t[1],
    );
    lerp(y0, y1, t[2])
}

/// The largest `|trilinear_L(p) − f(p)|` over every fine sample of a chunk.
fn residual(coarse: &[f64], n: usize, step: usize, fine: &[f64], fine_n: usize) -> f64 {
    let mut worst = 0.0f64;
    for k in 0..fine_n {
        for j in 0..fine_n {
            for i in 0..fine_n {
                let got = trilinear(coarse, n, step, [i, j, k]);
                let want = fine[i + fine_n * (j + fine_n * k)];
                let d = (got - want).abs();
                if d > worst {
                    worst = d;
                }
            }
        }
    }
    worst
}

/// One field's marginal-return curve and the metadata the CSV needs.
struct FieldCurve {
    /// The field's `ReferenceField` name.
    field: &'static str,
    /// Chunks per axis.
    per_axis: i32,
    /// Chunks visited.
    chunks_total: usize,
    /// Fine cell size.
    cell_size: f64,
    /// Camera position.
    camera: [f64; 3],
    /// One entry per **surface** chunk, in `ChunkId` order.
    curve: Vec<ChunkCurve>,
}

/// Sample a chunked world, build every surface chunk's error ladder.
///
/// `checked_nesting` gates the one-off bit-exactness assertions that the whole
/// level model rests on; they run on the first surface chunk of the first field
/// and nowhere else, because they are properties of the crate rather than of the
/// data.
fn measure_field<F>(
    field: &F,
    name: &'static str,
    per_axis: i32,
    checked_nesting: &mut bool,
) -> FieldCurve
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(per_axis * CHUNK_CELLS as i32);
    let layout = ChunkLayout::<f64>::new(CHUNK_CELLS, cell_size, lo)
        .expect("a positive cell size and a non-zero cell count make a layout");
    let shape = layout
        .sample_shape()
        .expect("nine samples per axis fit the index space");
    let fine_n = CHUNK_CELLS as usize + 1;
    assert_eq!(
        shape.element_count(),
        fine_n * fine_n * fine_n,
        "P-166: {name} — a chunk's sample grid is not {fine_n}³"
    );

    // The level model, asserted against the crate rather than assumed.
    // `at_lod(L)` doubles the spacing L times (chunk.rs:165-171), which is the
    // spacing the level-L lattice is on.
    for l in 0..LEVELS as u32 {
        let want = cell_size * f64::from(1u32 << l);
        let got = layout
            .at_lod(l)
            .expect("doubling a finite cell size three times stays finite")
            .cell_size();
        assert_eq!(
            got.to_bits(),
            want.to_bits(),
            "P-166: {name} — at_lod({l}) gives spacing {got}, not {want}; the level \
             ladder is not the crate's"
        );
    }

    let span = cell_size * f64::from(CHUNK_CELLS);
    let camera = [
        hi[0] * CAMERA_STANDOFF,
        hi[1] * CAMERA_STANDOFF,
        hi[2] * CAMERA_STANDOFF,
    ];

    let mut fine = vec![0.0f64; fine_n * fine_n * fine_n];
    let mut curve: Vec<ChunkCurve> = Vec::new();
    let mut chunks_total = 0usize;

    for cz in 0..per_axis {
        for cy in 0..per_axis {
            for cx in 0..per_axis {
                chunks_total += 1;
                let id = ChunkId::new([cx, cy, cz]);
                for k in 0..fine_n {
                    for j in 0..fine_n {
                        for i in 0..fine_n {
                            let global = layout.global_sample(id, [i as u32, j as u32, k as u32]);
                            fine[i + fine_n * (j + fine_n * k)] =
                                field.sample(layout.world_of_sample(global));
                        }
                    }
                }

                // A surface chunk is one the zero set passes through. `cube.rs:159`
                // fixes the convention: a sample of exactly zero is **outside**.
                let inside = fine[0] < 0.0;
                if fine.iter().all(|v| (*v < 0.0) == inside) {
                    continue;
                }

                // The local Lipschitz constant, measured on the grid this chunk is
                // reconstructed from. `fbm_terrain`'s declared bound is
                // `Unbounded` (fields/mod.rs:104-109), so a declared constant is
                // not available for every field and a measured one is used for all
                // three — one path, not two.
                let mut lipschitz = 0.0f64;
                for k in 0..fine_n {
                    for j in 0..fine_n {
                        for i in 0..fine_n {
                            let here = fine[i + fine_n * (j + fine_n * k)];
                            let mut probe = |a: f64| {
                                let d = (here - a).abs() / cell_size;
                                if d > lipschitz {
                                    lipschitz = d;
                                }
                            };
                            if i + 1 < fine_n {
                                probe(fine[i + 1 + fine_n * (j + fine_n * k)]);
                            }
                            if j + 1 < fine_n {
                                probe(fine[i + fine_n * (j + 1 + fine_n * k)]);
                            }
                            if k + 1 < fine_n {
                                probe(fine[i + fine_n * (j + fine_n * (k + 1))]);
                            }
                        }
                    }
                }
                if lipschitz <= 0.0 {
                    continue;
                }

                let mut geometric = [0.0f64; LEVELS];
                let mut level = fine.clone();
                let mut level_shape = shape;
                for l in 1..LEVELS {
                    let (coarse, coarse_shape) =
                        downsample(&level, &level_shape, Downsample::Decimate)
                            .expect("a 2^k + 1 grid halves");
                    let n = coarse_shape.size()[0] as usize;
                    let step = 1usize << l;
                    if !*checked_nesting {
                        // `Decimate`'s kernel is `[0, 1, 0]` (lod.rs:74) and the
                        // grids are sample-centred and nested (chunk.rs:152-157),
                        // so a coarse sample must be **the same bits** as the fine
                        // sample `2^L` times its index. If this ever failed, every
                        // level below would be measuring a filtered field rather
                        // than a coarser sampling of the same one.
                        for cq in 0..n {
                            for bq in 0..n {
                                for aq in 0..n {
                                    let got = coarse[aq + n * (bq + n * cq)];
                                    let want =
                                        fine[aq * step + fine_n * (bq * step + fine_n * cq * step)];
                                    assert_eq!(
                                        got.to_bits(),
                                        want.to_bits(),
                                        "P-166: decimation to level {l} moved a sample: \
                                         {got} against {want}"
                                    );
                                }
                            }
                        }
                    }
                    geometric[l] = residual(&coarse, n, step, &fine, fine_n) / lipschitz;
                    level = coarse;
                    level_shape = coarse_shape;
                }
                *checked_nesting = true;

                let chunk_lo = layout.sample_origin(id);
                let mut sum = 0.0f64;
                for a in 0..3 {
                    let d = (chunk_lo[a] - camera[a])
                        .max(camera[a] - (chunk_lo[a] + span))
                        .max(0.0);
                    sum += d * d;
                }
                let distance = sum.sqrt().max(cell_size);
                let factor = SCREEN_HEIGHT_PX / (2.0 * distance * TAN_HALF_FOV);

                let mut screen = [0.0f64; LEVELS];
                for l in 0..LEVELS {
                    screen[l] = geometric[l] * factor;
                }
                let mut gain = [0.0f64; STEPS];
                for j in 0..STEPS {
                    gain[j] = screen[TOP - j] - screen[TOP - j - 1];
                }

                curve.push(ChunkCurve {
                    coords: [cx, cy, cz],
                    distance,
                    lipschitz,
                    geometric,
                    screen,
                    gain,
                });
            }
        }
    }

    FieldCurve {
        field: name,
        per_axis,
        chunks_total,
        cell_size,
        camera,
        curve,
    }
}

/// The exact optimum of `max f(S)` subject to `|S| ≤ t`, for every `t ≤ kmax`.
///
/// A bounded knapsack over the per-chunk `value` tables: each chunk contributes
/// one item of size `0..=STEPS`. Reverse iteration over the budget keeps `dp` one
/// chunk behind on the right-hand side, so no second array is needed and no chunk
/// can be spent twice.
fn optimum_curve(values: &[[f64; LEVELS]], kmax: usize) -> Vec<f64> {
    let mut dp = vec![0.0f64; kmax + 1];
    for row in values {
        for t in (0..=kmax).rev() {
            let mut best = dp[t];
            for u in 1..LEVELS {
                if u <= t {
                    let cand = dp[t - u] + row[u];
                    if cand.total_cmp(&best) == Ordering::Greater {
                        best = cand;
                    }
                }
            }
            dp[t] = best;
        }
    }
    dp
}

/// Section 4's greedy: take the largest marginal still available, `k` times.
///
/// Faithful to the survey's pseudocode, which returns only when the valid set is
/// empty — so a non-positive marginal is still taken rather than skipped. That
/// matters only when `f` is not monotone, which `monotone` records.
fn greedy_value(gains: &[[f64; STEPS]], k: usize) -> f64 {
    let mut counts = vec![0usize; gains.len()];
    let mut total = 0.0f64;
    for _ in 0..k {
        let mut best: Option<(usize, f64)> = None;
        for (c, g) in gains.iter().enumerate() {
            if counts[c] >= STEPS {
                continue;
            }
            let marginal = g[counts[c]];
            let better = match best {
                None => true,
                Some((_, bm)) => marginal.total_cmp(&bm) == Ordering::Greater,
            };
            if better {
                best = Some((c, marginal));
            }
        }
        let Some((c, marginal)) = best else { break };
        counts[c] += 1;
        total += marginal;
    }
    total
}

/// The distance-only heuristic: refine the nearest chunk all the way down, then
/// the next. The SHARE's baseline.
fn distance_value(values: &[[f64; LEVELS]], order: &[usize], k: usize) -> f64 {
    let mut spent = 0usize;
    let mut total = 0.0f64;
    for &c in order {
        if spent >= k {
            break;
        }
        let take = STEPS.min(k - spent);
        total += values[c][take];
        spent += take;
    }
    total
}

/// Definition 2.1 on random subset pairs: `f(A) + f(B) ≥ f(A∪B) + f(A∩B)`.
///
/// A chunk's token subset is three bits, so the union and intersection are a
/// bitwise `|` and `&` and `f` is read off the `value` table by population count.
/// Deliberately the union/intersection form and not Proposition 2.3's marginal
/// form: the two are equivalent, so implementing both and asserting they agree
/// makes each a check on the other.
fn submodular_scan(values: &[[f64; LEVELS]], rng: &mut Rng, pairs: usize) -> (u64, f64) {
    let mut violations = 0u64;
    let mut worst = 0.0f64;
    for _ in 0..pairs {
        let (mut fa, mut fb, mut fu, mut fi) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for row in values {
            let word = rng.next_u64();
            let a = (word & 7) as u32;
            let b = (word >> 3 & 7) as u32;
            fa += row[a.count_ones() as usize];
            fb += row[b.count_ones() as usize];
            fu += row[(a | b).count_ones() as usize];
            fi += row[(a & b).count_ones() as usize];
        }
        let slack = (fa + fb) - (fu + fi);
        if slack < -SUBMODULAR_REL_TOL * (fa + fb).abs().max(1.0) {
            violations += 1;
        }
        if slack < worst {
            worst = slack;
        }
    }
    (violations, worst)
}

/// Definition 2.1 at the **constructed** extremal witness for a violating chunk.
///
/// With `A` the first `j+1` tokens of the chunk and `B` its first `j` tokens plus
/// token `j+2`, both have `j+1` elements, `A∪B` has `j+2` and `A∩B` has `j`. So
/// (2.1) reduces to `2·value(j+1) ≥ value(j+2) + value(j)`, which is exactly
/// Proposition 2.3's `gain[j] ≥ gain[j+1]`. Returns whether (2.1) is violated
/// there — deterministic, where the random scan is not.
fn witness_violates(row: &[f64; LEVELS], j: usize) -> bool {
    2.0 * row[j + 1] < row[j + 2] + row[j]
}

/// A gain ladder reversed step for step: the reversed-real control.
fn reverse_gains(curve: &[ChunkCurve]) -> Vec<[f64; STEPS]> {
    curve
        .iter()
        .map(|c| {
            let mut g = [0.0f64; STEPS];
            for j in 0..STEPS {
                g[j] = c.gain[STEPS - 1 - j];
            }
            g
        })
        .collect()
}

/// Prefix sums of a gain ladder.
fn values_of(gains: &[[f64; STEPS]]) -> Vec<[f64; LEVELS]> {
    gains
        .iter()
        .map(|g| {
            let mut v = [0.0f64; LEVELS];
            for t in 1..LEVELS {
                v[t] = v[t - 1] + g[t - 1];
            }
            v
        })
        .collect()
}

/// Chunks whose marginal returns increase, strictly and materially.
fn count_violations(gains: &[[f64; STEPS]]) -> (u64, u64) {
    let mut strict = 0u64;
    let mut material = 0u64;
    for g in gains {
        if (0..STEPS - 1).any(|j| g[j + 1] > g[j]) {
            strict += 1;
        }
        if (0..STEPS - 1).any(|j| g[j + 1] > g[j] * 1.01 + f64::EPSILON) {
            material += 1;
        }
    }
    (strict, material)
}

// ─── rows ───────────────────────────────────────────────────────────────────

/// One CSV row's nine registered columns, so every arm emits the same shape.
struct Row {
    /// The arm.
    problem: String,
    /// Definition 3.1(3) for that arm's constraint system.
    is_matroid: bool,
    /// Definition 2.1 for that arm's objective.
    is_submodular: bool,
    /// Proposition 2.3, per chunk where the arm has chunks.
    marginal_returns_monotone: bool,
    /// Achieved over best achievable, in `(0, 1]`.
    greedy_ratio: f64,
    /// Whether Theorem 4.2's hypotheses are met.
    bound_applies: bool,
    /// Chunks with increasing marginal returns.
    violating_chunks: u64,
}

impl Row {
    /// The registered columns, in registration order, with the global verdicts.
    fn registered(&self, c1: bool, c2: bool) -> Vec<(&'static str, String)> {
        vec![
            ("problem", self.problem.clone()),
            ("is_matroid", self.is_matroid.to_string()),
            ("is_submodular", self.is_submodular.to_string()),
            (
                "marginal_returns_monotone",
                self.marginal_returns_monotone.to_string(),
            ),
            ("greedy_ratio", format!("{:.6}", self.greedy_ratio)),
            ("bound_applies", self.bound_applies.to_string()),
            ("violating_chunks", self.violating_chunks.to_string()),
            ("c1_holds", c1.to_string()),
            ("c2_holds", c2.to_string()),
        ]
    }
}

/// A ratio that is `1.0` when both sides are zero, rather than a NaN.
fn ratio(achieved: f64, best: f64) -> f64 {
    if best == 0.0 { 1.0 } else { achieved / best }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-166");

    common::experiment::run(prereg, |run| {
        let started = Instant::now();

        // ── C1: the fixtures, the shipped mesher, and the two controls ──────
        let fixtures: Vec<Fixture> = FIXTURES
            .iter()
            .map(|(name, art)| measure_fixture(name, art))
            .collect();

        for f in &fixtures {
            println!(
                "{:>12}  art {:<12} cells {:>2}  |X| {:>3}  bases {:>6}{}  min {:>2}  max {:>2}  \
                 greedy {:>2}  disjoint-from-min {}  axioms {}",
                f.name,
                f.art,
                f.cells,
                f.ground,
                f.parts.count,
                if f.parts.complete { " " } else { "+" },
                f.parts.min,
                f.parts.max,
                f.greedy,
                f.disjoint_from_small,
                if f.axioms.is_some() {
                    "exhaustive"
                } else {
                    "witness"
                },
            );
        }

        // Example 3.5's `U(3,8)`: the tester must report the exchange property
        // HOLDING, or every `false` it reports elsewhere is worthless.
        let uniform = check_axioms(8, |m| m.count_ones() <= 3);
        // Section 3.3's published failure, verbatim: X = {a, b, c}, w = (2, 1, 1),
        // B = 2. Element 0 is `a`, 1 is `b`, 2 is `c`.
        let knapsack_weights = [2u32, 1, 1];
        let knapsack = check_axioms(3, |m| {
            let mut total = 0u32;
            for i in 0..3 {
                if m >> i & 1 == 1 {
                    total += knapsack_weights[i];
                }
            }
            total <= 2
        });
        println!(
            "\ncontrol U(3,8): sets {} augmentation_holds {}  |  control knapsack(B=2 w=2/1/1): \
             sets {} augmentation_holds {} failures {} witness {:?}",
            uniform.independent_sets,
            uniform.augmentation_holds,
            knapsack.independent_sets,
            knapsack.augmentation_holds,
            knapsack.augmentation_failures,
            knapsack
                .witness
                .map(|(a, b)| (elements(a), elements(b)))
                .unwrap_or_else(|| (String::from("none"), String::from("none"))),
        );

        let shipped: Vec<(&'static str, usize, usize)> = vec![
            (
                "box_exact",
                shipped_quads(&BoxExact::<f64>::canonical(), Merge::Off),
                shipped_quads(&BoxExact::<f64>::canonical(), Merge::Greedy),
            ),
            (
                "sphere",
                shipped_quads(&Sphere::<f64>::canonical(), Merge::Off),
                shipped_quads(&Sphere::<f64>::canonical(), Merge::Greedy),
            ),
        ];
        for (name, off, greedy) in &shipped {
            println!(
                "shipped GreedyQuads {name:>10} at {SHIPPED_SAMPLES}³: merge-off {off:>7} quads, \
                 merge-greedy {greedy:>7} quads, ratio {:>8.3}x",
                *off as f64 / *greedy as f64
            );
        }

        // ── C2: three chunked worlds ────────────────────────────────────────
        let mut checked_nesting = false;
        let fields = vec![
            measure_field(
                &Sphere::<f64>::canonical(),
                "sphere",
                12,
                &mut checked_nesting,
            ),
            measure_field(
                &capped_gyroid::<f64>(),
                "capped_gyroid",
                7,
                &mut checked_nesting,
            ),
            measure_field(
                &FbmTerrain::<f64>::canonical(),
                "fbm_terrain",
                10,
                &mut checked_nesting,
            ),
        ];

        let mut rng = Rng::new(RNG_SEED);
        struct Measured {
            field: &'static str,
            per_axis: i32,
            chunks_total: usize,
            cell_size: f64,
            camera: [f64; 3],
            chunks: usize,
            gains: Vec<[f64; STEPS]>,
            values: Vec<[f64; LEVELS]>,
            order: Vec<usize>,
            optimum: Vec<f64>,
            steps_total: usize,
            strict: u64,
            material: u64,
            monotone: bool,
            levels_differ: bool,
            pair_violations: u64,
            worst_slack: f64,
            worst_chunk: usize,
            worst_second: f64,
            worst_at: usize,
            reversed_strict: u64,
            reversed_ratio: f64,
        }

        let mut measured: Vec<Measured> = Vec::new();
        for fc in &fields {
            let chunks = fc.curve.len();
            let gains: Vec<[f64; STEPS]> = fc.curve.iter().map(|c| c.gain).collect();
            let values: Vec<[f64; LEVELS]> = fc.curve.iter().map(ChunkCurve::value).collect();
            let steps_total = STEPS * chunks;
            let optimum = optimum_curve(&values, steps_total);
            let mut order: Vec<usize> = (0..chunks).collect();
            order.sort_by(|a, b| {
                fc.curve[*a]
                    .distance
                    .total_cmp(&fc.curve[*b].distance)
                    .then(a.cmp(b))
            });
            let (strict, material) = count_violations(&gains);
            let monotone = gains.iter().all(|g| g.iter().all(|v| *v >= 0.0));
            let levels_differ = fc
                .curve
                .iter()
                .any(|c| c.screen[TOP] > c.screen[1] && c.screen[1] > 0.0);
            let (pair_violations, worst_slack) =
                submodular_scan(&values, &mut rng, SUBMODULAR_PAIRS);

            let mut worst_chunk = 0usize;
            let mut worst_second = f64::NEG_INFINITY;
            let mut worst_at = 0usize;
            for (i, c) in fc.curve.iter().enumerate() {
                let (d, at) = c.worst_second_difference();
                if d.total_cmp(&worst_second) == Ordering::Greater {
                    worst_second = d;
                    worst_chunk = i;
                    worst_at = at;
                }
            }

            let reversed = reverse_gains(&fc.curve);
            let reversed_values = values_of(&reversed);
            let (reversed_strict, _) = count_violations(&reversed);
            let reversed_k = (steps_total / 4).max(1);
            let reversed_opt = optimum_curve(&reversed_values, reversed_k);
            let reversed_ratio = ratio(
                greedy_value(&reversed, reversed_k),
                reversed_opt[reversed_k],
            );

            println!(
                "\n{:>14}: {:>5} chunks of {:>5} carry surface, cell {:.6}, camera [{:.2} {:.2} \
                 {:.2}], steps {:>5}\n{:>14}  violating {:>4} (material {:>4}) monotone {} \
                 levels_differ {} pair-violations {}/{} worst-slack {:.3e} worst-second {:.6e} \
                 at step {}  reversed-control violating {:>4} ratio {:.6}",
                fc.field,
                chunks,
                fc.chunks_total,
                fc.cell_size,
                fc.camera[0],
                fc.camera[1],
                fc.camera[2],
                steps_total,
                "",
                strict,
                material,
                monotone,
                levels_differ,
                pair_violations,
                SUBMODULAR_PAIRS,
                worst_slack,
                worst_second,
                worst_at,
                reversed_strict,
                reversed_ratio,
            );

            measured.push(Measured {
                field: fc.field,
                per_axis: fc.per_axis,
                chunks_total: fc.chunks_total,
                cell_size: fc.cell_size,
                camera: fc.camera,
                chunks,
                gains,
                values,
                order,
                optimum,
                steps_total,
                strict,
                material,
                monotone,
                levels_differ,
                pair_violations,
                worst_slack,
                worst_chunk,
                worst_second,
                worst_at,
                reversed_strict,
                reversed_ratio,
            });
        }

        // The synthetic convex control: two chunks, gains built to increase.
        let convex_gains: Vec<[f64; STEPS]> = vec![[1.0, 1.0, 10.0], [2.0, 0.0, 0.0]];
        let convex_values = values_of(&convex_gains);
        let convex_k = 3usize;
        let convex_opt = optimum_curve(&convex_values, convex_k);
        let convex_greedy = greedy_value(&convex_gains, convex_k);
        let convex_ratio = ratio(convex_greedy, convex_opt[convex_k]);
        let (convex_strict, _) = count_violations(&convex_gains);
        let (convex_pairs, _) = {
            let mut r = Rng::new(RNG_SEED ^ 0xC0FFEE);
            submodular_scan(&convex_values, &mut r, SUBMODULAR_PAIRS)
        };
        println!(
            "\ncontrol convex synthetic: greedy {convex_greedy:.3} of optimum {:.3} = \
             {convex_ratio:.6} against the 1-1/e bar {NWF_BOUND:.6}, violating {convex_strict}, \
             pair-violations {convex_pairs}/{SUBMODULAR_PAIRS}",
            convex_opt[convex_k],
        );

        // ── vacuity controls, every one before the first record (M-44) ───────

        // The ladder's own shape is checked at compile time by
        // `_LADDER_IS_LONG_ENOUGH`; what has to be checked at run time is that
        // the four levels it names measurably DIFFER, which is the per-field
        // `levels_differ` control below.
        assert!(
            uniform.empty_independent && uniform.downward_closed && uniform.augmentation_holds,
            "VOID: the exchange tester reported empty={} downward_closed={} \
             augmentation_holds={} on the uniform matroid U(3,8), which Example 3.5 of \
             arXiv:2606.10192 proves IS a matroid -- so a failure reported on the meshing \
             system would be a bug in this tester rather than a result (M-44)",
            uniform.empty_independent,
            uniform.downward_closed,
            uniform.augmentation_holds
        );
        assert!(
            knapsack.empty_independent && knapsack.downward_closed && !knapsack.augmentation_holds,
            "VOID: the exchange tester reported augmentation_holds={} on Section 3.3's \
             published knapsack X={{a,b,c}} w=(2,1,1) B=2, where {{b,c}} and {{a}} are both \
             feasible and neither b nor c fits beside a -- a tester that cannot reproduce a \
             failure someone else published cannot certify one of ours",
            knapsack.augmentation_holds
        );
        {
            let (a, b) = knapsack
                .witness
                .expect("the knapsack control has a failing pair");
            assert!(
                a.count_ones() == 2 && b.count_ones() == 1 && b == 1,
                "VOID: the knapsack control failed on ({}, {}) rather than on Section 3.3's \
                 ({{b,c}}, {{a}}), so the tester is failing for a reason the source does not \
                 name",
                elements(a),
                elements(b)
            );
        }
        for f in &fixtures {
            assert!(
                f.parts.min >= 1 && f.parts.max > f.parts.min,
                "VOID: {} has bases of {}..{} -- a region with one partition size has nothing \
                 to exchange, so `disjoint_from_small = 0` there would mean 'no larger basis \
                 exists' rather than 'no augmentation is possible' (M-44)",
                f.name,
                f.parts.min,
                f.parts.max
            );
            assert_eq!(
                f.parts.max, f.cells,
                "VOID: {}'s largest basis is {} against {} cells -- the all-unit-cells \
                 partition must be enumerated or the certificate's larger set is not the one \
                 the header claims",
                f.name, f.parts.max, f.cells
            );
            assert_eq!(
                f.disjoint_from_small, 0,
                "P-166: {} has {} rectangles disjoint from its minimum basis, so that basis \
                 does not cover the region and the exchange certificate is void",
                f.name, f.disjoint_from_small
            );
            if let Some(ax) = &f.axioms {
                assert!(
                    ax.empty_independent && ax.downward_closed,
                    "P-166: {}'s packing system is not an independence system (empty={} \
                     downward_closed={}), so Definition 3.1(3) is being asked about the wrong \
                     object",
                    f.name,
                    ax.empty_independent,
                    ax.downward_closed
                );
            }
        }
        for (name, off, greedy) in &shipped {
            assert!(
                off > greedy,
                "VOID: shipped GreedyQuads on {name} emitted {off} quads with merging off and \
                 {greedy} with it on -- if merging changes no count then Merge is inert here \
                 and the shipped arm witnesses no cardinality difference (M-44)"
            );
        }
        for m in &measured {
            assert!(
                m.chunks >= MIN_CURVE_CHUNKS,
                "VOID: {}'s marginal-return curve spans {} surface chunks of {} visited, below \
                 the registered {MIN_CURVE_CHUNKS} -- 'diminishing' would be a claim about a \
                 handful of points",
                m.field,
                m.chunks,
                m.chunks_total
            );
            assert!(
                m.levels_differ,
                "VOID: no chunk of {} has screen[3] > screen[1] > 0, so the four levels named \
                 do not measurably differ and 'spanning at least two LOD levels' is nominal \
                 rather than real",
                m.field
            );
            // Proposition 2.3 and Definition 2.1 are the same statement, so the
            // two instruments must agree at the extremal witness. A silent
            // disagreement would let either of them report a pass.
            if m.worst_second > 0.0 {
                assert!(
                    witness_violates(&m.values[m.worst_chunk], m.worst_at),
                    "P-166: {}'s chunk {:?} has gain[{}]-gain[{}] = {:.6e} > 0, so \
                     Proposition 2.3 fails there -- but Definition 2.1 at the constructed \
                     witness pair holds, and the two are equivalent, so one implementation is \
                     wrong",
                    m.field,
                    m.worst_chunk,
                    m.worst_at + 1,
                    m.worst_at,
                    m.worst_second
                );
            }
            assert!(
                m.reversed_strict > 0,
                "VOID: reversing every chunk ladder of {} produced {} chunks with increasing \
                 marginal returns -- the detector cannot be shown able to fire on this \
                 field's own magnitudes, so a zero on the measured arm is a silence (M-44)",
                m.field,
                m.reversed_strict
            );
        }
        assert!(
            convex_strict == 1 && convex_ratio < NWF_BOUND,
            "VOID: the synthetic convex control reported violating={convex_strict} and \
             greedy_ratio={convex_ratio:.6} against the 1-1/e bar {NWF_BOUND:.6}. It is built \
             from gains (1,1,10) and (2,0,0) at budget 3, where greedy takes 2+1+1 = 4 and the \
             optimum takes all three of the first chunk for 12, so exactly one chunk must be \
             flagged and the ratio must be 1/3. Anything else means the second-difference test \
             or the knapsack DP cannot report bad news (M-44)"
        );
        assert!(
            convex_pairs > 0,
            "VOID: Definition 2.1's random pair scan found no violation on the synthetic \
             convex control, whose second chunk is concave and whose first is not -- the \
             union/intersection instrument cannot be shown able to fire"
        );

        // ── global clause verdicts ──────────────────────────────────────────
        //
        // C1 holds if the exchange property is exhibited failing on greedy
        // meshing's own set system, with a certificate, and the tester is shown
        // able to report the property holding. C2 holds only if no chunk anywhere
        // has increasing marginal returns -- one such chunk is the registered
        // falsifier and is what this harness was built to look for.
        let c1_holds = fixtures
            .iter()
            .all(|f| f.parts.max > f.parts.min && f.disjoint_from_small == 0)
            && uniform.augmentation_holds
            && !knapsack.augmentation_holds
            && shipped.iter().all(|(_, off, greedy)| off > greedy);
        let total_violating: u64 = measured.iter().map(|m| m.strict).sum();
        let c2_holds = total_violating == 0;

        println!(
            "\nC1 {} -- the exchange property fails on greedy meshing's packing system on all \
             {} fixtures and on the shipped mesher, with U(3,8) confirming the tester can \
             report it holding\nC2 {} -- {} chunks of {} across three fields have increasing \
             marginal returns",
            if c1_holds { "HELD" } else { "FALSIFIED" },
            fixtures.len(),
            if c2_holds { "HELD" } else { "FALSIFIED" },
            total_violating,
            measured.iter().map(|m| m.chunks).sum::<usize>(),
        );

        let wall_ns = started.elapsed().as_nanos();

        // ── C1 rows ─────────────────────────────────────────────────────────
        //
        // `is_submodular` and `marginal_returns_monotone` are TRUE on every
        // meshing row and that is the finding, not filler: the objective is |S|,
        // which is modular, so its marginals are constant and constant is
        // non-increasing. Submodularity is not what greedy meshing lacks.
        // `bound_applies` is FALSE because Theorem 4.2 wants a monotone
        // submodular objective MAXIMISED under a cardinality constraint, and this
        // is a cardinality MINIMISATION under a covering constraint whose system
        // fails Definition 3.1(3). `violating_chunks` is 0 because this arm has no
        // chunks, which `arm_has_chunks` says on the row.
        for f in &fixtures {
            let (wa, wb) = f.axioms.as_ref().and_then(|ax| ax.witness).map_or_else(
                || (String::from("certificate"), String::from("certificate")),
                |(a, b)| (elements(a), elements(b)),
            );
            let mut row = Row {
                problem: format!("greedy_meshing/{}", f.name),
                is_matroid: false,
                is_submodular: true,
                marginal_returns_monotone: true,
                greedy_ratio: f.parts.min as f64 / f.greedy as f64,
                bound_applies: false,
                violating_chunks: 0,
            }
            .registered(c1_holds, c2_holds);
            // ── extras (M-273) ──
            row.extend([
                ("arm_kind", String::from("c1_fixture")),
                ("arm_verdict", c1_holds.to_string()),
                ("arm_has_chunks", String::from("false")),
                ("is_control", String::from("false")),
                ("fixture_art", f.art.to_string()),
                ("fixture_cells", f.cells.to_string()),
                ("ground_rectangles", f.ground.to_string()),
                ("bases_enumerated", f.parts.count.to_string()),
                ("bases_complete", f.parts.complete.to_string()),
                ("min_basis", f.parts.min.to_string()),
                ("min_basis_exact", String::from("true")),
                ("max_basis", f.parts.max.to_string()),
                ("greedy_basis", f.greedy.to_string()),
                (
                    "basis_size_ratio",
                    format!("{:.6}", f.parts.max as f64 / f.parts.min as f64),
                ),
                (
                    "p_system_lower_bound",
                    format!("{:.6}", f.parts.max as f64 / f.parts.min as f64),
                ),
                (
                    "min_basis_rects",
                    f.parts
                        .small
                        .iter()
                        .map(|r| r.tag())
                        .collect::<Vec<String>>()
                        .join("|"),
                ),
                ("disjoint_from_small", f.disjoint_from_small.to_string()),
                (
                    "axioms",
                    String::from(if f.axioms.is_some() {
                        "exhaustive"
                    } else {
                        "witness"
                    }),
                ),
                (
                    "independent_sets",
                    f.axioms.as_ref().map_or_else(
                        || String::from("uncounted"),
                        |ax| ax.independent_sets.to_string(),
                    ),
                ),
                (
                    "augmentation_failures",
                    f.axioms.as_ref().map_or_else(
                        || String::from("certificate"),
                        |ax| ax.augmentation_failures.to_string(),
                    ),
                ),
                ("witness_a", wa),
                ("witness_b", wb),
                ("wall_ns", wall_ns.to_string()),
            ]);
            run.record(&row);
        }

        for (name, off, greedy) in &shipped {
            let mut row = Row {
                problem: format!("greedy_meshing/shipped@{name}_{SHIPPED_SAMPLES}"),
                is_matroid: false,
                is_submodular: true,
                marginal_returns_monotone: true,
                // The best basis this arm found is greedy's own, so the ratio is
                // 1 against it. The exact minimum is R-165's row, which is what
                // `min_basis_exact = false` records.
                greedy_ratio: 1.0,
                bound_applies: false,
                violating_chunks: 0,
            }
            .registered(c1_holds, c2_holds);
            // ── extras (M-273) ──
            row.extend([
                ("arm_kind", String::from("c1_shipped")),
                ("arm_verdict", c1_holds.to_string()),
                ("arm_has_chunks", String::from("false")),
                ("is_control", String::from("false")),
                ("field", (*name).to_string()),
                ("resolution", SHIPPED_SAMPLES.to_string()),
                ("min_basis", greedy.to_string()),
                ("min_basis_exact", String::from("false")),
                ("max_basis", off.to_string()),
                ("greedy_basis", greedy.to_string()),
                ("quads_merge_off", off.to_string()),
                ("quads_merge_greedy", greedy.to_string()),
                (
                    "basis_size_ratio",
                    format!("{:.6}", *off as f64 / *greedy as f64),
                ),
                (
                    "p_system_lower_bound",
                    format!("{:.6}", *off as f64 / *greedy as f64),
                ),
                ("axioms", String::from("witness")),
                ("wall_ns", wall_ns.to_string()),
            ]);
            run.record(&row);
        }

        for (name, ax, matroid, applies) in [
            ("control/uniform_matroid_U3_8", &uniform, true, true),
            ("control/knapsack_B2_w211", &knapsack, false, false),
        ] {
            let (wa, wb) = ax.witness.map_or_else(
                || (String::from("none"), String::from("none")),
                |(a, b)| (elements(a), elements(b)),
            );
            let mut row = Row {
                problem: name.to_string(),
                is_matroid: matroid,
                is_submodular: true,
                marginal_returns_monotone: true,
                greedy_ratio: 1.0,
                bound_applies: applies,
                violating_chunks: 0,
            }
            .registered(c1_holds, c2_holds);
            // ── extras (M-273) ──
            row.extend([
                ("arm_kind", String::from("c1_control")),
                ("arm_verdict", ax.augmentation_holds.to_string()),
                ("arm_has_chunks", String::from("false")),
                ("is_control", String::from("true")),
                ("axioms", String::from("exhaustive")),
                ("independent_sets", ax.independent_sets.to_string()),
                ("empty_independent", ax.empty_independent.to_string()),
                ("downward_closed", ax.downward_closed.to_string()),
                (
                    "augmentation_failures",
                    ax.augmentation_failures.to_string(),
                ),
                ("witness_a", wa),
                ("witness_b", wb),
                ("wall_ns", wall_ns.to_string()),
            ]);
            run.record(&row);
        }

        // ── C2 rows ─────────────────────────────────────────────────────────
        for m in &measured {
            let arm_submodular = m.strict == 0 && m.pair_violations == 0;
            let arm_bound = m.monotone && arm_submodular;
            let arm_c2 = m.strict == 0;

            let mut budgets: Vec<usize> = vec![
                1,
                m.steps_total / 16,
                m.steps_total / 8,
                m.steps_total / 4,
                m.steps_total / 2,
                m.steps_total,
            ];
            budgets.retain(|k| *k >= 1);
            budgets.sort_unstable();
            budgets.dedup();

            for k in &budgets {
                let k = *k;
                let f_greedy = greedy_value(&m.gains, k);
                let f_opt = m.optimum[k];
                let f_distance = distance_value(&m.values, &m.order, k);
                let g_ratio = ratio(f_greedy, f_opt);
                let mut row = Row {
                    problem: format!("lod_budget/{}@k={k}", m.field),
                    // `|S| <= k` is Example 3.5's uniform matroid; the U(3,8)
                    // control row is where that tester is exercised.
                    is_matroid: true,
                    is_submodular: arm_submodular,
                    marginal_returns_monotone: m.strict == 0,
                    greedy_ratio: g_ratio,
                    bound_applies: arm_bound,
                    violating_chunks: m.strict,
                }
                .registered(c1_holds, c2_holds);
                // ── extras (M-273) ──
                row.extend([
                    ("arm_kind", String::from("c2_budget")),
                    ("arm_verdict", arm_c2.to_string()),
                    ("arm_has_chunks", String::from("true")),
                    ("is_control", String::from("false")),
                    ("field", m.field.to_string()),
                    ("chunks_per_axis", m.per_axis.to_string()),
                    ("chunks_visited", m.chunks_total.to_string()),
                    ("curve_chunks", m.chunks.to_string()),
                    ("levels_spanned", LEVELS.to_string()),
                    ("refinement_steps", STEPS.to_string()),
                    ("steps_total", m.steps_total.to_string()),
                    ("budget_k", k.to_string()),
                    ("cell_size", format!("{:.9}", m.cell_size)),
                    (
                        "camera",
                        format!("{:.4}|{:.4}|{:.4}", m.camera[0], m.camera[1], m.camera[2]),
                    ),
                    ("downsample_op", String::from(Downsample::Decimate.name())),
                    ("f_greedy", format!("{f_greedy:.6e}")),
                    ("f_optimal", format!("{f_opt:.6e}")),
                    ("f_distance_order", format!("{f_distance:.6e}")),
                    (
                        "share_over_distance_order",
                        format!("{:.6}", ratio(f_greedy - f_distance, f_opt)),
                    ),
                    (
                        "distance_order_ratio",
                        format!("{:.6}", ratio(f_distance, f_opt)),
                    ),
                    ("one_minus_1_over_e", format!("{NWF_BOUND:.6}")),
                    ("ratio_meets_nwf_bound", (g_ratio >= NWF_BOUND).to_string()),
                    ("greedy_is_optimal", (g_ratio >= 1.0).to_string()),
                    ("normalized", String::from("true")),
                    ("monotone", m.monotone.to_string()),
                    ("violating_chunks_material", m.material.to_string()),
                    ("submodular_pairs", SUBMODULAR_PAIRS.to_string()),
                    ("submodular_violations", m.pair_violations.to_string()),
                    ("submodular_worst_slack", format!("{:.6e}", m.worst_slack)),
                    (
                        "submodular_scan_saw_it",
                        (m.pair_violations > 0).to_string(),
                    ),
                    (
                        "witness_sees_violation",
                        (m.worst_second > 0.0
                            && witness_violates(&m.values[m.worst_chunk], m.worst_at))
                        .to_string(),
                    ),
                    ("worst_second_difference", format!("{:.6e}", m.worst_second)),
                    ("rng_seed", format!("{RNG_SEED:#x}")),
                    ("wall_ns", wall_ns.to_string()),
                ]);
                run.record(&row);
            }

            // The single most extremal chunk, named, so the CSV always exhibits
            // the closest call rather than only its aggregate.
            let fc = fields
                .iter()
                .find(|f| f.field == m.field)
                .expect("every measured arm came from a field curve");
            let c = &fc.curve[m.worst_chunk];
            let k = (m.steps_total / 4).max(1);
            let g_ratio = ratio(greedy_value(&m.gains, k), m.optimum[k]);
            let mut row = Row {
                problem: format!("lod_budget/{}@worst_chunk", m.field),
                is_matroid: true,
                is_submodular: arm_submodular,
                marginal_returns_monotone: m.strict == 0,
                greedy_ratio: g_ratio,
                bound_applies: arm_bound,
                violating_chunks: m.strict,
            }
            .registered(c1_holds, c2_holds);
            // ── extras (M-273) ──
            row.extend([
                ("arm_kind", String::from("c2_worst_chunk")),
                ("arm_verdict", arm_c2.to_string()),
                ("arm_has_chunks", String::from("true")),
                ("is_control", String::from("false")),
                ("field", m.field.to_string()),
                ("curve_chunks", m.chunks.to_string()),
                ("levels_spanned", LEVELS.to_string()),
                ("refinement_steps", STEPS.to_string()),
                ("steps_total", m.steps_total.to_string()),
                ("budget_k", k.to_string()),
                (
                    "worst_chunk_coords",
                    format!("{}|{}|{}", c.coords[0], c.coords[1], c.coords[2]),
                ),
                ("worst_chunk_distance", format!("{:.6}", c.distance)),
                ("worst_chunk_lipschitz", format!("{:.6e}", c.lipschitz)),
                (
                    "worst_chunk_geometric",
                    format!(
                        "{:.6e}|{:.6e}|{:.6e}|{:.6e}",
                        c.geometric[0], c.geometric[1], c.geometric[2], c.geometric[3]
                    ),
                ),
                (
                    "worst_chunk_screen",
                    format!(
                        "{:.6e}|{:.6e}|{:.6e}|{:.6e}",
                        c.screen[0], c.screen[1], c.screen[2], c.screen[3]
                    ),
                ),
                (
                    "worst_chunk_gains",
                    format!("{:.6e}|{:.6e}|{:.6e}", c.gain[0], c.gain[1], c.gain[2]),
                ),
                ("worst_second_difference", format!("{:.6e}", m.worst_second)),
                ("worst_second_at_step", (m.worst_at + 1).to_string()),
                ("worst_chunk_violates", c.violates().to_string()),
                (
                    "worst_chunk_violates_materially",
                    c.violates_materially().to_string(),
                ),
                (
                    "submodular_scan_saw_it",
                    (m.pair_violations > 0).to_string(),
                ),
                (
                    "witness_sees_violation",
                    (m.worst_second > 0.0
                        && witness_violates(&m.values[m.worst_chunk], m.worst_at))
                    .to_string(),
                ),
                ("violating_chunks_material", m.material.to_string()),
                ("wall_ns", wall_ns.to_string()),
            ]);
            run.record(&row);

            // The reversed-real control: the same magnitudes, increasing.
            let mut row = Row {
                problem: format!("control/lod_reversed_{}", m.field),
                is_matroid: true,
                is_submodular: false,
                marginal_returns_monotone: false,
                greedy_ratio: m.reversed_ratio,
                bound_applies: false,
                violating_chunks: m.reversed_strict,
            }
            .registered(c1_holds, c2_holds);
            // ── extras (M-273) ──
            row.extend([
                ("arm_kind", String::from("c2_control")),
                ("arm_verdict", String::from("false")),
                ("arm_has_chunks", String::from("true")),
                ("is_control", String::from("true")),
                ("field", m.field.to_string()),
                ("curve_chunks", m.chunks.to_string()),
                ("levels_spanned", LEVELS.to_string()),
                ("refinement_steps", STEPS.to_string()),
                ("steps_total", m.steps_total.to_string()),
                ("budget_k", (m.steps_total / 4).max(1).to_string()),
                ("one_minus_1_over_e", format!("{NWF_BOUND:.6}")),
                (
                    "ratio_meets_nwf_bound",
                    (m.reversed_ratio >= NWF_BOUND).to_string(),
                ),
                ("wall_ns", wall_ns.to_string()),
            ]);
            run.record(&row);
        }

        let mut row = Row {
            problem: String::from("control/lod_convex_synthetic"),
            is_matroid: true,
            is_submodular: false,
            marginal_returns_monotone: false,
            greedy_ratio: convex_ratio,
            bound_applies: false,
            violating_chunks: convex_strict,
        }
        .registered(c1_holds, c2_holds);
        // ── extras (M-273) ──
        row.extend([
            ("arm_kind", String::from("c2_control")),
            ("arm_verdict", String::from("false")),
            ("arm_has_chunks", String::from("true")),
            ("is_control", String::from("true")),
            ("field", String::from("synthetic")),
            ("curve_chunks", convex_gains.len().to_string()),
            ("levels_spanned", LEVELS.to_string()),
            ("refinement_steps", STEPS.to_string()),
            ("steps_total", (STEPS * convex_gains.len()).to_string()),
            ("budget_k", convex_k.to_string()),
            ("f_greedy", format!("{convex_greedy:.6e}")),
            ("f_optimal", format!("{:.6e}", convex_opt[convex_k])),
            (
                "worst_chunk_gains",
                format!(
                    "{:.6e}|{:.6e}|{:.6e}",
                    convex_gains[0][0], convex_gains[0][1], convex_gains[0][2]
                ),
            ),
            ("submodular_pairs", SUBMODULAR_PAIRS.to_string()),
            ("submodular_violations", convex_pairs.to_string()),
            ("one_minus_1_over_e", format!("{NWF_BOUND:.6}")),
            (
                "ratio_meets_nwf_bound",
                (convex_ratio >= NWF_BOUND).to_string(),
            ),
            ("wall_ns", wall_ns.to_string()),
        ]);
        run.record(&row);
    });
}
