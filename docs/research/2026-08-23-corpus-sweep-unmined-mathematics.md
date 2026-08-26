# The corpus sweep, and the parent theory that has been named since day two and never used

**Date:** 2026-08-23 · **Corpus:** home-still, 9,502 documents, 9,486 embedded, 289,980 chunks
**Supersedes** the acquisition judgements in `2026-08-23-unmined-mathematics-for-meshing.md`, which was
written web-only. Three of its recommendations were wrong about what you already own.

---

## 0. The measurement

For each candidate area: **is it in the corpus**, and **does the repository ever cite it?** That second
column is the one that matters — an area with a paper in the library and no reference anywhere in
`FINDINGS.md`, `BACKLOG.md`, `CLAUDE.md`, `docs/` or `crates/` is unmined *with the source already
paid for*.

| Paper | In corpus | Cited in repo |
|---|---|---|
| **Robins, Wood & Sheppard — Morse complexes from grayscale images** `10.1109/tpami.2011.95` | ✅ | **0 files** |
| **Delgado-Friedrichs, Robins & Sheppard — skeletonization via DMT** `10.1109/tpami.2014.2346172` | ✅ | **0 files** |
| **Saye — high-order distances to implicit surfaces** `10.2140/camcos.2014.9.107` | ✅ | **0 files** |
| **Bender, Kuszmaul, Teng & Wang — optimal cache-oblivious mesh layouts** `10.48550/arXiv.0705.1033` | ✅ | **0 files** |
| **Curry — Sheaves, Cosheaves and Applications** `10.48550/arXiv.1303.3255` | ✅ | **0 files** |
| **Reuter et al. — barycentric Bernstein polynomial-system solver** `10.1007/s00371-007-0184-x` | ✅ | **0 files** |
| **Lévy — L₂ semi-discrete optimal transport in 3D** `10.1051/m2an/2015055` | ✅ | **0 files** |
| **Practical box splines on the BCC lattice** `10.1109/tvcg.2007.70429` | ✅ (today) | **0 files** |
| **Klötzl et al. — local bilinear Jacobi sets** `10.1007/s00371-022-02557-4` | ✅ (today) | **0 files** |
| Blu, Thévenaz & Unser — shifted linear interpolation `10.1109/tip.2004.826093` | ✅ | 1 — a procurement list only |
| Ghrist & Riess — cellular sheaves of lattices `10.4310/hha.2022.v24.n1.a16` | ✅ | 1 — an acquisition list only |
| Carr et al. — data-parallel contour trees `10.1109/tvcg.2021.3064385` | ✅ | 1 — a research doc, never a ticket |
| Weber et al. — parallel peak pruning `10.1109/ldav.2016.7874312` | ✅ | 1 — same doc |
| Mourrain & Pavone — subdivision solvers `10.1016/j.jsc.2008.04.016` | ✅ | 1 — **my memo from this morning** |
| Shewchuk — adaptive precision predicates `10.1007/pl00009321` | ✅ | 5 — **mined** |
| Attene — indirect predicates `10.1016/j.cad.2020.102856` | ✅ | 4 — **mined** |
| Villar et al. — *Scalars are universal* `10.48550/arxiv.2106.06610` | ✅ | 2, incl. `dual_contouring/solve.rs` — **mined** |

---

## 1. The find: discrete Morse theory, and the sentence you wrote on 2026-08-12

`docs/research/2026-08-12-axes-and-vocabulary.md:117` says, in your own words:

> **Morse theory** | Study of a function's topology via its critical points | **The parent theory.
> *Ambiguity* is really "a critical point falls inside a cell"**

That is the correct diagnosis of the entire A-002 series — the asymptotic decider, MC33, the interior
test, body saddles, tunnels versus twelve-vertex contours, singular faces, the `[9,3]` cells — and in
the eleven days since, **nothing has cited a Morse-theoretic result.** The A-002 line instead
hand-derived quadratics in the monomial basis, and paid for it in M-207, M-219, M-221, ✗22, M-228 and
M-231.

**The constructive algorithm is in your library and is aimed at exactly your data model.**
Robins, Wood & Sheppard, *Theory and Algorithms for Constructing Discrete Morse Complexes from
Grayscale Digital Images*, **IEEE TPAMI 2011**, `10.1109/tpami.2011.95`. From the paper, verbatim:

> *"We present an algorithm for determining the Morse complex of a 2- or 3-dimensional grayscale
> digital image. Each cell in the Morse complex corresponds to a topological change in the level sets
> (i.e. a critical point) of the grayscale image. Since more than one critical point may be associated
> with a single image voxel, we model digital images by cubical complexes. A new homotopic algorithm
> is used to construct a discrete Morse function on the cubical complex that agrees with the digital
> image and **has exactly the number and type of critical cells necessary to characterize the
> topological changes in the level sets**."*

Four properties, each of which lands on a live problem in the ledger:

**It is per-voxel local.** The core routine is `ProcessLowerStars`, which operates on one voxel's
**lower star** — the cells whose values are at or below it. That is a `2×2×2`-neighbourhood
computation, so it is chunkable, and it is **edit-local in exactly the sense R-020/R-022 measured**
(M-311: 792 dirty cells constant across a 64× lattice). A discrete Morse complex maintained under
editing is the same shape of question as `Air`'s connectivity, and you have already built that.

**It carries a correctness proof, not an empirical claim.** Lemmas 7–10 plus a Mayer–Vietoris argument
establish a **one-to-one correspondence** between the critical cells the algorithm finds and the
homology changes of the lower level cuts. Compare P-52's monotone-edge gate, which I proposed this
morning from arXiv:2608.12142 — whose own abstract says *"correct with regards to critical points **in
our experiments**"*. **This is the theorem P-52 was groping for, it is stronger, and you already own
it.** P-52 should be withdrawn and re-registered against this.

**It is honest about its worst case, and the honesty is usable.** The V-path traversal is `O(N²)` on a
contrived 3D construction with `O(N)` critical 1-cells each meeting `O(N)` critical 2-cells — and the
authors say plainly: *"none of the images or objects that we have studied exhibit scaling behaviour
worse than O(N)."* That is a measurable claim on your eight reference fields plus `fuel` and `bonsai`,
and measuring it is a one-bench experiment.

**It generalises the two structures the ledger already reaches for.** The paper states the Morse
complex *"is a generalisation of the watershed and component tree"*. Your `Air` connectivity is a
component tree in disguise; R-021's blocked contour-tree ticket wants the other one.

**What it would buy, stated against specific rows.** A topological gate on `gyroid` and `fbm_terrain`,
where `CLAUDE.md` records that χ cannot be asserted at all — those two fields currently have no
topological check beyond manifoldness. A principled account of ambiguity: a cell is ambiguous exactly
when it contains a critical cell of the discrete Morse function, which is a *computed* fact rather
than a case-table lookup. And, per Kovalevsky as the paper cites him, *"a cubical complex is the only
topologically consistent model of a digital image"* — which is the same claim P-41/M-338's
well-composedness census is a special case of.

**Companion, also in corpus, also uncited:** Delgado-Friedrichs, Robins & Sheppard, *Skeletonization
and Partitioning of Digital Images Using Discrete Morse Theory*, `10.1109/tpami.2014.2346172`. That is
the medial-axis line — ✗27, M-324, M-325, O-19 — from the discrete side, with the λ-filter question
O-19 is blocked on posed combinatorially rather than as a continuous threshold.

---

## 2. Shifted linear interpolation — flagged in your own procurement doc, never acted on

`docs/research/2026-08-18-corpus-audit-and-procurement.md:128` lists Blu/Thévenaz/Unser under *vertex
placement and quality* as **present**. Nothing has used it.

**What the paper claims, read from the text rather than the title.** Standard piecewise-linear
interpolation is improved by shifting the sampling knots by a **fixed, signal-independent** amount
while enforcing the interpolation property. The optimal shift is *"nonzero and close to 1/5"*. The
measured consequences:

- a gain of *"almost 8 dB in the neighborhood of ω = 0"*, and the shifted method beats the unshifted
  one *"for signals with frequency range up to 3/4th of the sampling bandwidth"*;
- *"a quality that is similar to that of the computationally more costly 'high-quality' cubic
  convolution"*, at *"a similar computational cost"*;
- the overhead is a prefilter of support `W` over `N_in` samples — `2·W·N_in` operations, which they
  call *"usually a negligible overhead"* and note can be precomputed;
- and, directly relevant if you ever move to a tricubic filter: *"the same shifting trick may also be
  used for other kernels than Λ(x)… the optimal shift has then to be determined on a case-by-case
  basis."*

**Why this is aimed at your single most-executed operation.** Every Marching Cubes vertex is a linear
interpolation between two corner samples. M-12 measured the error falling like `h²` with a fitted
constant; F-007/M-250 measured that refining the crossing on the real field buys **13–15%** on curved
fields and nothing on the CSG one. A fixed shift that recovers cubic-convolution quality at linear cost
attacks the same constant from a completely different direction — approximation theory rather than
root-finding.

**Three objections that belong in the registration, not in a footnote.** The shift moves the surface,
so **every one of the 216 golden hashes changes** — this is a re-baseline, not a tweak. The prefilter is
a *recursive digital filter over the samples*, which makes reconstruction non-local and therefore
interacts badly with chunking and with M-32's seam arithmetic; whether a chunk-local approximation
preserves the gain is itself the experiment. And the entire A-002 apparatus is derived from the
**trilinear interpolant** — shifting changes the reconstruction, so the decider, the interior test and
the body-saddle algebra would all need re-deriving. The honest first experiment is therefore *how much
of the h² constant is recoverable*, measured on a 1-D edge in isolation, before anything touches an
extractor.

---

## 3. Three more corpus-resident areas with a hook into a measured row

**Cache-oblivious layout theory.** Bender, Kuszmaul, Teng & Wang, `10.48550/arXiv.0705.1033`,
*Optimal Cache-Oblivious Mesh Layouts*. The A-023/A-024 saga (M-285, M-286, M-287, O-11) fought cache
behaviour with `size[0] | 1` — a **cache-aware** fix, found empirically after a 3.37× penalty at 128³
was traced to a 64 KiB plane stride. Cache-oblivious theory says there is a provably optimal layout
that needs no knowledge of the cache parameters at all, which is the right property for a crate that
ships to unknown hardware and whose own M-281 says a millisecond is a property of the binary. The
measurable question: does a van Emde Boas or Hilbert-order sample layout beat the odd-stride flat array
on the cost-per-sample curve, and does it remove the residual `+11%` rise M-287 left across 16³…256³?

**High-order distance computation.** Saye, `10.2140/camcos.2014.9.107`, *High-order methods for
computing distances to implicitly defined surfaces*. Your entire S-001…S-009 construction line is first
or second order — exact transform quantised to the grid (M-251's one-spacing floor), Godunov fast
sweeping, fast marching, jump flooding. M-257 measured the approximate GPU method beating both exact
CPU ones *because the seeding, not the solver, carries the accuracy*. A high-order closest-point
projection is the seeding improved, and it is the natural next measurement on that line.

**Cellular sheaves for the seam question.** Curry, `10.48550/arXiv.1303.3255` (uncited) and Ghrist &
Riess, `10.4310/hha.2022.v24.n1.a16` (in an acquisition list only). M-32, M-307, M-278 and ✗18 are all
one question — *when do locally-meshed chunks agree on overlaps well enough to glue into a global
surface?* — which is literally the sheaf-gluing axiom. I flag this as the most speculative item here:
the formalism is exact and the algorithmic payoff is unproven. Worth reading before proposing, not
proposing before reading.

**And one bridge between two of my own findings.** Lévy, `10.1051/m2an/2015055`, *A Numerical Algorithm
for L₂ Semi-Discrete Optimal Transport in 3D* — uncited, in corpus. Semi-discrete optimal transport in
3D is computed via **power (Laguerre) diagrams**, which is precisely the object §3 of this morning's
memo identified as the one place the CSG tape becomes literally tropical. Optimal transport and the
power-diagram substitution are the same machinery approached from two directions, and Lévy has the 3-D
implementation.

---

## 4. Corrections I owe to this morning's memo

**Three "acquire" labels were wrong.** Mourrain & Pavone `10.1016/j.jsc.2008.04.016` is in the corpus —
so is Reuter et al.'s barycentric Bernstein solver `10.1007/s00371-007-0184-x`, which I did not know
existed. The certified-subdivision route of §4 needs **no acquisition at all**.

**One area I called unmined is mined.** I listed invariant theory as Tier-2 item E — "equivariance by
construction rather than by sorting". Villar et al., *Scalars are universal* `10.48550/arxiv.2106.06610`
is in the corpus **and is cited in `crates/isomesh/src/dual_contouring/solve.rs`**. The A-016 line
already knows about it. What remains open is narrower and should be stated that way: M-177 proved
reordering cannot buy negation equivariance, and the O(3) scalar-parameterisation gives a *construction*
that would — but nobody has written the vertex rule that way. That is a smaller claim than I made.

**Two negatives are now corpus-backed rather than web-backed, which makes them worth more.** Tropical
geometry: the best hit for *"tropical geometry min-plus algebra polyhedral complex"* across 9,486
embedded documents scores **0.570** and is a paper on surfaces in cusped 3-manifolds. Basu & Perrucci's
multi-affine hypersurface bound: **absent** — the top hit for the query is the same 3-manifolds paper.
The Lovász-extension Morse work (Jost et al.): **absent**. Additively-weighted Voronoi / Apollonius as a
dedicated source: **absent**. So the four things I flagged as the genuinely novel directions are novel
against your library too, not merely against a web search.

---

## 5. Acquisition list, ranked by what it unblocks

| Priority | Identifier | What it unblocks |
|---|---|---|
| 1 | `10.1016/j.aim.2023.108982` (arXiv:2204.01595) Basu & Perrucci | The `b₀ ≤ 2^(d−1)` multi-affine bound — turns the case-table budget from searched into derived |
| 2 | `10.1007/s00454-022-00461-1` (arXiv:2003.06021) Jost et al. | Lovász-extension Morse theory — the proof of why Marching Tetrahedra has no ambiguity |
| 3 | `10.1007/978-3-319-18720-4_47` Boutry/Géraud/Najman ISMM 2015 | **P-46's actual method.** The repair that ran and failed cited the *tutorial* (`10.1007/s10851-017-0769-6`, now in corpus) rather than the self-dual repair |
| 4 | `10.1145/3730868` Plateau-Holleville et al. 2025 | 3-D Apollonius on GPU — the computable form of the pure-`Add` seam |
| 5 | `10.1007/978-1-4614-4340-7` Scholtes, *Piecewise Differentiable Equations* | The "essentially active index set" — the analytic name for P-39's surviving brush set |
| 6 | `10.1145/49155.51123` Rossignac & Requicha, *Active zones in CSG* | The 1986 successor to Tilove that P-39 should cite alongside 1984 |
| 7 | arXiv:2604.00157 / `10.1145/3799902.3811116` Carrera et al. | P-49's tangency-energy vertex rule |
| 8 | Lipski / Ohtsuki, minimum rectangle partition | The optimum M-56's greedy ratio has never been measured against |

---

## 6. What I would register next, revised

**P-52 is withdrawn.** The monotone-edge gate was the weaker version of a theorem you already own.
Replace it with a discrete-Morse critical-cell census over the eight reference fields plus `fuel` and
`bonsai`: does the count of critical cells predict the ambiguous-cell count, and does the V-path
traversal scale as `O(N)` on real fields as the authors claim? Both halves are falsifiable, both are
counts rather than timings, and the second is the paper's own stated worst case tested against your
data. It also has a natural control — on `sphere`, where no cell is ambiguous, the critical cells
should be exactly the surface's own Morse structure and nothing else.

**P-48, P-50, P-51, P-53, P-54, P-55 stand as written.** P-49's source is still an acquisition.

**One new candidate, from §2:** the shifted-linear question, scoped to a 1-D edge before it touches any
extractor. The registered quantity is *what fraction of M-12's measured `h²` constant a fixed shift
recovers on a single grid edge*, with the golden-hash cost and the prefilter-locality objection stated
in the registration rather than discovered during it.

---

**Method note for the ledger.** This sweep found four papers in the corpus that the repository has
never cited and one that only a procurement list mentions — and it found them in about twenty
semantic queries. ✗4's rule is *"presence in the corpus is decided by `catalog_read`, never by
search"*; its converse is now also earned: **absence from the ledger is not evidence of absence from
the library**, and the gap between "acquired" and "used" is currently at least nine papers wide.
