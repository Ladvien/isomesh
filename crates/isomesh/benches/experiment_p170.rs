//! **P-170 — Marching Cubes has a 1985 ancestor, and it proves less than it is
//! tempting to claim.**
//!
//! Ticket: R-170. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p170
//! ```
//!
//! Writes `docs/experiments/p-170.csv`.
//!
//! # What was missing
//!
//! This is the phase's only desk-work row, and it is desk work because the gap
//! is bibliographic rather than numerical. The crate implements Plantinga &
//! Vegter's isotopy certificate at `validate/isotopy.rs:3`, marches Kuhn's /
//! Freudenthal's triangulation at `marching_tetrahedra/table.rs`, and has
//! measured the certificate twice — `M-378` (`P-62`) and `M-410` (`P-98`). Both
//! of those are *descendants*. The technique they descend from — following an
//! implicitly defined manifold through a simplicial subdivision, which is what
//! Marching Cubes does to a cube grid — was published in **1985**, two years
//! before Lorensen & Cline, and this repository has never cited it.
//!
//! The two counts below are the evidence for that sentence and they are
//! **re-counted from the working tree on every run** rather than transcribed,
//! because the registration's own preamble is a count and a count decays. It
//! already has: see *The preamble has drifted* below.
//!
//! # The resolution log
//!
//! `C1`'s vacuity control is the whole point of the row: *each DOI must be
//! resolved and its title recorded in the finding, not merely quoted from this
//! document*. Every DOI in the table below was resolved on **2026-08-30**
//! against **Crossref** (`api.crossref.org/works/<doi>`, the registration
//! agency, authoritative for title / volume / issue / pages) and independently
//! against the `home-still` corpus server (`paper_get`, which aggregates
//! Crossref + OpenAlex + Semantic Scholar). The titles, venues, volumes and
//! page ranges in `SOURCES` are transcribed from those responses and from no
//! other source.
//!
//! **Three resolution hazards were hit while doing it**, and they are recorded
//! because a verification that never met a failure mode has not been tested:
//!
//! - **The registered trap fires exactly as registered.** `10.1137/0722019` is
//!   *Differential Equations and the QR Algorithm* by T. Nanda, SIAM J. Numer.
//!   Anal. **22(2):310-321**, April 1985 — the article immediately *preceding*
//!   Allgower & Schmidt in the same issue. It is a row in the table, not a
//!   footnote, and `trap_gap_is_one` proves the adjacency by arithmetic.
//! - **A guessed DOI resolved to an unrelated paper.** Allgower & Georg's
//!   *Numerical Continuation Methods* is needed because it is what the
//!   corroborating source's quote actually points at. Guessing
//!   `10.1007/978-3-642-61257-2` from the Springer series and asking `paper_get`
//!   returned *Practical Ulam-Hyers-Rassias stability for nonlinear equations*
//!   by Fečkan & Wang — a merged aggregator record. The DOI is in fact correct;
//!   searching by title rather than trusting the merged record is what settled
//!   it. Recorded as `resolved_via = crossref_via_title_search`.
//! - **An aggregator reversed an author order.** `paper_get` on
//!   `10.1145/37402.37422` returns *Cline; Lorensen*; Crossref's `sequence`
//!   fields give **Lorensen** first. Author strings below follow Crossref.
//!
//! # What the sources actually prove, quoted rather than paraphrased
//!
//! **`C2`: AS85 proves a residual bound, not a homeomorphism.** Its abstract,
//! verbatim from Crossref: *"The method begins at a point `x₀` in the solution
//! set where the derivative `DH(x₀)` is of full rank. Given any `ε > 0`, a
//! piecewise linear manifold is constructed along which `‖H(x)‖_∞ < ε`."* That
//! is a full-rank hypothesis at the seed and a residual bound along the output.
//! Neither the title nor the abstract contains *homeomorphic*, *isotopic* or
//! *topologically equivalent*.
//!
//! **The full text of AS85 is paywalled and was not read.** `paper_download`
//! reports *"No open-access PDF found for DOI: 10.1137/0722020"*; Crossref lists
//! only an `epubs.siam.org` link marked `intended-application:
//! similarity-checking`. So `full_text_read = false` on that row and its
//! `evidence` column reads `abstract_verbatim_crossref`, not `full_text`. The
//! reading is of the abstract plus two corroborating primary sources, and the
//! row says so rather than implying a reading it does not have.
//!
//! **The corroboration was obtained and read in full.** Boissonnat &
//! Wintraecken, *The Topological Correctness of PL-Approximations of
//! Isomanifolds*, SoCG 2020, is open access at Dagstuhl and was read. Two
//! passages carry `C2`. Page 20:3, verbatim: *"The fact that the zero set of
//! `f_PL` is a manifold was proved (under strong condition) by Allgower and
//! Georg [3, Theorem 15.4.1], **without a homeomorphism with the zero set of
//! f**."* And page 20:2 on the whole older line: *"weaker results [2, 3] have
//! been known for a while, e.g. bounds on the one-sided Hausdorff distance, on
//! the approximation of tangent spaces, and manifoldness of the approximation
//! (under strong conditions)."*
//!
//! **One correction to the registration, which sharpens it rather than
//! weakening it.** The registration attributes that quote to AS85. It is
//! verbatim in B&W and it is about the Allgower line, but B&W's `[3]` is
//! **Allgower & Georg, *Numerical Continuation Methods*, Springer 1990** — their
//! `[1]`, `[2]` and `[3]` are all Allgower & Georg, and **AS85 is not in B&W's
//! bibliography at all**. So the quote is a fact about the textbook of the line,
//! not about AS85, and it is recorded on the `allgower_georg_1990_ncm` row where
//! it belongs. AS85's own non-homeomorphism status rests on its own abstract,
//! which is stronger evidence than a third party's remark about a different
//! paper. `C2` is therefore answered from two independent primary sources with
//! the attribution corrected.
//!
//! **`C3`: DLTW traces curves.** Abstract, verbatim: *"We present a method for
//! tracing a **curve** that is represented as the contour of a function in
//! Euclidean space of any dimension. … The algorithm does not fail in the
//! presence of high curvature of the contour; it accumulates essentially no
//! round-off error and has a well-defined integer test for detecting a loop."*
//! Robustness and loop detection. No topological-correctness claim, and the
//! object is a curve.
//!
//! # Two sources the registration did not name, and why they are rows anyway
//!
//! Both were found *inside* the reference lists of sources the registration did
//! name, which is the only reason they can be added without re-registering: they
//! are evidence for `C1`'s existing precedence clause, not a new clause.
//!
//! - **Allgower & Gnutzmann 1987** (`10.1137/0724033`), *An Algorithm for
//!   Piecewise Linear Approximation of Implicitly Defined Two-Dimensional
//!   Surfaces*, SIAM J. Numer. Anal. **24(2):452-469**, **April 1987** — four
//!   months before Marching Cubes, and it is the **two-dimensional surface**
//!   case, which is the case Marching Cubes solves. Its abstract again promises
//!   *"Error bounds on `‖H(x)‖_∞` … in terms of the mesh size of `T`"* from a
//!   *"regular point"*, so `C2`'s shape survives the specialisation to surfaces.
//!   Its Crossref reference `R2` is `10.1137/0722020`, publisher-asserted: the
//!   1985 → 1987 link is in the deposited metadata, not inferred.
//! - **Allgower & Schmidt 1984** (`10.1007/978-3-642-45567-4_26`), *Piecewise
//!   Linear Approximation of Solution Manifolds for Nonlinear Systems of
//!   Equations*, Lecture Notes in Econom. and Math. Systems **226**, pp.
//!   339-347, from the 8th International Conference on Operations Research,
//!   Karlsruhe 1983. Found as reference `R3` of Allgower & Gnutzmann 1987. It is
//!   **three years** before Lorensen & Cline and is the earliest member of the
//!   line this repository can name.
//!
//! And the descent into graphics is documented by the graphics paper itself:
//! DLTW's Crossref reference `e_1_2_1_2_2` reads *"ALLGOWER, E. L., AND SCHMIDT,
//! P.H. An algorithm for piecewise linear approximation of an implicitly defined
//! manifold. SIAM d. Numer. Anal. 22 (April 1985), 322-346."* — recorded per row
//! as `cited_by_dltw`.
//!
//! # The two counts, and exactly what corpus they are over
//!
//! `descendant_citations_in_repo` and `root_citations_in_repo` are **raw file
//! counts over one declared corpus**, global to the run and therefore identical
//! on every row (the header says so; the contract allows it). The corpus is
//! every file under the workspace root whose extension is in a fourteen-entry
//! text allowlist, minus four exclusions, each mechanical:
//!
//! - directories named `.git`, `target`, `node_modules`, `dist`, `.venv`,
//!   `__pycache__` — build and VCS artefacts.
//! - the `docs/experiments/` subtree — **generated output**, including this
//!   bench's own CSV. Counting it would make the measurement non-idempotent: a
//!   second run would find the titles this run wrote.
//! - the file named `experiment_p170.rs` — the scanner's own source, which
//!   mentions every pattern it searches for. `scanner_did_not_read_itself`
//!   asserts the exclusion took effect.
//! - nothing else. There is no judgement call in the corpus.
//!
//! A file counts once if any of its needles occurs as a lowercased substring.
//! The needles are declared, not derived: the descendant is `plantinga`; the
//! root is `allgower` **or** `10.1137/0722020`.
//!
//! **The root needle is `allgower`, not `schmidt`, and that is a measurement.**
//! Bare `schmidt` also matches `bevy_isomesh/examples/game_carve_seams.rs` and
//! `crates/isomesh/benches/experiment_p28.rs`, both of which say *Gram-Schmidt*.
//! `naive_schmidt_files` records the over-matching count and
//! `naive_pattern_over_matches` asserts it is strictly larger than the root
//! count, so the pattern choice is justified by the data rather than by taste.
//!
//! Two further cuts are reported as extras, because the registration's three
//! preamble numbers turn out not to have been counted the same way:
//!
//! - `*_excl_phase27_docs` removes the two `2026-08-29-phase-27-*.md` research
//!   documents — the corpus as it stood before this phase was written down.
//! - `*_excl_registration` additionally removes `FINDINGS.md`, `BACKLOG.md` and
//!   `crates/isomesh/src/experiment.rs` — the ledger-and-registration layer,
//!   which names a paper because an experiment *registered* it rather than
//!   because the crate *uses* it.
//!
//! **The corpus is the working tree at the stamped commit and it is growing
//! under this row.** Phase 27's fifty harnesses land into
//! `crates/isomesh/benches/` over the phase, and some of them cite
//! Plantinga-Vegter in their own headers — `experiment_p138.rs` already does —
//! so `descendant_citations_in_repo` **rises as the phase proceeds**, and the
//! figure in this CSV is a snapshot of that commit rather than a constant.
//! `phase27_harnesses_in_corpus` counts how many of `experiment_p127.rs` …
//! `experiment_p176.rs` were present, which is exactly the term that moves.
//!
//! No clause and no control depends on the value. `C1`, `C2` and `C3` are about
//! what the papers prove. The descendant count's only job is to be **non-zero**,
//! as the scanner's non-vacuity witness, and drift can only strengthen that. The
//! reading that carries the thesis — *"cites the parent theory in none"* — is
//! `root_citations_excl_registration`, and it is stable at **0** because the only
//! files in the tree that name the root are the registration layer itself.
//!
//! # The preamble has drifted, and only two thirds of it was ever right
//!
//! The registration's preamble says *"cites Plantinga-Vegter in 18 files and
//! Freudenthal/Kuhn in 6 including `marching_tetrahedra/table.rs`, and cites the
//! parent theory in none"*. Measured while authoring, on 2026-08-30, over a
//! 507-file corpus with none of Phase 27's harnesses yet landed — **the CSV is
//! authoritative and these are the figures that motivated the columns**:
//!
//! | claim | raw | excl. phase-27 docs | excl. registration | verdict |
//! |---|---|---|---|---|
//! | Plantinga-Vegter = 18 | 19 | **18** | 15 | right, under cut A |
//! | Freudenthal/Kuhn = 6 | 11 | 9 | 7 | **wrong under every cut** |
//! | parent theory = none | 4 | 2 | **0** | right, under cut B |
//!
//! So the descendant figure was exact for the corpus that existed before the
//! phase document itself, the root figure is exact once the registration layer
//! that introduced the name is removed — and the Kuhn figure is off by at least
//! three however it is cut. Narrowing the needle does not rescue it either:
//! `freudenthal` alone gives 9 raw / 8 excl-A. `marching_tetrahedra/table.rs`
//! *is* among the matches, as the registration says, and
//! `scanner_finds_the_named_file` asserts it — the named file is a known
//! positive and therefore the scanner's calibration.
//!
//! The three `doc_claim_*_reproduced` columns re-decide this on every run rather
//! than freezing the table above: `doc_claim_descendant_reproduced` compares cut
//! A against 18 and will read `false` once a Phase 27 harness cites the
//! descendant, which is drift rather than a defect and is why the raw figure is
//! reported beside it.
//!
//! None of this touches a clause verdict: `C1`, `C2` and `C3` are about what the
//! papers prove, not about the preamble's arithmetic. It is reported because a
//! preamble that has rotted and is not reported rots further.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `allgower_schmidt_1984` | earliest member of the line | no |
//! | `allgower_schmidt_1985` | **the root `C1` names** | no |
//! | `nanda_1985_trap` | **the registered DOI trap** | **yes** |
//! | `allgower_gnutzmann_1987` | the 2-surface case; April 1987 | no |
//! | `lorensen_cline_1987` | Marching Cubes; the precedence baseline | **yes** |
//! | `dobkin_levy_thurston_wilks_1990` | **`C3`'s subject** | no |
//! | `allgower_georg_1990_ncm` | the true referent of B&W's quote | no |
//! | `boissonnat_wintraecken_2020` | the isotopy result; the only `true` | no |
//!
//! Two controls. `nanda_1985_trap` proves the resolver discriminates: a
//! one-digit slip in the DOI returns a paper on the QR algorithm, and if it did
//! not, `doi_verified` would be measuring nothing. `lorensen_cline_1987` turns
//! *"predates"* from a word into `1985 < 1987`.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's `SHARE` line reads *"none — this is a citation repair and
//! a scope statement"*, and that is discharged rather than skipped: nothing here
//! moves a stage, changes a byte of `crates/isomesh/src/**` or proposes to. The
//! deliverable is the citation table and the scope sentence — that the crate's
//! PL-approximation technique has a 1984/1985 root which proves a residual bound
//! under a regularity hypothesis, and that the topological guarantee the crate
//! actually relies on dates from 2004 (Plantinga-Vegter) and 2020 (B&W), not
//! from the root.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record` and every panic begins `VOID: `.
//!
//! - **`every_doi_resolved`** — every row has a non-empty resolved title and a
//!   non-empty `resolved_via`, so `doi_verified` is true on all eight. Column:
//!   `doi_verified`.
//! - **`titles_are_distinct`** — eight rows, eight different titles. A
//!   copy-pasted title would make the table look resolved while carrying one
//!   answer eight times. Column: `title`.
//! - **`trap_gap_is_one`** — the two SIAM article numbers parse to integers
//!   differing by exactly `1`, on the same venue / volume / issue / year, with
//!   different titles. The trap is *computed*, not asserted. Columns: `doi`,
//!   `title`.
//! - **`homeomorphism_column_is_not_constant`** — at least one row `true` and at
//!   least one `false`. Eight `false` rows would be a column that could not have
//!   been true, which is not a measurement (`M-44`). Column:
//!   `is_homeomorphism`.
//! - **`precedence_is_arithmetic`** — `1984 < 1985 < 1987`, read off the `year`
//!   column rather than from the prose. Column: `year`.
//! - **`scanner_finds_the_named_file`** — `marching_tetrahedra/table.rs` is
//!   among the Kuhn matches and the descendant count is non-zero, both from the
//!   same walk that reports the root count. A zero root count from a scanner
//!   that finds nothing is not a citation gap. Columns:
//!   `descendant_citations_in_repo`, `kuhn_citations_in_repo`.
//! - **`naive_pattern_over_matches`** — bare `schmidt` matches strictly more
//!   files than the root needle, proving the discarded pattern really was
//!   unusable. Column: `naive_schmidt_files`.
//! - **`scanner_did_not_read_itself`** — no scanned path is
//!   `experiment_p170.rs` and none is under `docs/experiments/`, so the two
//!   idempotence exclusions took effect. Column: `files_scanned`.

#![allow(clippy::too_many_lines)]

mod common;

use std::path::{Path, PathBuf};

/// One resolved bibliographic record.
///
/// Every string field is transcribed from a Crossref or `home-still` response,
/// and every one is CSV-safe by construction: **no field may contain a comma**,
/// because `common::experiment::Run::record` refuses one rather than quoting.
/// Author lists are joined with `|`; clauses inside a `guarantee` are separated
/// with `;`.
struct Source {
    /// The `source` column: a stable snake_case row id.
    source: &'static str,
    year: u32,
    doi: &'static str,
    /// The resolved title. Empty would mean unresolved; none are.
    title: &'static str,
    authors: &'static str,
    venue: &'static str,
    /// `volume(issue)`, or just the volume where there is no issue.
    volume_issue: &'static str,
    pages: &'static str,
    claim_type: &'static str,
    /// The hypothesis the *paper* requires, in the paper's own terms.
    hypothesis: &'static str,
    /// The guarantee the *paper* delivers, in the paper's own terms.
    guarantee: &'static str,
    /// Does this work prove a homeomorphism (or stronger) with the true zero
    /// set? `false` for every residual bound. `true` only for B&W 2020, which
    /// proves the strictly stronger isotopy.
    is_homeomorphism: bool,
    /// Which lookups produced the record above.
    resolved_via: &'static str,
    /// Where the `guarantee` sentence was read from.
    evidence: &'static str,
    /// Was the full text obtained and read, as opposed to the abstract?
    full_text_read: bool,
    /// Is this work in DLTW 1990's own reference list?
    cited_by_dltw: bool,
    /// Lowercased needles used to count this work's mentions in the repository.
    needles: &'static [&'static str],
}

impl Source {
    /// A DOI is verified when a lookup returned a title for it.
    fn doi_verified(&self) -> bool {
        !self.title.is_empty() && !self.resolved_via.is_empty()
    }
}

/// The eight resolved sources, oldest first.
///
/// `nanda_1985_trap` is the registered trap and is a row rather than a remark:
/// recording the wrong paper's real title is what demonstrates that the right
/// paper's title was not merely copied out of the research document.
const SOURCES: &[Source] = &[
    Source {
        source: "allgower_schmidt_1984",
        year: 1984,
        doi: "10.1007/978-3-642-45567-4_26",
        title: "Piecewise Linear Approximation of Solution Manifolds for Nonlinear Systems of Equations",
        authors: "Eugene L. Allgower|Phillip H. Schmidt",
        venue: "Lecture Notes in Econom. and Math. Systems",
        volume_issue: "226",
        pages: "339-347",
        claim_type: "residual_bound",
        hypothesis: "a nonlinear system H(x) = 0 with a known starting point on the solution set",
        guarantee: "a piecewise linear approximation of the solution manifold by simplicial subdivision",
        is_homeomorphism: false,
        resolved_via: "home_still_paper_get|crossref_ref_R3_of_10.1137/0724033",
        evidence: "title_and_venue_only_no_abstract_deposited",
        full_text_read: false,
        cited_by_dltw: false,
        needles: &["10.1007/978-3-642-45567-4_26"],
    },
    Source {
        source: "allgower_schmidt_1985",
        year: 1985,
        doi: "10.1137/0722020",
        title: "An Algorithm for Piecewise-Linear Approximation of an Implicitly Defined Manifold",
        authors: "Eugene L. Allgower|Phillip H. Schmidt",
        venue: "SIAM J. Numer. Anal.",
        volume_issue: "22(2)",
        pages: "322-346",
        claim_type: "residual_bound",
        hypothesis: "H: R^(N+K) -> R^N smooth; the seed x0 lies in H^-1(0) and DH(x0) has full rank",
        guarantee: "for any eps > 0 a PL manifold is constructed along which ||H(x)||_inf < eps",
        is_homeomorphism: false,
        resolved_via: "crossref_api|home_still_paper_get",
        evidence: "abstract_verbatim_crossref",
        full_text_read: false,
        cited_by_dltw: true,
        needles: &["allgower", "10.1137/0722020"],
    },
    Source {
        source: "nanda_1985_trap",
        year: 1985,
        doi: "10.1137/0722019",
        title: "Differential Equations and the $QR$ Algorithm",
        authors: "T. Nanda",
        venue: "SIAM J. Numer. Anal.",
        volume_issue: "22(2)",
        pages: "310-321",
        claim_type: "unrelated_work_trap",
        hypothesis: "isospectral flows on n x n matrices arising from Lax pairs",
        guarantee: "asymptotics of those flows give a new method for the symmetric eigenvalue problem",
        is_homeomorphism: false,
        resolved_via: "crossref_api|home_still_paper_get",
        evidence: "abstract_verbatim_crossref",
        full_text_read: false,
        cited_by_dltw: false,
        needles: &["10.1137/0722019"],
    },
    Source {
        source: "allgower_gnutzmann_1987",
        year: 1987,
        doi: "10.1137/0724033",
        title: "An Algorithm for Piecewise Linear Approximation of Implicitly Defined Two-Dimensional Surfaces",
        authors: "Eugene L. Allgower|Stefan Gnutzmann",
        venue: "SIAM J. Numer. Anal.",
        volume_issue: "24(2)",
        pages: "452-469",
        claim_type: "residual_bound_2d_surfaces",
        hypothesis: "H: R^(N+2) -> R^N smooth and x0 in H^-1(0) a regular point of H",
        guarantee: "error bounds on ||H(x)||_inf over M_T in terms of the mesh size of T",
        is_homeomorphism: false,
        resolved_via: "crossref_api|home_still_paper_search",
        evidence: "abstract_verbatim_crossref",
        full_text_read: false,
        cited_by_dltw: false,
        needles: &["gnutzmann", "10.1137/0724033"],
    },
    Source {
        source: "lorensen_cline_1987",
        year: 1987,
        doi: "10.1145/37402.37422",
        title: "Marching cubes: A high resolution 3D surface construction algorithm",
        authors: "William E. Lorensen|Harvey E. Cline",
        venue: "SIGGRAPH Comput. Graph.",
        volume_issue: "21(4)",
        pages: "163-169",
        claim_type: "no_correctness_claim",
        hypothesis: "none stated",
        guarantee: "triangle models of constant-density surfaces from a case table plus linear interpolation; no topological guarantee stated",
        is_homeomorphism: false,
        resolved_via: "crossref_api|home_still_paper_get",
        evidence: "abstract_verbatim_crossref",
        full_text_read: false,
        cited_by_dltw: false,
        needles: &["lorensen", "10.1145/37402"],
    },
    Source {
        source: "dobkin_levy_thurston_wilks_1990",
        year: 1990,
        doi: "10.1145/88560.88575",
        title: "Contour tracing by piecewise linear approximations",
        authors: "David P. Dobkin|Allan R. Wilks|Silvio V. F. Levy|William P. Thurston",
        venue: "ACM Trans. Graph.",
        volume_issue: "9(4)",
        pages: "389-423",
        claim_type: "robustness_and_loop_detection",
        hypothesis: "a curve given as the contour of a function in Euclidean space of any dimension; a triangulation generated by reflections",
        guarantee: "traces a curve without failing at high curvature; accumulates essentially no round-off error; has a well-defined integer test for detecting a loop",
        is_homeomorphism: false,
        resolved_via: "crossref_api|home_still_paper_get",
        evidence: "abstract_verbatim_crossref",
        full_text_read: false,
        cited_by_dltw: false,
        needles: &["dobkin", "10.1145/88560"],
    },
    Source {
        source: "allgower_georg_1990_ncm",
        year: 1990,
        doi: "10.1007/978-3-642-61257-2",
        title: "Numerical Continuation Methods",
        authors: "Eugene L. Allgower|Kurt Georg",
        venue: "Springer Series in Computational Mathematics",
        volume_issue: "13",
        pages: "n/a-book",
        claim_type: "manifoldness_under_strong_conditions",
        hypothesis: "strong conditions incl. the zero set avoiding simplices of dimension below the codimension",
        guarantee: "Thm 15.4.1: the zero set of f_PL is a manifold -- explicitly without a homeomorphism with the zero set of f",
        is_homeomorphism: false,
        resolved_via: "crossref_via_title_search|guessed_doi_returned_a_merged_record",
        evidence: "boissonnat_wintraecken_socg_2020_p20_3_verbatim",
        full_text_read: false,
        cited_by_dltw: false,
        needles: &["10.1007/978-3-642-61257-2"],
    },
    Source {
        source: "boissonnat_wintraecken_2020",
        year: 2020,
        doi: "10.4230/LIPIcs.SoCG.2020.20",
        title: "The Topological Correctness of PL-Approximations of Isomanifolds",
        authors: "Jean-Daniel Boissonnat|Mathijs Wintraecken",
        venue: "LIPIcs (SoCG 2020)",
        volume_issue: "164",
        pages: "20:1-20:18",
        claim_type: "isotopy",
        hypothesis: "0 a regular value of f; longest edge D small enough at bounded simplex thickness T",
        guarantee: "Thm 20: the zero set of f_PL is a manifold isotopic to the zero set of f; plus a Frechet bound",
        is_homeomorphism: true,
        resolved_via: "home_still_paper_get|dagstuhl_open_access_pdf",
        evidence: "full_text_read_p20_1_to_20_18",
        full_text_read: true,
        cited_by_dltw: false,
        needles: &["wintraecken", "10.4230"],
    },
];

/// The descendant this repository does cite: Plantinga & Vegter.
const DESCENDANT_NEEDLES: &[&str] = &["plantinga"];
/// The 1985 root, matched by author surname or by DOI.
const ROOT_NEEDLES: &[&str] = &["allgower", "10.1137/0722020"];
/// The triangulation the registration also counts.
const KUHN_NEEDLES: &[&str] = &["freudenthal", "kuhn"];
/// The pattern deliberately **not** used for the root, kept as a control.
const NAIVE_ROOT_NEEDLES: &[&str] = &["schmidt"];

/// The file this bench lives in. Excluded from its own walk.
const SELF_FILE: &str = "experiment_p170.rs";
/// Generated output. Excluded so the count does not change on a second run.
const GENERATED_SUBTREE: &str = "docs/experiments";

/// Cut A: the corpus as it stood before this phase was written down.
const PHASE_27_DOCS: &[&str] = &[
    "docs/research/2026-08-29-phase-27-fifty-experiments-from-unmined-mathematics.md",
    "docs/research/2026-08-29-phase-27-axes-and-vocabulary-v2.md",
];
/// Cut B: cut A plus the ledger-and-registration layer.
const REGISTRATION_LAYER: &[&str] = &[
    "docs/research/2026-08-29-phase-27-fifty-experiments-from-unmined-mathematics.md",
    "docs/research/2026-08-29-phase-27-axes-and-vocabulary-v2.md",
    "FINDINGS.md",
    "BACKLOG.md",
    "crates/isomesh/src/experiment.rs",
];

/// The registration's own preamble figures, for the drift table.
const DOC_CLAIM_DESCENDANT: usize = 18;
const DOC_CLAIM_KUHN: usize = 6;
const DOC_CLAIM_ROOT: usize = 0;

/// One file of the scanned corpus, held lowercased for substring matching.
struct Scanned {
    /// Path relative to the workspace root, with `/` separators.
    rel: String,
    lowercased: String,
}

/// Is this a text file the corpus includes?
///
/// An allowlist rather than a denylist: a new binary format appearing in the
/// tree must not silently join the corpus and move a published count.
fn is_text(name: &str) -> bool {
    const EXTS: &[&str] = &[
        "rs", "md", "toml", "py", "sh", "wgsl", "html", "css", "js", "json", "yml", "yaml", "txt",
        "csv",
    ];
    name.rsplit_once('.')
        .is_some_and(|(_, ext)| EXTS.contains(&ext.to_ascii_lowercase().as_str()))
}

/// Walk the workspace and read every text file of the corpus.
///
/// Symlinked directories are not followed: `DirEntry::file_type` does not
/// dereference, so a link to a directory is neither descended nor read. A file
/// that is not valid UTF-8 is skipped, which cannot hide a citation because
/// every document in this repository is UTF-8 text.
fn scan(root: &Path) -> Vec<Scanned> {
    const SKIP_DIRS: &[&str] = &[
        ".git",
        "target",
        "node_modules",
        "dist",
        ".venv",
        "__pycache__",
    ];

    let mut out = Vec::new();
    let mut stack = vec![root.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = std::fs::read_dir(&dir) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
                continue;
            };
            let Ok(kind) = entry.file_type() else {
                continue;
            };
            if kind.is_dir() {
                if !SKIP_DIRS.contains(&name) {
                    stack.push(path);
                }
                continue;
            }
            if !kind.is_file() || !is_text(name) || name == SELF_FILE {
                continue;
            }
            let Some(rel) = path
                .strip_prefix(root)
                .ok()
                .and_then(Path::to_str)
                .map(str::to_owned)
            else {
                continue;
            };
            if rel.starts_with(GENERATED_SUBTREE) {
                continue;
            }
            let Ok(text) = std::fs::read_to_string(&path) else {
                continue;
            };
            out.push(Scanned {
                lowercased: text.to_lowercase(),
                rel,
            });
        }
    }
    out.sort_unstable_by(|a, b| a.rel.cmp(&b.rel));
    out
}

/// How many files of `corpus`, excluding `excluded`, mention any needle?
fn count(corpus: &[Scanned], needles: &[&str], excluded: &[&str]) -> usize {
    corpus
        .iter()
        .filter(|f| !excluded.contains(&f.rel.as_str()))
        .filter(|f| needles.iter().any(|n| f.lowercased.contains(n)))
        .count()
}

/// Which files of `corpus` mention any needle?
fn matching<'a>(corpus: &'a [Scanned], needles: &[&str]) -> Vec<&'a str> {
    corpus
        .iter()
        .filter(|f| needles.iter().any(|n| f.lowercased.contains(n)))
        .map(|f| f.rel.as_str())
        .collect()
}

/// How many of Phase 27's own harnesses — `experiment_p127.rs` through
/// `experiment_p176.rs` — are in the corpus?
///
/// This is the term that makes `descendant_citations_in_repo` move while the
/// phase is mid-flight: a sibling harness citing Plantinga-Vegter in its header
/// joins the corpus between one bench run and the next. Reporting it is what
/// turns that drift from an unexplained number into an accounted one.
fn phase27_harnesses(corpus: &[Scanned]) -> usize {
    corpus
        .iter()
        .filter(|f| {
            f.rel
                .strip_prefix("crates/isomesh/benches/experiment_p")
                .and_then(|rest| rest.strip_suffix(".rs"))
                .and_then(|id| id.parse::<u32>().ok())
                .is_some_and(|id| (127..=176).contains(&id))
        })
        .count()
}

/// The article number of a legacy SIAM DOI, e.g. `722020` for
/// `10.1137/0722020`.
///
/// Legacy SIAM suffixes encode journal, volume, issue and article position, so
/// two consecutive articles in one issue differ by one. That is what makes
/// `10.1137/0722019` a plausible slip and therefore a trap.
fn siam_article_number(doi: &str) -> u32 {
    doi.rsplit_once('/')
        .expect("a legacy SIAM DOI has exactly one slash")
        .1
        .parse()
        .expect("a legacy SIAM DOI suffix is all digits")
}

/// Look a row up by its `source` id.
fn row(source: &str) -> &'static Source {
    SOURCES
        .iter()
        .find(|s| s.source == source)
        .expect("every row referenced by a control is in SOURCES")
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-170");

    common::experiment::run(prereg, |run| {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .canonicalize()
            .expect("the workspace root two levels above the crate resolves");
        let corpus = scan(&root);

        let descendant = count(&corpus, DESCENDANT_NEEDLES, &[]);
        let descendant_a = count(&corpus, DESCENDANT_NEEDLES, PHASE_27_DOCS);
        let descendant_b = count(&corpus, DESCENDANT_NEEDLES, REGISTRATION_LAYER);
        let root_cites = count(&corpus, ROOT_NEEDLES, &[]);
        let root_a = count(&corpus, ROOT_NEEDLES, PHASE_27_DOCS);
        let root_b = count(&corpus, ROOT_NEEDLES, REGISTRATION_LAYER);
        let kuhn = count(&corpus, KUHN_NEEDLES, &[]);
        let kuhn_a = count(&corpus, KUHN_NEEDLES, PHASE_27_DOCS);
        let kuhn_b = count(&corpus, KUHN_NEEDLES, REGISTRATION_LAYER);
        let naive = count(&corpus, NAIVE_ROOT_NEEDLES, &[]);
        let harnesses = phase27_harnesses(&corpus);

        // ── vacuity controls, all of them before the first row ──

        for source in SOURCES {
            assert!(
                source.doi_verified(),
                "VOID: {} carries no resolved title, so `doi_verified` would record a DOI \
                 quoted from the research document rather than resolved against a registration \
                 agency, which is exactly what this row's vacuity control forbids",
                source.source
            );
        }

        for (i, a) in SOURCES.iter().enumerate() {
            for b in &SOURCES[i + 1..] {
                assert!(
                    a.title != b.title,
                    "VOID: `{}` and `{}` carry the same resolved title, so the `title` column \
                     holds one answer twice and the table only looks resolved",
                    a.source,
                    b.source
                );
            }
        }

        let as85 = row("allgower_schmidt_1985");
        let trap = row("nanda_1985_trap");
        let as84 = row("allgower_schmidt_1984");
        let mc = row("lorensen_cline_1987");
        let dltw = row("dobkin_levy_thurston_wilks_1990");

        let trap_gap = siam_article_number(as85.doi) - siam_article_number(trap.doi);
        assert!(
            trap_gap == 1
                && as85.venue == trap.venue
                && as85.volume_issue == trap.volume_issue
                && as85.year == trap.year
                && as85.title != trap.title,
            "VOID: the registered trap is not adjacent to the root -- gap {trap_gap}, venues \
             {}/{}, issues {}/{} -- so `10.1137/0722019` was never a plausible slip and \
             demonstrating it proves nothing about the resolution",
            as85.venue,
            trap.venue,
            as85.volume_issue,
            trap.volume_issue
        );

        let homeomorphisms = SOURCES.iter().filter(|s| s.is_homeomorphism).count();
        assert!(
            homeomorphisms > 0 && homeomorphisms < SOURCES.len(),
            "VOID: `is_homeomorphism` reads the same on all {} rows ({homeomorphisms} true), so \
             C2's false is a value that could not have been true and is not a measurement (M-44)",
            SOURCES.len()
        );

        assert!(
            as84.year < as85.year && as85.year < mc.year,
            "VOID: the precedence claim is not arithmetic on the recorded years -- {} then {} \
             then {} -- so `predates Lorensen & Cline` would be prose rather than a comparison",
            as84.year,
            as85.year,
            mc.year
        );

        let kuhn_files = matching(&corpus, KUHN_NEEDLES);
        let named = "crates/isomesh/src/marching_tetrahedra/table.rs";
        assert!(
            descendant > 0 && kuhn_files.contains(&named),
            "VOID: the scanner reports {descendant} descendant files and {} Kuhn files without \
             finding `{named}`, which the registration names as a known positive -- so a low \
             root count would be a broken walk rather than a citation gap",
            kuhn_files.len()
        );

        assert!(
            naive > root_cites,
            "VOID: the discarded `schmidt` needle matches {naive} files against the root \
             needle's {root_cites}, so it does not in fact over-match and choosing `allgower` \
             over it was taste rather than measurement"
        );

        assert!(
            corpus.len() > 300
                && !corpus
                    .iter()
                    .any(|f| f.rel.ends_with(SELF_FILE) || f.rel.starts_with(GENERATED_SUBTREE)),
            "VOID: the walk scanned {} files and did not exclude its own source or the generated \
             `{GENERATED_SUBTREE}` subtree, so the counts include this bench counting itself and \
             would change on a second run",
            corpus.len()
        );

        // ── the clause verdicts, global to the run and identical on every row ──

        // C1: the citation resolves as stated, and the trap resolves elsewhere.
        let c1 = as85.doi_verified()
            && as85.year == 1985
            && as85.venue == "SIAM J. Numer. Anal."
            && as85.volume_issue == "22(2)"
            && as85.pages == "322-346"
            && as85.authors == "Eugene L. Allgower|Phillip H. Schmidt"
            && as85.title
                == "An Algorithm for Piecewise-Linear Approximation of an Implicitly Defined Manifold"
            && trap.doi_verified()
            && trap.authors == "T. Nanda"
            && trap.claim_type == "unrelated_work_trap"
            && as85.year < mc.year;

        // C2: AS85 proves a residual bound under a regularity hypothesis, and
        // not a homeomorphism. Falsified only by AS85 proving one.
        let c2 = !as85.is_homeomorphism
            && as85.claim_type == "residual_bound"
            && as85.hypothesis.contains("full rank")
            && as85.guarantee.contains("||H(x)||_inf < eps")
            && homeomorphisms > 0;

        // C3: DLTW traces curves and claims robustness plus loop detection.
        // Falsified by DLTW covering surfaces.
        let c3 = !dltw.is_homeomorphism
            && dltw.claim_type == "robustness_and_loop_detection"
            && dltw.guarantee.contains("curve")
            && !dltw.guarantee.contains("surface")
            && dltw.hypothesis.contains("curve");

        for source in SOURCES {
            let mine = count(&corpus, source.needles, &[]);
            run.record(&[
                ("source", source.source.to_string()),
                ("year", source.year.to_string()),
                ("doi", source.doi.to_string()),
                ("doi_verified", source.doi_verified().to_string()),
                ("claim_type", source.claim_type.to_string()),
                ("hypothesis", source.hypothesis.to_string()),
                ("guarantee", source.guarantee.to_string()),
                ("is_homeomorphism", source.is_homeomorphism.to_string()),
                ("descendant_citations_in_repo", descendant.to_string()),
                ("root_citations_in_repo", root_cites.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // ── extras (M-273) ──
                ("title", source.title.to_string()),
                ("authors", source.authors.to_string()),
                ("venue", source.venue.to_string()),
                ("volume_issue", source.volume_issue.to_string()),
                ("pages", source.pages.to_string()),
                ("resolved_via", source.resolved_via.to_string()),
                ("evidence", source.evidence.to_string()),
                ("full_text_read", source.full_text_read.to_string()),
                ("cited_by_dltw", source.cited_by_dltw.to_string()),
                ("this_source_citations_in_repo", mine.to_string()),
                (
                    "descendant_citations_excl_phase27_docs",
                    descendant_a.to_string(),
                ),
                (
                    "descendant_citations_excl_registration",
                    descendant_b.to_string(),
                ),
                ("root_citations_excl_phase27_docs", root_a.to_string()),
                ("root_citations_excl_registration", root_b.to_string()),
                ("kuhn_citations_in_repo", kuhn.to_string()),
                ("kuhn_citations_excl_phase27_docs", kuhn_a.to_string()),
                ("kuhn_citations_excl_registration", kuhn_b.to_string()),
                ("naive_schmidt_files", naive.to_string()),
                ("doc_claim_descendant", DOC_CLAIM_DESCENDANT.to_string()),
                ("doc_claim_kuhn", DOC_CLAIM_KUHN.to_string()),
                ("doc_claim_root", DOC_CLAIM_ROOT.to_string()),
                (
                    "doc_claim_descendant_reproduced",
                    (descendant_a == DOC_CLAIM_DESCENDANT).to_string(),
                ),
                (
                    "doc_claim_kuhn_reproduced",
                    (kuhn == DOC_CLAIM_KUHN
                        || kuhn_a == DOC_CLAIM_KUHN
                        || kuhn_b == DOC_CLAIM_KUHN)
                        .to_string(),
                ),
                (
                    "doc_claim_root_reproduced",
                    (root_b == DOC_CLAIM_ROOT).to_string(),
                ),
                ("files_scanned", corpus.len().to_string()),
                ("phase27_harnesses_in_corpus", harnesses.to_string()),
                ("trap_gap", trap_gap.to_string()),
                (
                    "predates_marching_cubes",
                    (source.year < mc.year).to_string(),
                ),
            ]);
        }
    });
}
