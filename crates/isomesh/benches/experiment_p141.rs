//! **P-141 — deriving the case table by quantifier elimination, offline, once.**
//!
//! Ticket: R-141. Pre-registered before this harness existed.
//!
//! ```bash
//! Z3=~/.venvs/isomesh/bin/z3 cargo bench --bench experiment_p141
//! ```
//!
//! Writes `docs/experiments/p-141.csv`.
//!
//! **This is the one row in Phase 27 that needs an external binary**, and it is
//! not optional. The solver is the executable named by the `Z3` environment
//! variable, defaulting to `z3`. If it cannot be executed the harness panics
//! naming the variable and writes nothing: a CSV row for a solver that did not
//! run is the one artefact this row must never produce.
//!
//! # What was missing
//!
//! `docs/research/2026-08-12-axes-and-vocabulary.md` and the discovery dossier
//! dismiss cylindrical algebraic decomposition as *"global, doubly exponential,
//! needs arbitrary precision"*. All three are true and all three price the wrong
//! thing. A case table is computed **once in a project's lifetime**; the cost
//! this repository actually pays, over and over, is transcription. `M-207`,
//! `M-219`, `M-221`, `✗22`, `M-228` and `M-231` are six separate entries whose
//! subject is a symbol copied wrong — `M-219` is a reference implementation
//! storing `u_e1` where it compares `u_e2`, one character.
//!
//! The crate already refuses to transcribe the 256-case triangulation:
//! `table::CASES` is `build_cases()` evaluated during compilation
//! (marching_cubes/table.rs:180-194), and `AMBIGUOUS_FACES` is
//! `build_ambiguous_faces()` (:202-231). But *derived at compile time* is not
//! *proved*. Both are built from the same eight-corner sign convention the rest
//! of the module assumes, so what they establish is that the module agrees with
//! itself. Nothing in the tree compares either table against the object it is a
//! table **of** — the trilinear interpolant — over anything but finite samples:
//!
//! - `validate_table()` (marching_cubes/mod.rs:836) and
//!   `validate_decider_table()` (:857) check internal consistency — triangles on
//!   uncut edges, cut-edge mismatch, degeneracy. They never evaluate the
//!   interpolant.
//! - `M-40`'s ambiguous-face census is eight reference fields on 33³ grids.
//! - `M-165`'s 1,966-of-15,625 disagreement rate is a lattice of quantised
//!   corner values.
//! - `M-220` reports the singular face as `0 of 1,838` on eight fields and
//!   `0 of 299,215` over 400,000 random cells — a census that found the
//!   population empty.
//!
//! Every one of those is a sample. A statement about *every real corner
//! assignment* has never been made here, and it is exactly the statement a real
//! quantifier-elimination engine makes.
//!
//! **`M-220` is also the warning.** The tangency case — `S = 0` exactly, where
//! the bilinear's two hyperbola branches degenerate into crossing straight lines
//! — is measure zero, so a grid sweep is blind to it by construction. This
//! harness measured that directly while it was being written: an independent
//! 2,001-point sweep in `p` disagreed with the exact interval computation on
//! **162 of 98,304** sample assignments, and all 162 were exact zero-width
//! intersections, all in the same direction, every one a tangency. That is the
//! whole argument for a solver rather than a denser sweep, and it is why
//! `control_tangency_samples` is a vacuity control below rather than a curiosity.
//!
//! # The query
//!
//! Over eight `Real` corner parameters `f0..f7`, with the case's eight sign
//! constraints asserted (`fi < 0` where `cube::corner_inside` puts corner `i`
//! inside, `fi >= 0` otherwise — an exact zero is *outside*, cube.rs:171-173).
//!
//! On the face `(axis, side)`, `face_corners` (cube.rs:88-101) gives the four
//! corners counter-clockwise from outside; call their values `v0..v3` and place
//! them at `(0,0)`, `(1,0)`, `(1,1)`, `(0,1)` in a local `(p, q)` frame. The
//! trilinear interpolant restricted to that face is the bilinear
//!
//! ```text
//! B(p,q) = v0(1-p)(1-q) + v1 p(1-q) + v2 p q + v3 (1-p) q
//! ```
//!
//! and the topological question — *are the face's two inside corners joined
//! through the face?* — is stated without naming a formula:
//!
//! ```text
//! truth  <=>  for all p in [0,1] there is a q in [0,1] with B(p,q) < 0
//! ```
//!
//! **Why that is the connectivity statement, in both directions.** If the two
//! inside corners are joined, a path inside the negative region runs from
//! `(0,0)` to `(1,1)`; its `p`-coordinate is continuous and spans `[0,1]`, so by
//! the intermediate value theorem every `p` carries a negative point. Conversely
//! every fibre being non-empty forces one negative region rather than two, since
//! `B(p,·)` is affine and its negative set is therefore a single sub-interval,
//! continuous in `p`. Neither direction mentions a saddle, a diagonal or a
//! product — which is the point: the solver is not handed the answer.
//!
//! # Arms
//!
//! Two **encodings** of `truth`, and three **claims** to test it against. The
//! encodings are not two paths to one answer; the difference between them is
//! itself a measurement, and it is the one C1 is about.
//!
//! - `nested` states `truth` verbatim: `forall p . exists q . B(p,q) < 0`.
//! - `affine` discharges the inner existential analytically and nothing else.
//!   `B(p,·)` is **affine in `q`**, so its minimum over `[0,1]` is at an
//!   endpoint and `exists q in [0,1] . B(p,q) < 0` is exactly
//!   `B(p,0) < 0 or B(p,1) < 0`. `q = 0` is the ring edge `(v0,v1)` and `q = 1`
//!   is `(v3,v2)`, so the residual formula is a `forall` over one variable of a
//!   disjunction of two affine constraints with parametric coefficients. It is
//!   still nonlinear in the parameters and still needs real quantifier
//!   elimination; what it does not need is a model of a nested quantifier.
//!
//! | arm | encoding | claim under test | entries | is_control |
//! |---|---|---|---|---|
//! | `ambiguity_affine` | affine | `truth` is *constant* on the case's sign class — the claim `AMBIGUOUS_FACES` makes by leaving a bit clear | 1,536 = 256 × 6 | on the 1,344 unambiguous pairs |
//! | `decider_affine` | affine | `ambiguity::face_is_joined` is `truth` | 192 | no |
//! | `separate_affine` | affine | Marching Cubes proper's pairing is `truth` | 192 | yes |
//! | `decider_nested` | nested | `ambiguity::face_is_joined` is `truth` | 192 | no |
//! | `separate_nested` | nested | Marching Cubes proper's pairing is `truth` | 192 | yes |
//!
//! Every query has one shape — `(assert (not (= <claim> truth))) (check-sat)` —
//! so `unsat` means *the claim holds for every real corner assignment in this
//! case's sign class* and `sat` means a counterexample exists. There is no
//! second assertion form and no per-arm special case.
//!
//! `ambiguity_affine` is the derivation proper: `unsat` says `truth` does not
//! depend on the magnitudes, which **is** the definition of a non-ambiguous
//! face, so the arm re-derives `AMBIGUOUS_FACES` entry by entry by quantifier
//! elimination and compares it against the committed table. `separate_*` is the
//! positive control: `FaceAmbiguity::Separate` is the crate's default and
//! `CASES` at mask zero is its table, and the crate's own documentation says it
//! *"agrees with the field's own bilinear interpolant on the face only by luck"*
//! (ambiguity.rs:13-14). An arm that cannot come back `sat` has not shown the
//! instrument able to report bad news.
//!
//! # Budget, and what non-termination means
//!
//! `QUERY_TIMEOUT_S` is handed to the solver as its own `-T:<seconds>`, so no
//! query can hang and the sweep's worst case is bounded by construction.
//! `TOTAL_BUDGET` is a second, whole-run bound that **panics** rather than
//! truncating: a partial derivation is not a derivation, and a short CSV that
//! looks complete is worse than no CSV.
//!
//! Non-termination inside the budget is a **pre-registered outcome**, recorded
//! as `terminated=false` with the elapsed time, and `c1_holds` is that column —
//! C1 *is* the termination claim, so the two are the same predicate by
//! construction and the number worth reading is the aggregate.
//!
//! ## The asymmetry, measured before the sweep was designed
//!
//! The two encodings are in the file because of a fact worth stating in advance,
//! since it decided the design: **z3 refutes the nested `forall`/`exists` form
//! and cannot produce a model of it.** On the same entry, `(assert (not (=
//! decider truth))) (check-sat)` returns `unsat` in 59 ms, while flipping the
//! assertion to something satisfiable returns `timeout` at a 60 s cap, and the
//! `qsat` tactic likewise at 20 s. Refuting a `exists`/`forall`/`exists` prefix
//! reduces to a universal statement that nlsat handles; exhibiting a model of a
//! quantified nonlinear formula needs the elimination itself, and `nlqsat` on the
//! bare `truth` predicate did not finish in 90 s either.
//!
//! So a design that used only the nested form could never return `sat`, and an
//! instrument that can only produce one answer has not measured anything (`M-44`).
//! The affine encoding decides both directions, so the derivation runs on it and
//! the nested form is swept for what it is: C1's evidence, terminating on the
//! arms whose answer is a refutation and reported as **not terminating** on the
//! arm whose answer is a witness. Both are registered outcomes.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE line is *"none at runtime; this changes how the
//! table is produced, not what it costs"*, and that is discharged rather than
//! quoted: this harness calls no extractor, allocates no grid and touches
//! nothing on the extraction path. `CASES` and `AMBIGUOUS_FACES` are read, never
//! rebuilt. The cost moved here is a developer's, paid once, and the stage it
//! could move is the one `M-207`/`M-219`/`M-221` bill for — transcription — not
//! a per-cell cost. So there is no share to compute and no timed comparison that
//! could produce one.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record`, every panic begins `VOID: `,
//! and every one is also a column so a reader can check it without rerunning.
//!
//! - `control_table_mismatches` — this harness's own alternating-sign test
//!   against `AMBIGUOUS_FACES` on all 1,536 pairs. Must be 0, or the population
//!   being derived is not the table's.
//! - `control_ambiguous_pairs` / `control_ambiguous_cases` — must be 192 and
//!   120. 192 is `6 faces × 2 alternating patterns × 2^4 free corners`; 120 is
//!   `256 − 136` and `M-41` records the 136.
//! - `control_decider_transcription_mismatches` — the decider predicate as the
//!   SMT text states it, evaluated in Rust, against the shipped
//!   `ambiguity::face_is_joined`, over 12,288 assignments. Must be 0, or the
//!   solver is being asked about this file's transcription instead of the
//!   crate's rule, which is the failure `M-219` is.
//! - `control_forced_oracle_mismatches` — the exact interval evaluator against
//!   an **independent** oracle on all 1,344 unambiguous pairs: `truth` holds iff
//!   a whole ring edge at constant `q` is inside, which witnesses every `p`. Two
//!   derivations of one predicate, cross-checked; must be 0.
//! - `control_varying_ambiguous` — must equal 192. `truth` has to actually vary
//!   with the magnitudes on every ambiguous face, or the constancy question is
//!   asked of a constant and the `sat` answers are a zero that could not have
//!   been non-zero (`M-44`).
//! - `control_tangency_samples` — must be non-zero. The `all_unit` fixture puts
//!   `d_in == d_out` exactly on **every** ambiguous face, so the sample reaches
//!   the tie the crate documents at ambiguity.rs:45-48 by construction rather
//!   than by seed luck.
//! - `control_numeric_disagreements` — the decider against the interval
//!   evaluator over the same 12,288 assignments. Must be 0: a finite
//!   counterexample would make the solver's `unsat` impossible, and an
//!   instrument's first job is to agree with the other one where they overlap
//!   (`M-279`).
//! - `control_unambiguous_refuted` — **the registration's named control.** Must
//!   be 0: the pipeline re-derives all 1,344 known-correct sub-cases the table
//!   already contains, and a single `sat` there means the encoder disagrees with
//!   the table about a face the table says is not even a decision. Without this
//!   a clean run proves nothing about the solver's fidelity.
//! - `both answers seen` — the run must produce at least one `unsat` and at
//!   least one `sat`. A decision procedure that only ever returns one answer has
//!   not decided anything.

#![allow(
    clippy::cast_precision_loss,
    clippy::too_many_lines,
    reason = "sample counts are cast to f64 for reporting; the sweep is one long function"
)]

mod common;

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

use isomesh::marching_cubes::ambiguity::face_is_joined;
use isomesh::marching_cubes::table::{AMBIGUOUS_FACES, corner_inside, face_bit, face_corners};

/// Per-query cap, handed to the solver as `-T:<seconds>`.
///
/// Two seconds is **about 105× the terminating arms' measured cost**, so it
/// cannot be mistaken for a budget that merely happened to be too small.
/// Measured on this host over a 154-entry subset: the affine encoding decides an
/// entry in **6.7 ms** and the nested encoding *refutes* one in **19 ms**, while
/// the nested encoding's model search does not finish at any cap tried — 2 s
/// here, and 60 s and `qsat` at 20 s while this file was being written. The cap
/// is also what makes the sweep's worst case finite.
const QUERY_TIMEOUT_S: u64 = 2;

/// Whole-run solver budget. Exceeding it panics; it never truncates the sweep.
///
/// Fourteen minutes against a measured projection of about **6.7 minutes**:
/// 1,536 affine derivations at 6.7 ms, 192 + 192 rule queries at 6.7 ms, 192
/// nested refutations at 19 ms, and 192 nested model searches that each spend
/// the whole [`QUERY_TIMEOUT_S`] — 384 s, which is 96% of the run and is the
/// price of recording C1's non-termination branch over the full table rather
/// than over a single point.
const TOTAL_BUDGET: Duration = Duration::from_secs(840);

/// Corner parameter names, in the order `cube`'s corner indices run.
const SYMBOL: [&str; 8] = ["f0", "f1", "f2", "f3", "f4", "f5", "f6", "f7"];

/// Face names by `(axis, side)`, matching `face_bit`'s `axis * 2 + side`.
const FACE_NAME: [[&str; 2]; 3] = [["xlo", "xhi"], ["ylo", "yhi"], ["zlo", "zhi"]];

/// Assignments sampled per `(case, face)` pair: three fixed fixtures and 61
/// drawn ones.
const SAMPLES_PER_FACE: usize = 64;

/// Seed for the drawn assignments. Stated here because a control computed from
/// an unstated seed is not reproducible.
const SAMPLE_SEED: u64 = 0x0141_C0DE_5EED_1141;

/// Which part of the algorithm one of Custodio et al.'s corrections changes.
///
/// C3 asks whether a disagreement this harness finds is *already known*. That
/// question is a scope question, so the scope is the datum rather than a
/// hard-coded overlap count — if a `FaceRule`-scoped correction were ever added
/// to the list below, `custodio_covered` would start reporting it without
/// another line changing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Scope {
    /// Chernyaev's sweeping-plane interior test.
    InteriorTest,
    /// A reference implementation's missing subcases.
    Implementation,
    /// The corner label for a sample lying exactly on the isovalue.
    EqualCorner,
    /// The ambiguous-**face** rule — the asymptotic decider. Nothing Custodio et
    /// al. published is in this scope, which is the whole of C3's answer.
    FaceRule,
}

/// One published correction, as **this repository** recorded it.
///
/// The papers themselves are not read here. `docs/research/2026-08-10-meshing-\
/// library-target.md:82` marks Custodio 2013 acquired and `:85` Custodio,
/// Pesco & Silva 2019, and the four rows below are the repository's own verified
/// readings of them — `V-24`, `V-25`, `V-26` and `M-352`/`P-53` — not a
/// paraphrase of an abstract.
#[derive(Clone, Copy, Debug)]
struct Correction {
    /// The finding id this repository recorded it under.
    finding: &'static str,
    /// The DOI it came from.
    doi: &'static str,
    /// The stage it changes.
    scope: Scope,
}

/// Custodio et al.'s published corrections to Chernyaev, and their scopes.
///
/// Sources, all in-tree: `V-24` (FINDINGS.md:1256) re-derives §5.1's swept
/// saddle `F(t)/Δ(t)` and is quoted in marching_cubes/interior.rs:32-44;
/// `M-165` (FINDINGS.md:1259) measures its reach at 1,966 of 15,625; `V-25`
/// (:1257) records §5.4's missing cases 10 and 12; `V-26` (:1258) the 2019
/// extended triangulation; `M-352` and `P-53`'s registration
/// (experiment.rs:1800-1804) the `=`-corner label, whose §6.2 non-manifold-edge
/// mechanism is quoted at FINDINGS.md:915.
const CORRECTIONS: [Correction; 4] = [
    Correction {
        finding: "V-24",
        doi: "10.1016/j.cag.2013.04.004",
        scope: Scope::InteriorTest,
    },
    Correction {
        finding: "V-25",
        doi: "10.1016/j.cag.2013.04.004",
        scope: Scope::Implementation,
    },
    Correction {
        finding: "M-165",
        doi: "10.1016/j.cag.2013.04.004",
        scope: Scope::InteriorTest,
    },
    Correction {
        finding: "M-352",
        doi: "10.1186/s13173-019-0086-6",
        scope: Scope::EqualCorner,
    },
];

/// The findings this repository recorded the corrections under, as a `|` list.
///
/// Written into every row as `custodio_corrections_known`, so the file names the
/// source of its own overlap judgement instead of leaving it to prose.
fn correction_findings() -> String {
    CORRECTIONS
        .iter()
        .map(|c| c.finding)
        .collect::<Vec<_>>()
        .join("|")
}

/// How many of `disagreements` fall inside a correction Custodio et al. published.
///
/// Every disagreement this harness can find is about the ambiguous-**face**
/// rule, and no correction in [`CORRECTIONS`] is scoped there, so the answer is
/// zero — for a stated reason rather than by assumption.
fn custodio_covered(disagreements: u32) -> u32 {
    if CORRECTIONS.iter().any(|c| c.scope == Scope::FaceRule) {
        disagreements
    } else {
        0
    }
}

/// How `truth` is written in the emitted query.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Encoding {
    /// `forall p . exists q . B(p,q) < 0`, verbatim.
    Nested,
    /// The same, with the inner existential discharged by `B`'s affineness in
    /// `q`. See the module header.
    Affine,
}

impl Encoding {
    /// Column spelling.
    const fn name(self) -> &'static str {
        match self {
            Self::Nested => "nested",
            Self::Affine => "affine",
        }
    }
}

/// What `truth` is compared against.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Claim {
    /// `ambiguity::face_is_joined`: the inside diagonal's product exceeds the
    /// outside diagonal's.
    Decider,
    /// Marching Cubes proper — an ambiguous face's inside corners are always
    /// separated, so the claim is the constant `false`.
    Separate,
    /// The constant `truth` takes at the case's probe assignment. `unsat` then
    /// says `truth` is independent of the magnitudes, which is what
    /// `AMBIGUOUS_FACES` asserts by leaving the bit clear.
    ProbeConstant,
}

impl Claim {
    /// Column spelling.
    const fn name(self) -> &'static str {
        match self {
            Self::Decider => "asymptotic_decider",
            Self::Separate => "marching_cubes_separate",
            Self::ProbeConstant => "truth_is_magnitude_free",
        }
    }
}

/// The solver's answer to one query.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Answer {
    /// The claim holds for every real assignment in the sign class.
    Unsat,
    /// A counterexample exists.
    Sat,
    /// The solver gave up without deciding.
    Unknown,
    /// The solver hit its own wall-clock cap.
    Timeout,
}

impl Answer {
    /// The four words a solver prints, and nothing else.
    fn parse(line: &str) -> Option<Self> {
        match line {
            "unsat" => Some(Self::Unsat),
            "sat" => Some(Self::Sat),
            "unknown" => Some(Self::Unknown),
            "timeout" => Some(Self::Timeout),
            _ => None,
        }
    }

    /// Did the query decide? C1 is exactly this predicate.
    const fn terminated(self) -> bool {
        matches!(self, Self::Unsat | Self::Sat)
    }

    /// Column spelling.
    const fn name(self) -> &'static str {
        match self {
            Self::Unsat => "unsat",
            Self::Sat => "sat",
            Self::Unknown => "unknown",
            Self::Timeout => "timeout",
        }
    }
}

/// One arm of the sweep.
#[derive(Clone, Copy, Debug)]
struct Arm {
    /// Column spelling, and the certificate filename's stem.
    name: &'static str,
    /// How `truth` is written.
    encoding: Encoding,
    /// What it is compared against.
    claim: Claim,
    /// Does the arm issue a query for every `(case, face)` pair, or only for the
    /// ambiguous ones?
    every_face: bool,
}

impl Arm {
    /// The answer that would mean **the shipped table's claim survives** here.
    ///
    /// A different *terminated* answer is a disagreement between the machine
    /// derivation and the committed table, which is what C2 counts. Note this is
    /// not a prediction: for `Claim::Separate` it is the claim being refuted, and
    /// the arm exists because the crate's own documentation says it fails.
    const fn table_claims(self, face: &FaceRef) -> Answer {
        match self.claim {
            // `AMBIGUOUS_FACES` says "there is something to decide here"; a
            // decision exists exactly when `truth` is not constant on the class.
            Claim::ProbeConstant => {
                if face.ambiguous {
                    Answer::Sat
                } else {
                    Answer::Unsat
                }
            }
            // Both rules claim to be a pairing the interpolant agrees with, so
            // holding means no counterexample.
            Claim::Decider | Claim::Separate => Answer::Unsat,
        }
    }

    /// Does a refutation of this arm's claim indict the crate?
    ///
    /// `AsymptoticDecider` and `AMBIGUOUS_FACES` both claim to follow the
    /// interpolant, so a counterexample to either is a defect. `Separate` claims
    /// the opposite in the crate's own words — *"it agrees with the field's own
    /// bilinear interpolant on the face only by luck"*, ambiguity.rs:13-14 — so a
    /// counterexample to it is the documented design, already quantified by
    /// `M-40` (27 of 5,240 cells on `gyroid`) and `M-41` (88 of 256 cases change
    /// their Euler characteristic). This is C3's adjudication, mechanised.
    const fn refutation_is_a_defect(self) -> bool {
        !matches!(self.claim, Claim::Separate)
    }
}

/// The sweep, roots first: the derivation, then the two rules under both
/// encodings.
const ARMS: [Arm; 5] = [
    Arm {
        name: "ambiguity_affine",
        encoding: Encoding::Affine,
        claim: Claim::ProbeConstant,
        every_face: true,
    },
    Arm {
        name: "decider_affine",
        encoding: Encoding::Affine,
        claim: Claim::Decider,
        every_face: false,
    },
    Arm {
        name: "separate_affine",
        encoding: Encoding::Affine,
        claim: Claim::Separate,
        every_face: false,
    },
    Arm {
        name: "decider_nested",
        encoding: Encoding::Nested,
        claim: Claim::Decider,
        every_face: false,
    },
    Arm {
        name: "separate_nested",
        encoding: Encoding::Nested,
        claim: Claim::Separate,
        every_face: false,
    },
];

/// One `(case, face)` pair, resolved once so no loop re-derives it.
#[derive(Clone, Copy, Debug)]
struct FaceRef {
    /// The eight-corner sign configuration.
    case: u8,
    /// The face's normal axis.
    axis: usize,
    /// Near (`0`) or far (`1`) along that axis.
    side: u8,
    /// The four cube-corner indices, counter-clockwise as seen from outside.
    ring: [u8; 4],
    /// Whether each ring corner is inside for this case.
    inside: [bool; 4],
    /// Do the ring signs alternate? Cross-checked against `AMBIGUOUS_FACES`.
    ambiguous: bool,
}

impl FaceRef {
    /// Resolve `(case, axis, side)` through the crate's own primitives.
    fn of(case: u8, axis: usize, side: u8) -> Self {
        let ring = face_corners(axis, side);
        let inside = [
            corner_inside(case, ring[0]),
            corner_inside(case, ring[1]),
            corner_inside(case, ring[2]),
            corner_inside(case, ring[3]),
        ];
        // The same test `build_ambiguous_faces` runs (table.rs:214-220):
        // opposite corners agree, adjacent ones differ.
        let ambiguous = inside[0] == inside[2] && inside[1] == inside[3] && inside[0] != inside[1];
        Self {
            case,
            axis,
            side,
            ring,
            inside,
            ambiguous,
        }
    }

    /// Column spelling of the face.
    const fn name(&self) -> &'static str {
        FACE_NAME[self.axis][self.side as usize]
    }

    /// Is the bit set in the committed table?
    fn marked_by_table(&self) -> bool {
        AMBIGUOUS_FACES[self.case as usize] & face_bit(self.axis, self.side) != 0
    }

    /// The four ring values of an eight-corner assignment, in ring order — the
    /// order `face_is_joined` documents (ambiguity.rs:90-93).
    fn ring_values(&self, corner: &[f64; 8]) -> [f64; 4] {
        [
            corner[self.ring[0] as usize],
            corner[self.ring[1] as usize],
            corner[self.ring[2] as usize],
            corner[self.ring[3] as usize],
        ]
    }
}

/// `{p in [0,1] : (1-p)·a + p·b >= 0}`, as a closed interval, or `None` when
/// empty.
///
/// Closed on purpose. The constraint is `>= 0`, and a single-point intersection
/// is the tangency `S = 0` — the case ambiguity.rs:45-48 resolves to
/// *separated*. Reading it as open would silently adopt the other convention.
fn nonneg_interval(a: f64, b: f64) -> Option<(f64, f64)> {
    match (a >= 0.0, b >= 0.0) {
        (true, true) => Some((0.0, 1.0)),
        (false, false) => None,
        // One root, and it lies in [0,1] because |a| <= |b - a| when the signs
        // differ. `a / (a - b)` rather than `-a / (b - a)`: one subtraction.
        (false, true) => Some((a / (a - b), 1.0)),
        (true, false) => Some((0.0, a / (a - b))),
    }
}

/// `truth`, evaluated exactly for one ring assignment.
///
/// `truth` fails exactly when some `p` has `B(p,0) >= 0` **and** `B(p,1) >= 0`,
/// and both are affine in `p`, so each closed half-plane meets `[0,1]` in a
/// closed sub-interval and the two either intersect or they do not. No sampling,
/// no tolerance, and no blindness at the tangency.
fn truth_at(v: [f64; 4]) -> bool {
    let (Some(lower), Some(upper)) = (nonneg_interval(v[0], v[1]), nonneg_interval(v[3], v[2]))
    else {
        return true;
    };
    let lo = lower.0.max(upper.0);
    let hi = lower.1.min(upper.1);
    lo > hi
}

/// Does the intersection of the two `>= 0` regions collapse to a single point?
///
/// That is `S = 0` exactly. Compared through `total_cmp` rather than `==` so the
/// house rule against bare float equality holds without weakening the test.
fn is_tangency(v: [f64; 4]) -> bool {
    let (Some(lower), Some(upper)) = (nonneg_interval(v[0], v[1]), nonneg_interval(v[3], v[2]))
    else {
        return false;
    };
    let lo = lower.0.max(upper.0);
    let hi = lower.1.min(upper.1);
    lo.total_cmp(&hi) == std::cmp::Ordering::Equal
}

/// The connectivity answer a **non-ambiguous** face forces, from the ring signs
/// alone.
///
/// `q = 0` and `q = 1` are the ring edges `(v0,v1)` and `(v3,v2)`. If either has
/// both ends inside, that constant `q` witnesses every `p` and `truth` holds.
/// Over the twelve non-ambiguous sign patterns the converse holds too: when
/// neither does, the opposite `p`-edge is entirely outside and no `q` works
/// there. Verified against [`truth_at`] on all 1,344 unambiguous pairs — that is
/// `control_forced_oracle_mismatches`.
///
/// **Only valid off the ambiguous patterns.** On an alternating ring the answer
/// depends on the magnitudes, which is the whole subject of this row.
const fn forced_by_ring(inside: [bool; 4]) -> bool {
    (inside[0] && inside[1]) || (inside[3] && inside[2])
}

/// The decider predicate exactly as the emitted SMT text states it.
///
/// Kept beside the shipped `face_is_joined` and checked against it
/// (`control_decider_transcription_mismatches`), because a solver run against a
/// mistyped rule proves something about the typo.
fn emitted_decider(v: [f64; 4], inside_first: bool) -> bool {
    let d02 = v[0] * v[2];
    let d13 = v[1] * v[3];
    if inside_first { d02 > d13 } else { d13 > d02 }
}

/// The probe assignment for a case: corner `i` at `-(1 + i/8)` inside, `1 + i/8`
/// outside.
///
/// Any member of the sign class serves [`Claim::ProbeConstant`], whose question
/// is whether `truth` depends on the magnitudes at all. The offsets keep the
/// eight values distinct so the two diagonal products cannot tie, which would
/// put the probe on the boundary of its own class instead of inside it.
fn probe(case: u8) -> [f64; 8] {
    let mut out = [0.0f64; 8];
    for (i, slot) in out.iter_mut().enumerate() {
        let magnitude = 1.0 + i as f64 * 0.125;
        *slot = if corner_inside(case, i as u8) {
            -magnitude
        } else {
            magnitude
        };
    }
    out
}

/// The three fixed fixtures every sample set opens with, in order.
///
/// They make the vacuity controls hold **by construction** rather than by seed
/// luck, and each is a member of every case's sign class:
///
/// 1. `all_unit` — inside `-1`, outside `+1`. On any ambiguous face
///    `d_in = d_out = 1`, so this is an exact tangency on all 192 of them.
/// 2. `joined` — inside `-1`, outside `0`. `d_in = 1 > 0 = d_out`: `truth` holds.
///    Zero is in the class because `cube::is_inside` puts an exact zero outside.
/// 3. `separated` — inside `-0.25`, outside `1.75`. `d_in = 0.0625` against
///    `d_out = 3.0625`: `truth` fails.
///
/// So `truth` provably takes both values on every ambiguous face, and the tie is
/// provably reached.
fn fixtures(case: u8) -> [[f64; 8]; 3] {
    let mut out = [[0.0f64; 8]; 3];
    for (inside_value, outside_value, row) in
        [(-1.0, 1.0, 0usize), (-1.0, 0.0, 1), (-0.25, 1.75, 2)]
    {
        for i in 0u8..8 {
            out[row][i as usize] = if corner_inside(case, i) {
                inside_value
            } else {
                outside_value
            };
        }
    }
    out
}

/// A drawn assignment inside `case`'s sign class.
///
/// Inside magnitudes come from `{0.25, .., 2.0}` and outside from
/// `{0, 0.25, .., 1.75}`. Zero is in the outside set deliberately: the class
/// contains it (cube.rs:171-173) and a sample that never reaches it never tests
/// a vanishing diagonal product.
fn drawn(case: u8, rng: &mut common::poly::Rng) -> [f64; 8] {
    let mut out = [0.0f64; 8];
    for i in 0u8..8 {
        let step = rng.next_i64_in(0, 8) as f64 * 0.25;
        out[i as usize] = if corner_inside(case, i) {
            -(step + 0.25)
        } else {
            step
        };
    }
    out
}

/// The `q = <side>` ring edge of the bilinear, as an affine expression in `p`.
fn edge_expr(a: &str, b: &str) -> String {
    format!("(+ (* {a} (- 1.0 p)) (* {b} p))")
}

/// The bilinear itself, as an expression in `p` and `q`.
fn bilinear_expr(v: [&str; 4]) -> String {
    format!(
        "(+ (* {} (- 1.0 p) (- 1.0 q)) (* {} p (- 1.0 q)) (* {} p q) (* {} (- 1.0 p) q))",
        v[0], v[1], v[2], v[3]
    )
}

/// The SMT-LIB2 text for one query.
///
/// Self-describing: the `;` preamble names the case, the face, the ring and the
/// claim, so the emitted file is a certificate a reader can regenerate and
/// re-check without this harness. That is C3's answer.
fn emit(face: &FaceRef, arm: Arm, probe_truth: bool) -> String {
    let v = [
        SYMBOL[face.ring[0] as usize],
        SYMBOL[face.ring[1] as usize],
        SYMBOL[face.ring[2] as usize],
        SYMBOL[face.ring[3] as usize],
    ];
    let mut s = String::with_capacity(1024);
    let _ = writeln!(
        s,
        "; isomesh P-141 / R-141 -- one entry of the case table, derived."
    );
    let _ = writeln!(s, "; arm       {}", arm.name);
    let _ = writeln!(s, "; case      {}", face.case);
    let _ = writeln!(
        s,
        "; face      {} (axis {} side {}, bit {})",
        face.name(),
        face.axis,
        face.side,
        face_bit(face.axis, face.side)
    );
    let _ = writeln!(
        s,
        "; ring      {} {} {} {}  at (0,0) (1,0) (1,1) (0,1), ccw from outside",
        v[0], v[1], v[2], v[3]
    );
    let _ = writeln!(s, "; ambiguous {}", face.ambiguous);
    let _ = writeln!(s, "; encoding  {}", arm.encoding.name());
    let _ = writeln!(s, "; claim     {}", arm.claim.name());
    let _ = writeln!(s, ";");
    let _ = writeln!(
        s,
        "; unsat = the claim holds for EVERY real corner assignment in this sign class."
    );
    let _ = writeln!(s, "; sat   = a counterexample exists.");
    for name in SYMBOL {
        let _ = writeln!(s, "(declare-const {name} Real)");
    }
    for i in 0u8..8 {
        let name = SYMBOL[i as usize];
        if corner_inside(face.case, i) {
            let _ = writeln!(s, "(assert (< {name} 0.0))");
        } else {
            let _ = writeln!(s, "(assert (>= {name} 0.0))");
        }
    }
    let _ = writeln!(s, "(define-fun truth () Bool");
    let _ = writeln!(s, "  (forall ((p Real))");
    let _ = writeln!(s, "    (=> (and (<= 0.0 p) (<= p 1.0))");
    match arm.encoding {
        Encoding::Nested => {
            let _ = writeln!(s, "      (exists ((q Real))");
            let _ = writeln!(s, "        (and (<= 0.0 q) (<= q 1.0)");
            let _ = writeln!(s, "          (< {} 0.0))))))", bilinear_expr(v));
        }
        Encoding::Affine => {
            let _ = writeln!(s, "      (or (< {} 0.0)", edge_expr(v[0], v[1]));
            let _ = writeln!(s, "          (< {} 0.0)))))", edge_expr(v[3], v[2]));
        }
    }
    let claim = match arm.claim {
        Claim::Decider => {
            // `face_is_joined`: whichever diagonal is inside supplies `d_in`
            // (ambiguity.rs:108-116). The signs are fixed by the case, so the
            // branch is resolved here and the emitted text carries no condition.
            let (d_in, d_out) = if face.inside[0] {
                ((v[0], v[2]), (v[1], v[3]))
            } else {
                ((v[1], v[3]), (v[0], v[2]))
            };
            format!("(> (* {} {}) (* {} {}))", d_in.0, d_in.1, d_out.0, d_out.1)
        }
        Claim::Separate => String::from("false"),
        Claim::ProbeConstant => probe_truth.to_string(),
    };
    let _ = writeln!(s, "(assert (not (= {claim} truth)))");
    let _ = writeln!(s, "(check-sat)");
    s
}

/// FNV-1a over the query bytes — the same hash family `validate::mesh_hash` uses.
fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// The solver binary: `Z3` if set, `z3` otherwise.
fn solver_path() -> PathBuf {
    std::env::var_os("Z3").map_or_else(|| PathBuf::from("z3"), PathBuf::from)
}

/// `(column value, raw version string)`, or a panic naming the binary.
///
/// # Panics
///
/// If the solver cannot be executed or its `--version` fails. Deliberate, and
/// the one place this harness refuses to degrade: R-141's whole content is what
/// an external decision procedure says, and a row written without one would be a
/// row about nothing.
fn solver_id(solver: &Path) -> (String, String) {
    let raw = match Command::new(solver).arg("--version").output() {
        Ok(out) if out.status.success() => String::from_utf8_lossy(&out.stdout).trim().to_string(),
        Ok(out) => panic!(
            "P-141: `{} --version` exited with {} instead of reporting a version. R-141 needs a \
             real-quantifier-elimination solver and no CSV row is written for a solver that did \
             not run; point the `Z3` environment variable at a working z3 executable.",
            solver.display(),
            out.status
        ),
        Err(err) => panic!(
            "P-141: cannot execute the real-quantifier-elimination solver at `{}` ({err}). R-141 \
             is the one row in Phase 27 that needs an external binary and it is not optional: set \
             the `Z3` environment variable to a z3 executable, or put `z3` on PATH. No CSV row is \
             written for a solver that did not run.",
            solver.display()
        ),
    };
    let name = raw
        .split_whitespace()
        .next()
        .unwrap_or("solver")
        .to_ascii_lowercase();
    let version = raw
        .split_whitespace()
        .find(|token| token.starts_with(|c: char| c.is_ascii_digit()) && token.contains('.'))
        .unwrap_or("unknown");
    (format!("{name}-{version}"), raw.replace(' ', "_"))
}

/// Ask the solver one question.
///
/// The cap is the solver's own `-T:`, so the process cannot outlive it. z3 exits
/// `0` for `sat`, `unsat` and `timeout` alike, so the answer is read from stdout
/// and a malformed query is a panic rather than a data point.
///
/// # Panics
///
/// If the solver cannot be spawned, rejects the query, or prints no verdict.
fn ask(solver: &Path, file: &Path) -> Answer {
    let out = Command::new(solver)
        .arg("-smt2")
        .arg(format!("-T:{QUERY_TIMEOUT_S}"))
        .arg(file)
        .output()
        .unwrap_or_else(|err| {
            panic!(
                "P-141: the solver at `{}` could not be spawned for {} ({err})",
                solver.display(),
                file.display()
            )
        });
    let text = String::from_utf8_lossy(&out.stdout);
    let mut answer = None;
    for raw in text.lines() {
        let line = raw.trim();
        assert!(
            !line.starts_with("(error"),
            "P-141: the solver rejected {}: {line}\nThe emitted query is the certificate, so a \
             parse error here is a bug in this harness's emitter, not a measurement.",
            file.display()
        );
        if let Some(parsed) = Answer::parse(line) {
            answer = Some(parsed);
        }
    }
    answer.unwrap_or_else(|| {
        panic!(
            "P-141: the solver printed no verdict for {}; stdout was {:?}",
            file.display(),
            text.trim()
        )
    })
}

/// Everything one issued query produced.
struct Outcome {
    /// Index into [`ARMS`].
    arm: usize,
    /// The entry it decided.
    face: FaceRef,
    /// Which constant `Claim::ProbeConstant` asserted.
    probe_truth: bool,
    /// What came back.
    answer: Answer,
    /// Wall clock for the whole spawn-and-decide, in seconds.
    wall: f64,
    /// Certificate size.
    bytes: usize,
    /// Certificate hash.
    hash: u64,
    /// Certificate filename, inside the printed directory.
    file: String,
}

/// One arm's totals over the entries it issued.
#[derive(Clone, Copy, Debug, Default)]
struct Totals {
    /// Queries issued.
    queries: u32,
    /// Queries that decided. C1's aggregate.
    derived: u32,
    /// Decided queries whose answer is the one the committed table needs.
    matching: u32,
    /// Decided queries that contradict the committed table. C2's aggregate.
    disagreeing: u32,
    /// Solver wall clock, in seconds.
    solver_s: f64,
}

/// The Rust-side cross-checks, all computed before a single query is issued.
#[derive(Clone, Copy, Debug, Default)]
struct Controls {
    /// `(case, face)` pairs walked.
    pairs: u32,
    /// Pairs this harness calls ambiguous. Must be 192.
    ambiguous_pairs: u32,
    /// Cases with at least one ambiguous face. Must be 120.
    ambiguous_cases: u32,
    /// Pairs where this harness and `AMBIGUOUS_FACES` disagree. Must be 0.
    table_mismatches: u32,
    /// Assignments where the emitted decider text and `face_is_joined` differ.
    decider_transcription_mismatches: u32,
    /// Assignments where `truth_at` and the independent ring oracle differ.
    forced_oracle_mismatches: u32,
    /// Ambiguous faces where `truth` takes both values over the sample.
    varying_ambiguous: u32,
    /// Assignments sitting exactly on the tangency `S = 0`.
    tangency_samples: u32,
    /// Assignments where the decider and `truth_at` differ. Must be 0.
    numeric_disagreements: u32,
    /// Assignments examined.
    samples: u32,
}

/// Walk every `(case, face)` pair and cross-check the three evaluators.
fn measure_controls(faces: &[FaceRef]) -> Controls {
    let mut c = Controls::default();
    let mut rng = common::poly::Rng::new(SAMPLE_SEED);
    let mut case_has_ambiguous = [false; 256];

    for face in faces {
        c.pairs += 1;
        if face.ambiguous != face.marked_by_table() {
            c.table_mismatches += 1;
        }
        if face.ambiguous {
            c.ambiguous_pairs += 1;
            case_has_ambiguous[face.case as usize] = true;
        }

        let mut saw_true = false;
        let mut saw_false = false;
        let fixed = fixtures(face.case);
        for index in 0..SAMPLES_PER_FACE {
            let corner = if index < fixed.len() {
                fixed[index]
            } else {
                drawn(face.case, &mut rng)
            };
            let v = face.ring_values(&corner);
            let truth = truth_at(v);
            c.samples += 1;
            saw_true |= truth;
            saw_false |= !truth;
            if is_tangency(v) {
                c.tangency_samples += 1;
            }
            if face.ambiguous {
                let emitted = emitted_decider(v, face.inside[0]);
                if emitted != face_is_joined(v) {
                    c.decider_transcription_mismatches += 1;
                }
                if emitted != truth {
                    c.numeric_disagreements += 1;
                }
            } else if forced_by_ring(face.inside) != truth {
                c.forced_oracle_mismatches += 1;
            }
        }
        if face.ambiguous && saw_true && saw_false {
            c.varying_ambiguous += 1;
        }
    }

    c.ambiguous_cases = case_has_ambiguous.iter().filter(|has| **has).count() as u32;
    c
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-141");

    common::experiment::run(prereg, |run| {
        // ── the solver, before anything else ─────────────────────────────────
        let solver = solver_path();
        let (solver_column, solver_raw) = solver_id(&solver);
        println!("solver: {solver_raw} at {}", solver.display());
        println!("per-query cap: {QUERY_TIMEOUT_S}s (the solver's own -T:)\n");
        println!(
            "Custodio et al.'s corrections, as this repository verified them -- C3's reference set:"
        );
        for correction in CORRECTIONS {
            println!(
                "  {:<6} {:<28} scope {:?}",
                correction.finding, correction.doi, correction.scope
            );
        }
        println!();

        // ── the population, resolved through the crate's own primitives ──────
        let mut faces: Vec<FaceRef> = Vec::with_capacity(1536);
        for case in 0u8..=255 {
            for axis in 0..3usize {
                for side in 0u8..2 {
                    faces.push(FaceRef::of(case, axis, side));
                }
            }
        }

        // ── Rust-side controls, before a single query ────────────────────────
        let controls = measure_controls(&faces);
        assert_eq!(
            controls.table_mismatches, 0,
            "VOID: this harness's alternating-sign test disagrees with AMBIGUOUS_FACES on {} of \
             {} pairs, so the population being derived is not the committed table's and no \
             agreement or disagreement below is about the crate",
            controls.table_mismatches, controls.pairs
        );
        assert_eq!(
            controls.ambiguous_pairs, 192,
            "VOID: {} ambiguous (case, face) pairs, not 192 = 6 faces x 2 alternating patterns x \
             2^4 free corners. The arms' denominators are wrong and every rate below is over the \
             wrong population",
            controls.ambiguous_pairs
        );
        assert_eq!(
            controls.ambiguous_cases, 120,
            "VOID: {} cases carry an ambiguous face, not 120 = 256 - 136, and M-41 records the \
             136 cases with none",
            controls.ambiguous_cases
        );
        assert_eq!(
            controls.decider_transcription_mismatches, 0,
            "VOID: the decider predicate as the emitted SMT text states it disagrees with the \
             shipped ambiguity::face_is_joined on {} of {} assignments, so the solver would be \
             deciding a question about this file's transcription rather than about the crate's \
             rule -- which is exactly the defect M-219 is",
            controls.decider_transcription_mismatches, controls.samples
        );
        assert_eq!(
            controls.forced_oracle_mismatches, 0,
            "VOID: the exact interval evaluator disagrees with the independent ring oracle on {} \
             assignments over the 1,344 unambiguous pairs. The two are separate derivations of one \
             predicate and a clean solver run against an unchecked evaluator proves nothing",
            controls.forced_oracle_mismatches
        );
        assert_eq!(
            controls.varying_ambiguous, controls.ambiguous_pairs,
            "VOID: truth varies with the magnitudes on only {} of {} ambiguous faces, so on the \
             rest the constancy question is asked of a constant and its answer is a zero that \
             could not have been non-zero (M-44). The `joined` and `separated` fixtures are \
             supposed to make this hold by construction",
            controls.varying_ambiguous, controls.ambiguous_pairs
        );
        assert!(
            controls.tangency_samples > 0,
            "VOID: no sampled assignment reaches the tangency S = 0, so this run never tests the \
             tie ambiguity.rs:45-48 resolves to separated -- and the tangency is precisely the \
             measure-zero set a grid sweep cannot see, which is the whole reason this row needs a \
             solver. The `all_unit` fixture is supposed to guarantee one per ambiguous face"
        );
        assert_eq!(
            controls.numeric_disagreements, 0,
            "VOID: the decider and the interval evaluator already disagree on {} of {} finite \
             assignments, so any `unsat` the solver returns below would be impossible and one of \
             the two instruments is wrong (M-279: a new instrument's first job is to agree with \
             the old one where they overlap)",
            controls.numeric_disagreements, controls.samples
        );
        println!(
            "controls: {} pairs ({} ambiguous over {} cases), {} assignments, {} tangencies, \
             0 mismatches on all four cross-checks",
            controls.pairs,
            controls.ambiguous_pairs,
            controls.ambiguous_cases,
            controls.samples,
            controls.tangency_samples
        );

        // ── the sweep ────────────────────────────────────────────────────────
        let dir = std::env::temp_dir().join("isomesh-p141");
        std::fs::create_dir_all(&dir)
            .expect("P-141: a writable temporary directory for the certificates");
        println!("certificates: {}\n", dir.display());

        let started = Instant::now();
        let mut outcomes: Vec<Outcome> = Vec::new();
        let mut totals = [Totals::default(); ARMS.len()];

        for (index, arm) in ARMS.iter().enumerate() {
            let arm_started = Instant::now();
            for face in &faces {
                if !arm.every_face && !face.ambiguous {
                    continue;
                }
                assert!(
                    started.elapsed() < TOTAL_BUDGET,
                    "P-141: the whole-run solver budget of {}s is spent with arm `{}` unfinished. \
                     A truncated sweep is not a derivation, so this aborts rather than writing a \
                     short CSV that looks complete",
                    TOTAL_BUDGET.as_secs(),
                    arm.name
                );

                let probe_truth = truth_at(face.ring_values(&probe(face.case)));
                let text = emit(face, *arm, probe_truth);
                let name = format!("{}__case{:03}__{}.smt2", arm.name, face.case, face.name());
                let path = dir.join(&name);
                std::fs::write(&path, text.as_bytes())
                    .expect("P-141: the emitted query must be writable; it is the certificate");

                let query_started = Instant::now();
                let answer = ask(&solver, &path);
                let wall = query_started.elapsed().as_secs_f64();

                let tally = &mut totals[index];
                tally.queries += 1;
                tally.solver_s += wall;
                if answer.terminated() {
                    tally.derived += 1;
                    if answer == arm.table_claims(face) {
                        tally.matching += 1;
                    } else {
                        tally.disagreeing += 1;
                    }
                }
                outcomes.push(Outcome {
                    arm: index,
                    face: *face,
                    probe_truth,
                    answer,
                    wall,
                    bytes: text.len(),
                    hash: fnv1a64(text.as_bytes()),
                    file: name,
                });
            }
            let t = totals[index];
            println!(
                "{:<18} {:>5} queries  {:>5} decided  {:>5} matching  {:>5} disagreeing  {:>7.2}s",
                arm.name,
                t.queries,
                t.derived,
                t.matching,
                t.disagreeing,
                arm_started.elapsed().as_secs_f64()
            );
        }

        // ── the registration's named control, over the solver's own answers ──
        let unambiguous_refuted = outcomes
            .iter()
            .filter(|o| {
                matches!(ARMS[o.arm].claim, Claim::ProbeConstant)
                    && !o.face.ambiguous
                    && o.answer != Answer::Unsat
            })
            .count() as u32;
        assert_eq!(
            unambiguous_refuted, 0,
            "VOID: the pipeline failed to re-derive {} of the 1,344 known-correct sub-cases the \
             table already contains -- faces AMBIGUOUS_FACES says carry no decision at all. The \
             registration's vacuity control is exactly this: without it a clean run proves \
             nothing about the solver's fidelity",
            unambiguous_refuted
        );
        let total_unsat = outcomes
            .iter()
            .filter(|o| o.answer == Answer::Unsat)
            .count();
        let total_sat = outcomes.iter().filter(|o| o.answer == Answer::Sat).count();
        assert!(
            total_unsat > 0 && total_sat > 0,
            "VOID: the run returned {total_unsat} unsat and {total_sat} sat. A decision procedure \
             that only ever produces one answer has not decided anything, and every column below \
             would be a property of the harness rather than of the table"
        );

        // ── the verdicts ─────────────────────────────────────────────────────
        let total_queries: u32 = totals.iter().map(|t| t.queries).sum();
        let total_derived: u32 = totals.iter().map(|t| t.derived).sum();
        let total_disagreeing: u32 = totals.iter().map(|t| t.disagreeing).sum();
        let total_solver_s: f64 = totals.iter().map(|t| t.solver_s).sum();
        let adjudicated_defects: u32 = ARMS
            .iter()
            .zip(&totals)
            .filter(|(arm, _)| arm.refutation_is_a_defect())
            .map(|(_, t)| t.disagreeing)
            .sum();
        let defects_outside_custodio: u32 = ARMS
            .iter()
            .zip(&totals)
            .filter(|(arm, _)| arm.refutation_is_a_defect())
            .map(|(_, t)| t.disagreeing - custodio_covered(t.disagreeing))
            .sum();

        // C2 is falsified by zero disagreements between the derived taxonomy and
        // the committed table; the `separate_*` arms compare against `CASES` at
        // mask zero, which is the crate's shipped default.
        let c2 = total_disagreeing > 0;
        // C3 is falsified by disagreements outside Custodio's set that SURVIVE
        // adjudication. `Arm::refutation_is_a_defect` is the adjudication.
        let c3 = defects_outside_custodio == 0;
        println!(
            "\n{total_derived} of {total_queries} entries decided in {total_solver_s:.2}s; \
             {total_disagreeing} disagree with the committed table, {adjudicated_defects} of them \
             indict it, {defects_outside_custodio} outside Custodio's published set"
        );
        println!(
            "C2 {}   C3 {}",
            if c2 { "HELD" } else { "FALSIFIED" },
            if c3 { "HELD" } else { "FALSIFIED" }
        );

        // ── rows ─────────────────────────────────────────────────────────────
        for out in &outcomes {
            let arm = ARMS[out.arm];
            let t = totals[out.arm];
            let terminated = out.answer.terminated();
            let claims = arm.table_claims(&out.face);
            let agrees = terminated && out.answer == claims;
            // A control row is one whose answer is fixed in advance: the
            // 1,344 unambiguous re-derivations, and both `separate` arms.
            let is_control = matches!(arm.claim, Claim::Separate)
                || (matches!(arm.claim, Claim::ProbeConstant) && !out.face.ambiguous);
            run.record(&[
                ("solver", solver_column.clone()),
                (
                    "configuration",
                    format!("{}|case{:03}|{}", arm.name, out.face.case, out.face.name()),
                ),
                ("wall_clock_s", format!("{:.6}", out.wall)),
                ("terminated", terminated.to_string()),
                ("cases_derived", t.derived.to_string()),
                ("cases_matching_table", t.matching.to_string()),
                ("cases_disagreeing", t.disagreeing.to_string()),
                (
                    "custodio_disagreement_overlap",
                    format!("{}_of_{}", custodio_covered(t.disagreeing), t.disagreeing),
                ),
                (
                    "machine_checkable_certificate",
                    format!("{}b_{:016x}", out.bytes, out.hash),
                ),
                // C1 *is* the termination claim, so this column is `terminated`
                // by construction; the number worth reading is the aggregate.
                ("c1_holds", terminated.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // ── extras (M-273) ──
                ("arm", arm.name.to_string()),
                ("encoding", arm.encoding.name().to_string()),
                ("rule", arm.claim.name().to_string()),
                ("is_control", is_control.to_string()),
                ("case_index", out.face.case.to_string()),
                ("face", out.face.name().to_string()),
                ("face_ambiguous", out.face.ambiguous.to_string()),
                ("table_claims", claims.name().to_string()),
                ("query_status", out.answer.name().to_string()),
                ("agrees_with_table", agrees.to_string()),
                (
                    "refutation_is_a_defect",
                    arm.refutation_is_a_defect().to_string(),
                ),
                ("probe_truth", out.probe_truth.to_string()),
                ("query_timeout_s", QUERY_TIMEOUT_S.to_string()),
                ("arm_queries", t.queries.to_string()),
                ("arm_solver_s", format!("{:.3}", t.solver_s)),
                ("smt2_bytes", out.bytes.to_string()),
                ("smt2_hash", format!("{:016x}", out.hash)),
                ("smt2_file", out.file.clone()),
                ("total_solver_s", format!("{total_solver_s:.3}")),
                ("total_disagreeing", total_disagreeing.to_string()),
                ("adjudicated_defects", adjudicated_defects.to_string()),
                (
                    "defects_outside_custodio",
                    defects_outside_custodio.to_string(),
                ),
                ("solver_raw", solver_raw.clone()),
                ("control_samples", controls.samples.to_string()),
                (
                    "control_ambiguous_pairs",
                    controls.ambiguous_pairs.to_string(),
                ),
                (
                    "control_ambiguous_cases",
                    controls.ambiguous_cases.to_string(),
                ),
                (
                    "control_table_mismatches",
                    controls.table_mismatches.to_string(),
                ),
                (
                    "control_decider_transcription_mismatches",
                    controls.decider_transcription_mismatches.to_string(),
                ),
                (
                    "control_forced_oracle_mismatches",
                    controls.forced_oracle_mismatches.to_string(),
                ),
                (
                    "control_varying_ambiguous",
                    controls.varying_ambiguous.to_string(),
                ),
                (
                    "control_tangency_samples",
                    controls.tangency_samples.to_string(),
                ),
                (
                    "control_numeric_disagreements",
                    controls.numeric_disagreements.to_string(),
                ),
                (
                    "control_unambiguous_refuted",
                    unambiguous_refuted.to_string(),
                ),
                ("custodio_corrections_known", correction_findings()),
            ]);
        }
    });
}
