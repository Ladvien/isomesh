//! Pre-registration, as a compile error.
//!
//! Ticket: R-000. *"The feedback loop is currently a discipline; make it a
//! compile error."*
//!
//! # What was wrong with the discipline
//!
//! `FINDINGS.md` Part 4 holds seven pre-registered predictions and the record is
//! good — P-6 was falsified and the entry says so, P-7 held. The weakness is
//! that nothing stopped an experiment from being *written first and predicted
//! afterwards*, which is the failure the whole practice exists to prevent and is
//! invisible in the artefact: a prediction recorded after the numbers came in
//! reads exactly like one recorded before.
//!
//! So the prediction moves here, where the compiler can see it, and
//! [`experiment!`](macro@crate::experiment) refuses to build an experiment whose id is
//! not in [`PREREGISTERED`]. Registering one *is* a commit, so the git history
//! carries the ordering that prose cannot.
//!
//! # This file is the source and `FINDINGS.md` quotes it
//!
//! Two copies of a hypothesis drift, and the drift is undetectable — a
//! pre-registration that quietly acquired a clause is worse than none.
//! `scripts/backlog_gate.sh` checks that every id here appears in `FINDINGS.md`,
//! so the prose can elaborate and cannot contradict.
//!
//! # Registering is cheap; changing a registration is not
//!
//! Amending a [`Preregistration`] after its experiment has run is a **rewrite of
//! the prediction**, and the only honest way to do it is to register a new id
//! and record why the first one was inadequate. Nothing here enforces that —
//! git does.
//!
//! ```
//! use isomesh::experiment;
//!
//! // Fails to compile if "P-8" is not in `PREREGISTERED`.
//! let p = experiment!("P-8");
//! assert!(p.hypothesis.starts_with("A weld gated"));
//! ```

#[cfg(test)]
mod tests;

/// A prediction made before the experiment that tests it.
///
/// `#[non_exhaustive]`, so a consumer cannot build one: the only instances are
/// the constants in [`PREREGISTERED`], which is what makes
/// [`experiment!`](macro@crate::experiment) the sole way to obtain one and therefore
/// the sole gate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct Preregistration {
    /// `P-n`, unique and never reused.
    pub id: &'static str,
    /// The ticket that runs it.
    pub ticket: &'static str,
    /// What is predicted, in one sentence, stated so it can fail.
    pub hypothesis: &'static str,
    /// The observation that would refute it. **A hypothesis with no falsifier is
    /// not registered** — there is no way to express one here without this
    /// field, which is the point of the field existing.
    pub falsified_by: &'static str,
    /// Column names the run must emit. The harness writes exactly these, so a
    /// metric quietly dropped between prediction and run is a missing column
    /// rather than a silence.
    pub records: &'static [&'static str],
}

/// Every prediction registered so far.
///
/// Ordered by id. Ids are never reused, including for a prediction that was
/// withdrawn — a gap in the sequence is information.
pub const PREREGISTERED: &[Preregistration] = &[
    Preregistration {
        id: "P-8",
        ticket: "R-001",
        hypothesis: "A weld gated on Lk u ∩ Lk v = ∅, leaving rejected pairs \
                     split, yields exactly 0 non-manifold edges and 0 \
                     non-manifold vertices on all eight fields × all extractors, \
                     where the unconditional weld yields N > 0.",
        falsified_by: "The gated weld still producing non-manifold output — \
                       which would prove the surface link condition insufficient \
                       for index-buffer realisation, and is the more interesting \
                       result.",
        records: &[
            "non_manifold_edges_ungated",
            "non_manifold_edges_gated",
            "non_manifold_vertices_ungated",
            "non_manifold_vertices_gated",
            "rejected_merges",
            "vertex_delta",
            "weld_ms_ungated",
            "weld_ms_gated",
        ],
    },
    Preregistration {
        id: "P-9",
        ticket: "R-002",
        hypothesis: "For buckets of ≥3 coincident vertices, at least one \
                     reference field yields ≥2 distinct outputs across P seeded \
                     permutations of within-bucket merge order.",
        falsified_by: "All P permutations byte-identical on every field — \
                       meaning the k-way weld is confluent and no canonical \
                       order is needed.",
        records: &[
            "distinct_outputs",
            "vertex_count_spread",
            "buckets_of_three_or_more",
        ],
    },
    Preregistration {
        id: "P-10",
        ticket: "R-003",
        hypothesis: "Vertex inflation from gated-weld-plus-split is < 1%, and \
                     self-intersections per 1k are unchanged from the \
                     unconditional weld.",
        falsified_by: "Inflation > 1%, meaning a real merge/split trade-off \
                       needing a stated policy; or self-intersections rising, \
                       meaning M-93's duplication artefact returns and the \
                       metric must be defined on welded output only.",
        records: &[
            "vertex_inflation_pct",
            "self_intersections_per_1k_gated",
            "self_intersections_per_1k_ungated",
        ],
    },
    Preregistration {
        id: "P-11",
        ticket: "R-004",
        hypothesis: "With one canonical world_of_sample rather than an \
                     offset-and-add, seam cracks fall to 0 for all cell sizes — \
                     not only powers of two — and M-73's hairline disappears with \
                     no change to the transition-cell construction.",
        falsified_by: "Cracks surviving canonical reconstruction, which localises \
                       the defect back in Transvoxel and makes it a different \
                       ticket.",
        records: &[
            "cell_size",
            "lod_pair",
            "crack_count",
            "max_discontinuity",
            "arithmetic",
        ],
    },
    Preregistration {
        id: "P-12",
        ticket: "R-005",
        hypothesis: "The dual's superlinear cost is the four-cells-around-a-\
                     crossed-edge gather at stride n²: cache-miss count per \
                     sample rises with n for Surface Nets and stays flat for \
                     Marching Cubes.",
        falsified_by: "Flat miss rates, pointing at branch misprediction or \
                       allocation instead.",
        records: &[
            "samples",
            "extractor",
            "cache_misses_per_sample",
            "ns_per_sample",
        ],
    },
    Preregistration {
        id: "P-13",
        ticket: "R-006",
        hypothesis: "M-66's non-convergent angle is bounded below by the \
                     dihedral angle of the feature, so it is a property of sharp \
                     edges rather than of resolution and is predictable from the \
                     field.",
        falsified_by: "The angle failing to track the dihedral prediction, which \
                       makes it a defect with a location rather than a property.",
        records: &[
            "dihedral_deg",
            "samples",
            "measured_angle_deg",
            "predicted_angle_deg",
        ],
    },
    Preregistration {
        id: "P-14",
        ticket: "R-003",
        hypothesis: "The residual non-manifold vertices under Surface Nets and \
                     Dual Contouring are the one-vertex-per-cell rule meeting a \
                     cell that contains more than one surface component: almost \
                     all of them sit in cells where Manifold Dual Contouring \
                     emits more than one vertex, and MDC's own count is strictly \
                     lower on every field where either is non-zero.",
        falsified_by: "A substantial share of non-manifold vertices sitting in \
                       cells where MDC emits exactly one vertex — which would \
                       mean the defect is not the one-vertex-per-cell rule and \
                       is somewhere else entirely.",
        records: &[
            "non_manifold_vertices",
            "non_manifold_edges",
            "multi_vertex_cells",
            "nm_vertices_in_single_vertex_cells",
            "worst_link_components",
        ],
    },
    Preregistration {
        id: "P-15",
        ticket: "R-007",
        hypothesis: "More than half of the dual mesher's cycles per sample are \
                     spent in emit_quads, which is three unconditional O(n³) \
                     sweeps over the sample grid rather than work proportional \
                     to the surface.",
        falsified_by: "emit_quads accounting for half or less, which puts the \
                       cost in sample or place_vertices — work Marching Cubes \
                       does too, at four times the IPC — and means the dual's \
                       IPC is lost to something other than its extra traversal.",
        records: &[
            "stage",
            "cycles_per_sample",
            "instructions_per_sample",
            "ipc",
            "samples",
        ],
    },
    Preregistration {
        id: "P-16",
        ticket: "R-008",
        hypothesis: "Every vertex whose area-weighted normal is more than 90° \
                     from the field gradient lies on a grid edge at least one of \
                     whose incident cells straddles the crease — its eight \
                     corners do not all have the same nearer plane — so the \
                     phenomenon is two faces meeting inside one cell rather than \
                     a winding or ordering defect.",
        falsified_by: "More than 5% of past-90° vertices whose incident cells \
                       all lie on one side of the crease, which would make it a \
                       defect with a location and a fix rather than a property \
                       of one vertex per crossed edge.",
        records: &[
            "dihedral_deg",
            "samples",
            "past90_vertices",
            "past90_in_straddling_cell",
            "offending_faces_per_past90_vertex",
        ],
    },
    Preregistration {
        id: "P-17",
        ticket: "A-025",
        hypothesis: "Manifold Dual Contouring's residual non-manifold edges are \
                     an interior ambiguity the face decider cannot see: the \
                     cells sharing an offending ambiguous face report \
                     Interior::Joined on at least one axis sweep at a rate far \
                     above the same measurement over non-offending \
                     ambiguous-face pairs.",
        falsified_by: "Offending pairs reporting Interior::Joined at about the \
                       same rate as the control — which would mean the residue \
                       is not the interior ambiguity and the one-vertex-per-cycle \
                       rule fails for a reason nothing in the A-002 series has \
                       named.",
        records: &[
            "samples",
            "face_rule",
            "offending_pairs",
            "offending_with_interior_join",
            "control_with_interior_join",
        ],
    },
    Preregistration {
        id: "P-18",
        ticket: "R-011",
        hypothesis: "Every precondition a published convex decomposition method \
                     requires of its input mesh is already reported by \
                     ColliderReadiness, so a caller holding a readiness report \
                     has everything it needs to decide whether a mesh can be \
                     handed to a decomposer.",
        falsified_by: "Any precondition required by a method in the audit that \
                       no ColliderReadiness field reports. \
                       Self-intersection-freedom is the standing candidate and \
                       this hypothesis is expected to die on it, since \
                       SelfIntersectionReport is a separate type the readiness \
                       report does not fold in — registered as a prediction \
                       anyway, because an expectation that is not written down \
                       before the count is taken is not evidence.",
        records: &[
            "method",
            "precondition",
            "required",
            "readiness_field",
            "covered",
        ],
    },
    Preregistration {
        id: "P-19",
        ticket: "S-009",
        hypothesis: "The on-demand-versus-batch crossover for the generalized \
                     winding number is set by the batch path's one-ray-per-grid-row \
                     sharing, not by its point count: on a nearly-closed mesh the \
                     crossover query count Q* is of order N² rather than N³, \
                     landing within 0.5x to 4x of N² for grids from 17³ to 65³. \
                     The naive expectation -- that an on-demand field wins below \
                     N³ queries, since that is how many the batch path answers -- \
                     is registered here as the thing this is expected to beat.",
        falsified_by: "Q* reaching a tenth of N³ or more on a nearly-closed mesh, \
                       which would put the cost in the per-point boundary-edge \
                       correction rather than in the per-row ray, and mean the \
                       row sharing M-299 identified is not what sets the \
                       crossover. Also falsified if no crossover exists at all \
                       above the point where batching is trivially cheaper -- a \
                       single query -- which would say the on-demand field has no \
                       regime and S-009 should be closed rather than built.",
        records: &[
            "samples_per_axis",
            "triangles",
            "boundary_edges",
            "batch_total_ns",
            "on_demand_ns_per_query",
            "crossover_queries",
            "crossover_over_n_squared",
            "crossover_over_n_cubed",
        ],
    },
    Preregistration {
        id: "P-20",
        ticket: "R-010",
        hypothesis: "Splitting a weld on a caller-supplied key moves no topology \
                     metric relative to the unconditional weld beyond the splits \
                     the key itself names: with a constant key every metric is \
                     identical, and with a varying key the only change is the \
                     vertex count rising by the number of sub-classes the key \
                     creates. Across all eight reference fields and every \
                     extractor.",
        falsified_by: "Any topology metric -- non-manifold edges, non-manifold \
                       vertices, boundary edges -- moving where the key is \
                       constant, which would mean the hook itself rather than \
                       the key is doing something. Or, with a varying key, a \
                       non-manifold vertex count that rises by more than the \
                       split count, which would be E*4's failure reappearing: a \
                       partial refusal within one coincidence class leaving its \
                       representative a bowtie.",
        records: &[
            "field",
            "extractor",
            "key",
            "vertices_after",
            "splits",
            "non_manifold_edges",
            "non_manifold_vertices",
            "boundary_edges",
        ],
    },
    Preregistration {
        id: "P-21",
        ticket: "R-024",
        hypothesis: "A freshly extracted mesh separates exactly the sample pairs \
                     the field's own sign separates. For every 6-adjacent pair of \
                     grid samples, the mesh crosses the segment between them an \
                     ODD number of times when the two samples straddle the \
                     surface and an EVEN number when they do not -- so the \
                     connected components of the air sublevel set and the \
                     components of the same samples under mesh-cut adjacency \
                     agree in count and in partition. Marching Cubes achieves \
                     this on all eight reference fields, because its vertex is \
                     the root of the interpolant along the very edge being \
                     probed. At least one dual method does not, because it \
                     places its vertex by solve.",
        falsified_by: "Universal agreement -- every extractor sealing every field \
                       at every resolution, which would be a stronger \
                       correctness statement than this crate currently makes and \
                       is worth saying so. Separately falsified, and more \
                       interestingly, by Marching Cubes disagreeing: that would \
                       put the defect in the primal path, where the crossing is \
                       on the interpolant by construction, and would mean the \
                       failure is in triangulation rather than in vertex \
                       placement.",
        records: &[
            "field",
            "extractor",
            "samples_per_axis",
            "field_air_components",
            "mesh_air_components",
            "unsealed_walls",
            "spurious_walls",
            "mixed_regions",
        ],
    },
    Preregistration {
        id: "P-22",
        ticket: "T-026",
        hypothesis: "Two clauses about Grosso & Zint's mean-ratio triangle \
                     quality, q = 4*sqrt(3)*A / sum(l_i^2). (1) \
                     `marching_cubes` and `marching_cubes+decider` measure the \
                     IDENTICAL mean ratio on every reference field, because both \
                     place their vertices at the root of the interpolant along a \
                     grid edge and a face rule changes which crossings are \
                     joined rather than where they are. That is the paper's own \
                     explanation for its MC and TMC columns agreeing to two \
                     decimals on all seven of its rows. (2) This crate's \
                     Marching Cubes lands inside 0.65 to 0.71 on a smooth \
                     analytic field, the band their MC occupies, whose gen2 \
                     figure is resolution-independent at 64, 128 and 256 cubed.",
        falsified_by: "The two Marching Cubes entries differing at all, which \
                       would mean the face rule moves geometry and not only \
                       connectivity -- and would contradict the mechanism the \
                       paper states for its own two columns. Or a mean ratio \
                       outside 0.65-0.71 on a smooth analytic field, which would \
                       mean the metric is not implementation-independent and the \
                       published baseline cannot be compared against at all. The \
                       second is a real possibility rather than a formality: \
                       their MC is somebody else's code measured on somebody \
                       else's fields, and every cross-source comparison this \
                       repo has attempted so far has needed an amendment.",
        records: &[
            "field",
            "extractor",
            "samples_per_axis",
            "mean_ratio",
            "irregular_vertices",
            "referenced_vertices",
            "triangles",
        ],
    },
    Preregistration {
        id: "P-23",
        ticket: "R-022a",
        hypothesis: "Repairing the air sublevel set's connectivity after a \
                     brush dig costs work proportional to the DIRTY SET and not \
                     to the lattice. Concretely: at a fixed brush radius, the \
                     number of union operations an incremental insertion-only \
                     update performs is constant as the lattice grows through \
                     33, 65 and 129 cubed, so unions divided by newly-air \
                     samples stays inside a narrow band and never exceeds 6 -- \
                     the degree of the lattice -- while a full rebuild's union \
                     count grows as n cubed. Digging only ever inserts, and an \
                     insert joins at most two trees with no replacement-edge \
                     search, so a union-find is the entire structure.",
        falsified_by: "The incremental union count growing with n cubed at a \
                       fixed brush size. That is R-022's own stated falsifier \
                       and it would mean edit-proportional repair is \
                       unavailable even in the easy direction, which closes the \
                       whole direction rather than only this half. Separately \
                       falsified by unions per dirty sample exceeding 6, which \
                       would mean the harness is visiting something other than \
                       the six incident lattice edges and is measuring its own \
                       traversal rather than the repair.",
        records: &[
            "samples_per_axis",
            "dirty_samples",
            "incremental_unions",
            "rebuild_unions",
            "unions_per_dirty",
        ],
    },
    Preregistration {
        id: "P-24",
        ticket: "R-023",
        hypothesis: "The trilinear body-saddle value F(s) is a DECISION MARGIN \
                     for the interior ambiguity, not merely correlated with it: \
                     sign(F(s)) agrees with Interior::test() on every ambiguous \
                     cell whose sweep has NO pole in (0, 1), and every \
                     disagreement between them has a pole inside the sweep. A \
                     pole is the term Chernyaev's test drops and is what \
                     Custodio's Figure 6 counterexample is built from, so that \
                     is where the published algorithms part company and is the \
                     only place H permits a difference. Thresholding |F(s)| at \
                     epsilon is then a one-parameter family whose epsilon = 0 \
                     member is the published decider exactly.",
        falsified_by: "A disagreement on an ambiguous cell whose sweep has no \
                       pole in (0, 1). Chernyaev's quadratic is exact there, so \
                       a difference would mean the body-saddle value is not the \
                       quantity the interior test is a sign of, and the whole \
                       reframing collapses rather than needing a tolerance. \
                       Separately falsified if the two agree on every cell \
                       INCLUDING the poled ones, which would mean Custodio's \
                       correction is unreachable on this crate's fields and the \
                       census cannot discriminate -- a null that says the \
                       fixture is wrong rather than the hypothesis.",
        records: &[
            "field",
            "samples_per_axis",
            "ambiguous_cells",
            "agreements",
            "disagreements",
            "disagreements_with_pole",
            "disagreements_without_pole",
        ],
    },
    Preregistration {
        id: "P-25",
        ticket: "R-022b",
        hypothesis: "Repairing the air sublevel set's connectivity after a \
                     brush FILL costs work proportional to the SHED VOLUME -- \
                     the air samples that leave their component -- and not to \
                     the surviving component, nor to the lattice. Concretely: \
                     at a fixed brush radius the voxels the replacement search \
                     visits stays flat as the lattice grows through 33, 49 and \
                     65 cubed, while a rebuild's visit count grows as n cubed. \
                     Two mechanisms, both already measured. A union-find CAN \
                     absorb deletion, because a parent pointer only has to \
                     reach the right root and a filled sample sitting mid-tree \
                     is never queried -- so only the shed pieces are re-rooted \
                     and the surviving side is never walked. And the shed \
                     pieces are tiny: M-320 measured the smaller side of a \
                     split at ONE voxel at the median and 120 at the observed \
                     maximum, against 227,567 air samples. The levelled HDT \
                     scheme is therefore unnecessary here; lockstep search \
                     outward from every seed, stopping when all but one \
                     frontier exhausts, is enough.",
        falsified_by: "Visited voxels growing as n cubed at a fixed brush size, \
                       which would mean the search is exploring the SURVIVING \
                       component rather than the shed pieces -- i.e. the \
                       lockstep stop condition is wrong and the structure is \
                       walking the thing it was built to avoid walking. \
                       Separately and more seriously falsified by any \
                       disagreement between the incrementally maintained \
                       components and a full rebuild over the same values: \
                       component count, or any connected() answer. That would \
                       mean the structure is fast and WRONG, which is worse \
                       than slow and right, and it is the failure a \
                       measurement of cost alone cannot see.",
        records: &[
            "samples_per_axis",
            "fills",
            "dirty_samples",
            "seeds",
            "visited",
            "splits",
            "shed_components",
            "vanished_components",
            "rebuild_visited",
        ],
    },
];

/// `a == b`, in a const context.
///
/// `str`'s `PartialEq` is not const, and neither is `<[u8]>::eq`. A byte loop is
/// — and it only has to be correct, not fast, since it runs at compile time over
/// a handful of six-character ids.
const fn str_eq(a: &str, b: &str) -> bool {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    if a.len() != b.len() {
        return false;
    }
    let mut i = 0;
    while i < a.len() {
        if a[i] != b[i] {
            return false;
        }
        i += 1;
    }
    true
}

/// Is `id` registered?
///
/// `const`, because [`experiment!`](macro@crate::experiment) asserts on it at compile
/// time — that assertion is the entire ticket.
#[must_use]
pub const fn is_preregistered(id: &str) -> bool {
    let mut i = 0;
    while i < PREREGISTERED.len() {
        if str_eq(PREREGISTERED[i].id, id) {
            return true;
        }
        i += 1;
    }
    false
}

/// The registration for `id`.
///
/// # Panics
///
/// If `id` is not registered. In the `const` context
/// [`experiment!`](macro@crate::experiment) uses, that panic is a **compile error**,
/// which is the ticket's acceptance; called at run time it is a programming
/// error with the same cause.
#[must_use]
pub const fn preregistration(id: &str) -> &'static Preregistration {
    let mut i = 0;
    while i < PREREGISTERED.len() {
        if str_eq(PREREGISTERED[i].id, id) {
            return &PREREGISTERED[i];
        }
        i += 1;
    }
    panic!("no such pre-registration; add it to experiment::PREREGISTERED first")
}

/// The pre-registration for `$id`, or a **compile error** if there is none.
///
/// R-000's acceptance, in one line: *"an experiment without a pre-registered
/// `P-` fails to build."*
///
/// ```compile_fail
/// // "P-999" is not registered, so this does not compile.
/// let _ = isomesh::experiment!("P-999");
/// ```
#[macro_export]
macro_rules! experiment {
    ($id:literal) => {{
        // The gate. A `const` block is evaluated at compile time whether or not
        // the value is used, so an unregistered id cannot reach a test run.
        const _CHECK: () = assert!(
            $crate::experiment::is_preregistered($id),
            "this experiment id is not pre-registered — add it to \
             experiment::PREREGISTERED, in its own commit, before running \
             anything"
        );
        $crate::experiment::preregistration($id)
    }};
}
