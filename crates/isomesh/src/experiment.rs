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
    Preregistration {
        id: "P-26",
        ticket: "R-022b",
        hypothesis: "P-25 with its MECHANISM clause replaced and its COST \
                     clause unchanged; see the falsification in FINDINGS as \
                     cross-26. Unchanged: repairing the air sublevel set after \
                     a brush fill costs work proportional to the SHED VOLUME, \
                     not to the surviving component nor to the lattice, so \
                     visited voxels stay flat through 33, 49 and 65 cubed at a \
                     fixed brush radius. Replaced: the structure is not a \
                     union-find but a FLAT label array -- every sample carries \
                     its component id directly, so re-rooting a shed piece is \
                     one write per member and no surviving sample can route \
                     through it. Flat labels fix the REPRESENTATION, not the \
                     SEARCH: the lockstep replacement search was always \
                     required and the union-find merely promised falsely that \
                     it could be skipped. Added: lockstep bounds work by the \
                     SECOND-LARGEST piece, so M-320's one-voxel median is a \
                     property of the edit distribution rather than of the \
                     structure. Bisecting a tunnel between two equal caverns \
                     makes both frontiers huge and visited then grows with n at \
                     a fixed brush size. So the prediction is stated per \
                     fixture: FLAT on the measured distribution, GROWING on a \
                     deliberate bisect.",
        falsified_by: "Visited voxels growing as n cubed ON THE MEASURED \
                       DISTRIBUTION at a fixed brush size. Growth on the bisect \
                       fixture is predicted and is not falsifying; a structure \
                       that came out flat on BOTH would instead mean the bisect \
                       fixture is not adversarial and needs rebuilding, which \
                       is a fixture failure rather than a result. Separately \
                       and more seriously falsified, as in P-25, by any \
                       disagreement between the maintained components and a \
                       full rebuild over the same values -- component count, or \
                       any connected() answer. Fast and wrong is worse than \
                       slow and right.",
        records: &[
            "samples_per_axis",
            "fixture",
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
    Preregistration {
        id: "P-27",
        ticket: "R-029",
        hypothesis: "Offsetting a slab's mid-plane by half a voxel takes the \
                     count of lattice samples whose default SDF gradient is \
                     exactly [0,0,0] from one full lattice plane (65-squared = \
                     4,225 at 65 samples per axis) to ZERO, while the \
                     medial-stability population -- sub-voxel probes with \
                     voxel-step central-difference gradient magnitude under \
                     0.1 -- changes by less than 5 percent. Two populations \
                     on purpose: the sub-threshold set is a band of thickness \
                     0.2h, thinner than the lattice pitch, so counted on the \
                     voxel lattice itself clause two would be vacuously false \
                     under ANY misalignment. On the registered probe lattice \
                     of pitch h/200 the worst change a rigid offset can \
                     produce is 1/40 = 2.5 percent -- derived, half the \
                     registered bound, and the margin between them is the \
                     room for the discrete band to misbehave.",
        falsified_by: "Exact zeros surviving the half-voxel offset -- the \
                       discrete gradient is doing something the continuous \
                       identity does not describe, and R-030 inherits that \
                       question before anything else runs. Clause two is \
                       separately falsified by the band count moving more \
                       than 5 percent, which the derived 2.5 percent \
                       instrument bound says can only happen if the discrete \
                       band is not the rigid 0.2h ramp the arithmetic \
                       assumes.",
        records: &[
            "arm",
            "offset_voxels",
            "exact_zeros",
            "band_count",
            "probe_pitch_voxels",
        ],
    },
    Preregistration {
        id: "P-28",
        ticket: "R-030",
        hypothesis: "The matched-analytic re-instrumentation of the identity \
                     r = rho * sqrt(1 - |grad rho|^2), after V-46 showed the \
                     ticket's MEB-oracle shape cannot discriminate. Three \
                     clauses. C1, form: on slab, wedge and triangular-prism \
                     air fixtures in generic position -- piecewise-linear \
                     fields, so the voxel-step mollified truth is derivable \
                     exactly -- the measured r matches the derived closed \
                     form within 1e-9 of the gap on 100 percent of \
                     medial-band samples, exercising the two-point AND \
                     three-point closest-point cases. C2, curvature: on a \
                     capsule at generic off-axis samples, where the true \
                     inscribed radius is zero, the formula's own noise floor \
                     has a world-unit median that HALVES per resolution \
                     doubling -- 33-to-129 end-to-end ratio at most 0.35 -- \
                     the O(h) sqrt-amplification of the O(h squared) \
                     curvature error in the discrete gradient. C3, clearance: \
                     for slab gaps of 3, 6 and 10 voxels across 8 sub-voxel \
                     phases, the band-max r sits inside the derived envelope \
                     [sqrt(3)/2 * (W - h/2), W] on 24 of 24 rows.",
        falsified_by: "C2's end-to-end ratio at or above 0.7 -- an \
                       h-independent noise floor, meaning the discrete score \
                       cannot separate medial signal from curvature noise at \
                       any fixed world scale, and Calibre, the throat metric \
                       and handholds die together as the dossier said. C1 \
                       failing instead is an implementation or transcription \
                       finding, not a verdict on the identity. C3 failing \
                       means the clearance-envelope derivation is wrong and \
                       the lambda test loses its accuracy claim. And the \
                       wrong-form inversion -- rho times (1 - |grad rho|) -- \
                       must fail C1 on at least 30 percent of mid-band \
                       samples by more than a tenth of the gap, or this \
                       instrument has not been shown able to go red and V-46 \
                       applies to it too.",
        records: &[
            "fixture",
            "samples_per_axis",
            "band_samples",
            "within_tol_pct",
            "band_median_residual_world",
            "clearance_true_voxels",
            "clearance_est_voxels",
            "clamped",
        ],
    },
    Preregistration {
        id: "P-29",
        ticket: "R-031",
        hypothesis: "Dreybrodt and Gabrovsek's wormhole competition reproduces \
                     on a 64x64 fracture lattice with their own constants, \
                     read from the converted primary source rather than \
                     tuned: cubic-law resistance, linear kinetics \
                     F = k(1 - c/ceq) with the composite k(a), \
                     k1 = 4e-11 mol/cm2/s, ceq = 1e-6 mol/cm3, \
                     D = 1e-5 cm2/s, penetration-length transport, \
                     da/dt = 2*gamma*F. Clause one: under constant-head \
                     boundaries the post-breakthrough aperture distribution \
                     is BIMODAL by the central-gap statistic -- max gap in \
                     ln(a) with 1 percent tails dropped, at least 0.2, the \
                     registered initial log-sd -- on both the \
                     seeded-homogeneous and lognormal-heterogeneous nets, \
                     while the fixed-flux recharge-limited arm stays \
                     UNIMODAL at matched cumulative dissolved volume, per \
                     Perne, Covington and Gabrovsek's limited-recharge \
                     suppression. Clause two: past breakthrough, more than \
                     90 percent of dissolution flux concentrates in fewer \
                     than 10 percent of edges. Recorded, not registered: \
                     heterogeneous breakthrough earlier than \
                     seeded-homogeneous (the paper's 560 against 1890 \
                     years, as an ordering, not a magnitude), per-tick \
                     cost, and the Gini-over-flow series.",
        falsified_by: "No bimodal split on either constant-head arm, or \
                       dissolution flux failing to concentrate past \
                       breakthrough -- the positive-feedback premise of the \
                       mechanic dies with the kinetics. Separately falsified \
                       as an INSTRUMENT if the recharge-limited arm goes \
                       bimodal while the detector also calls the t-zero \
                       lognormal apertures bimodal -- both-arms-bimodal \
                       indicts the statistic, both-arms-unimodal indicts \
                       the kinetics. The detector's own red and green -- a \
                       synthetic half-shifted sample it must call bimodal, \
                       the t-zero sample it must call unimodal -- run before \
                       any verdict is read.",
        records: &[
            "arm",
            "ticks",
            "years",
            "breakthrough_years",
            "max_gap_ln",
            "bimodal",
            "flux_top10_pct",
            "gini_flow",
            "max_da_over_a_pct",
            "tick_ms_median",
        ],
    },
    Preregistration {
        id: "P-30",
        ticket: "R-031",
        hypothesis: "P-29's clause one with its INSTRUMENT replaced and its \
                     prediction unchanged; the falsification is M-326. The \
                     central-gap detector trimmed 1 percent tails -- 81 of \
                     8,128 edges -- while a winning wormhole path is at most \
                     one input-to-output chain of 63 edges, so the registered \
                     statistic ate the mode it was looking for. Replaced: no \
                     trim; a split counts as bimodal only when the gap in \
                     sorted ln(a) is at least 0.2 AND both sides hold at \
                     least 8 edges -- 0.1 percent, eight-fold below the \
                     smallest possible winner mode (the lattice's own path \
                     length, 63, knowable before any run) and eight-fold \
                     above single-edge outliers. Unchanged: under constant \
                     heads the post-breakthrough aperture distribution is \
                     bimodal on both the seeded-homogeneous and \
                     lognormal-heterogeneous nets, and the recharge-limited \
                     arm stays unimodal at the same at-breakthrough \
                     dissolved volume. Recorded beside it: the water-flux \
                     top-10-percent share, next to the dissolution share \
                     M-326 closed at 78 percent.",
        falsified_by: "No qualifying gap on either constant-head arm. With \
                       the mode guard in place that verdict would mean the \
                       aperture histogram is genuinely not gapped-bimodal -- \
                       the dossier's tier-R bimodality reading dies, and the \
                       competition claim rests on the concentration columns \
                       M-326 already measured (Gini 0.976 against 0.560). \
                       The positive-feedback premise itself is convicted \
                       only if those concentration columns also regress, \
                       which M-326 shows they do not. Instrument-falsified, \
                       as before, if the recharge arm reads bimodal while \
                       the t-zero green also fails; the synthetic red and \
                       t-zero green run again before any verdict is read.",
        records: &[
            "arm",
            "ticks",
            "years",
            "breakthrough_years",
            "max_gap_ln",
            "guarded_gap_ln",
            "guarded_bimodal",
            "flux_top10_pct",
            "flux_water_top10_pct",
            "gini_flow",
            "max_da_over_a_pct",
            "tick_ms_median",
        ],
    },
    Preregistration {
        id: "P-31",
        ticket: "R-032",
        hypothesis: "On 20 seeded dug scenes -- a solid block carved by 12 \
                     random overlapping capsule brushes through the crate's \
                     own BrushStack composition -- the weak feature size, \
                     measured as the minimum air-side distance-to-boundary \
                     over discrete critical points (voxel-step \
                     central-difference gradient magnitude under 0.5, the \
                     dossier's own theta-above-120-degrees filter constant, \
                     non-maximum-suppressed over 26-neighbourhoods), is \
                     below 2 voxels on MORE THAN 80 PERCENT of scenes -- so \
                     the homotopy certificate lambda < wfs essentially \
                     never holds at brush scale and the lambda-medial line \
                     rests on Hausdorff stability instead. A single-cavity \
                     control scene, one 20-voxel sphere in generic \
                     position, must report wfs of at least 10 voxels, \
                     demonstrating the instrument can call the certificate \
                     AVAILABLE before it is trusted calling it absent.",
        falsified_by: "wfs at or above 2 voxels on half or more of the dug \
                       scenes -- the certificate comfortably available at \
                       brush scale, which would make the stronger homotopy \
                       guarantee live and the Hausdorff fallback \
                       unnecessary. Instrument-falsified by any dug scene \
                       reporting zero critical points (a minimum over an \
                       empty set is not a measurement), or by the control \
                       cavity failing its 10-voxel floor.",
        records: &[
            "scene",
            "air_samples",
            "critical_points",
            "wfs_voxels",
            "epsilon",
            "samples_per_axis",
        ],
    },
    Preregistration {
        id: "P-32",
        ticket: "R-033",
        hypothesis: "A one-voxel interior edit does not move a carved \
                     pillar's fundamental eigenvalue audibly. The pitch JND \
                     is 0.6 percent on FREQUENCY and f is proportional to \
                     sqrt(lambda), so the audibility threshold on lambda-one \
                     is 1.2 percent -- the conversion is registered because \
                     applying 0.6 to lambda directly would be a silent \
                     factor of two. Clause one: all 8 seeded \
                     strictly-interior one-voxel digs move lambda-one by \
                     less than 1.2 percent. Clause two: at least 1 of 4 digs \
                     adjacent to the fixture's two-cell web exceeds 1.2 \
                     percent -- the ticket's 'only edits near a thin feature \
                     are audible', both sides of it. Clause three, the \
                     null's reachability: carving a 20-percent-volume cavity \
                     moves lambda-one by more than 15 percent. Instrument: \
                     hexahedral FEM assembled directly on the occupancy grid \
                     (trilinear, 2x2x2 Gauss, E=1 nu=0.3 rho=1, lumped mass, \
                     base layer fixed), matrix-free inverse power iteration \
                     (48 outer) over Jacobi-preconditioned CG (256 inner, \
                     warm-started), and every reported eigenvalue carries an \
                     a-posteriori certificate -- the residual Kx minus \
                     lambda Mx in the M-inverse norm, at or below 5e-4 of \
                     lambda -- so fixed iteration counts can never mean \
                     silently wrong.",
        falsified_by: "Any strictly-interior dig moving lambda-one at or \
                       above 1.2 percent -- per-edit modal audio earns its \
                       ticket, bounded in advance by Picard's O(m^2.8) \
                       scaling at about 40 modes in a 5 ms budget. Clause \
                       two failing the other way -- no web-adjacent dig \
                       audible -- closes the modal direction entirely at the \
                       cost of a day, the outcome the dossier priced. \
                       Instrument-falsified by the certificate exceeding its \
                       bound, the two deterministic starts disagreeing \
                       beyond twice the certificate, the control cavity \
                       moving lambda-one by less than 15 percent, or the \
                       free single-element spectrum holding other than \
                       exactly six rigid modes.",
        records: &[
            "edit",
            "cells",
            "dof",
            "lambda1_base",
            "lambda1_edited",
            "delta_pct",
            "audible",
            "certificate_rel",
        ],
    },
    Preregistration {
        id: "P-33",
        ticket: "R-034a",
        hypothesis: "A pure-Rust reimplementation of Whiting, Ochsendorf and \
                     Durand's rigid-block feasibility program -- per-vertex \
                     interface forces, equilibrium A f = -w, friction cone \
                     at mu = 0.7, compression-only -- reproduces the only \
                     external ground truths in the dossier. Bisecting \
                     thickness over centerline radius on a 100-block \
                     semicircular arch (the paper's own tessellation) finds \
                     the infeasibility threshold at 0.1075 plus or minus \
                     0.0010 -- Milankovitch's 1907 analytic value, which \
                     Whiting's solver hit at 0.10746 -- and bisecting ground \
                     tilt at t/r = 0.20 finds 15.84 degrees plus or minus \
                     0.05, Ochsendorf's value. The solver is alternating \
                     projection between the equilibrium affine set -- exact \
                     per iteration via one prefactored Cholesky of A times \
                     A-transpose -- and the per-vertex friction cones, whose \
                     projection is closed-form; compression is exact in the \
                     cone rather than penalized; 20,000 fixed iterations; \
                     the decision reads the CONE-side iterate's equilibrium \
                     residual per unit weight, 1e-5 feasible and 1e-4 \
                     infeasible, the band between asserted never hit, and \
                     the thresholds themselves are checked by the golden \
                     values -- a tuned-wrong threshold cannot hit 0.1075 \
                     from both sides. Block weights act at exact \
                     annular-sector centroids, because centerline weights \
                     reproduce Heyman's 0.106 rather than Milankovitch's \
                     0.1075 and the third decimal is the whole point. \
                     Recorded, not registered: the threshold at 25, 50, 100 \
                     and 200 blocks -- the paper's own warning that coarser \
                     blocks over-estimate stability, so coarse thresholds \
                     should sit BELOW fine ones.",
        falsified_by: "Missing either golden value -- the solver is wrong \
                       somewhere between the formulation and the arithmetic, \
                       and nothing structural may be built on it. \
                       Instrument-aborts before any verdict: doubling \
                       gravity moving any feasibility decision (the program \
                       is scale-invariant or it is wrong), or any bisection \
                       step landing in the undecided residual band. The \
                       game-facing rule and the warm-start economics are \
                       R-034b's, deliberately not here.",
        records: &[
            "test",
            "blocks",
            "value",
            "target",
            "abs_error",
            "within_tolerance",
            "residual_feasible",
            "residual_infeasible",
        ],
    },
    Preregistration {
        id: "P-34",
        ticket: "R-034b",
        hypothesis: "Warm-started re-solves of the M-330-validated feasibility \
                     program are BIMODAL over an edit corpus, in the ticket's \
                     two classes. Fixture: a running-bond masonry wall (8 \
                     courses, ~96 blocks, bed and head joints, mu = 0.7) -- \
                     chosen over the arch because redundancy is what lets a \
                     severing edit leave a standing structure, and M-330 \
                     showed the arch too simple to discriminate anything \
                     block-structural. Cost is COUNTED, not timed: iterations \
                     of the alternating projection until the cone-side \
                     residual first crosses the 1e-5 feasibility line, cold \
                     (from zero) versus warm (from the pre-edit solution, \
                     mapped by interface identity). Clause one: over 10 \
                     non-severing edits (single-block weight nudges, small \
                     gravity tilts) the median warm-to-cold ratio is at most \
                     0.15 -- under 15 percent of cold, the ticket's number. \
                     Clause two: over 10 severing edits (an interior block \
                     removed, forces rerouted around the hole) the median \
                     ratio is at least 0.5 -- under 2x speedup. The \
                     bimodality IS the prediction: the two class medians \
                     separated by more than 3x.",
        falsified_by: "A unimodal ratio distribution -- class medians within \
                       3x of each other -- which kills the cheap-incremental \
                       story and demotes the admissibility gate to \
                       background-budget-only, exactly as the original \
                       R-034 registered. Instrument notes: any edit that \
                       classifies infeasible at the 20,000-iteration cap is \
                       recorded as collapsed and excluded with its count \
                       printed -- a corpus that mostly collapses is a \
                       fixture failure and aborts; and every surviving edit \
                       must reach the decision line within the cap or the \
                       count, not the clock, has failed to decide.",
        records: &[
            "edit",
            "class",
            "iters_cold",
            "iters_warm",
            "ratio",
            "feasible",
        ],
    },
    Preregistration {
        id: "P-35",
        ticket: "R-034b",
        hypothesis: "P-34 with its instrument floor removed and its corpus \
                     unmixed; the falsification is M-331. Corrections, all \
                     derivable before running: the decision is probed every \
                     ITERATION (the floor becomes 1/cold instead of \
                     10/cold); the wall grows to 20 courses of 24 blocks \
                     (~470 blocks) so the cold count has three digits of \
                     dynamic range; and tilts leave the registered corpus \
                     -- a gravity rotation moves every interface force and \
                     is recorded as its own class, not averaged into local \
                     edits. Unchanged prediction, the ticket's own: over 10 \
                     single-block weight nudges the median warm-to-cold \
                     ratio is at most 0.15; over 10 interior removals it is \
                     at least 0.5; the class medians separate by more than \
                     3x. M-331's floor-pinned weight nudges and 0.8 removal \
                     median already point this way at 4x, undecidably.",
        falsified_by: "The weight-nudge median still above 0.15 with the \
                       floor at 1/cold -- then the 0.400 was never the \
                       instrument's, the cheap-incremental story dies on \
                       merit, and the admissibility gate is \
                       background-budget-only as the original R-034 \
                       registered. Or medians within 3x -- unimodal, same \
                       consequence. Collapsed removals excluded on the \
                       record as in P-34, aborting past three.",
        records: &[
            "edit",
            "class",
            "iters_cold",
            "iters_warm",
            "ratio",
            "feasible",
        ],
    },
    Preregistration {
        id: "P-36",
        ticket: "R-035b",
        hypothesis: "On M-333's verified substrate -- the 64-cubed gyroid \
                     chunk, 12,615 Surface Nets vertices, heat operator at \
                     t = h-bar squared, nested-dissection ordering -- a \
                     radius-4-voxel surface perturbation changes at most 400 \
                     vertex slots (M-318's 346 of 15,706, keyed by cell \
                     identity; the 1-ring operator-row halo is counted \
                     beside it, expected 3 to 6 times more), and a partial \
                     refactorization over the elimination-tree ancestor \
                     closure of the changed columns re-establishes a valid \
                     factorisation at least 20 TIMES cheaper than a full \
                     refactorisation by wall time, with the flop ratio at \
                     least 10 alongside -- a large time ratio over a small \
                     flop ratio is an implementation artifact, not a result. \
                     The absolutes (the ticket's under-5-ms and over-100-ms) \
                     are recorded, not load-bearing: M-333 already measured \
                     the full refactor at 87.7 ms. Validity is asserted, not \
                     assumed: the updated factor holds the same Frobenius \
                     residual bound as a fresh one, its solve agrees with \
                     the refactored solve within 1e-8 relative, and a \
                     deliberately skipped closure column must push the \
                     residual past its bound before any verdict is read. \
                     Timing is interleaved in both orders, 11 repetitions, \
                     medians.",
        falsified_by: "The update under 10 times cheaper by wall or under \
                       the flop-ratio floor -- the prefactored family is \
                       dead for live carving and everything \
                       surface-intrinsic routes to the Closest Point Method \
                       instead, the routing decision this ticket exists to \
                       make. Separately: the slot count exceeding 400 \
                       re-scopes M-318's extrapolation; the slot SET \
                       changing at all aborts the fixture (the experiment \
                       is about value updates on a stable pattern, and says \
                       so); and the skipped-column inversion failing to go \
                       red voids the validity oracle and the run with it.",
        records: &[
            "rep",
            "order",
            "update_ms",
            "refactor_ms",
            "update_flops",
            "refactor_flops",
            "changed_slots",
            "changed_rows",
            "closure_rows",
        ],
    },
    Preregistration {
        id: "P-37",
        ticket: "R-036",
        hypothesis: "The ticket's premise is half false, and the correction \
                     is the first finding: the tracker maintains component \
                     volume (Air::component_size) and does NOT maintain \
                     boundary surface area -- nothing in the crate does. \
                     With the accumulator added (per-label air-solid face \
                     counts, delta-maintained through build, dig's blob \
                     growth, fill's retirement, merge transfer and split \
                     hand-off; domain-boundary faces count as solid, the \
                     sealed-box convention; AirWorld roll-up deliberately \
                     out of scope), clause one: a Sabine RT60 for the \
                     breach-frame component -- 0.161 times volume over \
                     absorption times area, two accumulator reads and a \
                     divide -- costs under 0.1 ms, and structurally so. \
                     Clause two: a Planeverb-style 2D FDTD re-bake of a \
                     64x64 slice (Rosen, Godin and Raghuvanshi, \
                     10.1111/cgf.14099, public C++ reference; leapfrog \
                     pressure-velocity, 1,000 steps -- about half a second \
                     of audio at the CFL step, the length a decay \
                     measurement needs -- damped edges, cost is the claim \
                     and acoustic fidelity is not) completes in under 30 ms \
                     single-threaded. Recorded, not registered: dig and \
                     fill costs with the accumulator in place, and the \
                     split rate against M-319's one-in-six -- a divergence \
                     THERE would be news; the clauses holding is not.",
        falsified_by: "Either figure exceeding its bound by 3x, the \
                       ticket's own falsifier. Instrument notes: the area \
                       invariant -- a full recount equals the maintained \
                       counts, and the label-free global face total equals \
                       their sum -- is asserted in the crate's own tests \
                       over synchronous op sequences, with a deliberate \
                       corruption shown to turn the checker red; \
                       budget-truncated ops leave area conservatively stale \
                       exactly as labels already are, scoped and documented \
                       rather than solved.",
        records: &["quantity", "value", "unit", "bound", "held"],
    },
    Preregistration {
        id: "P-38",
        ticket: "R-037",
        hypothesis: "A-024's odd-row-stride remedy reached DualMesher and \
                     never reached Marching Cubes, and the defect is still \
                     there: MarchingCubes indexes its private `values` buffer \
                     through `shape.linearize`, so at 128 samples per axis the \
                     row stride is 512 bytes and the plane stride exactly \
                     65,536, and its `edge_vertices` cache -- 3 u32 per sample, \
                     the buffer with the scattered access -- has a plane \
                     stride of three times that. Clause one: measured on a \
                     sphere, Marching Cubes at 128 costs more than 1.5x the \
                     mean of its 127 and 129 neighbours per sample, where the \
                     dual path (already fixed) does not. Clause two: giving \
                     both private buffers the same `size[0] | 1` stride takes \
                     that ratio under 1.1 and is worth at least 1.25x at 128 \
                     alone. Clause three: every one of the 216 golden hashes \
                     is bit-identical afterwards, structurally -- the change \
                     permutes where floats are stored, not which floats are \
                     computed, nor the order cells are visited, nor the order \
                     vertices are created.",
        falsified_by: "The 128 ratio coming in under 1.5x before the change, \
                       which would mean Marching Cubes' access pattern does \
                       not alias and A-024's remedy is specific to the dual \
                       path -- a narrower result than the repo currently \
                       believes, and the more interesting one. A moved golden \
                       hash falsifies clause three outright and makes this a \
                       behaviour change rather than a layout change.",
        records: &[
            "extractor",
            "samples_per_axis",
            "ns_per_sample_before",
            "ns_per_sample_after",
            "neighbour_ratio_before",
            "neighbour_ratio_after",
            "golden_unchanged",
        ],
    },
    Preregistration {
        id: "P-39",
        ticket: "R-038",
        hypothesis: "A brush that provably cannot win the min/max chain inside \
                     a chunk can be deleted from the tape before the chunk is \
                     meshed, and deleted BIT-EXACTLY. The bound is one sample \
                     per brush per chunk: a shape with declared Lipschitz \
                     constant l varies by at most l*r over a box of \
                     circumradius r, so f(centre) +/- l*r encloses it. The \
                     exactness is a selection argument, not an approximation \
                     one -- `apply(Add)` is IEEE min and `apply(Subtract)` is \
                     IEEE max, both of which SELECT an operand rather than \
                     computing a new value, and negation is exact, so dropping \
                     a provably-losing Add or Subtract moves the result by \
                     zero ULP. The asymmetry is registered rather than \
                     discovered: smooth_min is exactly dominant at h == 0 and \
                     is NOT bit-exactly prunable in the losing direction, \
                     because at h == 1 it returns b + (a - b), which is not \
                     bit-identical to a. Clause one: on a 64-brush stack of \
                     Add and Subtract spheres and capsules scattered over a \
                     4x4x4 chunk world, the median surviving-brush fraction \
                     per chunk is under 0.5. Clause two: meshing a chunk \
                     against the pruned tape is at least 1.25x meshing it \
                     against the full one. Clause three: every chunk's mesh is \
                     byte-identical between the two.",
        falsified_by: "A survivor fraction at or near 1.0, meaning there is \
                       nothing to prune on a scattered stack and the mechanism \
                       has no purchase on this workload whatever the paper \
                       reports. Or -- the result that would matter far more -- \
                       any byte difference between pruned and full output on \
                       an Add/Subtract-only stack, which would refute the IEEE \
                       selection lemma the whole mechanism rests on.",
        records: &[
            "brushes",
            "chunks",
            "survivor_fraction_median",
            "survivor_fraction_max",
            "ns_per_sample_full",
            "ns_per_sample_pruned",
            "speedup",
            "mesh_identical",
        ],
    },
    Preregistration {
        id: "P-40",
        ticket: "R-039",
        hypothesis: "The active-cell test is one bit per sample, so 64 cells \
                     decide at once. Packing `value < 0` into a u64 bitmap \
                     along x and fusing the four rows that bound a cell row -- \
                     any = OR of (w | w>>1), all = AND of (w & w>>1), active = \
                     any & !all -- replaces DualMesher::place_vertices' \
                     eight-corner gather on the ~97% of cells that produce \
                     nothing. Clause one: on a surface-free 128-cubed field, \
                     where the active path never runs, the traversal stage is \
                     at least 2x faster. Clause two: on sphere at 128-cubed the \
                     WHOLE extractor is at least 1.25x faster. Clause three: \
                     the mesh is byte-identical, because the bit is the same \
                     comparison and the set-bit walk visits cells in the same \
                     lexicographic order, so vertex creation order is \
                     unchanged. The sign bit is deliberately NOT used as the \
                     inside bit: -0.0 has it set while -0.0 < 0.0 is false, and \
                     box_exact is exactly zero across its whole boundary, so \
                     signed zeros are reachable in this crate's own fixtures.",
        falsified_by: "The stage ratio under 2x on the surface-free field, \
                       which would mean the scalar gather was not the cost. Or \
                       the whole-extractor ratio under 1.25x on sphere at 128, \
                       which puts the mechanism under the Amdahl floor that \
                       M-296's field-evaluation dominance sets and makes the \
                       bitmap not worth its own bookkeeping. Any mesh \
                       difference falsifies the ordering argument outright.",
        records: &[
            "field",
            "samples_per_axis",
            "active_fraction",
            "stage_ns_scalar",
            "stage_ns_bitmap",
            "stage_ratio",
            "extract_ns_scalar",
            "extract_ns_bitmap",
            "extract_ratio",
            "mesh_identical",
        ],
    },
    Preregistration {
        id: "P-41",
        ticket: "R-040",
        hypothesis: "The sign lattices this crate meshes are not \
                     well-composed, and that is where its dual extractors go \
                     non-manifold. A digital set is well-composed exactly when \
                     its boundary is a 2-manifold, which Latecki characterises \
                     by the absence of two critical configurations: a diagonal \
                     pair sharing only an edge (2D-critical) and a diagonal \
                     pair sharing only a vertex (3D-critical). Clause one: \
                     over the eight reference fields at 65 samples per axis, \
                     the count of cells whose 2x2x2 sign neighbourhood hosts a \
                     critical configuration is non-zero on at least \
                     noise_cavity. Clause two: at least 90% of the \
                     non-manifold edge and vertex incidents that \
                     validate::MeshReport already reports for DualContouring \
                     and SurfaceNets occur in cells the census flagged \
                     critical. This registers the DETECTOR only. The repair is \
                     deliberately not registered: it moves the surface by up \
                     to a cell and breaks every golden hash, so it is only \
                     worth designing if clause two holds.",
        falsified_by: "A critical-configuration count of zero on all eight \
                       fields, meaning the lattice is already well-composed \
                       and the non-manifold output comes from somewhere else \
                       entirely -- the QEF vertex escaping its cell being the \
                       obvious other suspect. Or co-location below 90%, \
                       meaning the sign lattice is not the cause and a repair \
                       would be treating the wrong object.",
        records: &[
            "field",
            "samples_per_axis",
            "cells",
            "critical_2d_cells",
            "critical_3d_cells",
            "extractor",
            "non_manifold_edges",
            "non_manifold_vertices",
            "incidents_in_critical_cells",
            "colocation_fraction",
        ],
    },
    Preregistration {
        id: "P-42",
        ticket: "R-041",
        hypothesis: "Curvature computed as a normal-cycle MEASURE is additive \
                     over chunks and carries an error bound the crate can \
                     compute from its own output, which is a class of claim it \
                     has never made -- validate::accuracy measures distance and \
                     validate::isotopy certifies topology, and neither states a \
                     bound it derived rather than sampled. On a triangle mesh \
                     the Gaussian measure is the vertex angle defect and the \
                     mean measure is edge length times signed dihedral angle; \
                     Cohen-Steiner and Morvan's Theorem 6 bounds the deviation \
                     from the smooth surface by C*K*eps with K a sum of \
                     triangle circumradii and eps their maximum. Clause one: on \
                     sphere at 33, 65 and 129 samples per axis, the residual \
                     |sum of angle defects - 4*pi| falls inside the computed \
                     bound at all three. Clause two: that residual falls at \
                     least linearly as h halves. Clause three: on torus the \
                     defect sum is zero to within the same bound, which an \
                     accumulator that lost the sign would fail loudly. \
                     Recorded beside them, not registered: the mean-curvature \
                     total, and the chi recovered from the defect sum against \
                     the chi MeshReport computes independently.",
        falsified_by: "The residual exceeding its own computed bound at any of \
                       the three resolutions, or failing to fall as h halves. \
                       Either means Theorem 6's closely-inscribed hypothesis \
                       does not hold for a marching-cubes mesh -- whose \
                       vertices lie on the trilinear interpolant's zero set \
                       rather than on the field's -- and that is a better \
                       result than the bound holding, because it says the \
                       bound must be stated against the interpolant, which is \
                       the move isotopy.rs already makes.",
        records: &[
            "field",
            "samples_per_axis",
            "gaussian_total",
            "gaussian_expected",
            "residual",
            "bound",
            "within_bound",
            "mean_curvature_total",
            "chi_from_defect",
            "chi_from_report",
        ],
    },
    Preregistration {
        id: "P-43",
        ticket: "R-042",
        hypothesis: "One field evaluation at the cell centre, compared against \
                     the trilinear interpolant of the eight corners and \
                     normalised by cell size, is a usable witness that a chunk \
                     is under-sampled -- a statement the crate cannot make \
                     today, since validate::field_bound samples the gradient \
                     and is explicit that a sampled maximum is only a lower \
                     bound on a supremum. The witness is one-sided in the safe \
                     direction: it can prove a chunk inadequate and can never \
                     prove it adequate. Clause one: across noise_cavity and \
                     gyroid at 17, 33, 65 and 129 samples per axis, the \
                     per-grid maximum normalised centre residual correlates \
                     with the symmetric Hausdorff distance validate::accuracy \
                     already reports, at Pearson r of at least 0.7. Clause \
                     two: the extra evaluations are under 15% of the corner \
                     evaluations, which is the structural 1/8 plus slack.",
        falsified_by: "Pearson r below 0.7, meaning the centre residual does \
                       not witness the error it is supposed to predict. The \
                       failure mode is named in advance rather than \
                       rationalised afterwards: a feature that passes cleanly \
                       through the cell centre without perturbing it gives a \
                       residual near zero while the reconstruction is badly \
                       wrong, and thin_plate is constructed to be exactly \
                       that, so it is measured beside the two registered \
                       fields as the adversary.",
        records: &[
            "field",
            "samples_per_axis",
            "centre_residual_max",
            "centre_residual_mean",
            "symmetric_hausdorff",
            "pearson_r",
            "extra_eval_fraction",
        ],
    },
    Preregistration {
        id: "P-44",
        ticket: "R-042a",
        hypothesis: "P-43's mechanism survived its own falsification and its \
                     order statistic did not: the MAXIMUM normalised centre \
                     residual is pinned at every resolution by a single C1 \
                     crease cell, because a crease makes |f(centre) - \
                     mean(corners)| an O(h) quantity and dividing by h leaves a \
                     constant. The MEAN is not, and on P-43's own three fields \
                     it correlated at r = 0.983, 0.984 and 0.9998 with the \
                     symmetric Hausdorff, with log-log decay exponents agreeing \
                     to within 0.05. That is a hypothesis read off the data that \
                     killed its predecessor, so it is worthless until it is \
                     tested somewhere it can fail. This registers it on the four \
                     reference fields P-43 never touched -- sphere, torus, \
                     box_exact and csg_difference -- at 17, 33, 65 and 129 \
                     samples per axis. Clause one: the Pearson r between the \
                     per-grid MEAN normalised centre residual and the symmetric \
                     Hausdorff is at least 0.9 on each of the four. Clause two: \
                     the log-log decay exponent of the mean residual matches \
                     that of the symmetric Hausdorff to within 0.15 on each of \
                     the four. Clause three, the cost clause P-43 got wrong and \
                     which is restated with the right arithmetic: the witness \
                     costs one extra field evaluation per CELL against one per \
                     SAMPLE for the corners -- (n-1)^3/n^3, which rises toward 1 \
                     and is nowhere near an eighth, because this crate prefills \
                     one shared sample grid and every corner is evaluated once \
                     -- so the claim is on wall clock instead: computing the \
                     witness costs under 0.5x a Marching Cubes extraction of the \
                     same grid.",
        falsified_by: "r below 0.9 on any of the four fields, or an exponent gap \
                       above 0.15 on any of them -- either of which would mean \
                       the mean-residual agreement seen on P-43's three fields \
                       is a property of those fields' shared crease structure \
                       rather than of the witness, and the whole line dies with \
                       a second null rather than limping on. Clause three fails \
                       at 0.5x or worse, which would price the witness above the \
                       extraction it is diagnosing.",
        records: &[
            "field",
            "samples_per_axis",
            "centre_residual_mean",
            "symmetric_hausdorff",
            "pearson_r_mean",
            "decay_exponent_mean",
            "decay_exponent_hausdorff",
            "exponent_gap",
            "witness_ns",
            "extract_ns",
            "witness_cost_ratio",
        ],
    },
    Preregistration {
        id: "P-45",
        ticket: "R-041a",
        hypothesis: "P-42 took B to be the whole closed surface, which made its \
                     Gaussian clause an identity -- 3F = 2E, so the defect sum \
                     is 2*pi*chi combinatorially and the residual was one f64 \
                     epsilon per vertex with zero geometric content. The \
                     property this crate actually wants was therefore never \
                     tested: ADDITIVITY, N(A union B) = N(A) + N(B) - N(A \
                     intersect B), which is what lets a per-chunk curvature \
                     measure compose into a per-world one with no global pass. \
                     Partition an extracted closed mesh into a 4x4x4 grid of \
                     spatial chunks by triangle centroid, so every chunk is a \
                     patch with a real boundary. With the geodesic-curvature \
                     boundary term (pi minus the incident angle sum at a \
                     boundary vertex) and edges on a chunk boundary weighted one \
                     half: clause one, the sum over chunks of the Gaussian \
                     measure equals the global 2*pi*chi to within 1e-9 absolute \
                     on sphere and torus at 65 samples per axis; clause two, the \
                     sum over chunks of the mean measure equals the global sum \
                     of l*beta to within 1e-9 absolute on both; clause three, \
                     each chunk's value recomputed IN ISOLATION from its own \
                     triangles alone -- boundary detected from that set, nothing \
                     read from a neighbour -- reproduces its in-context value \
                     bit for bit, which is the difference between a measure that \
                     is chunk-local and one that merely sums.",
        falsified_by: "Any of the three exceeding its tolerance. A failure of \
                       clause one or two means the boundary term as transcribed \
                       is wrong; a failure of clause three means the measure \
                       needs a one-ring halo, which would make it composable in \
                       principle and not chunk-local in practice -- and this \
                       crate meshes chunks independently, so that is the \
                       difference between usable and not. Recorded beside them, \
                       not registered: half the mean measure against the \
                       analytic integral of H on box_exact, where the surface is \
                       grid-aligned and the estimator may be exact rather than \
                       convergent.",
        records: &[
            "field",
            "samples_per_axis",
            "chunks",
            "gaussian_global",
            "gaussian_chunk_sum",
            "gaussian_gap",
            "mean_global",
            "mean_chunk_sum",
            "mean_gap",
            "isolated_chunks_bit_identical",
        ],
    },
    Preregistration {
        id: "P-46",
        ticket: "R-040a",
        hypothesis: "P-41 measured a bijection -- critical cells == \
                     non-manifold vertices == critical cells hosting one, \
                     exactly, on every affected field -- which gives the repair \
                     a target it did not have when R-040 deliberately withheld \
                     it. Repair the SIGN LATTICE rather than the mesh, and do it \
                     by minimal value perturbation rather than by flipping a \
                     sign arbitrarily: for each critical cell, move the corner \
                     of smallest |value| across zero by the smallest \
                     representable step, which is Boutry's self-dual repair \
                     applied to the grey-level function rather than to the \
                     binary set. Clause one: on the three fields with a non-zero \
                     census -- noise_cavity, gyroid, fbm_terrain at 65 samples \
                     per axis -- DualContouring and SurfaceNets emit exactly 0 \
                     non-manifold edges and exactly 0 non-manifold vertices \
                     after the repair, where before they emitted 322/602, \
                     69/141 and 29/58. Clause two: the repair reaches a fixpoint \
                     in at most TWO sweeps over the chunk, because a \
                     perturbation can create a critical configuration in a \
                     neighbour and an unbounded cascade would rule the mechanism \
                     out of a frame budget outright; the residual critical count \
                     after two sweeps is recorded either way. Clause three, the \
                     price: the symmetric Hausdorff distance of the repaired \
                     mesh against the unrepaired field rises by less than 10% on \
                     each of the three fields -- the repair moves geometry, and a \
                     manifold mesh of the wrong surface is not a bargain.",
        falsified_by: "Any non-manifold output surviving the repair, which would \
                       mean well-composedness of the sign lattice is not \
                       sufficient for a manifold dual mesh even though P-41 \
                       measured it to be necessary -- and that gap is the more \
                       interesting result, because ManifoldDualContouring \
                       already pays extra vertices for the same guarantee. Or a \
                       residual critical count above zero after two sweeps, \
                       which prices the cascade out. Or Hausdorff rising 10% or \
                       more, which says the manifoldness was bought with \
                       geometry the caller did not agree to sell.",
        records: &[
            "field",
            "samples_per_axis",
            "extractor",
            "critical_before",
            "critical_after",
            "sweeps",
            "non_manifold_edges_before",
            "non_manifold_edges_after",
            "non_manifold_vertices_before",
            "non_manifold_vertices_after",
            "hausdorff_before",
            "hausdorff_after",
            "hausdorff_ratio",
        ],
    },
    Preregistration {
        id: "P-47",
        ticket: "R-043",
        hypothesis: "Every one of the eight reference fields overrides \
                     Sdf::gradient analytically, and the crate's docs say the \
                     central-difference default is never used by one. But \
                     BrushStack does not override gradient and neither does \
                     Capsule -- so the moment a caller composes anything, which \
                     is the entire point of a sculpting tool, every vertex \
                     normal silently falls back to the six-sample central \
                     difference at O(h^2), costing six evaluations of the WHOLE \
                     TAPE per normal. Carrying a forward-mode dual number \
                     (value, [dx, dy, dz]) through the fold fixes both halves in \
                     one traversal: min and max propagate the selected branch's \
                     derivative, smooth_min propagates the analytic derivative \
                     of its polynomial, and the tie at a CSG seam is broken \
                     deterministically by taking the lower index so the result \
                     stays a pure function of the inputs. Clause one, the \
                     accuracy hole, which the crate cannot currently state at \
                     all: on a 64-brush BrushStack the mean angular error of the \
                     central-difference normal against the exact dual-number \
                     normal exceeds 0.1 degrees and the maximum exceeds 5 \
                     degrees. Clause two, the speed: computing normals for a \
                     meshed chunk through the dual-number tape is at least 2x \
                     faster than through the six-sample central difference, \
                     since it is one traversal against six. Clause three, the \
                     control: on a bare Sphere, where the analytic gradient is \
                     already exact and known, the dual-number normal agrees with \
                     p/|p| to within 1e-12 -- an instrument that cannot \
                     reproduce a closed form has no business measuring a tape.",
        falsified_by: "A mean angular error at or under 0.1 degrees, which would \
                       mean the central difference is already good enough on \
                       composed fields and the hole is theoretical -- a real \
                       possibility, since DIFF_STEP is scaled by |p| and these \
                       shapes are smooth away from their seams. Or a speedup \
                       under 2x, which would mean the tape traversal is not the \
                       cost and the six evaluations were being amortised by \
                       something. Clause three failing means the dual arithmetic \
                       is wrong and nothing else in the row can be believed.",
        records: &[
            "fixture",
            "brushes",
            "vertices",
            "mean_angular_error_deg",
            "max_angular_error_deg",
            "central_ns_per_normal",
            "dual_ns_per_normal",
            "speedup",
            "sphere_control_max_error",
        ],
    },
    Preregistration {
        id: "P-48",
        ticket: "R-044",
        hypothesis: "validate/isotopy.rs names its own gap in its own header: \
                     'The general form needs interval arithmetic over an \
                     arbitrary F, which this crate has no way to do -- an Sdf \
                     hands back point values', so it settles for certifying the \
                     TRILINEAR INTERPOLANT and says plainly that it does not \
                     certify the analytic field against it. A compositional \
                     inclusion function over the crate's own field types closes \
                     that gap with no dependency and no interval library: (lo, \
                     hi) pairs and six operations, widened by one ULP per \
                     operation because core has no directed rounding, which \
                     keeps the enclosure sound in the only direction that \
                     matters. Clause one, SOUNDNESS, and it is fatal: over the \
                     eight reference fields at 33 samples per axis, every cell \
                     the certificate calls surface-free is surface-free under \
                     4096-point dense sampling -- zero unsound certifications, \
                     no tolerance. Clause two, NON-VACUITY: on the fields built \
                     from exact distance primitives -- sphere, box_exact, torus \
                     and csg_difference -- the certificate proves at least 90% \
                     of the cells that dense sampling shows to be surface-free. \
                     Clause three, REACH: the certified fraction is reported for \
                     all eight fields and is strictly greater than zero on at \
                     least six of them, so the mechanism is not confined to the \
                     eikonal cases.",
        falsified_by: "Any unsound certification whatever -- a cell called \
                       surface-free that dense sampling shows is not. That is \
                       fatal rather than a threshold miss: a certificate that \
                       can be wrong is not a certificate, and the failure would \
                       be in the one direction isotopy.rs says a predicate must \
                       never err. Or the certified share falling under 90% on \
                       the exact fields, which would mean the enclosure is too \
                       loose to be worth evaluating; or fewer than six fields \
                       above zero, which would confine it to the eikonal cases \
                       and leave gyroid, noise_cavity and fbm_terrain -- the \
                       three that actually go wrong -- outside its reach.",
        records: &[
            "field",
            "samples_per_axis",
            "cells",
            "cells_surface_free_sampled",
            "cells_certified_empty",
            "certified_fraction",
            "unsound_certifications",
            "undecided_fraction",
            "certified_vs_trilinear",
        ],
    },
    Preregistration {
        id: "P-49",
        ticket: "R-045",
        hypothesis: "The crate can say whether two places are connected and \
                     cannot say how big a thing fits between them. R-022a's Air \
                     tracker answers 'is this sealed'; a game asks 'can the \
                     player get through', and that is a bottleneck value rather \
                     than a boolean: for a pair of chunk faces, the maximum over \
                     air paths of the minimum distance-to-solid along the path. \
                     It is computable by a monotone union-find over air voxels \
                     processed in DESCENDING (field value, grid index), which is \
                     a total order and therefore deterministic with no PRNG, no \
                     atomics and no HashMap -- the same discipline the crate \
                     already applies to its weld. The output is a 6x6 symmetric \
                     matrix of face-to-face apertures plus a reachability mask, \
                     which is the composable boundary summary: neighbouring \
                     chunks combine their matrices with no global solve. Clause \
                     one, against exact ground truth: on a BoxExact slab with a \
                     Capsule of radius r subtracted along x, for r of 2, 4 and 8 \
                     cells, the reported -x/+x aperture equals r to within one \
                     cell size, and every other face pair reports unreachable. \
                     Clause two, on a real field: the capped gyroid at 65 \
                     samples per axis reports all six faces mutually reachable \
                     -- it is a bicontinuous triply periodic surface, so a \
                     disconnected pair is a proven bug in the instrument, not a \
                     property of the field -- with a positive aperture on all 15 \
                     pairs. Clause three, cost: the whole 6x6 computation costs \
                     under 2x a Marching Cubes extraction of the same grid.",
        falsified_by: "An aperture off by more than one cell on the drilled \
                       channel, where the answer is known exactly and the \
                       instrument has nowhere to hide -- that would mean the \
                       union-find is not computing the bottleneck it claims. Or \
                       any face pair reported reachable on the slab that is not, \
                       which is the unsound direction and is fatal for a \
                       clearance claim a game would gate movement on. Or a \
                       gyroid pair reported unreachable. Or cost at or above 2x, \
                       which prices a per-chunk summary above the mesh it \
                       summarises.",
        records: &[
            "fixture",
            "samples_per_axis",
            "channel_radius_cells",
            "aperture_reported_cells",
            "aperture_error_cells",
            "reachable_pairs",
            "expected_reachable_pairs",
            "false_reachable_pairs",
            "aperture_ns",
            "extract_ns",
            "cost_ratio",
        ],
    },
    Preregistration {
        id: "P-50",
        ticket: "R-039a",
        hypothesis: "P-40's C2 was a wall-clock ratio registered as a threshold, \
                     and it failed its own artefact: three quiet-machine runs \
                     read 1.022x, 1.184x and 1.177x against a registered 1.25x \
                     (M-348). Marching Cubes 24 already earned the rule -- gate \
                     the COUNT the ratio samples -- and the count exists. The \
                     bitmap prepass does not make an eight-corner gather faster; \
                     it REMOVES gathers, and how many run is an integer that is \
                     identical on every machine, under every governor, at every \
                     load. Clause one, as an EQUALITY and not a ratio: the \
                     scalar predicate performs exactly `cells` gathers and the \
                     bitmap predicate exactly `active_cells`, on every field at \
                     every resolution, with no tolerance. Clause two: the \
                     bitmap's own cost is exactly `sample_count` comparisons \
                     plus `cells_x.div_ceil(64) * cells_y * cells_z` fused word \
                     groups -- note CELLS and not SAMPLES, which is the defect \
                     E-307 found and which cost about 30% of the stage at 64k+1 \
                     grids -- so the prepass is O(n^3/64) word groups against the \
                     O(n^3) gathers it replaces, and the two counts predict the \
                     crossover without a clock. Clause three: the active set the \
                     two predicates name is the same ORDERED list, element for \
                     element, on every field and resolution, which is what makes \
                     the mesh identical rather than merely equal in count. \
                     Timing is recorded beside all three and gates nothing.",
        falsified_by: "Clause one failing as an exact equality. A gather \
                       performed on an inactive cell means the mask admits a \
                       cell it should not, and a mask wrong in the other \
                       direction drops a cell and punches a hole -- so this is a \
                       correctness gate wearing a performance gate's clothes, \
                       and a mismatch of one is a failure. Clause three failing \
                       would mean the set-bit walk does not reproduce \
                       lexicographic order, which changes every vertex index \
                       downstream of the first disagreement. Clause two failing \
                       would mean the word-group count is not what reading the \
                       loop says it is.",
        records: &[
            "field",
            "samples_per_axis",
            "cells",
            "active_cells",
            "gathers_scalar",
            "gathers_bitmap",
            "gathers_equal_cells",
            "gathers_equal_active",
            "bitmap_comparisons",
            "bitmap_word_groups",
            "word_groups_predicted",
            "same_ordered_list",
            "ns_per_cell_scalar",
            "ns_per_cell_bitmap",
        ],
    },
    Preregistration {
        id: "P-51",
        ticket: "R-046",
        hypothesis: "Sellan, Batty & Stein (10.1145/3610548.3618196) state the \
            full constraint a signed distance sample carries, and it has TWO \
            halves: the surface 'must be tangent to every sphere at least once \
            while strictly containing every sphere with negative value and \
            excluding every positive value one'. Every extractor in this crate \
            reads d as a number to interpolate and discards both halves. \
            Measured as integer counts over the output, at 65 samples per axis, \
            on the five fields declaring FieldBound::Exact -- sphere, torus, \
            box_exact, thin_plate, csg_difference. PIERCING, the exclusion half: \
            violation(v) = max over samples p within one cell of v of \
            (|d(p)| - ||v - p||) normalised by cell size, counting vertices \
            above 0.05 cells. (C1) marching_cubes pierces on fewer than 1 vertex \
            per 1,000, because an MC vertex is the root of the interpolant on a \
            grid edge and cannot lie deeply inside a neighbour's ball. \
            (C2) dual_contouring's rate is at least 20 TIMES marching_cubes' -- a \
            RATIO, so the clause cannot be cleared by a bar set below the ~150 \
            per 1,000 M-27 already measured (152 of 1,016 box_exact vertices \
            moving 0.35-0.57 cells). TOUCHING, the tangency half, which no \
            extractor here even attempts: (C3) the number of samples whose \
            sphere is never touched -- min over mesh vertices of \
            | ||v - p|| - |d(p)| | exceeding 0.05 cells -- is NON-ZERO for every \
            extractor on every field, and is strictly larger for marching_cubes \
            than for dual_contouring, because MC's vertices are confined to grid \
            edges and cannot reach a tangent point in a cell's interior. The \
            0.05-cell floor is derived, not chosen: M-12 measured h^2 \
            convergence and sphere at 65^3 has mean error 6.5e-4 against \
            h = 0.0635, i.e. 1.0% of a cell, so the gate sits five times above \
            the honest discretisation floor and far below M-27's 0.35 cells.",
        falsified_by: "C2's ratio under 20x, which says the exclusion half is \
            already respected by construction; or C1 exceeding C2's absolute \
            rate, which inverts the mechanism; or C3 finding ZERO untouched \
            spheres, which would say this crate's extractors already saturate \
            the tangency half and the four 2025-2026 papers built on recovering \
            it have nothing to offer here -- the most valuable of the three \
            outcomes to get wrong. The control that must pass first: on \
            box_exact, whose surface is planar and axis-aligned and where M-27 \
            measured every Dual Contouring vertex on a planar patch landing \
            exactly on the centroid, ALL THREE extractors must report zero \
            PIERCING. An instrument that finds piercing there is measuring its \
            own tolerance. Scope: the constraint needs the field to BE a \
            distance, and F-004/M-247 measured a voxel game's field degrading \
            from 0.577 to 0.004 over 256 strokes -- so gyroid, fbm_terrain and \
            noise_cavity are excluded by construction and a violation there \
            would not be a defect. Source honesty: the paper credits the \
            tangent-sphere interpretation to Batty 2011 and Kobbelt 2001 and \
            contributes the energy and flow; nothing here reuses its Table 1, \
            which measures a global remeshing gradient flow this crate does not \
            run.",
        records: &[
            "field",
            "extractor",
            "samples_per_axis",
            "vertices",
            "samples",
            "pierced",
            "pierced_per_1k",
            "worst_piercing_cells",
            "dc_over_mc_ratio",
            "untouched",
            "untouched_per_1k",
            "worst_untouched_cells",
            "samples_probed_per_vertex",
            "control_box_exact_zero",
            "threshold_cells",
        ],
    },
    Preregistration {
        id: "P-52",
        ticket: "R-047",
        hypothesis: "A third VertexRule that moves each dual vertex toward the \
            tangent points its cell's samples imply, using only Eq. (8) of \
            Sellan, Batty & Stein -- t_i = p_i + sigma_i |s_i| (c_i - p_i) / \
            ||c_i - p_i||, one normalize and one fma per sample, allocation-free \
            and Real-generic -- iterated TWICE and clamped to the cell, against \
            Qef through the identical classification and quad walk X-002 built. \
            THIS IS NOT THEIR ALGORITHM, which is a global sparse solve with \
            per-iteration remeshing over hundreds of iterations, and whose Fig. \
            17 ablation measures clamping away far spheres as progressive detail \
            loss; none of their reported accuracy is claimed here and the \
            baseline is this crate's own. (C1) On box_exact and thin_plate at \
            65^3 the symmetric Hausdorff improves by at least 1.25x over Qef. \
            (C2) On sphere and torus it is within +/-10% of Qef. (C3) Vertex, \
            triangle and non-manifold-edge counts are IDENTICAL to Qef on every \
            field, which is what says the rule changed placement and nothing \
            else -- M-237 established that property for the Qef/Centroid swap, \
            with byte-identical index buffers and all 680 positions different. \
            (C4) The mechanism M-315 predicts is visible: the vertex term of the \
            Hausdorff improves while the CENTROID term worsens, on at least 3 of \
            the 4 fields. M-315 measured that Dual Contouring's Hausdorff is \
            vertex-dominated on 8 of 8 rows AND that its centroid error is \
            already BETTER than the perfect-placement floor by 2.9-3.6x on \
            sphere, because the QEF minimises distance to tangent planes and \
            buys better-centred facets at the cost of worse-placed vertices. A \
            rule that pulls vertices onto spheres spends exactly that trade in \
            reverse.",
        falsified_by: "C1 under 1.25x, which says the tangency geometry does not \
            survive being clamped to its cell and iterated twice -- the honest \
            reading then being that this is offline CAD and not a chunk budget, \
            recorded as a null with a number. Or C3 failing, which means the arm \
            is not controlled and no number in it means anything. C4 is the \
            clause that can teach something whatever C1 does: if the vertex term \
            improves and the centroid term does NOT worsen, then M-315's \
            tangent-plane trade is not what the QEF is doing and a placement \
            rule is not zero-sum. The 1.25x bar is this crate's own: M-315 \
            measured the ceiling on ANY placement rule at 1.5-21.5% for SMOOTH \
            fields, which is why C1 is asked only of the two SHARP fields, where \
            no ceiling has been measured and where M-66 recorded box_exact's \
            worst normal error at 35.796 degrees identically at every \
            resolution -- a corner does not soften with h, so a placement rule \
            is the only thing that can move it.",
        records: &[
            "field",
            "samples_per_axis",
            "rule",
            "iterations",
            "symmetric_hausdorff",
            "hausdorff_ratio_vs_qef",
            "vertex_term",
            "vertex_term_ratio_vs_qef",
            "centroid_term",
            "centroid_term_ratio_vs_qef",
            "self_intersections_per_1k",
            "vertices",
            "triangles",
            "non_manifold_edges",
            "counts_identical_to_qef",
            "ns_per_sample",
        ],
    },
    Preregistration {
        id: "P-53",
        ticket: "R-048",
        hypothesis: "Custodio, Pesco & Silva (10.1186/s13173-019-0086-6) observe \
            that MC33 classifies a corner whose value EQUALS the isovalue as \
            inside, which marks all three incident edges cut and emits triangles \
            with coincident vertices. Their remedy is a third corner label. The \
            label assignment is a PURE PRE-PASS over the eight-corner sign \
            classification -- one strict comparison per corner instead of one \
            non-strict -- which is the half a bench can reproduce exactly; their \
            triangulator is a per-cube convex hull with cross-cell face dedup \
            and is NOT reproduced here and NOT claimed. On fuel (64^3) and \
            bonsai (256^3) at an INTEGER isovalue of 32: (C1) degenerate \
            triangles under marching_cubes are non-zero and at least 80% of them \
            are attributable to cells having a corner exactly equal to the \
            isovalue, measured by TAGGING each degenerate triangle with its \
            cell's equal-corner count, not inferred from a correlation. \
            (C2) with the ternary label and coincident-vertex collapse, \
            degenerate triangles fall by at least 10x while euler_characteristic, \
            non_manifold_edges and boundary_edges are UNCHANGED. (C3) at a \
            HALF-OFFSET isovalue of 32.5, where an integer sample cannot equal \
            the isosurface, the two paths produce BYTE-IDENTICAL meshes. \
            M-316 measured 16,284 of 529,508 bonsai surface-cell corners exactly \
            on the isosurface -- 3% -- and M-232 measured 20 singular faces per \
            400,000 cells at u8 density against 0 in continuous data, so this is \
            aimed at a defect already measured on this crate's real CT input.",
        falsified_by: "C1 under 80%, which would mean the degenerate count is \
            dominated by ordinary near-tangency slivers -- the thing CLAUDE.md \
            refuses to gate on -- and the paper solves a problem this crate does \
            not have. Or C2 changing the Euler characteristic, which would mean \
            the label is not topology-preserving as implemented here and rule 5 \
            applies. C3 is the control and is not decoration: M-317's own \
            guidance is to contour at a half-integer precisely because it is \
            unattainable by integer data, so if the two paths differ THERE, the \
            label is doing something beyond the exact-equality case and every \
            number in the row is suspect. Source honesty: the paper reports NO \
            count of degenerate triangles removed on any dataset -- only \
            radii-ratio histograms, Betti numbers and blocked-cube percentages -- \
            so C2's 10x is this crate's own bar and is not a reproduction of any \
            published figure.",
        records: &[
            "volume",
            "isovalue",
            "label_rule",
            "cells",
            "surface_cell_corners",
            "equal_corners",
            "degenerate_triangles",
            "degenerate_from_equal_corners",
            "degenerate_attributable_fraction",
            "degenerate_ratio",
            "triangles",
            "euler_characteristic",
            "non_manifold_edges",
            "boundary_edges",
            "mesh_hash",
            "half_offset_identical",
        ],
    },
    Preregistration {
        id: "P-54",
        ticket: "R-049",
        hypothesis: "M-248 measured empty-cell rejection by Hart's Lipschitz \
            bound at 16.8% of cells on gyroid against 80.6-95.1% on every other \
            field, and M-306 identified the cause: gyroid declares Lipschitz \
            l = 2*sqrt(3) = 3.464, derived correctly at M-244, while M-267 \
            measured its actual gradient supremum converging to 1.731. A revised \
            affine form (Fryazinov, Pasko & Comninos, 10.1016/j.cag.2010.07.003) \
            carries three noise symbols plus one accumulator -- FIVE stored \
            reals, fixed size, never growing, no heap -- and bounds an \
            expression over a BOX rather than a ball. The mechanism is \
            correlation: sin(x)cos(y) + sin(y)cos(z) + sin(z)cos(x) cannot have \
            all three terms extremal at once, and a per-term interval bound \
            throws that away. The prediction is therefore NOT uniform, and this \
            is the point: (C1) on gyroid at 17^3 the rejected-cell COUNT rises \
            from 688 of 4,096 by at least 1.5x, because gyroid is built only \
            from smooth trig with shared arguments and is exactly where \
            correlation lives. (C2) on box_exact and csg_difference the rejected \
            count rises by LESS THAN 5%, because both are built from min/max, \
            for which the source paper gives no affine rule at all and the only \
            sound treatment collapses the form to an interval and destroys every \
            correlation. (C3) the mesh is BYTE-IDENTICAL on every field, the one \
            property a rejection test must have, since a wrong rejection makes a \
            hole and a hole is invisible to every validity gate this crate has.",
        falsified_by: "C1 under 1.5x, which says the correlation slack is not \
            where the looseness lives and M-267's 2x gap is genuinely attainable \
            by the gradient rather than an artefact of the ball. C2 failing \
            UPWARD is the more interesting outcome and must be reported as such: \
            a tighter bound on a field whose constant is already tight would mean \
            the BALL geometry, not the constant, is what costs -- and that \
            generalises to every field. COUNTED, NOT TIMED: rejected cells are \
            integers, identical on every machine; evaluation cost is printed \
            beside them and gates nothing. Two things this registration owns \
            rather than cites: the source paper contains NO correlation argument \
            and NO quantified tightening figure -- only end-to-end wall-clock \
            tables -- so the 1.5x is derived here from M-267's measured \
            supremum, not reproduced; and the paper gives no min/max or abs rule, \
            so C2's mechanism is this crate's own reasoning about what the \
            absence of such a rule forces.",
        records: &[
            "field",
            "samples_per_axis",
            "bound",
            "cells",
            "rejected_cells",
            "rejected_fraction",
            "rejected_ratio_vs_lipschitz",
            "has_min_max",
            "mesh_identical",
            "mesh_hash",
            "bound_evals",
            "bound_ns_per_cell",
        ],
    },
    Preregistration {
        id: "P-55",
        ticket: "R-050",
        hypothesis: "validate checks manifoldness, orientation, Euler \
            characteristic, self-intersection, isotopy and Hausdorff accuracy. \
            NOTHING checks that the mesh's critical-point structure matches the \
            field's, and on gyroid and fbm_terrain the crate cannot even assert \
            chi, so those two fields have no topological gate beyond \
            manifoldness. Finken, Li, Wang, Guo & Levine (arXiv:2608.12142) \
            prove Theorem 1: a PL function monotonic with respect to a Morse f \
            has no spurious critical points. Their theorem is 2D and its \
            pigeonhole step -- 'since a triangle has only three edges' -- has no \
            hexahedral analogue, so what is tested here is a 3D PORT and is \
            labelled as such. Sampling k = max(2, ceil(||e||/w) + 1) points along \
            each mesh edge and calling it non-monotone when two directional \
            derivatives disagree in sign: (C1) marching_cubes on sphere, torus \
            and box_exact at 65^3 has ZERO non-monotone edges. (C2) gyroid and \
            fbm_terrain have a NON-ZERO count that FALLS monotonically across \
            17^3/33^3/65^3/129^3, making it a resolution witness rather than a \
            defect. (C3) noise_cavity has the highest per-1k count of the eight \
            fields, because it is the field with interior ambiguity (M-208).",
        falsified_by: "C1 non-zero, which would mean the gate is measuring the \
            sampling of the edge rather than the mesh and k is the problem; or \
            C2 flat or rising in resolution, which would make it a property of \
            the field rather than of the grid and useless as a witness. TWO \
            THINGS THIS CRATE OWNS, NOT THE PAPER: the paper gives a bare \
            sign-disagreement predicate with NO epsilon, NO relative tolerance \
            and NO flat-region guard, so the rule used here -- discard steps \
            under 1e-12 * (|f(a)| + |f(b)|), fixed in this registration before \
            the harness exists, and recorded at 1e-14 and 1e-10 beside it so the \
            answer's sensitivity is visible -- is isomesh's invention and must \
            not be attributed to Finken et al.; and the paper obtains gradients \
            by autodiff on a neural field while this uses the crate's central \
            difference, which changes the noise story in a way the paper never \
            analyses. The theorem is also inapplicable to the trilinear \
            interpolant, under which interior critical points genuinely exist -- \
            that is the origin of the ambiguous-face problem -- so a HELD C1 is \
            evidence about this crate's meshes, never a proof transported to 3D.",
        records: &[
            "field",
            "extractor",
            "samples_per_axis",
            "edges",
            "non_monotone_edges",
            "non_monotone_per_1k",
            "k_samples",
            "tolerance",
            "non_monotone_at_1e14",
            "non_monotone_at_1e10",
            "worst_reversal",
            "falls_with_resolution",
        ],
    },
    Preregistration {
        id: "P-56",
        ticket: "R-051",
        hypothesis: "P-47's accuracy clause died by three orders of magnitude -- \
            mean 7.6e-5 degrees against a registered 0.1 -- and its own artefact \
            says why: bulk_mean_angular_error_deg is 1.9e-8 while one vertex in \
            57,470 carries 4.365 degrees, vertices_over_1deg is 1, and \
            worst_stencil_straddles_seam is true from 32 brushes upward. The \
            surviving hypothesis is narrower and mechanical. At a min/max CSG \
            seam the field is C0 and not C1, so a central difference whose \
            six-sample stencil straddles the seam averages two different \
            gradients and the returned direction lies in the cone the two \
            branches span. The error is therefore at most HALF the angle between \
            them -- M-283's (180 - theta)/2 in a second setting -- and it does \
            NOT shrink with h, because the stencil step is DIFF_STEP * |p|, \
            independent of the grid. Over a swept family of two-sphere Subtract \
            fixtures with seam dihedral theta controlled from 30 to 175 degrees: \
            (C1) every vertex whose stencil straddles the seam has angular error \
            against the analytic gradient bounded by (180 - theta)/2, on every \
            fixture, with no exceptions. (C2) the count of such vertices scales \
            like the seam's length in cells, so it is O(n) on an n^3 grid rather \
            than O(n^2), which is why one vertex in 57,470 is the expected order \
            and not a fluke -- measured as a fitted exponent against n, required \
            to be under 1.5. (C3) vertices whose stencil does NOT straddle the \
            seam have mean error under 1e-6 degrees, the control that says the \
            effect is the seam and not the tape.",
        falsified_by: "C1 exceeding the bound on any fixture, which would mean \
            the error is not the two-branch average and M-283's mechanism does \
            not transfer; or C2's fitted exponent reaching 1.5 or above, which \
            would make it a surface-wide effect rather than a seam effect and \
            change what a consumer should do about it. The measured 4.365 \
            degrees predicts a seam dihedral near 171 degrees, so a fixture swept \
            down to 30 degrees must show the bound widening proportionally -- a \
            bound that holds only because it is loose everywhere is not evidence, \
            which is why tightness is recorded per fixture as \
            worst_over_bound_ratio and a median under 0.1 across the sweep would \
            be reported as a vacuous pass.",
        records: &[
            "dihedral_deg",
            "samples_per_axis",
            "seam_cells",
            "vertices",
            "straddling_vertices",
            "straddling_max_error_deg",
            "predicted_bound_deg",
            "worst_over_bound_ratio",
            "within_bound",
            "non_straddling_mean_error_deg",
            "scaling_exponent",
        ],
    },
    Preregistration {
        id: "P-57",
        ticket: "R-055",
        hypothesis: "Every element of the 48-element octahedral group is a \
            signed coordinate permutation and therefore exact in f64, so \
            mesh(g*f) and g*mesh(f) are comparable bit-for-bit. Compared as \
            SORTED VERTEX-POSITION MULTISETS, not as index buffers: table.rs \
            picks safe_apex by lowest edge index, which is not invariant under \
            axis relabelling, so a triangle-level relation is known-false in \
            advance. (C1) the four primal extractors (marching_cubes, \
            marching_cubes+decider, marching_tetrahedra, \
            subgrid_marching_tetrahedra) are bit-exactly equivariant on ALL 48 \
            elements, on all eight reference fields, because a primal vertex is \
            a/(a-b) on one grid edge -- two values, no accumulation. (C2) the \
            three dual extractors (surface_nets, dual_contouring, \
            manifold_dual_contouring) are bit-exact on STRICTLY FEWER than 48, \
            because M-177 established that a sum of position components is not \
            bit-exactly equivariant by ordering alone; the count and identity \
            of failing elements is the new number. (C3) a triangle-level \
            relation fails where the vertex-level one holds, on at least one \
            extractor, and the count of triangle-level mismatches is reported \
            so the dossier's \"2,688 false positives\" warning is quantified \
            rather than repeated.",
        falsified_by: "A primal failure -- which is an axis-dependent bug and \
            the thing this exists to catch; or all 48 holding on the dual path, \
            which would contradict M-177 and be the more interesting result; or \
            C3 finding zero triangle-level mismatches, which would mean \
            safe_apex is invariant after all.",
        records: &[
            "field",
            "extractor",
            "family",
            "samples_per_axis",
            "vertices",
            "elements_tested",
            "elements_vertex_exact",
            "elements_triangle_exact",
            "first_failing_element",
            "first_failing_det",
            "worst_component_ulp",
            "fixture_can_fail",
        ],
    },
    Preregistration {
        id: "P-58",
        ticket: "R-056",
        hypothesis: "Robins, Wood & Sheppard's ProcessLowerStars \
            (10.1109/tpami.2011.95) builds a discrete Morse function from a \
            sampled grid by pairing cells within each voxel's LOWER STAR -- at \
            most 27 cells in 3D -- which is per-voxel local and, by the paper's \
            own claim, produces a critical-cell census independent of the \
            processing order. The paper requires distinct values and perturbs \
            ties with a GLOBAL ramp that depends on the image dimensions I, J, \
            K; this crate's fields tie exactly, so this registration fixes a \
            CHUNK-LOCAL EXACT TIE-BREAK instead: order by (value, \
            linear_index) lexicographically, comparing values by \
            f64::total_cmp, which is a total order, deterministic, \
            allocation-free, and perturbs no sample. (C1) the census by \
            dimension (critical_0 .. critical_3) is IDENTICAL under the \
            registered tie-break and under its reverse (value, \
            Reverse(linear_index)), on all eight fields -- the paper's \
            ordering-independence claim tested on data it has never seen. (C2) \
            every cell the asymptotic decider flags ambiguous CONTAINS AT LEAST \
            ONE CRITICAL CELL -- stated as containment, not set equality, \
            because a Morse census can be non-empty where no MC ambiguity \
            exists and equality would indict the instrument; the excess \
            critical_cells_outside_ambiguous is reported per field either way. \
            (C3) the census is resolution-stable on simple topology and grows \
            on complex: critical_total on sphere and torus changes by less than \
            2x across 17^3/33^3/65^3 while on noise_cavity it grows by more \
            than 4x.",
        falsified_by: "C1 differing, which would mean the paper's \
            ordering-independence does not survive exact ties and the tie-break \
            is doing work it must not; or C2 finding an ambiguous cell with no \
            critical cell, which would break the containment the whole framing \
            rests on; or C3's sphere census growing like the grid, which would \
            make it an artefact count rather than a topological signature. \
            PROVENANCE CORRECTED BEFORE THE HARNESS EXISTED: the review that \
            proposed this registration asserted the corpus markdown terminates \
            mid-section-4 before Theorem 11 and that Theorem 11 is therefore \
            not transcribable. THAT IS FALSE, checked line by line -- section 4 \
            is complete (Lemmas 7-10), Theorem 11 is stated in full with both \
            of its arms, and so are Theorem 3, Theorem 6, Propositions 4-5, \
            Lemma 12 and the pseudocode of Algorithms 1, 2 and 3. The \
            ordering-independence C1 tests is therefore the paper's own \
            sentence, quoted: 'the results in Section 4 show that for 2D and 3D \
            complexes the number and type of critical cells found by \
            ProcessLowerStars are independent of this ordering'. What the paper \
            does NOT have is any test of that claim against exact ties, because \
            Eq. (8) assumes them away.",
        records: &[
            "field",
            "samples_per_axis",
            "voxels",
            "critical_0",
            "critical_1",
            "critical_2",
            "critical_3",
            "critical_total",
            "max_lower_star_cells",
            "census_matches_reverse_order",
            "ambiguous_cells",
            "ambiguous_with_critical",
            "ambiguous_containment_holds",
            "critical_cells_outside_ambiguous",
            "ns_per_voxel",
        ],
    },
    Preregistration {
        id: "P-59",
        ticket: "R-057",
        hypothesis: "M-341 measured that a Lipschitz interval bound prunes a \
            64-brush tape to a median of 19 survivors per chunk with the mesh \
            byte-identical on 64 of 64 chunks. Nothing measures how many \
            survivors are NECESSARY. Leave-one-out ablation decides it: for \
            each chunk, for each surviving brush, remove that one brush, \
            re-mesh, and compare mesh_hash. (C1) SOUNDNESS CONTROL, REPORTED \
            FIRST -- removing ALL non-survivors together changes no chunk's \
            mesh_hash, 64 of 64, re-confirming M-341's C3 on this harness; if \
            it fails, the bound is unsound and every other number here is void. \
            (C2) the bound OVER-KEEPS: the median over chunks of necessary / \
            survivors is AT MOST 0.75, i.e. at least a quarter of surviving \
            brushes can be dropped individually with the mesh bit-identical. \
            (C3) the over-keep has a nameable cause rather than being noise: of \
            the brushes that survive but prove unnecessary, at least 90% have \
            an interval over the chunk whose distance from zero exceeds one \
            cell size -- they win the min/max chain only where the field is far \
            from the surface.",
        falsified_by: "C2's median exceeding 0.75, which says the interval \
            bound is already close to tight and closes the direction with a \
            distribution rather than a hunch -- a null worth having; or C3 \
            under 90%, which means the over-keep is not explained by \
            distance-from-surface and the mechanism is unnamed.",
        records: &[
            "chunk",
            "brushes",
            "survivors",
            "necessary",
            "necessary_fraction",
            "non_survivors_removed",
            "control_hash_unchanged",
            "unnecessary_far_from_surface",
            "unnecessary_far_fraction",
            "mesh_hash",
            "remeshes",
            "ns_per_remesh",
        ],
    },
    Preregistration {
        id: "P-60",
        ticket: "R-058",
        hypothesis: "Every Marching Cubes vertex is a linear interpolation \
            between two corner samples. Blu, Thevenaz & Unser \
            (10.1109/tip.2004.826093) show that shifting the sampling knots by \
            a fixed, signal-independent tau_opt = (1 - sqrt(3)/3)/2 ~= 0.21 and \
            enforcing the interpolation property recovers \"about 8 dB \
            asymptotically\" over standard linear reconstruction, for \
            w < 3*pi/4. SCOPED TO A SINGLE GRID LINE, TOUCHING NO EXTRACTOR, \
            because the prefilter is a causal one-pole IIR and the \
            reconstruction change would re-derive the whole A-002 apparatus. \
            Along an axis-aligned line through each reference field, sampled at \
            65 points, with the exact root known from the field's analytic \
            distance: (C1) on the four smooth fields (sphere, torus, gyroid, \
            fbm_terrain) the shifted reconstruction's median root error is at \
            least 30% lower than standard linear at matched sample count, using \
            tau = 1/5, at which the recursion c_n = -2^-2 * c_{n-1} + (1 + \
            2^-2) * f_n is multiplication-free. (C2) THE PRE-REGISTERED \
            FAILURE -- on the two fields with a step-like restriction \
            (box_exact, csg_difference) the shifted method is WORSE OR EQUAL, \
            because the paper states a Gibbs phenomenon on a step and a sharp \
            CSG boundary is one; this is expected, and finding it absent would \
            be the surprise. (C3) the prefilter's non-locality is bounded: \
            computing it over a truncated window of k preceding samples instead \
            of the whole line changes the recovered root by less than 1e-6 \
            cells for k >= 10, consistent with the (tau/(1-tau))^k = (1/4)^k \
            decay.",
        falsified_by: "C1 under 30%, which says the reconstruction gain does \
            not transfer to root position and the direction closes for a \
            mesher; or C3 failing at k = 10, which would mean a chunked mesher \
            cannot have this at any affordable guard band. NO GOLDEN HASH MOVES \
            -- this experiment reads a 1-D sample line and never calls an \
            extractor.",
        records: &[
            "field",
            "samples",
            "tau",
            "root_error_standard",
            "root_error_shifted",
            "error_ratio",
            "median_error_ratio",
            "is_step_like",
            "gibbs_overshoot",
            "guard_band_k",
            "guard_band_delta_cells",
            "guard_band_converged",
        ],
    },
    Preregistration {
        id: "P-61",
        ticket: "R-059",
        hypothesis: "M-356 found a mesh bit-exactly equivariant under axis \
            RELABELLING and not under REFLECTION, and attributed it to a/(a-b) \
            and b/(b-a) being two divisions of the same two values. That is \
            right and incomplete. The subtraction is innocent: IEEE \
            round-to-nearest is sign-symmetric, so fl(b-a) = -fl(a-b) EXACTLY. \
            What breaks is the ANCHOR -- cube::edge_crossing returns a \
            parameter measured from the LOWER corner and the placements put the \
            vertex at lower + t, so a reflection swaps which corner is lower \
            and the correct reflected parameter is 1 - t, which is not \
            representable as b/(b-a). Store the crossing instead as a signed \
            offset from the edge MIDPOINT: d = ((a + b) / 2) / (a - b) in \
            [-1/2, +1/2], position = edge_midpoint + h * d. This is exactly \
            antisymmetric under the simultaneous endpoint-and-sign swap by four \
            IEEE 754 guarantees rather than by observation: fl(a+b) = fl(b+a) \
            because addition is commutative; halving is exact because 2 is a \
            power of two; fl(b-a) = -fl(a-b) because round-to-nearest is odd; \
            and fl(S / -D) = -fl(S / D) for the same reason. The [0,1] \
            parameter frame cannot host this because reflection acts there as \
            0 <-> 1, an AFFINE map, and floating point respects sign flips \
            exactly and affine maps only approximately. THIS IS A src/ CHANGE, \
            registered as one: cube::edge_crossing becomes cube::edge_offset \
            and all five placements move to the centred frame -- \
            marching_cubes/mod.rs edge_position (which computes a/(a-b) \
            inline), marching_cubes/trilinear.rs local_crossing, \
            surface_nets.rs, hermite.rs and transvoxel/cell.rs. (C1) The \
            octahedral relation is re-measured with P-57's own fixtures and \
            group, and elements_vertex_exact reads 48 on EVERY row that \
            fixture_can_fail marks true, in both the primal and the dual \
            family, at both 33 and 25 samples per axis. POPULATION, counted \
            from docs/experiments/p-57.csv before this harness was written: 98 \
            of 112 rows have fixture_can_fail = true -- 56 primal and 42 dual \
            -- and ZERO of those 98 currently reach 48, so the clause has 98 \
            rows on which it can fire and all 98 currently fail. The 14 \
            box_exact rows are excluded by fixture_can_fail because that \
            field's zero set lies on dyadic planes and its 8 rows at 48 of 48 \
            are the fixture rather than the extractor. The harness asserts the \
            fixture columns cut_edges, order_sensitive_edges, grid_symmetric \
            and fixture_can_fail against p-57.csv row for row before reporting \
            anything, so a drift in the re-implemented group or fixture is a \
            failure rather than a new number. (C2) A COST CLAUSE, NOT A BENEFIT \
            CLAUSE -- no reference field's golden hash is unchanged. The vertex \
            positions genuinely move, T-007's 216 golden hashes are rebaselined \
            in the same commit as the source change, and a claim that they do \
            not move would mean the centred form is not on the path that \
            produces vertices. Measured per fixture as edges_moved: the number \
            of cut grid edges whose centred world position differs in bits from \
            the lower-corner one, over all eight fields at 33 and 65 samples \
            per axis. POPULATION: every field has cut edges at both \
            resolutions, p-57.csv reading 1350 for box_exact alone at 33, so \
            the clause fires on 16 of 16 rows. (C3) Symmetric Hausdorff changes \
            by less than 1 percent in either direction on all eight reference \
            fields at 33 and 65 samples per axis, and the \
            self-intersections-per-1000 count changes by less than 1 percent. \
            SCOPED TO marching_cubes, which the doc's clause leaves unstated: \
            the sign classification is untouched by a placement change, so the \
            topology and the index buffer are bit-identical and the two arms \
            are the same mesh with substituted positions, which is an exact \
            comparison. A dual vertex is a QEF solve over the crossings and its \
            lower-corner counterpart is not recoverable by substitution, so it \
            is not claimed. (C4) REGISTERED AS MEASUREMENT WITH A STATED RISK \
            OF A NULL. M-32 says chunk seams are bit-exact only when the cell \
            size is a power of two, and the mechanism it names is \
            world_of_sample rather than the crossing, so the honest prediction \
            is that C4 changes nothing. But the off-repo pre-measurement moves \
            from 75.0 percent mismatches to 0 at h = 0.1, so the seam sweep is \
            re-run at h = 0.1 and h = 3/32 as well as h = 0.125 and the answer \
            recorded either way. The appendix pre-measurement is re-run inside \
            this harness in Real for BOTH scalars rather than quoted: 2000000 \
            random straddling pairs per row, cell-local and world frames, both \
            forms.",
        falsified_by: "C1: any row below 48 that fixture_can_fail marks true. \
            Two mechanisms are already known that this change cannot reach and \
            a falsification on them is expected rather than surprising -- the \
            two tetrahedral extractors fail because a six-tetrahedron \
            decomposition of a cell is not octahedrally invariant, which M-356 \
            measured at 6 to 12 of 48 even on box_exact where \
            order_sensitive_edges is 0, and the duals accumulate into A^T A in \
            an order that axis relabelling permutes. C2: any field whose golden \
            hash survives, or any fixture with edges_moved = 0, either of which \
            would mean the centred form is not on the path that produces \
            vertices. C3: any field moving more than 1 percent on either \
            metric -- in which case the roughly 3-ulp worst-case tail the \
            off-repo measurement recorded is buying a real geometric cost and \
            the trade is a decision rather than a fix, to be raised and not \
            merged. C4: NOTHING. That clause is registered as measurement and a \
            null is the expected outcome; it is reported either way. And the \
            pre-measurement itself is falsified by any non-zero mismatch count \
            for the centred form in either scalar, which would mean the \
            three-line IEEE argument is wrong.",
        records: &[
            "block",
            "field",
            "extractor",
            "samples_per_axis",
            "scalar",
            "cell_size",
            "elements_vertex_exact",
            "elements_vertex_exact_p57",
            "elements_triangle_exact",
            "pure_permutation_exact",
            "pure_sign_flip_exact",
            "fixture_can_fail",
            "cut_edges",
            "edges_moved",
            "worst_move_ulp",
            "hausdorff_lower_corner",
            "hausdorff_centred",
            "hausdorff_ratio",
            "self_intersections_per_1k_lower_corner",
            "self_intersections_per_1k_centred",
            "self_intersections_ratio",
            "seam_vertices",
            "seam_mismatches_lower_corner",
            "seam_mismatches_centred",
            "pairs",
            "mismatches_lower_corner",
            "mismatches_centred",
            "out_of_cell_centred",
            "c1_population",
            "c1_rows_at_48",
            "c1_holds",
            "c2_fixtures_with_moved_edges",
            "c2_holds",
            "c3_worst_hausdorff_ratio",
            "c3_worst_self_intersection_ratio",
            "c3_holds",
            "c4_seam_mismatch_delta",
            "p57_fixture_columns_match",
        ],
    },
    Preregistration {
        id: "P-62",
        ticket: "R-060",
        hypothesis: "P-48 gave this crate a certificate that a cell is EMPTY \
            (M-347: zero unsound over 1.07e9 evaluations) and P-54 a tighter \
            one via affine arithmetic (M-354: 3.85x more rejections on \
            gyroid). There is no certificate in the other direction: nothing \
            can say 'this cell's surface patch has no hidden topology', which \
            is the difference between a mesher that is correct and one that can \
            STATE WHERE it is correct, and it is the CAD half of the mandate. \
            THE PREDICATE IS ALREADY IN THE TREE, so this is a measurement \
            rather than a build: validate::isotopy::cell_is_certified shipped \
            under T-015 and is exactly the registered form, \
            0 not-in box-F(C) OR <box-grad-F(C), box-grad-F(C)> > 0, from \
            Plantinga & Vegter, Isotopic approximation of implicit curves and \
            surfaces, SGP 2004, 10.1145/1057432.1057465. Both bounds are EXACT \
            rather than interval approximations, because the surface Marching \
            Cubes approximates is the trilinear interpolant: F is a convex \
            combination of the eight corner values, so clause one is exactly \
            'all eight corners share a sign', and each partial derivative is \
            bilinear in the other two coordinates, so its exact range is the \
            min and max of four corner differences. The cell size cancels -- \
            the predicate tests the SIGN of a sum of three squares, so h^2 \
            factors out. WHAT HAS NEVER BEEN DONE IS THE SOUNDNESS CHECK, and \
            this crate owns a ground truth the PV literature does not: A-020's \
            classifier counts tunnels and twelve-vertex contours from the \
            trilinear itself, M-214 recorded 2,053 and 173 in 396,000 cells, \
            and M-222 established that chi falls by exactly two per tunnel. A \
            cell containing a tunnel is a cell whose patch is NOT a graph, so a \
            certificate on such a cell is unsound. (C1) SOUNDNESS, ONE-SIDED, \
            THE KILL-SHOT: over eight reference fields at 17^3 / 33^3 / 65^3, \
            ZERO cells are C1-certified while the A-020 classifier reports a \
            tunnel or a twelve-vertex contour in them. (C2) YIELD: the \
            certified fraction of SURFACE cells is above 50% on sphere, torus \
            and box_exact at 33^3, and rises monotonically with resolution on \
            all eight fields. (C3) COST: the predicate costs under 5% of \
            extraction wall time on marching_cubes at 65^3. REGISTERED CAVEAT, \
            not a discovery: C1 guarantees the patch is a GRAPH, not that its \
            planar domain is connected -- a graph over a disconnected planar \
            region still has several components. PV close that globally with a \
            BALANCED octree this crate does not have, so the honest claim is \
            'no hidden topology in this cell' and not 'exactly one component', \
            and the entry must say so. Lin & Yap document the same gap, \
            10.1007/s00454-011-9345-9.",
        falsified_by: "C1 by ONE certified cell that the A-020 classifier calls \
            a tunnel or a twelve-vertex contour. A single unsound certificate \
            kills the direction, and M-214's counts of 2,053 and 173 prove the \
            fixture can produce the configuration -- this is not an M-44 pass \
            over an unreached case, and the population of tunnel cells is \
            reported per field so a zero is shown to have been able to be \
            non-zero. C2 by a certified fraction below 50% on any of sphere, \
            torus or box_exact at 33^3, or by a non-monotone sequence on any \
            field -- the latter would mean the predicate is measuring the \
            arithmetic's slack rather than the field's geometry. C3 by above \
            5%, in which case it is a debug gate rather than a shippable \
            capability. NO GOLDEN HASH MOVES: the predicate is read-only and \
            no extractor calls it.",
        records: &[
            "field",
            "samples_per_axis",
            "cells",
            "surface_cells",
            "certified_cells",
            "certified_surface_cells",
            "certified_surface_fraction",
            "tunnel_cells",
            "twelve_vertex_cells",
            "unsound_certificates",
            "monotone_in_resolution",
            "predicate_ms",
            "extract_ms",
            "predicate_share",
            "c1_holds",
            "c2_holds",
            "c3_holds",
        ],
    },
    Preregistration {
        id: "P-63",
        ticket: "R-061",
        hypothesis: "O-12 is the oldest open question in the ledger -- 'is \
            Marching Cubes unconditionally manifold now?' -- and its own text \
            says what would settle it: an exhaustive search over \
            configurations spanning more than two cells, or a proof that a \
            cell-local cycle triangulation plus shared face segments cannot \
            produce a non-manifold VERTEX. The search space is much smaller \
            than it looks. In this crate every Marching Cubes vertex sits on a \
            grid EDGE, so every face incident to an edge vertex comes from one \
            of the FOUR cells sharing that grid edge. Those four cells span a \
            3 x 3 x 2 block of grid corners -- 18 corners, 2^18 = 262144 sign \
            patterns. That is not a sample, it is the whole space, and it runs \
            in seconds. This is a proof by exhaustion of the vertex-link case \
            for Marching Cubes, and it is the case Chernikov & Xu's Coq work \
            does not cover: their 2013 IMR proof enumerates all 2^8 \
            single-cube configurations and proves cohesion and water-tightness, \
            then composes to a grid via FACE-local consistency, which is \
            exactly the argument that does not reach a vertex link. (C1) Over \
            all 262144 patterns, meshing the four cells, welding, and walking \
            the connected components of the incident-face link of the shared \
            edge vertex yields ZERO non-manifold vertices, with the \
            interior-ambiguity rule both off and on. POPULATION, derived before \
            the harness: the shared edge is the block's central z-edge and is \
            cut exactly when its two endpoint corners differ in sign, which is \
            HALF the patterns -- 131072 per interior rule, 262144 link walks in \
            total, and the clause can fire on every one of them. The \
            shared-edge vertex is ALSO the only vertex in the block whose link \
            is COMPLETE: every other edge vertex has cells missing outside the \
            block, so a defect there could be an artefact of the truncation and \
            is reported separately rather than counted. A cell's INTERIOR \
            vertex is complete too, because all of its faces come from that one \
            cell, so it is counted. (C2) THE FIXTURE CAN FAIL, and the defect \
            is the one known to exist: injecting the pre-fix single-apex fan \
            into the same sweep produces a NON-ZERO count. Reproduced \
            bench-locally as a vertex identification rather than by editing \
            src/ -- merging all of one cell's interior apexes into a single \
            vertex IS the pre-fix topology, same triangles and one shared apex, \
            which is precisely what the per-ring apex fix undid. POPULATION: \
            patterns in which some cell fans two or more rings, measured and \
            reported as fan_patterns; a zero there voids C2 rather than \
            passing it. (C3) THE DUAL FAMILY IS WHERE IT IS INTERESTING. The \
            same sweep over surface_nets, dual_contouring and \
            manifold_dual_contouring produces a non-zero count, and that count \
            is a FUNCTION of the well-composedness census in the sense of \
            M-338's bijection: the number of link-defective dual vertices \
            equals the number of critical sign configurations in the block. \
            POPULATION: 128 of the 256 possible cell sign bytes are critical -- \
            120 by a checkerboard 2x2 face and 8 by a main-diagonal inside pair \
            or its complement -- so a critical cell is reachable on the great \
            majority of the 262144 patterns and the clause cannot be vacuous. \
            SCOPE NOTE, WHICH THE ENTRY MUST NOT EXCEED: C1 and C2 are \
            COMPLETE for Marching Cubes. For the dual family a vertex lives at \
            a cell CENTRE and its link involves the cell's 26 neighbours, which \
            is 4^3 = 64 corners and out of reach of this block, so C3 is a \
            NECESSARY-CONDITION sweep on the same 18 corners and nothing more. \
            The full dual sweep at 2^27 over a 3 x 3 x 3 corner block -- \
            134217728 patterns, an estimated 4 to 45 minutes single-threaded -- \
            is a nightly gate and a separate ticket, deliberately not \
            registered here.",
        falsified_by: "C1: one pattern whose shared-edge vertex link has more \
            than one component, which is the THIRD MECHANISM O-12 asks about \
            and takes the next free falsified id. (The research doc's own text \
            names '49' for this outcome, written before Phase 23 ran; that id \
            went to P-61's falsification, and ids are assigned when used and \
            never reused.) C2: ZERO -- which would mean the link walk cannot \
            see the one defect known to exist and the sweep proves nothing, \
            voiding C1 rather than confirming it. C3: a count that is non-zero \
            and NOT equal to the critical-configuration census, which is the \
            more interesting outcome and would mean M-338's bijection is \
            cell-local and does not extend to a block. NO GOLDEN HASH MOVES: \
            this experiment builds its own 3 x 3 x 2 sign lattice and never \
            touches a reference field or a src/ path.",
        records: &[
            "arm",
            "extractor",
            "interior_rule",
            "patterns",
            "patterns_shared_edge_cut",
            "shared_edge_vertices",
            "link_defective_shared_edge",
            "link_defective_interior",
            "link_defective_truncated",
            "worst_link_components",
            "first_defective_pattern",
            "fan_patterns",
            "critical_cells",
            "critical_patterns",
            "defective_equals_critical",
            "c1_holds",
            "c2_holds",
            "c3_holds",
            "wall_ms",
        ],
    },
    Preregistration {
        id: "P-64",
        ticket: "R-062",
        hypothesis: "Bounded model checking proves the combinatorics that \
            property tests only sample. THE SPLIT THAT MAKES THIS TRACTABLE: \
            bit-blasting IEEE 754 to SAT is the adversarial case for a model \
            checker, and a harness over eight nondeterministic f32 corner \
            values is 256 bits of unconstrained float -- precisely the shape \
            Kani is worst at. But this crate's correctness risk is not in the \
            arithmetic; CLAUDE.md rule 5 names it exactly: 'wrong case tables \
            produce meshes that look fine and are subtly non-manifold.' That is \
            COMBINATORICS OVER EIGHT SIGN BITS, which is 256 states and trivial \
            for BMC. So: verify the combinatorics, keep testing the \
            arithmetic. Kani, arXiv:2607.01504, ranges over the eight sign bits \
            and proves, for all 256 patterns, that no case-table index goes out \
            of range, no emitted index is at or past the vertex count, no \
            triangle carries two equal indices, and nothing panics. Both are \
            dev tools with no runtime footprint, so hard rule 3 is not engaged, \
            and NO PUBLISHED USE OF EITHER ON GEOMETRY OR GRAPHICS CODE WAS \
            FOUND, which makes this novel as well as useful. (C1) Kani proves \
            all four properties over all 256 sign patterns for marching_cubes \
            with the interior rule off, in under 10 minutes. (C2) It finds \
            nothing the existing suite does not already cover. This is \
            registered as the EXPECTED outcome and is still worth the run: a \
            proof and a passing property test are different objects, and \
            M-208 to M-213 is five pre-registered claims that were true on \
            seven fields and false on the eighth. (C3) Turning the \
            interior-ambiguity rule on keeps C1 under 30 minutes. SCOPE NOTE: \
            neither tool touches vertex placement, and the registration must \
            say so -- placement stays under proptest and golden hashes. The \
            honest scope is 'the table cannot be indexed wrongly', not 'the \
            mesh is correct'.",
        falsified_by: "C1 by a timeout, or by a property that cannot be \
            expressed against the sign abstraction -- the second outcome is the \
            more informative one and means the abstraction is wrong, not the \
            tool. C2 by Kani finding a reachable violation, which is a \
            falsification entry and the most valuable outcome available here. \
            C3 by a blow-up, which localises the state explosion to the \
            interior rule and is itself a finding about that rule's branching. \
            NO GOLDEN HASH MOVES: the proofs are read-only and no extractor \
            calls them. VOID if the solver reports zero checks for a property, \
            since a proof over an empty check set is M-44's vacuous zero in \
            formal clothing -- each harness must report its check count.",
        records: &[
            "harness",
            "property",
            "interior_rule",
            "patterns",
            "checks",
            "failed_checks",
            "status",
            "solver_seconds",
            "kani_version",
            "c1_holds",
            "c2_holds",
            "c3_holds",
        ],
    },
    Preregistration {
        id: "P-66",
        ticket: "R-064",
        hypothesis: "THE LINE THIS REPLACES DIED TWICE, and the third attempt \
            has to be a different mechanism. P-43 / x29 tried one evaluation at \
            the cell centre as an under-sampling witness and was falsified on \
            both clauses; P-44 / x31 tried the mean residual instead and was \
            falsified out of sample. Both were VALUE witnesses: they asked \
            whether the trilinear's value at a point disagrees with the \
            field's. The failure mode they were chasing is not a value \
            disagreement, it is a MISSED ROOT -- an edge whose two endpoints \
            have the same sign while the field crosses zero twice between them, \
            or opposite signs while it crosses three times. THE WITNESS IS A \
            DERIVATIVE SIGN TEST, from Finken, Li, Wang, Guo & Levine, \
            Topology-Preserving Meshing of Implicit Scalar Fields via \
            Monotonicity Constraints, arXiv:2608.12142, IEEE Vis 2026 short \
            paper. Their central statement: if every edge of a PL mesh is \
            monotonic with respect to f, the PL approximation is topologically \
            consistent with f's critical points. The test itself is one line -- \
            sample the directional derivative grad-f dot e-hat at k points \
            along the edge and declare it non-monotonic when any two sampled \
            projections disagree in sign. DO NOT PORT THE PAPER: it is \
            explicitly 2D and the authors say so; its sampling-density \
            argument, its Theorem 1 case analysis (3D Morse theory has four \
            critical-point types and a spherical link, not three and a circle) \
            and its separatrix refinement all fail to generalise, and it wants \
            a Hessian this crate's field trait does not expose. Take the EDGE \
            TEST ALONE, as a diagnostic rather than an extraction rule. WHAT \
            MAKES THIS REGISTRABLE HERE RATHER THAN ANYWHERE ELSE IS THE \
            ORACLE: subgrid::roots::all_roots already finds ALL roots along an \
            edge -- M-94 resolved a slab at 1/1000 of the edge length, M-168 \
            gave each crossing an identity, and M-169 established that \
            identity-based sharing is complete exactly when no root lands on a \
            grid sample point -- so the number of roots per edge is a KNOWN \
            QUANTITY IN THIS REPOSITORY and the monotone-edge test can be \
            scored against it rather than against a hunch. (C1) SOUNDNESS, \
            ONE-SIDED: on eight reference fields at 17^3 / 33^3 / 65^3, with \
            k = 5 samples per edge, every edge the subgrid root finder reports \
            with more than one root is flagged non-monotonic. ZERO FALSE \
            NEGATIVES. (C2) YIELD: the false-positive rate -- edges flagged \
            non-monotonic that carry exactly one root -- is below 20% at k = 5 \
            on the six smooth fields, and the rate falls as k rises. (C3) IT IS \
            A RESOLUTION WITNESS, WHICH IS THE POINT: the non-monotonic edge \
            fraction falls monotonically with resolution on all eight fields, \
            and is highest on thin_plate -- the field whose sub-cell features \
            Marching Cubes structurally cannot see (M-100). CONSEQUENCE IF IT \
            HOLDS: a chunk can report 'this grid under-resolves this field, \
            here' as a number, per chunk, cheaply -- which is the missing input \
            to an LOD decision that M-121's 3.14-cell surface pop and M-72's \
            aliasing both want and neither has.",
        falsified_by: "C1 by ONE multi-root edge the test calls monotonic. C2 \
            by a false-positive rate above 20% on any of the six smooth fields, \
            or a rate that does not fall with k, which would mean the test is \
            measuring sampling noise rather than the field. C3 by a \
            non-monotone sequence in resolution on any field, or by thin_plate \
            not ranking first. NO GOLDEN HASH MOVES: the test is a diagnostic \
            and no extractor calls it. VOID if the oracle reports no multi-root \
            edge at all, since C1's zero false negatives would then be a zero \
            over an empty population -- the multi-root edge count is a recorded \
            column and must be non-zero somewhere in the sweep. The oracle's \
            own resolution is a stated limit rather than an assumption: it \
            divides each edge into a fixed number of intervals and cannot see a \
            root pair closer together than one of them, so the oracle sample \
            count is recorded per row.",
        records: &[
            "field",
            "samples_per_axis",
            "k",
            "oracle_samples",
            "edges",
            "single_root_edges",
            "multi_root_edges",
            "flagged_non_monotonic",
            "false_negatives",
            "false_positives",
            "false_positive_rate",
            "non_monotonic_fraction",
            "falls_with_k",
            "falls_with_resolution",
            "thin_plate_ranks_first",
            "c1_holds",
            "c2_holds",
            "c3_holds",
        ],
    },
    Preregistration {
        id: "P-69",
        ticket: "R-067",
        hypothesis: "core::simd is nightly and is staying nightly -- the \
            LaneCount<N>: SupportedLaneCount bound and the mask-element-type \
            mismatch are unresolved (rust-lang/portable-simd#364) and the \
            maintainers' own 2025 summary is 'nightly-only and will remain such \
            for the foreseeable future'. So the lever is AUTOVECTORISATION, and \
            the measured prior says autovectorisation is enough: Wilcox's \
            AArch64/NEON study on 100k f32 samples measured scalar 77.67 us, \
            hand-written intrinsics 25.78 us, and autovectorised safe Rust \
            25.54 us -- safe code matched intrinsics. The patterns that decide \
            it are shape, not machinery: struct-of-fields rather than index \
            arithmetic, pre-slicing once outside the loop so LLVM can prove the \
            bound, and chunks_exact / zip iterators. dual.rs's sample() pushes \
            into a Vec inside a triple loop, so the bound is re-proved per \
            element and the store is not a provable contiguous write. THE FLOAT \
            CAVEAT CUTS THE RIGHT WAY: the blanket claim that autovectorisation \
            fails on floats is about REDUCTIONS, because LLVM will not \
            reassociate float adds without fast-math and stable Rust does not \
            expose it. Elementwise float map and zip vectorise fine. This \
            crate's field evaluation is elementwise over independent samples and \
            its accumulations -- active-cell popcounts, vertex counts -- are \
            INTEGER. THIS IS A src/ CHANGE, registered as one. (C1) \
            Restructuring the sample loop to a pre-sliced contiguous write with \
            the bound hoisted gives at least 2x on the marginal f32 cost. \
            POPULATION AND TWO CORRECTIONS TO THE CLAUSE'S OWN TERMS, both \
            established from source before this harness was written. FIRST, the \
            fields: the doc names sphere and gyroid. libm 0.2.16's sqrtf carries \
            a select_implementation on target_feature = 'sse2' (and on \
            aarch64+neon), so on x86-64 it is a hardware instruction and \
            sphere's body -- sqrt(x^2+y^2+z^2) - r -- CAN vectorise. libm's sinf \
            and cosf carry NO arch selection at all: they are pure software with \
            argument-reduction branches, so gyroid's body, which is six of them \
            per sample, CANNOT vectorise at any loop shape while libm is the \
            transcendental path (CLAUDE.md rule 4). So C1 can fire on sphere and \
            is expected to fail on gyroid for a structural reason named here \
            rather than discovered afterwards, and the clause is registered as \
            the doc states it rather than softened to sphere alone. SECOND, the \
            machine: the doc's threshold is M-20's 4.75 ns/sample falling below \
            2.4, which is an Apple M5 figure. The M5 is reachable and CONTENDED \
            -- Spotlight, WindowServer, Messages and loginwindow at a combined \
            ~76% of a core, load average 1.65 to 1.87, 13 days up -- which is \
            exactly the contention M-005 is blocked on, and M-005 exists because \
            a memory-bound single-threaded timing taken beside a persistent \
            competitor is not a figure. So C1 is measured on the Zen 3 / Ryzen 9 \
            5900X host against THIS repository's own committed baseline: \
            docs/measurements/resolution_sweep-ryzen9-5900x.csv fits a marginal \
            13.1892 ns/sample over 9 rows from 16^3 to 256^3 for \
            marching_cubes/f32/sphere, so 2x is below 6.5946 ns/sample. The \
            ratio is measured WITHIN ONE BINARY AND ONE RUN, both loop shapes \
            compiled into the harness, because M-281 says a millisecond is a \
            property of the binary; and every row carries cycles and ghz because \
            M-280 says a nanosecond is not a unit on a governed CPU. (C2) THE \
            GATE, AND THE CLAUSE MOST LIKELY TO KILL THIS. All 216 golden hashes \
            are UNCHANGED. Vectorisation must not move one bit: IEEE elementwise \
            operations are exact per lane, so a hash movement means LLVM \
            reassociated something, and the change is REJECTED rather than \
            rebaselined. Unlike P-61's C2, here a movement is a defect and not a \
            cost. POPULATION: 216 rows, and the instrument is proven able to \
            report the bad news rather than assumed to be -- P-61 moved 135 of \
            these same 216 four commits ago. (C3) The f64 gain is at most half \
            the f32 gain. The doc's reason is that NEON is 2-wide at f64 against \
            4-wide at f32; on this host's AVX2 it is 4-wide against 8-wide, so \
            the ratio the clause asserts is the same factor of two and the \
            clause survives the machine change unaltered. POPULATION: both \
            scalars are measured on every field, so the clause fires on every \
            row pair. VERIFICATION REQUIREMENT, STATED AT REGISTRATION: \
            cargo-show-asm output for the monomorphised f32 instance goes in the \
            ticket. The crate is generic over Real and LLVM vectorises the \
            monomorphised instance or does not; that must be inspected per \
            instantiation, and a Criterion delta alone cannot distinguish a \
            vectorised loop from a lucky one.",
        falsified_by: "C1 under 2x, in which case the ceiling is the vector \
            width and the honest number is smaller than the prior suggests. On \
            gyroid a falsification is EXPECTED and its cause is named above; a \
            gyroid speedup at or above 2x would be the surprise and would mean \
            libm's sinf is being inlined and vectorised, which is a finding \
            about libm rather than about this loop. C2: ANY hash movement, which \
            rejects the change outright -- this is the one clause here whose \
            failure is not a measurement but a veto. C3: an f64 gain above half \
            the f32 gain, which would mean the f32 path was not the vector path \
            and C1's number came from something else. And the whole experiment \
            is void if the asm dump shows no vector instruction in the \
            monomorphised f32 sphere instance while C1 reports a speedup, \
            because then the speedup is the loop's bookkeeping and not its \
            arithmetic.",
        records: &[
            "arm",
            "field",
            "scalar",
            "samples_per_axis",
            "samples",
            "loop_shape",
            "ns_per_sample",
            "cycles_per_sample",
            "ghz",
            "speedup_vs_push",
            "marginal_ns_per_sample",
            "baseline_marginal_ns_per_sample",
            "bit_identical_to_push",
            "golden_hashes_unchanged",
            "vectorisable_body",
            "c1_speedup_f32",
            "c1_holds",
            "c2_holds",
            "c3_f64_over_f32_gain",
            "c3_holds",
            "machine",
            "wall_ms",
        ],
    },
    Preregistration {
        id: "P-71",
        ticket: "R-069",
        hypothesis: "M-167 is the largest single number this project owns about \
            its own GPU path: synchronisation was 83% of an extraction. M-159 \
            localised it -- the last four bytes cost 0.033 ms to move and 0.375 \
            ms to WAIT FOR, because poll(Wait) with no submission index drains \
            every dispatch queued before it -- and M-160 showed what removing it \
            buys: CPU time flat at about 0.17 ms from 33^3 to 129^3. What wgpu \
            29 gives splits the two targets and that split is the reason this is \
            registered as one experiment rather than two: PollType::Poll is \
            'check the device for a single time without blocking', PollType::Wait \
            is 'block until the given submission has completed execution', and \
            verbatim from the docs, 'On WebGPU, this has no effect. Callbacks \
            are invoked from the window event loop.' So native Bevy has a real \
            CPU stall to design away and the web build has no blocking \
            primitive at all -- meaning any code shaped around Wait is \
            native-only scaffolding, and the restructuring must not become a \
            cfg fork, for the same one-path reason the libm justification \
            already gives. TIMESTAMP_QUERY is in FeaturesWebGPU and is \
            supported on Vulkan, DX12, Metal, OpenGL and WebGPU -- the one \
            feature here that behaves identically on both targets, which V-48 \
            establishes and which is why C1 can be stated for both. TWO OF THE \
            THREE MECHANISMS ALREADY EXIST IN THE TREE, so this is partly a \
            measurement of shipped code rather than a build: extract_buffers \
            waits once, for the four bytes of the triangle count, and \
            extract_indirect waits NOT AT ALL, sizing the geometry from a budget \
            and turning the total into indirect draw arguments on the device. \
            (C1) INSTRUMENT FIRST. Timestamp queries attribute M-167's 83% into \
            submit / execute / map / copy, and the largest single component is \
            MAP-WAIT, not execute. POPULATION, probed on this host before this \
            registration was written: the adapter is an NVIDIA RTX 3090 on \
            Vulkan and advertises TIMESTAMP_QUERY, \
            TIMESTAMP_QUERY_INSIDE_PASSES and TIMESTAMP_QUERY_INSIDE_ENCODERS, \
            all three true, so the instrument exists on the machine the numbers \
            will come from. It is not currently reachable: \
            isomesh_gpu::headless::Gpu requests Features::empty(), and \
            ExtractTimings' own doc says in as many words that timestamp \
            attribution 'need[s] a device feature this crate does not request'. \
            Enabling it is a src/ change to isomesh-gpu and is registered as \
            one, shaped as a capability check that REFUSES LOUDLY on an adapter \
            without the feature rather than a fallback that silently reports \
            CPU-side numbers under a GPU-side column name -- GPU-007's pattern. \
            (C2) Feeding the vertex and index counts into indirect draw \
            arguments from a GPU buffer -- so the CPU never learns the count -- \
            removes at least 60% of the measured synchronisation at 129^3. \
            POPULATION: this is the difference between two entry points that \
            both already ship, extract_buffers against extract_indirect, so the \
            clause fires on every resolution both accept and the comparison is \
            WITHIN ONE BINARY AND ONE RUN, which M-281 requires. The \
            denominator is stated rather than assumed: 'the measured \
            synchronisation' is the wait component C1 attributes, not the whole \
            extraction, and if C1 finds map-wait is not the largest component \
            then C2's own denominator is smaller than M-167 suggests and the \
            entry must say so. (C3) An N-frame-delayed double-buffered staging \
            ring for the paths that genuinely need CPU-side data -- collider \
            generation -- holds the amortised per-frame cost within one chunk of \
            the budget across a 320x range, i.e. M-124's property survives the \
            added latency. POPULATION: M-124's own sweep, re-run with the ring \
            in place, so the clause fires on exactly the rows M-124 has. This \
            is the one genuinely new capability here. A DESIGN QUESTION THAT IS \
            THE OWNER'S, NOT THE HARNESS'S, AND IS REGISTERED AS A QUESTION: C3 \
            costs one to two frames of latency on collision. For a voxel game \
            that is invisible; for a CAD tool it is a decision. The \
            registration records the question rather than assuming an answer, \
            per CLAUDE.md's rule about design decisions, and the entry must \
            surface it rather than pick.",
        falsified_by: "C1: execute being the largest component, which would mean \
            the arithmetic did move after all and M-167's 'the arithmetic never \
            moved and was never the point' needs re-tiering. C2: under 60% of \
            the measured synchronisation removed. C3: the amortised cost \
            drifting outside one chunk of the budget, which means the ring \
            traded a stall for a queue. And the whole experiment is void if the \
            timestamp period reads zero or the resolved query set comes back \
            monotonically non-increasing, because then the attribution is a \
            column that was named and not measured -- the harness asserts on \
            both rather than reporting them.",
        records: &[
            "arm",
            "entry_point",
            "samples_per_axis",
            "cells",
            "triangles",
            "wall_ms",
            "submit_ms",
            "execute_ms",
            "map_wait_ms",
            "copy_ms",
            "largest_component",
            "synchronisation_ms",
            "synchronisation_share",
            "synchronisation_removed_share",
            "timestamp_feature",
            "timestamp_period_ns",
            "amortised_ms_per_frame",
            "budget_chunks",
            "within_one_chunk",
            "ring_frames_delay",
            "c1_holds",
            "c2_holds",
            "c3_holds",
            "adapter",
        ],
    },
    Preregistration {
        id: "P-72",
        ticket: "R-070",
        hypothesis: "The granularity of the active-cell structure is a \
            first-class parameter and it has never been swept. THE MEASURED \
            PRIOR IS A 256x SPREAD FROM ONE KNOB: Hoetzlein, GVDB, HPG 2016, \
            all timings on a Quadro M6000 -- tree build at 2048^3, <3,3,3,3> \
            616,444 bricks in 461 ms, <3,3,3,4> 83,218 bricks in 69.8 ms, \
            <3,3,3,6> 2,036 bricks in 1.8 ms. Same resolution, same data, 256x \
            apart, and the paper's own conclusion is 'larger brick sizes \
            produce a fewer number of bricks resulting in faster tree \
            changes.' Also measured there: octrees were 30-40% slower on node \
            insertion than N^3-trees. P-40 chose 64 cells per word and never \
            asked whether 64 was right; M-337 measured the stage at 5.5x and \
            12/12 bit-identical, which settles that the bitmap works and says \
            nothing about its granularity. P-39's Lipschitz brush pruning \
            (M-341, 3.36x median) is the direct analogue of GVDB's topology \
            cull, and GVDB's result says its yield should likewise be \
            granularity-dependent. WHAT THE CUBIC 8^3-64^3 KNOB IS IN THIS \
            CRATE, read from the source before registering rather than assumed: \
            it is the CHUNK, not the word. build_inside_bits packs 64 cells per \
            u64 along X ONLY -- a flat per-row word array with no block or \
            brick layer -- so there is no cubic granularity inside the bitmap \
            to sweep, and u8/u16/u32/u64 would be a 1-D sweep of a packing \
            width rather than GVDB's knob. The unit that is REBUILT ON AN EDIT \
            is the chunk: mark_edit marks chunks, DirtySet holds chunks, \
            mesh_dirty re-meshes chunks. That is GVDB's leaf brick with the \
            same semantics, and C1 is denominated in EDIT-PLUS-REMESH time, \
            which only the chunk path has. THE TRADEOFF IS ARITHMETIC, which is \
            why an optimum should exist at all: total world cells are held \
            fixed, so a chunk of c cells re-samples its shared corner planes at \
            ((c+1)/c)^3 of the field evaluations -- 1.42x at c = 8 against \
            1.05x at c = 64 -- while a finer dirty set re-meshes fewer cells \
            per edit. Small chunks pay in duplicated samples and save in wasted \
            remesh. (C1) Sweeping the chunk granularity across 8^3, 16^3, 32^3 \
            and 64^3 cells per pruning unit on a live edit trace, at a FIXED \
            total cell count, produces a PRONOUNCED OPTIMUM: the best and worst \
            edit-plus-remesh times differ by at least 2x. (C2) The optimum is \
            FIELD-DEPENDENT, and specifically differs between gyroid (surface \
            everywhere) and fbm_terrain (surface on a sheet). (C3) The spread is \
            SMALLER than GVDB's 256x, and predicted below 4x: GVDB's figure is \
            a tree rebuild on a 2048^3 SPH volume, far larger than a chunk \
            here, and its level-set numbers (5-6x over CPU) are consistently \
            weaker than its volume numbers (60x).",
        falsified_by: "C1 by a flat curve -- best-to-worst under 2x -- which is \
            a genuine null and means M-337's granularity was already at a \
            plateau; that null is worth having and is the expected outcome on \
            the smooth fields. C2 by one granularity winning on BOTH fields, \
            which would make it a constant rather than a parameter. C3 by a \
            spread above 4x, which would be the more valuable outcome and would \
            say the chunk-size regime is not the damping factor it looks like. \
            NO GOLDEN HASH MOVES: this experiment changes no extractor and no \
            placement -- it varies only how the same total cell count is \
            partitioned, and a partition that changed the mesh would be a \
            seam defect and is asserted against.",
        records: &[
            "field",
            "chunk_cells",
            "world_cells",
            "chunks",
            "samples_per_chunk",
            "sample_duplication",
            "edits",
            "dirty_chunks_total",
            "remeshed_cells_total",
            "mark_ms",
            "remesh_ms",
            "total_ms",
            "ms_per_edit",
            "vertices",
            "triangles",
            "distinct_surface_points",
            "vertex_duplication",
            "best_chunk_cells",
            "worst_chunk_cells",
            "spread",
            "c1_holds",
            "c2_holds",
            "c3_holds",
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
