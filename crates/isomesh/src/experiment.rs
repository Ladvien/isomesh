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
