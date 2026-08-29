//! Which coincidence collapses weld pieces the mesh kept apart.
//!
//! Ticket: R-053, hypothesis `P-125`. Registered **before** this module was
//! written — the `P-61` and `P-69` way — because a landing that was not
//! registered in advance is exactly what `V-45` cost.
//!
//! # What was missing
//!
//! `M-352` measured the `=`-corner repair on two real CT volumes. Custodio,
//! Pesco & Silva (`10.1186/s13173-019-0086-6`) name a third corner label for a
//! sample lying exactly on the isosurface, and what this crate needs from them
//! is that label's *consequence*: one shared vertex at the corner. On a cut edge
//! with an exactly-zero endpoint the interpolation parameter `t = a/(a − b)` is
//! exactly `0` or exactly `1`, so the crossing lands **on** the sample, every
//! cut edge meeting there places its own vertex at the same point, and the
//! vertex cache — keyed on the grid *edge* — shares none of them.
//!
//! The repair removed **every** degenerate triangle on both volumes and moved
//! **no geometry whatsoever** (`max_snap_distance` exactly `0`: the vertices are
//! already at the corner, so only indices merge). On `fuel` 64³ the topology was
//! untouched — χ 19 → 19, non-manifold edges 0 → 0, boundary edges 24 → 24. On
//! `bonsai` 256³ the identical repair took χ 517 → 585, non-manifold edges
//! **0 → 561** and boundary edges 4,366 → 3,716.
//!
//! So the repair is safe on one volume and changes the topology of the other,
//! and **which** is a property of the data rather than of the algorithm. The
//! deciding property is a graph property of the *baseline* mesh, computable
//! before any repair is applied — and until this module the crate had no way to
//! ask it. A caller could only be shipped a repair that is safe on somebody
//! else's scan, and would discover the difference as 520 welded pieces and a
//! moved Euler characteristic.
//!
//! # The predicate
//!
//! Two vertices a collapse is about to identify are in one of exactly two
//! situations, and everything here turns on which:
//!
//! * They **already share a triangle**. That triangle is one of the degenerates
//!   — two of its corners are the same point — so merging them flattens a fold.
//!   The triangle goes and no edge, boundary or component can move.
//! * They **share no triangle**. Then they lie on pieces of the surface that
//!   meet only at that point, the isosurface genuinely touches itself there, and
//!   identifying them **welds** those pieces together. That is a change of
//!   topology no relabelling avoids.
//!
//! Take the transitive closure of "shares a triangle with" *within* one
//! coincidence group and call the classes its **sharing clusters**. A group of
//! one cluster is a fold. A group of two or more is a **pinch**, and it joins
//! `clusters − 1` pieces. `M-352` counted **516 pinches of 17,201** collapse
//! groups on `bonsai`, joining **520** pieces, against **0 of 50** on `fuel`.
//!
//! [`PinchReport::is_pinch_free`] is therefore the precondition a caller can
//! test on its own data: when it holds, the collapse removes degenerate
//! triangles and provably cannot move a component, a boundary or an edge count.
//!
//! # A group is an ε-connected component, and the alternatives are both wrong
//!
//! ε-closeness is **not transitive**, so "every pair within ε" is not a
//! partition at all and there would be nothing to count groups over — the same
//! obstacle [`weld`](crate::weld) states in its own rule. Its **transitive
//! closure** is a partition, and that is what a group is here: vertices joined
//! by a chain of pairs each within
//! [`ValidateConfig::weld_epsilon`](super::ValidateConfig::weld_epsilon).
//!
//! The other candidate was equality on the lattice cell key —
//! [`MeshReport::weld_buckets`](super::MeshReport::weld_buckets)' own notion, and
//! also an equivalence relation. It was **measured and rejected**. On `sphere` at
//! 25³, where `M-48` recorded [`Welder::weld`](crate::weld::Welder::weld)
//! removing **48 vertices and 96 triangles** at this very epsilon, the bucket
//! reading finds 17 groups over 47 vertices — **30 of the 48 merges, and 60 of
//! the 96 folding faces** — because a coincidence class straddling a cell face is
//! split between two buckets. The closure finds 24 groups over 72 vertices:
//! **48 and 96, exactly the weld**. A census that under-reports a collapse is a
//! **false clearance**, and issuing clearances is this report's entire job;
//! `tests::the_census_predicts_the_weld_on_a_real_extraction` is where both
//! readings were measured.
//!
//! The closure errs the other way, which is the safe way. The welder joins a
//! vertex to the *lowest-indexed representative* within ε and therefore stops a
//! chain — `a ~ b ~ c` with `a` and `c` further than ε apart yields two
//! representatives — so the welder's classes **refine** the closure. A
//! pinch-free census is a genuine clearance for a weld at that epsilon, and a
//! non-zero pinch count is an upper bound on what a weld could reach.
//!
//! The lattice and the 27-cell neighbourhood are `crate::weld::Lattice` and
//! `validate`'s own duplicate-vertex probe, **not copies of them**. Two
//! quantisers that disagreed about which cell a position falls in, or two probes
//! that disagreed about which cells are candidates, would make this census
//! describe a neighbourhood no weld ever visited. A non-finite coordinate
//! saturates into the origin cell and then fails the exact distance test, so it
//! is never joined to anything — the same behaviour, for the same reason, as
//! `duplicate_vertices`.
//!
//! # The cluster label is canonical, and that is why the union-find is local
//!
//! `super::Dsu` unions by size, and its determinism rests on *its* callers
//! requesting unions in sorted order. The unions in the second phase here arrive
//! in **face order**, which `P-125`'s C2 permutes 128 ways on purpose — that is
//! `✗26`'s objection asked in advance rather than after a landing. Under union
//! by size the surviving root depends on the order the unions arrived in, so
//! [`PinchGroups::clusters`] would carry a face-order leak into a public
//! artefact, and `validate`'s whole no-hash-map argument is that a report is a
//! pure function of the mesh's values.
//!
//! So the union-find here keeps the **lower** index as the root. Every set's
//! root is then its least member, the second phase only ever joins two members
//! of one group, and the label is a canonical vertex index whichever order the
//! faces arrived in. It also makes the cluster count free: a cluster has exactly
//! one root, that root is a member, so the number of clusters in a group is the
//! number of members that are their own label — no per-group allocation, which
//! is the other half of C2's falsifier.
//!
//! # What this does not claim
//!
//! [`pieces_joined`](PinchReport::pieces_joined) counts pieces **the group's own
//! members lie on**, not connected components of the whole mesh. Two sharing
//! clusters that meet at a pinch may be joined elsewhere by a long path of
//! triangles, in which case identifying them welds nothing globally. That gap is
//! `P-125`'s C3 and is measured by the harness, not asserted here; this report is
//! deliberately local, because a global component count is extra work a caller
//! may not want and is already
//! [`MeshReport::components`](super::MeshReport::components) — which counts
//! **edge**-connected components of the *face* set and is therefore a different
//! relation from the vertex-level one the predicate uses.
//!
//! Nothing here claims anything about position. `M-352`'s repair moved no
//! geometry, so this is a pure connectivity decision; a caller who wants the
//! geometric half asks
//! [`groups_moving_geometry`](PinchReport::groups_moving_geometry), which is zero
//! exactly when every group's members are bit-identical.
//!
//! # Cost
//!
//! `O(V log V + V·k + F)` time and `O(V)` memory, for `k` candidates in the 27
//! probed cells. The `log V` is the sort, and it is the module's own convention
//! rather than an accident: `validate`'s header states why every structural pass
//! here is *sort a flat `Vec`, then scan runs of equal keys* — one allocation
//! instead of a node per entry, cache-linear, and deterministic by construction
//! because each sort key ends in a vertex index and no two entries ever compare
//! equal. `P-125`'s hypothesis text says `O(V + F)`; that is the sort-free
//! reading and this module does not have it. Its falsifier gates on face-order
//! dependence and on per-group allocation, neither of which the sort introduces.
//!
//! **Seven buffers, all `O(V)` or smaller, none per group**: the sorted lattice
//! keys, the union-find parents (reset between the two phases rather than
//! reallocated), the roots to group by, the group id per vertex, and the three
//! the caller gets back. Every one is `Vec::with_capacity`-reserved to a size
//! known before its first push — the member and group counts come from one scan
//! of the sorted roots — so no buffer ever grows and the count does not move with
//! the number of groups. `crates/isomesh/Cargo.toml` forbids `unsafe_code`
//! workspace-wide, so a counting allocator cannot exist in this crate's own
//! benches; what stands in its place is that `len` equals `capacity` on all three
//! returned buffers, which is only true of a `Vec` that never reallocated.

#[cfg(test)]
mod tests;

use alloc::vec::Vec;
use core::cmp::Ordering;
use core::fmt;

use super::ValidateConfig;
use crate::Real;
use crate::weld::Lattice;

/// What a coincidence collapse would do to a mesh's connectivity.
///
/// Produced by [`pinch_census`]. Every field is a count, in keeping with the
/// rest of [`validate`](super) — [`is_pinch_free`](Self::is_pinch_free) and
/// [`moves_no_geometry`](Self::moves_no_geometry) are the two derived predicates
/// and both are opt-in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PinchReport {
    /// Vertices in the buffer, referenced or not.
    pub vertices: u64,
    /// Whole triangles in the index buffer.
    pub triangles: u64,
    /// Indices past the last whole triangle. Non-zero is a malformed buffer,
    /// reported rather than rejected.
    pub trailing_indices: u64,
    /// Faces excluded because an index is out of range or repeated.
    ///
    /// A face that cannot be dereferenced has no edges to contribute, and one
    /// naming a vertex twice has no well-defined edge set — the same intake rule
    /// [`validate_features`](super::validate_features) applies.
    pub faces_skipped: u64,

    /// Coincidence groups: ε-connected components holding **two or more**
    /// vertices.
    ///
    /// These are the collapses a weld at this epsilon could perform, and the
    /// denominator every figure below is out of. `M-352` measured **17,201** on
    /// `bonsai` 256³ at isovalue 32 and **50** on `fuel` 64³.
    pub collapse_groups: u64,
    /// Vertices in those groups. `collapsing_vertices − collapse_groups` is how
    /// many vertices the collapse removes.
    pub collapsing_vertices: u64,

    /// Groups whose members do **not** all already share triangles.
    ///
    /// **The number the safety of a collapse turns on.** Each of these joins
    /// pieces of the surface that the baseline mesh kept apart. `M-352` measured
    /// **516** on `bonsai` against **0** on `fuel`, at 3.0% and 0.0% of their
    /// collapse groups.
    pub pinch_groups: u64,
    /// Summed `clusters − 1` over the pinch groups: how many pieces the collapse
    /// joins in total.
    ///
    /// Larger than [`pinch_groups`](Self::pinch_groups) exactly when some group
    /// spans three or more sharing clusters — **520 against 516** on `bonsai`,
    /// so four of its groups do. Local to the groups: see the module docs for
    /// what this is not.
    pub pieces_joined: u64,

    /// Faces the collapse would drop, having two corners in one group.
    ///
    /// These are the degenerate triangles the repair exists to remove. `M-352`
    /// measured `triangles_dropped` **164** on `fuel` and **58,097** on
    /// `bonsai`, taking `degenerate_triangles` to 0 on both.
    pub folding_faces: u64,
    /// Face-side edges joining two members of one group.
    ///
    /// The unions the sharing clusters were built from. **Read this before
    /// believing a zero in [`pinch_groups`](Self::pinch_groups)**: a mesh whose
    /// groups are all folds has to have found sharing edges to say so, and zero
    /// here with non-zero `collapse_groups` means every group came out a pinch
    /// by default.
    pub sharing_edges: u64,

    /// Groups whose members are not all at bit-identical positions.
    ///
    /// Zero means the collapse is **purely combinatorial** — it merges indices
    /// and moves nothing — which is what `M-352` measured as `max_snap_distance`
    /// exactly `0` on both volumes, and is why the pinch question is a decision
    /// rather than an approximation. Non-zero means the group spans a
    /// neighbourhood rather than a point, and the caller is moving geometry.
    pub groups_moving_geometry: u64,
}

impl PinchReport {
    /// No collapse joins pieces the mesh kept apart.
    ///
    /// **The precondition for applying an `=`-corner or coincidence repair
    /// blind.** When this holds, the collapse removes
    /// [`folding_faces`](Self::folding_faces) degenerate triangles and cannot
    /// change a component count, a boundary edge or an Euler characteristic.
    /// When it does not, it welds [`pieces_joined`](Self::pieces_joined) pieces
    /// and the caller has to decide whether that is a repair or a defect —
    /// which is a question about the caller's data, not about the algorithm.
    #[must_use]
    pub const fn is_pinch_free(&self) -> bool {
        self.pinch_groups == 0
    }

    /// The collapse merges indices and moves no vertex.
    ///
    /// Every group's members are bit-identical, so nothing is snapped anywhere.
    /// Vacuously true of a mesh with no coincidences at all, which is why it is
    /// read beside [`collapse_groups`](Self::collapse_groups) rather than alone.
    #[must_use]
    pub const fn moves_no_geometry(&self) -> bool {
        self.groups_moving_geometry == 0
    }
}

impl fmt::Display for PinchReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "pinch census: {} vertices, {} triangles, {} coincidence groups over {} vertices",
            self.vertices, self.triangles, self.collapse_groups, self.collapsing_vertices
        )?;
        writeln!(
            f,
            "  pinches {} of {}   pieces joined {}   folds drop {} faces   sharing edges {}",
            self.pinch_groups,
            self.collapse_groups,
            self.pieces_joined,
            self.folding_faces,
            self.sharing_edges
        )?;
        write!(
            f,
            "  groups moving geometry {}   skipped faces {} (+{} trailing indices)   -> {}",
            self.groups_moving_geometry,
            self.faces_skipped,
            self.trailing_indices,
            if self.is_pinch_free() {
                "PINCH-FREE"
            } else {
                "WELDS PIECES"
            }
        )
    }
}

/// **Which** vertices pinch, rather than only how many.
///
/// Produced by [`pinch_features`]. A count cannot be drawn, and a caller told
/// its scan has 516 pinches needs to know where they are before it can decide
/// anything — the same argument
/// [`NonManifoldFeatures`](super::NonManifoldFeatures) exists for.
///
/// Compressed-row storage over the coincidence groups. Groups appear in
/// ascending order of their least member and members ascend within a group, so
/// this is a pure function of the mesh's values: no hash seed, no insertion
/// order, and no dependence on the order the faces arrived in.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PinchGroups {
    /// Every coincidence group's member vertices, concatenated.
    pub vertices: Vec<u32>,
    /// Group `g` occupies `vertices[starts[g] .. starts[g + 1]]`.
    ///
    /// One entry per group plus a terminator, so this is empty exactly when
    /// there are no groups.
    pub starts: Vec<u32>,
    /// Per entry of [`vertices`](Self::vertices), the label of its sharing
    /// cluster: the **least** vertex index among the members it already shares
    /// triangles with, transitively.
    ///
    /// Two members of one group carrying the same label are already connected
    /// through the faces; a group carrying two distinct labels is a pinch. The
    /// label is canonical rather than an arbitrary union-find root — see the
    /// module docs for why that is load-bearing.
    pub clusters: Vec<u32>,
}

impl PinchGroups {
    /// How many coincidence groups there are.
    #[must_use]
    pub fn group_count(&self) -> usize {
        self.starts.len().saturating_sub(1)
    }

    /// Group `g`'s member vertices, ascending.
    ///
    /// # Panics
    ///
    /// If `g` is not below [`group_count`](Self::group_count).
    #[must_use]
    pub fn members(&self, g: usize) -> &[u32] {
        let (lo, hi) = self.span(g);
        &self.vertices[lo..hi]
    }

    /// Group `g`'s sharing-cluster labels, aligned with
    /// [`members`](Self::members).
    ///
    /// # Panics
    ///
    /// If `g` is not below [`group_count`](Self::group_count).
    #[must_use]
    pub fn clusters_of(&self, g: usize) -> &[u32] {
        let (lo, hi) = self.span(g);
        &self.clusters[lo..hi]
    }

    /// How many sharing clusters group `g` spans.
    ///
    /// `1` is a fold, `2` or more is a pinch joining that many pieces less one.
    /// Counted as the members that are their own label, which is exact because a
    /// cluster's root is its least member and every root is a member.
    ///
    /// # Panics
    ///
    /// If `g` is not below [`group_count`](Self::group_count).
    #[must_use]
    pub fn clusters_in(&self, g: usize) -> usize {
        self.members(g)
            .iter()
            .zip(self.clusters_of(g))
            .filter(|(v, label)| v == label)
            .count()
    }

    /// Group `g`'s members do not all already share triangles.
    ///
    /// # Panics
    ///
    /// If `g` is not below [`group_count`](Self::group_count).
    #[must_use]
    pub fn clusters_are_split(&self, g: usize) -> bool {
        self.clusters_in(g) > 1
    }

    fn span(&self, g: usize) -> (usize, usize) {
        assert!(
            g < self.group_count(),
            "group {g} of {}",
            self.group_count()
        );
        (self.starts[g] as usize, self.starts[g + 1] as usize)
    }
}

/// Census the coincidence collapses a weld at this epsilon could perform, and
/// say which of them weld pieces the mesh kept apart.
///
/// Takes slices rather than a [`MeshBuffer`](crate::MeshBuffer) for the same
/// reason [`validate_indexed`](super::validate_indexed) does: the callers that
/// matter do not have one. Nothing is required of the mesh — it need not be
/// welded, closed, manifold or even well-formed — and nothing here panics on a
/// malformed one; an unusable face is skipped and counted.
///
/// The epsilon is
/// [`ValidateConfig::weld_epsilon`](super::ValidateConfig::weld_epsilon), so the
/// same config that governs the validator's duplicate count and
/// [`weld`](crate::weld)'s tolerance governs this. There is no invalid value to
/// guard against: a `ValidateConfig` that exists is a valid one.
///
/// # Reading it
///
/// [`PinchReport::is_pinch_free`] is the whole question for a caller about to
/// apply a repair. Everything else is the accounting behind it, and
/// [`pinch_features`] hands back the vertices so the answer can be drawn.
#[must_use]
pub fn pinch_census<R: Real>(
    positions: &[[R; 3]],
    indices: &[u32],
    cfg: &ValidateConfig,
) -> PinchReport {
    pinch_features(positions, indices, cfg).0
}

/// [`pinch_census`], plus **which** vertices are in each group and which
/// cluster.
///
/// Not a second implementation. It **is** the implementation; [`pinch_census`]
/// calls it and drops the second half, exactly as
/// [`validate_indexed`](super::validate_indexed) does to
/// [`validate_features`](super::validate_features). The lists are collected
/// unconditionally, which costs nothing on the common path: a mesh with no
/// coincidences leaves all three empty, and an empty [`Vec`] does not allocate.
#[must_use]
pub fn pinch_features<R: Real>(
    positions: &[[R; 3]],
    indices: &[u32],
    cfg: &ValidateConfig,
) -> (PinchReport, PinchGroups) {
    let n = positions.len();
    let whole = indices.len() - indices.len() % 3;
    let mut report = PinchReport {
        vertices: n as u64,
        triangles: (whole / 3) as u64,
        trailing_indices: (indices.len() % 3) as u64,
        faces_skipped: 0,
        collapse_groups: 0,
        collapsing_vertices: 0,
        pinch_groups: 0,
        pieces_joined: 0,
        folding_faces: 0,
        sharing_edges: 0,
        groups_moving_geometry: 0,
    };
    if n == 0 {
        return (report, PinchGroups::default());
    }

    // ── phase 1: the ε-graph, closed ────────────────────────────────────────
    //
    // Bucket onto the lattice `weld` welds on, sort, then probe the 27
    // neighbouring cells and union every pair that is genuinely within epsilon.
    // Quantising sidesteps `R: !Ord` entirely — the key is `[i64; 3]`, which has
    // a total order — and `-0.0` and `+0.0` land in the same cell. A non-finite
    // coordinate saturates into the origin cell and then fails the distance
    // test, so it is never joined to anything.
    let eps = R::from_f64(cfg.weld_epsilon());
    let eps_sq = eps * eps;
    let lattice = Lattice::new(positions, eps);
    let mut cells: Vec<([i64; 3], u32)> = Vec::with_capacity(n);
    for (v, p) in positions.iter().enumerate() {
        cells.push((lattice.key_of(*p), v as u32));
    }
    // The vertex index is part of the key and is unique, so no two entries
    // compare equal and an unstable sort is a deterministic one.
    cells.sort_unstable();

    let mut parent: Vec<u32> = (0..n as u32).collect();
    for &(base, v) in &cells {
        let p = positions[v as usize];
        for dz in -1..=1i64 {
            for dy in -1..=1i64 {
                for dx in -1..=1i64 {
                    let key = [base[0] + dx, base[1] + dy, base[2] + dz];
                    let from = cells.partition_point(|(k, _)| *k < key);
                    for &(k, u) in &cells[from..] {
                        if k != key {
                            break;
                        }
                        // Each unordered pair once, and never a vertex against
                        // itself.
                        if u >= v {
                            continue;
                        }
                        let q = positions[u as usize];
                        let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
                        if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= eps_sq {
                            union_to_lower(&mut parent, u, v);
                        }
                    }
                }
            }
        }
    }

    // ── the groups, counted before anything is reserved ─────────────────────
    //
    // One scan of the sorted roots for the sizes, one to fill: that is what lets
    // every buffer below be reserved exactly once rather than growing per group.
    // Sorting by `(root, vertex)` puts each group's members together and
    // ascending, and the root is the least member, so groups come out ordered by
    // their least member.
    let mut roots: Vec<(u32, u32)> = Vec::with_capacity(n);
    for v in 0..n as u32 {
        roots.push((find(&mut parent, v), v));
    }
    roots.sort_unstable();

    let mut group_count = 0usize;
    let mut member_count = 0usize;
    let mut scan = Runs::new(&roots);
    while let Some((lo, hi)) = scan.next_group() {
        group_count += 1;
        member_count += hi - lo;
    }
    report.collapse_groups = group_count as u64;
    report.collapsing_vertices = member_count as u64;

    let mut group_of = alloc::vec![u32::MAX; n];
    let mut members: Vec<u32> = Vec::with_capacity(member_count);
    let mut starts: Vec<u32> = Vec::with_capacity(group_count + 1);
    if group_count > 0 {
        starts.push(0);
    }
    let mut scan = Runs::new(&roots);
    while let Some((lo, hi)) = scan.next_group() {
        let g = starts.len() as u32 - 1;
        for &(_, v) in &roots[lo..hi] {
            group_of[v as usize] = g;
            members.push(v);
        }
        starts.push(members.len() as u32);
    }

    // ── phase 2: the sharing clusters, over the baseline faces ──────────────
    //
    // The same parent array, reset rather than reallocated. Only edges whose two
    // ends are in the *same* group are unioned, so every set is confined to one
    // group and its root is a member of it. Group membership is read off
    // `group_of` rather than recomputed, which is what keeps this pass linear in
    // the faces.
    parent.clear();
    parent.extend(0..n as u32);
    for tri in indices[..whole].as_chunks::<3>().0 {
        let unusable = tri.iter().any(|&i| i as usize >= n)
            || tri[0] == tri[1]
            || tri[1] == tri[2]
            || tri[0] == tri[2];
        if unusable {
            report.faces_skipped += 1;
            continue;
        }
        let mut folds = false;
        for (a, b) in [(0usize, 1usize), (1, 2), (2, 0)] {
            let (u, w) = (tri[a], tri[b]);
            let g = group_of[u as usize];
            if g != u32::MAX && g == group_of[w as usize] {
                report.sharing_edges += 1;
                folds = true;
                union_to_lower(&mut parent, u, w);
            }
        }
        if folds {
            report.folding_faces += 1;
        }
    }

    let mut clusters: Vec<u32> = Vec::with_capacity(member_count);
    for &v in &members {
        clusters.push(find(&mut parent, v));
    }
    let groups = PinchGroups {
        vertices: members,
        starts,
        clusters,
    };

    // ── the verdict per group ───────────────────────────────────────────────
    for g in 0..groups.group_count() {
        let spanned = groups.clusters_in(g);
        if spanned > 1 {
            report.pinch_groups += 1;
            report.pieces_joined += spanned as u64 - 1;
        }
        let group = groups.members(g);
        let at = positions[group[0] as usize];
        if group[1..]
            .iter()
            .any(|&v| !same_point(positions[v as usize], at))
        {
            report.groups_moving_geometry += 1;
        }
    }

    (report, groups)
}

/// Runs of equal first element in a sorted slice, yielding only the ones with
/// two or more members.
///
/// A cursor rather than a closure so that the counting scan and the filling scan
/// are the *same* traversal rule, written once. Two copies of a run-detection
/// loop that disagreed by one would put a group in the census that is not in the
/// feature list.
struct Runs<'a> {
    of: &'a [(u32, u32)],
    at: usize,
}

impl<'a> Runs<'a> {
    const fn new(of: &'a [(u32, u32)]) -> Self {
        Self { of, at: 0 }
    }

    fn next_group(&mut self) -> Option<(usize, usize)> {
        while self.at < self.of.len() {
            let lo = self.at;
            let mut hi = lo + 1;
            while hi < self.of.len() && self.of[hi].0 == self.of[lo].0 {
                hi += 1;
            }
            self.at = hi;
            if hi - lo > 1 {
                return Some((lo, hi));
            }
        }
        None
    }
}

/// Root of `x`, halving the path on the way.
fn find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        let grand = parent[parent[x as usize] as usize];
        parent[x as usize] = grand;
        x = grand;
    }
    x
}

/// Join two sets, keeping the **lower** index as the root.
///
/// Not union by size, and the module docs say why: the root is reported, so it
/// has to be canonical under any order the unions arrive in.
fn union_to_lower(parent: &mut [u32], a: u32, b: u32) {
    let (ra, rb) = (find(parent, a), find(parent, b));
    if ra == rb {
        return;
    }
    let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
    parent[hi as usize] = lo;
}

/// Bit-identical, componentwise.
///
/// `Real::total_cmp` rather than `==`, so signed zeros are distinguished and a
/// NaN coordinate compares equal only to itself — the question is whether the
/// collapse moves anything, and a comparison that answered `false` for a
/// position against itself would report movement that is not there.
fn same_point<R: Real>(a: [R; 3], b: [R; 3]) -> bool {
    a.iter()
        .zip(&b)
        .all(|(x, y)| x.total_cmp(y) == Ordering::Equal)
}
