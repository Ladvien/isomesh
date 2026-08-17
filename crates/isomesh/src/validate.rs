//! Mesh validity metrics.
//!
//! This module **reports**; it does not judge. Every metric is a count, nothing
//! panics on a malformed mesh, and [`validate_indexed`] never returns an error.
//! That is deliberate: a `Result` forces a branch at every call site, and `?`
//! would throw away the partial report exactly when it is most useful — a mesh
//! with one bad index is a mesh you still want the manifoldness of.
//!
//! Loudness comes from [`MeshReport::has_structural_errors`], from the banner
//! [`MeshReport`]'s `Display` prints, and from the opt-in
//! [`MeshReport::panic_if_invalid`].
//!
//! # Why there is no hash map here
//!
//! `no_std` plus `alloc` offers `BTreeMap` and no `HashMap`, and the crate's
//! dependency policy rules out `hashbrown`. But a B-tree is the wrong shape
//! anyway. Every structural pass below is *sort a flat `Vec`, then scan runs of
//! equal keys*: one allocation instead of a node per entry, cache-linear instead
//! of pointer-chasing, and — the part that matters — **deterministic by
//! construction**. Each sort key ends in a face or vertex index, so no two
//! entries ever compare equal, which means there is no iteration order to leak:
//! not a hash seed, not insertion order, not even the instability of
//! `sort_unstable`. A report is a pure function of the mesh's values, which is
//! what makes it safe to feed a golden hash.
//!
//! # Cost
//!
//! `O(F log F)` time and `O(F)` memory. Validation is explicitly not on the
//! per-chunk re-meshing hot path, so it allocates freely and is not built around
//! a reusable scratch buffer. If a benchmark ever shows otherwise, that change
//! arrives with the benchmark.

use alloc::vec::Vec;
use core::fmt;

use crate::{MeshBuffer, Real};

mod accuracy;
mod determinism;
mod field_bound;
mod isotopy;
mod sealing;
mod self_intersection;
mod tri_grid;

pub use accuracy::{AccuracyConfig, AccuracyReport, DistanceStats, accuracy};
pub use determinism::{DeterminismReport, Divergence, RunPair, check_determinism};
pub use field_bound::{EIKONAL_TOLERANCE, FieldBoundReport, field_bound_report};
pub use isotopy::{IsotopyReport, cell_is_certified, isotopy_report};
pub use sealing::{SealingReport, sealing};
pub use self_intersection::{SelfIntersectionReport, self_intersections};

/// Thresholds for the two metrics that have units.
///
/// There is deliberately **no `Default`**. An absolute area or weld distance is
/// meaningless without a length scale: a threshold tuned at 64³ silently
/// misfires at 256³, in the direction of reporting nothing. The caller states
/// the grid spacing and both thresholds follow from it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValidateConfig {
    cell_size: f64,
    weld_epsilon: f64,
    area_epsilon_rel: f64,
}

impl ValidateConfig {
    /// `weld_epsilon / cell_size`.
    ///
    /// Vertices that *should* coincide are the ones computed twice at a chunk
    /// seam. Computed by identical code from identical inputs they are
    /// bit-identical; computed from mathematically equal but differently ordered
    /// expressions they differ by about an ulp relative, which at a coordinate
    /// magnitude of a few dozen cells is around `2e-6·h`. This sits comfortably
    /// above that and comfortably below any real feature.
    pub const WELD_EPSILON_REL: f64 = 1e-4;

    /// Degenerate-area threshold, relative to `cell_size²`.
    ///
    /// Relative rather than absolute precisely because area has units. The value
    /// puts the threshold about thirty times above the `f32` rounding floor for
    /// a cross product of magnitude `h²` (~`6e-8·h²`), so it means "numerically
    /// zero" without firing on ordinary rounding. An order of magnitude smaller
    /// would sit *below* the noise floor and therefore mean nothing at all.
    pub const AREA_EPSILON_REL: f64 = 1e-6;

    /// The only constructor. Derives both thresholds from the grid spacing.
    ///
    /// Fields are private and this is the sole way in, so a `ValidateConfig`
    /// that exists is a valid one. That is why [`validate_indexed`] needs no
    /// runtime check of its own: the invalid state is unrepresentable rather
    /// than merely reported.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidCellSize`](crate::Error::InvalidCellSize) if `cell_size`
    /// is not finite and positive. Every threshold here is relative to the
    /// spacing, so a meaningless spacing makes them all meaningless.
    pub fn from_cell_size(cell_size: f64) -> crate::Result<Self> {
        if !cell_size.is_finite() || cell_size <= 0.0 {
            return Err(crate::Error::InvalidCellSize { value: cell_size });
        }
        Ok(Self {
            cell_size,
            weld_epsilon: cell_size * Self::WELD_EPSILON_REL,
            area_epsilon_rel: Self::AREA_EPSILON_REL,
        })
    }

    /// Grid spacing `h`, in world units. Every threshold is relative to it.
    #[must_use]
    pub fn cell_size(&self) -> f64 {
        self.cell_size
    }

    /// Two vertices closer together than this count as duplicates.
    #[must_use]
    pub fn weld_epsilon(&self) -> f64 {
        self.weld_epsilon
    }

    /// A triangle is degenerate when `area <= area_epsilon_rel · cell_size²`.
    #[must_use]
    pub fn area_epsilon_rel(&self) -> f64 {
        self.area_epsilon_rel
    }
}

/// A report on an indexed triangle mesh. Every field is a count; nothing here
/// asserts.
///
/// Counts are `u64` rather than `usize` so that the `Display` block is identical
/// on 32- and 64-bit targets — the block is stable enough to hash.
///
/// # Reading it
///
/// "Valid" is not one thing, so there are three predicates rather than one:
/// [`has_structural_errors`](Self::has_structural_errors) for malformed input,
/// [`is_manifold`](Self::is_manifold) for a surface that may have boundary, and
/// [`is_closed`](Self::is_closed) for one that may not. Which of the last two
/// applies is a property of the field being meshed, not of the mesher.
///
/// **If you are asking "did this pass", use [`satisfies`](Self::satisfies) and
/// name the [`SurfaceGate`] your artefact earns** — that is the rule, and
/// picking a predicate by intuition is how a correct mesh comes to read as
/// broken (✗22). Reaching for the three predicates directly is right when you
/// want to *describe* a mesh rather than judge it; `manifold_check` does exactly
/// that to print what it found.
///
/// # Solids and surfaces are measured differently, on purpose
///
/// A closed solid must have no boundary edges. An open surface — an open field,
/// one chunk of a larger extraction, a render mesh that is a subset of some
/// body — is **supposed** to have them, and its
/// [`boundary_edges`](Self::boundary_edges) count is a recorded number rather
/// than a failure. The same split runs through
/// [`violations`](Self::violations), which deliberately excludes the metrics a
/// correct extractor produces non-zero counts of.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeshReport {
    // ── census ──────────────────────────────────────────────────────────────
    /// `positions.len()`.
    pub vertices: u64,
    /// Vertices referenced by at least one considered face.
    pub referenced_vertices: u64,
    /// Distinct undirected edges `{u, v}` over considered faces.
    pub edges: u64,
    /// Faces included in the topological metrics.
    pub faces: u64,
    /// Faces excluded from them: an out-of-range or a repeated index.
    pub faces_skipped: u64,

    // ── topology ────────────────────────────────────────────────────────────
    /// `referenced_vertices − edges + faces`.
    ///
    /// Uses *referenced* vertices, not `positions.len()`. A vertex no face
    /// mentions is not part of the surface, and counting it would inflate `χ` by
    /// one each — so a mesh with a stale vertex pool would report `χ = 3` for a
    /// sphere and read as a topology bug when it is an allocation bug.
    /// [`unreferenced_vertices`](Self::unreferenced_vertices) is reported
    /// separately so the total is still recoverable.
    pub euler_characteristic: i64,
    /// Edge-connected components of the face set.
    pub components: u64,
    /// Connected components of the boundary-edge subgraph.
    ///
    /// This is what makes `χ` readable. `χ = 1` alone is ambiguous; `χ = 1` with
    /// one boundary loop is unambiguously a disk. A single count of boundary
    /// edges could not distinguish one loop from three, and it is the number of
    /// loops that enters `χ = 2 − 2g − b`.
    pub boundary_loops: u64,
    /// `(2 − χ − boundary_loops) / 2`, when that formula applies.
    ///
    /// `Some` only for a single, consistently oriented, manifold component with
    /// no skipped faces. The orientation precondition is what makes it sound: a
    /// closed manifold whose every edge is traversed oppositely by its two faces
    /// is orientable, so the formula holds and the division is exact. Without
    /// that check it would silently produce a half-integer on a Klein bottle.
    pub genus: Option<i64>,

    // ── violations ──────────────────────────────────────────────────────────
    /// Undirected edges incident to **three or more** faces.
    ///
    /// Note the split from [`boundary_edges`](Self::boundary_edges). Counting
    /// every edge with "not exactly two" faces would double-count and would make
    /// zero unachievable for any open mesh — including every individual chunk,
    /// which is the case the game path cares most about.
    pub non_manifold_edges: u64,
    /// Vertices whose incident faces do not form a single connected fan.
    ///
    /// Catches both the bowtie (two cones sharing an apex) and umbrella
    /// branching. See the module source for why the cheap "incident faces equals
    /// incident edges" test does not.
    pub non_manifold_vertices: u64,
    /// Undirected edges incident to exactly one face.
    pub boundary_edges: u64,
    /// Edges whose two incident faces traverse them in the *same* direction.
    ///
    /// Not on T-001's original list, and worth its keep: a transcribed case
    /// table with one flipped triangle passes the Euler check, edge
    /// manifoldness *and* vertex manifoldness while being inside out.
    pub inconsistently_oriented_edges: u64,
    /// Triangles with `area <= area_epsilon_rel · cell_size²`.
    ///
    /// A recorded metric rather than a violation to gate on. Marching cubes
    /// genuinely emits slivers whenever a grid corner value sits near zero;
    /// that is the algorithm, not a defect.
    pub degenerate_triangles: u64,
    /// Triangles with two or three equal indices. Also counted in
    /// [`faces_skipped`](Self::faces_skipped).
    ///
    /// Separate from [`degenerate_triangles`](Self::degenerate_triangles)
    /// because they are different bugs: this one is a case table emitting a
    /// collapsed triangle, that one is two edge crossings converging.
    pub repeated_index_triangles: u64,
    /// Distinct cells the weld lattice put vertices in.
    ///
    /// A diagnostic rather than a violation, and it exists because the duplicate
    /// scan's *cost* has a failure mode its *answer* does not. The lattice
    /// quantises through `as_f32`, which stops distinguishing consecutive
    /// integers above `2²⁴`; if the scale ever pushed coordinates past that,
    /// cells would silently merge, every vertex would land in a handful of
    /// buckets and the 27-cell probe would degrade toward comparing everything
    /// with everything — while still returning the right count, because the
    /// exact distance test runs regardless.
    ///
    /// Quantising relative to the mesh's own bounds is what prevents that
    /// (T-008, M-18), and this is the number that proves it: on a well-spread
    /// mesh it tracks the vertex count, and a collapse shows up here and
    /// nowhere else.
    pub weld_buckets: u64,

    /// Vertices having an earlier vertex within `weld_epsilon`.
    ///
    /// Equivalently, how many a first-fit welder could remove. Defined against
    /// *earlier* vertices rather than as equivalence classes because
    /// epsilon-closeness is not transitive, so classes are not well defined.
    pub duplicate_vertices: u64,
    /// Vertices no considered face references.
    pub unreferenced_vertices: u64,

    // ── malformed input ─────────────────────────────────────────────────────
    /// Individual index values `>= positions.len()`.
    pub out_of_range_indices: u64,
    /// `indices.len() % 3`: a partial trailing triangle.
    pub trailing_indices: u64,
    /// Positions with a non-finite component.
    ///
    /// Present because such a vertex would otherwise pass silently: it
    /// quantises into one bucket and never matches as a duplicate, and
    /// `NaN <= threshold` is false so it does not register as degenerate either.
    pub non_finite_positions: u64,
    /// `normals.len() != positions.len()`. Only [`validate`] can set this.
    pub normal_count_mismatch: bool,

    /// Echoed so `Display` can print the thresholds that produced the counts.
    pub config: ValidateConfig,
}

impl MeshReport {
    /// Malformed input: a bad index, a partial triangle, a non-finite position,
    /// or an attribute-length mismatch.
    ///
    /// Always a bug, whatever the field looks like. When this is true the
    /// derived metrics describe only the valid subset, and `Display` says so.
    #[must_use]
    pub fn has_structural_errors(&self) -> bool {
        self.out_of_range_indices > 0
            || self.trailing_indices > 0
            || self.repeated_index_triangles > 0
            || self.non_finite_positions > 0
            || self.normal_count_mismatch
    }

    /// A 2-manifold, possibly with boundary, consistently oriented.
    ///
    /// **The gate for open fields.**
    #[must_use]
    pub fn is_manifold(&self) -> bool {
        !self.has_structural_errors()
            && self.non_manifold_edges == 0
            && self.non_manifold_vertices == 0
            && self.inconsistently_oriented_edges == 0
    }

    /// [`is_manifold`](Self::is_manifold) and closed.
    ///
    /// **The gate for closed fields.** For any closed orientable surface
    /// `χ = 2 − 2g`, so `χ` must also be even; that is checked here rather than
    /// left to the caller, since it holds for every closed field regardless of
    /// whether the genus itself is known.
    #[must_use]
    pub fn is_closed(&self) -> bool {
        self.is_manifold() && self.boundary_edges == 0 && self.euler_characteristic % 2 == 0
    }

    /// Whether this mesh meets the gate its field and algorithm earn.
    ///
    /// **This is the rule, and it is not the caller's to guess.** "Valid" is not
    /// one thing: a closed solid must have no boundary, an open surface is
    /// *supposed* to, and a grid too coarse to resolve its field is permitted to
    /// be non-manifold without that being a mesher defect. Applying the wrong one
    /// makes a correct mesh read as broken — which is exactly what happened
    /// downstream before this was public (✗22).
    ///
    /// Every arm additionally requires
    /// [`has_structural_errors`](Self::has_structural_errors) to be false: a
    /// malformed mesh satisfies nothing, including the loosest gate.
    ///
    /// # Which artefact gets which
    ///
    /// - A **solid** — a closed proxy, a watertight body, a chunk-free extraction
    ///   of a closed field — takes [`SurfaceGate::Closed`].
    /// - A **surface** — an open field, a single chunk, a render mesh that is a
    ///   subset of some larger body — takes [`SurfaceGate::Manifold`]. Its open
    ///   edges are a **recorded number, not a failure**.
    /// - A solid on a grid that may not resolve it takes
    ///   [`SurfaceGate::ClosedAllowingUnresolvedTopology`].
    #[must_use]
    pub fn satisfies(&self, gate: SurfaceGate) -> bool {
        if self.has_structural_errors() {
            return false;
        }
        match gate {
            SurfaceGate::Closed => self.is_closed(),
            SurfaceGate::Manifold => self.is_manifold(),
            SurfaceGate::ClosedAllowingUnresolvedTopology => {
                self.boundary_edges == 0 && self.inconsistently_oriented_edges == 0
            }
        }
    }

    /// Total of every violation counter that gates correctness.
    ///
    /// [`degenerate_triangles`](Self::degenerate_triangles),
    /// [`duplicate_vertices`](Self::duplicate_vertices) and
    /// [`unreferenced_vertices`](Self::unreferenced_vertices) are excluded: they
    /// are recorded metrics, and a correct extractor produces non-zero counts of
    /// all three for perfectly ordinary reasons.
    #[must_use]
    pub fn violations(&self) -> u64 {
        self.non_manifold_edges
            + self.non_manifold_vertices
            + self.inconsistently_oriented_edges
            + self.repeated_index_triangles
            + self.out_of_range_indices
            + self.trailing_indices
            + self.non_finite_positions
            + u64::from(self.normal_count_mismatch)
    }

    /// The loud path, opt in.
    ///
    /// # Panics
    ///
    /// If the mesh fails [`is_closed`](Self::is_closed) (when `expect_closed`) or
    /// [`is_manifold`](Self::is_manifold) (otherwise). The message is the whole
    /// `Display` block, so the failure explains itself.
    pub fn panic_if_invalid(&self, expect_closed: bool) {
        let ok = if expect_closed {
            self.is_closed()
        } else {
            self.is_manifold()
        };
        assert!(ok, "{self}");
    }
}

impl fmt::Display for MeshReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let partial = self.faces_skipped > 0 || self.trailing_indices > 0;
        let rule = "  ------------------------------------------------------------";

        writeln!(f, "mesh report")?;
        writeln!(
            f,
            "  vertices                 {:8}  ({} referenced, {} unreferenced)",
            self.vertices, self.referenced_vertices, self.unreferenced_vertices
        )?;
        writeln!(f, "  edges                    {:8}", self.edges)?;
        if partial {
            writeln!(
                f,
                "  faces                    {:8}  ({} considered, {} skipped)",
                self.faces + self.faces_skipped,
                self.faces,
                self.faces_skipped
            )?;
            writeln!(
                f,
                "  euler characteristic     {:8}  (over the valid subset only)",
                self.euler_characteristic
            )?;
        } else {
            writeln!(f, "  faces                    {:8}", self.faces)?;
            writeln!(
                f,
                "  euler characteristic     {:8}",
                self.euler_characteristic
            )?;
        }
        writeln!(f, "  components               {:8}", self.components)?;
        writeln!(f, "  boundary loops           {:8}", self.boundary_loops)?;
        match self.genus {
            Some(g) => writeln!(f, "  genus                    {g:8}")?,
            None => writeln!(
                f,
                "  genus                           -  (not a single oriented manifold component)"
            )?,
        }

        writeln!(f, "{rule}")?;
        writeln!(
            f,
            "  non-manifold edges       {:8}",
            self.non_manifold_edges
        )?;
        writeln!(
            f,
            "  non-manifold vertices    {:8}",
            self.non_manifold_vertices
        )?;
        writeln!(f, "  boundary edges           {:8}", self.boundary_edges)?;
        writeln!(
            f,
            "  inconsistently oriented  {:8}",
            self.inconsistently_oriented_edges
        )?;
        writeln!(
            f,
            "  degenerate triangles     {:8}  (area <= {:e} * h^2, h = {})",
            self.degenerate_triangles, self.config.area_epsilon_rel, self.config.cell_size
        )?;
        writeln!(
            f,
            "  repeated-index triangles {:8}",
            self.repeated_index_triangles
        )?;
        writeln!(
            f,
            "  duplicate vertices       {:8}  (within {:e})",
            self.duplicate_vertices, self.config.weld_epsilon
        )?;
        writeln!(
            f,
            "  out-of-range indices     {:8}",
            self.out_of_range_indices
        )?;
        writeln!(f, "  trailing indices         {:8}", self.trailing_indices)?;
        writeln!(
            f,
            "  non-finite positions     {:8}",
            self.non_finite_positions
        )?;
        if self.normal_count_mismatch {
            writeln!(f, "  normal count mismatch         yes")?;
        }

        writeln!(f, "{rule}")?;
        if self.has_structural_errors() {
            writeln!(
                f,
                "  !! STRUCTURAL ERRORS - derived metrics cover the valid subset only"
            )?;
        }
        if self.is_closed() {
            write!(f, "  MANIFOLD, CLOSED")
        } else if self.is_manifold() {
            write!(f, "  MANIFOLD, WITH BOUNDARY")
        } else {
            write!(f, "  INVALID: {} violations", self.violations())
        }
    }
}

/// Disjoint-set union over `u32` labels, with path halving and union by size.
///
/// Deterministic: the result depends only on which unions were requested, and
/// the callers below request them in sorted order.
struct Dsu {
    parent: Vec<u32>,
    size: Vec<u32>,
}

impl Dsu {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n as u32).collect(),
            size: alloc::vec![1; n],
        }
    }

    /// Reuse the allocation for the next vertex's link, so the per-vertex walk
    /// does not allocate.
    fn reset(&mut self, n: usize) {
        self.parent.clear();
        self.parent.extend(0..n as u32);
        self.size.clear();
        self.size.resize(n, 1);
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            let grand = self.parent[self.parent[x as usize] as usize];
            self.parent[x as usize] = grand;
            x = grand;
        }
        x
    }

    fn union(&mut self, a: u32, b: u32) {
        let (mut ra, mut rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        if self.size[ra as usize] < self.size[rb as usize] {
            core::mem::swap(&mut ra, &mut rb);
        }
        self.parent[rb as usize] = ra;
        self.size[ra as usize] += self.size[rb as usize];
    }
}

/// Which features of a mesh offend, as vertex indices rather than as counts.
///
/// Produced by [`validate_features`]. Every list is in the same order the
/// validator visits them -- edges ascending by `[lo, hi]`, vertices ascending --
/// so this is a pure function of the mesh's values, exactly like [`MeshReport`].
///
/// Indices, not positions: the caller already has the position buffer, and
/// handing back copies would mean a consumer could hold geometry that no longer
/// matches the mesh it came from.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NonManifoldFeatures {
    /// Edges used by **three or more** faces, as `[lo, hi]` vertex indices.
    /// Length equals [`MeshReport::non_manifold_edges`].
    pub edges: Vec<[u32; 2]>,
    /// Vertices whose incident-face link is not a single fan -- the bowtie that
    /// only the link walk catches. Length equals
    /// [`MeshReport::non_manifold_vertices`].
    pub vertices: Vec<u32>,
    /// Edges used by **exactly one** face: the mesh's boundary. Not a defect on
    /// an open field, which is why it is kept separate from `edges`. Length
    /// equals [`MeshReport::boundary_edges`].
    pub boundary_edges: Vec<[u32; 2]>,
    /// Edges whose two faces traverse them the *same* way round. Length equals
    /// [`MeshReport::inconsistently_oriented_edges`].
    pub inconsistently_oriented_edges: Vec<[u32; 2]>,
}

/// Validate an indexed triangle mesh.
///
/// Never panics on a malformed mesh, never returns an error, and never
/// short-circuits: every pass runs to completion over the valid subset, so a
/// mesh with one bad index still yields a complete manifoldness report.
///
/// Takes slices rather than a [`MeshBuffer`] because the callers that matter do
/// not have one: a golden-hash fixture hashes raw arrays, a property test
/// generates them, a GPU path reads them back from a buffer, and an engine
/// wrapper holds them as attribute arrays.
///
/// Cannot fail: a [`ValidateConfig`] can only be built through its checked
/// constructor, so there is no invalid threshold to guard against here.
#[must_use]
pub fn validate_indexed<R: Real>(
    positions: &[[R; 3]],
    indices: &[u32],
    cfg: &ValidateConfig,
) -> MeshReport {
    validate_features(positions, indices, cfg).0
}

/// [`validate_indexed`], plus **which** features offend rather than only how
/// many.
///
/// A count cannot be drawn. `E-111` renders non-manifold edges as red lines and
/// non-manifold vertices as red spheres, and for that to mean anything the
/// geometry on screen and the number in the HUD have to come from the same pass
/// -- otherwise the picture and the caption can disagree and nobody can tell
/// which is wrong.
///
/// So this is not a second implementation. It **is** the implementation;
/// [`validate_indexed`] calls it and drops the second half. The lists are
/// collected unconditionally rather than behind a flag, which costs nothing on
/// the common path: a clean mesh leaves every one of them empty, and an empty
/// [`Vec`] does not allocate.
///
/// The same argument as [`self_intersections`]: report the offenders, not just
/// the tally, because the caller is the one who knows what to do with them.
pub fn validate_features<R: Real>(
    positions: &[[R; 3]],
    indices: &[u32],
    cfg: &ValidateConfig,
) -> (MeshReport, NonManifoldFeatures) {
    let mut features = NonManifoldFeatures::default();
    let vertex_count = positions.len();
    let mut report = MeshReport {
        vertices: vertex_count as u64,
        referenced_vertices: 0,
        edges: 0,
        faces: 0,
        faces_skipped: 0,
        euler_characteristic: 0,
        components: 0,
        boundary_loops: 0,
        genus: None,
        non_manifold_edges: 0,
        non_manifold_vertices: 0,
        boundary_edges: 0,
        inconsistently_oriented_edges: 0,
        degenerate_triangles: 0,
        repeated_index_triangles: 0,
        duplicate_vertices: 0,
        weld_buckets: 0,
        unreferenced_vertices: 0,
        out_of_range_indices: 0,
        trailing_indices: (indices.len() % 3) as u64,
        non_finite_positions: 0,
        normal_count_mismatch: false,
        config: *cfg,
    };

    for p in positions {
        if !(p[0].is_finite() && p[1].is_finite() && p[2].is_finite()) {
            report.non_finite_positions += 1;
        }
    }

    // ── pass 0: which faces are usable at all ───────────────────────────────
    //
    // A face with an out-of-range index cannot be dereferenced, and one with a
    // repeated index has no well-defined edge set, so including it would corrupt
    // the edge count. Both are skipped and counted; nothing short-circuits.
    let whole = indices.len() - indices.len() % 3;
    let mut faces: Vec<[u32; 3]> = Vec::with_capacity(whole / 3);
    for tri in indices[..whole].chunks_exact(3) {
        let mut out_of_range = 0u64;
        for &i in tri {
            if i as usize >= vertex_count {
                out_of_range += 1;
            }
        }
        let repeated = tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2];
        if repeated {
            report.repeated_index_triangles += 1;
        }
        if out_of_range > 0 || repeated {
            report.out_of_range_indices += out_of_range;
            report.faces_skipped += 1;
            continue;
        }
        faces.push([tri[0], tri[1], tri[2]]);
    }
    report.faces = faces.len() as u64;

    // ── vertex references ───────────────────────────────────────────────────
    let mut referenced = alloc::vec![false; vertex_count];
    for f in &faces {
        for &v in f {
            referenced[v as usize] = true;
        }
    }
    report.referenced_vertices = referenced.iter().filter(|r| **r).count() as u64;
    report.unreferenced_vertices = report.vertices - report.referenced_vertices;

    // ── pass 1: edges ───────────────────────────────────────────────────────
    //
    // One entry per directed half-edge, canonicalised to (min, max) with the
    // original direction kept as a flag. The face index is part of the key, so
    // no two entries compare equal and the sort is a total order on values.
    let mut half_edges: Vec<(u32, u32, u32, u8)> = Vec::with_capacity(faces.len() * 3);
    for (fi, f) in faces.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (f[k], f[(k + 1) % 3]);
            let (lo, hi, forward) = if a < b { (a, b, 1) } else { (b, a, 0) };
            half_edges.push((lo, hi, fi as u32, forward));
        }
    }
    half_edges.sort_unstable();

    let mut face_dsu = Dsu::new(faces.len());
    let mut boundary_dsu = Dsu::new(vertex_count);
    let mut on_boundary = alloc::vec![false; vertex_count];

    let mut start = 0usize;
    while start < half_edges.len() {
        let (lo, hi, _, _) = half_edges[start];
        let mut end = start + 1;
        while end < half_edges.len() && half_edges[end].0 == lo && half_edges[end].1 == hi {
            end += 1;
        }
        let run = &half_edges[start..end];
        report.edges += 1;

        match run.len() {
            1 => {
                report.boundary_edges += 1;
                features.boundary_edges.push([lo, hi]);
                on_boundary[lo as usize] = true;
                on_boundary[hi as usize] = true;
                boundary_dsu.union(lo, hi);
            }
            2 => {
                if run[0].3 == run[1].3 {
                    report.inconsistently_oriented_edges += 1;
                    features.inconsistently_oriented_edges.push([lo, hi]);
                }
            }
            _ => {
                report.non_manifold_edges += 1;
                features.edges.push([lo, hi]);
            }
        }

        // Faces meeting along this edge are in the same component regardless of
        // how many of them there are.
        for w in run.windows(2) {
            face_dsu.union(w[0].2, w[1].2);
        }

        start = end;
    }

    // ── components and boundary loops ───────────────────────────────────────
    {
        let mut seen: Vec<u32> = (0..faces.len() as u32).map(|f| face_dsu.find(f)).collect();
        seen.sort_unstable();
        seen.dedup();
        report.components = seen.len() as u64;

        let mut loops: Vec<u32> = (0..vertex_count as u32)
            .filter(|v| on_boundary[*v as usize])
            .map(|v| boundary_dsu.find(v))
            .collect();
        loops.sort_unstable();
        loops.dedup();
        report.boundary_loops = loops.len() as u64;
    }

    // ── pass 2: vertex links ────────────────────────────────────────────────
    //
    // A vertex is non-manifold when the faces around it do not form one
    // connected fan. The cheap test -- "incident faces equals incident edges" --
    // reports a bowtie as clean: two cones sharing an apex have 2k faces and 2k
    // wing edges, every edge has exactly two faces, and chi can come out right.
    // Nothing else in this module would catch it, which is what pays for the
    // walk.
    {
        let mut vertex_faces: Vec<(u32, u32)> = Vec::with_capacity(faces.len() * 3);
        for (fi, f) in faces.iter().enumerate() {
            for &v in f {
                vertex_faces.push((v, fi as u32));
            }
        }
        vertex_faces.sort_unstable();

        let mut wings: Vec<(u32, u32)> = Vec::new();
        let mut link = Dsu::new(0);

        let mut start = 0usize;
        while start < vertex_faces.len() {
            let v = vertex_faces[start].0;
            let mut end = start + 1;
            while end < vertex_faces.len() && vertex_faces[end].0 == v {
                end += 1;
            }
            let degree = end - start;

            // The two "wing" vertices of each incident face, tagged with that
            // face's position in this vertex's own list.
            wings.clear();
            for (local, &(_, fi)) in vertex_faces[start..end].iter().enumerate() {
                for &w in &faces[fi as usize] {
                    if w != v {
                        wings.push((w, local as u32));
                    }
                }
            }
            wings.sort_unstable();

            link.reset(degree);
            let mut branching = false;
            let mut i = 0usize;
            while i < wings.len() {
                let mut j = i + 1;
                while j < wings.len() && wings[j].0 == wings[i].0 {
                    j += 1;
                }
                if j - i >= 3 {
                    // Three faces share one edge at this vertex: an umbrella,
                    // not a fan.
                    branching = true;
                }
                for k in i + 1..j {
                    link.union(wings[i].1, wings[k].1);
                }
                i = j;
            }

            let mut roots: Vec<u32> = (0..degree as u32).map(|d| link.find(d)).collect();
            roots.sort_unstable();
            roots.dedup();
            if branching || roots.len() > 1 {
                report.non_manifold_vertices += 1;
                features.vertices.push(v);
            }

            start = end;
        }
    }

    // ── degenerate triangles ────────────────────────────────────────────────
    //
    // Compared squared, so there is no square root and this works unchanged in
    // no_std. A NaN position makes the comparison false rather than true, which
    // is why `non_finite_positions` is counted separately.
    {
        let two_area_limit =
            R::from_f64(2.0 * cfg.area_epsilon_rel * cfg.cell_size * cfg.cell_size);
        let limit_sq = two_area_limit * two_area_limit;
        for f in &faces {
            let a = positions[f[0] as usize];
            let b = positions[f[1] as usize];
            let c = positions[f[2] as usize];
            let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let cross = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let len_sq = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
            if len_sq <= limit_sq {
                report.degenerate_triangles += 1;
            }
        }
    }

    // ── duplicate vertices ──────────────────────────────────────────────────
    //
    // Bucket onto the same integer lattice `crate::weld` welds on, sort, then
    // probe the 27 neighbouring cells. Quantising sidesteps `f32: !Ord`
    // entirely: the key is `[i64; 3]`, which has a genuine total order. `-0.0`
    // and `+0.0` land in the same bucket. A NaN coordinate saturates to zero and
    // then fails the distance comparison, so it is never reported as a
    // duplicate.
    //
    // **The lattice is `weld::Lattice`, not a copy of it.** If the two disagreed
    // about which cells are probed, this count would describe a different
    // neighbourhood than the weld it is reporting on. What they do differ in is
    // the question asked: this asks whether *any* earlier vertex is within
    // epsilon, the welder asks for the lowest-indexed *kept* one, so this is an
    // upper bound on what a weld removes rather than a prediction of it — see
    // `weld::tests::the_validator_bounds_the_weld_rather_than_predicting_it`,
    // which measures both.
    {
        let eps = R::from_f64(cfg.weld_epsilon);
        let eps_sq = eps * eps;
        let lattice = crate::weld::Lattice::new(positions, eps);
        let key_of = |p: [R; 3]| lattice.key_of(p);

        let mut cells: Vec<([i64; 3], u32)> = (0..vertex_count)
            .map(|i| (key_of(positions[i]), i as u32))
            .collect();
        cells.sort_unstable();

        {
            let mut distinct = 0u64;
            let mut previous: Option<[i64; 3]> = None;
            for (k, _) in &cells {
                if previous != Some(*k) {
                    distinct += 1;
                    previous = Some(*k);
                }
            }
            report.weld_buckets = distinct;
        }

        for v in 0..vertex_count as u32 {
            let p = positions[v as usize];
            let base = key_of(p);
            let mut found = false;
            'probe: for dz in -1..=1i64 {
                for dy in -1..=1i64 {
                    for dx in -1..=1i64 {
                        let key = [base[0] + dx, base[1] + dy, base[2] + dz];
                        let lo = cells.partition_point(|(k, _)| *k < key);
                        for &(k, u) in &cells[lo..] {
                            if k != key {
                                break;
                            }
                            if u >= v {
                                continue;
                            }
                            let q = positions[u as usize];
                            let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
                            if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= eps_sq {
                                found = true;
                                break 'probe;
                            }
                        }
                    }
                }
            }
            if found {
                report.duplicate_vertices += 1;
            }
        }
    }

    // ── derived ─────────────────────────────────────────────────────────────
    report.euler_characteristic =
        report.referenced_vertices as i64 - report.edges as i64 + report.faces as i64;

    if report.components == 1
        && report.faces_skipped == 0
        && report.non_manifold_edges == 0
        && report.non_manifold_vertices == 0
        && report.inconsistently_oriented_edges == 0
    {
        let doubled = 2 - report.euler_characteristic - report.boundary_loops as i64;
        if doubled % 2 == 0 {
            report.genus = Some(doubled / 2);
        }
    }

    (report, features)
}

/// [`validate_indexed`] for the default sink, plus the one check that needs the
/// normals: that there is exactly one per vertex.
///
/// A plain forwarder — same code, same behaviour, no second path.
///
/// # Panics
///
/// As [`validate_indexed`].
#[must_use]
pub fn validate<R: Real>(mesh: &MeshBuffer<R>, cfg: &ValidateConfig) -> MeshReport {
    let mut report = validate_indexed(&mesh.positions, &mesh.indices, cfg);
    report.normal_count_mismatch = mesh.normals.len() != mesh.positions.len();
    report
}

/// Which validity gate a mesh is held to.
///
/// Deliberately an enum rather than a `bool`, and deliberately three cases: a
/// blanket gate is unsatisfiable for at least one field *and* at least one
/// algorithm, so the caller has to name which one applies and why.
///
/// # Why this is public (T-023, ✗22)
///
/// [`MeshReport`] offers three predicates and, until T-023, no reachable
/// statement of **which one applies to what**. That rule existed and was
/// correct, but lived in a `#[cfg(test)]` module — so it was compiled out of
/// every shipped build, and consumers re-derived it. They got it wrong in the
/// obvious way: calling [`is_closed`](MeshReport::is_closed) on a *render* mesh,
/// which was never a solid, and reading the failure as a mesher defect.
///
/// This type is the data half of the rule and
/// [`MeshReport::satisfies`] is the policy half. Both ship, because a tag with
/// no policy makes every consumer reinvent it, and a policy with no tag makes a
/// different one impossible.
///
/// # Choosing one
///
/// **Not from intuition.** The gate is a property of the *field and the
/// algorithm*, not of the caller's expectations — a closed field on a grid that
/// resolves it earns [`Closed`](Self::Closed), an open field earns
/// [`Manifold`](Self::Manifold), and a grid that may not resolve the field's
/// topology earns
/// [`ClosedAllowingUnresolvedTopology`](Self::ClosedAllowingUnresolvedTopology).
/// A field that knows whether it is closed in its own domain should say so; see
/// `ReferenceField::closed_in_domain` for how this crate's own fields do it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum SurfaceGate {
    /// A closed, oriented 2-manifold. **The gate for a closed field on a grid
    /// that resolves it** — the seven reference fields at their own resolutions,
    /// meshed by Marching Cubes.
    ///
    /// Note the condition is on the *pair*, not on the algorithm. Marching Cubes
    /// was believed to earn this unconditionally, by placing vertices on grid
    /// edges rather than one per cell. ✗15 shows it does not: refine far enough
    /// and it does, but a surface that pinches inside one cell defeats it.
    Closed,
    /// A 2-manifold, possibly with boundary. **The gate for an open field**, and
    /// for any single chunk once G-001 lands.
    Manifold,
    /// Closed, correctly oriented and wholly inside the grid, but *permitted* to
    /// be non-manifold.
    ///
    /// **The gate for a grid that may not resolve the field's topology.** Note
    /// what this is *not* keyed on: it started life as "the gate for
    /// one-vertex-per-cell methods", and that was wrong twice over.
    ///
    /// Two distinct mechanisms land here, which is why the name is about the
    /// grid rather than the algorithm:
    ///
    /// - **Surface Nets and plain Dual Contouring** place one vertex per cell, so
    ///   two sheets passing through one cell must share it. The literature calls
    ///   this DC's *"actual structural defect"*; A-010 fixes it architecturally
    ///   by vertex splitting. M-4 measured it on `gyroid` and `fbm_terrain` and
    ///   read it as a high-genus/open-field effect; M-15 corrected that — a
    ///   generated **convex body** does it too, so it is about resolution.
    /// - **Marching Cubes**, which was believed unconditionally manifold, does it
    ///   too where the surface *pinches* inside a single cell: the shared grid
    ///   edge ends up carrying four faces. See ✗15 and
    ///   `an_under_resolved_pinch_makes_marching_cubes_non_manifold`, which pins
    ///   the exact counts at `h = 2/3` and their disappearance by `h = 1/2`.
    ///
    /// The strict [`Closed`](Self::Closed) claim is still tested where it is
    /// actually true — the seven reference fields at their own resolutions, in
    /// `mc/tests.rs`. It is only the *generated* fields, which are adversarial by
    /// construction and go as coarse as `h = 2/3`, that need this.
    ///
    /// What is still asserted is everything unrelated to unresolved topology: no
    /// structural errors, no boundary (the surface did not leave the grid), and
    /// consistent winding.
    ///
    /// **The even-`χ` parity check is deliberately *not* asserted here**, and
    /// that is not an oversight. `χ = 2 − 2g` — hence `χ` even — holds for a
    /// closed *orientable manifold*, so parity is a corollary of manifoldness
    /// rather than an independent check. Waiving manifoldness and keeping the
    /// parity check is incoherent, and measurably so: Surface Nets on a generated
    /// convex body produces `χ = 1` with one non-manifold edge and zero boundary
    /// edges. A-010 is where this becomes assertable again.
    ClosedAllowingUnresolvedTopology,
}

#[cfg(test)]
mod tests;
