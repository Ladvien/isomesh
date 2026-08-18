//! Does the mesh separate what the field separates?
//!
//! Ticket: R-024, hypothesis P-21. The instrument is **not ours** — see the
//! *Provenance* section below.
//!
//! # The claim nothing else in this crate checks
//!
//! Everything else here validates a mesh against **itself** — manifoldness,
//! orientation, Euler characteristic, self-intersection — or against the field's
//! **geometry**, in [`accuracy`](super::accuracy). None of it asks whether the
//! mesh partitions *space* the way the field's sign does. Those are different
//! claims, and neither implies the other: a mesh can be closed, manifold,
//! correctly oriented and Hausdorff-close while sealing a passage the field
//! leaves open, or opening one the field seals.
//!
//! For a game that is the whole question. *Is this cave sealed? Did I just break
//! through?* are questions about the connected components of the air region, and
//! they are asked of the **mesh** — the collider — while the authority is the
//! **field**. If the two disagree, a player walks through a wall or into one
//! that is not there.
//!
//! # The probe is a grid edge, and that is what makes it work on open fields
//!
//! An extractor never sees the field anywhere but at the samples, so the air
//! sublevel set *as the extractor could know it* is the 6-connected graph on
//! samples with a non-negative value. Probing along the **grid edges** keeps
//! both sides of the comparison on the same lattice, and keeps the test
//! **local** — which is what lets it run on `gyroid` and `fbm_terrain`, whose
//! surfaces leave the domain and for which any global inside/outside test is
//! undefined at the boundary.
//!
//! For each 6-adjacent pair of samples the two sides say:
//!
//! - **the field**: the pair straddles the surface iff their signs differ;
//! - **the mesh**: the pair is separated iff the segment between them is crossed
//!   an **odd** number of times.
//!
//! Disagreements come in two kinds and they are not the same defect.
//! [`unsealed_walls`](SealingReport::unsealed_walls) is the field saying
//! *separated* and the mesh saying *connected* — a hole. [`spurious_walls`]
//! (SealingReport::spurious_walls) is the reverse — a membrane across open air.
//!
//! # Counting *points*, not triangle hits, and why that is load-bearing
//!
//! A Marching Cubes vertex lies **exactly on** the segment being probed — that
//! is what Marching Cubes *is*, the root of the interpolant along a grid edge —
//! and every triangle in the fan around it contains that point. Counting
//! triangle hits would therefore report 4 where the surface crosses once, and
//! **every sign-changing edge in the crate's flagship extractor would read
//! even**: a harness that got this wrong would report total failure while
//! measuring nothing but its own convention.
//!
//! So crossings are counted as **distinct points**, deduplicated along the
//! segment at the crate's own weld tolerance, and
//! [`merged_crossings`](SealingReport::merged_crossings) reports how many raw
//! hits collapsed. On a dual method that number is near zero — its vertices are
//! in cell interiors and its quads cross grid edges transversally — so the
//! column doubles as evidence that the mechanism is real and not a fudge.
//!
//! # Provenance — the test is Wojtan et al.'s, the audit is not
//!
//! Wojtan, Thürey, Gross & Turk, *Physics-inspired topology changes for thin
//! fluid features*, SIGGRAPH 2010 (`10.1145/1778765.1778787`) define this
//! exactly, as their **complex edge test**: *"we determine the complexity of a
//! cell edge by counting the number and orientation of its intersections with
//! triangles in the surface mesh and comparing the result with a single line
//! segment (0-sphere)."* They run it over every edge of the signed-distance grid
//! and bound its coverage — it *"is guaranteed to identify any topological flaws
//! that are well-resolved in at least two dimensions"*, with face and cell tests
//! above it for thinner defects.
//!
//! What they use it for is a different question: their mesh is *advected by a
//! velocity field* and the distance field rebuilt around it, so disagreement is
//! expected and is the thing being repaired. Running it as a **correctness audit
//! of extraction itself** — is a mesh freshly extracted from a field faithful to
//! that field's partition — is what R-024 adds. See V-37.
//!
//! # Coverage, stated rather than implied
//!
//! Inherited from the same source: this sees any disagreement that is resolved
//! in at least two dimensions. A defect thinner than a cell in all three — a
//! spindle or a void wholly inside one cell, touching no grid edge — is
//! invisible to it, and would need the face and cell tests Wojtan et al. put
//! above this one. Neither is implemented here.
//!
//! # Cost
//!
//! `O(T + E·k)` for `T` triangles, `E ≈ 3n³` grid edges and `k` candidate
//! triangles per edge, plus `O(n³)` for the two union-finds. Triangles are
//! binned into extraction cells once, in CSR form. This measures a mesh; it does
//! not produce one, and it is not on any hot path.

use alloc::vec::Vec;
use core::fmt;

use crate::cube::is_inside;
use crate::{Real, Sdf, Shape3, vec3};

/// Tolerance for calling two crossings the same point, as a fraction of the
/// probe segment.
///
/// The segment is one cell long, so this is
/// [`ValidateConfig::WELD_EPSILON_REL`](super::ValidateConfig::WELD_EPSILON_REL)
/// in the crate's usual units: two hits closer than a ten-thousandth of a cell
/// are the same place, which is the same judgement
/// [`weld`](crate::weld) makes about two vertices.
const SAME_POINT_REL: f64 = super::ValidateConfig::WELD_EPSILON_REL;

/// Relative tolerance for calling a probe parallel to a triangle's plane.
///
/// Scaled by `|d| · |e₁| · |e₂|`, so it is a test on the *sine* of the angle
/// between the segment and the plane rather than on a raw determinant, and does
/// not change meaning with the grid spacing.
const PARALLEL_REL: f64 = 1e-12;

/// How the mesh's partition of the sample lattice compares with the field's.
///
/// Produced by [`sealing`]. Every field is a count; nothing here judges, in
/// keeping with the rest of [`validate`](super) — [`agrees`](Self::agrees) is
/// the one derived predicate and it is opt-in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SealingReport {
    /// Grid samples, `sx · sy · sz`.
    pub samples: u64,
    /// Samples the field calls air — value `>= 0`, the complement of
    /// `cube::is_inside`.
    pub air_samples: u64,
    /// 6-adjacent sample pairs probed. Every grid edge, once.
    pub probes: u64,

    /// Probes whose endpoints straddle the surface, by the field's sign alone.
    ///
    /// What the mesh is *supposed* to separate.
    pub field_walls: u64,
    /// Probes the mesh separates: an **odd** number of distinct crossings.
    pub mesh_walls: u64,

    /// **A hole.** The field separates the pair and the mesh does not.
    ///
    /// Air leaks into solid across a boundary the field says is closed. For a
    /// collider this is the defect a player falls through.
    pub unsealed_walls: u64,
    /// Of the holes, how many touch a face of the sampled domain.
    ///
    /// **This separates two different defects that look identical in the
    /// total.** A dual method emits one quad per sign-changing grid edge, and
    /// that quad needs all **four** cells around the edge; on a domain face only
    /// one or two exist, so no quad is emitted and the wall is left open. That
    /// is a property of where the grid stops, not of the field — and for a
    /// chunked world it is the chunk seam. A hole strictly inside the domain is
    /// a different animal and has nothing to do with clipping.
    ///
    /// Measured on `fbm_terrain`, the one reference field whose surface leaves
    /// through the sides: **92 of 92** at 17³, and the same at 25³ and 33³, for
    /// all three dual methods.
    pub unsealed_on_domain_face: u64,
    /// **A membrane.** The mesh separates the pair and the field does not.
    ///
    /// A wall across open air, or across solid rock. Less dramatic than a hole
    /// and the same class of error: the mesh is asserting a boundary the field
    /// does not have.
    pub spurious_walls: u64,

    /// 6-connected components of the air sublevel set, from the field alone.
    pub field_air_components: u64,
    /// The same samples' components once the mesh cuts the adjacency.
    ///
    /// Cannot be **less** than [`field_air_components`](Self::field_air_components)
    /// unless [`unsealed_walls`](Self::unsealed_walls) is non-zero, since the
    /// mesh can only ever remove air-to-air adjacency it did not add.
    pub mesh_air_components: u64,
    /// Components of the **whole** lattice under mesh-cut adjacency, both phases.
    pub mesh_regions: u64,
    /// Of those, the ones holding samples of **both** signs.
    ///
    /// A region the mesh leaves connected while the field puts solid and air in
    /// it: the same defect as [`unsealed_walls`](Self::unsealed_walls) seen at
    /// component scale rather than at edge scale, and the one that survives when
    /// several holes conspire.
    pub mixed_regions: u64,

    /// Raw triangle hits that collapsed into an already-counted point.
    ///
    /// **Read this before believing a zero anywhere else.** It is large for
    /// Marching Cubes by construction — its vertex is *on* the probe and its
    /// whole fan contains it — and near zero for a dual method. See the module
    /// docs.
    pub merged_crossings: u64,
    /// Probes skipped against a triangle whose plane contains them.
    ///
    /// A crossing count is undefined there. Non-generic, reported so that a
    /// non-zero cannot hide inside an agreement.
    pub coplanar_probes: u64,
    /// Triangles with no area, excluded before binning.
    ///
    /// A degenerate triangle separates nothing, and including it is not
    /// conservative but *wrong*: its normal is a cancellation residue, so it
    /// answers "parallel" to every probe and lands in
    /// [`coplanar_probes`](Self::coplanar_probes) once per probe it is binned
    /// near. Measured before the exclusion existed: Marching Tetrahedra's **36**
    /// slivers on `sphere` at 17³ produced **6,624** coplanar events between
    /// them, and every other extractor produced none.
    pub degenerate_triangles: u64,

    /// Samples whose field value is **exactly** the iso value.
    ///
    /// The surface passes through them, so they are in neither open phase and
    /// `is_inside`'s tie-break — not the geometry — decides which side they
    /// count as. See [`degenerate_probes`](Self::degenerate_probes).
    pub boundary_samples: u64,
    /// Probes set aside because a crossing landed **on** an endpoint.
    ///
    /// The field-side test for a boundary sample is `value == 0` **exactly**,
    /// and an exact test undercounts this phenomenon — the same lesson A-002i
    /// records for singular faces. Measured: on `sphere` at 25³ a sample carries
    /// `−1.11e−16`, which is solid by `is_inside` and is the surface to within
    /// one ulp. A primal method puts its vertex at `t = a/(a − b)`, so `a ≈ 0`
    /// lands it essentially **at** the sample, and its fan then touches every
    /// same-sign probe there at `t = 1`. One rounding away from the exact case
    /// and geometrically identical to it.
    ///
    /// So the degeneracy is detected from the **mesh** side as well: a crossing
    /// within the same-point tolerance of `t = 0` or `t = 1` means the surface
    /// passes through a sample, and the parity of a path that *starts* on the
    /// surface has no defined value.
    pub endpoint_crossings: u64,
    /// Probes set aside as undecidable, for either reason.
    ///
    /// The field-side reason is a [boundary sample](Self::boundary_samples); the
    /// mesh-side one is an [endpoint crossing](Self::endpoint_crossings). Both
    /// are the same geometry — the surface passes through a sample — seen from
    /// the two sides, and the two counts overlap.
    ///
    /// **Not a defect, and not a convenience either.** A crossing count answers
    /// "does the path from A to B meet the surface", and it has no answer when
    /// the surface passes through A. Both directions fail, symmetrically: toward
    /// air the surface *touches* the segment at its endpoint and a parity count
    /// reads that as a crossing; toward solid it genuinely separates and the
    /// crossing sits at the closed endpoint where an open-interval count cannot
    /// see it.
    ///
    /// It is a **primal** phenomenon, measured. A method that places its vertex
    /// on a grid edge at `t = a/(a − b)` puts it exactly at the sample when
    /// `a = 0`, and its whole triangle fan then touches all five same-sign edges
    /// there. A dual method never places a vertex on a grid edge and never
    /// reaches this. Excluded from every wall column and from both graphs, so
    /// each side loses the same edges and any surviving difference is a real
    /// one.
    pub degenerate_probes: u64,
}

impl SealingReport {
    /// The mesh partitions the lattice exactly as the field's sign does.
    ///
    /// Every probe agrees, no region mixes phases, and the two air-component
    /// counts match. The last is implied by the first — equal edge sets give
    /// equal components — and is checked anyway, because a derived predicate
    /// that restates its own inputs catches a harness bug that a tighter one
    /// would not.
    #[must_use]
    pub const fn agrees(&self) -> bool {
        self.unsealed_walls == 0
            && self.spurious_walls == 0
            && self.mixed_regions == 0
            && self.field_air_components == self.mesh_air_components
    }
}

impl fmt::Display for SealingReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "sealing: {} probes, {} field walls, {} mesh walls",
            self.probes, self.field_walls, self.mesh_walls
        )?;
        writeln!(
            f,
            "  holes (unsealed) {} ({} on a domain face)   membranes (spurious) {}   \
             mixed regions {}",
            self.unsealed_walls,
            self.unsealed_on_domain_face,
            self.spurious_walls,
            self.mixed_regions
        )?;
        writeln!(
            f,
            "  air components: field {} / mesh {}   ({} regions over {} samples)",
            self.field_air_components, self.mesh_air_components, self.mesh_regions, self.samples
        )?;
        writeln!(
            f,
            "  merged crossings {}   coplanar probes {}   zero-area triangles {}",
            self.merged_crossings, self.coplanar_probes, self.degenerate_triangles
        )?;
        write!(
            f,
            "  {} samples on the surface, {} endpoint crossings, {} probes set aside   -> {}",
            self.boundary_samples,
            self.endpoint_crossings,
            self.degenerate_probes,
            if self.agrees() { "SEALED" } else { "DISAGREES" }
        )
    }
}

/// Compare the mesh's partition of the sample lattice with the field's.
///
/// The grid arguments are the ones handed to
/// [`Extractor::extract_into`](crate::extractor::Extractor::extract_into), so a
/// caller who has just extracted passes the same four values back. `positions`
/// and `indices` are that extraction's output; the mesh is **not** required to
/// be welded, closed or manifold, and nothing here assumes it is.
///
/// Sample positions are computed as `origin + cell_size · index`, character for
/// character the expression every extractor uses, so the classification this
/// compares against is the one the extractor saw.
///
/// # Sign convention
///
/// Negative is inside the solid, matching [`Sdf::sample`] and
/// `cube::is_inside`. **Air is `value >= 0`**, so a sample of exactly zero is
/// air — the same tie-break the case tables use, which is what keeps a
/// zero-valued corner from being classified one way by the extractor and the
/// other way here.
///
/// # Empty input
///
/// A mesh with no triangles is not an error: every probe reports zero crossings,
/// so the report says the mesh separates nothing and
/// [`unsealed_walls`](SealingReport::unsealed_walls) equals
/// [`field_walls`](SealingReport::field_walls). That is the honest reading of an
/// empty mesh and is more useful than a `Result`.
#[must_use]
pub fn sealing<R, S>(
    field: &S,
    shape: &impl Shape3,
    origin: [R; 3],
    cell_size: R,
    positions: &[[R; 3]],
    indices: &[u32],
) -> SealingReport
where
    R: Real,
    S: Sdf<Scalar = R>,
{
    let size = shape.size();
    let count = shape.element_count();

    let at = |i: [u32; 3]| {
        [
            origin[0] + cell_size * R::from_f64(f64::from(i[0])),
            origin[1] + cell_size * R::from_f64(f64::from(i[1])),
            origin[2] + cell_size * R::from_f64(f64::from(i[2])),
        ]
    };

    // ── the field's own classification, on the extractor's sample positions ──
    //
    // `on_surface` is the exact tie, and the exactness is the point: a sample
    // whose value is exactly the iso value is in neither open phase, so
    // `is_inside`'s convention rather than the geometry decides its side. See
    // `degenerate_probes`.
    let mut air = Vec::with_capacity(count);
    let mut on_surface = Vec::with_capacity(count);
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let v = field.sample(at([x, y, z]));
                air.push(!is_inside(v));
                #[allow(
                    clippy::float_cmp,
                    reason = "the exact tie is the condition being detected"
                )]
                on_surface.push(v == R::ZERO);
            }
        }
    }

    let bins = CellBins::build(size, origin, cell_size, positions, indices);

    let mut field_uf = UnionFind::new(count);
    let mut mesh_uf = UnionFind::new(count);
    let mut air_uf = UnionFind::new(count);

    let mut report = SealingReport {
        samples: count as u64,
        air_samples: air.iter().filter(|a| **a).count() as u64,
        probes: 0,
        field_walls: 0,
        mesh_walls: 0,
        unsealed_walls: 0,
        unsealed_on_domain_face: 0,
        spurious_walls: 0,
        field_air_components: 0,
        mesh_air_components: 0,
        mesh_regions: 0,
        mixed_regions: 0,
        merged_crossings: 0,
        coplanar_probes: 0,
        degenerate_triangles: bins.degenerate,
        boundary_samples: on_surface.iter().filter(|b| **b).count() as u64,
        endpoint_crossings: 0,
        degenerate_probes: 0,
    };

    let mut probe = Probe::new(bins.triangles);
    let same_point = R::from_f64(SAME_POINT_REL);

    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let from = [x, y, z];
                let i = shape.linearize(from) as usize;
                for axis in 0..3 {
                    let mut to = from;
                    to[axis] += 1;
                    if to[axis] >= size[axis] {
                        continue;
                    }
                    let j = shape.linearize(to) as usize;
                    report.probes += 1;

                    let (crossings, merged, coplanar, on_endpoint) = probe
                        .crossings(&bins, positions, indices, size, from, axis, at, same_point);
                    report.merged_crossings += merged;
                    report.coplanar_probes += coplanar;
                    if on_endpoint {
                        report.endpoint_crossings += 1;
                    }
                    if on_endpoint
                        || on_surface.get(i) == Some(&true)
                        || on_surface.get(j) == Some(&true)
                    {
                        report.degenerate_probes += 1;
                        continue;
                    }

                    let mesh_separates = crossings % 2 == 1;
                    let field_separates = air[i] != air[j];

                    if field_separates {
                        report.field_walls += 1;
                    }
                    if mesh_separates {
                        report.mesh_walls += 1;
                    }
                    match (field_separates, mesh_separates) {
                        (true, false) => {
                            report.unsealed_walls += 1;
                            if (0..3).any(|a| from[a] == 0 || to[a] == size[a] - 1) {
                                report.unsealed_on_domain_face += 1;
                            }
                        }
                        (false, true) => report.spurious_walls += 1,
                        _ => {}
                    }

                    if !field_separates && air[i] {
                        field_uf.union(i, j);
                    }
                    if !mesh_separates {
                        mesh_uf.union(i, j);
                        if air[i] && air[j] {
                            air_uf.union(i, j);
                        }
                    }
                }
            }
        }
    }

    for (i, &is_air) in air.iter().enumerate() {
        if is_air && field_uf.find(i) == i {
            report.field_air_components += 1;
        }
        if is_air && air_uf.find(i) == i {
            report.mesh_air_components += 1;
        }
    }

    // A region mixes phases when its members are not all one sign. Recorded per
    // root in one pass rather than by grouping, which would need a sort.
    let mut seen_air = alloc::vec![false; count];
    let mut seen_solid = alloc::vec![false; count];
    for (i, &is_air) in air.iter().enumerate() {
        let root = mesh_uf.find(i);
        let side = if is_air {
            &mut seen_air
        } else {
            &mut seen_solid
        };
        if let Some(slot) = side.get_mut(root) {
            *slot = true;
        }
    }
    for i in 0..count {
        if mesh_uf.find(i) == i {
            report.mesh_regions += 1;
            if seen_air.get(i) == Some(&true) && seen_solid.get(i) == Some(&true) {
                report.mixed_regions += 1;
            }
        }
    }

    report
}

/// Triangles binned into extraction cells, in CSR form.
///
/// One allocation pair instead of a `Vec` per cell. The bin is the **cell**
/// rather than a bounding-box hierarchy because the probes are grid edges, so
/// the candidate set for a probe is exactly the triangles in the at most four
/// cells that edge borders.
struct CellBins {
    cells: [u32; 3],
    starts: Vec<u32>,
    items: Vec<u32>,
    triangles: usize,
    degenerate: u64,
}

impl CellBins {
    fn build<R: Real>(
        size: [u32; 3],
        origin: [R; 3],
        cell_size: R,
        positions: &[[R; 3]],
        indices: &[u32],
    ) -> Self {
        let cells = [
            size[0].saturating_sub(1).max(1),
            size[1].saturating_sub(1).max(1),
            size[2].saturating_sub(1).max(1),
        ];
        let bin_count = cells[0] as usize * cells[1] as usize * cells[2] as usize;
        let triangles = indices.len() / 3;

        // A triangle's cell range, expanded by one cell each way. The expansion
        // is deliberate slack: a vertex sitting exactly on a cell boundary
        // belongs to the cell either side of it, and floor() picks one. Paying a
        // few extra narrow-phase tests is cheaper than reasoning about which.
        let range = |t: usize| -> Option<([u32; 3], [u32; 3])> {
            let tri = [
                *positions.get(*indices.get(t * 3)? as usize)?,
                *positions.get(*indices.get(t * 3 + 1)? as usize)?,
                *positions.get(*indices.get(t * 3 + 2)? as usize)?,
            ];
            if has_no_area(tri) {
                return None;
            }
            let mut lo = [0u32; 3];
            let mut hi = [0u32; 3];
            for a in 0..3 {
                let min = tri[0][a].min(tri[1][a]).min(tri[2][a]);
                let max = tri[0][a].max(tri[1][a]).max(tri[2][a]);
                if !min.is_finite() || !max.is_finite() {
                    return None;
                }
                let to_cell = |v: R| {
                    let s = ((v - origin[a]) / cell_size).floor().as_f64();
                    // Clamped in f64 before narrowing: `as` on an out-of-range
                    // float saturates in Rust, but the clamp says so out loud
                    // and covers a mesh far outside the grid.
                    s.clamp(0.0, f64::from(cells[a] - 1)) as u32
                };
                lo[a] = to_cell(min).saturating_sub(1);
                hi[a] = (to_cell(max) + 1).min(cells[a] - 1);
            }
            Some((lo, hi))
        };

        let index = |c: [u32; 3]| {
            (c[2] as usize * cells[1] as usize + c[1] as usize) * cells[0] as usize + c[0] as usize
        };

        let mut counts = alloc::vec![0u32; bin_count + 1];
        for t in 0..triangles {
            if let Some((lo, hi)) = range(t) {
                for cz in lo[2]..=hi[2] {
                    for cy in lo[1]..=hi[1] {
                        for cx in lo[0]..=hi[0] {
                            counts[index([cx, cy, cz]) + 1] += 1;
                        }
                    }
                }
            }
        }
        for i in 1..counts.len() {
            counts[i] += counts[i - 1];
        }
        let starts = counts.clone();
        let total = *counts.last().unwrap_or(&0) as usize;
        let mut items = alloc::vec![0u32; total];
        let mut cursor = counts;
        for t in 0..triangles {
            if let Some((lo, hi)) = range(t) {
                for cz in lo[2]..=hi[2] {
                    for cy in lo[1]..=hi[1] {
                        for cx in lo[0]..=hi[0] {
                            let b = index([cx, cy, cz]);
                            if let Some(slot) = items.get_mut(cursor[b] as usize) {
                                *slot = t as u32;
                            }
                            cursor[b] += 1;
                        }
                    }
                }
            }
        }

        let degenerate = (0..triangles).filter(|t| range(*t).is_none()).count() as u64;

        Self {
            cells,
            starts,
            items,
            triangles,
            degenerate,
        }
    }

    fn cell(&self, c: [u32; 3]) -> &[u32] {
        let b = (c[2] as usize * self.cells[1] as usize + c[1] as usize) * self.cells[0] as usize
            + c[0] as usize;
        let (Some(&lo), Some(&hi)) = (self.starts.get(b), self.starts.get(b + 1)) else {
            return &[];
        };
        self.items.get(lo as usize..hi as usize).unwrap_or(&[])
    }
}

/// Scratch reused across probes.
///
/// A stamp per triangle deduplicates the four incident cells' candidate lists in
/// `O(1)` per candidate without allocating; the parameter list is sorted per
/// probe and is at most a handful of entries.
struct Probe {
    stamp: Vec<u64>,
    epoch: u64,
    hits: Vec<f64>,
}

impl Probe {
    fn new(triangles: usize) -> Self {
        Self {
            stamp: alloc::vec![0u64; triangles],
            epoch: 0,
            hits: Vec::new(),
        }
    }

    /// Distinct crossings of one grid edge, how many hits merged, how many
    /// triangles were coplanar with it, and whether any crossing sat on an
    /// endpoint.
    #[allow(clippy::too_many_arguments, reason = "one probe, fully specified")]
    fn crossings<R: Real>(
        &mut self,
        bins: &CellBins,
        positions: &[[R; 3]],
        indices: &[u32],
        size: [u32; 3],
        from: [u32; 3],
        axis: usize,
        at: impl Fn([u32; 3]) -> [R; 3],
        same_point: R,
    ) -> (u64, u64, u64, bool) {
        self.epoch += 1;
        self.hits.clear();

        let mut to = from;
        to[axis] += 1;
        let p0 = at(from);
        let p1 = at(to);
        let d = vec3::sub(p1, p0);

        let mut coplanar = 0u64;

        // The at most four cells this edge borders: it runs the length of one
        // cell along `axis` and sits on a corner of the cells either side of it
        // in the other two.
        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
        let mut cell = [0u32; 3];
        cell[axis] = from[axis].min(bins.cells[axis] - 1);
        for du in 0..2u32 {
            let Some(cu) = (from[u] + du).checked_sub(1) else {
                continue;
            };
            if cu >= bins.cells[u] || from[u] + du > size[u] {
                continue;
            }
            for dv in 0..2u32 {
                let Some(cv) = (from[v] + dv).checked_sub(1) else {
                    continue;
                };
                if cv >= bins.cells[v] || from[v] + dv > size[v] {
                    continue;
                }
                cell[u] = cu;
                cell[v] = cv;
                for &t in bins.cell(cell) {
                    let slot = t as usize;
                    if self.stamp.get(slot).copied() == Some(self.epoch) {
                        continue;
                    }
                    if let Some(s) = self.stamp.get_mut(slot) {
                        *s = self.epoch;
                    }
                    match segment_triangle(p0, d, positions, indices, slot) {
                        Hit::At(t) => self.hits.push(t.as_f64()),
                        Hit::Coplanar => coplanar += 1,
                        Hit::Miss => {}
                    }
                }
            }
        }

        let raw = self.hits.len() as u64;
        self.hits.sort_by(f64::total_cmp);
        let tol = same_point.as_f64();
        let mut distinct = 0u64;
        let mut last = f64::NEG_INFINITY;
        let mut on_endpoint = false;
        for &t in &self.hits {
            if distinct == 0 || t - last > tol {
                distinct += 1;
                last = t;
                on_endpoint |= t <= tol || t >= 1.0 - tol;
            }
        }
        (distinct, raw - distinct, coplanar, on_endpoint)
    }
}

/// A triangle whose two edge vectors are parallel, so it has no area.
///
/// Tested on `|e₁ × e₂|` against `|e₁| · |e₂|` — the **sine** of the angle
/// between the edges — so it means the same thing at every grid spacing and does
/// not confuse a small triangle with a thin one.
fn has_no_area<R: Real>(tri: [[R; 3]; 3]) -> bool {
    let e1 = vec3::sub(tri[1], tri[0]);
    let e2 = vec3::sub(tri[2], tri[0]);
    vec3::length(vec3::cross(e1, e2))
        <= R::from_f64(PARALLEL_REL) * vec3::length(e1) * vec3::length(e2)
}

/// Outcome of one segment against one triangle.
enum Hit<R> {
    /// Crosses at this fraction along the segment.
    At(R),
    /// The segment lies in the triangle's plane; a crossing count is undefined.
    Coplanar,
    Miss,
}

/// Möller–Trumbore, with the parallel branch split into *coplanar* and *miss*.
///
/// Möller & Trumbore, *Fast, minimum storage ray/triangle intersection*, JGT
/// 1997 (`10.1080/10867651.1997.10487468`). The barycentric bounds are inclusive
/// so that a crossing exactly on a shared triangle edge is found by **both**
/// incident triangles — which is the case that then merges in
/// [`Probe::crossings`], and is the whole reason the merge exists.
fn segment_triangle<R: Real>(
    p0: [R; 3],
    d: [R; 3],
    positions: &[[R; 3]],
    indices: &[u32],
    t: usize,
) -> Hit<R> {
    let (Some(&ia), Some(&ib), Some(&ic)) = (
        indices.get(t * 3),
        indices.get(t * 3 + 1),
        indices.get(t * 3 + 2),
    ) else {
        return Hit::Miss;
    };
    let (Some(&a), Some(&b), Some(&c)) = (
        positions.get(ia as usize),
        positions.get(ib as usize),
        positions.get(ic as usize),
    ) else {
        return Hit::Miss;
    };

    let e1 = vec3::sub(b, a);
    let e2 = vec3::sub(c, a);
    let h = vec3::cross(d, e2);
    let det = vec3::dot(e1, h);

    // Relative to the three lengths, so this tests the sine of the angle
    // between segment and plane rather than a raw determinant, and means the
    // same thing at every grid spacing.
    let scale = vec3::length(d) * vec3::length(e1) * vec3::length(e2);
    if det.abs() <= R::from_f64(PARALLEL_REL) * scale {
        let s = vec3::sub(p0, a);
        let n = vec3::cross(e1, e2);
        // In-plane as well as parallel is coplanar; parallel and offset misses.
        //
        // The comparison is on the **distance** `|s·n| / |n|`, against the
        // segment length. An earlier version compared `|s·n|` against
        // `|n| · |d| · |e₁| · |e₂|`, which is a length⁵ threshold on a length³
        // quantity: it called 6,840 of Marching Tetrahedra's probes coplanar at
        // one grid spacing and would have said something different at another.
        return if vec3::dot(s, n).abs()
            <= R::from_f64(PARALLEL_REL) * vec3::length(n) * vec3::length(d)
        {
            Hit::Coplanar
        } else {
            Hit::Miss
        };
    }

    let inv = det.recip();
    let s = vec3::sub(p0, a);
    let u = vec3::dot(s, h) * inv;
    if u < R::ZERO || u > R::ONE {
        return Hit::Miss;
    }
    let q = vec3::cross(s, e1);
    let v = vec3::dot(d, q) * inv;
    if v < R::ZERO || u + v > R::ONE {
        return Hit::Miss;
    }
    let t = vec3::dot(e2, q) * inv;
    if t < R::ZERO || t > R::ONE {
        return Hit::Miss;
    }
    Hit::At(t)
}

/// Union-find with path halving and union by size.
///
/// Deterministic: the representative is fixed by the union order, and the union
/// order is the lattice scan, which is fixed. Component *counts* would not need
/// that; keeping it means a future caller can key on the representative.
struct UnionFind {
    parent: Vec<u32>,
    size: Vec<u32>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n as u32).collect(),
            size: alloc::vec![1u32; n],
        }
    }

    fn find(&mut self, mut i: usize) -> usize {
        while self.parent.get(i).copied().map(|p| p as usize) != Some(i) {
            let Some(p) = self.parent.get(i).copied() else {
                return i;
            };
            let grand = self.parent.get(p as usize).copied().unwrap_or(p);
            if let Some(slot) = self.parent.get_mut(i) {
                *slot = grand;
            }
            i = grand as usize;
        }
        i
    }

    fn union(&mut self, a: usize, b: usize) {
        let (mut ra, mut rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        if self.size.get(ra).copied().unwrap_or(0) < self.size.get(rb).copied().unwrap_or(0) {
            core::mem::swap(&mut ra, &mut rb);
        }
        if let Some(slot) = self.parent.get_mut(rb) {
            *slot = ra as u32;
        }
        let grew = self.size.get(rb).copied().unwrap_or(0);
        if let Some(slot) = self.size.get_mut(ra) {
            *slot += grew;
        }
    }
}

#[cfg(test)]
mod tests;
