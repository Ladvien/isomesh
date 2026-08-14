//! Counting triangle pairs that pass through each other.
//!
//! # What is reported, and why
//!
//! The literature's usual statistic is `p`, the *fraction of meshes with at
//! least one self-intersection*. For a chunked engine that number is actively
//! misleading: with an intersection rate `λ` per triangle and `T` triangles per
//! chunk, `p = 1 − e^{−λT}`, so a large enough chunk reads `p ≈ 1` for any
//! `λ > 0`. This module reports **`λ`** —
//! [`per_thousand_triangles`](SelfIntersectionReport::per_thousand_triangles).
//!
//! It also returns the intersecting **pairs**, not just their count. The
//! research is explicit that every detected pair should be bucketed as
//! intra-cell, inter-cell within a chunk, or cross-chunk-seam, because
//! cross-seam intersections are a stitching bug with nothing to do with
//! contouring and will otherwise dominate `λ` and send you chasing a problem you
//! do not have. Chunks do not exist yet, so rather than invent a chunk concept
//! here, the pairs come back and the bucketing becomes a caller's job the moment
//! there is something to bucket by.
//!
//! # Cost
//!
//! A uniform grid broadphase and an exact narrow phase. Not fast, and not meant
//! to be: this measures a mesh, it does not produce one.

use alloc::vec::Vec;
use core::fmt;

use crate::Real;
use crate::vec3;

/// Tolerance for deciding two triangles are coplanar, relative to `cell_size`.
const COPLANAR_EPSILON_REL: f64 = 1e-6;

use super::tri_grid::{MAX_CELLS_PER_TRIANGLE, cell_of};

/// Which triangles pass through which.
#[derive(Clone, Debug, PartialEq)]
pub struct SelfIntersectionReport {
    /// Triangles considered.
    pub triangles: u64,
    /// Intersecting pairs, as triangle indices, each `[low, high]`.
    ///
    /// Sorted ascending, so the report is a pure function of the mesh. Returned
    /// rather than merely counted so that a caller can bucket them once chunks
    /// exist — see the module docs.
    pub pairs: Vec<[u32; 2]>,
    /// Pairs that reached the exact test. A measure of how well the grid
    /// separated the mesh, not a defect.
    pub tested_pairs: u64,
    /// Pairs rejected for sharing at least one vertex index.
    ///
    /// Adjacent triangles touch along their shared edge or at their shared
    /// vertex by construction. Counting that as a self-intersection would make
    /// every well-formed mesh look catastrophic — this is the trap the ticket
    /// names, and this counter is how you can see the filter is working.
    pub adjacent_pairs_skipped: u64,
    /// Triangles with no well-defined plane, excluded from the exact test.
    ///
    /// A zero-area triangle has no normal to intersect against. They are counted
    /// here and reported by the validity harness as degenerate; this module does
    /// not silently pretend they were tested.
    pub degenerate_triangles: u64,
    /// The grid spacing the report was produced with.
    pub cell_size: f64,
}

impl SelfIntersectionReport {
    /// Number of intersecting pairs.
    #[must_use]
    pub fn count(&self) -> u64 {
        self.pairs.len() as u64
    }

    /// Intersecting pairs per 1,000 triangles — the rate `λ`.
    ///
    /// This is the number to record and compare. Zero triangles gives zero.
    #[must_use]
    pub fn per_thousand_triangles(&self) -> f64 {
        if self.triangles == 0 {
            0.0
        } else {
            1000.0 * self.count() as f64 / self.triangles as f64
        }
    }

    /// `true` when no pair intersects.
    #[must_use]
    pub fn is_intersection_free(&self) -> bool {
        self.pairs.is_empty()
    }
}

impl fmt::Display for SelfIntersectionReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "self-intersection report")?;
        writeln!(f, "  triangles                {:8}", self.triangles)?;
        writeln!(f, "  intersecting pairs       {:8}", self.count())?;
        writeln!(
            f,
            "  per 1000 triangles       {:8.3}",
            self.per_thousand_triangles()
        )?;
        writeln!(f, "  narrow-phase tests       {:8}", self.tested_pairs)?;
        writeln!(
            f,
            "  adjacent pairs skipped   {:8}",
            self.adjacent_pairs_skipped
        )?;
        writeln!(
            f,
            "  degenerate triangles     {:8}",
            self.degenerate_triangles
        )?;
        write!(f, "  cell size                {:8}", self.cell_size)
    }
}

/// Count triangle pairs that pass through each other.
///
/// `cell_size` is the broadphase grid spacing, normally the spacing of the grid
/// the mesh was extracted from. It is required rather than guessed, for the same
/// reason the validity harness requires one: a geometric tolerance without a
/// length scale is meaningless.
///
/// # What counts
///
/// A pair counts when the triangles overlap in a set of positive measure —
/// either crossing transversely, or overlapping while coplanar. A *tangential*
/// contact, where they meet exactly along a line or at a point, does not: that
/// is measure zero, it is what correctly-stitched neighbours do, and treating it
/// as an intersection would make every closed mesh report a defect.
///
/// Pairs sharing at least one vertex index are excluded outright. Adjacent
/// triangles necessarily touch, so including them would swamp the rate. The
/// documented cost of that choice: a fold that pinches *exactly* at a shared
/// vertex is not counted. The alternative — testing vertex-adjacent pairs for a
/// proper crossing — would report a false positive for every pair in every
/// vertex fan unless it also handled the touching case exactly, and that
/// trade-off is not worth it for a measurement.
///
/// Triangles with malformed indices are skipped, exactly as in the validity
/// harness, so a broken mesh yields a report rather than a panic.
///
/// # Errors
///
/// [`Error::InvalidCellSize`](crate::Error::InvalidCellSize) if `cell_size` is
/// not finite and positive, and
/// [`Error::CellSizeMismatch`](crate::Error::CellSizeMismatch) if a single
/// triangle spans more grid cells than the guard allows — which means the
/// spacing does not describe this mesh, and the broadphase would otherwise grow
/// until it exhausted memory.
pub fn self_intersections<R: Real>(
    positions: &[[R; 3]],
    indices: &[u32],
    cell_size: f64,
) -> crate::Result<SelfIntersectionReport> {
    if !cell_size.is_finite() || cell_size <= 0.0 {
        return Err(crate::Error::InvalidCellSize { value: cell_size });
    }

    // Same face filter as the validity harness: an out-of-range index cannot be
    // dereferenced and a repeated index has no plane.
    let whole = indices.len() - indices.len() % 3;
    let mut tris: Vec<[u32; 3]> = Vec::with_capacity(whole / 3);
    for tri in indices[..whole].chunks_exact(3) {
        let in_range = tri.iter().all(|&i| (i as usize) < positions.len());
        let distinct = tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2];
        if in_range && distinct {
            tris.push([tri[0], tri[1], tri[2]]);
        }
    }

    let mut report = SelfIntersectionReport {
        triangles: tris.len() as u64,
        pairs: Vec::new(),
        tested_pairs: 0,
        adjacent_pairs_skipped: 0,
        degenerate_triangles: 0,
        cell_size,
    };

    let eps = R::from_f64(COPLANAR_EPSILON_REL * cell_size);
    let inv_cell = R::from_f64(cell_size).recip();

    // Unit normals up front: a triangle with no well-defined plane is excluded
    // once here rather than rediscovered inside every pair test.
    let mut normals: Vec<Option<[R; 3]>> = Vec::with_capacity(tris.len());
    for t in &tris {
        let (a, b, c) = corners(positions, *t);
        let n = vec3::cross(vec3::sub(b, a), vec3::sub(c, a));
        let len = vec3::length(n);
        if len > R::ZERO && len.is_finite() {
            normals.push(Some(vec3::scale(n, len.recip())));
        } else {
            normals.push(None);
            report.degenerate_triangles += 1;
        }
    }

    // ── broadphase: a uniform grid, as a sorted Vec ─────────────────────────
    //
    // Same sort-and-scan shape as the validity harness, and deterministic for
    // the same reason: the key ends in the triangle index, so no two entries
    // compare equal.
    let mut buckets: Vec<([i64; 3], u32)> = Vec::new();
    for (ti, t) in tris.iter().enumerate() {
        if normals[ti].is_none() {
            continue;
        }
        let (a, b, c) = corners(positions, *t);
        let lo = [
            cell_of(a[0].min(b[0]).min(c[0]) * inv_cell),
            cell_of(a[1].min(b[1]).min(c[1]) * inv_cell),
            cell_of(a[2].min(b[2]).min(c[2]) * inv_cell),
        ];
        let hi = [
            cell_of(a[0].max(b[0]).max(c[0]) * inv_cell),
            cell_of(a[1].max(b[1]).max(c[1]) * inv_cell),
            cell_of(a[2].max(b[2]).max(c[2]) * inv_cell),
        ];
        let span =
            (hi[0] - lo[0] + 1) as u128 * (hi[1] - lo[1] + 1) as u128 * (hi[2] - lo[2] + 1) as u128;
        if span > MAX_CELLS_PER_TRIANGLE as u128 {
            return Err(crate::Error::CellSizeMismatch {
                triangle: ti as u64,
                cells: span,
                cell_size,
            });
        }
        for z in lo[2]..=hi[2] {
            for y in lo[1]..=hi[1] {
                for x in lo[0]..=hi[0] {
                    buckets.push(([x, y, z], ti as u32));
                }
            }
        }
    }
    buckets.sort_unstable();

    // Candidate pairs, deduplicated by sorting rather than by a hash set: a pair
    // sharing several cells is found several times.
    let mut candidates: Vec<[u32; 2]> = Vec::new();
    let mut start = 0usize;
    while start < buckets.len() {
        let key = buckets[start].0;
        let mut end = start + 1;
        while end < buckets.len() && buckets[end].0 == key {
            end += 1;
        }
        let run = &buckets[start..end];
        for i in 0..run.len() {
            for j in i + 1..run.len() {
                let (a, b) = (run[i].1, run[j].1);
                candidates.push(if a < b { [a, b] } else { [b, a] });
            }
        }
        start = end;
    }
    candidates.sort_unstable();
    candidates.dedup();

    // ── narrow phase ────────────────────────────────────────────────────────
    for [i, j] in candidates {
        let (ti, tj) = (tris[i as usize], tris[j as usize]);
        if shares_a_vertex(ti, tj) {
            report.adjacent_pairs_skipped += 1;
            continue;
        }
        let (Some(ni), Some(nj)) = (normals[i as usize], normals[j as usize]) else {
            continue;
        };
        report.tested_pairs += 1;

        let a = triangle(positions, ti);
        let b = triangle(positions, tj);
        if triangles_overlap(a, ni, b, nj, eps) {
            report.pairs.push([i, j]);
        }
    }

    Ok(report)
}

#[inline]
fn corners<R: Real>(positions: &[[R; 3]], t: [u32; 3]) -> ([R; 3], [R; 3], [R; 3]) {
    (
        positions[t[0] as usize],
        positions[t[1] as usize],
        positions[t[2] as usize],
    )
}

#[inline]
fn triangle<R: Real>(positions: &[[R; 3]], t: [u32; 3]) -> [[R; 3]; 3] {
    [
        positions[t[0] as usize],
        positions[t[1] as usize],
        positions[t[2] as usize],
    ]
}

#[inline]
fn shares_a_vertex(a: [u32; 3], b: [u32; 3]) -> bool {
    a.iter().any(|x| b.contains(x))
}

/// Exact triangle-triangle overlap.
///
/// Both normals are unit, so a plane distance is a true distance and `eps` has
/// length units.
///
/// The non-coplanar case is the standard interval argument: if the triangles
/// meet at all they both meet the line where their planes cross, so each one's
/// intersection with that line is an interval and the triangles overlap exactly
/// when the intervals do. The intervals are compared along the coordinate axis
/// most aligned with the line, which is a monotone reparametrisation of it and
/// therefore order-preserving.
fn triangles_overlap<R: Real>(
    a: [[R; 3]; 3],
    na: [R; 3],
    b: [[R; 3]; 3],
    nb: [R; 3],
    eps: R,
) -> bool {
    let da = [
        vec3::dot(nb, vec3::sub(a[0], b[0])),
        vec3::dot(nb, vec3::sub(a[1], b[0])),
        vec3::dot(nb, vec3::sub(a[2], b[0])),
    ];
    if da[0].abs() <= eps && da[1].abs() <= eps && da[2].abs() <= eps {
        return coplanar_overlap(a, b, na);
    }
    // A positive-measure transverse crossing needs each triangle strictly on
    // *both* sides of the other's plane. "Not entirely on one side" is weaker:
    // a triangle with a vertex or an edge in the plane and the rest on one
    // side merely *touches* it, its interval on the crossing line degenerates
    // to that contact, and counting it would report exactly the tangential
    // point/line contact the module contract excludes.
    if !straddles(da, eps) {
        return false;
    }

    let db = [
        vec3::dot(na, vec3::sub(b[0], a[0])),
        vec3::dot(na, vec3::sub(b[1], a[0])),
        vec3::dot(na, vec3::sub(b[2], a[0])),
    ];
    if !straddles(db, eps) {
        return false;
    }

    let axis = vec3::dominant_axis(vec3::cross(na, nb));
    let (Some((lo_a, hi_a)), Some((lo_b, hi_b))) = (
        plane_interval(a, da, axis, eps),
        plane_interval(b, db, axis, eps),
    ) else {
        return false;
    };

    // Strict: intervals that merely touch describe a tangential contact.
    hi_a > lo_b && hi_b > lo_a
}

/// Whether the triangle has vertices strictly on both sides of the plane.
#[inline]
fn straddles<R: Real>(d: [R; 3], eps: R) -> bool {
    let above = d[0] > eps || d[1] > eps || d[2] > eps;
    let below = d[0] < -eps || d[1] < -eps || d[2] < -eps;
    above && below
}

/// Where a triangle meets the other plane, projected onto `axis`.
///
/// Collects the vertices lying in the plane and the crossings of edges whose
/// endpoints are strictly on opposite sides. Handling both together is what
/// makes a vertex exactly on the plane an ordinary case rather than one needing
/// its own branch — which is where hand-written versions of this usually go
/// wrong.
fn plane_interval<R: Real>(t: [[R; 3]; 3], d: [R; 3], axis: usize, eps: R) -> Option<(R, R)> {
    let mut lo = R::INFINITY;
    let mut hi = R::NEG_INFINITY;
    let mut any = false;

    let mut record = |v: R| {
        lo = lo.min(v);
        hi = hi.max(v);
        any = true;
    };

    for k in 0..3 {
        if d[k].abs() <= eps {
            record(t[k][axis]);
        }
    }
    for k in 0..3 {
        let m = (k + 1) % 3;
        let (dk, dm) = (d[k], d[m]);
        if (dk > eps && dm < -eps) || (dk < -eps && dm > eps) {
            let s = dk / (dk - dm);
            record(t[k][axis] + (t[m][axis] - t[k][axis]) * s);
        }
    }

    if any { Some((lo, hi)) } else { None }
}

/// Two coplanar triangles, by separating axis.
///
/// Projected onto the coordinate plane the shared normal is least aligned with,
/// which keeps the projection non-degenerate. Six candidate axes, one per edge;
/// convexity makes that sufficient. Touching counts as separated, matching the
/// transverse case.
fn coplanar_overlap<R: Real>(a: [[R; 3]; 3], b: [[R; 3]; 3], n: [R; 3]) -> bool {
    let drop = vec3::dominant_axis(n);
    let (u, v) = match drop {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };
    let flat = |t: [[R; 3]; 3]| [[t[0][u], t[0][v]], [t[1][u], t[1][v]], [t[2][u], t[2][v]]];
    let (p, q) = (flat(a), flat(b));

    for tri in [p, q] {
        for k in 0..3 {
            let m = (k + 1) % 3;
            // Outward-facing candidate axis for this edge.
            let axis = [-(tri[m][1] - tri[k][1]), tri[m][0] - tri[k][0]];
            let (pa, pb) = (project2(p, axis), project2(q, axis));
            if pa.1 <= pb.0 || pb.1 <= pa.0 {
                return false;
            }
        }
    }
    true
}

#[inline]
fn project2<R: Real>(t: [[R; 2]; 3], axis: [R; 2]) -> (R, R) {
    let d = [
        t[0][0] * axis[0] + t[0][1] * axis[1],
        t[1][0] * axis[0] + t[1][1] * axis[1],
        t[2][0] * axis[0] + t[2][1] * axis[1],
    ];
    (d[0].min(d[1]).min(d[2]), d[0].max(d[1]).max(d[2]))
}

#[cfg(test)]
mod tests;
