//! A-013 — welding coincident vertices, and the lattice that finds them.
//!
//! # What this is for, and what it is not for
//!
//! Within a single volume *almost* nothing needs welding: Marching Cubes keys
//! its vertices on the grid edge they sit on and the dual methods key theirs on
//! the cell, so a shared vertex is shared by construction. Dong et al. (2018,
//! `10.1109/icra.2018.8463157`) reach the same design from the other end — their
//! whole contribution is distributing vertices on discrete edges so duplicates
//! cannot be allocated in the first place.
//!
//! The word "almost" is measured (M-48) and replaced a flat "nothing" that was
//! written here and was wrong. The edge cache shares vertices between cells
//! meeting on a grid **edge**, and that is all it can do. When a grid **sample**
//! lands on the isosurface, `t` is 0 or 1 and the crossing sits *at that
//! sample*, so every cut edge meeting there places its own vertex at the same
//! point and nothing shares them. On `sphere` at 25³ that is 48 vertices and 96
//! collapsed triangles — and the 96 is exactly the degenerate-sliver count A-001
//! recorded for that resolution, so **welding also removes that class of
//! sliver**, which nobody predicted.
//!
//! **The bulk of the duplicates are still at chunk seams.** Two chunks meshed independently each
//! compute the vertices on their shared plane, so every seam vertex exists
//! twice, and no amount of care inside one extraction removes it. That is what
//! this module is for, and it is why [`crate::chunk`] had to land first.
//!
//! # Why an epsilon rather than equality
//!
//! Because equality does not hold. **M-32:** two adjacent chunks compute the
//! same seam plane as `(o + h·cn) + h·n` and as `o + h·(c+1)n`, equal by algebra
//! and not by IEEE, and 22% of random `(origin, h, cells, chunk)` combinations
//! disagree by an ulp. Only a power-of-two cell size makes them bit-identical.
//! A weld keyed on exact equality would therefore work on `h = 0.125` and
//! silently leave a seam un-welded on `h = 4/35`, which is the worst kind of
//! bug: correct on the fixture and wrong in the field.
//!
//! # The rule, stated because it is a choice
//!
//! Epsilon-closeness is **not transitive**, so "weld everything within ε" does
//! not define equivalence classes and cannot be implemented as stated. What is
//! well defined is first fit against the vertices already kept:
//!
//! > Walk vertices in input index order. A vertex joins the **lowest-indexed
//! > representative** within `ε` of it, or becomes a representative itself.
//!
//! Two consequences worth stating rather than discovering. A vertex is compared
//! against representatives only, never against a vertex that has itself been
//! welded away — so a chain `a ~ b ~ c` with `a` and `c` further than `ε` apart
//! yields two representatives, not one, and no vertex is ever moved further than
//! `ε` from where the extractor put it. And the result depends on the input
//! order, which is why the order is the buffer's own and why T-004 covers this:
//! weld ordering is the classic determinism leak, and
//! [`crate::MeshSink`]'s contract already demands that a welding sink document
//! its epsilon and its tie-break.
//!
//! The kept vertex keeps its own position and normal. Averaging would move a
//! vertex the extractor placed deliberately, and at a seam the two candidates
//! differ by about an ulp anyway.
//!
//! # The lattice
//!
//! The `Lattice` this welds on is shared with [`crate::validate`] — one type,
//! not two copies — deliberately: the validator's
//! `duplicate_vertices` scan and this welder must agree about which cells are
//! probed, or the count would describe a different neighbourhood than the weld.
//! They differ in exactly one thing — the validator asks whether *any* earlier
//! vertex is within `ε`, the welder asks for the lowest-indexed *representative*
//! — and that difference is measured rather than assumed.

#[cfg(test)]
mod tests;

use alloc::vec::Vec;

use crate::{Error, MeshBuffer, Real, Result};

/// Floor to an `i64` lattice coordinate, saturating a non-finite input to zero.
///
/// A NaN coordinate lands in the origin cell and then fails the exact distance
/// test below, so it is never welded to anything.
pub(crate) fn quantise<R: Real>(scaled: R) -> i64 {
    let f = scaled.floor();
    if f.is_finite() { f.as_f32() as i64 } else { 0 }
}

/// The integer lattice a set of positions is bucketed onto, of side `epsilon`.
///
/// # Why the anchor is the mesh's own minimum, not the world origin
///
/// [`quantise`] narrows through `as_f32`, which stops distinguishing consecutive
/// integers above `2²⁴`. The scale is `1/epsilon`, which for the default
/// `h · 1e-4` is a factor of 160,000 at `h = 0.0625` — so an *absolute*
/// coordinate crosses that ceiling at about **105 world units** (M-18, measured;
/// refined at T-008 to be gradual rather than a cliff). A chunked world passes
/// 105 units almost immediately.
///
/// Anchoring at the minimum makes the scale depend on the mesh's *extent*
/// instead of its position, and an extent that large is a mesh nobody could weld
/// meaningfully anyway. Same reasoning as `validate::tri_grid`.
pub(crate) struct Lattice<R: Real> {
    anchor: [R; 3],
    inv_epsilon: R,
}

impl<R: Real> Lattice<R> {
    /// Anchored on the componentwise minimum of the finite positions.
    pub(crate) fn new(positions: &[[R; 3]], epsilon: R) -> Self {
        let mut anchor = [R::INFINITY; 3];
        for p in positions {
            for (a, slot) in anchor.iter_mut().enumerate() {
                if p[a].is_finite() && p[a] < *slot {
                    *slot = p[a];
                }
            }
        }
        for slot in &mut anchor {
            if !slot.is_finite() {
                *slot = R::ZERO;
            }
        }
        Self {
            anchor,
            inv_epsilon: epsilon.recip(),
        }
    }

    /// Which cell a position falls in.
    pub(crate) fn key_of(&self, p: [R; 3]) -> [i64; 3] {
        [
            quantise((p[0] - self.anchor[0]) * self.inv_epsilon),
            quantise((p[1] - self.anchor[1]) * self.inv_epsilon),
            quantise((p[2] - self.anchor[2]) * self.inv_epsilon),
        ]
    }
}

/// What a weld did.
///
/// Counts rather than a bare success, because "how much was there to weld" is
/// the number that says whether the seams were where you thought they were.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct WeldReport {
    /// Vertices before.
    pub vertices_before: usize,
    /// Vertices after.
    pub vertices_after: usize,
    /// Triangles before.
    pub triangles_before: usize,
    /// Triangles after.
    pub triangles_after: usize,
    /// Triangles dropped because welding made two of their corners the same
    /// vertex.
    ///
    /// Such a triangle has no area left. Dropping it is not a repair of
    /// something that went wrong — it is the correct consequence of merging two
    /// of its corners, and the count is here so a surprising number is visible
    /// rather than silent.
    pub triangles_collapsed: usize,
}

impl WeldReport {
    /// How many vertices the weld removed.
    #[must_use]
    pub const fn vertices_removed(&self) -> usize {
        self.vertices_before - self.vertices_after
    }

    /// `true` when nothing coincided and the mesh is unchanged.
    #[must_use]
    pub const fn is_noop(&self) -> bool {
        self.vertices_before == self.vertices_after && self.triangles_collapsed == 0
    }
}

/// Welds coincident vertices in place, reusing its scratch across calls.
///
/// Owns its buffers for the same reason the extractors do: a chunked world welds
/// thousands of meshes and should not allocate thousands of times.
///
/// # Example
///
/// ```
/// use isomesh::weld::Welder;
/// use isomesh::MeshBuffer;
///
/// // Two triangles that share an edge, written as two separate triangles.
/// let mut mesh = MeshBuffer::<f64>::new();
/// mesh.positions = vec![
///     [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
///     [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0],
/// ];
/// mesh.normals = vec![[0.0, 0.0, 1.0]; 6];
/// mesh.indices = vec![0, 1, 2, 3, 4, 5];
///
/// let mut welder = Welder::<f64>::new();
/// let report = welder.weld(&mut mesh, 1e-6)?;
///
/// assert_eq!(report.vertices_removed(), 2);
/// assert_eq!(mesh.vertex_count(), 4);
/// assert_eq!(mesh.triangle_count(), 2);
/// # Ok::<(), isomesh::Error>(())
/// ```
#[derive(Debug)]
pub struct Welder<R: Real> {
    /// `(lattice cell, vertex index)`, sorted. The broadphase.
    cells: Vec<([i64; 3], u32)>,
    /// Output index of every input vertex.
    remap: Vec<u32>,
    /// Whether each input vertex was kept.
    kept: Vec<bool>,
    _scalar: core::marker::PhantomData<R>,
}

impl<R: Real> Welder<R> {
    /// A welder that has allocated nothing yet.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            cells: Vec::new(),
            remap: Vec::new(),
            kept: Vec::new(),
            _scalar: core::marker::PhantomData,
        }
    }

    /// The output index each input vertex was mapped to, from the last
    /// [`weld`](Self::weld).
    ///
    /// Exposed because a caller with parallel per-vertex data — colours, texture
    /// coordinates, material ids — has to move it the same way, and could not
    /// otherwise.
    #[must_use]
    pub fn remap(&self) -> &[u32] {
        &self.remap
    }

    /// Weld `mesh` in place.
    ///
    /// Vertices within `epsilon` of a kept vertex are merged into it by the rule
    /// in the [module docs](self); indices are rewritten; triangles left with two
    /// equal corners are dropped and counted. Positions, normals and indices are
    /// compacted without reallocating, since the output is never longer than the
    /// input.
    ///
    /// Unreferenced vertices are **not** removed. Welding is not garbage
    /// collection, and a vertex no triangle mentions is the caller's to explain;
    /// [`crate::validate`] already counts them.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidWeldEpsilon`] if `epsilon` is not finite and positive. A
    /// zero epsilon would weld only bit-identical positions, which M-32 says is
    /// exactly the case that fails at a seam, so it is rejected rather than
    /// quietly accepted.
    ///
    /// # Panics
    ///
    /// In debug builds, if `mesh.normals` is not the same length as
    /// `mesh.positions`. That is [`MeshBuffer`]'s own invariant — every
    /// [`MeshSink::vertex`](crate::MeshSink::vertex) writes both — so a mismatch
    /// is a caller that hand-built a malformed buffer.
    pub fn weld(&mut self, mesh: &mut MeshBuffer<R>, epsilon: R) -> Result<WeldReport> {
        if !epsilon.is_finite() || epsilon <= R::ZERO {
            // Narrowed only to carry the value in the message; the decision above
            // is made at full width.
            return Err(Error::InvalidWeldEpsilon {
                value: f64::from(epsilon.as_f32()),
            });
        }
        debug_assert_eq!(
            mesh.normals.len(),
            mesh.positions.len(),
            "a MeshBuffer carries one normal per vertex"
        );

        let n = mesh.positions.len();
        let mut report = WeldReport {
            vertices_before: n,
            vertices_after: n,
            triangles_before: mesh.indices.len() / 3,
            triangles_after: mesh.indices.len() / 3,
            triangles_collapsed: 0,
        };
        if n == 0 {
            return Ok(report);
        }

        // ── broadphase ──────────────────────────────────────────────────────
        let lattice = Lattice::new(&mesh.positions, epsilon);
        self.cells.clear();
        self.cells.reserve(n);
        for (i, p) in mesh.positions.iter().enumerate() {
            self.cells.push((lattice.key_of(*p), i as u32));
        }
        // The vertex index is part of the key and is unique, so no two entries
        // compare equal and an unstable sort is still a deterministic one.
        self.cells.sort_unstable();

        // ── elect representatives ───────────────────────────────────────────
        self.remap.clear();
        self.remap.resize(n, u32::MAX);
        self.kept.clear();
        self.kept.resize(n, false);

        let eps_sq = epsilon * epsilon;
        let mut next_output = 0u32;

        for v in 0..n {
            let p = mesh.positions[v];
            let base = lattice.key_of(p);

            // The *lowest-indexed* representative within epsilon, not the first
            // one the probe happens to reach — so the answer does not depend on
            // the order the 27 cells are visited in.
            let mut best: Option<u32> = None;
            for dz in -1..=1i64 {
                for dy in -1..=1i64 {
                    for dx in -1..=1i64 {
                        let key = [base[0] + dx, base[1] + dy, base[2] + dz];
                        let lo = self.cells.partition_point(|(k, _)| *k < key);
                        for &(k, u) in &self.cells[lo..] {
                            if k != key {
                                break;
                            }
                            if u as usize >= v || !self.kept[u as usize] {
                                continue;
                            }
                            if best.is_some_and(|b| u >= b) {
                                continue;
                            }
                            let q = mesh.positions[u as usize];
                            let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
                            if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= eps_sq {
                                best = Some(u);
                            }
                        }
                    }
                }
            }

            match best {
                Some(u) => self.remap[v] = self.remap[u as usize],
                None => {
                    self.kept[v] = true;
                    self.remap[v] = next_output;
                    next_output += 1;
                }
            }
        }

        report.vertices_after = next_output as usize;
        if report.vertices_after == n {
            // Nothing coincided, so no index can have moved and no triangle can
            // have collapsed. Returning here is not a second path — the loops
            // below would run and change nothing.
            return Ok(report);
        }

        // ── compact, in place ───────────────────────────────────────────────
        //
        // A kept vertex's output index is never greater than its input index, so
        // the write cursor never overtakes the read cursor.
        let mut w = 0usize;
        for v in 0..n {
            if self.kept[v] {
                mesh.positions[w] = mesh.positions[v];
                mesh.normals[w] = mesh.normals[v];
                w += 1;
            }
        }
        debug_assert_eq!(w, report.vertices_after);
        mesh.positions.truncate(w);
        mesh.normals.truncate(w);

        let mut t = 0usize;
        for tri in 0..report.triangles_before {
            let a = self.remap[mesh.indices[tri * 3] as usize];
            let b = self.remap[mesh.indices[tri * 3 + 1] as usize];
            let c = self.remap[mesh.indices[tri * 3 + 2] as usize];
            if a == b || b == c || a == c {
                report.triangles_collapsed += 1;
                continue;
            }
            mesh.indices[t * 3] = a;
            mesh.indices[t * 3 + 1] = b;
            mesh.indices[t * 3 + 2] = c;
            t += 1;
        }
        mesh.indices.truncate(t * 3);
        report.triangles_after = t;

        Ok(report)
    }
}

impl<R: Real> Default for Welder<R> {
    fn default() -> Self {
        Self::new()
    }
}
