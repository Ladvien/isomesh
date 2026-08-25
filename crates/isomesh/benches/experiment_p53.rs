//! **P-53 — the third corner label, and the slivers an exactly-equal corner makes.**
//!
//! Ticket: R-048. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p53
//! ```
//!
//! Writes `docs/experiments/p-53.csv`, eight rows: two volumes × two isovalues ×
//! two label rules.
//!
//! # What is being measured
//!
//! This crate classifies a corner with [`is_inside`], which is `value < 0`, so a
//! sample sitting *exactly* on the isosurface is outside. Either global choice is
//! legitimate (Lengyel, 2010 §3.1.1) and neither avoids the consequence: on a cut
//! edge whose outside endpoint is that exactly-zero corner, `t = a / (a − b)` is
//! `0 / (0 − b)` = **exactly 0**, so the crossing lands **on the corner**. Every
//! cut edge meeting there does the same thing, the vertex cache keys on the grid
//! edge rather than on the point, and the cell emits a triangle whose corners are
//! all the same point. Custodio, Pesco & Silva
//! (`10.1186/s13173-019-0086-6`) name the fix: a **third corner label**, so the
//! equal case is not silently folded into one of the two sides.
//!
//! The input class is not hypothetical. M-316 measured **16,284 of 529,508**
//! `bonsai` surface-cell corners exactly on the isosurface at an integer
//! isovalue — 3% — and M-232 measured 20 singular faces per 400,000 cells at
//! `u8` density against **0** in continuous data.
//!
//! # What is bench-local, and what is the crate's
//!
//! `src/**` is read-only for this experiment, so nothing here adds a label to
//! `is_inside`. Three things are rebuilt in this file and nothing else is:
//!
//! 1. **The march is replayed for provenance only.** Cell iteration order, the
//!    eight-corner case index, [`CASES`] and the vertex cache keyed on
//!    `(lower sample, axis)` — no position and no normal. That recovers, for
//!    every triangle the crate emitted, *which cell emitted it*, and for every
//!    vertex, *which grid edge it sits on*. The replay's index buffer is then
//!    compared against the crate's, element for element (`replay_matches_crate`);
//!    the tagging is only licensed because that comparison holds.
//! 2. **The ternary label**, as a pre-pass over the corners: `−` strictly
//!    inside, `+` strictly outside, `=` exactly on the isosurface. It is decided
//!    on the **raw byte** — `raw == 32` — not on `iso − raw` against `0.0`, which
//!    is the same answer for this data with one fewer float comparison in the
//!    gate, and is why the half-offset row is a control rather than decoration: a
//!    `u8` cannot equal 32.5 at all.
//! 3. **Snap and collapse.** Every vertex on an edge incident to an `=` corner is
//!    moved to that corner and the vertices sharing an `=` corner become one, at
//!    the lowest of their indices — the tie-break [`isomesh::weld`] already uses.
//!    Triangles that then name one vertex twice are dropped.
//!
//! The mesh, its positions and its normals are the crate's own
//! [`MarchingCubes`] output. The treatment arm is that mesh rewritten, so both
//! arms differ in the label rule and in nothing else.
//!
//! # The degenerate test carries no epsilon
//!
//! A triangle is degenerate here when **two of its three indices are the same
//! vertex**, or **two of its three positions are bit-identical**, or **the
//! doubled-area cross product is exactly the zero vector**. All three are exact:
//! integers, bit patterns, and a float compared against zero rather than against
//! a tuned threshold — so the count is the same on every machine, which a
//! relative-area threshold is not obliged to be. `degenerate_validator_epsilon`
//! records what `isomesh::validate` counts at its own
//! `area <= 1e-6 · h²` beside it, so the two definitions can be read against each
//! other rather than confused. `h = 1` here, one world unit per voxel, so that
//! threshold is `1e-6` absolute.
//!
//! # The snap is measured, not assumed
//!
//! With `origin = 0`, `h = 1` and integer sample values, `t` is exactly 0 or
//! exactly 1 on an edge incident to an `=` corner, and `lo + (hi − lo)·t` is
//! then bit-exactly a corner position. So the snap **should** move nothing and
//! the label's whole effect should be the collapse. `max_snap_distance` reports
//! the largest move actually made; if it is not exactly zero, the reasoning above
//! is wrong and every number in the row is suspect. The normal of a snapped
//! vertex is recomputed from the field at the snapped position by the crate's own
//! rule, so it cannot silently disagree with the position.
//!
//! `unexplained_coincident_vertices` is the other half of the same question:
//! vertices whose position is bit-identical to an earlier vertex's and which the
//! `=` label does **not** reach. That is the part of the coincidence the paper's
//! mechanism does not explain, and it is reported rather than welded away — a
//! blanket epsilon weld is `isomesh::weld`'s job (M-48) and is a different claim.
//!
//! # Why a collapse can move the topology, and the column that says when
//!
//! Two vertices snapped to the same `=` corner are in one of two situations, and
//! C2 turns entirely on which. If they already share a triangle, that triangle is
//! one of the degenerate ones and merging them flattens a fold: the triangle goes
//! and no edge, boundary or component can move. If they share no triangle, they
//! are on **different pieces of the surface** that happen to meet at that sample —
//! the isosurface genuinely touches itself there — and identifying the point is a
//! change of topology no relabelling can avoid. `pinch_groups` counts the second
//! kind and `pinch_excess_components` counts the pieces welded, over a
//! disjoint-set built from the baseline mesh's own triangles. When C2's three
//! topological columns move, this is the number that says why.
//!
//! # Source honesty
//!
//! **The source paper reports no degenerate-triangle count on any dataset** —
//! only radii-ratio histograms, Betti numbers and blocked-cube percentages — so
//! C2's 10× is this crate's own bar and is not a reproduction of a published
//! figure. Their triangulator is a per-cube convex hull with cross-cell face
//! dedup; it is **not** reproduced here and nothing here is a claim about it.
//!
//! `half_offset_identical` is recorded on all eight rows. At 32.5 it is C3 and it
//! is a gate. At 32 it is the same comparison and is expected to be false —
//! reported because "the two arms differ where equal corners exist" is worth
//! seeing beside "they do not differ where none do".
//!
//! Counted, not timed: every clause is an integer, an equality or a ratio of
//! integers. `extract_ms` is printed beside them and gates nothing.

#![allow(
    clippy::float_cmp,
    reason = "exact equality is the phenomenon: a corner on the isosurface, a \
              vertex on that corner, a triangle of exactly zero area"
)]

mod common;

use std::path::{Path, PathBuf};
use std::time::Instant;

use isomesh::construct::SampledField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_cubes::table::{self, CASES, EDGE_AXIS, EDGE_CORNERS, is_inside};
use isomesh::validate::{ValidateConfig, mesh_hash, validate_indexed};
use isomesh::{MeshBuffer, MeshSink, RuntimeShape3, Sdf, Shape3};

/// One world unit per voxel: all three datasets declare `spacing 1x1x1`.
const CELL: f64 = 1.0;

/// Sample `[0, 0, 0]` sits at the world origin, so a corner's world position is
/// its integer grid coordinate exactly.
const ORIGIN: [f64; 3] = [0.0; 3];

/// A volume to read.
struct Volume {
    file: &'static str,
    short: &'static str,
    /// Samples per axis; the files are cubes.
    n: u32,
}

const VOLUMES: [Volume; 2] = [
    Volume {
        file: "fuel_64x64x64_uint8.raw",
        short: "fuel",
        n: 64,
    },
    Volume {
        file: "bonsai_256x256x256_uint8.raw",
        short: "bonsai",
        n: 256,
    },
];

/// An isovalue, and the raw byte that can sit exactly on it.
///
/// The `=` label is an integer comparison against `equal`, which is why 32.5 is a
/// control: `None` means no sample can be on the isosurface, so the label has
/// nothing to act on and the two arms must come out identical.
#[derive(Clone, Copy)]
struct Iso {
    value: f64,
    label: &'static str,
    equal: Option<u8>,
}

const ISOS: [Iso; 2] = [
    Iso {
        value: 32.0,
        label: "32",
        equal: Some(32),
    },
    Iso {
        value: 32.5,
        label: "32.5",
        equal: None,
    },
];

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/measurements/volumes")
}

/// Read a `uint8` volume, raw. The length is checked against the dimensions
/// rather than trusted from the filename, as `benches/volumes.rs` does.
fn read_u8(path: &Path, n: u32) -> Result<Vec<u8>, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let want = (n as usize).pow(3);
    if bytes.len() != want {
        return Err(format!(
            "{}: {} bytes, expected {want} for {n}³",
            path.display(),
            bytes.len()
        ));
    }
    Ok(bytes)
}

/// The local offset of cube corner `i`, as grid steps: bit 0 is x, bit 1 is y,
/// bit 2 is z.
///
/// `crate::cube::corner_offset` is `pub(crate)`. This is its documented
/// definition rather than a second convention.
const fn corner_offset(corner: u8) -> [u32; 3] {
    [
        (corner & 1) as u32,
        ((corner >> 1) & 1) as u32,
        ((corner >> 2) & 1) as u32,
    ]
}

/// Marching Cubes' own vertex numbering and cell attribution, recovered without
/// recomputing a single position.
struct Provenance {
    /// The index buffer in emission order, to be compared against the crate's.
    indices: Vec<u32>,
    /// Per triangle: the linear sample index of its cell's base corner.
    cell: Vec<u32>,
    /// Per vertex: the lower sample of its grid edge, or `u32::MAX` for a
    /// cell-local cycle centroid.
    edge_lo: Vec<u32>,
    /// Per vertex: the axis of that grid edge, or `3` for a centroid.
    edge_axis: Vec<u8>,
    /// Cells with at least one corner inside and at least one not — M-316's
    /// `lo < 0 && hi >= 0`, in case-index form.
    surface_cells: u64,
    /// Cells whose triangulation needed a cycle centroid. Plain Marching Cubes
    /// tops out at cycle length 7 and `safe_apex` covers 3..=7, so this should be
    /// zero; it is counted rather than assumed, because a centroid vertex is
    /// allocated on a path the replay has to follow to stay in step.
    centroid_cells: u64,
    /// Distinct samples that are a corner of at least one surface cell.
    surface_cell_corners: u64,
    /// How many of those are exactly on the isosurface.
    equal_corners: u64,
    /// Surface cells whose **base** corner is exactly on the isosurface —
    /// M-316's own narrower census, reproduced to check this loader against the
    /// record (16,284 of 529,508 on `bonsai` at iso 32).
    m316_equal_base_corners: u64,
}

/// Replay the march for provenance, reading the same values the crate read.
fn replay(values: &[f64], raw: &[u8], equal: Option<u8>, shape: &RuntimeShape3) -> Provenance {
    let size = shape.size();
    let samples = shape.element_count();
    let mut p = Provenance {
        indices: Vec::new(),
        cell: Vec::new(),
        edge_lo: Vec::new(),
        edge_axis: Vec::new(),
        surface_cells: 0,
        centroid_cells: 0,
        surface_cell_corners: 0,
        equal_corners: 0,
        m316_equal_base_corners: 0,
    };
    // The crate's cache: one slot per grid edge, keyed on the lower sample plus
    // the axis, so the key is the same whichever cell arrives first.
    let mut edge_vertex = vec![u32::MAX; samples * 3];
    let mut corner_seen = vec![0u64; samples.div_ceil(64)];
    let mut next = 0u32;

    for z in 0..size[2] - 1 {
        for y in 0..size[1] - 1 {
            for x in 0..size[0] - 1 {
                let base = [x, y, z];
                let mut case = 0u8;
                let mut sample = [0u32; 8];
                for (c, slot) in sample.iter_mut().enumerate() {
                    let o = corner_offset(c as u8);
                    let s = shape.linearize([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
                    *slot = s;
                    if is_inside(values[s as usize]) {
                        case |= 1 << c;
                    }
                }
                if case == 0 || case == u8::MAX {
                    continue;
                }
                p.surface_cells += 1;
                if equal == Some(raw[sample[0] as usize]) {
                    p.m316_equal_base_corners += 1;
                }
                for s in sample {
                    corner_seen[s as usize / 64] |= 1 << (s % 64);
                }

                let entry = CASES[case as usize];
                if entry.count == 0 {
                    continue;
                }
                // Cycle centroids are allocated before any triangle of the cell
                // and are cell-local, so no grid edge names them.
                let mut centroid = [u32::MAX; table::MAX_CENTROIDS];
                if entry.centroids > 0 {
                    p.centroid_cells += 1;
                }
                for slot in centroid.iter_mut().take(entry.centroids as usize) {
                    *slot = next;
                    next += 1;
                    p.edge_lo.push(u32::MAX);
                    p.edge_axis.push(3);
                }

                let cell = shape.linearize(base);
                for tri in &entry.triangles[..entry.count as usize] {
                    for &code in tri {
                        let index = if table::is_centroid(code) {
                            centroid[(code - table::CENTROID_BASE) as usize]
                        } else {
                            let axis = EDGE_AXIS[code as usize];
                            let lo = sample[EDGE_CORNERS[code as usize][0] as usize];
                            let key = lo as usize * 3 + axis as usize;
                            if edge_vertex[key] == u32::MAX {
                                edge_vertex[key] = next;
                                next += 1;
                                p.edge_lo.push(lo);
                                p.edge_axis.push(axis);
                            }
                            edge_vertex[key]
                        };
                        p.indices.push(index);
                    }
                    p.cell.push(cell);
                }
            }
        }
    }

    for (w, word) in corner_seen.iter().enumerate() {
        let mut bits = *word;
        while bits != 0 {
            let s = w * 64 + bits.trailing_zeros() as usize;
            bits &= bits - 1;
            p.surface_cell_corners += 1;
            if equal == Some(raw[s]) {
                p.equal_corners += 1;
            }
        }
    }

    p
}

/// Bit pattern of a position, so coincidence is bitwise rather than approximate.
fn bits(p: [f64; 3]) -> [u64; 3] {
    [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()]
}

/// Doubled area, as a vector: `(b − a) × (c − a)`.
fn doubled_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

/// The crate's normal rule, restated over the same public [`Sdf::gradient`]:
/// `marching_cubes::unit_gradient` is private, and a snapped vertex has to get
/// the normal the crate would have given it or the arms differ in two things.
fn unit_gradient<S: Sdf<Scalar = f64>>(field: &S, p: [f64; 3]) -> [f64; 3] {
    let g = field.gradient(p);
    let inv = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt().recip();
    [g[0] * inv, g[1] * inv, g[2] * inv]
}

/// Degenerate triangles, and how many are attributable to an `=` corner.
#[derive(Clone, Copy, Default)]
struct Degeneracy {
    total: u64,
    repeated_index: u64,
    coincident_position: u64,
    zero_area_only: u64,
    from_equal_corners: u64,
}

impl Degeneracy {
    /// Share of degenerate triangles whose cell has a corner exactly on the
    /// isosurface. Zero when there is nothing to attribute.
    fn attributable(self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            self.from_equal_corners as f64 / self.total as f64
        }
    }
}

/// Does this cell have a corner whose raw value is exactly the isovalue?
fn cell_has_equal_corner(shape: &RuntimeShape3, raw: &[u8], equal: Option<u8>, cell: u32) -> bool {
    let Some(target) = equal else {
        return false;
    };
    let base = shape.delinearize(cell);
    (0..8u8).any(|c| {
        let o = corner_offset(c);
        let s = shape.linearize([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
        raw[s as usize] == target
    })
}

/// Count degenerate triangles and tag each one with its own cell's corners.
///
/// `cell[t]` is the base sample of the cell that emitted triangle `t`, so the
/// tag is read off that cell rather than inferred from a correlation.
fn census(
    mesh: &MeshBuffer<f64>,
    cell: &[u32],
    raw: &[u8],
    equal: Option<u8>,
    shape: &RuntimeShape3,
) -> Degeneracy {
    let mut d = Degeneracy::default();
    for (t, tri) in mesh.indices.as_chunks::<3>().0.iter().enumerate() {
        let idx = [tri[0], tri[1], tri[2]];
        let p = [
            mesh.positions[idx[0] as usize],
            mesh.positions[idx[1] as usize],
            mesh.positions[idx[2] as usize],
        ];
        let repeated = idx[0] == idx[1] || idx[1] == idx[2] || idx[2] == idx[0];
        let coincident =
            bits(p[0]) == bits(p[1]) || bits(p[1]) == bits(p[2]) || bits(p[2]) == bits(p[0]);
        let flat = doubled_area(p[0], p[1], p[2]) == [0.0; 3];
        if !(repeated || coincident || flat) {
            continue;
        }
        d.total += 1;
        if repeated {
            d.repeated_index += 1;
        }
        if coincident {
            d.coincident_position += 1;
        }
        if flat && !repeated && !coincident {
            d.zero_area_only += 1;
        }
        if cell_has_equal_corner(shape, raw, equal, cell[t]) {
            d.from_equal_corners += 1;
        }
    }
    d
}

/// Disjoint-set over vertex indices, path halving, unioned to the **lower**
/// root.
///
/// `validate::Dsu` is private and this needs the same property it has: the
/// result depends only on which unions were requested, never on the order they
/// arrived in, so the component count below is a pure function of the mesh.
struct Dsu(Vec<u32>);

impl Dsu {
    fn new(n: usize) -> Self {
        Self((0..n as u32).collect())
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.0[x as usize] != x {
            let parent = self.0[x as usize];
            self.0[x as usize] = self.0[parent as usize];
            x = self.0[x as usize];
        }
        x
    }

    fn union(&mut self, a: u32, b: u32) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
        self.0[hi as usize] = lo;
    }
}

/// What the ternary label did, and the mesh it produced.
struct Ternary {
    mesh: MeshBuffer<f64>,
    /// Per surviving triangle: its cell's base sample, as in [`Provenance`].
    cell: Vec<u32>,
    snapped_vertices: u64,
    collapsed_groups: u64,
    vertices_removed: u64,
    triangles_dropped: u64,
    max_snap: f64,
    /// Vertices whose position is bit-identical to an earlier vertex's and which
    /// the `=` label does not reach.
    unexplained_coincident: u64,
    /// Collapse groups whose members did **not** all already share triangles —
    /// so the merge joins mesh pieces that were connected to nothing but each
    /// other's position.
    ///
    /// This is the number C2 turns on. A group whose members all co-occur in
    /// triangles is a fold being flattened: the degenerate triangles go and the
    /// topology cannot move. A group spanning two pieces is a **pinch** — the
    /// surface genuinely touches itself at that sample, and identifying the
    /// point is a change of topology whatever the label rule is called.
    pinch_groups: u64,
    /// Summed `components − 1` over those groups: how many separate pieces the
    /// collapse welded together in total.
    pinch_excess_components: u64,
}

/// Label, snap, collapse.
///
/// The label is the pre-pass: a ternary sign per corner, read off the raw bytes.
/// Every vertex on an edge incident to an `=` corner snaps to that corner, the
/// vertices sharing one become the lowest of their indices, and a triangle that
/// then names a vertex twice is dropped. Nothing else moves: with an empty label
/// set this returns the input mesh unchanged, vertex for vertex, which is what
/// makes the half-offset row a control.
fn ternary(
    field: &SampledField<'_, f64, RuntimeShape3>,
    base: &MeshBuffer<f64>,
    prov: &Provenance,
    raw: &[u8],
    equal: Option<u8>,
    shape: &RuntimeShape3,
) -> Ternary {
    let size = shape.size();
    let stride = [1, size[0], size[0] * size[1]];
    let n = base.positions.len();

    // ── the label, and what each vertex snaps to ────────────────────────────
    let mut target = vec![u32::MAX; n];
    let mut pairs: Vec<(u32, u32)> = Vec::new();
    if let Some(t) = equal {
        for (v, slot) in target.iter_mut().enumerate() {
            let lo = prov.edge_lo[v];
            if lo == u32::MAX {
                continue;
            }
            let hi = lo + stride[prov.edge_axis[v] as usize];
            // A cut edge has one endpoint strictly inside, and `=` is on the
            // outside of `is_inside`, so at most one of the two can be `=`.
            let s = if raw[lo as usize] == t {
                lo
            } else if raw[hi as usize] == t {
                hi
            } else {
                continue;
            };
            *slot = s;
            pairs.push((s, v as u32));
        }
    }
    // Sorted, then scanned in runs: no map, no iteration order to leak. The key
    // ends in the vertex index, so no two entries compare equal.
    pairs.sort_unstable();

    let mut remap: Vec<u32> = (0..n as u32).collect();
    let mut collapsed_groups = 0u64;
    let mut groups = 0u64;
    let mut i = 0;
    while i < pairs.len() {
        let mut j = i + 1;
        while j < pairs.len() && pairs[j].0 == pairs[i].0 {
            j += 1;
        }
        groups += 1;
        if j - i > 1 {
            collapsed_groups += 1;
        }
        // The lowest-indexed member represents the group, which is the tie-break
        // `isomesh::weld` documents.
        let rep = pairs[i].1;
        for &(_, v) in &pairs[i..j] {
            remap[v as usize] = rep;
        }
        i = j;
    }

    // ── the snap, measured ──────────────────────────────────────────────────
    let corner_position = |s: u32| {
        let c = shape.delinearize(s);
        [
            ORIGIN[0] + CELL * f64::from(c[0]),
            ORIGIN[1] + CELL * f64::from(c[1]),
            ORIGIN[2] + CELL * f64::from(c[2]),
        ]
    };
    let mut max_snap = 0.0f64;
    for (v, &s) in target.iter().enumerate() {
        if s == u32::MAX {
            continue;
        }
        let to = corner_position(s);
        let from = base.positions[v];
        for a in 0..3 {
            max_snap = max_snap.max((to[a] - from[a]).abs());
        }
    }

    // ── the collapsed mesh ──────────────────────────────────────────────────
    let mut mesh = MeshBuffer::<f64>::new();
    let mut new_index = vec![u32::MAX; n];
    for (v, slot) in new_index.iter_mut().enumerate() {
        if remap[v] as usize != v {
            continue;
        }
        let (position, normal) = if target[v] == u32::MAX {
            (base.positions[v], base.normals[v])
        } else {
            let p = corner_position(target[v]);
            (p, unit_gradient(field, p))
        };
        *slot = mesh.vertex(position, normal);
    }

    let mut cell = Vec::new();
    let mut triangles_dropped = 0u64;
    for (t, tri) in base.indices.as_chunks::<3>().0.iter().enumerate() {
        let a = new_index[remap[tri[0] as usize] as usize];
        let b = new_index[remap[tri[1] as usize] as usize];
        let c = new_index[remap[tri[2] as usize] as usize];
        if a == b || b == c || c == a {
            triangles_dropped += 1;
            continue;
        }
        mesh.triangle(a, b, c);
        cell.push(prov.cell[t]);
    }

    // ── coincidence the label does not explain ──────────────────────────────
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_unstable_by(|&a, &b| {
        let (p, q) = (base.positions[a as usize], base.positions[b as usize]);
        // `total_cmp`, so a NaN coordinate sorts into view instead of vanishing
        // through a partial comparison.
        p[0].total_cmp(&q[0])
            .then(p[1].total_cmp(&q[1]))
            .then(p[2].total_cmp(&q[2]))
            .then(a.cmp(&b))
    });
    let mut unexplained_coincident = 0u64;
    let mut i = 0;
    while i < order.len() {
        let key = bits(base.positions[order[i] as usize]);
        let mut j = i + 1;
        while j < order.len() && bits(base.positions[order[j] as usize]) == key {
            j += 1;
        }
        if j - i > 1 {
            for &v in &order[i..j] {
                if target[v as usize] == u32::MAX {
                    unexplained_coincident += 1;
                }
            }
        }
        i = j;
    }

    // ── which collapses join pieces, rather than flatten folds ──────────────
    //
    // Two vertices of one group that already share a triangle are two corners of
    // a triangle about to be dropped: the fold flattens and nothing moves. Two
    // that share no triangle are on different pieces of the surface, and
    // identifying their point is a pinch. Counted rather than argued, because
    // this is the difference between C2 holding and C2 failing.
    let mut dsu = Dsu::new(n);
    for tri in base.indices.as_chunks::<3>().0 {
        for (a, b) in [(0, 1), (1, 2), (2, 0)] {
            let (u, w) = (tri[a], tri[b]);
            if u != w && remap[u as usize] == remap[w as usize] {
                dsu.union(u, w);
            }
        }
    }
    let mut pinch_groups = 0u64;
    let mut pinch_excess_components = 0u64;
    let mut g = 0;
    while g < pairs.len() {
        let mut j = g + 1;
        while j < pairs.len() && pairs[j].0 == pairs[g].0 {
            j += 1;
        }
        if j - g > 1 {
            let mut roots: Vec<u32> = pairs[g..j].iter().map(|&(_, v)| dsu.find(v)).collect();
            roots.sort_unstable();
            roots.dedup();
            if roots.len() > 1 {
                pinch_groups += 1;
                pinch_excess_components += roots.len() as u64 - 1;
            }
        }
        g = j;
    }

    Ternary {
        mesh,
        cell,
        snapped_vertices: pairs.len() as u64,
        collapsed_groups,
        vertices_removed: pairs.len() as u64 - groups,
        triangles_dropped,
        max_snap,
        unexplained_coincident,
        pinch_groups,
        pinch_excess_components,
    }
}

fn main() {
    let prereg = isomesh::experiment!("P-53");

    common::experiment::run(prereg, |run| {
        let dir = dir();
        for v in &VOLUMES {
            let raw = match read_u8(&dir.join(v.file), v.n) {
                Ok(raw) => raw,
                Err(e) => {
                    println!("::error:: {e}");
                    std::process::exit(1);
                }
            };
            let shape = match RuntimeShape3::new([v.n; 3]) {
                Ok(shape) => shape,
                Err(e) => {
                    println!("::error:: {}: {e}", v.file);
                    std::process::exit(1);
                }
            };
            let cells = u64::from(v.n - 1).pow(3);

            for iso in ISOS {
                // `iso - value`, so a dense voxel is negative and the crate's
                // sign convention holds unchanged.
                let values: Vec<f64> = raw.iter().map(|b| iso.value - f64::from(*b)).collect();
                let field = match SampledField::new(&values, &shape, ORIGIN, CELL) {
                    Ok(field) => field,
                    Err(e) => {
                        println!("::error:: {}: {e}", v.file);
                        std::process::exit(1);
                    }
                };

                // ── baseline: the crate's own Marching Cubes ────────────────
                let mut base = MeshBuffer::<f64>::new();
                let t0 = Instant::now();
                if let Err(e) =
                    MarchingCubes::<f64>::new().extract(&field, &shape, ORIGIN, CELL, &mut base)
                {
                    println!("::error:: {} at iso {}: {e}", v.file, iso.label);
                    std::process::exit(1);
                }
                let extract_ms = t0.elapsed().as_secs_f64() * 1e3;

                let prov = replay(&values, &raw, iso.equal, &shape);
                let matches =
                    prov.indices == base.indices && prov.edge_lo.len() == base.positions.len();
                if !matches {
                    println!(
                        "::error:: {} at iso {}: the replay is not the crate's march \
                         ({} vs {} indices, {} vs {} vertices) — every attribution \
                         below would be a guess",
                        v.short,
                        iso.label,
                        prov.indices.len(),
                        base.indices.len(),
                        prov.edge_lo.len(),
                        base.positions.len()
                    );
                    std::process::exit(1);
                }

                let cfg = match ValidateConfig::from_cell_size(CELL) {
                    Ok(cfg) => cfg,
                    Err(e) => {
                        println!("::error:: {e}");
                        std::process::exit(1);
                    }
                };

                let base_deg = census(&base, &prov.cell, &raw, iso.equal, &shape);
                let base_report = validate_indexed(&base.positions, &base.indices, &cfg);
                let base_hash = mesh_hash(&base);

                let tern = ternary(&field, &base, &prov, &raw, iso.equal, &shape);
                let tern_deg = census(&tern.mesh, &tern.cell, &raw, iso.equal, &shape);
                let tern_report = validate_indexed(&tern.mesh.positions, &tern.mesh.indices, &cfg);
                let tern_hash = mesh_hash(&tern.mesh);
                let identical = base_hash == tern_hash;

                // `0 / 0` is undefined, not infinite: at 32.5 there was nothing
                // to remove, and writing `inf` there would read as an infinite
                // improvement over nothing.
                let ratio = match (base_deg.total, tern_deg.total) {
                    (0, 0) => String::from("n/a"),
                    (_, 0) => String::from("inf"),
                    (b, t) => format!("{:.4}", b as f64 / t as f64),
                };

                println!(
                    "{:>7} iso {:<5} cells {cells:>10}  surface-cell corners {:>8} \
                     ({} exactly equal)",
                    v.short, iso.label, prov.surface_cell_corners, prov.equal_corners
                );
                println!(
                    "        binary   tris {:>8}  degenerate {:>6} ({} from equal corners, \
                     {:.4})  chi {:>7}  nm_e {:>5}  bnd {:>5}",
                    base.triangle_count(),
                    base_deg.total,
                    base_deg.from_equal_corners,
                    base_deg.attributable(),
                    base_report.euler_characteristic,
                    base_report.non_manifold_edges,
                    base_report.boundary_edges
                );
                println!(
                    "        ternary  tris {:>8}  degenerate {:>6} ({} from equal corners, \
                     {:.4})  chi {:>7}  nm_e {:>5}  bnd {:>5}  ratio {ratio}",
                    tern.mesh.triangle_count(),
                    tern_deg.total,
                    tern_deg.from_equal_corners,
                    tern_deg.attributable(),
                    tern_report.euler_characteristic,
                    tern_report.non_manifold_edges,
                    tern_report.boundary_edges
                );
                println!(
                    "        snapped {} vertices into {} groups, removed {}, dropped {} \
                     triangles, max snap {:e}, unexplained coincidence {}, identical {identical}",
                    tern.snapped_vertices,
                    tern.collapsed_groups,
                    tern.vertices_removed,
                    tern.triangles_dropped,
                    tern.max_snap,
                    tern.unexplained_coincident
                );
                println!(
                    "        of those groups {} are pinches, welding {} separate pieces \
                     together in total",
                    tern.pinch_groups, tern.pinch_excess_components
                );

                let shared = |rule: &'static str| {
                    vec![
                        ("volume", v.short.to_string()),
                        ("isovalue", iso.label.to_string()),
                        ("label_rule", rule.to_string()),
                        ("cells", cells.to_string()),
                        ("samples_per_axis", v.n.to_string()),
                        ("surface_cells", prov.surface_cells.to_string()),
                        (
                            "surface_cell_corners",
                            prov.surface_cell_corners.to_string(),
                        ),
                        ("equal_corners", prov.equal_corners.to_string()),
                        (
                            "m316_equal_base_corners",
                            prov.m316_equal_base_corners.to_string(),
                        ),
                        ("centroid_cells", prov.centroid_cells.to_string()),
                        ("replay_matches_crate", matches.to_string()),
                        ("half_offset_identical", identical.to_string()),
                        ("extract_ms", format!("{extract_ms:.3}")),
                    ]
                };

                let mut binary = shared("binary");
                binary.extend([
                    ("triangles", base.triangle_count().to_string()),
                    ("vertices", base.vertex_count().to_string()),
                    ("degenerate_triangles", base_deg.total.to_string()),
                    (
                        "degenerate_from_equal_corners",
                        base_deg.from_equal_corners.to_string(),
                    ),
                    (
                        "degenerate_attributable_fraction",
                        format!("{:.6}", base_deg.attributable()),
                    ),
                    ("degenerate_ratio", String::from("1")),
                    (
                        "degenerate_repeated_index",
                        base_deg.repeated_index.to_string(),
                    ),
                    (
                        "degenerate_coincident_position",
                        base_deg.coincident_position.to_string(),
                    ),
                    (
                        "degenerate_zero_area_only",
                        base_deg.zero_area_only.to_string(),
                    ),
                    (
                        "degenerate_validator_epsilon",
                        base_report.degenerate_triangles.to_string(),
                    ),
                    (
                        "euler_characteristic",
                        base_report.euler_characteristic.to_string(),
                    ),
                    (
                        "non_manifold_edges",
                        base_report.non_manifold_edges.to_string(),
                    ),
                    ("boundary_edges", base_report.boundary_edges.to_string()),
                    ("mesh_hash", base_hash.to_string()),
                    ("snapped_vertices", String::from("0")),
                    ("collapsed_groups", String::from("0")),
                    ("vertices_removed", String::from("0")),
                    ("triangles_dropped", String::from("0")),
                    ("max_snap_distance", String::from("0")),
                    (
                        "unexplained_coincident_vertices",
                        tern.unexplained_coincident.to_string(),
                    ),
                    ("pinch_groups", String::from("0")),
                    ("pinch_excess_components", String::from("0")),
                ]);
                run.record(&binary);

                let mut ternary_row = shared("ternary");
                ternary_row.extend([
                    ("triangles", tern.mesh.triangle_count().to_string()),
                    ("vertices", tern.mesh.vertex_count().to_string()),
                    ("degenerate_triangles", tern_deg.total.to_string()),
                    (
                        "degenerate_from_equal_corners",
                        tern_deg.from_equal_corners.to_string(),
                    ),
                    (
                        "degenerate_attributable_fraction",
                        format!("{:.6}", tern_deg.attributable()),
                    ),
                    ("degenerate_ratio", ratio.clone()),
                    (
                        "degenerate_repeated_index",
                        tern_deg.repeated_index.to_string(),
                    ),
                    (
                        "degenerate_coincident_position",
                        tern_deg.coincident_position.to_string(),
                    ),
                    (
                        "degenerate_zero_area_only",
                        tern_deg.zero_area_only.to_string(),
                    ),
                    (
                        "degenerate_validator_epsilon",
                        tern_report.degenerate_triangles.to_string(),
                    ),
                    (
                        "euler_characteristic",
                        tern_report.euler_characteristic.to_string(),
                    ),
                    (
                        "non_manifold_edges",
                        tern_report.non_manifold_edges.to_string(),
                    ),
                    ("boundary_edges", tern_report.boundary_edges.to_string()),
                    ("mesh_hash", tern_hash.to_string()),
                    ("snapped_vertices", tern.snapped_vertices.to_string()),
                    ("collapsed_groups", tern.collapsed_groups.to_string()),
                    ("vertices_removed", tern.vertices_removed.to_string()),
                    ("triangles_dropped", tern.triangles_dropped.to_string()),
                    ("max_snap_distance", format!("{:e}", tern.max_snap)),
                    (
                        "unexplained_coincident_vertices",
                        tern.unexplained_coincident.to_string(),
                    ),
                    ("pinch_groups", tern.pinch_groups.to_string()),
                    (
                        "pinch_excess_components",
                        tern.pinch_excess_components.to_string(),
                    ),
                ]);
                run.record(&ternary_row);
            }
        }
    });
}
