//! **P-41 — is the sign lattice well-composed, and is that where the dual
//! extractors go non-manifold?**
//!
//! Ticket: R-040. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p41
//! ```
//!
//! Writes `docs/experiments/p-41.csv`.
//!
//! # What "critical configuration" means here, and where it comes from
//!
//! A binary digital set is **well-composed** exactly when its boundary is a
//! 2-manifold. Latecki, Eckhardt & Rosenfeld, *Well-Composed Sets*, CVIU 1995
//! (`10.1006/cviu.1995.1006`) and Latecki, *3D Well-Composed Pictures*,
//! GMIP 1997 (`10.1006/gmip.1997.0422`) characterise that by the **absence** of
//! two critical configurations in every `2×2×2` block of samples — which, on
//! this crate's grid, is exactly one **cell**, since a cell's corners are a
//! `2×2×2` block. Boutry, Géraud & Najman's tutorial
//! (`10.1007/s10851-017-0769-6`) states the same two in the same form.
//!
//! - **2D-critical** — within one of the cell's six `2×2` faces the two inside
//!   corners are one diagonal of that face and the two outside corners are the
//!   other. A checkerboard face; the pair shares only an edge of the cell.
//! - **3D-critical** — the cell's inside set is exactly two corners that differ
//!   in all three coordinates, so they share only a cell vertex; or the
//!   complementary case, six inside corners whose two outside corners are such a
//!   pair.
//!
//! **The classification is enumerated, not transcribed.** [`classify`] walks all
//! `256` sign bytes and decides each one combinatorially from the definitions
//! above — face sets built from the axis they fix, diagonality from an XOR of
//! corner indices. Rule 5 of `CLAUDE.md` forbids writing a case table from
//! memory, and a 256-entry table of critical configurations is precisely the
//! kind of thing that would be wrong in one entry and unfalsifiable.
//!
//! Corner bit layout: bit `i` is corner `(x, y, z)` with `i = x + 2y + 4z`, so
//! two corners are cell-diagonal exactly when `i ^ j == 0b111`, and two corners
//! of the face fixing axis `a` are face-diagonal exactly when
//! `i ^ j == 0b111 ^ (1 << a)`.
//!
//! # The sign test, and the `-0.0` trap
//!
//! Inside is `value < 0.0`, copied from `cube.rs::is_inside` rather than
//! reinvented — the census has to partition the samples the *extractor's* way or
//! it is a census of a different lattice. Note that `-0.0 < 0.0` is **false**,
//! so a negative zero is outside; `box_exact` is exactly zero across its whole
//! boundary, so this is reachable rather than theoretical.
//!
//! # The non-manifold incidents are recomputed here, and cross-checked
//!
//! [`recompute`] mirrors `validate.rs`'s two manifoldness passes **exactly**:
//!
//! - a **non-manifold edge** is an undirected edge `[lo, hi]` used by three or
//!   more faces (`validate.rs`'s run-length scan over canonicalised half-edges,
//!   `run.len()` in the `_ =>` arm);
//! - a **non-manifold vertex** is a vertex whose incident-face link is not one
//!   fan — either three or more incident faces share one wing vertex
//!   (`branching`, an umbrella) or the faces fall into more than one component
//!   under "shares a wing vertex" (`roots.len() > 1`, a bowtie).
//!
//! Same pass-0 skip rule (out-of-range or repeated indices), same
//! canonicalisation, same visit order, so the lists come out in the same order.
//! The one deliberate difference is union-by-index instead of the crate's
//! union-by-size, which changes which node becomes a root and cannot change the
//! partition.
//!
//! The crate *also* exposes the locations directly, through
//! `validate::validate_features` — so this is not a substitute for the crate's
//! answer, it is a second opinion against it. Every row carries both counts and
//! a `locations_agree` flag comparing the full lists element by element. A
//! disagreement is a finding about `validate.rs`, not a licence to prefer
//! whichever number is more convenient.
//!
//! # Two interpretations this harness had to choose, stated
//!
//! 1. **An edge incident spans two cells.** A dual method places one vertex per
//!    cell, so a non-manifold *vertex* names one cell and there is nothing to
//!    decide. A non-manifold *edge* joins two dual vertices and therefore two
//!    cells. The registered `incidents_in_critical_cells` uses the **union**
//!    rule — the incident counts as co-located if *either* endpoint's cell is
//!    flagged — because that is the reading under which "this incident occurred
//!    in a critical cell" is true. The stricter both-endpoints rule is recorded
//!    beside it as `incidents_in_critical_both_ends` so the choice is visible
//!    and reversible.
//! 2. **`colocation_fraction` is per row; the clause is decided pooled.** A row
//!    with zero incidents has no fraction, and inventing `1.0` for it would let
//!    seven clean fields outvote the one dirty one. Such a row records `n/a`,
//!    and the columns `pooled_*` (per extractor) and `overall_*` (both) carry
//!    the numbers clause two is actually decided on, repeated on every row.
//!
//! `null_colocation_fraction` is the control that makes 90% mean something: the
//! fraction of *vertex-hosting* cells that are critical. If the critical cells
//! are most of the active grid, a high co-location is arithmetic rather than
//! mechanism.
//!
//! # Not the banked topological safety gate
//!
//! That gates a vertex reposition, after extraction, on the mesh. This is a
//! census of the sign field, before extraction. Nothing here repairs anything:
//! the registration is the detector, and the repair moves the surface by up to a
//! cell and breaks every golden hash, so it is only worth designing if clause
//! two holds.

mod common;

use isomesh::dual_contouring::DualContouring;
use isomesh::surface_nets::SurfaceNets;
use isomesh::validate::{ValidateConfig, validate_features};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// The registered resolution. `65` samples span `64` cells per axis.
const SAMPLES: u32 = 65;

/// Inside, the way the extractors decide it.
///
/// Copied from `cube.rs::is_inside`: strictly negative. `-0.0` is **not** less
/// than `0.0`, so a negative zero is outside here exactly as it is there.
fn is_inside(value: f64) -> bool {
    value < 0.0
}

// ─── the 256-byte classification, enumerated ────────────────────────────────

/// Which of the 256 possible cell sign bytes host each critical configuration.
struct Critical {
    /// A checkerboard `2×2` face: the inside pair shares only a cell edge.
    two_d: [bool; 256],
    /// A main-diagonal inside pair, or its complement: the pair shares only a
    /// cell vertex.
    three_d: [bool; 256],
}

/// The inside corners of `byte`, and how many there are.
fn inside_corners(byte: u32) -> ([u32; 8], usize) {
    let mut out = [0u32; 8];
    let mut n = 0;
    for corner in 0..8u32 {
        if (byte >> corner) & 1 == 1 {
            out[n] = corner;
            n += 1;
        }
    }
    (out, n)
}

/// Exactly two inside corners, differing in all three coordinates.
///
/// `i ^ j == 0b111` is "differs in x, y and z", which for two corners of a cell
/// is "shares only a vertex".
fn is_vertex_diagonal_pair(byte: u32) -> bool {
    let (corners, n) = inside_corners(byte);
    n == 2 && (corners[0] ^ corners[1]) == 0b111
}

/// Some `2×2` face of the cell is a checkerboard.
///
/// The face fixing axis `a` at side `s` is the four corners with bit `a` equal
/// to `s`; its two diagonals are the corner pairs with `i ^ j == 0b111 ^
/// (1 << a)`. A checkerboard is two inside corners forming one of those
/// diagonals, which forces the other two — the outside pair — to be the other.
fn has_checkerboard_face(byte: u32) -> bool {
    for axis in 0..3u32 {
        let diagonal = 0b111 ^ (1 << axis);
        for side in 0..2u32 {
            let mut inside = [0u32; 4];
            let mut n = 0;
            for corner in 0..8u32 {
                if (corner >> axis) & 1 == side && (byte >> corner) & 1 == 1 {
                    inside[n] = corner;
                    n += 1;
                }
            }
            if n == 2 && (inside[0] ^ inside[1]) == diagonal {
                return true;
            }
        }
    }
    false
}

/// Decide all 256 sign bytes from the definitions, once.
fn classify() -> Critical {
    let mut two_d = [false; 256];
    let mut three_d = [false; 256];
    for byte in 0..256u32 {
        two_d[byte as usize] = has_checkerboard_face(byte);
        // The complementary case is the same configuration seen from the other
        // side, and Latecki lists both: six inside corners whose two outside
        // corners share only a vertex is just as non-well-composed.
        three_d[byte as usize] =
            is_vertex_diagonal_pair(byte) || is_vertex_diagonal_pair(!byte & 0xFF);
    }
    Critical { two_d, three_d }
}

// ─── the grid ───────────────────────────────────────────────────────────────

/// The grid a field is sampled and meshed on. One definition, so the census and
/// the extraction cannot disagree about which lattice they looked at.
struct Grid {
    /// World position of sample `[0, 0, 0]`.
    origin: [f64; 3],
    /// Spacing.
    cell_size: f64,
    /// Samples per axis.
    samples: u32,
    /// Cells per axis: `samples - 1`.
    cells: u32,
}

impl Grid {
    /// `i = x + y·sx + z·sx·sy`, the crate's order.
    fn sample_index(&self, x: usize, y: usize, z: usize) -> usize {
        let n = self.samples as usize;
        x + y * n + z * n * n
    }

    /// The same order over cells rather than samples.
    fn cell_index(&self, cell: [u32; 3]) -> usize {
        let c = self.cells as usize;
        cell[0] as usize + cell[1] as usize * c + cell[2] as usize * c * c
    }

    /// World position of sample `(x, y, z)`.
    fn point(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        [
            self.origin[0] + self.cell_size * x as f64,
            self.origin[1] + self.cell_size * y as f64,
            self.origin[2] + self.cell_size * z as f64,
        ]
    }

    /// The cell a dual vertex belongs to.
    ///
    /// Both dual extractors clamp their vertex strictly inside its own cell —
    /// `dual_contouring::apply_clamp` with `Clamp::ToCell` for the QEF path,
    /// the same function for the relaxation path in `dual.rs` — so the floor is
    /// exact rather than a guess, and the clamp is `(1 − ε)` about the cell
    /// centre so no vertex lands on a face. The clamp to `0..cells-1` is
    /// belt-and-braces for a vertex that escaped anyway; `escaped_cell` counts
    /// those instead of hiding them.
    fn cell_of(&self, p: [f64; 3]) -> ([u32; 3], bool) {
        let mut cell = [0u32; 3];
        let mut escaped = false;
        let last = self.cells - 1;
        for (axis, slot) in cell.iter_mut().enumerate() {
            let t = ((p[axis] - self.origin[axis]) / self.cell_size).floor();
            if t < 0.0 {
                escaped = true;
                *slot = 0;
            } else if t > f64::from(last) {
                escaped = true;
                *slot = last;
            } else {
                *slot = t as u32;
            }
        }
        (cell, escaped)
    }
}

// ─── the census ─────────────────────────────────────────────────────────────

/// The critical-configuration census of one field's sign lattice.
struct Census {
    /// One flag per cell, `true` when the cell hosts either configuration.
    critical: Vec<bool>,
    /// The cells that ought to host a dual vertex, ascending: those whose eight
    /// corner signs are not all equal, so at least one of the twelve cell edges
    /// is cut. A dual method places one vertex per such cell, and comparing this
    /// set against the set the vertex positions map back to is what turns "the
    /// vertex is in its cell" from an assumption into a measurement — the
    /// registration names the escaping vertex as the rival suspect, so it is
    /// checked rather than trusted.
    active: Vec<usize>,
    /// Cells hosting a checkerboard face.
    two_d: u64,
    /// Cells hosting a main-diagonal pair or its complement.
    three_d: u64,
    /// Cells hosting either. Not the sum: the classes are disjoint by
    /// construction, and this asserts nothing about that.
    either: u64,
}

/// Sample the field on the extraction grid, take signs, classify every cell.
fn census<F>(field: &F, grid: &Grid, table: &Critical) -> Census
where
    F: Sdf<Scalar = f64>,
{
    let n = grid.samples as usize;
    let mut inside = vec![false; n * n * n];
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                inside[grid.sample_index(x, y, z)] = is_inside(field.sample(grid.point(x, y, z)));
            }
        }
    }

    let c = grid.cells as usize;
    let mut critical = vec![false; c * c * c];
    let mut active = Vec::new();
    let mut two_d = 0u64;
    let mut three_d = 0u64;
    let mut either = 0u64;
    for cz in 0..c {
        for cy in 0..c {
            for cx in 0..c {
                let mut byte = 0u32;
                for corner in 0..8usize {
                    let dx = corner & 1;
                    let dy = (corner >> 1) & 1;
                    let dz = (corner >> 2) & 1;
                    if inside[grid.sample_index(cx + dx, cy + dy, cz + dz)] {
                        byte |= 1 << corner;
                    }
                }
                let flat = table.two_d[byte as usize];
                let solid = table.three_d[byte as usize];
                if flat {
                    two_d += 1;
                }
                if solid {
                    three_d += 1;
                }
                let index = grid.cell_index([cx as u32, cy as u32, cz as u32]);
                // Pushed in ascending index order, because `cx` is the fastest
                // loop and the index is `cx + cy·c + cz·c²`.
                if byte != 0x00 && byte != 0xFF {
                    active.push(index);
                }
                if flat || solid {
                    either += 1;
                    critical[index] = true;
                }
            }
        }
    }

    Census {
        critical,
        active,
        two_d,
        three_d,
        either,
    }
}

// ─── the non-manifold recomputation ─────────────────────────────────────────

/// Union-find over face indices. Union by index rather than by size, which
/// changes which node is a root and cannot change the partition.
#[derive(Default)]
struct Dsu {
    parent: Vec<u32>,
}

impl Dsu {
    fn reset(&mut self, n: usize) {
        self.parent.clear();
        self.parent.extend(0..n as u32);
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
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[rb as usize] = ra;
        }
    }
}

/// Where the mesh is non-manifold, recomputed from `positions`/`indices` by the
/// definitions in `validate.rs`.
struct Recomputed {
    /// Undirected edges used by three or more faces, ascending by `[lo, hi]`.
    edges: Vec<[u32; 2]>,
    /// Vertices whose incident-face link is not one fan, ascending.
    vertices: Vec<u32>,
}

/// `validate.rs`'s pass 0: the faces that can be dereferenced at all.
fn usable_faces(vertex_count: usize, indices: &[u32]) -> Vec<[u32; 3]> {
    let whole = indices.len() - indices.len() % 3;
    let mut faces = Vec::with_capacity(whole / 3);
    for tri in indices[..whole].chunks_exact(3) {
        let out_of_range = tri.iter().any(|&i| i as usize >= vertex_count);
        let repeated = tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2];
        if out_of_range || repeated {
            continue;
        }
        faces.push([tri[0], tri[1], tri[2]]);
    }
    faces
}

/// `validate.rs`'s pass 1, kept only for the `>= 3` arm: an undirected edge used
/// by three or more faces.
fn non_manifold_edges(faces: &[[u32; 3]]) -> Vec<[u32; 2]> {
    let mut half: Vec<(u32, u32, u32)> = Vec::with_capacity(faces.len() * 3);
    for (fi, f) in faces.iter().enumerate() {
        for k in 0..3usize {
            let (a, b) = (f[k], f[(k + 1) % 3]);
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            half.push((lo, hi, fi as u32));
        }
    }
    half.sort_unstable();

    let mut found = Vec::new();
    let mut start = 0usize;
    while start < half.len() {
        let (lo, hi, _) = half[start];
        let mut end = start + 1;
        while end < half.len() && half[end].0 == lo && half[end].1 == hi {
            end += 1;
        }
        if end - start >= 3 {
            found.push([lo, hi]);
        }
        start = end;
    }
    found
}

/// `validate.rs`'s pass 2: a vertex whose incident faces do not form one fan.
///
/// Two ways to fail, and both are counted there so both are counted here: three
/// or more incident faces sharing one wing vertex (an umbrella), or the incident
/// faces falling into more than one component under "shares a wing vertex" (a
/// bowtie, which the cheap faces-equals-edges test reports as clean).
fn non_manifold_vertices(faces: &[[u32; 3]]) -> Vec<u32> {
    let mut vertex_faces: Vec<(u32, u32)> = Vec::with_capacity(faces.len() * 3);
    for (fi, f) in faces.iter().enumerate() {
        for &v in f {
            vertex_faces.push((v, fi as u32));
        }
    }
    vertex_faces.sort_unstable();

    let mut found = Vec::new();
    let mut wings: Vec<(u32, u32)> = Vec::new();
    let mut link = Dsu::default();

    let mut start = 0usize;
    while start < vertex_faces.len() {
        let v = vertex_faces[start].0;
        let mut end = start + 1;
        while end < vertex_faces.len() && vertex_faces[end].0 == v {
            end += 1;
        }
        let degree = end - start;

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
            found.push(v);
        }

        start = end;
    }
    found
}

/// Both passes over one mesh.
fn recompute(vertex_count: usize, indices: &[u32]) -> Recomputed {
    let faces = usable_faces(vertex_count, indices);
    Recomputed {
        edges: non_manifold_edges(&faces),
        vertices: non_manifold_vertices(&faces),
    }
}

// ─── the arms ───────────────────────────────────────────────────────────────

/// Which dual extractor.
#[derive(Clone, Copy)]
enum Which {
    DualContouring,
    SurfaceNets,
}

impl Which {
    const fn name(self) -> &'static str {
        match self {
            Self::DualContouring => "dual_contouring",
            Self::SurfaceNets => "surface_nets",
        }
    }
}

/// The two dual extractors, in the order they are reported.
const ARMS: [Which; 2] = [Which::DualContouring, Which::SurfaceNets];

/// Everything one `(field, extractor)` pair produced.
struct Arm {
    field: &'static str,
    extractor: &'static str,
    cells: u64,
    critical_2d: u64,
    critical_3d: u64,
    critical: u64,
    mesh_vertices: u64,
    mesh_triangles: u64,
    report_edges: u64,
    report_vertices: u64,
    recomputed_edges: u64,
    recomputed_vertices: u64,
    locations_agree: bool,
    edge_incidents_critical: u64,
    vertex_incidents_critical: u64,
    incidents_both_ends: u64,
    active_cells: u64,
    critical_active_cells: u64,
    vertex_cells_unique: bool,
    escaped_cell: u64,
    /// Cells whose eight corner signs are not all equal: the cells that ought
    /// to host a dual vertex, one each.
    sign_active_cells: u64,
    /// `true` when the cells the vertex positions map back to are *exactly* the
    /// sign-active cells. False means at least one vertex is not in the cell
    /// that created it, which is the registration's rival suspect caught in the
    /// act rather than argued about.
    vertex_cells_match_sign_active: bool,
    /// Vertices whose mapped cell has no sign change at all, so it cannot be
    /// the cell that created them.
    vertices_in_inactive_cells: u64,
    /// Vertices sitting exactly on a cell boundary on some axis, where `floor`
    /// has to pick one of two cells. The explanation for
    /// `vertices_in_inactive_cells`, and the reason it is a mapping ambiguity
    /// rather than a vertex that escaped its cell.
    vertices_on_cell_boundary: u64,
    /// The same, restricted to the non-manifold vertices — the only ones whose
    /// attribution the co-location number depends on. Zero here means the
    /// ambiguity above cannot have moved a single incident.
    nm_vertices_on_cell_boundary: u64,
    /// Critical cells that host at least one non-manifold vertex. Equal to both
    /// `critical` and `recomputed_vertices` is a bijection rather than a
    /// coincidence, and that is worth being able to read off the row.
    critical_cells_with_nm_vertex: u64,
}

impl Arm {
    fn incidents(&self) -> u64 {
        self.recomputed_edges + self.recomputed_vertices
    }

    fn in_critical(&self) -> u64 {
        self.edge_incidents_critical + self.vertex_incidents_critical
    }
}

/// A fraction, or `n/a` when there is nothing to divide.
///
/// Never `1.0` for an empty numerator and denominator: a row with no incidents
/// has no co-location, and writing perfect agreement there would let the clean
/// fields outvote the dirty one in any later average.
fn fraction(part: u64, whole: u64) -> String {
    if whole == 0 {
        String::from("n/a")
    } else {
        format!("{:.6}", part as f64 / whole as f64)
    }
}

/// Extract, recompute, cross-tabulate against the census.
fn measure<F>(
    field_name: &'static str,
    field: &F,
    grid: &Grid,
    counted: &Census,
    which: Which,
) -> Arm
where
    F: Sdf<Scalar = f64>,
{
    let shape = RuntimeShape3::new([grid.samples; 3]).expect("census grid fits u32");
    let mut out = MeshBuffer::<f64>::new();
    match which {
        Which::DualContouring => DualContouring::<f64>::new()
            .extract(field, &shape, grid.origin, grid.cell_size, &mut out)
            .expect("dual contouring extraction"),
        Which::SurfaceNets => SurfaceNets::<f64>::new()
            .extract(field, &shape, grid.origin, grid.cell_size, &mut out)
            .expect("surface nets extraction"),
    }

    let cfg = ValidateConfig::from_cell_size(grid.cell_size).expect("positive cell size");
    let (report, features) = validate_features(&out.positions, &out.indices, &cfg);
    let mine = recompute(out.positions.len(), &out.indices);
    let locations_agree = mine.edges == features.edges && mine.vertices == features.vertices;

    // One vertex per cell is the property that makes this mapping meaningful, so
    // it is measured rather than assumed.
    let mut escaped_cell = 0u64;
    let cell_of_vertex: Vec<usize> = out
        .positions
        .iter()
        .map(|p| {
            let (cell, escaped) = grid.cell_of(*p);
            if escaped {
                escaped_cell += 1;
            }
            grid.cell_index(cell)
        })
        .collect();
    let mut distinct = cell_of_vertex.clone();
    distinct.sort_unstable();
    distinct.dedup();
    let critical_active_cells = distinct.iter().filter(|c| counted.critical[**c]).count() as u64;

    let vertex_incidents_critical = mine
        .vertices
        .iter()
        .filter(|v| counted.critical[cell_of_vertex[**v as usize]])
        .count() as u64;

    let mut edge_incidents_critical = 0u64;
    let mut edge_incidents_both = 0u64;
    for e in &mine.edges {
        let a = counted.critical[cell_of_vertex[e[0] as usize]];
        let b = counted.critical[cell_of_vertex[e[1] as usize]];
        if a || b {
            edge_incidents_critical += 1;
        }
        if a && b {
            edge_incidents_both += 1;
        }
    }

    // The cells the vertices landed in, against the cells that *created* them.
    // Equal sets mean every dual vertex is in its own cell, so the
    // position-to-cell map is a bijection onto the active grid.
    //
    // They are not always equal, and the reason is worth a column rather than a
    // guess. `Centroid` averages the cut points on a cell's twelve edges, and
    // with `smoothing_passes` at its default zero nothing clamps that average —
    // so on a field whose crossings sit exactly on cell corners (`box_exact` is
    // exactly zero across its whole boundary) the centroid can land exactly on
    // the cell's *upper* face, where `floor` names the neighbour. That is an
    // ambiguity in the mapping, not a vertex that ran away, and
    // `vertices_on_cell_boundary` is what tells the two apart. Dual Contouring
    // cannot hit it: `Clamp::ToCell` insets by `(1 − ε)` about the cell centre.
    let vertices_in_inactive_cells = cell_of_vertex
        .iter()
        .filter(|c| counted.active.binary_search(c).is_err())
        .count() as u64;
    let on_boundary = |p: &[f64; 3]| {
        (0..3usize).any(|axis| {
            let t = (p[axis] - grid.origin[axis]) / grid.cell_size;
            // `t - t.floor()` is never negative for a finite `t`, so `<= 0.0` is
            // "exactly integral" without an `==` on floats.
            t - t.floor() <= 0.0
        })
    };
    let vertices_on_cell_boundary = out.positions.iter().filter(|p| on_boundary(p)).count() as u64;
    let nm_vertices_on_cell_boundary = mine
        .vertices
        .iter()
        .filter(|v| on_boundary(&out.positions[**v as usize]))
        .count() as u64;

    // How many *distinct* critical cells host a non-manifold vertex. Together
    // with `critical` and `recomputed_vertices` this says whether the relation
    // is one-to-one.
    let mut nm_vertex_cells: Vec<usize> = mine
        .vertices
        .iter()
        .map(|v| cell_of_vertex[*v as usize])
        .filter(|c| counted.critical[*c])
        .collect();
    nm_vertex_cells.sort_unstable();
    nm_vertex_cells.dedup();

    Arm {
        field: field_name,
        extractor: which.name(),
        cells: counted.critical.len() as u64,
        critical_2d: counted.two_d,
        critical_3d: counted.three_d,
        critical: counted.either,
        mesh_vertices: out.vertex_count() as u64,
        mesh_triangles: out.triangle_count() as u64,
        report_edges: report.non_manifold_edges,
        report_vertices: report.non_manifold_vertices,
        recomputed_edges: mine.edges.len() as u64,
        recomputed_vertices: mine.vertices.len() as u64,
        locations_agree,
        edge_incidents_critical,
        // A vertex incident names one cell, so the union and both-ends rules
        // agree on it and only the edges differ.
        vertex_incidents_critical,
        incidents_both_ends: edge_incidents_both + vertex_incidents_critical,
        active_cells: distinct.len() as u64,
        critical_active_cells,
        vertex_cells_unique: distinct.len() == out.positions.len(),
        escaped_cell,
        sign_active_cells: counted.active.len() as u64,
        vertex_cells_match_sign_active: distinct == counted.active,
        vertices_in_inactive_cells,
        vertices_on_cell_boundary,
        nm_vertices_on_cell_boundary,
        critical_cells_with_nm_vertex: nm_vertex_cells.len() as u64,
    }
}

// ─── main ───────────────────────────────────────────────────────────────────

fn main() {
    let prereg = isomesh::experiment!("P-41");

    common::experiment::run(prereg, |run| {
        let table = classify();
        let two_d_bytes = table.two_d.iter().filter(|b| **b).count();
        let three_d_bytes = table.three_d.iter().filter(|b| **b).count();
        println!(
            "critical sign bytes, enumerated over all 256: \
             2d-critical {two_d_bytes}, 3d-critical {three_d_bytes}\n"
        );

        let mut arms: Vec<Arm> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            // No early return in here: the macro inlines its body once per
            // field, so a `return` would leave the sweep silently short.
            let (_, origin, cell_size) = common::grid::<f64, _>(&field, SAMPLES);
            let grid = Grid {
                origin,
                cell_size,
                samples: SAMPLES,
                cells: SAMPLES - 1,
            };
            let counted = census(&field, &grid, &table);
            println!(
                "{name:>14}  cells {:>7}  2d-critical {:>7}  3d-critical {:>7}  \
                 either {:>7} ({:.4}%)",
                counted.critical.len(),
                counted.two_d,
                counted.three_d,
                counted.either,
                100.0 * counted.either as f64 / counted.critical.len() as f64,
            );
            for which in ARMS {
                let arm = measure(name, &field, &grid, &counted, which);
                println!(
                    "{:>14}  {:>16}  nm-edges {:>5}  nm-vertices {:>5}  \
                     in-critical {:>5}  ({})  null {}",
                    "",
                    arm.extractor,
                    arm.recomputed_edges,
                    arm.recomputed_vertices,
                    arm.in_critical(),
                    fraction(arm.in_critical(), arm.incidents()),
                    fraction(arm.critical_active_cells, arm.active_cells),
                );
                arms.push(arm);
            }
        });

        // Clause two is a statement about the incidents, not about the fields,
        // so it is decided pooled. Per extractor and over both, repeated on
        // every row so the CSV can be read one line at a time.
        let overall_incidents: u64 = arms.iter().map(Arm::incidents).sum();
        let overall_critical: u64 = arms.iter().map(Arm::in_critical).sum();
        let overall_both: u64 = arms.iter().map(|a| a.incidents_both_ends).sum();

        println!(
            "\npooled over both extractors: {overall_critical}/{overall_incidents} incidents in \
             critical cells ({}), both-ends {}",
            fraction(overall_critical, overall_incidents),
            fraction(overall_both, overall_incidents),
        );

        for which in ARMS {
            let pooled = || arms.iter().filter(|a| a.extractor == which.name());
            let incidents: u64 = pooled().map(Arm::incidents).sum();
            let critical: u64 = pooled().map(Arm::in_critical).sum();
            println!(
                "{:>16}: {critical}/{incidents} ({})",
                which.name(),
                fraction(critical, incidents),
            );
        }

        for arm in &arms {
            let pooled = || arms.iter().filter(|a| a.extractor == arm.extractor);
            let pooled_incidents: u64 = pooled().map(Arm::incidents).sum();
            let pooled_critical: u64 = pooled().map(Arm::in_critical).sum();

            run.record(&[
                ("field", arm.field.to_string()),
                ("samples_per_axis", SAMPLES.to_string()),
                ("cells", arm.cells.to_string()),
                ("critical_2d_cells", arm.critical_2d.to_string()),
                ("critical_3d_cells", arm.critical_3d.to_string()),
                ("extractor", arm.extractor.to_string()),
                ("non_manifold_edges", arm.report_edges.to_string()),
                ("non_manifold_vertices", arm.report_vertices.to_string()),
                ("incidents_in_critical_cells", arm.in_critical().to_string()),
                (
                    "colocation_fraction",
                    fraction(arm.in_critical(), arm.incidents()),
                ),
                ("critical_cells", arm.critical.to_string()),
                ("critical_cell_fraction", fraction(arm.critical, arm.cells)),
                ("incidents_total", arm.incidents().to_string()),
                (
                    "edge_incidents_in_critical",
                    arm.edge_incidents_critical.to_string(),
                ),
                (
                    "vertex_incidents_in_critical",
                    arm.vertex_incidents_critical.to_string(),
                ),
                (
                    "incidents_in_critical_both_ends",
                    arm.incidents_both_ends.to_string(),
                ),
                (
                    "colocation_fraction_both_ends",
                    fraction(arm.incidents_both_ends, arm.incidents()),
                ),
                (
                    "null_colocation_fraction",
                    fraction(arm.critical_active_cells, arm.active_cells),
                ),
                ("active_cells", arm.active_cells.to_string()),
                ("sign_active_cells", arm.sign_active_cells.to_string()),
                (
                    "vertex_cells_match_sign_active",
                    arm.vertex_cells_match_sign_active.to_string(),
                ),
                (
                    "vertices_in_inactive_cells",
                    arm.vertices_in_inactive_cells.to_string(),
                ),
                (
                    "vertices_on_cell_boundary",
                    arm.vertices_on_cell_boundary.to_string(),
                ),
                (
                    "nm_vertices_on_cell_boundary",
                    arm.nm_vertices_on_cell_boundary.to_string(),
                ),
                (
                    "critical_cells_with_nm_vertex",
                    arm.critical_cells_with_nm_vertex.to_string(),
                ),
                (
                    "critical_active_cells",
                    arm.critical_active_cells.to_string(),
                ),
                (
                    "recomputed_non_manifold_edges",
                    arm.recomputed_edges.to_string(),
                ),
                (
                    "recomputed_non_manifold_vertices",
                    arm.recomputed_vertices.to_string(),
                ),
                (
                    "counts_agree",
                    (arm.report_edges == arm.recomputed_edges
                        && arm.report_vertices == arm.recomputed_vertices)
                        .to_string(),
                ),
                ("locations_agree", arm.locations_agree.to_string()),
                ("vertex_cells_unique", arm.vertex_cells_unique.to_string()),
                ("vertices_escaping_grid", arm.escaped_cell.to_string()),
                ("mesh_vertices", arm.mesh_vertices.to_string()),
                ("mesh_triangles", arm.mesh_triangles.to_string()),
                ("pooled_incidents", pooled_incidents.to_string()),
                ("pooled_in_critical_cells", pooled_critical.to_string()),
                (
                    "pooled_colocation_fraction",
                    fraction(pooled_critical, pooled_incidents),
                ),
                ("overall_incidents", overall_incidents.to_string()),
                ("overall_in_critical_cells", overall_critical.to_string()),
                (
                    "overall_colocation_fraction",
                    fraction(overall_critical, overall_incidents),
                ),
            ]);
        }
    });
}
