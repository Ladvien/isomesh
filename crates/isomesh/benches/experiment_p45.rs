//! **P-45 — additivity: does a chunk-local curvature measure compose?**
//!
//! Ticket: R-041a. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p45
//! ```
//!
//! Writes `docs/experiments/p-45.csv`.
//!
//! # What P-42 could not reach
//!
//! P-42 took `B` to be the whole closed surface, which made `∂B` empty on every
//! row and its Gaussian clause an identity: `3F = 2E`, so `Σ_v (2π − α_v) =
//! 2πχ` combinatorially and the residual was one f64 epsilon per vertex. The
//! property the crate actually wants was never touched — **additivity**,
//! `N(A ∪ B) = N(A) + N(B) − N(A ∩ B)`, the law that lets a per-chunk number
//! compose into a per-world one with no global pass.
//!
//! So here every `B` is a genuine patch: the mesh's triangles are partitioned
//! into a 4×4×4 grid of spatial chunks by centroid, and each chunk carries the
//! geodesic-curvature boundary term that an open patch needs.
//!
//! # The formulas, transcribed rather than remembered
//!
//! Sun & Morvan, *Curvature measures, normal cycles and asymptotic cones*,
//! `10.5802/acirm.50`:
//!
//! - **Theorem 5 (2)(b)** — `Φ_P^H(B) = Σ_{e ∈ E ∩ B} l(e ∩ B) β(e)`, with
//!   `β(e) ∈ [−π, π]` the angle between the normals of the two triangles
//!   incident on `e`, **positive if `e` is convex and negative if concave**.
//!   That sign rule is why the dihedral is computed with `atan2` against the
//!   edge direction and not with `acos`.
//! - **§3, (3.2)** — the additivity law itself, `N(A ∪ B) = N(A) + N(B) −
//!   N(A ∩ B)`. Note the third term: it is what this experiment turns out to be
//!   about.
//! - **Theorem 5 (2)(a)** prints the polyhedral Gaussian measure as
//!   `Σ_{v ∈ V ∩ B} α_v`, where Definition 1 (1) has just defined `α_v` as the
//!   *sum* of incident corner angles. That total sits near `2π|V|`, not near
//!   `∫ G da`; the quantity the theory produces is the **defect** `2π − α_v`,
//!   and this file uses the defect. The same printed discrepancy was recorded in
//!   P-42 and is carried forward here unchanged.
//! - **Definition 3 / Theorem 2 (3)** equate `Σ l(e) ∠(e)` with `∫_S H da` under
//!   `H = ½ trace A` (§1.1). Those differ by a factor of two and the cube
//!   settles it: for `[−1, 1]³` the Steiner `ε²` coefficient is twelve quarter
//!   cylinders, `12 · (π/4) · 2 = 6π`, while `Σ l(e) ∠(e) = 12 · 2 · (π/2) =
//!   12π`. So `∫_S H da = ½ Σ l(e) β(e)`. `mean_global` and `mean_chunk_sum` are
//!   the literal `Σ l(e) β(e)` the registration names; `mean_curvature_half` and
//!   `mean_smooth_integral_h` sit beside them, which is where the `box_exact`
//!   check lands: a grid-aligned `[−1, 1]³` box has twelve edges of length two
//!   at exactly `π/2`, so `½ Σ l β` should be `6π` **exactly** rather than
//!   convergently.
//!
//! # The boundary term, as registered
//!
//! A vertex of a patch is a *boundary vertex* when some edge of the patch has
//! exactly one patch face on it. The registered rule is `π − α_v` there and
//! `2π − α_v` at an interior vertex — the standard discrete Gauss–Bonnet for a
//! disc, whose interior-plus-boundary sum is `2πχ(disc) = 2π`. Edges on a chunk
//! boundary are weighted one half, so an edge split between two chunks is paid
//! for once in total.
//!
//! # Two things the caller arranges, not the measure
//!
//! **The chunk grid spans the mesh's own bounding box, not the field's domain.**
//! The registration says "a 4×4×4 grid of spatial chunks" without saying over
//! what, and over the domain it would be a poor test: the canonical `[-2, 2]³`
//! domain divided in four puts the unit sphere inside eight of the sixty-four
//! boxes, so most chunks would be empty and most of the surface would never see
//! a seam. Over the bounding box every box is a real patch. `chunk_boxes` and
//! `chunks` are both recorded so the count is never ambiguous.
//!
//! **The isolated and in-context arms differ in exactly one input.** Clause
//! three asks whether a chunk recomputed from its own triangles alone reproduces
//! its in-context value. Both arms therefore traverse the *same* patch, built
//! with its vertices in ascending global-index order so the patch-local order is
//! the global relative order and no summation order changes between them. The
//! only difference is that the in-context accumulator may consult the global
//! dihedral table for a seam edge and the isolated one may not — because a patch
//! holding one of an edge's two faces has no second normal to turn against. A
//! bit difference can therefore only come from a seam term, which is what makes
//! the comparison a measurement of locality rather than of float ordering. The
//! Gaussian arm is checked harder still: its in-context value is recomputed by a
//! separate routine that never builds a patch at all, gathering angle sums from
//! the global face list and boundary status from global edge incidence.
//!
//! # Nothing here is timed
//!
//! Every column is a deterministic function of the extracted mesh. Accumulation
//! orders are fixed on purpose — vertices in index order, edges in sorted-key
//! order, chunks in box order — so the float sums reproduce run to run.

mod common;

use common::experiment::Run;
use isomesh::fields::{BoxExact, ReferenceField, Sphere, Torus};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{ValidateConfig, validate};
use isomesh::{MeshBuffer, Sdf};

/// The registered resolution.
const SAMPLES: u32 = 65;

/// Chunks per axis. `4 × 4 × 4`, as registered.
const SIDE: usize = 4;

/// Boxes in the grid, empty ones included.
const BOXES: usize = SIDE * SIDE * SIDE;

const TAU: f64 = std::f64::consts::TAU;
const PI: f64 = std::f64::consts::PI;

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Unit vector, and `[0, 0, 0]` for the zero vector.
///
/// A triangle of exactly zero area has no normal, so its dihedral contribution
/// is `atan2(0, 0) = 0` rather than a NaN. Degenerate triangles are counted in
/// `degenerate_triangles` rather than hidden by this.
fn unit(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    if n > 0.0 {
        [a[0] / n, a[1] / n, a[2] / n]
    } else {
        [0.0; 3]
    }
}

/// The unsigned angle at `apex` between the rays to `u` and `v`.
///
/// `atan2(|a × b|, a · b)` rather than `acos(â · b̂)`: the same value, but
/// conditioned for the sliver angles marching cubes produces, where the dot
/// product of two nearly parallel unit vectors rounds to exactly one.
fn corner_angle(apex: [f64; 3], u: [f64; 3], v: [f64; 3]) -> f64 {
    let a = sub(u, apex);
    let b = sub(v, apex);
    norm(cross(a, b)).atan2(dot(a, b))
}

/// One triangle's geometry.
struct Face {
    /// Vertex indices as the input gave them — the traversal order the dihedral
    /// sign depends on.
    verts: [u32; 3],
    /// Outward unit normal, from the counter-clockwise-seen-from-outside winding
    /// this crate guarantees.
    normal: [f64; 3],
    centroid: [f64; 3],
    /// `A(t) <= area_epsilon_rel · cell_size²`, the crate's own definition.
    degenerate: bool,
}

/// One `Face`, from three positions.
///
/// The single constructor, used for the global mesh **and** for every patch, so
/// a patch's normals are bit-identical to the global mesh's: same positions,
/// same expression, same order.
fn make_face(verts: [u32; 3], pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], area_floor: f64) -> Face {
    let twice = cross(sub(pb, pa), sub(pc, pa));
    Face {
        verts,
        normal: unit(twice),
        centroid: [
            (pa[0] + pb[0] + pc[0]) / 3.0,
            (pa[1] + pb[1] + pc[1]) / 3.0,
            (pa[2] + pb[2] + pc[2]) / 3.0,
        ],
        degenerate: norm(twice) / 2.0 <= area_floor,
    }
}

/// Faces with three distinct in-range indices, and the count of the rest.
///
/// The predicate is the one [`isomesh::validate`] uses for `faces_skipped`, so
/// this accumulator and `MeshReport` describe the same face set.
fn faces_from(positions: &[[f64; 3]], indices: &[u32], area_floor: f64) -> (Vec<Face>, u64) {
    let nv = positions.len();
    let mut faces = Vec::with_capacity(indices.len() / 3);
    let mut skipped = 0u64;
    for tri in indices.chunks_exact(3) {
        let (a, b, c) = (tri[0], tri[1], tri[2]);
        let in_range = (a as usize) < nv && (b as usize) < nv && (c as usize) < nv;
        if !in_range || a == b || b == c || c == a {
            skipped += 1;
            continue;
        }
        faces.push(make_face(
            [a, b, c],
            positions[a as usize],
            positions[b as usize],
            positions[c as usize],
            area_floor,
        ));
    }
    (faces, skipped)
}

/// An occurrence of an undirected edge in one face.
struct EdgeRef {
    lo: u32,
    hi: u32,
    face: u32,
    /// `true` when this face traverses the edge `lo → hi`.
    forward: bool,
}

/// Every `EdgeRef` of every face, sorted so grouping is a linear scan and every
/// sum accumulates in a fixed order.
fn sorted_edge_refs(faces: &[Face]) -> Vec<EdgeRef> {
    let mut refs = Vec::with_capacity(faces.len() * 3);
    for (fi, face) in faces.iter().enumerate() {
        let [a, b, c] = face.verts;
        for (u, v) in [(a, b), (b, c), (c, a)] {
            refs.push(EdgeRef {
                lo: u.min(v),
                hi: u.max(v),
                face: fi as u32,
                forward: u < v,
            });
        }
    }
    refs.sort_unstable_by_key(|e| (e.lo, e.hi, e.face));
    refs
}

/// One past the end of the run of `refs` sharing `refs[i]`'s endpoints.
fn group_end(refs: &[EdgeRef], i: usize) -> usize {
    let mut j = i + 1;
    while j < refs.len() && refs[j].lo == refs[i].lo && refs[j].hi == refs[i].hi {
        j += 1;
    }
    j
}

/// The signed dihedral of Theorem 5 (2)(b): `+` convex, `−` concave.
///
/// `atan2((n₁ × n₂) · ê, n₁ · n₂)` with `ê` the edge direction **as the first
/// face traverses it**. Swapping which face is first flips `ê` and `n₁ × n₂`
/// together, so the value does not depend on the order the group was sorted
/// into. Checked against the cube: top face `n₁ = +z` meets side face `n₂ = +y`
/// along an edge the top traverses in `−x`, giving `atan2(1, 0) = +π/2` —
/// convex, positive, as the sign rule requires.
fn signed_dihedral(first: &Face, second: &Face, edge_dir: [f64; 3]) -> f64 {
    let s = dot(cross(first.normal, second.normal), edge_dir);
    let c = dot(first.normal, second.normal);
    s.atan2(c)
}

/// The edge direction as `refs[k]`'s face traverses it, unit length.
fn traversal_dir(positions: &[[f64; 3]], e: &EdgeRef) -> [f64; 3] {
    let (lo, hi) = (positions[e.lo as usize], positions[e.hi as usize]);
    if e.forward {
        unit(sub(hi, lo))
    } else {
        unit(sub(lo, hi))
    }
}

/// An undirected edge of the global mesh, carrying the dihedral only the global
/// mesh can compute. The table is sorted by `(lo, hi)`, so a lookup is a binary
/// search.
struct GlobalEdge {
    lo: u32,
    hi: u32,
    length: f64,
    /// `β(e)`. NaN when the edge does not have exactly two incident faces, so
    /// that an anomalous mesh poisons the total visibly instead of quietly.
    beta: f64,
    /// The first two incident faces in sorted order.
    faces: [u32; 2],
    face_count: u32,
}

/// The whole mesh: its faces, its edge table, and both global measures.
struct GlobalMesh {
    faces: Vec<Face>,
    edges: Vec<GlobalEdge>,
    /// `Σ_v (2π − α_v)`, vertices in index order.
    gaussian_defect_sum: f64,
    /// `Σ_e l(e) β(e)`, edges in sorted-key order.
    mean: f64,
    /// `Σ l β` split by dihedral class, so a piecewise-planar surface can be
    /// read: `flat` is `|β| <= 1e-12`, `right` is `||β| − π/2| <= 1e-9`, `other`
    /// is everything left. On `box_exact` the whole measure should live in
    /// `right`, twelve edges of total length twenty-four at exactly `π/2`, and
    /// what leaks into `other` is what the extractor did to the corners.
    mean_flat: f64,
    mean_right: f64,
    mean_other: f64,
    flat_edges: u64,
    right_edges: u64,
    other_edges: u64,
    /// `Σ l` over the `right` class. `24` for an exact `[−1, 1]³`.
    right_length: f64,
    vertices: usize,
    referenced_vertices: u64,
    skipped_faces: u64,
    degenerate_faces: u64,
    boundary_edges: u64,
    non_manifold_edges: u64,
}

fn analyse_global(mesh: &MeshBuffer<f64>, area_floor: f64) -> GlobalMesh {
    let (faces, skipped_faces) = faces_from(&mesh.positions, &mesh.indices, area_floor);
    let nv = mesh.positions.len();
    let mut angle_sum = vec![0.0_f64; nv];
    let mut referenced = vec![false; nv];
    for face in &faces {
        let [a, b, c] = face.verts;
        let (pa, pb, pc) = (
            mesh.positions[a as usize],
            mesh.positions[b as usize],
            mesh.positions[c as usize],
        );
        for (v, apex, u, w) in [(a, pa, pb, pc), (b, pb, pc, pa), (c, pc, pa, pb)] {
            angle_sum[v as usize] += corner_angle(apex, u, w);
            referenced[v as usize] = true;
        }
    }

    let refs = sorted_edge_refs(&faces);
    let mut edges = Vec::with_capacity(refs.len() / 2);
    let mut mean = 0.0_f64;
    let mut classes = EdgeClasses::default();
    let mut boundary_edges = 0u64;
    let mut non_manifold_edges = 0u64;
    let mut i = 0usize;
    while i < refs.len() {
        let j = group_end(&refs, i);
        let group = &refs[i..j];
        let (lo, hi) = (group[0].lo, group[0].hi);
        let length = norm(sub(
            mesh.positions[hi as usize],
            mesh.positions[lo as usize],
        ));
        let beta = if group.len() == 2 {
            signed_dihedral(
                &faces[group[0].face as usize],
                &faces[group[1].face as usize],
                traversal_dir(&mesh.positions, &group[0]),
            )
        } else {
            if group.len() == 1 {
                boundary_edges += 1;
            } else {
                non_manifold_edges += 1;
            }
            f64::NAN
        };
        if group.len() == 2 {
            mean += length * beta;
            classes.add(length, beta);
        }
        edges.push(GlobalEdge {
            lo,
            hi,
            length,
            beta,
            faces: [group[0].face, group[group.len().min(2) - 1].face],
            face_count: group.len() as u32,
        });
        i = j;
    }

    let mut gaussian_defect_sum = 0.0_f64;
    let mut referenced_vertices = 0u64;
    for (&sum, &seen) in angle_sum.iter().zip(&referenced) {
        if seen {
            gaussian_defect_sum += TAU - sum;
            referenced_vertices += 1;
        }
    }

    GlobalMesh {
        degenerate_faces: faces.iter().filter(|f| f.degenerate).count() as u64,
        faces,
        edges,
        gaussian_defect_sum,
        mean,
        mean_flat: classes.mean_flat,
        mean_right: classes.mean_right,
        mean_other: classes.mean_other,
        flat_edges: classes.flat_edges,
        right_edges: classes.right_edges,
        other_edges: classes.other_edges,
        right_length: classes.right_length,
        vertices: nv,
        referenced_vertices,
        skipped_faces,
        boundary_edges,
        non_manifold_edges,
    }
}

/// `Σ l β` split three ways by dihedral class.
///
/// A grid-aligned box should have every non-zero term in `right`: flat faces
/// contribute `β = 0` exactly, and the twelve real edges contribute exactly
/// `π/2`. Anything in `other` is turning the extractor put somewhere the
/// geometry does not have it.
#[derive(Default)]
struct EdgeClasses {
    mean_flat: f64,
    mean_right: f64,
    mean_other: f64,
    flat_edges: u64,
    right_edges: u64,
    other_edges: u64,
    right_length: f64,
}

impl EdgeClasses {
    fn add(&mut self, length: f64, beta: f64) {
        let term = length * beta;
        if beta.abs() <= 1e-12 {
            self.mean_flat += term;
            self.flat_edges += 1;
        } else if (beta.abs() - PI / 2.0).abs() <= 1e-9 {
            self.mean_right += term;
            self.right_edges += 1;
            self.right_length += length;
        } else {
            self.mean_other += term;
            self.other_edges += 1;
        }
    }
}

/// `β(e)` for a global edge named by its endpoints.
///
/// The only thing an isolated patch is not allowed to call.
fn global_beta(edges: &[GlobalEdge], lo: u32, hi: u32) -> f64 {
    match edges.binary_search_by(|e| (e.lo, e.hi).cmp(&(lo, hi))) {
        Ok(k) => edges[k].beta,
        // An edge of a patch that is not an edge of the mesh the patch came from
        // is not a rounding problem, so it does not get a rounding-sized answer.
        Err(_) => f64::NAN,
    }
}

/// One chunk, as a standalone mesh.
///
/// `global` is ascending, so patch-local index order **is** the global relative
/// order and no sum changes order when a chunk is lifted out of its mesh.
struct Patch {
    global: Vec<u32>,
    positions: Vec<[f64; 3]>,
    indices: Vec<u32>,
}

fn build_patch(mesh: &MeshBuffer<f64>, faces: &[Face], chunk_of_face: &[u8], chunk: u8) -> Patch {
    let mut global: Vec<u32> = Vec::new();
    for (fi, face) in faces.iter().enumerate() {
        if chunk_of_face[fi] == chunk {
            global.extend_from_slice(&face.verts);
        }
    }
    global.sort_unstable();
    global.dedup();

    let positions = global
        .iter()
        .map(|&v| mesh.positions[v as usize])
        .collect::<Vec<_>>();
    let mut indices = Vec::new();
    for (fi, face) in faces.iter().enumerate() {
        if chunk_of_face[fi] != chunk {
            continue;
        }
        for v in face.verts {
            let local = global
                .binary_search(&v)
                .expect("a patch vertex is in the patch's own vertex list");
            indices.push(local as u32);
        }
    }
    Patch {
        global,
        positions,
        indices,
    }
}

/// One chunk's measures, both arms.
struct PatchOut {
    /// `Σ_interior (2π − α_v) + Σ_boundary (π − α_v)`, from the patch alone.
    gaussian: f64,
    /// `Σ_e w(e) l(e) β(e)` with the seam terms **omitted**, because a patch
    /// holding one of an edge's two faces has no second normal to turn against.
    mean_isolated: f64,
    /// The same traversal, with each seam edge contributing `½ l(e) β_global(e)`.
    mean_in_context: f64,
    interior_vertices: u64,
    boundary_vertices: u64,
    interior_edges: u64,
    seam_edges: u64,
    faces: u64,
}

fn patch_measures(patch: &Patch, global_edges: &[GlobalEdge], area_floor: f64) -> PatchOut {
    let (faces, _) = faces_from(&patch.positions, &patch.indices, area_floor);
    let nv = patch.positions.len();
    let mut angle_sum = vec![0.0_f64; nv];
    let mut referenced = vec![false; nv];
    for face in &faces {
        let [a, b, c] = face.verts;
        let (pa, pb, pc) = (
            patch.positions[a as usize],
            patch.positions[b as usize],
            patch.positions[c as usize],
        );
        for (v, apex, u, w) in [(a, pa, pb, pc), (b, pb, pc, pa), (c, pc, pa, pb)] {
            angle_sum[v as usize] += corner_angle(apex, u, w);
            referenced[v as usize] = true;
        }
    }

    let refs = sorted_edge_refs(&faces);
    let mut on_boundary = vec![false; nv];
    let mut mean_isolated = 0.0_f64;
    let mut mean_in_context = 0.0_f64;
    let mut interior_edges = 0u64;
    let mut seam_edges = 0u64;
    let mut i = 0usize;
    while i < refs.len() {
        let j = group_end(&refs, i);
        let group = &refs[i..j];
        let (lo, hi) = (group[0].lo as usize, group[0].hi as usize);
        let length = norm(sub(patch.positions[hi], patch.positions[lo]));
        if group.len() == 2 {
            let beta = signed_dihedral(
                &faces[group[0].face as usize],
                &faces[group[1].face as usize],
                traversal_dir(&patch.positions, &group[0]),
            );
            mean_isolated += length * beta;
            mean_in_context += length * beta;
            interior_edges += 1;
        } else {
            // Exactly one patch face on this edge: a chunk seam. The isolated
            // arm gets nothing here — that omission is the measurement.
            seam_edges += 1;
            on_boundary[lo] = true;
            on_boundary[hi] = true;
            let beta = global_beta(global_edges, patch.global[lo], patch.global[hi]);
            mean_in_context += 0.5 * length * beta;
        }
        i = j;
    }

    let mut gaussian = 0.0_f64;
    let mut interior_vertices = 0u64;
    let mut boundary_vertices = 0u64;
    for ((&sum, &seen), &edge) in angle_sum.iter().zip(&referenced).zip(&on_boundary) {
        if !seen {
            continue;
        }
        if edge {
            gaussian += PI - sum;
            boundary_vertices += 1;
        } else {
            gaussian += TAU - sum;
            interior_vertices += 1;
        }
    }

    PatchOut {
        gaussian,
        mean_isolated,
        mean_in_context,
        interior_vertices,
        boundary_vertices,
        interior_edges,
        seam_edges,
        faces: faces.len() as u64,
    }
}

/// One chunk's Gaussian measure computed **without building a patch**.
///
/// Angle sums are gathered from the global face list filtered to `chunk`, and
/// boundary status from global edge incidence: an edge straddling the chunk
/// boundary makes both its endpoints boundary vertices of the chunk. This exists
/// so clause three's Gaussian arm compares two independent routines rather than
/// a value against itself.
fn in_context_gaussian(
    g: &GlobalMesh,
    chunk_of_face: &[u8],
    chunk: u8,
    positions: &[[f64; 3]],
) -> f64 {
    let nv = g.vertices;
    let mut angle_sum = vec![0.0_f64; nv];
    let mut present = vec![false; nv];
    for (fi, face) in g.faces.iter().enumerate() {
        if chunk_of_face[fi] != chunk {
            continue;
        }
        let [a, b, c] = face.verts;
        let (pa, pb, pc) = (
            positions[a as usize],
            positions[b as usize],
            positions[c as usize],
        );
        for (v, apex, u, w) in [(a, pa, pb, pc), (b, pb, pc, pa), (c, pc, pa, pb)] {
            angle_sum[v as usize] += corner_angle(apex, u, w);
            present[v as usize] = true;
        }
    }

    let mut on_boundary = vec![false; nv];
    for e in &g.edges {
        let inside_first = chunk_of_face[e.faces[0] as usize] == chunk;
        let inside_second = chunk_of_face[e.faces[1] as usize] == chunk;
        if inside_first != inside_second {
            on_boundary[e.lo as usize] = true;
            on_boundary[e.hi as usize] = true;
        }
    }

    let mut gaussian = 0.0_f64;
    for ((&sum, &seen), &edge) in angle_sum.iter().zip(&present).zip(&on_boundary) {
        if seen {
            gaussian += if edge { PI - sum } else { TAU - sum };
        }
    }
    gaussian
}

/// The `4 × 4 × 4` box a point falls in, `x` slowest.
///
/// Half-open per axis with the top box closed, so every point of the mesh lands
/// in exactly one box and the boxes are a partition.
fn box_of(p: [f64; 3], lo: [f64; 3], inv_width: [f64; 3]) -> usize {
    let mut idx = 0usize;
    for (v, (l, iw)) in p.iter().zip(lo.iter().zip(inv_width.iter())) {
        let t = ((v - l) * iw).floor();
        let i = if t < 0.0 {
            0
        } else if t > (SIDE - 1) as f64 {
            SIDE - 1
        } else {
            t as usize
        };
        idx = idx * SIDE + i;
    }
    idx
}

/// One chunk's pair of arms, kept so the aggregate and the per-chunk comparison
/// read the same numbers.
struct ChunkRow {
    gaussian_patch: f64,
    gaussian_in_context: f64,
    mean_isolated: f64,
    mean_in_context: f64,
    seam_edges: u64,
    interior_edges: u64,
    boundary_vertices: u64,
    interior_vertices: u64,
    faces: u64,
}

/// Everything the row needs that is not a single global number.
struct Aggregate {
    chunks: usize,
    gaussian_chunk_sum: f64,
    mean_chunk_sum: f64,
    mean_isolated_sum: f64,
    gaussian_bit_identical: bool,
    mean_bit_identical: bool,
    gaussian_mismatches: u64,
    mean_mismatches: u64,
    seam_edges: u64,
    interior_edges: u64,
    boundary_vertices: u64,
    interior_vertices: u64,
    max_faces: u64,
    min_faces: u64,
}

fn aggregate(rows: &[ChunkRow]) -> Aggregate {
    let mut out = Aggregate {
        chunks: rows.len(),
        gaussian_chunk_sum: 0.0,
        mean_chunk_sum: 0.0,
        mean_isolated_sum: 0.0,
        gaussian_bit_identical: true,
        mean_bit_identical: true,
        gaussian_mismatches: 0,
        mean_mismatches: 0,
        seam_edges: 0,
        interior_edges: 0,
        boundary_vertices: 0,
        interior_vertices: 0,
        max_faces: 0,
        min_faces: u64::MAX,
    };
    for r in rows {
        out.gaussian_chunk_sum += r.gaussian_patch;
        out.mean_chunk_sum += r.mean_in_context;
        out.mean_isolated_sum += r.mean_isolated;
        out.seam_edges += r.seam_edges;
        out.interior_edges += r.interior_edges;
        out.boundary_vertices += r.boundary_vertices;
        out.interior_vertices += r.interior_vertices;
        out.max_faces = out.max_faces.max(r.faces);
        out.min_faces = out.min_faces.min(r.faces);
        if r.gaussian_patch.to_bits() != r.gaussian_in_context.to_bits() {
            out.gaussian_mismatches += 1;
            out.gaussian_bit_identical = false;
        }
        if r.mean_isolated.to_bits() != r.mean_in_context.to_bits() {
            out.mean_mismatches += 1;
            out.mean_bit_identical = false;
        }
    }
    out
}

/// How many chunks each vertex has a face in, and the excess over two.
///
/// `Σ_v max(k(v) − 2, 0)` is the diagnostic for clause one: the registered rule
/// pays `π` once per chunk that a boundary vertex appears in, so a vertex shared
/// by `k` chunks is charged `kπ` where `2π` is owed.
struct Sharing {
    excess: i64,
    k_hist: [u64; 5],
}

fn sharing(faces: &[Face], chunk_of_face: &[u8], nv: usize) -> Sharing {
    let mut last = vec![u8::MAX; nv];
    let mut k = vec![0u32; nv];
    for chunk in 0..BOXES {
        let c = chunk as u8;
        for (fi, face) in faces.iter().enumerate() {
            if chunk_of_face[fi] != c {
                continue;
            }
            for v in face.verts {
                if last[v as usize] != c {
                    last[v as usize] = c;
                    k[v as usize] += 1;
                }
            }
        }
    }
    let mut out = Sharing {
        excess: 0,
        k_hist: [0; 5],
    };
    for &kv in &k {
        if kv == 0 {
            continue;
        }
        out.excess += i64::from(kv.max(2) - 2);
        out.k_hist[(kv as usize).min(4)] += 1;
    }
    out
}

/// Each global term assigned whole to the box holding its vertex, or its edge's
/// midpoint.
///
/// This is the *measure-theoretic* decomposition, and it is exact by
/// construction: every term is assigned exactly once, so the boxes' totals sum
/// to the global totals up to summation order. It is recorded because that
/// exactness is the point — `α_v` here is the sum over **all** faces at `v` and
/// `β(e)` needs **both** faces of `e`, so every term is global data wearing a
/// chunk's label. The paper's `l(e ∩ B)` would clip an edge across a box face
/// rather than assign it whole; that changes how the total is distributed, not
/// that the distribution is exact.
fn borel_sums(
    g: &GlobalMesh,
    positions: &[[f64; 3]],
    lo: [f64; 3],
    inv_width: [f64; 3],
) -> (f64, f64) {
    let mut boxes_g = [0.0_f64; BOXES];
    let mut boxes_h = [0.0_f64; BOXES];

    let mut angle_sum = vec![0.0_f64; g.vertices];
    let mut referenced = vec![false; g.vertices];
    for face in &g.faces {
        let [a, b, c] = face.verts;
        let (pa, pb, pc) = (
            positions[a as usize],
            positions[b as usize],
            positions[c as usize],
        );
        for (v, apex, u, w) in [(a, pa, pb, pc), (b, pb, pc, pa), (c, pc, pa, pb)] {
            angle_sum[v as usize] += corner_angle(apex, u, w);
            referenced[v as usize] = true;
        }
    }
    for (v, (&sum, &seen)) in angle_sum.iter().zip(&referenced).enumerate() {
        if seen {
            boxes_g[box_of(positions[v], lo, inv_width)] += TAU - sum;
        }
    }
    for e in &g.edges {
        if e.face_count != 2 {
            continue;
        }
        let (p, q) = (positions[e.lo as usize], positions[e.hi as usize]);
        let mid = [
            (p[0] + q[0]) / 2.0,
            (p[1] + q[1]) / 2.0,
            (p[2] + q[2]) / 2.0,
        ];
        boxes_h[box_of(mid, lo, inv_width)] += e.length * e.beta;
    }
    (boxes_g.iter().sum(), boxes_h.iter().sum())
}

/// One field: extract, partition, and record.
fn sweep<F>(run: &mut Run, field: &F, smooth_mean_h: f64)
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, cell) = common::grid(field, SAMPLES);
    let mut mc = MarchingCubes::new();
    let mut out = MeshBuffer::new();
    mc.extract(field, &shape, origin, cell, &mut out)
        .expect("marching cubes on a closed reference field");

    let cfg = ValidateConfig::from_cell_size(cell).expect("cell size is positive");
    let report = validate(&out, &cfg);
    let area_floor = cfg.area_epsilon_rel() * cfg.cell_size() * cfg.cell_size();
    let g = analyse_global(&out, area_floor);

    // The chunk grid spans the mesh's bounding box, so every box is a real
    // patch. See the module docs for why not the field's domain.
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in &out.positions {
        for (a, &v) in p.iter().enumerate() {
            lo[a] = lo[a].min(v);
            hi[a] = hi[a].max(v);
        }
    }
    let mut inv_width = [0.0_f64; 3];
    for (a, iw) in inv_width.iter_mut().enumerate() {
        *iw = (SIDE as f64) / (hi[a] - lo[a]);
    }

    let chunk_of_face: Vec<u8> = g
        .faces
        .iter()
        .map(|f| box_of(f.centroid, lo, inv_width) as u8)
        .collect();

    let mut rows: Vec<ChunkRow> = Vec::new();
    for chunk in 0..BOXES {
        let c = chunk as u8;
        if !chunk_of_face.contains(&c) {
            continue;
        }
        let patch = build_patch(&out, &g.faces, &chunk_of_face, c);
        let p = patch_measures(&patch, &g.edges, area_floor);
        rows.push(ChunkRow {
            gaussian_patch: p.gaussian,
            gaussian_in_context: in_context_gaussian(&g, &chunk_of_face, c, &out.positions),
            mean_isolated: p.mean_isolated,
            mean_in_context: p.mean_in_context,
            seam_edges: p.seam_edges,
            interior_edges: p.interior_edges,
            boundary_vertices: p.boundary_vertices,
            interior_vertices: p.interior_vertices,
            faces: p.faces,
        });
    }

    let a = aggregate(&rows);
    let share = sharing(&g.faces, &chunk_of_face, g.vertices);
    let (borel_g, borel_h) = borel_sums(&g, &out.positions, lo, inv_width);

    let chi = report.euler_characteristic;
    let gaussian_global = TAU * chi as f64;
    let gaussian_gap = (a.gaussian_chunk_sum - gaussian_global).abs();
    let mean_gap = (a.mean_chunk_sum - g.mean).abs();
    let bit_identical = a.gaussian_bit_identical && a.mean_bit_identical;

    println!(
        "{:>10}  chunks {:>3}  ΣΦG {:>14.6} vs 2πχ {:>12.6}  gap {gaussian_gap:.4e} \
         (= {:.4} π)  ΣΦH {:.9} vs {:.9}  gap {mean_gap:.4e}  isolated bits: G {} H {} \
         (seams {})",
        F::NAME,
        a.chunks,
        a.gaussian_chunk_sum,
        gaussian_global,
        gaussian_gap / PI,
        a.mean_chunk_sum,
        g.mean,
        a.gaussian_bit_identical,
        a.mean_bit_identical,
        a.seam_edges,
    );

    run.record(&[
        ("field", F::NAME.to_string()),
        ("samples_per_axis", SAMPLES.to_string()),
        ("chunks", a.chunks.to_string()),
        ("gaussian_global", format!("{gaussian_global:.12}")),
        (
            "gaussian_chunk_sum",
            format!("{:.12}", a.gaussian_chunk_sum),
        ),
        ("gaussian_gap", format!("{gaussian_gap:.6e}")),
        ("mean_global", format!("{:.12}", g.mean)),
        ("mean_chunk_sum", format!("{:.12}", a.mean_chunk_sum)),
        ("mean_gap", format!("{mean_gap:.6e}")),
        ("isolated_chunks_bit_identical", bit_identical.to_string()),
        // ── clause three, split by which measure ────────────────────────────
        (
            "isolated_gaussian_bit_identical",
            a.gaussian_bit_identical.to_string(),
        ),
        (
            "isolated_mean_bit_identical",
            a.mean_bit_identical.to_string(),
        ),
        (
            "isolated_gaussian_mismatched_chunks",
            a.gaussian_mismatches.to_string(),
        ),
        (
            "isolated_mean_mismatched_chunks",
            a.mean_mismatches.to_string(),
        ),
        ("mean_isolated_sum", format!("{:.12}", a.mean_isolated_sum)),
        (
            "mean_isolation_gap",
            format!("{:.6e}", (a.mean_isolated_sum - g.mean).abs()),
        ),
        ("seam_edges_total", a.seam_edges.to_string()),
        (
            "chunk_boundary_vertices_total",
            a.boundary_vertices.to_string(),
        ),
        (
            "chunk_interior_vertices_total",
            a.interior_vertices.to_string(),
        ),
        ("chunk_interior_edges_total", a.interior_edges.to_string()),
        // ── clause one, diagnosed ──────────────────────────────────────────
        ("gaussian_gap_over_pi", format!("{:.9}", gaussian_gap / PI)),
        ("excess_chunk_incidence", share.excess.to_string()),
        (
            "gaussian_gap_matches_excess",
            (((gaussian_gap / PI) - share.excess as f64).abs() < 1e-6).to_string(),
        ),
        ("vertices_in_1_chunk", share.k_hist[1].to_string()),
        ("vertices_in_2_chunks", share.k_hist[2].to_string()),
        ("vertices_in_3_chunks", share.k_hist[3].to_string()),
        ("vertices_in_4_or_more_chunks", share.k_hist[4].to_string()),
        // ── the measure-theoretic decomposition, for contrast ───────────────
        ("borel_gaussian_sum", format!("{borel_g:.12}")),
        (
            "borel_gaussian_gap",
            format!("{:.6e}", (borel_g - gaussian_global).abs()),
        ),
        ("borel_mean_sum", format!("{borel_h:.12}")),
        (
            "borel_mean_gap",
            format!("{:.6e}", (borel_h - g.mean).abs()),
        ),
        // ── global census ──────────────────────────────────────────────────
        ("chunk_boxes", BOXES.to_string()),
        ("cell_size", format!("{cell:.9}")),
        ("chi_from_report", chi.to_string()),
        ("genus", format!("{:?}", report.genus)),
        (
            "gaussian_defect_sum_global",
            format!("{:.12}", g.gaussian_defect_sum),
        ),
        ("triangles", out.triangle_count().to_string()),
        ("vertices", out.vertex_count().to_string()),
        ("referenced_vertices", g.referenced_vertices.to_string()),
        ("global_edges", g.edges.len().to_string()),
        ("faces_skipped", g.skipped_faces.to_string()),
        ("degenerate_triangles", g.degenerate_faces.to_string()),
        ("global_boundary_edges", g.boundary_edges.to_string()),
        (
            "global_non_manifold_edges",
            g.non_manifold_edges.to_string(),
        ),
        (
            "report_non_manifold_vertices",
            report.non_manifold_vertices.to_string(),
        ),
        ("chunk_min_faces", a.min_faces.to_string()),
        ("chunk_max_faces", a.max_faces.to_string()),
        // ── the mean measure against the analytic integral ──────────────────
        ("mean_curvature_half", format!("{:.12}", g.mean / 2.0)),
        ("mean_smooth_integral_h", format!("{smooth_mean_h:.12}")),
        (
            "mean_half_relative_error",
            format!(
                "{:.6e}",
                (g.mean / 2.0 - smooth_mean_h).abs() / smooth_mean_h
            ),
        ),
        (
            "mean_half_exact",
            ((g.mean / 2.0 - smooth_mean_h).abs() < 1e-9).to_string(),
        ),
        // ── where the mean measure's turning actually lives ─────────────────
        ("mean_from_flat_edges", format!("{:.12}", g.mean_flat)),
        (
            "mean_from_right_angle_edges",
            format!("{:.12}", g.mean_right),
        ),
        ("mean_from_other_edges", format!("{:.12}", g.mean_other)),
        ("flat_edges", g.flat_edges.to_string()),
        ("right_angle_edges", g.right_edges.to_string()),
        ("other_edges", g.other_edges.to_string()),
        ("right_angle_edge_length", format!("{:.9}", g.right_length)),
    ]);
}

fn main() {
    let prereg = isomesh::experiment!("P-45");

    let sphere = Sphere::<f64>::canonical();
    let torus = Torus::<f64>::canonical();
    let cube = BoxExact::<f64>::canonical();

    // `∫_S H da` for each, from the canonical parameters rather than a literal.
    // Sphere of radius `r`: `H = 1/r` and area `4πr²`, so the integral is `4πr`.
    // Torus of radii `(R, a)`: `2π²R`, independent of the tube radius — the
    // `cos θ / (R + a cos θ)` principal curvature integrates to zero around the
    // tube, leaving `½ · 2π · 2πR`.
    // Box of half-extents `(a, b, c)`: twelve edges at exactly `π/2`, four of
    // each length, so `Σ l β = 4π(a + b + c)` and half of it is `2π(a + b + c)`
    // — `6π` for the canonical unit cube, and **exact** rather than convergent
    // if the extraction really is grid-aligned.
    let sphere_mean = 4.0 * PI * sphere.radius;
    let torus_mean = 2.0 * PI * PI * torus.major;
    let cube_mean = 2.0 * PI * (cube.half_extents[0] + cube.half_extents[1] + cube.half_extents[2]);

    common::experiment::run(prereg, |run| {
        sweep(run, &sphere, sphere_mean);
        sweep(run, &torus, torus_mean);
        sweep(run, &cube, cube_mean);
    });
}
