//! **P-42 — curvature as a normal-cycle measure, with the bound the source states.**
//!
//! Ticket: R-041. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p42
//! ```
//!
//! Writes `docs/experiments/p-42.csv`.
//!
//! # The formulas, transcribed rather than remembered
//!
//! From Sun & Morvan, *Curvature measures, normal cycles and asymptotic cones*,
//! `10.5802/acirm.50`, the in-corpus companion to Cohen-Steiner & Morvan's SoCG
//! paper, carrying the same statements:
//!
//! - **Theorem 5 (2)(b)** — the mean curvature measure of a closed polyhedron
//!   `P` is `Φ_P^H(B) = Σ_{e ∈ E ∩ B} l(e ∩ B) β(e)`, with `β(e) ∈ [−π, π]` the
//!   angle between the normals of the two triangles incident on `e`, **positive
//!   if `e` is convex and negative if concave**. That is the sign rule this file
//!   implements, and it is why the dihedral is computed with `atan2` against the
//!   edge direction rather than with `acos`, which cannot tell convex from
//!   concave.
//! - **Theorem 6** — for `P` *closely inscribed* in a smooth `S`,
//!   `|Φ_P^G(B) − Φ_S^G(pr(B))| ≤ C_S · K · ε` with
//!   `K = Σ_{t ⊂ B} cr(t)² + Σ_{t ⊂ B, t ∩ ∂B ≠ ∅} cr(t)` and
//!   `ε = max{cr(t), t ∈ T ∩ B}`, `cr` the circumradius.
//! - **Theorem 11** adds the hypothesis a marching-cubes mesh is most likely to
//!   break: the *fatness* `A(t) / l_max(t)²` of the triangulation must be
//!   uniformly bounded below by a positive constant. `min_triangle_fatness` is
//!   recorded for exactly that reason.
//!
//! ## `C_S` is not stated in a form any mesh can evaluate
//!
//! The paper says only that `C_S` is "a constant depending on the geometry of
//! `S`". No value, no formula, no dependence on any quantity a mesh carries. So
//! this harness does **not** invent one. The `bound` column is `K · ε` —
//! Theorem 6's bound with the unstated factor left out — and
//! `bound_constant_c` says so in the CSV. The falsifiable content that survives
//! is stronger than a single comparison anyway: `c_required = residual / (K · ε)`
//! is the smallest `C_S` for which clause one would hold, and Theorem 6 asserts
//! *one* constant for all resolutions. So a `c_required` that grows without
//! limit as `h` shrinks falsifies the theorem's form, and a `c_required` that
//! collapses towards zero says the bound is vacuous. Both are readable off the
//! CSV; neither needs a number the source declined to give.
//!
//! ## Where the source is internally inconsistent, and which reading is used
//!
//! Two places, both recorded rather than papered over.
//!
//! 1. **Theorem 5 (2)(a)** prints the Gaussian measure of a polyhedron as
//!    `Σ_{v ∈ V ∩ B} α_v`, where Definition 1 (1) has just defined `α_v` as the
//!    *sum* of incident corner angles. That total sits near `2π|V|`, not near
//!    `∫ G da`; the quantity the theory produces — and the one the
//!    pre-registration names — is the **defect** `2π − α_v`. This file uses the
//!    defect.
//! 2. **Definition 3 / Theorem 2 (3)** give the global mean curvature of a
//!    convex polyhedron as `Σ l(e) ∠(e)`, matched against `∫_S H da` with
//!    `H = ½ trace A` from §1.1. Those differ by a factor of two, and the cube
//!    settles it: for `[−1, 1]³` the Steiner `ε²` coefficient is twelve quarter
//!    cylinders, `12 · (π/4) · 2 = 6π`, while `Σ l(e) ∠(e) = 12 · 2 · (π/2) =
//!    12π`. So `∫_S H da = ½ Σ l(e) β(e)`. `mean_curvature_total` is the literal
//!    `Σ l(e) β(e)` of Theorem 5 (2)(b), because that is what the
//!    pre-registration names ("edge length times signed dihedral angle");
//!    `mean_curvature_half` and `mean_smooth_integral_h` sit beside it so the
//!    factor is visible in the artefact instead of chosen in silence.
//!
//! # The boundary term, and why it is a column rather than a branch
//!
//! On an open mesh the plain defect sum is not `2πχ`: the boundary vertices
//! carry a geodesic-curvature term this file does not compute. `sphere` and
//! `torus` are closed, which is why they are the registered fields, and
//! `boundary_vertices` / `boundary_edges` are recorded so a later reader
//! pointing this at an open field cannot mistake a missing term for a result.
//! Both are zero here; if they were not, `chi_from_defect` and `chi_from_report`
//! would separate and the assertion in [`sweep`] would fire.
//!
//! # Nothing here is timed
//!
//! Every column is a deterministic function of the extracted mesh, so there is
//! no median-of-five to report and no warm-up to do. The accumulation orders are
//! fixed on purpose — vertices in index order, edges in sorted-key order — so
//! the float sums reproduce run to run.

mod common;

use common::experiment::Run;
use isomesh::fields::{ReferenceField, Sphere, Torus};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{ValidateConfig, validate};
use isomesh::{MeshBuffer, Sdf};

/// Samples per axis. Each halves `h` against the one before it, which is what
/// clause two is stated in terms of.
const RESOLUTIONS: [u32; 3] = [33, 65, 129];

/// `2π`, the full turn a flat vertex's incident angles use up.
const TAU: f64 = std::f64::consts::TAU;

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
/// is exactly zero rather than a NaN. Not a silent substitution: such triangles
/// are counted in `degenerate_triangles` and `zero_length_sides`, and their
/// circumradius is `0/0`, which propagates into `bound` where it cannot be
/// missed.
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
/// conditioned for the sliver angles marching cubes actually produces, where the
/// dot product of two nearly parallel unit vectors rounds to exactly one and
/// `acos` returns a hard zero.
fn corner_angle(apex: [f64; 3], u: [f64; 3], v: [f64; 3]) -> f64 {
    let a = sub(u, apex);
    let b = sub(v, apex);
    norm(cross(a, b)).atan2(dot(a, b))
}

/// One considered face's geometry, computed once.
struct Face {
    /// Vertex indices as the input gave them — the traversal order the dihedral
    /// sign depends on.
    verts: [u32; 3],
    /// Outward unit normal, from the counter-clockwise-seen-from-outside winding
    /// this crate guarantees.
    normal: [f64; 3],
    /// Circumradius `abc / 4A`. Theorem 6's `cr(t)`.
    circumradius: f64,
    /// `A(t) / l_max(t)²`. Theorem 11's fatness.
    fatness: f64,
    /// `true` when `A(t) <= area_epsilon_rel · cell_size²`, the crate's own
    /// definition of degenerate.
    degenerate: bool,
}

/// The two curvature measures, and the census that makes them checkable.
struct Measures {
    /// `Σ_v (2π − α_v)`, vertices in index order.
    gaussian: f64,
    /// `Σ_v |2π − α_v|`. Not a curvature measure — it is what the Gaussian
    /// measure degenerates into if the sign of the defect is dropped, and it is
    /// recorded so that clause three's "an accumulator that lost the sign would
    /// fail loudly" is a number rather than an assurance.
    gaussian_abs: f64,
    /// `Σ_e l(e) β(e)`, edges in sorted-key order.
    mean: f64,
    /// `Σ_e l(e) |β(e)|`. The same sign check, for the mean measure.
    mean_abs: f64,
    /// Largest `|β(e)|` seen. Bounded by `π` by construction.
    max_abs_dihedral: f64,
    referenced_vertices: u64,
    considered_faces: u64,
    skipped_faces: u64,
    edges: u64,
    boundary_edges: u64,
    non_manifold_edges: u64,
    boundary_vertices: u64,
    /// Considered faces with a side of exactly zero length. Non-zero here is
    /// what would cost `corner_angle` the per-face `Σ α = π` identity, so it is
    /// the first thing to read if the two `chi` columns separate.
    zero_length_sides: u64,
}

/// Theorem 6's `K` and `ε`, and the fatness statistics that say whether the
/// theorem's hypotheses hold at all.
struct Bound {
    /// `Σ cr(t)²` over considered faces.
    k_area: f64,
    /// `Σ cr(t)` over faces incident to a boundary edge. Zero on a closed mesh.
    k_boundary: f64,
    /// `max cr(t)`.
    eps: f64,
    mean_circumradius: f64,
    min_fatness: f64,
    degenerate_faces: u64,
    /// Faces whose `abc / 4A` is not finite: a genuinely zero-area triangle.
    non_finite_circumradius: u64,
    /// The same three quantities over the faces the crate does *not* call
    /// degenerate. A separately named measure, not a fallback: both are always
    /// computed and both are always written.
    k_area_fat: f64,
    k_boundary_fat: f64,
    eps_fat: f64,
}

/// An occurrence of an undirected edge in one face.
struct EdgeRef {
    lo: u32,
    hi: u32,
    face: u32,
    /// `true` when this face traverses the edge `lo → hi`.
    forward: bool,
}

/// Faces with three distinct in-range indices, with normals and circumradii.
///
/// The predicate is deliberately the same one [`isomesh::validate`] uses for
/// `faces_skipped`, so this accumulator and `MeshReport` describe the same face
/// set and their two Euler characteristics are comparable.
fn faces_of(mesh: &MeshBuffer<f64>, cfg: &ValidateConfig) -> (Vec<Face>, u64, u64) {
    let nv = mesh.positions.len();
    let area_floor = cfg.area_epsilon_rel() * cfg.cell_size() * cfg.cell_size();
    let mut faces = Vec::with_capacity(mesh.indices.len() / 3);
    let mut skipped = 0u64;
    let mut zero_sides = 0u64;

    for tri in mesh.indices.chunks_exact(3) {
        let (a, b, c) = (tri[0], tri[1], tri[2]);
        let in_range = (a as usize) < nv && (b as usize) < nv && (c as usize) < nv;
        if !in_range || a == b || b == c || c == a {
            skipped += 1;
            continue;
        }
        let (pa, pb, pc) = (
            mesh.positions[a as usize],
            mesh.positions[b as usize],
            mesh.positions[c as usize],
        );
        let twice = cross(sub(pb, pa), sub(pc, pa));
        let area = norm(twice) / 2.0;
        let (la, lb, lc) = (norm(sub(pb, pc)), norm(sub(pc, pa)), norm(sub(pa, pb)));
        let l_max = la.max(lb).max(lc);
        if la.min(lb).min(lc) <= 0.0 {
            zero_sides += 1;
        }
        faces.push(Face {
            verts: [a, b, c],
            normal: unit(twice),
            circumradius: la * lb * lc / (4.0 * area),
            fatness: area / (l_max * l_max),
            degenerate: area <= area_floor,
        });
    }
    (faces, skipped, zero_sides)
}

/// Every `EdgeRef` of every considered face, sorted so that grouping is a linear
/// scan and the mean-curvature sum accumulates in a fixed order.
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

/// The signed dihedral of Theorem 5 (2)(b): `+` convex, `−` concave.
///
/// `atan2((n₁ × n₂) · ê, n₁ · n₂)` with `ê` the edge direction **as the first
/// face traverses it**. Swapping which face is first flips `ê` and `n₁ × n₂`
/// together, so the value does not depend on the order the group was sorted
/// into. Checked against the cube: top face `n₁ = +z` meets side face `n₂ = +y`
/// along an edge the top traverses in `−x`, giving `atan2(1, 0) = +π/2` —
/// convex, positive, as the sign rule requires. Both arguments carry the same
/// positive factor `|n₁||n₂|`, so `atan2` does not care that the normals are
/// unit; they are normalised anyway so that a zero-area face contributes
/// `atan2(0, 0) = 0` instead of a NaN.
fn signed_dihedral(first: &Face, second: &Face, edge_dir: [f64; 3]) -> f64 {
    let s = dot(cross(first.normal, second.normal), edge_dir);
    let c = dot(first.normal, second.normal);
    s.atan2(c)
}

/// Both measures, the census and Theorem 6's `K` and `ε`, from one sort of the
/// edge list.
///
/// The two are computed together because both need the same grouping: the mean
/// measure needs the interior edges and the bound's second term needs the
/// boundary ones, and sorting the 115,368 edge references of the 129³ sphere
/// twice to keep them in separate functions would be a pointless copy.
fn analyse(
    mesh: &MeshBuffer<f64>,
    faces: &[Face],
    skipped: u64,
    zero_sides: u64,
) -> (Measures, Bound) {
    let nv = mesh.positions.len();
    let mut angle_sum = vec![0.0_f64; nv];
    let mut referenced = vec![false; nv];

    for face in faces {
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

    let refs = sorted_edge_refs(faces);
    let mut on_boundary = vec![false; nv];
    let mut boundary_face = vec![false; faces.len()];
    let mut mean = 0.0_f64;
    let mut mean_abs = 0.0_f64;
    let mut max_abs_dihedral = 0.0_f64;
    let mut edges = 0u64;
    let mut boundary_edges = 0u64;
    let mut non_manifold_edges = 0u64;

    let mut i = 0usize;
    while i < refs.len() {
        let mut j = i + 1;
        while j < refs.len() && refs[j].lo == refs[i].lo && refs[j].hi == refs[i].hi {
            j += 1;
        }
        let group = &refs[i..j];
        edges += 1;
        let (lo, hi) = (group[0].lo as usize, group[0].hi as usize);
        match group.len() {
            1 => {
                boundary_edges += 1;
                on_boundary[lo] = true;
                on_boundary[hi] = true;
                boundary_face[group[0].face as usize] = true;
            }
            2 => {
                let first = &faces[group[0].face as usize];
                let second = &faces[group[1].face as usize];
                let along = sub(mesh.positions[hi], mesh.positions[lo]);
                let dir = if group[0].forward {
                    unit(along)
                } else {
                    unit(sub(mesh.positions[lo], mesh.positions[hi]))
                };
                let beta = signed_dihedral(first, second, dir);
                let contribution = norm(along) * beta;
                mean += contribution;
                mean_abs += contribution.abs();
                max_abs_dihedral = max_abs_dihedral.max(beta.abs());
            }
            _ => non_manifold_edges += 1,
        }
        i = j;
    }

    let mut gaussian = 0.0_f64;
    let mut gaussian_abs = 0.0_f64;
    let mut referenced_vertices = 0u64;
    for (&sum, &seen) in angle_sum.iter().zip(&referenced) {
        if seen {
            gaussian += TAU - sum;
            gaussian_abs += (TAU - sum).abs();
            referenced_vertices += 1;
        }
    }

    let measures = Measures {
        gaussian,
        gaussian_abs,
        mean,
        mean_abs,
        max_abs_dihedral,
        referenced_vertices,
        considered_faces: faces.len() as u64,
        skipped_faces: skipped,
        edges,
        boundary_edges,
        non_manifold_edges,
        boundary_vertices: on_boundary
            .iter()
            .zip(&referenced)
            .filter(|&(&b, &r)| b && r)
            .count() as u64,
        zero_length_sides: zero_sides,
    };
    (measures, bound_of(faces, &boundary_face))
}

/// Theorem 6's `K` and `ε`, plus the same over the non-degenerate subset.
///
/// `B` is the whole surface here, so `∂B` is the *mesh* boundary and
/// `boundary_face` is the incidence [`analyse`] already worked out.
fn bound_of(faces: &[Face], boundary_face: &[bool]) -> Bound {
    let mut out = Bound {
        k_area: 0.0,
        k_boundary: 0.0,
        eps: 0.0,
        mean_circumradius: 0.0,
        min_fatness: f64::INFINITY,
        degenerate_faces: 0,
        non_finite_circumradius: 0,
        k_area_fat: 0.0,
        k_boundary_fat: 0.0,
        eps_fat: 0.0,
    };
    let mut cr_sum = 0.0_f64;

    for (fi, face) in faces.iter().enumerate() {
        let cr = face.circumradius;
        cr_sum += cr;
        out.k_area += cr * cr;
        if cr > out.eps {
            out.eps = cr;
        }
        if boundary_face[fi] {
            out.k_boundary += cr;
        }
        if !cr.is_finite() {
            out.non_finite_circumradius += 1;
        }
        if face.degenerate {
            out.degenerate_faces += 1;
        } else {
            out.k_area_fat += cr * cr;
            if cr > out.eps_fat {
                out.eps_fat = cr;
            }
            if boundary_face[fi] {
                out.k_boundary_fat += cr;
            }
        }
        out.min_fatness = out.min_fatness.min(face.fatness);
    }
    out.mean_circumradius = cr_sum / faces.len() as f64;
    out
}

/// One field, three resolutions, three rows.
fn sweep<F>(run: &mut Run, field: &F, smooth_mean_h: f64)
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let chi_expected = field
        .expected_euler()
        .expect("sphere and torus both have an analytically known chi");
    let gaussian_expected = TAU * chi_expected as f64;

    let mut mc = MarchingCubes::new();
    let mut out = MeshBuffer::new();
    let mut previous_residual: Option<f64> = None;

    for n in RESOLUTIONS {
        let (shape, origin, cell) = common::grid(field, n);
        out.reset();
        mc.extract(field, &shape, origin, cell, &mut out)
            .expect("marching cubes on a closed reference field");

        let cfg = ValidateConfig::from_cell_size(cell).expect("cell size is positive");
        let report = validate(&out, &cfg);
        let (faces, skipped, zero_sides) = faces_of(&out, &cfg);
        let (m, b) = analyse(&out, &faces, skipped, zero_sides);

        let chi_report = report.euler_characteristic;
        let chi_from_defect = m.gaussian / TAU;
        assert!(
            (chi_from_defect - chi_report as f64).abs() < 1e-6,
            "{} at {n}³: chi from the defect sum is {chi_from_defect}, MeshReport \
             says {chi_report}. Either the accumulator is wrong, or the mesh has \
             boundary ({} boundary edges) or non-manifold edges ({}) — in which \
             case the plain defect sum is not 2πχ and the geodesic-curvature \
             boundary term this file does not compute is the missing piece.",
            F::NAME,
            m.boundary_edges,
            m.non_manifold_edges,
        );

        let residual = (m.gaussian - gaussian_expected).abs();
        let k = b.k_area + b.k_boundary;
        let bound = k * b.eps;
        let bound_fat = (b.k_area_fat + b.k_boundary_fat) * b.eps_fat;
        let chi_total = TAU * chi_report as f64;
        let ratio = previous_residual.map_or(f64::NAN, |prev| residual / prev);

        println!(
            "{:>6} {n:>4}³  Σdefect {:+.12}  residual {residual:.4e}  K·ε {bound:.4e}  \
             C_req {:.4e}  ratio {ratio:.4}  Σlβ {:.6}  ε {:.4e}  min fatness {:.4e}",
            F::NAME,
            m.gaussian,
            residual / bound,
            m.mean,
            b.eps,
            b.min_fatness,
        );

        run.record(&[
            ("field", F::NAME.to_string()),
            ("samples_per_axis", n.to_string()),
            ("gaussian_total", format!("{:.12}", m.gaussian)),
            ("gaussian_expected", format!("{gaussian_expected:.12}")),
            ("residual", format!("{residual:.6e}")),
            ("bound", format!("{bound:.6e}")),
            ("within_bound", (residual <= bound).to_string()),
            ("mean_curvature_total", format!("{:.9}", m.mean)),
            ("chi_from_defect", format!("{chi_from_defect:.9}")),
            ("chi_from_report", chi_report.to_string()),
            // ── the bound, decomposed and labelled ──────────────────────────
            ("bound_constant_c", "unstated_in_source".to_string()),
            ("bound_k", format!("{k:.9e}")),
            ("bound_k_area_term", format!("{:.9e}", b.k_area)),
            ("bound_k_boundary_term", format!("{:.9e}", b.k_boundary)),
            ("bound_eps", format!("{:.9e}", b.eps)),
            ("c_required", format!("{:.6e}", residual / bound)),
            ("bound_over_residual", format!("{:.6e}", bound / residual)),
            ("bound_excl_degenerate", format!("{bound_fat:.6e}")),
            (
                "within_bound_excl_degenerate",
                (residual <= bound_fat).to_string(),
            ),
            // ── where the residual comes from ───────────────────────────────
            (
                "residual_topological",
                format!("{:.6e}", (chi_total - gaussian_expected).abs()),
            ),
            (
                "residual_rounding",
                format!("{:.6e}", (m.gaussian - chi_total).abs()),
            ),
            ("residual_ratio_prev", format!("{ratio:.6}")),
            (
                "chi_from_defect_combinatorial",
                format!(
                    "{:.9}",
                    m.referenced_vertices as f64 - m.considered_faces as f64 / 2.0
                ),
            ),
            // ── mesh census ─────────────────────────────────────────────────
            ("cell_size", format!("{cell:.9}")),
            ("triangles", out.triangle_count().to_string()),
            ("vertices", out.vertex_count().to_string()),
            ("referenced_vertices", m.referenced_vertices.to_string()),
            ("edges", m.edges.to_string()),
            ("faces_considered", m.considered_faces.to_string()),
            ("faces_skipped", m.skipped_faces.to_string()),
            ("boundary_edges", m.boundary_edges.to_string()),
            ("boundary_vertices", m.boundary_vertices.to_string()),
            ("non_manifold_edges", m.non_manifold_edges.to_string()),
            ("zero_length_sides", m.zero_length_sides.to_string()),
            ("genus", format!("{:?}", report.genus)),
            // ── the shape of the triangles the bound is built from ──────────
            ("max_circumradius", format!("{:.9e}", b.eps)),
            ("mean_circumradius", format!("{:.9e}", b.mean_circumradius)),
            ("min_triangle_fatness", format!("{:.9e}", b.min_fatness)),
            ("degenerate_triangles", b.degenerate_faces.to_string()),
            (
                "non_finite_circumradius",
                b.non_finite_circumradius.to_string(),
            ),
            // ── the mean measure against the smooth integral ────────────────
            ("mean_curvature_half", format!("{:.9}", m.mean / 2.0)),
            ("mean_smooth_integral_h", format!("{smooth_mean_h:.9}")),
            (
                "mean_half_relative_error",
                format!(
                    "{:.6e}",
                    (m.mean / 2.0 - smooth_mean_h).abs() / smooth_mean_h
                ),
            ),
            ("max_abs_dihedral", format!("{:.9}", m.max_abs_dihedral)),
            // ── the sign check clause three exists to make ──────────────────
            ("gaussian_abs_defect_sum", format!("{:.9}", m.gaussian_abs)),
            ("mean_abs_measure", format!("{:.9}", m.mean_abs)),
        ]);

        previous_residual = Some(residual);
    }
}

fn main() {
    let prereg = isomesh::experiment!("P-42");

    let sphere = Sphere::<f64>::canonical();
    let torus = Torus::<f64>::canonical();

    // `∫_S H da` for each, from the canonical parameters rather than a literal.
    // Sphere of radius `r`: `H = 1/r` and area `4πr²`, so the integral is `4πr`.
    // Torus of radii `(R, a)`: `∫ H da = 2π²R`, independent of the tube radius —
    // the `cos θ / (R + a cos θ)` principal curvature integrates to zero around
    // the tube, leaving `½ · 2π · 2πR`.
    let sphere_mean = 4.0 * std::f64::consts::PI * sphere.radius;
    let torus_mean = 2.0 * std::f64::consts::PI * std::f64::consts::PI * torus.major;

    common::experiment::run(prereg, |run| {
        sweep(run, &sphere, sphere_mean);
        sweep(run, &torus, torus_mean);
    });
}
