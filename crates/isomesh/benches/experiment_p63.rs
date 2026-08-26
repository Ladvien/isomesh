//! **P-63 — the vertex-link case for Marching Cubes, by exhaustion at 2¹⁸.**
//!
//! Ticket: R-061. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p63
//! ```
//!
//! Writes `docs/experiments/p-63.csv`.
//!
//! # Why 18 corners is the whole space
//!
//! Every Marching Cubes vertex in this crate sits on a grid **edge**, and every
//! face incident to an edge vertex belongs to a cell containing that edge. Four
//! cells share a grid edge. Lay the block out with corners indexed
//! `x ∈ {0,1,2}`, `y ∈ {0,1,2}`, `z ∈ {0,1}` and the shared edge is the central
//! `z`-edge from `(1,1,0)` to `(1,1,1)`; the four cells are
//! `[x₀..x₀+1] × [y₀..y₀+1] × [0..1]` for `x₀, y₀ ∈ {0,1}`, and each of them
//! contains that edge. **3 × 3 × 2 = 18 corners, 2¹⁸ = 262,144 sign patterns**,
//! and the harness asserts all four cells contain the shared edge rather than
//! trusting the arithmetic above.
//!
//! # Which vertices this can honestly judge
//!
//! The block has a boundary, so most vertices have cells missing and a truncated
//! link. Three classes, by how many coordinates land exactly on the sample
//! lattice — and `✗43`'s own correction is why the count alone is not the test:
//!
//! - **2 or 3 integral** — on a grid edge. If it is `x = 1, y = 1` it is the
//!   **shared edge** and its link is complete: this is C1's vertex. Any other
//!   edge vertex is `link_defective_truncated`, reported and not counted.
//! - **0 integral** — inside a cell. A Marching Cubes interior apex or a dual
//!   vertex. Its faces all come from its own cell for the primal case, so the
//!   link is complete and it is counted; for a dual extractor it is **not**, and
//!   C3 says so.
//! - **1 integral** — on a cell face. Reported with the truncated set.
//!
//! # The field is the trilinear interpolant of the pattern, not a step function
//!
//! `sample` returns `±1` at the lattice and the **trilinear interpolation**
//! between, and `gradient` returns that interpolant's analytic gradient. Both
//! matter: `FaceAmbiguity::AsymptoticDecider` and `InteriorAmbiguity::Trilinear`
//! are statements *about the trilinear*, so a field that is not one would be
//! measuring a different object than the rules were derived for.
//!
//! # C2 reproduces the pre-fix defect without touching `src/`
//!
//! `✗43`'s defect was `Contours::fan` naming one shared `INTERIOR` apex for
//! every ring of a cell, so two rings longer than three glued into a bowtie. The
//! fix names ring `r`'s apex `INTERIOR + r`. **Merging all of one cell's
//! interior apexes back into a single vertex is exactly the pre-fix topology** —
//! the same triangles with one shared apex — because a vertex identification is
//! the operation the fix reversed. So the control needs no second extractor and
//! no edit to the library.

#![allow(clippy::float_cmp)]

mod common;

use std::time::Instant;

use isomesh::dual_contouring::DualContouring;
use isomesh::extractor::Extractor;
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::marching_cubes::{FaceAmbiguity, InteriorAmbiguity, MarchingCubes};
use isomesh::surface_nets::SurfaceNets;
use isomesh::weld::{Welder, epsilon_for};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Corners per axis: 3, 3, 2.
const NX: usize = 3;
const NY: usize = 3;
const NZ: usize = 2;
/// 18.
const CORNERS: usize = NX * NY * NZ;
/// 262,144.
const PATTERNS: u32 = 1 << CORNERS;

/// `(x, y, z)` to the pattern's bit index.
#[inline]
const fn corner_index(x: usize, y: usize, z: usize) -> usize {
    x + NX * (y + NY * z)
}

/// The shared edge's two endpoints: the block's central `z`-edge.
const SHARED_LO: usize = corner_index(1, 1, 0);
const SHARED_HI: usize = corner_index(1, 1, 1);

/// The four cells, as their eight corner bit indices in `corner_offset` order
/// (`x` fastest, then `y`, then `z`) — the same order `cube::corner_offset` uses.
fn cells() -> [[usize; 8]; 4] {
    let mut out = [[0usize; 8]; 4];
    let mut n = 0;
    for y0 in 0..2 {
        for x0 in 0..2 {
            for k in 0..8u32 {
                let (dx, dy, dz) = (
                    (k & 1) as usize,
                    ((k >> 1) & 1) as usize,
                    ((k >> 2) & 1) as usize,
                );
                out[n][k as usize] = corner_index(x0 + dx, y0 + dy, dz);
            }
            n += 1;
        }
    }
    out
}

// ─── the field ──────────────────────────────────────────────────────────────

/// The trilinear interpolant of one sign pattern over the 3 × 3 × 2 block.
///
/// Origin `0`, cell size `1`, so sample `(i, j, k)` is at position `(i, j, k)`
/// exactly and the containing cell is the floor of the position.
struct Lattice {
    value: [f64; CORNERS],
}

/// How a sign pattern becomes corner *values*.
///
/// **This is the fixture decision, and the first version of this harness got it
/// wrong in a way `M-44` caught.** `Unit` gives every corner `±1`, which is the
/// most symmetric magnitude assignment there is — and the trilinear's saddles
/// are then symmetric too, so `has_inner_hexagon`'s strict `0 < x < 1` test
/// rejects and **the interior-ambiguity rule never fires at all**. The "interior
/// rule on" arm was therefore identical to "off", `✗43`'s pre-fix fan had no
/// apex to merge, and C2's control could not fire. Measured, not assumed:
/// `interior_vertices` is a column.
///
/// `Generic` draws a magnitude per corner from SplitMix64 seeded on
/// `(pattern, corner)`, so the sweep stays **exhaustive over signs** — which is
/// what C1 is about — while the magnitudes are in general position. Both arms
/// run, because the `Unit` result is itself a finding about the interior rule.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Magnitudes {
    Unit,
    Generic,
}

impl Magnitudes {
    fn label(self) -> &'static str {
        match self {
            Self::Unit => "unit",
            Self::Generic => "generic",
        }
    }
}

/// SplitMix64, so a magnitude is a pure function of `(pattern, corner)` and the
/// sweep is byte-identical on every machine and every run.
#[inline]
fn magnitude(pattern: u32, corner: usize) -> f64 {
    let mut z = u64::from(pattern)
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(corner as u64 + 1)
        .wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    // In [1/4, 5/4): bounded away from zero, so no corner sits on the surface
    // and the sign pattern the sweep enumerates is the sign pattern meshed.
    0.25 + (z >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
}

impl Lattice {
    fn new(pattern: u32, magnitudes: Magnitudes) -> Self {
        let mut value = [0.0f64; CORNERS];
        for (bit, slot) in value.iter_mut().enumerate() {
            let sign = if pattern >> bit & 1 == 1 { -1.0 } else { 1.0 };
            *slot = match magnitudes {
                Magnitudes::Unit => sign,
                Magnitudes::Generic => sign * magnitude(pattern, bit),
            };
        }
        Self { value }
    }

    #[inline]
    fn at(&self, x: usize, y: usize, z: usize) -> f64 {
        self.value[corner_index(x, y, z)]
    }

    /// The containing cell's base index and the fraction inside it, clamped to
    /// the block. Clamping is a boundary condition, not a fallback: the extractor
    /// only ever asks inside, and the decider's saddle solve can ask on a face.
    #[inline]
    fn locate(&self, p: [f64; 3]) -> ([usize; 3], [f64; 3]) {
        let limit = [NX - 2, NY - 2, NZ - 2];
        let mut base = [0usize; 3];
        let mut frac = [0.0f64; 3];
        for k in 0..3 {
            let f = p[k].floor();
            let i = if f < 0.0 {
                0
            } else {
                (f as usize).min(limit[k])
            };
            base[k] = i;
            frac[k] = p[k] - i as f64;
        }
        (base, frac)
    }
}

impl Sdf for Lattice {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let (b, f) = self.locate(p);
        let mut acc = 0.0;
        for k in 0..8u32 {
            let (dx, dy, dz) = (
                (k & 1) as usize,
                ((k >> 1) & 1) as usize,
                ((k >> 2) & 1) as usize,
            );
            let wx = if dx == 1 { f[0] } else { 1.0 - f[0] };
            let wy = if dy == 1 { f[1] } else { 1.0 - f[1] };
            let wz = if dz == 1 { f[2] } else { 1.0 - f[2] };
            acc += wx * wy * wz * self.at(b[0] + dx, b[1] + dy, b[2] + dz);
        }
        acc
    }

    /// The trilinear's own analytic gradient. Not a central difference: the
    /// stencil would reach outside the block and measure the clamp.
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let (b, f) = self.locate(p);
        let mut g = [0.0f64; 3];
        for k in 0..8u32 {
            let d = [
                (k & 1) as usize,
                ((k >> 1) & 1) as usize,
                ((k >> 2) & 1) as usize,
            ];
            let v = self.at(b[0] + d[0], b[1] + d[1], b[2] + d[2]);
            for axis in 0..3 {
                let mut w = if d[axis] == 1 { 1.0 } else { -1.0 };
                for other in 0..3 {
                    if other != axis {
                        w *= if d[other] == 1 {
                            f[other]
                        } else {
                            1.0 - f[other]
                        };
                    }
                }
                g[axis] += w * v;
            }
        }
        g
    }
}

// ─── the link walk ──────────────────────────────────────────────────────────

/// How many coordinates of `p` land exactly on the sample lattice.
#[inline]
fn integral_coordinates(p: [f64; 3]) -> usize {
    (0..3).filter(|&k| p[k] == p[k].round()).count()
}

/// Is `p` the crossing on the block's shared edge?
///
/// The shared edge runs `(1,1,0)`–`(1,1,1)`, so its crossing has `x = 1` and
/// `y = 1` exactly and `z` strictly inside. The coordinate count alone would not
/// decide this — `✗43`'s own correction is that a **face** vertex also has one
/// integral coordinate and an edge vertex two, so the axis values are checked
/// rather than inferred.
#[inline]
fn on_shared_edge(p: [f64; 3]) -> bool {
    p[0] == 1.0 && p[1] == 1.0 && p[2] > 0.0 && p[2] < 1.0
}

/// Connected components of one vertex's incident-face link.
///
/// Two incident faces are adjacent when they share an **edge through the
/// vertex** — that is, the vertex plus one other corner. Counting components of
/// that relation is exactly what an edge census cannot see: two cones glued at a
/// point give every edge exactly two faces and χ its consequence, which is
/// `M-301`'s signature.
fn link_components(vertex: u32, incident: &[[u32; 3]]) -> usize {
    if incident.is_empty() {
        return 0;
    }
    // The two other corners of each incident face.
    let others: Vec<[u32; 2]> = incident
        .iter()
        .map(|t| {
            let mut o = [u32::MAX; 2];
            let mut n = 0;
            for &v in t {
                if v != vertex && n < 2 {
                    o[n] = v;
                    n += 1;
                }
            }
            o
        })
        .collect();

    let n = incident.len();
    let mut seen = vec![false; n];
    let mut components = 0usize;
    let mut stack = Vec::new();
    for start in 0..n {
        if seen[start] {
            continue;
        }
        components += 1;
        seen[start] = true;
        stack.push(start);
        while let Some(i) = stack.pop() {
            for j in 0..n {
                if seen[j] {
                    continue;
                }
                let shares = others[i].iter().any(|a| others[j].contains(a));
                if shares {
                    seen[j] = true;
                    stack.push(j);
                }
            }
        }
    }
    components
}

/// What one welded 4-cell mesh says about its vertex links.
#[derive(Default)]
struct LinkCensus {
    /// Vertices with 0 integral coordinates: a Marching Cubes interior apex, or
    /// a dual vertex. The reachability counter for the interior rule.
    interior_vertices: usize,
    /// The most incident faces any single vertex carried.
    max_incident_faces: usize,
    /// The shared-edge vertex exists and its link has one component.
    shared_edge_present: bool,
    /// The shared-edge vertex's link has more than one component. **C1.**
    shared_edge_defective: bool,
    /// Interior vertices (0 integral coordinates) with a split link.
    interior_defective: usize,
    /// Everything else with a split link: truncated links, reported not counted.
    truncated_defective: usize,
    /// The largest component count seen on any vertex.
    worst_components: usize,
}

fn census(positions: &[[f64; 3]], indices: &[u32]) -> LinkCensus {
    let mut out = LinkCensus::default();
    let mut incident: Vec<Vec<[u32; 3]>> = vec![Vec::new(); positions.len()];
    for tri in indices.as_chunks::<3>().0 {
        if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
            continue;
        }
        for &v in tri {
            incident[v as usize].push([tri[0], tri[1], tri[2]]);
        }
    }
    for (v, faces) in incident.iter().enumerate() {
        if faces.is_empty() {
            continue;
        }
        let components = link_components(v as u32, faces);
        out.worst_components = out.worst_components.max(components);
        out.max_incident_faces = out.max_incident_faces.max(faces.len());
        let p = positions[v];
        let split = components > 1;
        if on_shared_edge(p) {
            out.shared_edge_present = true;
            if split {
                out.shared_edge_defective = true;
            }
        } else if integral_coordinates(p) == 0 {
            out.interior_vertices += 1;
            out.interior_defective += usize::from(split);
        } else {
            out.truncated_defective += usize::from(split);
        }
    }
    out
}

/// Merge every cell's interior vertices into one, reproducing `✗43`'s pre-fix
/// single-apex fan. Returns the number of merges performed.
fn merge_interior_apexes(mesh: &mut MeshBuffer<f64>) -> usize {
    // Cell of an interior vertex is the floor of its position.
    let mut first: std::collections::HashMap<[i64; 3], u32> = std::collections::HashMap::new();
    let mut remap: Vec<u32> = (0..mesh.positions.len() as u32).collect();
    let mut merges = 0usize;
    for (v, p) in mesh.positions.iter().enumerate() {
        if integral_coordinates(*p) != 0 {
            continue;
        }
        let cell: [i64; 3] = std::array::from_fn(|k| p[k].floor() as i64);
        match first.get(&cell) {
            Some(&keep) => {
                remap[v] = keep;
                merges += 1;
            }
            None => {
                first.insert(cell, v as u32);
            }
        }
    }
    if merges > 0 {
        for i in &mut mesh.indices {
            *i = remap[*i as usize];
        }
    }
    merges
}

// ─── the well-composedness census, for C3 ───────────────────────────────────

/// The inside corners of a cell sign byte.
fn inside_corners(byte: u32) -> ([u32; 8], usize) {
    let mut out = [0u32; 8];
    let mut n = 0;
    for corner in 0..8u32 {
        if byte >> corner & 1 == 1 {
            out[n] = corner;
            n += 1;
        }
    }
    (out, n)
}

/// Exactly two inside corners sharing only a cell vertex.
fn is_vertex_diagonal_pair(byte: u32) -> bool {
    let (corners, n) = inside_corners(byte);
    n == 2 && (corners[0] ^ corners[1]) == 0b111
}

/// Some 2 × 2 face of the cell is a checkerboard.
fn has_checkerboard_face(byte: u32) -> bool {
    for axis in 0..3u32 {
        let diagonal = 0b111 ^ (1 << axis);
        for side in 0..2u32 {
            let mut inside = [0u32; 4];
            let mut n = 0;
            for corner in 0..8u32 {
                if corner >> axis & 1 == side && byte >> corner & 1 == 1 {
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

/// `P-41`'s classification, transcribed rather than imported: 120 sign bytes are
/// 2D-critical, 8 are 3D-critical, 128 are one or the other.
fn critical_bytes() -> [bool; 256] {
    let mut out = [false; 256];
    for byte in 0..256u32 {
        out[byte as usize] = has_checkerboard_face(byte)
            || is_vertex_diagonal_pair(byte)
            || is_vertex_diagonal_pair(!byte & 0xFF);
    }
    out
}

/// The cell sign byte of `cell` under `pattern`, in `corner_offset` order.
#[inline]
fn cell_byte(pattern: u32, cell: &[usize; 8]) -> u32 {
    let mut byte = 0u32;
    for (k, &bit) in cell.iter().enumerate() {
        if pattern >> bit & 1 == 1 {
            byte |= 1 << k;
        }
    }
    byte
}

// ─── one arm ────────────────────────────────────────────────────────────────

struct Arm {
    name: &'static str,
    extractor: &'static str,
    interior_rule: &'static str,
    merge_apexes: bool,
    magnitudes: Magnitudes,
}

#[allow(clippy::too_many_lines)]
fn run_arm(
    arm: &Arm,
    cells: &[[usize; 8]; 4],
    critical: &[bool; 256],
) -> Vec<(&'static str, String)> {
    let shape = RuntimeShape3::new([NX as u32, NY as u32, NZ as u32]).expect("3x3x2 fits");
    let origin = [0.0f64; 3];
    let mut mesh = MeshBuffer::<f64>::new();
    let mut welder = Welder::<f64>::new();
    let epsilon = epsilon_for(1.0f64);

    let mut mc = MarchingCubes::<f64>::new();
    mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
    if arm.interior_rule == "trilinear" {
        mc.set_interior_ambiguity(InteriorAmbiguity::Trilinear);
    }
    let mut sn = SurfaceNets::<f64>::new();
    let mut dc = DualContouring::<f64>::new();
    let mut mdc = ManifoldDualContouring::<f64>::new();

    let started = Instant::now();
    let mut patterns_cut = 0u64;
    let mut shared_edge_vertices = 0u64;
    let mut defective_shared = 0u64;
    let mut defective_interior = 0u64;
    let mut defective_truncated = 0u64;
    let mut worst_components = 0usize;
    let mut first_defective = String::from("none");
    let mut fan_patterns = 0u64;
    let mut critical_cells = 0u64;
    let mut critical_patterns = 0u64;
    let mut interior_vertices = 0u64;
    let mut max_incident_faces = 0usize;

    for pattern in 0..PATTERNS {
        let field = Lattice::new(pattern, arm.magnitudes);

        let cut = (pattern >> SHARED_LO & 1) != (pattern >> SHARED_HI & 1);
        if cut {
            patterns_cut += 1;
        }

        let mut crit_here = 0u64;
        for cell in cells {
            if critical[cell_byte(pattern, cell) as usize] {
                crit_here += 1;
            }
        }
        critical_cells += crit_here;
        if crit_here > 0 {
            critical_patterns += 1;
        }

        mesh.reset();
        match arm.extractor {
            "marching_cubes" => mc.extract_into(&field, &shape, origin, 1.0, &mut mesh),
            "surface_nets" => sn.extract_into(&field, &shape, origin, 1.0, &mut mesh),
            "dual_contouring" => dc.extract_into(&field, &shape, origin, 1.0, &mut mesh),
            "manifold_dual_contouring" => mdc.extract_into(&field, &shape, origin, 1.0, &mut mesh),
            other => panic!("{other} is not one of the four the registration names"),
        }
        .expect("extraction on an 18-corner block");

        if arm.merge_apexes && merge_interior_apexes(&mut mesh) > 0 {
            fan_patterns += 1;
        }

        welder
            .weld(&mut mesh, epsilon)
            .expect("weld on an 18-corner block");

        let c = census(&mesh.positions, &mesh.indices);
        if c.shared_edge_present {
            shared_edge_vertices += 1;
        }
        worst_components = worst_components.max(c.worst_components);
        max_incident_faces = max_incident_faces.max(c.max_incident_faces);
        interior_vertices += c.interior_vertices as u64;
        defective_interior += c.interior_defective as u64;
        defective_truncated += c.truncated_defective as u64;
        if c.shared_edge_defective {
            defective_shared += 1;
            if first_defective == "none" {
                first_defective = format!("{pattern:#07x}");
            }
        }
        if first_defective == "none" && c.interior_defective > 0 {
            first_defective = format!("interior@{pattern:#07x}");
        }
    }

    let wall_ms = started.elapsed().as_millis();
    let is_dual = arm.extractor != "marching_cubes";
    let defective_total = defective_shared + defective_interior + defective_truncated;
    let c1 = !is_dual && !arm.merge_apexes && defective_shared == 0 && defective_interior == 0;
    let c2 = arm.merge_apexes && defective_total > 0;
    let c3 = is_dual && defective_total > 0 && defective_total == critical_cells;

    println!(
        "{:<30} {:<8} {:<9} {:>8} {:>9} {:>8} {:>8} {:>8} {:>4} {:>7}",
        arm.name,
        arm.magnitudes.label(),
        arm.interior_rule,
        shared_edge_vertices,
        interior_vertices,
        defective_shared,
        defective_interior,
        defective_truncated,
        worst_components,
        wall_ms
    );

    vec![
        ("arm", arm.name.to_string()),
        ("extractor", arm.extractor.to_string()),
        ("interior_rule", arm.interior_rule.to_string()),
        ("magnitudes", arm.magnitudes.label().to_string()),
        ("interior_vertices", interior_vertices.to_string()),
        ("max_incident_faces", max_incident_faces.to_string()),
        ("patterns", PATTERNS.to_string()),
        ("patterns_shared_edge_cut", patterns_cut.to_string()),
        ("shared_edge_vertices", shared_edge_vertices.to_string()),
        ("link_defective_shared_edge", defective_shared.to_string()),
        ("link_defective_interior", defective_interior.to_string()),
        ("link_defective_truncated", defective_truncated.to_string()),
        ("worst_link_components", worst_components.to_string()),
        ("first_defective_pattern", first_defective),
        ("fan_patterns", fan_patterns.to_string()),
        ("critical_cells", critical_cells.to_string()),
        ("critical_patterns", critical_patterns.to_string()),
        (
            "defective_equals_critical",
            (defective_total == critical_cells).to_string(),
        ),
        ("link_defective_total", defective_total.to_string()),
        ("c1_holds", c1.to_string()),
        ("c2_holds", c2.to_string()),
        ("c3_holds", c3.to_string()),
        ("wall_ms", wall_ms.to_string()),
    ]
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-63");
    common::experiment::run(prereg, |run| {
        let cells = cells();
        // The block's geometry, asserted rather than trusted.
        assert_eq!(CORNERS, 18, "the block is 3 x 3 x 2 corners");
        assert_eq!(PATTERNS, 262_144, "2^18 sign patterns");
        for cell in &cells {
            assert!(
                cell.contains(&SHARED_LO) && cell.contains(&SHARED_HI),
                "every one of the four cells must contain the shared edge, or the \
                 block is not the four cells around it"
            );
            let mut sorted = *cell;
            sorted.sort_unstable();
            for pair in sorted.windows(2) {
                assert_ne!(pair[0], pair[1], "a cell names a corner twice");
            }
        }
        let critical = critical_bytes();
        let n_critical = critical.iter().filter(|c| **c).count();
        assert_eq!(
            n_critical, 128,
            "P-41's classification puts 128 of 256 sign bytes critical"
        );
        println!(
            "block: {CORNERS} corners, {PATTERNS} patterns, 4 cells around the \
             central z-edge ({SHARED_LO} -> {SHARED_HI}); {n_critical} of 256 \
             cell sign bytes critical\n"
        );
        println!(
            "{:<30} {:<8} {:<9} {:>8} {:>9} {:>8} {:>8} {:>8} {:>4} {:>7}",
            "arm",
            "mags",
            "interior",
            "sharedV",
            "interiorV",
            "defShare",
            "defInt",
            "defTrunc",
            "cmp",
            "ms"
        );

        let mc_arm =
            |name: &'static str, interior_rule: &'static str, merge_apexes: bool, magnitudes| Arm {
                name,
                extractor: "marching_cubes",
                interior_rule,
                merge_apexes,
                magnitudes,
            };
        let dual_arm = |name: &'static str, extractor: &'static str, magnitudes| Arm {
            name,
            extractor,
            interior_rule: "n/a",
            merge_apexes: false,
            magnitudes,
        };
        let arms = [
            // Signs exhaustive at the most symmetric magnitude there is. The
            // pure combinatorial statement, and the arm that shows whether the
            // interior rule can fire at all here.
            mc_arm("mc/off/unit", "ignore", false, Magnitudes::Unit),
            mc_arm("mc/on/unit", "trilinear", false, Magnitudes::Unit),
            // Signs exhaustive with magnitudes in general position, which is
            // what reaches the interior rule and therefore C2.
            mc_arm("mc/off/generic", "ignore", false, Magnitudes::Generic),
            mc_arm("mc/on/generic", "trilinear", false, Magnitudes::Generic),
            mc_arm(
                "mc/pre_fix_apex/generic",
                "trilinear",
                true,
                Magnitudes::Generic,
            ),
            dual_arm("surface_nets/generic", "surface_nets", Magnitudes::Generic),
            dual_arm(
                "dual_contouring/generic",
                "dual_contouring",
                Magnitudes::Generic,
            ),
            dual_arm(
                "manifold_dual_contouring/generic",
                "manifold_dual_contouring",
                Magnitudes::Generic,
            ),
        ];

        let mut rows = Vec::new();
        for arm in &arms {
            rows.push((arm.name, run_arm(arm, &cells, &critical)));
        }

        let value = |name: &str, key: &str| -> String {
            rows.iter()
                .find(|(n, _)| *n == name)
                .and_then(|(_, r)| r.iter().find(|(k, _)| *k == key))
                .map(|(_, v)| v.clone())
                .unwrap_or_default()
        };
        let num = |name: &str, key: &str| -> u64 { value(name, key).parse().unwrap_or(0) };

        // **C1's population, asserted.** Half the patterns cut the shared edge,
        // and the shared-edge vertex must exist on every one of them or the
        // sweep is not walking the link it claims to.
        let cut = num("mc/on/generic", "patterns_shared_edge_cut");
        assert_eq!(
            cut,
            u64::from(PATTERNS) / 2,
            "the shared edge should be cut on exactly half the patterns"
        );
        assert_eq!(
            num("mc/on/generic", "shared_edge_vertices"),
            cut,
            "every pattern that cuts the shared edge must put a vertex on it"
        );

        // **The reachability control, and it is why this harness has two
        // magnitude arms.** The interior rule has to be shown firing before
        // "with the interior rule on" means anything, and C2 has to be shown
        // able to return non-zero before C1's zero means anything (M-44).
        let interior_unit = num("mc/on/unit", "interior_vertices");
        let interior_generic = num("mc/on/generic", "interior_vertices");
        let fan = num("mc/pre_fix_apex/generic", "fan_patterns");
        let fan_defects = num("mc/pre_fix_apex/generic", "link_defective_total");
        println!("\ninterior vertices: unit {interior_unit}, generic {interior_generic}");
        println!(
            "C1 shared-edge link defects, off/on, unit:    {} / {}",
            num("mc/off/unit", "link_defective_shared_edge"),
            num("mc/on/unit", "link_defective_shared_edge")
        );
        println!(
            "C1 shared-edge link defects, off/on, generic: {} / {}",
            num("mc/off/generic", "link_defective_shared_edge"),
            num("mc/on/generic", "link_defective_shared_edge")
        );
        println!(
            "   interior-apex link defects, off/on, generic: {} / {}",
            num("mc/off/generic", "link_defective_interior"),
            num("mc/on/generic", "link_defective_interior")
        );
        println!("C2 pre-fix arm: {fan} patterns fanned, {fan_defects} link defects");
        println!(
            "C3 critical cells {} vs dual link defects sn {} / dc {} / mdc {}; max \
             incident faces on any dual vertex {}",
            num("surface_nets/generic", "critical_cells"),
            num("surface_nets/generic", "link_defective_total"),
            num("dual_contouring/generic", "link_defective_total"),
            num("manifold_dual_contouring/generic", "link_defective_total"),
            num("dual_contouring/generic", "max_incident_faces")
        );
        assert!(
            interior_generic > 0,
            "VOID: the interior rule produced no apex on any of {PATTERNS} \
             generic-magnitude patterns, so 'with the interior rule on' is the \
             same arm as 'off' and neither C1 nor C2 means anything"
        );
        assert!(
            fan > 0,
            "VOID: the pre-fix arm merged no apexes on any of {PATTERNS} \
             patterns, so C2 cannot fire and C1's zero proves nothing"
        );

        for (_, row) in rows {
            run.record(&row);
        }
    });
}
