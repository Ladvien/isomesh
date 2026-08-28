//! **P-73 - the angle-weighted pseudonormal against the field's gradient.**
//!
//! Ticket: R-073. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p73
//! ```
//!
//! Writes `docs/experiments/p-73.csv`.
//!
//! # What is being decided
//!
//! The mechanics dossier proposes an **angle-weighted pseudonormal** as a ~40
//! line addition to `isomesh::normals::NormalStrategy`. Jin, Lewis & West,
//! *A comparison of algorithms for vertex normal computation*, The Visual
//! Computer 21(1-2), 2005 (`10.1007/s00371-004-0271-1`) scored six connectivity
//! weightings against the analytic normal and reported a **median discrepancy of
//! 5-20 degrees for every one of them on marching-tetrahedra output**, at the
//! highest resolution they tested. This harness puts that to the crate's own
//! eight reference fields on trilinear Marching Cubes output. **Nothing in
//! `src/` is changed and no `NormalStrategy` variant is added: the six
//! weightings live here, in the bench, which is the whole point of a row that
//! exists to stop a change rather than make one.**
//!
//! The weightings are Jin, Lewis & West's own six, in the registration's order:
//! MWE (equal), MWA (incident angle), MWSELR (sine over edge-length product),
//! MWAAT (adjacent triangle area), MWELR (reciprocal edge-length product),
//! MWRELR (reciprocal square root of the edge-length product). Each sums
//! `w * unit_face_normal` over the incident triangles and normalises.
//!
//! # SHARE, recomputed before the harness was written
//!
//! The registration's line is *"C1 and C2 move the whole normals stage; C3 moves
//! nothing and is a correctness clause."* Recomputing it:
//!
//! - **C1's ratio is not an Amdahl ratio and has no share cap.** Numerator and
//!   denominator are two angular-error distributions over *the same* vertex
//!   population, each denominated in itself; `ratio_to_gradient` can take any
//!   value in `[0, inf)`, so `1/(1 - share/factor)` does not bind and the
//!   registered 3x bar is arithmetically reachable. What *can* make it
//!   unreachable is a **floor on the denominator**, and that floor is not small
//!   on one of the two fields the registration names as worst:
//!   `CentralDifference` at the cell size differences the field over a step `h`,
//!   and `thin_plate`'s **full thickness is 0.025** while `h` is 0.125 at 33^3
//!   and 0.0625 at 65^3 - **5.0x and 2.5x the plate's whole thickness**. The
//!   stencil steps clean through the plate and out the other side, so
//!   `gradient_median_deg` there is not a small reference error, it is a second
//!   broken instrument. **C1 is therefore predicted to fail on `thin_plate` for
//!   a reason that has nothing to do with connectivity weighting**, and that is
//!   stated here before the run rather than discovered in the numbers. The same
//!   arithmetic is benign elsewhere: `sphere` has curvature radius 1 against
//!   `h = 0.125`, `gyroid` a period of `2*pi` against `h = 0.4375`.
//!
//!   **Post-run correction, kept beside the prediction rather than replacing
//!   it.** The conclusion was right and the mechanism was wrong.
//!   `gradient_median_deg` on `thin_plate` is not large, it is **exactly zero**:
//!   a slab is axis-aligned, so the differencing error is confined to the one
//!   axis the normal already points along and the two tangential differences are
//!   exactly zero either way. The step really is 5x the plate's thickness and
//!   `cd_magnitude_ratio_median` measures that - but a normal is a direction, and
//!   the direction survives. C1 is unreachable on `thin_plate`, as predicted, by
//!   `0/0` rather than by a large denominator.
//! - **C2's rank correlation has no share cap either** - Spearman's rho is
//!   bounded to `[-1, 1]` and its standard error over the tens of thousands of
//!   vertices per case is under 0.005, so the registered 0.7 is far outside
//!   sampling noise. Its vacuity risk is at the other end: if `|f(v)|` were
//!   identically zero the ranks would all tie and rho would not exist. That is
//!   not hypothetical: on `box_exact` **every** Marching Cubes vertex lands
//!   exactly on the surface, `canary_mean_abs_f` is exactly zero and
//!   `canary_zero_vertices == scored_vertices`, so C2 is vacuous there and
//!   `case_c2_verdict` says so.
//! - **C3 moves nothing and is a correctness clause**, as registered.
//!
//! # M-289 is directly on this path, and here is the answer to it
//!
//! M-289: a reference gradient that normalises a **cancellation residue** is a
//! random unit vector at exactly the points being measured, and it falsified two
//! true hypotheses. Every number below is an angle against
//! `field.gradient(vertex)`, so that failure is available again. Three
//! independent checks, all on the row and all asserted:
//!
//! 1. **`reference_matches_extraction`** - the reference normals computed here
//!    are asserted **bit-identical** to the ones the extractor already stored,
//!    so this harness cannot be measuring a second, private idea of the
//!    gradient.
//! 2. **`ref_closed_form_median_deg` / `_max_deg`** - a closed form written from
//!    the geometry rather than from `fields/mod.rs`, for the five fields that
//!    have one. `sphere` -> `(p - c)/|p - c|`. `torus` -> the revolved-circle
//!    normal. `box_exact` and `thin_plate` -> **must be exactly one of the six
//!    axis directions**, which is a property a residue-normalising bug cannot
//!    satisfy and which is independent of how the crate computes it.
//!    `csg_difference` -> an axis direction *or* the cut sphere's inward normal,
//!    whichever is nearer. `-1` where there is no closed form.
//! 3. **`ref_fd_median_deg` / `_p99_deg` / `_over_half_deg`** - the angle to a
//!    **fourth-order** central difference of `sample` at `1e-6` of the domain
//!    extent, on all eight fields. Round-off there is `~1e-10` and truncation
//!    `~1e-24`, so it is an independent gradient to ten digits. It disagrees
//!    legitimately at a kink, which is why the **median** is what is asserted:
//!    M-289's bug hit about half the vertices, so it moves a median; a CSG seam
//!    is a percent or two and moves only the tail, which is reported and not
//!    asserted.
//!
//! # Vacuity control
//!
//! `worst_decile_triangles`, and the registration's wording is load-bearing:
//! *"the fixture must contain triangles in the bottom decile of the aspect-ratio
//! distribution on every field"*. Both readings of that sentence are emitted,
//! because they answer different questions and the first run showed the
//! difference matters.
//!
//! - `worst_decile_triangles` is the **registered column** under the per-field
//!   reading: the count of this case's triangles at or below **its own** 10th
//!   percentile of the radius ratio `16*A^2 / ((a+b+c) * a*b*c)` (1 for
//!   equilateral, 0 for degenerate). It is `n/10` for any non-empty mesh, so the
//!   assertion on it fires on exactly one thing - a case that meshed nothing -
//!   and that is stated rather than dressed up.
//! - `worst_decile_max_radius_ratio` is the number that makes the registered
//!   one mean something: the decile's upper edge. A field whose worst tenth is
//!   at 0.5 has no badly shaped triangles at all.
//! - `worst_decile_triangles_below_pooled_p10` is the **strict** reading - one
//!   aspect-ratio distribution, the fixture's, pooled over all sixteen cases,
//!   its 10th percentile in `radius_ratio_p10_pooled` - and it **can be zero**.
//! - `triangles_below_0p15` is the dossier's own radius-ratio bar for *"over 15
//!   degrees"*, and `min_radius_ratio` is the single worst triangle.
//!
//! A case with no triangle below 0.15 cannot exhibit the effect the proposal
//! claims, so it is scored **vacuous for C1** in `case_c1_verdict` rather than
//! panicked over: a field with uniformly good triangles does not invalidate the
//! harness, it invalidates the clause on that field, which is what the
//! registration's own instruction to *"score it VACUOUS and say why"* asks for.
//! What **is** asserted is that the experiment as a whole has an instrument:
//! `c1_testable_cases` must be non-zero.
//!
//! # C1's denominator can be exactly zero, and on three fields it is
//!
//! `ratio_to_gradient` is `NaN` wherever `gradient_median_deg` is `0.000000`,
//! and `case_c1_testable` says so. That is not an instrument fault: an
//! axis-aligned exact-box SDF is **piecewise linear** near its faces, Marching
//! Cubes' linear interpolation is exact on a linear field, and a central
//! difference of a linear function is its derivative - so on `box_exact`,
//! `thin_plate` and `csg_difference` the median error of *every* method is
//! exactly zero and the ratio is `0/0`. `zero_error_vertices` counts the
//! bit-exact agreements per arm, `cd_magnitude_ratio_median` shows that a
//! central difference can have the exactly right *direction* with a badly wrong
//! *magnitude*, and `p99_angle_error_deg` with `ratio_p99_to_gradient` is where
//! those three fields still carry a signal.
//!
//! # One population, seven arms
//!
//! M-281: every arm is scored over **the same vertex set**, in one build and one
//! run. A vertex is dropped from all seven arms, and counted in a column, if
//!
//! - it touches a mesh edge used by other than exactly two triangles - an
//!   incomplete fan is the connectivity route scored on a fan it does not have,
//!   and `fbm_terrain` exits through the domain wall by construction
//!   (`closed_in_domain() == false`);
//! - its reference gradient, its central difference, or any of the six weighted
//!   sums has zero or non-finite length.
//!
//! `scored_vertices` is on the row and asserted non-zero.
//!
//! # No clock appears in this experiment
//!
//! M-280 asks for cycles or ratios rather than nanoseconds on a governed CPU.
//! Nothing here is timed: every column is an angle, an exact-arithmetic count, a
//! rank correlation or a 64-bit hash. The connectivity route's *cost* is not in
//! question - the paper's finding and this harness's are both about accuracy -
//! so no timing column exists to be misread.
//!
//! # C3 is a cross-machine clause, answered by argument *and* measurement
//!
//! The registered premise is *"angle weighting needs `acos`, libm's `acos` has no
//! architecture selection, and a connectivity-weighted normal is therefore a
//! golden-hash liability"*. **The crate's own `src/real.rs` already contradicts
//! the conclusion**: `Real::acos` exists, was added at S-006 *"the angle-weighted
//! pseudonormal weights each face normal by the incident angle, and `libm` is the
//! only source of that angle that stays bit-identical across platforms - which is
//! what T-007's golden hashes depend on"*, and is wired to `libm::acos`. Reading
//! libm 0.2.16: `acos` is generic pure-Rust f64 arithmetic whose only helper is
//! `sqrt`, and `src/math/arch/mod.rs` selects `sqrt` per architecture - `fsqrt`
//! on aarch64+neon, `sqrtsd` on x86-64+sse2 - both of which IEEE-754 requires to
//! be correctly rounded. So the prediction here is **zero** hash movement on the
//! connectivity route, which is C3's registered falsifier.
//!
//! That argument is not a measurement, so the measurement is made. A single-file
//! program with no dependencies (`p73_route.rs.in`) is written into
//! `target/p73_crossmachine`, compiled with a bare `rustc`, and run here and on
//! the LAN's Apple M5 (`ssh mac_air`). Per fixture it hashes:
//!
//! - the six connectivity normal arrays, using libm 0.2.16's `acos` transcribed
//!   verbatim into the program;
//! - the angle-weighted array again using **`std`'s** `acos` - the platform's own
//!   libm, which is the arm that *could* move and is what the registration's
//!   premise would require the crate to have used;
//! - the gradient-route normal array, with the fixture's gradient transcribed
//!   from `fields/mod.rs` so the whole route runs on both machines. All three
//!   fixtures were chosen because their gradients are transcendental-free -
//!   `sqrt`, divide, `signum`, `max` - so nothing in the arm predicted not to
//!   move depends on `sin`, `cos` or a Perlin lattice;
//! - `acos` itself over a 1,000,001-point sweep of `[-1, 1]`, both ways, which
//!   isolates the transcendental from the mesh entirely.
//!
//! Three controls make the program admissible as evidence.
//! `port_acos_matches_libm` asserts the transcription's sweep hash equals
//! `libm::acos`'s computed in-process here. `program_matches_bench` asserts all
//! twenty-three of the program's hashes equal the ones this bench computed
//! through the real crate. `hashes_moved_codegen` re-runs the local program at
//! `opt-level=0`, because the two machines are not on the same rustc and a
//! codegen-sensitive hash could not be attributed to an architecture.
//!
//! **If the second machine cannot be reached or has no toolchain, C3 is
//! `blocked`** - `second_machine` reads `absent`, the reason is in
//! `second_machine_note`, and the moved counts carry the local-only value of
//! zero, which must **not** be read as C3's falsifier. Guessing is the one
//! unacceptable outcome.
//!
//! # Assertion order
//!
//! Every control asserts **after** the CSV is written. A tripped control should
//! leave the evidence that tripped it on disk; `x35` and `x52` are both entries
//! whose numbers outlived the file that held them.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::needless_range_loop,
    clippy::print_literal,
    clippy::too_many_lines
)]

mod common;

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use isomesh::extractor::Extractor;
use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::normals::{NormalStrategy, central_difference, recompute};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Samples per axis, as registered.
const RESOLUTIONS: [u32; 2] = [33, 65];

/// The six connectivity weightings, in the registration's order.
const WEIGHTINGS: [&str; 6] = [
    "equal",
    "angle",
    "sine_edge_length_recip",
    "adjacent_triangle_area",
    "edge_length_recip",
    "sqrt_edge_length_recip",
];

/// Index of the angle weighting. The dossier's own proposal, and the arm C2 is
/// scored on.
const ANGLE: usize = 1;

/// Index of the area weighting. It is the crate's existing
/// `NormalStrategy::AreaWeightedFaces`, which is what makes it a control.
const AREA: usize = 3;

/// The seventh arm's name in the `weighting` column.
const GRADIENT_ROW: &str = "central_difference";

/// C1's registered bar.
const C1_RATIO_BAR: f64 = 3.0;

/// C2's registered bar, and the count of fields that must clear it.
const C2_RHO_BAR: f64 = 0.7;
const C2_FIELDS_NEEDED: usize = 6;

/// The dossier's radius-ratio bar for *"over 15 degrees"*.
const RADIUS_RATIO_BAR: f64 = 0.15;

/// Fourth-order probe step, as a fraction of the domain extent.
const FD_PROBE_FRACTION: f64 = 1e-6;

/// How far the reference gradient's median may sit from an independent one
/// before the reference is not trusted. M-289's bug moved about half the
/// vertices, so it moved a median.
const REF_MEDIAN_BAR_DEG: f64 = 0.05;

/// Fixtures the cross-machine hash comparison runs on.
const HASH_FIXTURES: [&str; 3] = ["sphere", "torus", "box_exact"];

/// Resolution the hash fixtures are meshed at. 33 is in T-007's own golden set.
const HASH_SAMPLES: u32 = 33;

/// Points in the `acos` sweep, so a one-ulp disagreement anywhere in `[-1, 1]`
/// is unlikely to be missed. Two million calls; nothing here is timed.
const ACOS_SWEEP: u64 = 1_000_000;

/// The second machine. An `~/.ssh/config` alias; Apple M5, aarch64, macOS.
const REMOTE_HOST: &str = "mac_air";

// ── FNV-1a, the crate's golden-hash function ────────────────────────────────
//
// Mirrors `isomesh::validate::mesh_hash`'s hasher exactly - published constants,
// bit patterns rather than values - because the question C3 asks is about the
// crate's committed hashes. Only the normals are hashed: positions and indices
// are shipped to the second machine, so they are inputs rather than outputs.

struct Fnv(u64);

impl Fnv {
    const OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;

    const fn new() -> Self {
        Self(Self::OFFSET_BASIS)
    }

    fn write_u64(&mut self, v: u64) {
        for b in v.to_le_bytes() {
            self.0 ^= u64::from(b);
            self.0 = self.0.wrapping_mul(Self::PRIME);
        }
    }

    fn write_f64(&mut self, v: f64) {
        self.write_u64(v.to_bits());
    }
}

fn hash_normals(normals: &[[f64; 3]]) -> u64 {
    let mut h = Fnv::new();
    h.write_u64(normals.len() as u64);
    for n in normals {
        for c in n {
            h.write_f64(*c);
        }
    }
    h.0
}

fn hash_acos(f: fn(f64) -> f64) -> u64 {
    let mut h = Fnv::new();
    for i in 0..=ACOS_SWEEP {
        let x = -1.0 + 2.0 * (i as f64) / (ACOS_SWEEP as f64);
        h.write_f64(f(x));
    }
    h.0
}

// ── vector helpers, mirroring `src/vec3.rs` operation for operation ─────────

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
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

fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// `src/normals.rs::from_field`'s normalise, step for step: length, then scale
/// by its reciprocal. `None` exactly where that function returns
/// `Error::DegenerateNormal`.
fn unit(g: [f64; 3]) -> Option<[f64; 3]> {
    let len = length(g);
    if len > 0.0 && len.is_finite() {
        Some(scale(g, len.recip()))
    } else {
        None
    }
}

/// A degenerate sum hashes as `[0, 0, 0]` rather than being dropped, so the
/// array length is the vertex count on both machines and a hash cannot move
/// because one side skipped a vertex.
fn unit_or_zero(sums: &[[f64; 3]]) -> Vec<[f64; 3]> {
    sums.iter().map(|s| unit(*s).unwrap_or([0.0; 3])).collect()
}

fn angle_deg(a: [f64; 3], b: [f64; 3]) -> f64 {
    dot(a, b).clamp(-1.0, 1.0).acos().to_degrees()
}

// ── order statistics ────────────────────────────────────────────────────────

fn sorted(v: &[f64]) -> Vec<f64> {
    let mut s = v.to_vec();
    s.sort_by(f64::total_cmp);
    s
}

fn median_of_sorted(s: &[f64]) -> f64 {
    if s.is_empty() {
        return 0.0;
    }
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        0.5 * (s[n / 2 - 1] + s[n / 2])
    }
}

/// Nearest-rank quantile. No interpolation: a tail figure that invents a value
/// between two measurements is a value nothing measured.
fn quantile_of_sorted(s: &[f64], q: f64) -> f64 {
    if s.is_empty() {
        return 0.0;
    }
    let i = ((s.len() - 1) as f64 * q).round() as usize;
    s[i.min(s.len() - 1)]
}

/// Mid-ranks, so ties are handled the way Spearman's rho requires.
fn ranks(v: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..v.len()).collect();
    order.sort_by(|&a, &b| v[a].total_cmp(&v[b]));
    let mut r = vec![0.0; v.len()];
    let mut i = 0;
    while i < order.len() {
        let mut j = i;
        while j + 1 < order.len() && v[order[j + 1]] == v[order[i]] {
            j += 1;
        }
        let mid = (i + j) as f64 / 2.0 + 1.0;
        for k in i..=j {
            r[order[k]] = mid;
        }
        i = j + 1;
    }
    r
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    if x.len() < 2 {
        return 0.0;
    }
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    if sxx <= 0.0 || syy <= 0.0 {
        // Zero variance on one side: every value tied. Not a correlation of
        // zero - a correlation that does not exist, which `canary_zero_vertices`
        // and `canary_max_abs_f` are on the row to expose.
        return 0.0;
    }
    sxy / (sxx * syy).sqrt()
}

fn spearman(x: &[f64], y: &[f64]) -> f64 {
    pearson(&ranks(x), &ranks(y))
}

/// Radius ratio `2r/R`, normalised to 1 for equilateral and 0 for degenerate:
/// `16*A^2 / ((a+b+c) * a*b*c)`.
fn radius_ratio(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let la = length(sub(b, a));
    let lb = length(sub(c, b));
    let lc = length(sub(a, c));
    let area = 0.5 * length(cross(sub(b, a), sub(c, a)));
    let den = (la + lb + lc) * la * lb * lc;
    if den > 0.0 && area.is_finite() {
        (16.0 * area * area / den).min(1.0)
    } else {
        0.0
    }
}

// ── the six weightings ──────────────────────────────────────────────────────

/// Accumulate every weighting's unnormalised sum in one pass over the triangles.
///
/// `acos` is a parameter so the same code serves `libm::acos` and `std`'s.
///
/// The area arm accumulates the **raw cross product** rather than
/// `area * unit_normal`, which is what `src/normals.rs::area_weighted` does and
/// is why `mwaat_crate_status` can be a bit-exact assertion: the cross product's
/// magnitude is twice the area, so summing raw cross products *is* the
/// area-weighted sum, and dividing each by its own length first would throw the
/// weighting away and add a rounding.
fn weighted_sums(
    positions: &[[f64; 3]],
    indices: &[u32],
    acos: fn(f64) -> f64,
) -> ([Vec<[f64; 3]>; 6], usize) {
    let mut acc: [Vec<[f64; 3]>; 6] = std::array::from_fn(|_| vec![[0.0; 3]; positions.len()]);
    let mut skipped_corners = 0usize;

    for tri in indices.as_chunks::<3>().0 {
        let (i, j, k) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
        if i >= positions.len() || j >= positions.len() || k >= positions.len() {
            continue;
        }
        let (a, b, c) = (positions[i], positions[j], positions[k]);
        let cr = cross(sub(b, a), sub(c, a));

        // MWAAT, exactly `area_weighted`'s loop and in its order.
        for vertex in [i, j, k] {
            for axis in 0..3 {
                acc[AREA][vertex][axis] += cr[axis];
            }
        }

        let twice_area = length(cr);
        // `<= 0.0` is false for NaN, which `!is_finite()` then catches: the two
        // together are exactly "not a usable positive length".
        if twice_area <= 0.0 || !twice_area.is_finite() {
            // A zero-area triangle has no unit face normal. It contributes to
            // MWAAT - nothing, correctly - and to nothing else.
            skipped_corners += 3;
            continue;
        }
        let face = scale(cr, twice_area.recip());

        for &(p0, p1, p2, vertex) in &[(a, b, c, i), (b, c, a, j), (c, a, b, k)] {
            let e1 = sub(p1, p0);
            let e2 = sub(p2, p0);
            let l1 = length(e1);
            let l2 = length(e2);
            let lprod = l1 * l2;
            if lprod <= 0.0 || !lprod.is_finite() {
                skipped_corners += 1;
                continue;
            }
            // sin(alpha) = |e1 x e2| / (l1 l2) = twice_area / lprod, so MWSELR
            // needs no second transcendental. Only MWA does, which is C3's hook.
            let weights = [
                1.0,
                acos((dot(e1, e2) / lprod).clamp(-1.0, 1.0)),
                (twice_area / lprod) / lprod,
                0.0, // MWAAT, already accumulated above
                lprod.recip(),
                lprod.sqrt().recip(),
            ];
            for w in [0, 1, 2, 4, 5] {
                for axis in 0..3 {
                    acc[w][vertex][axis] += weights[w] * face[axis];
                }
            }
        }
    }
    (acc, skipped_corners)
}

// ── the reference gradient, checked three ways ──────────────────────────────

/// Fourth-order central difference: `(-f(+2s) + 8f(+s) - 8f(-s) + f(-2s)) / 12s`.
fn high_order_gradient<S: Sdf<Scalar = f64>>(sdf: &S, p: [f64; 3], s: f64) -> [f64; 3] {
    let mut g = [0.0; 3];
    for axis in 0..3 {
        let mut q = p;
        q[axis] = p[axis] + 2.0 * s;
        let f2 = sdf.sample(q);
        q[axis] = p[axis] + s;
        let f1 = sdf.sample(q);
        q[axis] = p[axis] - s;
        let m1 = sdf.sample(q);
        q[axis] = p[axis] - 2.0 * s;
        let m2 = sdf.sample(q);
        g[axis] = (-f2 + 8.0 * f1 - 8.0 * m1 + m2) / (12.0 * s);
    }
    g
}

/// Angle from `n` to the nearest of the six axis directions.
///
/// The independent property for the two box-derived fields: an exact signed
/// distance to an axis-aligned box has an axis-aligned gradient wherever its
/// nearest feature is a face, and a residue-normalising bug returns a direction
/// with no relation to the geometry.
fn axis_deviation_deg(n: [f64; 3]) -> f64 {
    let mut best = 180.0f64;
    for axis in 0..3 {
        for sign in [1.0, -1.0] {
            let mut e = [0.0; 3];
            e[axis] = sign;
            best = best.min(angle_deg(n, e));
        }
    }
    best
}

/// A closed form written from the geometry, not from `fields/mod.rs`.
///
/// `None` for `gyroid`, `fbm_terrain` and `noise_cavity`, whose level sets have
/// no closed-form normal; those three rest on the fourth-order check alone.
fn closed_form_deviation_deg(name: &str, p: [f64; 3], n: [f64; 3]) -> Option<f64> {
    match name {
        // Unit sphere at the origin.
        "sphere" => unit(p).map(|r| angle_deg(n, r)),
        // Ring in the xz-plane, major 1, minor 0.3: the normal points away from
        // the nearest point of the core circle.
        "torus" => {
            let s = (p[0] * p[0] + p[2] * p[2]).sqrt();
            if s == 0.0 {
                return None;
            }
            let radial = s - 1.0;
            unit([radial * (p[0] / s), p[1], radial * (p[2] / s)]).map(|r| angle_deg(n, r))
        }
        "box_exact" | "thin_plate" => Some(axis_deviation_deg(n)),
        // `[-1,1]^3` minus a sphere of radius 0.75 at (0.6, 0.6, 0.6): either a
        // box face's outward normal or the cut sphere's *inward* one.
        "csg_difference" => {
            let d = sub(p, [0.6; 3]);
            let sphere = unit(d).map_or(180.0, |r| angle_deg(n, scale(r, -1.0)));
            Some(axis_deviation_deg(n).min(sphere))
        }
        _ => None,
    }
}

// ── one (field, resolution) case ────────────────────────────────────────────

struct Case {
    field: &'static str,
    samples: u32,
    cell_size: f64,
    vertices: usize,
    triangles: usize,
    scored: usize,
    boundary_excluded: usize,
    nonmanifold_edges: usize,
    ref_degenerate: usize,
    cd_degenerate: usize,
    weight_degenerate: usize,
    skipped_corners: usize,
    /// Per-vertex angular error over the scored population: six weightings then
    /// the central difference.
    errors: [Vec<f64>; 7],
    canary: Vec<f64>,
    canary_zero: usize,
    radius_ratios: Vec<f64>,
    below_bar: usize,
    /// Per-arm count of vertices whose angular error is **exactly** zero. The
    /// sharpest statement of "this case cannot falsify anything".
    zero_error: [usize; 7],
    /// Median of `|central difference| / |analytic gradient|`. Exactly `1` when
    /// the difference is exact; `sin(h)/h` on a trigonometric field, where the
    /// scalar cancels under normalisation; far from `1` on a slab whose
    /// thickness is under the step.
    cd_mag_ratio_median: f64,
    /// This case's own 10th percentile of the radius ratio: the upper edge of
    /// the decile `worst_decile_triangles` counts.
    field_p10: f64,
    /// The single worst triangle in this case.
    min_radius_ratio: f64,
    ref_fd_median: f64,
    ref_fd_p99: f64,
    ref_fd_over_half: usize,
    ref_closed_median: f64,
    ref_closed_max: f64,
    reference_mismatch: usize,
    mwaat_mismatch: usize,
    mwaat_status: &'static str,
}

fn measure<F>(name: &'static str, field: &F, samples: u32) -> Case
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("the reference grid fits u32");
    let mut mesh = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract_into(field, &shape, lo, cell_size, &mut mesh)
        .expect("marching cubes extraction of a reference field");

    let vertices = mesh.vertex_count();
    let triangles = mesh.indices.len() / 3;

    // ── mesh-boundary and non-manifold vertices ──────────────────────────────
    let mut edge_use: HashMap<(u32, u32), u32> = HashMap::new();
    for tri in mesh.indices.as_chunks::<3>().0 {
        for pair in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            let key = if pair.0 <= pair.1 {
                pair
            } else {
                (pair.1, pair.0)
            };
            *edge_use.entry(key).or_insert(0) += 1;
        }
    }
    let mut incomplete = vec![false; vertices];
    let mut nonmanifold_edges = 0usize;
    for (&(u, v), &count) in &edge_use {
        if count != 2 {
            if count > 2 {
                nonmanifold_edges += 1;
            }
            incomplete[u as usize] = true;
            incomplete[v as usize] = true;
        }
    }

    // ── the seven candidate normals per vertex ───────────────────────────────
    let (sums, skipped_corners) = weighted_sums(&mesh.positions, &mesh.indices, libm::acos);

    let mut reference: Vec<Option<[f64; 3]>> = Vec::with_capacity(vertices);
    let mut cd: Vec<Option<[f64; 3]>> = Vec::with_capacity(vertices);
    // Magnitudes, kept because `cd_magnitude_ratio_median` is what explains a
    // central difference that has the exactly right direction and the wrong
    // length - `sin(h)/h` on a trigonometric field, `thickness/h` on a slab.
    let mut mag_ratio: Vec<f64> = Vec::with_capacity(vertices);
    for &p in &mesh.positions {
        let g = field.gradient(p);
        let d = central_difference(field, p, cell_size);
        let gl = length(g);
        mag_ratio.push(if gl > 0.0 { length(d) / gl } else { 0.0 });
        reference.push(unit(g));
        cd.push(unit(d));
    }

    // Control 1: the reference is the extractor's own normal, bit for bit.
    let mut reference_mismatch = 0usize;
    for (v, r) in reference.iter().enumerate() {
        if let Some(r) = r {
            let stored = mesh.normals[v];
            if r[0].to_bits() != stored[0].to_bits()
                || r[1].to_bits() != stored[1].to_bits()
                || r[2].to_bits() != stored[2].to_bits()
            {
                reference_mismatch += 1;
            }
        }
    }

    // Control 2: the area arm is the crate's `AreaWeightedFaces`, bit for bit.
    let mut crate_area = mesh.clone();
    let (mwaat_status, mwaat_mismatch) =
        if recompute(&mut crate_area, field, NormalStrategy::AreaWeightedFaces).is_ok() {
            let mut bad = 0usize;
            for v in 0..vertices {
                match unit(sums[AREA][v]) {
                    Some(n) => {
                        let c = crate_area.normals[v];
                        if n[0].to_bits() != c[0].to_bits()
                            || n[1].to_bits() != c[1].to_bits()
                            || n[2].to_bits() != c[2].to_bits()
                        {
                            bad += 1;
                        }
                    }
                    None => bad += 1,
                }
            }
            (if bad == 0 { "match" } else { "mismatch" }, bad)
        } else {
            ("recompute_error", 0)
        };

    // ── the scored population ────────────────────────────────────────────────
    let mut boundary_excluded = 0usize;
    let mut ref_degenerate = 0usize;
    let mut cd_degenerate = 0usize;
    let mut weight_degenerate = 0usize;
    let mut scored_ids: Vec<usize> = Vec::new();
    for v in 0..vertices {
        if incomplete[v] {
            boundary_excluded += 1;
        } else if reference[v].is_none() {
            ref_degenerate += 1;
        } else if cd[v].is_none() {
            cd_degenerate += 1;
        } else if (0..6).any(|w| unit(sums[w][v]).is_none()) {
            weight_degenerate += 1;
        } else {
            scored_ids.push(v);
        }
    }

    // ── the seven arms, the canary, and the reference's own checks ───────────
    let probe = FD_PROBE_FRACTION * (hi[0] - lo[0]);
    let mut errors: [Vec<f64>; 7] = std::array::from_fn(|_| Vec::with_capacity(scored_ids.len()));
    let mut canary = Vec::with_capacity(scored_ids.len());
    let mut canary_zero = 0usize;
    let mut zero_error = [0usize; 7];
    let mut ref_fd = Vec::with_capacity(scored_ids.len());
    let mut ref_closed = Vec::new();
    let mut mag_scored = Vec::with_capacity(scored_ids.len());
    for &v in &scored_ids {
        let p = mesh.positions[v];
        let r = reference[v].expect("a scored vertex has a reference normal");
        for w in 0..6 {
            let n = unit(sums[w][v]).expect("a scored vertex has every weighted normal");
            let e = angle_deg(n, r);
            if e == 0.0 {
                zero_error[w] += 1;
            }
            errors[w].push(e);
        }
        let n = cd[v].expect("a scored vertex has a central difference");
        let e = angle_deg(n, r);
        if e == 0.0 {
            zero_error[6] += 1;
        }
        errors[6].push(e);
        mag_scored.push(mag_ratio[v]);

        let f = field.sample(p).abs();
        if f == 0.0 {
            canary_zero += 1;
        }
        canary.push(f);

        match unit(high_order_gradient(field, p, probe)) {
            Some(fd) => ref_fd.push(angle_deg(r, fd)),
            None => ref_fd.push(180.0),
        }
        if let Some(dev) = closed_form_deviation_deg(name, p, r) {
            ref_closed.push(dev);
        }
    }

    let ref_fd_sorted = sorted(&ref_fd);
    let ref_closed_sorted = sorted(&ref_closed);
    let mag_sorted = sorted(&mag_scored);

    // ── triangle shape ───────────────────────────────────────────────────────
    let mut radius_ratios = Vec::with_capacity(triangles);
    let mut below_bar = 0usize;
    for tri in mesh.indices.as_chunks::<3>().0 {
        let (i, j, k) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
        if i >= vertices || j >= vertices || k >= vertices {
            continue;
        }
        let q = radius_ratio(mesh.positions[i], mesh.positions[j], mesh.positions[k]);
        if q < RADIUS_RATIO_BAR {
            below_bar += 1;
        }
        radius_ratios.push(q);
    }
    let ratio_sorted = sorted(&radius_ratios);

    Case {
        field: name,
        samples,
        cell_size,
        vertices,
        triangles,
        scored: scored_ids.len(),
        boundary_excluded,
        nonmanifold_edges,
        ref_degenerate,
        cd_degenerate,
        weight_degenerate,
        skipped_corners,
        errors,
        canary,
        canary_zero,
        radius_ratios,
        below_bar,
        zero_error,
        cd_mag_ratio_median: median_of_sorted(&mag_sorted),
        field_p10: quantile_of_sorted(&ratio_sorted, 0.10),
        min_radius_ratio: ratio_sorted.first().copied().unwrap_or(0.0),
        ref_fd_median: median_of_sorted(&ref_fd_sorted),
        ref_fd_p99: quantile_of_sorted(&ref_fd_sorted, 0.99),
        ref_fd_over_half: ref_fd.iter().filter(|d| **d > 0.5).count(),
        ref_closed_median: if ref_closed_sorted.is_empty() {
            -1.0
        } else {
            median_of_sorted(&ref_closed_sorted)
        },
        ref_closed_max: ref_closed_sorted.last().copied().unwrap_or(-1.0),
        reference_mismatch,
        mwaat_mismatch,
        mwaat_status,
    }
}

// ── C3: the cross-machine hash comparison ───────────────────────────────────

/// One machine's whole hash set: six connectivity arrays per fixture, the
/// std-`acos` angle array per fixture, the gradient array per fixture, and
/// `acos` itself over the sweep, twice.
#[derive(Clone, PartialEq, Eq)]
struct Hashes {
    connectivity: Vec<u64>,
    std_acos: Vec<u64>,
    gradient: Vec<u64>,
    sweep_libm: u64,
    sweep_std: u64,
}

fn count_moved(a: &[u64], b: &[u64]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

struct Crossmachine {
    local: Hashes,
    reached_second: bool,
    second_machine: String,
    second_toolchain: String,
    note: String,
    local_toolchain: String,
    moved_connectivity: usize,
    moved_gradient: usize,
    moved_std_acos: usize,
    moved_std_acos_fixtures: String,
    moved_sweep_libm: usize,
    moved_sweep_std: usize,
    moved_codegen: usize,
    codegen_control_ran: bool,
    program_matches_bench: bool,
    port_matches_libm: bool,
    remote_sweep_libm: u64,
    remote_sweep_std: u64,
}

fn short(program: &str, args: &[&str]) -> String {
    Command::new(program)
        .args(args)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_default()
}

/// `rustc 1.98.0 (88d9e12ae 2026-08-18)` -> `1.98.0`.
fn version_only(v: &str) -> String {
    v.split_whitespace().nth(1).unwrap_or("unknown").to_string()
}

fn slug(v: &str) -> String {
    let s: String = v
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() {
                c.to_ascii_lowercase()
            } else {
                '-'
            }
        })
        .collect();
    let trimmed = s.trim_matches('-').to_string();
    if trimmed.is_empty() {
        String::from("unknown")
    } else {
        trimmed
    }
}

/// Parse the program's `label hex` lines. `fixture <name>` does not parse as
/// hex and is skipped, which is what makes it usable as a block separator.
fn parse_hashes(out: &str) -> HashMap<String, u64> {
    let mut map = HashMap::new();
    for line in out.lines() {
        let mut parts = line.split_whitespace();
        let parsed = parts
            .next()
            .zip(parts.next())
            .and_then(|(label, value)| u64::from_str_radix(value, 16).ok().map(|v| (label, v)));
        if let Some((label, v)) = parsed {
            map.insert(label.to_string(), v);
        }
    }
    map
}

fn hashes_from_output(outputs: &[String]) -> Option<Hashes> {
    let mut connectivity = Vec::new();
    let mut std_acos = Vec::new();
    let mut gradient = Vec::new();
    let mut sweep_libm = None;
    let mut sweep_std = None;
    for out in outputs {
        let map = parse_hashes(out);
        for w in WEIGHTINGS {
            connectivity.push(map.get(&format!("conn_{w}")).copied()?);
        }
        std_acos.push(map.get("conn_angle_std_acos").copied()?);
        gradient.push(map.get("grad").copied()?);
        let l = map.get("acos_sweep_libm").copied()?;
        let s = map.get("acos_sweep_std").copied()?;
        // The sweep is the same computation in every run; a machine that
        // disagreed with itself would make everything else unattributable.
        if sweep_libm.is_some_and(|p| p != l) || sweep_std.is_some_and(|p| p != s) {
            return None;
        }
        sweep_libm = Some(l);
        sweep_std = Some(s);
    }
    Some(Hashes {
        connectivity,
        std_acos,
        gradient,
        sweep_libm: sweep_libm?,
        sweep_std: sweep_std?,
    })
}

/// The remote runs all three fixtures in one shell, so the concatenated output
/// is split back on the `fixture` line each run prints first.
fn split_fixture_outputs(text: &str) -> Option<Hashes> {
    let mut blocks: Vec<String> = Vec::new();
    let mut current = String::new();
    for line in text.lines() {
        if line.starts_with("fixture ") && !current.is_empty() {
            blocks.push(std::mem::take(&mut current));
        }
        current.push_str(line);
        current.push('\n');
    }
    if !current.is_empty() {
        blocks.push(current);
    }
    if blocks.len() != HASH_FIXTURES.len() {
        return None;
    }
    hashes_from_output(&blocks)
}

/// The generated single-file program.
const ROUTE_PROGRAM: &str = include_str!("p73_route.rs.in");

/// Mesh one hash fixture, write it as exact bit patterns, and return the six
/// connectivity hashes, the std-`acos` angle hash and the gradient hash as this
/// bench computes them through the real crate.
fn write_fixture_mesh(dir: &Path, name: &str) -> Option<(Vec<u64>, u64, u64)> {
    let mut text = String::new();
    let mut conn = Vec::new();
    let std_conn: u64;
    let grad: u64;

    macro_rules! build {
        ($field:expr) => {{
            let field = $field;
            let (lo, hi) = field.domain();
            let cell_size = (hi[0] - lo[0]) / f64::from(HASH_SAMPLES - 1);
            let shape = RuntimeShape3::new([HASH_SAMPLES; 3]).expect("the hash grid fits u32");
            let mut mesh = MeshBuffer::<f64>::new();
            MarchingCubes::<f64>::new()
                .extract_into(&field, &shape, lo, cell_size, &mut mesh)
                .expect("marching cubes extraction of a hash fixture");
            text.push_str(&format!(
                "{} {} {} {}\n",
                name,
                HASH_SAMPLES,
                mesh.positions.len(),
                mesh.indices.len() / 3
            ));
            for p in &mesh.positions {
                text.push_str(&format!(
                    "{:016x} {:016x} {:016x}\n",
                    p[0].to_bits(),
                    p[1].to_bits(),
                    p[2].to_bits()
                ));
            }
            for tri in mesh.indices.as_chunks::<3>().0 {
                text.push_str(&format!("{} {} {}\n", tri[0], tri[1], tri[2]));
            }

            let (libm_sums, _) = weighted_sums(&mesh.positions, &mesh.indices, libm::acos);
            let (std_sums, _) = weighted_sums(&mesh.positions, &mesh.indices, f64::acos);
            for w in 0..6 {
                conn.push(hash_normals(&unit_or_zero(&libm_sums[w])));
            }
            std_conn = hash_normals(&unit_or_zero(&std_sums[ANGLE]));
            let g: Vec<[f64; 3]> = mesh
                .positions
                .iter()
                .map(|p| unit(field.gradient(*p)).unwrap_or([0.0; 3]))
                .collect();
            grad = hash_normals(&g);
        }};
    }

    match name {
        "sphere" => build!(isomesh::fields::Sphere::<f64>::canonical()),
        "torus" => build!(isomesh::fields::Torus::<f64>::canonical()),
        "box_exact" => build!(isomesh::fields::BoxExact::<f64>::canonical()),
        _ => return None,
    }
    std::fs::write(dir.join(format!("mesh_{name}.txt")), text).ok()?;
    Some((conn, std_conn, grad))
}

/// One ssh round trip. The script is wrapped in `bash -lc '...'` as a **single**
/// argv element, because ssh space-joins its remaining arguments and hands the
/// result to the remote login shell to re-split. No script here contains a
/// single quote.
fn ssh(script: &str) -> Option<String> {
    let wrapped = format!("bash -lc '{script}'");
    let out = Command::new("ssh")
        .args([
            "-o",
            "BatchMode=yes",
            "-o",
            "ConnectTimeout=10",
            REMOTE_HOST,
            &wrapped,
        ])
        .output()
        .ok()?;
    if out.status.success() {
        String::from_utf8(out.stdout).ok()
    } else {
        None
    }
}

fn send_payload(payload: &[u8]) -> bool {
    use std::io::Write as _;
    let wrapped = String::from(
        "bash -lc 'rm -rf ~/p73_crossmachine && mkdir -p ~/p73_crossmachine && \
         tar -xz -C ~/p73_crossmachine'",
    );
    let child = Command::new("ssh")
        .args([
            "-o",
            "BatchMode=yes",
            "-o",
            "ConnectTimeout=10",
            REMOTE_HOST,
            &wrapped,
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn();
    match child {
        Ok(mut c) => {
            let written = c
                .stdin
                .as_mut()
                .is_some_and(|stdin| stdin.write_all(payload).is_ok());
            drop(c.stdin.take());
            written && c.wait().map(|s| s.success()).unwrap_or(false)
        }
        Err(_) => false,
    }
}

fn crossmachine(dir: &Path) -> Crossmachine {
    let local_toolchain = version_only(&short("rustc", &["--version"]));
    let mut note = String::new();

    let _ = std::fs::create_dir_all(dir);
    let program = dir.join("p73_route.rs");
    let mut bench_conn: Vec<u64> = Vec::new();
    let mut bench_std: Vec<u64> = Vec::new();
    let mut bench_grad: Vec<u64> = Vec::new();
    let mut wrote = true;
    for name in HASH_FIXTURES {
        match write_fixture_mesh(dir, name) {
            Some((conn, std_conn, grad)) => {
                bench_conn.extend(conn);
                bench_std.push(std_conn);
                bench_grad.push(grad);
            }
            None => wrote = false,
        }
    }
    let write_ok = std::fs::write(&program, ROUTE_PROGRAM).is_ok() && wrote;
    if !write_ok {
        note.push_str("fixture-write-failed+");
    }

    let bench = Hashes {
        connectivity: bench_conn,
        std_acos: bench_std,
        gradient: bench_grad,
        sweep_libm: hash_acos(libm::acos),
        sweep_std: hash_acos(f64::acos),
    };

    let build = |opt: &str, out: &str| -> bool {
        write_ok
            && Command::new("rustc")
                .args([
                    "--edition",
                    "2021",
                    "-C",
                    opt,
                    program.to_str().unwrap_or(""),
                    "-o",
                    dir.join(out).to_str().unwrap_or(""),
                ])
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status()
                .map(|s| s.success())
                .unwrap_or(false)
    };
    let opt3 = build("opt-level=3", "route_o3");
    let opt0 = build("opt-level=0", "route_o0");

    let run_local = |binary: &str| -> Option<Vec<String>> {
        HASH_FIXTURES
            .iter()
            .map(|name| {
                let out = Command::new(dir.join(binary))
                    .arg(dir.join(format!("mesh_{name}.txt")))
                    .output()
                    .ok()?;
                if out.status.success() {
                    String::from_utf8(out.stdout).ok()
                } else {
                    None
                }
            })
            .collect()
    };

    let local_o3 = if opt3 {
        run_local("route_o3").as_deref().and_then(hashes_from_output)
    } else {
        None
    };
    let local_o0 = if opt0 {
        run_local("route_o0").as_deref().and_then(hashes_from_output)
    } else {
        None
    };
    if local_o3.is_none() {
        note.push_str("local-program-build-or-run-failed+");
    }

    let program_matches_bench = local_o3.as_ref() == Some(&bench);
    let port_matches_libm = local_o3
        .as_ref()
        .is_some_and(|l| l.sweep_libm == bench.sweep_libm);
    let codegen_control_ran = local_o3.is_some() && local_o0.is_some();
    let moved_codegen = match (&local_o3, &local_o0) {
        (Some(a), Some(b)) => {
            count_moved(&a.connectivity, &b.connectivity)
                + count_moved(&a.std_acos, &b.std_acos)
                + count_moved(&a.gradient, &b.gradient)
                + usize::from(a.sweep_libm != b.sweep_libm)
                + usize::from(a.sweep_std != b.sweep_std)
        }
        _ => 0,
    };
    let local = local_o3.unwrap_or_else(|| bench.clone());

    // ── the second machine ───────────────────────────────────────────────────
    let probe = ssh(
        "uname -m; sysctl -n machdep.cpu.brand_string 2>/dev/null || uname -s; rustc --version",
    );
    let (second_machine, second_toolchain, remote) = match probe {
        None => {
            note.push_str("ssh-unreachable+");
            (String::from("absent"), String::from("absent"), None)
        }
        Some(text) => {
            let lines: Vec<&str> = text.lines().collect();
            let arch = lines.first().copied().unwrap_or("unknown");
            let brand = lines.get(1).copied().unwrap_or("unknown");
            let rustc = lines.get(2).copied().unwrap_or("");
            let machine = format!("{}-{}", slug(brand), slug(arch));
            if !rustc.starts_with("rustc") {
                note.push_str("no-rust-toolchain-on-second-machine+");
                (machine, String::from("absent"), None)
            } else {
                let toolchain = version_only(rustc);
                let mut args: Vec<String> = vec![
                    String::from("-cz"),
                    String::from("-C"),
                    dir.to_string_lossy().to_string(),
                    String::from("p73_route.rs"),
                ];
                args.extend(HASH_FIXTURES.iter().map(|n| format!("mesh_{n}.txt")));
                let payload = Command::new("tar")
                    .args(&args)
                    .output()
                    .ok()
                    .filter(|o| o.status.success())
                    .map(|o| o.stdout)
                    .unwrap_or_default();
                if payload.is_empty() {
                    note.push_str("local-tar-failed+");
                    (machine, toolchain, None)
                } else if send_payload(&payload) {
                    let runs = HASH_FIXTURES
                        .iter()
                        .map(|n| format!("./route_o3 mesh_{n}.txt"))
                        .collect::<Vec<_>>()
                        .join(" && ");
                    let script = format!(
                        "cd ~/p73_crossmachine && rustc --edition 2021 -C opt-level=3 \
                         p73_route.rs -o route_o3 >/dev/null 2>&1 && {runs}"
                    );
                    let parsed = ssh(&script).as_deref().and_then(split_fixture_outputs);
                    if parsed.is_none() {
                        note.push_str("remote-build-or-run-failed+");
                    }
                    let _ = ssh("rm -rf ~/p73_crossmachine");
                    (machine, toolchain, parsed)
                } else {
                    note.push_str("payload-transfer-failed+");
                    (machine, toolchain, None)
                }
            }
        }
    };

    let (moved_connectivity, moved_gradient, moved_std_acos, moved_sweep_libm, moved_sweep_std) =
        match &remote {
            Some(r) => (
                count_moved(&local.connectivity, &r.connectivity),
                count_moved(&local.gradient, &r.gradient),
                count_moved(&local.std_acos, &r.std_acos),
                usize::from(local.sweep_libm != r.sweep_libm),
                usize::from(local.sweep_std != r.sweep_std),
            ),
            None => (0, 0, 0, 0, 0),
        };
    let (remote_sweep_libm, remote_sweep_std) =
        remote.as_ref().map_or((0, 0), |r| (r.sweep_libm, r.sweep_std));
    // Which fixtures the platform-libm arm moved on. The count alone does not
    // say whether the disagreement is one adapter's rounding on one triangle or
    // a systematic difference, and the fixture names are the cheapest handle on
    // that question.
    let moved_std_acos_fixtures = remote.as_ref().map_or_else(
        || String::from("none"),
        |r| {
            let names: Vec<&str> = HASH_FIXTURES
                .iter()
                .enumerate()
                .filter(|(i, _)| local.std_acos.get(*i) != r.std_acos.get(*i))
                .map(|(_, n)| *n)
                .collect();
            if names.is_empty() {
                String::from("none")
            } else {
                names.join("+")
            }
        },
    );

    Crossmachine {
        local,
        reached_second: remote.is_some(),
        second_machine,
        second_toolchain,
        note: if note.is_empty() {
            String::from("ok")
        } else {
            note.trim_end_matches('+').to_string()
        },
        local_toolchain,
        moved_connectivity,
        moved_gradient,
        moved_std_acos,
        moved_std_acos_fixtures,
        moved_sweep_libm,
        moved_sweep_std,
        moved_codegen,
        codegen_control_ran,
        program_matches_bench,
        port_matches_libm,
        remote_sweep_libm,
        remote_sweep_std,
    }
}

// ── the run ─────────────────────────────────────────────────────────────────

/// Per-case aggregates, computed once the pooled decile is known.
struct Agg {
    medians: [f64; 7],
    p99: [f64; 7],
    rho: [f64; 7],
    canary_mean: f64,
    canary_max: f64,
    /// The registered vacuity column: this case's triangles at or below its own
    /// 10th percentile of the radius ratio.
    worst_decile: usize,
    /// The strict reading: at or below the whole fixture's 10th percentile.
    worst_decile_pooled: usize,
    mean_connectivity_median: f64,
    /// C1 can be scored here at all: the denominator is non-zero **and** the
    /// case contains triangles the proposal's claim is about.
    c1_testable: bool,
    /// C2 can be scored here at all: the canary is not identically zero.
    c2_testable: bool,
}

type Row = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-73");

    // ── phase 1: every (field, resolution) ───────────────────────────────────
    println!(
        "{:>14} {:>5} {:>8} {:>8} {:>8} {:>9} {:>9} {:>9}",
        "field", "res", "verts", "tris", "scored", "angle", "cd", "ratio"
    );
    let mut cases: Vec<Case> = Vec::new();
    isomesh::for_each_reference_field!(f64, |name, field| {
        for samples in RESOLUTIONS {
            let case = measure(name, &field, samples);
            let angle = median_of_sorted(&sorted(&case.errors[ANGLE]));
            let cd = median_of_sorted(&sorted(&case.errors[6]));
            println!(
                "{:>14} {:>5} {:>8} {:>8} {:>8} {angle:>9.4} {cd:>9.4} {:>9.3}",
                case.field,
                case.samples,
                case.vertices,
                case.triangles,
                case.scored,
                angle / cd,
            );
            cases.push(case);
        }
    });

    // ── phase 2: the pooled aspect-ratio decile ──────────────────────────────
    let mut pooled: Vec<f64> = Vec::new();
    for case in &cases {
        pooled.extend_from_slice(&case.radius_ratios);
    }
    let pooled_sorted = sorted(&pooled);
    let p10 = quantile_of_sorted(&pooled_sorted, 0.10);
    let pooled_median = median_of_sorted(&pooled_sorted);
    println!(
        "\npooled radius ratio over {} triangles: p10 {p10:.6}, median {pooled_median:.6}",
        pooled_sorted.len()
    );

    // ── phase 3: the cross-machine hashes ────────────────────────────────────
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/p73_crossmachine");
    let cm = crossmachine(&dir);
    println!(
        "\nC3: second machine {} (rustc {}), local rustc {}, note {}",
        cm.second_machine, cm.second_toolchain, cm.local_toolchain, cm.note
    );
    println!(
        "    local  acos sweep: libm {:016x}  std {:016x}",
        cm.local.sweep_libm, cm.local.sweep_std
    );
    if cm.reached_second {
        println!(
            "    remote acos sweep: libm {:016x}  std {:016x}",
            cm.remote_sweep_libm, cm.remote_sweep_std
        );
    }
    println!(
        "    moved: connectivity {}/{}, gradient {}/{}, std-acos {}/{}, \
         sweep libm {}, sweep std {}, codegen {} (ran {})",
        cm.moved_connectivity,
        cm.local.connectivity.len(),
        cm.moved_gradient,
        cm.local.gradient.len(),
        cm.moved_std_acos,
        cm.local.std_acos.len(),
        cm.moved_sweep_libm,
        cm.moved_sweep_std,
        cm.moved_codegen,
        cm.codegen_control_ran
    );

    // ── phase 4: aggregates and clause verdicts ──────────────────────────────
    let aggs: Vec<Agg> = cases
        .iter()
        .map(|case| {
            let mut medians = [0.0; 7];
            let mut p99 = [0.0; 7];
            for arm in 0..7 {
                let s = sorted(&case.errors[arm]);
                medians[arm] = median_of_sorted(&s);
                p99[arm] = quantile_of_sorted(&s, 0.99);
            }
            let rho = std::array::from_fn(|arm| spearman(&case.canary, &case.errors[arm]));
            let n = case.canary.len().max(1) as f64;
            Agg {
                canary_mean: case.canary.iter().sum::<f64>() / n,
                canary_max: case.canary.iter().copied().fold(0.0, f64::max),
                worst_decile: case
                    .radius_ratios
                    .iter()
                    .filter(|q| **q <= case.field_p10)
                    .count(),
                worst_decile_pooled: case.radius_ratios.iter().filter(|q| **q <= p10).count(),
                mean_connectivity_median: medians[0..6].iter().sum::<f64>() / 6.0,
                // A zero denominator makes the ratio 0/0, and a case with no
                // triangle below the dossier's own bar cannot show the effect
                // the proposal claims. Either way the clause is not scored
                // there, and pretending otherwise is P-70's C3.
                c1_testable: medians[6] > 0.0 && case.below_bar > 0,
                c2_testable: case.canary.iter().copied().fold(0.0, f64::max) > 0.0,
                medians,
                p99,
                rho,
            }
        })
        .collect();

    // C1: every weighting on every **testable** case at or above 3x. An
    // untestable case is counted separately and scored VACUOUS rather than
    // allowed to falsify a clause it cannot measure.
    let mut c1_rows_at_bar = 0usize;
    let mut c1_rows_testable = 0usize;
    let mut c1_rows_total = 0usize;
    for agg in &aggs {
        for w in 0..6 {
            c1_rows_total += 1;
            if agg.c1_testable {
                c1_rows_testable += 1;
                if agg.medians[w] / agg.medians[6] >= C1_RATIO_BAR {
                    c1_rows_at_bar += 1;
                }
            }
        }
    }
    let c1_testable_cases = aggs.iter().filter(|a| a.c1_testable).count();
    let c1_holds = c1_rows_testable > 0 && c1_rows_at_bar == c1_rows_testable;

    // "thin_plate and noise_cavity are the worst": rank the eight fields at each
    // resolution by the mean over the six weightings.
    let mut rank_of: HashMap<(&str, u32), usize> = HashMap::new();
    let mut worst_two_ok: HashMap<u32, bool> = HashMap::new();
    for samples in RESOLUTIONS {
        let mut order: Vec<(&str, f64)> = cases
            .iter()
            .zip(aggs.iter())
            .filter(|(c, _)| c.samples == samples)
            .map(|(c, a)| (c.field, a.mean_connectivity_median))
            .collect();
        order.sort_by(|a, b| b.1.total_cmp(&a.1));
        for (i, (name, _)) in order.iter().enumerate() {
            rank_of.insert((name, samples), i + 1);
        }
        let top: Vec<&str> = order.iter().take(2).map(|(n, _)| *n).collect();
        worst_two_ok.insert(
            samples,
            top.contains(&"thin_plate") && top.contains(&"noise_cavity"),
        );
    }

    // C2: the angle weighting's rank correlation, six of eight fields per
    // resolution. An untestable case cannot clear the bar and is counted as a
    // field that did not, which is the direction that cannot flatter the clause.
    let mut c2_fields_at_res: HashMap<u32, usize> = HashMap::new();
    let mut c2_testable_at_res: HashMap<u32, usize> = HashMap::new();
    for samples in RESOLUTIONS {
        let n = cases
            .iter()
            .zip(aggs.iter())
            .filter(|(c, a)| c.samples == samples && a.c2_testable && a.rho[ANGLE] > C2_RHO_BAR)
            .count();
        c2_fields_at_res.insert(samples, n);
        c2_testable_at_res.insert(
            samples,
            cases
                .iter()
                .zip(aggs.iter())
                .filter(|(c, a)| c.samples == samples && a.c2_testable)
                .count(),
        );
    }
    let c2_holds = RESOLUTIONS
        .iter()
        .all(|s| c2_fields_at_res[s] >= C2_FIELDS_NEEDED);

    // C3. `blocked` is not a verdict, it is the absence of one.
    let c3_holds = if cm.reached_second {
        if cm.moved_connectivity >= 1 && cm.moved_gradient == 0 {
            "true"
        } else {
            "false"
        }
    } else {
        "blocked"
    };

    println!(
        "\nC1: {c1_rows_at_bar} of {c1_rows_testable} testable rows at or above {C1_RATIO_BAR}x \
         ({c1_rows_total} rows total, {c1_testable_cases} of 16 cases testable) -> {c1_holds}"
    );
    for samples in RESOLUTIONS {
        println!(
            "C2 at {samples}^3: {} of {} testable fields above rho {C2_RHO_BAR} (bar is \
             {C2_FIELDS_NEEDED} of 8); worst two as registered: {}",
            c2_fields_at_res[&samples], c2_testable_at_res[&samples], worst_two_ok[&samples]
        );
    }
    println!("C3 -> {c3_holds}");

    // ── phase 5: rows ────────────────────────────────────────────────────────
    let mut rows: Vec<Row> = Vec::new();
    for (case, agg) in cases.iter().zip(aggs.iter()) {
        let arms: Vec<(usize, &str)> = (0..6)
            .map(|w| (w, WEIGHTINGS[w]))
            .chain(std::iter::once((6usize, GRADIENT_ROW)))
            .collect();
        for (arm, label) in arms {
            let ratio = agg.medians[arm] / agg.medians[6];
            let mut row: Row = vec![
                ("field", case.field.to_string()),
                ("resolution", case.samples.to_string()),
                ("weighting", label.to_string()),
                (
                    "median_angle_error_deg",
                    format!("{:.6}", agg.medians[arm]),
                ),
                ("p99_angle_error_deg", format!("{:.6}", agg.p99[arm])),
                ("gradient_median_deg", format!("{:.6}", agg.medians[6])),
                ("ratio_to_gradient", format!("{ratio:.6}")),
                ("canary_mean_abs_f", format!("{:.6e}", agg.canary_mean)),
                (
                    "canary_rank_correlation",
                    format!("{:.6}", agg.rho[arm]),
                ),
                ("worst_decile_triangles", agg.worst_decile.to_string()),
                (
                    "hashes_moved_connectivity",
                    cm.moved_connectivity.to_string(),
                ),
                ("hashes_moved_gradient", cm.moved_gradient.to_string()),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_holds.to_string()),
                ("c3_holds", c3_holds.to_string()),
                // ── extras ───────────────────────────────────────────────────
                ("cell_size", format!("{:.9}", case.cell_size)),
                ("vertices", case.vertices.to_string()),
                ("triangles", case.triangles.to_string()),
                ("scored_vertices", case.scored.to_string()),
                (
                    "boundary_vertices_excluded",
                    case.boundary_excluded.to_string(),
                ),
                ("nonmanifold_edges", case.nonmanifold_edges.to_string()),
                ("ref_degenerate_excluded", case.ref_degenerate.to_string()),
                ("cd_degenerate_excluded", case.cd_degenerate.to_string()),
                (
                    "weight_degenerate_excluded",
                    case.weight_degenerate.to_string(),
                ),
                ("skipped_corners", case.skipped_corners.to_string()),
                (
                    "row_ratio_at_least_3",
                    (ratio >= C1_RATIO_BAR).to_string(),
                ),
                ("radius_ratio_p10_pooled", format!("{p10:.6}")),
                ("radius_ratio_median_pooled", format!("{pooled_median:.6}")),
                ("triangles_below_0p15", case.below_bar.to_string()),
                ("ref_fd_median_deg", format!("{:.9}", case.ref_fd_median)),
                ("ref_fd_p99_deg", format!("{:.9}", case.ref_fd_p99)),
                ("ref_fd_over_half_deg", case.ref_fd_over_half.to_string()),
                (
                    "ref_closed_form_median_deg",
                    format!("{:.9}", case.ref_closed_median),
                ),
                (
                    "ref_closed_form_max_deg",
                    format!("{:.9}", case.ref_closed_max),
                ),
                (
                    "reference_matches_extraction",
                    (case.reference_mismatch == 0).to_string(),
                ),
                ("mwaat_crate_status", case.mwaat_status.to_string()),
                (
                    "mwaat_mismatch_vertices",
                    case.mwaat_mismatch.to_string(),
                ),
                ("canary_zero_vertices", case.canary_zero.to_string()),
                ("canary_max_abs_f", format!("{:.6e}", agg.canary_max)),
                (
                    "field_rank_by_error",
                    rank_of
                        .get(&(case.field, case.samples))
                        .copied()
                        .unwrap_or(0)
                        .to_string(),
                ),
                (
                    "field_mean_connectivity_median_deg",
                    format!("{:.6}", agg.mean_connectivity_median),
                ),
                ("c1_rows_at_least_3", c1_rows_at_bar.to_string()),
                ("c1_rows_total", c1_rows_total.to_string()),
                (
                    "c1_worst_two_as_registered",
                    worst_two_ok[&case.samples].to_string(),
                ),
                (
                    "c2_fields_above_bar_at_res",
                    c2_fields_at_res[&case.samples].to_string(),
                ),
                (
                    "c2_testable_fields_at_res",
                    c2_testable_at_res[&case.samples].to_string(),
                ),
                ("c1_rows_testable", c1_rows_testable.to_string()),
                ("c1_testable_cases", c1_testable_cases.to_string()),
                ("case_c1_testable", agg.c1_testable.to_string()),
                ("case_c2_testable", agg.c2_testable.to_string()),
                (
                    "case_c1_verdict",
                    String::from(if !agg.c1_testable {
                        "vacuous"
                    } else if ratio >= C1_RATIO_BAR {
                        "held"
                    } else {
                        "falsified"
                    }),
                ),
                (
                    "case_c2_verdict",
                    String::from(if !agg.c2_testable {
                        "vacuous"
                    } else if agg.rho[ANGLE] > C2_RHO_BAR {
                        "held"
                    } else {
                        "falsified"
                    }),
                ),
                (
                    "ratio_p99_to_gradient",
                    format!("{:.6}", agg.p99[arm] / agg.p99[6]),
                ),
                ("zero_error_vertices", case.zero_error[arm].to_string()),
                (
                    "cd_magnitude_ratio_median",
                    format!("{:.9}", case.cd_mag_ratio_median),
                ),
                (
                    "worst_decile_max_radius_ratio",
                    format!("{:.6}", case.field_p10),
                ),
                (
                    "worst_decile_triangles_below_pooled_p10",
                    agg.worst_decile_pooled.to_string(),
                ),
                ("min_radius_ratio", format!("{:.6}", case.min_radius_ratio)),
                (
                    "hashes_moved_std_acos_fixtures",
                    cm.moved_std_acos_fixtures.clone(),
                ),
                ("second_machine", cm.second_machine.clone()),
                ("second_machine_toolchain", cm.second_toolchain.clone()),
                ("second_machine_note", cm.note.clone()),
                ("local_toolchain", cm.local_toolchain.clone()),
                ("hash_fixtures", HASH_FIXTURES.join("+")),
                (
                    "hashes_connectivity_total",
                    cm.local.connectivity.len().to_string(),
                ),
                (
                    "hashes_gradient_total",
                    cm.local.gradient.len().to_string(),
                ),
                ("hashes_moved_std_acos", cm.moved_std_acos.to_string()),
                (
                    "hashes_std_acos_total",
                    cm.local.std_acos.len().to_string(),
                ),
                ("acos_sweep_moved_libm", cm.moved_sweep_libm.to_string()),
                ("acos_sweep_moved_std", cm.moved_sweep_std.to_string()),
                ("hashes_moved_codegen", cm.moved_codegen.to_string()),
                (
                    "codegen_control_ran",
                    cm.codegen_control_ran.to_string(),
                ),
                (
                    "port_acos_matches_libm",
                    cm.port_matches_libm.to_string(),
                ),
                (
                    "program_matches_bench",
                    cm.program_matches_bench.to_string(),
                ),
            ];
            row.push(("acos_sweep_libm_local", format!("{:016x}", cm.local.sweep_libm)));
            row.push(("acos_sweep_std_local", format!("{:016x}", cm.local.sweep_std)));
            row.push((
                "acos_sweep_libm_remote",
                format!("{:016x}", cm.remote_sweep_libm),
            ));
            row.push((
                "acos_sweep_std_remote",
                format!("{:016x}", cm.remote_sweep_std),
            ));
            rows.push(row);
        }
    }

    common::experiment::run(prereg, |run| {
        for row in &rows {
            run.record(row);
        }
    });

    // ── controls, after the artefact is on disk ───────────────────────────────
    for (case, agg) in cases.iter().zip(aggs.iter()) {
        let at = format!("{} at {}^3", case.field, case.samples);
        assert!(
            case.scored > 0,
            "VOID: {at} scored no vertex, so every arm's median is over an empty set"
        );
        // The registered vacuity control. Under the per-field reading it is
        // `n/10`, so the only thing it can catch is a case that meshed nothing -
        // which is worth catching and is said plainly rather than dressed up as
        // a shape check. The substantive question is one column over.
        assert!(
            agg.worst_decile > 0,
            "VACUOUS: {at} has an empty bottom aspect-ratio decile, so C1 there is over no \
             triangles at all"
        );
        // Not an assertion, deliberately: a field with no badly shaped triangle
        // and a field whose two instruments are both exact do not invalidate the
        // harness, they invalidate the clause on that field. `case_c1_verdict`
        // and `case_c2_verdict` carry that, and the registration's own
        // instruction is to score it VACUOUS and say why.
        if !agg.c1_testable {
            println!(
                "VACUOUS for C1: {at} - gradient median {:.6} deg, {} triangles below \
                 {RADIUS_RATIO_BAR}, worst triangle {:.6}, {} of {} vertices at exactly zero \
                 angular error on the angle arm",
                agg.medians[6],
                case.below_bar,
                case.min_radius_ratio,
                case.zero_error[ANGLE],
                case.scored
            );
        }
        if !agg.c2_testable {
            println!(
                "VACUOUS for C2: {at} - max |f(v)| is exactly zero over all {} scored \
                 vertices, so the canary has no ranks",
                case.scored
            );
        }
        assert_eq!(
            case.reference_mismatch, 0,
            "M-289: {at} computed a reference normal that is not the one the extractor stored, \
             so this harness is measuring a private idea of the gradient"
        );
        assert_ne!(
            case.mwaat_status, "mismatch",
            "{at}: the area-weighted arm is not the crate's AreaWeightedFaces ({} vertices \
             differ), so the six weightings are anchored to nothing the crate already ships",
            case.mwaat_mismatch
        );
        assert!(
            case.ref_fd_median < REF_MEDIAN_BAR_DEG,
            "M-289: {at} has a reference gradient whose MEDIAN angle to an independent \
             fourth-order difference is {:.6} degrees, above {REF_MEDIAN_BAR_DEG}. Half the \
             population disagreeing is M-289's signature, and no error number in this CSV is \
             trustworthy until it is explained",
            case.ref_fd_median
        );
        assert!(
            case.ref_closed_median < 0.0 || case.ref_closed_median < REF_MEDIAN_BAR_DEG,
            "M-289: {at} has a reference gradient whose MEDIAN angle to an independently \
             written closed form is {:.6} degrees",
            case.ref_closed_median
        );
    }
    // The experiment as a whole must have an instrument. A per-field vacuity is
    // a finding; sixteen of them would be M-44's vacuous zero with a CSV
    // attached.
    assert!(
        c1_testable_cases > 0,
        "VACUOUS: not one of the sixteen cases can score C1 - every case either has a \
         central-difference median of exactly zero or contains no triangle below \
         {RADIUS_RATIO_BAR}, so the whole clause is measuring well-shaped triangles only"
    );
    assert!(
        RESOLUTIONS
            .iter()
            .all(|s| c2_testable_at_res[s] >= C2_FIELDS_NEEDED),
        "VACUOUS: fewer than {C2_FIELDS_NEEDED} fields per resolution have a canary with any \
         variance at all, so C2's own bar of {C2_FIELDS_NEEDED} of eight cannot be reached \
         whatever the correlations are"
    );
    assert!(
        cm.port_matches_libm,
        "C3 is not measurable as designed: the generated program's transcribed acos does not \
         reproduce libm::acos over the sweep, so it does not compute what the crate would"
    );
    assert!(
        cm.program_matches_bench,
        "C3 is not measurable as designed: the generated program's hashes disagree with the \
         ones this bench computed through the crate, so the two implementations of the six \
         weightings have drifted and a cross-machine difference could not be attributed"
    );
    assert_eq!(
        cm.moved_codegen, 0,
        "C3's codegen control moved: the same program at opt-level=0 and opt-level=3 on THIS \
         machine produced different hashes, so a difference against the second machine cannot \
         be attributed to its architecture"
    );
    assert!(
        cm.codegen_control_ran,
        "C3's codegen control did not run, so a hash that agreed across machines agreed with \
         no evidence that it was sensitive to codegen at all"
    );
}
