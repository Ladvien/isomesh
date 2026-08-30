//! **P-83 — mass, centre of mass and inertia from the triangles already emitted.**
//!
//! Ticket: R-083. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p83
//! ```
//!
//! Writes `docs/experiments/p-83.csv`.
//!
//! # What is being compared, and against what
//!
//! [`isomesh::mass::mass_properties`] is Hartmann & Ewougsi Tekeu's surface
//! integral (`10.1007/s00707-025-04419-1`) over the triangles Marching Cubes
//! just produced. The reference is a **dense voxel integration of the field
//! itself**, at 4–16× the method's resolution. So the difference the harness
//! reports is *not* the integrator's quadrature error — that is exact for a
//! polyhedron — it is **how far the extracted polyhedron is from the solid the
//! field defines**. C1's `h²` is a statement about Marching Cubes' geometry seen
//! through a surface integral, which is why a convergence order below 1.5 would
//! falsify it: it would mean the surface integral is not tracking the geometry
//! the volume integral sees.
//!
//! # Every field is clipped to a box, and that is one path rather than a special
//! case for `fbm_terrain`
//!
//! Mass properties exist only for a **closed** solid, and
//! [`ReferenceField::closed_in_domain`] is `false` for `fbm_terrain`: *"a
//! heightfield exits through the sides"*. So every field is intersected with an
//! axis-aligned box inset inside its own domain. For the seven closed fields the
//! box is placed strictly outside the surface, so the intersection is
//! geometrically inert — and the harness **asserts** that, by extracting with
//! and without the clip and comparing the mesh hash (`clip_inert`). For
//! `fbm_terrain` the clip is what closes the solid, and `clip_inert` is `false`
//! there, which is the column proving the control could have failed.
//!
//! # Controls, each of which could have failed
//!
//! - **`boundary_edges` must be zero.** The divergence theorem does not apply to
//!   a surface with a hole; a leaking mesh would make every number here a
//!   different question's answer. Asserted, per row.
//! - **`clip_inert`** — above. Seven `true`, one `false`.
//! - **`reference_order` is the registered vacuity control**, and it is reported
//!   beside **`reference_order_naive`**, which is the same integration with the
//!   boundary voxels classified by their centre's sign instead of by the
//!   tangent-plane volume fraction. Lattice-point counting of a smooth body
//!   converges at roughly `h^1.5` with oscillation, so the naive column is a
//!   *demonstration that this control can report a bad reference* — an `h²` fit
//!   next to a column that does not fit `h²` on the same data is evidence; an
//!   `h²` fit on its own is a hope.
//! - **Every relative error must be strictly positive.** A zero would mean the
//!   method and the reference are the same computation wearing two names.
//!
//! # The reference integrator
//!
//! Dense voxels over the same clipped field and the same domain. Each voxel
//! contributes `fraction × (exact monomial integral over the whole cube)`, so a
//! voxel entirely inside contributes *exactly* — the `h²/12` term is the cube's
//! own second moment, and dropping it would leave a midpoint-rule error on the
//! interior that dominates the inertia comparison at every affordable
//! resolution. The fraction comes from the tangent plane `f(c) + ∇f(c)·(x − c)`,
//! clipped against the cube in closed form, which makes the reference second-
//! order rather than the `O(h)` of counting whole voxels. Boundary voxels are
//! subdivided [`REFINE`]³, which costs almost nothing — they are an `O(N²)`
//! subset — and divides the geometric error by `REFINE²`.
//!
//! # `M-280` and `M-281`
//!
//! `share` is a ratio of two timings taken in the same run of the same build, so
//! it survives a governed CPU; the clock and governor are on every row anyway.
//! The absolute `mass_props_ms` and `extract_ms` are medians and are *not* the
//! result.

// `needless_range_loop` fires on every one of the 3×3 tensor walks here, and the
// suggested `iter_mut().enumerate()` rewrite is worse on a fixed 3×3: the index
// pair `(i, j)` is what the mathematics is written in, and two of the loops read
// `[j][j]` and `[k][k]` from rows the iterator is not on. The crate's own
// `mass.rs` and `dual_contouring/solve.rs` make the same call.
#![allow(clippy::needless_range_loop, clippy::too_many_lines)]

mod common;

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use isomesh::fields::{BoxExact, Intersection, ReferenceField};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::mass::{MassProperties, mass_properties};
use isomesh::validate::{ValidateConfig, validate_indexed};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Samples per axis for the extracted mesh. `33` is the registered resolution
/// for C1's `1e-4`; the other two exist so the `h²` fit has three points.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// Base voxels per axis for the reference, coarsest first. Three so the
/// reference's own convergence order can be estimated from it rather than
/// asserted about it.
const REFERENCE_CELLS: [usize; 3] = [64, 128, 256];

/// Boundary voxels are subdivided this many times per axis.
///
/// The reference's error is entirely in the boundary layer once interior cells
/// integrate exactly, and that layer is `O(N²)` cells — so refining it is the
/// cheapest accuracy available. `4` puts the finest arm at an effective 1024³
/// for about 1.3× the cost of 256³.
const REFINE: usize = 4;

/// Slabs the reference grid is cut into for threading.
///
/// Fixed rather than derived from the core count, so the *summation order* is a
/// property of the harness and not of the machine it ran on: each slab is
/// summed independently and the slabs are folded in index order.
const SLABS: usize = 96;

/// Extractions timed per row, median taken.
const EXTRACT_REPS: usize = 5;

/// Mass-property passes timed per row, median taken. More than the extraction
/// because it is two orders of magnitude shorter and the clock's granularity is
/// a bigger share of it.
const MASS_REPS: usize = 25;

/// The ten moments of the solid about the origin, at unit density.
#[derive(Clone, Copy, Debug, Default)]
struct Moments {
    volume: f64,
    first: [f64; 3],
    second: [[f64; 3]; 3],
}

impl Moments {
    fn add(&mut self, other: &Self) {
        self.volume += other.volume;
        for i in 0..3 {
            self.first[i] += other.first[i];
            for j in 0..3 {
                self.second[i][j] += other.second[i][j];
            }
        }
    }

    /// Accumulate one axis-aligned cube of side `h` centred at `c`, of which the
    /// fraction `f` is inside the solid.
    ///
    /// The monomial integrals are the cube's own and are exact: `∫x dV = h³cₓ`
    /// and `∫x² dV = h³(cₓ² + h²/12)`. A voxel wholly inside therefore
    /// contributes with no quadrature error at all, which is what leaves the
    /// reference's whole error in the boundary layer where `REFINE` can reach
    /// it.
    #[inline]
    fn cube(&mut self, c: [f64; 3], h: f64, f: f64) {
        let w = f * h * h * h;
        let k = h * h / 12.0;
        self.volume += w;
        for i in 0..3 {
            self.first[i] += w * c[i];
            for j in 0..3 {
                let m = if i == j { c[i] * c[i] + k } else { c[i] * c[j] };
                self.second[i][j] += w * m;
            }
        }
    }

    /// Centre of mass.
    fn centroid(&self) -> [f64; 3] {
        [
            self.first[0] / self.volume,
            self.first[1] / self.volume,
            self.first[2] / self.volume,
        ]
    }

    /// Inertia tensor about the origin, in the crate's convention:
    /// `Θᵢᵢ = ∫(xⱼ² + xₖ²) dV`, `Θᵢⱼ = −∫xᵢxⱼ dV`.
    fn inertia_about_origin(&self) -> [[f64; 3]; 3] {
        let mut out = [[0.0; 3]; 3];
        for i in 0..3 {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            out[i][i] = self.second[j][j] + self.second[k][k];
            out[i][j] = -self.second[i][j];
            out[j][i] = -self.second[j][i];
        }
        out
    }
}

/// Fraction of the cube of side `h` centred where `f` and `∇f` were sampled that
/// lies in `{f < 0}`, from the tangent plane there.
///
/// With `aᵢ = |∂ᵢf|·h` and `α = −f + ½Σaᵢ`, the wanted volume is
/// `vol{t ∈ [0,1]³ : a·t ≤ α}`, which is the CDF of `Σ aᵢUᵢ` for independent
/// uniforms — an inclusion–exclusion over the subsets of the axes:
///
/// ```text
/// F = ( Σ_{T ⊆ S} (−1)^|T| max(0, α − Σ_T a)^k ) / ( k! ∏_S a )
/// ```
///
/// `S` is the set of axes with a non-negligible `aᵢ` and `k = |S|`. Dropping the
/// negligible ones is not an approximation dodge: an axis-aligned face — every
/// face of `box_exact`, every wall of the clip box — has two components exactly
/// zero, the formula's `k = 3` form is then `0/0`, and the reduced-dimension
/// form is its exact limit. The cut is relative to the largest component, so the
/// error it can introduce is bounded by that ratio.
#[inline]
fn plane_fraction(f: f64, g: [f64; 3], h: f64) -> f64 {
    let a = [g[0].abs() * h, g[1].abs() * h, g[2].abs() * h];
    let total = a[0] + a[1] + a[2];
    let alpha = -f + 0.5 * total;
    if alpha <= 0.0 {
        return 0.0;
    }
    if alpha >= total {
        return 1.0;
    }
    let cut = 1e-9 * a[0].max(a[1]).max(a[2]);
    let mut kept = [0.0f64; 3];
    let mut k = 0usize;
    for &ai in &a {
        if ai > cut {
            kept[k] = ai;
            k += 1;
        }
    }
    let mut num = 0.0;
    for mask in 0..(1u32 << k) {
        let mut offset = 0.0;
        let mut sign = 1.0;
        for (bit, &ai) in kept.iter().take(k).enumerate() {
            if mask >> bit & 1 == 1 {
                offset += ai;
                sign = -sign;
            }
        }
        let d = alpha - offset;
        if d > 0.0 {
            num += sign * d.powi(k as i32);
        }
    }
    let mut den = match k {
        1 => 1.0,
        2 => 2.0,
        _ => 6.0,
    };
    for &ai in kept.iter().take(k) {
        den *= ai;
    }
    (num / den).clamp(0.0, 1.0)
}

/// Dense voxel integration of `{field < 0}` over `[lo, hi]³`.
///
/// Returns `(fraction, naive)`: the tangent-plane volume fraction and, from the
/// same sample, the whole-voxel classification by the centre's sign. The second
/// is the control on the first.
fn reference_moments<S: Sdf<Scalar = f64> + Sync>(
    field: &S,
    lo: f64,
    hi: f64,
    cells: usize,
    boundary_voxels: &AtomicUsize,
) -> (Moments, Moments) {
    let h = (hi - lo) / cells as f64;
    let sub = h / REFINE as f64;
    let next = AtomicUsize::new(0);
    let done: Mutex<BTreeMap<usize, (Moments, Moments)>> = Mutex::new(BTreeMap::new());
    let workers = std::thread::available_parallelism().map_or(1, std::num::NonZero::get);

    std::thread::scope(|scope| {
        for _ in 0..workers {
            let (next, done) = (&next, &done);
            scope.spawn(move || {
                loop {
                    let slab = next.fetch_add(1, Ordering::Relaxed);
                    if slab >= SLABS {
                        break;
                    }
                    let z0 = slab * cells / SLABS;
                    let z1 = (slab + 1) * cells / SLABS;
                    let mut fine = Moments::default();
                    let mut naive = Moments::default();
                    let mut boundary = 0usize;
                    for iz in z0..z1 {
                        for iy in 0..cells {
                            for ix in 0..cells {
                                let c = [
                                    lo + (ix as f64 + 0.5) * h,
                                    lo + (iy as f64 + 0.5) * h,
                                    lo + (iz as f64 + 0.5) * h,
                                ];
                                let v = field.sample(c);
                                naive.cube(c, h, if v < 0.0 { 1.0 } else { 0.0 });
                                let g = field.gradient(c);
                                let f = plane_fraction(v, g, h);
                                if f <= 0.0 || f >= 1.0 {
                                    fine.cube(c, h, f);
                                    continue;
                                }
                                boundary += 1;
                                for sz in 0..REFINE {
                                    for sy in 0..REFINE {
                                        for sx in 0..REFINE {
                                            let offset = |i: usize, s: usize| {
                                                c[i] + (s as f64 + 0.5 - REFINE as f64 * 0.5) * sub
                                            };
                                            let q = [offset(0, sx), offset(1, sy), offset(2, sz)];
                                            let vq = field.sample(q);
                                            let gq = field.gradient(q);
                                            fine.cube(q, sub, plane_fraction(vq, gq, sub));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    boundary_voxels.fetch_add(boundary, Ordering::Relaxed);
                    done.lock().expect("slab map").insert(slab, (fine, naive));
                }
            });
        }
    });

    // Folded in slab order, so the sum does not depend on how the threads raced.
    let mut fine = Moments::default();
    let mut naive = Moments::default();
    for (a, b) in done.into_inner().expect("slab map").values() {
        fine.add(a);
        naive.add(b);
    }
    (fine, naive)
}

/// The floor below which a series carries no convergence information.
///
/// **This constant is a fixture defect the first run's own numbers found.** Five
/// of the eight reference fields are centred on the origin, so their centre of
/// mass is *exactly* zero and the method reproduces it to `1e-55`. Fitting a
/// slope through three numbers that are all round-off gives a slope of zero, and
/// folding that through a `min` turned the best available outcome — the method
/// is exact on this quantity — into a reported convergence order of `0.00` and a
/// falsified C1 on `box_exact`. The same defect made every `reference_order`
/// negative, because the reference's centroid series is the same noise.
///
/// So a series whose values, or whose successive differences, never rise above
/// `1e-10` of their natural scale is **exact, not slow**, and is excluded from
/// the fit and named in its own column rather than silently averaged in.
const INFORMATIVE: f64 = 1e-10;

/// Least-squares slope of `ln(error)` against `ln(h)`: the observed convergence
/// order. `None` when the series is below [`INFORMATIVE`].
fn convergence_order(steps: &[f64], errors: &[f64]) -> Option<f64> {
    if errors.iter().fold(0.0f64, |m, e| m.max(*e)) <= INFORMATIVE {
        return None;
    }
    let n = steps.len() as f64;
    let xs: Vec<f64> = steps.iter().map(|h| h.ln()).collect();
    let ys: Vec<f64> = errors.iter().map(|e| e.ln()).collect();
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut den = 0.0;
    for (x, y) in xs.iter().zip(&ys) {
        num += (x - mx) * (y - my);
        den += (x - mx) * (x - mx);
    }
    Some(num / den)
}

/// Three-grid order estimate for a quantity nobody knows the exact value of:
/// with each grid twice the last, `p = log₂(‖q₁ − q₂‖ / ‖q₂ − q₃‖)`, and the
/// remaining error of the finest is `‖q₃ − q₂‖ / (2^p − 1)`.
///
/// **Over the whole ten-moment vector, not component by component**, and that is
/// the second fixture defect this harness's own output found: taking the
/// *minimum* order over ten scalar series gave `−3.91` on `fbm_terrain` and
/// `−5.75` on `noise_cavity`, because a single near-zero off-diagonal moment
/// changes sign between two grids and its two successive differences are then in
/// the wrong ratio. The norm cannot do that — it is the residual of the whole
/// integration, which is what the vacuity control is asking about.
///
/// `None` when the three grids agree to within [`INFORMATIVE`], which is not a
/// failure to converge but a reference that was already exact — as the voxel
/// integration is on an axis-aligned box.
fn self_convergence(q: &[[f64; 10]; 3]) -> Option<(f64, f64, f64)> {
    let gap = |a: &[f64; 10], b: &[f64; 10]| -> f64 {
        a.iter()
            .zip(b)
            .map(|(x, y)| (x - y) * (x - y))
            .sum::<f64>()
            .sqrt()
    };
    let coarse = gap(&q[0], &q[1]);
    let fine = gap(&q[1], &q[2]);
    if coarse <= INFORMATIVE || fine <= INFORMATIVE {
        return None;
    }
    let order = (coarse / fine).log2();
    // Richardson's `fine/(2^p − 1)` divides by something that goes to zero as
    // the observed order does, so on a sequence that has *stalled* it reports an
    // error larger than the quantity. `fbm_terrain` stalls at `p ≈ 0` and the
    // extrapolation read 1.6e-2 for a last-doubling movement of 1.4e-5. Both go
    // out: the extrapolation where the order supports it, and the raw gap, which
    // assumes nothing.
    let residual = fine / (2.0f64.powf(order) - 1.0).abs();
    Some((order, residual, fine))
}

/// The ten moments, each divided by its own natural scale, so that a
/// three-grid convergence estimate over them is comparing dimensionless numbers
/// of order one and [`INFORMATIVE`] means the same thing for all of them.
fn normalised(m: &Moments, volume: f64, length: f64) -> [f64; 10] {
    let v = volume;
    let a = v * length;
    let b = a * length;
    [
        m.volume / v,
        m.first[0] / a,
        m.first[1] / a,
        m.first[2] / a,
        m.second[0][0] / b,
        m.second[1][1] / b,
        m.second[2][2] / b,
        m.second[0][1] / b,
        m.second[1][2] / b,
        m.second[0][2] / b,
    ]
}

/// Volume, first and second moments of the polyhedron by the classical
/// signed-tetrahedron decomposition: the origin fanned to every triangle,
/// `V = det(a,b,c)/6` and
/// `∫xᵢxⱼ = (V/20)(aᵢaⱼ + bᵢbⱼ + cᵢcⱼ + sᵢsⱼ)` with `s = a + b + c`.
///
/// **This is the instrument that attributes C1's error.** It is a *volume*
/// integral, exact for a closed oriented triangle mesh, so any disagreement
/// with the surface integral is the surface integral's own — and any agreement
/// says the whole gap against the voxel reference is the polyhedron failing to
/// be the solid, which is precisely the distinction the registered falsifier
/// draws when it says "the surface integral is not seeing the geometry the
/// volume integral sees".
fn tetra_moments(positions: &[[f64; 3]], triangles: &[[u32; 3]]) -> Moments {
    let mut out = Moments::default();
    for triangle in triangles {
        let a = positions[triangle[0] as usize];
        let b = positions[triangle[1] as usize];
        let c = positions[triangle[2] as usize];
        let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
        let v = det / 6.0;
        let s = [a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2]];
        // Hoisted out of the 3x3, as the surface integral's own factors are: nine
        // divisions per triangle here would have made the speed ratio below a
        // measurement of this control's sloppiness rather than of the two
        // formulations.
        let quarter = v * 0.25;
        let twentieth = v * 0.05;
        out.volume += v;
        for i in 0..3 {
            out.first[i] += quarter * s[i];
            for j in 0..3 {
                out.second[i][j] +=
                    twentieth * (a[i] * a[j] + b[i] * b[j] + c[i] * c[j] + s[i] * s[j]);
            }
        }
    }
    out
}

/// FNV-1a over the bit patterns, so a tensor can be compared across machines
/// without a hashing dependency and without printing eighty digits.
fn digest(values: &[f64]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for v in values {
        for byte in v.to_bits().to_le_bytes() {
            h ^= u64::from(byte);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    h
}

fn mesh_digest(mesh: &MeshBuffer<f64>) -> u64 {
    let mut flat: Vec<f64> = Vec::with_capacity(mesh.positions.len() * 3 + mesh.indices.len());
    for p in &mesh.positions {
        flat.extend_from_slice(p);
    }
    for i in &mesh.indices {
        flat.push(f64::from(*i));
    }
    digest(&flat)
}

fn tensor_digest(props: &MassProperties<f64>) -> u64 {
    let mut flat = vec![props.volume];
    flat.extend_from_slice(&props.center_of_mass);
    for row in &props.inertia {
        flat.extend_from_slice(row);
    }
    for row in &props.inertia_about_origin {
        flat.extend_from_slice(row);
    }
    digest(&flat)
}

/// Largest absolute entry of a 3×3.
fn tensor_scale(t: &[[f64; 3]; 3]) -> f64 {
    let mut m = 0.0f64;
    for row in t {
        for v in row {
            m = m.max(v.abs());
        }
    }
    m
}

fn tensor_gap(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> f64 {
    let mut m = 0.0f64;
    for i in 0..3 {
        for j in 0..3 {
            m = m.max((a[i][j] - b[i][j]).abs());
        }
    }
    m
}

fn median(mut xs: Vec<f64>) -> f64 {
    xs.sort_by(f64::total_cmp);
    xs[xs.len() / 2]
}

/// The clip box's half extent, per field.
///
/// Data, not a code path: every field takes the same `Intersection`, and this is
/// the one number that differs. Each is between the field's own surface and the
/// wall of its domain, except `fbm_terrain`, which has no surface to stay
/// outside of — its solid *is* the box, lidded by the terrain.
fn clip_half_extent(name: &str) -> f64 {
    match name {
        // domain 2, surfaces inside |x| ≤ 1.3
        "sphere" | "torus" | "box_exact" | "csg_difference" | "thin_plate" => 1.5,
        // domain 2, capped by a sphere of radius 1.5
        "noise_cavity" => 1.75,
        // domain 7, capped by a sphere of radius 6
        "gyroid" => 6.5,
        // domain 8, terrain amplitude bounded well inside 5.875, and 5.875 is on
        // no grid plane of a 16-, 32- or 64-cell axis, so no sample lands
        // exactly on the lid.
        "fbm_terrain" => 5.875,
        other => panic!("P-83: no clip half extent registered for reference field {other}"),
    }
}

/// Machine clock, as `M-280` requires on the row rather than in someone's
/// memory.
fn clock() -> (String, String) {
    let read = |p: &str| std::fs::read_to_string(p).map(|s| s.trim().to_string());
    let mhz = read("/sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .map_or_else(
            || String::from("unknown"),
            |khz| format!("{:.0}", khz / 1000.0),
        );
    let governor = read("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor")
        .unwrap_or_else(|_| String::from("unknown"));
    (mhz, governor)
}

/// Hashes measured on the second machine, if they have been brought over.
///
/// Lines are `field,resolution,inertia_hash,mesh_hash`. Absent file means the
/// cross-machine half of C2 was not measured on this run, and the column says
/// `no_peer` rather than guessing.
fn peer_hashes() -> BTreeMap<(String, u32), (String, String)> {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/experiments/p-83-m5-hashes.txt");
    let mut out = BTreeMap::new();
    let Ok(text) = std::fs::read_to_string(&path) else {
        return out;
    };
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() != 4 {
            continue;
        }
        let Ok(resolution) = parts[1].parse::<u32>() else {
            continue;
        };
        out.insert(
            (parts[0].to_string(), resolution),
            (parts[2].to_string(), parts[3].to_string()),
        );
    }
    out
}

type Row = Vec<(&'static str, String)>;

/// Everything measured for one field at one resolution, before the per-field
/// fits are folded in.
struct Arm {
    resolution: u32,
    step: f64,
    triangles: usize,
    vertices: usize,
    boundary_edges: u64,
    non_manifold_edges: u64,
    clip_inert: bool,
    props: MassProperties<f64>,
    volume_rel: f64,
    com_rel: f64,
    inertia_rel: f64,
    inertia_com_rel: f64,
    asymmetry_rel: f64,
    asymmetry_post: f64,
    extract_ms: f64,
    mass_ms: f64,
    inertia_hash: u64,
    mesh_hash: u64,
    tetra_rel: f64,
    tetra_ms: f64,
}

fn main() {
    if !std::env::args().any(|a| a == "--bench") {
        return;
    }

    let peers = peer_hashes();
    let (mhz, governor) = clock();

    common::experiment::run(isomesh::experiment!("P-83"), |run| {
        let mut rows: Vec<Row> = Vec::new();

        isomesh::for_each_reference_field!(f64, |name, field| {
            let (lo, hi) = field.domain();
            let closed = field.closed_in_domain();
            let clip = clip_half_extent(name);
            let solid = Intersection {
                a: &field,
                b: BoxExact {
                    center: [0.0; 3],
                    half_extents: [clip; 3],
                },
            };

            // ── the reference, and its own convergence ──────────────────────
            let boundary_voxels = AtomicUsize::new(0);
            let mut fine = Vec::new();
            let mut naive = Vec::new();
            for cells in REFERENCE_CELLS {
                let (a, b) = reference_moments(&solid, lo[0], hi[0], cells, &boundary_voxels);
                fine.push(a);
                naive.push(b);
            }
            // The reference must have had a boundary to get wrong. A zero here
            // would mean the plane-fraction path never ran and the "dense voxel
            // integration" was whole-voxel counting under another name.
            assert!(
                boundary_voxels.load(Ordering::Relaxed) > 0,
                "P-83/{name}: no partial voxels, so the reference's fraction path never ran"
            );

            let reference = fine[REFERENCE_CELLS.len() - 1];
            let reference_centroid = reference.centroid();
            let reference_inertia = reference.inertia_about_origin();
            // The characteristic length the centroid error is relative to. The
            // centroid itself will not do: five of the eight fields are centred
            // on the origin and a relative error against zero is not a number.
            let length = reference.volume.cbrt();

            // The registered vacuity control, and its own control beside it: the
            // same ten moments from the same samples, classified by the tangent
            // plane and by the voxel centre's sign.
            let fit = |grids: &[Moments]| -> (Option<f64>, f64, f64) {
                let rows = [
                    normalised(&grids[0], reference.volume, length),
                    normalised(&grids[1], reference.volume, length),
                    normalised(&grids[2], reference.volume, length),
                ];
                self_convergence(&rows).map_or((None, 0.0, 0.0), |(p, r, g)| (Some(p), r, g))
            };
            let (order_ref, residual_ref, gap_ref) = fit(&fine);
            let (order_naive, residual_naive, gap_naive) = fit(&naive);

            // An independent check on the whole reference, where a closed form
            // exists. `reference_residual_rel` is the integration's estimate of
            // its own error; this is its actual error, and the two agreeing is
            // what makes the estimate believable on the five fields that have no
            // closed form.
            let analytic = match name {
                "sphere" => Some(4.0 * core::f64::consts::PI / 3.0),
                "torus" => Some(2.0 * core::f64::consts::PI * core::f64::consts::PI * 0.09),
                "box_exact" => Some(8.0),
                _ => None,
            };
            let analytic_rel = analytic.map_or_else(
                || String::from("na"),
                |v| format!("{:.6e}", (reference.volume - v).abs() / v),
            );

            // ── the method, at three resolutions ────────────────────────────
            let mut arms: Vec<Arm> = Vec::new();
            let mut mesh = MeshBuffer::<f64>::new();
            let mut bare = MeshBuffer::<f64>::new();
            for resolution in RESOLUTIONS {
                let shape = RuntimeShape3::new([resolution; 3]).expect("grid");
                let cells = f64::from(resolution - 1);
                let step = (hi[0] - lo[0]) / cells;
                let mut mc = MarchingCubes::<f64>::new();

                let mut extract_ms = Vec::with_capacity(EXTRACT_REPS);
                for _ in 0..EXTRACT_REPS {
                    mesh.reset();
                    let start = Instant::now();
                    mc.extract(&solid, &shape, lo, step, &mut mesh)
                        .expect("extract");
                    extract_ms.push(start.elapsed().as_secs_f64() * 1e3);
                }

                // The clip is inert on a field that was already closed, and the
                // proof is the same mesh, not the same triangle count.
                bare.reset();
                mc.extract(&field, &shape, lo, step, &mut bare)
                    .expect("extract");
                let clip_inert = mesh_digest(&mesh) == mesh_digest(&bare);
                assert_eq!(
                    clip_inert, closed,
                    "P-83/{name} at {resolution}: clip changed the mesh iff the field was open"
                );

                let cfg = ValidateConfig::from_cell_size(step).expect("cell size");
                let report = validate_indexed(&mesh.positions, &mesh.indices, &cfg);
                assert_eq!(
                    report.boundary_edges, 0,
                    "P-83/{name} at {resolution}: the surface leaks, so the divergence \
                     theorem does not apply to it"
                );

                let (triangles, rest) = mesh.indices.as_chunks::<3>();
                assert!(rest.is_empty(), "P-83/{name}: ragged index buffer");

                let props = mass_properties(&mesh.positions, triangles).expect("mass properties");
                let mut mass_ms = Vec::with_capacity(MASS_REPS);
                for _ in 0..MASS_REPS {
                    let start = Instant::now();
                    let p = mass_properties(std::hint::black_box(&mesh.positions), triangles)
                        .expect("mass properties");
                    std::hint::black_box(p.volume);
                    mass_ms.push(start.elapsed().as_secs_f64() * 1e3);
                }

                // The attribution instrument: the same polyhedron, integrated
                // over its volume instead of its surface. Exact, so the gap is
                // the surface integral's own error and nothing else's.
                let tetra = tetra_moments(&mesh.positions, triangles);
                let mut tetra_ms = Vec::with_capacity(MASS_REPS);
                for _ in 0..MASS_REPS {
                    let start = Instant::now();
                    let m = tetra_moments(std::hint::black_box(&mesh.positions), triangles);
                    std::hint::black_box(m.volume);
                    tetra_ms.push(start.elapsed().as_secs_f64() * 1e3);
                }
                let tetra_inertia = tetra.inertia_about_origin();
                let tetra_centroid = tetra.centroid();
                let mut tetra_rel = (props.volume - tetra.volume).abs() / tetra.volume.abs();
                for axis in 0..3 {
                    let e = (props.center_of_mass[axis] - tetra_centroid[axis]).abs();
                    tetra_rel = tetra_rel.max(e / length);
                }
                tetra_rel = tetra_rel.max(
                    tensor_gap(&props.inertia_about_origin, &tetra_inertia)
                        / tensor_scale(&tetra_inertia),
                );
                assert!(
                    tetra_rel < 1e-6,
                    "P-83/{name} at {resolution}: the surface integral disagrees with an \
                     exact volume integral over the same polyhedron by {tetra_rel:e}, which \
                     is the integrator's own error and not the mesh's"
                );

                let volume_rel = (props.volume - reference.volume).abs() / reference.volume.abs();
                let com_rel = {
                    let mut d = 0.0;
                    for axis in 0..3 {
                        let e = props.center_of_mass[axis] - reference_centroid[axis];
                        d += e * e;
                    }
                    d.sqrt() / length
                };
                let scale = tensor_scale(&reference_inertia);
                let inertia_rel =
                    tensor_gap(&props.inertia_about_origin, &reference_inertia) / scale;

                // The centred tensor, for the record: it compounds the centroid
                // error and is what a physics engine actually consumes.
                let mut reference_centred = reference_inertia;
                let c = reference_centroid;
                let v = reference.volume;
                for i in 0..3 {
                    let j = (i + 1) % 3;
                    let k = (i + 2) % 3;
                    reference_centred[i][i] -= v * (c[j] * c[j] + c[k] * c[k]);
                    reference_centred[i][j] += v * c[i] * c[j];
                    reference_centred[j][i] += v * c[j] * c[i];
                }
                let inertia_com_rel = tensor_gap(&props.inertia, &reference_centred)
                    / tensor_scale(&reference_centred);

                // Every error must be able to be non-zero, and is.
                assert!(
                    volume_rel > 0.0 && com_rel > 0.0 && inertia_rel > 0.0,
                    "P-83/{name} at {resolution}: an error came out exactly zero, so the \
                     method and the reference are not independent computations"
                );

                let post = {
                    let t = &props.inertia_about_origin;
                    let mut m = 0.0f64;
                    for i in 0..3 {
                        for j in (i + 1)..3 {
                            m = m.max((t[i][j] - t[j][i]).abs());
                        }
                    }
                    m
                };

                arms.push(Arm {
                    resolution,
                    step,
                    triangles: mesh.triangle_count(),
                    vertices: mesh.vertex_count(),
                    boundary_edges: report.boundary_edges,
                    non_manifold_edges: report.non_manifold_edges,
                    clip_inert,
                    props,
                    volume_rel,
                    com_rel,
                    inertia_rel,
                    inertia_com_rel,
                    asymmetry_rel: props.asymmetry / tensor_scale(&props.inertia_about_origin),
                    asymmetry_post: post,
                    extract_ms: median(extract_ms),
                    mass_ms: median(mass_ms),
                    inertia_hash: tensor_digest(&props),
                    mesh_hash: mesh_digest(&mesh),
                    tetra_rel,
                    tetra_ms: median(tetra_ms),
                });
            }

            let steps: Vec<f64> = arms.iter().map(|a| a.step).collect();
            let order_volume = convergence_order(
                &steps,
                &arms.iter().map(|a| a.volume_rel).collect::<Vec<_>>(),
            );
            let order_com =
                convergence_order(&steps, &arms.iter().map(|a| a.com_rel).collect::<Vec<_>>());
            let order_inertia = convergence_order(
                &steps,
                &arms.iter().map(|a| a.inertia_rel).collect::<Vec<_>>(),
            );
            // Only series that carry information. `box_exact`'s centroid is
            // exactly zero and `thin_plate`'s is 1e-15; a slope fitted through
            // those is a slope through round-off, and taking its minimum turned
            // "the method is exact here" into "the method does not converge".
            let order = [order_volume, order_com, order_inertia]
                .into_iter()
                .flatten()
                .fold(f64::INFINITY, f64::min);

            // A series that was already exact is reported as `exact`, not as a
            // zero slope: the two are opposite results and a reader who grabs
            // the obvious column must not be told the wrong one (C5's rule).
            let slope =
                |o: Option<f64>| o.map_or_else(|| String::from("exact"), |v| format!("{v:.4}"));

            for arm in &arms {
                let share = arm.mass_ms / arm.extract_ms;
                let peer = peers.get(&(name.to_string(), arm.resolution));
                let identical = match peer {
                    None => String::from("no_peer"),
                    Some((inertia, mesh_hash)) => format!(
                        "{}",
                        *inertia == format!("{:016x}", arm.inertia_hash)
                            && *mesh_hash == format!("{:016x}", arm.mesh_hash)
                    ),
                };
                let c1 = arm.volume_rel <= 1e-4
                    && arm.com_rel <= 1e-4
                    && arm.inertia_rel <= 1e-4
                    && order >= 1.5;
                let c2 = if identical == "no_peer" {
                    String::from("blocked")
                } else {
                    format!("{}", arm.props.asymmetry > 0.0 && identical == "true")
                };
                rows.push(vec![
                    ("field", name.to_string()),
                    ("resolution", arm.resolution.to_string()),
                    ("volume", format!("{:.12e}", arm.props.volume)),
                    ("volume_reference", format!("{:.12e}", reference.volume)),
                    ("volume_rel_error", format!("{:.6e}", arm.volume_rel)),
                    ("com_rel_error", format!("{:.6e}", arm.com_rel)),
                    ("inertia_rel_error", format!("{:.6e}", arm.inertia_rel)),
                    ("convergence_order", format!("{order:.4}")),
                    ("asymmetry_pre", format!("{:.6e}", arm.asymmetry_rel)),
                    ("asymmetry_post", format!("{:.6e}", arm.asymmetry_post)),
                    ("bit_identical_across_machines", identical),
                    ("mass_props_ms", format!("{:.6}", arm.mass_ms)),
                    ("extract_ms", format!("{:.6}", arm.extract_ms)),
                    ("share", format!("{share:.6}")),
                    ("c1_holds", c1.to_string()),
                    ("c2_holds", c2),
                    ("c3_holds", (share < 0.02).to_string()),
                    // ── extras ──────────────────────────────────────────────
                    ("asymmetry_pre_abs", format!("{:.6e}", arm.props.asymmetry)),
                    ("boundary_edges", arm.boundary_edges.to_string()),
                    ("cell_size", format!("{:.8}", arm.step)),
                    ("clip_half_extent", format!("{clip:.6}")),
                    ("clip_inert", arm.clip_inert.to_string()),
                    ("clock_governor", governor.clone()),
                    ("clock_mhz", mhz.clone()),
                    (
                        "inertia_com_rel_error",
                        format!("{:.6e}", arm.inertia_com_rel),
                    ),
                    ("inertia_hash", format!("{:016x}", arm.inertia_hash)),
                    (
                        "inertia_scale",
                        format!("{:.12e}", tensor_scale(&reference_inertia)),
                    ),
                    ("mesh_hash", format!("{:016x}", arm.mesh_hash)),
                    ("non_manifold_edges", arm.non_manifold_edges.to_string()),
                    ("order_com", slope(order_com)),
                    ("order_inertia", slope(order_inertia)),
                    ("order_volume", slope(order_volume)),
                    ("reference_cells", REFERENCE_CELLS[2].to_string()),
                    ("reference_order", slope(order_ref)),
                    ("reference_order_naive", slope(order_naive)),
                    ("reference_refine", REFINE.to_string()),
                    ("reference_gap_rel", format!("{gap_ref:.6e}")),
                    ("reference_gap_rel_naive", format!("{gap_naive:.6e}")),
                    ("reference_residual_rel", format!("{residual_ref:.6e}")),
                    ("reference_vs_analytic_rel", analytic_rel.clone()),
                    (
                        "reference_residual_rel_naive",
                        format!("{residual_naive:.6e}"),
                    ),
                    ("tetra_ms", format!("{:.6}", arm.tetra_ms)),
                    ("tetra_rel_error", format!("{:.6e}", arm.tetra_rel)),
                    (
                        "surface_over_tetra_ms",
                        format!("{:.4}", arm.mass_ms / arm.tetra_ms),
                    ),
                    ("triangles", arm.triangles.to_string()),
                    ("vertices", arm.vertices.to_string()),
                ]);
            }
            println!(
                "{name:>14}  clip {clip:>6.3}  ref order {:>6} (naive {:>6}) gap {gap_ref:.2e} \
                 (naive {gap_naive:.2e})  analytic {analytic_rel}  method order {order:5.2}",
                slope(order_ref),
                slope(order_naive)
            );
            for arm in &arms {
                println!(
                    "                {:>3}³  V {:.3e}  com {:.3e}  I {:.3e}  tetra {:.2e}  \
                     asym {:.2e}  share {:.4}  tris {}",
                    arm.resolution,
                    arm.volume_rel,
                    arm.com_rel,
                    arm.inertia_rel,
                    arm.tetra_rel,
                    arm.asymmetry_rel,
                    arm.mass_ms / arm.extract_ms,
                    arm.triangles
                );
            }
        });

        for row in rows {
            run.record(&row);
        }
    });
}
