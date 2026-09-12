//! **P-180 — the Gaussian kinematic formula grades `noise_cavity` once the
//! sphere's cap is cut out of the census box.**
//!
//! Ticket: `R-185`. Pre-registered before this harness existed; minted by the
//! discovery loop (`docs/research/2026-09-12-discovery-loop-prompt.md`,
//! iteration 6) from `M-490`'s one falsified clause.
//!
//! Writes `docs/experiments/p-180.csv`.
//!
//! # What `M-490` left open
//!
//! `M-490` (`P-178`) put `fbm_terrain` inside the formula's matched two-sigma
//! band at 0.25 σ and `noise_cavity` outside it at 6.77 σ with the **sign** of χ
//! wrong: **15** measured against **−6.602285203** predicted. Both fields are
//! hash noise, so Gaussianity alone does not separate them. The formula has a
//! second hypothesis: it is stated for a field that is *stationary* over the set
//! it is integrated on, and `noise_cavity` is
//! `Intersection(NoiseVolume, Sphere(radius 1.5))` on `cube_domain(2.0)` — so
//! outside radius `1.5` the sampled field is the sphere's distance function and
//! not noise at all. The sphere's inscribed cube has half-side
//! `1.5 / √3 = 0.866025404`; inside it the cap term is at most `−0.633974596`
//! and the intersection's `max` is decided by the noise.
//!
//! This row censuses that inscribed, cap-free sub-box with the sub-box's own
//! `λ`, the sub-box's own Lipschitz–Killing curvatures, and a matched Gaussian
//! ensemble on the same sub-box, and asks whether the miss was the cap.
//!
//! # The instrument is `M-490`'s, copied
//!
//! Every function below `constants` up to `sample_field` is `experiment_p178.rs`'s,
//! **verbatim**: the Ohser–Nagel–Schladitz block census, the Hermite recurrence,
//! the closed-form `ρ_j`, the Steiner fit of a box's `L_j`, the exactly-Gaussian
//! Fourier-shell fixture, and the central-difference `μ₂`. A re-derivation would
//! measure a different instrument. Two things are *not* copied:
//!
//! - `rho_quadrature` and the assert that the closed form agrees with it to
//!   `1e-10` — `M-490` already asserted that on the same four `ρ_j` and the
//!   same rungs, and this file reads only `u = 0`;
//! - `sample_field` takes an explicit box rather than the field's own
//!   `domain()`, because the sub-box is the whole point.
//!
//! The copy is not trusted. Two arms re-take known numbers before any
//! prediction is read: the **calibration** arm reads surface χ **2** on `sphere`
//! and **0** on `torus` (`M-458`'s numbers), and the **full-box** arm must
//! reproduce `M-490`'s own rows — `noise_cavity` χ **15** and `fbm_terrain` χ
//! **1** at 65³ — or the copy has drifted and nothing else in the file may be
//! read. Both are `VOID:` asserts.
//!
//! # Arms
//!
//! | arm | box | which clause |
//! |---|---|---|
//! | `full_box` | each field's own `domain()`, 65³ | reproduces `M-490`; half of C3 |
//! | `sub_box` | the inscribed cap-free cube, at 29³ (`h`-matched to the full box) and 65³ (finest) | C1, C2, other half of C3 |
//! | `calibration` | `sphere`, `torus` on their own domains | anti-drift |
//!
//! `fbm_terrain` has no cap, so its "own inscribed sub-box" is not defined by the
//! registration. It is taken as the **same fraction** of its half-extent as
//! `noise_cavity`'s — `0.866025404 / 2.0 = 0.433012702` of `8.0`, half-side
//! `3.464101616` — so that both fields lose the same share of their box.
//!
//! # The band
//!
//! As `M-490`'s `matched_gaussian` arm: a Gaussian ensemble on the *same* box,
//! digitised the same bounded way, from the integer mode shell whose closed-form
//! `λ` on that box is nearest the field's measured one; `REALISATIONS`
//! independent draws; `band_sigma` is the ensemble's standard deviation of χ.
//! On a small box the smallest shell (`|k|² = 1`) already has a `λ` above the
//! field's, so the match can be coarse; `matched_ksq` and `matched_lambda` are
//! recorded beside every row so the reader can see how coarse.
//!
//! # The vacuity control, as registered
//!
//! The sub-box band must reject a prediction built from `4λ` **at `u = 0`, and
//! only at `u = 0`** — `M-490` is why the rung is named: at `u = ±1` the
//! formula's `(u² − 1)` factor is zero for every `λ`. Column
//! `band_rejects_wrong_lambda`; the verdict is read on the `sub_box`
//! `noise_cavity` row at 65³.
//!
//! # Determinism
//!
//! One `Rng` (SplitMix64, `common::poly::Rng`) seeded by a constant, drawn in a
//! fixed order — full box, then sub-box, `noise_cavity` before `fbm_terrain`,
//! 29³ before 65³. The shipped fields are deterministic hash noise. Re-running
//! produces byte-identical rows.

#![allow(
    // These loops index parallel arrays by the same integer, and iterators over
    // one of them would hide the correspondence.
    clippy::needless_range_loop,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use isomesh::Sdf;
use isomesh::fields::{FbmTerrain, NoiseVolume, ReferenceField, Sphere, Torus, noise_cavity};

use common::poly::Rng;

// ─── constants ───────────────────────────────────────────────────────────────

/// Common denominator of every Ohser–Nagel–Schladitz weight, as in `P-145`.
const WEIGHT_DENOMINATOR: i64 = 8;

/// The voxel offsets of one `2×2×2` block, in bit order. `P-145`'s constant.
const BLOCK_OFFSETS: [[usize; 3]; 8] = [
    [0, 0, 0],
    [0, 0, 1],
    [0, 1, 0],
    [0, 1, 1],
    [1, 0, 0],
    [1, 0, 1],
    [1, 1, 0],
    [1, 1, 1],
];

/// Independent draws behind every band. `M-490`'s count.
const REALISATIONS: usize = 32;

/// Seed for every draw in this file.
const SEED: u64 = 0x0180_5f2a_9c1d_7e43;

/// Half-side of the cube inscribed in `noise_cavity`'s cap: `1.5 / √3`.
const CAVITY_SUB_HALF: f64 = 0.866_025_404;

/// The share of a half-extent the sub-box keeps: `CAVITY_SUB_HALF / 2.0`.
const SUB_FRACTION: f64 = CAVITY_SUB_HALF / 2.0;

/// `fbm_terrain`'s sub-box half-side: `8.0 · SUB_FRACTION = 3.464101616`
/// (exact — a power-of-two scaling of the fraction).
const TERRAIN_SUB_HALF: f64 = 8.0 * SUB_FRACTION;

/// Samples per axis on the full box.
const FULL_SAMPLES: usize = 65;

/// Samples per axis on the sub-box at the full box's spacing: side `1.732` at
/// `h = 0.0625` is 28 cells.
const SUB_SAMPLES_H_MATCHED: usize = 29;

/// Samples per axis on the sub-box at the finest rung.
const SUB_SAMPLES_FINEST: usize = 65;

// ─── the Ohser–Nagel–Schladitz weights (P-145's derivation, copied) ──────────

/// One cell of the cubical complex, as it appears in a `2×2×2` block.
#[derive(Clone, Copy, Debug)]
struct LocalCell {
    /// Which of the block's eight voxels the cell is incident to.
    voxels: u8,
    /// The cell's dimension; the Euler sum takes `(-1)^dimension`.
    dimension: u32,
    /// How many blocks of the whole lattice share this cell.
    shared: i64,
}

/// The bits of [`BLOCK_OFFSETS`] whose `axis` offset is `side`.
fn mask_where(axis: usize, side: usize) -> u8 {
    let mut mask = 0u8;
    for (bit, offset) in BLOCK_OFFSETS.iter().enumerate() {
        if offset[axis] == side {
            mask |= 1 << bit;
        }
    }
    mask
}

/// The cells of the union-of-closed-voxels complex a block owns `1/shared` of.
///
/// Foreground 26-connected, background 6-connected: Etiene et al.'s model, and
/// the one `M-458` reports `chi_solid_26` from.
fn cells_26() -> Vec<LocalCell> {
    let mut out = Vec::with_capacity(1 + 6 + 12 + 8);
    out.push(LocalCell {
        voxels: u8::MAX,
        dimension: 0,
        shared: 1,
    });
    for axis in 0..3 {
        for side in 0..2 {
            out.push(LocalCell {
                voxels: mask_where(axis, side),
                dimension: 1,
                shared: 2,
            });
        }
    }
    for normal in 0..3 {
        let u = (normal + 1) % 3;
        let v = (normal + 2) % 3;
        for su in 0..2 {
            for sv in 0..2 {
                out.push(LocalCell {
                    voxels: mask_where(u, su) & mask_where(v, sv),
                    dimension: 2,
                    shared: 4,
                });
            }
        }
    }
    for bit in 0..8 {
        out.push(LocalCell {
            voxels: 1 << bit,
            dimension: 3,
            shared: 8,
        });
    }
    out
}

/// The 256 weights, as numerators over [`WEIGHT_DENOMINATOR`]. Nothing quoted.
fn weight_numerators() -> [i64; 256] {
    let cells = cells_26();
    let mut out = [0i64; 256];
    for (configuration, weight) in out.iter_mut().enumerate() {
        let occupied = configuration as u8;
        for cell in &cells {
            if occupied & cell.voxels != 0 {
                let sign = if cell.dimension % 2 == 0 { 1 } else { -1 };
                *weight += sign * (WEIGHT_DENOMINATOR / cell.shared);
            }
        }
    }
    out
}

// ─── the digital object ──────────────────────────────────────────────────────

/// What the sampled box is topologically.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Domain {
    /// A finite box, padded by one empty voxel per side.
    Bounded,
    /// The 3-torus: the last sample on each axis is the first one.
    Periodic,
}

/// `χ` of the digital solid `{occupied}`, by the block census.
///
/// `occupied` is `samples³` in the crate's own `x`-fastest order. The array is
/// extended by one voxel per axis so every lattice vertex is the base of an
/// in-range `2×2×2` block; periodic domains glue the far layer modularly from the
/// *occupancy* rather than from the values, because a sample sitting on the level
/// can flip sign between the two ends of one period.
fn chi_solid(occupied_in: &[bool], samples: usize, domain: Domain, weights: &[i64; 256]) -> i64 {
    let voxels = match domain {
        Domain::Bounded => samples,
        Domain::Periodic => samples - 1,
    };
    let pad = match domain {
        Domain::Bounded => 2usize,
        Domain::Periodic => 1,
    };
    let extent = voxels + pad;
    let blocks = extent - 1;
    let mut occupied = vec![false; extent * extent * extent];
    let base = usize::from(domain == Domain::Bounded);
    for z in 0..voxels {
        for y in 0..voxels {
            for x in 0..voxels {
                if !occupied_in[x + samples * (y + samples * z)] {
                    continue;
                }
                let index = (x + base) + extent * ((y + base) + extent * (z + base));
                occupied[index] = true;
            }
        }
    }
    if domain == Domain::Periodic {
        for z in 0..extent {
            for y in 0..extent {
                for x in 0..extent {
                    if x < voxels && y < voxels && z < voxels {
                        continue;
                    }
                    let from = (x % voxels) + extent * ((y % voxels) + extent * (z % voxels));
                    occupied[x + extent * (y + extent * z)] = occupied[from];
                }
            }
        }
    }

    let stride = [1usize, extent, extent * extent];
    let offsets: [usize; 8] = std::array::from_fn(|bit| {
        let o = BLOCK_OFFSETS[bit];
        o[0] * stride[0] + o[1] * stride[1] + o[2] * stride[2]
    });
    let mut counts = [0u64; 256];
    for z in 0..blocks {
        for y in 0..blocks {
            let row = y * stride[1] + z * stride[2];
            for x in 0..blocks {
                let cell = row + x;
                let mut configuration = 0u8;
                for (bit, offset) in offsets.iter().enumerate() {
                    configuration |= u8::from(occupied[cell + offset]) << bit;
                }
                counts[configuration as usize] += 1;
            }
        }
    }
    let total: i64 = counts
        .iter()
        .zip(weights.iter())
        .map(|(count, weight)| *count as i64 * weight)
        .sum();
    assert_eq!(
        total % WEIGHT_DENOMINATOR,
        0,
        "VOID: the weighted block census came to {total}/{WEIGHT_DENOMINATOR}, which is not an \
         integer chi -- the cell-sharing multiplicities do not account for every cell once"
    );
    total / WEIGHT_DENOMINATOR
}

// ─── the Gaussian kinematic formula ──────────────────────────────────────────

/// The `j`th probabilists' Hermite polynomial at `x`, by its recurrence.
fn hermite(j: usize, x: f64) -> f64 {
    let mut previous = 1.0;
    let mut current = x;
    if j == 0 {
        return previous;
    }
    for n in 1..j {
        let next = (n as f64).mul_add(-previous, x * current);
        previous = current;
        current = next;
    }
    current
}

/// `ρ_j(u)` in closed form, from `d/dx [H_{j-1} e^(-x²/2)] = -H_j e^(-x²/2)`.
///
/// For `j ≥ 1` the integral telescopes to `H_{j-1}(u) e^(-u²/2)`; for `j = 0` it
/// is the Gaussian tail, written here as the complementary error function.
/// `M-490` asserted this against quadrature of the paper's own integral to
/// `1e-10` on every rung this file uses.
fn rho_closed_form(j: usize, u: f64) -> f64 {
    let scale = (2.0 * std::f64::consts::PI).powf(-((j as f64) + 1.0) / 2.0);
    if j == 0 {
        // (2π)^(-1/2) ∫_u^∞ e^(-x²/2) dx = ½ erfc(u/√2).
        return 0.5 * libm::erfc(u / std::f64::consts::SQRT_2);
    }
    scale * hermite(j - 1, u) * (-0.5 * u * u).exp()
}

/// Euclidean Lipschitz–Killing curvatures of a box, fitted from Steiner's tube.
///
/// The paper states `λ(T(M,r)) = Σ_j M_j(M) r^j / j!` with
/// `M_{k-j}(M)/(k-j)! = L_j(M) ω_{k-j}`. For a box of sides `(a,b,c)` the tube
/// volume is the exact polynomial
/// `abc + 2(ab+bc+ca)r + π(a+b+c)r² + (4/3)πr³`, so the four `L_j` are read off
/// its coefficients. The coefficients are **fitted** from four evaluations of the
/// tube volume rather than written down, and the fit is asserted exact.
fn box_lipschitz_killing(sides: [f64; 3]) -> [f64; 4] {
    let [a, b, c] = sides;
    let tube = |r: f64| {
        let linear = 2.0 * (a * b + b * c + c * a) * r;
        let square = std::f64::consts::PI * (a + b + c) * r * r;
        let cube = 4.0 / 3.0 * std::f64::consts::PI * r * r * r;
        a * b * c + linear + square + cube
    };
    // Solve the 4×4 Vandermonde on r = 0, 1, 2, 3 for the tube polynomial's
    // coefficients, then divide each by the unit-ball volume ω_{3-j}.
    let samples = [tube(0.0), tube(1.0), tube(2.0), tube(3.0)];
    // Newton's divided differences give the monomial coefficients exactly for a
    // cubic sampled at four equally spaced points.
    let d1 = [
        samples[1] - samples[0],
        samples[2] - samples[1],
        samples[3] - samples[2],
    ];
    let d2 = [d1[1] - d1[0], d1[2] - d1[1]];
    let d3 = d2[1] - d2[0];
    let coefficients = [
        samples[0],
        d1[0] - d2[0] / 2.0 + d3 / 3.0,
        d2[0] / 2.0 - d3 / 2.0,
        d3 / 6.0,
    ];
    let omega = [
        1.0,
        2.0,
        std::f64::consts::PI,
        4.0 / 3.0 * std::f64::consts::PI,
    ];
    // coefficient of r^(3-j) is L_j * ω_{3-j}.
    let lk = [
        coefficients[3] / omega[3],
        coefficients[2] / omega[2],
        coefficients[1] / omega[1],
        coefficients[0] / omega[0],
    ];
    let expected = [1.0, a + b + c, a * b + b * c + c * a, a * b * c];
    for j in 0..4 {
        assert!(
            (lk[j] - expected[j]).abs() <= 1e-9 * expected[j].max(1.0),
            "VOID: the Steiner fit gives L_{j} = {} where the box's own geometry gives {} -- the \
             tube polynomial and the Minkowski normalisation disagree, so every prediction below \
             is built on the wrong curvatures",
            lk[j],
            expected[j]
        );
    }
    lk
}

/// `E[χ]` of `{f ≥ u}` for a unit-variance Gaussian field with spectral moment
/// `lambda` on a set whose Euclidean Lipschitz–Killing curvatures are `lk`.
///
/// `L_j` in the induced metric is `λ^(j/2)` times the Euclidean one, which is
/// what makes this one line rather than a second geometry.
fn expected_chi(lk: [f64; 4], lambda: f64, u: f64) -> f64 {
    let mut total = 0.0;
    for j in 0..4 {
        total += lk[j] * lambda.powf(j as f64 / 2.0) * rho_closed_form(j, u);
    }
    total
}

// ─── the exactly-Gaussian fixture ────────────────────────────────────────────

/// The integer mode shell, every `k` with `1 ≤ |k|² ≤ max_k_sq`.
fn mode_shell(max_k_sq: i32) -> Vec<[i32; 3]> {
    let limit = (max_k_sq as f64).sqrt().floor() as i32;
    let mut out = Vec::new();
    for x in -limit..=limit {
        for y in -limit..=limit {
            for z in -limit..=limit {
                let square = x * x + y * y + z * z;
                if square >= 1 && square <= max_k_sq {
                    out.push([x, y, z]);
                }
            }
        }
    }
    out
}

/// A periodic, exactly Gaussian, zero-mean unit-variance field on `[0, side]³`.
///
/// `f(x) = N^(-1/2) Σ_j [ξ_j cos(2π k_j·x / side) + η_j sin(2π k_j·x / side)]`
/// with `ξ, η` i.i.d. standard normal. Each mode contributes
/// `E[(ξ cos + η sin)²] = 1`, so the `N^(-1/2)` makes the variance exactly one,
/// and differentiating term by term gives
/// `μ₂ = Var(∂f/∂x₁) = N^(-1) Σ_j (2π k_{j,1} / side)²`. The shell is closed
/// under the cube group, so `Σ k₁² = Σ|k|²/3` exactly and the field is isotropic
/// on the lattice.
struct GaussianField {
    /// The mode set.
    modes: Vec<[i32; 3]>,
    /// Cosine coefficients, one per mode.
    xi: Vec<f64>,
    /// Sine coefficients, one per mode.
    eta: Vec<f64>,
    /// Period.
    side: f64,
}

impl GaussianField {
    /// Draw one realisation.
    fn draw(modes: &[[i32; 3]], side: f64, rng: &mut Rng) -> Self {
        let n = modes.len();
        let mut xi = Vec::with_capacity(n);
        let mut eta = Vec::with_capacity(n);
        for _ in 0..n {
            let (a, b) = box_muller(rng);
            xi.push(a);
            eta.push(b);
        }
        Self {
            modes: modes.to_vec(),
            xi,
            eta,
            side,
        }
    }

    /// The closed-form spectral moment `μ₂`.
    fn lambda(modes: &[[i32; 3]], side: f64) -> f64 {
        let n = modes.len() as f64;
        let sum: f64 = modes
            .iter()
            .map(|k| f64::from(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]))
            .sum();
        let omega = 2.0 * std::f64::consts::PI / side;
        omega * omega * sum / (3.0 * n)
    }

    /// Sample the field on an `n³` grid spanning one period.
    fn sample(&self, n: usize) -> Vec<f64> {
        let scale = 1.0 / (self.modes.len() as f64).sqrt();
        let step = self.side / (n - 1) as f64;
        let mut out = vec![0.0f64; n * n * n];
        for (mode, (xi, eta)) in self.modes.iter().zip(self.xi.iter().zip(self.eta.iter())) {
            let w = [
                2.0 * std::f64::consts::PI * f64::from(mode[0]) / self.side,
                2.0 * std::f64::consts::PI * f64::from(mode[1]) / self.side,
                2.0 * std::f64::consts::PI * f64::from(mode[2]) / self.side,
            ];
            for z in 0..n {
                let pz = w[2] * (z as f64 * step);
                for y in 0..n {
                    let py = w[1] * (y as f64 * step);
                    let row = n * (y + n * z);
                    for x in 0..n {
                        let phase = w[0] * (x as f64 * step) + py + pz;
                        out[row + x] += scale * (xi * phase.cos() + eta * phase.sin());
                    }
                }
            }
        }
        out
    }
}

/// Two independent standard normals, by Box–Muller on the shared SplitMix64.
fn box_muller(rng: &mut Rng) -> (f64, f64) {
    // `next_f64_unit` is `[-1, 1)`; map to `(0, 1]` so the logarithm is finite.
    let u1 = (rng.next_f64_unit() + 1.0) / 2.0;
    let u1 = if u1 <= f64::MIN_POSITIVE { 1.0 } else { u1 };
    let u2 = (rng.next_f64_unit() + 1.0) / 2.0;
    let radius = (-2.0 * u1.ln()).sqrt();
    let angle = 2.0 * std::f64::consts::PI * u2;
    (radius * angle.cos(), radius * angle.sin())
}

// ─── sampling the shipped fields ─────────────────────────────────────────────

/// Sample an SDF on the box `[lo, hi]` at `n³`, returning the values, the
/// spacing and the box's sides.
///
/// `experiment_p178.rs`'s `sample_field` with the box passed in rather than
/// read from `ReferenceField::domain`; the body is otherwise its.
fn sample_field<S: Sdf<Scalar = f64>>(
    field: &S,
    n: usize,
    lo: [f64; 3],
    hi: [f64; 3],
) -> (Vec<f64>, f64, [f64; 3]) {
    let step = (hi[0] - lo[0]) / (n - 1) as f64;
    let sides = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
    let mut out = vec![0.0f64; n * n * n];
    for z in 0..n {
        for y in 0..n {
            let row = n * (y + n * z);
            for x in 0..n {
                let p = [
                    step.mul_add(x as f64, lo[0]),
                    step.mul_add(y as f64, lo[1]),
                    step.mul_add(z as f64, lo[2]),
                ];
                out[row + x] = field.sample(p);
            }
        }
    }
    (out, step, sides)
}

/// Mean and standard deviation of a sample grid.
fn moments(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let variance = values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n;
    (mean, variance.sqrt())
}

/// `μ₂` estimated from central differences of a standardised grid.
///
/// Averaged over the three axes, over interior samples only. This is the paper's
/// own definition — *"the variance of the first-order partial derivatives"* —
/// applied to what was actually sampled, which is the only way to reach it for a
/// field with no closed-form spectrum.
fn lambda_central_difference(values: &[f64], n: usize, step: f64) -> f64 {
    let at = |x: usize, y: usize, z: usize| values[x + n * (y + n * z)];
    let mut total = 0.0;
    let mut count = 0u64;
    for axis in 0..3 {
        for z in 1..n - 1 {
            for y in 1..n - 1 {
                for x in 1..n - 1 {
                    let (mut lo, mut hi) = ([x, y, z], [x, y, z]);
                    lo[axis] -= 1;
                    hi[axis] += 1;
                    let d = (at(hi[0], hi[1], hi[2]) - at(lo[0], lo[1], lo[2])) / (2.0 * step);
                    total += d * d;
                    count += 1;
                }
            }
        }
    }
    total / count as f64
}

// ─── what this row adds ──────────────────────────────────────────────────────

/// Skewness `m₃ / sd³` and excess kurtosis `m₄ / sd⁴ − 3` of a grid.
///
/// Population moments about the grid's own mean, over every sample. On a
/// standardised grid `sd` is `1` up to rounding, and the division is kept so the
/// two numbers are the dimensionless ones C2's bars are written on regardless.
fn higher_moments(values: &[f64]) -> (f64, f64) {
    let (mean, sd) = moments(values);
    let n = values.len() as f64;
    let mut m3 = 0.0;
    let mut m4 = 0.0;
    for v in values {
        let d = v - mean;
        let d2 = d * d;
        m3 += d2 * d;
        m4 += d2 * d2;
    }
    m3 /= n;
    m4 /= n;
    (m3 / (sd * sd * sd), m4 / (sd * sd * sd * sd) - 3.0)
}

/// The share of samples on the box at which `noise_cavity`'s `max` picks the
/// sphere over the noise — where the sampled field is the cap's distance
/// function rather than noise.
///
/// This is the number that says whether the sub-box is what the registration
/// claims: near zero on the inscribed cube, large on the full box.
fn cap_dominated_share(n: usize, lo: [f64; 3], hi: [f64; 3]) -> f64 {
    let noise = NoiseVolume::<f64>::canonical();
    let cap = Sphere::<f64> {
        center: [0.0; 3],
        radius: 1.5,
    };
    let step = (hi[0] - lo[0]) / (n - 1) as f64;
    let mut dominated = 0u64;
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let p = [
                    step.mul_add(x as f64, lo[0]),
                    step.mul_add(y as f64, lo[1]),
                    step.mul_add(z as f64, lo[2]),
                ];
                if cap.sample(p) >= noise.sample(p) {
                    dominated += 1;
                }
            }
        }
    }
    dominated as f64 / (n * n * n) as f64
}

/// A matched Gaussian band on one box, as `M-490`'s `matched_gaussian` arm.
struct Band {
    /// The mode shell's `|k|²` ceiling.
    ksq: i32,
    /// The shell's closed-form `λ` on this box.
    lambda: f64,
    /// Ensemble mean of χ at `u = 0`.
    mean: f64,
    /// Ensemble standard deviation of χ; the band's half-width is twice this.
    sigma: f64,
    /// The formula's prediction for the ensemble's own `λ`.
    predicted: f64,
    /// The formula's prediction from `4λ` — the control's wrong answer.
    wrong: f64,
}

impl Band {
    /// Whether the band tells the ensemble's `λ` from `4λ` at `u = 0`.
    fn rejects_wrong_lambda(&self) -> bool {
        (self.mean - self.wrong).abs() > 2.0 * self.sigma
    }
}

/// Draw `REALISATIONS` Gaussian fields on the box, from the shell whose
/// closed-form `λ` is nearest `lambda`, and census each one bounded at `u = 0`.
fn matched_band(
    lambda: f64,
    sides: [f64; 3],
    lk: [f64; 4],
    n: usize,
    rng: &mut Rng,
    weights: &[i64; 256],
) -> Band {
    let mut best = (1, f64::INFINITY);
    for candidate in 1..=64 {
        let shell = mode_shell(candidate);
        if shell.is_empty() {
            continue;
        }
        let l = GaussianField::lambda(&shell, sides[0]);
        let gap = (l - lambda).abs();
        if gap < best.1 {
            best = (candidate, gap);
        }
    }
    let shell = mode_shell(best.0);
    let matched_lambda = GaussianField::lambda(&shell, sides[0]);
    let mut matched: Vec<f64> = Vec::with_capacity(REALISATIONS);
    for _ in 0..REALISATIONS {
        let field = GaussianField::draw(&shell, sides[0], rng);
        let sample = field.sample(n);
        let occupied: Vec<bool> = sample.iter().map(|v| *v > 0.0).collect();
        matched.push(chi_solid(&occupied, n, Domain::Bounded, weights) as f64);
    }
    let count = matched.len() as f64;
    let mean = matched.iter().sum::<f64>() / count;
    let variance = matched.iter().map(|c| (c - mean) * (c - mean)).sum::<f64>() / (count - 1.0);
    Band {
        ksq: best.0,
        lambda: matched_lambda,
        mean,
        sigma: variance.sqrt(),
        predicted: expected_chi(lk, matched_lambda, 0.0),
        wrong: expected_chi(lk, 4.0 * matched_lambda, 0.0),
    }
}

// ─── rows ────────────────────────────────────────────────────────────────────

/// Which shipped field a row is about.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Field {
    NoiseCavity,
    FbmTerrain,
}

impl Field {
    fn name(self) -> &'static str {
        match self {
            Self::NoiseCavity => "noise_cavity",
            Self::FbmTerrain => "fbm_terrain",
        }
    }

    /// The field's own `domain()`.
    fn full_half(self) -> f64 {
        match self {
            Self::NoiseCavity => 2.0,
            Self::FbmTerrain => 8.0,
        }
    }

    /// The inscribed sub-box's half-side.
    fn sub_half(self) -> f64 {
        match self {
            Self::NoiseCavity => CAVITY_SUB_HALF,
            Self::FbmTerrain => TERRAIN_SUB_HALF,
        }
    }

    fn sample(self, n: usize, half: f64) -> (Vec<f64>, f64, [f64; 3]) {
        let (lo, hi) = ([-half; 3], [half; 3]);
        match self {
            Self::NoiseCavity => sample_field(&noise_cavity::<f64>(), n, lo, hi),
            Self::FbmTerrain => sample_field(&FbmTerrain::<f64>::canonical(), n, lo, hi),
        }
    }
}

/// Everything one CSV row needs that is not a verdict.
struct Row {
    /// Which arm produced it.
    arm: &'static str,
    /// Which field.
    field: &'static str,
    /// `full`, `inscribed`, or `own` for the calibration fields.
    region: &'static str,
    /// Samples per axis.
    resolution: usize,
    /// Draws behind `band_sigma`; `1` where there is no band.
    realisations: usize,
    /// The field's own spectral moment, by central differences.
    lambda: f64,
    /// χ of the digital solid at the standardised zero level.
    chi: i64,
    /// The formula's prediction from the field's own `λ` and box.
    predicted: f64,
    /// The matched ensemble's standard deviation of χ.
    sigma: f64,
    /// `|chi − predicted| ≤ 2σ`.
    within: bool,
    /// Third and fourth standardised moments, where measured.
    skewness: Option<f64>,
    excess_kurtosis: Option<f64>,
    /// Share of samples the cap decides, on `noise_cavity`.
    cap_share: Option<f64>,
    /// The band, where there is one.
    band: Option<Band>,
    /// Whether the control is scored on this row.
    control_scored: bool,
}

impl Row {
    fn sigma_distance(&self) -> Option<f64> {
        (self.sigma > 0.0).then(|| (self.chi as f64 - self.predicted).abs() / self.sigma)
    }
}

/// `f64` for the CSV: enough digits to reproduce, no separators.
fn num(value: f64) -> String {
    format!("{value:.9}")
}

/// `Option<f64>` for the CSV.
fn opt(value: Option<f64>) -> String {
    value.map_or_else(|| String::from("na"), num)
}

/// One field on one box: sample, standardise, census, band.
///
/// The moments are taken only where C2 reads them — the inscribed
/// `noise_cavity` grid at the finest rung — and the control is scored on every
/// inscribed row, as registered.
fn measure(
    field: Field,
    region: &'static str,
    half: f64,
    n: usize,
    rng: &mut Rng,
    weights: &[i64; 256],
) -> Row {
    let inscribed = region == "inscribed";
    let with_moments = inscribed && field == Field::NoiseCavity && n == SUB_SAMPLES_FINEST;
    let (values, step, sides) = field.sample(n, half);
    assert!(
        values.iter().all(|v| v.is_finite()),
        "VOID: {} sampled a non-finite value on the {region} box at {n}^3, and every comparison \
         below would read as a cleanly empty excursion set rather than as broken",
        field.name()
    );
    let (mean, sd) = moments(&values);
    let standard: Vec<f64> = values.iter().map(|v| (v - mean) / sd).collect();
    let lambda = lambda_central_difference(&standard, n, step);
    let lk = box_lipschitz_killing(sides);
    let band = matched_band(lambda, sides, lk, n, rng, weights);
    let occupied: Vec<bool> = standard.iter().map(|v| *v > 0.0).collect();
    let chi = chi_solid(&occupied, n, Domain::Bounded, weights);
    let predicted = expected_chi(lk, lambda, 0.0);
    let within = (chi as f64 - predicted).abs() <= 2.0 * band.sigma;
    let (skewness, excess_kurtosis) = if with_moments {
        let (s, k) = higher_moments(&standard);
        (Some(s), Some(k))
    } else {
        (None, None)
    };
    let cap_share = (field == Field::NoiseCavity).then(|| {
        let (lo, hi) = ([-half; 3], [half; 3]);
        cap_dominated_share(n, lo, hi)
    });
    Row {
        arm: if inscribed { "sub_box" } else { "full_box" },
        field: field.name(),
        region,
        resolution: n,
        realisations: REALISATIONS,
        lambda,
        chi,
        predicted,
        sigma: band.sigma,
        within,
        skewness,
        excess_kurtosis,
        cap_share,
        band: Some(band),
        control_scored: inscribed,
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-180");
    let weights = weight_numerators();

    common::experiment::run(prereg, |run| {
        let mut rng = Rng::new(SEED);
        let mut rows: Vec<Row> = Vec::new();

        // ── calibration: M-458's own numbers, through this file's copy ─────
        let sphere = Sphere::<f64>::canonical();
        let torus = Torus::<f64>::canonical();
        let mut calibration = Vec::new();
        for (name, values) in [
            ("sphere", {
                let (lo, hi) = sphere.domain();
                sample_field(&sphere, FULL_SAMPLES, lo, hi).0
            }),
            ("torus", {
                let (lo, hi) = torus.domain();
                sample_field(&torus, FULL_SAMPLES, lo, hi).0
            }),
        ] {
            let occupied: Vec<bool> = values.iter().map(|v| -*v > 0.0).collect();
            let solid = chi_solid(&occupied, FULL_SAMPLES, Domain::Bounded, &weights);
            calibration.push((name, 2 * solid));
        }
        let chi_sphere = calibration[0].1;
        let chi_torus = calibration[1].1;
        assert_eq!(
            (chi_sphere, chi_torus),
            (2, 0),
            "VOID: the copied oracle reads surface chi {chi_sphere} on sphere and {chi_torus} on \
             torus where M-458 reads 2 and 0 -- the copy has drifted from the instrument this row \
             claims to be using"
        );

        // ── full box: M-490's rows, re-taken ───────────────────────────────
        for field in [Field::NoiseCavity, Field::FbmTerrain] {
            rows.push(measure(
                field,
                "full",
                field.full_half(),
                FULL_SAMPLES,
                &mut rng,
                &weights,
            ));
        }
        let full_chi: Vec<i64> = rows.iter().map(|r| r.chi).collect();
        assert_eq!(
            full_chi,
            vec![15, 1],
            "VOID: on the full box this copy reads chi {full_chi:?} for (noise_cavity, \
             fbm_terrain) where M-490 reads (15, 1) -- the instrument has drifted from the one \
             whose miss this row exists to explain, and nothing below may be read"
        );

        // ── sub-box: the cap-free inscribed cube, h-matched and finest ─────
        for field in [Field::NoiseCavity, Field::FbmTerrain] {
            for n in [SUB_SAMPLES_H_MATCHED, SUB_SAMPLES_FINEST] {
                rows.push(measure(
                    field,
                    "inscribed",
                    field.sub_half(),
                    n,
                    &mut rng,
                    &weights,
                ));
            }
        }

        // ── verdicts, exactly as registered ────────────────────────────────
        let find = |field: Field, region: &str, n: usize| {
            rows.iter()
                .find(|r| r.field == field.name() && r.region == region && r.resolution == n)
                .unwrap_or_else(|| panic!("VOID: no row for {} {region} {n}", field.name()))
        };
        let cavity_finest = find(Field::NoiseCavity, "inscribed", SUB_SAMPLES_FINEST);
        let c1 = cavity_finest.within;
        let c2 = match (cavity_finest.skewness, cavity_finest.excess_kurtosis) {
            (Some(s), Some(k)) => s.abs() <= 0.1 && k.abs() <= 0.2,
            _ => panic!("VOID: the finest inscribed noise_cavity row carries no moments"),
        };
        let terrain_finest = find(Field::FbmTerrain, "inscribed", SUB_SAMPLES_FINEST);
        let cavity_full = find(Field::NoiseCavity, "full", FULL_SAMPLES);
        let c3 = terrain_finest.within && !cavity_full.within;
        let rejects_wrong = cavity_finest
            .band
            .as_ref()
            .is_some_and(Band::rejects_wrong_lambda);

        // ── the calibration rows, written last so the verdicts are known ───
        for (name, chi) in calibration {
            rows.push(Row {
                arm: "calibration",
                field: name,
                region: "own",
                resolution: FULL_SAMPLES,
                realisations: 1,
                lambda: f64::NAN,
                chi,
                predicted: f64::NAN,
                sigma: 0.0,
                within: false,
                skewness: None,
                excess_kurtosis: None,
                cap_share: None,
                band: None,
                control_scored: false,
            });
        }

        println!(
            "\nC1 sub-box noise_cavity within band {c1}   C2 |skew| <= 0.1 and |kurt| <= 0.2 \
             {c2}   C3 fbm_terrain inscribed within and full-box noise_cavity outside {c3}\n\
             VACUITY: 4-lambda rejected at u = 0 on the finest inscribed noise_cavity band \
             {rejects_wrong}"
        );
        for row in &rows {
            println!(
                "  {:<13} {:<9} {:>3}  chi {:>4}  predicted {:>14}  sigma_distance {}",
                row.field,
                row.region,
                row.resolution,
                row.chi,
                num(row.predicted),
                opt(row.sigma_distance())
            );
        }

        for row in &rows {
            let band = row.band.as_ref();
            run.record(&[
                ("arm", String::from(row.arm)),
                ("field", String::from(row.field)),
                ("region", String::from(row.region)),
                ("resolution", row.resolution.to_string()),
                ("isovalue", num(0.0)),
                ("realisations", row.realisations.to_string()),
                ("lambda_spectral", num(row.lambda)),
                ("chi_digital", row.chi.to_string()),
                ("chi_gkf_predicted", num(row.predicted)),
                ("band_sigma", num(row.sigma)),
                ("gkf_within_band", row.within.to_string()),
                ("sigma_distance", opt(row.sigma_distance())),
                ("skewness", opt(row.skewness)),
                ("excess_kurtosis", opt(row.excess_kurtosis)),
                ("cap_dominated_share", opt(row.cap_share)),
                (
                    "band_rejects_wrong_lambda",
                    match band {
                        Some(b) if row.control_scored => b.rejects_wrong_lambda().to_string(),
                        _ => String::from("na"),
                    },
                ),
                (
                    "matched_ksq",
                    band.map_or_else(|| String::from("na"), |b| b.ksq.to_string()),
                ),
                ("matched_lambda", opt(band.map(|b| b.lambda))),
                ("matched_chi_mean", opt(band.map(|b| b.mean))),
                ("matched_chi_predicted", opt(band.map(|b| b.predicted))),
                (
                    "wrong_lambda_predicted",
                    opt(band.filter(|_| row.control_scored).map(|b| b.wrong)),
                ),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
            ]);
        }

        // The registered VACUITY CONTROL is scored, not enforced by a panic: an
        // abort here would delete the CSV that is the evidence for WHY it could
        // not be satisfied (✗126's rule).
        if !rejects_wrong {
            println!(
                "VACUOUS as registered: the finest inscribed noise_cavity band does not reject \
                 4*lambda at u = 0, so the comparison cannot discriminate and C1 is unreadable \
                 whatever it says."
            );
        }
    });
}
