//! **P-178 — the Gaussian kinematic formula gives an excursion set an expected χ.**
//!
//! Ticket: `R-183`. Pre-registered before this harness existed; minted by the
//! discovery loop (`docs/research/2026-09-12-discovery-loop-prompt.md`,
//! iteration 2) from the question ledger's `Q4`.
//!
//! Writes `docs/experiments/p-178.csv`.
//!
//! # The source, and what is taken from it
//!
//! Taylor, *A Gaussian kinematic formula*, `10.1214/009117905000000594` — in the
//! corpus, converted, and uncited in this repository before this row
//! (`Lipschitz-Killing`, `Hermite polynomial` and `spectral moment` were all
//! zero under `git grep`). Two sentences are used, both quoted in the
//! registration:
//!
//! - the real-valued Gaussian case,
//!   `E[χ(M ∩ f⁻¹[u, ∞))] = Σ_j L_j(M) · ρ_j(u)`, the `L_j` being taken in *"the
//!   metric induced by `f`"*;
//! - the EC densities, `ρ_j(u) = (2π)^(-(j+1)/2) ∫_u^∞ H_j(x) e^(-x²/2) dx`,
//!   with `H_j` the `j`th Hermite polynomial and the spectral parameter `μ₂`
//!   *"just the variance of the first-order partial derivatives of `f`"*.
//!
//! **Nothing is transcribed from memory.** The four `ρ_j` this row needs are
//! computed two ways and asserted equal: by Gauss–Legendre quadrature of the
//! Hermite integral above, and by the closed forms that follow from
//! `d/dx [H_{j-1}(x) e^(-x²/2)] = -H_j(x) e^(-x²/2)`. The box's
//! Lipschitz–Killing curvatures are likewise *fitted* from the Steiner tube
//! polynomial the same paper states — `λ(T(M,r)) = Σ_j M_j(M) r^j / j!` — rather
//! than quoted, so a wrong `L_j` is a failed assert and not a wrong number.
//!
//! Under the induced metric of an isotropic field with spectral moment `λ`,
//! lengths scale by `√λ`, so `L_j(M, g) = λ^(j/2) · L_j(M, euclidean)`. On the
//! flat 3-torus every Lipschitz–Killing curvature except the volume vanishes,
//! which leaves the single term `E[χ] = λ^(3/2) · Vol · (2π)^(-2) (u² - 1)
//! e^(-u²/2)` — the cleanest available form, and the one C1 is scored on.
//!
//! # Why this row exists
//!
//! `M-458` reports `noise_cavity`'s χ as **82** at 33³ and **−152** at 65³ and
//! calls both *ground truth*, and `CLAUDE.md`'s validity gate records χ rather
//! than asserting it on every field whose `expected_euler()` is `None`. The
//! formula above is an analytic expectation for exactly that quantity. Either it
//! grades these fields — and the gate gains an oracle for two of the four fields
//! `M-459` found ungradeable — or it does not, and `relative_departure` is the
//! number saying how far this crate's hash noise sits from a Gaussian field.
//!
//! # The oracle, and why it is copied rather than shared
//!
//! χ is read off the field's own signs with **no mesh involved**, by the
//! Ohser–Nagel–Schladitz block census `M-458` measured: digitise, count the 256
//! `2×2×2` configurations, weight each by `(-1)^dim / shared` summed over the
//! cells it activates, divide by eight. `experiment_p145.rs` owns that oracle and
//! is a landed dataset; this file carries its own copy so that running this row
//! cannot perturb that one. The copy is not trusted — the **calibration arm**
//! re-takes `M-458`'s own numbers, `χ` **2** on `sphere` and **0** on `torus`
//! (surface χ, which is `2 · χ(solid)`), and a drifted copy fails there before
//! any prediction is read.
//!
//! # Arms
//!
//! | arm | what it is | which clause |
//! |---|---|---|
//! | `gaussian_periodic` | an *exactly* Gaussian field on the unit 3-torus, `REALISATIONS` independent draws, λ in closed form | C1 |
//! | `wrong_lambda` | the same measurements against a prediction built from `4λ` | vacuity control |
//! | `matched_gaussian` | a Gaussian ensemble on the shipped field's own box, mode shell chosen so its closed-form λ is nearest the shipped field's measured one, digitised bounded | supplies C2's band |
//! | `shipped` | `noise_cavity` and `fbm_terrain`, standardised, λ from their own sampled derivatives | C2 |
//! | `calibration` | `sphere` and `torus` through the same oracle | anti-drift |
//! | `convergence` | `noise_cavity` at three rungs, one fixed realisation | C3 |
//!
//! # The vacuity controls, and what each would catch
//!
//! - **The band must be able to reject something.** `wrong_lambda` predicts with
//!   `4λ`, which scales the torus term by `4^(3/2) = 8`. If that also lands
//!   inside the band, the comparison cannot discriminate and the row is vacuous
//!   whatever the clauses say. Column: `band_rejects_wrong_lambda`.
//! - **The fixture must really be Gaussian with the λ it claims.** The
//!   closed-form λ is checked against a central-difference estimate on the
//!   sampled grid, and the field's sample variance against 1.
//! - **The oracle must be `M-458`'s.** `calibration_chi_sphere` **2** and
//!   `calibration_chi_torus` **0**, or the copy has drifted.
//! - **The digitisation must be the crate's own inside test.** On the
//!   calibration arm the predicate `-f > 0` is asserted equal to
//!   `marching_cubes::table::is_inside(f)` sample by sample.
//!
//! # Determinism
//!
//! One `Rng` (SplitMix64, `common::poly::Rng`) seeded by a constant, drawn in a
//! fixed order; the shipped fields are deterministic hash noise; every sum is
//! over a fixed index range. Re-running produces byte-identical rows.

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
use isomesh::fields::{FbmTerrain, ReferenceField, Sphere, Torus, noise_cavity};
use isomesh::marching_cubes::table::is_inside;

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

/// Independent draws of the Gaussian fixture.
///
/// The registration says *"at least 32"*; the standard error falls as `1/√R` and
/// 32 is what keeps the whole row inside its `M` budget.
const REALISATIONS: usize = 32;

/// The isovalue ladder, in standard deviations of the field.
///
/// `±1` are the theorem's **own zeros** — `(u² − 1)` vanishes there — so two of
/// the five rungs are a parameter-free prediction of exactly zero that the
/// ensemble has to match without any scale to hide behind.
const RUNGS: [f64; 5] = [-2.0, -1.0, 0.0, 1.0, 2.0];

/// Samples per axis for the Gaussian and shipped arms.
const SAMPLES: usize = 65;

/// The C3 ladder, one fixed realisation of `noise_cavity`.
const CONVERGENCE_RUNGS: [usize; 3] = [33, 65, 129];

/// Largest `|k|²` in the unit-torus mode shell.
///
/// The shell is every integer triple with `1 ≤ |k|² ≤ MAX_K_SQ`, which is closed
/// under the 48-element cube group, so the field is isotropic on the lattice and
/// `Σ k₁² = Σ|k|²/3` exactly rather than approximately.
const MAX_K_SQ: i32 = 14;

/// Seed for every draw in this file.
const SEED: u64 = 0x0178_9d5e_c43a_1f07;

/// Quadrature nodes for the Hermite integral, per half-decade of tail.
const QUAD_NODES: usize = 4096;

/// Where the quadrature truncates the upper tail, in standard deviations.
const QUAD_UPPER: f64 = 12.0;

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

/// `ρ_j(u)` by quadrature of the paper's own integral.
///
/// Composite Simpson on `[u, QUAD_UPPER]`. The integrand decays like
/// `e^(-x²/2)` times a degree-`j` polynomial, so truncating at twelve standard
/// deviations costs less than `1e-30` — smaller than the `1e-12` the closed-form
/// agreement is asserted to.
fn rho_quadrature(j: usize, u: f64) -> f64 {
    let hi = QUAD_UPPER.max(u + 1.0);
    let n = QUAD_NODES * 2;
    let step = (hi - u) / n as f64;
    let mut sum = 0.0;
    for i in 0..=n {
        let x = step.mul_add(i as f64, u);
        let value = hermite(j, x) * (-0.5 * x * x).exp();
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        sum += weight * value;
    }
    let integral = sum * step / 3.0;
    let scale = (2.0 * std::f64::consts::PI).powf(-((j as f64) + 1.0) / 2.0);
    scale * integral
}

/// `ρ_j(u)` in closed form, from `d/dx [H_{j-1} e^(-x²/2)] = -H_j e^(-x²/2)`.
///
/// For `j ≥ 1` the integral telescopes to `H_{j-1}(u) e^(-u²/2)`; for `j = 0` it
/// is the Gaussian tail, written here as the complementary error function.
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

/// Sample an SDF on its own domain at `n³`, returning the values and the spacing.
fn sample_field<S: Sdf<Scalar = f64> + ReferenceField>(
    field: &S,
    n: usize,
) -> (Vec<f64>, f64, [f64; 3]) {
    let (lo, hi) = field.domain();
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

// ─── rows ────────────────────────────────────────────────────────────────────

/// Everything one CSV row needs that is not a verdict.
struct Row {
    /// Which arm produced it.
    arm: &'static str,
    /// Which field, or the fixture's name.
    field: String,
    /// Samples per axis.
    resolution: usize,
    /// The level, in standard deviations.
    isovalue: f64,
    /// How many independent draws the mean is over; `1` for a fixed field.
    realisations: usize,
    /// The spectral moment used for the prediction.
    lambda: f64,
    /// How it was obtained.
    lambda_source: &'static str,
    /// Mean measured χ of the solid.
    chi_mean: f64,
    /// Standard error of that mean; `0` for a single realisation.
    chi_sem: f64,
    /// The formula's prediction.
    predicted: f64,
    /// Whether the mean is inside the two-sigma band.
    within: bool,
    /// `|measured − predicted| / max(1, |predicted|)`.
    departure: f64,
    /// Whether a `4λ` prediction is outside the same band.
    rejects_wrong: String,
    /// χ of one fixed realisation, for the convergence arm.
    fixed: String,
    /// Relative change from the previous rung.
    convergence_gap: String,
}

/// `f64` for the CSV: enough digits to reproduce, no separators.
fn num(value: f64) -> String {
    format!("{value:.9}")
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-178");
    let weights = weight_numerators();

    // ── the ρ_j, derived twice and asserted equal ──────────────────────────
    for j in 0..4 {
        for u in RUNGS {
            let (a, b) = (rho_quadrature(j, u), rho_closed_form(j, u));
            assert!(
                // Absolute floor of 1e-10: the closed form is EXACTLY zero at
                // u = ±1 for j = 3, where H_2(u) = u² − 1 vanishes, so a purely
                // relative tolerance compares quadrature noise against nothing.
                (a - b).abs() <= 1e-10 * a.abs().max(b.abs()).max(1.0),
                "VOID: rho_{j}({u}) is {a} by quadrature of the paper's integral and {b} by the \
                 closed form, so one of the two derivations is wrong and every prediction below \
                 is unreadable"
            );
        }
    }

    common::experiment::run(prereg, |run| {
        let mut rng = Rng::new(SEED);
        let mut rows: Vec<Row> = Vec::new();

        // ── calibration: M-458's own numbers, through this file's copy ─────
        let sphere = Sphere::<f64>::canonical();
        let torus = Torus::<f64>::canonical();
        let mut calibration = Vec::new();
        for (name, values) in [
            ("sphere", sample_field(&sphere, SAMPLES).0),
            ("torus", sample_field(&torus, SAMPLES).0),
        ] {
            let occupied: Vec<bool> = values.iter().map(|v| -*v > 0.0).collect();
            for (index, value) in values.iter().enumerate() {
                assert_eq!(
                    occupied[index],
                    is_inside(*value),
                    "VOID: on {name} the digitisation predicate disagrees with the crate's own \
                     is_inside at sample {index} (value {value}), so this oracle is not reading \
                     the same solid M-458 read"
                );
            }
            let solid = chi_solid(&occupied, SAMPLES, Domain::Bounded, &weights);
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

        // ── C1: the exactly-Gaussian fixture on the unit 3-torus ───────────
        let modes = mode_shell(MAX_K_SQ);
        let torus_lambda = GaussianField::lambda(&modes, 1.0);
        // On a flat 3-torus every Lipschitz-Killing curvature but the volume
        // vanishes, so only the j = 3 term survives.
        let torus_lk = [0.0, 0.0, 0.0, 1.0];
        let mut per_rung: Vec<Vec<f64>> = vec![Vec::new(); RUNGS.len()];
        let mut variance_seen = 0.0;
        let mut lambda_seen = 0.0;
        for _ in 0..REALISATIONS {
            let field = GaussianField::draw(&modes, 1.0, &mut rng);
            let values = field.sample(SAMPLES);
            let (_, sd) = moments(&values);
            variance_seen += sd * sd;
            lambda_seen += lambda_central_difference(&values, SAMPLES, 1.0 / (SAMPLES - 1) as f64);
            for (index, u) in RUNGS.iter().enumerate() {
                let occupied: Vec<bool> = values.iter().map(|v| *v > *u).collect();
                per_rung[index]
                    .push(chi_solid(&occupied, SAMPLES, Domain::Periodic, &weights) as f64);
            }
        }
        variance_seen /= REALISATIONS as f64;
        lambda_seen /= REALISATIONS as f64;
        assert!(
            (variance_seen - 1.0).abs() <= 0.05,
            "VOID: the fixture's sample variance is {variance_seen}, not 1, so it is not the \
             unit-variance field the formula is stated for"
        );
        assert!(
            (lambda_seen - torus_lambda).abs() <= 0.05 * torus_lambda,
            "VOID: the fixture's closed-form lambda is {torus_lambda} and its sampled derivatives \
             say {lambda_seen} -- the two disagree by more than 5%, so the prediction's spectral \
             moment is not the field's"
        );

        let mut rungs_within = 0usize;
        let mut rejects_all = true;
        for (index, u) in RUNGS.iter().enumerate() {
            let draws = &per_rung[index];
            let n = draws.len() as f64;
            let mean = draws.iter().sum::<f64>() / n;
            let variance = draws.iter().map(|c| (c - mean) * (c - mean)).sum::<f64>() / (n - 1.0);
            let sem = (variance / n).sqrt();
            let predicted = expected_chi(torus_lk, torus_lambda, *u);
            let wrong = expected_chi(torus_lk, 4.0 * torus_lambda, *u);
            let within = (mean - predicted).abs() <= 2.0 * sem;
            let rejects = (mean - wrong).abs() > 2.0 * sem;
            if within {
                rungs_within += 1;
            }
            if !rejects {
                rejects_all = false;
            }
            rows.push(Row {
                arm: "gaussian_periodic",
                field: format!("fourier_shell_ksq{MAX_K_SQ}"),
                resolution: SAMPLES,
                isovalue: *u,
                realisations: REALISATIONS,
                lambda: torus_lambda,
                lambda_source: "closed_form",
                chi_mean: mean,
                chi_sem: sem,
                predicted,
                within,
                departure: (mean - predicted).abs() / predicted.abs().max(1.0),
                rejects_wrong: rejects.to_string(),
                fixed: String::from("na"),
                convergence_gap: String::from("na"),
            });
            rows.push(Row {
                arm: "wrong_lambda",
                field: format!("fourier_shell_ksq{MAX_K_SQ}"),
                resolution: SAMPLES,
                isovalue: *u,
                realisations: REALISATIONS,
                lambda: 4.0 * torus_lambda,
                lambda_source: "closed_form_x4",
                chi_mean: mean,
                chi_sem: sem,
                predicted: wrong,
                within: !rejects,
                departure: (mean - wrong).abs() / wrong.abs().max(1.0),
                rejects_wrong: rejects.to_string(),
                fixed: String::from("na"),
                convergence_gap: String::from("na"),
            });
        }
        let c1 = rungs_within >= 4;

        // ── C2: the shipped fields, and a matched Gaussian band ────────────
        let cavity = noise_cavity::<f64>();
        let terrain = FbmTerrain::<f64>::canonical();
        let mut shipped_within = Vec::new();
        for (name, values, sides, step) in [
            {
                let (v, h, s) = sample_field(&cavity, SAMPLES);
                ("noise_cavity", v, s, h)
            },
            {
                let (v, h, s) = sample_field(&terrain, SAMPLES);
                ("fbm_terrain", v, s, h)
            },
        ] {
            assert!(
                values.iter().all(|v| v.is_finite()),
                "VOID: {name} sampled a non-finite value, and every comparison below would read \
                 as a cleanly empty excursion set rather than as broken"
            );
            let (mean, sd) = moments(&values);
            let standard: Vec<f64> = values.iter().map(|v| (v - mean) / sd).collect();
            let lambda = lambda_central_difference(&standard, SAMPLES, step);
            let lk = box_lipschitz_killing(sides);

            // The band: a Gaussian ensemble on this same box, digitised the same
            // bounded way, with the mode shell whose closed-form lambda is
            // nearest this field's. A band taken from any other geometry would
            // be a band for a different question.
            let mut best = (MAX_K_SQ, f64::INFINITY);
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
                let field = GaussianField::draw(&shell, sides[0], &mut rng);
                let sample = field.sample(SAMPLES);
                let occupied: Vec<bool> = sample.iter().map(|v| *v > 0.0).collect();
                matched.push(chi_solid(&occupied, SAMPLES, Domain::Bounded, &weights) as f64);
            }
            let n = matched.len() as f64;
            let matched_mean = matched.iter().sum::<f64>() / n;
            let matched_var = matched
                .iter()
                .map(|c| (c - matched_mean) * (c - matched_mean))
                .sum::<f64>()
                / (n - 1.0);
            let sigma = matched_var.sqrt();
            let matched_predicted = expected_chi(lk, matched_lambda, 0.0);
            rows.push(Row {
                arm: "matched_gaussian",
                field: format!("{name}_matched_ksq{}", best.0),
                resolution: SAMPLES,
                isovalue: 0.0,
                realisations: REALISATIONS,
                lambda: matched_lambda,
                lambda_source: "closed_form_matched",
                chi_mean: matched_mean,
                chi_sem: sigma / n.sqrt(),
                predicted: matched_predicted,
                within: (matched_mean - matched_predicted).abs() <= 2.0 * sigma,
                departure: (matched_mean - matched_predicted).abs()
                    / matched_predicted.abs().max(1.0),
                rejects_wrong: String::from("na"),
                fixed: String::from("na"),
                convergence_gap: String::from("na"),
            });

            // The shipped field itself, one realisation, at its own zero level.
            let occupied: Vec<bool> = standard.iter().map(|v| *v > 0.0).collect();
            let chi = chi_solid(&occupied, SAMPLES, Domain::Bounded, &weights) as f64;
            let predicted = expected_chi(lk, lambda, 0.0);
            let within = (chi - predicted).abs() <= 2.0 * sigma;
            shipped_within.push(within);
            rows.push(Row {
                arm: "shipped",
                field: String::from(name),
                resolution: SAMPLES,
                isovalue: 0.0,
                realisations: 1,
                lambda,
                lambda_source: "central_difference",
                chi_mean: chi,
                chi_sem: sigma,
                predicted,
                within,
                departure: (chi - predicted).abs() / predicted.abs().max(1.0),
                rejects_wrong: String::from("na"),
                fixed: num(chi),
                convergence_gap: String::from("na"),
            });
        }
        let c2 = shipped_within.iter().all(|w| *w);

        // ── C3: one fixed realisation of noise_cavity, three rungs ─────────
        //
        // Two levels, because they are different sets and only one of them is
        // M-458's. `noise_cavity_raw_iso0` is the crate's own inside test,
        // `{f < 0}`, which is exactly the solid M-458 read and whose surface χ
        // it reports as 82 at 33³ and −152 at 65³; `noise_cavity_standardised`
        // is `{g > 0}` for the standardised field, which is the level the
        // Gaussian arm's u = 0 means. C3 is scored on the raw one, because the
        // clause names M-458's numbers.
        let mut previous: Option<f64> = None;
        let mut raw_previous: Option<f64> = None;
        let mut worst_gap = 0.0f64;
        let mut raw_worst_gap = 0.0f64;
        let mut raw_surface: Vec<(usize, i64)> = Vec::new();
        for n in CONVERGENCE_RUNGS {
            let (values, step, sides) = sample_field(&cavity, n);
            let (mean, sd) = moments(&values);
            let standard: Vec<f64> = values.iter().map(|v| (v - mean) / sd).collect();
            let lambda = lambda_central_difference(&standard, n, step);
            let lk = box_lipschitz_killing(sides);
            let occupied: Vec<bool> = standard.iter().map(|v| *v > 0.0).collect();
            let chi = chi_solid(&occupied, n, Domain::Bounded, &weights) as f64;
            let gap = previous.map(|p| (chi - p).abs() / p.abs().max(1.0));
            if let Some(g) = gap {
                worst_gap = worst_gap.max(g);
            }
            previous = Some(chi);
            let predicted = expected_chi(lk, lambda, 0.0);
            rows.push(Row {
                arm: "convergence",
                field: String::from("noise_cavity"),
                resolution: n,
                isovalue: 0.0,
                realisations: 1,
                lambda,
                lambda_source: "central_difference",
                chi_mean: chi,
                chi_sem: 0.0,
                predicted,
                within: false,
                departure: (chi - predicted).abs() / predicted.abs().max(1.0),
                rejects_wrong: String::from("na"),
                fixed: num(chi),
                convergence_gap: gap.map_or_else(|| String::from("na"), num),
            });

            // The same rung, the crate's own inside test, undisplaced by any
            // standardisation: M-458's comparand.
            let raw_occupied: Vec<bool> = values.iter().map(|v| -*v > 0.0).collect();
            let raw_solid = chi_solid(&raw_occupied, n, Domain::Bounded, &weights);
            raw_surface.push((n, 2 * raw_solid));
            let raw_chi = f64::from(raw_solid as i32);
            let raw_gap = raw_previous.map(|p: f64| (raw_chi - p).abs() / p.abs().max(1.0));
            if let Some(g) = raw_gap {
                raw_worst_gap = raw_worst_gap.max(g);
            }
            raw_previous = Some(raw_chi);
            rows.push(Row {
                arm: "convergence",
                field: String::from("noise_cavity_raw_iso0"),
                resolution: n,
                isovalue: 0.0,
                realisations: 1,
                lambda,
                lambda_source: "central_difference",
                chi_mean: f64::from(2 * raw_solid as i32),
                chi_sem: 0.0,
                predicted: f64::NAN,
                within: false,
                departure: f64::NAN,
                rejects_wrong: String::from("na"),
                fixed: (2 * raw_solid).to_string(),
                convergence_gap: raw_gap.map_or_else(|| String::from("na"), num),
            });
        }
        let c3 = raw_worst_gap <= 0.10;

        // ── the calibration rows, written last so the verdicts are known ───
        for (name, chi) in calibration {
            rows.push(Row {
                arm: "calibration",
                field: String::from(name),
                resolution: SAMPLES,
                isovalue: 0.0,
                realisations: 1,
                lambda: f64::NAN,
                lambda_source: "not_applicable",
                chi_mean: f64::from(chi as i32),
                chi_sem: 0.0,
                predicted: f64::NAN,
                within: false,
                departure: f64::NAN,
                rejects_wrong: String::from("na"),
                fixed: chi.to_string(),
                convergence_gap: String::from("na"),
            });
        }

        let rejecting: usize = rows
            .iter()
            .filter(|r| r.arm == "gaussian_periodic" && r.rejects_wrong == "true")
            .count();
        println!(
            "\nC1 rungs within band {rungs_within} of {}   C2 shipped within {:?}   \
             C3 raw worst gap {raw_worst_gap:.6} (standardised {worst_gap:.6})   \
             M-458 comparand {raw_surface:?}\nVACUITY: 4-lambda rejected on \
             {rejecting} of {} rungs, every-rung {rejects_all}",
            RUNGS.len(),
            shipped_within,
            RUNGS.len()
        );

        for row in &rows {
            run.record(&[
                ("arm", String::from(row.arm)),
                ("field", row.field.clone()),
                ("resolution", row.resolution.to_string()),
                ("isovalue", num(row.isovalue)),
                ("realisations", row.realisations.to_string()),
                ("lambda_spectral", num(row.lambda)),
                ("lambda_source", String::from(row.lambda_source)),
                ("chi_digital_mean", num(row.chi_mean)),
                ("chi_digital_sem", num(row.chi_sem)),
                ("chi_gkf_predicted", num(row.predicted)),
                ("gkf_within_band", row.within.to_string()),
                ("relative_departure", num(row.departure)),
                ("band_rejects_wrong_lambda", row.rejects_wrong.clone()),
                ("chi_fixed_realisation", row.fixed.clone()),
                ("chi_convergence_gap", row.convergence_gap.clone()),
                (
                    "rungs_within_band",
                    format!("{rungs_within}of{}", RUNGS.len()),
                ),
                ("calibration_chi_sphere", chi_sphere.to_string()),
                ("calibration_chi_torus", chi_torus.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
            ]);
        }

        // The registered VACUITY CONTROL is scored, not enforced by a panic: an
        // abort here would delete the CSV that is the evidence for WHY it could
        // not be satisfied. `band_rejects_wrong_lambda` carries the verdict per
        // rung and the FINDINGS entry carries the tier.
        if !rejects_all {
            println!(
                "VACUOUS as registered: the band does not reject 4*lambda on every rung. \
                 At u = +/-1 the prediction's (u^2 - 1) factor is zero for ANY lambda, so the \
                 right and wrong predictions are the same number there and no band can separate \
                 them. The control was written without noticing that two of its five rungs are \
                 lambda-blind."
            );
        }
    });
}
