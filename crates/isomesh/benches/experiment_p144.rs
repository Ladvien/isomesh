//! **P-144 — whether a periodic value-noise terrain admits a `chi` oracle at all.**
//!
//! Ticket: R-144. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p144
//! ```
//!
//! Writes `docs/experiments/p-144.csv`.
//!
//! # What was missing
//!
//! `P-142` closes the gyroid: with a periodic-conforming grid and the seam
//! identification in `common::tpms`, `chi` is `-8*N^3` and the prediction is a
//! theorem about the space group rather than a fitted constant. That leaves the
//! two *noise* reference fields with nothing at all. `fbm_terrain` declares
//! `closed_in_domain() == false` — *"a heightfield exits through the sides"*
//! (`fields/mod.rs:1384`) — and `expected_euler() == None` — *"not closed, so
//! there is nothing to assert"* (`fields/mod.rs:1387`). The consequence is that
//! every topology gate in the crate simply steps over it: the manifold-dual
//! contouring suite guards its check with `if let Some(chi) = field.expected_euler()`
//! (`manifold_dual_contouring/tests.rs:256`), and the marching-cubes suite says in
//! as many words that the number *"is recorded rather than asserted — exactly what
//! `expected_euler() == None` means"* (`marching_cubes/tests.rs:348`).
//!
//! So three of the eight reference fields — `gyroid`, `fbm_terrain`,
//! `noise_cavity` — have never had their topology gated by anything. `P-142`
//! takes one of them. This row asks whether the remaining kind can be taken at
//! all, and it is registered as an **empirical** question because there is no
//! formula to check: a lattice-noise field with a fixed seed is a deterministic
//! function with a definite `chi`, and either the extraction resolves it or it
//! does not.
//!
//! # Why the field here is volumetric, and why that is part of the answer
//!
//! The shipped `fbm_terrain` is `f(p) = p.y - (base + amplitude * fbm(p.x, p.z))`
//! (`fields/mod.rs:1280`). **It cannot be wrapped in `y` by construction** — it is
//! monotone in `y`, so the sign at `y = lo` and the sign at `y = hi` are opposite
//! everywhere and no identification of the top and bottom faces exists. A
//! heightfield therefore has no `chi` on the 3-torus, ever, at any resolution;
//! that is not a measurement, it is arithmetic, and it is recorded here as the
//! first half of C1's answer.
//!
//! What *can* have a `chi` is the volumetric sibling: the same lattice value
//! noise thresholded as a solid, which is the shape `noise_cavity`
//! (`fields/mod.rs:1226`) and every cave generator in the wild actually meshes.
//! Made exactly periodic — the corner hash taken **modulo the lattice period**,
//! so `F(p + P*e_k) == F(p)` to the last bit the arithmetic allows — it is a
//! closed surface in `T^3` and the wrap of `common::tpms` closes it. That is the
//! field this row measures, and it is the most favourable case available: if the
//! oracle is unreachable here it is unreachable for the terrain too.
//!
//! # The field, stated so it can be reproduced without this file
//!
//! One period is `PERIOD = 4` world units per axis and the extraction box is
//! exactly `[0, 4]^3`. Octave `o` samples an integer lattice of
//! `4 * 2^o` cells per axis; corner values are `2*u - 1` for the top 53 bits of a
//! Murmur-style finalizer mix of `(seed ^ octave salt, cell mod 4*2^o)`,
//! interpolated by the smoothstep-weighted trilinear blend. Amplitudes halve per
//! octave and the sum is normalised, so the field lands in `[-1, 1]`. The
//! isovalue is `0.0` and `sample(p) = 0.0 - fbm(p)`, negative inside the solid,
//! which is the crate's sign convention (`lib.rs:56`).
//!
//! The modulus is the whole mechanism. Without it a lattice-noise field has no
//! period at all and the seam identification is meaningless; with it the far face
//! of the box samples the same lattice corners as the near face, which is the
//! condition `NodalTpms::periodic_grid` exists to guarantee for the TPMS and
//! which the value-noise field satisfies at **every** resolution, because the box
//! is one whole period on the nose.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `wrap_mode = periodic` | opposite faces identified by the period translation, via `common::tpms::wrap_seams` | no |
//! | `wrap_mode = none` | the identical extraction, left open in the box | **yes** |
//! | `octaves = 1` | one lattice, `4^3` cells over the box | no |
//! | `octaves = 3` | three lattices, finest `16^3` cells over the box | no |
//! | `seed` (five of them) | the hash seed, with everything else held fixed | no |
//! | resolution `17 / 33 / 49 / 65 / 81 / 97 / 129` samples per axis | how finely the same field is sampled | no |
//!
//! Ten fields (`2 octaves x 5 seeds`), seven rungs each, two arms per rung:
//! **140 rows** off **70** extractions, because the two `wrap_mode` arms read the
//! same extracted buffer and differ only in whether the seam is identified before
//! `chi` is counted. Sharing the extraction is not a shortcut — it is the only
//! way the control arm is a control rather than a second experiment.
//!
//! The ladder carries **three** rungs above the `65^3` C1 names, and that is the
//! point of its shape. One rung above would let a `chi` that happened to repeat
//! once read as converged — measured, exactly that happens: at three octaves one
//! seed's `chi` is `-32` at both `65` and `81` and then moves again at `97`. At
//! `129^3` the grid holds `128^3` cells against `64^3`, **eight times** as many,
//! so a `chi` unchanged across all four upper rungs is a `chi` that has genuinely
//! stopped moving. The whole sweep is a few tens of seconds, well inside the
//! phase's two-minute budget, so there is no reason to buy the weaker ladder.
//!
//! # The three ladder columns, and why none of them can be trivially satisfied
//!
//! - **`chi_stable`** is true for a rung when `chi` there equals `chi` at every
//!   **higher** rung, and there is at least one higher rung. The top rung is
//!   therefore `false` by construction; a column that read `true` for a
//!   single-point ladder would be measuring nothing.
//! - **`resolution_convergence`** is the lowest rung whose `chi_stable` is true,
//!   as a bare integer, or the token `never`. It is a property of the whole
//!   (seed, octaves, wrap) ladder and repeats on each of that ladder's rows.
//! - **`oracle_exists`** additionally demands the surface actually be closed:
//!   `wrap_mode = periodic`, converged at or below `65`, and `boundary_edges == 0`
//!   on every rung of the ladder. It is `false` on every control row by
//!   definition — the `chi` of a surface with a boundary is a number about the
//!   box, not about the field.
//!
//! `chi_variance_across_seeds` is the population variance of `chi` over the five
//! seeds sharing this row's `(octaves, wrap_mode, resolution)`, so it is
//! identical on those five rows and is the quantity C2 is decided on.
//!
//! `c1_holds` and `c2_holds` are **global** verdicts and carry the same value on
//! every row: C1 is a claim about every seed at once ("converges and is stable to
//! at least `65^3`"), and C2 is a claim about the spread across seeds, which no
//! single row can hold. The per-row evidence is in `chi_stable`,
//! `resolution_convergence` and `chi_variance_across_seeds`.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE line is **`none`**, and it is right to be. A null
//! here moves nothing: an oracle that needs more than `65^3` to stabilise is an
//! oracle no gate in this repository can afford, since the golden fixture's
//! largest resolution is far below that and a validity gate has to run on every
//! push. A positive result would not move a stage either — it would open a
//! *fixture* ticket (a ninth reference field with a known `chi`), which is
//! Phase 28's business and is priced there, not here.
//!
//! # Vacuity controls
//!
//! - **Three seeds at least, and distinct.** The registered control, checked on
//!   the de-duplicated seed set rather than on the array length, so a copy-paste
//!   that repeated a seed is caught. Without it `chi_variance_across_seeds` is a
//!   zero that could not have been non-zero (`M-44`) and C2 cannot fire. The seed
//!   travels on every row as the `seed` extra.
//! - **The field really is periodic.** `max |F(p + P*e_k) - F(p)|` over 1331
//!   off-lattice probes and all three axes, for every seed and octave count, must
//!   be under `1e-12`. Recorded as `max_periodicity_residual`. If the field were
//!   not periodic the wrap would be identifying faces that do not match and every
//!   `chi` in the file would be an artefact.
//! - **The seed is wired in.** The smallest over seed pairs of
//!   `max_p |F_a(p) - F_b(p)|` must be strictly positive. Recorded as
//!   `min_seed_separation`. A seed that did not reach the hash would give a
//!   seed-independent `chi` and *falsify C2 for the wrong reason*.
//! - **`octaves` is not a decorative column.** The two octave counts must produce
//!   different fields on the same probes.
//! - **The two `wrap_mode` arms are different measurements.** Every control row
//!   must have `boundary_edges > 0` (the open extraction is genuinely open) and
//!   every wrapped row must have identified at least one seam pair. `common::tpms`
//!   measured that the non-wrapped arm is recognised by its boundary edges and
//!   **not** by its `chi`, which agreed with the prediction on one field in three
//!   by coincidence.
//! - **Every extraction produced a surface.** The smallest triangle count over
//!   all 50 extractions must be positive, or a `chi` of zero would be the `chi` of
//!   the empty mesh.
//! - **The ladder can falsify C1.** `65` must be on the ladder and at least one
//!   rung must sit above it, or "stable to at least `65^3`" is unfalsifiable.

// The periodicity control asserts a residual against an absolute tolerance
// rather than against zero, so no exact float comparison is needed anywhere.

mod common;

use std::collections::BTreeSet;
use std::time::Instant;

use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

use crate::common::tpms::{EulerCount, euler, wrap_seams};

/// Lattice cells per axis at octave 0, which is also the field's period in world
/// units and the side of the extraction box.
///
/// Four rather than one so that a single period already carries a few dozen
/// lattice cells of topology; a period of one would make `chi` a number about a
/// single blob and the convergence question would be about nothing.
const PERIOD: u32 = 4;

/// The isovalue. The solid is `fbm > 0`, so a little under half the volume.
///
/// Zero rather than a tuned offset: the corner values are 53-bit fractions, the
/// chance any of them lands exactly on the isovalue is `2^-53` apiece, and
/// `non_manifold_edges` is recorded so a hit would be visible rather than
/// silently costing one from `chi` (the mechanism `common::tpms` measured on
/// Schwarz D).
const ISO: f64 = 0.0;

/// Octave counts swept. One lattice, and three.
///
/// With `PERIOD = 4` the finest lattice is `4^3` at one octave and `16^3` at
/// three, so the two arms differ by a factor of four in feature size — which is
/// the whole variable C1 is sensitive to.
const OCTAVES: [u32; 2] = [1, 3];

/// The seeds. Five, where the registration's vacuity control demands three.
const SEEDS: [u64; 5] = [
    0x5EED_0000_0000_0001,
    0x5EED_0000_0000_0002,
    0x0BAD_C0FF_EE0D_D00D,
    0x1234_5678_9ABC_DEF0,
    0xF00D_FACE_CAFE_B0BA,
];

/// Samples per axis. `n` samples span `n - 1` cells (`shape.rs:11`), and the box
/// is exactly one period, so every rung is periodic-conforming automatically.
///
/// Three rungs sit above [`C1_RUNG`], so "stable to at least `65^3`" is a claim
/// about `65`, `81`, `97` and `129` agreeing rather than about one lucky repeat.
const LADDER: [u32; 7] = [17, 33, 49, 65, 81, 97, 129];

/// The rung C1 names. There must be a rung above it or C1 cannot be falsified.
const C1_RUNG: u32 = 65;

/// Probes per axis for the vacuity controls. `11^3 = 1331` points.
const PROBE: u32 = 11;

/// Per-axis probe offsets, in units of one probe step.
///
/// The same trick `common::tpms::shift_residuals` uses: `1/3`, `1/5` and `1/7`
/// keep every probe strictly inside a lattice cell, so the periodicity control
/// is evaluated where interpolation happens rather than only on the corners,
/// where it would be an identity on the hash alone.
const PROBE_OFFSETS: [f64; 3] = [1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0];

/// Largest tolerated `|F(p + P*e_k) - F(p)|`.
///
/// The corner hashes at `p` and `p + P` are the *same* `u64` by construction, so
/// the only difference is the rounding of `p + P` itself: one ulp of a number
/// under 8, amplified by the interpolation weights. That is `1e-15` territory,
/// and a genuinely aperiodic field would be `O(1)`.
const PERIODICITY_TOLERANCE: f64 = 1e-12;

/// `2^53`, the denominator that turns a hash's top 53 bits into a fraction.
const TWO_POW_53: f64 = 9_007_199_254_740_992.0;

/// Per-octave seed salt, the golden-ratio odd constant.
const OCTAVE_SALT: u64 = 0x9E37_79B9_7F4A_7C15;

/// Per-axis additive salts, so that `(1, 0, 0)` and `(0, 1, 0)` cannot collide.
const AXIS_SALT: [u64; 3] = [
    0xA076_1D64_78BD_642F,
    0xE703_7ED1_A0B4_28DB,
    0x8EBC_6AF0_9C88_C6E3,
];

/// Column token for the wrapped arm.
const WRAP_PERIODIC: &str = "periodic";

/// Column token for the open control arm.
const WRAP_NONE: &str = "none";

/// Index of the control arm within a rung's pair.
const ARM_NONE: usize = 0;

/// Index of the wrapped arm within a rung's pair.
const ARM_PERIODIC: usize = 1;

/// The finalizer mix. Deterministic, seeded, and written here rather than taken
/// from `rand`, which is not a dependency and would not be reproducible across
/// versions anyway.
///
/// Two rounds of the Murmur3 64-bit finalizer per axis, folded into a running
/// state, so every bit of every coordinate reaches every bit of the output.
fn mix(seed: u64, cell: [u64; 3]) -> u64 {
    let mut h = seed;
    for (coordinate, salt) in cell.iter().zip(AXIS_SALT.iter()) {
        h ^= coordinate.wrapping_add(*salt);
        h = h.wrapping_mul(0xFF51_AFD7_ED55_8CCD);
        h ^= h >> 33;
        h = h.wrapping_mul(0xC4CE_B9FE_1A85_EC53);
        h ^= h >> 29;
    }
    h
}

/// The value at one lattice corner, in `[-1, 1)`.
fn corner_value(seed: u64, cell: [u64; 3]) -> f64 {
    let fraction = (mix(seed, cell) >> 11) as f64 / TWO_POW_53;
    2.0 * fraction - 1.0
}

/// The smoothstep weight. `C1` at the cell boundaries, which keeps the isosurface
/// from carrying the lattice's own creases as spurious topology.
fn smoothstep(t: f64) -> f64 {
    t * t * (3.0 - 2.0 * t)
}

/// One octave of value noise on a lattice of `modulus` cells per axis, sampled
/// at `q` in lattice units.
///
/// The corner index is reduced `mod modulus` **before** it is hashed, which is
/// the entire periodicity mechanism: the corner at `modulus` is the corner at
/// `0`, so `F` repeats exactly every `modulus` lattice cells.
fn lattice_value(q: [f64; 3], modulus: u32, seed: u64) -> f64 {
    let base = [q[0].floor(), q[1].floor(), q[2].floor()];
    let m = i64::from(modulus);
    let lo = [
        (base[0] as i64).rem_euclid(m) as u64,
        (base[1] as i64).rem_euclid(m) as u64,
        (base[2] as i64).rem_euclid(m) as u64,
    ];
    let modulus = u64::from(modulus);
    let hi = [
        (lo[0] + 1) % modulus,
        (lo[1] + 1) % modulus,
        (lo[2] + 1) % modulus,
    ];
    let t = [
        smoothstep(q[0] - base[0]),
        smoothstep(q[1] - base[1]),
        smoothstep(q[2] - base[2]),
    ];

    let mut acc = 0.0;
    for corner in 0..8u32 {
        let mut cell = [0u64; 3];
        let mut weight = 1.0;
        for axis in 0..3 {
            if corner >> axis & 1 == 0 {
                cell[axis] = lo[axis];
                weight *= 1.0 - t[axis];
            } else {
                cell[axis] = hi[axis];
                weight *= t[axis];
            }
        }
        acc += weight * corner_value(seed, cell);
    }
    acc
}

/// Exactly periodic fractal value noise, thresholded into a solid.
#[derive(Clone, Copy, Debug)]
struct PeriodicValueNoise {
    /// Hash seed. Everything else about the field is a constant of this file.
    seed: u64,
    /// Octaves summed, each at twice the previous frequency.
    octaves: u32,
}

impl PeriodicValueNoise {
    /// The normalised fractal sum, in `[-1, 1]`.
    fn fbm(&self, p: [f64; 3]) -> f64 {
        let mut total = 0.0;
        let mut norm = 0.0;
        let mut amplitude = 1.0;
        for octave in 0..self.octaves {
            let frequency = f64::from(1u32 << octave);
            let modulus = PERIOD << octave;
            let q = [p[0] * frequency, p[1] * frequency, p[2] * frequency];
            let seed = self.seed ^ OCTAVE_SALT.wrapping_mul(u64::from(octave) + 1);
            total += amplitude * lattice_value(q, modulus, seed);
            norm += amplitude;
            amplitude *= 0.5;
        }
        total / norm
    }

    /// Cells per axis of this field's finest lattice.
    fn finest_lattice(&self) -> u32 {
        PERIOD << (self.octaves - 1)
    }
}

impl Sdf for PeriodicValueNoise {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        ISO - self.fbm(p)
    }
}

/// One `chi` reading: a `(seed, octaves, resolution, wrap_mode)` cell of the
/// sweep.
#[derive(Clone, Copy, Debug)]
struct Arm {
    /// What `common::tpms::euler` counted.
    count: EulerCount,
    /// Seam vertex pairs identified. Zero on the control arm by construction.
    identified: u64,
    /// Non-degenerate triangles the arm's own buffer carries.
    triangles: usize,
    /// Wall clock of the extraction the two arms share.
    extract_ms: f64,
}

/// Extract once at `samples` per axis, then read `chi` twice: open, and wrapped.
///
/// Returns `[control, wrapped]`, indexed by [`ARM_NONE`] and [`ARM_PERIODIC`].
fn measure(field: &PeriodicValueNoise, samples: u32, mc: &mut MarchingCubes<f64>) -> [Arm; 2] {
    let cell_size = f64::from(PERIOD) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("the value-noise grid fits u32");
    let mut mesh = MeshBuffer::<f64>::new();

    let started = Instant::now();
    mc.extract(field, &shape, [0.0; 3], cell_size, &mut mesh)
        .expect("marching cubes over a finite field on a grid of at least two samples");
    let extract_ms = started.elapsed().as_secs_f64() * 1e3;

    let tol = isomesh::weld::epsilon_for(cell_size);
    let open = Arm {
        count: euler(&mesh.positions, &mesh.indices, tol),
        identified: 0,
        triangles: mesh.triangle_count(),
        extract_ms,
    };

    let mut wrapped = mesh;
    let identified = wrap_seams(&mut wrapped, [0.0; 3], [f64::from(PERIOD); 3], tol);
    let closed = Arm {
        count: euler(&wrapped.positions, &wrapped.indices, tol),
        identified,
        triangles: wrapped.triangle_count(),
        extract_ms,
    };

    [open, closed]
}

/// Flat index of one `(octaves, seed, rung, arm)` cell.
const fn cell_index(octave: usize, seed: usize, rung: usize, arm: usize) -> usize {
    ((octave * SEEDS.len() + seed) * LADDER.len() + rung) * 2 + arm
}

/// The 1331 off-lattice probe points the vacuity controls are evaluated on.
fn probe_points() -> Vec<[f64; 3]> {
    let span = f64::from(PERIOD);
    let steps = f64::from(PROBE);
    let mut points = Vec::new();
    for iz in 0..PROBE {
        for iy in 0..PROBE {
            for ix in 0..PROBE {
                points.push([
                    span * (f64::from(ix) + PROBE_OFFSETS[0]) / steps,
                    span * (f64::from(iy) + PROBE_OFFSETS[1]) / steps,
                    span * (f64::from(iz) + PROBE_OFFSETS[2]) / steps,
                ]);
            }
        }
    }
    points
}

/// `max |F(p + P*e_k) - F(p)|` over the probes and all three axes.
fn periodicity_residual(field: &PeriodicValueNoise, probes: &[[f64; 3]]) -> f64 {
    let span = f64::from(PERIOD);
    let mut worst = 0.0f64;
    for p in probes {
        let here = field.sample(*p);
        for axis in 0..3 {
            let mut shifted = *p;
            shifted[axis] += span;
            worst = worst.max((field.sample(shifted) - here).abs());
        }
    }
    worst
}

/// `max |A(p) - B(p)|` over the probes.
fn separation(a: &PeriodicValueNoise, b: &PeriodicValueNoise, probes: &[[f64; 3]]) -> f64 {
    probes
        .iter()
        .map(|p| (a.sample(*p) - b.sample(*p)).abs())
        .fold(0.0f64, f64::max)
}

/// Population variance of a ladder of integers.
fn variance(values: &[i64]) -> f64 {
    let n = values.len() as f64;
    let mean = values.iter().map(|v| *v as f64).sum::<f64>() / n;
    values
        .iter()
        .map(|v| {
            let d = *v as f64 - mean;
            d * d
        })
        .sum::<f64>()
        / n
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-144");

    common::experiment::run(prereg, |run| {
        let probes = probe_points();

        // ── vacuity control 1: three distinct seeds, or C2 cannot fire ───────
        let distinct: BTreeSet<u64> = SEEDS.iter().copied().collect();
        assert!(
            distinct.len() >= 3,
            "VOID: {} distinct seeds, and C2 needs at least three or \
             chi_variance_across_seeds is a zero that could not have been non-zero (M-44)",
            distinct.len()
        );

        // ── vacuity control 2: the ladder can falsify C1 ─────────────────────
        let above_c1 = LADDER.iter().filter(|rung| **rung > C1_RUNG).count();
        assert!(
            LADDER.contains(&C1_RUNG) && above_c1 >= 1,
            "VOID: the ladder must contain {C1_RUNG} and at least one rung above it, \
             or 'stable to at least 65^3' is unfalsifiable; it has {above_c1} above"
        );

        let fields: Vec<PeriodicValueNoise> = OCTAVES
            .iter()
            .flat_map(|octaves| {
                SEEDS.iter().map(move |seed| PeriodicValueNoise {
                    seed: *seed,
                    octaves: *octaves,
                })
            })
            .collect();

        // ── vacuity control 3: the field is exactly periodic on the box ──────
        let mut max_periodicity_residual = 0.0f64;
        for field in &fields {
            let residual = periodicity_residual(field, &probes);
            assert!(
                residual < PERIODICITY_TOLERANCE,
                "VOID: seed {:#018x} at {} octaves is not periodic on the extraction box \
                 (residual {residual:.3e} >= {PERIODICITY_TOLERANCE:.0e}), so wrap_seams would \
                 identify faces that do not match and every chi in this file would be an artefact",
                field.seed,
                field.octaves
            );
            max_periodicity_residual = max_periodicity_residual.max(residual);
        }

        // ── vacuity control 4: the seed reaches the hash ─────────────────────
        let mut min_seed_separation = f64::INFINITY;
        for octaves in OCTAVES {
            for (i, left) in SEEDS.iter().enumerate() {
                for right in SEEDS.iter().skip(i + 1) {
                    let a = PeriodicValueNoise {
                        seed: *left,
                        octaves,
                    };
                    let b = PeriodicValueNoise {
                        seed: *right,
                        octaves,
                    };
                    min_seed_separation = min_seed_separation.min(separation(&a, &b, &probes));
                }
            }
        }
        assert!(
            min_seed_separation > 0.0,
            "VOID: two seeds produce the same field on all {} probes, so the seed does not \
             reach the hash and a seed-independent chi would falsify C2 for the wrong reason",
            probes.len()
        );

        // ── vacuity control 5: `octaves` is not a decorative column ──────────
        let octave_separation = separation(
            &PeriodicValueNoise {
                seed: SEEDS[0],
                octaves: OCTAVES[0],
            },
            &PeriodicValueNoise {
                seed: SEEDS[0],
                octaves: OCTAVES[1],
            },
            &probes,
        );
        assert!(
            octave_separation > 0.0,
            "VOID: {} and {} octaves give the same field, so the octaves column varies nothing",
            OCTAVES[0],
            OCTAVES[1]
        );

        // ── the sweep: one extraction per (octaves, seed, rung), two arms ────
        let mut mc = MarchingCubes::<f64>::new();
        let mut cells: Vec<Arm> =
            Vec::with_capacity(OCTAVES.len() * SEEDS.len() * LADDER.len() * 2);
        for field in &fields {
            for samples in LADDER {
                let pair = measure(field, samples, &mut mc);
                cells.push(pair[ARM_NONE]);
                cells.push(pair[ARM_PERIODIC]);
            }
        }

        // ── vacuity control 6: every extraction produced a surface ───────────
        let min_triangles = cells
            .iter()
            .map(|arm| arm.triangles)
            .min()
            .expect("the sweep ran at least one configuration");
        assert!(
            min_triangles > 0,
            "VOID: some extraction produced an empty mesh, and the chi of an empty mesh is a \
             zero that could not have been non-zero (M-44)"
        );

        // ── vacuity control 7: the two wrap_mode arms are different readings ─
        for octave in 0..OCTAVES.len() {
            for seed in 0..SEEDS.len() {
                for rung in 0..LADDER.len() {
                    let open = cells[cell_index(octave, seed, rung, ARM_NONE)];
                    let closed = cells[cell_index(octave, seed, rung, ARM_PERIODIC)];
                    assert!(
                        open.count.boundary_edges > 0,
                        "VOID: the control arm at {} octaves, seed {:#018x}, {}^3 has no boundary \
                         edge, so the open extraction is not open and the two wrap_mode arms are \
                         the same measurement",
                        OCTAVES[octave],
                        SEEDS[seed],
                        LADDER[rung]
                    );
                    assert!(
                        closed.identified > 0,
                        "VOID: wrap_seams identified nothing at {} octaves, seed {:#018x}, {}^3, \
                         so the periodic arm did not wrap anything",
                        OCTAVES[octave],
                        SEEDS[seed],
                        LADDER[rung]
                    );
                }
            }
        }

        // ── the ladder reading, per (octaves, seed, arm) ─────────────────────
        let chi_of = |octave: usize, seed: usize, rung: usize, arm: usize| -> i64 {
            cells[cell_index(octave, seed, rung, arm)].count.chi
        };
        // `stable[rung]`: chi here equals chi at every HIGHER rung, and there is
        // one. The top rung is false by construction.
        let stability = |octave: usize, seed: usize, arm: usize| -> Vec<bool> {
            let top = chi_of(octave, seed, LADDER.len() - 1, arm);
            let mut flags = vec![false; LADDER.len()];
            let mut holds = true;
            for rung in (0..LADDER.len()).rev() {
                holds = holds && chi_of(octave, seed, rung, arm) == top;
                flags[rung] = holds && rung + 1 < LADDER.len();
            }
            flags
        };
        let converged = |octave: usize, seed: usize, arm: usize| -> Option<u32> {
            stability(octave, seed, arm)
                .iter()
                .position(|stable| *stable)
                .map(|rung| LADDER[rung])
        };

        // C1 is global: every seed, every octave count, on the wrapped arm.
        let mut c1 = true;
        for octave in 0..OCTAVES.len() {
            for seed in 0..SEEDS.len() {
                c1 =
                    c1 && converged(octave, seed, ARM_PERIODIC).is_some_and(|rung| rung <= C1_RUNG);
            }
        }

        // C2 is global: at the finest rung, on the wrapped arm, the seeds must
        // disagree — for both octave counts, since "seed-independent chi" is a
        // claim about the field family and not about one member of it.
        let top_rung = LADDER.len() - 1;
        let mut c2 = true;
        for octave in 0..OCTAVES.len() {
            let chis: Vec<i64> = (0..SEEDS.len())
                .map(|seed| chi_of(octave, seed, top_rung, ARM_PERIODIC))
                .collect();
            c2 = c2 && chis.iter().any(|chi| *chi != chis[0]);
        }

        // ── the rows ─────────────────────────────────────────────────────────
        for octave in 0..OCTAVES.len() {
            let octaves = OCTAVES[octave];
            let finest = PeriodicValueNoise {
                seed: SEEDS[0],
                octaves,
            }
            .finest_lattice();

            for arm in [ARM_NONE, ARM_PERIODIC] {
                let wrap_mode = if arm == ARM_PERIODIC {
                    WRAP_PERIODIC
                } else {
                    WRAP_NONE
                };

                for seed in 0..SEEDS.len() {
                    let flags = stability(octave, seed, arm);
                    let convergence = converged(octave, seed, arm);
                    let all_closed = (0..LADDER.len()).all(|rung| {
                        cells[cell_index(octave, seed, rung, arm)]
                            .count
                            .boundary_edges
                            == 0
                    });
                    let oracle = arm == ARM_PERIODIC
                        && all_closed
                        && convergence.is_some_and(|rung| rung <= C1_RUNG);
                    let chi_top = chi_of(octave, seed, top_rung, arm);

                    for rung in 0..LADDER.len() {
                        let samples = LADDER[rung];
                        let cell = cells[cell_index(octave, seed, rung, arm)];
                        let across: Vec<i64> = (0..SEEDS.len())
                            .map(|other| chi_of(octave, other, rung, arm))
                            .collect();
                        let span =
                            across.iter().max().unwrap_or(&0) - across.iter().min().unwrap_or(&0);
                        let cells_per_axis = samples - 1;
                        let cell_size = f64::from(PERIOD) / f64::from(cells_per_axis);

                        run.record(&[
                            ("field", "periodic_value_noise".to_string()),
                            ("octaves", octaves.to_string()),
                            ("period", PERIOD.to_string()),
                            ("wrap_mode", wrap_mode.to_string()),
                            ("chi_measured", cell.count.chi.to_string()),
                            (
                                "chi_variance_across_seeds",
                                format!("{:.6}", variance(&across)),
                            ),
                            ("chi_stable", flags[rung].to_string()),
                            ("oracle_exists", oracle.to_string()),
                            (
                                "resolution_convergence",
                                convergence
                                    .map_or_else(|| "never".to_string(), |rung| rung.to_string()),
                            ),
                            ("c1_holds", c1.to_string()),
                            ("c2_holds", c2.to_string()),
                            // ── extras (M-273) ──
                            ("boundary_edges", cell.count.boundary_edges.to_string()),
                            ("cell_size", format!("{cell_size:.9}")),
                            ("cells_per_axis", cells_per_axis.to_string()),
                            ("chi_delta_from_top", (cell.count.chi - chi_top).to_string()),
                            ("chi_span_across_seeds", span.to_string()),
                            ("chi_top_rung", chi_top.to_string()),
                            ("edges", cell.count.edges.to_string()),
                            ("extract_ms", format!("{:.3}", cell.extract_ms)),
                            ("faces", cell.count.faces.to_string()),
                            ("finest_lattice_cells_per_axis", finest.to_string()),
                            ("identified_pairs", cell.identified.to_string()),
                            ("iso", format!("{ISO:.3}")),
                            (
                                "max_periodicity_residual",
                                format!("{max_periodicity_residual:.3e}"),
                            ),
                            ("min_seed_separation", format!("{min_seed_separation:.6}")),
                            (
                                "non_manifold_edges",
                                cell.count.non_manifold_edges.to_string(),
                            ),
                            ("samples_per_axis", samples.to_string()),
                            ("seed", format!("{:#018x}", SEEDS[seed])),
                            ("triangles", cell.triangles.to_string()),
                            ("vertices", cell.count.vertices.to_string()),
                            (
                                "voxels_per_finest_lattice_cell",
                                format!("{:.3}", f64::from(cells_per_axis) / f64::from(finest)),
                            ),
                        ]);
                    }
                }
            }
        }
    });
}
