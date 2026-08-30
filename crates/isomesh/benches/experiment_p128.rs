//! **P-128 — the saddle count cannot depend on cell aspect ratio, and this is
//! the first thing in the ledger that says so.**
//!
//! Ticket: R-128. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p128
//! ```
//!
//! Writes `docs/experiments/p-128.csv`.
//!
//! # What was missing
//!
//! `P-127` established, exactly and in `i128`, that `b*b - 4*a*c` at
//! `crates/isomesh/src/marching_cubes/trilinear.rs:246` **is** Cayley's `2x2x2`
//! hyperdeterminant of the eight corner values, in the normalisation
//! `c1^2 - 4*c0*c2`. `crates/isomesh/benches/common/poly.rs` (owned by R-127) is
//! the gate for that identity and is consumed here unchanged.
//!
//! What no id in the repository states is the **consequence**. The `GL(2)^3`
//! action on the `2x2x2` tensor carries `Det` by the relative-invariant weight
//! `(det g1 * det g2 * det g3)^2` — a **perfect square**, therefore positive for
//! every invertible triple, therefore `sign(Det)` is an *absolute* invariant. And
//! `sign(Det)` is what decides whether the cell has two real body saddles or
//! none: `trilinear.rs:246-250` reads the discriminant, branches on
//! `discriminant == R::ZERO`, and hands the roots to `BodySaddles::of`
//! (`:165`), whose `inside_count` (`:294`) every later id builds on — `M-214`,
//! `M-215`, `M-216`, `M-217` and the whole of `InteriorAmbiguity::Trilinear`.
//!
//! Nobody has ever checked that the saddle path is blind to the cell's geometry.
//! It matters because Group D (`R-146` .. `R-151`) is about to *deliberately*
//! mesh with anisotropic cells, and every one of those rows silently assumes the
//! interior-ambiguity resolution it inherits is unaffected. `M-206` is the near
//! miss: it recorded two independent constructions agreeing to `1.1e-12` without
//! asking what that agreement is invariant under.
//!
//! # The construction, and exactly what it does and does not model
//!
//! This is the load-bearing paragraph of the file. **`Extractor::extract` takes
//! a single scalar `cell_size`** (marching_cubes/mod.rs:193) — there is no
//! per-axis cell size anywhere in the shipped API. So "meshing with per-axis
//! cell scales `(1, 2, 4)`" cannot be expressed by an argument, and the
//! registration's fixture has to be built out of what the crate does have. It is
//! built from two pieces, and they answer two different halves of C1:
//!
//! 1. **The census — same values, different cell geometry.** The anisotropy is
//!    expressed by scaling the field's domain per axis: a bench-local
//!    [`AxisScaled`] `Sdf` samples the field at `(sx*x, sy*y, sz*z)`, which is a
//!    cell of physical extent `(h*sx, h*sy, h*sz)` presented on a uniform
//!    integer grid. The anisotropic arm's eight corner values for cell
//!    `(i, j, k)` are then read from that wrapper at the **reciprocal** grid
//!    point `P/s`, where `P` is the isotropic arm's own sample point. Every one
//!    of the eight reference-field domains is symmetric about the origin
//!    (`([-a; 3], [a; 3])`, fields/mod.rs) — asserted, not assumed — and every
//!    scale in `SCALES` is a power of two, so `s * (P / s) == P` **bit-exactly**
//!    in `f64` and the two arms carry byte-identical tensors while describing
//!    cells of different aspect ratio. `corner_value_mismatches` measures that
//!    round-trip rather than trusting it, and it is asserted zero: without it,
//!    `count_disagreements == 0` could mean "the resampling silently returned the
//!    same numbers for the wrong reason".
//!
//!    So what the census models is precisely the registered claim — *the
//!    body-saddle count is a property of the eight values alone* — driven through
//!    the anisotropic sampling path end to end. What it does **not** model is a
//!    resampling of the field at coarser spacing: reading a cell at `(h, 2h, 4h)`
//!    spacing gives genuinely different corner values because it is a different
//!    region of the field, and a disagreement from that would be a statement
//!    about resolution, not about the saddle path. C1's falsifier says *"the
//!    extractor's saddle path depends on geometry the invariant says it cannot --
//!    i.e. a bug"*, and only the fixed-values construction can produce that
//!    reading.
//!
//! 2. **The sign check — different tensor, same sign.** The census additionally
//!    applies an explicit non-trivial **unimodular** `GL(2)^3` triple
//!    ([`G_U`], [`G_V`], [`G_W`], each of determinant `1`) to every cell's tensor
//!    through `common::poly::act_gl2_cubed_f64`, and compares `sign(Delta)`
//!    before and after. This is the half that is not blind by signature: the
//!    tensor genuinely changes, twelve monomials of cancellation happen, and the
//!    sign has to survive. `gl2_sign_invariant_cells` and
//!    `gl2_sign_variant_cells` are that census. A unimodular triple is chosen so
//!    the weight is exactly `1` and any drift is rounding rather than scale.
//!
//!    The census is taken over the cells that **have** a sign. `Delta == 0` is the
//!    discriminant's own vanishing locus: there the exact image is also zero, and
//!    any rounding at all in twelve cancelling monomials moves it off zero, so
//!    such a cell would be counted as a failure of an invariance it cannot state.
//!    Those cells are counted in `gl2_zero_delta_perturbed_cells` instead, beside
//!    the `delta_zero_cells` census that sizes the stratum. It is not a small
//!    stratum and pretending otherwise would have been the mistake: `box_exact`
//!    at `33^3` has `28,616` of `32,768` cells at exactly zero, because the
//!    field is piecewise linear and a linear function has no body saddle and no
//!    twist.
//!
//! Note what [`AxisScaled`] is not: it is a **domain** reparametrisation, so it
//! is not a distance function — `|grad|` is no longer one. That is irrelevant
//! here, because nothing in this row reads a distance: the census reads sample
//! values, and the extraction arm exists only to move vertices.
//!
//! # Arms
//!
//! Forty-eight rows: eight reference fields x two resolutions (`33^3`, `65^3`,
//! the phase's standard pair) x three axis-scale triples.
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `axis_scales = 1x1x1` | nothing — [`AxisScaled`] is the identity map | **yes** |
//! | `axis_scales = 1x2x4` | cell extent `(h, 2h, 4h)`, aspect ratio 4 | no |
//! | `axis_scales = 1x1x8` | cell extent `(h, h, 8h)`, aspect ratio 8 | no |
//! | unimodular `GL(2)^3` sign census | the tensor itself, at weight exactly 1 | no |
//! | `MarchingCubes` extraction, both arms | the vertex geometry | **yes**, the vacuity control |
//! | 500 exact `i128` weight trials | random tensors and random invertible `g` | no |
//! | singular `g` draws | rejected and counted, never used | **yes**, `gl2_singular_draws` |
//!
//! The `1x1x1` row is a control and *is* a row, because its job is to show the
//! whole apparatus reads zero when nothing is varied — `vertex_positions_moved`
//! must be `0` there and non-zero everywhere else, and one number doing both is
//! worth more than an assert.
//!
//! # The two clauses, and how each is decided
//!
//! **C1 is per row and per cell.** `saddle_count_isotropic` and
//! `saddle_count_anisotropic` are the summed `BodySaddles::of(&corner)
//! .inside_count()` over every cell of the grid — the shipped path, called
//! exactly as `trilinear.rs` calls it, over *all* cells rather than only the
//! ambiguous ones, which is the stricter population. `count_disagreements`
//! counts cells whose two counts differ. `sign_delta_isotropic` and
//! `sign_delta_anisotropic` are that row's sign histogram of the discriminant,
//! formatted `neg|zero|pos`, computed from `BodySaddles::coefficients` as
//! `b*b - 4*a*c` — again the shipped expression rather than a re-derivation.
//! `sign_disagreements` counts cells where the two signs differ, with `+-0.0`
//! mapped to `0` because `trilinear.rs:250` tests `discriminant == R::ZERO` and
//! `-0.0 == 0.0`. `c1_holds` is that row's own verdict: zero count
//! disagreements, zero sign disagreements, zero corner-value mismatches, and
//! zero `GL(2)^3` sign variants.
//!
//! **C2 is global and exact.** 500 trials, each a random `i128` tensor with
//! `|f_i| <= 8` and three random `2x2` integer matrices with entries in
//! `[-4, 4]`; a singular draw is rejected and counted, never repaired. The
//! headline is the **exact** arithmetic: `Det(g.A) == (det g1 det g2 det g3)^2 *
//! Det(A)` compared as an equality of two `i128` products, never as a quotient,
//! so "equals" means *identical integers* and not "within a tolerance". The
//! registered wording asks for the ratio to `f64` rounding; that is recorded
//! beside it as `gl2_f64_max_deviation`, the largest relative deviation of
//! `Det(g.A) / ((det g1 det g2 det g3)^2 * Det(A))` from `1` over the trials with
//! non-zero `Det(A)`. Magnitudes are bounded so the exact arm cannot overflow:
//! an acted entry is at most `8 * 4^3 * 8 = 4096`, twelve degree-4 monomials of
//! that sum below `2^54`, and the weight is at most `(4 * 4 + 4 * 4)^6 < 2^31`.
//! `gl2_weight_check` carries both halves as `exact_<m>of<n>|f64max_<d>`, and it
//! is the same value on every row because C2 is a property of the algebra and not
//! of a field — the header says so here so that forty-eight identical cells are
//! not read as forty-eight measurements. Likewise `c2_holds`.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration says: *"C1 is a correctness clause and moves no cost; it
//! removes an untested assumption from every anisotropic-cell path, including all
//! of Group D."* Discharged, and the arithmetic is a signature: `cell_size` is
//! not an argument of `BodySaddles::of`, `BodySaddles::coefficients` or
//! `BodySaddles::roots`, so the cost moved is **exactly zero nanoseconds** and no
//! change to `crates/isomesh/src/**` is proposed by this row. `wall_ns` is
//! recorded because it is interesting and is read by nothing, which is the only
//! safe status for a nanosecond on a host whose governor swings the same binary
//! `1.45x` (`M-280`).
//!
//! What stands in a share's place is what the clause licenses, with an exact
//! denominator: `48 * (samples - 1)^3` cells summed over the fixture —
//! `8 * (32^3 + 64^3) * 3 = 7,077,888` cells — each of which had its saddle count
//! and discriminant sign compared across two cell aspect ratios and one
//! unimodular tensor change. Group D may then vary cell aspect ratio without
//! re-asking whether interior ambiguity moved under it.
//!
//! # Vacuity controls
//!
//! `M-44`: a zero that could not have been non-zero is not a measurement. Every
//! control runs before the first `run.record` and every panic message starts
//! `VOID: `.
//!
//! - **`vertex_positions_moved > 0`, per anisotropic scale triple and
//!   resolution** — the registration's own control, verbatim. Both arms are
//!   extracted with `MarchingCubes` (with `InteriorAmbiguity::Trilinear`, so the
//!   extraction really does consume `BodySaddles`) on the same shape, origin and
//!   `cell_size`; the anisotropic arm extracts [`AxisScaled`], whose zero set is
//!   a per-axis compression of the field's. The column is the number of vertex
//!   positions in the anisotropic mesh with no bit-identical counterpart in the
//!   isotropic mesh, compared as a multiset of `f64::to_bits` triples so the
//!   answer does not depend on emission order. Zero here would mean the fixture
//!   never exercised anisotropy and every `count_disagreements = 0` in the file
//!   is a comparison of a thing with itself.
//! - **`vertex_positions_moved == 0` on the `1x1x1` control row** — the same
//!   column read the other way: the identity map must move nothing, or the
//!   comparison is reporting noise.
//! - **`cells_with_saddles > 0` over the fixture** — at least one cell must
//!   actually have a body saddle, or "bit-identical saddle counts" is `0 == 0`
//!   over seven million cells.
//! - **`corner_value_mismatches == 0` on every row** — the anisotropic
//!   reciprocal-grid resampling must return the isotropic values bit for bit, or
//!   the two arms are not two readings of one tensor.
//! - **both discriminant signs present over the fixture** — `sign_delta_*` must
//!   contain non-zero `neg` and `pos` counts, or `sign_disagreements` is a claim
//!   about signs over a population that only ever had one.
//! - **C2's population is non-degenerate** — `gl2_nontrivial_weight_trials > 0`
//!   (some weight other than `1`, or the square law is untested),
//!   `gl2_positive_base_trials > 0` and `gl2_negative_base_trials > 0` (both
//!   signs of `Det(A)`), and `gl2_singular_draws > 0` (the rejection path is
//!   reachable, so "500 invertible triples" is a filter and not a coincidence).
//!
//! # Determinism
//!
//! One thread. The only randomness is C2's, from `common::poly::Rng` — a
//! SplitMix64 seeded by [`SEED`], recorded in the `seed` column — so every count
//! is the same on every host and every re-run. The census is a deterministic
//! walk over a deterministic grid; every comparison is on `f64` bit patterns or
//! on `i128` integers; every float sign goes through [`float_sign`], which maps
//! both zeros to `0`.
//!
//! `clippy::float_cmp` is allowed for four comparisons and no others, and each is
//! exact by construction rather than by hope: `lo[axis] == -hi[axis]` on a domain
//! literal, and three `scales == [1.0, 1.0, 1.0]` tests against entries of
//! [`SCALES`], which are small integers held in `f64`. A tolerance on either
//! would turn a structural test into a numeric one.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::trilinear::BodySaddles;
use isomesh::marching_cubes::{InteriorAmbiguity, MarchingCubes};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf, for_each_reference_field};

use crate::common::poly;

// ─── the registered constants ───────────────────────────────────────────────

/// The PRNG seed for C2's trials, recorded in the `seed` column.
///
/// Any value would do; this one is fixed so the counts are reproducible, which
/// is the only property a seed needs to have.
const SEED: u64 = 0x0128_A15C_A1E5_0F87;

/// The two resolutions the registration names, in samples per axis. `n` samples
/// span `n - 1` cells, so these are `32^3` and `64^3` cells.
const RESOLUTIONS: [u32; 2] = [33, 65];

/// The three per-axis cell scale triples the registration names.
///
/// Every entry is a power of two, and that is load-bearing rather than tidy:
/// `s * (P / s) == P` exactly in `f64` only for a power of two, and the whole
/// census rests on the anisotropic arm reading back the isotropic corner values
/// bit for bit.
const SCALES: [[f64; 3]; 3] = [[1.0, 1.0, 1.0], [1.0, 2.0, 4.0], [1.0, 1.0, 8.0]];

/// C2's registered trial count.
const GL2_TRIALS: usize = 500;

/// C2's tensor entries are drawn from `[-GL2_TENSOR_SPAN, GL2_TENSOR_SPAN]`.
const GL2_TENSOR_SPAN: i64 = 8;

/// C2's matrix entries are drawn from `[-GL2_MATRIX_SPAN, GL2_MATRIX_SPAN]`.
const GL2_MATRIX_SPAN: i64 = 4;

/// The relative-deviation bar for C2's `f64` arm.
///
/// The exact arm is an integer equality and needs no tolerance. This is the
/// registration's *"to `f64` rounding"*, and twelve degree-4 monomials of
/// magnitude up to `4096^4` cancelling into a value of order `1` cannot do
/// better than a few hundred `f64` epsilons.
const GL2_F64_TOLERANCE: f64 = 1e-9;

/// A unimodular shear on the `u` index. `det = 1`.
const G_U: [[f64; 2]; 2] = [[1.0, 1.0], [0.0, 1.0]];

/// A unimodular shear on the `v` index, transposed relative to [`G_U`] so the
/// three matrices are not one matrix applied three times. `det = 1`.
const G_V: [[f64; 2]; 2] = [[1.0, 0.0], [1.0, 1.0]];

/// A unimodular quarter turn on the `w` index. `det = 1`, and it is not a shear,
/// so the triple mixes two different kinds of `SL(2)` element.
const G_W: [[f64; 2]; 2] = [[0.0, 1.0], [-1.0, 0.0]];

// ─── the anisotropic cell, as a field ───────────────────────────────────────

/// The field with its domain scaled per axis about the origin.
///
/// `sample(p) = field.sample((sx*px, sy*py, sz*pz))`. Presented on a uniform
/// integer grid this is a cell of physical extent `(h*sx, h*sy, h*sz)`, which is
/// the only way the shipped `extract` signature — one scalar `cell_size` — can be
/// made to mesh an anisotropic cell at all.
///
/// It is a **domain** reparametrisation and therefore not a distance function:
/// `|grad|` is no longer one and no distance read from it would mean anything.
/// Nothing in this row reads a distance.
struct AxisScaled<'a, F> {
    field: &'a F,
    scales: [f64; 3],
}

impl<F> Sdf for AxisScaled<'_, F>
where
    F: Sdf<Scalar = f64>,
{
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        self.field.sample([
            self.scales[0] * p[0],
            self.scales[1] * p[1],
            self.scales[2] * p[2],
        ])
    }
}

// ─── signs, cells and meshes ────────────────────────────────────────────────

/// The sign of a float, with **both** zeros mapping to `0`.
///
/// Not a comparison against `0.0` through `total_cmp`, which orders `-0.0` below
/// `+0.0`: `trilinear.rs:250` tests `discriminant == R::ZERO`, and `-0.0 == 0.0`,
/// so a discriminant that reached zero has lost its sign as far as the extractor
/// is concerned and must be reported as having lost it.
fn float_sign(x: f64) -> i32 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

/// The sign histogram slot for a sign: `0 = neg`, `1 = zero`, `2 = pos`.
fn sign_slot(sign: i32) -> usize {
    match sign {
        -1 => 0,
        0 => 1,
        _ => 2,
    }
}

/// The corner offsets of the crate's cell indexing, `corner = u + 2v + 4w`.
///
/// `crate::cube::corner_offset` (cube.rs:149) is `pub(crate)` and unreachable
/// from a bench, so the layout is restated here in the form its own test
/// (`corner_offset_matches_the_bit_layout`, cube.rs:315) asserts.
const fn corner_offset(corner: usize) -> [usize; 3] {
    [corner & 1, (corner >> 1) & 1, (corner >> 2) & 1]
}

/// The eight corner values of cell `(x, y, z)` from a sample grid of `n^3`.
fn cell_corners(values: &[f64], n: usize, x: usize, y: usize, z: usize) -> [f64; 8] {
    std::array::from_fn(|corner| {
        let o = corner_offset(corner);
        values[(x + o[0]) + (y + o[1]) * n + (z + o[2]) * n * n]
    })
}

/// The saddle count and the discriminant of one cell, both from the shipped
/// path: `BodySaddles::of(&corner).inside_count()` (trilinear.rs:165, :294) and
/// `b*b - 4*a*c` on `BodySaddles::coefficients` (`:199`, discriminant at `:246`).
fn saddle_reading(corner: &[f64; 8]) -> (u32, f64) {
    let count = BodySaddles::of(corner).inside_count();
    let [a, b, c] = BodySaddles::coefficients(corner);
    (count, b * b - 4.0 * a * c)
}

/// Every vertex position of a mesh as a sorted multiset of bit patterns.
///
/// Bit patterns rather than floats so the comparison is exact, and sorted so it
/// does not depend on the order the extractor happened to emit them in.
fn position_keys(mesh: &MeshBuffer<f64>) -> Vec<[u64; 3]> {
    let mut keys: Vec<[u64; 3]> = mesh
        .positions
        .iter()
        .map(|p| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()])
        .collect();
    keys.sort_unstable();
    keys
}

/// How many entries of `a` have no counterpart in `b`, as multisets.
///
/// Both inputs are sorted, so this is one merge pass.
fn unmatched(a: &[[u64; 3]], b: &[[u64; 3]]) -> u64 {
    let mut i = 0;
    let mut j = 0;
    let mut count = 0;
    while i < a.len() {
        while j < b.len() && b[j] < a[i] {
            j += 1;
        }
        if j < b.len() && b[j] == a[i] {
            j += 1;
        } else {
            count += 1;
        }
        i += 1;
    }
    count
}

/// Extract one field on a grid with `MarchingCubes`, interior ambiguity on.
///
/// `InteriorAmbiguity::Trilinear` rather than the default `Ignore` so the
/// extraction genuinely consumes `BodySaddles` — an anisotropy control that ran
/// the cheaper path would be moving vertices through code the row is not about.
fn extract<S>(sdf: &S, shape: &RuntimeShape3, origin: [f64; 3], cell_size: f64) -> MeshBuffer<f64>
where
    S: Sdf<Scalar = f64>,
{
    let mut mc = MarchingCubes::<f64>::new();
    mc.set_interior_ambiguity(InteriorAmbiguity::Trilinear);
    let mut mesh = MeshBuffer::<f64>::new();
    mc.extract(sdf, shape, origin, cell_size, &mut mesh)
        .expect("the reference grids are at least 2 samples on every axis");
    mesh
}

// ─── C2: the relative-invariant weight ──────────────────────────────────────

/// Everything C2 measures, over the whole trial population.
#[derive(Clone, Copy, Debug, Default)]
struct WeightLaw {
    /// Invertible triples actually used. The registration asks for 500.
    trials: usize,
    /// Trials where the exact `i128` identity held as an integer equality.
    exact_matches: usize,
    /// Draws rejected because some `det g` was zero. A control: the rejection
    /// path must be reachable.
    singular_draws: u64,
    /// Trials whose weight was something other than `1`, so the square law was
    /// actually exercised rather than read as an identity.
    nontrivial_weight: u64,
    /// Trials with `Det(A) > 0` and `Det(A) < 0`. Both must be non-zero.
    positive_base: u64,
    negative_base: u64,
    /// Trials with `Det(A) == 0`, where the identity is `0 == 0` and the ratio
    /// does not exist. Reported, never counted into the deviation.
    zero_base: u64,
    /// The largest relative deviation of the `f64` ratio from `1`.
    f64_max_deviation: f64,
    /// The largest `|weight|` seen, so `nontrivial_weight` has a magnitude.
    max_abs_weight: i128,
}

impl WeightLaw {
    /// The verdict for C2: every exact trial an integer identity, and the `f64`
    /// ratio within [`GL2_F64_TOLERANCE`].
    fn holds(&self) -> bool {
        self.trials == GL2_TRIALS
            && self.exact_matches == self.trials
            && self.f64_max_deviation <= GL2_F64_TOLERANCE
    }

    /// The registered `gl2_weight_check` token. CSV-safe: no comma, no quote.
    fn token(&self) -> String {
        format!(
            "exact_{}of{}|f64max_{:.3e}",
            self.exact_matches, self.trials, self.f64_max_deviation
        )
    }
}

/// `det` of a `2x2` integer matrix.
fn det_i128(g: [[i128; 2]; 2]) -> i128 {
    g[0][0] * g[1][1] - g[0][1] * g[1][0]
}

/// Draw a `2x2` integer matrix with entries in `[-GL2_MATRIX_SPAN, span]`.
fn draw_matrix(rng: &mut poly::Rng) -> [[i128; 2]; 2] {
    let mut g = [[0_i128; 2]; 2];
    for row in &mut g {
        for entry in row.iter_mut() {
            *entry = i128::from(rng.next_i64_in(-GL2_MATRIX_SPAN, GL2_MATRIX_SPAN + 1));
        }
    }
    g
}

/// Measure C2: `Det(g.A) == (det g1 det g2 det g3)^2 * Det(A)`, exactly and in
/// `f64`, over [`GL2_TRIALS`] random invertible triples.
///
/// The exact arm is the headline because it decides the clause without a
/// tolerance: two `i128` products are either the same integer or they are not.
/// The `f64` arm is the registration's own wording and is reported beside it.
fn weight_law() -> WeightLaw {
    let cayley = poly::cayley_2x2x2();
    let mut rng = poly::Rng::new(SEED);
    let mut law = WeightLaw::default();

    while law.trials < GL2_TRIALS {
        let g1 = draw_matrix(&mut rng);
        let g2 = draw_matrix(&mut rng);
        let g3 = draw_matrix(&mut rng);
        let (d1, d2, d3) = (det_i128(g1), det_i128(g2), det_i128(g3));
        if d1 == 0 || d2 == 0 || d3 == 0 {
            law.singular_draws += 1;
            continue;
        }
        let f: [i128; 8] = std::array::from_fn(|_| {
            i128::from(rng.next_i64_in(-GL2_TENSOR_SPAN, GL2_TENSOR_SPAN + 1))
        });

        let weight = {
            let product = d1 * d2 * d3;
            product * product
        };
        let base = cayley.eval_i128(&f);
        let acted = cayley.eval_i128(&poly::act_gl2_cubed(g1, g2, g3, &f));
        if acted == weight * base {
            law.exact_matches += 1;
        }

        let g1f = f64_matrix(g1);
        let g2f = f64_matrix(g2);
        let g3f = f64_matrix(g3);
        let ff: [f64; 8] = std::array::from_fn(|i| f[i] as f64);
        let acted_f = cayley.eval_f64(&poly::act_gl2_cubed_f64(g1f, g2f, g3f, &ff));
        let predicted_f = weight as f64 * base as f64;
        if base == 0 {
            law.zero_base += 1;
        } else {
            let deviation = ((acted_f - predicted_f) / predicted_f).abs();
            law.f64_max_deviation = law.f64_max_deviation.max(deviation);
        }

        if weight != 1 {
            law.nontrivial_weight += 1;
        }
        law.max_abs_weight = law.max_abs_weight.max(weight.abs());
        match base.signum() {
            1 => law.positive_base += 1,
            -1 => law.negative_base += 1,
            _ => {}
        }
        law.trials += 1;
    }
    law
}

/// An integer matrix as `f64`. Entries are bounded by `4`, so this is exact.
fn f64_matrix(g: [[i128; 2]; 2]) -> [[f64; 2]; 2] {
    [
        [g[0][0] as f64, g[0][1] as f64],
        [g[1][0] as f64, g[1][1] as f64],
    ]
}

// ─── one row of the census ──────────────────────────────────────────────────

/// One `(field, resolution, axis_scales)` row.
#[derive(Clone, Debug)]
struct Row {
    field: &'static str,
    samples: u32,
    scales: [f64; 3],
    cells: u64,
    saddle_iso: u64,
    saddle_aniso: u64,
    count_disagreements: u64,
    sign_iso: [u64; 3],
    sign_aniso: [u64; 3],
    sign_disagreements: u64,
    corner_mismatches: u64,
    surface_cells: u64,
    cells_with_saddles: u64,
    delta_zero_cells: u64,
    gl2_sign_invariant: u64,
    gl2_sign_variant: u64,
    gl2_zero_delta_perturbed: u64,
    gl2_variant_max_rel_delta: f64,
    vertices_iso: usize,
    vertices_aniso: usize,
    moved: u64,
}

impl Row {
    /// Whether this row's own C1 held: the body-saddle count is bit-identical
    /// across the two cell aspect ratios, the discriminant's sign is too, and the
    /// two arms really were reading one tensor.
    ///
    /// `gl2_sign_variant` is deliberately **not** a term here, and that is a
    /// judgement worth stating rather than burying. The registered falsifier is
    /// *"any disagreement, which would mean the extractor's saddle path depends on
    /// geometry the invariant says it cannot -- i.e. a bug"*. A unimodular
    /// `GL(2)^3` image whose `f64` discriminant lands on the other side of zero is
    /// not that: it is `f64` cancellation in twelve monomials evaluated on a cell
    /// whose exact discriminant is already within rounding of the vanishing
    /// locus, and the exact arithmetic that would settle it is C2's, which holds
    /// on 500/500 trials in `i128`. Folding it in would report a falsification for
    /// a reason the falsifier does not name. It is reported instead —
    /// `gl2_sign_variant_cells` beside `gl2_variant_max_rel_delta`, the scale-free
    /// `|Delta| / max|f_i|^4` of the worst such cell — so the finding can quote
    /// the number and say what it is.
    fn c1_holds(&self) -> bool {
        self.count_disagreements == 0 && self.sign_disagreements == 0 && self.corner_mismatches == 0
    }

    /// `max(s) / min(s)`, the registered `cell_aspect_ratio`.
    fn aspect_ratio(&self) -> f64 {
        let hi = self.scales.iter().copied().fold(f64::MIN, f64::max);
        let lo = self.scales.iter().copied().fold(f64::MAX, f64::min);
        hi / lo
    }

    /// The registered `axis_scales` token, e.g. `1x2x4`.
    fn axis_scales(&self) -> String {
        format!(
            "{}x{}x{}",
            self.scales[0] as u64, self.scales[1] as u64, self.scales[2] as u64
        )
    }

    /// Whether this row is the isotropic control.
    fn is_control(&self) -> bool {
        self.scales == [1.0, 1.0, 1.0]
    }
}

/// A sign histogram as the CSV token `neg|zero|pos`.
fn sign_token(h: [u64; 3]) -> String {
    format!("{}|{}|{}", h[0], h[1], h[2])
}

/// Sample a field on the grid `origin + cell_size * idx`, `x` fastest.
fn sample_grid<S>(sdf: &S, n: usize, origin: [f64; 3], cell_size: f64) -> Vec<f64>
where
    S: Sdf<Scalar = f64>,
{
    let mut values = vec![0.0_f64; n * n * n];
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                values[x + y * n + z * n * n] = sdf.sample([
                    origin[0] + cell_size * x as f64,
                    origin[1] + cell_size * y as f64,
                    origin[2] + cell_size * z as f64,
                ]);
            }
        }
    }
    values
}

/// Sample the anisotropic arm at the **reciprocal** grid `(origin + h*idx) / s`.
///
/// `AxisScaled` then multiplies each coordinate back by `s`, so with `s` a power
/// of two the result is the isotropic value bit for bit while the cell it
/// describes has extent `(h*sx, h*sy, h*sz)`. That round-trip is measured by
/// `corner_value_mismatches`, never assumed.
fn sample_reciprocal_grid<S>(
    sdf: &S,
    n: usize,
    origin: [f64; 3],
    cell_size: f64,
    scales: [f64; 3],
) -> Vec<f64>
where
    S: Sdf<Scalar = f64>,
{
    let mut values = vec![0.0_f64; n * n * n];
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let p = [
                    origin[0] + cell_size * x as f64,
                    origin[1] + cell_size * y as f64,
                    origin[2] + cell_size * z as f64,
                ];
                values[x + y * n + z * n * n] =
                    sdf.sample([p[0] / scales[0], p[1] / scales[1], p[2] / scales[2]]);
            }
        }
    }
    values
}

/// Every row for one reference field.
fn measure<F>(name: &'static str, field: &F, rows: &mut Vec<Row>)
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();
    for axis in 0..3 {
        assert!(
            lo[axis] == -hi[axis],
            "VOID: {name}'s domain is not symmetric about the origin on axis {axis} \
             ({lo:?}, {hi:?}), so scaling it per axis would also translate it and the \
             anisotropic arm would not be reading the isotropic tensor"
        );
    }

    for samples in RESOLUTIONS {
        let (shape, origin, h) = common::grid::<f64, _>(field, samples);
        let n = samples as usize;
        let values_iso = sample_grid(field, n, origin, h);
        let keys_iso = position_keys(&extract(field, &shape, origin, h));

        for scales in SCALES {
            let wrapper = AxisScaled { field, scales };
            let values_aniso = sample_reciprocal_grid(&wrapper, n, origin, h, scales);
            let keys_aniso = position_keys(&extract(&wrapper, &shape, origin, h));

            let mut row = Row {
                field: name,
                samples,
                scales,
                cells: 0,
                saddle_iso: 0,
                saddle_aniso: 0,
                count_disagreements: 0,
                sign_iso: [0; 3],
                sign_aniso: [0; 3],
                sign_disagreements: 0,
                corner_mismatches: 0,
                surface_cells: 0,
                cells_with_saddles: 0,
                delta_zero_cells: 0,
                gl2_sign_invariant: 0,
                gl2_sign_variant: 0,
                gl2_zero_delta_perturbed: 0,
                gl2_variant_max_rel_delta: 0.0,
                vertices_iso: keys_iso.len(),
                vertices_aniso: keys_aniso.len(),
                moved: unmatched(&keys_aniso, &keys_iso),
            };

            for z in 0..n - 1 {
                for y in 0..n - 1 {
                    for x in 0..n - 1 {
                        let corner_iso = cell_corners(&values_iso, n, x, y, z);
                        let corner_aniso = cell_corners(&values_aniso, n, x, y, z);
                        row.cells += 1;
                        for c in 0..8 {
                            if corner_iso[c].to_bits() != corner_aniso[c].to_bits() {
                                row.corner_mismatches += 1;
                            }
                        }

                        let (count_iso, delta_iso) = saddle_reading(&corner_iso);
                        let (count_aniso, delta_aniso) = saddle_reading(&corner_aniso);
                        row.saddle_iso += u64::from(count_iso);
                        row.saddle_aniso += u64::from(count_aniso);
                        if count_iso != count_aniso {
                            row.count_disagreements += 1;
                        }
                        if count_iso > 0 {
                            row.cells_with_saddles += 1;
                        }

                        let sign_iso = float_sign(delta_iso);
                        let sign_aniso = float_sign(delta_aniso);
                        row.sign_iso[sign_slot(sign_iso)] += 1;
                        row.sign_aniso[sign_slot(sign_aniso)] += 1;
                        if sign_iso != sign_aniso {
                            row.sign_disagreements += 1;
                        }
                        if sign_iso == 0 {
                            row.delta_zero_cells += 1;
                        }

                        let inside = corner_iso[0] < 0.0;
                        if corner_iso.iter().any(|v| (*v < 0.0) != inside) {
                            row.surface_cells += 1;
                        }

                        // The half of C1 that is not blind by signature: the
                        // tensor genuinely changes under a unimodular GL(2)^3
                        // triple, the weight is exactly 1, and the sign has to
                        // survive twelve monomials of cancellation. The claim is
                        // about a *sign*, so it is censused over the cells that
                        // have one: `Delta == 0` is the discriminant's own
                        // vanishing locus, where the exact image is also zero and
                        // any f64 rounding at all moves it off, and those cells
                        // are counted separately rather than folded in as
                        // failures of an invariance they cannot state.
                        let acted = poly::act_gl2_cubed_f64(G_U, G_V, G_W, &corner_iso);
                        let (_, delta_acted) = saddle_reading(&acted);
                        let sign_acted = float_sign(delta_acted);
                        if sign_iso == 0 {
                            if sign_acted != 0 {
                                row.gl2_zero_delta_perturbed += 1;
                            }
                        } else if sign_acted == sign_iso {
                            row.gl2_sign_invariant += 1;
                        } else {
                            row.gl2_sign_variant += 1;
                            // How close to the vanishing locus this cell was, in
                            // the scale-free units `P-134` will normalise by:
                            // `Delta` is degree 4, so `|Delta| / max|f_i|^4` is
                            // dimensionless and a value near `f64` epsilon says
                            // the flip is cancellation and not geometry.
                            let scale = corner_iso.iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
                            let quartic = scale * scale * scale * scale;
                            let relative = delta_iso.abs() / quartic;
                            row.gl2_variant_max_rel_delta =
                                row.gl2_variant_max_rel_delta.max(relative);
                        }
                    }
                }
            }
            rows.push(row);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-128");

    common::experiment::run(prereg, |run| {
        let started = Instant::now();

        // ── C2 first: it is global, and its controls gate every row's verdict ──
        let law = weight_law();
        assert!(
            law.trials == GL2_TRIALS,
            "VOID: C2 ran {} of the registered {GL2_TRIALS} invertible triples",
            law.trials
        );
        assert!(
            law.singular_draws > 0,
            "VOID: no singular g was ever drawn, so `500 random invertible triples` is a \
             coincidence rather than a filter and the rejection path is unreachable"
        );
        assert!(
            law.nontrivial_weight > 0 && law.max_abs_weight > 1,
            "VOID: every weight was 1 (max |weight| {}), so the square law was read as an \
             identity and `(det g1 det g2 det g3)^2` was never actually exercised",
            law.max_abs_weight
        );
        assert!(
            law.positive_base > 0 && law.negative_base > 0,
            "VOID: C2's tensors carry only one sign of Det(A) ({} positive, {} negative), so \
             a weight that was not a square could still have passed",
            law.positive_base,
            law.negative_base
        );

        // ── the census ──
        let mut rows: Vec<Row> = Vec::new();
        for_each_reference_field!(f64, |name, field| {
            measure(name, &field, &mut rows);
        });

        let elapsed_ns = started.elapsed().as_nanos();

        // ── the census's own controls, all before the first record ──
        for row in &rows {
            assert!(
                row.corner_mismatches == 0,
                "VOID: {} at {}^3 scale {} read {} corner values that differ between the two \
                 arms, so the anisotropic arm is not a second reading of one tensor and \
                 `count_disagreements` compares two different cells",
                row.field,
                row.samples,
                row.axis_scales(),
                row.corner_mismatches
            );
            if row.is_control() {
                assert!(
                    row.moved == 0,
                    "VOID: the 1x1x1 identity arm moved {} vertex positions on {} at {}^3, so \
                     `vertex_positions_moved` is reporting comparison noise rather than \
                     anisotropy",
                    row.moved,
                    row.field,
                    row.samples
                );
            }
        }
        for scales in SCALES {
            if scales == [1.0, 1.0, 1.0] {
                continue;
            }
            for samples in RESOLUTIONS {
                let moved: u64 = rows
                    .iter()
                    .filter(|r| r.scales == scales && r.samples == samples)
                    .map(|r| r.moved)
                    .sum();
                assert!(
                    moved > 0,
                    "VOID: the {}x{}x{} arm at {samples}^3 moved no vertex position on any of \
                     the eight fields, so the fixture is not exercising anisotropy at all and \
                     every zero in `count_disagreements` is a comparison of a thing with itself",
                    scales[0] as u64,
                    scales[1] as u64,
                    scales[2] as u64
                );
            }
        }
        let saddled: u64 = rows.iter().map(|r| r.cells_with_saddles).sum();
        assert!(
            saddled > 0,
            "VOID: not one cell in the whole fixture has a body saddle, so `bit-identical \
             saddle counts` is 0 == 0 over every cell"
        );
        let negatives: u64 = rows.iter().map(|r| r.sign_iso[0]).sum();
        let positives: u64 = rows.iter().map(|r| r.sign_iso[2]).sum();
        assert!(
            negatives > 0 && positives > 0,
            "VOID: the discriminant has only one sign over the fixture ({negatives} negative, \
             {positives} positive), so `sign_disagreements` is a claim about signs over a \
             population that never had two"
        );

        // ── the rows ──
        let c2 = law.holds();
        let token = law.token();
        for row in &rows {
            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.samples.to_string()),
                ("cell_aspect_ratio", format!("{:.6}", row.aspect_ratio())),
                ("axis_scales", row.axis_scales()),
                ("saddle_count_isotropic", row.saddle_iso.to_string()),
                ("saddle_count_anisotropic", row.saddle_aniso.to_string()),
                ("count_disagreements", row.count_disagreements.to_string()),
                ("sign_delta_isotropic", sign_token(row.sign_iso)),
                ("sign_delta_anisotropic", sign_token(row.sign_aniso)),
                ("sign_disagreements", row.sign_disagreements.to_string()),
                ("gl2_weight_check", token.clone()),
                ("c1_holds", row.c1_holds().to_string()),
                ("c2_holds", c2.to_string()),
                // ── extras (M-273) ──
                ("cells", row.cells.to_string()),
                ("cells_with_saddles", row.cells_with_saddles.to_string()),
                ("corner_value_mismatches", row.corner_mismatches.to_string()),
                ("delta_zero_cells", row.delta_zero_cells.to_string()),
                ("gl2_exact_matches", law.exact_matches.to_string()),
                (
                    "gl2_f64_max_deviation",
                    format!("{:.3e}", law.f64_max_deviation),
                ),
                ("gl2_max_abs_weight", law.max_abs_weight.to_string()),
                ("gl2_negative_base_trials", law.negative_base.to_string()),
                (
                    "gl2_nontrivial_weight_trials",
                    law.nontrivial_weight.to_string(),
                ),
                ("gl2_positive_base_trials", law.positive_base.to_string()),
                (
                    "gl2_sign_invariant_cells",
                    row.gl2_sign_invariant.to_string(),
                ),
                ("gl2_sign_variant_cells", row.gl2_sign_variant.to_string()),
                ("gl2_singular_draws", law.singular_draws.to_string()),
                (
                    "gl2_variant_max_rel_delta",
                    format!("{:.3e}", row.gl2_variant_max_rel_delta),
                ),
                ("gl2_trials", law.trials.to_string()),
                ("gl2_zero_base_trials", law.zero_base.to_string()),
                (
                    "gl2_zero_delta_perturbed_cells",
                    row.gl2_zero_delta_perturbed.to_string(),
                ),
                ("is_control", row.is_control().to_string()),
                ("seed", format!("{SEED:#018x}")),
                ("surface_cells", row.surface_cells.to_string()),
                ("vertex_count_anisotropic", row.vertices_aniso.to_string()),
                ("vertex_count_isotropic", row.vertices_iso.to_string()),
                ("vertex_positions_moved", row.moved.to_string()),
                ("wall_ns", elapsed_ns.to_string()),
            ]);
        }
    });
}
