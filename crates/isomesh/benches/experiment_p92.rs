//! **P-92 — regenerate, do not transmit.**
//!
//! Ticket: R-092. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p92
//! ```
//!
//! Writes `docs/experiments/p-92.csv`.
//!
//! # The bar, and where it comes from
//!
//! meshoptimizer decodes sponza — 184k vertices, 262k triangles — in **1.92 ms
//! on a Core i7-8650U at 2 GHz**, against Draco's 169 ms. `1.92e6 ns / 262e3
//! triangles` is **7.33 ns per triangle, single-threaded, to decode**, and
//! decoding produces geometry and nothing else. The registration rounds it to
//! **7.3** and that is the number C1 is scored against, because a bar that moves
//! with the harness is not a bar.
//!
//! # SHARE, recomputed before the harness was written
//!
//! `✗51`'s rule: a clause stated against a cost has to name the share of that
//! cost it can move, and the check happens *before* the run. C1 is not a
//! speedup ratio, so the arithmetic is a reachability check on the two readings
//! the registration's own text contains — and **they do not agree**.
//!
//! A 33³ chunk is **35,937 samples** and **32,768 cells**. Two committed
//! per-sample marginals for `marching_cubes` on this Zen 3: **10.68 ns/sample**
//! (`✗51`'s extraction marginal) and **13.1892 ns/sample**
//! (`docs/measurements/resolution_sweep-ryzen9-5900x.csv`, 9 rows, 16³–256³).
//! So one 33³ chunk extraction costs **384 µs to 474 µs**.
//!
//! - **The TOTAL reading** — which is what the registration's gloss *"re-extracting
//!   is cheaper than decoding an encoded mesh of the same output"* actually means,
//!   because a receiving endpoint that decodes does **not** sample the field —
//!   needs `total / triangles ≤ 7.3 ns`, i.e. **T ≥ 384,000/7.3 = 52,576
//!   triangles** from 32,768 cells: **1.60 triangles per cell averaged over the
//!   whole volume**. A marching-cubes surface is `O(n²)`: a few thousand active
//!   cells of 32,768, at most 5 triangles each, so the *ceiling* is ≈16,500 and a
//!   realistic surface is ≈4,000. **The total reading is unreachable at 33³ on
//!   this machine by 3.2× at the theoretical ceiling and ≈13× in practice, and
//!   that is arithmetic, not a measurement.** It is written down here before the
//!   run because that is the rule.
//! - **The MARGINAL reading** — the registration's literal wording, *"Marginal
//!   extraction cost per triangle"* — puts the grid sampling and the cell-
//!   classification walk in the fitted fixed term and leaves only the emit stage
//!   in `b`. That is reachable, its value is genuinely unknown, and it is what
//!   this harness fits. `triangle_term_share` reports `b·T/(a + b·T)` so the
//!   reader can see how much of a chunk extraction the clause is talking about.
//!
//! Both numbers are reported on every C1 row: `ns_per_triangle` is the fitted
//! marginal (the registered clause), `ns_per_triangle_total` is
//! `extract_ms·1e6/triangles`.
//!
//! # How the marginal is fitted, and the two traps M-19/M-21 name
//!
//! **Not by dividing a total by a count.** Within one family the grid is fixed at
//! 33³ and the *per-sample* cost is held bit-identical, so the only thing that
//! moves is the triangle count, and `t = a + b·T` separates the two terms. The
//! knob is chosen per family precisely so the sample loop's instruction sequence
//! does not change:
//!
//! | family | knob | why the per-sample cost is constant |
//! |---|---|---|
//! | `box_exact` | half-extent | `min`/`max`/`abs`, no `libm` call at all — `✗51` classified this as the crate's cheapest field, which is what makes it the **lower bound** |
//! | `sphere` | radius | one `sqrt` whatever the radius |
//! | `gyroid` | `scale` | the six `libm` calls run at every scale; scaling the field is scaling the surface, so triangles-per-active-cell moves as little as any knob can |
//! | `fbm_terrain` | `amplitude` | a trailing multiply; the same four octaves and the same lattice either way |
//! | `gyroid_dug30` | `scale` | the thirty subtract-spheres are walked at every sample regardless |
//!
//! **The knob was not free to choose, and the first attempt got it wrong.** The
//! original sweep used gyroid's `iso` and fbm's `frequency`, both of which change
//! the surface's *topology* rather than only its area — and triangles-per-active-
//! cell then moves with the knob, so time is not linear in triangles. The fits
//! came back at `r² = 0.19` and `0.37` and `gyroid_dug30/f64` produced a fitted
//! fixed cost of **−9.32 ms**. That is `M-21`'s signature and it is why the
//! knobs above are area knobs.
//!
//! Two traps, both from `FINDINGS.md`:
//!
//! - **M-19**: a fitted coefficient means nothing until it is compared to the
//!   data's own range. `fit_a_ms` is reported and so is its share of the largest
//!   and smallest run in the family.
//! - **M-21**: a physically impossible fitted parameter is the model telling you
//!   it is wrong. A **negative** `a` here would mean the fixed cost of sampling
//!   35,937 grid points came out below zero, which is nonsense. So every family
//!   carries `fit_sound` — `a > 0` **and** `r² ≥ 0.95` — and the harness asserts
//!   that **at least one of the four cheap arms** (`box_exact` and `sphere`, both
//!   precisions) is sound, because the verdict rests on a lower bound and
//!   without one there is nothing to rest on. Every other family's soundness is
//!   reported rather than asserted: this machine is shared with a dozen sibling
//!   agents running `cargo`, and *which* of the expensive arms keeps its fit
//!   moves run to run. A slope with `fit_sound = false` beside it is not a
//!   marginal cost and must not be quoted as one.
//!
//! # The vacuity control the registration names
//!
//! > *the encoded-size arm must use a real encoder on the real triangles, not an
//! > estimate, and its decode time must be measured on this machine rather than
//! > quoted from the blog post.*
//!
//! So `meshopt` 0.6.2 — the Rust FFI binding to zeux's own meshoptimizer C
//! library — is a **dev**-dependency, and every row on this CSV runs the real
//! codec over that row's real triangles:
//!
//! - `meshopt_roundtrip_differing` — vertices that do not survive
//!   encode→decode, plus triangles that do not survive **up to rotation** (the
//!   index codec normalises each triangle's rotation, which is documented and is
//!   why a naïve `Vec<u32>` equality fails). **Asserted zero.**
//! - `meshopt_control_differing` — the same comparison after bumping **one**
//!   quantised position component by **one** unit in the last place. **Asserted
//!   exactly 1**: a zero that could not have been non-zero is not a
//!   measurement (`M-44`), and one that could not have been *exactly* one is a
//!   comparator that would not notice a formula error either.
//! - `meshopt_bytes_per_triangle` varies row to row. An estimate cannot do that,
//!   and the harness asserts the spread is non-zero.
//! - `meshopt_decode_ms` and `meshopt_ns_per_triangle` are **measured here**, on
//!   this machine, over `DECODE_REPS` decodes into fresh buffers, and `c1_holds`
//!   is scored against the *published* 7.3 while `c1_holds_vs_measured` is scored
//!   against the local figure, which on a 4-plus-GHz Zen 3 is the harsher of the
//!   two.
//!
//! # The other three controls
//!
//! - **`regen_control_differing`** — the field-plus-log byte string is not a size
//!   estimate either: it is decoded back, the chunk is re-extracted from the
//!   decoded description, and `bytes_differing` counts the mesh bytes that moved.
//!   **Asserted zero**, with a control that flips one bit of the encoded
//!   `cell_size`'s exponent and asserts the difference is **non-zero**. That is
//!   the whole claim of this experiment — *the bytes are sufficient to regenerate
//!   the chunk* — instrumented on every row rather than argued.
//! - **`brushes_biting`** — every brush in a C2 log must have changed the mesh.
//!   A log padded with no-ops would make the size arm compare a fat description
//!   against geometry it did not produce. **Asserted equal to `brushes`**, and
//!   earned by construction: the log is drawn by rejection and `brushes_drawn`
//!   records how many candidates it took. The first version drew blind and this
//!   control fired — *44 of 45* — which is `P-94`'s collapse showing up
//!   uninvited.
//! - **`replay_min_triangles`** — the smallest triangle count over the 10,000
//!   extractions of the C3 replay. A digest stream of 10,000 hashes of the empty
//!   mesh would agree perfectly across two machines and mean nothing.
//!   **Asserted positive.**
//!
//! # C3, and why the comparison set is never empty
//!
//! The replay is **10,000 edits**: 200 chunks × 50 subtract-spheres, each brush
//! placed strictly inside its own chunk so a chunk's log is exactly its own
//! brushes and the window needs no approximation. After every edit the chunk is
//! re-extracted in `f64` and hashed with **`isomesh::validate::mesh_hash`** —
//! `M-31`'s own instrument, the one the 216 golden hashes are taken with — and
//! the 8 bytes are appended to a digest stream. The stream is written to
//! `target/experiment_p92_replay-<machine-slug>.bin`.
//!
//! There is one code path: **compare the local stream against every other stream
//! in that directory.** Zero peers is a comparison over an empty set, not a
//! branch. And the set is never empty, because the harness also writes
//! `experiment_p92_replay-control-onebyteflipped.bin` — its own stream with one
//! byte flipped — so a live comparator is proved on every run:
//!
//! - the control peer's `bytes_differing` is **asserted non-zero**;
//! - every other peer's `bytes_differing` is **asserted zero**, which is C3.
//!
//! The cross-machine row therefore appears only on the machine that has the
//! other machine's stream on disk. Copy `target/experiment_p92_replay-apple-m5.bin`
//! next to the local one and re-run; `peer` names which stream each row compared
//! against.
//!
//! On a C3 row the geometry columns — `triangles`, `meshopt_bytes`,
//! `field_plus_log_bytes`, `size_ratio` — describe the **final chunk**, the one
//! whose triangles were encoded and whose description was serialised, so that
//! four columns are four facts about one mesh. `extract_ms` is the **mean
//! per-edit** re-extraction over the trace, `ns_per_triangle` is the replay's own
//! `replay_ms · 1e6 / replay_total_triangles`, and the two aggregates are on the
//! row under their own names. The first version put the 148-million-triangle
//! replay total in `triangles` beside the final chunk's 17 kB encoding and
//! produced a `meshopt_bytes_per_triangle` of **0.0005**, which is the `M-377`
//! defect — a column carrying one arm's number on another arm's row.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::print_literal,
    clippy::too_many_lines
)]

mod common;

use std::path::PathBuf;
use std::time::Instant;

use isomesh::brush::{Brush, BrushStack};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::extractor::Extractor;
use isomesh::fields::{BoxExact, FbmTerrain, Gyroid, Sphere};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};
use meshopt::encoding::{
    decode_index_buffer, decode_vertex_buffer, encode_index_buffer, encode_vertex_buffer,
};
use meshopt::optimize::{optimize_vertex_cache, optimize_vertex_fetch};

// ─── the fixture ────────────────────────────────────────────────────────────

/// Cells per chunk axis. C1 names 33³ chunks, which is 32 cells.
const CHUNK_CELLS: u32 = 32;

/// Samples per chunk axis: `CHUNK_CELLS + 1`.
const CHUNK_SAMPLES: u32 = CHUNK_CELLS + 1;

/// Cell size. A power of two, per `M-32`: chunk seams are bit-exact only when
/// `origin + h·local` is exact, and the C3 replay walks 200 chunks.
const CELL_SIZE: f64 = 0.125;

/// One chunk's world extent: `CHUNK_CELLS · CELL_SIZE`.
const CHUNK_EXTENT: f64 = CHUNK_CELLS as f64 * CELL_SIZE;

/// meshoptimizer's own sponza figure as a per-triangle decode cost, and the bar
/// C1 is scored against: 262k triangles in 1.92 ms is 7.33 ns/triangle, which
/// the registration rounds to 7.3.
const PUBLISHED_DECODE_NS_PER_TRIANGLE: f64 = 7.3;

/// Timed extractions per point **per pass**. The fastest of all
/// `REPS · PASSES` is the point's time.
///
/// **The minimum, not the median, and spread over passes.** The finding is a
/// *slope*, and a slope over eight points is far more fragile than a single
/// time. `amd-pstate-epp` on `powersave` spans 1.96–5.62 GHz on this host
/// (`M-280`), a dozen sibling agents run `cargo` on it, and every excursion is
/// one-sided: interference can only make an extraction slower, never faster. A
/// median still carries that noise on a 250 µs run, and the first version of
/// this harness measured it — `sphere/f32` came back with `r² = 0.658` where the
/// identical `f64` fixture gave `0.99963`, and `gyroid_dug30/f64` produced a
/// **negative** fitted fixed cost. The minimum is the standard estimator for
/// "this cost, without a competitor". Every row still carries `ghz` and
/// `worst_ratio` so the clock is on the artefact either way.
const REPS: usize = 5;

/// Round-robin passes over a family's knob points.
///
/// Fifteen, so every point is sampled seventy-five times at fifteen separate
/// moments. See [`sweep`] for the runs that made this necessary: at seven passes
/// `box_exact/f64` still lost its fit to a neighbour's compile (`r² = 0.80096`)
/// while `box_exact/f32` reached `0.99999` in the same run.
const PASSES: usize = 15;

/// Decodes per row, median taken. Cheap enough to afford many, and the decode is
/// the number the registration insists be measured here rather than quoted.
const DECODE_REPS: usize = 31;

/// `M-50`'s four log buckets, at their **upper** ends — the worst case for C2,
/// since a longer log is a bigger description against the same geometry.
const LOG_BUCKETS: [u32; 4] = [15, 30, 45, 60];

/// C2's bar.
const C2_SIZE_RATIO_BAR: f64 = 20.0;

/// Brushes in the C1 `gyroid_dug30` family. Inside `M-50`'s 16–30 bucket, and
/// held fixed so the family's per-sample cost does not move with the knob.
const C1_DUG_BRUSHES: u32 = 30;

/// Chunks in the C3 replay.
const REPLAY_CHUNKS: usize = 200;

/// Edits per chunk in the C3 replay.
const REPLAY_EDITS_PER_CHUNK: usize = 50;

/// The registered 10⁴-edit trace.
const REPLAY_EDITS: usize = REPLAY_CHUNKS * REPLAY_EDITS_PER_CHUNK;

/// The C3 base field's frequency.
///
/// 2.0 gives a period of `π ≈ 3.14`, strictly under the 4.0-unit chunk, so
/// **every** chunk in the replay contains gyroid surface. That is what makes
/// `replay_min_triangles > 0` a property of the fixture rather than a hope.
const REPLAY_GYROID_SCALE: f64 = 2.0;

/// Brush radius in the replay and in the C2 logs, in cells.
const BRUSH_CELLS: f64 = 3.0;

// ─── a scalar that can be put on a wire ─────────────────────────────────────

/// The two scalars this crate extracts in, as bytes.
///
/// Needed twice: the field-plus-log description is a byte string whose size is
/// C2's whole left-hand side, and the mesh comparison that proves the
/// description is sufficient counts **bytes**, not floats.
trait Wire: Real {
    /// Bytes one scalar occupies on the wire.
    const WIDTH: usize;
    /// Append the little-endian bits.
    fn put(self, out: &mut Vec<u8>);
    /// Read the little-endian bits at `at`, advancing it.
    fn get(bytes: &[u8], at: &mut usize) -> Self;
}

impl Wire for f32 {
    const WIDTH: usize = 4;

    fn put(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.to_le_bytes());
    }

    fn get(bytes: &[u8], at: &mut usize) -> Self {
        let mut b = [0u8; 4];
        b.copy_from_slice(&bytes[*at..*at + 4]);
        *at += 4;
        Self::from_le_bytes(b)
    }
}

impl Wire for f64 {
    const WIDTH: usize = 8;

    fn put(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.to_le_bytes());
    }

    fn get(bytes: &[u8], at: &mut usize) -> Self {
        let mut b = [0u8; 8];
        b.copy_from_slice(&bytes[*at..*at + 8]);
        *at += 8;
        Self::from_le_bytes(b)
    }
}

// ─── the thing that is transmitted instead of geometry ──────────────────────

/// A chunk's regenerable description: which field, where, and what was dug.
///
/// This is the left-hand side of C2 and it is a real byte string, not a size
/// estimate — [`encode_desc`] writes it and [`decode_desc`] reads it back, and
/// the row's `bytes_differing` is what the re-extracted mesh differs from the
/// original by.
#[derive(Clone, Debug, PartialEq)]
struct Desc<R: Real> {
    /// Which base field. 0 sphere, 1 gyroid, 2 fbm terrain, 3 box.
    tag: u8,
    /// The base field's integer parameters: fbm's seed and octave count. Written
    /// unconditionally, which costs the other two families eight bytes each and
    /// is the conservative direction for C2.
    seed: u32,
    /// Octaves, as above.
    octaves: u32,
    /// The base field's scalar parameters, in a family-fixed order.
    params: Vec<R>,
    /// Cells per chunk axis.
    cells: u32,
    /// The chunk's sample origin.
    origin: [R; 3],
    /// The grid spacing.
    cell_size: R,
    /// The edit log: subtract-spheres, in order.
    brushes: Vec<Sphere<R>>,
}

impl<R: Wire> Desc<R> {
    /// Bytes this description occupies on the wire.
    ///
    /// Derived from the encoder rather than asserted against it: `encoded_len`
    /// and `encode_desc` disagreeing is the `M-236` failure — a header that
    /// states a count the body does not produce — so the harness asserts they
    /// agree instead of trusting either.
    fn encoded_len(&self) -> usize {
        1 + 4
            + 4
            + 4
            + 2
            + R::WIDTH * (self.params.len() + 4)
            + 2
            + self.brushes.len() * (1 + 4 * R::WIDTH)
    }
}

/// Write a description, and say where the `cell_size` scalar landed.
///
/// The offset is the control's target: flipping one bit of that scalar's
/// exponent has to change the regenerated mesh, or the byte string is not what
/// the regeneration read.
fn encode_desc<R: Wire>(d: &Desc<R>) -> (Vec<u8>, usize) {
    let mut out = Vec::with_capacity(d.encoded_len());
    out.push(d.tag);
    out.extend_from_slice(&d.seed.to_le_bytes());
    out.extend_from_slice(&d.octaves.to_le_bytes());
    out.extend_from_slice(&d.cells.to_le_bytes());
    out.extend_from_slice(&(d.params.len() as u16).to_le_bytes());
    for p in &d.params {
        p.put(&mut out);
    }
    for o in d.origin {
        o.put(&mut out);
    }
    let cell_size_at = out.len();
    d.cell_size.put(&mut out);
    out.extend_from_slice(&(d.brushes.len() as u16).to_le_bytes());
    for b in &d.brushes {
        // The op byte. Every brush in every trace here is a subtract, and the
        // byte is written anyway: a log that could only hold one operation is
        // not the log this crate has.
        out.push(1u8);
        for c in b.center {
            c.put(&mut out);
        }
        b.radius.put(&mut out);
    }
    // The exponent's low byte is the last one of a little-endian float.
    (out, cell_size_at + R::WIDTH - 1)
}

/// Read a description back.
fn decode_desc<R: Wire>(bytes: &[u8]) -> Desc<R> {
    let mut at = 0usize;
    let tag = bytes[at];
    at += 1;
    let u32_at = |at: &mut usize| {
        let mut b = [0u8; 4];
        b.copy_from_slice(&bytes[*at..*at + 4]);
        *at += 4;
        u32::from_le_bytes(b)
    };
    let seed = u32_at(&mut at);
    let octaves = u32_at(&mut at);
    let cells = u32_at(&mut at);
    let u16_at = |at: &mut usize| {
        let mut b = [0u8; 2];
        b.copy_from_slice(&bytes[*at..*at + 2]);
        *at += 2;
        u16::from_le_bytes(b)
    };
    let param_count = u16_at(&mut at) as usize;
    let params: Vec<R> = (0..param_count).map(|_| R::get(bytes, &mut at)).collect();
    let origin = [
        R::get(bytes, &mut at),
        R::get(bytes, &mut at),
        R::get(bytes, &mut at),
    ];
    let cell_size = R::get(bytes, &mut at);
    let brush_count = u16_at(&mut at) as usize;
    let brushes: Vec<Sphere<R>> = (0..brush_count)
        .map(|_| {
            let op = bytes[at];
            at += 1;
            assert_eq!(op, 1u8, "the only operation this trace writes is subtract");
            let center = [
                R::get(bytes, &mut at),
                R::get(bytes, &mut at),
                R::get(bytes, &mut at),
            ];
            let radius = R::get(bytes, &mut at);
            Sphere { center, radius }
        })
        .collect();
    Desc {
        tag,
        seed,
        octaves,
        params,
        cells,
        origin,
        cell_size,
        brushes,
    }
}

// ─── the fields, rebuilt from a description ─────────────────────────────────

/// Where a chunk is and how finely it is sampled.
///
/// Bundled because every call that needs one needs all three, and eight loose
/// arguments is the shape a mistake hides in.
struct Grid<R: Real> {
    /// Samples per axis.
    shape: RuntimeShape3,
    /// World position of sample `[0, 0, 0]`.
    origin: [R; 3],
    /// Spacing between adjacent samples.
    cell_size: R,
}

/// The grid a description names.
fn grid_of<R: Wire>(desc: &Desc<R>) -> Grid<R> {
    Grid {
        shape: RuntimeShape3::new([desc.cells + 1; 3]).expect("chunk grid fits u32"),
        origin: desc.origin,
        cell_size: desc.cell_size,
    }
}

/// A base field with a log of subtract-spheres applied.
///
/// Written out rather than composed through [`BrushStack`] for the base cases so
/// that a family with an empty log pays no walk at all — the sample loop's
/// instruction sequence has to be constant *within* a family, and it is, but it
/// must also be the crate's real one for the number to mean anything.
fn extract_with<R: Real, S: Sdf<Scalar = R>>(
    mc: &mut MarchingCubes<R>,
    base: &S,
    brushes: &[Brush<Sphere<R>>],
    grid: &Grid<R>,
    out: &mut MeshBuffer<R>,
) {
    out.positions.clear();
    out.normals.clear();
    out.indices.clear();
    if brushes.is_empty() {
        mc.extract_into(base, &grid.shape, grid.origin, grid.cell_size, out)
            .expect("33^3 chunk extracts");
    } else {
        let stack = BrushStack { base, brushes };
        mc.extract_into(&stack, &grid.shape, grid.origin, grid.cell_size, out)
            .expect("33^3 chunk extracts");
    }
}

/// Every byte of a mesh, in the order `mesh_hash` reads them.
fn mesh_bytes<R: Wire>(m: &MeshBuffer<R>) -> Vec<u8> {
    let mut out = Vec::with_capacity((m.positions.len() * 6) * R::WIDTH + m.indices.len() * 4);
    for p in &m.positions {
        for v in p {
            v.put(&mut out);
        }
    }
    for n in &m.normals {
        for v in n {
            v.put(&mut out);
        }
    }
    for i in &m.indices {
        out.extend_from_slice(&i.to_le_bytes());
    }
    out
}

/// Bytes that differ between two byte strings, counting a length mismatch as a
/// difference in every byte the shorter one does not have.
fn bytes_differing(a: &[u8], b: &[u8]) -> usize {
    let common = a.len().min(b.len());
    let mismatched = a[..common]
        .iter()
        .zip(&b[..common])
        .filter(|(x, y)| x != y)
        .count();
    mismatched + (a.len().max(b.len()) - common)
}

// ─── the opponent: a real encoder on the real triangles ─────────────────────

/// Position quantised to 16 bits over the chunk's own box, normal octahedral in
/// two bytes. Eight bytes a vertex.
///
/// This is zeux's own `PackedVertexOct` with the two UV components dropped,
/// because this crate emits none and charging meshopt four bytes a vertex for
/// texture coordinates it was not given would be charging the opponent for a
/// column it did not play. The smaller the encoded arm, the harder C2's 20×.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Packed {
    /// Position, 16-bit fixed point over the chunk box.
    p: [u16; 3],
    /// Octahedral normal.
    n: [u8; 2],
}

/// Quantise a mesh into the packed vertex stream a real pipeline would encode.
fn pack<R: Real>(m: &MeshBuffer<R>, lo: [f32; 3], extent: f32) -> Vec<Packed> {
    let inv = 65535.0f32 / extent;
    m.positions
        .iter()
        .zip(&m.normals)
        .map(|(p, n)| {
            let q = |i: usize| {
                let v = (p[i].as_f32() - lo[i]) * inv;
                v.clamp(0.0, 65535.0).round() as u16
            };
            let (nx, ny, nz) = (n[0].as_f32(), n[1].as_f32(), n[2].as_f32());
            let sum = nx.abs() + ny.abs() + nz.abs();
            let (ax, ay) = (nx / sum, ny / sum);
            let (u, v) = if nz >= 0.0 {
                (ax, ay)
            } else {
                (
                    (1.0 - ay.abs()) * if ax >= 0.0 { 1.0 } else { -1.0 },
                    (1.0 - ax.abs()) * if ay >= 0.0 { 1.0 } else { -1.0 },
                )
            };
            let s = |x: f32| (x.clamp(-1.0, 1.0) * 127.0).round() as i8 as u8;
            Packed {
                p: [q(0), q(1), q(2)],
                n: [s(u), s(v)],
            }
        })
        .collect()
}

/// What the real encoder did to one chunk's real triangles.
struct Encoded {
    /// Encoded vertex stream bytes.
    vertex_bytes: usize,
    /// Encoded index stream bytes.
    index_bytes: usize,
    /// The same mesh with `f32` positions and normals, uncompressed by
    /// quantisation, encoded — reported so the size arm cannot be accused of
    /// having won by throwing precision away.
    f32_bytes: usize,
    /// Median decode of both streams, measured here, in milliseconds.
    decode_ms: f64,
    /// Vertices plus triangles that did not survive encode → decode. Zero.
    roundtrip_differing: usize,
    /// The same, after one quantised position component is moved one unit in the
    /// last place. Exactly one.
    control_differing: usize,
}

/// Encode a mesh, time its decode here, and prove both with a control.
fn encode_and_time<R: Real>(m: &MeshBuffer<R>, lo: [f32; 3], extent: f32) -> Encoded {
    let verts = pack(m, lo, extent);
    let count = verts.len();
    assert!(count > 0, "VOID: nothing to encode");

    // Cache-then-fetch, as meshoptimizer's own documentation requires for the
    // index codec to reach its published density. The opponent plays its best
    // line.
    let mut indices = optimize_vertex_cache(&m.indices, count);
    let verts = optimize_vertex_fetch(&mut indices, &verts);

    let vb = encode_vertex_buffer(&verts).expect("vertex encode");
    let ib = encode_index_buffer(&indices, count).expect("index encode");

    let mut f32_stream: Vec<[f32; 6]> = Vec::with_capacity(count);
    for (p, n) in m.positions.iter().zip(&m.normals) {
        f32_stream.push([
            p[0].as_f32(),
            p[1].as_f32(),
            p[2].as_f32(),
            n[0].as_f32(),
            n[1].as_f32(),
            n[2].as_f32(),
        ]);
    }
    let f32_vb = encode_vertex_buffer(&f32_stream).expect("f32 vertex encode");

    let mut times = Vec::with_capacity(DECODE_REPS);
    for _ in 0..DECODE_REPS {
        let t = Instant::now();
        let dv: Vec<Packed> = decode_vertex_buffer(&vb, count).expect("vertex decode");
        let di: Vec<u32> = decode_index_buffer(&ib, indices.len()).expect("index decode");
        times.push(t.elapsed().as_nanos() as f64);
        std::hint::black_box((&dv[0], di[0]));
    }
    times.sort_unstable_by(f64::total_cmp);

    let dv: Vec<Packed> = decode_vertex_buffer(&vb, count).expect("vertex decode");
    let di: Vec<u32> = decode_index_buffer(&ib, indices.len()).expect("index decode");
    let roundtrip_differing = differing_vertices(&verts, &dv) + differing_triangles(&indices, &di);

    // The control. One unit in the last place on one component of one vertex,
    // through the real encoder and the real decoder, compared by the same
    // comparator. Exactly one vertex must move.
    let mut nudged = verts.clone();
    nudged[0].p[0] = nudged[0].p[0].wrapping_add(1);
    let nvb = encode_vertex_buffer(&nudged).expect("control encode");
    let ndv: Vec<Packed> = decode_vertex_buffer(&nvb, count).expect("control decode");
    let control_differing = differing_vertices(&verts, &ndv);

    Encoded {
        vertex_bytes: vb.len(),
        index_bytes: ib.len(),
        f32_bytes: f32_vb.len() + ib.len(),
        decode_ms: times[times.len() / 2] / 1.0e6,
        roundtrip_differing,
        control_differing,
    }
}

/// Vertices that differ.
fn differing_vertices(a: &[Packed], b: &[Packed]) -> usize {
    assert_eq!(a.len(), b.len(), "decode returned a different vertex count");
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

/// Triangles that differ **up to rotation**.
///
/// meshoptimizer's index codec normalises each triangle's rotation, which is
/// documented behaviour and not a defect: a plain `Vec<u32>` equality fails on
/// a mesh that round-tripped perfectly. Measured here: 0 of 3,200 triangles
/// differ up to rotation on a fixture where 7,500 of 9,600 indices differ
/// literally.
fn differing_triangles(a: &[u32], b: &[u32]) -> usize {
    assert_eq!(a.len(), b.len(), "decode returned a different index count");
    let (a3, _) = a.as_chunks::<3>();
    let (b3, _) = b.as_chunks::<3>();
    a3.iter()
        .zip(b3)
        .filter(|(x, y)| {
            !((x[0] == y[0] && x[1] == y[1] && x[2] == y[2])
                || (x[1] == y[0] && x[2] == y[1] && x[0] == y[2])
                || (x[2] == y[0] && x[0] == y[1] && x[1] == y[2]))
        })
        .count()
}

// ─── the clock, and the counters where they exist ───────────────────────────

/// One timed window's clock evidence.
#[derive(Clone, Copy)]
struct Window {
    /// Cycles over the whole window.
    cycles: f64,
    /// Cycles ÷ nanoseconds. `M-280`: on a governed CPU a nanosecond is not a
    /// unit, so the clock goes on the row.
    ghz: f64,
    /// The worst counter scheduling ratio, so a multiplexed count cannot be
    /// mistaken for a measurement.
    worst_ratio: f64,
}

/// Hardware counters where `perf_event_open` exists.
#[cfg(target_os = "linux")]
struct Meter(common::counters::Probe);

/// No counters. The columns say `unavailable` rather than inventing a number —
/// the precedent is `benches/family.rs`.
#[cfg(not(target_os = "linux"))]
struct Meter;

#[cfg(target_os = "linux")]
impl Meter {
    fn open() -> Self {
        Self(common::counters::Probe::open())
    }

    fn start(&mut self) {
        self.0.reset_and_enable();
    }

    fn finish(&mut self, total_ns: f64) -> Option<Window> {
        self.0.disable();
        let counts = self.0.read();
        let cycles = counts.cycles.count as f64;
        Some(Window {
            cycles,
            ghz: cycles / total_ns,
            worst_ratio: counts.worst_ratio(),
        })
    }
}

#[cfg(not(target_os = "linux"))]
impl Meter {
    fn open() -> Self {
        Self
    }

    fn start(&mut self) {}

    fn finish(&mut self, _total_ns: f64) -> Option<Window> {
        None
    }
}

/// A number, or the word that says it was not measured.
fn num(v: Option<f64>, precision: usize) -> String {
    match v {
        Some(x) => format!("{x:.*}", precision),
        None => String::from("unavailable"),
    }
}

// ─── one measured point ─────────────────────────────────────────────────────

/// What one knob setting produced.
struct Point {
    /// The knob's value, for the row.
    knob: f64,
    /// Triangles the extraction emitted.
    triangles: usize,
    /// Vertices it emitted.
    vertices: usize,
    /// Median extraction, milliseconds.
    ms: f64,
    /// Clock evidence over the rep set.
    window: Option<Window>,
    /// The mesh, kept for the encoder arm and the regeneration control.
    mesh_bytes: Vec<u8>,
    /// What the real encoder did to it.
    encoded: Encoded,
    /// The description's byte count.
    desc_bytes: usize,
    /// Brushes in the log.
    brushes: u32,
    /// Mesh bytes that differ after regenerating from the decoded description.
    regen_differing: usize,
    /// The same, after one bit of the encoded `cell_size` is flipped.
    regen_control_differing: usize,
}

/// One point's clock, kept as the best extraction seen anywhere.
#[derive(Clone, Copy)]
struct Timing {
    /// The fastest single extraction, in nanoseconds.
    ns: f64,
    /// Counters from the pass that produced it.
    window: Option<Window>,
}

impl Timing {
    /// A clock that nothing has beaten yet.
    const fn unset() -> Self {
        Self {
            ns: f64::MAX,
            window: None,
        }
    }

    /// Keep whichever of the two saw the faster extraction.
    fn best(self, other: Self) -> Self {
        if other.ns < self.ns { other } else { self }
    }
}

/// Time `REPS` extractions of one point and keep the fastest.
fn time_point<R, S>(
    meter: &mut Meter,
    mc: &mut MarchingCubes<R>,
    base: &S,
    brushes: &[Brush<Sphere<R>>],
    grid: &Grid<R>,
    out: &mut MeshBuffer<R>,
) -> Timing
where
    R: Real,
    S: Sdf<Scalar = R>,
{
    let mut fastest = f64::MAX;
    meter.start();
    let all = Instant::now();
    for _ in 0..REPS {
        let one = Instant::now();
        extract_with(mc, base, brushes, grid, out);
        fastest = fastest.min(one.elapsed().as_nanos() as f64);
        std::hint::black_box(out.positions.len());
    }
    let window = meter.finish(all.elapsed().as_nanos() as f64);
    Timing {
        ns: fastest,
        window,
    }
}

/// Encode one point's mesh, and regenerate it from its own description bytes.
///
/// Separated from the timing because it runs **once** per point while the timing
/// runs `PASSES` times, and because none of it belongs inside a timed window.
fn characterise<R, S>(
    mc: &mut MarchingCubes<R>,
    base: &S,
    brushes: &[Brush<Sphere<R>>],
    desc: &Desc<R>,
    rebuild: &dyn Fn(&Desc<R>) -> S,
    timing: (f64, Timing),
) -> Point
where
    R: Wire,
    S: Sdf<Scalar = R>,
{
    let (knob, clock) = timing;
    let grid = grid_of(desc);
    let mut out = MeshBuffer::<R>::new();
    extract_with(mc, base, brushes, &grid, &mut out);

    let triangles = out.indices.len() / 3;
    assert!(
        triangles > 0,
        "VOID: knob {knob} emitted no triangle, so its time measures an empty surface"
    );
    let original = mesh_bytes(&out);

    let lo = [
        desc.origin[0].as_f32(),
        desc.origin[1].as_f32(),
        desc.origin[2].as_f32(),
    ];
    let extent = desc.cell_size.as_f32() * desc.cells as f32;
    let encoded = encode_and_time(&out, lo, extent);

    // ── the description is a byte string, and it is sufficient ──────────────
    let (bytes, cell_size_high_byte) = encode_desc(desc);
    assert_eq!(
        bytes.len(),
        desc.encoded_len(),
        "the encoder and its own length function disagree"
    );
    let back: Desc<R> = decode_desc(&bytes);
    assert_eq!(&back, desc, "the description did not survive its own codec");
    let regen_brushes: Vec<Brush<Sphere<R>>> =
        back.brushes.iter().copied().map(Brush::subtract).collect();
    let regen_base = rebuild(&back);
    let mut regen = MeshBuffer::<R>::new();
    extract_with(mc, &regen_base, &regen_brushes, &grid_of(&back), &mut regen);
    let regen_differing = bytes_differing(&original, &mesh_bytes(&regen));

    // ── the control: one bit of the encoded cell_size's exponent ────────────
    let mut broken = bytes.clone();
    broken[cell_size_high_byte] ^= 0x01;
    let bad: Desc<R> = decode_desc(&broken);
    let bad_brushes: Vec<Brush<Sphere<R>>> =
        bad.brushes.iter().copied().map(Brush::subtract).collect();
    let bad_base = rebuild(&bad);
    let mut bad_mesh = MeshBuffer::<R>::new();
    extract_with(mc, &bad_base, &bad_brushes, &grid_of(&bad), &mut bad_mesh);
    let regen_control_differing = bytes_differing(&original, &mesh_bytes(&bad_mesh));

    Point {
        knob,
        triangles,
        vertices: out.positions.len(),
        ms: clock.ns / 1.0e6,
        window: clock.window,
        mesh_bytes: original,
        encoded,
        desc_bytes: bytes.len(),
        brushes: desc.brushes.len() as u32,
        regen_differing,
        regen_control_differing,
    }
}

// ─── the two-term fit ───────────────────────────────────────────────────────

/// Least squares on `y = a + b·x`, with the coefficient of determination.
fn fit(points: &[(f64, f64)]) -> (f64, f64, f64) {
    let n = points.len() as f64;
    let sx: f64 = points.iter().map(|p| p.0).sum();
    let sy: f64 = points.iter().map(|p| p.1).sum();
    let sxx: f64 = points.iter().map(|p| p.0 * p.0).sum();
    let sxy: f64 = points.iter().map(|p| p.0 * p.1).sum();
    let b = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    let a = (sy - b * sx) / n;
    let mean = sy / n;
    let ss_tot: f64 = points.iter().map(|p| (p.1 - mean) * (p.1 - mean)).sum();
    let ss_res: f64 = points
        .iter()
        .map(|p| {
            let r = p.1 - (a + b * p.0);
            r * r
        })
        .sum();
    (a, b, 1.0 - ss_res / ss_tot)
}

// ─── deterministic scatter ──────────────────────────────────────────────────

/// A 64-bit LCG, so the brush path is the same on both machines and in both
/// runs. Numerical Recipes' constants.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        // The top 53 bits, so the value is exactly representable and identical
        // on any machine with IEEE doubles.
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

/// Brush centres scattered strictly inside one chunk.
///
/// The inset is `radius + one cell`, so a brush can never reach a neighbouring
/// chunk. That is what makes a chunk's log **exactly** its own brushes in the
/// replay, with no windowing approximation to get wrong.
fn scatter<R: Real>(origin: [R; 3], radius: f64, count: u32, seed: u64) -> Vec<Sphere<R>> {
    let inset = radius + CELL_SIZE;
    let span = CHUNK_EXTENT - 2.0 * inset;
    assert!(span > 0.0, "the brush does not fit inside a chunk");
    let mut lcg = Lcg(seed);
    (0..count)
        .map(|_| {
            let mut c = [R::ZERO; 3];
            for (slot, o) in c.iter_mut().zip(origin) {
                *slot = o + R::from_f64(inset + span * lcg.next_unit());
            }
            Sphere {
                center: c,
                radius: R::from_f64(radius),
            }
        })
        .collect()
}

/// `M-50`'s bucket a log of `n` brushes falls in.
fn bucket(n: u32) -> &'static str {
    match n {
        1..=15 => "1-15",
        16..=30 => "16-30",
        31..=45 => "31-45",
        46..=60 => "46-60",
        _ => panic!("{n} brushes is outside M-50's four buckets"),
    }
}

// ─── the row ────────────────────────────────────────────────────────────────

/// Every column, so that no row is ragged.
struct Emit<'a> {
    arm: &'a str,
    field: &'a str,
    precision: &'a str,
    knob: f64,
    point: &'a Point,
    ns_per_triangle_marginal: f64,
    cycles_per_triangle_marginal: Option<f64>,
    fit_a_ms: f64,
    fit_r2: f64,
    fit_points: usize,
    fit_sound: &'a str,
    brushes_biting: String,
    brushes_drawn: String,
    replay_edits: usize,
    replay_min_triangles: String,
    replay_total_triangles: String,
    replay_ms: String,
    bytes_differing: usize,
    bytes_differing_scope: &'a str,
    peer: &'a str,
    c1: &'a str,
    c2: &'a str,
    c3: &'a str,
}

type Row = Vec<(&'static str, String)>;

impl Emit<'_> {
    fn row(&self) -> Row {
        let p = self.point;
        let tri = p.triangles as f64;
        let meshopt_bytes = p.encoded.vertex_bytes + p.encoded.index_bytes;
        let size_ratio = meshopt_bytes as f64 / p.desc_bytes as f64;
        let total_ns = p.ms * 1.0e6;
        let marginal_share = self.ns_per_triangle_marginal * tri / total_ns;
        vec![
            ("field", self.field.to_string()),
            ("chunk_cells", CHUNK_CELLS.to_string()),
            ("triangles", p.triangles.to_string()),
            ("extract_ms", format!("{:.6}", p.ms)),
            (
                "ns_per_triangle",
                format!("{:.4}", self.ns_per_triangle_marginal),
            ),
            ("meshopt_bytes", meshopt_bytes.to_string()),
            ("meshopt_decode_ms", format!("{:.6}", p.encoded.decode_ms)),
            (
                "meshopt_ns_per_triangle",
                format!("{:.4}", p.encoded.decode_ms * 1.0e6 / tri),
            ),
            ("field_plus_log_bytes", p.desc_bytes.to_string()),
            ("size_ratio", format!("{size_ratio:.4}")),
            ("log_bucket", bucket_or_none(p.brushes).to_string()),
            ("replay_edits", self.replay_edits.to_string()),
            ("bytes_differing", self.bytes_differing.to_string()),
            ("c1_holds", self.c1.to_string()),
            ("c2_holds", self.c2.to_string()),
            ("c3_holds", self.c3.to_string()),
            // ── extras ──────────────────────────────────────────────────────
            ("arm", self.arm.to_string()),
            ("precision", self.precision.to_string()),
            ("knob", format!("{:.4}", self.knob)),
            ("samples", (CHUNK_SAMPLES.pow(3)).to_string()),
            ("vertices", p.vertices.to_string()),
            ("brushes", p.brushes.to_string()),
            ("brushes_biting", self.brushes_biting.clone()),
            ("brushes_drawn", self.brushes_drawn.clone()),
            ("ns_per_triangle_total", format!("{:.4}", total_ns / tri)),
            (
                "cycles_per_triangle",
                num(self.cycles_per_triangle_marginal, 4),
            ),
            (
                "cycles_per_triangle_total",
                num(self.point.window.map(|w| w.cycles / (tri * REPS as f64)), 4),
            ),
            ("ghz", num(p.window.map(|w| w.ghz), 4)),
            ("worst_ratio", num(p.window.map(|w| w.worst_ratio), 4)),
            ("fit_a_ms", format!("{:.6}", self.fit_a_ms)),
            ("fit_r2", format!("{:.6}", self.fit_r2)),
            ("fit_points", self.fit_points.to_string()),
            ("fit_sound", self.fit_sound.to_string()),
            ("triangle_term_share", format!("{marginal_share:.4}")),
            (
                "published_ns_per_triangle",
                format!("{PUBLISHED_DECODE_NS_PER_TRIANGLE:.1}"),
            ),
            (
                "c1_holds_vs_measured",
                (self.ns_per_triangle_marginal < p.encoded.decode_ms * 1.0e6 / tri).to_string(),
            ),
            ("meshopt_vertex_bytes", p.encoded.vertex_bytes.to_string()),
            ("meshopt_index_bytes", p.encoded.index_bytes.to_string()),
            ("meshopt_bytes_f32", p.encoded.f32_bytes.to_string()),
            (
                "meshopt_bytes_per_triangle",
                format!("{:.4}", meshopt_bytes as f64 / tri),
            ),
            (
                "meshopt_roundtrip_differing",
                p.encoded.roundtrip_differing.to_string(),
            ),
            (
                "meshopt_control_differing",
                p.encoded.control_differing.to_string(),
            ),
            (
                "size_ratio_f32",
                format!("{:.4}", p.encoded.f32_bytes as f64 / p.desc_bytes as f64),
            ),
            ("mesh_bytes", p.mesh_bytes.len().to_string()),
            ("regen_differing", p.regen_differing.to_string()),
            (
                "regen_control_differing",
                p.regen_control_differing.to_string(),
            ),
            (
                "bytes_differing_scope",
                self.bytes_differing_scope.to_string(),
            ),
            ("peer", self.peer.to_string()),
            ("replay_min_triangles", self.replay_min_triangles.clone()),
            (
                "replay_total_triangles",
                self.replay_total_triangles.clone(),
            ),
            ("replay_ms", self.replay_ms.clone()),
        ]
    }
}

/// The bucket, or `none` for a row with no log at all.
fn bucket_or_none(n: u32) -> &'static str {
    if n == 0 { "none" } else { bucket(n) }
}

// ─── C1: one family, swept ──────────────────────────────────────────────────

/// The coefficient of determination below which the two-term model is not
/// describing the data and its slope is not a marginal cost.
const FIT_R2_FLOOR: f64 = 0.95;

/// A family's fitted marginal and the points behind it.
struct Family {
    points: Vec<Point>,
    a_ms: f64,
    ns_per_triangle: f64,
    cycles_per_triangle: Option<f64>,
    r2: f64,
    /// `a > 0` **and** `r² ≥ FIT_R2_FLOOR`. False means the slope beside it is
    /// not a marginal cost and must not be quoted as one (`M-21`).
    sound: bool,
}

/// Sweep one family's knob at fixed 33³, then fit `t = a + b·T`.
///
/// **The knob points are timed round-robin, `PASSES` times over.** Timing every
/// rep of one point before moving to the next makes a burst of interference
/// land entirely on that point, which is a slope error rather than a level
/// error — and this machine is not quiet, because a dozen sibling agents run
/// `cargo` on it. Measured: the same `box_exact/f64` fixture gave `r² = 0.99994`
/// on one run and **`0.32081`** on the next with all fifteen reps of a point
/// taken together. Round-robin spreads each point's samples across the whole
/// sweep, so load that comes and goes moves every point rather than one.
fn sweep<R, S>(
    meter: &mut Meter,
    mc: &mut MarchingCubes<R>,
    knobs: &[f64],
    build: &dyn Fn(f64) -> (S, Desc<R>),
    rebuild: &dyn Fn(&Desc<R>) -> S,
    brushes: &[Brush<Sphere<R>>],
) -> Family
where
    R: Wire,
    S: Sdf<Scalar = R>,
{
    // One shared output buffer for the whole family, warmed on every knob
    // before anything is timed: the first-touch page faults and the capacity
    // growth belong to no point's arithmetic.
    let mut out = MeshBuffer::<R>::new();
    let mut grids = Vec::with_capacity(knobs.len());
    for &k in knobs {
        let (base, desc) = build(k);
        let grid = grid_of(&desc);
        extract_with(mc, &base, brushes, &grid, &mut out);
        grids.push(grid);
    }

    let mut clocks = vec![Timing::unset(); knobs.len()];
    for _pass in 0..PASSES {
        for (i, &k) in knobs.iter().enumerate() {
            let (base, _) = build(k);
            let seen = time_point(meter, mc, &base, brushes, &grids[i], &mut out);
            clocks[i] = clocks[i].best(seen);
        }
    }

    let mut points = Vec::with_capacity(knobs.len());
    for (i, &k) in knobs.iter().enumerate() {
        let (base, desc) = build(k);
        points.push(characterise(
            mc,
            &base,
            brushes,
            &desc,
            rebuild,
            (k, clocks[i]),
        ));
    }
    let ms: Vec<(f64, f64)> = points.iter().map(|p| (p.triangles as f64, p.ms)).collect();
    let (a_ms, b_ms, r2) = fit(&ms);
    let cycles_per_triangle = if points.iter().all(|p| p.window.is_some()) {
        let cy: Vec<(f64, f64)> = points
            .iter()
            .map(|p| {
                (
                    p.triangles as f64,
                    p.window.expect("checked").cycles / REPS as f64,
                )
            })
            .collect();
        Some(fit(&cy).1)
    } else {
        None
    };
    Family {
        points,
        a_ms,
        ns_per_triangle: b_ms * 1.0e6,
        cycles_per_triangle,
        r2,
        sound: a_ms > 0.0 && r2 >= FIT_R2_FLOOR,
    }
}

// ─── C3: the replay ─────────────────────────────────────────────────────────

/// What the 10⁴-edit replay produced.
struct Replay {
    /// One 8-byte `mesh_hash` per edit, in order.
    digests: Vec<u8>,
    /// The smallest triangle count over all 10,000 extractions.
    min_triangles: usize,
    /// Triangles summed over the replay.
    total_triangles: usize,
    /// Wall time of the whole replay, milliseconds.
    ms: f64,
    /// The last chunk's final description, so the row's size columns are real.
    final_desc: Desc<f64>,
    /// The last chunk's final mesh.
    final_mesh: MeshBuffer<f64>,
}

/// Replay 10,000 edits over 200 chunks and hash every intermediate mesh.
///
/// `mesh_hash` is `M-31`'s own instrument — the function the 216 committed golden
/// hashes are taken with — so C3 is that gate extended to a long trace rather
/// than a second opinion about it.
fn replay() -> Replay {
    let layout = ChunkLayout::<f64>::new(CHUNK_CELLS, CELL_SIZE, [0.0; 3]).expect("layout");
    let base = Gyroid::<f64> {
        scale: REPLAY_GYROID_SCALE,
        iso: 0.0,
    };
    let radius = BRUSH_CELLS * CELL_SIZE;

    let mut mc = MarchingCubes::<f64>::new();
    let mut digests = Vec::with_capacity(REPLAY_EDITS * 8);
    let mut min_triangles = usize::MAX;
    let mut total_triangles = 0usize;
    let mut out = MeshBuffer::<f64>::new();
    let mut final_desc = None;

    let start = Instant::now();
    for c in 0..REPLAY_CHUNKS {
        // A 10 x 10 x 2 lattice of chunks, so the trace crosses 200 distinct
        // origins and every one of them has its own log.
        let id = ChunkId::new([(c % 10) as i32, ((c / 10) % 10) as i32, (c / 100) as i32]);
        let origin = layout.sample_origin(id);
        let spheres = scatter::<f64>(
            origin,
            radius,
            REPLAY_EDITS_PER_CHUNK as u32,
            0x5EED_0000 ^ c as u64,
        );
        let brushes: Vec<Brush<Sphere<f64>>> =
            spheres.iter().copied().map(Brush::subtract).collect();
        let grid = Grid {
            shape: layout.sample_shape().expect("shape"),
            origin,
            cell_size: CELL_SIZE,
        };
        for k in 0..REPLAY_EDITS_PER_CHUNK {
            extract_with(&mut mc, &base, &brushes[..=k], &grid, &mut out);
            let triangles = out.indices.len() / 3;
            min_triangles = min_triangles.min(triangles);
            total_triangles += triangles;
            digests.extend_from_slice(&mesh_hash(&out).to_le_bytes());
        }
        if c == REPLAY_CHUNKS - 1 {
            final_desc = Some(Desc {
                tag: 1,
                seed: 0,
                octaves: 0,
                params: vec![base.scale, base.iso],
                cells: CHUNK_CELLS,
                origin,
                cell_size: CELL_SIZE,
                brushes: spheres,
            });
        }
    }
    let ms = start.elapsed().as_nanos() as f64 / 1.0e6;

    let mut final_mesh = MeshBuffer::<f64>::new();
    std::mem::swap(&mut final_mesh, &mut out);
    Replay {
        digests,
        min_triangles,
        total_triangles,
        ms,
        final_desc: final_desc.expect("the last chunk ran"),
        final_mesh,
    }
}

/// The workspace root, as `benches/common/experiment.rs` computes it.
fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// This machine's slug, from the one script that decides what a machine is
/// called.
fn machine_slug() -> String {
    let out = std::process::Command::new(repo_root().join("scripts/machine.sh"))
        .arg("--slug")
        .output()
        .expect("scripts/machine.sh runs");
    String::from_utf8(out.stdout)
        .expect("the slug is utf-8")
        .trim()
        .to_string()
}

// ─── main ───────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let mut rows: Vec<Row> = Vec::new();
    let mut meter = Meter::open();

    // ── the SHARE line, recomputed and printed before anything is measured ──
    let samples = f64::from(CHUNK_SAMPLES.pow(3));
    let cells = f64::from(CHUNK_CELLS.pow(3));
    for (label, per_sample) in [
        ("x51 extraction marginal", 10.68),
        ("committed sweep", 13.1892),
    ] {
        let total_ns = samples * per_sample;
        println!(
            "SHARE  {label}: {per_sample} ns/sample x {samples:.0} samples = {:.1} us per 33^3 \
             chunk, so the TOTAL reading of C1 needs {:.0} triangles from {cells:.0} cells = \
             {:.2} triangles per cell",
            total_ns / 1.0e3,
            total_ns / PUBLISHED_DECODE_NS_PER_TRIANGLE,
            total_ns / PUBLISHED_DECODE_NS_PER_TRIANGLE / cells
        );
    }
    println!(
        "SHARE  the marginal reading is what the registration's wording asks for and what is \
         fitted below; triangle_term_share reports how much of a chunk extraction it covers\n"
    );

    // ── C1 ─────────────────────────────────────────────────────────────────
    let origin32 = [-(CHUNK_EXTENT as f32) * 0.5; 3];
    let origin64 = [-CHUNK_EXTENT * 0.5; 3];
    let cell32 = CELL_SIZE as f32;

    // Area knobs, all of them. The half-extents and radii are deliberately not
    // multiples of `CELL_SIZE`, so no face or pole lands exactly on a sample
    // plane and the triangle count is not decided by a tie-break.
    let box_knobs = [0.31, 0.53, 0.77, 1.01, 1.23, 1.47, 1.69, 1.83];
    let sphere_knobs = [0.21, 0.43, 0.65, 0.87, 1.09, 1.31, 1.53, 1.81];
    // Gyroid `scale`: the surface area inside a fixed box grows with it, and
    // scaling the field is scaling the surface rather than changing what kind of
    // surface it is.
    let gyroid_knobs = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
    // fbm `amplitude` at the canonical frequency: a steeper sheet is a larger
    // sheet, and it is still one sheet.
    let fbm_knobs = [0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5];

    let dug32 = scatter::<f32>(
        origin32,
        BRUSH_CELLS * CELL_SIZE,
        C1_DUG_BRUSHES,
        0x0092_0001,
    );
    let dug64 = scatter::<f64>(
        origin64,
        BRUSH_CELLS * CELL_SIZE,
        C1_DUG_BRUSHES,
        0x0092_0001,
    );
    let dug32_brushes: Vec<Brush<Sphere<f32>>> =
        dug32.iter().copied().map(Brush::subtract).collect();
    let dug64_brushes: Vec<Brush<Sphere<f64>>> =
        dug64.iter().copied().map(Brush::subtract).collect();

    let mut mc32 = MarchingCubes::<f32>::new();
    let mut mc64 = MarchingCubes::<f64>::new();

    let mut families: Vec<(&'static str, &'static str, Family)> = Vec::new();

    // `box_exact` first, and it is the family that carries the verdict. `✗51`
    // classified the crate's field bodies by what `libm` does with them and this
    // one has "no calls at all" -- so its per-triangle marginal is a **lower
    // bound** on any field's, and a falsification here cannot be blamed on an
    // expensive field.
    families.push((
        "box_exact",
        "f32",
        sweep::<f32, BoxExact<f32>>(
            &mut meter,
            &mut mc32,
            &box_knobs,
            &|e| {
                (
                    BoxExact {
                        center: [0.0; 3],
                        half_extents: [e as f32; 3],
                    },
                    Desc {
                        tag: 3,
                        seed: 0,
                        octaves: 0,
                        params: vec![0.0, 0.0, 0.0, e as f32, e as f32, e as f32],
                        cells: CHUNK_CELLS,
                        origin: origin32,
                        cell_size: cell32,
                        brushes: Vec::new(),
                    },
                )
            },
            &|d| BoxExact {
                center: [d.params[0], d.params[1], d.params[2]],
                half_extents: [d.params[3], d.params[4], d.params[5]],
            },
            &[],
        ),
    ));
    families.push((
        "box_exact",
        "f64",
        sweep::<f64, BoxExact<f64>>(
            &mut meter,
            &mut mc64,
            &box_knobs,
            &|e| {
                (
                    BoxExact {
                        center: [0.0; 3],
                        half_extents: [e; 3],
                    },
                    Desc {
                        tag: 3,
                        seed: 0,
                        octaves: 0,
                        params: vec![0.0, 0.0, 0.0, e, e, e],
                        cells: CHUNK_CELLS,
                        origin: origin64,
                        cell_size: CELL_SIZE,
                        brushes: Vec::new(),
                    },
                )
            },
            &|d| BoxExact {
                center: [d.params[0], d.params[1], d.params[2]],
                half_extents: [d.params[3], d.params[4], d.params[5]],
            },
            &[],
        ),
    ));
    families.push((
        "sphere",
        "f32",
        sweep::<f32, Sphere<f32>>(
            &mut meter,
            &mut mc32,
            &sphere_knobs,
            &|r| {
                (
                    Sphere {
                        center: [0.0; 3],
                        radius: r as f32,
                    },
                    Desc {
                        tag: 0,
                        seed: 0,
                        octaves: 0,
                        params: vec![0.0, 0.0, 0.0, r as f32],
                        cells: CHUNK_CELLS,
                        origin: origin32,
                        cell_size: cell32,
                        brushes: Vec::new(),
                    },
                )
            },
            &|d| Sphere {
                center: [d.params[0], d.params[1], d.params[2]],
                radius: d.params[3],
            },
            &[],
        ),
    ));
    families.push((
        "sphere",
        "f64",
        sweep::<f64, Sphere<f64>>(
            &mut meter,
            &mut mc64,
            &sphere_knobs,
            &|r| {
                (
                    Sphere {
                        center: [0.0; 3],
                        radius: r,
                    },
                    Desc {
                        tag: 0,
                        seed: 0,
                        octaves: 0,
                        params: vec![0.0, 0.0, 0.0, r],
                        cells: CHUNK_CELLS,
                        origin: origin64,
                        cell_size: CELL_SIZE,
                        brushes: Vec::new(),
                    },
                )
            },
            &|d| Sphere {
                center: [d.params[0], d.params[1], d.params[2]],
                radius: d.params[3],
            },
            &[],
        ),
    ));
    families.push((
        "gyroid",
        "f32",
        sweep::<f32, Gyroid<f32>>(
            &mut meter,
            &mut mc32,
            &gyroid_knobs,
            &|scale| {
                (
                    Gyroid {
                        scale: scale as f32,
                        iso: 0.0,
                    },
                    Desc {
                        tag: 1,
                        seed: 0,
                        octaves: 0,
                        params: vec![scale as f32, 0.0f32],
                        cells: CHUNK_CELLS,
                        origin: origin32,
                        cell_size: cell32,
                        brushes: Vec::new(),
                    },
                )
            },
            &|d| Gyroid {
                scale: d.params[0],
                iso: d.params[1],
            },
            &[],
        ),
    ));
    families.push((
        "gyroid",
        "f64",
        sweep::<f64, Gyroid<f64>>(
            &mut meter,
            &mut mc64,
            &gyroid_knobs,
            &|scale| {
                (
                    Gyroid { scale, iso: 0.0 },
                    Desc {
                        tag: 1,
                        seed: 0,
                        octaves: 0,
                        params: vec![scale, 0.0f64],
                        cells: CHUNK_CELLS,
                        origin: origin64,
                        cell_size: CELL_SIZE,
                        brushes: Vec::new(),
                    },
                )
            },
            &|d| Gyroid {
                scale: d.params[0],
                iso: d.params[1],
            },
            &[],
        ),
    ));
    families.push((
        "fbm_terrain",
        "f32",
        sweep::<f32, FbmTerrain<f32>>(
            &mut meter,
            &mut mc32,
            &fbm_knobs,
            &|f| {
                let mut field = FbmTerrain::<f32>::canonical();
                field.amplitude = f as f32;
                (field, fbm_desc(&field, origin32, cell32))
            },
            &|d| fbm_from(d),
            &[],
        ),
    ));
    families.push((
        "fbm_terrain",
        "f64",
        sweep::<f64, FbmTerrain<f64>>(
            &mut meter,
            &mut mc64,
            &fbm_knobs,
            &|f| {
                let mut field = FbmTerrain::<f64>::canonical();
                field.amplitude = f;
                (field, fbm_desc(&field, origin64, CELL_SIZE))
            },
            &|d| fbm_from(d),
            &[],
        ),
    ));
    families.push((
        "gyroid_dug30",
        "f32",
        sweep::<f32, Gyroid<f32>>(
            &mut meter,
            &mut mc32,
            &gyroid_knobs,
            &|scale| {
                (
                    Gyroid {
                        scale: scale as f32,
                        iso: 0.0,
                    },
                    Desc {
                        tag: 1,
                        seed: 0,
                        octaves: 0,
                        params: vec![scale as f32, 0.0f32],
                        cells: CHUNK_CELLS,
                        origin: origin32,
                        cell_size: cell32,
                        brushes: dug32.clone(),
                    },
                )
            },
            &|d| Gyroid {
                scale: d.params[0],
                iso: d.params[1],
            },
            &dug32_brushes,
        ),
    ));
    families.push((
        "gyroid_dug30",
        "f64",
        sweep::<f64, Gyroid<f64>>(
            &mut meter,
            &mut mc64,
            &gyroid_knobs,
            &|scale| {
                (
                    Gyroid { scale, iso: 0.0 },
                    Desc {
                        tag: 1,
                        seed: 0,
                        octaves: 0,
                        params: vec![scale, 0.0f64],
                        cells: CHUNK_CELLS,
                        origin: origin64,
                        cell_size: CELL_SIZE,
                        brushes: dug64.clone(),
                    },
                )
            },
            &|d| Gyroid {
                scale: d.params[0],
                iso: d.params[1],
            },
            &dug64_brushes,
        ),
    ));

    println!(
        "{:>13} {:>4} {:>8} {:>9} {:>10} {:>10} {:>9} {:>8} {:>6} {:>9} {:>10}",
        "family",
        "prec",
        "tri lo",
        "tri hi",
        "a (ms)",
        "ns/tri",
        "cyc/tri",
        "r2",
        "sound",
        "share",
        "verdict"
    );
    let mut sound_families = 0usize;
    let mut cheap_arms = 0usize;
    let mut cheap_sound = 0usize;
    for (name, precision, family) in &families {
        let lo = family
            .points
            .iter()
            .map(|p| p.triangles)
            .min()
            .expect("swept");
        let hi = family
            .points
            .iter()
            .map(|p| p.triangles)
            .max()
            .expect("swept");
        let ms_lo = family.points.iter().map(|p| p.ms).fold(f64::MAX, f64::min);
        let ms_hi = family.points.iter().map(|p| p.ms).fold(0.0f64, f64::max);

        // ── M-21: a negative fixed cost is the model saying it is wrong ─────
        //
        // Recorded rather than fatal per family, because *which* families the
        // two-term model describes is itself a result — and because this machine
        // is shared, so an individual arm can lose its fit to a neighbour's
        // compile. What is fatal is the *group*: at least one of the four cheap
        // arms (`box_exact` and `sphere`, in both precisions) has to have a sound
        // fit, because the verdict rests on a lower bound and without one there
        // is no lower bound to rest on. Asserted after the loop.
        if *name == "box_exact" || *name == "sphere" {
            cheap_arms += 1;
            if family.sound {
                cheap_sound += 1;
            }
        }
        if family.sound {
            sound_families += 1;
        }
        // ── M-19: the coefficient against both ends of its own data ─────────
        println!(
            "  {name}/{precision}: a = {:.6} ms = {:.1}% of the largest run ({ms_hi:.6} ms) and \
             {:.1}% of the smallest ({ms_lo:.6} ms)",
            family.a_ms,
            100.0 * family.a_ms / ms_hi,
            100.0 * family.a_ms / ms_lo
        );

        let holds = family.ns_per_triangle < PUBLISHED_DECODE_NS_PER_TRIANGLE;
        let share_at_hi = family
            .points
            .iter()
            .map(|p| family.ns_per_triangle * p.triangles as f64 / (p.ms * 1.0e6))
            .fold(0.0f64, f64::max);
        println!(
            "{name:>13} {precision:>4} {lo:>8} {hi:>9} {:>10.6} {:>10.4} {:>9} {:>8.5} \
             {:>6} {:>9.4} {:>10}",
            family.a_ms,
            family.ns_per_triangle,
            num(family.cycles_per_triangle, 2),
            family.r2,
            family.sound,
            share_at_hi,
            if !family.sound {
                "MODEL-BAD"
            } else if holds {
                "HELD"
            } else {
                "FALSIFIED"
            }
        );

        for point in &family.points {
            assert_eq!(
                point.encoded.roundtrip_differing, 0,
                "{name}/{precision} knob {}: the real encoder lost {} vertices-or-triangles, so \
                 the encoded size is not an encoding of these triangles",
                point.knob, point.encoded.roundtrip_differing
            );
            assert_eq!(
                point.encoded.control_differing, 1,
                "{name}/{precision} knob {}: moving one quantised component by one ULP changed \
                 {} vertices, so the round-trip comparator is not the instrument it claims (M-44)",
                point.knob, point.encoded.control_differing
            );
            assert_eq!(
                point.regen_differing, 0,
                "{name}/{precision} knob {}: regenerating from the field-plus-log bytes moved {} \
                 mesh bytes, so the description is not sufficient and C2 is comparing sizes of \
                 things that are not equivalent",
                point.knob, point.regen_differing
            );
            assert!(
                point.regen_control_differing > 0,
                "{name}/{precision} knob {}: flipping one bit of the encoded cell_size changed \
                 nothing, so the regeneration is not reading the bytes (M-44)",
                point.knob
            );
            rows.push(
                Emit {
                    arm: "c1_marginal",
                    field: name,
                    precision,
                    knob: point.knob,
                    point,
                    ns_per_triangle_marginal: family.ns_per_triangle,
                    cycles_per_triangle_marginal: family.cycles_per_triangle,
                    fit_a_ms: family.a_ms,
                    fit_r2: family.r2,
                    fit_points: family.points.len(),
                    fit_sound: if family.sound { "true" } else { "false" },
                    brushes_biting: String::from("n/a"),
                    brushes_drawn: String::from("n/a"),
                    replay_edits: 0,
                    replay_min_triangles: String::from("n/a"),
                    replay_total_triangles: String::from("n/a"),
                    replay_ms: String::from("n/a"),
                    bytes_differing: point.regen_differing,
                    bytes_differing_scope: "regen-from-description",
                    peer: "n/a",
                    c1: if holds { "true" } else { "false" },
                    c2: "n/a",
                    c3: "n/a",
                }
                .row(),
            );
        }
    }
    assert!(
        cheap_sound > 0,
        "none of the {cheap_arms} cheap arms (box_exact, sphere) has a two-term model that \
         describes it, so there is no lower bound on the per-triangle cost and C1 cannot be \
         scored (M-21)"
    );
    // The count below is reported, not asserted. On a machine a dozen sibling
    // agents are compiling on, *which* families keep a sound fit moves run to
    // run — measured: the same `gyroid/f32` fixture gave `r² = 0.98488` on a
    // quiet minute and `0.76258` on a busy one. The cheap arms are what the
    // verdict needs and they are asserted; the rest is reported with `fit_sound`
    // beside it so nobody reads a slope that has no model behind it.
    println!(
        "\n{sound_families} of {} families have a sound two-term fit (a > 0 and r2 >= \
         {FIT_R2_FLOOR}), including {cheap_sound} of the {cheap_arms} cheap arms that carry the \
         lower bound",
        families.len()
    );

    // ── the encoder is not an estimate: bytes/triangle has to move ──────────
    let bpt: Vec<f64> = families
        .iter()
        .flat_map(|(_, _, f)| f.points.iter())
        .map(|p| (p.encoded.vertex_bytes + p.encoded.index_bytes) as f64 / p.triangles as f64)
        .collect();
    let bpt_lo = bpt.iter().copied().fold(f64::MAX, f64::min);
    let bpt_hi = bpt.iter().copied().fold(0.0f64, f64::max);
    assert!(
        bpt_hi - bpt_lo > 0.0,
        "meshopt_bytes_per_triangle is identical on all {} rows, which is what an estimate looks \
         like rather than an encoder",
        bpt.len()
    );
    println!(
        "\nencoder is real: bytes/triangle spans {bpt_lo:.4} to {bpt_hi:.4} over {} meshes",
        bpt.len()
    );

    // ── C2: the four log buckets ───────────────────────────────────────────
    //
    // **One 60-brush log, and its four prefixes are the buckets.** The first
    // version drew each bucket independently from the LCG and the biting
    // control fired on the third: *"only 44 of 45 brushes changed the mesh"*.
    // That is a real property of a dig -- `P-94`'s C2 is the ticket for exactly
    // this collapse -- but a swallowed brush in *this* log is 17 bytes of
    // description with no geometry behind it, which is a fixture that answers a
    // slightly different question than the one C2 asks. So the log is built by
    // rejection: a candidate is kept only if adding it changed the mesh, and
    // `brushes_drawn` records how many candidates it took, which is a measurement
    // of the collapse rather than a defect hidden by a looser gate.
    //
    // Prefixes rather than four independent logs, because "brush `k` bit when
    // added to brushes `0..k`" is exactly the prefix condition, so every bucket
    // inherits the property instead of re-earning it.
    let base32 = Gyroid::<f32> {
        scale: REPLAY_GYROID_SCALE as f32,
        iso: 0.0,
    };
    let c2_grid = Grid {
        shape: RuntimeShape3::new([CHUNK_SAMPLES; 3]).expect("chunk grid fits u32"),
        origin: origin32,
        cell_size: cell32,
    };
    let biting_radius = BRUSH_CELLS * CELL_SIZE;
    let inset = biting_radius + CELL_SIZE;
    let span = CHUNK_EXTENT - 2.0 * inset;
    let mut lcg = Lcg(0x0092_0002);
    let mut log: Vec<Sphere<f32>> = Vec::new();
    let mut brushes_drawn = 0u32;
    let mut kept = MeshBuffer::<f32>::new();
    let mut candidate_mesh = MeshBuffer::<f32>::new();
    extract_with(&mut mc32, &base32, &[], &c2_grid, &mut kept);
    while log.len() < *LOG_BUCKETS.last().expect("four buckets") as usize {
        let mut c = [0.0f32; 3];
        for (slot, o) in c.iter_mut().zip(origin32) {
            *slot = o + (inset + span * lcg.next_unit()) as f32;
        }
        log.push(Sphere {
            center: c,
            radius: biting_radius as f32,
        });
        brushes_drawn += 1;
        let trial: Vec<Brush<Sphere<f32>>> = log.iter().copied().map(Brush::subtract).collect();
        extract_with(&mut mc32, &base32, &trial, &c2_grid, &mut candidate_mesh);
        if bytes_differing(&mesh_bytes(&kept), &mesh_bytes(&candidate_mesh)) > 0 {
            std::mem::swap(&mut kept, &mut candidate_mesh);
        } else {
            log.pop();
        }
    }
    println!(
        "\nC2 log: {} biting brushes from {brushes_drawn} candidates -- {} swallowed on arrival",
        log.len(),
        brushes_drawn - log.len() as u32
    );

    println!(
        "\n{:>10} {:>8} {:>9} {:>8} {:>12} {:>10} {:>10} {:>9}",
        "bucket", "brushes", "biting", "tri", "meshopt B", "desc B", "ratio", "verdict"
    );
    let mut c2_rows: Vec<Row> = Vec::new();
    let mut c2_all_hold = true;
    let mut c2_ratios: Vec<(u32, f64)> = Vec::new();
    for &n in &LOG_BUCKETS {
        let spheres: Vec<Sphere<f32>> = log[..n as usize].to_vec();
        let brushes: Vec<Brush<Sphere<f32>>> =
            spheres.iter().copied().map(Brush::subtract).collect();
        let base = base32;

        let desc = Desc {
            tag: 1,
            seed: 0,
            octaves: 0,
            params: vec![REPLAY_GYROID_SCALE as f32, 0.0f32],
            cells: CHUNK_CELLS,
            origin: origin32,
            cell_size: cell32,
            brushes: spheres,
        };

        // ── control: the prefix property, verified rather than inherited ────
        let mut biting = 0u32;
        let mut before = MeshBuffer::<f32>::new();
        let mut after = MeshBuffer::<f32>::new();
        for k in 0..n as usize {
            extract_with(&mut mc32, &base, &brushes[..k], &c2_grid, &mut before);
            extract_with(&mut mc32, &base, &brushes[..=k], &c2_grid, &mut after);
            if bytes_differing(&mesh_bytes(&before), &mesh_bytes(&after)) > 0 {
                biting += 1;
            }
        }
        assert_eq!(
            biting, n,
            "only {biting} of {n} brushes changed the mesh, so the log is padded with no-ops and \
             the size arm would compare a fat description against geometry it did not produce"
        );

        // Timed the same way every C1 point is: `PASSES` round-robin visits of
        // `REPS` extractions, fastest kept. One point, so the round-robin is
        // degenerate, but the estimator has to be the same one or `M-281`'s rule
        // about comparing within one run is broken by the harness itself.
        let c2_rebuild = |d: &Desc<f32>| Gyroid {
            scale: d.params[0],
            iso: d.params[1],
        };
        let mut clock = Timing::unset();
        let mut timed = MeshBuffer::<f32>::new();
        extract_with(&mut mc32, &base, &brushes, &c2_grid, &mut timed);
        for _pass in 0..PASSES {
            clock = clock.best(time_point(
                &mut meter, &mut mc32, &base, &brushes, &c2_grid, &mut timed,
            ));
        }
        let point = characterise(
            &mut mc32,
            &base,
            &brushes,
            &desc,
            &c2_rebuild,
            (f64::from(n), clock),
        );
        assert_eq!(
            point.encoded.roundtrip_differing, 0,
            "C2 encoder round-trip"
        );
        assert_eq!(point.encoded.control_differing, 1, "C2 encoder control");
        assert_eq!(point.regen_differing, 0, "C2 regeneration");
        assert!(point.regen_control_differing > 0, "C2 regeneration control");

        let meshopt_bytes = point.encoded.vertex_bytes + point.encoded.index_bytes;
        let ratio = meshopt_bytes as f64 / point.desc_bytes as f64;
        let holds = ratio >= C2_SIZE_RATIO_BAR;
        c2_all_hold &= holds;
        c2_ratios.push((n, ratio));
        println!(
            "{:>10} {n:>8} {biting:>9} {:>8} {meshopt_bytes:>12} {:>10} {ratio:>10.4} {:>9}",
            bucket(n),
            point.triangles,
            point.desc_bytes,
            if holds { "HELD" } else { "FALSIFIED" }
        );

        c2_rows.push(
            Emit {
                arm: "c2_size",
                field: "gyroid_dug",
                precision: "f32",
                knob: f64::from(n),
                point: &point,
                ns_per_triangle_marginal: point.ms * 1.0e6 / point.triangles as f64,
                cycles_per_triangle_marginal: point
                    .window
                    .map(|w| w.cycles / (point.triangles as f64 * REPS as f64)),
                fit_a_ms: 0.0,
                fit_r2: 0.0,
                fit_points: 1,
                fit_sound: "n/a",
                brushes_biting: biting.to_string(),
                brushes_drawn: brushes_drawn.to_string(),
                replay_edits: 0,
                replay_min_triangles: String::from("n/a"),
                replay_total_triangles: String::from("n/a"),
                replay_ms: String::from("n/a"),
                bytes_differing: point.regen_differing,
                bytes_differing_scope: "regen-from-description",
                peer: "n/a",
                c1: "n/a",
                c2: if holds { "true" } else { "false" },
                c3: "n/a",
            }
            .row(),
        );
    }
    rows.extend(c2_rows);
    println!(
        "C2 {}: ratios {:?} against a bar of {C2_SIZE_RATIO_BAR}",
        if c2_all_hold { "HELD" } else { "FALSIFIED" },
        c2_ratios
            .iter()
            .map(|(n, r)| format!("{n}:{r:.2}"))
            .collect::<Vec<_>>()
    );

    // ── C3: the 10⁴-edit replay, and the peers on disk ─────────────────────
    println!("\nreplaying {REPLAY_EDITS} edits over {REPLAY_CHUNKS} chunks ...");
    let rep = replay();
    assert!(
        rep.min_triangles > 0,
        "VOID: some extraction in the replay emitted no triangle, so its digest is the hash of an \
         empty mesh and would agree across machines for no reason (M-44)"
    );
    assert_eq!(
        rep.digests.len(),
        REPLAY_EDITS * 8,
        "the digest stream is not one hash per edit"
    );
    println!(
        "  {:.1} s, {} triangles, min {} per extraction, {} digest bytes",
        rep.ms / 1.0e3,
        rep.total_triangles,
        rep.min_triangles,
        rep.digests.len()
    );

    let slug = machine_slug();
    let dir = repo_root().join("target");
    std::fs::create_dir_all(&dir).expect("target/ exists");
    let local_name = format!("experiment_p92_replay-{slug}.bin");
    std::fs::write(dir.join(&local_name), &rep.digests).expect("write the local stream");

    // The control peer: this machine's own stream with one byte flipped. It is
    // written every run so the comparison set is never empty and a live
    // comparator is proved rather than assumed.
    let control_name = String::from("experiment_p92_replay-control-onebyteflipped.bin");
    let mut flipped = rep.digests.clone();
    let at = flipped.len() / 2;
    flipped[at] ^= 0x01;
    std::fs::write(dir.join(&control_name), &flipped).expect("write the control stream");

    let mut peers: Vec<String> = std::fs::read_dir(&dir)
        .expect("read target/")
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .filter(|n| {
            n.starts_with("experiment_p92_replay-") && n.ends_with(".bin") && *n != local_name
        })
        .collect();
    peers.sort();

    let final_desc = rep.final_desc.clone();
    let final_lo = [
        final_desc.origin[0] as f32,
        final_desc.origin[1] as f32,
        final_desc.origin[2] as f32,
    ];
    let final_extent = (final_desc.cell_size * f64::from(final_desc.cells)) as f32;
    let final_encoded = encode_and_time(&rep.final_mesh, final_lo, final_extent);
    assert_eq!(
        final_encoded.roundtrip_differing, 0,
        "C3 encoder round-trip"
    );
    assert_eq!(final_encoded.control_differing, 1, "C3 encoder control");
    let (final_bytes, _) = encode_desc(&final_desc);
    // The geometry columns describe the **final chunk** -- the one whose
    // triangles were encoded and whose description was serialised -- so that
    // `triangles`, `meshopt_bytes`, `field_plus_log_bytes` and `size_ratio` are
    // four facts about one mesh. The replay's own aggregates go in
    // `replay_total_triangles` and `replay_ms`, and `extract_ms` is the mean
    // per-edit re-extraction, which is the number a game would budget.
    let final_point = Point {
        knob: REPLAY_EDITS as f64,
        triangles: rep.final_mesh.indices.len() / 3,
        vertices: rep.final_mesh.positions.len(),
        ms: rep.ms / REPLAY_EDITS as f64,
        window: None,
        mesh_bytes: mesh_bytes(&rep.final_mesh),
        encoded: final_encoded,
        desc_bytes: final_bytes.len(),
        brushes: final_desc.brushes.len() as u32,
        regen_differing: 0,
        regen_control_differing: 0,
    };

    println!(
        "\n{:>46} {:>12} {:>12} {:>9}",
        "peer", "diff bytes", "diff edits", "verdict"
    );
    let mut real_peers = 0usize;
    let mut c3_holds = true;
    for peer in &peers {
        let other = std::fs::read(dir.join(peer)).expect("read a peer stream");
        let differing = bytes_differing(&rep.digests, &other);
        let edits_differing = rep
            .digests
            .chunks(8)
            .zip(other.chunks(8))
            .filter(|(a, b)| a != b)
            .count()
            + rep.digests.len().abs_diff(other.len()).div_ceil(8);
        let is_control = peer.contains("control-onebyteflipped");
        if is_control {
            assert!(
                differing > 0,
                "the deliberately flipped control stream compares equal, so the comparator is \
                 blind and every zero beside it is worthless (M-44)"
            );
        } else {
            real_peers += 1;
            c3_holds &= differing == 0;
        }
        println!(
            "{peer:>46} {differing:>12} {edits_differing:>12} {:>9}",
            if is_control {
                "CONTROL"
            } else if differing == 0 {
                "IDENTICAL"
            } else {
                "DIVERGED"
            }
        );
        rows.push(
            Emit {
                arm: if is_control {
                    "c3_control"
                } else {
                    "c3_replay"
                },
                field: "gyroid_dug_replay",
                precision: "f64",
                knob: REPLAY_EDITS as f64,
                point: &final_point,
                ns_per_triangle_marginal: rep.ms * 1.0e6 / rep.total_triangles as f64,
                cycles_per_triangle_marginal: None,
                fit_a_ms: 0.0,
                fit_r2: 0.0,
                fit_points: 1,
                fit_sound: "n/a",
                brushes_biting: String::from("n/a"),
                brushes_drawn: String::from("n/a"),
                replay_edits: REPLAY_EDITS,
                replay_min_triangles: rep.min_triangles.to_string(),
                replay_total_triangles: rep.total_triangles.to_string(),
                replay_ms: format!("{:.3}", rep.ms),
                bytes_differing: differing,
                bytes_differing_scope: if is_control {
                    "control-onebyteflipped"
                } else {
                    "cross-machine-digest-stream"
                },
                peer,
                c1: "n/a",
                c2: "n/a",
                c3: if is_control {
                    "control"
                } else if differing == 0 {
                    "true"
                } else {
                    "false"
                },
            }
            .row(),
        );
    }
    println!(
        "C3: {real_peers} real peer stream(s) on disk beside {local_name}; {}",
        if real_peers == 0 {
            String::from(
                "BLOCKED -- no second machine's stream present, so nothing cross-machine \
                          was compared",
            )
        } else if c3_holds {
            String::from("HELD -- byte-identical")
        } else {
            String::from("FALSIFIED -- a byte differs")
        }
    );

    common::experiment::run(isomesh::experiment!("P-92"), |run| {
        for row in &rows {
            run.record(row);
        }
    });
}

/// An fbm description, so the two precisions do not spell it out twice.
fn fbm_desc<R: Wire>(field: &FbmTerrain<R>, origin: [R; 3], cell_size: R) -> Desc<R> {
    Desc {
        tag: 2,
        seed: field.seed,
        octaves: field.octaves,
        params: vec![
            field.lacunarity,
            field.gain,
            field.frequency,
            field.amplitude,
            field.base_height,
        ],
        cells: CHUNK_CELLS,
        origin,
        cell_size,
        brushes: Vec::new(),
    }
}

/// An fbm field, rebuilt from a decoded description.
fn fbm_from<R: Wire>(d: &Desc<R>) -> FbmTerrain<R> {
    FbmTerrain {
        seed: d.seed,
        octaves: d.octaves,
        lacunarity: d.params[0],
        gain: d.params[1],
        frequency: d.params[2],
        amplitude: d.params[3],
        base_height: d.params[4],
    }
}
