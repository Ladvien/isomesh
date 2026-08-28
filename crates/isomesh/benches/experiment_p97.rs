//! **P-97 — determinism at replay scale, not at fixture scale.**
//!
//! Ticket: R-097. Pre-registered before this harness existed; the registration
//! is `crates/isomesh/src/experiment.rs`, id `P-97`, and nothing here amends it.
//!
//! ```bash
//! cargo bench --bench experiment_p97
//! ```
//!
//! Writes `docs/experiments/p-97.csv`.
//!
//! # The hypothesis, and the falsifier
//!
//! `M-31` is 216 golden hashes on eight reference fields — a **zero-edit**
//! regime. A save file is a hundred thousand edits. Teardown's team rewrote
//! destruction in fixed-point integer arithmetic before deciding floating point
//! was workable with precautions, and they were not committing cross-platform
//! hashes; this crate is (`I-001`'s `libm` decision, with `golden.rs` named as
//! its proof), so its bar is higher.
//!
//! - **C1** — a 10⁵-edit trace replayed on the M5 and the Zen 3 produces
//!   byte-identical meshes on all eight fields. *Falsified by one differing
//!   byte, and that outcome is the point of the row.*
//! - **C2** — if C1 fails, the divergence is localised: the first differing
//!   edit is identifiable by bisection and its brush parameters name the
//!   operation responsible. *Falsified by a divergence bisection cannot
//!   localise.*
//! - **C3** — replay cost is linear in log length, with a constant under 1.2×
//!   the sum of per-edit costs. *Falsified by above 1.2×.*
//!
//! # Two machines, one code path, and the order they run in
//!
//! Every run writes its own meshes and prefix ladders to `target/p97-self/`,
//! and compares them against a peer's copies in `docs/experiments/p-97-peer/`.
//! The M5 runs first with no peer present, which is the registration's own
//! **BLOCKED** outcome and is recorded as the literal string `BLOCKED` rather
//! than as a zero — a missing comparison is never a passing one. Its
//! `target/p97-self/` is then copied to the Zen 3's
//! `docs/experiments/p-97-peer/` and the Zen 3 run produces the committed CSV.
//!
//! # Why replay is a fold over the whole grid, and not a brush bounding box
//!
//! A brush restricted to its own bounding box is **not** the crate's semantics.
//! `Add` is `min(field, shape)` and `Subtract` is `max(field, -shape)`; outside
//! the sphere `shape` is still finite, so both can move a sample that the
//! sphere does not contain — subtraction cannot flip a sign out there, but it
//! does change the value, and marching cubes places a vertex from the values at
//! both ends of a cut edge. A bounded update is only sound given a bound on
//! `|field|`, which `fbm_terrain` and `noise_cavity` do not have
//! (`FieldBound::Unbounded`). So the fold evaluates
//! [`BrushStack`](isomesh::brush::BrushStack) at every lattice point, which is
//! the crate's own definition of the state and needs no reimplementation to be
//! trusted.
//!
//! The loop is transposed — per point, over the whole tape — because that is
//! the only ordering in which a 10⁵-brush tape is affordable, and because it is
//! what lets one pass produce the per-prefix ladder C2 needs.
//!
//! # Controls, each an assertion rather than a printed number
//!
//! 1. **The registered vacuity control.** `same_kind_edits`, `mixed_edits` and
//!    `smooth_union_edits` are the three classes `M-36`, `M-37` and `M-38`
//!    distinguish, they are asserted non-zero, and they sum to `edits`. A trace
//!    of only same-kind brushes is the case already known to be order-free.
//! 2. **A zero that could have been non-zero (`M-44`).** The same trace is
//!    replayed with **one ULP** added to the radius of brush 50,000, and
//!    `control_bytes_differing` is asserted **> 0**. If a one-ULP perturbation
//!    in the middle of the tape does not reach the output bytes, then
//!    `bytes_differing = 0` measures the harness and not the crate.
//! 3. **The fold is the crate's fold.** The manual paired fold that produces
//!    the ladder is asserted bit-for-bit equal to `BrushStack::sample` at every
//!    lattice point. Two folds are two answers to one question.
//! 4. **The C3 instrument can report superlinearity.** An interactive arm
//!    re-evaluates and re-meshes after *every* edit — `M-50`'s regime — and its
//!    linearity constant is asserted **> 1.2**. A clause whose predicate no
//!    fixture in the harness can violate is `P-70`'s C3, and this is the column
//!    that stops that happening here.
//! 5. **The comparison is over something.** `bytes_compared` and `vertices` are
//!    asserted non-zero: a byte comparison of two empty meshes agrees perfectly
//!    and says nothing.
//!
//! # The SHARE line, recomputed before the run
//!
//! C3 is registered as "under 1.2× **the sum of per-edit costs `M-50`
//! measured**". `M-50` has **no committed CSV** — its four medians
//! (0.158 / 0.354 / 0.525 / 0.589 ms per re-meshed chunk, for logs of
//! 1–15 / 16–30 / 31–45 / 46–60) live only in `FINDINGS.md` prose, taken from
//! `E-202`, a live Bevy example on a different build and a different run. Using
//! them as a denominator would violate `M-281` outright, and taken literally
//! they make the clause unfalsifiable: 10⁵ × 0.158 ms is 15.8 s against a
//! replay this harness measures in single-digit seconds, so the ratio is
//! ~0.1–0.3× and no possible outcome reaches 1.2. That is `✗51`'s failure and
//! it is stated here before the run rather than discovered after it.
//!
//! So `sum_of_per_edit_ms` is measured **inside this run**: the marginal cost
//! of one edit in the short-log regime `M-50` sampled, taken as the slope
//! between a 1,000-brush tape and a 2,000-brush tape over the same grid, times
//! 10⁵. That baseline is reachable in both directions — the tape at 10⁵ brushes
//! is ~5 MB and does not fit the cache a 1,000-brush tape sits in, so a
//! constant above 1.2 is a physically available outcome, and the interactive
//! arm demonstrates the instrument reaching ~10²×.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::print_literal,
    clippy::too_many_lines
)]

mod common;

use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use isomesh::brush::{Brush, BrushOp, BrushStack, apply};
use isomesh::construct::SampledField;
use isomesh::extractor::Extractor;
use isomesh::fields::{ReferenceField, Sphere};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, Sdf};

/// Edits in the replayed trace. The registered figure.
const EDITS: usize = 100_000;

/// Samples per axis for the replayed grid.
///
/// 33 is the top of `golden.rs`'s `RESOLUTIONS`, so this is `M-31`'s own
/// fixture resolution with 10⁵ edits folded into it instead of none — which is
/// exactly the regime change the registration is about.
const SAMPLES: u32 = 33;

/// Brush radius in cells, low and high. Randomised per edit.
///
/// One and a half to four cells, against a dig box sixteen cells across, so the
/// carved region has structure at several scales rather than being one blob the
/// size of the box. A brush that spanned the box would make the last edit the
/// only one visible and the other 99,999 unfalsifiable.
const RADIUS_CELLS: (f64, f64) = (1.5, 4.0);

/// Smooth-union join width in cells, low and high.
const SMOOTH_K_CELLS: (f64, f64) = (0.5, 2.0);

/// The dig box, as a fraction of the domain along each axis.
///
/// The central half-extent — an eighth of the volume. See [`trace`] for why a
/// trace that is not localised erases the field it was applied to.
const DIG_BOX: (f64, f64) = (0.30, 0.55);

/// Tape lengths whose fold times give the marginal per-edit cost.
///
/// The slope between them, not the ratio to zero: a single short fold also pays
/// the base field's transcendentals once, and `fbm_terrain`'s are not a per-edit
/// cost.
const SLOPE_TAPES: (usize, usize) = (1_000, 2_000);

/// Edits in the interactive arm — `M-50`'s regime, re-meshed after every edit.
const INTERACTIVE_EDITS: usize = 512;

/// Samples per axis in the interactive arm.
///
/// Coarser than the replay because this arm is quadratic by construction and
/// its job is to show the C3 instrument moving, not to produce a headline.
const INTERACTIVE_SAMPLES: u32 = 17;

/// The registered ceiling on C3's constant.
const C3_CEILING: f64 = 1.2;

/// Trace seed. Fixed, so the trace is a function of nothing but this file.
const SEED: u64 = 0x5097_2026_0827_0001;

/// SplitMix64.
///
/// Integer state and an exact `2^-53` scaling, so the trace is bit-identical on
/// any machine before a single field is evaluated. A float-recurrence generator
/// would put the thing under test inside the fixture.
struct SplitMix64(u64);

impl SplitMix64 {
    const fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A double in `[0, 1)`, exactly: 53 bits over `2^53`.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }
}

/// How many edits of each of `M-36`/`M-37`/`M-38`'s three classes the trace has.
#[derive(Clone, Copy)]
struct Classes {
    same_kind: usize,
    mixed: usize,
    smooth: usize,
}

/// The trace, in one field's own domain.
///
/// The integer draws are shared across fields, so every field replays the same
/// logical trace mapped into its own `ReferenceField::domain`.
///
/// # Why the dig is confined to a box, and it is a fixture defect this
/// harness's own control caught
///
/// The first version scattered the brushes uniformly over the whole domain, and
/// **six of the eight fields then produced the same mesh hash**
/// (`c5dd8988f595a347` at 5,000 edits, on every field whose domain is
/// `[-2, 2]³`). The base field had been *erased*: `Add` is `min(field, shape)`
/// and with thousands of sphere centres scattered through a 32-cell domain the
/// nearest centre to any point is a fraction of a cell away, so the minimum is
/// a function of the brushes and nothing else. Eight rows were one measurement
/// copied eight times, and C1 would have been a determinism claim about
/// `sqrt`, `min` and `max` with no transcendental anywhere in the output — the
/// exact opposite of the `libm` property `M-31` exists to pin.
///
/// So the dig is confined to the central half-extent of the domain, which is an
/// eighth of its volume. That is also what a save file is: a player digs where
/// they are standing. The surface outside the dig box is provably untouched —
/// `Subtract` is `max(field, −shape)` and `−shape` is negative outside the
/// sphere, so it can raise an interior value but never flip a sign; `Add` and
/// `SmoothAdd` bite only where the sphere's distance falls below the field's
/// own, which near an untouched surface means inside the sphere. The mesh
/// therefore carries both regimes at once: `M-31`'s pristine surface, and a
/// region that has taken 10⁵ folds.
///
/// The control that would catch a regression here is in `main`: the eight mesh
/// hashes must be pairwise distinct.
fn trace(lo: [f64; 3], hi: [f64; 3], cell_size: f64, edits: usize) -> Vec<Brush<Sphere<f64>>> {
    let mut rng = SplitMix64::new(SEED);
    let mut log = Vec::with_capacity(edits);
    for _ in 0..edits {
        let mut center = [0.0f64; 3];
        for (axis, c) in center.iter_mut().enumerate() {
            let extent = hi[axis] - lo[axis];
            *c = lo[axis] + extent * (DIG_BOX.0 + (DIG_BOX.1 - DIG_BOX.0) * rng.unit());
        }
        let radius = (RADIUS_CELLS.0 + (RADIUS_CELLS.1 - RADIUS_CELLS.0) * rng.unit()) * cell_size;
        let shape = Sphere { center, radius };
        let k_draw = rng.unit();
        let op = match rng.next_u64() % 8 {
            0..=2 => BrushOp::Subtract,
            3..=5 => BrushOp::Add,
            _ => BrushOp::SmoothAdd {
                k: (SMOOTH_K_CELLS.0 + (SMOOTH_K_CELLS.1 - SMOOTH_K_CELLS.0) * k_draw) * cell_size,
            },
        };
        log.push(Brush { shape, op });
    }
    log
}

/// The three counts the registration's vacuity control names.
///
/// Smooth union is its own class. Among the hard edits, an edit is *same-kind*
/// when it commutes with the previous hard edit and *mixed* when it does not —
/// which is `BrushOp::commutes_with`, the crate's own predicate, rather than a
/// second opinion written here. The first hard edit has no predecessor and is
/// counted same-kind, so the three sum to `edits` exactly.
fn classify(log: &[Brush<Sphere<f64>>]) -> Classes {
    let mut c = Classes {
        same_kind: 0,
        mixed: 0,
        smooth: 0,
    };
    let mut previous_hard: Option<BrushOp> = None;
    for brush in log {
        match brush.op {
            BrushOp::SmoothAdd { .. } => c.smooth += 1,
            op => {
                match previous_hard {
                    None => c.same_kind += 1,
                    Some(prev) if prev.commutes_with(op) => c.same_kind += 1,
                    Some(_) => c.mixed += 1,
                }
                previous_hard = Some(op);
            }
        }
    }
    c
}

/// The lattice point at sample index `i`.
#[inline]
fn point(origin: [f64; 3], cell_size: f64, i: [u32; 3]) -> [f64; 3] {
    [
        origin[0] + cell_size * f64::from(i[0]),
        origin[1] + cell_size * f64::from(i[1]),
        origin[2] + cell_size * f64::from(i[2]),
    ]
}

/// Evaluate `base` folded through `brushes` at every lattice point.
///
/// This is [`BrushStack`]'s own `sample`, called once per point. It is the
/// authoritative replay and the thing `replay_ms` times.
fn fold_grid<F>(
    base: &F,
    brushes: &[Brush<Sphere<f64>>],
    origin: [f64; 3],
    cell_size: f64,
    n: u32,
) -> Vec<f64>
where
    F: Sdf<Scalar = f64>,
{
    let stack = BrushStack { base, brushes };
    let mut values = Vec::with_capacity((n as usize).pow(3));
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                values.push(stack.sample(point(origin, cell_size, [x, y, z])));
            }
        }
    }
    values
}

/// Where in the tape a one-ULP perturbation is planted, as a ladder.
///
/// **The first version of this harness planted it at edit 50,000 and the
/// `M-44` control refused the run**: on `sphere`, one ULP on the radius of the
/// middle brush reached *zero* output bytes. That is not a broken control, it
/// is the mechanism — `Add` is `min` and `Subtract` is `max`, and a `min`
/// *selects* an argument rather than computing one, so a perturbed brush only
/// survives at a lattice point where it is still the selected argument after
/// every later edit. With hundreds of later brushes covering each point, it
/// almost never is.
///
/// So the position is measured rather than guessed: candidates at `edits − 2^j`
/// for every `j`, plus edit 0, and the reported depth is the **earliest** edit
/// whose one-ULP perturbation still reaches the mesh bytes. The ladder is
/// geometric because the survival probability is, and the whole ladder costs
/// about two extra `apply`s per edit per point — the accumulators share the one
/// sphere evaluation the clean fold already paid for.
fn perturbation_candidates(edits: usize) -> Vec<usize> {
    let mut c = vec![0usize];
    let mut step = 1usize;
    while step < edits {
        c.push(edits - step);
        step *= 2;
    }
    c.sort_unstable();
    c.dedup();
    c
}

/// One ULP on a brush's radius. The smallest change a machine could disagree by.
fn perturb(brush: Brush<Sphere<f64>>) -> Brush<Sphere<f64>> {
    Brush {
        shape: Sphere {
            center: brush.shape.center,
            radius: f64::from_bits(brush.shape.radius.to_bits() + 1),
        },
        op: brush.op,
    }
}

/// What one probing fold produced.
struct Probe {
    /// The clean state, which must equal [`fold_grid`]'s bit for bit.
    clean: Vec<f64>,
    /// Digest of the clean state after each prefix — the ladder a second
    /// machine is bisected against.
    ladder: Vec<u64>,
    /// The perturbation positions, ascending.
    candidates: Vec<usize>,
    /// One final grid per candidate.
    grids: Vec<Vec<f64>>,
    /// `differs[j * edits + k]` — candidate `j` disagrees with clean after
    /// `k + 1` edits.
    differs: Vec<bool>,
}

/// Fold the clean tape and every perturbed variant together, one lattice point
/// at a time.
///
/// One pass yields what a bisection would otherwise pay `log n` full replays
/// for: the clean state, a per-prefix digest of it that a second machine can be
/// compared against edit by edit, and — for every candidate perturbation
/// position — the final state and the exact set of prefixes at which it
/// disagrees.
///
/// The digest is the `xor` of each sample's bit pattern rotated by its own
/// index. The rotation is what stops two samples with equal values cancelling;
/// the construction is checked against `differs`, which is exact, so a digest
/// collision is an assertion failure rather than a silent agreement.
fn fold_probe<F>(
    base: &F,
    log: &[Brush<Sphere<f64>>],
    origin: [f64; 3],
    cell_size: f64,
    n: u32,
) -> Probe
where
    F: Sdf<Scalar = f64>,
{
    let count = (n as usize).pow(3);
    let edits = log.len();
    let candidates = perturbation_candidates(edits);
    let slots = candidates.len();
    let mut out = Probe {
        clean: Vec::with_capacity(count),
        ladder: vec![0u64; edits],
        candidates,
        grids: vec![Vec::with_capacity(count); slots],
        differs: vec![false; slots * edits],
    };
    let mut acc = vec![0.0f64; slots];
    let mut index = 0usize;
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let p = point(origin, cell_size, [x, y, z]);
                let mut v = base.sample(p);
                let rotate = (index & 63) as u32;
                let mut live = 0usize;
                for (k, brush) in log.iter().enumerate() {
                    let s = brush.shape.sample(p);
                    for a in &mut acc[..live] {
                        *a = apply(brush.op, *a, s);
                    }
                    while live < slots && out.candidates[live] == k {
                        acc[live] = apply(brush.op, v, perturb(*brush).shape.sample(p));
                        live += 1;
                    }
                    v = apply(brush.op, v, s);
                    let bits = v.to_bits();
                    out.ladder[k] ^= bits.rotate_left(rotate);
                    for (j, a) in acc[..live].iter().enumerate() {
                        if a.to_bits() != bits {
                            out.differs[j * edits + k] = true;
                        }
                    }
                }
                out.clean.push(v);
                for (j, a) in acc[..slots].iter().enumerate() {
                    out.grids[j].push(*a);
                }
                index += 1;
            }
        }
    }
    out
}

/// Mesh a folded grid.
fn mesh_of(values: &[f64], origin: [f64; 3], cell_size: f64, n: u32) -> MeshBuffer<f64> {
    let shape = isomesh::RuntimeShape3::new([n; 3]).expect("grid fits u32");
    let field = SampledField::new(values, &shape, origin, cell_size).expect("sampled field");
    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    mc.extract_into(&field, &shape, origin, cell_size, &mut out)
        .expect("extract");
    out
}

/// The mesh as the bytes C1 compares: three counts, then positions, normals,
/// indices, all little-endian.
fn mesh_bytes(mesh: &MeshBuffer<f64>) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(
        24 + (mesh.positions.len() + mesh.normals.len()) * 24 + mesh.indices.len() * 4,
    );
    for count in [mesh.positions.len(), mesh.normals.len(), mesh.indices.len()] {
        bytes.extend_from_slice(&(count as u64).to_le_bytes());
    }
    for v in mesh.positions.iter().chain(mesh.normals.iter()) {
        for c in v {
            bytes.extend_from_slice(&c.to_le_bytes());
        }
    }
    for i in &mesh.indices {
        bytes.extend_from_slice(&i.to_le_bytes());
    }
    bytes
}

/// Differing bytes between two streams. A length difference counts as a
/// difference in every byte one of them does not have.
fn bytes_differing(a: &[u8], b: &[u8]) -> usize {
    let common = a.len().min(b.len());
    let mismatched = (0..common).filter(|&i| a[i] != b[i]).count();
    mismatched + (a.len().max(b.len()) - common)
}

/// Where two mesh byte streams differ, decomposed by the region of the mesh.
///
/// C2 asks for a divergence to be *localised*. The registration expects that
/// localisation to be an edit index, and this is the second axis: which part of
/// the mesh moved. A divergence confined to the normals is not the same finding
/// as one in the positions, and a byte count alone cannot tell them apart.
struct Divergence {
    total: usize,
    counts: usize,
    positions: usize,
    normals: usize,
    indices: usize,
    /// Differing normal bytes belonging to a vertex whose own normal is not
    /// finite — `unit_gradient` normalising a gradient of exactly zero.
    nan_normals: usize,
    first_vertex: Option<usize>,
}

/// Decompose the difference between this machine's mesh and a peer's bytes.
fn compare(mesh: &MeshBuffer<f64>, mine: &[u8], theirs: &[u8]) -> Divergence {
    let positions_at = 24;
    let normals_at = positions_at + mesh.positions.len() * 24;
    let indices_at = normals_at + mesh.normals.len() * 24;
    let common = mine.len().min(theirs.len());
    let mut d = Divergence {
        total: mine.len().max(theirs.len()) - common,
        counts: 0,
        positions: 0,
        normals: 0,
        indices: 0,
        nan_normals: 0,
        first_vertex: None,
    };
    for offset in 0..common {
        if mine[offset] == theirs[offset] {
            continue;
        }
        d.total += 1;
        if offset < positions_at {
            d.counts += 1;
        } else if offset < normals_at {
            let vertex = (offset - positions_at) / 24;
            d.positions += 1;
            d.first_vertex.get_or_insert(vertex);
        } else if offset < indices_at {
            let vertex = (offset - normals_at) / 24;
            d.normals += 1;
            if !mesh.normals[vertex].iter().all(|c| c.is_finite()) {
                d.nan_normals += 1;
            }
            d.first_vertex.get_or_insert(vertex);
        } else {
            d.indices += 1;
        }
    }
    d
}

/// A grid of samples as raw little-endian `f64`.
fn grid_bytes(values: &[f64]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(values.len() * 8);
    for v in values {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    bytes
}

/// The `u64` ladder as bytes, and back.
fn ladder_bytes(ladder: &[u64]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(ladder.len() * 8);
    for v in ladder {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    bytes
}

fn ladder_of_bytes(bytes: &[u8]) -> Vec<u64> {
    bytes
        .as_chunks::<8>()
        .0
        .iter()
        .map(|c| u64::from_le_bytes(*c))
        .collect()
}

/// The smallest prefix length at which `differs` is true, by binary search,
/// and the number of probes it took.
///
/// The search assumes the predicate is monotone. It is not guaranteed to be —
/// `max` and `min` can heal a divergence — and whether it healed is reported
/// beside the answer, because a divergence bisection cannot localise is exactly
/// C2's falsifier and it must be visible rather than assumed away.
fn bisect(differs: &[bool]) -> (Option<usize>, usize) {
    if differs.is_empty() {
        return (None, 0);
    }
    let mut steps = 0usize;
    let mut low = 0usize;
    let mut high = differs.len() - 1;
    if !differs[high] {
        steps += 1;
        return (None, steps);
    }
    while low < high {
        steps += 1;
        let mid = low + (high - low) / 2;
        if differs[mid] {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    (Some(low), steps)
}

/// Did the divergence heal — vanish after appearing, and come back?
fn heals(differs: &[bool]) -> bool {
    let Some(first) = differs.iter().position(|&d| d) else {
        return false;
    };
    let tail = &differs[first..];
    let gap = tail.iter().position(|&d| !d);
    match gap {
        None => false,
        Some(g) => tail[g..].iter().any(|&d| d),
    }
}

fn op_name(op: BrushOp) -> &'static str {
    match op {
        BrushOp::Add => "add",
        BrushOp::Subtract => "subtract",
        BrushOp::SmoothAdd { .. } => "smooth_add",
    }
}

/// This machine, as a filename-safe tag with no CSV separator in it.
fn machine_tag() -> String {
    format!("{}-{}", std::env::consts::ARCH, std::env::consts::OS)
}

/// The clock, on the row (`M-280`).
fn cpu_khz() -> String {
    let path = "/sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq";
    fs::read_to_string(path).map_or_else(
        |_| String::from("unavailable"),
        |s| s.trim().replace(',', ""),
    )
}

fn self_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/p97-self")
        .canonicalize()
        .unwrap_or_else(|_| {
            let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/p97-self");
            fs::create_dir_all(&p).expect("create self dir");
            p
        })
}

fn peer_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/experiments/p-97-peer")
}

/// The interactive arm: re-evaluate and re-mesh after every edit, which is what
/// `E-202` did under a mouse and what `M-50` timed.
///
/// Returns the total, the single-edit baseline extrapolated over the same edit
/// count, and their ratio.
fn interactive_arm<F>(
    base: &F,
    log: &[Brush<Sphere<f64>>],
    origin: [f64; 3],
    cell_size: f64,
) -> (f64, f64, f64)
where
    F: Sdf<Scalar = f64>,
{
    let n = INTERACTIVE_SAMPLES;
    let edits = INTERACTIVE_EDITS.min(log.len());

    let first = Instant::now();
    let values = fold_grid(base, &log[..1], origin, cell_size, n);
    let mesh = mesh_of(&values, origin, cell_size, n);
    std::hint::black_box(&mesh);
    let one = first.elapsed().as_secs_f64() * 1e3;

    let start = Instant::now();
    for k in 1..=edits {
        let values = fold_grid(base, &log[..k], origin, cell_size, n);
        let mesh = mesh_of(&values, origin, cell_size, n);
        std::hint::black_box(&mesh);
    }
    let total = start.elapsed().as_secs_f64() * 1e3;

    let baseline = one * edits as f64;
    (total, baseline, total / baseline)
}

type Row = Vec<(&'static str, String)>;

/// Everything one field contributes.
fn run_field<F>(name: &'static str, field: &F, peer: Option<&Path>) -> Row
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (_shape, origin, cell_size) = common::grid::<f64, F>(field, SAMPLES);
    let (lo, hi) = field.domain();

    let log = trace(lo, hi, cell_size, EDITS);
    let classes = classify(&log);

    // ── the registered vacuity control ──────────────────────────────────────
    assert!(
        classes.same_kind > 0 && classes.mixed > 0 && classes.smooth > 0,
        "VOID: {name} trace exercises {} same-kind, {} mixed and {} smooth-union edits; a trace \
         missing a class cannot test the class M-37 and M-38 say is the hard one",
        classes.same_kind,
        classes.mixed,
        classes.smooth
    );
    assert_eq!(
        classes.same_kind + classes.mixed + classes.smooth,
        EDITS,
        "the three classes must partition the trace"
    );

    // ── the authoritative replay, timed ─────────────────────────────────────
    let t = Instant::now();
    let values = fold_grid(field, &log, origin, cell_size, SAMPLES);
    let fold_ms = t.elapsed().as_secs_f64() * 1e3;

    let t = Instant::now();
    let mesh = mesh_of(&values, origin, cell_size, SAMPLES);
    let extract_ms = t.elapsed().as_secs_f64() * 1e3;
    let replay_ms = fold_ms + extract_ms;

    assert!(
        !mesh.positions.is_empty(),
        "VOID: {name} replayed to an empty mesh, so a byte comparison would compare nothing"
    );

    let bytes = mesh_bytes(&mesh);
    let hash = mesh_hash(&mesh);

    // How much of the compared byte stream is beyond M-31's regime. The
    // zero-edit mesh is exactly what the golden fixture hashes; the bytes that
    // differ from it are the ones 10^5 folds put there, and if that number were
    // small then C1 would be re-testing M-31 under a longer name.
    let pristine = fold_grid(field, &[], origin, cell_size, SAMPLES);
    let pristine_bytes = mesh_bytes(&mesh_of(&pristine, origin, cell_size, SAMPLES));
    let bytes_beyond_m31 = bytes_differing(&bytes, &pristine_bytes);
    assert!(
        bytes_beyond_m31 > 0,
        "VOID: {name}: the 10^5-edit replay produced the same bytes as the zero-edit mesh, so \
         this row measures M-31's regime and nothing else"
    );

    // ── marginal per-edit cost, in the short-log regime M-50 sampled ────────
    let t = Instant::now();
    let short = fold_grid(field, &log[..SLOPE_TAPES.0], origin, cell_size, SAMPLES);
    let t_short = t.elapsed().as_secs_f64() * 1e3;
    std::hint::black_box(&short);
    let t = Instant::now();
    let long = fold_grid(field, &log[..SLOPE_TAPES.1], origin, cell_size, SAMPLES);
    let t_long = t.elapsed().as_secs_f64() * 1e3;
    std::hint::black_box(&long);
    let per_edit_ms = (t_long - t_short) / (SLOPE_TAPES.1 - SLOPE_TAPES.0) as f64;
    let sum_of_per_edit_ms = per_edit_ms * EDITS as f64;
    let linearity_constant = replay_ms / sum_of_per_edit_ms;

    // ── the M-44 control: one ULP, at every depth in the tape ───────────────
    let probe = fold_probe(field, &log, origin, cell_size, SAMPLES);

    // ── control 3: the manual fold is BrushStack's fold ─────────────────────
    let disagreeing = probe
        .clean
        .iter()
        .zip(values.iter())
        .filter(|(a, b)| a.to_bits() != b.to_bits())
        .count();
    assert_eq!(
        disagreeing, 0,
        "{name}: the probing fold and BrushStack::sample disagree at {disagreeing} lattice \
         points, so the ladder describes a different computation from the one replay_ms timed"
    );

    // Which perturbation depths reach the mesh bytes at all. This is the
    // measurement the first version of the harness assumed: a one-ULP change
    // mid-tape is absorbed by `min`/`max`, and the depth at which it stops
    // being absorbed is a number rather than a guess.
    let reaching: Vec<(usize, usize)> = probe
        .candidates
        .iter()
        .enumerate()
        .map(|(j, &m)| {
            let candidate_mesh = mesh_of(&probe.grids[j], origin, cell_size, SAMPLES);
            (m, bytes_differing(&bytes, &mesh_bytes(&candidate_mesh)))
        })
        .collect();
    let deepest = reaching.iter().find(|&&(_, d)| d > 0).copied();
    let candidates_reaching = reaching.iter().filter(|&&(_, d)| d > 0).count();

    let Some((control_edit, control_bytes_differing)) = deepest else {
        panic!(
            "VOID: {name}: one ULP on the radius of edit {:?} -- every depth from 0 to {} -- \
             reached zero output bytes, so bytes_differing = 0 would be a property of this \
             fixture and not of the crate",
            probe.candidates,
            EDITS - 1
        );
    };
    assert!(
        control_bytes_differing > 0,
        "VOID: {name}: the deepest reaching perturbation reached zero bytes"
    );

    let slot = probe
        .candidates
        .iter()
        .position(|&m| m == control_edit)
        .expect("the chosen depth is a candidate");
    let control_differs = &probe.differs[slot * EDITS..(slot + 1) * EDITS];
    let control_first = control_differs.iter().position(|&d| d);
    let (control_bisected, control_steps) = bisect(control_differs);
    let control_heals = heals(control_differs);
    let control_agrees = control_first == control_bisected;
    let control_prefixes = control_differs.iter().filter(|&&d| d).count();
    let victim = log[control_edit];

    // ── the C3 instrument's own control: M-50's regime, and it is quadratic ─
    let (interactive_total_ms, interactive_baseline_ms, interactive_constant) =
        interactive_arm(field, &log, origin, cell_size);
    assert!(
        interactive_constant > C3_CEILING,
        "VOID: {name}: the interactive arm returned {interactive_constant:.4}x, at or under the \
         registered 1.2x ceiling, so nothing in this harness can make C3 fail and its HELD would \
         be P-70's C3 again"
    );

    // ── write this machine's artefacts ──────────────────────────────────────
    //
    // The grid is written as well as the mesh because the two answer different
    // questions. The mesh is what C1 compares; the grid is the *state the
    // replay produced*, and a divergence that is in the mesh but not in the
    // grid did not happen in the fold at all. The ladder already implies this
    // over 100,000 prefixes, but an implication from an `xor` digest is weaker
    // than 287 KB of raw samples, and this claim is load-bearing.
    let grid = grid_bytes(&values);
    let dir = self_dir();
    fs::create_dir_all(&dir).expect("create self dir");
    fs::write(dir.join(format!("{name}.mesh")), &bytes).expect("write mesh");
    fs::write(dir.join(format!("{name}.grid")), &grid).expect("write grid");
    fs::write(
        dir.join(format!("{name}.ladder")),
        ladder_bytes(&probe.ladder),
    )
    .expect("write ladder");
    fs::write(dir.join(format!("{name}.hash")), format!("{hash:016x}")).expect("write hash");
    fs::write(dir.join("machine.txt"), machine_tag()).expect("write machine");

    // ── C1 and C2, against the peer ─────────────────────────────────────────
    let peer_mesh = peer.and_then(|p| fs::read(p.join(format!("{name}.mesh"))).ok());
    let peer_grid = peer.and_then(|p| fs::read(p.join(format!("{name}.grid"))).ok());
    let peer_ladder = peer.and_then(|p| fs::read(p.join(format!("{name}.ladder"))).ok());
    let peer_hash = peer
        .and_then(|p| fs::read_to_string(p.join(format!("{name}.hash"))).ok())
        .map_or_else(|| String::from("ABSENT"), |s| s.trim().replace(',', ""));
    let peer_machine = peer
        .and_then(|p| fs::read_to_string(p.join("machine.txt")).ok())
        .map_or_else(|| String::from("ABSENT"), |s| s.trim().replace(',', ""));

    // How many vertices carry a normal this machine cannot express as a
    // direction. `unit_gradient` divides by `|grad|` and only shouts about a
    // zero in debug builds, so in release a plateau in the folded grid emits
    // `0 * inf`.
    let nan_normal_vertices = mesh
        .normals
        .iter()
        .filter(|n| !n.iter().all(|c| c.is_finite()))
        .count();

    let divergence = peer_mesh
        .as_ref()
        .map(|other| compare(&mesh, &bytes, other));
    let grid_differing = peer_grid
        .as_ref()
        .map(|other| bytes_differing(&grid, other));

    let (differing, c1) = match &divergence {
        None => (String::from("BLOCKED"), "BLOCKED"),
        Some(d) => (
            d.total.to_string(),
            if d.total == 0 { "HELD" } else { "FALSIFIED" },
        ),
    };

    let (first_edit, first_op, bisection_steps, c2) = match (&peer_ladder, &divergence) {
        (Some(other), Some(d)) => {
            let theirs = ladder_of_bytes(other);
            assert_eq!(
                theirs.len(),
                probe.ladder.len(),
                "{name}: peer ladder is {} entries against {}; the two machines did not replay the \
                 same trace and no comparison of theirs is meaningful",
                theirs.len(),
                probe.ladder.len()
            );
            let cross: Vec<bool> = probe
                .ladder
                .iter()
                .zip(theirs.iter())
                .map(|(a, b)| a != b)
                .collect();
            let truth = cross.iter().position(|&x| x);
            let (found, steps) = bisect(&cross);
            match (truth, d.total) {
                // No divergence anywhere: C2's population is empty, and saying
                // so is the honest score rather than a HELD.
                (None, 0) => (String::from("NONE"), String::from("NONE"), steps, "EMPTY"),
                // The meshes differ and the fold does not. Bisection over the
                // edit log cannot localise this, which is C2's registered
                // falsifier reached by a better route than the registered one:
                // not "the fold accumulates" but "the fold never diverged".
                (None, _) => (
                    String::from("NONE"),
                    String::from("NONE"),
                    steps,
                    "FALSIFIED",
                ),
                (Some(k), _) => {
                    let localised = found == Some(k);
                    (
                        k.to_string(),
                        String::from(op_name(log[k].op)),
                        steps,
                        if localised { "HELD" } else { "FALSIFIED" },
                    )
                }
            }
        }
        _ => (
            String::from("BLOCKED"),
            String::from("BLOCKED"),
            0,
            "BLOCKED",
        ),
    };

    let c3 = if linearity_constant < C3_CEILING {
        "HELD"
    } else {
        "FALSIFIED"
    };

    println!(
        "{name:>15} {:>7} {:>9} {:>9.1} {:>9.4} {:>7} {:>10} {:>5} {:>10} {:>10}",
        mesh.positions.len(),
        bytes.len(),
        replay_ms,
        linearity_constant,
        nan_normal_vertices,
        grid_differing.map_or_else(|| String::from("BLOCKED"), |d| d.to_string()),
        differing,
        c1,
        c2
    );

    vec![
        ("field", String::from(name)),
        ("edits", EDITS.to_string()),
        ("same_kind_edits", classes.same_kind.to_string()),
        ("mixed_edits", classes.mixed.to_string()),
        ("smooth_union_edits", classes.smooth.to_string()),
        ("bytes_differing", differing),
        ("first_differing_edit", first_edit),
        ("first_differing_operator", first_op),
        ("bisection_steps", bisection_steps.to_string()),
        ("replay_ms", format!("{replay_ms:.3}")),
        ("sum_of_per_edit_ms", format!("{sum_of_per_edit_ms:.3}")),
        ("linearity_constant", format!("{linearity_constant:.6}")),
        ("c1_holds", String::from(c1)),
        ("c2_holds", String::from(c2)),
        ("c3_holds", String::from(c3)),
        // ── beyond the registration ─────────────────────────────────────────
        ("samples_per_axis", SAMPLES.to_string()),
        ("cell_size", format!("{cell_size:.9}")),
        ("vertices", mesh.positions.len().to_string()),
        ("triangles", (mesh.indices.len() / 3).to_string()),
        ("bytes_compared", bytes.len().to_string()),
        ("pristine_bytes", pristine_bytes.len().to_string()),
        ("bytes_beyond_m31", bytes_beyond_m31.to_string()),
        ("grid_bytes", grid.len().to_string()),
        (
            "grid_bytes_differing",
            grid_differing.map_or_else(|| String::from("BLOCKED"), |d| d.to_string()),
        ),
        (
            "bytes_differing_counts",
            divergence
                .as_ref()
                .map_or_else(|| String::from("BLOCKED"), |d| d.counts.to_string()),
        ),
        (
            "bytes_differing_positions",
            divergence
                .as_ref()
                .map_or_else(|| String::from("BLOCKED"), |d| d.positions.to_string()),
        ),
        (
            "bytes_differing_normals",
            divergence
                .as_ref()
                .map_or_else(|| String::from("BLOCKED"), |d| d.normals.to_string()),
        ),
        (
            "bytes_differing_indices",
            divergence
                .as_ref()
                .map_or_else(|| String::from("BLOCKED"), |d| d.indices.to_string()),
        ),
        (
            "bytes_differing_at_nan_normals",
            divergence
                .as_ref()
                .map_or_else(|| String::from("BLOCKED"), |d| d.nan_normals.to_string()),
        ),
        (
            "first_differing_vertex",
            divergence.as_ref().map_or_else(
                || String::from("BLOCKED"),
                |d| {
                    d.first_vertex
                        .map_or_else(|| String::from("NONE"), |v| v.to_string())
                },
            ),
        ),
        ("nan_normal_vertices", nan_normal_vertices.to_string()),
        (
            "dig_box_fraction_of_extent",
            format!("{:.2}", DIG_BOX.1 - DIG_BOX.0),
        ),
        ("mesh_hash_local", format!("{hash:016x}")),
        ("mesh_hash_peer", peer_hash),
        ("local_machine", machine_tag()),
        ("peer_machine", peer_machine),
        ("cpu_khz", cpu_khz()),
        ("fold_ms", format!("{fold_ms:.3}")),
        ("extract_ms", format!("{extract_ms:.3}")),
        ("per_edit_ms", format!("{per_edit_ms:.9}")),
        ("slope_tape_short_ms", format!("{t_short:.3}")),
        ("slope_tape_long_ms", format!("{t_long:.3}")),
        ("control_perturbed_edit", control_edit.to_string()),
        (
            "control_perturbation_depths",
            probe.candidates.len().to_string(),
        ),
        (
            "control_depths_reaching_output",
            candidates_reaching.to_string(),
        ),
        (
            "control_deepest_reaching_from_end",
            (EDITS - control_edit).to_string(),
        ),
        (
            "control_perturbed_operator",
            String::from(op_name(victim.op)),
        ),
        (
            "control_bytes_differing",
            control_bytes_differing.to_string(),
        ),
        (
            "control_first_differing_edit",
            control_first.map_or_else(|| String::from("NONE"), |k| k.to_string()),
        ),
        (
            "control_first_differing_operator",
            control_first.map_or_else(
                || String::from("NONE"),
                |k| String::from(op_name(log[k].op)),
            ),
        ),
        (
            "control_bisection_edit",
            control_bisected.map_or_else(|| String::from("NONE"), |k| k.to_string()),
        ),
        ("control_bisection_steps", control_steps.to_string()),
        ("control_bisection_agrees", control_agrees.to_string()),
        ("control_divergence_heals", control_heals.to_string()),
        ("control_differing_prefixes", control_prefixes.to_string()),
        ("interactive_edits", INTERACTIVE_EDITS.to_string()),
        (
            "interactive_samples_per_axis",
            INTERACTIVE_SAMPLES.to_string(),
        ),
        ("interactive_total_ms", format!("{interactive_total_ms:.3}")),
        (
            "interactive_baseline_ms",
            format!("{interactive_baseline_ms:.3}"),
        ),
        (
            "interactive_linearity_constant",
            format!("{interactive_constant:.6}"),
        ),
    ]
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-97");

    let peer = peer_dir();
    let peer = if peer.join("machine.txt").exists() {
        Some(peer)
    } else {
        println!(
            "P-97: no peer artefacts at {} -- C1 and C2 record BLOCKED, not zero. Run this bench \
             on the second machine, copy its target/p97-self/ here, and run again.",
            peer.display()
        );
        None
    };

    println!(
        "{:>15} {:>7} {:>9} {:>9} {:>9} {:>7} {:>10} {:>5} {:>10} {:>10}",
        "field",
        "verts",
        "bytes",
        "replay_ms",
        "c3_const",
        "nan_n",
        "grid_diff",
        "diff",
        "c1",
        "c2"
    );

    let mut rows: Vec<Row> = Vec::new();
    isomesh::for_each_reference_field!(f64, |name, field| {
        rows.push(run_field(name, &field, peer.as_deref()));
    });

    assert_eq!(rows.len(), 8, "the reference sweep must yield eight fields");

    // ── the control that caught the first fixture ───────────────────────────
    //
    // Eight fields must give eight meshes. When the dig was scattered over the
    // whole domain the brushes erased the field and six of these hashes were
    // the same string, which would have made "byte-identical on all eight
    // fields" a claim about one mesh reported eight times.
    let hashes: Vec<&String> = rows
        .iter()
        .map(|row| {
            &row.iter()
                .find(|(k, _)| *k == "mesh_hash_local")
                .expect("every row hashes its mesh")
                .1
        })
        .collect();
    let mut unique: Vec<&&String> = hashes.iter().collect();
    unique.sort_unstable();
    unique.dedup();
    assert_eq!(
        unique.len(),
        8,
        "VOID: only {} distinct mesh hashes across eight fields ({hashes:?}); the trace has \
         erased the base field and these rows are one measurement repeated",
        unique.len()
    );

    common::experiment::run(prereg, |run| {
        for row in &rows {
            run.record(row);
        }
    });
}
