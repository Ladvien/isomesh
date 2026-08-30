//! **P-117 — FMA contraction as a latent golden-hash divergence — a risk audit, not an optimisation.**
//!
//! Ticket: R-117. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p117
//! ```
//!
//! Writes `docs/experiments/p-117.csv`.
//!
//! # What was missing
//!
//! Three separate things, and the crate had an **assertion** standing where each
//! of them should have had a measurement.
//!
//! 1. **Nobody had enumerated the multiply-add shapes on the extraction path.**
//!    `real.rs:45` says `mul_add` is deliberately absent from `isomesh::Real`,
//!    `fields/noise.rs:17` says *"Rust does not contract `a * b + c` into a fused
//!    multiply-add on its own, and nothing here asks for one"*, and
//!    `predicates.rs:28` says the same and then adds *"but this is asserted
//!    rather than assumed"* — pointing at `two_product_is_exact`, which is the
//!    only check in the repository that would notice contraction, and it covers
//!    one function in a module the extractors never call. An absent `mul_add` is
//!    a statement about what the source **requests**. It is not a statement
//!    about what the arithmetic is **sensitive to**. This harness names every
//!    multiply-add shape on the path, evaluates each one both ways on inputs
//!    harvested from real extraction, and puts the ULP gap on a row.
//!
//! 2. **Nobody had looked at the machine code.** `scripts/p69_asm.sh` exists
//!    because P-69 registered exactly this discipline for vectorisation: *"A
//!    Criterion delta cannot distinguish a vectorised loop from a lucky one."*
//!    The same holds here — a golden hash that has not moved cannot distinguish
//!    a compiler that refuses to contract from a compiler that contracts
//!    harmlessly. So this harness disassembles **its own running binary** with
//!    `objdump -d`, counts fused-multiply-add mnemonics per symbol, and reports
//!    the count for a probe that must fuse beside the count for the shapes that
//!    must not.
//!
//! 3. **`M-31`'s cross-platform claim had never been measured in one sitting.**
//!    `golden.rs:49-54` states the portability the fixture asserts and names
//!    `CLAUDE.md` as where the decision is recorded. CI runs Linux and macOS
//!    against one committed fixture, which is evidence — but it is evidence
//!    spread across two jobs and two moments. This harness recomputes all 216
//!    hashes on both machines from the same source revision, and the peer's 216
//!    are committed as an **input** (`docs/experiments/p-117-m5-fma.txt`), the
//!    same shape as `docs/experiments/p-83-m5-hashes.txt`.
//!
//! `M-31`'s own text says 63 hashes. That was `T-007`'s count. The fixture is
//! **216** today — 8 reference fields × 9 algorithms × 3 resolutions {17, 25, 33}
//! (`golden.rs:73`, `:122`, `fields/mod.rs:212`) — and the two must not be
//! conflated. The gate is proven able to fire: `P-61` moved 135 of these 216.
//!
//! Every mechanism this file needs out of a private module is **copied here**,
//! with the source line it came from on the row that uses it.
//! `crates/isomesh/src/` is read-only for the whole of Phase 25, and a copy
//! whose line number is a column is auditable in a way a `pub` would not be.
//!
//! # SHARE
//!
//! **None, and that is registered rather than discovered.** This row moves no
//! time and claims no ratio, so no clause is a fraction of a total and there is
//! no Amdahl ceiling to compute. What stands in a share's place is one column
//! per clause, named here so a reader can check each verdict against the number
//! that produced it:
//!
//! * **C1** — at least one expression in `crates/isomesh/src/` is
//!   contraction-sensitive. **Column: `ulp_difference`**, per expression site,
//!   maximised over every reachable input. C1 holds when any site whose `file`
//!   is under `crates/isomesh/src/` reads above zero. `inputs_tried` and
//!   `inputs_differing` are the denominator and the numerator behind it, and
//!   `worst_input` carries the exact bits so the number can be reproduced by
//!   hand.
//! * **C2** — the 216 golden hashes are identical on the two targets. **Column:
//!   `hashes_identical`, over `golden_rows`.** `golden_hashes_differing` is the
//!   count behind the boolean and `golden_matches_committed_fixture` says
//!   separately whether each arm reproduces `crates/isomesh/golden_hashes.json`.
//! * **C3** — a divergence is attributed to a named expression and reproduced by
//!   an isolated probe. **Column: `probe_reproduced`**, beside
//!   `expression_site`, `file` and `line`. `contracted_in_codegen` is the
//!   machine-code half of the same question, read from that site's own
//!   `#[inline(never)]` probe.
//!
//! **The vacuity controls, and both of them are two-sided.**
//!
//! * *The disassembly reader.* `fmadd_in_known_fused_probe` must be true — a
//!   probe written with an explicit `mul_add` must be seen fusing. A reader that
//!   cannot see an FMA where one must exist cannot claim to see its absence
//!   anywhere else. `known_fused_probe_evidence` names the mnemonic actually
//!   observed and `known_fused_dispatch_symbol` names the symbol it was observed
//!   in, rather than asserting the category.
//! * *The ULP comparator.* Two control rows, one in each direction.
//!   `control::exact_products` is a **real** site — `perlin`'s `g · d`, whose
//!   gradient components are exactly `0` and `±1` so every product is exact —
//!   and it must read `ulp_difference = 0`, because fusing an exact product
//!   cannot change a result. A non-zero there means the comparator cries wolf.
//!   `control::searched_separator` is a triple **searched** for with SplitMix64
//!   rather than chosen by hand — `cube.rs:349-356`'s precedent, and its stated
//!   reason: the first version of that test picked its separating pair by hand
//!   and landed on values where both forms agree exactly — and it must read
//!   above zero, or the comparator is blind. Both are `assert!`ed before a row
//!   is written.

mod common;

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

use isomesh::dual_contouring::DualContouring;
use isomesh::fields::{FbmTerrain, ReferenceField};
use isomesh::greedy_quads::GreedyQuads;
use isomesh::hermite::HermiteCell;
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::marching_cubes::{FaceAmbiguity, InteriorAmbiguity, MarchingCubes};
use isomesh::marching_tetrahedra::MarchingTetrahedra;
use isomesh::subgrid::extract::SubgridMarchingTetrahedra;
use isomesh::surface_nets::SurfaceNets;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

// ── the two targets the registration names, and nothing else ────────────────

/// The Linux arm.
const X86: &str = "x86_64-unknown-linux-gnu";
/// The aarch64 arm, measured on `mac_air`.
const ARM: &str = "aarch64-apple-darwin";

/// Which of the two this binary is.
///
/// A third machine is not this row's fixture. `M-31`'s claim is about these two
/// targets, and a number from a third filed under a `target` column the
/// registration never names is how a fixture quietly becomes a different one.
const HOST_TARGET: &str = if cfg!(all(target_arch = "x86_64", target_os = "linux")) {
    X86
} else if cfg!(all(target_arch = "aarch64", target_os = "macos")) {
    ARM
} else {
    "unregistered"
};

/// Where the peer arm's measurement is committed, as an **input**.
const PEER_ARTEFACT: &str = "docs/experiments/p-117-m5-fma.txt";

// ── mechanisms copied out of private modules, with their source lines ───────

/// `cube.rs:171`. **A sample of exactly zero is outside.**
#[inline]
fn is_inside(value: f64) -> bool {
    value < 0.0
}

/// `cube.rs:26`, expanded: the twelve edges, lower corner first, grouped by axis.
const EDGE_CORNERS: [[u8; 2]; 12] = {
    let mut out = [[0u8; 2]; 12];
    let mut n = 0;
    let mut axis = 0u8;
    while axis < 3 {
        let mut c = 0u8;
        while c < 8 {
            if c & (1 << axis) == 0 {
                out[n] = [c, c | (1 << axis)];
                n += 1;
            }
            c += 1;
        }
        axis += 1;
    }
    out
};

/// `cube.rs:149`.
#[inline]
const fn corner_offset(corner: u8) -> [u32; 3] {
    [
        (corner & 1) as u32,
        ((corner >> 1) & 1) as u32,
        ((corner >> 2) & 1) as u32,
    ]
}

/// `cube.rs:222`. The centred offset, `((a + b)/2)/(a − b)`.
#[inline]
fn edge_offset(a: f64, b: f64) -> f64 {
    ((a + b) * 0.5) / (a - b)
}

/// `cube.rs:234`. **A site.**
#[inline]
fn place(lo: f64, hi: f64, d: f64) -> f64 {
    (lo + hi) * 0.5 + (hi - lo) * d
}

/// `equivariant.rs:80`.
#[inline]
fn precedes(a: f64, b: f64) -> bool {
    match a.abs().total_cmp(&b.abs()) {
        core::cmp::Ordering::Less => true,
        core::cmp::Ordering::Greater => false,
        core::cmp::Ordering::Equal => a.total_cmp(&b) == core::cmp::Ordering::Less,
    }
}

/// The permutation `equivariant.rs`'s sorting network produces, as an insertion
/// sort so one function covers the three-, five- and twelve-element widths this
/// file needs. `precedes` is a total order, so the two agree.
fn sort_by_magnitude(t: &mut [f64]) {
    for i in 1..t.len() {
        let mut j = i;
        while j > 0 && precedes(t[j], t[j - 1]) {
            t.swap(j, j - 1);
            j -= 1;
        }
    }
}

/// `equivariant.rs:92`. Sorts in place: the caller owns the buffer, so nothing
/// on the harvest's hot path allocates.
fn sum_equivariant(t: &mut [f64]) -> f64 {
    sort_by_magnitude(t);
    let mut acc = 0.0;
    for value in t {
        acc += *value;
    }
    acc
}

/// `dual_contouring/solve.rs:114`.
fn mul_equivariant(mut t: [f64; 3]) -> f64 {
    sort_by_magnitude(&mut t);
    (t[0] * t[1]) * t[2]
}

/// `real.rs`'s `SPLITTER` for `f64`: `2^27 + 1`, Dekker–Veltkamp.
const SPLITTER: f64 = 134_217_729.0;

/// `predicates.rs:85`.
#[inline]
fn split(a: f64) -> (f64, f64) {
    let c = SPLITTER * a;
    let a_big = c - a;
    let a_hi = c - a_big;
    (a_hi, a - a_hi)
}

/// `dual_contouring/solve.rs:85`.
const LAMBDA: f64 = 0.01;

/// `cube.rs`'s `EDGE_COUNT`.
const EDGE_COUNT: usize = 12;

// ── the ULP comparator ──────────────────────────────────────────────────────

/// A monotone `u64` key over the `f64` total order, so a gap is a count of
/// representable numbers rather than a subtraction of two floats.
#[inline]
fn ordering_key(x: f64) -> u64 {
    let b = x.to_bits();
    if b & 0x8000_0000_0000_0000 == 0 {
        b | 0x8000_0000_0000_0000
    } else {
        !b
    }
}

/// How many representable `f64`s lie between `a` and `b`.
#[inline]
fn ulp_gap(a: f64, b: f64) -> u64 {
    ordering_key(a).abs_diff(ordering_key(b))
}

// ── the isolated probes, and the disassembly arm's subjects ─────────────────
//
// `#[inline(never)]` for `scripts/p69_asm.sh`'s stated reason: with thin LTO
// these otherwise inline into their callers and leave no symbol to inspect, so a
// dump without them shows nothing and looks like evidence. They are **not**
// `#[no_mangle]`, because `unsafe_code = "forbid"` is a workspace lint and
// `#[unsafe(no_mangle)]` trips it; the mangled symbol carries the function name
// as a substring, which is all `symbol_body` needs.
//
// They serve twice. First as C3's isolated probe: the worst input recorded on a
// row is pushed back through the probe and the result must agree bit for bit
// with the harness's own evaluation, which is what makes a divergence
// *attributable* rather than merely observed. Second as what gets disassembled.
//
// There is no separate probe for `origin + cell_size * n`: it is the same
// expression as `acc + m * v`, LLVM's function merger folds the two into one
// symbol, and two names for one body would mean one of them read
// `symbol-not-found` and looked like a failed measurement.

/// A shape that **must** fuse. The vacuity control for the whole reader.
#[inline(never)]
fn p117_probe_known_fused(a: f64, b: f64, c: f64) -> f64 {
    a.mul_add(b, c)
}

/// The question, in its simplest form: does this build contract this?
#[inline(never)]
fn p117_probe_written(a: f64, b: f64, c: f64) -> f64 {
    a * b + c
}

/// `cube.rs:234`.
#[inline(never)]
fn p117_probe_place(lo: f64, hi: f64, d: f64) -> f64 {
    (lo + hi) * 0.5 + (hi - lo) * d
}

/// `vec3.rs:29`.
#[inline(never)]
fn p117_probe_dot(a0: f64, b0: f64, a1: f64, b1: f64, a2: f64, b2: f64) -> f64 {
    a0 * b0 + a1 * b1 + a2 * b2
}

/// `dual_contouring/solve.rs:104` — the same dot product through the
/// magnitude-sorted reduction A-016 put there.
#[inline(never)]
fn p117_probe_dot_equivariant(a0: f64, b0: f64, a1: f64, b1: f64, a2: f64, b2: f64) -> f64 {
    sum_equivariant(&mut [a0 * b0, a1 * b1, a2 * b2])
}

/// `vec3.rs:33` and `dual_contouring/solve.rs:186` — the 2×2 determinant.
#[inline(never)]
fn p117_probe_cross(a: f64, b: f64, c: f64, d: f64) -> f64 {
    a * b - c * d
}

/// `fields/noise.rs:96` — Perlin's quintic fade, in Horner form.
#[inline(never)]
fn p117_probe_fade(t: f64) -> f64 {
    t * (t * 6.0 - 15.0) + 10.0
}

/// `fields/noise.rs:168`, `:243`, `fields/mod.rs:1360`,
/// `marching_cubes/mod.rs:761`, `dual.rs:616`, `hermite.rs:99` — accumulate one
/// product onto one addend.
#[inline(never)]
fn p117_probe_accumulate(acc: f64, m: f64, v: f64) -> f64 {
    acc + m * v
}

/// `fields/noise.rs:237` — two products either side of one add.
#[inline(never)]
fn p117_probe_two_products(p: f64, freq: f64, k: f64, off: f64) -> f64 {
    p * freq + k * off
}

/// `predicates.rs:97-105` split in two, so the site's `written`, its two fused
/// variants and the `#[inline(never)]` probe all read the same five values and
/// differ only in how the four subtractions round.
fn two_product_parts(a: f64, b: f64) -> (f64, f64, f64, f64, f64) {
    let x = a * b;
    let (a_hi, a_lo) = split(a);
    let (b_hi, b_lo) = split(b);
    (x, a_hi, a_lo, b_hi, b_lo)
}

/// `predicates.rs:104` — Shewchuk's roundoff, as written.
fn two_product_roundoff(a: f64, b: f64) -> f64 {
    let (x, a_hi, a_lo, b_hi, b_lo) = two_product_parts(a, b);
    let err1 = x - (a_hi * b_hi);
    let err2 = err1 - (a_lo * b_hi);
    let err3 = err2 - (a_hi * b_lo);
    (a_lo * b_lo) - err3
}

/// `predicates.rs:104`, isolated.
#[inline(never)]
fn p117_probe_two_product_err(a: f64, b: f64) -> f64 {
    two_product_roundoff(a, b)
}

/// Every probe the disassembly arm reads, so a merged or missing symbol is
/// visible as data rather than inferred from a blank.
const PROBES: [&str; 10] = [
    "p117_probe_known_fused",
    "p117_probe_written",
    "p117_probe_place",
    "p117_probe_dot",
    "p117_probe_dot_equivariant",
    "p117_probe_cross",
    "p117_probe_fade",
    "p117_probe_accumulate",
    "p117_probe_two_products",
    "p117_probe_two_product_err",
];

// ── one expression site ─────────────────────────────────────────────────────

/// One evaluation of a six-slot input word.
type Evaluate = fn(&[f64; 6]) -> f64;

/// The named fusions of one written expression.
type Variants = &'static [(&'static str, Evaluate)];

/// A multiply-add shape in `crates/isomesh/src/`, with every fusion a compiler
/// could plausibly choose written out by hand.
struct Site {
    name: &'static str,
    file: &'static str,
    line: u32,
    /// The expression, spelled the way the source spells it. No commas: the CSV
    /// writer does not quote.
    shape: &'static str,
    /// Can the crate reach this expression at all, from any entry point?
    reachable: bool,
    /// Is it reached by `extract` on the eight reference fields?
    on_path: bool,
    /// Does the product also feed a comparison? A-016's magnitude sort reads
    /// every product before summing it, and an FMA wants the multiply to feed
    /// only the add — so a `true` here is a structural shield the source
    /// acquired for an entirely different reason.
    product_also_compared: bool,
    /// Which `#[inline(never)]` probe carries this shape verbatim.
    probe: &'static str,
    /// The expression as written.
    written: Evaluate,
    /// Every fusion a contracting compiler could pick, named.
    variants: Variants,
    /// The probe, invoked on a recorded input, for C3's reproduction.
    reproduce: Evaluate,
}

/// What one site read, maximised over every reachable input.
#[derive(Clone, Default)]
struct Reading {
    tried: u64,
    differing: u64,
    worst_input: [f64; 6],
    worst_written: f64,
    worst_fused: f64,
    worst_ulp: u64,
    worst_variant: &'static str,
    seen: bool,
}

impl Reading {
    /// Evaluate one reachable input at one site, keeping the worst.
    fn feed(&mut self, site: &Site, input: [f64; 6]) {
        let written = (site.written)(&input);
        if !written.is_finite() {
            return;
        }
        self.tried += 1;
        let mut best_ulp = 0u64;
        let mut best_value = written;
        let mut best_name = site.variants[0].0;
        for (name, variant) in site.variants {
            let fused = variant(&input);
            if !fused.is_finite() {
                continue;
            }
            let gap = ulp_gap(written, fused);
            if gap >= best_ulp {
                best_ulp = gap;
                best_value = fused;
                best_name = name;
            }
        }
        if best_ulp > 0 {
            self.differing += 1;
        }
        if !self.seen || best_ulp > self.worst_ulp {
            self.seen = true;
            self.worst_ulp = best_ulp;
            self.worst_input = input;
            self.worst_written = written;
            self.worst_fused = best_value;
            self.worst_variant = best_name;
        }
    }
}

/// Pack a short argument list into the fixed six-slot input word.
#[inline]
fn word(v: &[f64]) -> [f64; 6] {
    let mut out = [0.0f64; 6];
    out[..v.len()].copy_from_slice(v);
    out
}

/// Format a six-slot input word as exact bits, `;`-joined so the CSV writer's
/// no-comma rule holds and the number on the row can be reproduced by hand.
fn format_word(w: &[f64; 6]) -> String {
    let mut out = String::with_capacity(6 * 19);
    for (i, v) in w.iter().enumerate() {
        if i > 0 {
            out.push(';');
        }
        let _ = write!(out, "{:#018x}", v.to_bits());
    }
    out
}

fn parse_word(text: &str) -> [f64; 6] {
    let mut out = [0.0f64; 6];
    for (slot, part) in out.iter_mut().zip(text.split(';')) {
        let hex = part.trim_start_matches("0x");
        *slot = f64::from_bits(u64::from_str_radix(hex, 16).expect("a hex bit pattern"));
    }
    out
}

// ── the site table ──────────────────────────────────────────────────────────
//
// The indices are named so the harvest can stream into a site without a lookup,
// and `check_site_indices` asserts every one still points at the site it is
// named for.

const S_PLACE: usize = 0;
const S_MC_CORNER: usize = 1;
const S_DUAL_PLACE: usize = 2;
const S_HERMITE: usize = 3;
const S_VEC3_DOT: usize = 4;
const S_VEC3_CROSS: usize = 5;
const S_ADJUGATE: usize = 6;
const S_DETERMINANT: usize = 7;
const S_DOT_EQUIVARIANT: usize = 8;
const S_TWO_PRODUCT: usize = 9;
const S_FADE: usize = 10;
const S_PERLIN_ACC: usize = 11;
const S_FBM_ACC: usize = 12;
const S_FBM_Q: usize = 13;
const S_TERRAIN: usize = 14;
const S_CONTROL_EXACT: usize = 15;
const S_CONTROL_SEPARATOR: usize = 16;

/// The plausible fusions of `a * b - c * d`. LLVM's FMA former can take either
/// product into the subtraction, so both are tried and the worse one is reported.
const CROSS_VARIANTS: Variants = &[
    ("fuse-first-product", |w| w[0].mul_add(w[1], -(w[2] * w[3]))),
    ("fuse-second-product", |w| {
        (-w[2]).mul_add(w[3], w[0] * w[1])
    }),
];

/// The plausible fusions of `acc + m * v`: there is exactly one product.
const ACCUMULATE_VARIANTS: Variants = &[("fuse-product", |w| w[1].mul_add(w[2], w[0]))];

/// The plausible fusions of `a0*b0 + a1*b1 + a2*b2`, which Rust associates left.
const DOT_VARIANTS: Variants = &[
    ("fuse-both-adds", |w| {
        w[4].mul_add(w[5], w[2].mul_add(w[3], w[0] * w[1]))
    }),
    ("fuse-last-add", |w| {
        w[4].mul_add(w[5], w[0] * w[1] + w[2] * w[3])
    }),
];

fn written_accumulate(w: &[f64; 6]) -> f64 {
    w[0] + w[1] * w[2]
}

fn written_cross(w: &[f64; 6]) -> f64 {
    w[0] * w[1] - w[2] * w[3]
}

fn written_dot(w: &[f64; 6]) -> f64 {
    w[0] * w[1] + w[2] * w[3] + w[4] * w[5]
}

fn reproduce_accumulate(w: &[f64; 6]) -> f64 {
    p117_probe_accumulate(w[0], w[1], w[2])
}

fn reproduce_cross(w: &[f64; 6]) -> f64 {
    p117_probe_cross(w[0], w[1], w[2], w[3])
}

fn reproduce_dot(w: &[f64; 6]) -> f64 {
    p117_probe_dot(w[0], w[1], w[2], w[3], w[4], w[5])
}

const SITES: &[Site] = &[
    Site {
        name: "cube::place",
        file: "crates/isomesh/src/cube.rs",
        line: 234,
        shape: "(lo + hi) * HALF + (hi - lo) * d",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_place",
        written: |w| (w[0] + w[1]) * 0.5 + (w[1] - w[0]) * w[2],
        variants: &[
            ("fuse-offset-term", |w| {
                (w[1] - w[0]).mul_add(w[2], (w[0] + w[1]) * 0.5)
            }),
            ("fuse-midpoint-term", |w| {
                (w[0] + w[1]).mul_add(0.5, (w[1] - w[0]) * w[2])
            }),
        ],
        reproduce: |w| p117_probe_place(w[0], w[1], w[2]),
    },
    Site {
        name: "marching_cubes::corner_position",
        file: "crates/isomesh/src/marching_cubes/mod.rs",
        line: 761,
        shape: "origin[0] + cell_size * R::from_f64(f64::from(base[0] + o[0]))",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_accumulate",
        written: written_accumulate,
        variants: ACCUMULATE_VARIANTS,
        reproduce: reproduce_accumulate,
    },
    Site {
        name: "dual::place_vertices",
        file: "crates/isomesh/src/dual.rs",
        line: 616,
        shape: "origin[0] + cell_size * R::from_f64(f64::from(x))",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_accumulate",
        written: written_accumulate,
        variants: ACCUMULATE_VARIANTS,
        reproduce: reproduce_accumulate,
    },
    Site {
        name: "hermite::HermiteCell::from_corners",
        file: "crates/isomesh/src/hermite.rs",
        line: 99,
        shape: "cell_origin[axis] + cell_size * place(from; to; d)",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_accumulate",
        written: written_accumulate,
        variants: ACCUMULATE_VARIANTS,
        reproduce: reproduce_accumulate,
    },
    Site {
        name: "vec3::dot",
        file: "crates/isomesh/src/vec3.rs",
        line: 29,
        shape: "a[0] * b[0] + a[1] * b[1] + a[2] * b[2]",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_dot",
        written: written_dot,
        variants: DOT_VARIANTS,
        reproduce: reproduce_dot,
    },
    Site {
        name: "vec3::cross",
        file: "crates/isomesh/src/vec3.rs",
        line: 33,
        shape: "a[1] * b[2] - a[2] * b[1]",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_cross",
        written: written_cross,
        variants: CROSS_VARIANTS,
        reproduce: reproduce_cross,
    },
    Site {
        name: "solve::Symmetric3::adjugate",
        file: "crates/isomesh/src/dual_contouring/solve.rs",
        line: 186,
        shape: "self.yy * self.zz - self.yz * self.yz",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_cross",
        written: written_cross,
        variants: CROSS_VARIANTS,
        reproduce: reproduce_cross,
    },
    Site {
        name: "solve::Symmetric3::determinant",
        file: "crates/isomesh/src/dual_contouring/solve.rs",
        line: 217,
        shape: "acc += two * mul_equivariant([xy; yz; xz])",
        reachable: true,
        on_path: true,
        product_also_compared: true,
        probe: "p117_probe_accumulate",
        written: written_accumulate,
        variants: ACCUMULATE_VARIANTS,
        reproduce: reproduce_accumulate,
    },
    Site {
        name: "solve::dot_equivariant",
        file: "crates/isomesh/src/dual_contouring/solve.rs",
        line: 104,
        shape: "sum_equivariant([a[0]*b[0]; a[1]*b[1]; a[2]*b[2]])",
        reachable: true,
        on_path: true,
        product_also_compared: true,
        probe: "p117_probe_dot_equivariant",
        written: |w| sum_equivariant(&mut [w[0] * w[1], w[2] * w[3], w[4] * w[5]]),
        variants: &[("fuse-into-sorted-sum", |w| {
            // The contracted form of the same reduction: sort the products by
            // magnitude exactly as `sum_equivariant` does, then accumulate each
            // one with an FMA of the pair that produced it.
            let mut pairs = [(w[0], w[1]), (w[2], w[3]), (w[4], w[5])];
            pairs.sort_by(|p, q| {
                let (a, b) = (p.0 * p.1, q.0 * q.1);
                if precedes(a, b) {
                    core::cmp::Ordering::Less
                } else if precedes(b, a) {
                    core::cmp::Ordering::Greater
                } else {
                    core::cmp::Ordering::Equal
                }
            });
            let mut acc = 0.0;
            for (a, b) in pairs {
                acc = a.mul_add(b, acc);
            }
            acc
        })],
        reproduce: |w| p117_probe_dot_equivariant(w[0], w[1], w[2], w[3], w[4], w[5]),
    },
    Site {
        name: "predicates::two_product",
        file: "crates/isomesh/src/predicates.rs",
        line: 104,
        shape: "y = (a_lo * b_lo) - err3 -- the whole Shewchuk roundoff",
        // `orient2d` and `incircle` are `pub` and reach it; no extractor does.
        reachable: true,
        on_path: false,
        product_also_compared: false,
        probe: "p117_probe_two_product_err",
        // The whole roundoff rather than one of its four subtractions, because
        // `predicates.rs:28` asks about the roundoff and not about a step: the
        // Dekker split makes `a_hi*b_hi`, `a_lo*b_hi` and `a_hi*b_lo` each
        // individually exact -- Shewchuk Theorem 17, and the stated reason the
        // split is there -- so fusing any of those three cannot change a bit.
        // `a_lo*b_lo` is the one product that can round, and it is the one the
        // return value subtracts.
        written: |w| two_product_roundoff(w[0], w[1]),
        variants: &[
            ("fuse-final-product", |w| {
                let (x, a_hi, a_lo, b_hi, b_lo) = two_product_parts(w[0], w[1]);
                let err1 = x - (a_hi * b_hi);
                let err2 = err1 - (a_lo * b_hi);
                let err3 = err2 - (a_hi * b_lo);
                a_lo.mul_add(b_lo, -err3)
            }),
            ("fuse-every-cross-product", |w| {
                let (x, a_hi, a_lo, b_hi, b_lo) = two_product_parts(w[0], w[1]);
                let err1 = (-a_hi).mul_add(b_hi, x);
                let err2 = (-a_lo).mul_add(b_hi, err1);
                let err3 = (-a_hi).mul_add(b_lo, err2);
                a_lo.mul_add(b_lo, -err3)
            }),
        ],
        reproduce: |w| p117_probe_two_product_err(w[0], w[1]),
    },
    Site {
        name: "noise::fade",
        file: "crates/isomesh/src/fields/noise.rs",
        line: 96,
        shape: "t3 * (t * (t * six - fifteen) + ten)",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_fade",
        written: |w| w[0] * (w[0] * 6.0 - 15.0) + 10.0,
        variants: &[
            ("fuse-outer", |w| w[0].mul_add(w[0] * 6.0 - 15.0, 10.0)),
            ("fuse-inner", |w| w[0] * w[0].mul_add(6.0, -15.0) + 10.0),
            ("fuse-both", |w| {
                w[0].mul_add(w[0].mul_add(6.0, -15.0), 10.0)
            }),
        ],
        reproduce: |w| p117_probe_fade(w[0]),
    },
    Site {
        name: "noise::perlin value accumulation",
        file: "crates/isomesh/src/fields/noise.rs",
        line: 168,
        shape: "value += w * dot",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_accumulate",
        written: written_accumulate,
        variants: ACCUMULATE_VARIANTS,
        reproduce: reproduce_accumulate,
    },
    Site {
        name: "noise::fbm octave accumulation",
        file: "crates/isomesh/src/fields/noise.rs",
        line: 243,
        shape: "value += amp * v",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_accumulate",
        written: written_accumulate,
        variants: ACCUMULATE_VARIANTS,
        reproduce: reproduce_accumulate,
    },
    Site {
        name: "noise::fbm lattice offset",
        file: "crates/isomesh/src/fields/noise.rs",
        line: 237,
        shape: "p[0] * freq + k * R::from_f64(OCTAVE_OFFSET[0])",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_two_products",
        written: |w| w[0] * w[1] + w[2] * w[3],
        variants: &[
            ("fuse-frequency", |w| w[0].mul_add(w[1], w[2] * w[3])),
            ("fuse-offset", |w| w[2].mul_add(w[3], w[0] * w[1])),
        ],
        reproduce: |w| p117_probe_two_products(w[0], w[1], w[2], w[3]),
    },
    Site {
        name: "FbmTerrain::sample",
        file: "crates/isomesh/src/fields/mod.rs",
        line: 1360,
        shape: "p[1] - (self.base_height + self.amplitude * n)",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_accumulate",
        written: written_accumulate,
        variants: ACCUMULATE_VARIANTS,
        reproduce: reproduce_accumulate,
    },
    Site {
        name: "control::exact_products",
        file: "crates/isomesh/src/fields/noise.rs",
        line: 145,
        shape: "gv[0] * d[0] + gv[1] * d[1] + gv[2] * d[2]",
        reachable: true,
        on_path: true,
        product_also_compared: false,
        probe: "p117_probe_dot",
        written: written_dot,
        variants: DOT_VARIANTS,
        reproduce: reproduce_dot,
    },
    Site {
        name: "control::searched_separator",
        file: "crates/isomesh/benches/experiment_p117.rs",
        line: 0,
        shape: "a * b + c on a SplitMix64-searched triple",
        reachable: false,
        on_path: false,
        product_also_compared: false,
        probe: "p117_probe_written",
        written: |w| w[0] * w[1] + w[2],
        variants: &[("fuse", |w| w[0].mul_add(w[1], w[2]))],
        reproduce: |w| p117_probe_written(w[0], w[1], w[2]),
    },
];

/// The named indices still point at the sites they are named for.
fn check_site_indices() {
    let expected = [
        (S_PLACE, "cube::place"),
        (S_MC_CORNER, "marching_cubes::corner_position"),
        (S_DUAL_PLACE, "dual::place_vertices"),
        (S_HERMITE, "hermite::HermiteCell::from_corners"),
        (S_VEC3_DOT, "vec3::dot"),
        (S_VEC3_CROSS, "vec3::cross"),
        (S_ADJUGATE, "solve::Symmetric3::adjugate"),
        (S_DETERMINANT, "solve::Symmetric3::determinant"),
        (S_DOT_EQUIVARIANT, "solve::dot_equivariant"),
        (S_TWO_PRODUCT, "predicates::two_product"),
        (S_FADE, "noise::fade"),
        (S_PERLIN_ACC, "noise::perlin value accumulation"),
        (S_FBM_ACC, "noise::fbm octave accumulation"),
        (S_FBM_Q, "noise::fbm lattice offset"),
        (S_TERRAIN, "FbmTerrain::sample"),
        (S_CONTROL_EXACT, "control::exact_products"),
        (S_CONTROL_SEPARATOR, "control::searched_separator"),
    ];
    assert_eq!(expected.len(), SITES.len(), "every site is named");
    for (index, name) in expected {
        assert_eq!(SITES[index].name, name, "site index {index}");
    }
}

// ── the harvest: inputs the crate actually reaches ──────────────────────────

/// The golden fixture's own three resolutions (`golden.rs:73`).
const RESOLUTIONS: [u32; 3] = [17, 25, 33];

/// One reading per site, fed by the harvest.
struct Readings(Vec<Reading>);

impl Readings {
    fn new() -> Self {
        Self(vec![Reading::default(); SITES.len()])
    }

    #[inline]
    fn feed(&mut self, index: usize, input: [f64; 6]) {
        self.0[index].feed(&SITES[index], input);
    }
}

/// Harvest every reachable input for every geometry site out of one field on one
/// grid.
///
/// This walks the same cells `extract` walks, on the same grid the golden
/// fixture uses, and hands each site the values that expression is actually
/// evaluated on. That is what makes `reachable` a fact rather than a claim: the
/// numbers on the rows below came out of the reference fields, not out of a
/// random generator.
fn harvest<F>(field: &F, samples: u32, readings: &mut Readings)
where
    F: Sdf<Scalar = f64> + ReferenceField,
{
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let n = samples as usize;

    // The whole grid, sampled once, on the same lattice `sample_grid` uses.
    let mut values = vec![0.0f64; n * n * n];
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let p = [
                    lo[0] + cell_size * (x as f64),
                    lo[1] + cell_size * (y as f64),
                    lo[2] + cell_size * (z as f64),
                ];
                values[x + n * (y + n * z)] = field.sample(p);
            }
        }
    }

    let cells = n - 1;
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let base = [x, y, z];
                let mut corner_value = [0.0f64; 8];
                for (corner, slot) in corner_value.iter_mut().enumerate() {
                    let o = corner_offset(corner as u8);
                    let (cx, cy, cz) = (x + o[0] as usize, y + o[1] as usize, z + o[2] as usize);
                    *slot = values[cx + n * (cy + n * cz)];
                }

                // `dual.rs:616` — the dual vertex's own cell index.
                for (axis, origin) in lo.iter().enumerate() {
                    readings.feed(S_DUAL_PLACE, word(&[*origin, cell_size, base[axis] as f64]));
                }

                let mut cut = false;
                for [c_lo, c_hi] in EDGE_CORNERS {
                    let (a, b) = (corner_value[c_lo as usize], corner_value[c_hi as usize]);
                    if is_inside(a) == is_inside(b) {
                        continue;
                    }
                    cut = true;
                    let d = edge_offset(a, b);

                    // `marching_cubes/mod.rs:761` — the two corner positions
                    // `edge_position` asks for, then `cube.rs:234` on each axis.
                    let o_lo = corner_offset(c_lo);
                    let o_hi = corner_offset(c_hi);
                    for (axis, origin) in lo.iter().enumerate() {
                        let i_lo = f64::from(base[axis] as u32 + o_lo[axis]);
                        let i_hi = f64::from(base[axis] as u32 + o_hi[axis]);
                        readings.feed(S_MC_CORNER, word(&[*origin, cell_size, i_lo]));
                        let p_lo = origin + cell_size * i_lo;
                        let p_hi = origin + cell_size * i_hi;
                        readings.feed(S_PLACE, word(&[p_lo, p_hi, d]));

                        // `hermite.rs:99` — the same crossing, reached through
                        // the cell-local `place` on `0`..`1` instead.
                        let from = f64::from(o_lo[axis]);
                        let to = f64::from(o_hi[axis]);
                        let cell_origin = origin + cell_size * f64::from(base[axis] as u32);
                        readings.feed(
                            S_HERMITE,
                            word(&[cell_origin, cell_size, place(from, to, d)]),
                        );
                    }
                }

                if !cut {
                    continue;
                }

                // A real Hermite cell, through the public constructor, so the
                // normals and crossing positions are the crate's own.
                let cell_origin = [
                    lo[0] + cell_size * (x as f64),
                    lo[1] + cell_size * (y as f64),
                    lo[2] + cell_size * (z as f64),
                ];
                let cell = HermiteCell::from_corners(field, &corner_value, cell_origin, cell_size);
                let Some(centroid) = cell.centroid() else {
                    continue;
                };

                let mut m_terms = [[0.0f64; EDGE_COUNT]; 6];
                let mut count = 0usize;
                let mut first_normal = [0.0f64; 3];
                let mut last_normal = [0.0f64; 3];
                for crossing in cell.iter() {
                    let nrm = crossing.normal;
                    // `vec3::dot` — the gradient normalisation every crossing
                    // runs, as `vec3::length_squared`'s `dot(a, a)`.
                    readings.feed(
                        S_VEC3_DOT,
                        word(&[nrm[0], nrm[0], nrm[1], nrm[1], nrm[2], nrm[2]]),
                    );
                    if count == 0 {
                        first_normal = nrm;
                    }
                    last_normal = nrm;

                    let rel = [
                        crossing.position[0] - centroid[0],
                        crossing.position[1] - centroid[1],
                        crossing.position[2] - centroid[2],
                    ];
                    // `solve.rs:289` — the plane distance, through the sorted dot.
                    readings.feed(
                        S_DOT_EQUIVARIANT,
                        word(&[nrm[0], rel[0], nrm[1], rel[1], nrm[2], rel[2]]),
                    );

                    let outer = [
                        nrm[0] * nrm[0],
                        nrm[0] * nrm[1],
                        nrm[0] * nrm[2],
                        nrm[1] * nrm[1],
                        nrm[1] * nrm[2],
                        nrm[2] * nrm[2],
                    ];
                    for (slot, value) in m_terms.iter_mut().zip(outer) {
                        slot[count] = value;
                    }
                    count += 1;
                }

                // `vec3::cross` — two crossing normals of one cell, the same 2×2
                // determinant `normals::area_weighted` evaluates on two triangle
                // edge vectors.
                if count >= 2 {
                    readings.feed(
                        S_VEC3_CROSS,
                        word(&[
                            first_normal[1],
                            last_normal[2],
                            first_normal[2],
                            last_normal[1],
                        ]),
                    );
                }

                // `solve.rs:186` and `:217` — the real regularized matrix.
                let mut e = [0.0f64; 6];
                for (slot, terms) in e.iter_mut().zip(m_terms.iter_mut()) {
                    *slot = sum_equivariant(&mut terms[..count]);
                }
                let (xx, xy, xz, yy, yz, zz) = (
                    e[0] + LAMBDA,
                    e[1],
                    e[2],
                    e[3] + LAMBDA,
                    e[4],
                    e[5] + LAMBDA,
                );
                readings.feed(S_ADJUGATE, word(&[yy, zz, yz, yz]));

                // The doubled term of the determinant, and the partial sum
                // standing at the moment `sum_equivariant` reaches it.
                let doubled_factor = mul_equivariant([xy, yz, xz]);
                let doubled = 2.0 * doubled_factor;
                let mut terms = [
                    mul_equivariant([xx, yy, zz]),
                    doubled,
                    -mul_equivariant([xx, yz, yz]),
                    -mul_equivariant([yy, xz, xz]),
                    -mul_equivariant([zz, xy, xy]),
                ];
                sort_by_magnitude(&mut terms);
                let mut acc = 0.0;
                for value in terms {
                    if value.to_bits() == doubled.to_bits() {
                        readings.feed(S_DETERMINANT, word(&[acc, 2.0, doubled_factor]));
                    }
                    acc += value;
                }

                // `predicates::two_product` on a real crossing coordinate pair.
                if let Some(crossing) = cell.iter().next() {
                    readings.feed(
                        S_TWO_PRODUCT,
                        word(&[crossing.position[0], crossing.position[2]]),
                    );
                }
            }
        }
    }
}

/// Harvest the noise sites by replaying `fbm`'s own loop on the real grid
/// positions `fbm_terrain` is sampled at, with that field's canonical
/// parameters read off `fields/mod.rs:1320-1328`.
fn harvest_noise(samples: u32, readings: &mut Readings) {
    let field = FbmTerrain::<f64>::canonical();
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let n = samples as usize;

    // `FbmTerrain::canonical`: octaves 4, lacunarity 2, gain 1/2, frequency 1/4,
    // amplitude 2, base_height 0.
    const OCTAVES: u32 = 4;
    const LACUNARITY: f64 = 2.0;
    const GAIN: f64 = 0.5;
    const FREQUENCY: f64 = 0.25;
    const AMPLITUDE: f64 = 2.0;
    const BASE_HEIGHT: f64 = 0.0;
    /// `noise.rs:206` — `octave × (1/φ, 1/φ², 1/φ³)`.
    const OCTAVE_OFFSET: [f64; 3] = [
        0.618_033_988_749_894_8,
        0.381_966_011_250_105_15,
        0.236_067_977_499_789_7,
    ];
    /// `noise.rs:44` — the twelve gradients, two `±1` components each.
    const GRAD12: [[i32; 3]; 12] = [
        [1, 1, 0],
        [-1, 1, 0],
        [1, -1, 0],
        [-1, -1, 0],
        [1, 0, 1],
        [-1, 0, 1],
        [1, 0, -1],
        [-1, 0, -1],
        [0, 1, 1],
        [0, -1, 1],
        [0, 1, -1],
        [0, -1, -1],
    ];

    for z in 0..n {
        for x in 0..n {
            // `FbmTerrain::sample` flattens `y` to zero before the noise.
            let p = [
                lo[0] + cell_size * (x as f64),
                0.0,
                lo[2] + cell_size * (z as f64),
            ];
            let mut value = 0.0f64;
            let mut freq = FREQUENCY;
            let mut amp = 1.0f64;
            for octave in 0..OCTAVES {
                let k = f64::from(octave);
                // `noise.rs:237` — two products either side of one add.
                for (axis, offset) in OCTAVE_OFFSET.iter().enumerate() {
                    readings.feed(S_FBM_Q, word(&[p[axis], freq, k, *offset]));
                }
                let q = [
                    p[0] * freq + k * OCTAVE_OFFSET[0],
                    p[1] * freq + k * OCTAVE_OFFSET[1],
                    p[2] * freq + k * OCTAVE_OFFSET[2],
                ];

                // `noise.rs:96` — the quintic fade on each fractional coordinate.
                let t = [
                    q[0] - q[0].floor(),
                    q[1] - q[1].floor(),
                    q[2] - q[2].floor(),
                ];
                let mut fade_u = [0.0f64; 3];
                for (slot, ti) in fade_u.iter_mut().zip(t) {
                    readings.feed(S_FADE, word(&[ti]));
                    *slot = (ti * ti * ti) * (ti * (ti * 6.0 - 15.0) + 10.0);
                }

                // `noise.rs:145` and `:168` — the gradient dot with its exactly
                // representable products, and the weighted accumulation.
                let mut v = 0.0f64;
                for corner in 0..8u32 {
                    let (cx, cy, cz) = (corner & 1, (corner >> 1) & 1, (corner >> 2) & 1);
                    let g = GRAD12[(corner as usize * 5 + x + z) % 12];
                    let gv = [f64::from(g[0]), f64::from(g[1]), f64::from(g[2])];
                    let dv = [
                        t[0] - f64::from(cx),
                        t[1] - f64::from(cy),
                        t[2] - f64::from(cz),
                    ];
                    readings.feed(
                        S_CONTROL_EXACT,
                        word(&[gv[0], dv[0], gv[1], dv[1], gv[2], dv[2]]),
                    );
                    let dot = gv[0] * dv[0] + gv[1] * dv[1] + gv[2] * dv[2];
                    let wx = if cx == 1 { fade_u[0] } else { 1.0 - fade_u[0] };
                    let wy = if cy == 1 { fade_u[1] } else { 1.0 - fade_u[1] };
                    let wz = if cz == 1 { fade_u[2] } else { 1.0 - fade_u[2] };
                    let weight = wx * wy * wz;
                    readings.feed(S_PERLIN_ACC, word(&[v, weight, dot]));
                    v += weight * dot;
                }

                readings.feed(S_FBM_ACC, word(&[value, amp, v]));
                value += amp * v;
                freq *= LACUNARITY;
                amp *= GAIN;
            }
            readings.feed(S_TERRAIN, word(&[BASE_HEIGHT, AMPLITUDE, value]));
        }
    }
}

/// The searched control: a triple on which the two forms **must** differ.
///
/// `cube.rs:349-356`'s precedent and its stated reason — the first version of
/// that test picked its separating pair by hand and the values turned out to be
/// ones where both forms agree exactly. So this searches, deterministically,
/// with the same SplitMix64 the crate's own sweep uses.
fn harvest_searched_separator(readings: &mut Readings) {
    let mut state = 0x2026_u64;
    let mut next = || {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        (z >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    };
    for _ in 0..4096 {
        let (a, b, c) = (next(), next(), next());
        readings.feed(S_CONTROL_SEPARATOR, word(&[a, b, c]));
    }
}

// ── the mul_add census ──────────────────────────────────────────────────────

/// What `grep mul_add crates/isomesh/src/` finds, split into calls and prose.
struct Census {
    call_sites: usize,
    mentions: usize,
    /// `<path>:<lines mentioning it>`, for the run's own log.
    where_mentioned: Vec<String>,
}

fn census(root: &Path, src: &Path) -> Census {
    fn walk(dir: &Path, out: &mut Vec<PathBuf>) {
        let Ok(entries) = std::fs::read_dir(dir) else {
            return;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                walk(&path, out);
            } else if path.extension().is_some_and(|e| e == "rs") {
                out.push(path);
            }
        }
    }
    let mut files = Vec::new();
    walk(src, &mut files);
    files.sort();

    let mut call_sites = 0usize;
    let mut mentions = 0usize;
    let mut where_mentioned = Vec::new();
    for file in &files {
        let Ok(text) = std::fs::read_to_string(file) else {
            continue;
        };
        let calls = text.matches(".mul_add(").count();
        let lines = text.lines().filter(|l| l.contains("mul_add")).count();
        call_sites += calls;
        mentions += lines;
        if lines > 0 {
            let short = file.strip_prefix(root).unwrap_or(file);
            where_mentioned.push(format!(
                "{} {lines} mention(s) {calls} call(s)",
                short.display()
            ));
        }
    }
    Census {
        call_sites,
        mentions,
        where_mentioned,
    }
}

// ── the disassembly arm ─────────────────────────────────────────────────────

/// The disassembled body of every symbol in this binary, keyed by symbol name.
fn disassemble() -> BTreeMap<String, String> {
    let exe = std::env::current_exe().expect("this binary has a path");
    let out = Command::new("objdump")
        .arg("-d")
        .arg("--no-show-raw-insn")
        .arg(&exe)
        .output()
        .unwrap_or_else(|e| {
            panic!(
                "P-117 refuses to run without a disassembler: `objdump -d {}` failed ({e}). \
                 The registered vacuity control is that a KNOWN-FUSED probe is seen fusing, \
                 and a harness that cannot read machine code cannot claim to have looked.",
                exe.display()
            )
        });
    assert!(
        out.status.success(),
        "P-117: objdump exited {:?} on {}",
        out.status.code(),
        exe.display()
    );
    let text = String::from_utf8_lossy(&out.stdout);

    let mut symbols = BTreeMap::new();
    let mut current: Option<(String, String)> = None;
    for line in text.lines() {
        if let Some(open) = line.find(" <") {
            let address = &line[..open];
            if line.ends_with(">:")
                && !address.is_empty()
                && address.bytes().all(|b| b.is_ascii_hexdigit())
            {
                if let Some((name, body)) = current.take() {
                    symbols.insert(name, body);
                }
                current = Some((line[open + 2..line.len() - 2].to_string(), String::new()));
                continue;
            }
        }
        if let Some((_, body)) = current.as_mut() {
            body.push_str(line);
            body.push('\n');
        }
    }
    if let Some((name, body)) = current.take() {
        symbols.insert(name, body);
    }
    assert!(
        symbols.len() > 100,
        "P-117: objdump produced {} symbols, which is too few to be this binary",
        symbols.len()
    );
    symbols
}

/// Is this mnemonic a hardware fused multiply-add?
///
/// x86-64: FMA3's `vfmadd`/`vfmsub`/`vfnmadd`/`vfnmsub` families, in every
/// operand ordering (`132`/`213`/`231`) and both widths. aarch64: the scalar
/// `fmadd`/`fmsub`/`fnmadd`/`fnmsub` and the vector `fmla`/`fmls`.
fn is_fma_mnemonic(m: &str) -> bool {
    m.starts_with("vfmadd")
        || m.starts_with("vfmsub")
        || m.starts_with("vfnmadd")
        || m.starts_with("vfnmsub")
        || m == "fmadd"
        || m == "fmsub"
        || m == "fnmadd"
        || m == "fnmsub"
        || m.starts_with("fmla")
        || m.starts_with("fmls")
}

/// The mnemonics of one disassembled body.
fn mnemonics(body: &str) -> impl Iterator<Item = &str> {
    body.lines().filter_map(|line| {
        let (_, rest) = line.split_once(':')?;
        rest.split_whitespace().next()
    })
}

/// How many hardware fused multiply-adds one body carries, and the first
/// mnemonic that decided it.
fn fusion(body: &str) -> (usize, String) {
    let mut count = 0usize;
    let mut evidence = String::new();
    for m in mnemonics(body) {
        if is_fma_mnemonic(m) {
            count += 1;
            if evidence.is_empty() {
                evidence = m.to_string();
            }
        }
    }
    if evidence.is_empty() {
        evidence = String::from("none");
    }
    (count, evidence)
}

/// Does this body do floating-point multiply/add arithmetic of its own, as
/// opposed to delegating through a call or a tail jump?
fn does_arithmetic(body: &str) -> bool {
    mnemonics(body).any(|m| {
        let m = m.strip_prefix('v').unwrap_or(m);
        m.starts_with("mul")
            || m.starts_with("add")
            || m.starts_with("sub")
            || m.starts_with("fmul")
    })
}

/// Find the probe symbol named exactly `needle`.
///
/// Both manglings in play here — rustc's `v0` on ELF and the legacy `_ZN…` on
/// Mach-O — length-prefix every path segment, so searching for
/// `"{len}{needle}"` is an exact segment match rather than a substring one. A
/// bare `contains` is ambiguous in a way that matters: `p117_probe_dot` is a
/// substring of `p117_probe_dot_equivariant`, and `p117_probe_two_product_err`
/// contains `p117_probe_two_product`. Which of the two a `contains` walk found
/// depended on the disambiguator hash sorting first, which is a different answer
/// on a different build of the same source.
fn symbol_body<'a>(symbols: &'a BTreeMap<String, String>, needle: &str) -> Option<&'a str> {
    let segment = format!("{}{needle}", needle.len());
    symbols
        .iter()
        .find(|(name, _)| name.contains(&segment))
        .map(|(_, body)| body.as_str())
}

/// What the disassembly arm found.
struct Codegen {
    known_fused_count: usize,
    known_fused_evidence: String,
    known_fused_symbol: String,
    /// Per probe: hardware FMA count in that probe's own body.
    probes: BTreeMap<&'static str, (usize, String)>,
    isomesh_symbols: usize,
    isomesh_fma: usize,
    total_fma: usize,
    extract_symbols: usize,
    extract_fma: usize,
}

fn read_codegen() -> Codegen {
    let symbols = disassemble();

    let known = symbol_body(&symbols, "p117_probe_known_fused").unwrap_or_else(|| {
        panic!(
            "P-117: the known-fused probe has no symbol in the disassembly of this binary. \
             It is `#[inline(never)]` and is called through `black_box`, so its absence \
             means the reader is looking at the wrong object -- and a reader that cannot \
             find the probe that must fuse cannot report on the shapes that must not."
        )
    });

    // The probe's own body first. On aarch64 `mul_add` is the baseline `fmadd`
    // and lands here. On x86-64 the FMA3 instruction is **not** in the baseline
    // ISA, so rustc lowers `mul_add` to a tail call into `compiler_builtins`'
    // runtime-dispatched `fma`, whose hardware variant is a real `vfmadd` --
    // still fusion, one rounding where the written form has two. Following that
    // one level is what keeps the vacuity control honest on a machine whose
    // baseline has no FMA, instead of reporting a false negative that would then
    // excuse every other zero on the sheet.
    let (mut count, mut evidence) = fusion(known);
    let mut symbol = String::from("p117_probe_known_fused");
    if count == 0 && !does_arithmetic(known) {
        // `symbols` is a `BTreeMap`, so this walk is in symbol order and the
        // choice is the same on every run of the same binary. It is deliberately
        // the FIRST fma-named symbol that carries a hardware FMA rather than the
        // one with the most: `max_by_key` broke ties arbitrarily and picked
        // `fma_with_fma4` -- the AMD FMA4 variant, which a Zen 3 does not have
        // and the runtime dispatch would never select -- over `fma_with_fma`,
        // the FMA3 one it actually takes. Both counted 1, so the count was right
        // and the mnemonic on the row named an instruction this CPU cannot run.
        let dispatched = symbols
            .iter()
            .filter(|(name, _)| name.contains("fma"))
            .map(|(name, body)| {
                let (c, e) = fusion(body);
                (name.as_str(), c, e)
            })
            .find(|(_, c, _)| *c > 0);
        if let Some((name, c, e)) = dispatched {
            count = c;
            evidence = format!("dispatch-to-fma:{e}");
            symbol = name.chars().take(120).collect();
        }
    }

    let mut probes = BTreeMap::new();
    for probe in PROBES {
        let reading = symbol_body(&symbols, probe)
            .map_or_else(|| (0usize, String::from("symbol-not-found")), fusion);
        probes.insert(probe, reading);
    }

    let mut isomesh_symbols = 0usize;
    let mut isomesh_fma = 0usize;
    let mut total_fma = 0usize;
    let mut extract_symbols = 0usize;
    let mut extract_fma = 0usize;
    for (name, body) in &symbols {
        let count = mnemonics(body).filter(|m| is_fma_mnemonic(m)).count();
        total_fma += count;
        if name.contains("isomesh") {
            isomesh_symbols += 1;
            isomesh_fma += count;
            if name.contains("extract") {
                extract_symbols += 1;
                extract_fma += count;
            }
        }
    }

    Codegen {
        known_fused_count: count,
        known_fused_evidence: evidence,
        known_fused_symbol: symbol,
        probes,
        isomesh_symbols,
        isomesh_fma,
        total_fma,
        extract_symbols,
        extract_fma,
    }
}

// ── the 216-hash fixture ────────────────────────────────────────────────────

/// The nine extractors, by the name `golden.rs:148` gives them.
const ALGORITHMS: [&str; 9] = [
    "greedy_quads",
    "marching_cubes",
    "marching_cubes+decider",
    "marching_cubes+trilinear",
    "marching_tetrahedra",
    "surface_nets",
    "dual_contouring",
    "manifold_dual_contouring",
    "subgrid_marching_tetrahedra",
];

/// `golden.rs:119` — the subgrid fixture's nailed-down sampling resolution.
const SUBGRID_SAMPLES: u32 = 16;

/// One row of the fixture: the key, and what the mesh hashed to.
#[derive(Clone, PartialEq, Eq)]
struct Golden {
    key: String,
    vertices: usize,
    triangles: usize,
    hash: String,
}

/// `golden.rs:159`, replayed through the public extractors.
fn golden_extract<F>(algorithm: &str, field: &F, samples: u32) -> MeshBuffer<f64>
where
    F: Sdf<Scalar = f64> + ReferenceField,
{
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let mut out = MeshBuffer::<f64>::new();
    match algorithm {
        "greedy_quads" => GreedyQuads::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        "marching_cubes" => MarchingCubes::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        "marching_cubes+decider" => {
            let mut mc = MarchingCubes::<f64>::new();
            mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
            mc.extract(field, &shape, lo, cell_size, &mut out)
                .expect("extraction");
        }
        "marching_cubes+trilinear" => {
            let mut mc = MarchingCubes::<f64>::new();
            mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
            mc.set_interior_ambiguity(InteriorAmbiguity::Trilinear);
            mc.extract(field, &shape, lo, cell_size, &mut out)
                .expect("extraction");
        }
        "marching_tetrahedra" => MarchingTetrahedra::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        "surface_nets" => SurfaceNets::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        "dual_contouring" => DualContouring::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        "manifold_dual_contouring" => ManifoldDualContouring::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        "subgrid_marching_tetrahedra" => {
            SubgridMarchingTetrahedra::<f64>::new(SUBGRID_SAMPLES)
                .expect("a positive sampling resolution")
                .extract(field, &shape, lo, cell_size, &mut out)
                .expect("extraction");
        }
        other => panic!("P-117: no such golden algorithm: {other}"),
    }
    out
}

/// All 216, in `golden.rs:211`'s order: field, then algorithm, then resolution.
fn compute_golden() -> Vec<Golden> {
    let mut entries = Vec::new();
    isomesh::for_each_reference_field!(f64, |name, field| {
        for algorithm in ALGORITHMS {
            for samples in RESOLUTIONS {
                let mesh = golden_extract(algorithm, &field, samples);
                entries.push(Golden {
                    key: format!("{algorithm}/{name}/{samples}"),
                    vertices: mesh.vertex_count(),
                    triangles: mesh.triangle_count(),
                    hash: format!("{:016x}", mesh_hash(&mesh)),
                });
            }
        }
    });
    entries
}

/// The committed fixture, read with `golden.rs:253`'s one-line scanner.
fn read_committed_golden(root: &Path) -> Vec<Golden> {
    let path = root.join("crates/isomesh/golden_hashes.json");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("P-117: no golden fixture at {}: {e}", path.display()));

    fn field_of<'a>(line: &'a str, key: &str) -> Option<&'a str> {
        let at = line.find(&format!("\"{key}\":"))? + key.len() + 3;
        let rest = &line[at..];
        let rest = rest.strip_prefix('"').unwrap_or(rest);
        let end = rest.find(['"', ',', '}'])?;
        Some(&rest[..end])
    }

    let mut out = Vec::new();
    for line in text.lines() {
        let Some(algorithm) = field_of(line, "algorithm") else {
            continue;
        };
        let field = field_of(line, "field").expect("field");
        let samples = field_of(line, "samples").expect("samples");
        out.push(Golden {
            key: format!("{algorithm}/{field}/{samples}"),
            vertices: field_of(line, "vertices")
                .expect("vertices")
                .parse()
                .expect("count"),
            triangles: field_of(line, "triangles")
                .expect("triangles")
                .parse()
                .expect("count"),
            hash: field_of(line, "hash").expect("hash").to_string(),
        });
    }
    out
}

// ── the peer artefact ───────────────────────────────────────────────────────

/// One site as the peer machine measured it.
#[derive(Clone)]
struct PeerSite {
    tried: u64,
    differing: u64,
    ulp: u64,
    written_bits: String,
    fused_bits: String,
    variant: String,
    reproduced: bool,
    worst_input: String,
}

/// Everything the aarch64 arm measured, read back from the committed input.
struct Peer {
    target: String,
    machine: String,
    rustc: String,
    mul_add_call_sites: usize,
    mul_add_mentions: usize,
    known_fused_count: usize,
    known_fused_evidence: String,
    known_fused_symbol: String,
    isomesh_symbols: usize,
    isomesh_fma: usize,
    total_fma: usize,
    extract_symbols: usize,
    extract_fma: usize,
    probes: BTreeMap<String, usize>,
    sites: BTreeMap<String, PeerSite>,
    golden: Vec<Golden>,
}

fn parse_peer(text: &str) -> Peer {
    let mut kv: BTreeMap<String, String> = BTreeMap::new();
    let mut probes = BTreeMap::new();
    let mut sites = BTreeMap::new();
    let mut golden = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let Some(head) = parts.next() else { continue };
        match head {
            "probe" => {
                let name = parts.next().expect("probe name").to_string();
                let count: usize = parts.next().expect("probe count").parse().expect("count");
                probes.insert(name, count);
            }
            "site" => {
                let name = parts.next().expect("site name").replace('~', " ");
                sites.insert(
                    name,
                    PeerSite {
                        tried: parts.next().expect("tried").parse().expect("u64"),
                        differing: parts.next().expect("differing").parse().expect("u64"),
                        ulp: parts.next().expect("ulp").parse().expect("u64"),
                        written_bits: parts.next().expect("written").to_string(),
                        fused_bits: parts.next().expect("fused").to_string(),
                        variant: parts.next().expect("variant").to_string(),
                        reproduced: parts.next().expect("reproduced") == "true",
                        worst_input: parts.next().expect("worst_input").to_string(),
                    },
                );
            }
            "golden" => {
                let algorithm = parts.next().expect("algorithm");
                let field = parts.next().expect("field");
                let samples = parts.next().expect("samples");
                golden.push(Golden {
                    key: format!("{algorithm}/{field}/{samples}"),
                    vertices: parts.next().expect("vertices").parse().expect("count"),
                    triangles: parts.next().expect("triangles").parse().expect("count"),
                    hash: parts.next().expect("hash").to_string(),
                });
            }
            other => {
                if let Some((key, value)) = other.split_once('=') {
                    kv.insert(key.to_string(), value.to_string());
                }
            }
        }
    }
    let get = |k: &str| {
        kv.get(k)
            .cloned()
            .unwrap_or_else(|| String::from("unknown"))
    };
    let num = |k: &str| kv.get(k).and_then(|v| v.parse().ok()).unwrap_or(0);
    Peer {
        target: get("target"),
        machine: get("machine"),
        rustc: get("rustc"),
        mul_add_call_sites: num("mul_add_call_sites"),
        mul_add_mentions: num("mul_add_mentions"),
        known_fused_count: num("fma_in_known_fused_probe"),
        known_fused_evidence: get("known_fused_probe_evidence"),
        known_fused_symbol: get("known_fused_dispatch_symbol"),
        isomesh_symbols: num("isomesh_symbols_scanned"),
        isomesh_fma: num("fma_instructions_isomesh_text"),
        total_fma: num("fma_instructions_binary_total"),
        extract_symbols: num("isomesh_extract_symbols_scanned"),
        extract_fma: num("fma_instructions_isomesh_extract"),
        probes,
        sites,
        golden,
    }
}

/// Write this machine's arm out as the peer input the other machine reads.
fn write_peer_artefact(
    root: &Path,
    census: &Census,
    codegen: &Codegen,
    readings: &Readings,
    golden: &[Golden],
) {
    let ask = |program: &str, args: &[&str]| -> String {
        Command::new(program)
            .args(args)
            .output()
            .ok()
            .filter(|o| o.status.success())
            .and_then(|o| String::from_utf8(o.stdout).ok())
            .and_then(|s| s.lines().next().map(|l| l.trim().replace(' ', "~")))
            .unwrap_or_else(|| String::from("unknown"))
    };
    let machine = Command::new(root.join("scripts/machine.sh"))
        .arg("--slug")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| String::from("unknown"));

    let mut out = String::new();
    out.push_str(
        "# P-117 C2 cross-machine fixture, and the aarch64 half of the disassembly arm.\n\
         #\n\
         # An INPUT, not a result. Produced by\n\
         #   ISOMESH_P117_PEER=1 cargo bench -p isomesh --bench experiment_p117\n\
         # on the peer machine, and read back by the same bench on the other target to fill\n\
         # `hashes_identical`, the aarch64 `fmadd_in_known_fused_probe` and every aarch64\n\
         # row. Without it C2 is recorded BLOCKED with the blocker named on the row -- never\n\
         # a zero-filled peer table, which is what R-081 shipped and had deleted.\n\
         #\n\
         # THREE THINGS ARE IN HERE. (1) The 216 golden hashes recomputed on this machine\n\
         # from the same source revision -- 8 reference fields x 9 algorithms x 3 resolutions\n\
         # {17, 25, 33}, `golden.rs:73`, `:122`, `fields/mod.rs:212`. (2) The disassembly\n\
         # reading of this machine's own bench binary, including the KNOWN-FUSED probe the\n\
         # registered vacuity control requires to be seen fusing on aarch64. (3) The\n\
         # per-site ULP gap between each multiply-add shape as written and the same shape\n\
         # fused, on inputs harvested from real extraction on this machine.\n\
         #\n\
         # `~` stands in for a space: this file is whitespace-delimited.\n\
         #\n",
    );
    let _ = writeln!(out, "target={HOST_TARGET}");
    let _ = writeln!(out, "machine={machine}");
    let _ = writeln!(out, "rustc={}", ask("rustc", &["--version"]));
    let _ = writeln!(out, "objdump={}", ask("objdump", &["--version"]));
    let _ = writeln!(out, "mul_add_call_sites={}", census.call_sites);
    let _ = writeln!(out, "mul_add_mentions={}", census.mentions);
    let _ = writeln!(out, "golden_rows={}", golden.len());
    let _ = writeln!(
        out,
        "fma_in_known_fused_probe={}",
        codegen.known_fused_count
    );
    let _ = writeln!(
        out,
        "known_fused_probe_evidence={}",
        codegen.known_fused_evidence
    );
    let _ = writeln!(
        out,
        "known_fused_dispatch_symbol={}",
        codegen.known_fused_symbol
    );
    let _ = writeln!(out, "fma_instructions_binary_total={}", codegen.total_fma);
    let _ = writeln!(out, "isomesh_symbols_scanned={}", codegen.isomesh_symbols);
    let _ = writeln!(out, "fma_instructions_isomesh_text={}", codegen.isomesh_fma);
    let _ = writeln!(
        out,
        "isomesh_extract_symbols_scanned={}",
        codegen.extract_symbols
    );
    let _ = writeln!(
        out,
        "fma_instructions_isomesh_extract={}",
        codegen.extract_fma
    );

    out.push_str("#\n# probe <name> <hardware_fma_instructions> <evidence>\n");
    for (name, (count, evidence)) in &codegen.probes {
        let _ = writeln!(out, "probe {name} {count} {evidence}");
    }

    out.push_str(
        "#\n# site <name> <tried> <differing> <max_ulp> <written_bits> <fused_bits> \
         <variant> <probe_reproduced> <worst_input>\n",
    );
    for (site, reading) in SITES.iter().zip(&readings.0) {
        let reproduced =
            (site.reproduce)(&reading.worst_input).to_bits() == reading.worst_written.to_bits();
        let _ = writeln!(
            out,
            "site {} {} {} {} {:#018x} {:#018x} {} {reproduced} {}",
            site.name.replace(' ', "~"),
            reading.tried,
            reading.differing,
            reading.worst_ulp,
            reading.worst_written.to_bits(),
            reading.worst_fused.to_bits(),
            reading.worst_variant,
            format_word(&reading.worst_input),
        );
    }

    out.push_str("#\n# golden <algorithm> <field> <samples> <vertices> <triangles> <hash>\n");
    for entry in golden {
        let _ = writeln!(
            out,
            "golden {} {} {} {}",
            entry.key.replace('/', " "),
            entry.vertices,
            entry.triangles,
            entry.hash
        );
    }

    let path = root.join(PEER_ARTEFACT);
    let _ = std::fs::create_dir_all(path.parent().expect("a parent"));
    std::fs::write(&path, &out).expect("write the peer artefact");
    println!("peer artefact → {}", path.display());
}

// ── main ────────────────────────────────────────────────────────────────────

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    assert_ne!(
        HOST_TARGET, "unregistered",
        "P-117's fixture is {X86} and {ARM}. This machine is neither, and a third \
         platform's numbers filed under a `target` column the registration never names \
         is how a fixture quietly becomes a different fixture."
    );
    check_site_indices();

    let root = repo_root();
    let src = root.join("crates/isomesh/src");

    println!("P-117: the mul_add census over crates/isomesh/src");
    let census = census(&root, &src);
    println!(
        "  `.mul_add(` call sites: {}   lines mentioning mul_add: {}",
        census.call_sites, census.mentions
    );
    for entry in &census.where_mentioned {
        println!("    {entry}");
    }

    println!("\nP-117: harvesting reachable inputs from the eight reference fields");
    let mut readings = Readings::new();
    isomesh::for_each_reference_field!(f64, |name, field| {
        for samples in RESOLUTIONS {
            harvest(&field, samples, &mut readings);
        }
        println!("  harvested {name}");
    });
    for samples in RESOLUTIONS {
        harvest_noise(samples, &mut readings);
    }
    harvest_searched_separator(&mut readings);

    // Keep every probe alive and honest: `#[inline(never)]` stops the inliner,
    // `black_box` stops constant folding, and the results are consumed so
    // nothing is dead.
    let mut alive = 0.0f64;
    for (site, reading) in SITES.iter().zip(&readings.0) {
        let w = std::hint::black_box(reading.worst_input);
        alive += (site.reproduce)(&w);
    }
    alive += p117_probe_known_fused(
        std::hint::black_box(1.0),
        std::hint::black_box(3.0),
        std::hint::black_box(7.0),
    );
    alive += p117_probe_written(
        std::hint::black_box(1.0),
        std::hint::black_box(3.0),
        std::hint::black_box(7.0),
    );
    assert!(alive.is_finite(), "every probe ran");

    println!("\nP-117: reading this binary's own machine code");
    let codegen = read_codegen();
    println!(
        "  known-fused probe: {} fused instruction(s), evidence `{}` in `{}`",
        codegen.known_fused_count, codegen.known_fused_evidence, codegen.known_fused_symbol
    );
    for (name, (count, evidence)) in &codegen.probes {
        println!("  {name:32} {count} {evidence}");
    }
    println!(
        "  {} isomesh symbols carrying {} fused instructions; {} extractor symbols carrying \
         {}; {} in the whole binary",
        codegen.isomesh_symbols,
        codegen.isomesh_fma,
        codegen.extract_symbols,
        codegen.extract_fma,
        codegen.total_fma
    );

    println!("\nP-117: recomputing the 216 golden hashes on {HOST_TARGET}");
    let golden = compute_golden();
    println!("  {} rows", golden.len());

    if std::env::var_os("ISOMESH_P117_PEER").is_some() {
        write_peer_artefact(&root, &census, &codegen, &readings, &golden);
        return;
    }

    let committed = read_committed_golden(&root);
    assert_eq!(
        committed.len(),
        golden.len(),
        "P-117: the committed fixture has {} rows and this run computed {}",
        committed.len(),
        golden.len()
    );
    let local_matches_committed = committed == golden;

    let peer = std::fs::read_to_string(root.join(PEER_ARTEFACT))
        .ok()
        .map(|t| parse_peer(&t));

    // C1 is scored EXACTLY as the clause is worded: at least one expression in
    // `crates/isomesh/src/` whose fused and unfused evaluations differ in the
    // last bit on an input the crate can reach. Whether the compiler *takes* the
    // fusion is a different question, it is the falsifier's stated reason, and
    // it is reported beside the verdict in `contracted_in_codegen` and
    // `fma_instructions_isomesh_extract` rather than folded into it.
    let sensitive: Vec<&str> = SITES
        .iter()
        .zip(&readings.0)
        .filter(|(site, r)| r.worst_ulp > 0 && site.file.starts_with("crates/isomesh/src/"))
        .map(|(site, _)| site.name)
        .collect();
    let c1 = !sensitive.is_empty();

    // C3: every divergence **in the crate** names a file and a line and is
    // reproduced by the isolated `#[inline(never)]` probe carrying that
    // expression verbatim. The filter is `crates/isomesh/src/`, the same one C1
    // uses, and it is load-bearing: `control::searched_separator` is a
    // bench-local control with `line: 0` by construction, so scoring C3 over
    // every row scored it FALSIFIED on the presence of its own control.
    let mut c3 = true;
    let mut unattributed: Vec<&str> = Vec::new();
    for (site, reading) in SITES.iter().zip(&readings.0) {
        if reading.worst_ulp == 0 || !site.file.starts_with("crates/isomesh/src/") {
            continue;
        }
        let reproduced =
            (site.reproduce)(&reading.worst_input).to_bits() == reading.worst_written.to_bits();
        if !reproduced || site.line == 0 {
            c3 = false;
            unattributed.push(site.name);
        }
    }
    if let Some(p) = &peer {
        for site in SITES {
            if !site.file.starts_with("crates/isomesh/src/") {
                continue;
            }
            let Some(s) = p.sites.get(site.name) else {
                continue;
            };
            if s.ulp > 0 && !s.reproduced {
                c3 = false;
                unattributed.push(site.name);
            }
        }
    }

    // The vacuity controls, in both directions, before a row is written.
    let exact_control = readings.0[S_CONTROL_EXACT].worst_ulp;
    let separator_control = readings.0[S_CONTROL_SEPARATOR].worst_ulp;
    assert!(
        readings.0[S_CONTROL_EXACT].tried > 0,
        "P-117: the exact-products control saw no input"
    );
    assert_eq!(
        exact_control, 0,
        "P-117 vacuity control: `control::exact_products` read {exact_control} ULP. Every \
         product there is exact -- the Perlin gradients have components of exactly 0 and \
         +/-1 -- so fusing cannot change the result. A non-zero means the ULP comparator \
         cries wolf and every other row on this sheet is suspect."
    );
    assert!(
        separator_control > 0,
        "P-117 vacuity control: `control::searched_separator` read 0 ULP over 4096 \
         SplitMix64 triples. The comparator cannot see a difference that must be there, \
         so a zero anywhere else on this sheet means nothing."
    );
    assert!(
        codegen.known_fused_count > 0,
        "P-117 vacuity control: the known-fused probe disassembled to no fused \
         instruction, in its own body or one level of dispatch away. A reader that cannot \
         see an FMA where one must exist cannot claim to see its absence elsewhere."
    );
    if let Some(p) = &peer {
        assert!(
            p.known_fused_count > 0,
            "P-117 vacuity control: the peer artefact reports {} fused instructions in the \
             known-fused probe on {}. The registration requires the control to fire on the \
             aarch64 arm specifically.",
            p.known_fused_count,
            p.target
        );
    }

    let prereg = isomesh::experiment!("P-117");
    common::experiment::run(prereg, |run| {
        let golden_rows = golden.len();

        let (c2, cross_identical, differing, blocker) = match &peer {
            Some(p) => {
                assert_eq!(
                    p.target, ARM,
                    "P-117: the peer artefact says target={} and the registration's peer \
                     arm is {ARM}",
                    p.target
                );
                assert_eq!(
                    p.golden.len(),
                    golden_rows,
                    "P-117: the peer measured {} golden rows and this machine {golden_rows}",
                    p.golden.len()
                );
                let differing = golden
                    .iter()
                    .zip(&p.golden)
                    .filter(|(a, b)| a.key != b.key || a.hash != b.hash)
                    .count();
                let identical = differing == 0;
                (
                    String::from(if identical { "true" } else { "false" }),
                    identical,
                    differing,
                    String::from("-"),
                )
            }
            None => (
                String::from("BLOCKED"),
                false,
                0,
                String::from("no-peer-artefact-committed"),
            ),
        };
        let peer_matches_committed = peer.as_ref().is_some_and(|p| p.golden == committed);

        for target in [X86, ARM] {
            let on_peer = target == ARM;
            if on_peer && peer.is_none() {
                continue;
            }
            let p = peer.as_ref();
            let (known_fused_count, known_fused_evidence, known_fused_symbol) = if on_peer {
                let p = p.expect("peer");
                (
                    p.known_fused_count,
                    p.known_fused_evidence.clone(),
                    p.known_fused_symbol.clone(),
                )
            } else {
                (
                    codegen.known_fused_count,
                    codegen.known_fused_evidence.clone(),
                    codegen.known_fused_symbol.clone(),
                )
            };
            let (isomesh_symbols, isomesh_fma, total_fma, extract_symbols, extract_fma) = if on_peer
            {
                let p = p.expect("peer");
                (
                    p.isomesh_symbols,
                    p.isomesh_fma,
                    p.total_fma,
                    p.extract_symbols,
                    p.extract_fma,
                )
            } else {
                (
                    codegen.isomesh_symbols,
                    codegen.isomesh_fma,
                    codegen.total_fma,
                    codegen.extract_symbols,
                    codegen.extract_fma,
                )
            };
            let (mul_add_call_sites, mul_add_mentions) = if on_peer {
                let p = p.expect("peer");
                (p.mul_add_call_sites, p.mul_add_mentions)
            } else {
                (census.call_sites, census.mentions)
            };
            let matches_committed = if on_peer {
                peer_matches_committed
            } else {
                local_matches_committed
            };

            for (index, site) in SITES.iter().enumerate() {
                let local = &readings.0[index];
                let (
                    tried,
                    differing_inputs,
                    ulp,
                    written_bits,
                    fused_bits,
                    variant,
                    reproduced,
                    worst_input,
                ) = if on_peer {
                    let Some(s) = p.expect("peer").sites.get(site.name) else {
                        panic!(
                            "P-117: the peer artefact has no row for `{}`. A missing site \
                             is a column the peer promised and did not deliver.",
                            site.name
                        )
                    };
                    (
                        s.tried,
                        s.differing,
                        s.ulp,
                        s.written_bits.clone(),
                        s.fused_bits.clone(),
                        s.variant.clone(),
                        s.reproduced,
                        s.worst_input.clone(),
                    )
                } else {
                    let reproduced = (site.reproduce)(&local.worst_input).to_bits()
                        == local.worst_written.to_bits();
                    (
                        local.tried,
                        local.differing,
                        local.worst_ulp,
                        format!("{:#018x}", local.worst_written.to_bits()),
                        format!("{:#018x}", local.worst_fused.to_bits()),
                        local.worst_variant.to_string(),
                        reproduced,
                        format_word(&local.worst_input),
                    )
                };

                let probe_fma = if on_peer {
                    p.expect("peer").probes.get(site.probe).copied()
                } else {
                    codegen.probes.get(site.probe).map(|(c, _)| *c)
                };
                let contracted = match probe_fma {
                    Some(count) => (count > 0).to_string(),
                    None => String::from("probe-not-read"),
                };

                // `ulp_difference` is a count of representable doubles and is
                // the registered column, but on a near-cancelling input it runs
                // to ~4.4e18 -- which is honest and unreadable. These two make
                // the row legible without replacing it: `relative_difference`
                // is 1.0 exactly when the written form cancels to zero and the
                // fused one does not, and `worst_case_cancels` says so outright.
                let (unfused_value, fused_value) = (
                    f64::from_bits(
                        u64::from_str_radix(written_bits.trim_start_matches("0x"), 16)
                            .expect("recorded bits"),
                    ),
                    f64::from_bits(
                        u64::from_str_radix(fused_bits.trim_start_matches("0x"), 16)
                            .expect("recorded bits"),
                    ),
                );
                let scale = unfused_value.abs().max(fused_value.abs());
                let relative = if scale > 0.0 {
                    (unfused_value - fused_value).abs() / scale
                } else {
                    0.0
                };
                let cancels = unfused_value == 0.0
                    || fused_value == 0.0
                    || unfused_value.is_sign_negative() != fused_value.is_sign_negative();

                run.record(&[
                    ("expression_site", site.name.replace(' ', "_")),
                    ("file", site.file.to_string()),
                    ("line", site.line.to_string()),
                    ("shape", site.shape.replace(", ", "; ")),
                    ("fused_result", fused_bits),
                    ("unfused_result", written_bits),
                    ("ulp_difference", ulp.to_string()),
                    ("reachable", site.reachable.to_string()),
                    ("target", target.to_string()),
                    (
                        "fmadd_in_known_fused_probe",
                        (known_fused_count > 0).to_string(),
                    ),
                    ("mul_add_call_sites", mul_add_call_sites.to_string()),
                    ("hashes_identical", c2.clone()),
                    ("golden_rows", golden_rows.to_string()),
                    ("c1_holds", c1.to_string()),
                    ("c2_holds", c2.clone()),
                    ("c3_holds", c3.to_string()),
                    // ── extra columns, M-273 ──
                    ("scalar", String::from("f64")),
                    ("on_extraction_path", site.on_path.to_string()),
                    (
                        "product_also_compared",
                        site.product_also_compared.to_string(),
                    ),
                    ("relative_difference", format!("{relative:.17e}")),
                    ("worst_case_cancels", cancels.to_string()),
                    ("unfused_value", format!("{unfused_value:.17e}")),
                    ("fused_value", format!("{fused_value:.17e}")),
                    ("fusion_variant", variant),
                    ("inputs_tried", tried.to_string()),
                    ("inputs_differing", differing_inputs.to_string()),
                    ("worst_input", worst_input),
                    ("probe", site.probe.to_string()),
                    ("probe_reproduced", reproduced.to_string()),
                    ("contracted_in_codegen", contracted),
                    ("known_fused_probe_evidence", known_fused_evidence.clone()),
                    ("known_fused_dispatch_symbol", known_fused_symbol.clone()),
                    (
                        "fma_instructions_in_known_fused_probe",
                        known_fused_count.to_string(),
                    ),
                    ("fma_instructions_binary_total", total_fma.to_string()),
                    ("isomesh_symbols_scanned", isomesh_symbols.to_string()),
                    ("fma_instructions_isomesh_text", isomesh_fma.to_string()),
                    (
                        "isomesh_extract_symbols_scanned",
                        extract_symbols.to_string(),
                    ),
                    ("fma_instructions_isomesh_extract", extract_fma.to_string()),
                    ("mul_add_mentions", mul_add_mentions.to_string()),
                    ("golden_hashes_differing", differing.to_string()),
                    (
                        "golden_matches_committed_fixture",
                        matches_committed.to_string(),
                    ),
                    (
                        "golden_cross_machine_identical",
                        cross_identical.to_string(),
                    ),
                    ("c2_blocker", blocker.clone()),
                    ("control_exact_products_ulp", exact_control.to_string()),
                    (
                        "control_searched_separator_ulp",
                        separator_control.to_string(),
                    ),
                ]);
            }
        }

        println!("\nP-117 verdicts");
        println!(
            "  C1 {} — {} contraction-sensitive expression(s) under crates/isomesh/src/: {}",
            if c1 { "HELD" } else { "FALSIFIED" },
            sensitive.len(),
            sensitive.join(" ")
        );
        println!(
            "  C2 {c2} — {golden_rows} golden rows; {differing} differing across the two \
             targets; local matches the committed fixture: {local_matches_committed}; \
             blocker {blocker}"
        );
        println!(
            "  C3 {} — unattributed divergences: {}",
            if c3 { "HELD" } else { "FALSIFIED" },
            if unattributed.is_empty() {
                String::from("none")
            } else {
                unattributed.join(" ")
            }
        );
        println!(
            "  controls: exact_products {exact_control} ULP (must be 0); \
             searched_separator {separator_control} ULP (must be > 0); known-fused probe \
             {} `{}`",
            codegen.known_fused_count, codegen.known_fused_evidence
        );
        for (site, reading) in SITES.iter().zip(&readings.0) {
            println!(
                "  {:44} {:>4} ULP over {:>9} inputs ({:>9} differing) via {}",
                site.name,
                reading.worst_ulp,
                reading.tried,
                reading.differing,
                reading.worst_variant
            );
        }
        if let Some(p) = &peer {
            println!(
                "  peer: {} on {} with {} — known-fused {} `{}`",
                p.target, p.machine, p.rustc, p.known_fused_count, p.known_fused_evidence
            );
        }
        // `parse_word` exists so the recorded bits are a round trip rather than a
        // one-way format: a reader can take `worst_input` off a row and get the
        // same six floats back. Compared as BITS, not as floats, because bit
        // identity is the property being asserted -- `-0.0 == 0.0` is true and
        // `golden.rs:35-38` is on record that the difference between them is
        // exactly the class of change this whole file is about.
        let round_trip = parse_word(&format_word(&readings.0[S_PLACE].worst_input));
        assert_eq!(
            round_trip.map(f64::to_bits),
            readings.0[S_PLACE].worst_input.map(f64::to_bits),
            "P-117: worst_input does not round trip, so the bits on the row are not the \
             bits that were measured"
        );
    });
}
