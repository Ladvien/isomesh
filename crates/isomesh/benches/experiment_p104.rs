//! **P-104 — the interleaved layout for the active-cell bitmap, which packs samples and not cells.**
//!
//! Ticket: R-104. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p104
//! ```
//!
//! Writes `docs/experiments/p-104.csv`. **Linux only**, `experiment_p12`'s
//! precedent: the registered columns are cycles per sample, instructions per
//! cell and nanoseconds per vertex, and on a governed CPU a nanosecond is not a
//! unit (`M-280`, `M-281`). Every quantity here is taken from
//! `common::counters::Probe`, `ghz` is on every row so a later reader can see
//! what clock the cycles were taken at, and off Linux the bench refuses rather
//! than recording a fabricated zero.
//!
//! # What was missing, and the correction this row carries
//!
//! **The active-cell bitmap packs *samples*, not cells.** `dual.rs:359-381`
//! sets one bit per **sample**, 64 to a `u64`, **along `x` only**, with
//! `bit_row = size[0].div_ceil(64)`. The research doc and
//! `experiment.rs:3192` both say "64 cells per `u64`" and **that is wrong by
//! one**: the cell row is one word *shorter* and is handled separately, by
//! `cell_words = cells_x.div_ceil(64)` (`dual.rs:484`) and `cell_mask`
//! (`:445`). `dual.rs:472-484` documents the asymmetry and prices it —
//! `E-307` measured the packed stage going 0.451–0.470 ns/cell at 128 samples
//! to 0.599–0.637 at 129, about **30% of the stage for a word that could not
//! contain a cell**. So this row is a **sample-plane relayout with a
//! cell/word asymmetry underneath it**, and the asymmetry is why the two arms
//! below cannot simply be "64 cells per word" against "64 cells per block".
//!
//! **The second constraint is load-bearing and is a registered column.**
//! `dual.rs:489-497` walks set bits with `trailing_zeros` and the
//! `active &= active - 1` clear, which visits the active cells of a row in
//! **ascending `x`**, and the comment at `:491-494` says why: *"the same order
//! the scalar loop did, which is what keeps vertex creation order and therefore
//! every index unchanged."* A 4×4×4-interleaved word does **not** enumerate in
//! that order, so the relayout must restore it or C3 fails by construction. The
//! restore is measured — `resort_ns_per_vertex`, by prefix difference against
//! the unordered walk — not footnoted.
//!
//! **The shipped crate can no longer produce the layout C1 is about.**
//! `M-287`'s tax is a property of `DualMesher::values`, whose row stride is
//! `row_stride(size) = size[0] | 1` floats since `A-024` (`dual.rs:333`). At
//! `f32` a 64 KiB plane stride needs `row · size[1] = 16384 = 2¹⁴`, and an odd
//! `row` greater than one can never divide a power of two — so **no caller can
//! reach the taxed layout through the public API any more**. The unpadded
//! layout is therefore reconstructed **bench-local** (`row = n`, plane stride
//! `n² · 4`), which at `n = 128` is exactly the 512-byte row and 64 KiB plane
//! `M-287` named. The padded layout is reconstructed beside it, from the same
//! field, as the control that says the fix still holds.
//!
//! # The two layouts
//!
//! Both arms are mirrors, because `crates/isomesh/src/**` is read-only for
//! Phase 25. Every line of the `along_x` arm has a named origin and the mirror
//! is **asserted to reproduce the shipped extractor's own traversal** in the
//! same run (`R-120`, `R-121`, `M-279`): the ordered active-cell list must
//! equal the eight-corner scalar definition *and* the list the shipped
//! `DualContouring` actually visits, recovered through a recording
//! [`experiment::Recording`] `VertexRule` injected with the public
//! `DualContouring::with_rule`.
//!
//! - **`along_x`** — `build_along_x` is `dual.rs:359-381`, `active_word` is
//!   `dual.rs:405-436` (`any & !all` over four `(y, z)` rows), `cell_mask` is
//!   `dual.rs:445`, and the walk is `dual.rs:489-497`. Bit `k` of word `w` of
//!   row `(y, z)` is `is_inside(sample[64w + k, y, z])`, and `is_inside` is
//!   inlined as `value < 0.0` because `cube.rs:171-173` is `value < R::ZERO` —
//!   **exact zero is outside**, which matters on `box_exact`.
//! - **`interleaved`** — bit `k = bx + 4·by + 16·bz` of the word for sample
//!   block `(BX, BY, BZ)` is `is_inside(sample[4BX+bx, 4BY+by, 4BZ+bz])`. One
//!   `u64` is a 4×4×4 sample block, and one `u64` of *cell* answers is a
//!   4×4×4 cell block — which is the packing the registration names, and the
//!   one place the interleaved layout is *better* structured than the shipped
//!   one, because a 4×4×4 cell block is self-contained in cell index space
//!   while a 64-cell `x` run is not.
//!
//!   The fold is the interleaved counterpart of `any & !all`, built from three
//!   lane-crossing shifts with carry-in from the neighbouring block:
//!
//!   ```text
//!   shift_x(w, wn) = ((w >>  1) & !0x8888…) | ((wn <<  3) & 0x8888…)
//!   shift_y(w, wn) = ((w >>  4) & !0xF000F000F000F000) | ((wn << 12) & 0xF000…)
//!   shift_z(w, wn) = ((w >> 16) & !0xFFFF000000000000) | ((wn << 48) & 0xFFFF…)
//!   ```
//!
//!   A shift by one inside a 4-wide lane crosses the lane, so the mask is not
//!   decoration: without it slot `(3, by, bz)` reads `(0, by+1, bz)` — the
//!   interleaved form of the hole `dual.rs:391-393` names for the along-`x`
//!   form. The eight corner words come from the 2×2×2 neighbourhood of block
//!   words, collapsed `z` then `y` then `x`, and the final fold is the same
//!   four-row `any |= a | b; all &= a & b; any & !all` shape as
//!   `dual.rs:424-436`.
//!
//! # Where `M-287`'s tax actually lives, and why the control needs its own arm
//!
//! `M-287` measured a whole surface-free `SurfaceNets` extraction: 33.10 /
//! 108.51 / 31.39 cycles per sample at 127³ / 128³ / 129³, a **3.37×** tax,
//! with a pad-`z` control at 107.89 proving it is the stride. That extraction
//! is `sample` + `build_inside_bits` + the `place_vertices` word walk +
//! `emit_quads`. Reading the loop nests:
//!
//! - `build_inside_bits` streams `values` **sequentially**;
//! - the `place_vertices` walk reads the *bitmap*, whose row stride at 128 is
//!   `bit_row · 8` = 16 bytes;
//! - `emit_quad_axis` (`dual.rs:697-713`) re-reads `values` with the axis
//!   constant and the **innermost loop over `cells[v]`** — so `AXIS = 0` walks
//!   `z` innermost at a `row · size[1] · 4` = **64 KiB** stride, and
//!   `AXIS = 2` walks `y` innermost at a `row · 4` = **512-byte** stride.
//!
//! So the two aliasing periods `M-287` separated are both in `emit_quads` and
//! **neither is in the bitmap**. The registered VACUITY CONTROL is scored
//! exactly where it is registered — `tax_vs_127` on the *unmodified* bitmap
//! layout, which must reproduce 3.37× to within 10% — and beside it the
//! `m287_*` columns carry an arm that mirrors the whole scaffolding
//! (`fill` + `build_along_x` + walk + the three `emit_quad_axis` sign passes)
//! at both strides, on `f32`, at 127/128/129, so the run can say whether the
//! tax was reproduced *anywhere* rather than only whether the bitmap has it.
//! `sphere_surface_free` — `experiment_p40.rs:84-88`'s field, the canonical
//! sphere sampled from `[10, 10, 10]` so no corner is ever inside — is
//! included as a **ninth field** for that arm's sake, because `M-287`'s own
//! fixture is surface-free and a tax measured on a field with a surface is a
//! different measurement.
//!
//! # Counter windows are siblings, and the stages are prefix differences
//!
//! Zen 3 has six general-purpose counters and `Probe` opens six plus a
//! software event, so two **nested** windows multiplex and
//! `Counts::worst_ratio` refuses. `R-121` paid for that discovery. Every
//! window here is a sibling window over one whole arm, and the stages come out
//! as prefix differences over four nested arms per layout:
//!
//! | arm | contents | difference |
//! |---|---|---|
//! | `build` | build the bitmap from `values` | — |
//! | `fold` | build + the `any & !all` fold, accumulated by XOR | `fold − build` = **the fold** |
//! | `stage` | build + fused fold and set-bit walk, packed ids | `stage − fold` = the walk |
//! | `ordered` | build + walk + restore ascending-`x` order | `ordered − stage` = **the restore** |
//!
//! The fold arm accumulates with `acc ^= word` and **never popcounts**. That
//! is deliberate: this build emits **zero** `popcnt` instructions (there is no
//! `.cargo/config.toml` and no `target-cpu`, so `cfg!(target_feature =
//! "popcnt")` is false and `u64::count_ones` lowers to a ~12-instruction SWAR
//! sequence), and a popcount in the fold arm would put twelve instructions per
//! word into a column that is supposed to measure four. `count_ones` is called
//! nowhere inside a counted window, `target_feature_popcnt` and
//! `count_ones_calls_per_counted_window` are columns saying so, and **no
//! verdict here is contingent on the popcount question**.
//!
//! All arms of a batch share one `inner`, sized on the slowest, so the
//! ~28 `perf_event` system calls a window costs are identical across the arms
//! and cancel in a difference. Each repetition takes one window per arm back
//! to back, interleaved, so a drifting clock moves every arm together
//! (`R-105`, `R-114`).
//!
//! # SHARE
//!
//! Each clause's reachable share is a column, measured in this run rather than
//! quoted:
//!
//! - **C1's share is the tax itself.** `M-287`'s 3.37× is a whole-extraction
//!   number, so what C1 moves is bounded by
//!   `bitmap_share_of_m287_scaffold` — the `along_x` bitmap stage's cycles per
//!   sample over the `m287` scaffold's, on the same field at the same size, in
//!   the same run — and by `c1_ceiling_at_128`, `1/(1 − share·(1 − 1/tax))`
//!   computed from this row's own numbers. C1 is scored on `tax_vs_127` at
//!   128³, as registered, and **only if the VACUITY CONTROL reproduced**;
//!   otherwise `c1_holds` is `vacuous`, because a harness that cannot see the
//!   tax cannot measure its removal.
//! - **C2's share is `instructions_per_cell`**, and it is the whole clause: the
//!   interleaved fold's instructions per cell against the `along_x` fold's,
//!   inside ±20%. `instruction_ratio_fold` is the number and
//!   `c2_holds` reads it. Instruction counts are **deterministic**; the cycle
//!   form is beside it as `cycle_ratio_fold` and does not carry the verdict.
//! - **C3 has no share; it is an equality.** Four ordered active-cell lists
//!   must be the same list — the eight-corner scalar definition, the `along_x`
//!   walk, the restored interleaved walk, and the list the *shipped* extractor
//!   traverses. The first, second and fourth are **asserted**: they are claims
//!   about the mirror, and a mirror that does not reproduce the shipped
//!   traversal licenses nothing (`M-279`). The third is C3 and is recorded.
//!
//! # Which form carries the verdict
//!
//! `cycles_per_sample` and therefore `tax_vs_127` are **cycle** quantities,
//! because `M-287` is a cycle quantity and C1 is denominated in it. Cycle
//! counts on this machine are not reproducible — `R-105` watched an identical
//! binary's cycle ratio band move from 0.984 to 1.035 across three runs while
//! its instruction counts held to four figures — so `instructions_per_sample`,
//! `tax_vs_127_instructions` and `tax_vs_127_l1d_misses` are beside every
//! cycle column. **For C1 the cycle form is the only form that can carry the
//! verdict**, because a cache-set aliasing tax is by construction invisible to
//! an instruction count: the taxed and untaxed layouts execute the *same
//! instructions*. `tax_vs_127_instructions` is reported precisely so a reader
//! can see it sitting at 1.0 while the cycle form moves — that is the
//! signature of a stride effect, and its absence would mean the arm is
//! measuring work rather than layout. C2's verdict, by contrast, reads the
//! instruction form, which is the reproducible one.

mod common;

/// `M-287`'s own fixture — 127³, 128³, 129³ — plus 64³.
///
/// 128 is the spike: at `f32` with an unpadded `row = 128` the row stride is
/// 512 bytes and the plane stride exactly 64 KiB = 2¹⁶. 127 is the tax
/// baseline and 129 the other neighbour; 64³ is the small-grid control, where
/// the whole `values` array is 1 MB and no stride can alias out of L2.
const RESOLUTIONS: [u32; 4] = [64, 127, 128, 129];

/// The resolution every `tax_vs_127` is divided by.
const TAX_BASELINE: u32 = 127;

/// The resolution C1 is scored at.
const TAX_SPIKE: u32 = 128;

/// C1's bar: the transposed packing must bring the 128³ tax below this.
const C1_BAR: f64 = 1.5;

/// C2's bar: the fold's instructions per cell must survive the relayout
/// within ±20%.
const C2_TOLERANCE: f64 = 0.20;

/// `M-287`'s measured tax, which the VACUITY CONTROL must reproduce.
const M287_TAX: f64 = 3.37;

/// How close the control has to get. Registered as "to within 10%".
const M287_TOLERANCE: f64 = 0.10;

/// Windows per quantity per arm. Odd, so a median is a reading rather than a
/// mean of two.
const REPS: usize = 5;

/// The field name carrying `M-287`'s own surface-free fixture.
const CONTROL_FIELD: &str = "sphere_surface_free";

#[cfg(target_os = "linux")]
mod experiment {
    use std::cell::RefCell;
    use std::hint::black_box;
    use std::rc::Rc;
    use std::time::Instant;

    use isomesh::dual::{CellVertices, VertexRule};
    use isomesh::dual_contouring::DualContouring;
    use isomesh::fields::{ReferenceField, Sphere};
    use isomesh::surface_nets::Centroid;
    use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};
    use crate::{
        C1_BAR, C2_TOLERANCE, CONTROL_FIELD, M287_TAX, M287_TOLERANCE, REPS, RESOLUTIONS,
        TAX_BASELINE, TAX_SPIKE,
    };

    /// About this long per counter window, so the `perf_event` round trip is
    /// negligible against what it measures.
    const TARGET_BATCH_NS: f64 = 10_000_000.0;

    /// Ceiling on the batch factor, so a grid that got faster than expected
    /// cannot turn one window into a minute.
    const MAX_INNER: usize = 1 << 20;

    /// Bits whose `bx` is 3, the lane that must take its `+x` neighbour's bit 0.
    const MASK_BX3: u64 = 0x8888_8888_8888_8888;
    /// Bits whose `by` is 3.
    const MASK_BY3: u64 = 0xF000_F000_F000_F000;
    /// Bits whose `bz` is 3.
    const MASK_BZ3: u64 = 0xFFFF_0000_0000_0000;

    // ------------------------------------------------------------ the values

    /// `sample_grid` (`sdf.rs:167-194`), mirrored with a caller-chosen row
    /// stride.
    ///
    /// `row` is a parameter rather than `size[0] | 1` because the whole point
    /// of this row is the layout `A-024` removed: `row = n` at `n = 128` gives
    /// the 512-byte row and 64 KiB plane `M-287` measured, and the shipped
    /// `row_stride` can no longer produce it. The excess slots are left at
    /// zero exactly as `sample_grid` leaves them.
    fn fill<F: Sdf<Scalar = f32>>(
        values: &mut Vec<f32>,
        field: &F,
        n: u32,
        origin: [f32; 3],
        cell: f32,
        row: usize,
    ) {
        let nx = n as usize;
        let pad = row - nx;
        values.clear();
        values.reserve(row * nx * nx);
        for z in 0..n {
            let fz = origin[2] + cell * (z as f32);
            for y in 0..n {
                let fy = origin[1] + cell * (y as f32);
                for x in 0..n {
                    values.push(field.sample([origin[0] + cell * (x as f32), fy, fz]));
                }
                for _ in 0..pad {
                    values.push(0.0);
                }
            }
        }
    }

    /// One `(field, resolution)`'s sample array, at the unpadded stride.
    ///
    /// Immutable while the arms run, so every arm can borrow it at once and
    /// only its own bitmap scratch is private to it.
    struct Grid {
        values: Vec<f32>,
        /// Samples per row. `n`, deliberately unpadded.
        row: usize,
        /// Samples per axis.
        n: usize,
        /// Cells per axis: `n` samples span `n − 1` cells.
        c: usize,
    }

    // ------------------------------------------------------- layout: along x

    /// `dual.rs:359-381`, verbatim apart from the destination. Returns
    /// `bit_row`.
    ///
    /// One bit per **sample**, 64 to a `u64`, along `x` only. `is_inside` is
    /// inlined as `value < 0.0` because `cube.rs:171-173` is `value < R::ZERO`
    /// — exact zero is **outside**.
    fn build_along_x(dst: &mut Vec<u64>, values: &[f32], row: usize, n: usize) -> usize {
        let rows = n * n;
        let bit_row = n.div_ceil(64);
        dst.clear();
        dst.resize(bit_row * rows, 0);
        for r in 0..rows {
            let src = row * r;
            let base_dst = bit_row * r;
            for w in 0..bit_row {
                let base = w * 64;
                let count = (n - base).min(64);
                let mut word = 0u64;
                for k in 0..count {
                    // Branchless on purpose: the whole point of the bitmap is
                    // that the sign test stops being a branch per corner.
                    word |= u64::from(values[src + base + k] < 0.0) << k;
                }
                dst[base_dst + w] = word;
            }
        }
        bit_row
    }

    /// `dual.rs:385`.
    #[inline]
    fn inside_word(inside: &[u64], bit_row: usize, n: usize, w: usize, y: usize, z: usize) -> u64 {
        inside[bit_row * (y + n * z) + w]
    }

    /// `dual.rs:395`. The high bit comes from the next word, or the cell
    /// straddling a word boundary reads its `+x` corner as outside.
    #[inline]
    fn inside_word_shifted(
        inside: &[u64],
        bit_row: usize,
        n: usize,
        w: usize,
        y: usize,
        z: usize,
    ) -> u64 {
        let lo = inside_word(inside, bit_row, n, w, y, z);
        let hi = if w + 1 < bit_row {
            inside_word(inside, bit_row, n, w + 1, y, z)
        } else {
            0
        };
        (lo >> 1) | (hi << 63)
    }

    /// `dual.rs:424-436`. Sixty-four active-cell answers in four fused word
    /// operations per row.
    #[inline]
    fn active_word(inside: &[u64], bit_row: usize, n: usize, w: usize, y: usize, z: usize) -> u64 {
        let mut any = 0u64;
        let mut all = !0u64;
        for dz in 0..2usize {
            for dy in 0..2usize {
                let a = inside_word(inside, bit_row, n, w, y + dy, z + dz);
                let b = inside_word_shifted(inside, bit_row, n, w, y + dy, z + dz);
                any |= a | b;
                all &= a & b;
            }
        }
        any & !all
    }

    /// `dual.rs:445`. `1u64 << 64` is undefined, so the full-word case is
    /// named rather than computed.
    #[inline]
    fn cell_mask(w: usize, cells_x: usize) -> u64 {
        let base = w * 64;
        let remaining = cells_x.saturating_sub(base);
        if remaining >= 64 {
            !0
        } else {
            (1u64 << remaining) - 1
        }
    }

    /// The fold alone, accumulated by XOR so no `count_ones` enters the window.
    fn fold_along_x(inside: &[u64], bit_row: usize, n: usize, c: usize) -> u64 {
        let cell_words = c.div_ceil(64);
        let mut acc = 0u64;
        for z in 0..c {
            for y in 0..c {
                for w in 0..cell_words {
                    acc ^= active_word(inside, bit_row, n, w, y, z) & cell_mask(w, c);
                }
            }
        }
        acc
    }

    /// `dual.rs:487-497`. Emits the **packed** id `word_index · 64 + bit`,
    /// which is monotone in grid order.
    fn walk_along_x(inside: &[u64], bit_row: usize, n: usize, c: usize, out: &mut Vec<u32>) {
        let cell_words = c.div_ceil(64);
        out.clear();
        for z in 0..c {
            for y in 0..c {
                let row_base = ((y + c * z) * cell_words) as u32;
                for w in 0..cell_words {
                    let mut active = active_word(inside, bit_row, n, w, y, z) & cell_mask(w, c);
                    let base = (row_base + w as u32) << 6;
                    // `active &= active - 1` clears the lowest set bit, so this
                    // visits the active cells of the row in ascending `x`.
                    while active != 0 {
                        out.push(base + active.trailing_zeros());
                        active &= active - 1;
                    }
                }
            }
        }
    }

    /// The same walk emitting **grid** ids `x + c·(y + c·z)` directly.
    ///
    /// The along-`x` layout's order restore is the identity, so this arm does
    /// the same work as [`walk_along_x`] with the row base hoisted out of the
    /// inner loop — one add per set bit, exactly `dual.rs:496`. Its prefix
    /// difference against [`walk_along_x`] is this layout's
    /// `resort_ns_per_vertex`, and it is a genuine measurement rather than a
    /// zero that could not have been non-zero.
    fn walk_along_x_ordered(
        inside: &[u64],
        bit_row: usize,
        n: usize,
        c: usize,
        out: &mut Vec<u32>,
    ) {
        let cell_words = c.div_ceil(64);
        out.clear();
        for z in 0..c {
            for y in 0..c {
                let grid_base = (y + c * z) * c;
                for w in 0..cell_words {
                    let mut active = active_word(inside, bit_row, n, w, y, z) & cell_mask(w, c);
                    let base = (grid_base + w * 64) as u32;
                    while active != 0 {
                        out.push(base + active.trailing_zeros());
                        active &= active - 1;
                    }
                }
            }
        }
    }

    // --------------------------------------------------- layout: interleaved

    /// One `u64` per 4×4×4 **sample** block, bit `k = bx + 4·by + 16·bz`.
    ///
    /// Returns the blocks per axis, `n.div_ceil(4)`. Samples past the grid
    /// leave their bit at zero — outside — which is the same convention the
    /// along-`x` tail uses and the same one `dual.rs:397-401` uses for the
    /// missing `+x` word.
    fn build_interleaved(dst: &mut Vec<u64>, values: &[f32], row: usize, n: usize) -> usize {
        let nb = n.div_ceil(4);
        dst.clear();
        dst.resize(nb * nb * nb, 0);
        for bz in 0..nb {
            let zs = (n - 4 * bz).min(4);
            for by in 0..nb {
                let ys = (n - 4 * by).min(4);
                for bx in 0..nb {
                    let xs = (n - 4 * bx).min(4);
                    let mut word = 0u64;
                    for k in 0..zs {
                        for j in 0..ys {
                            let src = row * ((4 * by + j) + n * (4 * bz + k)) + 4 * bx;
                            let shift = 4 * j + 16 * k;
                            for i in 0..xs {
                                word |= u64::from(values[src + i] < 0.0) << (shift + i);
                            }
                        }
                    }
                    dst[bx + nb * (by + nb * bz)] = word;
                }
            }
        }
        nb
    }

    /// The block word, or zero outside the block grid — `dual.rs:397-401`'s
    /// convention in three axes instead of one.
    #[inline]
    fn block_word(blocks: &[u64], nb: usize, bx: usize, by: usize, bz: usize) -> u64 {
        if bx >= nb || by >= nb || bz >= nb {
            0
        } else {
            blocks[bx + nb * (by + nb * bz)]
        }
    }

    /// Slot `(bx, by, bz)` ← the sample at `bx + 1`.
    ///
    /// The mask is not decoration. `w >> 1` at slot `(3, by, bz)` holds bit
    /// `4(by+1) + 16bz`, which is the *next row's* first sample, so the
    /// `bx = 3` lane must instead take bit `4by + 16bz` of the `+x`
    /// neighbouring block — `wn << 3`. Dropping the mask is the interleaved
    /// form of the hole `dual.rs:391-393` names.
    #[inline]
    fn shift_x(w: u64, wn: u64) -> u64 {
        ((w >> 1) & !MASK_BX3) | ((wn << 3) & MASK_BX3)
    }

    /// Slot `(bx, by, bz)` ← the sample at `by + 1`.
    #[inline]
    fn shift_y(w: u64, wn: u64) -> u64 {
        ((w >> 4) & !MASK_BY3) | ((wn << 12) & MASK_BY3)
    }

    /// Slot `(bx, by, bz)` ← the sample at `bz + 1`.
    #[inline]
    fn shift_z(w: u64, wn: u64) -> u64 {
        ((w >> 16) & !MASK_BZ3) | ((wn << 48) & MASK_BZ3)
    }

    /// Sixty-four active-cell answers for the 4×4×4 cell block at
    /// `(bx, by, bz)`.
    ///
    /// `dual.rs:424-436`'s `any & !all` over the eight corners, with the eight
    /// corner words derived from the 2×2×2 neighbourhood of block words by
    /// collapsing `z`, then `y`, then `x`. Eight loads, the same as
    /// `active_word`'s four rows of two.
    #[inline]
    fn active_block(blocks: &[u64], nb: usize, bx: usize, by: usize, bz: usize) -> u64 {
        // The `z` collapse first. `lo[dy][dx]` is the block word at offset
        // `(dx, dy)`; `hi[dy][dx]` is that word with its `z` slots advanced by
        // one. Two arrays rather than one `[ez][dy][dx]`, so no index is a loop
        // variable and the eight loads stay eight loads.
        let mut lo = [[0u64; 2]; 2];
        let mut hi = [[0u64; 2]; 2];
        for (dy, (lo_row, hi_row)) in lo.iter_mut().zip(hi.iter_mut()).enumerate() {
            for (dx, (lo_slot, hi_slot)) in lo_row.iter_mut().zip(hi_row.iter_mut()).enumerate() {
                let w = block_word(blocks, nb, bx + dx, by + dy, bz);
                *lo_slot = w;
                *hi_slot = shift_z(w, block_word(blocks, nb, bx + dx, by + dy, bz + 1));
            }
        }
        // Then `y`, then `x`, folded as it goes — the same four-pair shape as
        // `active_word`'s four `(y, z)` rows. For each `ez` plane, `ey = 0` is
        // the `dy = 0` row itself and `ey = 1` folds it with the `dy = 1` row.
        let mut any = 0u64;
        let mut all = !0u64;
        for plane in [&lo, &hi] {
            let (row0, row1) = (plane[0], plane[1]);
            let folded = [shift_y(row0[0], row1[0]), shift_y(row0[1], row1[1])];
            for row in [row0, folded] {
                let a = row[0];
                let b = shift_x(a, row[1]);
                any |= a | b;
                all &= a & b;
            }
        }
        any & !all
    }

    /// `(1 << span) - 1` in every nibble: which `bx` lanes of a block are real
    /// cells.
    #[inline]
    fn nibble_mask(span: u32) -> u64 {
        ((1u64 << span) - 1) * 0x1111_1111_1111_1111
    }

    /// `cell_mask`'s counterpart: the `(by, bz)` half of a block's cell mask,
    /// hoisted out of the `bx` loop the way a real implementation would hoist
    /// it.
    #[inline]
    fn block_row_template(by: usize, bz: usize, c: usize) -> u64 {
        let ys = (c - 4 * by).min(4);
        let zs = (c - 4 * bz).min(4);
        let mut m = 0u64;
        for k in 0..zs {
            for j in 0..ys {
                m |= 0xFu64 << (4 * j + 16 * k);
            }
        }
        m
    }

    /// The `bx` half, branching on the full-block case exactly as `cell_mask`
    /// does.
    #[inline]
    fn block_cell_mask(template: u64, bx: usize, c: usize) -> u64 {
        let base = 4 * bx;
        if base + 4 <= c {
            template
        } else {
            template & nibble_mask((c - base) as u32)
        }
    }

    /// The interleaved fold alone, accumulated by XOR.
    fn fold_interleaved(blocks: &[u64], nb: usize, cb: usize, c: usize) -> u64 {
        let mut acc = 0u64;
        for bz in 0..cb {
            for by in 0..cb {
                let template = block_row_template(by, bz, c);
                for bx in 0..cb {
                    acc ^= active_block(blocks, nb, bx, by, bz) & block_cell_mask(template, bx, c);
                }
            }
        }
        acc
    }

    /// The interleaved walk. Emits the **packed** id `block_index · 64 + k`,
    /// which is *not* monotone in grid order — that is the whole point, and
    /// [`restore_interleaved`] is what it costs to fix.
    fn walk_interleaved(blocks: &[u64], nb: usize, cb: usize, c: usize, out: &mut Vec<u32>) {
        out.clear();
        for bz in 0..cb {
            for by in 0..cb {
                let template = block_row_template(by, bz, c);
                for bx in 0..cb {
                    let mut active =
                        active_block(blocks, nb, bx, by, bz) & block_cell_mask(template, bx, c);
                    let base = ((bx + cb * (by + cb * bz)) as u32) << 6;
                    while active != 0 {
                        out.push(base + active.trailing_zeros());
                        active &= active - 1;
                    }
                }
            }
        }
    }

    /// Per-`(bz, y)` row buffers for the order restore, retained across
    /// blocks so the restore does not pay an allocation per slab.
    struct Buckets {
        rows: Vec<Vec<u32>>,
    }

    impl Buckets {
        fn new() -> Self {
            Self { rows: Vec::new() }
        }

        /// `4 · c` buffers: one per `(bz-local z, global y)` inside a
        /// four-deep `z` slab.
        fn fit(&mut self, c: usize) {
            self.rows.resize_with(4 * c, Vec::new);
        }
    }

    /// The interleaved walk **plus** the restore to ascending-`x`,
    /// row-major order, emitting grid ids.
    ///
    /// The restore has to span a whole four-deep `z` slab, not one block:
    /// inside a slab the blocks are visited `(by, bx)`, so block `(bz, by)`
    /// emits `z = 4bz+1` before block `(bz, by+1)` emits `z = 4bz`, and a
    /// per-block regroup would leave the list out of order. With `4·c` row
    /// buffers per slab the walk is one pass, the concatenation is `z`-major
    /// then `y`, and within one buffer the ids arrive in ascending `bx` and
    /// then ascending `bx`-local `x` — so ascending `x`. That is the order
    /// `dual.rs:491-494` requires, restored in `O(active cells)` with no sort.
    fn restore_interleaved(
        blocks: &[u64],
        nb: usize,
        cb: usize,
        c: usize,
        buckets: &mut Buckets,
        out: &mut Vec<u32>,
    ) {
        out.clear();
        for bz in 0..cb {
            let z0 = 4 * bz;
            let zs = (c - z0).min(4);
            for row in buckets.rows.iter_mut() {
                row.clear();
            }
            for by in 0..cb {
                let template = block_row_template(by, bz, c);
                for bx in 0..cb {
                    let mut active =
                        active_block(blocks, nb, bx, by, bz) & block_cell_mask(template, bx, c);
                    let x0 = 4 * bx;
                    let y0 = 4 * by;
                    while active != 0 {
                        let k = active.trailing_zeros() as usize;
                        active &= active - 1;
                        let y = y0 + ((k >> 2) & 3);
                        let lz = k >> 4;
                        let x = x0 + (k & 3);
                        buckets.rows[lz * c + y].push((x + c * (y + c * (z0 + lz))) as u32);
                    }
                }
            }
            for lz in 0..zs {
                for y in 0..c {
                    out.extend_from_slice(&buckets.rows[lz * c + y]);
                }
            }
        }
    }

    // ----------------------------------------------------- the M-287 scaffold

    /// `emit_quad_axis` (`dual.rs:682-713`)'s sign re-read, with the axis a
    /// constant exactly as the shipped code has it (`M-285`).
    ///
    /// This is the loop nest that carries `M-287`'s tax. `AXIS = 0` has `z`
    /// innermost, so its stride is `row · size[1] · 4` — 64 KiB at an
    /// unpadded 128 — and `AXIS = 2` has `y` innermost, so its stride is
    /// `row · 4` = 512 bytes. The crossing count is returned so the pass
    /// cannot be elided.
    fn emit_scan<const AXIS: usize>(values: &[f32], row: usize, n: usize, c: usize) -> u64 {
        let axis = AXIS;
        let u = (AXIS + 1) % 3;
        let v = (AXIS + 2) % 3;
        let mut crossings = 0u64;
        let mut p = [0usize; 3];
        for a in 0..n - 1 {
            for b in 1..c {
                for d in 1..c {
                    p[axis] = a;
                    p[u] = b;
                    p[v] = d;
                    let s0 = p[0] + row * (p[1] + n * p[2]);
                    let mut q = p;
                    q[axis] += 1;
                    let s1 = q[0] + row * (q[1] + n * q[2]);
                    crossings += u64::from((values[s0] < 0.0) != (values[s1] < 0.0));
                }
            }
        }
        crossings
    }

    /// Everything `DualMesher::extract` does on a field with no surface, at a
    /// caller-chosen row stride.
    ///
    /// `sample` + `build_inside_bits` + the `place_vertices` word walk + the
    /// three `emit_quad_axis` sign passes. `smooth` and `emit_vertices` are
    /// omitted because with no vertices they do nothing, and the quad bodies
    /// are omitted because with no crossing they never run — so on
    /// `sphere_surface_free` this **is** the extraction `M-287` measured, and
    /// on a field with a surface it is the same scaffolding plus the crossings
    /// it would have emitted.
    fn m287_scaffold<F: Sdf<Scalar = f32>>(
        scratch: &mut M287Scratch,
        field: &F,
        n: u32,
        origin: [f32; 3],
        cell: f32,
        row: usize,
    ) -> u64 {
        let nn = n as usize;
        let c = nn - 1;
        fill(&mut scratch.values, field, n, origin, cell, row);
        let bit_row = build_along_x(&mut scratch.inside, &scratch.values, row, nn);
        walk_along_x(&scratch.inside, bit_row, nn, c, &mut scratch.out);
        let mut acc = scratch.out.len() as u64;
        acc ^= emit_scan::<0>(&scratch.values, row, nn, c);
        acc ^= emit_scan::<1>(&scratch.values, row, nn, c);
        acc ^= emit_scan::<2>(&scratch.values, row, nn, c);
        acc
    }

    /// Buffers the `M-287` scaffold owns, one set per stride so the two arms
    /// can be measured back to back.
    struct M287Scratch {
        values: Vec<f32>,
        inside: Vec<u64>,
        out: Vec<u32>,
    }

    impl M287Scratch {
        fn new() -> Self {
            Self {
                values: Vec::new(),
                inside: Vec::new(),
                out: Vec::new(),
            }
        }
    }

    // ------------------------------------------------------------- the clock

    /// Cycles, instructions, L1D read misses and nanoseconds from one window.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
        l1d_misses: f64,
        nanos: f64,
    }

    /// One counter window over `inner` repetitions of `body`, divided by
    /// `inner`.
    ///
    /// Sibling, never nested: Zen 3 has six general-purpose counters and
    /// `Probe` opens six, so a window inside a window multiplexes and
    /// `worst_ratio` refuses (`R-121`).
    fn window(probe: &mut Probe, inner: usize, body: &mut dyn FnMut()) -> Counted {
        probe.reset_and_enable();
        let started = Instant::now();
        for _ in 0..inner {
            body();
        }
        let nanos = started.elapsed().as_nanos() as f64;
        probe.disable();
        let counts = probe.read();
        assert!(
            counts.worst_ratio() >= MIN_TIME_RATIO,
            "a counter ran only {:.1}% of the time it was enabled, so its value is an \
             extrapolation rather than a measurement",
            counts.worst_ratio() * 100.0
        );
        let scale = 1.0 / inner as f64;
        Counted {
            cycles: counts.cycles.count as f64 * scale,
            instructions: counts.instructions.count as f64 * scale,
            l1d_misses: counts.l1d_read_misses.count as f64 * scale,
            nanos: nanos * scale,
        }
    }

    fn median(mut values: Vec<f64>) -> f64 {
        values.sort_by(f64::total_cmp);
        values[values.len() / 2]
    }

    /// Repetitions to make one window about [`TARGET_BATCH_NS`] long.
    fn pick_inner(pass_ns: f64) -> usize {
        ((TARGET_BATCH_NS / pass_ns.max(1.0)).ceil() as usize).clamp(1, MAX_INNER)
    }

    /// Time one pass, to size the batch.
    fn probe_pass_ns(body: &mut dyn FnMut()) -> f64 {
        body();
        let started = Instant::now();
        body();
        started.elapsed().as_nanos() as f64
    }

    /// Every arm of a batch, measured **interleaved**: one window per arm per
    /// repetition, back to back.
    ///
    /// `R-105` and `R-114`'s lesson. Five windows of arm one and then five of
    /// arm two gives arm one the clock the machine happened to have at the
    /// start of the row, and on a governed CPU that is not the clock arm two
    /// gets. Interleaving makes a drifting clock move every arm together, so
    /// it largely cancels in a ratio or a difference.
    ///
    /// All arms share one `inner`, sized on the **slowest**, so the batching
    /// overhead is identical across the arms and cancels in a prefix
    /// difference.
    fn measure_arms(probe: &mut Probe, arms: &mut [&mut dyn FnMut()]) -> Vec<Counted> {
        let slowest = arms
            .iter_mut()
            .map(|body| probe_pass_ns(&mut **body))
            .fold(0.0f64, f64::max);
        let inner = pick_inner(slowest);
        let mut counted: Vec<Vec<Counted>> =
            (0..arms.len()).map(|_| Vec::with_capacity(REPS)).collect();
        for _ in 0..REPS {
            for (slot, body) in counted.iter_mut().zip(arms.iter_mut()) {
                slot.push(window(probe, inner, &mut **body));
            }
        }
        counted
            .into_iter()
            .map(|runs| Counted {
                cycles: median(runs.iter().map(|c| c.cycles).collect()),
                instructions: median(runs.iter().map(|c| c.instructions).collect()),
                l1d_misses: median(runs.iter().map(|c| c.l1d_misses).collect()),
                nanos: median(runs.iter().map(|c| c.nanos).collect()),
            })
            .collect()
    }

    /// A difference of two prefix windows, floored at zero.
    ///
    /// A negative stage cost is noise, not a negative cost; it is floored and
    /// the raw prefixes are recorded beside it so the flooring is visible.
    #[inline]
    fn stage(outer: Counted, inner: Counted) -> Counted {
        Counted {
            cycles: (outer.cycles - inner.cycles).max(0.0),
            instructions: (outer.instructions - inner.instructions).max(0.0),
            l1d_misses: (outer.l1d_misses - inner.l1d_misses).max(0.0),
            nanos: (outer.nanos - inner.nanos).max(0.0),
        }
    }

    // ------------------------------------------------------- the shipped list

    /// The shipped extractor's own ordered active-cell list.
    ///
    /// `place_vertices` calls `rule.place` exactly once per active cell in
    /// traversal order, **before** it knows whether the rule will produce
    /// anything (`dual.rs:516`), so wrapping the rule recovers the traversal
    /// itself rather than a filtered view of it. `Centroid` is the inner rule
    /// because it is public and because the recorded list is rule-independent.
    pub(crate) struct Recording<V> {
        inner: V,
        log: Rc<RefCell<Vec<u32>>>,
        cells: u32,
    }

    impl<R: Real, V: VertexRule<R>> VertexRule<R> for Recording<V> {
        fn place<S: Sdf<Scalar = R>>(
            &self,
            sdf: &S,
            corner: &[R; 8],
            base: [u32; 3],
            origin: [R; 3],
            cell_size: R,
            out: &mut CellVertices<R>,
        ) {
            self.log
                .borrow_mut()
                .push(base[0] + self.cells * (base[1] + self.cells * base[2]));
            self.inner.place(sdf, corner, base, origin, cell_size, out);
        }
    }

    /// Are two meshes bit-identical in every position, normal and index?
    ///
    /// Transcribed from `experiment_p101.rs:1260` at `f32`.
    fn bit_identical(a: &MeshBuffer<f32>, b: &MeshBuffer<f32>) -> bool {
        a.indices == b.indices
            && a.positions.len() == b.positions.len()
            && a.normals.len() == b.normals.len()
            && a.positions
                .iter()
                .zip(&b.positions)
                .all(|(p, q)| (0..3).all(|k| p[k].to_bits() == q[k].to_bits()))
            && a.normals
                .iter()
                .zip(&b.normals)
                .all(|(p, q)| (0..3).all(|k| p[k].to_bits() == q[k].to_bits()))
    }

    /// The shipped traversal's ordered list, plus whether the recording
    /// wrapper left the mesh alone.
    fn shipped_list<F: Sdf<Scalar = f32>>(
        field: &F,
        n: u32,
        origin: [f32; 3],
        cell: f32,
    ) -> (Vec<u32>, bool) {
        let shape = RuntimeShape3::new([n; 3]).expect("the fixture fits u32");
        let mut plain = MeshBuffer::<f32>::new();
        DualContouring::<f32, _>::with_rule(Centroid)
            .extract(field, &shape, origin, cell, &mut plain)
            .expect("extraction");

        let log = Rc::new(RefCell::new(Vec::new()));
        let rule = Recording {
            inner: Centroid,
            log: Rc::clone(&log),
            cells: n - 1,
        };
        let mut wrapped = MeshBuffer::<f32>::new();
        DualContouring::<f32, _>::with_rule(rule)
            .extract(field, &shape, origin, cell, &mut wrapped)
            .expect("extraction");

        // `validate::mesh_hash` takes `&MeshBuffer<f64>`
        // (`src/validate/mesh_hash.rs:96`) and this fixture is `f32`, so
        // transparency is checked bit for bit here instead — a strictly
        // stronger reading than hash equality, and it needs no second grid.
        let transparent = bit_identical(&plain, &wrapped);
        (
            Rc::try_unwrap(log).expect("sole owner").into_inner(),
            transparent,
        )
    }

    /// The definition: eight loads and eight comparisons per cell, in grid
    /// order.
    fn active_scalar(grid: &Grid) -> Vec<u32> {
        let mut out = Vec::new();
        for z in 0..grid.c {
            for y in 0..grid.c {
                for x in 0..grid.c {
                    let mut inside = 0u32;
                    for corner in 0..8usize {
                        let s = (x + (corner & 1))
                            + grid.row
                                * ((y + ((corner >> 1) & 1)) + grid.n * (z + ((corner >> 2) & 1)));
                        if grid.values[s] < 0.0 {
                            inside += 1;
                        }
                    }
                    if inside != 0 && inside != 8 {
                        out.push((x + grid.c * (y + grid.c * z)) as u32);
                    }
                }
            }
        }
        out
    }

    /// Packed along-`x` ids → grid ids. Monotone, so it does not change the
    /// order; run outside every counted window.
    fn along_x_packed_to_grid(packed: &[u32], cell_words: usize, c: usize) -> Vec<u32> {
        packed
            .iter()
            .map(|p| {
                let i = (*p >> 6) as usize;
                let bit = (*p & 63) as usize;
                let rest = i / cell_words;
                let x = (i % cell_words) * 64 + bit;
                (x + c * rest) as u32
            })
            .collect()
    }

    // ------------------------------------------------------------ one fixture

    /// One layout's four prefix windows and the stages they imply.
    #[derive(Clone, Copy)]
    struct Layout {
        build: Counted,
        fold: Counted,
        stage_total: Counted,
        ordered: Counted,
    }

    impl Layout {
        /// The `any & !all` fold alone: `fold − build`.
        fn fold_only(&self) -> Counted {
            stage(self.fold, self.build)
        }

        /// The order restore alone: `ordered − stage_total`.
        fn restore(&self) -> Counted {
            stage(self.ordered, self.stage_total)
        }
    }

    /// Everything one `(field, resolution)` produced.
    struct Measured {
        field: &'static str,
        resolution: u32,
        samples: usize,
        cells: usize,
        active_cells: usize,
        along_x: Layout,
        interleaved: Layout,
        along_x_words: usize,
        interleaved_words: usize,
        m287_unpadded: Counted,
        m287_padded: Counted,
        m287_row: usize,
        mesh_identical: bool,
        matches_shipped: bool,
        matches_scalar: bool,
        wrapper_transparent: bool,
    }

    #[allow(
        clippy::too_many_lines,
        reason = "one fixture is one function: ten arms, four list checks and their assertions \
                  belong together, and splitting them would hide which borrow feeds which window"
    )]
    fn measure<F: Sdf<Scalar = f32>>(
        name: &'static str,
        field: &F,
        n: u32,
        origin: [f32; 3],
        cell: f32,
        probe: &mut Probe,
    ) -> Measured {
        let nn = n as usize;
        let c = nn - 1;
        let cb = c.div_ceil(4);
        let cell_words = c.div_ceil(64);

        let mut grid = Grid {
            values: Vec::new(),
            row: nn,
            n: nn,
            c,
        };
        fill(&mut grid.values, field, n, origin, cell, nn);

        // ---- the lists, outside every counted window.
        let mut inside = Vec::new();
        let bit_row = build_along_x(&mut inside, &grid.values, grid.row, nn);
        let mut blocks = Vec::new();
        let nb = build_interleaved(&mut blocks, &grid.values, grid.row, nn);
        assert_eq!(
            bit_row,
            nn.div_ceil(64),
            "bit_row is div_ceil(64) of samples"
        );
        assert!(cb <= nb, "a cell block must have a sample block over it");

        let mut packed = Vec::new();
        walk_along_x(&inside, bit_row, nn, c, &mut packed);
        let along_x_ids = along_x_packed_to_grid(&packed, cell_words, c);
        let mut ordered_x = Vec::new();
        walk_along_x_ordered(&inside, bit_row, nn, c, &mut ordered_x);
        let mut buckets = Buckets::new();
        buckets.fit(c);
        let mut restored = Vec::new();
        restore_interleaved(&blocks, nb, cb, c, &mut buckets, &mut restored);

        let scalar_ids = active_scalar(&grid);
        let (shipped_ids, wrapper_transparent) = shipped_list(field, n, origin, cell);

        // The mirror has to reproduce the shipped structure in the same run or
        // it licenses nothing (`M-279`, `R-120`, `R-121`). These three are
        // claims about the harness, not clauses, so they are assertions.
        assert_eq!(
            along_x_ids, ordered_x,
            "{name} {n}: the packed and grid-id walks of the same layout disagree"
        );
        assert_eq!(
            scalar_ids, along_x_ids,
            "{name} {n}: the along-x mirror disagrees with the eight-corner definition"
        );
        assert_eq!(
            shipped_ids, along_x_ids,
            "{name} {n}: the mirror's traversal is not the shipped extractor's"
        );
        assert!(
            wrapper_transparent,
            "{name} {n}: the recording wrapper moved the mesh"
        );

        // C3, recorded rather than asserted: it is a registered clause and a
        // clause that panics cannot be reported.
        let mesh_identical = restored == along_x_ids;
        let active_cells = along_x_ids.len();

        // ---- the ten arms. Each owns its own scratch, so all ten can borrow
        // `grid.values` at once and none can see another's bitmap.
        let mut s0: (Vec<u64>, Vec<u32>) = (Vec::new(), Vec::new());
        let mut s1: (Vec<u64>, Vec<u32>) = (Vec::new(), Vec::new());
        let mut s2: (Vec<u64>, Vec<u32>) = (Vec::new(), Vec::new());
        let mut s3: (Vec<u64>, Vec<u32>) = (Vec::new(), Vec::new());
        let mut s4: (Vec<u64>, Vec<u32>) = (Vec::new(), Vec::new());
        let mut s5: (Vec<u64>, Vec<u32>) = (Vec::new(), Vec::new());
        let mut s6: (Vec<u64>, Vec<u32>) = (Vec::new(), Vec::new());
        let mut s7: (Vec<u64>, Vec<u32>, Buckets) = (Vec::new(), Vec::new(), Buckets::new());
        s7.2.fit(c);
        let values = &grid.values;

        let counted = measure_arms(
            probe,
            &mut [
                // along_x: build
                &mut || {
                    black_box(build_along_x(&mut s0.0, values, nn, nn));
                    black_box(&s0.1);
                },
                // along_x: build + fold
                &mut || {
                    let br = build_along_x(&mut s1.0, values, nn, nn);
                    black_box(fold_along_x(&s1.0, br, nn, c));
                    black_box(&s1.1);
                },
                // along_x: build + fused fold and walk
                &mut || {
                    let br = build_along_x(&mut s2.0, values, nn, nn);
                    walk_along_x(&s2.0, br, nn, c, &mut s2.1);
                    black_box(&s2.1);
                },
                // along_x: build + walk + order restore (the identity)
                &mut || {
                    let br = build_along_x(&mut s3.0, values, nn, nn);
                    walk_along_x_ordered(&s3.0, br, nn, c, &mut s3.1);
                    black_box(&s3.1);
                },
                // interleaved: build
                &mut || {
                    black_box(build_interleaved(&mut s4.0, values, nn, nn));
                    black_box(&s4.1);
                },
                // interleaved: build + fold
                &mut || {
                    let bn = build_interleaved(&mut s5.0, values, nn, nn);
                    black_box(fold_interleaved(&s5.0, bn, cb, c));
                    black_box(&s5.1);
                },
                // interleaved: build + fused fold and walk
                &mut || {
                    let bn = build_interleaved(&mut s6.0, values, nn, nn);
                    walk_interleaved(&s6.0, bn, cb, c, &mut s6.1);
                    black_box(&s6.1);
                },
                // interleaved: build + walk + order restore
                &mut || {
                    let bn = build_interleaved(&mut s7.0, values, nn, nn);
                    restore_interleaved(&s7.0, bn, cb, c, &mut s7.2, &mut s7.1);
                    black_box(&s7.1);
                },
            ],
        );

        // ---- the M-287 control, its own batch: it re-samples the field, so a
        // pass costs orders of magnitude more than a bitmap pass and one shared
        // `inner` across all twelve arms would starve the cheap ones.
        let mut unpadded = M287Scratch::new();
        let mut padded = M287Scratch::new();
        let m287_row = nn | 1;
        let m287 = measure_arms(
            probe,
            &mut [
                &mut || {
                    black_box(m287_scaffold(&mut unpadded, field, n, origin, cell, nn));
                },
                &mut || {
                    black_box(m287_scaffold(&mut padded, field, n, origin, cell, m287_row));
                },
            ],
        );

        Measured {
            field: name,
            resolution: n,
            samples: nn * nn * nn,
            cells: c * c * c,
            active_cells,
            along_x: Layout {
                build: counted[0],
                fold: counted[1],
                stage_total: counted[2],
                ordered: counted[3],
            },
            interleaved: Layout {
                build: counted[4],
                fold: counted[5],
                stage_total: counted[6],
                ordered: counted[7],
            },
            along_x_words: inside.len(),
            interleaved_words: blocks.len(),
            m287_unpadded: m287[0],
            m287_padded: m287[1],
            m287_row,
            mesh_identical,
            matches_shipped: shipped_ids == along_x_ids,
            matches_scalar: scalar_ids == along_x_ids,
            wrapper_transparent,
        }
    }

    // ------------------------------------------------------------------- run

    /// A tax: this row's cycles per sample over the 127³ row's, same field,
    /// same arm.
    fn tax(here: f64, baseline: f64) -> f64 {
        if baseline > 0.0 { here / baseline } else { 0.0 }
    }

    #[allow(
        clippy::too_many_lines,
        reason = "one registered row is one `record` call, and the extras M-273 encourages are \
                  what make the verdicts auditable"
    )]
    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        let mut probe = Probe::open();
        let mut measured: Vec<Measured> = Vec::new();

        isomesh::for_each_reference_field!(f32, |name, field| {
            for n in RESOLUTIONS {
                let (_, origin, cell) = crate::common::grid::<f32, _>(&field, n);
                let row = measure(name, &field, n, origin, cell, &mut probe);
                println!(
                    "{name:>20} {n:>4}  active {:>7.4}%  along_x {:>7.3} c/s   \
                     interleaved {:>7.3} c/s   m287 {:>8.3} c/s",
                    row.active_cells as f64 / row.cells as f64 * 100.0,
                    row.along_x.stage_total.cycles / row.samples as f64,
                    row.interleaved.stage_total.cycles / row.samples as f64,
                    row.m287_unpadded.cycles / row.samples as f64,
                );
                measured.push(row);
            }
        });

        // `M-287`'s own fixture, added as a ninth field because the control is
        // registered against a **surface-free** measurement and a tax read off
        // a field with a surface is a different quantity.
        // `experiment_p40.rs:84-89`: the canonical sphere sampled from
        // `[10, 10, 10]`, so no corner is ever inside.
        {
            let sphere = Sphere::<f32>::canonical();
            let (lo, hi) = sphere.domain();
            for n in RESOLUTIONS {
                let cell = (hi[0] - lo[0]) / ((n - 1) as f32);
                let row = measure(CONTROL_FIELD, &sphere, n, [10.0; 3], cell, &mut probe);
                assert_eq!(
                    row.active_cells, 0,
                    "the surface-free control has a surface, so it is not M-287's fixture"
                );
                println!(
                    "{CONTROL_FIELD:>20} {n:>4}  active  0.0000%  along_x {:>7.3} c/s   \
                     interleaved {:>7.3} c/s   m287 {:>8.3} c/s",
                    row.along_x.stage_total.cycles / row.samples as f64,
                    row.interleaved.stage_total.cycles / row.samples as f64,
                    row.m287_unpadded.cycles / row.samples as f64,
                );
                measured.push(row);
            }
        }

        // Baselines for every tax: the 127³ row of the same field.
        let baseline_of = |field: &str| -> &Measured {
            measured
                .iter()
                .find(|m| m.field == field && m.resolution == TAX_BASELINE)
                .expect("every field has a 127 row")
        };

        // ---- the VACUITY CONTROL, scored exactly where it is registered: the
        // unmodified layout's `tax_vs_127` at 128, which must reproduce
        // `M-287`'s 3.37x to within 10%.
        let spike = |field: &str| -> &Measured {
            measured
                .iter()
                .find(|m| m.field == field && m.resolution == TAX_SPIKE)
                .expect("every field has a 128 row")
        };
        let control = spike(CONTROL_FIELD);
        let control_base = baseline_of(CONTROL_FIELD);
        let bitmap_control_tax = tax(
            control.along_x.stage_total.cycles / control.samples as f64,
            control_base.along_x.stage_total.cycles / control_base.samples as f64,
        );
        let scaffold_control_tax = tax(
            control.m287_unpadded.cycles / control.samples as f64,
            control_base.m287_unpadded.cycles / control_base.samples as f64,
        );
        let scaffold_padded_tax = tax(
            control.m287_padded.cycles / control.samples as f64,
            control_base.m287_padded.cycles / control_base.samples as f64,
        );
        let reproduced = |t: f64| (t - M287_TAX).abs() <= M287_TAX * M287_TOLERANCE;
        let control_reproduced = reproduced(bitmap_control_tax);

        println!(
            "\nVACUITY CONTROL on {CONTROL_FIELD} — M-287's own fixture, f32, unpadded row:\n  \
             bitmap stage   tax_vs_127 at 128 = {bitmap_control_tax:.4}x  \
             (M-287: {M287_TAX:.2}x, {} within 10%)\n  \
             m287 scaffold  tax_vs_127 at 128 = {scaffold_control_tax:.4}x  \
             ({} within 10%)\n  \
             m287 scaffold, padded row {} = {scaffold_padded_tax:.4}x\n",
            if control_reproduced { "IS" } else { "is NOT" },
            if reproduced(scaffold_control_tax) {
                "IS"
            } else {
                "is NOT"
            },
            control.m287_row,
        );

        // C1: the interleaved layout's tax at 128 on M-287's own fixture,
        // scored only if the control saw the tax it is supposed to remove.
        let c1_tax = tax(
            control.interleaved.stage_total.cycles / control.samples as f64,
            control_base.interleaved.stage_total.cycles / control_base.samples as f64,
        );
        let c1_verdict = if control_reproduced {
            (c1_tax < C1_BAR).to_string()
        } else {
            String::from("vacuous")
        };

        // C2: the fold's instructions per cell, interleaved against along_x,
        // on the same fixture C1 is scored at. Instruction counts are the
        // deterministic form and carry this verdict.
        let c2_ratio_of = |m: &Measured| -> f64 {
            let along = m.along_x.fold_only().instructions / m.cells as f64;
            let inter = m.interleaved.fold_only().instructions / m.cells as f64;
            if along > 0.0 { inter / along } else { 0.0 }
        };
        let c2_ratio = c2_ratio_of(control);
        let c2_verdict = (c2_ratio - 1.0).abs() <= C2_TOLERANCE;

        // C3: every row's restored interleaved list must be the along_x list.
        let c3_verdict = measured.iter().all(|m| m.mesh_identical);

        println!(
            "C1 tax_vs_127 at 128, interleaved, on {CONTROL_FIELD} = {c1_tax:.4}x \
             (bar {C1_BAR}) -> c1_holds = {c1_verdict}\n\
             C2 fold instruction ratio = {c2_ratio:.4} (bar 1 +/- {C2_TOLERANCE}) \
             -> c2_holds = {c2_verdict}\n\
             C3 restored interleaved list == along_x list on every row -> \
             c3_holds = {c3_verdict}\n"
        );

        for m in &measured {
            let base = baseline_of(m.field);
            let samples = m.samples as f64;
            let base_samples = base.samples as f64;
            let cells = m.cells as f64;

            for (layout, here, there) in [
                ("along_x", &m.along_x, &base.along_x),
                ("interleaved", &m.interleaved, &base.interleaved),
            ] {
                let cps = here.stage_total.cycles / samples;
                let ips = here.stage_total.instructions / samples;
                let l1ps = here.stage_total.l1d_misses / samples;
                let nps = here.stage_total.nanos / samples;
                let fold = here.fold_only();
                let restore = here.restore();
                let base_cps = there.stage_total.cycles / base_samples;

                // `resort_ns_per_vertex` has no denominator on a field with no
                // vertices, and inventing one would be a fabricated number
                // (`M-44`). `unavailable` is `common/mod.rs:15`'s precedent.
                let per_vertex = |value: f64| -> String {
                    if m.active_cells == 0 {
                        String::from("unavailable")
                    } else {
                        format!("{:.6}", value / m.active_cells as f64)
                    }
                };

                let bitmap_share = if m.m287_unpadded.cycles > 0.0 {
                    here.stage_total.cycles / m.m287_unpadded.cycles
                } else {
                    0.0
                };
                let this_tax = tax(cps, base_cps);
                let saved = bitmap_share * (1.0 - 1.0 / this_tax.max(1e-9));
                let ceiling = if this_tax > 1.0 && saved < 1.0 {
                    1.0 / (1.0 - saved)
                } else {
                    1.0
                };

                run.record(&[
                    // ---------------------------------------- registered
                    ("field", m.field.to_string()),
                    ("resolution", m.resolution.to_string()),
                    ("layout", layout.to_string()),
                    ("cycles_per_sample", format!("{cps:.6}")),
                    ("tax_vs_127", format!("{this_tax:.4}")),
                    (
                        "instructions_per_cell",
                        format!("{:.6}", fold.instructions / cells),
                    ),
                    ("mesh_identical", m.mesh_identical.to_string()),
                    ("resort_ns_per_vertex", per_vertex(restore.nanos)),
                    ("c1_holds", c1_verdict.clone()),
                    ("c2_holds", c2_verdict.to_string()),
                    ("c3_holds", c3_verdict.to_string()),
                    // ---------------------------------------- M-280: the clock
                    ("ghz", format!("{:.4}", cps / nps.max(1e-12))),
                    ("ns_per_sample", format!("{nps:.6}")),
                    ("instructions_per_sample", format!("{ips:.6}")),
                    ("l1d_read_misses_per_sample", format!("{l1ps:.6}")),
                    (
                        "tax_vs_127_instructions",
                        format!(
                            "{:.4}",
                            tax(ips, there.stage_total.instructions / base_samples)
                        ),
                    ),
                    (
                        "tax_vs_127_l1d_misses",
                        format!(
                            "{:.4}",
                            tax(l1ps, there.stage_total.l1d_misses / base_samples)
                        ),
                    ),
                    (
                        "tax_vs_127_ns",
                        format!("{:.4}", tax(nps, there.stage_total.nanos / base_samples)),
                    ),
                    // ------------------------------- the prefix decomposition
                    (
                        "cycles_per_sample_build",
                        format!("{:.6}", here.build.cycles / samples),
                    ),
                    (
                        "cycles_per_sample_build_and_fold",
                        format!("{:.6}", here.fold.cycles / samples),
                    ),
                    (
                        "cycles_per_sample_ordered",
                        format!("{:.6}", here.ordered.cycles / samples),
                    ),
                    (
                        "instructions_per_sample_build",
                        format!("{:.6}", here.build.instructions / samples),
                    ),
                    (
                        "cycles_per_cell_fold",
                        format!("{:.6}", fold.cycles / cells),
                    ),
                    (
                        "l1d_read_misses_per_cell_fold",
                        format!("{:.6}", fold.l1d_misses / cells),
                    ),
                    (
                        "cycles_per_cell_walk",
                        format!("{:.6}", stage(here.stage_total, here.fold).cycles / cells),
                    ),
                    (
                        "instructions_per_cell_walk",
                        format!(
                            "{:.6}",
                            stage(here.stage_total, here.fold).instructions / cells
                        ),
                    ),
                    (
                        "instructions_per_cell_stage",
                        format!("{:.6}", here.stage_total.instructions / cells),
                    ),
                    ("resort_cycles_per_vertex", per_vertex(restore.cycles)),
                    (
                        "resort_instructions_per_vertex",
                        per_vertex(restore.instructions),
                    ),
                    // ------------------------------------------- the structure
                    ("samples", m.samples.to_string()),
                    ("cells", m.cells.to_string()),
                    ("active_cells", m.active_cells.to_string()),
                    (
                        "active_fraction",
                        format!("{:.8}", m.active_cells as f64 / cells),
                    ),
                    (
                        "bitmap_words",
                        if layout == "along_x" {
                            m.along_x_words.to_string()
                        } else {
                            m.interleaved_words.to_string()
                        },
                    ),
                    (
                        "bitmap_bytes",
                        if layout == "along_x" {
                            (m.along_x_words * 8).to_string()
                        } else {
                            (m.interleaved_words * 8).to_string()
                        },
                    ),
                    ("values_row_stride", m.samples.to_string().len().to_string()),
                    ("values_row_samples", m.resolution.to_string()),
                    (
                        "values_plane_bytes",
                        (m.resolution as usize * m.resolution as usize * 4).to_string(),
                    ),
                    // ------------------------------------ the M-287 control
                    (
                        "m287_cycles_per_sample_unpadded",
                        format!("{:.6}", m.m287_unpadded.cycles / samples),
                    ),
                    (
                        "m287_cycles_per_sample_padded",
                        format!("{:.6}", m.m287_padded.cycles / samples),
                    ),
                    (
                        "m287_tax_vs_127_unpadded",
                        format!(
                            "{:.4}",
                            tax(
                                m.m287_unpadded.cycles / samples,
                                base.m287_unpadded.cycles / base_samples
                            )
                        ),
                    ),
                    (
                        "m287_tax_vs_127_padded",
                        format!(
                            "{:.4}",
                            tax(
                                m.m287_padded.cycles / samples,
                                base.m287_padded.cycles / base_samples
                            )
                        ),
                    ),
                    (
                        "m287_l1d_misses_per_sample_unpadded",
                        format!("{:.6}", m.m287_unpadded.l1d_misses / samples),
                    ),
                    (
                        "m287_l1d_misses_per_sample_padded",
                        format!("{:.6}", m.m287_padded.l1d_misses / samples),
                    ),
                    ("m287_padded_row_samples", m.m287_row.to_string()),
                    (
                        "bitmap_share_of_m287_scaffold",
                        format!("{bitmap_share:.6}"),
                    ),
                    ("c1_ceiling_at_128", format!("{ceiling:.5}")),
                    // -------------------------------------- the verdict inputs
                    ("c1_bar", C1_BAR.to_string()),
                    ("c1_field", CONTROL_FIELD.to_string()),
                    ("c1_tax_vs_127_interleaved", format!("{c1_tax:.4}")),
                    ("c1_vacuity_control_tax", format!("{bitmap_control_tax:.4}")),
                    (
                        "c1_vacuity_control_reproduced",
                        control_reproduced.to_string(),
                    ),
                    ("c1_vacuity_control_target", M287_TAX.to_string()),
                    (
                        "m287_scaffold_control_tax",
                        format!("{scaffold_control_tax:.4}"),
                    ),
                    (
                        "m287_scaffold_control_reproduced",
                        reproduced(scaffold_control_tax).to_string(),
                    ),
                    (
                        "m287_scaffold_control_tax_padded",
                        format!("{scaffold_padded_tax:.4}"),
                    ),
                    ("c2_tolerance", C2_TOLERANCE.to_string()),
                    ("c2_instruction_ratio_fold", format!("{c2_ratio:.4}")),
                    (
                        "c2_instruction_ratio_fold_this_row",
                        format!("{:.4}", c2_ratio_of(m)),
                    ),
                    (
                        "c2_cycle_ratio_fold_this_row",
                        format!("{:.4}", {
                            let a = m.along_x.fold_only().cycles;
                            let i = m.interleaved.fold_only().cycles;
                            if a > 0.0 { i / a } else { 0.0 }
                        }),
                    ),
                    ("c3_list_matches_shipped", m.matches_shipped.to_string()),
                    ("c3_list_matches_scalar", m.matches_scalar.to_string()),
                    (
                        "c3_recording_wrapper_transparent",
                        m.wrapper_transparent.to_string(),
                    ),
                    // ------------------------------------------- the popcount
                    (
                        "target_feature_popcnt",
                        cfg!(target_feature = "popcnt").to_string(),
                    ),
                    ("count_ones_calls_per_counted_window", 0.to_string()),
                    ("reps", REPS.to_string()),
                ]);
            }
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-104");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent. Every quantity here is a cycle, an
    // instruction or an L1D miss, and all three come from `perf_event_open`.
    // C1 is a cache-set aliasing tax, which an instruction count cannot see
    // and a wall clock cannot carry a verdict on (`M-280`, `M-281`), so off
    // Linux there is nothing to degrade to and a recorded zero would be a
    // fabricated measurement.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} scores a cache-set aliasing tax in cycles per sample, which needs \
             `perf_event_open` and which this platform does not have.",
            prereg.id
        );
        std::process::exit(1);
    }
}
