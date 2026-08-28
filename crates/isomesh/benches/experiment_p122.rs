//! **P-122 — Stream VByte's control/data split, applied to the case stream.**
//!
//! Ticket: R-122. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p122
//! ```
//!
//! Writes `docs/experiments/p-122.csv`. **Linux only**, for `experiment_p12`'s
//! reason and the registration's own: C1 is a branch-misprediction rate, the
//! only instrument that can see one is `perf_event_open`, and a recorded zero
//! off Linux would be a fabricated measurement rather than a missing one.
//!
//! # What was missing
//!
//! Stream VByte (Lemire, Kurz & Rupp, `10.1016/j.ipl.2017.09.011`, acquired in
//! `M-415`) splits a stream of integers into a **dense control stream** — two
//! bits per integer, four integers to a control byte — and a **variable-length
//! data stream**. The decoder reads one control byte, looks its total length up
//! in a table, and advances the data pointer by that length. The branch on
//! *how long this record is* becomes a table lookup, so the decoder never
//! branches on payload.
//!
//! Its headline figure is **4 billion integers per second on a 3.4 GHz
//! Haswell**, and that is the **SIMD** path — a `_mm_shuffle_epi8` gather over
//! sixteen bytes at a time. This harness has no SIMD arm, and
//! `crates/isomesh/src/**` contains no `core::simd`, no intrinsics and no
//! `unsafe`. **That figure is therefore not a comparand and nothing here is
//! scored against it.** Only the *layout* transfers.
//!
//! The layout maps onto Marching Cubes exactly. The shipped cell body
//! (`marching_cubes/mod.rs:254-380`) computes the eight-sign case index and
//! then, in the same iteration of the same loop, reads the case table, probes
//! the edge cache, interpolates the crossings and pushes the triangles.
//! Classification and payload are one interleaved stream. The case index *is* a
//! control byte — 256 values, one per cell, dense, computed for every cell
//! whether or not it produces anything — and the triangles are the
//! variable-length payload keyed by it.
//!
//! Nobody had measured what the interleaving costs in mispredictions. `M-279`
//! reports instructions per sample for whole extractions; `P-121` reports cycles
//! and instructions per stage; neither reads `BRANCH_MISSES`, and there was no
//! per-cell misprediction number anywhere in the repository. This harness
//! produces one, for both layouts, on the same grid in the same run.
//!
//! # The two arms
//!
//! Both arms are bench-local — `crates/isomesh/src/**` is read-only for Phase
//! 25 — and both call the *same two functions* for the two halves of the work,
//! so the only difference between them is the loop structure and the stream
//! layout:
//!
//! - **`single_stream`** is the shipped shape. One triple loop over cells;
//!   per cell [`classify_cell`] then [`March::emit_cell`], the second inlined
//!   into the first's loop body exactly as `mod.rs:259-377` fuses them. The
//!   case byte never leaves a register.
//! - **`split_stream`** is Stream VByte's shape. Pass one runs
//!   [`classify_cell`] over every cell and writes the case byte to a dense
//!   `Vec<u8>` — the control stream, one byte per cell, 256 KiB at 65³ and 2
//!   MiB at 129³. Pass two scans that stream, reads each case byte's payload
//!   length out of a 256-byte table ([`March::payload_len`], Stream VByte's
//!   control-byte length lookup, derived once from `table::CASES`), skips the
//!   cell if the length is zero, and otherwise calls the same
//!   [`March::emit_cell`] to append to the payload stream.
//!
//! The two passes are separated by a [`std::hint::black_box`] on the control
//! stream, so the split is a real two-pass structure and not something the
//! optimiser may fuse back into one.
//!
//! Vertex creation order is identical by construction: both arms visit cells in
//! `z, y, x` order and both allocate a cell's centroids before its triangle
//! vertices, so `positions[i]` is the same vertex in both. That is what makes
//! C3 a bit-for-bit equality rather than a tolerance.
//!
//! # The corner gather is held constant, and three columns prove it
//!
//! `P-121` measured the shipped `MarchingCubes` at **15–24% more** than a
//! bit-identical classify-then-compact mirror, and named the mechanism:
//! `mod.rs:259-268` gathers all eight corner **values** into `[R; 8]` for every
//! cell before it knows whether the cell produces anything, at a ~1.8% active
//! fraction. The split arm here is a near relative of that mirror. A harness
//! whose control pass read only the eight *signs* would harvest `P-121`'s saving
//! and report it as the control/data split's — the row would measure the wrong
//! mechanism and the number would still look plausible.
//!
//! Three things stop that, and the third is checkable from the CSV alone:
//!
//! 1. **One gather, one function, both arms.** [`classify_cell`] returns
//!    `(u8, [R; 8])` and is the only place in this file where a corner value is
//!    read for classification. There is no sign-only path to take. The split
//!    arm's control pass discards the eight values; it does not decline to load
//!    them.
//! 2. **Both arms gather over every cell.** The control pass walks all `(n−1)³`
//!    cells, not a compacted list, exactly as the shipped loop does.
//! 3. **The split arm gathers strictly *more* often, and the count is a
//!    column.** `gather_calls_single` is `cells`; `gather_calls_split` is
//!    `cells + emitting_cells`, because the payload pass re-reads the eight
//!    corner values of each cell it is about to emit rather than carrying an
//!    `[R; 8]` per cell across the pass boundary (which would be 8 MiB of
//!    scratch at 65³ and a different experiment). Both are exact integers,
//!    `gather_calls_ratio > 1` is asserted on every row, and a mechanism cannot
//!    be credited with a saving on a stage it runs more times.
//!
//! One honest caveat rather than a hidden one: the control pass is now a
//! straight loop, so the optimiser is free to do things to it that it cannot do
//! to the fused loop. Any such win is a *consequence of the split* — it is the
//! mechanism, "the classifier never branches on payload" — and it shows up in
//! `instructions_per_cell_split`, which C2 bounds. It is not a corner-gather
//! saving, because clause 3 above holds.
//!
//! # Where the window boundary falls, and why
//!
//! - **Field sampling is outside both windows.** The values array is filled once
//!   per row, before anything is counted, and both arms read the same one.
//!   `P-121` measured field evaluation at **94% of extraction on fbm_terrain**;
//!   counting it would drive every instruction ratio to 1.00 and C2's 10% bar
//!   would then hold by dilution rather than by measurement. The registration's
//!   "total instructions per cell" is a rate over the march, which is what the
//!   split restructures.
//! - **The march is inside, whole.** Edge-cache preparation, classification,
//!   table lookup, cache probe, crossing interpolation, `unit_gradient` normals,
//!   triangle push — the shipped `extract` minus `sample_grid`. No stage of it
//!   is hoisted out to flatter either arm.
//! - `edge_vertices` is cleared and resized inside both windows, as
//!   `mod.rs:250-251` does on every extract. It is a 3.3 MB fill at 65³ and 25.8
//!   MB at 129³, which is real time — but as a *branch-free vectorised fill* it
//!   is under half an instruction per cell and no mispredictions at all, so it
//!   neither dilutes C2 nor flatters C1.
//! - **Windows are siblings, never nested.** `R-121` paid for that discovery:
//!   Zen 3 has six general-purpose counters and `Probe` opens six, so two nested
//!   windows multiplex and `Probe::worst_ratio` refuses. One window per arm per
//!   repetition, and the arms **alternate order by repetition parity** so
//!   neither is permanently the one that runs second.
//!
//! # SHARE
//!
//! Each clause's reachable share, as a column:
//!
//! - **C2 is reported first, because it is the deterministic clause.** `M-280`
//!   and `M-281`: on a governed CPU a nanosecond is not a unit, and `R-105`
//!   watched one binary's cycle ratio band drift from 0.984 to 1.035 across
//!   three runs while its instruction counts held to four figures. C2's share is
//!   `instructions_per_cell_single` — the whole march's instruction rate, the
//!   thing a 10% rise is 10% of — and the bar is `instruction_ratio ≤ 1.10`.
//!   Scored on the **worst row**, which is the tight reading: dilution by an
//!   expensive field can only push a ratio toward 1.0, never above the maximum.
//!   `worst_instruction_ratio` is a column, and so is
//!   `instruction_ratio_rep_spread`, which *demonstrates* the determinism
//!   instead of asserting it.
//! - **C1's share is `branch_misses_per_cell_single`**, the misprediction rate
//!   of the shipped shape, which is the entire quantity available to remove. It
//!   must be non-zero — that is the registered vacuity control, and it is
//!   asserted per row rather than merely recorded, because a ratio against a
//!   floor of zero is not a measurement. C1 is scored **per row**, and
//!   `c1_holds_all_rows` is the aggregate: a mechanism that removes
//!   mispredictions on some fields and not others is a split verdict, not a
//!   held one, and the CSV should be able to say so.
//!
//!   **A branch counter is not an instruction counter.** The compact fields
//!   mispredict about 0.008–0.029 times per cell, so the median of [`REPS`]
//!   repetitions carries the number and four columns report how firm it is:
//!   `branch_miss_ratio_rep_best` and `branch_miss_ratio_rep_worst` for the
//!   spread, and `branch_miss_ratio_single_first` /
//!   `branch_miss_ratio_split_first` for the median **within each ordering** of
//!   the two sibling windows. The second pair is the one that matters: the arm
//!   that runs second inherits the other's cache and predictor state, and a row
//!   whose two orderings disagree about the sign of the effect is a row nobody
//!   should quote as a firm reduction. The printed summary counts them.
//!
//!   C1's falsifier is *"no branch-misprediction reduction, measured with
//!   counters"*, and a measurement that only says "no" is worth less than one
//!   that says where the branches went. So `extra_branch_misses_per_pass` is
//!   the signed difference and `extra_branch_misses_per_x_row` divides it by
//!   `x_rows`, the innermost-loop exits in one traversal of the cell space. The
//!   split arm traverses that space **twice** — once to write the control
//!   stream, once to read it — so it pays `x_rows` extra loop exits, and a
//!   column reading near 1.0 attributes the whole difference to them rather
//!   than leaving it a mystery.
//! - **C3's share is the mesh**, and it is an equality rather than a share:
//!   positions, normals and indices compared as bit patterns, in creation order,
//!   between the two arms. `mesh_identical_to_shipped` is beside it — the
//!   single-stream arm against `MarchingCubes::extract` on the same grid in the
//!   same run — which is `M-279`'s rule that a new instrument's first job is to
//!   agree with the old one where they overlap, and the licence for calling the
//!   single-stream arm "the shipped shape" at all. Both are asserted.
//! - **No time claim is made.** `ns_per_cell` and `ghz` are on every row as
//!   provenance so a later reader can see what clock the counts were taken at,
//!   and no clause consults either.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::hint::black_box;
    use std::time::Instant;

    use isomesh::marching_cubes::MarchingCubes;
    use isomesh::marching_cubes::table::{
        CASES, CENTROID_BASE, EDGE_AXIS, EDGE_CORNERS, MAX_CENTROIDS, edge_offset, is_centroid,
        is_inside,
    };
    use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    // ─── the fixture ───────────────────────────────────────────────────────

    /// The registered resolutions, in samples per axis.
    const RESOLUTIONS: [u32; 2] = [65, 129];
    /// Measured repetitions per arm per row, and **even on purpose**: the arms
    /// alternate which of them runs second (see [`measure`]), so an even count
    /// gives each ordering exactly half the repetitions and the median ratio
    /// cannot be decided by whichever ordering happened to get one more.
    ///
    /// The median per quantity carries the verdict; the best and worst
    /// per-repetition branch-miss ratios, and the median within each ordering,
    /// are reported beside it. That matters here: at 65³ on the compact fields
    /// the single-stream arm mispredicts about **0.010–0.029 times per cell**,
    /// so one window holds a few tens of thousands of events and a branch
    /// counter — unlike an instruction counter — is not deterministic.
    const REPS: usize = 16;
    /// Untimed passes of both arms before anything is counted, so that the
    /// buffers are at final capacity and the pages are faulted in.
    const WARMUP: usize = 2;
    /// How long one counter window should last, in nanoseconds.
    ///
    /// Longer than `experiment_p121`'s 30 ms, and for a reason that is specific
    /// to this row rather than a matter of taste: the quantity C1 reads is
    /// two orders of magnitude rarer than a cycle, so the window has to hold
    /// enough mispredictions for the ratio to mean something. At 30 ms the
    /// per-repetition band on `box_exact 65³` spanned 0.995 to 2.296.
    const TARGET_BATCH_NS: f64 = 100_000_000.0;
    /// Ceiling on the batch, so a cheap row cannot take minutes.
    const MAX_INNER: usize = 1024;
    /// C2's bar, from the registration: instructions per cell may rise by 10%.
    const INSTRUCTION_CEILING: f64 = 1.10;

    // ─── private crate mechanisms, copied rather than made `pub` ────────────

    /// `cube::corner_offset`. Private, and `src/**` is read-only this phase.
    #[inline]
    const fn corner_offset(corner: u8) -> [u32; 3] {
        [
            (corner & 1) as u32,
            ((corner >> 1) & 1) as u32,
            ((corner >> 2) & 1) as u32,
        ]
    }

    /// `cube::place`: the centred frame, spelled once in the crate and once
    /// here.
    #[inline]
    fn place<R: Real>(lo: R, hi: R, d: R) -> R {
        (lo + hi) * R::HALF + (hi - lo) * d
    }

    /// `vec3::length`. Left-to-right summation, which is what makes the mirror's
    /// normals bit-identical to the crate's rather than merely close.
    #[inline]
    fn length<R: Real>(a: [R; 3]) -> R {
        (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
    }

    /// `marching_cubes::unit_gradient`, over `vec3::scale`.
    #[inline]
    fn unit_gradient<R: Real, S: Sdf<Scalar = R>>(sdf: &S, position: [R; 3]) -> [R; 3] {
        let g = sdf.gradient(position);
        let inv = length(g).recip();
        [g[0] * inv, g[1] * inv, g[2] * inv]
    }

    /// `marching_cubes::corner_position`.
    #[inline]
    fn corner_position<R: Real>(
        base: [u32; 3],
        corner: u8,
        origin: [R; 3],
        cell_size: R,
    ) -> [R; 3] {
        let o = corner_offset(corner);
        [
            origin[0] + cell_size * R::from_f64(f64::from(base[0] + o[0])),
            origin[1] + cell_size * R::from_f64(f64::from(base[1] + o[1])),
            origin[2] + cell_size * R::from_f64(f64::from(base[2] + o[2])),
        ]
    }

    /// `marching_cubes::edge_position` with `crossing_refinement == 0`, which is
    /// `MarchingCubes::new`'s default: `refine_crossing` returns `d0` unchanged
    /// at zero steps, so there is no field evaluation on this path.
    #[inline]
    fn edge_position<R: Real>(
        base: [u32; 3],
        edge: u8,
        corner_value: &[R; 8],
        origin: [R; 3],
        cell_size: R,
    ) -> [R; 3] {
        let [lo_corner, hi_corner] = EDGE_CORNERS[edge as usize];
        let d = edge_offset(
            corner_value[lo_corner as usize],
            corner_value[hi_corner as usize],
        );
        let lo_pos = corner_position(base, lo_corner, origin, cell_size);
        let hi_pos = corner_position(base, hi_corner, origin, cell_size);
        [
            place(lo_pos[0], hi_pos[0], d),
            place(lo_pos[1], hi_pos[1], d),
            place(lo_pos[2], hi_pos[2], d),
        ]
    }

    // ─── the grid ──────────────────────────────────────────────────────────

    /// One cubic grid: `n` samples per axis, spanning `n − 1` cells.
    #[derive(Clone, Copy)]
    struct Grid<R: Real> {
        n: u32,
        origin: [R; 3],
        cell_size: R,
    }

    impl<R: Real> Grid<R> {
        #[inline]
        fn samples(self) -> usize {
            let n = self.n as usize;
            n * n * n
        }

        #[inline]
        fn cells(self) -> u32 {
            self.n - 1
        }

        #[inline]
        fn cell_count(self) -> usize {
            let c = self.cells() as usize;
            c * c * c
        }

        /// Innermost-loop exits in one traversal of the cell space: one per
        /// `x`-row. Both arms' loops are `z, y, x`, so the split arm — which
        /// traverses twice — pays this many extra loop exits, and
        /// `extra_branch_misses_per_x_row` is what says whether the predictor
        /// absorbs them.
        #[inline]
        fn x_rows(self) -> usize {
            let c = self.cells() as usize;
            c * c
        }

        /// `RuntimeShape3::linearize`, which is the layout Marching Cubes
        /// samples into: no row padding, so the stride is the row itself.
        #[inline]
        fn sample_index(self, p: [u32; 3]) -> usize {
            let n = self.n as usize;
            p[0] as usize + n * (p[1] as usize + n * p[2] as usize)
        }
    }

    // ─── the gather, shared by both arms ───────────────────────────────────

    /// The eight-sign case index and the eight corner values, for one cell.
    ///
    /// `mod.rs:259-268` verbatim, in **one** place, called by both arms — the
    /// whole point being that neither arm has a cheaper way to classify a cell
    /// than the other. `is_inside` is `value < 0`, so an exact zero is outside.
    #[inline]
    fn classify_cell<R: Real>(values: &[R], g: Grid<R>, base: [u32; 3]) -> (u8, [R; 8]) {
        let mut case = 0u8;
        let mut corner_value = [R::ZERO; 8];
        for (c, slot) in corner_value.iter_mut().enumerate() {
            let o = corner_offset(c as u8);
            let v = values[g.sample_index([base[0] + o[0], base[1] + o[1], base[2] + o[2]])];
            *slot = v;
            if is_inside(v) {
                case |= 1 << c;
            }
        }
        (case, corner_value)
    }

    // ─── the march, in two layouts ─────────────────────────────────────────

    /// Both arms, sharing every buffer so that neither gets a different
    /// allocation shape, and sharing [`Self::emit_cell`] so that the payload
    /// code is literally the same instructions in both.
    struct March<R: Real> {
        /// Filled once per row, **outside** every counter window; both arms read
        /// this one array.
        values: Vec<R>,
        /// **The control stream.** One dense case byte per cell. Sized once per
        /// row and overwritten in place, which is what an implementation would
        /// do; the split arm's cost is the write and the read-back, and that
        /// cost belongs on the row.
        cases: Vec<u8>,
        /// **Stream VByte's control-byte length table**: triangles per case, 256
        /// bytes against `CASES`'s several kilobytes, so the control-stream scan
        /// touches four cache lines rather than the whole case table. Derived
        /// once from `table::CASES`, which is a `static` and so cannot be read
        /// in a `const`.
        payload_len: [u8; 256],
        /// `MarchingCubes::edge_vertices`: one `u32` slot per (sample, axis).
        edge_vertices: Vec<u32>,
        /// **The payload stream**, in three parts, exactly as `MeshBuffer` holds
        /// it.
        positions: Vec<[R; 3]>,
        normals: Vec<[R; 3]>,
        indices: Vec<u32>,
    }

    impl<R: Real> March<R> {
        fn new() -> Self {
            let mut payload_len = [0u8; 256];
            for (case, slot) in payload_len.iter_mut().enumerate() {
                *slot = CASES[case].count;
            }
            Self {
                values: Vec::new(),
                cases: Vec::new(),
                payload_len,
                edge_vertices: Vec::new(),
                positions: Vec::new(),
                normals: Vec::new(),
                indices: Vec::new(),
            }
        }

        /// `sdf::sample_grid` with `row_stride == size[0]`, which is what
        /// `MarchingCubes::extract` passes. Outside every window.
        fn sample<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            self.values.clear();
            self.values.reserve(g.samples());
            for z in 0..g.n {
                for y in 0..g.n {
                    for x in 0..g.n {
                        self.values.push(sdf.sample([
                            g.origin[0] + g.cell_size * R::from_f64(f64::from(x)),
                            g.origin[1] + g.cell_size * R::from_f64(f64::from(y)),
                            g.origin[2] + g.cell_size * R::from_f64(f64::from(z)),
                        ]));
                    }
                }
            }
            self.cases.clear();
            self.cases.resize(g.cell_count(), 0);
        }

        /// `mod.rs:250-251`, and it runs inside both arms' windows because it
        /// runs on every shipped extract.
        #[inline]
        fn prepare(&mut self, g: Grid<R>) {
            self.edge_vertices.clear();
            self.edge_vertices.resize(g.samples() * 3, u32::MAX);
            self.positions.clear();
            self.normals.clear();
            self.indices.clear();
        }

        /// `MeshSink::vertex`: append to the payload stream, return its index.
        #[inline]
        fn vertex<S: Sdf<Scalar = R>>(&mut self, sdf: &S, position: [R; 3]) -> u32 {
            let index = self.positions.len() as u32;
            self.positions.push(position);
            self.normals.push(unit_gradient(sdf, position));
            index
        }

        /// `MarchingCubes::vertex_on_edge` — the edge-cache probe, and the only
        /// place a shared vertex is created.
        #[inline]
        fn vertex_on_edge<S: Sdf<Scalar = R>>(
            &mut self,
            sdf: &S,
            g: Grid<R>,
            base: [u32; 3],
            edge: u8,
            corner_value: &[R; 8],
        ) -> u32 {
            let axis = EDGE_AXIS[edge as usize] as usize;
            let o = corner_offset(EDGE_CORNERS[edge as usize][0]);
            let lo_sample = g.sample_index([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
            let key = lo_sample * 3 + axis;
            let cached = self.edge_vertices[key];
            if cached != u32::MAX {
                return cached;
            }
            let position = edge_position(base, edge, corner_value, g.origin, g.cell_size);
            let index = self.vertex(sdf, position);
            self.edge_vertices[key] = index;
            index
        }

        /// **The payload for one cell**, `mod.rs:307-377` under
        /// `MarchingCubes::new`'s defaults — `FaceAmbiguity::Separate`,
        /// `InteriorAmbiguity::Ignore`, `crossing_refinement == 0` — so
        /// `ambiguous` is zero, `mask` is zero and the entry is the derived
        /// table's.
        ///
        /// Called by **both** arms, so the payload code is the same
        /// instructions in both and the only difference between them is where
        /// the case byte came from.
        #[inline]
        fn emit_cell<S: Sdf<Scalar = R>>(
            &mut self,
            sdf: &S,
            g: Grid<R>,
            base: [u32; 3],
            case: u8,
            corner_value: &[R; 8],
        ) {
            let entry = &CASES[case as usize];
            if entry.count == 0 {
                return;
            }
            // Cycle centroids first, because a triangle naming one needs every
            // edge vertex of that cycle averaged before it can be emitted
            // (A-015). Cell-local, so never cached.
            let mut centroid = [0u32; MAX_CENTROIDS];
            for (c, slot) in centroid
                .iter_mut()
                .enumerate()
                .take(entry.centroids as usize)
            {
                let code = CENTROID_BASE + c as u8;
                let mut sum = [R::ZERO; 3];
                let mut n = 0u32;
                for tri in &entry.triangles[..entry.count as usize] {
                    if tri[0] != code {
                        continue;
                    }
                    let p = edge_position(base, tri[1], corner_value, g.origin, g.cell_size);
                    sum = [sum[0] + p[0], sum[1] + p[1], sum[2] + p[2]];
                    n += 1;
                }
                let scale = R::from_f64(f64::from(n)).recip();
                let position = [sum[0] * scale, sum[1] * scale, sum[2] * scale];
                *slot = self.vertex(sdf, position);
            }
            for tri in &entry.triangles[..entry.count as usize] {
                let mut idx = [0u32; 3];
                for (k, &code) in tri.iter().enumerate() {
                    idx[k] = if is_centroid(code) {
                        centroid[(code - CENTROID_BASE) as usize]
                    } else {
                        self.vertex_on_edge(sdf, g, base, code, corner_value)
                    };
                }
                self.indices.extend_from_slice(&idx);
            }
        }

        /// **The single-stream arm.** The shipped shape: classification and
        /// payload interleaved in one loop body, the case byte never leaving a
        /// register.
        fn march_single<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            self.prepare(g);
            let c = g.cells();
            for z in 0..c {
                for y in 0..c {
                    for x in 0..c {
                        let base = [x, y, z];
                        let (case, corner_value) = classify_cell(&self.values, g, base);
                        self.emit_cell(sdf, g, base, case, &corner_value);
                    }
                }
            }
            black_box(&self.indices);
        }

        /// **The split-stream arm.** Stream VByte's shape: the whole dense
        /// control stream first, then a scan of it emitting payload.
        fn march_split<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            self.prepare(g);
            let c = g.cells();

            // ── the control stream: one case byte per cell, no payload ───────
            let mut i = 0usize;
            for z in 0..c {
                for y in 0..c {
                    for x in 0..c {
                        // The eight values are gathered and discarded, not
                        // declined — see the module docs on holding the gather
                        // constant.
                        let (case, _corner_value) = classify_cell(&self.values, g, [x, y, z]);
                        self.cases[i] = case;
                        i += 1;
                    }
                }
            }
            // The split is a two-pass structure, and this is what says so to the
            // optimiser as well as to the reader.
            black_box(&self.cases);

            // ── the data stream: the payload, walked off the control stream ──
            let mut i = 0usize;
            for z in 0..c {
                for y in 0..c {
                    for x in 0..c {
                        let case = self.cases[i];
                        i += 1;
                        // Stream VByte's control-byte length lookup: 256 bytes,
                        // four cache lines, and the only thing the scan reads on
                        // the 98% of cells that emit nothing.
                        if self.payload_len[case as usize] == 0 {
                            continue;
                        }
                        let base = [x, y, z];
                        let (_, corner_value) = classify_cell(&self.values, g, base);
                        self.emit_cell(sdf, g, base, case, &corner_value);
                    }
                }
            }
            black_box(&self.indices);
        }

        /// How many cells the case table gives triangles for. Outside every
        /// window; the denominator of `active_fraction` and the second term of
        /// `gather_calls_split`.
        fn emitting_cells(&self, g: Grid<R>) -> usize {
            let c = g.cells();
            let mut n = 0usize;
            for z in 0..c {
                for y in 0..c {
                    for x in 0..c {
                        let (case, _) = classify_cell(&self.values, g, [x, y, z]);
                        if self.payload_len[case as usize] != 0 {
                            n += 1;
                        }
                    }
                }
            }
            n
        }
    }

    // ─── counting ──────────────────────────────────────────────────────────

    /// What one or more counter windows read.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
        branch_misses: f64,
        cache_misses: f64,
        l1d_read_misses: f64,
    }

    impl Counted {
        fn scaled(self, by: f64) -> Self {
            Self {
                cycles: self.cycles * by,
                instructions: self.instructions * by,
                branch_misses: self.branch_misses * by,
                cache_misses: self.cache_misses * by,
                l1d_read_misses: self.l1d_read_misses * by,
            }
        }
    }

    /// One counter window over `inner` repetitions of `body`, divided by
    /// `inner`.
    ///
    /// The `perf_event` system calls are all **outside** the counted region.
    /// Windows are never nested: `Probe` opens six hardware counters and Zen 3
    /// has six, so a nested pair would multiplex and `worst_ratio` refuses —
    /// `R-121` paid for that discovery.
    fn window(probe: &mut Probe, inner: usize, mut body: impl FnMut()) -> (Counted, f64) {
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
        (
            Counted {
                cycles: counts.cycles.count as f64,
                instructions: counts.instructions.count as f64,
                branch_misses: counts.branch_misses.count as f64,
                cache_misses: counts.cache_misses.count as f64,
                l1d_read_misses: counts.l1d_read_misses.count as f64,
            }
            .scaled(scale),
            nanos * scale,
        )
    }

    /// One repetition: one sibling window per arm, in one of the two orderings.
    #[derive(Clone, Copy, Default)]
    struct Rep {
        single: Counted,
        split: Counted,
        single_ns: f64,
        split_ns: f64,
        /// Which arm ran first. Recorded rather than assumed away: the second
        /// arm of a repetition inherits the first's cache and predictor state,
        /// and a ratio that depends on the ordering is a ratio a reader has to
        /// be told about.
        single_first: bool,
    }

    impl Rep {
        fn branch_miss_ratio(&self) -> f64 {
            self.split.branch_misses / self.single.branch_misses
        }

        fn instruction_ratio(&self) -> f64 {
            self.split.instructions / self.single.instructions
        }
    }

    /// The median of a set of readings, taken **per quantity** rather than per
    /// repetition: one repetition disturbed by another process on the machine
    /// should move one number, not a whole row.
    fn median(pick: &dyn Fn(&Rep) -> f64, reps: &[Rep]) -> f64 {
        let mut values: Vec<f64> = reps.iter().map(pick).collect();
        values.sort_by(|a, b| a.total_cmp(b));
        values[values.len() / 2]
    }

    /// [`median`] of a derived quantity over a subset of the repetitions.
    fn median_of(values: &mut [f64]) -> f64 {
        values.sort_by(|a, b| a.total_cmp(b));
        values[values.len() / 2]
    }

    fn median_counted(pick: &dyn Fn(&Rep) -> Counted, reps: &[Rep]) -> Counted {
        Counted {
            cycles: median(&|r| pick(r).cycles, reps),
            instructions: median(&|r| pick(r).instructions, reps),
            branch_misses: median(&|r| pick(r).branch_misses, reps),
            cache_misses: median(&|r| pick(r).cache_misses, reps),
            l1d_read_misses: median(&|r| pick(r).l1d_read_misses, reps),
        }
    }

    // ─── one row ───────────────────────────────────────────────────────────

    /// Which layout a row reports.
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Arm {
        Single,
        Split,
    }

    impl Arm {
        fn name(self) -> &'static str {
            match self {
                Self::Single => "single_stream",
                Self::Split => "split_stream",
            }
        }
    }

    /// One measured `(field, resolution)`: both arms, and the comparison.
    struct Measured {
        field: &'static str,
        resolution: u32,
        inner: usize,
        cells: usize,
        emitting_cells: usize,
        vertices: usize,
        triangles: usize,
        mesh_identical: bool,
        mesh_identical_to_shipped: bool,
        single: Counted,
        split: Counted,
        single_ns: f64,
        split_ns: f64,
        branch_ratio_best: f64,
        branch_ratio_worst: f64,
        /// The median branch-miss ratio within each ordering. If these two
        /// disagree, the row's verdict depends on which arm ran second and the
        /// reader needs to know before quoting the median.
        branch_ratio_single_first: f64,
        branch_ratio_split_first: f64,
        instruction_ratio_spread: f64,
        x_rows: usize,
    }

    impl Measured {
        fn counted(&self, arm: Arm) -> Counted {
            match arm {
                Arm::Single => self.single,
                Arm::Split => self.split,
            }
        }

        fn nanos(&self, arm: Arm) -> f64 {
            match arm {
                Arm::Single => self.single_ns,
                Arm::Split => self.split_ns,
            }
        }

        fn per_cell(&self, arm: Arm, pick: impl Fn(Counted) -> f64) -> f64 {
            pick(self.counted(arm)) / self.cells as f64
        }

        /// Below 1 is a reduction, which is what C1 asks for.
        fn branch_miss_ratio(&self) -> f64 {
            self.split.branch_misses / self.single.branch_misses
        }

        /// C2's bar is 1.10 on this.
        fn instruction_ratio(&self) -> f64 {
            self.split.instructions / self.single.instructions
        }

        /// What the split arm's extra mispredictions cost per pass. Signed: a
        /// reduction would be negative, and a clause about reducing something
        /// should be able to report the reduction it found.
        fn extra_branch_misses(&self) -> f64 {
            self.split.branch_misses - self.single.branch_misses
        }

        /// [`Self::extra_branch_misses`] per innermost-loop exit of one
        /// traversal. The split arm walks the cell space **twice** — once to
        /// write the control stream and once to read it — so it pays
        /// [`Grid::x_rows`] extra loop exits, and this column says whether the
        /// branch predictor absorbs them or charges for them. Near 1.0 means
        /// the second traversal's own loop exits *are* the difference, which is
        /// a mechanism rather than a mystery.
        fn extra_branch_misses_per_x_row(&self) -> f64 {
            self.extra_branch_misses() / self.x_rows as f64
        }

        fn gather_calls_single(&self) -> usize {
            self.cells
        }

        /// Strictly more than [`Self::gather_calls_single`], which is what stops
        /// this row being credited with `P-121`'s corner-gather saving.
        fn gather_calls_split(&self) -> usize {
            self.cells + self.emitting_cells
        }

        fn c1_holds(&self) -> bool {
            self.branch_miss_ratio() < 1.0
        }

        fn c2_holds(&self) -> bool {
            self.instruction_ratio() <= INSTRUCTION_CEILING
        }

        fn c3_holds(&self) -> bool {
            self.mesh_identical && self.mesh_identical_to_shipped
        }
    }

    /// Bit for bit, as bit patterns rather than as values.
    fn same<R: Real>(a: &[[R; 3]], b: &[[R; 3]]) -> bool {
        a.len() == b.len()
            && a.iter()
                .zip(b)
                .all(|(p, q)| (0..3).all(|k| p[k].as_f64().to_bits() == q[k].as_f64().to_bits()))
    }

    /// Measure one `(field, resolution)`.
    fn measure<R, S>(field: &'static str, n: u32, sdf: &S, origin: [R; 3], cell_size: R) -> Measured
    where
        R: Real,
        S: Sdf<Scalar = R>,
    {
        let g = Grid {
            n,
            origin,
            cell_size,
        };
        let shape = RuntimeShape3::new([n; 3]).expect("the fixture fits u32");

        let mut m = March::<R>::new();
        m.sample(sdf, g);
        let emitting = m.emitting_cells(g);

        for _ in 0..WARMUP {
            m.march_single(sdf, g);
            m.march_split(sdf, g);
        }

        // ── the batch, chosen from a timed pass ──────────────────────────────
        let started = Instant::now();
        m.march_single(sdf, g);
        let pass_ns = started.elapsed().as_nanos() as f64;
        let inner = ((TARGET_BATCH_NS / pass_ns.max(1.0)).ceil() as usize).clamp(1, MAX_INNER);

        // ── C3, and M-279's agreement check ─────────────────────────────────
        m.march_single(sdf, g);
        let single_positions = m.positions.clone();
        let single_normals = m.normals.clone();
        let single_indices = m.indices.clone();
        m.march_split(sdf, g);
        let mesh_identical = same(&single_positions, &m.positions)
            && same(&single_normals, &m.normals)
            && single_indices.as_slice() == m.indices.as_slice();
        assert!(
            mesh_identical,
            "{field} {n}^3: the split-stream arm's mesh differs from the single-stream arm's \
             ({} vs {} vertices, {} vs {} indices) — C3 is an equality and the two arms are \
             supposed to be one algorithm in two layouts",
            single_positions.len(),
            m.positions.len(),
            single_indices.len(),
            m.indices.len()
        );

        let mut shipped = MarchingCubes::<R>::new();
        let mut out = MeshBuffer::<R>::new();
        shipped
            .extract(sdf, &shape, origin, cell_size, &mut out)
            .expect("extraction");
        let mesh_identical_to_shipped = same(&single_positions, &out.positions)
            && same(&single_normals, &out.normals)
            && single_indices.as_slice() == out.indices.as_slice();
        assert!(
            mesh_identical_to_shipped,
            "{field} {n}^3: the single-stream arm's mesh differs from MarchingCubes::extract's \
             ({} vs {} vertices, {} vs {} indices) — then it is not the shipped shape and the \
             comparison below is between two things neither of which ships",
            single_positions.len(),
            out.positions.len(),
            single_indices.len(),
            out.indices.len()
        );

        // ── REPS repetitions, one sibling window per arm ─────────────────────
        let mut probe = Probe::open();
        let mut reps: Vec<Rep> = Vec::with_capacity(REPS);
        for rep in 0..REPS {
            // Alternate which arm runs second, so neither is permanently the
            // one that inherits the other's cache and predictor state.
            let mut r = Rep {
                single_first: rep % 2 == 0,
                ..Rep::default()
            };
            if r.single_first {
                let (c, ns) = window(&mut probe, inner, || m.march_single(sdf, g));
                r.single = c;
                r.single_ns = ns;
                let (c, ns) = window(&mut probe, inner, || m.march_split(sdf, g));
                r.split = c;
                r.split_ns = ns;
            } else {
                let (c, ns) = window(&mut probe, inner, || m.march_split(sdf, g));
                r.split = c;
                r.split_ns = ns;
                let (c, ns) = window(&mut probe, inner, || m.march_single(sdf, g));
                r.single = c;
                r.single_ns = ns;
            }
            reps.push(r);
        }

        let branch_ratios: Vec<f64> = reps.iter().map(Rep::branch_miss_ratio).collect();
        let instruction_ratios: Vec<f64> = reps.iter().map(Rep::instruction_ratio).collect();
        let branch_ratio_best = branch_ratios.iter().copied().fold(f64::MAX, f64::min);
        let branch_ratio_worst = branch_ratios.iter().copied().fold(0.0f64, f64::max);
        let instruction_ratio_spread = instruction_ratios.iter().copied().fold(0.0f64, f64::max)
            - instruction_ratios.iter().copied().fold(f64::MAX, f64::min);
        let mut single_first: Vec<f64> = reps
            .iter()
            .filter(|r| r.single_first)
            .map(Rep::branch_miss_ratio)
            .collect();
        let mut split_first: Vec<f64> = reps
            .iter()
            .filter(|r| !r.single_first)
            .map(Rep::branch_miss_ratio)
            .collect();

        let single = median_counted(&|r| r.single, &reps);
        let split = median_counted(&|r| r.split, &reps);

        // The registration's vacuity control, asserted rather than merely
        // recorded: a ratio against a floor of zero is not a measurement.
        assert!(
            single.branch_misses > 0.0,
            "{field} {n}^3: the single-stream arm read zero branch mispredictions, so there is \
             nothing for the split to remove and `branch_miss_ratio` would be a division by a \
             floor"
        );

        Measured {
            field,
            resolution: n,
            inner,
            cells: g.cell_count(),
            emitting_cells: emitting,
            vertices: single_positions.len(),
            triangles: single_indices.len() / 3,
            mesh_identical,
            mesh_identical_to_shipped,
            single,
            split,
            single_ns: median(&|r| r.single_ns, &reps),
            split_ns: median(&|r| r.split_ns, &reps),
            branch_ratio_best,
            branch_ratio_worst,
            branch_ratio_single_first: median_of(&mut single_first),
            branch_ratio_split_first: median_of(&mut split_first),
            instruction_ratio_spread,
            x_rows: g.x_rows(),
        }
    }

    /// The registered fixture: eight reference fields × {65³, 129³}.
    fn sweep() -> Vec<Measured> {
        let mut rows = Vec::new();
        for n in RESOLUTIONS {
            isomesh::for_each_reference_field!(f32, |name, field| {
                let (_, origin, cell_size) = crate::common::grid(&field, n);
                rows.push(measure(name, n, &field, origin, cell_size));
            });
        }
        rows
    }

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        let rows = sweep();

        let worst_instruction_ratio = rows
            .iter()
            .map(Measured::instruction_ratio)
            .fold(0.0f64, f64::max);
        let worst_branch_miss_ratio = rows
            .iter()
            .map(Measured::branch_miss_ratio)
            .fold(0.0f64, f64::max);
        let best_branch_miss_ratio = rows
            .iter()
            .map(Measured::branch_miss_ratio)
            .fold(f64::MAX, f64::min);
        let reducing = rows.iter().filter(|r| r.c1_holds()).count();
        let c1_all = rows.iter().all(Measured::c1_holds);
        let c2_all = rows.iter().all(Measured::c2_holds);
        let c3_all = rows.iter().all(Measured::c3_holds);

        // ── C2 first: it is the deterministic clause ─────────────────────────
        println!(
            "{:<16} {:>4} {:>5} {:>10} {:>10} {:>7} {:>10} {:>10} {:>7} {:>7} {:>7} {:>7} \
             {:>7} {:>5}",
            "field",
            "n",
            "inner",
            "instr/cell",
            "instr/cell",
            "instr",
            "bmiss/cell",
            "bmiss/cell",
            "bmiss",
            "best",
            "worst",
            "1st=sgl",
            "1st=spl",
            "C1C2"
        );
        println!(
            "{:<16} {:>4} {:>5} {:>10} {:>10} {:>7} {:>10} {:>10} {:>7} {:>7} {:>7} {:>7} \
             {:>7} {:>5}",
            "",
            "",
            "",
            "single",
            "split",
            "ratio",
            "single",
            "split",
            "ratio",
            "rep",
            "rep",
            "ratio",
            "ratio",
            ""
        );
        for r in &rows {
            println!(
                "{:<16} {:>4} {:>5} {:>10.3} {:>10.3} {:>7.4} {:>10.5} {:>10.5} {:>7.4} \
                 {:>7.4} {:>7.4} {:>7.4} {:>7.4} {:>2}{:>2}",
                r.field,
                r.resolution,
                r.inner,
                r.per_cell(Arm::Single, |c| c.instructions),
                r.per_cell(Arm::Split, |c| c.instructions),
                r.instruction_ratio(),
                r.per_cell(Arm::Single, |c| c.branch_misses),
                r.per_cell(Arm::Split, |c| c.branch_misses),
                r.branch_miss_ratio(),
                r.branch_ratio_best,
                r.branch_ratio_worst,
                r.branch_ratio_single_first,
                r.branch_ratio_split_first,
                u8::from(r.c1_holds()),
                u8::from(r.c2_holds()),
            );
        }

        println!(
            "\nC2 instructions per cell: worst ratio over {} rows = {worst_instruction_ratio:.4} \
             (bar {INSTRUCTION_CEILING}) -> {}",
            rows.len(),
            if c2_all { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "    worst per-row spread of the instruction ratio over {REPS} repetitions: \
             {:.6} — M-279's determinism, demonstrated rather than asserted",
            rows.iter()
                .map(|r| r.instruction_ratio_spread)
                .fold(0.0f64, f64::max)
        );
        println!(
            "C1 branch mispredictions per cell: {reducing} of {} rows reduce; ratios span \
             {best_branch_miss_ratio:.4}..{worst_branch_miss_ratio:.4} -> {}",
            rows.len(),
            if c1_all {
                "HELD"
            } else if reducing == 0 {
                "FALSIFIED"
            } else {
                "SPLIT — reduces on some fields and not others"
            }
        );
        println!(
            "    rows whose verdict flips with the ordering (one median below 1 and the other \
             above): {} of {} — the two orderings' medians are columns, so a fragile row cannot \
             be quoted as a firm one",
            rows.iter()
                .filter(
                    |r| (r.branch_ratio_single_first < 1.0) != (r.branch_ratio_split_first < 1.0)
                )
                .count(),
            rows.len()
        );
        println!(
            "    where the difference goes: the split arm traverses the cell space twice, so it \
             pays x_rows extra innermost-loop exits. Extra mispredictions per x-row exit:"
        );
        for r in &rows {
            println!(
                "      {:<16} {:>4}  x_rows {:>7}  extra/pass {:>+9.0}  per x-row {:>+7.3}",
                r.field,
                r.resolution,
                r.x_rows,
                r.extra_branch_misses(),
                r.extra_branch_misses_per_x_row()
            );
        }
        println!(
            "C3 identical mesh: split == single AND single == MarchingCubes::extract, bit for \
             bit, on every row -> {}",
            if c3_all { "HELD" } else { "FALSIFIED" }
        );
        // The corner-gather control, asserted and then printed. `P-121`'s
        // 15-24% saving comes from *not* gathering eight corner values on a
        // cell that produces nothing; the split arm gathers on every cell and
        // again on every cell it emits, so it cannot be the mechanism behind
        // any number in this file.
        let gather_ratios: Vec<f64> = rows
            .iter()
            .map(|r| r.gather_calls_split() as f64 / r.gather_calls_single() as f64)
            .collect();
        for (r, ratio) in rows.iter().zip(&gather_ratios) {
            assert!(
                *ratio > 1.0,
                "{} {}^3: the split arm gathers {} times against the single arm's {}, which is \
                 not strictly more — then P-121's corner-gather saving is available to it and \
                 this row would measure the wrong mechanism",
                r.field,
                r.resolution,
                r.gather_calls_split(),
                r.gather_calls_single()
            );
        }
        println!(
            "\nCorner gather held constant: one `classify_cell` shared by both arms, called \
             `cells` times on the single arm and `cells + emitting_cells` times on the split \
             arm, so the split gathers MORE on every row (ratio {:.4}..{:.4}). P-121's 15-24% \
             corner-gather saving is therefore not available to it.",
            gather_ratios.iter().copied().fold(f64::MAX, f64::min),
            gather_ratios.iter().copied().fold(0.0f64, f64::max)
        );

        for r in &rows {
            for arm in [Arm::Single, Arm::Split] {
                let counted = r.counted(arm);
                let nanos = r.nanos(arm);
                run.record(&[
                    ("field", r.field.to_string()),
                    ("resolution", r.resolution.to_string()),
                    ("arm", arm.name().to_string()),
                    (
                        "branch_misses_per_cell_single",
                        format!("{:.6}", r.per_cell(Arm::Single, |c| c.branch_misses)),
                    ),
                    (
                        "branch_misses_per_cell_split",
                        format!("{:.6}", r.per_cell(Arm::Split, |c| c.branch_misses)),
                    ),
                    ("branch_miss_ratio", format!("{:.6}", r.branch_miss_ratio())),
                    (
                        "instructions_per_cell_single",
                        format!("{:.4}", r.per_cell(Arm::Single, |c| c.instructions)),
                    ),
                    (
                        "instructions_per_cell_split",
                        format!("{:.4}", r.per_cell(Arm::Split, |c| c.instructions)),
                    ),
                    ("instruction_ratio", format!("{:.6}", r.instruction_ratio())),
                    ("mesh_identical", r.mesh_identical.to_string()),
                    ("c1_holds", r.c1_holds().to_string()),
                    ("c2_holds", r.c2_holds().to_string()),
                    ("c3_holds", r.c3_holds().to_string()),
                    // ── extra columns (M-273) ──
                    ("scalar", "f32".to_string()),
                    ("extractor", "marching_cubes".to_string()),
                    (
                        "branch_misses_per_cell",
                        format!("{:.6}", counted.branch_misses / r.cells as f64),
                    ),
                    (
                        "instructions_per_cell",
                        format!("{:.4}", counted.instructions / r.cells as f64),
                    ),
                    (
                        "cycles_per_cell",
                        format!("{:.4}", counted.cycles / r.cells as f64),
                    ),
                    (
                        "cache_misses_per_cell",
                        format!("{:.6}", counted.cache_misses / r.cells as f64),
                    ),
                    (
                        "l1d_read_misses_per_cell",
                        format!("{:.6}", counted.l1d_read_misses / r.cells as f64),
                    ),
                    (
                        "ipc",
                        format!("{:.4}", counted.instructions / counted.cycles),
                    ),
                    ("ns_per_cell", format!("{:.6}", nanos / r.cells as f64)),
                    ("ghz", format!("{:.4}", counted.cycles / nanos)),
                    ("x_rows", r.x_rows.to_string()),
                    (
                        "extra_branch_misses_per_pass",
                        format!("{:.1}", r.extra_branch_misses()),
                    ),
                    (
                        "extra_branch_misses_per_x_row",
                        format!("{:.4}", r.extra_branch_misses_per_x_row()),
                    ),
                    (
                        "branch_miss_ratio_rep_best",
                        format!("{:.6}", r.branch_ratio_best),
                    ),
                    (
                        "branch_miss_ratio_rep_worst",
                        format!("{:.6}", r.branch_ratio_worst),
                    ),
                    (
                        "instruction_ratio_rep_spread",
                        format!("{:.6}", r.instruction_ratio_spread),
                    ),
                    (
                        "branch_miss_ratio_single_first",
                        format!("{:.6}", r.branch_ratio_single_first),
                    ),
                    (
                        "branch_miss_ratio_split_first",
                        format!("{:.6}", r.branch_ratio_split_first),
                    ),
                    ("cells", r.cells.to_string()),
                    ("emitting_cells", r.emitting_cells.to_string()),
                    (
                        "active_fraction",
                        format!("{:.6}", r.emitting_cells as f64 / r.cells as f64),
                    ),
                    ("vertices", r.vertices.to_string()),
                    ("triangles", r.triangles.to_string()),
                    ("gather_calls_single", r.gather_calls_single().to_string()),
                    ("gather_calls_split", r.gather_calls_split().to_string()),
                    (
                        "gather_calls_ratio",
                        format!(
                            "{:.6}",
                            r.gather_calls_split() as f64 / r.gather_calls_single() as f64
                        ),
                    ),
                    ("control_stream_bytes", r.cells.to_string()),
                    ("payload_length_table_bytes", "256".to_string()),
                    (
                        "mesh_identical_to_shipped",
                        r.mesh_identical_to_shipped.to_string(),
                    ),
                    ("inner_reps", r.inner.to_string()),
                    ("reps", REPS.to_string()),
                    (
                        "worst_instruction_ratio",
                        format!("{worst_instruction_ratio:.6}"),
                    ),
                    (
                        "worst_branch_miss_ratio",
                        format!("{worst_branch_miss_ratio:.6}"),
                    ),
                    (
                        "best_branch_miss_ratio",
                        format!("{best_branch_miss_ratio:.6}"),
                    ),
                    ("rows_reducing_branch_misses", reducing.to_string()),
                    ("rows_total", rows.len().to_string()),
                    ("c1_holds_all_rows", c1_all.to_string()),
                    ("c2_holds_all_rows", c2_all.to_string()),
                    ("c3_holds_all_rows", c3_all.to_string()),
                ]);
            }
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-122");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent, and the registration names it outright: the
    // instrument for C1 *is* `perf_event_open`. Off Linux there is nothing to
    // degrade to — a recorded zero would be a fabricated misprediction rate,
    // and M-281 forbids a clock carrying this verdict.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} measures branch mispredictions per cell, and this platform has no \
             `perf_event_open`. There is no clock substitute and a recorded zero would be a \
             fabrication.",
            prereg.id
        );
        std::process::exit(1);
    }
}
