//! **P-176 — caves connect in three dimensions in a way they cannot in two, and it is a theorem.**
//!
//! Ticket: R-176. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p176
//! ```
//!
//! Writes `docs/experiments/p-176.csv`.
//!
//! # What was missing
//!
//! `isomesh::connectivity::Air` has answered *"is this cave sealed?"* since
//! `R-022a`, and `M-311 / P-23` (`FINDINGS.md:5333`) measured its **repair
//! cost** — 4,872 unions flat across a 59.7× lattice growth, 5.27 per newly-air
//! sample, under the lattice degree of 6. That row's own closing paragraph names
//! the fixture it did not have: *"a field with many small pre-existing cavities
//! would exercise the find path harder without changing the union count, which is
//! a different measurement"* (`FINDINGS.md:5377`). **This is that fixture**, and
//! what it measures is not cost at all: it is the first time this repository has
//! asked what the air region *looks like* — how many components a generated world
//! has, how big the biggest one is, and how that changes with the isovalue.
//!
//! Nothing here has ever been measured. `connectivity.rs` has counted components
//! (`:554`) and sizes (`:584`) since Phase 8 and nothing has ever read the
//! *distribution*; `grep percolat FINDINGS.md` returns one hit and it is the
//! Phase 27 probe table recording the source as newly acquired
//! (`FINDINGS.md:25061`).
//!
//! Duminil-Copin, Rivera, Rodriguez & Vanneuville (`arXiv:2108.08008` /
//! `10.1214/22-aop1594`) prove that for smooth Gaussian fields in dimension
//! `d >= 3` the critical level `l_c(d)` is **strictly positive**: the excursion
//! set percolates at levels where the two-dimensional analogue — whose critical
//! level is zero — does not. In this crate's units the excursion set at level `l`
//! is exactly `Air`'s air region of the field shifted by `l`, because
//! `Air::build` takes air to be `value >= 0` (`connectivity.rs:218-219`), so
//! `Air::build(values - l)` **is** `{f >= l}` with no reinterpretation and no new
//! code inside `src/`.
//!
//! # The instruments, and why there are two
//!
//! **Instrument A, bench-local.** A three-forward-edge union-find with union by
//! size and path halving: scan in index order and, for each air sample, union it
//! with its `-x`, `-y`, `-z` neighbours where those are air, so each lattice edge
//! is seen once. This is deliberately *the retired algorithm* — M-311's own text
//! records that the old union-find build *"visits three forward edges per air
//! sample so each lattice edge is seen once"* (`FINDINGS.md:5364`) — because
//! `Air` is no longer a union-find at all. `✗26` replaced the parent pointers
//! with a **flat label array filled by breadth-first flood**
//! (`connectivity.rs:29-46`, `flood` at `:708`). The two instruments here
//! therefore share no data structure, no traversal order and no merge rule. They
//! share the six-neighbour adjacency and the `value >= 0` test, which is exactly
//! what C3 is entitled to assume.
//!
//! **Instrument B, the crate's.** `Air::build`, read through `components()`
//! (`:554`), `air_samples()` (`:560`), `label_count()` (`:577`) and
//! `component_size(l)` (`:584`). C3 compares the component count, the air-sample
//! count **and the full sorted multiset of component sizes** — not merely the
//! count, because two labellings can agree on how many components there are and
//! disagree about which samples are in them.
//!
//! Both are asked the identical question. The air test is written
//! `values[i] - iso >= 0.0` in instrument A, and instrument B is handed an array
//! whose entry is that same subtraction already performed, so the two cannot
//! disagree through rounding: they compare the same `f64`.
//!
//! **Connectivity is 6 in 3D and 4 in 2D**, and that is not this bench's choice.
//! `Air::neighbours` (`connectivity.rs:610-646`) is the six axis-aligned
//! neighbours, so instrument A must use six or C3 would be comparing two
//! different graphs; the 2D slice takes the same rule restricted to a plane,
//! which is 4. Both are written out as columns (`connectivity_3d`,
//! `connectivity_2d`) so no reader has to take this paragraph's word for it.
//!
//! # The measured region, and the control that makes it necessary
//!
//! `fbm_terrain` is swept over **its own domain**, `[-8, 8]^3`
//! (`fields/mod.rs:1381`), unmasked.
//!
//! `noise_cavity` is `NoiseVolume ∩ Sphere{r: 1.5}` over `[-2, 2]^3`
//! (`fields/mod.rs:1237-1245`), and over the **whole** domain the question stops
//! being about the field: outside the cap every sample has `sphere > 0`, hence
//! `max(noise, sphere) > 0`, hence air at every isovalue at or below zero — and
//! that shell is one connected object wrapping the solid. **Measured, as this
//! bench's first act**: at `iso = 0` the unmasked domain is `82.58%` air in 12
//! components whose largest holds `99.99%` of it. C1 would then be true by the
//! sampling box rather than by the field, which is `M-44`'s failure mode exactly.
//! Those three numbers go on every row (`outer_shell_*`) and are asserted, so the
//! mask below is load-bearing rather than decorative.
//!
//! So `noise_cavity`'s region is the field's **own** cap:
//! `{p : sphere_sdf(p) < -0.05}`, with the sweep floored at `iso = -0.05`. On
//! that region over that range the cap cannot contribute — `sphere_sdf` is below
//! `-0.05 <= iso` everywhere in the mask, so `{max(noise, sphere) >= iso}` is
//! exactly `{noise >= iso}`. `cap_max_sdf` records the largest `sphere_sdf`
//! inside the mask and a vacuity control asserts it is **strictly** below the
//! sweep floor. Masked-out samples are handed to instrument B as `-1.0`, i.e. as
//! solid, so both instruments see the identical region.
//!
//! # The ladder
//!
//! 41 isovalues spaced evenly across each field's **own sampled value range**,
//! swept from the top down, plus `iso = 0.0` inserted when the range straddles it
//! — which it does for both fields, giving 42 rungs each. Spanning the sampled
//! range is what makes the vacuity control reachable by construction: the top
//! rung admits exactly the argmax sample and the bottom rung admits everything
//! the mask allows.
//!
//! `iso = 0` is inserted rather than left to chance because the registration's
//! own sentence is about it: *"at isovalue 0 the excavated region should have one
//! giant component rather than many isolated pockets"*. `giant_at_zero` and
//! `largest_component_fraction_at_zero` answer that sentence directly.
//!
//! Resolution is a single `97^3` — 912,673 samples — and one resolution rather
//! than three because percolation is a property of the field's *feature count*,
//! not of the sampling rate: refining the grid resolves the same blobs better and
//! does not add blobs. 97 puts about 7 samples across `noise_cavity`'s
//! `1 / 3.45 ≈ 0.29` feature (`fields/mod.rs:1152`) and about 10 features across
//! the cap, and the whole bench runs in about 1.3 s.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `fbm_terrain` × 42 isovalues | the level, over `[-8, 8]^3` unmasked | no |
//! | `noise_cavity` × 42 isovalues | the level, over the field's own cap | no |
//! | the 2D slice, `z = 48`, 4-connectivity, every row | the **dimension**, at a matched level | **yes** — C2's control |
//! | the unmasked `noise_cavity` at `iso = 0` | the mask | **yes** — `outer_shell_*` |
//! | instrument B (`Air::build`) on every row | the algorithm | **yes** — C3's control |
//!
//! # The three clauses, and how each is decided
//!
//! **C1** — *a giant component appears, and the isovalue at which it appears is
//! reported.* `giant_component_exists` is the registered definition applied to
//! the row: one component holding **above 50%** of the air volume. That is also
//! true at the very top of the ladder, where the air is a single voxel, so
//! `percolation_isovalue` is the **persistent onset**: the highest swept isovalue
//! at or below which the giant component never disappears again. It is then
//! refined by 14 bisections inside the bracketing pair of rungs, and the bracket
//! is recorded (`percolation_rung_giant`, `percolation_rung_fragmented`) so the
//! refinement is checkable. `c1_holds` is per field — an onset exists **and** the
//! sweep also contains a fragmented rung, because an onset with nothing above it
//! is not a transition.
//!
//! **C2** — *the 3D behaviour differs qualitatively from a 2D slice.* Not "the
//! numbers differ", since every pair of numbers differs. The criterion is the
//! theorem's own shape: `c2_holds` iff the 3D set has a persistent
//! giant-component phase **and the 2D slice has none anywhere in the swept
//! range**. A sign-carrying reading comes free — an onset above `iso = 0` in 3D
//! and no onset at all down to a negative floor in 2D is `l_c(3) > 0 > l_c(2)`
//! measured in the crate's own field units.
//!
//! **C3** — per row: instrument B's component count, air-sample count and sorted
//! size multiset all equal instrument A's, and `label_count() == components()` so
//! no retired label is hiding a component. `air_uf_*` carry B's own readings, so
//! a disagreement is diagnosable from the CSV rather than from a rerun.
//!
//! `c1_holds` and `c2_holds` are properties of a whole sweep, so a field's
//! verdict is written on every one of that field's rows. `c3_holds` and
//! `air_union_find_agreement` are genuinely per row.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration says **`SHARE: none`** — *"this predicts a property of
//! generated worlds, and `M-311`'s union-find is the instrument"*. Recomputed and
//! unchanged: nothing here proposes a source change or moves a stage's time. It
//! changes what a level designer may assume, not what the extractor does. The
//! `cave_percolation` demo consumes these rows, which is why the
//! giant-component fraction is recorded at **every** swept isovalue and not only
//! at the transition.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record`, and every panic starts `VOID: `.
//!
//! - **Both regimes are inside the swept range**, per field, which is the
//!   registration's own control verbatim. `fragmented_rows` counts rungs with at
//!   least 8 components and no giant one; `giant_rows_with_real_air` counts rungs
//!   with a giant component over at least 5% air. The air-share floor is what
//!   stops the ladder's one-voxel top rung from answering for the whole
//!   one-large-component regime.
//! - **The 2D control is populated.** `slice_populated_rows` counts rungs whose
//!   slice holds at least 100 air pixels; without one, C2 compares a 3D census
//!   against an empty plane and its verdict is a default.
//! - **The cap cannot make air.** `cap_max_sdf < sweep_lo` for `noise_cavity`, or
//!   the components would be measuring the sphere.
//! - **The mask is necessary.** The unmasked `noise_cavity` at `iso = 0` must
//!   itself show a giant component over a more-than-half-air domain, or the
//!   masking above is a complication with nothing behind it.
//! - **Every sampled value is finite**, because `fbm_terrain` declares
//!   `FieldBound::Unbounded` (`fields/mod.rs:1393`) and a NaN would make
//!   `v - iso >= 0.0` false everywhere and read as a cleanly empty air set.
//!
//! # Determinism
//!
//! No RNG: both fields are deterministic hash noise, the ladder is arithmetic on
//! the sampled range, and every sort is on integer keys or `f64::total_cmp`.
//! `row_ms` sits beside the verdicts and gates nothing.

mod common;

use std::time::Instant;

use isomesh::connectivity::Air;
use isomesh::fields::{FbmTerrain, ReferenceField, Sphere, noise_cavity};
use isomesh::{RuntimeShape3, Sdf};

/// Samples per axis. One resolution; see the header for why.
const SAMPLES: u32 = 97;
/// Evenly spaced rungs across a field's sampled value range.
const RUNGS: usize = 41;
/// How far inside its own cap `noise_cavity` is measured, and the sweep floor.
const CAP_MARGIN: f64 = 0.05;
/// The cap's radius, from `fields/mod.rs:1242`.
const CAP_RADIUS: f64 = 1.5;
/// The registered giant-component threshold: *above* half the air volume.
const GIANT_SHARE: f64 = 0.5;
/// "Many" small components, for the fragmented-regime control.
const MANY_COMPONENTS: u64 = 8;
/// Air share of the region below which a giant component is a ladder artefact.
const REAL_AIR_SHARE: f64 = 0.05;
/// Air pixels a slice needs before its verdict is evidence about dimension.
const SLICE_FLOOR: u64 = 100;
/// Halvings used to refine a persistent onset inside its bracketing rungs.
const BISECTIONS: u32 = 14;
/// Component sizes written out individually in `component_size_distribution`.
const TOP_SIZES: usize = 5;
/// How close to zero a rung must be for `iso = 0` to count as already present.
const ZERO_RUNG: f64 = 1e-12;

/// One field, sampled once, with the region its sweep is confined to.
struct Region {
    /// The `ReferenceField::NAME` this field's rows report as.
    name: &'static str,
    /// The sampled values, `x` fastest.
    values: Vec<f64>,
    /// Which samples the sweep may call air at all.
    mask: Vec<bool>,
    /// Samples per axis, all three equal.
    dims: [u32; 3],
    /// The shape instrument B is handed.
    shape: RuntimeShape3,
    /// `mask`'s count of `true`: the denominator of every air share.
    region: u64,
    /// Inclusive index extent of `mask` per axis, for the span reading.
    extent: [(u32, u32); 3],
    /// Human-readable statement of what `mask` is.
    mask_rule: &'static str,
    /// Largest cap value inside the mask, or `-inf` where there is no cap.
    cap_max: f64,
    /// The sweep's floor.
    sweep_lo: f64,
    /// The sweep's ceiling.
    sweep_hi: f64,
}

/// One connected-component census.
struct Comp {
    /// Air samples found.
    air: u64,
    /// Components found.
    components: u64,
    /// Component sizes, descending.
    sizes: Vec<u32>,
    /// `sizes[0]`, or zero when there is no air.
    largest: u32,
    /// Slot index of the largest component's root, or `u32::MAX` when none.
    root: u32,
    /// Extent of the largest component over the region's extent, worst axis.
    span: f64,
}

impl Comp {
    /// The registered `largest_component_fraction`.
    fn fraction(&self) -> f64 {
        if self.air == 0 {
            0.0
        } else {
            f64::from(self.largest) / self.air as f64
        }
    }

    /// The registered giant-component test: *above* half the air volume.
    fn giant(&self) -> bool {
        self.fraction() > GIANT_SHARE
    }

    /// `"66c|26695|2885|854|611|470|+61"` — the count, the biggest few, the rest.
    fn distribution(&self) -> String {
        let mut parts = vec![format!("{}c", self.components)];
        parts.extend(self.sizes.iter().take(TOP_SIZES).map(u32::to_string));
        if self.sizes.len() > TOP_SIZES {
            parts.push(format!("+{}", self.sizes.len() - TOP_SIZES));
        }
        parts.join("|")
    }
}

/// One rung of the ladder: both dimensions and both instruments.
struct Row {
    /// The isovalue.
    iso: f64,
    /// Instrument A in 3D, 6-connectivity.
    three: Comp,
    /// Instrument A on the mid-`z` plane, 4-connectivity.
    two: Comp,
    /// Instrument B's component count.
    air_components: u64,
    /// Instrument B's air-sample count.
    air_air: u64,
    /// Instrument B's label-array length.
    air_labels: u64,
    /// Instrument B's sorted size multiset equals instrument A's.
    air_sizes_match: bool,
    /// Wall clock for the rung, gating nothing.
    ms: f64,
}

impl Row {
    /// C3 for this row: same count, same air, same sizes, no retired label.
    fn agrees(&self) -> bool {
        self.air_components == self.three.components
            && self.air_air == self.three.air
            && self.air_labels == self.air_components
            && self.air_sizes_match
    }
}

/// A union-find over the whole sample array; only air samples are made live.
///
/// Union by size with path halving. Reused across rungs — `reset` is two
/// `fill`s, so a 42-rung sweep allocates twice.
struct Uf {
    /// `u32::MAX` where the sample is not air; otherwise the parent index.
    parent: Vec<u32>,
    /// Set size, meaningful only at a root.
    size: Vec<u32>,
}

impl Uf {
    /// One slot per sample, all dead.
    fn with(n: usize) -> Self {
        Self {
            parent: vec![u32::MAX; n],
            size: vec![0; n],
        }
    }

    /// Back to all dead.
    fn reset(&mut self) {
        self.parent.fill(u32::MAX);
        self.size.fill(0);
    }

    /// Whether slot `i` has been made air on this rung.
    fn live(&self, i: usize) -> bool {
        self.parent[i] != u32::MAX
    }

    /// Start a singleton at `i`.
    fn make(&mut self, i: usize) {
        self.parent[i] = i as u32;
        self.size[i] = 1;
    }

    /// Root of `i`, halving the path on the way.
    fn find(&mut self, mut i: u32) -> u32 {
        while self.parent[i as usize] != i {
            let grand = self.parent[self.parent[i as usize] as usize];
            self.parent[i as usize] = grand;
            i = grand;
        }
        i
    }

    /// Join the sets of `a` and `b`, the larger set winning.
    fn union(&mut self, a: usize, b: usize) {
        let (mut ra, mut rb) = (self.find(a as u32), self.find(b as u32));
        if ra == rb {
            return;
        }
        if self.size[ra as usize] < self.size[rb as usize] {
            core::mem::swap(&mut ra, &mut rb);
        }
        self.parent[rb as usize] = ra;
        self.size[ra as usize] += self.size[rb as usize];
    }

    /// Roots and their sizes over the first `n` slots, descending.
    fn harvest(&self, n: usize) -> Comp {
        let mut sizes = Vec::new();
        let mut air = 0u64;
        let mut largest = 0u32;
        let mut root = u32::MAX;
        for (i, &p) in self.parent.iter().take(n).enumerate() {
            if p != u32::MAX {
                air += 1;
            }
            if p == i as u32 {
                let s = self.size[i];
                sizes.push(s);
                if s > largest {
                    largest = s;
                    root = i as u32;
                }
            }
        }
        sizes.sort_unstable_by(|a, b| b.cmp(a));
        Comp {
            air,
            components: sizes.len() as u64,
            sizes,
            largest,
            root,
            span: 0.0,
        }
    }
}

/// `field.sample` over `n^3` points from `lo` at spacing `h`, `x` fastest.
fn sample_grid<S: Sdf<Scalar = f64>>(field: &S, lo: [f64; 3], h: f64, n: u32) -> Vec<f64> {
    let mut out = Vec::with_capacity((n as usize).pow(3));
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                out.push(field.sample([
                    lo[0] + h * f64::from(i),
                    lo[1] + h * f64::from(j),
                    lo[2] + h * f64::from(k),
                ]));
            }
        }
    }
    out
}

/// Instrument A in 3D: 6-connectivity over `{mask && value - iso >= 0}`.
fn census_3d(f: &Region, iso: f64, uf: &mut Uf) -> Comp {
    let [nx, ny, nz] = f.dims;
    let sy = nx as usize;
    let sz = (nx as usize) * (ny as usize);
    uf.reset();
    for k in 0..nz as usize {
        for j in 0..ny as usize {
            for i in 0..nx as usize {
                let idx = i + j * sy + k * sz;
                let is_air = f.mask[idx] && f.values[idx] - iso >= 0.0;
                if !is_air {
                    continue;
                }
                uf.make(idx);
                if i > 0 && uf.live(idx - 1) {
                    uf.union(idx, idx - 1);
                }
                if j > 0 && uf.live(idx - sy) {
                    uf.union(idx, idx - sy);
                }
                if k > 0 && uf.live(idx - sz) {
                    uf.union(idx, idx - sz);
                }
            }
        }
    }
    let mut comp = uf.harvest(f.values.len());
    comp.span = span_of_component(f, uf, comp.root);
    comp
}

/// Index extent of the component rooted at `root` over the region's own extent,
/// worst axis. `1.0` means it reaches as far as the region does on every axis.
fn span_of_component(f: &Region, uf: &mut Uf, root: u32) -> f64 {
    if root == u32::MAX {
        return 0.0;
    }
    let [nx, ny, nz] = f.dims;
    let sy = nx as usize;
    let sz = (nx as usize) * (ny as usize);
    let mut lo = [u32::MAX; 3];
    let mut hi = [0u32; 3];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let idx = i as usize + j as usize * sy + k as usize * sz;
                if !uf.live(idx) || uf.find(idx as u32) != root {
                    continue;
                }
                for (axis, v) in [i, j, k].into_iter().enumerate() {
                    lo[axis] = lo[axis].min(v);
                    hi[axis] = hi[axis].max(v);
                }
            }
        }
    }
    (0..3)
        .map(|axis| {
            let (rlo, rhi) = f.extent[axis];
            let whole = f64::from(rhi - rlo) + 1.0;
            let mine = if lo[axis] == u32::MAX {
                0.0
            } else {
                f64::from(hi[axis] - lo[axis]) + 1.0
            };
            mine / whole
        })
        .fold(f64::INFINITY, f64::min)
}

/// Instrument A on the middle `z` plane: 4-connectivity, same air test.
fn census_2d(f: &Region, iso: f64, uf: &mut Uf) -> Comp {
    let [nx, ny, nz] = f.dims;
    let plane = (nz / 2) as usize * (nx as usize) * (ny as usize);
    uf.reset();
    for j in 0..ny as usize {
        for i in 0..nx as usize {
            let local = i + j * nx as usize;
            let idx = plane + local;
            let is_air = f.mask[idx] && f.values[idx] - iso >= 0.0;
            if !is_air {
                continue;
            }
            uf.make(local);
            if i > 0 && uf.live(local - 1) {
                uf.union(local, local - 1);
            }
            if j > 0 && uf.live(local - nx as usize) {
                uf.union(local, local - nx as usize);
            }
        }
    }
    uf.harvest((nx as usize) * (ny as usize))
}

/// Instrument B: `Air::build` on the shifted, masked array.
///
/// Returns `(components, air_samples, label_count, sizes descending)`.
fn census_air(f: &Region, iso: f64, shifted: &mut [f64]) -> (u64, u64, u64, Vec<u32>) {
    for (slot, (v, m)) in shifted.iter_mut().zip(f.values.iter().zip(f.mask.iter())) {
        *slot = if *m { *v - iso } else { -1.0 };
    }
    let (built, _repair) =
        Air::build(shifted, &f.shape).expect("P-176: Air::build over the bench grid");
    let mut sizes: Vec<u32> = (0..built.label_count() as u32)
        .map(|l| built.component_size(l))
        .filter(|s| *s > 0)
        .collect();
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    (
        built.components(),
        built.air_samples(),
        built.label_count() as u64,
        sizes,
    )
}

/// 41 rungs across `[lo, hi]`, descending, with `iso = 0` inserted when the
/// range straddles it and no rung already lands there.
fn ladder(lo: f64, hi: f64) -> Vec<f64> {
    let span = hi - lo;
    let last = (RUNGS - 1) as f64;
    let mut rungs: Vec<f64> = (0..RUNGS).map(|k| hi - span * (k as f64) / last).collect();
    if lo < 0.0 && hi > 0.0 && !rungs.iter().any(|v| v.abs() < ZERO_RUNG) {
        rungs.push(0.0);
    }
    rungs.sort_by(|a, b| b.total_cmp(a));
    rungs
}

/// The persistent onset: the lowest index from which `giant` holds on every
/// remaining rung. `None` when the bottom rung is not giant.
fn persistent_onset(rows: &[Row], giant: impl Fn(&Row) -> bool) -> Option<usize> {
    let mut onset = None;
    for i in (0..rows.len()).rev() {
        if !giant(&rows[i]) {
            break;
        }
        onset = Some(i);
    }
    onset
}

/// Refine an onset inside its bracketing rungs: `giant_at` carries a giant
/// component, `plain_at` does not, and the return is the highest isovalue shown
/// to carry one.
fn refine(f: &Region, giant_at: f64, plain_at: f64, uf: &mut Uf, slice: bool) -> f64 {
    let (mut lo, mut hi) = (giant_at, plain_at);
    for _ in 0..BISECTIONS {
        let mid = 0.5 * (lo + hi);
        let comp = if slice {
            census_2d(f, mid, uf)
        } else {
            census_3d(f, mid, uf)
        };
        if comp.giant() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

/// An isovalue as a CSV value, or `none` where there is not one.
fn onset_value(iso: Option<f64>) -> String {
    match iso {
        Some(v) => format!("{v:.6}"),
        None => String::from("none"),
    }
}

/// Inclusive index extent of a mask per axis.
fn mask_extent(mask: &[bool], dims: [u32; 3]) -> [(u32, u32); 3] {
    let [nx, ny, nz] = dims;
    let sy = nx as usize;
    let sz = (nx as usize) * (ny as usize);
    let mut lo = [u32::MAX; 3];
    let mut hi = [0u32; 3];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                if !mask[i as usize + j as usize * sy + k as usize * sz] {
                    continue;
                }
                for (axis, v) in [i, j, k].into_iter().enumerate() {
                    lo[axis] = lo[axis].min(v);
                    hi[axis] = hi[axis].max(v);
                }
            }
        }
    }
    [
        (lo[0].min(hi[0]), hi[0]),
        (lo[1].min(hi[1]), hi[1]),
        (lo[2].min(hi[2]), hi[2]),
    ]
}

/// The unmasked-domain control: `noise_cavity` over all of `[-2, 2]^3` at
/// `iso = 0`, which is what C1 would have measured without the cap mask.
struct Shell {
    /// Air share of the whole sampling box.
    air_fraction: f64,
    /// Components in it.
    components: u64,
    /// Share of that air held by the largest of them.
    largest_fraction: f64,
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-176");

    let n = SAMPLES;
    let count = (n as usize).pow(3);
    let shape = RuntimeShape3::new([n; 3]).expect("P-176: 97^3 fits u32");

    // ── fbm_terrain, its own domain, unmasked ──────────────────────────────
    let terrain_field = FbmTerrain::<f64>::canonical();
    let (tlo, thi) = terrain_field.domain();
    let th = (thi[0] - tlo[0]) / f64::from(n - 1);
    let tvalues = sample_grid(&terrain_field, tlo, th, n);
    let tmin = tvalues.iter().copied().fold(f64::INFINITY, f64::min);
    let tmax = tvalues.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let terrain = Region {
        name: "fbm_terrain",
        values: tvalues,
        mask: vec![true; count],
        dims: [n; 3],
        shape,
        region: count as u64,
        extent: [(0, n - 1); 3],
        mask_rule: "full_domain",
        cap_max: f64::NEG_INFINITY,
        sweep_lo: tmin,
        sweep_hi: tmax,
    };

    // ── noise_cavity, masked to the inside of its own cap ──────────────────
    let cavity_field = noise_cavity::<f64>();
    let (clo, chi) = cavity_field.domain();
    let ch = (chi[0] - clo[0]) / f64::from(n - 1);
    let cvalues = sample_grid(&cavity_field, clo, ch, n);
    let cap = Sphere::<f64> {
        center: [0.0; 3],
        radius: CAP_RADIUS,
    };
    let capvalues = sample_grid(&cap, clo, ch, n);
    let cmask: Vec<bool> = capvalues.iter().map(|c| *c < -CAP_MARGIN).collect();
    let cap_max = capvalues
        .iter()
        .zip(&cmask)
        .filter(|(_, m)| **m)
        .map(|(c, _)| *c)
        .fold(f64::NEG_INFINITY, f64::max);
    let cmax = cvalues
        .iter()
        .zip(&cmask)
        .filter(|(_, m)| **m)
        .map(|(v, _)| *v)
        .fold(f64::NEG_INFINITY, f64::max);
    let cregion = cmask.iter().filter(|m| **m).count() as u64;
    let cextent = mask_extent(&cmask, [n; 3]);
    let cavity = Region {
        name: "noise_cavity",
        values: cvalues,
        mask: cmask,
        dims: [n; 3],
        shape,
        region: cregion,
        extent: cextent,
        mask_rule: "cap_sdf_lt_-0.05",
        cap_max,
        sweep_lo: -CAP_MARGIN,
        sweep_hi: cmax,
    };

    let mut uf = Uf::with(count);
    let mut shifted = vec![0.0f64; count];

    // ── the unmasked control, measured before anything is recorded ─────────
    let shell = {
        let unmasked = Region {
            name: "noise_cavity_unmasked",
            values: cavity.values.clone(),
            mask: vec![true; count],
            dims: [n; 3],
            shape,
            region: count as u64,
            extent: [(0, n - 1); 3],
            mask_rule: "full_domain",
            cap_max: f64::NEG_INFINITY,
            sweep_lo: cavity.sweep_lo,
            sweep_hi: cavity.sweep_hi,
        };
        let comp = census_3d(&unmasked, 0.0, &mut uf);
        Shell {
            air_fraction: comp.air as f64 / count as f64,
            components: comp.components,
            largest_fraction: comp.fraction(),
        }
    };

    common::experiment::run(prereg, |run| {
        // ── vacuity controls, all of them, before any row ──────────────────
        for f in [&terrain, &cavity] {
            assert!(
                f.values.iter().all(|v| v.is_finite()),
                "VOID: {} sampled a non-finite value, and `v - iso >= 0.0` is false \
                 for a NaN, so the air set would read as cleanly empty instead of \
                 as broken",
                f.name
            );
            assert!(
                f.region > 0,
                "VOID: {}'s measured region is empty, so every component count \
                 below is zero by construction",
                f.name
            );
            assert!(
                f.sweep_hi > f.sweep_lo,
                "VOID: {}'s sweep range is degenerate ({} to {}), so the ladder is \
                 one isovalue repeated {} times",
                f.name,
                f.sweep_lo,
                f.sweep_hi,
                RUNGS
            );
        }
        assert!(
            cavity.cap_max < cavity.sweep_lo,
            "VOID: noise_cavity's cap reaches {} inside the mask while the sweep \
             floor is {}, so at the bottom rungs the cap alone makes samples air \
             and the components would be measuring the sphere rather than the cave",
            cavity.cap_max,
            cavity.sweep_lo
        );
        assert!(
            shell.largest_fraction > GIANT_SHARE && shell.air_fraction > GIANT_SHARE,
            "VOID: the unmasked noise_cavity at iso 0 is {:.4} air in {} components \
             whose largest holds {:.4} — the cap mask exists to remove a trivially \
             giant outer shell, and with no such shell the mask is a complication \
             with nothing behind it",
            shell.air_fraction,
            shell.components,
            shell.largest_fraction
        );

        for f in [&terrain, &cavity] {
            let rungs = ladder(f.sweep_lo, f.sweep_hi);
            let step = (f.sweep_hi - f.sweep_lo) / (RUNGS - 1) as f64;
            let mut rows: Vec<Row> = Vec::with_capacity(rungs.len());
            for iso in rungs {
                let started = Instant::now();
                let three = census_3d(f, iso, &mut uf);
                let (air_components, air_air, air_labels, air_sizes) =
                    census_air(f, iso, &mut shifted);
                let air_sizes_match = air_sizes == three.sizes;
                let two = census_2d(f, iso, &mut uf);
                rows.push(Row {
                    iso,
                    three,
                    two,
                    air_components,
                    air_air,
                    air_labels,
                    air_sizes_match,
                    ms: started.elapsed().as_secs_f64() * 1e3,
                });
            }

            // ── the sweep's shape, in five counts ──────────────────────────
            let fragmented = rows
                .iter()
                .filter(|r| r.three.components >= MANY_COMPONENTS && !r.three.giant())
                .count() as u64;
            let real_giant = rows
                .iter()
                .filter(|r| {
                    r.three.giant() && r.three.air as f64 / f.region as f64 >= REAL_AIR_SHARE
                })
                .count() as u64;
            let single = rows
                .iter()
                .filter(|r| r.three.components == 1 && r.three.air > 0)
                .count() as u64;
            let slice_populated = rows.iter().filter(|r| r.two.air >= SLICE_FLOOR).count() as u64;
            let dimension_gap = rows
                .iter()
                .filter(|r| r.three.giant() && !r.two.giant() && r.two.air >= SLICE_FLOOR)
                .count() as u64;

            assert!(
                fragmented > 0,
                "VOID: {}: not one of the {} swept isovalues produced at least {} \
                 components without a giant one, so the many-small-components \
                 regime is outside the swept range and the transition is not \
                 bracketed",
                f.name,
                rows.len(),
                MANY_COMPONENTS
            );
            assert!(
                real_giant > 0,
                "VOID: {}: not one swept isovalue produced a giant component over \
                 at least {:.0}% air, so the one-large-component regime is outside \
                 the swept range and any giant reading is the ladder's top rung \
                 answering with a single voxel",
                f.name,
                REAL_AIR_SHARE * 100.0
            );
            assert!(
                slice_populated > 0,
                "VOID: {}: no swept isovalue gave the 2D slice {} air pixels, so C2 \
                 would be comparing a 3D census against an empty plane",
                f.name,
                SLICE_FLOOR
            );

            // ── the onsets, and the brackets they were refined inside ──────
            let onset3 = persistent_onset(&rows, |r| r.three.giant());
            let onset2 = persistent_onset(&rows, |r| r.two.giant());
            let bracket3 = onset3
                .filter(|i| *i > 0)
                .map(|i| (rows[i].iso, rows[i - 1].iso));
            let bracket2 = onset2
                .filter(|i| *i > 0)
                .map(|i| (rows[i].iso, rows[i - 1].iso));
            let refined3 = match bracket3 {
                Some((giant_at, plain_at)) => Some(refine(f, giant_at, plain_at, &mut uf, false)),
                None => onset3.map(|i| rows[i].iso),
            };
            let refined2 = match bracket2 {
                Some((giant_at, plain_at)) => Some(refine(f, giant_at, plain_at, &mut uf, true)),
                None => onset2.map(|i| rows[i].iso),
            };

            let c1 = refined3.is_some() && fragmented > 0;
            let c2 = refined3.is_some() && refined2.is_none();

            let at_zero = rows.iter().find(|r| r.iso.abs() < ZERO_RUNG);
            let giant_at_zero = at_zero.is_some_and(|r| r.three.giant());
            let fraction_at_zero = at_zero.map_or(0.0, |r| r.three.fraction());

            for r in &rows {
                run.record(&[
                    ("field", f.name.to_string()),
                    ("isovalue", format!("{:.6}", r.iso)),
                    ("components", r.three.components.to_string()),
                    (
                        "largest_component_fraction",
                        format!("{:.6}", r.three.fraction()),
                    ),
                    ("component_size_distribution", r.three.distribution()),
                    ("giant_component_exists", r.three.giant().to_string()),
                    ("percolation_isovalue", onset_value(refined3)),
                    (
                        "two_d_slice_comparison",
                        format!(
                            "3d:{}@{:.4}|2d:{}@{:.4}",
                            r.three.components,
                            r.three.fraction(),
                            r.two.components,
                            r.two.fraction()
                        ),
                    ),
                    ("air_union_find_agreement", r.agrees().to_string()),
                    ("c1_holds", c1.to_string()),
                    ("c2_holds", c2.to_string()),
                    ("c3_holds", r.agrees().to_string()),
                    // ── extras (M-273) ──
                    //
                    // the grid, the region, and the rule that defined it
                    ("resolution", SAMPLES.to_string()),
                    ("connectivity_3d", "6".to_string()),
                    ("connectivity_2d", "4".to_string()),
                    ("slice_plane", format!("z={}", SAMPLES / 2)),
                    ("mask_rule", f.mask_rule.to_string()),
                    ("region_samples", f.region.to_string()),
                    ("cap_max_sdf", format!("{:.6}", f.cap_max)),
                    // this rung's 3D census, in full
                    ("air_samples", r.three.air.to_string()),
                    (
                        "air_fraction_of_region",
                        format!("{:.6}", r.three.air as f64 / f.region as f64),
                    ),
                    ("largest_component_size", r.three.largest.to_string()),
                    (
                        "largest_component_span_fraction",
                        format!("{:.6}", r.three.span),
                    ),
                    // C2's control arm, on the same rung
                    ("components_2d", r.two.components.to_string()),
                    ("air_samples_2d", r.two.air.to_string()),
                    (
                        "largest_component_fraction_2d",
                        format!("{:.6}", r.two.fraction()),
                    ),
                    ("giant_component_2d", r.two.giant().to_string()),
                    ("percolation_isovalue_2d", onset_value(refined2)),
                    // the onset, and the two rungs it was refined between
                    (
                        "percolation_rung_giant",
                        onset_value(onset3.map(|i| rows[i].iso)),
                    ),
                    (
                        "percolation_rung_fragmented",
                        onset_value(bracket3.map(|(_, plain_at)| plain_at)),
                    ),
                    // the ladder
                    ("sweep_lo", format!("{:.6}", f.sweep_lo)),
                    ("sweep_hi", format!("{:.6}", f.sweep_hi)),
                    ("sweep_rungs", rows.len().to_string()),
                    ("sweep_step", format!("{step:.6}")),
                    // the five counts, so a global verdict is checkable from one row
                    ("fragmented_rows", fragmented.to_string()),
                    ("giant_rows_with_real_air", real_giant.to_string()),
                    ("single_component_rows", single.to_string()),
                    ("slice_populated_rows", slice_populated.to_string()),
                    ("dimension_gap_rows", dimension_gap.to_string()),
                    // the registration's own sentence, answered directly
                    ("giant_at_zero", giant_at_zero.to_string()),
                    (
                        "largest_component_fraction_at_zero",
                        format!("{fraction_at_zero:.6}"),
                    ),
                    // instrument B's own readings, so a false C3 is diagnosable
                    ("air_uf_components", r.air_components.to_string()),
                    ("air_uf_air_samples", r.air_air.to_string()),
                    ("air_uf_labels", r.air_labels.to_string()),
                    ("air_uf_sizes_match", r.air_sizes_match.to_string()),
                    // the unmasked control, identical on every row
                    (
                        "outer_shell_air_fraction",
                        format!("{:.6}", shell.air_fraction),
                    ),
                    ("outer_shell_components", shell.components.to_string()),
                    (
                        "outer_shell_largest_fraction",
                        format!("{:.6}", shell.largest_fraction),
                    ),
                    // time, beside the verdict, gating nothing
                    ("row_ms", format!("{:.1}", r.ms)),
                ]);
            }
        }
    });
}
