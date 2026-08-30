//! **P-163 — a null registered on purpose: FCC against BCC is `0.011 dB` and
//! should be under the harness's own scatter.**
//!
//! Ticket: R-163. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p163
//! ```
//!
//! Writes `docs/experiments/p-163.csv`.
//!
//! # What was missing
//!
//! **P-162 measured BCC against the cubic grid and both of its measurable
//! clauses failed.** `docs/experiments/p-162.csv` is committed and says so in
//! numbers: BCC improved symmetric Hausdorff on **3 of 8** reference fields
//! against a bar of five, so C1 was FALSIFIED; the measured gain did not land
//! within a factor of two of the predicted `0.2571 dB` on five fields, so C2 was
//! FALSIFIED; C3 held, the BCC case table being 16 entries against the cubic
//! 256. The per-field gains it recorded run from `−8.585 dB` (`thin_plate`) to
//! `+0.727 dB` (`gyroid`) — **nine decibels of spread around a prediction of a
//! quarter of one**. That is the baseline this row is built on and it is quoted
//! from the CSV, not from a summary.
//!
//! Against that spread, `G(D₃) = 2^(−11/3) = 0.078745066` versus
//! `G(A₃*) = 19/(192·∛2) = 0.078543281` is a gap of
//!
//! ```text
//!   10·log₁₀(G(D₃) / G(A₃*)) = 0.011143 dB
//! ```
//!
//! — 4.3% of the cubic gap P-162 could not resolve. `Lattice::g` asserts both
//! decimals against their closed forms on every call
//! (`benches/common/lattice.rs:234-252`, residuals `2.17e-10` and `3.82e-10`),
//! so the prediction cannot rot into a transposed digit. What was missing is not
//! the prediction. What was missing is a **denominator**: nothing in this
//! repository knew how reproducible its own Hausdorff number is, so "0.011 dB is
//! too small to see" was an opinion. This row makes it arithmetic.
//!
//! **There is no FCC reconstruction filter and there deliberately never was
//! one.** `benches/common/lattice.rs:41-48` states it: FCC's Delaunay complex is
//! the tetrahedral-octahedral honeycomb, two cell shapes rather than one, and the
//! reconstruction question this phase registered is BCC's. The same paragraph
//! states what R-163 needs instead — *"only a zero-set point sample from each
//! lattice, which `zero_set_hausdorff` consumes without caring how it was
//! produced"*. So this harness does not reconstruct at all, and §1 below is the
//! sample it uses instead.
//!
//! Everything here drives `common::lattice` and the public `isomesh` API.
//! `crates/isomesh/src/**` is untouched, no reference field is added, and no
//! golden hash can move.
//!
//! # Arms
//!
//! Two per field, eight fields, sixteen rows. `(field, lattice)` is the primary
//! key. The per-field comparison columns are stamped on **both** of a field's
//! rows — see §4 for why that differs from P-162's layout on purpose.
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `A3*` — BCC sites, its own Delaunay bond graph | nothing: the lattice P-162 already measured and the one `G` says is optimal among 3D lattices | **yes** |
//! | `D3` — FCC sites, its own split-Delaunay bond graph | the sampling lattice, at matched point density and matched bond count | no |
//!
//! The cubic lattice is built on every field and then **not measured**. It exists
//! only to supply the anchor count, because `lattice_grid`'s documented protocol
//! (`benches/common/lattice.rs:341-364`) is that the cubic arm — which anchored
//! on the box centre can realise only an odd number of sites per axis, with gaps
//! of 30% between attainable totals — must be built first and the other lattices
//! asked for *the count it realised*. Anchoring on `49³ = 117,649` is P-162's own
//! `TARGET_POINTS`, so the BCC arm here sits at exactly the density P-162
//! measured it at and the two CSVs are talking about the same grid.
//!
//! # Method
//!
//! ## 1. Where the zero-set point sample comes from, with no filter anywhere
//!
//! Each lattice is contoured on **its own Delaunay bond graph**: for every pair
//! of Delaunay-adjacent sites whose sampled values disagree in sign, the crossing
//! is placed on the bond by linear interpolation of the two *site values*. No
//! reconstruction filter is evaluated, no auxiliary grid is introduced, and
//! nothing lattice-specific enters beyond the bond set itself.
//!
//! That is not a convenience. It is exactly the vertex set a simplicial extractor
//! produces: `marching_tetrahedra.rs:249-255` places one vertex per cut edge and
//! says *"on a cut edge one endpoint is strictly negative and the other is `>= 0`,
//! so `a - b` is never zero and no epsilon guard is wanted"*. The offset used
//! here is the crate's own **centred** form,
//! `d = ((a + b)/2) / (a − b) ∈ [−0.5, 0.5]`, with the crossing at
//! `mid + d·(hi − lo)` — `cube.rs:221-225` and `cube.rs:232-234`, reproduced
//! bench-locally because `mod cube` is private (`lib.rs:143`) and R-059's reason
//! for centring it is numerical, not cosmetic.
//!
//! **The bond sets, and the integer equality that makes the comparison fair.**
//! Both lattices are enumerated in the integer coordinates
//! `Lattice::generator` generates (`benches/common/lattice.rs:176-189`): BCC is
//! `{k : k₀ ≡ k₁ ≡ k₂ (mod 2)}` at `scale/∛4` per integer step, FCC is
//! `{k : k₀+k₁+k₂ even}` at `scale/∛2`.
//!
//! | lattice | Voronoi cell | bonds per site | positive-half offsets |
//! |---|---|---|---|
//! | `A3*` | truncated octahedron, 14 facets | **7** | 4 × `(±1,±1,±1)` at `√3`, 3 × `(2,0,0)`-type at `2` |
//! | `D3` | rhombic dodecahedron, 12 facets | **7** | 6 × `(±1,±1,0)`-type at `√2`, 1 × `(0,0,2)` at `2` |
//!
//! BCC's Delaunay edges are its 14 Voronoi facet normals, halved. FCC's are its
//! 12, halved — **plus one**, because FCC's honeycomb contains one octahedron per
//! site and `benches/common/lattice.rs:638-648` splits each along *one fixed body
//! diagonal* to make the complex all-tetrahedral. That diagonal is an edge of the
//! split complex, it joins `s` to `s + (0,0,2)` through the octahedron centred at
//! the odd point `s + (0,0,1)`, and there is exactly one per site. So both
//! lattices carry **exactly seven Delaunay bonds per site**, and at matched point
//! density that means matched *bond* density: neither arm gets a denser crossing
//! set than the other. The equality is an integer one and is asserted, not hoped
//! for.
//!
//! The z-axis choice of diagonal is a genuine anisotropy of the FCC arm and it is
//! the module's documented choice, not this bench's invention; `case_table` is
//! insensitive to it, this bond set is not, and saying so is cheaper than
//! pretending otherwise.
//!
//! ## 2. The dB convention, taken from P-162 unchanged
//!
//! `experiment_p162.rs:111-128` argues it and this row must not re-argue it, or
//! the two rows' decibels stop being commensurable. `G` is a **second moment** —
//! a mean *squared* error — so the registered `0.011 dB` lives in the power
//! convention, `10·log₁₀` of a ratio of squared errors. Hausdorff distance is a
//! *linear* distance, an amplitude, so it enters as
//!
//! ```text
//!   measured_gap_db = 20·log₁₀(h_D3 / h_A3*) = 10·log₁₀((h_D3 / h_A3*)²)
//! ```
//!
//! — the same number twice, the second form making the commensurability obvious.
//! `20·log₁₀` on every ratio of distances on every row, including the scatter.
//!
//! The **sign** is the module's, and `benches/common/lattice.rs:142-164`
//! documents the trap: lower `G` is better, so
//! `Lattice::Fcc.gain_db_over(Lattice::Bcc)` is `+0.011143 dB` and reads *gain
//! available by moving from FCC to BCC*. Putting the worse lattice's Hausdorff in
//! the numerator aligns the measurement with it: positive means BCC came out
//! better, which is the direction the prediction is written in.
//!
//! ## 3. What varies between the repeats, and what does not
//!
//! **What does not vary:** the field, the box, the lattice, the site count, the
//! sample values, the bond set, the crossing set, and the seed. Every one is
//! computed once per arm and reused, so no repeat can differ in the geometry
//! being measured. A deterministic computation has zero scatter and a zero
//! denominator would make C1 trivially false, so a repeat that changed nothing
//! would be worthless.
//!
//! **What varies: the probe set.** `zero_set_hausdorff`'s
//! truth-to-reconstruction direction seeds `probes` points from a SplitMix64
//! stream on the fixed seed `0x1362_A3B5_D1E7_9F11`
//! (`benches/common/lattice.rs:1128-1130`), Newton-projects each onto the true
//! zero set, and takes the maximum distance to the nearest crossing. That
//! maximum is a **Monte-Carlo extreme-value statistic over a probe sample**, and
//! the module says so itself: *"a Hausdorff distance quoted without its probe
//! count is not reproducible"* (`benches/common/lattice.rs:1168-1171`). The seed
//! is `const` and this bench may not edit the module, so the lever on the probe
//! set is `probes`, and the repeats walk it: **2600, 2800, 3000, 3200, 3400**,
//! five distinct probe sets, identical schedule on both arms.
//!
//! **The limitation, stated because it decides the direction of the error.**
//! Those five sets are nested prefixes of one stream, so the spread is a monotone
//! sensitivity to the probe budget rather than a symmetric re-draw of an
//! independent sample. A full re-draw would scatter at least as much. So
//! `measurement_scatter_db` is a **lower bound** on the harness's scatter — which
//! makes C1 *harder* to satisfy, not easier, and a null that survives a
//! deliberately tight denominator is worth more than one that needed a loose one.
//!
//! The recon-to-truth direction is probe-independent, so if it ever dominated the
//! maximum the scatter would be exactly zero. `scatter_is_zero` is recorded for
//! that reason and is deliberately **not** asserted: a zero there falsifies C1
//! with its own arithmetic, which is a registered outcome and not a broken
//! fixture.
//!
//! ## 4. How the five repeats become the four registered numbers
//!
//! Per lattice `L`, `h_L(r)` for the five schedule entries `r`. Then
//!
//! ```text
//!   arm_scatter_db(L)      = 20·log₁₀( max_r h_L(r) / min_r h_L(r) )
//!   measurement_scatter_db = max over the two arms of arm_scatter_db
//!   gap_db(r)              = 20·log₁₀( h_D3(r) / h_A3*(r) )
//!   measured_gap_db        = median_r gap_db(r)
//!   gap_below_scatter      = |measured_gap_db| < measurement_scatter_db
//! ```
//!
//! `measurement_scatter_db` is the worse arm's, because the harness's scatter is
//! at least its worst arm's; each arm's own figure is the extra `arm_scatter_db`.
//! The gap is the **median of the five paired gaps**, not the ratio of two
//! medians: pairing the arms at a matched probe budget removes the common-mode
//! part of the budget trend, which is the sharper estimator and the honest one,
//! since the trend is shared machinery and not a property of either lattice.
//!
//! **Unlike P-162, the comparison columns are not relative-to-control per row.**
//! P-162's cubic row carried `measured_gain_db = 0` by construction. Doing that
//! here would put a zero in `measured_gap_db` on the `A3*` row, and
//! `|0| < scatter` is trivially true — a control row that votes for C1 for free.
//! So the single per-field FCC-versus-BCC comparison is stamped identically on
//! both of that field's rows, and `G`, `samples`, `points` and the error columns
//! are what distinguish them.
//!
//! **C1's decision rule over eight fields, fixed before the run.** C1 claims the
//! difference is *unresolvable*. One field on which it resolves is a field on
//! which somebody could distinguish the two lattices, so `c1_holds` is the
//! **conjunction** over the fields measured: every field must have
//! `gap_below_scatter`. It is a global verdict, stamped identically on all
//! sixteen rows, and `fields_below_scatter` / `fields_measured` are recorded so a
//! reader can apply any other rule to the same numbers.
//!
//! There is no C2 and no C3: the registration has one clause.
//!
//! ## 5. What the mean-square columns are for
//!
//! Hausdorff is a **maximum**; `G` predicts a **mean square**. So the same
//! crossing set is also reduced to an RMS of the linearised distance
//! `|f(p)| / ‖∇f(p)‖`, exactly as P-162 does (`experiment_p162.rs:499-518`) and
//! for the same reason: a first-order distance whose common second-order bias
//! cancels in the *ratio* between two arms measured on the same distance scale.
//! `rms_gap_db` is that ratio in the same `20·log₁₀` convention. It carries no
//! scatter column because it does not depend on the probe stream — which is
//! itself the point: the mean-square form of this measurement is exactly
//! reproducible and the max form is not.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration says **`SHARE: none`**, and it is right for a reason worth
//! writing down rather than skipping. A null here moves no stage of the pipeline
//! because it *forecloses* one: it says the FCC-versus-BCC question does not
//! repay a week of anybody's time. `crates/isomesh/src/**` is unchanged, no
//! reference field is added, no golden hash can move, and there is no landing
//! ticket to register from either verdict — a positive C1 lands nothing by
//! construction, and a negative C1 lands nothing either, because the thing it
//! would report is that this *harness* is more precise than expected, which is a
//! fact about `zero_set_hausdorff` and not about the crate.
//!
//! # A prediction, written down before the run
//!
//! At matched point density and matched bond count, the two bond sets do **not**
//! have matched bond *lengths*. In units of `scale`:
//!
//! ```text
//!   A3*:  (4·√3 + 3·2) / 7 / ∛4 = 1.16355 · scale
//!   D3:   (6·√2 + 1·2) / 7 / ∛2 = 1.18894 · scale
//! ```
//!
//! FCC's mean Delaunay bond is **2.18% longer**. The reported Hausdorff is
//! dominated by its truth-to-reconstruction direction — the covering radius of
//! the crossing set on the true surface, which scales with bond length — while
//! the recon-to-truth direction is a placement error two orders of magnitude
//! smaller. So bond length alone predicts
//!
//! ```text
//!   20·log₁₀(1.18894 / 1.16355) = +0.1873 dB
//! ```
//!
//! which is **seventeen times** the `0.011143 dB` the registration predicts, in
//! the same direction. `bond_length_gap_db` records it per field so the CSV
//! carries both predictions and a reader can see which one the measurement
//! matched. If the harness's scatter comes out under `0.19 dB`, C1 is FALSIFIED —
//! and falsified by *both* branches of its own falsifier at once: the harness is
//! more precise than expected, **and** Hausdorff error responds to something `G`
//! does not capture, namely the edge length of the Delaunay complex at fixed
//! point density. That is P-162 C2's question answered from the other side, which
//! is exactly what the falsifier says a resolvable difference would mean.
//!
//! # Vacuity controls
//!
//! All six run before the first `run.record`, and every panic message starts
//! `VOID: `.
//!
//! - **The registration's own control: at least five repeated runs of the same
//!   lattice.** The schedule must hold at least five entries and they must be
//!   strictly increasing, so the five probe sets genuinely differ and "below
//!   scatter" has a denominator. Columns: `repeats`, `probe_counts`,
//!   `hausdorff_by_probes`.
//! - **The prediction under test is the registered one.**
//!   `Lattice::Fcc.gain_db_over(Lattice::Bcc)` must agree with the registered
//!   `0.011 dB` to `5e-4`, or C1 would be scored against a number nobody
//!   registered. Columns: `predicted_gap_db`, `G`.
//! - **Matched point density.** `|samples_D3 − samples_A3*| / samples_A3*` must
//!   be at most 5%, both counts recorded — P-162's tolerance, so the two rows
//!   agree about what "matched" means. Columns: `samples`, `density_mismatch`.
//! - **Matched bond density.** Seven positive-half offsets on both lattices, an
//!   integer equality; and the *realised* bonds per site, after the box clipped
//!   the sites, must also agree within 5%. Without it the arm with more crossings
//!   wins a covering-radius comparison for a reason that is not its lattice.
//!   Columns: `bonds`, `bonds_per_site`, `bond_mismatch`.
//! - **Both arms measured a real surface.** At least 64 crossings, and a
//!   Hausdorff and an RMS that are strictly positive and finite: a maximum over a
//!   handful of points is not a surface and a ratio of two zeros is not a gap
//!   (M-44). Columns: `points`, `hausdorff`, `rms_error`.
//! - **The integer decoding is exact.** Every site's decoded coordinate must
//!   satisfy its lattice's parity condition and reproduce the site to `1e-6` of
//!   one integer step. A wrong decode produces a wrong bond set silently, and a
//!   silently wrong bond set is a wrong covering radius. Column:
//!   `coord_residual`.

mod common;

use common::lattice::{Lattice, LatticeGrid, lattice_grid, zero_set_hausdorff};
use isomesh::Sdf;
use isomesh::fields::ReferenceField;

// ─── the configuration, all of it derived in the header ─────────────────────

/// Sites the **cubic** anchor is asked for: `49³`.
///
/// P-162's own `TARGET_POINTS`, so the BCC arm here sits at exactly the density
/// `docs/experiments/p-162.csv` measured it at. The cubic grid is built for this
/// number and then discarded; the protocol reason is in the header's Arms
/// section.
const TARGET_POINTS: usize = 117_649;

/// The probe budgets of the five repeats, in order.
///
/// Five distinct probe sets, centred on 3000, identical on both arms. Nested
/// prefixes of one fixed stream — the header says why that makes the resulting
/// scatter a lower bound and why a lower bound is the conservative direction.
const PROBE_SCHEDULE: [usize; 5] = [2_600, 2_800, 3_000, 3_200, 3_400];

/// Repeats, which is the schedule's length.
const REPEATS: usize = PROBE_SCHEDULE.len();

/// The registration's vacuity control: *"at least five repeated runs of the same
/// lattice"*.
const MIN_REPEATS: usize = 5;

/// Largest site-count gap between the two arms this comparison will accept.
///
/// P-162's `DENSITY_TOLERANCE`, unchanged, so both rows mean the same thing by
/// "matched point density".
const DENSITY_TOLERANCE: f64 = 0.05;

/// Fewest crossings an arm may report and still be describing a surface.
///
/// P-162's `MIN_CROSSINGS`.
const MIN_CROSSINGS: usize = 64;

/// Delaunay bonds per site, and it is **7 on both lattices**.
///
/// BCC: its truncated-octahedral Voronoi cell has 14 facets, so 7 bonds per
/// site. FCC: its rhombic dodecahedron has 12, so 6, plus the one split diagonal
/// of the one octahedron per site. This integer equality is what makes a
/// covering-radius comparison at matched point density also a comparison at
/// matched bond density.
const DELAUNAY_BONDS_PER_SITE: usize = 7;

/// The dB figure the registration predicts, quoted from it.
const REGISTERED_GAP_DB: f64 = 0.011;

/// Agreement required between the registration's rounded dB and the module's
/// computed `0.011143`.
const GAP_DB_TOLERANCE: f64 = 5e-4;

/// The amplitude dB constant, `20`.
///
/// Not `10`, and P-162 argues it at `experiment_p162.rs:111-128`: `G` is a mean
/// *squared* error, Hausdorff is a linear distance, so
/// `20·log₁₀(h₁/h₂) = 10·log₁₀((h₁/h₂)²)` is the commensurable form.
const AMPLITUDE_DB: f64 = 20.0;

/// How far a decoded integer coordinate may sit from its site, in integer steps.
const COORD_TOLERANCE: f64 = 1e-6;

/// Marker for "no site at this integer coordinate" in the dense lookup.
const ABSENT: u32 = u32::MAX;

// ─── the lattice geometry this bench needs and the module does not export ───

/// The positive half of one lattice's Delaunay bond set, in integer coordinates.
///
/// Positive half means one offset per bond: the first non-zero component is
/// positive, so walking every site against every offset visits each bond exactly
/// once. Derived from each lattice's Voronoi cell — a Delaunay edge is a pair of
/// sites whose Voronoi cells share a facet — plus, for FCC, the one split
/// diagonal per octahedron that `benches/common/lattice.rs:638-648` chooses.
///
/// The cubic entry is the cube's three axis directions and is correct rather than
/// merely present: this bench never contours the cubic lattice, but the enum has
/// three variants and an arm that is wrong by omission is worse than one that is
/// right and unused.
fn bond_offsets(lattice: Lattice) -> &'static [[i64; 3]] {
    match lattice {
        Lattice::Cubic => &[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        Lattice::Bcc => &[
            [1, 1, 1],
            [1, 1, -1],
            [1, -1, 1],
            [1, -1, -1],
            [2, 0, 0],
            [0, 2, 0],
            [0, 0, 2],
        ],
        Lattice::Fcc => &[
            [1, 1, 0],
            [1, -1, 0],
            [1, 0, 1],
            [1, 0, -1],
            [0, 1, 1],
            [0, 1, -1],
            [0, 0, 2],
        ],
    }
}

/// World distance of one step of the lattice's integer coordinate.
///
/// `Lattice::generator`'s integer bases carry divisors `1`, `∛4` and `∛2`
/// (`benches/common/lattice.rs:176-189`), and the rows are then multiplied by
/// `LatticeGrid::scale`, so one integer step is `scale / divisor`.
fn integer_unit(lattice: Lattice, scale: f64) -> f64 {
    match lattice {
        Lattice::Cubic => scale,
        Lattice::Bcc => scale / 4f64.cbrt(),
        Lattice::Fcc => scale / 2f64.cbrt(),
    }
}

/// Does `k` lie on the lattice?
///
/// BCC is `{k : k₀ ≡ k₁ ≡ k₂ (mod 2)}` and FCC is `{k : k₀+k₁+k₂ even}`. The
/// remainder is Euclidean rather than `%`, because `(-3) % 2` is `-1` in Rust and
/// a parity test written with `%` accepts the wrong half of the lattice on the
/// negative side of the box.
fn on_lattice(lattice: Lattice, k: [i64; 3]) -> bool {
    match lattice {
        Lattice::Cubic => true,
        Lattice::Bcc => (k[0] - k[1]).rem_euclid(2) == 0 && (k[1] - k[2]).rem_euclid(2) == 0,
        Lattice::Fcc => (k[0] + k[1] + k[2]).rem_euclid(2) == 0,
    }
}

/// Mean Delaunay bond length of a lattice, in world units.
///
/// A lattice constant, not a measurement: the mean of `|offset|` over the
/// positive-half bond set, times one integer step. `1.16355·scale` for BCC and
/// `1.18894·scale` for FCC, which is the header's `+0.1873 dB` prediction.
fn mean_bond_length(lattice: Lattice, unit: f64) -> f64 {
    let offsets = bond_offsets(lattice);
    let total: f64 = offsets
        .iter()
        .map(|o| {
            let s = (o[0] * o[0] + o[1] * o[1] + o[2] * o[2]) as f64;
            s.sqrt()
        })
        .sum();
    unit * total / offsets.len() as f64
}

/// The integer coordinates of a grid's sites, and an `O(1)` lookup from an
/// integer coordinate back to the site.
///
/// The table is dense because it can be: the integer span is set by the point
/// density and not by the domain, so it is about `78³` entries for BCC and `62³`
/// for FCC whatever field is being sampled, a couple of megabytes against
/// 823,543 bond lookups per arm.
#[derive(Debug)]
struct Sites {
    /// One integer coordinate per site, parallel to `LatticeGrid::sites`.
    coords: Vec<[i64; 3]>,
    /// Lowest integer coordinate on each axis.
    origin: [i64; 3],
    /// Extent of the table on each axis.
    dims: [usize; 3],
    /// Site index per integer coordinate, or [`ABSENT`].
    table: Vec<u32>,
    /// Largest distance, in integer steps, between a site and the position its
    /// decoded coordinate names.
    residual: f64,
}

impl Sites {
    /// The site at integer coordinate `k`, or `None` if the box clipped it away.
    fn at(&self, k: [i64; 3]) -> Option<usize> {
        let d0 = k[0] - self.origin[0];
        let d1 = k[1] - self.origin[1];
        let d2 = k[2] - self.origin[2];
        if d0 < 0 || d1 < 0 || d2 < 0 {
            return None;
        }
        let (d0, d1, d2) = (d0 as usize, d1 as usize, d2 as usize);
        if d0 >= self.dims[0] || d1 >= self.dims[1] || d2 >= self.dims[2] {
            return None;
        }
        let found = self.table[d0 + self.dims[0] * (d1 + self.dims[1] * d2)];
        (found != ABSENT).then_some(found as usize)
    }
}

/// Decode every site of `grid` into its integer coordinate and build the lookup.
///
/// The grid is anchored on the box centre, which is always a lattice site
/// (`benches/common/lattice.rs:272-285`), so the coordinate is
/// `round((site − centre) / unit)` and the rounding residual is a check on the
/// whole decode rather than a tolerance being spent.
///
/// # Panics
///
/// If a decoded coordinate misses its site by more than [`COORD_TOLERANCE`] of an
/// integer step, if a coordinate is off the lattice, or if two sites decode to
/// the same coordinate. All three are the same fault — a wrong `unit` or a wrong
/// parity — and all three would silently produce a wrong bond set.
fn decode(grid: &LatticeGrid, unit: f64) -> Sites {
    assert!(
        grid.sites.len() < ABSENT as usize,
        "VOID: {} holds {} sites, which does not fit the dense lookup's u32 index",
        grid.lattice.name(),
        grid.sites.len()
    );
    let centre = [
        f64::midpoint(grid.lo[0], grid.hi[0]),
        f64::midpoint(grid.lo[1], grid.hi[1]),
        f64::midpoint(grid.lo[2], grid.hi[2]),
    ];

    let mut coords: Vec<[i64; 3]> = Vec::with_capacity(grid.sites.len());
    let mut residual = 0.0f64;
    for site in &grid.sites {
        let mut k = [0i64; 3];
        for (axis, slot) in k.iter_mut().enumerate() {
            let exact = (site[axis] - centre[axis]) / unit;
            let rounded = exact.round();
            residual = residual.max((exact - rounded).abs());
            *slot = rounded as i64;
        }
        assert!(
            on_lattice(grid.lattice, k),
            "VOID: site {site:?} of {} decodes to {k:?}, which is not on that lattice, so the \
             bond set built from these coordinates is not the lattice's Delaunay complex",
            grid.lattice.name()
        );
        coords.push(k);
    }
    assert!(
        residual <= COORD_TOLERANCE,
        "VOID: decoding {}'s sites at one integer step of {unit} leaves a rounding residual of \
         {residual}, over the {COORD_TOLERANCE} this decode allows — the integer coordinate \
         system is wrong and every bond derived from it is wrong with it",
        grid.lattice.name()
    );

    let mut origin = [i64::MAX; 3];
    let mut top = [i64::MIN; 3];
    for k in &coords {
        for ((slot, cap), value) in origin.iter_mut().zip(top.iter_mut()).zip(k.iter()) {
            *slot = (*slot).min(*value);
            *cap = (*cap).max(*value);
        }
    }
    let dims = [
        (top[0] - origin[0] + 1) as usize,
        (top[1] - origin[1] + 1) as usize,
        (top[2] - origin[2] + 1) as usize,
    ];

    let mut table = vec![ABSENT; dims[0] * dims[1] * dims[2]];
    for (i, k) in coords.iter().enumerate() {
        let at = (k[0] - origin[0]) as usize
            + dims[0] * ((k[1] - origin[1]) as usize + dims[1] * (k[2] - origin[2]) as usize);
        assert_eq!(
            table[at],
            ABSENT,
            "VOID: sites {} and {i} of {} decode to the same integer coordinate {k:?}, so the \
             lookup would hide one of them and the bond count would be wrong",
            table[at],
            grid.lattice.name()
        );
        table[at] = i as u32;
    }

    Sites {
        coords,
        origin,
        dims,
        table,
        residual,
    }
}

// ─── the measurement ────────────────────────────────────────────────────────

/// One lattice's arm on one field.
#[derive(Debug)]
struct Arm {
    /// Which lattice this arm sampled on.
    lattice: Lattice,
    /// Sites the lattice realised in the box — the matched-density number, read
    /// from `LatticeGrid::sites` and never from [`TARGET_POINTS`].
    samples: usize,
    /// The factor the unit-volume generator rows were multiplied by.
    scale: f64,
    /// World distance of one integer step.
    unit: f64,
    /// Mean Delaunay bond length, a lattice constant.
    bond_length: f64,
    /// Largest decoding residual over the sites, in integer steps.
    coord_residual: f64,
    /// Delaunay bonds with both endpoints inside the box.
    bonds: u64,
    /// Those bonds per site — under seven by the boundary deficit.
    bonds_per_site: f64,
    /// Crossings of this arm's bond graph: one per bond whose endpoints disagree
    /// in sign.
    points: usize,
    /// Hausdorff per repeat, in [`PROBE_SCHEDULE`] order.
    hausdorff: [f64; REPEATS],
    /// Median of them, the headline.
    h_median: f64,
    /// Smallest of them.
    h_min: f64,
    /// Largest of them.
    h_max: f64,
    /// `20·log₁₀(h_max / h_min)` — this arm's own scatter.
    scatter_db: f64,
    /// RMS of the linearised distance over the crossings — the mean-square form
    /// `G` actually predicts.
    rms: f64,
    /// Mean of the same linearised distance.
    mean: f64,
    /// Largest of it.
    worst: f64,
}

/// One field: both arms, and the comparison between them.
#[derive(Debug)]
struct FieldRow {
    /// The reference field's name.
    field: &'static str,
    /// `|samples_fcc − samples_bcc| / samples_bcc`.
    mismatch: f64,
    /// The same for realised bonds per site.
    bond_mismatch: f64,
    /// `h_fcc / h_bcc` on the medians.
    ratio: f64,
    /// The five paired gaps, in schedule order.
    gaps: [f64; REPEATS],
    /// Median of them — `measured_gap_db`.
    gap_median: f64,
    /// Smallest of them.
    gap_min: f64,
    /// Largest of them.
    gap_max: f64,
    /// The harness's scatter: the worse arm's `scatter_db`.
    scatter_db: f64,
    /// `|gap_median| < scatter_db`.
    below: bool,
    /// The bond-length-only prediction, `20·log₁₀(L_fcc / L_bcc)`.
    bond_length_gap_db: f64,
    /// The mean-square form of the same comparison.
    rms_gap_db: f64,
    /// The control arm.
    bcc: Arm,
    /// The challenger.
    fcc: Arm,
}

/// Linearised distance from `p` to the field's zero set: `|f(p)| / ‖∇f(p)‖`.
///
/// First order, and that is enough: it is used only inside a *ratio* between two
/// arms whose crossings sit at the same distance scale, where the common
/// second-order bias cancels. P-162 uses the identical form
/// (`experiment_p162.rs:499-518`), so the two rows' RMS columns are the same
/// quantity.
fn linear_residual<F>(field: &F, p: [f64; 3]) -> f64
where
    F: Sdf<Scalar = f64>,
{
    let value = field.sample(p);
    let g = field.gradient(p);
    let norm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    assert!(
        norm > 0.0 && norm.is_finite(),
        "the field's gradient is {g:?} at {p:?}, so the linearised distance to its zero set is \
         undefined there"
    );
    value.abs() / norm
}

/// Contour one lattice's Delaunay bond graph, returning the crossings and the
/// number of bonds that had both endpoints.
///
/// One crossing per bond whose endpoint values disagree in sign, placed by the
/// crate's own centred offset. A bond whose partner was clipped away by the box
/// is not a bond of this grid and is not counted; that is the same rule on both
/// arms and `bond_mismatch` is the control on it.
fn crossings(grid: &LatticeGrid, sites: &Sites, values: &[f64]) -> (Vec<[f64; 3]>, u64) {
    let offsets = bond_offsets(grid.lattice);
    let mut bonds = 0u64;
    let mut out: Vec<[f64; 3]> = Vec::new();
    for (i, k) in sites.coords.iter().enumerate() {
        let a = grid.sites[i];
        let fa = values[i];
        for o in offsets {
            let Some(j) = sites.at([k[0] + o[0], k[1] + o[1], k[2] + o[2]]) else {
                continue;
            };
            bonds += 1;
            let fb = values[j];
            if (fa < 0.0) == (fb < 0.0) {
                continue;
            }
            // `cube.rs:221-225`: one endpoint is strictly negative and the other
            // is `>= 0`, so `fa - fb` is never zero and no epsilon guard is
            // wanted. Centred on the bond's midpoint, which is R-059's frame.
            let d = f64::midpoint(fa, fb) / (fa - fb);
            assert!(
                (-0.5..=0.5).contains(&d),
                "the centred crossing offset on a cut bond is {d}, outside the half-edge it is \
                 defined on: values {fa} and {fb}"
            );
            let b = grid.sites[j];
            out.push([
                f64::midpoint(a[0], b[0]) + (b[0] - a[0]) * d,
                f64::midpoint(a[1], b[1]) + (b[1] - a[1]) * d,
                f64::midpoint(a[2], b[2]) + (b[2] - a[2]) * d,
            ]);
        }
    }
    (out, bonds)
}

/// Measure one arm: decode, sample, contour, then five Hausdorffs on five probe
/// sets.
///
/// The geometry is computed once and reused across the repeats, which is what
/// makes the probe set the only thing that varies.
fn arm<F>(field: &F, grid: &LatticeGrid) -> Arm
where
    F: Sdf<Scalar = f64>,
{
    let unit = integer_unit(grid.lattice, grid.scale);
    let sites = decode(grid, unit);

    let values: Vec<f64> = grid.sites.iter().map(|s| field.sample(*s)).collect();
    for (i, v) in values.iter().enumerate() {
        assert!(
            v.is_finite(),
            "the field samples {v} at site {i} of {}, {:?} — a non-finite value has no sign and \
             every bond touching it would be silently uncut",
            grid.lattice.name(),
            grid.sites[i]
        );
    }

    let (points, bonds) = crossings(grid, &sites, &values);

    let mut hausdorff = [0.0f64; REPEATS];
    for (slot, probes) in hausdorff.iter_mut().zip(PROBE_SCHEDULE) {
        *slot = zero_set_hausdorff(field, &points, probes);
    }
    let mut sorted = hausdorff;
    sorted.sort_unstable_by(f64::total_cmp);

    let mut sum = 0.0f64;
    let mut sum_sq = 0.0f64;
    let mut worst = 0.0f64;
    for p in &points {
        let r = linear_residual(field, *p);
        sum += r;
        sum_sq += r * r;
        worst = worst.max(r);
    }
    let count = points.len() as f64;

    Arm {
        lattice: grid.lattice,
        samples: grid.sites.len(),
        scale: grid.scale,
        unit,
        bond_length: mean_bond_length(grid.lattice, unit),
        coord_residual: sites.residual,
        bonds,
        bonds_per_site: bonds as f64 / grid.sites.len() as f64,
        points: points.len(),
        hausdorff,
        h_median: sorted[REPEATS / 2],
        h_min: sorted[0],
        h_max: sorted[REPEATS - 1],
        scatter_db: AMPLITUDE_DB * (sorted[REPEATS - 1] / sorted[0]).log10(),
        rms: (sum_sq / count).sqrt(),
        mean: sum / count,
        worst,
    }
}

/// Measure both arms on one field and derive the comparison.
fn measure<F>(name: &'static str, field: &F) -> FieldRow
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();

    // The cubic grid is built for its realised count and then dropped. The
    // module's protocol is that the coarse-grained lattice anchors the density
    // (`benches/common/lattice.rs:341-364`), and `49³` is P-162's own target.
    let anchor = lattice_grid(Lattice::Cubic, lo, hi, TARGET_POINTS)
        .sites
        .len();
    let bcc = lattice_grid(Lattice::Bcc, lo, hi, anchor);
    let fcc = lattice_grid(Lattice::Fcc, lo, hi, anchor);

    let bcc_arm = arm(field, &bcc);
    let fcc_arm = arm(field, &fcc);

    let mismatch = (fcc_arm.samples as f64 - bcc_arm.samples as f64).abs() / bcc_arm.samples as f64;
    let bond_mismatch =
        (fcc_arm.bonds_per_site - bcc_arm.bonds_per_site).abs() / bcc_arm.bonds_per_site;

    let mut gaps = [0.0f64; REPEATS];
    for (slot, (h_fcc, h_bcc)) in gaps
        .iter_mut()
        .zip(fcc_arm.hausdorff.iter().zip(bcc_arm.hausdorff.iter()))
    {
        *slot = AMPLITUDE_DB * (h_fcc / h_bcc).log10();
    }
    let mut sorted = gaps;
    sorted.sort_unstable_by(f64::total_cmp);

    let scatter_db = bcc_arm.scatter_db.max(fcc_arm.scatter_db);
    let gap_median = sorted[REPEATS / 2];

    FieldRow {
        field: name,
        mismatch,
        bond_mismatch,
        ratio: fcc_arm.h_median / bcc_arm.h_median,
        gaps,
        gap_median,
        gap_min: sorted[0],
        gap_max: sorted[REPEATS - 1],
        scatter_db,
        below: gap_median.abs() < scatter_db,
        bond_length_gap_db: AMPLITUDE_DB * (fcc_arm.bond_length / bcc_arm.bond_length).log10(),
        rms_gap_db: AMPLITUDE_DB * (fcc_arm.rms / bcc_arm.rms).log10(),
        bcc: bcc_arm,
        fcc: fcc_arm,
    }
}

/// Five floats joined by `|`, because `Run::record` refuses a comma.
fn joined(values: &[f64]) -> String {
    values
        .iter()
        .map(|v| format!("{v:.9}"))
        .collect::<Vec<String>>()
        .join("|")
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-163");

    common::experiment::run(prereg, |run| {
        let predicted = Lattice::Fcc.gain_db_over(Lattice::Bcc);
        println!(
            "predicted gap   G(D3) = {:.9}  G(A3*) = {:.9}  ->  {predicted:+.6} dB  \
             (MSE excess {:.4}%)",
            Lattice::Fcc.g(),
            Lattice::Bcc.g(),
            // `G(D3)/G(A3*) - 1`: FCC's excess over the lattice it loses to,
            // which is the reading the positive dB figure is the log of.
            100.0 * (Lattice::Fcc.g() / Lattice::Bcc.g() - 1.0)
        );
        println!(
            "bond sets       A3*: {} offsets, mean |o| {:.6}  |  D3: {} offsets, mean |o| {:.6}",
            bond_offsets(Lattice::Bcc).len(),
            mean_bond_length(Lattice::Bcc, 1.0),
            bond_offsets(Lattice::Fcc).len(),
            mean_bond_length(Lattice::Fcc, 1.0)
        );
        println!(
            "probe schedule  {:?} over {REPEATS} repeats, identical on both arms",
            PROBE_SCHEDULE
        );
        println!();

        let mut rows: Vec<FieldRow> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            // An inline block per field, not a closure: a `return` in here would
            // return from `main` and the run would stop at the first field
            // (M-253).
            let row = measure(name, &field);
            println!(
                "{:>15}  samples {:>7} / {:>7} ({:.3}% apart)  bonds/site {:.4} / {:.4}  \
                 points {:>6} / {:>6}",
                row.field,
                row.bcc.samples,
                row.fcc.samples,
                row.mismatch * 100.0,
                row.bcc.bonds_per_site,
                row.fcc.bonds_per_site,
                row.bcc.points,
                row.fcc.points
            );
            println!(
                "{:>15}  h {:.6e} / {:.6e}  gap {:+.4} dB  scatter {:.4} dB  \
                 (arms {:.4} / {:.4})  below {}",
                "",
                row.bcc.h_median,
                row.fcc.h_median,
                row.gap_median,
                row.scatter_db,
                row.bcc.scatter_db,
                row.fcc.scatter_db,
                row.below
            );
            println!(
                "{:>15}  rms {:.6e} / {:.6e}  {:+.4} dB   bond-length prediction {:+.4} dB   \
                 gaps {:+.4}..{:+.4} dB",
                "",
                row.bcc.rms,
                row.fcc.rms,
                row.rms_gap_db,
                row.bond_length_gap_db,
                row.gap_min,
                row.gap_max
            );
            rows.push(row);
        });
        println!();

        // ── the vacuity controls, before any verdict is reported ─────────────

        // 1. The registration's own control: at least five repeated runs of the
        //    same lattice, on genuinely different probe sets. The probe sets are
        //    *counted* rather than assumed from the schedule's length: two equal
        //    entries would be one probe set measured twice, and a scatter
        //    estimated from that is a zero that could not have been non-zero
        //    (M-44).
        let distinct_probe_sets = PROBE_SCHEDULE
            .iter()
            .enumerate()
            .filter(|(seen, probes)| !PROBE_SCHEDULE[..*seen].contains(probes))
            .count();
        assert!(
            distinct_probe_sets >= MIN_REPEATS,
            "VOID: the schedule {PROBE_SCHEDULE:?} holds only {distinct_probe_sets} distinct \
             probe counts against the {MIN_REPEATS} repeated runs of the same lattice the \
             registration requires — 'below scatter' would have no denominator"
        );
        assert!(
            PROBE_SCHEDULE.windows(2).all(|w| w[0] < w[1]),
            "VOID: the probe schedule {PROBE_SCHEDULE:?} is not strictly increasing, so two \
             repeats would draw the identical probe set and the scatter would be a zero that \
             could not have been non-zero (M-44)"
        );

        // 2. The number C1 is scored against is the number that was registered.
        assert!(
            (predicted - REGISTERED_GAP_DB).abs() < GAP_DB_TOLERANCE,
            "VOID: the module computes {predicted:.9} dB from G(D3)/G(A3*) while the \
             registration predicts {REGISTERED_GAP_DB} dB — C1 would be scored against a \
             prediction nobody registered"
        );

        // 3. Matched point density, reported as counts — P-162's tolerance, so
        //    both rows mean the same thing by "matched".
        for row in &rows {
            assert!(
                row.mismatch <= DENSITY_TOLERANCE,
                "VOID: {}: the arms hold {} BCC sites against {} FCC sites, a {:.3}% gap \
                 against the {:.1}% this comparison allows — at that gap the row is a \
                 resolution change wearing a lattice's name",
                row.field,
                row.bcc.samples,
                row.fcc.samples,
                row.mismatch * 100.0,
                DENSITY_TOLERANCE * 100.0
            );
        }

        // 4. Matched bond density: the integer equality first, then the realised
        //    one after the box clipped the sites.
        for lattice in [Lattice::Bcc, Lattice::Fcc] {
            assert_eq!(
                bond_offsets(lattice).len(),
                DELAUNAY_BONDS_PER_SITE,
                "VOID: {} carries {} positive-half Delaunay offsets against the \
                 {DELAUNAY_BONDS_PER_SITE} both lattices must have — without that equality the \
                 arm with more bonds wins a covering-radius comparison for a reason that is not \
                 its lattice",
                lattice.name(),
                bond_offsets(lattice).len()
            );
        }
        for row in &rows {
            assert!(
                row.bond_mismatch <= DENSITY_TOLERANCE,
                "VOID: {}: {:.4} realised bonds per BCC site against {:.4} per FCC site, a \
                 {:.3}% gap against the {:.1}% this comparison allows — the crossing sets are \
                 not at matched density",
                row.field,
                row.bcc.bonds_per_site,
                row.fcc.bonds_per_site,
                row.bond_mismatch * 100.0,
                DENSITY_TOLERANCE * 100.0
            );
        }

        // 5. Both arms measured a real surface, so no ratio here is a ratio of
        //    two zeros (M-44).
        for row in &rows {
            for side in [&row.bcc, &row.fcc] {
                assert!(
                    side.points >= MIN_CROSSINGS,
                    "VOID: {} on {}: only {} crossings, under the {MIN_CROSSINGS} this harness \
                     will call a surface — a Hausdorff maximum over a handful of points is not \
                     an error measurement",
                    row.field,
                    side.lattice.name(),
                    side.points
                );
                for (probes, h) in PROBE_SCHEDULE.iter().zip(side.hausdorff.iter()) {
                    assert!(
                        *h > 0.0 && h.is_finite(),
                        "VOID: {} on {} at {probes} probes: Hausdorff {h} — a zero or \
                         non-finite error makes the gap and the scatter meaningless",
                        row.field,
                        side.lattice.name()
                    );
                }
                assert!(
                    side.rms > 0.0 && side.rms.is_finite(),
                    "VOID: {} on {}: RMS linearised distance {} — the mean-square reading of \
                     the same comparison would be evaluated on a zero",
                    row.field,
                    side.lattice.name(),
                    side.rms
                );
            }
        }

        // 6. The integer decoding is exact, so the bond set really is the
        //    lattice's Delaunay complex. `decode` asserts per site; this is the
        //    aggregate, so a reader of the CSV can see it was asked.
        for row in &rows {
            for side in [&row.bcc, &row.fcc] {
                assert!(
                    side.coord_residual <= COORD_TOLERANCE,
                    "VOID: {} on {}: decoding residual {} over {COORD_TOLERANCE} integer steps",
                    row.field,
                    side.lattice.name(),
                    side.coord_residual
                );
            }
        }

        // ── the verdict, global, with its arithmetic ─────────────────────────

        let below = rows.iter().filter(|r| r.below).count();

        // C1: the measured FCC-BCC Hausdorff difference is below the measurement
        // scatter of the harness itself. The conjunction over the fields, fixed
        // before the run: the claim is unresolvability, so one field on which the
        // gap resolves is a field on which somebody could tell the two lattices
        // apart.
        let c1 = below == rows.len();

        println!(
            "C1  |gap| under the harness scatter on {below} of {} fields (needs all) -> {c1}",
            rows.len()
        );
        for row in &rows {
            println!(
                "      {:>15}  |{:+.4}| dB  vs scatter {:.4} dB  -> {}   \
                 (bond-length prediction {:+.4} dB, registered {predicted:+.6} dB)",
                row.field, row.gap_median, row.scatter_db, row.below, row.bond_length_gap_db
            );
        }
        println!();

        // ── the rows ────────────────────────────────────────────────────────

        for row in &rows {
            for side in [&row.bcc, &row.fcc] {
                let is_control = matches!(side.lattice, Lattice::Bcc);
                let over_scatter = if row.scatter_db > 0.0 {
                    row.gap_median.abs() / row.scatter_db
                } else {
                    f64::INFINITY
                };

                run.record(&[
                    ("lattice", side.lattice.name().to_string()),
                    ("G", format!("{:.9}", side.lattice.g())),
                    ("predicted_gap_db", format!("{predicted:.6}")),
                    ("measured_gap_db", format!("{:.6}", row.gap_median)),
                    ("measurement_scatter_db", format!("{:.6}", row.scatter_db)),
                    ("gap_below_scatter", row.below.to_string()),
                    ("c1_holds", c1.to_string()),
                    // ── extras (M-273) ──
                    ("arm_scatter_db", format!("{:.6}", side.scatter_db)),
                    (
                        "bond_length_gap_db",
                        format!("{:.6}", row.bond_length_gap_db),
                    ),
                    ("bond_mismatch", format!("{:.6}", row.bond_mismatch)),
                    ("bonds", side.bonds.to_string()),
                    ("bonds_per_site", format!("{:.6}", side.bonds_per_site)),
                    ("coord_residual", format!("{:.3e}", side.coord_residual)),
                    ("density_mismatch", format!("{:.6}", row.mismatch)),
                    ("field", row.field.to_string()),
                    ("fields_below_scatter", below.to_string()),
                    ("fields_measured", rows.len().to_string()),
                    ("gap_db_by_repeat", joined(&row.gaps)),
                    ("gap_db_max", format!("{:.6}", row.gap_max)),
                    ("gap_db_min", format!("{:.6}", row.gap_min)),
                    ("gap_over_scatter", format!("{over_scatter:.6}")),
                    ("gap_range_db", format!("{:.6}", row.gap_max - row.gap_min)),
                    ("hausdorff", format!("{:.9}", side.h_median)),
                    ("hausdorff_by_probes", joined(&side.hausdorff)),
                    ("hausdorff_max", format!("{:.9}", side.h_max)),
                    ("hausdorff_min", format!("{:.9}", side.h_min)),
                    ("hausdorff_ratio", format!("{:.6}", row.ratio)),
                    ("is_control", is_control.to_string()),
                    ("lattice_scale", format!("{:.9}", side.scale)),
                    ("lattice_unit", format!("{:.9}", side.unit)),
                    ("mean_bond_length", format!("{:.9}", side.bond_length)),
                    ("mean_linear_error", format!("{:.9}", side.mean)),
                    ("points", side.points.to_string()),
                    (
                        "probe_counts",
                        PROBE_SCHEDULE
                            .iter()
                            .map(usize::to_string)
                            .collect::<Vec<String>>()
                            .join("|"),
                    ),
                    ("probes_max", PROBE_SCHEDULE[REPEATS - 1].to_string()),
                    ("probes_min", PROBE_SCHEDULE[0].to_string()),
                    ("repeats", REPEATS.to_string()),
                    ("rms_error", format!("{:.9}", side.rms)),
                    ("rms_gap_db", format!("{:.6}", row.rms_gap_db)),
                    ("samples", side.samples.to_string()),
                    ("scatter_is_zero", (row.scatter_db <= 0.0).to_string()),
                    ("target_points", TARGET_POINTS.to_string()),
                    ("worst_linear_error", format!("{:.9}", side.worst)),
                ]);
            }
        }
    });
}
