//! **P-58 — a discrete Morse census under exact ties.**
//!
//! Ticket: R-056. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p58
//! ```
//!
//! Writes `docs/experiments/p-58.csv`.
//!
//! # What is transcribed, and from where
//!
//! Robins, Wood & Sheppard, *Theory and Algorithms for Constructing Discrete
//! Morse Complexes from Grayscale Digital Images*, `10.1109/tpami.2011.95`. The
//! corpus copy is complete, so everything below is quoted rather than
//! reconstructed, and cited by number:
//!
//! - **§2.1**, the model. Voxels are the **0-cells**; higher cells are the unit
//!   edges, unit squares and unit cubes. `α < β` iff `dim α < dim β` and `α`'s
//!   vertices are a subset of `β`'s. That is the whole cell structure, and it is
//!   what `STAR_VERTS` and `STAR_FACES` below encode.
//! - **§3.1, Eq. (7)**, the lower star:
//!   `L(x) = { α ∈ K | x ∈ α and g(x) = max_{y ∈ α} g(y) }`, and *"the lower
//!   stars of all voxels `x ∈ D` form a disjoint partition of `K`"*. That
//!   partition is `cells_total` below, and it is asserted rather than trusted.
//! - **§3.1**, the ordering: `G(α)` for `α ∈ L(x)` with vertices
//!   `{x, y_1 … y_k}` is `(g(x), g(y_1), …, g(y_k))` listed in **decreasing**
//!   `g`-order and compared **lexicographically**. Both priority queues order by
//!   `G`. The paper's own worked lower star in Fig. 1 lists
//!   `6, 61, 62, 63, 6431, 6532` in that order, which pins the convention: a
//!   shorter tuple that is a prefix of a longer one is the smaller.
//! - **§3.1, Algorithm 1**, `ProcessLowerStars`, all 27 lines, transcribed line
//!   for line in `Star::process` with the pseudocode's line numbers in comments.
//!   `num_unpaired_faces(α)` counts the faces of `α` **that are in `L(x)`** and
//!   are in neither `V` nor `C` yet; inside a star those are exactly the
//!   codimension-1 faces that still contain `x`, of which a `p`-cell has `p`.
//! - **§3.1, Proposition 5**: *"the maximum number of cells in the lower star of
//!   a voxel in a 3D lattice is 27"* — the full star is 8 cubes, 12 squares, 6
//!   edges and the voxel. `STAR_CELLS` is that 27 and `max_lower_star_cells` is
//!   asserted against it on every row.
//! - **§3.1, Proposition 4**: *"Each cell in `L(x)` will be paired and included
//!   in `V` or inserted into `C`"*. Asserted per voxel and again in aggregate.
//! - **§4, Lemma 10**: *"A critical 2-cell occurs only when `R_i` is the entire
//!   octahedron"*, with the converse proved in the same paragraph. Through the
//!   §4 bijection `φ : L_i \ x_i → R_i`, which drops dimension by one, that is a
//!   statement about critical **3-cells** of the lower star, and it is the
//!   `local_maxima` cross-check.
//! - **§3.1, Eq. (8)**, the tie-break `g'(i,j,k) = g + η(i + Ij + IJk)/(3IJK)`,
//!   quoted **in order to reject it**: it depends on the whole image dimensions
//!   `I, J, K`, so it is chunk-dependent and hash-breaking. It is not used here.
//! - **§4, Theorem 11** and the ordering-independence sentence C1 tests —
//!   *"the results in Section 4 show that for 2D and 3D complexes the number and
//!   type of critical cells found by ProcessLowerStars are independent of this
//!   ordering"* — are the claim under test, not an input to the code.
//!
//! # What this crate owns rather than the paper: the tie-break
//!
//! Algorithm 1 wants distinct voxel values, and this crate's reference fields
//! tie *exactly* — `box_exact` is exactly zero across whole faces of its
//! boundary. Eq. (8) is unusable for the reason above, so the registration fixed
//! a **chunk-local exact** order instead: rank voxels by `(value, linear_index)`
//! lexicographically, with values compared by `f64::total_cmp`. That is a total
//! order on any grid, needs no `η`, and perturbs no sample.
//!
//! It is realised once per grid as a dense `rank: Vec<u32>`, and then **every**
//! comparison in Algorithm 1 — the lower-star membership test and the whole `G`
//! ordering — runs on ranks, which are integers. There is no floating-point
//! comparison anywhere inside the algorithm, so "exact" is structural rather
//! than argued.
//!
//! C1 is the same census a second time under `(value, Reverse(linear_index))`:
//! the same `total_cmp` on values, ties broken by *descending* linear index.
//! Both censuses are on every row, as `critical_0..3` and
//! `critical_0_reverse..3_reverse`, so a reader can see both without rerunning.
//!
//! # Why the priority queue's order is what makes line 16 well-defined
//!
//! Algorithm 1 line 16 is `remove pair(α) from PQzero`, unconditionally. Nothing
//! in the pseudocode says why `pair(α)` must be *in* `PQzero` — and it is worth
//! saying, because the assertion that checks it is the sharpest control here.
//!
//! `G` ranks a cell after its faces: a face's descending rank tuple is a strict
//! subsequence of its coface's, so at the first index where they differ the face
//! carries the smaller rank, and otherwise it is a shorter prefix. `PQone` pops
//! the `G`-minimum. So if `pair(α)` were itself in `PQone`, it would be popped
//! **before** `α`, and on being popped it has zero unpaired faces and is moved
//! to `PQzero` by line 13. A cell reaches one unpaired face only when a face of
//! it is placed, and every such placement — lines 6, 15 and 22 — is immediately
//! followed by a push of the newly-eligible cofaces at lines 8, 17 and 23, so
//! nothing eligible is ever missed. `Star::remove_from_zero` encodes that
//! argument as an assertion; if it fires, the transcription is wrong.
//!
//! # The ambiguous-cell set is rebuilt from the table, not read from an extractor
//!
//! No public API exposes per-cell ambiguity, so C2's population is reconstructed
//! from the grid exactly as `experiment_p17` does: gather the eight corner
//! values in this crate's corner numbering, form the 8-bit case with the crate's
//! own `is_inside`, and index the public `AMBIGUOUS_FACES[case]`. Non-zero means
//! the cell has an ambiguous face. The corner numbering — bit `k` is axis `k` —
//! is **asserted** in `main` from `EDGE_CORNERS` and `EDGE_AXIS` rather than
//! assumed, because getting it wrong silently relabels the whole population.
//!
//! Containment is over the cubical complex, not over voxels: a grid cell
//! contains 27 cells of `K` — 8 vertices, 12 edges, 6 squares and 1 cube — and
//! it *contains a critical cell* when any of those 27 is critical. C2 is
//! containment and deliberately **not** set equality: a Morse census can be
//! non-empty where no Marching Cubes ambiguity exists, so
//! `critical_cells_outside_ambiguous` is a reported excess rather than a
//! failure. Equality would indict the instrument instead of testing the claim.
//!
//! # Which columns decide which clause
//!
//! - **C1** is the conjunction of `census_matches_reverse_order` over all 24
//!   rows. The four `critical_*_reverse` columns say *by how much* if it fails.
//! - **C2** is the conjunction of `ambiguous_containment_holds`, which is
//!   `ambiguous_with_critical == ambiguous_cells`, over all 24 rows.
//! - **C3** is `critical_total` read across `17³/33³/65³`: below 2× change on
//!   `sphere` and `torus`, above 4× growth on `noise_cavity`.
//!
//! `ns_per_voxel` is reported and **gates nothing**. M-348 is the incident where
//! a discovery was demoted for resting on a wall clock; a census is a count.
//!
//! # The controls, and what each would catch
//!
//! - `max_lower_star_cells ≤ 27`, per voxel and per row. Above 27 means the star
//!   enumeration is wrong, not that Proposition 5 is.
//! - Proposition 4, twice: every member of every lower star ends in `V` or `C`,
//!   and `2·pairs + critical_total == cells_total == (2n−1)³`. That closed form
//!   is `n³ + 3n²(n−1) + 3n(n−1)² + (n−1)³` collapsed, and it is the partition
//!   claim of §3.1 stated as an integer identity.
//! - Every pair `(α, β) ∈ V` satisfies `α < β` and `dim β = dim α + 1`, and no
//!   cell is placed twice or placed in both `V` and `C`.
//! - Lemma 10: `critical_3 == local_maxima`, where `local_maxima` counts voxels
//!   whose lower star holds all eight cubes — the reduced lower star being the
//!   entire octahedron. A census that drifts from that has invented a 3-cell.
//! - The critical bitmask rebuilt for containment must have exactly
//!   `critical_total` bits set, which cross-checks the star-code to canonical
//!   anchor conversion the containment test depends on.

mod common;

use std::time::Instant;

use isomesh::Sdf;
use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::table::{AMBIGUOUS_FACES, EDGE_AXIS, EDGE_CORNERS, is_inside};

type Scalar = f64;

/// Samples per axis. C3 needs three grids to read a growth ratio from.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// A star cell is coded `(axes << 3) | dirs`, so codes span `0..64`.
///
/// `axes` is the set of axes the cell spans out of the voxel; a set `dirs` bit
/// says the cell extends in the **negative** direction along that axis. Only
/// codes with `dirs & !axes == 0` are cells, which is 27 of the 64.
const CODES: usize = 64;

/// Cells in the full star of a voxel of a 3-D lattice: `1 + 6 + 12 + 8`.
///
/// Proposition 5 of `10.1109/tpami.2011.95`: *"the maximum number of cells in
/// the lower star of a voxel in a 3D lattice is 27"*.
const STAR_CELLS: usize = 27;

/// Strides of the 3×3×3 neighbourhood, one per axis.
const POW3: [i32; 3] = [1, 3, 9];

/// The centre of the 3×3×3 neighbourhood — the voxel itself.
const SELF_NB: usize = 13;

/// Is this code a cell? `dirs` may only name axes the cell actually spans.
const fn is_cell(code: usize) -> bool {
    let axes = (code >> 3) as u8;
    let dirs = (code & 7) as u8;
    dirs & !axes == 0
}

/// The dimension of star cell `code` — the number of axes it spans.
const fn cell_dim(code: usize) -> usize {
    ((code >> 3) as u8).count_ones() as usize
}

/// Index into the 3×3×3 neighbourhood of `x + Σ_{a ∈ subset} s_a e_a`.
const fn nb_index(subset: u8, dirs: u8) -> u8 {
    let mut idx = SELF_NB as i32;
    let mut k = 0usize;
    while k < 3 {
        if subset & (1 << k) != 0 {
            idx += if dirs & (1 << k) != 0 {
                -POW3[k]
            } else {
                POW3[k]
            };
        }
        k += 1;
    }
    idx as u8
}

/// The non-centre vertices of each star cell, as 3×3×3 neighbourhood indices.
///
/// §2.1: a cell spanning axes `A` out of `x` has vertex set
/// `{ x + Σ_{a ∈ T} s_a e_a : T ⊆ A }`, so its `2^|A| − 1` non-centre vertices
/// are indexed by the non-empty subsets of `A`.
const STAR_VERTS: [[u8; 8]; CODES] = build_star_verts();

const fn build_star_verts() -> [[u8; 8]; CODES] {
    let mut out = [[0u8; 8]; CODES];
    let mut code = 0usize;
    while code < CODES {
        if is_cell(code) {
            let axes = (code >> 3) as u8;
            let dirs = (code & 7) as u8;
            let mut n = 0usize;
            let mut subset = 1u8;
            while subset < 8 {
                if subset & !axes == 0 {
                    out[code][n] = nb_index(subset, dirs);
                    n += 1;
                }
                subset += 1;
            }
        }
        code += 1;
    }
    out
}

/// How many non-centre vertices each star cell has: `2^dim − 1`.
const STAR_NVERTS: [u8; CODES] = build_star_nverts();

const fn build_star_nverts() -> [u8; CODES] {
    let mut out = [0u8; CODES];
    let mut code = 0usize;
    while code < CODES {
        if is_cell(code) {
            out[code] = ((1usize << cell_dim(code)) - 1) as u8;
        }
        code += 1;
    }
    out
}

/// The codimension-1 faces of each star cell that are themselves star cells.
///
/// Dropping one spanned axis and keeping the directions on the rest. The other
/// codimension-1 faces of the cell do not contain `x`, so §3.1's
/// `num_unpaired_faces` — which counts faces *in `L(x)`* — never sees them.
const STAR_FACES: [[u8; 3]; CODES] = build_star_faces();

const fn build_star_faces() -> [[u8; 3]; CODES] {
    let mut out = [[0u8; 3]; CODES];
    let mut code = 0usize;
    while code < CODES {
        if is_cell(code) {
            let axes = (code >> 3) as u8;
            let dirs = (code & 7) as u8;
            let mut n = 0usize;
            let mut k = 0usize;
            while k < 3 {
                if axes & (1 << k) != 0 {
                    let sub_axes = axes & !(1 << k);
                    out[code][n] = (sub_axes << 3) | (dirs & sub_axes);
                    n += 1;
                }
                k += 1;
            }
        }
        code += 1;
    }
    out
}

/// How many codimension-1 faces each star cell has inside the star: its `dim`.
const STAR_NFACES: [u8; CODES] = build_star_nfaces();

const fn build_star_nfaces() -> [u8; CODES] {
    let mut out = [0u8; CODES];
    let mut code = 0usize;
    while code < CODES {
        if is_cell(code) {
            out[code] = cell_dim(code) as u8;
        }
        code += 1;
    }
    out
}

/// Neighbourhood index of each star cell's canonical anchor.
///
/// A cell of `K` has one description independent of which voxel's star names it:
/// the lowest-coordinate vertex plus the set of spanned axes. That anchor is the
/// vertex reached by taking every negative direction, i.e. subset `dirs`.
const STAR_ANCHOR: [u8; CODES] = build_star_anchor();

const fn build_star_anchor() -> [u8; CODES] {
    let mut out = [0u8; CODES];
    let mut code = 0usize;
    while code < CODES {
        if is_cell(code) {
            let dirs = (code & 7) as u8;
            out[code] = nb_index(dirs, dirs);
        }
        code += 1;
    }
    out
}

/// Corner `c`'s offset in this crate's cube numbering — asserted in `main`.
const fn corner_offset(c: u8) -> [usize; 3] {
    [
        (c & 1) as usize,
        ((c >> 1) & 1) as usize,
        ((c >> 2) & 1) as usize,
    ]
}

/// Where a cell of the lower star ended up: `V`, `C`, or neither yet.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Placed {
    /// In neither `V` nor `C`, so it still counts toward `num_unpaired_faces`.
    No,
    /// In `V`, as either member of a pair.
    Paired,
    /// In `C`.
    Critical,
}

/// The census of one grid under one voxel ordering.
struct Census {
    /// Critical cells by dimension.
    critical: [u64; 4],
    /// Pairs in `V`, over every lower star.
    pairs: u64,
    /// The largest `|L(x)|` seen — Proposition 5's 27 is the ceiling.
    max_star: usize,
    /// Voxels whose reduced lower star is the whole octahedron (Lemma 10).
    local_maxima: u64,
    /// Bit `axes` set on the anchor voxel when that cell of `K` is critical.
    mask: Vec<u8>,
}

impl Census {
    /// Critical cells of every dimension.
    fn total(&self) -> u64 {
        self.critical.iter().sum()
    }

    /// How many bits the critical mask carries — must equal `total`.
    fn mask_bits(&self) -> u64 {
        self.mask.iter().map(|&m| u64::from(m.count_ones())).sum()
    }

    /// `c0 − c1 + c2 − c3`, the Euler characteristic the census implies.
    fn euler(&self) -> i64 {
        self.critical[0] as i64 - self.critical[1] as i64 + self.critical[2] as i64
            - self.critical[3] as i64
    }
}

/// Scratch for one lower star, reused across every voxel of every grid.
///
/// Everything is a fixed-size array indexed by star code, so `ProcessLowerStar`
/// allocates nothing: at 65³ it runs 274,625 times per ordering.
struct Star {
    /// Is this code a cell of `L(x)`?
    member: [bool; CODES],
    /// `V`/`C` membership.
    placed: [Placed; CODES],
    /// Is this code currently in `PQzero`?
    in_zero: [bool; CODES],
    /// Has this code ever been pushed to `PQone`? §3.3: *"at most once"*.
    pushed_one: [bool; CODES],
    /// `G(α)`: the cell's vertex ranks in decreasing order.
    g: [[u32; 8]; CODES],
    /// How many entries of `g` are live.
    glen: [u8; CODES],
    /// `PQzero`, unordered; `pop_front` scans for the `G`-minimum.
    pq_zero: [u8; STAR_CELLS],
    /// Live length of `pq_zero`.
    nz: usize,
    /// `PQone`, unordered; `pop_front` scans for the `G`-minimum.
    pq_one: [u8; STAR_CELLS],
    /// Live length of `pq_one`.
    no: usize,
    /// Linear voxel index of each 3×3×3 neighbour, `u32::MAX` off-grid.
    nb_lin: [u32; 27],
    /// Rank of each 3×3×3 neighbour, `u32::MAX` off-grid so it is never lower.
    nb_rank: [u32; 27],
}

impl Star {
    /// A star with nothing in it.
    fn new() -> Self {
        Self {
            member: [false; CODES],
            placed: [Placed::No; CODES],
            in_zero: [false; CODES],
            pushed_one: [false; CODES],
            g: [[0u32; 8]; CODES],
            glen: [0u8; CODES],
            pq_zero: [0u8; STAR_CELLS],
            nz: 0,
            pq_one: [0u8; STAR_CELLS],
            no: 0,
            nb_lin: [u32::MAX; 27],
            nb_rank: [u32::MAX; 27],
        }
    }

    /// Load the 3×3×3 neighbourhood of voxel `at` and build `L(x)`.
    ///
    /// Eq. (7): a star cell is in the lower star exactly when every one of its
    /// other vertices is in the grid and ranks below `x`. Off-grid neighbours
    /// carry `u32::MAX`, which is never below a real rank, so the two conditions
    /// collapse into one comparison.
    fn load(&mut self, rank: &[u32], n: u32, at: [u32; 3]) -> usize {
        let nu = n as usize;
        let here = at[0] as usize + nu * at[1] as usize + nu * nu * at[2] as usize;
        let rx = rank[here];

        for dz in -1i32..=1 {
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let slot =
                        (SELF_NB as i32 + dx * POW3[0] + dy * POW3[1] + dz * POW3[2]) as usize;
                    let p = [at[0] as i32 + dx, at[1] as i32 + dy, at[2] as i32 + dz];
                    if p.iter().any(|&c| c < 0 || c >= n as i32) {
                        self.nb_lin[slot] = u32::MAX;
                        self.nb_rank[slot] = u32::MAX;
                    } else {
                        let lin = p[0] as usize + nu * p[1] as usize + nu * nu * p[2] as usize;
                        self.nb_lin[slot] = lin as u32;
                        self.nb_rank[slot] = rank[lin];
                    }
                }
            }
        }

        let mut star = 0usize;
        for code in 0..CODES {
            self.placed[code] = Placed::No;
            self.in_zero[code] = false;
            self.pushed_one[code] = false;
            self.member[code] = false;
            if !is_cell(code) {
                continue;
            }
            let nv = STAR_NVERTS[code] as usize;
            let lower = STAR_VERTS[code][..nv]
                .iter()
                .all(|&v| self.nb_rank[v as usize] < rx);
            if !lower {
                continue;
            }
            self.member[code] = true;
            star += 1;
            // §3.1: `G(α) = (g(x), g(y_1), …, g(y_k))` with the values in
            // decreasing order, compared lexicographically. On ranks.
            let mut buf = [0u32; 8];
            buf[0] = rx;
            for (slot, &v) in buf[1..=nv].iter_mut().zip(&STAR_VERTS[code][..nv]) {
                *slot = self.nb_rank[v as usize];
            }
            buf[..=nv].sort_unstable_by(|a, b| b.cmp(a));
            self.g[code] = buf;
            self.glen[code] = (nv + 1) as u8;
        }
        self.nz = 0;
        self.no = 0;
        star
    }

    /// Is `G(a)` below `G(b)`? Lexicographic, shorter-prefix-first.
    fn g_less(&self, a: usize, b: usize) -> bool {
        self.g[a][..self.glen[a] as usize] < self.g[b][..self.glen[b] as usize]
    }

    /// §3.1's `num_unpaired_faces(α)`: faces of `α` in `L(x)` not yet in `V`/`C`.
    fn unpaired_faces(&self, code: usize) -> usize {
        STAR_FACES[code][..STAR_NFACES[code] as usize]
            .iter()
            .filter(|&&f| self.member[f as usize] && self.placed[f as usize] == Placed::No)
            .count()
    }

    /// §3.1's `pair(α)`: the single available face, when there is exactly one.
    fn sole_unpaired_face(&self, code: usize) -> usize {
        let mut found = usize::MAX;
        for &f in &STAR_FACES[code][..STAR_NFACES[code] as usize] {
            let f = f as usize;
            if self.member[f] && self.placed[f] == Placed::No {
                assert!(found == usize::MAX, "pair(alpha) is not unique");
                found = f;
            }
        }
        assert!(found != usize::MAX, "pair(alpha) does not exist");
        found
    }

    /// Push to `PQzero`, once.
    fn push_zero(&mut self, code: usize) {
        assert!(!self.in_zero[code], "a cell entered PQzero twice");
        assert!(self.nz < STAR_CELLS, "PQzero exceeded the 27-cell star");
        self.pq_zero[self.nz] = code as u8;
        self.nz += 1;
        self.in_zero[code] = true;
    }

    /// Push to `PQone`, once — §3.3: *"each cell is inserted in PQone at most
    /// once"*.
    fn push_one(&mut self, code: usize) {
        assert!(self.no < STAR_CELLS, "PQone exceeded the 27-cell star");
        self.pq_one[self.no] = code as u8;
        self.no += 1;
        self.pushed_one[code] = true;
    }

    /// `PQzero.pop_front`, i.e. the `G`-minimum.
    fn pop_zero(&mut self) -> usize {
        let mut best = 0usize;
        for i in 1..self.nz {
            if self.g_less(self.pq_zero[i] as usize, self.pq_zero[best] as usize) {
                best = i;
            }
        }
        let code = self.pq_zero[best] as usize;
        self.nz -= 1;
        self.pq_zero[best] = self.pq_zero[self.nz];
        self.in_zero[code] = false;
        code
    }

    /// `PQone.pop_front`, i.e. the `G`-minimum.
    fn pop_one(&mut self) -> usize {
        let mut best = 0usize;
        for i in 1..self.no {
            if self.g_less(self.pq_one[i] as usize, self.pq_one[best] as usize) {
                best = i;
            }
        }
        let code = self.pq_one[best] as usize;
        self.no -= 1;
        self.pq_one[best] = self.pq_one[self.no];
        code
    }

    /// Drop `code` out of `PQzero` — Algorithm 1 line 16.
    ///
    /// The precondition is asserted, not assumed: see the module header on why
    /// the `G` order forces `pair(α)` to be in `PQzero` and nowhere else.
    fn remove_from_zero(&mut self, code: usize) {
        assert!(
            self.in_zero[code],
            "Algorithm 1 line 16: pair(alpha) was not in PQzero, so the G \
             ordering is not ranking faces before cofaces"
        );
        let at = (0..self.nz)
            .find(|&i| self.pq_zero[i] as usize == code)
            .expect("in_zero and pq_zero disagree");
        self.nz -= 1;
        self.pq_zero[at] = self.pq_zero[self.nz];
        self.in_zero[code] = false;
    }

    /// `β > α` inside one star: `α`'s axes are a proper subset, directions agree.
    fn is_coface(beta: usize, alpha: usize) -> bool {
        let (ab, db) = ((beta >> 3) as u8, (beta & 7) as u8);
        let (aa, da) = ((alpha >> 3) as u8, (alpha & 7) as u8);
        aa & !ab == 0 && aa != ab && db & aa == da
    }

    /// Record `V[lower] := upper`, with the audit Forman's definition wants.
    fn pair(&mut self, lower: usize, upper: usize) {
        assert!(
            self.placed[lower] == Placed::No && self.placed[upper] == Placed::No,
            "a cell was placed twice"
        );
        assert!(
            Self::is_coface(upper, lower) && cell_dim(upper) == cell_dim(lower) + 1,
            "V must pair alpha < beta with dim beta = dim alpha + 1"
        );
        assert!(
            self.member[lower] && self.member[upper],
            "V must pair cells of the same lower star"
        );
        self.placed[lower] = Placed::Paired;
        self.placed[upper] = Placed::Paired;
    }

    /// Push every star cell that a placement has just made eligible.
    ///
    /// Algorithm 1 lines 8, 17 and 23 are the same statement with a different
    /// set of triggers, so they are the same call with one or two of them.
    fn push_newly_eligible(&mut self, triggers: [usize; 2]) {
        for code in 0..CODES {
            if !self.member[code] || self.placed[code] != Placed::No || self.pushed_one[code] {
                continue;
            }
            let over = triggers
                .iter()
                .any(|&t| t != usize::MAX && Self::is_coface(code, t));
            if over && self.unpaired_faces(code) == 1 {
                self.push_one(code);
            }
        }
    }

    /// Algorithm 1, lines 2–26, for one voxel.
    ///
    /// Returns `(pairs, critical_by_dimension)` for this lower star.
    fn process(&mut self, star: usize) -> (u64, [u64; 4]) {
        let mut pairs = 0u64;
        let mut critical = [0u64; 4];

        // Line 2: `if L(x) = {x}` — code 0 is `x` itself.
        if star == 1 {
            // Line 3.
            self.placed[0] = Placed::Critical;
            critical[0] += 1;
            return (pairs, critical);
        }

        // Line 5: δ is the G-minimal 1-cell of L(x).
        let mut delta = usize::MAX;
        for code in 0..CODES {
            let one_cell = self.member[code] && cell_dim(code) == 1;
            if one_cell && (delta == usize::MAX || self.g_less(code, delta)) {
                delta = code;
            }
        }
        assert!(
            delta != usize::MAX,
            "a lower star with more than one cell has no 1-cell, which \
             contradicts §3.1: faces of a star cell are star cells"
        );

        // Line 6: V[x] := δ.
        self.pair(0, delta);
        pairs += 1;

        // Line 7: all other 1-cells to PQzero.
        for code in 0..CODES {
            if self.member[code] && cell_dim(code) == 1 && code != delta {
                self.push_zero(code);
            }
        }

        // Line 8.
        self.push_newly_eligible([delta, usize::MAX]);

        // Line 9.
        while self.no > 0 || self.nz > 0 {
            // Line 10.
            while self.no > 0 {
                // Line 11.
                let alpha = self.pop_one();
                if self.unpaired_faces(alpha) == 0 {
                    // Lines 12–13.
                    self.push_zero(alpha);
                } else {
                    // Lines 15–17.
                    let mate = self.sole_unpaired_face(alpha);
                    self.pair(mate, alpha);
                    pairs += 1;
                    self.remove_from_zero(mate);
                    self.push_newly_eligible([alpha, mate]);
                }
            }
            // Lines 20–23.
            if self.nz > 0 {
                let gamma = self.pop_zero();
                self.placed[gamma] = Placed::Critical;
                critical[cell_dim(gamma)] += 1;
                self.push_newly_eligible([gamma, usize::MAX]);
            }
        }

        (pairs, critical)
    }
}

/// Dense ranks of the voxels under the registered order, or its reverse.
///
/// `(value, linear_index)` lexicographically with `f64::total_cmp` on values —
/// a total order that perturbs nothing. `reverse` breaks ties by *descending*
/// linear index instead, which is C1's second ordering.
fn ranks(values: &[Scalar], reverse: bool) -> Vec<u32> {
    let mut order: Vec<u32> = (0..values.len() as u32).collect();
    order.sort_unstable_by(|&a, &b| {
        values[a as usize]
            .total_cmp(&values[b as usize])
            .then_with(|| if reverse { b.cmp(&a) } else { a.cmp(&b) })
    });
    let mut rank = vec![0u32; values.len()];
    for (r, &v) in order.iter().enumerate() {
        rank[v as usize] = r as u32;
    }
    rank
}

/// `ProcessLowerStars` over the whole grid, under one ordering.
fn process_lower_stars(rank: &[u32], n: u32, star: &mut Star) -> Census {
    let nu = n as usize;
    let mut out = Census {
        critical: [0; 4],
        pairs: 0,
        max_star: 0,
        local_maxima: 0,
        mask: vec![0u8; nu * nu * nu],
    };

    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let size = star.load(rank, n, [x, y, z]);
                assert!(
                    size <= STAR_CELLS,
                    "|L(x)| = {size} exceeds Proposition 5's bound of {STAR_CELLS}"
                );
                out.max_star = out.max_star.max(size);

                // Lemma 10 through φ: a critical 3-cell of the lower star needs
                // the reduced lower star to be the entire octahedron, i.e. all
                // eight cubes of the full star present.
                if (0..8).all(|d| star.member[(7 << 3) | d]) {
                    out.local_maxima += 1;
                }

                let (pairs, critical) = star.process(size);
                out.pairs += pairs;
                for (slot, add) in out.critical.iter_mut().zip(critical) {
                    *slot += add;
                }

                // Proposition 4, per star.
                let mut placed_cells = 0usize;
                for (code, ((&member, &placed), &anchor_nb)) in star
                    .member
                    .iter()
                    .zip(&star.placed)
                    .zip(&STAR_ANCHOR)
                    .enumerate()
                {
                    if !member {
                        continue;
                    }
                    assert!(
                        placed != Placed::No,
                        "Proposition 4: a cell of L(x) reached neither V nor C"
                    );
                    placed_cells += 1;
                    if placed == Placed::Critical {
                        let anchor = star.nb_lin[anchor_nb as usize] as usize;
                        out.mask[anchor] |= 1 << (code >> 3);
                    }
                }
                assert_eq!(
                    placed_cells, size,
                    "the member count moved during processing"
                );
                assert_eq!(
                    2 * pairs + critical.iter().sum::<u64>(),
                    size as u64,
                    "the star's cells are not exactly its pairs plus its criticals"
                );
            }
        }
    }
    out
}

/// The ambiguous-cell population, and C2's containment counts.
///
/// Returns `(ambiguous_cells, ambiguous_with_critical, critical_outside)`.
fn containment(values: &[Scalar], n: u32, critical: &[u8]) -> (u64, u64, u64) {
    let nu = n as usize;
    let cells = nu - 1;
    let mut covered = vec![0u8; nu * nu * nu];
    let mut ambiguous = 0u64;
    let mut with_critical = 0u64;

    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let mut case = 0u8;
                for c in 0..8u8 {
                    let o = corner_offset(c);
                    let lin = (x + o[0]) + nu * (y + o[1]) + nu * nu * (z + o[2]);
                    if is_inside(values[lin]) {
                        case |= 1 << c;
                    }
                }
                if AMBIGUOUS_FACES[case as usize] == 0 {
                    continue;
                }
                ambiguous += 1;
                // The 27 cells of `K` inside this grid cell: choose which axes
                // the sub-cell spans, and for each remaining axis whether it
                // sits at the low or the high face.
                let mut hit = false;
                for axes in 0..8usize {
                    for high in 0..8usize {
                        if high & axes != 0 {
                            continue;
                        }
                        let lin = (x + (high & 1))
                            + nu * (y + ((high >> 1) & 1))
                            + nu * nu * (z + ((high >> 2) & 1));
                        covered[lin] |= 1 << axes;
                        hit |= critical[lin] & (1 << axes) != 0;
                    }
                }
                if hit {
                    with_critical += 1;
                }
            }
        }
    }

    let outside = critical
        .iter()
        .zip(&covered)
        .map(|(&c, &v)| u64::from((c & !v).count_ones()))
        .sum();
    (ambiguous, with_critical, outside)
}

/// The corner numbering and star coding this experiment depends on, checked.
fn assert_conventions() {
    // Bit `k` is axis `k`, from the crate's own edge tables rather than assumed.
    for (e, corners) in EDGE_CORNERS.iter().enumerate() {
        assert_eq!(
            corners[0] ^ corners[1],
            1 << EDGE_AXIS[e],
            "edge {e} does not join corners differing in the bit of its own axis"
        );
    }
    for c in 0..8u8 {
        for (k, &v) in corner_offset(c).iter().enumerate() {
            assert_eq!(
                v,
                usize::from((c >> k) & 1),
                "corner bit {k} is not axis {k}"
            );
        }
    }
    // Proposition 5's 27, and the full star's census by dimension.
    let live = (0..CODES).filter(|&c| is_cell(c)).count();
    assert_eq!(live, STAR_CELLS, "the star coding does not have 27 cells");
    for (d, want) in [1usize, 6, 12, 8].into_iter().enumerate() {
        let got = (0..CODES)
            .filter(|&c| is_cell(c) && cell_dim(c) == d)
            .count();
        assert_eq!(
            got, want,
            "the full star should have {want} cells of dimension {d}"
        );
    }
    // A face of a star cell is a star cell of one lower dimension, and the
    // coface test agrees with it. This is §2.1's `α < β` in both directions.
    for code in 0..CODES {
        if !is_cell(code) {
            continue;
        }
        for &f in &STAR_FACES[code][..STAR_NFACES[code] as usize] {
            let f = f as usize;
            assert!(is_cell(f), "face {f} of {code} is not a star cell");
            assert_eq!(cell_dim(f) + 1, cell_dim(code), "face dimension is wrong");
            assert!(Star::is_coface(code, f), "the coface test rejects a face");
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    assert_conventions();

    let prereg = isomesh::experiment!("P-58");
    common::experiment::run(prereg, |run| {
        println!(
            "{:<15} {:>4} {:>7} {:>7} {:>7} {:>6} {:>8} {:>5} {:>6} {:>8} {:>9} {:>8} {:>9}",
            "field",
            "n",
            "crit0",
            "crit1",
            "crit2",
            "crit3",
            "total",
            "max*",
            "rev?",
            "ambig",
            "amb+crit",
            "outside",
            "ns/voxel"
        );

        let mut star = Star::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            // Inline block, so no `return` in here (M-199 / M-253).
            let (lo, hi) = field.domain();
            for &n in &RESOLUTIONS {
                let nu = n as usize;
                let h = (hi[0] - lo[0]) / f64::from(n - 1);
                let voxels = nu * nu * nu;

                let mut values = Vec::with_capacity(voxels);
                for z in 0..nu {
                    for y in 0..nu {
                        for x in 0..nu {
                            values.push(field.sample([
                                lo[0] + h * x as f64,
                                lo[1] + h * y as f64,
                                lo[2] + h * z as f64,
                            ]));
                        }
                    }
                }

                // How much work the tie-break can possibly be doing. If a grid
                // has no ties the two orderings are the *same* ordering and C1
                // is vacuous on that row, so the count belongs on the row.
                let mut sorted = values.clone();
                sorted.sort_unstable_by(|a, b| a.total_cmp(b));
                let distinct = 1 + sorted
                    .windows(2)
                    .filter(|w| w[0].total_cmp(&w[1]).is_ne())
                    .count();
                let tied = voxels - distinct;

                let forward = ranks(&values, false);
                let reverse = ranks(&values, true);

                let started = Instant::now();
                let fwd = process_lower_stars(&forward, n, &mut star);
                let elapsed = started.elapsed();
                let rev = process_lower_stars(&reverse, n, &mut star);

                // §3.1: the lower stars partition K, and Proposition 4 places
                // every cell. `(2n−1)³` is the collapsed closed form of
                // `n³ + 3n²(n−1) + 3n(n−1)² + (n−1)³`.
                let cells_total = (2 * nu - 1).pow(3) as u64;
                for census in [&fwd, &rev] {
                    assert_eq!(
                        2 * census.pairs + census.total(),
                        cells_total,
                        "{name} at {n}³: V and C do not exhaust the cubical complex"
                    );
                    assert!(
                        census.max_star <= STAR_CELLS,
                        "{name} at {n}³: max lower star {} exceeds 27",
                        census.max_star
                    );
                    // Lemma 10, both arms: a critical 3-cell exactly where the
                    // reduced lower star is the whole octahedron.
                    assert_eq!(
                        census.critical[3], census.local_maxima,
                        "{name} at {n}³: Lemma 10 broken — {} critical 3-cells \
                         against {} full-octahedron stars",
                        census.critical[3], census.local_maxima
                    );
                    assert_eq!(
                        census.mask_bits(),
                        census.total(),
                        "{name} at {n}³: the critical bitmask lost cells, so the \
                         star-code to anchor conversion is wrong"
                    );
                    // Theorem 3 with Theorem 6: the Morse chain complex
                    // computes `H_*(K)`. `K` here is the full cubical complex of
                    // a box, which is contractible, so its Euler
                    // characteristic is 1 and the alternating sum of the
                    // critical counts has to be 1 as well. Nothing in the
                    // transcription uses this, so it is a global check on the
                    // whole pairing rather than a restatement of one.
                    assert_eq!(
                        census.euler(),
                        1,
                        "{name} at {n}³: χ = {} from ({}, {}, {}, {}), but the \
                         cubical complex of a box is contractible, so Theorem 6 \
                         demands χ = 1",
                        census.euler(),
                        census.critical[0],
                        census.critical[1],
                        census.critical[2],
                        census.critical[3],
                    );
                }

                let (ambiguous, with_critical, outside) = containment(&values, n, &fwd.mask);
                let matches = fwd.critical == rev.critical;
                let holds = with_critical == ambiguous;
                let ns_per_voxel = elapsed.as_secs_f64() * 1e9 / voxels as f64;

                println!(
                    "{name:<15} {n:>4} {:>7} {:>7} {:>7} {:>6} {:>8} {:>5} {:>6} {:>8} {:>9} \
                     {:>8} {ns_per_voxel:>9.1}",
                    fwd.critical[0],
                    fwd.critical[1],
                    fwd.critical[2],
                    fwd.critical[3],
                    fwd.total(),
                    fwd.max_star,
                    matches,
                    ambiguous,
                    with_critical,
                    outside,
                );

                run.record(&[
                    ("field", name.to_string()),
                    ("samples_per_axis", n.to_string()),
                    ("voxels", voxels.to_string()),
                    ("critical_0", fwd.critical[0].to_string()),
                    ("critical_1", fwd.critical[1].to_string()),
                    ("critical_2", fwd.critical[2].to_string()),
                    ("critical_3", fwd.critical[3].to_string()),
                    ("critical_total", fwd.total().to_string()),
                    ("max_lower_star_cells", fwd.max_star.to_string()),
                    ("census_matches_reverse_order", matches.to_string()),
                    ("ambiguous_cells", ambiguous.to_string()),
                    ("ambiguous_with_critical", with_critical.to_string()),
                    ("ambiguous_containment_holds", holds.to_string()),
                    ("critical_cells_outside_ambiguous", outside.to_string()),
                    ("ns_per_voxel", format!("{ns_per_voxel:.3}")),
                    ("critical_0_reverse", rev.critical[0].to_string()),
                    ("critical_1_reverse", rev.critical[1].to_string()),
                    ("critical_2_reverse", rev.critical[2].to_string()),
                    ("critical_3_reverse", rev.critical[3].to_string()),
                    ("critical_total_reverse", rev.total().to_string()),
                    ("cells_total", cells_total.to_string()),
                    ("pairs", fwd.pairs.to_string()),
                    ("local_maxima", fwd.local_maxima.to_string()),
                    ("grid_cells", (nu - 1).pow(3).to_string()),
                    ("distinct_values", distinct.to_string()),
                    ("tied_voxels", tied.to_string()),
                    ("euler_characteristic", fwd.euler().to_string()),
                    ("euler_characteristic_reverse", rev.euler().to_string()),
                    (
                        "algorithm_source",
                        String::from(
                            "10.1109/tpami.2011.95 Alg.1 §3.1 + Eq.(7) + Prop.4-5 + \
                             Lemma 10 + Thm 11; Eq.(8) deliberately unused",
                        ),
                    ),
                ]);
            }
        });
    });
}
