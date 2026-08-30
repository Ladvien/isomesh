//! Sample lattices, their generated case tables, and the BCC box spline.
//!
//! Ticket: R-162, which owns this module. Consumed unchanged by R-163 (FCC
//! against BCC, a null registered on purpose) and R-164 (the box-spline half of
//! the same question). A second copy of any of this would be two paths to one
//! answer.
//!
//! # The quantity being predicted
//!
//! A sampling lattice `L` is judged by its **dimensionless second moment**
//!
//! ```text
//!            1        1
//!   G(L) = ------ · ------- · ∫    ‖x‖² dx
//!          V^(5/3)     3      V(L)
//! ```
//!
//! where `V(L)` is the Voronoi cell and `V = vol V(L)`. `G` is invariant under
//! scaling and rotation, so it ranks lattices independently of how finely they
//! are sampled — which is exactly what makes a *matched point density*
//! comparison meaningful. Barnes & Sloane (`10.1137/0604005`) prove `A₃*`, the
//! body-centred cubic lattice, minimises `G` **among three-dimensional
//! lattices**; the bracketing bounds over all quantisers remain unproven, and
//! this module claims nothing about them.
//!
//! Lower `G` is better. That fact is the whole reason [`Lattice::gain_db_over`]
//! reads the way it does, and the sign trap is documented on that method rather
//! than left for a bench author to trip over.
//!
//! # Why the case tables are generated
//!
//! `crates/isomesh/src/marching_cubes/table.rs:1-84` states the repository rule:
//! a case table is *constructed*, never typed from a paper, because a wrong
//! entry produces a mesh that looks fine and is silently non-manifold. This
//! module obeys the same rule for the same reason, and then goes one step
//! further — [`case_table`] for the cubic cell is checked, entry by entry,
//! against the shipped `isomesh::marching_cubes::table::CASES`. That check is
//! the calibration: it says the enumeration here understands the cube before
//! anything is believed about the tetrahedron.
//!
//! # What is deliberately absent
//!
//! There is no FCC reconstruction filter. FCC's Delaunay complex is the
//! tetrahedral-octahedral honeycomb — two cell shapes, not one — and the
//! reconstruction question this phase registered is BCC's
//! (`10.1109/tvcg.2007.70429`). R-163 needs only a zero-set point sample from
//! each lattice, which [`zero_set_hausdorff`] consumes without caring how it was
//! produced.

use isomesh::marching_cubes::table::{
    CASES, EDGE_CORNERS, corner_inside, edge_index, face_corners,
};

/// Corners of a cube.
const CUBE_CORNERS: usize = 8;
/// Edges of a cube.
const CUBE_EDGES: usize = 12;
/// Sign configurations of a cube's corners.
const CUBE_CASES: usize = 1 << CUBE_CORNERS;
/// Corners of a tetrahedron — the natural cell of both non-cubic lattices here.
const SIMPLEX_CORNERS: usize = 4;
/// Sign configurations of a tetrahedron's corners.
const SIMPLEX_CASES: usize = 1 << SIMPLEX_CORNERS;
/// Marker for "no incident segment yet" in the cycle graph.
const NO_LINK: u8 = u8::MAX;

/// Approximation order of the linear BCC box spline: **2**, the same as the
/// trilinear's.
///
/// A box spline with `s` direction vectors in `d` dimensions has polynomial
/// degree `s − d`; here `s = 4`, `d = 3`, so the spline is piecewise **linear**
/// and reproduces polynomials up to degree 1. Approximation order is degree plus
/// one. This is `P-164` C1's prediction, derived rather than measured, and the
/// bench measures it independently.
pub(crate) const BCC_BOX_SPLINE_ORDER: usize = 2;

/// Lattice sites carrying a non-zero weight at a generic point: **4**.
///
/// The linear box spline is the simplex ("Courant") element of the BCC Delaunay
/// complex, so at a point interior to one Delaunay tetrahedron exactly that
/// tetrahedron's four vertices contribute. Its support is a rhombic
/// dodecahedron of volume 16 in the integer lattice coordinates of
/// [`bcc_box_spline`], which is **four** fundamental cells — against the
/// trilinear's eight. Same order, half the footprint: that is the whole content
/// of `10.1109/tvcg.2007.70429` for the linear case.
pub(crate) const BCC_BOX_SPLINE_STENCIL: usize = 4;

/// Lattice sites carrying a non-zero weight in a trilinear evaluation: **8**.
pub(crate) const TRILINEAR_STENCIL: usize = 8;

/// The three 3D sample lattices this phase compares.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Lattice {
    /// The integer lattice `Z³` — what every extractor in this crate samples on.
    Cubic,
    /// Body-centred cubic, `A₃*`: the optimal three-dimensional lattice
    /// quantiser.
    Bcc,
    /// Face-centred cubic, `D₃`: the optimal three-dimensional *sphere packing*,
    /// and — this is the point of R-163 — very nearly as good a quantiser as
    /// `A₃*`.
    Fcc,
}

impl Lattice {
    /// All three, in the order every table and CSV in this phase uses.
    pub(crate) const ALL: [Lattice; 3] = [Lattice::Cubic, Lattice::Bcc, Lattice::Fcc];

    /// The lattice's name in the crystallographic notation the registrations use.
    pub(crate) fn name(self) -> &'static str {
        match self {
            Lattice::Cubic => "Z3",
            Lattice::Bcc => "A3*",
            Lattice::Fcc => "D3",
        }
    }

    /// The dimensionless second moment `G`.
    ///
    /// Closed forms, each of which the literal below is asserted against to
    /// `1e-9` on every call — see [`check_g_literals`]:
    ///
    /// | lattice | closed form | value |
    /// |---|---|---|
    /// | `Z³` | `1/12` | `0.083333333…` |
    /// | `A₃*` | `19 / (192·∛2)` | `0.078543281…` |
    /// | `D₃` | `2^(−11/3)` = `1 / (8·∛4)` | `0.078745066…` |
    ///
    /// `Z³` is returned as `1.0 / 12.0` rather than the registration's rounded
    /// `0.0833333`, because the registration writes `1/12 = 0.0833333` and the
    /// fraction is the claim. The two differ in the eighth decimal and by
    /// about `2·10⁻⁷` dB, far under anything this phase measures.
    pub(crate) fn g(self) -> f64 {
        check_g_literals();
        match self {
            Lattice::Cubic => 1.0 / 12.0,
            Lattice::Bcc => 0.078_543_281,
            Lattice::Fcc => 0.078_745_066,
        }
    }

    /// `10 · log10(G(self) / G(other))` — the quantisation gain in dB of `self`
    /// over `other`.
    ///
    /// # Read the sign before quoting a number
    ///
    /// Lower `G` is better, so this ratio is **positive when `self` is the worse
    /// lattice** and the dB figure is the gain *available by moving from `self`
    /// to `other`*. That is the direction both headline predictions are written
    /// in:
    ///
    /// - `Lattice::Cubic.gain_db_over(Lattice::Bcc)` = `+0.2571` dB — `P-162`'s
    ///   `predicted_gain_db`, the 5.748% MSE reduction of BCC over cubic.
    /// - `Lattice::Fcc.gain_db_over(Lattice::Bcc)` = `+0.0111` dB — `P-163`'s
    ///   `predicted_gap_db`, BCC's slender lead over FCC, about 4.3% of the
    ///   cubic gap.
    ///
    /// The reversed call `Lattice::Bcc.gain_db_over(Lattice::Fcc)` is
    /// `−0.0111` dB and is *not* a bug: BCC cannot gain over the lattice it
    /// already beats. A bench that wants a positive magnitude must name the
    /// worse lattice as `self`, or take `abs()` and say so in its column.
    pub(crate) fn gain_db_over(self, other: Lattice) -> f64 {
        10.0 * (self.g() / other.g()).log10()
    }

    /// Generator matrix, rows are the basis vectors, scaled so the fundamental
    /// cell has unit volume.
    ///
    /// Unit determinant is what makes [`lattice_grid`] able to solve for a
    /// point density with one cube root and no per-lattice special case: a
    /// lattice whose fundamental cell has volume `1` and which is then scaled by
    /// `s` holds exactly `1/s³` points per unit volume, whatever its shape.
    /// Every geometric difference between the three lattices survives in the
    /// *shape* of the rows; only the density is normalised away.
    ///
    /// The integer bases before normalisation, and the divisor applied to each:
    ///
    /// | lattice | integer basis rows | `det` | divisor |
    /// |---|---|---|---|
    /// | `Z³` | `(1,0,0) (0,1,0) (0,0,1)` | 1 | 1 |
    /// | `A₃*` | `(2,0,0) (0,2,0) (1,1,1)` | 4 | `∛4` |
    /// | `D₃` | `(1,1,0) (0,1,1) (1,0,1)` | 2 | `∛2` |
    ///
    /// The BCC basis generates `{k ∈ Z³ : k₀ ≡ k₁ ≡ k₂ (mod 2)}` — a cube of
    /// side 2 with its body centre — which is the integer coordinate system
    /// [`bcc_box_spline`] is written in. The FCC basis generates
    /// `{k ∈ Z³ : k₀+k₁+k₂ even}`. Row order is chosen so every determinant is
    /// **positive**, so the sign of a solved lattice coordinate is never a
    /// surprise.
    pub(crate) fn generator(self) -> [[f64; 3]; 3] {
        match self {
            Lattice::Cubic => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            Lattice::Bcc => {
                let d = 4f64.cbrt();
                [
                    [2.0 / d, 0.0, 0.0],
                    [0.0, 2.0 / d, 0.0],
                    [1.0 / d, 1.0 / d, 1.0 / d],
                ]
            }
            Lattice::Fcc => {
                let d = 2f64.cbrt();
                [
                    [1.0 / d, 1.0 / d, 0.0],
                    [0.0, 1.0 / d, 1.0 / d],
                    [1.0 / d, 0.0, 1.0 / d],
                ]
            }
        }
    }

    /// Points per unit volume for the **unscaled** conventional description —
    /// the crystallographic cubic cell of side 1.
    ///
    /// `Z³` has one point per unit cube, BCC two (corner plus body centre), FCC
    /// four (corner plus three face centres). This is the constant a reader
    /// expects to see beside a lattice's name and the one [`lattice_grid`]
    /// reports its density against in prose; it is **not** what solves the
    /// scale, because [`Lattice::generator`] has already normalised the
    /// fundamental cell to unit volume. Keeping both numbers visible is
    /// deliberate: `points_per_cell` is how the lattice is described in the
    /// literature, unit determinant is how it is *measured* here, and confusing
    /// the two is how a matched-density comparison turns into a resolution
    /// change wearing a lattice's name.
    pub(crate) fn points_per_cell(self) -> f64 {
        match self {
            Lattice::Cubic => 1.0,
            Lattice::Bcc => 2.0,
            Lattice::Fcc => 4.0,
        }
    }
}

/// Assert every `G` literal against its closed form.
///
/// `19 / (192·∛2)` and `2^(−11/3)` are computed here and compared to the decimals
/// [`Lattice::g`] returns. The literals are what the registrations quote and what
/// the CSVs will carry, so they are the values the code must use; this is the
/// guard that stops a transposed digit from becoming a published prediction.
/// Measured residuals: `2.17·10⁻¹⁰` for `A₃*` and `3.82·10⁻¹⁰` for `D₃`.
fn check_g_literals() {
    let bcc_closed = 19.0 / (192.0 * 2f64.cbrt());
    let fcc_closed = 2f64.powf(-11.0 / 3.0);
    assert!(
        (0.078_543_281 - bcc_closed).abs() < 1e-9,
        "G(A3*) literal disagrees with 19/(192*cbrt(2)): {bcc_closed}"
    );
    assert!(
        (0.078_745_066 - fcc_closed).abs() < 1e-9,
        "G(D3) literal disagrees with 2^(-11/3): {fcc_closed}"
    );
}

/// A concrete set of lattice sample sites inside a box, at a chosen point
/// density.
#[derive(Clone, Debug)]
pub(crate) struct LatticeGrid {
    /// Which lattice these sites came from.
    pub(crate) lattice: Lattice,
    /// The sites, in world coordinates, ordered lexicographically by their
    /// integer coordinate in [`Lattice::generator`]'s basis.
    pub(crate) sites: Vec<[f64; 3]>,
    /// The factor [`Lattice::generator`]'s unit-volume rows were multiplied by.
    /// Point density is `1 / scale³` per unit volume, for all three lattices.
    pub(crate) scale: f64,
    /// Low corner of the box the sites were clipped to.
    pub(crate) lo: [f64; 3],
    /// High corner of the box the sites were clipped to.
    pub(crate) hi: [f64; 3],
}

impl LatticeGrid {
    /// The box centre, which is always a lattice site.
    ///
    /// The lattice is anchored on the centre rather than on `lo` so that no arm
    /// of the comparison gets a corner-aligned advantage: all three lattices see
    /// the same window onto an infinite periodic point set, and the field being
    /// sampled is centred in that window.
    fn centre(&self) -> [f64; 3] {
        [
            0.5 * (self.lo[0] + self.hi[0]),
            0.5 * (self.lo[1] + self.hi[1]),
            0.5 * (self.lo[2] + self.hi[2]),
        ]
    }

    /// The integer coordinate, in [`Lattice::generator`]'s basis, of a world
    /// point that is known to be a site of this grid.
    fn index_of(&self, inv_scaled: &[[f64; 3]; 3], site: [f64; 3]) -> [i64; 3] {
        let c = self.centre();
        let d = [site[0] - c[0], site[1] - c[1], site[2] - c[2]];
        let v = row_times_matrix(d, inv_scaled);
        [
            v[0].round() as i64,
            v[1].round() as i64,
            v[2].round() as i64,
        ]
    }

    /// Position of the site with integer coordinate `ijk`, or `None` if that
    /// site was clipped away by the box.
    ///
    /// `sites` is built in lexicographic order of `ijk`, so this is a binary
    /// search rather than a scan — which is what keeps
    /// [`bcc_reconstruct`] and [`trilinear_reconstruct`] logarithmic in the
    /// site count instead of linear in it.
    fn find(&self, inv_scaled: &[[f64; 3]; 3], ijk: [i64; 3]) -> Option<usize> {
        let at = self
            .sites
            .partition_point(|s| self.index_of(inv_scaled, *s) < ijk);
        let found = self.sites.get(at)?;
        if self.index_of(inv_scaled, *found) == ijk {
            Some(at)
        } else {
            None
        }
    }
}

/// Build a lattice grid inside `(lo, hi)` holding as close as possible to
/// `target_points` sites, so two lattices can be compared at **matched point
/// density** rather than at matched resolution.
///
/// # How the scale is solved for
///
/// [`Lattice::generator`]'s rows span a fundamental cell of volume exactly 1, so
/// a lattice scaled by `s` holds `1/s³` points per unit volume — identically for
/// all three lattices. With box volume `V` the seed is therefore
/// `s₀ = ∛(V / target)`, one cube root and no lattice-specific arithmetic.
///
/// The seed is only a seed: the realised count is the number of lattice points
/// that survive clipping to the box, which is a **step function** of `s`,
/// non-increasing, and generally misses `target` by a few. So the seed is
/// bracketed (`count(s_lo) ≥ target ≥ count(s_hi)`), bisected 48 times, and the
/// endpoint whose count is nearer `target` is kept; a tie keeps `s_lo`, the
/// denser grid. The result is the closest count this lattice can actually
/// realise in this box, and `sites.len()` — not `target_points` — is what a
/// caller must report. `P-162`'s vacuity control is exactly this: *"both arms
/// must be at genuinely matched point density, reported as counts"*.
///
/// # The protocol a caller must follow, and why it is not symmetric
///
/// The two endpoints of the bracket **are** the two realisable counts either
/// side of `target`, so the count returned is provably the closest this lattice
/// can reach. That does not make the granularity fine. Anchored on the box
/// centre, the cubic lattice can only realise an odd number of sites per axis,
/// so its attainable totals are `15³ = 3375`, `17³ = 4913`, `31³ = 29791`,
/// `33³ = 35937` … — gaps of 30% and more. BCC and FCC interleave two cubic
/// sub-lattices and are far finer grained.
///
/// So the arms must be matched in one direction only: **build the cubic grid
/// first, then target its realised count.**
///
/// ```text
/// let cubic = lattice_grid(Lattice::Cubic, lo, hi, wanted);
/// let bcc   = lattice_grid(Lattice::Bcc, lo, hi, cubic.sites.len());
/// ```
///
/// Measured on `[-2, 2]³`: anchoring on `31³ = 29,791` lands BCC at 29,449
/// (−1.15%) and FCC at 29,659 (−0.44%); anchoring on `63³ = 250,047` lands BCC
/// at 242,649 (−2.96%) and FCC at 246,519 (−1.41%). Asking all three for the
/// same round number instead lands the cubic arm 9% away and matches nothing.
/// The residual mismatch is real and is why `P-162` records `samples` per row
/// rather than one figure for the experiment.
///
/// # Ordering
///
/// Sites are emitted with the first basis coordinate outermost and the third
/// innermost, so `sites` is sorted lexicographically by integer coordinate.
/// Deterministic, and it is what [`LatticeGrid::find`] binary-searches.
pub(crate) fn lattice_grid(
    lattice: Lattice,
    lo: [f64; 3],
    hi: [f64; 3],
    target_points: usize,
) -> LatticeGrid {
    assert!(
        target_points > 0,
        "a lattice grid of zero points is not a grid"
    );
    assert!(
        hi[0] > lo[0] && hi[1] > lo[1] && hi[2] > lo[2],
        "lattice_grid needs a non-degenerate box, got {lo:?}..{hi:?}"
    );

    let volume = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
    let count_at = |s: f64| {
        let mut n = 0usize;
        for_each_site(lattice, lo, hi, s, &mut |_| n += 1);
        n
    };

    let seed = (volume / target_points as f64).cbrt();
    let mut s_lo = seed;
    let mut s_hi = seed;
    let mut guard = 0u32;
    while count_at(s_lo) < target_points {
        s_lo *= 0.75;
        guard += 1;
        assert!(
            guard < 256,
            "no scale dense enough for {target_points} points"
        );
    }
    guard = 0;
    while count_at(s_hi) > target_points {
        s_hi *= 1.25;
        guard += 1;
        assert!(
            guard < 256,
            "no scale sparse enough for {target_points} points"
        );
    }
    for _ in 0..48 {
        let mid = 0.5 * (s_lo + s_hi);
        if count_at(mid) >= target_points {
            s_lo = mid;
        } else {
            s_hi = mid;
        }
    }

    let n_lo = count_at(s_lo);
    let n_hi = count_at(s_hi);
    let scale = if n_lo.abs_diff(target_points) <= n_hi.abs_diff(target_points) {
        s_lo
    } else {
        s_hi
    };

    let mut sites = Vec::with_capacity(n_lo.max(n_hi));
    for_each_site(lattice, lo, hi, scale, &mut |p| sites.push(p));
    LatticeGrid {
        lattice,
        sites,
        scale,
        lo,
        hi,
    }
}

/// Visit every site of `lattice`, scaled by `scale` and anchored on the centre of
/// `(lo, hi)`, that lies inside the closed box.
///
/// The integer range is derived by mapping the box's eight corners through the
/// inverse generator and taking a one-cell margin, so the enumeration covers the
/// box whatever the basis skew and never depends on a hand-picked bound.
fn for_each_site(
    lattice: Lattice,
    lo: [f64; 3],
    hi: [f64; 3],
    scale: f64,
    visit: &mut impl FnMut([f64; 3]),
) {
    let basis = scaled_generator(lattice, scale);
    let inv = invert3(&basis);
    let centre = [
        0.5 * (lo[0] + hi[0]),
        0.5 * (lo[1] + hi[1]),
        0.5 * (lo[2] + hi[2]),
    ];

    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for corner in 0..8u8 {
        let p = [
            if corner & 1 == 0 { lo[0] } else { hi[0] },
            if corner & 2 == 0 { lo[1] } else { hi[1] },
            if corner & 4 == 0 { lo[2] } else { hi[2] },
        ];
        let v = row_times_matrix([p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]], &inv);
        min = [min[0].min(v[0]), min[1].min(v[1]), min[2].min(v[2])];
        max = [max[0].max(v[0]), max[1].max(v[1]), max[2].max(v[2])];
    }
    let a = [
        min[0].floor() as i64 - 1,
        min[1].floor() as i64 - 1,
        min[2].floor() as i64 - 1,
    ];
    let b = [
        max[0].ceil() as i64 + 1,
        max[1].ceil() as i64 + 1,
        max[2].ceil() as i64 + 1,
    ];

    for i in a[0]..=b[0] {
        for j in a[1]..=b[1] {
            for k in a[2]..=b[2] {
                let v = [i as f64, j as f64, k as f64];
                let d = row_times_matrix(v, &basis);
                let p = [centre[0] + d[0], centre[1] + d[1], centre[2] + d[2]];
                if p[0] >= lo[0]
                    && p[0] <= hi[0]
                    && p[1] >= lo[1]
                    && p[1] <= hi[1]
                    && p[2] >= lo[2]
                    && p[2] <= hi[2]
                {
                    visit(p);
                }
            }
        }
    }
}

/// [`Lattice::generator`] with every row multiplied by `scale`.
fn scaled_generator(lattice: Lattice, scale: f64) -> [[f64; 3]; 3] {
    let g = lattice.generator();
    [
        [g[0][0] * scale, g[0][1] * scale, g[0][2] * scale],
        [g[1][0] * scale, g[1][1] * scale, g[1][2] * scale],
        [g[2][0] * scale, g[2][1] * scale, g[2][2] * scale],
    ]
}

/// The row vector `v` times the matrix `m`: `(v·m)_j = Σ_a v_a · m[a][j]`.
///
/// Row-vector convention throughout, because [`Lattice::generator`]'s rows are
/// the basis vectors and a lattice point is `(i,j,k) · B`.
fn row_times_matrix(v: [f64; 3], m: &[[f64; 3]; 3]) -> [f64; 3] {
    [
        v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0],
        v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1],
        v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2],
    ]
}

/// Inverse of a 3×3 matrix by cofactors.
///
/// Every matrix this is asked to invert is a scaled lattice generator, whose
/// determinant is `scale³` by construction — so a zero here means the caller
/// passed a degenerate scale, which is a programming error and not a case to
/// route around.
fn invert3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    let c01 = m[1][2] * m[2][0] - m[1][0] * m[2][2];
    let c02 = m[1][0] * m[2][1] - m[1][1] * m[2][0];
    let det = m[0][0] * c00 + m[0][1] * c01 + m[0][2] * c02;
    assert!(
        det.abs() > 0.0,
        "a lattice generator must be non-singular, got det = {det}"
    );
    let inv_det = 1.0 / det;
    [
        [
            c00 * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            c01 * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            c02 * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}

/// The generated case table for one lattice's natural cell: every distinct sign
/// configuration of the cell's corners, and the triangles it emits.
///
/// Generated by enumeration, never transcribed. `P-162` C3 asks what a
/// non-cubic extractor costs, and these five numbers are the answer in the only
/// currency that does not depend on a machine.
#[derive(Clone, Debug)]
pub(crate) struct CaseTable {
    /// Corners of the natural cell: 8 for the cube, 4 for either tetrahedral
    /// complex.
    pub(crate) corners_per_cell: usize,
    /// `2^corners_per_cell` — the number of sign configurations enumerated.
    pub(crate) cases: usize,
    /// Orbits of those configurations under the cell's **combinatorial**
    /// automorphism group: the 48 signed axis permutations for the cube, the 24
    /// corner permutations `S₄` for a simplex. See [`case_table`] for why
    /// combinatorial and not metric, and why the complement is not folded in.
    pub(crate) distinct_up_to_symmetry: usize,
    /// The largest triangle count any single case emits.
    pub(crate) max_triangles_per_case: usize,
    /// Triangles summed over all `cases` configurations — the table's total
    /// size, and the number a memory budget is set from.
    pub(crate) total_triangles: usize,
}

/// Build, by enumeration, the case table for the given lattice's natural cell.
///
/// # The construction, which is one rule for all three lattices
///
/// A cell edge is **cut** when its two corners are classified differently. The
/// cut edges of one cell link into closed polygons, and a `k`-gon fans into
/// `k − 2` triangles. That is the whole rule; everything lattice-specific is in
/// what the cell's edges are and how the cut edges link.
///
/// - **Cubic (8 corners, 256 cases).** Cut edges are linked face by face: walk
///   each face's four corners counter-clockwise as seen from outside, and pair
///   each entry transition with the exit that closes its own run of inside
///   corners — the *separating* rule, which is what
///   `isomesh::marching_cubes::table::CASES` uses with a zero ambiguity mask.
///   Every cut edge lies on exactly two faces and is an entry on one and an exit
///   on the other, so the segments close into cycles with nothing left over;
///   that is asserted here rather than assumed.
/// - **Tetrahedral (4 corners, 16 cases).** A plane separating a simplex's
///   corners cuts either three edges or four and nothing else, so the polygon is
///   a triangle or a quadrilateral. Asserted, not assumed.
///
/// # The calibration
///
/// For the cubic cell the generated per-case triangle count is compared, entry
/// by entry, against [`shipped_cubic_triangle_counts`]. A mismatch panics. This
/// is the only external check available — the shipped table is itself generated
/// (`marching_cubes/table.rs:177-194`), by a *different* construction that
/// builds directed segment links and picks a chord-safe fan apex — so agreement
/// across all 256 cases is two independent enumerations reaching the same
/// numbers, which is what licenses believing the tetrahedral tables that have
/// nothing to be checked against.
///
/// # The BCC cell decomposition
///
/// BCC's natural cell is its **Delaunay tetrahedron**: two sites from one cubic
/// sub-lattice and two from the other. In the integer coordinates of
/// [`bcc_box_spline`] one representative is
/// `(0,0,0) (2,0,0) (1,1,1) (1,−1,1)` — an edge of the even sub-lattice, an edge
/// of the odd one, skew to it, and the four cross edges. Its volume is `2/3`
/// against a volume-per-site of 4, so there are exactly **6 congruent
/// tetrahedra per lattice site** and twelve per cube of the coarse sub-lattice.
/// This is the decomposition chosen, for three reasons: every circumsphere is
/// empty, so the cell is the *natural* one and not a convention; all four
/// corners are real sample sites, with no interpolated cube centre; and it is
/// the complex the linear box spline of [`bcc_box_spline`] is piecewise linear
/// over, so sampling and reconstruction are described by the same cells. The
/// alternative — six pyramids from a cube centre, each split in two — has the
/// same twelve-per-cube count and the same 4 corners, but its cells are not
/// Delaunay and one corner is not a lattice site.
///
/// # The FCC cell decomposition
///
/// FCC's Delaunay complex is the tetrahedral-octahedral honeycomb, which has
/// **two** cell shapes. To get a single-cell case table each octahedron is split
/// into four tetrahedra along one fixed body diagonal; every octahedron in the
/// honeycomb is a translate of every other, so one global choice of diagonal is
/// consistent everywhere and no seam can disagree. The split is a genuine
/// choice — there are three diagonals — and it changes the piecewise-linear
/// interpolant, but not the case count, because the resulting complex is
/// all-tetrahedral and a tetrahedron has 16 sign configurations however it was
/// obtained.
///
/// # Why the symmetry count is combinatorial
///
/// `distinct_up_to_symmetry` counts orbits under the cell's combinatorial
/// automorphism group — 48 signed axis permutations for the cube, all 24 corner
/// permutations for a simplex — because a case table maps *corner signs* to
/// *edges*, and any relabelling that preserves the edge structure gives an
/// identical entry. The BCC Delaunay tetrahedron's **metric** symmetry group is
/// smaller, order 8, since its six edges come in lengths `2, 2` and four of
/// `√3`; counting orbits under that instead gives 6 rather than 5. The
/// complement `S ↦ Sᶜ` is deliberately **not** folded in: it is not a symmetry
/// of the cell, the shipped table does not treat a case and its complement as
/// one entry, and folding it is what turns the cube's count into the
/// often-quoted 15.
pub(crate) fn case_table(lattice: Lattice) -> CaseTable {
    check_g_literals();
    match lattice {
        Lattice::Cubic => cubic_case_table(),
        Lattice::Bcc | Lattice::Fcc => simplex_case_table(),
    }
}

/// Per-case triangle count of the **shipped** cubic table, read from
/// `isomesh::marching_cubes::table::CASES`.
///
/// The crate's own table, so [`case_table`] has something external to be checked
/// against. Read, never edited; this is the one place in the module where a
/// number comes from outside it.
pub(crate) fn shipped_cubic_triangle_counts() -> [u8; 256] {
    let mut out = [0u8; 256];
    for (slot, case) in out.iter_mut().zip(CASES.iter()) {
        *slot = case.count;
    }
    out
}

/// The 256-case cubic table, enumerated and then calibrated against the shipped
/// one. See [`case_table`].
fn cubic_case_table() -> CaseTable {
    let shipped = shipped_cubic_triangle_counts();
    let mut max = 0usize;
    let mut total = 0usize;
    for (case, expected) in shipped.iter().enumerate() {
        let generated = cubic_triangles_for_case(case as u8);
        assert_eq!(
            generated, *expected as usize,
            "generated cubic case table disagrees with isomesh::marching_cubes::table::CASES \
             at case {case}: {generated} triangles against the shipped {expected}"
        );
        max = max.max(generated);
        total += generated;
    }
    CaseTable {
        corners_per_cell: CUBE_CORNERS,
        cases: CUBE_CASES,
        distinct_up_to_symmetry: orbit_count(CUBE_CORNERS, &cube_symmetries()),
        max_triangles_per_case: max,
        total_triangles: total,
    }
}

/// Triangles one cubic sign configuration emits, by linking cut edges face by
/// face and fanning each closed polygon.
fn cubic_triangles_for_case(case: u8) -> usize {
    let mut cut = [false; CUBE_EDGES];
    for (edge, corners) in EDGE_CORNERS.iter().enumerate() {
        cut[edge] = corner_inside(case, corners[0]) != corner_inside(case, corners[1]);
    }

    // Two incident segments per cut edge, filled in as the faces are walked.
    let mut link = [[NO_LINK; 2]; CUBE_EDGES];
    let mut degree = [0usize; CUBE_EDGES];
    let join =
        |a: u8, b: u8, link: &mut [[u8; 2]; CUBE_EDGES], degree: &mut [usize; CUBE_EDGES]| {
            assert!(
                degree[a as usize] < 2 && degree[b as usize] < 2,
                "a cut edge received a third segment on case {case}"
            );
            link[a as usize][degree[a as usize]] = b;
            degree[a as usize] += 1;
            link[b as usize][degree[b as usize]] = a;
            degree[b as usize] += 1;
        };

    for axis in 0..3usize {
        for side in 0..2u8 {
            let ring = face_corners(axis, side);
            // Start the walk on an outside corner, so the first transition seen
            // is an entry into the solid. An all-inside face is not cut.
            let Some(start) = (0..4).find(|k| !corner_inside(case, ring[*k])) else {
                continue;
            };
            let mut pending: Option<u8> = None;
            for step in 0..4usize {
                let from = (start + step) % 4;
                let to = (from + 1) % 4;
                let inside_from = corner_inside(case, ring[from]);
                let inside_to = corner_inside(case, ring[to]);
                if inside_from == inside_to {
                    continue;
                }
                let edge = edge_index(ring[from], ring[to]);
                if inside_to {
                    assert!(
                        pending.is_none(),
                        "two entries without an exit on case {case}"
                    );
                    pending = Some(edge);
                } else {
                    let entry = pending.take().expect("an exit must follow an entry");
                    join(entry, edge, &mut link, &mut degree);
                }
            }
            assert!(
                pending.is_none(),
                "a face run was never closed on case {case}"
            );
        }
    }

    for edge in 0..CUBE_EDGES {
        assert_eq!(
            degree[edge],
            if cut[edge] { 2 } else { 0 },
            "cut edge {edge} of case {case} does not carry exactly two segments"
        );
    }

    let mut visited = [false; CUBE_EDGES];
    let mut triangles = 0usize;
    for edge in 0..CUBE_EDGES {
        if !cut[edge] || visited[edge] {
            continue;
        }
        let mut length = 0usize;
        let mut previous = NO_LINK;
        let mut current = edge as u8;
        loop {
            visited[current as usize] = true;
            length += 1;
            assert!(length <= CUBE_EDGES, "a cycle did not close on case {case}");
            let next = if link[current as usize][0] == previous {
                link[current as usize][1]
            } else {
                link[current as usize][0]
            };
            if next == edge as u8 {
                break;
            }
            previous = current;
            current = next;
        }
        assert!(length >= 3, "a cycle of length {length} on case {case}");
        triangles += length - 2;
    }
    triangles
}

/// The 16-case table of a tetrahedral cell, shared by BCC and FCC. See
/// [`case_table`].
fn simplex_case_table() -> CaseTable {
    // The six edges of a tetrahedron, lower corner first.
    const EDGES: [[usize; 2]; 6] = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
    let mut max = 0usize;
    let mut total = 0usize;
    for case in 0..SIMPLEX_CASES {
        let inside = |corner: usize| case & (1 << corner) != 0;
        let cuts = EDGES
            .iter()
            .filter(|e| inside(e[0]) != inside(e[1]))
            .count();
        assert!(
            cuts == 0 || cuts == 3 || cuts == 4,
            "a plane through a simplex cuts 0, 3 or 4 edges, not {cuts} (case {case})"
        );
        let triangles = cuts.saturating_sub(2);
        max = max.max(triangles);
        total += triangles;
    }
    CaseTable {
        corners_per_cell: SIMPLEX_CORNERS,
        cases: SIMPLEX_CASES,
        distinct_up_to_symmetry: orbit_count(SIMPLEX_CORNERS, &simplex_symmetries()),
        max_triangles_per_case: max,
        total_triangles: total,
    }
}

/// The 48 corner permutations induced by the cube's combinatorial automorphisms:
/// a permutation of the three axes composed with an independent flip of each.
///
/// Generated, in a fixed order, so the orbit count is deterministic.
fn cube_symmetries() -> Vec<Vec<u8>> {
    const AXIS_PERMUTATIONS: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut out = Vec::with_capacity(48);
    for sigma in AXIS_PERMUTATIONS {
        for flips in 0..8u8 {
            let mut map = Vec::with_capacity(CUBE_CORNERS);
            for corner in 0..CUBE_CORNERS as u8 {
                let mut image = 0u8;
                for (axis, source) in sigma.iter().enumerate() {
                    let value = ((corner >> source) & 1) ^ ((flips >> axis) & 1);
                    image |= value << axis;
                }
                map.push(image);
            }
            out.push(map);
        }
    }
    out
}

/// The 24 corner permutations of a tetrahedron, `S₄`, in lexicographic order.
fn simplex_symmetries() -> Vec<Vec<u8>> {
    let mut out = Vec::with_capacity(24);
    for a in 0..4u8 {
        for b in 0..4u8 {
            for c in 0..4u8 {
                for d in 0..4u8 {
                    if a != b && a != c && a != d && b != c && b != d && c != d {
                        out.push(vec![a, b, c, d]);
                    }
                }
            }
        }
    }
    out
}

/// Orbits of the `2^corners` sign configurations under `group`, counted by
/// canonical representative.
///
/// The representative is the numerically smallest image, and the canonical set
/// is a [`std::collections::BTreeSet`] — no hash iteration order anywhere near
/// a reported number.
fn orbit_count(corners: usize, group: &[Vec<u8>]) -> usize {
    let mut canonical = std::collections::BTreeSet::new();
    for case in 0..(1usize << corners) {
        let mut best = case;
        for map in group {
            let mut image = 0usize;
            for (corner, target) in map.iter().enumerate() {
                if case & (1 << corner) != 0 {
                    image |= 1 << *target as usize;
                }
            }
            best = best.min(image);
        }
        canonical.insert(best);
    }
    canonical.len()
}

/// The linear box spline on the BCC lattice (Entezari, Van De Ville & Möller,
/// `10.1109/tvcg.2007.70429`): the four-direction box spline whose directions
/// are the four body diagonals of the cube. Evaluated at `p` in **lattice
/// coordinates**.
///
/// # Lattice coordinates
///
/// `L = {k ∈ Z³ : k₀ ≡ k₁ ≡ k₂ (mod 2)}` — a cubic lattice of side 2 together
/// with its body centres. That is exactly the lattice
/// `Lattice::Bcc.generator()`'s integer basis generates, so a world point maps
/// in with one division: `p_lattice = (p_world − centre) / (scale / ∛4)`.
///
/// # The construction
///
/// With `Ξ = [ξ₁ ξ₂ ξ₃ ξ₄]`, the four body diagonals
/// `(1,1,1) (−1,1,1) (1,−1,1) (1,1,−1)` as columns, the box spline is
/// `M(x) = ∫_{[0,1]⁴} δ(x − Ξt) dt`. Because there is exactly one direction more
/// than there are dimensions, the fibre `{t ∈ [0,1]⁴ : Ξt = x}` is a **line
/// segment** and `M` is proportional to its length — so the spline needs no
/// recursion and no case table, just four `min`s and four `max`s.
///
/// `Ξ`'s null space is spanned by `u = (−1,1,1,1)`, and `t₀ = (0, (x₁+x₂)/2,
/// (x₀+x₂)/2, (x₀+x₁)/2)` solves `Ξt₀ = x`. The fibre is `t₀ + λu`, and the box
/// constraints `0 ≤ t₀ + λu ≤ 1` reduce to `λ ∈ [L, U]` with
/// `L = max(−1, −b, −c, −d)` and `U = min(0, 1−b, 1−c, 1−d)`. The Jacobian of
/// `(x, λ) ↦ t` is `1/4`, so `M(x) = (U − L)/4` integrates to 1.
///
/// # Normalisation and centring
///
/// `Ξ[0,1]⁴` is centred on `Ξ·(½,½,½,½) = (1,1,1)`, a lattice site, so this
/// function evaluates `M(p + (1,1,1))` to put the peak at the origin. It returns
/// `4·M`, not `M`: the lattice's fundamental cell has volume 4, so `4·M` is the
/// scaling under which `Σ_{k ∈ L} φ(p − k) = 1`. That makes `φ(0) = 1` and
/// `φ(k) = 0` at every other site — the linear box spline **interpolates** — and
/// it is the normalisation [`bcc_reconstruct`] needs.
///
/// # Support and order
///
/// Support is the zonotope of the four diagonals, a **rhombic dodecahedron** of
/// volume 16 in these coordinates — four fundamental cells, against the
/// trilinear's eight — with its 8 vertices at `(±1,±1,±1)` and its 6 at
/// `(±2,0,0)`. Degree is `4 − 3 = 1`, so approximation order is
/// [`BCC_BOX_SPLINE_ORDER`] = 2 and the stencil at a generic point is
/// [`BCC_BOX_SPLINE_STENCIL`] = 4.
pub(crate) fn bcc_box_spline(p: [f64; 3]) -> f64 {
    let x = [p[0] + 1.0, p[1] + 1.0, p[2] + 1.0];
    let b = 0.5 * (x[1] + x[2]);
    let c = 0.5 * (x[0] + x[2]);
    let d = 0.5 * (x[0] + x[1]);
    let low = (-1.0f64).max(-b).max(-c).max(-d);
    let high = 0.0f64.min(1.0 - b).min(1.0 - c).min(1.0 - d);
    (high - low).max(0.0)
}

/// Reconstruct a value at `p` from samples on a BCC grid using
/// [`bcc_box_spline`].
///
/// `values` is parallel to `grid.sites`. Only the sites within the spline's
/// support are visited: the candidate integer coordinates are derived from `p`
/// and looked up by binary search, so the cost is
/// `O(BCC_BOX_SPLINE_STENCIL · log sites)` and not a scan.
///
/// # Panics
///
/// If the weights do not sum to 1 within `1e-9`. That is the exact condition for
/// *"`p` is in the interior where this reconstruction is defined"*: the box
/// spline is a partition of unity over the infinite lattice, so a deficit means
/// a needed site was clipped away by the box. A caller sweeping a field must
/// inset its probes by the support radius — 2 lattice units — rather than accept
/// a quietly degraded value near the boundary.
pub(crate) fn bcc_reconstruct(grid: &LatticeGrid, values: &[f64], p: [f64; 3]) -> f64 {
    assert_eq!(
        grid.lattice,
        Lattice::Bcc,
        "bcc_reconstruct is defined on the BCC lattice, not {}",
        grid.lattice.name()
    );
    assert_eq!(
        values.len(),
        grid.sites.len(),
        "one value per site: {} values against {} sites",
        values.len(),
        grid.sites.len()
    );

    let basis = scaled_generator(grid.lattice, grid.scale);
    let inv = invert3(&basis);
    let centre = grid.centre();
    // One lattice-coordinate unit in world distance.
    let unit = grid.scale / 4f64.cbrt();
    let q = [
        (p[0] - centre[0]) / unit,
        (p[1] - centre[1]) / unit,
        (p[2] - centre[2]) / unit,
    ];

    // A site's lattice coordinate is m = (2i + k, 2j + k, k) for basis
    // coordinate (i, j, k), and the support has radius 2 in the sup norm.
    let mut acc = 0.0;
    let mut weight = 0.0;
    let k_lo = (q[2] - 2.0).ceil() as i64;
    let k_hi = (q[2] + 2.0).floor() as i64;
    for k in k_lo..=k_hi {
        let kf = k as f64;
        let i_lo = ((q[0] - 2.0 - kf) / 2.0).ceil() as i64;
        let i_hi = ((q[0] + 2.0 - kf) / 2.0).floor() as i64;
        let j_lo = ((q[1] - 2.0 - kf) / 2.0).ceil() as i64;
        let j_hi = ((q[1] + 2.0 - kf) / 2.0).floor() as i64;
        for i in i_lo..=i_hi {
            for j in j_lo..=j_hi {
                let m = [(2 * i + k) as f64, (2 * j + k) as f64, kf];
                let w = bcc_box_spline([q[0] - m[0], q[1] - m[1], q[2] - m[2]]);
                if w <= 0.0 {
                    continue;
                }
                if let Some(at) = grid.find(&inv, [i, j, k]) {
                    acc += w * values[at];
                    weight += w;
                }
            }
        }
    }

    assert!(
        (weight - 1.0).abs() < 1e-9,
        "bcc_reconstruct at {p:?} is outside the reconstructible interior: \
         the box spline's weights sum to {weight}, not 1"
    );
    acc
}

/// Trilinear reconstruction on a cubic grid, for the matched baseline.
///
/// `values` is parallel to `grid.sites`. The eight corners of the containing
/// cell are looked up by binary search, exactly as [`bcc_reconstruct`] looks up
/// its four — so the two arms differ in the filter and in nothing else, which is
/// what `P-164` C2 needs in order to attribute a difference to the filter.
///
/// # Panics
///
/// If any of the eight corners was clipped away by the box. The trilinear
/// weights sum to 1 by construction, so unlike [`bcc_reconstruct`] there is no
/// weight deficit to detect and the missing corner is named directly.
pub(crate) fn trilinear_reconstruct(grid: &LatticeGrid, values: &[f64], p: [f64; 3]) -> f64 {
    assert_eq!(
        grid.lattice,
        Lattice::Cubic,
        "trilinear_reconstruct is defined on the cubic lattice, not {}",
        grid.lattice.name()
    );
    assert_eq!(
        values.len(),
        grid.sites.len(),
        "one value per site: {} values against {} sites",
        values.len(),
        grid.sites.len()
    );

    let basis = scaled_generator(grid.lattice, grid.scale);
    let inv = invert3(&basis);
    let centre = grid.centre();
    let q = [
        (p[0] - centre[0]) / grid.scale,
        (p[1] - centre[1]) / grid.scale,
        (p[2] - centre[2]) / grid.scale,
    ];
    let base = [
        q[0].floor() as i64,
        q[1].floor() as i64,
        q[2].floor() as i64,
    ];
    let t = [
        q[0] - base[0] as f64,
        q[1] - base[1] as f64,
        q[2] - base[2] as f64,
    ];

    let mut acc = 0.0;
    for corner in 0..TRILINEAR_STENCIL {
        let d = [corner & 1, (corner >> 1) & 1, (corner >> 2) & 1];
        let w = (if d[0] == 0 { 1.0 - t[0] } else { t[0] })
            * (if d[1] == 0 { 1.0 - t[1] } else { t[1] })
            * (if d[2] == 0 { 1.0 - t[2] } else { t[2] });
        let ijk = [
            base[0] + d[0] as i64,
            base[1] + d[1] as i64,
            base[2] + d[2] as i64,
        ];
        let at = grid.find(&inv, ijk).unwrap_or_else(|| {
            panic!(
                "trilinear_reconstruct at {p:?} needs the clipped-away site {ijk:?}: \
                 inset the probes by one cell"
            )
        });
        acc += w * values[at];
    }
    acc
}

/// Newton iterations allowed when projecting a point onto a field's zero set.
///
/// **1024, not 64, and the difference is a measured verdict rather than a
/// margin.** Gradient-Newton converges linearly, not quadratically, on a CSG
/// crease where the field is a `max` of two smooth pieces: R-162 instrumented
/// the cap crease of `noise_cavity` at `[-1.405105, 0.0, 0.523810]`, `|p| =
/// 1.4996` against the `r = 1.5` sphere, and measured a contraction of about
/// `0.965` per step — 64 steps take `|f|` only from `4.48e-3` to `4.74e-4`, so
/// the projection aborted and no CSV was produced at all. The budget was also
/// silently dropping legitimate probes on every CSG field, which
/// **under-reports** Hausdorff in the truth-to-reconstruction direction:
/// `csg_difference`'s cubic Hausdorff moves `6.135e-2` to `6.570e-2` when the
/// budget rises, and its lattice verdict flips from `-0.177 dB` to `+0.417 dB`.
/// At 1024 the stall count falls from 240 to 86, and 80 of those 86 are the
/// identical value `5.7180e-2` — a genuine positive local minimum of the `max`
/// where no zero lies along any descent path, which is the legitimate drop this
/// module documents. Verdicts are identical at 512 and 1024, so the numbers are
/// insensitive to the budget once it is large enough, and the whole 49^3 sweep
/// costs 37.0 s against 35.2 s — inside this host's governor noise.
const PROJECTION_STEPS: usize = 1024;
/// Seed for the probe stream in [`zero_set_hausdorff`]. Fixed, so the number is
/// reproducible; changing it changes the measurement.
const PROBE_SEED: u64 = 0x1362_A3B5_D1E7_9F11;

/// Symmetric Hausdorff between a point sample of a reconstructed zero set and
/// the true zero set of `sdf`, both restricted to the box.
///
/// The box is the axis-aligned bounding box of `points` — the region the
/// reconstruction actually covered, which is the only box this signature can
/// derive. It is a **seeding region, not a clipping region**, and the difference
/// is load-bearing; see below.
///
/// # The two directions, and why they fail differently
///
/// - **Reconstruction → truth.** Every point of `points` is *claimed* to lie on
///   the surface, so each is Newton-projected onto `sdf`'s zero set and the
///   distance it moved is its error. A projection that does not converge means
///   the claim was false, so it **panics** rather than being skipped: a skipped
///   point would make the reported maximum a lie.
/// - **Truth → reconstruction.** `probes` seeds from a fixed SplitMix64 stream
///   are drawn uniformly in the box and Newton-projected; a seed that fails to
///   converge is simply **not a sample** of the true zero set and is dropped.
///   Each survivor's distance to the nearest point of `points` is measured by
///   exhaustive search, which is `O(probes · points)` and deliberately not
///   accelerated — an approximate nearest neighbour would put a second error
///   term inside an error measurement.
///
/// # Why a landing outside the box is kept
///
/// The obvious refinement — discard a probe whose projection leaves the box — is
/// wrong, and wrong in the direction that flatters the reconstruction. The box
/// comes from `points`, so a reconstruction that **misses part of the surface**
/// shrinks the box to exclude exactly the region it missed, and every probe that
/// would have found the hole is then thrown away as out of bounds. Measured on
/// the unit sphere with the cap `z > 0.6` deleted from a 40,000-point sample:
/// rejecting out-of-box landings reads `0.0159`, indistinguishable from the
/// intact sample's `0.0120`, while keeping them reads `0.9`. So a probe seeded
/// inside the box is followed wherever the field's own gradient takes it. The
/// projection step is clamped, so it cannot take it far.
///
/// The result is the larger of the two directed maxima. `probes` is a parameter
/// precisely so a bench can report it: this direction is a *sample* of an
/// uncountable set, and a Hausdorff distance quoted without its probe count is
/// not reproducible.
///
/// # Panics
///
/// If `points` is empty, if `probes` is zero, if a point of `points` cannot be
/// projected onto the zero set, or if no probe at all survives — the last being
/// the vacuity guard, since a maximum over an empty set would silently read as
/// zero error.
pub(crate) fn zero_set_hausdorff<S: isomesh::Sdf<Scalar = f64>>(
    sdf: &S,
    points: &[[f64; 3]],
    probes: usize,
) -> f64 {
    assert!(
        !points.is_empty(),
        "an empty reconstructed zero set has no Hausdorff distance to anything"
    );
    assert!(probes > 0, "the true zero set needs at least one probe");

    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in points {
        lo = [lo[0].min(p[0]), lo[1].min(p[1]), lo[2].min(p[2])];
        hi = [hi[0].max(p[0]), hi[1].max(p[1]), hi[2].max(p[2])];
    }
    let diagonal = ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2))
        .sqrt()
        .max(f64::MIN_POSITIVE);
    let tolerance = 1e-10 * diagonal;

    let mut worst = 0.0f64;
    for p in points {
        let landed = project(sdf, *p, tolerance, diagonal).unwrap_or_else(|| {
            panic!(
                "a point of the reconstructed zero set, {p:?}, does not project onto the \
                 field's zero set: residual {}",
                sdf.sample(*p)
            )
        });
        worst = worst.max(distance(*p, landed));
    }

    let mut rng = SplitMix64::new(PROBE_SEED);
    let mut survivors = 0usize;
    for _ in 0..probes {
        let seed = [
            lo[0] + (hi[0] - lo[0]) * rng.next_f64(),
            lo[1] + (hi[1] - lo[1]) * rng.next_f64(),
            lo[2] + (hi[2] - lo[2]) * rng.next_f64(),
        ];
        let Some(on_surface) = project(sdf, seed, tolerance, diagonal) else {
            continue;
        };
        survivors += 1;
        let mut nearest = f64::INFINITY;
        for q in points {
            nearest = nearest.min(distance(on_surface, *q));
        }
        worst = worst.max(nearest);
    }
    assert!(
        survivors > 0,
        "none of {probes} probes landed on the field's zero set, so the \
         truth-to-reconstruction direction was measured over an empty set"
    );
    worst
}

/// Newton-project `p` onto `sdf`'s zero set along the gradient.
///
/// The step is clamped to a quarter of the box diagonal so a near-flat gradient
/// cannot fling the iterate across the domain, which is the one way this
/// iteration diverges in practice. `None` means it did not reach `tolerance`.
fn project<S: isomesh::Sdf<Scalar = f64>>(
    sdf: &S,
    p: [f64; 3],
    tolerance: f64,
    diagonal: f64,
) -> Option<[f64; 3]> {
    let limit = 0.25 * diagonal;
    let mut at = p;
    for _ in 0..PROJECTION_STEPS {
        let value = sdf.sample(at);
        if value.abs() <= tolerance {
            return Some(at);
        }
        let g = sdf.gradient(at);
        let square = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
        // A NaN square fails `is_finite`, so both a flat and a broken gradient
        // land here.
        if square <= 0.0 || !square.is_finite() {
            return None;
        }
        let mut step = [
            value * g[0] / square,
            value * g[1] / square,
            value * g[2] / square,
        ];
        let length = (step[0] * step[0] + step[1] * step[1] + step[2] * step[2]).sqrt();
        if length > limit {
            let shrink = limit / length;
            step = [step[0] * shrink, step[1] * shrink, step[2] * shrink];
        }
        at = [at[0] - step[0], at[1] - step[1], at[2] - step[2]];
        if !at[0].is_finite() || !at[1].is_finite() || !at[2].is_finite() {
            return None;
        }
    }
    if sdf.sample(at).abs() <= tolerance {
        Some(at)
    } else {
        None
    }
}

/// Euclidean distance between two points.
fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// SplitMix64, the ten-line seeded generator this module needs and the reason it
/// needs no dependency.
///
/// Vigna's finaliser, from the reference C in the public-domain `splitmix64.c`.
/// Used only to place probe seeds in a box, where the requirement is
/// reproducibility rather than statistical quality.
#[derive(Clone, Debug)]
struct SplitMix64 {
    /// The 64-bit state, advanced by the golden-ratio increment.
    state: u64,
}

impl SplitMix64 {
    /// A stream from `seed`.
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// The next 64 bits.
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// The next value in `[0, 1)`, from the top 53 bits — the only bits an `f64`
    /// mantissa can hold without rounding.
    fn next_f64(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) * (1.0 / 9_007_199_254_740_992.0)
    }
}
