//! Triply periodic minimal surfaces as nodal fields, the periodic-wrap
//! identification, and the Euler reader the wrap is judged by.
//!
//! Ticket: R-142, which owns this module. Consumed **unchanged** by R-143
//! (Schwarz P and D as the second and third fields with exactly known
//! topology), R-144 (whether a periodic value-noise terrain admits a `chi`
//! oracle at all) and R-145 (stratified-Morse ground truth, cross-checked
//! against R-142's analytic prediction on the gyroid).
//!
//! # The three surfaces and the three predictions
//!
//! With the domain scaled so that one period is `2*pi` on every axis, the
//! standard nodal (Gandy–Klinowski) approximations are
//!
//! ```text
//! gyroid     F_G = sin x cos y + sin y cos z + sin z cos x
//! Schwarz P  F_P = cos x + cos y + cos z
//! Schwarz D  F_D = sin x sin y sin z + sin x cos y cos z
//!                + cos x sin y cos z + cos x cos y sin z
//! ```
//!
//! All three have **genus 3 in their own primitive translational cell**, so
//! `chi = 2 - 2*3 = -4` there. What differs between them is how many primitive
//! cells fit inside the *conventional cubic cell* of side one period, and that
//! is the whole content of the `-8 / -4 / -16` spread:
//!
//! | surface | space group | translation lattice | primitive cells per cubic cell | `chi` per cubic cell | genus per cubic cell |
//! |---|---|---|---|---|---|
//! | gyroid | `Ia-3d` | body-centred cubic | 2 | **-8** | 5 |
//! | Schwarz P | `Im-3m` | simple cubic | 1 | **-4** | 3 |
//! | Schwarz D | `Pn-3m` | face-centred cubic | 4 | **-16** | 9 |
//!
//! `chi` is additive over a disjoint decomposition, so `N` periods per axis
//! give `chi = N^3 * chi_cell`: the `-8*N^3` R-142 registered, and the `-4*N^3`
//! and `-16*N^3` R-143 registered.
//!
//! # The lattice index is an identity, not a lookup
//!
//! The extra translations in the middle column are exactly the shifts `t` for
//! which the nodal function is *invariant* — `F(p + t) = F(p)` — as opposed to
//! merely *negated*, `F(p + t) = -F(p)`. A negating shift maps the zero set to
//! itself but exchanges the two labyrinths, so it is a symmetry of the surface
//! and **not** a translation of the labelled structure, and it therefore does
//! not shrink the translational cell. Substituting the body-centring shift
//! `(x, y, z) -> (x + pi, y + pi, z + pi)`, where every `sin` and every `cos`
//! flips sign:
//!
//! ```text
//! F_G(p + (pi,pi,pi)) = (-sin x)(-cos y) + (-sin y)(-cos z) + (-sin z)(-cos x)
//!                     = + F_G(p)        two flips per term  -> INVARIANT
//!
//! F_P(p + (pi,pi,pi)) = (-cos x) + (-cos y) + (-cos z)
//!                     = - F_P(p)        one flip per term   -> NEGATED
//!
//! F_D(p + (pi,pi,pi)) = each term is a product of THREE trig factors,
//!                       so each term picks up (-1)^3
//!                     = - F_D(p)                            -> NEGATED
//! ```
//!
//! That is P-143 C2's mechanism, and it is why the gyroid's cubic cell holds
//! two primitive cells while P's holds one. The face-centring shift
//! `(pi, pi, 0)` flips two of the three factors and so completes the picture:
//! `F_D` is invariant under it (giving Schwarz D the face-centred lattice and
//! its four primitive cells), while `F_G` and `F_P` are neither invariant nor
//! negated by it. [`shift_residuals`] evaluates any of these claims
//! numerically; [`Tpms::lattice_shifts`] is the list the counts are derived
//! from, so [`Tpms::chi_per_cubic_cell`] cannot drift away from its own
//! justification.
//!
//! # Sign convention, and what it costs
//!
//! [`Sdf::sample`] returns the nodal value **directly**: negative inside one
//! labyrinth, positive in the other, zero on the surface. This is a level-set
//! function, **not** a signed distance field — `|grad F|` is nowhere near 1
//! (it vanishes on the whole singular skeleton). So
//! `isomesh::validate::accuracy`, `field_bound_report` and anything else that
//! reads the field value as a distance are **meaningless on this field** and
//! must not be run on it. `chi` needs only the sign of the field, which is why
//! it is the invariant this phase can gate on.
//!
//! For the same reason [`NodalTpms`] deliberately does **not** implement
//! `isomesh::fields::ReferenceField`: `common::grid` derives its cell size
//! from `(hi - lo) / (samples - 1)` for a caller-chosen `samples`, which is
//! precisely the non-periodic-conforming grid P-142 names as the defect ("the
//! extraction box is not periodic-conforming"). Use
//! [`NodalTpms::periodic_grid`], which takes voxels *per period* and can only
//! produce a conforming grid.
//!
//! # What `wrap_seams` does, and what it destroys
//!
//! A closed surface on the 3-torus cannot be embedded in a box: identifying
//! opposite faces is a **combinatorial** operation, and afterwards some
//! triangles have one corner at one side of the domain and another corner at
//! the far side. That is not a defect, it is what "closes on the torus" means.
//! So after [`wrap_seams`] the buffer is a valid *simplicial complex* and an
//! invalid *geometric mesh*: connectivity readings (`chi`, components, genus,
//! boundary and non-manifold edge counts) are exact, and every metric reading
//! (area, mean ratio, Hausdorff distance, self-intersections, degenerate
//! triangles) is nonsense. Take the non-wrapped arm for those.
//!
//! # Measured: which sample grids the oracle survives
//!
//! Everything below was measured through this module with `MarchingCubes`
//! (default `FaceAmbiguity::Separate`; the asymptotic decider and the trilinear
//! interior test change nothing here) at `tol = isomesh::weld::epsilon_for(h)`,
//! over 168 configurations: `voxels_per_period` from 8 to 48 at `N = 1`, plus
//! `{32, 33, 56, 64, 65, 96, 97, 128, 129}` at `N = 1` and `{32, 33}` at
//! `N` in `{1, 2, 3}`.
//!
//! - **`boundary_edges` was 0 in every single wrapped run.** The wrap closes the
//!   mesh at every resolution and every period count tested, on all three
//!   surfaces. That is the claim [`wrap_seams`] exists to support.
//! - **The gyroid gave `chi = -8*N^3` in every run**, and Schwarz P gave
//!   `-4*N^3` in all but one.
//! - **Schwarz D failed at exactly the resolutions where `voxels_per_period` is
//!   a multiple of 8**, and only there: `-12` instead of `-16` at 32 and 56,
//!   `-9` at 64, `-7` at 96, `+1` at 128 — while 15, 17, 31, 33, 65, 97, 129
//!   and every even non-multiple of 8 were exact. The mechanism is M-48's, not
//!   the wrap's: a multiple of 8 puts samples on the `pi/4` lattice, where
//!   `F_D`'s four terms are equal in magnitude and cancel to **exactly** zero
//!   (e.g. at `(pi/4, pi/4, 3*pi/4)`), so the crossing parameter is 0 or 1, the
//!   cut edges of one cell place coincident vertices, and the weld turns them
//!   into a pinch. Schwarz P shows the same failure once, at `N = 3, v = 33`,
//!   from ordinary floating-point cancellation rather than an exact lattice.
//! - **In all 12 pinching runs, `chi_measured - chi_predicted` equalled
//!   `non_manifold_edges` exactly** — 4, 4, 7, 9, 17, 32, 135, 5, 6 and so on,
//!   never off by one. Each pinch merges two sheets and costs exactly one from
//!   `chi`.
//!
//! So: **choose an odd `voxels_per_period` (33, 65, 97, 129) and record
//! `non_manifold_edges` beside `chi`.** A non-zero count says "a sample landed
//! on the isosurface", names how much of the `chi` gap it accounts for, and is
//! a statement about the extractor's degenerate-crossing handling rather than
//! about the `-8/-4/-16` prediction. The registration prose for P-142 names 32,
//! 64, 96 and 128 voxels per period — every one of them a multiple of 8, and so
//! every one of them in Schwarz D's failing family. That is a fact about the
//! grid, reportable as such, and the analytic predictions are unaffected: they
//! reproduce exactly on 118 clean runs.
//!
//! One more control-arm caveat, measured: the non-wrapped arm is recognised by
//! `boundary_edges > 0` (hundreds to thousands), **not** by its `chi`. The
//! non-wrapped gyroid reads `-3` against `-8` and non-wrapped Schwarz D reads
//! `-11` against `-16`, but non-wrapped Schwarz P reads `-4` at `N = 1` and
//! `-32` at `N = 2` — its own prediction, by coincidence of the caps the box
//! cuts. A vacuity control that only compares `chi` would pass the wrong arm on
//! one field in three.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use isomesh::weld::Welder;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Residual below which a shift identity is called exact.
///
/// The nodal functions are sums of at most four products of unit-magnitude
/// terms, so an exact identity evaluates to a few ulp of `4.0`, i.e. under
/// `1e-15`; a broken one is `O(1)`. Any threshold in between decides the same
/// way, and `1e-12` is three orders clear of both.
pub(crate) const SHIFT_RESIDUAL_TOLERANCE: f64 = 1e-12;

/// Samples per axis in the deterministic grid the shift identities are checked
/// on.
const SHIFT_GRID: u32 = 17;

/// The three nodal triply periodic minimal surfaces this phase measures.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Tpms {
    /// Schoen's gyroid, `Ia-3d`, the field the crate already ships capped and
    /// the one `CLAUDE.md` records as having no `chi` gate.
    Gyroid,
    /// Schwarz' primitive surface, `Im-3m`.
    SchwarzP,
    /// Schwarz' diamond surface, `Pn-3m`.
    SchwarzD,
}

impl Tpms {
    /// All three, in the order R-142 and R-143 report them.
    pub(crate) const ALL: [Tpms; 3] = [Tpms::Gyroid, Tpms::SchwarzP, Tpms::SchwarzD];

    /// The CSV `field` column: `gyroid`, `schwarz_p`, `schwarz_d`.
    pub(crate) fn name(self) -> &'static str {
        match self {
            Tpms::Gyroid => "gyroid",
            Tpms::SchwarzP => "schwarz_p",
            Tpms::SchwarzD => "schwarz_d",
        }
    }

    /// The crystallographic space group of the surface.
    pub(crate) fn space_group(self) -> &'static str {
        match self {
            Tpms::Gyroid => "Ia-3d",
            Tpms::SchwarzP => "Im-3m",
            Tpms::SchwarzD => "Pn-3m",
        }
    }

    /// The Bravais lattice of the surface's **translation** group, as a short
    /// tag for the CSV: `bcc`, `simple_cubic`, `fcc`.
    ///
    /// This is the lattice generated by the three period vectors together with
    /// [`lattice_shifts`](Self::lattice_shifts) — the shifts that leave the
    /// nodal function *invariant*. It is not read off the space-group symbol:
    /// `Im-3m` is body-centred, but the body-centring operation of Schwarz P
    /// **negates** `F_P` and so exchanges its two labyrinths instead of
    /// translating the structure, which is exactly why P's cubic cell holds one
    /// primitive cell and carries `chi = -4` rather than `-8`. Reporting `bcc`
    /// here beside `chi_per_cubic_cell = -4` would put a contradiction in one
    /// CSV row and would falsify P-143 C2's own explanation of the factor of
    /// two.
    pub(crate) fn primitive_lattice(self) -> &'static str {
        match self {
            Tpms::Gyroid => "bcc",
            Tpms::SchwarzP => "simple_cubic",
            Tpms::SchwarzD => "fcc",
        }
    }

    /// The non-trivial centring translations, in units where one period is
    /// `2*pi`, under which the nodal function is invariant.
    ///
    /// Together with the identity these are the cosets of the simple cubic
    /// lattice inside the surface's own translation lattice, so their count
    /// plus one is [`primitive_cells_per_cubic_cell`](Self::primitive_cells_per_cubic_cell).
    /// Every entry is checkable with [`shift_residuals`], and every vector
    /// *absent* from this list must fail that check — which is what makes the
    /// count a measurement rather than a transcription.
    pub(crate) fn lattice_shifts(self) -> &'static [[f64; 3]] {
        match self {
            // Body centring: the half diagonal.
            Tpms::Gyroid => &[[PI, PI, PI]],
            // Nothing beyond the periods themselves.
            Tpms::SchwarzP => &[],
            // Face centring: the three half face diagonals.
            Tpms::SchwarzD => &[[PI, PI, 0.0], [PI, 0.0, PI], [0.0, PI, PI]],
        }
    }

    /// How many primitive translational cells tile the conventional cubic cell:
    /// 2 (bcc), 1 (simple cubic), 4 (fcc).
    pub(crate) fn primitive_cells_per_cubic_cell(self) -> i64 {
        1 + self.lattice_shifts().len() as i64
    }

    /// Euler characteristic of the surface inside one conventional cubic cell of
    /// side one period, on the 3-torus: `-8` (gyroid), `-4` (Schwarz P), `-16`
    /// (Schwarz D).
    ///
    /// Genus 3 per primitive cell gives `chi = 2 - 2*3 = -4` there, and `chi`
    /// is additive over the
    /// [`primitive_cells_per_cubic_cell`](Self::primitive_cells_per_cubic_cell)
    /// copies that tile the cubic cell. Written as that product rather than as
    /// three literals so the number and its explanation cannot part company.
    pub(crate) fn chi_per_cubic_cell(self) -> i64 {
        -4 * self.primitive_cells_per_cubic_cell()
    }

    /// Is the nodal function invariant (rather than negated) under the
    /// body-centring shift `(pi, pi, pi)`?
    ///
    /// True for the gyroid, false for P and D — this is P-143 C2's mechanism and
    /// the reason the gyroid's conventional cubic cell holds two primitive
    /// cells. The claim is *asserted* here and *checked* by
    /// [`body_centring_check`].
    pub(crate) fn body_centring_invariant(self) -> bool {
        match self {
            Tpms::Gyroid => true,
            Tpms::SchwarzP | Tpms::SchwarzD => false,
        }
    }
}

/// The nodal function of `kind` at `p`, with one period equal to `2*pi` per
/// axis.
///
/// Negative inside one labyrinth, positive in the other. See the module header
/// for why this is not a distance.
pub(crate) fn nodal(kind: Tpms, p: [f64; 3]) -> f64 {
    let (sx, cx) = (p[0].sin(), p[0].cos());
    let (sy, cy) = (p[1].sin(), p[1].cos());
    let (sz, cz) = (p[2].sin(), p[2].cos());
    match kind {
        Tpms::Gyroid => sx * cy + sy * cz + sz * cx,
        Tpms::SchwarzP => cx + cy + cz,
        Tpms::SchwarzD => sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz,
    }
}

/// `(max |F(p + shift) - F(p)|, max |F(p + shift) + F(p)|)` over a fixed grid.
///
/// The first component is the residual of the *invariance* claim, the second of
/// the *negation* claim, and reporting both is what makes either non-vacuous: a
/// grid on which `F` happened to vanish would return two zeros, and a reader
/// can see that instead of being told "invariant" twice.
///
/// # The grid, and why it is offset
///
/// `SHIFT_GRID^3 = 4913` points at `2*pi*(i + o_k)/17` with per-axis offsets
/// `1/3`, `1/5`, `1/7`. Those offsets are what keeps every sample **off** the
/// symmetry planes: a coordinate is a multiple of `pi/2` only if
/// `4*(i + o_k)/17` is an integer, and `4/3`, `4/5` and `4/7` are not integers,
/// so no `sin` and no `cos` at any sample is zero. On a naive grid (`o_k = 0`)
/// a large fraction of samples sit on mirror planes where both residuals
/// collapse and the check stops discriminating.
pub(crate) fn shift_residuals(kind: Tpms, shift: [f64; 3]) -> (f64, f64) {
    let offsets = [1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0];
    let coord =
        |i: u32, axis: usize| 2.0 * PI * (f64::from(i) + offsets[axis]) / f64::from(SHIFT_GRID);
    let mut symmetric = 0.0f64;
    let mut antisymmetric = 0.0f64;
    for i in 0..SHIFT_GRID {
        let x = coord(i, 0);
        for j in 0..SHIFT_GRID {
            let y = coord(j, 1);
            for k in 0..SHIFT_GRID {
                let p = [x, y, coord(k, 2)];
                let here = nodal(kind, p);
                let there = nodal(kind, [p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]]);
                symmetric = symmetric.max((there - here).abs());
                antisymmetric = antisymmetric.max((there + here).abs());
            }
        }
    }
    (symmetric, antisymmetric)
}

/// [`shift_residuals`] at the body-centring shift `(pi, pi, pi)`.
pub(crate) fn body_centring_residuals(kind: Tpms) -> (f64, f64) {
    shift_residuals(kind, [PI, PI, PI])
}

/// The residual of the body-centring relation `kind` claims, and whether the
/// claim held.
///
/// The residual is `max |F(p + (pi,pi,pi)) - F(p)|` when
/// [`Tpms::body_centring_invariant`] is true and
/// `max |F(p + (pi,pi,pi)) + F(p)|` when it is false, over
/// [`shift_residuals`]'s grid.
///
/// The verdict requires **both** that the claimed relation holds within
/// [`SHIFT_RESIDUAL_TOLERANCE`] and that the opposite relation does not. The
/// second half is the vacuity guard: `F = 0` satisfies invariance and negation
/// at once, so a residual of zero on its own establishes nothing.
pub(crate) fn body_centring_check(kind: Tpms) -> (f64, bool) {
    let (symmetric, antisymmetric) = body_centring_residuals(kind);
    let (claimed, opposite) = if kind.body_centring_invariant() {
        (symmetric, antisymmetric)
    } else {
        (antisymmetric, symmetric)
    };
    (
        claimed,
        claimed <= SHIFT_RESIDUAL_TOLERANCE && opposite > SHIFT_RESIDUAL_TOLERANCE,
    )
}

/// A nodal TPMS over `periods` periods per axis on the domain
/// `[0, 2*pi*periods]^3`.
///
/// The domain is exactly `periods^3` conventional cubic cells, so the predicted
/// Euler characteristic under periodic wrap is
/// [`chi_predicted`](Self::chi_predicted)`= periods^3 * kind.chi_per_cubic_cell()`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct NodalTpms {
    /// Which surface.
    pub(crate) kind: Tpms,
    /// Periods per axis; the domain grows with it, the function does not scale.
    pub(crate) periods: u32,
}

impl NodalTpms {
    /// # Panics
    ///
    /// If `periods` is zero, which would give an empty domain and a predicted
    /// `chi` of zero — a number that would look like a measurement.
    pub(crate) fn new(kind: Tpms, periods: u32) -> Self {
        assert!(periods >= 1, "a TPMS domain needs at least one full period");
        Self { kind, periods }
    }

    /// `(lo, hi)` of the extraction box: exactly `periods` full periods per
    /// axis, from the origin.
    pub(crate) fn domain(&self) -> ([f64; 3], [f64; 3]) {
        let span = 2.0 * PI * f64::from(self.periods);
        ([0.0; 3], [span; 3])
    }

    /// `periods^3 * kind.chi_per_cubic_cell()`, the number R-142 and R-143
    /// record as `chi_predicted`.
    pub(crate) fn chi_predicted(&self) -> i64 {
        let n = i64::from(self.periods);
        self.kind.chi_per_cubic_cell() * n * n * n
    }

    /// The **periodic-conforming** sample grid: `(shape, origin, cell_size)`.
    ///
    /// `voxels_per_period` cells span one period, so `cell_size = 2*pi /
    /// voxels_per_period` divides the period exactly and the sample at `hi` is
    /// the same period point as the sample at `lo`. That identity is what makes
    /// the two opposite boundary faces carry the same sign configuration and
    /// hence the same cut, which is the only reason [`wrap_seams`] can close the
    /// mesh. A grid whose spacing does not divide the period cannot be wrapped
    /// no matter how the seam is welded, and P-142 names that as the defect the
    /// project mistook for "gyroid has no `chi`".
    ///
    /// Sample count is `voxels_per_period * periods + 1` per axis, because
    /// `Shape3::size` counts samples and `n` samples span `n - 1` cells.
    ///
    /// # Panics
    ///
    /// If `voxels_per_period` is zero, or if the sample grid does not fit `u32`.
    pub(crate) fn periodic_grid(&self, voxels_per_period: u32) -> (RuntimeShape3, [f64; 3], f64) {
        assert!(
            voxels_per_period >= 1,
            "a period needs at least one voxel to be sampled at all"
        );
        let samples = voxels_per_period * self.periods + 1;
        let shape = RuntimeShape3::new([samples; 3]).expect("TPMS periodic grid fits u32");
        let cell_size = 2.0 * PI / f64::from(voxels_per_period);
        (shape, self.domain().0, cell_size)
    }
}

impl Sdf for NodalTpms {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        nodal(self.kind, p)
    }
}

/// Identify opposite boundary faces of a periodic extraction by the period
/// translation, welding coincident boundary vertices so the mesh closes on the
/// 3-torus. Returns the number of vertex pairs identified.
///
/// Each identification merges one vertex into another, so the return value is
/// both the number of pairwise identifications performed and the number of
/// vertices the buffer lost. A corner where three identifications meet holds
/// eight coincident vertices and contributes seven.
///
/// # How
///
/// Every coordinate within `tol` of `hi[k]` is folded down to `lo[k]` **in a
/// matching key**, not in the mesh, and the keys are then welded with
/// [`isomesh::weld::Welder`]: vertices are walked in buffer order and each
/// joins the lowest-indexed representative within `tol`, on a sorted integer
/// lattice of side `tol`. The surviving buffer keeps each representative's
/// **original** position and normal, and indices are rewritten with triangles
/// left with two equal corners dropped.
///
/// The fold lives in the key because folding the mesh would move a boundary
/// vertex a whole period away from the triangle that owns it, which is a
/// different mesh, not a wrapped one. Either way the combinatorics — and so
/// every number [`euler`] reports — are identical; see the module header for
/// what the wrap costs geometrically.
///
/// The weld is the crate's shipped welder rather than a bench-local copy, so
/// what R-142 measures is the mechanism the project would land, and the
/// determinism argument (a sorted `Vec` broadphase, lowest-indexed
/// representative, no hash iteration anywhere) is the one already stated in
/// `isomesh::weld`'s module docs. `tol` should be
/// `isomesh::weld::epsilon_for(cell_size)`: the broadphase narrows its lattice
/// key through `f32`, so the usable domain extent is about `2^24 * tol`, which
/// a period-scaled tolerance never approaches and an absolute `1e-9` would.
///
/// # Panics
///
/// If `tol` is not finite and positive, if any axis of the box is not longer
/// than `2 * tol`, if `mesh.normals` is not one per vertex, or if an index is
/// out of range.
pub(crate) fn wrap_seams(mesh: &mut MeshBuffer<f64>, lo: [f64; 3], hi: [f64; 3], tol: f64) -> u64 {
    assert!(
        tol.is_finite() && tol > 0.0,
        "the seam tolerance must be finite and positive"
    );
    for (l, h) in lo.iter().zip(hi.iter()) {
        assert!(
            h - l > 2.0 * tol,
            "the wrap box must be longer than the seam tolerance on every axis"
        );
    }
    assert_eq!(
        mesh.normals.len(),
        mesh.positions.len(),
        "a MeshBuffer carries one normal per vertex"
    );

    let count = mesh.positions.len();
    let mut keys = MeshBuffer::<f64>::new();
    keys.positions = mesh
        .positions
        .iter()
        .map(|p| fold_to_lo(*p, lo, hi, tol))
        .collect();
    // The welder reads normals only to compact them, and this buffer's normals
    // are discarded; the real ones travel with the representatives below.
    keys.normals = vec![[0.0; 3]; count];
    keys.indices = mesh.indices.clone();

    let mut welder = Welder::<f64>::new();
    let report = welder
        .weld(&mut keys, tol)
        .expect("the seam tolerance is positive and every index is in range");
    let remap = welder.remap();

    let survivors = keys.positions.len();
    let mut positions = vec![[0.0f64; 3]; survivors];
    let mut normals = vec![[0.0f64; 3]; survivors];
    let mut written = vec![false; survivors];
    // Ascending input order, and the welder gives a representative an output
    // index no greater than its input index, so the first vertex to reach a
    // given output slot is that slot's representative.
    for (input, &output) in remap.iter().enumerate() {
        let output = output as usize;
        if !written[output] {
            written[output] = true;
            positions[output] = mesh.positions[input];
            normals[output] = mesh.normals[input];
        }
    }

    mesh.positions = positions;
    mesh.normals = normals;
    mesh.indices = keys.indices;
    report.vertices_removed() as u64
}

/// One position with every coordinate within `tol` of the far face folded onto
/// the near one.
///
/// This is the period translation, written per axis so that an edge of the box
/// folds twice and a corner three times — the multi-axis identifications P-142
/// needs and the ones a single "is this the `x = hi` face" test would miss.
fn fold_to_lo(p: [f64; 3], lo: [f64; 3], hi: [f64; 3], tol: f64) -> [f64; 3] {
    let mut folded = p;
    for ((slot, &l), &h) in folded.iter_mut().zip(lo.iter()).zip(hi.iter()) {
        if (h - *slot).abs() <= tol {
            *slot = l;
        }
    }
    folded
}

/// What [`euler`] counted.
///
/// `boundary_edges` and `non_manifold_edges` are here so a caller can tell a
/// closed surface from an open one *before* believing `chi`: on a
/// periodic-conforming extraction they must both be zero after
/// [`wrap_seams`], and the non-wrapped control arm is recognised by a
/// `boundary_edges` in the thousands rather than by its `chi` alone.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct EulerCount {
    /// `vertices - edges + faces`.
    pub(crate) chi: i64,
    /// Welded vertices a kept triangle actually names.
    pub(crate) vertices: u64,
    /// Distinct undirected edges.
    pub(crate) edges: u64,
    /// Non-degenerate triangles.
    pub(crate) faces: u64,
    /// Edges used by exactly one triangle.
    pub(crate) boundary_edges: u64,
    /// Edges used by three or more triangles.
    pub(crate) non_manifold_edges: u64,
}

/// `V - E + F` of an indexed triangle mesh, counting each undirected edge once
/// and welding positions within `tol` first.
///
/// # The three counting decisions, stated because each is a choice
///
/// **Weld first.** Marching Cubes shares a vertex only between cells meeting on
/// a grid *edge*, so a sample landing on the isosurface leaves genuinely
/// coincident duplicates (M-48), and `chi` on an unwelded buffer counts one
/// surface as several. The weld is the same [`isomesh::weld::Welder`] and the
/// same rule [`wrap_seams`] uses — one mechanism, so the two readings cannot
/// disagree about what "coincident" means.
///
/// **`vertices` counts referenced vertices**, not buffer length: welding is not
/// garbage collection and leaves unreferenced vertices in place, and a vertex
/// no triangle names is not part of the surface whose `chi` this is. This is
/// also exactly how the crate defines the field it will be cross-checked
/// against — `isomesh::validate`'s `euler_characteristic` is
/// `referenced_vertices - edges + faces`.
///
/// **Degenerate triangles do not count as faces.** A triangle with two equal
/// corners contributes no face and no edge; it is a triangle the weld removed
/// the area from, and counting it would break `V - E + F` rather than record
/// anything.
///
/// Components and genus are deliberately not returned:
/// `isomesh::validate::validate_indexed` already computes them, and a second
/// union-find here would be a second answer to one question.
///
/// # Panics
///
/// If `tol` is not finite and positive, or if an index is out of range.
pub(crate) fn euler(positions: &[[f64; 3]], indices: &[u32], tol: f64) -> EulerCount {
    assert!(
        tol.is_finite() && tol > 0.0,
        "the weld tolerance must be finite and positive"
    );

    let mut welded = MeshBuffer::<f64>::new();
    welded.positions = positions.to_vec();
    welded.normals = vec![[0.0; 3]; positions.len()];
    welded.indices = indices.to_vec();
    let mut welder = Welder::<f64>::new();
    welder
        .weld(&mut welded, tol)
        .expect("the weld tolerance is positive and every index is in range");

    let mut referenced = vec![false; welded.positions.len()];
    let mut edges: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    let mut faces = 0u64;
    // `as_chunks` over `chunks_exact` is the crate's own convention for walking
    // an index buffer (collider.rs:67, normals.rs:167, validate.rs:701); it
    // drops a ragged tail the same way.
    for tri in welded.indices.as_chunks::<3>().0 {
        let (a, b, c) = (tri[0], tri[1], tri[2]);
        // The welder drops collapsed triangles only when it merged something,
        // so a caller's own repeated index can still arrive here.
        if a == b || b == c || a == c {
            continue;
        }
        faces += 1;
        for &v in &[a, b, c] {
            referenced[v as usize] = true;
        }
        for &(u, v) in &[(a, b), (b, c), (c, a)] {
            let key = if u < v { (u, v) } else { (v, u) };
            *edges.entry(key).or_insert(0) += 1;
        }
    }

    let vertices = referenced.iter().filter(|seen| **seen).count() as u64;
    let mut boundary_edges = 0u64;
    let mut non_manifold_edges = 0u64;
    for &uses in edges.values() {
        if uses == 1 {
            boundary_edges += 1;
        } else if uses >= 3 {
            non_manifold_edges += 1;
        }
    }
    let edge_count = edges.len() as u64;

    EulerCount {
        chi: vertices as i64 - edge_count as i64 + faces as i64,
        vertices,
        edges: edge_count,
        faces,
        boundary_edges,
        non_manifold_edges,
    }
}
