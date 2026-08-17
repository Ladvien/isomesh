//! The extractor — [`roots`](super::roots) into [`surface`](super::surface),
//! over a grid.
//!
//! This is where the two halves of subgrid Marching Tetrahedra meet: 1D root
//! finding replaces the sign test, and §3.2's reconstruction triangulates
//! whatever that finds, however many crossings an edge carries.
//!
//! # What it does that the others cannot
//!
//! [`MarchingTetrahedra`](crate::marching_tetrahedra) asks one question per edge
//! and gets one bit. This asks for every zero along that edge. The difference is
//! not incremental: M-67 measured that a sign test cannot distinguish **95.6%**
//! of the configurations a tetrahedron can be in, and A-005 measured
//! `thin_plate` — a plate 0.4 cells thick — returning **zero triangles** from
//! greedy quads, because no cell centre is inside it.
//!
//! # Conventions
//!
//! Identical to [`marching_cubes`](crate::marching_cubes) and
//! [`marching_tetrahedra`](crate::marching_tetrahedra), deliberately: sign
//! negative-inside, zero counts as outside, normals the field's own gradient,
//! winding counter-clockwise seen from outside the solid. The last of those is
//! imposed here rather than inherited — see
//! [`extract`](SubgridMarchingTetrahedra::extract).
//!
//! # Vertices are shared, and the key is a property of the grid
//!
//! A crossing is the `index`-th root along a tetrahedron edge, and
//! [`FacePoint`] already counts `index` from that
//! edge's lower-numbered corner for exactly this reason. Lift it out of the
//! tetrahedron and it becomes a **global** name: the edge's two grid points and
//! the root's ordinal.
//!
//! The direction is the part that has to hold, and it does, for a reason that
//! does not depend on the cell. `TETS[t]` orders its corners by *inclusion*, so
//! a tet edge runs from the corner whose offset bits are a subset to the one
//! whose bits are a superset. Offsets are the grid point minus the cell origin,
//! and subtracting the same origin from both preserves the componentwise
//! comparison — so inclusion is equivalent to `P ≤ Q` componentwise, which names
//! the same direction from whichever cell the edge is viewed. Two tetrahedra in
//! two cells therefore derive the same key, and — since their corner positions
//! are bit-identical and [`all_roots`] is deterministic — the same position for
//! it.
//!
//! **What this buys is not speed.** It is that a vertex has an *identity* rather
//! than merely a location, so a later stage can move one and have every triangle
//! standing on it move too. Welding by position cannot do that: after the move
//! the copies no longer coincide, which is precisely the tearing M-101 measured
//! twice and M-162 traced to 150 triangles in neighbouring cells.
//!
//! Steiner points are deliberately **not** shared. A centroid is a property of
//! one tetrahedron, and giving it a key would assert an identity it does not
//! have.
//!
//! # Why every cell recomputes its neighbours' edges
//!
//! A grid edge shared by two cells has its roots found twice, and a tetrahedron
//! edge shared inside a cell likewise. That is deliberate for now and it is what
//! makes the result correct without a cache: both calls pass **bit-identical**
//! endpoints, because a corner's position is always `origin + cell_size · index`
//! computed the same way, and [`super::roots::all_roots`] is
//! deterministic for identical arguments. A cache keyed on the grid edge would
//! be faster and is the obvious optimisation, but it is an optimisation with a
//! correctness precondition, and this ticket owes a working extractor before a
//! fast one.

use alloc::collections::BTreeMap;
use alloc::vec::Vec;

use crate::cube::corner_offset;
use crate::marching_tetrahedra::table::{TET_EDGE_COUNT, TET_EDGES, TETS};
use crate::mesh::MeshSink;
use crate::real::Real;
use crate::sdf::Sdf;
use crate::shape::Shape3;
use crate::vec3;

use super::curves::FacePoint;
use super::roots::all_roots;
use super::surface::{TetCrossings, TetPatch, Unfilled, fill};

/// A crossing's identity, independent of which tetrahedron found it.
///
/// Two rules, because a crossing sits in one of two structurally different
/// places and only one of them is named by an edge. See this module's docs for
/// why the edge's direction is the same from either cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum CrossingKey {
    /// A root in the interior of a tetrahedron edge: the edge's two grid points,
    /// componentwise smaller first, and which root along it this is.
    OnEdge {
        lo: [u32; 3],
        hi: [u32; 3],
        index: u32,
    },
    /// A root that lies **on** a grid sample point, named by that point.
    ///
    /// Up to 24 tetrahedron edges meet at a grid point, so such a root has a
    /// different `(edge, index)` on every one of them — one point wearing many
    /// names, which no correct sharing under [`OnEdge`] can merge (M-169). A grid
    /// point has a single coordinate triple and needs no ordinal, which makes
    /// this rule *simpler* than the edge one rather than a special case bolted
    /// onto it: it needs no inclusion argument at all, because an absolute grid
    /// coordinate is cell-independent outright.
    ///
    /// [`OnEdge`]: CrossingKey::OnEdge
    OnGridPoint { at: [u32; 3] },
}

/// Which of a tetrahedron's four corners the field is exactly zero at.
///
/// The **exact** test — no tolerance — because it is the definition of the
/// surface passing through a grid point, and because it must give the same
/// answer in every cell that shares that point. `f(G)` depends on `G` alone,
/// where anything derived from an edge does not.
///
/// Note this samples the corner directly rather than reusing
/// [`all_roots`]'s own endpoint evaluations: that function evaluates `t = 1` as
/// `a + (b − a)·1`, which is *not* bit-identical to `b`, so its answer there is a
/// property of the edge and not of the grid point (M-183).
fn corners_on_surface<R: Real, S: Sdf<Scalar = R>>(sdf: &S, corners: &[[R; 3]; 4]) -> [bool; 4] {
    // `partial_cmp` rather than `== 0`: `Some(Equal)` accepts both zeros, which
    // the sign convention already treats alike, and `None` rejects NaN, which a
    // bare equality would silently answer `false` to without saying why.
    corners.map(|c| {
        matches!(
            sdf.sample(c).partial_cmp(&R::ZERO),
            Some(core::cmp::Ordering::Equal)
        )
    })
}

impl CrossingKey {
    /// The two grid points of a tetrahedron edge, in the order `TETS` gives.
    fn edge_grid_points(cell: [u32; 3], tet: usize, edge: u8) -> [[u32; 3]; 2] {
        let [a, b] = TET_EDGES[edge as usize];
        [a, b].map(|corner| {
            let offset = corner_offset(TETS[tet][corner as usize]);
            [
                cell[0] + offset[0],
                cell[1] + offset[1],
                cell[2] + offset[2],
            ]
        })
    }

    /// Lift a tetrahedron-local crossing to its global name.
    ///
    /// `endpoint` names the edge corner this root sits *on*, when it sits on one.
    fn of(cell: [u32; 3], tet: usize, point: FacePoint, endpoint: Option<u8>) -> Self {
        let [lo, hi] = Self::edge_grid_points(cell, tet, point.edge);
        match endpoint {
            Some(0) => Self::OnGridPoint { at: lo },
            Some(_) => Self::OnGridPoint { at: hi },
            None => Self::OnEdge {
                lo,
                hi,
                index: point.index,
            },
        }
    }
}

/// Relative floor for calling a gradient too small to be trusted.
///
/// Scaled by the tetrahedron's own slope — its corner-value spread over the cell
/// size — so it means the same thing at every field magnitude and every grid
/// spacing. Same reasoning as
/// [`ValidateConfig`](crate::validate::ValidateConfig)'s thresholds: an absolute
/// floor on a gradient is meaningless without a scale, and picks a different
/// answer on the same shape at a different resolution.
const NORMAL_CONDITION_REL: f64 = 1e-6;

/// Why a tetrahedron could not be given normals.
///
/// **These are different problems and the report keeps them apart.** Folding
/// them into one count would let an ill-conditioning bug hide inside a
/// degeneracy total, which is the shape of failure this crate has caught
/// repeatedly (M-279's rule: a falsifier must separate the hypothesis from its
/// rivals).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NormalCause {
    /// The gradient is **exactly zero**: a critical point of the field.
    ///
    /// If the surface also passes through that point, the level set is singular
    /// there and **no normal exists** — not a missing one, an absent one. No
    /// arithmetic recovers it, and widening the stencil would return the
    /// gradient of a *smoothed* field, which is a different field.
    Degenerate,
    /// The gradient is non-zero but below the conditioning floor.
    ///
    /// A normal exists and is not trustworthy: its direction is dominated by
    /// rounding. Distinct from [`Degenerate`](Self::Degenerate) because the
    /// remedy is different — this one is about precision, that one about
    /// topology.
    IllConditioned,
}

/// One place the extractor declined to emit, and what was true there.
///
/// **Positions, not just a count.** A count can stay at 33 while the cells move,
/// so a count alone is not a regression test; and a caller repairing the mesh
/// needs to know *where*.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalSite<R> {
    /// The grid cell.
    pub cell: [u32; 3],
    /// Which of the six tetrahedra of that cell.
    pub tet: u8,
    /// Where the normal was wanted, in world coordinates.
    pub position: [R; 3],
    /// `|∇f|` there. Exactly zero for [`NormalCause::Degenerate`].
    pub gradient_length: R,
    /// Whether the position is a **cell corner the surface passes exactly
    /// through**.
    ///
    /// The diagnostic signal, and on quantised data it is usually `true`: `u8`
    /// samples against an integer isovalue land *on* the surface constantly —
    /// **3% of surface-cell corners** on `bonsai` (M-316). Contouring at a
    /// half-offset isovalue removes that case entirely; see
    /// [`SubgridMarchingTetrahedra`]'s docs.
    pub on_surface_corner: bool,
    /// Which problem it was.
    pub cause: NormalCause,
}

/// What one extraction could not do, and where.
///
/// **Recorded, not asserted** — the same standing `MeshReport`'s metrics have,
/// and for the same reason: a non-zero here can be a property of the *input*
/// rather than a defect. A field with a critical point on its isosurface has no
/// normal there, and refusing to invent one is correct.
///
/// **It is not silence, either.** The mesh has a hole wherever a site is listed,
/// and the count comes back with the mesh rather than going to a log, so a caller
/// cannot fail to see it.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct NormalReport<R> {
    /// Tetrahedra skipped for want of a normal. Each leaves a hole.
    pub skipped_tetrahedra: u64,
    /// Of those, [`NormalCause::Degenerate`].
    pub degenerate: u64,
    /// Of those, [`NormalCause::IllConditioned`].
    pub ill_conditioned: u64,
    /// Where each was. One entry per skipped tetrahedron.
    pub sites: Vec<NormalSite<R>>,
}

impl<R> NormalReport<R> {
    /// Nothing was skipped: every tetrahedron got its normals.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.skipped_tetrahedra == 0
    }

    fn clear(&mut self) {
        self.skipped_tetrahedra = 0;
        self.degenerate = 0;
        self.ill_conditioned = 0;
        self.sites.clear();
    }

    fn record(&mut self, site: NormalSite<R>) {
        self.skipped_tetrahedra += 1;
        match site.cause {
            NormalCause::Degenerate => self.degenerate += 1,
            NormalCause::IllConditioned => self.ill_conditioned += 1,
        }
        self.sites.push(site);
    }
}

/// A patch vertex, resolved to what emission needs.
#[derive(Clone, Copy, Debug)]
enum Prepared<R> {
    /// Already emitted by another tetrahedron; reuse its index.
    Shared(u32),
    /// Needs emitting, with this position and normal.
    Fresh {
        position: [R; 3],
        normal: [R; 3],
        key: Option<CrossingKey>,
    },
    /// The same vertex as an earlier entry **in this tetrahedron**.
    ///
    /// A-014h's case: where the surface passes exactly through a grid point, the
    /// two tetrahedron edges meeting there each carry a root *at* it, and naming
    /// both by the point makes them one vertex. The shared table cannot serve
    /// this — it is not written until the tetrahedron is known to be emittable —
    /// so the repeat is resolved here instead. **Losing this collapses the two
    /// into distinct vertices, and the zero-area sliver between them stops
    /// looking degenerate and gets emitted**, which is the index buffer M-185
    /// says this crate's own validator calls invalid.
    Repeat(usize),
}

/// Subgrid Marching Tetrahedra — Baktash, Gillespie & Crane,
/// `10.48550/arXiv.2606.00454`.
///
/// Holds its working buffers so a repeated extraction allocates nothing new,
/// per `CLAUDE.md` rule 6.
///
/// # Normals, and the one input this cannot mesh (A-028)
///
/// Each vertex takes the field's own gradient as its normal. Where the gradient
/// **vanishes**, no normal exists — and if the isosurface passes through that
/// point, the level set is genuinely singular there rather than merely awkward.
/// Such a tetrahedron is **skipped**, leaving a hole its size, and
/// [`report`](Self::report) says how many and where. Nothing is substituted: a
/// wider stencil would return the gradient of a *smoothed* field, which is a
/// different field, and at a saddle there may be no correct normal at all.
///
/// # If your data is integer, contour at a half-offset isovalue
///
/// **This removes the problem rather than reporting it, and costs one line.**
/// `u8` or `u16` samples against an **integer** isovalue land *exactly* on the
/// surface constantly — **3% of surface-cell corners** on the `bonsai` CT volume
/// (M-316) — and a corner the surface passes exactly through is where this
/// extractor asks for a normal at a grid point, which is where the degeneracy
/// lives.
///
/// Contour at `127.5` rather than `127` and **no sample can sit on the
/// isosurface at all**, because a half-integer is not attainable by integer
/// data. Standard practice in volume rendering, for this reason.
///
/// ```text
/// let iso = 127.5;                       // not 127
/// let values: Vec<f64> = bytes.iter().map(|b| iso - f64::from(*b)).collect();
/// ```
///
/// It does **not** remove ill-conditioning near a critical point, and the
/// critical points are still in the field — it removes the *exact* zero that
/// this extractor trips over, which is the case that actually fires.
#[derive(Clone, Debug)]
pub struct SubgridMarchingTetrahedra<R: Real> {
    samples: u32,
    /// The field's Lipschitz constant, when the caller knows one.
    ///
    /// Enables the empty-cell rejection described on
    /// [`set_lipschitz`](Self::set_lipschitz). `None` means unknown, and every
    /// cell is subdivided.
    lipschitz: Option<R>,
    along: [Vec<R>; TET_EDGE_COUNT],
    patch: TetPatch<R>,
    index: Vec<u32>,
    /// Which sink vertex each crossing was emitted as, for the extraction in
    /// progress. Cleared per [`extract`](Self::extract), so it never carries a
    /// stale index from a previous grid.
    shared: BTreeMap<CrossingKey, u32>,
    /// Per-vertex resolution for the tetrahedron in progress. A field rather
    /// than a local so a repeated extraction allocates nothing (rule 6).
    prepared: Vec<Prepared<R>>,
    /// What the extraction in progress could not do. Cleared per
    /// [`extract`](Self::extract).
    report: NormalReport<R>,
}

impl<R: Real> SubgridMarchingTetrahedra<R> {
    /// A new extractor sampling each tetrahedron edge `samples` times.
    ///
    /// `samples` is the 1D marching resolution, and it is the knob that decides
    /// which features exist: a pair of crossings closer together than
    /// `1 / samples` of an edge is invisible, exactly as a pair closer than the
    /// grid spacing is invisible to a sign-based method (§1.3). It is **not** a
    /// quality setting on the same axis as grid resolution — raising it resolves
    /// thinner features at the same triangle count, which is the entire point of
    /// the method.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidCellSize`](crate::Error::InvalidCellSize) if `samples` is
    /// zero, which would make every edge unsampled and every mesh empty.
    pub fn new(samples: u32) -> crate::Result<Self> {
        if samples == 0 {
            return Err(crate::Error::InvalidCellSize { value: 0.0 });
        }
        Ok(Self {
            samples,
            lipschitz: None,
            along: Default::default(),
            patch: TetPatch::new(),
            index: Vec::new(),
            shared: BTreeMap::new(),
            prepared: Vec::new(),
            report: NormalReport::default(),
        })
    }

    /// The 1D sampling resolution this extractor was built with.
    #[must_use]
    pub const fn samples(&self) -> u32 {
        self.samples
    }
    /// Tell the extractor the field's Lipschitz constant, enabling empty-cell
    /// rejection.
    ///
    /// Ticket: F-005. Source: Hart, *Sphere tracing* (`10.1007/s003710050084`).
    ///
    /// # What it buys, and why this extractor in particular
    ///
    /// M-98 measured this extractor at **70× the cost of Marching Cubes**, and
    /// found the constant is the whole story: `6 tets × 6 edges × 16 samples =
    /// 576` field evaluations per cell against Marching Cubes' 8. The prediction
    /// was ~72× and the measurement was 70×, so there is nothing else to blame.
    ///
    /// A field with Lipschitz constant `l` cannot change by more than `l · d`
    /// over a distance `d`. So if `|f(centre)| > l · r`, where `r` is the
    /// circumradius of the cell — half its diagonal, `(√3/2)·h` — then `f` cannot
    /// reach zero anywhere inside it, and **one evaluation replaces 576**.
    ///
    /// # It cannot change the output, and that is asserted rather than argued
    ///
    /// This is a short-circuit, not an alternative algorithm: a rejected cell is
    /// one the full path would have found empty. `rejection_does_not_change_the_mesh`
    /// runs both ways on all eight reference fields and compares **bit for bit**,
    /// so the claim is checked rather than reasoned about. That is also why this
    /// does not violate the one-path rule — there is one path, and a proof that
    /// part of it need not be walked.
    ///
    /// # Pass the constant, not a guess
    ///
    /// [`FieldBound::lipschitz`](crate::fields::FieldBound::lipschitz) is where
    /// a reference field's constant comes from. **A value smaller than the
    /// field's true constant will reject cells that contain surface**, which is
    /// a hole in the mesh rather than a slow path — so `None`, the default, is
    /// the only safe answer when it is not known. M-244 is the incident: a
    /// hand-reasoned constant was wrong by 3× on the first try.
    /// What the last [`extract`](Self::extract) could not do, and where.
    ///
    /// Empty after a clean run. A non-zero
    /// [`skipped_tetrahedra`](NormalReport::skipped_tetrahedra) means the mesh
    /// has that many holes, each one tetrahedron in size, and
    /// [`sites`](NormalReport::sites) says where.
    ///
    /// **Read it.** A caller that ignores this gets a mesh with unannounced
    /// holes, which is exactly the silence the report exists to prevent.
    #[must_use]
    pub fn report(&self) -> &NormalReport<R> {
        &self.report
    }

    /// Tell the extractor the field's Lipschitz constant, enabling empty-cell
    /// rejection.
    ///
    /// Ticket: F-005. Source: Hart, *Sphere tracing* (`10.1007/s003710050084`).
    ///
    /// # What it buys, and why this extractor in particular
    ///
    /// M-98 measured this extractor at **70× the cost of Marching Cubes**, and
    /// found the constant is the whole story: `6 tets × 6 edges × 16 samples =
    /// 576` field evaluations per cell against Marching Cubes' 8. The prediction
    /// was ~72× and the measurement was 70×, so there is nothing else to blame.
    ///
    /// A field with Lipschitz constant `l` cannot change by more than `l · d`
    /// over a distance `d`. So if `|f(centre)| > l · r`, where `r` is the
    /// circumradius of the cell — half its diagonal, `(√3/2)·h` — then `f` cannot
    /// reach zero anywhere inside it, and **one evaluation replaces 576**.
    ///
    /// # It cannot change the output, and that is asserted rather than argued
    ///
    /// This is a short-circuit, not an alternative algorithm: a rejected cell is
    /// one the full path would have found empty. `rejection_does_not_change_the_mesh`
    /// runs both ways on all eight reference fields and compares **bit for bit**,
    /// so the claim is checked rather than reasoned about. That is also why this
    /// does not violate the one-path rule — there is one path, and a proof that
    /// part of it need not be walked.
    ///
    /// # Pass the constant, not a guess
    ///
    /// [`FieldBound::lipschitz`](crate::fields::FieldBound::lipschitz) is where
    /// a reference field's constant comes from. **A value smaller than the
    /// field's true constant will reject cells that contain surface**, which is
    /// a hole in the mesh rather than a slow path — so `None`, the default, is
    /// the only safe answer when it is not known. M-244 is the incident: a
    /// hand-reasoned constant was wrong by 3× on the first try.
    pub fn set_lipschitz(&mut self, l: Option<R>) {
        self.lipschitz = l;
    }

    /// Can one evaluation prove this cell contains no surface?
    ///
    /// `false` whenever no constant is known, which is the default — see
    /// [`set_lipschitz`](Self::set_lipschitz).
    fn cell_is_provably_empty<S>(
        &self,
        sdf: &S,
        origin: [R; 3],
        cell_size: R,
        cell: [u32; 3],
    ) -> bool
    where
        S: Sdf<Scalar = R>,
    {
        let Some(l) = self.lipschitz else {
            return false;
        };
        let half = cell_size * R::HALF;
        let centre = [
            origin[0] + cell_size * R::from_f64(f64::from(cell[0])) + half,
            origin[1] + cell_size * R::from_f64(f64::from(cell[1])) + half,
            origin[2] + cell_size * R::from_f64(f64::from(cell[2])) + half,
        ];
        // Circumradius of the cell: half the space diagonal.
        let radius = half * R::from_f64(1.732_050_807_568_877_2);
        // Strict, so a value exactly on the bound subdivides. The surface can
        // touch a corner at equality, and a hole is worse than a wasted cell.
        sdf.sample(centre).abs() > l * radius
    }

    /// Extract the zero level set into `out`.
    ///
    /// `shape` counts **samples**, so `[n; 3]` spans `n - 1` cells per axis, and
    /// `origin` is the world position of sample `[0, 0, 0]`. Note that unlike
    /// every other extractor here, the grid is used only for its *geometry* —
    /// the field is never sampled at the grid nodes, because node values are
    /// exactly the information this method replaces.
    ///
    /// # Winding, and why the output must be welded before it is judged
    ///
    /// Counter-clockwise seen from outside the solid, matching every other
    /// extractor here. That is **imposed rather than inherited**: §3.2 fixes
    /// each polygon's vertex order from its own boundary curve, which is
    /// consistent within a tetrahedron and carries no relationship to which side
    /// the field calls inside. Each triangle is therefore flipped, if needed, to
    /// agree with the gradient at its own centroid — per triangle and not per
    /// patch, because a sheet thinner than a cell puts two oppositely-facing
    /// surfaces inside one tetrahedron.
    ///
    /// # Vertices *are* shared, and welding this output is unnecessary and can
    /// damage it
    ///
    /// This said the opposite until A-018 measured it: *"vertices are emitted per
    /// tetrahedron and are not shared … before welding, the output is a triangle
    /// soup"*. That was true when M-93 and M-96 were written and **A-014h ended
    /// it** — the extractor gives every crossing a global identity and emits it
    /// once, so the raw output already has `boundary_edges == 0` on every closed
    /// reference field, with no weld at all.
    ///
    /// A positional weld afterwards is therefore a no-op at best. On seven of the
    /// eight reference fields it is exactly that — same vertex count, same
    /// topology. On `noise_cavity` it merges **one** pair of vertices that are
    /// coincident *by position* and distinct *by identity*, fusing two sheets and
    /// **adding two non-manifold edges and three non-manifold vertices**
    /// (M-226). Sharing by identity is strictly finer than sharing by position,
    /// and once it is complete the coarser rule has nothing left to contribute
    /// except mistakes.
    ///
    /// So do not weld this output to make it a surface; it already is one. Weld
    /// only to *join it to other geometry*, and expect that to cost the pair above
    /// wherever two sheets pass within the tolerance of each other.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if any axis has fewer
    /// than two samples. [`Error::SubgridUnfilled`](crate::Error::SubgridUnfilled)
    /// if a tetrahedron could not be triangulated — a defect rather than an
    /// unsupported input, since every case §3.2 defines is implemented.
    /// [`Error::DegenerateNormal`](crate::Error::DegenerateNormal) if the field's
    /// gradient at a crossing is zero or non-finite.
    /// [`Error::IndexSpaceExhausted`](crate::Error::IndexSpaceExhausted) if the
    /// extraction reaches `u32::MAX` emitted vertices — checked as it emits
    /// rather than up front, because an edge can carry any number of roots, so
    /// no a-priori bound exists.
    pub fn extract<S, M>(
        &mut self,
        sdf: &S,
        shape: &impl Shape3,
        origin: [R; 3],
        cell_size: R,
        out: &mut M,
    ) -> crate::Result<()>
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        let size = shape.size();
        if size[0] < 2 || size[1] < 2 || size[2] < 2 {
            return Err(crate::Error::GridTooSmall { size });
        }
        // Indices in here name vertices of `out`, so a table held across two
        // extractions would hand the second one the first's numbering.
        self.shared.clear();
        self.report.clear();

        // Unlike every other extractor, the vertex count here has no a-priori
        // bound -- an edge can carry any number of roots -- so the u32 index
        // space is guarded at emission rather than up front.
        let mut vertices: u64 = 0;

        for z in 0..size[2] - 1 {
            for y in 0..size[1] - 1 {
                for x in 0..size[0] - 1 {
                    let cell = [x, y, z];
                    if self.cell_is_provably_empty(sdf, origin, cell_size, cell) {
                        continue;
                    }
                    for t in 0..TETS.len() {
                        self.cell_tet(sdf, origin, cell_size, cell, t, &mut vertices, out)?;
                    }
                }
            }
        }
        Ok(())
    }

    /// Which tetrahedron-edge corner this crossing sits on, if it sits on one.
    ///
    /// `0` for the edge's lower corner, `1` for its upper, in `TET_EDGES` order.
    ///
    /// # Why the ordinal decides it rather than the parameter
    ///
    /// The obvious test — is the root's parameter 0 or 1 — describes nothing:
    /// **no root anywhere reports `t == 0`** and almost none reports `t == 1`
    /// (M-179), because [`refine`](super::roots) returns the *upper* end of its
    /// final bracket so it can keep the ascending-and-distinct contract when a
    /// root sits on a sample.
    ///
    /// What is exact is `f(G) == 0`, and `all_roots` reports at most one root per
    /// bracketing interval in ascending order. So if the surface passes through
    /// the lower corner and a root exists in the *first* interval, that root is
    /// the one at the corner, and it is `index == 0`. Symmetrically at the top.
    /// Both facts are properties of the grid point and the edge's own sampling,
    /// so every tetrahedron meeting there reaches the same answer.
    fn endpoint_of(&self, point: FacePoint, on_surface: &[bool; 4]) -> Option<u8> {
        let [a, b] = TET_EDGES[point.edge as usize];
        let roots = &self.along[point.edge as usize];
        let count = roots.len() as u32;
        if count == 0 {
            return None;
        }
        let step = R::ONE / R::from_f64(f64::from(self.samples));
        let last_interval_start = R::from_f64(f64::from(self.samples - 1)) * step;
        let t = *roots.get(point.index as usize)?;

        let at_lo = on_surface[a as usize] && point.index == 0 && t <= step;
        let at_hi = on_surface[b as usize] && point.index == count - 1 && t >= last_interval_start;
        match (at_lo, at_hi) {
            // Both ends on the surface with a single root between them: the
            // parameter is the only thing left that can say which one it is.
            (true, true) => Some(u8::from(t + t >= R::ONE)),
            (true, false) => Some(0),
            (false, true) => Some(1),
            (false, false) => None,
        }
    }

    /// One tetrahedron of one cell.
    #[allow(clippy::too_many_arguments)]
    fn cell_tet<S, M>(
        &mut self,
        sdf: &S,
        origin: [R; 3],
        cell_size: R,
        cell: [u32; 3],
        t: usize,
        vertices: &mut u64,
        out: &mut M,
    ) -> crate::Result<()>
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        // Corner positions. Written as `origin + cell_size · index` with the
        // index formed once, so the two cells sharing a face compute the same
        // world position by the same expression — M-32's caveat is that equal
        // by algebra is not equal by IEEE, and only the same expression is safe.
        let mut corners = [[R::ZERO; 3]; 4];
        for (c, slot) in corners.iter_mut().enumerate() {
            let offset = corner_offset(TETS[t][c]);
            for axis in 0..3 {
                let index = f64::from(cell[axis]) + f64::from(offset[axis]);
                slot[axis] = origin[axis] + cell_size * R::from_f64(index);
            }
        }

        // Every zero along every edge. `TETS[t]` is ordered by inclusion, so a
        // tet edge always runs from the lower cube-corner index to the higher —
        // which makes the traversal direction a property of the grid rather than
        // of the tetrahedron, and is what lets two tetrahedra sharing an edge
        // agree bit-for-bit without consulting each other.
        let mut total = 0usize;
        for (e, slot) in self.along.iter_mut().enumerate() {
            slot.clear();
            let [lo, hi] = TET_EDGES[e];
            all_roots(
                corners[lo as usize],
                corners[hi as usize],
                sdf,
                self.samples,
                slot,
            );
            total += slot.len();
        }
        if total == 0 {
            return Ok(());
        }

        let mut borrowed: [&[R]; TET_EDGE_COUNT] = [&[]; TET_EDGE_COUNT];
        for (slot, v) in borrowed.iter_mut().zip(self.along.iter()) {
            *slot = v.as_slice();
        }
        let crossings = TetCrossings {
            corners,
            along: borrowed,
        };

        let unfilled = fill(&crossings, &mut self.patch).map_err(|_| {
            // `check` can only fail on unsorted or out-of-range parameters, and
            // `all_roots` produces neither — so this is this crate's bug, and it
            // is reported as one rather than swallowed.
            crate::Error::SubgridUnfilled {
                cell,
                tet: t as u8,
                reason: "crossings rejected as malformed",
            }
        })?;
        if unfilled != Unfilled::None {
            return Err(crate::Error::SubgridUnfilled {
                cell,
                tet: t as u8,
                reason: match unfilled {
                    Unfilled::SingleLoop => "single loop",
                    Unfilled::Subdivision => "subdivision",
                    Unfilled::NonNormalLoop => "non-normal loop",
                    Unfilled::NoPattern => "residual is not a (d1, d2) pattern",
                    Unfilled::Inconsistent => "curves disagree with the crossings",
                    Unfilled::None => unreachable!(),
                },
            });
        }

        // Which corners the surface passes exactly through. A root adjacent to
        // one of those is *at* it, and is named and placed by the grid point
        // rather than by the edge it was found along (A-014h).
        let on_surface = corners_on_surface(sdf, &corners);

        // ── pass one: resolve every vertex before emitting any of them ─────
        //
        // **The order is what makes a skip clean.** If a normal is unavailable
        // the whole tetrahedron is declined, and declining after emitting some
        // of its vertices would leave them in the sink with no triangle
        // referencing them — vertices `validate_indexed` counts as unreferenced
        // and a caller cannot explain. Nothing reaches the sink until every
        // vertex of this tetrahedron is known to be emittable.
        // The conditioning floor, relative to this tetrahedron's own slope. An
        // absolute floor would mean a different thing at every field magnitude
        // and grid spacing. Computed once per tetrahedron: it is a property of
        // the tetrahedron, not of the vertex being resolved.
        let mut lo = R::ZERO;
        let mut hi = R::ZERO;
        for (k, corner) in corners.iter().enumerate() {
            let v = sdf.sample(*corner);
            if k == 0 || v < lo {
                lo = v;
            }
            if k == 0 || v > hi {
                hi = v;
            }
        }
        let floor = R::from_f64(NORMAL_CONDITION_REL) * ((hi - lo) / cell_size);

        let mut prepared = core::mem::take(&mut self.prepared);
        prepared.clear();
        prepared.reserve(self.patch.positions.len());
        let mut refused: Option<NormalSite<R>> = None;

        for (at, position) in self.patch.positions.iter().enumerate() {
            // A crossing carries a global name; a Steiner point does not, and
            // gets a fresh vertex because it is shared with nothing.
            let crossing = self.patch.crossings.get(at).copied();
            let endpoint = crossing.and_then(|point| self.endpoint_of(point, &on_surface));
            let key = crossing.map(|point| CrossingKey::of(cell, t, point, endpoint));
            if let Some(existing) = key.and_then(|k| self.shared.get(&k)) {
                prepared.push(Prepared::Shared(*existing));
                continue;
            }
            // Already resolved earlier in *this* tetrahedron? The shared table
            // is not written until the tetrahedron survives, so the within-tet
            // repeat has to be found here. The patch has a handful of vertices,
            // so a scan is cheaper than a second map.
            if let Some(k) = key {
                let earlier = prepared.iter().position(
                    |p| matches!(p, Prepared::Fresh { key: Some(other), .. } if *other == k),
                );
                if let Some(index) = earlier {
                    prepared.push(Prepared::Repeat(index));
                    continue;
                }
            }

            // A root on a grid point is emitted **at** that point, by the same
            // `origin + cell_size · index` expression every cell computes
            // identically (M-32) -- not at the lerp's answer, which differs in
            // the last bits between the edges that meet there and is what left
            // 690 of `box_exact`'s duplicates unmergeable (M-180).
            let position = &match (crossing, endpoint) {
                (Some(point), Some(end)) => {
                    let [a, b] = TET_EDGES[point.edge as usize];
                    corners[if end == 0 { a } else { b } as usize]
                }
                _ => *position,
            };

            let g = sdf.gradient(*position);
            let length = vec3::length(g);

            // `!is_finite() || <= 0` rather than a negated `>`: NaN is excluded by
            // the finiteness test first, so the comparison is total by the time
            // it runs.
            let cause = if !length.is_finite() || length <= R::ZERO {
                Some(NormalCause::Degenerate)
            } else if floor > R::ZERO && length < floor {
                Some(NormalCause::IllConditioned)
            } else {
                None
            };
            if let Some(cause) = cause {
                refused = Some(NormalSite {
                    cell,
                    tet: t as u8,
                    position: *position,
                    gradient_length: length,
                    on_surface_corner: endpoint.is_some(),
                    cause,
                });
                break;
            }

            prepared.push(Prepared::Fresh {
                position: *position,
                normal: vec3::scale(g, length.recip()),
                key,
            });
        }

        // ── the decision, before anything reaches the sink ─────────────────
        if let Some(site) = refused {
            self.report.record(site);
            self.prepared = prepared;
            // Skipping this tetrahedron leaves a hole exactly its size, and the
            // report says where. **Not a fallback**: nothing is substituted for
            // the geometry, and the omission is returned with the mesh rather
            // than logged.
            return Ok(());
        }

        // ── pass two: emit ─────────────────────────────────────────────────
        self.index.clear();
        self.index.reserve(prepared.len());
        for step in &prepared {
            match step {
                Prepared::Shared(existing) => self.index.push(*existing),
                Prepared::Repeat(earlier) => {
                    let Some(index) = self.index.get(*earlier).copied() else {
                        self.prepared = prepared;
                        return Err(crate::Error::SubgridUnfilled {
                            cell,
                            tet: t as u8,
                            reason: "a repeated vertex referred to one not yet emitted",
                        });
                    };
                    self.index.push(index);
                }
                Prepared::Fresh {
                    position,
                    normal,
                    key,
                } => {
                    // The sink's index space is u32 and `vertex` has no way to
                    // report exhaustion, so the count is enforced here: after
                    // `u32::MAX` emissions the next index a non-welding sink
                    // hands back would be `u32::MAX` itself, one past the last
                    // addressable vertex.
                    if *vertices >= u64::from(u32::MAX) {
                        self.prepared = prepared;
                        return Err(crate::Error::IndexSpaceExhausted {
                            needed: *vertices + 1,
                        });
                    }
                    let emitted = out.vertex(*position, *normal);
                    *vertices += 1;
                    if let Some(key) = key {
                        self.shared.insert(*key, emitted);
                    }
                    self.index.push(emitted);
                }
            }
        }
        self.prepared = prepared;
        for tri in &self.patch.triangles {
            let (a, b, c) = (
                self.index.get(tri[0] as usize),
                self.index.get(tri[1] as usize),
                self.index.get(tri[2] as usize),
            );
            let (Some(a), Some(b), Some(c)) = (a, b, c) else {
                return Err(crate::Error::SubgridUnfilled {
                    cell,
                    tet: t as u8,
                    reason: "a triangle indexed a vertex the patch does not have",
                });
            };

            // Two of this triangle's corners are the same vertex, so it has no
            // area and is not a triangle. It arises only from A-014h: when the
            // surface passes exactly through a grid point, the two tetrahedron
            // edges meeting there each carry a root *at* it, and naming both by
            // the point makes them one vertex. The sliver between them was always
            // zero-area; it is now zero-area and says so in its indices, which
            // `validate_indexed` counts as a structural error rather than as a
            // recorded metric. Declining to emit it is not dropping geometry —
            // there is no geometry to drop — and the alternative is shipping an
            // index buffer this crate's own validator calls invalid (M-185).
            if a == b || b == c || a == c {
                continue;
            }

            // Orientation. §3.2 fixes each polygon's vertex order from its own
            // boundary curve, which is consistent within a tetrahedron and
            // carries no relation to which side the field calls inside — so the
            // winding has to be imposed here, against the only thing that knows:
            // the gradient, which points away from the solid.
            //
            // Per triangle rather than per patch, because one tetrahedron can
            // carry sheets facing opposite ways — `thin_plate`'s two faces are
            // 0.4 cells apart and routinely land in the same cell.
            let (pa, pb, pc) = (
                self.patch.positions[tri[0] as usize],
                self.patch.positions[tri[1] as usize],
                self.patch.positions[tri[2] as usize],
            );
            let face = vec3::cross(vec3::sub(pb, pa), vec3::sub(pc, pa));
            let third = R::ONE / R::from_f64(3.0);
            let centroid = [
                (pa[0] + pb[0] + pc[0]) * third,
                (pa[1] + pb[1] + pc[1]) * third,
                (pa[2] + pb[2] + pc[2]) * third,
            ];
            let outward = vec3::dot(face, sdf.gradient(centroid));

            if outward < R::ZERO {
                out.triangle(*a, *c, *b);
            } else {
                // Includes the exactly-zero case, which is a triangle with no
                // area — §3.2's boundary disks emit those by construction
                // (V-21) and there is no orientation to choose for one. Left in
                // its original order rather than dropped, so the connectivity
                // §3.2 built stays intact.
                out.triangle(*a, *b, *c);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
