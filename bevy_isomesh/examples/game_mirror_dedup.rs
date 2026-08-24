//! E-314 — you mirrored the boulder to reuse it, and the mesh cache stopped
//! hitting.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_mirror_dedup --release
//! ```
//!
//! **Always `--release`.** Startup runs 8 fields × 49 extractions of a 33³ grid
//! to reproduce every number this demo quotes, and a debug build meshes 20-50x
//! slower.
//!
//! It runs itself and loops through all 48 elements of the octahedral group.
//! `Right`/`Left` step the element by hand, `A` hands it back to the tour,
//! `H` toggles the marker on every vertex that moved, `Space` freezes the clock,
//! `W` puts the wireframe over the rock and `N` draws the normals.
//! Digits `1`-`7` pick the field; `ISOMESH_FIELD=0..7` pins one for a still
//! without a keyboard, and `7` (`fbm_terrain`) is only reachable that way
//! because the shared harness only binds seven digits.
//!
//! ```bash
//! ISOMESH_VIEW=nogrid ISOMESH_CAPTURE_FRAMES=120 ISOMESH_CAPTURE_EVERY=2 \
//!   FPS=10 WIDTH=1000 DISPLAY=:77 \
//!   ./scripts/record_gif.sh game_mirror_dedup docs/gifs/mirrored-is-not-the-same-mesh.gif
//! ```
//!
//! Demonstrates **✗39 / M-356 / P-57** (`docs/experiments/p-57.csv`) as a studio
//! meets it: not as a pass rate over a swept fixture, but as two chunks that
//! render identically and hash differently.
//!
//! # The finding, in one paragraph
//!
//! There are exactly **48** ways to reorient a chunk on the cubic lattice — the
//! signed coordinate permutations, three axes reordered and each optionally
//! negated. Every one of them is exact in `f64`: a permutation moves bits and a
//! negation flips one, and neither rounds. So `mesh(g·f)` and `g·mesh(f)` ought
//! to be the same `f64` mesh for all 48. Measured on `marching_cubes`, **six
//! are**, and they are precisely the six that **never negate a coordinate**.
//! [`isomesh::marching_cubes::table`]'s `EDGE_CORNERS` orients every grid edge
//! along *increasing grid index*, so negating an axis swaps which endpoint is
//! `lo`: the extractor computes `b/(b−a)` where it computed `a/(a−b)` and
//! anchors the lerp at the other corner. Those two are `1 − t` and `t` of the
//! same edge — the same point to infinite precision, and **not bit-reciprocal**.
//!
//! # What a studio meets
//!
//! Mirroring an asset to reuse it is the oldest trick in level art. Anything
//! keyed on a **content hash** — GPU instancing, a collision-mesh cache, a chunk
//! mesh cache, a network delta, an asset-dedup pass in a cooker — silently
//! misses on the mirrored copy and hits on some of the rotated ones. The two
//! meshes are visually identical, so nothing looks wrong; you just pay twice for
//! one rock, and the cache-hit graph in your telemetry is worse than it should
//! be for a reason nobody can name.
//!
//! # It is not handedness, and the ledger's own columns say so
//!
//! The obvious story is "reflections break it, rotations do not". That story is
//! **false in both directions**, and this demo exists partly to kill it. The six
//! elements that survive split `det = +1` three ways and `det = −1` three ways —
//! `det_plus_vertex_exact = 3` and `det_minus_vertex_exact = 3` on six of the
//! eight `marching_cubes` rows in `p-57.csv`. Concretely:
//!
//! - `102/+++` swaps `x` and `y`. That is a **reflection** in the plane `x = y`,
//!   `det = −1`, and it is bit-exact. A mirrored asset *can* dedup.
//! - `012/--+` negates `x` and `y`. That is a **180° rotation about `z`**,
//!   `det = +1`, and it **moves**. A rotated asset can fail to dedup.
//!
//! So the predicate is not orientation. It is: *does the element negate any
//! coordinate?* Six of the 48 do not. The tour walks past both counterexamples
//! and the HUD calls each of them out by name when it lands on one.
//!
//! # The fixture is load-bearing, and one wrong number would make this lie
//!
//! `mesh(g·f)` versus `g·mesh(f)` is a statement about the **extractor** only if
//! the sample grid is itself exactly closed under `g`. A grid `origin + i·h`
//! mirrors bit-exactly only when `origin = −((n−1)/2)·h` with `n` odd **and**
//! `h` is a binary fraction, so that `i·h` is exact for every `i`.
//!
//! This demo uses the field's own symmetric domain at **33³**, where the crate's
//! spacing `2L/(n−1)` is `L/16` — dyadic on every reference field. The crate's
//! own 25³ spacing is `2L/24 = L/12`, which has a factor of three in the
//! denominator and **fails the bit-exact mirror test on 16 of 25 coordinates per
//! axis**. Running there would report a falsification that belonged to the grid,
//! not to `marching_cubes`. [`GridFacts::grid_symmetric`] is checked for the
//! selected field before anything is meshed and a failure goes on the HUD as
//! `SELF CHECK FAILED`, because every number below it would be void.
//!
//! `SAMPLES` is deliberately **not** settable from the environment, unlike every
//! other resolution in this repo: the fixture *is* the finding.
//!
//! # Signed zero is a representation, and it is folded on both sides
//!
//! 33 is odd and the grid is centred, so the coordinate `0.0` is on the grid —
//! and `−(0.0)` is `−0.0`, whose bit pattern is not `0.0`'s. Left raw, every one
//! of the 42 sign-flipping elements "fails" on every field for a reason that has
//! nothing to do with any extractor: the two sides agree exactly and disagree
//! only about which encoding of zero got written down. So [`key`] folds `−0.0`
//! onto `0.0`, identically on both sides, and touches no other value — a
//! one-ULP difference anywhere else is still a failure.
//!
//! That fold is not cosmetic and the ledger measures what it buys:
//! `elements_vertex_exact_raw` is **6 on all eight fields**, while the folded
//! `elements_vertex_exact` is 6 on six of them, **24** on `thin_plate` and
//! **48** on `box_exact`. Both numbers are on the startup cross-check table.
//!
//! # The number that closes the mechanism
//!
//! Before any extractor runs, [`grid_facts`] walks every axis-aligned grid edge
//! that changes sign and interpolates it **twice** — forward as
//! `p_lo + (p_hi − p_lo)·(a/(a−b))` and from the far end as
//! `p_hi + (p_lo − p_hi)·(b/(b−a))`. An edge whose two answers differ by a bit
//! is an edge where *the order the endpoints are visited in decides the vertex*.
//! On all sixteen `marching_cubes` rows of `p-57.csv`,
//!
//! ```text
//! worst_differing_vertices == order_sensitive_edges
//! ```
//!
//! **exactly**. A quantity computed from the grid with no extractor in sight
//! predicts, to the unit, how many vertices move. That is the mechanism closing,
//! not a correlation, and it is the pair this demo reproduces:
//!
//! | field | cut edges | order-sensitive | moved | exact of 48 |
//! |---|---:|---:|---:|---:|
//! | `noise_cavity` | 6,522 | 643 | 643 | 6 |
//! | `csg_difference` | 1,386 | 50 | 50 | 6 |
//! | `box_exact` | 1,350 | **0** | **0** | **48** |
//! | `sphere` | 1,158 | 72 | 72 | 6 |
//! | `torus` | 1,128 | 152 | 152 | 6 |
//! | `thin_plate` | 510 | 450 | 450 | **24** |
//! | `gyroid` | 5,292 | 532 | 532 | 6 |
//! | `fbm_terrain` | 2,069 | 291 | 291 | 6 |
//!
//! `cut edges == vertices` on every row because Marching Cubes puts exactly one
//! vertex on each cut grid edge, so the "moved" column is also a share of the
//! whole mesh: 9.9% on `noise_cavity`, 88% on `thin_plate`.
//!
//! # `box_exact` is the control, and it is the one field where mirroring is safe
//!
//! All 1,350 of its cut edges are order-**in**sensitive: its zero set lies on the
//! planes `|x| = 1`, `|y| = 1`, `|z| = 1`, which any dyadic grid hits exactly, so
//! `a/(a−b)` and `b/(b−a)` are both exactly representable and agree. It reaches
//! **48 of 48**, and it is also the field where `fixture_can_fail` is `false` —
//! a pass there is not evidence about the extractor, which is why it is offered
//! as a control and not as the headline. `thin_plate` sits in between at 24: its
//! plate is normal to one axis, so the four sign patterns that leave that axis
//! alone are exact and the four that negate it are not.
//!
//! # What is on screen
//!
//! Three chunks of the same field, side by side, all the same colour, all the
//! same size:
//!
//! 1. **`mesh(f)`** — the asset as authored. The reference every comparison is
//!    against.
//! 2. **`mesh(r·f)`** with `r = 120/+++`, a 120° rotation about `(1,1,1)` that
//!    relabels the axes and negates nothing. Bit-identical, on every field, in
//!    every run.
//! 3. **`mesh(g·f)`** for the element under test, which the tour walks through
//!    all 48 of. Every vertex with no counterpart in `g·mesh(f)` is tinted hot
//!    and carries a marker.
//!
//! Nothing about the picture distinguishes a chunk that dedups from one that
//! does not — that is the point, and it is why the markers exist. The comparison
//! is on the **`f64`** positions the extractor produced; the `f32` cast and the
//! uniform display scale happen after it and are display only.
//!
//! # Which numbers decide it
//!
//! - **`order-sensitive` against `moved`.** Equal is the mechanism. Unequal on
//!   `marching_cubes` would mean this demo has the wrong grid or the wrong key.
//! - **`tally`.** How many of the 48 walked so far were bit-identical, and how
//!   they split by determinant. On six of the eight fields it ends at 6 = 3 + 3.
//! - **`p-57.csv` against `this run`.** Eight numbers per field, checked at
//!   startup and printed as a table. Any disagreement is loud.
//!
//! Every live figure is measured in this process. The only numbers read from
//! `p-57.csv` are the committed ones, quoted for comparison and never mixed in.

mod common;

use bevy::asset::RenderAssetUsages;
use bevy::mesh::{Indices, PrimitiveTopology};
use bevy::platform::time::Instant;
use bevy::prelude::*;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{
    BoxExact, CappedGyroid, CsgDifference, FbmTerrain, NoiseCavity, ReferenceField, Sphere,
    ThinPlate, Torus, capped_gyroid, csg_difference, noise_cavity,
};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

// ─── the group ──────────────────────────────────────────────────────────────

/// The 6 axis permutations. Crossed with 8 sign patterns this is all 48.
///
/// The same order the harness uses, so element indices printed here and bit
/// positions in `p-57.csv`'s `vertex_failing_mask` are the same numbering.
const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

/// The order of the octahedral group.
const GROUP_ORDER: usize = 48;

/// The element the second slot is pinned to: `120/+++`.
///
/// Permutation outer, sign inner, so this is `3 · 8 + 0`. A 120° rotation about
/// the `(1, 1, 1)` diagonal — it relabels all three axes and negates none, which
/// is exactly the property that makes it exact.
const ROTATION: usize = 24;

/// The element the tour opens and closes its intro on: `012/-++`.
///
/// A mirror in `x`, and `p-57.csv`'s `first_failing_element` on six of the eight
/// fields.
const MIRROR: usize = 1;

/// Probe points for the inverse round-trip check.
///
/// Both zeros deliberately: `−(−0.0)` is `0.0`, so a signed permutation
/// round-trips a zero bit-exactly even though it does not preserve it.
const PROBES: [[f64; 3]; 4] = [
    [0.3, -1.7, 2.9],
    [0.0, -0.0, 1.5],
    [-2.25, 0.656_25, -0.187_5],
    [1.0, 1.0, 1.0],
];

/// One element of the octahedral group, as a signed axis permutation.
///
/// `apply(p)[k] = sign[k] · p[perm[k]]`. Signs are `i8` and applied by negation,
/// so the action is bit-exact by construction rather than by an argument about
/// multiplying by `±1.0`.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Element {
    /// Which source axis each output component reads.
    perm: [usize; 3],
    /// `+1` or `−1` per output component.
    sign: [i8; 3],
}

impl Element {
    /// `g·p`. Permute, then negate where the sign says so.
    #[inline]
    fn apply(self, p: [f64; 3]) -> [f64; 3] {
        std::array::from_fn(|k| {
            let v = p[self.perm[k]];
            if self.sign[k] < 0 { -v } else { v }
        })
    }

    /// `g⁻¹`, exactly.
    ///
    /// From `q[k] = sign[k]·p[perm[k]]`: with `j = perm[k]`,
    /// `p[j] = sign[k]·q[k]`, so the inverse reads axis `k` into slot `j` with
    /// the same sign.
    fn inverse(self) -> Self {
        let mut perm = [0usize; 3];
        let mut sign = [0i8; 3];
        for k in 0..3 {
            let j = self.perm[k];
            perm[j] = k;
            sign[j] = self.sign[k];
        }
        Self { perm, sign }
    }

    /// The determinant, `+1` for a rotation and `−1` for a reflection.
    ///
    /// Computed on the integer matrix `m[k][perm[k]] = sign[k]` by the cofactor
    /// formula, so it is exactly `±1` with no float in sight.
    fn det(self) -> i32 {
        let mut m = [[0i32; 3]; 3];
        for k in 0..3 {
            m[k][self.perm[k]] = i32::from(self.sign[k]);
        }
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    /// How many output components are negated. **The predicate that decides
    /// bit-exactness**, which is why it has a name of its own.
    fn negations(self) -> usize {
        self.sign.iter().filter(|s| **s < 0).count()
    }

    /// `perm=012 sign=-++`, for the HUD.
    fn label(self) -> String {
        let mut s = String::from("perm=");
        for k in 0..3 {
            s.push(char::from(b'0' + self.perm[k] as u8));
        }
        s.push_str(" sign=");
        for k in 0..3 {
            s.push(if self.sign[k] < 0 { '-' } else { '+' });
        }
        s
    }

    /// `012/-++` — the compact form `p-57.csv`'s `vertex_failing_labels` uses.
    fn short(self) -> String {
        let mut s = String::new();
        for k in 0..3 {
            s.push(char::from(b'0' + self.perm[k] as u8));
        }
        s.push('/');
        for k in 0..3 {
            s.push(if self.sign[k] < 0 { '-' } else { '+' });
        }
        s
    }

    /// What a level artist would call this reorientation.
    fn what(self) -> &'static str {
        match (self.perm == [0, 1, 2], self.negations()) {
            (true, 0) => "the identity",
            (true, 1) => "a mirror in one axis",
            (true, 2) => "a 180 deg turn about one axis",
            (true, _) => "a point inversion",
            (false, 0) => "a pure axis relabelling",
            (false, 1) => "relabel + one mirror",
            (false, 2) => "relabel + 180 deg turn",
            (false, _) => "relabel + point inversion",
        }
    }
}

/// All 48, permutation outer and sign pattern inner, so element 0 is the
/// identity and the printed index is stable across runs.
///
/// **Not filtered on determinant.** The 24 reflections are exactly as exact as
/// the 24 rotations, and half the point of this demo is that the split does not
/// fall along that line.
fn group() -> Vec<Element> {
    let mut out = Vec::with_capacity(GROUP_ORDER);
    for perm in PERMS {
        for bits in 0..8u8 {
            let sign = std::array::from_fn(|k| if bits & (1 << k) == 0 { 1i8 } else { -1i8 });
            out.push(Element { perm, sign });
        }
    }
    out
}

/// The group is checked before it is used. A wrong group would make every
/// verdict below it meaningless, and this is arithmetic this file owns.
///
/// # Panics
///
/// If the 48 are not distinct, not 24 rotations and 24 reflections, or if any
/// element fails to round-trip a probe point bit-exactly through its own
/// inverse.
fn verify_group(g: &[Element]) {
    assert_eq!(g.len(), GROUP_ORDER, "the octahedral group has 48 elements");
    assert!(
        g[0] == Element {
            perm: [0, 1, 2],
            sign: [1, 1, 1]
        },
        "element 0 must be the identity"
    );
    assert_eq!(g[ROTATION].short(), "120/+++", "ROTATION indexes 120/+++");
    assert_eq!(g[MIRROR].short(), "012/-++", "MIRROR indexes 012/-++");
    for (i, a) in g.iter().enumerate() {
        for b in &g[i + 1..] {
            assert!(a != b, "duplicate group element at {i}: {}", a.label());
        }
    }
    let mut rotations = 0;
    let mut reflections = 0;
    for e in g {
        match e.det() {
            1 => rotations += 1,
            -1 => reflections += 1,
            d => panic!("{} has determinant {d}, not +-1", e.label()),
        }
        let inv = e.inverse();
        for p in PROBES {
            let round = inv.apply(e.apply(p));
            for k in 0..3 {
                assert_eq!(
                    round[k].to_bits(),
                    p[k].to_bits(),
                    "{} does not round-trip {p:?} bit-exactly",
                    e.label()
                );
            }
        }
    }
    assert_eq!(rotations, 24, "24 elements must have det = +1");
    assert_eq!(reflections, 24, "24 elements must have det = -1");
}

// ─── the fixture ────────────────────────────────────────────────────────────

/// Samples per axis.
///
/// **Not settable from the environment, unlike every other resolution in this
/// repo, and that is deliberate.** At 33 over `[−L, L]` the spacing is `L/16`, a
/// binary fraction on every reference field, so the grid is closed under a sign
/// flip *to the bit*. At the crate's own 25³ the spacing is `L/12` and 16 of the
/// 25 coordinates per axis fail that test — a run there would report a property
/// of the grid as a property of `marching_cubes`. See the module header.
const SAMPLES: u32 = 33;

/// Grid spacing for a field of half-extent `L`. `SAMPLES − 1 == 32`, and `L/16`
/// is what `2L/32` reduces to.
fn cell_size(half_extent: f64) -> f64 {
    half_extent / 16.0
}

// ─── the rock ───────────────────────────────────────────────────────────────

/// How many fields the demo offers.
const ROCKS: usize = 8;

/// One of the crate's eight reference fields, chosen at runtime.
///
/// An enum rather than a generic, because the field is picked by a digit key and
/// the whole scene has to change with it. The `Sdf` impl below dispatches to the
/// concrete field, so what gets meshed is bit-for-bit what the harness meshed.
enum Rock {
    /// Perlin noise capped to a sphere — a boulder full of cavities, and the
    /// closest thing in the crate to a carved chunk of terrain.
    NoiseCavity(NoiseCavity<f64>),
    /// The `[−1, 1]³` block with a `0.75` sphere scooped out of its `+++`
    /// corner: a stone block a tool has bitten.
    CsgDifference(CsgDifference<f64>),
    /// The `[−1, 1]³` block. **The control** — its zero set is on planes any
    /// dyadic grid hits exactly, so it is the one field where mirroring is
    /// bit-exact, and the one field whose fixture cannot fail.
    BoxExact(BoxExact<f64>),
    /// A sphere. Fully symmetric as a *field*, and still only 6 of 48 as a mesh.
    Sphere(Sphere<f64>),
    /// A torus, whose axis a relabelling visibly moves.
    Torus(Torus<f64>),
    /// A thin plate. **24 of 48**, because its plate is normal to one axis.
    ThinPlate(ThinPlate<f64>),
    /// A capped gyroid — the densest surface here, 5,292 vertices at 33³.
    Gyroid(CappedGyroid<f64>),
    /// A fractal heightfield. Not closed in its own domain, which is why it is
    /// last and needs `ISOMESH_FIELD=7`.
    FbmTerrain(FbmTerrain<f64>),
}

/// Run `$body` against whichever concrete field `$rock` holds.
///
/// The body is expanded once per variant, so each instance type-checks against
/// its own field type. Three uses — `sample`, `gradient` and `domain` — and
/// writing them out would be 24 near-identical arms.
macro_rules! on_rock {
    ($rock:expr, |$f:ident| $body:expr) => {
        match $rock {
            Rock::NoiseCavity($f) => $body,
            Rock::CsgDifference($f) => $body,
            Rock::BoxExact($f) => $body,
            Rock::Sphere($f) => $body,
            Rock::Torus($f) => $body,
            Rock::ThinPlate($f) => $body,
            Rock::Gyroid($f) => $body,
            Rock::FbmTerrain($f) => $body,
        }
    };
}

impl Rock {
    /// The field at `index`, ordered for this demo rather than for the crate.
    ///
    /// Index 0 is the boulder, because that is the asset a level artist mirrors;
    /// 1 is the carved block; 2 is the control. Digits `1`-`7` in the shared
    /// harness map to indices 0-6, so index 7 needs `ISOMESH_FIELD=7`.
    fn at(index: usize) -> Self {
        match index % ROCKS {
            0 => Self::NoiseCavity(noise_cavity()),
            1 => Self::CsgDifference(csg_difference()),
            2 => Self::BoxExact(BoxExact::canonical()),
            3 => Self::Sphere(Sphere::canonical()),
            4 => Self::Torus(Torus::canonical()),
            5 => Self::ThinPlate(ThinPlate::canonical()),
            6 => Self::Gyroid(capped_gyroid()),
            _ => Self::FbmTerrain(FbmTerrain::canonical()),
        }
    }

    /// The `ReferenceField::NAME` the CSV keys on.
    fn name(&self) -> &'static str {
        on_rock!(self, |f| field_name(f))
    }

    /// Half the side of the field's own cubic domain, `L`.
    ///
    /// Read from the field rather than tabulated: `ReferenceField::domain` is
    /// what the harness read, and a second table of the same eight numbers is a
    /// second thing to keep in step.
    fn half_extent(&self) -> f64 {
        on_rock!(self, |f| ReferenceField::domain(f).1[0])
    }

    /// The domain is the symmetric cube both the harness and this demo assume.
    fn domain_is_symmetric(&self) -> bool {
        on_rock!(self, |f| {
            let (lo, hi) = ReferenceField::domain(f);
            (0..3).all(|k| lo[k].to_bits() == (-hi[k]).to_bits())
        })
    }
}

/// `F::NAME`, as a free function so [`on_rock`] can reach an associated const.
fn field_name<F: ReferenceField>(_field: &F) -> &'static str {
    F::NAME
}

impl Sdf for Rock {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        on_rock!(self, |f| f.sample(p))
    }

    #[inline]
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        on_rock!(self, |f| f.gradient(p))
    }
}

/// `g·f`, the field pushed forward by `g`: `(g·f)(p) = f(g⁻¹·p)`.
///
/// Both `g` and `g⁻¹` are stored because both are needed and neither costs
/// anything: the sample point goes in through `g⁻¹` and the gradient comes back
/// out through `g`.
struct Rotated<'a> {
    /// The field being reoriented.
    field: &'a Rock,
    /// The element.
    g: Element,
    /// Its inverse, precomputed.
    g_inv: Element,
}

impl Sdf for Rotated<'_> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        self.field.sample(self.g_inv.apply(p))
    }

    /// `∇(f∘g⁻¹)(p) = (g⁻¹)ᵀ·(∇f)(g⁻¹·p)`, and `(g⁻¹)ᵀ = g` because `g` is
    /// orthogonal.
    ///
    /// Overridden rather than inherited. `Sdf`'s default is a six-sample central
    /// difference, and letting that in here would put the differencing stencil
    /// into the normals rather than the field's own gradient.
    #[inline]
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        self.g.apply(self.field.gradient(self.g_inv.apply(p)))
    }
}

// ─── bit keys ───────────────────────────────────────────────────────────────

/// The bit pattern of `−0.0`.
const NEGATIVE_ZERO: u64 = 1u64 << 63;

/// A position component as a comparison key, with `−0.0` folded onto `0.0`.
///
/// Nothing else is touched — see the module header. This is not a tolerance: a
/// one-ULP difference on any other value is still a difference.
#[inline]
fn key(v: f64) -> u64 {
    let b = v.to_bits();
    if b == NEGATIVE_ZERO { 0 } else { b }
}

/// The sign-magnitude-ordered integer image of an `f64`'s bits.
///
/// Monotone in the value and continuous across zero, so a difference of these is
/// a ULP count that means something for a pair straddling zero.
#[inline]
fn monotone(bits: u64) -> i128 {
    if bits & NEGATIVE_ZERO == 0 {
        i128::from(bits)
    } else {
        -i128::from(bits & !NEGATIVE_ZERO)
    }
}

/// ULP distance between two `f64` bit patterns.
#[inline]
fn ulp_distance(a: u64, b: u64) -> u128 {
    (monotone(a) - monotone(b)).unsigned_abs()
}

/// Sorted multiset of vertex positions as bit triples, optionally mapped through
/// a group element first.
fn vertex_keys(positions: &[[f64; 3]], g: Option<Element>) -> Vec<[u64; 3]> {
    let mut out: Vec<[u64; 3]> = positions
        .iter()
        .map(|p| {
            let q = match g {
                Some(e) => e.apply(*p),
                None => *p,
            };
            std::array::from_fn(|k| key(q[k]))
        })
        .collect();
    out.sort_unstable();
    out
}

/// What is left after cancelling every vertex the two meshes agree on.
struct Residue {
    /// Keys present in `mesh(g·f)` and absent from `g·mesh(f)`. Sorted.
    only_got: Vec<[u64; 3]>,
    /// How many keys went the other way. Equal to `only_got.len()` whenever the
    /// vertex counts match, which on `marching_cubes` they always do.
    only_want: usize,
    /// Largest ULP gap over the paired residues.
    worst_ulp: u128,
}

impl Residue {
    /// Bit-identical means *nothing* left over on either side.
    fn exact(&self) -> bool {
        self.only_got.is_empty() && self.only_want == 0
    }
}

/// The multiset symmetric difference of `mesh(g·f)` against `g·mesh(f)`.
///
/// A **sorted merge**, not a positional walk over the two sorted lists, and that
/// distinction is not cosmetic: diffing entry `i` against entry `i` reports
/// nonsense the moment the lists differ by an *insertion* rather than a
/// perturbation, because every later entry is then compared against its
/// neighbour. The harness's first run reported a `9.2e18` ULP gap — the distance
/// from `+2` to `−2` — for meshes that differ by one bit on 72 edges, for exactly
/// that reason. The merge cancels every agreeing vertex and leaves the ones that
/// moved, and pairing the two residues in sorted order pairs each moved vertex
/// with the position it moved from.
fn residue(rotated: &[[f64; 3]], reference: &[[f64; 3]], g: Element) -> Residue {
    let got = vertex_keys(rotated, None);
    let want = vertex_keys(reference, Some(g));

    let mut only_got: Vec<[u64; 3]> = Vec::new();
    let mut only_want: Vec<[u64; 3]> = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < got.len() && j < want.len() {
        match got[i].cmp(&want[j]) {
            std::cmp::Ordering::Less => {
                only_got.push(got[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                only_want.push(want[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    only_got.extend_from_slice(&got[i..]);
    only_want.extend_from_slice(&want[j..]);

    let mut worst_ulp = 0u128;
    for (a, b) in only_got.iter().zip(only_want.iter()) {
        for k in 0..3 {
            if a[k] != b[k] {
                worst_ulp = worst_ulp.max(ulp_distance(a[k], b[k]));
            }
        }
    }
    Residue {
        only_got,
        only_want: only_want.len(),
        worst_ulp,
    }
}

/// Which vertices of the displayed mesh are the ones that moved.
///
/// `only_got` is sorted by construction, so this is a binary search per vertex
/// rather than a hash set — and the residue is the set the count came from, so
/// the markers cannot disagree with the number on the HUD.
fn moved_flags(positions: &[[f64; 3]], residue: &Residue) -> Vec<bool> {
    positions
        .iter()
        .map(|p| {
            let k: [u64; 3] = std::array::from_fn(|c| key(p[c]));
            residue.only_got.binary_search(&k).is_ok()
        })
        .collect()
}

// ─── the grid, before any extractor runs ────────────────────────────────────

/// What is true of a `(field, 33³)` pair with no extractor involved.
struct GridFacts {
    /// `pos[i] == −pos[n−1−i]` bit-exactly, on every `i` and every axis. **If
    /// this is false the whole demo is void**, which is why it is on the HUD.
    grid_symmetric: bool,
    /// Grid coordinates that are exactly zero. One, for a centred odd grid, and
    /// the reason [`key`] folds the sign of zero at all.
    zero_coordinates: usize,
    /// Axis-aligned grid edges that change sign. On `marching_cubes` this is
    /// also the vertex count.
    cut_edges: usize,
    /// Of those, the ones whose crossing coordinate depends on which end it was
    /// interpolated from. **The number that predicts how many vertices move.**
    order_sensitive_edges: usize,
}

/// The axis coordinate list, built exactly as the extractor builds it —
/// `origin + cell_size · i`, per `marching_cubes/mod.rs:229`.
fn axis_coords(origin: f64, cell: f64) -> Vec<f64> {
    (0..SAMPLES).map(|i| origin + cell * f64::from(i)).collect()
}

/// Walk the grid and every cut edge on it, twice.
///
/// One pass over `33³` samples and one over the `3 · 33² · 32` axis-aligned
/// edges. Nothing here calls an extractor, which is the whole point: the
/// `order_sensitive_edges` it returns is what `marching_cubes` is about to be
/// held against.
fn grid_facts(field: &Rock, origin: f64, cell: f64) -> GridFacts {
    let n = SAMPLES as usize;
    let coords = axis_coords(origin, cell);

    let mut grid_symmetric = true;
    let mut zero_coordinates = 0;
    for i in 0..n {
        if key(coords[i]) != key(-coords[n - 1 - i]) {
            grid_symmetric = false;
        }
        if key(coords[i]) == 0 {
            zero_coordinates += 1;
        }
    }

    let mut values = Vec::with_capacity(n * n * n);
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                values.push(field.sample([coords[x], coords[y], coords[z]]));
            }
        }
    }
    let at = |x: usize, y: usize, z: usize| values[x + n * (y + n * z)];

    let mut cut_edges = 0usize;
    let mut order_sensitive_edges = 0usize;
    for axis in 0..3 {
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let mut hi = [x, y, z];
                    hi[axis] += 1;
                    if hi[axis] >= n {
                        continue;
                    }
                    let va = at(x, y, z);
                    let vb = at(hi[0], hi[1], hi[2]);
                    if (va < 0.0) == (vb < 0.0) {
                        continue;
                    }
                    cut_edges += 1;
                    let lo_c = coords[[x, y, z][axis]];
                    let hi_c = coords[hi[axis]];
                    // The two interpolations `EDGE_CORNERS` chooses between.
                    let fwd = lo_c + (hi_c - lo_c) * (va / (va - vb));
                    let rev = hi_c + (lo_c - hi_c) * (vb / (vb - va));
                    if key(fwd) != key(rev) {
                        order_sensitive_edges += 1;
                    }
                }
            }
        }
    }

    GridFacts {
        grid_symmetric,
        zero_coordinates,
        cut_edges,
        order_sensitive_edges,
    }
}

// ─── the sweep ──────────────────────────────────────────────────────────────

/// Everything one field's 48-element sweep measured, plus the reference mesh it
/// measured against.
struct Sweep {
    /// The grid census, taken before the first extraction.
    facts: GridFacts,
    /// Vertices in `mesh(f)`.
    vertices: usize,
    /// Triangles in `mesh(f)`.
    triangles: usize,
    /// Bit `i` set = element `i` produced a bit-identical vertex multiset.
    exact_mask: u64,
    /// Largest number of vertices that moved, over all 48.
    worst_differing_vertices: usize,
    /// Largest ULP gap anywhere in the sweep.
    worst_component_ulp: u128,
    /// Of the 6 elements that negate nothing, how many were exact.
    pure_permutation_exact: usize,
    /// Of the 8 with `perm = 012`, how many were exact.
    pure_sign_flip_exact: usize,
    /// Rotations (`det = +1`) that were exact, out of 24.
    det_plus_exact: usize,
    /// Reflections (`det = −1`) that were exact, out of 24.
    det_minus_exact: usize,
    /// Lowest-indexed failure, for the CSV's `first_failing_element`.
    first_failing: Option<usize>,
    /// The identity element reproduced the reference exactly. **The negative
    /// control**: it goes through [`Rotated`] like every other element, so it
    /// checks the extractor is deterministic and the wrapper is transparent at
    /// once.
    identity_exact: bool,
    /// `mesh(f)`'s vertex positions, in `f64`, kept because every comparison is
    /// against them and re-extracting would be a second path to the same mesh.
    reference_positions: Vec<[f64; 3]>,
    /// `mesh(f)`'s normals, for the display.
    reference_normals: Vec<[f64; 3]>,
    /// `mesh(f)`'s index buffer, for the display.
    reference_indices: Vec<u32>,
    /// World units per display unit: `DISPLAY_HALF / max|position|`, so every
    /// field renders at the same size. **Display only** — the comparison never
    /// sees it.
    display_scale: f32,
    /// Wall time for the 49 extractions. Gates nothing; M-348 is the incident
    /// where a discovery was demoted for resting on a wall clock.
    wall_ms: f64,
}

impl Sweep {
    /// How many of the first `seen` elements of the group were exact.
    fn exact_among(&self, seen: usize) -> usize {
        let mask = if seen >= GROUP_ORDER {
            u64::MAX
        } else {
            (1u64 << seen) - 1
        };
        (self.exact_mask & mask).count_ones() as usize
    }

    /// Total exact, out of 48.
    fn exact_total(&self) -> usize {
        self.exact_mask.count_ones() as usize
    }
}

/// Mesh a field under all 48 elements and compare each against the reference.
///
/// One extraction of `mesh(f)` and 48 of `mesh(g·f)`, on the 33³ fixture. This
/// is the same arithmetic the per-frame rebuild uses — [`residue`] — so a number
/// in the startup table and a number on the HUD cannot be produced by two
/// different pieces of code.
fn sweep(field: &Rock, elements: &[Element], mc: &mut MarchingCubes<f64>) -> Sweep {
    let l = field.half_extent();
    let cell = cell_size(l);
    let origin = [-l; 3];
    let shape = RuntimeShape3::new([SAMPLES; 3]).expect("a 33 cube fits u32");
    let facts = grid_facts(field, -l, cell);

    let started = Instant::now();
    let mut reference = MeshBuffer::<f64>::new();
    let mut rotated = MeshBuffer::<f64>::new();
    if let Err(e) = mc.extract(field, &shape, origin, cell, &mut reference) {
        error!("{}: reference extraction failed: {e}", field.name());
    }

    let mut exact_mask = 0u64;
    let mut worst_differing_vertices = 0usize;
    let mut worst_component_ulp = 0u128;
    let mut pure_permutation_exact = 0usize;
    let mut pure_sign_flip_exact = 0usize;
    let mut det_plus_exact = 0usize;
    let mut det_minus_exact = 0usize;
    let mut first_failing = None;
    let mut identity_exact = false;

    for (index, &g) in elements.iter().enumerate() {
        rotated.reset();
        let wrapped = Rotated {
            field,
            g,
            g_inv: g.inverse(),
        };
        if let Err(e) = mc.extract(&wrapped, &shape, origin, cell, &mut rotated) {
            error!("{}: {} extraction failed: {e}", field.name(), g.short());
            continue;
        }
        let res = residue(&rotated.positions, &reference.positions, g);
        worst_differing_vertices = worst_differing_vertices.max(res.only_got.len());
        worst_component_ulp = worst_component_ulp.max(res.worst_ulp);
        if res.exact() {
            exact_mask |= 1u64 << index;
            if g.negations() == 0 {
                pure_permutation_exact += 1;
            }
            if g.perm == [0, 1, 2] {
                pure_sign_flip_exact += 1;
            }
            if g.det() > 0 {
                det_plus_exact += 1;
            } else {
                det_minus_exact += 1;
            }
            if index == 0 {
                identity_exact = true;
            }
        } else if first_failing.is_none() {
            first_failing = Some(index);
        }
    }

    let scale = reference
        .positions
        .iter()
        .flat_map(|p| p.iter().map(|v| v.abs()))
        .fold(0.0f64, f64::max);
    Sweep {
        facts,
        vertices: reference.positions.len(),
        triangles: reference.triangle_count(),
        exact_mask,
        worst_differing_vertices,
        worst_component_ulp,
        pure_permutation_exact,
        pure_sign_flip_exact,
        det_plus_exact,
        det_minus_exact,
        first_failing,
        identity_exact,
        display_scale: if scale > 0.0 {
            DISPLAY_HALF / scale as f32
        } else {
            1.0
        },
        reference_positions: reference.positions.clone(),
        reference_normals: reference.normals.clone(),
        reference_indices: reference.indices.clone(),
        wall_ms: started.elapsed().as_secs_f64() * 1000.0,
    }
}

// ─── the ledger, compiled in ────────────────────────────────────────────────

/// P-57's committed artefact, embedded at compile time.
///
/// `include_str!` rather than transcribed constants: the path resolves against
/// this source file so no working directory can break it, and a number that
/// lived only here could drift away from the CSV with nothing to say so.
const LEDGER_CSV: &str = include_str!("../../docs/experiments/p-57.csv");

/// The one row family this demo can honestly be held against.
const LEDGER_EXTRACTOR: &str = "marching_cubes";

/// One committed `marching_cubes` row of `p-57.csv` at 33³.
#[derive(Clone, Copy)]
struct LedgerRow {
    /// `vertices`.
    vertices: usize,
    /// `triangles`.
    triangles: usize,
    /// `cut_edges`.
    cut_edges: usize,
    /// `order_sensitive_edges` — the headline cross-check.
    order_sensitive_edges: usize,
    /// `worst_differing_vertices` — the other half of the pair.
    worst_differing_vertices: usize,
    /// `elements_vertex_exact`, out of 48.
    elements_vertex_exact: usize,
    /// `elements_vertex_exact_raw`, i.e. without the signed-zero fold.
    elements_vertex_exact_raw: usize,
    /// `pure_permutation_exact`, out of 6.
    pure_permutation_exact: usize,
    /// `worst_component_ulp`.
    worst_component_ulp: u128,
    /// `grid_zero_coordinates`.
    grid_zero_coordinates: usize,
    /// `fixture_can_fail`.
    fixture_can_fail: bool,
}

impl LedgerRow {
    /// Pull the `marching_cubes`, 33³ row for `field` out of the CSV by header
    /// name, so a reordered column cannot silently shift a value.
    fn load(field: &str) -> Option<Self> {
        let mut lines = LEDGER_CSV.lines().filter(|l| !l.starts_with('#'));
        let header: Vec<&str> = lines.next()?.split(',').collect();
        let column = |name: &str| header.iter().position(|h| *h == name);
        let c_field = column("field")?;
        let c_extractor = column("extractor")?;
        let c_samples = column("samples_per_axis")?;
        let want = [
            column("vertices")?,
            column("triangles")?,
            column("cut_edges")?,
            column("order_sensitive_edges")?,
            column("worst_differing_vertices")?,
            column("elements_vertex_exact")?,
            column("elements_vertex_exact_raw")?,
            column("pure_permutation_exact")?,
            column("worst_component_ulp")?,
            column("grid_zero_coordinates")?,
        ];
        let c_can_fail = column("fixture_can_fail")?;
        for line in lines {
            let cells: Vec<&str> = line.split(',').collect();
            if cells.len() != header.len()
                || cells[c_field] != field
                || cells[c_extractor] != LEDGER_EXTRACTOR
                || cells[c_samples] != "33"
            {
                continue;
            }
            let n: Vec<u128> = want
                .iter()
                .map(|c| cells[*c].parse::<u128>().unwrap_or(u128::MAX))
                .collect();
            if n.contains(&u128::MAX) {
                return None;
            }
            return Some(Self {
                vertices: n[0] as usize,
                triangles: n[1] as usize,
                cut_edges: n[2] as usize,
                order_sensitive_edges: n[3] as usize,
                worst_differing_vertices: n[4] as usize,
                elements_vertex_exact: n[5] as usize,
                elements_vertex_exact_raw: n[6] as usize,
                pure_permutation_exact: n[7] as usize,
                worst_component_ulp: n[8],
                grid_zero_coordinates: n[9] as usize,
                fixture_can_fail: cells[c_can_fail] == "true",
            });
        }
        None
    }

    /// Every number this run can be held against, and whether it matches.
    ///
    /// A `Vec` of `(name, expected, measured)` rather than one boolean, because
    /// "agrees" is worth nothing without saying what agreed.
    fn against(&self, sweep: &Sweep) -> Vec<(&'static str, u128, u128)> {
        vec![
            (
                "order_sensitive_edges",
                self.order_sensitive_edges as u128,
                sweep.facts.order_sensitive_edges as u128,
            ),
            (
                "worst_differing_vertices",
                self.worst_differing_vertices as u128,
                sweep.worst_differing_vertices as u128,
            ),
            ("vertices", self.vertices as u128, sweep.vertices as u128),
            ("triangles", self.triangles as u128, sweep.triangles as u128),
            (
                "cut_edges",
                self.cut_edges as u128,
                sweep.facts.cut_edges as u128,
            ),
            (
                "elements_vertex_exact",
                self.elements_vertex_exact as u128,
                sweep.exact_total() as u128,
            ),
            (
                "pure_permutation_exact",
                self.pure_permutation_exact as u128,
                sweep.pure_permutation_exact as u128,
            ),
            (
                "worst_component_ulp",
                self.worst_component_ulp,
                sweep.worst_component_ulp,
            ),
            (
                "grid_zero_coordinates",
                self.grid_zero_coordinates as u128,
                sweep.facts.zero_coordinates as u128,
            ),
        ]
    }
}

// ─── layout, framing and colour ─────────────────────────────────────────────

/// Half-extent every field is scaled to for display, in world units.
///
/// The three slots are laid out in these units, so the framing does not change
/// when the field does — `gyroid`'s domain is 3.5x `sphere`'s and a fixed camera
/// would show one of them as a speck.
const DISPLAY_HALF: f32 = 1.0;

/// Distance between slot centres, in world units.
///
/// Constrained from below by the *labels*, not by the chunks: a slot label has
/// to say `MOVED -- 226 vertices differ`, which is 267 px at the 15 px font, so
/// the slot pitch cannot be under about 280 px. `3.0 · 92 = 276` px, and the
/// 184 px chunks then sit in it with a 92 px gutter.
const SLOT_GAP: f32 = 3.0;

/// Logical pixels per world unit on the 1280x720 capture, at
/// [`CAMERA_RADIUS`] and Bevy's default 45° vertical FOV.
///
/// Not used for rendering — it is what the slot labels are positioned by, and
/// what makes [`SLOT_GAP`], [`CAMERA_RADIUS`] and [`FOCUS_LIFT`] checkable
/// rather than tuned blind. `720 / 92 = 7.826` units of height, so
/// `radius = 3.913 / tan(22.5°) = 9.447`.
const PIXELS_PER_UNIT: f32 = 92.0;

/// Orbit radius, in world units. See [`PIXELS_PER_UNIT`] for where it comes
/// from.
const CAMERA_RADIUS: f32 = 9.447;

/// Radians of pitch. Small: three chunks compared side by side want to be seen
/// from the same angle, and a steep look-down turns the outer two into
/// foreshortened ellipses.
const CAMERA_PITCH: f32 = 0.20;

/// How far the camera looks *above* the row, in world units.
///
/// `(496/720 − 0.5) · 7.826`, which drops the 184 px row of chunks to
/// **404-588 px** down a 720 px frame: clear of the HUD panel above it and of
/// the labels below it. Every one of those three numbers has to hold at once,
/// which is why they are derived rather than nudged.
const FOCUS_LIFT: f32 = 1.478;

/// Width and height of the backdrop the HUD is read against, in logical pixels.
///
/// Measured on a 1280x720 capture rather than chosen. The HUD reaches 25 lines
/// at the harness's 13 px font, whose pitch is 15.6 px, so
/// `10 + 25 · 15.6 = 400`; and its longest line is 78 characters at 8.12 px
/// each, so `12 + 78 · 8.12 = 645`. It has to stop clear of the chunks at
/// 404 px — which is why [`report`] keeps its block to eleven lines and none of
/// them wide, and why the punchline about the two counterexamples lives in the
/// caption instead.
const HUD_PANEL: Vec2 = Vec2::new(672.0, 400.0);

/// Where the slot labels sit, in logical pixels from the top.
///
/// The chunks end at 588 px and the two-line caption box starts at 652, so a
/// two-line label has 64 px to live in and this puts it in the middle of them.
const LABEL_TOP: f32 = 594.0;

/// Width of one slot label, in logical pixels.
const LABEL_WIDTH: f32 = 250.0;

/// Gap between slot labels, chosen so their centres land on the slot centres:
/// `LABEL_WIDTH + gap == SLOT_GAP · PIXELS_PER_UNIT`.
const LABEL_GAP: f32 = SLOT_GAP * PIXELS_PER_UNIT - LABEL_WIDTH;

/// Unpainted rock. **The same colour on all three slots** — the whole claim is
/// that the picture does not distinguish them.
const ROCK_SRGB: [f32; 4] = [0.47, 0.44, 0.40, 1.0];

/// A vertex with no counterpart in `g·mesh(f)`.
const MOVED_SRGB: [f32; 4] = [1.0, 0.22, 0.08, 1.0];

/// Half-length of the three-line cross drawn on a moved vertex, in world units.
///
/// A cross rather than `Gizmos::sphere`: 643 spheres is 61,000 line segments at
/// the default resolution, and at this size a sphere and a cross photograph the
/// same. `0.030` is 3 px of arm at [`PIXELS_PER_UNIT`], and the clip is scaled to
/// `WIDTH=1000` afterwards, so it survives at about 2 px — the smallest a mark
/// can be and still read as deliberate.
const MARKER_ARM: f32 = 0.030;

/// Seconds for one pass through the tour, when nobody is capturing.
const TOUR_SECONDS: f32 = 30.0;

/// Where the opening hold on the mirror ends, as a fraction of the tour.
const INTRO_END: f32 = 0.10;

/// Where the walk through all 48 ends and the closing hold begins.
///
/// `0.90 − 0.10 = 0.80` of the clip for 48 elements, which at the recording's
/// 120 captured frames is **exactly two frames each**.
const WALK_END: f32 = 0.90;

/// sRGB as a human picks it into the linear RGBA [`Mesh::ATTRIBUTE_COLOR`]
/// wants. Feeding sRGB in raw renders it washed out (E-208).
fn linear(srgb: [f32; 4]) -> [f32; 4] {
    Color::srgba(srgb[0], srgb[1], srgb[2], srgb[3])
        .to_linear()
        .to_f32_array()
}

/// `n` with a thousands separator, for a HUD a designer reads.
fn commas(n: usize) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

// ─── resources ──────────────────────────────────────────────────────────────

/// The 48 elements, verified once at startup.
#[derive(Resource)]
struct Group(Vec<Element>);

/// One [`Sweep`] per field, measured before the window opens.
#[derive(Resource)]
struct Sweeps(Vec<Sweep>);

/// One [`LedgerRow`] per field, or `None` if the CSV does not carry it.
#[derive(Resource)]
struct Ledger(Vec<Option<LedgerRow>>);

/// What this frame is showing.
#[derive(Resource, Clone, Copy, PartialEq, Eq)]
struct Shot {
    /// Index into [`Rock::at`].
    field: usize,
    /// Index into the group.
    element: usize,
}

impl Default for Shot {
    fn default() -> Self {
        Self {
            field: 0,
            element: MIRROR,
        }
    }
}

/// The steering a human has done.
#[derive(Resource, Default)]
struct Steer {
    /// An element pinned by hand, which stops the tour until `A`.
    pinned: Option<usize>,
    /// Whether the per-vertex crosses are drawn.
    markers: bool,
}

/// The extractor, its buffer, and the three assets the chunks are drawn from.
#[derive(Resource)]
struct Rig {
    /// Kept across frames so a step does not allocate a fresh 6,500-vertex
    /// buffer.
    mc: MarchingCubes<f64>,
    /// Scratch for the slot being rebuilt.
    buffer: MeshBuffer<f64>,
    /// One asset per slot, **replaced in place** rather than swapped.
    ///
    /// Adding a fresh asset and dropping the old one is the obvious shape and it
    /// makes `bevy_render` log `Use-after-free: attempted to copy element data
    /// for an unallocated key` once per slot: the render world had already
    /// queued the handle this frame. `Assets::insert` on the same id fires
    /// `AssetEvent::Modified`, re-uploads, and never frees anything a queued
    /// draw still names.
    meshes: [Handle<Mesh>; 3],
    /// Linear rock colour.
    rock: [f32; 4],
    /// Linear moved-vertex colour.
    moved: [f32; 4],
}

/// What the current `(field, element)` measured.
#[derive(Resource, Default)]
struct Live {
    /// Vertices in `mesh(g·f)`.
    vertices: usize,
    /// Triangles in `mesh(g·f)`.
    triangles: usize,
    /// Vertices of `mesh(g·f)` with no counterpart in `g·mesh(f)`.
    moved: usize,
    /// Largest ULP gap over the paired residues.
    worst_ulp: u128,
    /// Time the last extraction took.
    extract_ms: f64,
    /// Whether the vertex multisets were bit-identical.
    exact: bool,
    /// World-space marker positions for the moved vertices of slot 2.
    markers: Vec<Vec3>,
    /// The `Shot` the assets were built from, so a frame that changes nothing
    /// costs nothing.
    built: Option<Shot>,
}

/// The startup self-check, so a wrong fixture is visible on screen rather than
/// only in the log.
#[derive(Resource)]
struct SelfCheck {
    /// Fields whose 33³ grid is not bit-exactly closed under a sign flip. Must
    /// be zero, or every verdict is about the grid rather than the extractor.
    asymmetric_grids: usize,
    /// Fields where the identity element did not reproduce the reference. Must
    /// be zero, or `marching_cubes` is not deterministic here.
    nondeterministic: usize,
    /// `(field, column, expected, measured)` for every cross-check that
    /// disagreed with `p-57.csv`.
    disagreements: Vec<(&'static str, &'static str, u128, u128)>,
    /// Fields with no `marching_cubes` 33³ row in the CSV at all.
    missing_rows: usize,
}

/// A slot label, by index.
#[derive(Component)]
struct SlotLabel(usize);

/// The bottom caption — the line a viewer reads instead of the HUD.
#[derive(Component)]
struct Caption;

// ─── app ────────────────────────────────────────────────────────────────────

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-314 mirrored is not the same mesh".into(),
                // Web only, inert on native: bind to the 1280x720 canvas the
                // page supplies rather than letting Bevy append its own. The HUD
                // panels are laid out in pixels for that size, so the canvas is
                // fixed and CSS scales it -- `fit_canvas_to_parent` stays at its
                // `false` default for the same reason.
                canvas: Some("#isomesh-canvas".into()),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        // Not black: three grey chunks against a void read as floating polygons.
        .insert_resource(ClearColor(Color::srgb(0.10, 0.12, 0.16)))
        .init_resource::<Shot>()
        .init_resource::<Live>()
        .add_systems(Startup, setup)
        // `PreUpdate`, not `Update`, and it is a correctness fix rather than a
        // preference: the harness's `update_hud` lives in `Update` and system
        // order within a schedule is unspecified, so a caption written there
        // would sometimes disagree with the HUD above it by one frame (E-312).
        .add_systems(PreUpdate, (steer, advance, rebuild, report).chain())
        .add_systems(Update, draw_markers)
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut camera: Query<&mut OrbitCamera>,
    flags: Res<ViewFlags>,
) {
    for mut orbit in &mut camera {
        orbit.focus = Vec3::new(0.0, FOCUS_LIFT, 0.0);
        orbit.radius = CAMERA_RADIUS;
        // Straight down `−z`, so all three chunks are seen from the same angle.
        // A comparison shot whose three subjects are lit and turned differently
        // is not a comparison.
        orbit.yaw = std::f32::consts::FRAC_PI_2;
        orbit.pitch = CAMERA_PITCH;
    }

    // White base colour, because every chunk's colour arrives per vertex: the
    // same attribute carries the moved-vertex tint, so there is one path onto a
    // vertex rather than a base colour and an override.
    let material = materials.add(StandardMaterial {
        base_color: Color::WHITE,
        perceptual_roughness: 0.92,
        metallic: 0.0,
        ..default()
    });
    // One asset per slot, created empty and then written through for the rest of
    // the run. See `Rig::meshes`.
    let slot_meshes: [Handle<Mesh>; 3] = std::array::from_fn(|_| {
        meshes.add(Mesh::new(
            PrimitiveTopology::TriangleList,
            RenderAssetUsages::default(),
        ))
    });
    for (slot, handle) in slot_meshes.iter().enumerate() {
        commands.spawn((
            Mesh3d(handle.clone()),
            MeshMaterial3d(material.clone()),
            DemoMesh,
            Transform::from_xyz(slot_x(slot), 0.0, 0.0),
        ));
    }

    // Behind the harness HUD, which `CommonPlugin` spawns at the default z.
    commands.spawn((
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(6.0),
            left: Val::Px(6.0),
            width: Val::Px(HUD_PANEL.x),
            height: Val::Px(HUD_PANEL.y),
            ..default()
        },
        BackgroundColor(Color::srgba(0.04, 0.05, 0.07, 0.60)),
        GlobalZIndex(-1),
    ));

    // One label per slot, in a centred row whose child pitch is exactly the slot
    // pitch, so a label sits under its own chunk on the 1280x720 capture.
    commands
        .spawn((
            Node {
                position_type: PositionType::Absolute,
                left: Val::Px(0.0),
                right: Val::Px(0.0),
                top: Val::Px(LABEL_TOP),
                justify_content: JustifyContent::Center,
                column_gap: Val::Px(LABEL_GAP),
                ..default()
            },
            GlobalZIndex(4),
        ))
        .with_children(|parent| {
            for slot in 0..3usize {
                parent.spawn((
                    Text::new(""),
                    TextFont {
                        font_size: FontSize::Px(15.0),
                        ..default()
                    },
                    TextColor(Color::srgb(0.96, 0.94, 0.90)),
                    TextLayout {
                        justify: Justify::Center,
                        ..default()
                    },
                    BackgroundColor(Color::srgba(0.03, 0.03, 0.05, 0.80)),
                    Node {
                        width: Val::Px(LABEL_WIDTH),
                        padding: UiRect::axes(Val::Px(6.0), Val::Px(5.0)),
                        ..default()
                    },
                    SlotLabel(slot),
                ));
            }
        });

    commands
        .spawn((
            Node {
                position_type: PositionType::Absolute,
                left: Val::Px(0.0),
                right: Val::Px(0.0),
                bottom: Val::Px(6.0),
                justify_content: JustifyContent::Center,
                ..default()
            },
            GlobalZIndex(4),
        ))
        .with_children(|parent| {
            parent.spawn((
                Text::new(""),
                TextFont {
                    font_size: FontSize::Px(19.0),
                    ..default()
                },
                // `NoWrap`: in a centring flex row the measure is handed the
                // container's whole width but the node's height resolves before
                // the wrap, so a soft wrap pushes the second line off frame.
                TextLayout {
                    linebreak: bevy::text::LineBreak::NoWrap,
                    justify: Justify::Center,
                },
                TextColor(Color::srgb(0.97, 0.94, 0.90)),
                BackgroundColor(Color::srgba(0.03, 0.03, 0.05, 0.84)),
                Node {
                    padding: UiRect::axes(Val::Px(14.0), Val::Px(7.0)),
                    ..default()
                },
                Caption,
            ));
        });

    let elements = group();
    verify_group(&elements);
    let mut mc = MarchingCubes::<f64>::new();

    // The whole table, before the window opens. 8 fields x 49 extractions of a
    // 33 cube, which the harness's own `wall_ms` column puts at 1.4 s in total.
    // Measuring it up front is what lets a pinned still carry the argument: a
    // sweep that only visited the elements the clip walked would report a
    // different number for `ISOMESH_FIELD=2` than for the same field reached by
    // pressing `3`.
    info!(
        "p-57 cross-check, marching_cubes at {SAMPLES}^3 -- expected (p-57.csv) / measured (this run)"
    );
    info!(
        "{:<15} {:>11} {:>11} {:>9} {:>9} {:>7} {:>9} {:>6}",
        "field", "cut_edges", "osEdges", "moved", "vex/48", "p6/6", "ulp", "ms"
    );
    let mut check = SelfCheck {
        asymmetric_grids: 0,
        nondeterministic: 0,
        disagreements: Vec::new(),
        missing_rows: 0,
    };
    let mut sweeps = Vec::with_capacity(ROCKS);
    let mut ledger = Vec::with_capacity(ROCKS);
    for index in 0..ROCKS {
        let rock = Rock::at(index);
        let name = rock.name();
        assert!(
            rock.domain_is_symmetric(),
            "{name}: the domain is not the symmetric cube this fixture assumes"
        );
        let done = sweep(&rock, &elements, &mut mc);
        if !done.facts.grid_symmetric {
            check.asymmetric_grids += 1;
            error!("{name}: the {SAMPLES}^3 grid does not mirror bit-exactly");
        }
        if !done.identity_exact {
            check.nondeterministic += 1;
            error!("{name}: the identity element is not exact");
        }
        let row = LedgerRow::load(name);
        match row {
            Some(r) => {
                for (column, expected, measured) in r.against(&done) {
                    if expected != measured {
                        check.disagreements.push((name, column, expected, measured));
                        error!("{name}: {column} expected {expected}, measured {measured}");
                    }
                }
                info!(
                    "{:<15} {:>5}/{:<5} {:>5}/{:<5} {:>4}/{:<4} {:>4}/{:<4} {:>3}/{:<3} {:>4}/{:<4} {:>6.0}  raw {} can_fail {}",
                    name,
                    r.cut_edges,
                    done.facts.cut_edges,
                    r.order_sensitive_edges,
                    done.facts.order_sensitive_edges,
                    r.worst_differing_vertices,
                    done.worst_differing_vertices,
                    r.elements_vertex_exact,
                    done.exact_total(),
                    r.pure_permutation_exact,
                    done.pure_permutation_exact,
                    r.worst_component_ulp,
                    done.worst_component_ulp,
                    done.wall_ms,
                    r.elements_vertex_exact_raw,
                    r.fixture_can_fail,
                );
            }
            None => {
                check.missing_rows += 1;
                error!("{name}: p-57.csv has no marching_cubes 33^3 row");
            }
        }
        info!(
            "{name}: exact {} of 48 -- det+ {} of 24, det- {} of 24, negate-nothing {} of 6, \
             perm=012 {} of 8, first failure {}",
            done.exact_total(),
            done.det_plus_exact,
            done.det_minus_exact,
            done.pure_permutation_exact,
            done.pure_sign_flip_exact,
            done.first_failing
                .map_or_else(|| String::from("none"), |i| elements[i].short())
        );
        sweeps.push(done);
        ledger.push(row);
    }
    if check.disagreements.is_empty() && check.missing_rows == 0 {
        info!("every committed number reproduced: {ROCKS} fields x 9 columns of p-57.csv");
    }

    commands.insert_resource(Rig {
        mc,
        buffer: MeshBuffer::<f64>::new(),
        meshes: slot_meshes,
        rock: linear(ROCK_SRGB),
        moved: linear(MOVED_SRGB),
    });
    commands.insert_resource(Group(elements));
    commands.insert_resource(Sweeps(sweeps));
    commands.insert_resource(Ledger(ledger));
    commands.insert_resource(check);
    commands.insert_resource(Steer {
        pinned: None,
        markers: true,
    });
    commands.insert_resource(Shot {
        field: flags.field % ROCKS,
        element: MIRROR,
    });
}

/// World `x` of a slot centre.
fn slot_x(slot: usize) -> f32 {
    (slot as f32 - 1.0) * SLOT_GAP
}

/// The keys this example owns, on top of the harness's.
fn steer(keys: Res<ButtonInput<KeyCode>>, shot: Res<Shot>, mut steer: ResMut<Steer>) {
    let from = steer.pinned.unwrap_or(shot.element) as i32;
    let step = |by: i32| Some((from + by).rem_euclid(GROUP_ORDER as i32) as usize);
    if keys.just_pressed(KeyCode::ArrowRight) {
        steer.pinned = step(1);
    }
    if keys.just_pressed(KeyCode::ArrowLeft) {
        steer.pinned = step(-1);
    }
    if keys.just_pressed(KeyCode::KeyA) {
        steer.pinned = None;
    }
    if keys.just_pressed(KeyCode::KeyH) {
        steer.markers = !steer.markers;
    }
}

/// `ISOMESH_CAPTURE_FRAMES`, or the harness default.
///
/// Read from the environment rather than from [`Capture`], which keeps its
/// length private, because pacing the tour off the captured-frame count is what
/// stops a six-frame smoke test and a 120-frame clip from both being a still.
fn capture_frames() -> u32 {
    std::env::var("ISOMESH_CAPTURE_FRAMES")
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .filter(|n| *n > 0)
        .unwrap_or(60)
}

/// Decide which field and which element this frame is about.
///
/// Under capture the tour advances with the captured frame count, so a clip of
/// any length is the whole tour rather than a slice of it — the mistake
/// `scripts/record_all_gifs.sh` documents three times. Interactively it runs on
/// wall-clock time and loops, and `Right`/`Left` pin an element.
fn advance(
    time: Res<Time>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    steer: Res<Steer>,
    mut shot: ResMut<Shot>,
    mut elapsed: Local<f32>,
) {
    let element = match steer.pinned {
        Some(pinned) => pinned,
        None => {
            let phase = if capture.is_active() {
                f32::from(u16::try_from(capture.taken).unwrap_or(u16::MAX))
                    / f32::from(u16::try_from(capture_frames()).unwrap_or(1).max(1))
            } else {
                if !flags.paused {
                    *elapsed += time.delta_secs();
                }
                (*elapsed / TOUR_SECONDS).fract()
            }
            .clamp(0.0, 1.0);
            if phase < INTRO_END {
                // Open on the mirror: the case a studio actually hits.
                MIRROR
            } else if phase < WALK_END {
                let walked = (phase - INTRO_END) / (WALK_END - INTRO_END);
                ((walked * GROUP_ORDER as f32) as usize).min(GROUP_ORDER - 1)
            } else {
                // Close on `210/---`, which is a proper rotation and still
                // moves. The tally is at 48 seen by then.
                GROUP_ORDER - 1
            }
        }
    };
    *shot = Shot {
        field: flags.field % ROCKS,
        element,
    };
}

/// Mesh whatever the shot asks for — only when the answer would change.
///
/// Slot 0 is `mesh(f)` and comes straight out of the sweep, so it is never
/// re-extracted. Slot 1 is pinned to [`ROTATION`] and is re-extracted only when
/// the field changes. Slot 2 is the element under test and is the only thing a
/// step costs.
fn rebuild(
    shot: Res<Shot>,
    group: Res<Group>,
    sweeps: Res<Sweeps>,
    mut live: ResMut<Live>,
    mut rig: ResMut<Rig>,
    mut meshes: ResMut<Assets<Mesh>>,
) {
    let field_changed = live.built.is_none_or(|built| built.field != shot.field);
    if !field_changed
        && live
            .built
            .is_some_and(|built| built.element == shot.element)
    {
        return;
    }

    let rock = Rock::at(shot.field);
    let sweep = &sweeps.0[shot.field];
    let scale = sweep.display_scale;

    if field_changed {
        let no_movement = vec![false; sweep.reference_positions.len()];
        let mesh = build_mesh(
            &sweep.reference_positions,
            &sweep.reference_normals,
            &sweep.reference_indices,
            &no_movement,
            scale,
            rig.rock,
            rig.moved,
        );
        install(&mut meshes, &rig.meshes[0], 0, mesh);

        let g = group.0[ROTATION];
        // The extraction time is reported for slot 2, the one a step pays for.
        extract(&rock, &mut rig, g);
        let res = residue(&rig.buffer.positions, &sweep.reference_positions, g);
        let flags = moved_flags(&rig.buffer.positions, &res);
        let mesh = build_mesh(
            &rig.buffer.positions,
            &rig.buffer.normals,
            &rig.buffer.indices,
            &flags,
            scale,
            rig.rock,
            rig.moved,
        );
        install(&mut meshes, &rig.meshes[1], 1, mesh);
        // Slot 1's own residue is not reported anywhere, so a silent failure
        // there would show only as a tinted chunk. Say it out loud instead.
        if !res.exact() {
            error!(
                "{}: {} is not bit-exact -- p-57.csv says every negate-nothing element is",
                rock.name(),
                g.short()
            );
        }
    }

    let g = group.0[shot.element];
    let extract_ms = extract(&rock, &mut rig, g);
    let res = residue(&rig.buffer.positions, &sweep.reference_positions, g);
    let flags = moved_flags(&rig.buffer.positions, &res);
    live.vertices = rig.buffer.positions.len();
    live.triangles = rig.buffer.triangle_count();
    live.moved = res.only_got.len();
    live.worst_ulp = res.worst_ulp;
    live.exact = res.exact();
    live.extract_ms = extract_ms;
    live.markers.clear();
    live.markers.extend(
        rig.buffer
            .positions
            .iter()
            .zip(&flags)
            .filter(|(_, moved)| **moved)
            .map(|(p, _)| {
                Vec3::new(
                    p[0] as f32 * scale + slot_x(2),
                    p[1] as f32 * scale,
                    p[2] as f32 * scale,
                )
            }),
    );
    let mesh = build_mesh(
        &rig.buffer.positions,
        &rig.buffer.normals,
        &rig.buffer.indices,
        &flags,
        scale,
        rig.rock,
        rig.moved,
    );
    install(&mut meshes, &rig.meshes[2], 2, mesh);

    live.built = Some(*shot);
}

/// Extract `mesh(g·f)` into `rig.buffer`, and report how long it took.
fn extract(rock: &Rock, rig: &mut Rig, g: Element) -> f64 {
    let l = rock.half_extent();
    let cell = cell_size(l);
    let shape = RuntimeShape3::new([SAMPLES; 3]).expect("a 33 cube fits u32");
    let wrapped = Rotated {
        field: rock,
        g,
        g_inv: g.inverse(),
    };
    rig.buffer.reset();
    let started = Instant::now();
    if let Err(e) = rig
        .mc
        .extract(&wrapped, &shape, [-l; 3], cell, &mut rig.buffer)
    {
        error!("{}: {} extraction failed: {e}", rock.name(), g.short());
    }
    started.elapsed().as_secs_f64() * 1000.0
}

/// Write a slot's mesh into the asset it already owns.
///
/// `Assets::insert` returns `Err` only when the id's generation is stale, i.e.
/// the asset has been dropped — which cannot happen while [`Rig`] holds a strong
/// handle to it. It is reported rather than discarded because the failure mode
/// of swallowing it is a chunk that silently stops updating, which reads as a
/// wrong verdict rather than as a broken asset.
fn install(meshes: &mut Assets<Mesh>, handle: &Handle<Mesh>, slot: usize, mesh: Mesh) {
    if let Err(e) = meshes.insert(handle, mesh) {
        error!("slot {slot}: the mesh asset could not be replaced: {e}");
    }
}

/// The `f64` extraction as a Bevy mesh, uniformly scaled for display.
///
/// Cast and scaled rather than re-extracted in `f32`: the comparison this whole
/// example is about happened on the `f64` positions, and a second extraction at
/// a different precision would not be the mesh the verdict describes.
fn build_mesh(
    positions: &[[f64; 3]],
    normals: &[[f64; 3]],
    indices: &[u32],
    moved: &[bool],
    scale: f32,
    rock: [f32; 4],
    hot: [f32; 4],
) -> Mesh {
    let mut mesh = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::default(),
    );
    mesh.insert_attribute(
        Mesh::ATTRIBUTE_POSITION,
        positions
            .iter()
            .map(|p| {
                [
                    p[0] as f32 * scale,
                    p[1] as f32 * scale,
                    p[2] as f32 * scale,
                ]
            })
            .collect::<Vec<[f32; 3]>>(),
    );
    mesh.insert_attribute(
        Mesh::ATTRIBUTE_NORMAL,
        normals
            .iter()
            .map(|n| [n[0] as f32, n[1] as f32, n[2] as f32])
            .collect::<Vec<[f32; 3]>>(),
    );
    mesh.insert_attribute(
        Mesh::ATTRIBUTE_COLOR,
        moved
            .iter()
            .map(|m| if *m { hot } else { rock })
            .collect::<Vec<[f32; 4]>>(),
    );
    mesh.insert_indices(Indices::U32(indices.to_vec()));
    mesh
}

/// A three-line cross on every vertex that moved.
///
/// Gizmos draw on top, which is what makes a marker on the far side of the chunk
/// visible — the population is a share of the whole surface, not of its
/// silhouette, and hiding half of it would understate the count the HUD reports.
fn draw_markers(steer: Res<Steer>, live: Res<Live>, mut gizmos: Gizmos) {
    if !steer.markers {
        return;
    }
    let colour = Color::srgb(MOVED_SRGB[0], MOVED_SRGB[1], MOVED_SRGB[2]);
    for p in &live.markers {
        gizmos.line(*p - Vec3::X * MARKER_ARM, *p + Vec3::X * MARKER_ARM, colour);
        gizmos.line(*p - Vec3::Y * MARKER_ARM, *p + Vec3::Y * MARKER_ARM, colour);
        gizmos.line(*p - Vec3::Z * MARKER_ARM, *p + Vec3::Z * MARKER_ARM, colour);
    }
}

// ─── what is on screen ──────────────────────────────────────────────────────

/// The HUD, the three slot labels and the caption, all from one read of the
/// state, so two numbers on screen cannot disagree by a frame.
#[allow(clippy::too_many_arguments)]
fn report(
    live: Res<Live>,
    shot: Res<Shot>,
    group: Res<Group>,
    sweeps: Res<Sweeps>,
    ledger: Res<Ledger>,
    check: Res<SelfCheck>,
    steer: Res<Steer>,
    mut stats: ResMut<DemoStats>,
    mut labels: Query<(&SlotLabel, &mut Text), Without<Caption>>,
    mut caption: Query<&mut Text, With<Caption>>,
) {
    let rock = Rock::at(shot.field);
    let sweep = &sweeps.0[shot.field];
    let g = group.0[shot.element];
    let l = rock.half_extent();
    let seen = shot.element + 1;
    let exact_seen = sweep.exact_among(seen);

    stats.title = String::from("E-314 game mirror dedup -- M-356 / P-57, in an asset cache");
    stats.vertices = live.vertices;
    stats.triangles = live.triangles;
    stats.extract_ms = live.extract_ms;

    // Eleven lines and three blanks, and the count is a constraint rather than a
    // taste: `HUD_PANEL` is sized for 25 rendered lines and the chunks start
    // 4 px below it.
    let verdict = if live.exact { "identical" } else { "MOVED" };
    let mut extra = vec![
        format!(
            "field       {:<15} {SAMPLES}^3  cell {:.9}  origin {:.9}",
            rock.name(),
            cell_size(l),
            -l
        ),
        format!(
            "cut edges   {:>7}  order-sensitive {:>5}  h = L/16, no extractor",
            commas(sweep.facts.cut_edges),
            commas(sweep.facts.order_sensitive_edges)
        ),
        format!(
            "\nelement {:>2} of 48   {:<18} det {:+}   {}",
            shot.element + 1,
            g.label(),
            g.det(),
            g.what()
        ),
        format!(
            "verdict     {verdict:<10} {:>6} of {} vertices   worst gap {} ULP",
            commas(live.moved),
            commas(live.vertices),
            live.worst_ulp
        ),
    ];

    // The mechanism, stated as the equality it is. Quoted against the sweep's
    // worst element rather than this one, because a negate-nothing element moves
    // nothing and `0 == 643` is not a disagreement.
    extra.push(format!(
        "mechanism   worst moved {} == order-sensitive edges {}, {}",
        sweep.worst_differing_vertices,
        sweep.facts.order_sensitive_edges,
        if sweep.worst_differing_vertices == sweep.facts.order_sensitive_edges {
            "to the unit"
        } else {
            "WHICH DISAGREE"
        }
    ));
    extra.push(format!(
        "\ntally       {exact_seen} of {seen} walked bit-identical   {} of 48 over the group",
        sweep.exact_total()
    ));
    extra.push(format!(
        "            det +1  {} of 24        det -1  {} of 24",
        sweep.det_plus_exact, sweep.det_minus_exact
    ));
    // The split that says the predicate is negation and not handedness. The two
    // named counterexamples are in the caption, where a viewer reads them.
    extra.push(format!(
        "            negate nothing  {} of 6       perm=012  {} of 8",
        sweep.pure_permutation_exact, sweep.pure_sign_flip_exact
    ));

    match ledger.0[shot.field] {
        Some(row) => {
            let mismatches = row
                .against(sweep)
                .into_iter()
                .filter(|(_, e, m)| e != m)
                .count();
            extra.push(format!(
                "\np-57.csv / run   order_sensitive_edges {:>5} / {:<5}  exact {:>2} / {}",
                row.order_sensitive_edges,
                sweep.facts.order_sensitive_edges,
                row.elements_vertex_exact,
                sweep.exact_total()
            ));
            extra.push(format!(
                "                 worst_differing_verts {:>5} / {:<5}  {}",
                row.worst_differing_vertices,
                sweep.worst_differing_vertices,
                if mismatches == 0 {
                    "agrees on all 9 columns"
                } else {
                    "DISAGREES"
                }
            ));
            if !row.fixture_can_fail {
                extra.push(String::from(
                    "                 fixture_can_fail = false: 48 of 48 is a control here",
                ));
            }
        }
        None => extra.push(String::from(
            "\nSELF CHECK FAILED: p-57.csv has no marching_cubes 33^3 row here",
        )),
    }

    if check.asymmetric_grids > 0 {
        extra.push(format!(
            "SELF CHECK FAILED: {} fields have a grid that does not mirror bit-exactly",
            check.asymmetric_grids
        ));
    }
    if check.nondeterministic > 0 {
        extra.push(format!(
            "SELF CHECK FAILED: the identity element moved on {} fields",
            check.nondeterministic
        ));
    }
    if !check.disagreements.is_empty() {
        extra.push(format!(
            "SELF CHECK FAILED: {} committed numbers did not reproduce -- see the log",
            check.disagreements.len()
        ));
    }
    // No blank line before this one: the panel is sized for 25 rendered lines
    // and `fixture_can_fail = false` adds a twelfth on `box_exact`.
    extra.push(format!(
        "[Right/Left] step element   [A] auto   [H] markers {}",
        if steer.markers { "on" } else { "off" }
    ));
    stats.extra = extra;

    let rotation = group.0[ROTATION];
    let slot_text = [
        format!("original\nmesh(f) -- {} vertices", commas(sweep.vertices)),
        format!("reused: {}\nidentical -- hash matches", rotation.short()),
        format!(
            "reused: {}\n{}",
            g.short(),
            if live.exact {
                String::from("identical -- hash matches")
            } else {
                format!("MOVED -- {} differ", commas(live.moved))
            }
        ),
    ];
    for (slot, mut text) in &mut labels {
        if let Some(line) = slot_text.get(slot.0) {
            text.0.clone_from(line);
        }
    }

    // Four captions, because the four quadrants of (determinant, verdict) are
    // four different sentences and two of them are the counterexamples that kill
    // "reflections break it, rotations do not".
    let line = match (g.det() > 0, live.exact) {
        (false, true) => format!(
            "{} is a REFLECTION and the mesh is bit-identical -- the cache hits.\n\
             Handedness is not the test: {} of the {} that dedup are reflections.",
            g.short(),
            sweep.det_minus_exact,
            sweep.exact_total()
        ),
        (true, true) => format!(
            "{} reorients the chunk and the mesh is bit-identical -- the cache hits.\n\
             {} of the {seen} walked so far are, and every one of them negates no axis.",
            g.short(),
            exact_seen
        ),
        (true, false) => format!(
            "{} is a proper ROTATION and it still breaks the hash: {} vertices moved.\n\
             Only {} of 48 dedup -- the test is negation, not handedness.",
            g.short(),
            commas(live.moved),
            sweep.exact_total()
        ),
        (false, false) => format!(
            "{} looks like the same rock and hashes differently: {} of {} moved.\n\
             Only {} of 48 orientations dedup: exactly the {} that negate no axis.",
            g.short(),
            commas(live.moved),
            commas(live.vertices),
            sweep.exact_total(),
            sweep.pure_permutation_exact
        ),
    };
    for mut target in &mut caption {
        target.0.clone_from(&line);
    }
}
