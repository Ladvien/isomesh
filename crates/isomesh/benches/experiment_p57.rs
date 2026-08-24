//! **P-57 — octahedral equivariance of the seven extractors, to the bit.**
//!
//! Ticket: R-055. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p57
//! ```
//!
//! Writes `docs/experiments/p-57.csv`.
//!
//! # The relation, and the level it is stated at
//!
//! Every element of the 48-element octahedral group `B₃` is a **signed
//! coordinate permutation**: three components reordered and each optionally
//! negated. In `f64` that is exact — a permutation moves bits and a negation
//! flips one bit, and neither rounds. So `mesh(g·f)` and `g·mesh(f)` are two
//! `f64` meshes that a *correct, axis-independent* extractor must produce
//! bit-for-bit identically, and any difference between them is arithmetic that
//! depended on which axis happens to be called `x`.
//!
//! The comparison is at **vertex level**: sorted multisets of vertex positions,
//! each component reduced to its `f64` bit pattern. Not index buffers, and not
//! triangles as the primary statement, because
//! `docs/research/2026-08-23-discovery-dossier.md:267` already worked out why a
//! triangle-level relation is known-false in advance —
//! `marching_cubes/table.rs` picks its `safe_apex` by **lowest edge index**, and
//! edge indices are a naming of the axes rather than a property of the geometry.
//! Stating the relation at triangle level therefore "manufactures 2,688 false
//! positives". C3 exists to *quantify* that warning rather than repeat it, so
//! the triangle-level count is measured beside the vertex-level one on every
//! row and neither is allowed to stand in for the other.
//!
//! # THE FIXTURE DECISION — the single most important paragraph here
//!
//! The relation is only about the extractor if the **sample grid is exactly
//! closed under `g`**. A grid `origin + i·h` is closed under a sign flip only
//! when `origin + i·h == −(origin + (n−1−i)·h)` **bit-exactly**, which needs
//! `origin = −((n−1)/2)·h` with `n` odd, and needs `h` to be a binary fraction
//! so that `i·h` is exact for every `i`.
//!
//! The crate's own spacing is `2L/(n−1)`. That is a binary fraction at
//! `n = 17, 33, 65` — `L/8`, `L/16`, `L/32` — and it is **not** one at
//! `n = 25`, where `2L/24 = L/12` has a factor of three in the denominator and
//! **16 of the 25 coordinates per axis fail the bit-exact mirror test**. Running
//! at the crate's 25³ spacing would report a falsification of the hypothesis
//! that was really a property of the fixture, which is the exact inversion this
//! experiment exists to avoid. So the two fixtures are:
//!
//! - **`samples_per_axis = 33`** — `cell_size = L/16`, `origin = −L`. The
//!   crate's own domain, its own spacing, unaltered.
//! - **`samples_per_axis = 25`** — `cell_size = 3L/32`,
//!   `origin = −12·cell_size = −1.125·L`. A box 12.5% wider than the field's
//!   own. `3L/32` is a binary fraction for every `L` here (`0.1875` at `L = 2`,
//!   `0.65625` at `L = 7`, `0.75` at `L = 8`) so the grid is exactly mirror
//!   symmetric — **and** it is not of the form `L/2^k`, so crossing
//!   coordinates are generically non-representable.
//!
//! That second property is the whole reason the 25³ fixture exists. M-178
//! recorded the vacuous-fixture trap: at 17³ and 33³ the crate's spacing makes
//! `box_exact` and `thin_plate` crossings land on exactly representable
//! coordinates, so those two fields *cannot move* and a pass there means
//! nothing. `grid_symmetric` is asserted `true` on every row and
//! `fixture_can_fail` is asserted `true` somewhere, and a run where no fixture
//! can fail prints a VOID banner and is worthless regardless of its pass rate.
//!
//! `fixture_can_fail` is measured independently of any extractor. Every
//! axis-aligned grid edge that changes sign is interpolated **twice** — forward
//! as `p_lo + (p_hi − p_lo)·(a/(a−b))` and from the far end as
//! `p_hi + (p_lo − p_hi)·(b/(b−a))` — and an edge whose two answers differ in a
//! single bit is an edge where the *order the endpoints are visited in* decides
//! the vertex. `order_sensitive_edges` counts them and `cut_edges` is the
//! denominator.
//!
//! # Signed zero is a representation, not an error, and it is normalised
//!
//! Both fixtures have an odd sample count and are centred, so the coordinate
//! `0.0` is on the grid — and `−(0.0)` is `−0.0`, whose bit pattern is not
//! `0.0`'s. Left alone, every sign-flipping element would "fail" on every field
//! for a reason that has nothing to do with any extractor: the two sides agree
//! to infinite precision and disagree only about which of the two encodings of
//! zero they wrote down.
//!
//! So the key normalises `−0.0` to `0.0`, identically on both sides, and says
//! so. This is not a tolerance and not a fallback — no other value is touched,
//! and a one-ULP difference anywhere else is still a failure. To keep the
//! decision visible rather than buried, the un-normalised counts are recorded
//! beside the registered ones as `elements_vertex_exact_raw` and
//! `elements_triangle_exact_raw`: a reader can see exactly what the
//! normalisation bought, on every row, and `grid_zero_coordinates` says how many
//! grid coordinates it applies to.
//!
//! # What this crate owns rather than the source
//!
//! - **The group.** All 48 elements, generated as the 6 permutations of
//!   `[0, 1, 2]` crossed with the 8 sign patterns, and **not** filtered on
//!   determinant. The generator at `dual_contouring/solve/tests.rs:130-165`
//!   keeps only the 24 with `det = +1`, because a QEF vertex rule is a statement
//!   about rotations; bit-exactness of `a/(a−b)` is a statement about
//!   arithmetic, and a reflection is exactly as exact as a rotation. Signs are
//!   stored as `i8` and applied by negation rather than by multiplying by
//!   `±1.0`, so no float arithmetic enters the group action at all.
//! - **The inverse.** `apply(p)[k] = ±p[perm[k]]` inverts to
//!   `inv.perm[perm[k]] = k`, `inv.sign[perm[k]] = sign[k]`. Checked, not
//!   assumed: `g⁻¹(g(p))` must be bit-identical to `p` on a probe set that
//!   includes both zeros, and the 48 elements must be pairwise distinct.
//! - **The rotated field.** `Rotated` is `g·f`, i.e.
//!   `sample(p) = f(g⁻¹·p)`, and it forwards `gradient` as well —
//!   `∇(f∘g⁻¹)(p) = g·(∇f)(g⁻¹·p)` for orthogonal `g`. Forwarding only
//!   `sample` would silently substitute `Sdf`'s central-difference default and
//!   change what the three dual extractors see, which would measure the
//!   differencing stencil rather than the vertex rule.
//!
//! # Which columns decide which clause
//!
//! - **C1** — `elements_vertex_exact == 48` on all 64 rows with
//!   `family = primal`.
//! - **C2** — `elements_vertex_exact < 48` on the 48 rows with `family = dual`,
//!   with the count and the identity of the failures in
//!   `48 − elements_vertex_exact`, `vertex_failing_labels` and
//!   `vertex_failing_mask` (a 48-bit hex mask, bit `i` = element `i` of the
//!   printed index table).
//! - **C3** — `elements_vertex_exact_triangle_fail`, the number of elements that
//!   are vertex-exact and **not** triangle-exact. Non-zero anywhere is C3
//!   holding; zero everywhere would mean `safe_apex` is invariant after all.
//!
//! # Sizing a failure, and naming its mechanism
//!
//! A count of failing elements says nothing about how badly they failed, so two
//! columns size it. `worst_differing_vertices` is the largest number of
//! vertices that moved, over all 48 elements, and `worst_component_ulp` is the
//! largest ULP gap across the pairing of what moved with what it replaced —
//! measured on the sign-magnitude-ordered integer image of the bits, so it is
//! meaningful across zero.
//!
//! Both are read off the **multiset symmetric difference** of the two sorted
//! key lists rather than a positional walk over them, and that distinction is
//! not cosmetic: the first run of this harness walked positionally and reported
//! `worst_component_ulp = 9.2e18` — the gap between `+2` and `−2` — for meshes
//! that differ by one bit on 72 edges, because a single *inserted* vertex
//! shifts every later entry against its neighbour. The merge that fixes it is
//! `multiset_difference` below.
//!
//! **`worst_component_ulp` must be read next to `worst_differing_vertices`,
//! and it is not a tolerance either way.** When the residue is small against
//! `vertices`, the two residues line up one-to-one and the number is the
//! perturbation: the marching-cubes rows have
//! `worst_differing_vertices == order_sensitive_edges` exactly, and the ULP gap
//! there is the size of the `a/(a−b)` versus `b/(b−a)` disagreement on those
//! edges. When the residue is most of the mesh — which is what the
//! tetrahedral extractors do, because a six-tetrahedron cell decomposition
//! is not octahedrally symmetric and its diagonals cut different edges — the
//! two meshes have different *vertex sets* rather than perturbed ones, the
//! sorted pairing is between unrelated vertices, and the number degenerates to
//! the domain width in ULP. That is a "different mesh" flag, not a rounding
//! measurement, and the neighbouring count is what says which of the two it is.
//!
//! Two more columns split the mechanism. `pure_permutation_exact` counts the 6
//! elements with `sign = +++`, which relabel the axes and nothing else, and
//! `pure_sign_flip_exact` counts the 8 with `perm = 012`, which additionally
//! reverse the order every grid edge is interpolated in. A rule that fails only
//! the sign flips is an edge-interpolation asymmetry; a rule that fails
//! relabellings is an asymmetric cell decomposition. `det_plus_vertex_exact`
//! and `det_minus_vertex_exact` split the same 48 the other way.
//!
//! # The negative control, and the one that is not negotiable
//!
//! The identity element is **in** the group and is element 0. It must be
//! vertex-exact and triangle-exact on every row, and it is asserted. It goes
//! through the `Rotated` wrapper like every other element, so it checks three
//! things at once: that the extractor is deterministic, that the wrapper is
//! transparent, and that the comparison is not accidentally comparing a mesh to
//! itself by some other route. A failure there voids the row rather than
//! falsifying anything, which is why it is an assertion and not a column.
//!
//! `fixture_can_fail` is the other control and points the other way: it is the
//! check that a pass *could* have been a failure. `boundary_sign_constant` is
//! recorded but not asserted — it is `false` for `fbm_terrain` by construction,
//! which exits through the sides of its own domain.
//!
//! # Counted, not timed
//!
//! Every registered metric is an integer or a bit pattern and is identical on
//! every machine. `wall_ms` sits beside them and **gates nothing**: M-348 is the
//! incident where a discovery was demoted for resting on a wall clock.
//!
//! One more format note, because it is load-bearing: `common::experiment` writes
//! CSV without quoting, so no recorded value may contain a comma.
//! `first_failing_element` is therefore `perm=201;sign=+-+` rather than the
//! comma form, and label lists are joined with `|`.

mod common;

use std::time::Instant;

use isomesh::extractor::Extractor;
use isomesh::fields::ReferenceField;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// The 6 axis permutations. Crossed with 8 sign patterns this is all 48.
const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

/// The order of the octahedral group, and the value `elements_tested` must
/// report on every row.
const GROUP_ORDER: usize = 48;

/// Probe points for the inverse round-trip check.
///
/// Includes both zeros deliberately: `−(−0.0)` is `0.0`, so a signed
/// permutation round-trips a zero bit-exactly even though it does not preserve
/// it, and that is worth pinning rather than assuming.
const PROBES: [[f64; 3]; 4] = [
    [0.3, -1.7, 2.9],
    [0.0, -0.0, 1.5],
    [-2.25, 0.656_25, -0.187_5],
    [1.0, 1.0, 1.0],
];

// ─── the group ──────────────────────────────────────────────────────────────

/// One element of the octahedral group, as a signed axis permutation.
///
/// `apply(p)[k] = sign[k] · p[perm[k]]`. Signs are `i8` and applied by
/// negation, so the action is bit-exact by construction rather than by an
/// argument about multiplying by `±1.0`.
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

    /// `perm=201;sign=+-+`. Semicolons, because the CSV writer does not quote.
    fn label(self) -> String {
        let mut s = String::from("perm=");
        for k in 0..3 {
            s.push(char::from(b'0' + self.perm[k] as u8));
        }
        s.push_str(";sign=");
        for k in 0..3 {
            s.push(if self.sign[k] < 0 { '-' } else { '+' });
        }
        s
    }

    /// `201/+-+`, for packing many labels into one cell.
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
}

/// All 48, in a fixed order: permutation outer, sign pattern inner, so element
/// 0 is the identity and the printed index table is stable across runs.
///
/// **Not filtered on determinant.** The 24 reflections are as exact as the 24
/// rotations and are half the group the hypothesis names.
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

/// The group is checked before it is used.
fn verify_group(g: &[Element]) {
    assert_eq!(g.len(), GROUP_ORDER, "the octahedral group has 48 elements");
    assert!(
        g[0] == Element {
            perm: [0, 1, 2],
            sign: [1, 1, 1]
        },
        "element 0 must be the identity: it is the negative control"
    );
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
            d => panic!("{} has determinant {d}, not ±1", e.label()),
        }
        let inv = e.inverse();
        for p in PROBES {
            let round = inv.apply(e.apply(p));
            for k in 0..3 {
                assert_eq!(
                    round[k].to_bits(),
                    p[k].to_bits(),
                    "{} does not round-trip {p:?} bit-exactly: got {round:?}",
                    e.label()
                );
            }
        }
    }
    assert_eq!(rotations, 24, "24 elements must have det = +1");
    assert_eq!(reflections, 24, "24 elements must have det = -1");
}

// ─── the rotated field ──────────────────────────────────────────────────────

/// `g·f`, i.e. the field `f` pushed forward by `g`: `(g·f)(p) = f(g⁻¹·p)`.
///
/// Both `g` and `g⁻¹` are stored because both are needed and neither costs
/// anything: the sample point goes in through `g⁻¹` and the gradient comes back
/// out through `g`.
struct Rotated<'a, S> {
    /// The field being rotated.
    field: &'a S,
    /// The element.
    g: Element,
    /// Its inverse, precomputed.
    g_inv: Element,
}

impl<S: Sdf<Scalar = f64>> Sdf for Rotated<'_, S> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        self.field.sample(self.g_inv.apply(p))
    }

    /// `∇(f∘g⁻¹)(p) = (g⁻¹)ᵀ·(∇f)(g⁻¹·p)`, and `(g⁻¹)ᵀ = g` because `g` is
    /// orthogonal.
    ///
    /// Overridden rather than inherited. The default is six `sample` calls, and
    /// letting the dual extractors see a central difference here would measure
    /// the stencil instead of the vertex rule.
    #[inline]
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        self.g.apply(self.field.gradient(self.g_inv.apply(p)))
    }
}

// ─── bit keys ───────────────────────────────────────────────────────────────

/// The bit pattern of `−0.0`.
const NEGATIVE_ZERO: u64 = 1u64 << 63;

/// A position component as a comparison key.
///
/// With `normalise`, `−0.0` folds onto `0.0` and nothing else moves — see the
/// module header. Without it, the raw bits.
#[inline]
fn key(v: f64, normalise: bool) -> u64 {
    let b = v.to_bits();
    if normalise && b == NEGATIVE_ZERO {
        0
    } else {
        b
    }
}

/// The sign-magnitude-ordered integer image of an `f64`'s bits.
///
/// Monotone in the value and continuous across zero, so a difference of these
/// is a ULP count that means something for a pair straddling zero. Both zeros
/// map to `0`.
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

/// The multiset symmetric difference of two **sorted** key lists, as
/// `(only_in_a, only_in_b)`.
///
/// This exists because the obvious thing is wrong. Walking two sorted lists
/// positionally and diffing entry `i` against entry `i` reports nonsense the
/// moment the lists differ by an *insertion* rather than a perturbation: every
/// entry after the first difference is compared against its neighbour, and the
/// first run of this harness duly reported `worst_component_ulp` of
/// `9.2e18` — the gap between `+2` and `−2` — for meshes that actually differ
/// by one bit on a handful of vertices.
///
/// A sorted merge instead cancels every vertex the two meshes agree on and
/// leaves exactly the ones that moved. Those two residues have the same length
/// whenever the vertex counts match, and pairing them in sorted order pairs
/// each moved vertex with the position it moved from, because a one-ULP
/// perturbation cannot reorder a vertex past an unrelated one.
fn multiset_difference(a: &[[u64; 3]], b: &[[u64; 3]]) -> (Vec<[u64; 3]>, Vec<[u64; 3]>) {
    let mut only_a = Vec::new();
    let mut only_b = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                only_a.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                only_b.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    only_a.extend_from_slice(&a[i..]);
    only_b.extend_from_slice(&b[j..]);
    (only_a, only_b)
}

/// Sorted multiset of vertex positions as bit triples, optionally mapped
/// through a group element first.
fn vertex_keys(positions: &[[f64; 3]], g: Option<Element>, normalise: bool) -> Vec<[u64; 3]> {
    let mut out: Vec<[u64; 3]> = positions
        .iter()
        .map(|p| {
            let q = match g {
                Some(e) => e.apply(*p),
                None => *p,
            };
            std::array::from_fn(|k| key(q[k], normalise))
        })
        .collect();
    out.sort_unstable();
    out
}

/// Sorted list of triangles, each as the bit triple of its three **positions**
/// rotated to its lexicographically smallest rotation.
///
/// Rotations only, never a full sort of the three corners: a rotation preserves
/// winding, so a triangle whose normal has been flipped stays distinguishable
/// from one that has not. Indices are deliberately not part of the key — an
/// index buffer is a naming and this relation is about geometry.
fn triangle_keys(
    positions: &[[f64; 3]],
    indices: &[u32],
    g: Option<Element>,
    normalise: bool,
) -> Vec<[[u64; 3]; 3]> {
    let mapped = |i: u32| -> [u64; 3] {
        let p = positions[i as usize];
        let q = match g {
            Some(e) => e.apply(p),
            None => p,
        };
        std::array::from_fn(|k| key(q[k], normalise))
    };
    let mut out: Vec<[[u64; 3]; 3]> = Vec::with_capacity(indices.len() / 3);
    for tri in indices.chunks_exact(3) {
        let c = [mapped(tri[0]), mapped(tri[1]), mapped(tri[2])];
        let rotations = [[c[0], c[1], c[2]], [c[1], c[2], c[0]], [c[2], c[0], c[1]]];
        let best = rotations
            .into_iter()
            .min()
            .expect("three rotations is not empty");
        out.push(best);
    }
    out.sort_unstable();
    out
}

// ─── fixtures ───────────────────────────────────────────────────────────────

/// A grid: samples per axis, spacing, and the common origin component.
#[derive(Clone, Copy)]
struct Fixture {
    /// Samples per axis. Odd, so the grid can be centred.
    samples: u32,
    /// Spacing. A binary fraction on every field, which is the point.
    cell_size: f64,
    /// `origin[0] == origin[1] == origin[2]`, because the grid is a cube
    /// centred on the origin.
    origin: f64,
}

/// The two fixtures for a field of half-extent `L`. See the module header —
/// these two values are the experiment's most load-bearing choice.
fn fixtures(l: f64) -> [Fixture; 2] {
    let h33 = l / 16.0;
    let h25 = 3.0 * l / 32.0;
    [
        Fixture {
            samples: 33,
            cell_size: h33,
            origin: -l,
        },
        Fixture {
            samples: 25,
            cell_size: h25,
            origin: -12.0 * h25,
        },
    ]
}

/// What is true of a `(field, fixture)` pair before any extractor runs.
struct FixtureFacts {
    /// `pos[i] == −pos[n−1−i]` bit-exactly, for every `i` and every axis.
    grid_symmetric: bool,
    /// Grid coordinates that are exactly zero. One, for a centred odd grid.
    zero_coordinates: usize,
    /// Axis-aligned grid edges that change sign.
    cut_edges: usize,
    /// Of those, the ones whose crossing coordinate depends on which end it is
    /// interpolated from.
    order_sensitive_edges: usize,
    /// Whether every sample on the grid box's boundary has the same sign.
    boundary_sign_constant: bool,
}

/// The axis coordinate list, built exactly as the extractors build it —
/// `origin + cell_size · i`, per `marching_cubes/mod.rs:229`.
fn axis_coords(fx: &Fixture) -> Vec<f64> {
    (0..fx.samples)
        .map(|i| fx.origin + fx.cell_size * f64::from(i))
        .collect()
}

/// Measure a `(field, fixture)` pair: grid symmetry, and whether the fixture is
/// able to fail at all.
fn fixture_facts<S: Sdf<Scalar = f64>>(field: &S, fx: &Fixture) -> FixtureFacts {
    let n = fx.samples as usize;
    let coords = axis_coords(fx);

    let mut grid_symmetric = true;
    let mut zero_coordinates = 0;
    for i in 0..n {
        if key(coords[i], true) != key(-coords[n - 1 - i], true) {
            grid_symmetric = false;
        }
        if key(coords[i], true) == 0 {
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

    // Boundary sign constancy. Recorded, never asserted: `fbm_terrain` exits
    // through the sides of its own domain and this is `false` for it by
    // construction.
    let mut boundary_sign_constant = true;
    let outside = at(0, 0, 0) < 0.0;
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let on_boundary =
                    x == 0 || y == 0 || z == 0 || x == n - 1 || y == n - 1 || z == n - 1;
                if on_boundary && (at(x, y, z) < 0.0) != outside {
                    boundary_sign_constant = false;
                }
            }
        }
    }

    // The order-sensitivity census. One pass over every axis-aligned grid edge.
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
                    let fwd = lo_c + (hi_c - lo_c) * (va / (va - vb));
                    let rev = hi_c + (lo_c - hi_c) * (vb / (vb - va));
                    if key(fwd, true) != key(rev, true) {
                        order_sensitive_edges += 1;
                    }
                }
            }
        }
    }

    FixtureFacts {
        grid_symmetric,
        zero_coordinates,
        cut_edges,
        order_sensitive_edges,
        boundary_sign_constant,
    }
}

// ─── the measurement ────────────────────────────────────────────────────────

/// Everything one `(field, extractor, fixture)` row reports.
struct Measured {
    /// Vertices in the reference mesh.
    vertices: usize,
    /// Triangles in the reference mesh.
    triangles: usize,
    /// Elements whose vertex multiset matched, with `−0.0` normalised.
    vertex_exact: usize,
    /// The same without the normalisation.
    vertex_exact_raw: usize,
    /// Elements whose canonical triangle list matched.
    triangle_exact: usize,
    /// The same without the normalisation.
    triangle_exact_raw: usize,
    /// Elements that were vertex-exact and not triangle-exact. **C3.**
    vertex_exact_triangle_fail: usize,
    /// Elements whose mesh had a different vertex count entirely.
    vertex_count_mismatches: usize,
    /// Rotations (`det = +1`) that were vertex-exact, out of 24.
    det_plus_vertex_exact: usize,
    /// Reflections (`det = −1`) that were vertex-exact, out of 24.
    det_minus_vertex_exact: usize,
    /// Pure axis relabellings (`sign = +++`) that were vertex-exact, out of 6.
    ///
    /// Split out from the sign flips because the two mechanisms are different:
    /// a relabelling changes which axis a tetrahedral decomposition cuts along,
    /// while a sign flip additionally reverses the order every grid edge is
    /// interpolated in.
    pure_permutation_exact: usize,
    /// Pure sign flips (`perm = 012`) that were vertex-exact, out of 8.
    pure_sign_flip_exact: usize,
    /// The largest number of vertices that moved, over all 48 elements.
    worst_differing_vertices: usize,
    /// Lowest-indexed vertex-level failure, or `"none"`.
    first_failing_element: String,
    /// Its determinant, or `0`.
    first_failing_det: i32,
    /// Largest ULP gap between a mismatching component pair.
    worst_component_ulp: u128,
    /// Bit `i` set = element `i` failed at vertex level.
    vertex_failing_mask: u64,
    /// Bit `i` set = element `i` failed at triangle level.
    triangle_failing_mask: u64,
    /// Short labels of the vertex-level failures, joined with `|`.
    vertex_failing_labels: String,
    /// Wall time for the 49 extractions. Gates nothing.
    wall_ms: u128,
}

/// Extract the reference once, then all 48 rotated meshes, and compare.
///
/// # Panics
///
/// On any extraction error, and on the identity element failing — which means
/// the extractor is not deterministic and voids the row rather than falsifying
/// anything.
fn measure<S, E>(
    field: &S,
    extractor: &mut E,
    fx: &Fixture,
    elements: &[Element],
    reference: &mut MeshBuffer<f64>,
    rotated: &mut MeshBuffer<f64>,
) -> Measured
where
    S: Sdf<Scalar = f64>,
    E: Extractor<f64>,
{
    let shape = RuntimeShape3::new([fx.samples; 3]).expect("fixture grid fits u32");
    let origin = [fx.origin; 3];
    let started = Instant::now();

    reference.reset();
    extractor
        .extract_into(field, &shape, origin, fx.cell_size, reference)
        .expect("reference extraction");

    let mut m = Measured {
        vertices: reference.positions.len(),
        triangles: reference.triangle_count(),
        vertex_exact: 0,
        vertex_exact_raw: 0,
        triangle_exact: 0,
        triangle_exact_raw: 0,
        vertex_exact_triangle_fail: 0,
        vertex_count_mismatches: 0,
        det_plus_vertex_exact: 0,
        det_minus_vertex_exact: 0,
        pure_permutation_exact: 0,
        pure_sign_flip_exact: 0,
        worst_differing_vertices: 0,
        first_failing_element: String::from("none"),
        first_failing_det: 0,
        worst_component_ulp: 0,
        vertex_failing_mask: 0,
        triangle_failing_mask: 0,
        vertex_failing_labels: String::new(),
        wall_ms: 0,
    };
    let mut failing: Vec<String> = Vec::new();

    for (index, &g) in elements.iter().enumerate() {
        let wrapped = Rotated {
            field,
            g,
            g_inv: g.inverse(),
        };
        rotated.reset();
        extractor
            .extract_into(&wrapped, &shape, origin, fx.cell_size, rotated)
            .expect("rotated extraction");

        let got = vertex_keys(&rotated.positions, None, true);
        let want = vertex_keys(&reference.positions, Some(g), true);
        let vertex_ok = got == want;

        let got_raw = vertex_keys(&rotated.positions, None, false);
        let want_raw = vertex_keys(&reference.positions, Some(g), false);
        if got_raw == want_raw {
            m.vertex_exact_raw += 1;
        }

        let got_tri = triangle_keys(&rotated.positions, &rotated.indices, None, true);
        let want_tri = triangle_keys(&reference.positions, &reference.indices, Some(g), true);
        let triangle_ok = got_tri == want_tri;

        let got_tri_raw = triangle_keys(&rotated.positions, &rotated.indices, None, false);
        let want_tri_raw = triangle_keys(&reference.positions, &reference.indices, Some(g), false);
        if got_tri_raw == want_tri_raw {
            m.triangle_exact_raw += 1;
        }

        if got.len() != want.len() {
            m.vertex_count_mismatches += 1;
        }
        let pure_permutation = g.sign == [1, 1, 1];
        let pure_sign_flip = g.perm == [0, 1, 2];
        if vertex_ok {
            m.vertex_exact += 1;
            if g.det() > 0 {
                m.det_plus_vertex_exact += 1;
            } else {
                m.det_minus_vertex_exact += 1;
            }
            if pure_permutation {
                m.pure_permutation_exact += 1;
            }
            if pure_sign_flip {
                m.pure_sign_flip_exact += 1;
            }
        } else {
            // Only the vertices that actually moved, paired with where they
            // moved from. See `multiset_difference` for why a positional walk
            // over the full sorted lists is the wrong instrument.
            let (only_got, only_want) = multiset_difference(&got, &want);
            m.worst_differing_vertices = m.worst_differing_vertices.max(only_got.len());
            for (a, b) in only_got.iter().zip(only_want.iter()) {
                for k in 0..3 {
                    if a[k] != b[k] {
                        m.worst_component_ulp = m.worst_component_ulp.max(ulp_distance(a[k], b[k]));
                    }
                }
            }
            if m.vertex_failing_mask == 0 {
                m.first_failing_element = g.label();
                m.first_failing_det = g.det();
            }
            m.vertex_failing_mask |= 1u64 << index;
            failing.push(g.short());
        }
        if triangle_ok {
            m.triangle_exact += 1;
        } else {
            m.triangle_failing_mask |= 1u64 << index;
            if vertex_ok {
                m.vertex_exact_triangle_fail += 1;
            }
        }

        // **The control.** Element 0 is the identity, and the identity going
        // through the wrapper must reproduce the reference exactly. A failure
        // is non-determinism in the extractor and voids the row.
        assert!(
            index != 0 || (vertex_ok && triangle_ok),
            "the identity element is not exact: the extractor is not \
             deterministic, so this row measures nothing"
        );
    }

    m.vertex_failing_labels = failing.join("|");
    m.wall_ms = started.elapsed().as_millis();
    m
}

/// `primal` or `dual`, as the registration splits them.
fn family(extractor: &str) -> &'static str {
    match extractor {
        "marching_cubes"
        | "marching_cubes+decider"
        | "marching_tetrahedra"
        | "subgrid_marching_tetrahedra" => "primal",
        "surface_nets" | "dual_contouring" | "manifold_dual_contouring" => "dual",
        other => panic!("{other} is in neither family: the registration names all seven"),
    }
}

/// What the running table's abbreviated columns mean.
const LEGEND: &str = "vex=elements_vertex_exact tex=elements_triangle_exact \
                      raw=elements_vertex_exact_raw \
                      v!t=elements_vertex_exact_triangle_fail \
                      ulp=worst_component_ulp dv=worst_differing_vertices \
                      p6=pure_permutation_exact/6 f8=pure_sign_flip_exact/8 \
                      osEdges=order_sensitive_edges";

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-57");
    common::experiment::run(prereg, |run| {
        let elements = group();
        verify_group(&elements);

        println!("the group, by index — bit i of the failing masks is element i:");
        for (i, e) in elements.iter().enumerate() {
            print!("  {i:>2} {} det={:+}", e.short(), e.det());
            if i % 4 == 3 {
                println!();
            }
        }
        println!("\n{LEGEND}");
        println!(
            "{:<15} {:<28} {:>3} {:>7} {:>4} {:>4} {:>4} {:>4} {:>6} {:>6} {:>3} {:>3} {:>7} {:>7}",
            "field",
            "extractor",
            "n",
            "verts",
            "vex",
            "tex",
            "raw",
            "v!t",
            "ulp",
            "dv",
            "p6",
            "f8",
            "osEdges",
            "ms"
        );

        let mut any_fixture_can_fail = false;
        let mut reference = MeshBuffer::<f64>::new();
        let mut rotated = MeshBuffer::<f64>::new();

        isomesh::for_each_reference_field!(f64, |field_name, field| {
            let (lo, hi) = field.domain();
            let l = hi[0];
            for k in 0..3 {
                assert_eq!(
                    lo[k].to_bits(),
                    (-hi[k]).to_bits(),
                    "{field_name}: domain is not the symmetric cube this \
                     experiment's fixtures assume"
                );
            }

            for fx in fixtures(l) {
                let facts = fixture_facts(&field, &fx);
                assert!(
                    facts.grid_symmetric,
                    "{field_name} at {}³: the grid is not bit-exactly closed \
                     under a sign flip, so the relation would be falsified by \
                     the fixture rather than by the extractor",
                    fx.samples
                );
                if facts.order_sensitive_edges > 0 {
                    any_fixture_can_fail = true;
                }

                isomesh::for_each_extractor!(f64, |extractor_name, extractor| {
                    let m = measure(
                        &field,
                        &mut extractor,
                        &fx,
                        &elements,
                        &mut reference,
                        &mut rotated,
                    );
                    println!(
                        "{:<15} {:<28} {:>3} {:>7} {:>4} {:>4} {:>4} {:>4} {:>6} {:>6} {:>3} {:>3} {:>7} {:>7}",
                        field_name,
                        extractor_name,
                        fx.samples,
                        m.vertices,
                        m.vertex_exact,
                        m.triangle_exact,
                        m.vertex_exact_raw,
                        m.vertex_exact_triangle_fail,
                        m.worst_component_ulp,
                        m.worst_differing_vertices,
                        m.pure_permutation_exact,
                        m.pure_sign_flip_exact,
                        facts.order_sensitive_edges,
                        m.wall_ms
                    );
                    run.record(&[
                        ("field", field_name.to_string()),
                        ("extractor", extractor_name.to_string()),
                        ("family", family(extractor_name).to_string()),
                        ("samples_per_axis", fx.samples.to_string()),
                        ("vertices", m.vertices.to_string()),
                        ("elements_tested", GROUP_ORDER.to_string()),
                        ("elements_vertex_exact", m.vertex_exact.to_string()),
                        ("elements_triangle_exact", m.triangle_exact.to_string()),
                        ("first_failing_element", m.first_failing_element.clone()),
                        ("first_failing_det", format!("{:+}", m.first_failing_det)),
                        ("worst_component_ulp", m.worst_component_ulp.to_string()),
                        (
                            "fixture_can_fail",
                            (facts.order_sensitive_edges > 0).to_string(),
                        ),
                        ("triangles", m.triangles.to_string()),
                        ("grid_symmetric", facts.grid_symmetric.to_string()),
                        ("cell_size", format!("{:.9}", fx.cell_size)),
                        ("origin", format!("{:.9}", fx.origin)),
                        ("cut_edges", facts.cut_edges.to_string()),
                        (
                            "order_sensitive_edges",
                            facts.order_sensitive_edges.to_string(),
                        ),
                        (
                            "boundary_sign_constant",
                            facts.boundary_sign_constant.to_string(),
                        ),
                        ("grid_zero_coordinates", facts.zero_coordinates.to_string()),
                        ("elements_vertex_exact_raw", m.vertex_exact_raw.to_string()),
                        (
                            "elements_triangle_exact_raw",
                            m.triangle_exact_raw.to_string(),
                        ),
                        (
                            "elements_vertex_exact_triangle_fail",
                            m.vertex_exact_triangle_fail.to_string(),
                        ),
                        (
                            "vertex_count_mismatches",
                            m.vertex_count_mismatches.to_string(),
                        ),
                        ("det_plus_vertex_exact", m.det_plus_vertex_exact.to_string()),
                        (
                            "det_minus_vertex_exact",
                            m.det_minus_vertex_exact.to_string(),
                        ),
                        (
                            "pure_permutation_exact",
                            m.pure_permutation_exact.to_string(),
                        ),
                        ("pure_sign_flip_exact", m.pure_sign_flip_exact.to_string()),
                        (
                            "worst_differing_vertices",
                            m.worst_differing_vertices.to_string(),
                        ),
                        (
                            "vertex_failing_mask",
                            format!("{:012x}", m.vertex_failing_mask),
                        ),
                        (
                            "triangle_failing_mask",
                            format!("{:012x}", m.triangle_failing_mask),
                        ),
                        ("vertex_failing_labels", m.vertex_failing_labels.clone()),
                        ("wall_ms", m.wall_ms.to_string()),
                    ]);
                });
            }
        });

        assert!(
            any_fixture_can_fail,
            "VOID: not one fixture is order-sensitive, so every crossing is \
             exactly representable and no pass rate here means anything \
             (M-178's vacuous-fixture trap)"
        );
    });
}
