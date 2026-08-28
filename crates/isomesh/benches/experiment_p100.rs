//! **P-100 — is there a cell decomposition that is octahedrally invariant?**
//!
//! Ticket: R-100. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p100
//! ```
//!
//! Writes `docs/experiments/p-100.csv`.
//!
//! # Hypothesis
//!
//! `M-372` found Marching Tetrahedra stuck at a flat **12 of 48** octahedral
//! elements and named the obstruction structural: a six-tetrahedron
//! decomposition of a cell is not octahedrally invariant, because all six
//! tetrahedra share the main diagonal `0–7` and a relabelling that moves that
//! diagonal moves the whole decomposition. That is a statement about the
//! *decomposition*, and invariant decompositions exist. The one measured here
//! adds the six face centres and the cell centre and splits the cell into
//! **24 tetrahedra** — one per (face, face-edge) pair — a set the full
//! octahedral group permutes onto itself.
//!
//! **C1** the 24-tet decomposition reaches 48 of 48 on all eight fields at both
//! resolutions, where the six-tet split reaches 12. Falsified by under 48.
//!
//! **C2** it costs at most 2× the six-tet split in triangle count, and the
//! combined penalty against Marching Cubes stays under 6×. Falsified by above
//! 2×.
//!
//! **C3** it still tiles across chunk seams with zero open edges. Falsified by
//! any open edge.
//!
//! # The registered vacuity control
//!
//! *"The six-tet arm must reproduce `p-61.csv`'s flat 12 of 48 before the 24-tet
//! arm is believed, reported side by side."* So `p-61.csv` is read, not quoted:
//! [`p61_baseline`] parses the committed file and every `six_tet_crate` row
//! asserts its `elements_vertex_exact`, `cut_edges` and `triangles` against it.
//! The 12 is a number in a file this harness opens.
//!
//! That control has two halves and both matter. An instrument that always said
//! 48 would make C1 meaningless, and an instrument that could never say 48 would
//! make it unreachable — so the run brackets itself: `marching_cubes` must read
//! **48** and `six_tet_crate` must read **12**, in the same instrument, on the
//! same rows, before any 24-tet number is looked at.
//!
//! # Five arms, one instrument
//!
//! | `decomposition` | `tetrahedra_per_cell` | what it is for |
//! |---|---|---|
//! | `marching_cubes` | 0 | the triangle denominator, and the instrument's 48 |
//! | `six_tet_crate` | 6 | `p-61.csv`'s row, reproduced — the vacuity control |
//! | `six_tet_bench` | 6 | the *same* generic marcher as the 24-tet arm, on Kuhn's six |
//! | `barycentric_24` | 24 | **the registered arm** |
//! | `barycentric_24_field_sampled` | 24 | the same split with the field sampled at the new points |
//!
//! `six_tet_bench` is the arm that makes the comparison mean something. It runs
//! the identical runtime tetrahedron contour, the identical vertex cache and the
//! identical crossing placement as `barycentric_24`, and differs **only** in the
//! tetrahedron list. If it reads 12 while `barycentric_24` reads 48, the
//! difference is the decomposition and cannot be the marcher.
//!
//! # The instrument is P-61's, checked rather than trusted
//!
//! The octahedral group, [`Element`], [`vertex_keys`], [`multiset_difference`],
//! [`ulp_distance`], the two fixtures and the identity-element control are
//! `experiment_p61.rs`'s, unchanged. A bench cannot import another bench, so the
//! copy is **held against the committed artefact** the way P-61 held itself
//! against `p-57.csv`: `cut_edges` and `order_sensitive_edges` are properties of
//! the grid and the field values alone, no extractor touches them, and a
//! mismatch means this is not P-61's instrument. `p61_columns_match` records the
//! outcome per row.
//!
//! # Where the 24 tetrahedra get their values
//!
//! A mesher over a sampled grid has eight numbers per cell, and the face centres
//! and the cell centre are not among them. `barycentric_24` therefore takes the
//! **multilinear** value there — the mean of the 4 or 8 surrounding corners —
//! which keeps the arm on exactly the same data as the six-tet arm and makes the
//! triangle ratio a property of the decomposition rather than of a finer sample
//! set.
//!
//! The mean is summed **in sorted order**, and that is not a detail: floating
//! addition is commutative but not associative, so a fixed corner order would
//! make the face-centre value depend on the axis labelling and C1 would be
//! falsified by the averaging rather than by the decomposition. Sorting by
//! [`f64::total_cmp`] makes the sum a function of the *multiset* of corner
//! values, which an octahedral relabelling preserves exactly.
//!
//! `barycentric_24_field_sampled` is the other choice — ask the field — and is
//! carried as an extra arm to price the averaging.
//!
//! # C3's own control
//!
//! `open_edges` is a difference, not a census: the two-chunk mesh and the
//! single-pass mesh over the same union grid are both welded on exact position
//! bits, and an open edge is one that is unpaired in the two-chunk mesh and
//! **not** unpaired in the single pass. That subtraction is what removes the
//! domain boundary — `gyroid` and `fbm_terrain` are open surfaces and their
//! boundary is not a seam crack.
//!
//! A difference that is always zero proves nothing, so the run also meshes a
//! deliberately broken seam: `seam_control_open_edges` gives the second chunk
//! Kuhn's six tetrahedra anchored on corner `1` instead of corner `0`, which is
//! `P-3`'s crack — the two chunks then split their shared face on different
//! diagonals. That column must be non-zero or the counter is not measuring
//! anything, and `seam_configurations` counts the distinct cut sign patterns the
//! seam plane actually carried.
//!
//! # Falsified by
//!
//! C1 by any `barycentric_24` row under 48. C2 by `triangles_vs_six_tet` above
//! 2 on any row, or `triangles_vs_marching_cubes` above 6. C3 by any non-zero
//! `open_edges`.

// Exact comparisons throughout: this harness is about bit patterns.
#![allow(clippy::float_cmp)]

mod common;

use std::collections::HashMap;
use std::path::PathBuf;
use std::time::Instant;

use isomesh::extractor::Extractor;
use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_tetrahedra::MarchingTetrahedra;
use isomesh::validate::{AccuracyConfig, accuracy};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

// ─── the crossing, in the one frame the crate uses ───────────────────────────

/// `cube::edge_offset`, spelled out because `cube` is private.
///
/// Identical expression to the shipped one — `((a + b)/2)/(a − b)`, R-059's
/// centred offset. A bench comparing decompositions has to place crossings the
/// way the crate does or it is measuring a placement change instead.
#[inline]
fn edge_offset(a: f64, b: f64) -> f64 {
    ((a + b) * 0.5) / (a - b)
}

/// `cube::place`, likewise. `(lo + hi)/2 + (hi − lo)·d`.
///
/// P-61's finding is load-bearing here: this form is invariant under swapping
/// the endpoints, because `fl(b − a) = −fl(a − b)` exactly and `d` negates
/// exactly with them. So the marcher below never has to canonicalise which end
/// of a tetrahedron edge is "first", which is what lets a relabelled cell
/// produce bit-identical vertices.
#[inline]
fn place(lo: f64, hi: f64, d: f64) -> f64 {
    (lo + hi) * 0.5 + (hi - lo) * d
}

/// The crate's sign convention: negative inside, zero counts as outside.
#[inline]
fn is_inside(v: f64) -> bool {
    v < 0.0
}

// ─── the group: P-61's, unchanged ───────────────────────────────────────────

/// The 6 axis permutations. Crossed with 8 sign patterns this is all 48.
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

/// Probe points for the inverse round-trip check. Both zeros deliberately.
const PROBES: [[f64; 3]; 4] = [
    [0.3, -1.7, 2.9],
    [0.0, -0.0, 1.5],
    [-2.25, 0.656_25, -0.187_5],
    [1.0, 1.0, 1.0],
];

/// One element of the octahedral group, as a signed axis permutation.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Element {
    perm: [usize; 3],
    sign: [i8; 3],
}

impl Element {
    #[inline]
    fn apply(self, p: [f64; 3]) -> [f64; 3] {
        std::array::from_fn(|k| {
            let v = p[self.perm[k]];
            if self.sign[k] < 0 { -v } else { v }
        })
    }

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

    fn det(self) -> i32 {
        let mut m = [[0i32; 3]; 3];
        for k in 0..3 {
            m[k][self.perm[k]] = i32::from(self.sign[k]);
        }
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

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

    /// Does this element fix the cell's main diagonal `0–7`?
    ///
    /// The diagonal runs along `(1,1,1)`, so a signed permutation maps it to
    /// itself exactly when every sign agrees. Six permutations times the two
    /// uniform sign patterns is **twelve**, and that is the mechanism claim C1
    /// is testing: the six-tet split's 12 should be precisely this subgroup.
    fn fixes_main_diagonal(self) -> bool {
        self.sign == [1, 1, 1] || self.sign == [-1, -1, -1]
    }
}

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
    let (mut rotations, mut reflections) = (0, 0);
    let mut diagonal_fixing = 0;
    for e in g {
        match e.det() {
            1 => rotations += 1,
            -1 => reflections += 1,
            d => panic!("{} has determinant {d}, not +/-1", e.label()),
        }
        if e.fixes_main_diagonal() {
            diagonal_fixing += 1;
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
    assert_eq!(
        diagonal_fixing, 12,
        "the stabiliser of the main diagonal must have order 12: that is the \
         number the six-tet arm is predicted to reach"
    );
}

/// `g·f`, i.e. `(g·f)(p) = f(g⁻¹·p)`.
struct Rotated<'a, S> {
    field: &'a S,
    g: Element,
    g_inv: Element,
}

impl<S: Sdf<Scalar = f64>> Sdf for Rotated<'_, S> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        self.field.sample(self.g_inv.apply(p))
    }

    #[inline]
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        self.g.apply(self.field.gradient(self.g_inv.apply(p)))
    }
}

// ─── bit keys: P-61's, unchanged ────────────────────────────────────────────

const NEGATIVE_ZERO: u64 = 1u64 << 63;

#[inline]
fn key(v: f64, normalise: bool) -> u64 {
    let b = v.to_bits();
    if normalise && b == NEGATIVE_ZERO { 0 } else { b }
}

#[inline]
fn monotone(bits: u64) -> i128 {
    if bits & NEGATIVE_ZERO == 0 {
        i128::from(bits)
    } else {
        -i128::from(bits & !NEGATIVE_ZERO)
    }
}

#[inline]
fn ulp_distance(a: u64, b: u64) -> u128 {
    (monotone(a) - monotone(b)).unsigned_abs()
}

/// Multiset symmetric difference of two **sorted** key lists.
fn multiset_difference(a: &[[u64; 3]], b: &[[u64; 3]]) -> (Vec<[u64; 3]>, Vec<[u64; 3]>) {
    let (mut only_a, mut only_b) = (Vec::new(), Vec::new());
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

/// Triangles as **unordered** vertex-key triples, sorted.
///
/// Unordered on purpose: half the group is a reflection, which reverses
/// winding, so a winding-sensitive key would report 24 failures on an extractor
/// that is perfectly equivariant. C1 is a statement about the vertex set;
/// this column is the corroborating statement about which triples are joined.
fn triangle_keys(
    positions: &[[f64; 3]],
    indices: &[u32],
    g: Option<Element>,
    normalise: bool,
) -> Vec<[[u64; 3]; 3]> {
    let mut out = Vec::with_capacity(indices.len() / 3);
    for tri in indices.as_chunks::<3>().0 {
        let mut t: [[u64; 3]; 3] = std::array::from_fn(|k| {
            let p = positions[tri[k] as usize];
            let q = match g {
                Some(e) => e.apply(p),
                None => p,
            };
            std::array::from_fn(|c| key(q[c], normalise))
        });
        t.sort_unstable();
        out.push(t);
    }
    out.sort_unstable();
    out
}

// ─── the decompositions ─────────────────────────────────────────────────────

/// A point of a decomposed cell.
///
/// Eight corners, six face centres, one cell centre — the most any arm here
/// needs, and the six-tet arms simply never touch slots 8..15.
const CELL_POINTS: usize = 15;

/// Slot of the face centre normal to `axis` on side `side` (`0` = low).
#[inline]
const fn face_slot(axis: usize, side: usize) -> usize {
    8 + 2 * axis + side
}

/// Slot of the cell centre.
const CENTRE_SLOT: usize = 14;

/// The two corners of each edge of a tetrahedron, as indices into its corner
/// list. Lower index first; the crate's `TET_EDGES`.
const TET_EDGES: [[usize; 2]; 6] = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];

/// Kuhn's six tetrahedra, anchored on cube corner `anchor`.
///
/// The six monotone paths from `anchor` to its antipode along cube edges, one
/// per ordering of the three axes — the crate's `build_tets`, generalised in the
/// one way this harness needs: `anchor = 0` is the shipped decomposition, and a
/// different anchor is the same tiling on a **different main diagonal**, which
/// is `P-3`'s crack and the positive control for the seam counter.
///
/// # Not every other anchor is a crack, and the control caught it
///
/// The first version of this control used `anchor = 1` and reported **zero**
/// open edges — a clean pass on a fixture that cannot fail, which is `M-44`
/// exactly. The reason is worth writing down: Kuhn's diagonal on the face
/// normal to axis `a` joins the face corner where the *other two* axes both
/// carry the anchor's bits to the corner where both are flipped, so it depends
/// only on the anchor's bits on those two axes. `anchor = 1` differs from
/// `anchor = 0` on axis `x` alone, so both split every `x`-face on `0–6`/`1–7`
/// and they tile across an `x`-seam perfectly — they differ only *inside* the
/// cell. `anchor = 7` is the same story on all six faces: face-compatible with
/// `anchor = 0` everywhere.
///
/// An `x`-seam therefore needs an anchor whose `y` or `z` bit differs:
/// `anchor = 2` splits every `x`-face on `2–4`/`3–5` instead, and that is the
/// mismatch `P-3` describes.
fn six_tets(anchor: u8) -> Vec<[u8; 4]> {
    let mut out = Vec::with_capacity(6);
    for order in PERMS {
        let mut tet = [anchor; 4];
        let mut c = anchor;
        for (step, axis) in order.iter().enumerate() {
            c ^= 1 << axis;
            tet[step + 1] = c;
        }
        out.push(tet);
    }
    out
}

/// The 24-tetrahedron barycentric split: cell centre, face centre, and one edge
/// of that face.
///
/// Each of the six faces is cut into four triangles by its centre, and each
/// triangle is coned to the cell centre. The construction never mentions an
/// axis order or a diagonal, so the octahedral group permutes the 24 onto
/// themselves: an element sends a face to a face and a face edge to a face
/// edge, and the stabiliser of one tetrahedron is the order-2 reflection that
/// swaps its two cube corners — 48/2 = 24, one orbit.
fn barycentric_24_tets() -> Vec<[u8; 4]> {
    let mut out = Vec::with_capacity(24);
    for axis in 0..3usize {
        let (b, c) = ((axis + 1) % 3, (axis + 2) % 3);
        for side in 0..2usize {
            // The face's four corners, walked round the ring so consecutive
            // pairs are face edges.
            let ring = [[0, 0], [1, 0], [1, 1], [0, 1]];
            let corner_of = |uv: [usize; 2]| -> u8 {
                let mut bits = 0u8;
                if side == 1 {
                    bits |= 1 << axis;
                }
                if uv[0] == 1 {
                    bits |= 1 << b;
                }
                if uv[1] == 1 {
                    bits |= 1 << c;
                }
                bits
            };
            for k in 0..4usize {
                out.push([
                    CENTRE_SLOT as u8,
                    face_slot(axis, side) as u8,
                    corner_of(ring[k]),
                    corner_of(ring[(k + 1) % 4]),
                ]);
            }
        }
    }
    out
}

/// Which decomposition an arm marches, and where the extra points get values.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Decomp {
    /// Kuhn's six, on the main diagonal through `anchor`.
    SixTet { anchor: u8 },
    /// The barycentric 24. `sampled` asks the field at the face and cell
    /// centres instead of averaging the surrounding corners.
    Barycentric24 { sampled: bool },
}

impl Decomp {
    fn tets(self) -> Vec<[u8; 4]> {
        match self {
            Self::SixTet { anchor } => six_tets(anchor),
            Self::Barycentric24 { .. } => barycentric_24_tets(),
        }
    }

    fn uses_extra_points(self) -> bool {
        matches!(self, Self::Barycentric24 { .. })
    }
}

/// The mean of a corner multiset, summed in **sorted** order.
///
/// Sorted so the result is a function of the multiset and not of the labelling:
/// floating addition is commutative but not associative, and a fixed corner
/// order would leak the axis labels into the face-centre value. `total_cmp` is a
/// total order on the bits, so an octahedral relabelling — which permutes the
/// corners and leaves every value bit-identical — cannot move this sum.
///
/// The scale is a power of two, so the division is exact.
#[inline]
fn sorted_mean<const N: usize>(values: [f64; N]) -> f64 {
    let mut v = values;
    v.sort_unstable_by(f64::total_cmp);
    let mut sum = 0.0f64;
    for x in v {
        sum += x;
    }
    sum * (1.0 / N as f64)
}

// ─── the marcher ────────────────────────────────────────────────────────────

/// Positions and indices. No normals: nothing downstream of here reads them.
struct Mesh {
    positions: Vec<[f64; 3]>,
    indices: Vec<u32>,
}

impl Mesh {
    fn triangles(&self) -> usize {
        self.indices.len() / 3
    }
}

/// A tetrahedral marcher whose decomposition is a parameter.
///
/// One code path for every tetrahedral arm — the case classification, the quad
/// ordering, the crossing placement and the vertex cache are shared, so a
/// difference between two arms is a difference between two tetrahedron lists
/// and can be nothing else. That is what makes `six_tet_bench`'s 12 evidence
/// about `barycentric_24`'s 48.
struct TetMarcher {
    tets: Vec<[u8; 4]>,
    decomp: Decomp,
    values: Vec<f64>,
    cache: HashMap<u64, u32>,
}

impl TetMarcher {
    fn new(decomp: Decomp) -> Self {
        Self {
            tets: decomp.tets(),
            decomp,
            values: Vec::new(),
            cache: HashMap::new(),
        }
    }

    fn mesh<S: Sdf<Scalar = f64>>(
        &mut self,
        field: &S,
        samples: [u32; 3],
        origin: [f64; 3],
        h: f64,
    ) -> Mesh {
        let n: [usize; 3] = std::array::from_fn(|k| samples[k] as usize);
        let coords: [Vec<f64>; 3] =
            std::array::from_fn(|k| (0..n[k]).map(|i| origin[k] + h * i as f64).collect());

        self.values.clear();
        self.values.reserve(n[0] * n[1] * n[2]);
        for z in 0..n[2] {
            for y in 0..n[1] {
                for x in 0..n[0] {
                    self.values
                        .push(field.sample([coords[0][x], coords[1][y], coords[2][z]]));
                }
            }
        }
        let stride = |x: usize, y: usize, z: usize| x + n[0] * (y + n[1] * z);

        self.cache.clear();
        let mut mesh = Mesh {
            positions: Vec::new(),
            indices: Vec::new(),
        };

        let sample_total = n[0] * n[1] * n[2];
        let extra = self.decomp.uses_extra_points();
        let mut pos = [[0.0f64; 3]; CELL_POINTS];
        let mut val = [0.0f64; CELL_POINTS];
        let mut pid = [0u64; CELL_POINTS];

        for z in 0..n[2] - 1 {
            for y in 0..n[1] - 1 {
                for x in 0..n[0] - 1 {
                    let base = [x, y, z];
                    let mut all_same = true;
                    let first = is_inside(self.values[stride(x, y, z)]);
                    for c in 0..8usize {
                        let off = [c & 1, (c >> 1) & 1, (c >> 2) & 1];
                        let s = stride(x + off[0], y + off[1], z + off[2]);
                        val[c] = self.values[s];
                        pos[c] = [
                            coords[0][x + off[0]],
                            coords[1][y + off[1]],
                            coords[2][z + off[2]],
                        ];
                        // Corner point identity: the global sample it sits on.
                        pid[c] = (s as u64) * 4;
                        if is_inside(val[c]) != first {
                            all_same = false;
                        }
                    }

                    if extra {
                        for axis in 0..3usize {
                            let (b, c) = ((axis + 1) % 3, (axis + 2) % 3);
                            for side in 0..2usize {
                                let slot = face_slot(axis, side);
                                let mut p = [0.0f64; 3];
                                p[axis] = coords[axis][base[axis] + side];
                                p[b] = (coords[b][base[b]] + coords[b][base[b] + 1]) * 0.5;
                                p[c] = (coords[c][base[c]] + coords[c][base[c] + 1]) * 0.5;
                                pos[slot] = p;
                                let corners: [f64; 4] = std::array::from_fn(|k| {
                                    let mut off = [0usize; 3];
                                    off[axis] = side;
                                    off[b] = k & 1;
                                    off[c] = (k >> 1) & 1;
                                    val[off[0] | (off[1] << 1) | (off[2] << 2)]
                                });
                                val[slot] = if let Decomp::Barycentric24 { sampled: true } =
                                    self.decomp
                                {
                                    field.sample(p)
                                } else {
                                    sorted_mean(corners)
                                };
                                // Face identity: the sample at the low corner of
                                // the face's own plane, so both cells across a
                                // face agree on the key.
                                let mut fb = base;
                                fb[axis] += side;
                                pid[slot] = (stride(fb[0], fb[1], fb[2]) as u64 * 3 + axis as u64)
                                    * 4
                                    + 1;
                            }
                        }
                        let centre: [f64; 3] = std::array::from_fn(|k| {
                            (coords[k][base[k]] + coords[k][base[k] + 1]) * 0.5
                        });
                        pos[CENTRE_SLOT] = centre;
                        val[CENTRE_SLOT] =
                            if let Decomp::Barycentric24 { sampled: true } = self.decomp {
                                field.sample(centre)
                            } else {
                                sorted_mean([
                                    val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7],
                                ])
                            };
                        pid[CENTRE_SLOT] = (stride(x, y, z) as u64) * 4 + 2;
                    }

                    // The averaged arm cannot find a crossing in a cell whose
                    // eight corners agree, because a mean of like-signed values
                    // keeps the sign. The sampled arm can, so it never skips.
                    let sampled = matches!(self.decomp, Decomp::Barycentric24 { sampled: true });
                    if all_same && !sampled {
                        continue;
                    }

                    for t in 0..self.tets.len() {
                        self.contour_tet(t, &pos, &val, &pid, sample_total, &mut mesh);
                    }
                }
            }
        }

        mesh
    }

    /// Contour one tetrahedron. Sixteen sign cases, no ambiguity.
    fn contour_tet(
        &mut self,
        t: usize,
        pos: &[[f64; 3]; CELL_POINTS],
        val: &[f64; CELL_POINTS],
        pid: &[u64; CELL_POINTS],
        sample_total: usize,
        mesh: &mut Mesh,
    ) {
        let tet: [usize; 4] = std::array::from_fn(|k| self.tets[t][k] as usize);
        let mut mask = 0u8;
        for i in 0..4 {
            if is_inside(val[tet[i]]) {
                mask |= 1 << i;
            }
        }
        if mask == 0 || mask == 15 {
            return;
        }

        let mut cut = [0usize; 4];
        let mut cuts = 0usize;
        for (e, [i, j]) in TET_EDGES.iter().enumerate() {
            if (mask >> i) & 1 != (mask >> j) & 1 {
                cut[cuts] = e;
                cuts += 1;
            }
        }
        assert!(
            cuts == 3 || cuts == 4,
            "a tetrahedron with a sign change cuts three or four edges, not {cuts}"
        );

        // Inside centroid, for winding. Cheap and exact enough: the test is a
        // sign, and a triangle whose plane contains the inside centroid has
        // zero area anyway.
        let mut inside_c = [0.0f64; 3];
        let mut inside_n = 0.0f64;
        for i in 0..4 {
            if (mask >> i) & 1 == 1 {
                for k in 0..3 {
                    inside_c[k] += pos[tet[i]][k];
                }
                inside_n += 1.0;
            }
        }
        for c in &mut inside_c {
            *c /= inside_n;
        }

        let mut ring = [0usize; 4];
        let count = if cuts == 3 {
            ring[..3].copy_from_slice(&cut[..3]);
            3
        } else {
            // Two cut edges are adjacent on the quad when they share a corner
            // of the tetrahedron; the pair sharing none is the diagonal.
            ring[0] = cut[0];
            let mut placed = [false; 4];
            placed[0] = true;
            let mut adj = Vec::with_capacity(2);
            let mut diag = usize::MAX;
            for (i, &e) in cut.iter().enumerate().skip(1) {
                if edges_share_a_corner(cut[0], e) {
                    adj.push(e);
                } else {
                    diag = e;
                    placed[i] = true;
                }
            }
            assert_eq!(adj.len(), 2, "a quad's ring has exactly two neighbours");
            assert_ne!(diag, usize::MAX, "a quad has exactly one diagonal");
            ring[1] = adj[0];
            ring[2] = diag;
            ring[3] = adj[1];
            4
        };

        let mut ids = [0u32; 4];
        for (k, &e) in ring[..count].iter().enumerate() {
            let [a, b] = TET_EDGES[e];
            ids[k] = self.vertex_on(tet[a], tet[b], pos, val, pid, sample_total, mesh);
        }

        let tri_of = |mesh: &mut Mesh, a: u32, b: u32, c: u32| {
            let p = [
                mesh.positions[a as usize],
                mesh.positions[b as usize],
                mesh.positions[c as usize],
            ];
            let u: [f64; 3] = std::array::from_fn(|k| p[1][k] - p[0][k]);
            let v: [f64; 3] = std::array::from_fn(|k| p[2][k] - p[0][k]);
            let cr = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let to_inside: [f64; 3] = std::array::from_fn(|k| inside_c[k] - p[0][k]);
            let dot = cr[0] * to_inside[0] + cr[1] * to_inside[1] + cr[2] * to_inside[2];
            // Normals point away from the solid, as everywhere else in the crate.
            if dot > 0.0 {
                mesh.indices.extend_from_slice(&[a, c, b]);
            } else {
                mesh.indices.extend_from_slice(&[a, b, c]);
            }
        };

        if count == 3 {
            tri_of(mesh, ids[0], ids[1], ids[2]);
        } else {
            tri_of(mesh, ids[0], ids[1], ids[2]);
            tri_of(mesh, ids[0], ids[2], ids[3]);
        }
    }

    /// The vertex on one cut edge, created once per edge of the lattice.
    ///
    /// Keyed on the two **point identities** — global sample for a corner,
    /// (plane sample, axis) for a face centre, cell for the centre — so every
    /// cell that contains an edge agrees about its vertex and the traversal
    /// order cannot show up in the output.
    #[allow(clippy::too_many_arguments)]
    fn vertex_on(
        &mut self,
        a: usize,
        b: usize,
        pos: &[[f64; 3]; CELL_POINTS],
        val: &[f64; CELL_POINTS],
        pid: &[u64; CELL_POINTS],
        sample_total: usize,
        mesh: &mut Mesh,
    ) -> u32 {
        let span = (sample_total as u64) * 12 + 4;
        debug_assert!(pid[a] < span && pid[b] < span, "point identity overflows");
        let (lo_id, hi_id) = if pid[a] < pid[b] {
            (pid[a], pid[b])
        } else {
            (pid[b], pid[a])
        };
        let ekey = lo_id * span + hi_id;
        if let Some(&i) = self.cache.get(&ekey) {
            return i;
        }
        let d = edge_offset(val[a], val[b]);
        let p: [f64; 3] = std::array::from_fn(|k| place(pos[a][k], pos[b][k], d));
        let index = u32::try_from(mesh.positions.len()).expect("vertex count fits u32");
        mesh.positions.push(p);
        self.cache.insert(ekey, index);
        index
    }
}

/// Do these two tetrahedron edges meet at a corner?
fn edges_share_a_corner(x: usize, y: usize) -> bool {
    let a = TET_EDGES[x];
    let b = TET_EDGES[y];
    a[0] == b[0] || a[0] == b[1] || a[1] == b[0] || a[1] == b[1]
}

// ─── the arms ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Arm {
    MarchingCubes,
    SixTetCrate,
    SixTetBench,
    Barycentric24,
    Barycentric24Sampled,
}

impl Arm {
    fn name(self) -> &'static str {
        match self {
            Self::MarchingCubes => "marching_cubes",
            Self::SixTetCrate => "six_tet_crate",
            Self::SixTetBench => "six_tet_bench",
            Self::Barycentric24 => "barycentric_24",
            Self::Barycentric24Sampled => "barycentric_24_field_sampled",
        }
    }

    fn tets_per_cell(self) -> usize {
        match self {
            Self::MarchingCubes => 0,
            Self::SixTetCrate | Self::SixTetBench => 6,
            Self::Barycentric24 | Self::Barycentric24Sampled => 24,
        }
    }

    fn decomp(self) -> Option<Decomp> {
        match self {
            Self::MarchingCubes | Self::SixTetCrate => None,
            Self::SixTetBench => Some(Decomp::SixTet { anchor: 0 }),
            Self::Barycentric24 => Some(Decomp::Barycentric24 { sampled: false }),
            Self::Barycentric24Sampled => Some(Decomp::Barycentric24 { sampled: true }),
        }
    }
}

const ARMS: [Arm; 5] = [
    Arm::MarchingCubes,
    Arm::SixTetCrate,
    Arm::SixTetBench,
    Arm::Barycentric24,
    Arm::Barycentric24Sampled,
];

fn mesh_of<S: Sdf<Scalar = f64>>(
    arm: Arm,
    field: &S,
    samples: [u32; 3],
    origin: [f64; 3],
    h: f64,
) -> Mesh {
    let shape = RuntimeShape3::new(samples).expect("grid fits u32");
    match arm {
        Arm::MarchingCubes => {
            let mut mc = MarchingCubes::<f64>::new();
            let mut out = MeshBuffer::<f64>::new();
            mc.extract_into(field, &shape, origin, h, &mut out)
                .expect("marching cubes extraction");
            Mesh {
                positions: out.positions,
                indices: out.indices,
            }
        }
        Arm::SixTetCrate => {
            let mut mt = MarchingTetrahedra::<f64>::new();
            let mut out = MeshBuffer::<f64>::new();
            mt.extract_into(field, &shape, origin, h, &mut out)
                .expect("marching tetrahedra extraction");
            Mesh {
                positions: out.positions,
                indices: out.indices,
            }
        }
        other => {
            let decomp = other.decomp().expect("a tetrahedral arm has a decomposition");
            TetMarcher::new(decomp).mesh(field, samples, origin, h)
        }
    }
}

// ─── fixtures: P-61's two, unchanged ────────────────────────────────────────

#[derive(Clone, Copy)]
struct Fixture {
    samples: u32,
    cell_size: f64,
    origin: f64,
}

/// P-57's and P-61's two fixtures. The 25³ arm uses `3L/32` rather than `2L/24`,
/// because `L/12` is not dyadic and its grid does not mirror.
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

fn axis_coords(fx: &Fixture) -> Vec<f64> {
    (0..fx.samples)
        .map(|i| fx.origin + fx.cell_size * f64::from(i))
        .collect()
}

/// Grid properties P-61 also reports, recomputed so they can be asserted
/// against the committed file.
struct FixtureFacts {
    grid_symmetric: bool,
    cut_edges: usize,
    order_sensitive_edges: usize,
}

fn fixture_facts<S: Sdf<Scalar = f64>>(field: &S, fx: &Fixture) -> FixtureFacts {
    let n = fx.samples as usize;
    let coords = axis_coords(fx);

    let mut grid_symmetric = true;
    for i in 0..n {
        if key(coords[i], true) != key(-coords[n - 1 - i], true) {
            grid_symmetric = false;
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
                    if is_inside(va) == is_inside(vb) {
                        continue;
                    }
                    cut_edges += 1;
                    let lo_c = coords[[x, y, z][axis]];
                    let hi_c = coords[hi[axis]];
                    // P-61's census, still in the lower-corner frame: it is what
                    // `fixture_can_fail` is derived from and redefining it would
                    // silently redefine the population.
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
        cut_edges,
        order_sensitive_edges,
    }
}

// ─── P-61's committed artefact, as the registered vacuity control ───────────

/// One `equivariance` row of `docs/experiments/p-61.csv`.
struct Baseline {
    elements_vertex_exact: usize,
    cut_edges: usize,
    order_sensitive_edges: usize,
    triangles: usize,
}

/// Read `p-61.csv`. **Not optional**: the registered vacuity control is that the
/// six-tet arm reproduces this file's flat 12 of 48, so the 12 has to come out
/// of the file rather than out of a comment.
fn p61_baseline() -> HashMap<String, Baseline> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/experiments/p-61.csv");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("P-100 needs {} as its control: {e}", path.display()));
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty());
    let header: Vec<&str> = lines
        .next()
        .expect("p-61.csv has a header")
        .split(',')
        .collect();
    let col = |name: &str| {
        header
            .iter()
            .position(|h| *h == name)
            .unwrap_or_else(|| panic!("p-61.csv has no `{name}` column"))
    };
    let (c_block, c_field, c_extractor, c_samples) = (
        col("block"),
        col("field"),
        col("extractor"),
        col("samples_per_axis"),
    );
    let c_vex = col("elements_vertex_exact");
    let c_cut = col("cut_edges");
    let c_ose = col("order_sensitive_edges");
    let c_tris = col("triangles");

    let mut out = HashMap::new();
    for line in lines {
        let f: Vec<&str> = line.split(',').collect();
        if f[c_block] != "equivariance" {
            continue;
        }
        let k = format!("{}/{}/{}", f[c_field], f[c_extractor], f[c_samples]);
        out.insert(
            k,
            Baseline {
                elements_vertex_exact: f[c_vex].parse().expect("integer"),
                cut_edges: f[c_cut].parse().expect("integer"),
                order_sensitive_edges: f[c_ose].parse().expect("integer"),
                triangles: f[c_tris].parse().expect("integer"),
            },
        );
    }
    assert_eq!(
        out.len(),
        112,
        "p-61.csv's equivariance block should carry 8 fields x 2 grids x 7 extractors"
    );
    out
}

// ─── the equivariance measurement ───────────────────────────────────────────

struct Measured {
    vertices: usize,
    triangles: usize,
    vertex_exact: usize,
    triangle_exact: usize,
    worst_component_ulp: u128,
    first_failing_element: String,
    failing_labels: String,
    /// Are the exact elements **precisely** the stabiliser of the main
    /// diagonal? The mechanism claim behind the six-tet arm's 12.
    exact_are_diagonal_stabiliser: bool,
    wall_ms: u128,
}

fn measure<S: Sdf<Scalar = f64>>(
    arm: Arm,
    field: &S,
    fx: &Fixture,
    elements: &[Element],
) -> Measured {
    let samples = [fx.samples; 3];
    let origin = [fx.origin; 3];
    let started = Instant::now();
    let reference = mesh_of(arm, field, samples, origin, fx.cell_size);

    let mut m = Measured {
        vertices: reference.positions.len(),
        triangles: reference.triangles(),
        vertex_exact: 0,
        triangle_exact: 0,
        worst_component_ulp: 0,
        first_failing_element: String::from("none"),
        failing_labels: String::new(),
        exact_are_diagonal_stabiliser: true,
        wall_ms: 0,
    };
    let mut failing: Vec<String> = Vec::new();

    for (index, &g) in elements.iter().enumerate() {
        let wrapped = Rotated {
            field,
            g,
            g_inv: g.inverse(),
        };
        let rotated = mesh_of(arm, &wrapped, samples, origin, fx.cell_size);

        let got = vertex_keys(&rotated.positions, None, true);
        let want = vertex_keys(&reference.positions, Some(g), true);
        let vertex_ok = got == want;

        let got_tri = triangle_keys(&rotated.positions, &rotated.indices, None, true);
        let want_tri = triangle_keys(&reference.positions, &reference.indices, Some(g), true);
        let triangle_ok = got_tri == want_tri;

        if vertex_ok {
            m.vertex_exact += 1;
        } else {
            let (only_got, only_want) = multiset_difference(&got, &want);
            for (a, b) in only_got.iter().zip(only_want.iter()) {
                for k in 0..3 {
                    if a[k] != b[k] {
                        m.worst_component_ulp = m.worst_component_ulp.max(ulp_distance(a[k], b[k]));
                    }
                }
            }
            if failing.is_empty() {
                m.first_failing_element = g.label();
            }
            failing.push(g.short());
        }
        if vertex_ok != g.fixes_main_diagonal() {
            m.exact_are_diagonal_stabiliser = false;
        }
        if triangle_ok {
            m.triangle_exact += 1;
        }

        // **The control**: element 0 is the identity and must reproduce the
        // reference exactly, or the arm is not deterministic and the row
        // measures nothing.
        assert!(
            index != 0 || (vertex_ok && triangle_ok),
            "{}: the identity element is not exact, so this row measures nothing",
            arm.name()
        );
    }

    m.failing_labels = failing.join("|");
    m.wall_ms = started.elapsed().as_millis();
    m
}

// ─── C3: the two-chunk seam ─────────────────────────────────────────────────

/// Undirected edges of a mesh, welded on exact position bits, with their use
/// counts. Winding is irrelevant to whether a mesh is closed, so the key is the
/// unordered pair.
type EdgeCensus = HashMap<([u64; 3], [u64; 3]), usize>;

fn edge_census(mesh: &Mesh) -> EdgeCensus {
    let mut out: EdgeCensus = HashMap::new();
    let welded: Vec<[u64; 3]> = mesh
        .positions
        .iter()
        .map(|p| std::array::from_fn(|k| key(p[k], true)))
        .collect();
    for tri in mesh.indices.as_chunks::<3>().0 {
        for k in 0..3 {
            let a = welded[tri[k] as usize];
            let b = welded[tri[(k + 1) % 3] as usize];
            let e = if a <= b { (a, b) } else { (b, a) };
            *out.entry(e).or_insert(0) += 1;
        }
    }
    out
}

fn unpaired(census: &EdgeCensus) -> Vec<([u64; 3], [u64; 3])> {
    census
        .iter()
        .filter(|&(_, &c)| c != 2)
        .map(|(k, _)| *k)
        .collect()
}

struct Seam {
    /// Unpaired edges the split introduced that the single pass does not have.
    open_edges: usize,
    /// The single pass's own unpaired edges: the domain boundary. Non-zero on
    /// the open fields, which is what proves the counter counts.
    boundary_open_edges: usize,
    /// Is the two-chunk mesh the *same mesh* as the single pass?
    identical: bool,
    /// Distinct cut sign patterns on the shared plane's faces.
    configurations: usize,
    /// Were the two chunks' shared coordinates bit-identical? M-32's question,
    /// asked so a mismatch is attributed to arithmetic and not to the split.
    coords_exact: bool,
    /// The positive control: chunk B marched on a **different main diagonal**,
    /// which is P-3's crack. Must be non-zero.
    control_open_edges: usize,
}

fn seam_census<S: Sdf<Scalar = f64>>(arm: Arm, field: &S, fx: &Fixture) -> Seam {
    let n = fx.samples;
    let h = fx.cell_size;
    let coords = axis_coords(fx);
    let split = (n - 1) / 2;
    let origin = [fx.origin; 3];
    let seam_x = coords[split as usize];

    let mut coords_exact = true;
    for i in 0..(n - split) {
        let theirs = seam_x + h * f64::from(i);
        if theirs.to_bits() != coords[(split + i) as usize].to_bits() {
            coords_exact = false;
        }
    }

    let whole = mesh_of(arm, field, [n; 3], origin, h);
    let a = mesh_of(arm, field, [split + 1, n, n], origin, h);
    let b = mesh_of(arm, field, [n - split, n, n], [seam_x, origin[1], origin[2]], h);

    let two_chunk = Mesh {
        positions: a
            .positions
            .iter()
            .chain(b.positions.iter())
            .copied()
            .collect(),
        indices: a
            .indices
            .iter()
            .copied()
            .chain(
                b.indices
                    .iter()
                    .map(|i| i + u32::try_from(a.positions.len()).expect("fits u32")),
            )
            .collect(),
    };

    let whole_census = edge_census(&whole);
    let two_census = edge_census(&two_chunk);
    let whole_unpaired: std::collections::HashSet<_> = unpaired(&whole_census).into_iter().collect();
    let open_edges = unpaired(&two_census)
        .into_iter()
        .filter(|e| !whole_unpaired.contains(e))
        .count();

    let identical = {
        let wk = triangle_keys(&whole.positions, &whole.indices, None, true);
        let tk = triangle_keys(&two_chunk.positions, &two_chunk.indices, None, true);
        wk == tk
    };

    // The positive control: chunk B marches Kuhn's six anchored on corner `2`,
    // which splits every x-face on `2–4` instead of `0–6`, so the two chunks
    // disagree about the diagonal of their shared face. That is P-3's crack, it
    // is the one reading in this run that has to be non-zero, and it is what
    // makes every zero in `open_edges` a measurement rather than a silence.
    let control_open_edges = {
        let mut mirrored = TetMarcher::new(Decomp::SixTet { anchor: 2 });
        let bb = mirrored.mesh(field, [n - split, n, n], [seam_x, origin[1], origin[2]], h);
        let aa = mesh_of(
            Arm::SixTetBench,
            field,
            [split + 1, n, n],
            origin,
            h,
        );
        let mixed = Mesh {
            positions: aa
                .positions
                .iter()
                .chain(bb.positions.iter())
                .copied()
                .collect(),
            indices: aa
                .indices
                .iter()
                .copied()
                .chain(
                    bb.indices
                        .iter()
                        .map(|i| i + u32::try_from(aa.positions.len()).expect("fits u32")),
                )
                .collect(),
        };
        let mixed_census = edge_census(&mixed);
        unpaired(&mixed_census)
            .into_iter()
            .filter(|e| !whole_unpaired.contains(e))
            .count()
    };

    // The seam plane's own sign patterns, so a zero above is not a zero over an
    // uncut plane.
    let ns = n as usize;
    let sx = split as usize;
    let mut seen = [false; 16];
    for z in 0..ns - 1 {
        for y in 0..ns - 1 {
            let mut mask = 0u8;
            for (k, [dy, dz]) in [[0usize, 0usize], [1, 0], [0, 1], [1, 1]].iter().enumerate() {
                let v = field.sample([coords[sx], coords[y + dy], coords[z + dz]]);
                if is_inside(v) {
                    mask |= 1 << k;
                }
            }
            seen[mask as usize] = true;
        }
    }
    let configurations = (1..15).filter(|&m| seen[m]).count();

    Seam {
        open_edges,
        boundary_open_edges: whole_unpaired.len(),
        identical,
        configurations,
        coords_exact,
        control_open_edges,
    }
}

// ─── one row of the run ─────────────────────────────────────────────────────

struct RowData {
    field: &'static str,
    resolution: u32,
    cell_size: f64,
    arm: Arm,
    measured: Measured,
    seam: Seam,
    hausdorff: f64,
    hausdorff_mesh_to_field: f64,
    hausdorff_field_to_mesh: f64,
    accuracy_forward_samples: u64,
    accuracy_reverse_samples: u64,
    accuracy_triangles: u64,
    accuracy_degenerate_triangles: u64,
    facts_cut_edges: usize,
    facts_order_sensitive: usize,
    grid_symmetric: bool,
    baseline: Option<(usize, usize, usize)>,
    p61_columns_match: bool,
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-100");
    common::experiment::run(prereg, |run| {
        let elements = group();
        verify_group(&elements);
        let baseline = p61_baseline();

        // The decompositions themselves, checked before anything is meshed.
        let six = six_tets(0);
        assert_eq!(six.len(), 6, "Kuhn's decomposition has six tetrahedra");
        for tet in &six {
            assert_eq!(tet[0], 0, "every Kuhn tetrahedron starts at corner 0");
            assert_eq!(tet[3], 7, "every Kuhn tetrahedron ends at corner 7");
        }

        // The seam control's premise, verified rather than assumed: the two
        // anchors must split the shared x-face on **different** diagonals, or
        // the control is a fixture that cannot fail. `anchor = 1` does not, and
        // the first version of this harness used it.
        {
            let diagonal = |anchor: u8| -> [u8; 2] {
                let mut faces: Vec<Vec<u8>> = Vec::new();
                for tet in six_tets(anchor) {
                    let on: Vec<u8> = tet.iter().copied().filter(|c| c & 1 == 0).collect();
                    if on.len() == 3 {
                        faces.push(on);
                    }
                }
                assert_eq!(faces.len(), 2, "a cube face carries two Kuhn triangles");
                let mut shared: Vec<u8> = faces[0]
                    .iter()
                    .copied()
                    .filter(|c| faces[1].contains(c))
                    .collect();
                shared.sort_unstable();
                assert_eq!(shared.len(), 2, "the two triangles share a diagonal");
                [shared[0], shared[1]]
            };
            assert_eq!(diagonal(0), [0, 6], "the shipped x-face diagonal is 0-6");
            assert_eq!(diagonal(1), [0, 6], "anchor 1 is x-face-compatible with 0");
            assert_eq!(diagonal(2), [2, 4], "anchor 2 is the crack the control needs");
        }
        let bary = barycentric_24_tets();
        assert_eq!(bary.len(), 24, "the barycentric split has 24 tetrahedra");
        {
            // Every tetrahedron is (centre, face centre, two adjacent corners of
            // that face), and the 24 are distinct.
            let mut seen = std::collections::HashSet::new();
            for tet in &bary {
                assert_eq!(tet[0] as usize, CENTRE_SLOT, "coned to the cell centre");
                assert!((8..14).contains(&(tet[1] as usize)), "on a face centre");
                assert!(tet[2] < 8 && tet[3] < 8, "two cube corners");
                let shared = (tet[2] ^ tet[3]).count_ones();
                assert_eq!(shared, 1, "the two corners are a cube edge");
                let mut k = [tet[1], tet[2].min(tet[3]), tet[2].max(tet[3])];
                k.sort_unstable();
                assert!(seen.insert(k), "duplicate tetrahedron in the 24");
            }
        }

        println!("\n-- P-100: five decompositions through one equivariance instrument --");
        println!(
            "{:<15} {:>3} {:<30} {:>4} {:>4} {:>7} {:>6} {:>5} {:>6}",
            "field", "n", "decomposition", "vex", "tex", "tris", "open", "cfgs", "ms"
        );

        let mut rows: Vec<RowData> = Vec::new();

        isomesh::for_each_reference_field!(f64, |field_name, field| {
            let (lo, hi) = field.domain();
            let l = hi[0];
            for k in 0..3 {
                assert_eq!(
                    lo[k].to_bits(),
                    (-hi[k]).to_bits(),
                    "{field_name}: domain is not the symmetric cube the fixtures assume"
                );
            }

            for fx in fixtures(l) {
                let facts = fixture_facts(&field, &fx);
                assert!(
                    facts.grid_symmetric,
                    "{field_name} at {}: the grid is not bit-exactly closed under a \
                     sign flip, so the relation would be falsified by the fixture",
                    fx.samples
                );

                for arm in ARMS {
                    let measured = measure(arm, &field, &fx, &elements);
                    let seam = seam_census(arm, &field, &fx);

                    let shape = RuntimeShape3::new([fx.samples; 3]).expect("grid fits");
                    let mesh = mesh_of(arm, &field, [fx.samples; 3], [fx.origin; 3], fx.cell_size);
                    let cfg = AccuracyConfig::from_cell_size(fx.cell_size).expect("cell size");
                    let acc = accuracy(
                        &mesh.positions,
                        &mesh.indices,
                        &field,
                        &shape,
                        [fx.origin; 3],
                        &cfg,
                    )
                    .expect("accuracy");

                    // `hausdorff` reads exactly 0 on `box_exact` at 33³, and a
                    // zero has to prove it could have been non-zero: both
                    // directions of the symmetric distance must have been
                    // sampled, and the triangle set the index was built over
                    // must be the mesh's. Without this the column would be
                    // indistinguishable from an empty measurement.
                    assert!(
                        acc.mesh_to_field.samples > 0 && acc.field_to_mesh.samples > 0,
                        "{field_name}/{}/{}: accuracy sampled {} forward and {} \
                         reverse points, so its hausdorff is a zero over an empty \
                         set rather than a measurement",
                        fx.samples,
                        arm.name(),
                        acc.mesh_to_field.samples,
                        acc.field_to_mesh.samples
                    );
                    assert_eq!(
                        acc.triangles as usize + acc.degenerate_triangles as usize,
                        mesh.triangles(),
                        "{field_name}/{}/{}: accuracy did not see the whole mesh",
                        fx.samples,
                        arm.name()
                    );

                    // The registered vacuity control, per row: the six-tet arm
                    // against the committed file.
                    let p61_key = match arm {
                        Arm::MarchingCubes => Some(format!(
                            "{field_name}/marching_cubes/{}",
                            fx.samples
                        )),
                        Arm::SixTetCrate => Some(format!(
                            "{field_name}/marching_tetrahedra/{}",
                            fx.samples
                        )),
                        _ => None,
                    };
                    let (base, p61_columns_match) = match p61_key {
                        Some(k) => {
                            let b = baseline
                                .get(&k)
                                .unwrap_or_else(|| panic!("p-61.csv has no row {k}"));
                            let matched = b.cut_edges == facts.cut_edges
                                && b.order_sensitive_edges == facts.order_sensitive_edges
                                && b.elements_vertex_exact == measured.vertex_exact
                                && b.triangles == measured.triangles;
                            assert!(
                                b.cut_edges == facts.cut_edges
                                    && b.order_sensitive_edges == facts.order_sensitive_edges,
                                "{k}: the grid columns disagree with p-61.csv, so this is \
                                 not P-61's instrument ({} vs {}, {} vs {})",
                                facts.cut_edges,
                                b.cut_edges,
                                facts.order_sensitive_edges,
                                b.order_sensitive_edges
                            );
                            assert_eq!(
                                measured.vertex_exact, b.elements_vertex_exact,
                                "{k}: the arm does not reproduce p-61.csv's \
                                 elements_vertex_exact, so no 24-tet number here is \
                                 believable"
                            );
                            assert_eq!(
                                measured.triangles, b.triangles,
                                "{k}: triangle count disagrees with p-61.csv"
                            );
                            (
                                Some((
                                    b.elements_vertex_exact,
                                    b.cut_edges,
                                    b.triangles,
                                )),
                                matched,
                            )
                        }
                        None => (None, true),
                    };

                    println!(
                        "{:<15} {:>3} {:<30} {:>4} {:>4} {:>7} {:>6} {:>5} {:>6}",
                        field_name,
                        fx.samples,
                        arm.name(),
                        measured.vertex_exact,
                        measured.triangle_exact,
                        measured.triangles,
                        seam.open_edges,
                        seam.configurations,
                        measured.wall_ms
                    );

                    rows.push(RowData {
                        field: field_name,
                        resolution: fx.samples,
                        cell_size: fx.cell_size,
                        arm,
                        measured,
                        seam,
                        hausdorff: acc.symmetric_hausdorff(),
                        hausdorff_mesh_to_field: acc.mesh_to_field.max,
                        hausdorff_field_to_mesh: acc.field_to_mesh.max,
                        accuracy_forward_samples: acc.mesh_to_field.samples,
                        accuracy_reverse_samples: acc.field_to_mesh.samples,
                        accuracy_triangles: acc.triangles,
                        accuracy_degenerate_triangles: acc.degenerate_triangles,
                        facts_cut_edges: facts.cut_edges,
                        facts_order_sensitive: facts.order_sensitive_edges,
                        grid_symmetric: facts.grid_symmetric,
                        baseline: base,
                        p61_columns_match,
                    });
                }
            }
        });

        // ── the denominators, and the clause verdicts ───────────────────────
        let tri_of = |field: &str, res: u32, arm: Arm| -> usize {
            rows.iter()
                .find(|r| r.field == field && r.resolution == res && r.arm == arm)
                .map(|r| r.measured.triangles)
                .unwrap_or_else(|| panic!("no {} row for {field}/{res}", arm.name()))
        };

        let mut c1_rows_at_48 = 0usize;
        let mut c1_population = 0usize;
        let mut c1_worst = usize::MAX;
        let mut c2_worst_six = 0.0f64;
        let mut c2_worst_mc = 0.0f64;
        let mut c3_open_total = 0usize;
        let mut six_tet_rows_at_12 = 0usize;
        let mut mc_rows_at_48 = 0usize;
        let mut stabiliser_rows = 0usize;
        let mut control_min = usize::MAX;

        for r in &rows {
            let six = tri_of(r.field, r.resolution, Arm::SixTetCrate) as f64;
            let mc = tri_of(r.field, r.resolution, Arm::MarchingCubes) as f64;
            let ratio_six = r.measured.triangles as f64 / six;
            let ratio_mc = r.measured.triangles as f64 / mc;
            control_min = control_min.min(r.seam.control_open_edges);
            match r.arm {
                Arm::Barycentric24 => {
                    c1_population += 1;
                    if r.measured.vertex_exact == GROUP_ORDER {
                        c1_rows_at_48 += 1;
                    }
                    c1_worst = c1_worst.min(r.measured.vertex_exact);
                    c2_worst_six = c2_worst_six.max(ratio_six);
                    c2_worst_mc = c2_worst_mc.max(ratio_mc);
                    c3_open_total += r.seam.open_edges;
                }
                Arm::SixTetCrate | Arm::SixTetBench => {
                    if r.measured.vertex_exact == 12 {
                        six_tet_rows_at_12 += 1;
                    }
                    if r.measured.exact_are_diagonal_stabiliser {
                        stabiliser_rows += 1;
                    }
                }
                Arm::MarchingCubes => {
                    if r.measured.vertex_exact == GROUP_ORDER {
                        mc_rows_at_48 += 1;
                    }
                }
                Arm::Barycentric24Sampled => {}
            }
        }

        // The instrument brackets itself before any verdict is read.
        assert_eq!(
            mc_rows_at_48, 16,
            "marching_cubes must read 48 on all 16 rows or this instrument cannot \
             report a 48 at all and C1 is unreachable"
        );
        assert_eq!(
            six_tet_rows_at_12, 32,
            "both six-tet arms must read 12 on all 16 rows -- the registered \
             vacuity control -- or this instrument is not measuring what p-61.csv \
             measured"
        );
        assert!(
            control_min > 0,
            "the mismatched-diagonal seam control found no open edge, so the \
             open-edge counter cannot report a crack and every zero in C3 is vacuous"
        );

        let c1_holds = c1_rows_at_48 == c1_population;
        let c2_holds = c2_worst_six <= 2.0 && c2_worst_mc <= 6.0;
        let c3_holds = c3_open_total == 0;

        println!("\n-- verdicts --");
        println!(
            "C1  barycentric_24 at 48: {c1_rows_at_48}/{c1_population}, worst row \
             {c1_worst} of 48 -> {}",
            if c1_holds { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C2  worst triangles_vs_six_tet {c2_worst_six:.4} (bar 2.0), worst \
             triangles_vs_marching_cubes {c2_worst_mc:.4} (bar 6.0) -> {}",
            if c2_holds { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C3  open_edges total {c3_open_total}, control found {control_min}+ per \
             row -> {}",
            if c3_holds { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "vacuity  marching_cubes 48 on {mc_rows_at_48}/16, six-tet 12 on \
             {six_tet_rows_at_12}/32, diagonal-stabiliser rows {stabiliser_rows}/32"
        );

        for r in &rows {
            let six = tri_of(r.field, r.resolution, Arm::SixTetCrate) as f64;
            let mc = tri_of(r.field, r.resolution, Arm::MarchingCubes) as f64;
            let (b_vex, b_cut, b_tris) = match r.baseline {
                Some((v, c, t)) => (v.to_string(), c.to_string(), t.to_string()),
                None => (String::new(), String::new(), String::new()),
            };
            run.record(&[
                ("field", r.field.to_string()),
                ("resolution", r.resolution.to_string()),
                ("decomposition", r.arm.name().to_string()),
                ("tetrahedra_per_cell", r.arm.tets_per_cell().to_string()),
                (
                    "elements_vertex_exact",
                    r.measured.vertex_exact.to_string(),
                ),
                (
                    "worst_component_ulp",
                    r.measured.worst_component_ulp.to_string(),
                ),
                ("triangles", r.measured.triangles.to_string()),
                (
                    "triangles_vs_six_tet",
                    format!("{:.6}", r.measured.triangles as f64 / six),
                ),
                (
                    "triangles_vs_marching_cubes",
                    format!("{:.6}", r.measured.triangles as f64 / mc),
                ),
                ("open_edges", r.seam.open_edges.to_string()),
                (
                    "seam_configurations",
                    r.seam.configurations.to_string(),
                ),
                ("hausdorff", format!("{:.9}", r.hausdorff)),
                (
                    "hausdorff_mesh_to_field",
                    format!("{:.9}", r.hausdorff_mesh_to_field),
                ),
                (
                    "hausdorff_field_to_mesh",
                    format!("{:.9}", r.hausdorff_field_to_mesh),
                ),
                (
                    "accuracy_forward_samples",
                    r.accuracy_forward_samples.to_string(),
                ),
                (
                    "accuracy_reverse_samples",
                    r.accuracy_reverse_samples.to_string(),
                ),
                ("accuracy_triangles", r.accuracy_triangles.to_string()),
                (
                    "accuracy_degenerate_triangles",
                    r.accuracy_degenerate_triangles.to_string(),
                ),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_holds.to_string()),
                ("c3_holds", c3_holds.to_string()),
                // ── beyond the registration ──
                ("vertices", r.measured.vertices.to_string()),
                ("cell_size", format!("{:.9}", r.cell_size)),
                ("elements_tested", GROUP_ORDER.to_string()),
                (
                    "elements_triangle_exact",
                    r.measured.triangle_exact.to_string(),
                ),
                (
                    "first_failing_element",
                    r.measured.first_failing_element.clone(),
                ),
                (
                    "vertex_failing_labels",
                    r.measured.failing_labels.clone(),
                ),
                (
                    "exact_are_diagonal_stabiliser",
                    r.measured.exact_are_diagonal_stabiliser.to_string(),
                ),
                ("cut_edges", r.facts_cut_edges.to_string()),
                (
                    "order_sensitive_edges",
                    r.facts_order_sensitive.to_string(),
                ),
                (
                    "fixture_can_fail",
                    (r.facts_order_sensitive > 0).to_string(),
                ),
                ("grid_symmetric", r.grid_symmetric.to_string()),
                ("elements_vertex_exact_p61", b_vex),
                ("cut_edges_p61", b_cut),
                ("triangles_p61", b_tris),
                ("p61_columns_match", r.p61_columns_match.to_string()),
                (
                    "boundary_open_edges",
                    r.seam.boundary_open_edges.to_string(),
                ),
                (
                    "seam_control_open_edges",
                    r.seam.control_open_edges.to_string(),
                ),
                ("seam_mesh_identical", r.seam.identical.to_string()),
                ("seam_coords_exact", r.seam.coords_exact.to_string()),
                ("wall_ms", r.measured.wall_ms.to_string()),
                // ── run aggregates, repeated on every row ──
                ("c1_population", c1_population.to_string()),
                ("c1_rows_at_48", c1_rows_at_48.to_string()),
                ("c1_worst_row", c1_worst.to_string()),
                ("c2_worst_ratio_six_tet", format!("{c2_worst_six:.6}")),
                (
                    "c2_worst_ratio_marching_cubes",
                    format!("{c2_worst_mc:.6}"),
                ),
                ("c3_open_edges_total", c3_open_total.to_string()),
                (
                    "vacuity_six_tet_rows_at_12",
                    six_tet_rows_at_12.to_string(),
                ),
                ("vacuity_mc_rows_at_48", mc_rows_at_48.to_string()),
                (
                    "vacuity_stabiliser_rows",
                    stabiliser_rows.to_string(),
                ),
            ]);
        }
    });
}
