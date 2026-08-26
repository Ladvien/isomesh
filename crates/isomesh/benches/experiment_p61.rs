//! **P-61 — the crossing as a signed offset from the edge midpoint.**
//!
//! Ticket: R-059. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p61
//! ```
//!
//! Writes `docs/experiments/p-61.csv`.
//!
//! # Five blocks, one CSV, and the `block` column says which
//!
//! | `block` | clause | what it measures |
//! |---|---|---|
//! | `premeasure` | the appendix | the off-repo mirror count, re-run in `Real` for both scalars |
//! | `equivariance` | **C1** | `✗39`'s harness, re-run, against `p-57.csv` row for row |
//! | `edges` | **C2** | cut grid edges whose world position moved, per fixture |
//! | `geometry` | **C3** | Hausdorff and self-intersections, both placements, same topology |
//! | `seam` | **C4** | the two-chunk seam census at three spacings |
//!
//! The registered aggregates are repeated on **every** row, which is why the
//! rows are buffered and written at the end: a clause verdict is a property of
//! the whole run and a reader should not have to filter to find it. That is
//! `p-59.csv`'s shape.
//!
//! # The instrument is held against the old one where they overlap
//!
//! The equivariance block re-implements P-57's octahedral group, its two
//! fixtures and its comparison keys. A second copy of an instrument drifts, so
//! this one is **checked rather than trusted**: `cut_edges`,
//! `order_sensitive_edges`, `grid_symmetric` and `fixture_can_fail` are asserted
//! against the committed `docs/experiments/p-57.csv`, row for row, before any new
//! number is reported. Those four are properties of the grid and the field
//! values alone — no extractor touches them — so the placement change cannot
//! move them and a mismatch means this harness is not P-57's.
//!
//! `p57_fixture_columns_match` records the outcome so the CSV carries it too.
//!
//! # Why C3 compares two position arrays rather than two runs
//!
//! A placement change cannot touch the sign classification, so the case index,
//! the triangle list and the index buffer are bit-identical under both forms.
//! The two arms are therefore **the same mesh with substituted positions**, and
//! the substitution is exact: walk the grid, and for every cut edge compute both
//! `lo + (hi − lo)·t` and `(lo + hi)/2 + (hi − lo)·d`, keyed on the shipped
//! position's bits. Running an old build would compare two binaries, which
//! M-281 says is a comparison of layouts and not of formulas.
//!
//! The substitution covers vertices that sit **on a grid edge**, which is
//! Marching Cubes' whole vertex set apart from the interior-ambiguity apex. So
//! C3 is scoped to `marching_cubes` with the interior rule off, and the
//! registration says so.

#![allow(clippy::float_cmp)]

mod common;

use std::collections::HashMap;
use std::path::PathBuf;
use std::time::Instant;

use isomesh::extractor::Extractor;
use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{AccuracyConfig, accuracy, self_intersections};
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

// ─── the two placements, side by side, in one place ──────────────────────────

/// The lower-corner parameter this crate shipped until R-059: `t = a/(a − b)`.
///
/// Kept **here** and nowhere in `src/`, because a second placement in the
/// library would be the two-paths defect the crate's own rules forbid. A bench
/// comparing an old formula against a new one is the one place the old formula
/// belongs, which is how P-48 and P-54 carried their reference arithmetic.
#[inline]
fn lower_corner_t(a: f64, b: f64) -> f64 {
    a / (a - b)
}

/// `lo + (hi − lo)·t`, the placement that went with it.
#[inline]
fn place_lower_corner(lo: f64, hi: f64, t: f64) -> f64 {
    lo + (hi - lo) * t
}

/// The shipped form, spelled out generically so the pre-measurement can run it
/// in `f32` as well as `f64`. Identical expression to `cube::edge_offset`.
#[inline]
fn centred_d<R: Real>(a: R, b: R) -> R {
    ((a + b) * R::HALF) / (a - b)
}

/// The shipped placement. Identical expression to `cube::place`.
#[inline]
fn place_centred<R: Real>(lo: R, hi: R, d: R) -> R {
    (lo + hi) * R::HALF + (hi - lo) * d
}

/// The lower-corner pair, generic, for the pre-measurement's other arm.
#[inline]
fn generic_t<R: Real>(a: R, b: R) -> R {
    a / (a - b)
}

#[inline]
fn generic_place_lower<R: Real>(lo: R, hi: R, t: R) -> R {
    lo + (hi - lo) * t
}

// ─── a deterministic PRNG, because the numbers have to be re-runnable ────────

/// SplitMix64. Ten lines, no dependency, and a fixed seed, so the mismatch
/// counts below are the same on every machine and in every run — which is what
/// makes them a measurement rather than a sample.
struct SplitMix(u64);

impl SplitMix {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `(0, 1]`, never zero — a zero endpoint is not a straddling
    /// pair and the appendix's `random.uniform(0, 1)` excludes it in practice.
    fn unit(&mut self) -> f64 {
        let v = (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0);
        if v == 0.0 { 1.0 } else { v }
    }
}

// ─── the group ──────────────────────────────────────────────────────────────

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
    for e in g {
        match e.det() {
            1 => rotations += 1,
            -1 => reflections += 1,
            d => panic!("{} has determinant {d}, not +/-1", e.label()),
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

// ─── bit keys ───────────────────────────────────────────────────────────────

const NEGATIVE_ZERO: u64 = 1u64 << 63;

#[inline]
fn key(v: f64, normalise: bool) -> u64 {
    let b = v.to_bits();
    if normalise && b == NEGATIVE_ZERO {
        0
    } else {
        b
    }
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

/// Multiset symmetric difference of two **sorted** key lists. See P-57's own
/// note: a positional walk reports the `+x`/`−x` gap the moment the lists differ
/// by an insertion rather than a perturbation.
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
    for tri in indices.as_chunks::<3>().0 {
        let c = [mapped(tri[0]), mapped(tri[1]), mapped(tri[2])];
        let rotations = [[c[0], c[1], c[2]], [c[1], c[2], c[0]], [c[2], c[0], c[1]]];
        out.push(rotations.into_iter().min().expect("three rotations"));
    }
    out.sort_unstable();
    out
}

// ─── fixtures ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct Fixture {
    samples: u32,
    cell_size: f64,
    origin: f64,
}

/// P-57's two fixtures, unchanged. The 25³ arm uses `3L/32` rather than the
/// crate's `2L/24`, because `L/12` is not dyadic and its grid does not mirror.
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

struct FixtureFacts {
    grid_symmetric: bool,
    zero_coordinates: usize,
    cut_edges: usize,
    order_sensitive_edges: usize,
    boundary_sign_constant: bool,
    /// Cut edges whose **centred** world coordinate differs in bits from the
    /// lower-corner one. **C2's population, per fixture.**
    edges_moved: usize,
    /// The largest ULP gap over those.
    worst_move_ulp: u128,
}

fn axis_coords(fx: &Fixture) -> Vec<f64> {
    (0..fx.samples)
        .map(|i| fx.origin + fx.cell_size * f64::from(i))
        .collect()
}

/// Grid symmetry, order-sensitivity and the placement delta, in one pass.
///
/// `order_sensitive_edges` is P-57's own census and is deliberately still stated
/// in the **lower-corner** frame: it is a property of the grid coordinates and
/// the field values, it is what `fixture_can_fail` is derived from, and changing
/// its definition would silently redefine the population C1 is measured over.
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

    let mut cut_edges = 0usize;
    let mut order_sensitive_edges = 0usize;
    let mut edges_moved = 0usize;
    let mut worst_move_ulp = 0u128;
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
                    let fwd = place_lower_corner(lo_c, hi_c, lower_corner_t(va, vb));
                    let rev = place_lower_corner(hi_c, lo_c, lower_corner_t(vb, va));
                    if key(fwd, true) != key(rev, true) {
                        order_sensitive_edges += 1;
                    }
                    let centred = place_centred(lo_c, hi_c, centred_d(va, vb));
                    if key(centred, true) != key(fwd, true) {
                        edges_moved += 1;
                        worst_move_ulp =
                            worst_move_ulp.max(ulp_distance(centred.to_bits(), fwd.to_bits()));
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
        edges_moved,
        worst_move_ulp,
    }
}

// ─── P-57's committed artefact, as the baseline ──────────────────────────────

/// One row of `docs/experiments/p-57.csv`, keyed `field/extractor/samples`.
struct Baseline {
    elements_vertex_exact: usize,
    cut_edges: usize,
    order_sensitive_edges: usize,
    grid_symmetric: bool,
    fixture_can_fail: bool,
}

/// Read `p-57.csv`. **Not optional**: without it there is no before-arm for C1
/// and no check that this harness is P-57's instrument.
fn p57_baseline() -> HashMap<String, Baseline> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/experiments/p-57.csv");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("P-61 needs {} as its baseline: {e}", path.display()));
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty());
    let header: Vec<&str> = lines
        .next()
        .expect("p-57.csv has a header")
        .split(',')
        .collect();
    let col = |name: &str| {
        header
            .iter()
            .position(|h| *h == name)
            .unwrap_or_else(|| panic!("p-57.csv has no `{name}` column"))
    };
    let (c_field, c_extractor, c_samples) =
        (col("field"), col("extractor"), col("samples_per_axis"));
    let c_vex = col("elements_vertex_exact");
    let c_cut = col("cut_edges");
    let c_ose = col("order_sensitive_edges");
    let c_sym = col("grid_symmetric");
    let c_canfail = col("fixture_can_fail");

    let mut out = HashMap::new();
    for line in lines {
        let f: Vec<&str> = line.split(',').collect();
        let k = format!("{}/{}/{}", f[c_field], f[c_extractor], f[c_samples]);
        out.insert(
            k,
            Baseline {
                elements_vertex_exact: f[c_vex].parse().expect("integer"),
                cut_edges: f[c_cut].parse().expect("integer"),
                order_sensitive_edges: f[c_ose].parse().expect("integer"),
                grid_symmetric: f[c_sym] == "true",
                fixture_can_fail: f[c_canfail] == "true",
            },
        );
    }
    assert_eq!(out.len(), 112, "p-57.csv should carry 112 rows");
    out
}

// ─── the equivariance measurement ───────────────────────────────────────────

struct Measured {
    vertices: usize,
    triangles: usize,
    vertex_exact: usize,
    triangle_exact: usize,
    vertex_exact_triangle_fail: usize,
    pure_permutation_exact: usize,
    pure_sign_flip_exact: usize,
    worst_differing_vertices: usize,
    first_failing_element: String,
    worst_component_ulp: u128,
    vertex_failing_labels: String,
    wall_ms: u128,
}

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
        triangle_exact: 0,
        vertex_exact_triangle_fail: 0,
        pure_permutation_exact: 0,
        pure_sign_flip_exact: 0,
        worst_differing_vertices: 0,
        first_failing_element: String::from("none"),
        worst_component_ulp: 0,
        vertex_failing_labels: String::new(),
        wall_ms: 0,
    };
    let mut failing: Vec<String> = Vec::new();
    let mut any_failed = false;

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

        let got_tri = triangle_keys(&rotated.positions, &rotated.indices, None, true);
        let want_tri = triangle_keys(&reference.positions, &reference.indices, Some(g), true);
        let triangle_ok = got_tri == want_tri;

        if vertex_ok {
            m.vertex_exact += 1;
            if g.sign == [1, 1, 1] {
                m.pure_permutation_exact += 1;
            }
            if g.perm == [0, 1, 2] {
                m.pure_sign_flip_exact += 1;
            }
        } else {
            let (only_got, only_want) = multiset_difference(&got, &want);
            m.worst_differing_vertices = m.worst_differing_vertices.max(only_got.len());
            for (a, b) in only_got.iter().zip(only_want.iter()) {
                for k in 0..3 {
                    if a[k] != b[k] {
                        m.worst_component_ulp = m.worst_component_ulp.max(ulp_distance(a[k], b[k]));
                    }
                }
            }
            if !any_failed {
                m.first_failing_element = g.label();
                any_failed = true;
            }
            failing.push(g.short());
        }
        if triangle_ok {
            m.triangle_exact += 1;
        } else if vertex_ok {
            m.vertex_exact_triangle_fail += 1;
        }

        // **The control**: element 0 is the identity and must reproduce the
        // reference exactly, or the extractor is not deterministic and the row
        // measures nothing.
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

// ─── the appendix, re-run in `Real` ─────────────────────────────────────────

/// One row of the pre-measurement.
struct PreRow {
    scalar: &'static str,
    frame: &'static str,
    cell_size: f64,
    pairs: u64,
    mismatches_lower_corner: u64,
    mismatches_centred: u64,
    out_of_cell_centred: u64,
}

/// The appendix script, in `Real`.
///
/// `cell_local` compares the two forms inside a unit cell, `world` inside cell
/// `i` of a grid at `origin + h·i` — which is the frame that matters, because a
/// mesher never sees a bare parameter. The mirrored arm swaps the endpoint roles
/// **and** negates the coordinates, which is what a reflection of the grid does.
fn premeasure<R: Real>(scalar: &'static str, frame: &'static str, h: f64, pairs: u64) -> PreRow {
    let mut rng = SplitMix(2026);
    let (mut bad_lower, mut bad_centred, mut out_of_cell) = (0u64, 0u64, 0u64);
    let origin = R::from_f64(-2.0);
    let hr = R::from_f64(h);
    for _ in 0..pairs {
        let i = (rng.next_u64() % 32) as u32;
        let a = R::from_f64(rng.unit());
        let b = -R::from_f64(rng.unit());

        let (lo, up) = if frame == "world" {
            (
                origin + hr * R::from_f64(f64::from(i)),
                origin + hr * R::from_f64(f64::from(i + 1)),
            )
        } else {
            (R::ZERO, R::ONE)
        };

        // Forward, and the mirrored cell: the corner roles swap and every
        // coordinate negates.
        let x = generic_place_lower(lo, up, generic_t(a, b));
        let x2 = generic_place_lower(-up, -lo, generic_t(b, a));
        if x2 != -x {
            bad_lower += 1;
        }

        let y = place_centred(lo, up, centred_d(a, b));
        let y2 = place_centred(-up, -lo, centred_d(b, a));
        if y2 != -y {
            bad_centred += 1;
        }
        if y < lo || y > up {
            out_of_cell += 1;
        }
    }
    PreRow {
        scalar,
        frame,
        cell_size: h,
        pairs,
        mismatches_lower_corner: bad_lower,
        mismatches_centred: bad_centred,
        out_of_cell_centred: out_of_cell,
    }
}

// ─── C3: the same mesh, two position arrays ─────────────────────────────────

/// Positions under both placements for one `(field, fixture)` pair, keyed by the
/// shipped position's bits.
///
/// Every Marching Cubes edge vertex is `place(lo, hi, d)` for exactly one grid
/// edge, so walking the grid rebuilds the whole map. The interior-ambiguity apex
/// is not on a grid edge, which is why C3 runs with the interior rule **off**.
fn placement_map<S: Sdf<Scalar = f64>>(field: &S, fx: &Fixture) -> HashMap<[u64; 3], [f64; 3]> {
    let n = fx.samples as usize;
    let coords = axis_coords(fx);
    let mut values = Vec::with_capacity(n * n * n);
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                values.push(field.sample([coords[x], coords[y], coords[z]]));
            }
        }
    }
    let at = |x: usize, y: usize, z: usize| values[x + n * (y + n * z)];

    let mut map = HashMap::new();
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
                    let lo_w = [coords[x], coords[y], coords[z]];
                    let mut hi_w = lo_w;
                    hi_w[axis] = coords[hi[axis]];

                    let d = centred_d(va, vb);
                    let t = lower_corner_t(va, vb);
                    let new: [f64; 3] = std::array::from_fn(|k| place_centred(lo_w[k], hi_w[k], d));
                    let old: [f64; 3] =
                        std::array::from_fn(|k| place_lower_corner(lo_w[k], hi_w[k], t));
                    map.insert(std::array::from_fn(|k| new[k].to_bits()), old);
                }
            }
        }
    }
    map
}

/// What one C3 row says.
struct Geometry {
    hausdorff_new: f64,
    hausdorff_old: f64,
    si_new: f64,
    si_old: f64,
    vertices: usize,
    substituted: usize,
}

fn geometry<S: Sdf<Scalar = f64>>(field: &S, fx: &Fixture) -> Geometry {
    let shape = RuntimeShape3::new([fx.samples; 3]).expect("grid fits");
    let origin = [fx.origin; 3];
    let mut mc = MarchingCubes::<f64>::new();
    let mut mesh = MeshBuffer::<f64>::new();
    mc.extract_into(field, &shape, origin, fx.cell_size, &mut mesh)
        .expect("extraction");

    let map = placement_map(field, fx);
    let mut old_positions = mesh.positions.clone();
    let mut substituted = 0usize;
    for p in &mut old_positions {
        let k: [u64; 3] = std::array::from_fn(|i| p[i].to_bits());
        if let Some(q) = map.get(&k) {
            *p = *q;
            substituted += 1;
        }
    }

    let cfg = AccuracyConfig::from_cell_size(fx.cell_size).expect("valid cell size");
    let acc_new = accuracy(&mesh.positions, &mesh.indices, field, &shape, origin, &cfg)
        .expect("accuracy, centred");
    let acc_old = accuracy(&old_positions, &mesh.indices, field, &shape, origin, &cfg)
        .expect("accuracy, lower corner");

    let per_1k = |report: isomesh::validate::SelfIntersectionReport| {
        if report.triangles == 0 {
            0.0
        } else {
            1000.0 * report.pairs.len() as f64 / report.triangles as f64
        }
    };
    let si_new = per_1k(
        self_intersections(&mesh.positions, &mesh.indices, fx.cell_size)
            .expect("self-intersections, centred"),
    );
    let si_old = per_1k(
        self_intersections(&old_positions, &mesh.indices, fx.cell_size)
            .expect("self-intersections, lower corner"),
    );

    Geometry {
        hausdorff_new: acc_new.symmetric_hausdorff(),
        hausdorff_old: acc_old.symmetric_hausdorff(),
        si_new,
        si_old,
        vertices: mesh.positions.len(),
        substituted,
    }
}

// ─── C4: the two-chunk seam census ──────────────────────────────────────────

/// One seam row.
struct Seam {
    cell_size: f64,
    vertices: usize,
    mismatches_lower_corner: usize,
    mismatches_centred: usize,
}

/// Two neighbouring chunks compute the seam plane's crossings from their own
/// origins, which is `M-32`'s mechanism, and this counts the ones that disagree
/// in bits — under both placements.
///
/// The two arms differ only in the placement, so a difference here would mean
/// the crossing form participates in the seam and a zero means it does not. The
/// registration expects the latter: `M-32` names `world_of_sample`.
fn seam_census<S: Sdf<Scalar = f64>>(field: &S, h: f64, cells: u32) -> Seam {
    let n = cells as usize + 1;
    // Chunk A spans samples `0..=cells` from origin `-2`, chunk B the next span,
    // each reconstructing its own sample positions the way an extractor does:
    // `chunk_origin + h * local`.
    let origin_a = -2.0f64;
    let origin_b = origin_a + h * f64::from(cells);
    let coord = |origin: f64, i: usize| origin + h * (i as f64);

    let mut vertices = 0usize;
    let (mut bad_lower, mut bad_centred) = (0usize, 0usize);
    // The seam is chunk A's last sample plane and chunk B's first: A reaches it
    // as local index `cells`, B as local index 0.
    for y in 0..n {
        for z in 0..n {
            // The y-axis edge on the seam plane, from each side.
            if y + 1 >= n {
                continue;
            }
            let pa = [
                coord(origin_a, cells as usize),
                coord(origin_a, y),
                coord(origin_a, z),
            ];
            let pa_hi = [pa[0], coord(origin_a, y + 1), pa[2]];
            let pb = [coord(origin_b, 0), coord(origin_b, y), coord(origin_b, z)];
            let pb_hi = [pb[0], coord(origin_b, y + 1), pb[2]];

            let va = field.sample(pa);
            let vb = field.sample(pa_hi);
            if (va < 0.0) == (vb < 0.0) {
                continue;
            }
            // The same edge as B sees it. Its field values are the same numbers
            // only if the coordinates agree bit-for-bit, which is the property
            // under test, so each side samples its own reconstruction.
            let wa = field.sample(pb);
            let wb = field.sample(pb_hi);
            if (wa < 0.0) == (wb < 0.0) {
                continue;
            }
            vertices += 1;

            let ta = lower_corner_t(va, vb);
            let tb = lower_corner_t(wa, wb);
            for k in 0..3 {
                if key(place_lower_corner(pa[k], pa_hi[k], ta), true)
                    != key(place_lower_corner(pb[k], pb_hi[k], tb), true)
                {
                    bad_lower += 1;
                    break;
                }
            }
            let da = centred_d(va, vb);
            let db = centred_d(wa, wb);
            for k in 0..3 {
                if key(place_centred(pa[k], pa_hi[k], da), true)
                    != key(place_centred(pb[k], pb_hi[k], db), true)
                {
                    bad_centred += 1;
                    break;
                }
            }
        }
    }
    Seam {
        cell_size: h,
        vertices,
        mismatches_lower_corner: bad_lower,
        mismatches_centred: bad_centred,
    }
}

// ─── the run ────────────────────────────────────────────────────────────────

/// Every row, before the aggregates are known.
type Row = Vec<(&'static str, String)>;

const NA: &str = "";

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-61");
    common::experiment::run(prereg, |run| {
        let elements = group();
        verify_group(&elements);
        let baseline = p57_baseline();

        let mut rows: Vec<Row> = Vec::new();

        // ── block: premeasure ───────────────────────────────────────────────
        println!("\n-- premeasure: the appendix, in Real, both scalars --");
        println!(
            "{:<7} {:<11} {:>10} {:>10} {:>12} {:>10} {:>8}",
            "scalar", "frame", "h", "pairs", "lower_corner", "centred", "outside"
        );
        let mut pre_lower_max = 0u64;
        let mut pre_centred_max = 0u64;
        let mut pre_rows: Vec<PreRow> = Vec::new();
        for (frame, h, pairs) in [
            ("cell_local", 1.0, 2_000_000u64),
            ("world", 0.125, 300_000),
            ("world", 0.1, 300_000),
            ("world", 3.0 / 32.0, 300_000),
        ] {
            pre_rows.push(premeasure::<f64>("f64", frame, h, pairs));
            pre_rows.push(premeasure::<f32>("f32", frame, h, pairs));
        }
        for p in &pre_rows {
            println!(
                "{:<7} {:<11} {:>10.6} {:>10} {:>12} {:>10} {:>8}",
                p.scalar,
                p.frame,
                p.cell_size,
                p.pairs,
                p.mismatches_lower_corner,
                p.mismatches_centred,
                p.out_of_cell_centred
            );
            pre_lower_max = pre_lower_max.max(p.mismatches_lower_corner);
            pre_centred_max = pre_centred_max.max(p.mismatches_centred);
            rows.push(vec![
                ("block", "premeasure".to_string()),
                ("field", NA.to_string()),
                ("extractor", NA.to_string()),
                ("samples_per_axis", NA.to_string()),
                ("scalar", p.scalar.to_string()),
                ("cell_size", format!("{:.9}", p.cell_size)),
                ("pairs", p.pairs.to_string()),
                (
                    "mismatches_lower_corner",
                    p.mismatches_lower_corner.to_string(),
                ),
                ("mismatches_centred", p.mismatches_centred.to_string()),
                ("out_of_cell_centred", p.out_of_cell_centred.to_string()),
                ("frame", p.frame.to_string()),
            ]);
        }

        // ── block: equivariance, edges, geometry ────────────────────────────
        println!("\n-- equivariance: P-57's harness, re-run --");
        println!(
            "{:<15} {:<28} {:>3} {:>5} {:>5} {:>4} {:>3} {:>3} {:>8} {:>7}",
            "field", "extractor", "n", "vex", "p57", "tex", "p6", "f8", "ulp", "ms"
        );
        let mut reference = MeshBuffer::<f64>::new();
        let mut rotated = MeshBuffer::<f64>::new();
        let mut c1_population = 0usize;
        let mut c1_rows_at_48 = 0usize;
        let mut c1_can_fail_at_48 = 0usize;
        let mut fixture_columns_match = true;
        let mut c2_fixtures_with_moved_edges = 0usize;
        let mut c2_fixtures = 0usize;
        let mut worst_hausdorff_ratio = 0.0f64;
        let mut worst_si_ratio = 0.0f64;
        let mut edge_rows: Vec<Row> = Vec::new();
        let mut geometry_rows: Vec<Row> = Vec::new();

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

                // C2's row, per fixture.
                c2_fixtures += 1;
                if facts.edges_moved > 0 {
                    c2_fixtures_with_moved_edges += 1;
                }
                edge_rows.push(vec![
                    ("block", "edges".to_string()),
                    ("field", field_name.to_string()),
                    ("extractor", NA.to_string()),
                    ("samples_per_axis", fx.samples.to_string()),
                    ("scalar", "f64".to_string()),
                    ("cell_size", format!("{:.9}", fx.cell_size)),
                    ("cut_edges", facts.cut_edges.to_string()),
                    ("edges_moved", facts.edges_moved.to_string()),
                    ("worst_move_ulp", facts.worst_move_ulp.to_string()),
                    (
                        "order_sensitive_edges",
                        facts.order_sensitive_edges.to_string(),
                    ),
                    (
                        "edges_moved_fraction",
                        format!(
                            "{:.6}",
                            if facts.cut_edges == 0 {
                                0.0
                            } else {
                                facts.edges_moved as f64 / facts.cut_edges as f64
                            }
                        ),
                    ),
                ]);

                // C3's row, per fixture, on marching_cubes.
                let g = geometry(&field, &fx);
                let h_ratio = if g.hausdorff_old == 0.0 {
                    1.0
                } else {
                    g.hausdorff_new / g.hausdorff_old
                };
                let si_ratio = if g.si_old == 0.0 {
                    if g.si_new == 0.0 { 1.0 } else { f64::INFINITY }
                } else {
                    g.si_new / g.si_old
                };
                worst_hausdorff_ratio = worst_hausdorff_ratio.max((h_ratio - 1.0).abs());
                if si_ratio.is_finite() {
                    worst_si_ratio = worst_si_ratio.max((si_ratio - 1.0).abs());
                }
                geometry_rows.push(vec![
                    ("block", "geometry".to_string()),
                    ("field", field_name.to_string()),
                    ("extractor", "marching_cubes".to_string()),
                    ("samples_per_axis", fx.samples.to_string()),
                    ("scalar", "f64".to_string()),
                    ("cell_size", format!("{:.9}", fx.cell_size)),
                    (
                        "hausdorff_lower_corner",
                        format!("{:.12e}", g.hausdorff_old),
                    ),
                    ("hausdorff_centred", format!("{:.12e}", g.hausdorff_new)),
                    ("hausdorff_ratio", format!("{h_ratio:.9}")),
                    (
                        "self_intersections_per_1k_lower_corner",
                        format!("{:.6}", g.si_old),
                    ),
                    (
                        "self_intersections_per_1k_centred",
                        format!("{:.6}", g.si_new),
                    ),
                    ("self_intersections_ratio", format!("{si_ratio:.9}")),
                    ("vertices", g.vertices.to_string()),
                    ("vertices_substituted", g.substituted.to_string()),
                ]);

                isomesh::for_each_extractor!(f64, |extractor_name, extractor| {
                    let m = measure(
                        &field,
                        &mut extractor,
                        &fx,
                        &elements,
                        &mut reference,
                        &mut rotated,
                    );
                    let k = format!("{field_name}/{extractor_name}/{}", fx.samples);
                    let b = baseline
                        .get(&k)
                        .unwrap_or_else(|| panic!("p-57.csv has no row {k}"));
                    // The four fixture columns are properties of the grid and
                    // the field values, so the placement change cannot move
                    // them. A mismatch means this harness is not P-57's.
                    let matched = b.cut_edges == facts.cut_edges
                        && b.order_sensitive_edges == facts.order_sensitive_edges
                        && b.grid_symmetric == facts.grid_symmetric
                        && b.fixture_can_fail == (facts.order_sensitive_edges > 0);
                    assert!(
                        matched,
                        "{k}: this harness disagrees with p-57.csv about the FIXTURE \
                         (cut_edges {} vs {}, order_sensitive {} vs {}) -- the \
                         instrument drifted and no new number here means anything",
                        facts.cut_edges,
                        b.cut_edges,
                        facts.order_sensitive_edges,
                        b.order_sensitive_edges
                    );
                    fixture_columns_match &= matched;

                    if b.fixture_can_fail {
                        c1_population += 1;
                        if m.vertex_exact == GROUP_ORDER {
                            c1_can_fail_at_48 += 1;
                        }
                    }
                    if m.vertex_exact == GROUP_ORDER {
                        c1_rows_at_48 += 1;
                    }

                    println!(
                        "{:<15} {:<28} {:>3} {:>5} {:>5} {:>4} {:>3} {:>3} {:>8} {:>7}",
                        field_name,
                        extractor_name,
                        fx.samples,
                        m.vertex_exact,
                        b.elements_vertex_exact,
                        m.triangle_exact,
                        m.pure_permutation_exact,
                        m.pure_sign_flip_exact,
                        m.worst_component_ulp,
                        m.wall_ms
                    );
                    rows.push(vec![
                        ("block", "equivariance".to_string()),
                        ("field", field_name.to_string()),
                        ("extractor", extractor_name.to_string()),
                        ("samples_per_axis", fx.samples.to_string()),
                        ("scalar", "f64".to_string()),
                        ("cell_size", format!("{:.9}", fx.cell_size)),
                        ("family", family(extractor_name).to_string()),
                        ("elements_tested", GROUP_ORDER.to_string()),
                        ("elements_vertex_exact", m.vertex_exact.to_string()),
                        (
                            "elements_vertex_exact_p57",
                            b.elements_vertex_exact.to_string(),
                        ),
                        ("elements_triangle_exact", m.triangle_exact.to_string()),
                        (
                            "elements_vertex_exact_triangle_fail",
                            m.vertex_exact_triangle_fail.to_string(),
                        ),
                        (
                            "pure_permutation_exact",
                            m.pure_permutation_exact.to_string(),
                        ),
                        ("pure_sign_flip_exact", m.pure_sign_flip_exact.to_string()),
                        ("fixture_can_fail", b.fixture_can_fail.to_string()),
                        ("cut_edges", facts.cut_edges.to_string()),
                        (
                            "order_sensitive_edges",
                            facts.order_sensitive_edges.to_string(),
                        ),
                        ("edges_moved", facts.edges_moved.to_string()),
                        ("worst_move_ulp", facts.worst_move_ulp.to_string()),
                        ("vertices", m.vertices.to_string()),
                        ("triangles", m.triangles.to_string()),
                        (
                            "worst_differing_vertices",
                            m.worst_differing_vertices.to_string(),
                        ),
                        ("first_failing_element", m.first_failing_element.clone()),
                        ("worst_component_ulp", m.worst_component_ulp.to_string()),
                        ("vertex_failing_labels", m.vertex_failing_labels.clone()),
                        (
                            "boundary_sign_constant",
                            facts.boundary_sign_constant.to_string(),
                        ),
                        ("grid_symmetric", facts.grid_symmetric.to_string()),
                        ("grid_zero_coordinates", facts.zero_coordinates.to_string()),
                        ("wall_ms", m.wall_ms.to_string()),
                    ]);
                });
            }
        });
        rows.extend(edge_rows);
        rows.extend(geometry_rows);

        // ── block: seam ─────────────────────────────────────────────────────
        println!("\n-- seam: two chunks, three spacings, both placements --");
        println!(
            "{:>10} {:>10} {:>16} {:>10}",
            "h", "vertices", "lower_corner", "centred"
        );
        let mut c4_delta = 0i64;
        {
            let field = isomesh::fields::Gyroid::<f64>::canonical();
            for h in [0.125f64, 0.1, 3.0 / 32.0] {
                let s = seam_census(&field, h, 16);
                println!(
                    "{:>10.6} {:>10} {:>16} {:>10}",
                    s.cell_size, s.vertices, s.mismatches_lower_corner, s.mismatches_centred
                );
                c4_delta += s.mismatches_centred as i64 - s.mismatches_lower_corner as i64;
                rows.push(vec![
                    ("block", "seam".to_string()),
                    ("field", "gyroid".to_string()),
                    ("extractor", NA.to_string()),
                    ("samples_per_axis", "17".to_string()),
                    ("scalar", "f64".to_string()),
                    ("cell_size", format!("{:.9}", s.cell_size)),
                    ("seam_vertices", s.vertices.to_string()),
                    (
                        "seam_mismatches_lower_corner",
                        s.mismatches_lower_corner.to_string(),
                    ),
                    ("seam_mismatches_centred", s.mismatches_centred.to_string()),
                ]);
            }
        }

        // ── the aggregates, and the verdicts ────────────────────────────────
        let c1_holds = c1_can_fail_at_48 == c1_population;
        let c2_holds = c2_fixtures_with_moved_edges == c2_fixtures;
        let c3_holds = worst_hausdorff_ratio < 0.01 && worst_si_ratio < 0.01;

        assert!(
            c1_population > 0,
            "VOID: no fixture-can-fail row, so C1's population is empty and its \
             verdict would be vacuous (the audit's central finding)"
        );
        assert!(
            pre_lower_max > 0,
            "VOID: the lower-corner form mismatched nothing in the pre-measurement, \
             so the instrument cannot report the bad news it exists to report"
        );

        println!(
            "\nC1 {}/{} fixture-can-fail rows at 48 -> {}",
            c1_can_fail_at_48,
            c1_population,
            if c1_holds { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C2 {}/{} fixtures with moved edges -> {}",
            c2_fixtures_with_moved_edges,
            c2_fixtures,
            if c2_holds { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C3 worst |hausdorff ratio - 1| {:.6}, worst |self-intersection ratio - 1| {:.6} -> {}",
            worst_hausdorff_ratio,
            worst_si_ratio,
            if c3_holds { "HELD" } else { "FALSIFIED" }
        );
        println!("C4 seam mismatch delta {c4_delta} (a null is the registered expectation)");
        println!(
            "pre-measurement: worst lower-corner {pre_lower_max}, worst centred {pre_centred_max}"
        );

        let aggregates: Row = vec![
            ("c1_population", c1_population.to_string()),
            ("c1_rows_at_48", c1_rows_at_48.to_string()),
            ("c1_can_fail_rows_at_48", c1_can_fail_at_48.to_string()),
            ("c1_holds", c1_holds.to_string()),
            (
                "c2_fixtures_with_moved_edges",
                c2_fixtures_with_moved_edges.to_string(),
            ),
            ("c2_fixtures", c2_fixtures.to_string()),
            ("c2_holds", c2_holds.to_string()),
            (
                "c3_worst_hausdorff_ratio",
                format!("{worst_hausdorff_ratio:.9}"),
            ),
            (
                "c3_worst_self_intersection_ratio",
                format!("{worst_si_ratio:.9}"),
            ),
            ("c3_holds", c3_holds.to_string()),
            ("c4_seam_mismatch_delta", c4_delta.to_string()),
            (
                "p57_fixture_columns_match",
                fixture_columns_match.to_string(),
            ),
            ("pre_worst_lower_corner", pre_lower_max.to_string()),
            ("pre_worst_centred", pre_centred_max.to_string()),
        ];

        // Every registered column on every row: a reader should not have to
        // filter a CSV to find a clause verdict, and `Run::record` requires it.
        let registered: [&str; 38] = [
            "block",
            "field",
            "extractor",
            "samples_per_axis",
            "scalar",
            "cell_size",
            "elements_vertex_exact",
            "elements_vertex_exact_p57",
            "elements_triangle_exact",
            "pure_permutation_exact",
            "pure_sign_flip_exact",
            "fixture_can_fail",
            "cut_edges",
            "edges_moved",
            "worst_move_ulp",
            "hausdorff_lower_corner",
            "hausdorff_centred",
            "hausdorff_ratio",
            "self_intersections_per_1k_lower_corner",
            "self_intersections_per_1k_centred",
            "self_intersections_ratio",
            "seam_vertices",
            "seam_mismatches_lower_corner",
            "seam_mismatches_centred",
            "pairs",
            "mismatches_lower_corner",
            "mismatches_centred",
            "out_of_cell_centred",
            "c1_population",
            "c1_rows_at_48",
            "c1_holds",
            "c2_fixtures_with_moved_edges",
            "c2_holds",
            "c3_worst_hausdorff_ratio",
            "c3_worst_self_intersection_ratio",
            "c3_holds",
            "c4_seam_mismatch_delta",
            "p57_fixture_columns_match",
        ];
        for mut row in rows {
            row.extend(aggregates.iter().cloned());
            for name in registered {
                if !row.iter().any(|(k, _)| *k == name) {
                    row.push((name, NA.to_string()));
                }
            }
            run.record(&row);
        }
    });
}
