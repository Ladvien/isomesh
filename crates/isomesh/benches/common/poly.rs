//! Exact multivariate polynomial arithmetic, and the `2x2x2` hyperdeterminant
//! algebra that Group A of Phase 27 is built on.
//!
//! Ticket: R-127, which owns this module. Consumed unchanged by R-128, R-129,
//! R-135, R-137, R-138, R-143 and R-155 — a second copy of this algebra would
//! be two paths to one answer, and the answer is an identity.
//!
//! # The mechanism
//!
//! A cell's eight corner values `f[u + 2v + 4w]` are a `2x2x2` tensor. Cayley's
//! hyperdeterminant `Det(A)` is the unique (up to scale) `SL(2)^3`-invariant of
//! that tensor: a degree-4 form in the eight entries, twelve monomials, and the
//! discriminant of the quadratic-in-`lambda` pencil `det(A0 + lambda*A1)` for
//! *any* of the three ways of slicing the cube into two opposite faces. The
//! claim of P-127 is that `b*b - 4*a*c` at
//! `crates/isomesh/src/marching_cubes/trilinear.rs:246`, whose `a`, `b`, `c`
//! come from `BodySaddles::coefficients` at :199-214, **is** that invariant —
//! the crate has been shipping a classical object under a transcribed name.
//!
//! Everything here is exact. Coefficients are `i128`, monomial exponents are
//! `[u8; 8]`, and the term map is a `BTreeMap`, so two runs of the same
//! construction produce byte-identical output and a difference of zero means
//! zero rather than "smaller than a tolerance". Every arithmetic step is
//! checked: an overflow is a panic, never a wrapped answer wearing the word
//! *exact*.
//!
//! # Why `cayley_2x2x2` is written out rather than derived from the pencil
//!
//! [`cayley_2x2x2`] is the standard explicit twelve-term normalisation, the one
//! `docs/research/2026-08-29-phase-27-hyperdeterminant-identity.py:45-67`
//! writes in sympy. It could instead be *defined* as `pencil_discriminant(0)`,
//! and that is two lines shorter — but then P-127's C2, "the same polynomial
//! equals `disc(det(A0 + lambda*A1))` for all three axis pairings", would be
//! true by construction for one of the three pairings and the clause would be
//! measuring two facts while reporting three. The explicit form is the
//! independent witness; all three pencils are then genuine checks against it.
//!
//! The sign is not a free choice and was not chosen: with `Det` normalised as
//! `c1^2 - 4*c0*c2` (the plain discriminant of the pencil, no leading minus),
//! `repo_discriminant().sub(&cayley_2x2x2())` is **identically zero** — not
//! zero up to sign. The sympy script above is the independent cross-check that
//! produced the claim; this module is the one that gates it, because sympy is
//! not installed in `~/.venvs/isomesh` and is wired into no gate.

use std::collections::BTreeMap;
use std::fmt;

/// The number of variables. Eight, because a cell has eight corners; a
/// consumer using this for something else (R-137's `grad f = 0` system lives in
/// `x`, `y`, `z`) simply leaves the upper variables at exponent zero.
pub(crate) const VARS: usize = 8;

/// A multivariate polynomial in 8 variables `f0..f7` with exact `i128`
/// coefficients.
///
/// Monomial exponents are `[u8; 8]`; the map is ordered, so iteration,
/// [`Display`](fmt::Display) output and any CSV column derived from it are
/// deterministic. Zero coefficients are pruned on every write, so
/// [`Poly::terms`] is the count of genuinely non-zero terms and
/// [`Poly::is_zero`] is an emptiness test rather than a scan.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Poly {
    /// Exponent vector to coefficient. Never holds a zero coefficient.
    coeffs: BTreeMap<[u8; VARS], i128>,
}

/// Exact addition. An overflow here is a bug in the caller's magnitudes, not a
/// number to carry forward.
fn chk_add(a: i128, b: i128) -> i128 {
    a.checked_add(b)
        .expect("exact i128 addition overflowed; the inputs exceed this module's range")
}

/// Exact multiplication, with the same contract as [`chk_add`].
fn chk_mul(a: i128, b: i128) -> i128 {
    a.checked_mul(b)
        .expect("exact i128 multiplication overflowed; the inputs exceed this module's range")
}

/// `base^exponent` in checked `i128`.
fn ipow(base: i128, exponent: u32) -> i128 {
    let mut acc: i128 = 1;
    for _ in 0..exponent {
        acc = chk_mul(acc, base);
    }
    acc
}

/// Euclid's algorithm on magnitudes.
const fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Reduce `n/d` to lowest terms with a positive denominator. `0` reduces to
/// `(0, 1)` so that two zero results compare equal as pairs.
fn reduce(n: i128, d: i128) -> (i128, i128) {
    assert!(d != 0, "reduce: zero denominator");
    let (n, d) = if d < 0 {
        (
            n.checked_neg().expect("negating an exact numerator"),
            d.checked_neg().expect("negating an exact denominator"),
        )
    } else {
        (n, d)
    };
    if n == 0 {
        return (0, 1);
    }
    let g = gcd(n.unsigned_abs(), d.unsigned_abs());
    let g = i128::try_from(g).expect("a gcd of two i128 magnitudes fits in i128");
    (n / g, d / g)
}

impl Poly {
    /// The zero polynomial: no terms, [`Poly::total_degree`] zero.
    pub(crate) fn zero() -> Self {
        Self {
            coeffs: BTreeMap::new(),
        }
    }

    /// The constant polynomial `c`. `constant(0)` is [`Poly::zero`].
    pub(crate) fn constant(c: i128) -> Self {
        Self::monomial([0; VARS], c)
    }

    /// The variable `f_i`, for `i` in `0..8`.
    pub(crate) fn var(i: usize) -> Self {
        assert!(i < VARS, "variable index {i} is out of range 0..{VARS}");
        let mut exp = [0u8; VARS];
        exp[i] = 1;
        Self::monomial(exp, 1)
    }

    /// The single monomial `c * prod f_i^exp[i]`. A zero coefficient yields the
    /// zero polynomial rather than a stored zero term.
    pub(crate) fn monomial(exp: [u8; VARS], c: i128) -> Self {
        let mut coeffs = BTreeMap::new();
        if c != 0 {
            coeffs.insert(exp, c);
        }
        Self { coeffs }
    }

    /// Insert `c * exp` into the map, pruning the entry if it cancels.
    fn accumulate(&mut self, exp: [u8; VARS], c: i128) {
        if c == 0 {
            return;
        }
        match self.coeffs.get_mut(&exp) {
            Some(slot) => {
                let sum = chk_add(*slot, c);
                if sum == 0 {
                    self.coeffs.remove(&exp);
                } else {
                    *slot = sum;
                }
            }
            None => {
                self.coeffs.insert(exp, c);
            }
        }
    }

    /// `self + other`.
    pub(crate) fn add(&self, other: &Self) -> Self {
        let mut out = self.clone();
        for (exp, c) in &other.coeffs {
            out.accumulate(*exp, *c);
        }
        out
    }

    /// `self - other`. The zero-polynomial test on this is P-127's C1.
    pub(crate) fn sub(&self, other: &Self) -> Self {
        let mut out = self.clone();
        for (exp, c) in &other.coeffs {
            out.accumulate(
                *exp,
                c.checked_neg().expect("negating an exact coefficient"),
            );
        }
        out
    }

    /// `self * other`, expanded. Exponents are added with a checked `u8` add,
    /// so a degree above 255 in one variable panics rather than wrapping.
    pub(crate) fn mul(&self, other: &Self) -> Self {
        let mut out = Self::zero();
        for (ea, ca) in &self.coeffs {
            for (eb, cb) in &other.coeffs {
                let exp: [u8; VARS] = std::array::from_fn(|i| {
                    ea[i]
                        .checked_add(eb[i])
                        .expect("monomial exponent exceeded u8")
                });
                out.accumulate(exp, chk_mul(*ca, *cb));
            }
        }
        out
    }

    /// `k * self`. `scale(0)` is [`Poly::zero`].
    pub(crate) fn scale(&self, k: i128) -> Self {
        if k == 0 {
            return Self::zero();
        }
        Self {
            coeffs: self
                .coeffs
                .iter()
                .map(|(exp, c)| (*exp, chk_mul(*c, k)))
                .collect(),
        }
    }

    /// `self^n`, by repeated multiplication. `pow(0)` is the constant `1`.
    pub(crate) fn pow(&self, n: u32) -> Self {
        let mut acc = Self::constant(1);
        for _ in 0..n {
            acc = acc.mul(self);
        }
        acc
    }

    /// The number of non-zero terms. Twelve, for both sides of P-127's C1.
    pub(crate) fn terms(&self) -> usize {
        self.coeffs.len()
    }

    /// The largest total degree of any term; `0` for the zero polynomial.
    pub(crate) fn total_degree(&self) -> u32 {
        self.coeffs
            .keys()
            .map(|exp| exp.iter().map(|e| u32::from(*e)).sum::<u32>())
            .max()
            .unwrap_or(0)
    }

    /// The largest exponent of `f_i` in any term; `0` if `f_i` does not occur.
    pub(crate) fn degree_in(&self, i: usize) -> u32 {
        assert!(i < VARS, "variable index {i} is out of range 0..{VARS}");
        self.coeffs
            .keys()
            .map(|exp| u32::from(exp[i]))
            .max()
            .unwrap_or(0)
    }

    /// Degree at most one in every variable — the hypothesis every claim in
    /// Group A rests on, which R-135 records rather than assumes.
    pub(crate) fn is_multi_affine(&self) -> bool {
        self.coeffs.keys().all(|exp| exp.iter().all(|e| *e <= 1))
    }

    /// Whether every coefficient cancelled.
    pub(crate) fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// The coefficient of one monomial, `0` if absent.
    pub(crate) fn coefficient(&self, exp: [u8; VARS]) -> i128 {
        self.coeffs.get(&exp).copied().unwrap_or(0)
    }

    /// Every `(exponent, coefficient)` pair, ascending by exponent vector.
    /// Deterministic, which is why the map is a `BTreeMap`: a residual reported
    /// term-by-term (R-135's `symbolic_residual_terms`) must read the same on
    /// every run and every host.
    pub(crate) fn monomials(&self) -> impl Iterator<Item = ([u8; VARS], i128)> + '_ {
        self.coeffs.iter().map(|(exp, c)| (*exp, *c))
    }

    /// Exact evaluation at integer corner values.
    pub(crate) fn eval_i128(&self, f: &[i128; VARS]) -> i128 {
        let mut total: i128 = 0;
        for (exp, c) in &self.coeffs {
            let mut term = *c;
            for (v, e) in f.iter().zip(exp.iter()) {
                term = chk_mul(term, ipow(*v, u32::from(*e)));
            }
            total = chk_add(total, term);
        }
        total
    }

    /// Evaluation in `f64`, term by term over the *expanded* form.
    ///
    /// This is deliberately not the nested expression the crate evaluates. The
    /// gap between the two is the rounding P-127's C3 counts.
    pub(crate) fn eval_f64(&self, f: &[f64; VARS]) -> f64 {
        let mut total = 0.0_f64;
        for (exp, c) in &self.coeffs {
            let mut term = *c as f64;
            for (v, e) in f.iter().zip(exp.iter()) {
                for _ in 0..*e {
                    term *= *v;
                }
            }
            total += term;
        }
        total
    }

    /// Evaluation in `f32`, the precision `BodySaddles` actually runs in when
    /// the crate is instantiated at `f32`. Same expanded form as
    /// [`Poly::eval_f64`], so the two differ only in precision.
    pub(crate) fn eval_f32(&self, f: &[f32; VARS]) -> f32 {
        let mut total = 0.0_f32;
        for (exp, c) in &self.coeffs {
            let mut term = *c as f32;
            for (v, e) in f.iter().zip(exp.iter()) {
                for _ in 0..*e {
                    term *= *v;
                }
            }
            total += term;
        }
        total
    }

    /// Exact evaluation at the rationals `f_i = num[i] / den[i]`, returned as a
    /// reduced fraction with a positive denominator.
    ///
    /// The common denominator is `prod den[i]^degree_in(i)` — the smallest one
    /// this polynomial's own shape forces — so the fraction is exact before it
    /// is reduced and no division is ever inexact. Two calls whose reduced
    /// pairs are equal are equal as rationals, which is how C3's "ratio exactly
    /// 1" is checked without ever forming a ratio.
    ///
    /// # Panics
    ///
    /// If any `den[i]` is zero, or if the exact numerator overflows `i128`.
    pub(crate) fn eval_ratio(&self, num: &[i64; VARS], den: &[i64; VARS]) -> (i128, i128) {
        for (i, d) in den.iter().enumerate() {
            assert!(*d != 0, "eval_ratio: denominator {i} is zero");
        }
        let dmax: [u32; VARS] = std::array::from_fn(|i| self.degree_in(i));
        let mut common: i128 = 1;
        for (d, e) in den.iter().zip(dmax.iter()) {
            common = chk_mul(common, ipow(i128::from(*d), *e));
        }
        let mut total: i128 = 0;
        for (exp, c) in &self.coeffs {
            let mut term = *c;
            for i in 0..VARS {
                let e = u32::from(exp[i]);
                term = chk_mul(term, ipow(i128::from(num[i]), e));
                term = chk_mul(term, ipow(i128::from(den[i]), dmax[i] - e));
            }
            total = chk_add(total, term);
        }
        reduce(total, common)
    }

    /// `self` with `f_i` replaced by `g`, expanded.
    ///
    /// Powers of `g` are built once and reused, so substituting into a form of
    /// degree `d` in `f_i` costs `d` polynomial products rather than one per
    /// term. R-135 uses this to push a non-multi-affine reconstruction through
    /// the same construction and read the residual off.
    pub(crate) fn substitute(&self, i: usize, g: &Poly) -> Self {
        assert!(i < VARS, "variable index {i} is out of range 0..{VARS}");
        let mut powers = vec![Self::constant(1)];
        let mut out = Self::zero();
        for (exp, c) in &self.coeffs {
            let n = usize::from(exp[i]);
            while powers.len() <= n {
                let next = powers
                    .last()
                    .expect("the power table is seeded with g^0")
                    .mul(g);
                powers.push(next);
            }
            let mut rest = *exp;
            rest[i] = 0;
            out = out.add(&Self::monomial(rest, *c).mul(&powers[n]));
        }
        out
    }
}

impl fmt::Display for Poly {
    /// Ascending by exponent vector, `*` between factors, `^` for a power, and
    /// **no commas, quotes or newlines** — `common::experiment::Run::record`
    /// refuses all three, and P-127 records this string as its `expression`
    /// column.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.coeffs.is_empty() {
            return f.write_str("0");
        }
        let mut first = true;
        for (exp, c) in &self.coeffs {
            if first {
                if *c < 0 {
                    f.write_str("-")?;
                }
            } else {
                f.write_str(if *c < 0 { " - " } else { " + " })?;
            }
            first = false;

            let magnitude = c.unsigned_abs();
            let bare = exp.iter().all(|e| *e == 0);
            if magnitude != 1 || bare {
                write!(f, "{magnitude}")?;
                if !bare {
                    f.write_str("*")?;
                }
            }
            let mut between = false;
            for (i, e) in exp.iter().enumerate() {
                if *e == 0 {
                    continue;
                }
                if between {
                    f.write_str("*")?;
                }
                between = true;
                if *e == 1 {
                    write!(f, "f{i}")?;
                } else {
                    write!(f, "f{i}^{e}")?;
                }
            }
        }
        Ok(())
    }
}

/// The monomial `coefficient * prod f[indices[k]]`, with repeats becoming
/// powers. Writing Cayley's form as corner-index quadruples keeps it readable
/// beside the `a(i, j, k)` notation of the sympy script.
fn corner_product(indices: [usize; 4], coefficient: i128) -> Poly {
    let mut exp = [0u8; VARS];
    for i in indices {
        exp[i] += 1;
    }
    Poly::monomial(exp, coefficient)
}

/// Cayley's `2x2x2` hyperdeterminant of the eight corner values, in the crate's
/// corner indexing `f[u + 2v + 4w]`.
///
/// Twelve terms, total degree 4: four squares `a_ijk^2 * a_(1-i)(1-j)(1-k)^2`,
/// six antipodal-pair products at `-2`, and two "even/odd tetrahedron" products
/// at `+4`. Under `f[u + 2v + 4w]` the tensor entry `a(i, j, k)` is `f[i + 2j +
/// 4k]`, so `a(0,0,0) = f0`, `a(0,0,1) = f4`, `a(1,1,1) = f7`, and so on.
///
/// The normalisation is the pencil discriminant's, `c1^2 - 4*c0*c2`, with no
/// leading sign — see the module header for why that is a fact rather than a
/// convention here.
pub(crate) fn cayley_2x2x2() -> Poly {
    // a(0,0,0)^2 a(1,1,1)^2 + a(0,0,1)^2 a(1,1,0)^2
    //   + a(0,1,0)^2 a(1,0,1)^2 + a(1,0,0)^2 a(0,1,1)^2
    let squares = [[0, 0, 7, 7], [4, 4, 3, 3], [2, 2, 5, 5], [1, 1, 6, 6]];
    // The six products of two antipodal corner pairs.
    let pairs = [
        [0, 4, 3, 7],
        [0, 2, 5, 7],
        [0, 1, 6, 7],
        [4, 2, 5, 3],
        [4, 1, 6, 3],
        [2, 1, 6, 5],
    ];
    // The two products of one whole inscribed tetrahedron.
    let tetrahedra = [[0, 6, 5, 3], [4, 2, 1, 7]];

    let mut out = Poly::zero();
    for indices in squares {
        out = out.add(&corner_product(indices, 1));
    }
    for indices in pairs {
        out = out.add(&corner_product(indices, -2));
    }
    for indices in tetrahedra {
        out = out.add(&corner_product(indices, 4));
    }
    out
}

/// `b*b - 4*a*c` built symbolically from `BodySaddles::coefficients`.
///
/// Transcribed line for line from
/// `crates/isomesh/src/marching_cubes/trilinear.rs:202-213`, with the
/// discriminant from :246. Nothing is simplified on the way in: the two twists,
/// the four edge differences and the three coefficients are formed exactly as
/// the crate forms them, so a future edit to `coefficients` that changes the
/// polynomial will show up here as a non-zero difference rather than as a
/// comment that has quietly stopped being true.
pub(crate) fn repo_discriminant() -> Poly {
    let f: [Poly; VARS] = std::array::from_fn(Poly::var);

    let twist_lo = f[0].add(&f[3]).sub(&f[1].add(&f[2]));
    let twist_hi = f[4].add(&f[7]).sub(&f[5].add(&f[6]));
    let du_lo = f[1].sub(&f[0]);
    let du_hi = f[5].sub(&f[4]);
    let dv_lo = f[2].sub(&f[0]);
    let dv_hi = f[6].sub(&f[4]);

    let a = du_hi.mul(&twist_lo).sub(&du_lo.mul(&twist_hi));
    // `i0 = 0`, so `(i0 - f0)` is `-f0` and `(i0 - f4)` is `-f4`.
    let b = f[4]
        .mul(&twist_lo)
        .sub(&f[0].mul(&twist_hi))
        .add(&du_hi.mul(&dv_lo).sub(&du_lo.mul(&dv_hi)));
    let c = f[2].mul(&f[4]).sub(&f[0].mul(&f[6]));

    b.mul(&b).sub(&a.mul(&c).scale(4))
}

/// The two opposite-face corner index sets of one axis pairing, each read
/// row-major into a `2x2` matrix.
///
/// `0 => 0123|4567` splits along `w`, `1 => 0145|2367` along `v`, and
/// `2 => 0246|1357` along `u`. Each set lists the four corners of one face in
/// the order the remaining two bits count up.
const PAIRINGS: [([usize; 4], [usize; 4]); 3] = [
    ([0, 1, 2, 3], [4, 5, 6, 7]),
    ([0, 1, 4, 5], [2, 3, 6, 7]),
    ([0, 2, 4, 6], [1, 3, 5, 7]),
];

/// The human-readable name of each pairing, for the `pencil_axis_pairing`
/// column.
pub(crate) const PAIRING_NAMES: [&str; 3] = [
    "w-slices 0123|4567",
    "v-slices 0145|2367",
    "u-slices 0246|1357",
];

/// `disc(det(A0 + lambda*A1))` for one of the three axis pairings.
///
/// `det(A0 + lambda*A1)` is a quadratic `c2*lambda^2 + c1*lambda + c0` in the
/// pencil parameter; its discriminant is `c1^2 - 4*c0*c2`. That the three
/// pairings give the *same* polynomial is not obvious from the construction —
/// it is the `SL(2)^3` invariance, and it is the mechanism behind `M-206`.
///
/// # Panics
///
/// If `pairing` is not in `0..3`.
pub(crate) fn pencil_discriminant(pairing: usize) -> Poly {
    assert!(pairing < 3, "axis pairing {pairing} is out of range 0..3");
    let (s0, s1) = PAIRINGS[pairing];
    let m: [Poly; 4] = std::array::from_fn(|k| Poly::var(s0[k]));
    let n: [Poly; 4] = std::array::from_fn(|k| Poly::var(s1[k]));

    // det([[m0, m1], [m2, m3]] + lambda * [[n0, n1], [n2, n3]]).
    let c0 = m[0].mul(&m[3]).sub(&m[1].mul(&m[2]));
    let c1 = m[0]
        .mul(&n[3])
        .add(&n[0].mul(&m[3]))
        .sub(&m[1].mul(&n[2]))
        .sub(&n[1].mul(&m[2]));
    let c2 = n[0].mul(&n[3]).sub(&n[1].mul(&n[2]));

    c1.mul(&c1).sub(&c0.mul(&c2).scale(4))
}

/// The 48 octahedral relabellings of the eight corners, as permutations of
/// `0..8`.
///
/// Generated, never transcribed: the octahedral group in its full form is the
/// `3! = 6` axis permutations times the `2^3 = 8` axis flips, i.e. the signed
/// permutation matrices of the cube. Each `(pi, flips)` pair relabels a corner
/// by permuting its three coordinate bits and then flipping them, and the
/// corner index is rebuilt as `u + 2v + 4w`. The result is applied with
/// [`relabel`], which reads `out[i] = f[perm[i]]`.
///
/// # Panics
///
/// If the 48 are not distinct, or are not closed under composition. Both are
/// plain `assert!`s and not `debug_assert!`s: benches run in the release
/// profile, and a consumer that silently received a 47-element group with a
/// duplicate would report an invariance result about the wrong set.
pub(crate) fn octahedral_relabellings() -> [[u8; VARS]; 48] {
    // The six permutations of the three axes, enumerated rather than listed.
    let mut axis_perms = [[0u8; 3]; 6];
    let mut found = 0;
    for a in 0..3u8 {
        for b in 0..3u8 {
            for c in 0..3u8 {
                if a != b && b != c && a != c {
                    axis_perms[found] = [a, b, c];
                    found += 1;
                }
            }
        }
    }
    assert!(found == 6, "there are exactly 3! axis permutations");

    let mut out = [[0u8; VARS]; 48];
    let mut slot = 0;
    for pi in axis_perms {
        for flips in 0..8u8 {
            let mut perm = [0u8; VARS];
            for idx in 0..8u8 {
                let bit = |t: u8| (idx >> t) & 1;
                let mut moved = 0u8;
                for t in 0..3u8 {
                    let b = bit(pi[usize::from(t)]) ^ ((flips >> t) & 1);
                    moved |= b << t;
                }
                perm[usize::from(moved)] = idx;
            }
            out[slot] = perm;
            slot += 1;
        }
    }
    assert!(slot == 48, "6 axis permutations times 8 flips is 48");

    for i in 0..48 {
        for j in (i + 1)..48 {
            assert!(
                out[i] != out[j],
                "octahedral relabellings {i} and {j} coincide; the table is not a 48-element set"
            );
        }
    }
    for p in &out {
        for q in &out {
            let composed: [u8; VARS] = std::array::from_fn(|i| q[usize::from(p[i])]);
            assert!(
                out.contains(&composed),
                "the octahedral relabellings are not closed under composition"
            );
        }
    }
    out
}

/// Apply a relabelling to a corner tuple: `out[i] = f[perm[i]]`.
///
/// The same convention the sympy script applies its permutations with, so an
/// index into [`octahedral_relabellings`] names the same relabelling in both.
pub(crate) fn relabel<T: Copy>(perm: &[u8; VARS], f: &[T; VARS]) -> [T; VARS] {
    std::array::from_fn(|i| f[usize::from(perm[i])])
}

/// The `GL(2)^3` action on the `2x2x2` tensor: `g1` on the `u` index, `g2` on
/// `v`, `g3` on `w`. Each `g` is `[[a, b], [c, d]]` row-major. Exact over
/// `i128`.
///
/// `out[i + 2j + 4k] = sum_{p,q,r} g1[i][p] * g2[j][q] * g3[k][r] * f[p + 2q +
/// 4r]`. Cayley's `Det` is a *relative* invariant of this action with weight
/// `(det g1 * det g2 * det g3)^2` — a perfect square, which is exactly why the
/// **sign** of the discriminant, and hence the body-saddle count, cannot depend
/// on a per-axis affine reparametrisation of the cell. That is R-128's claim,
/// and it is checked with this function rather than argued for.
pub(crate) fn act_gl2_cubed(
    g1: [[i128; 2]; 2],
    g2: [[i128; 2]; 2],
    g3: [[i128; 2]; 2],
    f: &[i128; VARS],
) -> [i128; VARS] {
    let mut out = [0i128; VARS];
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..2 {
                let mut acc: i128 = 0;
                for r in 0..2 {
                    for q in 0..2 {
                        for p in 0..2 {
                            let w = chk_mul(chk_mul(g1[i][p], g2[j][q]), g3[k][r]);
                            acc = chk_add(acc, chk_mul(w, f[p + 2 * q + 4 * r]));
                        }
                    }
                }
                out[i + 2 * j + 4 * k] = acc;
            }
        }
    }
    out
}

/// `f64` version of [`act_gl2_cubed`], for the numeric weight check.
///
/// Same index arithmetic, no checked operations: overflow in `f64` is an
/// infinity that propagates visibly, and the point of this arm is to measure
/// how the exact weight law degrades in floating point rather than to be exact.
pub(crate) fn act_gl2_cubed_f64(
    g1: [[f64; 2]; 2],
    g2: [[f64; 2]; 2],
    g3: [[f64; 2]; 2],
    f: &[f64; VARS],
) -> [f64; VARS] {
    let mut out = [0.0f64; VARS];
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..2 {
                let mut acc = 0.0f64;
                for r in 0..2 {
                    for q in 0..2 {
                        for p in 0..2 {
                            acc += g1[i][p] * g2[j][q] * g3[k][r] * f[p + 2 * q + 4 * r];
                        }
                    }
                }
                out[i + 2 * j + 4 * k] = acc;
            }
        }
    }
    out
}

/// A deterministic SplitMix64, so every bench that draws random 8-tuples draws
/// the same ones.
///
/// Vigna's mixing function on a 64-bit Weyl sequence: ten lines, no dependency,
/// and the same stream on every host and every re-run. That last property is
/// the reason it exists — P-127's C3 reports a count of `f32` sign
/// disagreements over 3,000 tuples, and a count that changes between runs is
/// not a measurement.
#[derive(Clone, Debug)]
pub(crate) struct Rng {
    /// The Weyl-sequence state; advanced by the golden-ratio odd constant.
    state: u64,
}

impl Rng {
    /// Seed the generator. Every seed is valid, including zero.
    pub(crate) fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// The next 64 bits of the stream.
    pub(crate) fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A value in `[lo, hi)` — inclusive `lo`, exclusive `hi`.
    ///
    /// Lemire's multiply-shift rather than a modulo: the high half of a
    /// `u128` product spreads the bias across the whole range instead of
    /// concentrating it in the first `2^64 mod span` values, which matters when
    /// the range is tiny (`randint(1, 5)` for a denominator) and the bias would
    /// otherwise land on one specific denominator.
    ///
    /// # Panics
    ///
    /// If `hi <= lo`.
    pub(crate) fn next_i64_in(&mut self, lo: i64, hi: i64) -> i64 {
        assert!(hi > lo, "next_i64_in: empty range [{lo}, {hi})");
        let span = (i128::from(hi) - i128::from(lo)) as u128;
        let scaled = (u128::from(self.next_u64()) * span) >> 64;
        let offset = i128::try_from(scaled).expect("a value below span fits i128");
        i64::try_from(i128::from(lo) + offset).expect("lo + offset is inside [lo, hi)")
    }

    /// A value in `[-1, 1)`, from the top 53 bits so every representable
    /// multiple of `2^-52` in the range is reachable and none is favoured.
    pub(crate) fn next_f64_unit(&mut self) -> f64 {
        let unit = (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0);
        2.0f64.mul_add(unit, -1.0)
    }
}
