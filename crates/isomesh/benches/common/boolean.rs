//! Analysis of Boolean functions, applied to the 256-case Marching Cubes table.
//!
//! Ticket: R-167, which owns this module. Consumed unchanged by R-168 (noise
//! stability against the wrong noise model) and R-169 (per-corner influence as
//! an independent check on `validate_table()`).
//!
//! # The object
//!
//! A Marching Cubes cell is eight sign bits and a lookup. That makes the table a
//! map `{0,1}^8 -> N`, and once an output integer is fixed, each of its bits is
//! an honest Boolean function `f: {0,1}^8 -> {0,1}` — exactly the object
//! O'Donnell's *Analysis of Boolean Functions* (`arXiv:2105.10386`) is about.
//! The whole apparatus below is that book's Chapter 1 and Chapter 2, at `n = 8`.
//!
//! Two encodings are in play and keeping them apart is the only subtlety here:
//!
//! - the **ANF convention**, `f: {0,1}^8 -> {0,1}`, in which the natural algebra
//!   is `GF(2)` and the natural transform is Moebius inversion, giving a
//!   multilinear polynomial over `GF(2)` — the algebraic normal form;
//! - the **Fourier convention**, `chi = 1 - 2f: {0,1}^8 -> {-1,+1}`, in which the
//!   natural algebra is `R` and the natural transform is Walsh-Hadamard, giving a
//!   multilinear polynomial over the reals whose coefficients are orthonormal.
//!
//! Storage is the `0/1` table; `fourier()` converts on the way in. The two
//! transforms answer different questions and neither substitutes for the other:
//! the ANF term count is a *sparsity* statement about a branchless evaluation,
//! and the Fourier weight distribution is a *degree* statement about how much of
//! the function a low-order approximation could ever capture.
//!
//! # What the transform buys, concretely
//!
//! Three facts do all the work and each one is checked here rather than assumed:
//!
//! 1. **Parseval.** For a `+/-1`-valued function the Fourier coefficients are the
//!    coordinates of a unit vector: `sum_S fhat(S)^2 = E[chi^2] = 1`. A transform
//!    that has lost its `1/256` normalisation still looks plausible — every
//!    weight is merely scaled — so every clause in Group J would silently read a
//!    wrong concentration. [`Bool8::fourier`] asserts Parseval before returning.
//! 2. **Influence is a spectral quantity.** `Inf_i(f) = Pr[f(x) != f(x + e_i)]`
//!    counted combinatorially equals `sum_{S ni i} fhat(S)^2` read off the
//!    spectrum. Two independent computations of the same number, so
//!    [`Bool8::influence`] runs both and asserts they agree — a total check on the
//!    transform that costs nothing.
//! 3. **Noise stability is the spectrum, reweighted.** `Stab_rho(f) = sum_S
//!    rho^|S| fhat(S)^2`, which is why R-168 can predict a flip rate without ever
//!    sampling noise: `Pr[flip] = (1 - Stab_rho)/2`.
//!
//! # Instrument validation
//!
//! [`self_check`] runs two functions whose spectra are known in closed form
//! through the same code path the table uses. Parity is `chi(x) = (-1)^{|x|}`,
//! which *is* the character `chi_{[8]}`, so its entire weight sits at degree 8 —
//! and that is asserted, not merely reported. Majority-of-8 with ties broken to
//! `1` has `fhat({i}) = 70/256` for every `i` and `fhat(empty) = -70/256`, giving
//! `W^1 = 8 * (70/256)^2 = 0.598...` and `Inf_i = C(7,3)/2^7 = 35/128`, because
//! flipping one input changes the majority exactly when the other seven split
//! `3-4`. Total influence is therefore `8 * 35/128 = 2.1875` exactly.

use isomesh::marching_cubes::table::{CASES, CENTROID_BASE, EDGE_CORNERS, corner_inside};

/// Below this magnitude a Fourier coefficient is treated as absent.
///
/// The coefficients of a `+/-1` table on 256 points are exact multiples of
/// `1/256` — integer sums of at most 256 terms of `+/-1`, all representable — so
/// a genuine coefficient is never smaller than `1/256` and a genuine zero is
/// never larger than a few multiples of `f64::EPSILON`. Any threshold between
/// those two extremes gives the same answer; `1e-12` is stated once so degree,
/// concentration and the Parseval assertion cannot drift apart.
pub(crate) const NEGLIGIBLE: f64 = 1e-12;

/// A Boolean function on 8 variables, as a truth table of 256 entries valued
/// `-1`/`+1` in the Fourier convention, or `0`/`1` in the ANF convention. Stored
/// as the `0/1` table.
///
/// Input `x` is the index: bit `i` of the index is the value of variable `i`.
/// For the case table that means bit `i` is the sign of corner `i`, matching
/// `corner_inside(case, corner) == case & (1 << corner) != 0`, so a `Bool8` built
/// from a case table has its variables numbered the way the cube numbers its
/// corners — corner `i` at local coordinate `(i & 1, (i >> 1) & 1, (i >> 2) & 1)`.
///
/// Every entry must be `0` or `1`. The field is visible so a consumer can read a
/// truth table out, not so it can put arbitrary integers in; [`Bool8::fourier`]'s
/// Parseval assertion is what catches a violation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Bool8(pub(crate) [u8; 256]);

impl Bool8 {
    /// Bit `bit` of the shipped table's output, as a Boolean function of the
    /// eight corner signs. `values` supplies the per-case output integer.
    ///
    /// This is the per-bit treatment P-167's C1 falsifier names. A table whose
    /// output is an integer is a *vector*-valued Boolean function, and the
    /// Walsh-Hadamard transform is defined for a scalar one; the standard
    /// resolution is to analyse each output bit separately and report the
    /// spectrum per bit, which is what `output_bit` is a column for.
    ///
    /// # Panics
    ///
    /// If `bit >= 32`, which would shift a `u32` out of range.
    pub(crate) fn from_values(values: &[u32; 256], bit: u32) -> Self {
        assert!(bit < 32, "output bit {bit} does not exist in a u32");
        let mut table = [0u8; 256];
        for (slot, &v) in table.iter_mut().zip(values.iter()) {
            *slot = ((v >> bit) & 1) as u8;
        }
        Self(table)
    }

    /// Parity of all eight inputs — a known degree-8 function, for instrument
    /// validation.
    ///
    /// In the `+/-1` encoding this is literally the character `chi_{[8]}(x) =
    /// (-1)^{x_0 + ... + x_7}`, so its transform is the indicator of the full set
    /// and nothing else. It is the sharpest possible calibration: any lost
    /// normalisation, any transposed butterfly stride and any sign error moves
    /// weight off degree 8, where it is trivially visible.
    pub(crate) fn parity() -> Self {
        let mut table = [0u8; 256];
        for (x, slot) in table.iter_mut().enumerate() {
            *slot = (x.count_ones() & 1) as u8;
        }
        Self(table)
    }

    /// Majority of the eight inputs, ties broken to `1` — a known
    /// low-degree-heavy function.
    ///
    /// Eight is even, so `Maj_8` is not a total function without a tie rule;
    /// breaking to `1` makes it `[popcount >= 4]`. The cost of the convention is
    /// that the function is no longer odd, so it carries a degree-0 term
    /// (`fhat(empty) = -70/256`, a bias of `70/256` towards `1`) and even-degree
    /// terms generally. What survives is the point of the calibration: nearly
    /// 60% of its weight sits at degree 1, against parity's 0%, so a transform
    /// that cannot tell the two apart is broken.
    pub(crate) fn majority() -> Self {
        let mut table = [0u8; 256];
        for (x, slot) in table.iter_mut().enumerate() {
            *slot = u8::from(x.count_ones() >= 4);
        }
        Self(table)
    }

    /// The Walsh-Hadamard transform: 256 Fourier coefficients `fhat(S)`, indexed
    /// by the subset `S` as an 8-bit mask. Uses the `+/-1` convention, so the sum
    /// of squares is 1.
    ///
    /// The in-place butterfly, `O(n log n)` — eight passes of 128 add/subtract
    /// pairs, 2048 flops against the 65,536 of the defining double sum. Stride
    /// `len` doubles from 1 to 128 and pass `k` is exactly the one-variable
    /// transform in variable `k`, which is why the output index is the subset
    /// mask with no bit reversal: the transform is a tensor product of eight
    /// identical `2x2` blocks and the natural ordering is already the right one.
    ///
    /// # Panics
    ///
    /// If Parseval fails to `1e-12`, which means either the normalisation is
    /// wrong or the stored table was not `0/1`-valued.
    pub(crate) fn fourier(&self) -> [f64; 256] {
        let mut a = [0.0f64; 256];
        for (slot, &b) in a.iter_mut().zip(self.0.iter()) {
            *slot = 1.0 - 2.0 * f64::from(b);
        }

        let mut len = 1usize;
        while len < 256 {
            let mut base = 0usize;
            while base < 256 {
                for j in base..base + len {
                    let u = a[j];
                    let v = a[j + len];
                    a[j] = u + v;
                    a[j + len] = u - v;
                }
                base += len << 1;
            }
            len <<= 1;
        }

        for c in &mut a {
            *c /= 256.0;
        }

        let parseval: f64 = a.iter().map(|c| c * c).sum();
        assert!(
            (parseval - 1.0).abs() <= NEGLIGIBLE,
            "Walsh-Hadamard lost Parseval: sum of squares is {parseval}, not 1 — \
             the transform is unnormalised or the truth table is not 0/1-valued"
        );
        a
    }

    /// Fourier weight by degree: `w[k]` is the sum of `fhat(S)^2` over `|S| = k`.
    /// Sums to 1.
    ///
    /// This is the whole answer to "is the table low-degree": `w` is a probability
    /// distribution over `{0..8}` and the question is where its mass is. Parity
    /// puts all of it on 8; majority puts 60% on 1.
    pub(crate) fn weight_by_degree(&self) -> [f64; 9] {
        let f = self.fourier();
        let mut w = [0.0f64; 9];
        for (s, c) in f.iter().enumerate() {
            w[s.count_ones() as usize] += c * c;
        }
        w
    }

    /// The largest `k` with non-negligible weight — the function's degree.
    ///
    /// Read off the individual coefficients rather than the per-degree weights,
    /// because a single surviving coefficient at degree 8 is a degree-8 function
    /// however small its weight, and summing first could hide it under a
    /// threshold. A constant function has degree 0.
    pub(crate) fn max_degree(&self) -> u32 {
        let f = self.fourier();
        let mut deg = 0u32;
        for (s, c) in f.iter().enumerate() {
            if c.abs() > NEGLIGIBLE {
                deg = deg.max(s.count_ones());
            }
        }
        deg
    }

    /// The fraction of Fourier weight on degrees `<= k`.
    ///
    /// `k >= 8` returns 1 by Parseval. This is the number P-167's C2 is about:
    /// "concentrated on low degrees" means this is close to 1 for small `k`.
    pub(crate) fn concentration_up_to(&self, k: u32) -> f64 {
        let w = self.weight_by_degree();
        let top = (k as usize).min(8);
        w.iter().take(top + 1).sum()
    }

    /// Algebraic normal form over `GF(2)`: the 256 ANF coefficients, and the
    /// count of non-zero ones. Computed by the Moebius transform.
    ///
    /// The ANF is the unique multilinear polynomial over `GF(2)` agreeing with
    /// the function, `f(x) = sum_S a_S prod_{i in S} x_i`, and `a_S = sum_{T
    /// subset S} f(T)` — the `GF(2)` Moebius inversion, which is the same
    /// butterfly as Walsh-Hadamard with `xor` in place of `+/-`. The non-zero
    /// count is the term count of a branchless evaluation: `sparse` in P-167's
    /// record list, and the quantitative form of its negative prediction.
    ///
    /// Note that ANF degree and Fourier degree are different numbers. An edge
    /// occupancy bit `x_a xor x_b` has ANF degree 1 and two terms, but Fourier
    /// degree 2, because `chi_a chi_b` is its `+/-1` form.
    pub(crate) fn anf(&self) -> ([u8; 256], usize) {
        let mut a = self.0;
        for i in 0..8u32 {
            let bit = 1usize << i;
            for x in 0..256usize {
                if x & bit != 0 {
                    let lower = a[x ^ bit];
                    a[x] ^= lower;
                }
            }
        }
        let terms = a.iter().filter(|&&c| c != 0).count();
        (a, terms)
    }

    /// Influence of variable `i`: the probability a uniformly random input has
    /// its output changed by flipping bit `i`. Equals the sum of `fhat(S)^2` over
    /// `S` containing `i`.
    ///
    /// Computed **both** ways and asserted to agree. The identity is exact:
    /// the discrete derivative `D_i f` has Fourier coefficients `fhat(S)` for
    /// `S ni i` and zero elsewhere, and `Pr[flip]` is its second moment. So the
    /// combinatorial count and the spectral sum are the same rational number,
    /// and disagreement can only mean the transform is wrong.
    ///
    /// # Panics
    ///
    /// If `i >= 8`, or if the two computations disagree by more than `1e-12`.
    pub(crate) fn influence(&self, i: usize) -> f64 {
        assert!(
            i < 8,
            "variable {i} does not exist in an 8-variable function"
        );
        let f = self.fourier();
        let bit = 1usize << i;
        let mut spectral = 0.0f64;
        for (s, c) in f.iter().enumerate() {
            if s & bit != 0 {
                spectral += c * c;
            }
        }
        let combinatorial = self.influence_combinatorial(i);
        assert!(
            (spectral - combinatorial).abs() <= NEGLIGIBLE,
            "influence of variable {i} disagrees between the spectrum \
             ({spectral}) and the flip count ({combinatorial}) — the \
             Walsh-Hadamard transform is wrong"
        );
        spectral
    }

    /// Influence of variable `i`, counted directly: the fraction of the 256
    /// inputs whose output changes when bit `i` is flipped.
    ///
    /// The naive definition, kept as an independent instrument rather than an
    /// optimisation — 256 comparisons against a 2048-flop transform is not a
    /// speed question, it is the second opinion [`Bool8::influence`] checks
    /// against.
    ///
    /// # Panics
    ///
    /// If `i >= 8`.
    pub(crate) fn influence_combinatorial(&self, i: usize) -> f64 {
        assert!(
            i < 8,
            "variable {i} does not exist in an 8-variable function"
        );
        let bit = 1usize << i;
        let mut flips = 0u32;
        for x in 0..256usize {
            if self.0[x] != self.0[x ^ bit] {
                flips += 1;
            }
        }
        f64::from(flips) / 256.0
    }

    /// Total influence, the sum over all eight variables. Equals
    /// `sum_S |S| fhat(S)^2`.
    ///
    /// Also known as average sensitivity: the expected number of pivotal
    /// coordinates at a random input. For the table it answers "how many of the
    /// eight corners is a random cell's output actually resting on", which is
    /// the quantity a refinement heuristic would be trying to spend on.
    ///
    /// # Panics
    ///
    /// If the spectral form and the sum of the eight combinatorial influences
    /// disagree by more than `1e-12`.
    pub(crate) fn total_influence(&self) -> f64 {
        let f = self.fourier();
        let spectral: f64 = f
            .iter()
            .enumerate()
            .map(|(s, c)| f64::from(s.count_ones()) * c * c)
            .sum();
        let combinatorial: f64 = (0..8).map(|i| self.influence_combinatorial(i)).sum();
        assert!(
            (spectral - combinatorial).abs() <= NEGLIGIBLE,
            "total influence disagrees between sum |S| fhat(S)^2 ({spectral}) \
             and the eight flip counts ({combinatorial})"
        );
        spectral
    }

    /// Noise stability at correlation `rho`: `sum_S rho^|S| fhat(S)^2`.
    ///
    /// `Stab_rho(f) = E[chi(x) chi(y)]` where `y` is `x` with each bit
    /// independently re-randomised so that `Pr[y_i = x_i] = (1 + rho)/2`. For a
    /// `+/-1`-valued function the product is `+/-1`, so
    /// `Pr[chi(y) != chi(x)] = (1 - Stab_rho)/2` — the predicted flip rate R-168
    /// compares against `f32` rounding, and the comparison is deliberately
    /// against the wrong noise model: rounding is deterministic and spatially
    /// correlated, not eight independent coin flips.
    ///
    /// `rho = 1` gives 1 by Parseval, `rho = 0` gives `fhat(empty)^2`, and
    /// negative `rho` is meaningful (anti-correlated inputs) since only integer
    /// powers up to 8 are taken.
    pub(crate) fn noise_stability(&self, rho: f64) -> f64 {
        let mut pow = [1.0f64; 9];
        for k in 1..9usize {
            pow[k] = pow[k - 1] * rho;
        }
        let f = self.fourier();
        f.iter()
            .enumerate()
            .map(|(s, c)| pow[s.count_ones() as usize] * c * c)
            .sum()
    }
}

/// Validate the instrument on two functions with known spectra, and return
/// (parity weight at degree 8, majority weight at degree 1, majority total
/// influence) so a bench can record the calibration as columns.
///
/// Parity's numbers are *asserted*, because they are exact and structural: the
/// transform of `chi_{[8]}` is the indicator of `[8]`, so degree-8 weight is 1,
/// degree is 8, and every influence is 1 (flipping any input always flips the
/// parity), giving total influence 8. Majority's numbers are returned rather than
/// asserted, so that the bench rather than the module owns the comparison — but
/// they are equally exact: `W^1 = 8 * (70/256)^2 = 0.598144531250` and total
/// influence `= 8 * C(7,3)/2^7 = 2.1875`.
///
/// This is P-167's registered vacuity control: "a known low-degree function
/// (parity, majority) must be run through the same transform and reproduce its
/// known spectrum, or the instrument is unvalidated."
///
/// # Panics
///
/// If parity's spectrum is not the single degree-8 coefficient.
pub(crate) fn self_check() -> (f64, f64, f64) {
    let parity = Bool8::parity();
    let parity_weight = parity.weight_by_degree();
    assert!(
        (parity_weight[8] - 1.0).abs() <= NEGLIGIBLE,
        "VOID: parity does not have all its Fourier weight at degree 8 \
         (got {}), so the Walsh-Hadamard transform is not measuring what it \
         claims and no Group J number means anything",
        parity_weight[8]
    );
    assert_eq!(parity.max_degree(), 8, "VOID: parity is not degree 8");
    assert!(
        (parity.total_influence() - 8.0).abs() <= NEGLIGIBLE,
        "VOID: parity's total influence is not 8"
    );

    let majority = Bool8::majority();
    let majority_weight = majority.weight_by_degree();
    (
        parity_weight[8],
        majority_weight[1],
        majority.total_influence(),
    )
}

/// The eight corners' octahedral symmetry classes under the 48-element cube
/// group. Two corners are in the same class iff some cube symmetry maps one to
/// the other — which for the eight corners of a cube is ALL of them, one class.
/// Returns the class index per corner, and the class count. Generated, not
/// asserted.
///
/// # How the 48 are built
///
/// The full octahedral group `O_h` is the signed permutations of the three axes:
/// six orderings times eight sign patterns. Acting on the unit cube it is
/// affine — conjugate by the translation that puts the centre at the origin — so
/// on corner coordinates `c in {0,1}^3` a symmetry `(p, s)` sends
/// `c'[a] = if s[a] { 1 - c[p[a]] } else { c[p[a]] }`, and the corner index is
/// `index = u + 2v + 4w`, i.e. bit `a` of the index is coordinate `a`. Each of
/// the 48 is therefore a permutation of `0..8`, and this function computes the
/// orbits of that permutation group by flood fill.
///
/// The answer is one orbit — the cube group is transitive on corners, since the
/// three sign flips alone already carry corner 0 to all eight — so "equal within
/// each octahedral symmetry class" in P-169's C1 is the strongest reading it
/// could have: **all eight influences must be equal to each other**. That is
/// worth generating rather than assuming, because the clause's meaning depends
/// entirely on it: had the corners split into two classes, C1 would only have
/// constrained four-against-four and a real table defect could have hidden in
/// the gap.
///
/// # Panics
///
/// If the 48 constructions are not 48 distinct permutations, which would mean the
/// generator is degenerate and the orbit count is not the cube group's.
pub(crate) fn corner_symmetry_classes() -> ([u8; 8], usize) {
    let maps = cube_symmetries();

    let mut distinct: Vec<[u8; 8]> = maps.to_vec();
    distinct.sort_unstable();
    distinct.dedup();
    assert_eq!(
        distinct.len(),
        48,
        "the cube symmetry generator produced {} distinct corner permutations, \
         not 48 — the octahedral group has order 48",
        distinct.len()
    );

    let mut class = [u8::MAX; 8];
    let mut count = 0usize;
    for start in 0..8u8 {
        if class[start as usize] != u8::MAX {
            continue;
        }
        let id = count as u8;
        count += 1;
        class[start as usize] = id;
        let mut frontier = vec![start];
        while let Some(c) = frontier.pop() {
            for m in &maps {
                let d = m[c as usize];
                if class[d as usize] == u8::MAX {
                    class[d as usize] = id;
                    frontier.push(d);
                }
            }
        }
    }
    (class, count)
}

/// The 48 elements of the cube group, each as the permutation it induces on the
/// eight corner indices.
///
/// Enumerated in a fixed order — the six axis orderings outermost, the eight sign
/// patterns innermost — so the output is byte-identical on every run and on every
/// platform. Nothing downstream depends on the order, but a shared module whose
/// output can permute is a shared module that can make two benches disagree.
fn cube_symmetries() -> [[u8; 8]; 48] {
    const PERMS: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];

    let mut out = [[0u8; 8]; 48];
    let mut g = 0usize;
    for perm in PERMS {
        for signs in 0..8u8 {
            for corner in 0..8u8 {
                let mut image = 0u8;
                for (axis, &source) in perm.iter().enumerate() {
                    let mut bit = (corner >> source) & 1;
                    if (signs >> axis) & 1 == 1 {
                        bit ^= 1;
                    }
                    image |= bit << axis;
                }
                out[g][corner as usize] = image;
            }
            g += 1;
        }
    }
    out
}

/// The shipped 256-case table's triangle count per case, read from
/// `isomesh::marching_cubes::table::CASES`.
///
/// **The recommended primary reading, for both P-167 and P-169.** Two properties
/// pick it out, and only it:
///
/// - It is where the table's triangulation *decisions* live. `centroids` is
///   identically zero across all 256 cases of the plain table (the centroid path
///   is only reached by the joined pairing of the asymptotic decider), and
///   `edge_masks` is by construction twelve independent two-corner parities. Both
///   describe the cube; only the count describes the table.
/// - It is **octahedrally invariant**: relabelling the corners by any of the 48
///   cube symmetries permutes the case index and leaves the count unchanged —
///   verified here at 0 violations over all 48 x 256 pairs. P-169's C1 is only a
///   meaningful check on a reading with that property, since "influence equal
///   within a symmetry class" is a *consequence* of equivariance and an unequal
///   influence therefore localises a real defect.
///
/// The count is a small integer — measured maximum 5 over the plain table, with
/// histogram `[2, 16, 50, 80, 76, 32]` over counts `0..=5` — so it needs four
/// bits and the per-bit analysis is four functions rather than thirty-two. Two of
/// those four are degenerate and a consumer must expect it: **bit 3 is the
/// constant zero** (no case reaches 8 triangles), and **bit 0 is exactly
/// [`Bool8::parity`]** — the triangle count is odd precisely when an odd number
/// of corners are inside, so that bit carries all its Fourier weight at degree 8
/// and has every influence equal to 1. Bits 1 and 2 are the informative ones.
pub(crate) fn shipped_triangle_counts() -> [u32; 256] {
    let mut out = [0u32; 256];
    for (slot, case) in out.iter_mut().zip(CASES.iter()) {
        *slot = u32::from(case.count);
    }
    out
}

/// The shipped 256-case table's edge-occupancy bitmask per case: bit `e` is set
/// when the case's triangles place a vertex on cube edge `e`.
///
/// Read from the emitted triangles, then checked against the sign rule. Edge `e`
/// is *cut* when its two corners have opposite signs, which for case index `x` is
/// `x_a xor x_b` with `[a, b] = EDGE_CORNERS[e]`; Marching Cubes places exactly
/// one vertex on each cut edge and none elsewhere, so the triangle-derived mask
/// and the sign-derived mask must be identical. That is precisely what
/// `validate_table()`'s `triangles_on_uncut_edges` and `cut_edge_mismatch` gate,
/// and re-deriving it here gives R-169 its "cheap independent check" for free.
///
/// Because each bit is a two-corner parity, this reading's spectrum is known
/// before it is computed and measurement confirms it exactly: every bit is
/// `chi_a chi_b`, a single Fourier coefficient at degree 2, two ANF terms, degree
/// 2, total influence 2, influence **1** on each of its own two corners and `0`
/// on the other six. That makes it a *third* calibration alongside parity and
/// majority — a function whose whole spectrum is one coefficient at a degree
/// neither of the other two occupies.
///
/// It is a **poor primary, and for P-169 an actively misleading one.** A single
/// edge bit is attached to a specific pair of corners, so it is not invariant
/// under the cube group even though the mask as a whole is equivariant (the
/// symmetries permute the edges). Its eight influences are therefore `(1, 1, 0,
/// 0, 0, 0, 0, 0)` on the *correct* shipped table, and P-169's C1 read against
/// this reading would report "unequal within class" for a table with no defect
/// whatsoever. The verified octahedral invariance of
/// [`shipped_triangle_counts`] is what makes C1's equality a real test.
///
/// # Panics
///
/// If any case's triangles reference an edge the sign rule says is not cut, or
/// leave a cut edge unreferenced.
pub(crate) fn shipped_edge_masks() -> [u32; 256] {
    let mut out = [0u32; 256];
    for (case_index, (slot, case)) in out.iter_mut().zip(CASES.iter()).enumerate() {
        let mut emitted = 0u32;
        for tri in case.triangles.iter().take(case.count as usize) {
            for &code in tri {
                if code < CENTROID_BASE {
                    emitted |= 1u32 << code;
                }
            }
        }

        let mut cut = 0u32;
        for (e, corners) in EDGE_CORNERS.iter().enumerate() {
            let inside_lo = corner_inside(case_index as u8, corners[0]);
            let inside_hi = corner_inside(case_index as u8, corners[1]);
            if inside_lo != inside_hi {
                cut |= 1u32 << e;
            }
        }

        assert_eq!(
            emitted, cut,
            "case {case_index}: the shipped triangles occupy edges {emitted:#014b} \
             but the corner signs cut edges {cut:#014b} — the table and the sign \
             rule disagree, which validate_table() should already have caught"
        );
        *slot = emitted;
    }
    out
}

/// The shipped 256-case table's cycle-centroid count per case, read from
/// `isomesh::marching_cubes::table::CASES`.
///
/// A centroid is the interior vertex a cycle gets when no chord-safe apex exists.
/// `CASES` is built at `segment_links(case, 0)` — the all-separate resolution,
/// which is Marching Cubes proper — and plain Marching Cubes tops out at cycles
/// of length seven, below the length at which the centroid path is reached. So
/// this reading is the **constant zero function**, and measurement confirms it:
/// degree 0, all weight at degree 0, **zero** ANF terms (the empty polynomial,
/// not a constant-1 term), every influence zero, noise stability 1 at every
/// `rho`.
///
/// It is provided because the task of choosing "the output" should be settled by
/// measurement rather than by assertion, and a constant is a legitimate answer
/// that says something real — the centroid machinery is dead code for the
/// unambiguous table. It is a *bad* primary for exactly that reason: a constant
/// function's spectrum cannot falsify C2 either way.
pub(crate) fn shipped_centroid_counts() -> [u32; 256] {
    let mut out = [0u32; 256];
    for (slot, case) in out.iter_mut().zip(CASES.iter()) {
        *slot = u32::from(case.centroids);
    }
    out
}

/// A deliberately corrupted copy of a table: one case's output altered, for the
/// vacuity control that says the instrument can report bad news. `which` selects
/// the case, `delta` the alteration. Deterministic.
///
/// The alteration is `xor`, not addition, and that is a choice with a reason: at
/// the level a Boolean function sees, the output is a bundle of bits, and `xor`
/// flips exactly the bits named by `delta` at exactly one point of the domain.
/// Addition would carry, so `delta = 1` on a case whose count is 3 would move two
/// bits and a `delta` chosen to hit a specific bit would be case-dependent.
///
/// # What a single-case flip does to the influences, exactly
///
/// Flipping one point `x0` of one bit-function changes `Inf_i` for every `i` by
/// exactly `+/- 2/256`: only the pair `(x0, x0 xor e_i)` is touched, and it goes
/// from unequal to equal (`-2/256`) or from equal to unequal (`+2/256`)
/// depending on the neighbour. So the influence check detects the corruption
/// **iff the eight neighbours of `x0` do not all agree with each other about
/// whether they matched `x0`** — otherwise all eight influences move the same
/// way and equality survives.
///
/// That is not a caveat to bury; it decides which corruption P-169's vacuity
/// control may use. Measured over all 256 single-case flips of
/// [`shipped_triangle_counts`], per output bit:
///
/// | bit | what it is | flips detected |
/// |---|---|---|
/// | 0 | exactly [`Bool8::parity`] | **0 / 256** |
/// | 1 | informative | 236 / 256 |
/// | 2 | informative | 204 / 256 |
/// | 3 | constant zero | **0 / 256** |
///
/// Bits 0 and 3 are structurally undetectable by this instrument — for parity
/// every neighbour always differs, and for a constant none ever does, so in both
/// cases all eight influences move together. A bench that corrupts bit 0 and
/// concludes "C1 cannot detect a defect" has measured the wrong bit. Corrupt bit
/// 1 or bit 2, at a case whose neighbourhood is mixed — `delta = 2` at case 37
/// and `delta = 4` at case 7 are two that work — and the eight influences split.
///
/// # Panics
///
/// If `which >= 256`.
pub(crate) fn corrupt(values: &[u32; 256], which: usize, delta: u32) -> [u32; 256] {
    assert!(
        which < 256,
        "case {which} does not exist in a 256-case table"
    );
    let mut out = *values;
    out[which] ^= delta;
    out
}
