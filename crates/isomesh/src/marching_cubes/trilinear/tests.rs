//! Tests for the body-saddle classifier.
//!
//! The load-bearing one is
//! [`the_body_saddle_heights_agree_with_the_swept_saddle_roots`]: this module and
//! [`crate::marching_cubes::interior`] locate the same points by constructions
//! that share no arithmetic, so their agreement is evidence rather than a
//! tautology.

// Several tests here compare floats exactly on purpose: a root that must be
// reproduced bit-for-bit, a coordinate that must be untouched by a relabelling.
// An approximate comparison would be a different test, not a weaker one.
#![allow(clippy::float_cmp)]

use super::{
    ALL_INSIDE, BodySaddles, Contours, INTERIOR, MAX_CONTOURS, MAX_INTERIOR_VERTICES, Pair,
    Topology,
};
use crate::Sdf;
use crate::cube::{corner_offset, is_inside};
use crate::fields::ReferenceField;
use crate::marching_cubes::interior::SweptFaces;

/// A deterministic generator, so a census is reproducible and a failure is a
/// fixture rather than a rumour. Numerical Recipes' 64-bit LCG constants.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    /// A value in `[-1, 1)`.
    fn signed(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        // Top 53 bits, so the mantissa is filled and the low-order weakness of an
        // LCG stays out of the result.
        let unit = (self.0 >> 11) as f64 / (1u64 << 53) as f64;
        unit * 2.0 - 1.0
    }

    fn corners(&mut self) -> [f64; 8] {
        let mut f = [0.0; 8];
        for slot in &mut f {
            *slot = self.signed();
        }
        f
    }
}

/// The trilinear interpolant on the unit cell, written out independently of
/// anything in the module under test.
fn trilinear(f: &[f64; 8], p: [f64; 3]) -> f64 {
    let mut sum = 0.0;
    for (i, &value) in f.iter().enumerate() {
        let o = corner_offset(i as u8);
        let mut weight = value;
        for axis in 0..3 {
            weight *= if o[axis] == 1 { p[axis] } else { 1.0 - p[axis] };
        }
        sum += weight;
    }
    sum
}

/// How many connected components the cell's **inside corners** fall into.
///
/// The inside (negative) region within the cell is sampled on an `n³` grid and
/// flood-filled with 6-connectivity; the answer is how many distinct components
/// the cell's own inside corners land in. This is the independent oracle for
/// "is there a tunnel": a tunnel joins two same-signed corner groups *through*
/// the interior, so its corners land in **one** component, and two separate
/// blobs land in two.
///
/// It shares no arithmetic with [`super::Contours::topology`] or with
/// [`super::BodySaddles`] — it evaluates the interpolant and nothing else, which
/// is the only reason its agreement with either is evidence (M-229).
fn corner_groups(f: &[f64; 8], n: usize) -> usize {
    let at = |i: usize, j: usize, k: usize| {
        let p = [
            i as f64 / (n - 1) as f64,
            j as f64 / (n - 1) as f64,
            k as f64 / (n - 1) as f64,
        ];
        trilinear(f, p) < 0.0
    };
    let idx = |i: usize, j: usize, k: usize| (k * n + j) * n + i;
    let mut label = alloc::vec![usize::MAX; n * n * n];
    let mut next = 0usize;
    let mut stack = alloc::vec::Vec::new();
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                if !at(i, j, k) || label[idx(i, j, k)] != usize::MAX {
                    continue;
                }
                let id = next;
                next += 1;
                stack.push((i, j, k));
                label[idx(i, j, k)] = id;
                while let Some((x, y, z)) = stack.pop() {
                    let push = |a: usize,
                                b: usize,
                                c: usize,
                                st: &mut alloc::vec::Vec<_>,
                                lb: &mut alloc::vec::Vec<usize>| {
                        if at(a, b, c) && lb[idx(a, b, c)] == usize::MAX {
                            lb[idx(a, b, c)] = id;
                            st.push((a, b, c));
                        }
                    };
                    if x > 0 {
                        push(x - 1, y, z, &mut stack, &mut label);
                    }
                    if x + 1 < n {
                        push(x + 1, y, z, &mut stack, &mut label);
                    }
                    if y > 0 {
                        push(x, y - 1, z, &mut stack, &mut label);
                    }
                    if y + 1 < n {
                        push(x, y + 1, z, &mut stack, &mut label);
                    }
                    if z > 0 {
                        push(x, y, z - 1, &mut stack, &mut label);
                    }
                    if z + 1 < n {
                        push(x, y, z + 1, &mut stack, &mut label);
                    }
                }
            }
        }
    }
    // How many distinct components do the *inside corners* fall into?
    let mut seen = alloc::collections::BTreeSet::new();
    for (c, &value) in f.iter().enumerate() {
        if !is_inside(value) {
            continue;
        }
        let g = |bit: usize| if (c >> bit) & 1 == 1 { n - 1 } else { 0 };
        let l = label[idx(g(0), g(1), g(2))];
        if l != usize::MAX {
            seen.insert(l);
        }
    }
    seen.len()
}

/// Grosso's `v0 = (0,0,0)`, `v1 = (1,0,0)`, `v2 = (0,1,0)`, `v3 = (1,1,0)`,
/// `v4 = (0,0,1)`, `v5 = (1,0,1)`, `v6 = (0,1,1)`, `v7 = (1,1,1)`.
///
/// Every formula in the module is written directly in his indices, which is only
/// sound because they are also ours. It is a coincidence, not a design, so it is
/// pinned: if `cube.rs` ever renumbers, this fails first and loudly rather
/// than the meshes going subtly wrong.
#[test]
fn grosso_corner_numbering_is_ours() {
    let grosso: [[u32; 3]; 8] = [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [1, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 1],
    ];
    for (i, expected) in grosso.iter().enumerate() {
        assert_eq!(
            corner_offset(i as u8),
            *expected,
            "corner {i} is not where Grosso puts it"
        );
    }
}

/// The three coefficients against the expression they were derived from.
///
/// `P(u) = (i₀ − g₀)(g̃₁ − g̃₀) − (i₀ − g̃₀)(g₁ − g₀)` at `i₀ = 0`, evaluated
/// directly from the four edge interpolations, must equal `a·u² + b·u + c`. The
/// two share no term ordering, so a mis-grouped coefficient shows up here.
#[test]
fn the_coefficients_reproduce_the_face_hyperbola_difference() {
    let mut rng = Lcg::new(0x0000_A002_D000_0001);
    let mut worst: f64 = 0.0;
    for _ in 0..20_000 {
        let f = rng.corners();
        let [a, b, c] = BodySaddles::coefficients(&f);
        for step in 0..=8 {
            let u = f64::from(step) / 8.0;
            let s = 1.0 - u;
            let g0 = f[0] * s + f[1] * u;
            let g1 = f[2] * s + f[3] * u;
            let h0 = f[4] * s + f[5] * u;
            let h1 = f[6] * s + f[7] * u;
            let direct = (-g0) * (h1 - h0) - (-h0) * (g1 - g0);
            let expanded = (a * u + b) * u + c;
            worst = worst.max((direct - expanded).abs());
        }
    }
    assert!(worst < 1e-12, "coefficients disagree by up to {worst:e}");
    std::println!("measured: worst coefficient disagreement {worst:e} over 20,000 cells");
}

/// Every saddle the classifier reports inside the cell is **on** the surface.
///
/// This is the geometric claim the whole construction rests on — Grosso's
/// equation (5), that the lines joining solutions at opposite faces lie entirely
/// on the level set. Checked against an independent evaluation of the interpolant.
#[test]
fn every_body_saddle_lies_on_the_level_set() {
    let mut rng = Lcg::new(0x0000_A002_D000_0002);
    let mut checked = 0usize;
    let mut worst: f64 = 0.0;
    for _ in 0..200_000 {
        let f = rng.corners();
        let saddles = BodySaddles::of(&f);
        let Some(hexagon) = saddles.inner_hexagon() else {
            continue;
        };
        for vertex in hexagon {
            worst = worst.max(trilinear(&f, vertex).abs());
            checked += 1;
        }
    }
    assert!(
        checked > 0,
        "no cell in the sweep had six saddles, so this measured nothing"
    );
    assert!(
        worst < 1e-9,
        "a hexagon vertex is {worst:e} off the surface"
    );
    std::println!(
        "measured: {checked} hexagon vertices, worst distance from the level set {worst:e}"
    );
}

/// **The cross-check.** This module and [`SweptFaces`] find the same body saddles.
///
/// [`SweptFaces`] sweeps a plane between two opposite faces and looks for heights
/// where the plane's bilinear saddle sits on the surface. This module intersects
/// the two faces' hyperbolas and reads a coordinate off. Nothing is shared: one
/// solves a quadratic in the sweep height `t`, the other a quadratic in the face
/// coordinate `u` and then interpolates.
///
/// They agree. Recorded rather than assumed — the relationship was not obvious in
/// advance and is now measured (M-206).
#[test]
fn the_body_saddle_heights_agree_with_the_swept_saddle_roots() {
    let mut rng = Lcg::new(0x0000_A002_D000_0003);
    let mut compared = 0usize;
    let mut worst: f64 = 0.0;

    for _ in 0..200_000 {
        let f = rng.corners();
        let saddles = BodySaddles::of(&f);
        if !saddles.has_inner_hexagon() {
            continue;
        }

        // The `w` pair: the low face is `z = 0`, whose corners in the cyclic order
        // `SweptFaces` wants — `A`/`C` one diagonal, `B`/`D` the other, `lo[k]` and
        // `hi[k]` the ends of a cell edge — are 0, 1, 3, 2.
        let Ok(swept) = SweptFaces::new([f[0], f[1], f[3], f[2]], [f[4], f[5], f[7], f[6]]) else {
            continue;
        };

        let mut theirs = [0.0f64; 2];
        let mut count = 0;
        for t in swept.numerator_roots() {
            if t > 0.0 && t < 1.0 && count < 2 {
                theirs[count] = t;
                count += 1;
            }
        }
        if count != 2 {
            continue;
        }

        let mut ours = saddles.axis(2);
        if ours[0] > ours[1] {
            ours.swap(0, 1);
        }
        if theirs[0] > theirs[1] {
            theirs.swap(0, 1);
        }

        for k in 0..2 {
            worst = worst.max((ours[k] - theirs[k]).abs());
        }
        compared += 1;
    }

    assert!(
        compared > 100,
        "only {compared} cells reached the comparison, which measures nothing"
    );
    assert!(
        worst < 1e-9,
        "the two constructions disagree by up to {worst:e}"
    );
    std::println!(
        "measured: {compared} cells cross-checked against the swept saddle, worst gap {worst:e}"
    );
}

/// The hexagon is a closed ring whose every edge is parallel to an axis.
///
/// Grosso's Proposition 2. Consecutive vertices must differ in exactly one
/// coordinate — including the wrap from the sixth back to the first, which is the
/// half a corrupted listing would get wrong.
#[test]
fn hexagon_edges_are_axis_parallel_and_close() {
    let mut rng = Lcg::new(0x0000_A002_D000_0004);
    let mut rings = 0usize;
    for _ in 0..200_000 {
        let f = rng.corners();
        let Some(hexagon) = BodySaddles::of(&f).inner_hexagon() else {
            continue;
        };
        let mut per_axis = [0usize; 3];
        for k in 0..6 {
            let a = hexagon[k];
            let b = hexagon[(k + 1) % 6];
            let mut differing = 0;
            for (axis, count) in per_axis.iter_mut().enumerate() {
                if a[axis] != b[axis] {
                    *count += 1;
                    differing += 1;
                }
            }
            assert_eq!(
                differing,
                1,
                "hexagon edge {k}→{} changes {differing} coordinates, not one",
                (k + 1) % 6
            );
        }
        // Each axis must supply exactly two of the six edges, or the ring folds
        // back on itself instead of enclosing the tunnel.
        for (axis, &used) in per_axis.iter().enumerate() {
            assert_eq!(
                used, 2,
                "axis {axis} supplies {used} hexagon edges, not two"
            );
        }
        rings += 1;
    }
    assert!(rings > 0, "no hexagon was built, so this measured nothing");
    std::println!("measured: {rings} inner hexagons, all closed and axis-parallel");
}

/// Relabelling the axes cannot change how many saddles a cell has.
///
/// The classifier solves for `u` and derives `v` and `w`, so it is deliberately
/// asymmetric in the axes. The *geometry* is not: the body saddles are a property
/// of the interpolant, and a cyclic relabelling of the axes has to leave their
/// count alone. M-204 is the reason this is checked rather than assumed — that
/// entry records a rotation property that held by algebra and failed by IEEE.
#[test]
fn the_classification_is_invariant_under_axis_relabelling() {
    /// New axis 0 is old axis 2, new 1 is old 0, new 2 is old 1.
    const PERM: [usize; 3] = [2, 0, 1];

    fn relabel(f: &[f64; 8]) -> [f64; 8] {
        let mut out = [0.0; 8];
        for (n, slot) in out.iter_mut().enumerate() {
            let mut old = 0usize;
            for (axis, &target) in PERM.iter().enumerate() {
                old |= ((n >> axis) & 1) << target;
            }
            *slot = f[old];
        }
        out
    }

    let mut rng = Lcg::new(0x0000_A002_D000_0005);
    let mut disagreements = 0usize;
    let mut hexagons = 0usize;
    let mut cells = 0usize;

    for _ in 0..200_000 {
        let f = rng.corners();
        let g = relabel(&f);

        // The relabelling is only meaningful if it really is the same field seen
        // through permuted coordinates. Checked, not assumed.
        // The relabelled field at `probe` is the original at `mapped`, where
        // `mapped[PERM[axis]] = probe[axis]` — the inverse permutation, which is
        // what carries new coordinates back to old ones.
        let probe = [0.31, 0.57, 0.13];
        let mut mapped = [0.0f64; 3];
        for (axis, &target) in PERM.iter().enumerate() {
            mapped[target] = probe[axis];
        }
        assert!((trilinear(&g, probe) - trilinear(&f, mapped)).abs() < 1e-12);

        let a = BodySaddles::of(&f);
        let b = BodySaddles::of(&g);
        if a.inside_count() != b.inside_count() {
            disagreements += 1;
        }
        if a.has_inner_hexagon() {
            hexagons += 1;
        }
        cells += 1;
    }

    assert!(
        hexagons > 0,
        "no hexagon was reached, so this measured nothing"
    );
    assert_eq!(
        disagreements, 0,
        "{disagreements} of {cells} cells classified differently after relabelling"
    );
    std::println!(
        "measured: {cells} cells relabelled, {hexagons} with six saddles, 0 disagreements"
    );
}

/// A tangency is one intersection point, not two.
///
/// When the discriminant vanishes the two face hyperbolas touch rather than cross.
/// Reporting a double root as two solutions would let a degenerate, zero-area
/// "hexagon" claim six saddles and be meshed as a tunnel. The authors'
/// implementation drops the point entirely; this keeps it, once.
#[test]
fn a_double_root_is_one_point_not_two() {
    // A tangency needs `b² = 4ac` **exactly**, which a bisection search cannot
    // deliver — it lands near the crossing, where the discriminant is a small
    // non-zero number and the code correctly reports two nearly-equal roots. So
    // the fixture is constructed backwards from the coefficients instead, out of
    // exact binary fractions, so that the discriminant is exactly `+0.0`:
    //
    //   twist_lo = 2, twist_hi = 0, f₀ = 0
    //     ⟹ a = du_hi·2 = 1        (du_hi = 0.5)
    //       b = f₄·2 + du_hi·dv_lo − du_lo·dv_hi = −1
    //       c = f₂·f₄ = 0.25
    //   and b² − 4ac = 1 − 1 = 0, with the double root at u = −b/2a = 0.5.
    let f = [0.0f64, 1.5, 0.5, 4.0, 0.5, 1.0, 2.0, 2.5];

    let [a, b, c] = BodySaddles::coefficients(&f);
    assert_eq!([a, b, c], [1.0, -1.0, 0.25], "the fixture drifted");
    let discriminant = b * b - 4.0 * a * c;
    assert_eq!(
        discriminant, 0.0,
        "the fixture is not an exact tangency; it measures nothing"
    );

    let saddles = BodySaddles::of(&f);
    // The double root is at 0.5, inside the cell — so this is not passing merely
    // because both roots fell out of range.
    assert!(
        saddles.inside_mask() & 0b01 != 0,
        "the tangency's own root was not reported at all"
    );
    assert!(
        saddles.inside_mask() & 0b10 == 0,
        "a tangency was reported as two distinct intersection points"
    );
    assert!(
        !saddles.has_inner_hexagon(),
        "a tangency was reported as a full hexagon"
    );
    std::println!(
        "measured: exact tangency, double root {}, mask {:#08b}",
        saddles.axis(0)[0],
        saddles.inside_mask()
    );
}

/// A vanishing leading coefficient is a smaller polynomial, not an absence.
///
/// With `a == 0` the equation is linear and still has a root. The textbook
/// quadratic formula divides by `2a` and loses it; the authors' implementation
/// inherits that and reports no solution. Keeping the root matters because
/// Grosso §5.3 chooses a cell's interior vertex by counting face pairs with a
/// *single* solution — a count the textbook form cannot produce.
#[test]
fn a_linear_equation_keeps_its_root() {
    // `a = du_hi·twist_lo − du_lo·twist_hi`, so making **both** twists zero makes
    // `a` zero whatever the edge differences are. `f₃` and `f₇` are then forced:
    //   twist_lo = (f₀ + f₃) − (f₁ + f₂) = 0  ⟹  f₃ = −1.5
    //   twist_hi = (f₄ + f₇) − (f₅ + f₆) = 0  ⟹  f₇ = −2.5
    // The remaining six are chosen so that `b` and `c` stay non-trivial and the
    // single root lands inside the cell rather than out of it, which is what
    // makes the mask assertion below say something.
    let f = [0.25f64, -0.5, -0.75, -1.5, 0.5, 0.0, -2.0, -2.5];

    let [a, b, c] = BodySaddles::coefficients(&f);
    assert_eq!(a, 0.0, "the fixture does not have a vanishing leading term");
    assert_ne!(b, 0.0, "the fixture is not linear, it is constant");
    assert_eq!([b, c], [-1.375, 0.125], "the fixture drifted");

    let saddles = BodySaddles::of(&f);
    let root = saddles.axis(0)[0];
    assert!(
        (b * root + c).abs() < 1e-12,
        "the linear root does not satisfy the equation"
    );

    // And it is reported, not silently dropped.
    let expected_inside = root > 0.0 && root < 1.0;
    assert_eq!(
        saddles.inside_mask() & 1 != 0,
        expected_inside,
        "the linear root's in-range bit does not match its value {root}"
    );
    std::println!("measured: linear case root {root:.12}, in range {expected_inside}");
}

/// How often a cell has all six body saddles.
///
/// Pinned in both directions: it is the population A-002e has to find a field to
/// reach, and a change in either direction is a change in the classifier.
#[test]
fn how_often_a_cell_has_six_body_saddles() {
    let mut rng = Lcg::new(0x0000_A002_D000_0006);
    let mut histogram = [0usize; 7];
    const CELLS: usize = 200_000;
    for _ in 0..CELLS {
        let saddles = BodySaddles::of(&rng.corners());
        histogram[saddles.inside_count() as usize] += 1;
    }
    assert_eq!(histogram.iter().sum::<usize>(), CELLS);
    assert!(
        histogram[6] > 0,
        "no cell reached six saddles, so the census measured nothing"
    );
    std::println!("measured: saddle-count histogram over {CELLS} random cells: {histogram:?}");
}

/// The mask constant and the predicate cannot drift apart.
#[test]
fn the_hexagon_predicate_is_the_full_mask() {
    let mut rng = Lcg::new(0x0000_A002_D000_0007);
    for _ in 0..50_000 {
        let saddles = BodySaddles::of(&rng.corners());
        assert_eq!(
            saddles.has_inner_hexagon(),
            saddles.inside_mask() == ALL_INSIDE
        );
        assert_eq!(
            saddles.has_inner_hexagon(),
            saddles.inside_count() == 6,
            "count and mask disagree"
        );
    }
}

/// Both widths classify, and the same call twice gives the same bits.
#[test]
fn the_classification_is_deterministic_in_both_widths() {
    let mut rng = Lcg::new(0x0000_A002_D000_0008);
    for _ in 0..20_000 {
        let wide = rng.corners();
        let mut narrow = [0.0f32; 8];
        for (slot, &value) in narrow.iter_mut().zip(wide.iter()) {
            *slot = value as f32;
        }

        let a = BodySaddles::of(&wide);
        let b = BodySaddles::of(&wide);
        assert_eq!(a.inside_mask(), b.inside_mask());
        for axis in 0..3 {
            assert_eq!(
                a.axis(axis),
                b.axis(axis),
                "axis {axis} is not reproducible"
            );
        }

        let n = BodySaddles::of(&narrow);
        assert_eq!(n.inside_mask(), BodySaddles::of(&narrow).inside_mask());
    }
}

/// How often the seven reference fields reach six body saddles.
///
/// **This is the measurement A-002e turns on.** M-40 established that five of the
/// seven never produce an ambiguous *face* at all. Interior ambiguity is rarer
/// still — 0.57% of uniformly random cells, and a signed-distance field is very
/// far from uniformly random — so if no reference field reaches six saddles then
/// the T-001 field gates and the golden fixture cannot exercise a tunnel, and an
/// eighth field has to be built for them to.
///
/// Recorded rather than gated: a zero here is a fact about the fields, not a
/// failure. What *is* asserted is that the sweep visited surface cells at all, so
/// a zero cannot be an empty loop reported as a result.
#[test]
fn how_often_the_reference_fields_reach_six_body_saddles() {
    let mut reached = 0u64;
    let mut surface_total = 0u64;

    crate::for_each_reference_field!(f64, |name, field| {
        // No `return` in this body: `for_each_reference_field!` is a macro, not a
        // closure, and one would exit the whole test at `sphere` (M-199).
        for samples in [17u32, 33, 65] {
            let (lo, hi) = field.domain();
            let h = (hi[0] - lo[0]) / f64::from(samples - 1);
            let mut surface = 0u64;
            let mut histogram = [0u64; 7];

            for z in 0..samples - 1 {
                for y in 0..samples - 1 {
                    for x in 0..samples - 1 {
                        let mut corner = [0.0f64; 8];
                        let mut case = 0u8;
                        for (c, slot) in corner.iter_mut().enumerate() {
                            let o = corner_offset(c as u8);
                            *slot = field.sample([
                                lo[0] + h * f64::from(x + o[0]),
                                lo[1] + h * f64::from(y + o[1]),
                                lo[2] + h * f64::from(z + o[2]),
                            ]);
                            if is_inside(*slot) {
                                case |= 1 << c;
                            }
                        }
                        if case == 0 || case == 255 {
                            continue;
                        }
                        surface += 1;
                        histogram[BodySaddles::of(&corner).inside_count() as usize] += 1;
                    }
                }
            }

            surface_total += surface;
            reached += histogram[6];
            std::println!(
                "measured: {name} at {samples}^3 -> {surface} surface cells, \
                 saddle histogram {histogram:?}"
            );
        }
    });

    assert!(
        surface_total > 0,
        "the sweep visited no surface cell, so it measured nothing"
    );
    std::println!(
        "measured: {reached} of {surface_total} reference-field surface cells have six body saddles"
    );
}

/// The eighth reference field really does contain the configuration it exists for.
///
/// **A-002e's acceptance, and M-44's rule as a test.** `noise_cavity` was added for
/// one reason: every other reference field has an interior-ambiguity rate of
/// exactly zero (M-208), so without it the tunnel case — the thing MC33's interior
/// rule exists for — is unreachable by this crate's own suite and the T-001 gates
/// and golden fixture can never exercise it.
///
/// A field that reached no six-saddle cell would be a gate that measures nothing,
/// so the count is asserted **non-zero at every golden resolution** rather than
/// merely recorded. Pinned in both directions: an increase is as much a change to
/// the field as a decrease.
#[test]
fn the_tunnel_field_actually_contains_tunnels() {
    let field = crate::fields::noise_cavity::<f64>();
    let (lo, hi) = field.domain();

    let mut counts = alloc::vec::Vec::new();
    for samples in [17u32, 25, 33] {
        let h = (hi[0] - lo[0]) / f64::from(samples - 1);
        let mut six = 0u64;
        for z in 0..samples - 1 {
            for y in 0..samples - 1 {
                for x in 0..samples - 1 {
                    let mut corner = [0.0f64; 8];
                    let mut case = 0u8;
                    for (c, slot) in corner.iter_mut().enumerate() {
                        let o = corner_offset(c as u8);
                        *slot = field.sample([
                            lo[0] + h * f64::from(x + o[0]),
                            lo[1] + h * f64::from(y + o[1]),
                            lo[2] + h * f64::from(z + o[2]),
                        ]);
                        if is_inside(*slot) {
                            case |= 1 << c;
                        }
                    }
                    if case == 0 || case == 255 {
                        continue;
                    }
                    if BodySaddles::of(&corner).inside_count() == 6 {
                        six += 1;
                    }
                }
            }
        }
        assert!(
            six > 0,
            "noise_cavity at {samples}^3 has no six-saddle cell, so every gate \
             that rests on it measures nothing"
        );
        counts.push(six);
    }

    assert_eq!(
        counts,
        alloc::vec![3, 4, 4],
        "the tunnel population moved; the field or the classifier changed"
    );
    std::println!("measured: noise_cavity six-saddle cells at 17/25/33 -> {counts:?}");
}

// ─── A-002f: contours ───────────────────────────────────────────────────────

/// Every cut edge lies on exactly one ring, and the rings close.
///
/// Exhaustive over all 256 cases × all 64 face-resolution masks, in the shape
/// `validate_decider_table` uses — 16,384 combinations, not the handful some
/// field happens to produce.
#[test]
fn every_cut_edge_lies_on_exactly_one_closed_ring() {
    use crate::cube::{EDGE_CORNERS, corner_inside};
    use crate::marching_cubes::table::segment_links;

    let mut worst_rings = 0usize;
    let mut worst_len = 0usize;
    for case in 0..=255u8 {
        for mask in 0..64u8 {
            let contours = Contours::of(case, mask);

            let mut cut = 0u16;
            for (edge, [lo, hi]) in EDGE_CORNERS.into_iter().enumerate() {
                if corner_inside(case, lo) != corner_inside(case, hi) {
                    cut |= 1 << edge;
                }
            }

            let mut seen = 0u16;
            let next = segment_links(case, mask);
            for r in 0..contours.count() {
                let ring = contours.ring(r);
                assert!(
                    ring.len() >= 3,
                    "case {case:#010b} mask {mask:#08b}: ring of {}",
                    ring.len()
                );
                for (k, &e) in ring.iter().enumerate() {
                    assert_eq!(seen & (1 << e), 0, "edge {e} is on two rings");
                    seen |= 1 << e;
                    // The ring really is the link walk, closed at the wrap.
                    let follows = ring[(k + 1) % ring.len()];
                    assert_eq!(
                        next[e as usize], follows,
                        "case {case:#010b} mask {mask:#08b}: ring does not follow the links"
                    );
                }
                worst_len = worst_len.max(ring.len());
            }
            assert_eq!(
                seen, cut,
                "case {case:#010b} mask {mask:#08b}: rings do not cover the cut edges"
            );
            worst_rings = worst_rings.max(contours.count());
        }
    }
    // Four rings is the tetrahedral case; twelve is the single contour through
    // every cut edge. Pinned so an over-sized buffer is as visible as an
    // under-sized one.
    assert_eq!(worst_rings, MAX_CONTOURS);
    assert_eq!(worst_len, 12);
    std::println!(
        "measured: 16,384 combinations, at most {worst_rings} rings, longest {worst_len}"
    );
}

/// Grosso lists the ring lengths his construction can produce; ours produces
/// exactly that set.
///
/// *"Contours can have 3, 4, 5, 6, 7, 8, 9 or 12 vertices."* Ten and eleven are
/// **absent**, which is the interesting half — a ring of eleven would mean one
/// cut edge left over from a twelve-edge cell, and the parity of the face walk
/// forbids it.
#[test]
fn ring_lengths_are_the_set_the_paper_lists() {
    let mut seen = [false; 13];
    for case in 0..=255u8 {
        for mask in 0..64u8 {
            let contours = Contours::of(case, mask);
            for r in 0..contours.count() {
                seen[contours.ring(r).len()] = true;
            }
        }
    }
    let produced: alloc::vec::Vec<usize> = (0..13).filter(|&k| seen[k]).collect();
    assert_eq!(
        produced,
        alloc::vec![3, 4, 5, 6, 7, 8, 9, 12],
        "the ring-length set moved"
    );
    std::println!("measured: ring lengths produced = {produced:?}");
}

/// The tunnel and the twelve-vertex contour are told apart by the ring count,
/// and both are reachable.
///
/// **A-002d deliberately could not make this call** — six body saddles mean one
/// or the other and the saddles do not say which. This is the test that the ring
/// count does say, and that neither answer is dead code.
#[test]
fn six_saddles_split_into_tunnels_and_twelve_vertex_contours() {
    use crate::cube::is_inside as inside;
    use crate::marching_cubes::ambiguity::joined_mask;
    use crate::marching_cubes::table::AMBIGUOUS_FACES;

    let mut rng = Lcg::new(0x0000_A002_F000_0001);
    let mut tunnels = 0usize;
    let mut twelves = 0usize;
    let mut disks = 0usize;
    let mut separate = 0usize;

    for _ in 0..400_000 {
        let f = rng.corners();
        let mut case = 0u8;
        for (c, &v) in f.iter().enumerate() {
            if inside(v) {
                case |= 1 << c;
            }
        }
        if case == 0 || case == 255 {
            continue;
        }
        let mask = joined_mask(&f, AMBIGUOUS_FACES[case as usize]);
        let contours = Contours::of(case, mask);
        match contours.topology(&BodySaddles::of(&f)) {
            Topology::Disks => disks += 1,
            Topology::Tunnel => tunnels += 1,
            Topology::TwelveVertexContour => twelves += 1,
            Topology::SeparateDisks => separate += 1,
        }
    }

    assert!(
        tunnels > 0,
        "no tunnel was reached, so the branch is untested"
    );
    assert!(
        twelves > 0,
        "no twelve-vertex contour was reached, so that branch is untested"
    );
    // Recorded as zero rather than asserted absent: `[9,3]` is reachable but not
    // by uniform sampling, which is exactly the fixture trap M-228 fell into.
    // `a_nine_and_three_cell_is_refused` carries the searched configuration.
    assert_eq!(
        separate, 0,
        "uniform sampling reached a separate-disks cell, which no sweep has before"
    );
    std::println!(
        "measured: {disks} disks, {tunnels} tunnels, {twelves} twelve-vertex contours, \
         {separate} separate-disk cells"
    );
}

/// Grosso's Corollary 4, checked against the rings rather than assumed.
///
/// *"If the three quadratic equations have less than two solutions within the
/// range, or if only one of the equations has two solutions, then the contour has
/// 3, 4, 5 or 6 vertices and does not intersect an ambiguous face."*
///
/// Recorded as a census rather than asserted as a bare implication, because the
/// interesting part is the shape of the exceptions if there are any.
#[test]
fn few_saddles_mean_short_rings() {
    use crate::cube::is_inside as inside;
    use crate::marching_cubes::ambiguity::joined_mask;
    use crate::marching_cubes::table::AMBIGUOUS_FACES;

    let mut rng = Lcg::new(0x0000_A002_F000_0002);
    let mut checked = 0usize;
    let mut violations = 0usize;
    let mut worst = 0usize;

    for _ in 0..400_000 {
        let f = rng.corners();
        let mut case = 0u8;
        for (c, &v) in f.iter().enumerate() {
            if inside(v) {
                case |= 1 << c;
            }
        }
        if case == 0 || case == 255 {
            continue;
        }
        let saddles = BodySaddles::of(&f);
        if saddles.inside_count() >= 2 {
            continue;
        }
        let mask = joined_mask(&f, AMBIGUOUS_FACES[case as usize]);
        let contours = Contours::of(case, mask);
        let longest = contours.longest();
        worst = worst.max(longest);
        if longest > 6 {
            violations += 1;
        }
        checked += 1;
    }

    assert!(checked > 1000, "only {checked} cells reached the corollary");
    assert_eq!(
        violations, 0,
        "{violations} of {checked} cells with fewer than two saddles had a ring longer than six"
    );
    std::println!(
        "measured: {checked} cells with fewer than two body saddles, longest ring {worst}"
    );
}

/// The line inventory is checked against the inner hexagon, not against the
/// implementation it agrees with.
///
/// **This is the test that earns the crossed `u` pairing (M-215).** The obvious
/// reading — solution `k` gives the `u`-line at `(v[k], w[k])` — is wrong, and
/// only the hexagon says so: a pair of ring vertices differing in exactly one
/// coordinate *is* a line, so reading the ring off gives the pairing directly.
/// Where all six saddles are inside, every pair must report both its lines, and
/// each line's two fixing coordinates must match a hexagon edge.
#[test]
fn the_line_inventory_agrees_with_the_inner_hexagon() {
    let mut rng = Lcg::new(0x0000_A002_6000_0001);
    let mut checked = 0usize;
    for _ in 0..200_000 {
        let f = rng.corners();
        let saddles = BodySaddles::of(&f);
        let Some(hexagon) = saddles.inner_hexagon() else {
            continue;
        };

        // Every pair has both lines when the hexagon exists.
        for pair in Pair::ALL {
            assert_eq!(
                saddles.line_count(pair),
                2,
                "a full hexagon must give every pair both lines"
            );
        }

        // Each hexagon edge changes exactly one coordinate; the axis it changes
        // names the pair, and the two coordinates it holds fixed must be the ones
        // `lines` requires for one of that pair's two lines.
        for k in 0..6 {
            let a = hexagon[k];
            let b = hexagon[(k + 1) % 6];
            let axis = (0..3)
                .find(|&ax| a[ax] != b[ax])
                .expect("a hexagon edge changes a coordinate");
            let pair = Pair::ALL[axis];

            // Recover the line by asking `point_on_line` for both of the pair's
            // lines and seeing which reproduces the held coordinates.
            let mut matched = 0;
            for line in 0..2 {
                let probe = saddles.point_on_line_for_test(pair, line, a[axis]);
                if probe == a {
                    matched += 1;
                }
            }
            assert_eq!(
                matched, 1,
                "hexagon edge {k} along axis {axis} matches {matched} of the \
                 pair's two lines, not exactly one — the pairing is wrong"
            );
        }
        checked += 1;
    }
    assert!(checked > 100, "only {checked} hexagons reached the check");
    std::println!("measured: {checked} hexagons agree with the line inventory");
}

/// **The interior vertex lies on the level set.** This is what verifies §5.3.
///
/// The branch selection in `interior_vertex` is transcribed from the authors'
/// program, because the paper determines neither the four-line case nor which of
/// three lines' three pairwise intersections its "two saddle points" are. A
/// transcription cannot be checked by re-reading it; it can be checked by the
/// property it must have. The point is a **saddle of the trilinear interpolant**,
/// so it lies *on* the surface — and an index picked from the wrong axis lands it
/// somewhere the interpolant is not zero.
#[test]
fn the_interior_vertex_lies_on_the_level_set() {
    let mut rng = Lcg::new(0x0000_A002_6000_0002);
    let mut checked = 0usize;
    let mut worst: f64 = 0.0;
    let mut by_total = [0usize; 5];

    for _ in 0..400_000 {
        let f = rng.corners();
        let saddles = BodySaddles::of(&f);
        if saddles.has_inner_hexagon() {
            continue; // A-002h's territory
        }
        let Some(p) = saddles.interior_vertex() else {
            continue;
        };
        assert!(
            p.iter().all(|c| c.is_finite()),
            "interior vertex is not finite: {p:?}"
        );
        let lines: u32 = Pair::ALL.iter().map(|&q| saddles.line_count(q)).sum();
        by_total[lines as usize] += 1;
        worst = worst.max(trilinear(&f, p).abs());
        checked += 1;
    }

    assert!(
        checked > 1000,
        "only {checked} interior vertices were built"
    );
    assert!(
        worst < 1e-9,
        "an interior vertex is {worst:e} off the level set — the branch selection is wrong"
    );
    // Every branch that can fire must have fired, or one of them is untested.
    assert!(
        by_total[2] > 0 && by_total[3] > 0 && by_total[4] > 0,
        "{by_total:?}"
    );
    std::println!(
        "measured: {checked} interior vertices, worst distance from the level set {worst:e}, \
         by line count {by_total:?}"
    );
}

/// The disk path never reaches five or six lines.
///
/// `interior_vertex` returns `None` there, and that arm exists because the paper
/// and the reference both stop at four. Asserting it is unreachable from the disk
/// path is what stops the `None` being a silent wrong answer.
#[test]
fn the_disk_path_never_sees_five_or_six_lines() {
    let mut rng = Lcg::new(0x0000_A002_6000_0003);
    let mut seen_five_or_six = 0usize;
    let mut checked = 0usize;
    for _ in 0..400_000 {
        let f = rng.corners();
        let saddles = BodySaddles::of(&f);
        if saddles.has_inner_hexagon() {
            continue;
        }
        let lines: u32 = Pair::ALL.iter().map(|&q| saddles.line_count(q)).sum();
        if lines >= 5 {
            seen_five_or_six += 1;
        }
        checked += 1;
    }
    assert!(checked > 1000);
    assert_eq!(
        seen_five_or_six, 0,
        "the disk path reached {seen_five_or_six} cells with five or six lines, \
         which `interior_vertex` answers `None` to"
    );
    std::println!("measured: {checked} disk cells, none with five or six lines");
}

/// The fan uses every ring edge exactly twice and never repeats an index.
///
/// Twice is the manifold condition inside one cell: each ring edge is shared by
/// the two triangles either side of it, except that a fan from a ring vertex
/// leaves the ring's own boundary once. Checked over all 16,384 combinations,
/// both with and without an interior vertex.
#[test]
fn the_fan_covers_each_ring_once_and_repeats_no_index() {
    for case in 0..=255u8 {
        for mask in 0..64u8 {
            let contours = Contours::of(case, mask);
            for interior in [false, true] {
                let mut count = 0usize;
                let mut used = [0u32; 13];
                contours.fan(interior, |t| {
                    assert!(
                        t[0] != t[1] && t[1] != t[2] && t[0] != t[2],
                        "case {case:#010b} mask {mask:#08b}: degenerate triangle {t:?}"
                    );
                    for &c in &t {
                        used[c as usize] += 1;
                    }
                    count += 1;
                });
                assert_eq!(
                    count,
                    contours.triangle_count(interior),
                    "case {case:#010b} mask {mask:#08b}: triangle_count disagrees with fan"
                );
                // Every cut edge the rings carry appears at least once, and the
                // interior vertex appears only when one was asked for.
                for r in 0..contours.count() {
                    for &e in contours.ring(r) {
                        assert!(used[e as usize] > 0, "edge {e} was never emitted");
                    }
                }
                if !interior {
                    assert_eq!(
                        used[INTERIOR as usize], 0,
                        "an interior vertex was used without one"
                    );
                }
            }
        }
    }
    std::println!("measured: fan checked over 16,384 combinations, both interior modes");
}

/// A fan from an interior vertex costs `k` triangles against `k − 2`, and the
/// worst case is what the per-cell budget has to cover.
#[test]
fn the_worst_case_triangle_count_is_pinned() {
    let mut worst_plain = 0usize;
    let mut worst_interior = 0usize;
    for case in 0..=255u8 {
        for mask in 0..64u8 {
            let c = Contours::of(case, mask);
            worst_plain = worst_plain.max(c.triangle_count(false));
            worst_interior = worst_interior.max(c.triangle_count(true));
        }
    }
    // Twelve either way: the single twelve-vertex ring gives 10 fanned from a
    // vertex and 12 from an interior one, and four separate rings of three give
    // four. `MAX_TRIANGLES` is 12 and does not have to move.
    assert_eq!((worst_plain, worst_interior), (10, 12));
    assert!(
        worst_interior <= crate::marching_cubes::table::MAX_TRIANGLES,
        "the fan exceeds MAX_TRIANGLES"
    );
    std::println!(
        "measured: worst triangle count {worst_plain} without an interior vertex, \
         {worst_interior} with one (MAX_TRIANGLES is {})",
        crate::marching_cubes::table::MAX_TRIANGLES
    );
}

// ─── A-002h: tunnels ────────────────────────────────────────────────────────

/// Build one cell's tunnel patch and return its triangles.
fn tunnel_patch(f: &[f64; 8]) -> Option<(alloc::vec::Vec<[u8; 3]>, usize)> {
    use crate::cube::is_inside as inside;
    use crate::marching_cubes::ambiguity::joined_mask;
    use crate::marching_cubes::table::AMBIGUOUS_FACES;

    let mut case = 0u8;
    for (c, &v) in f.iter().enumerate() {
        if inside(v) {
            case |= 1 << c;
        }
    }
    if case == 0 || case == 255 {
        return None;
    }
    let saddles = BodySaddles::of(f);
    if !saddles.has_inner_hexagon() {
        return None;
    }
    let contours = Contours::of(case, joined_mask(f, AMBIGUOUS_FACES[case as usize]));
    let mut tris = alloc::vec::Vec::new();
    let unresolved = contours.fan_tunnel(&saddles, f, |t| tris.push(t));
    Some((tris, unresolved))
}

/// **The tunnel patch is manifold inside the cell.** This is what §5.1 is for.
///
/// Every edge of the emitted patch that is *not* a contour edge must be shared by
/// exactly two triangles, and every contour edge by exactly one — because a
/// contour edge lies on a cell face and its second triangle belongs to the
/// neighbouring cell. That is the whole manifoldness argument, and it is the
/// property Chernyaev's tunnel triangulation fails: his lays part of the tunnel
/// *on* the ambiguous face, so two neighbours both claim it and the shared edge
/// ends up with four faces.
#[test]
fn the_tunnel_patch_is_manifold_inside_the_cell() {
    let mut rng = Lcg::new(0x0000_A002_8000_0001);
    let mut cells = 0usize;
    let mut unresolved_total = 0usize;
    let mut worst_tris = 0usize;
    let mut worst_interior = 0usize;

    for _ in 0..400_000 {
        let f = rng.corners();
        let Some((tris, unresolved)) = tunnel_patch(&f) else {
            continue;
        };
        unresolved_total += unresolved;
        worst_tris = worst_tris.max(tris.len());

        // Which codes are contour vertices, and which consecutive pairs are
        // contour edges.
        let mut ring_edge = [[false; 20]; 20];
        {
            use crate::cube::is_inside as inside;
            use crate::marching_cubes::ambiguity::joined_mask;
            use crate::marching_cubes::table::AMBIGUOUS_FACES;
            let mut case = 0u8;
            for (c, &v) in f.iter().enumerate() {
                if inside(v) {
                    case |= 1 << c;
                }
            }
            let contours = Contours::of(case, joined_mask(&f, AMBIGUOUS_FACES[case as usize]));
            for r in 0..contours.count() {
                let ring = contours.ring(r);
                for k in 0..ring.len() {
                    let (a, b) = (ring[k] as usize, ring[(k + 1) % ring.len()] as usize);
                    ring_edge[a.min(b)][a.max(b)] = true;
                }
            }
        }

        // Directed, so this pins the winding as well as the count: a
        // consistently oriented patch traverses every interior edge once in each
        // direction, and every boundary edge once.
        let mut directed = [[0u32; 20]; 20];
        let mut used = [[0u32; 20]; 20];
        let mut interior_seen = [false; 20];
        for t in &tris {
            assert!(
                t[0] != t[1] && t[1] != t[2] && t[0] != t[2],
                "degenerate tunnel triangle {t:?}"
            );
            for &c in t {
                if c >= INTERIOR {
                    interior_seen[c as usize] = true;
                }
            }
            for k in 0..3 {
                let (a, b) = (t[k] as usize, t[(k + 1) % 3] as usize);
                directed[a][b] += 1;
                used[a.min(b)][a.max(b)] += 1;
            }
        }
        for (a, row) in directed.iter().enumerate() {
            for (b, &count) in row.iter().enumerate() {
                assert!(
                    count <= 1,
                    "directed edge ({a},{b}) traversed {count} times — the patch \
                     is not consistently wound"
                );
            }
        }
        worst_interior = worst_interior.max(interior_seen.iter().filter(|&&s| s).count());

        for a in 0..20 {
            for b in a..20 {
                let n = used[a][b];
                if n == 0 {
                    continue;
                }
                let want = if ring_edge[a][b] { 1 } else { 2 };
                assert_eq!(
                    n, want,
                    "edge ({a},{b}) carries {n} faces, wanted {want} — the tunnel \
                     patch is not manifold inside the cell. corners {f:?}"
                );
            }
        }
        cells += 1;
    }

    assert!(cells > 500, "only {cells} tunnel cells were built");
    // **Zero here, and A-020 explained why (M-230).** A contour edge spanning
    // three hexagon steps has no rule and emits nothing. Every configuration ever
    // found to produce one was a nine-and-three case-13 cell, which Corollary 6's
    // bound now classifies as `SeparateDisks` and keeps out of this sweep
    // entirely — see `a_nine_and_three_cell_is_refused_before_triangulation`. So
    // this is a live guard on a case nothing classified as a tunnel has reached,
    // kept because that is a measurement over a sample rather than a proof.
    assert_eq!(
        unresolved_total, 0,
        "{unresolved_total} ring edges spanned three hexagon steps, so the patches \
         above were judged outside the construction's defined domain"
    );
    std::println!(
        "measured: {cells} tunnel patches, all manifold inside the cell; \
         worst {worst_tris} triangles and {worst_interior} interior vertices"
    );
}

/// Every hexagon vertex the tunnel names is on the surface, and the codes stay in
/// range.
#[test]
fn tunnel_interior_vertices_are_hexagon_saddles() {
    let mut rng = Lcg::new(0x0000_A002_8000_0002);
    let mut checked = 0usize;
    let mut worst: f64 = 0.0;
    for _ in 0..200_000 {
        let f = rng.corners();
        let Some((tris, _)) = tunnel_patch(&f) else {
            continue;
        };
        let hexagon = BodySaddles::of(&f)
            .inner_hexagon()
            .expect("a tunnel patch has a hexagon");
        for t in &tris {
            for &c in t {
                if c >= INTERIOR {
                    let k = (c - INTERIOR) as usize;
                    assert!(k < MAX_INTERIOR_VERTICES, "interior code {c} out of range");
                    worst = worst.max(trilinear(&f, hexagon[k]).abs());
                    checked += 1;
                }
            }
        }
    }
    assert!(checked > 1000);
    assert!(
        worst < 1e-9,
        "a tunnel interior vertex is {worst:e} off the surface"
    );
    std::println!("measured: {checked} tunnel interior-vertex references, worst {worst:e} off");
}

/// What a tunnel costs, measured, because it is more than the disk path and the
/// per-cell budgets have to cover it.
///
/// **This is where A-002g's "no budget moves" stops being true (M-217 amended).**
/// A fanned disk adds one cell-local vertex and at most 12 triangles, both
/// already covered by `MAX_CENTROIDS` and `MAX_TRIANGLES`. A tunnel names all six
/// hexagon vertices and can reach 22 triangles, so A-002b has to raise both.
#[test]
fn the_worst_case_tunnel_triangle_count_is_pinned() {
    let mut rng = Lcg::new(0x0000_A002_8000_0003);
    let mut worst_tris = 0usize;
    let mut worst_interior = 0usize;
    let mut cells = 0usize;
    for _ in 0..400_000 {
        let f = rng.corners();
        let Some((tris, _)) = tunnel_patch(&f) else {
            continue;
        };
        worst_tris = worst_tris.max(tris.len());
        let mut seen = [false; 20];
        for t in &tris {
            for &c in t {
                if c >= INTERIOR {
                    seen[c as usize] = true;
                }
            }
        }
        worst_interior = worst_interior.max(seen.iter().filter(|&&s| s).count());
        cells += 1;
    }
    assert!(cells > 500);
    assert_eq!(
        (worst_tris, worst_interior),
        (22, MAX_INTERIOR_VERTICES),
        "the tunnel's worst case moved"
    );
    std::println!(
        "measured: {cells} tunnels, worst {worst_tris} triangles and \
         {worst_interior} interior vertices"
    );
}

/// How often the asymptotic decider has no answer: never, on continuous data.
///
/// **A-002i's reachability, measured before it was implemented (M-220).** A face
/// is *singular* when its bilinear saddle sits exactly on the level set — the two
/// hyperbola branches degenerate into crossing straight lines and the decider's
/// binary choice is between two answers that are both wrong, because the surface
/// passes through the saddle. Grosso 2017 exists for this case and counts **8, 58
/// and 20** of them per 512²×~700 CT volume.
///
/// Here it is **zero**, on every reference field and on 400,000 random cells. The
/// reason is the data, not the code: a singular face needs `v₀·v₂` and `v₁·v₃` to
/// be *bit-identical* `f64`s, which quantised CT voxels collide into readily and
/// a continuous field essentially never does.
///
/// Recorded rather than gated — a zero here is a fact about the fields. What is
/// asserted is that ambiguous faces were reached at all, so the zero cannot be an
/// empty loop reported as a result.
#[test]
fn how_often_a_face_is_singular() {
    use crate::cube::{face_corners, is_inside as inside};
    use crate::marching_cubes::table::AMBIGUOUS_FACES;

    // A face is singular when the bilinear saddle sits exactly on the level set,
    // i.e. the two diagonal products are bit-identical.
    let singular = |v: [f64; 4]| v[0] * v[2] == v[1] * v[3];

    let mut ambiguous = 0u64;
    let mut singular_count = 0u64;
    crate::for_each_reference_field!(f64, |name, field| {
        for samples in [17u32, 33, 65] {
            let (lo, hi) = field.domain();
            let h = (hi[0] - lo[0]) / f64::from(samples - 1);
            let mut amb = 0u64;
            let mut sing = 0u64;
            for z in 0..samples - 1 {
                for y in 0..samples - 1 {
                    for x in 0..samples - 1 {
                        let mut corner = [0.0f64; 8];
                        let mut case = 0u8;
                        for (c, slot) in corner.iter_mut().enumerate() {
                            let o = corner_offset(c as u8);
                            *slot = field.sample([
                                lo[0] + h * f64::from(x + o[0]),
                                lo[1] + h * f64::from(y + o[1]),
                                lo[2] + h * f64::from(z + o[2]),
                            ]);
                            if inside(*slot) {
                                case |= 1 << c;
                            }
                        }
                        let mask = AMBIGUOUS_FACES[case as usize];
                        if mask == 0 {
                            continue;
                        }
                        for axis in 0..3usize {
                            for side in 0..2u8 {
                                let bit = 1u8 << (axis * 2 + side as usize);
                                if mask & bit == 0 {
                                    continue;
                                }
                                let c = face_corners(axis, side);
                                let v = [
                                    corner[c[0] as usize],
                                    corner[c[1] as usize],
                                    corner[c[2] as usize],
                                    corner[c[3] as usize],
                                ];
                                amb += 1;
                                if singular(v) {
                                    sing += 1;
                                }
                            }
                        }
                    }
                }
            }
            ambiguous += amb;
            singular_count += sing;
            if sing > 0 {
                std::println!("  {name} at {samples}^3: {sing} singular of {amb} ambiguous faces");
            }
        }
    });
    std::println!(
        "measured: {singular_count} singular of {ambiguous} ambiguous faces, eight fields"
    );

    // And over random cells, where a tie is far likelier than on a smooth field.
    let mut rng = Lcg::new(0x0000_A002_9000_0001);
    let mut ramb = 0u64;
    let mut rsing = 0u64;
    for _ in 0..400_000 {
        let f = rng.corners();
        let mut case = 0u8;
        for (c, &v) in f.iter().enumerate() {
            if inside(v) {
                case |= 1 << c;
            }
        }
        let mask = AMBIGUOUS_FACES[case as usize];
        for axis in 0..3usize {
            for side in 0..2u8 {
                if mask & (1u8 << (axis * 2 + side as usize)) == 0 {
                    continue;
                }
                let c = face_corners(axis, side);
                let v = [
                    f[c[0] as usize],
                    f[c[1] as usize],
                    f[c[2] as usize],
                    f[c[3] as usize],
                ];
                ramb += 1;
                if singular(v) {
                    rsing += 1;
                }
            }
        }
    }
    std::println!("measured: {rsing} singular of {ramb} ambiguous faces, 400,000 random cells");

    assert!(
        ambiguous > 1000 && ramb > 100_000,
        "the sweep reached {ambiguous} and {ramb} ambiguous faces, too few to mean anything"
    );
    assert_eq!(
        (singular_count, rsing),
        (0, 0),
        "a singular face appeared — A-002i's premise has changed and it is now reachable"
    );
}

/// **The nine-and-three cell, and the classification that now stops it before any
/// triangulation runs (M-230).**
///
/// Marching Cubes' **case 13** — the four alternating corners, the only case with
/// all six faces ambiguous — at particular face resolutions gives a cell with an
/// inner hexagon and contours of **nine and three** vertices. The ring count
/// called it a tunnel, `fan_tunnel` then met a contour edge whose endpoints land
/// three steps apart on the hexagon, and Grosso's construction has no rule for
/// that: it emits nothing and leaves a hole (M-228).
///
/// **The hole was a symptom.** Corollary 6 bounds a tunnel's contours at six, so
/// this cell was never a tunnel, and
/// [`a_nine_and_three_cell_is_not_one_connected_tunnel`] shows independently that
/// its inside region is two blobs. [`Contours::topology`] now says
/// [`Topology::SeparateDisks`] and `extract` refuses there, before a vertex is
/// emitted.
///
/// Three things are pinned: the configuration still produces the shape, the
/// classifier still excludes it, and **the three-step rule is no longer reached**
/// — `fan_tunnel` is not called for such a cell at all.
#[test]
fn a_nine_and_three_cell_is_refused_before_triangulation() {
    use crate::cube::is_inside as inside;
    use crate::marching_cubes::ambiguity::joined_mask;
    use crate::marching_cubes::table::AMBIGUOUS_FACES;

    // Searched, not invented: a sweep over rounded corner values looking for a
    // tunnel whose contours are not Corollary 6's shape.
    let f = [-0.8f64, 0.8, 0.6, -0.8, 0.8, -0.8, -0.2, 0.7];

    let mut case = 0u8;
    for (c, &v) in f.iter().enumerate() {
        if inside(v) {
            case |= 1 << c;
        }
    }
    assert_eq!(case, 0b0110_1001, "the fixture is no longer case 13");

    let saddles = BodySaddles::of(&f);
    assert!(saddles.has_inner_hexagon());
    let mask = joined_mask(&f, AMBIGUOUS_FACES[case as usize]);
    let contours = Contours::of(case, mask);

    let mut sizes: alloc::vec::Vec<usize> = (0..contours.count())
        .map(|r| contours.ring(r).len())
        .collect();
    sizes.sort_unstable();
    assert_eq!(
        sizes,
        alloc::vec![3, 9],
        "the fixture no longer produces the nine-and-three shape Corollary 6 excludes"
    );
    assert_eq!(
        contours.topology(&saddles),
        Topology::SeparateDisks,
        "the length bound no longer excludes the nine-and-three cell"
    );

    // What the old classification led to, kept as the record of why the bound
    // exists: sent to the tunnel rule, this cell leaves two contour edges with no
    // triangulation and therefore a hole.
    let mut tris = 0usize;
    let unresolved = contours.fan_tunnel(&saddles, &f, |_| tris += 1);
    assert_eq!(
        unresolved, 2,
        "the three-step count moved — the construction's gap has changed shape"
    );
    std::println!(
        "measured: case 13 mask {mask:#08b} gives rings {sizes:?}; classified \
         SeparateDisks. Under the old tunnel rule it would emit {tris} triangles and \
         leave {unresolved} contour edges with no rule"
    );
}

/// **The ring-count discriminator misclassifies case 13, and that is exactly
/// where the missing triangulation bites (M-229).**
///
/// `Contours::topology` calls a six-saddle cell a tunnel when it has two or more
/// contours and a twelve-vertex disk when it has one — which is what the authors'
/// implementation branches on, and what V-31 recorded as making Proposition 1's
/// asymptote-side predicate unnecessary. **It is not unnecessary.**
///
/// Measured by flood-filling the cell's inside region on a 96³ grid and counting
/// how many components its inside *corners* fall into — a computation that shares
/// nothing with the classifier:
///
/// - a genuine tunnel joins its same-signed corners through the interior, so they
///   land in **one** component. Both shipped tunnel fixtures do.
/// - the `[9,3]` case-13 cells land in **two**. Their inside region is two
///   separate blobs, which is not one cylinder — a different case-13 subcase, and
///   Grosso's Corollary 6 excludes them from the tunnel case for that reason.
///
/// And those are precisely the cells whose contour edges span three hexagon steps
/// and hit the rule the construction does not define (A-020). So the gap is a
/// **classification** failure first and a triangulation gap second, which is what
/// A-020 has to settle before it derives anything.
#[test]
fn a_nine_and_three_cell_is_not_one_connected_tunnel() {
    use crate::cube::is_inside as inside;
    use crate::marching_cubes::ambiguity::joined_mask;
    use crate::marching_cubes::table::AMBIGUOUS_FACES;

    for (label, f) in [
        (
            "[9,3] case 13",
            [-0.8f64, 0.8, 0.6, -0.8, 0.8, -0.8, -0.2, 0.7],
        ),
        (
            "[9,3] case 13 b",
            [-0.3, 0.2, 0.3, -0.4, 0.6, -0.5, -0.4, 1.0],
        ),
        (
            "[6,3] tunnel",
            [-0.2, -0.9, 0.7, -0.3, 0.8, -0.4, -0.9, 0.6],
        ),
        ("[4,3] tunnel", [0.2, -0.6, 0.2, -0.2, 0.4, 0.1, -0.9, 1.0]),
    ] {
        let mut case = 0u8;
        for (c, &v) in f.iter().enumerate() {
            if inside(v) {
                case |= 1 << c;
            }
        }
        let saddles = BodySaddles::of(&f);
        let contours = Contours::of(case, joined_mask(&f, AMBIGUOUS_FACES[case as usize]));
        let sizes: alloc::vec::Vec<usize> = (0..contours.count())
            .map(|r| contours.ring(r).len())
            .collect();
        let inside_count = f.iter().filter(|v| inside(**v)).count();
        let parts = corner_groups(&f, 96);
        let is_nine_and_three = sizes.contains(&9);
        assert_eq!(
            parts,
            if is_nine_and_three { 2 } else { 1 },
            "{label}: inside corners fall into {parts} components, which contradicts \
             the classification this test exists to check"
        );
        // The ring count called every one of these a tunnel; Corollary 6's bound
        // splits them, and it splits them the way the flood fill does.
        assert_eq!(
            contours.topology(&saddles),
            if is_nine_and_three {
                Topology::SeparateDisks
            } else {
                Topology::Tunnel
            },
            "{label}"
        );
        std::println!(
            "measured: {label:<16} case {case:#010b} rings {sizes:?} | inside corners \
             {inside_count} in {parts} component(s)"
        );
    }
}

/// **Proposition 1's asymptote-side predicate is Corollary 1, and it is not the
/// tunnel test (M-230).**
///
/// [`super::BodySaddles::same_asymptote_side`] is *derived* from the paper's
/// normal form rather than transcribed, which by rule 5 makes it a claim until
/// something independent grades it. The grader here is the contour count, which
/// shares no arithmetic with it at all: the predicate is trigonometry on the
/// quadratic's roots, the count is a walk over the face segment links.
///
/// Corollary 1 says a six-solution cell that is **not** a tunnel has a contour of
/// twelve vertices, and that its solutions are then *not* all on the same side of
/// the asymptotes. So the two must agree exactly — and they do, on every cell
/// measured, which is what makes the derivation trustworthy.
///
/// **It also shows what the predicate does not decide.** Every multi-contour
/// shape passes it, the `[9,3]` cells included, so Proposition 1 is not what
/// separates them. Corollary 6's length bound is.
#[test]
fn the_asymptote_side_predicate_is_corollary_1() {
    use crate::cube::is_inside as inside;
    use crate::marching_cubes::ambiguity::joined_mask;
    use crate::marching_cubes::table::AMBIGUOUS_FACES;

    let mut rng = Lcg::new(0x0000_A020_F000_0001);
    let mut same = 0usize;
    let mut opposite = 0usize;
    let mut disagreements = 0usize;

    for _ in 0..400_000 {
        let f = rng.corners();
        let mut case = 0u8;
        for (c, &v) in f.iter().enumerate() {
            if inside(v) {
                case |= 1 << c;
            }
        }
        if case == 0 || case == 255 {
            continue;
        }
        let saddles = BodySaddles::of(&f);
        if !saddles.has_inner_hexagon() {
            continue;
        }

        let predicate = saddles.same_asymptote_side(&f);
        let contours = Contours::of(case, joined_mask(&f, AMBIGUOUS_FACES[case as usize]));
        // Corollary 1: not all on the same side <=> the single twelve-vertex contour.
        if predicate == (contours.count() == 1) {
            disagreements += 1;
        }
        if predicate {
            same += 1;
        } else {
            opposite += 1;
        }
    }

    assert_eq!(
        disagreements, 0,
        "the derived predicate and the contour count disagree, so one of them is wrong"
    );
    assert!(same > 0, "no same-side cell was reached");
    assert!(opposite > 0, "no opposite-side cell was reached");
    std::println!(
        "measured: {same} same-side and {opposite} opposite-side six-saddle cells, \
         0 disagreements with Corollary 1"
    );
}

/// **Corollary 6's length bound is the tunnel test, and the flood fill grades it
/// (M-230).**
///
/// For a cell with exactly **two** contours the question "is this one cylinder or
/// two disks" has an unambiguous independent answer: [`corner_groups`] flood-fills
/// the inside region and counts the components its inside corners land in. A
/// cylinder joins them, so one component; two disks cap two blobs, so two.
///
/// The restriction to two contours is load-bearing and is why M-229's oracle
/// could not simply be applied everywhere. A **three**-contour tunnel has a
/// detached ring capping its own corner group, so two components there is the
/// correct answer *for a tunnel* — the census below records 14 such cells. And a
/// twelve-vertex contour has one component while being a disk, so the flood fill
/// says nothing about that split either. It is decisive only here.
#[test]
fn corollary_6s_length_bound_agrees_with_the_flood_fill() {
    use crate::cube::is_inside as inside;
    use crate::marching_cubes::ambiguity::joined_mask;
    use crate::marching_cubes::table::AMBIGUOUS_FACES;

    // **Both constants are measured, not guessed.** The oracle's resolution has a
    // floor: a narrow tunnel closes up on a coarse grid and reads as two
    // components, which would fail this test for a reason that has nothing to do
    // with the classifier. Swept against a 96³ reference over the same 400 cells,
    // 24³ and 32³ each miss one; **48³ and 64³ miss none**, as do all four on the
    // decisive fixtures. 48 is therefore the floor plus nothing, and the cell
    // count is what keeps a debug run in seconds rather than minutes — 96³ over
    // 400 cells took 415 s.
    const GRID: usize = 48;
    const CELLS: usize = 150;

    let mut rng = Lcg::new(0x0000_A020_F000_0002);
    let mut tunnels = 0usize;
    let mut seen = 0usize;
    let mut disagreements = alloc::vec::Vec::new();

    while seen < CELLS {
        let f = rng.corners();
        let mut case = 0u8;
        for (c, &v) in f.iter().enumerate() {
            if inside(v) {
                case |= 1 << c;
            }
        }
        if case == 0 || case == 255 {
            continue;
        }
        let saddles = BodySaddles::of(&f);
        if !saddles.has_inner_hexagon() {
            continue;
        }
        let contours = Contours::of(case, joined_mask(&f, AMBIGUOUS_FACES[case as usize]));
        if contours.count() != 2 {
            continue;
        }
        seen += 1;

        let is_tunnel = contours.topology(&saddles) == Topology::Tunnel;
        if is_tunnel {
            tunnels += 1;
        }
        if is_tunnel != (corner_groups(&f, GRID) == 1) {
            disagreements.push(f);
        }
    }

    assert!(
        disagreements.is_empty(),
        "the classifier and the {GRID}^3 flood fill disagree on {} of {seen} two-contour \
         six-saddle cells; first is {:?}",
        disagreements.len(),
        disagreements.first()
    );
    assert_eq!(
        tunnels, seen,
        "the random sweep is expected to reach only genuine tunnels here; the \
         nine-and-three shape is found by search, not by sampling"
    );
    std::println!(
        "measured: {seen} two-contour six-saddle cells, all classified tunnels, \
         0 disagreements with the {GRID}^3 flood fill"
    );
}

/// **What a tunnel's contours actually look like, and which half of Corollary 6
/// survives measurement (M-230).**
///
/// Corollary 6 reads *"one of the contours can have at most 6 vertices and the
/// other 3 vertices."* Only the **bound** is used by
/// [`Contours::topology`], and this is why: `[4,4]` and `[3,3,6]` tunnels are
/// common, and neither has a second contour of three vertices. The corollary's
/// second half is a description of the typical case, not a property.
///
/// The census is pinned as a whole so a new shape is as visible as a vanished
/// one, and so that [`super::MAX_TUNNEL_CONTOUR`] is a measured bound rather than
/// a transcribed one.
#[test]
fn the_tunnel_contour_shapes_are_pinned() {
    use crate::cube::is_inside as inside;
    use crate::marching_cubes::ambiguity::joined_mask;
    use crate::marching_cubes::table::AMBIGUOUS_FACES;

    let mut rng = Lcg::new(0x0000_A002_F000_0001);
    let mut shapes: alloc::collections::BTreeMap<alloc::vec::Vec<usize>, usize> =
        alloc::collections::BTreeMap::new();

    for _ in 0..400_000 {
        let f = rng.corners();
        let mut case = 0u8;
        for (c, &v) in f.iter().enumerate() {
            if inside(v) {
                case |= 1 << c;
            }
        }
        if case == 0 || case == 255 {
            continue;
        }
        let saddles = BodySaddles::of(&f);
        let contours = Contours::of(case, joined_mask(&f, AMBIGUOUS_FACES[case as usize]));
        if contours.topology(&saddles) != Topology::Tunnel {
            continue;
        }
        let mut sizes: alloc::vec::Vec<usize> = (0..contours.count())
            .map(|r| contours.ring(r).len())
            .collect();
        sizes.sort_unstable();
        *shapes.entry(sizes).or_insert(0) += 1;
    }

    for (sizes, n) in &shapes {
        std::println!("measured: tunnel contours {sizes:?} x{n}");
        assert!(
            sizes.iter().all(|&s| s <= super::MAX_TUNNEL_CONTOUR),
            "a tunnel with a contour past Corollary 6's bound reached the census"
        );
    }
    // The half of Corollary 6 that measurement falsifies: a second contour of
    // three vertices is not a property of tunnels.
    assert!(
        shapes.keys().any(|s| !s.contains(&3)),
        "no tunnel without a three-vertex contour was reached, so the falsification \
         of Corollary 6's second half is untested"
    );
}

/// **Every separate-disks cell has a body saddle sitting on a cell face, and
/// continuous corner values produce none at all (M-231).**
///
/// This is A-020b's premise falling over, and it is the reason that ticket asks
/// for no new triangulation. [`Topology::SeparateDisks`] looked like a topological
/// subcase of case 13 that Grosso's construction forgot. It is not. It is the
/// numerical signature of a **singular face** — a saddle of the trilinear
/// interpolant lying *on* a face of the cell, which is Grosso 2017 §4.2's case and
/// this repository's A-002i — surviving [`BodySaddles::has_inner_hexagon`]'s
/// strict `0 < x < 1` test because floating point puts the root a few ulps inside.
///
/// Two measurements say so, and the second is what makes the first mean something:
///
/// - **Continuous corner values produce zero** separate-disks cells. The shipped
///   fixtures were found by a search over *rounded* values, and rounding is
///   precisely what makes a root land on a face.
/// - **Every one that quantisation does produce has a root within `1e-12` of a
///   face** — with no exceptions across every quantum tried — while the same
///   degeneracy among other six-saddle cells runs at a rate that swings with the
///   quantum. The invariant is not "rounding makes cells degenerate"; it is this
///   configuration specifically.
///
/// **Not proof, and the test is written to say so.** Zero over a sample is a
/// sample, so the continuous arm asserts zero and would fail loudly if a
/// non-degenerate one ever appeared — which is exactly the event that would put
/// A-020b back in scope.
#[test]
fn every_separate_disks_cell_has_a_saddle_on_a_face() {
    use crate::cube::is_inside as inside;
    use crate::marching_cubes::ambiguity::joined_mask;
    use crate::marching_cubes::table::AMBIGUOUS_FACES;

    // A root this close to 0 or 1 is a saddle on the face; the strict interior
    // test admits it because the arithmetic misses the face by a few ulps.
    const ON_FACE: f64 = 1e-12;
    const CELLS: usize = 500_000;

    let classify = |f: &[f64; 8]| {
        let mut case = 0u8;
        for (c, &val) in f.iter().enumerate() {
            if inside(val) {
                case |= 1 << c;
            }
        }
        if case == 0 || case == 255 {
            return None;
        }
        let saddles = BodySaddles::of(f);
        if !saddles.has_inner_hexagon() {
            return None;
        }
        let contours = Contours::of(case, joined_mask(f, AMBIGUOUS_FACES[case as usize]));
        Some((saddles, contours.topology(&saddles)))
    };
    let on_face = |saddles: &BodySaddles<f64>| {
        (0..3).any(|axis| saddles.axis(axis).iter().any(|&x| x.min(1.0 - x) < ON_FACE))
    };

    // Continuous: none exist.
    let mut rng = Lcg::new(0x0000_A020_B000_0003);
    let mut six = 0usize;
    let mut separate = 0usize;
    for _ in 0..CELLS {
        let f = rng.corners();
        if let Some((_, topology)) = classify(&f) {
            six += 1;
            if topology == Topology::SeparateDisks {
                separate += 1;
            }
        }
    }
    assert!(
        six > 1000,
        "only {six} six-saddle cells, too few to mean anything"
    );
    assert_eq!(
        separate, 0,
        "a separate-disks cell appeared without quantisation — A-020b's premise is \
         back in scope and the singular-face explanation is not the whole story"
    );
    std::println!(
        "measured: {six} six-saddle cells from continuous values, {separate} separate-disks"
    );

    // Quantised: they exist, and every one has a saddle on a face.
    for (label, quantum) in [("0.1", 0.1f64), ("1/3", 1.0 / 3.0)] {
        let mut rng = Lcg::new(0x0000_A020_B000_0003);
        let mut found = 0usize;
        let mut degenerate = 0usize;
        for _ in 0..CELLS {
            let mut f = rng.corners();
            for x in &mut f {
                *x = (*x / quantum).round() * quantum;
            }
            if let Some((saddles, topology)) = classify(&f)
                && topology == Topology::SeparateDisks
            {
                found += 1;
                if on_face(&saddles) {
                    degenerate += 1;
                }
            }
        }
        assert!(
            found > 0,
            "quantum {label} produced no separate-disks cell, so this asserted nothing"
        );
        assert_eq!(
            found,
            degenerate,
            "quantum {label}: {} of {found} separate-disks cells have no saddle on a face, \
             so the singular-face explanation does not cover them",
            found - degenerate
        );
        std::println!(
            "measured: quantum {label} gives {found} separate-disks cells, all {degenerate} \
             with a saddle on a cell face"
        );
    }
}
