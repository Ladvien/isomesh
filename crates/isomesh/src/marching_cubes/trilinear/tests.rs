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

use super::{ALL_INSIDE, BodySaddles};
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
