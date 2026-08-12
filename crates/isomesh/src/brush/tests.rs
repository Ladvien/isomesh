//! The load-bearing test is the acceptance criterion: eight brushes in all
//! `8! = 40,320` orderings, counting distinct results.

use alloc::vec::Vec;

use super::{Brush, BrushOp, BrushStack, Capsule, smooth_min};
use crate::fields::{BoxExact, Sphere};
use crate::{Real, Sdf};

/// A brush shape, as one type so a stack can hold a mixture.
#[derive(Clone, Copy, Debug)]
enum Shape {
    Sphere(Sphere<f64>),
    Cube(BoxExact<f64>),
    Capsule(Capsule<f64>),
}

impl Sdf for Shape {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Sphere(s) => s.sample(p),
            Self::Cube(b) => b.sample(p),
            Self::Capsule(c) => c.sample(p),
        }
    }
}

/// Eight brushes, all three shapes represented, deliberately overlapping so the
/// order could matter.
fn eight(op_of: impl Fn(usize) -> BrushOp) -> Vec<Brush<Shape>> {
    let centres = [
        [0.30, 0.10, -0.20],
        [-0.25, 0.35, 0.15],
        [0.05, -0.30, 0.25],
        [-0.15, -0.10, -0.35],
        [0.40, 0.25, 0.05],
        [-0.35, 0.05, 0.30],
        [0.20, -0.40, -0.10],
        [-0.05, 0.20, -0.30],
    ];
    centres
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let shape = match i % 3 {
                0 => Shape::Sphere(Sphere {
                    center: *c,
                    radius: 0.30 + 0.02 * i as f64,
                }),
                1 => Shape::Cube(BoxExact {
                    center: *c,
                    half_extents: [0.22 + 0.01 * i as f64; 3],
                }),
                _ => Shape::Capsule(Capsule {
                    a: *c,
                    b: [c[0] + 0.25, c[1] - 0.15, c[2] + 0.1],
                    radius: 0.16,
                }),
            };
            Brush {
                shape,
                op: op_of(i),
            }
        })
        .collect()
}

/// Sample points the stacks are compared at.
///
/// A fixed lattice rather than random, so the comparison is reproducible, and
/// spread across the region the brushes occupy so a difference anywhere shows up.
fn probes() -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..4 {
                out.push([
                    -0.6 + 0.4 * f64::from(x),
                    -0.6 + 0.4 * f64::from(y),
                    -0.6 + 0.4 * f64::from(z),
                ]);
            }
        }
    }
    out
}

/// The field's values at the probes, as raw bits, so `+0.0` and `-0.0` are
/// distinct and two stacks count as equal only if they agree exactly.
fn signature(brushes: &[Brush<Shape>], probes: &[[f64; 3]]) -> Vec<u64> {
    let stack = BrushStack {
        base: BoxExact::<f64>::canonical(),
        brushes,
    };
    probes.iter().map(|p| stack.sample(*p).to_bits()).collect()
}

/// Every permutation of `0..n`, by Heap's algorithm.
///
/// Generated rather than sampled: the acceptance criterion says all 40,320, and
/// "we tried a thousand random orderings" is a different and weaker claim.
fn permutations(n: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut items: Vec<usize> = (0..n).collect();
    let mut counters = alloc::vec![0usize; n];
    out.push(items.clone());
    let mut i = 0;
    while i < n {
        if counters[i] < i {
            if i % 2 == 0 {
                items.swap(0, i);
            } else {
                items.swap(counters[i], i);
            }
            out.push(items.clone());
            counters[i] += 1;
            i = 0;
        } else {
            counters[i] = 0;
            i += 1;
        }
    }
    out
}

/// How many distinct results the orderings produce.
fn distinct_results(brushes: &[Brush<Shape>]) -> usize {
    let probes = probes();
    let orders = permutations(brushes.len());
    let mut seen: Vec<Vec<u64>> = Vec::new();
    for order in &orders {
        let permuted: Vec<Brush<Shape>> = order.iter().map(|&i| brushes[i]).collect();
        let sig = signature(&permuted, &probes);
        if let Err(at) = seen.binary_search(&sig) {
            seen.insert(at, sig);
        }
    }
    seen.len()
}

#[test]
fn there_really_are_40320_orderings() {
    let orders = permutations(8);
    assert_eq!(orders.len(), 40_320);
    let mut sorted = orders.clone();
    sorted.sort();
    sorted.dedup();
    assert_eq!(sorted.len(), 40_320, "permutations must be distinct");
}

/// **The acceptance criterion.** Eight adds, all 40,320 orderings, one result.
///
/// `min` is commutative *and* associative in IEEE, and it introduces no rounding
/// at all — it selects one of its arguments rather than computing a new value —
/// so this is exact rather than nearly so. Concurrent clients may apply a run of
/// additions in any order and converge on the same solid, bit for bit.
#[test]
fn eight_additions_commute_exactly_over_all_orderings() {
    let brushes = eight(|_| BrushOp::Add);
    let distinct = distinct_results(&brushes);
    std::println!(
        "measured: G-003 eight Add brushes over 40,320 orderings -> {distinct} distinct result(s)"
    );
    assert_eq!(
        distinct, 1,
        "additions must commute, or concurrent editing dies here"
    );
}

/// The same for subtraction, which is `max` and has the same algebra.
#[test]
fn eight_subtractions_commute_exactly_over_all_orderings() {
    let brushes = eight(|_| BrushOp::Subtract);
    let distinct = distinct_results(&brushes);
    std::println!(
        "measured: G-003 eight Subtract brushes over 40,320 orderings -> {distinct} distinct result(s)"
    );
    assert_eq!(distinct, 1);
}

/// **Mixed adds and subtracts do not commute**, and the acceptance criterion's
/// "expect 1" does not apply to them.
///
/// This is geometry rather than floating point: carving a hole and then filling
/// it is a different solid from filling and then carving. No storage format or
/// arithmetic fixes it, and a concurrent-editing protocol has to order edits
/// across an add/subtract boundary rather than merging them freely.
///
/// The count is reported rather than pinned — it is a property of these eight
/// brushes, not a constant of the algebra — but it must be greater than one, or
/// the fixture is not overlapping enough to be testing anything.
#[test]
fn mixing_adds_and_subtracts_does_not_commute() {
    let brushes = eight(|i| {
        if i % 2 == 0 {
            BrushOp::Add
        } else {
            BrushOp::Subtract
        }
    });
    let distinct = distinct_results(&brushes);
    std::println!(
        "measured: G-003 four Add + four Subtract over 40,320 orderings -> {distinct} distinct results"
    );
    assert!(
        distinct > 1,
        "if these commuted the fixture would not be overlapping enough to test anything"
    );
    assert!(
        !BrushOp::Add.commutes_with(BrushOp::Subtract),
        "the API must not claim they commute"
    );
}

/// **Smooth union does not commute either**, and for a different reason.
///
/// Smooth-min is commutative in its two arguments but **not associative**, so a
/// stack of smooth adds depends on the order it was folded in even though every
/// operation is the same kind. That is worth separating from the mixed case: it
/// means "all the same operation" is not sufficient for reordering, only "all the
/// same *hard* operation".
#[test]
fn smooth_union_does_not_commute() {
    let brushes = eight(|_| BrushOp::SmoothAdd { k: 0.15 });
    let distinct = distinct_results(&brushes);
    std::println!(
        "measured: G-003 eight SmoothAdd(k=0.15) over 40,320 orderings -> {distinct} distinct results"
    );
    assert!(
        distinct > 1,
        "smooth-min is not associative; this should differ"
    );
    assert!(!BrushOp::SmoothAdd { k: 0.15 }.commutes_with(BrushOp::SmoothAdd { k: 0.15 }));
}

/// Smooth-min fails reordering in **two** independent ways, and separating them
/// matters because they have different fixes.
///
/// It is not associative in exact arithmetic, which no storage format repairs.
/// And it is not *bit*-commutative — swapping the arguments evaluates a
/// different expression that agrees only to rounding, which a fixed-point
/// representation would repair. This module's docs claimed a flat "commutative"
/// until this test disagreed.
#[test]
fn smooth_min_fails_reordering_in_two_separate_ways() {
    let cases: [(f64, f64, f64); 4] = [
        (0.3, -0.2, 0.25),
        (1.0, 1.0, 0.5),
        (-0.75, 0.4, 0.1),
        (0.05, 0.06, 0.5),
    ];
    // Commutative to rounding, not to the bit.
    let mut worst_ulp = 0i64;
    for (a, b, k) in cases {
        let forward = smooth_min(a, b, k);
        let reverse = smooth_min(b, a, k);
        assert!(
            (forward - reverse).abs() <= 4.0 * f64::EPSILON * forward.abs().max(1.0),
            "smooth_min({a}, {b}, {k}) differs by more than rounding: {forward} vs {reverse}"
        );
        let gap = (forward.to_bits() as i64 - reverse.to_bits() as i64).abs();
        worst_ulp = worst_ulp.max(gap);
    }
    std::println!(
        "measured: G-003 smooth-min commutativity -- worst disagreement {worst_ulp} ulp when the arguments are swapped"
    );
    assert!(
        worst_ulp > 0,
        "if swapping is bit-exact, the module docs should say so plainly"
    );

    // And the associativity failure, which is what breaks reordering.
    //
    // The values are **close together relative to `k`**, and that is essential
    // rather than incidental: once `|a − b| ≥ k`, `h` saturates, the `k·h·(1−h)`
    // term vanishes and smooth-min degenerates to an ordinary `min` — which *is*
    // associative. A fixture chosen for looking irregular rather than for
    // exercising the smooth region passes this test while proving nothing, which
    // is G-001's method rule arriving a second time. These come from a search.
    let (a, b, c, k): (f64, f64, f64, f64) = (-0.15, 0.2, 0.2, 0.5);
    let left = smooth_min(smooth_min(a, b, k), c, k);
    let right = smooth_min(a, smooth_min(b, c, k), k);
    assert_ne!(
        left.to_bits(),
        right.to_bits(),
        "smooth-min is not associative; if this passes the claim above is wrong"
    );
    // Not a rounding difference: this is 14 orders of magnitude above an ulp,
    // which is why no storage format repairs it.
    let gap = (left - right).abs();
    std::println!(
        "measured: G-003 smooth-min associativity gap -> {left} vs {right}, {gap:.3e} apart"
    );
    assert!(
        gap > 1e-3,
        "the fixture has drifted into the saturated region where smooth-min is just min"
    );
}

/// A zero join width must degenerate to an ordinary `min` rather than dividing
/// by zero — and then it commutes and associates like one.
#[test]
fn a_zero_join_width_is_an_ordinary_min() {
    for (a, b) in [(0.3f64, -0.2f64), (1.0, 1.0), (-0.75, 0.4)] {
        assert_eq!(smooth_min(a, b, 0.0).to_bits(), a.min(b).to_bits());
        assert_eq!(smooth_min(a, b, -1.0).to_bits(), a.min(b).to_bits());
    }
}

// ─── the shapes ─────────────────────────────────────────────────────────────

/// A capsule is an exact distance field, which is what makes it usable as a
/// brush: `|∇f| = 1`, so combining it does not distort its neighbourhood.
#[test]
fn the_capsule_is_an_exact_distance_field() {
    let c = Capsule {
        a: [-0.4, 0.0, 0.0],
        b: [0.4, 0.0, 0.0],
        radius: 0.25,
    };
    // On the segment's axis, distance is radius-relative and exact.
    assert!((c.sample([0.0, 0.0, 0.0]) + 0.25).abs() < 1e-15);
    assert!((c.sample([0.0, 0.25, 0.0])).abs() < 1e-15);
    // Past an end cap it is the distance to that end, so the capsule is a
    // stadium rather than an infinite cylinder.
    assert!((c.sample([0.9, 0.0, 0.0]) - 0.25).abs() < 1e-15);

    // Unit gradient, by central difference, away from the axis singularity.
    let eps = 1e-6;
    for p in [[0.0, 0.5, 0.1], [0.7, 0.2, -0.3], [-0.6, -0.4, 0.2]] {
        let mut grad = [0.0f64; 3];
        for axis in 0..3 {
            let mut lo = p;
            let mut hi = p;
            lo[axis] -= eps;
            hi[axis] += eps;
            grad[axis] = (c.sample(hi) - c.sample(lo)) / (2.0 * eps);
        }
        let len = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-6, "|grad| = {len} at {p:?}");
    }
}

/// A zero-length capsule is a sphere, which is the right answer rather than a
/// case to reject.
#[test]
fn a_degenerate_capsule_is_a_sphere() {
    let c = Capsule {
        a: [0.1, 0.2, 0.3],
        b: [0.1, 0.2, 0.3],
        radius: 0.4,
    };
    let s = Sphere {
        center: [0.1, 0.2, 0.3],
        radius: 0.4,
    };
    for p in [[0.0; 3], [1.0, 0.0, 0.0], [-0.5, 0.5, 0.25]] {
        assert!((c.sample(p) - s.sample(p)).abs() < 1e-15, "at {p:?}");
    }
}

/// The stack really applies its brushes, and in order.
#[test]
fn a_stack_applies_its_brushes() {
    let inside = [0.0, 0.0, 0.0];
    let base = BoxExact::<f64>::canonical();
    assert!(
        base.sample(inside) < 0.0,
        "the origin starts inside the box"
    );

    let carve = [Brush::subtract(Shape::Sphere(Sphere {
        center: [0.0; 3],
        radius: 0.5,
    }))];
    let stack = BrushStack {
        base,
        brushes: &carve,
    };
    assert!(
        stack.sample(inside) > 0.0,
        "the carve should have removed it"
    );

    let refill = [
        carve[0],
        Brush::add(Shape::Sphere(Sphere {
            center: [0.0; 3],
            radius: 0.5,
        })),
    ];
    let stack = BrushStack {
        base,
        brushes: &refill,
    };
    assert!(
        stack.sample(inside) < 0.0,
        "adding it back should restore it"
    );
}
