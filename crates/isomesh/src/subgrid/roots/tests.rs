//! §4.3.2's tests.
//!
//! The one that matters is
//! `a_slab_thinner_than_the_grid_still_has_two_roots`: it is the whole premise
//! of the subgrid track, measured directly rather than inferred from a mesh.

use super::*;
use alloc::vec::Vec;

/// A field defined by a closure, so a test can state the function it means.
struct Field<F>(F);

impl<F: Fn([f64; 3]) -> f64> Sdf for Field<F> {
    type Scalar = f64;
    fn sample(&self, p: [f64; 3]) -> f64 {
        (self.0)(p)
    }
}

/// Roots along the unit segment on the x axis.
fn roots_of(f: impl Fn(f64) -> f64, samples: u32) -> Vec<f64> {
    let mut out = Vec::new();
    all_roots(
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        &Field(|p: [f64; 3]| f(p[0])),
        samples,
        &mut out,
    );
    out
}

#[test]
fn a_monotone_crossing_gives_one_root_at_the_right_place() {
    // f(x) = x - 1/3, so the root is at t = 1/3.
    let roots = roots_of(|x| x - 1.0 / 3.0, 8);
    assert_eq!(roots.len(), 1);
    assert!(
        (roots[0] - 1.0 / 3.0).abs() < 1e-15,
        "root at {} rather than 1/3",
        roots[0]
    );
}

#[test]
fn no_sign_change_gives_no_roots() {
    assert!(roots_of(|_| 1.0, 16).is_empty());
    assert!(roots_of(|_| -1.0, 16).is_empty());
}

#[test]
fn a_slab_thinner_than_the_grid_still_has_two_roots() {
    // The premise of the whole subgrid track, stated at its simplest. A slab of
    // half-width w centred at 1/2 is inside on (1/2 - w, 1/2 + w) and outside
    // elsewhere, so the edge carries exactly two roots however thin it gets --
    // as long as the 1D sampling resolves it, which is the bound §1.3 states.
    //
    // A sign test at the endpoints sees "outside, outside" and reports nothing,
    // which is A-005's zero-triangle result at its root.
    for w in [0.25, 0.05, 0.01, 0.001] {
        let slab = move |x: f64| (x - 0.5).abs() - w;
        // Endpoints are both outside: a sign test learns nothing here.
        assert!(slab(0.0) > 0.0 && slab(1.0) > 0.0);

        // Enough samples to bracket a slab of this width.
        let samples = (4.0 / w) as u32;
        let roots = roots_of(slab, samples);
        assert_eq!(roots.len(), 2, "w = {w} gave {roots:?}");
        assert!((roots[0] - (0.5 - w)).abs() < 1e-12, "w = {w}: {roots:?}");
        assert!((roots[1] - (0.5 + w)).abs() < 1e-12, "w = {w}: {roots:?}");
    }
}

#[test]
fn a_slab_thinner_than_the_sampling_is_missed_and_that_is_the_stated_bound() {
    // §1.3: "1D marching can of course miss intersections, [but] we are no worse
    // off than classic marching." This is that limitation, asserted rather than
    // left implicit -- the failure is silent, so it belongs in a test where it
    // can be read.
    // The centre is deliberately *not* on the sample lattice. Centring it at
    // 0.5 makes this test pass for the wrong reason: with 8 samples one lands
    // exactly on 0.5, reads inside, and the slab is found after all. Fourth
    // occurrence of the fixture trap in this repo (M-32, M-38, M-44), and the
    // first where the fixture was mine.
    let slab = |x: f64| (x - 0.5137).abs() - 1e-4;
    assert!(
        roots_of(slab, 8).is_empty(),
        "8 samples should step over a slab 2e-4 wide"
    );
    // And it is a resolution question, not a defect: enough samples find it.
    assert_eq!(roots_of(slab, 100_000).len(), 2);
}

#[test]
fn roots_come_back_ascending_and_distinct() {
    // Four crossings from a quartic with roots at 0.2, 0.4, 0.6, 0.8.
    let f = |x: f64| (x - 0.2) * (x - 0.4) * (x - 0.6) * (x - 0.8);
    let roots = roots_of(f, 32);
    assert_eq!(roots.len(), 4);
    for w in roots.windows(2) {
        assert!(w[0] < w[1], "not ascending: {roots:?}");
    }
    for (found, expected) in roots.iter().zip([0.2, 0.4, 0.6, 0.8]) {
        assert!((found - expected).abs() < 1e-12, "{roots:?}");
    }
}

#[test]
fn zero_counts_as_outside_so_a_touch_is_not_a_crossing() {
    // f(x) = (x - 1/2)^2 touches zero and does not pass through it. Inside is
    // `f < 0`, so nothing here is ever inside and there is no crossing -- the
    // same convention marching cubes and marching tetrahedra use, and the place
    // an inconsistency in it would first show.
    assert!(roots_of(|x| (x - 0.5) * (x - 0.5), 64).is_empty());

    // The mirrored case does cross, twice.
    assert_eq!(roots_of(|x| -((x - 0.5) * (x - 0.5)) + 0.04, 64).len(), 2);
}

#[test]
fn the_same_edge_from_either_direction_is_a_caller_obligation_not_a_promise() {
    // Walking an edge backwards gives the mirrored parameters, and they are the
    // same points in space -- but only to within the arithmetic. This asserts
    // what is actually true, because the docs promise determinism for identical
    // arguments and nothing more, and a test claiming symmetry would be claiming
    // something the extractor must not rely on.
    let f = |p: [f64; 3]| p[0] - 1.0 / 3.0;
    let mut forward = Vec::new();
    all_roots([0.0; 3], [1.0, 0.0, 0.0], &Field(f), 8, &mut forward);
    let mut backward = Vec::new();
    all_roots([1.0, 0.0, 0.0], [0.0; 3], &Field(f), 8, &mut backward);

    assert_eq!(forward.len(), 1);
    assert_eq!(backward.len(), 1);
    let forward_x = forward[0];
    let backward_x = 1.0 - backward[0];
    assert!(
        (forward_x - backward_x).abs() < 1e-15,
        "{forward_x} vs {backward_x}"
    );
}

#[test]
fn the_same_call_twice_is_bit_identical() {
    // Two tetrahedra sharing an edge call this independently, so anything less
    // than bit-equality would put a crack in the mesh that no combinatorial
    // guarantee could close.
    let f = |p: [f64; 3]| (p[0] - 0.3) * (p[0] - 0.7);
    let mut a = Vec::new();
    let mut b = Vec::new();
    all_roots([0.0; 3], [1.0, 0.0, 0.0], &Field(f), 16, &mut a);
    all_roots([0.0; 3], [1.0, 0.0, 0.0], &Field(f), 16, &mut b);
    assert_eq!(a, b);
    assert_eq!(a.len(), 2);
}

// The sentinel is compared exactly on purpose: the point is that `all_roots`
// left it *untouched*, which is bit-equality and not approximation.
#[expect(clippy::float_cmp, reason = "asserting a value was not written to")]
#[test]
fn results_append_rather_than_replace() {
    // The buffer is caller-provided and reused per edge, per CLAUDE.md rule 6.
    let mut out = alloc::vec![-1.0f64];
    all_roots(
        [0.0; 3],
        [1.0, 0.0, 0.0],
        &Field(|p: [f64; 3]| p[0] - 0.5),
        4,
        &mut out,
    );
    assert_eq!(out.len(), 2);
    assert_eq!(out[0], -1.0);
}

#[test]
fn zero_samples_finds_nothing_rather_than_dividing_by_zero() {
    assert!(roots_of(|x| x - 0.5, 0).is_empty());
}
