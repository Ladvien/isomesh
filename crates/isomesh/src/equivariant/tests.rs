//! The reductions do what the module claims, including the parts that are
//! easiest to assume rather than check.

use super::{sort_by_magnitude, sum_equivariant};

/// Compare sorted output by bit pattern: these tests are *about* which exact
/// value landed where, including which zero, so `==` on the floats would both
/// trip `float_cmp` and be the weaker check.
fn bits<const N: usize>(t: [f64; N]) -> [u64; N] {
    t.map(f64::to_bits)
}

/// Three terms chosen so the naive left-to-right sum is *wrong*: the two large
/// ones cancel, and adding the small one first is what preserves it.
const CANCELLING: [f64; 3] = [1.0e-16, 1.0, -1.0];

#[test]
fn the_sum_does_not_depend_on_the_order_the_terms_arrive_in() {
    let mut seen = None;
    for permutation in [
        [CANCELLING[0], CANCELLING[1], CANCELLING[2]],
        [CANCELLING[0], CANCELLING[2], CANCELLING[1]],
        [CANCELLING[1], CANCELLING[0], CANCELLING[2]],
        [CANCELLING[1], CANCELLING[2], CANCELLING[0]],
        [CANCELLING[2], CANCELLING[0], CANCELLING[1]],
        [CANCELLING[2], CANCELLING[1], CANCELLING[0]],
    ] {
        let got = sum_equivariant(permutation);
        match seen {
            None => seen = Some(got),
            Some(first) => assert_eq!(
                got.to_bits(),
                first.to_bits(),
                "permutation {permutation:?} summed to {got:e}, not {first:e}"
            ),
        }
    }
}

/// The negative control. Without this the test above could pass on a reduction
/// that is order-independent because the inputs are too tame to tell, and the
/// whole exercise would be decoration.
#[test]
fn an_unsorted_sum_of_the_same_terms_really_does_disagree() {
    let naive = |t: [f64; 3]| t[0] + t[1] + t[2];
    assert_ne!(
        naive([CANCELLING[0], CANCELLING[1], CANCELLING[2]]).to_bits(),
        naive([CANCELLING[1], CANCELLING[2], CANCELLING[0]]).to_bits(),
        "the fixture cannot detect an accumulation-order defect"
    );
}

/// The second negative control, and the one that found M-175: sorting on the
/// magnitudes *alone* is not enough. `+1` and `−1` tie, so a stable sort leaves
/// their order to the order they arrived in — which is precisely the labelling a
/// lattice rotation permutes.
#[test]
fn sorting_on_magnitude_alone_would_still_depend_on_the_arrival_order() {
    let sum = |t: [f64; 3]| t.iter().fold(0.0f64, |a, b| a + b);
    // Both are correctly sorted smallest-magnitude-first. They disagree.
    assert_ne!(
        sum([1.0e-16, 1.0, -1.0]).to_bits(),
        sum([1.0e-16, -1.0, 1.0]).to_bits(),
        "a magnitude-only sort cannot be what makes the reduction equivariant"
    );
    // The tie-break picks the second of those unconditionally: −c before +c.
    assert_eq!(
        sum_equivariant(CANCELLING).to_bits(),
        sum([1.0e-16, -1.0, 1.0]).to_bits()
    );
}

/// The padding claim the callers rely on: a twelve-slot reduction over three
/// real terms must give the same bits as a three-slot one, or every cell with
/// fewer than twelve crossings would answer differently for no reason.
#[test]
fn zero_padding_a_reduction_does_not_change_a_single_bit() {
    let mut padded = [0.0f64; 12];
    padded[4] = CANCELLING[0];
    padded[9] = CANCELLING[1];
    padded[1] = CANCELLING[2];
    assert_eq!(
        sum_equivariant(padded).to_bits(),
        sum_equivariant(CANCELLING).to_bits()
    );
}

/// Negative zero is the obvious place to expect padding to leak, and it does
/// not: the accumulator seeds at `+0.0` and `+0.0 + (−0.0)` is `+0.0`, so an
/// all-negative-zero reduction returns `+0.0` either way (M-176). Pinned because
/// the claim was made in the opposite direction before it was checked.
#[test]
fn negative_zero_survives_padding_because_the_accumulator_starts_positive() {
    assert_eq!(sum_equivariant([-0.0f64; 3]).to_bits(), 0.0f64.to_bits());

    let mut padded = [0.0f64; 12];
    padded[0] = -0.0;
    padded[5] = -0.0;
    padded[11] = -0.0;
    assert_eq!(
        sum_equivariant(padded).to_bits(),
        sum_equivariant([-0.0f64; 3]).to_bits()
    );
}

#[test]
fn sorting_is_by_magnitude_and_keeps_signs() {
    let mut t = [3.0f64, -1.0, 2.0, -0.5];
    sort_by_magnitude(&mut t);
    assert_eq!(bits(t), bits([-0.5, -1.0, 2.0, 3.0]));
}

/// The tie-break itself: equal magnitudes come out negative-first, from any
/// starting arrangement.
#[test]
fn equal_magnitudes_come_out_negative_first() {
    for start in [[1.0f64, -1.0], [-1.0, 1.0]] {
        let mut t = start;
        sort_by_magnitude(&mut t);
        assert_eq!(bits(t), bits([-1.0, 1.0]), "from {start:?}");
    }
}

/// **The limit of what reordering can buy, pinned as a fact rather than left as
/// an assumption (M-177).** Permutation invariance is achievable; negation
/// *equivariance* — `φ(−S) = −φ(S)` — is not, by any ordering rule.
///
/// The obstruction is structural. A magnitude tie group holds `m` copies of `−c`
/// and `n` of `+c`; negating swaps those counts, so no ordering rule that is a
/// function of the multiset can map the group's order onto its own reverse
/// unless `m == n`. Here is a witness, and it is not a corner case: the sum of
/// `[1e-16, +1, −1]` is `1.11e-16`, and of its negation `0`.
///
/// This matters because a lattice rotation *can* negate a component, so a sum of
/// position or normal components is not bit-exactly equivariant under the full
/// octahedral group by this route. The three-term dot product is unaffected and
/// that is why ✗12 measured 0/9600: `(Ra)ᵢ(Rb)ᵢ = aⱼbⱼ` for a signed permutation
/// — the two sign flips cancel inside each product, so its terms permute without
/// negating.
#[test]
fn negation_equivariance_is_not_achievable_by_ordering_and_here_is_the_witness() {
    let terms = [1.0e-16f64, 1.0, -1.0];
    let negated = terms.map(|t| -t);
    assert_ne!(
        sum_equivariant(negated).to_bits(),
        (-sum_equivariant(terms)).to_bits(),
        "if this ever passes, revisit A-016's successor -- the obstruction moved"
    );

    // What *does* hold, and is what the callers need: both are functions of
    // their own multiset, whichever order the terms are presented in.
    for permuted in [[1.0f64, -1.0, 1.0e-16], [-1.0, 1.0e-16, 1.0]] {
        assert_eq!(
            sum_equivariant(permuted).to_bits(),
            sum_equivariant(terms).to_bits()
        );
    }
}

/// `total_cmp` rather than `<` is load-bearing: a NaN makes every `<`
/// comparison false, which leaves an insertion sort's output dependent on the
/// starting arrangement. The sum is NaN either way — what must not vary is
/// whether the *sort* terminated somewhere reproducible.
#[test]
fn a_nan_does_not_make_the_order_depend_on_the_arrangement() {
    let mut a = [f64::NAN, 1.0, -2.0];
    let mut b = [-2.0, f64::NAN, 1.0];
    sort_by_magnitude(&mut a);
    sort_by_magnitude(&mut b);
    assert_eq!(bits(a), bits(b));
}

#[test]
fn f32_reduces_the_same_way() {
    let terms = [1.0e-8f32, 1.0, -1.0];
    assert_eq!(
        sum_equivariant([terms[2], terms[0], terms[1]]).to_bits(),
        sum_equivariant(terms).to_bits()
    );
}
