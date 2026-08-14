//! Reductions whose result is a function of the *set* of terms, not of the
//! order they were written in.
//!
//! # Why this exists
//!
//! IEEE addition and multiplication are **commutative but not associative**:
//! `(a + b) + c` and `(b + c) + a` differ in the last bits. That is harmless
//! until something permutes the order, and in this crate something does — a
//! lattice rotation permutes the axis labels of a vector and the edge labels of
//! a cell, so a running sum taken in label order is a function of *how the cell
//! is oriented on the grid* rather than of the cell.
//!
//! The audit measured it directly on the three-term case: **4328 of 9600**
//! lattice-symmetry trials disagreed with an unsorted dot product, and **0 of
//! 9600** with a sorted one (✗12).
//!
//! Sorting by magnitude first fixes it, because magnitude order is invariant
//! under permutation of the terms. The sum then depends only on which values are
//! present.
//!
//! # Why the reductions are fixed-size
//!
//! Every caller reduces over a slot per *label* — three axes, or twelve cell
//! edges — with absent labels left at `R::ZERO`. That is deliberate rather than
//! wasteful: padding with `+0.0` sorts first and adds exactly, so a twelve-slot
//! reduction over three real terms gives the same bits as a three-slot one.
//!
//! **Including for negative zero**, which is worth stating because it is the
//! obvious place to expect an exception and there isn't one. `+0.0 + (−0.0)` is
//! `+0.0` under round-to-nearest, and [`sum_equivariant`] seeds its accumulator
//! at `R::ZERO`, so a reduction whose every term is `−0.0` returns `+0.0`
//! whether it was padded or not — the padding changes nothing a golden hash
//! could see (M-176).

use crate::Real;

/// The order the reductions sum in: ascending by magnitude, ties broken by the
/// signed value.
///
/// `total_cmp` is a total order even across NaN and signed zero, so the network
/// cannot produce an order that depends on which comparison ran first.
///
/// # The tie-break is load-bearing, and the comment this replaced said it was not
///
/// The original form sorted on the absolute values alone and argued that ties
/// need no tie-break, *"since IEEE addition and multiplication are both
/// commutative — only associativity fails — so two terms of equal magnitude give
/// the same answer either way round."* That is true of `a + b` against `b + a`
/// and false of everything longer, which is the only case that occurs here.
/// Measured (M-175): summing `[1e-16, +1, −1]` smallest-first gives `0`, and
/// `[1e-16, −1, +1]` gives `1.11e-16` — the pair ties in magnitude, so a *stable*
/// sort left their order to the order they arrived in, which for this crate is
/// the axis or edge labelling that a lattice rotation permutes.
///
/// Comparing the signed values on a tie removes it: `−c` sorts before `+c` for
/// any `c`, so the sequence is a function of the multiset of terms alone. It
/// survives negation, which is what a rotation can do to a component: negating
/// every term maps each magnitude group onto itself (a group is either repeats
/// of one value, or a `±c` pair), and `−(a + b)` is exactly `(−a) + (−b)`.
///
/// Insertion sort, because `N` is 3, 6 or 12 here and it is branch-predictable
/// and allocation-free, which a comparison sort over a slice would not be in
/// `no_std`. At the largest of those it is 66 comparisons worst case, against a
/// reduction that would otherwise be 11 additions — the cost is real and it buys
/// the only property that makes the result a function of the cell.
#[inline]
pub(crate) fn sort_by_magnitude<R: Real, const N: usize>(t: &mut [R; N]) {
    let mut i = 1;
    while i < N {
        let mut j = i;
        while j > 0 && precedes(t[j], t[j - 1]) {
            t.swap(j - 1, j);
            j -= 1;
        }
        i += 1;
    }
}

/// Whether `a` sums before `b`: smaller magnitude first, then smaller value.
#[inline]
fn precedes<R: Real>(a: R, b: R) -> bool {
    match a.abs().total_cmp(&b.abs()) {
        core::cmp::Ordering::Less => true,
        core::cmp::Ordering::Greater => false,
        core::cmp::Ordering::Equal => a.total_cmp(&b) == core::cmp::Ordering::Less,
    }
}

/// Sum smallest-magnitude-first, so the result is a function of the *set* of
/// terms rather than of the order they were written in.
#[inline]
pub(crate) fn sum_equivariant<R: Real, const N: usize>(mut t: [R; N]) -> R {
    sort_by_magnitude(&mut t);
    let mut acc = R::ZERO;
    let mut i = 0;
    while i < N {
        acc += t[i];
        i += 1;
    }
    acc
}

#[cfg(test)]
mod tests;
