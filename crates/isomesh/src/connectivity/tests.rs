use super::*;
use crate::RuntimeShape3;

/// All-solid values on an `n³` lattice.
fn solid(n: u32) -> (Vec<f64>, RuntimeShape3) {
    let shape = RuntimeShape3::new([n; 3]).expect("valid shape");
    (alloc::vec![-1.0_f64; shape.element_count()], shape)
}

/// A solid block, so every sample starts non-air.
fn air_of(n: u32) -> Air {
    let (values, shape) = solid(n);
    Air::build(&values, &shape).expect("build").0
}

#[test]
fn a_solid_block_has_no_air_and_no_components() {
    let mut a = air_of(8);
    assert_eq!(a.air_samples(), 0);
    assert_eq!(a.components(), 0);
    assert!(!a.connected([1, 1, 1], [2, 1, 1]), "solid is not connected");
}

/// **Two digs that do not touch stay two components; the dig that bridges them
/// merges both.**
///
/// This is the mechanic R-022 exists for — *did I break through* — reduced to
/// the smallest fixture that can show it.
#[test]
fn digging_through_merges_two_caves_and_the_merge_is_the_event() {
    let mut a = air_of(8);

    a.dig(&[[1, 1, 1], [2, 1, 1]]);
    a.dig(&[[5, 1, 1], [6, 1, 1]]);
    assert_eq!(a.components(), 2, "two separate caves");
    assert!(!a.connected([1, 1, 1], [6, 1, 1]));

    // Still two: this widens one cave without reaching the other.
    let r = a.dig(&[[3, 1, 1]]);
    assert_eq!(r.dirty, 1);
    assert_eq!(r.merges, 1, "joined to the left cave only");
    assert_eq!(a.components(), 2);

    // The breakthrough.
    let r = a.dig(&[[4, 1, 1]]);
    assert_eq!(r.dirty, 1);
    assert_eq!(
        r.merges, 2,
        "one join each side, and the second is the event"
    );
    assert_eq!(a.components(), 1);
    assert!(a.connected([1, 1, 1], [6, 1, 1]));
}

/// **Digging where you have already dug costs nothing.**
///
/// A brush applied twice has an empty dirty set the second time. Without
/// `already_air` a harness could not tell that from a brush that missed.
#[test]
fn the_same_brush_twice_has_an_empty_dirty_set() {
    let mut a = air_of(8);
    let cells = [[1, 1, 1], [2, 1, 1], [3, 1, 1]];

    let first = a.dig(&cells);
    assert_eq!(first.dirty, 3);
    assert_eq!(first.already_air, 0);

    let second = a.dig(&cells);
    assert_eq!(second.dirty, 0);
    assert_eq!(second.already_air, 3);
    assert_eq!(second.merges, 0, "nothing left to merge");
}

/// **P-23's second falsifier, as a test: at most six unions per newly-air
/// sample.**
///
/// Six is the lattice degree. A harness walking a 26-neighbourhood, or visiting
/// each edge from both ends, would exceed it — and would still produce the flat
/// curve against `n³` that H predicts, which is why this is asserted separately
/// from the flatness. Flat is the answer we want, so flat is where a bug hides.
#[test]
fn union_calls_never_exceed_the_lattice_degree() {
    for n in [8u32, 12, 16] {
        let mut a = air_of(n);
        // A solid slab dug in one batch: interior samples have all six
        // neighbours air, which is the worst case for this bound.
        let mut cells = alloc::vec::Vec::new();
        for z in 2..n - 2 {
            for y in 2..n - 2 {
                for x in 2..n - 2 {
                    cells.push([x, y, z]);
                }
            }
        }
        let r = a.dig(&cells);
        assert!(r.dirty > 0);
        assert!(
            r.unions <= 6 * r.dirty,
            "n={n}: {} unions for {} dirty is more than six per sample",
            r.unions,
            r.dirty
        );
    }
}

/// **A batch joins samples that are newly air together**, which is what the
/// two-pass order in `dig` buys.
///
/// Digging a line in one call must give one component. If the phase field were
/// updated and linked in the same pass, the answer would depend on the order the
/// slice happened to be in — so this is run forwards and backwards and required
/// to agree.
#[test]
fn a_batch_is_order_independent() {
    let line: Vec<[u32; 3]> = (1..7).map(|x| [x, 3, 3]).collect();
    let mut forward = air_of(8);
    let f = forward.dig(&line);

    let mut reversed: Vec<[u32; 3]> = line.clone();
    reversed.reverse();
    let mut backward = air_of(8);
    let b = backward.dig(&reversed);

    assert_eq!(forward.components(), 1);
    assert_eq!(backward.components(), 1);
    assert_eq!(f, b, "the same batch in a different order cost the same");
}

/// Building from values and digging the same samples reach the same components.
///
/// The incremental path is only worth having if it agrees with the batch one,
/// and this is the check that says so rather than assuming it.
#[test]
fn incremental_digging_agrees_with_a_rebuild() {
    const N: u32 = 10;
    let shape = RuntimeShape3::new([N; 3]).expect("valid shape");

    // Two disjoint boxes, then a bridge.
    let mut cells = alloc::vec::Vec::new();
    for x in 1..4u32 {
        cells.push([x, 5, 5]);
    }
    for x in 6..9u32 {
        cells.push([x, 5, 5]);
    }
    cells.push([4, 5, 5]);
    cells.push([5, 5, 5]);

    let mut incremental = air_of(N);
    for c in &cells {
        incremental.dig(&[*c]);
    }

    let mut values = alloc::vec![-1.0_f64; shape.element_count()];
    for c in &cells {
        let i = (c[2] as usize * N as usize + c[1] as usize) * N as usize + c[0] as usize;
        values[i] = 1.0;
    }
    let (mut rebuilt, _) = Air::build(&values, &shape).expect("build");

    assert_eq!(incremental.air_samples(), rebuilt.air_samples());
    assert_eq!(incremental.components(), rebuilt.components());
    assert_eq!(incremental.components(), 1);
}

/// A value of exactly zero is air, matching `cube::is_inside` and therefore
/// matching what every extractor here decides about the same sample.
#[test]
fn exactly_zero_is_air_like_everywhere_else() {
    let shape = RuntimeShape3::new([4; 3]).expect("valid shape");
    let mut values = alloc::vec![-1.0_f64; shape.element_count()];
    values[0] = 0.0;
    let (a, _) = Air::build(&values, &shape).expect("build");
    assert_eq!(a.air_samples(), 1);
}

/// A brush running off the grid is ordinary, not an error.
#[test]
fn a_brush_over_the_edge_ignores_what_is_not_there() {
    let mut a = air_of(4);
    let r = a.dig(&[[1, 1, 1], [99, 1, 1], [1, 99, 1]]);
    assert_eq!(r.dirty, 1);
    assert_eq!(a.air_samples(), 1);
}

/// A mismatched value slice is refused rather than truncated.
#[test]
fn the_wrong_number_of_values_is_an_error() {
    let shape = RuntimeShape3::new([4; 3]).expect("valid shape");
    let values = alloc::vec![-1.0_f64; 5];
    assert!(Air::build(&values, &shape).is_err());
}
