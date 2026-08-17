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
    let a = air_of(8);
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

    a.dig(&[[1, 1, 1], [2, 1, 1]], || true);
    a.dig(&[[5, 1, 1], [6, 1, 1]], || true);
    assert_eq!(a.components(), 2, "two separate caves");
    assert!(!a.connected([1, 1, 1], [6, 1, 1]));

    // Still two: this widens one cave without reaching the other.
    let r = a.dig(&[[3, 1, 1]], || true);
    assert_eq!(r.dirty, 1);
    assert_eq!(r.merges, 1, "joined to the left cave only");
    assert_eq!(a.components(), 2);

    // The breakthrough.
    let r = a.dig(&[[4, 1, 1]], || true);
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

    let first = a.dig(&cells, || true);
    assert_eq!(first.dirty, 3);
    assert_eq!(first.already_air, 0);

    let second = a.dig(&cells, || true);
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
fn label_writes_never_exceed_the_lattice_degree() {
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
        let r = a.dig(&cells, || true);
        assert!(r.dirty > 0);
        // Was `unions <= 6 * dirty` when this was a union-find. Under flat
        // labels (✗26) there are no union calls, and the instrument is label
        // writes -- so the invariant is restated rather than dropped: every
        // newly-air sample is written at least once, and the repair must not
        // touch more than the lattice degree per sample.
        assert!(
            r.relabels >= r.dirty,
            "n={n}: {} relabels for {} dirty leaves a sample unwritten",
            r.relabels,
            r.dirty
        );
        assert!(
            r.relabels <= 6 * r.dirty,
            "n={n}: {} relabels for {} dirty is more than six per sample",
            r.relabels,
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
    let f = forward.dig(&line, || true);

    let mut reversed: Vec<[u32; 3]> = line.clone();
    reversed.reverse();
    let mut backward = air_of(8);
    let b = backward.dig(&reversed, || true);

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
        incremental.dig(&[*c], || true);
    }

    let mut values = alloc::vec![-1.0_f64; shape.element_count()];
    for c in &cells {
        let i = (c[2] as usize * N as usize + c[1] as usize) * N as usize + c[0] as usize;
        values[i] = 1.0;
    }
    let (rebuilt, _) = Air::build(&values, &shape).expect("build");

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
    let r = a.dig(&[[1, 1, 1], [99, 1, 1], [1, 99, 1]], || true);
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

/// **Filling the midpoint of a tunnel severs the two caverns it joined.**
///
/// The smallest form of P-26's adversarial fixture, and the sealed-volume
/// mechanic in miniature: this is the edit the connectivity layer exists for.
#[test]
fn filling_a_tunnel_midpoint_severs_the_two_caverns() {
    let mut a = air_of(9);
    let line: Vec<[u32; 3]> = (1..8).map(|x| [x, 4, 4]).collect();
    a.dig(&line, || true);
    assert_eq!(a.components(), 1);
    assert!(a.connected([1, 4, 4], [7, 4, 4]));

    let f = a.fill(&[[4, 4, 4]], || true);
    assert_eq!(f.dirty, 1);
    assert_eq!(f.splits, 1, "one component was severed");
    assert_eq!(f.shed, 1, "into two pieces, so one is new");
    assert_eq!(f.vanished, 0);
    assert_eq!(a.components(), 2);
    assert!(!a.connected([1, 4, 4], [7, 4, 4]), "the passage is sealed");
    assert!(a.connected([1, 4, 4], [3, 4, 4]), "the left side is intact");
    assert!(a.connected([5, 4, 4], [7, 4, 4]), "so is the right");
}

/// **A fill that does not disconnect anything reports no split.**
///
/// Most filling is this — M-319 measured five fills in six changing nothing.
#[test]
fn filling_the_end_of_a_tunnel_shortens_it_without_splitting() {
    let mut a = air_of(9);
    let line: Vec<[u32; 3]> = (1..8).map(|x| [x, 4, 4]).collect();
    a.dig(&line, || true);

    let f = a.fill(&[[7, 4, 4]], || true);
    assert_eq!(f.dirty, 1);
    assert_eq!(f.splits, 0);
    assert_eq!(f.shed, 0);
    assert_eq!(a.components(), 1);
    assert!(a.connected([1, 4, 4], [6, 4, 4]));
}

/// **A component consumed outright vanishes, and needs no search.**
///
/// The other half of M-319's split-or-vanish pair: 27 severed against 5 consumed
/// at 65³. A vanished component has no surviving piece to relabel.
#[test]
fn consuming_a_pocket_vanishes_it() {
    let mut a = air_of(9);
    a.dig(&[[1, 1, 1], [2, 1, 1]], || true);
    a.dig(&[[6, 6, 6]], || true);
    assert_eq!(a.components(), 2);

    let f = a.fill(&[[6, 6, 6]], || true);
    assert_eq!(f.vanished, 1);
    assert_eq!(f.splits, 0);
    assert_eq!(a.components(), 1);
    assert!(!a.connected([6, 6, 6], [1, 1, 1]));
}

/// **Filling agrees with a rebuild — P-26's serious falsifier.**
///
/// Cost measurement cannot see a structure that is fast and wrong, so this is an
/// assertion rather than a benchmark. Every air sample's component membership
/// must match a from-scratch flood fill, not merely the component *count*: two
/// wrong labellings can agree on the total.
#[test]
fn filling_agrees_with_a_rebuild() {
    let n = 12;
    let (mut values, shape) = solid(n);
    // A cross of tunnels, so fills have something to sever.
    let mut air_cells: Vec<[u32; 3]> = Vec::new();
    for t in 1..n - 1 {
        air_cells.push([t, 5, 5]);
        air_cells.push([5, t, 5]);
        air_cells.push([5, 5, t]);
    }
    let (mut a, _) = Air::build(&values, &shape).expect("build");
    a.dig(&air_cells, || true);
    for c in &air_cells {
        let i = (c[2] as usize * n as usize + c[1] as usize) * n as usize + c[0] as usize;
        if let Some(v) = values.get_mut(i) {
            *v = 1.0;
        }
    }

    // Fill one arm's midpoint at a time; each is a potential severance.
    for cut in [[3, 5, 5], [5, 3, 5], [5, 5, 3], [8, 5, 5], [5, 8, 5]] {
        a.fill(&[cut], || true);
        let i = (cut[2] as usize * n as usize + cut[1] as usize) * n as usize + cut[0] as usize;
        if let Some(v) = values.get_mut(i) {
            *v = -1.0;
        }
        let (rebuilt, _) = Air::build(&values, &shape).expect("rebuild");

        assert_eq!(
            a.components(),
            rebuilt.components(),
            "component count diverged after filling {cut:?}"
        );
        assert_eq!(a.air_samples(), rebuilt.air_samples());

        // Membership, not just the count. Every pair of air samples must agree.
        for p in &air_cells {
            for q in &air_cells {
                assert_eq!(
                    a.connected(*p, *q),
                    rebuilt.connected(*p, *q),
                    "after filling {cut:?}, connected({p:?}, {q:?}) diverged"
                );
            }
        }
    }
}

/// **A budget of zero defers the repair, and the deferred answer is the safe
/// one.**
///
/// Before a `fill` repair completes the severed pieces still share a label, so
/// `connected` says yes and a caller reads *"not sealed yet"*. Water leaking for
/// three frames is recoverable; water not leaking out of a room the engine
/// wrongly believes is sealed is a broken game rule.
#[test]
fn an_exhausted_budget_defers_and_reads_conservatively() {
    let mut a = air_of(9);
    let line: Vec<[u32; 3]> = (1..8).map(|x| [x, 4, 4]).collect();
    a.dig(&line, || true);

    let f = a.fill(&[[4, 4, 4]], || false);
    assert_eq!(f.dirty, 1, "the removal itself always happens");
    assert_eq!(f.splits, 0, "but the search did not run");
    assert_eq!(f.pending, 1);
    assert_eq!(a.pending(), 1);
    assert!(
        a.connected([1, 4, 4], [7, 4, 4]),
        "stale reads as still-connected, which is the safe direction"
    );

    // Draining finishes the job and the answer corrects itself.
    let drained = a.repair(&mut || true);
    assert_eq!(drained.splits, 1);
    assert_eq!(a.pending(), 0);
    assert!(!a.connected([1, 4, 4], [7, 4, 4]));
    assert_eq!(a.components(), 2);
}

/// **`label_of` exposes what a chunk-stitching layer needs.**
///
/// Two samples share a label exactly when they are connected, which is what lets
/// a higher layer join `Air`s across a shared face without this type knowing
/// anything about chunks.
#[test]
fn label_of_names_the_component_for_stitching() {
    let mut a = air_of(9);
    a.dig(&[[1, 1, 1], [2, 1, 1]], || true);
    a.dig(&[[6, 6, 6]], || true);

    let left = a.label_of([1, 1, 1]).expect("air");
    assert_eq!(a.label_of([2, 1, 1]), Some(left), "same cave, same label");
    assert_ne!(a.label_of([6, 6, 6]), Some(left), "other cave, other label");
    assert_eq!(a.label_of([0, 0, 0]), None, "solid has no label");
    assert_eq!(a.label_of([99, 0, 0]), None, "off the lattice has no label");
}

/// The R-036 area invariant, checked from outside: a full recount of
/// air–solid faces per label must equal the maintained counts, and a
/// label-free global face total must equal their sum.
fn assert_area_invariant(a: &Air) {
    assert_eq!(a.pending(), 0, "the invariant is a drained-state contract");
    let mut recount = alloc::vec![0u32; a.label_count()];
    let mut global = 0u64;
    let mut nb = [0usize; 6];
    let count = a.air.len();
    for i in 0..count {
        if a.air.get(i) != Some(&true) {
            continue;
        }
        let used = a.neighbours(i, &mut nb);
        let air_faces = nb
            .iter()
            .take(used)
            .filter(|&&j| a.air.get(j) == Some(&true))
            .count() as u32;
        let solid = 6 - air_faces;
        global += u64::from(solid);
        let Some(&l) = a.label.get(i) else { continue };
        if l == NONE {
            continue;
        }
        if let Some(slot) = recount.get_mut(l as usize) {
            *slot += solid;
        }
    }
    let mut maintained_sum = 0u64;
    for l in 0..a.label_count() as u32 {
        assert_eq!(
            a.component_area(l),
            recount.get(l as usize).copied().unwrap_or(0),
            "maintained area for label {l} diverges from the recount"
        );
        maintained_sum += u64::from(a.component_area(l));
    }
    assert_eq!(
        maintained_sum, global,
        "label-free global face total diverges from the per-label sum"
    );
}

/// Every synchronous op keeps the area invariant: build, dig with merges,
/// fill with splits and vanishes, on the P-26-shaped bisect fixture.
#[test]
fn the_area_accumulator_survives_dig_merge_fill_split_and_vanish() {
    const N: u32 = 10;
    let mut a = air_of(N);
    assert_area_invariant(&a);

    // Two caverns and a tunnel — the bisect shape.
    let mut cells = alloc::vec::Vec::new();
    for x in 1..4u32 {
        for y in 4..7u32 {
            for z in 4..7u32 {
                cells.push([x, y, z]);
            }
        }
    }
    for x in 6..9u32 {
        for y in 4..7u32 {
            for z in 4..7u32 {
                cells.push([x, y, z]);
            }
        }
    }
    a.dig(&cells, || true);
    assert_area_invariant(&a);
    assert_eq!(a.components(), 2);

    // The tunnel merges them — and the merged component's area is the
    // recount, not the sum of stale halves.
    a.dig(&[[4, 5, 5], [5, 5, 5]], || true);
    assert_area_invariant(&a);
    assert_eq!(a.components(), 1);

    // Filling the tunnel midpoint severs them again (the split hand-off).
    a.fill(&[[4, 5, 5]], || true);
    assert_area_invariant(&a);
    assert_eq!(a.components(), 2);

    // Consume one cavern entirely (retire must leave its area at zero).
    let mut right = alloc::vec::Vec::new();
    right.push([5, 5, 5]);
    for x in 6..9u32 {
        for y in 4..7u32 {
            for z in 4..7u32 {
                right.push([x, y, z]);
            }
        }
    }
    a.fill(&right, || true);
    assert_area_invariant(&a);
    assert_eq!(a.components(), 1);
}

/// A seeded random op sequence holds the invariant at every drained step —
/// merges, splits, re-digs and vanishes included.
#[test]
fn the_area_invariant_holds_across_a_seeded_op_sequence() {
    const N: u32 = 12;
    let mut a = air_of(N);
    let mut state = 0x00C0_FFEE_u64;
    let mut next = move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 33) as u32
    };
    for step in 0..60 {
        let cx = 1 + next() % (N - 2);
        let cy = 1 + next() % (N - 2);
        let cz = 1 + next() % (N - 2);
        let r = 1 + next() % 3;
        let mut ball = alloc::vec::Vec::new();
        for x in cx.saturating_sub(r)..(cx + r + 1).min(N) {
            for y in cy.saturating_sub(r)..(cy + r + 1).min(N) {
                for z in cz.saturating_sub(r)..(cz + r + 1).min(N) {
                    let d = (i64::from(x) - i64::from(cx)).pow(2)
                        + (i64::from(y) - i64::from(cy)).pow(2)
                        + (i64::from(z) - i64::from(cz)).pow(2);
                    if d <= i64::from(r * r) {
                        ball.push([x, y, z]);
                    }
                }
            }
        }
        if step % 3 == 2 {
            a.fill(&ball, || true);
        } else {
            a.dig(&ball, || true);
        }
        assert_area_invariant(&a);
    }
}

/// The oracle itself can go red: a deliberate one-count corruption must be
/// caught. (`should_panic` is the inversion — an oracle that cannot fail is
/// not checking anything.)
#[test]
#[should_panic(expected = "diverges from the recount")]
fn the_area_oracle_goes_red_on_a_corrupted_count() {
    const N: u32 = 8;
    let mut a = air_of(N);
    a.dig(&[[3, 3, 3], [4, 3, 3]], || true);
    let l = a.label_of([3, 3, 3]).expect("dug air has a label");
    if let Some(slot) = a.area.get_mut(l as usize) {
        *slot += 1;
    }
    assert_area_invariant(&a);
}

/// The split hands both sides a non-zero area, and the Sabine read is two
/// lookups — the structural half of R-036's first clause.
#[test]
fn a_severed_tunnel_leaves_both_components_with_boundary_area() {
    const N: u32 = 10;
    let mut a = air_of(N);
    let mut cells = alloc::vec::Vec::new();
    for x in 1..9u32 {
        cells.push([x, 5, 5]);
    }
    a.dig(&cells, || true);
    a.fill(&[[4, 5, 5]], || true);
    assert_eq!(a.components(), 2);
    let left = a.label_of([1, 5, 5]).expect("left survives");
    let right = a.label_of([8, 5, 5]).expect("right survives");
    assert!(a.component_area(left) > 0, "left side has boundary faces");
    assert!(a.component_area(right) > 0, "right side has boundary faces");
    // A 1×1×3 tube sealed at both ends: 3 samples × 6 faces − 2·2 shared.
    assert_eq!(a.component_size(left), 3);
    assert_eq!(a.component_area(left), 14);
}
