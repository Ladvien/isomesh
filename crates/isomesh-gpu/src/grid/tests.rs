//! `GridParams` needs no GPU to be wrong, so none of this asks for one.
//!
//! Exact float comparison throughout, deliberately: the packing is a copy, not
//! arithmetic, so a tolerance here would absorb exactly the layout bugs these
//! tests exist to catch.
#![allow(clippy::float_cmp)]

use super::GridParams;
use crate::Error;

#[test]
fn a_grid_needs_two_samples_on_every_axis() {
    for samples in [[1, 33, 33], [33, 1, 33], [33, 33, 1], [0, 0, 0]] {
        assert_eq!(
            GridParams::new(samples, [0.0; 3], 0.125),
            Err(Error::DegenerateGrid { samples }),
            "accepted {samples:?}"
        );
    }
    assert!(GridParams::new([2; 3], [0.0; 3], 0.125).is_ok());
}

#[test]
fn geometry_must_be_finite_and_positive() {
    for bad in [0.0, -0.125, f32::NAN, f32::INFINITY] {
        assert_eq!(
            GridParams::new([33; 3], [0.0; 3], bad),
            Err(Error::InvalidCellSize),
            "accepted cell size {bad}"
        );
    }
    assert_eq!(
        GridParams::new([33; 3], [0.0, f32::NAN, 0.0], 0.125),
        Err(Error::InvalidOrigin)
    );
}

/// The overflow is rejected at construction so every buffer size downstream can
/// be plain arithmetic rather than a checked chain.
#[test]
fn a_grid_too_large_to_address_is_refused() {
    let samples = [u32::MAX, u32::MAX, u32::MAX];
    assert_eq!(
        GridParams::new(samples, [0.0; 3], 0.125),
        Err(Error::GridTooLarge { samples })
    );
}

#[test]
fn counts_are_samples_and_cells_respectively() {
    let grid = GridParams::new([33, 17, 9], [0.0; 3], 0.125).expect("valid grid");
    assert_eq!(grid.sample_count(), 33 * 17 * 9);
    assert_eq!(grid.cell_count(), 32 * 16 * 8);
    assert_eq!(grid.field_buffer_size(), 33 * 17 * 9 * 4);
}

/// M-70 and M-73 both record cracks caused by accumulating a sample cursor
/// instead of multiplying the index. This pins the multiply.
#[test]
fn a_sample_position_is_the_index_multiplied_not_accumulated() {
    // 4/35 is deliberately not a power of two -- the spacing where the two
    // expressions disagree in the last bit.
    let h = 4.0f32 / 35.0;
    let grid = GridParams::new([64, 2, 2], [-2.0, 0.0, 0.0], h).expect("valid grid");

    let mut accumulated = -2.0f32;
    let mut disagreements = 0;
    for i in 0..64u32 {
        let multiplied = grid.sample_position([i, 0, 0])[0];
        assert_eq!(multiplied, -2.0 + h * i as f32);
        if multiplied.to_bits() != accumulated.to_bits() {
            disagreements += 1;
        }
        accumulated += h;
    }
    // And the accumulation really does drift, or the test above proves nothing.
    assert!(
        disagreements > 0,
        "accumulation agreed everywhere -- this fixture no longer exercises the trap"
    );
}

#[test]
fn the_uniform_layout_is_two_vec4s_little_endian() {
    let grid = GridParams::new([33, 17, 9], [-2.0, 0.5, 1.0], 0.125).expect("valid grid");
    let bytes = grid.to_std140();
    assert_eq!(bytes.len() as u64, GridParams::UNIFORM_SIZE);

    let word = |i: usize| u32::from_le_bytes([bytes[i], bytes[i + 1], bytes[i + 2], bytes[i + 3]]);
    assert_eq!(word(0), 33);
    assert_eq!(word(4), 17);
    assert_eq!(word(8), 9);
    assert_eq!(word(12), 0, "padding word must be zero");
    assert_eq!(f32::from_bits(word(16)), -2.0);
    assert_eq!(f32::from_bits(word(20)), 0.5);
    assert_eq!(f32::from_bits(word(24)), 1.0);
    assert_eq!(f32::from_bits(word(28)), 0.125);
}

#[test]
fn packing_is_deterministic() {
    let grid = GridParams::new([65; 3], [-4.0; 3], 1.0 / 3.0).expect("valid grid");
    assert_eq!(grid.to_std140(), grid.to_std140());
}

/// The expression must be `origin + h * i` as **two** operations, matching
/// `isomesh`'s `corner_position`, and not the fused `mul_add`.
///
/// This is the guard for M-143. `sample_position` decides where the field is
/// sampled before upload, so a fused form makes the GPU read a field evaluated
/// at different points from the CPU's, and every downstream comparison then
/// measures that rather than the algorithm.
///
/// The fixture spacing is `0.1`, deliberately. **At a power of two the two
/// forms agree bit for bit and this test cannot fail** — which is exactly why
/// the bug survived a full GPU-004 test suite run at `h = 0.125`.
#[test]
fn a_sample_position_is_not_a_fused_multiply_add() {
    let h = 0.1f32;
    let origin = -2.0f32;
    let grid = GridParams::new([64, 2, 2], [origin, 0.0, 0.0], h).expect("valid grid");

    let mut fused_differs = 0;
    for i in 0..64u32 {
        let separate = origin + h * i as f32;
        let fused = h.mul_add(i as f32, origin);
        assert_eq!(
            grid.sample_position([i, 0, 0])[0].to_bits(),
            separate.to_bits(),
            "sample {i} is not `origin + h * i` evaluated as two operations"
        );
        if separate.to_bits() != fused.to_bits() {
            fused_differs += 1;
        }
    }
    assert!(
        fused_differs > 0,
        "this spacing cannot tell the two forms apart, so the assertion above is vacuous"
    );
}
