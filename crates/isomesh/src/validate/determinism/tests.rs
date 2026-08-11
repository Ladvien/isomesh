//! Fixtures for the determinism harness.
//!
//! The two that matter are the signed-zero and NaN cases: they are where float
//! equality gives the wrong answer in each direction, and a harness that used
//! `==` would pass one and fail the other.

use super::*;
use crate::MeshSink;
use core::cell::Cell;

fn triangle(out: &mut MeshBuffer<f64>) {
    let a = out.vertex([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let b = out.vertex([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let c = out.vertex([0.0, 0.0, 1.0], [0.0, 1.0, 0.0]);
    out.triangle(a, b, c);
}

#[test]
fn a_deterministic_extractor_passes() {
    let report = check_determinism(triangle);
    assert!(report.is_deterministic(), "{report}");
    assert_eq!(report.vertices, 3);
    assert_eq!(report.triangles, 1);
    report.panic_if_divergent();
}

#[test]
fn an_empty_extractor_is_deterministic() {
    let report = check_determinism(|_: &mut MeshBuffer<f32>| {});
    assert!(report.is_deterministic());
    assert_eq!(report.vertices, 0);
    assert_eq!(report.triangles, 0);
}

/// The case float equality gets wrong in one direction: `+0.0 == -0.0` is true,
/// but the bit patterns differ. A sign flip on a zero coordinate is exactly what
/// a reordered summation produces, and a harness comparing with `==` would call
/// this deterministic.
#[test]
fn signed_zero_is_a_divergence() {
    let run = Cell::new(0u32);
    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        let n = run.get();
        run.set(n + 1);
        let z = if n == 0 { 0.0 } else { -0.0 };
        out.vertex([z, 0.0, 0.0], [0.0, 1.0, 0.0]);
    });

    // `0.0 == -0.0` is true, so a harness comparing with `==` would call this
    // deterministic. Bit comparison does not.
    assert!(!report.is_deterministic(), "{report}");
    match report.divergence {
        Some((RunPair::RepeatedCall, Divergence::Position { vertex, axis, .. })) => {
            assert_eq!((vertex, axis), (0, 0));
        }
        other => panic!("expected a position divergence, got {other:?}"),
    }
}

/// The case float equality gets wrong in the other direction: `NaN != NaN`, so
/// `==` would report a divergence where the bits are identical. Meshing a field
/// that produces a NaN is a bug, but it is not a *determinism* bug, and this
/// harness must not claim it is.
#[test]
fn identical_nan_bits_are_not_a_divergence() {
    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        out.vertex([f64::NAN, 0.0, 0.0], [0.0, 1.0, 0.0]);
    });

    // `NaN != NaN` is true, so a harness comparing with `==` would report a
    // divergence here. Bit comparison does not.
    assert!(report.is_deterministic(), "{report}");
}

#[test]
fn a_differing_position_is_located_exactly() {
    let run = Cell::new(0u32);
    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        let n = run.get();
        run.set(n + 1);
        out.vertex([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        out.vertex([1.0, if n == 0 { 2.0 } else { 3.0 }, 0.0], [0.0, 1.0, 0.0]);
    });

    match report.divergence {
        Some((
            RunPair::RepeatedCall,
            Divergence::Position {
                vertex,
                axis,
                first,
                second,
            },
        )) => {
            assert_eq!((vertex, axis), (1, 1));
            assert_eq!((first, second), (2.0, 3.0));
        }
        other => panic!("expected a position divergence, got {other:?}"),
    }
}

/// The bug class the ticket names: iteration order reaching vertex order. The
/// positions are the same set both times, in a different order.
#[test]
fn reordered_output_is_a_divergence() {
    let run = Cell::new(0u32);
    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        let n = run.get();
        run.set(n + 1);
        let order: [f64; 2] = if n == 0 { [1.0, 2.0] } else { [2.0, 1.0] };
        for x in order {
            out.vertex([x, 0.0, 0.0], [0.0, 1.0, 0.0]);
        }
    });
    assert!(!report.is_deterministic(), "{report}");
}

#[test]
fn a_differing_vertex_count_is_reported_as_such() {
    let run = Cell::new(0u32);
    let report = check_determinism(|out: &mut MeshBuffer<f32>| {
        let n = run.get();
        run.set(n + 1);
        for _ in 0..=n {
            out.vertex([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        }
    });

    match report.divergence {
        Some((RunPair::RepeatedCall, Divergence::VertexCount { first, second })) => {
            assert_eq!((first, second), (1, 2));
        }
        other => panic!("expected a vertex-count divergence, got {other:?}"),
    }
}

#[test]
fn a_differing_index_is_located_exactly() {
    let run = Cell::new(0u32);
    let report = check_determinism(|out: &mut MeshBuffer<f32>| {
        let n = run.get();
        run.set(n + 1);
        for _ in 0..3 {
            out.vertex([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        }
        if n == 0 {
            out.triangle(0, 1, 2);
        } else {
            out.triangle(0, 2, 1);
        }
    });

    match report.divergence {
        Some((RunPair::RepeatedCall, Divergence::Index { at, first, second })) => {
            assert_eq!((at, first, second), (1, 1, 2));
        }
        other => panic!("expected an index divergence, got {other:?}"),
    }
}

/// The third run's reason for existing: output that depends on the output
/// buffer's prior state rather than on the input.
///
/// The trigger here is deliberately blunt — reading the buffer's capacity, which
/// is zero when fresh and non-zero after a reset. A real extractor would reach
/// this through a cache or a scratch buffer it forgot to clear. What matters is
/// that such a thing is caught at all, since nothing else in the suite drives an
/// extractor the way the API intends it to be driven.
#[test]
fn output_depending_on_a_reused_buffer_is_a_divergence() {
    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        let reused = out.positions.capacity() > 0;
        out.vertex([if reused { 1.0 } else { 0.0 }, 0.0, 0.0], [0.0, 1.0, 0.0]);
    });

    assert!(!report.is_deterministic(), "{report}");
    let (pair, _) = report.divergence.expect("a divergence");
    assert_eq!(
        pair,
        RunPair::ReusedBuffer,
        "the first two runs both used fresh buffers and agreed"
    );
}

#[test]
#[should_panic(expected = "non-deterministic output between")]
fn panic_if_divergent_names_the_run_pair() {
    let run = Cell::new(0u32);
    check_determinism(|out: &mut MeshBuffer<f32>| {
        let n = run.get();
        run.set(n + 1);
        out.vertex([n as f32, 0.0, 0.0], [0.0, 1.0, 0.0]);
    })
    .panic_if_divergent();
}

#[test]
fn a_normal_divergence_is_distinguished_from_a_position_one() {
    let run = Cell::new(0u32);
    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        let n = run.get();
        run.set(n + 1);
        out.vertex([0.0, 0.0, 0.0], [0.0, if n == 0 { 1.0 } else { -1.0 }, 0.0]);
    });

    match report.divergence {
        Some((RunPair::RepeatedCall, Divergence::Normal { vertex, axis, .. })) => {
            assert_eq!((vertex, axis), (0, 1));
        }
        other => panic!("expected a normal divergence, got {other:?}"),
    }
}
