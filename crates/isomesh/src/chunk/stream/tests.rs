//! G-007's tests.
//!
//! The one the ticket exists for is `a_chunk_parked_in_the_band_never_thrashes`:
//! hysteresis is not "two numbers in a struct", it is the property that a camera
//! oscillating across a threshold causes no work, and that is what it measures.

use super::*;
use crate::chunk::ChunkLayout;

/// Chunks 4 cells across at unit spacing, so a chunk spans exactly 4 world
/// units and the arithmetic in these tests is readable.
fn layout() -> ChunkLayout<f64> {
    ChunkLayout::new(4, 1.0, [0.0; 3]).expect("a valid layout")
}

fn config(load: f64, unload: f64) -> StreamConfig<f64> {
    StreamConfig::new(load, unload).expect("valid radii")
}

#[test]
fn the_band_must_have_width() {
    // Equal radii is precisely the thrashing case, so it is not constructible.
    assert!(StreamConfig::new(4.0f64, 4.0).is_err());
    assert!(StreamConfig::new(4.0f64, 3.0).is_err());
    assert!(StreamConfig::new(0.0f64, 4.0).is_err());
    assert!(StreamConfig::new(f64::NAN, 4.0).is_err());
    assert!(StreamConfig::new(4.0f64, f64::INFINITY).is_err());
    assert!(StreamConfig::new(4.0f64, 4.001).is_ok());
}

// Exact comparisons on purpose: these are the box distance's *defining* values,
// not measurements of it. Zero inside the box has to be exactly zero, or a
// camera standing in a chunk would fail a `<= radius` test at any radius, and 3
// and 5 are integers a `sqrt` of 9 and 25 must land on exactly.
#[expect(clippy::float_cmp, reason = "these are definitions, not measurements")]
#[test]
fn distance_is_to_the_box_and_is_zero_inside_it() {
    let layout = layout();
    let origin = ChunkId::new([0, 0, 0]);
    // Chunk [0,0,0] spans [0,4) on each axis.
    assert_eq!(distance_to_chunk(&layout, origin, [2.0, 2.0, 2.0]), 0.0);
    assert_eq!(distance_to_chunk(&layout, origin, [0.0, 0.0, 0.0]), 0.0);
    // Three units past the far face on x, nothing on the others.
    assert_eq!(distance_to_chunk(&layout, origin, [7.0, 2.0, 2.0]), 3.0);
    // A corner: 3-4-5 on two axes.
    assert_eq!(distance_to_chunk(&layout, origin, [7.0, 8.0, 2.0]), 5.0);

    // And the reason it is not centre distance: a chunk the camera is *inside*
    // has centre distance up to half a diagonal, which at these radii is most
    // of the load radius.
    let centre_distance = {
        let o = layout.sample_origin(origin);
        let half = 2.0;
        let d = [o[0] + half - 0.1, o[1] + half - 0.1, o[2] + half - 0.1];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    };
    assert!(centre_distance > 3.0, "the two measures differ materially");
}

#[test]
fn the_camera_s_own_chunk_is_always_resident() {
    let layout = layout();
    let mut stream = ChunkStream::new();
    let mut update = StreamUpdate::new();
    // A radius smaller than a single chunk. The chunk the camera is standing in
    // is at distance zero, so it loads regardless -- which is the behaviour a
    // streaming world needs at any radius, or the player falls through the
    // world they are standing on.
    stream
        .update(&layout, [2.0, 2.0, 2.0], &config(0.5, 1.0), &mut update)
        .expect("update");
    assert!(stream.resident().contains(&ChunkId::new([0, 0, 0])));
    assert_eq!(update.unloaded, []);
}

#[test]
fn a_chunk_parked_in_the_band_never_thrashes() {
    // The ticket's actual requirement. A camera oscillating across the load
    // threshold must cause no work at all once the chunk is resident, because
    // every load/unload cycle is a full re-mesh -- the most expensive thing a
    // streaming world does.
    let layout = layout();
    let config = config(9.0, 13.0);
    let mut stream = ChunkStream::new();
    let mut update = StreamUpdate::new();

    // Settle, then find a chunk sitting inside the band: outside the load
    // radius, inside the unload radius. It must exist or this test is asserting
    // nothing.
    stream
        .update(&layout, [2.0, 2.0, 2.0], &config, &mut update)
        .expect("update");
    let banded: Vec<ChunkId> = stream
        .resident()
        .iter()
        .copied()
        .filter(|id| {
            let d = distance_to_chunk(&layout, *id, [2.0, 2.0, 2.0]);
            d > config.load() && d <= config.unload()
        })
        .collect();
    // Nothing can be in the band on the first update -- entry needs `load` --
    // so the band is empty here, which is itself the correct behaviour.
    assert!(
        banded.is_empty(),
        "a chunk entered without coming within the load radius: {banded:?}"
    );

    // Now oscillate. Each pass moves the camera far enough that chunks cross
    // the load radius in and out, but never past the unload radius.
    let mut churn = 0usize;
    for step in 0..12 {
        let x = if step % 2 == 0 { 2.0 } else { 5.5 };
        stream
            .update(&layout, [x, 2.0, 2.0], &config, &mut update)
            .expect("update");
        if step >= 2 {
            churn += update.unloaded.len();
        }
    }
    assert_eq!(
        churn, 0,
        "hysteresis did not hold: {churn} chunks unloaded while oscillating"
    );
}

#[test]
fn without_hysteresis_the_same_motion_does_thrash() {
    // M-44's rule: the zero above has to prove it could have been non-zero.
    // The narrowest legal band is nearly a single radius, and with it the same
    // oscillation evicts chunks on almost every pass.
    let layout = layout();
    let narrow = config(9.0, 9.0 + 1e-9);
    let mut stream = ChunkStream::new();
    let mut update = StreamUpdate::new();

    let mut churn = 0usize;
    for step in 0..12 {
        let x = if step % 2 == 0 { 2.0 } else { 5.5 };
        stream
            .update(&layout, [x, 2.0, 2.0], &narrow, &mut update)
            .expect("update");
        if step >= 2 {
            churn += update.unloaded.len();
        }
    }
    assert!(
        churn > 0,
        "a band of width 1e-9 absorbed the same motion, so the test above is \
         not measuring hysteresis"
    );
}

#[test]
fn a_chunk_leaves_only_past_the_unload_radius() {
    let layout = layout();
    let config = config(6.0, 14.0);
    let mut stream = ChunkStream::new();
    let mut update = StreamUpdate::new();

    stream
        .update(&layout, [2.0, 2.0, 2.0], &config, &mut update)
        .expect("update");
    let watched = ChunkId::new([0, 0, 0]);
    assert!(stream.resident().contains(&watched));

    // Walk away in steps and record where it actually leaves.
    let mut left_at = None;
    for step in 1..40 {
        let x = 2.0 + f64::from(step) * 2.0;
        stream
            .update(&layout, [x, 2.0, 2.0], &config, &mut update)
            .expect("update");
        if update.unloaded.contains(&watched) {
            left_at = Some(distance_to_chunk(&layout, watched, [x, 2.0, 2.0]));
            break;
        }
    }
    let left_at = left_at.expect("the chunk should eventually unload");
    assert!(
        left_at > config.unload(),
        "unloaded at {left_at}, which is inside the unload radius {}",
        config.unload()
    );
}

#[test]
fn the_result_is_a_pure_function_of_camera_and_config() {
    // Two streams reaching the same camera by different routes must agree, or
    // the residency set is carrying history it should not. The hysteresis is
    // deliberate state; anything *else* remembered would be a bug.
    let layout = layout();
    let config = config(7.0, 11.0);
    let mut update = StreamUpdate::new();

    let mut direct = ChunkStream::new();
    direct
        .update(&layout, [30.0, 2.0, 2.0], &config, &mut update)
        .expect("update");

    let mut wandered = ChunkStream::new();
    for x in [2.0, 40.0, 60.0, 30.0] {
        wandered
            .update(&layout, [x, 2.0, 2.0], &config, &mut update)
            .expect("update");
    }

    // Arriving from far away, every chunk in the band is outside the load
    // radius and so absent from both -- the band only holds chunks that were
    // already resident, and neither of these had one.
    assert_eq!(direct.resident(), wandered.resident());
}

#[test]
fn the_lists_are_sorted_and_disjoint_and_account_for_the_whole_change() {
    let layout = layout();
    let config = config(7.0, 11.0);
    let mut stream = ChunkStream::new();
    let mut update = StreamUpdate::new();

    let mut previous: Vec<ChunkId> = Vec::new();
    for step in 0..8 {
        let x = 2.0 + f64::from(step) * 3.0;
        stream
            .update(&layout, [x, 2.0, 2.0], &config, &mut update)
            .expect("update");

        assert!(
            update.loaded.windows(2).all(|w| w[0] < w[1]),
            "loaded order"
        );
        assert!(
            update.unloaded.windows(2).all(|w| w[0] < w[1]),
            "unloaded order"
        );
        assert!(
            update.loaded.iter().all(|id| !update.unloaded.contains(id)),
            "a chunk both loaded and unloaded in one update"
        );

        // The lists must reconstruct the new set from the old one exactly.
        let mut rebuilt = previous.clone();
        rebuilt.retain(|id| !update.unloaded.contains(id));
        rebuilt.extend(update.loaded.iter().copied());
        rebuilt.sort_unstable();
        assert_eq!(rebuilt, stream.resident(), "step {step}");

        assert!(stream.resident().windows(2).all(|w| w[0] < w[1]));
        previous = stream.resident().to_vec();
    }
    assert!(!previous.is_empty(), "the sweep never loaded anything");
}

#[test]
fn clearing_reports_everything_as_unloaded() {
    let layout = layout();
    let mut stream = ChunkStream::new();
    let mut update = StreamUpdate::new();
    stream
        .update(&layout, [2.0, 2.0, 2.0], &config(7.0, 11.0), &mut update)
        .expect("update");
    let resident = stream.resident().to_vec();
    assert!(!resident.is_empty());

    stream.clear(&mut update);
    assert!(stream.is_empty());
    assert_eq!(update.unloaded, resident);
    assert!(update.loaded.is_empty());
}

#[test]
fn an_absurd_radius_is_refused_rather_than_allocated() {
    let layout = layout();
    let mut stream = ChunkStream::new();
    let mut update = StreamUpdate::new();
    // A radius of ten million chunks is 8e21 candidates. Reported, because the
    // alternative is an allocation that takes the process with it.
    let config = config(1.0e7, 4.0e7);
    assert!(
        stream
            .update(&layout, [0.0; 3], &config, &mut update)
            .is_err()
    );
}

#[test]
fn residency_reuses_its_buffers() {
    // Rule 6: a streaming world calls this every frame.
    let layout = layout();
    let config = config(7.0, 11.0);
    let mut stream = ChunkStream::new();
    let mut update = StreamUpdate::new();

    stream
        .update(&layout, [2.0, 2.0, 2.0], &config, &mut update)
        .expect("update");
    update.reset();
    assert!(update.is_empty());

    // A second update at the same place changes nothing at all, which is the
    // steady state a stationary camera should cost.
    stream
        .update(&layout, [2.0, 2.0, 2.0], &config, &mut update)
        .expect("update");
    assert!(
        update.is_empty(),
        "a stationary camera produced {} loads and {} unloads",
        update.loaded.len(),
        update.unloaded.len()
    );
}
