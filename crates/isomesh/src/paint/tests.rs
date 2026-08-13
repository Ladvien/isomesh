//! The load-bearing test is E-208's acceptance criterion, stated as an
//! equality rather than a picture: spray a wall, blow a hole through it, and
//! the colour at every surviving point is **bit-identical** to what it was
//! before the carve.
//!
//! "Exactly where you sprayed it" admits a stronger test than a screenshot, and
//! this is it. An L²-nearest transfer between two meshes could only ever pass a
//! tolerance.
//!
//! Exact float comparison is the point throughout — an unpainted surface is
//! *the background colour*, not something near it, and the acceptance test
//! above is an equality of bit patterns.
#![allow(clippy::float_cmp)]

use alloc::vec::Vec;

use super::{Edit, PaintStack, Splat, ramp, shade};
use crate::Sdf;
use crate::brush::{Brush, BrushStack};
use crate::fields::{BoxExact, Sphere};

/// Every log here carves and sprays with spheres, so one alias names both
/// shape parameters. A log holding *only* sprays cannot infer the carve shape
/// on its own, which is why these appear explicitly below.
type Log = Edit<Sphere<f64>, Sphere<f64>, f64>;
type Log32 = Edit<Sphere<f32>, Sphere<f32>, f32>;

/// A wall in the `z ≈ 0` plane, 4×4 across and half a unit thick.
fn wall() -> BoxExact<f64> {
    BoxExact {
        center: [0.0; 3],
        half_extents: [2.0, 2.0, 0.25],
    }
}

const RED: [f64; 4] = [1.0, 0.0, 0.0, 1.0];
const GREY: [f64; 4] = [0.5, 0.5, 0.5, 1.0];

/// Red sprayed onto the front face, centred on the origin.
fn spray() -> Splat<Sphere<f64>, f64> {
    Splat {
        shape: Sphere {
            center: [0.0, 0.0, -0.25],
            radius: 0.75,
        },
        color: RED,
        softness: 0.1,
        depth: 0.05,
    }
}

/// A hole punched through the wall at `x = 1.2`, clear of the paint.
fn hole() -> Brush<Sphere<f64>> {
    Brush::subtract(Sphere {
        center: [1.2, 0.0, 0.0],
        radius: 0.5,
    })
}

/// Points on the painted front face that the hole does not reach.
const SURVIVING: [[f64; 3]; 5] = [
    [0.00, 0.00, -0.25],
    [0.20, 0.10, -0.25],
    [-0.30, 0.25, -0.25],
    [0.15, -0.40, -0.25],
    [-0.45, -0.20, -0.25],
];

/// E-208's acceptance. Not "close" — identical, because the carve changed the
/// surface and the paint was never on the surface.
#[test]
fn paint_on_the_surviving_wall_does_not_move_when_a_hole_is_carved() {
    let before_log: [Log; 1] = [Edit::Spray(spray())];
    let after_log = [Edit::Spray(spray()), Edit::Carve(hole())];

    let before = PaintStack {
        base: wall(),
        edits: &before_log,
        background: GREY,
    };
    let after = PaintStack {
        base: wall(),
        edits: &after_log,
        background: GREY,
    };

    for p in SURVIVING {
        let b = before.color_at(p);
        let a = after.color_at(p);
        assert_eq!(
            b.map(f64::to_bits),
            a.map(f64::to_bits),
            "colour moved at {p:?}: {b:?} -> {a:?}"
        );
        // And it is actually painted, or the equality above is vacuous.
        assert!(b[0] > 0.9, "expected red at {p:?}, got {b:?}");
    }
}

/// The thin-shell half of the answer, and the half that makes it graffiti
/// rather than dyed material.
#[test]
fn the_interior_a_carve_exposes_is_bare() {
    // A hole straight through the middle of the painted patch this time.
    let log = [
        Edit::Spray(spray()),
        Edit::Carve(Brush::subtract(Sphere {
            center: [0.0, 0.0, 0.0],
            radius: 0.35,
        })),
    ];
    let world = PaintStack {
        base: wall(),
        edits: &log,
        background: GREY,
    };

    // Inside the wall's thickness, on the rim the carve just opened: more than
    // `depth` from where the front face used to be, so it was never painted.
    for z in [-0.10, 0.0, 0.10] {
        let p = [0.35, 0.0, z];
        let c = world.color_at(p);
        assert_eq!(c, GREY, "newly exposed interior is painted at {p:?}");
    }
}

/// A splat is confined to its shape as well as to the shell.
#[test]
fn paint_stops_at_the_edge_of_the_spray() {
    let log: [Log; 1] = [Edit::Spray(spray())];
    let world = PaintStack {
        base: wall(),
        edits: &log,
        background: GREY,
    };

    // Radius 0.75 plus 0.1 of softness: 1.0 out is clear of both.
    assert_eq!(world.color_at([1.0, 0.0, -0.25]), GREY);
}

/// The splat's shape is the only thing that decides *reach*, and that has a
/// consequence worth pinning rather than discovering: a shape that passes
/// through a thin wall paints the far face too, because both faces are
/// surface and both are inside the shape.
///
/// This is the design, not a defect — the shape is an arbitrary [`Sdf`], so it
/// already carries whatever directionality a caller wants, and a `direction`
/// field on [`Splat`] would duplicate it. A spray that must not reach through
/// gets a shape that does not reach through.
#[test]
fn a_shape_that_reaches_through_a_wall_paints_both_faces() {
    // Radius 0.75 from z = −0.25 reaches the back face at z = +0.25.
    let reaching: [Log; 1] = [Edit::Spray(spray())];
    let both = PaintStack {
        base: wall(),
        edits: &reaching,
        background: GREY,
    };
    assert!(
        both.color_at([0.0, 0.0, 0.25])[0] > 0.9,
        "a shape spanning the wall paints the far face"
    );

    // Radius 0.3 from the same centre stops inside the wall's thickness.
    let shallow: [Log; 1] = [Edit::Spray(Splat {
        shape: Sphere {
            center: [0.0, 0.0, -0.25],
            radius: 0.3,
        },
        ..spray()
    })];
    let one_face = PaintStack {
        base: wall(),
        edits: &shallow,
        background: GREY,
    };
    assert!(
        one_face.color_at([0.0, 0.0, -0.25])[0] > 0.9,
        "front is painted"
    );
    assert_eq!(
        one_face.color_at([0.0, 0.0, 0.25]),
        GREY,
        "back face is out of the shape's reach"
    );
}

/// Rule 1's real risk is two paths that agree until they don't. This pins that
/// the painted field and the unpainted one are the same arithmetic.
#[test]
fn a_log_without_sprays_samples_bit_identically_to_a_brush_stack() {
    let brushes = [
        hole(),
        Brush::add(Sphere {
            center: [-1.0, 0.5, 0.0],
            radius: 0.4,
        }),
        Brush::smooth_add(
            Sphere {
                center: [0.8, -0.6, 0.1],
                radius: 0.3,
            },
            0.15,
        ),
    ];
    let edits: Vec<Edit<Sphere<f64>, Sphere<f64>, f64>> =
        brushes.iter().copied().map(Edit::Carve).collect();

    let stack = BrushStack {
        base: wall(),
        brushes: &brushes,
    };
    let painted = PaintStack {
        base: wall(),
        edits: &edits,
        background: GREY,
    };

    for i in 0..40 {
        let t = f64::from(i) * 0.1 - 2.0;
        let p = [t, t * 0.5 - 0.3, t * 0.25];
        assert_eq!(
            stack.sample(p).to_bits(),
            painted.sample(p).to_bits(),
            "field diverged at {p:?}"
        );
    }
}

/// Sprays cost the field walk nothing but a `match` arm — asserted through
/// behaviour, since a spray that affected `sample` would move the surface.
#[test]
fn sprays_do_not_move_the_surface() {
    let carve_only = [Edit::<Sphere<f64>, Sphere<f64>, f64>::Carve(hole())];
    let with_spray = [Edit::Spray(spray()), Edit::Carve(hole())];

    let a = PaintStack {
        base: wall(),
        edits: &carve_only,
        background: GREY,
    };
    let b = PaintStack {
        base: wall(),
        edits: &with_spray,
        background: GREY,
    };

    for i in 0..40 {
        let t = f64::from(i) * 0.1 - 2.0;
        let p = [t, t * 0.5 - 0.3, t * 0.25];
        assert_eq!(a.sample(p).to_bits(), b.sample(p).to_bits());
    }
}

/// Later sprays composite over earlier ones.
#[test]
fn sprays_layer_in_log_order() {
    let blue = Splat {
        color: [0.0, 0.0, 1.0, 1.0],
        ..spray()
    };
    let red_then_blue: [Log; 2] = [Edit::Spray(spray()), Edit::Spray(blue)];
    let blue_then_red: [Log; 2] = [Edit::Spray(blue), Edit::Spray(spray())];

    let p = [0.0, 0.0, -0.25];
    let first = PaintStack {
        base: wall(),
        edits: &red_then_blue,
        background: GREY,
    }
    .color_at(p);
    let second = PaintStack {
        base: wall(),
        edits: &blue_then_red,
        background: GREY,
    }
    .color_at(p);

    assert!(
        first[2] > 0.9 && first[0] < 0.1,
        "blue should be on top: {first:?}"
    );
    assert!(
        second[0] > 0.9 && second[2] < 0.1,
        "red should be on top: {second:?}"
    );
}

/// Alpha is coverage, consumed into the weight, not written out.
#[test]
fn a_half_covered_splat_blends_and_leaves_alpha_alone() {
    let half = Splat {
        color: [1.0, 0.0, 0.0, 0.5],
        ..spray()
    };
    let log: [Log; 1] = [Edit::Spray(half)];
    let c = PaintStack {
        base: wall(),
        edits: &log,
        background: GREY,
    }
    .color_at([0.0, 0.0, -0.25]);

    assert!((c[0] - 0.75).abs() < 1e-12, "half of the way to red: {c:?}");
    assert_eq!(c[3], GREY[3], "alpha is the background's");
}

/// The degenerate widths have answers rather than divisions by zero — the same
/// treatment `smooth_min` gives `k <= 0`.
#[test]
fn a_zero_width_ramp_is_a_step() {
    assert_eq!(ramp(-1.0_f64, 0.0), 1.0);
    assert_eq!(ramp(0.0_f64, 0.0), 1.0);
    assert_eq!(ramp(1e-300_f64, 0.0), 0.0);
    // And a negative width does not invert it.
    assert_eq!(ramp(1.0_f64, -1.0), 0.0);
}

#[test]
fn a_positive_width_ramp_is_linear_and_clamped() {
    assert_eq!(ramp(-1.0_f64, 2.0), 1.0);
    assert_eq!(ramp(0.0_f64, 2.0), 1.0);
    assert_eq!(ramp(1.0_f64, 2.0), 0.5);
    assert_eq!(ramp(2.0_f64, 2.0), 0.0);
    assert_eq!(ramp(3.0_f64, 2.0), 0.0);
}

/// Rule 6: the buffer is the caller's and survives being refilled.
#[test]
fn shade_refills_without_releasing_capacity() {
    let log: [Log; 1] = [Edit::Spray(spray())];
    let world = PaintStack {
        base: wall(),
        edits: &log,
        background: GREY,
    };
    let positions: Vec<[f64; 3]> = SURVIVING.into_iter().collect();

    let mut out = Vec::new();
    shade(&positions, &world, &mut out);
    assert_eq!(out.len(), positions.len());
    let capacity = out.capacity();
    let first = out.clone();

    shade(&positions, &world, &mut out);
    assert_eq!(out, first, "shade is not deterministic");
    assert_eq!(out.capacity(), capacity, "shade reallocated");

    shade(&[], &world, &mut out);
    assert!(out.is_empty());
    assert_eq!(out.capacity(), capacity, "shade released capacity");
}

/// `shade` and `color_at` are the same function, so a caller can mix them.
#[test]
fn shade_agrees_with_color_at() {
    let log = [Edit::Spray(spray()), Edit::Carve(hole())];
    let world = PaintStack {
        base: wall(),
        edits: &log,
        background: GREY,
    };
    let positions: Vec<[f64; 3]> = SURVIVING.into_iter().collect();

    let mut out = Vec::new();
    shade(&positions, &world, &mut out);

    for (p, c) in positions.iter().zip(&out) {
        assert_eq!(world.color_at(*p).map(f64::to_bits), c.map(f64::to_bits));
    }
}

/// The crate is generic over `Real`; so is this.
#[test]
fn it_works_in_f32() {
    let log: [Log32; 1] = [Edit::Spray(Splat {
        shape: Sphere {
            center: [0.0_f32, 0.0, -0.25],
            radius: 0.75,
        },
        color: [1.0_f32, 0.0, 0.0, 1.0],
        softness: 0.1,
        depth: 0.05,
    })];
    let world = PaintStack {
        base: BoxExact {
            center: [0.0_f32; 3],
            half_extents: [2.0, 2.0, 0.25],
        },
        edits: &log,
        background: [0.5_f32, 0.5, 0.5, 1.0],
    };
    assert!(world.color_at([0.0, 0.0, -0.25])[0] > 0.9);
    assert_eq!(world.color_at([1.5, 0.0, -0.25]), [0.5, 0.5, 0.5, 1.0]);
}
