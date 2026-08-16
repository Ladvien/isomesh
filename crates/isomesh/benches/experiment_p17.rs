//! **P-17 — is Manifold Dual Contouring's residue an interior ambiguity?**
//!
//! Ticket: A-025. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p17
//! ```
//!
//! Writes `docs/experiments/p-17.csv`.
//!
//! # The residue, and why the interior is the suspect
//!
//! ✗19: the paper claims the uniform-grid dual *"is always a manifold"*, and it
//! is not. M-290 measured 143 non-manifold edges with the crate's default table
//! and 114 with the decider-modified one the paper specifies, against Marching
//! Cubes' **0** — so the premise holds and *"the dual preserves the topology"*
//! is what fails.
//!
//! **Every one of them is on `noise_cavity`.** That is the field A-002e added
//! precisely because none of the other seven produces a cell with an **interior**
//! ambiguity (M-208), and a *face* decider cannot see an interior one by
//! construction.
//!
//! # The predicate is the crate's, not a new one
//!
//! `the_defect_count_is_predicted_from_the_grid_alone` already computes the
//! offending pairs from the grid with no mesh involved:
//!
//! ```text
//! non_manifold_edges  ==  shared ambiguous faces whose four cut edges
//!                         lie in one cycle on both sides
//! ```
//!
//! This reuses it exactly — same `Contours::of`, same `joined_mask`, same
//! `edge_on_face` — and adds two things: the interior test on each cell of each
//! pair, and **the same measurement over the ambiguous-face pairs that are not
//! offending**. Without that control the finding would be unfalsifiable:
//! ambiguous faces are where interior ambiguities live, so "offenders have one"
//! could easily be true of every cell in the population.
//!
//! # One convention is asserted rather than assumed
//!
//! [`SweptFaces`] wants `[A, B, C, D]` cyclic with `A`/`C` one diagonal, and
//! building those faces needs the crate's corner-bit order — bit `k` is axis
//! `k`. `corner_offset` is not re-exported, so the order is **checked** instead:
//! every edge must join two corners differing exactly in the bit of its own
//! axis, which is true of the crate's numbering and of no other. Rule 5.

mod common;

use isomesh::Sdf;
use isomesh::fields::{ReferenceField, noise_cavity};
use isomesh::marching_cubes::FaceAmbiguity;
use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::interior::{Interior, SweptFaces};
use isomesh::marching_cubes::table::{
    AMBIGUOUS_FACES, EDGE_AXIS, EDGE_CORNERS, edge_on_face, is_inside,
};
use isomesh::marching_cubes::trilinear::Contours;

type Scalar = f64;

/// Samples per axis. 17 and 33 are the resolutions the crate's own prediction
/// test pins (30/64 under `Separate`, 8/40 under the decider), so they double as
/// a check that this reuses the predicate correctly.
const RESOLUTIONS: [u32; 4] = [17, 33, 49, 65];

/// Corner `c`'s offset, with the crate's bit order — asserted in `main`.
fn corner_offset(c: u8) -> [u32; 3] {
    [
        u32::from(c & 1),
        u32::from((c >> 1) & 1),
        u32::from((c >> 2) & 1),
    ]
}

/// Does any axis sweep of this cell report the positive regions joined through
/// the interior?
///
/// Returns `(joined, degenerate)` — the second counts sweeps whose bilinear
/// denominator vanishes, where the criterion is undefined and
/// [`SweptFaces::new`] refuses rather than guessing.
fn interior_join(corner: &[Scalar; 8]) -> (bool, usize) {
    let mut joined = false;
    let mut degenerate = 0usize;
    for w in 0..3usize {
        if let (Some(j), _) = interior_join_axis(corner, w) {
            joined |= j;
        } else {
            degenerate += 1;
        }
    }
    (joined, degenerate)
}

/// The same test on **one** sweep: the pair of faces perpendicular to `w`.
///
/// **Added after the first run, and it is the sharper question.** Taking the
/// disjunction over all three sweeps reports `Joined` on 100% of ambiguous-face
/// pairs — offenders and control alike — which falsifies P-17 by its own stated
/// falsifier and says nothing about the world, because a disjunction over three
/// chances is nearly always true. Custodio's criterion is about a *specific*
/// pair of opposite faces, and for an offending pair the one that matters is the
/// sweep whose endpoint faces include the **shared** face.
fn interior_join_axis(corner: &[Scalar; 8], w: usize) -> (Option<bool>, usize) {
    let u = (w + 1) % 3;
    let v = (w + 2) % 3;
    // `[A, B, C, D]` cyclic, `A`/`C` one diagonal: (0,0), (1,0), (1,1), (0,1).
    let face = |t: u8| -> [Scalar; 4] {
        let at = |du: u8, dv: u8| corner[usize::from((du << u) | (dv << v) | (t << w))];
        [at(0, 0), at(1, 0), at(1, 1), at(0, 1)]
    };
    match SweptFaces::new(face(0), face(1)) {
        Ok(sweep) => (Some(sweep.test() == Interior::Joined), 0),
        Err(_) => (None, 1),
    }
}

#[allow(dead_code, reason = "kept so the first run's predicate stays readable")]
fn interior_join_all_axes(corner: &[Scalar; 8]) -> (bool, usize) {
    let mut joined = false;
    let mut degenerate = 0usize;
    for w in 0..3usize {
        let u = (w + 1) % 3;
        let v = (w + 2) % 3;
        // `[A, B, C, D]` cyclic, `A`/`C` one diagonal: (0,0), (1,0), (1,1), (0,1).
        let face = |t: u8| -> [Scalar; 4] {
            let at = |du: u8, dv: u8| corner[usize::from((du << u) | (dv << v) | (t << w))];
            [at(0, 0), at(1, 0), at(1, 1), at(0, 1)]
        };
        match SweptFaces::new(face(0), face(1)) {
            Ok(sweep) => {
                if sweep.test() == Interior::Joined {
                    joined = true;
                }
            }
            Err(_) => degenerate += 1,
        }
    }
    (joined, degenerate)
}

/// One row.
struct Row {
    samples: u32,
    face_rule: &'static str,
    offending_pairs: usize,
    offending_with_interior_join: usize,
    control_pairs: usize,
    control_with_interior_join: usize,
    /// The same counts using only the sweep across the **shared** face.
    offending_with_axis_join: usize,
    control_with_axis_join: usize,
    degenerate_sweeps: usize,
}

fn measure(samples: u32, rule: FaceAmbiguity, rule_name: &'static str) -> Row {
    let field = noise_cavity::<Scalar>();
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);

    let cell_of = |x: u32, y: u32, z: u32| {
        let mut corner = [0.0f64; 8];
        let mut case = 0u8;
        for (c, slot) in corner.iter_mut().enumerate() {
            let o = corner_offset(c as u8);
            *slot = field.sample([
                lo[0] + h * f64::from(x + o[0]),
                lo[1] + h * f64::from(y + o[1]),
                lo[2] + h * f64::from(z + o[2]),
            ]);
            if is_inside(*slot) {
                case |= 1 << c;
            }
        }
        (case, corner)
    };
    let rings_of = |case: u8, corner: &[f64; 8]| {
        let mask = match rule {
            FaceAmbiguity::Separate => 0,
            FaceAmbiguity::AsymptoticDecider => joined_mask(corner, AMBIGUOUS_FACES[case as usize]),
        };
        let contours = Contours::of(case, mask);
        let mut owner = [255u8; 12];
        for r in 0..contours.count() {
            for &e in contours.ring(r) {
                owner[e as usize] = r as u8;
            }
        }
        owner
    };

    let mut row = Row {
        samples,
        face_rule: rule_name,
        offending_pairs: 0,
        offending_with_interior_join: 0,
        control_pairs: 0,
        control_with_interior_join: 0,
        offending_with_axis_join: 0,
        control_with_axis_join: 0,
        degenerate_sweeps: 0,
    };

    for axis in 0..3usize {
        for z in 0..samples - 1 {
            for y in 0..samples - 1 {
                for x in 0..samples - 1 {
                    let mut n = [x, y, z];
                    n[axis] += 1;
                    if n[axis] >= samples - 1 {
                        continue;
                    }
                    let (ca, va) = cell_of(x, y, z);
                    // Only pairs sharing an *ambiguous* face are in the
                    // population at all — that is what makes the control a
                    // control rather than a comparison against empty space.
                    if AMBIGUOUS_FACES[ca as usize] & (1u8 << (axis * 2 + 1)) == 0 {
                        continue;
                    }
                    let (cb, vb) = cell_of(n[0], n[1], n[2]);
                    let oa = rings_of(ca, &va);
                    let ob = rings_of(cb, &vb);
                    let cut = |owner: &[u8; 12], side: u8| -> Vec<u8> {
                        (0..12u8)
                            .filter(|&e| edge_on_face(e, axis, side) && owner[e as usize] != 255)
                            .collect()
                    };
                    let (cut_a, cut_b) = (cut(&oa, 1), cut(&ob, 0));
                    if cut_a.len() != 4 || cut_b.len() != 4 {
                        continue;
                    }
                    let one_cycle = |owner: &[u8; 12], edges: &[u8]| {
                        edges
                            .iter()
                            .all(|&e| owner[e as usize] == owner[edges[0] as usize])
                    };
                    let offending = one_cycle(&oa, &cut_a) && one_cycle(&ob, &cut_b);

                    let (ja, da) = interior_join(&va);
                    let (jb, db) = interior_join(&vb);
                    row.degenerate_sweeps += da + db;
                    // The sweep across the shared face, for each cell.
                    let (sa, _) = interior_join_axis(&va, axis);
                    let (sb, _) = interior_join_axis(&vb, axis);
                    let axis_join = sa.unwrap_or(false) || sb.unwrap_or(false);
                    if offending {
                        row.offending_pairs += 1;
                        if ja || jb {
                            row.offending_with_interior_join += 1;
                        }
                        if axis_join {
                            row.offending_with_axis_join += 1;
                        }
                    } else {
                        row.control_pairs += 1;
                        if ja || jb {
                            row.control_with_interior_join += 1;
                        }
                        if axis_join {
                            row.control_with_axis_join += 1;
                        }
                    }
                }
            }
        }
    }
    row
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    // The corner-bit order, asserted rather than assumed. Every edge joins two
    // corners differing exactly in the bit of its own axis.
    for e in 0..12usize {
        let [a, b] = EDGE_CORNERS[e];
        assert_eq!(
            b,
            a | (1 << EDGE_AXIS[e]),
            "edge {e} runs along axis {} but its corners are {a} and {b}, so bit k is not axis k \
             and every face this bench builds is transposed",
            EDGE_AXIS[e]
        );
    }

    let prereg = isomesh::experiment!("P-17");
    common::experiment::run(prereg, |run| {
        println!(
            "{:>5} {:<12} {:>10} {:>12} {:>12} {:>10} {:>12} {:>12} {:>11}",
            "n",
            "face rule",
            "offending",
            "any-axis",
            "shared-axis",
            "control",
            "any-axis",
            "shared-axis",
            "degenerate"
        );
        for samples in RESOLUTIONS {
            for (rule, name) in [
                (FaceAmbiguity::Separate, "separate"),
                (FaceAmbiguity::AsymptoticDecider, "decider"),
            ] {
                let r = measure(samples, rule, name);
                let share = |a: usize, b: usize| if b == 0 { 0.0 } else { a as f64 / b as f64 };
                println!(
                    "{:>5} {:<12} {:>10} {:>12.4} {:>12.4} {:>10} {:>12.4} {:>12.4} {:>11}",
                    r.samples,
                    r.face_rule,
                    r.offending_pairs,
                    share(r.offending_with_interior_join, r.offending_pairs),
                    share(r.offending_with_axis_join, r.offending_pairs),
                    r.control_pairs,
                    share(r.control_with_interior_join, r.control_pairs),
                    share(r.control_with_axis_join, r.control_pairs),
                    r.degenerate_sweeps
                );
                run.record(&[
                    ("samples", r.samples.to_string()),
                    ("face_rule", r.face_rule.to_string()),
                    ("offending_pairs", r.offending_pairs.to_string()),
                    (
                        "offending_with_interior_join",
                        r.offending_with_interior_join.to_string(),
                    ),
                    (
                        "control_with_interior_join",
                        r.control_with_interior_join.to_string(),
                    ),
                    ("control_pairs", r.control_pairs.to_string()),
                    (
                        "offending_share",
                        format!(
                            "{:.6}",
                            share(r.offending_with_interior_join, r.offending_pairs)
                        ),
                    ),
                    (
                        "control_share",
                        format!(
                            "{:.6}",
                            share(r.control_with_interior_join, r.control_pairs)
                        ),
                    ),
                    (
                        "offending_shared_axis_share",
                        format!(
                            "{:.6}",
                            share(r.offending_with_axis_join, r.offending_pairs)
                        ),
                    ),
                    (
                        "control_shared_axis_share",
                        format!("{:.6}", share(r.control_with_axis_join, r.control_pairs)),
                    ),
                    ("degenerate_sweeps", r.degenerate_sweeps.to_string()),
                ]);
            }
        }
        println!(
            "\n`offending` is shared ambiguous faces whose four cut edges lie in one cycle on both \
             sides —\nthe crate's own predicate, computed from the grid with no mesh involved. \
             `control` is every\nother ambiguous-face pair, which is what says whether the interior \
             join means anything."
        );
    });
}
