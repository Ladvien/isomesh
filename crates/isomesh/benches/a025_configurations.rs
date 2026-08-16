//! **A-025 — every two-cell configuration that makes Manifold Dual Contouring
//! non-manifold, enumerated rather than sampled.**
//!
//! ```bash
//! cargo bench --bench a025_configurations
//! ```
//!
//! Writes `docs/measurements/a025-configurations.csv`.
//!
//! # Why exhaustive and not another census
//!
//! P-17 excluded the interior ambiguity (M-291) and left the mechanism unnamed.
//! The remaining question — *which* pairs — is finite and does not need a field:
//! two cells stacked along one axis share a face and have **twelve samples**
//! between them, so all `2¹² = 4,096` sign patterns fit in a loop. ✗17 settled
//! its own attribution exactly this way, and A-021's lesson is that a census
//! beats a minimisation when the defect is a local *predicate*.
//!
//! Sampling `noise_cavity` can only ever say which configurations that field
//! happens to reach. This says which exist.
//!
//! # The joined masks are enumerated too, so the decider arm is covered
//!
//! `joined_mask` needs corner *magnitudes*, not signs, so a sign enumeration
//! cannot evaluate the asymptotic decider directly. It does not need to: the
//! decider's only output is one bit per ambiguous face, so enumerating **every
//! reachable mask** on each cell covers every answer the decider could give,
//! and more besides. A configuration that offends under some reachable mask is
//! one the decider can reach; one that offends under none cannot.
//!
//! # The predicate is the crate's
//!
//! ```text
//! offending  ==  the shared face is ambiguous, and its four cut edges
//!                lie in one cycle on both sides
//! ```
//!
//! — `the_defect_count_is_predicted_from_the_grid_alone`, which M-291 already
//! reproduced against a mesh-derived count from the other direction.

mod common;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use isomesh::marching_cubes::table::{AMBIGUOUS_FACES, edge_on_face};
use isomesh::marching_cubes::trilinear::Contours;

/// The two cells are stacked along `z`, so the shared face is `z = 1`.
const AXIS: usize = 2;

/// Which cycle owns each of the twelve edges, or 255.
fn owners(case: u8, joined: u8) -> [u8; 12] {
    let contours = Contours::of(case, joined);
    let mut owner = [255u8; 12];
    for r in 0..contours.count() {
        for &e in contours.ring(r) {
            owner[e as usize] = r as u8;
        }
    }
    owner
}

/// Every mask the decider could produce for this case: one bit per ambiguous
/// face, all combinations.
fn reachable_masks(case: u8) -> Vec<u8> {
    let ambiguous = AMBIGUOUS_FACES[case as usize];
    let bits: Vec<u8> = (0..6u8).filter(|f| ambiguous & (1 << f) != 0).collect();
    let mut out = Vec::with_capacity(1 << bits.len());
    for combination in 0..(1u32 << bits.len()) {
        let mut mask = 0u8;
        for (k, bit) in bits.iter().enumerate() {
            if combination & (1 << k) != 0 {
                mask |= 1 << bit;
            }
        }
        out.push(mask);
    }
    out
}

/// Do the shared face's cut edges all lie in one cycle, and how many are there?
fn shared_face_cycle(owner: &[u8; 12], side: u8) -> Option<usize> {
    let cut: Vec<u8> = (0..12u8)
        .filter(|&e| edge_on_face(e, AXIS, side) && owner[e as usize] != 255)
        .collect();
    if cut.len() != 4 {
        return None;
    }
    let first = owner[cut[0] as usize];
    if cut.iter().all(|&e| owner[e as usize] == first) {
        Some(cut.len())
    } else {
        None
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    // Twelve samples: bits 0..3 are the z = 0 plane, 4..7 the shared z = 1
    // plane, 8..11 the z = 2 plane. Within a plane the bit is `(dy << 1) | dx`,
    // matching the crate's corner order with z fixed.
    let case_of = |pattern: u32, plane: u32| -> u8 {
        let mut case = 0u8;
        for c in 0..8u8 {
            let dz = u32::from((c >> 2) & 1);
            let low = u32::from(c & 3);
            if pattern & (1 << ((plane + dz) * 4 + low)) != 0 {
                case |= 1 << c;
            }
        }
        case
    };

    let mut offending: BTreeSet<(u8, u8)> = BTreeSet::new();
    let mut ambiguous_pairs: BTreeSet<(u8, u8)> = BTreeSet::new();
    let mut offending_patterns = 0usize;
    let mut ambiguous_patterns = 0usize;
    let mut separate_patterns = 0usize;
    let mut unavoidable_patterns = 0usize;
    let mut unavoidable: BTreeSet<(u8, u8)> = BTreeSet::new();
    let mut total_masks = 0usize;
    // How many cycles each cell has, for offenders and for the rest.
    let mut cycles_offending: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    let mut cycles_other: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    // Is the face opposite the shared one also ambiguous?
    let mut opposite_offending = [0usize; 4];
    let mut opposite_other = [0usize; 4];

    for pattern in 0..4096u32 {
        let ca = case_of(pattern, 0);
        let cb = case_of(pattern, 1);
        // The shared face is cell A's `+z` and cell B's `−z`.
        let shared_bit_a = 1u8 << (AXIS * 2 + 1);
        if AMBIGUOUS_FACES[ca as usize] & shared_bit_a == 0 {
            continue;
        }
        ambiguous_patterns += 1;
        ambiguous_pairs.insert((ca, cb));

        // Three questions, not one. `mask 0` is exactly what
        // `FaceAmbiguity::Separate` does and is sign-determined, so it is the
        // default's own answer. `any` is what the configuration space permits —
        // an upper bound on any face rule. `all` is the one that matters most:
        // a configuration offending under **every** mask cannot be fixed by any
        // face-disambiguation rule whatsoever, because the mask is the whole of
        // what such a rule chooses.
        //
        // **The two cells must agree about the shared face, and the first
        // version of this did not require that.** A face rule evaluates the
        // *face*, and both cells see the same four corners, so the bit is one
        // bit and not two — the asymptotic decider's whole point is that
        // neighbouring cells cannot disagree. Enumerating the masks
        // independently admits assignments no rule can produce, and "0 of 512
        // are unavoidable" measured that way is a statement about impossible
        // rules.
        let shared_a = 1u8 << (AXIS * 2 + 1);
        let shared_b = 1u8 << (AXIS * 2);
        let mut any = false;
        let mut all = true;
        let mut masks = 0usize;
        for ma in reachable_masks(ca) {
            for mb in reachable_masks(cb) {
                if (ma & shared_a != 0) != (mb & shared_b != 0) {
                    continue;
                }
                let oa = owners(ca, ma);
                let ob = owners(cb, mb);
                let bad =
                    shared_face_cycle(&oa, 1).is_some() && shared_face_cycle(&ob, 0).is_some();
                any |= bad;
                all &= bad;
                masks += 1;
            }
        }
        let zero = {
            let oa = owners(ca, 0);
            let ob = owners(cb, 0);
            shared_face_cycle(&oa, 1).is_some() && shared_face_cycle(&ob, 0).is_some()
        };
        if zero {
            separate_patterns += 1;
        }
        if all {
            unavoidable_patterns += 1;
            unavoidable.insert((ca, cb));
        }
        total_masks += masks;

        let na = Contours::of(ca, 0).count();
        let nb = Contours::of(cb, 0).count();
        let opposite_a = u8::from(AMBIGUOUS_FACES[ca as usize] & (1 << (AXIS * 2)) != 0);
        let opposite_b = u8::from(AMBIGUOUS_FACES[cb as usize] & (1 << (AXIS * 2 + 1)) != 0);
        let opposite = usize::from(opposite_a + opposite_b * 2);
        if any {
            offending_patterns += 1;
            offending.insert((ca, cb));
            *cycles_offending.entry((na, nb)).or_insert(0) += 1;
            opposite_offending[opposite] += 1;
        } else {
            *cycles_other.entry((na, nb)).or_insert(0) += 1;
            opposite_other[opposite] += 1;
        }
    }

    println!(
        "over all 4,096 two-cell sign patterns:\n  \
         {ambiguous_patterns} share an ambiguous face, of which **{offending_patterns}** can offend \
         under some reachable mask"
    );
    println!(
        "  distinct (case_a, case_b) pairs: {} sharing an ambiguous face, {} that can offend",
        ambiguous_pairs.len(),
        offending.len()
    );
    println!(
        "\n  under mask 0 — exactly `FaceAmbiguity::Separate`:     {separate_patterns} of \
         {ambiguous_patterns}"
    );
    println!(
        "  under **every** reachable mask — no face rule can help: {unavoidable_patterns} of \
         {ambiguous_patterns}  ({} distinct case pairs)",
        unavoidable.len()
    );
    println!(
        "  consistent mask combinations examined: {total_masks} — the two cells are required to \
         agree\n  about the shared face, because a face rule evaluates the face and both cells see \
         the same four corners"
    );

    println!("\ncycle counts (cell A, cell B) — offending against the rest:");
    let keys: BTreeSet<(usize, usize)> = cycles_offending
        .keys()
        .chain(cycles_other.keys())
        .copied()
        .collect();
    for k in keys {
        let o = cycles_offending.get(&k).copied().unwrap_or(0);
        let n = cycles_other.get(&k).copied().unwrap_or(0);
        println!(
            "  ({}, {}) {:>6} offending  {:>6} not   share {:.4}",
            k.0,
            k.1,
            o,
            n,
            if o + n == 0 {
                0.0
            } else {
                o as f64 / (o + n) as f64
            }
        );
    }

    println!("\nis the *opposite* face ambiguous too? (bit 0 = cell A's, bit 1 = cell B's)");
    for i in 0..4 {
        let (o, n) = (opposite_offending[i], opposite_other[i]);
        println!(
            "  {i} {:>6} offending  {:>6} not   share {:.4}",
            o,
            n,
            if o + n == 0 {
                0.0
            } else {
                o as f64 / (o + n) as f64
            }
        );
    }

    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/measurements");
    fs::create_dir_all(&dir).expect("create docs/measurements");
    let mut csv =
        String::from("# A-025: every two-cell configuration that can make MDC non-manifold\n");
    let _ = writeln!(csv, "case_a,case_b,cycles_a,cycles_b,offending,unavoidable");
    for (ca, cb) in &ambiguous_pairs {
        let _ = writeln!(
            csv,
            "{},{},{},{},{},{}",
            ca,
            cb,
            Contours::of(*ca, 0).count(),
            Contours::of(*cb, 0).count(),
            u8::from(offending.contains(&(*ca, *cb))),
            u8::from(unavoidable.contains(&(*ca, *cb)))
        );
    }
    let path = dir.join("a025-configurations.csv");
    fs::write(&path, csv).expect("write csv");
    println!("\n{} rows → {}", ambiguous_pairs.len(), path.display());
}
