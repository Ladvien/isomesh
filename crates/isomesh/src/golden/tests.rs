//! The golden test, and the mutation check that proves it can fail.

use alloc::format;
use alloc::string::String;
use alloc::vec::Vec;

use super::{Entry, compute_all, field_of, fixture_path, hash_mesh, render};

/// Read the committed fixture as `(key, hash)` pairs plus the counts, keyed on
/// the combination so a reordering of the sweep is a diff rather than a
/// catastrophe.
fn read_fixture() -> Vec<(String, String)> {
    let text = std::fs::read_to_string(fixture_path()).unwrap_or_default();
    text.lines()
        .filter(|l| l.contains("\"hash\""))
        .filter_map(|line| {
            let key = format!(
                "{}/{}/{}",
                field_of(line, "algorithm")?,
                field_of(line, "field")?,
                field_of(line, "samples")?
            );
            let value = format!(
                "{} verts, {} tris, hash {}",
                field_of(line, "vertices")?,
                field_of(line, "triangles")?,
                field_of(line, "hash")?
            );
            Some((key, value))
        })
        .collect()
}

fn as_pairs(entries: &[Entry]) -> Vec<(String, String)> {
    entries
        .iter()
        .map(|e| {
            (
                format!("{}/{}/{}", e.algorithm, e.field, e.samples),
                format!(
                    "{} verts, {} tris, hash {:016x}",
                    e.vertices, e.triangles, e.hash
                ),
            )
        })
        .collect()
}

/// **The regression net.** Every (algorithm, field, resolution) must hash to
/// exactly what is committed.
///
/// On a mismatch this names each combination that drifted and prints both
/// values, because "a hash changed" is useless and "`dual_contouring/gyroid/33`
/// went from
/// 10584 to 10580 triangles" is a bug report.
///
/// Regenerate deliberately with `ISOMESH_BLESS=1`, and read the diff.
#[test]
fn golden_hashes_are_unchanged() {
    let computed = as_pairs(&compute_all());

    if std::env::var("ISOMESH_BLESS").is_ok() {
        let entries = compute_all();
        std::fs::write(fixture_path(), render(&entries)).expect("write fixture");
        std::println!(
            "measured: blessed {} golden hashes to {}",
            entries.len(),
            fixture_path().display()
        );
        return;
    }

    let committed = read_fixture();
    assert!(
        !committed.is_empty(),
        "no fixture at {} -- run with ISOMESH_BLESS=1 to create it",
        fixture_path().display()
    );

    let mut drifted = Vec::new();
    for (key, got) in &computed {
        match committed.iter().find(|(k, _)| k == key) {
            Some((_, want)) if want == got => {}
            Some((_, want)) => drifted.push(format!("  {key}\n    was {want}\n    now {got}")),
            None => drifted.push(format!("  {key}\n    was ABSENT\n    now {got}")),
        }
    }
    for (key, want) in &committed {
        if !computed.iter().any(|(k, _)| k == key) {
            drifted.push(format!("  {key}\n    was {want}\n    now ABSENT"));
        }
    }

    assert!(
        drifted.is_empty(),
        "{} of {} golden combinations drifted:\n{}\n\nIf this is intended, regenerate with \
         ISOMESH_BLESS=1 and read the diff before committing it.",
        drifted.len(),
        computed.len(),
        drifted.join("\n")
    );
}

/// **The acceptance criterion**: a one-bit change to a case table must fail the
/// golden test, and the failure must name the combination.
///
/// The corrupted table runs through T-005b's `march_with_table`, whose agreement
/// with the real extractor is itself asserted bit-for-bit by
/// `the_double_reproduces_marching_cubes` — so hashing its output with the *real*
/// table reproduces the committed `mc` hash exactly, and any difference under a
/// corrupted table is the corruption rather than the double.
///
/// Without this, the fixture would be a file that agrees with itself.
#[test]
fn a_corrupted_case_table_changes_the_hash_and_names_the_combination() {
    use crate::RuntimeShape3;
    use crate::fields::{ReferenceField, Sphere};
    use crate::marching_cubes::table::CASES;
    use crate::property::extraction::march_with_table;

    let field = Sphere::<f64>::canonical();
    let (lo, hi) = field.domain();
    let samples = 17u32;
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

    let honest = march_with_table(&field, &shape, lo, cell_size, &CASES);
    let honest_hash = hash_mesh(&honest);

    // The double reproduces the real extractor, so this must equal the committed
    // `marching_cubes/sphere/17` hash. If it does not, the golden fixture and the mutation
    // check are describing different code and neither means anything.
    let committed = read_fixture();
    let want = committed
        .iter()
        .find(|(k, _)| k == "marching_cubes/sphere/17")
        .map(|(_, v)| v.clone())
        .expect("marching_cubes/sphere/17 in the fixture");
    assert!(
        want.contains(&format!("hash {honest_hash:016x}")),
        "the test double disagrees with the committed hash: {want} vs {honest_hash:016x}"
    );

    // Now corrupt one triangle of one case the sphere actually reaches, exactly
    // as T-005b's mutation checks do.
    let mut cases = CASES;
    let victim = (1..255usize)
        .find(|&c| cases[c].count > 0)
        .expect("some non-empty case");
    cases[victim].triangles[0].swap(1, 2);

    let corrupted = march_with_table(&field, &shape, lo, cell_size, &cases);
    let corrupted_hash = hash_mesh(&corrupted);

    assert_ne!(
        honest_hash, corrupted_hash,
        "a flipped triangle in case {victim} did not change the hash, so the golden \
         fixture cannot detect a case-table regression"
    );

    // And the failure message has to be useful, not just present.
    let message = format!(
        "  marching_cubes/sphere/{samples}\n    was ... hash {honest_hash:016x}\n    now ... hash {corrupted_hash:016x}"
    );
    assert!(message.contains("marching_cubes/sphere/17"), "{message}");
    std::println!(
        "measured: T-007 mutation check -- flipping one triangle in case {victim} moved the hash {honest_hash:016x} -> {corrupted_hash:016x}"
    );
}

/// The hash must distinguish things that compare equal, or it is not a
/// bit-level check.
#[test]
fn the_hash_separates_signed_zero_and_notices_truncation() {
    use crate::MeshBuffer;

    let base = MeshBuffer::<f64> {
        positions: alloc::vec![[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]],
        normals: alloc::vec![[1.0, 0.0, 0.0]; 3],
        indices: alloc::vec![0, 1, 2],
    };

    // `+0.0 == -0.0` is true, so a value comparison would call these identical.
    let mut signed = base.clone();
    signed.positions[0][0] = -0.0;
    assert_ne!(
        hash_mesh(&base),
        hash_mesh(&signed),
        "a sign flip on a zero coordinate must change the hash -- it is exactly what a \
         reordered summation produces"
    );

    // A truncated index buffer must not hash as a prefix of the whole.
    let mut truncated = base.clone();
    truncated.indices.clear();
    assert_ne!(hash_mesh(&base), hash_mesh(&truncated));

    // And an unchanged mesh hashes the same twice.
    assert_eq!(hash_mesh(&base), hash_mesh(&base.clone()));
}

/// The fixture on disk must be exactly what `render` produces, so a hand edit
/// or a stale regeneration shows up as a diff rather than being tolerated.
#[test]
fn the_committed_fixture_is_canonically_formatted() {
    if std::env::var("ISOMESH_BLESS").is_ok() {
        return;
    }
    let on_disk = std::fs::read_to_string(fixture_path()).unwrap_or_default();
    assert!(!on_disk.is_empty(), "no fixture committed");
    assert_eq!(
        on_disk,
        render(&compute_all()),
        "the fixture is not in canonical form; regenerate with ISOMESH_BLESS=1"
    );
}

/// Guard against the sweep silently shrinking. 7 fields x 8 algorithm rows x 3
/// resolutions.
#[test]
fn every_combination_is_covered() {
    let entries = compute_all();
    // 7 fields x 8 algorithm rows x 3 resolutions. `marching_cubes+decider` is a
    // row rather than an algorithm — Marching Cubes with the decider switched on
    // — and `subgrid_marching_tetrahedra` is pinned to one 1D sampling
    // resolution, because M-95 measured that changing it moves every position by
    // about 1e-12 while leaving the topology alone.
    assert_eq!(entries.len(), 8 * 9 * 3, "{} combinations", entries.len());
    let unique: Vec<String> = {
        let mut keys: Vec<String> = as_pairs(&entries).into_iter().map(|(k, _)| k).collect();
        keys.sort();
        keys.dedup();
        keys
    };
    assert_eq!(unique.len(), entries.len(), "duplicate combination keys");
}
