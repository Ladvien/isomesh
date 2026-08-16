//! Tests for the extractor registry.
//!
//! The load-bearing one is [`every_extractor_impl_is_registered_or_excused`]:
//! X-001's whole point is that adding a ninth algorithm is one edit, and the
//! failure it prevents is silent — an `impl` written and never registered
//! compiles, warns about nothing, and leaves a bench measuring seven algorithms
//! while its header says eight.

use super::{ALL_EXTRACTORS, Extractor, UNREGISTERED};
use crate::fields::Sphere;
use crate::{MeshBuffer, RuntimeShape3};

/// **Every `Extractor` impl is either in the registry or excused by name.**
///
/// This is X-001's acceptance, and it is a *source* check because Rust has no
/// runtime reflection: there is no way to ask "how many types implement this
/// trait". Both halves live in `extractor.rs` — the impls come from one
/// `forward_extractor!` invocation and the registry from one
/// `for_each_extractor!` definition — so this reads the file it lives beside and
/// compares the two lists.
///
/// **A count alone would pass for the wrong reason**, which is why this matches
/// names rather than lengths. The registry has seven entries and there are seven
/// impls, and those two sevens are unrelated: `MarchingCubes` appears in the
/// registry *twice* — once bare, once with the asymptotic decider — and
/// `GreedyQuads` appears not at all. Two errors cancelling is exactly what an
/// `assert_eq!(7, 7)` cannot see.
#[test]
fn every_extractor_impl_is_registered_or_excused() {
    let source = include_str!("../extractor.rs");

    // The types given to `forward_extractor!`, which is the complete set of
    // `Extractor` impls by construction — there is no other way to get one.
    let (_, after) = source
        .split_once("forward_extractor!(")
        .expect("the impl list moved");
    let (list, _) = after
        .split_once(");")
        .expect("the impl list is unterminated");
    let impls: alloc::vec::Vec<&str> = list
        .lines()
        .map(str::trim)
        .filter(|line| line.starts_with("crate::"))
        .map(|line| {
            // `crate::marching_cubes::MarchingCubes<R>,` -> `MarchingCubes`
            let path = line.trim_end_matches(',');
            let leaf = path.rsplit("::").next().unwrap_or(path);
            leaf.split('<').next().unwrap_or(leaf)
        })
        .collect();
    assert!(
        impls.len() >= 7,
        "only {} impls parsed, so the parse broke rather than the code: {impls:?}",
        impls.len()
    );

    // The registry's body, where a registered type is named by its constructor.
    //
    // **Bounded at both ends on purpose.** Taking everything after the macro
    // swallows `UNREGISTERED`, whose whole job is to name the types that are
    // *not* registered — so every excused type appeared registered and the
    // check passed while testing nothing. Caught by running it.
    let (_, after_macro) = source
        .split_once("macro_rules! for_each_extractor")
        .expect("the registry moved");
    let (macro_body, _) = after_macro
        .split_once("pub const ALL_EXTRACTORS")
        .expect("the registry's end marker moved");

    for name in &impls {
        let registered = macro_body.contains(name);
        let excused = UNREGISTERED.contains(name);
        assert!(
            registered ^ excused,
            "{name} is {} — an Extractor impl must be either in for_each_extractor! \
             or in UNREGISTERED with a reason, and never both",
            if registered {
                "both registered and excused"
            } else {
                "neither registered nor excused"
            }
        );
    }
    std::println!(
        "measured: {} Extractor impls, {} registered, {} excused",
        impls.len(),
        impls.len() - UNREGISTERED.len(),
        UNREGISTERED.len()
    );
}

/// **The macro and the name list visit the same entries, in the same order.**
///
/// `ALL_EXTRACTORS` exists so a reader can see the list without expanding a
/// macro. A list that drifts from the thing it describes is worse than no list,
/// so it is checked rather than trusted.
#[test]
fn the_registry_and_the_macro_agree() {
    let mut visited = alloc::vec::Vec::new();
    crate::for_each_extractor!(f64, |name, extractor| {
        // Touch it so the binding is not merely declared: a registry entry that
        // does not construct is a registry entry that does not run.
        let _ = &mut extractor;
        visited.push(name);
    });
    assert_eq!(
        visited.as_slice(),
        ALL_EXTRACTORS.as_slice(),
        "the macro and ALL_EXTRACTORS disagree"
    );
}

/// **Every registered extractor meshes a sphere to something.**
///
/// The cheapest possible property, and it is here because the registry's
/// entries are *expressions*: one that fails to configure, or configures itself
/// into producing nothing, is a silent hole in every sweep that enumerates from
/// this list. The subgrid entry is the one that can genuinely fail to
/// construct — its resolution has to be positive — so this is where that is
/// caught rather than in whichever bench runs first.
#[test]
fn every_registered_extractor_meshes_a_sphere() {
    let field = Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new([17; 3]).expect("valid shape");
    let mut out = MeshBuffer::<f64>::new();
    let mut rows = alloc::vec::Vec::new();

    crate::for_each_extractor!(f64, |name, extractor| {
        out.reset();
        extractor
            .extract_into(&field, &shape, [-2.0; 3], 0.25, &mut out)
            .expect("extraction");
        assert!(
            out.triangle_count() > 0,
            "{name} meshed a sphere to nothing"
        );
        rows.push((name, out.triangle_count()));
    });

    assert_eq!(rows.len(), ALL_EXTRACTORS.len());
    std::println!("measured: sphere at 17^3 -> {rows:?}");
}
