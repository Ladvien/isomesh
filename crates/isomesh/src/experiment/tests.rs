//! R-000's registry, checked. The *gate* is proved by the `compile_fail`
//! doctest on [`experiment!`](crate::experiment); these cover the properties a
//! doctest cannot see.

extern crate std;

use super::{PREREGISTERED, is_preregistered, preregistration};

/// Ids are unique, and in order.
///
/// Uniqueness because a duplicate would make `preregistration` return whichever
/// came first and silently attach an experiment to the wrong hypothesis. Order
/// because the file is read by humans and a list that drifts out of order gets
/// appended to in the wrong place.
#[test]
fn ids_are_unique_and_sorted() {
    let mut seen: std::vec::Vec<u32> = std::vec::Vec::new();
    for p in PREREGISTERED {
        let n: u32 =
            p.id.strip_prefix("P-")
                .and_then(|s| s.parse().ok())
                .unwrap_or_else(|| panic!("{} is not of the form P-<number>", p.id));
        assert!(!seen.contains(&n), "{} is registered twice", p.id);
        if let Some(last) = seen.last() {
            assert!(n > *last, "{} is out of order", p.id);
        }
        seen.push(n);
    }
}

/// Every registration is complete enough to be falsifiable.
///
/// **The falsifier is the field that matters.** A hypothesis with no stated
/// refutation is not a prediction, it is a description with a future tense, and
/// the whole reason `Preregistration` has a `falsified_by` field is to make
/// omitting one impossible rather than merely discouraged. The length floors
/// exist because `""` type-checks.
#[test]
fn every_registration_states_what_would_refute_it() {
    for p in PREREGISTERED {
        assert!(
            p.hypothesis.len() > 40,
            "{}: the hypothesis is too short to be a hypothesis",
            p.id
        );
        assert!(
            p.falsified_by.len() > 40,
            "{}: no falsifier worth the name",
            p.id
        );
        assert!(
            !p.records.is_empty(),
            "{}: an experiment that records nothing cannot be checked",
            p.id
        );
        // The *shape* of a ticket id, not a list of series that run experiments.
        // This was `starts_with('R') || 'A' || 'T'` and went stale the first time
        // an `S` ticket registered one (S-009/P-19) -- Part 5's rule about
        // hand-written bounds, third instance. The shape here is the one
        // `scripts/backlog_gate.sh` uses: `[A-Za-z]+-[0-9]+[a-z]?`.
        let well_formed = match p.ticket.split_once('-') {
            Some((series, number)) => {
                !series.is_empty()
                    && series.chars().all(|c| c.is_ascii_alphabetic())
                    && number.chars().next().is_some_and(|c| c.is_ascii_digit())
                    && number.chars().all(|c| c.is_ascii_alphanumeric())
            }
            None => false,
        };
        assert!(well_formed, "{}: `{}` is not a ticket id", p.id, p.ticket);
    }
}

/// Column names are unique within a registration.
///
/// A duplicate would make `Run::record`'s completeness check pass while one of
/// the two columns went unwritten — the exact silence that check exists to
/// prevent.
#[test]
fn records_are_unique_within_a_registration() {
    for p in PREREGISTERED {
        for (i, a) in p.records.iter().enumerate() {
            for b in &p.records[i + 1..] {
                assert_ne!(a, b, "{}: `{a}` is listed twice", p.id);
            }
        }
    }
}

/// Lookup agrees with membership, both ways.
#[test]
fn lookup_and_membership_agree() {
    for p in PREREGISTERED {
        assert!(is_preregistered(p.id));
        assert_eq!(preregistration(p.id).hypothesis, p.hypothesis);
    }
    assert!(!is_preregistered("P-999"));
    assert!(!is_preregistered("p-8"), "the check must be case sensitive");
    assert!(!is_preregistered("P-8 "), "trailing space must not match");
    assert!(!is_preregistered(""), "the empty id must not match");
}

/// `is_preregistered` really is const-evaluable.
///
/// The macro's whole gate is a `const` assertion, so if this ever stopped
/// folding at compile time the check would silently become a run-time one — and
/// a run-time check on a registry that is a compile-time constant would never
/// fire, because the experiment would already be running.
#[test]
fn the_check_folds_at_compile_time() {
    // **The assertion is the compilation.** These two items are evaluated by
    // the compiler; if `is_preregistered` ever stopped being const-evaluable
    // they would not build, and if either answer were wrong the `assert!` would
    // fire during evaluation. Nothing at run time is needed or wanted -- a
    // run-time check on a compile-time constant would never fire, because by
    // then the experiment is already running.
    const _YES: () = assert!(is_preregistered("P-8"));
    const _NO: () = assert!(!is_preregistered("P-999"));
}
