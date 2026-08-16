//! **Field quality as first-class recorded numbers.**
//!
//! Ticket: T-017. This crate measures its *output* exhaustively — Hausdorff,
//! manifoldness, orientation, self-intersections, golden hashes — and until now
//! measured its *input* not at all. Every one of those output numbers is a
//! statement about an extractor applied to a field, and half of the pair went
//! unrecorded.
//!
//! ```bash
//! cargo bench --bench field_quality
//! ```
//!
//! Writes `docs/measurements/field_quality.csv`, which
//! `scripts/regress.sh field_quality` then diffs against its committed baseline
//! — so a field that silently degrades fails CI the same way a slower extractor
//! does.
//!
//! # The four numbers, and what each one catches
//!
//! - **`sup_grad`** — the largest `‖∇f‖` sampled over the field's own domain.
//!   Against a declared [`Lipschitz`](isomesh::fields::FieldBound::Lipschitz)
//!   constant this is the number that must not exceed it, and F-002 exists
//!   because a first draft declared `noise_cavity` at `2.598` when its gradient
//!   reaches `7.73` (M-244).
//! - **`eikonal_pct`** — share of samples with `‖∇f‖ ≈ 1`, the differential form
//!   of *"the value is the distance"*. An exact field should be near 100 and a
//!   gyroid nowhere near it.
//! - **`bound_gap`** — `sup_grad` divided by whatever the declaration allows.
//!   Above 1 means the declaration is **provably** wrong; well below 1 means it
//!   is loose, which costs a sphere tracer step count rather than correctness.
//! - **`certified_pct`** — T-015's share of active cells carrying Plantinga &
//!   Vegter's isotopy guarantee. It belongs here rather than with the mesh
//!   metrics because it is a property of the *field on this grid*: no extractor
//!   changes it.
//!
//! # Direction of error, stated once
//!
//! `sup_grad` is a **sampled maximum**, so it is a lower bound on a supremum. It
//! can prove a declaration wrong and can never prove one right — the same
//! asymmetry F-002's [`violates`](isomesh::validate::FieldBoundReport::violates)
//! encodes, and the reason `noise_cavity` is declared
//! [`Unbounded`](isomesh::fields::FieldBound::Unbounded) rather than given its
//! measured figure.
//!
//! # Why the columns are compared exactly
//!
//! Every one is a **max, a min, or a ratio of counts** over deterministic
//! `libm` arithmetic — never a sum, so no accumulation order to differ across
//! architectures. `regress.sh` therefore gets `None` for all of them: a
//! tolerance would only hide the change this file exists to catch.

use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use isomesh::fields::{FieldBound, ReferenceField};
use isomesh::validate::{field_bound_report, isotopy_report};
use isomesh::{RuntimeShape3, Sdf, Shape3};

/// Samples per axis for the gradient sweep. `n³` gradient evaluations.
const GRADIENT_SAMPLES: u32 = 32;

/// Grid sizes the isotopy certificate is evaluated at.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// Densities the gradient supremum is re-sampled at, to show it moving.
///
/// The point of the sweep is that `sup_grad` is a **sampled maximum**, so it can
/// only ever rise with density. A field whose figure is still climbing at the
/// last column has not been bounded, it has been under-sampled — which is a
/// claim worth making with four numbers instead of prose.
const DENSITY_SWEEP: [u32; 4] = [8, 16, 32, 64];

struct Row {
    field: &'static str,
    samples: u32,
    declared: &'static str,
    declared_bound: f64,
    sup_grad: f64,
    inf_grad: f64,
    eikonal_pct: f64,
    bound_gap: f64,
    certified_pct: f64,
}

/// The name and the numeric ceiling a declaration implies.
///
/// `Exact` and `Underestimate` both imply a Lipschitz constant of exactly one —
/// Corollary 1 of Bálint, Valasek & Gergó makes `1` the *smallest* Lipschitz
/// constant of a true distance, and this crate's underestimates are built from
/// exact operands by `min`/`max`, which preserve the constant while destroying
/// exactness. `Unbounded` claims nothing, so its ceiling is infinite and its gap
/// is meaningless rather than zero.
fn ceiling(bound: FieldBound) -> (&'static str, f64) {
    match bound {
        FieldBound::Exact => ("exact", 1.0),
        FieldBound::Underestimate { q } => {
            let _ = q;
            ("underestimate", 1.0)
        }
        FieldBound::Lipschitz { l } => ("lipschitz", l),
        FieldBound::Unbounded => ("unbounded", f64::INFINITY),
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!("field_quality — the input half of the pipeline\n");
    println!(
        "{:<16} {:>7} {:>14} {:>10} {:>10} {:>9} {:>11} {:>10} {:>13}",
        "field",
        "samples",
        "declared",
        "bound",
        "sup_grad",
        "inf_grad",
        "eikonal_%",
        "gap",
        "certified_%"
    );

    let mut rows: Vec<Row> = Vec::new();
    let mut density: Vec<(&'static str, [f64; 4])> = Vec::new();

    isomesh::for_each_reference_field!(f64, |name, field| {
        // Inline block, so no `return` in here (M-253).
        let mut sweep = [0.0f64; 4];
        for (slot, n) in sweep.iter_mut().zip(DENSITY_SWEEP) {
            *slot = field_bound_report(&field, n).sup;
        }
        density.push((name, sweep));

        let report = field_bound_report(&field, GRADIENT_SAMPLES);
        let (declared, bound) = ceiling(field.bound());
        let gap = report.sup / bound;

        let (lo, hi) = field.domain();
        for samples in RESOLUTIONS {
            let h = (hi[0] - lo[0]) / f64::from(samples - 1);
            let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
            let mut grid = Vec::with_capacity(shape.element_count());
            for z in 0..samples {
                for y in 0..samples {
                    for x in 0..samples {
                        grid.push(field.sample([
                            lo[0] + h * f64::from(x),
                            lo[1] + h * f64::from(y),
                            lo[2] + h * f64::from(z),
                        ]));
                    }
                }
            }
            let isotopy = isotopy_report(&grid, &shape).expect("isotopy report");
            let certified = 100.0 * isotopy.certified_fraction();

            println!(
                "{name:<16} {samples:>7} {declared:>14} {bound:>10.4} {:>10.5} \
                 {:>9.5} {:>11.3} {gap:>10.4} {certified:>12.3}%",
                report.sup,
                report.inf,
                100.0 * report.eikonal_fraction
            );
            rows.push(Row {
                field: name,
                samples,
                declared,
                declared_bound: bound,
                sup_grad: report.sup,
                inf_grad: report.inf,
                eikonal_pct: 100.0 * report.eikonal_fraction,
                bound_gap: gap,
                certified_pct: certified,
            });
        }

        // The gate, and it is one-sided on purpose. A sampled maximum above a
        // declared constant settles the question; one below it settles nothing.
        assert!(
            !report.violates(1e-6),
            "{name} declares {declared} and its gradient reaches {} — the \
             declaration is provably wrong",
            report.sup
        );
    });

    println!("\nsup ‖∇f‖ against sampling density — a maximum that is still climbing");
    println!("has not been bounded, it has been under-sampled:");
    println!(
        "{:<16} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "field", "n=8", "n=16", "n=32", "n=64", "still ↑?"
    );
    for (name, sweep) in &density {
        // "Still climbing" means the last step gained more than a thousandth,
        // which is far below any real feature and far above float noise on a
        // max of deterministic values.
        let climbing = sweep[3] > sweep[2] * (1.0 + 1e-3);
        println!(
            "{name:<16} {:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>10}",
            sweep[0],
            sweep[1],
            sweep[2],
            sweep[3],
            if climbing { "yes" } else { "no" }
        );
    }

    let mut csv = String::from(
        "field,samples,declared,declared_bound,sup_grad,inf_grad,eikonal_pct,\
         bound_gap,certified_pct\n",
    );
    for r in &rows {
        let _ = writeln!(
            csv,
            "{},{},{},{},{:.8},{:.8},{:.4},{},{:.4}",
            r.field,
            r.samples,
            r.declared,
            // An `Unbounded` declaration has no ceiling and therefore no gap.
            // Empty cells, not zeros: `regress.sh` skips an empty value, and a
            // literal `0` would read as "the gradient is nowhere near the
            // bound" — the opposite of "there is no bound".
            if r.declared_bound.is_finite() {
                std::format!("{:.6}", r.declared_bound)
            } else {
                String::new()
            },
            r.sup_grad,
            r.inf_grad,
            r.eikonal_pct,
            if r.declared_bound.is_finite() {
                std::format!("{:.6}", r.bound_gap)
            } else {
                String::new()
            },
            r.certified_pct
        );
    }
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements")
        .join("field_quality.csv");
    if let Some(dir) = path.parent() {
        let _ = fs::create_dir_all(dir);
    }
    let _ = fs::write(&path, csv);
    println!("\nwrote {}", path.display());

    report(&rows);
}

/// The reading, stated rather than left to whoever opens the CSV.
fn report(rows: &[Row]) {
    println!("\nwhat the numbers say:");

    let mut loose = 0;
    let mut seen = Vec::new();
    for r in rows {
        if seen.contains(&r.field) {
            continue;
        }
        seen.push(r.field);
        if r.declared_bound.is_finite() && r.bound_gap < 0.5 {
            loose += 1;
            println!(
                "  {:<16} declares {} but only reaches {:.3}× it — loose, which \
                 costs sphere-tracing steps rather than correctness",
                r.field, r.declared, r.bound_gap
            );
        }
    }
    if loose == 0 {
        println!("  every finite declaration is within 2× of its measured supremum");
    }

    // The certificate's trend with resolution, which is the field-quality
    // question the mesh metrics cannot answer: a field whose certified share
    // climbs is being resolved, one that is already at 100 has nothing to
    // resolve, and one that stalls has a feature no spacing will smooth.
    println!("\n  certified share, 17³ → 65³:");
    for field in &seen {
        let mine: Vec<&Row> = rows.iter().filter(|r| r.field == *field).collect();
        if mine.len() < 3 {
            continue;
        }
        let verdict = if mine[0].certified_pct >= 100.0 {
            "smooth at every spacing"
        } else if mine[2].certified_pct > mine[0].certified_pct {
            "resolving"
        } else {
            "stalled"
        };
        println!(
            "  {field:<16} {:>7.2}% → {:>7.2}% → {:>7.2}%  {verdict}",
            mine[0].certified_pct, mine[1].certified_pct, mine[2].certified_pct
        );
    }
}
