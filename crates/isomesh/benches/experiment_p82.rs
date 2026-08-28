//! **P-82 — tunnelling through a wall you just made thin.**
//!
//! Ticket: R-082. Pre-registered before this harness existed; the registration
//! in `isomesh::experiment` is the contract and is not edited here.
//!
//! ```bash
//! cargo bench --bench experiment_p82
//! ```
//!
//! Writes `docs/experiments/p-82.csv`.
//!
//! # The source, read rather than paraphrased
//!
//! Pelletier-Guénette, Mercier-Aubin and Andrews, *Real-Time Triangle-SDF
//! Continuous Collision Detection*, `10.1145/3747862` (PACMCGIT 8(4) Article 49,
//! SCA 2025). The method implemented below is their **Algorithm 2, FWGSS**:
//! a modified Frank-Wolfe over the barycentric coordinates *and* time, with a
//! golden-section line search (their Algorithm 1, after Kiefer 1953) solving each
//! one-dimensional sub-problem, and `φ(x_t) ≤ 0` used as a new **temporal upper
//! bound** rather than as a termination. Their §4.6 broad phase — a bounding
//! sphere marched against the field — is here too, and it is load-bearing:
//!
//! **The registration's two quoted figures only reconcile through it.** Table 4
//! reads `Shuriken / FWGSS: Total 0.4 ms, Tri. Tested 36, Mean 11.13 µs`. The
//! shuriken *mesh* has 888 triangles; **36** of them reached the narrow phase.
//! `36 × 11.13 µs = 0.40 ms` exactly, and `0.4 ms / 888 = 0.45 µs` — the naive
//! reading of the registration's sentence — is not a figure in the paper. The
//! `0.96 µs/tri` at 100K is Table 3's per-triangle mean on the armadillo, whose
//! **total is 86.85 ms**. So this harness reports cost per element *tested* (the
//! paper's own unit) and the broad-phase cull beside it, and never divides a
//! total by a mesh size.
//!
//! Two more numbers from the same tables, because they set the bar C2 is
//! measured against: FWGSS's per-triangle mean is **25.57 µs on Skate**, 7.93 on
//! Bunny, 11.13 on Shuriken, against Macklin et al. 2020's DCD at 30.69–67.97.
//! The paper's own method is *over* the registered 25 µs bar on one of its three
//! scenes. And Table 3 is **sub-linear in triangle count** — 6.47, 24.18, 86.85
//! ms for 1K, 10K, 100K, so 100× the triangles for 13.4× the time — because their
//! time interval shrinks whenever an earlier contact is found, which makes every
//! later triangle a cheaper problem. That is an inter-element optimisation, so it
//! is measured here as its own arm rather than folded in silently.
//!
//! # SHARE, recomputed before the harness was written
//!
//! The registration's SHARE line is: *"C2 moves a single dynamic body's collision
//! budget, not the whole frame."*
//!
//! - **C2's 25 µs is absolute, not a ratio**, so `✗51`'s Amdahl obstruction
//!   cannot apply to it. What the share does is bound the *element count*: 25 µs
//!   per element inside a 10% slice of `game_dig`'s own 16 ms frame buys
//!   **64 elements**. `game_dig`'s body is 4 spheres (0.625% of the frame at the
//!   bar); an 8-sphere capsule proxy is 1.25%; a 200-triangle capsule mesh is
//!   **31.25% of the whole frame** and is not a single body's budget at all.
//! - **C3's 10× is a ratio between two whole measurements**, so its share is 1.0
//!   and `1/(1 − 1/10) = 1.11` — no obstruction. The element counts alone give
//!   200/8 = **25×** before any per-element difference, so the clause has
//!   headroom and could only fail if a triangle test were 2.5× *cheaper* than a
//!   sphere test.
//! - **C1 is a count, not a ratio.**
//!
//! ## But C1 is arithmetically unreachable at `game_dig`'s own body, and that is
//! ## stated here before the run
//!
//! The registration writes the tunnelling threshold as `t < v·Δt`, which is the
//! **point-particle** threshold. A moving *sphere* of radius `R` is caught by a
//! discrete test whenever either endpoint of a frame lies in the overlap window,
//! and that window is `2R + t` wide, so the discrete path tunnels only when
//!
//! ```text
//! v·Δt > 2R + t.
//! ```
//!
//! `game_dig`'s body radius is `BODY_RADIUS = 0.25` = **2 cells**, so `2R` alone
//! is **4 cells**, while the fastest thing in the demo moves **2.459 cells** in a
//! 16 ms frame (sprint 9.0 u/s composed with a full-sandbox-height fall,
//! `sqrt(9² + 2·18·8) = 19.209 u/s`). **The discrete arm cannot tunnel at
//! `game_dig`'s speeds with `game_dig`'s body, by arithmetic, at any wall
//! thickness in the registered sweep.** That is the registration's own falsifier
//! for C1 — "the discrete path not tunnelling, which would mean the game's speeds
//! are below the threshold and this is premature" — and it is computed rather than
//! discovered, the `✗54` pattern.
//!
//! The run happens anyway, and the sweep is widened in the two directions that
//! turn a predicted zero into a located threshold: **element radius** (a point
//! projectile, a quarter-cell, one cell, and `game_dig`'s two cells) and **speed**
//! (four speeds the demo actually reaches, plus 40, 80 and 160 u/s, which no
//! entity in the demo has). Those extra rows are what make the registered
//! VACUITY CONTROL non-vacuous, and they answer the engineering question the
//! clause was reaching for: *above what projectile speed does discrete collision
//! stop working for this body?*
//!
//! # The vacuity control, quoted from the registration
//!
//! > VACUITY CONTROL: the discrete arm must tunnel at least once, reported as a
//! > count, or C1 is comparing two methods that both work.
//!
//! Emitted as `tunnels_discrete` per row and `tunnels_discrete_total` on every
//! row, and asserted `> 0` after the sweep.
//!
//! Two further controls, because a tunnel count is only meaningful if the shots
//! were aimed at solid rock:
//!
//! - **`shots_crossing_solid` must equal `shots`, asserted.** Every shot is
//!   generated so that its centre path passes through the wall's mid-plane, and
//!   the field is sampled *there* to prove the point is inside the solid. A
//!   tunnel over a shot that missed the wall is `M-44`'s zero wearing a one's
//!   clothes.
//! - **`prediction_mismatches` must be 0 on the exact-slab arm, asserted.** For a
//!   slab the discrete predicate is `|x| − t/2 < R`, so which frames detect is
//!   closed-form. The harness computes that from the geometry and compares it to
//!   what sampling the field actually said. This is `P-72`'s duplication-factor
//!   control in a different costume: a harness that asserts arithmetic against
//!   itself has asserted nothing, so the arithmetic is asserted against the
//!   instrument.
//!
//! # Fixtures, and the two constants lifted from the demo
//!
//! **`bevy_isomesh` is not a dependency of this crate** (`CLAUDE.md` hard rule 2),
//! so the demo's numbers are **restated as constants with the source line named**
//! rather than imported. They are: `CELL_SIZE = 0.125`
//! (`bevy_isomesh/examples/game_dig.rs:127`), the fixed frame
//! `Duration::from_millis(16)` (`:3060`, the demo's own test frame — the shipped
//! app reads `time.delta_secs()`), walk `2.5` and sprint `9.0` (`:2487-2491`),
//! `GRAVITY = 18.0` (`:304`), `BODY_RADIUS = 0.25` (`:289`),
//! `BODY_OFFSETS = [0.25, 0.65, 1.05, 1.45]` (`:301`, so the body's sphere
//! centres span `1.2`), `WALL_THICKNESS = 0.5` (`:148`) and the default brush
//! radius `0.25` (`:287-289`).
//!
//! **There is no projectile in `game_dig`.** The demo edits the field through a
//! ray-cast aim, and its only moving body is the player. So "the speed a
//! `game_dig` projectile travels" is read as the fastest speed a `game_dig` body
//! reaches, and the two composed maxima are computed rather than guessed.
//!
//! Three fields, all exact-or-conservative and all **1-Lipschitz**, which is what
//! makes conservative advancement sound with `L = 1` and no fudge factor:
//!
//! - `slab` — the registered fixture, thickness `t`, exact distance.
//! - `dug` — `game_dig`'s own 0.5-unit sandbox wall with a default-radius brush
//!   subtracted from behind so the residual thickness on the axis is exactly `t`.
//!   `max(dA, −dB)` is the intersection of two exact fields, so outside the solid
//!   it *under*-states the distance, which is the safe direction for advancement.
//! - `rim` — a convex nub of 2 cells, the shape of a dug pit's edge. C3 needs it:
//!   see below.
//!
//! # C3's clause is nearly unfalsifiable against a flat wall, and why
//!
//! The support function of a union of spheres is the max of the spheres' support
//! functions. The 8 bead centres include **both segment endpoints**, so
//!
//! ```text
//! sup_beads(n̂) = |h·(â·n̂)| + r = sup_capsule(n̂)   for every n̂,
//! ```
//!
//! exactly. **A swept-sphere proxy is support-exact against a plane**, so on the
//! registered wall fixture the ToI disagreement is not the proxy's error at all —
//! it is the *triangle mesh's* inscribed-polyhedron deficit, `r(1 − cos(π/10)) =
//! 0.0979 cells` for 10 longitudes. Asserting "within one cell" there is
//! asserting a tenth of a cell of faceting, and the proxy could not have failed.
//! That is exactly `P-70`'s C3 — a HELD with no instrument — so the wall arm is
//! reported *with* the analytic support deficits that prove the mechanism, and a
//! second fixture is added where the clause can fail: a **convex rim**, where the
//! scallop troughs between beads are real gaps. The capsule's half-height is swept
//! from `game_dig`'s own 0.6 up to 3.6, which crosses `h = 7r = 1.75` — above
//! that the 8 beads are **disjoint** and the proxy has holes a rim can enter.
//!
//! # Deviations from the paper, stated
//!
//! - **Fixed iteration counts.** Algorithm 1 loops on a tolerance and Algorithm 2
//!   breaks on convergence; both are given fixed counts here, so cost is
//!   data-independent and the numbers are comparable across rows. The cost of the
//!   deviation is measured: `solver_error_cells` is the solver's ToI error against
//!   the closed-form answer on the slab, and the iteration counts are swept.
//! - **No rotation.** `game_dig`'s body does not rotate, so trajectories are
//!   linear, `v_t` of Equation 14 is the constant step, and §4.6's non-linear
//!   padding is zero. The broad phase is therefore *exact* here, which flatters
//!   the cull rate relative to the paper's rigid-body scenes.
//! - **No adaptive triangle subdivision (§4.5).** It generates *additional*
//!   simultaneous contacts on a coarse mesh; this experiment asks only for the
//!   first time of impact, which is one contact.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::print_literal,
    clippy::too_many_lines,
    clippy::similar_names
)]

mod common;

fn main() {
    // A bench target is also built by `cargo test`, which runs it with no
    // arguments. Doing the whole sweep there would make every test run pay for
    // it, so the work happens only under `cargo bench`.
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-82");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    #[cfg(not(target_os = "linux"))]
    {
        // M-280: on a governed CPU a nanosecond is not a unit, so every timed row
        // here carries `ghz` computed as cycles ÷ nanoseconds. That needs
        // `perf_event_open`, which is Linux-only. Refusing beats inventing.
        eprintln!(
            "{} reports cycles and a measured clock on every timed row (M-280), and that needs\n\
             hardware performance counters. This platform has no `perf_event_open`.\n\
             Run it on Linux; `perf_event_paranoid = 2` is permissive enough and no root is needed.",
            prereg.id
        );
        std::process::exit(1);
    }
}

#[cfg(target_os = "linux")]
mod experiment {
    use std::hint::black_box;
    use std::time::Instant;

    use isomesh::Sdf;

    use crate::common::counters::{MIN_TIME_RATIO, Probe};
    use crate::common::experiment::Run;

    // ── constants restated from bevy_isomesh/examples/game_dig.rs ─────────────
    // `bevy_isomesh` is outside this crate's dependency graph, so these are
    // restated with their source line rather than imported.

    /// `game_dig.rs:127` — `CELL_SIZE`. A power of two, which is what makes the
    /// demo's chunk seams bit-exact (`M-32`).
    const CELL: f64 = 0.125;
    /// `game_dig.rs:3060` — `FRAME: Duration = Duration::from_millis(16)`, the
    /// demo's own fixed test frame. The shipped app integrates
    /// `time.delta_secs()`, so this is the demo's stated frame rather than a
    /// measured one.
    const DT: f64 = 0.016;
    /// `game_dig.rs:2490` — the walk speed.
    const WALK: f64 = 2.5;
    /// `game_dig.rs:2488` — the sprint speed, with its own comment: "a 16-unit
    /// box takes 6.4 s to cross at the walk speed and 1.8 s at this one".
    const SPRINT: f64 = 9.0;
    /// `game_dig.rs:304` — `GRAVITY`, "roughly twice Earth's".
    const GRAVITY: f64 = 18.0;
    /// `game_dig.rs:137` — `EXTENT = [8, 4, 8]` chunks of
    /// `CHUNK_CELLS * CELL_SIZE = 2.0` units, so the sandbox is 8 units tall and
    /// a body falling its full height arrives at `sqrt(2 g H)`.
    const SANDBOX_HEIGHT: f64 = 8.0;
    /// `game_dig.rs:289` — `BODY_RADIUS`. Two cells.
    const BODY_RADIUS: f64 = 0.25;
    /// `game_dig.rs:301` — `BODY_OFFSETS = [0.25, 0.65, 1.05, 1.45]`, so the four
    /// sphere centres span `1.45 − 0.25 = 1.2` units. Half of that is the
    /// capsule half-height C3's proxy has to cover.
    const BODY_SPAN: f64 = 1.2;
    /// `game_dig.rs:148` — `WALL_THICKNESS`, the five slabs lining the sandbox.
    /// Four cells, and the thing a digger is standing in front of.
    const WALL: f64 = 0.5;
    /// `game_dig.rs:287-289` — the default brush radius, two cells.
    const BRUSH: f64 = 0.25;

    // ── the registered sweeps ────────────────────────────────────────────────

    /// Wall thickness in cells: from the registered top of two cells down to
    /// `t/h = 0.05`, which is `subgrid`'s floor — `M-95` measured a slab
    /// **1/20 of a cell thick** returning 0 triangles at 2 samples per edge and a
    /// mesh at 256.
    const THICKNESS_CELLS: [f64; 8] = [2.0, 1.5, 1.0, 0.75, 0.5, 0.25, 0.125, 0.05];

    /// Radius of the moving element, in cells.
    ///
    /// `0.0` is the point particle the registration's `t < v·Δt` threshold
    /// implicitly assumes; `2.0` is `game_dig`'s own `BODY_RADIUS`. The two in
    /// between locate the threshold rather than bracketing it.
    const RADII_CELLS: [f64; 4] = [0.0, 0.25, 1.0, 2.0];

    /// Randomised shots per row. The registered `10^4`.
    const SHOTS: u32 = 10_000;

    /// Half-angle of the cone the shot direction is drawn from, radians.
    ///
    /// 45°, not 90°: past that the shot is more parallel to the wall than across
    /// it, and a "shot at a wall" that mostly travels along the wall is not the
    /// configuration the registration describes.
    const CONE: f64 = std::f64::consts::FRAC_PI_4;

    /// Lateral spread of the shots about the wall's axis, in cells.
    ///
    /// Half a cell, so that on the `dug` field every shot sees a residual
    /// thickness within `0.063` cells of the nominal `t` — the pocket is a sphere
    /// and the residual grows off-axis.
    const RHO_CELLS: f64 = 0.5;

    /// Interval-halving depth of the sphere solver's branch and bound. The leaf
    /// resolution is `|step| · 2^-40`, which is `2.3e-12` world units at the
    /// fastest speed swept.
    const MAX_DEPTH: u32 = 40;

    /// Hard cap on nodes visited by one sphere query. Asserted, never used as a
    /// fallback: a query that reaches it is a query whose geometry the Lipschitz
    /// bound could not separate, and that has to stop the experiment rather than
    /// round to "no collision".
    const MAX_NODES: u32 = 4096;

    /// Contact tolerance for the sphere solver, world units. `1e-12` is a
    /// thousandth of a millionth of a cell and well inside `f64`'s reach on
    /// quantities of order 1.
    const CONTACT_TOL: f64 = 1e-12;

    /// Element counts for C2's sweep. `2 · longitudes · rings` for the triangle
    /// arm, so every count is exactly realisable: `200 = 2·10·10` is C3's capsule
    /// and `888 = 2·12·37` is the paper's shuriken.
    const ELEMENT_COUNTS: [usize; 8] = [8, 16, 40, 100, 200, 400, 888, 1600];

    /// `(longitudes, rings)` giving each of [`ELEMENT_COUNTS`].
    const TESSELLATIONS: [(usize, usize); 8] = [
        (4, 1),
        (4, 2),
        (5, 4),
        (5, 10),
        (10, 10),
        (10, 20),
        (12, 37),
        (20, 40),
    ];

    /// FWGSS iteration budgets: `(Algorithm 2 iterations, Algorithm 1 iterations)`.
    const SOLVERS: [Solver; 3] = [
        Solver {
            fw: 4,
            gss: 6,
            name: "fw4_gss6",
        },
        Solver {
            fw: 8,
            gss: 10,
            name: "fw8_gss10",
        },
        Solver {
            fw: 16,
            gss: 20,
            name: "fw16_gss20",
        },
    ];

    /// The bar a solver setting must clear before its cost is allowed to answer
    /// C2: a hundredth of a cell of ToI error against the closed form. Two orders
    /// inside C3's one-cell tolerance, so an accepted setting cannot be trading
    /// accuracy for the 25 µs bar.
    const ACCURACY_BAR_CELLS: f64 = 0.01;

    /// Capsule half-heights for C3, world units. `0.6` is `game_dig`'s own body;
    /// the sweep crosses `7 · BODY_RADIUS = 1.75`, above which eight beads
    /// spaced `2h/7` apart no longer overlap.
    const HALF_HEIGHTS: [f64; 6] = [0.6, 1.2, 1.8, 2.4, 3.0, 3.6];

    /// Capsule orientations per C3 row-group: three canonical, thirteen random.
    const ORIENTATIONS: usize = 16;

    /// Where in the step the **true capsule** first touches, in every C3 fixture.
    ///
    /// A tenth, not a half. The disagreement C3 measures is a difference of
    /// support, and a representation that falls short of the capsule by more than
    /// the *remaining* travel cannot touch the wall inside one step at all — so the
    /// fixture's dynamic range is `(1 − TAU_PLACE)` steps, and it has to exceed
    /// C3's one-cell bar or the clause cannot be falsified.
    const TAU_PLACE: f64 = 0.1;

    /// Step length of C3's second wall arm, in cells.
    ///
    /// `game_dig`'s sprint step is 1.152 cells, so `0.9 × 1.152 = 1.037` cells is
    /// all the range the game's own speed leaves — a one-cell bar with 3.7% of
    /// headroom. C3 names no speed, so a second arm at eight cells per step gives
    /// the bar **7.2 cells** of room to be exceeded.
    const WIDE_STEP_CELLS: f64 = 8.0;

    /// Wall thickness of the wide-step arm, in cells. Larger than the step, so the
    /// body cannot pass clean through and the unsigned-distance line search has
    /// only one crossing to find — the paper's Figure 2 failure, foreclosed by the
    /// fixture rather than by hope.
    const WIDE_WALL_CELLS: f64 = 16.0;

    /// Timing repetitions; the median is taken. Three is `P-72`'s figure and the
    /// reason is the same — the quantity being compared is a ratio, and one
    /// sample of it on a governed CPU is `M-337`'s mistake.
    const REPS: usize = 5;

    /// The sentinel for a column this row does not measure. `P-71`'s convention:
    /// `record` requires every registered key to be present, and an empty value
    /// says "not this arm" without inventing a number.
    const NA: &str = "";

    type Row = Vec<(&'static str, String)>;

    /// `(element kind, solver budget, the paper's inter-element interval shrink)`
    /// — the group a linearity fit belongs to.
    type FitKey = (&'static str, &'static str, bool);

    // ── fields ───────────────────────────────────────────────────────────────

    /// The three fixtures. Every variant is 1-Lipschitz, which is the whole
    /// reason conservative advancement below can use `L = 1` with no safety
    /// factor: `max` of exact distance fields is 1-Lipschitz, and outside the
    /// solid it under-states the true distance, which is the safe direction.
    #[derive(Clone, Copy)]
    enum Scene {
        /// The registered fixture: the slab `[lo, hi]` on x, unbounded in y and z.
        Slab { lo: f64, hi: f64 },
        /// `game_dig`'s own sandbox wall, dug from behind by a default-radius
        /// brush until the residual thickness on the axis is `hi − lo`'s target.
        Dug {
            lo: f64,
            hi: f64,
            c: [f64; 3],
            r: f64,
        },
        /// A convex nub: the rim of a dug pit, or a rock standing proud of one.
        Rim { c: [f64; 3], r: f64 },
    }

    fn dist(p: [f64; 3], c: [f64; 3]) -> f64 {
        let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    }

    impl Sdf for Scene {
        type Scalar = f64;

        fn sample(&self, p: [f64; 3]) -> f64 {
            match *self {
                Self::Slab { lo, hi } => (lo - p[0]).max(p[0] - hi),
                Self::Dug { lo, hi, c, r } => {
                    let slab = (lo - p[0]).max(p[0] - hi);
                    slab.max(-(dist(p, c) - r))
                }
                Self::Rim { c, r } => dist(p, c) - r,
            }
        }

        /// Analytic, not central differences.
        ///
        /// Two reasons and both matter here: the FWGSS line searches query the
        /// gradient once per Frank-Wolfe iteration, so six extra samples each
        /// would put the cost of a finite-difference stencil into C2's answer;
        /// and a differenced gradient across the `max` seam is a blend of two
        /// faces rather than either, which is the discontinuity the paper's §6
        /// warns gradient methods about.
        ///
        /// On the sphere's centre — the medial axis — the direction is `+x̂`.
        /// `M-172` measured `BrushStack::gradient` returning exactly `[0, 0, 0]`
        /// there, and a zero direction is not a direction; a fixed one is at
        /// least deterministic.
        fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
            let radial = |c: [f64; 3], sign: f64| {
                let d = dist(p, c);
                if d <= 0.0 {
                    [sign, 0.0, 0.0]
                } else {
                    [
                        sign * (p[0] - c[0]) / d,
                        sign * (p[1] - c[1]) / d,
                        sign * (p[2] - c[2]) / d,
                    ]
                }
            };
            let slab_grad = |lo: f64, hi: f64| {
                if (lo - p[0]) >= (p[0] - hi) {
                    [-1.0, 0.0, 0.0]
                } else {
                    [1.0, 0.0, 0.0]
                }
            };
            match *self {
                Self::Slab { lo, hi } => slab_grad(lo, hi),
                Self::Dug { lo, hi, c, r } => {
                    let slab = (lo - p[0]).max(p[0] - hi);
                    if slab >= -(dist(p, c) - r) {
                        slab_grad(lo, hi)
                    } else {
                        radial(c, -1.0)
                    }
                }
                Self::Rim { c, r } => {
                    let _ = r;
                    radial(c, 1.0)
                }
            }
        }
    }

    impl Scene {
        /// The exact slab, thickness `t`, centred on `x = 0`.
        fn slab(t: f64) -> Self {
            Self::Slab {
                lo: -0.5 * t,
                hi: 0.5 * t,
            }
        }

        /// `game_dig`'s 0.5-unit wall dug from behind so the residual thickness on
        /// the axis is exactly `t`, and the residual is centred on `x = 0`.
        ///
        /// The brush centre is `BRUSH + t/2`, so the brush's near face lands at
        /// `+t/2` and the residual is `[−t/2, +t/2]`. Off-axis the residual grows
        /// to `t + BRUSH − sqrt(BRUSH² − ρ²)`, which is why [`RHO_CELLS`] is half
        /// a cell.
        fn dug(t: f64) -> Self {
            Self::Dug {
                lo: -0.5 * t,
                hi: -0.5 * t + WALL,
                c: [BRUSH + 0.5 * t, 0.0, 0.0],
                r: BRUSH,
            }
        }

        /// The x-interval a discrete frame has to land in for the wall to be seen
        /// at all, widened by the element radius. Used only to bound the frame
        /// range that gets tested; testing extra frames is harmless.
        fn span(&self) -> (f64, f64) {
            match *self {
                Self::Slab { lo, hi } | Self::Dug { lo, hi, .. } => (lo, hi),
                Self::Rim { c, r } => (c[0] - r, c[0] + r),
            }
        }
    }

    // ── a deterministic stream of shots ──────────────────────────────────────

    /// SplitMix64. Seeded per row, so every row's shot set is reproducible from
    /// the row's own identity rather than from execution order.
    struct Rng(u64);

    impl Rng {
        fn new(seed: u64) -> Self {
            Self(seed)
        }

        fn bits(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }

        /// Uniform on `[0, 1)`.
        fn unit(&mut self) -> f64 {
            (self.bits() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
        }
    }

    // ── vector helpers ───────────────────────────────────────────────────────

    fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
    }

    fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }

    fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
        [a[0] * s, a[1] * s, a[2] * s]
    }

    fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    fn norm(a: [f64; 3]) -> f64 {
        dot(a, a).sqrt()
    }

    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    /// A right-handed orthonormal frame whose **second** axis is `axis`.
    ///
    /// The helper vector switches on `|a.x| < 0.9` so the cross product is never
    /// near-degenerate; the resulting azimuth of `u` about `axis` is arbitrary,
    /// which is fine — C3 randomises the capsule's axis, not its roll, because a
    /// capsule is a surface of revolution and its roll is not observable.
    fn frame(axis: [f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let a = scale(axis, 1.0 / norm(axis));
        let helper = if a[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 0.0, 1.0]
        };
        let u = cross(a, helper);
        let u = scale(u, 1.0 / norm(u));
        let w = cross(a, u);
        (u, a, scale(w, 1.0 / norm(w)))
    }

    // ── the sphere path: a sound Lipschitz branch and bound ─────────────────

    /// First time of impact of a sphere swept linearly from `p0` by `step`, over
    /// `τ ∈ [0, 1]`.
    ///
    /// **This started as plain conservative advancement (Mirtich 1996; the
    /// temporal half of Macklin, Erleben, Müller, Chentanez, Jeschke and Corse,
    /// `10.1145/3384538`), and the first run of this harness stalled it.** With
    /// `L = 1`, a clearance `c` licenses advancing `τ` by `c / |step|`, and the
    /// clearance then shrinks by the factor `1 − cos θ` where `θ` is the angle
    /// between the motion and the surface normal. On the `slab` face that factor
    /// is at most `0.293` and twenty iterations reach `f64`. On the **`dug`**
    /// field it is not: the brush leaves a concave seam where the pocket meets the
    /// slab, and a sphere passing that seam has a motion nearly *tangential* to
    /// the pocket's own surface, so `cos θ → 0`, the factor `→ 1`, and the march
    /// never converges. The harness hit its 512-iteration assertion at
    /// `τ = 0.1824`. That is the paper's own §6 limitation — "our ray-marching
    /// approach is cheap and avoids missing contacts, but requires a wide
    /// padding" — arriving as a stall rather than as a missed contact.
    ///
    /// So the march is replaced by a branch and bound over time intervals, which
    /// resolves a graze in `O(log)` rather than never. For a 1-Lipschitz `φ` and a
    /// straight segment of length `L`, over `[a, b]` with midpoint `m`:
    ///
    /// ```text
    /// min φ(p(t)) ≥ φ(p(m)) − L·(b − a)/2.
    /// t ∈ [a,b]
    /// ```
    ///
    /// An interval whose bound exceeds the sphere's radius **provably** contains
    /// no contact and is dropped; otherwise it is halved and the *earlier* half is
    /// examined first, so the first leaf reached is the earliest contact. Pruning
    /// is sound, so no contact is ever skipped — which is the property C1's CCD
    /// arm is claiming, and it is a property of the algorithm rather than of the
    /// fixture.
    ///
    /// The leaf resolution is `L · 2^-MAX_DEPTH`, which at the fastest speed in
    /// the sweep is `2.3e-12` world units — six orders below `CONTACT_TOL` — so a
    /// leaf that could not be pruned is a genuine approach to within a distance no
    /// `f64` position can distinguish from contact.
    ///
    /// # Panics
    ///
    /// If [`MAX_NODES`] is exhausted. There is deliberately no "give up and say no
    /// collision" branch: that is the tunnelling this experiment measures,
    /// arriving through the instrument instead of through the geometry.
    fn toi_sphere(
        scene: &Scene,
        p0: [f64; 3],
        step: [f64; 3],
        radius: f64,
        nodes_high_water: &mut u32,
    ) -> Option<f64> {
        let len = norm(step);
        assert!(len > 0.0, "a zero step has no time of impact to find");
        if scene.sample(p0) - radius <= CONTACT_TOL {
            *nodes_high_water = (*nodes_high_water).max(1);
            return Some(0.0);
        }
        // Depth-first, earliest interval first: push the later half, then the
        // earlier, so the earlier one pops next.
        let mut stack: Vec<(f64, f64, u32)> = Vec::with_capacity(2 * MAX_DEPTH as usize);
        stack.push((0.0, 1.0, 0));
        let mut nodes = 0u32;
        while let Some((a, b, depth)) = stack.pop() {
            nodes += 1;
            assert!(
                nodes <= MAX_NODES,
                "the Lipschitz branch and bound exhausted {MAX_NODES} nodes on [{a}, {b}] at \
                 depth {depth}: the trajectory grazes the isosurface over a region the bound \
                 cannot separate, and rounding that to `no collision` is the failure this \
                 experiment exists to measure"
            );
            let mid = 0.5 * (a + b);
            let clearance = scene.sample(add(p0, scale(step, mid))) - radius;
            if clearance > len * (b - a) * 0.5 {
                continue;
            }
            if depth >= MAX_DEPTH {
                *nodes_high_water = (*nodes_high_water).max(nodes);
                return Some(a);
            }
            stack.push((mid, b, depth + 1));
            stack.push((a, mid, depth + 1));
        }
        *nodes_high_water = (*nodes_high_water).max(nodes);
        None
    }

    // ── the triangle path: the paper's Algorithm 2 ───────────────────────────

    /// Iteration budgets for FWGSS, and the label the CSV carries them under.
    #[derive(Clone, Copy)]
    struct Solver {
        /// Algorithm 2's `max iterations`.
        fw: u32,
        /// Algorithm 1's iteration count, fixed instead of tolerance-driven.
        gss: u32,
        /// The `solver` column's value, so a row names the budget it was run at.
        name: &'static str,
    }

    /// Algorithm 1, golden-section search, with a fixed iteration count.
    ///
    /// Four maintained points at `α = 0, 1 − 1/φ, 1/φ, 1`, comparing
    /// `min(f0, f1)` against `min(f2, f3)` and keeping the two intervals on the
    /// winning side, exactly as the paper's Algorithm 1. The tolerance test is
    /// replaced by a count so that cost is data-independent; the price of that is
    /// measured as `solver_error_cells`.
    fn gss(iters: u32, mut f: impl FnMut(f64) -> f64) -> f64 {
        const INV_PHI: f64 = 0.618_033_988_749_894_9;
        let (mut a0, mut a3) = (0.0_f64, 1.0_f64);
        let mut a1 = a3 - (a3 - a0) * INV_PHI;
        let mut a2 = a0 + (a3 - a0) * INV_PHI;
        let (mut f0, mut f1) = (f(a0), f(a1));
        let (mut f2, mut f3) = (f(a2), f(a3));
        for _ in 0..iters {
            if f0.min(f1) < f2.min(f3) {
                a3 = a2;
                f3 = f2;
                a2 = a1;
                f2 = f1;
                a1 = a3 - (a3 - a0) * INV_PHI;
                f1 = f(a1);
            } else {
                a0 = a1;
                f0 = f1;
                a1 = a2;
                f1 = f2;
                a2 = a0 + (a3 - a0) * INV_PHI;
                f2 = f(a2);
            }
        }
        let mut best = (f0, a0);
        for (v, a) in [(f1, a1), (f2, a2), (f3, a3)] {
            if v < best.0 {
                best = (v, a);
            }
        }
        best.1
    }

    /// Algorithm 2, FWGSS: first time of impact between a linearly translating
    /// triangle and the field, over `[t_start, t_end] = [0, hi]`.
    ///
    /// Structure, line for line with the paper:
    ///
    /// - the starting iterate is §4.4's vertex, the one minimising `v · ∇φ(s)`
    ///   rather than Macklin et al.'s closest vertex;
    /// - the temporal sub-problem golden-sections the *unsigned* distance
    ///   backward over `[t_start, t_i]` when `φ ≤ 0`, and the signed distance
    ///   over `[t_start, t_i]` or `[t_i, t_end]` according to Equation 17's
    ///   `d_i = sign(∇φ · v)` otherwise;
    /// - every `φ ≤ 0` lowers `t_end` (Equation 9's constraint as an upper bound,
    ///   which is what culls later temporal minima);
    /// - the spatial sub-problem golden-sections `φ` along the segment from the
    ///   current point to Equation 15's support vertex.
    ///
    /// Returns the final `t_end` if any penetration was seen, `None` otherwise.
    fn toi_triangle(
        scene: &Scene,
        p0: [f64; 3],
        step: [f64; 3],
        tri: [[f64; 3]; 3],
        solver: Solver,
        hi: f64,
    ) -> Option<f64> {
        let at = |bary: [f64; 3], t: f64| {
            let q = [
                bary[0] * tri[0][0] + bary[1] * tri[1][0] + bary[2] * tri[2][0],
                bary[0] * tri[0][1] + bary[1] * tri[1][1] + bary[2] * tri[2][1],
                bary[0] * tri[0][2] + bary[1] * tri[1][2] + bary[2] * tri[2][2],
            ];
            add(add(p0, scale(step, t)), q)
        };

        // §4.4's starting iterate: the vertex most likely to collide is the one
        // minimising v · ∇φ(s), not the one nearest the surface. The paper found
        // Macklin et al.'s nearest-vertex heuristic *hindered* convergence,
        // because the vertex closest to the SDF at the start of the interval is
        // not necessarily the one moving toward it.
        let mut bary = {
            let mut best = (f64::INFINITY, 0usize);
            for (i, v) in tri.iter().enumerate() {
                let score = dot(step, scene.gradient(add(p0, *v)));
                if score < best.0 {
                    best = (score, i);
                }
            }
            let mut b = [0.0, 0.0, 0.0];
            b[best.1] = 1.0;
            b
        };

        let mut t_end = hi;
        let mut t_i = 0.0_f64;
        // The earliest time at which contact has been *established*.
        //
        // **This is `t_i`'s crossing, not `t_end`, and the first run of this
        // harness read the wrong one.** `t_end` is Equation 9's constraint used as
        // a culling device: it is lowered only where a *penetrating* sample was
        // seen, and the unsigned-distance line search deliberately lands ON the
        // isosurface, where `φ` is `±ε`. When it landed on the `+ε` side, `t_end`
        // was never lowered below the penetrating time and the harness reported a
        // ToI of exactly `1.0` against a closed-form `0.5` — a 0.576-cell error,
        // systematic per orientation because the sign of that `ε` is deterministic.
        // The paper's own §4.3 describes the same behaviour from the other side:
        // "the candidate contact point will often *wiggle* around a solution,
        // getting continually closer, but may never quite reach its target". Its
        // Table 1 measures ΔTOI against ground truth, so the reported figure is
        // the converged crossing; `t_end` is the bound that culls later minima.
        let mut toi: Option<f64> = None;

        for _ in 0..solver.fw {
            // ── temporal sub-problem ─────────────────────────────────────────
            let x = at(bary, t_i);
            let phi = scene.sample(x);
            let grad = scene.gradient(x);
            let t_next = if phi <= 0.0 {
                // Equation 9 violated: a new temporal upper bound, and a valid if
                // late upper bound on the ToI itself.
                t_end = t_end.min(t_i);
                toi = Some(toi.map_or(t_i, |b: f64| b.min(t_i)));
                let span = t_i;
                let a = gss(solver.gss, |g| scene.sample(at(bary, span * g)).abs());
                let cross = span * a;
                // The crossing is a contact time: `φ` penetrates at `t_i ≥ cross`
                // and the search minimised `|φ|` on the way there.
                toi = Some(toi.map_or(cross, |b: f64| b.min(cross)));
                cross
            } else {
                // Equation 17. `v_t` is the constant `step` because the motion is
                // a translation; the paper's rigid-body form projects ω as well.
                // `d_i = −sign(∇φ · v)` is the descent direction in time, so a
                // *receding* point searches backward and an approaching one
                // forward.
                let (lo, span) = if dot(grad, step) > 0.0 {
                    (0.0, t_i)
                } else {
                    (t_i, t_end - t_i)
                };
                let a = gss(solver.gss, |g| scene.sample(at(bary, lo + span * g)));
                lo + span * a
            };

            // ── spatial sub-problem ──────────────────────────────────────────
            let x = at(bary, t_next);
            let phi = scene.sample(x);
            if phi <= 0.0 {
                t_end = t_end.min(t_next);
                toi = Some(toi.map_or(t_next, |b: f64| b.min(t_next)));
            }
            let grad = scene.gradient(x);
            // Equation 15: the support vertex of the triangle against ∇φ.
            let mut support = (f64::INFINITY, 0usize);
            for (i, v) in tri.iter().enumerate() {
                let s = dot(grad, add(add(p0, scale(step, t_next)), *v));
                if s < support.0 {
                    support = (s, i);
                }
            }
            let mut target = [0.0, 0.0, 0.0];
            target[support.1] = 1.0;
            let towards = |g: f64| {
                [
                    bary[0] + (target[0] - bary[0]) * g,
                    bary[1] + (target[1] - bary[1]) * g,
                    bary[2] + (target[2] - bary[2]) * g,
                ]
            };
            let a = gss(solver.gss, |g| scene.sample(at(towards(g), t_next)));
            bary = towards(a);
            t_i = t_next;

            if scene.sample(at(bary, t_i)) <= 0.0 {
                t_end = t_end.min(t_i);
                toi = Some(toi.map_or(t_i, |b: f64| b.min(t_i)));
            }
        }

        toi.map(|t| t.clamp(0.0, hi))
    }

    /// §4.6's broad phase, and the reason the paper's 888 became 36.
    ///
    /// The triangle's bounding sphere, marched against the field by
    /// [`toi_sphere`]. Because `φ` is 1-Lipschitz, `min over the triangle of φ ≥
    /// φ(centroid) − R`, so a bounding sphere that never reaches the isosurface
    /// proves the triangle never does — the cull is **sound**, not heuristic.
    /// The paper needs a padding term for rotation (their Figure 4); these
    /// trajectories are pure translations, so the padding is zero and this cull is
    /// exact. That flatters the cull rate relative to their rigid-body scenes and
    /// is said out loud rather than left implicit.
    fn broad_phase(
        scene: &Scene,
        p0: [f64; 3],
        step: [f64; 3],
        tri: [[f64; 3]; 3],
        iters: &mut u32,
    ) -> Option<f64> {
        let g = scale(add(add(tri[0], tri[1]), tri[2]), 1.0 / 3.0);
        let r = tri
            .iter()
            .map(|v| norm(sub(*v, g)))
            .fold(0.0_f64, f64::max);
        toi_sphere(scene, add(p0, g), step, r, iters)
    }

    /// The whole narrow-phase query for a triangle body: broad phase, then FWGSS
    /// on the survivors.
    ///
    /// `shrink` is the paper's own inter-element optimisation — "the time interval
    /// for CCD shrinks whenever an earlier contact is found ... this makes
    /// subsequent triangle tests more efficient" — and is what makes their Table 3
    /// sub-linear in triangle count. It is a switch here because C2's clause is
    /// about cost *per element*, and an optimisation that makes element `n + 1`
    /// cheaper than element `n` is exactly the thing a linearity fit would report
    /// as a failure while being an improvement.
    fn query_triangles(
        scene: &Scene,
        p0: [f64; 3],
        step: [f64; 3],
        tris: &[[[f64; 3]; 3]],
        solver: Solver,
        shrink: bool,
    ) -> (Option<f64>, usize, u32) {
        let mut best: Option<f64> = None;
        let mut tested = 0usize;
        let mut iters = 0u32;
        for tri in tris {
            let hi = if shrink { best.unwrap_or(1.0) } else { 1.0 };
            if hi <= 0.0 {
                continue;
            }
            // The broad phase is asked about the full interval regardless of
            // `shrink`: a cull is a cull, and narrowing its interval too would
            // conflate the two effects.
            if broad_phase(scene, p0, step, *tri, &mut iters).is_none() {
                continue;
            }
            tested += 1;
            if let Some(t) = toi_triangle(scene, p0, step, *tri, solver, hi) {
                best = Some(best.map_or(t, |b: f64| b.min(t)));
            }
        }
        (best, tested, iters)
    }

    /// The same query for a body of swept spheres.
    fn query_spheres(
        scene: &Scene,
        p0: [f64; 3],
        step: [f64; 3],
        beads: &[([f64; 3], f64)],
        shrink: bool,
    ) -> (Option<f64>, usize, u32) {
        let mut best: Option<f64> = None;
        let mut iters = 0u32;
        for (c, r) in beads {
            if let Some(t) = toi_sphere(scene, add(p0, *c), step, *r, &mut iters) {
                if shrink && best.is_some_and(|b: f64| b <= t) {
                    continue;
                }
                best = Some(best.map_or(t, |b: f64| b.min(t)));
            }
        }
        (best, beads.len(), iters)
    }

    // ── capsule representations ──────────────────────────────────────────────

    /// Eight sphere centres on the capsule's segment, endpoints included.
    ///
    /// Endpoints included is the whole mechanism of C3's null: the support
    /// function of a union of spheres is the max of theirs, and with `±h·â` among
    /// the centres that max is `|h (â·n̂)| + r`, which is the capsule's support
    /// **exactly, in every direction**. So against a plane this proxy has no
    /// error at all, and the ToI disagreement C3 measures there is entirely the
    /// triangle mesh's.
    fn capsule_beads(radius: f64, half_height: f64, axis: [f64; 3], n: usize) -> Vec<([f64; 3], f64)> {
        assert!(n >= 2, "a bead chain needs both endpoints");
        (0..n)
            .map(|k| {
                let s = -half_height + 2.0 * half_height * (k as f64) / ((n - 1) as f64);
                (scale(axis, s), radius)
            })
            .collect()
    }

    /// A capsule mesh with exactly `2 · longitudes · rings` triangles.
    ///
    /// Rings are placed at **uniform arc length along the profile** — quarter
    /// circle, cylinder, quarter circle — rather than at uniform polar angle,
    /// which would collapse the cylindrical band to a single quad however long the
    /// capsule is. `(10, 10)` is the registered 200; `(12, 37)` is the paper's 888.
    fn capsule_mesh(
        radius: f64,
        half_height: f64,
        axis: [f64; 3],
        longitudes: usize,
        rings: usize,
    ) -> Vec<[[f64; 3]; 3]> {
        let (u, a, w) = frame(axis);
        let cap = std::f64::consts::FRAC_PI_2 * radius;
        let profile = 2.0 * cap + 2.0 * half_height;
        let place = |s: f64| -> (f64, f64) {
            if s <= cap {
                let ang = s / radius;
                (-half_height - radius * ang.cos(), radius * ang.sin())
            } else if s <= cap + 2.0 * half_height {
                (-half_height + (s - cap), radius)
            } else {
                let ang = (s - cap - 2.0 * half_height) / radius;
                (half_height + radius * ang.sin(), radius * ang.cos())
            }
        };
        let point = |axial: f64, radial: f64, theta: f64| {
            add(
                scale(a, axial),
                add(scale(u, radial * theta.cos()), scale(w, radial * theta.sin())),
            )
        };
        let ring: Vec<Vec<[f64; 3]>> = (0..rings)
            .map(|i| {
                let s = profile * ((i + 1) as f64) / ((rings + 1) as f64);
                let (axial, radial) = place(s);
                (0..longitudes)
                    .map(|j| {
                        let theta =
                            std::f64::consts::TAU * (j as f64) / (longitudes as f64);
                        point(axial, radial, theta)
                    })
                    .collect()
            })
            .collect();
        let bottom = scale(a, -half_height - radius);
        let top = scale(a, half_height + radius);
        let mut tris = Vec::with_capacity(2 * longitudes * rings);
        for j in 0..longitudes {
            let k = (j + 1) % longitudes;
            tris.push([bottom, ring[0][j], ring[0][k]]);
            tris.push([top, ring[rings - 1][k], ring[rings - 1][j]]);
        }
        for i in 0..rings.saturating_sub(1) {
            for j in 0..longitudes {
                let k = (j + 1) % longitudes;
                tris.push([ring[i][j], ring[i + 1][j], ring[i + 1][k]]);
                tris.push([ring[i][j], ring[i + 1][k], ring[i][k]]);
            }
        }
        assert_eq!(
            tris.len(),
            2 * longitudes * rings,
            "the capsule tessellation must hit its registered triangle count exactly"
        );
        tris
    }

    /// Support of a triangle soup in direction `+x̂`. The capsule mesh is convex,
    /// so its support is attained at a vertex.
    fn support_x_mesh(tris: &[[[f64; 3]; 3]]) -> f64 {
        tris.iter()
            .flat_map(|t| t.iter())
            .map(|v| v[0])
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Support of a bead chain in direction `+x̂`.
    fn support_x_beads(beads: &[([f64; 3], f64)]) -> f64 {
        beads
            .iter()
            .map(|(c, r)| c[0] + r)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Exact distance from `q` to a triangle. Ericson, *Real-Time Collision
    /// Detection*, §5.1.5: the closest point is in the interior, on an edge, or at
    /// a vertex, and the barycentric region tests decide which without iteration.
    fn tri_point_distance(tri: [[f64; 3]; 3], q: [f64; 3]) -> f64 {
        let ab = sub(tri[1], tri[0]);
        let ac = sub(tri[2], tri[0]);
        let aq = sub(q, tri[0]);
        let d1 = dot(ab, aq);
        let d2 = dot(ac, aq);
        if d1 <= 0.0 && d2 <= 0.0 {
            return norm(aq);
        }
        let bq = sub(q, tri[1]);
        let d3 = dot(ab, bq);
        let d4 = dot(ac, bq);
        if d3 >= 0.0 && d4 <= d3 {
            return norm(bq);
        }
        let vc = d1 * d4 - d3 * d2;
        if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
            let v = d1 / (d1 - d3);
            return norm(sub(q, add(tri[0], scale(ab, v))));
        }
        let cq = sub(q, tri[2]);
        let d5 = dot(ab, cq);
        let d6 = dot(ac, cq);
        if d6 >= 0.0 && d5 <= d6 {
            return norm(cq);
        }
        let vb = d5 * d2 - d1 * d6;
        if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
            let w = d2 / (d2 - d6);
            return norm(sub(q, add(tri[0], scale(ac, w))));
        }
        let va = d3 * d6 - d5 * d4;
        if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
            let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return norm(sub(q, add(tri[1], scale(sub(tri[2], tri[1]), w))));
        }
        let denom = 1.0 / (va + vb + vc);
        let v = vb * denom;
        let w = vc * denom;
        norm(sub(q, add(tri[0], add(scale(ab, v), scale(ac, w)))))
    }

    /// The **exact** minimum of `φ` over the whole triangle soup at time `τ`.
    ///
    /// Both fixtures admit a closed form, which is what makes an independent
    /// reference possible at all:
    ///
    /// - on a slab `φ` depends only on `x`, the mesh's `x` range is
    ///   `[inf, sup]` over its vertices, every intermediate `x` is attained
    ///   because the surface is connected, so the minimum is the one-dimensional
    ///   `min over x ∈ [inf, sup] of max(lo − x, x − hi)`;
    /// - against a sphere `φ = |p − c| − r`, so the minimum is the exact
    ///   point-to-triangle-soup distance minus `r`.
    fn mesh_min_phi(scene: &Scene, origin: [f64; 3], tris: &[[[f64; 3]; 3]]) -> f64 {
        match *scene {
            Scene::Slab { lo, hi } => {
                let (mut inf, mut sup) = (f64::INFINITY, f64::NEG_INFINITY);
                for t in tris {
                    for v in t {
                        let x = origin[0] + v[0];
                        inf = inf.min(x);
                        sup = sup.max(x);
                    }
                }
                let x = (0.5 * (lo + hi)).clamp(inf, sup);
                (lo - x).max(x - hi)
            }
            Scene::Rim { c, r } => {
                let mut best = f64::INFINITY;
                for t in tris {
                    let moved = [add(origin, t[0]), add(origin, t[1]), add(origin, t[2])];
                    best = best.min(tri_point_distance(moved, c));
                }
                best - r
            }
            // `Dug` is C1's field and no C3 fixture uses it; a reference for a
            // `max(a, −b)` composition would need the composition's own closest
            // point, which is not a closed form. Reaching here is a wiring error
            // rather than a limitation, so it says so.
            Scene::Dug { .. } => unreachable!("C3's reference is defined for `Slab` and `Rim`"),
        }
    }

    /// The exact first time of impact of the whole mesh, by the same sound
    /// Lipschitz branch and bound [`toi_sphere`] uses — but with the **exact**
    /// per-time minimum over the mesh instead of one sphere's clearance.
    ///
    /// This is the reference `toi_capsule_mesh` is scored against, and it exists
    /// because the first run could not tell a **solver false negative** from a mesh
    /// that genuinely does not reach the wall inside one step. The paper measures
    /// the same distinction against a brute-force ground truth and reports 0.94%
    /// false negatives and 0.98% false positives over 2048 configurations, so the
    /// distinction is not hypothetical.
    fn mesh_reference_toi(
        scene: &Scene,
        p0: [f64; 3],
        step: [f64; 3],
        tris: &[[[f64; 3]; 3]],
    ) -> Option<f64> {
        let len = norm(step);
        if mesh_min_phi(scene, p0, tris) <= CONTACT_TOL {
            return Some(0.0);
        }
        let mut stack: Vec<(f64, f64, u32)> = vec![(0.0, 1.0, 0)];
        while let Some((a, b, depth)) = stack.pop() {
            let mid = 0.5 * (a + b);
            let clearance = mesh_min_phi(scene, add(p0, scale(step, mid)), tris);
            if clearance > len * (b - a) * 0.5 {
                continue;
            }
            if depth >= MAX_DEPTH {
                return Some(a);
            }
            stack.push((mid, b, depth + 1));
            stack.push((a, mid, depth + 1));
        }
        None
    }

    // ── timing ───────────────────────────────────────────────────────────────

    /// Median nanoseconds and median cycles over [`REPS`] runs of `body`.
    ///
    /// Cycles because of `M-280`: the same binary reported Marching Cubes at 48³
    /// as 8.13 and 14.66 ns/sample with cycles/sample unchanged at ~34, and
    /// nothing on the face of either number said which clock it was. Every timed
    /// row here carries `ghz = cycles / nanoseconds` so the artefact states the
    /// clock rather than inviting the inference.
    fn timed<T>(probe: &mut Probe, mut body: impl FnMut() -> T) -> (f64, f64) {
        black_box(body());
        let mut runs: Vec<(f64, f64)> = Vec::with_capacity(REPS);
        for _ in 0..REPS {
            probe.reset_and_enable();
            let start = Instant::now();
            black_box(body());
            let nanos = start.elapsed().as_nanos() as f64;
            probe.disable();
            let counts = probe.read();
            assert!(
                counts.worst_ratio() >= MIN_TIME_RATIO,
                "a counter was multiplexed ({:.4}); its value is an extrapolation, not a measurement",
                counts.worst_ratio()
            );
            runs.push((nanos, counts.cycles.count as f64));
        }
        runs.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        runs[REPS / 2]
    }

    /// Least squares `y = a + b·n`, returning `(b, r2)`.
    fn fit(points: &[(f64, f64)]) -> (f64, f64) {
        let n = points.len() as f64;
        let mean_x = points.iter().map(|p| p.0).sum::<f64>() / n;
        let mean_y = points.iter().map(|p| p.1).sum::<f64>() / n;
        let sxx: f64 = points.iter().map(|p| (p.0 - mean_x).powi(2)).sum();
        let sxy: f64 = points
            .iter()
            .map(|p| (p.0 - mean_x) * (p.1 - mean_y))
            .sum();
        let slope = sxy / sxx;
        let intercept = mean_y - slope * mean_x;
        let ss_tot: f64 = points.iter().map(|p| (p.1 - mean_y).powi(2)).sum();
        let ss_res: f64 = points
            .iter()
            .map(|p| (p.1 - (intercept + slope * p.0)).powi(2))
            .sum();
        (slope, 1.0 - ss_res / ss_tot)
    }

    fn f(v: f64) -> String {
        format!("{v:.6}")
    }

    // ── C1 ───────────────────────────────────────────────────────────────────

    struct C1Row {
        field: &'static str,
        thickness_cells: f64,
        radius_cells: f64,
        speed_name: &'static str,
        speed: f64,
        step_cells: f64,
        tunnels_discrete: u32,
        tunnels_ccd: u32,
        crossing: u32,
        mismatches: u32,
        predicted: u32,
        bnb_nodes: u32,
        worst_toi_error_cells: f64,
        /// How many times the closed form was actually compared against the
        /// solver. A zero error over zero comparisons is not a measurement.
        toi_checks: u32,
    }

    /// One `(field, thickness, radius, speed)` row: [`SHOTS`] randomised shots,
    /// each one a straight trajectory across the wall, tested by a discrete
    /// end-of-frame probe and by CCD **over the same frames**.
    fn c1_row(
        field: &'static str,
        thickness_cells: f64,
        radius_cells: f64,
        speed_name: &'static str,
        speed: f64,
        seed: u64,
    ) -> C1Row {
        let t = thickness_cells * CELL;
        let radius = radius_cells * CELL;
        let scene = if field == "slab" {
            Scene::slab(t)
        } else {
            Scene::dug(t)
        };
        let (lo, hi) = scene.span();
        let step_len = speed * DT;
        let mut rng = Rng::new(seed);
        let mut row = C1Row {
            field,
            thickness_cells,
            radius_cells,
            speed_name,
            speed,
            step_cells: step_len / CELL,
            tunnels_discrete: 0,
            tunnels_ccd: 0,
            crossing: 0,
            mismatches: 0,
            predicted: 0,
            bnb_nodes: 0,
            worst_toi_error_cells: 0.0,
            toi_checks: 0,
        };

        for _ in 0..SHOTS {
            // Direction: uniform in a cone of half-angle CONE about +x.
            let cos_t = 1.0 - rng.unit() * (1.0 - CONE.cos());
            let sin_t = (1.0 - cos_t * cos_t).max(0.0).sqrt();
            let phi = std::f64::consts::TAU * rng.unit();
            let dir = [cos_t, sin_t * phi.cos(), sin_t * phi.sin()];
            // Lateral offset of the crossing point, uniform on a disc.
            let rho = RHO_CELLS * CELL * rng.unit().sqrt();
            let psi = std::f64::consts::TAU * rng.unit();
            let cross = [0.0, rho * psi.cos(), rho * psi.sin()];
            // Sub-frame phase: where the frame boundaries fall relative to the
            // wall. This is the only quantity that decides whether a discrete
            // probe lands in the overlap window, so it is the one that must be
            // uniform.
            let u = rng.unit();
            let step = scale(dir, step_len);
            let p0 = sub(cross, scale(step, u));
            let dx = step[0];

            // Non-vacuity: the centre path is inside the solid at the crossing.
            if scene.sample(cross) < 0.0 {
                row.crossing += 1;
            }

            // Frames that could possibly see the wall. Extra frames outside the
            // window are harmless; missing one would be a false tunnel.
            let margin = radius + 2.0 * CELL;
            let k_lo = ((lo - margin - p0[0]) / dx).floor() as i64;
            let k_hi = ((hi + margin - p0[0]) / dx).ceil() as i64;

            let mut discrete = false;
            let mut ccd = false;
            let mut predicted = false;
            for k in k_lo..=k_hi {
                let p = add(p0, scale(step, k as f64));
                // The discrete path, faithful to `resolve_body`: sample the field
                // at the body's position and call it a contact when the sphere
                // overlaps. No sweep, no history.
                if scene.sample(p) < radius {
                    discrete = true;
                }
                // The closed form, for the exact slab only: |x| − t/2 < radius.
                if field == "slab" && p[0].abs() - 0.5 * t < radius {
                    predicted = true;
                }
                if let Some(toi) = toi_sphere(&scene, p, step, radius, &mut row.bnb_nodes) {
                    ccd = true;
                    if field == "slab" {
                        // Closed-form ToI for a plane face, as a check on the
                        // solver rather than on the clause: contact when the
                        // sphere's leading point reaches the near face. The
                        // *count* of comparisons is carried too — a zero error
                        // over a population that was never reached is `M-44`'s
                        // vacuous zero, and this instrument is as exposed to that
                        // as anything it is checking.
                        let exact = (lo - radius - p[0]) / dx;
                        if (0.0..=1.0).contains(&exact) {
                            let err = (toi - exact).abs() * step_len / CELL;
                            row.worst_toi_error_cells = row.worst_toi_error_cells.max(err);
                            row.toi_checks += 1;
                        }
                    }
                }
            }
            if !discrete {
                row.tunnels_discrete += 1;
            }
            if !ccd {
                row.tunnels_ccd += 1;
            }
            if !predicted {
                row.predicted += 1;
            }
            if field == "slab" && discrete != predicted {
                row.mismatches += 1;
            }
        }
        row
    }

    // ── C2 ───────────────────────────────────────────────────────────────────

    struct C2Row {
        kind: &'static str,
        solver_name: &'static str,
        shrink: bool,
        elements: usize,
        tested: usize,
        total_ns: f64,
        total_cycles: f64,
        toi: Option<f64>,
        bnb_nodes: u32,
    }

    // ── C3 ───────────────────────────────────────────────────────────────────

    struct C3Row {
        fixture: &'static str,
        thickness_cells: Option<f64>,
        half_height: f64,
        orientation: usize,
        axis: [f64; 3],
        toi_spheres: Option<f64>,
        toi_mesh: Option<f64>,
        /// The exact mesh ToI from [`mesh_reference_toi`], so a `None` in
        /// `toi_mesh` can be told apart from a mesh that never reaches the wall.
        toi_mesh_reference: Option<f64>,
        /// Reference says contact, FWGSS says none.
        mesh_false_negative: bool,
        /// FWGSS says contact, reference says none.
        mesh_false_positive: bool,
        /// A lower bound on the disagreement when one arm found nothing: the
        /// travel remaining after the other arm's contact, in cells.
        disagreement_floor_cells: Option<f64>,
        disagreement_cells: Option<f64>,
        solver_error_cells: Option<f64>,
        deficit_spheres_cells: Option<f64>,
        deficit_mesh_cells: Option<f64>,
        ns_spheres: f64,
        ns_mesh: f64,
        cycles_spheres: f64,
        cycles_mesh: f64,
        tested_mesh: usize,
        step_cells: f64,
    }

    /// Bisect the launch position so that the **true capsule** just touches the
    /// rim at `TAU_PLACE`, whatever the orientation and however long the capsule is.
    ///
    /// **The first version of this placed by the bead proxy's own analytic time of
    /// impact, and that was a fixture defect its own output caught.** At
    /// `h = 3.6` the eight beads are spaced `2h/7 = 1.029` apart, so the bead
    /// nearest the capsule's mid-point sits `h/7 = 0.514` off the axis — further
    /// than `r_rim + r_bead = 0.5`, which means **no bead can ever touch the rim**
    /// for a broadside strike. The bead ToI was therefore `+∞` at every launch
    /// position, the bisection ran its bracket up to `x₀ = 0`, and the body
    /// started centred *inside* the rim: sixteen rows reported `toi = 0.00000` on
    /// both arms, which is a comparison of two methods that were both already
    /// penetrating. `M-44` in a new place, and the reason the reference has to be
    /// the geometry both proxies approximate rather than either proxy.
    ///
    /// The true capsule's distance to a point is the point's distance to the
    /// segment, minus the radius, and that is finite for every orientation — so
    /// this placement always exists.
    fn rim_launch(
        radius: f64,
        half_height: f64,
        axis: [f64; 3],
        rim: ([f64; 3], f64),
        step_len: f64,
    ) -> f64 {
        // Surface gap between the true capsule at `TAU_PLACE` and the rim.
        let gap_at = |x0: f64| -> f64 {
            let centre = [x0 + TAU_PLACE * step_len - rim.0[0], -rim.0[1], -rim.0[2]];
            let s = (-dot(centre, axis)).clamp(-half_height, half_height);
            norm(add(centre, scale(axis, s))) - radius - rim.1
        };
        // Monotone decreasing in x0 over the approach: moving the launch closer
        // can only shrink the gap.
        let mut lo = -(rim.1 + radius + half_height + 2.0 * step_len);
        let mut hi = rim.0[0];
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if gap_at(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    pub(crate) fn run(run: &mut Run) {
        let mut probe = Probe::open();
        let mut rows: Vec<Row> = Vec::new();

        let fall = (2.0 * GRAVITY * SANDBOX_HEIGHT).sqrt();
        let speeds: [(&str, f64); 7] = [
            ("walk", WALK),
            ("sprint", SPRINT),
            ("fall_full_height", fall),
            ("sprint_plus_fall", (SPRINT * SPRINT + 2.0 * GRAVITY * SANDBOX_HEIGHT).sqrt()),
            ("projectile_40", 40.0),
            ("projectile_80", 80.0),
            ("projectile_160", 160.0),
        ];
        // The four speeds a `game_dig` body actually reaches. C1 is scored on
        // these; the three above them locate a threshold and are reported apart.
        const GAME_SPEEDS: usize = 4;

        println!("-- SHARE, recomputed before the run --");
        println!(
            "  h = {CELL}, dt = {DT} s -> a frame of travel in cells: walk {:.4}, sprint {:.4}, \
             full-height fall {:.4}, both {:.4}",
            WALK * DT / CELL,
            SPRINT * DT / CELL,
            fall * DT / CELL,
            speeds[3].1 * DT / CELL
        );
        println!(
            "  discrete tunnelling needs v*dt > 2R + t; game_dig's 2R = {:.4} cells, and the \
             fastest game_dig body moves {:.4} cells per frame.",
            2.0 * BODY_RADIUS / CELL,
            speeds[3].1 * DT / CELL
        );
        println!(
            "  => C1's discrete arm is UNREACHABLE at game_dig's own body radius, by arithmetic, \
             at every registered thickness. Running anyway."
        );
        println!(
            "  C2: 25 us/element inside 10% of a 16 ms frame buys {:.0} elements; a 200-triangle \
             capsule at the bar is {:.2}% of the whole frame.",
            0.10 * DT * 1e6 / 25.0,
            200.0 * 25.0 / (DT * 1e6) * 100.0
        );
        println!(
            "  C3: the ratio's share is 1.0 (two whole measurements), so 1/(1 - 1/10) = {:.4} and \
             the element counts alone give {:.1}x before any per-element difference.",
            1.0 / (1.0 - 0.1),
            200.0 / 8.0
        );

        // ── C1 ───────────────────────────────────────────────────────────────
        println!("\n-- C1: 10^4 shots per row, discrete against CCD over the same frames --");
        println!(
            "{:>6} {:>7} {:>7} {:>18} {:>9} {:>9} {:>9} {:>7}",
            "field", "t/h", "R/h", "speed", "step/h", "tun_dcd", "tun_ccd", "nodes"
        );
        let mut c1_rows: Vec<C1Row> = Vec::new();
        let mut seed = 0x5EED_0082_u64;
        for field in ["slab", "dug"] {
            for &radius_cells in &RADII_CELLS {
                for &(speed_name, speed) in &speeds {
                    for &thickness_cells in &THICKNESS_CELLS {
                        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                        let row =
                            c1_row(field, thickness_cells, radius_cells, speed_name, speed, seed);
                        assert_eq!(
                            row.crossing, SHOTS,
                            "{field}: {} of {SHOTS} shots did not put the element's centre path \
                             inside the solid, so a tunnel over them would be M-44's zero",
                            SHOTS - row.crossing
                        );
                        if field == "slab" {
                            assert_eq!(
                                row.mismatches, 0,
                                "slab t/h = {thickness_cells}, R/h = {radius_cells}: the discrete \
                                 arm disagreed with the closed-form predicate |x| - t/2 < R on {} \
                                 of {SHOTS} shots, so the instrument is not measuring the geometry",
                                row.mismatches
                            );
                            assert!(
                                row.toi_checks > 0,
                                "slab t/h = {thickness_cells}, R/h = {radius_cells}: the ToI \
                                 accuracy check made zero comparisons, so its zero error is \
                                 M-44's vacuous zero rather than a measurement"
                            );
                            assert!(
                                row.worst_toi_error_cells < 1e-6,
                                "slab: the branch and bound's ToI is {} cells from the closed \
                                 form over {} comparisons; the solver is wrong before the clause \
                                 is asked",
                                row.worst_toi_error_cells,
                                row.toi_checks
                            );
                        }
                        c1_rows.push(row);
                    }
                }
            }
        }
        let tunnels_discrete_total: u64 =
            c1_rows.iter().map(|r| u64::from(r.tunnels_discrete)).sum();
        let tunnels_ccd_total: u64 = c1_rows.iter().map(|r| u64::from(r.tunnels_ccd)).sum();
        for r in &c1_rows {
            // Printed selectively: 448 rows is a CSV, not a table. The two
            // endpoints of the thickness sweep are the interesting ones — two
            // cells, where a wall is still a wall, and `subgrid`'s 0.05 floor.
            let endpoint = [THICKNESS_CELLS[0], THICKNESS_CELLS[THICKNESS_CELLS.len() - 1]]
                .iter()
                .any(|t| (r.thickness_cells - t).abs() < 1e-12);
            if endpoint {
                println!(
                    "{:>6} {:>7.3} {:>7.2} {:>18} {:>9.4} {:>9} {:>9} {:>7}",
                    r.field,
                    r.thickness_cells,
                    r.radius_cells,
                    r.speed_name,
                    r.step_cells,
                    r.tunnels_discrete,
                    r.tunnels_ccd,
                    r.bnb_nodes
                );
            }
        }
        println!(
            "  VACUITY CONTROL: tunnels_discrete summed over the sweep = {tunnels_discrete_total} \
             (must be > 0); tunnels_ccd = {tunnels_ccd_total}"
        );
        // The registration's vacuity control, verbatim: "the discrete arm must
        // tunnel at least once, reported as a count, or C1 is comparing two
        // methods that both work."
        assert!(
            tunnels_discrete_total > 0,
            "VACUITY CONTROL FAILED: the discrete arm never tunnelled anywhere in the sweep, so \
             the fixture could not have distinguished the two methods"
        );

        // C1 is scored where the registration puts it: game_dig's speeds,
        // game_dig's body, and thicknesses below v*dt.
        let mut c1_scored = 0usize;
        let mut c1_holds = true;
        for r in &c1_rows {
            let game_speed = speeds[..GAME_SPEEDS].iter().any(|s| s.0 == r.speed_name);
            if game_speed
                && (r.radius_cells - BODY_RADIUS / CELL).abs() < 1e-12
                && r.thickness_cells < r.step_cells
            {
                c1_scored += 1;
                if r.tunnels_discrete == 0 || r.tunnels_ccd != 0 {
                    c1_holds = false;
                }
            }
        }
        assert!(
            c1_scored > 0,
            "C1 has no rows to score: no game_dig speed reaches a thickness in the registered sweep"
        );
        println!(
            "  C1 scored over {c1_scored} rows at game_dig's own body radius and speeds: {}",
            if c1_holds { "HELD" } else { "FALSIFIED" }
        );

        // ── C2 ───────────────────────────────────────────────────────────────
        println!("\n-- C2: cost per element tested, swept over element count --");
        // The fixture is the two-cell wall at sprint speed, and the body's extent
        // is held FIXED while the element count grows. P-72's lesson: a sweep that
        // let the body grow with the count would be measuring approach distance
        // wearing a linearity sweep's name, because conservative advancement's
        // iteration count grows with the initial clearance.
        let wall = Scene::slab(2.0 * CELL);
        let (wall_lo, _) = wall.span();
        let c2_step = [SPRINT * DT, 0.0, 0.0];
        let c2_step_len = norm(c2_step);
        let c2_axis = [0.0, 1.0, 0.0];
        let mut c2_rows: Vec<C2Row> = Vec::new();

        println!(
            "{:>16} {:>7} {:>8} {:>8} {:>11} {:>11} {:>10} {:>7}",
            "solver", "shrink", "elems", "tested", "us/tested", "us/elem", "ghz", "toi"
        );
        for shrink in [false, true] {
            for solver in SOLVERS {
                for (idx, &(longitudes, rings)) in TESSELLATIONS.iter().enumerate() {
                    let tris = capsule_mesh(BODY_RADIUS, 0.5 * BODY_SPAN, c2_axis, longitudes, rings);
                    assert_eq!(tris.len(), ELEMENT_COUNTS[idx]);
                    let sup = support_x_mesh(&tris);
                    let p0 = [wall_lo - sup - 0.5 * c2_step_len, 0.0, 0.0];
                    let (ns, cycles) = timed(&mut probe, || {
                        query_triangles(&wall, p0, c2_step, &tris, solver, shrink)
                    });
                    let (toi, tested, iters) =
                        query_triangles(&wall, p0, c2_step, &tris, solver, shrink);
                    c2_rows.push(C2Row {
                        kind: "triangle_fwgss",
                        solver_name: solver.name,
                        shrink,
                        elements: tris.len(),
                        tested,
                        total_ns: ns,
                        total_cycles: cycles,
                        toi,
                        bnb_nodes: iters,
                    });
                }
            }
            // The sphere arm: the same body, the same extent, as a bead chain.
            for &n in &ELEMENT_COUNTS {
                let beads = capsule_beads(BODY_RADIUS, 0.5 * BODY_SPAN, c2_axis, n);
                let sup = support_x_beads(&beads);
                let p0 = [wall_lo - sup - 0.5 * c2_step_len, 0.0, 0.0];
                let (ns, cycles) =
                    timed(&mut probe, || query_spheres(&wall, p0, c2_step, &beads, shrink));
                let (toi, tested, iters) = query_spheres(&wall, p0, c2_step, &beads, shrink);
                c2_rows.push(C2Row {
                    kind: "sphere_ca",
                    solver_name: "exact",
                    shrink,
                    elements: n,
                    tested,
                    total_ns: ns,
                    total_cycles: cycles,
                    toi,
                    bnb_nodes: iters,
                });
            }
        }

        // Linearity, fitted per (kind, solver, shrink) group over element count.
        let mut fits: Vec<(FitKey, (f64, f64))> = Vec::new();
        for kind in ["triangle_fwgss", "sphere_ca"] {
            for solver_name in ["fw4_gss6", "fw8_gss10", "fw16_gss20", "exact"] {
                for shrink in [false, true] {
                    let points: Vec<(f64, f64)> = c2_rows
                        .iter()
                        .filter(|r| {
                            r.kind == kind && r.solver_name == solver_name && r.shrink == shrink
                        })
                        .map(|r| (r.elements as f64, r.total_ns / 1000.0))
                        .collect();
                    if points.len() == ELEMENT_COUNTS.len() {
                        fits.push(((kind, solver_name, shrink), fit(&points)));
                    }
                }
            }
        }
        let fit_of = |kind: &str, solver_name: &str, shrink: bool| -> (f64, f64) {
            fits.iter()
                .find(|(k, _)| k.0 == kind && k.1 == solver_name && k.2 == shrink)
                .map_or((f64::NAN, f64::NAN), |(_, v)| *v)
        };
        for r in &c2_rows {
            let per_tested = if r.tested == 0 {
                f64::NAN
            } else {
                r.total_ns / 1000.0 / r.tested as f64
            };
            println!(
                "{:>16} {:>7} {:>8} {:>8} {:>11.4} {:>11.4} {:>10.4} {:>7}",
                r.solver_name,
                r.shrink,
                r.elements,
                r.tested,
                per_tested,
                r.total_ns / 1000.0 / r.elements as f64,
                r.total_cycles / r.total_ns,
                r.toi.map_or_else(|| "none".to_string(), |t| format!("{t:.4}"))
            );
        }
        for ((kind, solver_name, shrink), (slope, r2)) in &fits {
            println!(
                "  fit {kind}/{solver_name}/shrink={shrink}: marginal {slope:.4} us per element, \
                 r2 {r2:.6}"
            );
        }

        // ── C3 ───────────────────────────────────────────────────────────────
        println!("\n-- C3: 8 swept spheres against a 200-triangle capsule mesh --");
        let rim = ([0.0, 0.0, 0.0], 2.0 * CELL);
        let mut c3_rows: Vec<C3Row> = Vec::new();
        let mut orient_rng = Rng::new(0xC3_0082);
        // Three canonical orientations - end-on, broadside, and 45° in the plane
        // of motion - then uniform on the sphere. A capsule is a surface of
        // revolution, so only the axis matters and the roll does not.
        let mut axes: Vec<[f64; 3]> = vec![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [
                std::f64::consts::FRAC_1_SQRT_2,
                std::f64::consts::FRAC_1_SQRT_2,
                0.0,
            ],
        ];
        while axes.len() < ORIENTATIONS {
            let z = 2.0 * orient_rng.unit() - 1.0;
            let r = (1.0 - z * z).max(0.0).sqrt();
            let phi = std::f64::consts::TAU * orient_rng.unit();
            axes.push([r * phi.cos(), z, r * phi.sin()]);
        }

        // C3's fixtures: `(name, step in cells, wall thickness in cells or None
        // for the convex rim)`.
        //
        // `wall` is the registered fixture at `game_dig`'s own sprint step.
        // `wall_wide` is the same geometry with the range the one-cell bar needs.
        // `rim` is the convex nub — a dug pit's edge — and it is here because on a
        // *plane* the bead chain is support-exact and C3's clause could not fail;
        // see this module's header.
        let c3_fixtures: [(&str, f64, Option<f64>); 3] = [
            ("wall", SPRINT * DT / CELL, Some(2.0)),
            ("wall_wide", WIDE_STEP_CELLS, Some(WIDE_WALL_CELLS)),
            ("rim", SPRINT * DT / CELL, None),
        ];

        // The solver budget is CHOSEN BY MEASUREMENT, over exactly the population
        // whose numbers C3 then reports, against [`mesh_reference_toi`] — the
        // exact, independent answer. The rule is stated once here: **the cheapest
        // budget whose worst error over that population clears
        // `ACCURACY_BAR_CELLS` with no false negative and no false positive, and
        // the most accurate budget available if none does**, with
        // `solver_cleared_accuracy_bar` on every row saying which happened.
        // Picking the budget after seeing the cost would be tuning the instrument
        // to the answer.
        //
        // The population is every C3 fixture, not just the wall: the rim is where
        // the clause can fail, so a solver that is only accurate on planes is not
        // good enough to report it. And the error is taken against the reference
        // rather than against the closed-form *support*, because on the narrow arm
        // a mesh whose support deficit exceeds the remaining travel genuinely does
        // not touch the wall — the reference says so, and scoring that as solver
        // error is scoring the fixture.
        let solver_accuracy = |solver: Solver| -> (f64, u32, u32) {
            let mut worst = 0.0_f64;
            let (mut false_neg, mut false_pos) = (0u32, 0u32);
            for &(_, step_cells, thickness) in &c3_fixtures {
                let step = [step_cells * CELL, 0.0, 0.0];
                let step_len = norm(step);
                for &half_height in &HALF_HEIGHTS {
                    for &axis in &axes {
                        let tris = capsule_mesh(BODY_RADIUS, half_height, axis, 10, 10);
                        let (scene, p0) = if let Some(thickness) = thickness {
                            let scene = Scene::slab(thickness * CELL);
                            let (lo, _) = scene.span();
                            let sup = (half_height * axis[0]).abs() + BODY_RADIUS;
                            (scene, [lo - sup - TAU_PLACE * step_len, 0.0, 0.0])
                        } else {
                            let scene = Scene::Rim { c: rim.0, r: rim.1 };
                            let x0 = rim_launch(BODY_RADIUS, half_height, axis, rim, step_len);
                            (scene, [x0, 0.0, 0.0])
                        };
                        let want = mesh_reference_toi(&scene, p0, step, &tris);
                        let got = query_triangles(&scene, p0, step, &tris, solver, true).0;
                        match (got, want) {
                            (Some(a), Some(b)) => {
                                worst = worst.max((a - b).abs() * step_len / CELL);
                            }
                            (None, Some(_)) => false_neg += 1,
                            (Some(_), None) => false_pos += 1,
                            (None, None) => {}
                        }
                    }
                }
            }
            (worst, false_neg, false_pos)
        };
        println!(
            "  solver calibration over all {} C3 configurations, against the exact reference:",
            c3_fixtures.len() * HALF_HEIGHTS.len() * axes.len()
        );
        let mut calibration: Vec<(Solver, f64, u32, u32)> = Vec::new();
        for solver in SOLVERS {
            let (err, fneg, fpos) = solver_accuracy(solver);
            let clears = err <= ACCURACY_BAR_CELLS && fneg == 0 && fpos == 0;
            println!(
                "    {:>12}: worst {:>12.8} cells, {fneg} false negative, {fpos} false positive  {}",
                solver.name,
                err,
                if clears { "clears" } else { "over the bar" }
            );
            calibration.push((solver, err, fneg, fpos));
        }
        let cleared = calibration
            .iter()
            .find(|(_, e, n, p)| *e <= ACCURACY_BAR_CELLS && *n == 0 && *p == 0)
            .copied();
        let c3_solver = cleared.map_or_else(
            || {
                calibration
                    .iter()
                    .copied()
                    .min_by(|a, b| {
                        (a.2 + a.3)
                            .cmp(&(b.2 + b.3))
                            .then_with(|| a.1.total_cmp(&b.1))
                    })
                    .expect("SOLVERS is non-empty")
                    .0
            },
            |(s, ..)| s,
        );
        let solver_cleared_bar = cleared.is_some();
        println!(
            "  chosen: {} ({})",
            c3_solver.name,
            if solver_cleared_bar {
                "cheapest budget clearing the bar"
            } else {
                "NO budget cleared the bar; fewest misses then most accurate"
            }
        );

        println!(
            "{:>8} {:>7} {:>5} {:>10} {:>10} {:>12} {:>12} {:>10}",
            "fixture", "h", "orient", "toi_sph", "toi_mesh", "disagree/h", "cost_ratio", "err/h"
        );
        for &(fixture, step_cells, thickness_cells) in &c3_fixtures {
            for &half_height in &HALF_HEIGHTS {
                for (orientation, &axis) in axes.iter().enumerate() {
                    let beads = capsule_beads(BODY_RADIUS, half_height, axis, 8);
                    let tris = capsule_mesh(BODY_RADIUS, half_height, axis, 10, 10);
                    assert_eq!(tris.len(), 200, "C3's mesh arm is registered at 200 triangles");
                    let step = [step_cells * CELL, 0.0, 0.0];
                    let step_len = norm(step);

                    let (scene, p0, exact) = if let Some(thickness) = thickness_cells {
                        let scene = Scene::slab(thickness * CELL);
                        let (lo, _) = scene.span();
                        // Place the TRUE CAPSULE's contact at `TAU_PLACE`. Its
                        // support in `+x̂` is `|h (â·x̂)| + r`, which is also the
                        // bead chain's, exactly — which is C3's whole mechanism.
                        //
                        // `TAU_PLACE` is 0.1 and was 0.5, and that was a fixture
                        // range defect its own output caught. A representation
                        // whose support falls short of the capsule's by more than
                        // the remaining travel cannot contact the wall *at all*
                        // inside one step, so the fixture reported `none` rather
                        // than a disagreement — and at half a step of sprint that
                        // ceiling is 0.576 cells, **below C3's own one-cell bar**.
                        // A clause cannot be falsified by a fixture whose dynamic
                        // range stops short of its threshold, which is `✗51`'s rule
                        // wearing a geometric costume.
                        let sup_capsule = (half_height * axis[0]).abs() + BODY_RADIUS;
                        let p0 = [lo - sup_capsule - TAU_PLACE * step_len, 0.0, 0.0];
                        let sup_beads = support_x_beads(&beads);
                        let sup_mesh = support_x_mesh(&tris);
                        (
                            scene,
                            p0,
                            Some((
                                (lo - sup_beads - p0[0]) / step[0],
                                (lo - sup_mesh - p0[0]) / step[0],
                                (sup_capsule - sup_beads) / CELL,
                                (sup_capsule - sup_mesh) / CELL,
                            )),
                        )
                    } else {
                        let scene = Scene::Rim { c: rim.0, r: rim.1 };
                        let x0 = rim_launch(BODY_RADIUS, half_height, axis, rim, step_len);
                        (scene, [x0, 0.0, 0.0], None)
                    };

                    let (ns_s, cyc_s) =
                        timed(&mut probe, || query_spheres(&scene, p0, step, &beads, true));
                    let (ns_m, cyc_m) = timed(&mut probe, || {
                        query_triangles(&scene, p0, step, &tris, c3_solver, true)
                    });
                    let (toi_s, _, _) = query_spheres(&scene, p0, step, &beads, true);
                    let (toi_m, tested_m, _) =
                        query_triangles(&scene, p0, step, &tris, c3_solver, true);

                    // The reference: exact, independent, and the instrument that
                    // tells a solver false negative from a mesh that genuinely
                    // does not reach the wall inside one step. The bead arm needs
                    // none — `query_spheres` *is* the exact branch and bound, to
                    // `2e-12` world units — so a `None` there is a geometric fact
                    // about the proxy, which is what C3's falsifier is about.
                    let reference = mesh_reference_toi(&scene, p0, step, &tris);
                    let disagreement = match (toi_s, toi_m) {
                        (Some(a), Some(b)) => Some((a - b).abs() * step_len / CELL),
                        _ => None,
                    };
                    // When one arm found nothing, the disagreement is at least the
                    // travel left after the other arm's contact. Reported, so a
                    // miss is a number rather than a blank.
                    let floor = match (toi_s, toi_m) {
                        (Some(a), None) => Some((1.0 - a) * step_len / CELL),
                        (None, Some(b)) => Some((1.0 - b) * step_len / CELL),
                        _ => None,
                    };
                    let solver_error = match (toi_m, reference) {
                        (Some(got), Some(want)) => Some((got - want).abs() * step_len / CELL),
                        _ => None,
                    };
                    c3_rows.push(C3Row {
                        thickness_cells,
                        fixture,
                        half_height,
                        orientation,
                        axis,
                        toi_spheres: toi_s,
                        toi_mesh: toi_m,
                        toi_mesh_reference: reference,
                        mesh_false_negative: toi_m.is_none() && reference.is_some(),
                        mesh_false_positive: toi_m.is_some() && reference.is_none(),
                        disagreement_floor_cells: floor,
                        disagreement_cells: disagreement,
                        solver_error_cells: solver_error,
                        deficit_spheres_cells: exact.map(|e| e.2),
                        deficit_mesh_cells: exact.map(|e| e.3),
                        ns_spheres: ns_s,
                        ns_mesh: ns_m,
                        cycles_spheres: cyc_s,
                        cycles_mesh: cyc_m,
                        tested_mesh: tested_m,
                        step_cells: step_len / CELL,
                    });
                }
            }
        }
        // The mechanism, asserted rather than argued: a bead chain whose endpoints
        // are the segment's endpoints has the capsule's support function exactly,
        // so its deficit against a plane is zero to floating point.
        let worst_bead_deficit = c3_rows
            .iter()
            .filter_map(|r| r.deficit_spheres_cells)
            .fold(0.0_f64, |a, b| a.max(b.abs()));
        assert!(
            worst_bead_deficit < 1e-12,
            "the bead chain's support deficit against a plane is {worst_bead_deficit} cells, not \
             zero - then the analytic mechanism behind C3's wall arm is wrong"
        );
        let worst_solver_error = c3_rows
            .iter()
            .filter_map(|r| r.solver_error_cells)
            .fold(0.0_f64, f64::max);
        let mesh_false_negatives = c3_rows.iter().filter(|r| r.mesh_false_negative).count();
        let mesh_false_positives = c3_rows.iter().filter(|r| r.mesh_false_positive).count();
        let proxy_true_misses = c3_rows
            .iter()
            .filter(|r| r.toi_spheres.is_none() && r.toi_mesh_reference.is_some())
            .count();
        println!(
            "  bead support deficit against a plane: {worst_bead_deficit:.3e} cells (exactly zero \
             by construction). Worst FWGSS ToI error against the exact reference: \
             {worst_solver_error:.6} cells, bar {ACCURACY_BAR_CELLS}. FWGSS misses: \
             {mesh_false_negatives} false negative, {mesh_false_positives} false positive of {}. \
             8-bead proxy misses a contact the mesh really has: {proxy_true_misses}",
            c3_rows.len()
        );
        for r in &c3_rows {
            if r.orientation < 3 {
                println!(
                    "{:>8} {:>7.2} {:>5} {:>10} {:>10} {:>12} {:>12} {:>10}",
                    r.fixture,
                    r.half_height,
                    r.orientation,
                    r.toi_spheres
                        .map_or_else(|| "none".to_string(), |t| format!("{t:.5}")),
                    r.toi_mesh
                        .map_or_else(|| "none".to_string(), |t| format!("{t:.5}")),
                    r.disagreement_cells
                        .map_or_else(|| "MISS".to_string(), |d| format!("{d:.5}")),
                    format!("{:.2}", r.ns_mesh / r.ns_spheres),
                    r.solver_error_cells
                        .map_or_else(|| "-".to_string(), |e| format!("{e:.6}"))
                );
            }
        }

        // ── verdicts ─────────────────────────────────────────────────────────
        // C2 is scored on the paper's own unit - cost per element TESTED - at the
        // solver budget the calibration above chose by measurement, with the
        // paper's inter-element interval shrink OFF, because C2's clause is about
        // cost per element and the shrink makes later elements cheaper than
        // earlier ones.
        let chosen_name = c3_solver.name;
        let (c2_marginal, c2_r2) = fit_of("triangle_fwgss", chosen_name, false);
        let c2_worst_per_tested = c2_rows
            .iter()
            .filter(|r| r.kind == "triangle_fwgss" && r.solver_name == chosen_name && !r.shrink)
            .map(|r| r.total_ns / 1000.0 / r.tested.max(1) as f64)
            .fold(0.0_f64, f64::max);
        let c2_holds = solver_cleared_bar
            && c2_marginal < 25.0
            && c2_worst_per_tested < 25.0
            && c2_r2 >= 0.99;

        // C3 is scored where it is registered — a **wall**, at `game_dig`'s own
        // body — over both wall arms, because a clause is only tested by the arm
        // whose range can exceed its bar. The rim arm and the aspect sweep are
        // reported beside it and are where the clause can actually fail.
        let c3_wall_game: Vec<&C3Row> = c3_rows
            .iter()
            .filter(|r| {
                r.thickness_cells.is_some() && (r.half_height - 0.5 * BODY_SPAN).abs() < 1e-12
            })
            .collect();
        let c3_worst = c3_wall_game
            .iter()
            .map(|r| r.disagreement_cells.unwrap_or(f64::INFINITY))
            .fold(0.0_f64, f64::max);
        let c3_ratio = c3_wall_game
            .iter()
            .map(|r| r.ns_mesh / r.ns_spheres)
            .fold(f64::INFINITY, f64::min);
        let c3_holds = c3_worst <= 1.0 && c3_ratio >= 10.0;
        // The aspect ratio at which eight beads stop representing the body, taken
        // from the rim arm where a scallop trough is a real gap rather than an
        // invisible one. Reported, not scored: C3 registers `game_dig`'s capsule.
        let c3_rim_breaks: Option<f64> = c3_rows
            .iter()
            .filter(|r| {
                r.thickness_cells.is_none()
                    && (r.toi_spheres.is_none()
                        || r.disagreement_cells.is_some_and(|d| d > 1.0))
            })
            .map(|r| r.half_height)
            .fold(None, |acc: Option<f64>, h| {
                Some(acc.map_or(h, |a: f64| a.min(h)))
            });

        println!("\n-- verdicts --");
        println!("  C1 {}", if c1_holds { "HELD" } else { "FALSIFIED" });
        println!(
            "  C2 {} (solver {chosen_name}, marginal {c2_marginal:.4} us/element, worst \
             {c2_worst_per_tested:.4} us per element tested, r2 {c2_r2:.6})",
            if c2_holds { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "  C3 {} (worst disagreement {c3_worst:.6} cells over both wall arms at game_dig's \
             own capsule, worst-case cost ratio {c3_ratio:.2}x)",
            if c3_holds { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "  C3, reported not scored: the smallest capsule half-height at which the 8-bead \
             proxy misses the rim or disagrees by more than a cell is {}",
            c3_rim_breaks.map_or_else(
                || "never in the swept range".to_string(),
                |h| format!(
                    "{h} units = {:.2} cells (aspect {:.2})",
                    h / CELL,
                    (h + BODY_RADIUS) / BODY_RADIUS
                )
            )
        );

        // ── rows ─────────────────────────────────────────────────────────────
        for r in &c1_rows {
            rows.push(vec![
                ("arm", "c1_tunnelling".to_string()),
                ("field", r.field.to_string()),
                ("wall_thickness_cells", f(r.thickness_cells)),
                ("body_radius_cells", f(r.radius_cells)),
                ("speed_name", r.speed_name.to_string()),
                ("speed_units_per_second", f(r.speed)),
                ("speed_cells_per_step", f(r.step_cells)),
                ("frame_seconds", f(DT)),
                ("cell_size", f(CELL)),
                ("shots", SHOTS.to_string()),
                ("tunnels_discrete", r.tunnels_discrete.to_string()),
                ("tunnels_ccd", r.tunnels_ccd.to_string()),
                ("tunnels_discrete_predicted", r.predicted.to_string()),
                ("prediction_mismatches", r.mismatches.to_string()),
                ("shots_crossing_solid", r.crossing.to_string()),
                (
                    "discrete_threshold_cells",
                    f(2.0 * r.radius_cells + r.thickness_cells),
                ),
                ("sphere_bnb_nodes_high_water", r.bnb_nodes.to_string()),
                ("ccd_toi_error_cells", f(r.worst_toi_error_cells)),
                ("ccd_toi_checks", r.toi_checks.to_string()),
            ]);
        }
        for r in &c2_rows {
            let (slope, r2) = fit_of(r.kind, r.solver_name, r.shrink);
            rows.push(vec![
                ("arm", "c2_scaling".to_string()),
                ("field", "slab".to_string()),
                ("wall_thickness_cells", f(2.0)),
                ("element_kind", r.kind.to_string()),
                ("solver", r.solver_name.to_string()),
                ("interval_shrink", r.shrink.to_string()),
                ("elements", r.elements.to_string()),
                ("elements_tested", r.tested.to_string()),
                (
                    "broad_phase_cull",
                    f(1.0 - r.tested as f64 / r.elements as f64),
                ),
                ("speed_cells_per_step", f(SPRINT * DT / CELL)),
                ("total_us", f(r.total_ns / 1000.0)),
                (
                    "us_per_element",
                    f(r.total_ns / 1000.0 / r.tested.max(1) as f64),
                ),
                (
                    "us_per_element_in_body",
                    f(r.total_ns / 1000.0 / r.elements as f64),
                ),
                ("marginal_us_per_element", f(slope)),
                ("linear_fit_r2", f(r2)),
                ("cycles", f(r.total_cycles)),
                (
                    "cycles_per_element",
                    f(r.total_cycles / r.tested.max(1) as f64),
                ),
                ("ghz", f(r.total_cycles / r.total_ns)),
                ("sphere_bnb_nodes_high_water", r.bnb_nodes.to_string()),
                (
                    "toi",
                    r.toi.map_or_else(|| NA.to_string(), f),
                ),
            ]);
        }
        for r in &c3_rows {
            let ratio = r.ns_mesh / r.ns_spheres;
            rows.push(vec![
                ("arm", "c3_proxy".to_string()),
                ("field", r.fixture.to_string()),
                (
                    "wall_thickness_cells",
                    r.thickness_cells.map_or_else(|| NA.to_string(), f),
                ),
                ("capsule_half_height_cells", f(r.half_height / CELL)),
                (
                    "capsule_aspect",
                    f((r.half_height + BODY_RADIUS) / BODY_RADIUS),
                ),
                (
                    "bead_spacing_cells",
                    f(2.0 * r.half_height / 7.0 / CELL),
                ),
                (
                    "beads_disjoint",
                    (2.0 * r.half_height / 7.0 > 2.0 * BODY_RADIUS).to_string(),
                ),
                ("orientation_index", r.orientation.to_string()),
                ("axis_x", f(r.axis[0])),
                ("axis_y", f(r.axis[1])),
                ("axis_z", f(r.axis[2])),
                ("speed_cells_per_step", f(r.step_cells)),
                (
                    "toi_capsule_spheres",
                    r.toi_spheres.map_or_else(|| NA.to_string(), f),
                ),
                (
                    "toi_capsule_mesh",
                    r.toi_mesh.map_or_else(|| NA.to_string(), f),
                ),
                (
                    "toi_capsule_mesh_reference",
                    r.toi_mesh_reference
                        .map_or_else(|| NA.to_string(), f),
                ),
                (
                    "toi_disagreement_cells",
                    r.disagreement_cells
                        .map_or_else(|| NA.to_string(), f),
                ),
                (
                    "toi_disagreement_floor_cells",
                    r.disagreement_floor_cells
                        .map_or_else(|| NA.to_string(), f),
                ),
                (
                    "proxy_missed_impact",
                    (r.toi_spheres.is_none() && r.toi_mesh_reference.is_some()).to_string(),
                ),
                (
                    "mesh_missed_impact",
                    (r.toi_mesh_reference.is_none() && r.toi_spheres.is_some()).to_string(),
                ),
                ("mesh_false_negative", r.mesh_false_negative.to_string()),
                ("mesh_false_positive", r.mesh_false_positive.to_string()),
                (
                    "solver_error_cells",
                    r.solver_error_cells
                        .map_or_else(|| NA.to_string(), f),
                ),
                (
                    "support_deficit_spheres_cells",
                    r.deficit_spheres_cells
                        .map_or_else(|| NA.to_string(), f),
                ),
                (
                    "support_deficit_mesh_cells",
                    r.deficit_mesh_cells
                        .map_or_else(|| NA.to_string(), f),
                ),
                ("cost_ratio_spheres_vs_mesh", f(ratio)),
                ("total_us", f(r.ns_mesh / 1000.0)),
                ("us_per_element", f(r.ns_mesh / 1000.0 / 200.0)),
                ("elements", "200".to_string()),
                ("elements_tested", r.tested_mesh.to_string()),
                (
                    "broad_phase_cull",
                    f(1.0 - r.tested_mesh as f64 / 200.0),
                ),
                ("cycles", f(r.cycles_mesh)),
                ("ghz", f(r.cycles_mesh / r.ns_mesh)),
                ("spheres_us", f(r.ns_spheres / 1000.0)),
                ("spheres_ghz", f(r.cycles_spheres / r.ns_spheres)),
                ("solver", c3_solver.name.to_string()),
                ("interval_shrink", "true".to_string()),
            ]);
        }

        let aggregates: Row = vec![
            ("c1_holds", c1_holds.to_string()),
            ("c2_holds", c2_holds.to_string()),
            ("c3_holds", c3_holds.to_string()),
            ("tunnels_discrete_total", tunnels_discrete_total.to_string()),
            ("tunnels_ccd_total", tunnels_ccd_total.to_string()),
            ("solver_chosen", chosen_name.to_string()),
            (
                "solver_cleared_accuracy_bar",
                solver_cleared_bar.to_string(),
            ),
            ("solver_accuracy_bar_cells", f(ACCURACY_BAR_CELLS)),
        ];
        let registered: [&str; 15] = [
            "wall_thickness_cells",
            "speed_cells_per_step",
            "shots",
            "tunnels_discrete",
            "tunnels_ccd",
            "us_per_element",
            "elements",
            "linear_fit_r2",
            "toi_capsule_spheres",
            "toi_capsule_mesh",
            "toi_disagreement_cells",
            "cost_ratio_spheres_vs_mesh",
            "c1_holds",
            "c2_holds",
            "c3_holds",
        ];
        for mut row in rows {
            row.extend(aggregates.iter().cloned());
            for name in registered {
                if !row.iter().any(|(k, _)| *k == name) {
                    row.push((name, NA.to_string()));
                }
            }
            run.record(&row);
        }
    }
}
