# Experiments that held, and what they buy you

Moved from the root README on 2026-10-31, with its link paths updated and stale Bevy versions corrected, when the README was brought to the [org documentation guide](https://github.com/bevy-coven/.github/blob/main/DOCUMENTATION.md).

Twenty times, this project wrote down what it expected to measure — and what result would prove it
wrong — **before running the measurement**. Thirteen of those are enforced by the compiler:

```rust,ignore
// This line does not compile if "P-8" is not in `PREREGISTERED`.
let p = isomesh::experiment!("P-8");
```

An unregistered id is a **compile error**, not a test failure. `Preregistration` is `#[non_exhaustive]`
so the macro is the only way to obtain one. And `falsified_by` is not optional — there is no way to
express a registration without it, which is the point of the field existing.

Registering is a commit, so git carries the ordering that prose cannot.

| prediction | verdict | what it means if you depend on this |
|---|---|---|
| one canonical `world_of_sample` closes seam cracks at every cell size | **held** | Chunk seams are exact **only at power-of-two cell sizes** today — 0 unmatched edges there, 63–348 at `0.1`. The weld hides it; an unwelded consumer gets it in full |
| over half the dual mesher's cycles are in `emit_quads` | **held**, `r²` = 0.995 | 45% of instructions, **82% of the time**. This is the measurement that made the 4.26× below possible |
| the k-way weld's order changes its output | **held** | It is not confluent, so vertex order is pinned deliberately. Your meshes reproduce byte-for-byte |
| sharp-feature angles are a property of the field, not a defect with a location | **held**, after a bug reversed two wrong verdicts | Nothing to fix; something to know |
| a link-condition-gated weld reaches zero non-manifold output | **falsified**, both clauses | The gate makes it **worse**. The crate does not ship one, and now you know why rather than wondering |
| the dual's superlinearity is a cache-miss effect | **falsified** … then true | A null measured under a dominant confound is a statement about the confound |
| Manifold Dual Contouring's residue is an interior ambiguity | **falsified** | Bounded exhaustively instead, then reduced to a 48-sample fixture |

Five of the ten held, four were falsified, one never ran. A section called "experiments that held"
which quietly dropped the failures would be exactly the failure this machinery exists to prevent, so
they are all here.

**Three that are worth your time:**

- **The 4.26×, and not one triangle changed.** Making the dual's loop axis a `const` generic and
  forcing its sample row length odd took Surface Nets at 256³ from **693.8 ms to 162.7 ms**, IPC from
  **1.20 to 4.09**, and `SN/MC` from **5.43× to 1.26×**. Golden hashes untouched. It also falsified
  ✗14, this repository's own claim that Surface Nets *"never wins on Zen 3, at any resolution"* — it
  now wins at 16³, 24³ and 32³.
- **One bit of a row length cost 3.4× at 128³**, the canonical voxel chunk size, because a 64 KiB plane
  stride is a cache-set aliasing period. Diagnosed before any code changed, with a `z`-axis control
  that adds the same work and keeps the whole penalty.
- **A reference gradient was returning noise**, and it had recorded two true hypotheses as falsified.
  Past-90° vertices went 6,959 → **472** when it was fixed. Both falsified entries are still in the
  ledger under a banner, because deleting them would delete the evidence that the instrument can be
  the broken thing.

Every number in these docs carries a tier: **M** measured here · **V** verified from a primary source ·
**R** reported · **F** folklore · **✗** falsified. Never cite an R or F claim as justification for a
design decision without saying which it is.

**[The full record →](../experiments.md)** · [`FINDINGS.md`](../../FINDINGS.md), 363 entries
