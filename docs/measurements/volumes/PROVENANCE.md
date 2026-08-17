# Volumes

**The `.raw` files in this directory are not committed.** `scripts/fetch_volumes.sh` downloads them and
verifies each against the publisher's own SHA-512. This file is committed, and it is what makes the
figures traceable when the bytes are not in git.

```bash
./scripts/fetch_volumes.sh          # fetch anything missing, verify all
./scripts/fetch_volumes.sh --check  # verify what is present, fetch nothing
```

Benchmarks that read these **skip cleanly when they are absent**, so a clean clone with no network still
builds and tests. M-006.

## Why real volumes at all

Every reference field in this crate is **analytic**. That is what makes the accuracy harness exact — a
Hausdorff distance against a closed-form surface — and it is exactly why the timings cannot be set beside
a published isosurfacing benchmark, because those are run on CT and simulation data.

There is a second reason and it is the sharper one. **Quantised data is what makes Grosso 2017's singular
faces reachable.** A singular face needs `v₀·v₂` and `v₁·v₃` bit-identical, which a continuous `f64` field
essentially never produces — measured at **0 of 299,215** over 400,000 random cells (M-220) — while `u8`
voxels collide readily, and M-232 reproduced the order of magnitude synthetically at quantum `1/255`.
Grosso counts **8, 58 and 20** singular faces per 512²×~700 CT volume. Both files here are `uint8`
deliberately, so they are also the fixture **A-002i** and **A-020b** have been waiting for.

## Source

[Open Scientific Visualization Datasets](http://klacansky.com/open-scivis-datasets/), curated by Pavol
Klacansky.

**The site serves over HTTP only.** Port 443 refuses the connection outright — not a timeout, not a
certificate problem, a `connection refused` from `openssl s_client` — so there is no HTTPS URL to prefer
and nothing to fall back from (V-40). Integrity is carried by the **published SHA-512**, which is the
right guarantee for content-addressed data regardless of transport: a wrong or tampered file fails the
hash, is deleted, and the script exits non-zero.

## Format

Stated by the site, and both statements matter here:

> *"All datasets are in little-endian byte order."*
>
> *"Dimensions are width x height x depth (e.g., `array[depth][height][width]` in C)."*

So **x is the fastest-varying axis**, which is exactly the layout
[`isomesh::construct::SampledField`](../../../crates/isomesh/src/construct.rs) documents — *"laid out
x-fastest over `shape`"* — and the same layout `Shape3::linearize` defines. No transpose is needed and
none is done.

`uint8` files are one byte per sample, so the file size is exactly `width · height · depth` and the
reader checks that rather than trusting the filename.

## Files

| file | dimensions | bytes | type |
|---|---|---:|---|
| `fuel_64x64x64_uint8.raw` | 64 × 64 × 64 | 262,144 | `uint8` |
| `bonsai_256x256x256_uint8.raw` | 256 × 256 × 256 | 16,777,216 | `uint8` |

### `fuel`

*Simulation of fuel injection into a combustion chamber.* **The higher the density value, the less
presence of air** — so high is solid, which is why the adapter forms `f = iso − value` to keep this
crate's convention that negative is inside.

**Acknowledgement:** volvis.org and SFB 382 of the German Research Council (DFG).

```
url     http://klacansky.com/open-scivis-datasets/fuel/fuel_64x64x64_uint8.raw
sha512  77fdd7c657da1946bafc84e88c6b8a03ae104a79a5bdec3c7db9257480ef4bf7
        2551a08d22fd237c8e387dd2571b575f1a1a11f5f32b1fa4d4ef385d9fe1d613
```

Chosen for size: 256 kB fetches in a second and every extractor can run it at full resolution, including
subgrid Marching Tetrahedra at ~200× Marching Cubes (M-308).

### `bonsai`

*CT scan of a bonsai tree.* The volume published comparisons actually use, which is the entire point of
M-006.

**Acknowledgement:** volvis.org and S. Roettger, VIS, University of Stuttgart.

```
url     http://klacansky.com/open-scivis-datasets/bonsai/bonsai_256x256x256_uint8.raw
sha512  b34156a0ffc80ffaf84d069f3d05a40fdd999a35f05492829a2b0c13403a3147
        e73712b1d10c2cc34da66a59540a1632dae6adc96f3ebf3efa5d4d6c10598997
```

## Licence and citation

The site does not publish a blanket licence. What it asks is explicit:

> *"Please, cite individual datasets to support authors (some have BibTeX)."*

So each file carries its acknowledgement above, taken verbatim from its dataset page, and **that is why
nothing here is redistributed**: the files are fetched from the source by the person running the
benchmark rather than vendored into this repository, which keeps the question of redistribution rights
from arising at all. If a volume is ever committed, its licence has to be established first — the
per-dataset acknowledgement is an attribution request and is **not** a grant.
