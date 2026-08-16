# Bevy Assets submission

`bevy_isomesh.toml` is the entry for [`bevyengine/bevy-assets`](https://github.com/bevyengine/bevy-assets),
staged here rather than submitted. **Opening the pull request is a decision for the maintainer**, not
something a build step should do: it is outward-facing, it happens under a personal GitHub account, and
listing a crate publicly is a claim that it is ready for other people to use.

## The schema is verified, not guessed

Checked against the live repository on 2026-08-15 rather than written from memory, which is B-011's
whole risk — a submission in the wrong shape is a rejected PR and a wasted round trip.

- Directory layout is `Assets/<Category>/`, and `Assets/` carries `2D`, `3D`, `Animation`, `Physics`,
  `Shapes` and twenty-odd others.
- Each entry is one `.toml` file, optionally beside an image of the same stem.
- Fields, from the repository's own README: `name`, `description` (under 100 characters, no
  formatting) and `link` are required; `image` and `crate` are optional.

Two real entries, read to confirm the shape:

```toml
name = "bevy_mod_outline"
description = "A plugin for drawing outlines around meshes"
link = "https://github.com/komadori/bevy_mod_outline"
image = "bevy_mod_outline.png"
crate = "bevy_mod_outline"
```

## To submit

1. Fork `bevyengine/bevy-assets`.
2. Copy `bevy_isomesh.toml` to `Assets/3D/bevy_isomesh.toml`. **`3D` is the judgement call** — the
   crate meshes 3D volumes, and `Shapes` is closer to primitive geometry helpers than to extraction.
   `bevy_hanabi`, `bevy_mod_outline` and `bevy-hikari` all sit in `3D`, which is the company this
   belongs in.
3. Optionally add `bevy_isomesh.png` beside it — 16:9, at most 600 px wide and 2 MB. **There is no
   image yet**: E-214 is blocked on a capture environment, so the crate has no committed screenshot
   of its own that fits. The entry is valid without one and an image can be added later.
4. Open the pull request.

## Keeping it true

The description is under 100 characters, as required. If the crate's scope changes, this file and the
listing both need updating — the listing is a copy, not a link, so it rots independently of this
repository and nothing here can detect that.
