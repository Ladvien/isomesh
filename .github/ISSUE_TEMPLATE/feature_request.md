---
name: Feature request
about: Something the crates should do and don't
title: ""
labels: enhancement
---

**The problem, before the solution**

What are you building, and where does the current API stop short? The repo's tickets are phrased as
problems with acceptance criteria (see `BACKLOG.md`), and requests in that shape travel fastest.

**Which crate**

`isomesh` (core, `no_std`), `isomesh-gpu` (wgpu), or `bevy_isomesh` (the plugin)? Note the
boundaries that are load-bearing: no math library types in the core's public API, no engine types
in `isomesh-gpu`, no `bevy_render` dependency in `bevy_isomesh`.

**Already checked**

`BACKLOG.md` lists open tickets and, at the bottom, the things deliberately out of scope with
reasons (networked editing, Nanite-style simplification, …). If your request is one of those, say
why the recorded reason doesn't cover your case — that is genuinely useful information.
