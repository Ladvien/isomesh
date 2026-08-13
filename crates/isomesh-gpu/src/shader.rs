//! WGSL composition: `#include` and `#ifdef`, and nothing else.
//!
//! # Why this and not `naga_oil`
//!
//! `docs/research/2026-08-11-meshing-crate-architecture.md` §5 ranks three
//! options and puts a small preprocessor first. `naga_oil` is Bevy-owned and
//! was 14 months stale at the time of that review, and Bevy has said WESL is
//! where it is going — so depending on it couples this crate to the *outgoing*
//! path of an engine it is not supposed to know about. WESL is the answer when
//! the shader tree gets big, and its README is candid that the features a
//! meshing kernel would most want are still in the planned column.
//!
//! What is actually needed is `#include` and a boolean. That is this file.
//!
//! # The whole language
//!
//! | directive | meaning |
//! |---|---|
//! | `#include <name>` | paste module `name`, **at most once per composition** |
//! | `#ifdef NAME` / `#ifndef NAME` | keep the following lines if `NAME` is (not) defined |
//! | `#else` | invert, once per `#ifdef` |
//! | `#endif` | close |
//!
//! Directives must be the first non-space text on their line. Everything else
//! is passed through untouched, so the input stays valid WGSL that an editor
//! can highlight and `naga` can parse once composed.
//!
//! **`#include` pastes a module at most once.** WGSL has no forward
//! declarations and a duplicated function is a hard error, so "include the
//! shared header from both of the two files that need it" has to work — and the
//! only semantics under which it does is include-once. A *cycle* is still an
//! error rather than being silently absorbed by that rule, because a cycle is a
//! bug about which module owns what and quietly resolving it would hide the
//! question.
//!
//! ```
//! use isomesh_gpu::Composer;
//!
//! # fn main() -> Result<(), isomesh_gpu::Error> {
//! let mut composer = Composer::new();
//! composer.insert("maths", "fn double(x: f32) -> f32 { return x * 2.0; }\n");
//! composer.insert(
//!     "kernel",
//!     "#include <maths>\n#ifdef WIDE\nconst WIDTH: u32 = 256u;\n#else\nconst WIDTH: u32 = 64u;\n#endif\n",
//! );
//!
//! let narrow = composer.compose("kernel", &[])?;
//! assert!(narrow.contains("64u") && !narrow.contains("256u"));
//! assert!(narrow.contains("fn double"));
//!
//! let wide = composer.compose("kernel", &["WIDE"])?;
//! assert!(wide.contains("256u") && !wide.contains("64u"));
//! # Ok(())
//! # }
//! ```

use std::collections::{BTreeMap, BTreeSet};

use crate::{Error, Result};

/// The shader-side counterpart of [`GridParams`](crate::GridParams).
///
/// Included by every kernel here. Its layout has to agree with
/// [`GridParams::to_std140`](crate::GridParams::to_std140) byte for byte, which
/// is why the two live next to each other in the repository and why the WGSL
/// says so in its own comments.
pub const GRID_WGSL: &str = include_str!("shaders/grid.wgsl");

/// A registry of WGSL modules, and the preprocessor over them.
///
/// Modules are `&'static str` because they come from `include_str!` — compiled
/// in, so a shipped binary cannot be missing one, and no file is read at run
/// time.
#[derive(Clone, Debug, Default)]
pub struct Composer {
    modules: BTreeMap<String, &'static str>,
}

/// One `#ifdef` region being tracked.
struct Region {
    /// Whether output was enabled *outside* this region.
    outer: bool,
    /// Whether any branch of this region has been taken yet.
    taken: bool,
    /// Whether `#else` has already been seen, so a second one is an error.
    closed: bool,
}

impl Composer {
    /// An empty registry.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// A registry holding this crate's own modules, ready for a kernel to
    /// `#include <grid>`.
    #[must_use]
    pub fn with_builtins() -> Self {
        let mut composer = Self::new();
        composer.insert("grid", GRID_WGSL);
        composer
    }

    /// Register `source` under `name`, replacing any previous entry.
    pub fn insert(&mut self, name: &str, source: &'static str) {
        self.modules.insert(String::from(name), source);
    }

    /// Module names currently registered, in sorted order.
    ///
    /// Sorted because it feeds GPU-003's permutation sweep, and a validation
    /// run that visits modules in a different order each time reports its
    /// failures in a different order each time.
    #[must_use]
    pub fn module_names(&self) -> Vec<&str> {
        self.modules.keys().map(String::as_str).collect()
    }

    /// Expand `name` with `defines` set.
    ///
    /// # Errors
    ///
    /// [`Error::ShaderModuleMissing`] for an unregistered name — including the
    /// root — [`Error::ShaderCircularInclude`] for a module that includes
    /// itself through any chain, and [`Error::ShaderDirective`] for a malformed
    /// or unbalanced directive.
    pub fn compose(&self, name: &str, defines: &[&str]) -> Result<String> {
        let defines: BTreeSet<&str> = defines.iter().copied().collect();
        let mut out = String::new();
        let mut included = BTreeSet::new();
        let mut in_progress = Vec::new();
        self.expand(name, &defines, &mut out, &mut included, &mut in_progress)?;
        Ok(out)
    }

    fn expand(
        &self,
        name: &str,
        defines: &BTreeSet<&str>,
        out: &mut String,
        included: &mut BTreeSet<String>,
        in_progress: &mut Vec<String>,
    ) -> Result<()> {
        let Some(source) = self.modules.get(name) else {
            return Err(Error::ShaderModuleMissing {
                name: String::from(name),
            });
        };
        in_progress.push(String::from(name));
        included.insert(String::from(name));

        let mut regions: Vec<Region> = Vec::new();
        let mut emitting = true;

        for (number, line) in source.lines().enumerate() {
            let trimmed = line.trim_start();
            let here = || Error::ShaderDirective {
                module: String::from(name),
                line: number + 1,
            };

            // Split `#keyword rest` without requiring the space, so a bare
            // `#ifdef` is a *malformed directive* rather than a line of text
            // that happens to start with a hash. An unrecognised `#word` is
            // passed through untouched — WGSL has no preprocessor of its own
            // today, and inventing an error for syntax it might grow is not
            // this file's business.
            let mut handled = false;
            if let Some(rest) = trimmed.strip_prefix('#') {
                let split = rest.find(char::is_whitespace).unwrap_or(rest.len());
                let (keyword, argument) = rest.split_at(split);
                let argument = argument.trim();
                handled = true;
                match keyword {
                    "ifdef" | "ifndef" => {
                        if argument.is_empty() || argument.split_whitespace().count() != 1 {
                            return Err(here());
                        }
                        let hit = defines.contains(argument) == (keyword == "ifdef");
                        regions.push(Region {
                            outer: emitting,
                            taken: emitting && hit,
                            closed: false,
                        });
                        emitting = emitting && hit;
                    }
                    "else" => {
                        if !argument.is_empty() {
                            return Err(here());
                        }
                        let Some(region) = regions.last_mut() else {
                            return Err(here());
                        };
                        if region.closed {
                            return Err(here());
                        }
                        region.closed = true;
                        emitting = region.outer && !region.taken;
                        region.taken |= emitting;
                    }
                    "endif" => {
                        if !argument.is_empty() {
                            return Err(here());
                        }
                        let Some(region) = regions.pop() else {
                            return Err(here());
                        };
                        emitting = region.outer;
                    }
                    "include" => {
                        let target = argument
                            .strip_prefix('<')
                            .and_then(|r| r.strip_suffix('>'))
                            .ok_or_else(here)?
                            .trim();
                        if target.is_empty() {
                            return Err(here());
                        }
                        if emitting {
                            // The cycle check comes **before** the include-once
                            // skip, and the order is the whole point: with it
                            // the other way round, `a` includes `b` includes
                            // `a` finds `a` already included, skips, and
                            // reports success on a genuinely circular graph.
                            // A test caught exactly that.
                            if in_progress.iter().any(|n| n == target) {
                                return Err(Error::ShaderCircularInclude {
                                    name: String::from(target),
                                });
                            }
                            if !included.contains(target) {
                                self.expand(target, defines, out, included, in_progress)?;
                            }
                        }
                    }
                    _ => handled = false,
                }
            }

            if !handled && emitting {
                out.push_str(line);
                out.push('\n');
            }
        }

        if !regions.is_empty() {
            return Err(Error::ShaderDirective {
                module: String::from(name),
                // The count of lines, i.e. one past the end: the `#endif` that
                // is missing would have gone here.
                line: source.lines().count() + 1,
            });
        }
        in_progress.pop();
        Ok(())
    }
}

#[cfg(test)]
mod tests;
