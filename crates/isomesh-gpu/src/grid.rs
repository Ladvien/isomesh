//! The grid a GPU extraction samples, and the bytes a shader reads it as.

use crate::{Error, Result};

/// A sampling grid: how many samples, where they start, how far apart.
///
/// Private fields and one checked constructor, for the same reason
/// `isomesh`'s `ValidateConfig` has them — a `GridParams` that exists is one
/// whose arithmetic is meaningful, so nothing downstream needs a runtime guard.
///
/// # Conventions, inherited from `isomesh` and repeated here because a
/// mismatch across this boundary produces plausible garbage
///
/// - **Index order.** `x` varies fastest: `i = x + y·sx + z·sx·sy`.
/// - **Sample position.** `origin + cell_size · [x, y, z]`, computed from the
///   grid origin and the index rather than by accumulation. That is not
///   fussiness: `isomesh`'s M-70 and M-73 both record cracks caused by
///   `(origin + h·i) + h ≠ origin + h·(i + 1)` at a spacing that is not a power
///   of two, and a shader that walks a sample cursor reintroduces exactly that.
/// - **Sign.** Negative inside the solid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridParams {
    samples: [u32; 3],
    origin: [f32; 3],
    cell_size: f32,
}

impl GridParams {
    /// Bytes this occupies in a uniform buffer.
    ///
    /// Two `vec4`s, which is what the layout below packs to. Chosen so std140
    /// and std430 agree and there is no padding rule to get wrong: a shader
    /// declaring `vec4<u32>` then `vec4<f32>` reads exactly these bytes under
    /// either.
    pub const UNIFORM_SIZE: u64 = 32;

    /// A grid, or the reason it is not one.
    ///
    /// # Errors
    ///
    /// [`Error::DegenerateGrid`] if any axis has fewer than two samples — one
    /// sample spans no cell, so there is nothing to extract along that axis.
    /// [`Error::InvalidCellSize`] and [`Error::InvalidOrigin`] for non-finite
    /// or non-positive geometry. [`Error::GridTooLarge`] if the sample count
    /// overflows the byte arithmetic every buffer size here depends on.
    pub fn new(samples: [u32; 3], origin: [f32; 3], cell_size: f32) -> Result<Self> {
        if samples.iter().any(|&s| s < 2) {
            return Err(Error::DegenerateGrid { samples });
        }
        if !cell_size.is_finite() || cell_size <= 0.0 {
            return Err(Error::InvalidCellSize);
        }
        if origin.iter().any(|c| !c.is_finite()) {
            return Err(Error::InvalidOrigin);
        }
        // Checked here so `sample_count` and every buffer size derived from it
        // can be plain arithmetic. Four bytes per sample, so the count itself
        // must leave room for the multiply.
        let count = u64::from(samples[0])
            .checked_mul(u64::from(samples[1]))
            .and_then(|n| n.checked_mul(u64::from(samples[2])))
            .and_then(|n| n.checked_mul(4));
        if count.is_none() {
            return Err(Error::GridTooLarge { samples });
        }
        Ok(Self {
            samples,
            origin,
            cell_size,
        })
    }

    /// Samples along each axis.
    #[must_use]
    pub const fn samples(&self) -> [u32; 3] {
        self.samples
    }

    /// World position of the sample at index `[0, 0, 0]`.
    #[must_use]
    pub const fn origin(&self) -> [f32; 3] {
        self.origin
    }

    /// Distance between adjacent samples.
    #[must_use]
    pub const fn cell_size(&self) -> f32 {
        self.cell_size
    }

    /// Total samples. Never overflows — [`new`](Self::new) rejected the grids
    /// where it would.
    #[must_use]
    pub fn sample_count(&self) -> u64 {
        u64::from(self.samples[0]) * u64::from(self.samples[1]) * u64::from(self.samples[2])
    }

    /// Total cells, which is one fewer than the samples on each axis.
    #[must_use]
    pub fn cell_count(&self) -> u64 {
        u64::from(self.samples[0] - 1)
            * u64::from(self.samples[1] - 1)
            * u64::from(self.samples[2] - 1)
    }

    /// Bytes an `f32` sample buffer for this grid needs.
    #[must_use]
    pub fn field_buffer_size(&self) -> u64 {
        self.sample_count() * 4
    }

    /// World position of the sample at `index`.
    ///
    /// The reference implementation of the rule a shader must follow: multiply
    /// the index, never accumulate. A CPU-side caller uses this to check a
    /// shader against it.
    ///
    /// # `origin + h * i`, as two operations, and not `mul_add`
    ///
    /// This must be the *same expression* `isomesh`'s own extractors evaluate —
    /// `marching_cubes::corner_position` is `origin + cell_size * index`, a
    /// multiply and then an add, rounding twice. `mul_add` rounds **once**, and
    /// the two disagree in the last bit at any spacing where `h * i` is not
    /// exact.
    ///
    /// That is not a stylistic difference. This function decides where the
    /// field is *sampled* before upload, so using the fused form makes the GPU
    /// read a field evaluated at slightly different points from the CPU's, and
    /// every downstream comparison then measures that instead of the algorithm.
    ///
    /// **It was written with `mul_add` and the error was invisible at `h =
    /// 0.125`**, where `h * i` is exact and the two forms agree bit for bit. It
    /// surfaced the moment E-301 ran at `h = 0.1` (M-143). Fifth instance in
    /// this repository of an algebraic identity IEEE does not honour, after
    /// M-32, M-49, M-70 and M-73.
    #[must_use]
    pub fn sample_position(&self, index: [u32; 3]) -> [f32; 3] {
        [
            self.origin[0] + self.cell_size * index[0] as f32,
            self.origin[1] + self.cell_size * index[1] as f32,
            self.origin[2] + self.cell_size * index[2] as f32,
        ]
    }

    /// The bytes a shader reads this as.
    ///
    /// ```wgsl
    /// struct GridParams {
    ///     samples:   vec4<u32>,   // xyz used, w padding
    ///     placement: vec4<f32>,   // xyz origin, w cell size
    /// }
    /// ```
    ///
    /// Packed by hand rather than through a derive, because it is 32 bytes with
    /// no conditional padding and a layout crate would be a dependency carrying
    /// a rule this does not need. Little-endian explicitly: every backend wgpu
    /// targets is little-endian, and saying so beats inheriting the host's.
    #[must_use]
    pub fn to_std140(&self) -> [u8; Self::UNIFORM_SIZE as usize] {
        let mut out = [0u8; Self::UNIFORM_SIZE as usize];
        let words = [
            self.samples[0].to_le_bytes(),
            self.samples[1].to_le_bytes(),
            self.samples[2].to_le_bytes(),
            0u32.to_le_bytes(),
            self.origin[0].to_le_bytes(),
            self.origin[1].to_le_bytes(),
            self.origin[2].to_le_bytes(),
            self.cell_size.to_le_bytes(),
        ];
        for (slot, word) in out.as_chunks_mut::<4>().0.iter_mut().zip(words) {
            slot.copy_from_slice(&word);
        }
        out
    }
}

#[cfg(test)]
mod tests;
