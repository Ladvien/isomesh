//! Getting a field onto the GPU, and results back off it.

use isomesh::Sdf;

use crate::{Error, GridParams, Result};

/// A grid of `f32` field samples in GPU memory, with the grid that describes it.
///
/// The two travel together because they are useless apart: a buffer of floats
/// with no grid cannot be indexed, and a grid with no buffer cannot be sampled.
/// A shader binds this as `array<f32>` alongside
/// [`GridParams::to_std140`] as its uniform.
#[derive(Debug)]
pub struct FieldBuffer {
    buffer: wgpu::Buffer,
    params: GridParams,
}

impl FieldBuffer {
    /// Allocate the samples for `params` without writing them.
    ///
    /// Usable as a compute *output* — a shader that evaluates the field on the
    /// GPU writes here — as well as an input. Usage is `STORAGE | COPY_DST |
    /// COPY_SRC` for exactly that reason: the same buffer is a destination when
    /// the CPU fills it and a source when [`read_buffer`] takes it back.
    #[must_use]
    pub fn new(device: &wgpu::Device, params: GridParams) -> Self {
        let buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("isomesh field samples"),
            size: params.field_buffer_size(),
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_DST
                | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        Self { buffer, params }
    }

    /// Allocate and write `samples`, in `x`-fastest order.
    ///
    /// # Errors
    ///
    /// [`Error::SampleCountMismatch`] if the slice is not exactly
    /// [`GridParams::sample_count`] long. A short slice is rejected rather than
    /// zero-filled and a long one rather than truncated: both repairs produce a
    /// buffer that looks meshable and describes a different surface.
    pub fn uploaded(
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        params: GridParams,
        samples: &[f32],
    ) -> Result<Self> {
        let expected = params.sample_count();
        let got = samples.len() as u64;
        if got != expected {
            return Err(Error::SampleCountMismatch { expected, got });
        }
        let field = Self::new(device, params);
        let mut bytes = Vec::with_capacity(samples.len() * 4);
        for s in samples {
            bytes.extend_from_slice(&s.to_le_bytes());
        }
        queue.write_buffer(&field.buffer, 0, &bytes);
        Ok(field)
    }

    /// Sample an `isomesh` field on the CPU and upload the result.
    ///
    /// The bridge between the two crates, and the thing that makes a GPU
    /// extraction comparable to a CPU one: both read the same numbers, so a
    /// difference in the output is the *algorithm*, which is what GPU-005's
    /// bit-identity acceptance is about.
    ///
    /// # Errors
    ///
    /// Propagates [`uploaded`](Self::uploaded)'s, which cannot actually fire
    /// here — the slice is built from `params` — and is returned rather than
    /// asserted so this stays one code path with the caller-supplied case.
    pub fn sampled<F>(
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        params: GridParams,
        field: &F,
    ) -> Result<Self>
    where
        F: Sdf<Scalar = f32>,
    {
        let [sx, sy, sz] = params.samples();
        let mut samples = Vec::with_capacity(params.sample_count() as usize);
        // x fastest, matching GridParams' documented index order.
        for z in 0..sz {
            for y in 0..sy {
                for x in 0..sx {
                    samples.push(field.sample(params.sample_position([x, y, z])));
                }
            }
        }
        Self::uploaded(device, queue, params, &samples)
    }

    /// The underlying buffer, for binding into a pipeline.
    #[must_use]
    pub const fn buffer(&self) -> &wgpu::Buffer {
        &self.buffer
    }

    /// The grid these samples describe.
    #[must_use]
    pub const fn params(&self) -> GridParams {
        self.params
    }
}

/// Copy `bytes` from `source` back to the CPU as `f32`s.
///
/// Blocks until the copy completes. That is the honest signature for what this
/// does — a caller wanting it off the frame thread runs it off the frame
/// thread — and it is why nothing here is `async`: read-back needs the device
/// polled, not an executor, and hiding the wait behind a future would suggest
/// otherwise.
///
/// `source` needs `COPY_SRC`; a staging buffer is created and dropped per call,
/// which is right for the validation and comparison this is for and wrong for a
/// per-frame path. When a per-frame path exists it gets a persistent staging
/// ring, not a flag on this function.
///
/// # Errors
///
/// [`Error::UnalignedReadback`] if `bytes` is not a multiple of 4,
/// [`Error::MapFailed`] if the map is refused, [`Error::DeviceLost`] if the
/// submission never completes.
pub fn read_buffer(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    source: &wgpu::Buffer,
    bytes: u64,
) -> Result<Vec<f32>> {
    if bytes % 4 != 0 {
        return Err(Error::UnalignedReadback { bytes, stride: 4 });
    }
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("isomesh readback"),
        size: bytes,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("isomesh readback copy"),
    });
    encoder.copy_buffer_to_buffer(source, 0, &staging, 0, bytes);
    queue.submit(Some(encoder.finish()));

    let (sender, receiver) = std::sync::mpsc::channel();
    staging
        .slice(..)
        .map_async(wgpu::MapMode::Read, move |result| {
            // A closed channel means the caller is gone, which cannot happen
            // while this function is on the stack. Dropping the result is the
            // only thing left to do with it and there is nobody to report to.
            let _ = sender.send(result);
        });
    // `Wait` with no submission index waits for everything queued, and with no
    // timeout waits indefinitely -- a read-back that gives up early would
    // return a partially written buffer, which is worse than blocking.
    device
        .poll(wgpu::PollType::Wait {
            submission_index: None,
            timeout: None,
        })
        .map_err(|_| Error::DeviceLost)?;

    match receiver.recv() {
        Ok(Ok(())) => {}
        Ok(Err(_)) => return Err(Error::MapFailed),
        Err(_) => return Err(Error::DeviceLost),
    }

    let out = {
        let view = staging.slice(..).get_mapped_range();
        view.chunks_exact(4)
            .map(|word| f32::from_le_bytes([word[0], word[1], word[2], word[3]]))
            .collect()
    };
    staging.unmap();
    Ok(out)
}

#[cfg(test)]
mod tests;
