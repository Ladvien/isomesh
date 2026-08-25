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
    /// # The obvious optimisation here is a pessimisation, measured
    ///
    /// This builds a `Vec<u8>` and hands it to `write_buffer`, which looks like
    /// one copy too many. Three variants were measured at 129³ (M-153) and this
    /// one wins:
    ///
    /// | | upload |
    /// |---|---:|
    /// | this: `Vec<u8>` + `write_buffer` | **8.40 ms** |
    /// | `mapped_at_creation` + `write_iter`, no intermediate | 13.47 ms |
    /// | `mapped_at_creation` + bulk `copy_from_slice` | 8.62 ms |
    ///
    /// Per-element writes into a mapping cost **1.6×**, because mapped memory
    /// may be write-combining — which is exactly why `BufferViewMut` refuses to
    /// deref to `[u8]`. A bulk memcpy into a mapping ties, so `write_buffer` was
    /// already doing the efficient thing and there is nothing to reclaim here.
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
        let mut bytes = Vec::with_capacity(samples.len() * 4);
        for sample in samples {
            bytes.extend_from_slice(&sample.to_le_bytes());
        }
        Ok(Self::from_bytes(device, queue, params, &bytes))
    }

    /// Allocate and write bytes that are already in the shader's layout.
    ///
    /// The one place a `FieldBuffer` is filled, so the descriptor exists once.
    fn from_bytes(
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        params: GridParams,
        bytes: &[u8],
    ) -> Self {
        let field = Self::new(device, params);
        queue.write_buffer(&field.buffer, 0, bytes);
        field
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
        // Bytes directly, rather than a `Vec<f32>` that is converted afterwards.
        // That conversion was a second full pass over 8.4 MB at 129³, 1.14 ms of
        // a 4.65 ms upload (M-152). Evaluation already touches every sample, so
        // writing it in the shader's layout here costs nothing extra.
        let mut bytes = Vec::with_capacity(params.sample_count() as usize * 4);
        // x fastest, matching GridParams' documented index order.
        for z in 0..sz {
            for y in 0..sy {
                for x in 0..sx {
                    let sample = field.sample(params.sample_position([x, y, z]));
                    bytes.extend_from_slice(&sample.to_le_bytes());
                }
            }
        }
        Ok(Self::from_bytes(device, queue, params, &bytes))
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
    let raw = read_bytes(device, queue, source, bytes)?;
    Ok(raw
        .as_chunks::<4>()
        .0
        .iter()
        .map(|w| f32::from_le_bytes([w[0], w[1], w[2], w[3]]))
        .collect())
}

/// The same, as `u32`.
///
/// A separate function rather than a generic or a flag: there are exactly two
/// element types crossing this boundary and naming them is shorter than
/// abstracting over them.
///
/// # Errors
///
/// As [`read_buffer`].
pub fn read_buffer_u32(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    source: &wgpu::Buffer,
    bytes: u64,
) -> Result<Vec<u32>> {
    let raw = read_bytes(device, queue, source, bytes)?;
    Ok(raw
        .as_chunks::<4>()
        .0
        .iter()
        .map(|w| u32::from_le_bytes([w[0], w[1], w[2], w[3]]))
        .collect())
}

/// Copy `bytes` from `source` back to the CPU, untyped.
///
/// One request through [`read_bytes_many`], which is where the staging/map/poll
/// dance actually lives.
///
/// # Errors
///
/// As [`read_buffer`].
pub fn read_bytes(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    source: &wgpu::Buffer,
    bytes: u64,
) -> Result<Vec<u8>> {
    let mut out = read_bytes_many(device, queue, &[(source, bytes)])?;
    match out.pop() {
        Some(only) => Ok(only),
        // `read_bytes_many` returns one result per request and it was given one.
        None => Err(Error::DeviceLost),
    }
}

/// Copy several buffers back to the CPU in **one submission and one device
/// wait**.
///
/// The one place that touches a staging buffer, so the map/poll/unmap dance
/// exists once.
///
/// # Why the plural version is the primitive
///
/// `device.poll(Wait { submission_index: None })` **drains the entire queue**,
/// not just the copy that preceded it — so the cost of a read-back is a full
/// device synchronisation whether it moves four bytes or forty megabytes, and
/// two read-backs in a row cost two of them for no reason. M-167 measured
/// synchronisation at 82.6% of this crate's GPU time, which makes the second
/// drain the expensive part of the operation rather than an detail of it.
///
/// Batching is only legal when nothing between the reads consumes a result. A
/// read-back whose value sizes the next dispatch **cannot** join the batch, and
/// `extract_buffers` has exactly one of those — the triangle total that sizes
/// the geometry buffers.
///
/// # Errors
///
/// [`Error::UnalignedReadback`] if any `bytes` is not a multiple of 4,
/// [`Error::MapFailed`] if a map is refused, [`Error::DeviceLost`] if the
/// submission never completes.
pub fn read_bytes_many(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    requests: &[(&wgpu::Buffer, u64)],
) -> Result<Vec<Vec<u8>>> {
    for (_, bytes) in requests {
        if !bytes.is_multiple_of(4) {
            return Err(Error::UnalignedReadback {
                bytes: *bytes,
                stride: 4,
            });
        }
    }
    if requests.is_empty() {
        return Ok(Vec::new());
    }

    // **One staging buffer, not one per request.** Peak staging memory is the
    // same either way -- every copy in a batch has to be resident until the
    // single wait completes -- but this is one allocation and one map instead of
    // `n` of each, and the device allocator is the thing being avoided.
    let total: u64 = requests.iter().map(|(_, bytes)| *bytes).sum();
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("isomesh readback"),
        size: total,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("isomesh readback copy"),
    });
    let mut at = 0u64;
    for (source, bytes) in requests {
        // Every size is a multiple of 4, checked above, so each offset satisfies
        // `COPY_BUFFER_ALIGNMENT` by construction rather than by rounding.
        encoder.copy_buffer_to_buffer(source, 0, &staging, at, *bytes);
        at += *bytes;
    }
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

    let view = staging.slice(..).get_mapped_range();
    let mut out = Vec::with_capacity(requests.len());
    let mut at = 0usize;
    for (_, bytes) in requests {
        let len = usize::try_from(*bytes).map_err(|_| Error::DeviceLost)?;
        let end = at.checked_add(len).ok_or(Error::DeviceLost)?;
        let slice = view.get(at..end).ok_or(Error::DeviceLost)?;
        out.push(slice.to_vec());
        at = end;
    }
    drop(view);
    staging.unmap();
    Ok(out)
}

/// A read-back in flight, with no wait in it.
///
/// [`read_bytes_many`] waits on an mpsc channel, which is correct on a native
/// thread and fatal on the browser's: under WebGPU `Device::poll` is a
/// documented no-op and `map_async` completes only when control returns to the
/// event loop, so a thread that blocks on `recv()` is the thread that would have
/// run the event loop and the tab deadlocks. This is the same staging copy with
/// the wait taken out -- submit, then ask again next frame.
///
/// One frame of latency at best, more under load, and that is the honest shape
/// of a GPU read-back. [`read_bytes_many`] stays as it is: a native test wants
/// the answer on the next line, and this is its deferred sibling rather than a
/// replacement.
#[derive(Debug)]
pub struct Readback {
    /// `None` for an empty request list. A zero-size `wgpu::Buffer` is invalid,
    /// so the empty case is represented by having no buffer at all rather than
    /// by a buffer that validation would reject.
    staging: Option<wgpu::Buffer>,
    /// 0 pending, 1 mapped, 2 refused. Written once by the `map_async` callback,
    /// which wgpu may run on any thread, and read by [`Readback::ready`] and
    /// [`Readback::take`] -- hence the atomic rather than a `Cell`.
    state: std::sync::Arc<std::sync::atomic::AtomicU8>,
    /// Byte length per request, in request order, so [`Readback::take`] can
    /// split the single mapped range back into the buffers that were asked for.
    spans: Vec<usize>,
}

/// Copy `requests` into one staging buffer and start mapping it.
///
/// The staging allocation, the alignment guard and the copy loop are
/// [`read_bytes_many`]'s; what is missing is its `poll(Wait)` and its
/// `recv()`. The caller asks [`Readback::ready`] once a frame and calls
/// [`Readback::take`] when it says yes.
///
/// # Errors
///
/// [`Error::UnalignedReadback`] if any `bytes` is not a multiple of 4. Nothing
/// else can fail here: a map that is refused is reported by
/// [`Readback::take`], because that is where the caller is asking for the bytes.
pub fn read_bytes_many_deferred(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    requests: &[(&wgpu::Buffer, u64)],
) -> Result<Readback> {
    for (_, bytes) in requests {
        if !bytes.is_multiple_of(4) {
            return Err(Error::UnalignedReadback {
                bytes: *bytes,
                stride: 4,
            });
        }
    }

    let state = std::sync::Arc::new(std::sync::atomic::AtomicU8::new(0));
    if requests.is_empty() {
        // Nothing to map, so nothing to wait for: mark it done here and the
        // caller's `ready` loop terminates on its first look, matching
        // `read_bytes_many`'s empty-input `Ok(Vec::new())`.
        state.store(1, std::sync::atomic::Ordering::Release);
        return Ok(Readback {
            staging: None,
            state,
            spans: Vec::new(),
        });
    }

    let mut spans = Vec::with_capacity(requests.len());
    for (_, bytes) in requests {
        spans.push(usize::try_from(*bytes).map_err(|_| Error::DeviceLost)?);
    }

    // One staging buffer for the batch, for the reason `read_bytes_many` gives:
    // peak residency is the same either way and this is one device allocation
    // and one map instead of `n` of each.
    let total: u64 = requests.iter().map(|(_, bytes)| *bytes).sum();
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("isomesh deferred readback"),
        size: total,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("isomesh deferred readback copy"),
    });
    let mut at = 0u64;
    for (source, bytes) in requests {
        // Every size is a multiple of 4, checked above, so each offset satisfies
        // `COPY_BUFFER_ALIGNMENT` by construction rather than by rounding.
        encoder.copy_buffer_to_buffer(source, 0, &staging, at, *bytes);
        at += *bytes;
    }
    queue.submit(Some(encoder.finish()));

    let signal = std::sync::Arc::clone(&state);
    staging
        .slice(..)
        .map_async(wgpu::MapMode::Read, move |result| {
            // `Release` so that a thread which loads 1 or 2 with `Acquire` also
            // sees the mapped memory the driver wrote before calling this.
            let code = if result.is_ok() { 1 } else { 2 };
            signal.store(code, std::sync::atomic::Ordering::Release);
        });

    Ok(Readback {
        staging: Some(staging),
        state,
        spans,
    })
}

impl Readback {
    /// Whether the mapping has completed. Call once a frame.
    ///
    /// Polls with [`wgpu::PollType::Poll`], which is one call site for both
    /// targets and no `cfg`: on a native backend it is the pump that lets the
    /// `map_async` callback fire, and under WebGPU it is a documented no-op
    /// returning `PollStatus::QueueEmpty` because the browser polls for us and
    /// the callback arrives from the event loop.
    ///
    /// The poll result is dropped deliberately. A device that failed to poll
    /// shows up here as a read-back that never becomes ready, and
    /// [`Readback::take`] is the one place that reports an error, because it is
    /// the one place a caller is asking for bytes.
    #[must_use]
    pub fn ready(&self, device: &wgpu::Device) -> bool {
        let _ = device.poll(wgpu::PollType::Poll);
        self.state.load(std::sync::atomic::Ordering::Acquire) != 0
    }

    /// The bytes, one `Vec` per request, in request order. Unmaps the staging
    /// buffer.
    ///
    /// # Errors
    ///
    /// [`Error::MapFailed`] if the map was refused, [`Error::DeviceLost`] if the
    /// mapping has not completed -- calling this before [`Readback::ready`]
    /// returns `true` is the caller's bug, and returning half a buffer would
    /// hide it.
    pub fn take(self) -> Result<Vec<Vec<u8>>> {
        match self.state.load(std::sync::atomic::Ordering::Acquire) {
            1 => {}
            2 => return Err(Error::MapFailed),
            _ => return Err(Error::DeviceLost),
        }
        let Some(staging) = self.staging else {
            // The empty request list, which never had a buffer to map.
            return Ok(Vec::new());
        };

        let view = staging.slice(..).get_mapped_range();
        let mut out = Vec::with_capacity(self.spans.len());
        let mut at = 0usize;
        for len in &self.spans {
            let end = at.checked_add(*len).ok_or(Error::DeviceLost)?;
            let slice = view.get(at..end).ok_or(Error::DeviceLost)?;
            out.push(slice.to_vec());
            at = end;
        }
        drop(view);
        staging.unmap();
        Ok(out)
    }
}

#[cfg(test)]
mod tests;
