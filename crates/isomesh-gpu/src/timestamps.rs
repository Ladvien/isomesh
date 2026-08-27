//! GPU-side spans for the compute passes, so "execute" is a measurement.
//!
//! Ticket: R-069 (P-71). `M-167` put synchronisation at **83%** of an
//! extraction, and `ExtractTimings`' own documentation says why that number
//! cannot be split any finer from the CPU: *"Wall-clock from the CPU's side.
//! `poll(Wait)` inside the read-backs is what makes these meaningful at all:
//! without a wait the submit returns immediately and every millisecond lands on
//! whichever call happens to block first. Timestamp queries would attribute
//! GPU-side time more precisely and need a device feature this crate does not
//! request."*
//!
//! This is that feature, requested.
//!
//! # One path, and the switch is a constructor rather than a branch
//!
//! [`MarchingCubesGpu::with_timestamps`](crate::MarchingCubesGpu::with_timestamps)
//! builds an extractor that carries a query set; the plain constructor builds
//! one that does not. **The dispatch code is the same code either way** —
//! `timestamp_writes` is a field of `ComputePassDescriptor` and passing `None`
//! or `Some(..)` is passing *data*, not taking a second path. Nothing in the
//! extraction reads a flag and does something different.
//!
//! # Why it refuses rather than degrades
//!
//! An adapter without `TIMESTAMP_QUERY` gets [`Error::TimestampsUnsupported`],
//! not a silent fall back to CPU wall-clock under a GPU-side column name. That
//! is `GPU-007`'s pattern: *a capability check that refuses loudly*. A harness
//! that reported `execute_ms` from `Instant::now()` on a device that cannot
//! measure it would be reporting a column it did not measure, which is the one
//! failure `Run::record`'s missing-column panic exists to prevent, arriving by a
//! different door.
//!
//! # Two timestamps per pass, and the resolve is separate
//!
//! WebGPU's model is: write a `u64` tick into a `QuerySet` slot at the beginning
//! and end of a pass, `resolve_query_set` into a `QUERY_RESOLVE` buffer on the
//! GPU, then copy that to a mappable buffer. The ticks are **not** nanoseconds;
//! `Queue::get_timestamp_period` returns the multiplier, and a period of zero
//! means the driver does not actually implement it — which this module treats as
//! a hard error rather than reporting zero-length spans.

use crate::{Error, Result};

/// The most passes one extraction can time.
///
/// Two per `dispatch` call — `count` and `emit` — and the four entry points use
/// at most two dispatches each. Sized at eight pairs so a caller can run several
/// extractions before resolving, and asserted rather than assumed: `writes`
/// returns `None` once the set is full, so an over-run loses a span instead of
/// corrupting one, and [`Spans::complete`] reports it.
pub const MAX_PASSES: u32 = 8;

/// A query set, its resolve chain, and the cursor into it.
///
/// Not `Clone`: two collectors sharing one query set would interleave their
/// slots and neither would know.
#[derive(Debug)]
pub struct StageTimestamps {
    set: wgpu::QuerySet,
    /// `QUERY_RESOLVE | COPY_SRC`, GPU-only.
    resolve: wgpu::Buffer,
    /// `MAP_READ | COPY_DST`, the only mappable one.
    readback: wgpu::Buffer,
    /// Nanoseconds per tick, from `Queue::get_timestamp_period`.
    period_ns: f32,
    /// Next free slot. Two are consumed per pass.
    ///
    /// **Atomic, not `Cell`, and that is not gold-plating.** `MarchingCubesGpu`
    /// carries one of these and `bevy_isomesh` holds a `MarchingCubesGpu` as a
    /// Bevy `Resource`, which requires `Send + Sync`. A `Cell` here compiles
    /// fine in this crate and breaks that crate's examples — which is exactly
    /// what it did, and `bevy: check --all-targets` is the gate that said so
    /// (M-293's step, earning its keep).
    next: core::sync::atomic::AtomicU32,
    /// Labels, in the order the passes were opened.
    labels: std::sync::Mutex<Vec<&'static str>>,
}

/// One pass's GPU-side span.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Span {
    /// What the pass was called.
    pub label: &'static str,
    /// End minus begin, in milliseconds.
    pub ms: f64,
}

/// Every span resolved from one collector.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Spans {
    /// In the order the passes ran.
    pub spans: Vec<Span>,
    /// Nanoseconds per tick on this device.
    pub period_ns: f64,
    /// `false` if a pass was opened after the set filled up, so its span is
    /// missing from `spans` rather than silently merged into a neighbour's.
    pub complete: bool,
}

impl Spans {
    /// Every span summed.
    #[must_use]
    pub fn total_ms(&self) -> f64 {
        self.spans.iter().map(|s| s.ms).sum()
    }
}

impl StageTimestamps {
    /// Open a query set on `device`.
    ///
    /// # Errors
    ///
    /// [`Error::TimestampsUnsupported`] if the device was not created with
    /// [`wgpu::Features::TIMESTAMP_QUERY`], or if the queue reports a timestamp
    /// period of zero — which means the driver advertises the feature and does
    /// not implement it, and is indistinguishable from a working device that
    /// measures every pass as instantaneous.
    pub fn new(device: &wgpu::Device, queue: &wgpu::Queue) -> Result<Self> {
        if !device.features().contains(wgpu::Features::TIMESTAMP_QUERY) {
            return Err(Error::TimestampsUnsupported);
        }
        let period_ns = queue.get_timestamp_period();
        // `period_ns > 0.0` written as a match on the ordering, so a NaN period
        // — a driver returning garbage rather than zero — is also refused rather
        // than sliding through a negated comparison.
        if !matches!(
            period_ns.partial_cmp(&0.0),
            Some(core::cmp::Ordering::Greater)
        ) {
            return Err(Error::TimestampsUnsupported);
        }
        let count = MAX_PASSES * 2;
        let bytes = u64::from(count) * 8;
        Ok(Self {
            set: device.create_query_set(&wgpu::QuerySetDescriptor {
                label: Some("isomesh stage timestamps"),
                ty: wgpu::QueryType::Timestamp,
                count,
            }),
            resolve: device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("isomesh timestamp resolve"),
                size: bytes,
                usage: wgpu::BufferUsages::QUERY_RESOLVE | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            }),
            readback: device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("isomesh timestamp readback"),
                size: bytes,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }),
            period_ns,
            next: core::sync::atomic::AtomicU32::new(0),
            labels: std::sync::Mutex::new(Vec::new()),
        })
    }

    /// The `timestamp_writes` for the next pass, or `None` if the set is full.
    ///
    /// Returning `None` rather than panicking is deliberate and is not a
    /// fallback: the caller passes it straight into `ComputePassDescriptor`,
    /// which takes an `Option` anyway, so a full set costs a missing span and
    /// [`Spans::complete`] says so. A panic here would take down an extraction
    /// to protect a measurement.
    pub fn writes(&self, label: &'static str) -> Option<wgpu::ComputePassTimestampWrites<'_>> {
        use core::sync::atomic::Ordering;
        // `fetch_update` rather than a load-then-store: two threads recording
        // passes on one extractor would otherwise hand out the same slot pair
        // and neither would know. The extraction path is single-threaded today
        // and this costs nothing; a silently shared slot would cost a wrong
        // number that looks like a right one.
        let at = self
            .next
            .fetch_update(Ordering::SeqCst, Ordering::SeqCst, |at| {
                (at + 2 <= MAX_PASSES * 2).then_some(at + 2)
            })
            .ok()?;
        self.labels.lock().expect("timestamp labels").push(label);
        Some(wgpu::ComputePassTimestampWrites {
            query_set: &self.set,
            beginning_of_pass_write_index: Some(at),
            end_of_pass_write_index: Some(at + 1),
        })
    }

    /// Resolve every recorded pass and read the ticks back, then reset.
    ///
    /// One submission and one wait, deliberately: this is a measurement path and
    /// it runs once per extraction, so it is not the thing being optimised. It
    /// is also **after** everything it measures, so the wait it costs is not
    /// inside any span it reports.
    ///
    /// # Errors
    ///
    /// [`Error::MapFailed`] if the readback map is refused, [`Error::DeviceLost`]
    /// if the submission never completes, and [`Error::TimestampsUnsupported`]
    /// if a resolved pair is non-increasing — a span that ends before it begins
    /// is a driver that is not measuring, and reporting it as a negative
    /// millisecond would put a number in a column that means nothing.
    pub fn resolve(&self, device: &wgpu::Device, queue: &wgpu::Queue) -> Result<Spans> {
        let passes = self.next.load(core::sync::atomic::Ordering::SeqCst) / 2;
        if passes == 0 {
            return Ok(Spans {
                spans: Vec::new(),
                period_ns: f64::from(self.period_ns),
                complete: true,
            });
        }
        let used = passes * 2;
        let bytes = u64::from(used) * 8;

        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("isomesh timestamp resolve"),
        });
        encoder.resolve_query_set(&self.set, 0..used, &self.resolve, 0);
        encoder.copy_buffer_to_buffer(&self.resolve, 0, &self.readback, 0, bytes);
        queue.submit(Some(encoder.finish()));

        let (sender, receiver) = std::sync::mpsc::channel();
        self.readback
            .slice(..bytes)
            .map_async(wgpu::MapMode::Read, move |r| {
                let _ = sender.send(r);
            });
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

        let ticks: Vec<u64> = {
            let view = self.readback.slice(..bytes).get_mapped_range();
            view.as_chunks::<8>()
                .0
                .iter()
                .map(|c| u64::from_le_bytes(*c))
                .collect()
        };
        self.readback.unmap();

        let labels = self.labels.lock().expect("timestamp labels").clone();
        let mut spans = Vec::with_capacity(passes as usize);
        for (index, label) in labels.iter().enumerate() {
            let (begin, end) = (ticks[index * 2], ticks[index * 2 + 1]);
            if end < begin {
                return Err(Error::TimestampsUnsupported);
            }
            spans.push(Span {
                label,
                ms: (end - begin) as f64 * f64::from(self.period_ns) / 1.0e6,
            });
        }
        let complete = passes <= MAX_PASSES;
        self.next.store(0, core::sync::atomic::Ordering::SeqCst);
        self.labels.lock().expect("timestamp labels").clear();
        Ok(Spans {
            spans,
            period_ns: f64::from(self.period_ns),
            complete,
        })
    }
}
