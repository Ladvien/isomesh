//! A device with no window, and a report of what it turned out to be.
//!
//! Two callers need this and neither is a game. A **CAD tool** has no renderer
//! and still wants GPU extraction. A **test** has no display server at all —
//! GPU-004's ticket is explicit that the harness comes first, *"no Bevy in the
//! room; if it can't run against raw wgpu, the abstraction leaked."*
//!
//! A Bevy consumer does **not** use this. It already has a device and passes it
//! in, which is the whole point of the API rule in the crate docs.

use crate::{Error, Result};

/// An instance, adapter, device and queue, created without a surface.
///
/// Owns them, so a caller can keep one alive for the length of a session.
/// Everything else in this crate borrows a `&Device` and does not care where it
/// came from.
#[derive(Debug)]
pub struct Gpu {
    device: wgpu::Device,
    queue: wgpu::Queue,
    report: AdapterReport,
}

/// What the adapter actually turned out to be.
///
/// Recorded because a GPU measurement without it is not reproducible, and
/// because "which backend did this run on" is the first question of any
/// disagreement between two machines. This crate's own tests print it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AdapterReport {
    /// Driver's name for the device.
    pub name: String,
    /// Backend the device was created on.
    pub backend: wgpu::Backend,
    /// Discrete, integrated, virtual, CPU, or unknown.
    pub device_type: wgpu::DeviceType,
    /// Driver identification, as reported.
    pub driver: String,
    /// Largest single storage buffer binding this device will accept, in bytes.
    ///
    /// The limit that decides how big a grid one dispatch can cover, which is
    /// why it is the one surfaced here rather than the whole `Limits`.
    pub max_storage_buffer_binding_size: u64,
    /// Largest total workgroups per dispatch, per axis.
    pub max_compute_workgroups_per_dimension: u32,
}

impl Gpu {
    /// Ask for a high-performance adapter and open a device on it.
    ///
    /// # No software fallback, deliberately
    ///
    /// `fallback_adapter` stays `false` and a missing adapter is
    /// [`Error::NoAdapter`] rather than a CPU reference driver. A benchmark that
    /// silently ran on lavapipe reports numbers three orders of magnitude off
    /// and looks exactly like a slow GPU — the whole reason to have this crate
    /// is to compare against the CPU path, and a fallback would make that
    /// comparison meaningless without saying so.
    ///
    /// # Errors
    ///
    /// [`Error::NoAdapter`] when nothing matches, [`Error::DeviceUnavailable`]
    /// when the adapter refuses the request.
    pub fn new() -> Result<Self> {
        // No display handle: this is the whole point of the module. Backend
        // selection still honours `WGPU_BACKEND` through the `_from_env`
        // variant, so a machine with two backends can be told which to use
        // without this crate inventing a preference.
        let instance =
            wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle_from_env());
        let adapter =
            pollster_lite::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                force_fallback_adapter: false,
                compatible_surface: None,
            }))
            .map_err(|_| Error::NoAdapter)?;

        let info = adapter.get_info();
        let limits = adapter.limits();
        let report = AdapterReport {
            name: info.name,
            backend: info.backend,
            device_type: info.device_type,
            driver: info.driver,
            max_storage_buffer_binding_size: limits.max_storage_buffer_binding_size,
            max_compute_workgroups_per_dimension: limits.max_compute_workgroups_per_dimension,
        };

        let (device, queue) =
            pollster_lite::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
                label: Some("isomesh headless"),
                required_features: wgpu::Features::empty(),
                // The adapter's own limits rather than downlevel defaults: this
                // is a native compute path, and asking for less would cap the
                // grid size for the benefit of a target this crate does not
                // support.
                required_limits: limits,
                memory_hints: wgpu::MemoryHints::Performance,
                // Nothing here needs an experimental feature, and asking for
                // one that a driver half-implements is how a GPU path starts
                // producing results that differ per machine.
                experimental_features: wgpu::ExperimentalFeatures::disabled(),
                trace: wgpu::Trace::Off,
            }))
            .map_err(|_| Error::DeviceUnavailable)?;

        Ok(Self {
            device,
            queue,
            report,
        })
    }

    /// The device, to hand to anything in this crate.
    #[must_use]
    pub const fn device(&self) -> &wgpu::Device {
        &self.device
    }

    /// The queue.
    #[must_use]
    pub const fn queue(&self) -> &wgpu::Queue {
        &self.queue
    }

    /// What this turned out to run on.
    #[must_use]
    pub const fn report(&self) -> &AdapterReport {
        &self.report
    }
}

/// Blocking on a future, in the twenty lines it actually takes.
///
/// `wgpu`'s adapter and device requests are futures and nothing else in this
/// crate is async. Pulling in an executor to await two calls made once per
/// session would be a dependency for a thread park, so this is the thread park.
mod pollster_lite {
    use core::future::Future;
    use core::pin::pin;
    use core::task::{Context, Poll, Waker};
    use std::sync::Arc;
    use std::task::Wake;

    /// A waker that unparks the thread doing the blocking.
    ///
    /// `park`/`unpark` are permit-based rather than edge-triggered — an unpark
    /// arriving before the park is remembered — so a wake landing between
    /// `poll` and `park` cannot be lost, and no mutex is needed to close that
    /// window.
    struct ThreadWaker(std::thread::Thread);

    impl Wake for ThreadWaker {
        fn wake(self: Arc<Self>) {
            self.0.unpark();
        }

        fn wake_by_ref(self: &Arc<Self>) {
            self.0.unpark();
        }
    }

    /// Drive `future` on this thread until it completes.
    ///
    /// Built on [`Wake`] rather than a raw vtable because the workspace forbids
    /// `unsafe_code`, and because the safe version is shorter.
    pub(super) fn block_on<F: Future>(future: F) -> F::Output {
        let mut future = pin!(future);
        let waker = Waker::from(Arc::new(ThreadWaker(std::thread::current())));
        let mut context = Context::from_waker(&waker);
        loop {
            match future.as_mut().poll(&mut context) {
                Poll::Ready(value) => return value,
                Poll::Pending => std::thread::park(),
            }
        }
    }
}

#[cfg(test)]
mod tests;
