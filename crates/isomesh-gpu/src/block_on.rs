//! Blocking on a future, in the twenty lines it actually takes.
//!
//! `wgpu` hands out futures in three places this crate touches — adapter
//! request, device request, and popping an error scope — and nothing else here
//! is async. Pulling in an executor to await three calls would be a dependency
//! for a thread park, so this is the thread park.
//!
//! Built on [`Wake`] rather than a raw waker vtable because the workspace sets
//! `unsafe_code = "forbid"`, and because the safe version is shorter.

use core::future::Future;
use core::pin::pin;
use core::task::{Context, Poll, Waker};
use std::sync::Arc;
use std::task::Wake;

/// A waker that unparks the thread doing the blocking.
///
/// `park`/`unpark` are permit-based rather than edge-triggered — an unpark
/// arriving before the park is remembered — so a wake landing between `poll`
/// and `park` cannot be lost, and no mutex is needed to close that window.
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
pub(crate) fn block_on<F: Future>(future: F) -> F::Output {
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
