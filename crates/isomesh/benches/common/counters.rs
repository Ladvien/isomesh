//! Hardware performance counters, for the benches that need them.
//!
//! Ticket: M-001, extracted from R-005's `experiment_p12` so there is one
//! implementation rather than two. Two copies of a counter set drift, and the
//! drift is invisible: a bench reporting `cache_misses_per_sample` from its own
//! private probe cannot be compared with one reporting it from another.
//!
//! # Linux only, and the callers say so out loud
//!
//! These come from `perf_event_open`, which is a Linux system call. It needs no
//! privileges at `perf_event_paranoid = 2` — user-space counting of one's own
//! process is permitted — and it has no macOS equivalent a bench can call. This
//! module is therefore `#[cfg(target_os = "linux")]` and every caller either
//! refuses on other platforms (`experiment_p12`) or records the columns as
//! `unavailable` (`family`). Neither invents a number.
//!
//! # Why the set is fixed
//!
//! Zen 3 has six general-purpose counters and this opens exactly six plus one
//! software event, so nothing is multiplexed — [`Counts::worst_ratio`] is what
//! says so rather than hoping. A caller wanting a seventh has to give one up,
//! which is a decision worth making explicitly.
//!
//! `STALLED_CYCLES_BACKEND` would be the seventh and the one that says where the
//! cycles go; **AMD does not map it** and `perf_event_open` answers ENOENT.

use perf_event::events::{Cache, CacheOp, CacheResult, Hardware, Software};
use perf_event::{Builder, Counter};

/// Below this, a counter was multiplexed and its value is an extrapolation.
///
/// Zen 3 has six general-purpose counters and this opens exactly six plus one
/// software event, so nothing should be scheduled out. Asserted rather than
/// hoped: a multiplexed count is a scaled estimate, and reporting one as a
/// measurement is the kind of thing this file exists to prevent.
pub(crate) const MIN_TIME_RATIO: f64 = 0.99;

/// One counter, its label, and what it read.
pub(crate) struct Reading {
    pub(crate) count: u64,
    /// `time_running / time_enabled`. Below 1 means the kernel had to share
    /// the counter and scaled the result.
    pub(crate) ratio: f64,
}

/// The six hardware events and one software event, opened together.
///
/// Zen 3 has six general-purpose counters, so this is exactly full and
/// nothing should be multiplexed; [`MIN_TIME_RATIO`] is what says so rather
/// than hoping. `STALLED_CYCLES_BACKEND` would be the seventh and the one
/// that would settle where the cycles go — it is **not available on this
/// machine**, `perf_event_open` answering ENOENT, because AMD does not map
/// the generic event.
pub(crate) struct Probe {
    cycles: Counter,
    instructions: Counter,
    cache_misses: Counter,
    l1d_read_misses: Counter,
    /// Transparent huge pages are `always` here, so a 67 MB array is ~34
    /// pages and this should stay near zero. Measured rather than assumed,
    /// because "the working set grew" and "the page walk grew" are different
    /// mechanisms with the same symptom.
    dtlb_read_misses: Counter,
    branch_misses: Counter,
    page_faults: Counter,
}

/// What one counted run produced.
pub(crate) struct Counts {
    pub(crate) cycles: Reading,
    pub(crate) instructions: Reading,
    pub(crate) cache_misses: Reading,
    pub(crate) l1d_read_misses: Reading,
    pub(crate) dtlb_read_misses: Reading,
    pub(crate) branch_misses: Reading,
    pub(crate) page_faults: Reading,
}

impl Probe {
    /// Open every counter, or say which one the kernel refused.
    ///
    /// # Panics
    ///
    /// If any counter cannot be opened. An experiment that silently drops an
    /// event reports a column it did not measure.
    pub(crate) fn open() -> Self {
        let hardware = |kind: Hardware| {
            Builder::new()
                .kind(kind)
                .build()
                .unwrap_or_else(|e| panic!("perf_event_open for {kind:?}: {e}"))
        };
        Self {
            cycles: hardware(Hardware::CPU_CYCLES),
            instructions: hardware(Hardware::INSTRUCTIONS),
            cache_misses: hardware(Hardware::CACHE_MISSES),
            l1d_read_misses: Builder::new()
                .kind(Cache {
                    which: perf_event::events::WhichCache::L1D,
                    operation: CacheOp::READ,
                    result: CacheResult::MISS,
                })
                .build()
                .unwrap_or_else(|e| panic!("perf_event_open for L1D read miss: {e}")),
            dtlb_read_misses: Builder::new()
                .kind(Cache {
                    which: perf_event::events::WhichCache::DTLB,
                    operation: CacheOp::READ,
                    result: CacheResult::MISS,
                })
                .build()
                .unwrap_or_else(|e| panic!("perf_event_open for dTLB read miss: {e}")),
            branch_misses: hardware(Hardware::BRANCH_MISSES),
            page_faults: Builder::new()
                .kind(Software::PAGE_FAULTS)
                .build()
                .unwrap_or_else(|e| panic!("perf_event_open for page faults: {e}")),
        }
    }

    fn each(&mut self) -> [&mut Counter; 7] {
        [
            &mut self.cycles,
            &mut self.instructions,
            &mut self.cache_misses,
            &mut self.l1d_read_misses,
            &mut self.dtlb_read_misses,
            &mut self.branch_misses,
            &mut self.page_faults,
        ]
    }

    pub(crate) fn reset_and_enable(&mut self) {
        for counter in self.each() {
            counter.reset().expect("reset");
            counter.enable().expect("enable");
        }
    }

    pub(crate) fn disable(&mut self) {
        for counter in self.each() {
            counter.disable().expect("disable");
        }
    }

    pub(crate) fn read(&mut self) -> Counts {
        fn one(counter: &mut Counter) -> Reading {
            let read = counter.read_count_and_time().expect("read");
            Reading {
                count: read.count,
                ratio: if read.time_enabled == 0 {
                    1.0
                } else {
                    read.time_running as f64 / read.time_enabled as f64
                },
            }
        }
        Counts {
            cycles: one(&mut self.cycles),
            instructions: one(&mut self.instructions),
            cache_misses: one(&mut self.cache_misses),
            l1d_read_misses: one(&mut self.l1d_read_misses),
            dtlb_read_misses: one(&mut self.dtlb_read_misses),
            branch_misses: one(&mut self.branch_misses),
            page_faults: one(&mut self.page_faults),
        }
    }
}

impl Counts {
    /// The worst scheduling ratio of the six.
    pub(crate) fn worst_ratio(&self) -> f64 {
        [
            self.cycles.ratio,
            self.instructions.ratio,
            self.cache_misses.ratio,
            self.l1d_read_misses.ratio,
            self.dtlb_read_misses.ratio,
            self.branch_misses.ratio,
            self.page_faults.ratio,
        ]
        .into_iter()
        .fold(1.0f64, f64::min)
    }
}
