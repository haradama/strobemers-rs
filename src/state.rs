//! Everything [`MinStrobes`](crate::MinStrobes) and
//! [`RandStrobes`](crate::RandStrobes) have in common.
//!
//! The two iterators differ only in how they pick the next strobe inside a
//! window. Where the windows lie, when iteration stops, how many items remain
//! and the whole public surface are identical, and live here.

use crate::{
    Result, StrobeError, constants::DEFAULT_PRIME_NUMBER, hashes::HashBuf, util::roundup64,
};

/// Parameters, k-mer hashes and position of a strobemer iterator.
#[derive(Debug, Clone)]
pub(crate) struct State<'a> {
    /// Order of the strobemer: 2 or 3.
    pub(crate) n: u8,
    pub(crate) w_min: usize,
    pub(crate) w_max: usize,
    /// Offset from `m1` to the end of the last search window when nothing
    /// clips it: `(n - 1) * w_max`.
    pub(crate) reach: usize,

    pub(crate) hashes: HashBuf<'a>,

    /// Index of the next first strobe (m1).
    pub(crate) idx: usize,
    /// One past the last index that still yields a strobemer; see
    /// [`State::update_stop`].
    stop: usize,
    /// Last index at which `n` strobes fit, or `None` when the sequence cannot
    /// hold them at all. The Go reference stores this in a signed `int` and
    /// lets it go negative; `usize` would wrap, so that case is explicit.
    end_idx: Option<usize>,
    /// Index of the last k-mer.
    end_hash: usize,

    /// Positions of the second and third strobe of the last item.
    pub(crate) idx2: usize,
    pub(crate) idx3: usize,

    /// Mask applied when candidates are compared by `(base + hash) & prime`.
    pub(crate) prime: u64,
    /// Whether windows may be clipped at the end of the sequence.
    shrink: bool,
}

impl<'a> State<'a> {
    /// `seq_len` and the parameters must already have passed
    /// [`validate_params`](crate::util::validate_params).
    pub(crate) fn new(
        seq_len: usize,
        n: u8,
        k: usize,
        w_min: usize,
        w_max: usize,
        hashes: HashBuf<'a>,
    ) -> Self {
        let mut state = Self {
            n,
            w_min,
            w_max,
            reach: (n as usize - 1) * w_max,
            hashes,
            idx: 0,
            stop: 0,
            // `None` when the sequence is shorter than `n * k`.
            end_idx: seq_len.checked_sub(n as usize * k),
            end_hash: seq_len - k,
            idx2: 0,
            idx3: 0,
            prime: DEFAULT_PRIME_NUMBER,
            shrink: true,
        };
        state.update_stop();
        state
    }

    /// Recomputes where iteration ends.
    ///
    /// A position yields a strobemer while all `n` strobes fit and its last
    /// search window is usable: non-empty when windows may shrink, entirely
    /// inside the sequence otherwise. Both bounds only tighten as the index
    /// grows, so a single cut-off describes them, and the remaining length is
    /// exact.
    ///
    /// - order 2 needs `idx + w_min ≤ end_hash` when shrinking, `idx + w_max ≤ end_hash` otherwise
    /// - order 3 needs `idx + w_max + w_min ≤ end_hash` when shrinking, `idx + 2·w_max ≤ end_hash` otherwise
    ///
    /// When `w_min > (n - 1) * k` the window bound is the tighter one. The Go
    /// reference keeps going there and silently reuses the previous strobe
    /// position; stopping means every reported position is a real k-mer.
    fn update_stop(&mut self) {
        let needed = if self.shrink {
            self.reach - self.w_max + self.w_min
        } else {
            self.reach
        };
        self.stop = match (self.end_idx, self.end_hash.checked_sub(needed)) {
            (Some(end_idx), Some(limit)) => end_idx.min(limit) + 1,
            _ => 0,
        };
    }

    /// Starts the next item.
    ///
    /// Returns the offset from `m1` to the end of its last search window —
    /// `reach`, or less if the sequence ends first — after making sure every
    /// hash up to there is buffered. `None` once iteration is over.
    #[inline(always)]
    pub(crate) fn begin(&mut self) -> Option<usize> {
        if self.idx >= self.stop {
            return None;
        }
        let hi = self.reach.min(self.end_hash - self.idx);
        self.hashes.ensure(self.idx, self.idx + hi);
        Some(hi)
    }

    #[inline]
    pub(crate) fn remaining(&self) -> usize {
        self.stop.saturating_sub(self.idx)
    }

    pub(crate) fn set_prime(&mut self, q: u64) -> Result<()> {
        if q < 256 {
            return Err(StrobeError::PrimeNumberTooSmall);
        }
        // Round up to the next power of two, subtract one → an all-ones mask.
        // `roundup64` wraps to 0 for `q > 2^63`, which the reference turns into
        // an all-ones mask, so subtract with wrap-around.
        self.prime = roundup64(q).wrapping_sub(1);
        Ok(())
    }

    pub(crate) fn set_window_shrink(&mut self, shrink: bool) {
        self.shrink = shrink;
        self.update_stop();
    }
}

/// Offset within `window` of the earliest element minimising
/// `(base + window[offset]) & prime`. `window` must not be empty.
#[inline(always)]
pub(crate) fn masked_argmin(window: &[u64], base: u64, prime: u64) -> usize {
    // Fast path: when the mask fits in 32 bits — as the default one does —
    // the candidate and its offset pack into a single `u64`, and the argmin
    // collapses into a branchless `min` reduction.
    if prime <= u32::MAX as u64 && window.len() <= u32::MAX as usize {
        return (packed_argmin(window, base, prime) & 0xFFFF_FFFF) as usize;
    }

    // General path for a caller-supplied mask wider than 32 bits.
    let mut best_off = 0usize;
    let mut best_val = u64::MAX;
    for (off, &h) in window.iter().enumerate() {
        let cand = base.wrapping_add(h) & prime;
        if cand < best_val {
            best_val = cand;
            best_off = off;
        }
    }
    best_off
}

/// Returns `(masked_candidate << 32) | offset` for the best position in `window`,
/// where `masked_candidate` is `(base + window[offset]) & prime`.
///
/// Packing the offset into the low bits turns the argmin into a plain `min`
/// reduction and keeps the earliest position among ties — the same choice the
/// branchy form's strict `<` makes. Callers must ensure `prime <= u32::MAX` and
/// `window.len() <= u32::MAX`, and `window` must not be empty.
///
/// AVX-512 has a 64-bit unsigned SIMD minimum (`vpminuq`), so there the plain
/// loop vectorises on its own and beats any hand-unrolling.
#[cfg(target_feature = "avx512f")]
#[inline(always)]
fn packed_argmin(window: &[u64], base: u64, prime: u64) -> u64 {
    let mut best = u64::MAX;
    for (off, &h) in window.iter().enumerate() {
        best = best.min(((base.wrapping_add(h) & prime) << 32) | off as u64);
    }
    best
}

/// See the AVX-512 variant above for the packing scheme.
///
/// Without a 64-bit SIMD minimum a single `min` chain is bound by the latency
/// of each `min` feeding the next, so four independent accumulators are used to
/// let the CPU overlap four reductions.
#[cfg(not(target_feature = "avx512f"))]
#[inline(always)]
// The lint only exists in recent clippy releases.
#[allow(unknown_lints, clippy::chunks_exact_to_as_chunks)]
fn packed_argmin(window: &[u64], base: u64, prime: u64) -> u64 {
    // `as_chunks::<4>()` (what clippy suggests here) measured slower on this
    // loop, so the runtime-length form is kept deliberately.
    let mut chunks = window.chunks_exact(4);
    let mut acc = [u64::MAX; 4];
    let mut off = 0u64;

    for c in &mut chunks {
        for (j, &h) in c.iter().enumerate() {
            acc[j] = acc[j].min(((base.wrapping_add(h) & prime) << 32) | (off + j as u64));
        }
        off += 4;
    }
    for (j, &h) in chunks.remainder().iter().enumerate() {
        acc[0] = acc[0].min(((base.wrapping_add(h) & prime) << 32) | (off + j as u64));
    }

    acc[0].min(acc[1]).min(acc[2]).min(acc[3])
}

/// Implements the public surface shared by both iterators on top of a
/// `state: State` field and an inherent `advance(&mut self) -> Option<u64>`.
macro_rules! strobe_iterator {
    ($ty:ident) => {
        impl $ty<'_> {
            /// Sets the mask used when comparing candidate strobes by
            /// `(base_hash + candidate_hash) & mask`.
            ///
            /// The name follows the reference implementations, which call this
            /// value a "prime". It is used as a bit mask: `q` is rounded up to
            /// the next power of two and decremented by one, giving `2^k - 1`,
            /// which is generally not a prime number. The default is
            /// [`DEFAULT_PRIME_NUMBER`](crate::DEFAULT_PRIME_NUMBER).
            ///
            /// # Errors
            ///
            /// [`StrobeError::PrimeNumberTooSmall`](crate::StrobeError::PrimeNumberTooSmall)
            /// if `q < 256`; the mask is left unchanged.
            pub fn set_prime(&mut self, q: u64) -> $crate::Result<()> {
                self.state.set_prime(q)
            }

            /// Enables or disables window shrinking at the sequence end.
            ///
            /// When `shrink = true`, terminal windows may be smaller than `w_max`.
            /// When `shrink = false`, iteration stops if a full window cannot be formed.
            pub fn set_window_shrink(&mut self, shrink: bool) {
                self.state.set_window_shrink(shrink)
            }

            /// Returns the index of the last returned first-strobe (m1).
            ///
            /// If no strobe has been generated yet, returns `None`.
            pub fn index(&self) -> Option<usize> {
                self.state.idx.checked_sub(1)
            }

            /// Returns the indices of the most recently generated strobes: [m1, m2, m3].
            ///
            /// The third entry is only meaningful for order 3. If no strobe has been
            /// generated yet, returns `[0, 0, 0]`.
            ///
            /// This reads state left behind by the last `next()` call, so it cannot be
            /// combined with iterator adapters such as `take` or `enumerate`. Prefer
            /// [`with_positions`](Self::with_positions), which yields the positions
            /// together with each hash.
            pub fn indexes(&self) -> [usize; 3] {
                [
                    self.index().unwrap_or(0),
                    self.state.idx2,
                    self.state.idx3,
                ]
            }

            /// Advances the iterator and returns the next strobemer with the start
            /// position of each of its strobes.
            ///
            /// Equivalent to calling `next()` and then [`indexes`](Self::indexes).
            #[inline]
            pub fn next_strobemer(&mut self) -> Option<$crate::Strobemer> {
                let hash = self.advance()?;
                Some($crate::Strobemer {
                    hash,
                    m1: self.state.idx - 1,
                    m2: self.state.idx2,
                    m3: (self.state.n == 3).then_some(self.state.idx3),
                })
            }

            /// Converts this iterator into one yielding
            /// [`Strobemer`](crate::Strobemer)s: each hash together with its
            /// strobe positions.
            ///
            /// # Example
            /// ```
            #[doc = concat!("use strobemers_rs::", stringify!($ty), ";")]
            ///
            #[doc = concat!("let it = ", stringify!($ty), "::new(b\"ACGATCTGGTACCTAG\", 3, 3, 3, 5)?;")]
            /// for (i, s) in it.with_positions().take(5).enumerate() {
            ///     println!("{i}: m1={} m2={} m3={:?} hash={}", s.m1, s.m2, s.m3, s.hash);
            /// }
            /// # Ok::<(), strobemers_rs::StrobeError>(())
            /// ```
            pub fn with_positions(self) -> $crate::WithPositions<Self> {
                $crate::WithPositions::new(self)
            }
        }

        impl Iterator for $ty<'_> {
            type Item = u64;

            /// Returns the next strobemer hash value.
            #[inline]
            fn next(&mut self) -> Option<u64> {
                self.advance()
            }

            /// Exact remaining length, so `collect()` allocates exactly once.
            #[inline]
            fn size_hint(&self) -> (usize, Option<usize>) {
                let n = self.state.remaining();
                (n, Some(n))
            }
        }

        /// Exact via [`Iterator::size_hint`].
        impl ExactSizeIterator for $ty<'_> {}

        /// Once exhausted, `next()` keeps returning `None`: nothing advances
        /// past the end.
        impl core::iter::FusedIterator for $ty<'_> {}
    };
}
