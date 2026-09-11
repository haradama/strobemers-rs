use crate::{Result, StrobeError};
use nthash_rs::kmer::NtHash;

pub trait KmerHasher: Send + Sync + 'static {
    fn hash_all(&self, seq: &[u8], k: usize) -> Result<Vec<u64>>;
}

pub struct NtHash64;
impl Default for NtHash64 {
    fn default() -> Self {
        Self
    }
}

impl KmerHasher for NtHash64 {
    fn hash_all(&self, seq: &[u8], k: usize) -> Result<Vec<u64>> {
        if !(1..=64).contains(&k) {
            return Err(StrobeError::StrobeLengthTooSmall);
        }
        if seq.len() < k {
            return Err(StrobeError::SequenceTooShort);
        }

        // `NtHashIter` yields `(pos, Vec<u64>)`, so iterating it heap-allocates a
        // one-element vector per k-mer. Drive the roller directly instead and
        // read the canonical hash out of its reusable buffer.
        let mut roller = NtHash::new(seq, k as u16, 1, 0).map_err(StrobeError::from)?;

        let mut out = Vec::with_capacity(seq.len() - k + 1);
        while roller.roll() {
            out.push(roller.hashes()[0]);
        }
        Ok(out)
    }
}

/// Verifies that a [`KmerHasher`] returned exactly one hash per `k`-mer.
///
/// Both strobemer iterators treat the index into the hash vector as the position
/// of the corresponding `k`-mer in the sequence. `nthash-rs` **skips** windows
/// containing `N` (or any non-ACGT byte), which silently breaks that mapping and
/// would lead to out-of-range strobe positions, so the length is checked here.
/// The Go reference performs the same check and reports `ErrIncompleteHashValues`.
///
/// # Parameters
/// - `hashes`: the hash values returned by the hasher
/// - `seq_len`: length of the hashed sequence
/// - `k`: the `k`-mer length
///
/// # Returns
/// `Ok(())` when `hashes.len() == seq_len - k + 1`, otherwise
/// [`StrobeError::IncompleteHashValues`].
pub(crate) fn check_hash_count(hashes: &[u64], seq_len: usize, k: usize) -> Result<()> {
    // `validate_params!` guarantees `seq_len >= k`, so this cannot underflow.
    if hashes.len() != seq_len - k + 1 {
        return Err(StrobeError::IncompleteHashValues);
    }
    Ok(())
}

/// Online sliding-window minimum over a fixed window width.
///
/// Positions are absorbed strictly in order. Alongside the monotonic deque the
/// window keeps the minima of the last `lag + 1` absorbed positions, so a caller
/// that needs two minima a fixed distance apart — as order-3 `MinStrobes` does,
/// at `idx + w_max` and `idx + 2 * w_max` — gets both from a single pass instead
/// of running two deques over the same data.
///
/// This is the same computation [`compute_min_hashes`] runs in bulk, but it
/// keeps `O(w + lag)` entries rather than one pair per k-mer.
#[derive(Debug, Clone)]
pub(crate) struct MinWindow {
    /// Positions of the candidates, increasing from `head`.
    idx_q: Box<[usize]>,
    /// Their hash values, also increasing from `head` (monotonic deque).
    val_q: Box<[u64]>,
    head: usize,
    len: usize,
    w: usize,
    /// Next position the window expects to absorb.
    next_pos: usize,
    /// Ring of `(position, value)` minima, indexed by `position & hist_mask`.
    /// Its length is a power of two so the wrap is a mask rather than a `div`.
    hist: Box<[(usize, u64)]>,
    hist_mask: usize,
}

impl MinWindow {
    /// Creates an empty window of width `w` (must be ≥ 1) that can still report
    /// the minimum `lag` positions behind the one most recently absorbed.
    pub(crate) fn new(w: usize, lag: usize) -> Self {
        assert!(w >= 1, "window size must be ≥ 1");
        let hist_len = (lag + 1).next_power_of_two();
        Self {
            idx_q: vec![0usize; w].into_boxed_slice(),
            val_q: vec![0u64; w].into_boxed_slice(),
            head: 0,
            len: 0,
            w,
            next_pos: 0,
            hist: vec![(0usize, u64::MAX); hist_len].into_boxed_slice(),
            hist_mask: hist_len - 1,
        }
    }

    /// `head + len <= 2 * w` always holds, so wrapping the deque never needs a
    /// division: one conditional subtraction is enough. `%` compiles to a real
    /// `div` whenever `w` is not a power of two, and it sits in the inner loop.
    #[inline(always)]
    fn wrap(&self, x: usize) -> usize {
        if x >= self.w { x - self.w } else { x }
    }

    /// Absorbs the value at position `i`, dropping anything that has fallen out
    /// of the window or can no longer be the minimum.
    #[inline]
    fn push(&mut self, i: usize, h: u64) {
        let window_start = i.saturating_sub(self.w - 1);
        while self.len > 0 && self.idx_q[self.head] < window_start {
            self.head = self.wrap(self.head + 1);
            self.len -= 1;
        }

        // Strictly-greater comparison: equal values are kept, so the front of the
        // deque holds the *earliest* position among ties. The Go reference keeps
        // a sorted buffer that inserts after equal elements, giving the same
        // choice; using `>=` here would report the latest position instead.
        while self.len > 0 && self.val_q[self.wrap(self.head + self.len - 1)] > h {
            self.len -= 1;
        }

        let tail = self.wrap(self.head + self.len);
        self.idx_q[tail] = i;
        self.val_q[tail] = h;
        self.len += 1;
    }

    /// `(position, value)` of the minimum over the currently absorbed window.
    #[inline]
    fn min(&self) -> (usize, u64) {
        (self.idx_q[self.head], self.val_q[self.head])
    }

    /// The minimum of the window that ended at `p`.
    ///
    /// Valid for the `lag + 1` most recently absorbed positions.
    #[inline]
    pub(crate) fn at(&self, p: usize) -> (usize, u64) {
        self.hist[p & self.hist_mask]
    }

    /// Absorbs every position up to and including `p`, then reports the minimum
    /// of the window ending at `p`.
    ///
    /// `p` must not go backwards by more than `lag` between calls, and
    /// `hashes[..=p]` must be valid.
    #[inline]
    pub(crate) fn advance_to(&mut self, hashes: &[u64], p: usize) -> (usize, u64) {
        // With no lag the caller only ever asks for the position just absorbed,
        // so the history is pure overhead — one extra 16-byte store per k-mer.
        if self.hist_mask == 0 {
            while self.next_pos <= p {
                self.push(self.next_pos, hashes[self.next_pos]);
                self.next_pos += 1;
            }
            return self.min();
        }

        while self.next_pos <= p {
            let i = self.next_pos;
            self.push(i, hashes[i]);
            self.hist[i & self.hist_mask] = self.min();
            self.next_pos += 1;
        }
        self.at(p)
    }
}

/// For a sliding window of width `w` over the given slice of hash values,
/// computes the index and value of the minimum hash in each window.
///
/// Materialises both series. [`MinStrobes`](crate::MinStrobes) uses this for
/// short sequences and [`MinWindow`] for long ones; the deque below is
/// deliberately written out rather than driven through [`MinWindow`], because in
/// a bulk pass the hot state has to live in locals — behind `&mut self` it costs
/// about 18% here. `minima_strategies_agree` in `minstrobes` pins the two
/// against each other.
///
/// # Parameters
/// - `hashes`: slice of u64 hash values to slide over
/// - `w`: the size of the sliding window (must be ≥ 1)
///
/// # Returns
/// A tuple `(locs, mins)` where:
/// - `locs[i]` is the index of the minimum hash in the window ending at position `i`
/// - `mins[i]` is the minimum hash value in that same window
///
/// Only valid when `i ≥ w - 1`; for indices `< w - 1`, the values in `locs` and `mins` remain default (0 and `u64::MAX`).
///
pub fn compute_min_hashes(hashes: &[u64], w: usize) -> (Vec<usize>, Vec<u64>) {
    assert!(w >= 1, "window size must be ≥ 1");
    let n = hashes.len();

    if w == 1 {
        return ((0..n).collect(), hashes.to_vec());
    }

    let mut locs = vec![0usize; n];
    let mut mins = vec![u64::MAX; n];

    let mut idx_q = vec![0usize; w];
    let mut val_q = vec![0u64; w];
    let mut head = 0usize;
    let mut len = 0usize;

    // `head + len <= 2 * w` always holds, so wrapping the ring buffer never needs
    // a division: one conditional subtraction is enough. `%` compiles to a real
    // `div` whenever `w` is not a power of two, and it sits in the inner loop.
    #[inline(always)]
    fn wrap(x: usize, w: usize) -> usize {
        if x >= w { x - w } else { x }
    }

    for (i, &h) in hashes.iter().enumerate() {
        let window_start = i.saturating_sub(w - 1);
        while len > 0 && idx_q[head] < window_start {
            head = wrap(head + 1, w);
            len -= 1;
        }

        // Strictly-greater comparison: equal values are kept, so the front of the
        // deque holds the *earliest* position among ties. The Go reference keeps
        // a sorted buffer that inserts after equal elements, giving the same
        // choice; using `>=` here would report the latest position instead.
        while len > 0 && val_q[wrap(head + len - 1, w)] > h {
            len -= 1;
        }

        let tail = wrap(head + len, w);
        idx_q[tail] = i;
        val_q[tail] = h;
        len += 1;

        if i >= w - 1 {
            locs[i] = idx_q[head];
            mins[i] = val_q[head];
        }
    }
    (locs, mins)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slide_min_window_ties_pick_earliest() {
        // Values chosen so that every window width exercises a tie.
        // Expected locations are the output of the Go reference implementation
        // (shenwei356/strobemers, computeMinHashes).
        let v = [5u64, 3, 3, 1, 1, 1, 4, 2, 2, 9, 9, 0];
        for (w, want) in [
            (2usize, &[1usize, 1, 3, 3, 4, 5, 7, 7, 8, 9, 11][..]),
            (3, &[1, 3, 3, 3, 4, 5, 7, 7, 8, 11][..]),
            (4, &[3, 3, 3, 3, 4, 5, 7, 7, 11][..]),
        ] {
            let (locs, mins) = compute_min_hashes(&v, w);
            assert_eq!(&locs[w - 1..], want, "locations for w={w}");
            // The reported location must actually carry the window minimum.
            for (i, &loc) in locs.iter().enumerate().skip(w - 1) {
                assert_eq!(v[loc], mins[i], "value at reported loc for w={w}, i={i}");
            }
        }
    }

    /// [`MinWindow`] and [`compute_min_hashes`] hold the tie-breaking rule
    /// separately, so pin the online form against the same Go-derived
    /// expectations — and against the bulk form, position by position.
    #[test]
    fn min_window_matches_bulk_on_ties() {
        let v = [5u64, 3, 3, 1, 1, 1, 4, 2, 2, 9, 9, 0];
        for (w, want) in [
            (2usize, &[1usize, 1, 3, 3, 4, 5, 7, 7, 8, 9, 11][..]),
            (3, &[1, 3, 3, 3, 4, 5, 7, 7, 8, 11][..]),
            (4, &[3, 3, 3, 3, 4, 5, 7, 7, 11][..]),
        ] {
            let (bulk_locs, bulk_mins) = compute_min_hashes(&v, w);
            // `lag = 2` so the history ring is exercised rather than bypassed.
            let mut win = MinWindow::new(w, 2);
            for (offset, &want_loc) in want.iter().enumerate() {
                let p = w - 1 + offset;
                assert_eq!(
                    win.advance_to(&v, p),
                    (want_loc, v[want_loc]),
                    "w={w} p={p}"
                );
                assert_eq!(win.at(p), (bulk_locs[p], bulk_mins[p]), "w={w} p={p}");
                // Positions still inside the lag must remain readable.
                // `compute_min_hashes` only fills indices from `w - 1` on, so
                // stay inside the range where both forms are defined.
                if p > w {
                    assert_eq!(win.at(p - 2), (bulk_locs[p - 2], bulk_mins[p - 2]));
                }
            }
        }
    }

    #[test]
    fn slide_min_window_three() {
        // Test vector: [5, 3, 6, 1, 4]
        // Windows of size 3:
        //  - [5, 3, 6] → min=3 at index=1
        //  - [3, 6, 1] → min=1 at index=3
        //  - [6, 1, 4] → min=1 at index=3
        let v = [5, 3, 6, 1, 4];
        let (locs, mins) = compute_min_hashes(&v, 3);
        // We only care about positions ≥ w - 1 (i.e., indices 2, 3, 4)
        assert_eq!(&mins[2..], &[3, 1, 1]);
        assert_eq!(&locs[2..], &[1, 3, 3]);
    }
}
