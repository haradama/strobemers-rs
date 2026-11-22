use crate::{Result, StrobeError, util::is_supported_base};
use nthash_rs::kmer::NtHash;

/// Source of k-mer hashes for [`MinStrobes`](crate::MinStrobes) and
/// [`RandStrobes`](crate::RandStrobes).
///
/// Implement this to replace the default ntHash-based [`NtHash64`], and pass
/// the implementation to `with_hasher`.
pub trait KmerHasher {
    /// Hashes every `k`-mer of `seq`, in order.
    ///
    /// Must return exactly `seq.len() - k + 1` values, where element `i` is the
    /// hash of `seq[i..i + k]`; the iterators use the index as the k-mer's
    /// position, and reject any other length with
    /// [`StrobeError::IncompleteHashValues`].
    ///
    /// The constructors only call this with a non-empty ASCII `seq` and
    /// `1 <= k <= min(64, seq.len())`.
    fn hash_all(&self, seq: &[u8], k: usize) -> Result<Vec<u64>>;
}

/// The default [`KmerHasher`]: canonical 64-bit ntHash, via `nthash-rs`.
///
/// Accepts `A`, `C`, `G` and `T` in either case. Any other base — including
/// `U`, so RNA must be converted first — makes
/// [`hash_all`](KmerHasher::hash_all) fail with
/// [`StrobeError::UnsupportedBase`]; see [`acgt_runs`](crate::acgt_runs).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NtHash64;

/// The error for the first base of `seq` that ntHash cannot hash, if any.
pub(crate) fn unsupported_base(seq: &[u8]) -> Option<StrobeError> {
    // Checking a block at a time without an early exit lets the common,
    // all-valid case vectorise; the position is only looked up on failure.
    if seq
        .chunks(4096)
        .all(|block| block.iter().fold(true, |ok, &b| ok & is_supported_base(b)))
    {
        return None;
    }
    seq.iter()
        .position(|&b| !is_supported_base(b))
        .map(|position| StrobeError::UnsupportedBase {
            position,
            byte: seq[position],
        })
}

/// Appends the canonical ntHash of every `k`-mer of `seq` to `out`.
///
/// `seq` must hold at least one `k`-mer and `k` must be in `1..=64`.
fn nthash_into(seq: &[u8], k: usize, out: &mut Vec<u64>) -> Result<()> {
    let expected = seq.len() - k + 1;
    let before = out.len();
    out.reserve(expected);

    // `NtHashIter` yields `(pos, Vec<u64>)`, so iterating it heap-allocates a
    // one-element vector per k-mer. Drive the roller directly instead and
    // read the canonical hash out of its reusable buffer.
    let mut roller = match NtHash::new(seq, k as u16, 1, 0) {
        Ok(roller) => roller,
        Err(e) => return Err(unsupported_base(seq).unwrap_or_else(|| e.into())),
    };
    while roller.roll() {
        out.push(roller.hashes()[0]);
    }

    // `nthash-rs` skips every window containing a base it has no seed for,
    // which leaves fewer hashes than k-mers. Name the culprit rather than
    // letting the caller see a bare length mismatch.
    if out.len() - before != expected {
        return Err(unsupported_base(seq).unwrap_or(StrobeError::IncompleteHashValues));
    }
    Ok(())
}

impl KmerHasher for NtHash64 {
    fn hash_all(&self, seq: &[u8], k: usize) -> Result<Vec<u64>> {
        if !(1..=64).contains(&k) {
            return Err(StrobeError::InvalidStrobeLength);
        }
        if seq.len() < k {
            return Err(StrobeError::SequenceTooShort);
        }
        let mut out = Vec::new();
        nthash_into(seq, k, &mut out)?;
        Ok(out)
    }
}

/// How many k-mers to hash ahead of the iterator at a time.
///
/// Large enough that sliding the buffer forward is rare and the hashing loop
/// runs in long stretches; small enough that the buffer stays cache-resident.
/// The `scaling` benchmark group shows the effect, so re-run it before changing
/// this value.
pub(crate) const CHUNK_KMERS: usize = 1 << 14;

/// The k-mer hashes an iterator can currently see.
///
/// The iterators only ever look a bounded distance ahead of their position, so
/// the hashes do not all have to exist at once. This holds a contiguous run of
/// them, `buf[i]` being the hash of the k-mer at `base + i`, and — when built
/// with [`HashBuf::streaming`] — hashes the next stretch of the sequence as the
/// iterator approaches the end of the run.
///
/// Keeping the run contiguous, rather than a ring, lets the callers slice their
/// search windows directly; refilling in bulk keeps the hashing loop tight.
#[derive(Debug, Clone)]
pub(crate) struct HashBuf<'a> {
    buf: Vec<u64>,
    /// Position of the k-mer `buf[0]` belongs to.
    base: usize,
    /// Most k-mers ever held at once.
    cap: usize,
    /// Number of k-mers in the whole sequence.
    total: usize,
    /// The sequence still to be hashed from; unused once `buf` reaches the end.
    seq: &'a [u8],
    k: usize,
}

impl HashBuf<'static> {
    /// Wraps hashes that already cover the whole sequence.
    pub(crate) fn complete(hashes: Vec<u64>) -> Self {
        Self {
            base: 0,
            cap: hashes.len(),
            total: hashes.len(),
            buf: hashes,
            seq: &[],
            k: 0,
        }
    }
}

impl<'a> HashBuf<'a> {
    /// Hashes `seq` with ntHash on demand.
    ///
    /// `lookahead` is the furthest a caller reads past the position it passes
    /// to [`ensure`](Self::ensure); `chunk` is how many further k-mers to hash
    /// per refill. `seq` must be at least `k` long and contain only bases
    /// ntHash supports, which is what makes the later refills infallible.
    pub(crate) fn streaming(
        seq: &'a [u8],
        k: usize,
        lookahead: usize,
        chunk: usize,
    ) -> Result<Self> {
        let total = seq.len() - k + 1;
        let cap = total.min(lookahead.saturating_add(chunk).saturating_add(1));
        let mut buf = Vec::new();
        nthash_into(&seq[..cap + k - 1], k, &mut buf)?;
        Ok(Self {
            buf,
            base: 0,
            cap,
            total,
            seq,
            k,
        })
    }

    /// Whether every k-mer of the sequence is in the buffer.
    pub(crate) fn is_complete(&self) -> bool {
        self.buf.len() == self.total
    }

    /// The buffered hashes and the position of the first one.
    #[inline(always)]
    pub(crate) fn view(&self) -> (&[u64], usize) {
        (&self.buf, self.base)
    }

    /// Makes positions `idx..=hi` available, discarding everything before
    /// `idx`. Requires `hi - idx <= lookahead`, and `idx` must not decrease.
    #[inline(always)]
    pub(crate) fn ensure(&mut self, idx: usize, hi: usize) {
        if hi - self.base >= self.buf.len() {
            self.refill(idx);
        }
    }

    /// Moves the run forward so that it starts at `idx`.
    #[cold]
    #[inline(never)]
    fn refill(&mut self, idx: usize) {
        let keep_from = idx - self.base;
        self.buf.copy_within(keep_from.., 0);
        self.buf.truncate(self.buf.len() - keep_from);
        self.base = idx;

        let from = idx + self.buf.len();
        let to = self.total.min(idx + self.cap);
        if from < to {
            nthash_into(&self.seq[from..to + self.k - 1], self.k, &mut self.buf)
                .expect("the sequence was validated when the iterator was built");
        }
    }
}

/// Verifies that a [`KmerHasher`] returned exactly one hash per `k`-mer.
///
/// Both strobemer iterators treat the index into the hash vector as the position
/// of the corresponding `k`-mer in the sequence. `nthash-rs` **skips** windows
/// containing `N` (or any non-ACGT byte) — [`NtHash64`] reports that itself — and
/// a custom hasher may do the same, which silently breaks that mapping and would
/// lead to out-of-range strobe positions, so the length is checked here.
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
    // `validate_params` guarantees `seq_len >= k`, so this cannot underflow.
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
    /// `buf[i]` is the hash at position `base + i`. `p` must not go backwards
    /// by more than `lag` between calls, and every position not yet absorbed up
    /// to `p` must lie inside `buf`.
    #[inline]
    pub(crate) fn advance_to(&mut self, buf: &[u64], base: usize, p: usize) -> (usize, u64) {
        // With no lag the caller only ever asks for the position just absorbed,
        // so the history is pure overhead — one extra store per k-mer.
        if self.hist_mask == 0 {
            while self.next_pos <= p {
                self.push(self.next_pos, buf[self.next_pos - base]);
                self.next_pos += 1;
            }
            return self.min();
        }

        while self.next_pos <= p {
            let i = self.next_pos;
            self.push(i, buf[i - base]);
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
/// a bulk pass the hot state has to live in locals — behind `&mut self` it is
/// measurably slower. `minima_strategies_agree` in `minstrobes` pins the two
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
pub(crate) fn compute_min_hashes(hashes: &[u64], w: usize) -> (Vec<usize>, Vec<u64>) {
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
                    win.advance_to(&v, 0, p),
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
