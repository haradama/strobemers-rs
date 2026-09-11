use crate::{
    Result, StrobeError,
    constants::DEFAULT_PRIME_NUMBER,
    hashes::{KmerHasher, MinWindow, NtHash64, check_hash_count, compute_min_hashes},
    util::{last_start_index, roundup64},
};

/// Largest sliding-window-minima table this module will precompute, in bytes.
///
/// The table is a `usize` and a `u64` per k-mer, so 16 bytes per base on top of
/// the hash array itself. Below this the bulk loop that fills it is the faster
/// way to get the minima; above it the extra memory traffic costs more than
/// recomputing them as the iterator advances.
///
/// The threshold sits where the two were measured to cross (Xeon Platinum 8480+,
/// k=31, w=20..50): at 500 kbp the table still wins slightly, at 1 Mbp the online
/// window is about 20% faster, and it stays that way as the sequence grows.
const MAX_MINIMA_TABLE_BYTES: usize = 8 * 1024 * 1024;

/// Bytes of table per k-mer: one `usize` location plus one `u64` value.
const MINIMA_BYTES_PER_KMER: usize = size_of::<usize>() + size_of::<u64>();

/// Where the sliding-window minima come from.
///
/// Both variants answer the same question — the minimum of the window ending at
/// a given position — and produce identical results; they differ only in whether
/// the answers are materialised up front.
#[derive(Debug, Clone)]
enum Minima {
    /// Precomputed for both series. Fastest while the arrays stay cache-resident.
    Table { loc: Vec<usize>, val: Vec<u64> },
    /// Computed as the iterator advances, keeping `O(w_max)` entries instead of
    /// one pair per k-mer. About 20% faster once the table would no longer fit
    /// in cache, and it removes 16 bytes per base of peak memory.
    Online(MinWindow),
}

impl Minima {
    /// Chooses a strategy for a sequence with `hashes.len()` k-mers.
    ///
    /// `lag` is how far behind the most recent position the caller may still
    /// read; order 3 needs `w_max`, order 2 needs none.
    fn new(hashes: &[u64], window: usize, lag: usize) -> Self {
        if hashes.len().saturating_mul(MINIMA_BYTES_PER_KMER) <= MAX_MINIMA_TABLE_BYTES {
            let (loc, val) = compute_min_hashes(hashes, window);
            Self::Table { loc, val }
        } else {
            Self::Online(MinWindow::new(window, lag))
        }
    }

    /// `(position, value)` of the minimum over the window ending at `p`.
    ///
    /// For the online variant `p` must not go backwards by more than `lag`
    /// between calls, which both call sites respect.
    #[inline]
    fn at(&mut self, hashes: &[u64], p: usize) -> (usize, u64) {
        match self {
            Self::Table { loc, val } => (loc[p], val[p]),
            Self::Online(win) => win.advance_to(hashes, p),
        }
    }
}

/// Iterator for generating MinStrobes of order 2 or 3 from a DNA/RNA sequence.
///
/// A MinStrobe is a concatenation of k-mers selected based on minimum hash
/// values within sliding windows. This struct precomputes k-mer hashes and
/// window minima to efficiently produce strobemer hash values.
///
#[derive(Debug, Clone)]
pub struct MinStrobes {
    // Parameters controlling strobemer generation
    n: u8,        // Order of strobemer: 2 or 3
    w_min: usize, // Minimum window offset
    w_max: usize, // Maximum window offset

    // Precomputed data
    hashes: Vec<u64>, // Hash values for each k-mer in the sequence

    minima: Minima, // Sliding-window minima, precomputed or computed online

    // Iteration state
    idx: usize, // Current index of the first k-mer (m1)
    // Last index at which a complete strobemer can start, or `None` when the
    // sequence cannot hold `n` strobes at all. The Go reference stores this in a
    // signed `int` and lets it go negative; `usize` would wrap, so the
    // "no strobemer fits" case is modelled explicitly.
    end_idx: Option<usize>,
    end_hash: usize, // Last index in `hashes` (i.e., sequence length minus k)

    // Strobe indices for current item
    idx2: usize, // Index of second k-mer (m2)
    idx3: usize, // Index of third k-mer (m3) if order = 3

    // Prime number and shrink-window flag
    prime: u64,   // Used for combining hash values in order 3
    shrink: bool, // Whether to shrink windows near sequence end

    // Working registers for hash values
    h1: u64, // Hash of first k-mer (m1)
    h2: u64, // Combined hash after selecting m2
    h3: u64, // Combined hash after selecting m3 (order 3 only)
}

impl MinStrobes {
    /// Constructs a new [`MinStrobes`] iterator using the default hash function (`NtHash64`).
    ///
    /// Internally delegates to [`MinStrobes::with_hasher`] with the standard ntHash implementation.
    /// This function performs all necessary validation and preprocessing to enable
    /// efficient strobemer generation.
    ///
    /// # Arguments
    ///
    /// * `seq` – Input nucleotide sequence as a byte slice (DNA/RNA, ASCII only).
    /// * `n` – Order of the strobemer (must be 2 or 3).
    /// * `k` – Length of each strobe segment (k-mer); must be in `[1, 64]`.
    /// * `w_min` – Minimum offset (in bases) between strobes.
    /// * `w_max` – Maximum offset (inclusive); must satisfy `w_min ≤ w_max`.
    ///
    /// # Returns
    ///
    /// * `Ok(MinStrobes)` on success.
    /// * `Err(StrobeError)` if parameters are invalid or the sequence is too short.
    ///
    /// # Example
    /// ```
    /// use strobemers_rs::MinStrobes;
    /// let ms = MinStrobes::new(b"ACGTACGTACGT", 2, 3, 1, 4).unwrap();
    /// ```
    pub fn new(seq: &[u8], n: u8, k: usize, w_min: usize, w_max: usize) -> Result<Self> {
        Self::with_hasher(seq, n, k, w_min, w_max, &NtHash64)
    }

    /// Constructs a new [`MinStrobes`] iterator with a user-defined hash function.
    ///
    /// This method accepts any implementation of the [`KmerHasher`] trait,
    /// allowing for custom k-mer hashing strategies, such as XOR-based, cryptographic,
    /// or rolling hashes optimized for performance or reproducibility.
    ///
    /// Precomputes:
    /// - `k`-mer hashes from the sequence using `hasher`
    /// - Sliding window minima (position and value) for efficient downstream selection
    ///
    /// # Arguments
    ///
    /// * `seq` – Input DNA/RNA sequence as bytes (e.g., `b"ACGT..."`).
    /// * `n` – Strobemer order (only 2 or 3 are supported).
    /// * `k` – Length of each strobe (k-mer), must be `1..=64`.
    /// * `w_min` – Minimum window offset after the first strobe.
    /// * `w_max` – Maximum window offset after the first strobe.
    /// * `hasher` – A reference to a type implementing the [`KmerHasher`] trait.
    ///
    /// # Returns
    ///
    /// * `Ok(MinStrobes)` – Ready-to-use iterator for strobemers.
    /// * `Err(StrobeError)` – On invalid parameters or hash failure.
    ///
    /// # Example
    /// ```
    /// use strobemers_rs::{MinStrobes, KmerHasher};
    ///
    /// struct DummyHasher;
    /// impl KmerHasher for DummyHasher {
    ///     fn hash_all(&self, seq: &[u8], k: usize) -> strobemers_rs::Result<Vec<u64>> {
    ///         Ok(seq.windows(k).map(|w| w.iter().map(|b| *b as u64).sum()).collect())
    ///     }
    /// }
    ///
    /// let hasher = DummyHasher;
    /// let ms = MinStrobes::with_hasher(b"ACGTACGT", 2, 3, 1, 4, &hasher).unwrap();
    /// ```
    pub fn with_hasher<H>(
        seq: &[u8],
        n: u8,
        k: usize,
        w_min: usize,
        w_max: usize,
        hasher: &H,
    ) -> Result<Self>
    where
        H: KmerHasher,
    {
        // Check all preconditions
        validate_params!(seq, n, k, w_min, w_max);

        // Compute k-mer hash values via user-supplied hasher
        let hashes = hasher.hash_all(seq, k)?;
        // The iterator indexes `hashes` by sequence position, so the hasher must
        // have produced exactly one value per k-mer.
        check_hash_count(&hashes, seq.len(), k)?;

        // Window width shared by both strobe selections. Order 3 reads the
        // minima at two positions `w_max` apart; order 2 only ever reads the
        // position it just reached.
        let window = w_max - w_min + 1;
        let minima = Minima::new(&hashes, window, if n == 3 { w_max } else { 0 });

        // Define range bounds for m1 (starting point of each strobemer)
        let seq_len = seq.len();
        let end_hash = seq_len - k;
        // `None` when the sequence is shorter than `n * k`, in which case the
        // iterator yields nothing.
        let end_idx = seq_len.checked_sub(k + (n as usize - 1) * k);

        Ok(Self {
            n,
            w_min,
            w_max,
            hashes,
            minima,
            idx: 0,
            end_hash,
            end_idx,
            idx2: 0,
            idx3: 0,
            prime: DEFAULT_PRIME_NUMBER,
            shrink: true,
            h1: 0,
            h2: 0,
            h3: 0,
        })
    }

    /// Sets a new prime number for combining hash values in order-3 strobes.
    ///
    /// The provided `q` must be at least 256. Internally, the value is rounded up
    /// to the next power of two and then decremented by one to form a Mersenne prime.
    ///
    /// # Arguments
    ///
    /// * `q` – Candidate prime (will be rounded to Mersenne form: `2^k - 1`).
    ///
    /// # Returns
    ///
    /// * `Ok(())` – If `q` ≥ 256, updates `self.prime`.
    /// * `Err(StrobeError::PrimeNumberTooSmall)` – If `q` < 256.
    pub fn set_prime(&mut self, q: u64) -> Result<()> {
        if q < 256 {
            return Err(StrobeError::PrimeNumberTooSmall);
        }
        // Round up to next power of two, subtract one → Mersenne prime form.
        // `roundup64` wraps to 0 for `q > 2^63`, which the reference turns into
        // an all-ones mask, so subtract with wrap-around.
        self.prime = roundup64(q).wrapping_sub(1);
        Ok(())
    }

    /// Enables or disables window shrinking at the sequence end.
    ///
    /// When `shrink = true`, terminal windows may be smaller than `w_max`.
    /// When `shrink = false`, iteration stops if a full window cannot be formed.
    pub fn set_window_shrink(&mut self, s: bool) {
        self.shrink = s;
    }

    /// Returns the index of the last returned first-strobe (m1).
    ///
    /// If no strobe has been generated yet, returns `None`.
    pub fn index(&self) -> Option<usize> {
        self.idx.checked_sub(1)
    }

    /// Returns the indices of the most recently generated strobes: [m1, m2, (m3)].
    ///
    /// If no strobe has been generated yet, returns `[0, 0, 0]`.
    pub fn indexes(&self) -> [usize; 3] {
        [self.index().unwrap_or(0), self.idx2, self.idx3]
    }

    /// Number of strobemers this iterator can still produce. Exact — see
    /// [`last_start_index`].
    #[inline]
    fn remaining(&self) -> usize {
        let (Some(end_idx), Some(limit)) = (
            self.end_idx,
            last_start_index(self.n, self.shrink, self.w_min, self.w_max, self.end_hash),
        ) else {
            return 0;
        };
        (limit.min(end_idx) + 1).saturating_sub(self.idx)
    }

    /// Computes the next hash value for an order-2 MinStrobe.
    fn next_order2(&mut self) -> Option<u64> {
        // Stop if no more valid starting positions for m1
        if self.idx > self.end_idx? {
            return None;
        }

        // Define the search window range for m2
        let w_start = self.idx + self.w_min;
        let mut w_end = self.idx + self.w_max;

        // Hash of the first k-mer (m1)
        self.h1 = self.hashes[self.idx];

        // If window extends past last hash index, adjust or stop
        if w_end > self.end_hash {
            if !self.shrink {
                return None;
            }
            w_end = self.end_hash;
        }

        // If the full window fits, take the running minimum. Once the window
        // starts being clamped it stays clamped, so `win_m2` is never needed
        // again and leaving it behind is safe.
        if w_end == self.idx + self.w_max {
            let (loc, val) = self.minima.at(&self.hashes, w_end);
            self.idx2 = loc;
            // Combine h1 and the window minimum
            self.h2 = (self.h1 >> 1) + val / 3;
        } else {
            // Partial window: manually scan to find minimum.
            // An empty window (`w_start > w_end`, possible once `w_min > (n-1) * k`)
            // ends iteration; `w_start` only grows, so no later index can succeed.
            // The Go reference instead silently reuses the previous `idx2`.
            if w_start > w_end {
                return None;
            }
            let (mut best_hash, mut best_pos) = (u64::MAX, w_start);
            for pos in w_start..=w_end {
                let cand = self.hashes[pos];
                if cand < best_hash {
                    best_hash = cand;
                    best_pos = pos;
                }
            }
            self.idx2 = best_pos;
            self.h2 = self.h1 / 2 + best_hash / 3;
        }

        // Advance to next starting index for m1
        self.idx += 1;
        Some(self.h2)
    }

    /// Computes the next hash value for an order-3 MinStrobe.
    ///
    /// # Returns
    /// - `Some(u64)` – Combined hash value of m1, m2, and m3, if available.
    /// - `None` – When no further strobes can be formed.
    ///
    fn next_order3(&mut self) -> Option<u64> {
        // Stop if no more valid starting positions for m1
        if self.idx > self.end_idx? {
            return None;
        }

        // Window range for selecting m2
        let w_end = self.idx + self.w_max;
        // Window range for selecting m3 (after m2 block)
        let w2_start = self.idx + self.w_max + self.w_min;
        let mut w2_end = self.idx + (self.w_max << 1);

        // If there's no room for a third k-mer, stop
        if w2_start > self.end_hash {
            return None;
        }
        // If second window extends past end, adjust or stop
        if w2_end > self.end_hash {
            if !self.shrink {
                return None;
            }
            w2_end = self.end_hash;
        }

        // Compute m1 (first k-mer)
        self.h1 = self.hashes[self.idx];
        // Select m2 from the running minimum at the window end
        let (loc2, val2) = self.minima.at(&self.hashes, w_end);
        self.idx2 = loc2;
        self.h2 = self.h1 / 3 + (val2 >> 2);

        // Select m3
        if w2_end == self.idx + (self.w_max << 1) {
            // Full second window fits: take the running minimum
            let (loc3, val3) = self.minima.at(&self.hashes, w2_end);
            self.idx3 = loc3;
            self.h3 = self.h2 + val3 / 5;
        } else {
            // Partial second window near the end: manual scan
            let (mut best_hash, mut best_pos) = (u64::MAX, w2_start);
            for pos in w2_start..=w2_end {
                // Combine current h2 with candidate hash, then mask with prime.
                // `h2` is unreduced here, so the sum can exceed `u64::MAX`; the Go
                // reference relies on uint64 wrap-around, hence `wrapping_add`.
                let cand = self.h2.wrapping_add(self.hashes[pos]) & self.prime;
                if cand < best_hash {
                    best_hash = cand;
                    best_pos = pos;
                }
            }
            self.idx3 = best_pos;
            self.h3 = self.h2 + self.hashes[self.idx3] / 5;
        }

        // Advance to next starting index for m1
        self.idx += 1;
        Some(self.h3)
    }
}

impl Iterator for MinStrobes {
    type Item = u64;

    /// Exact remaining length, so `collect()` allocates exactly once.
    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.remaining();
        (n, Some(n))
    }

    /// Advances the iterator, returning the next strobemer hash value.
    ///
    /// Dispatches to `next_order2` or `next_order3` based on `self.n`.
    /// If `n` is not 2 or 3, returns `None`.
    fn next(&mut self) -> Option<Self::Item> {
        match self.n {
            2 => self.next_order2(),
            3 => self.next_order3(),
            _ => None, // Should not occur due to prior validation
        }
    }
}

/// Exact via [`Iterator::size_hint`].
impl ExactSizeIterator for MinStrobes {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn order2_basic() {
        // Basic smoke test: sequence "ACGTACGTACGT", order=2, k=3, w_min=1, w_max=4
        let mut ms = MinStrobes::new("ACGTACGTACGT".as_bytes(), 2, 3, 1, 4).unwrap();
        // Expect at least one strobemer
        assert!(ms.next().is_some());
    }

    /// Deterministic pseudo-random DNA, long enough to cross
    /// `MAX_MINIMA_TABLE_BYTES` so that `Minima::new` picks the online strategy.
    fn long_seq(n: usize) -> Vec<u8> {
        let mut state: u64 = 0x2545F491_4F6CDD1D;
        const BASES: [u8; 4] = *b"ACGT";
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                BASES[(state >> 33) as usize & 3]
            })
            .collect()
    }

    /// The precomputed table and the online window must be interchangeable.
    ///
    /// Every other test in the crate uses short sequences and so only exercises
    /// the table; this one runs a sequence past the threshold both ways and
    /// compares the full output, strobe positions included.
    #[test]
    fn minima_strategies_agree() {
        let (k, w_min, w_max) = (21usize, 10usize, 25usize);
        // Just past the point where `Minima::new` stops building a table.
        let seq = long_seq(MAX_MINIMA_TABLE_BYTES / MINIMA_BYTES_PER_KMER + k + 1);

        for n in [2u8, 3] {
            let mut online = MinStrobes::new(&seq, n, k, w_min, w_max).unwrap();
            assert!(
                matches!(online.minima, Minima::Online(_)),
                "expected the online strategy for a {}-base sequence",
                seq.len()
            );

            // Same iterator, but with the minima materialised instead.
            let mut table = MinStrobes::new(&seq, n, k, w_min, w_max).unwrap();
            let (loc, val) = compute_min_hashes(&table.hashes, w_max - w_min + 1);
            table.minima = Minima::Table { loc, val };

            let mut a = Vec::new();
            while let Some(h) = online.next() {
                a.push((online.indexes(), h));
            }
            let mut b = Vec::new();
            while let Some(h) = table.next() {
                b.push((table.indexes(), h));
            }

            assert!(!a.is_empty(), "expected output for order {n}");
            assert_eq!(a, b, "strategies disagree for order {n}");
        }
    }

    /// The strategy is chosen purely by how large the table would be.
    #[test]
    fn strategy_switches_on_table_size() {
        let cutoff = MAX_MINIMA_TABLE_BYTES / MINIMA_BYTES_PER_KMER;
        let small = vec![0u64; cutoff];
        let large = vec![0u64; cutoff + 1];
        assert!(matches!(Minima::new(&small, 4, 0), Minima::Table { .. }));
        assert!(matches!(Minima::new(&large, 4, 0), Minima::Online(_)));
    }

    #[test]
    fn order3_basic() {
        // Basic smoke test: sequence repeated, order=3, k=3, w_min=1, w_max=4
        let seq = "ACGTACGTACGTACGTACGTACGT";
        let ms = MinStrobes::new(seq.as_bytes(), 3, 3, 1, 4).unwrap();
        // Take first 10 strobemers; expect exactly 10 values
        assert_eq!(ms.take(10).count(), 10);
    }
}
