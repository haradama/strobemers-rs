use crate::{
    constants::DEFAULT_PRIME_NUMBER,
    hashes::{compute_hashes, compute_min_hashes},
    util::roundup64,
    Result, StrobeError,
};

/// Iterator for generating MinStrobes of order 2 or 3 from a DNA/RNA sequence.
///
/// A MinStrobe is a concatenation of k-mers selected based on minimum hash
/// values within sliding windows. This struct precomputes k-mer hashes and
/// window minima to efficiently produce strobemer hash values.
///
#[derive(Debug, Clone)]
pub struct MinStrobes {
    // Parameters controlling strobemer generation
    n:      u8,       // Order of strobemer: 2 or 3
    w_min:  usize,    // Minimum window offset
    w_max:  usize,    // Maximum window offset

    // Precomputed data
    hashes: Vec<u64>,   // Hash values for each k-mer in the sequence
    minloc: Vec<usize>, // Location of the minimum hash within each sliding window
    minval: Vec<u64>,   // Minimum hash value within each sliding window

    // Iteration state
    idx:      usize, // Current index of the first k-mer (m1)
    end_idx:  usize, // Last index at which a complete strobemer can start
    end_hash: usize, // Last index in `hashes` (i.e., sequence length minus k)

    // Strobe indices for current item
    idx2: usize, // Index of second k-mer (m2)
    idx3: usize, // Index of third k-mer (m3) if order = 3

    // Prime number and shrink-window flag
    prime: u64,  // Used for combining hash values in order 3
    shrink: bool, // Whether to shrink windows near sequence end

    // Working registers for hash values
    h1: u64, // Hash of first k-mer (m1)
    h2: u64, // Combined hash after selecting m2
    h3: u64, // Combined hash after selecting m3 (order 3 only)
}

impl MinStrobes {
    /// Constructs a new `MinStrobes` iterator.
    ///
    /// Precomputes k-mer hash values and per-window minima for efficient strobemer generation.
    ///
    /// # Arguments
    ///
    /// * `seq` – Byte slice of the DNA/RNA sequence.
    /// * `n` – Order of strobemer (must be 2 or 3).
    /// * `l` – k-mer length for each strobe (must be between 1 and 64).
    /// * `w_min` – Minimum window offset (must be ≥ 1).
    /// * `w_max` – Maximum window offset (must be ≥ `w_min`).
    ///
    /// # Returns
    ///
    /// * `Ok(Self)` – Successfully initialized `MinStrobes`.
    /// * `Err(StrobeError)` – Parameter validation or hashing failure.
    ///
    pub fn new(seq: &[u8], n: u8, l: usize, w_min: usize, w_max: usize) -> Result<Self> {
        // Validate input parameters (sequence, order, window sizes, k-mer length)
        validate_params!(seq, n, l, w_min, w_max);

        // --- Precompute k-mer hashes and window minima ---
        let hashes = compute_hashes(seq, l)?;
        let (minloc, minval) = compute_min_hashes(&hashes, w_max - w_min + 1);

        // Compute the last valid index for selecting m1 such that all strobes fit
        let seq_len = seq.len();
        let end_hash = seq_len - l; // last index where a k-mer hash exists
        let end_idx = seq_len - l - (n as usize - 1) * l;

        Ok(Self {
            n,
            w_min,
            w_max,
            hashes,
            minloc,
            minval,
            idx: 0,
            end_hash,
            end_idx,
            idx2: 0,
            idx3: 0,
            prime: DEFAULT_PRIME_NUMBER, // default Mersenne prime (2^20 - 1)
            shrink: true,                // allow shrinking windows near the end
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
        // Round up to next power of two, subtract one → Mersenne prime form
        self.prime = roundup64(q) - 1;
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

    /// Computes the next hash value for an order-2 MinStrobe.
    fn next_order2(&mut self) -> Option<u64> {
        // Stop if no more valid starting positions for m1
        if self.idx > self.end_idx {
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

        // If full window fits, use precomputed minimum
        if w_end == self.idx + self.w_max {
            self.idx2 = self.minloc[w_end];
            // Combine h1 and precomputed minimum hash
            self.h2 = self.h1 / 2 + self.minval[w_end] / 3;
        } else {
            // Partial window: manually scan to find minimum
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
        if self.idx > self.end_idx {
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
        // Select m2 using precomputed minima at window end
        self.idx2 = self.minloc[w_end];
        self.h2 = self.h1 / 3 + self.minval[w_end] / 4;

        // Select m3
        if w2_end == self.idx + (self.w_max << 1) {
            // Full second window fits: use precomputed minima
            self.idx3 = self.minloc[w2_end];
            self.h3 = self.h2 + self.minval[w2_end] / 5;
        } else {
            // Partial second window near the end: manual scan
            let (mut best_hash, mut best_pos) = (u64::MAX, w2_start);
            for pos in w2_start..=w2_end {
                // Combine current h2 with candidate hash, then mask with prime
                let cand = (self.h2 + self.hashes[pos]) & self.prime;
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

    #[test]
    fn order3_basic() {
        // Basic smoke test: sequence repeated, order=3, k=3, w_min=1, w_max=4
        let seq = "ACGTACGTACGTACGTACGTACGT";
        let ms = MinStrobes::new(seq.as_bytes(), 3, 3, 1, 4).unwrap();
        // Take first 10 strobemers; expect exactly 10 values
        assert_eq!(ms.take(10).count(), 10);
    }
}
