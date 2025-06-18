use crate::{
    constants::DEFAULT_PRIME_NUMBER,
    hashes::{KmerHasher, NtHash64},
    util::roundup64,
    Result, StrobeError,
};

/// Iterator for generating RandStrobes of order 2 or 3 from a DNA/RNA sequence.
///
/// A RandStrobe is a strobemer that selects subsequent k-mers by choosing the
/// position that minimizes `(base_hash + candidate_hash) & prime`. This approach
/// provides a pseudo-random yet deterministic selection of k-mers within sliding windows.
///
#[derive(Debug, Clone)]
pub struct RandStrobes {
    // Parameters controlling strobemer generation
    n:      u8,      // Order of strobemer: 2 or 3
    _k:     usize,   // k-mer length (only needed during construction)
    w_min:  usize,   // Minimum window offset
    w_max:  usize,   // Maximum window offset

    // Precomputed data
    hashes: Vec<u64>, // Hash values for each k-mer in the sequence

    // Iteration state
    idx:      usize, // Current index of the first k-mer (m1)
    end_idx:  usize, // Last index at which a complete strobemer can start
    end_hash: usize, // Last index in `hashes` (i.e., sequence length minus k)

    // Strobe indices for current item
    idx2: usize, // Index of second k-mer (m2)
    idx3: usize, // Index of third k-mer (m3) if order = 3

    // Prime number and shrink-window flag
    prime: u64,  // Used for mask-based combination: `(base_hash + candidate_hash) & prime`
    shrink: bool, // Whether to shrink windows near the end if the full window does not fit

    // Working registers for hash values
    h1: u64, // Hash of first k-mer (m1)
    h2: u64, // Combined hash after selecting m2
    h3: u64, // Combined hash after selecting m3 (order 3 only)
}

impl RandStrobes {
    /// Constructs a new [`RandStrobes`] iterator using the default hash function (`NtHash64`).
    ///
    /// This method serves as a convenience wrapper for [`RandStrobes::with_hasher`],
    /// providing a standard ntHash-based setup for k-mer hashing.
    ///
    /// The generated iterator will produce strobemers using the **RandStrobe protocol**,
    /// where the second (and optionally third) k-mer is selected based on a minimum
    /// of a randomized hash function over a windowed region.
    ///
    /// # Arguments
    ///
    /// * `seq` – Nucleotide sequence as a byte slice (e.g., `b"ACGT..."`). Must be ASCII.
    /// * `n` – Strobemer order (2 or 3 only).
    /// * `k` – k-mer length for each strobe. Must be between 1 and 64 (inclusive).
    /// * `w_min` – Minimum window offset for selecting the next strobe.
    /// * `w_max` – Maximum window offset (inclusive); must satisfy `w_min ≤ w_max`.
    ///
    /// # Returns
    ///
    /// * `Ok(RandStrobes)` – Ready-to-use iterator for random strobemers.
    /// * `Err(StrobeError)` – Returned if parameters are invalid or the sequence is too short.
    ///
    /// # Example
    /// ```
    /// use strobemers_rs::RandStrobes;
    /// let rs = RandStrobes::new(b"ACGTACGTACGT", 2, 3, 1, 4).unwrap();
    /// for h in rs.take(5) {
    ///     println!("{}", h);
    /// }
    /// ```
    pub fn new(seq: &[u8], n: u8, k: usize, w_min: usize, w_max: usize) -> Result<Self> {
        Self::with_hasher(seq, n, k, w_min, w_max, &NtHash64)
    }

    /// Constructs a new [`RandStrobes`] iterator using a user-defined k-mer hash function.
    ///
    /// This method enables **dependency injection** of the hashing algorithm via the [`KmerHasher`] trait.
    /// It allows experimentation with custom hashers (e.g. fast XOR, cryptographic, locality-aware) for advanced use cases.
    ///
    /// The resulting iterator emits strobemer hashes using the **RandStrobe method**:
    /// - The first k-mer is fixed at position `i`
    /// - The next k-mer is chosen within a window `[i + w_min ..= i + w_max]`
    ///   to **minimize a masked combination** `(h₁ + h₂) & prime`
    /// - If `n = 3`, the third k-mer is chosen similarly after `w_max + w_min`
    ///
    /// # Arguments
    ///
    /// * `seq` – Input DNA/RNA sequence as ASCII bytes.
    /// * `n` – Order of the strobemer (must be 2 or 3).
    /// * `k` – Length of each strobe (k-mer), within the inclusive range [1, 64].
    /// * `w_min` – Minimum offset for the search window (must be ≥ 1).
    /// * `w_max` – Maximum offset (inclusive); must satisfy `w_min ≤ w_max`.
    /// * `hasher` – Reference to a [`KmerHasher`] implementation for computing all k-mer hashes.
    ///
    /// # Returns
    ///
    /// * `Ok(RandStrobes)` – If input and hashes are valid.
    /// * `Err(StrobeError)` – On invalid input, hashing errors, or insufficient sequence length.
    ///
    /// # Example
    /// ```
    /// use strobemers_rs::{RandStrobes, KmerHasher};
    ///
    /// struct DummyHasher;
    /// impl KmerHasher for DummyHasher {
    ///     fn hash_all(&self, seq: &[u8], k: usize) -> strobemers_rs::Result<Vec<u64>> {
    ///         Ok(seq.windows(k).map(|w| w.iter().map(|b| *b as u64).sum()).collect())
    ///     }
    /// }
    ///
    /// let rs = RandStrobes::with_hasher(b"ACGTACGT", 2, 3, 1, 4, &DummyHasher).unwrap();
    /// for h in rs.take(3) {
    ///     println!("strobemer hash: {}", h);
    /// }
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
        // Ensure all parameters are valid before proceeding
        validate_params!(seq, n, k, w_min, w_max);

        // Precompute hash values for all valid k-mers
        let hashes = hasher.hash_all(seq, k)?;

        // Calculate the valid iteration bounds
        let end_hash = seq.len().saturating_sub(k); // maximum hash index
        let end_idx = seq.len().saturating_sub(k + (n as usize - 1) * k); // max starting index for m₁

        Ok(Self {
            n,
            _k: k,
            w_min,
            w_max,
            hashes,
            idx: 0,
            end_idx,
            end_hash,
            idx2: 0,
            idx3: 0,
            prime: DEFAULT_PRIME_NUMBER,
            shrink: true,
            h1: 0,
            h2: 0,
            h3: 0,
        })
    }

    /// Sets a new prime number for combining hash values.
    ///
    /// The formula used is `(base_hash + candidate_hash) & prime`. The provided `q` must
    /// be at least 256. Internally, `q` is rounded up to the next power of two and then
    /// decremented by one to form a Mersenne prime.
    ///
    /// # Arguments
    ///
    /// * `q` – Candidate prime (will be rounded to Mersenne form: `2^k − 1`).
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

    /// Chooses the position within `range` that minimizes `(base_hash + hashes[pos]) & prime`.
    ///
    /// # Arguments
    ///
    /// * `base_hash` – The hash value of the previous strobe (m1 or m2).
    /// * `range` – Inclusive range of indices to consider for the next strobe.
    ///
    /// # Returns
    ///
    /// *(best_pos, best_val)* – Index of the chosen k-mer and the resulting combined hash value.
    ///
    fn choose_min(&self, base_hash: u64, range: std::ops::RangeInclusive<usize>) -> (usize, u64) {
        let mut best_pos = *range.start();
        let mut best_val = u64::MAX;

        for pos in range {
            // Wrap-around addition, then bitwise AND with prime (Mersenne prime mask)
            let cand = base_hash
                .wrapping_add(self.hashes[pos])
                & self.prime;
            if cand < best_val {
                best_val = cand;
                best_pos = pos;
            }
        }
        (best_pos, best_val)
    }

    // -------------------- order-specific next ---------------------------- //

    /// Computes the next RandStrobe hash value for order 2.
    ///
    /// # Returns
    /// - `Some(u64)` – Combined hash of m1 and m2, if available.
    /// - `None` – When `idx > end_idx` (no more valid strobes).
    ///
    fn next_order2(&mut self) -> Option<u64> {
        if self.idx > self.end_idx {
            return None;
        }

        // Define the search window for m2
        let w_start = self.idx + self.w_min;
        let mut w_end = self.idx + self.w_max;
        if w_end > self.end_hash {
            if !self.shrink {
                return None;
            }
            w_end = self.end_hash;
        }

        // Hash of the first k-mer (m1)
        self.h1 = self.hashes[self.idx];
        // Choose m2 by minimizing `(h1 + hash[m2]) & prime`
        let (pos2, _) = self.choose_min(self.h1, w_start..=w_end);
        self.idx2 = pos2;
        // Combine h1 and second k-mer’s hash
        self.h2 = self.h1 / 2 + self.hashes[pos2] / 3;

        // Advance to next starting index for m1
        self.idx += 1;
        Some(self.h2)
    }

    /// Computes the next RandStrobe hash value for order 3.
    ///
    /// # Returns
    /// - `Some(u64)` – Combined hash of m1, m2, and m3, if available.
    /// - `None` – When no further strobes can be formed.
    ///
    fn next_order3(&mut self) -> Option<u64> {
        if self.idx > self.end_idx {
            return None;
        }

        // First window range for selecting m2
        let w1_start = self.idx + self.w_min;
        let w1_end = self.idx + self.w_max;

        // Second window range for selecting m3
        let w2_start = self.idx + self.w_max + self.w_min;
        let mut w2_end = self.idx + (self.w_max << 1);
        if w2_start > self.end_hash {
            return None;
        }
        if w2_end > self.end_hash {
            if !self.shrink {
                return None;
            }
            w2_end = self.end_hash;
        }

        // Compute m1 (first k-mer)
        self.h1 = self.hashes[self.idx];
        // Select m2
        let (pos2, _) = self.choose_min(self.h1, w1_start..=w1_end);
        self.idx2 = pos2;
        self.h2 = self.h1 / 3 + self.hashes[pos2] / 4;

        // Select m3
        let (pos3, _) = self.choose_min(self.h2, w2_start..=w2_end);
        self.idx3 = pos3;
        self.h3 = self.h2 + self.hashes[pos3] / 5;

        // Advance to next starting index for m1
        self.idx += 1;
        Some(self.h3)
    }
}

impl Iterator for RandStrobes {
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
        let mut rs = RandStrobes::new("ACGTACGTACGT".as_bytes(), 2, 3, 1, 4).unwrap();
        // Expect at least one strobemer
        assert!(rs.next().is_some());
    }

    #[test]
    fn order3_basic() {
        // Basic smoke test: sequence repeated, order=3, k=3, w_min=1, w_max=4
        let seq = "ACGTACGTACGTACGTACGTACGT";
        let rs = RandStrobes::new(seq.as_bytes(), 3, 3, 1, 4).unwrap();
        // Take first 10 strobemers; expect exactly 10 values
        assert_eq!(rs.take(10).count(), 10);
    }
}
