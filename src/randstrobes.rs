use crate::{
    constants::DEFAULT_PRIME_NUMBER,
    hashes::compute_hashes,
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
    _l:     usize,   // k-mer length (only needed during construction)
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
    /// Constructs a new `RandStrobes` iterator.
    ///
    /// Precomputes k-mer hashes and initializes state for pseudo-random strobemer selection.
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
    /// * `Ok(Self)` – Successfully initialized `RandStrobes`.
    /// * `Err(StrobeError)` – Parameter validation or hashing failure.
    ///
    pub fn new(seq: &[u8], n: u8, l: usize, w_min: usize, w_max: usize) -> Result<Self> {
        // Validate input parameters (sequence, order, window sizes, k-mer length)
        validate_params!(seq, n, l, w_min, w_max);

        // Precompute k-mer hash values
        let hashes = compute_hashes(seq, l)?;

        // Compute last valid hash index (sequence length minus k-mer length)
        let end_hash = seq.len().saturating_sub(l);
        // Compute last valid starting index for a full strobemer (n*k-mer segments)
        let end_idx = seq.len().saturating_sub(l + (n as usize - 1) * l);

        Ok(Self {
            n,
            _l: l,
            w_min,
            w_max,
            hashes,
            idx: 0,
            end_idx,
            end_hash,
            idx2: 0,
            idx3: 0,
            prime: DEFAULT_PRIME_NUMBER, // default Mersenne prime (2^20 − 1)
            shrink: true,                // allow shrinking windows near the end
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
