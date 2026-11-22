use crate::{
    Result,
    hashes::{CHUNK_KMERS, HashBuf, KmerHasher, check_hash_count, unsupported_base},
    state::{State, masked_argmin},
    util::validate_params,
};

/// Iterator for generating RandStrobes of order 2 or 3 from a DNA sequence.
///
/// A RandStrobe is a strobemer that selects subsequent k-mers by choosing the
/// position that minimizes `(base_hash + candidate_hash) & prime`. This approach
/// provides a pseudo-random yet deterministic selection of k-mers within sliding windows.
///
/// Built with [`new`](Self::new), the iterator borrows the sequence and hashes
/// it as it goes, so its memory use does not grow with the sequence. Built with
/// [`with_hasher`](Self::with_hasher), every k-mer is hashed up front and the
/// iterator owns the result.
#[derive(Debug, Clone)]
pub struct RandStrobes<'a> {
    state: State<'a>,
}

impl<'a> RandStrobes<'a> {
    /// Constructs a new [`RandStrobes`] iterator using the default hash function (`NtHash64`).
    ///
    /// The sequence is validated here, but hashed lazily while iterating, which
    /// is why the iterator borrows it.
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
    pub fn new(seq: &'a [u8], n: u8, k: usize, w_min: usize, w_max: usize) -> Result<Self> {
        Self::with_chunk(seq, n, k, w_min, w_max, CHUNK_KMERS)
    }

    /// [`new`](Self::new) with an explicit refill size, so that tests can put
    /// buffer boundaries wherever they like.
    fn with_chunk(
        seq: &'a [u8],
        n: u8,
        k: usize,
        w_min: usize,
        w_max: usize,
        chunk: usize,
    ) -> Result<Self> {
        validate_params(seq, n, k, w_min, w_max)?;
        // Hashing is deferred, so anything it would trip over has to be
        // reported now.
        if let Some(e) = unsupported_base(seq) {
            return Err(e);
        }
        let hashes = HashBuf::streaming(seq, k, (n as usize - 1) * w_max, chunk)?;
        Ok(Self {
            state: State::new(seq.len(), n, k, w_min, w_max, hashes),
        })
    }
}

impl RandStrobes<'static> {
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
    /// All k-mers are hashed before this returns, so a hasher that fails or
    /// returns the wrong number of values is reported here, and the iterator
    /// keeps no reference to the sequence or the hasher.
    ///
    /// # Arguments
    ///
    /// * `seq` – Input DNA sequence as ASCII bytes.
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
        H: KmerHasher + ?Sized,
    {
        validate_params(seq, n, k, w_min, w_max)?;

        let hashes = hasher.hash_all(seq, k)?;
        // The iterator indexes the hashes by sequence position, so the hasher
        // must have produced exactly one value per k-mer.
        check_hash_count(&hashes, seq.len(), k)?;

        Ok(Self {
            state: State::new(seq.len(), n, k, w_min, w_max, HashBuf::complete(hashes)),
        })
    }
}

impl RandStrobes<'_> {
    fn advance(&mut self) -> Option<u64> {
        match self.state.n {
            2 => self.next_order2(),
            _ => self.next_order3(),
        }
    }

    /// Computes the next RandStrobe hash value for order 2.
    fn next_order2(&mut self) -> Option<u64> {
        let s = &mut self.state;
        let hi = s.begin()?;
        let (buf, base) = s.hashes.view();
        // `hs[d]` is the hash of the k-mer `d` positions after m1.
        let hs = &buf[s.idx - base..][..=hi];

        let h1 = hs[0];
        // Choose m2 by minimizing `(h1 + hash[m2]) & prime`.
        let m2 = s.w_min + masked_argmin(&hs[s.w_min..], h1, s.prime);
        let hash = (h1 >> 1) + hs[m2] / 3;

        s.idx2 = s.idx + m2;
        s.idx += 1;
        Some(hash)
    }

    /// Computes the next RandStrobe hash value for order 3.
    ///
    /// Both windows are anchored on m1: m2 comes from `w_min..=w_max` after it,
    /// m3 from `w_max + w_min..=2 * w_max`, the latter clipped by the end of
    /// the sequence.
    fn next_order3(&mut self) -> Option<u64> {
        let s = &mut self.state;
        let hi = s.begin()?;
        let (buf, base) = s.hashes.view();
        let hs = &buf[s.idx - base..][..=hi];

        let h1 = hs[0];
        let m2 = s.w_min + masked_argmin(&hs[s.w_min..=s.w_max], h1, s.prime);
        let h2 = h1 / 3 + (hs[m2] >> 2);

        let w2_start = s.w_max + s.w_min;
        let m3 = w2_start + masked_argmin(&hs[w2_start..], h2, s.prime);
        let hash = h2 + hs[m3] / 5;

        s.idx2 = s.idx + m2;
        s.idx3 = s.idx + m3;
        s.idx += 1;
        Some(hash)
    }
}

strobe_iterator!(RandStrobes);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{NtHash64, Strobemer};

    #[test]
    fn order2_basic() {
        let mut rs = RandStrobes::new("ACGTACGTACGT".as_bytes(), 2, 3, 1, 4).unwrap();
        assert!(rs.next().is_some());
    }

    #[test]
    fn order3_basic() {
        let seq = "ACGTACGTACGTACGTACGTACGT";
        let rs = RandStrobes::new(seq.as_bytes(), 3, 3, 1, 4).unwrap();
        // The sequence yields more strobemers than are taken here.
        assert_eq!(rs.take(10).count(), 10);
    }

    /// Deterministic pseudo-random DNA.
    fn random_seq(n: usize) -> Vec<u8> {
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

    /// Refilling the hash buffer must be invisible, wherever the boundaries
    /// fall relative to the search windows and whatever the refill size.
    /// Hashing everything up front is the reference.
    #[test]
    fn chunk_boundaries_are_invisible() {
        let seq = random_seq(700);
        for n in [2u8, 3] {
            for (k, w_min, w_max) in [(3usize, 1usize, 1usize), (5, 2, 9), (7, 4, 4), (4, 1, 30)] {
                for shrink in [true, false] {
                    let mut want =
                        RandStrobes::with_hasher(&seq, n, k, w_min, w_max, &NtHash64).unwrap();
                    assert!(want.state.hashes.is_complete());
                    want.set_window_shrink(shrink);
                    let want: Vec<Strobemer> = want.with_positions().collect();
                    assert!(!want.is_empty());

                    for chunk in [1usize, 2, 3, 17, 64, 699, 5000] {
                        let mut got =
                            RandStrobes::with_chunk(&seq, n, k, w_min, w_max, chunk).unwrap();
                        got.set_window_shrink(shrink);
                        let got: Vec<Strobemer> = got.with_positions().collect();
                        assert_eq!(
                            got, want,
                            "n={n} k={k} w=({w_min},{w_max}) shrink={shrink} chunk={chunk}"
                        );
                    }
                }
            }
        }
    }

    /// The default constructor streams once the sequence outgrows one refill,
    /// and still agrees with hashing everything up front.
    #[test]
    fn long_sequences_are_streamed() {
        let (k, w_min, w_max) = (21usize, 10usize, 25usize);
        let seq = random_seq(4 * CHUNK_KMERS + 123);
        for n in [2u8, 3] {
            let streamed = RandStrobes::new(&seq, n, k, w_min, w_max).unwrap();
            assert!(!streamed.state.hashes.is_complete());
            let a: Vec<Strobemer> = streamed.with_positions().collect();
            let b: Vec<Strobemer> = RandStrobes::with_hasher(&seq, n, k, w_min, w_max, &NtHash64)
                .unwrap()
                .with_positions()
                .collect();
            assert!(!a.is_empty());
            assert!(a == b, "streaming disagrees for order {n}");
        }
    }
}
