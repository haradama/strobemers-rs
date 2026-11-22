use crate::{
    Result,
    hashes::{
        CHUNK_KMERS, HashBuf, KmerHasher, MinWindow, check_hash_count, compute_min_hashes,
        unsupported_base,
    },
    state::{State, masked_argmin},
    util::validate_params,
};

/// Largest sliding-window-minima table this module will precompute, in bytes.
///
/// The table holds a `usize` and a `u64` per k-mer, on top of the hash array
/// itself. Below this the bulk loop that fills it is the faster way to get the
/// minima; above it the extra memory traffic costs more than recomputing them
/// as the iterator advances.
///
/// The threshold was chosen where the two strategies were measured to cross;
/// the `scaling` benchmark group shows the crossover, so re-run it before
/// changing this value.
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
    /// one pair per k-mer. Faster once the table would no longer fit in cache,
    /// and its memory use does not grow with the sequence.
    Online(MinWindow),
}

impl Minima {
    /// Whether a sequence with `kmers` k-mers is short enough for the table.
    fn fits_table(kmers: usize) -> bool {
        kmers.saturating_mul(MINIMA_BYTES_PER_KMER) <= MAX_MINIMA_TABLE_BYTES
    }

    /// Chooses a strategy for the given hashes.
    ///
    /// The table needs every hash up front, so it is only an option when they
    /// are all buffered. `lag` is how far behind the most recent position the
    /// caller may still read; order 3 needs `w_max`, order 2 needs none.
    fn new(hashes: &HashBuf<'_>, window: usize, lag: usize) -> Self {
        let (buf, _) = hashes.view();
        if hashes.is_complete() && Self::fits_table(buf.len()) {
            let (loc, val) = compute_min_hashes(buf, window);
            Self::Table { loc, val }
        } else {
            Self::Online(MinWindow::new(window, lag))
        }
    }

    /// `(position, value)` of the minimum over the window ending at `p`,
    /// where `buf[i]` is the hash at position `base + i`.
    ///
    /// For the online variant `p` must not go backwards by more than `lag`
    /// between calls, which both call sites respect.
    #[inline]
    fn at(&mut self, buf: &[u64], base: usize, p: usize) -> (usize, u64) {
        match self {
            Self::Table { loc, val } => (loc[p], val[p]),
            Self::Online(win) => win.advance_to(buf, base, p),
        }
    }
}

/// Iterator for generating MinStrobes of order 2 or 3 from a DNA sequence.
///
/// A MinStrobe is a concatenation of k-mers selected based on minimum hash
/// values within sliding windows.
///
/// Built with [`new`](Self::new), the iterator borrows the sequence; for long
/// sequences it hashes as it goes, so its memory use does not grow with the
/// sequence. Built with [`with_hasher`](Self::with_hasher), every k-mer is
/// hashed up front and the iterator owns the result.
#[derive(Debug, Clone)]
pub struct MinStrobes<'a> {
    state: State<'a>,
    /// Sliding-window minima, precomputed or computed online.
    minima: Minima,
}

/// Window width shared by both strobe selections, and how far back the minima
/// are read: order 3 reads them at two positions `w_max` apart, order 2 only
/// ever at the position it just reached.
fn minima_for(hashes: &HashBuf<'_>, n: u8, w_min: usize, w_max: usize) -> Minima {
    Minima::new(hashes, w_max - w_min + 1, if n == 3 { w_max } else { 0 })
}

impl<'a> MinStrobes<'a> {
    /// Constructs a new [`MinStrobes`] iterator using the default hash function (`NtHash64`).
    ///
    /// The sequence is validated here. Short sequences are hashed here too;
    /// long ones are hashed lazily while iterating, which is why the iterator
    /// borrows the sequence.
    ///
    /// # Arguments
    ///
    /// * `seq` – Input nucleotide sequence as a byte slice (DNA, ASCII only).
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
    pub fn new(seq: &'a [u8], n: u8, k: usize, w_min: usize, w_max: usize) -> Result<Self> {
        validate_params(seq, n, k, w_min, w_max)?;
        // While the minima table is worth building, all hashes are needed up
        // front anyway; beyond that, hash in chunks.
        let kmers = seq.len() - k + 1;
        let chunk = if Minima::fits_table(kmers) {
            kmers
        } else {
            CHUNK_KMERS
        };
        Self::with_chunk(seq, n, k, w_min, w_max, chunk)
    }

    /// Builds the iterator with an explicit refill size, so that tests can put
    /// buffer boundaries wherever they like. The parameters must already be
    /// validated.
    fn with_chunk(
        seq: &'a [u8],
        n: u8,
        k: usize,
        w_min: usize,
        w_max: usize,
        chunk: usize,
    ) -> Result<Self> {
        // Hashing may be deferred, so anything it would trip over has to be
        // reported now.
        if let Some(e) = unsupported_base(seq) {
            return Err(e);
        }
        let hashes = HashBuf::streaming(seq, k, (n as usize - 1) * w_max, chunk)?;
        let minima = minima_for(&hashes, n, w_min, w_max);
        Ok(Self {
            state: State::new(seq.len(), n, k, w_min, w_max, hashes),
            minima,
        })
    }
}

impl MinStrobes<'static> {
    /// Constructs a new [`MinStrobes`] iterator with a user-defined hash function.
    ///
    /// This method accepts any implementation of the [`KmerHasher`] trait,
    /// allowing for custom k-mer hashing strategies, such as XOR-based, cryptographic,
    /// or rolling hashes optimized for performance or reproducibility.
    ///
    /// All k-mers are hashed before this returns, so a hasher that fails or
    /// returns the wrong number of values is reported here, and the iterator
    /// keeps no reference to the sequence or the hasher.
    ///
    /// # Arguments
    ///
    /// * `seq` – Input DNA sequence as bytes (e.g., `b"ACGT..."`).
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
        H: KmerHasher + ?Sized,
    {
        validate_params(seq, n, k, w_min, w_max)?;

        let hashes = hasher.hash_all(seq, k)?;
        // The iterator indexes the hashes by sequence position, so the hasher
        // must have produced exactly one value per k-mer.
        check_hash_count(&hashes, seq.len(), k)?;

        let hashes = HashBuf::complete(hashes);
        let minima = minima_for(&hashes, n, w_min, w_max);
        Ok(Self {
            state: State::new(seq.len(), n, k, w_min, w_max, hashes),
            minima,
        })
    }
}

impl MinStrobes<'_> {
    fn advance(&mut self) -> Option<u64> {
        match self.state.n {
            2 => self.next_order2(),
            _ => self.next_order3(),
        }
    }

    /// Computes the next hash value for an order-2 MinStrobe.
    fn next_order2(&mut self) -> Option<u64> {
        let s = &mut self.state;
        let hi = s.begin()?;
        let (buf, base) = s.hashes.view();
        let h1 = buf[s.idx - base];

        let (m2, val) = if hi == s.reach {
            // The full window fits: take the running minimum. Once the window
            // starts being clipped it stays clipped, so the minima are never
            // needed again and leaving them behind is safe.
            self.minima.at(buf, base, s.idx + hi)
        } else {
            // Clipped window: scan what is left of it.
            let start = s.idx + s.w_min;
            let (off, val) = earliest_min(&buf[start - base..=s.idx + hi - base]);
            (start + off, val)
        };

        s.idx2 = m2;
        s.idx += 1;
        Some((h1 >> 1) + val / 3)
    }

    /// Computes the next hash value for an order-3 MinStrobe.
    ///
    /// Both windows are anchored on m1: m2 comes from `w_min..=w_max` after it,
    /// m3 from `w_max + w_min..=2 * w_max`, the latter clipped by the end of
    /// the sequence.
    fn next_order3(&mut self) -> Option<u64> {
        let s = &mut self.state;
        let hi = s.begin()?;
        let (buf, base) = s.hashes.view();
        let h1 = buf[s.idx - base];

        // The first window always fits when the second one is not empty.
        let (m2, val2) = self.minima.at(buf, base, s.idx + s.w_max);
        let h2 = h1 / 3 + (val2 >> 2);

        let (m3, val3) = if hi == s.reach {
            self.minima.at(buf, base, s.idx + hi)
        } else {
            // Clipped second window. As in the reference, this case does not
            // take the plain minimum but the RandStrobes-style masked one.
            let start = s.idx + s.w_max + s.w_min;
            let window = &buf[start - base..=s.idx + hi - base];
            let off = masked_argmin(window, h2, s.prime);
            (start + off, window[off])
        };

        s.idx2 = m2;
        s.idx3 = m3;
        s.idx += 1;
        Some(h2 + val3 / 5)
    }
}

/// Offset and value of the earliest minimum of `window`, which must not be empty.
#[inline]
fn earliest_min(window: &[u64]) -> (usize, u64) {
    let (mut best_off, mut best_val) = (0, window[0]);
    for (off, &h) in window.iter().enumerate() {
        if h < best_val {
            best_val = h;
            best_off = off;
        }
    }
    (best_off, best_val)
}

strobe_iterator!(MinStrobes);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{NtHash64, Strobemer};

    #[test]
    fn order2_basic() {
        let mut ms = MinStrobes::new("ACGTACGTACGT".as_bytes(), 2, 3, 1, 4).unwrap();
        assert!(ms.next().is_some());
    }

    #[test]
    fn order3_basic() {
        let seq = "ACGTACGTACGTACGTACGTACGT";
        let ms = MinStrobes::new(seq.as_bytes(), 3, 3, 1, 4).unwrap();
        // The sequence yields more strobemers than are taken here.
        assert_eq!(ms.take(10).count(), 10);
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

    /// Every k-mer hashed up front, with the minima materialised as a table:
    /// the plainest configuration, used as the reference below.
    fn eager_table(seq: &[u8], n: u8, k: usize, w_min: usize, w_max: usize) -> MinStrobes<'static> {
        let mut ms = MinStrobes::with_hasher(seq, n, k, w_min, w_max, &NtHash64).unwrap();
        let (buf, base) = ms.state.hashes.view();
        assert_eq!(base, 0);
        let (loc, val) = compute_min_hashes(buf, w_max - w_min + 1);
        ms.minima = Minima::Table { loc, val };
        ms
    }

    /// The precomputed table and the online window must be interchangeable,
    /// and so must hashing up front and hashing in chunks.
    ///
    /// Every other test in the crate uses short sequences and so only exercises
    /// the table; this one runs a sequence past the threshold, where
    /// [`MinStrobes::new`] streams the hashes and computes the minima online,
    /// and compares the full output, strobe positions included.
    #[test]
    fn long_sequences_match_the_table() {
        let (k, w_min, w_max) = (21usize, 10usize, 25usize);
        // Just past the point where the table stops being built.
        let seq = random_seq(MAX_MINIMA_TABLE_BYTES / MINIMA_BYTES_PER_KMER + k + 1);

        for n in [2u8, 3] {
            let online = MinStrobes::new(&seq, n, k, w_min, w_max).unwrap();
            assert!(
                matches!(online.minima, Minima::Online(_)),
                "expected the online strategy for a {}-base sequence",
                seq.len()
            );
            assert!(
                !online.state.hashes.is_complete(),
                "expected the hashes to be streamed"
            );

            let a: Vec<Strobemer> = online.with_positions().collect();
            let b: Vec<Strobemer> = eager_table(&seq, n, k, w_min, w_max)
                .with_positions()
                .collect();

            assert!(!a.is_empty(), "expected output for order {n}");
            assert!(a == b, "strategies disagree for order {n}");
        }
    }

    /// Refilling the hash buffer must be invisible, wherever the boundaries
    /// fall relative to the search windows and whatever the refill size.
    #[test]
    fn chunk_boundaries_are_invisible() {
        let seq = random_seq(700);
        for n in [2u8, 3] {
            for (k, w_min, w_max) in [(3usize, 1usize, 1usize), (5, 2, 9), (7, 4, 4), (4, 1, 30)] {
                for shrink in [true, false] {
                    let mut want = eager_table(&seq, n, k, w_min, w_max);
                    want.set_window_shrink(shrink);
                    let want: Vec<Strobemer> = want.with_positions().collect();
                    assert!(!want.is_empty());

                    for chunk in [1usize, 2, 3, 17, 64, 699, 5000] {
                        let mut got =
                            MinStrobes::with_chunk(&seq, n, k, w_min, w_max, chunk).unwrap();
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

    /// The strategy is chosen purely by how large the table would be.
    #[test]
    fn strategy_switches_on_table_size() {
        let cutoff = MAX_MINIMA_TABLE_BYTES / MINIMA_BYTES_PER_KMER;
        let small = HashBuf::complete(vec![0u64; cutoff]);
        let large = HashBuf::complete(vec![0u64; cutoff + 1]);
        assert!(matches!(Minima::new(&small, 4, 0), Minima::Table { .. }));
        assert!(matches!(Minima::new(&large, 4, 0), Minima::Online(_)));
    }
}
