use crate::{Result, StrobeError};
use nthash_rs::kmer::NtHashBuilder;

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

        let it = NtHashBuilder::new(seq)
            .k(k as u16)
            .num_hashes(1)
            .finish()
            .map_err(StrobeError::from)?;

        let mut out = Vec::with_capacity(seq.len() - k + 1);
        for (_, h) in it {
            out.push(h[0]);
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

/// For a sliding window of width `w` over the given slice of hash values,
/// computes the index and value of the minimum hash in each window.
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

    for (i, &h) in hashes.iter().enumerate() {
        let window_start = i.saturating_sub(w - 1);
        while len > 0 && idx_q[head] < window_start {
            head = (head + 1) % w;
            len -= 1;
        }

        // Strictly-greater comparison: equal values are kept, so the front of the
        // deque holds the *earliest* position among ties. The Go reference keeps
        // a sorted buffer that inserts after equal elements, giving the same
        // choice; using `>=` here would report the latest position instead.
        while len > 0 && val_q[(head + len - 1) % w] > h {
            len -= 1;
        }

        let tail = (head + len) % w;
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
