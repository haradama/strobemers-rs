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

        while len > 0 && val_q[(head + len - 1) % w] >= h {
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
