use std::collections::VecDeque;
use nthash_rs::kmer::NtHashBuilder;
use crate::{Result, StrobeError};

/// Generates k-mer hash values from the given sequence `seq`, using exactly one hash per k-mer.
/// 
/// # Parameters
/// - `seq`: byte slice representing the DNA/RNA sequence
/// - `k`: length of each k-mer (must be between 1 and 64, inclusive)
///
/// # Returns
/// - `Ok(Vec<u64>)`: a vector of hash values, one per k-mer, in sequential order
/// - `Err(StrobeError::StrobeLengthTooSmall)`: if `k` is outside the valid range (1..=64)
/// - `Err(StrobeError::SequenceTooShort)`: if `seq.len() < k`
/// - `Err(StrobeError::IncompleteHashValues)`: if the number of hashes computed does not match `seq.len() - k + 1`
///
pub fn compute_hashes(seq: &[u8], k: usize) -> Result<Vec<u64>> {
    // Validate k is within [1, 64]
    if !(1..=64).contains(&k) {
        return Err(StrobeError::StrobeLengthTooSmall);
    }
    // Ensure the sequence length is sufficient for at least one k-mer
    if seq.len() < k {
        return Err(StrobeError::SequenceTooShort);
    }

    // Build the NtHash iterator: `k` specifies the k-mer length,
    // `num_hashes(1)` means we only care about the first hash for each k-mer
    let iter = NtHashBuilder::new(seq)
        .k(k as u16)
        .num_hashes(1)
        .finish()
        .map_err(StrobeError::from)?; // Convert NtHashError into StrobeError

    let capacity = seq.len() - k + 1;
    let mut hashes = Vec::with_capacity(capacity);

    // For each position, `iter` yields (start_index, [hash_values...])
    // We only take the first hash in the returned array.
    for (_, h_array) in iter {
        hashes.push(h_array[0]);
    }

    // Sanity check: we should have exactly one hash per k-mer position
    if hashes.len() != capacity {
        return Err(StrobeError::IncompleteHashValues);
    }
    Ok(hashes)
}

/// For a sliding window of width `w` over the given slice of hash values,
/// computes the index and value of the minimum hash in each window. Runs in O(N) time
/// by maintaining a monotonic (increasing) queue of candidates.
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

    // Trivial case: window size of 1 → each position is its own minimum
    if w == 1 {
        let locs: Vec<usize> = (0..n).collect();
        let mins = hashes.to_vec();
        return (locs, mins);
    }

    // Prepare output vectors: default locs=0, mins=MAX (for positions < w-1)
    let mut locs = vec![0; n];
    let mut mins = vec![u64::MAX; n];
    
    // Monotonic queue: stores (index, value) pairs in increasing order of value
    let mut dq: VecDeque<(usize, u64)> = VecDeque::with_capacity(w);

    for (idx, &h) in hashes.iter().enumerate() {
        // Pop from back while the back's hash is ≥ current hash `h`,
        // ensuring that only strictly smaller values remain earlier in queue.
        while matches!(dq.back(), Some(&(_, back_h)) if back_h >= h) {
            dq.pop_back();
        }
        // Push the current (idx, h) onto the back
        dq.push_back((idx, h));

        // Determine the start of the sliding window that ends at `idx`
        let window_start = idx.saturating_sub(w - 1);

        // Pop from front if that element's index is outside the current window
        while matches!(dq.front(), Some(&(front_idx, _)) if front_idx < window_start) {
            dq.pop_front();
        }

        // Once we've processed at least w elements, record the minimum for this window
        if idx >= w - 1 {
            // The front of the deque is the (index, value) of the minimum
            let (min_pos, min_val) = dq.front().copied().unwrap();
            locs[idx] = min_pos;
            mins[idx] = min_val;
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
