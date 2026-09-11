//! Regression tests for defects found while auditing this crate against the Go
//! reference implementation (<https://github.com/shenwei356/strobemers>).
//!
//! Each test pins one behaviour that previously panicked, produced an
//! out-of-range strobe position, or diverged from the reference.

use strobemers_rs::{KmerHasher, MinStrobes, RandStrobes, Result, StrobeError};

/// Hasher that returns values near `u64::MAX`, so that the order-3 combination
/// `h2 + hashes[pos]` overflows a `u64`.
struct HugeHasher;

impl KmerHasher for HugeHasher {
    fn hash_all(&self, seq: &[u8], k: usize) -> Result<Vec<u64>> {
        Ok((0..seq.len() - k + 1)
            .map(|i| u64::MAX - (i as u64 % 7))
            .collect())
    }
}

/// Hasher that returns one hash too few, mimicking a hasher that skips k-mers.
struct ShortHasher;

impl KmerHasher for ShortHasher {
    fn hash_all(&self, seq: &[u8], k: usize) -> Result<Vec<u64>> {
        let mut v: Vec<u64> = (0..seq.len() - k + 1).map(|i| i as u64).collect();
        v.pop();
        Ok(v)
    }
}

/// `nthash-rs` skips every k-mer window containing a non-ACGT base, so the hash
/// vector becomes shorter than the sequence and its indices no longer line up
/// with sequence positions. That used to yield strobes at shifted positions and
/// then panic with an out-of-bounds index; it must be reported as an error, as
/// the Go reference does with `ErrIncompleteHashValues`.
#[test]
fn ambiguous_base_is_rejected() {
    let with_n: &[u8] = b"ACGATCTGGTANCTAGGCTTACGATCAGT";
    let clean: &[u8] = b"ACGATCTGGTACCTAGGCTTACGATCAGT";

    for n in [2u8, 3] {
        assert_eq!(
            RandStrobes::new(with_n, n, 3, 3, 5).err(),
            Some(StrobeError::IncompleteHashValues),
            "RandStrobes must reject a sequence containing N (order {n})"
        );
        assert_eq!(
            MinStrobes::new(with_n, n, 3, 3, 5).err(),
            Some(StrobeError::IncompleteHashValues),
            "MinStrobes must reject a sequence containing N (order {n})"
        );
        // The same sequence without the ambiguous base still works.
        assert!(RandStrobes::new(clean, n, 3, 3, 5).is_ok());
        assert!(MinStrobes::new(clean, n, 3, 3, 5).is_ok());
    }
}

/// A custom `KmerHasher` returning the wrong number of hashes must be rejected
/// rather than allowed to drive the iterator out of bounds.
#[test]
fn hasher_returning_too_few_hashes_is_rejected() {
    let seq: &[u8] = b"ACGATCTGGTACCTAGGCTTACGATCAGT";
    assert_eq!(
        RandStrobes::with_hasher(seq, 2, 3, 1, 4, &ShortHasher).err(),
        Some(StrobeError::IncompleteHashValues)
    );
    assert_eq!(
        MinStrobes::with_hasher(seq, 2, 3, 1, 4, &ShortHasher).err(),
        Some(StrobeError::IncompleteHashValues)
    );
}

/// Once `w_min > (n - 1) * l`, the last few m1 positions have a search window
/// that starts past the final k-mer. Iteration must stop there. Previously
/// `RandStrobes` panicked and `MinStrobes` emitted a strobe whose `m2` index
/// pointed past the end of the sequence.
#[test]
fn empty_second_window_ends_iteration() -> Result<()> {
    let seq: &[u8] = b"ACGATCTGGTACCTAG"; // len 16
    let l = 2;
    let end_hash = seq.len() - l; // index of the last l-mer == 14

    for (w_min, w_max, expected) in [(3usize, 5usize, 12usize), (4, 4, 11), (5, 9, 10)] {
        let mut rand_indexes = Vec::new();
        let mut rs = RandStrobes::new(seq, 2, l, w_min, w_max)?;
        while rs.next().is_some() {
            rand_indexes.push(rs.indexes());
        }
        assert_eq!(
            rand_indexes.len(),
            expected,
            "RandStrobes count for w_min={w_min} w_max={w_max}"
        );

        let mut min_indexes = Vec::new();
        let mut ms = MinStrobes::new(seq, 2, l, w_min, w_max)?;
        while ms.next().is_some() {
            min_indexes.push(ms.indexes());
        }
        assert_eq!(
            min_indexes.len(),
            expected,
            "MinStrobes count for w_min={w_min} w_max={w_max}"
        );

        // Every reported strobe position must address a real l-mer.
        for [m1, m2, _] in rand_indexes.iter().chain(min_indexes.iter()) {
            assert!(
                *m1 <= end_hash && *m2 <= end_hash,
                "strobe out of range: {m1},{m2}"
            );
        }
    }
    Ok(())
}

/// When the sequence cannot hold `n` strobes, both iterators must yield nothing,
/// matching the Go reference (whose signed `endIdx` simply goes negative).
/// `MinStrobes` used to underflow its `usize` end index and panic.
#[test]
fn sequence_shorter_than_n_strobes_yields_nothing() -> Result<()> {
    // len 15 < n * l == 20
    let seq: &[u8] = b"ACGTACGTACGTACG";
    assert_eq!(RandStrobes::new(seq, 2, 10, 1, 2)?.count(), 0);
    assert_eq!(MinStrobes::new(seq, 2, 10, 1, 2)?.count(), 0);

    // len 20 < n * l == 21, order 3
    let seq3: &[u8] = b"ACGTACGTACGTACGTACGT";
    assert_eq!(RandStrobes::new(seq3, 3, 7, 1, 3)?.count(), 0);
    assert_eq!(MinStrobes::new(seq3, 3, 7, 1, 3)?.count(), 0);
    Ok(())
}

/// The order-3 `MinStrobes` end-of-sequence branch adds an unreduced running
/// hash to a raw k-mer hash. The Go reference relies on `uint64` wrap-around, so
/// the Rust port must wrap too instead of panicking in debug builds.
#[test]
fn order3_minstrobes_tail_does_not_overflow() -> Result<()> {
    let seq: &[u8] = b"ACGTACGTACGTACGTACGTACGT";
    let count = MinStrobes::with_hasher(seq, 3, 3, 1, 4, &HugeHasher)?.count();
    assert!(
        count > 0,
        "expected the shrinking tail branch to be exercised"
    );
    Ok(())
}

/// `n < 2` and `n > 3` are distinct failures in the Go reference.
#[test]
fn order_validation_matches_reference() {
    let seq: &[u8] = b"ACGATCTGGTACCTAG";
    assert_eq!(
        RandStrobes::new(seq, 1, 3, 1, 4).err(),
        Some(StrobeError::InvalidOrder)
    );
    assert_eq!(
        RandStrobes::new(seq, 4, 3, 1, 4).err(),
        Some(StrobeError::OrderNotSupported)
    );
    assert_eq!(
        MinStrobes::new(seq, 0, 3, 1, 4).err(),
        Some(StrobeError::InvalidOrder)
    );
    assert_eq!(
        MinStrobes::new(seq, 4, 3, 1, 4).err(),
        Some(StrobeError::OrderNotSupported)
    );
}

/// A sequence shorter than one strobe must be rejected rather than underflow
/// when computing the index of the last l-mer.
#[test]
fn sequence_shorter_than_one_strobe_is_rejected() {
    let seq: &[u8] = b"ACGT";
    assert_eq!(
        RandStrobes::new(seq, 2, 10, 1, 2).err(),
        Some(StrobeError::SequenceTooShort)
    );
    assert_eq!(
        MinStrobes::new(seq, 2, 10, 1, 2).err(),
        Some(StrobeError::SequenceTooShort)
    );
}

/// `set_prime` rounds up to a Mersenne-style mask. Very large values must not
/// panic on the `next_power_of_two` overflow; the reference wraps to `u64::MAX`.
#[test]
fn set_prime_accepts_large_values() -> Result<()> {
    let seq: &[u8] = b"ACGATCTGGTACCTAG";
    let mut rs = RandStrobes::new(seq, 2, 3, 1, 4)?;
    rs.set_prime(u64::MAX)?;
    assert!(rs.count() > 0);

    let mut ms = MinStrobes::new(seq, 2, 3, 1, 4)?;
    ms.set_prime(u64::MAX)?;
    assert!(ms.count() > 0);

    assert_eq!(
        RandStrobes::new(seq, 2, 3, 1, 4)?.set_prime(255).err(),
        Some(StrobeError::PrimeNumberTooSmall)
    );
    Ok(())
}
