use crate::constants::{COMPL_BASES, SEQ_NT4_TABLE};

/// Rounds up the given value `x` to the next power of two.
///
/// Examples:
/// - `roundup64(5)` -> 8
/// - `roundup64(8)` -> 8
///
/// # Arguments
///
/// * `x` – An unsigned 64-bit integer to round up.
///
/// # Returns
///
/// * The smallest power of two greater than or equal to `x`.
#[inline(always)]
pub const fn roundup64(x: u64) -> u64 {
    // `next_power_of_two` panics on overflow; the Go reference wraps to 0 for
    // `x > 2^63`, which `SetPrime` then turns into `u64::MAX`. Match that.
    match x.checked_next_power_of_two() {
        Some(v) => v,
        None => 0,
    }
}

/// Returns the complementary DNA/RNA base for the given ASCII byte.
///
/// Looks up the byte in the `COMPL_BASES` table, which maps:
/// `A ↔ T`, `C ↔ G` (uppercase and lowercase), `U/u → A`, and all others to `N`.
///
/// # Arguments
///
/// * `b` – An ASCII byte representing a nucleotide.
///
/// # Returns
///
/// * The ASCII byte for the complementary base, or `b'N'` if outside A/C/G/T/U.
#[inline(always)]
pub const fn complement(b: u8) -> u8 {
    COMPL_BASES[b as usize]
}

/// Encodes a nucleotide ASCII byte into its 2-bit code (0‒3), or 4 for invalid.
///
/// Uses the `SEQ_NT4_TABLE`, which assigns:
/// - A/a → 0
/// - C/c → 1
/// - G/g → 2
/// - T/t, U/u → 3
/// - Any other ASCII byte → 4
///
/// # Arguments
///
/// * `b` – An ASCII byte representing a nucleotide.
///
/// # Returns
///
/// * A 2-bit encoding (0..=3) for valid nucleotides, or 4 for any other byte.
#[inline(always)]
pub const fn nt4(b: u8) -> u8 {
    SEQ_NT4_TABLE[b as usize]
}

/// Largest `m1` index at which a strobemer iterator can still produce an item.
///
/// Both iterators stop at the first position whose search window runs past the
/// last `l`-mer, and that bound only tightens as the index grows, so this makes
/// the remaining length exact:
///
/// - order 2 needs `idx + w_min ≤ end_hash` when shrinking, `idx + w_max ≤ end_hash` otherwise
/// - order 3 needs `idx + w_max + w_min ≤ end_hash` when shrinking, `idx + 2·w_max ≤ end_hash` otherwise
///
/// # Arguments
///
/// * `n` – Strobemer order (2 or 3).
/// * `shrink` – Whether terminal windows may be shortened.
/// * `w_min`, `w_max` – Window offsets.
/// * `end_hash` – Index of the last `l`-mer (`seq.len() - l`).
///
/// # Returns
///
/// * `Some(idx)` – The last usable starting index.
/// * `None` – No position qualifies, or `n` is out of range.
#[inline]
pub(crate) const fn last_start_index(
    n: u8,
    shrink: bool,
    w_min: usize,
    w_max: usize,
    end_hash: usize,
) -> Option<usize> {
    match (n, shrink) {
        (2, true) => end_hash.checked_sub(w_min),
        (2, false) => end_hash.checked_sub(w_max),
        (3, true) => end_hash.checked_sub(w_max + w_min),
        (3, false) => end_hash.checked_sub(w_max << 1),
        _ => None,
    }
}

/// Validates parameters for strobemer construction and returns early on error.
///
/// This macro is intended to be invoked at the start of constructors or functions
/// that require:
/// - A non-empty, ASCII-only sequence slice (`$seq`)
/// - An order (`$n`) of either 2 or 3
/// - A strobe length (`$l`) between 1 and 64
/// - Window offsets (`$w_min`, `$w_max`) where both are > 0 and `w_min ≤ w_max`
/// - Sequence length sufficient to accommodate `(n - 1)` windows of size `(w_max + 1)`
/// - Sequence length of at least `l`, so that one strobe fits
///
/// Returns the corresponding `StrobeError` on any validation failure:
/// - `InvalidSequence` if the sequence is empty
/// - `InvalidOrder` if `n < 2`
/// - `OrderNotSupported` if `n > 3`
/// - `StrobeLengthTooSmall` if `l` is outside [1..=64]
/// - `InvalidWindowOffsets` if `w_min` or `w_max` are zero or `w_min > w_max`
/// - `SequenceTooShort` if `seq.len()` is too small for the given parameters
///
/// # Example
///
/// ```ignore
/// validate_params!(seq, n, l, w_min, w_max);
/// ```
macro_rules! validate_params {
    ($seq:expr, $n:expr, $l:expr, $w_min:expr, $w_max:expr) => {{
        // Sequence must be non-empty
        if $seq.is_empty() || !$seq.is_ascii() {
            return Err(StrobeError::InvalidSequence);
        }
        // Order must be exactly 2 or 3. The two failure modes are kept distinct,
        // matching the Go reference (ErrInvalidOrder vs ErrOrderNotSupported).
        if $n < 2 {
            return Err(StrobeError::InvalidOrder);
        }
        if $n > 3 {
            return Err(StrobeError::OrderNotSupported);
        }
        // Strobe length must be between 1 and 64 inclusive
        if !(1..=64).contains(&$l) {
            return Err(StrobeError::StrobeLengthTooSmall);
        }
        // Window offsets must be greater than zero and w_min ≤ w_max
        if $w_min == 0 || $w_max == 0 || $w_min > $w_max {
            return Err(StrobeError::InvalidWindowOffsets);
        }
        // Sequence must be long enough to fit (n − 1) windows of size (w_max + 1),
        // and to hold at least one strobe. The latter guarantees that
        // `seq.len() - l` (the index of the last l-mer) cannot underflow.
        if $seq.len() < ($n as usize - 1) * ($w_max + 1) || $seq.len() < $l {
            return Err(StrobeError::SequenceTooShort);
        }
    }};
}
