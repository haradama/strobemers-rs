use crate::constants::{COMPL_BASES, SEQ_NT4_TABLE};

/// Rounds up the given value `x` to the next power of two.
///
/// Internally:
/// 1. Ensures `x` is odd by OR’ing with 1.
/// 2. Subtracts 1 to handle exact powers-of-two correctly.
/// 3. Calls `.next_power_of_two()` to obtain the next power-of-two boundary.
///
/// Examples:
/// - `roundup64(5)` → 8
/// - `roundup64(8)` → 8
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
    x.next_power_of_two()
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

/// Validates parameters for strobemer construction and returns early on error.
///
/// This macro is intended to be invoked at the start of constructors or functions
/// that require:
/// - A non-empty, ASCII-only sequence slice (`$seq`)
/// - An order (`$n`) of either 2 or 3
/// - A strobe length (`$l`) between 1 and 64
/// - Window offsets (`$w_min`, `$w_max`) where both are > 0 and `w_min ≤ w_max`
/// - Sequence length sufficient to accommodate `(n - 1)` windows of size `(w_max + 1)`
///
/// Returns the corresponding `StrobeError` on any validation failure:
/// - `InvalidSequence` if the sequence is empty
/// - `OrderNotSupported` if `n` is not 2 or 3
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
        if $seq.is_empty() {
            return Err(StrobeError::InvalidSequence);
        }
        // Order must be exactly 2 or 3
        if !matches!($n, 2 | 3) {
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
        // Sequence must be long enough to fit (n − 1) windows of size (w_max + 1)
        if $seq.len() < ($n as usize - 1) * ($w_max + 1) {
            return Err(StrobeError::SequenceTooShort);
        }
    }};
}
