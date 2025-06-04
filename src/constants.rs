#![allow(clippy::unreadable_literal)]

// Global constants used throughout the product code.

// `DEFAULT_PRIME_NUMBER` – Default prime number (2²⁰ - 1) used for Stöber calculations.
pub const DEFAULT_PRIME_NUMBER: u64 = (1u64 << 20) - 1;

// `ASCII_SIZE` – Number of possible ASCII values (0..255).
pub const ASCII_SIZE: usize = 256;

/// Complement base lookup table. Maps ASCII nucleotide characters to their complementary DNA base:
/// `A ↔ T`, `C ↔ G`. Any other character defaults to `N`.
pub const COMPL_BASES: [u8; ASCII_SIZE] = {
    // Initialize all entries to 'N' (unknown nucleotide).
    let mut tbl = [b'N'; ASCII_SIZE];

    // DNA (uppercase)
    tbl[b'A' as usize] = b'T'; // A → T
    tbl[b'C' as usize] = b'G'; // C → G
    tbl[b'G' as usize] = b'C'; // G → C
    tbl[b'T' as usize] = b'A'; // T → A

    // DNA (lowercase)
    tbl[b'a' as usize] = b'T'; // a → T
    tbl[b'c' as usize] = b'G'; // c → G
    tbl[b'g' as usize] = b'C'; // g → C
    tbl[b't' as usize] = b'A'; // t → A

    // RNA characters map to DNA complement: U (or u) → A
    tbl[b'U' as usize] = b'A';
    tbl[b'u' as usize] = b'A';

    tbl
};

/// 2-bit encoding table for nucleotide sequences. Maps ASCII characters to:
/// A=0, C=1, G=2, T/U=3, any other character=4 (invalid).
pub const SEQ_NT4_TABLE: [u8; ASCII_SIZE] = {
    const INVALID: u8 = 4;
    // Initialize all entries to 4 (invalid).
    let mut t = [INVALID; ASCII_SIZE];

    // Assign valid encodings (uppercase and lowercase).
    t[b'A' as usize] = 0; // A → 00b
    t[b'a' as usize] = 0;

    t[b'C' as usize] = 1; // C → 01b
    t[b'c' as usize] = 1;

    t[b'G' as usize] = 2; // G → 10b
    t[b'g' as usize] = 2;

    t[b'T' as usize] = 3; // T → 11b
    t[b't' as usize] = 3;

    // RNA: U (or u) is treated as T.
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;

    t
};
