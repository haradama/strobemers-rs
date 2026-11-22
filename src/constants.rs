/// Default mask used when combining strobe hashes.
///
/// The reference implementations call this value a "prime", and the setters
/// are named accordingly, but it is used as a bit mask — candidates are compared
/// by `(h1 + h2) & mask` — and it is not a prime number.
pub const DEFAULT_PRIME_NUMBER: u64 = (1u64 << 20) - 1;
