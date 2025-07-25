mod constants;
#[macro_use]
mod util;
mod hashes;
mod minstrobes;
mod randstrobes;

pub use constants::*;
pub use hashes::{KmerHasher, compute_min_hashes};
pub use minstrobes::MinStrobes;
pub use randstrobes::RandStrobes;
pub use util::*;

use nthash_rs::NtHashError;

/// Common `Result` type for all library operations, using `StrobeError` for errors.
pub type Result<T, E = StrobeError> = core::result::Result<T, E>;

/// Error variants for strobemer generation and related operations.
///
/// This `enum` covers validation failures, invalid parameters, and errors
/// propagated from the `nthash-rs` crate.
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum StrobeError {
    /// Thrown when the requested strobemer order is not supported.
    /// Only orders 2 and 3 are allowed.
    #[error("strobemer order not supported (must be 2 or 3)")]
    OrderNotSupported,

    /// Thrown when the requested strobemer order is less than 2.
    #[error("strobemer order must be ≥ 2")]
    InvalidOrder,

    /// Thrown when the input sequence is empty or contains non-ASCII characters.
    #[error("invalid DNA sequence (empty or contains non-ASCII)")]
    InvalidSequence,

    /// Thrown when the sequence is too short to generate any strobemers
    /// given the provided parameters.
    #[error("sequence too short for given parameters")]
    SequenceTooShort,

    /// Thrown when the strobe (k-mer) length `l` is less than 1 or greater than 64.
    #[error("strobe length (l) must be ≥ 1 and ≤ 64")]
    StrobeLengthTooSmall,

    /// Thrown when window offsets are invalid (must be > 0 and `w_min ≤ w_max`).
    #[error("window offsets must be > 0 and w_min ≤ w_max")]
    InvalidWindowOffsets,

    /// Indicates that the precomputed k-mer hash values (via `nthash-rs`) were incomplete.
    /// This should not happen under normal circumstances.
    #[error("incomplete pre-computed hash values (nthash)")]
    IncompleteHashValues,

    /// Thrown when the provided prime number is too small (minimum allowed is 256).
    #[error("prime number too small (must be ≥ 256)")]
    PrimeNumberTooSmall,

    /// Wraps errors originating from the `nthash-rs` crate.
    #[error(transparent)]
    NtHashError(#[from] NtHashError),
}
