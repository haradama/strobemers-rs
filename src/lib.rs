#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod constants;
mod hashes;
#[macro_use]
mod state;
mod minstrobes;
mod randstrobes;
mod strobemer;
mod util;

pub use constants::DEFAULT_PRIME_NUMBER;
pub use hashes::{KmerHasher, NtHash64};
pub use minstrobes::MinStrobes;
pub use randstrobes::RandStrobes;
pub use strobemer::{Strobemer, WithPositions};
pub use util::{AcgtRuns, acgt_runs};

use nthash_rs::NtHashError;

/// Common `Result` type for all library operations, using `StrobeError` for errors.
pub type Result<T, E = StrobeError> = core::result::Result<T, E>;

/// Error variants for strobemer generation and related operations.
///
/// This `enum` covers validation failures, invalid parameters, and errors
/// propagated from the `nthash-rs` crate.
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
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

    /// Thrown by the default hasher when the sequence contains a base other
    /// than `A`, `C`, `G` or `T` (in either case), such as `N` or `U`.
    ///
    /// ntHash has no value for a k-mer containing such a base, so the sequence
    /// cannot be hashed position by position. Use [`acgt_runs`] to process the
    /// stretches between ambiguous bases separately.
    #[error("unsupported base {:?} at position {position}", char::from(*.byte))]
    UnsupportedBase {
        /// Zero-based position of the first unsupported base.
        position: usize,
        /// The offending byte.
        byte: u8,
    },

    /// Thrown when the sequence is too short to generate any strobemers
    /// given the provided parameters.
    #[error("sequence too short for given parameters")]
    SequenceTooShort,

    /// Thrown when the strobe (k-mer) length `l` is outside `1..=64`.
    #[error("strobe length (l) must be ≥ 1 and ≤ 64")]
    InvalidStrobeLength,

    /// Thrown when window offsets are invalid (must be > 0 and `w_min ≤ w_max`).
    #[error("window offsets must be > 0 and w_min ≤ w_max")]
    InvalidWindowOffsets,

    /// Thrown when a [`KmerHasher`] did not return exactly one hash per k-mer
    /// (`seq.len() - k + 1` values).
    #[error("k-mer hasher returned the wrong number of hash values")]
    IncompleteHashValues,

    /// Thrown when the provided prime number is too small (minimum allowed is 256).
    #[error("prime number too small (must be ≥ 256)")]
    PrimeNumberTooSmall,

    /// Wraps errors originating from the `nthash-rs` crate.
    #[error(transparent)]
    NtHashError(#[from] NtHashError),
}
