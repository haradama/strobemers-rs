use crate::{Result, StrobeError};

/// Rounds up the given value `x` to the next power of two.
///
/// Examples:
/// - `roundup64(5)` -> 8
/// - `roundup64(8)` -> 8
#[inline(always)]
pub(crate) const fn roundup64(x: u64) -> u64 {
    // `next_power_of_two` panics on overflow; the Go reference wraps to 0 for
    // `x > 2^63`, which `SetPrime` then turns into `u64::MAX`. Match that.
    match x.checked_next_power_of_two() {
        Some(v) => v,
        None => 0,
    }
}

/// Whether the default ntHash-based hasher can hash `b`: `A`, `C`, `G` or `T`,
/// in either case. `nthash-rs` has no seed for `U`.
#[inline(always)]
pub(crate) const fn is_supported_base(b: u8) -> bool {
    // Setting bit 5 folds the upper-case letters onto the lower-case ones and
    // maps no other byte onto them.
    matches!(b | 0x20, b'a' | b'c' | b'g' | b't')
}

/// Splits `seq` into its maximal runs of unambiguous bases.
///
/// The default hasher rejects any sequence containing a base other than `A`,
/// `C`, `G` or `T` (in either case) with
/// [`StrobeError::UnsupportedBase`], because a strobemer cannot span such a
/// base. Real assemblies contain stretches of `N`, so this helper yields each
/// stretch between them, together with its offset in `seq`, ready to be handed
/// to [`MinStrobes`](crate::MinStrobes) or [`RandStrobes`](crate::RandStrobes).
///
/// Runs too short for the chosen parameters make the constructors return
/// [`StrobeError::SequenceTooShort`]; skip those.
///
/// # Example
/// ```
/// use strobemers_rs::{MinStrobes, StrobeError, acgt_runs};
///
/// let seq = b"ACGATCTGGTACCTAGNNNNTTGACCTGATCGATCGGGCTNAC";
/// let mut positions = Vec::new();
/// for (offset, run) in acgt_runs(seq) {
///     let strobes = match MinStrobes::new(run, 2, 3, 3, 5) {
///         Ok(it) => it,
///         Err(StrobeError::SequenceTooShort) => continue,
///         Err(e) => return Err(e),
///     };
///     for s in strobes.with_positions() {
///         // Positions are relative to the run; add the offset back.
///         positions.push((offset + s.m1, offset + s.m2));
///     }
/// }
/// assert!(positions.iter().any(|&(m1, _)| m1 >= 20));
/// # Ok::<(), StrobeError>(())
/// ```
pub fn acgt_runs(seq: &[u8]) -> AcgtRuns<'_> {
    AcgtRuns { seq, pos: 0 }
}

/// Iterator over the unambiguous runs of a sequence, created by [`acgt_runs`].
///
/// Yields `(offset, run)` pairs, where `run == &seq[offset..offset + run.len()]`.
#[derive(Debug, Clone)]
pub struct AcgtRuns<'a> {
    seq: &'a [u8],
    pos: usize,
}

impl<'a> Iterator for AcgtRuns<'a> {
    type Item = (usize, &'a [u8]);

    fn next(&mut self) -> Option<Self::Item> {
        let rest = &self.seq[self.pos..];
        let skipped = rest.iter().position(|&b| is_supported_base(b))?;
        let start = self.pos + skipped;
        let rest = &self.seq[start..];
        let len = rest
            .iter()
            .position(|&b| !is_supported_base(b))
            .unwrap_or(rest.len());
        self.pos = start + len;
        Some((start, &rest[..len]))
    }
}

impl core::iter::FusedIterator for AcgtRuns<'_> {}

/// Validates the parameters shared by both strobemer constructors.
///
/// Requires:
/// - a non-empty, ASCII-only sequence
/// - an order `n` of either 2 or 3
/// - a strobe length `l` in `1..=64`
/// - window offsets that are both > 0 with `w_min ≤ w_max`
/// - a sequence long enough for `(n - 1)` windows of size `(w_max + 1)`, and
///   for at least one strobe
///
/// On success every index the iterators compute (`idx + 2 * w_max` at most) is
/// known to fit in a `usize`, because `w_max < seq.len()`.
///
/// # Errors
/// - `InvalidSequence` if the sequence is empty or not ASCII
/// - `InvalidOrder` if `n < 2`
/// - `OrderNotSupported` if `n > 3`
/// - `InvalidStrobeLength` if `l` is outside `1..=64`
/// - `InvalidWindowOffsets` if `w_min` or `w_max` are zero or `w_min > w_max`
/// - `SequenceTooShort` if `seq.len()` is too small for the given parameters
pub(crate) fn validate_params(
    seq: &[u8],
    n: u8,
    l: usize,
    w_min: usize,
    w_max: usize,
) -> Result<()> {
    if seq.is_empty() || !seq.is_ascii() {
        return Err(StrobeError::InvalidSequence);
    }
    // Order must be exactly 2 or 3. The two failure modes are kept distinct,
    // matching the Go reference (ErrInvalidOrder vs ErrOrderNotSupported).
    if n < 2 {
        return Err(StrobeError::InvalidOrder);
    }
    if n > 3 {
        return Err(StrobeError::OrderNotSupported);
    }
    if !(1..=64).contains(&l) {
        return Err(StrobeError::InvalidStrobeLength);
    }
    if w_min == 0 || w_max == 0 || w_min > w_max {
        return Err(StrobeError::InvalidWindowOffsets);
    }
    // Sequence must be long enough to fit (n − 1) windows of size (w_max + 1),
    // and to hold at least one strobe. The latter guarantees that
    // `seq.len() - l` (the index of the last l-mer) cannot underflow.
    //
    // A window span that overflows `usize` is longer than any sequence.
    let span = w_max
        .checked_add(1)
        .and_then(|w| w.checked_mul(n as usize - 1))
        .ok_or(StrobeError::SequenceTooShort)?;
    if seq.len() < span || seq.len() < l {
        return Err(StrobeError::SequenceTooShort);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn acgt_runs_splits_on_ambiguous_bases() {
        let runs: Vec<_> = acgt_runs(b"NNACGTNacgtUnT").collect();
        assert_eq!(
            runs,
            vec![(2, &b"ACGT"[..]), (7, &b"acgt"[..]), (13, &b"T"[..])]
        );
        assert_eq!(acgt_runs(b"").count(), 0);
        assert_eq!(acgt_runs(b"NNN").count(), 0);
        assert_eq!(
            acgt_runs(b"ACGT").collect::<Vec<_>>(),
            vec![(0, &b"ACGT"[..])]
        );
    }

    #[test]
    fn huge_window_is_an_error_not_a_panic() {
        let seq = b"ACGATCTGGTACCTAG";
        for n in [2u8, 3] {
            for w_max in [usize::MAX, usize::MAX / 2, usize::MAX / 2 + 1] {
                assert_eq!(
                    validate_params(seq, n, 3, 1, w_max),
                    Err(StrobeError::SequenceTooShort),
                    "n={n} w_max={w_max}"
                );
            }
        }
    }
}
