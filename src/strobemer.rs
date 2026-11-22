use core::iter::FusedIterator;

use crate::{MinStrobes, RandStrobes};

/// One strobemer: its hash together with where each strobe starts.
///
/// Produced by [`MinStrobes::with_positions`] and
/// [`RandStrobes::with_positions`]. Positions are zero-based offsets into the
/// sequence the iterator was built from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Strobemer {
    /// Combined hash of all strobes.
    pub hash: u64,
    /// Start of the first strobe.
    pub m1: usize,
    /// Start of the second strobe.
    pub m2: usize,
    /// Start of the third strobe; `None` for order 2.
    pub m3: Option<usize>,
}

/// Iterator adapter yielding [`Strobemer`]s — hash and strobe positions
/// together — instead of bare hashes.
///
/// Created by [`MinStrobes::with_positions`] and
/// [`RandStrobes::with_positions`].
#[derive(Debug, Clone)]
pub struct WithPositions<I> {
    inner: I,
}

impl<I> WithPositions<I> {
    pub(crate) fn new(inner: I) -> Self {
        Self { inner }
    }

    /// Returns the underlying hash iterator, positioned after the last item
    /// yielded so far.
    pub fn into_inner(self) -> I {
        self.inner
    }
}

macro_rules! impl_with_positions {
    ($ty:ident) => {
        impl Iterator for WithPositions<$ty<'_>> {
            type Item = Strobemer;

            #[inline]
            fn next(&mut self) -> Option<Strobemer> {
                self.inner.next_strobemer()
            }

            #[inline]
            fn size_hint(&self) -> (usize, Option<usize>) {
                self.inner.size_hint()
            }
        }

        impl ExactSizeIterator for WithPositions<$ty<'_>> {}
        impl FusedIterator for WithPositions<$ty<'_>> {}
    };
}

impl_with_positions!(MinStrobes);
impl_with_positions!(RandStrobes);
