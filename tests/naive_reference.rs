//! Randomised comparison against a deliberately naive implementation.
//!
//! `tests/go_conformance.rs` pins a handful of parameter sets against the Go
//! reference. This file covers the space between them: for a few thousand
//! random sequences and parameter sets it recomputes every strobemer the slow,
//! obvious way — rescanning each window from scratch — and requires the
//! iterators to agree on every hash and every strobe position. The fast paths
//! (sliding-window minima, the packed argmin) have no counterpart here, so a
//! mistake in them cannot cancel out.

use std::iter::FusedIterator;

use strobemers_rs::{DEFAULT_PRIME_NUMBER, KmerHasher, MinStrobes, RandStrobes, Result, Strobemer};

/// 64-bit FNV-1a over each k-mer, optionally reduced to a handful of distinct
/// values so that windows are full of ties.
struct Fnv1a {
    modulus: Option<u64>,
}

impl KmerHasher for Fnv1a {
    fn hash_all(&self, seq: &[u8], k: usize) -> Result<Vec<u64>> {
        Ok(seq
            .windows(k)
            .map(|w| {
                let mut h: u64 = 14695981039346656037;
                for &b in w {
                    h ^= b as u64;
                    h = h.wrapping_mul(1099511628211);
                }
                match self.modulus {
                    Some(m) => h % m,
                    None => h,
                }
            })
            .collect())
    }
}

/// xorshift64, so the cases are reproducible without a dependency.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    /// Uniform-enough value in `lo..=hi`.
    fn range(&mut self, lo: usize, hi: usize) -> usize {
        lo + (self.next() >> 33) as usize % (hi - lo + 1)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Kind {
    Min,
    Rand,
}

#[derive(Debug, Clone, Copy)]
struct Params {
    kind: Kind,
    n: u8,
    k: usize,
    w_min: usize,
    w_max: usize,
    shrink: bool,
    /// The mask, as stored after `set_prime`.
    mask: u64,
}

/// Earliest position in `lo..=hi` minimising `key`.
fn argmin(lo: usize, hi: usize, key: impl Fn(usize) -> u64) -> usize {
    let mut best = lo;
    for p in lo..=hi {
        if key(p) < key(best) {
            best = p;
        }
    }
    best
}

/// Every strobemer for `hashes` under `p`, computed without any shared state
/// between positions.
fn naive(hashes: &[u64], seq_len: usize, p: Params) -> Vec<Strobemer> {
    let mut out = Vec::new();
    let Some(end_idx) = seq_len.checked_sub(p.n as usize * p.k) else {
        return out;
    };
    let end_hash = seq_len - p.k;
    let masked = |base: u64, pos: usize| base.wrapping_add(hashes[pos]) & p.mask;

    for idx in 0..=end_idx {
        let h1 = hashes[idx];

        // The last window may run off the end: shrink it, or stop.
        let last_lo = idx + (p.n as usize - 2) * p.w_max + p.w_min;
        let mut last_hi = idx + (p.n as usize - 1) * p.w_max;
        let full = last_hi <= end_hash;
        if !full {
            if !p.shrink {
                break;
            }
            last_hi = end_hash;
        }
        if last_lo > last_hi {
            break;
        }

        if p.n == 2 {
            let m2 = match p.kind {
                Kind::Min => argmin(last_lo, last_hi, |q| hashes[q]),
                Kind::Rand => argmin(last_lo, last_hi, |q| masked(h1, q)),
            };
            out.push(Strobemer {
                hash: (h1 >> 1) + hashes[m2] / 3,
                m1: idx,
                m2,
                m3: None,
            });
            continue;
        }

        let (lo, hi) = (idx + p.w_min, idx + p.w_max);
        let m2 = match p.kind {
            Kind::Min => argmin(lo, hi, |q| hashes[q]),
            Kind::Rand => argmin(lo, hi, |q| masked(h1, q)),
        };
        let h2 = h1 / 3 + (hashes[m2] >> 2);
        let m3 = match p.kind {
            // As in the reference, a clipped final window falls back to the
            // masked comparison even for MinStrobes.
            Kind::Min if full => argmin(last_lo, last_hi, |q| hashes[q]),
            _ => argmin(last_lo, last_hi, |q| masked(h2, q)),
        };
        out.push(Strobemer {
            hash: h2 + hashes[m3] / 5,
            m1: idx,
            m2,
            m3: Some(m3),
        });
    }
    out
}

fn assert_fused<I: FusedIterator>(_: &I) {}

#[test]
fn iterators_match_naive_reference() -> Result<()> {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut compared = 0usize;
    let mut items = 0usize;

    for case in 0..4000 {
        let len = rng.range(4, 160);
        // A two-letter alphabet now and then, for long runs of equal k-mers.
        let alphabet: &[u8] = if case % 5 == 0 { b"AC" } else { b"ACGT" };
        let seq: Vec<u8> = (0..len)
            .map(|_| alphabet[rng.range(0, alphabet.len() - 1)])
            .collect();

        let w_min = rng.range(1, 7);
        let p = Params {
            kind: if rng.next() & 1 == 0 {
                Kind::Min
            } else {
                Kind::Rand
            },
            n: rng.range(2, 3) as u8,
            k: rng.range(1, 9),
            w_min,
            w_max: w_min + rng.range(0, 9),
            shrink: rng.next() & 1 == 0,
            mask: DEFAULT_PRIME_NUMBER,
        };
        // (argument to `set_prime`, resulting mask)
        let (prime, mask) = [
            (None, DEFAULT_PRIME_NUMBER),
            (Some(256), 255),
            (Some(1 << 40), (1 << 40) - 1),
            (Some(u64::MAX), u64::MAX),
        ][rng.range(0, 3)];
        let p = Params { mask, ..p };
        let hasher = Fnv1a {
            modulus: [None, Some(5), Some(1 << 30)][rng.range(0, 2)],
        };

        let label = format!("case {case}: {p:?} seq={}", String::from_utf8_lossy(&seq));

        let got: Vec<Strobemer> = match p.kind {
            Kind::Min => {
                let Ok(mut it) = MinStrobes::with_hasher(&seq, p.n, p.k, p.w_min, p.w_max, &hasher)
                else {
                    continue;
                };
                it.set_window_shrink(p.shrink);
                if let Some(q) = prime {
                    it.set_prime(q)?;
                }
                assert_fused(&it);
                let mut it = it.with_positions();
                let got: Vec<_> = it.by_ref().collect();
                assert_eq!(it.next(), None, "{label}");
                got
            }
            Kind::Rand => {
                let Ok(mut it) =
                    RandStrobes::with_hasher(&seq, p.n, p.k, p.w_min, p.w_max, &hasher)
                else {
                    continue;
                };
                it.set_window_shrink(p.shrink);
                if let Some(q) = prime {
                    it.set_prime(q)?;
                }
                assert_fused(&it);
                let mut it = it.with_positions();
                let got: Vec<_> = it.by_ref().collect();
                assert_eq!(it.next(), None, "{label}");
                got
            }
        };

        let hashes = hasher.hash_all(&seq, p.k)?;
        assert_eq!(got, naive(&hashes, seq.len(), p), "{label}");
        compared += 1;
        items += got.len();
    }

    assert!(compared > 2000, "only {compared} cases were valid");
    assert!(items > 50_000, "only {items} strobemers were compared");
    Ok(())
}

/// `with_positions` reports exactly what `next()` followed by `indexes()` does,
/// and keeps the exact length.
#[test]
fn with_positions_agrees_with_indexes() -> Result<()> {
    let seq: &[u8] =
        b"ACGATCTGGTACCTAGTTGACCTGATCGATCGGGCTTATCAGCTTACGGATCAGCTTAGCATCGAAGGCTTACGCT";

    for n in [2u8, 3] {
        let mut plain = MinStrobes::new(seq, n, 5, 3, 7)?;
        let positioned = MinStrobes::new(seq, n, 5, 3, 7)?.with_positions();
        assert_eq!(positioned.len(), plain.len());
        for s in positioned {
            assert_eq!(plain.next(), Some(s.hash));
            let [m1, m2, m3] = plain.indexes();
            assert_eq!((s.m1, s.m2), (m1, m2));
            assert_eq!(s.m3, (n == 3).then_some(m3));
        }
        assert_eq!(plain.next(), None);

        let mut plain = RandStrobes::new(seq, n, 5, 3, 7)?;
        let positioned = RandStrobes::new(seq, n, 5, 3, 7)?.with_positions();
        assert_eq!(positioned.len(), plain.len());
        for s in positioned {
            assert_eq!(plain.next(), Some(s.hash));
            let [m1, m2, m3] = plain.indexes();
            assert_eq!((s.m1, s.m2), (m1, m2));
            assert_eq!(s.m3, (n == 3).then_some(m3));
        }
        assert_eq!(plain.next(), None);
    }

    // Mixing the two styles on one iterator stays consistent.
    let mut it = RandStrobes::new(seq, 3, 5, 3, 7)?;
    it.next();
    let s = it.next_strobemer().expect("second strobemer");
    assert_eq!(s.m1, 1);
    assert_eq!(it.indexes(), [s.m1, s.m2, s.m3.unwrap()]);
    let mut rest = it.with_positions();
    assert_eq!(rest.next().map(|s| s.m1), Some(2));
    assert_eq!(rest.into_inner().index(), Some(2));
    Ok(())
}
