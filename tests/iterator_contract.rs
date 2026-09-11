//! Checks the `Iterator` contract that `MinStrobes` and `RandStrobes` promise.
//!
//! Both implement [`ExactSizeIterator`], so `size_hint` must report the exact
//! remaining length at every point of the iteration — including before the first
//! `next()` and after the last one — for every combination of order, strobe
//! length, window offsets and window-shrinking.

use strobemers_rs::{MinStrobes, RandStrobes};

const SEQS: [&str; 5] = [
    "ACGATCTGGTACCTAG",
    "ACGTACGTACGTACGTACGTACGT",
    "TTGACCTGATCGATCGGGCTTATCAGCTTACGGATCAGCTTAGCATCG",
    "AAAAAAAAAAAAAAAAAAAA",
    "ACGATCTGGTACCTAGTTGACCTGATCGATCGGGCTTATCAGCTTACGGATCAGCTTAGCATCGAAGGCTTACGCT",
];

const WINDOWS: [(usize, usize); 8] = [
    (1, 1),
    (1, 4),
    (2, 5),
    (3, 5),
    (4, 4),
    (1, 10),
    (5, 9),
    (7, 7),
];

/// Walks `it` to exhaustion, asserting that `size_hint` stays exact throughout.
fn assert_exact<I: Iterator<Item = u64>>(mut it: I, label: &str) {
    let total = it.size_hint();
    assert_eq!(
        total.1,
        Some(total.0),
        "upper bound must be present and equal to the lower bound ({label})"
    );

    let mut seen = 0usize;
    loop {
        let hint = it.size_hint();
        assert_eq!(
            hint.1,
            Some(hint.0),
            "upper bound diverged after {seen} items ({label})"
        );
        assert_eq!(
            hint.0 + seen,
            total.0,
            "remaining + consumed must stay constant after {seen} items ({label})"
        );
        if it.next().is_none() {
            break;
        }
        seen += 1;
    }

    assert_eq!(
        seen, total.0,
        "size_hint promised {} items, iterator produced {seen} ({label})",
        total.0
    );
}

#[test]
fn size_hint_is_exact() {
    let mut checked = 0usize;
    for seq in SEQS {
        for n in [2u8, 3] {
            for l in [2usize, 3, 5, 7] {
                for (w_min, w_max) in WINDOWS {
                    for shrink in [true, false] {
                        let label = format!(
                            "len={} n={n} l={l} w=({w_min},{w_max}) shrink={shrink}",
                            seq.len()
                        );

                        if let Ok(mut it) = RandStrobes::new(seq.as_bytes(), n, l, w_min, w_max) {
                            it.set_window_shrink(shrink);
                            assert_exact(it, &format!("RandStrobes {label}"));
                            checked += 1;
                        }
                        if let Ok(mut it) = MinStrobes::new(seq.as_bytes(), n, l, w_min, w_max) {
                            it.set_window_shrink(shrink);
                            assert_exact(it, &format!("MinStrobes {label}"));
                            checked += 1;
                        }
                    }
                }
            }
        }
    }
    assert!(
        checked > 1000,
        "expected broad coverage, only checked {checked}"
    );
}

/// `ExactSizeIterator::len` must agree with what the iterator actually yields.
#[test]
fn exact_size_len_matches_output() {
    let seq: &[u8] =
        b"ACGATCTGGTACCTAGTTGACCTGATCGATCGGGCTTATCAGCTTACGGATCAGCTTAGCATCGAAGGCTTACGCT";
    for n in [2u8, 3] {
        for (w_min, w_max) in [(1usize, 4usize), (3, 5), (5, 9)] {
            let rs = RandStrobes::new(seq, n, 5, w_min, w_max).unwrap();
            let len = rs.len();
            assert_eq!(len, rs.count(), "RandStrobes n={n} w=({w_min},{w_max})");

            let ms = MinStrobes::new(seq, n, 5, w_min, w_max).unwrap();
            let len = ms.len();
            assert_eq!(len, ms.count(), "MinStrobes n={n} w=({w_min},{w_max})");
        }
    }
}

/// An iterator that cannot produce anything must report zero rather than
/// over-promising to `collect()`.
#[test]
fn size_hint_zero_when_nothing_fits() {
    // Sequence shorter than `n * l`.
    let seq: &[u8] = b"ACGTACGTACGTACG";
    assert_eq!(
        RandStrobes::new(seq, 2, 10, 1, 2).unwrap().size_hint(),
        (0, Some(0))
    );
    assert_eq!(
        MinStrobes::new(seq, 2, 10, 1, 2).unwrap().size_hint(),
        (0, Some(0))
    );

    // Order-3 needs `idx + 2 * w_max <= end_hash`; here `end_hash` is 9 and
    // `2 * w_max` is 10, so without shrinking no position ever qualifies.
    let short: &[u8] = b"ACGATCTGGTAC"; // len 12
    let mut rs = RandStrobes::new(short, 3, 3, 2, 5).unwrap();
    rs.set_window_shrink(false);
    assert_eq!(rs.size_hint(), (0, Some(0)));
    assert_eq!(rs.count(), 0);

    let mut ms = MinStrobes::new(short, 3, 3, 2, 5).unwrap();
    ms.set_window_shrink(false);
    assert_eq!(ms.size_hint(), (0, Some(0)));
    assert_eq!(ms.count(), 0);

    // ...while shrinking (`idx + w_max + w_min <= end_hash`) still yields items.
    let rs = RandStrobes::new(short, 3, 3, 2, 5).unwrap();
    assert_eq!(rs.len(), 3);
    assert_eq!(rs.count(), 3);
}
