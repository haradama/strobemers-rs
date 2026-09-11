# strobemers-rs

[<img alt="github" src="https://img.shields.io/badge/github-haradama/strobemers__rs-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/haradama/strobemers-rs)
[<img alt="crates.io" src="https://img.shields.io/crates/v/strobemers-rs.svg?style=for-the-badge&color=fc8d62&logo=rust" height="20">](https://crates.io/crates/strobemers-rs)
[<img alt="docs.rs" src="https://img.shields.io/badge/docs.rs-strobemers__rs-66c2a5?style=for-the-badge&labelColor=555555&logo=docs.rs" height="20">](https://docs.rs/strobemers-rs)
[<img alt="build status" src="https://img.shields.io/github/actions/workflow/status/haradama/strobemers-rs/rust.yml?branch=main&style=for-the-badge" height="20">](https://github.com/haradama/strobemers-rs/actions)

Rust implementation of **strobemer** generation (MinStrobes, RandStrobes).  
Streaming iterators produce 64-bit strobemer hashes over DNA/RNA sequences.

> **Background**  
>
> - The original strobemer concept and reference implementation live at [ksahlin/strobemers](https://github.com/ksahlin/strobemers).  
> - This crate’s design was guided by the Go port found at [shenwei356/strobemers](https://github.com/shenwei356/strobemers), but the code here is a complete Rust‐native rewrite (zero `unsafe`).
> - Strobe selection is verified against that Go port value-by-value in `tests/go_conformance.rs`. Hash *values* differ between the two, because they use different ntHash libraries; the tests therefore drive both with the same injected k-mer hasher.
> - One deliberate deviation: when `w_min > (n-1) * k`, the search window for the last few positions starts past the end of the sequence. The Go code silently reuses the previous strobe there; this crate ends iteration instead, so every reported strobe position is a real k-mer.

## Overview

Strobemers are an alternative to contiguous k-mer hashing designed to improve sensitivity in long-read alignment and genome comparison.  
This crate provides:

- **MinStrobes (order 2/3)**  
  Choose subsequent strobes by selecting the minimum hash within a sliding window.

- **RandStrobes (order 2/3)**  
  Choose subsequent strobes by a pseudo-random but reproducible rule derived from hash values.

Sequences must consist of `A`, `C`, `G`, `T`/`U` only. The default hasher skips
`k`-mer windows containing any other base, which would break the mapping from
hash index to sequence position, so such sequences are rejected with
`StrobeError::IncompleteHashValues` (as in the Go reference).

## Installation

Add to your `Cargo.toml` via command line.

```shell
cargo add strobemers-rs
```

## Quick Start

```rust
use strobemers_rs::{MinStrobes, RandStrobes};

fn main() -> anyhow::Result<()> {
    // Example DNA sequence
    let seq = b"ACGATCTGGTACCTAG";
    // k-mer length
    let k = 3;
    // window offsets
    let w_min = 3;
    let w_max = 5;

    // MinStrobes order-2
    // Produces one 64-bit hash per position, selecting the minimum-hash strobe in [i+w_min .. i+w_max].
    let mut min2 = MinStrobes::new(seq, 2, k, w_min, w_max)?;
    while let Some(h) = min2.next() {
        println!("MinStrobes(order=2) hash: {}", h);
    }

    // MinStrobes order-3
    // Two-level selection: pick m2 in [i+w_min .. i+w_max], then m3 in
    // [i+w_max+w_min .. i+2*w_max]. Both windows are anchored on the position of
    // m1, not on where m2 happened to land.
    let mut min3 = MinStrobes::new(seq, 3, k, w_min, w_max)?;
    for (i, h) in min3.take(5).enumerate() {
        let [m1, m2, m3_pos] = min3.indexes();
        println!("MinStrobes(order=3)[{}] m1={} m2={} m3={} -> hash={}",
            i, m1, m2, m3_pos, h);
    }

    // RandStrobes order-2
    // Similar to MinStrobes, but the second strobe is chosen by a pseudo-random rule:
    let mut rand2 = RandStrobes::new(seq, 2, k, w_min, w_max)?;
    while let Some(h) = rand2.next() {
        println!("RandStrobes(order=2) hash: {}", h);
    }

    // RandStrobes order-3
    let mut rand3 = RandStrobes::new(seq, 3, k, w_min, w_max)?;
    for (i, h) in rand3.take(5).enumerate() {
        let [m1, m2, m3_pos] = rand3.indexes();
        println!("RandStrobes(order=3)[{}] m1={} m2={} m3={} -> hash={}",
            i, m1, m2, m3_pos, h);
    }

    Ok(())
}
```

## Performance notes

Both iterators precompute one 64-bit hash per k-mer up front, so peak scratch
memory is 8 bytes per base for `RandStrobes` and 24 for `MinStrobes`. On a
1 Mbp sequence with `k = 31, w_min = 20, w_max = 50` a full pass costs roughly
24 ms for `MinStrobes` and 31–55 ms for `RandStrobes` (Xeon Platinum 8480+).

`RandStrobes` picks its next strobe by minimising `(h1 + h2) & prime` over the
window. With the default prime that reduces to a branchless 64-bit `min`, which
vectorises into `vpminuq` when the target has AVX-512 — building with
`RUSTFLAGS="-C target-cpu=native"` (or `-C target-cpu=x86-64-v4`) is worth about
another 1.4x there. `MinStrobes` is unaffected, being bound by memory traffic
rather than by that loop.

Both iterators are [`ExactSizeIterator`], so `collect()` allocates exactly once.

[`ExactSizeIterator`]: https://doc.rust-lang.org/std/iter/trait.ExactSizeIterator.html

## Using a Custom Hash Function

By default, strobemers-rs uses [nthash-rs](https://github.com/haradama/nthash-rs) for k-mer hashing.
However, you can inject your own hash function by implementing the `KmerHasher` trait and passing it via the with_hasher method. See [the example](./examples/custom_hash.rs) for more details.

## License

This project is MIT‑licensed (see [LICENSE](LICENSE)).
