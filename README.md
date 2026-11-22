# strobemers-rs

[<img alt="github" src="https://img.shields.io/badge/github-haradama/strobemers__rs-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/haradama/strobemers-rs)
[<img alt="crates.io" src="https://img.shields.io/crates/v/strobemers-rs.svg?style=for-the-badge&color=fc8d62&logo=rust" height="20">](https://crates.io/crates/strobemers-rs)
[<img alt="docs.rs" src="https://img.shields.io/badge/docs.rs-strobemers__rs-66c2a5?style=for-the-badge&labelColor=555555&logo=docs.rs" height="20">](https://docs.rs/strobemers-rs)
[<img alt="build status" src="https://img.shields.io/github/actions/workflow/status/haradama/strobemers-rs/rust.yml?branch=main&style=for-the-badge" height="20">](https://github.com/haradama/strobemers-rs/actions)

Rust implementation of **strobemer** generation (MinStrobes, RandStrobes).  
Iterators produce 64-bit strobemer hashes, optionally with strobe positions, over DNA sequences.

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

Sequences must consist of `A`, `C`, `G`, `T` only (either case). The default
hasher has no value for a `k`-mer containing any other base — `N`, IUPAC codes,
or `U`, so convert RNA to `T` first — and such sequences are rejected
with `StrobeError::UnsupportedBase` rather than hashed with gaps — see
[Sequences containing `N`](#sequences-containing-n).

## Installation

Add to your `Cargo.toml` via command line.

```shell
cargo add strobemers-rs
```

## Quick Start

```rust
use strobemers_rs::{MinStrobes, RandStrobes};

fn main() -> strobemers_rs::Result<()> {
    // Example DNA sequence
    let seq = b"ACGATCTGGTACCTAG";
    // k-mer length
    let k = 3;
    // window offsets
    let w_min = 3;
    let w_max = 5;

    // MinStrobes order-2
    // Produces one 64-bit hash per position, selecting the minimum-hash strobe in [i+w_min .. i+w_max].
    for h in MinStrobes::new(seq, 2, k, w_min, w_max)? {
        println!("MinStrobes(order=2) hash: {h}");
    }

    // MinStrobes order-3
    // Two-level selection: pick m2 in [i+w_min .. i+w_max], then m3 in
    // [i+w_max+w_min .. i+2*w_max]. Both windows are anchored on the position of
    // m1, not on where m2 happened to land.
    //
    // `with_positions()` yields the start of each strobe along with the hash.
    let min3 = MinStrobes::new(seq, 3, k, w_min, w_max)?;
    for (i, s) in min3.with_positions().take(5).enumerate() {
        println!("MinStrobes(order=3)[{i}] m1={} m2={} m3={:?} -> hash={}",
            s.m1, s.m2, s.m3, s.hash);
    }

    // RandStrobes order-2
    // Similar to MinStrobes, but the second strobe is chosen by a pseudo-random rule:
    for h in RandStrobes::new(seq, 2, k, w_min, w_max)? {
        println!("RandStrobes(order=2) hash: {h}");
    }

    // RandStrobes order-3
    let rand3 = RandStrobes::new(seq, 3, k, w_min, w_max)?;
    for (i, s) in rand3.with_positions().take(5).enumerate() {
        println!("RandStrobes(order=3)[{i}] m1={} m2={} m3={:?} -> hash={}",
            s.m1, s.m2, s.m3, s.hash);
    }

    Ok(())
}
```

## Sequences containing `N`

A strobemer cannot span an ambiguous base, so the default hasher rejects the
whole sequence with `StrobeError::UnsupportedBase`, which reports where the
first such base is. `acgt_runs` splits a sequence into the stretches between
them, each with its offset:

```rust
use strobemers_rs::{RandStrobes, StrobeError, acgt_runs};

fn main() -> strobemers_rs::Result<()> {
    let seq = b"ACGATCTGGTACCTAGNNNNTTGACCTGATCGATCGGGCTNAC";

    for (offset, run) in acgt_runs(seq) {
        let strobes = match RandStrobes::new(run, 2, 3, 3, 5) {
            Ok(it) => it,
            // This run is too short for the chosen parameters.
            Err(StrobeError::SequenceTooShort) => continue,
            Err(e) => return Err(e),
        };
        for s in strobes.with_positions() {
            println!("m1={} m2={} hash={}", offset + s.m1, offset + s.m2, s.hash);
        }
    }
    Ok(())
}
```

## Performance notes

Iterators built with `new` borrow the sequence and hash it in chunks as they
advance, so their memory use does not grow with the sequence length. `MinStrobes`
also needs the minimum of a sliding window: for short sequences it hashes
everything up front and precomputes the minima as a table, which is faster
while the table stays cache-resident; for longer ones it computes them online,
keeping only `O(w_max)` entries.

Iterators built with `with_hasher` call the custom hasher once, for the whole
sequence, and own the resulting hashes — one 64-bit value per k-mer.

`RandStrobes` picks its next strobe by minimising `(h1 + h2) & prime` over the
window. With the default prime that reduces to a branchless 64-bit `min`, which
vectorises into `vpminuq` when the target has AVX-512, so building with
`RUSTFLAGS="-C target-cpu=native"` (or `-C target-cpu=x86-64-v4`) makes it
noticeably faster there. `MinStrobes` is unaffected, being bound by memory
traffic rather than by that loop.

Run `cargo bench` to measure on your own hardware.

Both iterators are [`ExactSizeIterator`], so `collect()` allocates exactly once.

[`ExactSizeIterator`]: https://doc.rust-lang.org/std/iter/trait.ExactSizeIterator.html

## Using a Custom Hash Function

By default, strobemers-rs uses [nthash-rs](https://github.com/haradama/nthash-rs) for k-mer hashing.
However, you can inject your own hash function by implementing the `KmerHasher` trait and passing it via the `with_hasher` method. See [the example](https://github.com/haradama/strobemers-rs/blob/main/examples/custom_hash.rs) for more details.

## License

This project is MIT‑licensed (see [LICENSE](https://github.com/haradama/strobemers-rs/blob/main/LICENSE)).
