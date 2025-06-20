# strobemers-rs

Rust implementation of **strobemer** generation (MinStrobes, RandStrobes).  
Streaming iterators produce 64-bit strobemer hashes over DNA/RNA sequences.

> **Background**  
>
> - The original strobemer concept and reference implementation live at [ksahlin/strobemers](https://github.com/ksahlin/strobemers).  
> - This crate’s design was guided by the Go port found at [shenwei356/strobemers](https://github.com/shenwei356/strobemers), but the code here is a complete Rust‐native rewrite (zero `unsafe`).  

## Overview

Strobemers are an alternative to contiguous k-mer hashing designed to improve sensitivity in long-read alignment and genome comparison.  
This crate provides:

- **MinStrobes (order 2/3)**  
  Choose subsequent strobes by selecting the minimum hash within a sliding window.

- **RandStrobes (order 2/3)**  
  Choose subsequent strobes by a pseudo-random but reproducible rule derived from hash values.

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
    // Two-level selection: first pick m2 in [i+w_min..i+w_max], then pick m3 in [m2+w_min..m2+w_max].
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

## License

This project is MIT‑licensed (see [LICENSE](LICENSE)).
