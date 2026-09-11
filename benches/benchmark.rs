use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use strobemers_rs::{MinStrobes, RandStrobes};

use rand::{Rng, SeedableRng};
use std::hint::black_box;

const L: usize = 31;
const W_MIN: usize = 20;
const W_MAX: usize = 50;

/// Default working length. Large enough that the precomputed hash array
/// (8 bytes per base, plus as much again for the window minima) no longer fits
/// in L2, which is where genome-scale inputs live.
const LEN: usize = 1_000_000;

/// Generate a reproducible pseudo-random DNA sequence of `n` bases.
fn make_seq(n: usize) -> Vec<u8> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);
    const BASES: [u8; 4] = *b"ACGT";
    (0..n).map(|_| BASES[rng.random_range(0..4)]).collect()
}

/// Both orders of both iterators, driven to exhaustion.
fn bench_iterators(c: &mut Criterion) {
    let seq = make_seq(LEN);
    let mut g = c.benchmark_group("strobemers");
    g.throughput(Throughput::Bytes(LEN as u64));

    for n in [2u8, 3] {
        g.bench_function(format!("MinStrobes/order-{n}"), |b| {
            b.iter(|| {
                let it = MinStrobes::new(&seq, n, L, W_MIN, W_MAX).unwrap();
                let _sum: u64 = black_box(it).sum();
            })
        });
        g.bench_function(format!("RandStrobes/order-{n}"), |b| {
            b.iter(|| {
                let it = RandStrobes::new(&seq, n, L, W_MIN, W_MAX).unwrap();
                let _sum: u64 = black_box(it).sum();
            })
        });
    }
    g.finish();
}

/// Isolates the k-mer hashing stage, which the iterators run up front and which
/// dominates the cheaper MinStrobes configurations.
fn bench_hashing(c: &mut Criterion) {
    let seq = make_seq(LEN);
    let mut g = c.benchmark_group("hashing");
    g.throughput(Throughput::Bytes(LEN as u64));
    g.bench_function("nthash", |b| {
        b.iter(|| {
            // Construction performs the hashing; order 2 adds the least on top.
            let it = RandStrobes::new(&seq, 2, L, W_MIN, W_MAX).unwrap();
            black_box(it);
        })
    });
    g.finish();
}

/// `collect()` rather than `sum()`, so the allocation path (and therefore
/// `size_hint`) is on the critical path.
fn bench_collect(c: &mut Criterion) {
    let seq = make_seq(LEN);
    let mut g = c.benchmark_group("collect");
    g.throughput(Throughput::Bytes(LEN as u64));
    g.bench_function("MinStrobes/order-2", |b| {
        b.iter(|| {
            let v: Vec<u64> = MinStrobes::new(&seq, 2, L, W_MIN, W_MAX).unwrap().collect();
            black_box(v);
        })
    });
    g.finish();
}

/// Per-base cost as the working set outgrows the caches. The iterators keep one
/// `u64` per base (MinStrobes keeps three), so throughput here shows how much of
/// the runtime is memory traffic rather than arithmetic.
fn bench_scaling(c: &mut Criterion) {
    let mut g = c.benchmark_group("scaling/MinStrobes-order-2");
    for len in [
        10_000usize,
        100_000,
        200_000,
        500_000,
        1_000_000,
        10_000_000,
    ] {
        let seq = make_seq(len);
        g.throughput(Throughput::Bytes(len as u64));
        g.bench_with_input(BenchmarkId::from_parameter(len), &seq, |b, seq| {
            b.iter(|| {
                let it = MinStrobes::new(seq, 2, L, W_MIN, W_MAX).unwrap();
                let _sum: u64 = black_box(it).sum();
            })
        });
    }
    g.finish();
}

criterion_group!(
    benches,
    bench_iterators,
    bench_hashing,
    bench_collect,
    bench_scaling
);
criterion_main!(benches);
