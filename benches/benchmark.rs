use criterion::{criterion_group, criterion_main, Criterion};
use strobemers_rs::{MinStrobes, RandStrobes};

use rand::{Rng, SeedableRng};
use std::hint::black_box;

const L: usize = 31;
const W_MIN: usize = 20;
const W_MAX: usize = 50;

/// Generate a reproducible 100-kbp pseudo-random DNA sequence.
fn make_seq() -> Vec<u8> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);
    const BASES: [u8; 4] = *b"ACGT";
    (0..100_000)
        .map(|_| BASES[rng.random_range(0..4)])
        .collect()
}

fn bench_minstrobes_iter(c: &mut Criterion) {
    let seq = make_seq();
    c.bench_function("MinStrobes order-3", |b| {
        b.iter(|| {
            let it = MinStrobes::new(&seq, 3, L, W_MIN, W_MAX).unwrap();
            let _sum: u64 = black_box(it).sum();
        })
    });
}

fn bench_randstrobes_iter(c: &mut Criterion) {
    let seq = make_seq();
    c.bench_function("RandStrobes order-3", |b| {
        b.iter(|| {
            let it = RandStrobes::new(&seq, 3, L, W_MIN, W_MAX).unwrap();
            let _sum: u64 = black_box(it).sum();
        })
    });
}

criterion_group!(
    benches,
    bench_minstrobes_iter,
    bench_randstrobes_iter
);
criterion_main!(benches);
