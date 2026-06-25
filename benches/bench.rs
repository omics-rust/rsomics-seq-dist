use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use rsomics_seq_dist::Metric;
use rsomics_seq_dist::matrix::compute;

fn synth(n: usize, len: usize) -> Vec<Vec<u8>> {
    let bases = b"ACGT";
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    (0..n)
        .map(|_| (0..len).map(|_| bases[(next() % 4) as usize]).collect())
        .collect()
}

fn bench(c: &mut Criterion) {
    let seqs = synth(300, 800);
    for (name, metric) in [
        ("hamming", Metric::Hamming),
        ("jc69", Metric::Jc69),
        ("k2p", Metric::K2p),
    ] {
        c.bench_function(name, |b| {
            b.iter(|| black_box(compute(metric, black_box(&seqs))));
        });
    }
}

criterion_group!(benches, bench);
criterion_main!(benches);
