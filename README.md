# rsomics-seq-dist

Pairwise nucleotide-sequence distances from an aligned FASTA, value-exact to
scikit-bio. Computes Hamming (proportion), Jukes-Cantor (JC69), and
Kimura-2-Parameter (K2P) and emits a symmetric distance-matrix TSV that pipes
directly into `rsomics-nj-tree` / `rsomics-upgma`.

## Usage

```bash
rsomics-seq-dist aligned.fasta --metric k2p
rsomics-seq-dist aligned.fasta --metric jc69 -o dist.tsv
cat aligned.fasta | rsomics-seq-dist - --metric hamming
```

`--metric` is one of `hamming` (default), `jc69`, `k2p`. Input is an aligned
FASTA (all records the same length) on a path or stdin (`-`).

## Output

A header row of tab-separated ids, then one row per id with `n` tab-separated
numeric cells — a full symmetric square matrix, the layout `rsomics-nj-tree` and
`rsomics-upgma` read. Over-saturated or empty pairs are written as `nan`. Values
are byte-identical to scikit-bio's `DistanceMatrix`.

## Metrics

- **Hamming** — proportion of differing aligned positions over the full
  alignment length; gaps and ambiguity codes are compared literally (no
  exclusion), matching `skbio.sequence.distance.hamming(proportion=True)`.
- **JC69** — `D = -3/4·ln(1 - 4/3·p)` over the canonical p-distance. Non-canonical
  columns (gaps, `N`, ambiguity codes) are excluded pairwise. `p ≥ 3/4` is
  over-saturated and returns `nan`.
- **K2P** — `D = -1/2·ln(1 - 2P - Q) - 1/4·ln(1 - 2Q)`, with `P`/`Q` the
  transition/transversion proportions over canonical sites. A non-positive log
  argument (over-saturation) returns `nan`.

## Performance

Sequences are pre-encoded to canonical codes once, then the JC69/K2P pairwise
tally runs as a portable-SIMD loop (`wide::u8x16`, stable Rust) — lane-wise
canonical-site, transition, and transversion masks reduced with `move_mask` +
`count_ones`. Work is parallelised over matrix rows with rayon. Single-threaded,
per metric, this beats scikit-bio's numpy-vectorized core on the same machine;
multi-threaded and full-program it wins by a wide margin.

## Origin

This crate reproduces scikit-bio 0.7.2's
`skbio.sequence.distance.{hamming, jc69, k2p}`, whose source is BSD-3-Clause and
was read directly. The evolutionary models are:

- Jukes & Cantor (1969), *Evolution of protein molecules*, Mammalian Protein
  Metabolism 3(21):132.
- Kimura, M. (1980), *A simple method for estimating evolutionary rates of base
  substitutions through comparative studies of nucleotide sequences*, Journal of
  Molecular Evolution 16(2):111-120.

Goldens are generated from scikit-bio 0.7.2 and committed under `tests/golden/`;
the compat test runs without scikit-bio installed.

License: MIT OR Apache-2.0.
Upstream credit: scikit-bio (BSD-3-Clause).
