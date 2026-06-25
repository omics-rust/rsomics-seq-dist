use std::io::{BufRead, Write};

use rsomics_common::Result;

pub mod cli;
pub mod fasta;
pub mod matrix;
pub mod metric;
pub mod simd;

pub use metric::Metric;

/// Read an aligned FASTA, compute the pairwise distance matrix under `metric`,
/// and write it as a distance-matrix TSV.
pub fn run(metric: Metric, input: &mut dyn BufRead, output: &mut dyn Write) -> Result<()> {
    let aln = fasta::read(input)?;
    let flat = matrix::compute(metric, &aln.seqs);
    matrix::write_tsv(output, &aln.ids, &flat)
}
