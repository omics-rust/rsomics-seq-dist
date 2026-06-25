//! Symmetric distance-matrix computation and TSV serialization.
//!
//! Output layout: a header row of tab-separated ids, then one row per id with
//! `n` tab-separated numeric cells — a full symmetric square matrix. This is the
//! distance-matrix TSV that `rsomics-nj-tree` and `rsomics-upgma` parse (header
//! = ids, body = `n × n` numbers, no row-id prefix), so the output pipes
//! directly into tree construction. The numeric values equal scikit-bio's
//! `DistanceMatrix`; the only layout difference from scikit-bio's own
//! `DistanceMatrix.write` is the absence of its leading-tab/row-label
//! decoration, which the rsomics consumers do not read.

use std::io::Write;

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};

use crate::metric::{Metric, encode, hamming, jc69_from_counts, k2p_from_counts};
use crate::simd::{tally_jc69, tally_k2p};

/// Compute the full symmetric distance matrix in row-major order. The diagonal
/// is exactly 0.0; off-diagonal `(i, j)` is the metric distance between
/// sequences `i` and `j`.
///
/// Work is parallelised over rows: row `i` computes its distances to every
/// `j < i` in one tight inner pass while sequence `i` stays cache-hot. JC69 and
/// K2P pre-encode every sequence to canonical codes once (not per pair), so the
/// inner loop is a branchless tally the compiler auto-vectorizes.
#[must_use]
pub fn compute(metric: Metric, seqs: &[Vec<u8>]) -> Vec<f64> {
    let n = seqs.len();
    let mut flat = vec![0.0f64; n * n];

    // Each row i owns its own output slice (lower triangle); the upper triangle
    // is mirrored afterwards. par_chunks_mut over the n*n buffer gives one row
    // per task with no cross-row aliasing.
    let codes: Vec<Vec<u8>> = match metric {
        Metric::Hamming => Vec::new(),
        Metric::Jc69 | Metric::K2p => seqs.par_iter().map(|s| encode(s)).collect(),
    };

    flat.par_chunks_mut(n)
        .enumerate()
        .for_each(|(i, row)| match metric {
            Metric::Hamming => {
                for (j, cell) in row.iter_mut().enumerate().take(i) {
                    *cell = hamming(&seqs[i], &seqs[j]);
                }
            }
            Metric::Jc69 => {
                let ci = &codes[i];
                for (j, cell) in row.iter_mut().enumerate().take(i) {
                    let (valid, differing) = tally_jc69(ci, &codes[j]);
                    *cell = jc69_from_counts(valid, differing);
                }
            }
            Metric::K2p => {
                let ci = &codes[i];
                for (j, cell) in row.iter_mut().enumerate().take(i) {
                    let (valid, ts, tv) = tally_k2p(ci, &codes[j]);
                    *cell = k2p_from_counts(valid, ts, tv);
                }
            }
        });

    for i in 0..n {
        for j in 0..i {
            flat[j * n + i] = flat[i * n + j];
        }
    }
    flat
}

/// Format a cell to match scikit-bio's float repr where finite, and emit
/// lowercase `nan` / `inf` (saturated / over-divergent pairs) instead of
/// panicking. Finite values use the shortest round-trip representation, which
/// reproduces scikit-bio's `repr(float)` output for the same value.
fn fmt_cell(v: f64) -> String {
    if v.is_nan() {
        "nan".to_string()
    } else if v.is_infinite() {
        if v.is_sign_negative() { "-inf" } else { "inf" }.to_string()
    } else {
        let mut s = format!("{v}");
        if !s.contains(['.', 'e', 'E']) {
            s.push_str(".0");
        }
        s
    }
}

pub fn write_tsv(out: &mut dyn Write, ids: &[String], flat: &[f64]) -> Result<()> {
    let n = ids.len();
    writeln!(out, "{}", ids.join("\t")).map_err(RsomicsError::Io)?;
    let mut line = String::new();
    for i in 0..n {
        line.clear();
        for j in 0..n {
            if j > 0 {
                line.push('\t');
            }
            line.push_str(&fmt_cell(flat[i * n + j]));
        }
        writeln!(out, "{line}").map_err(RsomicsError::Io)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagonal_is_zero_and_symmetric() {
        let seqs = vec![b"ACGT".to_vec(), b"AGGT".to_vec(), b"ACGA".to_vec()];
        let m = compute(Metric::Hamming, &seqs);
        let n = 3;
        for i in 0..n {
            assert_eq!(m[i * n + i], 0.0);
            for j in 0..n {
                assert_eq!(m[i * n + j], m[j * n + i]);
            }
        }
    }

    #[test]
    fn fmt_cell_matches_scikit_bio_repr() {
        assert_eq!(fmt_cell(0.0), "0.0");
        assert_eq!(fmt_cell(0.5), "0.5");
        assert_eq!(fmt_cell(0.188_485_821_210_679_53), "0.18848582121067953");
        assert_eq!(fmt_cell(f64::NAN), "nan");
        assert_eq!(fmt_cell(f64::INFINITY), "inf");
    }
}
